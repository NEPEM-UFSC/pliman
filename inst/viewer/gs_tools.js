/* =============================================================================
 * GS Studio: Gaussian-splatting tools for gs_viewer() (pliman)
 *
 * Loaded after the main viewer script; it uses the viewer globals (scene,
 * camera, renderer, controls, modelContainer, raw* arrays, ...) and plugs into
 * it through a few hooks:
 *   gsToolMask           Uint8Array, 1 = splat hidden (deleted / floater / section)
 *   gsToolsTint()        tints the selection inside updateSplatBuffers()
 *   gsToolsBeforeRender  drives animations / photo comparison every frame
 *   gsToolsPick()        picks the front-most *visible* splat (model coords)
 *   gsToolsOnReindex()   keeps the tool state in sync when the cloud is cropped
 *   gsFmtLen()           formats lengths in calibrated units (also used by the ruler)
 * ========================================================================== */
(function () {
  "use strict";

  // ---------------------------------------------------------------- helpers
  const $ = (id) => document.getElementById(id);
  const N = () => (typeof numGaussians === "number" ? numGaussians : 0);
  const activeCam = () => ((isOrtho && orthoCamera) ? orthoCamera : camera);
  const SH_C0 = 0.28209479177387814;
  const toast = (msg) => { if (typeof showToast === "function") showToast(msg); };
  const hasSplatBuffers = () => (typeof colorsAttr !== "undefined" && colorsAttr && alphasAttr && N() > 0);
  const refreshBuffers = () => { if (hasSplatBuffers()) updateSplatBuffers(); };
  const clamp = (v, a, b) => Math.min(b, Math.max(a, v));

  function download(blob, name) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 2000);
  }

  // ------------------------------------------------------------------ state
  const st = {
    ver: 0,
    deleted: null, delCount: 0,
    sel: null, selCount: 0,
    undo: [],
    sizes: null, sortedSizes: null,
    minOpac: 0, sizePct: 100,
    clip: { on: false, axis: "z", pos: 50, flip: false, normal: null, lo: 0, hi: 1, helper: null, showHelper: true },
    tool: null,          // rect | lasso | brush | polyline | angle | area | calib
    brushR: 30,
    drag: null,
    orig: null,
    anim: null,
    rec: null,
    cmp: { on: false, idx: 0, saved: null, mode: "blend", opacity: 0.5, split: 50 },
    meas: [], cur: null, measGroup: null, hover: null,
    calibPts: [], calibDist: null,
    hidden: null
  };

  function ensureArrays() {
    const n = N();
    if (!st.deleted || st.deleted.length !== n) { st.deleted = new Uint8Array(n); st.delCount = 0; st.undo = []; }
    if (!st.sel || st.sel.length !== n) { st.sel = new Uint8Array(n); st.selCount = 0; }
  }

  // ------------------------------------------------------- linear algebra
  // Jacobi eigen-decomposition of a symmetric 3x3 matrix given as
  // [a00, a01, a02, a11, a12, a22]; returns eigenvalues and eigenvectors (columns).
  function jacobi3(c) {
    const A = [[c[0], c[1], c[2]], [c[1], c[3], c[4]], [c[2], c[4], c[5]]];
    const V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 24; sweep++) {
      const off = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];
      if (off < 1e-30) break;
      for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
        const apq = A[p][q];
        if (Math.abs(apq) < 1e-300) continue;
        const theta = (A[q][q] - A[p][p]) / (2 * apq);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
        for (let k = 0; k < 3; k++) {           // A <- A J
          const akp = A[k][p], akq = A[k][q];
          A[k][p] = cs * akp - sn * akq;
          A[k][q] = sn * akp + cs * akq;
        }
        for (let k = 0; k < 3; k++) {           // A <- J^T A
          const apk = A[p][k], aqk = A[q][k];
          A[p][k] = cs * apk - sn * aqk;
          A[q][k] = sn * apk + cs * aqk;
        }
        for (let k = 0; k < 3; k++) {           // V <- V J
          const vkp = V[k][p], vkq = V[k][q];
          V[k][p] = cs * vkp - sn * vkq;
          V[k][q] = sn * vkp + cs * vkq;
        }
      }
    }
    return { vals: [A[0][0], A[1][1], A[2][2]], vecs: V };
  }

  // largest eigenvalue of a symmetric 3x3 (closed form)
  function maxEig(a00, a01, a02, a11, a12, a22) {
    const p1 = a01 * a01 + a02 * a02 + a12 * a12;
    if (p1 < 1e-30) return Math.max(a00, a11, a22);
    const q = (a00 + a11 + a22) / 3;
    const b00 = a00 - q, b11 = a11 - q, b22 = a22 - q;
    const p = Math.sqrt((b00 * b00 + b11 * b11 + b22 * b22 + 2 * p1) / 6);
    const r = (b00 * (b11 * b22 - a12 * a12) - a01 * (a01 * b22 - a12 * a02) + a02 * (a01 * a12 - b11 * a02)) / (2 * p * p * p);
    const phi = Math.acos(clamp(r, -1, 1)) / 3;
    return q + 2 * p * Math.cos(phi);
  }

  function computeSizes() {
    const n = N();
    const sz = new Float32Array(n);
    if (rawCov3D && rawCov3D.length === n * 6) {
      for (let i = 0; i < n; i++) {
        const c = 6 * i;
        sz[i] = Math.sqrt(Math.max(0, maxEig(rawCov3D[c], rawCov3D[c + 1], rawCov3D[c + 2], rawCov3D[c + 3], rawCov3D[c + 4], rawCov3D[c + 5])));
      }
    } else if (rawScales && rawScales.length === n) {
      sz.set(rawScales);
    }
    st.sizes = sz;
    st.sortedSizes = Float32Array.from(sz).sort();
  }

  // ---------------------------------------------------------- visibility mask
  function clipNormal() {
    const ax = st.clip.axis;
    if (ax === "x") return new THREE.Vector3(1, 0, 0);
    if (ax === "y") return new THREE.Vector3(0, 1, 0);
    if (ax === "z") return new THREE.Vector3(0, 0, 1);
    return st.clip.normal ? st.clip.normal.clone() : new THREE.Vector3(0, 0, 1);
  }

  function clipRange() {
    const n = N(), nv = clipNormal();
    const step = Math.max(1, Math.floor(n / 40000));
    const v = [];
    for (let i = 0; i < n; i += step) v.push(nv.x * rawPositions[3 * i] + nv.y * rawPositions[3 * i + 1] + nv.z * rawPositions[3 * i + 2]);
    v.sort((a, b) => a - b);
    st.clip.lo = v[Math.floor(v.length * 0.005)] || 0;
    st.clip.hi = v[Math.min(v.length - 1, Math.floor(v.length * 0.995))] || 1;
  }

  function rebuildMask() {
    if (N() === 0) return;
    ensureArrays();
    const n = N();
    const useOpac = st.minOpac > 0;
    const useSize = st.sizePct < 100;
    const useClip = st.clip.on;
    let cntF = 0, cntC = 0;
    if (!useOpac && !useSize && !useClip && st.delCount === 0) {
      gsToolMask = null;
    } else {
      if ((useSize) && (!st.sizes || st.sizes.length !== n)) computeSizes();
      const mask = new Uint8Array(n);
      const sizeCut = useSize ? st.sortedSizes[Math.min(n - 1, Math.floor(st.sizePct / 100 * n))] : Infinity;
      const nv = clipNormal();
      const d = st.clip.lo + (st.clip.hi - st.clip.lo) * st.clip.pos / 100;
      const flip = st.clip.flip;
      for (let i = 0; i < n; i++) {
        let m = st.deleted[i];
        if (!m && (useOpac || useSize)) {
          if ((useOpac && rawOpacities[i] < st.minOpac) || (useSize && st.sizes[i] > sizeCut)) { m = 1; cntF++; }
        }
        if (!m && useClip) {
          const dd = nv.x * rawPositions[3 * i] + nv.y * rawPositions[3 * i + 1] + nv.z * rawPositions[3 * i + 2] - d;
          if (flip ? dd < 0 : dd > 0) { m = 1; cntC++; }
        }
        mask[i] = m;
      }
      gsToolMask = mask;
    }
    updateClipHelper();
    refreshBuffers();
    const lf = $("gss-floater-stat");
    if (lf) lf.textContent = (useOpac || useSize) ? `${cntF.toLocaleString()} floaters hidden` : "off";
    const lc = $("gss-clip-stat");
    if (lc) lc.textContent = useClip ? `${cntC.toLocaleString()} splats cut` : "off";
    updateEditStats();
  }

  // ---------------------------------------------------------------- hooks
  window.gsToolsTint = function (col, al) {
    if (!st.sel || st.selCount === 0 || st.sel.length !== al.length) return;
    for (let i = 0; i < al.length; i++) {
      if (st.sel[i] && al[i] > 0) {
        col[3 * i] = col[3 * i] * 0.35 + 0.65;
        col[3 * i + 1] = col[3 * i + 1] * 0.35 + 0.33;
        col[3 * i + 2] = col[3 * i + 2] * 0.35 + 0.03;
      }
    }
  };

  window.gsToolsOnReindex = function (kept) {
    st.ver++;
    proj.key = "";
    st.sizes = null;
    gsToolMask = null;
    if (kept === "reset") {
      st.deleted = null; st.sel = null; st.orig = null; st.undo = []; st.delCount = 0; st.selCount = 0;
    } else if (kept === null) {
      if (st.orig) { st.deleted = st.orig.deleted; st.delCount = st.orig.delCount; }
      st.orig = null; st.sel = null; st.selCount = 0; st.undo = [];
    } else if (st.deleted && st.deleted.length > 0) {
      if (!st.orig) st.orig = { deleted: new Uint8Array(st.deleted), delCount: st.delCount };
      const nd = new Uint8Array(kept.length);
      const ns = new Uint8Array(kept.length);
      let dc = 0, sc = 0;
      for (let j = 0; j < kept.length; j++) {
        nd[j] = st.deleted[kept[j]]; dc += nd[j];
        ns[j] = st.sel ? st.sel[kept[j]] : 0; sc += ns[j];
      }
      st.deleted = nd; st.delCount = dc; st.sel = ns; st.selCount = sc; st.undo = [];
    }
    if (st.clip.on) setTimeout(() => { clipRange(); rebuildMask(); }, 0);
    else setTimeout(rebuildMask, 0);
  };

  // -------------------------------------------------- screen projection cache
  const proj = { key: "", sx: null, sy: null, sz: null, M: new THREE.Matrix4() };
  function project() {
    const n = N();
    const cam = activeCam();
    cam.updateMatrixWorld();
    modelContainer.updateMatrixWorld();
    proj.M.multiplyMatrices(cam.projectionMatrix, cam.matrixWorldInverse).multiply(modelContainer.matrixWorld);
    const e = proj.M.elements;
    const W = window.innerWidth, H = window.innerHeight;
    const key = e.join(",") + "|" + W + "|" + H + "|" + n + "|" + st.ver;
    if (key === proj.key && proj.sx && proj.sx.length === n) return proj;
    if (!proj.sx || proj.sx.length !== n) {
      proj.sx = new Float32Array(n); proj.sy = new Float32Array(n); proj.sz = new Float32Array(n);
    }
    const P = rawPositions;
    for (let i = 0; i < n; i++) {
      const x = P[3 * i], y = P[3 * i + 1], z = P[3 * i + 2];
      const w = e[3] * x + e[7] * y + e[11] * z + e[15];
      if (w <= 1e-9) { proj.sx[i] = NaN; continue; }
      proj.sx[i] = ((e[0] * x + e[4] * y + e[8] * z + e[12]) / w + 1) * 0.5 * W;
      proj.sy[i] = (1 - (e[1] * x + e[5] * y + e[9] * z + e[13]) / w) * 0.5 * H;
      proj.sz[i] = (e[2] * x + e[6] * y + e[10] * z + e[14]) / w;
    }
    proj.key = key;
    return proj;
  }

  // Front-most visible splat centre under the cursor (model coordinates)
  window.gsToolsPick = function (cx, cy) {
    if (N() === 0 || !hasSplatBuffers()) return null;
    const P = project();
    const al = alphasAttr.array;
    const n = N();
    const r2 = 64;
    let best = -1, bestZ = Infinity, best2 = -1, bestZ2 = Infinity;
    for (let i = 0; i < n; i++) {
      const a = al[i];
      if (a <= 0) continue;
      const dx = P.sx[i] - cx, dy = P.sy[i] - cy;
      if (!(dx * dx + dy * dy <= r2)) continue;
      const z = P.sz[i];
      if (z < -1 || z > 1) continue;
      if (a >= 0.3) { if (z < bestZ) { bestZ = z; best = i; } }
      else if (z < bestZ2) { bestZ2 = z; best2 = i; }
    }
    const b = best >= 0 ? best : best2;
    if (b < 0) return null;
    return new THREE.Vector3(rawPositions[3 * b], rawPositions[3 * b + 1], rawPositions[3 * b + 2]);
  };

  // ------------------------------------------------------------- units
  function fmtLen(d) {
    if (!(d >= 0)) return "-";
    if (gsUnitScale > 0) {
      const m = d * gsUnitScale;
      if (m < 0.01) return (m * 1000).toFixed(2) + " mm";
      if (m < 1) return (m * 100).toFixed(2) + " cm";
      return m.toFixed(3) + " m";
    }
    return d.toFixed(4) + " u";
  }
  function fmtArea(a) {
    if (gsUnitScale > 0) {
      const m2 = a * gsUnitScale * gsUnitScale;
      if (m2 < 1e-4) return (m2 * 1e6).toFixed(2) + " mm\u00b2";
      if (m2 < 1) return (m2 * 1e4).toFixed(2) + " cm\u00b2";
      return m2.toFixed(4) + " m\u00b2";
    }
    return a.toFixed(5) + " u\u00b2";
  }
  function fmtVol(v) {
    if (gsUnitScale > 0) {
      const m3 = v * Math.pow(gsUnitScale, 3);
      if (m3 < 1e-6) return (m3 * 1e9).toFixed(1) + " mm\u00b3";
      if (m3 < 0.1) return (m3 * 1e6).toFixed(2) + " cm\u00b3";
      return m3.toFixed(4) + " m\u00b3";
    }
    return v.toFixed(5) + " u\u00b3";
  }
  window.gsFmtLen = fmtLen;

  function refreshUnits() {
    const lbl = $("gss-unit-status");
    if (lbl) {
      lbl.innerHTML = gsUnitScale > 0
        ? `1 model unit = <b>${fmtLen(1)}</b>`
        : `<span style="color:#facc15">Unknown scale</span> (e.g. COLMAP): lengths in model units (u). Calibrate below.`;
    }
    // re-label the viewer ruler
    if (typeof rulerMeasurements !== "undefined" && Array.isArray(rulerMeasurements)) {
      rulerMeasurements.forEach((m) => {
        if (typeof m.distM !== "number") return;
        m.distStr = fmtLen(m.distM);
        m.dxStr = fmtLen(m.dx); m.dyStr = fmtLen(m.dy); m.dzStr = fmtLen(m.dz);
        if (m.sprite) {
          rulerGroup.remove(m.sprite);
          m.sprite = createMeasurementBadge(m.distStr, m.pA.clone().add(m.pB).multiplyScalar(0.5));
          rulerGroup.add(m.sprite);
        }
      });
      if (typeof updateMeasurementsListUI === "function") updateMeasurementsListUI();
    }
    st.meas.forEach(relabelMeasurement);
    renderMeasList();
  }

  // ================================================================ RENDER
  function setAdjust(key, v) {
    gsAdjust[key] = v;
    const o = $("gss-adj-" + key + "-v");
    if (o) o.textContent = (key === "exposure" ? (v >= 0 ? "+" : "") : "") + Number(v).toFixed(2);
  }

  function setBackground(hex) {
    if (scene && scene.background && scene.background.isColor) scene.background.set(hex);
    else if (scene) scene.background = new THREE.Color(hex);
    const ci = $("gss-bg-color");
    if (ci) ci.value = hex;
  }

  // ================================================================== EDIT
  function updateEditStats() {
    const s = $("gss-edit-stat");
    if (!s) return;
    let vis = 0;
    if (hasSplatBuffers()) { const al = alphasAttr.array; for (let i = 0; i < al.length; i++) if (al[i] > 0) vis++; }
    s.textContent = `Selected: ${st.selCount.toLocaleString()}   Deleted: ${st.delCount.toLocaleString()}   Visible: ${vis.toLocaleString()}`;
  }

  function setTool(t) {
    if (st.cur) {
      if (st.cur.pts.length >= minPts(st.cur.type) && st.cur.type !== "angle") finishMeasurement();
      else cancelMeasurement();
    }
    st.tool = (st.tool === t) ? null : t;
    document.querySelectorAll("[data-gss-tool]").forEach((b) => b.classList.toggle("active", b.getAttribute("data-gss-tool") === st.tool));
    if (st.tool && typeof isRulerActive !== "undefined" && isRulerActive && typeof toggleRulerTool === "function") toggleRulerTool();
    renderer.domElement.style.cursor = st.tool ? "crosshair" : "default";
    if (st.tool === "calib") { st.calibPts = []; clearCalibMarkers(); setStatus("gss-calib-status", "Click the 1st point of a known distance..."); }
    const hints = {
      rect: "Drag a rectangle. Shift = add, Alt = remove. Right button pans, wheel zooms.",
      lasso: "Draw a lasso. Shift = add, Alt = remove.",
      brush: "Paint to add to the selection (Alt = remove). Brush size below.",
      polyline: "Click points along the path. Finish: double-click, Enter, right-click or the Finish button.",
      angle: "Click 3 points: start, vertex, end.",
      area: "Click the polygon vertices. Close: double-click, Enter, right-click or the Finish button.",
      calib: "Click two points whose real distance you know."
    };
    setStatus("gss-hint", st.tool ? hints[st.tool] : "");
    drawOverlay();
  }

  function setStatus(id, txt) { const e = $(id); if (e) e.textContent = txt || ""; }

  function applySelection(idxTest, mode) {
    ensureArrays();
    const P = project();
    const al = alphasAttr.array;
    const n = N();
    if (mode === "replace") { st.sel.fill(0); }
    for (let i = 0; i < n; i++) {
      if (al[i] <= 0) continue;
      const x = P.sx[i];
      if (x !== x) continue;
      if (idxTest(x, P.sy[i])) st.sel[i] = (mode === "remove") ? 0 : 1;
    }
    let c = 0;
    for (let i = 0; i < n; i++) c += st.sel[i];
    st.selCount = c;
    refreshBuffers();
    updateEditStats();
  }

  function pointInPoly(x, y, poly) {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1];
      if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi)) inside = !inside;
    }
    return inside;
  }

  function opDelete() {
    ensureArrays();
    if (st.selCount === 0) { toast("Nothing selected"); return; }
    const idx = [];
    const al = alphasAttr.array;
    for (let i = 0; i < st.sel.length; i++) if (st.sel[i] && !st.deleted[i] && al[i] > 0) { st.deleted[i] = 1; idx.push(i); }
    st.undo.push(Uint32Array.from(idx));
    st.delCount += idx.length;
    st.sel.fill(0); st.selCount = 0;
    rebuildMask();
    toast(`Deleted ${idx.length.toLocaleString()} splats`);
  }

  function opIsolate() {
    ensureArrays();
    if (st.selCount === 0) { toast("Select what to keep first"); return; }
    const idx = [];
    const al = alphasAttr.array;
    for (let i = 0; i < st.sel.length; i++) if (!st.sel[i] && !st.deleted[i] && al[i] > 0) { st.deleted[i] = 1; idx.push(i); }
    st.undo.push(Uint32Array.from(idx));
    st.delCount += idx.length;
    st.sel.fill(0); st.selCount = 0;
    rebuildMask();
    toast(`Isolated selection (${idx.length.toLocaleString()} splats removed)`);
  }

  function opInvert() {
    ensureArrays();
    const al = alphasAttr.array;
    let c = 0;
    for (let i = 0; i < st.sel.length; i++) { st.sel[i] = (al[i] > 0 && !st.sel[i]) ? 1 : 0; c += st.sel[i]; }
    st.selCount = c;
    refreshBuffers(); updateEditStats();
  }

  function opSelectAll() {
    ensureArrays();
    const al = alphasAttr.array;
    let c = 0;
    for (let i = 0; i < st.sel.length; i++) { st.sel[i] = al[i] > 0 ? 1 : 0; c += st.sel[i]; }
    st.selCount = c;
    refreshBuffers(); updateEditStats();
  }

  function opClearSel() {
    ensureArrays();
    st.sel.fill(0); st.selCount = 0;
    refreshBuffers(); updateEditStats();
  }

  function opUndo() {
    const last = st.undo.pop();
    if (!last) { toast("Nothing to undo"); return; }
    for (let k = 0; k < last.length; k++) st.deleted[last[k]] = 0;
    st.delCount -= last.length;
    rebuildMask();
    toast(`Restored ${last.length.toLocaleString()} splats`);
  }

  function opResetEdits() {
    ensureArrays();
    st.deleted.fill(0); st.delCount = 0; st.undo = [];
    st.sel.fill(0); st.selCount = 0;
    rebuildMask();
    toast("All deletions reverted");
  }

  // ---- section plane
  function updateClipHelper() {
    if (!st.clip.on || !st.clip.showHelper) {
      if (st.clip.helper) st.clip.helper.visible = false;
      return;
    }
    if (!st.clip.helper) {
      const g = new THREE.PlaneGeometry(1, 1);
      const m = new THREE.MeshBasicMaterial({ color: 0x38bdf8, transparent: true, opacity: 0.13, side: THREE.DoubleSide, depthWrite: false });
      st.clip.helper = new THREE.Mesh(g, m);
      const edges = new THREE.LineSegments(new THREE.EdgesGeometry(g), new THREE.LineBasicMaterial({ color: 0x38bdf8 }));
      st.clip.helper.add(edges);
      modelContainer.add(st.clip.helper);
    }
    const nv = clipNormal();
    const d = st.clip.lo + (st.clip.hi - st.clip.lo) * st.clip.pos / 100;
    const size = Math.max(0.05, targetDim * 2.6);
    st.clip.helper.scale.set(size, size, 1);
    st.clip.helper.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), nv);
    st.clip.helper.position.copy(nv.clone().multiplyScalar(d));
    st.clip.helper.visible = true;
  }

  function setClipAxis(ax) {
    st.clip.axis = ax;
    if (ax === "view") {
      const cam = activeCam();
      const dir = new THREE.Vector3(0, 0, -1).applyQuaternion(cam.quaternion);
      const qInv = modelContainer.quaternion.clone().invert();
      st.clip.normal = dir.applyQuaternion(qInv).normalize();
    }
    document.querySelectorAll("[data-gss-axis]").forEach((b) => b.classList.toggle("active", b.getAttribute("data-gss-axis") === ax));
    clipRange();
    rebuildMask();
  }

  // ---- export (full 3DGS PLY of the visible splats)
  function exportPLY3DGS() {
    const n = N();
    if (!hasSplatBuffers()) { toast("No model loaded"); return; }
    if (!rawCov3D || rawCov3D.length !== n * 6) {
      if (typeof exportFilteredPLY === "function") exportFilteredPLY();
      return;
    }
    const al = alphasAttr.array;
    const idx = [];
    for (let i = 0; i < n; i++) if (al[i] > 0) idx.push(i);
    if (idx.length === 0) { toast("No visible splats to export"); return; }
    const deg = (gsSHData && gsSHTexels > 0 && gsSHData.length >= n * gsSHTexels * 4) ? gsSHDegree : 0;
    const Kr = (deg + 1) * (deg + 1) - 1;
    const props = ["x", "y", "z", "nx", "ny", "nz", "f_dc_0", "f_dc_1", "f_dc_2"];
    for (let k = 0; k < 3 * Kr; k++) props.push("f_rest_" + k);
    props.push("opacity", "scale_0", "scale_1", "scale_2", "rot_0", "rot_1", "rot_2", "rot_3");
    const header = "ply\nformat binary_little_endian 1.0\ncomment exported by pliman gs_viewer (GS Studio)\nelement vertex " + idx.length + "\n" +
      props.map((p) => "property float " + p).join("\n") + "\nend_header\n";
    const hb = new TextEncoder().encode(header);
    const np = props.length;
    const buf = new ArrayBuffer(hb.length + idx.length * np * 4);
    new Uint8Array(buf).set(hb, 0);
    const dv = new DataView(buf, hb.length);
    const T4 = gsSHTexels * 4;
    const bandScale = (j) => (j <= 3 ? gsSHScale[0] : (j <= 8 ? gsSHScale[1] : gsSHScale[2])) / 127;
    const m4 = new THREE.Matrix4(), q = new THREE.Quaternion();
    const v0 = new THREE.Vector3(), v1 = new THREE.Vector3(), v2 = new THREE.Vector3();
    let off = 0;
    const put = (v) => { dv.setFloat32(off, v, true); off += 4; };
    for (let t = 0; t < idx.length; t++) {
      const i = idx[t];
      put(rawPositions[3 * i] + centerOrigin[0]);
      put(rawPositions[3 * i + 1] + centerOrigin[1]);
      put(rawPositions[3 * i + 2] + centerOrigin[2]);
      put(0); put(0); put(0);
      for (let c = 0; c < 3; c++) put((rawColorsRGB[3 * i + c] - 0.5) / SH_C0);
      for (let c = 0; c < 3; c++) {
        for (let j = 1; j <= Kr; j++) {
          const b = gsSHData[i * T4 + (j - 1) * 3 + c];
          put((b - 128) * bandScale(j));
        }
      }
      const a = clamp(rawOpacities[i], 1e-4, 1 - 1e-4);
      put(Math.log(a / (1 - a)));
      const e = jacobi3([rawCov3D[6 * i], rawCov3D[6 * i + 1], rawCov3D[6 * i + 2], rawCov3D[6 * i + 3], rawCov3D[6 * i + 4], rawCov3D[6 * i + 5]]);
      for (let k = 0; k < 3; k++) put(0.5 * Math.log(Math.max(e.vals[k], 1e-20)));
      const V = e.vecs;
      v0.set(V[0][0], V[1][0], V[2][0]); v1.set(V[0][1], V[1][1], V[2][1]); v2.set(V[0][2], V[1][2], V[2][2]);
      if (v0.clone().cross(v1).dot(v2) < 0) v2.negate();
      m4.makeBasis(v0, v1, v2);
      q.setFromRotationMatrix(m4).normalize();
      put(q.w); put(q.x); put(q.y); put(q.z);
    }
    download(new Blob([buf], { type: "application/octet-stream" }), "gs_edited.ply");
    toast(`Exported ${idx.length.toLocaleString()} splats (3DGS PLY, SH degree ${deg})`);
  }

  // =============================================================== MEASURE
  function markerMesh(p, color) {
    const r = Math.max(0.004, targetDim * 0.011);
    const m = new THREE.Mesh(new THREE.SphereGeometry(r, 12, 12), new THREE.MeshBasicMaterial({ color: color, depthTest: false }));
    m.position.copy(p);
    m.renderOrder = 999;
    return m;
  }

  function measGroup() {
    if (!st.measGroup) { st.measGroup = new THREE.Group(); modelContainer.add(st.measGroup); }
    return st.measGroup;
  }

  function polyLine(pts, closed, color) {
    const arr = closed ? pts.concat([pts[0]]) : pts;
    const g = new THREE.BufferGeometry().setFromPoints(arr);
    const l = new THREE.Line(g, new THREE.LineBasicMaterial({ color: color, depthTest: false, transparent: true }));
    l.renderOrder = 998;
    return l;
  }

  // Filled polygon (salmon, alpha 0.5) for area measurements: the 3D polygon is
  // projected onto its best-fit plane (Newell normal) and triangulated.
  function fillMesh(pts) {
    const n = new THREE.Vector3();
    for (let i = 0; i < pts.length; i++) {
      const a = pts[i], b = pts[(i + 1) % pts.length];
      n.x += (a.y - b.y) * (a.z + b.z);
      n.y += (a.z - b.z) * (a.x + b.x);
      n.z += (a.x - b.x) * (a.y + b.y);
    }
    if (n.lengthSq() < 1e-30) return null;
    n.normalize();
    const u = new THREE.Vector3().subVectors(pts[1], pts[0]);
    u.sub(n.clone().multiplyScalar(u.dot(n)));
    if (u.lengthSq() < 1e-30) return null;
    u.normalize();
    const v = new THREE.Vector3().crossVectors(n, u);
    const contour = pts.map((p) => new THREE.Vector2(p.clone().sub(pts[0]).dot(u), p.clone().sub(pts[0]).dot(v)));
    const tris = THREE.ShapeUtils.triangulateShape(contour, []);
    if (!tris.length) return null;
    const g = new THREE.BufferGeometry().setFromPoints(pts);
    g.setIndex([].concat.apply([], tris));
    const m = new THREE.Mesh(g, new THREE.MeshBasicMaterial({
      color: 0xfa8072, transparent: true, opacity: 0.5, side: THREE.DoubleSide, depthTest: false, depthWrite: false
    }));
    m.renderOrder = 997;
    return m;
  }

  function measValue(m) {
    const p = m.pts;
    if (m.type === "polyline") {
      let L = 0;
      for (let i = 1; i < p.length; i++) L += p[i].distanceTo(p[i - 1]);
      return { text: fmtLen(L), detail: `${p.length} pts`, raw: L };
    }
    if (m.type === "angle") {
      const a = p[0].clone().sub(p[1]), b = p[2].clone().sub(p[1]);
      const ang = Math.acos(clamp(a.dot(b) / (a.length() * b.length() + 1e-30), -1, 1)) * 180 / Math.PI;
      return { text: ang.toFixed(1) + "\u00b0", detail: `arms ${fmtLen(a.length())} / ${fmtLen(b.length())}`, raw: ang };
    }
    // area (vector area of the 3D polygon) + perimeter
    const cr = new THREE.Vector3(), acc = new THREE.Vector3();
    for (let i = 1; i < p.length - 1; i++) {
      cr.crossVectors(p[i].clone().sub(p[0]), p[i + 1].clone().sub(p[0]));
      acc.add(cr);
    }
    const A = 0.5 * acc.length();
    let per = 0;
    for (let i = 0; i < p.length; i++) per += p[i].distanceTo(p[(i + 1) % p.length]);
    return { text: fmtArea(A), detail: `perimeter ${fmtLen(per)}`, raw: A };
  }

  function measAnchor(m) {
    if (m.type === "angle") return m.pts[1].clone();
    const c = new THREE.Vector3();
    m.pts.forEach((q) => c.add(q));
    return c.multiplyScalar(1 / m.pts.length);
  }

  function relabelMeasurement(m) {
    if (!m.done) return;
    const g = measGroup();
    if (m.badge) g.remove(m.badge);
    m.val = measValue(m);
    m.badge = createMeasurementBadge(m.val.text, measAnchor(m));
    g.add(m.badge);
  }

  function addMeasPoint(p) {
    const g = measGroup();
    if (!st.cur) st.cur = { type: st.tool, pts: [], objs: [], done: false };
    const cur = st.cur;
    if (cur.pts.length > 0 && cur.pts[cur.pts.length - 1].distanceTo(p) < targetDim * 1e-4) return;
    cur.pts.push(p.clone());
    const mk = markerMesh(p, cur.pts.length === 1 ? 0x22c55e : 0xfacc15);
    g.add(mk); cur.objs.push(mk);
    if (cur.line) { g.remove(cur.line); }
    if (cur.pts.length > 1) { cur.line = polyLine(cur.pts, cur.type === "area" && cur.pts.length > 2, 0xfacc15); g.add(cur.line); }
    if (cur.type === "area") {
      if (cur.fill) { g.remove(cur.fill); cur.fill = null; }
      if (cur.pts.length > 2) { cur.fill = fillMesh(cur.pts); if (cur.fill) g.add(cur.fill); }
    }
    if (cur.type === "angle" && cur.pts.length === 3) { finishMeasurement(); return; }
    updateLive();
  }

  function minPts(type) { return type === "polyline" ? 2 : 3; }

  // running value of the measurement being drawn (panel + 3D label)
  function updateLive() {
    const cur = st.cur;
    const box = $("gss-meas-live");
    if (!cur) { if (box) box.style.display = "none"; return; }
    const g = measGroup();
    if (cur.badge) { g.remove(cur.badge); cur.badge = null; }
    let txt = `${cur.type}: ${cur.pts.length} point${cur.pts.length === 1 ? "" : "s"}`;
    if (cur.pts.length >= minPts(cur.type)) {
      const v = measValue(cur);
      txt += ` \u2014 ${v.text} (${v.detail})`;
      cur.badge = createMeasurementBadge(v.text, measAnchor(cur));
      g.add(cur.badge);
    } else {
      txt += ` \u2014 need ${minPts(cur.type)}`;
    }
    if (box) {
      box.style.display = "flex";
      const t = $("gss-meas-live-txt");
      if (t) t.textContent = txt;
      const f = $("gss-meas-finish");
      if (f) f.disabled = cur.pts.length < minPts(cur.type) || cur.type === "angle";
    }
  }

  function finishMeasurement() {
    const cur = st.cur;
    if (!cur) return;
    const need = cur.type === "polyline" ? 2 : 3;
    if (cur.pts.length < need) { toast(`Need at least ${need} points`); return; }
    cur.done = true;
    cur.id = (st.meas.length ? st.meas[st.meas.length - 1].id : 0) + 1;
    relabelMeasurement(cur);
    st.meas.push(cur);
    st.cur = null;
    updateLive();
    renderMeasList();
    toast(`${cur.type}: ${cur.val.text}`);
  }

  function cancelMeasurement() {
    const cur = st.cur;
    if (!cur) return;
    const g = measGroup();
    cur.objs.forEach((o) => g.remove(o));
    if (cur.line) g.remove(cur.line);
    if (cur.badge) g.remove(cur.badge);
    if (cur.fill) g.remove(cur.fill);
    st.cur = null;
    updateLive();
  }

  function removeMeas(k) {
    const m = st.meas[k];
    if (!m) return;
    const g = measGroup();
    m.objs.forEach((o) => g.remove(o));
    if (m.line) g.remove(m.line);
    if (m.badge) g.remove(m.badge);
    if (m.fill) g.remove(m.fill);
    st.meas.splice(k, 1);
    renderMeasList();
  }

  function renderMeasList() {
    const L = $("gss-meas-list");
    if (!L) return;
    if (st.meas.length === 0) { L.innerHTML = `<div class="gss-note">No measurements yet.</div>`; return; }
    L.innerHTML = st.meas.map((m, k) =>
      `<div class="gss-item"><span>#${m.id} ${m.type}: <b>${m.val.text}</b> <span class="gss-note">${m.val.detail}</span></span><button data-gss-rm="${k}" title="Delete">&#10005;</button></div>`
    ).join("");
    L.querySelectorAll("[data-gss-rm]").forEach((b) => b.addEventListener("click", () => removeMeas(+b.getAttribute("data-gss-rm"))));
  }

  function copyMeasCSV() {
    const unit = gsUnitScale > 0 ? "calibrated" : "model units";
    const rows = ["id,type,value,detail,points"];
    st.meas.forEach((m) => rows.push(`${m.id},${m.type},"${m.val.text}","${m.val.detail}",${m.pts.length}`));
    (typeof rulerMeasurements !== "undefined" ? rulerMeasurements : []).forEach((m) => rows.push(`ruler${m.id},distance,"${m.distStr}","dx ${m.dxStr} dy ${m.dyStr} dz ${m.dzStr}",2`));
    const txt = rows.join("\n");
    if (navigator.clipboard) navigator.clipboard.writeText(txt).then(() => toast(`Measurements copied (${unit})`), () => download(new Blob([txt]), "measurements.csv"));
    else download(new Blob([txt]), "measurements.csv");
  }

  // ---- scale calibration
  function clearCalibMarkers() {
    const g = measGroup();
    (st.calibObjs || []).forEach((o) => g.remove(o));
    st.calibObjs = [];
  }

  function addCalibPoint(p) {
    const g = measGroup();
    st.calibPts.push(p.clone());
    const mk = markerMesh(p, 0x38bdf8);
    g.add(mk); st.calibObjs.push(mk);
    if (st.calibPts.length === 2) {
      const ln = polyLine(st.calibPts, false, 0x38bdf8);
      g.add(ln); st.calibObjs.push(ln);
      st.calibDist = st.calibPts[0].distanceTo(st.calibPts[1]);
      setStatus("gss-calib-status", `Model distance: ${st.calibDist.toFixed(5)} u. Enter the real length and press Apply.`);
      setTool("calib");   // toggles the tool off
    } else {
      setStatus("gss-calib-status", "Now click the 2nd point...");
    }
  }

  function applyCalibration() {
    const v = parseFloat(($("gss-calib-len") || {}).value);
    const unit = ($("gss-calib-unit") || {}).value || "cm";
    if (!(st.calibDist > 0)) { toast("Pick two points first"); return; }
    if (!(v > 0)) { toast("Enter the real length"); return; }
    const toM = { mm: 0.001, cm: 0.01, m: 1 }[unit];
    gsUnitScale = (v * toM) / st.calibDist;
    clearCalibMarkers();
    refreshUnits();
    setStatus("gss-calib-status", `Calibrated: 1 u = ${fmtLen(1)}`);
    toast("Scale calibrated");
  }

  // ---- volume / dimensions of the selection (or of all visible splats)
  function computeVolume() {
    if (!hasSplatBuffers()) return;
    const n = N(), al = alphasAttr.array;
    const useSel = st.sel && st.selCount > 0;
    const pts = [];
    for (let i = 0; i < n; i++) if (al[i] > 0 && (!useSel || st.sel[i])) pts.push(i);
    if (pts.length < 20) { setStatus("gss-vol-out", "Too few splats."); return; }
    const res = parseInt(($("gss-vol-res") || {}).value || "96", 10);
    const P = rawPositions;
    let mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
    let mean = [0, 0, 0];
    for (const i of pts) for (let k = 0; k < 3; k++) {
      const v = P[3 * i + k];
      if (v < mn[k]) mn[k] = v;
      if (v > mx[k]) mx[k] = v;
      mean[k] += v;
    }
    mean = mean.map((v) => v / pts.length);
    const ext = [mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]];
    const h = Math.max(ext[0], ext[1], ext[2]) / res;
    const nx = Math.ceil(ext[0] / h) + 3, ny = Math.ceil(ext[1] / h) + 3, nz = Math.ceil(ext[2] / h) + 3;
    const idx3 = (x, y, z) => x + nx * (y + ny * z);
    const occ = new Uint8Array(nx * ny * nz);
    for (const i of pts) {
      const x = Math.floor((P[3 * i] - mn[0]) / h) + 1, y = Math.floor((P[3 * i + 1] - mn[1]) / h) + 1, z = Math.floor((P[3 * i + 2] - mn[2]) / h) + 1;
      occ[idx3(x, y, z)] = 1;
    }
    // close small gaps (dilation), fill the interior (flood fill of the exterior), erode back
    const dil = new Uint8Array(occ.length);
    for (let z = 1; z < nz - 1; z++) for (let y = 1; y < ny - 1; y++) for (let x = 1; x < nx - 1; x++) {
      if (!occ[idx3(x, y, z)]) continue;
      for (let dz = -1; dz <= 1; dz++) for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) dil[idx3(x + dx, y + dy, z + dz)] = 1;
    }
    const ext3 = new Uint8Array(occ.length);
    const queue = new Int32Array(occ.length);
    let qh = 0, qt = 0;
    queue[qt++] = 0; ext3[0] = 1;
    while (qh < qt) {
      const c = queue[qh++];
      const x = c % nx, y = Math.floor(c / nx) % ny, z = Math.floor(c / (nx * ny));
      const nb = [[x - 1, y, z], [x + 1, y, z], [x, y - 1, z], [x, y + 1, z], [x, y, z - 1], [x, y, z + 1]];
      for (const [a, b, d] of nb) {
        if (a < 0 || b < 0 || d < 0 || a >= nx || b >= ny || d >= nz) continue;
        const k = idx3(a, b, d);
        if (!ext3[k] && !dil[k]) { ext3[k] = 1; queue[qt++] = k; }
      }
    }
    let solid = 0;
    for (let z = 1; z < nz - 1; z++) for (let y = 1; y < ny - 1; y++) for (let x = 1; x < nx - 1; x++) {
      if (ext3[idx3(x, y, z)]) continue;
      let interior = true;
      for (let dz = -1; dz <= 1 && interior; dz++) for (let dy = -1; dy <= 1 && interior; dy++) for (let dx = -1; dx <= 1; dx++) {
        if (ext3[idx3(x + dx, y + dy, z + dz)]) { interior = false; break; }
      }
      if (interior || occ[idx3(x, y, z)]) solid++;
    }
    const vol = solid * h * h * h;
    // principal dimensions (PCA)
    let c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
    for (const i of pts) {
      const dx = P[3 * i] - mean[0], dy = P[3 * i + 1] - mean[1], dz = P[3 * i + 2] - mean[2];
      c00 += dx * dx; c01 += dx * dy; c02 += dx * dz; c11 += dy * dy; c12 += dy * dz; c22 += dz * dz;
    }
    const e = jacobi3([c00, c01, c02, c11, c12, c22]);
    const dims = [0, 1, 2].map((k) => {
      const ax = [e.vecs[0][k], e.vecs[1][k], e.vecs[2][k]];
      let lo = Infinity, hi = -Infinity;
      for (const i of pts) {
        const v = ax[0] * P[3 * i] + ax[1] * P[3 * i + 1] + ax[2] * P[3 * i + 2];
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
      return hi - lo;
    }).sort((a, b) => b - a);
    setStatus("gss-vol-out",
      `${useSel ? "Selection" : "Visible splats"}: ${pts.length.toLocaleString()}\n` +
      `Volume (voxel fill): ${fmtVol(vol)}\n` +
      `Principal dims: ${dims.map(fmtLen).join(" x ")}\n` +
      `Box (X x Y x Z): ${ext.map(fmtLen).join(" x ")}\n` +
      `Voxel: ${fmtLen(h)}`);
  }

  // ================================================================ OUTPUT
  function overlayObjects() {
    return [typeof gridHelper !== "undefined" ? gridHelper : null,
      typeof camerasGroup !== "undefined" ? camerasGroup : null,
      typeof trackballGizmo !== "undefined" ? trackballGizmo : null,
      typeof rulerGroup !== "undefined" ? rulerGroup : null,
      st.measGroup, st.clip.helper,
      typeof clipBoxGroup !== "undefined" ? clipBoxGroup : null].filter(Boolean);
  }
  function hideOverlays() {
    if (st.hidden) return;
    st.hidden = overlayObjects().map((o) => [o, o.visible]);
    st.hidden.forEach(([o]) => { o.visible = false; });
    if (typeof transformControls !== "undefined" && transformControls) { st.hiddenTC = transformControls.visible; transformControls.visible = false; }
  }
  function restoreOverlays() {
    if (!st.hidden) return;
    st.hidden.forEach(([o, v]) => { o.visible = v; });
    st.hidden = null;
    if (typeof transformControls !== "undefined" && transformControls && st.hiddenTC !== undefined) transformControls.visible = st.hiddenTC;
  }

  // camera pose of training camera c (row-major world->camera R, centred position)
  function poseFromRig(c) {
    const R = c.R;
    const X = new THREE.Vector3(R[0], R[1], R[2]);
    const Y = new THREE.Vector3(-R[3], -R[4], -R[5]);
    const Z = new THREE.Vector3(-R[6], -R[7], -R[8]);
    const q = new THREE.Quaternion().setFromRotationMatrix(new THREE.Matrix4().makeBasis(X, Y, Z));
    modelContainer.updateMatrixWorld();
    const qm = new THREE.Quaternion();
    modelContainer.getWorldQuaternion(qm);
    const pos = modelContainer.localToWorld(new THREE.Vector3(c.pos[0], c.pos[1], c.pos[2]));
    const fwdLocal = new THREE.Vector3(R[6], R[7], R[8]);
    return {
      pos: pos,
      quat: qm.clone().multiply(q),
      fov: 2 * Math.atan((c.h || 1080) / (2 * (c.fy || c.fx || 1000))) * 180 / Math.PI,
      fwd: fwdLocal.applyQuaternion(qm)
    };
  }

  function saveCamState() {
    return {
      pos: camera.position.clone(), quat: camera.quaternion.clone(), up: camera.up.clone(),
      target: controls.target.clone(), fov: camera.fov
    };
  }
  function restoreCamState(s) {
    if (!s) return;
    camera.position.copy(s.pos); camera.quaternion.copy(s.quat); camera.up.copy(s.up);
    camera.fov = s.fov; camera.updateProjectionMatrix();
    controls.target.copy(s.target);
    controls.enabled = true;
    controls.update();
  }
  function applyPose(p) {
    camera.position.copy(p.pos);
    camera.quaternion.copy(p.quat);
    if (Math.abs(camera.fov - p.fov) > 1e-3) { camera.fov = p.fov; camera.updateProjectionMatrix(); }
    controls.target.copy(p.pos.clone().addScaledVector(p.fwd, Math.max(0.01, targetDim * 0.8)));
  }

  function ensurePerspective() {
    if (isOrtho && typeof toggleCameraProjection === "function") toggleCameraProjection();
  }

  // ---- animations (turntable / fly-through)
  function startAnim(kind, record) {
    if (st.anim) { stopAnim(); return; }
    if (st.cmp.on) exitCompare();
    ensurePerspective();
    const dur = Math.max(2, parseFloat(($("gss-anim-dur") || {}).value || "12"));
    const saved = saveCamState();
    let apply;
    if (kind !== "fly") {
      apply = orbitPath(kind);
    } else {
      if (!rawCameras || rawCameras.length < 2) { toast("Fly-through needs camera poses"); return; }
      const step = Math.max(1, parseInt(($("gss-fly-step") || {}).value || "4", 10));
      const cams = rawCameras.slice().sort((a, b) => String(a.label).localeCompare(String(b.label), undefined, { numeric: true }));
      const keys = [];
      for (let i = 0; i < cams.length; i += step) keys.push(poseFromRig(cams[i]));
      if (keys.length < 2) { toast("Not enough cameras for a fly-through"); return; }
      for (let i = 1; i < keys.length; i++) if (keys[i].quat.dot(keys[i - 1].quat) < 0) { const q = keys[i].quat; q.set(-q.x, -q.y, -q.z, -q.w); }
      const curve = new THREE.CatmullRomCurve3(keys.map((k) => k.pos), false, "centripetal");
      const fov = keys.reduce((s, k) => s + k.fov, 0) / keys.length;
      apply = (t) => {
        const u = t * (keys.length - 1);
        const i0 = Math.min(keys.length - 2, Math.floor(u));
        const f = u - i0;
        const fs = f * f * (3 - 2 * f);
        const pos = curve.getPoint(t);
        const quat = keys[i0].quat.clone().slerp(keys[i0 + 1].quat, fs);
        const fwd = new THREE.Vector3(0, 0, -1).applyQuaternion(quat);
        applyPose({ pos: pos, quat: quat, fov: fov, fwd: fwd });
      };
    }
    if (($("gss-clean") || {}).checked) hideOverlays();
    controls.enabled = false;
    st.anim = { t0: performance.now() + (record ? 300 : 0), dur: dur, apply: apply, saved: saved, kind: kind };
    if (record) startRecording("gs_" + kind);
    updateAnimButtons();
  }

  // Orbit-type camera paths around the orbit target (up = +z), in spherical
  // coordinates relative to the current view: azimuth, elevation and radius.
  const PATHS = {
    turntable: { turns: 1, height: 0, zoom: 0, note: "One or more full turns at constant height and distance." },
    spiral: { turns: 3, height: 35, zoom: 0, note: "Helix: orbits while the height oscillates by \u00b1Height degrees." },
    vortex: { turns: 2, height: 30, zoom: 60, note: "Spirals in towards the subject while rising (Zoom %), then back out." },
    zoom: { turns: 0, height: 0, zoom: 60, note: "Dolly in and out along the current view (Zoom %); add Turns to orbit too." },
    orbitzoom: { turns: 1, height: 0, zoom: 40, note: "Orbit with two zoom pulses per animation." },
    wiggle: { turns: 2, height: 20, zoom: 0, note: "Wobble: swings \u00b1Height degrees around the subject; Turns = number of swings." },
    fly: { turns: 0, height: 0, zoom: 0, note: "Follows the training cameras (sorted by file name), using every k-th camera." }
  };

  function orbitPath(kind) {
    const num = (id, def) => { const v = parseFloat(($(id) || {}).value); return isFinite(v) ? v : def; };
    const turns = num("gss-turns", 1);
    const A = num("gss-height", 0) * Math.PI / 180;
    const Z = clamp(num("gss-zoom", 0) / 100, 0, 0.95);
    const dir = (($("gss-anim-dir") || {}).value === "cw") ? -1 : 1;
    const ease = !!($("gss-ease") || {}).checked;
    const target = controls.target.clone();
    const off = camera.position.clone().sub(target);
    const r0 = Math.max(1e-6, off.length());
    const az0 = Math.atan2(off.y, off.x);
    const el0 = Math.asin(clamp(off.z / r0, -1, 1));
    const lim = 85 * Math.PI / 180;
    const TAU = 2 * Math.PI;
    return (t) => {
      const u = ease ? 0.5 - 0.5 * Math.cos(Math.PI * t) : t;
      let az = az0 + dir * TAU * turns * u, el = el0, r = r0;
      if (kind === "spiral") el = el0 + A * Math.sin(TAU * u);
      else if (kind === "vortex") { el = el0 + A * 0.5 * (1 - Math.cos(TAU * u)); r = r0 * (1 - Z * 0.5 * (1 - Math.cos(TAU * u))); }
      else if (kind === "zoom") r = r0 * (1 - Z * 0.5 * (1 - Math.cos(TAU * u)));
      else if (kind === "orbitzoom") r = r0 * (1 - Z * 0.5 * (1 - Math.cos(2 * TAU * u)));
      else if (kind === "wiggle") az = az0 + dir * A * Math.sin(TAU * Math.max(turns, 0.5) * u);
      el = clamp(el, -lim, lim);
      camera.position.set(
        target.x + r * Math.cos(el) * Math.cos(az),
        target.y + r * Math.cos(el) * Math.sin(az),
        target.z + r * Math.sin(el));
      camera.up.set(0, 0, 1);
      camera.lookAt(target);
      controls.target.copy(target);
    };
  }

  function setPathPreset(kind) {
    const pr = PATHS[kind] || PATHS.turntable;
    [["gss-turns", pr.turns], ["gss-height", pr.height], ["gss-zoom", pr.zoom]].forEach(([id, v]) => {
      const e = $(id);
      if (!e) return;
      e.value = v;
      e.disabled = kind === "fly";
      const o = $(id + "-v");
      if (o) o.textContent = id === "gss-zoom" ? v + "%" : (id === "gss-height" ? v + "\u00b0" : String(v));
    });
    const fs = $("gss-fly-row");
    if (fs) fs.style.display = kind === "fly" ? "flex" : "none";
    const nt = $("gss-path-note");
    if (nt) nt.textContent = pr.note;
  }

  function stopAnim() {
    const a = st.anim;
    if (!a) return;
    st.anim = null;
    if (st.rec) stopRecording();
    restoreOverlays();
    restoreCamState(a.saved);
    updateAnimButtons();
  }

  function updateAnimButtons() {
    const busy = !!st.anim;
    ["gss-anim-prev", "gss-anim-rec"].forEach((id) => {
      const b = $(id);
      if (!b) return;
      if (busy) { b.textContent = (id.endsWith("rec") && st.rec) ? "\u25a0 Stop recording" : "\u25a0 Stop"; }
      else { b.textContent = id.endsWith("rec") ? "\u25cf Record WebM" : "\u25b6 Preview"; }
    });
  }

  function startRecording(name) {
    if (typeof MediaRecorder === "undefined" || !renderer.domElement.captureStream) { toast("Video recording is not supported by this browser"); return; }
    const fps = parseInt(($("gss-fps") || {}).value || "30", 10);
    const stream = renderer.domElement.captureStream(fps);
    const mime = ["video/webm;codecs=vp9", "video/webm;codecs=vp8", "video/webm"].find((m) => MediaRecorder.isTypeSupported(m)) || "";
    const rec = new MediaRecorder(stream, mime ? { mimeType: mime, videoBitsPerSecond: 16e6 } : undefined);
    const chunks = [];
    rec.ondataavailable = (e) => { if (e.data && e.data.size > 0) chunks.push(e.data); };
    rec.onstop = () => {
      download(new Blob(chunks, { type: "video/webm" }), name + ".webm");
      const b = $("gss-rec-badge"); if (b) b.style.display = "none";
      toast("Video saved (" + name + ".webm)");
    };
    rec.start(250);
    st.rec = rec;
    const b = $("gss-rec-badge"); if (b) b.style.display = "block";
  }

  function stopRecording() {
    const r = st.rec;
    st.rec = null;
    if (r && r.state !== "inactive") r.stop();
  }

  // ---- high-resolution / transparent snapshot
  function snapshot() {
    const k = parseFloat(($("gss-snap-scale") || {}).value || "2");
    const transparent = !!($("gss-snap-transp") || {}).checked;
    const clean = !!($("gss-clean") || {}).checked;
    const cam = activeCam();
    const pr0 = renderer.getPixelRatio();
    const maxSz = (renderer.capabilities && renderer.capabilities.maxTextureSize) || 8192;
    const W = window.innerWidth, H = window.innerHeight;
    const pr = Math.min(pr0 * k, maxSz / Math.max(W, H));
    const bg = scene.background;
    if (clean) hideOverlays();
    if (transparent) { scene.background = null; renderer.setClearColor(0x000000, 0); }
    renderer.setPixelRatio(pr);
    if (typeof gsToolsBeforeRender === "function") gsToolsBeforeRender(cam);
    gsSortDirty = true;
    gsUpdateFrame(cam);
    renderer.render(scene, cam);
    const url = renderer.domElement.toDataURL("image/png");
    renderer.setPixelRatio(pr0);
    if (transparent) { scene.background = bg; }
    if (clean) restoreOverlays();
    gsSortDirty = true;
    const a = document.createElement("a");
    a.href = url;
    a.download = `gs_snapshot_${Math.round(W * pr)}x${Math.round(H * pr)}.png`;
    a.click();
    toast(`Snapshot ${Math.round(W * pr)} x ${Math.round(H * pr)} px saved`);
  }

  // ---- compare the render with the original photo
  function photoSrc(c) {
    let src = c.img || "";
    if (!src && typeof localPhotoStore !== "undefined") {
      const lbl = (c.label || "").toLowerCase();
      const base = lbl.split(/[\\/]/).pop();
      src = localPhotoStore[lbl] || localPhotoStore[base] || "";
    }
    return src;
  }

  function layoutCompare() {
    const c = rawCameras[st.cmp.idx];
    const img = $("gss-compare-img");
    if (!c || !img) return;
    const W = window.innerWidth, H = window.innerHeight;
    const s = H / (c.h || 1080);
    img.style.width = (c.w * s) + "px";
    img.style.height = H + "px";
    img.style.left = (W / 2 - (c.cx != null ? c.cx : c.w / 2) * s) + "px";
    img.style.top = (H / 2 - (c.cy != null ? c.cy : c.h / 2) * s) + "px";
    const wrap = $("gss-compare"), line = $("gss-compare-line");
    if (st.cmp.mode === "split") {
      img.style.opacity = 1;
      wrap.style.clipPath = `inset(0 ${100 - st.cmp.split}% 0 0)`;
      line.style.display = "block";
      line.style.left = (st.cmp.split / 100 * W - 1) + "px";
    } else {
      img.style.opacity = st.cmp.opacity;
      wrap.style.clipPath = "none";
      line.style.display = "none";
    }
  }

  function enterCompare(idx) {
    if (!rawCameras || rawCameras.length === 0) { toast("No camera poses available"); return; }
    if (st.anim) stopAnim();
    ensurePerspective();
    const n = rawCameras.length;
    st.cmp.idx = ((idx % n) + n) % n;
    const c = rawCameras[st.cmp.idx];
    if (!st.cmp.on) { st.cmp.saved = saveCamState(); hideOverlays(); }
    st.cmp.on = true;
    st.cmp.pose = poseFromRig(c);
    controls.enabled = false;
    applyPose(st.cmp.pose);
    const img = $("gss-compare-img");
    const src = photoSrc(c);
    img.style.display = src ? "block" : "none";
    if (src) img.src = src;
    $("gss-compare").style.display = "block";
    layoutCompare();
    const lab = $("gss-cmp-label");
    if (lab) lab.textContent = `${st.cmp.idx + 1}/${n}: ${c.label}${src ? "" : " (photo not found)"}`;
    document.querySelectorAll(".gss-cmp-on").forEach((e) => { e.disabled = false; });
  }

  function exitCompare() {
    if (!st.cmp.on) return;
    st.cmp.on = false;
    $("gss-compare").style.display = "none";
    $("gss-compare-line").style.display = "none";
    restoreOverlays();
    restoreCamState(st.cmp.saved);
    const lab = $("gss-cmp-label");
    if (lab) lab.textContent = "";
    document.querySelectorAll(".gss-cmp-on").forEach((e) => { e.disabled = true; });
  }

  // ------------------------------------------------------ per-frame hook
  window.gsToolsBeforeRender = function (cam) {
    if (st.anim) {
      const t = clamp((performance.now() - st.anim.t0) / (st.anim.dur * 1000), 0, 1);
      st.anim.apply(t);
      if (t >= 1) setTimeout(stopAnim, 120);
    } else if (st.cmp.on && st.cmp.pose) {
      applyPose(st.cmp.pose);
    }
  };

  // ============================================================ INTERACTION
  const overlay = document.createElement("canvas");
  overlay.id = "gss-overlay";
  const octx = overlay.getContext("2d");
  function sizeOverlay() { overlay.width = window.innerWidth; overlay.height = window.innerHeight; drawOverlay(); }

  function drawOverlay(mx, my) {
    octx.clearRect(0, 0, overlay.width, overlay.height);
    const d = st.drag;
    octx.lineWidth = 1.5;
    octx.strokeStyle = "#facc15";
    octx.fillStyle = "rgba(250, 204, 21, 0.08)";
    octx.setLineDash([6, 4]);
    if (d && st.tool === "rect") {
      octx.fillRect(d.x0, d.y0, d.x1 - d.x0, d.y1 - d.y0);
      octx.strokeRect(d.x0, d.y0, d.x1 - d.x0, d.y1 - d.y0);
    } else if (d && st.tool === "lasso" && d.poly.length > 1) {
      octx.beginPath();
      octx.moveTo(d.poly[0][0], d.poly[0][1]);
      d.poly.forEach((p) => octx.lineTo(p[0], p[1]));
      octx.closePath(); octx.fill(); octx.stroke();
    }
    if (st.tool === "brush" && mx !== undefined) {
      octx.setLineDash([]);
      octx.beginPath();
      octx.arc(mx, my, st.brushR, 0, 2 * Math.PI);
      octx.stroke();
    }
  }

  function selMode(e) { return e.altKey ? "remove" : (e.shiftKey ? "add" : "replace"); }

  function onDown(e) {
    onRightDown(e);
    if (!st.tool || e.target !== renderer.domElement || e.button !== 0) return;
    e.stopPropagation();
    e.preventDefault();
    const x = e.clientX, y = e.clientY;
    if (st.tool === "rect") st.drag = { x0: x, y0: y, x1: x, y1: y, mode: selMode(e) };
    else if (st.tool === "lasso") st.drag = { poly: [[x, y]], mode: selMode(e) };
    else if (st.tool === "brush") {
      st.drag = { mode: e.altKey ? "remove" : "add" };
      brushAt(x, y, st.drag.mode);
    } else if (st.tool === "calib" || st.tool === "polyline" || st.tool === "angle" || st.tool === "area") {
      const p = window.gsToolsPick(x, y);
      if (!p) { toast("No visible splat under the cursor"); return; }
      if (st.tool === "calib") addCalibPoint(p);
      else addMeasPoint(p);
    }
  }

  function brushAt(x, y, mode) {
    const r2 = st.brushR * st.brushR;
    applySelection((sx, sy) => (sx - x) * (sx - x) + (sy - y) * (sy - y) <= r2, mode);
  }

  let moveQueued = null;
  function onMove(e) {
    if (!st.tool) return;
    const x = e.clientX, y = e.clientY;
    if (st.drag) {
      e.stopPropagation();
      if (st.tool === "rect") { st.drag.x1 = x; st.drag.y1 = y; }
      else if (st.tool === "lasso") st.drag.poly.push([x, y]);
      else if (st.tool === "brush") brushAt(x, y, st.drag.mode);
    }
    if (st.tool === "brush" || st.drag) drawOverlay(x, y);
    if (["polyline", "angle", "area", "calib"].includes(st.tool) && !moveQueued) {
      moveQueued = requestAnimationFrame(() => {
        moveQueued = null;
        const p = window.gsToolsPick(x, y);
        if (!st.hover) { st.hover = markerMesh(new THREE.Vector3(), 0xffffff); st.hover.material.wireframe = true; measGroup().add(st.hover); }
        st.hover.visible = !!p;
        if (p) st.hover.position.copy(p);
      });
    }
  }

  let rDown = null;
  function onRightDown(e) {
    if (e.button === 2 && e.target === renderer.domElement) rDown = [e.clientX, e.clientY];
  }
  function onRightUp(e) {
    if (e.button !== 2 || !rDown) return;
    const moved = Math.hypot(e.clientX - rDown[0], e.clientY - rDown[1]);
    rDown = null;
    if (moved < 5 && st.cur && (st.cur.type === "polyline" || st.cur.type === "area")) finishMeasurement();
  }

  function onUp(e) {
    onRightUp(e);
    if (!st.tool) return;
    if (e.target === renderer.domElement && e.button === 0) e.stopPropagation();
    const d = st.drag;
    if (!d) return;
    st.drag = null;
    if (st.tool === "rect") {
      const x0 = Math.min(d.x0, d.x1), x1 = Math.max(d.x0, d.x1), y0 = Math.min(d.y0, d.y1), y1 = Math.max(d.y0, d.y1);
      if (x1 - x0 > 2 && y1 - y0 > 2) applySelection((sx, sy) => sx >= x0 && sx <= x1 && sy >= y0 && sy <= y1, d.mode);
    } else if (st.tool === "lasso" && d.poly.length > 2) {
      let bx0 = Infinity, bx1 = -Infinity, by0 = Infinity, by1 = -Infinity;
      d.poly.forEach(([px, py]) => { bx0 = Math.min(bx0, px); bx1 = Math.max(bx1, px); by0 = Math.min(by0, py); by1 = Math.max(by1, py); });
      const poly = d.poly;
      applySelection((sx, sy) => sx >= bx0 && sx <= bx1 && sy >= by0 && sy <= by1 && pointInPoly(sx, sy, poly), d.mode);
    }
    drawOverlay(e.clientX, e.clientY);
  }

  function onDbl(e) {
    if (!st.tool || e.target !== renderer.domElement) return;
    e.stopPropagation();
    if ((st.tool === "polyline" || st.tool === "area") && st.cur) finishMeasurement();
  }

  function onKey(e) {
    if (e.target && (e.target.tagName === "INPUT" || e.target.tagName === "TEXTAREA" || e.target.tagName === "SELECT")) return;
    if (e.key === "Escape") {
      if (st.anim) { stopAnim(); return; }
      if (st.cmp.on) { exitCompare(); return; }
      if (st.cur) { cancelMeasurement(); return; }
      if (st.tool) { setTool(st.tool); return; }
    }
    if (e.key === "Enter" && st.cur) { finishMeasurement(); return; }
    if (!panel.classList.contains("open")) return;
    if ((e.key === "Delete" || e.key === "Backspace") && st.selCount > 0) { e.preventDefault(); opDelete(); }
    if ((e.ctrlKey || e.metaKey) && (e.key === "z" || e.key === "Z")) { e.preventDefault(); opUndo(); }
    if (e.key === "[" ) { st.brushR = Math.max(4, st.brushR - 4); syncBrush(); }
    if (e.key === "]") { st.brushR = Math.min(200, st.brushR + 4); syncBrush(); }
  }

  function syncBrush() {
    const s = $("gss-brush"); if (s) s.value = st.brushR;
    const v = $("gss-brush-v"); if (v) v.textContent = st.brushR + " px";
    drawOverlay();
  }

  // ===================================================================== UI
  const panel = document.createElement("div");
  panel.id = "gss-panel";

  function slider(id, label, min, max, step, value, fmt) {
    return `<div class="gss-row"><label for="${id}">${label}</label><input type="range" id="${id}" min="${min}" max="${max}" step="${step}" value="${value}"><span class="gss-num" id="${id}-v">${fmt}</span></div>`;
  }

  function buildPanel() {
    const shAvail = gsSHDegree > 0 && gsSHData;
    const hasCams = rawCameras && rawCameras.length > 0;
    panel.innerHTML = `
      <div class="gss-head"><span>&#129520; GS Studio <small>Gaussian splat tools</small></span>
        <button class="gss-btn" id="gss-close" title="Close">&#10005;</button></div>
      <div class="gss-tabs">
        <button class="gss-tab active" data-gss-tab="render">Render</button>
        <button class="gss-tab" data-gss-tab="edit">Edit</button>
        <button class="gss-tab" data-gss-tab="measure">Measure</button>
        <button class="gss-tab" data-gss-tab="output">Output</button>
      </div>
      <div class="gss-body">

      <div class="gss-page active" data-gss-page="render">
        <div class="gss-sec"><h4>View-dependent color <span class="gss-val">${shAvail ? "SH degree " + gsSHDegree : "n/a"}</span></h4>
          <label class="gss-check"><input type="checkbox" id="gss-sh" ${shAvail ? "checked" : "disabled"}> Spherical harmonics (reflections and highlights change with the viewing angle)</label>
          <div class="gss-note">${shAvail ? "Applies to the RGB color mode in Splats mode." : "This model has base colors only (untrained cloud or SH degree 0)."}</div>
        </div>
        <div class="gss-sec"><h4>Floater filter <span class="gss-val" id="gss-floater-stat">off</span></h4>
          ${slider("gss-minopac", "Min opacity", 0, 0.95, 0.01, 0, "0.00")}
          ${slider("gss-maxsize", "Max size", 50, 100, 0.5, 100, "off")}
          <div class="gss-note">Hides translucent splats and the largest ones (size percentile), the haze that surrounds trained models. Non-destructive; exports respect it.</div>
        </div>
        <div class="gss-sec"><h4>Image adjustments <button class="gss-btn" id="gss-adj-reset">Reset</button></h4>
          ${slider("gss-adj-exposure", "Exposure", -3, 3, 0.05, 0, "+0.00")}
          ${slider("gss-adj-contrast", "Contrast", 0.5, 2, 0.01, 1, "1.00")}
          ${slider("gss-adj-saturation", "Saturation", 0, 2, 0.01, 1, "1.00")}
          ${slider("gss-adj-gamma", "Gamma", 0.5, 2.5, 0.01, 1, "1.00")}
        </div>
        <div class="gss-sec"><h4>Background</h4>
          <div class="gss-row">
            <button class="gss-btn" data-gss-bg="#090d16">Dark</button>
            <button class="gss-btn" data-gss-bg="#000000">Black</button>
            <button class="gss-btn" data-gss-bg="#ffffff">White</button>
            <button class="gss-btn" data-gss-bg="#808080">Gray</button>
            <button class="gss-btn" data-gss-bg="#00b140">Chroma</button>
            <input type="color" id="gss-bg-color" value="#090d16" title="Custom color">
          </div>
        </div>
      </div>

      <div class="gss-page" data-gss-page="edit">
        <div class="gss-sec"><h4>Select</h4>
          <div class="gss-grid3">
            <button class="gss-btn" data-gss-tool="rect">&#9634; Rectangle</button>
            <button class="gss-btn" data-gss-tool="lasso">&#10227; Lasso</button>
            <button class="gss-btn" data-gss-tool="brush">&#9679; Brush</button>
          </div>
          ${slider("gss-brush", "Brush size", 4, 200, 1, st.brushR, st.brushR + " px")}
          <div class="gss-grid4">
            <button class="gss-btn" id="gss-sel-all">All</button>
            <button class="gss-btn" id="gss-sel-inv">Invert</button>
            <button class="gss-btn" id="gss-sel-clear">None</button>
            <button class="gss-btn" id="gss-undo" title="Ctrl+Z">Undo</button>
          </div>
          <div class="gss-grid2">
            <button class="gss-btn danger" id="gss-del" title="Delete key">&#128465; Delete selected</button>
            <button class="gss-btn" id="gss-iso">&#9673; Keep selected only</button>
          </div>
          <div class="gss-stat" id="gss-edit-stat"></div>
          <div class="gss-status" id="gss-hint"></div>
          <div class="gss-note">Selections go through the model (all depths). Right button pans and the wheel zooms while a tool is active; Esc leaves the tool. [ and ] resize the brush.</div>
        </div>
        <div class="gss-sec"><h4>Section plane <span class="gss-val" id="gss-clip-stat">off</span></h4>
          <label class="gss-check"><input type="checkbox" id="gss-clip-on"> Cut the model with a plane</label>
          <div class="gss-grid4">
            <button class="gss-btn" data-gss-axis="x">X</button>
            <button class="gss-btn" data-gss-axis="y">Y</button>
            <button class="gss-btn active" data-gss-axis="z">Z</button>
            <button class="gss-btn" data-gss-axis="view" title="Plane facing the current view">View</button>
          </div>
          ${slider("gss-clip-pos", "Position", 0, 100, 0.1, 50, "50%")}
          <div class="gss-row">
            <label class="gss-check"><input type="checkbox" id="gss-clip-flip"> Flip side</label>
            <label class="gss-check"><input type="checkbox" id="gss-clip-helper" checked> Show plane</label>
          </div>
        </div>
        <div class="gss-sec"><h4>Export</h4>
          <button class="gss-btn primary" id="gss-export">&#128190; Export visible splats (3DGS .ply)</button>
          <div class="gss-note">Writes every visible splat (after deletions, floater filter, section plane, crop and noise filter) with SH, opacity, scales and rotations, readable by gs_import_metashape() and other 3DGS viewers. Colors and SH are stored at viewer precision (8-bit SH).</div>
        </div>
      </div>

      <div class="gss-page" data-gss-page="measure">
        <div class="gss-sec"><h4>Scale</h4>
          <div class="gss-note" id="gss-unit-status"></div>
          <button class="gss-btn" data-gss-tool="calib">&#128207; Calibrate: pick 2 points</button>
          <div class="gss-row"><label>Real length</label><input type="number" id="gss-calib-len" min="0" step="any" placeholder="e.g. 10">
            <select id="gss-calib-unit"><option value="mm">mm</option><option value="cm" selected>cm</option><option value="m">m</option></select>
            <button class="gss-btn primary" id="gss-calib-apply">Apply</button></div>
          <div class="gss-row"><label>or factor</label><input type="number" id="gss-unit-factor" min="0" step="any" placeholder="m per unit">
            <button class="gss-btn" id="gss-unit-set">Set</button><button class="gss-btn" id="gss-unit-reset">Unknown</button></div>
          <div class="gss-status" id="gss-calib-status"></div>
        </div>
        <div class="gss-sec"><h4>Measure</h4>
          <div class="gss-grid3">
            <button class="gss-btn" data-gss-tool="polyline">&#12316; Polyline</button>
            <button class="gss-btn" data-gss-tool="angle">&ang; Angle</button>
            <button class="gss-btn" data-gss-tool="area">&#9649; Area</button>
          </div>
          <div class="gss-row" id="gss-meas-live" style="display:none; background: rgba(250,204,21,0.1); border: 1px solid rgba(250,204,21,0.4); border-radius: 7px; padding: 5px 7px;">
            <span id="gss-meas-live-txt" class="gss-stat" style="flex:1; color:#facc15;"></span>
            <button class="gss-btn primary" id="gss-meas-finish">&#10003; Finish</button>
            <button class="gss-btn danger" id="gss-meas-cancel">Cancel</button>
          </div>
          <div class="gss-list" id="gss-meas-list"></div>
          <div class="gss-grid2"><button class="gss-btn" id="gss-meas-copy">Copy as CSV</button><button class="gss-btn danger" id="gss-meas-clear">Clear all</button></div>
          <div class="gss-note">Points snap to the front-most visible splat. The toolbar 3D Ruler uses the same units.</div>
        </div>
        <div class="gss-sec"><h4>Volume &amp; dimensions</h4>
          ${slider("gss-vol-res", "Resolution", 32, 192, 8, 96, "96")}
          <button class="gss-btn primary" id="gss-vol-run">Compute (selection, or all visible)</button>
          <div class="gss-stat" id="gss-vol-out"></div>
          <div class="gss-note">Voxel-filled volume of the splat centres (closed shapes), principal (PCA) and box dimensions. Select or crop the subject first to exclude the background.</div>
        </div>
      </div>

      <div class="gss-page" data-gss-page="output">
        <div class="gss-sec"><h4>Animation &amp; video</h4>
          <div class="gss-row"><label>Path</label><select id="gss-path" style="flex:1">
            <option value="turntable">Turntable 360&deg;</option>
            <option value="spiral">Spiral (helix)</option>
            <option value="vortex">Vortex (spiral zoom-in)</option>
            <option value="zoom">Zoom in-out</option>
            <option value="orbitzoom">Orbit + zoom pulse</option>
            <option value="wiggle">Wiggle (3D wobble)</option>
            <option value="fly" ${hasCams ? "" : "disabled"}>Fly-through (training cameras)</option>
          </select></div>
          ${slider("gss-turns", "Turns", 0, 6, 0.25, 1, "1")}
          ${slider("gss-height", "Height", 0, 80, 1, 0, "0&deg;")}
          ${slider("gss-zoom", "Zoom", 0, 90, 1, 0, "0%")}
          <div class="gss-row" id="gss-fly-row" style="display:none"><label>Every k-th cam</label><input type="number" id="gss-fly-step" value="4" min="1" step="1"></div>
          <div class="gss-row"><label>Duration (s)</label><input type="number" id="gss-anim-dur" value="12" min="2" step="1">
            <select id="gss-anim-dir"><option value="ccw">CCW</option><option value="cw">CW</option></select>
            <label class="gss-check"><input type="checkbox" id="gss-ease" checked> Ease</label></div>
          <div class="gss-row"><button class="gss-btn" id="gss-anim-prev" style="flex:1">&#9654; Preview</button><button class="gss-btn rec" id="gss-anim-rec" style="flex:1">&#9679; Record WebM</button>
            <label>FPS</label><select id="gss-fps"><option>24</option><option selected>30</option><option>60</option></select></div>
          <div class="gss-note" id="gss-path-note"></div>
          <label class="gss-check"><input type="checkbox" id="gss-clean" checked> Clean render (hide grid, cameras, gizmos, measurements)</label>
          <div class="gss-note">Paths start from the current view and orbit the rotation centre (double-click the model to set it). Video size = window size. Esc stops.</div>
        </div>
        <div class="gss-sec"><h4>Snapshot</h4>
          <div class="gss-row"><label>Scale</label><select id="gss-snap-scale"><option value="1">1x</option><option value="2" selected>2x</option><option value="3">3x</option><option value="4">4x</option></select>
            <label class="gss-check"><input type="checkbox" id="gss-snap-transp"> Transparent</label></div>
          <button class="gss-btn primary" id="gss-snap">&#128248; Save PNG</button>
        </div>
        <div class="gss-sec"><h4>Compare with photo</h4>
          <div class="gss-row"><button class="gss-btn primary" id="gss-cmp-start" ${hasCams ? "" : "disabled"}>Align to camera</button>
            <input type="number" id="gss-cmp-idx" value="1" min="1" step="1" ${hasCams ? "" : "disabled"}>
            <button class="gss-btn gss-cmp-on" id="gss-cmp-prev" disabled>&#9664;</button>
            <button class="gss-btn gss-cmp-on" id="gss-cmp-next" disabled>&#9654;</button>
            <button class="gss-btn gss-cmp-on danger" id="gss-cmp-exit" disabled>Exit</button></div>
          <div class="gss-row"><label>Mode</label><select id="gss-cmp-mode"><option value="blend">Blend</option><option value="split">Split</option></select></div>
          ${slider("gss-cmp-op", "Photo opacity", 0, 1, 0.01, 0.5, "0.50")}
          ${slider("gss-cmp-split", "Split", 0, 100, 0.5, 50, "50%")}
          <div class="gss-status" id="gss-cmp-label"></div>
          <div class="gss-note">Looks through a training camera (its focal length and principal point) and overlays the original photo to check the fidelity of the reconstruction.</div>
        </div>
      </div>
      </div>`;
  }

  function bind(id, ev, fn) { const e = $(id); if (e) e.addEventListener(ev, fn); }

  function bindPanel() {
    $("gss-close").addEventListener("click", togglePanel);
    panel.querySelectorAll("[data-gss-tab]").forEach((b) => b.addEventListener("click", () => {
      const t = b.getAttribute("data-gss-tab");
      panel.querySelectorAll("[data-gss-tab]").forEach((x) => x.classList.toggle("active", x === b));
      panel.querySelectorAll("[data-gss-page]").forEach((p) => p.classList.toggle("active", p.getAttribute("data-gss-page") === t));
      if (t === "edit") updateEditStats();
    }));
    // render
    bind("gss-sh", "change", (e) => { gsSHEnabled = e.target.checked; });
    bind("gss-minopac", "input", (e) => { st.minOpac = +e.target.value; $("gss-minopac-v").textContent = st.minOpac.toFixed(2); rebuildMask(); });
    bind("gss-maxsize", "input", (e) => { st.sizePct = +e.target.value; $("gss-maxsize-v").textContent = st.sizePct >= 100 ? "off" : st.sizePct.toFixed(1) + "%"; rebuildMask(); });
    ["exposure", "contrast", "saturation", "gamma"].forEach((k) => bind("gss-adj-" + k, "input", (e) => setAdjust(k, +e.target.value)));
    bind("gss-adj-reset", "click", () => {
      [["exposure", 0], ["contrast", 1], ["saturation", 1], ["gamma", 1]].forEach(([k, v]) => { $("gss-adj-" + k).value = v; setAdjust(k, v); });
    });
    panel.querySelectorAll("[data-gss-bg]").forEach((b) => b.addEventListener("click", () => setBackground(b.getAttribute("data-gss-bg"))));
    bind("gss-bg-color", "input", (e) => setBackground(e.target.value));
    // edit
    panel.querySelectorAll("[data-gss-tool]").forEach((b) => b.addEventListener("click", () => setTool(b.getAttribute("data-gss-tool"))));
    bind("gss-brush", "input", (e) => { st.brushR = +e.target.value; syncBrush(); });
    bind("gss-sel-all", "click", opSelectAll);
    bind("gss-sel-inv", "click", opInvert);
    bind("gss-sel-clear", "click", opClearSel);
    bind("gss-undo", "click", opUndo);
    bind("gss-del", "click", opDelete);
    bind("gss-iso", "click", opIsolate);
    bind("gss-clip-on", "change", (e) => { st.clip.on = e.target.checked; if (st.clip.on) clipRange(); rebuildMask(); });
    panel.querySelectorAll("[data-gss-axis]").forEach((b) => b.addEventListener("click", () => setClipAxis(b.getAttribute("data-gss-axis"))));
    bind("gss-clip-pos", "input", (e) => { st.clip.pos = +e.target.value; $("gss-clip-pos-v").textContent = st.clip.pos.toFixed(1) + "%"; if (st.clip.on) rebuildMask(); });
    bind("gss-clip-flip", "change", (e) => { st.clip.flip = e.target.checked; if (st.clip.on) rebuildMask(); });
    bind("gss-clip-helper", "change", (e) => { st.clip.showHelper = e.target.checked; updateClipHelper(); });
    bind("gss-export", "click", exportPLY3DGS);
    // measure
    bind("gss-calib-apply", "click", applyCalibration);
    bind("gss-unit-set", "click", () => {
      const v = parseFloat($("gss-unit-factor").value);
      if (v > 0) { gsUnitScale = v; refreshUnits(); toast("Scale factor set"); }
    });
    bind("gss-unit-reset", "click", () => { gsUnitScale = 0; refreshUnits(); });
    bind("gss-meas-copy", "click", copyMeasCSV);
    bind("gss-meas-finish", "click", finishMeasurement);
    bind("gss-meas-cancel", "click", cancelMeasurement);
    bind("gss-meas-clear", "click", () => { while (st.meas.length) removeMeas(0); });
    bind("gss-vol-res", "input", (e) => { $("gss-vol-res-v").textContent = e.target.value; });
    bind("gss-vol-run", "click", computeVolume);
    // output
    bind("gss-path", "change", (e) => setPathPreset(e.target.value));
    bind("gss-turns", "input", (e) => { $("gss-turns-v").textContent = e.target.value; });
    bind("gss-height", "input", (e) => { $("gss-height-v").textContent = e.target.value + "\u00b0"; });
    bind("gss-zoom", "input", (e) => { $("gss-zoom-v").textContent = e.target.value + "%"; });
    bind("gss-anim-prev", "click", () => startAnim(($("gss-path") || {}).value || "turntable", false));
    bind("gss-anim-rec", "click", () => startAnim(($("gss-path") || {}).value || "turntable", true));
    setPathPreset("turntable");
    bind("gss-snap", "click", snapshot);
    bind("gss-cmp-start", "click", () => enterCompare((parseInt($("gss-cmp-idx").value, 10) || 1) - 1));
    bind("gss-cmp-prev", "click", () => { enterCompare(st.cmp.idx - 1); $("gss-cmp-idx").value = st.cmp.idx + 1; });
    bind("gss-cmp-next", "click", () => { enterCompare(st.cmp.idx + 1); $("gss-cmp-idx").value = st.cmp.idx + 1; });
    bind("gss-cmp-exit", "click", exitCompare);
    bind("gss-cmp-mode", "change", (e) => { st.cmp.mode = e.target.value; layoutCompare(); });
    bind("gss-cmp-op", "input", (e) => { st.cmp.opacity = +e.target.value; $("gss-cmp-op-v").textContent = st.cmp.opacity.toFixed(2); layoutCompare(); });
    bind("gss-cmp-split", "input", (e) => { st.cmp.split = +e.target.value; $("gss-cmp-split-v").textContent = st.cmp.split.toFixed(0) + "%"; layoutCompare(); });
  }

  function applyDefaultMinOpacity() {
    const thr = (typeof gsMinOpacityDefault === "number") ? gsMinOpacityDefault : 0;
    const n = N();
    if (!(thr > 0) || n === 0 || !rawOpacities || rawOpacities.length !== n) return;
    let hidden = 0;
    for (let i = 0; i < n; i++) if (rawOpacities[i] < thr) hidden++;
    if (hidden > 0.9 * n) {
      toast(`Min opacity ${thr} would hide ${Math.round(100 * hidden / n)}% of the splats: floater filter left off`);
      return;
    }
    st.minOpac = thr;
    const sl = $("gss-minopac");
    if (sl) sl.value = thr;
    const v = $("gss-minopac-v");
    if (v) v.textContent = thr.toFixed(2);
    rebuildMask();
  }

  function togglePanel() {
    panel.classList.toggle("open");
    const b = $("tb-gss");
    if (b) b.classList.toggle("active", panel.classList.contains("open"));
    if (!panel.classList.contains("open") && st.tool) setTool(st.tool);
  }

  function init() {
    if (typeof renderer === "undefined" || !renderer || typeof modelContainer === "undefined" || !modelContainer) {
      setTimeout(init, 100);
      return;
    }
    document.body.appendChild(overlay);
    const cmp = document.createElement("div");
    cmp.id = "gss-compare";
    cmp.innerHTML = `<img id="gss-compare-img" alt="">`;
    document.body.appendChild(cmp);
    const line = document.createElement("div");
    line.id = "gss-compare-line";
    document.body.appendChild(line);
    const badge = document.createElement("div");
    badge.id = "gss-rec-badge";
    badge.textContent = "\u25cf REC";
    document.body.appendChild(badge);

    buildPanel();
    document.body.appendChild(panel);
    bindPanel();

    const tb = document.querySelector(".metashape-toolbar");
    if (tb) {
      const div = document.createElement("div");
      div.className = "tb-divider";
      const btn = document.createElement("button");
      btn.id = "tb-gss";
      btn.className = "tb-btn";
      btn.title = "GS Studio: Gaussian splat render, edit, measure and output tools";
      btn.innerHTML = "&#129520; GS Studio";
      btn.style.color = "#fbbf24";
      btn.style.fontWeight = "700";
      btn.addEventListener("click", togglePanel);
      tb.appendChild(div);
      tb.appendChild(btn);
    }

    window.addEventListener("pointerdown", onDown, true);
    window.addEventListener("pointermove", onMove, true);
    window.addEventListener("pointerup", onUp, true);
    window.addEventListener("dblclick", onDbl, true);
    window.addEventListener("keydown", onKey, true);
    window.addEventListener("resize", () => { sizeOverlay(); if (st.cmp.on) layoutCompare(); });
    sizeOverlay();
    applyDefaultMinOpacity();
    refreshUnits();
    renderMeasList();
    updateEditStats();
    // open by default for optimized Gaussian splats
    if (rawCov3D && rawCov3D.length === N() * 6 && N() > 0) togglePanel();
  }

  window.GSStudio = { state: st, rebuildMask: rebuildMask, exportPLY: exportPLY3DGS, snapshot: snapshot, enterCompare: enterCompare };
  if (document.readyState === "complete") setTimeout(init, 0);
  else window.addEventListener("load", () => setTimeout(init, 0));
})();
