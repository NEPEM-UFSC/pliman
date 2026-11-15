# Plano de Otimização — Inferência de Imagem/Vídeo e PerSAM / Grounded-SAM

> Branch analisada: `machine_vision` (commit `9521017`). Escopo: motor ONNX em C++
> (`src/onnx_inference.cpp`, `src/image_resize.cpp`), runners R (`R/utils_neural.R`,
> `R/utils_video.R`, `R/utils_dl_batch.R`, `R/utils_dl_fastpath.R`, `R/utils_yolo_train.R`).
> Os ganhos citados são **estimativas** a confirmar com o protocolo de medição da Seção 1.

---

## STATUS DE IMPLEMENTAÇÃO

Atualizado após a execução das Fases 0–1. Hardware das medições: Intel i7 (12 threads),
NVIDIA GTX 1650 (3935 MB VRAM, DirectML), Windows 11, ONNX Runtime 1.20.1, R 4.5.1.

### Implementado e validado

**Bugs funcionais (todos com teste de regressão em `tests/testthat/test-onnx-inference-fixes.R`)**

| # | Bug | Correção |
|---|-----|----------|
| 1 | Hash amostrado do cache do encoder | `hash_bytes64()` sobre o tensor inteiro; cache virou **LRU de K entradas** (`PLIMAN_SAM_CACHE`, padrão 4) |
| 2 | Modo denso sem kernels → 0 objetos | Kernels ZNCC serializados em `prototypes$dense` e recarregados; se ainda faltarem, cai para `feat_res = 256` com aviso |
| 3 | `resize(8)` descartava os cliques novos | Protótipos do clique atual entram **antes** dos acumulados |
| 4 | `omp_set_num_threads()` global | Removido; só a cláusula `num_threads(...)` do pragma |
| 5 | Prompt de caixa virava 2 cliques | Cantos (labels 2/3) enviados como **um** prompt `[1, 2, 2]`; pontos label 1 mantêm o laço |
| 6 | `text_threshold` ignorado; dims 900×256 fixas | Aplicado como piso no token vencedor; dims lidas do tensor de saída |
| 8 | OpenMP aninhado em `run_yolo_cpp` | Laço interno condicionado a `!par_inst` |
| 9 | Inconsistência de dispositivo do SAM | **Resolvida por medição** (ver abaixo) |

Bug 7 (distorção de aspecto por `letterbox = FALSE`) **não foi corrigido de propósito**:
alteraria o resultado de todos os usuários silenciosamente. Deve entrar junto com a
reescrita do pré-processamento, como opt-in.

**Desempenho**

- Evicção seletiva de sessão: novo `onnx_evict_session_cpp(model_path)`; o Grounded-SAM
  deixou de chamar `clear_onnx_sessions_cpp()` (item 1 da tabela de prioridades).
  Relevante nesta máquina: 3935 MB < 6000 MB, logo o flush disparava **a cada imagem**.
- Cópias de tensor: `as.numeric()` removido dos 7 pontos de chamada SAM/DINO
  (25 MB por chamada em 1024²).
- `image_data(curr_img, type = "numeric")` → `image_data(curr_img)` no autolabel:
  array `raw` em vez de `integer`, **4× menos memória**, tensor bit-idêntico (verificado).
- Projeção de rótulos do PerSAM com LUT de linha/coluna.
- NMS espacial dos picos: grade de células → **O(N)** no lugar de O(N_cand × N_aceitos).
- `pixel_mask` do DINO (5 MB int64) em buffer estático.
- `run_onnx_inference_cpp` usa os metadados da sessão (sem alocar nomes por chamada).
- Sigmoide do DINO em logit-space, com corte antecipado por query.

**PerSAM em GPU (opt-in)** — novo argumento `sam_device = c("cpu", "gpu", "auto")` em
`image_segment_dl()`, `.run_persam()` e `.compute_dl_mask()`. **Padrão `"cpu"`,
comportamento anterior preservado byte a byte.** Em GPU o decoder roda serialmente
(uma submissão por vez) em lugar do laço OpenMP.

**Novos utilitários:** `onnx_evict_session_cpp()`, `onnx_clear_encoder_cache_cpp()`,
`pliman_omp_max_threads_cpp()`.

### Medições obtidas

| Teste | Resultado |
|---|---|
| Encoder SAM 2.1 em DirectML, 40 iterações encode+decode | **40/40 estáveis**, 0,31 s/iteração, 2/2 sessões na GPU, sem fallback |
| PerSAM `sam_device="cpu"` vs `"gpu"` (mesma imagem/parâmetros) | 7,22 s → **2,58 s (2,8×)**, mesmos 4 objetos |
| Cache do encoder: mesmo tensor | 7,93 s (frio) → 1,88 s (hit) |
| Cache do encoder: tensor alterado | 4,13 s — **miss correto** (antes era hit falso) |
| Grounded-SAM, 3 chamadas seguidas (CPU) | 16,6 s → 11,2 s em regime; 3 sessões permanecem quentes |
| Grounding DINO `text_threshold` 0,05 → 0,95 | 4 → 0 caixas (antes: sem efeito) |

**Conclusão sobre o bug 9:** o travamento do DirectML alegado para justificar CPU
forçada **não se reproduz** nesta configuração. Como a falha original presumivelmente
ocorreu em algum driver/runtime, a GPU entrou como opt-in e não como padrão.

---

## STATUS — FASE 2

### Aviso sobre as medições da Fase 1

As medições registradas acima na seção da Fase 1 foram obtidas com
`pkgbuild::compile_dll()` no padrão (`debug = TRUE`), que compila com **`-O0`**. Os
arquivos alterados ficaram sem otimização enquanto o resto do pacote estava em `-O2`,
o que **inflou sistematicamente** os tempos dos laços em C++. Todos os números da Fase
2 abaixo vêm de build otimizado (`debug = FALSE`, `-O2`), e as conclusões qualitativas
da Fase 1 (CPU vs GPU, cache correto, evicção seletiva) se mantiveram ao refazer.

Para medir desempenho neste pacote: `pkgbuild::compile_dll(".", debug = FALSE)`.

### Implementado e validado

**§3.2-2 — Cache de encoder compartilhado.** Novo helper `sam2_encode_cached()` usado
pelos três runners (`run_sam2_inference_cpp`, `run_sam2_instances_cpp`,
`run_sam2_persam_cpp`). Antes só o PerSAM tinha cache; o Grounded-SAM reexecutava o
backbone a cada chamada.

| Cenário | Antes | Depois |
|---|---|---|
| Grounded-SAM, 2º prompt na mesma imagem | 17,54 s | **9,6–11,0 s** (1,70×) |
| PerSAM após Grounded-SAM na mesma imagem | 8,40 s | **5,78 s** (1,45×) |

**§2.3 — Decoder do SAM em lote.** O decoder exportado aceita lote de prompts
(`num_labels` dinâmico), então B prompts vão em um `Run`. Medido em imagem 2448²
com 176 picos, **instâncias idênticas em todos os tamanhos de lote**:

| | lote 1 | lote 8 | lote 16 | lote 32 | lote 64 |
|---|---|---|---|---|---|
| GPU | 3,49 s | 2,76 s | 2,07 s | **1,92 s** | 2,82 s |
| CPU (OpenMP sobre picos) | **8,05 s** | — | — | — | — |

Resultado contrário à expectativa do plano para CPU: lotes **pioram** lá (19 s),
porque o decoder é pequeno, o overhead por chamada já é baixo e o conjunto de trabalho
maior prejudica a localidade de cache. Mantive os **dois caminhos**, como o plano
previa: GPU usa lote 32 serial; CPU mantém o fan-out OpenMP com sessão de 1 thread.
Ajustável por `PLIMAN_SAM_DECODER_BATCH`.

**§2.4 — ZNCC denso.** Duas etapas, ambas exatas quanto às instâncias:

- *Opção C* — kernel reorganizado em **spans contíguos** (SoA) com `#pragma omp simd`,
  pesos `w·t_zero` pré-multiplicados. Transforma um gather sobre um disco em cargas
  sequenciais.
- *Opção A* — backend **FFT** novo (`src/dense_zncc.cpp`, `src/dense_zncc.h`), usando
  o FFTW já linkado. As seis somas ponderadas do ZNCC são correlações com kernel fixo,
  logo produtos de espectros; custo passa a ser independente do raio.

Medido em build otimizado, estágio isolado com encoder quente:

| Imagem | Raio | Direto | FFT | Ganho |
|---|---:|---:|---:|---|
| `coffee.jpg` (objeto grande) | 48 | 5,31 s | **3,37 s** | 1,57× |
| `soybean.jpg` | 15 | 0,37 s | — | FFT não ativa (correto) |
| `flaxgrain.jpg` (grãos) | 7 | ~0 s | — | FFT não ativa (correto) |

**Correção de expectativa:** o plano estimava 10–100× para o modo denso, baseado na
minha medição inicial de 18,27 s. Em `-O2` o passe denso custa **0,4 s ou menos** nos
raios pequenos para os quais o modo existe (grãos, sementes, esporos) e só chega a
5,3 s no teto de raio 48 — caso em que `feat_res = 256` seria a escolha certa. O
ganho real do FFT é 1,57× no pior caso, com as instâncias preservadas. O portão por
raio (`PLIMAN_ZNCC_FFT_MIN_RADIUS`, padrão 20) garante que abaixo disso o custo das
transformadas não seja pago.

**PerSAM ponta-a-ponta** (imagem 2448², 173/169 objetos, `-O2`):

| | `feat_res` | frio | re-execução (ajuste de limiar) |
|---|---|---:|---:|
| `sam_device = "cpu"` | 256 | 12,83 s | 8,57 s |
| `sam_device = "gpu"` | 256 | **7,44 s** | **2,30 s** |
| `sam_device = "cpu"` | 1024 | 12,91 s | 6,53 s |
| `sam_device = "gpu"` | 1024 | **7,71 s** | **1,73 s** |

**Diagnóstico novo:** `pliman_zncc_fft_stats_cpp()` conta tentativas e sucessos do
backend FFT — usado pelos testes para garantir que o caminho rápido é realmente
exercitado, e não silenciosamente ignorado.

**Testes:** `tests/testthat/test-persam-dense-zncc.R` compara FFT contra avaliação
direta (contagem, geometria das caixas, correlação do mapa > 0,9999) e verifica que o
portão troca de caminho de fato. Suíte completa: 23 asserções, todas passando.

### Não implementado (Fases 2–5)

Pendentes após a Fase 2:

- **§2.3 (parte)** — supressão antecipada de picos já cobertos por máscaras aceitas.
  Deixei de fora deliberadamente: ela decide por *centro de pico dentro de máscara*,
  critério mais agressivo que o NMS atual, que protege explicitamente objetos convexos
  encostados (grãos, toras) cujas caixas se sobrepõem mas cujas máscaras são
  disjuntas. Implementar como opt-in, com medição de recall em cena densa.
- **§2.5** — encoder INT8 (CPU) / FP16 (GPU), grafos pré-otimizados em disco, cache de
  embeddings em disco.
- **§3.2-1** — variante FP16 do Grounding DINO (o arquivo atual é INT8 dinâmico e
  provavelmente fragmenta a partição no DirectML; confirmar com o profiler do ORT).
- **§3.2-3** — pós-processamento do Grounded-SAM em C++: hoje retorna uma matriz
  `NumericMatrix` 256×256 em double por instância (≈100 MB para 200 caixas).
- **§4.1** — EPs CUDA/TensorRT/CoreML. **§4.2** — IoBinding.
- **§5.1, §6.3** — batch de tiles e de frames. **§6.1** — pipeline assíncrono de
  vídeo. **§6.2** — overhead R por frame. **§7.2** — `save_image` nativo no lote.
- **§8** — dispatch SIMD em tempo de execução (AVX2/FMA) para os kernels quentes.

---

## STATUS — FASES 3 A 5

Hardware: Intel i7 (12 threads), NVIDIA GTX 1650 (4 GB, **sem Tensor Cores**), Windows 11,
ORT 1.20.1 (build DirectML e build CUDA oficial), CUDA 12.4 / cuDNN 9 vindos do PyTorch,
build `-O2`. Testes: `tests/testthat/test-onnx-providers-pipeline.R` (48 asserções);
suíte completa passando.

### Fase 3 — GPU

**§4.1 Provedores de execução.** Toda sessão passa por `get_or_create_cached_session()`,
que agora anexa o provedor pedido em `PLIMAN_GPU_PROVIDER` (`auto | dml | cuda | tensorrt |
coreml | openvino`), com cadeia de fallback provedor pedido → CUDA → DirectML → CPU.
O padrão `auto` reproduz o comportamento anterior. API R: `pliman_gpu_provider()`
(`cuda_path` para as DLLs CUDA, cache de engines TensorRT, `trt_fp16`),
`onnx_available_providers()`, `onnx_install(engine = "cuda")` (build oficial do GitHub
em `bin/cuda/`). `pliman_session_info()` ganhou `provider`, `graph_cache`, `create_ms`.
Testado de fato: CUDA ativo; TensorRT sem as libs → CUDA; CUDA sem as DLLs → CPU, sem erro.
O `engine = c("gpu","cpu")` dos ~40 runners não mudou: o provedor é uma configuração
de processo, herdada por workers `mirai`.

**§1.2 Profiler** — `pliman_onnx_profiling()` + `onnx_profile_summary()`.

**§3.2-1 Grounding DINO.** O profiler **refutou a hipótese** do plano para o DirectML:

| | DML (INT8) | CUDA (INT8) |
|---|---|---|
| Tempo/chamada (regime) | **0,8 s** (CPU: 8,6 s) | 3,2 s |
| Nós na CPU | 39 triviais (0,2% do tempo) | 392 MatMulInteger + 291 DynamicQuantizeLinear (68%) + 693 cópias |

O FP16 público (onnx-community, irmão exato do arquivo INT8 distribuído) piorou tudo
(CUDA 7,3 s; DML 14,5 s): `GridSample` FP16 não tem kernel CUDA no ORT 1.20 (48 nós na
CPU, 69% do tempo). Com os `GridSample` envolvidos em casts FP32: **CUDA 1,75 s (1,8×
sobre o INT8)**, mesmas caixas. No DML a variante segue 14,5 s, por isso ela é escolhida
automaticamente só com CUDA/TensorRT (`.resolve_dino_model()`); o DML mantém o INT8.
O runner do DINO passou a aceitar E/S FP16 (§4.3).

**§4.2 IoBinding** no decoder SAM em GPU (features ligadas uma vez por imagem):
**sem ganho mensurável** (DML 0,61 s vs 0,61–0,66 s; CUDA 0,50 vs 0,48–0,50 s). A
reutilização dos buffers de saída pelo binding dava ~25% no DML, mas **perdia uma
instância** (resultado errado), então nunca é feita. Mantido opt-in
(`PLIMAN_SAM_IOBINDING=1`) para GPUs com barramento lento.

**Matriz YOLO** (`yolo_benchmark()`, 1600×900, média de 30):

| Modelo | CPU | DirectML | CUDA |
|---|---:|---:|---:|
| yolo26n | 95,7 ms | **13,3 ms** | 25,0 ms |
| yolo26s | 131,0 ms | **19,7 ms** | 22,7 ms |
| yolo26m | 394,3 ms | **41,7 ms** | 50,3 ms |

Nesta GPU o DirectML vence o CUDA em YOLO e empata no encoder SAM (0,36 vs 0,39 s).
O CUDA ganha no DINO (com a variante mista) e no lote de tiles.

### Fase 4 — Vídeo / lote

**§7.2 `save_image` nativo** (`.save_detections_native()`): decodifica direto do arquivo,
desenha com `draw_yolo_detections_bgr_cpp()` e grava com stb (JPEG/PNG/BMP; PNG novo em
`save_frame_bgr_cpp`). 8 imagens: 0,36 s vs 0,61–2,36 s pelo `grDevices`. Fallback
automático para outros formatos, keypoints ou `options(pliman.native_save = FALSE)`.

**§5.1 Lote de tiles SAHI** com exports batch-dinâmicos: um `Run` `[B,3,S,S]`, e cada
tile recebe uma *view* da sua fatia (`g_yolo_inject`), pós-processada pelo mesmo código
do caminho por tile. Caixas iguais a ≤ 2,4e-4 px. 24 tiles:

| | B=1 | B=4 | B=8 | B=16 |
|---|---:|---:|---:|---:|
| CUDA | 0,42 s | 0,36 s | **0,35 s** | 0,35 s |
| DirectML | **0,43 s** | 0,78 s | 0,56 s | 0,49 s |
| CPU | **1,30 s** | 1,50 s | 1,61 s | 1,67 s |

Padrão B=8 só em CUDA/TensorRT (`PLIMAN_YOLO_TILE_BATCH`). **Os modelos YOLO distribuídos
têm batch fixo = 1**, então o caminho só atua com `yolo_export(dynamic = TRUE)`.

**§6.1 Pipeline assíncrono de vídeo**: pré-busca (`yolo_prefetch_cpp`). O quadro *i+1*
é pré-processado e o `Run` dele vai para uma `std::thread`; o `run_yolo_cpp()` seguinte
consome as saídas se o quadro for byte-idêntico (hash completo). Detecções **idênticas**
(716/716 pré-buscas usadas). **Sem ganho de vazão**: yolo26n/DML 30,6 → 28,5 fps;
yolo26s/CPU dentro do ruído. Motivo medido: o trabalho R fora da inferência já é de 2–5
ms por quadro (tracker, desenho e codificação estão em C++), então quase todo o `Run`
continua exposto. Fica **opt-in** (`options(pliman.video_pipeline = TRUE)`).

**Não implementados, com justificativa:** §6.3 lote de quadros (mesmo bloqueio do batch
fixo, e o ganho em DML seria negativo como nos tiles); §6.2 overhead R por quadro (medido
em 0,11 ms/quadro em `rdata`, não há o que tirar); laço nativo completo.

### Fase 5 — Modelos

**Grafos pré-otimizados em disco** (`pliman_onnx_graph_cache()`, opt-in, só sessões CPU):

| Sessão (processo novo) | sem cache | com cache |
|---|---:|---:|
| SAM 2.1 encoder | 3.254 ms | **178 ms** |
| Grounding DINO | 1.700 ms | **690 ms** |
| SAM 2.1 decoder | 206 ms | 102 ms |

Caixas e rótulos idênticos. Chave = tamanho + mtime + versão da API; gravação atômica.
É opt-in porque duplica o espaço em disco dos modelos grandes (~360 MB aqui).

**Encoder SAM INT8/FP16** (`pliman_build_model_variant()`, `options(pliman.sam_encoder)`):

| Encoder | FP32 | INT8 | FP16 |
|---|---:|---:|---:|
| CPU | **2,31 s** | 2,31 s (IoU 0,995) | 5,75 s |
| DirectML | **0,36 s** | 38 s, **IoU 0** | 36,7 s |
| CUDA | **0,39 s** | 4,65 s | 1,52 s |

Nenhuma variante ganha nesta máquina (CPU sem VNNI, GPU sem Tensor Cores). O padrão
continua FP32, sem seleção automática; INT8 é recusado em sessões GPU (máscaras erradas
no DML). A conversão FP16 exigiu realinhar os `Cast` originais do modelo (o
`onnxconverter_common` deixava o grafo inválido).

**Não implementado:** encoders leves (EfficientSAM/EdgeSAM), que exigem novos modelos e
outra interface de decoder. Cache de embeddings em disco (§2.5-4).

---

## 0. Resumo executivo

### Os 12 itens de maior impacto (ordenados por ganho/esforço)

| # | Item | Onde | Tipo | Ganho esperado | Esforço |
|---|------|------|------|----------------|---------|
| 1 | `clear_onnx_sessions_cpp()` derruba **todas** as sessões a cada chamada do Grounded-SAM em GPUs com < 6 GB | `R/utils_neural.R:4265-4273` | Bug de desempenho | Elimina 5–30 s de recriação de sessão por imagem em lote | Baixo |
| 2 | Decoder do PerSAM chamado **1 prompt por `Run`**, sem supressão antecipada de picos já cobertos | `onnx_inference.cpp:2644-2802` | Gargalo | 2–5× na etapa de decoder | Médio |
| 3 | Modo denso (`feat_res = 1024`) faz ZNCC direto O(N·r²·P) sobre 1M pixels | `onnx_inference.cpp:2419-2533` | Gargalo | 10–100× (FFT ou 512 px) | Médio |
| 4 | Modo denso com protótipos pré-computados **não tem kernels** → mapa zerado → 0 detecções (e o loop de 1M pixels roda à toa) | `onnx_inference.cpp:1935-2290, 2453` | **Bug funcional** | Corrige lote de autolabel | Baixo |
| 5 | Cache do encoder SAM usa hash de **4096 amostras** → vídeo/time-lapse com câmera fixa reutiliza embeddings velhos | `onnx_inference.cpp:1785-1809` | **Bug funcional** | Corrige resultados errados | Baixo |
| 6 | Encoder e decoder do PerSAM forçados em CPU; decoder (pequeno, sem Hiera) poderia ir para GPU em lote | `onnx_inference.cpp:1764-1774`, `utils_neural.R:4677-4685` | Gargalo | 3–10× no decoder (GPU) | Médio |
| 7 | Grounding DINO é **INT8 dinâmico** (MatMulInteger/ConvInteger/DynamicQuantizeLinear) e está sendo executado no DirectML | modelo `groundingdino-tiny.onnx` | Gargalo provável | 2–5× em GPU com variante FP16 | Médio |
| 8 | Tensores trafegam R(double 8 B) → `as.numeric()` (cópia) → `std::vector<float>`: 3 cópias de 12–25 MB por chamada + pré-processamento single-thread em double | `image_resize.cpp:335`, `utils_neural.R:4763/4288/4129` | Overhead | 50–300 ms/imagem | Baixo |
| 9 | `image_data(curr_img, type = "numeric")` no loop de autolabel: imagem inteira em double (≈ 290 MB para 12 MP) por classe/iteração | `utils_yolo_train.R:882, 1015` | Overhead/memória | 0,3–1 s/imagem + GC | Baixo |
| 10 | Tiles SAHI e frames de vídeo rodam **1 por `Run`**, sem batch, mesmo quando o modelo tem batch dinâmico | `onnx_inference.cpp:4892-4903`, `utils_video.R:1856` | Gargalo (GPU) | 1,5–3× em GPU | Médio |
| 11 | Laço de vídeo síncrono: inferência, tracker, `data.frame()`s e desenho na mesma thread R | `utils_video.R:2984-3500` | Gargalo | 1,3–2× FPS em modelos nano | Alto |
| 12 | `save_image = TRUE` no lote renderiza via `grDevices::jpeg/png + plot()` | `utils_dl_batch.R:495-560` | Overhead | 0,5–2 s/imagem | Baixo |

### Bugs funcionais encontrados (corrigir independentemente de desempenho)

1. **PerSAM – hash amostrado do encoder** (item 5).
2. **PerSAM – modo denso com protótipos** sem kernels ZNCC → 0 instâncias (item 4).
3. **PerSAM – truncamento de protótipos descarta os novos exemplares**: os pré-computados entram primeiro (`:1912`) e os novos são anexados (`:2251`); `resize(8)` em `:2266` mantém só os antigos. No review interativo com 8 protótipos acumulados, novos cliques **não alteram** o mapa de similaridade.
4. **`omp_set_num_threads()` global** em `:2639`: após uma chamada PerSAM com `threads = 2`, **todas** as regiões OpenMP do pacote passam a usar 2 threads até o fim da sessão R.
5. **`run_sam2_inference_cpp` quebra prompts de caixa**: os cantos (labels 2 e 3) são enviados como **dois prompts de 1 ponto independentes** (`:1160-1270`) em vez de um prompt `[1, 2, 2]`. Caixa vira "dois cliques com labels inválidos".
6. **Grounding DINO – `text_threshold` é ignorado** (parâmetro em `:1295` nunca usado) e as dimensões `900` queries × `256` tokens estão fixas no código (`:1394-1403`); um export diferente lê memória fora do tensor.
7. **Distorção de aspecto**: SAM (1024) e DINO (800) recebem a imagem esticada (`letterbox = FALSE`, `utils_neural.R:4112, 4282, 4748`). Não é velocidade, mas degrada máscaras/caixas em imagens não quadradas — vale corrigir junto com a reescrita do pré-processamento.
8. **OpenMP aninhado** em `run_yolo_cpp` (`:4510` paralelo por instância contém `:4654` paralelo por linha). Com `OMP_MAX_ACTIVE_LEVELS > 1` gera sobre-assinatura.
9. **Inconsistência de dispositivo do SAM**: o PerSAM força CPU alegando travamento do DirectML (`DXGI_ERROR_DEVICE_HUNG`) no Hiera, mas o Grounded-SAM (`run_sam2_instances_cpp`) roda **o mesmo encoder** em GPU. Ou o risco existe no Grounded-SAM, ou a restrição do PerSAM é desnecessária. É preciso decidir com teste (Seção 2.3).

---

## 1. Protocolo de medição (fazer antes de qualquer mudança)

Sem linha de base, otimizações de inferência enganam (aquecimento, cache de sessão, turbo de CPU).

1. **Cronometragem por estágio** dentro do C++ com `std::chrono::steady_clock` e retorno opcional em `attr(res, "timings")` (habilitado por `options(pliman.profile = TRUE)` → flag passada ao C++). Estágios mínimos:
   - PerSAM: `preprocess`, `hash`, `encoder`, `prototypes`, `similarity`, `dense_zncc`, `peaks`, `decoder`, `nms`, `labels`, `R_post`.
   - Grounded-SAM: `dino_pre`, `dino_run`, `sam_pre`, `sam_enc`, `sam_dec`, `R_post`.
   - Vídeo: `decode_wait`, `infer`, `track`, `R_df`, `draw`, `encode_wait`, `display`.
2. **Profiler do ONNX Runtime**: expor `PLIMAN_ORT_PROFILE=<prefixo>` → `EnableProfiling(opts, prefix)` em `make_session_options()`. O JSON mostra **quais nós caíram para o CPU EP** quando a sessão é DirectML (essencial para o item 7).
3. **Conjunto fixo de testes**: 3 imagens (1 MP, 12 MP, 24 MP), 1 vídeo 1080p de 30 s, 1 dataset de 50 imagens; rodar 5× descartando a 1ª (aquecimento). Registrar CPU, GPU, versão do ORT, `threads`.
4. Usar `yolo_benchmark()` e `pliman_session_info()` já existentes; acrescentar `persam_benchmark()` (interno) com os estágios acima.
5. **Critério de aceitação de cada PR**: tempo mediano do estágio-alvo, saída idêntica (ou IoU ≥ 0,99 das máscaras / mesmas caixas) contra a linha de base.

---

## 2. PerSAM (SAM 2.1 + similaridade de protótipos) — foco principal

### 2.1 Como o tempo se distribui hoje (estimativa para 12 MP, CPU 8 núcleos)

| Estágio | Código | Custo estimado | Observação |
|---|---|---|---|
| `image_data(type="numeric")` no autolabel | `utils_yolo_train.R:882` | 0,2–0,6 s | ≈ 290 MB em double |
| `.preprocess_nchw()` → `preprocess_nchw_cpp` | `image_resize.cpp:335` | 0,15–0,4 s | single-thread, double, filtro com muitos taps |
| `as.numeric(sam_tensor)` + `sexp_to_float` | `utils_neural.R:4763`, `:1783` | 20–60 ms | 2 cópias extras |
| Encoder Hiera-Tiny 1024² (CPU fp32) | `:1833` | **0,8–2,5 s** | maior custo fixo |
| Protótipos + similaridade 64×64 | `:1942-2326` | < 20 ms | ok |
| Modo 256 (padrão) | `:2366-2417` | < 10 ms | ok |
| **Modo 1024 (ZNCC denso)** | `:2419-2533` | **2 s – vários minutos** | depende do raio (6–48) e de nº de exemplares (≤ 8) |
| Detecção de picos + NMS espacial | `:2558-2618` | ms – s | O(N_cand × N_aceitos) |
| **Decoder por pico** | `:2644-2802` | **0,3–3 s** | 10–30 ms/pico em 1 thread, ÷ nº threads |
| NMS de máscaras | `:2813-2890` | ms – s | O(n²) com bbox reject |
| Projeção para rótulos full-res | `:2897-2929` | 50–150 ms | `std::round` em double por pixel |
| `contour()` + `poly_measures()` | `utils_neural.R:4984-5007` | 50–300 ms | full-res |

### 2.2 P0 — correções e ganhos rápidos (1–3 dias)

**(a) Hash completo e cache LRU de embeddings** (bug 1)
- Trocar o FNV amostrado por hash de **todos** os bytes. Ideal: hashear a imagem **de origem** (raw/uint8, antes do resize) + dims + caminho do encoder — mais barato que hashear 12 MB de float e invariante ao pré-processamento.
- Implementar um FNV-1a/xxHash64 de 64 bits em blocos de 8 bytes (≈ 2–4 ms para 36 MB); sem dependência nova.
- Converter `g_cached_encoder` em **LRU de K entradas** (`PLIMAN_SAM_CACHE=4` por padrão): cobre multi-classe por imagem, review interativo e "voltar à imagem anterior".
- Guardar as features em **FP16** no cache (16 MB → 8 MB por entrada) quando K > 1; converter para float só ao criar o tensor.

```cpp
static uint64_t hash_bytes64(const void* p, size_t n, uint64_t seed) {
  const uint8_t* b = (const uint8_t*)p; uint64_t h = 1469598103934665603ULL ^ seed;
  size_t i = 0;
  for (; i + 8 <= n; i += 8) { uint64_t v; std::memcpy(&v, b + i, 8);
    h ^= v * 0x9E3779B97F4A7C15ULL; h = (h << 27 | h >> 37) * 0x100000001B3ULL; }
  for (; i < n; ++i) { h ^= b[i]; h *= 0x100000001B3ULL; }
  return h;
}
```

**(b) Remover `omp_set_num_threads(actual_omp_threads)`** (bug 4) — a cláusula `num_threads(...)` do `pragma` já basta.

**(c) Protótipos: novos exemplares primeiro + kernels serializados** (bugs 2 e 3)
- Ordem: protótipos dos cliques atuais **antes** dos pré-computados; compactação por Farthest-Point Sampling no C++ (hoje só existe em R, `.compact_prototypes`, `utils_neural.R:4576`).
- Serializar os kernels ZNCC (`dx, dy, w, t_zero[3]`, `tmpl_std`, `sum_w`) em `prototypes$dense` para que o modo 1024 funcione nas imagens 2..N do autolabel.
- Se `num_kernels == 0`, **cair para o modo 256** com `cli_alert_warning` em vez de rodar 1M pixels que retornam zero.

**(d) Entrada raw direto no C++, pré-processamento em float e paralelo**
- Nova assinatura (mantendo a antiga como alias para retrocompatibilidade): `run_sam2_persam_cpp(img = <array raw/double [W,H,C]>, ...)`. O C++ faz resize 1024 + normalização em **float**, com `#pragma omp parallel for` por linha e filtro bilinear de 2 taps (ou área quando o fator de redução > 2), escrevendo direto no buffer do tensor.
- No R: deixar de chamar `.preprocess_nchw()` e `as.numeric()`; no autolabel, trocar `image_data(curr_img, type = "numeric")` pelo array nativo (raw quando disponível).
- Mesmo helper usado por SAM-ponto, Grounded-SAM e DINO (ver 3.2) — um único `preprocess_to_nchw_f32()`.

**(e) Projeção de rótulos com LUT**
- Pré-computar `col_lut[x] = round((x+0.5)/W*256-0.5)` e `row_lut[y]` uma vez; laço por instância vira só lookups; paralelizar por linhas com "primeiro escreve vence" via ordem de score (como já feito no YOLO com `__atomic_compare_exchange_n`).

**(f) NMS espacial por grade** (`:2604-2618`)
- Grade de células de lado `grid_min_dist`; cada candidato só compara com aceitos nas 9 células vizinhas → O(N).

### 2.3 P1 — decoder em lote com supressão antecipada (maior ganho do PerSAM)

O decoder exportado aceita lote: `point_coords [num_labels, num_points, 2]`, `mask_input [num_labels, 1, 256, 256]`, `has_mask_input [num_labels]` (confirmado inspecionando `sam2.1.decoder.onnx`). Hoje cada pico é um `Run` separado e todos os picos são decodificados antes do NMS — inclusive vários picos **do mesmo objeto**.

**Algoritmo proposto** (substitui `:2620-2890`):

```
ordenar picos por similaridade (desc)
coverage = bitmap 256x256 zerado
while restam picos:
    lote = próximos B picos cujo (r,c) em coverage == 0     # B = 16..64
    se lote vazio: break
    Run(decoder, lote)                                      # 1 chamada
    para cada máscara do lote (em ordem de score):
        selecionar melhor das 3 (regras atuais)
        NMS contra instâncias aceitas (grade espacial)
        se aceita: OR da máscara em coverage
    descartar picos pendentes que caíram em coverage
```

- Picos dentro de objetos já segmentados **não chegam ao decoder**. Em objetos grandes ou `min_dist` pequeno isso costuma eliminar a maior parte das chamadas.
- Exemplares: reaproveitar a máscara já calculada na extração de protótipos (`:1973`) em vez de decodificar de novo como pico.
- Buffers `mask_input`/`has_mask` alocados **uma vez** no tamanho máximo do lote (hoje 16 MB por lote em `run_sam2_instances_cpp:1576`).
- **CPU**: comparar (i) lote B=32 em sessão com `intra_op = threads` vs (ii) o esquema atual OpenMP × sessão 1-thread. Manter ambos e escolher por heurística (`n_picos < 2·threads` → (ii)). Lote maior gera GEMMs maiores, mais eficientes, mas o decoder replica o embedding por prompt — limitar B pela memória (≈ 4 MB × B).
- **GPU**: aqui o lote é decisivo. O decoder não contém o backbone Hiera (que é o suspeito do travamento DirectML). Proposta: **encoder em CPU, decoder em DirectML**, com `B = 64`. O upload de 16 MB de features por `Run` fica amortizado no lote; na versão P2, usar IoBinding para subir as features uma vez por imagem.
- Expor `sam_device = c("auto", "cpu", "gpu")` e `encoder_device` separados (mantendo `engine` com o comportamento atual como padrão) e **documentar** que o ganho depende de GPU compatível.

**Teste de estabilidade do encoder em GPU** (resolve o bug 9): script de estresse com 200 encodes consecutivos no DirectML, em (a) ORT 1.20.1 atual e (b) versão mais recente do ORT/DirectML. Se estável → liberar `encoder_device = "gpu"` opt-in; se instável → forçar CPU também no Grounded-SAM.

### 2.4 P1 — ZNCC denso (`feat_res = 1024`) sem custo O(r²)

Para cada exemplar p, os termos são correlações com kernels fixos:

- `dot = Σ_c (I_c ⋆ (w·t_c))`, `s_c = I_c ⋆ w`, `Σsq = Σ_c (I_c² ⋆ w)`
- `var = Σsq − Σ_c s_c² / Σw`, `zncc = dot / (√var · tmpl_std)`

**Opção A — FFT (recomendada, FFTW3 já está em `PKG_LIBS`)**:
- FFT direta de `I_c` e `I_c²` (6 transformadas) **uma vez por imagem**; por exemplar, FFT de `w` e `w·t_c` (4) e 5 inversas (`dot` e `Σsq` somados no domínio da frequência, `s_c` × 3).
- Tamanho 1024 + 2·48 → pad para 1152 (fator suave para FFTW). Custo praticamente **independente do raio**; ≈ 9 FFTs por exemplar.
- Planos FFTW em cache estático (`FFTW_MEASURE` uma vez por sessão; usar `fftw_plan_with_nthreads` se `fftw3_threads` estiver disponível; caso contrário paralelizar exemplares com OpenMP — atenção: o *planner* do FFTW não é thread-safe, criar planos fora da região paralela).

**Opção B — resolução reduzida + refinamento**: ZNCC a 512² com kernel reduzido (raio/2) → ≈ 16× menos operações; refinar só em janelas 9×9 ao redor dos picos a 1024. Mais simples, pequena perda em grãos muito pequenos.

**Opção C — mínimo imediato**: reorganizar o kernel em **spans por linha** (SoA: `dx_ini, dx_fim, dy` + arrays de pesos/templates contíguos) para o laço interno ser contíguo e vetorizável (`#pragma omp simd`), e restringir o cálculo a ROIs com `vit_1024 ≥ vit_thresh` dilatadas. Ganho 3–8×.

Sugestão: C já (baixo risco), A em seguida, mantendo o caminho direto como referência de teste.

### 2.5 P2 — custo fixo do encoder

1. **Quantização do encoder para CPU** (INT8 dinâmico/estático via `onnxruntime.quantization`, ou QDQ estático com ~100 imagens agrícolas de calibração). Típico 1,5–2,5× em CPUs com AVX-VNNI; **validar mIoU** contra fp32. Distribuir como `sam2.1.encoder.int8.onnx` e escolher automaticamente quando `engine = "cpu"`.
2. **FP16 para GPU** (`keep_io_types=True`, assim as entradas/saídas continuam float32 e o C++ não muda). Documentar: ganho exige suporte nativo a FP16 (Tensor Cores/equivalente); em GPUs sem suporte pode ficar mais lento por casting.
3. **Modelo ORT pré-otimizado**: salvar o grafo otimizado na 1ª carga (`SetOptimizedModelFilePath` → `*.opt.onnx` em `pliman_model_dir()`), recarregar com `ORT_ENABLE_BASIC`. Corta o tempo de criação de sessão (relevante após `pliman_clear_sessions()` e em workers `mirai`).
4. **Cache em disco de embeddings** (opcional, `cache_dir =`): FP16, chave = hash do arquivo + mtime + encoder. Útil em datasets revisitados (review, re-export).
5. **Encoders mais leves como opção** (EfficientSAM / EdgeSAM / entrada 768): só para o mapa de similaridade; manter Hiera-Tiny para as máscaras finais quando a precisão importar.

### 2.6 P2 — lado R do PerSAM
- `sim_norm`/`as_image` ok; `contour()` + `poly_measures()` sobre a imagem inteira: gerar contornos **por instância no C++** a partir da máscara 256 reamostrada na janela do bbox (mesma ideia de `.instance_window_mask`), evitando varrer W×H.
- Quando `mask = TRUE` é chamado por `.compute_dl_mask` só para obter a máscara binária (`p_res$mask * 1.0`, `utils_neural.R:2896`), pular contornos/medidas.
- `p_res$mask * 1.0` cria um double W×H: retornar lógico/inteiro.

---

## 3. Grounded-SAM (Grounding DINO + SAM 2.1)

### 3.1 P0
1. **Remover o `clear_onnx_sessions_cpp()` incondicional** (`utils_neural.R:4265-4273`). Substituir por: (a) evicção **somente da sessão do DINO** após a detecção quando VRAM < limite (nova `onnx_evict_session_cpp(model_path)` — a função interna `evict_cached_session_by_key` já existe); (b) ou tentar e, em caso de erro de alocação, evictar e repetir. Com o código atual, todo lote/autolabel com Grounded-SAM recria DINO + SAM a cada imagem.
2. **`text_threshold`**: aplicar no nível de token (um box só é mantido se o token vencedor ≥ `text_threshold`), como no Grounding DINO original — ou documentar que é ignorado. Ler `num_queries`/`num_tokens` do shape de `logits` em vez de `900`/`256` fixos.
3. Sigmoide em float e só para `t ∈ [1, seq_len-2]` (já é) — trocar `double`/`std::exp` por `float`/`std::exp` e pular a query cedo quando `max_logit < logit(box_threshold)` (compara logits, sem sigmoide).
4. `pixel_mask` int64 800×800 (5 MB) alocado e zerado a cada chamada → buffer estático reutilizado.

### 3.2 P1
1. **Variante do DINO por EP**: o arquivo atual é INT8 dinâmico (392 `MatMulInteger`, 291 `DynamicQuantizeLinear`, 5 `ConvInteger`, 48 `GridSample`). Quantização dinâmica é pensada para CPU; no DirectML esses operadores podem não ter kernel e cair em partições CPU com cópias GPU↔CPU a cada camada. **Confirmar com o profiler** (Seção 1.2) e, se confirmado, publicar `groundingdino-tiny.fp16.onnx` para `engine = "gpu"` e manter o INT8 para CPU, escolhendo automaticamente.
2. **Encoder SAM compartilhado**: `run_sam2_instances_cpp` e `run_sam2_inference_cpp` não usam o cache de embeddings do PerSAM. Extrair `get_sam2_features(img, enc_path, device)` usado pelos três runners — Grounded-SAM → PerSAM na mesma imagem passa a custar só o decoder.
3. **Pós-processamento em C++**: hoje cada instância volta como `NumericMatrix` 256×256 double (512 KB; 100 MB para 200 caixas, `:1639-1648`) e o R faz `.instance_window_mask()` por instância. Retornar direto a matriz de rótulos (como o PerSAM faz) com reamostragem bilinear restrita ao bbox, e opcionalmente logits float só da ROI.
4. **Pré-processamento único**: DINO (800) e SAM (1024) partem da mesma imagem; gerar ambos num só passe C++ a partir do raw (ver 2.2-d) e, com isso, corrigir o letterbox (bug 7).

### 3.3 P2
- Para vários prompts de texto na mesma imagem, cachear o resultado do DINO por `(hash_imagem, prompt)` — útil no review interativo.
- Caixas muito grandes (fundo) podem ser descartadas antes do SAM (`filter_nested` hoje é um laço R O(n²) em `:4188-4216`; mover para `filter_nested_boxes_cpp`, que já existe).

---

## 4. Motor ONNX comum (`src/onnx_inference.cpp`)

### 4.1 Provedores de execução (maior alavanca de hardware)
- Hoje GPU = **somente DirectML no Windows** (`utils_neural.R:2131-2135, 2297-2310`); Linux/macOS sempre CPU.
- Adicionar carregamento opcional, pela mesma API C dinâmica, de:
  - **CUDA EP** (`SessionOptionsAppendExecutionProvider_CUDA_V2`) e **TensorRT EP** (`..._TensorRT_V2`, com cache de engine em `pliman_model_dir()`) para NVIDIA em Windows/Linux;
  - **CoreML EP** no macOS; **OpenVINO EP** para iGPU/NPU Intel.
- `engine = c("gpu", "cpu", "cuda", "tensorrt", "dml", "coreml", "openvino")` mantendo `"gpu"` como alias de "melhor disponível". Documentar explicitamente a dependência de hardware/driver e que FP16/TensorRT só aceleram com suporte nativo (AGENTS.md §4).
- Avaliar atualizar o ORT de **1.20.1** para uma versão mais recente (medir antes/depois; DirectML recebe correções de estabilidade a cada release).

### 4.2 Sessões
- **IoBinding** para modelos de shape estática usados em laço (vídeo, tiles): entrada/saída pré-alocadas, sem alocação por `Run`; no DirectML permite manter saídas na GPU quando o próximo modelo as consome (encoder → decoder SAM).
- **Aquecimento** hoje só para sessões GPU com shape totalmente estática (`:810`). Adicionar `warmup_shape` opcional (o vídeo já faz um warm-up manual; `image_detect_dl()` em lote não).
- **Grafo otimizado em disco** (ver 2.5-3) para todos os modelos grandes (DINO, SAM, RMBG/BiRefNet).
- `run_onnx_inference_cpp` pede `SessionGetInputName/OutputName` e aloca/libera nomes a cada chamada (`:919-923`); usar `SessionMeta` como o YOLO já faz.
- Saídas FP16→double em laços escalares (`half_to_double`) → converter em float com OpenMP e só no fim para `double` se o R precisar.
- **Concorrência ORT × OpenMP**: nos pipelines que alternam `Run` e laços OpenMP pesados (PerSAM: encoder → ZNCC → decoder), os workers do ORT ficam em *spin* e competem com o OpenMP. Para essas sessões, usar `session.intra_op.allow_spinning = 0` por sessão (não globalmente, já que no vídeo o spin ajuda).

### 4.3 FP16 / INT8 nos runners especializados
- Os runners SAM/DINO criam tensores **sempre `FLOAT`**; um modelo FP16 com I/O em FP16 falharia. Padronizar via `session_input_type()` (já existe) ou exigir exports com `keep_io_types=True`.

---

## 5. YOLO — detecção e segmentação em imagem

`run_yolo_cpp` já está bem otimizado (pré-processamento fundido, argmax por blocos, NMS por classe, máscara só no bbox). Pontos restantes:

1. **Batch real** para tiles SAHI e lotes de imagens quando `meta->batch_dynamic` (`:4892-4903`): montar `[N,3,S,S]` (N = 4–16), um `Run`, pós-processar cada fatia. Em GPU o ganho é grande para modelos n/s, cujo custo por chamada é dominado por overhead.
2. **Crop + pré-processamento dos tiles em paralelo** (threads OpenMP por tile para preencher o lote) e, quando `tile == target_size`, pular o resize (hoje sempre passa pelo bilinear).
3. **`static std::vector<float> input_vals`** (`:3804`): mover para buffer por sessão em `SessionMeta` (evita retenção de memória de shapes antigos e prepara uso multi-thread).
4. **Candidatos**: `YoloCand` carrega dois `std::vector` (máscara/keypoints) por candidato; com `conf` baixo são milhares de alocações. Separar caixas (POD) dos dados extras (já coletados só após o NMS).
5. **OpenMP aninhado** (bug 8): quando `par_inst`, desabilitar o `parallel for` interno (`if(!par_inst && ...)`).
6. **Rótulos full-res**: `labels` W×H int (33 MB em 4K) é alocado sempre que `need_masks`. Oferecer `mask_format = c("labels", "roi")` retornando máscaras recortadas por bbox; `extract_contours_cpp(lbl_mat)` passa a rodar por ROI.
7. **`.run_yolo()` no R** com morfologia (`utils_neural.R:5263-5308`): `df_boxes[k, ]` dentro do laço (acesso a linha de data.frame é lento) → usar a matriz `raw_res$boxes`.

---

## 6. Vídeo (`video_detect_dl`)

Já existe: leitor Media Foundation com *prefetch* de 2 slots e gravação em thread, tracker em C++, `infer_every`. O gargalo restante é a **thread R** fazendo tudo em série.

### 6.1 P1 — pipeline assíncrono de inferência
- Novo par `video_infer_submit_cpp(engine_ptr, frame_raw)` / `video_infer_collect_cpp(engine_ptr)`: uma `std::thread` em C++ executa pré-processamento + `Run` + pós-processamento do frame N+1 enquanto o R faz tracker/desenho/HUD do frame N. Regras: a thread **não toca objetos R** (copia o frame para `std::vector<uint8_t>`; resultados em structs C++, convertidos para `Rcpp::List` só no `collect`).
- Latência +1 frame, throughput ≈ max(inferência, resto) em vez de soma.

### 6.2 P1 — reduzir overhead R por frame
- Hoje há 2–3 `data.frame()` por frame (`utils_video.R:3063, 3115, 3170`), `ifelse` de labels, `.get_class_palette()` e `data.frame` de keypoints por objeto (`:3310-3320`). Guardar por frame **matrizes/listas simples** e montar o `data.frame` final uma vez com `.rbind_df_fast()`; cachear a paleta por classe/ID.
- Máscaras em vídeo: desenhar direto no buffer BGR a partir do protótipo (C++), sem materializar a matriz de rótulos full-res por frame.
- SAHI + segmentação ainda usa laço R por tile (`:1966-2006`); estender `run_yolo_tiles_cpp` para devolver `seg_stats` e usar o caminho único.

### 6.3 P2
- **Batch de frames** para arquivo (não câmera): decodificar 4–8 frames, um `Run` com batch, desenhar em ordem — exige export com batch dinâmico.
- **Laço nativo completo** (decodificar → inferir → rastrear → desenhar → codificar em C++) quando `show` não usa o dispositivo gráfico R; o R recebe agregados a cada K frames.
- `infer_every` adaptativo (pular inferência quando a diferença entre frames é baixa; o tracker já propaga caixas).
- Câmera usa `MF_READWRITE_ENABLE_HARDWARE_TRANSFORMS = FALSE` (`native_camera.cpp:111`) enquanto o leitor de arquivo usa `TRUE`; testar `TRUE` também na câmera.
- Suportar `as_job = TRUE` com `rstudioapi::jobRunScript()` para processamento offline longo (AGENTS.md §4), se ainda não coberto em todos os modos.

---

## 7. Lote (`image_*_dl_batch`, SAHI)

1. **GPU**: workers `mirai` em GPU duplicam sessões e VRAM (limitado a 2). Para GPU, preferir **um único processo** com pipeline: threads C++ decodificando JPEG (`stb_image` já incluso) e pré-processando o lote N+1 enquanto a GPU roda o lote N.
2. **`save_image`** (item 12): trocar `grDevices::jpeg/png` + `plot()` + `.plot_yolo_bboxes()` por `draw_yolo_detections_bgr_cpp()` + `save_frame_bgr_cpp()` (já usados no vídeo). Manter o caminho R como fallback para estilos não suportados.
3. Reimportação da imagem só para salvar (`utils_dl_batch.R:514-516`): reutilizar o buffer raw da inferência.
4. CPU com `mirai`: `.tune_parallel_dl` divide núcleos físicos por worker (correto); somar `allow_spinning = 0` nos workers para evitar *spin* cruzado entre processos.

---

## 8. Compilação e SIMD

- `src/Makevars*` só têm `$(SHLIB_OPENMP_CXXFLAGS)`; no Rtools isso significa `-O2` + SSE2. Os laços quentes (ZNCC, pré-processamento, decodificação de máscaras, conversões FP16) não usam AVX2/FMA.
- Abordagem compatível com CRAN: *dispatch* em tempo de execução para os 4–5 kernels quentes com `__attribute__((target("avx2,fma")))` + `__builtin_cpu_supports("avx2")`, mantendo versão genérica. Não alterar flags globais.
- Conversão float↔half: usar instruções F16C (`_mm256_cvtps_ph`) no caminho AVX2.

---

## 9. Roteiro por fases

| Fase | Duração | Conteúdo | Saída |
|---|---|---|---|
| **0 – Base** | 1–2 dias | Timers por estágio, profiler ORT, script de benchmark fixo | Tabela de linha de base |
| **1 – Bugs + quick wins** | 3–5 dias | Itens 1, 4, 5; bugs 3–6 e 8; 2.2 (a)–(f); 3.1; 5.3–5.5; 7.2 | PR "persam-fixes", "grounded-sam-sessions" |
| **2 – PerSAM rápido** | 1–2 semanas | Decoder em lote + supressão antecipada (2.3); ZNCC opção C e depois A (2.4); encoder compartilhado (3.2-2); pós-processamento C++ (3.2-3) | PerSAM ≥ 3× no modo 256, ≥ 10× no modo 1024 |
| **3 – GPU** | 1–2 semanas | Teste de estabilidade DML do encoder; decoder em GPU; DINO FP16; IoBinding; CUDA/TensorRT EP opcionais | Matriz CPU/DML/CUDA documentada |
| **4 – Vídeo/lote** | 2–3 semanas | Pipeline assíncrono (6.1); overhead R (6.2); batch de tiles/frames (5.1, 6.3); `save_image` nativo | FPS e imagens/s medidos |
| **5 – Modelos** | contínuo | Encoder INT8 (CPU) e FP16 (GPU), grafos pré-otimizados, encoders leves opcionais | Novos arquivos no repositório de modelos |

### Regras do repositório a respeitar durante a implementação
- `R/utils_neural.R` (415 KB) e `R/utils_yolo_train.R` (490 KB) já estão perto do limite de 512 KB citado no AGENTS.md: código novo em arquivos temáticos (`R/utils_sam.R`, `R/utils_dl_profile.R`) e C++ novo em `src/sam2_engine.cpp` / `src/dense_zncc.cpp`. Após editar arquivos grandes, checar tamanho e `parse()`.
- Nenhum argumento removido: novos parâmetros (`sam_device`, `encoder_device`, `mask_format`, `cache_dir`) entram com padrão que reproduz o comportamento atual; assinaturas antigas dos `*_cpp` mantidas como wrappers.
- Mensagens via `cli::`; documentação roxygen deixa explícito quando o ganho depende de hardware (FP16, TensorRT, DirectML).
- Convenção `dim(img)[1] = largura`: os novos caminhos C++ que leem raw `[W,H,C]` e escrevem NCHW (`H×W`) precisam de teste de transposição com imagem não quadrada (ex.: 640×360).

### Testes de regressão sugeridos (`tests/testthat/`)
- PerSAM: mesma contagem/IoU de máscaras antes/depois em 3 imagens; modo 1024 com `precomputed_prototypes` retorna > 0 instâncias; duas imagens que diferem em uma pequena região **não** compartilham cache.
- SAM caixa: prompt `c(xmin, ymin, xmax, ymax)` gera máscara contida na caixa.
- Grounding DINO: `text_threshold = 0.99` reduz detecções.
- Threads: `omp_get_max_threads()` igual antes/depois de `image_segment_dl(model = "persam", threads = 2)`.
- Imagem não quadrada em todos os runners (convenção W×H).

---

## Apêndice A — Referências de código

| Tema | Local |
|---|---|
| Opções de sessão | `src/onnx_inference.cpp:499-522` |
| Cache/LRU de sessões | `src/onnx_inference.cpp:673-822` |
| Warm-up (só shape estática em GPU) | `src/onnx_inference.cpp:575-608, 810` |
| SAM por ponto (bug da caixa) | `src/onnx_inference.cpp:1044-1279` |
| Grounding DINO | `src/onnx_inference.cpp:1288-1481` |
| SAM por caixas (Grounded-SAM) | `src/onnx_inference.cpp:1484-1729` |
| PerSAM | `src/onnx_inference.cpp:1747-2963` |
| YOLO | `src/onnx_inference.cpp:3572-4770` |
| Tiles YOLO | `src/onnx_inference.cpp:4821-4923` |
| Pré-processamento NCHW | `src/image_resize.cpp:335-418` |
| `.run_grounded_sam` | `R/utils_neural.R:4059-4342` |
| `.run_persam` | `R/utils_neural.R:4629-5009` |
| Laço de vídeo | `R/utils_video.R:2984-3690` |
| PerSAM no autolabel | `R/utils_yolo_train.R:771-1080` |
| `save_image` do lote | `R/utils_dl_batch.R:495-560` |

## Apêndice B — Inspeção dos modelos locais

| Modelo | Entradas | Observações |
|---|---|---|
| `groundingdino-tiny.onnx` (204 MB) | `pixel_values [1,3,800,800]` fixo, `input_ids [1,seq]`, `pixel_mask [1,800,800]` int64 | **INT8 dinâmico**: 392 MatMulInteger, 291 DynamicQuantizeLinear, 5 ConvInteger; 48 GridSample; saídas `[1,900,256]` e `[1,900,4]` |
| `sam2.1.encoder.onnx` (109 MB) | `image [1,3,1024,1024]` | FP32, opset 18, 2 nós `If` (podem fragmentar a partição no DirectML) |
| `sam2.1.decoder.onnx` (16,5 MB) | `point_coords [num_labels, num_points, 2]`, `mask_input [num_labels,1,256,256]` | **lote dinâmico** → permite o decoder em lote do PerSAM |
