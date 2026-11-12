---
name: pliman-dev-workflow
description: >-
  Runbook de verificação, geração de documentação roxygen2 e carregamento do pacote pliman em Windows com R 4.5.1.
  Use após editar código R ou C++ no repositório pliman para validar limites de tamanho, sintaxe, arquivos .Rd e evitar falhas com pliman.dll travada pelo RStudio.
---

# Fluxo de Verificação e Build do `pliman`

Execute estes passos após editar arquivos em `R/` ou `src/`:

## 1. Verificar Limite de Tamanho (512 KB) e Sintaxe R
Garanta que nenhum arquivo em `R/` ultrapassou 512.000 bytes e que todos fazem `parse()` sem erro:

```powershell
powershell -Command "Get-ChildItem R/*.R | Where-Object { $_.Length -ge 500000 } | Select-Object Name, Length"
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e "for (f in list.files('R', pattern = '\\.R$', full.names = TRUE)) parse(file = f)"
```

## 2. Atualizar Documentação (`roxygen2`)
Sempre que assinaturas, `@param` ou `@export` forem alterados:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e "roxygen2::roxygenise()"
```

## 3. Carregar e Testar Sem Conflito de DLL
No Windows, uma sessão ativa do RStudio frequentemente mantém `src/pliman.dll` travada em memória, fazendo `R CMD INSTALL` falhar com erro de permissão.

- **Se apenas arquivos `.R` foram modificados**: teste carregando em memória com `pkgload::load_all('.')`:
  ```powershell
  & 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e "pkgload::load_all('.'); # rode testes diretos aqui"
  ```
- **Se arquivos C++ (`src/*.cpp`, `src/*.h`) foram modificados**: tente compilar com `pkgbuild::compile_dll()` ou `R CMD INSTALL --no-multiarch --with-keep.source .`. Se falhar por bloqueio de `pliman.dll`, avise o usuário para reiniciar/fechar a sessão R no RStudio que está usando o pacote antes de recompilar.

## 4. Higiene do Pacote
Se novos arquivos de workspace fora do padrão CRAN forem criados na raiz, confirme que estão listados no `.Rbuildignore`.

