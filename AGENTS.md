# Diretrizes do Repositório `pliman`

## 1. Limite Crítico de Tamanho por Arquivo (512 KB)
- **REGRA ESTRITA**: Nenhum arquivo `.R` em `R/` pode atingir ou ultrapassar (700 KiB).
- Se um arquivo ultrapassar **modularize novas funções** em arquivos temáticos dedicados em vez de inflar o arquivo existente.
- Após qualquer edição em arquivos grandes em `R/`, verifique o tamanho em bytes e valide o parsing em R.

## 2. Convenções de Imagens, Matrizes e Coordenadas
- Objetos `Image` do `pliman` seguem o padrão `EBImage`: `dim(img)[1]` representa **largura** (eixo X / horizontal) e `dim(img)[2]` representa **altura** (eixo Y / vertical).
- Ao integrar R com C++ (`src/*.cpp`), OpenCV ou tensores ONNX/NCNN (geralmente `linha × coluna` / `H × W`), verifique sempre a transposição e os índices de largura vs. altura antes de finalizar.

## 3. Preservação de Funcionalidades e Estilo de Código
- **Nunca remova funcionalidades ou argumentos existentes** ao otimizar, refatorar ou corrigir uma função (mantenha retrocompatibilidade e aliases quando aplicável).
- Preserve comentários e documentação `roxygen2` existentes que não estejam ligados à mudança atual.
- Use o pacote `cli` (`cli::cli_alert_info()`, `cli::cli_alert_success()`, `cli::cli_alert_warning()`, `cli::cli_bullets()`) para mensagens de console consistentes.

## 4. Padrões de UX, Modelos DL e Pipelines Longos
- **Edição interativa**: Fluxos interativos sobre imagens (ex.: `yolo_dataset_review()`) devem permitir ações contínuas em lote (adicionar/remover múltiplas caixas/pontos em sequência até `<Esc>`, `<Enter>` ou clique direito), evitando retornar ao menu principal a cada item.
- **Gerenciamento de modelos**: Funções que salvam, convertem ou removem modelos (ex.: `pliman_remove_model()`, conversores FP16/TensorRT) devem usar `pliman_model_dir()` como local padrão e aceitar tanto nomes curtos registrados quanto caminhos completos de arquivo.
- **Dependências de hardware**: Ao implementar ou documentar otimizações dependentes de GPU/aceleradores (como FP16 e TensorRT), informe explicitamente na documentação que ganhos de velocidade dependem de hardware compatível (ex.: Tensor Cores), podendo haver perda de desempenho por casting em GPUs sem suporte nativo.
- **Jobs em segundo plano**: Pipelines demorados de dataset ou inferência em lote devem suportar `as_job = TRUE` usando `rstudioapi::jobRunScript()` quando executados no RStudio.

