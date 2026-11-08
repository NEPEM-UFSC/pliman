lines <- readLines("C:/Users/tiago/.gemini/antigravity-ide/brain/1d5856a2-4fc0-4b36-a08b-c77a309b6255/.system_generated/logs/transcript.jsonl")
user_idx <- grep('"type":"USER_INPUT"', lines)
for (i in user_idx) {
  txt <- lines[i]
  if (grepl("FPS|fps|infer|tempo|12", txt, ignore.case = TRUE)) {
    cat(substr(txt, 1, 300), "\n---\n")
  }
}
