MAX_ROUNDS = 1

LLM_CONFIG = {
    # "default_model": "llama3.1:8b",
    "default_model": "qwen3:8b",
    # "default_model": "gpt-oss:20b",
    "temperature": {
        "research": 0.3,
        "coding": 0.2,
        "critic": 0.4,
        "execution": 0.1,
        "review": 0.1,
    },
}