# backend/logging_config.py
import logging
import sys
from loguru import logger

# Удаляем стандартный логгер
logger.remove()

# Логирование в stdout
logger.add(sys.stdout, level="INFO", format="{time} {level} {message}")

# Логирование в файл
logger.add("logs/app.log", rotation="100 MB", level="DEBUG")

# Sentry (ошибки в продакшене)
try:
    import sentry_sdk
    sentry_sdk.init(
        dsn="https://your-sentry-dsn@app.sentry.io/project-id",
        traces_sample_rate=0.5,
        profiles_sample_rate=0.5,
    )
    logger.add(lambda msg: sentry_sdk.capture_message(msg), level="ERROR")
except ImportError:
    pass

# Prometheus (метрики)
from prometheus_client import start_http_server, Counter, Histogram

start_http_server(8001)  # /metrics

REQUEST_COUNT = Counter('http_requests_total', 'Total HTTP Requests', ['method', 'endpoint'])
REQUEST_LATENCY = Histogram('request_latency_seconds', 'Request latency', ['endpoint'])