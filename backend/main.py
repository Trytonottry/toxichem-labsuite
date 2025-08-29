# backend/main.py (обновлённый)
from fastapi import FastAPI
from api import toxins, docking, predict, pathways, pbpk
import models, database
from fastapi.openapi.docs import get_swagger_ui_html
# backend/main.py — добавь в начало
from logging_config import logger, REQUEST_COUNT, REQUEST_LATENCY
from fastapi.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware
import time

class MetricsMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request, call_next):
        start_time = time.time()
        response = await call_next(request)
        duration = time.time() - start_time

        REQUEST_COUNT.labels(method=request.method, endpoint=request.url.path).inc()
        REQUEST_LATENCY.labels(endpoint=request.url.path).observe(duration)

        logger.info(f"{request.method} {request.url.path} {response.status_code} {duration:.2f}s")
        return response

app.add_middleware(MetricsMiddleware)

app = FastAPI(title="ToxiChem LabSuite", version="2.0")

# Создание таблиц
models.Base.metadata.create_all(bind=database.engine)

app.include_router(toxins.router, prefix="/api", tags=["toxins"])
app.include_router(docking.router, prefix="/api", tags=["docking"])
app.include_router(predict.router, prefix="/api", tags=["predict"])
app.include_router(pathways.router, prefix="/api", tags=["pathways"])
app.include_router(pbpk.router, prefix="/api", tags=["pbpk"])

@app.get("/")
def root():
    return {"message": "Welcome to ToxiChem LabSuite! Docs: /docs"}

@app.get("/docs", include_in_schema=False)
async def custom_swagger_ui():
    return get_swagger_ui_html(
        openapi_url="/openapi.json",
        title="ToxiChem LabSuite API",
        swagger_favicon_url="https://raw.githubusercontent.com/3dmol/3Dmol.js/master/icons/favicon.ico",
        swagger_js_url="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5/swagger-ui-bundle.js",
        swagger_css_url="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5/swagger-ui.css",
        init_oauth={
            "clientId": "toxichem-client",
            "appName": "ToxiChem LabSuite"
        }
    )