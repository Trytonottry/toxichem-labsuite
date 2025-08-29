# backend/workers/celery_app.py
from celery import Celery
from .tasks import run_docking_task, simulate_pbk_task

celery_app = Celery(
    'toxichem',
    broker='redis://redis:6379',
    backend='redis://redis:6379'
)

celery_app.conf.update(
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
)

# Регистрация задач
celery_app.autodiscover_tasks(['workers.tasks'])