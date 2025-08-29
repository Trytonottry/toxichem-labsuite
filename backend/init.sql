-- init.sql
-- Установка RDKit (через контейнер с rdkit)
-- Примечание: стандартный PostgreSQL не включает RDKit.
-- Нужен образ: dockingbay/rdkit или mcs07/postgresql-rdkit

-- Если используешь образ с RDKit, выполнится автоматически:
CREATE EXTENSION IF NOT EXISTS rdkit;