// desktop/store.js
const Store = require('electron-store');
const path = require('path');
const fs = require('fs');

const schema = {
  lastProject: { type: 'string', nullable: true },
  recentToxins: {
    type: 'array',
    items: { type: 'string' },
    default: []
  },
  settings: {
    type: 'object',
    properties: {
      theme: { type: 'string', enum: ['light', 'dark'], default: 'dark' },
      language: { type: 'string', default: 'en' },
      offlineMode: { type: 'boolean', default: false }
    }
  }
};

const store = new Store({ schema });

// Локальная SQLite для молекул (если не запущен backend)
const dbPath = path.join(store.path, '../data/local.db');
const sqlite3 = require('sqlite3').verbose();
const db = new sqlite3.Database(dbPath);

db.serialize(() => {
  db.run(`
    CREATE TABLE IF NOT EXISTS toxins (
      id INTEGER PRIMARY KEY,
      name TEXT,
      smiles TEXT UNIQUE,
      data TEXT,
      created_at DATETIME DEFAULT CURRENT_TIMESTAMP
    )
  `);
});

module.exports = { store, db };