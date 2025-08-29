// desktop/main.js
const { autoUpdater } = require('electron-updater');

// После создания окна
if (!process.env.ELECTRON_IS_DEV) {
  autoUpdater.checkForUpdatesAndNotify();

  autoUpdater.on('update-available', () => {
    win.webContents.send('update-available');
  });

  autoUpdater.on('update-downloaded', () => {
    win.webContents.send('update-downloaded');
  });
}
const { app, BrowserWindow, Menu } = require('electron');

const { store, db } = require('./store');

// IPC для локальных данных
ipcMain.handle('db:save-toxin', (event, toxin) => {
  return new Promise((resolve, reject) => {
    db.run(
      `INSERT OR REPLACE INTO toxins (name, smiles, data) VALUES (?, ?, ?)`,
      [toxin.name, toxin.smiles, JSON.stringify(toxin)],
      function (err) {
        if (err) reject(err);
        else resolve({ id: this.lastID });
      }
    );
  });
});

ipcMain.handle('db:get-toxins', () => {
  return new Promise((resolve, reject) => {
    db.all(`SELECT * FROM toxins`, (err, rows) => {
      if (err) reject(err);
      else resolve(rows.map(r => ({ ...r, data: JSON.parse(r.data) })));
    });
  });
});

function createWindow () {
  const win = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      nodeIntegration: false
    }
  });

  win.loadURL('http://localhost:3000');
  win.webContents.openDevTools();
}

app.whenReady().then(() => {
  createWindow();
  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow();
  });
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});

// Меню
Menu.setApplicationMenu(Menu.buildFromTemplate([
  { role: 'fileMenu' },
  { role: 'editMenu' },
  { role: 'viewMenu' },
  { role: 'windowMenu' },
  { role: 'help' }
]));