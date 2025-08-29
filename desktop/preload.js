// desktop/preload.js
const { contextBridge, ipcRenderer } = require('electron');

// Безопасный API для renderer-процесса (фронтенда)
contextBridge.exposeInMainWorld('electronAPI', {
  // Пример: отправка асинхронного запроса
  invoke: (channel, data) => {
    // Ограниченный набор каналов
    const validChannels = ['docking:start', 'pbpk:simulate'];
    if (validChannels.includes(channel)) {
      return ipcRenderer.invoke(channel, data);
    }
  },

  // Подписка на события
  on: (channel, func) => {
    const validChannels = ['docking:complete', 'pbpk:updated'];
    if (validChannels.includes(channel)) {
      ipcRenderer.on(channel, (event, ...args) => func(...args));
    }
  },

  // Открытие внешних ссылок
  openExternal: (url) => {
    const allowedDomains = ['pubchem.ncbi.nlm.nih.gov', 'www.rcsb.org', 'kegg.jp'];
    if (allowedDomains.some(domain => url.includes(domain))) {
      ipcRenderer.send('open-external', url);
    }
  },

  // Получение версии приложения
  getAppVersion: async () => {
    return await ipcRenderer.invoke('app:version');
  }
});