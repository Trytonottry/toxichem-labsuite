// frontend/src/env.js
export const isElectron = () => !!window.electronAPI;
export const isPWA = () => window.matchMedia('(display-mode: standalone)').matches;