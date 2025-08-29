// frontend/vite.config.js
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig(({ mode }) => {
  return {
    plugins: [react()],
    define: {
      __IS_ELECTRON__: mode === 'electron',
      __IS_PWA__: mode === 'pwa'
    },
    build: {
      outDir: mode === 'electron' ? '../desktop/frontend/dist' : 'dist'
    }
  };
});