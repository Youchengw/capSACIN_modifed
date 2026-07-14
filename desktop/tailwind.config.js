/** @type {import('tailwindcss').Config} */
export default {
  content: ["./index.html", "./src/**/*.{js,ts,jsx,tsx}"],
  theme: {
    extend: {
      colors: {
        capsacin: {
          blue: "#1e3a5f",
          light: "#e8f0fe",
          accent: "#3b82f6",
          gold: "#eab308",
          green: "#22c55e",
        },
      },
    },
  },
  plugins: [],
};
