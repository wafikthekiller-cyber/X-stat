# X-Solver (lean version)
Run: `pip install -r requirements.txt && uvicorn main:app --host 0.0.0.0 --port 8000` then open http://localhost:8000
Deploy: any host that runs Python (Render, Railway, Fly.io): start command `uvicorn main:app --host 0.0.0.0 --port $PORT`.
Files: solver.py (SymPy engine, input whitelist), main.py (FastAPI + 15 s timeout), static/index.html (UI, KaTeX via CDN).
Input: `2x+3=7`, `x^2-5x+6=0`, inequalities (`<`,`>`,`<=`,`>=`), systems separated by `;`. Unknown must be `x` (or x,y,... in systems).
Not included yet: image OCR, other-variable equations, derivatives/limits/integrals/matrices, PDF file export (use Print > Save as PDF), accounts/DB.
