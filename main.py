from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from solver import solve, MathError
import concurrent.futures

app = FastAPI(title="X-Solver")
pool = concurrent.futures.ProcessPoolExecutor(max_workers=2)

class Q(BaseModel):
    expr: str

@app.post("/api/solve")
def api_solve(q: Q):
    try:
        return pool.submit(solve, q.expr).result(timeout=15)  # time limit protects the server
    except MathError as e:
        raise HTTPException(400, str(e))
    except concurrent.futures.TimeoutError:
        raise HTTPException(408, "Calcul trop long / العملية طويلة جداً")

@app.get("/")
def index():
    return FileResponse("static/index.html")

app.mount("/static", StaticFiles(directory="static"), name="static")
