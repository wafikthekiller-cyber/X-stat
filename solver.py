"""X-Solver core: pure SymPy, no network, no eval of raw user code."""
import re
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
from sympy.calculus.util import continuous_domain
from sympy.solvers.inequalities import solve_univariate_inequality

x = sp.Symbol('x', real=True)
ALLOWED = {'sqrt': sp.sqrt, 'ln': sp.log, 'log': lambda a, b=10: sp.log(a, b), 'exp': sp.exp, 'abs': sp.Abs,
           'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'pi': sp.pi, 'e': sp.E}
GLOB = {'Integer': sp.Integer, 'Float': sp.Float, 'Rational': sp.Rational, 'Symbol': sp.Symbol, 'Function': sp.Function}
TR = standard_transformations + (implicit_multiplication_application, convert_xor)
RELS = ['<=', '>=', '<', '>', '=']

class MathError(Exception):
    pass

def parse(s):
    s = s.replace('≤', '<=').replace('≥', '>=').replace('×', '*').replace('÷', '/').replace('π', 'pi').replace('√', 'sqrt')
    if len(s) > 200 or re.search(r'[^0-9a-zA-Z+\-*/^().,\s]', s):
        raise MathError("Caractères non autorisés / رموز غير مسموحة")
    for w in re.findall(r'[a-zA-Z]{2,}', s):
        if w not in ALLOWED:
            raise MathError(f"Fonction ou mot non reconnu : « {w} » / كلمة غير معروفة")
    loc = dict(ALLOWED, x=x)
    try:
        return parse_expr(s, local_dict=loc, global_dict=GLOB, transformations=TR)
    except Exception:
        raise MathError("Expression illisible. Vérifiez les parenthèses. / عبارة غير مفهومة")

def conditions(expr):
    out = []
    for n in sp.preorder_traversal(expr):
        if n.is_Pow and n.has(x):
            b, e = n.args
            if e.is_negative:
                out.append(sp.Ne(b, 0))
            if e.is_Rational and not e.is_integer and e.q % 2 == 0:
                out.append(sp.Ge(b, 0) if not e.is_negative else sp.Gt(b, 0))
        elif isinstance(n, sp.log) and n.args[0].has(x):
            out.append(sp.Gt(n.args[0], 0))
    seen, res = set(), []
    for c in out:
        if c not in seen and c != True:
            seen.add(c); res.append(c)
    return res

def step(fr, ar, tex=None, tfr=None, tar=None):
    return {'fr': fr, 'ar': ar, 'latex': tex}

def plot_points(exprs):
    out = []
    for ex in exprs:
        f = sp.lambdify(x, ex, 'math'); pts = []
        for i in range(-400, 401):
            xv = i / 40
            try:
                v = float(f(xv)); pts.append([xv, v] if abs(v) < 1e6 else None)
            except Exception:
                pts.append(None)
        out.append(pts)
    return out

def method(f):
    if f.is_polynomial(x):
        d = sp.degree(f, x)
        return {0: 'constante', 1: 'équation du 1er degré', 2: 'équation du 2e degré'}.get(d, f'polynôme de degré {d}')
    if f.has(sp.Abs): return 'valeur absolue'
    if f.has(sp.sin, sp.cos, sp.tan): return 'trigonométrique'
    if f.has(sp.log): return 'logarithmique'
    if f.has(sp.exp) or any(n.is_Pow and n.exp.has(x) for n in sp.preorder_traversal(f)): return 'exponentielle'
    if any(n.is_Pow and n.exp.is_Rational and not n.exp.is_integer for n in sp.preorder_traversal(f)): return 'radicale'
    return 'rationnelle'

def solve(text):
    text = text.strip()
    if ';' in text:
        return solve_system(text)
    rel = next((r for r in RELS if r in text), None)
    if rel is None:
        raise MathError("Entrez une équation (=) ou une inéquation (<, >, ≤, ≥). / أدخل معادلة أو متراجحة")
    lhs_s, rhs_s = text.split(rel, 1)
    lhs, rhs = parse(lhs_s), parse(rhs_s)
    syms = (lhs - rhs).free_symbols - {x}
    if syms:
        raise MathError("Une seule inconnue x est supportée dans cette version. / المجهول x فقط")
    f = sp.together(lhs - rhs)
    steps = []
    steps.append({'fr': "Énoncé", 'ar': "المعطى", 'latex': sp.latex(lhs) + {'=':'=','<':'<','>':'>','<=':r'\leq','>=':r'\geq'}[rel] + sp.latex(rhs)})
    conds = conditions(lhs) + conditions(rhs)
    dom = continuous_domain(lhs - rhs, x, sp.S.Reals)
    if conds:
        steps.append({'fr': "Conditions d'existence", 'ar': "نحدد شروط التعريف أولاً حتى لا نقسم على الصفر أو نأخذ جذراً لعدد سالب",
                      'latex': r',\ '.join(sp.latex(c) for c in conds)})
    steps.append({'fr': "Domaine de définition", 'ar': "مجموعة التعريف", 'latex': 'D=' + sp.latex(dom)})
    m = method(f)
    steps.append({'fr': f"Méthode : {m}. On ramène tout dans un membre.", 'ar': "ننقل كل الحدود إلى طرف واحد",
                  'latex': sp.latex(f) + {'=':'=0','<':'<0','>':'>0','<=':r'\leq0','>=':r'\geq0'}[rel]})
    num, den = sp.fraction(f)
    if den != 1:
        steps.append({'fr': "Un quotient est nul si son numérateur est nul (et son dénominateur non nul).",
                      'ar': "الكسر ينعدم عندما ينعدم بسطه ولا ينعدم مقامه", 'latex': sp.latex(sp.factor(num)) + '=0'})
    if rel == '=' and f.is_polynomial(x) and sp.degree(f, x) == 2:
        a, b, c = sp.Poly(f, x).all_coeffs()
        steps.append({'fr': "On calcule le discriminant Δ = b² − 4ac.", 'ar': "نحسب المميز Δ",
                      'latex': rf'\Delta={sp.latex(b**2-4*a*c)}'})
    elif rel == '=' and f.is_polynomial(x):
        steps.append({'fr': "Factorisation.", 'ar': "نحلل العبارة إلى جداء عوامل", 'latex': sp.latex(sp.factor(f)) + '=0'})
    # solving
    approx = False
    try:
        if rel == '=':
            sol = sp.solveset(f, x, domain=sp.S.Reals)
            if isinstance(sol, sp.ConditionSet):
                raise MathError("Pas de solution exacte trouvée (équation non supportée). / لم يتم إيجاد حل دقيق")
        else:
            r = {'<': f < 0, '>': f > 0, '<=': f <= 0, '>=': f >= 0}[rel]
            sol = solve_univariate_inequality(r, x, relational=False, domain=sp.S.Reals)
    except MathError:
        raise
    except Exception:
        raise MathError("Résolution impossible pour cette expression. / تعذر الحل")
    steps.append({'fr': "Résolution dans ℝ", 'ar': "نحل في مجموعة الأعداد الحقيقية", 'latex': sp.latex(sol)})
    final = sp.Intersection(sol, dom)
    if final != sol:
        removed = sp.Complement(sol, dom)
        steps.append({'fr': "On élimine les valeurs interdites (solutions étrangères).", 'ar': "نستبعد القيم الممنوعة (حلول دخيلة)", 'latex': sp.latex(removed)})
    # verification
    ver = []
    if rel == '=' and final.is_FiniteSet:
        for s in final:
            d = sp.simplify(lhs.subs(x, s) - rhs.subs(x, s))
            ok = d == 0 or abs(complex(sp.N(d, 30))) < 1e-15
            ver.append(rf'x={sp.latex(s)}:\ {sp.latex(sp.simplify(lhs.subs(x, s)))}={sp.latex(sp.simplify(rhs.subs(x, s)))}\ \checkmark' if ok else rf'x={sp.latex(s)}\ \text{{rejeté}}')
            if not ok: final = final - {s}
    elif rel != '=':
        for t in (sp.Float(-7.31), sp.Float(0.57), sp.Float(4.13)):
            if t in dom:
                try:
                    holds = bool(sp.simplify({'<': lhs < rhs, '>': lhs > rhs, '<=': lhs <= rhs, '>=': lhs >= rhs}[rel].subs(x, t)))
                    ver.append(rf'x={t}:\ {"vrai" if holds else "faux"}\ ({"∈" if t in final else "∉"}\ S)')
                except Exception: pass
    if ver:
        steps.append({'fr': "Vérification dans l'équation initiale" if rel == '=' else "Vérification par valeurs test", 'ar': "نتحقق بالتعويض في المعادلة الأصلية", 'latex': r'\\ '.join(ver)})
    if final.is_FiniteSet and any(not s.is_rational for s in final):
        approx_vals = ', '.join(str(sp.N(s, 6)) for s in final)
    else:
        approx_vals = None
    return {'kind': 'equation' if rel == '=' else 'inequation', 'steps': steps,
            'solution': ('S=' + sp.latex(final)) if final != sp.S.EmptySet else r'S=\varnothing', 'approx': approx_vals,
            'graph': plot_points([lhs, rhs]) if rel else None}

def solve_system(text):
    eqs = []
    for part in text.split(';'):
        if '=' not in part or any(r in part for r in ('<', '>')):
            raise MathError("Système : équations séparées par « ; ». / النظام: معادلات مفصولة بـ ;")
        a, b = part.split('=', 1)
        eqs.append(sp.Eq(parse(a), parse(b)))
    vars_ = sorted(set().union(*[e.free_symbols for e in eqs]), key=str)
    if not vars_ or len(vars_) > 4:
        raise MathError("Variables invalides. / متغيرات غير صالحة")
    sol = sp.solve(eqs, vars_, dict=True)
    steps = [{'fr': "Système", 'ar': "الجملة", 'latex': r'\begin{cases}' + r'\\'.join(sp.latex(e) for e in eqs) + r'\end{cases}'},
             {'fr': "Résolution (substitution/élimination)", 'ar': "نحل بالتعويض أو بالحذف", 'latex': sp.latex(sol)}]
    ver = [r'\\ '.join(sp.latex(sp.simplify(e.lhs.subs(s) - e.rhs.subs(s))) + '=0\\ \\checkmark' for e in eqs) for s in sol]
    if ver:
        steps.append({'fr': "Vérification", 'ar': "التحقق بالتعويض", 'latex': r'\\ '.join(ver)})
    return {'kind': 'system', 'steps': steps, 'solution': 'S=' + (sp.latex(sol) if sol else r'\varnothing'), 'approx': None, 'graph': None}
