import os
import numpy as np


# 1) Lecture robuste d'un fichier 2 colonnes (T, valeur)

def load_two_columns(path):
    """
    Charge un fichier 2 colonnes numériques (T, val).
    - Ignore lignes vides et commentaires (#).
    - Accepte ',' ou ';' comme séparateur.
    - Tolère une éventuelle ligne d'en-tête (si non-numérique).
    Retourne (T_sorted, Y_sorted) triés par T.
    """
    T, Y = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            s = s.replace(";", ",")
            parts = [p.strip() for p in s.split(",") if p.strip() != ""]
            if len(parts) < 2:
                parts = s.split()
            try:
                t = float(parts[0]); y = float(parts[1])
            except Exception:
                continue
            T.append(t); Y.append(y)

    if not T or len(T) != len(Y):
        raise ValueError(f"Fichier {path} vide/mal formé (attendu: 2 colonnes numériques).")

    T = np.asarray(T, float)
    Y = np.asarray(Y, float)
    idx = np.argsort(T)
    return T[idx], Y[idx]


# 1bis) Lecture d'un .txt contenant une valeur constante (rho)

def load_constant_txt(path):
    """
    Lit la première valeur numérique d'un .txt (ignore lignes vides/commentées).
    Retourne une 'table' plate (T=[0,1], Y=[rho,rho]) pour compatibilité avec np.interp.
    """
    val = None
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # prend le 1er token numérique rencontré
            for tok in s.replace(";", " ").replace(",", " ").split():
                try:
                    val = float(tok)
                    break
                except Exception:
                    continue
            if val is not None:
                break
    if val is None:
        raise ValueError(f"Aucune valeur numérique trouvée dans {path}")
    # Table 2 points pour que l'interpolation marche (saturation identique aux bornes)
    return np.array([0.0, 1.0], float), np.array([val, val], float)

# ----------------------------------------
# 2) Chargement des tables pour tous les matériaux présents
# ----------------------------------------
def chargement_tables_caracteristique(index_materiaux, dossier="data"):
    """
    index_materiaux: matrice 2 x N (ligne 0 = x/indices, ligne 1 = noms matériaux)
    Cherche:
      - k_<mat>.csv
      - Cp_<mat>.csv
      - rho_<mat>.csv OU rho_<mat>.txt
    Retourne: tables_k, tables_cp, tables_rho (dict {mat: (Ttab, Ytab)})
    """
    if index_materiaux.shape[0] < 2:
        raise ValueError("index_materiaux doit avoir au moins 2 lignes (x, matériau).")

    noms = np.asarray(index_materiaux[-1, :], dtype=object).astype(str)
    mats_uniques = np.unique(noms)

    def _path_first_existing(bases, mat, exts):
        """
        bases: liste de préfixes, ex: ["k_", "K_"] si besoin (ici 1 suffit)
        exts : liste d'extensions, ex: [".csv"] ou [".csv",".txt"]
        Essaie mat et mat.lower().
        """
        candidates = []
        for base in bases:
            for name in (mat, mat.lower()):
                for ext in exts:
                    candidates.append(os.path.join(dossier, f"{base}{name}{ext}"))
        for p in candidates:
            if os.path.isfile(p):
                return p
        return None

    tables_k, tables_cp, tables_rho = {}, {}, {}
    manquants = []

    for mat in mats_uniques:
        p_k   = _path_first_existing(["k_"],   mat, [".csv"])
        p_cp  = _path_first_existing(["Cp_"],  mat, [".csv"])
        # rho accepte .csv OU .txt
        p_rho_csv = _path_first_existing(["rho_"], mat, [".csv"])
        p_rho_txt = _path_first_existing(["rho_"], mat, [".txt"])
        p_rho = p_rho_csv if p_rho_csv is not None else p_rho_txt

        if p_k is None:      manquants.append(f"k_{mat}.csv")
        if p_cp is None:     manquants.append(f"Cp_{mat}.csv")
        if p_rho is None:    manquants.append(f"rho_{mat}.csv|txt")

        if (p_k is not None) and (p_cp is not None) and (p_rho is not None):
            tables_k[mat]  = load_two_columns(p_k)
            tables_cp[mat] = load_two_columns(p_cp)
            # rho: csv -> table, txt -> constante
            if p_rho.endswith(".txt"):
                tables_rho[mat] = load_constant_txt(p_rho)
            else:
                tables_rho[mat] = load_two_columns(p_rho)

    if manquants:
        # On ne lève l'erreur que si au moins un des 3 manque pour au moins un matériau
        # (et qu'on n'a pas réussi à tout charger pour tous)
        # On vérifie que tous les mats ont bien été ajoutés:
        if (set(mats_uniques) - set(tables_k.keys())) \
           or (set(mats_uniques) - set(tables_cp.keys())) \
           or (set(mats_uniques) - set(tables_rho.keys())):
            raise FileNotFoundError(
                f"Fichiers de tables manquants dans '{dossier}': {', '.join(manquants)}"
            )

    return tables_k, tables_cp, tables_rho

def interp_table(tables_prop, Tvals, mats):
    """
    Interpole une propriété (k, Cp ou rho) cellule par cellule en fonction du matériau.
    - tables_prop : dict {mat: (Ttab, Ytab)}   (ex: "AFRSI": (T_k, k_k))
    - Tvals       : array (Nx,)  températures par cellule
    - mats        : array (Nx,)  noms de matériau par cellule (strings)

    -> Retourne array (Nx,) de la propriété interpolée (linéaire, saturée aux bornes).
    """
    Tvals = np.asarray(Tvals, dtype=float)
    mats  = np.asarray(mats,  dtype=object).astype(str)

    if Tvals.shape != mats.shape:
        raise ValueError("Tvals et mats doivent avoir la même forme.")

    out = np.empty_like(Tvals, dtype=float)

    # dictionnaire pour tolérer la casse (lower-case lookup)
    lower_key = {k.lower(): k for k in tables_prop.keys()}

    for m in np.unique(mats):
        key = m if m in tables_prop else lower_key.get(m.lower(), None)
        if key is None:
            raise KeyError(f"Aucune table trouvée pour le matériau '{m}'.")
        Ttab, Ytab = tables_prop[key]
        mask = (mats == m)
        out[mask] = np.interp(Tvals[mask], Ttab, Ytab, left=Ytab[0], right=Ytab[-1])

    return out


# run_thermal_design.py
# Python minimal (sans classes) : Fourier analytique + FVM T-dépendant + tableau de résultats

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# 1) Fourier analytique (mur plan, gauche T=Ts, droite adiabatique)


def fourier_Tcentre_adiabatique(e, t, k, rho, Cp, Ts, Ti, nmax=500, tol=1e-12):
    """
    Température au point x=e (extrémité adiabatique == centre par symétrie) à l'instant t.
    Série : θ(L,t) = 2 θ_i * Σ_{n=0..∞} [ (-1)^n / ((2n+1)π) * exp(-α * ((2n+1)π/(2L))^2 * t) ]
    où θ = T - Ts, L = e, α = k/(ρCp), θ_i = Ti - Ts.
    """
    L = float(e)
    alpha = float(k) / (float(rho) * float(Cp))
    theta_i = float(Ti) - float(Ts)

    s = 0.0
    for n in range(nmax):
        m = 2*n + 1
        lam = (m*np.pi) / (2.0*L)        # (1/m)
        term = ((-1.0)**n) / (m*np.pi) * np.exp(-alpha * (lam**2) * t)
        s_new = s + term
        # critère d'arrêt
        if abs(term) < tol*max(1.0, abs(s_new)):
            s = s_new
            break
        s = s_new
    return float(Ts + 4.0 * theta_i * s)

def fourier_profile_adiabatique(e, t, k, rho, Cp, Ts, Ti, x, nmax=600):
    """
    Retourne T(x,t) (shape = (Nx,)) pour Dirichlet à x=0 (Ts) et Neumann à x=e.
    Formule série:
      T(x,t) = Ts + (4/π)*(Ti-Ts) * Σ_{n=0..∞} [ sin((2n+1)π x/(2e)) / (2n+1)
                                                * exp(-α ((2n+1)π/(2e))^2 t) ]
    """
    x = np.asarray(x, float)                 # (Nx,)
    n = np.arange(nmax)                      # (nmax,)
    m = 2*n + 1                              # (nmax,)
    alpha = k/(rho*Cp)
    lam   = (m*np.pi)/(2.0*e)                # (nmax,)
    coeff = (4.0*(Ti - Ts)/np.pi) * (1.0/m) * np.exp(-alpha*(lam**2)*t)  # (nmax,)
    S     = np.sin(np.outer(lam, x))         # (nmax, Nx)
    theta = (coeff[:, None] * S).sum(axis=0) # (Nx,)
    return Ts + theta                        # (Nx,)

def dichotomie_epaisseur(T_target, t, k, rho, Cp, Ts, Ti, e_lo, e_hi, tol=1e-5, max_iter=80):
    """
    Bisection sur e avec la solution analytique. On cherche Tcentre(e) - T_target <= 0.
    Hypothèse vraie ici : Tcentre(e) décroît quand e augmente.
    """
    def f(e):
        return fourier_Tcentre_adiabatique(e, t, k, rho, Cp, Ts, Ti) - T_target

    f_lo = f(e_lo)
    f_hi = f(e_hi)
    # élargit la borne haute si nécessaire
    it_guard = 0
    while f_hi > 0.0 and e_hi < 1.0 and it_guard < 40:
        e_hi *= 1.5
        f_hi = f(e_hi)
        it_guard += 1

    if f_lo <= 0.0:
        return e_lo
    if f_hi > 0.0:
        raise RuntimeError("Impossible de borner la solution : augmentez e_hi.")

    for _ in range(max_iter):
        e_mid = 0.5*(e_lo + e_hi)
        fm = f(e_mid)
        if abs(e_hi - e_lo) < tol:
            return e_mid
        if fm > 0.0:
            e_lo = e_mid
        else:
            e_hi = e_mid
    return 0.5*(e_lo + e_hi)

def mean(a, b):
    a = float(a); b = float(b)
    return (a+b)/2

def Initialisation(e, Nx, Tinitiale, index_Materiaux, tmax, dt):
    x = np.linspace(0.0, e, Nx)
    dx = e/Nx;
    To = Tinitiale*np.ones(Nx)
    nt = int(np.ceil(tmax / dt))
    T_xt = np.empty((nt + 1, Nx), dtype=float)
    T_xt[0, :] = To
    materiaux = Construction_Vecteur_Materiaux(index_Materiaux,Nx)
    return x,dx,T_xt,To,materiaux
    
def Construction_Vecteur_Materiaux(index_Materiaux, Nx):
    """
    Construit le vecteur des matériaux pour chaque cellule du maillage 1D.
    index_Materiaux : matrice 2xN, ex:
        [[0.1, 0.2],
         ["AFRSI", "FRCI12"]]
      => première ligne : positions limites en mètres
         deuxième ligne : noms des matériaux
    Nx : nombre total de cellules
    Retour : np.array (Nx,) de strings (matériau par cellule)
    """
    # coordonnées de centres de cellules
    e_total = float(index_Materiaux[0, -1])
    x_cells = np.linspace(0, e_total, Nx)
    materiaux = np.empty(Nx, dtype=object)

    # boucle sur les couches
    x_prev = 0.0
    for i in range(index_Materiaux.shape[1]):
        x_lim = float(index_Materiaux[0, i])
        nom = str(index_Materiaux[1, i])
        mask = (x_cells >= x_prev) & (x_cells <= x_lim)
        materiaux[mask] = nom
        x_prev = x_lim

    return materiaux

    
def ResolutionVolumesFinis(x, dx, Nx, dt, tmax, To, T_xt, materiaux,
                           méthode, k_defaut, rho_defaut, Cp_defaut,
                           CLdroit, ValCLDroit, CLgauche, ValClgauche):
    """
    Schéma FV explicite 1D.
    - T_xt est une matrice 2D (nt+1, Nx) avec T_xt[0,:] = To déjà posé.
    - méthode:
        * "constant": k_defaut/Cp_defaut/rho_defaut = dict {mat: valeur}
        * "dependant": k_defaut/Cp_defaut/rho_defaut = dict {mat: (Ttab,Ytab)}, utilisé via interp_table
    - CLgauche / CLdroit: "flux" ou "Temp"
    """
    nt = int(tmax / dt)

    # Préparer les "tables" si dépendant
    if méthode == "dependant":
        tables_k   = k_defaut
        tables_Cp  = Cp_defaut
        tables_rho = rho_defaut

    for i in range(nt):
        # État courant
        T = T_xt[i, :]

        # --- Propriétés au pas i (vecteurs Nx) ---
        if méthode == "constant":
            k_de_x   = np.array([k_defaut[str(m)]   for m in materiaux], dtype=float)
            Cp_de_x  = np.array([Cp_defaut[str(m)]  for m in materiaux], dtype=float)
            rho_de_x = np.array([rho_defaut[str(m)] for m in materiaux], dtype=float)
        else:  # "dependant"
            k_de_x   = interp_table(tables_k,   T, materiaux)
            Cp_de_x  = interp_table(tables_Cp,  T, materiaux)
            rho_de_x = interp_table(tables_rho, T, materiaux)

        # --- Boucle spatiale (indices corrigés) ---
        vecteurPassageT = np.empty_like(T)
        for j in range(Nx):
            # Bord gauche (j = 0) : ne jamais utiliser j-1 ici
            if j == 0:
                if CLgauche == "flux":
                    D = 0.5*(k_de_x[j+1] + k_de_x[j]) * (T[j+1] - T[j])/dx
                    G = -ValClgauche
                elif CLgauche == "Temp":
                    D = 0.5*(k_de_x[j+1] + k_de_x[j]) * (T[j+1] - T[j])/dx
                    G = -1*(k_de_x[j+1] + k_de_x[j]) * (T[j] - ValClgauche)/dx
                else:
                    raise ValueError("CLgauche doit être 'flux' ou 'Temp'.")

            # Bord droit (j = Nx-1) : j+1 n'existe pas
            elif j == Nx-1:
                if CLdroit == "flux":
                    D = -ValCLDroit
                    G = -0.5*(k_de_x[j] + k_de_x[j-1]) * (T[j] - T[j-1])/dx
                elif CLdroit == "Temp":
                    D = -1*(k_de_x[j] + k_de_x[j-1]) * (T[j] - ValCLDroit)/dx
                    G = -0.5*(k_de_x[j] + k_de_x[j-1]) * (T[j] - T[j-1])/dx
                else:
                    raise ValueError("CLdroit doit être 'flux' ou 'Temp'.")

            # Intérieur (1 .. Nx-2)
            else:
                D = 0.5*(k_de_x[j+1] + k_de_x[j]) * (T[j+1] - T[j])/dx
                G = -0.5*(k_de_x[j] + k_de_x[j-1]) * (T[j]   - T[j-1])/dx

            # Mise à jour cellule j (indices sur j, pas j-1)
            vecteurPassageT[j] = T[j] + dt/(rho_de_x[j]*Cp_de_x[j]*dx) * (D + G)

        # Écrire la ligne suivante (et pas la même)
        T_xt[i+1, :] = vecteurPassageT

    return T_xt


    
def thermal_stress_profile(x, T_x, alpha, E, nu, T_ref):
    """
    Contrainte thermo-élastique 1D sur l'épaisseur avec :
      - face externe (x=0) libre :      sigma(0,t) = 0
      - face arrière (x=h) bloquée :     epsilon(h,t) = 0
    Modèle membrane/poutre avec déformation linéaire : eps(x)=eps0 + kappa*(x - h/2)
    Loi : sigma = Ebar*(eps - eps_theta), avec eps_theta = alpha*(T - T_ref)
    Ebar = E/(1 - nu) (choix conservatif pour une tuile collée; tu peux utiliser E si souhaité).

    Paramètres
    ----------
    x : (N,) abscisses croissantes sur [0,h]
    T_x : (N,) profil de température (K) au temps considéré
    alpha : (1/K) CTE linéaire
    E : (Pa) module d'Young (thru-thickness)
    nu : coefficient de Poisson (adopte 0.1-0.2 si inconnu)
    T_ref : (K) température de référence des dilatations
    """
    x = np.asarray(x, float)
    T_x = np.asarray(T_x, float)
    h = float(x[-1] - x[0])

    # E "effectif" (voir commentaires ci-dessus)
    Ebar = E/(1.0 - nu)

    # Déformation thermique locale
    eps_th = alpha*(T_x - T_ref)

    # CLs → eps0 et kappa :
    #   sigma(0)=0  ⇒  eps(0)=eps_th(0)
    #   eps(h)=0
    # eps(x)=eps0 + kappa*(x - h/2)
    eps0   = 0.5*(eps_th[0] + eps_th[-1])
    kappa  = (eps_th[0] - eps_th[-1]) / h

    # Contrainte sur l'épaisseur
    sigma = Ebar*(eps0 + kappa*(x - 0.5*h) - eps_th)

    return sigma  # (Pa)