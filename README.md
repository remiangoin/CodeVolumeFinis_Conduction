#  CodeVolumeFinis_Conduction

Ce projet traite de la **conduction thermique généralisée en 1D** dans un matériau composite, avec prise en compte des **propriétés dépendantes de la température** (conductivité, capacité thermique, densité).  
Il compare les résultats d’un **modèle analytique (Fourier)** et d’une **méthode numérique par volumes finis (FVM)**, puis calcule les **contraintes thermo-élastiques internes**.

---

##  Structure du dépôt

|
|-- Fonctions.py # Fonctions de lecture et d'interpolation des propriétés matériaux
|-- main.py # Programme principal : simulation et post-traitement
|-- data/ # Dossier à créer contenant les propriétés thermiques
| |-- k_<mat>.csv # Conductivité thermique (W/m/K)
| |-- Cp_<mat>.csv # Capacité thermique (J/kg/K)
| |-- rho_<mat>.txt # Densité (kg/m³)
|-- README.md 


---

##  Fonctionnalités principales

### 1. Lecture des données matériaux
- Lecture robuste de fichiers à **2 colonnes (T, valeur)** avec `,` ou `;` comme séparateur.  
- Gestion des **lignes commentées (#)** et des **fichiers vides/mal formés**.  
- Lecture de valeurs constantes depuis un `.txt` pour la densité.  
- Chargement automatique des fichiers `k_`, `Cp_` et `rho_` pour chaque matériau.

### 2. Modèles de conduction
- **Modèle analytique (Fourier)** : conduction transitoire dans un mur plan avec condition adiabatique à droite et température imposée à gauche.  
- **Méthode numérique (FVM 1D explicite)** :  
  - Gestion des conditions aux limites (flux ou température).  
  - Propriétés thermiques constantes ou dépendantes de la température.  
  - Schéma explicite simple, lisible et modulaire.

### 3. Post-traitement
- Tracé des profils de température à différents instants.  
- Comparaison entre FVM constant, FVM variable et solution analytique.  
- Calcul et affichage de l’erreur relative (%) entre les modèles.  
- Calcul des **contraintes thermo-élastiques internes** au cours du temps.

---

##  Paramètres principaux

Les paramètres définis dans `main.py` contrôlent la simulation :

| Paramètre | Signification | Valeur par défaut |
|------------|----------------|-------------------|
| `e` | Épaisseur totale du matériau (m) | 0.10 |
| `Nx` | Nombre de volumes finis | 200 |
| `dt` | Pas de temps (s) | 0.02 |
| `tmax` | Durée totale (s) | 900 |
| `Tinitiale` | Température initiale uniforme (K) | 273.15 |
| `CLgauche` / `CLdroit` | Type de condition limite ("Temp" ou "flux") | `"Temp"` / `"flux"` |
| `ValClgauche` / `ValCLDroit` | Valeurs associées (K ou W/m²) | `1000.0` / `0.0` |

---

##  Exemple de configuration

```python
index_materiaux = np.array([
    [0.1],          # Limite de la couche (m)
    ["AFRSI"]       # Nom du matériau
], dtype=object)

## Dépendance
pip install numpy scipy matplotlib





