import Fonctions

def main():

    # 1) Données physiques et numériques

    e = 0.10                 # épaisseur (m) = 10 cm
    Nx = 200                 # nombre de volumes finis
    dt = 0.02                 # pas de temps (s)
    tmax = 900.0             # durée totale = 15 min
    Tinitiale = 273.15       # température initiale uniforme (K)


    # Propriétés constantes pour l'analytique (cohérentes FRCI-12 / placeholders à ajuster si tu as mieux)

    k_defaut   = {"AFRSI": 0.049} # W/m/K
    rho_defaut = {"AFRSI": 192.0} # kg/m^3
    Cp_defaut  = {"AFRSI": 711.0} # J/kg/K

    # Paramètres mécaniques pour contrainte (ordre de grandeur FRCI-12)
    alpha = 5e-7          # 1/K (CTE moyen demandé)
    E     = 1.0e8         # Pa  (module Young thru-thickness ~1e8 Pa)
    nu    = 0.15          # (-)  (Poisson ; ajuste si tu as une valeur)


    # 2) Définition du matériau et de sa matrice

    # 1 seule couche : AFRSI sur 10 cm
    index_materiaux = np.array([
        [0.1],          # limite de la couche (0.1 m)
        ["AFRSI"]       # matériau
    ], dtype=object)

    methode = "dependant"     # dépendant (propriétés dépendantes de T) ou "constant" (rentrées dans "valeur"_defaut)


    # 3) Conditions aux limites

    CLgauche = "Temp"         # température imposée
    ValClgauche = 1000.0      # K
    CLdroit = "flux"          # flux imposé (Neumann)
    ValCLDroit = 0.0          # adiabatique


    # 4) Chargement des tables matériaux

    # le dossier doit contenir k_AFRSI.csv, Cp_AFRSI.csv, rho_AFRSI.txt
    tables_k, tables_cp, tables_rho = chargement_tables_caracteristique(index_materiaux, dossier="data")


    # 5) Initialisation du domaine

    x, dx, T_xt, To, materiaux = Initialisation(
    e=e, Nx=Nx, Tinitiale=Tinitiale,
    index_Materiaux=index_materiaux,
    tmax=tmax, dt=dt
)



    # 6) Lancement des solutions

    T_xt_FVM_variable = ResolutionVolumesFinis(
        x, dx, Nx, dt, tmax,
        To, T_xt.copy(), materiaux,
        methode,
        k_defaut=tables_k,        # pour compatibilité du prototype, même si non utilisé ici
        rho_defaut=tables_rho,
        Cp_defaut=tables_cp,
        CLdroit=CLdroit, ValCLDroit=ValCLDroit,
        CLgauche=CLgauche, ValClgauche=ValClgauche
    )

    T_xt_FVM_constant = ResolutionVolumesFinis(
        x, dx, Nx, dt, tmax,
        To, T_xt.copy(), materiaux,
        "constant",
        k_defaut,        # pour compatibilité du prototype, même si non utilisé ici
        rho_defaut,
        Cp_defaut,
        CLdroit=CLdroit, ValCLDroit=ValCLDroit,
        CLgauche=CLgauche, ValClgauche=ValClgauche
    )
    Nt = int(np.ceil(tmax / dt))
    T_xt_analytique = np.zeros((Nt+1, Nx))   # (temps, espace)
    
    for j in range(Nt+1):
        tcur = j * dt
        # Profil analytique à cet instant tcur
        T_xt_analytique[j, :] = fourier_profile_adiabatique(
            e, tcur,
            k_defaut["AFRSI"], rho_defaut["AFRSI"], Cp_defaut["AFRSI"],
            ValClgauche, Tinitiale, x,
            nmax=500
        )

    erreur_FVM_constant = ((T_xt_FVM_constant - T_xt_analytique)/T_xt_analytique)*100
    erreur_FVM_variable = ((T_xt_FVM_variable - T_xt_analytique)/T_xt_analytique)*100

    sigma1 = thermal_stress_profile(x, T_xt_FVM_variable[int(tmax/(3*dt)), :], alpha, E, nu, Tinitiale)
    sigma2 = thermal_stress_profile(x, T_xt_FVM_variable[int(2*tmax/(3*dt)), :], alpha, E, nu, Tinitiale)
    sigma3 = thermal_stress_profile(x, T_xt_FVM_variable[int(tmax/dt), :], alpha, E, nu, Tinitiale)
    

    # 7) Post-traitement

    # Affichage de la température pour t = tmax/3, 2*tmax/3 et tmax sur une seule ligne
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    
    axes[0].plot(x, T_xt_FVM_variable[int(tmax/(3*dt))], label=f"T(t={tmax/3:.0f}s), FVM variable")
    axes[0].plot(x, T_xt_FVM_constant[int(tmax/(3*dt))],"--", label=f"T(t={tmax/3:.0f}s) FVM constant")
    axes[0].plot(x, T_xt_analytique[int(tmax/(3*dt))],":", label=f"T(t={tmax/3:.0f}s) FVM analytique")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("Température (K)")
    axes[0].legend()
    axes[0].grid(True)
    
    axes[1].plot(x, T_xt_FVM_variable[int(2*tmax/(3*dt))], label=f"T(t={2*tmax/3:.0f}s) FVM variable")
    axes[1].plot(x, T_xt_FVM_constant[int(2*tmax/(3*dt))], label=f"T(t={2*tmax/3:.0f}s) FVM constant")
    axes[1].plot(x, T_xt_analytique[int(2*tmax/(3*dt))], label=f"T(t={2*tmax/3:.0f}s) FVM analytique")
    axes[1].set_xlabel("x (m)")
    axes[1].legend()
    axes[1].grid(True)
    
    axes[2].plot(x, T_xt_FVM_variable[-1], label=f"T(t={tmax:.0f}s) FVM variable")
    axes[2].plot(x, T_xt_FVM_constant[-1], label=f"T(t={tmax:.0f}s) FVM constant")
    axes[2].plot(x, T_xt_analytique[-1], label=f"T(t={tmax:.0f}s) FVM analytique")
    axes[2].set_xlabel("x (m)")
    axes[2].legend()
    axes[2].grid(True)
    
    fig.suptitle("Profils de température", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()



    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    
    axes[0].plot(x, erreur_FVM_constant[int(tmax/(3*dt))], label=f"t={tmax/3:.0f}s, erreur constant")
    axes[0].plot(x, erreur_FVM_variable[int(tmax/(3*dt))],"--", label=f"t={tmax/3:.0f}s, erreur variable")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("erreur (%)")
    axes[0].legend()
    axes[0].grid(True)
    
    axes[1].plot(x, erreur_FVM_constant[int(2*tmax/(3*dt))], label=f"t={2*tmax/3:.0f}s, erreur constant")
    axes[1].plot(x, erreur_FVM_variable[int(2*tmax/(3*dt))], label=f"t={2*tmax/3:.0f}s, erreur variable")
    axes[1].set_xlabel("x (m)")
    axes[1].legend()
    axes[1].grid(True)
    
    axes[2].plot(x, erreur_FVM_constant[-1], label=f"t={tmax:.0f}s erreur constant")
    axes[2].plot(x, erreur_FVM_variable[-1], label=f"t={tmax:.0f}s erreur variable")
    axes[2].set_xlabel("x (m)")
    axes[2].legend()
    axes[2].grid(True)
    
    fig.suptitle("Erreur par rapport à l'analytique", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    
    axes[0].plot(x, sigma1, label=f"t={tmax/3:.0f}s, contrainte interne (cas FVM_variable)")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("Contrainte (Pa)")
    axes[0].legend()
    axes[0].grid(True)
    
    axes[1].plot(x, sigma2, label=f"t={2*tmax/3:.0f}s, contrainte interne (cas FVM_variable)")
    axes[1].set_xlabel("x (m)")
    axes[1].legend()
    axes[1].grid(True)
    
    axes[2].plot(x, sigma3, label=f"t={tmax:.0f}s contrainte interne (cas FVM_variable)")
    axes[2].set_xlabel("x (m)")
    axes[2].legend()
    axes[2].grid(True)
    
    fig.suptitle("Contraintes internes", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    
if __name__ == "__main__":
    main()
