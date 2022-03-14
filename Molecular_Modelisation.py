from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def write_file(rp, vp, np, it):
    """Écrit les valeurs de positions et de vitesse de chaque atome de la liste à chaque tour"""
    i = 0
    j = 0
    obfile = open("result.xyz.txt", "a")
    obfile.seek(0, 2)  # on place le curseur à la fin du fichier (pour ne pas écraser les données existantes

    obfile.write("Simulation step : " + str(it) + "\n")
    while i < np:
        obfile.write("Ar    ")
        while j < 3:
            obfile.write("%18.15f" % rp[i][j] + "  ")  # pour avoir un fichier bien ordonné en colonnes
            j += 1
        j = 0

        while j < 3:
            obfile.write(str("%18.15f" % vp[i][j]) + "     ")
            j += 1
        j = 0
        i += 1
        obfile.write(str(i) + "\n")

    obfile.write("\n\n")
    obfile.close()

def create_file(name):
    open(name, "w").close()  # crée le fichier s'il n'existe pas et efface tout son contenu s'il en a (clear)


def read_file(rp, vp, np):
    """Lit le fichier au format xyz et initialise
    dles variables de position et de vitesse"""
    i = 0

    obfile = open('lj_init.xyz.txt', 'r')
    obfile.readline()
    obfile.readline()

    while i < np:
        strin = obfile.readline()
        liste = strin.split()
        rp[i][0] = float(liste[1])
        rp[i][1] = float(liste[2])
        rp[i][2] = float(liste[3])
        vp[i][0] = float(liste[4])
        vp[i][1] = float(liste[5])
        vp[i][2] = float(liste[6])
        i += 1
    obfile.close()

    return rp, vp


def init(rp, vp, fp, tmp_ap, ap, np):
    """Initialise toutes les listes, pour
    en faire des listes à deux dimensions de la taille coordonnées * np"""
    ip = 0

    while ip < np:
        rp.append([0.0, 0.0, 0.0])
        vp.append([0.0, 0.0, 0.0])
        fp.append([0.0, 0.0, 0.0])
        tmp_ap.append([0.0, 0.0, 0.0])
        ap.append([0.0, 0.0, 0.0])
        ip += 1

    return rp, vp, fp, tmp_ap, ap

def movement(epsilon, sigma, np):
    """Calcule les positions et les vitesses d'une liste d'atomes dans le temps en fonction des paramètres"""

    ip = 0  # variable à incrémenter
    it = 0  # variable à incrémenter (pour le fichier)
    rp = []  # la liste des positions
    vp = []  # la liste des vitesses
    fp = []  # la liste des forces
    tmp_ap = []  # la liste temporaire des accélérations
    ap = []  # la liste des accélérations


    rp, vp, fp, tmp_ap, ap = init(rp, vp, fp, tmp_ap, ap, np)
    rp, vp = read_file(rp, vp, np)
    create_file("result.xyz.txt")

    r12 = [0, 0, 0]  # distance entre deux points sur les 3 axes x, y et z
    lbox = 25  # taille de la box en Angstrom
    ulj_tot = 0.0   # énergie totale
    m = 40  # masse de l'argon en ua
    u_cin = 0

    dt = 0.02  # pas de temps
    t_max = 3  # temps maximum
    t = 0  # variable du temps à incrémenter

    p = 0
    k = 0

    # initialisation des valeurs de fp et ap
    while ip < np:
        jp = 0
        while jp < ip:
            i_xyz = 0
            while i_xyz < 3:
                r12[i_xyz] = rp[ip][i_xyz] - rp[jp][i_xyz]
                if r12[i_xyz] > 0.5 * lbox:
                    r12[i_xyz] = r12[i_xyz] - lbox
                elif r12[i_xyz] < -0.5 * lbox:
                    r12[i_xyz] = r12[i_xyz] + lbox
                i_xyz += 1
            r12c = 0.0
            i_xyz = 0

            while i_xyz < 3:
                r12c += r12[i_xyz] ** 2
                i_xyz += 1
            r12a = sqrt(r12c)
            ulj = 4.0 * epsilon * ((sigma / r12a) ** 12 - (sigma / r12a) ** 6)
            ulj_tot += ulj
            flj = (24.0 * epsilon * (2.0 * (sigma / r12a) ** 12 - (sigma / r12a) ** 6) / r12a) / 300
            i_xyz = 0

            while i_xyz < 3:
                f12xyz = flj * r12[i_xyz] / r12a
                fp[ip][i_xyz] += f12xyz
                fp[jp][i_xyz] -= f12xyz
                i_xyz += 1

            jp += 1
        ip += 1

    ip = 0
    while ip < np:
        jp = 0
        while jp < np:
            i_xyz = 0
            while i_xyz < 3:
                ap[ip][i_xyz] = fp[ip][i_xyz] / m
                ap[jp][i_xyz] = fp[jp][i_xyz] / m
                i_xyz += 1
            jp += 1
        ip += 1

    # calculs
    while t < t_max:
        write_file(rp, vp, np, it)
        # calcul des nouvelles positions
        ip = 0
        while ip < np:
            i_xyz = 0
            while i_xyz < 3:
                rp[ip][i_xyz] += dt * vp[ip][i_xyz] + dt ** 2 * ap[ip][i_xyz] / 2
                i_xyz += 1
            ip += 1

        # calcul de la distance entre les atomes
        ip = 0
        while ip < np:
            jp = 0
            while jp < ip:
                i_xyz = 0
                while i_xyz < 3:
                    r12[i_xyz] = rp[ip][i_xyz] - rp[jp][i_xyz]
                    if r12[i_xyz] > 0.5 * lbox:
                        r12[i_xyz] = r12[i_xyz] - lbox
                    elif r12[i_xyz] < -0.5 * lbox:
                        r12[i_xyz] = r12[i_xyz] + lbox
                    # periodicité
                    if r12[i_xyz] > lbox / 2:
                        r12[i_xyz] -= lbox
                    if r12[i_xyz] > lbox / 2:
                        r12[i_xyz] -= lbox
                    if r12[i_xyz] > lbox / 2:
                        r12[i_xyz] -= lbox

                    if r12[i_xyz] < - lbox / 2:
                        r12[i_xyz] += lbox
                    if r12[i_xyz] < - lbox / 2:
                        r12[i_xyz] += lbox
                    if r12[i_xyz] < - lbox / 2:
                        r12[i_xyz] += lbox
                    i_xyz += 1

                # calcul des énergies
                r12c = 0.0
                i_xyz = 0
                while i_xyz < 3:
                    r12c += r12[i_xyz] ** 2
                    i_xyz += 1
                r12a = sqrt(r12c)
                ulj = 4.0 * epsilon * ((sigma / r12a) ** 12 - (sigma / r12a) ** 6)
                ulj_tot += ulj
                flj = (24.0 * epsilon * (2.0 * (sigma / r12a) ** 12 - (sigma / r12a) ** 6) / r12a) / 300.0

                # calcul des forces
                i_xyz = 0
                while i_xyz < 3:
                    f12xyz = flj * r12[i_xyz] / r12a
                    fp[ip][i_xyz] += f12xyz
                    fp[jp][i_xyz] -= f12xyz
                    i_xyz += 1

                # mise en tmp des accélérations
                i_xyz = 0
                while i_xyz < 3:
                    tmp_ap[ip][i_xyz] = ap[ip][i_xyz]
                    i_xyz += 1
                jp += 1
            ip += 1

        # calcul des accélérations
        ip = 0
        while ip < np:
            jp = 0
            while jp < ip:
                i_xyz = 0
                while i_xyz < 3:
                    ap[ip][i_xyz] = fp[ip][i_xyz] / m
                    ap[jp][i_xyz] = fp[jp][i_xyz] / m
                    i_xyz += 1
                jp += 1
            ip += 1

        # calcul des nouvelles vitesses
        ip = 0
        while ip < np:
            i_xyz = 0
            while i_xyz < 3:
                vp[ip][i_xyz] += dt * (ap[ip][i_xyz] + tmp_ap[ip][i_xyz]) / 2
                u_cin += m * vp[ip][i_xyz] ** 2 / 2
                i_xyz += 1
            ip += 1

        # affichage
        if p == 5:
            plt.close()
            fig = plt.figure()
            ax = Axes3D(fig)  # définition du repère 3D
            while k < np:
                ax.scatter(rp[k][0], rp[k][1], rp[k][2], c=[[k / 100, k / 100, k / 100]], marker="o")  # ajout des points avant affichage et définition des couleurs
                k += 1
            k = 0
            print("Affichage au temps t = : ", t)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_zlim(0, lbox)
            plt.xlim(0, lbox)
            plt.ylim(0, lbox)
            p = 0

        plt.draw()  # affichage
        plt.pause(1e-17)
        p += 1

        t += dt
        it += 1

    print(ulj_tot, u_cin)


epsilon = 112
sigma = 3.6
np = 100  # nombre de particules

movement(epsilon, sigma, np)
