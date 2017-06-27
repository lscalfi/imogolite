import numpy as np
import random
import math

#units: kcal/mol

diameter = 12
units = 5
p_ads = 0.4
p_des = 0.4
p_silanol = 0.2

def kron(a, b):
    if a == True and b == True:
        return -1
    return 0

def infi(a, b):
    if a == True and b == True:
        return 2.0e31
    return 0

def energy(sub_nt):
    H = 0
    J = 6
    assert(sub_nt.shape == (3, 3))
    if sub_nt[0, 0] == 0: #green type
        H += infi(sub_nt[1, 1] == 8, sub_nt[0, 1] == 4)
        H += J * kron(sub_nt[1, 1] == 8, sub_nt[2, 0] == 2)
        H += J * kron(sub_nt[1, 1] == 8, sub_nt[2, 2] == 6)

        H += J * kron(sub_nt[1, 1] == 9, sub_nt[0, 1] == 4)
        H += J * kron(sub_nt[1, 1] == 9, sub_nt[2, 0] == 2)
        H += infi(sub_nt[1, 1] == 9, sub_nt[2, 2] == 6)
 
        H += J * kron(sub_nt[1, 1] == 10, sub_nt[0, 1] == 4)
        H += infi(sub_nt[1, 1] == 10, sub_nt[2, 0] == 2)
        H += J * kron(sub_nt[1, 1] == 10, sub_nt[2, 2] == 6)
 
        H += infi(sub_nt[1, 1] == 11, sub_nt[0, 1] == 4)
        H += J * kron(sub_nt[1, 1] == 11, sub_nt[2, 0] == 2)
        H += infi(sub_nt[1, 1] == 11, sub_nt[2, 2] == 6)
 
        H += J * kron(sub_nt[1, 1] == 12, sub_nt[0, 1] == 4)
        H += infi(sub_nt[1, 1] == 12, sub_nt[2, 0] == 2)
        H += infi(sub_nt[1, 1] == 12, sub_nt[2, 2] == 6)
 
        H += infi(sub_nt[1, 1] == 13, sub_nt[0, 1] == 4)
        H += infi(sub_nt[1, 1] == 13, sub_nt[2, 0] == 2)
        H += J * kron(sub_nt[1, 1] == 13, sub_nt[2, 2] == 6)
 
    else : #red type
        H += infi(sub_nt[1, 1] == 8, sub_nt[2, 1] == 1)
        H += J * kron(sub_nt[1, 1] == 8, sub_nt[0, 0] == 3)
        H += J * kron(sub_nt[1, 1] == 8, sub_nt[0, 2] == 5)

        H += J * kron(sub_nt[1, 1] == 9, sub_nt[2, 1] == 1)
        H += J * kron(sub_nt[1, 1] == 9, sub_nt[0, 0] == 3)
        H += infi(sub_nt[1, 1] == 9, sub_nt[0, 2] == 5)
 
        H += J * kron(sub_nt[1, 1] == 10, sub_nt[2, 1] == 1)
        H += infi(sub_nt[1, 1] == 10, sub_nt[0, 0] == 3)
        H += J * kron(sub_nt[1, 1] == 10, sub_nt[0, 2] == 5)
 
        H += infi(sub_nt[1, 1] == 11, sub_nt[2, 1] == 1)
        H += J * kron(sub_nt[1, 1] == 11, sub_nt[0, 0] == 3)
        H += infi(sub_nt[1, 1] == 11, sub_nt[0, 2] == 5)
 
        H += J * kron(sub_nt[1, 1] == 12, sub_nt[2, 1] == 1)
        H += infi(sub_nt[1, 1] == 12, sub_nt[0, 0] == 3)
        H += infi(sub_nt[1, 1] == 12, sub_nt[0, 2] == 5)
 
        H += infi(sub_nt[1, 1] == 13, sub_nt[2, 1] == 1)
        H += infi(sub_nt[1, 1] == 13, sub_nt[0, 0] == 3)
        H += J * kron(sub_nt[1, 1] == 13, sub_nt[0, 2] == 5)

    return H

def print_nt(nt):
    for ring in range(0, 2 * 2 * units):
        for silanol in range(0, 2*diameter):
            out.write(str(int(nt[ring, silanol])))
            out.write("  ")
        out.write("\n")
    out.write("###########################################\n")

def initialize(nanotube):
    for ring in range(0, 2 * 2 * units, 4):
        for silanol in range(0, 2*diameter, 2):
            nanotube[ring, silanol] = random.randint(1,6)
    for ring in range(2, 2 * 2 * units, 4):
        for silanol in range(1, 2*diameter, 2):
            nanotube[ring, silanol] = random.randint(1,6)
    for ring in range(1, 2 * 2 * units, 2):
        for silanol in range(0, 2*diameter):
            nanotube[ring, silanol] = 7

def sub_nanotube(nanotube, i, j):
    if i < 1 or j < 1 or i > 2 * 2 * units - 3 or j > 2 * diameter - 3 : 
        pbc = np.zeros((2 * 2 * units + 4, 2 * diameter + 4))
        pbc[2 : -2, 2 : -2] = nanotube[:, :]
        pbc[0, 2 : -2] = nanotube[-2, :]
        pbc[1, 2 : -2] = nanotube[-1, :]
        pbc[-1, 2 : -2] = nanotube[1, :]
        pbc[-2, 2 : -2] = nanotube[0, :]
        pbc[2 : -2, 0] = nanotube[:, -2]
        pbc[2 : -2, 1] = nanotube[:, -1]
        pbc[2 : -2, -1] = nanotube[:, 1]
        pbc[2 : -2, -2] = nanotube[:, 0]
        pbc[0 : 2, 0 : 2] = nanotube[-2:, -2:]
        pbc[0 : 2, -2:] = nanotube[-2:, 0 : 2]
        pbc[-2:, 0 : 2] = nanotube[0 : 2, -2:]
        pbc[-2:, -2:] = nanotube[0 : 2, 0 : 2]
        i += 2
        j += 2
           
        return pbc[i - 1 : i + 2, j - 1 : j + 2]
    return nanotube[i - 1 : i + 2, j - 1 : j + 2]

def adsorption(tube):
    nanotube = np.copy(tube)
    #try adsorption or rotation of water
    print("Adsorption move")
    trial_ring = random.randrange(1, 2 * 2 * units, 2)
    trial_triangle = random.randint(0, 2 * diameter - 1)

    print(str(trial_ring) + ', ' + str(trial_triangle))

    ads = False
    if nanotube[trial_ring, trial_triangle] == 7:
        ads = True

    orientation = random.randint(8, 13)

    sub_nt = sub_nanotube(nanotube, trial_ring, trial_triangle)
    H1 = energy(sub_nt)
    nanotube[trial_ring, trial_triangle] = orientation
    sub_nt[1, 1] = orientation
    H2 = energy(sub_nt)
    dH = H2 - H1

    print("dH = " + str(dH))
    return nanotube, dH, ads

def desorption(tube):
    nanotube = np.copy(tube)
    print("Desorption move")
    trial_ring = random.randrange(1, 2 * 2 * units, 2)
    trial_triangle = random.randint(0, 2 * diameter - 1)
    
    print(str(trial_ring) + ', ' + str(trial_triangle))

    des = False
    dH = 2e+32

    if not (nanotube[trial_ring, trial_triangle] == 7) :
        sub_nt = sub_nanotube(nanotube, trial_ring, trial_triangle)
        dH = - energy(sub_nt)
        des = True

    nanotube[trial_ring, trial_triangle] = 7

    print("dH = " + str(dH))
    return nanotube, dH, des

def rotationSi(tube):
    nanotube = np.copy(tube)
    print("Rotation of silanol move")
    trial_ring = random.randrange(0, 2 * 2 * units, 2)
    trial_triangle = 2 * random.randint(0, diameter - 1)

    if trial_ring % 4 == 2:
        trial_triangle += 1

    print(str(trial_ring) + ', ' + str(trial_triangle))

    orientation = random.randint(1, 6)

    sub_nt1 = sub_nanotube(nanotube, trial_ring - 1, trial_triangle - 1)
    sub_nt2 = sub_nanotube(nanotube, trial_ring - 1, trial_triangle)
    sub_nt3 = sub_nanotube(nanotube, trial_ring - 1, trial_triangle + 1)
    sub_nt4 = sub_nanotube(nanotube, trial_ring + 1, trial_triangle - 1)
    sub_nt5 = sub_nanotube(nanotube, trial_ring + 1, trial_triangle)
    sub_nt6 = sub_nanotube(nanotube, trial_ring + 1, trial_triangle + 1)
    H1 = energy(sub_nt1) + energy(sub_nt2) + energy(sub_nt3) + energy(sub_nt4) + energy(sub_nt5) + energy(sub_nt6)
    nanotube[trial_ring, trial_triangle] = orientation
    sub_nt1 = sub_nanotube(nanotube, trial_ring - 1, trial_triangle - 1)
    sub_nt2 = sub_nanotube(nanotube, trial_ring - 1, trial_triangle)
    sub_nt3 = sub_nanotube(nanotube, trial_ring - 1, trial_triangle + 1)
    sub_nt4 = sub_nanotube(nanotube, trial_ring + 1, trial_triangle - 1)
    sub_nt5 = sub_nanotube(nanotube, trial_ring + 1, trial_triangle)
    sub_nt6 = sub_nanotube(nanotube, trial_ring + 1, trial_triangle + 1)
    H2 = energy(sub_nt1) + energy(sub_nt2) + energy(sub_nt3) + energy(sub_nt4) + energy(sub_nt5) + energy(sub_nt6)
    dH = H2 - H1

    print("dH = " + str(dH))
    return nanotube, dH

def accept(dH, beta):
    if dH < 0:
        return True
    elif dH > 2e30:
        return False
    else:
        ran = random.random()
        if ran < math.exp(-beta * dH):
            return True
    return False

def main(steps, temp):
    random.seed()
    H = 0
    beta = 1/(1.987e-3*temp)
    nt = np.zeros((2 * 2 * units, 2 * diameter))

    initialize(nt)
    n_ads = 0

    for n in range(steps):
        out.write("Step: " + str(n) + '\n')
        print("Step: " + str(n))
        trial = random.uniform(0, 1)

        if trial < p_ads :
            nt_post, dH, ads = adsorption(nt)
            if accept(dH, beta) :
                print("Accepted")
                nt = nt_post
                H += dH
                if ads :
                    n_ads +=1
            else :
                print("Rejected")
            print_nt(nt)
            print("N_ads = " + str(n_ads) + "\n")

        elif trial < p_ads + p_des :
            nt_post, dH, des = desorption(nt)
            if accept(dH, beta) :
                print("Accepted")
                nt = nt_post
                H += dH
                if des :
                    n_ads -=1
            else :
                print("Rejected")
            print_nt(nt)
            print("N_ads = " + str(n_ads) + "\n")

        else :
            nt_post, dH = rotationSi(nt)
            if accept(dH, beta) :
                print("Accepted")
                nt = nt_post
                H += dH
            else :
                print("Rejected")
            print_nt(nt)
            print("N_ads = " + str(n_ads) + "\n")
        print("Energy = " + str(H) + " kcal/mol\n")


out = open("conf.dat","w")
main(1000000, 300)
out.close()

