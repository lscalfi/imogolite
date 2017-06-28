import numpy as np
import math

unit = 2.2
height = 1.6

from chemfiles import Trajectory, UnitCell, Atom, Topology, Frame, Selection

def coordinate_si(ring, column, state):
    x = unit * column
    y = - unit * ring
    if state == 0 or state == 7:
        return 0, 0, 0, 0, 0, 0, 0
    if state == 1:
        return x, y, x, y, x, y + 1, height
    if state == 4:
        return x, y, x, y, x, y - 1, height
    if state == 2:
        return x, y, x, y, x + 0.866, y + 0.5, height
    if state == 3:
        return x, y, x, y, x + 0.866, y - 0.5, height
    if state == 5:
        return x, y, x, y, x - 0.866, y - 0.5, height
    if state == 6:
        return x, y, x, y, x - 0.866, y + 0.5, height
    return 0, 0, 0, 0, 0, 0, height

def coordinate_up(ring, column, state):
    x = unit * column
    y = - unit * ring
    if state == 7:
        return 0, 0, 0, 0, 0, 0, 0
    if state == 8:
        return x, y, x, y + 1, x, y, height + 1
    if state == 9:
        return x, y, x + 0.866, y - 0.5, x, y, height + 1
    if state == 10:
        return x, y, x - 0.866, y - 0.5, x, y, height + 1
    if state == 11:
        return x, y, x, y + 1, x + 0.866, y - 0.5, height
    if state == 12:
        return x, y, x - 0.866, y - 0.5, x + 0.866, y - 0.5, height
    if state == 13:
        return x, y, x - 0.866, y - 0.5, x, y + 1, height

def coordinate_down(ring, column, state):
    x = unit * column
    y = - unit * ring
    if state == 7:
        return 0, 0, 0, 0, 0, 0, 0
    if state == 8:
        return x, y, x, y - 1, x, y, height + 1
    if state == 9:
        return x, y, x + 0.866, y + 0.5, x, y, height + 1
    if state == 10:
        return x, y, x - 0.866, y + 0.5, x, y, height + 1
    if state == 11:
        return x, y, x, y - 1, x + 0.866, y + 0.5, height
    if state == 12:
        return x, y, x - 0.866, y + 0.5, x + 0.866, y + 0.5, height
    if state == 13:
        return x, y, x - 0.866, y + 0.5, x, y - 1, height

with Trajectory('conf.xyz','w') as output:
    with open('conf.dat', 'r') as conf:
        ring = 0
        skip = 0
        for line in conf:
            if line.startswith("Step:") and skip % 100 == 0:
                frame = Frame()
                ring = 0
            elif line.startswith("####"):
                if skip % 100 == 0:
                    output.write(frame)
                skip += 1
            elif skip % 100 == 0:
                row = map(int,line.split())
                for column, state in enumerate(row):
                    if state < 7:
                        x1, y1, x2, y2, x3, y3, z3 = coordinate_si(ring, column, state)
                    if state > 6 and ring % 4 == 1 and column % 2 == 0: #up (vert)
                        x1, y1, x2, y2, x3, y3, z3 = coordinate_up(ring, column, state)
                    if state > 6 and ring % 4 == 1 and column % 2 == 1: #down (rouge)
                        x1, y1, x2, y2, x3, y3, z3 = coordinate_down(ring, column, state)
                    if state > 6 and ring % 4 == 3 and column % 2 == 1: #up (vert)
                        x1, y1, x2, y2, x3, y3, z3 = coordinate_up(ring, column, state)
                    if state > 6 and ring % 4 == 3 and column % 2 == 0: #down (rouge)
                        x1, y1, x2, y2, x3, y3, z3 = coordinate_down(ring, column, state)
                    if state < 7 and state > 0:
                        atom1 = Atom("Si", "Si")
                        atom2 = Atom("O", "O")
                        atom3 = Atom("H", "H")
                        frame.add_atom(atom1, [x1, y1, 0])
                        frame.add_atom(atom2, [x2, y2, height])
                        frame.add_atom(atom3, [x3, y3, height])
                    if state > 7:
                        atom1 = Atom("Ow", "O")
                        atom2 = Atom("Hw", "H")
                        atom3 = Atom("Hw", "H")
                        frame.add_atom(atom1, [x1, y1, height])
                        frame.add_atom(atom2, [x2, y2, height])
                        frame.add_atom(atom3, [x3, y3, z3])
                ring += 1

