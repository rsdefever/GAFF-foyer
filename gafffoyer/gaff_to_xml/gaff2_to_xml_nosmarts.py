# coding: utf-8

import sys
import re
import math
from lxml import etree as ET

A_TO_NM = 0.1
KCAL_TO_KJ = 4.184

def main():

    if len(sys.argv) != 3:
        print('usage: python gaff2xml.py [path_to_gaff_dat_file] [path_to_output_xml_file]')
        exit(1)
    else:
        gaff_parm_path = sys.argv[1]
        xml_out_path = sys.argv[2]

    start_mass=1; end_mass=99
    start_lj=7522; end_lj=7618
    start_bond=101; end_bond=1036
    start_angle=1037; end_angle=6349
    start_torsion=6350; end_torsion=7478
    start_improper=7479; end_improper=7517

    data = []
    with open(gaff_parm_path) as f:
        for line in f:
            data.append(line.strip())

    # Extract gaff parameters
    mass_parms = {x.split()[0]: x.split()[1] for x in data[start_mass:end_mass]}
    lj_parms = {x.split()[0]: (x.split()[1:3]) for x in data[start_lj:end_lj]}
    bond_parms = [x for x in data[start_bond:end_bond]]
    angle_parms = [x for x in data[start_angle:end_angle]]
    torsion_parms = [x for x in data[start_torsion:end_torsion]]
    improper_parms = [x for x in data[start_improper:end_improper]]

    # Create our force field
    root = ET.Element("ForceField")

    atomtypes = ET.SubElement(root, 'AtomTypes')
    nonbonded = ET.SubElement(root, 'NonbondedForce')
    nonbonded.set('coulomb14scale', '0.833333333')
    nonbonded.set('lj14scale', '0.5')

    # Remove any params that don't have LJ params too
    keys_to_remove = []
    for key in mass_parms.keys():
        if key not in lj_parms.keys():
            keys_to_remove.append(key)
    remove_keys(mass_parms, keys_to_remove)


    elements_to_remove = set()
    for bond in bond_parms:
        classes = re.split('\s+-|-|\s+', bond[0:5])
        for iclass in classes[:2]:
            if iclass not in mass_parms.keys():
                print(iclass)
                elements_to_remove.add(bond)
    remove_elements(bond_parms, elements_to_remove)

    elements_to_remove = set()
    for angle in angle_parms:
        classes = re.split('\s+-|-|\s+', angle[0:8])
        for iclass in classes[:3]:
            if iclass not in mass_parms.keys():
                elements_to_remove.add(angle)
    remove_elements(angle_parms, elements_to_remove)

    elements_to_remove = set()
    for torsion in torsion_parms:
        classes = re.split('\s+-|-|\s+', torsion[0:11])
        for iclass in classes[:4]:
            if iclass not in mass_parms.keys() and iclass.upper() != "X":
                elements_to_remove.add(torsion)
    remove_elements(torsion_parms, elements_to_remove)

    elements_to_remove = set()
    for improper in improper_parms:
        classes = re.split('\s+-|-|\s+', improper[0:11])
        for iclass in classes[:4]:
            if iclass not in mass_parms.keys() and iclass.upper() != "X":
                elements_to_remove.add(improper)
    remove_elements(torsion_parms, elements_to_remove)



    # Atoms/nonbonded (charges set to zero)
    for atype in mass_parms.keys():
        atomtype = ET.SubElement(atomtypes, "Type")
        nb_force = ET.SubElement(nonbonded, "Atom")
        atomtype.set("name", atype)
        atomtype.set("class", atype)
        atomtype.set("element", determine_element(mass_parms[atype]))
        atomtype.set("mass", mass_parms[atype])
        atomtype.set("doi", "10.1021/ja9621760")
        nb_force.set("type", atype)
        nb_force.set("charge", "0.0")
        nb_force.set("sigma", convert_sigma(lj_parms[atype][0]))
        nb_force.set("epsilon", convert_epsilon(lj_parms[atype][1]))

    # Bonds
    bond_forces = ET.SubElement(root,'HarmonicBondForce')
    for bond in bond_parms:
        bond_force = ET.SubElement(bond_forces, 'Bond')
        classes = re.split('\s+-|-|\s+',bond[0:5])
        parms = bond[5:24].split()
        bond_force.set('class1', classes[0])
        bond_force.set('class2', classes[1])
        bond_force.set('length', convert_bondlength(parms[1]))
        bond_force.set('k', convert_bondk(parms[0]))

    # Angles
    angle_forces = ET.SubElement(root,'HarmonicAngleForce')
    for angle in angle_parms:
        angle_force = ET.SubElement(angle_forces, 'Angle')
        classes = re.split('\s+-|-|\s+',angle[0:8])
        parms = angle[8:30].split()
        angle_force.set('class1', classes[0])
        angle_force.set('class2', classes[1])
        angle_force.set('class3', classes[2])
        angle_force.set('angle', convert_theta(parms[1]))
        angle_force.set('k', convert_anglek(parms[0]))

    # Proper dihedrals
    torsion_forces = ET.SubElement(root,'PeriodicTorsionForce')
    continue_reading = False
    for torsion in torsion_parms:
        # Amber leap convention...if periodicity < 0 it means it is
        # there will follow additional torsional terms for the same set of atoms
        # See http://ambermd.org/FileFormats.php#parm.dat
        classes = re.split('\s+-|-|\s+',torsion[0:11])
        parms = torsion[11:54].split()
        if continue_reading == False:
            torsion_force = ET.SubElement(torsion_forces, 'Proper')
            torsion_ctr = 1
            if classes[0].upper() == 'X':
                torsion_force.set('class1', '')
            else:
                torsion_force.set('class1', classes[0])

            if classes[1].upper() == 'X':
                torsion_force.set('class2', '')
            else:
                torsion_force.set('class2', classes[1])

            if classes[2].upper() == 'X':
                torsion_force.set('class3', '')
            else:
                torsion_force.set('class3', classes[2])

            if classes[3].upper() == 'X':
                torsion_force.set('class4', '')
            else:
                torsion_force.set('class4', classes[3])
        else:
            torsion_ctr += 1
        if float(parms[3]) < 0.0:
                continue_reading = True
        else:
            continue_reading = False

        name = 'periodicity' + str(torsion_ctr)
        torsion_force.set(name, str(int(abs(float(parms[3])))))
        name = 'k' + str(torsion_ctr)
        torsion_force.set(name, convert_torsionk(parms[1],parms[0]))
        name = 'phase' + str(torsion_ctr)
        torsion_force.set(name, convert_theta(parms[2]))

    # Improper dihedrals
    continue_reading = False
    for torsion in improper_parms:
        classes = re.split('\s+-|-|\s+',torsion[0:11])
        parms = torsion[11:54].split()
        if continue_reading == False:
            torsion_force = ET.SubElement(torsion_forces, 'Improper')
            torsion_ctr = 1
            if classes[2].upper() == 'X':
                torsion_force.set('class1', '')
            else:
                torsion_force.set('class1', classes[2])

            if classes[0].upper() == 'X':
                torsion_force.set('class2', '')
            else:
                torsion_force.set('class2', classes[0])

            if classes[1].upper() == 'X':
                torsion_force.set('class3', '')
            else:
                torsion_force.set('class3', classes[1])

            if classes[3].upper() == 'X':
                torsion_force.set('class4', '')
            else:
                torsion_force.set('class4', classes[3])
        else:
            torsion_ctr += 1
        if float(parms[2]) < 0.0:
                continue_reading = True
        else:
            continue_reading = False

        name = 'periodicity' + str(torsion_ctr)
        torsion_force.set(name, str(int(abs(float(parms[2])))))
        name = 'k' + str(torsion_ctr)
        torsion_force.set(name, convert_improperk(parms[0]))
        name = 'phase' + str(torsion_ctr)
        torsion_force.set(name, convert_theta(parms[1]))

    # Write XML without SMARTS defs
    ET.ElementTree(root).write(xml_out_path, pretty_print=True)


def determine_element(mass):

    if mass == '1.008':
        element = 'H'
    elif mass == '12.01':
        element = 'C'
    elif mass == '14.01':
        element = 'N'
    elif mass == '16.00':
        element = 'O'
    elif mass == '19.00':
        element = 'F'
    elif mass == '30.97':
        element = 'P'
    elif mass == '32.06':
        element = 'S'
    elif mass == '35.45':
        element = 'Cl'
    elif mass == '63.010':
        element = 'Cu'
    elif mass == '79.90':
        element = 'Br'
    elif mass == '126.9':
        element = 'I'
    else:
        raise ValueError("Invalid element, mass = {}".format(mass))

    return element

def convert_sigma(rmin):
    sigma = 2.0*float(rmin)*A_TO_NM/(2**(1./6.))
    return str(sigma)

def convert_epsilon(amber_epsilon):
    epsilon = float(amber_epsilon)*KCAL_TO_KJ
    return str(epsilon)

def convert_bondlength(length):
    converted_length = float(length)*A_TO_NM
    return str(converted_length)

def convert_bondk(bondk):
    converted_bondk = 2.*float(bondk)*KCAL_TO_KJ/(A_TO_NM*A_TO_NM)
    return str(converted_bondk)

def convert_theta(angle_degrees):
    angle_radians = float(angle_degrees)*math.pi/180.0
    return str(angle_radians)

def convert_anglek(anglek):
    converted_anglek = 2.*float(anglek)*KCAL_TO_KJ
    return str(converted_anglek)

def convert_torsionk(torsionk,idivf):
    converted_torsionk = float(torsionk)*KCAL_TO_KJ/float(idivf)
    return str(converted_torsionk)

def convert_improperk(torsionk):
    converted_torsionk = float(torsionk)*KCAL_TO_KJ
    return str(converted_torsionk)

def remove_keys(dict_, keys):
    for key in keys:
        dict_.pop(key)
        print(f"Removed {key}")

def remove_elements(list_, elements):
    for element in elements:
        list_.remove(element)
        print(f"Removed {element}")


if __name__ == "__main__":
    main()

