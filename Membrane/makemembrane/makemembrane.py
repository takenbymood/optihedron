#hex lattice of flat membrane beads
membraneXYRatio = 15.0/26.0 #rhombus count Y/X
boxLength = 50
rhombusDimX = 1.0
outputForLAMMPS = 'rigidmembrane.data'

#lattice parameters
rhombusDimY = rhombusDimX*(3.0**0.5)

rhombusCountDimX = int(round(2*boxLength/rhombusDimX))
rhombusCountDimY = int(round(rhombusCountDimX*membraneXYRatio))

membraneLengthDimX = (rhombusCountDimX*rhombusDimX)/2.0
membraneLengthDimY = (rhombusCountDimY*rhombusDimY)/2.0
membraneLengthDimZ = 100*rhombusDimX

membraneBeadPositions = [] #[[x,y,z]]
for aRhombusAlongDimX in range(rhombusCountDimX):
    for aRhombusAlongDimY in range(rhombusCountDimY):
        membraneBeadPositions.append([(rhombusDimX*aRhombusAlongDimX)-membraneLengthDimX+(rhombusDimX/4.0), (rhombusDimY*aRhombusAlongDimY)-membraneLengthDimY+(rhombusDimY/4.0), 0.0])
        membraneBeadPositions.append([(rhombusDimX*aRhombusAlongDimX)-membraneLengthDimX+(3.0*rhombusDimX/4.0), (rhombusDimY*aRhombusAlongDimY)-membraneLengthDimY+(3.0*rhombusDimY/4.0), 0.0])

membraneBeadCount = len(membraneBeadPositions)

with open( outputForLAMMPS, 'w' ) as file:
    file.write('Hex lattice of flat membrane beads.\n\n')

    file.write('{0} atoms\n'.format(membraneBeadCount))
    file.write('1 atom types\n\n')    

    file.write('-{0} {0} xlo xhi\n'.format(membraneLengthDimX))
    file.write('-{0} {0} ylo yhi\n'.format(membraneLengthDimY))
    file.write('-{0} {0} zlo zhi\n\n'.format(membraneLengthDimZ))

    file.write('Masses\n\n')
    file.write('1 1\n\n')

    file.write('Atoms\n\n')

    for membraneBeadIndex, membraneBead in enumerate(membraneBeadPositions, 1):
        file.write('{0} 1 {1} {2} {3}   1 1 0   0 0 1   0 0 0\n'.format(membraneBeadIndex, membraneBead[0], membraneBead[1], membraneBead[2]))
    file.write('\n')

    file.write('Velocities\n\n')    
    for membraneBeadIndex in range(len(membraneBeadPositions)):
        file.write('{0} 0 0 0 0 0 0\n'.format(membraneBeadIndex+1))