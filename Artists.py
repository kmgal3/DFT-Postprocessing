import numpy as np

def BandWriter(s: str) -> list[int]:
    """Parse a string of ranges and values (separated by commas or whitespace)
    into a 0-based Python list.

    Args:
        s (str): eg. "4,5 7-9"

    Returns:
        list[int]: eg. [3,4,6,7,8]
    """
    bands = []
    
    # Replace commas with spaces, then split on spaces
    tokens = s.replace(",", " ").split()
    
    for token in tokens:
        if "-" in token:  # range
            start, end = map(int, token.split("-"))
            bands.extend(range(start - 1, end))  # 0-based
        else:  # single number
            num = int(token)
            bands.append(num - 1)
    
    return bands if bands != [] else ...


def ProjectionPainter(ions:list, num_ions:list, description:str) -> dict:
    """Translate an intuitive string to select proper elements of the projection matrix
    Args:
        ions (list): Names of ions
        num_ions (list): Number of each of the ions POSCAR style
        description (str): The string to translate, eg. "Cr down py pz"

    Returns:
        dict: {'direction':list, 'ion':list, 'orbital':list}
    """    
    orbital_names = {'s':[0], 'p':[1,2,3], 'd':[4,5,6,7,8], 'f':[i for i in range(9,16)],
                     'sp':[0,1,2,3,4], 'pd':[i for i in range(1,9)], 'spd':[i for i in range(0,9)],
                     'py':[1], 'pz':[2], 'px':[3],
                     'dxy':[4], 'dyz':[5], 'dz2':[6], 'dxz':[7], 'x2-y2':[8],
                     'fy3x2':[9], 'fxyz':[10], 'fyz2':[11], 'fz3':[12], 'fxz2':[13], 'fzx2':[14], 'fx3':[15]}
        
    direction_names = {'up':[0], 'down':[1], 
                    'tot':[0], 'x':[1], 'y':[2], 'z':[3]}

    # For each ion type note which indexes it occupies
    now = 0
    list_ions = []
    atom_indexes = {}
    for n in num_ions:
        list_ions.append([i for i in range(now, now+n)])
        now += n
    for name,pos in zip(ions,list_ions): atom_indexes[name] = pos

    projection = {'direction':[], 'ion':[], 'orbital':[]}

    description = description.replace(",", " ").split()

    iatoms = []
    for i in description:
        if i.isnumeric(): iatoms.append(i)

    # Gather the ion projection list
    for atom in atom_indexes.keys():
        if atom in description:
            if len(iatoms):
                for iatom in iatoms:
                    projection['ion'].append(atom_indexes[atom][int(iatom)-1])
            else: projection['ion'] = atom_indexes[atom]
    if projection['ion']==[]: projection['ion'] = ...

    # Gather the direction projection list
    for dir in direction_names.keys():
        if dir in description: projection['direction'] += direction_names[dir]
    if projection['direction']==[]: projection['direction'] = [0]

    # Gather the orbital projection list
    for orb in orbital_names.keys():
        if orb in description: projection['orbital'] += orbital_names[orb]
    if projection['orbital']==[]: projection['orbital'] = ...

    return projection


def SegmentArranger(nSeg:int, mat: np.ndarray, arrangement: str = ''):
    """Reordering a given matrix with dimention corresponding to nkpts,
    to suit any chosen combination of segments.

    Args:
        nSeg (int): the length of each segment (must be the same for each)
        mat (np.ndarray): any array with shape (nkpts, ...) 
        arrangement (str, optional): a sequence of segments. Might include negative numbers for reverse order. Defaults to ''.

    Returns:
        np.ndarray: Rearranged matrix
    """    
    if arrangement == '': return mat

    mat = np.split(mat, nSeg)
    arrangement = arrangement.replace(",", " ").split()

    newmat = []
    for i in arrangement:
        i = int(i)
        if i > 0: 
            # If the list is still empty set the current segment
            try: newmat = np.concatenate([newmat,mat[i-1]])
            except ValueError: newmat = mat[i-1]
        if i < 0: 
            try: newmat = np.concatenate([newmat,mat[np.abs(i)-1][::-1,...]])
            except ValueError: newmat = mat[np.abs(i)-1][::-1,...]

    return newmat


def MatrixSculptor(mat: np.ndarray, 
                   projection = {'direction':[0], 'ion':..., 'orbital':...}
                   ) -> np.ndarray:
    """Pruning last 3 dimentions by summation of chosen indexes.

    Args:
        mat (np.ndarray): any matrix with shape (..., ndirs, nions, norbits)
        projection (dict, optional): a ProjectionPainter output. Defaults to {'direction':[0], 'ion':..., 'orbital':...}.

    Returns:
        np.ndarray: an array with 3 less dimentions
    """
    try: mat = np.sum(mat[..., projection['orbital']], axis=-1)
    except IndexError: mat = np.sum(mat, axis=-1)
    try: mat = np.sum(mat[..., projection['ion']], axis=-1)
    except IndexError: mat = np.sum(mat, axis=-1)
    try: mat = np.sum(mat[..., projection['direction']], axis=-1)
    except IndexError: mat = np.sum(mat, axis=-1)
    return mat


def KlineEngineer(Lvec, Scal, Kpts, Klen, Klab, 
                  nSeg, 
                  arrangement: str = ''
                  ):
    """Takes care of preparing everything, so that the x-axis 
        of the band structure plot will look as it should

    Args:
        Lvec (_type_): _description_
        Scal (_type_): _description_
        Kpts (_type_): _description_
        Klen (_type_): _description_
        Klab (_type_): _description_
        nSeg (_type_): _description_
        arrangement (str, optional): _description_. Defaults to ''.

    Returns:
        _type_: _description_
    """
    Kpts = SegmentArranger(nSeg, Kpts, arrangement)
    Dvec = Kpts.copy()

    # calculate reciprocal lattice vectors and constants
    a1,a2,a3 = Lvec*Scal
    bM = np.zeros((3,3))
    V = np.dot(a1, np.cross(a2, a3))
    bM[0,...] = 2*np.pi*np.cross(a2, a3) / V
    bM[1,...] = 2*np.pi*np.cross(a3, a1) / V
    bM[2,...] = 2*np.pi*np.cross(a1, a2) / V
    Rval = np.sum(bM**2, axis = 1)**.5

    # calulate the difference beteen subsequent kpoints
    Dvec *= Rval
    Dvec[1:] -= Dvec[:-1]
    Dvec = np.cumsum(np.sum((Dvec)**2, axis = 1)**.5)

    # place k-point labels at the start and end of each segment
    Nkpt = len(Dvec)
    Ktik = []
    for i in range(0,Nkpt//Klen):
        Ktik.append(i*Klen)
        Ktik.append((i+1)*Klen-1)
    Ktik = Dvec[Ktik]
    Klab = SegmentArranger(nSeg, Klab, arrangement)
    Klab = [lab if lab!='GAMMA' else 'Γ' for lab in Klab]

    return Kpts, Dvec, Ktik, Klab