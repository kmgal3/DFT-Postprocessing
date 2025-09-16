class InformationCollector:
    """A class to collect all relevant DFT information
    """    
    def __init__(self):
        ## General
        self.description = ''   # e.g. '1,2 Mo z p dxy' means 1st and 2nd Mo atom, in the z direction, all p-orbitals, dxy-orbital 
        self.ions = []          # poscar/ion_types
        self.num_ions = []      # poscar/number_ion_types
        self.folder = ''

        ## Band Structure
        self.kpaths = ''    # e.g. '2 3 -4 -1'
        self.bands = ''     # e.g. '10,12,14-18 22-26'
        self.Pmat = []      # PROCAR matrix shape ( nkpts, nbands, ndirs, nions, norbits )
        self.Emat = []      # EIGENVAL matrix shape ( nkpts, nbands, 2 if SpinPolarized else 1  )
        self.Occp = []      # electron_eigenvalues/fermiweights shape ( nkpts, nbands )
        self.Lvec = []      # poscar/lattice_vectors
        self.Kpts = []      # kpoint_coords
        self.Dvec = []      # distances between kpoints
        self.Klab = []      # kpoints/labels_kpoints
        self.Ktik = []      # positions corresponding to Dvec of the kpoint labels
        self.nSeg = 1       # kpoints/number_line_segments
        self.Scal = 1       # poscar/scale
        self.Klen = 1       # kpoints/number_kpoints
        self.Evalmax = 0
        self.Econmin = 0
        self.Egap = 0
        
        # Density of States
        self.Vdos = []      # electron_dos/dos shape ( nedos, 1 )
        self.Pdos = []      # DOSCAR matrix shape ( nedos, ndirs, nions, norbits )
        self.Edos = []      # electron_dos/energies
        self.Efer = 0       # fermi energy
        self.Emin = 0
        self.Emax = 0
        self.Ediff = 0