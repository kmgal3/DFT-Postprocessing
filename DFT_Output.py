import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from tqdm import tqdm
from matplotlib import colormaps
from matplotlib.lines import Line2D

class PhysicalConstants:
    pi      = 3.141592653589793     # no unit
    c       = 299792458;            # m/s
    e       = 1.602176565E-19;      # As = C
    me      = 9.10938291E-31;       # kg
    k       = 1.3806488E-23;        # J/K
    h       = 6.62606957E-34;       # Js = Nms = VAs^2
    hbar    = 6.62606957E-34/2/pi;  # Js = Nms = VAs^2
    eps0    = 8.854187817E-12;      # dielectric constant [As/Vm = F/m]
    a0      = 0.52917721092E-10;    # Bohr radius [m]
    Ry      = 2.179872171E-18;      # Rydberg energy [J] 
    u       = 1.660538921E-27;      # atomic mass unit [kg]
    NA      = 6.02214129E23;        # Avogadro constant[1/mol]

class DFT_Output :
    """This class contains various methods for calculation and analysis of 
    parameters important in g-factor analysis, such as orbital angular momentum 
    and optical selection rules
    """    
    
    def __init__ (self, source, dir = '', name = ''):
        """_summary_

        Args:
            source (str): The program, which supplies the relevant DFT info
            dir (str, optional): The directory in which the files are contained. Defaults to ''.
            name (str, optional): An optional name for the mommat2up file. Defaults to ''.
        """        
        self.WAVEDER    = None
        self.EIGENVAL   = None
        self.WAVEDERF   = None
        self.PROCAR     = None
        self.DOSCAR     = None
        self.mommat2up  = None

        if source == 'VASP':
            self.WAVEDER    = dir + "WAVEDER"
            self.EIGENVAL   = dir + "EIGENVAL"
            self.WAVEDERF   = dir + "WAVEDERF"
            self.PROCAR     = dir + "PROCAR"
            self.DOSCAR     = dir + "DOSCAR"

        elif source == 'Wien2K' :
            self.mommat2up  = dir + name

        elif source == 'DFTB+' :
            self.bandout   = dir + "band.out"

        self.source = source
        self.name = name
        self.band_vecs = {}
        self.P_matrix  = np.array([])
        self.eval = np.array([])
        self.occ = np.array([])
        self.dir = dir

    def set_VASP_files(self, WAVEDER = None, EIGENVAL = None, WAVEDERF = None):
        if WAVEDER != None: self.WAVEDER  = WAVEDER
        if EIGENVAL != None: self.EIGENVAL = EIGENVAL
        if WAVEDERF != None: self.WAVEDERF = WAVEDERF

    def set_wien_file(self, file = None):
        self.mommat2up = file

    def read_P_matrix(self, file):
        self.P_matrix = np.load(file, mmap_mode='r')

    def __getitem__ (self, index):
        if self.P_matrix.size > 0 :
            return self.P_matrix[:,index-1,:,:]
        elif index in self.known_bands():
            return self.band_vecs[index]
        else : 
            self.read_bands([index])
            return self.band_vecs[index]
    
    def known_bands (self):
        bands = list(self.band_vecs.keys())
        bands.sort()
        return bands

    def VASP_ThinRead(self, bands):
        vecs = self.band_vecs
        for band in bands:
            vecs[band] = []
        with open(self.WAVEDERF,'r') as fh:
            Q = np.fromstring(fh.readline(), dtype=int, sep = ' ')
            nkpts = Q[1]
            nbands = Q[2]
            ibands = bands.copy()
            for line in fh:
                if len(ibands) == 0: break
                for band in ibands:
                    if f" {band} " in line[:10]:
                        vecs[band].append(np.fromstring(line, dtype=float, sep=' '))
                        if len(vecs[band]) == nkpts*nbands : ibands.remove(band)

        for band in bands:
            vec = np.array(vecs[band])
            vec = np.concatenate([vec[:,6:7]   + vec[:,7:8]*1j,
                                  vec[:,8:9]   + vec[:,9:10]*1j,
                                  vec[:,10:11] + vec[:,11:12]*1j,
                                  -vec[:,4:5]  + vec[:,1:2]     ], axis = 1)
            vec[:,:3] = vec[:,:3] * np.abs(vec[:,3:4]) 
            vecs[band] = vec.reshape((nkpts,nbands,4))

        return vecs
    
    def wien_ThinRead(self, bands):
        nkpnts, nbands = 0, 0
        nlines, skip   = 0, 0     # Optimisation variables
        vecs = self.band_vecs
        for band in bands:
            vecs[band] = []
        with open(self.mommat2up,'r') as fh:
            ibands = bands.copy()
            for line in fh:
                if skip: skip -= 1 ; continue
                if len(ibands) == 0: 
                    skip = nlines
                    ibands = bands.copy()
                if "KP:" in line:
                    nkpnts += 1
                    nbands = int(line.split()[6])
                    nlines = int(0.5*nbands*(nbands+1))
                for band in ibands:
                    if f" {band} " in line:
                        vecs[band].append(np.fromstring(line, dtype=float, sep=' '))
                        if len(vecs[band]) == nkpnts * nbands and band != 1: ibands.remove(band)
                nlines -= 1

        for band in bands:
            vec = np.array(vecs[band])
            vec = np.concatenate([vec[:,2:3]+vec[:,3:4]*1j,
                                vec[:,4:5]+vec[:,5:6]*1j,
                                vec[:,6:7]+vec[:,7:8]*1j,
                                vec[:,8:9]], axis = 1)
            vec[:band-1,0:3] = np.conj(vec[:band-1,0:3])
            vec[:band-1,3] = -vec[:band-1,3]
            vecs[band] = vec.reshape((nkpnts,nbands,4))

        return vecs
    
    def read_bands(self, bands):
        if self.mommat2up != None : return self.wien_ThinRead(bands)
        else : return self.VASP_ThinRead(bands)

    def Optical_Selection_Rules(self, kpoint, vb_max, num_vb, num_cb, save = '', txt = False):
        # calculate the optical matrix elements for circular and linear polarized 
        # light and determine the dipole strengths and degree of circular polarization
        #
        # the k-point is defined by the WAVEDERF file
        # the transitions are given in the vb_list and cb_list

        # find VBM

        vb_list = [vb_max - num_vb + i + 1 for i in range(num_vb)]
        cb_list = [vb_max + i + 1 for i in range(num_cb)]

        bands = vb_list + cb_list

        if self.P_matrix.size == 0 :
            bands_to_read = []
            for band in bands:
                if band not in self.known_bands():
                    bands_to_read.append(band)
            self.read_bands(bands_to_read)

        optical = []

        # polarization verctors
        e_plus  = 1/np.sqrt(2)*np.array([[1,  1j, 0]])
        e_minus = 1/np.sqrt(2)*np.array([[1, -1j, 0]])

        vecs = self.band_vecs
        e_flag = 0

        for vb in vb_list:
            for cb in cb_list:
                # momentum matrix element
                # complex valued and not uniquely defined by an arbitray phase

                if self.P_matrix.size > 0 :
                    p = np.array([[self.P_matrix[kpoint-1, cb-1, vb-1, 0], 
                                   self.P_matrix[kpoint-1, cb-1, vb-1, 1], 
                                   self.P_matrix[kpoint-1, cb-1, vb-1, 2]]])

                else :
                    p = np.array([[vecs[cb][kpoint-1, vb-1, 0], 
                                   vecs[cb][kpoint-1, vb-1, 1], 
                                   vecs[cb][kpoint-1, vb-1, 2]]])

                # circular polarization
                p_plus  = e_plus @ p.T  
                p_minus = e_minus @ p.T

                # oscillator strengths
                p_plus_squ  =   p_plus[0,0] * np.conj(p_plus[0,0])
                p_minus_squ =  p_minus[0,0] * np.conj(p_minus[0,0])
                p_x_squ     =        p[0,0] * np.conj(p[0,0])
                p_y_squ     =        p[0,1] * np.conj(p[0,1])
                p_z_squ     =        p[0,2] * np.conj(p[0,2])
                p_squ       = p_x_squ + p_y_squ + p_z_squ
                
                # degree of circular polarization
                polariz = (p_plus_squ - p_minus_squ)/(p_plus_squ + p_minus_squ)
                
                optical += [vb, cb, polariz, p_plus_squ, p_minus_squ, 
                            p_x_squ, p_y_squ, p_z_squ, p_squ]
                
                if self.eval.size > 0 : 
                    optical += [self.eval[cb-1,0] - self.eval[vb-1,0]]
                    e_flag = 1
                
        optical = np.real(np.array(optical)).reshape(-1, 9 + e_flag)
        header = "vb\t\t# cb\t\tcir_pol\t\t|P+|^2\t\t|P-|^2\t\t|Px|^2\t\t|Py|^2\t\t|Pz|^2\t\t|P|^2"
        if e_flag : header += "\t\t  dE"

        if save : 
            if txt : np.savetxt(save+'.txt', optical, fmt = '%.4f', delimiter='\t\t', header=header)
            else :   np.save(save, optical)
                
        return "  vb \tcb \tcir_pol\t|P+|^2\t|P-|^2\t|Px|^2\t|Py|^2\t|Pz|^2\t|P|^2\t  dE\n" + \
                np.array2string(optical, formatter={'float_kind':lambda x: "%.2f" % x}, separator='\t', suppress_small=True)
    

    def L_Plot(self, kpoint, direction, vb_max, num_vb, num_cb, plot = True, save = '', txt = False):
    # Calculates orbital angular momentum (L)
    # from matrix elements 
    # and band energies read from vasp's WAVEDERF
    # for selected k point number
    # direction: 1 - X, 2 - Y, 3 - Z
        if   direction == 1 : beta = 2; gamma = 3
        elif direction == 2 : beta = 1; gamma = 3
        elif direction == 3 : beta = 1; gamma = 2
        else : raise Exception('wrong direction')

        val = [vb_max - num_vb + i + 1 for i in range(num_vb)]
        con = [vb_max + i + 1 for i in range(num_cb)]

        bands_oi = val + con
        nbands_oi = len(bands_oi)

        if self.P_matrix.size == 0 :
            bands_to_read = []
            for band in bands_oi:
                if band not in self.known_bands():
                    bands_to_read.append(band)
            self.read_bands(bands_to_read)
            P_vecs = self.band_vecs
            nbands = P_vecs[vb_max].shape[1]
        else : nbands = self.P_matrix.shape[1]

        L  = np.zeros((nbands, nbands_oi), dtype = complex)

        for count, band in enumerate(bands_oi):
            if self.P_matrix.size > 0 : P_vec = self.P_matrix[kpoint-1, band-1,:,:].copy()
            else :                      P_vec = P_vecs[band][kpoint-1].copy()
            P_vec[abs(P_vec[:,3]) < .00000001] = [0,0,0,1]
            L_iter = 2*np.imag(P_vec[:,beta-1]*np.conj(P_vec[:,gamma-1]))/P_vec[:,3]
            L[:,count] = np.cumsum(L_iter)

        C = PhysicalConstants
        if self.source == 'VASP': 
            convert =      C.e * 1E-20 * C.me * 4 * C.pi**2 / C.h**2
        else: 
            convert = C.Ry * (C.a0)**2 * C.me * 4 * C.pi**2 / C.h**2

        L = np.real(L*convert)

        if plot:
            plt.figure()
            plt.plot(L, '-')
            dirs = {1: 'X', 2: 'Y', 3: 'Z'}
            plt.title('K-point: ' + str(kpoint+1) + ', ' \
                    + 'Direction: ' + dirs[direction])
            plt.xlabel('number of bands')
            plt.ylabel('L (hbar)')
            plt.legend(['v'+str(i) for i in val]+['c'+str(i) for i in con] )
            if plot == 'show': plt.show()
            if plot == 'save': plt.savefig('L_Plot.png')
            else: plt.savefig(plot)

        if save : 
            if txt : np.savetxt(save+'.txt', np.concatenate([np.array(bands_oi).reshape((1,-1)), L], axis = 0))
            else :   np.save(save, np.concatenate([np.array(bands_oi).reshape((1,-1)), L], axis = 0))

        return L

    def read_EIGENVAL(self):
        with open(self.EIGENVAL) as eig:
            # Get relevant data
            for i in range(5):
                eig.readline()
            data = np.fromstring(eig.readline(), 
                                dtype=int, 
                                sep = ' ')
            eig.readline() #skip next empty line
            vb_max, nkpts, nbands = data

            Q = np.empty((nbands,3,nkpts))
            for kpnt in range(nkpts):
                eig.readline()
                Q[:,:,kpnt] = np.fromfile(eig,
                                        count = 3*nbands,
                                        sep = ' ').reshape(-1,3)
        
        eval, occ = Q[:,1:2,:], Q[:,2:3,:]

        self.eval = eval
        self.occ  = occ

        return eval, occ

    def read_WAVEDER(self):
        with FortranFile(self.WAVEDER, 'r') as file:
            nb_tot, nbands_cder, nkpts, ispin = file.read_record(dtype= np.int32)
            nodesn_i_dielectric_function = file.read_record(dtype= float)
            wplasmon = file.read_record(dtype= float).reshape(3,3)
            cder = file.read_record(dtype= np.complex64).reshape(3, ispin,nkpts,nbands_cder,nb_tot)

        # Format cder matrix to resemble P from WAVEDERF
        P = np.moveaxis(cder, [2, 4, 3, 0], [0, 1, 2, 3])[:,:,:,:,0]
        P = np.concatenate([P, np.zeros_like(P[:,:,:,0:1])], axis = 3)
        eval, occ = self.read_EIGENVAL()
        for kpnt in range(nkpts):
            eval_t = eval[:,:,kpnt]
            z = (np.zeros((nbands_cder, nbands_cder)) - eval_t.T + eval_t).reshape(nbands_cder, nb_tot,1)
            P[kpnt,:,:,:] *= np.abs(z)
            P[kpnt,:,:,3:4] = z

        eval, occ = eval[:,:,-1], occ[:,:,-1]

        np.save(self.name+"_P_matrix.npy", P)
        self.P_matrix = np.load(self.name+"_P_matrix.npy", mmap_mode='r')
        self.eval = eval

    def read_mommat2up(self):
        with open(self.mommat2up, 'r') as fh:
            Q = np.array([])
            fh.readline() 
            fh.readline() # skip first two lines
            nkpnts = 0

            while True:
                try:
                    # Following is a header for a given k-point
                    a = fh.readline().split()
                    nbands = int(a[6])  # number of bands
                    nlines = int(0.5*nbands*(nbands+1))
                    fh.readline() # skip next line
                    Q = np.append(Q, np.fromfile(fh, count = 9*nlines, sep = ' '))
                    nkpnts += 1
                except:
                    break
                
        dat = Q.reshape((-1, 9))
        P = np.zeros((nkpnts, nbands, nbands, 4), dtype=complex)

        for kpnt in range(nkpnts):
            start =  kpnt   *nlines
            stop  = (kpnt+1)*nlines

            mat_vec = np.concatenate(
                    [dat[start:stop, 2:3] + 1j*dat[start:stop, 3:4],
                    dat[start:stop, 4:5] + 1j*dat[start:stop, 5:6],
                    dat[start:stop, 6:7] + 1j*dat[start:stop, 7:8],
                    dat[start:stop, 8:9]],
                    axis = 1 )
            
            indices = np.triu_indices_from(P[kpnt,:,:,0], k =  0)
            P[kpnt,:,:,0][indices] = mat_vec[:,0]
            P[kpnt,:,:,0] += P[kpnt,:,:,0].conj().T - np.conj(np.diag(np.diag(P[kpnt,:,:,0])))

            P[kpnt,:,:,1][indices] = mat_vec[:,1]
            P[kpnt,:,:,1] += P[kpnt,:,:,1].conj().T - np.conj(np.diag(np.diag(P[kpnt,:,:,1])))

            P[kpnt,:,:,2][indices] = mat_vec[:,2]
            P[kpnt,:,:,2] += P[kpnt,:,:,2].conj().T - np.conj(np.diag(np.diag(P[kpnt,:,:,2])))

            P[kpnt,:,:,3][indices] =  mat_vec[:,3]
            P[kpnt,:,:,3] -= P[kpnt,:,:,3].T + np.diag(np.diag(P[kpnt,:,:,3]))

        np.save(self.name+"_P_matrix.npy", P)
        self.P_matrix = np.load(self.name+"_P_matrix.npy", mmap_mode='r')


    def read_WAVEDERF(self):
        with open(self.WAVEDERF, 'r') as fh:

            Q = np.fromstring(fh.readline(), dtype=int, 
                            sep = ' ')

            if Q.shape[0] != 3 : 
                raise Exception('error reading file')

            nkpnts = Q[1]
            nbands = Q[2]

            Q = np.fromfile(fh, sep = ' ')

        if Q.shape[0] != (12*nbands*nbands*nkpnts) : 
            raise Exception('error reading file')

        dat = Q.reshape((-1, 12))

        P = np.zeros((nkpnts, nbands, nbands, 4), dtype = complex)

        for kpnt in range(nkpnts) :
            start =  kpnt   *nbands**2
            stop  = (kpnt+1)*nbands**2

            eval = dat[start:start+nbands, 4:5]
            occ  = dat[start:start+nbands, 5:6]

            mat_vec = np.concatenate(
                    [dat[start:stop,  6: 7] + 1j*dat[start:stop,  7: 8],
                    dat[start:stop,  8: 9] + 1j*dat[start:stop,  9:10],
                    dat[start:stop, 10:11] + 1j*dat[start:stop, 11:12]],
                    axis = 1 )

            # Introduce z matrix for analogous operation to WAVEDER
            z = np.zeros((nbands, nbands)) - eval.T + eval
            
            P[kpnt,:,:,0] = mat_vec[:,0:1].reshape(nbands, nbands) *abs(z)
            P[kpnt,:,:,1] = mat_vec[:,1:2].reshape(nbands, nbands) *abs(z)
            P[kpnt,:,:,2] = mat_vec[:,2:3].reshape(nbands, nbands) *abs(z)
            P[kpnt,:,:,3] = z

        np.save(self.name+"_P_matrix.npy", P)
        self.P_matrix = np.load(self.name+"_P_matrix.npy", mmap_mode='r')
        self.eval = eval


    def read_procar(self):
        try: 
            with open(self.PROCAR+'.npy', 'rb') as f:
                promat = np.load(f, allow_pickle=True)
                energs = np.load(f, allow_pickle=True)
                E_max  = np.load(f, allow_pickle=True)
                E_min  = np.load(f, allow_pickle=True)
                E_gap  = np.load(f, allow_pickle=True)
                kindexer = np.load(f, allow_pickle=True)
                print('Read from PROCAR.npy')
        except:
            print(f"Checking {self.PROCAR}")
            nlines = 0
            nprjct = 0
            with open(self.PROCAR) as procar:
                ion_flag = False
                while True:
                    data = procar.readline().split()
                    nlines += 1
                    if len(data) > 0 :
                        if data[0] == 'ion' and ion_flag == False:
                            norbitals = len(data)-2
                            ion_flag = True
                        elif data[0] == 'tot' and ion_flag == True:
                            nprjct += 1
                        elif data[0] == 'ion' and ion_flag == True:
                            break
                for line in procar:
                    if "# of k-points:" in line:
                        nprjct += 1
                    nlines += 1

            print(f"Reading {self.PROCAR}", end=' ')
            with open(self.PROCAR) as procar:
                procar.readline()
                data = procar.readline().split()
                nkpts, nbands, nions = int(data[3]), int(data[7]), int(data[11])
                kindexer = {}
                proj = 0
                promat = np.zeros((nkpts, nbands, nprjct, nions, norbitals), dtype=float)
                energs = np.zeros((nkpts, nbands), dtype=float)
                vbmax = 0
                for line in tqdm(procar, total=nlines - ((nions+1)*nprjct)*nbands*nkpts - 2):
                    if "# of k-points:" in line:
                        proj += 1
                    data = line.split()
                    if len(data) > 0:
                        if data[0] == 'k-point': 
                            kpoint = int(data[1])-1
                            kindexer[tuple(float(f'{float(num):.4f}') for num in data[3:6])] = kpoint
                        if data[0] == 'band':
                            band = int(data[1])-1
                            energs[kpoint, band] = float(data[4])
                            if vbmax==0 and int(float(data[7]))==0:
                                vbmax = band
                        if data[0] == 'ion':
                            if nprjct != 2:
                                for proj in range(nprjct):
                                    for ion in range(nions):
                                        mrow = np.fromstring(procar.readline(),
                                                            dtype=float, sep=' ')[1:-1]
                                        promat[kpoint, band, proj, ion, :] = mrow
                                    procar.readline()
                            else:
                                for ion in range(nions):
                                    mrow = np.fromstring(procar.readline(),
                                                        dtype=float, sep=' ')[1:-1]
                                    promat[kpoint, band, proj, ion, :] = mrow
                                procar.readline()
            print("PROCAR read, promat shape is:", promat.shape, f'vbmax is {vbmax}')
            E_max = np.max(energs[:,vbmax-1])
            E_min = np.min(energs[:,vbmax])
            E_gap = E_min - E_max
            try:
                with open(self.PROCAR+'.npy', 'wb') as f:
                    np.save(f, promat, allow_pickle=True)
                    np.save(f, energs, allow_pickle=True)
                    np.save(f, E_max, allow_pickle=True)
                    np.save(f, E_min, allow_pickle=True)
                    np.save(f, E_gap, allow_pickle=True)
                    np.save(f, kindexer, allow_pickle=True)
                    print('Saved to PROCAR.npy')
            except: pass
            return promat, energs, E_max, E_min, E_gap, kindexer


    def band_plot(self, ax=None,
                kpoints = [], bands = [], projections = [], ions = [], orbits = [], 
                fatbands = True, bandlines = False, sum_orbits = False, sum_ions = True, 
                label = '', bandcolor='black', interpolate = False):
        promat, energs, E_max, E_min, E_gap, kindexer = self.read_procar()
        energs -= E_max

        if ax == None: ax = plt.gca()

        klist = []
        for pos in kindexer.keys():
            if pos[0] == 0 and pos[2] == 0 and pos[1] >= 0 and pos[1] <= 0.5:
                klist.append(kindexer[pos])
        klist = klist[-1::-1]
        print(klist)
        for pos in kindexer.keys():
            if pos[1] == 0 and pos[2] == 0 and pos[0] > 0 and pos[0] <= 0.5:
                klist.append(kindexer[pos])
        kpoints = klist

        ebands = np.array([]).reshape(0, energs.shape[1])
        pbands = np.array([]).reshape(0, *promat.shape[1:])
        orbit_names = ['$s$', '$p_y$', '$p_z$', '$p_x$', '$d_{xy}$', 
                '$d_{yz}$', '$d_{z^2}$', '$d_{xz}$', '$d_{x^2-y^2}$']

        if orbits == []: orbits = np.arange(9)
        elif orbits == 's': orbits = [1]
        elif orbits == 'p': orbits = [2,3,4]
        elif orbits == 'd': orbits = [5,6,7,8,9]
        else: orbits = np.array(orbits)-1
        if ions == []: ions = np.arange(promat.shape[3])+1
        if projections == []: projections = np.arange(promat.shape[2])+1
        cmap = colormaps['Set1']
        colors = cmap(np.linspace(0, 1, 9))
        x_line = np.linspace(-0.5,0.5,len(kpoints))
        for kpoint in kpoints:
            ebands = np.concatenate([ebands, energs[kpoint:kpoint+1]], axis = 0)
            pbands = np.concatenate([pbands, promat[kpoint:kpoint+1]], axis = 0)
        if len(bands) == 0:
            if bandlines: ax.plot(x_line, ebands, color=bandcolor, linewidth = 0.5)
        elif bandlines: ax.plot(x_line, ebands[:,bands], color=bandcolor, linewidth = 0.5)
        if sum_ions:
            for ion in ions[1:]:
                pbands[:,:,:,ions[0]-1,:] += pbands[:,:,:,ion-1,:]
            ions = [ions[0]]
        if fatbands:
            legend_lines = []
            legend_labels= []
            for ion in ions:
                for proj in projections:
                    if sum_orbits:
                        ax.scatter(np.stack([x_line]*energs.shape[1], axis=1), 
                                    ebands, np.sum(np.abs([100*pbands[:,:,proj-1,ion-1,j] for j in orbits]), axis=0), 
                                    label = label, color = colors[1])
                    else:
                        for j in orbits:
                            legend_lines.append(Line2D([0], [0], color=colors[j], lw=4, label = orbit_names[j]))
                            legend_labels.append(orbit_names[j])
                            if interpolate:
                                for band in bands:
                                    x_int = np.linspace(x_line[0], x_line[-1], interpolate)
                                    e_int = np.interp(x_int, x_line, ebands[:,band])
                                    p_int = np.interp(x_int, x_line, 100*pbands[:,band,proj-1,ion-1,j])
                                    ax.scatter(x_int, e_int, p_int, 
                                            #    label = orbit_names[j], 
                                            color = colors[j], alpha=0.08)
                            else:
                                if len(bands) == 0:
                                    ax.scatter(np.stack([x_line]*energs.shape[1], axis=1), 
                                            ebands, np.abs(100*pbands[:,:,proj-1,ion-1,j]), 
                                            color = colors[j])
                                else:
                                    ax.scatter(np.stack([x_line]*len(bands), axis=1), 
                                            ebands[:,bands], 100*pbands[:,bands,proj-1,ion-1,j], 
                                            color = colors[j])
            ax.legend(legend_lines, legend_labels)
        ax.axvline(0, color="black", linestyle="--", linewidth = 0.5)
        ax.axhline(0, color="black", linestyle="--", linewidth = 0.5)


    def read_band_out(self):
        """
        Read band structure from DFTB+ band.out file.

        Format:
        KPT i spin j kweight w
            band_index   energy   occupation
            ...

        Stores:
            self.eval -> ndarray shape (nkpts, nbands)
            self.occ  -> ndarray shape (nkpts, nbands)
        """
        kpoints = []
        all_energies = []
        all_occs = []

        with open(self.bandout, "r") as f:
            lines = f.readlines()

        energies, occs = [], []
        for line in lines:
            if line.startswith("KPT"):
                # save previous block
                if energies:
                    all_energies.append(energies)
                    all_occs.append(occs)
                    energies, occs = [], []

                parts = line.split()
                ik = int(parts[1])
                spin = int(parts[3])
                weight = float(parts[5])
                kpoints.append({"ik": ik, "spin": spin, "weight": weight})

            else:
                vals = line.split()
                if len(vals) >= 3:
                    band_index = int(vals[0])   # not used but could check
                    energies.append(float(vals[1]))
                    occs.append(float(vals[2]))

        # append last block
        if energies:
            all_energies.append(energies)
            all_occs.append(occs)

        self.eval = np.array(all_energies)
        self.occ = np.array(all_occs)
        self.kpoints = kpoints

        return self.eval, self.occ


    def band_plot_DFTB(self, ax=None, fermi_level=None, color="black", label="DFTB+"):
        """
        Plot DFTB+ band structure from band.out (requires self.eval, self.occ).
        """
        if self.eval is None or self.eval.size == 0:
            raise Exception("No band data loaded. Run read_band_out() first.")

        energies = self.eval
        occs = self.occ
        nkpts, nbands = energies.shape
        x = np.arange(nkpts)

        # auto-guess Fermi level if not given
        if fermi_level is None:
            fermi_level = np.max(energies[occs > 0.5]) if np.any(occs > 0) else 0.0

        if ax is None:
            ax = plt.gca()

        for ib in range(nbands):
            ax.plot(x, energies[:, ib] - fermi_level, color=color, linewidth=0.8)

        ax.axhline(0, color="red", linestyle="--", linewidth=0.8, label="Fermi level")
        ax.set_xlabel("k-path index")
        ax.set_ylabel("Energy (eV)")
        ax.set_title(label)

        return ax



def addCHG(nums, folder):
    nums = [str(num) for num in nums]
    with open(folder+'CHG.'+nums[0]) as file:
        head = ''
        for line in file:
            head += line
            if len(line.split()) == 0:
                break
        dims1 = file.readline()
        vec1 = np.fromstring(file.read(), sep=' ')

    with open(folder+'CHG.'+nums[1]) as file:
        for line in file:
            if len(line.split()) == 0:
                break
        dims2 = file.readline()
        vec2 = np.fromstring(file.read(), sep=' ')

    with open(folder+'CHG.'+'-'.join(nums), 'w') as file:
        file.write(head)
        file.write(dims2)
        vec = vec1+vec2
        rest = len(vec)%10
        vec1 = vec[:-rest].reshape((-1,10))
        vec2 = vec[-rest:].reshape((1,rest))
        np.savetxt(file, vec1, fmt='%.5e')
        np.savetxt(file, vec2, fmt='%.5e')