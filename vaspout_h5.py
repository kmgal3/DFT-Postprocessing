import h5py
import numpy as np
import matplotlib.pyplot as plt
from Plotting import BandStructurePlot, DensityOfStatesPlot
from InformationCollector import InformationCollector

# This script can be used to produce the plots straight from the terminal.
# For the full description go to the bottom.

class vaspout_h5:
    """This class knows how to extract relevant information from vaspout.h5
    and put it into an InformationCollector object
    """
    def __init__(self, vaspout: str = 'vaspout.h5'):
        self.vaspout = h5py.File(vaspout, 'r')
        self.ic = InformationCollector()

        # Get PATH to directory for conveniance
        try: self.ic.folder = vaspout[:vaspout.rindex('/')+1]
        except ValueError: self.ic.folder = ''


    def extract_magmom(self):
        """Get MAGMOM info and put it into a separate file in a format accepted by INCAR
        """        
        vaspout = self.vaspout
        magmom_matrix: np.ndarray = vaspout['/intermediate/ion_dynamics/magnetism/spin_moments/values'][:]

        # Get the full MAGMOM matrix and sum over orbitals
        if magmom_matrix.shape[1]==1: # ISPIN = 1
            print("No magmom in this calculation")
        elif magmom_matrix.shape[1]==2: # ISPIN = 2
            magmom = np.sum(magmom_matrix[0,1,...], axis=1)
        else: # lnoncolinear = .true.
            x = np.sum(magmom_matrix[0,1,...], axis=1)[:,None]
            y = np.sum(magmom_matrix[0,2,...], axis=1)[:,None]
            z = np.sum(magmom_matrix[0,3,...], axis=1)[:,None]
            magmom = np.concatenate([x,y,z], axis=1)

        # write the MAGMOM file in the INCAR format
        with open(f"{self.ic.folder}MAGMOM", "w") as f:
            f.write("MAGMOM = ")
            for i, row in enumerate(magmom):
                # Format row with aligned spacing
                row_str = " ".join(f"{x:12.8f}" for x in row)
                if i == 0:
                    f.write(f"{row_str} \\\n")
                elif i == len(magmom)-1:
                    f.write(f"         {row_str}\n")
                else:
                    f.write(f"         {row_str} \\\n")


    def read_BS(self):
        """Put all relevant band structure information from vaspout to self.ic: InformationCollector
        """        
        vaspout = self.vaspout
        ic = self.ic
        ic.Pmat = vaspout['results/projectors/par'][:]
        ic.Emat = vaspout['results/electron_eigenvalues/eigenvalues'][:]
        ic.Occp = vaspout['results/electron_eigenvalues/fermiweights'][:]
        ic.nSeg = vaspout['input/kpoints/number_line_segments'][()]/2
        ic.Lvec = vaspout['input/poscar/lattice_vectors'][:]
        ic.Scal = vaspout['input/poscar/scale'][()]
        ic.Kpts = vaspout['results/electron_eigenvalues/kpoint_coords'][:]
        ic.Klen = vaspout['input/kpoints/number_kpoints'][()]
        ic.Klab = vaspout['input/kpoints/labels_kpoints'].asstr()[:]
        ic.Efer = vaspout['results/electron_dos/efermi'][()]
        ic.ions = vaspout['input/poscar/ion_types'].asstr()[:]
        ic.num_ions = vaspout['input/poscar/number_ion_types'][:]

        # Rearrange PROCAR and EIGENVAL matrices
        # into the desired format: ( ..., ndirs, nions, norbits )
        # and calculate bandgap
        ic.Pmat = np.moveaxis(ic.Pmat, [-2,-1], [0,1])
        ic.Evalmax = np.max(ic.Emat[ic.Occp>1e-10])
        ic.Econmin = np.min(ic.Emat[ic.Occp<=1e-10])
        ic.Egap = ic.Econmin-ic.Evalmax
        ic.Emat = np.moveaxis(ic.Emat, [-2,-1], [0,1])


    def read_DOS(self):
        """Put all relevant density of states information from vaspout to self.ic: InformationCollector
        """        
        vaspout = self.vaspout
        ic = self.ic

        ic.Vdos = vaspout['results/electron_dos/dos'][:]
        ic.Pdos = vaspout['results/electron_dos/dospar'][:]
        ic.Edos = vaspout['results/electron_dos/energies'][:]
        ic.Efer = vaspout['results/electron_dos/efermi'][()]
        ic.ions = vaspout['input/poscar/ion_types'].asstr()[:]
        ic.num_ions = vaspout['input/poscar/number_ion_types'][:]

        # Rearrange the matrices to desired format
        ic.Vdos = np.moveaxis(ic.Vdos, [-1], [0])
        ic.Pdos = np.moveaxis(ic.Pdos, [-1], [0])

        # Calculate energy range
        ic.Emin = np.min(ic.Edos)
        ic.Emax = np.max(ic.Edos)
        ic.Ediff = ic.Emax - ic.Emin
    

    def plot_BS(self,  
                description: str = '', 
                kpaths: str = '', 
                bands: str = '', 
                **kwargs
                ):
        """Make a band structue plot according to string specifiers 
        and other args that BandStructurePlot accepts

        Args:
            description (str, optional): the sum of which projections should be plotted. Defaults to ''.
            kpaths (str, optional): assuming line mode in band structure calculation 
                                    any ordering of those lines can be selected (including inversion). Defaults to ''.
            bands (str, optional): which bands are to be displayed in the form of an intuitive string. Defaults to ''.

        Kwargs:
            ax (plt.axes): to which axis should the plot be rendered. If None (default) uses plt.gca(). 
            
            E0 (float | str): what energy should be used as 0. str options are ('fermi'). If None (default) will look for valence band maximum. 
            
            folder (str): to which folder should the output(s) be saved. If '' (default) uses vaspout.h5's folder. 
            
            save (bool): should the plot be immediately saved. Defaults to False.
            
            gnuplot (bool): should the plot be written to a file in gnuplot format. Defaults to False, 
            
            bandnums (bool): should band numbers be displayed next to the plot. Defaults to False. 
            
            min_diff (float): if band numbers are to be displayed, how far apart they should be to not group them. Defaults to .1
            
            color (int | float | str | None): same as in all pyplot figures. Defaults to None 
            
            mult (float): the dot size multiplier. Defaults to 1
            
            alpha (float): the dot opacity level in range [0, 1]. Defaults to .1
        """        
        self.read_BS()
        ic = self.ic
        ic.description, ic.arrangement, ic.bands = description, kpaths, bands
        BandStructurePlot(ic, **kwargs)


    def plot_BS_default(self):
        """The default plot, gives overview of the entire band structure
        """        
        plt.figure(figsize=(4,8))
        self.plot_BS('clean', bandnums=True, E0='fermi')
        plt.tight_layout()
        plt.savefig(self.ic.folder+'BSplot')
    

    def plot_DOS(self, 
                 description='', 
                 ax=None, 
                 E0=None,
                 gnuplot = False,
                 **kwargs
                 ):
        """Make a density of states plot according to string specifiers 
        and other args that BandStructurePlot accepts

        Args:
            description (str, optional): the sum of which projections should be plotted. Defaults to ''.
            ax (plt.axes): to which axis should the plot be rendered. If None (default) uses plt.gca(). 
            E0 (float | str): what energy should be used as 0. str options are ('fermi'). If None (default) will look for valence band maximum. 
            gnuplot (bool): should the plot be written to a file in gnuplot format. Defaults to False, 
        """        
        self.read_DOS()
        self.ic.description = description
        DensityOfStatesPlot(self.ic, ax=ax, E0=E0, gnuplot=gnuplot, **kwargs)


    def plot_DOS_default(self):
        """The default plot, shows the entire density of states
        """        
        plt.figure(figsize=(4,8))
        self.plot_DOS('total',E0=0)
        plt.tight_layout()
        plt.savefig(self.ic.folder+'DOSplot')
    

def plot_BS_DOS(BS:vaspout_h5, 
                DOS:vaspout_h5, 
                description='', 
                arrangement='', 
                bands= ''
                ):
    """A simple function to produce an elegant band structure and density of states plot

    Args:
        BS (vaspout_h5): _description_
        DOS (vaspout_h5): _description_
        description (str, optional): _description_. Defaults to ''.
        arrangement (str, optional): _description_. Defaults to ''.
        bands (str, optional): _description_. Defaults to ''.
    """    
    fig, axes = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
    DOS.read_DOS()
    BS.plot_BS(description, arrangement, bands, 
               min_diff=DOS.ic.Ediff*1e-2, ax=axes[0], bandnums=True)
    DOS.plot_DOS(description, ax=axes[1], E0=BS.ic.Evalmax)
    plt.tight_layout()


# This script can be used to produce the plots straight from the terminal.
if __name__ == "__main__":
    import argparse

    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Make band stucture or density of states plot')
    
    # Required positional argument
    parser.add_argument('calc_type', type=str,
                        help='Either band or dos')

    parser.add_argument('description', type=str, nargs='?', default='',
                        help='The description of specific projections to be processed, e.g. "1,2 Mo x p dxy"')

    parser.add_argument('--kpaths', type=str, default='',
                        help='Order and direction of kpaths for band, e.g. "1 2 3 -4"')
    
    parser.add_argument('--bands', type=str, default='',
                        help='choice of bands for band')
    
    parser.add_argument("--min_diff", type=float, default=0.1,
                        help="Minimum difference threshold (default: 0.1).")

    parser.add_argument("--ax", default=None,
                        help="Matplotlib Axes object (default: None).")

    parser.add_argument("--E0", type=float, default=None,
                        help="Reference energy level (default: None).")

    parser.add_argument("--gnuplot", action="store_true",
                        help="Enable gnuplot output (default: False).")

    parser.add_argument("--folder", type=str, default="",
                        help="Output folder (default: '').")

    parser.add_argument("--save", action="store_true",
                        help="Save plots to file instead of showing (default: False).")

    parser.add_argument("--bandnums", action="store_true",
                        help="Show band numbers (default: False).")

    parser.add_argument("--color", type=str, default=None,
                        help="Color for plotting (default: None).")

    parser.add_argument("--mult", type=float, default=1,
                        help="Scaling multiplier (default: 1).")

    parser.add_argument("--alpha", type=float, default=0.1,
                        help="Transparency level (default: 0.1).")
    
    args = parser.parse_args()

    calc = vaspout_h5()
    plt.figure(figsize=(4,8))

    if args.calc_type == 'band':
        calc.plot_BS(**vars(args))

    if args.calc_type == 'dos':
        calc.plot_DOS(**vars(args))

    plt.tight_layout()
    plt.savefig(f'{calc.ic.folder}{args.calc_type}_plot')