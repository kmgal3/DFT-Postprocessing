import numpy as np
import matplotlib.pyplot as plt
from InformationCollector import InformationCollector
from Artists import MatrixSculptor, ProjectionPainter, SegmentArranger, KlineEngineer, BandWriter

def BandToGnuplot(Pmat, Emat, Kpts, Dvec, 
                  projection, bands, description, kpaths, 
                  folder=''
                  ):
    """Write a band structure plot to a gnuplot compatible file

    Args:
        Pmat (array_like): PROCAR matrix
        Emat (array_like): EIGENVAL matrix
        Kpts (array_like): kpoint_coords
        Dvec (array_like): distances between kpoints
        projection (dict): projection dictionary as output by ProjectionPainter
        bands (str): selected bands e.g. '10,12,14-18 22-26'
        description (str): the sum of which projections should be plotted e.g. '1,2 Mo z p dxy' means 1st and 2nd Mo atom, in the z direction, all p-orbitals, dxy-orbital 
        kpaths (str): assuming line mode in band structure calculation 
                      any ordering of those lines can be selected (including inversion).
        folder (str, optional): to which folder should the output be saved. If '' (default) uses vaspout.h5's folder. 
    """    
    orbit_names=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
    chosen_orbits = ''
    orbitals = projection['orbital'] if projection['orbital'] != ... else np.arange(Pmat.shape[-1])
    directions = projection['direction']
    ions = projection['ion']
    bands = bands if bands != ... else np.arange(Pmat.shape[1])

    for i in orbitals:
        chosen_orbits += ' '+orbit_names[i]
    Dvec = Dvec.copy()
    Dvec -= Dvec[0]
    with open(folder+description.replace(' ', '_'), 'w') as file:
        file.write(description+' '+kpaths+' '+(str(bands) if len(bands)<Pmat.shape[1] else '')+'\n')
        header = 'k point coordinates x y z             | Distance   | Energy     |'
        # file.write(' Dir 1 Ion 1 px py pz                 : Ion 2 px py pz\n')
        for band in bands:
            mat = np.concatenate([Kpts, Dvec[:,None], Emat[:,band][:,None]], axis = 1)
            for in_dir, direction in enumerate(directions):
                stub = f' Dir {in_dir+1} '
                for in_ion, ion in enumerate(ions):
                    stub += f': Ion {in_ion+1}' + chosen_orbits
                    stub += ' '*(len(orbitals)*13*(in_ion+1) - len(stub)-1)
                    mat = np.concatenate([mat, Pmat[:,band,direction,ion,orbitals]], axis = 1)
                if header: header += stub +'|'
            if header: file.write(header+'\n')
            header = ''
            np.savetxt(file, mat, '%12.8f')
            file.write('\n')


def DOSToGnuplot(Pdos, Edos, 
                  projection, description, 
                  folder=''
                  ):
    """Write a density of states plot to a gnuplot compatible file

    Args:
        Pdos (_type_): DOSCAR matrix
        Edos (_type_): electron_dos/energies
        projection (dict): projection dictionary as output by ProjectionPainter
        description (str): the sum of which projections should be plotted e.g. '1,2 Mo z p dxy' means 1st and 2nd Mo atom, in the z direction, all p-orbitals, dxy-orbital 
        folder (str, optional): to which folder should the output be saved. If '' (default) uses vaspout.h5's folder. 
    """    
    orbit_names=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
    chosen_orbits = ''
    orbitals = projection['orbital'] if projection['orbital'] != ... else np.arange(Pdos.shape[-1])
    directions = projection['direction']
    ions = projection['ion']

    for i in orbitals:
        chosen_orbits += ' '+orbit_names[i]
    with open(folder+description.replace(' ', '_'), 'w') as file:
        file.write(description+'\n')
        header = ' Energy     |'
        # file.write(' Dir 1 Ion 1 px py pz                 : Ion 2 px py pz\n')
        mat = Edos[:,None]
        for in_dir, direction in enumerate(directions):
            stub = f' Dir {in_dir+1} '
            for in_ion, ion in enumerate(ions):
                stub += f': Ion {in_ion+1}' + chosen_orbits
                stub += ' '*(len(orbitals)*13*(in_ion+1) - len(stub)-1)
                mat = np.concatenate([mat, Pdos[:,direction,ion,orbitals]], axis = 1)
            if header: header += stub +'|'
        if header: file.write(header+'\n')
        header = ''
        np.savetxt(file, mat, '%12.8f')
        file.write('\n')


def BandStructurePlot(ic: InformationCollector,
                      min_diff=.1, 
                      ax=None, 
                      E0=None, 
                      gnuplot=False, 
                      folder='', 
                      save=False,
                      bandnums=False, 
                      color = None, 
                      mult=1, 
                      alpha=.1,
                      **kwargs
                      ):
    """Make a band structure plot using data from InformationCollector object

    Args:
        ic (InformationCollector): _description_
        min_diff (float): if band numbers are to be displayed, how far apart they should be to not group them. Defaults to .1
        ax (plt.axes): to which axis should the plot be rendered. If None (default) uses plt.gca().     
        E0 (float | str): what energy should be used as 0. str options are ('fermi'). If None (default) will look for valence band maximum. 
        gnuplot (bool): should the plot be written to a file in gnuplot format. Defaults to False, 
        folder (str): to which folder should the output(s) be saved. If '' (default) uses vaspout.h5's folder. 
        save (bool): should the plot be immediately saved. Defaults to False.
        bandnums (bool): should band numbers be displayed next to the plot. Defaults to False. 
        color (int | float | str | None): same as in all pyplot figures. Defaults to None 
        mult (float): the dot size multiplier. Defaults to 1
        alpha (float): the dot opacity level in range [0, 1]. Defaults to .1

    Returns:
        _type_: _description_
    """    
    Pmat = SegmentArranger(ic.nSeg, ic.Pmat, ic.kpaths)
    Emat = SegmentArranger(ic.nSeg, ic.Emat, ic.kpaths)
    kline_args = (ic.Lvec, ic.Scal, ic.Kpts, ic.Klen, ic.Klab, ic.nSeg)
    Kpts, Dvec, Ktik, Klab = KlineEngineer(*kline_args, ic.kpaths)

    Evalmax  = ic.Evalmax
    Egap, ions, num_ions = ic.Egap, ic.ions, ic.num_ions
    description, kpaths, bands = ic.description, ic.kpaths, BandWriter(ic.bands)

    if not folder: folder = ic.folder

    if E0 == None: Emat -= Evalmax
    elif E0 == 'fermi': Emat -= ic.Efer
    else: Emat -= E0

    projection = ProjectionPainter(ions, num_ions, description)

    try: Emat = np.sum(Emat[..., projection['direction']], axis=-1)
    except IndexError: Emat = Emat[...,0]

    if gnuplot: BandToGnuplot(Pmat, Emat, Kpts, Dvec, projection, bands, description, kpaths,folder)

    Pmat = MatrixSculptor(Pmat, projection)
    if ax==None: ax = plt.gca()

    ### Band labelling
    if bandnums:
        diff = Emat[-1, 0]
        band_prev = 0
        for band in range(Pmat.shape[1]):
            diff = abs(Emat[-1, band] - diff)
            if diff > min_diff:
                if band-1 != band_prev:
                    ax.annotate(str(band_prev+1)+'-'+str(band-1+1), xy = (Dvec[-1], Emat[-1, band-1]))
                else: ax.annotate(str(band_prev+1), xy = (Dvec[-1], Emat[-1, band-1]))
                band_prev = band
            diff = Emat[-1, band]
        band = Pmat.shape[1]
        if band-1 != band_prev:
            ax.annotate(str(band_prev+1)+'-'+str(band-1+1), xy = (Dvec[-1], Emat[-1, band-1]))
        else: ax.annotate(str(band_prev+1), xy = (Dvec[-1], Emat[-1, band-1]))

    Emat = Emat[:,bands]
    Pmat = Pmat[:,bands]

    if 'clean' not in description:
        ax.plot(Dvec, Emat, color="black", linewidth=0.5, linestyle=':')
        ax.scatter(np.stack([Dvec]*Emat.shape[1], axis=1), Emat, np.abs(Pmat)*100*mult, 
                   alpha=alpha, label=description, c=color)
    else: ax.plot(Dvec, Emat, color=color, linewidth=1, label=description)

    ax.vlines(Ktik[1:-1], [Emat.min()]*(len(Ktik)-2), [Emat.max()]*(len(Ktik)-2),
            color="black", linestyle=":")
    ax.axhline(0, color="black", linestyle=":")
    if E0==None:
        ax.axhline(Egap, color="black", linestyle=":")
        ax.annotate(f'{Egap:.3f}', xy = (Dvec[0], Egap))

    ax.set_xticks(Ktik, Klab)
    ax.set_xlim(Dvec[0], Dvec[-1])
    if save: plt.savefig(folder+description.replace(' ', '_'))

    # return Evalmax, Egap


def DensityOfStatesPlot(ic: InformationCollector,
                        ax=None, E0=None, gnuplot=False,
                        **kwargs
                        ):
    """Make a density of states plot using data from InformationCollector object

    Args:
        ic (InformationCollector): _description_
        ax (plt.axes): to which axis should the plot be rendered. If None (default) uses plt.gca(). 
        E0 (float | str): what energy should be used as 0. str options are ('fermi'). If None (default) will look for valence band maximum. 
    """    
    Vdos, Pdos, Edos, Efer = ic.Vdos, ic.Pdos, ic.Edos, ic.Efer
    ions, num_ions, description = ic.ions, ic.num_ions, ic.description
    Emin, Emax = ic.Emin, ic.Emax

    if E0 == None: 
        Edos -= Efer
        Emin -= Efer
        Emax -= Efer
    else: 
        Edos -= E0
        Emin -= E0
        Emax -= E0

    projection = ProjectionPainter(ions, num_ions, description)

    if gnuplot: DOSToGnuplot(Pdos, Edos,
                             projection, description, ic.folder)

    Pdos = MatrixSculptor(Pdos, projection)

    Vdos = np.sum(Vdos[..., projection['direction']], axis=-1)

    if ax==None: ax = plt.gca()
    if description in ['total', 'clean']:
        ax.plot(Vdos, Edos, label='total')
        ax.set_xlim([np.percentile(Vdos, .1),np.percentile(Vdos, 99.9)])
    else:
        ax.plot(Pdos, Edos, label=description)
        ax.set_xlim([np.percentile(Pdos, .1),np.percentile(Pdos, 99.9)])

    
    ax.set_ylim([Emin, Emax])