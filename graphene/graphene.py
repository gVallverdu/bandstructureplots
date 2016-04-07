#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys

import numpy as np
from numpy import array as npa

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital


def rgbline(ax, k, e, red, green, blue, alpha=1.):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=2)
    ax.add_collection(lc)

if __name__ == "__main__":
    # read data
    # ---------

    # kpoints labels
    path = HighSymmKpath(mg.Structure.from_file("./opt/CONTCAR")).kpath["path"]
    labels = [r"$%s$" % lab for lab in path[0][0:4]]

    # bands
    bands = Vasprun("./bands/vasprun.xml").get_band_structure("./bands/KPOINTS", line_mode=True)

    # projected bands
    data = Procar("./bands/PROCAR").data

    # density of state
    dosrun = Vasprun("./dos/vasprun.xml")

    # set up matplotlib plot
    # ----------------------

    # general options for plot
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of State
    gs = GridSpec(1, 2, width_ratios=[2, 1])
    fig = plt.figure(figsize=(11.69, 8.27))
    fig.suptitle("Bands diagram of graphene")
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])  # , sharey=ax1)

    # set ylim for the plot
    # ---------------------
    emin = 1e100
    emax = -1e100
    for spin in bands.bands.keys():
        for b in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][b]))
            emax = max(emax, max(bands.bands[spin][b]))

    emin -= bands.efermi + 1
    emax -= bands.efermi - 1
    ax1.set_ylim(emin, emax)
    ax2.set_ylim(emin, emax)

    # Band Diagram
    # ------------

    # sum up contribution over carbon atoms
    data = data[Spin.up].sum(axis=2)

    # sum up px and py contributions and normalize contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            sc = data[k][b][Orbital.s.value]**2
            pxpyc = data[k][b][Orbital.px.value]**2 + \
                data[k][b][Orbital.py.value]**2
            pzc = data[k][b][Orbital.pz.value]**2
            tot = sc + pxpyc + pzc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pxpyc / tot
                contrib[b, k, 2] = pzc / tot

    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax1,
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0],
                contrib[b, :, 1],
                contrib[b, :, 2])

    # style
    ax1.set_xlabel("k-points")
    ax1.set_ylabel(r"$E - E_f$   /   eV")
    ax1.grid()

    # fermi level at 0
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=2)

    # labels
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        ax1.vlines(i * step, emin, emax, "k")
    ax1.set_xticks([i * step for i in range(nlabs)])
    ax1.set_xticklabels(labels)

    ax1.set_xlim(0, len(bands.kpoints))

    # Density of state
    # ----------------

    ax2.set_yticklabels([])
    ax2.grid()
    ax2.set_xticks(np.arange(0, 1.5, 0.4))
    ax2.set_xticklabels(np.arange(0, 1.5, 0.4))
    ax2.set_xlim(1e-6, 1.5)
    ax2.hlines(y=0, xmin=0, xmax=1.5, color="k", lw=2)
    ax2.set_xlabel("Density of State")

    # s contribution
    ax2.plot(npa(dosrun.pdos[0][Orbital.s][Spin.up]) +
             npa(dosrun.pdos[1][Orbital.s][Spin.up]),
             dosrun.tdos.energies - dosrun.efermi,
             "r-", label="s", linewidth=2)

    # px py contribution
    ax2.plot(npa(dosrun.pdos[0][Orbital.px][Spin.up]) +
             npa(dosrun.pdos[1][Orbital.px][Spin.up]) +
             npa(dosrun.pdos[0][Orbital.py][Spin.up]) +
             npa(dosrun.pdos[1][Orbital.py][Spin.up]),
             dosrun.tdos.energies - dosrun.efermi,
             "g-",
             label="(px, py)",
             linewidth=2)

    # pz contribution
    ax2.plot(npa(dosrun.pdos[0][Orbital.pz][Spin.up]) +
             npa(dosrun.pdos[1][Orbital.pz][Spin.up]),
             dosrun.tdos.energies - dosrun.efermi,
             "b-", label="pz", linewidth=2)

    # total dos
    ax2.fill_between(dosrun.tdos.densities[Spin.up],
                     0,
                     dosrun.tdos.energies - dosrun.efermi,
                     color=(0.7, 0.7, 0.7),
                     facecolor=(0.7, 0.7, 0.7))

    ax2.plot(dosrun.tdos.densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             color=(0.6, 0.6, 0.6),
             label="total DOS")

    # plot format style
    # -----------------

    ax2.legend(fancybox=True, shadow=True, prop={'size': 18})
    plt.subplots_adjust(wspace=0)

    # plt.show()
    plt.savefig(sys.argv[0].strip(".py") + ".pdf", format="pdf")
