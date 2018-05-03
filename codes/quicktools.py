#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Feb, 15 2018
Author: Aleksander Molak || aleksander.molak@gmail.com
"""
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def visualize(what, title, figsize=(5,5)):
    f = plt.figure(figsize=figsize)
    ax = plt.gca()
    im = plt.imshow(what, clim=[-0.5,0.5],cmap='coolwarm')
    im = plt.imshow(what, clim=[-1,1],cmap='coolwarm')
    plt.title(title, fontsize=20)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.show()
