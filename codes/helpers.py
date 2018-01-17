#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 14:19:30 2018

@author: dabrowska
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def figure(grid = [1, 1], size = [8., 6.], top=0.9, bottom=0.1, left=0.1, right=0.98, hspace=0.25, wspace = 0.25):
    '''
    Creates a figure and a gridspace. Allows for nesting subplots in grid objects defined this way. 
    
    grid - a list [rows, columns] indicating the number of subplots
         - [1, 1] gives a single subplot
    
    defaults: 
    figure size = [8., 6.] (inches)
    top = 0.9, 
    bottom = 0.1
    left = 0.1
    right = 0.98
    hspace = 0.25
    wspace = 0.25
    
    Returns figure and gridspace objects. 
    '''
    fig = plt.figure(figsize = size, dpi=110)
    gs = gridspec.GridSpec(grid[0], grid[1])
    gs.update(top = top, bottom = bottom, left = left, right = right, 
              hspace = hspace, wspace = wspace)
    
    
    return fig, gs

def labels(ax, fontsize = 12, title = None, xlabel = None, ylabel = None, 
           legend = False, color = 'k', xcolor = 'k', ycolor = 'k', 
           suptitle = None, fig = None):
    '''
    Generates all demanded labels with fixed size ratios. 
    
    title, suptitle, xlabel, ylabel - strings
    color, xcolor, ycolor - colors of the (sup)title and axes' labels
    legend - if True, places the legend in 'best' location (labels have to be already provided while plotting)

    To add a suptitle provide a figure handle 'fig'.
    '''
    if title: 
        ax.set_title(title, fontsize = fontsize + 1, color = color)

    if suptitle: 
        if fig:
            fig.suptitle(suptitle, fontsize = fontsize + 1, color = color)
        else:
            print 'Provide a figure handle to add a suptitle'
            
    if xlabel:
        ax.set_xlabel(xlabel, fontsize = fontsize, color = xcolor)
    
    if ylabel:
        ax.set_ylabel(ylabel, fontsize = fontsize, color = ycolor)
        
    if legend:
        ax.legend(loc = 'best')
    
    if xcolor == ycolor:
        ax.tick_params('both', labelsize = fontsize - 1, color = xcolor)
    else:
        ax.tick_params('x', labelsize = fontsize - 1, colors = xcolor)
        ax.tick_params('y', labelsize = fontsize - 1, colors = ycolor)
        # print 'OK but this does not work => fontsize has no effect here and above'

def savef(fig, fname, dpis, formats):
    '''
    Saves the figure with given name. 
    
    fname - name of a file (full directory)
    dpis - a list of dpi values to be used when saving the figure
    formats - a list of strings indicating all formats in which the figure 
        is to be saved (should work with png, pdf, ps, eps and svg, 
        depending on an active backend)
    '''
    for d in dpis:
        for f in formats:
            if f in ['png', 'pdf']:
                fig.savefig(fname + '_' + str(d) + 'dpi.' + f, bbox_inches = 'tight')
            elif f == 'eps':
                fig.savefig(fname + '_' + str(d) + 'dpi.' + f, bbox_inches = 'tight', frameon = False)
            else:
                fig.savefig(fname + '_' + str(d) + 'dpi.' + f)
