
# %%
import numpy as np
import visualization_prx as vis
from matplotlib import pyplot as plt
from labellines import labelLines


figurepath = 'figures/isoFisherInformation/'

# Parameters:
T = 0.2
# A_two_point_start = 1.
w_two_point = 0.1
# A_one_point_start = 0.45
w_one_point = 0.07

g = 1. # Not stored in file name

xi = 0.2
mean_act = 0.15
N_x = 100

FisherInfo_wishlist = np.array([0.01, 0.015, 0.02])

datapath = 'data/Theory/FisherInfo_vs_Aone_Atwo/'

A_one_point_collection = []
A_two_point_collection = []
FisherInfo_const_collection = []

for ii, FisherInfo_wish in enumerate(FisherInfo_wishlist):
    with open( datapath + 'const_longer_FisherInfowish=' + str(FisherInfo_wish) + '_T=' + str(T) + '_f=' + str(mean_act) 
                        + '_w_one_point=' + str(w_one_point) + '_w_two_point=' + str(w_two_point) 
                        + '_xi='+str(xi) + '_N_x=' + str(N_x), 'rb') as f:

        FisherInfo_wish = np.load(f)
        A_one_point_list = np.load(f)
        A_two_point_list = np.load(f)

    
    A_one_point_collection.append(A_one_point_list)
    A_two_point_collection.append(A_two_point_list)
    FisherInfo_const_collection.append(FisherInfo_wish*np.ones(len(A_one_point_list)))

#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

ax.set_title(r'\textbf{b}',loc='left')

for ii in range(len(A_one_point_collection)):
    ax.plot(A_one_point_collection[ii]/T, A_two_point_collection[ii]/T, label='${\cal I}/N=$' + str(FisherInfo_wishlist[ii]), color = 'black')
plt.xlim([0.013/T,0.026/T])
plt.xticks([0.014/T, 0.018/T, 0.022/T, 0.026/T])
plt.ylim([0.,9./T])
plt.xlabel('feedforward strength $\mathrm{U}_{\mathrm{inp}}$')
plt.ylabel('recurrence strength $K_{\mathrm{rec}}$', rotation = 90, labelpad=10)
labelLines(ax.get_lines(), zorder=2.5)

figurename = ('isoFisherInfo_T=' + str(T) + '_w_two_point=' + str(w_two_point) + '_w_one_point='  + str(w_one_point)
              + '_g=' + str(g) + '_xi=' + str(xi) + '_N_x=' + str(N_x) + '_f=' + str(mean_act))

plt.savefig(figurepath + figurename + '.pdf', format='pdf')
plt.savefig(figurepath + figurename + '.eps', format='eps')
plt.savefig(figurepath + figurename + '.jpg', format='jpg')

#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

for ii in range(len(A_one_point_collection)):
    ax.plot(A_one_point_collection[ii]/T, A_two_point_collection[ii]/T, label='${\cal I}/N=$' + str(FisherInfo_wishlist[ii]), color = 'black')
plt.xlim([0.013/T,0.026/T])
plt.xticks([0.014/T, 0.018/T, 0.022/T, 0.026/T])
plt.ylim([0.,9./T])
plt.xlabel('strength of feedforward input')
plt.ylabel('strength of recurrence', rotation = 90, labelpad=10)
labelLines(ax.get_lines(), zorder=2.5)

figurename = ('isoFisherInfo_talk_T=' + str(T) + '_w_two_point=' + str(w_two_point) + '_w_one_point='  + str(w_one_point)
              + '_g=' + str(g) + '_xi=' + str(xi) + '_N_x=' + str(N_x) + '_f=' + str(mean_act))

plt.savefig(figurepath + figurename + '.pdf', format='pdf')
plt.savefig(figurepath + figurename + '.eps', format='eps')
plt.savefig(figurepath + figurename + '.jpg', format='jpg')
# %%
