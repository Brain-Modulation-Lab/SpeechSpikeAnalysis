function [spk_h, mean_h, err_h] = plotRastersandIFRs(time, ifr_ax, ifr, ifr_err, spk_ax, spikes, spk_offset, pcolor)
% function plotRastersandIFRs(time, ifr_ax, ifr, ifr_err, spike_ax, spikes, spk_offset)
%
% This is a general use function that will plot a spike raster and the
% corresponding instantaneous firing rate. The calling function should
% create the axes ahead of time, and provide them as arguments.

spk_h = plotSpikeRasters(spk_ax, time, spikes, spk_offset);
set(spk_h, 'Color', pcolor, 'LineWidth', 1);

[mean_h, err_h] = plot_err_poly(ifr_ax, time, ifr, ifr_err, pcolor, (pcolor + [1 1 1])/2, .5);
%mean_h = line('Parent', ifr_ax,'Xdata', time, 'Ydata', ifr); 