#!/usr/bin/env python3
# -- encoding: utf-8 --
import numpy
import matplotlib.pyplot as plt

from rvlm.labhelpers.timesignals import argmax

probes_solver = numpy.loadtxt('solver_output.txt')

probe1_solver = probes_solver[:, [0,1]]
probe2_solver = probes_solver[:, [0,2]]

probe1_cst = numpy.loadtxt('cst_probe1.txt', skiprows=2)
probe2_cst = numpy.loadtxt('cst_probe2.txt', skiprows=2)

plt.plot(probes_solver[:, 0], probes_solver[:, 1],
         probe1_cst[:, 0],    10*probe1_cst[:, 1])

plt.savefig("probe1.png", format='png')

plt.clf()
plt.plot(probes_solver[:, 0], probes_solver[:, 2],
         probe2_cst[:, 0],    probe2_cst[:, 1])

plt.savefig("probe2.png", format='png')

print(argmax(probe1_solver) - argmax(probe1_cst))
print(argmax(probe2_solver) - argmax(probe2_cst))

