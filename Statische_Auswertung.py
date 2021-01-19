#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:49:45 2020

@author: valentin
"""

import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as sma
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import numpy as np
import matplotlib.pyplot as plt
import os

nBlades = 3
nEig = 10
nRPM = 10
StaticAnalysisType = "linear" #"non_linear"
cleanupfolders = True

if cleanupfolders:
    os.system("rm -r EigenResults*")
    os.system("rm -r StaticResults*")
    os.system("rm -r fig")

# Modify postprocessing of prestress to extract maximum stress
class StructuralMechanicsAnalysisStatic(StructuralMechanicsAnalysis):
    def Finalize(self):
        mpName = self.project_parameters["solver_settings"]["model_part_name"].GetString()
        mp = self.model[mpName]
        self.stressM = [[]]*len(mp.Elements)
        for ii in range(len(mp.Elements)):
            self.stressM[ii] = mp.Elements[ii+1].CalculateOnIntegrationPoints(sma.VON_MISES_STRESS,
                                                                                  mp.ProcessInfo)[0]

# Modify postprocessing to extract mass and eigenfrequencies (as eigenvalues, eigenfrequencies in rad/s and in Hz)
class StructuralMechanicsAnalysisModal(StructuralMechanicsAnalysis):
        def Finalize(self):
            super(StructuralMechanicsAnalysisModal, self).Finalize()
            mpName = self.project_parameters["solver_settings"]["model_part_name"].GetString()
            mp = self.model[mpName]
            # mass of structure via nodmal mass
            mass_process = sma.TotalStructuralMassProcess(mp)
            mass_process.Execute()
            self.mass = mp.ProcessInfo[km.NODAL_MASS]
            # modal results
            self.Lambda = np.abs(mp.ProcessInfo[sma.EIGENVALUE_VECTOR])
            self.Phi = mp.GetNode(1).GetValue(sma.EIGENVECTOR_MATRIX)
            self.omegan = np.sqrt(self.Lambda)
            self.fn = self.omegan/(2*np.pi)


def Analysis(Omega):
    # Static analysis for pre-stress
    print("\nbegin static analysis ------------------------------------------")
    with open("ProjectParametersPrestressed.json", "r") as parameter_file:
        StaticParameters = km.Parameters(parameter_file.read())
    model = km.Model()                                                                  #km 
    StaticParameters["solver_settings"]["analysis_type"].SetString(StaticAnalysisType)
    StaticParameters["processes"]["loads_process_list"][0]["Parameters"]["modulus"].SetString("((x**2+y**2)**0.5)*"+str(((Omega/60)*2*np.pi)**2))
    StaticParameters["processes"]["loads_process_list"][0]["Parameters"]["direction"][0].SetString("x/((x**2+y**2)**0.5)")
    StaticParameters["processes"]["loads_process_list"][0]["Parameters"]["direction"][1].SetString("y/((x**2+y**2)**0.5)")
    #StaticParameters["processes"]["loads_process_list"][0]["Parameters"]["modulus"].SetDouble(Omega)
    StaticSimulation = StructuralMechanicsAnalysisStatic(model, StaticParameters)
    StaticSimulation.Run()
    print("\nend static analysis --------------------------------------------")
    os.system("mv vtk_output StaticResults_"+str(Omega))

    # Modal analysis
    print("\nbegin modal analysis -------------------------------------------")
    with open("ProjectParametersModal.json", "r") as parameter_file:
        ModalParameters = km.Parameters(parameter_file.read())
    ModalParameters["solver_settings"]["eigensolver_settings"]["number_of_eigenvalues"].SetInt(nEig)
    ModalSimulation = StructuralMechanicsAnalysisModal(model, ModalParameters)
    ModalSimulation.Run()
    print("\nend modal analysis ---------------------------------------------")
    os.system("mv EigenResults EigenResults_"+str(Omega))
    return(StaticSimulation, ModalSimulation)

#RPMi = 12.1
RPMList = np.linspace(0, 1500, nRPM)
#RPMList = [0] #die Belastung nur ein Mal durchführen
fnList = []
sigmaMax = []
for i, RPMi in enumerate(RPMList):
    StaticSimulation, ModalSimulation = Analysis(RPMi)
    fnList.append(ModalSimulation.fn)
    sigmaMax.append(np.max(StaticSimulation.stressM))

print(sigmaMax)


folderFig = "fig"
try:
    os.system("mkdir "+folderFig)
except:
    pass

#LabelList = ["$\omega_{"+str(i)+"}$" for i in list(range(1, len(fnList)+1))]
for i in range(len(np.array(fnList)[0, :])):
    plt.plot(RPMList, np.array(fnList)[:, i], label="$\\omega_{"+str(i+1)+"}$")
plt.plot(RPMList, RPMList, label="$\\Omega$")
#plt.plot(RPMList, RPMList*nBlades, label="$"+str(nBlades)+"\\Omega$")
plt.ylabel("Eigenfrequency $\\omega_0$ [Hz]")
plt.xlabel("Rotor speed $\\Omega$ [$\\frac{1}{\\mathrm{min}}$]")
#plt.axvline(x=6.0, label='Cut-in rot speed')
plt.axvline(x=1500, label='Nominal rot speed')
plt.legend(frameon=False)
plt.ylim(0, np.nanmax(np.array(fnList)))
plt.xlim(0, RPMList[-1])
plt.savefig("fig/Campbell.svg")
plt.savefig("fig/Campbell.png")
plt.show()

plt.plot(RPMList, sigmaMax)
plt.ylabel("Stress after von Mises $\\sigma_M$ [MPa]")
plt.xlabel("Rotor speed $\\Omega$ [$\\frac{1}{\\mathrm{min}}$]")
plt.axvline(x=6.0, label='Cut-in rot speed')
plt.axvline(x=12.1, label='Nominal rot speed')
plt.axhline(y=200, label='AL Yield Strenght')
plt.xlim(0, RPMList[-1])
plt.savefig("fig/Stress.svg")
plt.savefig("fig/Stress.png")
plt.show()
#save file.


# Post-processing with pyvista for publication-ready figures
plotGen = True
if plotGen:
    import pyvista as pv
    #import vtk
    #from vtk.numpy_interface import dataset_adapter as dsa

    font = "arial"  #"times"
    colorMap = "jet"  #"winter" #"gnuplot" #"jet" #"bcyr" #"rainbow" #"inferno" #"plasma" #"viridis_r" #"bwr" #
    showMesh =False#True
    FactorStatic = 25
    FactorModal = 25
    for i, Omega in enumerate(RPMList):
        #Static
        Responses = ["DISPLACEMENT", "VON_MISES_STRESS_MIDDLE_SURFACE"]
        FileNames = ["Displacement", "Stress"]
        BarTitles = ["displacement [mm]", "Mises stress [MPa]"]
        inputStaticFile = "StaticResults_"+str(Omega)+os.sep+"Structure_0_1.vtk"
        for j in range(len(Responses)):
            pv.set_plot_theme("document")
            feMesh = pv.UnstructuredGrid(inputStaticFile)
            plotter = pv.Plotter()
            if showMesh:
                plotter.add_mesh(feMesh.copy(), color='k', style='wireframe')
            sargs = dict(height=0.5, vertical=True, position_x=0.95,
                         position_y=0.25, font_family=font,  bold=False,
                         n_colors=100, title_font_size=20, n_labels=11,
                         label_font_size=16)
            plotter.add_mesh(feMesh.warp_by_vector("DISPLACEMENT", FactorStatic),
                             show_scalar_bar=True, cmap=colorMap,
                             culling=True, scalars=Responses[j],
                             stitle=BarTitles[j]+"\n", scalar_bar_args=sargs)
            plotter.view_isometric()
            plotter.save_graphic(folderFig+os.sep+FileNames[j]+"_"+str(Omega)+".svg",
                                 raster=False, painter=True)
            plotter.save_graphic(folderFig+os.sep+FileNames[j]+"_"+str(Omega)+".tex",
                                 raster=False)
            plotter.show(screenshot=folderFig+os.sep+FileNames[j]+"_"+str(Omega)+".png",
                         window_size=[1024, 768])

        #Modal
        vtkFolder = "EigenResults_"+str(Omega)
        inputModalFile = "prestressed_eigen_EigenResults_1_0.vtk"
        BarTitles = ["mode "+str(ii+1) for ii in range(nEig)]
        for j in range(nEig):
            pv.set_plot_theme("document")
            feMesh = pv.UnstructuredGrid(vtkFolder+os.sep+inputModalFile)
            #plotter = pv.Plotter(line_smoothing=True)
            plotter = pv.Plotter()
            if showMesh:
                plotter.add_mesh(feMesh.copy(), color='k', style='wireframe')
            sargs = dict(height=0.5, vertical=True, position_x=0.95,
                         position_y=0.25, font_family=font,  bold=False,
                         n_colors=100, title_font_size=20, n_labels=11,
                         label_font_size=16)
            actor = plotter.add_mesh(feMesh.warp_by_vector(feMesh.array_names[j],
                                                           FactorModal),
                                     show_scalar_bar=True, cmap=colorMap,
                                     culling=True, scalars=feMesh.array_names[j],
                                     stitle=BarTitles[j]+" \n"+ "%.2f" % fnList[i][j] + " Hz \n",
                                     scalar_bar_args=sargs)
            plotter.view_isometric()
            #plotter.save_graphic(folderFig+os.sep+"Eigen"+str(j+1)+"_"+str(Omega)+".eps", raster=False, painter=True)
            #plotter.save_graphic(folderFig+os.sep+"Eigen"+str(j+1)+"_"+str(Omega)+".pdf", raster=False, painter=True)
            plotter.save_graphic(folderFig+os.sep+"Eigen"+str(Omega)+"_"+str(Omega)+".svg",
                                 raster=False, painter=True)
            plotter.save_graphic(folderFig+os.sep+"Eigen"+str(Omega)+"_"+str(Omega)+".tex",
                                 raster=False, painter=True)
            plotter.view_isometric()
            plotter.show(screenshot=folderFig+os.sep+"Eigen"+str(Omega)+"_"+str(Omega)+".png",
                         window_size=[1024, 768])


animate = False
if animate:
    filename = "sphere-shrinking.mp4"
    mesh = pv.Sphere()
    mesh.cell_arrays["data"] = np.random.random(mesh.n_cells)

    plotter = pv.Plotter()
    # Open a movie file
    plotter.open_movie(filename)

    # Add initial mesh
    plotter.add_mesh(mesh, scalars="data", clim=[0, 1])
    # Add outline for shrinking reference
    plotter.add_mesh(mesh.outline_corners())

    print('Orient the view, then press "q" to close window and produce movie')

    # Render and do NOT close
    plotter.show(auto_close=False)

    # Run through each frame
    plotter.write_frame()  # write initial data

    # Update scalars on each frame
    for i in range(100):
        random_points = np.random.random(mesh.points.shape)
        mesh.points = random_points * 0.01 + mesh.points * 0.99
        mesh.points -= mesh.points.mean(0)
        mesh.cell_arrays["data"] = np.random.random(mesh.n_cells)
        plotter.write_frame()  # Write this frame

    # Be sure to close the plotter when finished
    plotter.close()
