from paraview.simple import *
xHvtk = LegacyVTKReader(FileNames=['XHP_saltwatereos-enhanced_computations.vtk'])
renderView1 = GetActiveViewOrCreate('RenderView')
xHvtkDisplay = Show(xHvtk, renderView1)
xHvtkDisplay.Representation = 'Surface'
renderView1.AxesGrid.Visibility = 1
xHvtkDisplay.Scale = [1, 8.62069e-05, 0.000555556]
renderView1.AxesGrid.DataScale = [1, 8.62069e-05, 0.000555556]
renderView1.AxesGrid.DataBoundsInflateFactor = 0
renderView1.AxesGrid.XTitle = 'Salinity'
renderView1.AxesGrid.YTitle = 'Enthalpy (kJ/kg)'
renderView1.AxesGrid.ZTitle = 'Pressure (bar)'
renderView1.AxesGrid.XTitleFontSize = 16
renderView1.AxesGrid.XTitleBold = 1
renderView1.AxesGrid.YTitleFontSize = 16
renderView1.AxesGrid.YTitleBold = 1
renderView1.AxesGrid.ZTitleFontSize = 16
renderView1.AxesGrid.ZTitleBold = 1
#set default data source as PhaseRegion
paraview.simple._DisableFirstRenderCameraReset()
legacyVTKReader1 = GetActiveSource()
renderView1 = GetActiveViewOrCreate('RenderView')
legacyVTKReader1Display = GetDisplayProperties(legacyVTKReader1, view=renderView1)
ColorBy(legacyVTKReader1Display, ('POINTS', 'PhaseRegion'))
legacyVTKReader1Display.RescaleTransferFunctionToDataRange(True, False)
legacyVTKReader1Display.SetScalarBarVisibility(renderView1, True)
phaseRegionLUT = GetColorTransferFunction('PhaseRegion')

renderView1.ResetCamera()
