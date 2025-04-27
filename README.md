The code for "A Novel Algorithm and Software for Efficient Global Gravimetric Forward Modeling in the Spherical Coordinate System"
A detailed description is provided below:
Operating System: Windows 10/11, MATLAB

Data Description:
The provided text includes CRUST1.0 data, where rmap-bd1 to rmap-bd9 are boundary data, and rmap-ro1 to rmap-ro9 are density data.
These cover bathymetry (bd1, bd2, ro1), ice (bd2, bd3, ro2), sediment (upper: bd3, bd4, ro3; middle: bd4, bd5, ro4; lower: bd5, bd6, ro5), and crust (with the same selection method as sediment). In all datasets, the first column is longitude, the second is latitude, and the third is the corresponding value.

Code Description
(1)"Spherical_Coordinate_Fast_Gravity_Forward_Modeling" is the main program.
(2)In the Gravity Parameters section, three types of gravity parameters are available for selection.
(3)In the Calculation Methods section, two calculation methods are provided:Variable Upper and Lower Boundaries: The values of the upper and lower boundaries are variable. This is a basic forward modeling method based on the boundaries and densities.Constant Upper and Lower Boundaries (proposed in this paper): The model is layered so that the upper and lower boundaries are constant within each layer, allowing lateral density variations. The gravity is calculated layer by layer and then summed.
(4)In the Computation Area module, two options are available:Constructing a regular grid: Typically for global coverage, with a range of longitude -179.5° to 179.5° and latitude -89.5° to 89.5°, at 1° intervals.Higher resolutions, such as 0.5° × 0.5° (longitude -179.75° to 179.75°, latitude -89.75° to 89.75°), or finer grids are also supported by adjusting the range and intervals.Importing point data: Allows calculations based on a set of imported points.
(5)In the Input Parameters module, users must input the maximum spherical harmonic degree and computation height.If the "Constant Upper and Lower Boundaries" method is selected, the number of layers must also be specified.
(6)Clicking the Load Model Data button opens a sub-interface to import boundary and density data.
You can click Load Data multiple times to import multiple boundary and density datasets for cumulative calculations.
(7)After clicking the Calculate button, the computed gravity values and statistical results will be displayed in the output area on the right.Clicking the Image button will display the data visualization.
