# cc_sandy_matlab
 
This repository contains the code to compute much of the analysis presented in "Economic Damages from Hurricane Sandy Attributable to Sea Level Rise Caused by Anthropogenic Climate Change." It expects data from <link>, which includes simulation water height fields, in a subdirectory named "data."

Importantly, this repository includes the code to load and correct the ADCIRC simulations. To do this correction, run the scripts in order:

WarpSandy()
SandyCorrectionError()

This will load the table "sandyObservationData," and build the "sandyWarpedModelCompleteData" matrix variable, which contains latitude and longitude in the first two columns, and the corrected max water heights for simulations alt01-alt09 (alt06 being historical) in the remaining columns.

Errors can be assessed using the script:

EvaluateSandyError( sandyObservationData, sandyWarpedModelCompleteData, 9 )

This does not include data nor code to compute exposure of any variable, as these are part of Climate Central's proprietary analysis systems, which due to licensing restrictions, cannot be shared.
