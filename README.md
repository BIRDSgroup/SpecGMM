# SpecGMM

This repository contains codes for the analyses done for the manuscript “SpecGMM: .....” by Saish Jaiswal, Hema A Murthy, and Manikandan Narayanan.

## License preamble 

Copyright 2024 BIRDS Group, IIT Madras

SpecGMM is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

SpecGMM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with Demixer.  If not, see <https://www.gnu.org/licenses/>.

## Setup Instructions    

Note: Setup instructions are to help the user set up the SpecGMM project on a local machine.  

1) Place the downloaded folder SpecGMM in the MATLAB working directory and add it to the path (right-click on the SpecGMM folder and choose the option "Add to Path -> Selected Folders and SubFolders").  
2) Change MATLAB's current working directory to ~/SpecGMM/.

## Requirements
- MATLAB R2020a or above
- [Fathom Toolbox](https://www.usf.edu/marine-science/research/matlab-resources/fathom-toolbox-for-matlab.aspx)
- [MATLAB Progress Bar](https://github.com/JAAdrian/MatlabProgressBar)

Note: MATLAB R2023a on a Linux machine was used to run the experiments. Download the Fathom Toolbox and the MATLAB Progress Bar and put them in the SpecGMM folder and add them to the path.

## Instructions

- Modified baseline code
  
  Run the Baseline_Analysis.m script on the MATLAB terminal:
  ```
  run('Baseline_Analysis.m')
  ```

- SpecGMM Code
  
  Run the SpecGMM_Analysis.m script on the MATLAB terminal:
  ```
  run('SpecGMM_Analysis.m')
  ```

## Folder Description

- [Database](https://github.com/BIRDSgroup/SpecGMM/tree/main/DataBase) contains the datasets used in our analyses
- [Figures](https://github.com/BIRDSgroup/SpecGMM/tree/main/Figures) contains scripts to reproduce the figures in the manuscript
- [Supplementary Data Files](https://github.com/BIRDSgroup/SpecGMM/tree/main/Supplementary%20Data%20Files) contains supplementary files containing detailed analyses and results
