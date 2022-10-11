# Brain Fluencerate Simulation
Simulate the fluence rate and refelctance, pathlength of the 3D head model using MCX.

---

## Prepare
1. Install MCXLAB on the computer  
Make sure the simulation can run properly.
2. Prepare the subject MRI 3D model file  
Containing the voxel, EEG points and head model surface mesh.  Some example models can be downloaded [here](https://drive.google.com/drive/folders/1UCdTMV5SDQeSsa4-zapeVYGJhnCzXKyx?usp=sharing).
3. Prepare the subject's tissue parameter  
Containing the A,K and absorber concentration for scalp, skull and gray matter.  Preparing the absorption and scattering spectrum is also Okey.

---

## Usage
1. Run `Find_probe_pos_before_sim.m` to prepare the location and direction of the sources and probes.  The result will be store in `models` folder.
2. Make sure there are some sets of OPs to run in the `literature_OPs` folder.
3. Run the `main_simulation_literature_param.m` to simulate the fluence rate.
4. Run `main_cal_layer_eng_hist.m` to calculate the layer absorbed energy and it's histogram.


---

## Parameters
The parameters that can be setting.  All of them are in the `main_simulation.m`  
* Subjects:
    * subject_name_arr: The names of subjects to run simulation.  Which also need their segmented MRI voxel model in `models` directory.  
* Detectors:
    * num_SDS: The number of detectors.
    * SDS_x_arr: The SDS x displacement (cm) compares to source position.
    * SDS_z_arr: The SDS z displacement (cm) compares to source position.
    * reference_index_arr: The index of EEG points that close to the setting detectors.
    * detector_r: The radius of the detector fibers (mm).
    * detector_larger_r: The larger radius for better detect photon in MCX simulation (mm).
    * detector_NA: The numerical aperture of the detector fibers.
    * source_NA: The numerical aperture of the source fiber.
    * fiber_n: The refractive index of the fiber.
* Tissue:
    * n: The refractive index of tissue.
    * outer_n: The refractive index of outer air or porbe surface.
* simulation:
    * output_fluence_rate: Output the total fluence rate
    * output_pathlength: Output the pathlength of detected photons, also the reflectance
    * output_jacobian: Output the jacobian of each detector fiber.
    * num_photon: Number of photons to simulate.
    * output_folder: Output folder name
    * to_plot_figure: Plot the probe position and direction
    * lambda: The wavelength to simulate
    * do_probe_surface: Add a layer of probe surface outside the head model.
    * do_crop: Experimental function, crop the model for faster simulation.
    * crop_range: The range of model to crop.

---
## Output Files
In each subject's output directory, you may find the following files for each wavelength:  
(The * can be replaced by integers, that's the index of the simulated wavelength)  
* sim_summary_*.json: The summary of simulation of this wavelength.  
* average_fluence_*.mat: The average fluence rate for each voxel of this wavelength.  
* PL_*.mat: The pathlength of each detected photon of this wavelength.  


---

## Files
### For Simulation

#### Find_probe_pos_before_sim.m
Find the location of source and detectors, also the source location for laser/LED stimulation.  Run this script before run the main simulation script.

#### main_simulation.m
The main program that preparing model, OPs, settings to simulate.  Run this program to begin simulation. 

#### main_simulation_literature_param.m
The main function to use while you had prepared many sets of OPs to run.  

#### main_cal_layer_eng_hist.m
Calculate the layer absorbed energy and histogram of the simulated results.  

#### fun_MCX_sim_dist2axis.m
The kernel function for the simulation, called `by main_simulation.m`.

#### fun_cal_probe_position_and_direction_2D.m
Calculate the probe position and direction, called by `main_simulation.m`.

#### fun_assign_param.m
Assigh the subject's tissue parameter to the parameter data struct, called by `main_simulation.m`.

#### fun_generate_MCX_tissue_param.m
Generate the optical parameter given the tissue parameter data struct, called by `main_simulation.m`.

#### fun_param_to_mu_spec.m
Generate the optical parameter.  Called by `fun_generate_MCX_tissue_param.m`.

#### fun_add_probe_surface_to_vol.m
Add the probe surface to outside the model.  Called by `main_simulation.m`.

#### fun_cal_layer_energy_and_hist.m
A function calculate the layer absorbed energy (or fluence rate) and histogram.  Called by `main_cal_layer_eng_hist.m`.  

#### GPU_setting.txt
If your computer has more than 1 GPU, then edit this file can set which GPU(s) to use, and the workload between them.  
In matlab, Type `gpuDeviceCount` to check how many GPU do you have.  Type `gpuDevice(1)`, `gpuDevice(2)` and so on to check the index of each GPU.  
The workload for each GPU can be found by run the official mcx gpu benchmark.  

#### stop_flag.txt
If this file exist, and the content in it is `1`, then the simulation will stop after finish this subsim.

---

### After Simulation
#### calculate_layer_fluence.m
After the simulation, calculate the fluence rate for each layer.

#### plot_fluence_rate_result_single_sim.m
After the simulation, plot the fluence rate.

#### plot_fluence_rate_result_compare_wl.m
After the simulation, plot the fluence rate for all wavelength for comparision.

---

## Simulation Steps
1. Open the `main_simulation.m`
2. Setting the subjece name, proper model and parameter files.  Also set the probe settings.
3. Set the `to_plot_figure`(line 33) to 1, and pause at around line 56.  Check the 2D and 3D plot of the subject's model and the probes.  Make sure the probes are on the subject's head surface.  
![head model plot](https://i.imgur.com/PLQJPVB.png)
If the position or direction of the probes are wrong, try to fix the problem.
4. Start to run the MCX simulation by continuing to line 86.
5. After the simulation, use `plot_fluence_rate_result.m`, `plot_fluence_rate_result_compare_wl.m` or `calculate_layer_fluence.m` to check the result.

---

## Update
### 2020/03/16
* Add the voxel size option. The model containing voxel size can be download [here](https://drive.google.com/drive/folders/1UCdTMV5SDQeSsa4-zapeVYGJhnCzXKyx?usp=sharing).  
* Add checkpoint.  
* Set the fluence rate to the right scale.  
* Improve performance.  

### Version 2.01
* 2020/03/27
* Add the version file and checker.
* Set the right default parameter for simulate laser or LED stimulation.
* Improve the plot function.

### Version 3.11
* 2020/04/01
* Fix the setting bug of cone source simulation
    * No longer do 2 times of asin computation of the source angle

### Version 3.21
* 2020/04/04
* Fix the setting bug of cone source simulation
    * No longer gradually move the cone source away while simulate multiple time.

### Version 3.31
* 2020/05/04
* Using a function to calculate the position of cone source.
* Make calculate layer fluence rate as a function.
* Change the output of MCX from fluence rate to energy.

### Version 3.41
* 2020/05/06
* Updata the pre-calculate pos function
    * Use a fitting function to make sure the shined area of cone source is similar to disk source
    * Plot the disk source and cone source on the 3D figure.
* Change the angle of cone source from 60 degree to 40 degree
* Let the calculate probe position and direction prior to simulation. (remove calculate position, plot figures from main simulation script)

### Version 4.00
* 2020/05/16
* Update the simulation code for version 2 MRI models.
* Update the find position code, to make the cone source more close to disk source.

### Version 4.31
* 2020/06/13
* Add the check detpt CV code
* Slightly change the MCX function for load the detected photon seed if there are no photon detected in the first ckeckpoint


---

## TODO
1. Update the normal version (non-litOP) main function.