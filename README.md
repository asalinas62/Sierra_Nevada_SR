# Sierra_Nevada_SR

The programs are the following: 
1.	‘Reading_Fourier.py’. It reads the raw data in the time domain and obtains the amplitude spectrum in the calibrated band of the measurement system.
2.	‘Anthropo.ipynb’. When present, this program eliminates the anthropogenic noise in the amplitude spectrum of the previous step.
3.	‘Lorentz_Lmfit.py’. The program applies a non-linear fit using the ‘Lmfit’ module of Python to generate a Lorentzian fitted spectrum together with its corresponding error parameters.
4.	‘Packaging.py’. This program carries out a final storage of the relevant data. The program generates a file for each month and sensor that allows a direct reading and recovering of the structure (arrays of the Numpy package of Python) of all the processed data generated by the previous packages for subsequent analysis purposes.

The folder hierarchy of the data provided is as follows:

S_N_Data. This folder contains subfolders for all the years with time domain raw data from the station, named with the four digits corresponding to the specific year, e.g., 2014. Each year contains folders for each month, named with the four digits corresponding to year and month, e.g., 1412 stands for the raw data during December 2014. Within each month folder, the data and information files for each sensor and each hour are available. Data are stored in files containing measurements corresponding to a time period of approximately one hour. The filenames begin with a common part, “smplGRTU1_sensor_”, followed by a specific part to denote the sensor used (0 for the NS orientation and 1 for the EW orientation), the date and the initial time of the measurement recorded.  For example, smplGRTU1_sensor_0_1412010430, stands for data measured by the NS-oriented magnetometer, in the year 2014, month 12, day 01, hour 04, and starting minute, 30. The information file has the same name but ends with _info.txt. Since the initial time is not fixed, the filenames in each folder are not completely determined. For this reason, each month folder includes two files, 'ficheros0' and 'ficheros1', that contain the set of filenames for that month. Therefore, for each normally measured we have 24 data files and also 24 information files for each sensor. Each hour data file occupies 1.8 MB, so each month has roughly a data volume of 2.6 GB.

S_N_DF. This folder contains the intermediate results of the data processing, i.e., the outputs of the programs for stages 1 to 3 listed in the previous section, but not for the last stage, devoted to packaging the final information. It has the same structure of years and months as the S_N_Data folder. Within each month, output files for each step of the processing can be found. Note that, if the code is executed, the folder structure of years and months has to be previously created. 
The files associated to the programs used for each intermediate processing stage are:
1.	‘Reading_Fourier.py’: the amplitude spectra and related information generated by this program are stored in the following files (December 2014 is used as an example): 
•	SR1412_media_0,1: it contains the amplitude spectra within the calibrated band of the magnetometers (6 - 25 Hz), for the NS sensor (sensor 0) and for the EW sensor (sensor 1), respectively.
•	SR1412_satper_0,1: it stores the percentage of those 10 s intervals that show saturation within each ten-minute interval. There is a file for each sensor. 
•	SR1412_saturados_0,1: it contains the different times at which the measurements are saturated with respect to the +-10 V reference. The objective is to keep a record with the moments at which saturations occur for a possible subsequent study of correlation with the occurrence of high intensity lightning events. There is a file for each sensor.
2.	‘Anthropo.ipynb’: the anthropogenic noise detected in the measurements is eliminated from the spectrum. The following file is generated: 
•	SR1412_mediaNA_0,1: similar to SR1412_media_0,1 but with the anthropogenic noise filtered from the amplitude spectra. 
3.	‘Lorentz_Lmfit.py’: a Lorentzian fitting of the amplitude spectrum is made. The output files for this step are: 
•	SR1412_mediaLO_0,1: it contains the parameters of the Lorentzian adjustment of the spectra in amplitude for each ten-minute interval and sensor. 
•	SR1412_mediaLOC_0,1: it includes the fitting error parameter of the function  in each spectrum. 
•	SR1412_mediaLOE_0,1: it contains the errors that the fitting program assigns to each calculated fitting parameter.

S_N_Proc. This folder stores the final SR results of the raw data processing. There is a file for each month of measurements and for each sensor. The filename makes reference to the year, month and sensor used. The name for December 2014 and sensor 0 is SN_1412_0.npz, for instance. The format is Numpy’s .npz. Each file has a size of 43.3 MB, therefore, all the months of measurements can be stored in the same directory in order to facilitate the reading for a later analysis.
