This is a repository of a version of the ray trace simulation of a compact fourier transform spectrometer (FTS) written by Mira Liu over the years 2016-2019 while working with Dr. Stephan S. Meyer, Dr. Ritoban Basu Thakur, and Zhaodi Pan. All of this was written in python, and tested with jupyter notebooks. 

Included in this repository are three base .py files necessary to run the simulation. Also included are sample .ipynb jupyter notebook files to demonstrate the use of different functions. Within each base .py file there are also descriptions of the functions. There is also a labeled diagram of the FTS called 'CAD_FTS_labeled' that demonstrates all physical aspects included in this simulation. DetectorPlots.png is an example of a plot of the hit points of the rays on the detector. It only includes 8 mirror positions between 0 and the furthest +y location of the mirror. ALL hit points on the detector at that mirror position are colored respectively (red for mirror position 0, orange for +y/8, yellow for +y/4... black for +y).

1. CheckGeometric.ipynb is a jupyter notebook that statistically simulates the modulating envelope caused by loss of rays due to geometric scattering (simply fraction of launched rays that reach the detector). 
2. PlotSimulation.ipynb is a jupyter notebook that plots the paths of rays (n=50) through simulation (specified by a function that sets which of the 8 paths the rays will be taking). 
3. RunSimulation.ipynb is a jupyter notebook that calls functions and gives appropriate input to launch rays through the simulation and can implement different methods of power summation at the detector.
More notes are included in the notebooks.
3. TimeConstantFit is a jupyter notebook that explicitly goes through the analysis of the time constant (excluding the MCMC that was run that takes a very long time). The parameters returned by the MCMC are plugged into the mathematical model of a convolved interferogram, which is compared to a 30mm/s experimental interferogram. Those parameters are then used to deconvolve the 30mm/s interferogram, which is then compared to a 10mm/s experimental interferogram. 
4. GitAnalysisOverview is a jupyter notebook that is simply markdown copy of an old analysis explanation. 

Previous Runs contains the pickle files of both experimental data and simulated data.
20160921_1724_20160921_10mms_new_fts_all_band.pkl and 20170615_1729_30mms_reference.pkl were used in TimeConstantFit.

20180318_2211_150GHz_new_polarizer_35.pkl is the reference experimental data set run with the 144GHz source. It can be ploted to compare with the results of the simulation with Power1 = (d['sig0F']), Power = [x/Power1.max() for x in Power1], Delay = (d['delay0F']). 

There are pickle files of recently run simulations saved in the folder 'Previous Runs'. There is a note in the folder describing the names of the runs. There are many more pickle files of simulations of different versions of code available if wanted. 

Lastly there is a .pkl file that contains every step of the simulations of a 90Ghz light source with 500 rays, located at this link https://uchicago.box.com/s/vhxjh8tygffc7ips137com2ajwc9lrfs. ToPickle.txt included in this repository describes the format of the pickle file. 

All code is written by Mira Liu.
- 10/08/2019
