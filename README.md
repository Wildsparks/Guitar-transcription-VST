# Guitar-transcription-VST

First git-hub version of the waspaa paper implementation in c++.

keep in mind : some variables are choose a bit randomly :

_betares => une seul valeur possible pour l'instant

_fftorder for onset and pitch

_number of realisation for the model

_scope size

_number of string

_number of fret

_standar deviation of the noise

_size of the windows for onset

_threashold for onset

_size of mean windows and max windows for onset

_résolution for pitch detection

_number of zoom for pitch detection

_temperament considération

_range of frequence for pitch detection

_time between two prediction

What I have to do !

First, It was interesting to compute the classifier but it can’t be used if I don’t have the features. The features for this model are B and Fo which can be determine by some computation but if I have an audio segment to analyse. The audio segment can be collected by an onset detector. The idea is to detect the transient of the sound to get a pluck attack and then keep the sound for around 40 ms. 
This onset detector is done on matlab with a function. I must find a way to do my own in c++.
The vst could test this detector by using the microphone and some clap of hands.

edit : 10/06

C++ done for the onset detector. The plotting is until now not adapted... need to change that !

edit : 10/06 : 16.47

plotting perfect with the onset detector working in real time over guitar sound.

edit : 17/06 : 16.37

pitchcandidate implemented + graphical adaptation over string and fret estimation.
