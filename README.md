# Guitar-transcription-VST

First git-hub version of the waspaa paper implementation in c++.


What I have to do !

First, It was interesting to compute the classifier but it can’t be used if I don’t have the features. The features for this model are B and Fo which can be determine by some computation but if I have an audio segment to analyse. The audio segment can be collected by an onset detector. The idea is to detect the transient of the sound to get a pluck attack and then keep the sound for around 40 ms. 
This onset detector is done on matlab with a function. I must find a way to do my own in c++.
The vst could test this detector by using the microphone and some clap of hands.

edit : 10/06

C++ done for the onset detector. The plotting is until now not adapted... need to change that !
