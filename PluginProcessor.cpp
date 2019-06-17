/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <vector>
#include <iostream> 
#include <algorithm>

using namespace std;

#define M_PI 3.141592653589793238460
#define TIMEBUFFER 0.04
#define SCOPESIZE 4096
#define NBFRET 12
#define NBSTRING 6
struct Model;
struct result;

//==============================================================================
Jacode_iiiAudioProcessor::Jacode_iiiAudioProcessor() : //class audioprocessor constructor

	//input buffer//
	X(2),

	//sweep//
	trueString(6),
	trueFret(13),

	//classifier parameters//
	mu(2),
	C(2, 2),
	JMatrix(78),

	pitchReference(78),

	w0(6, 13), //= squeeze(est_f0);
	B(6, 13),  //= squeeze(BCoeff);
	P(0),

	prediction(0),

	//onset//
	counter(0),
	thresholdValue(0),
	onsetdetected(false),
	newSegment(1, 0),

	storageActual(1, 0),
	storagePast(1, 0),

	//moi
	afficheValue(0),
	DataLastBuffer(192000, 0),

	//affichage
	scopeDataIsPlayed(),
	scopeDataFret(),
	scopeDataString(),
	counterAnalyser(0),
	nextFFTBlockReady(false),

	//result
	stringPlayed(0),
	fretPlayed(0),
	timeOfOnset(1, 0),

#ifndef JucePlugin_PreferredChannelConfigurations
      AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
	classifier = Calcul_parametre();
	pitchReference.resize(78);
	for (int i = 0; i < 78; ++i)
	{
		pitchReference[i] = classifier.Mu(0, i);
	}
}

Jacode_iiiAudioProcessor::~Jacode_iiiAudioProcessor()                           //class audioprocessor destructor
{
}

//==============================================================================//
//=================================METHODE======================================//
//==============================================================================//

const String Jacode_iiiAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool Jacode_iiiAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool Jacode_iiiAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool Jacode_iiiAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double Jacode_iiiAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int Jacode_iiiAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int Jacode_iiiAudioProcessor::getCurrentProgram()
{
    return 0;
}

void Jacode_iiiAudioProcessor::setCurrentProgram (int index)
{
}

const String Jacode_iiiAudioProcessor::getProgramName (int index)
{
    return {};
}

void Jacode_iiiAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void Jacode_iiiAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
}

void Jacode_iiiAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool Jacode_iiiAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

//==============================================================================//
//===============================------------===================================//
//==============================================================================//


void Jacode_iiiAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{

    ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();
	
    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
    // Make sure to reset the state if your inner loop is processing
    // the samples and the outer loop is handling the channels.
    // Alternatively, you can process the samples with the channels
    // interleaved by keeping the same state.

	//main process loop

    for (int channel = 0; channel < totalNumInputChannels; ++channel)
    {
		if (channel == 0)
		{
			int segmentSize = TIMEBUFFER * getSampleRate(); //40ms en sample
			int frameSize   = 0.001 * getSampleRate();      //1ms en sample

			//onset detector/
			timeOfOnset.clear();

			Onset(buffer.getReadPointer(channel), buffer.getNumSamples(), getSampleRate(), timeOfOnset, counter, storageActual, storagePast, thresholdValue);

			for (int i = 0; i < timeOfOnset.size(); i++)
			{

				if (timeOfOnset[i] == -1) //detection of an Onset
				{
					if (onsetdetected == false)
					{
						onsetdetected = true;
						newSegment.clear();

						for (int j = 0; j < frameSize; ++j)//fill with the first frame
						{
							newSegment.push_back(storagePast[(double(i) * frameSize + j + double(segmentSize))]);
						}

						fretPlayed = 4;
						stringPlayed = 4;

					}
					else if(newSegment.size() < segmentSize)//if a new onset is detect earlier
					{

						afficheValue = PitchEstimator(newSegment);
						fretPlayed = 5;
						stringPlayed = 5;

						newSegment.clear(); //we create a new segment

						for (int j = 0; j < frameSize; ++j)
						{
							newSegment.push_back(storagePast[(double(i) * frameSize + j + double(segmentSize))]);
						}
						//we continue to fill the new segment

					}
					else // no new onset and the segment have its 40ms but this will never arrive normaly cause we will not detect between 40ms buffer
					{
						afficheValue = PitchEstimator(newSegment);
						fretPlayed = 6;
						stringPlayed = 6;
						onsetdetected = false;
					}
				}
				else 
				{
					
					if ((newSegment.size() < segmentSize) && (onsetdetected==true)) // no detection and no 40ms so just fill the segment
					{
						for (int j = 0; j < frameSize; ++j)
						{
							newSegment.push_back(storagePast[(double(i) * frameSize + j + double(segmentSize))]);
						}
						fretPlayed = 7;
						stringPlayed = 7;
						
					}
					else if (onsetdetected==true) //segment ready so predict and clear
					{
						
						afficheValue= PitchEstimator(newSegment);
						fretPlayed = 8;
						stringPlayed = 8;

						onsetdetected = false;
						newSegment.clear();
					}
					else 
					{
						//do nothing
						fretPlayed = 9;
						stringPlayed = 9;
					}
					
				}


				pushNextSampleIntoFifo(timeOfOnset[i], fretPlayed, stringPlayed, (4096));

			}
			//end//

			//featuresExtraction(w0, B, channelData, bufsize);

			/*for (int j = 0; j < trueString; j++)
			{
				for (int i = 0; i < trueFret; i++)
				{
					X(0) = w0(j, i);
					X(1) = B(j, i);

					for (int k = 0; k < (trueFret * trueString); k++)
					{

						mu(0) = classifier.Mu(0, k);
						mu(1) = classifier.Mu(1, k);
						C = classifier.Sigma[k];
						P = classifier.W(k);

						JMatrix(k) = JClassifier(mu, C, P, X);

					}
				}
			}*/

			//int position;
			//JMatrix.maxCoeff(&position);
			//afficheValue = counter;
		}
    }
}

//==============================================================================//
//===============================------------===================================//
//==============================================================================//

bool Jacode_iiiAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* Jacode_iiiAudioProcessor::createEditor()
{
    return new Jacode_iiiAudioProcessorEditor (*this);
}

//==============================================================================
void Jacode_iiiAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void Jacode_iiiAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================//
//=================================MOI==========================================//
//==============================================================================//

double Jacode_iiiAudioProcessor::getAfficheValue()
{
	return(afficheValue);
}

//==============================================================================//
//=================================Analyser=====================================//
//==============================================================================//

void Jacode_iiiAudioProcessor::pushNextSampleIntoFifo(bool isPlayed, int fretPlayed, int stringPlayed, int size) noexcept
{
	if (!nextFFTBlockReady)
	{
		if (DataScope.isPlayed.size() != size)
		{
			DataScope.Fret.clear();
			DataScope.String.clear();
			DataScope.isPlayed.clear();
			DataScopeFifo.Fret.clear();
			DataScopeFifo.String.clear();
			DataScopeFifo.isPlayed.clear();

			for (int i = 0; i < size;i++)
			{
				DataScope.Fret.push_back(0);
				DataScope.String.push_back(0);
				DataScope.isPlayed.push_back(0);
				DataScopeFifo.Fret.push_back(0);
				DataScopeFifo.String.push_back(0);
				DataScopeFifo.isPlayed.push_back(0);
			}
		}

		std::rotate(DataScopeFifo.Fret.begin(),     DataScopeFifo.Fret.begin() + 1,     DataScopeFifo.Fret.end());
		std::rotate(DataScopeFifo.String.begin(),   DataScopeFifo.String.begin() + 1,   DataScopeFifo.String.end());
		std::rotate(DataScopeFifo.isPlayed.begin(), DataScopeFifo.isPlayed.begin() + 1, DataScopeFifo.isPlayed.end());

		DataScopeFifo.isPlayed[size - 1.0] = isPlayed;
		DataScopeFifo.String[size - 1.0]   = stringPlayed;
		DataScopeFifo.Fret[size - 1.0]     = fretPlayed;

		DataScope.Fret     = DataScopeFifo.Fret;
		DataScope.String   = DataScopeFifo.String;
		DataScope.isPlayed = DataScopeFifo.isPlayed;
		nextFFTBlockReady = true;
	}
	else
	{
		std::rotate(DataScopeFifo.Fret.begin(), DataScopeFifo.Fret.begin() + 1, DataScopeFifo.Fret.end());
		std::rotate(DataScopeFifo.String.begin(), DataScopeFifo.String.begin() + 1, DataScopeFifo.String.end());
		std::rotate(DataScopeFifo.isPlayed.begin(), DataScopeFifo.isPlayed.begin() + 1, DataScopeFifo.isPlayed.end());

		DataScopeFifo.isPlayed[size - 1.0] = isPlayed;
		DataScopeFifo.String[size - 1.0]   = stringPlayed;
		DataScopeFifo.Fret[size - 1.0]     = fretPlayed;
	}
}

void Jacode_iiiAudioProcessor::drawNextFrameOfSpectrum()
{
	for (int i = 0; i < SCOPESIZE; i++)
	{
		scopeDataIsPlayed[i] = DataScope.isPlayed[floor(i * DataScope.isPlayed.size() / SCOPESIZE)];
		scopeDataFret[i]     = DataScope.String[floor(i * DataScope.String.size() / SCOPESIZE)];
		scopeDataString[i]   = DataScope.Fret[floor(i * DataScope.Fret.size() / SCOPESIZE)];
	}
}

bool Jacode_iiiAudioProcessor::getNextFFTBlockReady()
{
	return nextFFTBlockReady;
}

void Jacode_iiiAudioProcessor::setNextFFTBlockReady(bool value)
{
	nextFFTBlockReady = value;
}

void Jacode_iiiAudioProcessor::drawFrame(Graphics& g)
{
	for (int i = 1; i < SCOPESIZE; i++)
	{
		auto width = 1000;
		auto height = 600;

		if (scopeDataIsPlayed[i])
		{
			g.drawEllipse((float)jmap(i, 0, SCOPESIZE - 1, 65, width-60), (0.25 * height+63-(5.0- scopeDataString[i])*24), 12, 12, 1);
			g.drawSingleLineText(std::to_string(scopeDataFret[i]), (float)jmap(i, 0, SCOPESIZE - 1 -2, 65, width - 60), (0.25 * height + 63 + 10 - (5.0 - scopeDataString[i]) * 24));
		}
	}	
}

void Jacode_iiiAudioProcessor::setThresholdValue(int value)
{
	thresholdValue = value;
}

//==============================================================================//
//================================Prediction====================================//
//==============================================================================//

float Jacode_iiiAudioProcessor::PitchEstimator(std::vector<float> const& segment)
{
	std::vector<float> processSegment;
	processSegment = segment;
	Preprocessing(processSegment);
	return (HarmonicSummation(processSegment, f0LimitsInf, f0LimitsSup, nbOfHarmonicsInit, getSampleRate()));
};

float Jacode_iiiAudioProcessor::BetaEstimator()
{
	return 0;
};

void Jacode_iiiAudioProcessor::FretStringPrediction(float observedPitch , std::vector<float> pitchReference , int nbFret, Model classifier)
{

	std::vector<int> fretCandidate;
	std::vector<int> stringCandidate;
	PitchCandidate(fretCandidate, stringCandidate, observedPitch, pitchReference, NBSTRING);

};

float Jacode_iiiAudioProcessor::PluckPositionPrediction()
{
	return 0;
};

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Jacode_iiiAudioProcessor();
}

















//thread t1(&Jacode_iiiAudioProcessor::Prediction, std::ref(newSegment), std::ref(afficheValue));//, std::ref(f0LimitsInf), std::ref(f0LimitsSup), std::ref(nbOfHarmonicsInit), getSampleRate());
//auto myFuture = std::async(std::launch::async, Predict, std::ref(newSegment), std::ref(afficheValue));

