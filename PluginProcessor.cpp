/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Estimation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Features.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Onset.h"
#include <vector>

#define M_PI 3.141592653589793238460
struct Model;

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

	w0(6, 13), //= squeeze(est_f0);
	B(6, 13),  //= squeeze(BCoeff);
	P(0),

	prediction(0),

	//onset//
	counter(0),
	time(10,0),

	storageActual(10, 0),
	storagePast(10, 0),

	//moi
	afficheValue(0),
	DataLastBuffer(192000,0),

	//affichage
	DataScopefifo(10, 0),
	DataScope(10, 0),
	scopeData(),
	scopeData2(),
	counterAnalyser(0),
	nextFFTBlockReady(false),


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

	//classifier=(Calcul_parametre())

	//onset//
	

	

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
			int bufsize = buffer.getNumSamples();
			int sampleRate = getSampleRate();
			const float* BufferData = buffer.getReadPointer(channel);
			
			//onset detector/
			time.clear();
			Onset(BufferData, bufsize, sampleRate, time, counter, storageActual, storagePast);

			for (int i = 0; i < time.size(); i++)
			{
				pushNextSampleIntoFifo(time[i], (2048));
					/*if (time[i]==-1)
					{
						pushNextSampleIntoFifo(time[i], sampleRate*0.40);
					}
					else
					{
						pushNextSampleIntoFifo(time[i], sampleRate*0.40);
					}*/
					//pushNextSampleIntoFifo(time[i], time.size());
			}
		


			for (int i = 0; i < SCOPESIZE;i++)
			{
				scopeData[i] = 0;//channelData[int(round(i * bufsize / double(scopeSize)))];
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


//=================================MOI===========================================

double Jacode_iiiAudioProcessor::getAfficheValue()
{
	return(afficheValue);
}


//==============================================================================//
//=================================Analyser=====================================//
//==============================================================================//


void Jacode_iiiAudioProcessor::pushNextSampleIntoFifo(float sample, int size) noexcept
{
	if (!nextFFTBlockReady)
	{
		if (DataScope.size() != size)
		{
			DataScope.clear();
			DataScopefifo.clear();
			for (int i = 0; i < size;i++)
			{
				DataScope.push_back(0);
				DataScopefifo.push_back(0);
			}
		}

		std::rotate(DataScopefifo.begin(), DataScopefifo.begin() + 1, DataScopefifo.end());
		DataScopefifo[size - 1.0] = sample;

		DataScope=DataScopefifo;

		nextFFTBlockReady = true;
	}
	else
	{
		std::rotate(DataScopefifo.begin(), DataScopefifo.begin() + 1, DataScopefifo.end());
		DataScopefifo[size - 1.0] = sample;
	}
}

void Jacode_iiiAudioProcessor::drawNextFrameOfSpectrum()
{
	for (int i = 0; i < SCOPESIZE; i++)
	{

		float level = DataScope[floor(i * DataScope.size() / SCOPESIZE)];

		//float(cos(2*M_PI* afficheValue++ /1000));
		//afficheValue = (level);
		scopeData2[i] = level;//jlimit(-1.0f, +1.0f, level);
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
		auto width = 1200;
		auto height = 1000;
		g.drawLine({ (float)jmap(i - 1, 0, SCOPESIZE - 1, 0, width),
							  jmap(scopeData[i - 1], -1.0f, 1.0f, (float)height, 0.0f),
					 (float)jmap(i,    0, SCOPESIZE - 1, 0, width),
							  jmap(scopeData[i],     -1.0f, 1.0f, (float)height, 0.0f) });

		g.drawLine({ (float)jmap(i - 1, 0, SCOPESIZE - 1, 0, width),
					          jmap(scopeData2[i - 1], -1.0f, 1.0f, (float)height, 0.0f),
			         (float)jmap(i,     0, SCOPESIZE - 1, 0, width),
					          jmap(scopeData2[i],     -1.0f, 1.0f, (float)height, 0.0f) });
	}	
}



//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Jacode_iiiAudioProcessor();
}
