/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

using namespace std;

#define M_PI       3.141592653589793238460
#define TIMEBUFFER 0.04
#define TIMEFRAME  0.001
#define SCOPESIZE  4096
#define NBFRET     12
#define NBSTRING   6

struct Model;
struct result;

//==============================================================================////==============================================================================//
Jacode_iiiAudioProcessor::Jacode_iiiAudioProcessor() : //class audioprocessor constructor

	//-------first--step-------//

	//--------//onset//--------//
	
	counterOnset  (0),
	thresholdValue(0),

	storageSparse(120, 0),
	storageActual(1, 0),
	storagePast  (1, 0),
	newSegment   (1, 0),
	onsetScope   (1, 0),

	detectionDone(false),
	onsetdetected(false),

	timeOfOnset(1, 0),

	//-------second-step-------//

	//--------//pitch//--------//

	f0LimitsInf (75),
	f0LimitsSup (700),         // boundaries for f0 search grid in Hz.
	nbOfHarmonicsInit(5),      // number of harmonics for initial harmonic estimate(B = 0).








	//-------final-step-------//

	//--------//plot//--------//

	DataScope(),
	DataScopeFifo(),
	nextFFTBlockReady(false),








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
	
	





	//moi
	afficheValue(10),



	

	//result
	stringPlayed(0),
	fretPlayed(0),



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

	//==============================================================================//

	classifier = Calcul_parametre(featureMatrix);

	pitchReference.resize(78);

	for (int i = 0; i < 78; ++i)
	{
		pitchReference[i] = classifier.Mu(0, i);
	}

	for (int i = 0; i < SCOPESIZE; ++i)
	{
		scopeDataIsPlayed[i] = false;
		scopeDataFret[i] = 0;
		scopeDataString[i] = 0;
		scopeDataOnset[i] = 0.0f;//just for working
	}

	//==============================================================================//

}

Jacode_iiiAudioProcessor::~Jacode_iiiAudioProcessor()                           //class audioprocessor destructor
{

}

//==============================================================================////==============================================================================//
//=================================METHODE======================================////=================================METHODE======================================//
//==============================================================================////==============================================================================//

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

			//===================================//
			//=====say to Onset Thread to go=====//
			//===================================//

			int buffsize(buffer.getNumSamples());
			int counterOfBuffer(1);
			std::vector<float> data;
			const float* pointerData(buffer.getReadPointer(channel));

			while (buffsize > TIMEBUFFER * getSampleRate()) //if buffer > 40ms we are not able to predict for each time we need.
			{
				buffsize /= 2;
				counterOfBuffer++;
			}

			data.resize(buffsize);

			for (int i = 0; i < counterOfBuffer; ++i)
			{
				for (int j = 0; j < buffsize; ++j)
				{
					data[j] = pointerData[j + i * buffsize];
				}
				ProcessWithSpectrogramOnsetDetector(&data[0], buffsize, getSampleRate());
			}

			//===================================//
			//===================================//
			//===================================//

		}
		
    }//loop chanel 
}//function end

//==============================================================================////==============================================================================//
//=================================TEST=========================================////==============================================================================//
//==============================================================================////==============================================================================//

float Jacode_iiiAudioProcessor::getAfficheValue()
{
	return(afficheValue);
}

//==============================================================================//
//================================--Process--===================================//
//==============================================================================//

void Jacode_iiiAudioProcessor::ProcessWithSpectrogramOnsetDetector(const float* data, int bufSize, int sampleRate)
{
	int    segmentSize (TIMEBUFFER * sampleRate); //40ms in sample
	int    frameSize   (TIMEFRAME  * sampleRate); //1ms in sample
	int    Noteprediction  (0);
	float time        (-0.040 * 2.0);            //decay time because of storing

	Onset(data, bufSize, sampleRate, timeOfOnset, counterOnset, storageActual, storagePast, thresholdValue, onsetScope, detectionDone, storageSparse);
	
	if (detectionDone)//many buffer are just store => predictions aren't done every time we get buffer but every time we get 40ms.
	{

		int index = 40; //decay for take temporal sample from the time onset is detect

		for (int i = 0; i < timeOfOnset.size(); ++i)
		{
			time += TIMEFRAME; 
			if (timeOfOnset[i] == true) //detection of an Onset
			{

				if (onsetdetected == false)
				{
					onsetdetected = true;
					newSegment.clear();

					for (int j = 0; j < frameSize; ++j)//fill with the first frame
					{
						newSegment.push_back(storagePast[(((int)i + (int)index) * (int)frameSize + (int)j)]);
					}

				}
				else if (newSegment.size() < segmentSize)//if a new onset is detect earlier
				{

					Noteprediction = FretStringPrediction(newSegment);
					fretPlayed = (unsigned int)floor(Noteprediction / NBSTRING);
					stringPlayed = Noteprediction % NBSTRING;

					newSegment.clear(); //we create a new segment

					for (int j = 0; j < frameSize; ++j)//we continue to fill the new segment
					{
						newSegment.push_back(storagePast[(((int)i + (int)index) * (int)frameSize + (int)j)]);
					}

					

				}
				else // no new onset and the segment have its 40ms but this will never arrive normaly cause we will not detect between 40ms buffers
				{
					Noteprediction = FretStringPrediction(newSegment);
					fretPlayed = (unsigned int)floor(Noteprediction / NBSTRING);
					stringPlayed = Noteprediction % NBSTRING;

					newSegment.clear();
					onsetdetected = false;
				}
			}
			else
			{

				if ((newSegment.size() < segmentSize) && (onsetdetected == true)) // no detection and no 40ms : just fill the segment
				{
					for (int j = 0; j < frameSize; ++j)
					{
						newSegment.push_back(storagePast[(((int)i+(int)index) * (int)frameSize + (int)j)]);
					}

				}
				else if (onsetdetected == true) //segment ready : predict and clear
				{

					Noteprediction = FretStringPrediction(newSegment);
					stringPlayed = Noteprediction % NBSTRING;
					fretPlayed = (int)floor(Noteprediction / NBSTRING);

					onsetdetected = false;
					newSegment.clear();

				}
				else
				{
					//do nothing
				}

			}

			if (isnan(onsetScope[i])) //to prevent error
			{
				pushNextSampleIntoFifo(timeOfOnset[i], fretPlayed, stringPlayed, (4096), 0);
			}
			else
			{
				pushNextSampleIntoFifo(timeOfOnset[i], fretPlayed, stringPlayed, (4096), onsetScope[i]);
			}

		}  //----------------//end process//----------------//

	}      //----------------//endif onset//----------------//
	
}

//==============================================================================//
//=================================Analyser=====================================//
//==============================================================================//

void Jacode_iiiAudioProcessor::pushNextSampleIntoFifo(bool isPlayed, int fretPlayed, int stringPlayed, int size, float onsetValue) noexcept
{
	if (!nextFFTBlockReady)
	{
		if (DataScope.isPlayed.size() != size)
		{
			DataScope.Fret.clear();
			DataScope.String.clear();
			DataScope.isPlayed.clear();
			DataScope.Onset.clear();
			DataScopeFifo.Fret.clear();
			DataScopeFifo.String.clear();
			DataScopeFifo.isPlayed.clear();
			DataScopeFifo.Onset.clear();

			for (int i = 0; i < size;++i)
			{
				DataScope.Fret.push_back(0);
				DataScope.String.push_back(0);
				DataScope.isPlayed.push_back(false);
				DataScope.Onset.push_back(0.0f);
				DataScopeFifo.Fret.push_back(0);
				DataScopeFifo.String.push_back(0);
				DataScopeFifo.isPlayed.push_back(false);
				DataScopeFifo.Onset.push_back(0.0f);
			}
		}

		std::rotate(DataScopeFifo.Fret.begin(),     DataScopeFifo.Fret.begin() + 1,     DataScopeFifo.Fret.end());
		std::rotate(DataScopeFifo.String.begin(),   DataScopeFifo.String.begin() + 1,   DataScopeFifo.String.end());
		std::rotate(DataScopeFifo.isPlayed.begin(), DataScopeFifo.isPlayed.begin() + 1, DataScopeFifo.isPlayed.end());
		std::rotate(DataScopeFifo.Onset.begin(),    DataScopeFifo.Onset.begin() + 1,    DataScopeFifo.Onset.end());

		DataScopeFifo.isPlayed[(int)size - (int)1] = isPlayed;
		DataScopeFifo.String  [(int)size - (int)1] = stringPlayed;
		DataScopeFifo.Fret    [(int)size - (int)1] = fretPlayed;
		DataScopeFifo.Onset   [(int)size - (int)1] = onsetValue;

		DataScope.Fret     = DataScopeFifo.Fret;
		DataScope.String   = DataScopeFifo.String;
		DataScope.isPlayed = DataScopeFifo.isPlayed;
		DataScope.Onset    = DataScopeFifo.Onset;
		nextFFTBlockReady  = true;
	}
	else
	{
		std::rotate(DataScopeFifo.Fret.begin(),     DataScopeFifo.Fret.begin() + 1,     DataScopeFifo.Fret.end());
		std::rotate(DataScopeFifo.String.begin(),   DataScopeFifo.String.begin() + 1,   DataScopeFifo.String.end());
		std::rotate(DataScopeFifo.isPlayed.begin(), DataScopeFifo.isPlayed.begin() + 1, DataScopeFifo.isPlayed.end());
		std::rotate(DataScopeFifo.Onset.begin(),    DataScopeFifo.Onset.begin() + 1,    DataScopeFifo.Onset.end());

		DataScopeFifo.isPlayed[(int)size - (int)1] = isPlayed;
		DataScopeFifo.String  [(int)size - (int)1] = stringPlayed;
		DataScopeFifo.Fret    [(int)size - (int)1] = fretPlayed;
		DataScopeFifo.Onset   [(int)size - (int)1] = onsetValue;
	}
}

void Jacode_iiiAudioProcessor::drawNextFrameOfSpectrum()
{
	for (int i = 0; i < SCOPESIZE; ++i)
	{
		scopeDataIsPlayed[i]   = DataScope.isPlayed[(int)floor(i * DataScope.isPlayed.size() / SCOPESIZE)];
		scopeDataString[i]     = DataScope.String  [(int)floor(i * DataScope.String.size()   / SCOPESIZE)];
		scopeDataFret[i]       = DataScope.Fret    [(int)floor(i * DataScope.Fret.size()     / SCOPESIZE)];
		scopeDataOnset[i]      = DataScope.Onset   [(int)floor(i * DataScope.Onset.size()    / SCOPESIZE)];
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
	for (int i = 1; i < SCOPESIZE; ++i)
	{
		int width = 1600;
		int height = 600;

		if (scopeDataIsPlayed[i])
		{
			g.drawEllipse((float)jmap(i, 0, SCOPESIZE - 1, 90, width-80), (float)(0.25 * height + 58.0 - (5.0 - scopeDataString[i]) * 24.0), 15, 15, 1);
			g.drawSingleLineText(std::to_string(scopeDataFret[i]), jmap(i, 0, SCOPESIZE - 1 -5, 90, width - 80), (int)(0.25 * height + 58.0 + 10.0 + 2.0 - (5.0 - scopeDataString[i]) * 24.0));
			g.fillRect((int)(250.0+(afficheValue)*10.0), 270, 8,200);
		}
		afficheValue = scopeDataOnset[i - 1];
		g.drawLine({ (float)jmap(i - 1, 0, SCOPESIZE - 1, 720, width-140) , jmap(scopeDataOnset[i - 1], -1.0f, 1.0f, (float)height-103.0f, 265.0f),
				     (float)jmap(i ,    0, SCOPESIZE - 1, 720, width-140) , jmap(scopeDataOnset[i],     -1.0f, 1.0f, (float)height-103.0f, 265.0f) });
		afficheValue = scopeDataOnset[i - 1];
		g.drawLine(720, jmap(float(thresholdValue) / 100, 0.0f, 1.0f, (float)height - 218, 265.0f), (float)(width - 140.0), jmap(float(thresholdValue) / 100, 0.0f, 1.0f, (float)height - 218, 265.0f));


	}	
}

void Jacode_iiiAudioProcessor::setThresholdValue(int value)
{
	thresholdValue = value;
}

//==============================================================================//
//================================Prediction====================================//
//==============================================================================//
int Jacode_iiiAudioProcessor::FretStringPrediction(const std::vector<float>& segment)
{
	float f0Features(0);
	float betaFeatures(0);

	std::vector<float> processSegment;
	std::vector<std::complex<float>> hilberOutput;
	processSegment = segment;

	Preprocessing(processSegment, hilberOutput);

	float pitch = PitchEstimator(processSegment);

	std::vector<int> listOfCandidate = BetaEstimator(processSegment, pitch, pitchReference, f0Features, betaFeatures);

	X(0) = f0Features;
	X(1) = betaFeatures;

	int position(0);

	JMatrix.resize(listOfCandidate.size());
	for (int k = 0; k < listOfCandidate.size(); k++)
	{
		int indexCandidat = listOfCandidate[k];
		mu(0) = classifier.Mu(0, indexCandidat);
		mu(1) = classifier.Mu(1, indexCandidat);
		C = classifier.Sigma[indexCandidat];
		P = classifier.W(indexCandidat);

		JMatrix(k) = JClassifier(mu, C, P, X);

	}
	JMatrix.maxCoeff(&position);

	afficheValue = float(PluckPositionPrediction(hilberOutput, f0Features, betaFeatures, (int)floor(listOfCandidate[position] / NBSTRING)));

	return (listOfCandidate[position]);
};

float Jacode_iiiAudioProcessor::PitchEstimator(const std::vector<float>& segment)
{
	return (HarmonicSummation(segment, f0LimitsInf, f0LimitsSup, nbOfHarmonicsInit, getSampleRate()));
};

std::vector<int> Jacode_iiiAudioProcessor::BetaEstimator(const std::vector<float>& segment, float observedPitch, std::vector<float> pitchReference, float& f0Features, float& betaFeatures)
{
	std::vector<int> finalCandidate;

	PitchCandidate(finalCandidate, observedPitch, pitchReference, NBSTRING);

	std::vector<float> pitchEstimate;
	std::vector<float> BEstimate;
	std::vector<float> costFunctionMaxVal;

	for (int indexCandidate = 0; indexCandidate < finalCandidate.size(); ++indexCandidate)
	{
		double minBetaGrid = std::numeric_limits<double>::max();
		double maxBetaGrid = std::numeric_limits<double>::min();

		for (int k = 0; k < NOMBRE_ESTIMATION; ++k)
		{
			if (minBetaGrid > featureMatrix[k](1, finalCandidate[indexCandidate]))
			{
				minBetaGrid = featureMatrix[k](1, finalCandidate[indexCandidate]);
			}
			if (maxBetaGrid < featureMatrix[k](1, finalCandidate[indexCandidate]))
			{
				maxBetaGrid = featureMatrix[k](1, finalCandidate[indexCandidate]);
			}
		}

		pitchEstimate.push_back(0);
		BEstimate.push_back(0);
		costFunctionMaxVal.push_back(0);

		InharmonicSummation(segment, observedPitch, highNbOfHarmonics, (int)getSampleRate(), maxBetaGrid, minBetaGrid, betaRes, lengthFFT, pitchEstimate[indexCandidate], BEstimate[indexCandidate], costFunctionMaxVal[indexCandidate]);
	}

	int maxIndexCost(0);
	double maxvalueCost(0);

	for (int i = 0; i < costFunctionMaxVal.size(); ++i)
	{
		if (costFunctionMaxVal[i] > maxvalueCost)
		{
			maxvalueCost = costFunctionMaxVal[i];
			maxIndexCost = i;
		}
	}

	f0Features = pitchEstimate[maxIndexCost];
	betaFeatures = BEstimate[maxIndexCost];

	return(finalCandidate);
};

float Jacode_iiiAudioProcessor::PluckPositionPrediction(const std::vector<std::complex<float>>& hilberOutput, float f0Features, float betaFeatures, int fretPlayed)
{
	std::vector<double> amplitudesAbs;
	double pluckingPosition;
	double L;

	amplitudesAbs = AmplitudesEstimation(hilberOutput, f0Features, betaFeatures, hilberOutput.size(), (int)getSampleRate(), highNbOfHarmonics);
	L = LOpen * pow(2.0, (-(double)fretPlayed / 12.0));
	pluckingPosition = PluckingPositionEstimatorLSD(amplitudesAbs, L);

	return(pluckingPosition);
};

//==============================================================================////==============================================================================//
//===============================------------===================================////===============================------------===================================//
//==============================================================================////==============================================================================//

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

void Jacode_iiiAudioProcessor::setCurrentProgram(int index)
{
}

const String Jacode_iiiAudioProcessor::getProgramName(int index)
{
	return {};
}

void Jacode_iiiAudioProcessor::changeProgramName(int index, const String& newName)
{
}

//==============================================================================//
//==============================================================================//
//==============================================================================//

void Jacode_iiiAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
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
bool Jacode_iiiAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
#if JucePlugin_IsMidiEffect
	ignoreUnused(layouts);
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

bool Jacode_iiiAudioProcessor::hasEditor() const
{
	return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* Jacode_iiiAudioProcessor::createEditor()
{
	return new Jacode_iiiAudioProcessorEditor(*this);
}

//==============================================================================
void Jacode_iiiAudioProcessor::getStateInformation(MemoryBlock& destData)
{
	// You should use this method to store your parameters in the memory block.
	// You could do that either as raw data, or use the XML or ValueTree classes
	// as intermediaries to make it easy to save and load complex data.
}

void Jacode_iiiAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
	// You should use this method to restore your parameters from this memory block,
	// whose contents will have been created by the getStateInformation() call.
}

//==============================================================================//
//===============================------------===================================//
//==============================================================================//

//==============================================================================//
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
	return new Jacode_iiiAudioProcessor();
}
//==============================================================================//

//==============================================================================////==============================================================================//
//===============================----END-----===================================////===============================----END----===================================//
//==============================================================================////==============================================================================//

































/*




//thread t1(&Jacode_iiiAudioProcessor::Prediction, std::ref(newSegment), std::ref(afficheValue));//, std::ref(f0LimitsInf), std::ref(f0LimitsSup), std::ref(nbOfHarmonicsInit), getSampleRate());
//auto myFuture = std::async(std::launch::async, Predict, std::ref(newSegment), std::ref(afficheValue));


//std::thread t([&](Jacode_iiiAudioProcessor* P) { P->threadPredict(newSegment, prediction); }, this);

//thread//

tOnset.detach();
tPrediction.detach();
tDisplay.detach();


void Jacode_iiiAudioProcessor::threadPredict (std::vector<float> const& newSegment, int& prediction)
{
	std:vector<float> hereWeGo = newSegment;
	int prediction2 = 8;//FretStringPrediction(hereWeGo);
};

void Jacode_iiiAudioProcessor::OnsetThread()
{
	while (!EndProcess)
	{
		//std::unique_lock<std::mutex> lck(bufferReady);
		while (!bufferEmpty) {};// waitBufferReady.wait(lck); }

		//======================================//
		afficheValue = 54;
		//std::cout << "onset detection in process... " << std::endl;

		//======================================//

		/* Notify next threads it is its turn */
/*
while (segmentToProcess) {};
segmentToProcess = true;
//waitOnsetReady.notify_all();
bufferEmpty = false;
	}
}

void Jacode_iiiAudioProcessor::PredictionThread()
{
	while (!EndProcess)
	{
		//std::unique_lock<std::mutex> lck(onsetReady);
		while (!segmentToProcess) {};//waitOnsetReady.wait(lck); }

		//======================================//

		//std::cout << "Prediction in process... " << std::endl;

		//======================================//

		while (noteToProcess) {};
		noteToProcess = true;
		//waitNoteReady.notify_all();
		segmentToProcess = false;
	}
}

void Jacode_iiiAudioProcessor::DisplayThread()
{
	while (!EndProcess)
	{
		//std::unique_lock<std::mutex> lck(noteReady);
		while (!noteToProcess) {};// waitNoteReady.wait(lck);}


		//======================================//

		//std::cout << "diplay in process... " << std::endl;
		afficheValue = 85.0;

		//======================================//

		goToDisplay = true;
		noteToProcess = false;
	}
}

void Jacode_iiiAudioProcessor::RunThread()
{
	while (bufferEmpty) {};
	bufferEmpty = true;
	//waitBufferReady.notify_all();
}

			//std::thread t([&](Jacode_iiiAudioProcessor* P) { P->ProcessWithSpectrogramOnsetDetector(buffer.getReadPointer(channel), buffer.getNumSamples(), getSampleRate()); }, this);
			//t.detach();



	//thread
	ready (false),
	bufferEmpty (false),
	EndProcess (false),
	segmentToProcess (false),
	noteToProcess (false),
	goToDisplay (false),
	tOnset([&](Jacode_iiiAudioProcessor* P) { P->OnsetThread(); }, this),
	tPrediction([&](Jacode_iiiAudioProcessor* P) { P->PredictionThread(); }, this),
	tDisplay([&](Jacode_iiiAudioProcessor* P) { P->DisplayThread(); }, this),



	//thread//
	std::mutex bufferReady;
	std::mutex onsetReady;
	std::mutex noteReady;

	std::condition_variable waitBufferReady;
	std::condition_variable waitOnsetReady;
	std::condition_variable waitNoteReady;

	bool ready;
	bool bufferEmpty;
	bool EndProcess;
	bool segmentToProcess;
	bool noteToProcess;
	bool goToDisplay;

	std::thread tOnset;
	std::thread tPrediction;
	std::thread tDisplay;

	void Jacode_iiiAudioProcessor::ProcessWithLowCPUDetector(const float* data, int bufSize, int sampleRate)
{
	int segmentSize = TIMEBUFFER * sampleRate; //40ms en sample
	int prediction(0);

	//onset detector/

	//Onset(data, bufSize, sampleRate, timeOfOnset, counter, storageActual, storagePast, thresholdValue, onsetScope, detectionDone);

	if (detectionDone)//many buffer are just store => predictions aren't done every time we get buffer but every time we get 40ms.
	{

		for (int i = 0; i < timeOfOnset.size(); ++i)
		{
			if (timeOfOnset[i] == true) //detection of an Onset
			{
				if (onsetdetected == false)
				{
					onsetdetected = true;
					newSegment.clear();
					newSegment.push_back(storagePast[double(i) + segmentSize]);

				}
				else if (newSegment.size() < segmentSize)//if a new onset is detect earlier
				{

					prediction = FretStringPrediction(newSegment);
					fretPlayed = (unsigned int)floor(prediction / NBSTRING);
					stringPlayed = prediction % NBSTRING;

					newSegment.clear(); //we create a new segment
					newSegment.push_back(storagePast[double(i) + segmentSize]);

					//we continue to fill the new segment

				}
				else // no new onset and the segment have its 40ms but this will never arrive normaly cause we will not detect between 40ms buffer
				{
					prediction = FretStringPrediction(newSegment);
					fretPlayed = (unsigned int)floor(prediction / NBSTRING);
					stringPlayed = prediction % NBSTRING;
					onsetdetected = false;
				}
			}
			else
			{

				if ((newSegment.size() < segmentSize) && (onsetdetected == true)) // no detection and no 40ms so just fill the segment
				{

					newSegment.push_back(storagePast[double(i) + segmentSize]);

				}
				else if (onsetdetected == true) //segment ready so predict and clear
				{

					prediction = FretStringPrediction(newSegment);
					fretPlayed = (unsigned int)floor(prediction / NBSTRING);
					stringPlayed = prediction % NBSTRING;

					onsetdetected = false;
					newSegment.clear();
				}
				else
				{
					//do nothing
				}

			}

			pushNextSampleIntoFifo(timeOfOnset[i], fretPlayed, stringPlayed, (4096), onsetScope[i]);

		}//end process//

	}//if onset

}

*/