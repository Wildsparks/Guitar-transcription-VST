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
#define TIMEFRAME  0.004
#define SCOPESIZE  4096
#define NBFRET     21
#define NBSTRING   6

struct Model;
struct result;

//==============================================================================////==============================================================================//
Jacode_iiiAudioProcessor::Jacode_iiiAudioProcessor() : //class audioprocessor constructor

	//-------first--step-------//

	//--------//onset//--------//

	counterOnset(0),
	thresholdValue(0),

	storageSparse(120 * 10, 0),
	storageActual(1, 0),
	storagePast(1, 0),
	newSegment(1, 0),
	onsetScope(1, 0),

	detectionDone(false),
	onsetdetected(false),

	timeOfOnset(1, 0),

	//-------second-step-------//

	//------//prepocess//------//

	lengthFFT(pow(2.0, 17.0)),  // Length of  zero - padded FFT. //2^19 before opti 
	segmentDuration(40e-3),       // segment duration in seconds.

	//--------third-step-------//

	//--------//pitch//--------//

	f0LimitsInf(75),
	f0LimitsSup(700),          // boundaries for f0 search grid in Hz.
	nbOfHarmonicsInit(5),      // number of harmonics for initial harmonic estimate(B = 0).

	//--------fourth-step------//

	//------//Candidate//------//

	//--------fifth-step--------//

	//------//inharmonic//------//

	betaRes(1e-5),              // resolution of search grid for B. //1e-7 before opti
	highNbOfHarmonics(20),      //25 before opti

	//--------sixth-step--------//

	//------//Prediction//------//

	distribution(0.0, 1.0),
	X(2),
	trueString(6),
	trueFret(132),
	//classifier parameters//
	mu(2),
	C(2, 2),
	JMatrix(132),
	pitchReference(132),
	w0(6, 22), //= squeeze(est_f0);
	B(6, 22),  //= squeeze(BCoeff);
	P(0),
	prediction(0),
	pluckValue(0.0),

	//-------seventh-step-------//

	//----//pluck position//----//

	LOpen(64.3),
	distanceMaxPluck(LOpen / 2.0),
	resPluck(0.1),

	//-------final-step-------//

	//--------//plot//--------//

	DataScopeFifo(),
	nextBlockReady(false),
	stringPlayed(0),
	fretPlayed(0),
	isPlayed(false),

	//-----miscellaneous-----//
	PluckOnOff(false),
	afficheValue(10),
	//--------//---//--------//

	

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

	classifier = Calcul_parametre(featureMatrix, 0.010, 0.013, 0.0172, 0.026, 0.036, 0.046, 0.010, 0.013, 0.0172, 0.0146, 0.016, 0.018, 16.5, 16.0, 17.5, 18.5, 20, 17, 7950.0, 6000.0, 79.3e9, 2.27e11, generator, distribution);

	pitchReference.resize(132);

	for (int i = 0; i < 132; ++i)
	{
		pitchReference[i] = classifier.Mu(0, i);
	}

	for (int i = 0; i < SCOPESIZE; ++i)
	{
		scopeDataIsPlayed[i] = false;
		scopeDataFret[i]     = 0;
		scopeDataString[i]   = 0;
		scopeDataOnset[i]    = 0.0;  
		scopeDataPluck[i]    = 0.0;
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
			//=========say to Onset to go========//
			//===================================//

			int buffsize(buffer.getNumSamples());
			int counterOfBuffer(1);
			std::vector<double> data;
			const float* pointerData(buffer.getReadPointer(channel));

			while (buffsize > TIMEBUFFER * getSampleRate()) //if buffer > 40ms we are not able to predict for each time we need.
			{
				buffsize /= 2.0;
				counterOfBuffer*=2.0;
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
//==============================================================================////==============================================================================//
//==============================================================================////==============================================================================//

//==============================================================================//
//================================--Process--===================================//
//==============================================================================//

void Jacode_iiiAudioProcessor::ProcessWithSpectrogramOnsetDetector(const double* data, int bufSize, double sampleRate)
{
	int    segmentSize     (static_cast<int>(TIMEBUFFER * sampleRate)); //40ms in sample
	int    frameSize       (static_cast<int>(TIMEFRAME  * sampleRate)); //1ms in sample
	int    Noteprediction  (0);
	double time            (-0.040 * 2.0);            //decay time because of storing

	Onset(data, bufSize, sampleRate, timeOfOnset, storageActual, storagePast, thresholdValue, onsetScope, detectionDone, storageSparse);
	
	if (detectionDone)//many buffer are just store => predictions aren't done every time we get buffer but every time we get 40ms.
	{

		int index = (40.0+24.0)/4.0; //decay for take temporal sample from the time onset is detect

		for (int i = 0; i < timeOfOnset.size(); ++i)
		{
			time += TIMEFRAME;
			fretPlayed = -1;
			stringPlayed = -1;
			isPlayed = false;

			if (timeOfOnset[i] == true) //detection of an Onset
			{

				if (onsetdetected == false)
				{
					onsetdetected = true;
					newSegment.clear();

					for (int j = 0; j < frameSize; ++j)//fill with the first frame
					{
						newSegment.push_back(storagePast[(((__int64)i + (__int64)index) * (__int64)frameSize + (__int64)j)]);
					}

				}
				else if (newSegment.size() < segmentSize)//if a new onset is detect earlier
				{

					FretStringPrediction(newSegment, Noteprediction, pluckValue);
					fretPlayed = (int)floor(Noteprediction / NBSTRING);
					stringPlayed = Noteprediction % NBSTRING;
					isPlayed = true;

					newSegment.clear(); //we create a new segment

					for (int j = 0; j < frameSize; ++j)//we continue to fill the new segment
					{
						newSegment.push_back(storagePast[(((__int64)i + (__int64)index) * (__int64)frameSize + (__int64)j)]);
					}
				}
				else // no new onset and the segment have its 40ms but this will never arrive normaly cause we will not detect between 40ms buffers
				{
					FretStringPrediction(newSegment, Noteprediction, pluckValue);
					fretPlayed = (int)floor(Noteprediction / NBSTRING);
					stringPlayed = Noteprediction % NBSTRING;
					isPlayed = true;

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
						newSegment.push_back(storagePast[(((__int64)i + (__int64)index) * (__int64)frameSize + (__int64)j)]);
					}

				}
				else if (onsetdetected == true) //segment ready : predict and clear
				{

					FretStringPrediction(newSegment, Noteprediction, pluckValue);
					stringPlayed = Noteprediction % NBSTRING;
					fretPlayed = (int)floor(Noteprediction / NBSTRING);
					isPlayed = true;

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
				pushNextSampleIntoFifo(isPlayed, fretPlayed, stringPlayed, (4096), 0, pluckValue);
			}
			else
			{
				pushNextSampleIntoFifo(isPlayed, fretPlayed, stringPlayed, (4096), onsetScope[i], pluckValue);
			}

		}  //----------------//end process//----------------//

	}      //----------------//endif onset//----------------//
	
}

void Jacode_iiiAudioProcessor::ProcessWithLowCPUOnsetDetector(const double* data, int bufSize, double sampleRate)
{
	int    segmentSize(static_cast<int>(TIMEBUFFER * sampleRate)); //40ms in sample
	int    frameSize(static_cast<int>(TIMEFRAME * sampleRate)); //1ms in sample
	int    Noteprediction(0);
	double time(-0.040 * 2.0);            //decay time because of storing

	Onset(data, bufSize, sampleRate, timeOfOnset, storageActual, storagePast, thresholdValue, onsetScope, detectionDone, storageSparse);

	if (detectionDone)//many buffer are just store => predictions aren't done every time we get buffer but every time we get 40ms.
	{

		int index = 40; //decay for take temporal sample from the time onset is detect

		for (int i = 0; i < timeOfOnset.size(); ++i)
		{
			time += TIMEFRAME;
			fretPlayed = -1;
			stringPlayed = -1;
			isPlayed = false;

			if (timeOfOnset[i] == true) //detection of an Onset
			{

				if (onsetdetected == false)
				{
					onsetdetected = true;
					newSegment.clear();

					for (int j = 0; j < frameSize; ++j)//fill with the first frame
					{
						newSegment.push_back(storagePast[((__int64)i + (__int64)index * (__int64)frameSize + (__int64)j)]);
					}

				}
				else if (newSegment.size() < segmentSize)//if a new onset is detect earlier
				{

					FretStringPrediction(newSegment, Noteprediction, pluckValue);
					fretPlayed = (unsigned int)floor(Noteprediction / NBSTRING);
					stringPlayed = Noteprediction % NBSTRING;
					isPlayed = true;

					newSegment.clear(); //we create a new segment

					for (int j = 0; j < frameSize; ++j)//we continue to fill the new segment
					{
						newSegment.push_back(storagePast[((__int64)i + (__int64)index * (__int64)frameSize + (__int64)j)]);
					}



				}
				else // no new onset and the segment have its 40ms but this will never arrive normaly cause we will not detect between 40ms buffers
				{
					FretStringPrediction(newSegment, Noteprediction, pluckValue);
					fretPlayed = (unsigned int)floor(Noteprediction / NBSTRING);
					stringPlayed = Noteprediction % NBSTRING;
					isPlayed = true;

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
						newSegment.push_back(storagePast[((__int64)i + (__int64)index * (__int64)frameSize + (__int64)j)]);
					}

				}
				else if (onsetdetected == true) //segment ready : predict and clear
				{

					FretStringPrediction(newSegment, Noteprediction, pluckValue);
					stringPlayed = Noteprediction % NBSTRING;
					fretPlayed = (int)floor(Noteprediction / NBSTRING);
					isPlayed = true;

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
				pushNextSampleIntoFifo(isPlayed, fretPlayed, stringPlayed, (4096), 0, pluckValue);
			}
			else
			{
				pushNextSampleIntoFifo(isPlayed, fretPlayed, stringPlayed, (4096), onsetScope[i], pluckValue);
			}

		}  //----------------//end process//----------------//

	}      //----------------//endif onset//----------------//

}

//==============================================================================//
//=================================Analyser=====================================//
//==============================================================================//

void Jacode_iiiAudioProcessor::pushNextSampleIntoFifo(bool isPlayed, int fretPlayed, int stringPlayed, int size, double onsetValue, double pluckValue) noexcept
{

	nextBlockReady = false;

	if (DataScopeFifo.isPlayed.size() != size)
	{
		DataScopeFifo.Fret.resize(size);
		DataScopeFifo.String.resize(size);
		DataScopeFifo.isPlayed.resize(size);
		DataScopeFifo.Onset.resize(size);
		DataScopeFifo.Pluck.resize(size);
	}

	std::rotate(DataScopeFifo.Fret.begin(),     DataScopeFifo.Fret.begin() + 1,     DataScopeFifo.Fret.end());
	std::rotate(DataScopeFifo.String.begin(),   DataScopeFifo.String.begin() + 1,   DataScopeFifo.String.end());
	std::rotate(DataScopeFifo.isPlayed.begin(), DataScopeFifo.isPlayed.begin() + 1, DataScopeFifo.isPlayed.end());
	std::rotate(DataScopeFifo.Onset.begin(),    DataScopeFifo.Onset.begin() + 1,    DataScopeFifo.Onset.end());
	std::rotate(DataScopeFifo.Pluck.begin(),    DataScopeFifo.Pluck.begin() + 1,    DataScopeFifo.Pluck.end());

	DataScopeFifo.isPlayed[(__int64)size - (__int64)1] = isPlayed;
	DataScopeFifo.String  [(__int64)size - (__int64)1] = stringPlayed;
	DataScopeFifo.Fret    [(__int64)size - (__int64)1] = fretPlayed;
	DataScopeFifo.Onset   [(__int64)size - (__int64)1] = onsetValue;
	DataScopeFifo.Pluck   [(__int64)size - (__int64)1] = pluckValue;

	nextBlockReady = true;

}

void Jacode_iiiAudioProcessor::drawNextFrame()
{
	for (int i = 0; i < SCOPESIZE; ++i)
	{
		scopeDataIsPlayed[i]   = DataScopeFifo.isPlayed  [(int)round(i * DataScopeFifo.isPlayed.size() / SCOPESIZE)];
		scopeDataString[i]     = DataScopeFifo.String    [(int)round(i * DataScopeFifo.String.size()   / SCOPESIZE)];
		scopeDataFret[i]       = DataScopeFifo.Fret      [(int)round(i * DataScopeFifo.Fret.size()     / SCOPESIZE)];
		scopeDataOnset[i]      = DataScopeFifo.Onset     [(int)round(i * DataScopeFifo.Onset.size()    / SCOPESIZE)];
		scopeDataPluck[i]      = DataScopeFifo.Pluck     [(int)round(i * DataScopeFifo.Pluck.size()    / SCOPESIZE)];
	}
}

bool Jacode_iiiAudioProcessor::getNextBlockReady()
{
	return nextBlockReady;
}

void Jacode_iiiAudioProcessor::setNextBlockReady(bool value)
{
	nextBlockReady = value;
}

void Jacode_iiiAudioProcessor::drawFrame(Graphics& g)
{
	int width = 1600;
	int height = 700;

	if (scopeDataIsPlayed[0])
	{
		g.setFont(40.0f);
		g.drawSingleLineText(std::to_string(scopeDataFret[0]), jmap(0, 0, SCOPESIZE - 1 - 5, 60, width - 40), (int)(0.25 * height + 58.0 + 10.0 + 4.0 - (8.0 - (5.0 - scopeDataString[0])) * 20.0));
	}

	for (int i = 1; i < SCOPESIZE; ++i)
	{

		if (scopeDataIsPlayed[i])
		{
			g.setFont(40.0f);
			g.drawSingleLineText(std::to_string(scopeDataFret[i]), jmap(i, 0, SCOPESIZE - 1 -5, 60, width - 40), (int)(0.25 * height + 58.0 + 10.0 + 4.0 - (8.0 - (5.0-scopeDataString[i])) * 20.0));
		}
		//pluck position
		g.drawLine({ (float)jmap(i - 1, 0, SCOPESIZE - 1, 80, width - 40) , float(400.0 - (scopeDataPluck[i]) * 3.0),
			         (float)jmap(i ,    0, SCOPESIZE - 1, 80, width - 40) , float(400.0 - (scopeDataPluck[i]) * 3.0)});
		//onset
		g.drawLine({ (float)jmap(i - 1, 0, SCOPESIZE - 1, 40, width-40) , jmap((float)scopeDataOnset[i - 1], -1.0f, 100.0f, (float)height-60.0f, 480.0f),
				     (float)jmap(i ,    0, SCOPESIZE - 1, 40, width-40) , jmap((float)scopeDataOnset[i],     -1.0f, 100.0f, (float)height-60.0f, 480.0f) });
		//treshold
		g.drawLine(40, jmap(float(thresholdValue) / 100, 0.0f, 1.0f, (float)height - 140, 480.0f), (float)(width - 40.0), jmap(float(thresholdValue) / 100, 0.0f, 1.0f, (float)height - 140, 480.0f));

	}	
}


//==============================================================================//
//=================================accesors=====================================//
//==============================================================================//

void Jacode_iiiAudioProcessor::setThresholdValue(int value)
{
	thresholdValue = value;
}

void Jacode_iiiAudioProcessor::setnfftValue(double value)
{
	lengthFFT = pow(2.0, value);
}

void Jacode_iiiAudioProcessor::setbetaresValue(double value)
{
	betaRes = pow(10.0, -value);
}

void Jacode_iiiAudioProcessor::setPluckOnOff(bool value)
{
	PluckOnOff = value;
}

void Jacode_iiiAudioProcessor::setResPluck(double value)
{
	resPluck = std::pow(10.0,-value);
}

void Jacode_iiiAudioProcessor::setDistanceMaxPluck(double value)
{
	distanceMaxPluck = value; 
}

void Jacode_iiiAudioProcessor::setNBHarmonics(double value)
{
	highNbOfHarmonics = static_cast<int>(value);
}

void Jacode_iiiAudioProcessor::setNBHarmonicsF0(double value)
{
	nbOfHarmonicsInit = static_cast<int>(value);
}

double Jacode_iiiAudioProcessor::getAfficheValue()
{
	return(afficheValue);
}

void Jacode_iiiAudioProcessor::reloadModel(double val1, double val2, double val3, double val4, double val5,
	double val6, double val7, double val8, double val9, double val10,
	double val11, double val12, double val13, double val14, double val15,
	double val16, double val17, double val18, double val19, double val20,
	double val21, double val22)
{
	classifier = Calcul_parametre(featureMatrix, 
		 val1,   val2,   val3,   val4,   val5,
		 val6,   val7,   val8,   val9,   val10,
		 val11,  val12,  val13,  val14,  val15,
		 val16,  val17,  val18,  val19,  val20,
		 val21,  val22, generator, distribution);

	pitchReference.resize(132);

	for (int i = 0; i < 132; ++i)
	{
		pitchReference[i] = classifier.Mu(0, i);
	}
}

//==============================================================================//
//================================Prediction====================================//
//==============================================================================//
void Jacode_iiiAudioProcessor::FretStringPrediction(const std::vector<double>& segment, int& Noteprediction, double& pluckValue)
{
	double f0Features(0);
	double betaFeatures(0);

	std::vector<double> processSegment;
	std::vector<std::complex<double>> hilberOutput;
	processSegment = segment;

	Preprocessing(processSegment, hilberOutput, static_cast<int>(lengthFFT));

	double pitch = PitchEstimator(processSegment);

	std::vector<int> listOfCandidate = BetaEstimator(processSegment, pitch, pitchReference, f0Features, betaFeatures);

	X(0) = f0Features;
	X(1) = betaFeatures;

	int position(0);

	JMatrix.resize(listOfCandidate.size());
	for (int k = 0; k < listOfCandidate.size(); ++k)
	{
		int indexCandidat = listOfCandidate[k];
		mu(0) = classifier.Mu(0, indexCandidat);
		mu(1) = classifier.Mu(1, indexCandidat);
		C = classifier.Sigma[indexCandidat];
		P = classifier.W(indexCandidat);

		JMatrix(k) = JClassifier(mu, C, P, X);

	}
	JMatrix.maxCoeff(&position);

	Noteprediction = listOfCandidate[position];
	if (PluckOnOff)
	{
		pluckValue = PluckPositionPrediction(hilberOutput, f0Features, betaFeatures, (int)floor(listOfCandidate[position] / NBSTRING));
	}
};

double Jacode_iiiAudioProcessor::PitchEstimator(const std::vector<double>& segment)
{
	return (HarmonicSummation(segment, f0LimitsInf, f0LimitsSup, nbOfHarmonicsInit, getSampleRate()));
};

std::vector<int> Jacode_iiiAudioProcessor::BetaEstimator(const std::vector<double>& segment, double observedPitch, std::vector<double> pitchReference, double& f0Features, double& betaFeatures)
{
	std::vector<int> finalCandidate;

	PitchCandidate(finalCandidate, observedPitch, pitchReference, NBSTRING);

	std::vector<double> pitchEstimate;
	std::vector<double> BEstimate;
	std::vector<double> costFunctionMaxVal;

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

		InharmonicSummation(segment, observedPitch, highNbOfHarmonics, getSampleRate(), maxBetaGrid, minBetaGrid, betaRes, lengthFFT, pitchEstimate[indexCandidate], BEstimate[indexCandidate], costFunctionMaxVal[indexCandidate]);
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

double Jacode_iiiAudioProcessor::PluckPositionPrediction(const std::vector<std::complex<double>>& hilberOutput, double f0Features, double betaFeatures, int fretPlayed)
{
	std::vector<double> amplitudesAbs;
	double pluckingPosition;
	double L;

	amplitudesAbs = AmplitudesEstimation(hilberOutput, f0Features, betaFeatures, static_cast<int>(hilberOutput.size()), getSampleRate(), highNbOfHarmonics);
	L = LOpen * pow(2.0, (-(double)fretPlayed / 12.0));
	pluckingPosition = PluckingPositionEstimatorLSD(amplitudesAbs, L, distanceMaxPluck, resPluck);

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
