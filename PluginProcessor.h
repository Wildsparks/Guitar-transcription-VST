/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include <Eigen/Dense>

#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Estimation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Features.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Onset.h"

#define SCOPESIZE 4096
struct Model;

//==============================================================================
/**
*/
class Jacode_iiiAudioProcessor  : public AudioProcessor                             //héritage de audioprocessor
{
public:                                                                             //déclaration des méthodes et accesseur
    //==============================================================================
    Jacode_iiiAudioProcessor();                                                     //constructor
    ~Jacode_iiiAudioProcessor();                                                    //décontructor

    //=================================METHODE======================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;
/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include <Eigen/Dense>
#include <thread>
#include <mutex>

#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Estimation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Onset.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PreProcessing.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\Features.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\HarmonicSummation.h"
#include "D:\ecler\Documents\Cours\Ingenieur_4A\Stage\Jacode_III\Builds\VisualStudio2019\PitchCandidate.h"

#define SCOPESIZE 4096
#define NBFRET 12
#define NBSTRING 6

struct Model;

struct result
{
	std::vector<bool> isPlayed;
	std::vector<int> Fret;
	std::vector<int> String;
};

//==============================================================================
/**
*/
class Jacode_iiiAudioProcessor  : public AudioProcessor                             //héritage de audioprocessor
{
public:                                                                             //déclaration des méthodes et accesseur
    //==============================================================================
    Jacode_iiiAudioProcessor();                                                     //constructor
    ~Jacode_iiiAudioProcessor();                                                    //décontructor

    //=================================METHODE======================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

	//==============================================================================//
	//================================Osciloscope===================================//
	//==============================================================================//

	void pushNextSampleIntoFifo(bool isPlayed, int fretPlayed, int stringPlayed, int size) noexcept;
	void drawNextFrameOfSpectrum();
	bool getNextFFTBlockReady();
	void setNextFFTBlockReady(bool value);
	void drawFrame(Graphics& g);

	//==============================================================================//
	//==================================Onset=======================================//
	//==============================================================================//

	void setThresholdValue(int value);

	//==============================================================================//
	//================================Prediction====================================//
	//==============================================================================//
	
	float PitchEstimator(std::vector<float> const& segment);
	float Jacode_iiiAudioProcessor::BetaEstimator();
	void Jacode_iiiAudioProcessor::FretStringPrediction(float observedPitch, std::vector<float> pitchReference, int nbFret, Model classifier);
	float Jacode_iiiAudioProcessor::PluckPositionPrediction();

	//==============================================================================//
	//===================================MOI========================================//
	//==============================================================================//

	double getAfficheValue();

private:                                                                            //déclaration des Attributs
    //================================ATTRIBUT======================================
	double afficheValue;
	Model classifier;

	Eigen::MatrixXd w0; //= squeeze(est_f0);
	Eigen::MatrixXd B;  //= squeeze(BCoeff);

	//sweep//
	int trueString;
	int trueFret;

	//input buffer//
	Eigen::VectorXd X;

	//classifier parameters//
	Eigen::VectorXd mu;
	Eigen::MatrixXd C;
	double P;
	Eigen::VectorXd JMatrix;
	double prediction;
	std::vector<float> pitchReference;

	//Initialize estimator / classifier implementation constants
	float        segmentDuration   = 40e-3; // segment duration in seconds.
	float        LOpen             = 64.3;  // assumed length of all open strings.
	unsigned int highNbOfHarmonics = 25;    // assumed high number of harmonics(M >> 1).Used for inharmonic pitch estimationand for estimation of plucking position on the fretted string.
	unsigned int nbOfHarmonicsInit = 5;     // number of harmonics for initial harmonic estimate(B = 0).
	unsigned int f0LimitsInf       = 75;
	unsigned int f0LimitsSup       = 700;   // boundaries for f0 search grid in Hz.
	float        lengthFFT         = pow(2,19);   // Length of  zero - padded FFT.
	float        betaRes           = 1e-5;  // resolution of search grid for B.

	//onset
	int counter;
	std::vector<float> storageActual;
	std::vector<float> storagePast;
	int thresholdValue;
	bool onsetdetected;
	std::vector<float> newSegment;

	//analyser//

	result DataScope;
	result DataScopeFifo;

	bool nextFFTBlockReady = false;       //This temporary boolean tells us whether the next FFT block is ready to be rendered.
	bool newBuffer = true;

	bool scopeDataIsPlayed[SCOPESIZE];           //The scopeData float array of size 512 will contain the points to display on the screen.
	int scopeDataFret[SCOPESIZE];
	int scopeDataString[SCOPESIZE];

	int counterAnalyser;

	//result
	int stringPlayed;
	int fretPlayed;
	std::vector<bool> timeOfOnset;

	//tempo
	std::vector<float> DataLastBuffer;

	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Jacode_iiiAudioProcessor)
};





