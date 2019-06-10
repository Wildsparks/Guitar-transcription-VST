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

#define SCOPESIZE 40000
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

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

	//===================================MOI========================================

	double getAfficheValue();

	//================================ANALYSEUR=====================================

	void pushNextSampleIntoFifo(float sample, int size) noexcept;
	void drawNextFrameOfSpectrum();
	bool getNextFFTBlockReady();
	void setNextFFTBlockReady(bool value);
	void drawFrame(Graphics& g);

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

	//onset
	int counter;
	std::vector<float> storageActual;
	std::vector<float> storagePast;
	std::vector<double> time;

	//analyser//

	std::vector<float> DataScope;      
	std::vector<float> DataScopefifo;
	bool nextFFTBlockReady = false;       //This temporary boolean tells us whether the next FFT block is ready to be rendered.
	bool newBuffer = true;
	float scopeData[SCOPESIZE];           //The scopeData float array of size 512 will contain the points to display on the screen.
	float scopeData2[SCOPESIZE];
	int counterAnalyser;
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Jacode_iiiAudioProcessor)
};




