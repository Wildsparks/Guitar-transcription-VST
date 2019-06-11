/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class Jacode_iiiAudioProcessorEditor : public AudioProcessorEditor,                //hérite de audioprocessoreditor
	                                   public Timer,
	                                   public Slider::Listener
{
public:  
    Jacode_iiiAudioProcessorEditor (Jacode_iiiAudioProcessor&);                     //constructor
    ~Jacode_iiiAudioProcessorEditor();                                              //destructor 

	//=================================METHODE======================================
    void paint (Graphics&) override;                                                
    void resized() override;
	void timerCallback() override;
	void sliderValueChanged(Slider* slider) override;

private:
	

	//================================ATTRIBUT======================================
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    Jacode_iiiAudioProcessor& processor;

	//backgournd//
	Image background;

	//onset treashold//
	Slider thresholdSlider;
	Label thresholdLabel;

	Label labelvalue;
	double value;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Jacode_iiiAudioProcessorEditor)
};
