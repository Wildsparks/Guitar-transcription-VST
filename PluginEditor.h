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
class Jacode_iiiAudioProcessorEditor : public AudioProcessorEditor,                //hérite de audioprocessoreditor
	                                   public Timer,
	                                   public Slider::Listener,
									   public Button::Listener,
									   public TextEditor::Listener
{
public:  
    Jacode_iiiAudioProcessorEditor (Jacode_iiiAudioProcessor&);                     //constructor
    ~Jacode_iiiAudioProcessorEditor();                                              //destructor 

	//=================================METHODE======================================
    void paint (Graphics&) override;                                                
    void resized() override;
	void timerCallback() override;
	void sliderValueChanged(Slider* slider) override;
	void buttonClicked(Button* button) override;

private:
	

	//================================ATTRIBUT======================================
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    Jacode_iiiAudioProcessor& processor;

	//background//
	Image background;

	//onset treashold//
	Slider thresholdSlider;
	Label thresholdLabel;

	//precision
	Slider betaresSlider;
	Label betaresLabel;

	Slider nfftSlider;
	Label nfftLabel;
	
	Slider distanceMaxPluck;
	Label distanceMaxPluckLabel;

	Slider resPluck;
	Label resPluckLabel;

	Slider nbHarmoF0;
	Label nbHarmoF0Label;

	Slider nbHarmo;
	Label nbHarmoLabel;

	ToggleButton switchPluck;
	Label switchPluckLabel;

	TextButton reloadButton;

	Label DfullLabel;
	Label DcoreLabel;
	Label TensionLabel;
	Label rhoCore;
	Label rhoWrapping;
	Label steelShearModulus;
	Label steelYoungModulus;

	TextEditor string_1;
	TextEditor string_2;
	TextEditor string_3;
	TextEditor string_4;
	TextEditor string_5;
	TextEditor string_6;
	TextEditor string_7;
	TextEditor string_8;
	TextEditor string_9;
	TextEditor string_10;
	TextEditor string_11;
	TextEditor string_12;

	TextEditor string_13;
	TextEditor string_14;
	TextEditor string_15;
	TextEditor string_16;
	TextEditor string_17;
	TextEditor string_18;

	TextEditor string_19;
	TextEditor string_20;
	TextEditor string_21;
	TextEditor string_22;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Jacode_iiiAudioProcessorEditor)
};
