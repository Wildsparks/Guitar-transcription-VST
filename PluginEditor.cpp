/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <string> 

#define BPM 110

//==============================================================================
Jacode_iiiAudioProcessorEditor::Jacode_iiiAudioProcessorEditor (Jacode_iiiAudioProcessor& p) //constructor
    : AudioProcessorEditor (&p), processor (p), value(0)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.


	//size of the windows
    setSize (1600, 600);

	//FPS
	startTimerHz(120);  //interval en miliseconde divisé par 8 pour faire des 60000 / BPM / 8

	//background
	background = ImageCache::getFromMemory(BinaryData::backgroundImage_png, BinaryData::backgroundImage_pngSize);

	//--------------------------------//
	//--------onset sliders-----------//
	//--------------------------------//

	//--------------------------------//
	addAndMakeVisible(thresholdSlider);
	thresholdSlider.setRange(0, 100);            
	thresholdSlider.setTextValueSuffix(" %");    
	thresholdSlider.setValue(5.0);
	thresholdSlider.addListener(this);
	thresholdSlider.setColour(Slider::textBoxTextColourId, Colours::black);
	//--------------------------------//
	addAndMakeVisible(thresholdLabel);
	thresholdLabel.setText("Threshold", dontSendNotification);
	thresholdLabel.attachToComponent(&thresholdSlider, true); 
	thresholdLabel.setColour(Label::textColourId, Colours::black);
	//--------------------------------//



	labelvalue.setFont(60.0f);
	labelvalue.setText("load...", sendNotification);
	addAndMakeVisible(&labelvalue);

}

Jacode_iiiAudioProcessorEditor::~Jacode_iiiAudioProcessorEditor()                            //destructor 
{

}

//=================================METHODE======================================
void Jacode_iiiAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
    g.setColour (Colours::darkgrey);
    g.setFont (15.0f);
	
	g.drawImage(background,
		0,
		0,
		1600,
		600,
		0,
		0,
		1641,
		501
	);
	processor.drawFrame(g);


}

void Jacode_iiiAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	auto sliderLeft = 820;
	labelvalue.setBounds(sliderLeft, 400, getWidth() - sliderLeft, 100);
	thresholdSlider.setBounds(sliderLeft, 500, getWidth() - sliderLeft - 180, 80);
	
}

void Jacode_iiiAudioProcessorEditor::timerCallback()
{
	labelvalue.setText(std::to_string(processor.getAfficheValue()), sendNotification);
	if (processor.getNextFFTBlockReady())
	{
		labelvalue.setText(std::to_string(processor.getAfficheValue()), sendNotification);
		processor.setThresholdValue((int)thresholdSlider.getValue());
		processor.drawNextFrameOfSpectrum();
		repaint();
		processor.setNextFFTBlockReady(false);
		
	}

}

void Jacode_iiiAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
	if (slider == &thresholdSlider)
	{
		thresholdSlider.setValue(thresholdSlider.getValue(), dontSendNotification);
	}
		
}
