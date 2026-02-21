% Load the .NET speech assembly
NET.addAssembly('System.Speech');

% Create the speech synthesizer
speaker = System.Speech.Synthesis.SpeechSynthesizer;

% Set the voice to Microsoft Zira
speaker.SelectVoice('Microsoft Zira Desktop');

% Speak your custom message
Speak(speaker, 'Welcome to the O P C Calcium Analysis App. Where do you wish to begin?');
