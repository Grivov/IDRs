% Run script for generating Initial and In files for Chain450
% Correct sequence: Units -> IdealChain -> StickerSpacer -> Input files
clear; clc;

fprintf('\n=== Running Units Setup ===\n');
Units_nano_Chain450

fprintf('\n=== Running Ideal Chain Generation (Random Walks) ===\n');
InitialState_IdealChain_Chain450

fprintf('\n=== Running StickerSpacer Initial State Generation ===\n');
InitialState_StickerSpacer_Chain450

fprintf('\n=== Running Input File Generation ===\n');
In_StickerSpacer_Chain450_Equilibrium

fprintf('\n=== Complete ===\n');
