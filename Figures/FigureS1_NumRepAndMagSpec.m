close all;

% Sample DNA sequence (replace with actual sequence)
dnaSequence = 'GTGCTTAACACATGCAAGTCGAACGATGATCCCAGCTTGCTGGGGGATTAGTGGCGAACGGGTGAGTAACACGTGAGTAACC';

% Convert DNA sequence to numerical PP representation
% Purine (A,G) = -1, Pyrimidine (T,C) = 1
ppSignal = zeros(1, length(dnaSequence));
ppSignal(ismember(dnaSequence, 'AGag')) = -1;
ppSignal(ismember(dnaSequence, 'TCtc')) = 1;

% Compute the magnitude spectrum of the PP signal
% This is a placeholder - you would use your actual data and method here
magnitudeSpectrum = abs(fft(ppSignal));

% Define colors for each nucleotide
colors = containers.Map({'A', 'T', 'G', 'C'}, {'#c22010', '#082ac0', '#bcaf0a', '#00bc0c'}); % Red, Blue, Green, Magenta

% Create the figure
figure('units','normalized','outerposition',[0 0 1 1]);

% Subplot 1: DNA Sequence
subplot(3, 1, 1);
axis([1 length(dnaSequence) -1 1]); % Set axis limits to match the stem plot
set(gca, 'Visible', 'off'); % Turn off the axis
hold on;
% Place each nucleotide at a specific position
for i = 1:length(dnaSequence)
    text(i, 0, dnaSequence(i), 'HorizontalAlignment', 'center', 'FontName', 'Courier', 'FontWeight','bold', 'FontSize', 20, 'Color', colors(dnaSequence(i)));
end
hold off;

% Subplot 2: Genomic Signal (PP Numerical Representation)
subplot(3, 1, 2);
stem(1:length(ppSignal), ppSignal, 'LineWidth', 2, 'MarkerFaceColor', '#2c85c5');
xlim([1 length(ppSignal)]);
ylim([-1.2 1.2]);
xlabel('Position', 'FontSize', 18);  % X-axis label with font size 18
ylabel({'Assigned', 'Num. Value'}, 'FontSize', 18);  % Y-axis label with font size 18
title('PP Numerical Representation', 'FontSize', 18);  % Title of the plot with font size 18
yticks([-1 0 1]);
set(gca, 'FontSize', 18);

% Subplot 3: Magnitude Spectrum
subplot(3, 1, 3);
f = linspace(0, 2*pi, length(magnitudeSpectrum));
plot(f, magnitudeSpectrum, 'LineWidth', 2, 'MarkerFaceColor', '#2c85c5');
xlim([0 2*pi]);
xlabel('Frequency (rad/sample)', 'FontSize', 18);  % X-axis label with font size 18
ylabel('Magnitude', 'FontSize', 18);
title('Magnitude Spectrum', 'FontSize', 18);
xticks([0 3.173 2*pi]);
xticklabels({'0', '\pi', '2\pi'});
set(gca, 'FontSize', 18);
