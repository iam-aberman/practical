raw_accelerometer_data = load('ins-mems-aks1.log', '-ascii');
accelerometer_data = raw_accelerometer_data(:, 5:7);

figure;
hold on;
plot(accelerometer_data(:, 1), 'r');
plot(accelerometer_data(:, 2), 'g');
plot(accelerometer_data(:, 3), 'b');

% Calibration
positions_count = 6;

% Starts and ends of segments are chosen "empirically"
segments_starts = [1, 700, 1200, 1700, 2200, 2700];
segments_ends   = [500, 1100, 1500, 2100, 2500, 3200];

for i = 1:positions_count
  avg_s_1(i) = sum(accelerometer_data(segments_starts(i):segments_ends(i), 1)) / (segments_ends(i) - segments_starts(i) + 1);
  avg_s_2(i) = sum(accelerometer_data(segments_starts(i):segments_ends(i), 2)) / (segments_ends(i) - segments_starts(i) + 1); 
  avg_s_3(i) = sum(accelerometer_data(segments_starts(i):segments_ends(i), 3)) / (segments_ends(i) - segments_starts(i) + 1);
endfor

b_a(1) = 1 / 2 * (avg_s_1(1) + avg_s_1(2));  
b_a(2) = 1 / 2 * (avg_s_2(1) + avg_s_2(2));  
b_a(3) = 1 / 2 * (avg_s_3(1) + avg_s_3(2));

c = 1 / (2 * 9.8);
S_a(1, :) = c * [(avg_s_1(6) - avg_s_1(5)); (avg_s_1(4) - avg_s_1(3)); (avg_s_1(1) - avg_s_1(2))];
S_a(2, :) = c * [(avg_s_2(6) - avg_s_2(5)); (avg_s_2(4) - avg_s_2(3)); (avg_s_2(1) - avg_s_2(2))];
S_a(3, :) = c * [(avg_s_3(6) - avg_s_3(5)); (avg_s_3(4) - avg_s_3(3)); (avg_s_3(1) - avg_s_3(2))];

S_f = inv(-S_a);
b_f = -S_f * b_a';
f_z = (S_f * accelerometer_data' + b_f)'

figure;
hold on;
plot(f_z);
