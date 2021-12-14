raw_angular_rate_data = load('ins-mems-dus2.log', '-ascii');
angular_rate_data = raw_angular_rate_data(:, 8:10);

figure;
hold on;
plot(angular_rate_data(:, 1), 'r');
plot(angular_rate_data(:, 2), 'g');
plot(angular_rate_data(:, 3), 'b');

% Calibration
positions_count = 6;

% Starts and ends of segments are chosen "empirically"
segments_starts = [50, 1000, 2050, 2950, 4350, 5100];
segments_ends   = [600, 1450, 2600, 3400, 4750, 5500];

% Introducing additional position, where phi = 0
b_w = [sum(angular_rate_data(3500:4000, 1)); sum(angular_rate_data(3500:4000, 2)); sum(angular_rate_data(3500:4000, 3))] / 501;
for i = 1:3
  angular_rate_data(:, i) -= b_w(i);
endfor

for i = 1:positions_count
  dt = raw_angular_rate_data(segments_starts(i), 12) * 60; 
  dt += raw_angular_rate_data(segments_starts(i), 13);
  dt += raw_angular_rate_data(segments_starts(i), 14) * 0.001;
  dt -= raw_angular_rate_data(segments_ends(i), 12) * 60;
  dt -= raw_angular_rate_data(segments_ends(i), 13); 
  dt -= raw_angular_rate_data(segments_ends(i), 14) * 0.001;
  
  dt
  
  count = segments_ends(i) - segments_starts(i) + 1;
  phi_1(i) = dt * sum(angular_rate_data(segments_starts(i):segments_ends(i), 1)) / count;
  phi_2(i) = dt * sum(angular_rate_data(segments_starts(i):segments_ends(i), 2)) / count;
  phi_3(i) = dt * sum(angular_rate_data(segments_starts(i):segments_ends(i), 3)) / count;
end

S_w(1, :) = 1 / 360 * [(phi_1(6) - phi_1(5)); (phi_1(3) - phi_1(4)); (phi_1(1) - phi_1(2))];
S_w(2, :) = 1 / 360 * [(phi_2(6) - phi_2(5)); (phi_2(3) - phi_2(4)); (phi_2(1) - phi_2(2))];
S_w(3, :) = 1 / 360 * [(phi_3(6) - phi_3(5)); (phi_3(3) - phi_3(4)); (phi_3(1) - phi_3(2))];

phi = (-inv(S_w) * angular_rate_data');

figure;
hold on;
plot(phi(1, :), 'r');
plot(phi(2, :), 'g');
plot(phi(3, :), 'b');

