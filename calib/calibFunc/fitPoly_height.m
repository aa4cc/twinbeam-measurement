function p_coeff = fitPoly_height(height, posDif, voltages)

% p_coeff  = polyfit(posDif, height, 1);
p_coeff = posDif'\height';
% Identify the dependence of

figure; clf;
plot(voltages, [height', p_coeff*posDif'], '-*')
legend('Measured', 'Estimated', 'Location', 'SouthEast')
xlabel('Voltage [V]')
ylabel('Levitation Height [\mum]')
grid on

figure; clf;
plot(posDif, [height', p_coeff*posDif'], '-*')
legend('Measured', 'Estimated', 'Location', 'SouthEast')
xlabel('Diffraction patterns distance [um]')
ylabel('Levitation Height [\mum]')
grid on