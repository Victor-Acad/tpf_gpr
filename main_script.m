% ------------------------------------------------------------------------

% |--------------------------------------|
% |                                      |
% |     Electromagnetismo Aplicado       |
% |                                      |
% |    TPF Georradar - Simulaciones      |
% |                                      |
% |--------------------------------------|

% Limpieza.
clear; close all; clc;

% Parámetros físicos.
global mu_0 eps_0 c;

mu_0 = 4*pi*1e-7;
eps_0 = 8.854e-12;
c = 1/sqrt(mu_0*eps_0);

% Frecuencia de operación.
f = 1000e6;
omega = 2*pi*f;

% Amplitud del campo eléctrico de la onda incidente emitida por la antena.
E_i_0 = 1;

% ------------------------------------------------------------------------

% Medio 1 (espacio libre).
sigma_1 = 0;
mu_r_1 = 1;
eps_r_1 = 1;

% Medio 2 (suelo de propiedades conocidas).
sigma_2 = 0.01;
mu_r_2 = 1;
eps_r_2 = 9;

% Medio 3 (medio a identificar, a priori desconocido).
sigma_3 = 0.05;
mu_r_3 = 1;
eps_r_3 = 16;

% Espesores de los medios.
d_1 = 0.5; % Distancia a la cual se emite desde el espacio libre.
d_2 = 0.5; % Espesor del medio 1 (intermedio).

% Matriz de propiedades de los medios (filas = medios; columnas = propiedades).
props = [
  sigma_1, mu_r_1, eps_r_1;
  sigma_2, mu_r_2, eps_r_2;
  sigma_3, mu_r_3, eps_r_3
];

% ------------------------------------------------------------------------

% Cálculo de parámetros de los medios.
function params_vec = calc_params(sigma, mu_r, eps_r, omega)
  global mu_0 eps_0;

  % Constante de propagación.
  gamma = sqrt((1i*omega*mu_r*mu_0) .* (sigma + 1i*omega*eps_r*eps_0));

  alpha = real(gamma);
  beta = imag(gamma);

  % Impedancia compleja.
  eta = sqrt((1i*omega*mu_r*mu_0) ./ (sigma + 1i*omega*eps_r*eps_0));

  params_vec = [alpha, beta, gamma, eta];
endfunction

% Matriz de parámetros (filas = medios; columnas = parámetros).
params = zeros(3, 4);

for i = 1:3
  params(i, :) = calc_params(props(i, 1), props(i, 2), props(i, 3), omega);
endfor

% ------------------------------------------------------------------------

% Vectores de parámetros.
alphas = params(:, 1);
betas = params(:, 2);
gammas = params(:, 3);
etas = params(:, 4);

% Velocidades de propagación.
vels = omega ./ betas;

% Coeficientes de reflexión.
Gamma_12 = (etas(2) - etas(1)) / (etas(1) + etas(2));
Gamma_23 = (etas(3) - etas(2)) / (etas(2) + etas(3));

Gamma_21 = -Gamma_12;
Gamma_32 = -Gamma_23;

% Coeficientes de transmisión.
tau_12 = 2 * etas(2) / (etas(1) + etas(2));
tau_21 = 2 * etas(1) / (etas(1) + etas(2));

tau_23 = 2 * etas(3) / (etas(2) + etas(3));
tau_32 = 2 * etas(2) / (etas(2) + etas(3));

% Reflexión total hacía el medio 1.
Gamma_T_1 = Gamma_12 + (tau_12*tau_21*Gamma_23*exp(-2*gammas(2)*d_2)) / (1 - Gamma_21*Gamma_23*exp(-2*gammas(2)*d_2));

% Transmisión total hacía el medio 3.
tau_T_3 = (tau_12*tau_23*exp(-gammas(2)*d_2)) / (1 - Gamma_21*Gamma_23*exp(-2*gammas(2)*d_2));

% Factor interviniente en la onda progresiva total en el medio 2.
tau_T_2 = tau_12 / (1 - Gamma_21*Gamma_23*exp(-2*gammas(2)*d_2));

% Factor interviniente en la onda regresiva total en el medio 2.
Gamma_T_2 = Gamma_23 * exp(-2*gammas(2)*d_2) * tau_T_2;

% ------------------------------------------------------------------------
% ----- Simulación 1: foco isotrópico puntual de ondas sinusoidales ------
% ------------------------------------------------------------------------

% Vector espacial.
z = linspace(0, d_1 + d_2 + 0.5, 2000);

% Cálculo de la envolvente final.
E_env = abs(E_i_0 .* exp(-gammas(1) .* z) + E_i_0 .* Gamma_T_1 .* exp(gammas(1) .* z - 2*gammas(1)*d_1)) .* ((z >= 0) & (z <= d_1));

% Configuración de la figura.
figure();

set(gcf, 'defaultLineLineWidth', 1);

med_1 = fill([0 d_1 d_1 0], [-2 -2 2 2], [0.95 0.95 0.95]);

hold on;

med_2 = fill([d_1 d_1+d_2 d_1+d_2 d_1], [-2 -2 2 2], [0.9 0.9 1]);
med_3 = fill([d_1+d_2 d_1+d_2+0.5 d_1+d_2+0.5 d_1+d_2], [-2 -2 2 2], [0.9 1 0.9]);

env_arriba = plot(z, E_env, ':', 'color', 'k', 'linewidth', 0.75);
env_abajo = plot(z, -E_env, ':', 'color', 'k', 'linewidth', 0.75);

E_T = plot(z, zeros(size(z)), 'b', 'linewidth', 2);
plot([d_1 d_1], [-2, 2], '--k');
plot([d_1+d_2 d_1+d_2], [-2, 2], '--k');

hold off;

xlim([0, d_1 + d_2 + 0.5]);
ylim([-2, 2]);

legend([med_1 med_2 med_3 E_T env_arriba], ...
       {"Medio 1", "Medio 2", "Medio 3", "Campo eléctrico", "Envolvente final"});

xlabel('z [m]', 'fontsize', 14);
ylabel('E [V/m]', 'fontsize', 14);

title("Simulación de ondas en los 3 medios (solamente E)", 'fontsize', 14);

% ------------------------------------------------------------------------

% Tiempos clave.
t_1 = d_1 / vels(1); % La incidente llegó a la primera interfaz.
t_2 = t_1 + d_2 / vels(2); % La onda llegó a la segunda interfaz.
t_3 = t_2 + d_2 / vels(2); % La onda volvió a atravesar la primera interfaz.

% Tiempo para la animación.
t_max = t_3 * 1.75;
N = 1500; % Número de frames.

% Vector temporal.
t = linspace(0, t_max, N);

% ------------------------------------------------------------------------

% Animación temporal por fases.
for k = 1:N
  t_eff = t(k);
  E_eff = zeros(size(z));

  % Incidente (siempre presente).
  z_front = vels(1)*t_eff;
  mask = z <= min(z_front, d_1);
  E_eff(mask) = E_i_0 .* exp(-alphas(1)*z(mask)) .* cos(omega*t_eff - betas(1)*z(mask));

  % Fase 2.
  if t_eff >= t_1
    % Reflejada parcial.
    z_back = vels(1)*(t_eff - t_1);
    mask = (z <= d_1) & (z >= d_1 - z_back);
    E_eff(mask) += E_i_0 * Gamma_12 * cos(omega*t_eff + betas(1)*(z(mask) - 2*d_1));

    % Transmitida al medio 2.
    z2_prog = vels(2)*(t_eff - t_1);
    mask = (z >= d_1) & (z <= d_1 + min(z2_prog, d_2));
    z_2 = z(mask) - d_1;
    E_eff(mask) += E_i_0 * tau_12 .* exp(-alphas(2)*z_2) .* cos(omega*t_eff - betas(2)*z_2);
  endif

  % Fase 3.
  if t_eff >= t_2
    % Transmitida al medio 3.
    z3_prog = vels(3)*(t_eff - t_2);
    mask = (z >= d_1 + d_2) & (z <= d_1 + d_2 + z3_prog);
    z_3 = z(mask) - (d_1 + d_2);
    E_eff(mask) += E_i_0 * tau_T_3 .* exp(-alphas(3)*z_3) .* cos(omega*t_eff - betas(3)*z_3);

    % Reflejada desde la segunda interfaz al medio 2.
    z2_back = vels(2)*(t_eff - t_2);
    mask = (z >= d_1 + d_2 - z2_back) & (z <= d_1 + d_2);
    z_2 = z(mask) - d_1;
    z_2_from_end = d_2 - z_2;
    E_eff(mask) += E_i_0 * Gamma_T_2 .* exp(-alphas(2)*z_2_from_end) .* cos(omega*t_eff + betas(2)*(z_2 - d_2));
  endif

  % Fase 4.
  if t_eff >= t_3
    % Reflejada final hacía el medio 1.
    z_back = vels(1)*(t_eff - t_3);
    mask = (z <= d_1) & (z >= d_1 - z_back);
    E_eff(mask) += E_i_0 * (Gamma_T_1 - Gamma_12) * cos(omega*t_eff + betas(1)*(z(mask) - 2*d_1));
  endif

  % Actualización.
  set(E_T, 'YData', E_eff);

  drawnow;
  pause(0.01);
endfor

pause(1);

% ------------------------------------------------------------------------
% -------------------- Simulación 2: pulso gaussiano ---------------------
% ------------------------------------------------------------------------

% La idea es cambiar el enfoque y hacer una simulación FDTD 1D (Yee)
% que directamente resuelve las ecuaciones de Maxwell en pasos espaciales.

dz = 3e-3; % Tamaño de paso espacial (resolución espacial).
Nc = round((d_1 + d_2 + 1)/dz); % Número de celdas.

dt = 0.99 * dz/c; % Tamaño de paso temporal (resolución temporal).
N = 1800; % Número de frames.

% Parámetros de la fuente.
t_0 = 40; % Tiempo central del pulso en frames.
ancho = 12; % Factor de ancho del pulso.
A = 1; % Amplitud del pulso.
loc = 2; % Ubicación de la fuente (primer celda de z != 0).

gauss = @(k) 2*A*exp(-((k - t_0)/ancho).^2);

% Vector espacial.
z = (0:Nc-1) * dz;

% Índices de las interfaces.
idx_d_1 = round(d_1/dz);
idx_d_2 = round((d_1 + d_2)/dz);

% Vectores de propiedades de los medios.
sigma_z = zeros(1, Nc); mu_r_z = zeros(1, Nc); eps_r_z = zeros(1, Nc);

sigma_z(1:idx_d_1) = sigma_1;
sigma_z(idx_d_1+1:idx_d_2) = sigma_2;
sigma_z(idx_d_2+1:end) = sigma_3;

mu_r_z(1:idx_d_1) = mu_r_1;
mu_r_z(idx_d_1+1:idx_d_2) = mu_r_2;
mu_r_z(idx_d_2+1:end) = mu_r_3;

eps_r_z(1:idx_d_1) = eps_r_1;
eps_r_z(idx_d_1+1:idx_d_2) = eps_r_2;
eps_r_z(idx_d_2+1:end) = eps_r_3;

% Campos y coeficientes de actualización.
E = zeros(1, Nc); CE_E = zeros(1, Nc); CE_H = zeros(1, Nc);
H = zeros(1, Nc);

% Seteo de los coeficientes de actualización.
CE_E = (1 - sigma_z .* dt ./ (2 .* eps_0 .* eps_r_z)) ./ ((1 + sigma_z .* dt ./ (2 .* eps_0 .* eps_r_z)));
CE_H = (dt ./ (eps_0 .* eps_r_z .* dz)) ./ (1 + sigma_z .* dt ./ (2 .* eps_0 .* eps_r_z));
CH_H = ones(1, Nc-1);
CH_E = dt ./ ((mu_0 .* mu_r_z)(1:end - 1) .* dz);

% Variables que guardan el estado inmediatamente anterior del campo eléctrico.
E_left_prev = 0; E_right_prev = 0;

% ------------------------------------------------------------------------

% Configuración de la figura.

med_1 = fill([0 d_1 d_1 0], [-2 -2 2 2], [0.95 0.95 0.95]);

hold on;

med_2 = fill([d_1 d_1+d_2 d_1+d_2 d_1], [-2 -2 2 2], [0.9 0.9 1]);
med_3 = fill([d_1+d_2 d_1+d_2+0.5 d_1+d_2+0.5 d_1+d_2], [-2 -2 2 2], [0.9 1 0.9]);

E_T = plot(z, zeros(size(z)), 'b', 'linewidth', 2);
plot([d_1 d_1], [-2, 2], '--k');
plot([d_1+d_2 d_1+d_2], [-2, 2], '--k');

hold off;

xlim([0, d_1 + d_2 + 0.5]);

legend([med_1 med_2 med_3 E_T], ...
       {"Medio 1", "Medio 2", "Medio 3", "Campo eléctrico"});

xlabel('z [m]', 'fontsize', 14);
ylabel('E [V/m]', 'fontsize', 14);

title("Simulación de pulso gaussiano en los 3 medios (solamente E)", 'fontsize', 14);

% ------------------------------------------------------------------------

threshold = 0.0035 * A; % Umbral de detección de ecos, calibrado.

% Variables de detección de ecos.
over_threshold_prev = false;
t_eco = 0;
eco_count = 0;

% Animación temporal.
for k = 1:N
  % Paso 1 - actualizar H.
  for i = 1:Nc-1
    H(i) = H(i) + CH_E(i) * (E(i+1) - E(i));
  endfor

  % Paso 2 - actualizar E.
  for i=2:Nc-1
    E(i) = CE_E(i) * E(i) + CE_H(i) * (H(i) - H(i-1));
  endfor

  % Paso 3 - sumar el efecto del pulso gaussiano.
  E(loc) = E(loc) + gauss(k);

  % Paso 4 - Aplicar las condiciones de contorno.
  E(1) = E_left_prev + ((c*dt - dz) / (c*dt + dz) * (E(2) - E(1)));
  E_left_prev = E(2);

  E(end) = E_right_prev + ((c*dt - dz) / (c*dt + dz) * (E(end-1) - E(end)));
  E_right_prev = E(end-1);

  over_threshold = abs(E(2)) > threshold;

  % Medimos solo hasta el segundo eco, que es el relevante.
  if over_threshold && ~over_threshold_prev && eco_count < 3
    t_eco = k * dt; % Guardar el tiempo del eco.
    eco_count++;
  endif

  over_threshold_prev = over_threshold;

  % Actualización.
  set(E_T, 'YData', E);
  drawnow;
endfor

pause(1);
clf;

% ------------------------------------------------------------------------
% --------------------- Mostrar resultados numéricos ---------------------
% ------------------------------------------------------------------------

% Calculamos la distancia hasta el medio 3 en el subsuelo en base al tiempo del eco.
d_2_est = (t_eco - 2 * d_1 / vels(1)) * vels(2) / 2;

% Calculamos el coeficiente de reflexión total Gamma_T_1.
E_rec = E_i_0 * (Gamma_T_1 - Gamma_12) .* cos(omega*t);
Gamma_T_1_est = ((E_rec ./ cos(omega*t)) / E_i_0 + Gamma_12)(1);

% Despejamos el coeficiente de reflexión Gamma_23.
Gamma_23_est = ((Gamma_T_1 - Gamma_12) * exp(2*gammas(2)*d_2_est)) / ...
               (tau_12*tau_21 + (Gamma_T_1 - Gamma_12)*Gamma_21);

% Despejamos la impedancia compleja del medio desconocido.
eta_est = etas(2) * (1 + Gamma_23_est) / (1 - Gamma_23_est);

% Despejamos la conductividad y permitividad dieléctrica del medio desconocido.
sigma_eps_rel = (1i * omega * mu_0) / (eta_est.^2);

sigma_3_est = abs(real(sigma_eps_rel));
eps_r_3_est = imag(sigma_eps_rel) / (omega * eps_0);

% Calculamos gamma.
gamma_3_est = sqrt((1i*omega*mu_0) .* (sigma_3_est + 1i*omega*eps_r_3_est*eps_0));

% ------------------------------------------------------------------------

% Presentamos todos los resultados en una tabla, junto con sus errores porcentuales.

% Función auxiliar para el error relativo porcentual.
function error = er(v_est, v_real)
  error = abs((v_est - v_real) / v_real) * 100;
endfunction

axis([0 1 0 1]);
axis off;

title(["Resultados (parámetros hallados numéricamente) - f=", num2str(f), " Hz"], 'fontsize', 14);

% Bordes.
line([0.00 0.85], [0.85 0.85], 'linewidth', 2);
line([0.00 0.85], [0.10 0.10], 'linewidth', 2);
line([0.00 0.00], [0.10 0.85], 'linewidth', 2);
line([0.85 0.85], [0.10 0.85], 'linewidth', 2);

% Texto.
text(0.03, 0.80, "Variable", 'fontsize', 14, 'fontweight', 'bold');
line([0.20 0.20], [0.10 0.85], 'linewidth', 2);
text(0.22, 0.80, "Valor estimado", 'fontsize', 14, 'fontweight', 'bold');
line([0.48 0.48], [0.10 0.85], 'linewidth', 2);
text(0.50, 0.80, "Valor real", 'fontsize', 14, 'fontweight', 'bold');
line([0.68 0.68], [0.10 0.85], 'linewidth', 2);
text(0.70, 0.80, "Error %", 'fontsize', 14, 'fontweight', 'bold');

line([0.00 0.85], [0.75 0.75], 'linewidth', 2);

% Variables y valores.

% t_eco.
text(0.02, 0.70, "t_{eco}", 'fontsize', 16); text(0.22, 0.70, num2str(t_eco), 'fontsize', 16);
text(0.50, 0.70, num2str(2*t_2), 'fontsize', 16); text(0.70, 0.70, [num2str(er(t_eco, 2*t_2)), " %"], 'fontsize', 16);

% d_2.
text(0.02, 0.62, "d_2", 'fontsize', 16); text(0.22, 0.62, num2str(d_2_est), 'fontsize', 16);
text(0.50, 0.62, num2str(d_2), 'fontsize', 16); text(0.70, 0.62, [num2str(er(d_2_est, d_2)), " %"], 'fontsize', 16);

% Gamma_23.
text(0.02, 0.54, "Γ_{23}", 'fontsize', 16); text(0.22, 0.54, num2str(Gamma_23_est), 'fontsize', 16);
text(0.50, 0.54, num2str(Gamma_23), 'fontsize', 16); text(0.70, 0.54, [num2str(er(Gamma_23_est, Gamma_23)), " %"], 'fontsize', 16);

% eta_3.
text(0.02, 0.46, "η_3", 'fontsize', 16); text(0.22, 0.46, num2str(eta_est), 'fontsize', 16);
text(0.50, 0.46, num2str(etas(3)), 'fontsize', 16); text(0.70, 0.46, [num2str(er(eta_est, etas(3))), " %"], 'fontsize', 16);

% sigma_3.
text(0.02, 0.38, "σ_3", 'fontsize', 16); text(0.22, 0.38, num2str(sigma_3_est), 'fontsize', 16);
text(0.50, 0.38, num2str(sigma_3), 'fontsize', 16); text(0.70, 0.38, [num2str(er(sigma_3_est, sigma_3)), " %"], 'fontsize', 16);

% mu_r_3.
text(0.02, 0.30, "μr_3", 'fontsize', 16); text(0.22, 0.30, num2str(1), 'fontsize', 16);
text(0.50, 0.30, num2str(mu_r_3), 'fontsize', 16); text(0.70, 0.30, [num2str(er(1, mu_r_3)), " %"], 'fontsize', 16);

% eps_r_3.
text(0.02, 0.22, "εr_3", 'fontsize', 16); text(0.22, 0.22, num2str(eps_r_3_est), 'fontsize', 16);
text(0.50, 0.22, num2str(eps_r_3), 'fontsize', 16); text(0.70, 0.22, [num2str(er(eps_r_3_est, eps_r_3)), " %"], 'fontsize', 16);

% gamma_3.
text(0.02, 0.14, "γ_3", 'fontsize', 16); text(0.22, 0.14, num2str(gamma_3_est), 'fontsize', 16);
text(0.50, 0.14, num2str(gammas(3)), 'fontsize', 16); text(0.70, 0.14, [num2str(er(gamma_3_est, gammas(3))), " %"], 'fontsize', 16);

% ------------------------------------------------------------------------
% --------------------------------- FIN ----------------------------------
% ------------------------------------------------------------------------
