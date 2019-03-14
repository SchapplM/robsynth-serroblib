% Inverse Kinematik basierend auf Euler-Winkel-Residuum für
% S6RRRRRR1
%
% Eingabe:
% xE_soll
%   EE-Lage (Sollwert)
% q0
%   Anfangs-Gelenkwinkel für Algorithmus
% s
%   Struktur mit Eingabedaten. Felder, siehe Quelltext.
%
% Ausgabe:
% q
%   Lösung der IK
% Phi
%   Restfehler mit der IK-Lösung

% TODO: Nullraum-Optimierung ist noch nicht gut implementiert (z.B. Indizes)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-28 15:58
% Revision: 2bf3b907e1213de0593c9d1d0a7eb98ef6ddbfca (2019-02-28)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [q, Phi] = S6RRRRRR1_invkin_eulangresidual(xE_soll, q0, s)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),struct(
%$cgargs            'pkin', zeros(11,1),
%$cgargs          'sigmaJ', zeros(6,1),
%$cgargs             'NQJ', 0,
%$cgargs            'qlim', zeros(6,2),
%$cgargs            'I_EE', true(1,6),
%$cgargs     'phiconv_W_E', uint8(2),
%$cgargs        'I_EElink', uint8(0),
%$cgargs            'reci', true,
%$cgargs           'T_N_E', zeros(4,4),
%$cgargs        'task_red', false,
%$cgargs               'K', zeros(6,1),
%$cgargs              'Kn', zeros(6,1),
%$cgargs              'wn', 0,
%$cgargs           'n_min', 0,
%$cgargs           'n_max', 1000,
%$cgargs        'Phit_tol', 1.0000e-10,
%$cgargs        'Phir_tol', 1.0000e-10,
%$cgargs     'retry_limit', 100)}

%% Initialisierung

% Einstellungsvariablen aus Struktur herausholen
sigmaJ = s.sigmaJ; % Marker für Dreh-/Schubgelenk (in den Minimalkoordinaten)
phiconv_W_E = s.phiconv_W_E;
K = s.K;
Kn = s.Kn;
task_red = s.task_red;
n_min = s.n_min;
n_max = s.n_max;
wn = s.wn;
Phit_tol = s.Phit_tol;
Phir_tol = s.Phir_tol;
retry_limit = s.retry_limit;
reci = s.reci;
I_EElink = s.I_EElink;
T_N_E = s.T_N_E;
pkin = s.pkin;
qmin = s.qlim(:,1);
qmax = s.qlim(:,2);
n_Phi_t = sum(s.I_EE(1:3)); % Anzahl der translatorischen Zwangsbedingungen

if wn ~= 0
  nsoptim = true;
else
  % Keine zusätzlichen Optimierungskriterien
  nsoptim = false;
end
% Damit der Roboter einen Nullraum für Nebenoptimierungen hat, muss er min.
% 7FG für 6FG-Aufgaben und 6FG für 5FG-Aufgaben haben.
if nsoptim && s.NQJ < 6-task_red
  nsoptim = false;
end

% Indizes für Koordinaten festlegen
I_IK = 1:6;
if task_red
  xE_soll(6) = 0; % Dieser Wert hat keinen Einfluss auf die Berechnung, darf aber aufgrund der Implementierung nicht NaN sein.
  % Indizes zur Auswahl der berücksichtigten ZB-Komponenten
  I_IK = [1 2 3 5 6]; % Nehme an, dass immer die vierte Komponente des Fehlers weggelassen wird
elseif any(~s.I_EE)
  % Falls EE-FG nicht gefordert sind: Streiche die entsprechenden Zeilen in
  % den ZB und der Jacobi
  % TODO: Bessere Anpassung an Euler-Winkel
  % TODO: Funktioniert aktuell eigentlich nur für Methode 1
  I_IK = find(s.I_EE);
  task_red = true; % TODO: Eigener Marker hierfür
end
% Ausgabe belegen
Phi = NaN(length(I_IK),1);

success = false;
q1 = q0;
%% Iterative Berechnung der inversen Kinematik
for rr = 1:retry_limit % Schleife über Neu-Anfänge der Berechnung

  % Variablen zum Speichern der Zwischenergebnisse
  q1 = q0;
  for jj = 2:n_max % Schleife über iteratives Annähern mit Newton-Verfahren

    % Gradienten-Matrix
    Jdk_voll = S6RRRRRR1_constr2grad(q1, xE_soll, pkin, T_N_E, phiconv_W_E, I_EElink, reci);
    % Zwangsbedingungen
    Phi_voll = S6RRRRRR1_constr2(q1, xE_soll, pkin, T_N_E, phiconv_W_E, I_EElink, reci);

    %% Aufgabenredundanz
    if task_red
      % Aufgabenredundanz. Lasse den letzten Rotations-FG wegfallen
      Jdk = Jdk_voll(I_IK,:);
      Phi = Phi_voll(I_IK);
    else
      Jdk = Jdk_voll;
      Phi = Phi_voll;
    end
    %% Nullstellensuche für Positions- und Orientierungsfehler
    % (Optimierung der Aufgabe)
    % Normale Invertierung der Jacobi-Matrix der seriellen Kette
    delta_q_T = Jdk \ (-Phi);
    %% Optimierung der Nebenbedingungen (Nullraum)
    delta_q_N = zeros(size(delta_q_T));
    if nsoptim && jj < n_max-10 % die letzten Iterationen sind zum Ausgleich des Positionsfehlers (ohne Nullraum)
      % Berechne Gradienten der zusätzlichen Optimierungskriterien
      v = zeros(s.NQJ, 1);
      if wn(1) ~= 0
        [~, hdq] = invkin_optimcrit_limits1(q1, [qmin, qmax]);
        % [1], Gl. (25)
        v = v - hdq';
      end
      % [1], Gl. (24)
      delta_q_N = (eye(s.NQJ) - pinv(Jdk)* Jdk) * v;
    end

    % [1], Gl. (23)
    delta_q = K.*delta_q_T + Kn.*delta_q_N;
    q2 = q1 + delta_q;

    if any(isnan(q2)) || any(isinf(q2))
      break; % ab hier kann das Ergebnis nicht mehr besser werden wegen NaN/Inf
    end

    % Nachverarbeitung der Ergebnisse der Iteration
    q1 = q2;
    q1(sigmaJ==0) = normalize_angle(q1(sigmaJ==0)); % nur Winkel normalisieren

    % Abbruchbedingungen prüfen
    if jj > n_min ... % Mindestzahl Iterationen erfüllt
        && all(abs(Phi(1:n_Phi_t)) < Phit_tol) && all(abs(Phi(n_Phi_t+1:end)) < Phir_tol) && ... % Haupt-Bedingung ist erfüllt
        ( ~nsoptim || ... %  und keine Nebenoptimierung läuft
        nsoptim && all(abs(delta_q_N) < 1e-10) ) % oder die Nullraumoptimierung läuft noch
      success = true;
      break;
    end
  end
  if success
    break;
  end
  % Beim vorherigen Durchlauf kein Erfolg. Generiere neue Anfangswerte
  q0 = qmin + rand(s.NQJ,1).*(qmax-qmin);
end
q = q1;