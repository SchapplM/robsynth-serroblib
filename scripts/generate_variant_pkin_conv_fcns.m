% Generiere Matlab-Funktionen, um die Kinematik-Parameter von Varianten
% automatisch aus dem Allgemeinen Modell abzuleiten
% 
% Gehe die gesamte Datenbank durch und erzeuge die Funktionen für alle
% Modelle
% 
% Erzeugt Dateien: SRR...PRy_pkin_gen2var.m, SRR...PRy_pkin_var2gen.m
% im Ordner /mdl_xdof/SRR...PRy/var
% 
% Mögliche Modifikation: Überschreiben aller bestehender Dateien durch
% Setzen von `only_add_new` auf `false`

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-04
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Initialisierung
roblibpath=fileparts(which('serroblib_path_init.m'));
only_add_new = true;

%% Alle Modelle durchgehen
for N = 1:7
  fprintf('Prüfe Strukturen mit %d FG\n', N);
  % Alle Roboter aus Datenbank laden
  mdllistfile_Ndof = fullfile(roblibpath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof', 'BitArrays_EEdof0', 'AdditionalInfo');
  for j = 1:length(l.Names_Ndof)
    Name = l.Names_Ndof{j};
    isvariant = l.AdditionalInfo(j,2);
    if ~isvariant
      continue % Nur Varianten durchgehen
    end
    fprintf('[%d/%d] Generiere Funktionen für %s\n', j, length(l.Names_Ndof), Name);
    variantof = l.AdditionalInfo(j,3);
    Name_GenMdl = l.Names_Ndof{variantof};
    
    % Kinematik-Parameter auslesen:
    % ... für allgemeines Modell
    [~, PSG] = serroblib_create_robot_class(Name_GenMdl, [], true);
    % ... für Variante
    [~, PSV] = serroblib_create_robot_class(Name, [], true);
    
    % Unterschiede in den MDH-Parametern feststellen und daraus die
    % Funktionen erstellen
    serroblib_addtopath({Name_GenMdl});
    mdh2pkin_hdl = eval(sprintf('@%s_mdhparam2pkin', Name_GenMdl));
    % Berechnung des Kinematik-Parametervektors des allgemeinen Modells mit
    % den Parametern der Variante. Allgemein gelassene Parameter sind NaN
    pkinGV = mdh2pkin_hdl(PSV.beta, PSV.b, PSV.alpha, PSV.a, PSV.theta, PSV.d, PSV.offset);
    % Indizes für Parameter, die bei der Variante von der allgemeinen Form
    % übernommen werden können
    I_pkinV = find(isnan(pkinGV)); % Alle NaN sind weiterhin allgemeine Parameter
    % Erstelle daraus die Funktion pkin_gen2var
    var_dir = fullfile(roblibpath, sprintf('mdl_%ddof', N), Name_GenMdl, 'var');
    mkdirs(var_dir);
    filename_gen2var = fullfile(var_dir, sprintf('%s_pkin_gen2var.m', Name));
    if ~only_add_new || only_add_new && ~exist(filename_gen2var, 'file')
      fid = fopen(filename_gen2var, 'w');
      fprintf(fid, '%% Umwandlung der Kinematikparameter von %s zu %s\n', Name_GenMdl, Name);
      fprintf(fid, '%% Eingabe:\n');
      fprintf(fid, '%% pkin_gen (%dx1) double\n', length(pkinGV));
      fprintf(fid, '%%   Kinematikparameter (pkin) von %s\n', Name_GenMdl);
      fprintf(fid, '%% Ausgabe:\n');
      fprintf(fid, '%% pkin_var (%dx1) double\n', length(I_pkinV));
      fprintf(fid, '%%   Kinematikparameter (pkin) von %s\n', Name);
      fprintf(fid, '%% I_gv (%dx1)\n', length(I_pkinV));
      fprintf(fid, '%%   Vektor mit Indizes zur Selektion von Kinematikparametern\n');
      fprintf(fid, 'function [pkin_var, I_gv] = %s_pkin_gen2var(pkin_gen)\n', Name);
      fprintf(fid, 'I_gv = [%s];\n', disp_array(I_pkinV','%d'));
      fprintf(fid, 'pkin_var = pkin_gen(I_gv);\n');
      fclose(fid);
    else
      fprintf('Datei %s existiert schon.\n', filename_gen2var);
    end
    % Zahlenwerte für Kinematik-Parameter, die festgelegt werden
    I_pkinG_fix = find(~isnan(pkinGV)); % Parameter ungleich NaN sind auf Zahlenwert der Variante festgelegt
    values_fix = pkinGV(I_pkinG_fix);
    % Erstelle daraus Funktion pkin_var2gen
    filename_var2gen = fullfile(var_dir, sprintf('%s_pkin_var2gen.m', Name));
    if ~only_add_new || only_add_new && ~exist(filename_var2gen, 'file')
      fid = fopen(filename_var2gen, 'w');
      fprintf(fid, '%% Umwandlung der Kinematikparameter von %s zu %s\n', Name, Name_GenMdl);
      fprintf(fid, '%% Eingabe:\n');
      fprintf(fid, '%% pkin_var (%dx1) double\n', length(I_pkinV));
      fprintf(fid, '%%   Kinematikparameter (pkin) von %s\n', Name);
      fprintf(fid, '%% Ausgabe:\n');
      fprintf(fid, '%% pkin_gen (%dx1) double\n', length(pkinGV));
      fprintf(fid, '%%   Kinematikparameter (pkin) von %s\n', Name_GenMdl);
      fprintf(fid, 'function pkin_gen = %s_pkin_var2gen(pkin_var)\n', Name);
      fprintf(fid, 'pkin_gen = zeros(%d,1);\n', length(pkinGV));
      fprintf(fid, 'pkin_gen([%s]) = pkin_var;\n', disp_array(I_pkinV','%d'));
      for ii = 1:length(I_pkinG_fix)
        if     abs(values_fix(ii)-pi/2)<1e-10, value_str = 'pi/2'; 
        elseif abs(values_fix(ii)-(-pi/2))<1e-10, value_str = '-pi/2'; 
        elseif abs(values_fix(ii)-pi)<1e-10, value_str = 'pi'; 
        else,  value_str = sprintf('%1.1f', values_fix(ii));
        end
        fprintf(fid, 'pkin_gen(%d) = %s;\n', I_pkinG_fix(ii), value_str);
      end
      fclose(fid);
    else
      fprintf('Datei %s existiert schon.\n', filename_var2gen);
    end
  end
end