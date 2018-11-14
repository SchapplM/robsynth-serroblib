% Wandle eine Roboterstruktur aus gegebener CSV-Zeile in Bit-Array um
% Die csv-Zeile ist besser lesbar, das Bit-Array ist schneller
% Maschinen-Verarbeitbar
% 
% Eingabe:
% csvline
%   Cell-Array mit Spaltenweise Kinematik-Parameter (MDH) für alle Gelenke
%   csv-Format (aber nur die ersten Spalten, die noch Gelenk-Bezug haben.
%   Keine zusätzlichen Spalten, z.B. für EE-FG)
% 
% Ausgabe:
% BAJ [1xN] uint16
%   Bit-Array zur Kennzeichnung aller MDH-Kinematikparameter.
%   Jede Spalte des Arrays (2Byte, uint16) entspricht einer Gelenk-Transfo.
%   Bits:
%   01 (LSB):    Gelenktyp
%   02-05: beta
%   06:    b
%   07-10: alpha
%   ...
% 
% BAE [1x1] uint16
%   Bit-Array zur Kennzeichnung der EE-FG
%   Bits:
%   01 (LSB): vx0
%   ...
%   06: wz0
% 
% Siehe auch: serroblib_bits2csvline

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function [BAJ, BAE] = serroblib_csvline2bits(csvline)

%% Initialisierung
% Prüfe Eingabe
% Nicht Teil der Gelenk-Einträge: Name (1 Spalte), EE-FG (6 Spalten)
N1 = (length(csvline)-1-6);
if mod(N1,8) ~= 0
  error('falsche Anzahl Einträge in csvline');
end
N = N1/8; % 8 Spalten pro Gelenk. So herausfinden der Gelenkzahl

% Ausgabevariable vorbelegen
BAJ = uint16(zeros(1,N));
BAE = uint16(0);
%% Bit-Vektor für Gelenk-Parameter aus csv-Zeile gewinnen
c = 1;
for kk = 1:N % über alle Gelenk-FG
  % Inhalt der csv-Zeilen. Siehe auch serroblib_bits2csvline.m und
  % serroblib_add_robot.m
  descr_type = {'R', 'P'};
  descr_beta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('beta%d',kk)};
  descr_b = {'0', sprintf('b%d',kk)};
  descr_theta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('theta%d',kk)};
  descr_d = {'0', sprintf('d%d',kk)};
  descr_alpha = {'0', 'pi/2', 'pi', '-pi/2', sprintf('alpha%d',kk)};
  descr_a = {'0', sprintf('a%d',kk)};
  descr_offset = {'0', 'pi/2', 'pi', '-pi/2', sprintf('offset%d',kk)};

  % Finde die Nummer des Eintrages in den möglichen Optionen für die
  % csv-Zeile. Index ist Null-Basiert. Aus diesem Index (Bit_...) wird dann
  % das Bit-Array erzeugt.
  c=c+1; Bit_type   = uint16(find(strcmp(csvline{c},descr_type  ))-1);
  c=c+1; Bit_beta   = uint16(find(strcmp(csvline{c},descr_beta  ))-1);
  c=c+1; Bit_b      = uint16(find(strcmp(csvline{c},descr_b     ))-1);
  c=c+1; Bit_alpha  = uint16(find(strcmp(csvline{c},descr_alpha ))-1); 
  c=c+1; Bit_a      = uint16(find(strcmp(csvline{c},descr_a     ))-1);
  c=c+1; Bit_theta  = uint16(find(strcmp(csvline{c},descr_theta ))-1);
  c=c+1; Bit_d      = uint16(find(strcmp(csvline{c},descr_d     ))-1);
  c=c+1; Bit_offset = uint16(find(strcmp(csvline{c},descr_offset))-1);
  
  % Bit-Array aus den Bits für alle Parameter zusammenstellen
  b = 0; % Bit-Offset zur Verschiebung der Parameter-Bits in der Gesamtvariable
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_type,0)); b = b+1;
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_beta,b)); b = b+3;
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_b,b)); b = b+1;
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_alpha,b)); b = b+3;
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_a,b));b = b+1;
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_theta,b)); b = b+3;
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_d,b));b = b+1;
  BAJ(kk) = bitor( BAJ(kk), bitshift(Bit_offset,b));
end

%% Bit-Vektor für EE-FG aus csv-Zeile gewinnen
BAE = serroblib_csvline2bits_EE(csvline(c+1:end));

