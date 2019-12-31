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
%   02-04: beta
%   ... (siehe serroblib_bits2csvline)
% 
% BAR [1x1] uint16
%   Bit-Array zur Kennzeichnung der Rotation vom KS N zum KS E
%   In X-Y-Z-Euler-Winkeln (Die Bits kodieren diskrete Zustände der
%   einzelnen Winkel; 0,pi/2,pi,...; siehe Quelltext)
%   Bits:
%   01 (LSB) - 03: Erster Euler-Winkel
%   04-06: Zweiter Euler-Winkel
%   07-09: Dritter Euler-Winkel
% 
% BAE [1x1] uint16
%   Bit-Array zur Kennzeichnung der EE-FG
%   Bits:
%   01 (LSB): vx0
%   ...
%   06: wz0
%   07: phix0
%   08: phiy0
%   09: phiz0
%
% BAJVF [1xN] uint16
%   Bit-Array zum Filtern von Varianten eines Roboters
%   Dadurch kann bestimmt werden, welcher Roboter bereits durch ein
%   allgemeineres Robotermodell beschrieben werden kann (durch allgemeine
%   DH-Parameter, die beim spezifischen Modell auf einen Wert wie 0, pi/2
%   gesetzt werden)
%   Bits entsprechen denen aus BAJ, sind auf 0 gesetzt, wenn Parameter/Bits
%   nicht betrachtet werden müssen
% 
% BAO [1x1] uint16
%   Bit-Array zur Kennzeichnung der Modellherkunft
%   Bits:
%   01 (LSB): Spalte "Manuell" (Kette wurde vom Benutzer erstellt)
%   02: Spalte "Roboter" (Kette kommt aus Struktursynthese mit EE-FG)
%   03: Spalte "3T0R-PKM" (Beinkette für PKM mit EE-FG)
%   04: Spalte "3T1R-PKM"
% 
% Siehe auch: serroblib_bits2csvline

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [BAJ, BAR, BAE, BAJVF, BAO] = serroblib_csvline2bits(csvline)

%% Initialisierung
% Prüfe Eingabe
% Nicht Teil der Gelenk-Einträge: Name (1 Spalte), EE-Trafo (3), EE-FG (6
% Spalten), EE-FG Euler-Winkel (3), Nummer des Pos.-beeinfl. Gelenks (1),
% Gelenktypen (U/S/...) (1), Herkunft der Kinematik (5 Spalten)
N1 = (length(csvline)-1-3-6-3-2-5);
if mod(N1,8) ~= 0
  error('falsche Anzahl Einträge in csvline');
end
N = N1/8; % 8 Spalten pro Gelenk. So herausfinden der Gelenkzahl

% Ausgabevariable vorbelegen
BAJ = uint16(zeros(1,N));
BAJVF = uint16(Inf(1,N));
BAR = uint16(0);
BAO = uint16(0);
%% Bit-Vektor für Gelenk-Parameter aus csv-Zeile gewinnen
c = 1;
for kk = 1:N % über alle Gelenk-FG
  % Inhalt der csv-Zeilen. Siehe auch serroblib_bits2csvline.m und
  % serroblib_add_robot.m
  descr_type = {'R', 'P'};
  descr_beta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('beta%d',kk)};
  descr_b = {'0', sprintf('b%d',kk)};
  descr_theta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('theta%d',kk), 'q%d'};
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
  c=c+1;
  if Bit_type == 1
    % Schubgelenk: theta ist Parameter
    Bit_theta  = uint16(find(strcmp(csvline{c},descr_theta ))-1);
  else
    % Drehgelenk: theta ist Gelenkkoordinate. Setze auf 0
    Bit_theta = uint16(0);
  end
  c=c+1; 
  if Bit_type == 0
    % Drehgelenk, d ist Parameter
    Bit_d      = uint16(find(strcmp(csvline{c},descr_d     ))-1);
  else
    Bit_d = uint16(0);
  end
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
  
  % Bit-Array für Filter genauso zusammenstellen
  bits_disable = [];
  if Bit_beta == 4,  bits_disable = [bits_disable, 2:4]; end %#ok<AGROW>
  if Bit_b == 1,     bits_disable = [bits_disable, 5]; end %#ok<AGROW>
  if Bit_alpha == 4, bits_disable = [bits_disable, 6:8]; end %#ok<AGROW>
  if Bit_a == 1,     bits_disable = [bits_disable, 9]; end %#ok<AGROW>
  if Bit_theta == 4, bits_disable = [bits_disable, 10:12]; end %#ok<AGROW>
  if Bit_d == 1,     bits_disable = [bits_disable, 13]; end %#ok<AGROW>
  for i = bits_disable
    BAJVF(kk) = bitset(BAJVF(kk), i, 0);
  end
end
% Prüfen mit: `dec2bin(BAJ, 16)`, `dec2bin(BAJVF, 16)`
%% Bit-Vektor für EE-Transformation von letztem Körper
% Mögliche diskrete Zustände für die phi-Winkel in der csv-Tabelle. 
% "?" steht für einen noch undefinierten Winkel
descr_phi = {'0', 'pi/2', 'pi', '-pi/2', '?'};
b = 0; % Bit-Offset zur Verschiebung der Parameter-Bits in der Gesamtvariable
for kk = 1:3
  c=c+1; Bit_phi   = uint16(find( strcmp(csvline{c},descr_phi))-1);
  BAR = bitor( BAR, bitshift(Bit_phi,b)); b = b+3;
	% Prüfen mit: `dec2bin(Bit_phi)`
end
% Prüfen mit: `dec2bin(BAR)`
%% Bit-Vektor für EE-FG aus csv-Zeile gewinnen
BAE = serroblib_csvline2bits_EE(csvline(c+1:c+9));
c = c+9;
%% Weitere Spalten für zusätzliche Eigenschaften
% Spalte für Nummer des Positionsbeeinflussenden Gelenks überspringen
c = c+1;
% Spalte für Gelenkfolge (z.B. "UPS") überspringen
c = c+1;
%% Bit-Vektor für Modellherkunft
b = 0; % Bit-Offset zur Verschiebung der Parameter-Bits in der Gesamtvariable
for kk = 1:5 % 5 Tabellenspalten -> 5 Bits
  c=c+1; Bit_phi = uint16( strcmp(csvline{c},'1') );
  BAO = bitor( BAO, bitshift(Bit_phi,b)); b = b+1;
  % Prüfen mit: `dec2bin(Bit_phi)`
end
return
% Prüfen mit: `dec2bin(BAO)`
