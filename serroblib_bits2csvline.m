% Wandle eine Roboterstruktur aus gegebenem Bit-Array in eine CSV-Zeile um
% Die csv-Zeile ist besser lesbar, das Bit-Array ist schneller
% Maschinen-Verarbeitbar
% 
% Eingabe:
% BA [1xN] uint16
%   Bit-Array zur Kennzeichnung aller MDH-Kinematikparameter.
%   Jede Spalte des Arrays (2Byte, uint16) entspricht einer Gelenk-Transfo.
%   Bits:
%   01 (LSB): Gelenktyp (1 Bit)
%   02-04:    beta (3 Bit)
%   05:       b (1 Bit)
%   06-08:    alpha (3 Bit)
%   09:       a (1 Bit)
%   10-12:    theta (3 Bit)
%   13:       d (1 Bit)
%   14-16:    offset (3 Bit)
% 
% Ausgabe:
% csvline
%   Cell-Array mit Spaltenweise Kinematik-Parameter (MDH) für alle Gelenke
% csvbits
%   Spaltenweise Indizes für die Daten in csvline
% 
% Siehe auch: serroblib_csvline2bits

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function [csvline, csvbits] = serroblib_bits2csvline(BA)

N = length(BA);
csvline = {'mdlname'}; % Erste Spalte: Platzhalter für Modellnamen
csvbits = [];
c = 1;
for kk = 1:N
  % In den 2 Byte sind alle Kinematik-Parameter enthalten mit jeweils 1 bis
  % 3 Bit. Die Bits kodieren den Index auf die möglichen Zustände.
  % Die Reihenfolge ist überall gleich
  b = 0; % Zähler für Bit-Verschiebung
  Bit_type   = bitand( bitshift( BA(kk), -b), bin2dec('1'));   b = b+1;
  Bit_beta   = bitand( bitshift( BA(kk), -b), bin2dec('111')); b = b+3;
  Bit_b      = bitand( bitshift( BA(kk), -b), bin2dec('1'));   b = b+1;
  Bit_alpha  = bitand( bitshift( BA(kk), -b), bin2dec('111')); b = b+3;
  Bit_a      = bitand( bitshift( BA(kk), -b), bin2dec('1'));   b = b+1;
  Bit_theta  = bitand( bitshift( BA(kk), -b), bin2dec('111')); b = b+3;
  Bit_d      = bitand( bitshift( BA(kk), -b), bin2dec('1'));   b = b+1;
  Bit_offset = bitand( bitshift( BA(kk), -b), bin2dec('111'));
  % Prüfen: `dec2bin(BA(kk),16)`

  % Index aus Bit-Vektor wieder herausholen
  % Das hier ist doppelt mit serroblib_add_robot.m.
  % Die Definitionen werden aber an den Gelenk-Index `kk` angepasst. Daher
  % keine zentrale Speicherung
  descr_type = {'R', 'P'};
  descr_beta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('beta%d',kk)};
  descr_b = {'0', sprintf('b%d',kk)};
  descr_alpha = {'0', 'pi/2', 'pi', '-pi/2', sprintf('alpha%d',kk)};
  descr_a = {'0', sprintf('a%d',kk)};
  descr_theta = {'0', 'pi/2', 'pi', '-pi/2', sprintf('theta%d',kk)};
  descr_d = {'0', sprintf('d%d',kk)};
  descr_offset = {'0', 'pi/2', 'pi', '-pi/2', sprintf('offset%d',kk)};
  
  % csv-Zeile erstellen
  % Siehe auch serroblib_add_robot.m
  c=c+1; csvline{c} = descr_type{Bit_type+1};     csvbits(c) = Bit_type+1;  %#ok<AGROW>
  c=c+1; csvline{c} = descr_beta{Bit_beta+1};     csvbits(c) = Bit_beta+1;  %#ok<AGROW>
  c=c+1; csvline{c} = descr_b{Bit_b+1};           csvbits(c) = Bit_b+1;     %#ok<AGROW>
  c=c+1; csvline{c} = descr_alpha{Bit_alpha+1};   csvbits(c) = Bit_alpha+1; %#ok<AGROW>
  c=c+1; csvline{c} = descr_a{Bit_a+1};           csvbits(c) = Bit_a+1;     %#ok<AGROW>
  c=c+1; csvline{c} = descr_theta{Bit_theta+1};   csvbits(c) = Bit_theta+1; %#ok<AGROW>
  if Bit_type == 0 % Drehgelenk
    csvline{c} = sprintf('q%d', kk);
  end
  c=c+1; csvline{c} = descr_d{Bit_d+1};           csvbits(c) = Bit_d+1;     %#ok<AGROW>
  if Bit_type == 1 % Schubgelenk
    csvline{c} = sprintf('q%d', kk);
  end
  c=c+1; csvline{c} = descr_offset{Bit_offset+1}; csvbits(c) = Bit_offset+1;%#ok<AGROW>
end