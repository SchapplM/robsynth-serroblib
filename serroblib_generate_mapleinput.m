% Generiere die Eingabedateien zur Code-Generierung für gegebene
% Roboterstrukturen
% 
% Eingabe:
% Names
%   Cell-Array mit Liste der Namen der Roboterstrukturen, für die der Code
%   erzeugt werden soll. Der Name entspricht dem Schema "SxRRRyyy" mit
%   x=Anzahl Gelenk-FG und yyy laufende Nummer für "RRR".

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function serroblib_generate_mapleinput(Names)

%#ok<*AGROW>
for i = 1:length(Names)
  n = Names{i};
  N = str2double(n(2));
  
  % Datenbank laden
  repopath=fileparts(which('serroblib_path_init.m'));
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof', 'BitArrays_Ndof');
  
  
  % Index in Gesamt-Tabelle laden
  index_ges = find(strcmp(l.Names_Ndof, n));
  
  % MDH-Tabelle laden
  BA = l.BitArrays_Ndof(index_ges,:); %#ok<FNDSB>
  [~, csvbits] = serroblib_bits2csvline(BA);
  
  % Maple-Toolbox-Eingabe erzeugen
  mapleinputfile=fullfile(repopath, sprintf('mdl_%ddof', N), n, 'hd', sprintf('robot_env_%s', n));
  mkdirs(fileparts(mapleinputfile));
  fid = fopen(mapleinputfile, 'w');
  % Allgemeine Definitionen
  fprintf(fid, sprintf('robot_name := "%s":\n', n));
  fprintf(fid, sprintf('NJ := %d:\n', N));
  fprintf(fid, sprintf('NQJ := %d:\n', N));

  % Zeilen für alle MDH-Parameter
  line_v = 'v := Matrix(NJ,1,[';
  line_mu = 'mu := Matrix(NJ,1,[';
  line_sigma = 'sigma := Matrix(NJ,1,[';
  
  line_beta = 'beta := Matrix(NJ,1,[';
  line_b = 'b := Matrix(NJ,1,[';
  line_alpha = 'alpha := Matrix(NJ,1,[';
  line_a = 'a := Matrix(NJ,1,[';
  line_theta = 'theta := Matrix(NJ,1,[';
  line_d = 'd := Matrix(NJ,1,[';
  line_offset = 'qoffset := Matrix(NJ,1,[';
  
  % Hänge die einzelnen MDH-Parameter für alle Gelenk-Transformationen an
  % die Eingabevektoren an
  for kk = 1:N
    % Namensformat für die Maple-Dynamik-Toolbox
    % Die Indizes entsprechen der Reihenfolge für die csv-Zeile
    maple_type = {'0', '1'};
    maple_beta = {'0', 'Pi/2', 'Pi', '-Pi/2', sprintf('beta%d',kk)};
    maple_b = {'0', sprintf('b%d',kk)};
    maple_alpha = {'0', 'Pi/2', 'Pi', '-Pi/2', sprintf('alpha%d',kk)};
    maple_a = {'0', sprintf('a%d',kk)};
    maple_theta = {'0', 'Pi/2', 'Pi', '-Pi/2', sprintf('theta%d',kk)};
    maple_d = {'0', sprintf('d%d',kk)};
    maple_offset = {'0', 'Pi/2', 'Pi', '-Pi/2', sprintf('qoffset%d',kk)};
    
    
    line_v     = [line_v,     sprintf('%d,', kk-1)];
    line_mu    = [line_mu,    sprintf('%d,', 1)];
    % Benutze den Index aus csvbits zur Adressierung der Eigenschaften
    % aller MDH-Parameter
    line_sigma = [line_sigma, maple_type{  csvbits(2+8*(kk-1)) }, ',' ];
    line_beta  = [line_beta,  maple_beta{  csvbits(3+8*(kk-1)) }, ',' ];
    line_b     = [line_b,     maple_b{     csvbits(4+8*(kk-1)) }, ',' ];
    line_alpha = [line_alpha, maple_alpha{ csvbits(5+8*(kk-1)) }, ',' ];
    line_a     = [line_a,     maple_a{     csvbits(6+8*(kk-1)) }, ',' ];
    
    if csvbits(2+8*(kk-1)) == 1
      % Drehgelenk, nehme die Gelenkvariable
      line_theta = [line_theta, sprintf('qJ_t(%d,1) + %s', kk, maple_offset{ csvbits(9+8*(kk-1)) }), ',' ];
    else
      % Schubgelenk, nehme den konstanten Winkel theta
      line_theta = [line_theta, maple_theta{ csvbits(7+8*(kk-1)) }, ',' ];
    end
    if csvbits(2+8*(kk-1)) == 2
      % Schubgelenk, nehme die Gelenkvariable
      line_d = [line_d, sprintf('qJ_t(%d,1) + %s', kk, maple_offset{ csvbits(9+8*(kk-1)) }), ',' ];
    else
      % Drehgelenk, nehme die konstante Verschiebung d
      line_d = [line_d, maple_d{ csvbits(8+8*(kk-1)) }, ',' ];
    end
    
    line_offset= [line_offset, maple_offset{ csvbits(9+8*(kk-1)) }, ',' ];
  end
  
  % Eckige Klammer für Liste schließen, Zeilen beenden
  line_v      = [     line_v(1:end-1), ']):'];
  line_mu     = [    line_mu(1:end-1), ']):'];
  line_sigma  = [ line_sigma(1:end-1), ']):'];
  line_beta   = [  line_beta(1:end-1), ']):'];
  line_b      = [     line_b(1:end-1), ']):'];
  line_alpha  = [ line_alpha(1:end-1), ']):'];
  line_a      = [     line_a(1:end-1), ']):'];
  line_theta  = [ line_theta(1:end-1), ']):'];
  line_d      = [     line_d(1:end-1), ']):'];
  line_offset = [line_offset(1:end-1), ']):'];
  
  fwrite(fid, [line_v,      newline]);
  fwrite(fid, [line_mu,     newline]);
  fwrite(fid, [line_sigma,  newline]);
  fwrite(fid, [line_beta,   newline]);
  fwrite(fid, [line_b,      newline]);
  fwrite(fid, [line_alpha,  newline]);
  fwrite(fid, [line_a,      newline]);
  fwrite(fid, [line_theta,  newline]);
  fwrite(fid, [line_d,      newline]);
  fwrite(fid, [line_offset, newline]);
  
  fclose(fid);
end
