% Rufe alle Funktionen aller in der Bibliothek enthaltenen Roboter auf
% Dadurch wird sichergestellt, dass alles vorhanden ist und funktioniert.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-05
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc
repopath=fileparts(which('serroblib_path_init.m'));

%% Funktionen für alle Modelle kompilieren
for N = 1:7
  % Liste zusammenstellen
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof');
  
  % Debug: Modell finden:
  % II = find(strcmp(l.Names_Ndof, 'S6RRPRRR14'));
  % Alle Modelle nehmen
  II = 1:length(l.Names_Ndof);
  % Alle Funktionen ausführen
  for iFK = II
    %% Definition des Modells
    Name = l.Names_Ndof{iFK};
    fprintf('Modell %d/%d für %d FG: %s\n', iFK, length(l.Names_Ndof), N, Name);
    RS = serroblib_create_robot_class(Name);
    RS.gen_testsettings(true, true);
    q0 = rand(RS.NQJ,1);
    qD0 = rand(RS.NQJ,1);
    qDD0 = rand(RS.NQJ,1);
    T_E = RS.fkineEE(q0);
    xE = [T_E(1:3,4); r2eul(T_E(1:3,1:3), RS.phiconv_W_E)];

    %% Funktionen aufrufen
    for mex = [false, true]
      % Funktionen bei Bedarf kompilieren
      RS.fill_fcn_handles(mex, mex&true);
      % Debuggen: Bestimmte Funktionen löschen (falls sie alt sind)
      % delete(which(sprintf('%s_invkin_eulangresidual_mex', Name)));
      % RS.fill_fcn_handles(mex, mex&true);
      
      % Rufe Kinematik-Funktionen auf
      RS.fkine(q0);
      RS.fkineEE(q0);
      RS.jacobiR(q0);
      RS.jacobig(q0);
      RS.jacobit(q0);
      RS.jacobiw(q0);
      RS.jacobiwD(q0,qD0);
      RS.jacobia(q0);
      RS.jacobiaD(q0,qD0);
      RS.jtraf(q0);
      RS.jointvar(q0);

      % Rufe Dynamik-Funktionen mit Inertialparametern als Eingang auf
      RS.DynPar.mode = 2;
      RS.ekin(q0,qD0);
      RS.epot(q0);
      RS.gravload(q0);
      RS.inertia(q0);
      RS.corvec(q0,qD0);
      RS.cormat(q0,qD0);
      RS.invdyn(q0,qD0,qDD0);
      [w, wreg] = RS.internforce(q0,qD0,qDD0);
      f_test = reshape( wreg(1:3*RS.NL,:)*RS.DynPar.ipv_floatb, 3, RS.NL);
      f_test - w(1:3,:);
      m_test = reshape( wreg(3*RS.NL+1:end,:)*RS.DynPar.ipv_floatb, 3, RS.NL);
      m_test - w(4:6,:);
      RS.internforce_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RS.invdyn_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      
      % Rufe Dynamik-Funktionen mit Minimalparametervektor als Eingang auf
      RS.DynPar.mode = 4;
      [T,Treg]=RS.ekin(q0,qD0);
      [U,Ureg]=RS.epot(q0);
      [g,greg]=RS.gravload(q0);
      [M,Mreg]=RS.inertia(q0);
      [c,creg]=RS.corvec(q0,qD0);
      [C,Creg]=RS.cormat(q0,qD0);
      [t,treg]=RS.invdyn(q0,qD0,qDD0);
      RS.invdyn_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RS.invdyn2_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RV = RS.invdynregmat_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RS.invdyn3_traj(RV);

      % Rufe IK-Funktionen auf
      RS.constr1(q0, xE);
      RS.constr1grad(q0, xE);
      RS.constr2(q0, xE);
      RS.constr2grad(q0, xE);
      [q_ik1, Phi_ik1] = RS.invkin2(xE, q0+0.1*ones(RS.NJ,1), struct('reci', true));
      [q_ik2, Phi_ik2] = RS.invkin2(xE, q0+0.1*ones(RS.NJ,1), struct('reci', false));
    end
    % Roboter wieder aus Pfad entfernen
    serroblib_removefrompath({Name});
  end
end
