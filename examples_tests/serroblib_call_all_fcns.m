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
      Tcges = RS.fkine(q0);
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
      T2 = RS.ekin(q0,qD0);
      U2 = RS.epot(q0);
      g2 = RS.gravload(q0);
      M2 = RS.inertia(q0);
      c2 = RS.corvec(q0,qD0);
      C2 = RS.cormat(q0,qD0);
      t2 = RS.invdyn(q0,qD0,qDD0);

      [w, wreg] = RS.internforce(q0,qD0,qDD0);
      f_test = reshape( wreg(1:3*RS.NL,:)*RS.DynPar.ipv_floatb, 3, RS.NL);
      test_f = f_test - w(1:3,:);
      if any(abs(test_f(:)) > 1e-10)
        error('Schnittkräfte stimmen nicht zwischen zwei Rechenwegen');
      end
      m_test = reshape( wreg(3*RS.NL+1:end,:)*RS.DynPar.ipv_floatb, 3, RS.NL);
      test_m = m_test - w(4:6,:);
      if any(abs(test_m(:)) > 1e-10)
        error('Schnittmomente stimmen nicht zwischen zwei Rechenwegen');
      end
      
      RS.internforce_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RS.invdyn_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      
      % Rufe Dynamik-Funktionen mit Inertialparametern als Eingang auf und
      % benutze die Regressorform
      RS.DynPar.mode = 3;
      [T3,Treg3]=RS.ekin(q0,qD0);
      [U3,Ureg3]=RS.epot(q0);
      [g3,greg3]=RS.gravload(q0);
      [M3,Mreg3]=RS.inertia(q0);
      [c3,creg3]=RS.corvec(q0,qD0);
      [C3,Creg3]=RS.cormat(q0,qD0);
      [t3,treg3]=RS.invdyn(q0,qD0,qDD0);
      if abs(T3-Treg3*RS.DynPar.ipv) > 1e-10 || abs(T3-T2) > 1e-10
        error('Regressor für kinetische Energie stimmt nicht mit DynPar Methode 3');
      end
      if abs(U3-Ureg3*RS.DynPar.ipv) > 1e-10 || abs(U3-U2) > 1e-10
        error('Regressor für potentielle Energie stimmt nicht mit DynPar Methode 3');
      end
      if any(abs(g3-greg3*RS.DynPar.ipv) > 1e-10) || abs(U3-U2) > 1e-10
        error('Regressor für Gravitationsmoment stimmt nicht mit DynPar Methode 3');
      end
      if any(abs(c3-creg3*RS.DynPar.ipv) > 1e-10) || any(abs(c3-c2) > 1e-10)
        error('Regressor für Coriolismoment stimmt nicht mit DynPar Methode 3');
      end
      if any(abs(t3-treg3*RS.DynPar.ipv)> 1e-10)  || any(abs(t3-t2) > 1e-10)
        error('Regressor für Inversdynamik-Moment stimmt nicht mit DynPar Methode 3');
      end
      
      % Rufe Dynamik-Funktionen mit Minimalparametervektor als Eingang auf
      RS.DynPar.mode = 4;
      [T4,Treg4]=RS.ekin(q0,qD0);
      [U4,Ureg4]=RS.epot(q0);
      [g4,greg4]=RS.gravload(q0);
      [M4,Mreg4]=RS.inertia(q0);
      [c4,creg4]=RS.corvec(q0,qD0);
      [C4,Creg4]=RS.cormat(q0,qD0);
      [t4,treg4]=RS.invdyn(q0,qD0,qDD0);
      if abs(T4-Treg4*RS.DynPar.mpv) > 1e-10 || abs(T4-T2) > 1e-10
        error('Regressor für kinetische Energie stimmt nicht mit DynPar Methode 4');
      end
      if abs(U4-Ureg4*RS.DynPar.mpv) > 1e-10 || abs(U4-U2) > 1e-10
        error('Regressor für potentielle Energie stimmt nicht mit DynPar Methode 4');
      end
      if any(abs(g4-greg4*RS.DynPar.mpv) > 1e-10) || abs(U4-U2) > 1e-10
        error('Regressor für Gravitationsmoment stimmt nicht mit DynPar Methode 4');
      end
      if any(abs(c4-creg4*RS.DynPar.mpv) > 1e-10) || any(abs(c4-c2) > 1e-10)
        error('Regressor für Coriolismoment stimmt nicht mit DynPar Methode 4');
      end
      if any(abs(t4-treg4*RS.DynPar.mpv)> 1e-10)  || any(abs(t4-t2) > 1e-10)
        error('Regressor für Inversdynamik-Moment stimmt nicht mit DynPar Methode 4');
      end
      
      RS.invdyn_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RS.invdyn2_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RV = RS.invdynregmat_traj(repmat(q0',5,1),repmat(qD0',5,1),repmat(qDD0',5,1));
      RS.invdyn3_traj(RV);

      % Rufe IK-Funktionen auf
      RS.constr1(q0, xE);
      RS.constr1grad(q0, xE);
      Phi_fromclass = RS.constr2(q0, xE, true);
      eval(sprintf('[Phi_fromfcn,Tcges_stack] = %s_constr2(q0, xE, RS.pkin, RS.T_N_E, RS.phiconv_W_E, RS.I_EElink, true);', RS.mdlname));
      if any(abs(Phi_fromclass-Phi_fromfcn)>1e-14)
        error('Funktion constr2 stimmt nicht zwischen Klassenmethode und eigener Funktion');
      end
      % Prüfe einige Ergebnisse
      for kk = 1:RS.NJ+1
        Tc_kk = Tcges(:,:,kk);
        Tc_kk_from_stack = [Tcges_stack((kk-1)*3+1:kk*3,1:4); [0 0 0 1]];
        test_Tc_kk = Tc_kk-Tc_kk_from_stack;
        if any(abs(test_Tc_kk(:))>1e-14)
          error('Ausgegebene direkte Kinematik aus constr2 stimmt nicht gegen direkte Berechnung. Max Fehler: %1.1e', max(abs(test_Tc_kk(:))));
        end
      end
      RS.constr2grad(q0, xE);
      [q_ik1, Phi_ik1, Tcges_stack_IK1] = RS.invkin2(xE, q0+0.1*ones(RS.NJ,1), struct('reci', true, 'Phir_tol', 1e-3, 'Phit_tol', 1e-3));
      [q_ik2, Phi_ik2] = RS.invkin2(xE, q0+0.1*ones(RS.NJ,1), struct('reci', false));
      Tcges_IK1 = RS.fkine(q_ik1);
      for kk = 1:RS.NJ+1
        Tc_kk = Tcges_IK1(:,:,kk);
        Tc_kk_from_stack = [Tcges_stack_IK1((kk-1)*3+1:kk*3,1:4); [0 0 0 1]];
        test_Tc_kk = Tc_kk-Tc_kk_from_stack;
        if any(abs(test_Tc_kk(:))>1e-14)
          error('Ausgegebene direkte Kinematik aus invkin2 stimmt nicht gegen direkte Berechnung. Max Fehler: %1.1e', max(abs(test_Tc_kk(:))));
        end
      end
    end
    % Roboter wieder aus Pfad entfernen
    serroblib_removefrompath({Name});
  end
end
