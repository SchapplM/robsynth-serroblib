% Script to fix outdated functions in the robot model database.
% The fkine functions were updated to also allow symbolic computations.
% Earlier versions of the code did not support this.
% This script updates the database for all kinematic structures

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2026-07
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
repopath=fileparts(which('serroblib_path_init.m'));

% Define symbolic variables to test for the bug
syms q1 q2 q3 q4 q5 q6 q7
q_7dof = [q1;q2;q3;q4;q5;q6;q7];
assume(q_7dof, 'real');

%% Loop all robot kinematics
for N = 1:7
  mdllistfile_Ndof = fullfile(repopath, sprintf('mdl_%ddof', N), ...
    sprintf('S%d_list.mat',N));
  l = load(mdllistfile_Ndof, 'Names_Ndof');
  for iFK = 1:length(l.Names_Ndof)
    Name = l.Names_Ndof{iFK};
    fprintf('Model %d/%d for %d DoF: %s\n', ...
      iFK, length(l.Names_Ndof), N, Name);
    % Initialize robot class
    RS = serroblib_create_robot_class(Name);
    RS.fill_fcn_handles(false, false);
    % if ~strcmp(Name, 'S3RRR1'), continue; end % Debug: Only single robot
    pkin_sym = sym(zeros(length(RS.pkin_names),1));
    RS.update_mdh(pkin_sym);
    q0 = q_7dof(1:RS.NJ); % define joint vector for testing
    %% Test the kinematics function to trigger the bug
    for i_try = 1:2
      clear T0E_test
      try
        T0E_test = RS.fkineEE(q0);
      catch err
        fprintf('Error calling fkineEE.\n');
        if strcmp(err.identifier, 'MATLAB:UnableToConvert')
          if i_try == 2
            error('failed on second try. Something went wrong');
          else
            fprintf('Regenerate code for kinematics.\n');
            % expected error for outdated version. Recreate code.
            serroblib_generate_mapleinput({Name});
            serroblib_generate_code({Name}, true, false, 4);
            % Notify Matlab about updated functions
            rehash; rehash path; pause(2);
            continue % try running the code again
          end
        else
          error('Unexpected error');
        end
      end % try-catch
      fprintf('Erfolgreich getestet.\n');
    end % for i_try
    serroblib_removefrompath({Name});
  end
end
