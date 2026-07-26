% Calculate homogenous joint transformation matrices for
% S3PPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,theta1]';
% 
% Output:
% T_mdh [4x4x3]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(3+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2026-07-26 14:53
% Revision: 05cfe754a3b58d1958e6f60508538bf5b2996c81 (2026-06-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S3PPP1_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2026-07-26 14:53:29
% EndTime: 2026-07-26 14:53:29
% DurationCPUTime: 0.05s
% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
t9 = cos(pkin(3));
t8 = sin(pkin(3));
t1 = [t9, -t8, 0, 0; t8, t9, 0, 0; 0, 0, 1, qJ(1); 0, -1, 0, pkin(1); 0, 0, -1, -qJ(2); 1, 0, 0, 0; 1, 0, 0, pkin(2); 0, 0, -1, -qJ(3); 0, 1, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,3);             % numerisch
else,                         T_mdh = sym('xx', [4,4,3]); end % symbolisch

for i = 1:3
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
