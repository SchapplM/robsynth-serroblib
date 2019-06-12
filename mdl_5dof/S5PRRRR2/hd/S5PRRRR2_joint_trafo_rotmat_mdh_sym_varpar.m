% Calculate homogenous joint transformation matrices for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-03 15:11
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5PRRRR2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5PRRRR2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-03 15:11:16
% EndTime: 2019-06-03 15:11:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (16->8), ass. (0->9)
t29 = cos(qJ(2));
t28 = cos(qJ(3));
t27 = cos(qJ(4));
t26 = cos(qJ(5));
t25 = sin(qJ(2));
t24 = sin(qJ(3));
t23 = sin(qJ(4));
t22 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t29, -t25, 0, 0; 0, 0, -1, 0; t25, t29, 0, 0; 0, 0, 0, 1; t28, -t24, 0, 0; 0, 0, -1, 0; t24, t28, 0, 0; 0, 0, 0, 1; t27, -t23, 0, pkin(1); t23, t27, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t22, 0, 0; 0, 0, -1, 0; t22, t26, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
