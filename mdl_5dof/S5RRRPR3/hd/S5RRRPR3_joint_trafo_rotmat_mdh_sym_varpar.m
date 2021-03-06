% Calculate homogenous joint transformation matrices for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RRRPR3_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:44
% EndTime: 2020-01-03 12:08:44
% DurationCPUTime: 0.03s
% Computational Cost: add. (6->6), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t35 = cos(qJ(1));
t34 = cos(qJ(2));
t33 = cos(qJ(3));
t32 = cos(qJ(5));
t31 = sin(qJ(1));
t30 = sin(qJ(2));
t29 = sin(qJ(3));
t28 = sin(qJ(5));
t27 = cos(pkin(9));
t26 = sin(pkin(9));
t1 = [0, 0, 1, pkin(5); t31, t35, 0, 0; -t35, t31, 0, 0; 0, 0, 0, 1; t34, -t30, 0, pkin(1); t30, t34, 0, 0; 0, 0, 1, pkin(6); 0, 0, 0, 1; t33, -t29, 0, pkin(2); 0, 0, -1, -pkin(7); t29, t33, 0, 0; 0, 0, 0, 1; t27, -t26, 0, pkin(3); t26, t27, 0, 0; 0, 0, 1, qJ(4); 0, 0, 0, 1; t32, -t28, 0, pkin(4); t28, t32, 0, 0; 0, 0, 1, pkin(8); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
