% Calculate homogenous joint transformation matrices for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RRRRR2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:28
% EndTime: 2019-03-29 15:25:28
% DurationCPUTime: 0.03s
% Computational Cost: add. (5->5), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t37 = cos(qJ(1));
t36 = cos(qJ(2));
t35 = cos(qJ(3));
t34 = cos(qJ(4));
t33 = cos(qJ(5));
t32 = sin(qJ(1));
t31 = sin(qJ(2));
t30 = sin(qJ(3));
t29 = sin(qJ(4));
t28 = sin(qJ(5));
t1 = [t37, -t32, 0, 0; t32, t37, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t36, -t31, 0, pkin(1); t31, t36, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t30, 0, 0; 0, 0, -1, 0; t30, t35, 0, 0; 0, 0, 0, 1; t34, -t29, 0, pkin(2); t29, t34, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t28, 0, 0; 0, 0, -1, 0; t28, t33, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
