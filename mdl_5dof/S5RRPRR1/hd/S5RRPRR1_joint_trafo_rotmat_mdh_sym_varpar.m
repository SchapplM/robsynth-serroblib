% Calculate homogenous joint transformation matrices for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RRPRR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:36
% EndTime: 2019-07-18 17:20:36
% DurationCPUTime: 0.03s
% Computational Cost: add. (5->5), mult. (0->0), div. (0->0), fcn. (16->8), ass. (0->9)
t36 = cos(qJ(1));
t35 = cos(qJ(2));
t34 = cos(qJ(4));
t33 = cos(qJ(5));
t32 = sin(qJ(1));
t31 = sin(qJ(2));
t30 = sin(qJ(4));
t29 = sin(qJ(5));
t1 = [t36, -t32, 0, 0; t32, t36, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t31, 0, 0; 0, 0, -1, 0; t31, t35, 0, 0; 0, 0, 0, 1; 1, 0, 0, pkin(1); 0, 1, 0, 0; 0, 0, 1, qJ(3); 0, 0, 0, 1; t34, -t30, 0, pkin(2); t30, t34, 0, 0; 0, 0, 1, pkin(3); 0, 0, 0, 1; t33, -t29, 0, 0; 0, 0, -1, -pkin(4); t29, t33, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
