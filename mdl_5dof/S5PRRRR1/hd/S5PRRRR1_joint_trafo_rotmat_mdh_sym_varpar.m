% Calculate homogenous joint transformation matrices for
% S5PRRRR1
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5PRRRR1_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:16
% EndTime: 2019-07-18 13:28:16
% DurationCPUTime: 0.02s
% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (16->8), ass. (0->9)
t30 = cos(qJ(2));
t29 = cos(qJ(3));
t28 = cos(qJ(4));
t27 = cos(qJ(5));
t26 = sin(qJ(2));
t25 = sin(qJ(3));
t24 = sin(qJ(4));
t23 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t30, -t26, 0, pkin(1); 0, 0, -1, 0; t26, t30, 0, 0; 0, 0, 0, 1; t29, -t25, 0, 0; 0, 0, -1, 0; t25, t29, 0, 0; 0, 0, 0, 1; t28, -t24, 0, pkin(2); t24, t28, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t23, 0, 0; 0, 0, -1, 0; t23, t27, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
