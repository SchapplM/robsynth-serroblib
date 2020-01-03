% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:15
% EndTime: 2020-01-03 11:41:15
% DurationCPUTime: 0.19s
% Computational Cost: add. (124->61), mult. (117->64), div. (0->0), fcn. (168->10), ass. (0->34)
t18 = sin(qJ(3));
t27 = -t18 * pkin(3) - qJ(2);
t20 = cos(qJ(3));
t6 = t20 * pkin(3) + pkin(2);
t16 = cos(pkin(8));
t19 = sin(qJ(1));
t36 = t19 * t16;
t35 = t19 * t18;
t34 = t19 * t20;
t15 = sin(pkin(8));
t21 = cos(qJ(1));
t33 = t21 * t15;
t32 = t21 * t16;
t31 = t21 * t18;
t30 = t21 * t20;
t13 = qJ(3) + pkin(9);
t7 = sin(t13);
t29 = -pkin(4) * t7 + t27;
t17 = -qJ(4) - pkin(6);
t14 = pkin(5) + 0;
t28 = t19 * pkin(1) + 0;
t26 = -t19 * qJ(2) + 0;
t25 = pkin(2) * t16 + pkin(6) * t15;
t8 = cos(t13);
t1 = pkin(4) * t8 + t6;
t12 = -pkin(7) + t17;
t24 = t1 * t16 - t12 * t15;
t23 = -t15 * t17 + t16 * t6;
t22 = -t21 * qJ(2) + t28;
t9 = qJ(5) + t13;
t5 = cos(t9);
t4 = sin(t9);
t3 = t19 * t15;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t14; t19, t21, 0, 0; -t21, t19, 0, 0; 0, 0, 0, 1; t15, t16, 0, t14; t36, -t3, -t21, t22; -t32, t33, -t19, -t21 * pkin(1) + t26; 0, 0, 0, 1; t15 * t20, -t15 * t18, -t16, t15 * pkin(2) - t16 * pkin(6) + t14; t16 * t34 - t31, -t16 * t35 - t30, t3, t25 * t19 + t22; -t16 * t30 - t35, t16 * t31 - t34, -t33, (-pkin(1) - t25) * t21 + t26; 0, 0, 0, 1; t15 * t8, -t15 * t7, -t16, t15 * t6 + t16 * t17 + t14; -t21 * t7 + t8 * t36, -t21 * t8 - t7 * t36, t3, t19 * t23 + t21 * t27 + t28; -t19 * t7 - t8 * t32, -t19 * t8 + t7 * t32, -t33, 0 + t27 * t19 + (-pkin(1) - t23) * t21; 0, 0, 0, 1; t15 * t5, -t15 * t4, -t16, t15 * t1 + t16 * t12 + t14; -t21 * t4 + t5 * t36, -t21 * t5 - t4 * t36, t3, t19 * t24 + t21 * t29 + t28; -t19 * t4 - t5 * t32, -t19 * t5 + t4 * t32, -t33, 0 + t29 * t19 + (-pkin(1) - t24) * t21; 0, 0, 0, 1;];
T_ges = t2;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
