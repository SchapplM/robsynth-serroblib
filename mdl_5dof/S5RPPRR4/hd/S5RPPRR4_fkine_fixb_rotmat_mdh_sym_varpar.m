% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:43:54
% EndTime: 2019-12-05 17:43:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (124->60), mult. (117->64), div. (0->0), fcn. (168->10), ass. (0->30)
t19 = cos(pkin(9));
t6 = t19 * pkin(3) + pkin(2);
t17 = sin(pkin(9));
t22 = sin(qJ(1));
t33 = t22 * t17;
t18 = sin(pkin(8));
t32 = t22 * t18;
t20 = cos(pkin(8));
t31 = t22 * t20;
t23 = cos(qJ(1));
t30 = t23 * t17;
t29 = t23 * t20;
t21 = -pkin(6) - qJ(3);
t16 = pkin(5) + 0;
t15 = pkin(9) + qJ(4);
t28 = t23 * qJ(2) + 0;
t27 = t23 * pkin(1) + t22 * qJ(2) + 0;
t8 = cos(t15);
t1 = pkin(4) * t8 + t6;
t14 = -pkin(7) + t21;
t26 = t1 * t20 - t14 * t18;
t25 = -t18 * t21 + t20 * t6;
t24 = pkin(2) * t20 + qJ(3) * t18;
t9 = qJ(5) + t15;
t7 = sin(t15);
t5 = cos(t9);
t4 = sin(t9);
t3 = t23 * t18;
t2 = t17 * pkin(3) + pkin(4) * t7;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t16; -t22, -t23, 0, 0; t23, -t22, 0, 0; 0, 0, 0, 1; t18, t20, 0, t16; -t31, t32, t23, -t22 * pkin(1) + t28; t29, -t3, t22, t27; 0, 0, 0, 1; t18 * t19, -t18 * t17, -t20, t18 * pkin(2) - t20 * qJ(3) + t16; -t19 * t31 + t30, t17 * t31 + t23 * t19, -t32, (-pkin(1) - t24) * t22 + t28; t19 * t29 + t33, -t17 * t29 + t22 * t19, t3, t24 * t23 + t27; 0, 0, 0, 1; t18 * t8, -t18 * t7, -t20, t18 * t6 + t20 * t21 + t16; t23 * t7 - t8 * t31, t23 * t8 + t7 * t31, -t32, pkin(3) * t30 + (-pkin(1) - t25) * t22 + t28; t22 * t7 + t8 * t29, t22 * t8 - t7 * t29, t3, pkin(3) * t33 + t25 * t23 + t27; 0, 0, 0, 1; t18 * t5, -t18 * t4, -t20, t18 * t1 + t20 * t14 + t16; t23 * t4 - t5 * t31, t23 * t5 + t4 * t31, -t32, t23 * t2 + (-pkin(1) - t26) * t22 + t28; t22 * t4 + t5 * t29, t22 * t5 - t4 * t29, t3, t22 * t2 + t26 * t23 + t27; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
