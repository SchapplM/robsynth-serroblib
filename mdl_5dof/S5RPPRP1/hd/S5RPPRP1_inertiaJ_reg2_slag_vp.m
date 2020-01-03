% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t26 = sin(pkin(8));
t47 = -0.2e1 * t26;
t46 = 0.2e1 * t26;
t27 = sin(pkin(7));
t45 = t27 * pkin(1);
t29 = cos(pkin(7));
t44 = t29 * pkin(1);
t16 = qJ(3) + t45;
t30 = sin(qJ(4));
t43 = t16 * t30;
t22 = t26 ^ 2;
t31 = cos(qJ(4));
t42 = t22 * t31;
t41 = t30 * t26;
t28 = cos(pkin(8));
t18 = t30 * t28;
t20 = t31 * t26;
t40 = t31 * t28;
t23 = t28 ^ 2;
t39 = t22 + t23;
t24 = t30 ^ 2;
t25 = t31 ^ 2;
t14 = t24 + t25;
t38 = qJ(5) * t26;
t37 = t28 * t46;
t36 = t16 * t40;
t21 = -pkin(2) - t44;
t8 = -t28 * pkin(3) - t26 * pkin(6) + t21;
t5 = t31 * t8;
t35 = -t31 * t38 + t5;
t1 = (-pkin(4) - t43) * t28 + t35;
t2 = t36 + (t8 - t38) * t30;
t34 = t1 * t31 + t2 * t30;
t3 = -t16 * t18 + t5;
t4 = t30 * t8 + t36;
t33 = t3 * t31 + t4 * t30;
t19 = t25 * t22;
t17 = t24 * t22;
t15 = t16 ^ 2;
t13 = -0.2e1 * t30 * t42;
t12 = t40 * t47;
t11 = t30 * t37;
t10 = t22 * t15;
t9 = t14 * t26;
t7 = t19 + t17 + t23;
t6 = (pkin(4) * t30 + t16) * t26;
t32 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t45, 0, (t27 ^ 2 + t29 ^ 2) * pkin(1) ^ 2, t22, t37, 0, t23, 0, 0, -0.2e1 * t21 * t28, t21 * t46, 0.2e1 * t39 * t16, t23 * t15 + t21 ^ 2 + t10, t19, t13, t12, t17, t11, t23, 0.2e1 * t22 * t43 - 0.2e1 * t3 * t28, 0.2e1 * t16 * t42 + 0.2e1 * t4 * t28, t33 * t47, t3 ^ 2 + t4 ^ 2 + t10, t19, t13, t12, t17, t11, t23, -0.2e1 * t1 * t28 + 0.2e1 * t6 * t41, 0.2e1 * t2 * t28 + 0.2e1 * t6 * t20, t34 * t47, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t16 * t28 - t3 * t30 + t31 * t4) * t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t28 + (-t1 * t30 + t2 * t31) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t26, 0, t21, 0, 0, 0, 0, 0, 0, -t40, t18, -t9, t33, 0, 0, 0, 0, 0, 0, -t40, t18, -t9, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t41, -t28, t3, -t4, 0, 0, 0, 0, t20, 0, -t41, -t28, (-0.2e1 * pkin(4) - t43) * t28 + t35, -t2, -pkin(4) * t20, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t20, 0, -pkin(4) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, t31 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t20, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t32;
