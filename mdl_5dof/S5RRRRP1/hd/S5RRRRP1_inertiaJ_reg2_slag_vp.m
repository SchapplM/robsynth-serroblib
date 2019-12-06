% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t47 = cos(qJ(3));
t48 = cos(qJ(2));
t28 = -t44 * t45 + t47 * t48;
t29 = t44 * t48 + t47 * t45;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t14 = -t46 * t28 + t43 * t29;
t36 = -t48 * pkin(2) - pkin(1);
t21 = -t28 * pkin(3) + t36;
t6 = t14 * pkin(4) + t21;
t62 = 0.2e1 * t6;
t61 = 0.2e1 * t21;
t60 = 0.2e1 * t29;
t59 = 0.2e1 * t48;
t58 = -pkin(7) - pkin(6);
t57 = t43 * pkin(3);
t56 = t44 * pkin(2);
t39 = t46 * pkin(3);
t41 = t45 ^ 2;
t42 = t48 ^ 2;
t55 = t41 + t42;
t54 = t46 * t56;
t31 = t58 * t45;
t32 = t58 * t48;
t18 = t44 * t31 - t47 * t32;
t10 = t28 * pkin(8) + t18;
t17 = t47 * t31 + t44 * t32;
t9 = -t29 * pkin(8) + t17;
t3 = -t43 * t10 + t46 * t9;
t40 = t47 * pkin(2);
t35 = t40 + pkin(3);
t25 = t46 * t35 - t43 * t56;
t4 = t46 * t10 + t43 * t9;
t49 = 0.2e1 * pkin(4);
t53 = t25 + t49;
t51 = pkin(3) ^ 2;
t38 = t43 ^ 2 * t51;
t37 = -0.2e1 * t57;
t34 = t39 + pkin(4);
t26 = t43 * t35 + t54;
t24 = t26 ^ 2;
t23 = pkin(4) + t25;
t22 = 0.2e1 * t26;
t20 = t26 * t57;
t19 = -t54 + (-pkin(3) - t35) * t43;
t16 = t43 * t28 + t46 * t29;
t13 = t16 ^ 2;
t12 = t14 ^ 2;
t11 = t14 * t57;
t7 = t26 * t14;
t5 = -0.2e1 * t16 * t14;
t2 = -t14 * qJ(5) + t4;
t1 = -t16 * qJ(5) + t3;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t41, t45 * t59, 0, t42, 0, 0, pkin(1) * t59, -0.2e1 * pkin(1) * t45, 0.2e1 * t55 * pkin(6), t55 * pkin(6) ^ 2 + pkin(1) ^ 2, t29 ^ 2, t28 * t60, 0, t28 ^ 2, 0, 0, -0.2e1 * t36 * t28, t36 * t60, -0.2e1 * t17 * t29 + 0.2e1 * t18 * t28, t17 ^ 2 + t18 ^ 2 + t36 ^ 2, t13, t5, 0, t12, 0, 0, t14 * t61, t16 * t61, -0.2e1 * t4 * t14 - 0.2e1 * t3 * t16, t21 ^ 2 + t3 ^ 2 + t4 ^ 2, t13, t5, 0, t12, 0, 0, t14 * t62, t16 * t62, -0.2e1 * t1 * t16 - 0.2e1 * t2 * t14, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t48, 0, -t45 * pkin(6), -t48 * pkin(6), 0, 0, 0, 0, t29, 0, t28, 0, t17, -t18, (t28 * t44 - t29 * t47) * pkin(2), (t17 * t47 + t18 * t44) * pkin(2), 0, 0, t16, 0, -t14, 0, t3, -t4, -t25 * t16 - t7, t3 * t25 + t4 * t26, 0, 0, t16, 0, -t14, 0, t1, -t2, -t23 * t16 - t7, t1 * t23 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t56, 0, (t44 ^ 2 + t47 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -t22, 0, t25 ^ 2 + t24, 0, 0, 0, 0, 0, 1, 0.2e1 * t23, -t22, 0, t23 ^ 2 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t28, 0, t17, -t18, 0, 0, 0, 0, t16, 0, -t14, 0, t3, -t4, -t16 * t39 - t11, (t3 * t46 + t4 * t43) * pkin(3), 0, 0, t16, 0, -t14, 0, t1, -t2, -t34 * t16 - t11, t1 * t34 + t2 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t56, 0, 0, 0, 0, 0, 0, 0, 1, t25 + t39, t19, 0, t25 * t39 + t20, 0, 0, 0, 0, 0, 1, t39 + t53, t19, 0, t23 * t34 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t39, t37, 0, t46 ^ 2 * t51 + t38, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, t37, 0, t34 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t3, -t4, 0, 0, 0, 0, t16, 0, -t14, 0, t1, -t2, -t16 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t26, 0, 0, 0, 0, 0, 0, 0, 1, t53, -t26, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t39, -t57, 0, 0, 0, 0, 0, 0, 0, 1, t49 + t39, -t57, 0, t34 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t49, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
