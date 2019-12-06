% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t41 = sin(pkin(9));
t42 = cos(pkin(9));
t45 = sin(qJ(2));
t48 = cos(qJ(2));
t28 = t41 * t45 - t42 * t48;
t30 = t41 * t48 + t42 * t45;
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t17 = t47 * t28 + t44 * t30;
t37 = -t48 * pkin(2) - pkin(1);
t23 = t28 * pkin(3) + t37;
t10 = t17 * pkin(4) + t23;
t61 = 0.2e1 * t10;
t60 = 0.2e1 * t23;
t59 = 0.2e1 * t37;
t58 = 0.2e1 * t48;
t57 = t41 * pkin(2);
t56 = t42 * pkin(2);
t43 = sin(qJ(5));
t55 = t43 * pkin(4);
t36 = pkin(3) + t56;
t26 = t44 * t36 + t47 * t57;
t46 = cos(qJ(5));
t54 = t46 * t26;
t53 = -qJ(3) - pkin(6);
t39 = t45 ^ 2;
t40 = t48 ^ 2;
t52 = t39 + t40;
t32 = t53 * t45;
t33 = t53 * t48;
t20 = t42 * t32 + t41 * t33;
t14 = -t30 * pkin(7) + t20;
t21 = t41 * t32 - t42 * t33;
t15 = -t28 * pkin(7) + t21;
t5 = t47 * t14 - t44 * t15;
t25 = t47 * t36 - t44 * t57;
t24 = pkin(4) + t25;
t12 = t46 * t24 - t43 * t26;
t6 = t44 * t14 + t47 * t15;
t38 = t46 * pkin(4);
t19 = -t44 * t28 + t47 * t30;
t13 = t43 * t24 + t54;
t9 = -t43 * t17 + t46 * t19;
t7 = t46 * t17 + t43 * t19;
t4 = -t17 * pkin(8) + t6;
t3 = -t19 * pkin(8) + t5;
t2 = t43 * t3 + t46 * t4;
t1 = t46 * t3 - t43 * t4;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t39, t45 * t58, 0, t40, 0, 0, pkin(1) * t58, -0.2e1 * pkin(1) * t45, 0.2e1 * t52 * pkin(6), t52 * pkin(6) ^ 2 + pkin(1) ^ 2, t30 ^ 2, -0.2e1 * t30 * t28, 0, t28 ^ 2, 0, 0, t28 * t59, t30 * t59, -0.2e1 * t20 * t30 - 0.2e1 * t21 * t28, t20 ^ 2 + t21 ^ 2 + t37 ^ 2, t19 ^ 2, -0.2e1 * t19 * t17, 0, t17 ^ 2, 0, 0, t17 * t60, t19 * t60, -0.2e1 * t6 * t17 - 0.2e1 * t5 * t19, t23 ^ 2 + t5 ^ 2 + t6 ^ 2, t9 ^ 2, -0.2e1 * t9 * t7, 0, t7 ^ 2, 0, 0, t7 * t61, t9 * t61, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t48, 0, -t45 * pkin(6), -t48 * pkin(6), 0, 0, 0, 0, t30, 0, -t28, 0, t20, -t21, (-t28 * t41 - t30 * t42) * pkin(2), (t20 * t42 + t21 * t41) * pkin(2), 0, 0, t19, 0, -t17, 0, t5, -t6, -t26 * t17 - t25 * t19, t5 * t25 + t6 * t26, 0, 0, t9, 0, -t7, 0, t1, -t2, -t12 * t9 - t13 * t7, t1 * t12 + t2 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t57, 0, (t41 ^ 2 + t42 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t26, 0, t25 ^ 2 + t26 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t12, -0.2e1 * t13, 0, t12 ^ 2 + t13 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t30, 0, t37, 0, 0, 0, 0, 0, 0, t17, t19, 0, t23, 0, 0, 0, 0, 0, 0, t7, t9, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t5, -t6, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, (-t43 * t7 - t46 * t9) * pkin(4), (t1 * t46 + t2 * t43) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t26, 0, 0, 0, 0, 0, 0, 0, 1, t12 + t38, -t54 + (-pkin(4) - t24) * t43, 0, (t12 * t46 + t13 * t43) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t55, 0, (t43 ^ 2 + t46 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t12, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t55, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
