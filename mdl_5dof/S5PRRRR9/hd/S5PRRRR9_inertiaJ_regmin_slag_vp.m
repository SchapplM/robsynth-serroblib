% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(5));
t34 = sin(qJ(4));
t37 = cos(qJ(5));
t38 = cos(qJ(4));
t19 = t33 * t34 - t37 * t38;
t35 = sin(qJ(3));
t15 = t19 * t35;
t60 = 0.2e1 * t15;
t26 = -t38 * pkin(4) - pkin(3);
t59 = 0.2e1 * t26;
t58 = -0.2e1 * t35;
t39 = cos(qJ(3));
t57 = 0.2e1 * t39;
t56 = pkin(8) + pkin(9);
t55 = pkin(3) * t38;
t54 = pkin(7) * t34;
t53 = t33 * pkin(4);
t52 = t37 * pkin(4);
t22 = -t39 * pkin(3) - t35 * pkin(8) - pkin(2);
t43 = t38 * t39;
t41 = pkin(7) * t43;
t9 = t41 + (-pkin(9) * t35 + t22) * t34;
t51 = t37 * t9;
t50 = t39 * pkin(4);
t31 = sin(pkin(5));
t49 = t31 * sin(qJ(2));
t48 = t31 * cos(qJ(2));
t47 = t34 * t35;
t46 = t34 * t38;
t45 = t34 * t39;
t44 = t38 * t35;
t42 = t35 * t57;
t18 = t38 * t22;
t6 = -pkin(9) * t44 + t18 + (-pkin(4) - t54) * t39;
t1 = -t33 * t9 + t37 * t6;
t20 = t33 * t38 + t37 * t34;
t32 = cos(pkin(5));
t30 = t39 ^ 2;
t29 = t38 ^ 2;
t28 = t35 ^ 2;
t27 = t34 ^ 2;
t24 = t56 * t38;
t23 = t56 * t34;
t21 = (pkin(4) * t34 + pkin(7)) * t35;
t17 = t32 * t35 + t39 * t49;
t16 = -t32 * t39 + t35 * t49;
t14 = t20 * t35;
t13 = t34 * t22 + t41;
t12 = -pkin(7) * t45 + t18;
t11 = -t33 * t23 + t37 * t24;
t10 = -t37 * t23 - t33 * t24;
t8 = t17 * t38 - t34 * t48;
t7 = -t17 * t34 - t38 * t48;
t4 = t33 * t7 + t37 * t8;
t3 = -t33 * t8 + t37 * t7;
t2 = t33 * t6 + t51;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t48, -t49, 0, 0, 0, 0, 0, t39 * t48, -t35 * t48, 0, 0, 0, 0, 0, t16 * t47 - t7 * t39, t16 * t44 + t8 * t39, 0, 0, 0, 0, 0, t16 * t14 - t3 * t39, -t16 * t15 + t4 * t39; 0, 1, 0, 0, t28, t42, 0, 0, 0, pkin(2) * t57, pkin(2) * t58, t29 * t28, -0.2e1 * t28 * t46, t43 * t58, t34 * t42, t30, -0.2e1 * t12 * t39 + 0.2e1 * t28 * t54, 0.2e1 * t28 * pkin(7) * t38 + 0.2e1 * t13 * t39, t15 ^ 2, t14 * t60, t39 * t60, t14 * t57, t30, -0.2e1 * t1 * t39 + 0.2e1 * t21 * t14, -0.2e1 * t21 * t15 + 0.2e1 * t2 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, -t16 * t38, t16 * t34, 0, 0, 0, 0, 0, t16 * t19, t16 * t20; 0, 0, 0, 0, 0, 0, t35, t39, 0, -t35 * pkin(7), -t39 * pkin(7), t34 * t44, (-t27 + t29) * t35, -t45, -t43, 0, -pkin(7) * t44 + (-pkin(3) * t35 + pkin(8) * t39) * t34, pkin(8) * t43 + (t54 - t55) * t35, -t15 * t20, -t20 * t14 + t15 * t19, -t20 * t39, t19 * t39, 0, -t10 * t39 + t26 * t14 + t21 * t19, t11 * t39 - t26 * t15 + t21 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t27, 0.2e1 * t46, 0, 0, 0, 0.2e1 * t55, -0.2e1 * pkin(3) * t34, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t59, t20 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t47, -t39, t12, -t13, 0, 0, -t15, -t14, -t39, -t37 * t50 + t1, -t51 + (-t6 + t50) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t38, 0, -t34 * pkin(8), -t38 * pkin(8), 0, 0, t20, -t19, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t14, -t39, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
