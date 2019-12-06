% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t36 = sin(pkin(9));
t60 = -0.2e1 * t36;
t37 = cos(pkin(9));
t59 = 0.2e1 * t37;
t38 = sin(qJ(5));
t24 = t38 * t37;
t39 = sin(qJ(3));
t57 = t39 * pkin(2);
t27 = qJ(4) + t57;
t40 = cos(qJ(5));
t12 = -t37 * pkin(4) - t36 * pkin(7) - pkin(3);
t41 = cos(qJ(3));
t56 = t41 * pkin(2);
t7 = t12 - t56;
t2 = -t27 * t24 + t40 * t7;
t46 = qJ(4) * t37;
t5 = t40 * t12 - t38 * t46;
t58 = -t2 - t5;
t28 = -pkin(3) - t56;
t55 = pkin(3) - t28;
t50 = t40 * t37;
t3 = t27 * t50 + t38 * t7;
t32 = t36 ^ 2;
t52 = t32 * t40;
t54 = t27 * t52 + t3 * t37;
t30 = t32 * qJ(4);
t6 = t38 * t12 + t40 * t46;
t53 = t40 * t30 + t6 * t37;
t16 = t32 * t27;
t51 = t38 * t36;
t33 = t37 ^ 2;
t49 = t33 * t27 + t16;
t31 = t33 * qJ(4);
t48 = t31 + t30;
t34 = t38 ^ 2;
t35 = t40 ^ 2;
t47 = t34 + t35;
t21 = t36 * t59;
t45 = t2 * t40 + t3 * t38;
t44 = t6 * t38 + t5 * t40;
t42 = qJ(4) ^ 2;
t29 = t32 * t42;
t26 = t40 * t36;
t25 = t35 * t32;
t23 = t34 * t32;
t22 = t27 ^ 2;
t19 = t38 * t30;
t18 = -0.2e1 * t38 * t52;
t15 = t32 * t22;
t14 = t50 * t60;
t13 = t38 * t21;
t11 = t27 * t30;
t9 = t38 * t16;
t8 = t47 * t36;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 + t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 + t23 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t2 * t38 - t27 * t37 + t3 * t40) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t57, 0, (t39 ^ 2 + t41 ^ 2) * pkin(2) ^ 2, t32, t21, 0, t33, 0, 0, -0.2e1 * t28 * t37, 0.2e1 * t28 * t36, 0.2e1 * t49, t33 * t22 + t28 ^ 2 + t15, t25, t18, t14, t23, t13, t33, -0.2e1 * t2 * t37 + 0.2e1 * t9, 0.2e1 * t54, t45 * t60, t2 ^ 2 + t3 ^ 2 + t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t38 * t5 + t40 * t6 - t46) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t57, 0, 0, t32, t21, 0, t33, 0, 0, t55 * t37, -t55 * t36, t48 + t49, -t28 * pkin(3) + t27 * t31 + t11, t25, t18, t14, t23, t13, t33, t58 * t37 + t19 + t9, t53 + t54, (t58 * t40 + (-t3 - t6) * t38) * t36, t2 * t5 + t3 * t6 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t32, t21, 0, t33, 0, 0, pkin(3) * t59, pkin(3) * t60, 0.2e1 * t48, pkin(3) ^ 2 + t33 * t42 + t29, t25, t18, t14, t23, t13, t33, -0.2e1 * t5 * t37 + 0.2e1 * t19, 0.2e1 * t53, t44 * t60, t5 ^ 2 + t6 ^ 2 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, 0, t28, 0, 0, 0, 0, 0, 0, -t50, t24, -t8, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t50, t24, -t8, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t51, -t37, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t51, -t37, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
