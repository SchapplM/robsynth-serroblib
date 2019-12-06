% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t30 = sin(pkin(10));
t32 = cos(pkin(10));
t35 = sin(qJ(4));
t59 = cos(qJ(4));
t20 = t59 * t30 + t35 * t32;
t66 = -0.2e1 * t20;
t33 = cos(pkin(5));
t31 = sin(pkin(5));
t36 = sin(qJ(2));
t58 = t31 * t36;
t13 = -t30 * t58 + t33 * t32;
t14 = t33 * t30 + t32 * t58;
t5 = -t59 * t13 + t35 * t14;
t65 = t5 ^ 2;
t52 = pkin(7) + qJ(3);
t21 = t52 * t32;
t47 = t52 * t30;
t9 = t35 * t21 + t59 * t47;
t64 = t9 ^ 2;
t18 = t35 * t30 - t59 * t32;
t63 = t18 ^ 2;
t24 = -t32 * pkin(3) - pkin(2);
t62 = 0.2e1 * t24;
t61 = 0.2e1 * t32;
t60 = t5 * t9;
t38 = cos(qJ(2));
t57 = t31 * t38;
t34 = sin(qJ(5));
t56 = t34 * t18;
t55 = t34 * t20;
t37 = cos(qJ(5));
t54 = t34 * t37;
t53 = t37 * t20;
t25 = t30 ^ 2;
t27 = t32 ^ 2;
t51 = t25 + t27;
t28 = t34 ^ 2;
t29 = t37 ^ 2;
t50 = t28 + t29;
t49 = t18 * t66;
t48 = t34 * t53;
t46 = -pkin(4) * t20 - pkin(8) * t18;
t11 = t59 * t21 - t35 * t47;
t8 = t18 * pkin(4) - t20 * pkin(8) + t24;
t1 = -t34 * t11 + t37 * t8;
t2 = t37 * t11 + t34 * t8;
t45 = t1 * t37 + t2 * t34;
t44 = -t1 * t34 + t2 * t37;
t7 = t35 * t13 + t59 * t14;
t3 = -t34 * t7 - t37 * t57;
t4 = -t34 * t57 + t37 * t7;
t43 = t3 * t37 + t4 * t34;
t42 = -t3 * t34 + t4 * t37;
t41 = -t13 * t30 + t14 * t32;
t26 = t31 ^ 2;
t22 = t26 * t38 ^ 2;
t16 = t20 ^ 2;
t15 = t37 * t18;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t36 ^ 2 + t33 ^ 2 + t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t14 ^ 2 + t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t22 + t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t58, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t57, -t30 * t57, t41, pkin(2) * t57 + t41 * qJ(3), 0, 0, 0, 0, 0, 0, -t18 * t57, -t20 * t57, -t7 * t18 + t5 * t20, t7 * t11 - t24 * t57 + t60, 0, 0, 0, 0, 0, 0, t3 * t18 + t5 * t55, -t4 * t18 + t5 * t53, -t43 * t20, t3 * t1 + t4 * t2 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t25, t30 * t61, 0, t27, 0, 0, pkin(2) * t61, -0.2e1 * pkin(2) * t30, 0.2e1 * t51 * qJ(3), t51 * qJ(3) ^ 2 + pkin(2) ^ 2, t16, t49, 0, t63, 0, 0, t18 * t62, t20 * t62, -0.2e1 * t11 * t18 + 0.2e1 * t9 * t20, t11 ^ 2 + t24 ^ 2 + t64, t29 * t16, -0.2e1 * t16 * t54, 0.2e1 * t18 * t53, t28 * t16, t34 * t49, t63, 0.2e1 * t1 * t18 + 0.2e1 * t9 * t55, -0.2e1 * t2 * t18 + 0.2e1 * t9 * t53, t45 * t66, t1 ^ 2 + t2 ^ 2 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t18, t20, 0, t24, 0, 0, 0, 0, 0, 0, t15, -t56, -t50 * t20, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t37, t5 * t34, t42, -t5 * pkin(4) + t42 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, -t9, -t11, 0, 0, t48, (-t28 + t29) * t20, t56, -t48, t15, 0, t46 * t34 - t9 * t37, t9 * t34 + t46 * t37, t44, -t9 * pkin(4) + t44 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, 0.2e1 * t54, 0, t29, 0, 0, 0.2e1 * pkin(4) * t37, -0.2e1 * pkin(4) * t34, 0.2e1 * t50 * pkin(8), t50 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t55, t18, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t37, 0, -t34 * pkin(8), -t37 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
