% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:04:26
% EndTime: 2019-05-05 15:04:29
% DurationCPUTime: 1.21s
% Computational Cost: add. (807->98), mult. (1401->157), div. (0->0), fcn. (1608->6), ass. (0->76)
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t56 = cos(qJ(4));
t73 = t56 * t52;
t78 = sin(qJ(4));
t34 = -t78 * t51 + t73;
t87 = -0.2e1 * t34;
t68 = t78 * t52;
t32 = t56 * t51 + t68;
t29 = t32 ^ 2;
t30 = t34 ^ 2;
t86 = t29 + t30;
t53 = -pkin(1) - qJ(3);
t79 = -pkin(7) + t53;
t36 = t79 * t51;
t14 = t78 * t36 - t79 * t73;
t85 = t14 ^ 2;
t42 = t51 * pkin(3) + qJ(2);
t84 = 0.2e1 * t42;
t83 = 0.2e1 * qJ(2);
t82 = t32 * pkin(5);
t81 = t34 * pkin(4);
t55 = cos(qJ(5));
t80 = t55 * pkin(5);
t77 = t14 * t34;
t54 = sin(qJ(5));
t76 = t54 * t32;
t27 = t54 * t34;
t75 = t54 * t55;
t16 = t56 * t36 + t79 * t68;
t74 = t55 * t16;
t26 = t55 * t34;
t72 = -qJ(6) - pkin(8);
t46 = t51 ^ 2;
t47 = t52 ^ 2;
t39 = t46 + t47;
t49 = t54 ^ 2;
t50 = t55 ^ 2;
t40 = t49 + t50;
t71 = qJ(6) * t34;
t70 = t32 * t87;
t10 = t32 * pkin(4) - t34 * pkin(8) + t42;
t3 = t55 * t10 - t54 * t16;
t13 = t40 * t32;
t67 = -pkin(8) * t32 - t81;
t59 = -t55 * t71 + t3;
t1 = t59 + t82;
t2 = t74 + (t10 - t71) * t54;
t66 = t1 * t55 + t2 * t54;
t65 = -t1 * t54 + t2 * t55;
t4 = t54 * t10 + t74;
t64 = t3 * t55 + t4 * t54;
t63 = -t3 * t54 + t4 * t55;
t62 = -t16 * t32 + t77;
t37 = t72 * t54;
t38 = t72 * t55;
t61 = t55 * t37 - t54 * t38;
t60 = -t37 * t54 - t38 * t55;
t57 = qJ(2) ^ 2;
t44 = -pkin(4) - t80;
t41 = 0.2e1 * t75;
t31 = t39 * t53;
t25 = t55 * t32;
t24 = t50 * t30;
t21 = t49 * t30;
t20 = t54 * t26;
t19 = -0.2e1 * t30 * t75;
t18 = 0.2e1 * t32 * t26;
t17 = t54 * t70;
t12 = t40 * t34;
t11 = (-t49 + t50) * t34;
t8 = t86 * t54;
t7 = t40 * t29 + t30;
t6 = pkin(5) * t27 + t14;
t5 = t86 * t55;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t83 (pkin(1) ^ 2) + t57, t47, -0.2e1 * t52 * t51, 0, t46, 0, 0, t51 * t83, t52 * t83, -0.2e1 * t31, t39 * t53 ^ 2 + t57, t30, t70, 0, t29, 0, 0, t32 * t84, t34 * t84, 0.2e1 * t62, t16 ^ 2 + t42 ^ 2 + t85, t24, t19, t18, t21, t17, t29, 0.2e1 * t14 * t27 + 0.2e1 * t3 * t32, 0.2e1 * t14 * t26 - 0.2e1 * t4 * t32, t64 * t87, t3 ^ 2 + t4 ^ 2 + t85, t24, t19, t18, t21, t17, t29, 0.2e1 * t1 * t32 + 0.2e1 * t27 * t6, -0.2e1 * t2 * t32 + 0.2e1 * t26 * t6, t66 * t87, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t39, t31, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t62, 0, 0, 0, 0, 0, 0, -t8, -t5, 0, t32 * t63 - t77, 0, 0, 0, 0, 0, 0, -t8, -t5, 0, t32 * t65 - t6 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t52, 0, qJ(2), 0, 0, 0, 0, 0, 0, t32, t34, 0, t42, 0, 0, 0, 0, 0, 0, t25, -t76, -t12, t64, 0, 0, 0, 0, 0, 0, t25, -t76, -t12, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t32, 0, -t14, -t16, 0, 0, t20, t11, t76, -t20, t25, 0, -t14 * t55 + t54 * t67, t14 * t54 + t55 * t67, t63, -t14 * pkin(4) + pkin(8) * t63, t20, t11, t76, -t20, t25, 0, t27 * t44 + t37 * t32 - t6 * t55, t26 * t44 + t38 * t32 + t6 * t54, -t34 * t61 + t65, t1 * t37 - t2 * t38 + t6 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t27, t13, pkin(8) * t13 + t81, 0, 0, 0, 0, 0, 0, t26, -t27, t13, t32 * t60 - t34 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t49, t41, 0, t50, 0, 0, 0.2e1 * pkin(4) * t55, -0.2e1 * pkin(4) * t54, 0.2e1 * t40 * pkin(8), pkin(8) ^ 2 * t40 + pkin(4) ^ 2, t49, t41, 0, t50, 0, 0, -0.2e1 * t44 * t55, 0.2e1 * t44 * t54, 0.2e1 * t60, t37 ^ 2 + t38 ^ 2 + t44 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t27, t32, t3, -t4, 0, 0, 0, 0, t26, 0, -t27, t32, t59 + 0.2e1 * t82, -t2, -pkin(5) * t26, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t25, 0, -pkin(5) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t55, 0, -t54 * pkin(8), -t55 * pkin(8), 0, 0, 0, 0, t54, 0, t55, 0, t37, t38, -t54 * pkin(5), t37 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t54, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t9;
