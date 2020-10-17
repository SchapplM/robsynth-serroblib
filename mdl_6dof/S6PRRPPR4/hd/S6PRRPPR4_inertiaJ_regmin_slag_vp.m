% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:12:32
% EndTime: 2019-05-05 03:12:34
% DurationCPUTime: 0.63s
% Computational Cost: add. (405->107), mult. (942->190), div. (0->0), fcn. (1081->10), ass. (0->72)
t52 = sin(pkin(11));
t54 = cos(pkin(11));
t89 = t52 ^ 2 + t54 ^ 2;
t55 = sin(qJ(6));
t58 = cos(qJ(6));
t31 = t58 * t52 - t55 * t54;
t56 = sin(qJ(3));
t23 = t31 * t56;
t88 = 0.2e1 * t23;
t68 = t52 * qJ(5) + pkin(3);
t83 = pkin(4) + pkin(5);
t28 = t83 * t54 + t68;
t87 = 0.2e1 * t28;
t33 = -t54 * pkin(4) - t68;
t86 = -0.2e1 * t33;
t85 = 0.2e1 * t56;
t59 = cos(qJ(3));
t84 = 0.2e1 * t59;
t82 = pkin(3) * t52;
t81 = pkin(8) * t54;
t47 = t56 * pkin(8);
t80 = t59 * pkin(8);
t70 = cos(pkin(6));
t53 = sin(pkin(6));
t77 = t53 * sin(qJ(2));
t27 = t70 * t56 + t59 * t77;
t76 = t53 * cos(qJ(2));
t13 = t27 * t54 - t52 * t76;
t9 = t13 * t54;
t26 = t56 * t77 - t70 * t59;
t79 = t26 * t52;
t78 = t26 * t54;
t41 = t52 * t56;
t42 = t54 * t56;
t34 = -t59 * pkin(3) - t56 * qJ(4) - pkin(2);
t21 = t52 * t34 + t54 * t80;
t73 = t89 * qJ(4) ^ 2;
t72 = qJ(4) * t59;
t71 = qJ(5) * t54;
t45 = t52 * qJ(4);
t11 = t27 * t52 + t54 * t76;
t69 = t11 * t52 + t9;
t67 = t11 * t59 + t26 * t41;
t39 = t52 * t80;
t20 = t54 * t34 - t39;
t66 = t11 ^ 2 + t13 ^ 2 + t26 ^ 2;
t65 = qJ(4) * t9 + t11 * t45;
t17 = -t59 * qJ(5) + t21;
t48 = t59 * pkin(4);
t18 = -t20 + t48;
t64 = t17 * t54 + t18 * t52;
t63 = -t20 * t52 + t21 * t54;
t30 = t55 * t52 + t58 * t54;
t62 = (t11 * t54 - t13 * t52) * t56;
t51 = t56 ^ 2;
t38 = t59 * t45;
t36 = (-pkin(9) + qJ(4)) * t54;
t35 = -t52 * pkin(9) + t45;
t32 = 0.2e1 * t89 * qJ(4);
t24 = t30 * t56;
t22 = t47 + (pkin(4) * t52 - t71) * t56;
t16 = -t47 + (-t83 * t52 + t71) * t56;
t15 = t55 * t35 + t58 * t36;
t14 = t58 * t35 - t55 * t36;
t8 = pkin(9) * t41 + t17;
t6 = t59 * pkin(5) + t39 + t48 + (-pkin(9) * t56 - t34) * t54;
t5 = t13 * t59 + t26 * t42;
t4 = t11 * t55 + t13 * t58;
t3 = t11 * t58 - t13 * t55;
t2 = t55 * t6 + t58 * t8;
t1 = -t55 * t8 + t58 * t6;
t7 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0; 0, 0, t76, -t77, 0, 0, 0, 0, 0, t59 * t76, -t56 * t76, t67, t5, t62, -t11 * t20 + t13 * t21 + t26 * t47, t67, t62, -t5, t11 * t18 + t13 * t17 + t26 * t22, 0, 0, 0, 0, 0, t26 * t23 + t3 * t59, -t26 * t24 - t4 * t59; 0, 1, 0, 0, t51, t56 * t84, 0, 0, 0, pkin(2) * t84, -0.2e1 * pkin(2) * t56, 0.2e1 * t51 * pkin(8) * t52 - 0.2e1 * t20 * t59, 0.2e1 * t21 * t59 + 0.2e1 * t51 * t81 (-t20 * t54 - t21 * t52) * t85, t51 * pkin(8) ^ 2 + t20 ^ 2 + t21 ^ 2, 0.2e1 * t18 * t59 + 0.2e1 * t22 * t41 (-t17 * t52 + t18 * t54) * t85, -0.2e1 * t17 * t59 - 0.2e1 * t22 * t42, t17 ^ 2 + t18 ^ 2 + t22 ^ 2, t24 ^ 2, t24 * t88, t24 * t84, t59 * t88, t59 ^ 2, 0.2e1 * t1 * t59 - 0.2e1 * t16 * t23, 0.2e1 * t16 * t24 - 0.2e1 * t2 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, -t78, t79, t69, -t26 * pkin(3) + t65, -t78, t69, -t79, t26 * t33 + t65, 0, 0, 0, 0, 0, -t26 * t30, -t26 * t31; 0, 0, 0, 0, 0, 0, t56, t59, 0, -t47, -t80, t38 + (-t81 - t82) * t56, pkin(8) * t41 + (-pkin(3) * t56 + t72) * t54, t63, -pkin(3) * t47 + t63 * qJ(4), -t22 * t54 + t33 * t41 + t38, t64, -t22 * t52 + (-t33 * t56 - t72) * t54, t64 * qJ(4) + t22 * t33, t24 * t31, t31 * t23 - t24 * t30, t31 * t59, -t30 * t59, 0, t14 * t59 + t16 * t30 - t28 * t23, -t15 * t59 + t16 * t31 + t28 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t54, -0.2e1 * t82, t32, pkin(3) ^ 2 + t73, t54 * t86, t32, t52 * t86, t33 ^ 2 + t73, t31 ^ 2, -0.2e1 * t31 * t30, 0, 0, 0, t30 * t87, t31 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, 0, t47, t41, 0, -t42, t22, 0, 0, 0, 0, 0, t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t52, 0, -pkin(3), -t54, 0, -t52, t33, 0, 0, 0, 0, 0, -t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t42, 0, t18, 0, 0, 0, 0, 0, t58 * t59, -t55 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t23, t59, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t7;
