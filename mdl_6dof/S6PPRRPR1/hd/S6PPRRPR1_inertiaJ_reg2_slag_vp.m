% Calculate inertial parameters regressor of joint inertia matrix for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:05:13
% EndTime: 2019-05-04 20:05:17
% DurationCPUTime: 1.50s
% Computational Cost: add. (1054->160), mult. (2791->300), div. (0->0), fcn. (3499->14), ass. (0->85)
t59 = sin(pkin(13));
t63 = cos(pkin(13));
t67 = sin(qJ(6));
t70 = cos(qJ(6));
t102 = -t67 * t59 + t70 * t63;
t60 = sin(pkin(12));
t62 = sin(pkin(6));
t66 = cos(pkin(6));
t69 = sin(qJ(3));
t72 = cos(qJ(3));
t64 = cos(pkin(12));
t65 = cos(pkin(7));
t88 = t64 * t65;
t61 = sin(pkin(7));
t92 = t61 * t69;
t17 = t66 * t92 + (t60 * t72 + t69 * t88) * t62;
t29 = -t62 * t64 * t61 + t66 * t65;
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t11 = t17 * t68 - t29 * t71;
t10 = t11 ^ 2;
t91 = t61 * t72;
t15 = -t66 * t91 + (t60 * t69 - t72 * t88) * t62;
t101 = t15 ^ 2;
t31 = -t71 * t65 + t68 * t92;
t30 = t31 ^ 2;
t50 = -t63 * pkin(5) - pkin(4);
t100 = 0.2e1 * t50;
t99 = -0.2e1 * t71;
t57 = t68 ^ 2;
t98 = t57 * pkin(9);
t52 = t68 * pkin(9);
t5 = t11 * t31;
t97 = t11 * t68;
t96 = t31 * t68;
t95 = t59 * t63;
t94 = t59 * t68;
t93 = t59 * t71;
t90 = t63 * t68;
t89 = t63 * t71;
t85 = t71 * t68;
t84 = pkin(10) + qJ(5);
t41 = -t71 * pkin(4) - t68 * qJ(5) - pkin(3);
t24 = pkin(9) * t89 + t59 * t41;
t53 = t59 ^ 2;
t56 = t63 ^ 2;
t83 = t53 + t56;
t82 = 0.2e1 * t85;
t81 = t59 * t90;
t13 = t17 * t71 + t29 * t68;
t3 = -t13 * t59 + t15 * t63;
t4 = t13 * t63 + t15 * t59;
t80 = -t3 * t59 + t4 * t63;
t79 = -pkin(4) * t68 + qJ(5) * t71;
t78 = t13 * t71 + t97;
t33 = t68 * t65 + t71 * t92;
t18 = -t59 * t33 - t63 * t91;
t19 = t63 * t33 - t59 * t91;
t77 = -t18 * t59 + t19 * t63;
t35 = t63 * t41;
t23 = -pkin(9) * t93 + t35;
t76 = -t23 * t59 + t24 * t63;
t75 = t33 * t71 + t96;
t38 = t70 * t59 + t67 * t63;
t74 = pkin(9) ^ 2;
t58 = t71 ^ 2;
t54 = t61 ^ 2;
t51 = t57 * t74;
t47 = t54 * t72 ^ 2;
t43 = t84 * t63;
t42 = t84 * t59;
t40 = pkin(5) * t94 + t52;
t28 = t102 * t68;
t26 = t38 * t68;
t22 = -t67 * t42 + t70 * t43;
t21 = -t70 * t42 - t67 * t43;
t20 = -pkin(10) * t94 + t24;
t14 = -pkin(10) * t90 + t35 + (-pkin(9) * t59 - pkin(5)) * t71;
t9 = t67 * t18 + t70 * t19;
t8 = t70 * t18 - t67 * t19;
t7 = t67 * t14 + t70 * t20;
t6 = t70 * t14 - t67 * t20;
t2 = t67 * t3 + t70 * t4;
t1 = t70 * t3 - t67 * t4;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 ^ 2 + (t60 ^ 2 + t64 ^ 2) * t62 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 + t29 ^ 2 + t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t10 + t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t65 + (-t15 * t72 + t17 * t69) * t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t33 - t15 * t91 + t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t18 + t4 * t19 + t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t8 + t2 * t9 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t69 ^ 2 + t65 ^ 2 + t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 ^ 2 + t30 + t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 + t19 ^ 2 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 ^ 2 + t9 ^ 2 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t71, t15 * t68, t78, -t15 * pkin(3) + t78 * pkin(9), 0, 0, 0, 0, 0, 0, t11 * t94 - t3 * t71, t11 * t90 + t4 * t71 (-t3 * t63 - t4 * t59) * t68, pkin(9) * t97 + t3 * t23 + t4 * t24, 0, 0, 0, 0, 0, 0, -t1 * t71 + t11 * t26, t11 * t28 + t2 * t71, -t1 * t28 - t2 * t26, t1 * t6 + t11 * t40 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t92, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t91, -t68 * t91, t75, pkin(3) * t91 + t75 * pkin(9), 0, 0, 0, 0, 0, 0, -t18 * t71 + t31 * t94, t19 * t71 + t31 * t90 (-t18 * t63 - t19 * t59) * t68, pkin(9) * t96 + t18 * t23 + t19 * t24, 0, 0, 0, 0, 0, 0, t31 * t26 - t8 * t71, t31 * t28 + t9 * t71, -t9 * t26 - t8 * t28, t31 * t40 + t8 * t6 + t9 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t57, t82, 0, t58, 0, 0, 0.2e1 * pkin(3) * t71, -0.2e1 * pkin(3) * t68, 0.2e1 * (t57 + t58) * pkin(9), pkin(3) ^ 2 + t58 * t74 + t51, t56 * t57, -0.2e1 * t57 * t95, -0.2e1 * t63 * t85, t53 * t57, t59 * t82, t58, -0.2e1 * t23 * t71 + 0.2e1 * t59 * t98, 0.2e1 * t24 * t71 + 0.2e1 * t63 * t98, 0.2e1 * (-t23 * t63 - t24 * t59) * t68, t23 ^ 2 + t24 ^ 2 + t51, t28 ^ 2, -0.2e1 * t28 * t26, t28 * t99, t26 ^ 2, -t26 * t99, t58, 0.2e1 * t40 * t26 - 0.2e1 * t6 * t71, 0.2e1 * t40 * t28 + 0.2e1 * t7 * t71, -0.2e1 * t7 * t26 - 0.2e1 * t6 * t28, t40 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t13, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t63, t11 * t59, t80, -t11 * pkin(4) + t80 * qJ(5), 0, 0, 0, 0, 0, 0, -t11 * t102, t11 * t38, -t1 * t38 + t102 * t2, t1 * t21 + t11 * t50 + t2 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t63, t31 * t59, t77, -t31 * pkin(4) + t77 * qJ(5), 0, 0, 0, 0, 0, 0, -t31 * t102, t31 * t38, t102 * t9 - t8 * t38, t8 * t21 + t9 * t22 + t31 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, t71, 0, -t52, -t71 * pkin(9), 0, 0, t81 (-t53 + t56) * t68, -t93, -t81, -t89, 0, -pkin(9) * t90 + t79 * t59, pkin(9) * t94 + t79 * t63, t76, -pkin(4) * t52 + t76 * qJ(5), t28 * t38, t102 * t28 - t38 * t26, -t38 * t71, -t26 * t102, -t102 * t71, 0, -t102 * t40 - t21 * t71 + t50 * t26, t22 * t71 + t50 * t28 + t40 * t38, t102 * t7 - t21 * t28 - t22 * t26 - t6 * t38, t6 * t21 + t7 * t22 + t40 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t53, 0.2e1 * t95, 0, t56, 0, 0, 0.2e1 * pkin(4) * t63, -0.2e1 * pkin(4) * t59, 0.2e1 * t83 * qJ(5), t83 * qJ(5) ^ 2 + pkin(4) ^ 2, t38 ^ 2, 0.2e1 * t38 * t102, 0, t102 ^ 2, 0, 0, -t102 * t100, t38 * t100, 0.2e1 * t102 * t22 - 0.2e1 * t21 * t38, t21 ^ 2 + t22 ^ 2 + t50 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t90, 0, t52, 0, 0, 0, 0, 0, 0, t26, t28, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t59, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t102, t38, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, -t71, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t102, 0, t21, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t12;
