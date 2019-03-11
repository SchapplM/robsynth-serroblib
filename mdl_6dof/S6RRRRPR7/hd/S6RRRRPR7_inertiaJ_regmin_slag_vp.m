% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t82 = sin(pkin(6));
t92 = cos(qJ(2));
t115 = t82 * t92;
t88 = sin(qJ(2));
t116 = t82 * t88;
t84 = cos(pkin(6));
t87 = sin(qJ(3));
t91 = cos(qJ(3));
t61 = t87 * t116 - t84 * t91;
t62 = t91 * t116 + t84 * t87;
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t40 = t90 * t61 + t86 * t62;
t41 = -t86 * t61 + t90 * t62;
t81 = sin(pkin(12));
t83 = cos(pkin(12));
t29 = -t81 * t40 + t83 * t41;
t85 = sin(qJ(6));
t89 = cos(qJ(6));
t19 = t89 * t115 + t85 * t29;
t128 = -0.2e1 * t19;
t76 = -t91 * pkin(3) - pkin(2);
t127 = 0.2e1 * t76;
t126 = 0.2e1 * t85;
t125 = -0.2e1 * t89;
t124 = 0.2e1 * t91;
t123 = pkin(9) + pkin(10);
t122 = pkin(1) * t88;
t121 = pkin(1) * t92;
t102 = pkin(3) * t115;
t101 = pkin(8) * t115;
t54 = t101 + (pkin(9) + t122) * t84;
t55 = (-pkin(2) * t92 - pkin(9) * t88 - pkin(1)) * t82;
t35 = -t87 * t54 + t91 * t55;
t33 = -t62 * pkin(10) - t102 + t35;
t36 = t91 * t54 + t87 * t55;
t34 = -t61 * pkin(10) + t36;
t16 = t90 * t33 - t86 * t34;
t11 = -pkin(4) * t115 - t41 * qJ(5) + t16;
t106 = t90 * t34;
t17 = t86 * t33 + t106;
t15 = -t40 * qJ(5) + t17;
t6 = t83 * t11 - t81 * t15;
t4 = pkin(5) * t115 - t6;
t120 = t4 * t89;
t119 = t86 * pkin(3);
t20 = -t85 * t115 + t89 * t29;
t18 = t20 * t85;
t97 = t123 * t91;
t98 = t123 * t87;
t50 = -t86 * t98 + t90 * t97;
t65 = t86 * t87 - t90 * t91;
t39 = -t65 * qJ(5) + t50;
t49 = -t86 * t97 - t90 * t98;
t66 = t86 * t91 + t90 * t87;
t94 = -t66 * qJ(5) + t49;
t25 = t81 * t39 - t83 * t94;
t118 = t25 * t89;
t78 = t82 ^ 2;
t117 = t78 * t92;
t114 = t84 * t88;
t28 = t83 * t40 + t81 * t41;
t22 = t85 * t28;
t47 = t83 * t65 + t81 * t66;
t44 = t85 * t47;
t48 = -t81 * t65 + t83 * t66;
t113 = t85 * t48;
t77 = t90 * pkin(3);
t75 = t77 + pkin(4);
t59 = t83 * t119 + t81 * t75;
t57 = pkin(11) + t59;
t112 = t85 * t57;
t73 = t81 * pkin(4) + pkin(11);
t111 = t85 * t73;
t110 = t85 * t89;
t109 = t89 * t48;
t108 = t89 * t57;
t107 = t89 * t73;
t7 = t81 * t11 + t83 * t15;
t58 = -t81 * t119 + t83 * t75;
t56 = -pkin(5) - t58;
t74 = -t83 * pkin(4) - pkin(5);
t105 = t56 + t74;
t104 = -0.2e1 * t115;
t103 = 0.2e1 * t115;
t100 = t87 * t115;
t99 = t91 * t115;
t96 = -t47 * t57 + t48 * t56;
t95 = -t47 * t73 + t48 * t74;
t52 = t65 * pkin(4) + t76;
t68 = pkin(8) * t116;
t53 = t68 + (-pkin(2) - t121) * t84;
t43 = t61 * pkin(3) + t53;
t30 = t40 * pkin(4) + t43;
t80 = t89 ^ 2;
t79 = t85 ^ 2;
t71 = t78 * t92 ^ 2;
t70 = 0.2e1 * t110;
t64 = pkin(1) * t114 + t101;
t63 = t84 * t121 - t68;
t46 = t48 ^ 2;
t45 = t89 * t47;
t42 = t85 * t109;
t32 = (-t79 + t80) * t48;
t27 = t83 * t39 + t81 * t94;
t24 = t47 * pkin(5) - t48 * pkin(11) + t52;
t23 = t89 * t28;
t21 = t25 * t85;
t14 = t85 * t24 + t89 * t27;
t13 = t89 * t24 - t85 * t27;
t9 = -t85 * t19 + t20 * t89;
t8 = t28 * pkin(5) - t29 * pkin(11) + t30;
t5 = -pkin(11) * t115 + t7;
t3 = t4 * t85;
t2 = t89 * t5 + t85 * t8;
t1 = -t85 * t5 + t89 * t8;
t10 = [1, 0, 0, t78 * t88 ^ 2, 0.2e1 * t88 * t117, 0.2e1 * t82 * t114, t84 * t103, t84 ^ 2, 0.2e1 * pkin(1) * t117 + 0.2e1 * t63 * t84, -0.2e1 * t78 * t122 - 0.2e1 * t64 * t84, t62 ^ 2, -0.2e1 * t62 * t61, t62 * t104, t61 * t103, t71, -0.2e1 * t35 * t115 + 0.2e1 * t53 * t61, 0.2e1 * t36 * t115 + 0.2e1 * t53 * t62, t41 ^ 2, -0.2e1 * t41 * t40, t41 * t104, t40 * t103, t71, -0.2e1 * t115 * t16 + 0.2e1 * t43 * t40, 0.2e1 * t115 * t17 + 0.2e1 * t43 * t41, -0.2e1 * t7 * t28 - 0.2e1 * t6 * t29, t30 ^ 2 + t6 ^ 2 + t7 ^ 2, t20 ^ 2, t20 * t128, 0.2e1 * t20 * t28, t28 * t128, t28 ^ 2, 0.2e1 * t1 * t28 + 0.2e1 * t4 * t19, -0.2e1 * t2 * t28 + 0.2e1 * t4 * t20; 0, 0, 0, 0, 0, t116, t115, t84, t63, -t64, t62 * t87, -t87 * t61 + t62 * t91, -t100, -t99, 0, -pkin(2) * t61 + pkin(9) * t100 - t53 * t91, -pkin(2) * t62 + pkin(9) * t99 + t53 * t87, t41 * t66, -t66 * t40 - t41 * t65, -t66 * t115, t65 * t115, 0, -t115 * t49 + t76 * t40 + t43 * t65, t115 * t50 + t76 * t41 + t43 * t66, t25 * t29 - t27 * t28 - t7 * t47 - t6 * t48, -t6 * t25 + t7 * t27 + t30 * t52, t20 * t109 (-t19 * t89 - t18) * t48, t109 * t28 + t20 * t47, -t113 * t28 - t19 * t47, t28 * t47, t1 * t47 + t113 * t4 + t13 * t28 + t25 * t19, t109 * t4 - t14 * t28 - t2 * t47 + t25 * t20; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t87 ^ 2, t87 * t124, 0, 0, 0, pkin(2) * t124, -0.2e1 * pkin(2) * t87, t66 ^ 2, -0.2e1 * t66 * t65, 0, 0, 0, t65 * t127, t66 * t127, 0.2e1 * t25 * t48 - 0.2e1 * t27 * t47, t25 ^ 2 + t27 ^ 2 + t52 ^ 2, t80 * t46, -0.2e1 * t46 * t110, 0.2e1 * t47 * t109, -0.2e1 * t47 * t113, t47 ^ 2, 0.2e1 * t113 * t25 + 0.2e1 * t13 * t47, 0.2e1 * t109 * t25 - 0.2e1 * t14 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t61, -t115, t35, -t36, 0, 0, t41, -t40, -t115, -t102 * t90 + t16, -t106 + (-t33 + t102) * t86, -t59 * t28 - t58 * t29, t6 * t58 + t7 * t59, t18, t9, t22, t23, 0, -t112 * t28 + t56 * t19 - t120, -t108 * t28 + t56 * t20 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t91, 0, -t87 * pkin(9), -t91 * pkin(9), 0, 0, t66, -t65, 0, t49, -t50, -t59 * t47 - t58 * t48, -t25 * t58 + t27 * t59, t42, t32, t44, t45, 0, t85 * t96 - t118, t89 * t96 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t119, 0, t58 ^ 2 + t59 ^ 2, t79, t70, 0, 0, 0, t56 * t125, t56 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, -t115, t16, -t17 (-t28 * t81 - t29 * t83) * pkin(4) (t6 * t83 + t7 * t81) * pkin(4), t18, t9, t22, t23, 0, -t111 * t28 + t74 * t19 - t120, -t107 * t28 + t74 * t20 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t65, 0, t49, -t50 (-t47 * t81 - t48 * t83) * pkin(4) (-t25 * t83 + t27 * t81) * pkin(4), t42, t32, t44, t45, 0, t85 * t95 - t118, t89 * t95 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t77, -t119, 0 (t58 * t83 + t59 * t81) * pkin(4), t79, t70, 0, 0, 0, -t105 * t89, t105 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t81 ^ 2 + t83 ^ 2) * pkin(4) ^ 2, t79, t70, 0, 0, 0, t74 * t125, t74 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, t23, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, t45, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t28, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t113, t47, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t89, 0, -t112, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t89, 0, -t111, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
