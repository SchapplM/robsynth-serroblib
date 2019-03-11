% Calculate inertial parameters regressor of gravitation load for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t120 = cos(pkin(8));
t52 = sin(pkin(8));
t58 = cos(qJ(3));
t116 = sin(pkin(7));
t117 = sin(pkin(6));
t100 = t117 * t116;
t119 = cos(pkin(13));
t121 = cos(pkin(7));
t118 = cos(pkin(14));
t103 = t119 * t118;
t122 = cos(pkin(6));
t114 = sin(pkin(14));
t115 = sin(pkin(13));
t96 = t115 * t114;
t78 = t122 * t103 - t96;
t67 = -t119 * t100 + t78 * t121;
t127 = sin(qJ(3));
t101 = t119 * t114;
t97 = t115 * t118;
t79 = t122 * t101 + t97;
t73 = t79 * t127;
t61 = -t67 * t58 + t73;
t102 = t119 * t117;
t66 = -t121 * t102 - t78 * t116;
t132 = t61 * t120 - t66 * t52;
t80 = -t122 * t97 - t101;
t99 = t117 * t115;
t69 = t116 * t99 + t80 * t121;
t81 = -t122 * t96 + t103;
t74 = t81 * t127;
t62 = -t69 * t58 + t74;
t68 = -t80 * t116 + t121 * t99;
t131 = t62 * t120 - t68 * t52;
t77 = t118 * t121 * t117 + t122 * t116;
t98 = t117 * t114;
t85 = t127 * t98;
t70 = -t77 * t58 + t85;
t76 = -t118 * t100 + t122 * t121;
t130 = t70 * t120 - t76 * t52;
t129 = pkin(10) * t52;
t128 = cos(qJ(4));
t54 = sin(qJ(5));
t126 = t52 * t54;
t57 = cos(qJ(5));
t125 = t52 * t57;
t53 = sin(qJ(6));
t124 = t53 * t57;
t56 = cos(qJ(6));
t123 = t56 * t57;
t40 = t67 * t127 + t79 * t58;
t55 = sin(qJ(4));
t12 = t132 * t128 + t40 * t55;
t13 = t40 * t128 - t132 * t55;
t113 = -t12 * pkin(4) + pkin(11) * t13;
t41 = t69 * t127 + t81 * t58;
t14 = t131 * t128 + t41 * t55;
t15 = t41 * t128 - t131 * t55;
t112 = -t14 * pkin(4) + pkin(11) * t15;
t48 = t77 * t127 + t58 * t98;
t27 = t130 * t128 + t48 * t55;
t28 = t48 * t128 - t130 * t55;
t111 = -t27 * pkin(4) + pkin(11) * t28;
t110 = t55 * t120;
t109 = -t61 * pkin(3) + t40 * t129;
t108 = -t62 * pkin(3) + t41 * t129;
t107 = -t70 * pkin(3) + t48 * t129;
t106 = -pkin(5) * t57 - pkin(12) * t54;
t105 = t120 * t128;
t42 = t76 * t120 + t70 * t52;
t16 = -t28 * t54 + t42 * t57;
t29 = t66 * t120 + t61 * t52;
t2 = -t13 * t54 + t29 * t57;
t30 = t68 * t120 + t62 * t52;
t4 = -t15 * t54 + t30 * t57;
t95 = g(1) * t4 + g(2) * t2 + g(3) * t16;
t17 = t28 * t57 + t42 * t54;
t3 = t13 * t57 + t29 * t54;
t5 = t15 * t57 + t30 * t54;
t94 = g(1) * t5 + g(2) * t3 + g(3) * t17;
t33 = -t48 * t110 - t70 * t128;
t24 = -t48 * t125 + t33 * t54;
t21 = -t40 * t110 - t61 * t128;
t6 = -t40 * t125 + t21 * t54;
t23 = -t41 * t110 - t62 * t128;
t8 = -t41 * t125 + t23 * t54;
t93 = g(1) * t8 + g(2) * t6 + g(3) * t24;
t92 = g(1) * t14 + g(2) * t12 + g(3) * t27;
t91 = g(1) * t15 + g(2) * t13 + g(3) * t28;
t20 = t40 * t105 - t61 * t55;
t22 = t41 * t105 - t62 * t55;
t32 = t48 * t105 - t70 * t55;
t90 = g(1) * t22 + g(2) * t20 + g(3) * t32;
t89 = g(1) * t41 + g(2) * t40 + g(3) * t48;
t88 = t21 * pkin(4) + pkin(11) * t20 + t109;
t87 = t23 * pkin(4) + pkin(11) * t22 + t108;
t86 = t33 * pkin(4) + pkin(11) * t32 + t107;
t50 = -g(1) * t99 + g(2) * t102 - g(3) * t122;
t25 = t48 * t126 + t33 * t57;
t9 = t41 * t126 + t23 * t57;
t7 = t40 * t126 + t21 * t57;
t1 = t92 * t54;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t74 + g(2) * t73 + g(3) * t85 + (-g(1) * t69 - g(2) * t67 - g(3) * t77) * t58, t89, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t21 - g(3) * t33, t90, -t89 * t52, -g(1) * t108 - g(2) * t109 - g(3) * t107, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t25, t93, -t90, -g(1) * t87 - g(2) * t88 - g(3) * t86, 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t53 + t56 * t9) - g(2) * (t20 * t53 + t56 * t7) - g(3) * (t25 * t56 + t32 * t53) -g(1) * (t22 * t56 - t53 * t9) - g(2) * (t20 * t56 - t53 * t7) - g(3) * (-t25 * t53 + t32 * t56) -t93, -g(1) * (t9 * pkin(5) + t8 * pkin(12) + t87) - g(2) * (t7 * pkin(5) + t6 * pkin(12) + t88) - g(3) * (pkin(5) * t25 + pkin(12) * t24 + t86); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t91, 0, 0, 0, 0, 0, 0, 0, 0, t92 * t57, -t1, -t91, -g(1) * t112 - g(2) * t113 - g(3) * t111, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t123 + t15 * t53) - g(2) * (-t12 * t123 + t13 * t53) - g(3) * (-t27 * t123 + t28 * t53) -g(1) * (t14 * t124 + t15 * t56) - g(2) * (t12 * t124 + t13 * t56) - g(3) * (t27 * t124 + t28 * t56) t1, -g(1) * (t106 * t14 + t112) - g(2) * (t106 * t12 + t113) - g(3) * (t106 * t27 + t111); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t94, 0, 0, 0, 0, 0, 0, 0, 0, -t95 * t56, t95 * t53, -t94, -g(1) * (pkin(5) * t4 + pkin(12) * t5) - g(2) * (pkin(5) * t2 + pkin(12) * t3) - g(3) * (pkin(5) * t16 + pkin(12) * t17); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t56 - t5 * t53) - g(2) * (t12 * t56 - t3 * t53) - g(3) * (-t17 * t53 + t27 * t56) -g(1) * (-t14 * t53 - t5 * t56) - g(2) * (-t12 * t53 - t3 * t56) - g(3) * (-t17 * t56 - t27 * t53) 0, 0;];
taug_reg  = t10;
