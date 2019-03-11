% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t69 = sin(pkin(13));
t73 = cos(pkin(13));
t79 = sin(qJ(3));
t83 = cos(qJ(3));
t104 = t83 * t69 + t79 * t73;
t72 = sin(pkin(6));
t84 = cos(qJ(1));
t128 = t72 * t84;
t57 = t79 * t69 - t83 * t73;
t71 = sin(pkin(7));
t46 = t57 * t71;
t75 = cos(pkin(7));
t48 = t57 * t75;
t74 = cos(pkin(12));
t120 = t84 * t74;
t70 = sin(pkin(12));
t80 = sin(qJ(1));
t124 = t80 * t70;
t76 = cos(pkin(6));
t53 = -t120 * t76 + t124;
t121 = t84 * t70;
t123 = t80 * t74;
t54 = t121 * t76 + t123;
t19 = t104 * t54 - t128 * t46 - t53 * t48;
t36 = t128 * t75 - t53 * t71;
t78 = sin(qJ(5));
t82 = cos(qJ(5));
t47 = t104 * t71;
t49 = t104 * t75;
t90 = t128 * t47 + t53 * t49 + t54 * t57;
t7 = t36 * t78 + t82 * t90;
t77 = sin(qJ(6));
t81 = cos(qJ(6));
t143 = t19 * t81 + t7 * t77;
t142 = -t19 * t77 + t7 * t81;
t130 = t72 * t74;
t139 = t75 * t130 + t71 * t76;
t138 = -t36 * t82 + t78 * t90;
t129 = t72 * t80;
t55 = -t123 * t76 - t121;
t101 = t129 * t75 - t55 * t71;
t29 = -g(1) * t36 - g(2) * t101;
t28 = t76 * t47 + (t49 * t74 - t57 * t70) * t72;
t136 = t83 * pkin(3);
t135 = t54 * t79;
t56 = -t124 * t76 + t120;
t134 = t56 * t79;
t133 = t70 * t72;
t131 = t71 * t79;
t127 = t75 * t79;
t125 = t77 * t82;
t122 = t81 * t82;
t119 = pkin(9) + qJ(4);
t117 = qJ(2) * t72;
t118 = t84 * pkin(1) + t80 * t117;
t116 = t75 * t136;
t115 = t71 * t129;
t114 = t71 * t128;
t112 = -t80 * pkin(1) + t84 * t117;
t23 = t129 * t47 + t55 * t49 - t56 * t57;
t8 = -t101 * t82 + t23 * t78;
t111 = g(1) * t138 + g(2) * t8;
t50 = pkin(3) * t131 + t119 * t75;
t51 = pkin(3) * t127 - t119 * t71;
t67 = pkin(2) + t136;
t110 = t50 * t129 + t55 * t51 + t56 * t67 + t118;
t109 = pkin(5) * t82 + pkin(11) * t78;
t22 = -t104 * t56 - t129 * t46 - t55 * t48;
t108 = g(1) * t19 + g(2) * t22;
t107 = g(1) * t84 + g(2) * t80;
t106 = g(1) * t80 - g(2) * t84;
t103 = -pkin(3) * t134 + t115 * t136 + t55 * t116;
t102 = t53 * t75 + t114;
t100 = t55 * t75 + t115;
t99 = -pkin(3) * t133 * t79 + t139 * t136;
t98 = t50 * t128 + t53 * t51 - t54 * t67 + t112;
t52 = -t130 * t71 + t76 * t75;
t14 = -t28 * t78 + t52 * t82;
t97 = g(1) * t8 - g(2) * t138 - g(3) * t14;
t15 = t28 * t82 + t52 * t78;
t9 = t101 * t78 + t23 * t82;
t96 = g(1) * t9 - g(2) * t7 + g(3) * t15;
t95 = t79 * t114 + t127 * t53 - t54 * t83;
t94 = -g(1) * t23 + g(2) * t90 - g(3) * t28;
t27 = -t76 * t46 + (-t104 * t70 - t48 * t74) * t72;
t93 = g(1) * t22 - g(2) * t19 + g(3) * t27;
t91 = t23 * pkin(4) - t22 * pkin(10) + t110;
t89 = t22 * pkin(4) + pkin(10) * t23 + t103;
t88 = t27 * pkin(4) + pkin(10) * t28 + t99;
t87 = pkin(4) * t90 - t19 * pkin(10) + t98;
t86 = -t53 * t116 + (-t114 * t83 - t135) * pkin(3);
t85 = -pkin(4) * t19 - pkin(10) * t90 + t86;
t45 = -g(3) * t76 - t106 * t72;
t31 = t100 * t79 + t56 * t83;
t30 = t100 * t83 - t134;
t24 = -g(1) * t101 + g(2) * t36 - g(3) * t52;
t3 = -t22 * t77 + t9 * t81;
t2 = -t22 * t81 - t9 * t77;
t1 = t93 * t78;
t4 = [0, 0, 0, 0, 0, 0, t106, t107, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t56, -g(1) * t53 - g(2) * t55, -t107 * t72, -g(1) * t112 - g(2) * t118, 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t31, -g(1) * (t102 * t83 + t135) - g(2) * t30, t29, -g(1) * (-t54 * pkin(2) + t112) - g(2) * (t56 * pkin(2) + t118) + t29 * pkin(9), 0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t23, -t108, t29, -g(1) * t98 - g(2) * t110, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t111, t108, -g(1) * t87 - g(2) * t91, 0, 0, 0, 0, 0, 0, -g(1) * t142 - g(2) * t3, g(1) * t143 - g(2) * t2, -t111, -g(1) * (t7 * pkin(5) + pkin(11) * t138 + t87) - g(2) * (t9 * pkin(5) + t8 * pkin(11) + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t30 + (g(2) * t54 + g(3) * t133) * t79 + (g(2) * t102 - g(3) * t139) * t83, g(1) * t31 - g(2) * t95 - g(3) * (-t76 * t131 + (-t127 * t74 - t70 * t83) * t72) 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t94, 0, -g(1) * t103 - g(2) * t86 - g(3) * t99, 0, 0, 0, 0, 0, 0, -t93 * t82, t1, t94, -g(1) * t89 - g(2) * t85 - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (t122 * t22 + t23 * t77) - g(2) * (-t122 * t19 - t77 * t90) - g(3) * (t122 * t27 + t28 * t77) -g(1) * (-t125 * t22 + t23 * t81) - g(2) * (t125 * t19 - t81 * t90) - g(3) * (-t125 * t27 + t28 * t81) -t1, -g(1) * (t109 * t22 + t89) - g(2) * (-t109 * t19 + t85) - g(3) * (t109 * t27 + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t96, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t81, -t97 * t77, -t96, -g(1) * (-t8 * pkin(5) + t9 * pkin(11)) - g(2) * (pkin(5) * t138 - pkin(11) * t7) - g(3) * (t14 * pkin(5) + t15 * pkin(11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t143 - g(3) * (-t15 * t77 - t27 * t81) g(1) * t3 - g(2) * t142 - g(3) * (-t15 * t81 + t27 * t77) 0, 0;];
taug_reg  = t4;
