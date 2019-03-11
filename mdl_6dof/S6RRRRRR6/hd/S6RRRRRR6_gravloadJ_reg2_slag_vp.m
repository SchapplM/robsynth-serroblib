% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t107 = cos(pkin(6));
t71 = sin(pkin(6));
t74 = sin(qJ(2));
t120 = t71 * t74;
t73 = sin(qJ(3));
t77 = cos(qJ(3));
t143 = t107 * t77 - t73 * t120;
t118 = t71 * t77;
t134 = cos(qJ(1));
t78 = cos(qJ(2));
t75 = sin(qJ(1));
t96 = t75 * t107;
t48 = t134 * t78 - t74 * t96;
t30 = t75 * t118 - t48 * t73;
t103 = t71 * t134;
t90 = t107 * t134;
t46 = t74 * t90 + t75 * t78;
t70 = qJ(3) + qJ(4);
t65 = sin(t70);
t67 = cos(t70);
t25 = -t65 * t103 + t46 * t67;
t45 = t75 * t74 - t78 * t90;
t72 = sin(qJ(5));
t76 = cos(qJ(5));
t142 = t25 * t72 - t45 * t76;
t131 = t45 * t72;
t141 = t25 * t76 + t131;
t69 = qJ(5) + qJ(6);
t64 = sin(t69);
t66 = cos(t69);
t140 = t25 * t64 - t45 * t66;
t139 = t25 * t66 + t45 * t64;
t114 = t76 * t78;
t119 = t71 * t75;
t29 = t119 * t65 + t48 * t67;
t47 = t134 * t74 + t78 * t96;
t14 = -t29 * t72 + t47 * t76;
t40 = t107 * t65 + t120 * t67;
t138 = g(2) * t142 - g(3) * (-t114 * t71 - t40 * t72) - g(1) * t14;
t117 = t71 * t78;
t63 = t77 * pkin(3) + pkin(2);
t49 = t63 * t117;
t136 = g(3) * t49;
t135 = g(3) * t71;
t129 = t46 * t72;
t128 = t47 * t72;
t127 = t48 * t72;
t125 = t64 * t67;
t124 = t66 * t67;
t123 = t67 * t72;
t122 = t67 * t76;
t121 = t67 * t78;
t116 = t72 * t78;
t80 = -pkin(10) - pkin(9);
t115 = t74 * t80;
t62 = t76 * pkin(5) + pkin(4);
t79 = -pkin(12) - pkin(11);
t98 = -t67 * t103 - t46 * t65;
t113 = -t25 * t79 + t62 * t98;
t28 = -t119 * t67 + t48 * t65;
t112 = -t28 * t62 - t29 * t79;
t39 = t107 * t67 - t120 * t65;
t111 = t39 * t62 - t40 * t79;
t110 = -t45 * t63 - t46 * t80;
t109 = -t47 * t63 - t48 * t80;
t108 = t134 * pkin(1) + pkin(8) * t119;
t105 = t73 * t119;
t102 = -t75 * pkin(1) + pkin(8) * t103;
t101 = pkin(4) * t98 + t25 * pkin(11);
t100 = -t28 * pkin(4) + t29 * pkin(11);
t99 = t39 * pkin(4) + t40 * pkin(11);
t57 = t73 * t103;
t97 = t46 * t77 - t57;
t94 = t30 * pkin(3);
t93 = pkin(3) * t105 - t47 * t80 + t48 * t63 + t108;
t92 = pkin(4) * t67 + pkin(11) * t65;
t91 = g(1) * t98 + g(2) * t28;
t21 = g(1) * t45 - g(2) * t47;
t89 = t62 * t67 - t65 * t79;
t88 = t143 * pkin(3);
t87 = g(1) * t134 + g(2) * t75;
t86 = pkin(3) * t57 + t45 * t80 - t46 * t63 + t102;
t7 = g(1) * t28 - g(2) * t98 - g(3) * t39;
t9 = g(1) * t29 + g(2) * t25 + g(3) * t40;
t85 = t103 * t77 + t46 * t73;
t84 = -g(1) * t47 - g(2) * t45 + g(3) * t117;
t83 = g(1) * t48 + g(2) * t46 + g(3) * t120;
t82 = t85 * pkin(3);
t31 = t48 * t77 + t105;
t15 = t29 * t76 + t128;
t13 = t84 * t65;
t12 = t29 * t66 + t47 * t64;
t11 = -t29 * t64 + t47 * t66;
t6 = t7 * t76;
t5 = t7 * t72;
t4 = t7 * t66;
t3 = t7 * t64;
t2 = g(1) * t12 + g(2) * t139 - g(3) * (t117 * t64 - t40 * t66);
t1 = -g(1) * t11 + g(2) * t140 - g(3) * (-t117 * t66 - t40 * t64);
t8 = [0, 0, 0, 0, 0, 0, g(1) * t75 - g(2) * t134, t87, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t48, -t21, -t87 * t71, -g(1) * t102 - g(2) * t108, 0, 0, 0, 0, 0, 0, g(1) * t97 - g(2) * t31, -g(1) * t85 - g(2) * t30, t21, -g(1) * (-t46 * pkin(2) - t45 * pkin(9) + t102) - g(2) * (t48 * pkin(2) + t47 * pkin(9) + t108) 0, 0, 0, 0, 0, 0, g(1) * t25 - g(2) * t29, t91, t21, -g(1) * t86 - g(2) * t93, 0, 0, 0, 0, 0, 0, g(1) * t141 - g(2) * t15, -g(1) * t142 - g(2) * t14, -t91, -g(1) * (-pkin(4) * t25 + pkin(11) * t98 + t86) - g(2) * (t29 * pkin(4) + t28 * pkin(11) + t93) 0, 0, 0, 0, 0, 0, g(1) * t139 - g(2) * t12, -g(1) * t140 - g(2) * t11, -t91, -g(1) * (-pkin(5) * t131 - t25 * t62 - t79 * t98 + t86) - g(2) * (pkin(5) * t128 - t28 * t79 + t29 * t62 + t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t83, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t77, t84 * t73, -t83, -g(1) * (-t47 * pkin(2) + t48 * pkin(9)) - g(2) * (-t45 * pkin(2) + t46 * pkin(9)) - (pkin(2) * t78 + pkin(9) * t74) * t135, 0, 0, 0, 0, 0, 0, -t84 * t67, t13, -t83, -g(1) * t109 - g(2) * t110 - g(3) * (-t115 * t71 + t49) 0, 0, 0, 0, 0, 0, -g(1) * (-t122 * t47 + t127) - g(2) * (-t122 * t45 + t129) - (t114 * t67 + t72 * t74) * t135, -g(1) * (t123 * t47 + t48 * t76) - g(2) * (t123 * t45 + t46 * t76) - (-t116 * t67 + t74 * t76) * t135, -t13, -g(1) * (-t47 * t92 + t109) - g(2) * (-t45 * t92 + t110) - t136 - (t78 * t92 - t115) * t135, 0, 0, 0, 0, 0, 0, -g(1) * (-t124 * t47 + t48 * t64) - g(2) * (-t124 * t45 + t46 * t64) - (t121 * t66 + t64 * t74) * t135, -g(1) * (t125 * t47 + t48 * t66) - g(2) * (t125 * t45 + t46 * t66) - (-t121 * t64 + t66 * t74) * t135, -t13, -g(1) * (pkin(5) * t127 - t47 * t89 + t109) - g(2) * (pkin(5) * t129 - t45 * t89 + t110) - t136 - (t89 * t78 + (pkin(5) * t72 - t80) * t74) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t30 + g(2) * t85 - g(3) * t143, g(1) * t31 + g(2) * t97 - g(3) * (-t107 * t73 - t118 * t74) 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, -g(1) * t94 + g(2) * t82 - g(3) * t88, 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * (t100 + t94) - g(2) * (t101 - t82) - g(3) * (t88 + t99) 0, 0, 0, 0, 0, 0, t4, -t3, -t9, -g(1) * (t94 + t112) - g(2) * (-t82 + t113) - g(3) * (t88 + t111); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * t100 - g(2) * t101 - g(3) * t99, 0, 0, 0, 0, 0, 0, t4, -t3, -t9, -g(1) * t112 - g(2) * t113 - g(3) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, g(1) * t15 + g(2) * t141 - g(3) * (t116 * t71 - t40 * t76) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t138 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t8;
