% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t103 = cos(pkin(6));
t68 = sin(pkin(6));
t72 = sin(qJ(2));
t116 = t68 * t72;
t71 = sin(qJ(3));
t75 = cos(qJ(3));
t131 = t103 * t75 - t71 * t116;
t114 = t68 * t75;
t125 = cos(qJ(1));
t76 = cos(qJ(2));
t73 = sin(qJ(1));
t93 = t73 * t103;
t48 = t125 * t76 - t72 * t93;
t30 = t73 * t114 - t48 * t71;
t87 = t103 * t125;
t46 = t72 * t87 + t73 * t76;
t67 = qJ(3) + qJ(4);
t64 = sin(t67);
t65 = cos(t67);
t99 = t68 * t125;
t25 = t46 * t65 - t64 * t99;
t45 = t73 * t72 - t76 * t87;
t70 = sin(qJ(5));
t74 = cos(qJ(5));
t130 = t25 * t70 - t45 * t74;
t124 = t45 * t70;
t129 = t25 * t74 + t124;
t110 = t74 * t76;
t115 = t68 * t73;
t29 = t64 * t115 + t48 * t65;
t47 = t125 * t72 + t76 * t93;
t14 = -t29 * t70 + t47 * t74;
t40 = t103 * t64 + t65 * t116;
t1 = g(2) * t130 - g(3) * (-t68 * t110 - t40 * t70) - g(1) * t14;
t113 = t68 * t76;
t63 = t75 * pkin(3) + pkin(2);
t49 = t63 * t113;
t127 = g(3) * t49;
t126 = g(3) * t68;
t122 = t46 * t70;
t121 = t47 * t70;
t120 = t48 * t70;
t118 = t65 * t70;
t117 = t65 * t74;
t112 = t70 * t76;
t77 = -pkin(10) - pkin(9);
t111 = t72 * t77;
t24 = t46 * t64 + t65 * t99;
t62 = t74 * pkin(5) + pkin(4);
t69 = -qJ(6) - pkin(11);
t109 = -t24 * t62 - t25 * t69;
t28 = -t65 * t115 + t48 * t64;
t108 = -t28 * t62 - t29 * t69;
t39 = -t103 * t65 + t64 * t116;
t107 = -t39 * t62 - t40 * t69;
t106 = -t45 * t63 - t46 * t77;
t105 = -t47 * t63 - t48 * t77;
t104 = t125 * pkin(1) + pkin(8) * t115;
t101 = t71 * t115;
t98 = -t73 * pkin(1) + pkin(8) * t99;
t97 = -t24 * pkin(4) + t25 * pkin(11);
t96 = -t28 * pkin(4) + t29 * pkin(11);
t95 = -t39 * pkin(4) + t40 * pkin(11);
t56 = t71 * t99;
t94 = t46 * t75 - t56;
t91 = t30 * pkin(3);
t90 = pkin(3) * t101 - t47 * t77 + t48 * t63 + t104;
t89 = pkin(4) * t65 + pkin(11) * t64;
t88 = -g(1) * t24 + g(2) * t28;
t21 = g(1) * t45 - g(2) * t47;
t86 = t62 * t65 - t64 * t69;
t85 = t131 * pkin(3);
t84 = g(1) * t125 + g(2) * t73;
t83 = pkin(3) * t56 + t45 * t77 - t46 * t63 + t98;
t9 = g(1) * t28 + g(2) * t24 + g(3) * t39;
t11 = g(1) * t29 + g(2) * t25 + g(3) * t40;
t82 = t46 * t71 + t75 * t99;
t81 = -g(1) * t47 - g(2) * t45 + g(3) * t113;
t80 = g(1) * t48 + g(2) * t46 + g(3) * t116;
t79 = t82 * pkin(3);
t31 = t48 * t75 + t101;
t15 = t29 * t74 + t121;
t13 = t81 * t64;
t8 = t9 * t74;
t7 = t9 * t70;
t6 = g(1) * t129 - g(2) * t15;
t5 = -g(1) * t130 - g(2) * t14;
t4 = -g(1) * (-t47 * t117 + t120) - g(2) * (-t45 * t117 + t122) - (t65 * t110 + t70 * t72) * t126;
t3 = -g(1) * (t47 * t118 + t48 * t74) - g(2) * (t45 * t118 + t46 * t74) - (-t65 * t112 + t72 * t74) * t126;
t2 = g(1) * t15 + g(2) * t129 - g(3) * (t68 * t112 - t40 * t74);
t10 = [0, 0, 0, 0, 0, 0, g(1) * t73 - g(2) * t125, t84, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t48, -t21, -t84 * t68, -g(1) * t98 - g(2) * t104, 0, 0, 0, 0, 0, 0, g(1) * t94 - g(2) * t31, -g(1) * t82 - g(2) * t30, t21, -g(1) * (-t46 * pkin(2) - t45 * pkin(9) + t98) - g(2) * (t48 * pkin(2) + t47 * pkin(9) + t104) 0, 0, 0, 0, 0, 0, g(1) * t25 - g(2) * t29, t88, t21, -g(1) * t83 - g(2) * t90, 0, 0, 0, 0, 0, 0, t6, t5, -t88, -g(1) * (-pkin(4) * t25 - pkin(11) * t24 + t83) - g(2) * (t29 * pkin(4) + t28 * pkin(11) + t90) 0, 0, 0, 0, 0, 0, t6, t5, -t88, -g(1) * (-pkin(5) * t124 + t24 * t69 - t25 * t62 + t83) - g(2) * (pkin(5) * t121 - t28 * t69 + t29 * t62 + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t80, 0, 0, 0, 0, 0, 0, 0, 0, -t81 * t75, t81 * t71, -t80, -g(1) * (-t47 * pkin(2) + t48 * pkin(9)) - g(2) * (-t45 * pkin(2) + t46 * pkin(9)) - (pkin(2) * t76 + pkin(9) * t72) * t126, 0, 0, 0, 0, 0, 0, -t81 * t65, t13, -t80, -g(1) * t105 - g(2) * t106 - g(3) * (-t68 * t111 + t49) 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * (-t89 * t47 + t105) - g(2) * (-t89 * t45 + t106) - t127 - (t89 * t76 - t111) * t126, 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * (pkin(5) * t120 - t86 * t47 + t105) - g(2) * (pkin(5) * t122 - t86 * t45 + t106) - t127 - (t86 * t76 + (pkin(5) * t70 - t77) * t72) * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t30 + g(2) * t82 - g(3) * t131, g(1) * t31 + g(2) * t94 - g(3) * (-t103 * t71 - t72 * t114) 0, 0, 0, 0, 0, 0, 0, 0, t9, t11, 0, -g(1) * t91 + g(2) * t79 - g(3) * t85, 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * (t91 + t96) - g(2) * (-t79 + t97) - g(3) * (t85 + t95) 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * (t91 + t108) - g(2) * (-t79 + t109) - g(3) * (t85 + t107); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t11, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * t96 - g(2) * t97 - g(3) * t95, 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * t108 - g(2) * t109 - g(3) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9;];
taug_reg  = t10;
