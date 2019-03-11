% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP10
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t127 = cos(pkin(6));
t154 = cos(qJ(1));
t111 = t127 * t154;
t153 = sin(qJ(1));
t90 = sin(qJ(2));
t93 = cos(qJ(2));
t65 = -t93 * t111 + t153 * t90;
t91 = cos(qJ(4));
t147 = t65 * t91;
t87 = sin(pkin(6));
t124 = t87 * t154;
t66 = t90 * t111 + t153 * t93;
t89 = sin(qJ(3));
t92 = cos(qJ(3));
t37 = -t89 * t124 + t66 * t92;
t88 = sin(qJ(4));
t152 = t37 * t88;
t158 = t147 - t152;
t94 = -pkin(11) - pkin(10);
t135 = t89 * t94;
t82 = t91 * pkin(4) + pkin(3);
t157 = -t82 * t92 + t135;
t148 = t65 * t88;
t156 = t37 * t91 + t148;
t86 = qJ(4) + qJ(5);
t83 = sin(t86);
t84 = cos(t86);
t12 = t37 * t83 - t65 * t84;
t13 = t37 * t84 + t65 * t83;
t155 = g(3) * t87;
t123 = t87 * t153;
t110 = t127 * t153;
t68 = -t90 * t110 + t154 * t93;
t41 = t89 * t123 + t68 * t92;
t151 = t41 * t88;
t146 = t66 * t88;
t67 = t93 * t110 + t154 * t90;
t145 = t67 * t88;
t144 = t67 * t91;
t143 = t68 * t88;
t141 = t83 * t92;
t140 = t84 * t92;
t139 = t87 * t90;
t138 = t87 * t93;
t137 = t88 * t90;
t136 = t88 * t92;
t134 = t91 * t92;
t133 = t92 * t93;
t120 = -t92 * t124 - t66 * t89;
t132 = t120 * t82 - t37 * t94;
t40 = -t92 * t123 + t68 * t89;
t131 = -t40 * t82 - t41 * t94;
t63 = t127 * t92 - t89 * t139;
t64 = t127 * t89 + t92 * t139;
t130 = t63 * t82 - t64 * t94;
t129 = pkin(2) * t138 + pkin(9) * t139;
t128 = t154 * pkin(1) + pkin(8) * t123;
t126 = t83 * t138;
t125 = t91 * t138;
t122 = -t65 * pkin(2) + t66 * pkin(9);
t121 = -t67 * pkin(2) + t68 * pkin(9);
t119 = -t12 * pkin(5) + t13 * qJ(6);
t16 = t41 * t83 - t67 * t84;
t17 = t41 * t84 + t67 * t83;
t118 = -t16 * pkin(5) + t17 * qJ(6);
t30 = t84 * t138 + t64 * t83;
t31 = t64 * t84 - t126;
t117 = -t30 * pkin(5) + t31 * qJ(6);
t116 = -t153 * pkin(1) + pkin(8) * t124;
t115 = pkin(3) * t92 + pkin(10) * t89;
t114 = -g(1) * t12 + g(2) * t16;
t113 = g(1) * t120 + g(2) * t40;
t112 = g(1) * t65 - g(2) * t67;
t109 = pkin(5) * t84 + qJ(6) * t83;
t108 = t68 * pkin(2) + t67 * pkin(9) + t128;
t107 = -t64 * t88 - t125;
t1 = g(1) * t16 + g(2) * t12 + g(3) * t30;
t3 = g(1) * t17 + g(2) * t13 + g(3) * t31;
t21 = -t65 * t141 - t66 * t84;
t23 = -t67 * t141 - t68 * t84;
t42 = t92 * t126 - t84 * t139;
t106 = g(1) * t23 + g(2) * t21 + g(3) * t42;
t105 = g(1) * t40 - g(2) * t120 - g(3) * t63;
t104 = g(1) * t41 + g(2) * t37 + g(3) * t64;
t103 = g(1) * t154 + g(2) * t153;
t102 = -t66 * pkin(2) - t65 * pkin(9) + t116;
t101 = -g(1) * t67 - g(2) * t65 + g(3) * t138;
t100 = g(1) * t68 + g(2) * t66 + g(3) * t139;
t99 = pkin(4) * t145 - t40 * t94 + t41 * t82 + t108;
t98 = -t135 * t138 + t129 + (pkin(4) * t137 + t133 * t82) * t87;
t97 = pkin(4) * t146 + t157 * t65 + t122;
t96 = pkin(4) * t143 + t157 * t67 + t121;
t95 = -pkin(4) * t148 - t120 * t94 - t37 * t82 + t102;
t56 = pkin(4) * t144;
t52 = pkin(4) * t147;
t43 = (t84 * t133 + t83 * t90) * t87;
t24 = -t67 * t140 + t68 * t83;
t22 = -t65 * t140 + t66 * t83;
t20 = t101 * t89;
t19 = t41 * t91 + t145;
t18 = t144 - t151;
t7 = t105 * t84;
t6 = t105 * t83;
t5 = g(1) * t13 - g(2) * t17;
t4 = -g(1) * t24 - g(2) * t22 - g(3) * t43;
t2 = [0, 0, 0, 0, 0, 0, g(1) * t153 - g(2) * t154, t103, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t66 - g(2) * t68, -t112, -t103 * t87, -g(1) * t116 - g(2) * t128, 0, 0, 0, 0, 0, 0, g(1) * t37 - g(2) * t41, t113, t112, -g(1) * t102 - g(2) * t108, 0, 0, 0, 0, 0, 0, g(1) * t156 - g(2) * t19, g(1) * t158 - g(2) * t18, -t113, -g(1) * (-pkin(3) * t37 + pkin(10) * t120 + t102) - g(2) * (t41 * pkin(3) + t40 * pkin(10) + t108) 0, 0, 0, 0, 0, 0, t5, t114, -t113, -g(1) * t95 - g(2) * t99, 0, 0, 0, 0, 0, 0, t5, -t113, -t114, -g(1) * (-pkin(5) * t13 - qJ(6) * t12 + t95) - g(2) * (t17 * pkin(5) + t16 * qJ(6) + t99); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, t100, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t92, t20, -t100, -g(1) * t121 - g(2) * t122 - g(3) * t129, 0, 0, 0, 0, 0, 0, -g(1) * (-t67 * t134 + t143) - g(2) * (-t65 * t134 + t146) - (t91 * t133 + t137) * t155, -g(1) * (t67 * t136 + t68 * t91) - g(2) * (t65 * t136 + t66 * t91) - (-t88 * t133 + t90 * t91) * t155, -t20, -g(1) * (-t115 * t67 + t121) - g(2) * (-t115 * t65 + t122) - g(3) * (t115 * t138 + t129) 0, 0, 0, 0, 0, 0, t4, t106, -t20, -g(1) * t96 - g(2) * t97 - g(3) * t98, 0, 0, 0, 0, 0, 0, t4, -t20, -t106, -g(1) * (t24 * pkin(5) + t23 * qJ(6) + t96) - g(2) * (t22 * pkin(5) + t21 * qJ(6) + t97) - g(3) * (t43 * pkin(5) + t42 * qJ(6) + t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t104, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t91, -t105 * t88, -t104, -g(1) * (-t40 * pkin(3) + t41 * pkin(10)) - g(2) * (pkin(3) * t120 + t37 * pkin(10)) - g(3) * (t63 * pkin(3) + t64 * pkin(10)) 0, 0, 0, 0, 0, 0, t7, -t6, -t104, -g(1) * t131 - g(2) * t132 - g(3) * t130, 0, 0, 0, 0, 0, 0, t7, -t104, t6, -g(1) * (-t109 * t40 + t131) - g(2) * (t109 * t120 + t132) - g(3) * (t109 * t63 + t130); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t158 - g(3) * t107, g(1) * t19 + g(2) * t156 - g(3) * (t88 * t138 - t64 * t91) 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t56 - g(2) * t52 + (g(3) * t125 + t104 * t88) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t151 + t118 + t56) - g(2) * (-pkin(4) * t152 + t119 + t52) - g(3) * (pkin(4) * t107 + t117); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t118 - g(2) * t119 - g(3) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
