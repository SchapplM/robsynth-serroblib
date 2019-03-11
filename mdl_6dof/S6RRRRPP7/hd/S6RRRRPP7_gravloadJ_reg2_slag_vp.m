% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t120 = cos(pkin(6));
t147 = cos(qJ(1));
t108 = t120 * t147;
t146 = sin(qJ(1));
t87 = sin(qJ(2));
t90 = cos(qJ(2));
t61 = -t90 * t108 + t146 * t87;
t88 = cos(qJ(4));
t140 = t61 * t88;
t83 = sin(pkin(6));
t117 = t83 * t147;
t62 = t87 * t108 + t146 * t90;
t86 = sin(qJ(3));
t89 = cos(qJ(3));
t33 = -t86 * t117 + t62 * t89;
t85 = sin(qJ(4));
t145 = t33 * t85;
t151 = t140 - t145;
t84 = -qJ(5) - pkin(10);
t130 = t84 * t86;
t78 = t88 * pkin(4) + pkin(3);
t150 = -t78 * t89 + t130;
t141 = t61 * t85;
t149 = t33 * t88 + t141;
t82 = qJ(4) + pkin(11);
t79 = sin(t82);
t80 = cos(t82);
t9 = t33 * t79 - t61 * t80;
t10 = t33 * t80 + t61 * t79;
t148 = g(3) * t83;
t116 = t83 * t146;
t107 = t120 * t146;
t64 = -t87 * t107 + t147 * t90;
t37 = t86 * t116 + t64 * t89;
t144 = t37 * t85;
t139 = t62 * t85;
t63 = t90 * t107 + t147 * t87;
t138 = t63 * t85;
t137 = t63 * t88;
t136 = t64 * t85;
t134 = t79 * t89;
t133 = t80 * t89;
t132 = t83 * t87;
t131 = t83 * t90;
t129 = t85 * t87;
t128 = t85 * t89;
t127 = t88 * t89;
t126 = t89 * t90;
t32 = t89 * t117 + t62 * t86;
t125 = -t32 * t78 - t33 * t84;
t36 = -t89 * t116 + t64 * t86;
t124 = -t36 * t78 - t37 * t84;
t59 = -t120 * t89 + t86 * t132;
t60 = t120 * t86 + t89 * t132;
t123 = -t59 * t78 - t60 * t84;
t122 = pkin(2) * t131 + pkin(9) * t132;
t121 = t147 * pkin(1) + pkin(8) * t116;
t119 = t79 * t131;
t118 = t88 * t131;
t115 = -t61 * pkin(2) + t62 * pkin(9);
t114 = -t63 * pkin(2) + t64 * pkin(9);
t113 = -t146 * pkin(1) + pkin(8) * t117;
t112 = pkin(3) * t89 + pkin(10) * t86;
t13 = t37 * t79 - t63 * t80;
t111 = -g(1) * t9 + g(2) * t13;
t110 = -g(1) * t32 + g(2) * t36;
t109 = g(1) * t61 - g(2) * t63;
t106 = -pkin(5) * t80 - qJ(6) * t79;
t105 = t64 * pkin(2) + t63 * pkin(9) + t121;
t104 = -t60 * t85 - t118;
t26 = t80 * t131 + t60 * t79;
t1 = g(1) * t13 + g(2) * t9 + g(3) * t26;
t14 = t37 * t80 + t63 * t79;
t27 = t60 * t80 - t119;
t103 = g(1) * t14 + g(2) * t10 + g(3) * t27;
t18 = -t61 * t134 - t62 * t80;
t20 = -t63 * t134 - t64 * t80;
t38 = t89 * t119 - t80 * t132;
t102 = g(1) * t20 + g(2) * t18 + g(3) * t38;
t101 = g(1) * t36 + g(2) * t32 + g(3) * t59;
t100 = g(1) * t37 + g(2) * t33 + g(3) * t60;
t99 = g(1) * t147 + g(2) * t146;
t98 = -t62 * pkin(2) - t61 * pkin(9) + t113;
t97 = -g(1) * t63 - g(2) * t61 + g(3) * t131;
t96 = g(1) * t64 + g(2) * t62 + g(3) * t132;
t95 = pkin(4) * t138 - t36 * t84 + t37 * t78 + t105;
t94 = -t130 * t131 + t122 + (pkin(4) * t129 + t126 * t78) * t83;
t93 = pkin(4) * t139 + t150 * t61 + t115;
t92 = pkin(4) * t136 + t150 * t63 + t114;
t91 = -pkin(4) * t141 + t32 * t84 - t33 * t78 + t98;
t52 = pkin(4) * t137;
t48 = pkin(4) * t140;
t39 = (t80 * t126 + t79 * t87) * t83;
t21 = -t63 * t133 + t64 * t79;
t19 = -t61 * t133 + t62 * t79;
t17 = t97 * t86;
t16 = t37 * t88 + t138;
t15 = t137 - t144;
t5 = t101 * t80;
t4 = t101 * t79;
t3 = g(1) * t10 - g(2) * t14;
t2 = -g(1) * t21 - g(2) * t19 - g(3) * t39;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t146 - g(2) * t147, t99, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t62 - g(2) * t64, -t109, -t99 * t83, -g(1) * t113 - g(2) * t121, 0, 0, 0, 0, 0, 0, g(1) * t33 - g(2) * t37, t110, t109, -g(1) * t98 - g(2) * t105, 0, 0, 0, 0, 0, 0, g(1) * t149 - g(2) * t16, g(1) * t151 - g(2) * t15, -t110, -g(1) * (-pkin(3) * t33 - pkin(10) * t32 + t98) - g(2) * (t37 * pkin(3) + t36 * pkin(10) + t105) 0, 0, 0, 0, 0, 0, t3, t111, -t110, -g(1) * t91 - g(2) * t95, 0, 0, 0, 0, 0, 0, t3, -t110, -t111, -g(1) * (-pkin(5) * t10 - qJ(6) * t9 + t91) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t96, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t89, t17, -t96, -g(1) * t114 - g(2) * t115 - g(3) * t122, 0, 0, 0, 0, 0, 0, -g(1) * (-t63 * t127 + t136) - g(2) * (-t61 * t127 + t139) - (t88 * t126 + t129) * t148, -g(1) * (t63 * t128 + t64 * t88) - g(2) * (t61 * t128 + t62 * t88) - (-t85 * t126 + t87 * t88) * t148, -t17, -g(1) * (-t112 * t63 + t114) - g(2) * (-t112 * t61 + t115) - g(3) * (t112 * t131 + t122) 0, 0, 0, 0, 0, 0, t2, t102, -t17, -g(1) * t92 - g(2) * t93 - g(3) * t94, 0, 0, 0, 0, 0, 0, t2, -t17, -t102, -g(1) * (t21 * pkin(5) + t20 * qJ(6) + t92) - g(2) * (t19 * pkin(5) + t18 * qJ(6) + t93) - g(3) * (t39 * pkin(5) + t38 * qJ(6) + t94); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t100, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t88, -t101 * t85, -t100, -g(1) * (-t36 * pkin(3) + t37 * pkin(10)) - g(2) * (-t32 * pkin(3) + t33 * pkin(10)) - g(3) * (-t59 * pkin(3) + t60 * pkin(10)) 0, 0, 0, 0, 0, 0, t5, -t4, -t100, -g(1) * t124 - g(2) * t125 - g(3) * t123, 0, 0, 0, 0, 0, 0, t5, -t100, t4, -g(1) * (t106 * t36 + t124) - g(2) * (t106 * t32 + t125) - g(3) * (t106 * t59 + t123); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t151 - g(3) * t104, g(1) * t16 + g(2) * t149 - g(3) * (t85 * t131 - t60 * t88) 0, 0, 0, 0, 0, 0, 0, 0, t1, t103, 0, -g(1) * t52 - g(2) * t48 + (g(3) * t118 + t100 * t85) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t103, -g(1) * (-pkin(4) * t144 - t13 * pkin(5) + t14 * qJ(6) + t52) - g(2) * (-pkin(4) * t145 - t9 * pkin(5) + t10 * qJ(6) + t48) - g(3) * (t104 * pkin(4) - t26 * pkin(5) + t27 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
