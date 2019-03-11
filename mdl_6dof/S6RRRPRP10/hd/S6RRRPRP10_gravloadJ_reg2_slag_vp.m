% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t81 = cos(pkin(11));
t74 = t81 * pkin(4) + pkin(3);
t82 = -pkin(10) - qJ(4);
t83 = sin(qJ(3));
t85 = cos(qJ(3));
t139 = -t74 * t85 + t82 * t83;
t137 = cos(qJ(1));
t80 = sin(pkin(6));
t115 = t80 * t137;
t116 = cos(pkin(6));
t103 = t116 * t137;
t135 = sin(qJ(1));
t136 = cos(qJ(2));
t84 = sin(qJ(2));
t58 = t84 * t103 + t135 * t136;
t31 = -t83 * t115 + t58 * t85;
t57 = -t136 * t103 + t135 * t84;
t78 = pkin(11) + qJ(5);
t75 = sin(t78);
t76 = cos(t78);
t9 = t31 * t75 - t57 * t76;
t10 = t31 * t76 + t57 * t75;
t138 = g(3) * t80;
t79 = sin(pkin(11));
t132 = t57 * t79;
t131 = t58 * t79;
t102 = t116 * t135;
t59 = t136 * t102 + t137 * t84;
t130 = t59 * t79;
t60 = -t84 * t102 + t137 * t136;
t129 = t60 * t79;
t127 = t75 * t85;
t126 = t76 * t85;
t125 = t79 * t85;
t124 = t80 * t84;
t123 = t81 * t85;
t30 = t85 * t115 + t58 * t83;
t121 = -t30 * t74 - t31 * t82;
t113 = t80 * t135;
t34 = -t85 * t113 + t60 * t83;
t35 = t83 * t113 + t60 * t85;
t120 = -t34 * t74 - t35 * t82;
t55 = -t116 * t85 + t83 * t124;
t56 = t116 * t83 + t85 * t124;
t119 = -t55 * t74 - t56 * t82;
t114 = t80 * t136;
t118 = pkin(2) * t114 + pkin(9) * t124;
t117 = t137 * pkin(1) + pkin(8) * t113;
t112 = t83 * t136;
t111 = t85 * t136;
t110 = -t57 * pkin(2) + t58 * pkin(9);
t109 = -t59 * pkin(2) + t60 * pkin(9);
t108 = t75 * t114;
t107 = -t135 * pkin(1) + pkin(8) * t115;
t13 = t35 * t75 - t59 * t76;
t106 = -g(1) * t9 + g(2) * t13;
t105 = -g(1) * t30 + g(2) * t34;
t104 = g(1) * t57 - g(2) * t59;
t101 = -pkin(3) * t85 - qJ(4) * t83;
t100 = -pkin(5) * t76 - qJ(6) * t75;
t99 = t60 * pkin(2) + t59 * pkin(9) + t117;
t24 = t76 * t114 + t56 * t75;
t1 = g(1) * t13 + g(2) * t9 + g(3) * t24;
t14 = t35 * t76 + t59 * t75;
t25 = t56 * t76 - t108;
t98 = g(1) * t14 + g(2) * t10 + g(3) * t25;
t16 = -t57 * t127 - t58 * t76;
t18 = -t59 * t127 - t60 * t76;
t36 = t85 * t108 - t76 * t124;
t97 = g(1) * t18 + g(2) * t16 + g(3) * t36;
t96 = g(1) * t34 + g(2) * t30 + g(3) * t55;
t95 = g(1) * t35 + g(2) * t31 + g(3) * t56;
t94 = g(1) * t137 + g(2) * t135;
t93 = -t58 * pkin(2) - t57 * pkin(9) + t107;
t92 = g(1) * t60 + g(2) * t58 + g(3) * t124;
t91 = pkin(4) * t130 - t34 * t82 + t35 * t74 + t99;
t90 = pkin(4) * t131 + t139 * t57 + t110;
t89 = pkin(4) * t129 + t139 * t59 + t109;
t88 = -g(1) * t59 - g(2) * t57 + g(3) * t114;
t87 = t79 * pkin(4) * t124 + t118 + (t111 * t74 - t112 * t82) * t80;
t86 = -pkin(4) * t132 + t30 * t82 - t31 * t74 + t93;
t37 = (t76 * t111 + t75 * t84) * t80;
t19 = -t59 * t126 + t60 * t75;
t17 = -t57 * t126 + t58 * t75;
t15 = t88 * t83;
t5 = t96 * t76;
t4 = t96 * t75;
t3 = g(1) * t10 - g(2) * t14;
t2 = -g(1) * t19 - g(2) * t17 - g(3) * t37;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t135 - g(2) * t137, t94, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t60, -t104, -t94 * t80, -g(1) * t107 - g(2) * t117, 0, 0, 0, 0, 0, 0, g(1) * t31 - g(2) * t35, t105, t104, -g(1) * t93 - g(2) * t99, 0, 0, 0, 0, 0, 0, -g(1) * (-t31 * t81 - t132) - g(2) * (t35 * t81 + t130) -g(1) * (t31 * t79 - t57 * t81) - g(2) * (-t35 * t79 + t59 * t81) -t105, -g(1) * (-pkin(3) * t31 - qJ(4) * t30 + t93) - g(2) * (t35 * pkin(3) + t34 * qJ(4) + t99) 0, 0, 0, 0, 0, 0, t3, t106, -t105, -g(1) * t86 - g(2) * t91, 0, 0, 0, 0, 0, 0, t3, -t105, -t106, -g(1) * (-pkin(5) * t10 - qJ(6) * t9 + t86) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t92, 0, 0, 0, 0, 0, 0, 0, 0, -t88 * t85, t15, -t92, -g(1) * t109 - g(2) * t110 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (-t59 * t123 + t129) - g(2) * (-t57 * t123 + t131) - (t81 * t111 + t79 * t84) * t138, -g(1) * (t59 * t125 + t60 * t81) - g(2) * (t57 * t125 + t58 * t81) - (-t79 * t111 + t81 * t84) * t138, -t15, -g(1) * (t101 * t59 + t109) - g(2) * (t101 * t57 + t110) - g(3) * ((pkin(3) * t111 + qJ(4) * t112) * t80 + t118) 0, 0, 0, 0, 0, 0, t2, t97, -t15, -g(1) * t89 - g(2) * t90 - g(3) * t87, 0, 0, 0, 0, 0, 0, t2, -t15, -t97, -g(1) * (t19 * pkin(5) + t18 * qJ(6) + t89) - g(2) * (t17 * pkin(5) + t16 * qJ(6) + t90) - g(3) * (t37 * pkin(5) + t36 * qJ(6) + t87); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t95, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t81, -t96 * t79, -t95, -g(1) * (-t34 * pkin(3) + t35 * qJ(4)) - g(2) * (-t30 * pkin(3) + t31 * qJ(4)) - g(3) * (-t55 * pkin(3) + t56 * qJ(4)) 0, 0, 0, 0, 0, 0, t5, -t4, -t95, -g(1) * t120 - g(2) * t121 - g(3) * t119, 0, 0, 0, 0, 0, 0, t5, -t95, t4, -g(1) * (t100 * t34 + t120) - g(2) * (t100 * t30 + t121) - g(3) * (t100 * t55 + t119); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t98, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t98, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t24 * pkin(5) + t25 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
