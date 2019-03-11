% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t106 = cos(pkin(6));
t66 = sin(pkin(6));
t69 = sin(qJ(2));
t116 = t66 * t69;
t68 = sin(qJ(3));
t72 = cos(qJ(3));
t132 = t106 * t72 - t68 * t116;
t114 = t66 * t72;
t124 = cos(qJ(1));
t73 = cos(qJ(2));
t70 = sin(qJ(1));
t99 = t70 * t106;
t44 = t124 * t73 - t69 * t99;
t24 = t70 * t114 - t44 * t68;
t102 = t66 * t124;
t87 = t106 * t124;
t42 = t69 * t87 + t70 * t73;
t65 = qJ(3) + qJ(4);
t62 = sin(t65);
t63 = cos(t65);
t18 = t63 * t102 + t42 * t62;
t41 = t70 * t69 - t73 * t87;
t67 = sin(qJ(6));
t71 = cos(qJ(6));
t131 = t18 * t67 + t41 * t71;
t130 = t18 * t71 - t41 * t67;
t107 = qJ(5) * t62;
t113 = t66 * t73;
t129 = (pkin(4) * t63 + t107) * t113;
t128 = g(3) * t66;
t127 = t18 * pkin(11);
t115 = t66 * t70;
t22 = -t63 * t115 + t44 * t62;
t126 = t22 * pkin(11);
t35 = -t106 * t63 + t62 * t116;
t125 = t35 * pkin(11);
t123 = t41 * t63;
t43 = t124 * t69 + t73 * t99;
t120 = t43 * t63;
t118 = t62 * t67;
t117 = t62 * t71;
t112 = t67 * t73;
t111 = t71 * t73;
t61 = t72 * pkin(3) + pkin(2);
t74 = -pkin(10) - pkin(9);
t110 = -t41 * t61 - t42 * t74;
t109 = -t43 * t61 - t44 * t74;
t108 = t124 * pkin(1) + pkin(8) * t115;
t104 = t68 * t115;
t101 = -t70 * pkin(1) + pkin(8) * t102;
t19 = -t62 * t102 + t42 * t63;
t55 = t68 * t102;
t100 = t42 * t72 - t55;
t97 = -t18 * pkin(4) + t19 * qJ(5);
t23 = t62 * t115 + t44 * t63;
t96 = -t22 * pkin(4) + t23 * qJ(5);
t36 = t106 * t62 + t63 * t116;
t95 = -t35 * pkin(4) + t36 * qJ(5);
t94 = -pkin(4) * t123 - t41 * t107 + t110;
t93 = -pkin(4) * t120 - t43 * t107 + t109;
t92 = t24 * pkin(3);
t46 = t61 * t113;
t91 = -t74 * t116 + t46;
t90 = pkin(3) * t104 - t43 * t74 + t44 * t61 + t108;
t89 = -g(1) * t18 + g(2) * t22;
t88 = -g(1) * t19 + g(2) * t23;
t13 = g(1) * t41 - g(2) * t43;
t86 = t132 * pkin(3);
t85 = g(1) * t124 + g(2) * t70;
t84 = pkin(3) * t55 + t41 * t74 - t42 * t61 + t101;
t4 = g(1) * t22 + g(2) * t18 + g(3) * t35;
t6 = g(1) * t23 + g(2) * t19 + g(3) * t36;
t83 = t72 * t102 + t42 * t68;
t82 = t92 + t96;
t81 = -g(1) * t43 - g(2) * t41 + g(3) * t113;
t80 = g(1) * t44 + g(2) * t42 + g(3) * t116;
t79 = t23 * pkin(4) + t22 * qJ(5) + t90;
t78 = t83 * pkin(3);
t77 = t86 + t95;
t76 = -pkin(4) * t19 - qJ(5) * t18 + t84;
t75 = -t78 + t97;
t25 = t44 * t72 + t104;
t11 = t22 * t67 + t43 * t71;
t10 = t22 * t71 - t43 * t67;
t9 = t81 * t63;
t8 = t81 * t62;
t2 = t6 * t71;
t1 = t6 * t67;
t3 = [0, 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t124, t85, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t44, -t13, -t85 * t66, -g(1) * t101 - g(2) * t108, 0, 0, 0, 0, 0, 0, g(1) * t100 - g(2) * t25, -g(1) * t83 - g(2) * t24, t13, -g(1) * (-t42 * pkin(2) - t41 * pkin(9) + t101) - g(2) * (t44 * pkin(2) + t43 * pkin(9) + t108) 0, 0, 0, 0, 0, 0, -t88, t89, t13, -g(1) * t84 - g(2) * t90, 0, 0, 0, 0, 0, 0, t13, t88, -t89, -g(1) * t76 - g(2) * t79, 0, 0, 0, 0, 0, 0, g(1) * t131 - g(2) * t11, g(1) * t130 - g(2) * t10, -t88, -g(1) * (-t41 * pkin(5) - pkin(11) * t19 + t76) - g(2) * (t43 * pkin(5) + t23 * pkin(11) + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t80, 0, 0, 0, 0, 0, 0, 0, 0, -t81 * t72, t81 * t68, -t80, -g(1) * (-t43 * pkin(2) + t44 * pkin(9)) - g(2) * (-t41 * pkin(2) + t42 * pkin(9)) - (pkin(2) * t73 + pkin(9) * t69) * t128, 0, 0, 0, 0, 0, 0, -t9, t8, -t80, -g(1) * t109 - g(2) * t110 - g(3) * t91, 0, 0, 0, 0, 0, 0, -t80, t9, -t8, -g(1) * t93 - g(2) * t94 - g(3) * (t91 + t129) 0, 0, 0, 0, 0, 0, -g(1) * (-t43 * t118 + t44 * t71) - g(2) * (-t41 * t118 + t42 * t71) - (t62 * t112 + t69 * t71) * t128, -g(1) * (-t43 * t117 - t44 * t67) - g(2) * (-t41 * t117 - t42 * t67) - (t62 * t111 - t67 * t69) * t128, -t9, -g(1) * (t44 * pkin(5) - pkin(11) * t120 + t93) - g(2) * (t42 * pkin(5) - pkin(11) * t123 + t94) - g(3) * (t46 + t129) - (pkin(11) * t63 * t73 + (pkin(5) - t74) * t69) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t24 + g(2) * t83 - g(3) * t132, g(1) * t25 + g(2) * t100 - g(3) * (-t106 * t68 - t69 * t114) 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, -g(1) * t92 + g(2) * t78 - g(3) * t86, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, -g(1) * t82 - g(2) * t75 - g(3) * t77, 0, 0, 0, 0, 0, 0, -t1, -t2, t4, -g(1) * (t82 - t126) - g(2) * (t75 - t127) - g(3) * (t77 - t125); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, -g(1) * t96 - g(2) * t97 - g(3) * t95, 0, 0, 0, 0, 0, 0, -t1, -t2, t4, -g(1) * (t96 - t126) - g(2) * (t97 - t127) - g(3) * (t95 - t125); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t130 - g(3) * (t66 * t112 + t35 * t71) g(1) * t11 + g(2) * t131 - g(3) * (t66 * t111 - t35 * t67) 0, 0;];
taug_reg  = t3;
