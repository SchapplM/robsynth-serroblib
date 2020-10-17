% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:44:23
% EndTime: 2019-05-05 05:44:26
% DurationCPUTime: 0.88s
% Computational Cost: add. (951->172), mult. (2313->274), div. (0->0), fcn. (2947->16), ass. (0->94)
t107 = cos(pkin(12));
t67 = sin(pkin(6));
t100 = t67 * t107;
t108 = cos(pkin(7));
t66 = sin(pkin(7));
t106 = sin(pkin(12));
t120 = sin(qJ(2));
t122 = cos(qJ(2));
t109 = cos(pkin(6));
t86 = t109 * t107;
t75 = t106 * t120 - t122 * t86;
t128 = t66 * t100 + t75 * t108;
t85 = t109 * t106;
t76 = t107 * t120 + t122 * t85;
t99 = t67 * t106;
t127 = t76 * t108 - t66 * t99;
t50 = t106 * t122 + t120 * t86;
t51 = t107 * t122 - t120 * t85;
t126 = -g(1) * t51 - g(2) * t50;
t125 = pkin(9) * t66;
t121 = cos(qJ(3));
t64 = pkin(13) + qJ(5);
t62 = sin(t64);
t119 = t62 * t66;
t63 = cos(t64);
t118 = t63 * t66;
t70 = sin(qJ(6));
t117 = t63 * t70;
t72 = cos(qJ(6));
t116 = t63 * t72;
t65 = sin(pkin(13));
t115 = t65 * t66;
t68 = cos(pkin(13));
t114 = t66 * t68;
t71 = sin(qJ(3));
t19 = t128 * t121 + t50 * t71;
t20 = t50 * t121 - t128 * t71;
t61 = t68 * pkin(4) + pkin(3);
t69 = -pkin(10) - qJ(4);
t113 = -t19 * t61 - t20 * t69;
t21 = t127 * t121 + t51 * t71;
t22 = t51 * t121 - t127 * t71;
t112 = -t21 * t61 - t22 * t69;
t101 = t67 * t120;
t102 = t67 * t122;
t89 = t108 * t121;
t97 = t109 * t66;
t35 = t71 * t101 - t89 * t102 - t121 * t97;
t98 = t71 * t108;
t36 = t71 * t97 + (t120 * t121 + t122 * t98) * t67;
t111 = -t35 * t61 - t36 * t69;
t96 = t66 * t101;
t110 = pkin(2) * t102 + pkin(9) * t96;
t27 = t50 * t89 - t75 * t71;
t28 = -t75 * t121 - t50 * t98;
t47 = t75 * pkin(2);
t105 = -t27 * t69 + t28 * t61 - t47;
t29 = t51 * t89 - t76 * t71;
t30 = -t76 * t121 - t51 * t98;
t48 = t76 * pkin(2);
t104 = -t29 * t69 + t30 * t61 - t48;
t95 = t50 * t125 - t47;
t94 = t51 * t125 - t48;
t88 = t108 * t120;
t45 = (t121 * t88 + t122 * t71) * t67;
t46 = (t122 * t121 - t71 * t88) * t67;
t87 = t65 * t96;
t91 = pkin(4) * t87 - t45 * t69 + t46 * t61 + t110;
t90 = -pkin(5) * t63 - pkin(11) * t62;
t49 = -t66 * t102 + t109 * t108;
t15 = -t36 * t62 + t49 * t63;
t37 = -t108 * t100 + t75 * t66;
t5 = -t20 * t62 + t37 * t63;
t38 = t108 * t99 + t76 * t66;
t7 = -t22 * t62 + t38 * t63;
t84 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t16 = t36 * t63 + t49 * t62;
t6 = t20 * t63 + t37 * t62;
t8 = t22 * t63 + t38 * t62;
t83 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t11 = -t51 * t118 + t30 * t62;
t31 = t46 * t62 - t63 * t96;
t9 = -t50 * t118 + t28 * t62;
t82 = g(1) * t11 + g(2) * t9 + g(3) * t31;
t81 = g(1) * t21 + g(2) * t19 + g(3) * t35;
t80 = g(1) * t22 + g(2) * t20 + g(3) * t36;
t79 = g(1) * t29 + g(2) * t27 + g(3) * t45;
t78 = -g(3) * t101 + t126;
t77 = t126 * (pkin(4) * t65 + pkin(9)) * t66;
t32 = t46 * t63 + t62 * t96;
t12 = t51 * t119 + t30 * t63;
t10 = t50 * t119 + t28 * t63;
t1 = t81 * t62;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t76 + g(2) * t75 - g(3) * t102, -t78, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t28 - g(3) * t46, t79, t78 * t66, -g(1) * t94 - g(2) * t95 - g(3) * t110, 0, 0, 0, 0, 0, 0, -g(1) * (t51 * t115 + t30 * t68) - g(2) * (t50 * t115 + t28 * t68) - g(3) * (t46 * t68 + t87) -g(1) * (t51 * t114 - t30 * t65) - g(2) * (t50 * t114 - t28 * t65) - g(3) * (-t46 * t65 + t68 * t96) -t79, -g(1) * (t30 * pkin(3) + t29 * qJ(4) + t94) - g(2) * (t28 * pkin(3) + t27 * qJ(4) + t95) - g(3) * (t46 * pkin(3) + t45 * qJ(4) + t110) 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10 - g(3) * t32, t82, -t79, -g(1) * t104 - g(2) * t105 - g(3) * t91 + t77, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t72 + t29 * t70) - g(2) * (t10 * t72 + t27 * t70) - g(3) * (t32 * t72 + t45 * t70) -g(1) * (-t12 * t70 + t29 * t72) - g(2) * (-t10 * t70 + t27 * t72) - g(3) * (-t32 * t70 + t45 * t72) -t82, -g(1) * (t12 * pkin(5) + t11 * pkin(11) + t104) - g(2) * (t10 * pkin(5) + t9 * pkin(11) + t105) - g(3) * (t32 * pkin(5) + t31 * pkin(11) + t91) + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t80, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t68, -t81 * t65, -t80, -g(1) * (-t21 * pkin(3) + t22 * qJ(4)) - g(2) * (-t19 * pkin(3) + t20 * qJ(4)) - g(3) * (-t35 * pkin(3) + t36 * qJ(4)) 0, 0, 0, 0, 0, 0, t81 * t63, -t1, -t80, -g(1) * t112 - g(2) * t113 - g(3) * t111, 0, 0, 0, 0, 0, 0, -g(1) * (-t21 * t116 + t22 * t70) - g(2) * (-t19 * t116 + t20 * t70) - g(3) * (-t35 * t116 + t36 * t70) -g(1) * (t21 * t117 + t22 * t72) - g(2) * (t19 * t117 + t20 * t72) - g(3) * (t35 * t117 + t36 * t72) t1, -g(1) * (t90 * t21 + t112) - g(2) * (t90 * t19 + t113) - g(3) * (t90 * t35 + t111); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t83, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t72, t84 * t70, -t83, -g(1) * (t7 * pkin(5) + t8 * pkin(11)) - g(2) * (t5 * pkin(5) + t6 * pkin(11)) - g(3) * (t15 * pkin(5) + t16 * pkin(11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t21 * t72 - t8 * t70) - g(2) * (t19 * t72 - t6 * t70) - g(3) * (-t16 * t70 + t35 * t72) -g(1) * (-t21 * t70 - t8 * t72) - g(2) * (-t19 * t70 - t6 * t72) - g(3) * (-t16 * t72 - t35 * t70) 0, 0;];
taug_reg  = t2;
