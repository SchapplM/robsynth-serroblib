% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 18:43:27
% EndTime: 2019-05-08 18:43:44
% DurationCPUTime: 2.30s
% Computational Cost: add. (1508->172), mult. (4402->328), div. (0->0), fcn. (5838->18), ass. (0->110)
t73 = sin(pkin(6));
t85 = cos(qJ(1));
t114 = t73 * t85;
t124 = sin(qJ(1));
t80 = sin(qJ(2));
t84 = cos(qJ(2));
t104 = cos(pkin(6));
t97 = t85 * t104;
t64 = t124 * t80 - t84 * t97;
t72 = sin(pkin(7));
t75 = cos(pkin(7));
t60 = t75 * t114 - t64 * t72;
t71 = sin(pkin(8));
t123 = t60 * t71;
t125 = cos(qJ(4));
t102 = t72 * t114;
t65 = t124 * t84 + t80 * t97;
t79 = sin(qJ(3));
t83 = cos(qJ(3));
t46 = (t64 * t75 + t102) * t83 + t65 * t79;
t112 = t75 * t79;
t47 = t79 * t102 + t64 * t112 - t65 * t83;
t74 = cos(pkin(8));
t78 = sin(qJ(4));
t17 = t47 * t125 + (t46 * t74 + t123) * t78;
t32 = t46 * t71 - t60 * t74;
t77 = sin(qJ(5));
t82 = cos(qJ(5));
t5 = t17 * t82 - t32 * t77;
t76 = sin(qJ(6));
t134 = t5 * t76;
t81 = cos(qJ(6));
t133 = t5 * t81;
t132 = t17 * t77 + t32 * t82;
t99 = t74 * t125;
t128 = -t46 * t99 + t47 * t78;
t120 = t71 * t72;
t119 = t71 * t77;
t118 = t71 * t82;
t117 = t72 * t74;
t116 = t73 * t80;
t115 = t73 * t84;
t113 = t74 * t78;
t111 = t75 * t83;
t110 = t76 * t82;
t109 = t79 * t80;
t108 = t79 * t84;
t107 = t80 * t83;
t106 = t81 * t82;
t105 = t83 * t84;
t103 = t72 * t116;
t101 = t71 * t125;
t100 = t73 * t124;
t98 = t72 * t104;
t96 = t72 * t101;
t95 = t104 * t124;
t58 = t83 * t98 + (t75 * t105 - t109) * t73;
t59 = t79 * t98 + (t75 * t108 + t107) * t73;
t88 = t104 * t75 - t72 * t115;
t86 = t88 * t71;
t31 = t59 * t125 + (t58 * t74 + t86) * t78;
t43 = -t58 * t71 + t88 * t74;
t67 = -t80 * t95 + t85 * t84;
t66 = -t85 * t80 - t84 * t95;
t90 = t72 * t100 + t66 * t75;
t48 = -t67 * t79 + t90 * t83;
t49 = t67 * t83 + t90 * t79;
t91 = t75 * t100 - t66 * t72;
t87 = t91 * t71;
t19 = t49 * t125 + (t48 * t74 + t87) * t78;
t34 = -t48 * t71 + t91 * t74;
t6 = -t19 * t77 + t34 * t82;
t93 = g(1) * t6 + g(2) * t132 + g(3) * (-t31 * t77 + t43 * t82);
t14 = t123 * t125 - t128;
t18 = -t125 * t87 - t48 * t99 + t49 * t78;
t30 = -t125 * t86 - t58 * t99 + t59 * t78;
t92 = g(1) * t18 + g(2) * t14 + g(3) * t30;
t63 = (-t75 * t109 + t105) * t73;
t62 = (-t75 * t107 - t108) * t73;
t54 = t74 * t103 - t62 * t71;
t53 = -t67 * t112 + t66 * t83;
t52 = -t67 * t111 - t66 * t79;
t51 = -t65 * t112 - t64 * t83;
t50 = -t65 * t111 + t64 * t79;
t40 = t67 * t117 - t52 * t71;
t39 = t65 * t117 - t50 * t71;
t38 = t63 * t125 + (t71 * t103 + t62 * t74) * t78;
t37 = -t96 * t116 - t62 * t99 + t63 * t78;
t36 = -t59 * t113 + t58 * t125;
t35 = t58 * t78 + t59 * t99;
t29 = t59 * t119 + t36 * t82;
t28 = t38 * t82 + t54 * t77;
t27 = -t49 * t113 + t48 * t125;
t26 = t48 * t78 + t49 * t99;
t25 = t113 * t47 - t125 * t46;
t24 = -t46 * t78 - t47 * t99;
t23 = t53 * t125 + (t67 * t120 + t52 * t74) * t78;
t22 = -t52 * t99 + t53 * t78 - t67 * t96;
t21 = t51 * t125 + (t65 * t120 + t50 * t74) * t78;
t20 = -t50 * t99 + t51 * t78 - t65 * t96;
t16 = -t60 * t101 + t128;
t13 = t31 * t82 + t43 * t77;
t11 = t49 * t119 + t27 * t82;
t10 = -t119 * t47 + t25 * t82;
t9 = t23 * t82 + t40 * t77;
t8 = t21 * t82 + t39 * t77;
t7 = t19 * t82 + t34 * t77;
t2 = t18 * t76 + t7 * t81;
t1 = t18 * t81 - t7 * t76;
t3 = [0, g(1) * t124 - g(2) * t85, g(1) * t85 + g(2) * t124, 0, 0, 0, 0, 0, g(1) * t65 - g(2) * t67, -g(1) * t64 - g(2) * t66, 0, 0, 0, 0, 0, -g(1) * t47 - g(2) * t49, -g(1) * t46 - g(2) * t48, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, g(1) * t16 + g(2) * t18, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, g(1) * t132 - g(2) * t6, 0, 0, 0, 0, 0, -g(1) * (t16 * t76 + t133) - g(2) * t2, -g(1) * (t16 * t81 - t134) - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t66 + g(2) * t64 - g(3) * t115, g(1) * t67 + g(2) * t65 + g(3) * t116, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t51 - g(3) * t63, -g(1) * t52 - g(2) * t50 - g(3) * t62, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t21 - g(3) * t38, g(1) * t22 + g(2) * t20 + g(3) * t37, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t28, -g(1) * (-t23 * t77 + t40 * t82) - g(2) * (-t21 * t77 + t39 * t82) - g(3) * (-t38 * t77 + t54 * t82) 0, 0, 0, 0, 0, -g(1) * (t22 * t76 + t81 * t9) - g(2) * (t20 * t76 + t8 * t81) - g(3) * (t28 * t81 + t37 * t76) -g(1) * (t22 * t81 - t76 * t9) - g(2) * (t20 * t81 - t76 * t8) - g(3) * (-t28 * t76 + t37 * t81); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t48 + g(2) * t46 - g(3) * t58, g(1) * t49 - g(2) * t47 + g(3) * t59, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t25 - g(3) * t36, g(1) * t26 + g(2) * t24 + g(3) * t35, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t10 - g(3) * t29, -g(1) * (t49 * t118 - t27 * t77) - g(2) * (-t118 * t47 - t25 * t77) - g(3) * (t59 * t118 - t36 * t77) 0, 0, 0, 0, 0, -g(1) * (t11 * t81 + t26 * t76) - g(2) * (t10 * t81 + t24 * t76) - g(3) * (t29 * t81 + t35 * t76) -g(1) * (-t11 * t76 + t26 * t81) - g(2) * (-t10 * t76 + t24 * t81) - g(3) * (-t29 * t76 + t35 * t81); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, g(1) * t19 - g(2) * t17 + g(3) * t31, 0, 0, 0, 0, 0, t92 * t82, -t92 * t77, 0, 0, 0, 0, 0, -g(1) * (-t18 * t106 + t19 * t76) - g(2) * (-t14 * t106 - t17 * t76) - g(3) * (-t30 * t106 + t31 * t76) -g(1) * (t18 * t110 + t19 * t81) - g(2) * (t14 * t110 - t17 * t81) - g(3) * (t30 * t110 + t31 * t81); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, g(1) * t7 - g(2) * t5 + g(3) * t13, 0, 0, 0, 0, 0, -t93 * t81, t93 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * (t14 * t81 + t134) - g(3) * (-t13 * t76 + t30 * t81) g(1) * t2 - g(2) * (-t14 * t76 + t133) - g(3) * (-t13 * t81 - t30 * t76);];
taug_reg  = t3;
