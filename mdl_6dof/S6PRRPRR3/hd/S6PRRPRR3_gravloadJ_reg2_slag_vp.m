% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR3
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:50:19
% EndTime: 2019-05-05 04:50:23
% DurationCPUTime: 1.19s
% Computational Cost: add. (1117->181), mult. (3060->302), div. (0->0), fcn. (3967->16), ass. (0->102)
t110 = cos(pkin(13));
t68 = sin(pkin(13));
t77 = sin(qJ(3));
t81 = cos(qJ(3));
t59 = -t110 * t77 - t68 * t81;
t70 = sin(pkin(7));
t49 = t59 * t70;
t73 = cos(pkin(7));
t51 = t59 * t73;
t71 = sin(pkin(6));
t74 = cos(pkin(6));
t78 = sin(qJ(2));
t82 = cos(qJ(2));
t95 = t110 * t81 - t68 * t77;
t24 = -t74 * t49 + (-t51 * t82 + t78 * t95) * t71;
t121 = t71 * t82;
t106 = t73 * t121;
t72 = cos(pkin(12));
t124 = t71 * t72;
t69 = sin(pkin(12));
t128 = t69 * t71;
t116 = t74 * t82;
t54 = t116 * t72 - t69 * t78;
t56 = -t116 * t69 - t72 * t78;
t132 = g(1) * (t128 * t70 + t56 * t73) - g(2) * (t124 * t70 - t54 * t73) + g(3) * (t70 * t74 + t106);
t131 = g(3) * t71;
t130 = t70 * pkin(9);
t129 = t81 * pkin(3);
t76 = sin(qJ(5));
t127 = t70 * t76;
t80 = cos(qJ(5));
t126 = t70 * t80;
t125 = t70 * t81;
t123 = t71 * t73;
t122 = t71 * t78;
t120 = t73 * t77;
t119 = t73 * t81;
t117 = t74 * t78;
t75 = sin(qJ(6));
t115 = t75 * t80;
t114 = t77 * t78;
t79 = cos(qJ(6));
t113 = t79 * t80;
t52 = pkin(3) * t120 + (-pkin(9) - qJ(4)) * t70;
t55 = t117 * t72 + t69 * t82;
t67 = pkin(2) + t129;
t112 = -t52 * t55 + t54 * t67;
t57 = -t117 * t69 + t72 * t82;
t111 = -t52 * t57 + t56 * t67;
t109 = pkin(3) * t119;
t108 = t70 * t122;
t107 = t71 * t125;
t105 = t121 * t67 - t122 * t52;
t104 = pkin(5) * t80 + pkin(11) * t76;
t102 = t56 * t109 + (t107 * t69 - t57 * t77) * pkin(3);
t50 = t95 * t73;
t27 = -t50 * t55 + t54 * t59;
t28 = t51 * t55 + t54 * t95;
t101 = pkin(4) * t28 - pkin(10) * t27 + t112;
t29 = -t50 * t57 + t56 * t59;
t30 = t51 * t57 + t56 * t95;
t100 = pkin(4) * t30 - pkin(10) * t29 + t111;
t96 = t106 * t129 + (-t114 * t71 + t125 * t74) * pkin(3);
t53 = -t121 * t70 + t73 * t74;
t18 = -t24 * t76 + t53 * t80;
t12 = -t124 * t49 + t51 * t54 - t55 * t95;
t39 = -t123 * t72 - t54 * t70;
t2 = t12 * t76 + t39 * t80;
t17 = -t128 * t49 - t51 * t56 + t57 * t95;
t40 = t123 * t69 - t56 * t70;
t4 = -t17 * t76 + t40 * t80;
t94 = g(1) * t4 + g(2) * t2 + g(3) * t18;
t19 = t24 * t80 + t53 * t76;
t3 = -t12 * t80 + t39 * t76;
t5 = t17 * t80 + t40 * t76;
t93 = g(1) * t5 + g(2) * t3 + g(3) * t19;
t10 = -t126 * t57 + t30 * t76;
t36 = (t51 * t78 + t82 * t95) * t71;
t31 = -t108 * t80 + t36 * t76;
t8 = -t126 * t55 + t28 * t76;
t92 = g(1) * t10 + g(2) * t8 + g(3) * t31;
t91 = -g(1) * t17 + g(2) * t12 - g(3) * t24;
t48 = t95 * t70;
t13 = -t124 * t48 + t50 * t54 + t55 * t59;
t16 = t128 * t48 + t50 * t56 + t57 * t59;
t23 = t74 * t48 + (t50 * t82 + t59 * t78) * t71;
t90 = g(1) * t16 + g(2) * t13 + g(3) * t23;
t35 = t121 * t59 - t122 * t50;
t89 = g(1) * t29 + g(2) * t27 + g(3) * t35;
t88 = pkin(4) * t36 - pkin(10) * t35 + t105;
t87 = g(1) * t57 + g(2) * t55 + g(3) * t122;
t86 = pkin(4) * t16 + pkin(10) * t17 + t102;
t85 = pkin(4) * t23 + pkin(10) * t24 + t96;
t84 = t54 * t109 + (-t107 * t72 - t55 * t77) * pkin(3);
t83 = pkin(4) * t13 - t12 * pkin(10) + t84;
t33 = t87 * t70;
t32 = t108 * t76 + t36 * t80;
t20 = -g(1) * t40 - g(2) * t39 - g(3) * t53;
t11 = t127 * t57 + t30 * t80;
t9 = t127 * t55 + t28 * t80;
t1 = t90 * t76;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t56 - g(2) * t54 - g(3) * t121, t87, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t120 * t57 + t56 * t81) - g(2) * (-t120 * t55 + t54 * t81) - (-t114 * t73 + t81 * t82) * t131, -g(1) * (-t119 * t57 - t56 * t77) - g(2) * (-t119 * t55 - t54 * t77) - (-t119 * t78 - t77 * t82) * t131, -t33, -g(1) * (pkin(2) * t56 + t130 * t57) - g(2) * (pkin(2) * t54 + t130 * t55) - (pkin(2) * t82 + t130 * t78) * t131, 0, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t28 - g(3) * t36, -t89, -t33, -g(1) * t111 - g(2) * t112 - g(3) * t105, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t9 - g(3) * t32, t92, t89, -g(1) * t100 - g(2) * t101 - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t79 - t29 * t75) - g(2) * (-t27 * t75 + t79 * t9) - g(3) * (t32 * t79 - t35 * t75) -g(1) * (-t11 * t75 - t29 * t79) - g(2) * (-t27 * t79 - t75 * t9) - g(3) * (-t32 * t75 - t35 * t79) -t92, -g(1) * (pkin(5) * t11 + pkin(11) * t10 + t100) - g(2) * (pkin(5) * t9 + pkin(11) * t8 + t101) - g(3) * (pkin(5) * t32 + pkin(11) * t31 + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132 * t81 + t77 * t87, t132 * t77 + t81 * t87, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t91, 0, -g(1) * t102 - g(2) * t84 - g(3) * t96, 0, 0, 0, 0, 0, 0, -t90 * t80, t1, t91, -g(1) * t86 - g(2) * t83 - g(3) * t85, 0, 0, 0, 0, 0, 0, -g(1) * (t113 * t16 + t17 * t75) - g(2) * (t113 * t13 - t12 * t75) - g(3) * (t113 * t23 + t24 * t75) -g(1) * (-t115 * t16 + t17 * t79) - g(2) * (-t115 * t13 - t12 * t79) - g(3) * (-t115 * t23 + t24 * t79) -t1, -g(1) * (t104 * t16 + t86) - g(2) * (t104 * t13 + t83) - g(3) * (t104 * t23 + t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t93, 0, 0, 0, 0, 0, 0, 0, 0, -t94 * t79, t94 * t75, -t93, -g(1) * (pkin(5) * t4 + pkin(11) * t5) - g(2) * (pkin(5) * t2 + pkin(11) * t3) - g(3) * (pkin(5) * t18 + pkin(11) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t79 - t5 * t75) - g(2) * (-t13 * t79 - t3 * t75) - g(3) * (-t19 * t75 - t23 * t79) -g(1) * (t16 * t75 - t5 * t79) - g(2) * (t13 * t75 - t3 * t79) - g(3) * (-t19 * t79 + t23 * t75) 0, 0;];
taug_reg  = t6;
