% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t68 = sin(qJ(4));
t109 = qJ(5) * t68;
t65 = sin(pkin(11));
t70 = sin(qJ(2));
t74 = cos(qJ(2));
t107 = cos(pkin(11));
t108 = cos(pkin(6));
t83 = t108 * t107;
t51 = t65 * t74 + t70 * t83;
t69 = sin(qJ(3));
t73 = cos(qJ(3));
t66 = sin(pkin(6));
t96 = t66 * t107;
t28 = -t51 * t69 - t73 * t96;
t72 = cos(qJ(4));
t121 = t28 * t72;
t131 = pkin(4) * t121 + t28 * t109;
t115 = t66 * t73;
t97 = t65 * t108;
t53 = t107 * t74 - t70 * t97;
t30 = t65 * t115 - t53 * t69;
t120 = t30 * t72;
t130 = pkin(4) * t120 + t30 * t109;
t116 = t66 * t70;
t54 = t108 * t73 - t69 * t116;
t117 = t54 * t72;
t129 = pkin(4) * t117 + t54 * t109;
t80 = g(1) * t30 + g(2) * t28 + g(3) * t54;
t29 = t51 * t73 - t69 * t96;
t50 = t65 * t70 - t74 * t83;
t10 = t29 * t68 - t50 * t72;
t11 = t29 * t72 + t50 * t68;
t31 = t65 * t66 * t69 + t53 * t73;
t52 = t107 * t70 + t74 * t97;
t12 = t31 * t68 - t52 * t72;
t13 = t31 * t72 + t52 * t68;
t114 = t66 * t74;
t55 = t108 * t69 + t70 * t115;
t32 = t72 * t114 + t55 * t68;
t105 = t68 * t114;
t33 = t55 * t72 - t105;
t67 = sin(qJ(6));
t71 = cos(qJ(6));
t128 = g(1) * (t12 * t67 + t13 * t71) + g(2) * (t10 * t67 + t11 * t71) + g(3) * (t32 * t67 + t33 * t71);
t127 = g(1) * (t12 * t71 - t13 * t67) + g(2) * (t10 * t71 - t11 * t67) + g(3) * (t32 * t71 - t33 * t67);
t126 = pkin(9) - pkin(10);
t125 = pkin(3) * t73;
t119 = t50 * t69;
t118 = t52 * t69;
t113 = t68 * t73;
t112 = t72 * t73;
t111 = t73 * t74;
t110 = pkin(2) * t114 + pkin(8) * t116;
t106 = t69 * t114;
t104 = -t50 * pkin(2) + t51 * pkin(8);
t103 = -t52 * pkin(2) + t53 * pkin(8);
t25 = t28 * pkin(3);
t102 = t29 * pkin(9) + t25;
t26 = t30 * pkin(3);
t101 = t31 * pkin(9) + t26;
t49 = t54 * pkin(3);
t100 = t55 * pkin(9) + t49;
t99 = -t10 * pkin(4) + t11 * qJ(5);
t98 = -t12 * pkin(4) + t13 * qJ(5);
t95 = -t32 * pkin(4) + t33 * qJ(5);
t94 = t66 * pkin(3) * t111 + pkin(9) * t106 + t110;
t85 = -pkin(9) * t119 - t50 * t125 + t104;
t84 = -pkin(9) * t118 - t52 * t125 + t103;
t2 = g(1) * t12 + g(2) * t10 + g(3) * t32;
t82 = g(1) * t13 + g(2) * t11 + g(3) * t33;
t17 = -t50 * t113 - t51 * t72;
t19 = -t52 * t113 - t53 * t72;
t35 = t73 * t105 - t72 * t116;
t81 = g(1) * t19 + g(2) * t17 + g(3) * t35;
t7 = g(1) * t31 + g(2) * t29 + g(3) * t55;
t36 = (t72 * t111 + t68 * t70) * t66;
t79 = t36 * pkin(4) + t35 * qJ(5) + t94;
t78 = -g(1) * t52 - g(2) * t50 + g(3) * t114;
t77 = g(1) * t53 + g(2) * t51 + g(3) * t116;
t18 = -t50 * t112 + t51 * t68;
t76 = t18 * pkin(4) + t17 * qJ(5) + t85;
t20 = -t52 * t112 + t53 * t68;
t75 = t20 * pkin(4) + t19 * qJ(5) + t84;
t14 = t78 * t69;
t5 = t80 * t72;
t4 = t80 * t68;
t3 = -g(1) * t20 - g(2) * t18 - g(3) * t36;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t77, 0, 0, 0, 0, 0, 0, 0, 0, -t78 * t73, t14, -t77, -g(1) * t103 - g(2) * t104 - g(3) * t110, 0, 0, 0, 0, 0, 0, t3, t81, -t14, -g(1) * t84 - g(2) * t85 - g(3) * t94, 0, 0, 0, 0, 0, 0, t3, -t14, -t81, -g(1) * t75 - g(2) * t76 - g(3) * t79, 0, 0, 0, 0, 0, 0, -g(1) * (t19 * t67 + t20 * t71) - g(2) * (t17 * t67 + t18 * t71) - g(3) * (t35 * t67 + t36 * t71) -g(1) * (t19 * t71 - t20 * t67) - g(2) * (t17 * t71 - t18 * t67) - g(3) * (t35 * t71 - t36 * t67) t14, -g(1) * (t20 * pkin(5) + pkin(10) * t118 + t75) - g(2) * (t18 * pkin(5) + pkin(10) * t119 + t76) - g(3) * (t36 * pkin(5) - pkin(10) * t106 + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, -t7, -g(1) * t101 - g(2) * t102 - g(3) * t100, 0, 0, 0, 0, 0, 0, -t5, -t7, -t4, -g(1) * (t101 + t130) - g(2) * (t102 + t131) - g(3) * (t100 + t129) 0, 0, 0, 0, 0, 0, -t80 * (t67 * t68 + t71 * t72) t80 * (t67 * t72 - t68 * t71) t7, -g(1) * (pkin(5) * t120 + t126 * t31 + t130 + t26) - g(2) * (pkin(5) * t121 + t126 * t29 + t131 + t25) - g(3) * (pkin(5) * t117 + t126 * t55 + t129 + t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t82, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t82, -g(1) * t98 - g(2) * t99 - g(3) * t95, 0, 0, 0, 0, 0, 0, t127, -t128, 0, -g(1) * (-t12 * pkin(5) + t98) - g(2) * (-t10 * pkin(5) + t99) - g(3) * (-t32 * pkin(5) + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, t128, 0, 0;];
taug_reg  = t1;
