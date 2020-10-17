% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP6
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:08:01
% EndTime: 2019-05-07 08:08:04
% DurationCPUTime: 0.85s
% Computational Cost: add. (739->151), mult. (1355->226), div. (0->0), fcn. (1627->12), ass. (0->85)
t58 = sin(pkin(6));
t63 = sin(qJ(2));
t101 = t58 * t63;
t62 = sin(qJ(3));
t66 = cos(qJ(3));
t91 = cos(pkin(6));
t116 = -t101 * t62 + t66 * t91;
t110 = cos(qJ(1));
t67 = cos(qJ(2));
t64 = sin(qJ(1));
t84 = t64 * t91;
t38 = t110 * t67 - t63 * t84;
t99 = t58 * t66;
t23 = -t38 * t62 + t64 * t99;
t78 = t91 * t110;
t36 = t63 * t78 + t64 * t67;
t57 = qJ(3) + pkin(11);
t54 = sin(t57);
t55 = cos(t57);
t87 = t58 * t110;
t18 = t36 * t55 - t54 * t87;
t35 = t63 * t64 - t67 * t78;
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t115 = t18 * t61 - t35 * t65;
t109 = t35 * t61;
t114 = t18 * t65 + t109;
t100 = t58 * t64;
t22 = t100 * t54 + t38 * t55;
t37 = t110 * t63 + t67 * t84;
t12 = -t22 * t61 + t37 * t65;
t30 = t101 * t55 + t54 * t91;
t95 = t65 * t67;
t1 = g(2) * t115 - g(3) * (-t30 * t61 - t58 * t95) - g(1) * t12;
t53 = pkin(3) * t66 + pkin(2);
t98 = t58 * t67;
t39 = t53 * t98;
t112 = g(3) * t39;
t111 = g(3) * t58;
t107 = t36 * t61;
t106 = t37 * t61;
t105 = t38 * t61;
t103 = t55 * t61;
t102 = t55 * t65;
t60 = -qJ(4) - pkin(9);
t97 = t60 * t63;
t96 = t61 * t67;
t94 = -t35 * t53 - t36 * t60;
t93 = -t37 * t53 - t38 * t60;
t92 = pkin(1) * t110 + pkin(8) * t100;
t89 = t62 * t100;
t86 = -pkin(1) * t64 + pkin(8) * t87;
t46 = t62 * t87;
t85 = t36 * t66 - t46;
t82 = t23 * pkin(3);
t81 = pkin(3) * t89 - t37 * t60 + t38 * t53 + t92;
t80 = pkin(4) * t55 + pkin(10) * t54;
t17 = t36 * t54 + t55 * t87;
t21 = -t100 * t55 + t38 * t54;
t79 = -g(1) * t17 + g(2) * t21;
t16 = g(1) * t35 - g(2) * t37;
t52 = pkin(5) * t65 + pkin(4);
t59 = -qJ(6) - pkin(10);
t77 = t52 * t55 - t54 * t59;
t76 = t116 * pkin(3);
t75 = g(1) * t110 + g(2) * t64;
t74 = pkin(3) * t46 + t35 * t60 - t36 * t53 + t86;
t29 = t101 * t54 - t55 * t91;
t73 = g(1) * t21 + g(2) * t17 + g(3) * t29;
t72 = g(1) * t22 + g(2) * t18 + g(3) * t30;
t71 = t36 * t62 + t66 * t87;
t14 = -g(1) * t37 - g(2) * t35 + g(3) * t98;
t70 = g(1) * t38 + g(2) * t36 + g(3) * t101;
t69 = t71 * pkin(3);
t24 = t38 * t66 + t89;
t13 = t22 * t65 + t106;
t11 = t14 * t54;
t8 = t73 * t65;
t7 = t73 * t61;
t6 = -g(1) * (-t102 * t37 + t105) - g(2) * (-t102 * t35 + t107) - (t55 * t95 + t61 * t63) * t111;
t5 = -g(1) * (t103 * t37 + t38 * t65) - g(2) * (t103 * t35 + t36 * t65) - (-t55 * t96 + t63 * t65) * t111;
t4 = g(1) * t114 - g(2) * t13;
t3 = -g(1) * t115 - g(2) * t12;
t2 = g(1) * t13 + g(2) * t114 - g(3) * (-t30 * t65 + t58 * t96);
t9 = [0, 0, 0, 0, 0, 0, g(1) * t64 - g(2) * t110, t75, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t36 - g(2) * t38, -t16, -t75 * t58, -g(1) * t86 - g(2) * t92, 0, 0, 0, 0, 0, 0, g(1) * t85 - g(2) * t24, -g(1) * t71 - g(2) * t23, t16, -g(1) * (-pkin(2) * t36 - pkin(9) * t35 + t86) - g(2) * (pkin(2) * t38 + pkin(9) * t37 + t92) 0, 0, 0, 0, 0, 0, g(1) * t18 - g(2) * t22, t79, t16, -g(1) * t74 - g(2) * t81, 0, 0, 0, 0, 0, 0, t4, t3, -t79, -g(1) * (-pkin(4) * t18 - pkin(10) * t17 + t74) - g(2) * (pkin(4) * t22 + pkin(10) * t21 + t81) 0, 0, 0, 0, 0, 0, t4, t3, -t79, -g(1) * (-pkin(5) * t109 + t17 * t59 - t18 * t52 + t74) - g(2) * (pkin(5) * t106 - t21 * t59 + t22 * t52 + t81); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t70, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t66, t14 * t62, -t70, -g(1) * (-pkin(2) * t37 + pkin(9) * t38) - g(2) * (-pkin(2) * t35 + pkin(9) * t36) - (pkin(2) * t67 + pkin(9) * t63) * t111, 0, 0, 0, 0, 0, 0, -t14 * t55, t11, -t70, -g(1) * t93 - g(2) * t94 - g(3) * (-t58 * t97 + t39) 0, 0, 0, 0, 0, 0, t6, t5, -t11, -g(1) * (-t37 * t80 + t93) - g(2) * (-t35 * t80 + t94) - t112 - (t67 * t80 - t97) * t111, 0, 0, 0, 0, 0, 0, t6, t5, -t11, -g(1) * (pkin(5) * t105 - t37 * t77 + t93) - g(2) * (pkin(5) * t107 - t35 * t77 + t94) - t112 - (t77 * t67 + (pkin(5) * t61 - t60) * t63) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 + g(2) * t71 - g(3) * t116, g(1) * t24 + g(2) * t85 - g(3) * (-t62 * t91 - t63 * t99) 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, -g(1) * t82 + g(2) * t69 - g(3) * t76, 0, 0, 0, 0, 0, 0, t8, -t7, -t72, -g(1) * (-pkin(4) * t21 + pkin(10) * t22 + t82) - g(2) * (-t17 * pkin(4) + t18 * pkin(10) - t69) - g(3) * (-pkin(4) * t29 + pkin(10) * t30 + t76) 0, 0, 0, 0, 0, 0, t8, -t7, -t72, -g(1) * (-t21 * t52 - t22 * t59 + t82) - g(2) * (-t17 * t52 - t18 * t59 - t69) - g(3) * (-t29 * t52 - t30 * t59 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73;];
taug_reg  = t9;
