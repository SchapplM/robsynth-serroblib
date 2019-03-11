% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t74 = sin(pkin(6));
t82 = cos(qJ(1));
t117 = t74 * t82;
t77 = sin(qJ(2));
t78 = sin(qJ(1));
t81 = cos(qJ(2));
t109 = cos(pkin(6));
t96 = t82 * t109;
t45 = t77 * t96 + t78 * t81;
t73 = qJ(3) + qJ(4);
t67 = pkin(12) + t73;
t63 = sin(t67);
t64 = cos(t67);
t16 = -t63 * t117 + t45 * t64;
t44 = t78 * t77 - t81 * t96;
t75 = sin(qJ(6));
t79 = cos(qJ(6));
t131 = t16 * t75 - t44 * t79;
t130 = t16 * t79 + t44 * t75;
t121 = t74 * t77;
t80 = cos(qJ(3));
t119 = t74 * t80;
t97 = t78 * t109;
t47 = -t77 * t97 + t82 * t81;
t76 = sin(qJ(3));
t24 = t78 * t119 - t47 * t76;
t88 = t80 * t117 + t45 * t76;
t129 = -g(3) * (t109 * t80 - t76 * t121) + g(2) * t88 - g(1) * t24;
t83 = -pkin(10) - pkin(9);
t127 = g(3) * t74;
t68 = sin(t73);
t124 = t47 * t68;
t123 = t64 * t75;
t122 = t64 * t79;
t120 = t74 * t78;
t118 = t74 * t81;
t116 = t75 * t81;
t115 = t79 * t81;
t69 = cos(t73);
t70 = t80 * pkin(3);
t55 = pkin(4) * t69 + t70;
t52 = pkin(2) + t55;
t72 = -qJ(5) + t83;
t114 = -t44 * t52 - t45 * t72;
t46 = t82 * t77 + t81 * t97;
t113 = -t46 * t52 - t47 * t72;
t54 = t76 * pkin(3) + pkin(4) * t68;
t112 = t55 * t120 - t47 * t54;
t111 = t109 * t55 - t54 * t121;
t110 = t82 * pkin(1) + pkin(8) * t120;
t108 = t68 * t121;
t107 = t69 * t120;
t106 = t76 * t120;
t105 = t69 * t117;
t59 = t76 * t117;
t104 = -t78 * pkin(1) + pkin(8) * t117;
t100 = -t64 * t117 - t45 * t63;
t103 = pkin(5) * t100 + t16 * pkin(11);
t19 = -t64 * t120 + t47 * t63;
t20 = t63 * t120 + t47 * t64;
t102 = -t19 * pkin(5) + t20 * pkin(11);
t33 = t109 * t64 - t63 * t121;
t34 = t109 * t63 + t64 * t121;
t101 = t33 * pkin(5) + t34 * pkin(11);
t99 = -t68 * t117 + t45 * t69;
t98 = t45 * t80 - t59;
t95 = t109 * t69;
t94 = -t55 * t117 - t45 * t54;
t93 = t54 * t120 - t46 * t72 + t47 * t52 + t110;
t92 = pkin(5) * t64 + pkin(11) * t63;
t91 = g(1) * t100 + g(2) * t19;
t21 = g(1) * t44 - g(2) * t46;
t90 = g(1) * t82 + g(2) * t78;
t89 = t45 * t68 + t105;
t87 = t54 * t117 + t44 * t72 - t45 * t52 + t104;
t3 = g(1) * t19 - g(2) * t100 - g(3) * t33;
t5 = g(1) * t20 + g(2) * t16 + g(3) * t34;
t11 = -g(1) * t46 - g(2) * t44 + g(3) * t118;
t85 = g(1) * t47 + g(2) * t45 + g(3) * t121;
t66 = t70 + pkin(2);
t57 = pkin(4) * t95;
t53 = pkin(4) * t107;
t39 = t52 * t118;
t25 = t47 * t80 + t106;
t23 = t68 * t120 + t47 * t69;
t22 = t107 - t124;
t10 = t20 * t79 + t46 * t75;
t9 = -t20 * t75 + t46 * t79;
t8 = t11 * t63;
t7 = g(1) * t23 + g(2) * t99 - g(3) * (-t109 * t68 - t69 * t121);
t6 = -g(1) * t22 + g(2) * t89 - g(3) * (t95 - t108);
t2 = t3 * t79;
t1 = t3 * t75;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t78 - g(2) * t82, t90, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t45 - g(2) * t47, -t21, -t90 * t74, -g(1) * t104 - g(2) * t110, 0, 0, 0, 0, 0, 0, g(1) * t98 - g(2) * t25, -g(1) * t88 - g(2) * t24, t21, -g(1) * (-t45 * pkin(2) - t44 * pkin(9) + t104) - g(2) * (t47 * pkin(2) + t46 * pkin(9) + t110) 0, 0, 0, 0, 0, 0, g(1) * t99 - g(2) * t23, -g(1) * t89 - g(2) * t22, t21, -g(1) * (pkin(3) * t59 + t44 * t83 - t45 * t66 + t104) - g(2) * (pkin(3) * t106 - t46 * t83 + t47 * t66 + t110) 0, 0, 0, 0, 0, 0, g(1) * t16 - g(2) * t20, t91, t21, -g(1) * t87 - g(2) * t93, 0, 0, 0, 0, 0, 0, g(1) * t130 - g(2) * t10, -g(1) * t131 - g(2) * t9, -t91, -g(1) * (-pkin(5) * t16 + pkin(11) * t100 + t87) - g(2) * (t20 * pkin(5) + t19 * pkin(11) + t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t85, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t80, t11 * t76, -t85, -g(1) * (-t46 * pkin(2) + t47 * pkin(9)) - g(2) * (-t44 * pkin(2) + t45 * pkin(9)) - (pkin(2) * t81 + pkin(9) * t77) * t127, 0, 0, 0, 0, 0, 0, -t11 * t69, t11 * t68, -t85, -g(1) * (-t46 * t66 - t47 * t83) - g(2) * (-t44 * t66 - t45 * t83) - (t66 * t81 - t77 * t83) * t127, 0, 0, 0, 0, 0, 0, -t11 * t64, t8, -t85, -g(1) * t113 - g(2) * t114 - g(3) * (-t72 * t121 + t39) 0, 0, 0, 0, 0, 0, -g(1) * (-t46 * t122 + t47 * t75) - g(2) * (-t44 * t122 + t45 * t75) - (t64 * t115 + t75 * t77) * t127, -g(1) * (t46 * t123 + t47 * t79) - g(2) * (t44 * t123 + t45 * t79) - (-t64 * t116 + t77 * t79) * t127, -t8, -g(1) * (-t92 * t46 + t113) - g(2) * (-t92 * t44 + t114) - g(3) * t39 - (-t72 * t77 + t92 * t81) * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, g(1) * t25 + g(2) * t98 - g(3) * (-t109 * t76 - t77 * t119) 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t129 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t112 - g(2) * t94 - g(3) * t111, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (t102 + t112) - g(2) * (t103 + t94) - g(3) * (t101 + t111); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t53 - g(3) * t57 + (g(2) * t105 + t85 * t68) * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (-pkin(4) * t124 + t102 + t53) - g(2) * (-t89 * pkin(4) + t103) - g(3) * (-pkin(4) * t108 + t101 + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t131 - g(3) * (-t74 * t115 - t34 * t75) g(1) * t10 + g(2) * t130 - g(3) * (t74 * t116 - t34 * t79) 0, 0;];
taug_reg  = t4;
