% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t73 = sin(pkin(6));
t81 = cos(qJ(1));
t117 = t73 * t81;
t76 = sin(qJ(2));
t77 = sin(qJ(1));
t80 = cos(qJ(2));
t109 = cos(pkin(6));
t96 = t81 * t109;
t44 = t76 * t96 + t77 * t80;
t72 = qJ(3) + qJ(4);
t68 = qJ(5) + t72;
t63 = sin(t68);
t64 = cos(t68);
t16 = -t63 * t117 + t44 * t64;
t43 = t77 * t76 - t80 * t96;
t74 = sin(qJ(6));
t78 = cos(qJ(6));
t131 = t16 * t74 - t43 * t78;
t130 = t16 * t78 + t43 * t74;
t121 = t73 * t76;
t79 = cos(qJ(3));
t119 = t73 * t79;
t97 = t77 * t109;
t46 = -t76 * t97 + t81 * t80;
t75 = sin(qJ(3));
t23 = t119 * t77 - t46 * t75;
t88 = t117 * t79 + t44 * t75;
t129 = -g(3) * (t109 * t79 - t121 * t75) + g(2) * t88 - g(1) * t23;
t82 = -pkin(10) - pkin(9);
t127 = g(3) * t73;
t66 = sin(t72);
t124 = t46 * t66;
t123 = t64 * t74;
t122 = t64 * t78;
t120 = t73 * t77;
t118 = t73 * t80;
t116 = t74 * t80;
t115 = t78 * t80;
t67 = cos(t72);
t69 = t79 * pkin(3);
t54 = pkin(4) * t67 + t69;
t51 = pkin(2) + t54;
t71 = -pkin(11) + t82;
t114 = -t43 * t51 - t44 * t71;
t45 = t81 * t76 + t80 * t97;
t113 = -t45 * t51 - t46 * t71;
t53 = t75 * pkin(3) + pkin(4) * t66;
t112 = t54 * t120 - t46 * t53;
t111 = t109 * t54 - t53 * t121;
t110 = t81 * pkin(1) + pkin(8) * t120;
t108 = t66 * t121;
t107 = t67 * t120;
t106 = t75 * t120;
t105 = t67 * t117;
t58 = t75 * t117;
t104 = -t77 * pkin(1) + pkin(8) * t117;
t100 = -t64 * t117 - t44 * t63;
t103 = pkin(5) * t100 + t16 * pkin(12);
t19 = -t120 * t64 + t46 * t63;
t20 = t120 * t63 + t46 * t64;
t102 = -t19 * pkin(5) + t20 * pkin(12);
t32 = t109 * t64 - t121 * t63;
t33 = t109 * t63 + t121 * t64;
t101 = t32 * pkin(5) + t33 * pkin(12);
t99 = -t66 * t117 + t44 * t67;
t98 = t44 * t79 - t58;
t95 = t109 * t67;
t94 = -t117 * t54 - t44 * t53;
t93 = t53 * t120 - t45 * t71 + t46 * t51 + t110;
t92 = pkin(5) * t64 + pkin(12) * t63;
t91 = g(1) * t100 + g(2) * t19;
t14 = g(1) * t43 - g(2) * t45;
t90 = g(1) * t81 + g(2) * t77;
t89 = t44 * t66 + t105;
t87 = t53 * t117 + t43 * t71 - t44 * t51 + t104;
t3 = g(1) * t19 - g(2) * t100 - g(3) * t32;
t5 = g(1) * t20 + g(2) * t16 + g(3) * t33;
t85 = -g(1) * t45 - g(2) * t43 + g(3) * t118;
t84 = g(1) * t46 + g(2) * t44 + g(3) * t121;
t65 = t69 + pkin(2);
t56 = pkin(4) * t95;
t52 = pkin(4) * t107;
t38 = t51 * t118;
t24 = t46 * t79 + t106;
t22 = t120 * t66 + t46 * t67;
t21 = t107 - t124;
t10 = t20 * t78 + t45 * t74;
t9 = -t20 * t74 + t45 * t78;
t8 = t85 * t63;
t7 = g(1) * t22 + g(2) * t99 - g(3) * (-t109 * t66 - t121 * t67);
t6 = -g(1) * t21 + g(2) * t89 - g(3) * (t95 - t108);
t2 = t3 * t78;
t1 = t3 * t74;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t77 - g(2) * t81, t90, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t46, -t14, -t90 * t73, -g(1) * t104 - g(2) * t110, 0, 0, 0, 0, 0, 0, g(1) * t98 - g(2) * t24, -g(1) * t88 - g(2) * t23, t14, -g(1) * (-t44 * pkin(2) - t43 * pkin(9) + t104) - g(2) * (t46 * pkin(2) + t45 * pkin(9) + t110) 0, 0, 0, 0, 0, 0, g(1) * t99 - g(2) * t22, -g(1) * t89 - g(2) * t21, t14, -g(1) * (pkin(3) * t58 + t43 * t82 - t44 * t65 + t104) - g(2) * (pkin(3) * t106 - t45 * t82 + t46 * t65 + t110) 0, 0, 0, 0, 0, 0, g(1) * t16 - g(2) * t20, t91, t14, -g(1) * t87 - g(2) * t93, 0, 0, 0, 0, 0, 0, g(1) * t130 - g(2) * t10, -g(1) * t131 - g(2) * t9, -t91, -g(1) * (-pkin(5) * t16 + pkin(12) * t100 + t87) - g(2) * (t20 * pkin(5) + t19 * pkin(12) + t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t84, 0, 0, 0, 0, 0, 0, 0, 0, -t85 * t79, t85 * t75, -t84, -g(1) * (-t45 * pkin(2) + t46 * pkin(9)) - g(2) * (-t43 * pkin(2) + t44 * pkin(9)) - (pkin(2) * t80 + pkin(9) * t76) * t127, 0, 0, 0, 0, 0, 0, -t85 * t67, t85 * t66, -t84, -g(1) * (-t45 * t65 - t46 * t82) - g(2) * (-t43 * t65 - t44 * t82) - (t65 * t80 - t76 * t82) * t127, 0, 0, 0, 0, 0, 0, -t85 * t64, t8, -t84, -g(1) * t113 - g(2) * t114 - g(3) * (-t121 * t71 + t38) 0, 0, 0, 0, 0, 0, -g(1) * (-t122 * t45 + t46 * t74) - g(2) * (-t122 * t43 + t44 * t74) - (t115 * t64 + t74 * t76) * t127, -g(1) * (t123 * t45 + t46 * t78) - g(2) * (t123 * t43 + t44 * t78) - (-t116 * t64 + t76 * t78) * t127, -t8, -g(1) * (-t45 * t92 + t113) - g(2) * (-t43 * t92 + t114) - g(3) * t38 - (-t71 * t76 + t80 * t92) * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, g(1) * t24 + g(2) * t98 - g(3) * (-t109 * t75 - t119 * t76) 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t129 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t112 - g(2) * t94 - g(3) * t111, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (t102 + t112) - g(2) * (t103 + t94) - g(3) * (t101 + t111); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t52 - g(3) * t56 + (g(2) * t105 + t66 * t84) * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (-pkin(4) * t124 + t102 + t52) - g(2) * (-pkin(4) * t89 + t103) - g(3) * (-pkin(4) * t108 + t101 + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t102 - g(2) * t103 - g(3) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t131 - g(3) * (-t115 * t73 - t33 * t74) g(1) * t10 + g(2) * t130 - g(3) * (t116 * t73 - t33 * t78) 0, 0;];
taug_reg  = t4;
