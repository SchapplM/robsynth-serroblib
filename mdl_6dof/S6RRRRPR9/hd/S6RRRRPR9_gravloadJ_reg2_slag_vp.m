% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR9
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t104 = cos(pkin(6));
t69 = sin(pkin(6));
t73 = sin(qJ(2));
t114 = t69 * t73;
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t132 = t104 * t75 - t114 * t72;
t112 = t69 * t75;
t126 = cos(qJ(2));
t127 = cos(qJ(1));
t74 = sin(qJ(1));
t94 = t73 * t104;
t45 = t126 * t127 - t74 * t94;
t27 = t112 * t74 - t45 * t72;
t43 = t126 * t74 + t127 * t94;
t67 = qJ(3) + qJ(4);
t63 = sin(t67);
t64 = cos(t67);
t98 = t69 * t127;
t22 = t43 * t64 - t63 * t98;
t86 = t104 * t126;
t42 = -t127 * t86 + t73 * t74;
t66 = pkin(12) + qJ(6);
t61 = sin(t66);
t62 = cos(t66);
t131 = t22 * t61 - t42 * t62;
t130 = t22 * t62 + t42 * t61;
t60 = pkin(3) * t75 + pkin(2);
t97 = t69 * t126;
t46 = t60 * t97;
t129 = g(3) * t46;
t128 = g(3) * t69;
t68 = sin(pkin(12));
t123 = t42 * t68;
t122 = t43 * t68;
t44 = t127 * t73 + t74 * t86;
t121 = t44 * t68;
t120 = t45 * t68;
t118 = t61 * t64;
t117 = t62 * t64;
t116 = t64 * t68;
t70 = cos(pkin(12));
t115 = t64 * t70;
t113 = t69 * t74;
t76 = -pkin(10) - pkin(9);
t111 = t73 * t76;
t21 = t43 * t63 + t64 * t98;
t59 = pkin(5) * t70 + pkin(4);
t71 = -pkin(11) - qJ(5);
t110 = -t21 * t59 - t22 * t71;
t25 = -t113 * t64 + t45 * t63;
t26 = t113 * t63 + t45 * t64;
t109 = -t25 * t59 - t26 * t71;
t36 = -t104 * t64 + t114 * t63;
t37 = t104 * t63 + t114 * t64;
t108 = -t36 * t59 - t37 * t71;
t107 = -t42 * t60 - t43 * t76;
t106 = -t44 * t60 - t45 * t76;
t105 = pkin(1) * t127 + pkin(8) * t113;
t102 = t72 * t113;
t100 = t63 * t126;
t99 = t64 * t126;
t96 = -pkin(1) * t74 + pkin(8) * t98;
t53 = t72 * t98;
t95 = t43 * t75 - t53;
t92 = -pkin(4) * t21 + t22 * qJ(5);
t91 = -pkin(4) * t25 + qJ(5) * t26;
t90 = -pkin(4) * t36 + qJ(5) * t37;
t89 = t27 * pkin(3);
t88 = pkin(3) * t102 - t44 * t76 + t45 * t60 + t105;
t87 = -g(1) * t21 + g(2) * t25;
t18 = g(1) * t42 - g(2) * t44;
t85 = -pkin(4) * t64 - qJ(5) * t63;
t84 = -t59 * t64 + t63 * t71;
t83 = t132 * pkin(3);
t82 = g(1) * t127 + g(2) * t74;
t81 = pkin(3) * t53 + t42 * t76 - t43 * t60 + t96;
t6 = g(1) * t25 + g(2) * t21 + g(3) * t36;
t8 = g(1) * t26 + g(2) * t22 + g(3) * t37;
t80 = t43 * t72 + t75 * t98;
t79 = g(1) * t45 + g(2) * t43 + g(3) * t114;
t78 = t80 * pkin(3);
t77 = -g(1) * t44 - g(2) * t42 + g(3) * t97;
t28 = t45 * t75 + t102;
t12 = t77 * t63;
t11 = t26 * t62 + t44 * t61;
t10 = -t26 * t61 + t44 * t62;
t4 = t6 * t70;
t3 = t6 * t68;
t2 = t6 * t62;
t1 = t6 * t61;
t5 = [0, 0, 0, 0, 0, 0, g(1) * t74 - g(2) * t127, t82, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t43 - g(2) * t45, -t18, -t82 * t69, -g(1) * t96 - g(2) * t105, 0, 0, 0, 0, 0, 0, g(1) * t95 - g(2) * t28, -g(1) * t80 - g(2) * t27, t18, -g(1) * (-pkin(2) * t43 - pkin(9) * t42 + t96) - g(2) * (pkin(2) * t45 + pkin(9) * t44 + t105) 0, 0, 0, 0, 0, 0, g(1) * t22 - g(2) * t26, t87, t18, -g(1) * t81 - g(2) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (-t22 * t70 - t123) - g(2) * (t26 * t70 + t121) -g(1) * (t22 * t68 - t42 * t70) - g(2) * (-t26 * t68 + t44 * t70) -t87, -g(1) * (-pkin(4) * t22 - qJ(5) * t21 + t81) - g(2) * (pkin(4) * t26 + qJ(5) * t25 + t88) 0, 0, 0, 0, 0, 0, g(1) * t130 - g(2) * t11, -g(1) * t131 - g(2) * t10, -t87, -g(1) * (-pkin(5) * t123 + t21 * t71 - t22 * t59 + t81) - g(2) * (pkin(5) * t121 - t25 * t71 + t26 * t59 + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t79, 0, 0, 0, 0, 0, 0, 0, 0, -t77 * t75, t77 * t72, -t79, -g(1) * (-pkin(2) * t44 + pkin(9) * t45) - g(2) * (-pkin(2) * t42 + pkin(9) * t43) - (pkin(2) * t126 + pkin(9) * t73) * t128, 0, 0, 0, 0, 0, 0, -t77 * t64, t12, -t79, -g(1) * t106 - g(2) * t107 - g(3) * (-t111 * t69 + t46) 0, 0, 0, 0, 0, 0, -g(1) * (-t115 * t44 + t120) - g(2) * (-t115 * t42 + t122) - (t68 * t73 + t70 * t99) * t128, -g(1) * (t116 * t44 + t45 * t70) - g(2) * (t116 * t42 + t43 * t70) - (-t68 * t99 + t70 * t73) * t128, -t12, -g(1) * (t44 * t85 + t106) - g(2) * (t42 * t85 + t107) - t129 - (pkin(4) * t99 + qJ(5) * t100 - t111) * t128, 0, 0, 0, 0, 0, 0, -g(1) * (-t117 * t44 + t45 * t61) - g(2) * (-t117 * t42 + t43 * t61) - (t61 * t73 + t62 * t99) * t128, -g(1) * (t118 * t44 + t45 * t62) - g(2) * (t118 * t42 + t43 * t62) - (-t61 * t99 + t62 * t73) * t128, -t12, -g(1) * (pkin(5) * t120 + t44 * t84 + t106) - g(2) * (pkin(5) * t122 + t42 * t84 + t107) - t129 - (t59 * t99 - t71 * t100 + (pkin(5) * t68 - t76) * t73) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t27 + g(2) * t80 - g(3) * t132, g(1) * t28 + g(2) * t95 - g(3) * (-t104 * t72 - t112 * t73) 0, 0, 0, 0, 0, 0, 0, 0, t6, t8, 0, -g(1) * t89 + g(2) * t78 - g(3) * t83, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t89 + t91) - g(2) * (-t78 + t92) - g(3) * (t83 + t90) 0, 0, 0, 0, 0, 0, t2, -t1, -t8, -g(1) * (t89 + t109) - g(2) * (-t78 + t110) - g(3) * (t83 + t108); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t8, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t91 - g(2) * t92 - g(3) * t90, 0, 0, 0, 0, 0, 0, t2, -t1, -t8, -g(1) * t109 - g(2) * t110 - g(3) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t131 - g(3) * (-t37 * t61 - t62 * t97) g(1) * t11 + g(2) * t130 - g(3) * (-t37 * t62 + t61 * t97) 0, 0;];
taug_reg  = t5;
