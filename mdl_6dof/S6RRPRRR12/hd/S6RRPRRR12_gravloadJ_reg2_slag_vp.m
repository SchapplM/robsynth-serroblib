% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t64 = sin(pkin(6));
t72 = cos(qJ(1));
t109 = t64 * t72;
t67 = sin(qJ(2));
t68 = sin(qJ(1));
t71 = cos(qJ(2));
t102 = cos(pkin(6));
t92 = t72 * t102;
t42 = t67 * t68 - t71 * t92;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t80 = t109 * t66 + t42 * t70;
t111 = t64 * t68;
t93 = t68 * t102;
t44 = t67 * t72 + t71 * t93;
t20 = -t111 * t66 + t44 * t70;
t124 = pkin(4) * t66;
t43 = t67 * t92 + t68 * t71;
t73 = -pkin(10) - pkin(9);
t129 = t124 * t43 + t42 * t73;
t45 = -t67 * t93 + t71 * t72;
t128 = t124 * t45 + t44 * t73;
t127 = -g(1) * t45 - g(2) * t43;
t65 = sin(qJ(6));
t69 = cos(qJ(6));
t63 = qJ(4) + qJ(5);
t60 = sin(t63);
t61 = cos(t63);
t82 = t109 * t61 - t42 * t60;
t126 = t43 * t69 + t65 * t82;
t125 = -t43 * t65 + t69 * t82;
t121 = g(3) * t64;
t120 = t42 * t66;
t116 = t44 * t66;
t114 = t60 * t65;
t113 = t60 * t69;
t112 = t64 * t67;
t110 = t64 * t71;
t108 = t65 * t67;
t107 = t67 * t69;
t106 = t71 * t73;
t105 = t80 * pkin(4);
t104 = pkin(2) * t110 + qJ(3) * t112;
t103 = pkin(1) * t72 + pkin(8) * t111;
t99 = t112 * t124 + t104;
t98 = -pkin(1) * t68 + pkin(8) * t109;
t14 = t111 * t60 - t44 * t61;
t15 = t111 * t61 + t44 * t60;
t97 = -pkin(5) * t14 + pkin(11) * t15;
t81 = t109 * t60 + t42 * t61;
t96 = pkin(5) * t81 - pkin(11) * t82;
t26 = -t102 * t60 - t110 * t61;
t27 = t102 * t61 - t110 * t60;
t95 = pkin(5) * t26 + t27 * pkin(11);
t94 = -t109 * t70 + t120;
t38 = t42 * pkin(2);
t91 = qJ(3) * t43 - t38;
t40 = t44 * pkin(2);
t90 = qJ(3) * t45 - t40;
t89 = pkin(5) * t60 - pkin(11) * t61;
t88 = g(1) * t81 + g(2) * t14;
t87 = g(1) * t42 - g(2) * t44;
t11 = g(1) * t43 - g(2) * t45;
t86 = g(1) * t72 + g(2) * t68;
t85 = t20 * pkin(4);
t84 = pkin(2) * t45 + qJ(3) * t44 + t103;
t79 = -pkin(2) * t43 - qJ(3) * t42 + t98;
t3 = g(1) * t14 - g(2) * t81 - g(3) * t26;
t5 = g(1) * t15 - g(2) * t82 + g(3) * t27;
t78 = -t102 * t66 - t110 * t70;
t9 = -g(1) * t44 - g(2) * t42 + g(3) * t110;
t77 = g(3) * t112 - t127;
t59 = pkin(4) * t70 + pkin(3);
t76 = pkin(4) * t116 + t111 * t59 - t45 * t73 + t84;
t75 = t78 * pkin(4);
t74 = -pkin(4) * t120 + t109 * t59 + t43 * t73 + t79;
t46 = t86 * t64;
t21 = t111 * t70 + t116;
t8 = t15 * t69 + t45 * t65;
t7 = -t15 * t65 + t45 * t69;
t6 = t77 * t61;
t2 = t3 * t69;
t1 = t3 * t65;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t72, t86, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t87, -t46, -g(1) * t98 - g(2) * t103, 0, 0, 0, 0, 0, 0, -t46, -t11, t87, -g(1) * t79 - g(2) * t84, 0, 0, 0, 0, 0, 0, g(1) * t94 - g(2) * t21, g(1) * t80 - g(2) * t20, t11, -g(1) * (pkin(3) * t109 - pkin(9) * t43 + t79) - g(2) * (pkin(3) * t111 + pkin(9) * t45 + t84) 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t15, t88, t11, -g(1) * t74 - g(2) * t76, 0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t8, g(1) * t126 - g(2) * t7, -t88, -g(1) * (pkin(5) * t82 + pkin(11) * t81 + t74) - g(2) * (pkin(5) * t15 + pkin(11) * t14 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t77, -g(1) * t90 - g(2) * t91 - g(3) * t104, 0, 0, 0, 0, 0, 0, -t77 * t66, -t77 * t70, -t9, -g(1) * (-pkin(9) * t44 + t90) - g(2) * (-pkin(9) * t42 + t91) - g(3) * (pkin(9) * t110 + t104) 0, 0, 0, 0, 0, 0, -t77 * t60, -t6, -t9, -g(1) * (t90 + t128) - g(2) * (t91 + t129) - g(3) * (-t106 * t64 + t99) 0, 0, 0, 0, 0, 0, -g(1) * (t113 * t45 - t44 * t65) - g(2) * (t113 * t43 - t42 * t65) - (t107 * t60 + t65 * t71) * t121, -g(1) * (-t114 * t45 - t44 * t69) - g(2) * (-t114 * t43 - t42 * t69) - (-t108 * t60 + t69 * t71) * t121, t6, -g(1) * (-t40 + t128) - g(2) * (-t38 + t129) - g(3) * t99 - (t67 * t89 - t106) * t121 + t127 * (qJ(3) + t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t80 - g(3) * t78, g(1) * t21 + g(2) * t94 - g(3) * (-t102 * t70 + t110 * t66) 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t85 - g(2) * t105 - g(3) * t75, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (t85 + t97) - g(2) * (t96 + t105) - g(3) * (t75 + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t97 - g(2) * t96 - g(3) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t126 - g(3) * (t107 * t64 - t27 * t65) g(1) * t8 - g(2) * t125 - g(3) * (-t108 * t64 - t27 * t69) 0, 0;];
taug_reg  = t4;
