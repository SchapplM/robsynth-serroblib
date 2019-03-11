% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t73 = sin(pkin(6));
t82 = cos(qJ(1));
t120 = t73 * t82;
t113 = cos(pkin(12));
t72 = sin(pkin(12));
t77 = sin(qJ(2));
t81 = cos(qJ(2));
t53 = -t77 * t113 - t81 * t72;
t74 = cos(pkin(6));
t45 = t53 * t74;
t78 = sin(qJ(1));
t89 = t81 * t113 - t77 * t72;
t31 = -t82 * t45 + t78 * t89;
t71 = qJ(4) + qJ(5);
t69 = sin(t71);
t70 = cos(t71);
t14 = -t69 * t120 + t31 * t70;
t85 = t89 * t74;
t30 = t78 * t53 + t82 * t85;
t75 = sin(qJ(6));
t79 = cos(qJ(6));
t133 = t14 * t75 + t30 * t79;
t132 = t14 * t79 - t30 * t75;
t32 = -t78 * t45 - t82 * t89;
t44 = t53 * t73;
t7 = g(1) * t32 - g(2) * t31 + g(3) * t44;
t121 = t73 * t81;
t115 = t82 * t77;
t116 = t78 * t81;
t49 = -t74 * t116 - t115;
t131 = -g(1) * t49 - g(3) * t121;
t76 = sin(qJ(4));
t126 = t32 * t76;
t125 = t44 * t76;
t124 = t70 * t75;
t123 = t70 * t79;
t122 = t73 * t78;
t80 = cos(qJ(4));
t119 = t74 * t80;
t117 = t78 * t77;
t114 = t82 * t81;
t111 = t76 * t122;
t110 = t80 * t122;
t63 = t76 * t120;
t109 = t80 * t120;
t108 = t74 * t114;
t43 = t89 * t73;
t64 = pkin(2) * t121;
t67 = t80 * pkin(4) + pkin(3);
t83 = -pkin(10) - pkin(9);
t107 = t43 * t67 + t44 * t83 + t64;
t103 = -t70 * t120 - t31 * t69;
t106 = pkin(5) * t103 + t14 * pkin(11);
t17 = -t70 * t122 - t32 * t69;
t18 = t69 * t122 - t32 * t70;
t105 = -t17 * pkin(5) + t18 * pkin(11);
t36 = t44 * t69 + t74 * t70;
t37 = -t44 * t70 + t74 * t69;
t104 = t36 * pkin(5) + t37 * pkin(11);
t102 = t31 * t80 - t63;
t46 = t74 * t77 * pkin(2) + (-pkin(8) - qJ(3)) * t73;
t68 = t81 * pkin(2) + pkin(1);
t101 = -t78 * t46 + t82 * t68;
t60 = pkin(2) * t108;
t99 = -pkin(2) * t117 + t60;
t98 = pkin(5) * t70 + pkin(11) * t69;
t97 = g(1) * t103 + g(2) * t17;
t33 = t82 * t53 - t78 * t85;
t96 = g(1) * t30 - g(2) * t33;
t95 = g(1) * t82 + g(2) * t78;
t94 = g(1) * t78 - g(2) * t82;
t93 = -t82 * t46 - t78 * t68;
t92 = t31 * t76 + t109;
t91 = t30 * t67 - t31 * t83 + t99;
t90 = pkin(4) * t111 - t32 * t67 + t33 * t83 + t101;
t3 = g(1) * t17 - g(2) * t103 - g(3) * t36;
t5 = g(1) * t18 + g(2) * t14 + g(3) * t37;
t88 = g(1) * t33 + g(2) * t30 + g(3) * t43;
t87 = t49 * pkin(2);
t86 = pkin(4) * t63 - t30 * t83 - t31 * t67 + t93;
t84 = t32 * t83 + t33 * t67 + t87;
t65 = pkin(4) * t119;
t58 = pkin(4) * t110;
t51 = t95 * t73;
t50 = -t74 * t117 + t114;
t48 = -t74 * t115 - t116;
t47 = -t108 + t117;
t42 = -g(3) * t74 - t94 * t73;
t20 = -t32 * t80 + t111;
t19 = t110 + t126;
t9 = t18 * t79 - t33 * t75;
t8 = -t18 * t75 - t33 * t79;
t6 = t88 * t69;
t2 = t3 * t79;
t1 = t3 * t75;
t4 = [0, 0, 0, 0, 0, 0, t94, t95, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t48 - g(2) * t50, -g(1) * t47 - g(2) * t49, -t51, -g(1) * (-t78 * pkin(1) + pkin(8) * t120) - g(2) * (t82 * pkin(1) + pkin(8) * t122) 0, 0, 0, 0, 0, 0, g(1) * t31 + g(2) * t32, t96, -t51, -g(1) * t93 - g(2) * t101, 0, 0, 0, 0, 0, 0, g(1) * t102 - g(2) * t20, -g(1) * t92 - g(2) * t19, -t96, -g(1) * (-t31 * pkin(3) + t30 * pkin(9) + t93) - g(2) * (-pkin(3) * t32 - t33 * pkin(9) + t101) 0, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t97, -t96, -g(1) * t86 - g(2) * t90, 0, 0, 0, 0, 0, 0, g(1) * t132 - g(2) * t9, -g(1) * t133 - g(2) * t8, -t97, -g(1) * (-pkin(5) * t14 + pkin(11) * t103 + t86) - g(2) * (t18 * pkin(5) + t17 * pkin(11) + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t47 + t131, g(3) * t73 * t77 + g(1) * t50 - g(2) * t48, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t7, 0, -g(2) * t60 + (g(2) * t117 + t131) * pkin(2), 0, 0, 0, 0, 0, 0, -t88 * t80, t88 * t76, t7, -g(1) * (t33 * pkin(3) - t32 * pkin(9) + t87) - g(2) * (t30 * pkin(3) + pkin(9) * t31 + t99) - g(3) * (t43 * pkin(3) - t44 * pkin(9) + t64) 0, 0, 0, 0, 0, 0, -t88 * t70, t6, t7, -g(1) * t84 - g(2) * t91 - g(3) * t107, 0, 0, 0, 0, 0, 0, -g(1) * (t33 * t123 - t32 * t75) - g(2) * (t30 * t123 + t31 * t75) - g(3) * (t43 * t123 - t44 * t75) -g(1) * (-t33 * t124 - t32 * t79) - g(2) * (-t30 * t124 + t31 * t79) - g(3) * (-t43 * t124 - t44 * t79) -t6, -g(1) * (t98 * t33 + t84) - g(2) * (t98 * t30 + t91) - g(3) * (t98 * t43 + t107); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t19 + g(2) * t92 - g(3) * (t119 + t125) g(1) * t20 + g(2) * t102 - g(3) * (t44 * t80 - t74 * t76) 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t58 - g(3) * t65 + (g(2) * t109 - t7 * t76) * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (pkin(4) * t126 + t105 + t58) - g(2) * (-t92 * pkin(4) + t106) - g(3) * (pkin(4) * t125 + t104 + t65); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t105 - g(2) * t106 - g(3) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t133 - g(3) * (-t37 * t75 - t43 * t79) g(1) * t9 + g(2) * t132 - g(3) * (-t37 * t79 + t43 * t75) 0, 0;];
taug_reg  = t4;
