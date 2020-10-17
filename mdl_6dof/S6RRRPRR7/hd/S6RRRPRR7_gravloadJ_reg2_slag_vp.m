% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:37:59
% EndTime: 2019-05-07 11:38:05
% DurationCPUTime: 0.98s
% Computational Cost: add. (842->157), mult. (1199->246), div. (0->0), fcn. (1423->14), ass. (0->85)
t70 = sin(pkin(6));
t79 = cos(qJ(1));
t109 = t70 * t79;
t74 = sin(qJ(2));
t75 = sin(qJ(1));
t78 = cos(qJ(2));
t101 = cos(pkin(6));
t91 = t79 * t101;
t43 = t74 * t91 + t75 * t78;
t69 = qJ(3) + pkin(12);
t65 = qJ(5) + t69;
t60 = sin(t65);
t61 = cos(t65);
t14 = -t60 * t109 + t43 * t61;
t42 = t75 * t74 - t78 * t91;
t72 = sin(qJ(6));
t76 = cos(qJ(6));
t122 = t14 * t72 - t42 * t76;
t121 = t14 * t76 + t42 * t72;
t113 = t70 * t74;
t77 = cos(qJ(3));
t111 = t70 * t77;
t92 = t75 * t101;
t45 = -t74 * t92 + t79 * t78;
t73 = sin(qJ(3));
t22 = t75 * t111 - t45 * t73;
t84 = t77 * t109 + t43 * t73;
t120 = -g(3) * (t101 * t77 - t73 * t113) + g(2) * t84 - g(1) * t22;
t118 = g(3) * t70;
t115 = t61 * t72;
t114 = t61 * t76;
t112 = t70 * t75;
t110 = t70 * t78;
t108 = t72 * t78;
t107 = t76 * t78;
t71 = -qJ(4) - pkin(9);
t64 = cos(t69);
t66 = t77 * pkin(3);
t52 = pkin(4) * t64 + t66;
t50 = pkin(2) + t52;
t68 = -pkin(10) + t71;
t106 = -t42 * t50 - t43 * t68;
t44 = t79 * t74 + t78 * t92;
t105 = -t44 * t50 - t45 * t68;
t63 = sin(t69);
t51 = t73 * pkin(3) + pkin(4) * t63;
t104 = t52 * t112 - t45 * t51;
t103 = t101 * t52 - t51 * t113;
t102 = t79 * pkin(1) + pkin(8) * t112;
t100 = t73 * t112;
t55 = t73 * t109;
t99 = -t75 * pkin(1) + pkin(8) * t109;
t95 = -t61 * t109 - t43 * t60;
t98 = pkin(5) * t95 + t14 * pkin(11);
t17 = -t61 * t112 + t45 * t60;
t18 = t60 * t112 + t45 * t61;
t97 = -t17 * pkin(5) + t18 * pkin(11);
t31 = t101 * t61 - t60 * t113;
t32 = t101 * t60 + t61 * t113;
t96 = t31 * pkin(5) + t32 * pkin(11);
t94 = -t63 * t109 + t43 * t64;
t93 = t43 * t77 - t55;
t90 = -t52 * t109 - t43 * t51;
t89 = t51 * t112 - t44 * t68 + t45 * t50 + t102;
t88 = pkin(5) * t61 + pkin(11) * t60;
t87 = g(1) * t95 + g(2) * t17;
t19 = g(1) * t42 - g(2) * t44;
t86 = g(1) * t79 + g(2) * t75;
t85 = t64 * t109 + t43 * t63;
t83 = t51 * t109 + t42 * t68 - t43 * t50 + t99;
t3 = g(1) * t17 - g(2) * t95 - g(3) * t31;
t5 = g(1) * t18 + g(2) * t14 + g(3) * t32;
t9 = -g(1) * t44 - g(2) * t42 + g(3) * t110;
t81 = g(1) * t45 + g(2) * t43 + g(3) * t113;
t62 = t66 + pkin(2);
t37 = t50 * t110;
t23 = t45 * t77 + t100;
t21 = t63 * t112 + t45 * t64;
t20 = t64 * t112 - t45 * t63;
t8 = t18 * t76 + t44 * t72;
t7 = -t18 * t72 + t44 * t76;
t6 = t9 * t60;
t2 = t3 * t76;
t1 = t3 * t72;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t75 - g(2) * t79, t86, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t43 - g(2) * t45, -t19, -t86 * t70, -g(1) * t99 - g(2) * t102, 0, 0, 0, 0, 0, 0, g(1) * t93 - g(2) * t23, -g(1) * t84 - g(2) * t22, t19, -g(1) * (-t43 * pkin(2) - t42 * pkin(9) + t99) - g(2) * (t45 * pkin(2) + t44 * pkin(9) + t102) 0, 0, 0, 0, 0, 0, g(1) * t94 - g(2) * t21, -g(1) * t85 - g(2) * t20, t19, -g(1) * (pkin(3) * t55 + t42 * t71 - t43 * t62 + t99) - g(2) * (pkin(3) * t100 - t44 * t71 + t45 * t62 + t102) 0, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t87, t19, -g(1) * t83 - g(2) * t89, 0, 0, 0, 0, 0, 0, g(1) * t121 - g(2) * t8, -g(1) * t122 - g(2) * t7, -t87, -g(1) * (-pkin(5) * t14 + pkin(11) * t95 + t83) - g(2) * (t18 * pkin(5) + t17 * pkin(11) + t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t81, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t77, t9 * t73, -t81, -g(1) * (-t44 * pkin(2) + t45 * pkin(9)) - g(2) * (-t42 * pkin(2) + t43 * pkin(9)) - (pkin(2) * t78 + pkin(9) * t74) * t118, 0, 0, 0, 0, 0, 0, -t9 * t64, t9 * t63, -t81, -g(1) * (-t44 * t62 - t45 * t71) - g(2) * (-t42 * t62 - t43 * t71) - (t62 * t78 - t71 * t74) * t118, 0, 0, 0, 0, 0, 0, -t9 * t61, t6, -t81, -g(1) * t105 - g(2) * t106 - g(3) * (-t68 * t113 + t37) 0, 0, 0, 0, 0, 0, -g(1) * (-t44 * t114 + t45 * t72) - g(2) * (-t42 * t114 + t43 * t72) - (t61 * t107 + t72 * t74) * t118, -g(1) * (t44 * t115 + t45 * t76) - g(2) * (t42 * t115 + t43 * t76) - (-t61 * t108 + t74 * t76) * t118, -t6, -g(1) * (-t88 * t44 + t105) - g(2) * (-t88 * t42 + t106) - g(3) * t37 - (-t68 * t74 + t88 * t78) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, g(1) * t23 + g(2) * t93 - g(3) * (-t101 * t73 - t74 * t111) 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + g(2) * t85 - g(3) * (t101 * t64 - t63 * t113) g(1) * t21 + g(2) * t94 - g(3) * (-t101 * t63 - t64 * t113) 0, t120 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t104 - g(2) * t90 - g(3) * t103, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (t97 + t104) - g(2) * (t90 + t98) - g(3) * (t96 + t103); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t97 - g(2) * t98 - g(3) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t122 - g(3) * (-t70 * t107 - t32 * t72) g(1) * t8 + g(2) * t121 - g(3) * (t70 * t108 - t32 * t76) 0, 0;];
taug_reg  = t4;
