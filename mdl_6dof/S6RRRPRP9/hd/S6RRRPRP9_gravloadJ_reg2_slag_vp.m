% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t70 = cos(qJ(3));
t72 = cos(qJ(1));
t108 = t72 * t70;
t68 = sin(qJ(1));
t71 = cos(qJ(2));
t112 = t68 * t71;
t66 = sin(qJ(3));
t40 = t66 * t112 + t108;
t109 = t72 * t66;
t111 = t70 * t71;
t41 = t68 * t111 - t109;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t126 = -t40 * t69 + t41 * t65;
t67 = sin(qJ(2));
t114 = t67 * t70;
t116 = t66 * t69;
t28 = t65 * t114 - t67 * t116;
t42 = t71 * t109 - t68 * t70;
t43 = t71 * t108 + t68 * t66;
t89 = -t42 * t69 + t43 * t65;
t2 = g(1) * t89 + g(2) * t126 + g(3) * t28;
t16 = t42 * t65 + t43 * t69;
t131 = -pkin(5) * t89 + t16 * qJ(6);
t10 = t40 * t65 + t41 * t69;
t130 = -pkin(5) * t126 + t10 * qJ(6);
t88 = t65 * t66 + t69 * t70;
t127 = t88 * t67;
t81 = g(1) * t16 + g(2) * t10 + g(3) * t127;
t45 = g(1) * t72 + g(2) * t68;
t128 = t45 * t67;
t107 = t71 * pkin(2) + t67 * pkin(8);
t124 = -pkin(3) - pkin(4);
t123 = g(1) * t68;
t117 = t66 * t67;
t115 = t67 * t68;
t113 = t67 * t72;
t110 = t71 * t72;
t106 = t72 * pkin(1) + t68 * pkin(7);
t105 = qJ(4) * t66;
t104 = -pkin(2) - t105;
t103 = -t40 * pkin(3) + t41 * qJ(4);
t102 = -t42 * pkin(3) + t43 * qJ(4);
t101 = pkin(3) * t111 + t71 * t105 + t107;
t100 = pkin(2) * t110 + pkin(8) * t113 + t106;
t51 = pkin(8) * t112;
t99 = -pkin(9) * t112 + t51;
t55 = pkin(8) * t110;
t98 = -pkin(9) * t110 + t55;
t97 = -t40 * pkin(4) + t103;
t63 = t72 * pkin(7);
t96 = -t41 * pkin(3) - t40 * qJ(4) + t63;
t95 = -t42 * pkin(4) + t102;
t94 = pkin(4) * t111 + t101;
t93 = -g(1) * t126 + g(2) * t89;
t92 = g(1) * t40 - g(2) * t42;
t91 = -g(2) * t72 + t123;
t90 = -t28 * pkin(5) + qJ(6) * t127;
t48 = qJ(4) * t114;
t87 = t124 * t117 + t48;
t85 = -t41 * pkin(4) + pkin(9) * t115 + t96;
t84 = t67 * (t65 * t70 - t116);
t20 = t68 * t84;
t22 = t72 * t84;
t30 = t65 * t111 - t71 * t116;
t80 = g(1) * t22 + g(2) * t20 - g(3) * t30;
t78 = t43 * pkin(3) + t42 * qJ(4) + t100;
t77 = (-pkin(1) - t107) * t123;
t6 = g(1) * t42 + g(2) * t40 + g(3) * t117;
t76 = g(1) * t43 + g(2) * t41 + g(3) * t114;
t75 = -g(3) * t71 + t128;
t74 = t43 * pkin(4) - pkin(9) * t113 + t78;
t73 = (g(3) * pkin(9) + t45 * (-t124 * t70 - t104)) * t67;
t44 = g(1) * t115 - g(2) * t113;
t31 = t88 * t71;
t24 = g(3) * t67 + t45 * t71;
t23 = t72 * t127;
t21 = t68 * t127;
t19 = t75 * t70;
t18 = t75 * t66;
t17 = g(1) * t41 - g(2) * t43;
t4 = g(1) * t23 + g(2) * t21 - g(3) * t31;
t3 = g(1) * t10 - g(2) * t16;
t1 = [0, 0, 0, 0, 0, 0, t91, t45, 0, 0, 0, 0, 0, 0, 0, 0, t91 * t71, -t44, -t45, -g(1) * (-t68 * pkin(1) + t63) - g(2) * t106, 0, 0, 0, 0, 0, 0, t17, -t92, t44, -g(1) * t63 - g(2) * t100 - t77, 0, 0, 0, 0, 0, 0, t17, t44, t92, -g(1) * t96 - g(2) * t78 - t77, 0, 0, 0, 0, 0, 0, t3, t93, -t44, -g(1) * t85 - g(2) * t74 - t77, 0, 0, 0, 0, 0, 0, t3, -t44, -t93, -g(1) * (-pkin(5) * t10 - qJ(6) * t126 + t85) - g(2) * (t16 * pkin(5) + qJ(6) * t89 + t74) - t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t24, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, -t24, -g(1) * (-pkin(2) * t113 + t55) - g(2) * (-pkin(2) * t115 + t51) - g(3) * t107, 0, 0, 0, 0, 0, 0, t19, -t24, t18, -g(1) * t55 - g(2) * t51 - g(3) * t101 + (pkin(3) * t70 - t104) * t128, 0, 0, 0, 0, 0, 0, t4, -t80, t24, -g(1) * t98 - g(2) * t99 - g(3) * t94 + t73, 0, 0, 0, 0, 0, 0, t4, t24, t80, -g(1) * (-t23 * pkin(5) - t22 * qJ(6) + t98) - g(2) * (-t21 * pkin(5) - t20 * qJ(6) + t99) - g(3) * (t31 * pkin(5) + t30 * qJ(6) + t94) + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t76, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, -t76, -g(1) * t102 - g(2) * t103 - g(3) * (-pkin(3) * t117 + t48) 0, 0, 0, 0, 0, 0, -t2, -t81, 0, -g(1) * t95 - g(2) * t97 - g(3) * t87, 0, 0, 0, 0, 0, 0, -t2, 0, t81, -g(1) * (-t131 + t95) - g(2) * (-t130 + t97) - g(3) * (t87 - t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t81, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t81, -g(1) * t131 - g(2) * t130 - g(3) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg  = t1;
