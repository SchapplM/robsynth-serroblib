% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t64 = sin(qJ(2));
t65 = sin(qJ(1));
t68 = cos(qJ(2));
t116 = cos(qJ(1));
t95 = cos(pkin(6));
t78 = t95 * t116;
t44 = t64 * t78 + t65 * t68;
t63 = sin(qJ(3));
t67 = cos(qJ(3));
t61 = sin(pkin(6));
t92 = t61 * t116;
t23 = t44 * t63 + t67 * t92;
t43 = t65 * t64 - t68 * t78;
t62 = sin(qJ(6));
t66 = cos(qJ(6));
t121 = t23 * t62 + t43 * t66;
t120 = t23 * t66 - t43 * t62;
t113 = t43 * t67;
t99 = qJ(4) * t63;
t119 = -pkin(3) * t113 - t43 * t99;
t88 = t65 * t95;
t45 = t116 * t64 + t68 * t88;
t112 = t45 * t67;
t118 = -pkin(3) * t112 - t45 * t99;
t117 = g(3) * t61;
t111 = t61 * t64;
t110 = t61 * t65;
t109 = t61 * t67;
t108 = t61 * t68;
t107 = t62 * t63;
t106 = t62 * t68;
t105 = t63 * t66;
t104 = t66 * t68;
t103 = pkin(9) - qJ(5);
t102 = qJ(4) + pkin(5);
t101 = pkin(2) * t108 + pkin(9) * t111;
t100 = t116 * pkin(1) + pkin(8) * t110;
t98 = qJ(5) * t64;
t97 = t23 * qJ(4);
t46 = t116 * t68 - t64 * t88;
t27 = -t65 * t109 + t46 * t63;
t96 = t27 * qJ(4);
t94 = t67 * t108;
t93 = t46 * pkin(2) + t100;
t91 = -t65 * pkin(1) + pkin(8) * t92;
t37 = t43 * pkin(2);
t90 = t44 * pkin(9) - t37;
t39 = t45 * pkin(2);
t89 = t46 * pkin(9) - t39;
t24 = t44 * t67 - t63 * t92;
t16 = t23 * pkin(3);
t87 = t24 * qJ(4) - t16;
t20 = t27 * pkin(3);
t28 = t63 * t110 + t46 * t67;
t86 = t28 * qJ(4) - t20;
t41 = t63 * t111 - t95 * t67;
t36 = t41 * pkin(3);
t42 = t64 * t109 + t95 * t63;
t85 = t42 * qJ(4) - t36;
t84 = pkin(3) * t94 + t99 * t108 + t101;
t83 = -t44 * pkin(2) + t91;
t82 = pkin(4) * t94 + t84;
t81 = pkin(5) * t63 + pkin(10) * t67;
t80 = -g(1) * t23 + g(2) * t27;
t79 = -g(1) * t24 + g(2) * t28;
t14 = g(1) * t43 - g(2) * t45;
t77 = t45 * pkin(9) + t93;
t76 = g(1) * t116 + g(2) * t65;
t75 = -t43 * pkin(9) + t83;
t2 = g(1) * t27 + g(2) * t23 + g(3) * t41;
t74 = g(1) * t28 + g(2) * t24 + g(3) * t42;
t73 = -pkin(4) * t113 + t103 * t44 + t119 - t37;
t72 = -pkin(4) * t112 + t103 * t46 + t118 - t39;
t71 = -g(1) * t45 - g(2) * t43 + g(3) * t108;
t12 = g(1) * t46 + g(2) * t44 + g(3) * t111;
t22 = t28 * pkin(3);
t70 = t28 * pkin(4) + t103 * t45 + t22 + t93;
t18 = t24 * pkin(3);
t69 = -pkin(4) * t24 - t103 * t43 - t18 + t83;
t35 = t41 * pkin(4);
t19 = t27 * pkin(4);
t15 = t23 * pkin(4);
t9 = t71 * t67;
t8 = t71 * t63;
t7 = t27 * t66 - t45 * t62;
t6 = -t27 * t62 - t45 * t66;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t65 - g(2) * t116, t76, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t46, -t14, -t76 * t61, -g(1) * t91 - g(2) * t100, 0, 0, 0, 0, 0, 0, -t79, t80, t14, -g(1) * t75 - g(2) * t77, 0, 0, 0, 0, 0, 0, -t79, t14, -t80, -g(1) * (-t18 + t75 - t97) - g(2) * (t22 + t77 + t96) 0, 0, 0, 0, 0, 0, -t80, t79, -t14, -g(1) * (t69 - t97) - g(2) * (t70 + t96) 0, 0, 0, 0, 0, 0, g(1) * t120 - g(2) * t7, -g(1) * t121 - g(2) * t6, -t79, -g(1) * (-pkin(10) * t24 - t102 * t23 + t69) - g(2) * (t28 * pkin(10) + t102 * t27 + t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t12, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, -t12, -g(1) * t89 - g(2) * t90 - g(3) * t101, 0, 0, 0, 0, 0, 0, -t9, -t12, -t8, -g(1) * (t89 + t118) - g(2) * (t90 + t119) - g(3) * t84, 0, 0, 0, 0, 0, 0, -t8, t9, t12, -g(1) * t72 - g(2) * t73 - g(3) * (-t61 * t98 + t82) 0, 0, 0, 0, 0, 0, -g(1) * (-t45 * t105 - t46 * t62) - g(2) * (-t43 * t105 - t44 * t62) - (t63 * t104 - t62 * t64) * t117, -g(1) * (t45 * t107 - t46 * t66) - g(2) * (t43 * t107 - t44 * t66) - (-t63 * t106 - t64 * t66) * t117, -t9, -g(1) * (-t81 * t45 + t72) - g(2) * (-t81 * t43 + t73) - g(3) * t82 - (t81 * t68 - t98) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t74, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t74, -g(1) * t86 - g(2) * t87 - g(3) * t85, 0, 0, 0, 0, 0, 0, -t74, -t2, 0, -g(1) * (-t19 + t86) - g(2) * (-t15 + t87) - g(3) * (-t35 + t85) 0, 0, 0, 0, 0, 0, -t74 * t66, t74 * t62, t2, -g(1) * (-t27 * pkin(10) + t102 * t28 - t19 - t20) - g(2) * (-t23 * pkin(10) + t102 * t24 - t15 - t16) - g(3) * (-t41 * pkin(10) + t102 * t42 - t35 - t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t121 - g(3) * (t61 * t104 - t41 * t62) g(1) * t7 + g(2) * t120 - g(3) * (-t61 * t106 - t41 * t66) 0, 0;];
taug_reg  = t1;
