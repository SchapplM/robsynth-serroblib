% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t117 = cos(pkin(12));
t118 = cos(pkin(7));
t114 = sin(pkin(12));
t123 = sin(qJ(2));
t125 = cos(qJ(2));
t119 = cos(pkin(6));
t97 = t119 * t117;
t73 = t114 * t123 - t125 * t97;
t115 = sin(pkin(7));
t116 = sin(pkin(6));
t93 = t116 * t115;
t128 = t117 * t93 + t73 * t118;
t95 = t119 * t114;
t74 = t117 * t123 + t125 * t95;
t92 = t116 * t114;
t127 = -t115 * t92 + t74 * t118;
t94 = t118 * t116;
t126 = t119 * t115 + t125 * t94;
t124 = cos(qJ(3));
t63 = sin(qJ(5));
t67 = cos(qJ(4));
t122 = t63 * t67;
t66 = cos(qJ(5));
t121 = t66 * t67;
t101 = t116 * t125;
t79 = t123 * t93;
t120 = pkin(2) * t101 + pkin(9) * t79;
t65 = sin(qJ(3));
t81 = t123 * t94;
t47 = t124 * t101 - t65 * t81;
t113 = t47 * pkin(3) + t120;
t112 = pkin(5) * t63 + pkin(10);
t51 = t114 * t125 + t123 * t97;
t22 = t128 * t124 + t51 * t65;
t20 = t22 * pkin(3);
t23 = t51 * t124 - t128 * t65;
t111 = t23 * pkin(10) - t20;
t52 = t117 * t125 - t123 * t95;
t24 = t127 * t124 + t52 * t65;
t21 = t24 * pkin(3);
t25 = t52 * t124 - t127 * t65;
t110 = t25 * pkin(10) - t21;
t100 = t116 * t123;
t39 = t65 * t100 - t126 * t124;
t38 = t39 * pkin(3);
t40 = t124 * t100 + t126 * t65;
t109 = t40 * pkin(10) - t38;
t108 = t51 * t115;
t107 = t52 * t115;
t64 = sin(qJ(4));
t106 = t64 * t115;
t105 = t65 * t118;
t104 = t67 * t115;
t103 = -pkin(4) * t67 - pkin(11) * t64;
t102 = t118 * t124;
t61 = t66 * pkin(5) + pkin(4);
t62 = -qJ(6) - pkin(11);
t99 = -t61 * t67 + t62 * t64;
t46 = t65 * t101 + t124 * t81;
t98 = t46 * pkin(10) + t113;
t91 = -t73 * pkin(2) + pkin(9) * t108;
t90 = -t74 * pkin(2) + pkin(9) * t107;
t31 = -t51 * t105 - t73 * t124;
t89 = t31 * pkin(3) + t91;
t33 = -t52 * t105 - t74 * t124;
t88 = t33 * pkin(3) + t90;
t68 = t73 * t115 - t117 * t94;
t12 = t23 * t64 - t68 * t67;
t69 = t74 * t115 + t118 * t92;
t14 = t25 * t64 - t69 * t67;
t72 = t119 * t118 - t125 * t93;
t26 = t40 * t64 - t72 * t67;
t87 = g(1) * t14 + g(2) * t12 + g(3) * t26;
t13 = t23 * t67 + t68 * t64;
t15 = t25 * t67 + t69 * t64;
t27 = t40 * t67 + t72 * t64;
t86 = g(1) * t15 + g(2) * t13 + g(3) * t27;
t16 = -t51 * t104 + t31 * t64;
t18 = -t52 * t104 + t33 * t64;
t34 = t47 * t64 - t67 * t79;
t85 = g(1) * t18 + g(2) * t16 + g(3) * t34;
t84 = g(1) * t24 + g(2) * t22 + g(3) * t39;
t83 = g(1) * t25 + g(2) * t23 + g(3) * t40;
t30 = t51 * t102 - t73 * t65;
t32 = t52 * t102 - t74 * t65;
t82 = g(1) * t32 + g(2) * t30 + g(3) * t46;
t76 = t30 * pkin(10) + t89;
t75 = t32 * pkin(10) + t88;
t1 = -g(1) * (-t15 * t63 + t24 * t66) - g(2) * (-t13 * t63 + t22 * t66) - g(3) * (-t27 * t63 + t39 * t66);
t35 = t47 * t67 + t64 * t79;
t19 = t52 * t106 + t33 * t67;
t17 = t51 * t106 + t31 * t67;
t11 = t84 * t64;
t8 = t87 * t66;
t7 = t87 * t63;
t6 = -g(1) * (t19 * t66 + t32 * t63) - g(2) * (t17 * t66 + t30 * t63) - g(3) * (t35 * t66 + t46 * t63);
t5 = -g(1) * (-t19 * t63 + t32 * t66) - g(2) * (-t17 * t63 + t30 * t66) - g(3) * (-t35 * t63 + t46 * t66);
t4 = -g(1) * (-t24 * t121 + t25 * t63) - g(2) * (-t22 * t121 + t23 * t63) - g(3) * (-t39 * t121 + t40 * t63);
t3 = -g(1) * (t24 * t122 + t25 * t66) - g(2) * (t22 * t122 + t23 * t66) - g(3) * (t39 * t122 + t40 * t66);
t2 = -g(1) * (-t15 * t66 - t24 * t63) - g(2) * (-t13 * t66 - t22 * t63) - g(3) * (-t27 * t66 - t39 * t63);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t74 + g(2) * t73 - g(3) * t101, g(1) * t52 + g(2) * t51 + g(3) * t100, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t33 - g(2) * t31 - g(3) * t47, t82, -g(1) * t107 - g(2) * t108 - g(3) * t79, -g(1) * t90 - g(2) * t91 - g(3) * t120, 0, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t17 - g(3) * t35, t85, -t82, -g(1) * t75 - g(2) * t76 - g(3) * t98, 0, 0, 0, 0, 0, 0, t6, t5, -t85, -g(1) * (t19 * pkin(4) + t18 * pkin(11) + t75) - g(2) * (t17 * pkin(4) + t16 * pkin(11) + t76) - g(3) * (t35 * pkin(4) + t34 * pkin(11) + t98) 0, 0, 0, 0, 0, 0, t6, t5, -t85, -g(1) * (t112 * t32 - t18 * t62 + t19 * t61 + t88) - g(2) * (t112 * t30 - t16 * t62 + t17 * t61 + t89) - g(3) * (t112 * t46 - t34 * t62 + t35 * t61 + t113); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t83, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t67, -t11, -t83, -g(1) * t110 - g(2) * t111 - g(3) * t109, 0, 0, 0, 0, 0, 0, t4, t3, t11, -g(1) * (t103 * t24 + t110) - g(2) * (t103 * t22 + t111) - g(3) * (t103 * t39 + t109) 0, 0, 0, 0, 0, 0, t4, t3, t11, -g(1) * (t112 * t25 + t99 * t24 - t21) - g(2) * (t112 * t23 + t99 * t22 - t20) - g(3) * (t112 * t40 + t99 * t39 - t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t86, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t86, -g(1) * (-t14 * pkin(4) + t15 * pkin(11)) - g(2) * (-t12 * pkin(4) + t13 * pkin(11)) - g(3) * (-t26 * pkin(4) + t27 * pkin(11)) 0, 0, 0, 0, 0, 0, t8, -t7, -t86, -g(1) * (-t14 * t61 - t15 * t62) - g(2) * (-t12 * t61 - t13 * t62) - g(3) * (-t26 * t61 - t27 * t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87;];
taug_reg  = t9;
