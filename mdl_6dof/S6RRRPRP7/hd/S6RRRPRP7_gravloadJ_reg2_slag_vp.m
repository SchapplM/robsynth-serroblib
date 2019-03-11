% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP7
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t83 = qJ(3) + pkin(11);
t80 = sin(t83);
t81 = cos(t83);
t143 = pkin(4) * t81 + pkin(10) * t80;
t84 = sin(pkin(6));
t93 = cos(qJ(1));
t131 = t84 * t93;
t126 = cos(pkin(6));
t116 = t93 * t126;
t88 = sin(qJ(2));
t89 = sin(qJ(1));
t92 = cos(qJ(2));
t59 = t88 * t116 + t89 * t92;
t29 = -t80 * t131 + t59 * t81;
t58 = -t92 * t116 + t89 * t88;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t9 = t29 * t86 - t58 * t90;
t10 = t29 * t90 + t58 * t86;
t117 = t89 * t126;
t61 = -t88 * t117 + t93 * t92;
t87 = sin(qJ(3));
t138 = t61 * t87;
t137 = t81 * t86;
t136 = t81 * t90;
t135 = t84 * t88;
t134 = t84 * t89;
t91 = cos(qJ(3));
t133 = t84 * t91;
t132 = t84 * t92;
t130 = t90 * t92;
t79 = t91 * pkin(3) + pkin(2);
t85 = -qJ(4) - pkin(9);
t129 = -t58 * t79 - t59 * t85;
t60 = t92 * t117 + t93 * t88;
t128 = -t60 * t79 - t61 * t85;
t127 = t93 * pkin(1) + pkin(8) * t134;
t125 = t87 * t135;
t124 = t86 * t132;
t123 = t87 * t134;
t122 = t89 * t133;
t74 = t87 * t131;
t121 = t91 * t131;
t120 = -t89 * pkin(1) + pkin(8) * t131;
t119 = -t81 * t131 - t59 * t80;
t118 = t59 * t91 - t74;
t115 = t126 * t91;
t114 = -t143 * t58 + t129;
t113 = -t143 * t60 + t128;
t112 = t79 * t132 - t85 * t135;
t111 = pkin(3) * t123 - t60 * t85 + t61 * t79 + t127;
t33 = t80 * t134 + t61 * t81;
t13 = t33 * t86 - t60 * t90;
t110 = -g(1) * t9 + g(2) * t13;
t32 = -t81 * t134 + t61 * t80;
t109 = g(1) * t119 + g(2) * t32;
t21 = g(1) * t58 - g(2) * t60;
t108 = g(1) * t93 + g(2) * t89;
t107 = pkin(5) * t90 + qJ(6) * t86;
t106 = t59 * t87 + t121;
t105 = pkin(3) * t74 + t58 * t85 - t59 * t79 + t120;
t104 = t143 * t132 + t112;
t48 = t126 * t80 + t81 * t135;
t26 = t84 * t130 + t48 * t86;
t1 = g(1) * t13 + g(2) * t9 + g(3) * t26;
t14 = t33 * t90 + t60 * t86;
t27 = t48 * t90 - t124;
t103 = g(1) * t14 + g(2) * t10 + g(3) * t27;
t15 = -t58 * t137 - t59 * t90;
t17 = -t60 * t137 - t61 * t90;
t36 = t81 * t124 - t90 * t135;
t102 = g(1) * t17 + g(2) * t15 + g(3) * t36;
t47 = t126 * t81 - t80 * t135;
t101 = g(1) * t32 - g(2) * t119 - g(3) * t47;
t100 = g(1) * t33 + g(2) * t29 + g(3) * t48;
t69 = pkin(3) * t122;
t99 = -pkin(3) * t138 - t32 * pkin(4) + t33 * pkin(10) + t69;
t19 = -g(1) * t60 - g(2) * t58 + g(3) * t132;
t98 = g(1) * t61 + g(2) * t59 + g(3) * t135;
t97 = t33 * pkin(4) + t32 * pkin(10) + t111;
t78 = pkin(3) * t115;
t96 = -pkin(3) * t125 + t47 * pkin(4) + t48 * pkin(10) + t78;
t95 = -pkin(4) * t29 + pkin(10) * t119 + t105;
t94 = -t106 * pkin(3) + pkin(4) * t119 + t29 * pkin(10);
t37 = (t81 * t130 + t86 * t88) * t84;
t35 = t61 * t91 + t123;
t34 = t122 - t138;
t18 = -t60 * t136 + t61 * t86;
t16 = -t58 * t136 + t59 * t86;
t8 = t19 * t80;
t5 = t101 * t90;
t4 = t101 * t86;
t3 = -g(1) * t18 - g(2) * t16 - g(3) * t37;
t2 = g(1) * t10 - g(2) * t14;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t89 - g(2) * t93, t108, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * t61, -t21, -t108 * t84, -g(1) * t120 - g(2) * t127, 0, 0, 0, 0, 0, 0, g(1) * t118 - g(2) * t35, -g(1) * t106 - g(2) * t34, t21, -g(1) * (-t59 * pkin(2) - t58 * pkin(9) + t120) - g(2) * (t61 * pkin(2) + t60 * pkin(9) + t127) 0, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t33, t109, t21, -g(1) * t105 - g(2) * t111, 0, 0, 0, 0, 0, 0, t2, t110, -t109, -g(1) * t95 - g(2) * t97, 0, 0, 0, 0, 0, 0, t2, -t109, -t110, -g(1) * (-pkin(5) * t10 - qJ(6) * t9 + t95) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t97); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t98, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t91, t19 * t87, -t98, -g(1) * (-t60 * pkin(2) + t61 * pkin(9)) - g(2) * (-t58 * pkin(2) + t59 * pkin(9)) - g(3) * (pkin(2) * t92 + pkin(9) * t88) * t84, 0, 0, 0, 0, 0, 0, -t19 * t81, t8, -t98, -g(1) * t128 - g(2) * t129 - g(3) * t112, 0, 0, 0, 0, 0, 0, t3, t102, -t8, -g(1) * t113 - g(2) * t114 - g(3) * t104, 0, 0, 0, 0, 0, 0, t3, -t8, -t102, -g(1) * (t18 * pkin(5) + t17 * qJ(6) + t113) - g(2) * (t16 * pkin(5) + t15 * qJ(6) + t114) - g(3) * (t37 * pkin(5) + t36 * qJ(6) + t104); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t34 + g(2) * t106 - g(3) * (t115 - t125) g(1) * t35 + g(2) * t118 - g(3) * (-t126 * t87 - t88 * t133) 0, 0, 0, 0, 0, 0, 0, 0, t101, t100, 0, -g(1) * t69 - g(3) * t78 + (g(2) * t121 + t98 * t87) * pkin(3), 0, 0, 0, 0, 0, 0, t5, -t4, -t100, -g(1) * t99 - g(2) * t94 - g(3) * t96, 0, 0, 0, 0, 0, 0, t5, -t100, t4, -g(1) * (-t107 * t32 + t99) - g(2) * (t107 * t119 + t94) - g(3) * (t107 * t47 + t96); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t103, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t103, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t26 * pkin(5) + t27 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
