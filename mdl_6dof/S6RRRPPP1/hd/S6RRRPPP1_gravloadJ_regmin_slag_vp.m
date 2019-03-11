% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t81 = cos(pkin(6));
t83 = sin(qJ(2));
t142 = t83 * t81;
t79 = sin(pkin(6));
t86 = cos(qJ(2));
t147 = t79 * t86;
t82 = sin(qJ(3));
t102 = t82 * t142 + t147;
t85 = cos(qJ(3));
t140 = t83 * t85;
t78 = sin(pkin(10));
t80 = cos(pkin(10));
t157 = -t102 * t78 + t80 * t140;
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t109 = g(1) * t87 + g(2) * t84;
t158 = t109 * t83;
t93 = -g(3) * t86 + t158;
t133 = t87 * t85;
t52 = t86 * t133 + t84 * t82;
t156 = g(1) * t52;
t155 = g(1) * t84;
t134 = t87 * t82;
t137 = t85 * t86;
t50 = t84 * t137 - t134;
t153 = g(2) * t50;
t75 = t86 * pkin(2);
t150 = t78 * t81;
t149 = t78 * t82;
t148 = t79 * t83;
t146 = t80 * t81;
t145 = t80 * t82;
t144 = t81 * t85;
t143 = t82 * t83;
t141 = t83 * t84;
t139 = t83 * t87;
t138 = t84 * t86;
t136 = t86 * t81;
t135 = t86 * t87;
t128 = qJ(4) * t81;
t119 = t86 * t128;
t131 = pkin(9) * t138 + t84 * t119;
t130 = pkin(9) * t135 + t87 * t119;
t129 = qJ(4) * t79;
t127 = qJ(4) * t83;
t126 = g(3) * t140;
t124 = t78 * t148;
t123 = t80 * t148;
t120 = t82 * t129;
t118 = t50 * t129;
t117 = t52 * t129;
t66 = t81 * t127;
t115 = t85 * t79 * t127 - pkin(3) * t143;
t49 = t82 * t138 + t133;
t20 = t50 * t146 - t49 * t78;
t21 = -t50 * t150 - t49 * t80;
t44 = t49 * pkin(3);
t114 = t21 * pkin(4) + t20 * qJ(5) - t44;
t51 = -t86 * t134 + t84 * t85;
t22 = t52 * t146 + t51 * t78;
t23 = -t52 * t150 + t51 * t80;
t46 = t51 * pkin(3);
t113 = t23 * pkin(4) + t22 * qJ(5) + t46;
t112 = pkin(3) * t137 + t83 * pkin(9) + t86 * t120 + t66 + t75;
t13 = t84 * t123 - t49 * t146 - t50 * t78;
t15 = -t87 * t123 - t51 * t146 + t52 * t78;
t111 = g(1) * t13 + g(2) * t15;
t14 = -t84 * t124 + t49 * t150 - t50 * t80;
t16 = t52 * t80 + (t79 * t139 + t51 * t81) * t78;
t110 = g(1) * t14 + g(2) * t16;
t108 = -g(2) * t87 + t155;
t107 = -t50 * pkin(3) + t87 * pkin(8) - t49 * t129;
t90 = t102 * t80 + t78 * t140;
t26 = t90 * t84;
t27 = t157 * t84;
t106 = -t27 * pkin(4) - t26 * qJ(5) + t131;
t28 = t90 * t87;
t29 = t157 * t87;
t105 = -t29 * pkin(4) - t28 * qJ(5) + t130;
t104 = -t126 - t153;
t103 = t81 * t141 + t49 * t79;
t101 = t79 * t143 - t136;
t40 = (-t80 * t144 + t149) * t83;
t99 = g(1) * t22 + g(2) * t20 - g(3) * t40;
t41 = (t78 * t144 + t145) * t83;
t98 = g(1) * t23 + g(2) * t21 - g(3) * t41;
t31 = -t123 + (t81 * t145 + t78 * t85) * t86;
t2 = g(1) * t28 + g(2) * t26 - g(3) * t31;
t32 = -t136 * t149 + t80 * t137 + t124;
t3 = g(1) * t29 + g(2) * t27 - g(3) * t32;
t97 = -t41 * pkin(4) - t40 * qJ(5) + t115;
t96 = -t104 + t156;
t95 = t32 * pkin(4) + t31 * qJ(5) + t112;
t94 = t14 * pkin(4) + t13 * qJ(5) + t107;
t92 = pkin(2) * t135 + t52 * pkin(3) + t84 * pkin(8) + pkin(9) * t139 - t51 * t129 + (pkin(1) + t66) * t87;
t91 = (-t75 - pkin(1) + (-pkin(9) - t128) * t83) * t155;
t89 = t16 * pkin(4) + t15 * qJ(5) + t92;
t88 = (pkin(3) * t85 + pkin(2) + t120) * t158;
t48 = t82 * t147 + t142;
t43 = t101 * t87;
t42 = t101 * t84;
t34 = t81 * t139 - t51 * t79;
t17 = t96 * t79;
t10 = g(1) * t103 - g(2) * t34;
t9 = g(1) * t43 + g(2) * t42 - g(3) * t48;
t8 = -g(1) * t34 - g(2) * t103 - g(3) * t101;
t1 = -g(1) * t15 + t104 * t78 + (-g(2) * (-t79 * t141 + t49 * t81) - g(3) * t102) * t80;
t4 = [0, t108, t109, 0, 0, 0, 0, 0, t108 * t86, -t108 * t83, 0, 0, 0, 0, 0, g(1) * t50 - g(2) * t52, -g(1) * t49 - g(2) * t51, -t110, t111, t10, -g(1) * t107 - g(2) * t92 - t91, t10, t110, -t111, -g(1) * t94 - g(2) * t89 - t91, t10, -t111, -t110, -g(1) * (-pkin(5) * t103 + t14 * qJ(6) + t94) - g(2) * (t34 * pkin(5) + t16 * qJ(6) + t89) - t91; 0, 0, 0, 0, 0, 0, 0, 0, t93, g(3) * t83 + t109 * t86, 0, 0, 0, 0, 0, t93 * t85, -t93 * t82, t3, -t2, t9, -g(1) * t130 - g(2) * t131 - g(3) * t112 + t88, t9, -t3, t2, -g(1) * t105 - g(2) * t106 - g(3) * t95 + t88, t9, t2, t3, -g(1) * (-t43 * pkin(5) - t29 * qJ(6) + t105) - g(2) * (-t42 * pkin(5) - t27 * qJ(6) + t106) - g(3) * (t48 * pkin(5) + t32 * qJ(6) + t95) + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t51 + g(2) * t49 + g(3) * t143, t96, -t98, t99, -t17, -g(1) * (t46 + t117) - g(2) * (-t44 + t118) - g(3) * t115, -t17, t98, -t99, -g(1) * (t113 + t117) - g(2) * (t114 + t118) - g(3) * t97, -t17, -t99, -t98, -g(1) * (t23 * qJ(6) + t113) - g(2) * (t21 * qJ(6) + t114) - g(3) * (-t41 * qJ(6) + t97) + (-pkin(5) * t126 + (-t153 - t156) * (pkin(5) + qJ(4))) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, t8, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t14 - g(3) * t157;];
taug_reg  = t4;
