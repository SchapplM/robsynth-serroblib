% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t84 = cos(pkin(6));
t86 = sin(qJ(2));
t148 = t86 * t84;
t82 = sin(pkin(6));
t89 = cos(qJ(2));
t153 = t82 * t89;
t85 = sin(qJ(3));
t105 = t85 * t148 + t153;
t88 = cos(qJ(3));
t146 = t86 * t88;
t81 = sin(pkin(10));
t83 = cos(pkin(10));
t163 = -t105 * t81 + t83 * t146;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t68 = g(1) * t90 + g(2) * t87;
t164 = t68 * t86;
t96 = -g(3) * t89 + t164;
t139 = t90 * t88;
t53 = t139 * t89 + t87 * t85;
t162 = g(1) * t53;
t161 = g(1) * t87;
t140 = t90 * t85;
t143 = t88 * t89;
t51 = t143 * t87 - t140;
t159 = g(2) * t51;
t76 = t86 * pkin(9);
t78 = t89 * pkin(2);
t156 = t81 * t84;
t155 = t81 * t85;
t154 = t82 * t86;
t152 = t83 * t84;
t151 = t83 * t85;
t150 = t84 * t88;
t149 = t85 * t86;
t147 = t86 * t87;
t145 = t86 * t90;
t144 = t87 * t89;
t142 = t89 * t84;
t141 = t89 * t90;
t132 = qJ(4) * t84;
t123 = t89 * t132;
t70 = pkin(9) * t144;
t137 = t87 * t123 + t70;
t73 = pkin(9) * t141;
t136 = t90 * t123 + t73;
t135 = t78 + t76;
t134 = t90 * pkin(1) + t87 * pkin(8);
t133 = qJ(4) * t82;
t131 = qJ(4) * t86;
t130 = g(3) * t146;
t128 = t81 * t154;
t127 = t83 * t154;
t124 = t85 * t133;
t122 = t51 * t133;
t121 = t53 * t133;
t69 = t84 * t131;
t120 = -pkin(1) - t78;
t118 = pkin(2) * t141 + pkin(9) * t145 + t134;
t117 = t88 * t82 * t131 - pkin(3) * t149;
t50 = t144 * t85 + t139;
t20 = t152 * t51 - t50 * t81;
t21 = -t156 * t51 - t50 * t83;
t45 = t50 * pkin(3);
t116 = t21 * pkin(4) + t20 * qJ(5) - t45;
t52 = -t140 * t89 + t87 * t88;
t22 = t152 * t53 + t52 * t81;
t23 = -t156 * t53 + t52 * t83;
t47 = t52 * pkin(3);
t115 = t23 * pkin(4) + t22 * qJ(5) + t47;
t114 = pkin(3) * t143 + t89 * t124 + t135 + t69;
t13 = t127 * t87 - t152 * t50 - t51 * t81;
t15 = -t127 * t90 - t152 * t52 + t53 * t81;
t113 = g(1) * t13 + g(2) * t15;
t14 = -t87 * t128 + t156 * t50 - t51 * t83;
t16 = t53 * t83 + (t145 * t82 + t52 * t84) * t81;
t112 = g(1) * t14 + g(2) * t16;
t111 = -g(2) * t90 + t161;
t79 = t90 * pkin(8);
t110 = -t51 * pkin(3) - t133 * t50 + t79;
t93 = t105 * t83 + t146 * t81;
t26 = t93 * t87;
t27 = t163 * t87;
t109 = -t27 * pkin(4) - t26 * qJ(5) + t137;
t28 = t93 * t90;
t29 = t163 * t90;
t108 = -t29 * pkin(4) - t28 * qJ(5) + t136;
t107 = -t130 - t159;
t106 = t147 * t84 + t50 * t82;
t104 = t149 * t82 - t142;
t40 = (-t150 * t83 + t155) * t86;
t102 = g(1) * t22 + g(2) * t20 - g(3) * t40;
t41 = (t150 * t81 + t151) * t86;
t101 = g(1) * t23 + g(2) * t21 - g(3) * t41;
t31 = -t127 + (t151 * t84 + t81 * t88) * t89;
t2 = g(1) * t28 + g(2) * t26 - g(3) * t31;
t32 = -t142 * t155 + t143 * t83 + t128;
t3 = g(1) * t29 + g(2) * t27 - g(3) * t32;
t100 = -t41 * pkin(4) - t40 * qJ(5) + t117;
t99 = -t107 + t162;
t98 = t32 * pkin(4) + t31 * qJ(5) + t114;
t97 = t14 * pkin(4) + t13 * qJ(5) + t110;
t95 = t53 * pkin(3) - t133 * t52 + t90 * t69 + t118;
t94 = ((-pkin(9) - t132) * t86 + t120) * t161;
t92 = t16 * pkin(4) + t15 * qJ(5) + t95;
t91 = (pkin(3) * t88 + pkin(2) + t124) * t164;
t56 = t111 * t86;
t49 = t153 * t85 + t148;
t44 = g(3) * t86 + t68 * t89;
t43 = t104 * t90;
t42 = t104 * t87;
t34 = t145 * t84 - t52 * t82;
t17 = t99 * t82;
t10 = g(1) * t106 - g(2) * t34;
t9 = g(1) * t43 + g(2) * t42 - g(3) * t49;
t8 = -g(1) * t34 - g(2) * t106 - g(3) * t104;
t1 = -g(1) * t15 + t107 * t81 + (-g(2) * (-t147 * t82 + t50 * t84) - g(3) * t105) * t83;
t4 = [0, 0, 0, 0, 0, 0, t111, t68, 0, 0, 0, 0, 0, 0, 0, 0, t111 * t89, -t56, -t68, -g(1) * (-t87 * pkin(1) + t79) - g(2) * t134, 0, 0, 0, 0, 0, 0, g(1) * t51 - g(2) * t53, -g(1) * t50 - g(2) * t52, t56, -g(1) * t79 - g(2) * t118 - (t120 - t76) * t161, 0, 0, 0, 0, 0, 0, -t112, t113, t10, -g(1) * t110 - g(2) * t95 - t94, 0, 0, 0, 0, 0, 0, t10, t112, -t113, -g(1) * t97 - g(2) * t92 - t94, 0, 0, 0, 0, 0, 0, t10, -t113, -t112, -g(1) * (-pkin(5) * t106 + t14 * qJ(6) + t97) - g(2) * (t34 * pkin(5) + t16 * qJ(6) + t92) - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t44, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t88, -t96 * t85, -t44, -g(1) * (-pkin(2) * t145 + t73) - g(2) * (-pkin(2) * t147 + t70) - g(3) * t135, 0, 0, 0, 0, 0, 0, t3, -t2, t9, -g(1) * t136 - g(2) * t137 - g(3) * t114 + t91, 0, 0, 0, 0, 0, 0, t9, -t3, t2, -g(1) * t108 - g(2) * t109 - g(3) * t98 + t91, 0, 0, 0, 0, 0, 0, t9, t2, t3, -g(1) * (-t43 * pkin(5) - t29 * qJ(6) + t108) - g(2) * (-t42 * pkin(5) - t27 * qJ(6) + t109) - g(3) * (t49 * pkin(5) + t32 * qJ(6) + t98) + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t52 + g(2) * t50 + g(3) * t149, t99, 0, 0, 0, 0, 0, 0, 0, 0, -t101, t102, -t17, -g(1) * (t47 + t121) - g(2) * (-t45 + t122) - g(3) * t117, 0, 0, 0, 0, 0, 0, -t17, t101, -t102, -g(1) * (t115 + t121) - g(2) * (t116 + t122) - g(3) * t100, 0, 0, 0, 0, 0, 0, -t17, -t102, -t101, -g(1) * (t23 * qJ(6) + t115) - g(2) * (t21 * qJ(6) + t116) - g(3) * (-t41 * qJ(6) + t100) + (-pkin(5) * t130 + (-t159 - t162) * (pkin(5) + qJ(4))) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t14 - g(3) * t163;];
taug_reg  = t4;
