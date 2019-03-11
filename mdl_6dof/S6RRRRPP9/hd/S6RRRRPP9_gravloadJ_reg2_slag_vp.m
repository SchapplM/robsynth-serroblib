% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t156 = pkin(4) * t89 + qJ(5) * t86;
t147 = cos(qJ(1));
t85 = sin(pkin(6));
t129 = t85 * t147;
t132 = cos(pkin(6));
t112 = t132 * t147;
t146 = sin(qJ(1));
t88 = sin(qJ(2));
t91 = cos(qJ(2));
t68 = t88 * t112 + t146 * t91;
t87 = sin(qJ(3));
t90 = cos(qJ(3));
t41 = -t87 * t129 + t68 * t90;
t67 = -t91 * t112 + t146 * t88;
t16 = t41 * t86 - t67 * t89;
t17 = t41 * t89 + t67 * t86;
t122 = -t90 * t129 - t68 * t87;
t155 = t156 * t122;
t128 = t85 * t146;
t111 = t132 * t146;
t70 = -t88 * t111 + t147 * t91;
t44 = -t90 * t128 + t70 * t87;
t154 = t156 * t44;
t141 = t85 * t88;
t65 = t132 * t90 - t87 * t141;
t153 = t156 * t65;
t152 = pkin(5) + pkin(10);
t151 = pkin(3) * t90;
t149 = t122 * pkin(10);
t148 = t44 * pkin(10);
t144 = t67 * t87;
t69 = t91 * t111 + t147 * t88;
t142 = t69 * t87;
t140 = t85 * t91;
t139 = t86 * t90;
t138 = t89 * t90;
t137 = t90 * t91;
t136 = pkin(2) * t140 + pkin(9) * t141;
t135 = t147 * pkin(1) + pkin(8) * t128;
t133 = qJ(6) * t89;
t131 = t87 * t140;
t130 = t86 * t140;
t127 = -t67 * pkin(2) + t68 * pkin(9);
t126 = -t69 * pkin(2) + t70 * pkin(9);
t34 = t122 * pkin(3);
t125 = t41 * pkin(10) + t34;
t36 = t44 * pkin(3);
t45 = t87 * t128 + t70 * t90;
t124 = t45 * pkin(10) - t36;
t60 = t65 * pkin(3);
t66 = t132 * t87 + t90 * t141;
t123 = t66 * pkin(10) + t60;
t121 = -t16 * pkin(4) + t17 * qJ(5);
t20 = t45 * t86 - t69 * t89;
t21 = t45 * t89 + t69 * t86;
t120 = -t20 * pkin(4) + t21 * qJ(5);
t38 = t89 * t140 + t66 * t86;
t39 = t66 * t89 - t130;
t119 = -t38 * pkin(4) + t39 * qJ(5);
t118 = t85 * pkin(3) * t137 + pkin(10) * t131 + t136;
t117 = -t146 * pkin(1) + pkin(8) * t129;
t116 = -g(1) * t16 + g(2) * t20;
t115 = -g(1) * t17 + g(2) * t21;
t114 = g(1) * t122 + g(2) * t44;
t113 = g(1) * t67 - g(2) * t69;
t110 = -pkin(10) * t144 - t67 * t151 + t127;
t109 = -pkin(10) * t142 - t69 * t151 + t126;
t108 = t70 * pkin(2) + t69 * pkin(9) + t135;
t107 = t45 * pkin(3) + t108;
t2 = g(1) * t20 + g(2) * t16 + g(3) * t38;
t106 = g(1) * t21 + g(2) * t17 + g(3) * t39;
t25 = -t67 * t139 - t68 * t89;
t27 = -t69 * t139 - t70 * t89;
t47 = t90 * t130 - t89 * t141;
t105 = g(1) * t27 + g(2) * t25 + g(3) * t47;
t26 = -t67 * t138 + t68 * t86;
t28 = -t69 * t138 + t70 * t86;
t48 = (t89 * t137 + t86 * t88) * t85;
t104 = g(1) * t28 + g(2) * t26 + g(3) * t48;
t103 = g(1) * t44 - g(2) * t122 - g(3) * t65;
t102 = g(1) * t45 + g(2) * t41 + g(3) * t66;
t101 = t48 * pkin(4) + t47 * qJ(5) + t118;
t100 = g(1) * t147 + g(2) * t146;
t99 = -t68 * pkin(2) - t67 * pkin(9) + t117;
t98 = -g(1) * t69 - g(2) * t67 + g(3) * t140;
t97 = g(1) * t70 + g(2) * t68 + g(3) * t141;
t96 = t26 * pkin(4) + t25 * qJ(5) + t110;
t95 = t28 * pkin(4) + t27 * qJ(5) + t109;
t94 = -pkin(3) * t41 + t99;
t93 = t21 * pkin(4) + t20 * qJ(5) + t107;
t92 = -pkin(4) * t17 - qJ(5) * t16 + t94;
t22 = t98 * t87;
t9 = t103 * t89;
t8 = t103 * t86;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t146 - g(2) * t147, t100, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t70, -t113, -t100 * t85, -g(1) * t117 - g(2) * t135, 0, 0, 0, 0, 0, 0, g(1) * t41 - g(2) * t45, t114, t113, -g(1) * t99 - g(2) * t108, 0, 0, 0, 0, 0, 0, -t115, t116, -t114, -g(1) * (t94 + t149) - g(2) * (t107 + t148) 0, 0, 0, 0, 0, 0, -t114, t115, -t116, -g(1) * (t92 + t149) - g(2) * (t93 + t148) 0, 0, 0, 0, 0, 0, -t114, -t116, -t115, -g(1) * (-qJ(6) * t17 + t122 * t152 + t92) - g(2) * (t21 * qJ(6) + t152 * t44 + t93); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t97, 0, 0, 0, 0, 0, 0, 0, 0, -t98 * t90, t22, -t97, -g(1) * t126 - g(2) * t127 - g(3) * t136, 0, 0, 0, 0, 0, 0, -t104, t105, -t22, -g(1) * t109 - g(2) * t110 - g(3) * t118, 0, 0, 0, 0, 0, 0, -t22, t104, -t105, -g(1) * t95 - g(2) * t96 - g(3) * t101, 0, 0, 0, 0, 0, 0, -t22, -t105, -t104, -g(1) * (-pkin(5) * t142 + t28 * qJ(6) + t95) - g(2) * (-pkin(5) * t144 + t26 * qJ(6) + t96) - g(3) * (pkin(5) * t131 + t48 * qJ(6) + t101); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t102, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, -t102, -g(1) * t124 - g(2) * t125 - g(3) * t123, 0, 0, 0, 0, 0, 0, -t102, -t9, t8, -g(1) * (t124 - t154) - g(2) * (t125 + t155) - g(3) * (t123 + t153) 0, 0, 0, 0, 0, 0, -t102, t8, t9, -g(1) * (-t44 * t133 + t152 * t45 - t154 - t36) - g(2) * (t122 * t133 + t152 * t41 + t155 + t34) - g(3) * (t65 * t133 + t152 * t66 + t153 + t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t106, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t106, -g(1) * t120 - g(2) * t121 - g(3) * t119, 0, 0, 0, 0, 0, 0, 0, -t106, t2, -g(1) * (-t20 * qJ(6) + t120) - g(2) * (-t16 * qJ(6) + t121) - g(3) * (-t38 * qJ(6) + t119); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106;];
taug_reg  = t1;
