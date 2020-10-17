% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:47:46
% EndTime: 2019-05-07 15:47:53
% DurationCPUTime: 1.98s
% Computational Cost: add. (1315->212), mult. (3217->343), div. (0->0), fcn. (4099->16), ass. (0->110)
t143 = cos(pkin(7));
t161 = cos(qJ(3));
t120 = t143 * t161;
t92 = sin(pkin(6));
t98 = cos(qJ(1));
t150 = t92 * t98;
t91 = sin(pkin(7));
t142 = t91 * t150;
t144 = cos(pkin(6));
t162 = cos(qJ(2));
t122 = t144 * t162;
t159 = sin(qJ(2));
t160 = sin(qJ(1));
t68 = -t98 * t122 + t160 * t159;
t121 = t144 * t159;
t69 = t98 * t121 + t160 * t162;
t96 = sin(qJ(3));
t28 = t68 * t120 + t161 * t142 + t69 * t96;
t135 = t92 * t143;
t169 = t98 * t135 - t68 * t91;
t134 = t96 * t143;
t29 = -t68 * t134 - t96 * t142 + t69 * t161;
t89 = pkin(13) + qJ(5);
t86 = sin(t89);
t87 = cos(t89);
t9 = -t169 * t86 + t29 * t87;
t95 = sin(qJ(6));
t97 = cos(qJ(6));
t174 = -t28 * t97 + t9 * t95;
t173 = t28 * t95 + t9 * t97;
t107 = t160 * t122 + t98 * t159;
t170 = t107 * t91 + t160 * t135;
t168 = -t169 * t87 - t29 * t86;
t137 = t92 * t160;
t167 = t107 * t143 - t91 * t137;
t70 = -t160 * t121 + t98 * t162;
t166 = -g(1) * t70 - g(2) * t69;
t165 = pkin(10) * t91;
t90 = sin(pkin(13));
t158 = t169 * t90;
t156 = t86 * t91;
t155 = t87 * t91;
t154 = t87 * t95;
t153 = t87 * t97;
t152 = t90 * t91;
t93 = cos(pkin(13));
t151 = t91 * t93;
t84 = pkin(4) * t93 + pkin(3);
t94 = -pkin(11) - qJ(4);
t149 = -t28 * t84 - t29 * t94;
t32 = t167 * t161 + t70 * t96;
t33 = t70 * t161 - t167 * t96;
t148 = -t32 * t84 - t33 * t94;
t133 = t144 * t91;
t136 = t92 * t159;
t138 = t92 * t162;
t50 = -t120 * t138 - t161 * t133 + t96 * t136;
t51 = t96 * t133 + (t162 * t134 + t161 * t159) * t92;
t147 = -t50 * t84 - t51 * t94;
t132 = t91 * t136;
t146 = pkin(2) * t138 + pkin(10) * t132;
t145 = t98 * pkin(1) + pkin(9) * t137;
t40 = t69 * t120 - t68 * t96;
t41 = -t69 * t134 - t68 * t161;
t63 = t68 * pkin(2);
t141 = -t40 * t94 + t41 * t84 - t63;
t42 = -t107 * t96 + t70 * t120;
t43 = -t107 * t161 - t70 * t134;
t65 = t107 * pkin(2);
t140 = -t42 * t94 + t43 * t84 - t65;
t130 = t69 * t165 - t63;
t129 = t70 * t165 - t65;
t127 = -t160 * pkin(1) + pkin(9) * t150;
t118 = t90 * t132;
t119 = t143 * t159;
t61 = (t161 * t119 + t162 * t96) * t92;
t62 = (-t96 * t119 + t161 * t162) * t92;
t126 = pkin(4) * t118 - t61 * t94 + t62 * t84 + t146;
t125 = -pkin(5) * t87 - pkin(12) * t86;
t12 = -t170 * t87 + t33 * t86;
t124 = g(1) * t168 + g(2) * t12;
t123 = -g(1) * t28 + g(2) * t32;
t116 = -g(1) * t98 - g(2) * t160;
t67 = -t91 * t138 + t144 * t143;
t22 = -t51 * t86 + t67 * t87;
t115 = g(1) * t12 - g(2) * t168 - g(3) * t22;
t13 = t170 * t86 + t33 * t87;
t23 = t51 * t87 + t67 * t86;
t114 = g(1) * t13 + g(2) * t9 + g(3) * t23;
t14 = -t69 * t155 + t41 * t86;
t16 = -t70 * t155 + t43 * t86;
t38 = -t87 * t132 + t62 * t86;
t113 = g(1) * t16 + g(2) * t14 + g(3) * t38;
t112 = g(1) * t32 + g(2) * t28 + g(3) * t50;
t111 = g(1) * t33 + g(2) * t29 + g(3) * t51;
t110 = g(1) * t42 + g(2) * t40 + g(3) * t61;
t109 = -g(3) * t136 + t166;
t108 = -t69 * pkin(2) + t169 * pkin(10) + t127;
t106 = pkin(4) * t158 + t28 * t94 - t29 * t84 + t108;
t104 = t166 * (pkin(4) * t90 + pkin(10)) * t91;
t102 = t70 * pkin(2) + t170 * pkin(10) + t145;
t99 = t170 * t90;
t100 = pkin(4) * t99 - t32 * t94 + t33 * t84 + t102;
t39 = t86 * t132 + t62 * t87;
t17 = t70 * t156 + t43 * t87;
t15 = t69 * t156 + t41 * t87;
t3 = t13 * t97 + t32 * t95;
t2 = -t13 * t95 + t32 * t97;
t1 = t112 * t86;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t160 - g(2) * t98, -t116, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t69 - g(2) * t70, -g(1) * t68 + g(2) * t107, t116 * t92, -g(1) * t127 - g(2) * t145, 0, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t33, t123, -g(1) * t169 - g(2) * t170, -g(1) * t108 - g(2) * t102, 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t93 + t158) - g(2) * (t33 * t93 + t99) -g(1) * (t169 * t93 + t29 * t90) - g(2) * (t170 * t93 - t33 * t90) -t123, -g(1) * (-pkin(3) * t29 - qJ(4) * t28 + t108) - g(2) * (t33 * pkin(3) + t32 * qJ(4) + t102) 0, 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t13, t124, -t123, -g(1) * t106 - g(2) * t100, 0, 0, 0, 0, 0, 0, g(1) * t173 - g(2) * t3, -g(1) * t174 - g(2) * t2, -t124, -g(1) * (-pkin(5) * t9 + pkin(12) * t168 + t106) - g(2) * (t13 * pkin(5) + t12 * pkin(12) + t100); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t107 + g(2) * t68 - g(3) * t138, -t109, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t43 - g(2) * t41 - g(3) * t62, t110, t109 * t91, -g(1) * t129 - g(2) * t130 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t152 + t43 * t93) - g(2) * (t69 * t152 + t41 * t93) - g(3) * (t62 * t93 + t118) -g(1) * (t70 * t151 - t43 * t90) - g(2) * (t69 * t151 - t41 * t90) - g(3) * (t93 * t132 - t62 * t90) -t110, -g(1) * (pkin(3) * t43 + qJ(4) * t42 + t129) - g(2) * (pkin(3) * t41 + qJ(4) * t40 + t130) - g(3) * (pkin(3) * t62 + qJ(4) * t61 + t146) 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15 - g(3) * t39, t113, -t110, -g(1) * t140 - g(2) * t141 - g(3) * t126 + t104, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t97 + t42 * t95) - g(2) * (t15 * t97 + t40 * t95) - g(3) * (t39 * t97 + t61 * t95) -g(1) * (-t17 * t95 + t42 * t97) - g(2) * (-t15 * t95 + t40 * t97) - g(3) * (-t39 * t95 + t61 * t97) -t113, -g(1) * (pkin(5) * t17 + pkin(12) * t16 + t140) - g(2) * (pkin(5) * t15 + pkin(12) * t14 + t141) - g(3) * (pkin(5) * t39 + pkin(12) * t38 + t126) + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t111, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t93, -t112 * t90, -t111, -g(1) * (-pkin(3) * t32 + qJ(4) * t33) - g(2) * (-pkin(3) * t28 + qJ(4) * t29) - g(3) * (-pkin(3) * t50 + qJ(4) * t51) 0, 0, 0, 0, 0, 0, t112 * t87, -t1, -t111, -g(1) * t148 - g(2) * t149 - g(3) * t147, 0, 0, 0, 0, 0, 0, -g(1) * (-t32 * t153 + t33 * t95) - g(2) * (-t28 * t153 + t29 * t95) - g(3) * (-t50 * t153 + t51 * t95) -g(1) * (t32 * t154 + t33 * t97) - g(2) * (t28 * t154 + t29 * t97) - g(3) * (t50 * t154 + t51 * t97) t1, -g(1) * (t125 * t32 + t148) - g(2) * (t125 * t28 + t149) - g(3) * (t125 * t50 + t147); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t114, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t97, -t115 * t95, -t114, -g(1) * (-pkin(5) * t12 + pkin(12) * t13) - g(2) * (pkin(5) * t168 + pkin(12) * t9) - g(3) * (pkin(5) * t22 + pkin(12) * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t174 - g(3) * (-t23 * t95 + t50 * t97) g(1) * t3 + g(2) * t173 - g(3) * (-t23 * t97 - t50 * t95) 0, 0;];
taug_reg  = t4;
