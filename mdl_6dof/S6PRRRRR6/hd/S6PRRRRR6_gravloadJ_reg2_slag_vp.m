% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:41:51
% EndTime: 2019-05-05 12:41:56
% DurationCPUTime: 1.58s
% Computational Cost: add. (2357->237), mult. (6881->391), div. (0->0), fcn. (9003->18), ass. (0->143)
t101 = cos(qJ(3));
t102 = cos(qJ(2));
t168 = sin(pkin(14));
t173 = cos(pkin(6));
t146 = t173 * t168;
t170 = cos(pkin(14));
t181 = sin(qJ(2));
t118 = t102 * t146 + t170 * t181;
t169 = sin(pkin(6));
t143 = t169 * t168;
t172 = cos(pkin(7));
t94 = sin(pkin(7));
t112 = -t118 * t172 + t94 * t143;
t86 = t170 * t102 - t181 * t146;
t98 = sin(qJ(3));
t105 = -t112 * t101 + t86 * t98;
t113 = t118 * t94 + t172 * t143;
t171 = cos(pkin(8));
t93 = sin(pkin(8));
t190 = t105 * t171 - t113 * t93;
t147 = t173 * t170;
t119 = -t102 * t147 + t168 * t181;
t144 = t170 * t169;
t110 = -t119 * t172 - t94 * t144;
t85 = t168 * t102 + t181 * t147;
t106 = -t110 * t101 + t85 * t98;
t111 = t119 * t94 - t172 * t144;
t189 = t106 * t171 - t111 * t93;
t145 = t172 * t169;
t121 = t102 * t145 + t173 * t94;
t150 = t169 * t181;
t115 = -t121 * t101 + t98 * t150;
t160 = t102 * t169;
t120 = -t94 * t160 + t173 * t172;
t188 = t115 * t171 - t120 * t93;
t164 = t94 * t171;
t161 = t101 * t172;
t66 = t119 * t98 - t85 * t161;
t48 = t85 * t164 - t66 * t93;
t68 = t118 * t98 - t86 * t161;
t49 = t86 * t164 - t68 * t93;
t187 = -g(1) * t49 - g(2) * t48;
t186 = pkin(10) * t94;
t185 = pkin(11) * t93;
t182 = cos(qJ(4));
t126 = t181 * t145;
t81 = -t101 * t126 - t98 * t160;
t180 = t81 * t93;
t179 = t93 * t94;
t96 = sin(qJ(5));
t178 = t93 * t96;
t142 = t94 * t150;
t177 = pkin(2) * t160 + pkin(10) * t142;
t100 = cos(qJ(5));
t176 = t100 * t93;
t95 = sin(qJ(6));
t175 = t100 * t95;
t99 = cos(qJ(6));
t174 = t100 * t99;
t60 = t85 * t101 + t110 * t98;
t97 = sin(qJ(4));
t18 = t189 * t182 + t60 * t97;
t19 = t60 * t182 - t189 * t97;
t167 = -t18 * pkin(4) + t19 * pkin(12);
t61 = t86 * t101 + t112 * t98;
t20 = t190 * t182 + t61 * t97;
t21 = t61 * t182 - t190 * t97;
t166 = -t20 * pkin(4) + t21 * pkin(12);
t77 = t101 * t150 + t121 * t98;
t39 = t188 * t182 + t77 * t97;
t40 = t77 * t182 - t188 * t97;
t165 = -t39 * pkin(4) + t40 * pkin(12);
t163 = t97 * t171;
t162 = t98 * t172;
t123 = t171 * t142;
t82 = t101 * t160 - t98 * t126;
t159 = t82 * pkin(3) + pkin(11) * t123 + t177;
t158 = t182 * t179;
t157 = -t119 * pkin(2) + t85 * t186;
t156 = -t118 * pkin(2) + t86 * t186;
t155 = -t106 * pkin(3) + t60 * t185;
t154 = -t105 * pkin(3) + t61 * t185;
t153 = -t115 * pkin(3) + t77 * t185;
t152 = -pkin(5) * t100 - pkin(13) * t96;
t151 = t171 * t182;
t67 = -t119 * t101 - t85 * t162;
t149 = t67 * pkin(3) + t157;
t69 = -t118 * t101 - t86 * t162;
t148 = t69 * pkin(3) + t156;
t59 = t115 * t93 + t120 * t171;
t16 = t59 * t100 - t40 * t96;
t41 = t106 * t93 + t111 * t171;
t2 = t41 * t100 - t19 * t96;
t42 = t105 * t93 + t113 * t171;
t4 = t42 * t100 - t21 * t96;
t140 = g(1) * t4 + g(2) * t2 + g(3) * t16;
t17 = t40 * t100 + t59 * t96;
t3 = t19 * t100 + t41 * t96;
t5 = t21 * t100 + t42 * t96;
t139 = g(1) * t5 + g(2) * t3 + g(3) * t17;
t137 = t93 * t142;
t52 = t82 * t182 + (t171 * t81 + t137) * t97;
t71 = t123 - t180;
t34 = -t71 * t100 + t52 * t96;
t27 = t67 * t182 + (t171 * t66 + t85 * t179) * t97;
t6 = -t48 * t100 + t27 * t96;
t29 = t69 * t182 + (t171 * t68 + t86 * t179) * t97;
t8 = -t49 * t100 + t29 * t96;
t138 = g(1) * t8 + g(2) * t6 + g(3) * t34;
t31 = -t106 * t182 - t60 * t163;
t10 = -t60 * t176 + t31 * t96;
t33 = -t105 * t182 - t61 * t163;
t12 = -t61 * t176 + t33 * t96;
t47 = -t115 * t182 - t77 * t163;
t36 = -t77 * t176 + t47 * t96;
t136 = g(1) * t12 + g(2) * t10 + g(3) * t36;
t135 = g(1) * t20 + g(2) * t18 + g(3) * t39;
t134 = g(1) * t21 + g(2) * t19 + g(3) * t40;
t26 = -t66 * t151 - t85 * t158 + t67 * t97;
t28 = -t68 * t151 - t86 * t158 + t69 * t97;
t51 = -t182 * t137 - t81 * t151 + t82 * t97;
t133 = g(1) * t28 + g(2) * t26 + g(3) * t51;
t30 = -t106 * t97 + t60 * t151;
t32 = -t105 * t97 + t61 * t151;
t46 = -t115 * t97 + t77 * t151;
t132 = g(1) * t32 + g(2) * t30 + g(3) * t46;
t131 = g(1) * t61 + g(2) * t60 + g(3) * t77;
t130 = t52 * pkin(4) + t51 * pkin(12) + t159;
t129 = t31 * pkin(4) + t30 * pkin(12) + t155;
t128 = t33 * pkin(4) + t32 * pkin(12) + t154;
t127 = t47 * pkin(4) + t46 * pkin(12) + t153;
t125 = t27 * pkin(4) + t26 * pkin(12) + t149;
t124 = t29 * pkin(4) + t28 * pkin(12) + t148;
t122 = -g(1) * t86 - g(2) * t85 - g(3) * t150;
t109 = (g(3) * t180 + t187) * pkin(11);
t37 = t47 * t100 + t77 * t178;
t35 = t52 * t100 + t71 * t96;
t13 = t33 * t100 + t61 * t178;
t11 = t31 * t100 + t60 * t178;
t9 = t29 * t100 + t49 * t96;
t7 = t27 * t100 + t48 * t96;
t1 = t135 * t96;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t118 + g(2) * t119 - g(3) * t160, -t122, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t69 - g(2) * t67 - g(3) * t82, -g(1) * t68 - g(2) * t66 - g(3) * t81, t122 * t94, -g(1) * t156 - g(2) * t157 - g(3) * t177, 0, 0, 0, 0, 0, 0, -g(1) * t29 - g(2) * t27 - g(3) * t52, t133, -g(3) * t71 + t187, -g(1) * t148 - g(2) * t149 - g(3) * t159 + t109, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t35, t138, -t133, -g(1) * t124 - g(2) * t125 - g(3) * t130 + t109, 0, 0, 0, 0, 0, 0, -g(1) * (t28 * t95 + t9 * t99) - g(2) * (t26 * t95 + t7 * t99) - g(3) * (t35 * t99 + t51 * t95) -g(1) * (t28 * t99 - t9 * t95) - g(2) * (t26 * t99 - t7 * t95) - g(3) * (-t35 * t95 + t51 * t99) -t138, -g(1) * (t9 * pkin(5) + t8 * pkin(13) + t124) - g(2) * (t7 * pkin(5) + t6 * pkin(13) + t125) - g(3) * (t35 * pkin(5) + t34 * pkin(13) + t130) + t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122 * t98 + (-g(1) * t112 - g(2) * t110 - g(3) * t121) * t101, t131, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t33 - g(2) * t31 - g(3) * t47, t132, -t131 * t93, -g(1) * t154 - g(2) * t155 - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t37, t136, -t132, -g(1) * t128 - g(2) * t129 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t99 + t32 * t95) - g(2) * (t11 * t99 + t30 * t95) - g(3) * (t37 * t99 + t46 * t95) -g(1) * (-t13 * t95 + t32 * t99) - g(2) * (-t11 * t95 + t30 * t99) - g(3) * (-t37 * t95 + t46 * t99) -t136, -g(1) * (t13 * pkin(5) + t12 * pkin(13) + t128) - g(2) * (t11 * pkin(5) + t10 * pkin(13) + t129) - g(3) * (t37 * pkin(5) + t36 * pkin(13) + t127); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t134, 0, 0, 0, 0, 0, 0, 0, 0, t135 * t100, -t1, -t134, -g(1) * t166 - g(2) * t167 - g(3) * t165, 0, 0, 0, 0, 0, 0, -g(1) * (-t20 * t174 + t21 * t95) - g(2) * (-t18 * t174 + t19 * t95) - g(3) * (-t39 * t174 + t40 * t95) -g(1) * (t20 * t175 + t21 * t99) - g(2) * (t18 * t175 + t19 * t99) - g(3) * (t39 * t175 + t40 * t99) t1, -g(1) * (t152 * t20 + t166) - g(2) * (t152 * t18 + t167) - g(3) * (t152 * t39 + t165); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t139, 0, 0, 0, 0, 0, 0, 0, 0, -t140 * t99, t140 * t95, -t139, -g(1) * (t4 * pkin(5) + t5 * pkin(13)) - g(2) * (t2 * pkin(5) + t3 * pkin(13)) - g(3) * (t16 * pkin(5) + t17 * pkin(13)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t20 * t99 - t5 * t95) - g(2) * (t18 * t99 - t3 * t95) - g(3) * (-t17 * t95 + t39 * t99) -g(1) * (-t20 * t95 - t5 * t99) - g(2) * (-t18 * t95 - t3 * t99) - g(3) * (-t17 * t99 - t39 * t95) 0, 0;];
taug_reg  = t14;
