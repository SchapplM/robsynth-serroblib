% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR9
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t101 = cos(qJ(6));
t103 = cos(qJ(3));
t151 = cos(pkin(13));
t92 = sin(pkin(13));
t98 = sin(qJ(3));
t121 = t103 * t92 + t98 * t151;
t105 = cos(qJ(1));
t94 = sin(pkin(6));
t153 = t105 * t94;
t77 = -t103 * t151 + t92 * t98;
t93 = sin(pkin(7));
t66 = t77 * t93;
t95 = cos(pkin(7));
t68 = t77 * t95;
t100 = sin(qJ(1));
t104 = cos(qJ(2));
t152 = cos(pkin(6));
t140 = t105 * t152;
t99 = sin(qJ(2));
t73 = t100 * t99 - t104 * t140;
t74 = t100 * t104 + t99 * t140;
t23 = t121 * t74 - t66 * t153 - t73 * t68;
t102 = cos(qJ(5));
t67 = t121 * t93;
t69 = t121 * t95;
t110 = t67 * t153 + t73 * t69 + t74 * t77;
t54 = t95 * t153 - t73 * t93;
t97 = sin(qJ(5));
t7 = t102 * t110 + t54 * t97;
t96 = sin(qJ(6));
t180 = t101 * t23 + t7 * t96;
t179 = t101 * t7 - t23 * t96;
t142 = t152 * t93;
t155 = t104 * t94;
t176 = t95 * t155 + t142;
t175 = -t54 * t102 + t110 * t97;
t160 = t100 * t94;
t141 = t100 * t152;
t75 = -t104 * t141 - t105 * t99;
t128 = t95 * t160 - t75 * t93;
t39 = -g(1) * t54 - g(2) * t128;
t32 = (t104 * t69 - t77 * t99) * t94 + t152 * t67;
t174 = pkin(10) * t93;
t172 = g(3) * t94;
t171 = pkin(3) * t103;
t170 = t74 * t98;
t76 = t104 * t105 - t99 * t141;
t169 = t76 * t98;
t168 = t93 * t97;
t167 = t94 * t99;
t166 = t95 * t98;
t165 = t98 * t99;
t164 = pkin(10) + qJ(4);
t71 = pkin(3) * t166 - t164 * t93;
t90 = pkin(2) + t171;
t163 = -t74 * t71 - t73 * t90;
t162 = -t76 * t71 + t75 * t90;
t161 = t105 * pkin(1) + pkin(9) * t160;
t159 = t102 * t93;
t158 = t102 * t96;
t157 = t103 * t95;
t156 = t103 * t99;
t154 = t104 * t98;
t150 = t101 * t102;
t149 = pkin(3) * t157;
t148 = t93 * t167;
t147 = t93 * t160;
t145 = t93 * t153;
t144 = -t100 * pkin(1) + pkin(9) * t153;
t27 = t67 * t160 + t69 * t75 - t76 * t77;
t8 = -t128 * t102 + t27 * t97;
t139 = g(1) * t175 + g(2) * t8;
t138 = t90 * t155 - t71 * t167;
t70 = pkin(3) * t93 * t98 + t164 * t95;
t137 = t70 * t160 + t75 * t71 + t76 * t90 + t161;
t26 = -t121 * t76 - t66 * t160 - t68 * t75;
t136 = g(1) * t23 + g(2) * t26;
t135 = pkin(5) * t102 + pkin(12) * t97;
t134 = g(1) * t105 + g(2) * t100;
t132 = -pkin(3) * t169 + t147 * t171 + t75 * t149;
t35 = t121 * t73 + t68 * t74;
t36 = -t69 * t74 + t73 * t77;
t131 = t36 * pkin(4) - pkin(11) * t35 + t163;
t37 = -t121 * t75 + t68 * t76;
t38 = -t69 * t76 - t75 * t77;
t130 = t38 * pkin(4) - pkin(11) * t37 + t162;
t129 = g(2) * t74 + g(3) * t167;
t127 = t75 * t95 + t147;
t126 = t73 * t95 + t145;
t125 = -pkin(3) * t94 * t165 + t176 * t171;
t124 = t70 * t153 + t73 * t71 - t74 * t90 + t144;
t72 = t152 * t95 - t93 * t155;
t14 = t102 * t72 - t32 * t97;
t123 = g(1) * t8 - g(2) * t175 - g(3) * t14;
t15 = t102 * t32 + t72 * t97;
t9 = t27 * t102 + t128 * t97;
t122 = g(1) * t9 - g(2) * t7 + g(3) * t15;
t16 = -t74 * t159 + t36 * t97;
t18 = -t76 * t159 + t38 * t97;
t47 = (-t104 * t77 - t69 * t99) * t94;
t40 = -t102 * t148 + t47 * t97;
t120 = g(1) * t18 + g(2) * t16 + g(3) * t40;
t119 = -g(1) * t27 + g(2) * t110 - g(3) * t32;
t31 = -t152 * t66 + (-t104 * t68 - t121 * t99) * t94;
t118 = g(1) * t26 - g(2) * t23 + g(3) * t31;
t46 = -t121 * t155 + t68 * t167;
t117 = g(1) * t37 + g(2) * t35 + g(3) * t46;
t116 = -t103 * t74 + t98 * t145 + t73 * t166;
t114 = t47 * pkin(4) - pkin(11) * t46 + t138;
t113 = g(1) * t76 + t129;
t112 = t27 * pkin(4) - pkin(11) * t26 + t137;
t111 = t26 * pkin(4) + pkin(11) * t27 + t132;
t109 = t31 * pkin(4) + pkin(11) * t32 + t125;
t108 = pkin(4) * t110 - pkin(11) * t23 + t124;
t107 = -t73 * t149 + (-t103 * t145 - t170) * pkin(3);
t106 = -pkin(4) * t23 - pkin(11) * t110 + t107;
t44 = t113 * t93;
t43 = t103 * t76 + t127 * t98;
t42 = t127 * t103 - t169;
t41 = t102 * t47 + t97 * t148;
t28 = -g(1) * t128 + g(2) * t54 - g(3) * t72;
t19 = t102 * t38 + t76 * t168;
t17 = t102 * t36 + t74 * t168;
t3 = t101 * t9 - t26 * t96;
t2 = -t101 * t26 - t9 * t96;
t1 = t118 * t97;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t100 - g(2) * t105, t134, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t74 - g(2) * t76, -g(1) * t73 - g(2) * t75, -t134 * t94, -g(1) * t144 - g(2) * t161, 0, 0, 0, 0, 0, 0, -g(1) * t116 - g(2) * t43, -g(1) * (t126 * t103 + t170) - g(2) * t42, t39, -g(1) * (-t74 * pkin(2) + t144) - g(2) * (pkin(2) * t76 + t161) + t39 * pkin(10), 0, 0, 0, 0, 0, 0, -g(1) * t110 - g(2) * t27, -t136, t39, -g(1) * t124 - g(2) * t137, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t139, t136, -g(1) * t108 - g(2) * t112, 0, 0, 0, 0, 0, 0, -g(1) * t179 - g(2) * t3, g(1) * t180 - g(2) * t2, -t139, -g(1) * (pkin(5) * t7 + pkin(12) * t175 + t108) - g(2) * (pkin(5) * t9 + pkin(12) * t8 + t112); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t75 + g(2) * t73 - g(3) * t155, t113, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t75 - t76 * t166) - g(2) * (-t103 * t73 - t74 * t166) - (t103 * t104 - t95 * t165) * t172, -g(1) * (-t76 * t157 - t75 * t98) - g(2) * (-t74 * t157 + t73 * t98) - (-t95 * t156 - t154) * t172, -t44, -g(1) * (pkin(2) * t75 + t76 * t174) - g(2) * (-pkin(2) * t73 + t74 * t174) - (pkin(2) * t104 + t99 * t174) * t172, 0, 0, 0, 0, 0, 0, -g(1) * t38 - g(2) * t36 - g(3) * t47, -t117, -t44, -g(1) * t162 - g(2) * t163 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t17 - g(3) * t41, t120, t117, -g(1) * t130 - g(2) * t131 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t19 - t37 * t96) - g(2) * (t101 * t17 - t35 * t96) - g(3) * (t101 * t41 - t46 * t96) -g(1) * (-t101 * t37 - t19 * t96) - g(2) * (-t101 * t35 - t17 * t96) - g(3) * (-t101 * t46 - t41 * t96) -t120, -g(1) * (pkin(5) * t19 + pkin(12) * t18 + t130) - g(2) * (pkin(5) * t17 + pkin(12) * t16 + t131) - g(3) * (pkin(5) * t41 + pkin(12) * t40 + t114); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t42 + t129 * t98 + (g(2) * t126 - g(3) * t176) * t103, g(1) * t43 - g(2) * t116 - g(3) * (-t98 * t142 + (-t95 * t154 - t156) * t94) 0, 0, 0, 0, 0, 0, 0, 0, -t118, -t119, 0, -g(1) * t132 - g(2) * t107 - g(3) * t125, 0, 0, 0, 0, 0, 0, -t118 * t102, t1, t119, -g(1) * t111 - g(2) * t106 - g(3) * t109, 0, 0, 0, 0, 0, 0, -g(1) * (t150 * t26 + t27 * t96) - g(2) * (-t110 * t96 - t150 * t23) - g(3) * (t150 * t31 + t32 * t96) -g(1) * (t101 * t27 - t158 * t26) - g(2) * (-t101 * t110 + t158 * t23) - g(3) * (t101 * t32 - t158 * t31) -t1, -g(1) * (t135 * t26 + t111) - g(2) * (-t135 * t23 + t106) - g(3) * (t135 * t31 + t109); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t122, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t101, -t123 * t96, -t122, -g(1) * (-pkin(5) * t8 + pkin(12) * t9) - g(2) * (pkin(5) * t175 - pkin(12) * t7) - g(3) * (pkin(5) * t14 + pkin(12) * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t180 - g(3) * (-t101 * t31 - t15 * t96) g(1) * t3 - g(2) * t179 - g(3) * (-t101 * t15 + t31 * t96) 0, 0;];
taug_reg  = t4;
