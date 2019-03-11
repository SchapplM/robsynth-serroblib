% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t102 = cos(qJ(6));
t101 = sin(qJ(3));
t151 = cos(pkin(7));
t176 = cos(qJ(3));
t125 = t151 * t176;
t104 = cos(qJ(1));
t97 = sin(pkin(6));
t153 = t104 * t97;
t96 = sin(pkin(7));
t148 = t96 * t153;
t152 = cos(pkin(6));
t177 = cos(qJ(2));
t127 = t152 * t177;
t174 = sin(qJ(2));
t175 = sin(qJ(1));
t74 = -t104 * t127 + t175 * t174;
t126 = t152 * t174;
t75 = t104 * t126 + t175 * t177;
t30 = t75 * t101 + t74 * t125 + t176 * t148;
t143 = t97 * t151;
t187 = t104 * t143 - t74 * t96;
t141 = t101 * t151;
t31 = -t101 * t148 - t74 * t141 + t75 * t176;
t95 = qJ(4) + pkin(13);
t92 = sin(t95);
t93 = cos(t95);
t9 = -t187 * t92 + t31 * t93;
t99 = sin(qJ(6));
t194 = -t102 * t30 + t9 * t99;
t193 = t102 * t9 + t30 * t99;
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t192 = t100 * t31 + t103 * t187;
t161 = t100 * t187;
t191 = -t103 * t31 + t161;
t111 = t104 * t174 + t175 * t127;
t58 = t111 * t96 + t175 * t143;
t145 = t97 * t174;
t182 = t177 * t143 + t152 * t96;
t55 = t182 * t101 + t176 * t145;
t147 = t97 * t177;
t73 = -t96 * t147 + t152 * t151;
t186 = -t100 * t55 + t103 * t73;
t146 = t97 * t175;
t183 = t111 * t151 - t96 * t146;
t76 = t104 * t177 - t175 * t126;
t35 = -t183 * t101 + t76 * t176;
t14 = -t100 * t35 + t103 * t58;
t184 = -t187 * t93 - t31 * t92;
t181 = -g(1) * t76 - g(2) * t75;
t180 = pkin(10) * t96;
t172 = t92 * t96;
t171 = t93 * t96;
t170 = t93 * t99;
t91 = pkin(4) * t103 + pkin(3);
t98 = -qJ(5) - pkin(11);
t169 = -t30 * t91 - t31 * t98;
t34 = t101 * t76 + t183 * t176;
t168 = -t34 * t91 - t35 * t98;
t54 = t101 * t145 - t182 * t176;
t167 = -t54 * t91 - t55 * t98;
t140 = t96 * t145;
t166 = pkin(2) * t147 + pkin(10) * t140;
t165 = t104 * pkin(1) + pkin(9) * t146;
t160 = t100 * t58;
t159 = t100 * t96;
t158 = t102 * t93;
t154 = t103 * t96;
t42 = -t101 * t74 + t75 * t125;
t43 = -t75 * t141 - t74 * t176;
t69 = t74 * pkin(2);
t150 = -t42 * t98 + t43 * t91 - t69;
t44 = -t111 * t101 + t76 * t125;
t45 = -t111 * t176 - t76 * t141;
t71 = t111 * pkin(2);
t149 = -t44 * t98 + t45 * t91 - t71;
t138 = t75 * t180 - t69;
t137 = t76 * t180 - t71;
t136 = t192 * pkin(4);
t135 = t14 * pkin(4);
t134 = t186 * pkin(4);
t132 = -t175 * pkin(1) + pkin(9) * t153;
t123 = t100 * t140;
t124 = t151 * t174;
t66 = (t177 * t101 + t176 * t124) * t97;
t67 = (-t101 * t124 + t176 * t177) * t97;
t131 = pkin(4) * t123 - t66 * t98 + t67 * t91 + t166;
t130 = -pkin(5) * t93 - pkin(12) * t92;
t12 = t35 * t92 - t58 * t93;
t129 = g(1) * t184 + g(2) * t12;
t128 = -g(1) * t30 + g(2) * t34;
t120 = -g(1) * t104 - g(2) * t175;
t24 = -t55 * t92 + t73 * t93;
t119 = g(1) * t12 - g(2) * t184 - g(3) * t24;
t13 = t35 * t93 + t58 * t92;
t25 = t55 * t93 + t73 * t92;
t118 = g(1) * t13 + g(2) * t9 + g(3) * t25;
t16 = -t75 * t171 + t43 * t92;
t18 = -t76 * t171 + t45 * t92;
t40 = -t93 * t140 + t67 * t92;
t117 = g(1) * t18 + g(2) * t16 + g(3) * t40;
t116 = g(1) * t34 + g(2) * t30 + g(3) * t54;
t115 = g(1) * t35 + g(2) * t31 + g(3) * t55;
t114 = g(1) * t44 + g(2) * t42 + g(3) * t66;
t113 = -g(3) * t145 + t181;
t112 = -t75 * pkin(2) + t187 * pkin(10) + t132;
t110 = pkin(4) * t161 + t30 * t98 - t31 * t91 + t112;
t108 = t181 * (pkin(4) * t100 + pkin(10)) * t96;
t106 = t76 * pkin(2) + t58 * pkin(10) + t165;
t105 = pkin(4) * t160 - t34 * t98 + t35 * t91 + t106;
t41 = t92 * t140 + t67 * t93;
t19 = t76 * t172 + t45 * t93;
t17 = t75 * t172 + t43 * t93;
t15 = t103 * t35 + t160;
t3 = t102 * t13 + t34 * t99;
t2 = t102 * t34 - t13 * t99;
t1 = t116 * t92;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t175 - g(2) * t104, -t120, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t75 - g(2) * t76, -g(1) * t74 + g(2) * t111, t120 * t97, -g(1) * t132 - g(2) * t165, 0, 0, 0, 0, 0, 0, g(1) * t31 - g(2) * t35, t128, -g(1) * t187 - g(2) * t58, -g(1) * t112 - g(2) * t106, 0, 0, 0, 0, 0, 0, -g(1) * t191 - g(2) * t15, -g(1) * t192 - g(2) * t14, -t128, -g(1) * (-pkin(3) * t31 - pkin(11) * t30 + t112) - g(2) * (t35 * pkin(3) + t34 * pkin(11) + t106) 0, 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t13, t129, -t128, -g(1) * t110 - g(2) * t105, 0, 0, 0, 0, 0, 0, g(1) * t193 - g(2) * t3, -g(1) * t194 - g(2) * t2, -t129, -g(1) * (-pkin(5) * t9 + pkin(12) * t184 + t110) - g(2) * (t13 * pkin(5) + t12 * pkin(12) + t105); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t111 + g(2) * t74 - g(3) * t147, -t113, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t43 - g(3) * t67, t114, t113 * t96, -g(1) * t137 - g(2) * t138 - g(3) * t166, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t45 + t76 * t159) - g(2) * (t103 * t43 + t75 * t159) - g(3) * (t67 * t103 + t123) -g(1) * (-t100 * t45 + t76 * t154) - g(2) * (-t100 * t43 + t75 * t154) - g(3) * (-t67 * t100 + t103 * t140) -t114, -g(1) * (pkin(3) * t45 + pkin(11) * t44 + t137) - g(2) * (pkin(3) * t43 + pkin(11) * t42 + t138) - g(3) * (pkin(3) * t67 + pkin(11) * t66 + t166) 0, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t17 - g(3) * t41, t117, -t114, -g(1) * t149 - g(2) * t150 - g(3) * t131 + t108, 0, 0, 0, 0, 0, 0, -g(1) * (t102 * t19 + t44 * t99) - g(2) * (t102 * t17 + t42 * t99) - g(3) * (t102 * t41 + t66 * t99) -g(1) * (t102 * t44 - t19 * t99) - g(2) * (t102 * t42 - t17 * t99) - g(3) * (t102 * t66 - t41 * t99) -t117, -g(1) * (pkin(5) * t19 + pkin(12) * t18 + t149) - g(2) * (pkin(5) * t17 + pkin(12) * t16 + t150) - g(3) * (pkin(5) * t41 + pkin(12) * t40 + t131) + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t115, 0, 0, 0, 0, 0, 0, 0, 0, t116 * t103, -t116 * t100, -t115, -g(1) * (-pkin(3) * t34 + pkin(11) * t35) - g(2) * (-pkin(3) * t30 + pkin(11) * t31) - g(3) * (-pkin(3) * t54 + pkin(11) * t55) 0, 0, 0, 0, 0, 0, t116 * t93, -t1, -t115, -g(1) * t168 - g(2) * t169 - g(3) * t167, 0, 0, 0, 0, 0, 0, -g(1) * (-t34 * t158 + t35 * t99) - g(2) * (-t30 * t158 + t31 * t99) - g(3) * (-t54 * t158 + t55 * t99) -g(1) * (t102 * t35 + t34 * t170) - g(2) * (t102 * t31 + t30 * t170) - g(3) * (t102 * t55 + t54 * t170) t1, -g(1) * (t130 * t34 + t168) - g(2) * (t130 * t30 + t169) - g(3) * (t130 * t54 + t167); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 + g(2) * t192 - g(3) * t186, g(1) * t15 - g(2) * t191 - g(3) * (-t100 * t73 - t103 * t55) 0, 0, 0, 0, 0, 0, 0, 0, t119, t118, 0, -g(1) * t135 + g(2) * t136 - g(3) * t134, 0, 0, 0, 0, 0, 0, t119 * t102, -t119 * t99, -t118, -g(1) * (-t12 * pkin(5) + pkin(12) * t13 + t135) - g(2) * (pkin(5) * t184 + pkin(12) * t9 - t136) - g(3) * (pkin(5) * t24 + t25 * pkin(12) + t134); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t194 - g(3) * (t102 * t54 - t25 * t99) g(1) * t3 + g(2) * t193 - g(3) * (-t102 * t25 - t54 * t99) 0, 0;];
taug_reg  = t4;
