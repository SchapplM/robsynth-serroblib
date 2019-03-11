% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t173 = cos(qJ(4));
t161 = cos(pkin(8));
t159 = sin(pkin(6));
t162 = cos(pkin(7));
t135 = t162 * t159;
t158 = sin(pkin(7));
t174 = cos(qJ(1));
t160 = cos(pkin(14));
t163 = cos(pkin(6));
t138 = t163 * t160;
t157 = sin(pkin(14));
t172 = sin(qJ(1));
t72 = -t174 * t138 + t172 * t157;
t115 = -t174 * t135 + t72 * t158;
t81 = sin(pkin(8));
t169 = t115 * t81;
t134 = t159 * t158;
t114 = t174 * t134 + t162 * t72;
t171 = sin(qJ(3));
t136 = t163 * t157;
t73 = t174 * t136 + t172 * t160;
t156 = t73 * t171;
t87 = cos(qJ(3));
t176 = t114 * t87 + t156;
t177 = t161 * t176 - t169;
t54 = -t114 * t171 + t73 * t87;
t84 = sin(qJ(4));
t21 = t54 * t173 - t177 * t84;
t39 = t115 * t161 + t176 * t81;
t83 = sin(qJ(5));
t86 = cos(qJ(5));
t5 = t21 * t86 + t39 * t83;
t82 = sin(qJ(6));
t190 = t5 * t82;
t85 = cos(qJ(6));
t189 = t5 * t85;
t188 = -t21 * t83 + t39 * t86;
t184 = t54 * t84;
t112 = t172 * t138 + t174 * t157;
t180 = t112 * t158 + t172 * t135;
t113 = -t172 * t136 + t174 * t160;
t106 = t113 * t171;
t98 = -t112 * t162 + t172 * t134;
t93 = -t98 * t87 + t106;
t183 = t161 * t180 + t93 * t81;
t179 = t93 * t161 - t180 * t81;
t109 = -t160 * t134 + t163 * t162;
t110 = t160 * t135 + t163 * t158;
t133 = t159 * t157;
t120 = t171 * t133;
t99 = -t110 * t87 + t120;
t178 = t109 * t81 - t99 * t161;
t175 = pkin(11) * t81;
t168 = t81 * t83;
t167 = t81 * t86;
t166 = t82 * t86;
t165 = t85 * t86;
t139 = t159 * t172;
t164 = t174 * pkin(1) + qJ(2) * t139;
t20 = t173 * t177 + t184;
t155 = -t20 * pkin(4) + t21 * pkin(12);
t57 = t113 * t87 + t171 * t98;
t24 = t173 * t179 + t57 * t84;
t25 = t57 * t173 - t179 * t84;
t154 = -t24 * pkin(4) + t25 * pkin(12);
t66 = t110 * t171 + t87 * t133;
t35 = -t173 * t178 + t66 * t84;
t36 = t66 * t173 + t178 * t84;
t153 = -t35 * pkin(4) + t36 * pkin(12);
t150 = t84 * t161;
t148 = -pkin(3) * t176 + t54 * t175;
t147 = -t93 * pkin(3) + t57 * t175;
t146 = -t99 * pkin(3) + t66 * t175;
t8 = -t183 * t86 + t25 * t83;
t145 = g(1) * t188 + g(2) * t8;
t140 = t174 * t159;
t144 = -t172 * pkin(1) + qJ(2) * t140;
t143 = -pkin(5) * t86 - pkin(13) * t83;
t141 = t161 * t173;
t22 = -t141 * t176 + t173 * t169 - t184;
t142 = g(1) * t22 + g(2) * t24;
t49 = t109 * t161 + t99 * t81;
t14 = -t36 * t83 + t49 * t86;
t132 = g(1) * t8 - g(2) * t188 - g(3) * t14;
t15 = t36 * t86 + t49 * t83;
t9 = t183 * t83 + t25 * t86;
t131 = g(1) * t9 + g(2) * t5 + g(3) * t15;
t29 = -t54 * t150 - t173 * t176;
t10 = -t54 * t167 + t29 * t83;
t31 = -t57 * t150 - t93 * t173;
t12 = -t57 * t167 + t31 * t83;
t43 = -t66 * t150 - t99 * t173;
t32 = -t66 * t167 + t43 * t83;
t130 = g(1) * t12 + g(2) * t10 + g(3) * t32;
t129 = g(1) * t24 + g(2) * t20 + g(3) * t35;
t128 = g(1) * t25 + g(2) * t21 + g(3) * t36;
t28 = t54 * t141 - t176 * t84;
t30 = t57 * t141 - t93 * t84;
t42 = t66 * t141 - t99 * t84;
t127 = g(1) * t30 + g(2) * t28 + g(3) * t42;
t126 = g(1) * t57 + g(2) * t54 + g(3) * t66;
t125 = t29 * pkin(4) + t28 * pkin(12) + t148;
t124 = t31 * pkin(4) + t30 * pkin(12) + t147;
t123 = t43 * pkin(4) + t42 * pkin(12) + t146;
t116 = -t73 * pkin(2) - t115 * pkin(10) + t144;
t107 = -t54 * pkin(3) - pkin(11) * t39 + t116;
t102 = -pkin(4) * t21 + t22 * pkin(12) + t107;
t100 = t113 * pkin(2) + t180 * pkin(10) + t164;
t90 = t57 * pkin(3) + t183 * pkin(11) + t100;
t89 = t25 * pkin(4) + t24 * pkin(12) + t90;
t69 = -g(1) * t139 + g(2) * t140 - g(3) * t163;
t33 = t66 * t168 + t43 * t86;
t13 = t57 * t168 + t31 * t86;
t11 = t54 * t168 + t29 * t86;
t3 = t24 * t82 + t9 * t85;
t2 = t24 * t85 - t9 * t82;
t1 = t129 * t83;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t172 - g(2) * t174, g(1) * t174 + g(2) * t172, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t73 - g(2) * t113, -g(1) * t72 + g(2) * t112, -g(1) * t140 - g(2) * t139, -g(1) * t144 - g(2) * t164, 0, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t57, -g(1) * t176 + g(2) * t93, g(1) * t115 - g(2) * t180, -g(1) * t116 - g(2) * t100, 0, 0, 0, 0, 0, 0, g(1) * t21 - g(2) * t25, t142, g(1) * t39 - g(2) * t183, -g(1) * t107 - g(2) * t90, 0, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t145, -t142, -g(1) * t102 - g(2) * t89, 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t82 - t189) - g(2) * t3, -g(1) * (t22 * t85 + t190) - g(2) * t2, -t145, -g(1) * (-pkin(5) * t5 + pkin(13) * t188 + t102) - g(2) * (t9 * pkin(5) + t8 * pkin(13) + t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t106 + g(2) * t156 + g(3) * t120 + (-g(1) * t98 + g(2) * t114 - g(3) * t110) * t87, t126, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t29 - g(3) * t43, t127, -t126 * t81, -g(1) * t147 - g(2) * t148 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t33, t130, -t127, -g(1) * t124 - g(2) * t125 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t85 + t30 * t82) - g(2) * (t11 * t85 + t28 * t82) - g(3) * (t33 * t85 + t42 * t82) -g(1) * (-t13 * t82 + t30 * t85) - g(2) * (-t11 * t82 + t28 * t85) - g(3) * (-t33 * t82 + t42 * t85) -t130, -g(1) * (t13 * pkin(5) + t12 * pkin(13) + t124) - g(2) * (t11 * pkin(5) + t10 * pkin(13) + t125) - g(3) * (t33 * pkin(5) + t32 * pkin(13) + t123); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t128, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t86, -t1, -t128, -g(1) * t154 - g(2) * t155 - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t165 + t25 * t82) - g(2) * (-t20 * t165 + t21 * t82) - g(3) * (-t35 * t165 + t36 * t82) -g(1) * (t24 * t166 + t25 * t85) - g(2) * (t20 * t166 + t21 * t85) - g(3) * (t35 * t166 + t36 * t85) t1, -g(1) * (t143 * t24 + t154) - g(2) * (t143 * t20 + t155) - g(3) * (t143 * t35 + t153); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t131, 0, 0, 0, 0, 0, 0, 0, 0, t132 * t85, -t132 * t82, -t131, -g(1) * (-t8 * pkin(5) + t9 * pkin(13)) - g(2) * (pkin(5) * t188 + t5 * pkin(13)) - g(3) * (t14 * pkin(5) + t15 * pkin(13)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * (t20 * t85 - t190) - g(3) * (-t15 * t82 + t35 * t85) g(1) * t3 - g(2) * (-t20 * t82 - t189) - g(3) * (-t15 * t85 - t35 * t82) 0, 0;];
taug_reg  = t4;
