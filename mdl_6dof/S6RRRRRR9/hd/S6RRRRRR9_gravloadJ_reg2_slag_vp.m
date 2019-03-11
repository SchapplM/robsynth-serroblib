% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t152 = cos(pkin(6));
t162 = cos(qJ(2));
t133 = t152 * t162;
t159 = sin(qJ(2));
t160 = sin(qJ(1));
t163 = cos(qJ(1));
t104 = -t163 * t133 + t160 * t159;
t150 = sin(pkin(6));
t151 = cos(pkin(7));
t123 = t151 * t150;
t149 = sin(pkin(7));
t167 = t104 * t149 - t163 * t123;
t161 = cos(qJ(3));
t122 = t150 * t149;
t173 = t104 * t151 + t163 * t122;
t132 = t152 * t159;
t60 = t163 * t132 + t160 * t162;
t82 = sin(qJ(3));
t30 = -t60 * t161 + t173 * t82;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t14 = -t167 * t81 + t30 * t84;
t27 = t173 * t161 + t60 * t82;
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t181 = t14 * t80 + t27 * t83;
t180 = t14 * t83 - t27 * t80;
t79 = qJ(5) + qJ(6);
t76 = sin(t79);
t77 = cos(t79);
t179 = t14 * t76 + t27 * t77;
t178 = t14 * t77 - t27 * t76;
t172 = t167 * t84 + t30 * t81;
t97 = t160 * t133 + t163 * t159;
t86 = -t160 * t123 - t97 * t149;
t169 = t160 * t122 - t97 * t151;
t127 = t150 * t159;
t165 = t162 * t123 + t152 * t149;
t45 = t161 * t127 + t165 * t82;
t59 = -t162 * t122 + t152 * t151;
t26 = t45 * t84 + t59 * t81;
t44 = t82 * t127 - t165 * t161;
t61 = -t160 * t132 + t163 * t162;
t32 = t61 * t161 + t169 * t82;
t16 = t32 * t84 - t86 * t81;
t31 = -t169 * t161 + t61 * t82;
t9 = -t16 * t80 + t31 * t83;
t168 = -g(1) * t9 - g(2) * t181 - g(3) * (-t26 * t80 + t44 * t83);
t158 = t76 * t84;
t157 = t77 * t84;
t156 = t80 * t84;
t155 = t83 * t84;
t107 = t159 * t122;
t129 = t162 * t150;
t154 = pkin(2) * t129 + pkin(10) * t107;
t128 = t150 * t160;
t153 = t163 * pkin(1) + pkin(9) * t128;
t110 = t159 * t123;
t54 = -t82 * t110 + t161 * t129;
t148 = t54 * pkin(3) + t154;
t147 = t80 * pkin(5) + pkin(11);
t21 = t27 * pkin(3);
t146 = -pkin(11) * t30 - t21;
t23 = t31 * pkin(3);
t145 = t32 * pkin(11) - t23;
t43 = t44 * pkin(3);
t144 = t45 * pkin(11) - t43;
t143 = t60 * t149;
t142 = t61 * t149;
t141 = t81 * t149;
t140 = t82 * t151;
t139 = t84 * t149;
t130 = t163 * t150;
t138 = -t160 * pkin(1) + pkin(9) * t130;
t137 = -pkin(4) * t84 - pkin(12) * t81;
t15 = t32 * t81 + t86 * t84;
t136 = g(1) * t172 + g(2) * t15;
t135 = -g(1) * t27 + g(2) * t31;
t131 = t151 * t161;
t75 = t83 * pkin(5) + pkin(4);
t85 = -pkin(13) - pkin(12);
t126 = -t75 * t84 + t81 * t85;
t53 = t161 * t110 + t82 * t129;
t125 = t53 * pkin(11) + t148;
t121 = -t104 * pkin(2) + pkin(10) * t143;
t120 = -t97 * pkin(2) + pkin(10) * t142;
t36 = -t104 * t161 - t60 * t140;
t119 = t36 * pkin(3) + t121;
t38 = -t61 * t140 - t97 * t161;
t118 = t38 * pkin(3) + t120;
t25 = -t45 * t81 + t59 * t84;
t117 = g(1) * t15 - g(2) * t172 - g(3) * t25;
t116 = g(1) * t16 - g(2) * t14 + g(3) * t26;
t17 = -t60 * t139 + t36 * t81;
t19 = -t61 * t139 + t38 * t81;
t39 = -t84 * t107 + t54 * t81;
t115 = g(1) * t19 + g(2) * t17 + g(3) * t39;
t114 = g(1) * t31 + g(2) * t27 + g(3) * t44;
t113 = g(1) * t32 - g(2) * t30 + g(3) * t45;
t35 = -t104 * t82 + t60 * t131;
t37 = t61 * t131 - t97 * t82;
t112 = g(1) * t37 + g(2) * t35 + g(3) * t53;
t103 = t35 * pkin(11) + t119;
t102 = t37 * pkin(11) + t118;
t93 = -t60 * pkin(2) - t167 * pkin(10) + t138;
t91 = t30 * pkin(3) + t93;
t90 = t61 * pkin(2) - t86 * pkin(10) + t153;
t89 = t32 * pkin(3) + t90;
t88 = -pkin(11) * t27 + t91;
t87 = t31 * pkin(11) + t89;
t40 = t81 * t107 + t54 * t84;
t20 = t61 * t141 + t38 * t84;
t18 = t60 * t141 + t36 * t84;
t10 = t16 * t83 + t31 * t80;
t8 = t114 * t81;
t7 = t16 * t77 + t31 * t76;
t6 = -t16 * t76 + t31 * t77;
t2 = g(1) * t7 - g(2) * t178 - g(3) * (-t26 * t77 - t44 * t76);
t1 = -g(1) * t6 - g(2) * t179 - g(3) * (-t26 * t76 + t44 * t77);
t3 = [0, 0, 0, 0, 0, 0, g(1) * t160 - g(2) * t163, g(1) * t163 + g(2) * t160, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t61, -g(1) * t104 + g(2) * t97, -g(1) * t130 - g(2) * t128, -g(1) * t138 - g(2) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t32, t135, g(1) * t167 + g(2) * t86, -g(1) * t93 - g(2) * t90, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, t136, -t135, -g(1) * t88 - g(2) * t87, 0, 0, 0, 0, 0, 0, -g(1) * t180 - g(2) * t10, g(1) * t181 - g(2) * t9, -t136, -g(1) * (t14 * pkin(4) + pkin(12) * t172 + t88) - g(2) * (t16 * pkin(4) + t15 * pkin(12) + t87) 0, 0, 0, 0, 0, 0, -g(1) * t178 - g(2) * t7, g(1) * t179 - g(2) * t6, -t136, -g(1) * (t14 * t75 - t147 * t27 - t172 * t85 + t91) - g(2) * (t147 * t31 - t15 * t85 + t16 * t75 + t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t97 + g(2) * t104 - g(3) * t129, g(1) * t61 + g(2) * t60 + g(3) * t127, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t38 - g(2) * t36 - g(3) * t54, t112, -g(1) * t142 - g(2) * t143 - g(3) * t107, -g(1) * t120 - g(2) * t121 - g(3) * t154, 0, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t18 - g(3) * t40, t115, -t112, -g(1) * t102 - g(2) * t103 - g(3) * t125, 0, 0, 0, 0, 0, 0, -g(1) * (t20 * t83 + t37 * t80) - g(2) * (t18 * t83 + t35 * t80) - g(3) * (t40 * t83 + t53 * t80) -g(1) * (-t20 * t80 + t37 * t83) - g(2) * (-t18 * t80 + t35 * t83) - g(3) * (-t40 * t80 + t53 * t83) -t115, -g(1) * (t20 * pkin(4) + t19 * pkin(12) + t102) - g(2) * (t18 * pkin(4) + t17 * pkin(12) + t103) - g(3) * (t40 * pkin(4) + t39 * pkin(12) + t125) 0, 0, 0, 0, 0, 0, -g(1) * (t20 * t77 + t37 * t76) - g(2) * (t18 * t77 + t35 * t76) - g(3) * (t40 * t77 + t53 * t76) -g(1) * (-t20 * t76 + t37 * t77) - g(2) * (-t18 * t76 + t35 * t77) - g(3) * (-t40 * t76 + t53 * t77) -t115, -g(1) * (t147 * t37 - t19 * t85 + t20 * t75 + t118) - g(2) * (t147 * t35 - t17 * t85 + t18 * t75 + t119) - g(3) * (t147 * t53 - t39 * t85 + t40 * t75 + t148); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t113, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t84, -t8, -t113, -g(1) * t145 - g(2) * t146 - g(3) * t144, 0, 0, 0, 0, 0, 0, -g(1) * (-t31 * t155 + t32 * t80) - g(2) * (-t27 * t155 - t30 * t80) - g(3) * (-t44 * t155 + t45 * t80) -g(1) * (t31 * t156 + t32 * t83) - g(2) * (t27 * t156 - t30 * t83) - g(3) * (t44 * t156 + t45 * t83) t8, -g(1) * (t137 * t31 + t145) - g(2) * (t137 * t27 + t146) - g(3) * (t137 * t44 + t144) 0, 0, 0, 0, 0, 0, -g(1) * (-t31 * t157 + t32 * t76) - g(2) * (-t27 * t157 - t30 * t76) - g(3) * (-t44 * t157 + t45 * t76) -g(1) * (t31 * t158 + t32 * t77) - g(2) * (t27 * t158 - t30 * t77) - g(3) * (t44 * t158 + t45 * t77) t8, -g(1) * (t126 * t31 + t147 * t32 - t23) - g(2) * (t126 * t27 - t147 * t30 - t21) - g(3) * (t126 * t44 + t147 * t45 - t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t116, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t83, -t117 * t80, -t116, -g(1) * (-pkin(4) * t15 + pkin(12) * t16) - g(2) * (pkin(4) * t172 - pkin(12) * t14) - g(3) * (t25 * pkin(4) + t26 * pkin(12)) 0, 0, 0, 0, 0, 0, t117 * t77, -t117 * t76, -t116, -g(1) * (-t15 * t75 - t16 * t85) - g(2) * (t14 * t85 + t172 * t75) - g(3) * (t25 * t75 - t26 * t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, g(1) * t10 - g(2) * t180 - g(3) * (-t26 * t83 - t44 * t80) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t168 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
