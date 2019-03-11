% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t157 = cos(pkin(6));
t165 = cos(qJ(2));
t138 = t157 * t165;
t162 = sin(qJ(2));
t163 = sin(qJ(1));
t166 = cos(qJ(1));
t109 = -t166 * t138 + t163 * t162;
t155 = sin(pkin(6));
t156 = cos(pkin(7));
t128 = t156 * t155;
t154 = sin(pkin(7));
t171 = t109 * t154 - t166 * t128;
t164 = cos(qJ(3));
t127 = t155 * t154;
t174 = t109 * t156 + t166 * t127;
t137 = t157 * t162;
t66 = t166 * t137 + t163 * t165;
t86 = sin(qJ(3));
t36 = -t66 * t164 + t174 * t86;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t20 = -t171 * t85 + t36 * t88;
t33 = t174 * t164 + t66 * t86;
t84 = sin(qJ(5));
t87 = cos(qJ(5));
t178 = t20 * t84 + t33 * t87;
t177 = t20 * t87 - t33 * t84;
t19 = t171 * t88 + t36 * t85;
t102 = t163 * t138 + t166 * t162;
t89 = -t102 * t154 - t163 * t128;
t170 = -t102 * t156 + t163 * t127;
t67 = -t163 * t137 + t166 * t165;
t38 = t67 * t164 + t170 * t86;
t22 = t38 * t88 - t85 * t89;
t37 = -t170 * t164 + t67 * t86;
t15 = -t22 * t84 + t37 * t87;
t100 = -t165 * t127 + t157 * t156;
t132 = t155 * t162;
t168 = t165 * t128 + t157 * t154;
t52 = t164 * t132 + t168 * t86;
t32 = t100 * t85 + t52 * t88;
t51 = t86 * t132 - t168 * t164;
t1 = -g(1) * t15 - g(2) * t178 - g(3) * (-t32 * t84 + t51 * t87);
t161 = t84 * t88;
t160 = t87 * t88;
t112 = t162 * t127;
t134 = t165 * t155;
t159 = pkin(2) * t134 + pkin(10) * t112;
t133 = t155 * t163;
t158 = t166 * pkin(1) + pkin(9) * t133;
t115 = t162 * t128;
t60 = -t86 * t115 + t164 * t134;
t153 = t60 * pkin(3) + t159;
t152 = t84 * pkin(5) + pkin(11);
t27 = t33 * pkin(3);
t151 = -pkin(11) * t36 - t27;
t29 = t37 * pkin(3);
t150 = t38 * pkin(11) - t29;
t50 = t51 * pkin(3);
t149 = t52 * pkin(11) - t50;
t148 = t66 * t154;
t147 = t67 * t154;
t146 = t85 * t154;
t145 = t86 * t156;
t144 = t88 * t154;
t135 = t166 * t155;
t143 = -t163 * pkin(1) + pkin(9) * t135;
t142 = -pkin(4) * t88 - pkin(12) * t85;
t21 = t38 * t85 + t88 * t89;
t141 = g(1) * t19 + g(2) * t21;
t140 = -g(1) * t33 + g(2) * t37;
t136 = t156 * t164;
t81 = pkin(5) * t87 + pkin(4);
t83 = -qJ(6) - pkin(12);
t131 = -t81 * t88 + t83 * t85;
t59 = t164 * t115 + t86 * t134;
t130 = pkin(11) * t59 + t153;
t126 = -t109 * pkin(2) + pkin(10) * t148;
t125 = -t102 * pkin(2) + pkin(10) * t147;
t42 = -t109 * t164 - t66 * t145;
t124 = t42 * pkin(3) + t126;
t44 = -t102 * t164 - t67 * t145;
t123 = t44 * pkin(3) + t125;
t31 = -t100 * t88 + t52 * t85;
t122 = g(1) * t21 - g(2) * t19 + g(3) * t31;
t121 = g(1) * t22 - g(2) * t20 + g(3) * t32;
t23 = -t66 * t144 + t42 * t85;
t25 = -t67 * t144 + t44 * t85;
t45 = -t88 * t112 + t60 * t85;
t120 = g(1) * t25 + g(2) * t23 + g(3) * t45;
t119 = g(1) * t37 + g(2) * t33 + g(3) * t51;
t118 = g(1) * t38 - g(2) * t36 + g(3) * t52;
t41 = -t109 * t86 + t66 * t136;
t43 = -t102 * t86 + t67 * t136;
t117 = g(1) * t43 + g(2) * t41 + g(3) * t59;
t108 = t41 * pkin(11) + t124;
t107 = t43 * pkin(11) + t123;
t97 = -t66 * pkin(2) - t171 * pkin(10) + t143;
t94 = t36 * pkin(3) + t97;
t93 = t67 * pkin(2) - t89 * pkin(10) + t158;
t92 = t38 * pkin(3) + t93;
t91 = -pkin(11) * t33 + t94;
t90 = t37 * pkin(11) + t92;
t46 = t85 * t112 + t60 * t88;
t26 = t67 * t146 + t44 * t88;
t24 = t66 * t146 + t42 * t88;
t16 = t22 * t87 + t37 * t84;
t14 = t119 * t85;
t10 = t122 * t87;
t9 = t122 * t84;
t8 = -g(1) * (t26 * t87 + t43 * t84) - g(2) * (t24 * t87 + t41 * t84) - g(3) * (t46 * t87 + t59 * t84);
t7 = -g(1) * (-t26 * t84 + t43 * t87) - g(2) * (-t24 * t84 + t41 * t87) - g(3) * (-t46 * t84 + t59 * t87);
t6 = -g(1) * t177 - g(2) * t16;
t5 = g(1) * t178 - g(2) * t15;
t4 = -g(1) * (-t37 * t160 + t38 * t84) - g(2) * (-t33 * t160 - t36 * t84) - g(3) * (-t51 * t160 + t52 * t84);
t3 = -g(1) * (t37 * t161 + t38 * t87) - g(2) * (t33 * t161 - t36 * t87) - g(3) * (t51 * t161 + t52 * t87);
t2 = g(1) * t16 - g(2) * t177 - g(3) * (-t32 * t87 - t51 * t84);
t11 = [0, 0, 0, 0, 0, 0, g(1) * t163 - g(2) * t166, g(1) * t166 + g(2) * t163, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t66 - g(2) * t67, -g(1) * t109 + g(2) * t102, -g(1) * t135 - g(2) * t133, -g(1) * t143 - g(2) * t158, 0, 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t38, t140, g(1) * t171 + g(2) * t89, -g(1) * t97 - g(2) * t93, 0, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t22, t141, -t140, -g(1) * t91 - g(2) * t90, 0, 0, 0, 0, 0, 0, t6, t5, -t141, -g(1) * (t20 * pkin(4) + t19 * pkin(12) + t91) - g(2) * (t22 * pkin(4) + t21 * pkin(12) + t90) 0, 0, 0, 0, 0, 0, t6, t5, -t141, -g(1) * (-t152 * t33 - t19 * t83 + t20 * t81 + t94) - g(2) * (t152 * t37 - t21 * t83 + t22 * t81 + t92); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t102 + g(2) * t109 - g(3) * t134, g(1) * t67 + g(2) * t66 + g(3) * t132, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t44 - g(2) * t42 - g(3) * t60, t117, -g(1) * t147 - g(2) * t148 - g(3) * t112, -g(1) * t125 - g(2) * t126 - g(3) * t159, 0, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t24 - g(3) * t46, t120, -t117, -g(1) * t107 - g(2) * t108 - g(3) * t130, 0, 0, 0, 0, 0, 0, t8, t7, -t120, -g(1) * (t26 * pkin(4) + t25 * pkin(12) + t107) - g(2) * (t24 * pkin(4) + t23 * pkin(12) + t108) - g(3) * (pkin(4) * t46 + pkin(12) * t45 + t130) 0, 0, 0, 0, 0, 0, t8, t7, -t120, -g(1) * (t152 * t43 - t25 * t83 + t26 * t81 + t123) - g(2) * (t152 * t41 - t23 * t83 + t24 * t81 + t124) - g(3) * (t152 * t59 - t45 * t83 + t46 * t81 + t153); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t118, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t88, -t14, -t118, -g(1) * t150 - g(2) * t151 - g(3) * t149, 0, 0, 0, 0, 0, 0, t4, t3, t14, -g(1) * (t142 * t37 + t150) - g(2) * (t142 * t33 + t151) - g(3) * (t142 * t51 + t149) 0, 0, 0, 0, 0, 0, t4, t3, t14, -g(1) * (t131 * t37 + t152 * t38 - t29) - g(2) * (t131 * t33 - t152 * t36 - t27) - g(3) * (t131 * t51 + t152 * t52 - t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t121, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t121, -g(1) * (-pkin(4) * t21 + pkin(12) * t22) - g(2) * (pkin(4) * t19 - pkin(12) * t20) - g(3) * (-pkin(4) * t31 + pkin(12) * t32) 0, 0, 0, 0, 0, 0, t10, -t9, -t121, -g(1) * (-t21 * t81 - t22 * t83) - g(2) * (t19 * t81 + t20 * t83) - g(3) * (-t31 * t81 - t32 * t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122;];
taug_reg  = t11;
