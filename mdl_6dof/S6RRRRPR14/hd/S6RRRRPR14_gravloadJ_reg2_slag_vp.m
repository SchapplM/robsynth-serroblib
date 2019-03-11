% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR14
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t150 = cos(pkin(6));
t160 = cos(qJ(2));
t132 = t150 * t160;
t157 = sin(qJ(2));
t158 = sin(qJ(1));
t161 = cos(qJ(1));
t102 = -t161 * t132 + t158 * t157;
t148 = sin(pkin(6));
t149 = cos(pkin(7));
t121 = t149 * t148;
t147 = sin(pkin(7));
t164 = t102 * t147 - t161 * t121;
t159 = cos(qJ(3));
t120 = t148 * t147;
t169 = t102 * t149 + t161 * t120;
t131 = t150 * t157;
t57 = t161 * t131 + t158 * t160;
t81 = sin(qJ(3));
t27 = -t57 * t159 + t169 * t81;
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t11 = -t164 * t80 + t27 * t82;
t24 = t159 * t169 + t57 * t81;
t76 = pkin(13) + qJ(6);
t73 = sin(t76);
t74 = cos(t76);
t173 = t11 * t73 + t24 * t74;
t172 = t11 * t74 - t24 * t73;
t10 = t164 * t82 + t27 * t80;
t95 = t158 * t132 + t161 * t157;
t166 = t158 * t121 + t95 * t147;
t165 = t158 * t120 - t95 * t149;
t162 = t160 * t121 + t150 * t147;
t156 = t73 * t82;
t155 = t74 * t82;
t77 = sin(pkin(13));
t154 = t77 * t82;
t78 = cos(pkin(13));
t153 = t78 * t82;
t105 = t157 * t120;
t128 = t160 * t148;
t152 = pkin(2) * t128 + pkin(10) * t105;
t127 = t148 * t158;
t151 = t161 * pkin(1) + pkin(9) * t127;
t108 = t157 * t121;
t51 = -t81 * t108 + t159 * t128;
t146 = t51 * pkin(3) + t152;
t145 = t77 * pkin(5) + pkin(11);
t18 = t24 * pkin(3);
t144 = -pkin(11) * t27 - t18;
t58 = -t158 * t131 + t161 * t160;
t28 = -t165 * t159 + t58 * t81;
t20 = t28 * pkin(3);
t29 = t58 * t159 + t165 * t81;
t143 = t29 * pkin(11) - t20;
t126 = t148 * t157;
t42 = t81 * t126 - t162 * t159;
t41 = t42 * pkin(3);
t43 = t159 * t126 + t162 * t81;
t142 = t43 * pkin(11) - t41;
t141 = t57 * t147;
t140 = t58 * t147;
t139 = t80 * t147;
t138 = t81 * t149;
t137 = t82 * t147;
t129 = t161 * t148;
t136 = -t158 * pkin(1) + pkin(9) * t129;
t12 = -t166 * t82 + t29 * t80;
t135 = g(1) * t10 + g(2) * t12;
t134 = -g(1) * t24 + g(2) * t28;
t130 = t149 * t159;
t125 = -pkin(4) * t82 - qJ(5) * t80;
t71 = pkin(5) * t78 + pkin(4);
t79 = -pkin(12) - qJ(5);
t124 = -t71 * t82 + t79 * t80;
t50 = t159 * t108 + t81 * t128;
t123 = pkin(11) * t50 + t146;
t119 = -t102 * pkin(2) + pkin(10) * t141;
t118 = -t95 * pkin(2) + pkin(10) * t140;
t33 = -t102 * t159 - t57 * t138;
t117 = t33 * pkin(3) + t119;
t35 = -t58 * t138 - t95 * t159;
t116 = t35 * pkin(3) + t118;
t94 = -t160 * t120 + t150 * t149;
t22 = t43 * t80 - t94 * t82;
t115 = g(1) * t12 - g(2) * t10 + g(3) * t22;
t13 = t166 * t80 + t29 * t82;
t23 = t43 * t82 + t94 * t80;
t114 = g(1) * t13 - g(2) * t11 + g(3) * t23;
t14 = -t57 * t137 + t33 * t80;
t16 = -t58 * t137 + t35 * t80;
t36 = -t82 * t105 + t51 * t80;
t113 = g(1) * t16 + g(2) * t14 + g(3) * t36;
t112 = g(1) * t28 + g(2) * t24 + g(3) * t42;
t111 = g(1) * t29 - g(2) * t27 + g(3) * t43;
t32 = -t102 * t81 + t57 * t130;
t34 = t58 * t130 - t95 * t81;
t110 = g(1) * t34 + g(2) * t32 + g(3) * t50;
t101 = t32 * pkin(11) + t117;
t100 = t34 * pkin(11) + t116;
t91 = -t57 * pkin(2) - t164 * pkin(10) + t136;
t88 = t27 * pkin(3) + t91;
t87 = t58 * pkin(2) + t166 * pkin(10) + t151;
t86 = t29 * pkin(3) + t87;
t85 = -pkin(11) * t24 + t88;
t84 = t28 * pkin(11) + t86;
t37 = t80 * t105 + t51 * t82;
t17 = t58 * t139 + t35 * t82;
t15 = t57 * t139 + t33 * t82;
t7 = t112 * t80;
t6 = t13 * t74 + t28 * t73;
t5 = -t13 * t73 + t28 * t74;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t158 - g(2) * t161, g(1) * t161 + g(2) * t158, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t57 - g(2) * t58, -g(1) * t102 + g(2) * t95, -g(1) * t129 - g(2) * t127, -g(1) * t136 - g(2) * t151, 0, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t29, t134, g(1) * t164 - g(2) * t166, -g(1) * t91 - g(2) * t87, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, t135, -t134, -g(1) * t85 - g(2) * t84, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t78 - t24 * t77) - g(2) * (t13 * t78 + t28 * t77) -g(1) * (-t11 * t77 - t24 * t78) - g(2) * (-t13 * t77 + t28 * t78) -t135, -g(1) * (t11 * pkin(4) + t10 * qJ(5) + t85) - g(2) * (t13 * pkin(4) + t12 * qJ(5) + t84) 0, 0, 0, 0, 0, 0, -g(1) * t172 - g(2) * t6, g(1) * t173 - g(2) * t5, -t135, -g(1) * (-t10 * t79 + t11 * t71 - t145 * t24 + t88) - g(2) * (-t12 * t79 + t13 * t71 + t145 * t28 + t86); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t95 + g(2) * t102 - g(3) * t128, g(1) * t58 + g(2) * t57 + g(3) * t126, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t35 - g(2) * t33 - g(3) * t51, t110, -g(1) * t140 - g(2) * t141 - g(3) * t105, -g(1) * t118 - g(2) * t119 - g(3) * t152, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15 - g(3) * t37, t113, -t110, -g(1) * t100 - g(2) * t101 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t78 + t34 * t77) - g(2) * (t15 * t78 + t32 * t77) - g(3) * (t37 * t78 + t50 * t77) -g(1) * (-t17 * t77 + t34 * t78) - g(2) * (-t15 * t77 + t32 * t78) - g(3) * (-t37 * t77 + t50 * t78) -t113, -g(1) * (t17 * pkin(4) + t16 * qJ(5) + t100) - g(2) * (t15 * pkin(4) + t14 * qJ(5) + t101) - g(3) * (pkin(4) * t37 + qJ(5) * t36 + t123) 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t74 + t34 * t73) - g(2) * (t15 * t74 + t32 * t73) - g(3) * (t37 * t74 + t50 * t73) -g(1) * (-t17 * t73 + t34 * t74) - g(2) * (-t15 * t73 + t32 * t74) - g(3) * (-t37 * t73 + t50 * t74) -t113, -g(1) * (t145 * t34 - t16 * t79 + t17 * t71 + t116) - g(2) * (-t14 * t79 + t145 * t32 + t15 * t71 + t117) - g(3) * (t145 * t50 - t36 * t79 + t37 * t71 + t146); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t111, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t82, -t7, -t111, -g(1) * t143 - g(2) * t144 - g(3) * t142, 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t153 + t29 * t77) - g(2) * (-t24 * t153 - t27 * t77) - g(3) * (-t42 * t153 + t43 * t77) -g(1) * (t28 * t154 + t29 * t78) - g(2) * (t24 * t154 - t27 * t78) - g(3) * (t42 * t154 + t43 * t78) t7, -g(1) * (t125 * t28 + t143) - g(2) * (t125 * t24 + t144) - g(3) * (t125 * t42 + t142) 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t155 + t29 * t73) - g(2) * (-t24 * t155 - t27 * t73) - g(3) * (-t42 * t155 + t43 * t73) -g(1) * (t28 * t156 + t29 * t74) - g(2) * (t24 * t156 - t27 * t74) - g(3) * (t42 * t156 + t43 * t74) t7, -g(1) * (t124 * t28 + t145 * t29 - t20) - g(2) * (t124 * t24 - t145 * t27 - t18) - g(3) * (t124 * t42 + t145 * t43 - t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t114, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t78, -t115 * t77, -t114, -g(1) * (-pkin(4) * t12 + qJ(5) * t13) - g(2) * (pkin(4) * t10 - qJ(5) * t11) - g(3) * (-pkin(4) * t22 + qJ(5) * t23) 0, 0, 0, 0, 0, 0, t115 * t74, -t115 * t73, -t114, -g(1) * (-t12 * t71 - t13 * t79) - g(2) * (t10 * t71 + t11 * t79) - g(3) * (-t22 * t71 - t23 * t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t173 - g(3) * (-t23 * t73 + t42 * t74) g(1) * t6 - g(2) * t172 - g(3) * (-t23 * t74 - t42 * t73) 0, 0;];
taug_reg  = t1;
