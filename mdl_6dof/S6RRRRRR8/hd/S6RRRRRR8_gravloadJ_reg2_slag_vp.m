% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR8
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t105 = sin(qJ(6));
t108 = cos(qJ(6));
t102 = qJ(4) + qJ(5);
t100 = cos(t102);
t103 = sin(pkin(7));
t110 = cos(qJ(1));
t104 = sin(pkin(6));
t164 = cos(pkin(7));
t147 = t104 * t164;
t165 = cos(pkin(6));
t185 = cos(qJ(2));
t131 = t165 * t185;
t182 = sin(qJ(2));
t183 = sin(qJ(1));
t81 = -t110 * t131 + t183 * t182;
t193 = -t81 * t103 + t110 * t147;
t107 = sin(qJ(3));
t146 = t107 * t164;
t158 = t104 * t110;
t155 = t103 * t158;
t184 = cos(qJ(3));
t130 = t165 * t182;
t82 = t110 * t130 + t183 * t185;
t38 = -t107 * t155 - t81 * t146 + t82 * t184;
t99 = sin(t102);
t15 = t100 * t38 - t193 * t99;
t129 = t164 * t184;
t37 = t82 * t107 + t81 * t129 + t184 * t155;
t200 = t105 * t15 - t108 * t37;
t199 = t105 * t37 + t108 * t15;
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t198 = t106 * t38 + t109 * t193;
t172 = t106 * t193;
t197 = -t109 * t38 + t172;
t118 = t110 * t182 + t183 * t131;
t65 = t118 * t103 + t183 * t147;
t145 = t165 * t103;
t62 = t107 * t145 + (t185 * t146 + t184 * t182) * t104;
t150 = t104 * t185;
t80 = -t103 * t150 + t165 * t164;
t192 = -t106 * t62 + t109 * t80;
t149 = t104 * t183;
t189 = -t103 * t149 + t118 * t164;
t83 = t110 * t185 - t183 * t130;
t42 = -t189 * t107 + t83 * t184;
t20 = -t106 * t42 + t109 * t65;
t190 = -t100 * t193 - t38 * t99;
t188 = -g(1) * t83 - g(2) * t82;
t181 = pkin(10) * t103;
t111 = -pkin(12) - pkin(11);
t98 = pkin(4) * t109 + pkin(3);
t180 = -t38 * t111 - t37 * t98;
t41 = t107 * t83 + t189 * t184;
t179 = -t42 * t111 - t41 * t98;
t148 = t104 * t182;
t61 = t107 * t148 - t129 * t150 - t184 * t145;
t178 = -t62 * t111 - t61 * t98;
t139 = t103 * t148;
t177 = pkin(2) * t150 + pkin(10) * t139;
t176 = t103 * t99;
t171 = t106 * t65;
t166 = t110 * pkin(1) + pkin(9) * t149;
t163 = t100 * t103;
t162 = t100 * t105;
t161 = t100 * t108;
t160 = t103 * t106;
t159 = t103 * t109;
t47 = -t107 * t81 + t82 * t129;
t48 = -t82 * t146 - t81 * t184;
t76 = t81 * pkin(2);
t157 = -t47 * t111 + t48 * t98 - t76;
t49 = -t118 * t107 + t83 * t129;
t50 = -t118 * t184 - t83 * t146;
t78 = t118 * pkin(2);
t156 = -t49 * t111 + t50 * t98 - t78;
t153 = pkin(5) * t190 + pkin(13) * t15;
t18 = -t65 * t100 + t42 * t99;
t19 = t100 * t42 + t65 * t99;
t152 = -t18 * pkin(5) + pkin(13) * t19;
t31 = t100 * t80 - t62 * t99;
t32 = t100 * t62 + t80 * t99;
t151 = t31 * pkin(5) + pkin(13) * t32;
t144 = t198 * pkin(4);
t143 = t20 * pkin(4);
t142 = t192 * pkin(4);
t141 = t82 * t181 - t76;
t140 = t83 * t181 - t78;
t137 = -t183 * pkin(1) + pkin(9) * t158;
t127 = t106 * t139;
t128 = t164 * t182;
t73 = (t185 * t107 + t184 * t128) * t104;
t74 = (-t107 * t128 + t184 * t185) * t104;
t136 = pkin(4) * t127 - t73 * t111 + t74 * t98 + t177;
t135 = g(1) * t190 + g(2) * t18;
t134 = -g(1) * t37 + g(2) * t41;
t132 = -pkin(5) * t100 - pkin(13) * t99;
t125 = -g(1) * t110 - g(2) * t183;
t3 = g(1) * t18 - g(2) * t190 - g(3) * t31;
t5 = g(1) * t19 + g(2) * t15 + g(3) * t32;
t22 = -t82 * t163 + t48 * t99;
t24 = -t83 * t163 + t50 * t99;
t51 = -t100 * t139 + t74 * t99;
t124 = g(1) * t24 + g(2) * t22 + g(3) * t51;
t123 = g(1) * t41 + g(2) * t37 + g(3) * t61;
t122 = g(1) * t42 + g(2) * t38 + g(3) * t62;
t121 = g(1) * t49 + g(2) * t47 + g(3) * t73;
t120 = -g(3) * t148 + t188;
t119 = -t82 * pkin(2) + t193 * pkin(10) + t137;
t117 = pkin(4) * t172 + t111 * t37 - t38 * t98 + t119;
t115 = t188 * (pkin(4) * t106 + pkin(10)) * t103;
t113 = t83 * pkin(2) + t65 * pkin(10) + t166;
t112 = pkin(4) * t171 - t41 * t111 + t42 * t98 + t113;
t52 = t74 * t100 + t99 * t139;
t25 = t100 * t50 + t83 * t176;
t23 = t100 * t48 + t82 * t176;
t21 = t109 * t42 + t171;
t8 = t105 * t41 + t108 * t19;
t7 = -t105 * t19 + t108 * t41;
t6 = t123 * t99;
t2 = t3 * t108;
t1 = t3 * t105;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t183 - g(2) * t110, -t125, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t82 - g(2) * t83, -g(1) * t81 + g(2) * t118, t125 * t104, -g(1) * t137 - g(2) * t166, 0, 0, 0, 0, 0, 0, g(1) * t38 - g(2) * t42, t134, -g(1) * t193 - g(2) * t65, -g(1) * t119 - g(2) * t113, 0, 0, 0, 0, 0, 0, -g(1) * t197 - g(2) * t21, -g(1) * t198 - g(2) * t20, -t134, -g(1) * (-pkin(3) * t38 - pkin(11) * t37 + t119) - g(2) * (t42 * pkin(3) + t41 * pkin(11) + t113) 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t19, t135, -t134, -g(1) * t117 - g(2) * t112, 0, 0, 0, 0, 0, 0, g(1) * t199 - g(2) * t8, -g(1) * t200 - g(2) * t7, -t135, -g(1) * (-pkin(5) * t15 + pkin(13) * t190 + t117) - g(2) * (t19 * pkin(5) + t18 * pkin(13) + t112); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t118 + g(2) * t81 - g(3) * t150, -t120, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t50 - g(2) * t48 - g(3) * t74, t121, t120 * t103, -g(1) * t140 - g(2) * t141 - g(3) * t177, 0, 0, 0, 0, 0, 0, -g(1) * (t109 * t50 + t83 * t160) - g(2) * (t109 * t48 + t82 * t160) - g(3) * (t74 * t109 + t127) -g(1) * (-t106 * t50 + t83 * t159) - g(2) * (-t106 * t48 + t82 * t159) - g(3) * (-t74 * t106 + t109 * t139) -t121, -g(1) * (pkin(3) * t50 + pkin(11) * t49 + t140) - g(2) * (pkin(3) * t48 + pkin(11) * t47 + t141) - g(3) * (pkin(3) * t74 + pkin(11) * t73 + t177) 0, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t23 - g(3) * t52, t124, -t121, -g(1) * t156 - g(2) * t157 - g(3) * t136 + t115, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t49 + t108 * t25) - g(2) * (t105 * t47 + t108 * t23) - g(3) * (t105 * t73 + t108 * t52) -g(1) * (-t105 * t25 + t108 * t49) - g(2) * (-t105 * t23 + t108 * t47) - g(3) * (-t105 * t52 + t108 * t73) -t124, -g(1) * (pkin(5) * t25 + pkin(13) * t24 + t156) - g(2) * (pkin(5) * t23 + pkin(13) * t22 + t157) - g(3) * (pkin(5) * t52 + pkin(13) * t51 + t136) + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t122, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t109, -t123 * t106, -t122, -g(1) * (-pkin(3) * t41 + pkin(11) * t42) - g(2) * (-pkin(3) * t37 + pkin(11) * t38) - g(3) * (-pkin(3) * t61 + pkin(11) * t62) 0, 0, 0, 0, 0, 0, t123 * t100, -t6, -t122, -g(1) * t179 - g(2) * t180 - g(3) * t178, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t42 - t41 * t161) - g(2) * (t105 * t38 - t37 * t161) - g(3) * (t105 * t62 - t61 * t161) -g(1) * (t108 * t42 + t41 * t162) - g(2) * (t108 * t38 + t37 * t162) - g(3) * (t108 * t62 + t61 * t162) t6, -g(1) * (t132 * t41 + t179) - g(2) * (t132 * t37 + t180) - g(3) * (t132 * t61 + t178); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + g(2) * t198 - g(3) * t192, g(1) * t21 - g(2) * t197 - g(3) * (-t106 * t80 - t109 * t62) 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t143 + g(2) * t144 - g(3) * t142, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (t143 + t152) - g(2) * (-t144 + t153) - g(3) * (t142 + t151); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t152 - g(2) * t153 - g(3) * t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t200 - g(3) * (-t105 * t32 + t108 * t61) g(1) * t8 + g(2) * t199 - g(3) * (-t105 * t61 - t108 * t32) 0, 0;];
taug_reg  = t4;
