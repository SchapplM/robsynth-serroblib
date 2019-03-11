% Calculate Gravitation load on the joints for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:57:20
% EndTime: 2019-03-09 14:57:27
% DurationCPUTime: 3.13s
% Computational Cost: add. (2526->280), mult. (7133->426), div. (0->0), fcn. (9298->18), ass. (0->134)
t101 = sin(qJ(6));
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t103 = sin(qJ(4));
t185 = cos(qJ(4));
t100 = sin(pkin(7));
t171 = sin(pkin(6));
t174 = cos(pkin(7));
t150 = t174 * t171;
t187 = cos(qJ(1));
t137 = t187 * t150;
t175 = cos(pkin(6));
t186 = cos(qJ(2));
t158 = t175 * t186;
t183 = sin(qJ(2));
t184 = sin(qJ(1));
t85 = -t158 * t187 + t183 * t184;
t129 = t85 * t100 - t137;
t173 = cos(pkin(8));
t172 = cos(pkin(14));
t148 = t172 * t171;
t141 = t100 * t148;
t151 = t174 * t172;
t170 = sin(pkin(14));
t157 = t175 * t183;
t86 = t157 * t187 + t184 * t186;
t55 = t141 * t187 + t151 * t85 + t170 * t86;
t99 = sin(pkin(8));
t190 = -t129 * t99 + t173 * t55;
t147 = t171 * t170;
t140 = t100 * t147;
t149 = t174 * t170;
t54 = -t140 * t187 - t149 * t85 + t172 * t86;
t22 = t103 * t190 - t54 * t185;
t38 = t129 * t173 + t55 * t99;
t4 = t38 * t102 - t105 * t22;
t198 = t101 * t4;
t104 = cos(qJ(6));
t197 = t104 * t4;
t196 = t102 * t22 + t105 * t38;
t195 = t103 * t54;
t126 = t158 * t184 + t183 * t187;
t121 = t126 * t172;
t87 = -t157 * t184 + t186 * t187;
t111 = t121 * t174 - t141 * t184 + t170 * t87;
t117 = t126 * t100 + t184 * t150;
t106 = t111 * t99 + t117 * t173;
t188 = rSges(7,3) + pkin(13);
t192 = t111 * t173 - t117 * t99;
t135 = t174 * t148;
t163 = t100 * t175;
t114 = -t135 * t186 + t147 * t183 - t163 * t172;
t154 = t186 * t171;
t127 = -t100 * t154 + t174 * t175;
t191 = t114 * t173 - t127 * t99;
t189 = rSges(6,3) + pkin(12);
t182 = t105 * pkin(5);
t79 = -t135 * t183 - t147 * t186;
t180 = t79 * t99;
t152 = t171 * t183;
t143 = t100 * t152;
t179 = pkin(2) * t154 + qJ(3) * t143;
t153 = t171 * t184;
t178 = t187 * pkin(1) + pkin(10) * t153;
t177 = t100 * t99;
t169 = qJ(3) * t100;
t168 = t101 * t105;
t167 = t104 * t105;
t166 = t99 * t185;
t162 = t100 * t173;
t133 = t173 * t143;
t134 = t174 * t147;
t80 = -t134 * t183 + t148 * t186;
t161 = t80 * pkin(3) + pkin(11) * t133 + t179;
t160 = t100 * t166;
t155 = t187 * t171;
t159 = -pkin(1) * t184 + pkin(10) * t155;
t156 = t173 * t185;
t146 = -rSges(6,1) * t105 + rSges(6,2) * t102;
t63 = -t149 * t86 - t172 * t85;
t81 = t85 * pkin(2);
t145 = t63 * pkin(3) + t169 * t86 - t81;
t65 = -t149 * t87 - t121;
t83 = t126 * pkin(2);
t144 = t65 * pkin(3) + t169 * t87 - t83;
t142 = t104 * rSges(7,1) - t101 * rSges(7,2) + pkin(5);
t139 = t99 * t143;
t44 = t103 * t80 - t139 * t185 - t156 * t79;
t45 = t80 * t185 + (t173 * t79 + t139) * t103;
t138 = t45 * pkin(4) + t44 * pkin(12) + t161;
t62 = -t151 * t86 + t170 * t85;
t46 = t162 * t86 - t62 * t99;
t120 = t126 * t170;
t64 = -t151 * t87 + t120;
t47 = t162 * t87 - t64 * t99;
t27 = t63 * t103 - t156 * t62 - t160 * t86;
t28 = t63 * t185 + (t173 * t62 + t177 * t86) * t103;
t132 = t28 * pkin(4) + t27 * pkin(12) + t145;
t29 = t65 * t103 - t156 * t64 - t160 * t87;
t30 = t65 * t185 + (t173 * t64 + t177 * t87) * t103;
t131 = t30 * pkin(4) + t29 * pkin(12) + t144;
t130 = -t86 * pkin(2) + qJ(3) * t137 - t169 * t85 + t159;
t123 = -t54 * pkin(3) - pkin(11) * t38 + t130;
t122 = t22 * pkin(4) + t123;
t119 = t87 * pkin(2) + qJ(3) * t117 + t178;
t116 = (g(1) * t47 + g(2) * t46 - g(3) * t180) * pkin(11);
t57 = -t120 * t174 + t140 * t184 + t172 * t87;
t108 = t57 * pkin(3) + pkin(11) * t106 + t119;
t24 = -t103 * t192 + t57 * t185;
t107 = t24 * pkin(4) + t108;
t74 = t134 * t186 + t148 * t183 + t163 * t170;
t67 = t133 - t180;
t51 = t114 * t99 + t127 * t173;
t35 = -t103 * t191 + t74 * t185;
t34 = t103 * t74 + t185 * t191;
t33 = t34 * pkin(4);
t32 = t67 * t102 + t105 * t45;
t31 = t102 * t45 - t67 * t105;
t23 = t103 * t57 + t185 * t192;
t21 = t129 * t166 - t156 * t55 - t195;
t19 = t185 * t190 + t195;
t17 = t23 * pkin(4);
t15 = t19 * pkin(4);
t14 = t51 * t102 + t105 * t35;
t13 = -t35 * t102 + t105 * t51;
t12 = t47 * t102 + t105 * t30;
t11 = t102 * t30 - t47 * t105;
t10 = t46 * t102 + t105 * t28;
t9 = t102 * t28 - t46 * t105;
t8 = t102 * t106 + t24 * t105;
t7 = t102 * t24 - t105 * t106;
t2 = t101 * t23 + t104 * t8;
t1 = -t101 * t8 + t104 * t23;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t184 - rSges(2,2) * t187) + g(2) * (rSges(2,1) * t187 - rSges(2,2) * t184)) - m(3) * (g(1) * (-t86 * rSges(3,1) + t85 * rSges(3,2) + rSges(3,3) * t155 + t159) + g(2) * (t87 * rSges(3,1) - rSges(3,2) * t126 + rSges(3,3) * t153 + t178)) - m(4) * (g(1) * (-rSges(4,1) * t54 + t55 * rSges(4,2) - rSges(4,3) * t129 + t130) + g(2) * (t57 * rSges(4,1) - rSges(4,2) * t111 + rSges(4,3) * t117 + t119)) - m(5) * (g(1) * (t22 * rSges(5,1) - t21 * rSges(5,2) - rSges(5,3) * t38 + t123) + g(2) * (t24 * rSges(5,1) - t23 * rSges(5,2) + rSges(5,3) * t106 + t108)) - m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t196 + t189 * t21 + t122) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t189 * t23 + t107)) - m(7) * (g(1) * (-t4 * pkin(5) + t21 * pkin(12) + (t101 * t21 - t197) * rSges(7,1) + (t104 * t21 + t198) * rSges(7,2) + t188 * t196 + t122) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t8 * pkin(5) + t23 * pkin(12) + t188 * t7 + t107)) -m(3) * (g(1) * (-rSges(3,1) * t126 - t87 * rSges(3,2)) + g(2) * (-rSges(3,1) * t85 - rSges(3,2) * t86) + g(3) * (rSges(3,1) * t154 - rSges(3,2) * t152)) - m(4) * (g(1) * (rSges(4,1) * t65 + rSges(4,2) * t64 - t83) + g(2) * (rSges(4,1) * t63 + rSges(4,2) * t62 - t81) + g(3) * (t80 * rSges(4,1) + t79 * rSges(4,2) + t179) + (rSges(4,3) * g(3) * t152 + (g(1) * t87 + g(2) * t86) * (rSges(4,3) + qJ(3))) * t100) - m(5) * (g(1) * (t30 * rSges(5,1) - t29 * rSges(5,2) + t47 * rSges(5,3) + t144) + g(2) * (t28 * rSges(5,1) - t27 * rSges(5,2) + t46 * rSges(5,3) + t145) + g(3) * (rSges(5,1) * t45 - rSges(5,2) * t44 + rSges(5,3) * t67 + t161) + t116) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t29 * rSges(6,3) + t131) + g(2) * (t10 * rSges(6,1) - t9 * rSges(6,2) + t27 * rSges(6,3) + t132) + g(3) * (rSges(6,1) * t32 - rSges(6,2) * t31 + rSges(6,3) * t44 + t138) + t116) - m(7) * (g(1) * (t12 * pkin(5) + (t101 * t29 + t104 * t12) * rSges(7,1) + (-t101 * t12 + t104 * t29) * rSges(7,2) + t131 + t188 * t11) + g(2) * (t10 * pkin(5) + (t10 * t104 + t101 * t27) * rSges(7,1) + (-t10 * t101 + t104 * t27) * rSges(7,2) + t132 + t188 * t9) + g(3) * (t32 * pkin(5) + (t101 * t44 + t104 * t32) * rSges(7,1) + (-t101 * t32 + t104 * t44) * rSges(7,2) + t138 + t188 * t31) + t116) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t117 + g(2) * t129 + g(3) * t127) -m(5) * (g(1) * (-rSges(5,1) * t23 - rSges(5,2) * t24) + g(2) * (-rSges(5,1) * t19 + rSges(5,2) * t22) + g(3) * (-rSges(5,1) * t34 - rSges(5,2) * t35)) - m(6) * (g(1) * (t146 * t23 + t189 * t24 - t17) + g(2) * (t146 * t19 - t189 * t22 - t15) + g(3) * (t146 * t34 + t189 * t35 - t33)) + (-g(1) * (-t23 * t182 - t17 + t24 * pkin(12) + (t24 * t101 - t23 * t167) * rSges(7,1) + (t24 * t104 + t168 * t23) * rSges(7,2)) - g(2) * (-t19 * t182 - t15 - t22 * pkin(12) + (-t101 * t22 - t167 * t19) * rSges(7,1) + (-t104 * t22 + t168 * t19) * rSges(7,2)) - g(3) * (-t34 * t182 - t33 + t35 * pkin(12) + (t35 * t101 - t167 * t34) * rSges(7,1) + (t35 * t104 + t168 * t34) * rSges(7,2)) - (-g(1) * t23 - g(2) * t19 - g(3) * t34) * t102 * t188) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (rSges(6,1) * t196 - rSges(6,2) * t4) + g(3) * (rSges(6,1) * t13 - rSges(6,2) * t14)) - m(7) * (g(1) * (-t142 * t7 + t188 * t8) + (t142 * t13 + t188 * t14) * g(3) + (t142 * t196 + t188 * t4) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * ((t104 * t19 - t198) * rSges(7,1) + (-t101 * t19 - t197) * rSges(7,2)) + g(3) * ((-t101 * t14 + t104 * t34) * rSges(7,1) + (-t101 * t34 - t104 * t14) * rSges(7,2)))];
taug  = t3(:);
