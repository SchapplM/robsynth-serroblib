% Calculate Gravitation load on the joints for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:45
% EndTime: 2019-03-10 05:49:57
% DurationCPUTime: 4.30s
% Computational Cost: add. (3118->343), mult. (8861->521), div. (0->0), fcn. (11549->18), ass. (0->150)
t122 = sin(qJ(6));
t123 = sin(qJ(5));
t127 = cos(qJ(5));
t124 = sin(qJ(4));
t211 = cos(qJ(4));
t120 = sin(pkin(8));
t201 = cos(pkin(8));
t203 = cos(pkin(6));
t209 = sin(qJ(2));
t179 = t203 * t209;
t210 = sin(qJ(1));
t212 = cos(qJ(2));
t213 = cos(qJ(1));
t106 = t179 * t213 + t210 * t212;
t125 = sin(qJ(3));
t128 = cos(qJ(3));
t180 = t203 * t212;
t105 = -t180 * t213 + t209 * t210;
t121 = sin(pkin(7));
t200 = sin(pkin(6));
t177 = t213 * t200;
t202 = cos(pkin(7));
t222 = t202 * t105 + t121 * t177;
t218 = t106 * t125 + t128 * t222;
t173 = t202 * t200;
t223 = -t105 * t121 + t213 * t173;
t219 = t120 * t223 + t201 * t218;
t72 = t106 * t128 - t222 * t125;
t24 = -t124 * t219 + t211 * t72;
t51 = -t218 * t120 + t201 * t223;
t6 = t123 * t51 - t127 * t24;
t229 = t122 * t6;
t126 = cos(qJ(6));
t228 = t126 * t6;
t227 = -t123 * t24 - t51 * t127;
t224 = t124 * t72;
t107 = -t179 * t210 + t212 * t213;
t149 = t180 * t210 + t209 * t213;
t175 = t200 * t210;
t139 = t121 * t175 - t149 * t202;
t134 = t107 * t125 - t128 * t139;
t140 = t149 * t121 + t210 * t173;
t129 = t134 * t120 + t140 * t201;
t214 = rSges(7,3) + pkin(14);
t221 = -t120 * t140 + t134 * t201;
t152 = t121 * t203 + t173 * t212;
t174 = t200 * t209;
t141 = t125 * t174 - t128 * t152;
t176 = t212 * t200;
t151 = -t121 * t176 + t202 * t203;
t220 = -t151 * t120 + t141 * t201;
t215 = rSges(6,3) + pkin(13);
t208 = pkin(11) * t121;
t207 = pkin(12) * t120;
t206 = t127 * pkin(5);
t161 = t209 * t173;
t99 = -t125 * t176 - t128 * t161;
t204 = t99 * t120;
t197 = t120 * t121;
t196 = t120 * t123;
t195 = t120 * t127;
t194 = t122 * t127;
t193 = t126 * t127;
t166 = t121 * t174;
t192 = pkin(2) * t176 + pkin(11) * t166;
t191 = t213 * pkin(1) + pkin(10) * t175;
t190 = t120 * t211;
t188 = t121 * t201;
t187 = t124 * t201;
t186 = t125 * t202;
t185 = t128 * t202;
t100 = -t125 * t161 + t128 * t176;
t158 = t201 * t166;
t183 = t100 * pkin(3) + pkin(12) * t158 + t192;
t182 = t121 * t190;
t181 = -pkin(1) * t210 + pkin(10) * t177;
t178 = t201 * t211;
t172 = -rSges(6,1) * t127 + rSges(6,2) * t123;
t38 = -t187 * t72 - t211 * t218;
t68 = t218 * pkin(3);
t171 = t38 * pkin(4) + t207 * t72 - t68;
t75 = t107 * t128 + t125 * t139;
t40 = -t134 * t211 - t187 * t75;
t70 = t134 * pkin(3);
t170 = t40 * pkin(4) + t207 * t75 - t70;
t94 = t125 * t152 + t128 * t174;
t54 = -t141 * t211 - t187 * t94;
t93 = t141 * pkin(3);
t169 = t54 * pkin(4) + t207 * t94 - t93;
t101 = t105 * pkin(2);
t81 = -t105 * t128 - t106 * t186;
t168 = t81 * pkin(3) + t106 * t208 - t101;
t103 = t149 * pkin(2);
t83 = -t107 * t186 - t128 * t149;
t167 = t83 * pkin(3) + t107 * t208 - t103;
t164 = t126 * rSges(7,1) - t122 * rSges(7,2) + pkin(5);
t163 = t120 * t166;
t58 = t100 * t124 - t163 * t211 - t178 * t99;
t59 = t100 * t211 + (t201 * t99 + t163) * t124;
t162 = t59 * pkin(4) + t58 * pkin(13) + t183;
t80 = t105 * t125 - t106 * t185;
t60 = t106 * t188 - t80 * t120;
t82 = -t107 * t185 + t125 * t149;
t61 = t107 * t188 - t82 * t120;
t33 = -t106 * t182 + t81 * t124 - t178 * t80;
t34 = t81 * t211 + (t106 * t197 + t201 * t80) * t124;
t157 = t34 * pkin(4) + t33 * pkin(13) + t168;
t35 = -t107 * t182 + t83 * t124 - t178 * t82;
t36 = t83 * t211 + (t107 * t197 + t201 * t82) * t124;
t156 = t36 * pkin(4) + t35 * pkin(13) + t167;
t155 = -t106 * pkin(2) + t223 * pkin(11) + t181;
t146 = -t72 * pkin(3) + t51 * pkin(12) + t155;
t144 = -pkin(4) * t24 + t146;
t142 = t107 * pkin(2) + t140 * pkin(11) + t191;
t137 = (g(1) * t61 + g(2) * t60 - g(3) * t204) * pkin(12);
t131 = t75 * pkin(3) + t129 * pkin(12) + t142;
t27 = t124 * t75 + t211 * t221;
t28 = -t124 * t221 + t211 * t75;
t130 = t28 * pkin(4) + t27 * pkin(13) + t131;
t85 = t158 - t204;
t67 = t120 * t141 + t151 * t201;
t53 = -t124 * t141 + t178 * t94;
t47 = -t124 * t220 + t94 * t211;
t46 = t124 * t94 + t211 * t220;
t45 = t46 * pkin(4);
t44 = t127 * t54 + t196 * t94;
t43 = t123 * t54 - t195 * t94;
t42 = t123 * t85 + t127 * t59;
t41 = t123 * t59 - t85 * t127;
t39 = -t124 * t134 + t178 * t75;
t37 = -t124 * t218 + t178 * t72;
t25 = -t178 * t218 - t190 * t223 - t224;
t23 = t211 * t219 + t224;
t21 = t27 * pkin(4);
t19 = t23 * pkin(4);
t18 = t123 * t67 + t127 * t47;
t17 = -t123 * t47 + t127 * t67;
t16 = t127 * t40 + t196 * t75;
t15 = t123 * t40 - t195 * t75;
t14 = t127 * t38 + t196 * t72;
t13 = t123 * t38 - t195 * t72;
t12 = t123 * t61 + t127 * t36;
t11 = t123 * t36 - t61 * t127;
t10 = t123 * t60 + t127 * t34;
t9 = t123 * t34 - t60 * t127;
t8 = t123 * t129 + t28 * t127;
t7 = t123 * t28 - t127 * t129;
t2 = t122 * t27 + t126 * t8;
t1 = -t122 * t8 + t126 * t27;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t210 - rSges(2,2) * t213) + g(2) * (rSges(2,1) * t213 - rSges(2,2) * t210)) - m(3) * (g(1) * (-t106 * rSges(3,1) + t105 * rSges(3,2) + rSges(3,3) * t177 + t181) + g(2) * (t107 * rSges(3,1) - rSges(3,2) * t149 + rSges(3,3) * t175 + t191)) - m(4) * (g(1) * (-rSges(4,1) * t72 + rSges(4,2) * t218 + rSges(4,3) * t223 + t155) + g(2) * (t75 * rSges(4,1) - rSges(4,2) * t134 + rSges(4,3) * t140 + t142)) - m(5) * (g(1) * (-rSges(5,1) * t24 - t25 * rSges(5,2) + t51 * rSges(5,3) + t146) + g(2) * (t28 * rSges(5,1) - t27 * rSges(5,2) + rSges(5,3) * t129 + t131)) - m(6) * (g(1) * (t6 * rSges(6,1) - rSges(6,2) * t227 + t215 * t25 + t144) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t27 * rSges(6,3) + t130)) - m(7) * (g(1) * (t6 * pkin(5) + t25 * pkin(13) + (t122 * t25 + t228) * rSges(7,1) + (t126 * t25 - t229) * rSges(7,2) + t214 * t227 + t144) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t8 * pkin(5) + t214 * t7 + t130)) -m(3) * (g(1) * (-rSges(3,1) * t149 - t107 * rSges(3,2)) + g(2) * (-rSges(3,1) * t105 - rSges(3,2) * t106) + g(3) * (rSges(3,1) * t176 - rSges(3,2) * t174)) - m(4) * (g(1) * (rSges(4,1) * t83 + rSges(4,2) * t82 - t103) + g(2) * (rSges(4,1) * t81 + rSges(4,2) * t80 - t101) + g(3) * (t100 * rSges(4,1) + t99 * rSges(4,2) + t192) + (g(3) * rSges(4,3) * t174 + (g(1) * t107 + g(2) * t106) * (rSges(4,3) + pkin(11))) * t121) - m(5) * (g(1) * (t36 * rSges(5,1) - t35 * rSges(5,2) + t61 * rSges(5,3) + t167) + g(2) * (t34 * rSges(5,1) - t33 * rSges(5,2) + t60 * rSges(5,3) + t168) + g(3) * (rSges(5,1) * t59 - rSges(5,2) * t58 + rSges(5,3) * t85 + t183) + t137) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t35 * rSges(6,3) + t156) + g(2) * (t10 * rSges(6,1) - t9 * rSges(6,2) + t33 * rSges(6,3) + t157) + g(3) * (rSges(6,1) * t42 - rSges(6,2) * t41 + rSges(6,3) * t58 + t162) + t137) - m(7) * (g(1) * (t12 * pkin(5) + (t12 * t126 + t122 * t35) * rSges(7,1) + (-t12 * t122 + t126 * t35) * rSges(7,2) + t156 + t214 * t11) + g(2) * (t10 * pkin(5) + (t10 * t126 + t122 * t33) * rSges(7,1) + (-t10 * t122 + t126 * t33) * rSges(7,2) + t157 + t214 * t9) + g(3) * (t42 * pkin(5) + (t122 * t58 + t126 * t42) * rSges(7,1) + (-t122 * t42 + t126 * t58) * rSges(7,2) + t162 + t214 * t41) + t137) -m(4) * (g(1) * (-rSges(4,1) * t134 - t75 * rSges(4,2)) + g(2) * (-rSges(4,1) * t218 - t72 * rSges(4,2)) + g(3) * (-rSges(4,1) * t141 - t94 * rSges(4,2))) - m(6) * (g(1) * (rSges(6,1) * t16 - rSges(6,2) * t15 + t215 * t39 + t170) + g(2) * (rSges(6,1) * t14 - rSges(6,2) * t13 + t215 * t37 + t171) + g(3) * (rSges(6,1) * t44 - rSges(6,2) * t43 + t215 * t53 + t169)) - m(7) * (g(1) * (t16 * pkin(5) + t39 * pkin(13) + (t122 * t39 + t126 * t16) * rSges(7,1) + (-t122 * t16 + t126 * t39) * rSges(7,2) + t214 * t15 + t170) + g(2) * (t14 * pkin(5) + t37 * pkin(13) + (t122 * t37 + t126 * t14) * rSges(7,1) + (-t122 * t14 + t126 * t37) * rSges(7,2) + t214 * t13 + t171) + g(3) * (t44 * pkin(5) + t53 * pkin(13) + (t122 * t53 + t126 * t44) * rSges(7,1) + (-t122 * t44 + t126 * t53) * rSges(7,2) + t214 * t43 + t169)) + (-g(1) * (rSges(5,1) * t40 - rSges(5,2) * t39 - t70) - g(2) * (rSges(5,1) * t38 - rSges(5,2) * t37 - t68) - g(3) * (rSges(5,1) * t54 - rSges(5,2) * t53 - t93) - (g(1) * t75 + g(2) * t72 + g(3) * t94) * t120 * (rSges(5,3) + pkin(12))) * m(5), -m(5) * (g(1) * (-rSges(5,1) * t27 - rSges(5,2) * t28) + g(2) * (-rSges(5,1) * t23 - rSges(5,2) * t24) + g(3) * (-rSges(5,1) * t46 - rSges(5,2) * t47)) - m(6) * (g(1) * (t172 * t27 + t215 * t28 - t21) + g(2) * (t172 * t23 + t215 * t24 - t19) + g(3) * (t172 * t46 + t215 * t47 - t45)) + (-g(1) * (-t27 * t206 - t21 + t28 * pkin(13) + (t122 * t28 - t193 * t27) * rSges(7,1) + (t126 * t28 + t194 * t27) * rSges(7,2)) - g(2) * (-t23 * t206 - t19 + t24 * pkin(13) + (t122 * t24 - t193 * t23) * rSges(7,1) + (t126 * t24 + t194 * t23) * rSges(7,2)) - g(3) * (-t46 * t206 - t45 + t47 * pkin(13) + (t122 * t47 - t193 * t46) * rSges(7,1) + (t126 * t47 + t194 * t46) * rSges(7,2)) - (-g(1) * t27 - g(2) * t23 - g(3) * t46) * t123 * t214) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (rSges(6,1) * t227 + rSges(6,2) * t6) + g(3) * (rSges(6,1) * t17 - rSges(6,2) * t18)) - m(7) * (g(1) * (-t164 * t7 + t214 * t8) + (t164 * t17 + t214 * t18) * g(3) + (t164 * t227 - t214 * t6) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * ((t126 * t23 + t229) * rSges(7,1) + (-t122 * t23 + t228) * rSges(7,2)) + g(3) * ((-t122 * t18 + t126 * t46) * rSges(7,1) + (-t122 * t46 - t126 * t18) * rSges(7,2)))];
taug  = t3(:);
