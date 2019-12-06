% Calculate joint inertia matrix for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:19
% EndTime: 2019-12-05 17:23:32
% DurationCPUTime: 4.54s
% Computational Cost: add. (32112->470), mult. (89102->695), div. (0->0), fcn. (117424->14), ass. (0->220)
t184 = sin(pkin(11));
t186 = cos(pkin(11));
t191 = sin(qJ(2));
t187 = cos(pkin(5));
t193 = cos(qJ(2));
t218 = t193 * t187;
t177 = -t184 * t191 + t186 * t218;
t219 = t191 * t187;
t178 = t193 * t184 + t186 * t219;
t190 = sin(qJ(3));
t185 = sin(pkin(5));
t223 = sin(pkin(6));
t201 = t185 * t223;
t224 = cos(pkin(6));
t229 = cos(qJ(3));
t146 = t178 * t229 + (t177 * t224 - t186 * t201) * t190;
t202 = t185 * t224;
t166 = -t177 * t223 - t186 * t202;
t189 = sin(qJ(4));
t228 = cos(qJ(4));
t132 = t146 * t189 - t166 * t228;
t179 = -t184 * t218 - t186 * t191;
t180 = -t184 * t219 + t193 * t186;
t148 = t180 * t229 + (t179 * t224 + t184 * t201) * t190;
t167 = -t179 * t223 + t184 * t202;
t134 = t148 * t189 - t167 * t228;
t165 = t187 * t223 * t190 + (t190 * t193 * t224 + t191 * t229) * t185;
t176 = t187 * t224 - t193 * t201;
t149 = t165 * t189 - t176 * t228;
t133 = t146 * t228 + t166 * t189;
t196 = t229 * t223;
t194 = t185 * t196;
t197 = t224 * t229;
t145 = -t177 * t197 + t178 * t190 + t186 * t194;
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t104 = -t133 * t188 + t145 * t192;
t105 = t133 * t192 + t145 * t188;
t74 = Icges(6,5) * t105 + Icges(6,6) * t104 + Icges(6,3) * t132;
t76 = Icges(6,4) * t105 + Icges(6,2) * t104 + Icges(6,6) * t132;
t78 = Icges(6,1) * t105 + Icges(6,4) * t104 + Icges(6,5) * t132;
t24 = t104 * t76 + t105 * t78 + t132 * t74;
t135 = t148 * t228 + t167 * t189;
t147 = -t179 * t197 + t180 * t190 - t184 * t194;
t106 = -t135 * t188 + t147 * t192;
t107 = t135 * t192 + t147 * t188;
t75 = Icges(6,5) * t107 + Icges(6,6) * t106 + Icges(6,3) * t134;
t77 = Icges(6,4) * t107 + Icges(6,2) * t106 + Icges(6,6) * t134;
t79 = Icges(6,1) * t107 + Icges(6,4) * t106 + Icges(6,5) * t134;
t25 = t104 * t77 + t105 * t79 + t132 * t75;
t150 = t165 * t228 + t176 * t189;
t221 = t185 * t191;
t164 = -t185 * t193 * t197 - t187 * t196 + t190 * t221;
t130 = -t150 * t188 + t164 * t192;
t131 = t150 * t192 + t164 * t188;
t90 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t149;
t93 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t149;
t96 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t149;
t41 = t104 * t93 + t105 * t96 + t132 * t90;
t1 = t132 * t24 + t134 * t25 + t149 * t41;
t235 = t1 / 0.2e1;
t26 = t106 * t76 + t107 * t78 + t134 * t74;
t27 = t106 * t77 + t107 * t79 + t134 * t75;
t42 = t106 * t93 + t107 * t96 + t134 * t90;
t2 = t132 * t26 + t134 * t27 + t149 * t42;
t234 = t2 / 0.2e1;
t32 = t130 * t76 + t131 * t78 + t149 * t74;
t33 = t130 * t77 + t131 * t79 + t149 * t75;
t50 = t130 * t93 + t131 * t96 + t149 * t90;
t9 = t132 * t32 + t134 * t33 + t149 * t50;
t233 = t9 / 0.2e1;
t232 = t132 / 0.2e1;
t231 = t134 / 0.2e1;
t230 = t149 / 0.2e1;
t80 = rSges(6,1) * t105 + rSges(6,2) * t104 + rSges(6,3) * t132;
t227 = pkin(4) * t133 + pkin(10) * t132 + t80;
t81 = rSges(6,1) * t107 + rSges(6,2) * t106 + rSges(6,3) * t134;
t226 = pkin(4) * t135 + pkin(10) * t134 + t81;
t125 = pkin(3) * t146 + pkin(9) * t145;
t99 = rSges(5,1) * t133 - rSges(5,2) * t132 + rSges(5,3) * t145;
t225 = -t125 - t99;
t222 = t184 * t185;
t220 = t186 * t185;
t101 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t149;
t217 = pkin(4) * t150 + pkin(10) * t149 + t101;
t119 = rSges(5,1) * t150 - rSges(5,2) * t149 + rSges(5,3) * t164;
t142 = pkin(3) * t165 + pkin(9) * t164;
t216 = -t119 - t142;
t126 = pkin(3) * t148 + pkin(9) * t147;
t153 = t180 * pkin(2) + pkin(8) * t167;
t151 = t187 * t153;
t215 = t187 * t126 + t151;
t152 = t178 * pkin(2) + pkin(8) * t166;
t214 = t152 * t222 + t153 * t220;
t91 = Icges(5,5) * t133 - Icges(5,6) * t132 + Icges(5,3) * t145;
t94 = Icges(5,4) * t133 - Icges(5,2) * t132 + Icges(5,6) * t145;
t97 = Icges(5,1) * t133 - Icges(5,4) * t132 + Icges(5,5) * t145;
t46 = -t132 * t94 + t133 * t97 + t145 * t91;
t92 = Icges(5,5) * t135 - Icges(5,6) * t134 + Icges(5,3) * t147;
t95 = Icges(5,4) * t135 - Icges(5,2) * t134 + Icges(5,6) * t147;
t98 = Icges(5,1) * t135 - Icges(5,4) * t134 + Icges(5,5) * t147;
t47 = -t132 * t95 + t133 * t98 + t145 * t92;
t114 = Icges(5,5) * t150 - Icges(5,6) * t149 + Icges(5,3) * t164;
t115 = Icges(5,4) * t150 - Icges(5,2) * t149 + Icges(5,6) * t164;
t116 = Icges(5,1) * t150 - Icges(5,4) * t149 + Icges(5,5) * t164;
t57 = t114 * t145 - t115 * t132 + t116 * t133;
t13 = t145 * t46 + t147 * t47 + t164 * t57;
t3 = t145 * t24 + t147 * t25 + t164 * t41;
t213 = t3 / 0.2e1 + t13 / 0.2e1;
t48 = -t134 * t94 + t135 * t97 + t147 * t91;
t49 = -t134 * t95 + t135 * t98 + t147 * t92;
t58 = t114 * t147 - t115 * t134 + t116 * t135;
t14 = t145 * t48 + t147 * t49 + t164 * t58;
t4 = t145 * t26 + t147 * t27 + t164 * t42;
t212 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t166 * t46 + t167 * t47 + t176 * t57;
t5 = t166 * t24 + t167 * t25 + t176 * t41;
t211 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t166 * t48 + t167 * t49 + t176 * t58;
t6 = t166 * t26 + t167 * t27 + t176 * t42;
t210 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t187 * t57 + (t184 * t47 - t186 * t46) * t185;
t7 = t187 * t41 + (t184 * t25 - t186 * t24) * t185;
t209 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t187 * t58 + (t184 * t49 - t186 * t48) * t185;
t8 = t187 * t42 + (t184 * t27 - t186 * t26) * t185;
t208 = t8 / 0.2e1 + t18 / 0.2e1;
t10 = t145 * t32 + t147 * t33 + t164 * t50;
t51 = -t149 * t94 + t150 * t97 + t164 * t91;
t52 = -t149 * t95 + t150 * t98 + t164 * t92;
t68 = t114 * t164 - t115 * t149 + t116 * t150;
t19 = t145 * t51 + t147 * t52 + t164 * t68;
t207 = t10 / 0.2e1 + t19 / 0.2e1;
t11 = t166 * t32 + t167 * t33 + t176 * t50;
t20 = t166 * t51 + t167 * t52 + t176 * t68;
t206 = t11 / 0.2e1 + t20 / 0.2e1;
t12 = t187 * t50 + (t184 * t33 - t186 * t32) * t185;
t21 = t187 * t68 + (t184 * t52 - t186 * t51) * t185;
t205 = t12 / 0.2e1 + t21 / 0.2e1;
t204 = -t125 - t227;
t203 = -t142 - t217;
t139 = rSges(4,1) * t165 - rSges(4,2) * t164 + rSges(4,3) * t176;
t168 = pkin(2) * t221 + pkin(8) * t176;
t200 = (-t139 - t168) * t185;
t199 = t125 * t222 + t126 * t220 + t214;
t198 = (-t168 + t216) * t185;
t195 = (-t168 + t203) * t185;
t174 = t187 * rSges(3,3) + (rSges(3,1) * t191 + rSges(3,2) * t193) * t185;
t173 = Icges(3,5) * t187 + (Icges(3,1) * t191 + Icges(3,4) * t193) * t185;
t172 = Icges(3,6) * t187 + (Icges(3,4) * t191 + Icges(3,2) * t193) * t185;
t171 = Icges(3,3) * t187 + (Icges(3,5) * t191 + Icges(3,6) * t193) * t185;
t161 = rSges(3,1) * t180 + rSges(3,2) * t179 + rSges(3,3) * t222;
t160 = rSges(3,1) * t178 + rSges(3,2) * t177 - rSges(3,3) * t220;
t159 = Icges(3,1) * t180 + Icges(3,4) * t179 + Icges(3,5) * t222;
t158 = Icges(3,1) * t178 + Icges(3,4) * t177 - Icges(3,5) * t220;
t157 = Icges(3,4) * t180 + Icges(3,2) * t179 + Icges(3,6) * t222;
t156 = Icges(3,4) * t178 + Icges(3,2) * t177 - Icges(3,6) * t220;
t155 = Icges(3,5) * t180 + Icges(3,6) * t179 + Icges(3,3) * t222;
t154 = Icges(3,5) * t178 + Icges(3,6) * t177 - Icges(3,3) * t220;
t141 = -t160 * t187 - t174 * t220;
t140 = t161 * t187 - t174 * t222;
t138 = Icges(4,1) * t165 - Icges(4,4) * t164 + Icges(4,5) * t176;
t137 = Icges(4,4) * t165 - Icges(4,2) * t164 + Icges(4,6) * t176;
t136 = Icges(4,5) * t165 - Icges(4,6) * t164 + Icges(4,3) * t176;
t129 = t166 * t142;
t128 = (t160 * t184 + t161 * t186) * t185;
t121 = t176 * t126;
t120 = t167 * t125;
t118 = rSges(4,1) * t148 - rSges(4,2) * t147 + rSges(4,3) * t167;
t117 = rSges(4,1) * t146 - rSges(4,2) * t145 + rSges(4,3) * t166;
t113 = Icges(4,1) * t148 - Icges(4,4) * t147 + Icges(4,5) * t167;
t112 = Icges(4,1) * t146 - Icges(4,4) * t145 + Icges(4,5) * t166;
t111 = Icges(4,4) * t148 - Icges(4,2) * t147 + Icges(4,6) * t167;
t110 = Icges(4,4) * t146 - Icges(4,2) * t145 + Icges(4,6) * t166;
t109 = Icges(4,5) * t148 - Icges(4,6) * t147 + Icges(4,3) * t167;
t108 = Icges(4,5) * t146 - Icges(4,6) * t145 + Icges(4,3) * t166;
t100 = rSges(5,1) * t135 - rSges(5,2) * t134 + rSges(5,3) * t147;
t89 = t118 * t176 - t139 * t167;
t88 = -t117 * t176 + t139 * t166;
t87 = (-t117 - t152) * t187 + t186 * t200;
t86 = t118 * t187 + t184 * t200 + t151;
t85 = t136 * t176 - t137 * t164 + t138 * t165;
t84 = t117 * t167 - t118 * t166;
t83 = t136 * t167 - t137 * t147 + t138 * t148;
t82 = t136 * t166 - t137 * t145 + t138 * t146;
t73 = (t117 * t184 + t118 * t186) * t185 + t214;
t72 = t100 * t164 - t119 * t147;
t71 = t119 * t145 - t164 * t99;
t70 = t109 * t176 - t111 * t164 + t113 * t165;
t69 = t108 * t176 - t110 * t164 + t112 * t165;
t67 = t109 * t167 - t111 * t147 + t113 * t148;
t66 = t108 * t167 - t110 * t147 + t112 * t148;
t65 = t109 * t166 - t111 * t145 + t113 * t146;
t64 = t108 * t166 - t110 * t145 + t112 * t146;
t63 = -t100 * t145 + t147 * t99;
t62 = t100 * t176 + t167 * t216 + t121;
t61 = t119 * t166 + t176 * t225 + t129;
t60 = (-t152 + t225) * t187 + t186 * t198;
t59 = t100 * t187 + t184 * t198 + t215;
t56 = -t101 * t134 + t149 * t81;
t55 = t101 * t132 - t149 * t80;
t54 = t167 * t99 + t120 + (-t100 - t126) * t166;
t53 = (t100 * t186 + t184 * t99) * t185 + t199;
t45 = -t132 * t81 + t134 * t80;
t44 = -t147 * t217 + t164 * t226;
t43 = t145 * t217 - t164 * t227;
t40 = (-t152 + t204) * t187 + t186 * t195;
t39 = t184 * t195 + t187 * t226 + t215;
t38 = t167 * t203 + t176 * t226 + t121;
t37 = t166 * t217 + t176 * t204 + t129;
t36 = t187 * t85 + (t184 * t70 - t186 * t69) * t185;
t35 = -t145 * t226 + t147 * t227;
t34 = t166 * t69 + t167 * t70 + t176 * t85;
t31 = (t184 * t227 + t186 * t226) * t185 + t199;
t30 = t187 * t83 + (t184 * t67 - t186 * t66) * t185;
t29 = t187 * t82 + (t184 * t65 - t186 * t64) * t185;
t28 = t120 + t227 * t167 + (-t126 - t226) * t166;
t23 = t166 * t66 + t167 * t67 + t176 * t83;
t22 = t166 * t64 + t167 * t65 + t176 * t82;
t102 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t128 + m(4) * t73 + m(5) * t53 + m(6) * t31; m(6) * (t31 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t53 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t73 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(3) * (t128 ^ 2 + t140 ^ 2 + t141 ^ 2) + (t8 + t18 + t30 + (t155 * t222 + t157 * t179 + t159 * t180) * t222) * t222 + (-t7 - t17 - t29 + (-t154 * t220 + t156 * t177 + t158 * t178) * t220 + (-t154 * t222 + t155 * t220 - t156 * t179 - t157 * t177 - t158 * t180 - t159 * t178) * t222) * t220 + (t12 + t21 + t36 + ((t157 * t193 + t159 * t191) * t184 - (t156 * t193 + t158 * t191) * t186) * t185 ^ 2 - (-t171 * t220 + t172 * t177 + t173 * t178) * t220 + (t171 * t222 + t172 * t179 + t173 * t180) * t222 + ((-t154 * t186 + t155 * t184 + t172 * t193 + t173 * t191) * t185 + t187 * t171) * t187) * t187; m(4) * t84 + m(5) * t54 + m(6) * t28; (t34 / 0.2e1 + t206) * t187 + (t36 / 0.2e1 + t205) * t176 + (t30 / 0.2e1 + t208) * t167 + (t29 / 0.2e1 + t209) * t166 + m(6) * (t28 * t31 + t37 * t40 + t38 * t39) + m(5) * (t53 * t54 + t59 * t62 + t60 * t61) + m(4) * (t73 * t84 + t86 * t89 + t87 * t88) + ((-t22 / 0.2e1 - t211) * t186 + (t23 / 0.2e1 + t210) * t184) * t185; (t11 + t20 + t34) * t176 + (t6 + t16 + t23) * t167 + (t5 + t15 + t22) * t166 + m(6) * (t28 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t54 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(4) * (t84 ^ 2 + t88 ^ 2 + t89 ^ 2); m(5) * t63 + m(6) * t35; t207 * t187 + t205 * t164 + t208 * t147 + t209 * t145 + m(6) * (t31 * t35 + t39 * t44 + t40 * t43) + m(5) * (t53 * t63 + t59 * t72 + t60 * t71) + (t184 * t212 - t186 * t213) * t185; t207 * t176 + t212 * t167 + t213 * t166 + t206 * t164 + t210 * t147 + t211 * t145 + m(6) * (t28 * t35 + t37 * t43 + t38 * t44) + m(5) * (t54 * t63 + t61 * t71 + t62 * t72); (t10 + t19) * t164 + (t4 + t14) * t147 + (t3 + t13) * t145 + m(6) * (t35 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t63 ^ 2 + t71 ^ 2 + t72 ^ 2); m(6) * t45; t12 * t230 + m(6) * (t31 * t45 + t39 * t56 + t40 * t55) + t8 * t231 + t7 * t232 + t187 * t233 + (t184 * t234 - t186 * t1 / 0.2e1) * t185; m(6) * (t28 * t45 + t37 * t55 + t38 * t56) + t6 * t231 + t11 * t230 + t5 * t232 + t167 * t234 + t176 * t233 + t166 * t235; t164 * t233 + t3 * t232 + t4 * t231 + m(6) * (t35 * t45 + t43 * t55 + t44 * t56) + t147 * t234 + t145 * t235 + t10 * t230; m(6) * (t45 ^ 2 + t55 ^ 2 + t56 ^ 2) + t134 * t2 + t132 * t1 + t149 * t9;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t102(1), t102(2), t102(4), t102(7), t102(11); t102(2), t102(3), t102(5), t102(8), t102(12); t102(4), t102(5), t102(6), t102(9), t102(13); t102(7), t102(8), t102(9), t102(10), t102(14); t102(11), t102(12), t102(13), t102(14), t102(15);];
Mq = res;
