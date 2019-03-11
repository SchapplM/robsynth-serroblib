% Calculate joint inertia matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:09
% EndTime: 2019-03-08 21:04:16
% DurationCPUTime: 3.72s
% Computational Cost: add. (16240->541), mult. (28506->793), div. (0->0), fcn. (35503->12), ass. (0->246)
t212 = sin(pkin(6));
t271 = t212 ^ 2;
t241 = m(6) / 0.2e1 + m(7) / 0.2e1;
t270 = 0.2e1 * t241;
t211 = sin(pkin(10));
t213 = cos(pkin(10));
t221 = cos(qJ(2));
t214 = cos(pkin(6));
t218 = sin(qJ(2));
t255 = t214 * t218;
t200 = t211 * t221 + t213 * t255;
t242 = qJ(3) + pkin(11);
t233 = sin(t242);
t225 = t212 * t233;
t234 = cos(t242);
t183 = t200 * t234 - t213 * t225;
t202 = -t211 * t255 + t213 * t221;
t185 = t202 * t234 + t211 * t225;
t226 = t212 * t234;
t197 = t214 * t233 + t218 * t226;
t184 = t202 * t233 - t211 * t226;
t254 = t214 * t221;
t201 = t211 * t254 + t213 * t218;
t216 = sin(qJ(6));
t219 = cos(qJ(6));
t150 = t184 * t219 - t201 * t216;
t151 = t184 * t216 + t201 * t219;
t182 = t200 * t233 + t213 * t226;
t199 = t211 * t218 - t213 * t254;
t148 = t182 * t219 - t199 * t216;
t149 = t182 * t216 + t199 * t219;
t89 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t183;
t91 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t183;
t93 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t183;
t37 = t150 * t91 + t151 * t93 + t185 * t89;
t90 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t185;
t92 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t185;
t94 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t185;
t38 = t150 * t92 + t151 * t94 + t185 * t90;
t196 = -t214 * t234 + t218 * t225;
t257 = t212 * t221;
t186 = t196 * t219 + t216 * t257;
t187 = t196 * t216 - t219 * t257;
t105 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t197;
t106 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t197;
t107 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t197;
t48 = t105 * t185 + t106 * t150 + t107 * t151;
t2 = t183 * t37 + t185 * t38 + t197 * t48;
t269 = t2 / 0.2e1;
t268 = t183 / 0.2e1;
t267 = t185 / 0.2e1;
t266 = t197 / 0.2e1;
t220 = cos(qJ(3));
t265 = pkin(3) * t220;
t95 = rSges(7,1) * t149 + rSges(7,2) * t148 + rSges(7,3) * t183;
t263 = pkin(5) * t199 + pkin(9) * t183 + t95;
t96 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t185;
t262 = pkin(5) * t201 + pkin(9) * t185 + t96;
t261 = t211 * t212;
t260 = t212 * t213;
t217 = sin(qJ(3));
t259 = t212 * t217;
t258 = t212 * t220;
t256 = t214 * t217;
t239 = t213 * t259;
t128 = -pkin(3) * t239 + qJ(4) * t199 + t265 * t200;
t100 = t201 * t128;
t144 = pkin(4) * t183 + qJ(5) * t182;
t253 = t201 * t144 + t100;
t174 = pkin(3) * t256 + (-qJ(4) * t221 + t265 * t218) * t212;
t252 = t128 * t257 + t199 * t174;
t108 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t197;
t251 = -pkin(5) * t257 + t197 * pkin(9) + t108;
t240 = t211 * t259;
t129 = pkin(3) * t240 + qJ(4) * t201 + t265 * t202;
t181 = t202 * pkin(2) + t201 * pkin(8);
t179 = t214 * t181;
t250 = t214 * t129 + t179;
t126 = rSges(5,1) * t185 - rSges(5,2) * t184 + rSges(5,3) * t201;
t249 = -t126 - t129;
t180 = t200 * pkin(2) + t199 * pkin(8);
t248 = -t128 - t180;
t145 = pkin(4) * t185 + qJ(5) * t184;
t247 = -t129 - t145;
t161 = t197 * rSges(5,1) - t196 * rSges(5,2) - rSges(5,3) * t257;
t246 = -t161 - t174;
t162 = pkin(4) * t197 + qJ(5) * t196;
t245 = -t162 - t174;
t244 = t180 * t261 + t181 * t260;
t243 = -m(5) - m(6) - m(7);
t238 = t214 * t145 + t250;
t124 = rSges(6,1) * t201 - rSges(6,2) * t185 + rSges(6,3) * t184;
t237 = -t124 + t247;
t236 = -t144 + t248;
t160 = -rSges(6,1) * t257 - t197 * rSges(6,2) + t196 * rSges(6,3);
t235 = -t160 + t245;
t203 = t214 * t220 - t218 * t259;
t204 = t218 * t258 + t256;
t175 = t204 * rSges(4,1) + t203 * rSges(4,2) - rSges(4,3) * t257;
t205 = (pkin(2) * t218 - pkin(8) * t221) * t212;
t232 = (-t175 - t205) * t212;
t231 = t247 - t262;
t230 = t128 * t261 + t129 * t260 + t244;
t229 = t144 * t257 + t199 * t162 + t252;
t228 = t245 - t251;
t227 = (-t205 + t246) * t212;
t224 = (-t205 + t235) * t212;
t223 = t144 * t261 + t145 * t260 + t230;
t222 = (-t205 + t228) * t212;
t198 = t214 * rSges(3,3) + (rSges(3,1) * t218 + rSges(3,2) * t221) * t212;
t195 = Icges(3,5) * t214 + (Icges(3,1) * t218 + Icges(3,4) * t221) * t212;
t194 = Icges(3,6) * t214 + (Icges(3,4) * t218 + Icges(3,2) * t221) * t212;
t193 = Icges(3,3) * t214 + (Icges(3,5) * t218 + Icges(3,6) * t221) * t212;
t192 = t202 * t220 + t240;
t191 = -t202 * t217 + t211 * t258;
t190 = t200 * t220 - t239;
t189 = -t200 * t217 - t213 * t258;
t173 = Icges(4,1) * t204 + Icges(4,4) * t203 - Icges(4,5) * t257;
t172 = Icges(4,4) * t204 + Icges(4,2) * t203 - Icges(4,6) * t257;
t171 = Icges(4,5) * t204 + Icges(4,6) * t203 - Icges(4,3) * t257;
t170 = rSges(3,1) * t202 - rSges(3,2) * t201 + rSges(3,3) * t261;
t169 = rSges(3,1) * t200 - rSges(3,2) * t199 - rSges(3,3) * t260;
t168 = Icges(3,1) * t202 - Icges(3,4) * t201 + Icges(3,5) * t261;
t167 = Icges(3,1) * t200 - Icges(3,4) * t199 - Icges(3,5) * t260;
t166 = Icges(3,4) * t202 - Icges(3,2) * t201 + Icges(3,6) * t261;
t165 = Icges(3,4) * t200 - Icges(3,2) * t199 - Icges(3,6) * t260;
t164 = Icges(3,5) * t202 - Icges(3,6) * t201 + Icges(3,3) * t261;
t163 = Icges(3,5) * t200 - Icges(3,6) * t199 - Icges(3,3) * t260;
t159 = Icges(5,1) * t197 - Icges(5,4) * t196 - Icges(5,5) * t257;
t158 = Icges(5,4) * t197 - Icges(5,2) * t196 - Icges(5,6) * t257;
t157 = Icges(5,5) * t197 - Icges(5,6) * t196 - Icges(5,3) * t257;
t156 = -Icges(6,1) * t257 - Icges(6,4) * t197 + Icges(6,5) * t196;
t155 = -Icges(6,4) * t257 - Icges(6,2) * t197 + Icges(6,6) * t196;
t154 = -Icges(6,5) * t257 - Icges(6,6) * t197 + Icges(6,3) * t196;
t143 = -t169 * t214 - t198 * t260;
t142 = t170 * t214 - t198 * t261;
t137 = rSges(4,1) * t192 + rSges(4,2) * t191 + rSges(4,3) * t201;
t136 = rSges(4,1) * t190 + rSges(4,2) * t189 + rSges(4,3) * t199;
t135 = Icges(4,1) * t192 + Icges(4,4) * t191 + Icges(4,5) * t201;
t134 = Icges(4,1) * t190 + Icges(4,4) * t189 + Icges(4,5) * t199;
t133 = Icges(4,4) * t192 + Icges(4,2) * t191 + Icges(4,6) * t201;
t132 = Icges(4,4) * t190 + Icges(4,2) * t189 + Icges(4,6) * t199;
t131 = Icges(4,5) * t192 + Icges(4,6) * t191 + Icges(4,3) * t201;
t130 = Icges(4,5) * t190 + Icges(4,6) * t189 + Icges(4,3) * t199;
t125 = rSges(5,1) * t183 - rSges(5,2) * t182 + rSges(5,3) * t199;
t123 = rSges(6,1) * t199 - rSges(6,2) * t183 + rSges(6,3) * t182;
t122 = Icges(5,1) * t185 - Icges(5,4) * t184 + Icges(5,5) * t201;
t121 = Icges(5,1) * t183 - Icges(5,4) * t182 + Icges(5,5) * t199;
t120 = Icges(6,1) * t201 - Icges(6,4) * t185 + Icges(6,5) * t184;
t119 = Icges(6,1) * t199 - Icges(6,4) * t183 + Icges(6,5) * t182;
t118 = Icges(5,4) * t185 - Icges(5,2) * t184 + Icges(5,6) * t201;
t117 = Icges(5,4) * t183 - Icges(5,2) * t182 + Icges(5,6) * t199;
t116 = Icges(6,4) * t201 - Icges(6,2) * t185 + Icges(6,6) * t184;
t115 = Icges(6,4) * t199 - Icges(6,2) * t183 + Icges(6,6) * t182;
t114 = Icges(5,5) * t185 - Icges(5,6) * t184 + Icges(5,3) * t201;
t113 = Icges(5,5) * t183 - Icges(5,6) * t182 + Icges(5,3) * t199;
t112 = Icges(6,5) * t201 - Icges(6,6) * t185 + Icges(6,3) * t184;
t111 = Icges(6,5) * t199 - Icges(6,6) * t183 + Icges(6,3) * t182;
t101 = (t169 * t211 + t170 * t213) * t212;
t98 = -t137 * t257 - t201 * t175;
t97 = t136 * t257 + t199 * t175;
t88 = -t171 * t257 + t203 * t172 + t204 * t173;
t87 = t136 * t201 - t137 * t199;
t86 = (-t136 - t180) * t214 + t213 * t232;
t85 = t214 * t137 + t211 * t232 + t179;
t84 = -t157 * t257 - t196 * t158 + t197 * t159;
t83 = t196 * t154 - t197 * t155 - t156 * t257;
t82 = t171 * t201 + t172 * t191 + t173 * t192;
t81 = t171 * t199 + t172 * t189 + t173 * t190;
t80 = t157 * t201 - t158 * t184 + t159 * t185;
t79 = t157 * t199 - t158 * t182 + t159 * t183;
t78 = t154 * t184 - t155 * t185 + t156 * t201;
t77 = t154 * t182 - t155 * t183 + t156 * t199;
t76 = (t136 * t211 + t137 * t213) * t212 + t244;
t75 = -t131 * t257 + t203 * t133 + t204 * t135;
t74 = -t130 * t257 + t203 * t132 + t204 * t134;
t73 = -t108 * t185 + t197 * t96;
t72 = t108 * t183 - t197 * t95;
t71 = -t114 * t257 - t196 * t118 + t197 * t122;
t70 = -t113 * t257 - t196 * t117 + t197 * t121;
t69 = t196 * t112 - t197 * t116 - t120 * t257;
t68 = t196 * t111 - t197 * t115 - t119 * t257;
t67 = t246 * t201 + t249 * t257;
t66 = t125 * t257 + t199 * t161 + t252;
t65 = t131 * t201 + t133 * t191 + t135 * t192;
t64 = t130 * t201 + t132 * t191 + t134 * t192;
t63 = t131 * t199 + t133 * t189 + t135 * t190;
t62 = t130 * t199 + t132 * t189 + t134 * t190;
t61 = (-t125 + t248) * t214 + t213 * t227;
t60 = t214 * t126 + t211 * t227 + t250;
t59 = t114 * t201 - t118 * t184 + t122 * t185;
t58 = t113 * t201 - t117 * t184 + t121 * t185;
t57 = t114 * t199 - t118 * t182 + t122 * t183;
t56 = t113 * t199 - t117 * t182 + t121 * t183;
t55 = t112 * t184 - t116 * t185 + t120 * t201;
t54 = t111 * t184 - t115 * t185 + t119 * t201;
t53 = t112 * t182 - t116 * t183 + t120 * t199;
t52 = t111 * t182 - t115 * t183 + t119 * t199;
t51 = t105 * t197 + t106 * t186 + t107 * t187;
t50 = -t183 * t96 + t185 * t95;
t49 = t201 * t125 + t249 * t199 + t100;
t47 = t105 * t183 + t106 * t148 + t107 * t149;
t46 = t235 * t201 + t237 * t257;
t45 = t123 * t257 + t199 * t160 + t229;
t44 = (t125 * t211 + t126 * t213) * t212 + t230;
t43 = (-t123 + t236) * t214 + t213 * t224;
t42 = t214 * t124 + t211 * t224 + t238;
t41 = t186 * t92 + t187 * t94 + t197 * t90;
t40 = t186 * t91 + t187 * t93 + t197 * t89;
t39 = t201 * t123 + t237 * t199 + t253;
t36 = t148 * t92 + t149 * t94 + t183 * t90;
t35 = t148 * t91 + t149 * t93 + t183 * t89;
t34 = (t123 * t211 + t124 * t213) * t212 + t223;
t33 = t228 * t201 + t231 * t257;
t32 = t251 * t199 + t263 * t257 + t229;
t31 = (t236 - t263) * t214 + t213 * t222;
t30 = t211 * t222 + t262 * t214 + t238;
t29 = t88 * t214 + (t211 * t75 - t213 * t74) * t212;
t28 = t74 * t199 + t75 * t201 - t88 * t257;
t27 = t84 * t214 + (t211 * t71 - t213 * t70) * t212;
t26 = t83 * t214 + (t211 * t69 - t213 * t68) * t212;
t25 = t82 * t214 + (t211 * t65 - t213 * t64) * t212;
t24 = t81 * t214 + (t211 * t63 - t213 * t62) * t212;
t23 = t231 * t199 + t263 * t201 + t253;
t22 = (t263 * t211 + t262 * t213) * t212 + t223;
t21 = t70 * t199 + t71 * t201 - t84 * t257;
t20 = t68 * t199 + t69 * t201 - t83 * t257;
t19 = t64 * t199 + t65 * t201 - t82 * t257;
t18 = t62 * t199 + t63 * t201 - t81 * t257;
t17 = t80 * t214 + (t211 * t59 - t213 * t58) * t212;
t16 = t79 * t214 + (t211 * t57 - t213 * t56) * t212;
t15 = t78 * t214 + (t211 * t55 - t213 * t54) * t212;
t14 = t77 * t214 + (t211 * t53 - t213 * t52) * t212;
t13 = t58 * t199 + t59 * t201 - t80 * t257;
t12 = t56 * t199 + t57 * t201 - t79 * t257;
t11 = t54 * t199 + t55 * t201 - t78 * t257;
t10 = t52 * t199 + t53 * t201 - t77 * t257;
t9 = t51 * t214 + (t211 * t41 - t213 * t40) * t212;
t8 = t40 * t199 + t41 * t201 - t51 * t257;
t7 = t183 * t40 + t185 * t41 + t197 * t51;
t6 = t48 * t214 + (t211 * t38 - t213 * t37) * t212;
t5 = t47 * t214 + (t211 * t36 - t213 * t35) * t212;
t4 = t37 * t199 + t38 * t201 - t48 * t257;
t3 = t35 * t199 + t36 * t201 - t47 * t257;
t1 = t183 * t35 + t185 * t36 + t197 * t47;
t99 = [m(2) + m(3) + m(4) - t243; m(3) * t101 + m(4) * t76 + m(5) * t44 + m(6) * t34 + m(7) * t22; m(7) * (t22 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t34 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t44 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t76 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(3) * (t101 ^ 2 + t142 ^ 2 + t143 ^ 2) + (t6 + t25 + t15 + t17 + (t164 * t261 - t166 * t201 + t168 * t202) * t261) * t261 + (-t5 - t24 - t16 - t14 + (-t163 * t260 - t165 * t199 + t167 * t200) * t260 + (-t163 * t261 + t164 * t260 + t165 * t201 + t166 * t199 - t167 * t202 - t168 * t200) * t261) * t260 + (-(-t193 * t260 - t199 * t194 + t200 * t195) * t260 + (t193 * t261 - t201 * t194 + t202 * t195) * t261 + t9 + t27 + t26 + t29 + ((t166 * t221 + t168 * t218) * t211 - (t165 * t221 + t167 * t218) * t213) * t271 + ((-t163 * t213 + t164 * t211 + t194 * t221 + t195 * t218) * t212 + t214 * t193) * t214) * t214; m(4) * t87 + m(5) * t49 + m(6) * t39 + m(7) * t23; (t21 / 0.2e1 + t20 / 0.2e1 + t28 / 0.2e1 + t8 / 0.2e1) * t214 + (t17 / 0.2e1 + t15 / 0.2e1 + t25 / 0.2e1 + t6 / 0.2e1) * t201 + (t24 / 0.2e1 + t14 / 0.2e1 + t16 / 0.2e1 + t5 / 0.2e1) * t199 + m(7) * (t22 * t23 + t30 * t33 + t31 * t32) + m(6) * (t34 * t39 + t42 * t46 + t43 * t45) + m(5) * (t44 * t49 + t60 * t67 + t61 * t66) + m(4) * (t76 * t87 + t85 * t98 + t86 * t97) + ((-t26 / 0.2e1 - t27 / 0.2e1 - t29 / 0.2e1 - t9 / 0.2e1) * t221 + (-t3 / 0.2e1 - t18 / 0.2e1 - t10 / 0.2e1 - t12 / 0.2e1) * t213 + (t19 / 0.2e1 + t4 / 0.2e1 + t11 / 0.2e1 + t13 / 0.2e1) * t211) * t212; (-t20 - t21 - t28 - t8) * t257 + (t13 + t19 + t11 + t4) * t201 + (t10 + t12 + t18 + t3) * t199 + m(4) * (t87 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(7) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t49 ^ 2 + t66 ^ 2 + t67 ^ 2); t243 * t257; m(7) * (t199 * t30 + t201 * t31 - t22 * t257) + m(6) * (t199 * t42 + t201 * t43 - t34 * t257) + m(5) * (t199 * t60 + t201 * t61 - t44 * t257); m(7) * (t199 * t33 + t201 * t32 - t23 * t257) + m(6) * (t199 * t46 + t201 * t45 - t39 * t257) + m(5) * (t199 * t67 + t201 * t66 - t49 * t257); 0.2e1 * (m(5) / 0.2e1 + t241) * (t271 * t221 ^ 2 + t199 ^ 2 + t201 ^ 2); t196 * t270; m(7) * (t182 * t30 + t184 * t31 + t196 * t22) + m(6) * (t182 * t42 + t184 * t43 + t196 * t34); m(7) * (t182 * t33 + t184 * t32 + t196 * t23) + m(6) * (t182 * t46 + t184 * t45 + t196 * t39); (t182 * t199 + t184 * t201 - t196 * t257) * t270; (t182 ^ 2 + t184 ^ 2 + t196 ^ 2) * t270; m(7) * t50; t5 * t268 + t214 * t7 / 0.2e1 + t9 * t266 + t6 * t267 + m(7) * (t22 * t50 + t30 * t73 + t31 * t72) + (t211 * t269 - t213 * t1 / 0.2e1) * t212; m(7) * (t23 * t50 + t32 * t72 + t33 * t73) + t8 * t266 + t201 * t269 + t3 * t268 + t4 * t267 + t199 * t1 / 0.2e1 - t7 * t257 / 0.2e1; m(7) * (t73 * t199 + t72 * t201 - t50 * t257); m(7) * (t182 * t73 + t184 * t72 + t196 * t50); t185 * t2 + t183 * t1 + t197 * t7 + m(7) * (t50 ^ 2 + t72 ^ 2 + t73 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t99(1) t99(2) t99(4) t99(7) t99(11) t99(16); t99(2) t99(3) t99(5) t99(8) t99(12) t99(17); t99(4) t99(5) t99(6) t99(9) t99(13) t99(18); t99(7) t99(8) t99(9) t99(10) t99(14) t99(19); t99(11) t99(12) t99(13) t99(14) t99(15) t99(20); t99(16) t99(17) t99(18) t99(19) t99(20) t99(21);];
Mq  = res;
