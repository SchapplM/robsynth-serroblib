% Calculate joint inertia matrix for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:28
% EndTime: 2019-03-09 13:12:35
% DurationCPUTime: 3.04s
% Computational Cost: add. (10513->396), mult. (7597->581), div. (0->0), fcn. (7854->12), ass. (0->198)
t286 = Icges(3,3) + Icges(4,3);
t183 = qJ(2) + pkin(11);
t172 = sin(t183);
t173 = cos(t183);
t188 = sin(qJ(2));
t191 = cos(qJ(2));
t285 = Icges(3,5) * t191 + Icges(4,5) * t173 - Icges(3,6) * t188 - Icges(4,6) * t172;
t189 = sin(qJ(1));
t184 = t189 ^ 2;
t284 = t189 * pkin(7);
t192 = cos(qJ(1));
t174 = qJ(4) + t183;
t168 = sin(t174);
t169 = cos(t174);
t213 = Icges(5,5) * t169 - Icges(5,6) * t168;
t104 = -Icges(5,3) * t192 + t213 * t189;
t105 = Icges(5,3) * t189 + t213 * t192;
t185 = t192 ^ 2;
t254 = Icges(5,4) * t169;
t217 = -Icges(5,2) * t168 + t254;
t107 = Icges(5,6) * t189 + t217 * t192;
t255 = Icges(5,4) * t168;
t221 = Icges(5,1) * t169 - t255;
t109 = Icges(5,5) * t189 + t221 * t192;
t210 = -t107 * t168 + t109 * t169;
t106 = -Icges(5,6) * t192 + t217 * t189;
t108 = -Icges(5,5) * t192 + t221 * t189;
t211 = t106 * t168 - t108 * t169;
t170 = qJ(5) + t174;
t163 = sin(t170);
t164 = cos(t170);
t252 = Icges(6,4) * t164;
t216 = -Icges(6,2) * t163 + t252;
t95 = Icges(6,6) * t189 + t216 * t192;
t253 = Icges(6,4) * t163;
t220 = Icges(6,1) * t164 - t253;
t97 = Icges(6,5) * t189 + t220 * t192;
t224 = -t163 * t95 + t164 * t97;
t94 = -Icges(6,6) * t192 + t216 * t189;
t96 = -Icges(6,5) * t192 + t220 * t189;
t225 = t163 * t94 - t164 * t96;
t190 = cos(qJ(6));
t245 = t190 * t192;
t187 = sin(qJ(6));
t247 = t189 * t187;
t122 = -t164 * t247 - t245;
t246 = t189 * t190;
t248 = t187 * t192;
t123 = t164 * t246 - t248;
t251 = t163 * t189;
t62 = Icges(7,5) * t123 + Icges(7,6) * t122 + Icges(7,3) * t251;
t64 = Icges(7,4) * t123 + Icges(7,2) * t122 + Icges(7,6) * t251;
t66 = Icges(7,1) * t123 + Icges(7,4) * t122 + Icges(7,5) * t251;
t16 = t122 * t64 + t123 * t66 + t62 * t251;
t124 = -t164 * t248 + t246;
t125 = t164 * t245 + t247;
t250 = t163 * t192;
t63 = Icges(7,5) * t125 + Icges(7,6) * t124 + Icges(7,3) * t250;
t65 = Icges(7,4) * t125 + Icges(7,2) * t124 + Icges(7,6) * t250;
t67 = Icges(7,1) * t125 + Icges(7,4) * t124 + Icges(7,5) * t250;
t17 = t122 * t65 + t123 * t67 + t63 * t251;
t8 = -t16 * t192 + t17 * t189;
t212 = Icges(6,5) * t164 - Icges(6,6) * t163;
t92 = -Icges(6,3) * t192 + t212 * t189;
t93 = Icges(6,3) * t189 + t212 * t192;
t276 = -t185 * t92 - (t224 * t189 + (t225 - t93) * t192) * t189 - t8;
t283 = -t185 * t104 - (t210 * t189 + (-t105 + t211) * t192) * t189 + t276;
t249 = t164 * t192;
t69 = t125 * rSges(7,1) + t124 * rSges(7,2) + rSges(7,3) * t250;
t282 = pkin(5) * t249 + pkin(10) * t250 + t69;
t228 = rSges(5,1) * t169 - rSges(5,2) * t168;
t229 = rSges(4,1) * t173 - rSges(4,2) * t172;
t281 = -t285 * t189 + t286 * t192;
t280 = t286 * t189 + t285 * t192;
t279 = t189 / 0.2e1;
t278 = -t192 / 0.2e1;
t18 = t124 * t64 + t125 * t66 + t62 * t250;
t19 = t124 * t65 + t125 * t67 + t63 * t250;
t9 = -t18 * t192 + t19 * t189;
t277 = (t184 * t93 + t9 + (t225 * t192 + (t224 - t92) * t189) * t192) * t189;
t275 = pkin(2) * t188;
t274 = pkin(4) * t168;
t273 = pkin(5) * t164;
t186 = -qJ(3) - pkin(7);
t171 = t191 * pkin(2) + pkin(1);
t148 = pkin(3) * t173 + t171;
t121 = pkin(4) * t169 + t148;
t118 = t192 * t121;
t143 = t192 * t148;
t272 = t192 * (t118 - t143) + (t121 - t148) * t184;
t200 = rSges(6,1) * t249 - rSges(6,2) * t250 + t189 * rSges(6,3);
t227 = rSges(6,1) * t164 - rSges(6,2) * t163;
t52 = t189 * (-t192 * rSges(6,3) + t227 * t189) + t192 * t200;
t201 = t189 * rSges(5,3) + t228 * t192;
t55 = t189 * (-t192 * rSges(5,3) + t228 * t189) + t192 * t201;
t162 = t192 * t171;
t181 = t192 * pkin(7);
t271 = t189 * (t181 + (-pkin(1) + t171) * t189) + t192 * (-t192 * pkin(1) + t162 - t284);
t270 = rSges(3,1) * t191;
t267 = rSges(3,2) * t188;
t83 = -Icges(7,6) * t164 + (Icges(7,4) * t190 - Icges(7,2) * t187) * t163;
t264 = t187 * t83;
t263 = t192 * rSges(3,3);
t24 = -t164 * t62 + (-t187 * t64 + t190 * t66) * t163;
t262 = t24 * t192;
t25 = -t164 * t63 + (-t187 * t65 + t190 * t67) * t163;
t261 = t25 * t189;
t85 = -t164 * rSges(7,3) + (rSges(7,1) * t190 - rSges(7,2) * t187) * t163;
t260 = -pkin(5) * t163 + pkin(10) * t164 - t85;
t259 = Icges(3,4) * t188;
t258 = Icges(3,4) * t191;
t257 = Icges(4,4) * t172;
t256 = Icges(4,4) * t173;
t243 = t189 * rSges(3,3) + t192 * t270;
t240 = t184 + t185;
t182 = -pkin(8) + t186;
t239 = t189 * (t184 * t105 + (t211 * t192 + (-t104 + t210) * t189) * t192) + t277;
t238 = -rSges(4,1) * t172 - rSges(4,2) * t173 - t275;
t136 = rSges(6,1) * t163 + rSges(6,2) * t164;
t237 = -t136 - t274;
t27 = t52 + t272;
t236 = t192 * (t143 - t162) + t271 + (t148 - t171) * t184;
t175 = -pkin(9) + t182;
t235 = -t189 * t175 + t118;
t226 = -t123 * rSges(7,1) - t122 * rSges(7,2);
t68 = rSges(7,3) * t251 - t226;
t28 = t189 * t68 + t184 * (pkin(10) * t163 + t273) + t282 * t192;
t82 = -Icges(7,3) * t164 + (Icges(7,5) * t190 - Icges(7,6) * t187) * t163;
t84 = -Icges(7,5) * t164 + (Icges(7,1) * t190 - Icges(7,4) * t187) * t163;
t31 = t122 * t83 + t123 * t84 + t82 * t251;
t3 = -t31 * t164 + (t16 * t189 + t17 * t192) * t163;
t32 = t124 * t83 + t125 * t84 + t82 * t250;
t4 = -t32 * t164 + (t18 * t189 + t19 * t192) * t163;
t234 = t3 * t278 + t4 * t279 - t164 * (t261 - t262) / 0.2e1 + t8 * t251 / 0.2e1 + t9 * t250 / 0.2e1;
t233 = t260 - t274;
t231 = -pkin(3) * t172 - t275;
t230 = -t267 + t270;
t223 = Icges(3,1) * t191 - t259;
t222 = Icges(4,1) * t173 - t257;
t219 = -Icges(3,2) * t188 + t258;
t218 = -Icges(4,2) * t172 + t256;
t134 = Icges(6,2) * t164 + t253;
t135 = Icges(6,1) * t163 + t252;
t205 = -t134 * t163 + t135 * t164;
t140 = Icges(5,2) * t169 + t255;
t141 = Icges(5,1) * t168 + t254;
t204 = -t140 * t168 + t141 * t169;
t203 = t276 * t192 + t277;
t13 = t28 + t272;
t202 = t189 * rSges(4,3) + t229 * t192;
t142 = rSges(5,1) * t168 + rSges(5,2) * t169;
t199 = -t142 + t231;
t133 = Icges(6,5) * t163 + Icges(6,6) * t164;
t198 = -t262 / 0.2e1 + t261 / 0.2e1 + (t189 * t133 + t163 * t97 + t164 * t95 + t205 * t192 + t32) * t279 + (-t192 * t133 + t163 * t96 + t164 * t94 + t205 * t189 + t31) * t278;
t197 = t231 - t274;
t196 = t283 * t192 + t239;
t195 = -t136 + t197;
t194 = t197 + t260;
t139 = Icges(5,5) * t168 + Icges(5,6) * t169;
t193 = t198 + (t107 * t169 + t109 * t168 + t189 * t139 + t204 * t192) * t279 + (t106 * t169 + t108 * t168 - t192 * t139 + t204 * t189) * t278;
t160 = rSges(2,1) * t192 - t189 * rSges(2,2);
t159 = -t189 * rSges(2,1) - rSges(2,2) * t192;
t158 = rSges(3,1) * t188 + rSges(3,2) * t191;
t111 = t238 * t192;
t110 = t238 * t189;
t101 = t284 + (pkin(1) - t267) * t192 + t243;
t100 = t263 + t181 + (-pkin(1) - t230) * t189;
t87 = t237 * t192;
t86 = t237 * t189;
t81 = t199 * t192;
t80 = t199 * t189;
t79 = -t189 * t186 + t162 + t202;
t78 = (rSges(4,3) - t186) * t192 + (-t171 - t229) * t189;
t75 = t163 * t190 * t84;
t74 = t192 * (-t192 * t267 + t243) + (t230 * t189 - t263) * t189;
t73 = -t189 * t182 + t143 + t201;
t72 = (rSges(5,3) - t182) * t192 + (-t148 - t228) * t189;
t71 = t195 * t192;
t70 = t195 * t189;
t57 = t260 * t192;
t56 = t260 * t189;
t54 = t200 + t235;
t53 = (rSges(6,3) - t175) * t192 + (-t121 - t227) * t189;
t51 = t233 * t192;
t50 = t233 * t189;
t43 = t194 * t192;
t42 = t194 * t189;
t39 = -t164 * t69 - t85 * t250;
t38 = t164 * t68 + t85 * t251;
t37 = t192 * t202 + (-t192 * rSges(4,3) + t229 * t189) * t189 + t271;
t36 = t235 + t282;
t35 = -t192 * t175 + (-t273 - t121 + (-rSges(7,3) - pkin(10)) * t163) * t189 + t226;
t34 = -t163 * t264 - t164 * t82 + t75;
t33 = (-t189 * t69 + t192 * t68) * t163;
t26 = t236 + t55;
t12 = t236 + t27;
t10 = t13 + t236;
t1 = [t169 * t140 + t168 * t141 + t173 * (Icges(4,2) * t173 + t257) + t172 * (Icges(4,1) * t172 + t256) + t191 * (Icges(3,2) * t191 + t259) + t188 * (Icges(3,1) * t188 + t258) + Icges(2,3) + t75 + (-t82 + t134) * t164 + (t135 - t264) * t163 + m(7) * (t35 ^ 2 + t36 ^ 2) + m(6) * (t53 ^ 2 + t54 ^ 2) + m(5) * (t72 ^ 2 + t73 ^ 2) + m(4) * (t78 ^ 2 + t79 ^ 2) + m(3) * (t100 ^ 2 + t101 ^ 2) + m(2) * (t159 ^ 2 + t160 ^ 2); t193 + m(3) * (-t100 * t192 - t101 * t189) * t158 + m(4) * (t110 * t79 + t111 * t78) + m(5) * (t72 * t81 + t73 * t80) + m(6) * (t53 * t71 + t54 * t70) + m(7) * (t35 * t43 + t36 * t42) + (t173 * (Icges(4,6) * t189 + t218 * t192) + t172 * (Icges(4,5) * t189 + t222 * t192) + t191 * (Icges(3,6) * t189 + t219 * t192) + t188 * (Icges(3,5) * t189 + t223 * t192)) * t279 + (t173 * (-Icges(4,6) * t192 + t218 * t189) + t172 * (-Icges(4,5) * t192 + t222 * t189) + t191 * (-Icges(3,6) * t192 + t219 * t189) + t188 * (-Icges(3,5) * t192 + t223 * t189)) * t278 + (Icges(3,5) * t188 + Icges(4,5) * t172 + Icges(3,6) * t191 + Icges(4,6) * t173) * (t185 / 0.2e1 + t184 / 0.2e1); m(7) * (t10 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t12 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t26 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t110 ^ 2 + t111 ^ 2 + t37 ^ 2) + m(3) * (t240 * t158 ^ 2 + t74 ^ 2) + t239 + t280 * t189 * t184 + (t281 * t185 + (t189 * t281 + t192 * t280) * t189 + t283) * t192; m(7) * (t189 * t35 - t192 * t36) + m(6) * (t189 * t53 - t192 * t54) + m(5) * (t189 * t72 - t192 * t73) + m(4) * (t189 * t78 - t192 * t79); m(7) * (t189 * t43 - t192 * t42) + m(6) * (t189 * t71 - t192 * t70) + m(5) * (t189 * t81 - t192 * t80) + m(4) * (-t110 * t192 + t189 * t111); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t240; m(7) * (t35 * t51 + t36 * t50) + m(6) * (t53 * t87 + t54 * t86) + m(5) * (-t189 * t73 - t192 * t72) * t142 + t193; m(7) * (t10 * t13 + t42 * t50 + t43 * t51) + m(6) * (t12 * t27 + t70 * t86 + t71 * t87) + m(5) * (t55 * t26 + (-t189 * t80 - t192 * t81) * t142) + t196; m(6) * (t87 * t189 - t192 * t86) + m(7) * (t51 * t189 - t192 * t50); m(7) * (t13 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t27 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(5) * (t142 ^ 2 * t240 + t55 ^ 2) + t196; m(7) * (t35 * t57 + t36 * t56) + m(6) * (-t189 * t54 - t192 * t53) * t136 + t198; m(7) * (t10 * t28 + t42 * t56 + t43 * t57) + m(6) * (t52 * t12 + (-t189 * t70 - t192 * t71) * t136) + t203; m(7) * (t57 * t189 - t192 * t56); m(7) * (t13 * t28 + t50 * t56 + t51 * t57) + m(6) * (t52 * t27 + (-t189 * t86 - t192 * t87) * t136) + t203; m(6) * (t136 ^ 2 * t240 + t52 ^ 2) + m(7) * (t28 ^ 2 + t56 ^ 2 + t57 ^ 2) + t203; m(7) * (t35 * t38 + t36 * t39) - t34 * t164 + ((t25 / 0.2e1 + t32 / 0.2e1) * t192 + (t24 / 0.2e1 + t31 / 0.2e1) * t189) * t163; m(7) * (t10 * t33 + t38 * t43 + t39 * t42) + t234; m(7) * (t38 * t189 - t192 * t39); m(7) * (t13 * t33 + t38 * t51 + t39 * t50) + t234; m(7) * (t28 * t33 + t38 * t57 + t39 * t56) + t234; t164 ^ 2 * t34 + m(7) * (t33 ^ 2 + t38 ^ 2 + t39 ^ 2) + (t192 * t4 + t189 * t3 - t164 * (t189 * t24 + t192 * t25)) * t163;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
