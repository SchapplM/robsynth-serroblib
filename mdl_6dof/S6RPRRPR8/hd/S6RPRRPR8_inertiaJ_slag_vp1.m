% Calculate joint inertia matrix for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:42
% EndTime: 2019-03-09 05:22:49
% DurationCPUTime: 3.59s
% Computational Cost: add. (7754->441), mult. (9488->622), div. (0->0), fcn. (10117->10), ass. (0->225)
t213 = cos(qJ(3));
t301 = -0.2e1 * t213;
t300 = Icges(4,5) * t213;
t210 = sin(qJ(3));
t299 = Icges(4,6) * t210;
t298 = t300 / 0.2e1 - t299 / 0.2e1;
t204 = qJ(4) + pkin(10);
t194 = sin(t204);
t195 = cos(t204);
t128 = Icges(6,3) * t210 + (Icges(6,5) * t195 - Icges(6,6) * t194) * t213;
t130 = Icges(6,5) * t210 + (Icges(6,1) * t195 - Icges(6,4) * t194) * t213;
t209 = sin(qJ(4));
t212 = cos(qJ(4));
t145 = Icges(5,3) * t210 + (Icges(5,5) * t212 - Icges(5,6) * t209) * t213;
t151 = Icges(5,5) * t210 + (Icges(5,1) * t212 - Icges(5,4) * t209) * t213;
t297 = (t130 * t195 + t151 * t212) * t213 + (t128 + t145) * t210;
t129 = Icges(6,6) * t210 + (Icges(6,4) * t195 - Icges(6,2) * t194) * t213;
t148 = Icges(5,6) * t210 + (Icges(5,4) * t212 - Icges(5,2) * t209) * t213;
t296 = -t129 * t194 - t148 * t209;
t211 = sin(qJ(1));
t214 = cos(qJ(1));
t262 = t211 * t194;
t141 = t195 * t214 - t210 * t262;
t261 = t211 * t195;
t142 = t194 * t214 + t210 * t261;
t258 = t211 * t213;
t87 = Icges(6,5) * t142 + Icges(6,6) * t141 - Icges(6,3) * t258;
t89 = Icges(6,4) * t142 + Icges(6,2) * t141 - Icges(6,6) * t258;
t91 = Icges(6,1) * t142 + Icges(6,4) * t141 - Icges(6,5) * t258;
t30 = t141 * t89 + t142 * t91 - t258 * t87;
t265 = t210 * t214;
t143 = t194 * t265 + t261;
t144 = -t195 * t265 + t262;
t256 = t213 * t214;
t88 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t256;
t90 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t256;
t92 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t256;
t31 = t141 * t90 + t142 * t92 - t258 * t88;
t257 = t212 * t214;
t260 = t211 * t209;
t162 = -t210 * t260 + t257;
t259 = t211 * t212;
t267 = t209 * t214;
t163 = t210 * t259 + t267;
t100 = Icges(5,4) * t163 + Icges(5,2) * t162 - Icges(5,6) * t258;
t102 = Icges(5,1) * t163 + Icges(5,4) * t162 - Icges(5,5) * t258;
t98 = Icges(5,5) * t163 + Icges(5,6) * t162 - Icges(5,3) * t258;
t40 = t100 * t162 + t102 * t163 - t258 * t98;
t164 = t209 * t265 + t259;
t165 = -t210 * t257 + t260;
t101 = Icges(5,4) * t165 + Icges(5,2) * t164 + Icges(5,6) * t256;
t103 = Icges(5,1) * t165 + Icges(5,4) * t164 + Icges(5,5) * t256;
t99 = Icges(5,5) * t165 + Icges(5,6) * t164 + Icges(5,3) * t256;
t41 = t101 * t162 + t103 * t163 - t258 * t99;
t54 = -t128 * t258 + t129 * t141 + t130 * t142;
t61 = -t145 * t258 + t148 * t162 + t151 * t163;
t295 = ((t31 + t41) * t214 + (-t30 - t40) * t211) * t213 + (t61 + t54) * t210;
t32 = t143 * t89 + t144 * t91 + t256 * t87;
t33 = t143 * t90 + t144 * t92 + t256 * t88;
t42 = t100 * t164 + t102 * t165 + t256 * t98;
t43 = t101 * t164 + t103 * t165 + t256 * t99;
t55 = t128 * t256 + t129 * t143 + t130 * t144;
t62 = t145 * t256 + t148 * t164 + t151 * t165;
t294 = ((t33 + t43) * t214 + (-t32 - t42) * t211) * t213 + (t62 + t55) * t210;
t198 = t214 * qJ(2);
t193 = pkin(4) * t212 + pkin(3);
t208 = -qJ(5) - pkin(8);
t248 = t193 * t265 + t208 * t256;
t280 = pkin(4) * t209;
t288 = -pkin(1) - pkin(7);
t94 = t144 * rSges(6,1) + t143 * rSges(6,2) + rSges(6,3) * t256;
t67 = t198 + (-t280 + t288) * t211 - t94 + t248;
t247 = pkin(1) * t214 + qJ(2) * t211;
t239 = pkin(7) * t214 + t247;
t266 = t210 * t211;
t241 = pkin(4) * t267 + t193 * t266 + t208 * t258;
t93 = rSges(6,1) * t142 + rSges(6,2) * t141 - rSges(6,3) * t258;
t68 = t93 + t239 + t241;
t293 = m(6) * (t211 * t67 - t214 * t68);
t173 = pkin(5) * t194 + t280;
t203 = -pkin(9) + t208;
t196 = qJ(6) + t204;
t191 = sin(t196);
t192 = cos(t196);
t263 = t211 * t192;
t137 = t191 * t265 + t263;
t264 = t211 * t191;
t138 = -t192 * t265 + t264;
t229 = -t138 * rSges(7,1) - t137 * rSges(7,2);
t170 = pkin(5) * t195 + t193;
t268 = t170 * t210;
t58 = t198 + (t268 + (-rSges(7,3) + t203) * t213) * t214 + (-t173 + t288) * t211 + t229;
t242 = t170 * t266 + t173 * t214 + t203 * t258;
t135 = t192 * t214 - t210 * t264;
t136 = t191 * t214 + t210 * t263;
t85 = rSges(7,1) * t136 + rSges(7,2) * t135 - rSges(7,3) * t258;
t59 = t85 + t239 + t242;
t292 = m(7) * (t211 * t58 - t214 * t59);
t38 = t210 * t87 + (-t194 * t89 + t195 * t91) * t213;
t46 = t210 * t98 + (-t100 * t209 + t102 * t212) * t213;
t291 = t38 + t46;
t39 = t210 * t88 + (-t194 * t90 + t195 * t92) * t213;
t47 = t210 * t99 + (-t101 * t209 + t103 * t212) * t213;
t290 = t39 + t47;
t289 = (rSges(4,1) * t210 + rSges(4,2) * t213) * t214;
t205 = t211 ^ 2;
t207 = t214 ^ 2;
t79 = Icges(7,5) * t136 + Icges(7,6) * t135 - Icges(7,3) * t258;
t81 = Icges(7,4) * t136 + Icges(7,2) * t135 - Icges(7,6) * t258;
t83 = Icges(7,1) * t136 + Icges(7,4) * t135 - Icges(7,5) * t258;
t36 = t210 * t79 + (-t191 * t81 + t192 * t83) * t213;
t80 = Icges(7,5) * t138 + Icges(7,6) * t137 + Icges(7,3) * t256;
t82 = Icges(7,4) * t138 + Icges(7,2) * t137 + Icges(7,6) * t256;
t84 = Icges(7,1) * t138 + Icges(7,4) * t137 + Icges(7,5) * t256;
t37 = t210 * t80 + (-t191 * t82 + t192 * t84) * t213;
t27 = t137 * t81 + t138 * t83 + t256 * t79;
t28 = t137 * t82 + t138 * t84 + t256 * t80;
t122 = Icges(7,3) * t210 + (Icges(7,5) * t192 - Icges(7,6) * t191) * t213;
t123 = Icges(7,6) * t210 + (Icges(7,4) * t192 - Icges(7,2) * t191) * t213;
t124 = Icges(7,5) * t210 + (Icges(7,1) * t192 - Icges(7,4) * t191) * t213;
t51 = t122 * t256 + t123 * t137 + t124 * t138;
t5 = t51 * t210 + (-t211 * t27 + t214 * t28) * t213;
t254 = t124 * t192 * t213 + t122 * t210;
t271 = t123 * t191;
t60 = (-t213 * t271 + t254) * t210;
t287 = t5 * t256 + t210 * (t60 + (-t211 * t36 + t214 * t37) * t213);
t285 = t210 / 0.2e1;
t284 = t211 / 0.2e1;
t282 = t214 / 0.2e1;
t180 = rSges(4,1) * t213 - rSges(4,2) * t210;
t281 = m(4) * t180;
t279 = (t213 * t296 + t297) * t210;
t86 = rSges(7,3) * t256 - t229;
t278 = (-t203 * t213 - t268) * t214 + (t173 - t280) * t211 + t248 + t86;
t188 = pkin(3) * t266;
t167 = -pkin(8) * t258 + t188;
t104 = -t167 + t241;
t277 = -t104 - t93;
t190 = pkin(3) * t265;
t233 = pkin(8) * t256 - t190;
t105 = pkin(4) * t260 - t233 - t248;
t276 = -t105 - t94;
t125 = t210 * rSges(7,3) + (rSges(7,1) * t192 - rSges(7,2) * t191) * t213;
t65 = t125 * t258 + t210 * t85;
t127 = (-pkin(3) + t193) * t213 + (-pkin(8) - t208) * t210;
t275 = t104 * t210 + t127 * t258;
t158 = t214 * t233;
t274 = t105 * t214 + t158;
t108 = (t170 - t193) * t213 + (-t203 + t208) * t210;
t255 = t108 + t125;
t182 = t213 * pkin(3) + t210 * pkin(8);
t169 = t211 * t182;
t252 = t127 * t211 + t169;
t250 = -t127 - t182;
t249 = rSges(5,1) * t163 + rSges(5,2) * t162;
t246 = t205 + t207;
t245 = m(6) / 0.2e1 + m(7) / 0.2e1;
t72 = -t241 + t242;
t244 = -t104 - t72 - t85;
t243 = -t105 - t278;
t240 = rSges(4,1) * t266 + rSges(4,2) * t258 + rSges(4,3) * t214;
t238 = (-rSges(5,3) - pkin(8)) * t213;
t237 = -t258 / 0.2e1;
t236 = t256 / 0.2e1;
t25 = t135 * t81 + t136 * t83 - t258 * t79;
t26 = t135 * t82 + t136 * t84 - t258 * t80;
t11 = t211 * t26 + t214 * t25;
t12 = t211 * t28 + t214 * t27;
t50 = -t122 * t258 + t123 * t135 + t124 * t136;
t4 = t50 * t210 + (-t211 * t25 + t214 * t26) * t213;
t235 = t11 * t237 + t12 * t236 + t4 * t282 + t5 * t284 + (t37 * t211 + t36 * t214) * t285;
t234 = t60 + (t36 + t50) * t237 + (t37 + t51) * t236;
t232 = -t258 * t4 + t287;
t230 = -t165 * rSges(5,1) - t164 * rSges(5,2);
t22 = t108 * t258 + t210 * t72 + t275 + t65;
t116 = t125 * t256;
t118 = t127 * t256;
t23 = t108 * t256 + t210 * t243 + t116 + t118;
t228 = t211 * t23 - t214 * t22;
t133 = t210 * rSges(6,3) + (rSges(6,1) * t195 - rSges(6,2) * t194) * t213;
t44 = t133 * t258 + t210 * t93 + t275;
t45 = t133 * t256 + t210 * t276 + t118;
t227 = t211 * t45 - t214 * t44;
t56 = t211 * t255 + t252;
t57 = (t250 - t255) * t214;
t226 = t211 * t56 - t214 * t57;
t66 = -t210 * t86 + t116;
t224 = t211 * t66 - t214 * t65;
t76 = t133 * t211 + t252;
t77 = (-t133 + t250) * t214;
t222 = t211 * t76 - t214 * t77;
t219 = Icges(4,5) * t210 + Icges(4,6) * t213;
t216 = -t38 / 0.2e1 - t46 / 0.2e1 - t61 / 0.2e1 - t54 / 0.2e1;
t215 = t39 / 0.2e1 + t47 / 0.2e1 + t62 / 0.2e1 + t55 / 0.2e1;
t181 = rSges(2,1) * t214 - rSges(2,2) * t211;
t179 = -rSges(2,1) * t211 - rSges(2,2) * t214;
t175 = -t299 + t300;
t156 = -rSges(3,2) * t214 + rSges(3,3) * t211 + t247;
t155 = t214 * rSges(3,3) + t198 + (rSges(3,2) - pkin(1)) * t211;
t154 = t210 * rSges(5,3) + (rSges(5,1) * t212 - rSges(5,2) * t209) * t213;
t147 = Icges(4,3) * t211 - t214 * t219;
t146 = Icges(4,3) * t214 + t211 * t219;
t114 = t239 + t240;
t113 = t198 + t289 + (-rSges(4,3) + t288) * t211;
t111 = (-t154 - t182) * t214;
t110 = t154 * t211 + t169;
t107 = rSges(5,3) * t256 - t230;
t106 = -rSges(5,3) * t258 + t249;
t95 = -t211 * t240 + (t211 * rSges(4,3) - t289) * t214;
t75 = t211 * t238 + t188 + t239 + t249;
t74 = t211 * t288 + t214 * t238 + t190 + t198 + t230;
t71 = -t107 * t210 + t154 * t256;
t70 = t106 * t210 + t154 * t258;
t64 = (-t106 * t214 - t107 * t211) * t213;
t53 = (-t211 * t86 - t214 * t85) * t213;
t52 = t214 * t107 + t158 + (-t106 - t167) * t211;
t29 = (t211 * t276 + t214 * t277) * t213;
t24 = t214 * t94 + (-t167 + t277) * t211 + t274;
t21 = t211 * t43 + t214 * t42;
t20 = t211 * t41 + t214 * t40;
t19 = (t211 * t243 + t214 * t244) * t213;
t18 = t278 * t214 + (-t167 + t244) * t211 + t274;
t16 = t211 * t33 + t214 * t32;
t15 = t211 * t31 + t214 * t30;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t213 - t271 + t296) * t213 + m(7) * (t58 ^ 2 + t59 ^ 2) + m(6) * (t67 ^ 2 + t68 ^ 2) + m(5) * (t74 ^ 2 + t75 ^ 2) + m(4) * (t113 ^ 2 + t114 ^ 2) + m(3) * (t155 ^ 2 + t156 ^ 2) + m(2) * (t179 ^ 2 + t181 ^ 2) + t254 + t297 + (Icges(4,4) * t301 + Icges(4,2) * t210) * t210; t292 + t293 + m(5) * (t211 * t74 - t214 * t75) + m(4) * (t113 * t211 - t114 * t214) + m(3) * (t155 * t211 - t156 * t214); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t245) * t246; m(7) * (t56 * t58 + t57 * t59) + m(6) * (t67 * t76 + t68 * t77) + m(5) * (t110 * t74 + t111 * t75) + (t50 / 0.2e1 + t36 / 0.2e1 - t114 * t281 + t175 * t282 - t216 + t298 * t214) * t214 + (t51 / 0.2e1 + t37 / 0.2e1 + t113 * t281 + t175 * t284 + t215 + t298 * t211) * t211; m(5) * (t110 * t211 - t111 * t214) + m(6) * t222 + m(7) * t226 + t246 * t281; m(7) * (t18 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t24 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t110 ^ 2 + t111 ^ 2 + t52 ^ 2) + m(4) * (t180 ^ 2 * t246 + t95 ^ 2) + (t207 * t146 + t11 + t15 + t20) * t214 + (t205 * t147 + t12 + t16 + t21 + (t211 * t146 + t214 * t147) * t214) * t211; m(7) * (t22 * t59 + t23 * t58) + m(6) * (t44 * t68 + t45 * t67) + m(5) * (t70 * t75 + t71 * t74) + (t211 * t216 + t214 * t215) * t213 + t234 + t279; m(5) * (t211 * t71 - t214 * t70) + m(6) * t227 + m(7) * t228; m(7) * (t18 * t19 + t22 * t57 + t23 * t56) + m(6) * (t24 * t29 + t44 * t77 + t45 * t76) + m(5) * (t110 * t71 + t111 * t70 + t52 * t64) + ((t21 / 0.2e1 + t16 / 0.2e1) * t214 + (-t15 / 0.2e1 - t20 / 0.2e1) * t211) * t213 + t235 + (t211 * t290 + t214 * t291) * t285 + t294 * t284 + t295 * t282; t279 * t210 + m(7) * (t19 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(6) * (t29 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(5) * (t64 ^ 2 + t70 ^ 2 + t71 ^ 2) + ((t210 * t290 + t294) * t214 + (-t210 * t291 - t295 - t4) * t211) * t213 + t287; 0.2e1 * (-t292 / 0.2e1 - t293 / 0.2e1) * t213; t245 * t246 * t301; m(7) * (t210 * t18 - t213 * t226) + m(6) * (t210 * t24 - t213 * t222); m(7) * (t210 * t19 - t213 * t228) + m(6) * (t210 * t29 - t213 * t227); 0.2e1 * t245 * (t213 ^ 2 * t246 + t210 ^ 2); m(7) * (t58 * t66 + t59 * t65) + t234; m(7) * t224; m(7) * (t18 * t53 + t56 * t66 + t57 * t65) + t235; m(7) * (t19 * t53 + t22 * t65 + t23 * t66) + t232; m(7) * (t53 * t210 - t213 * t224); m(7) * (t53 ^ 2 + t65 ^ 2 + t66 ^ 2) + t232;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
