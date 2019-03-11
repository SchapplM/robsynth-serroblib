% Calculate joint inertia matrix for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:45
% EndTime: 2019-03-09 16:30:51
% DurationCPUTime: 4.57s
% Computational Cost: add. (10214->430), mult. (9124->611), div. (0->0), fcn. (9596->10), ass. (0->214)
t215 = qJ(2) + qJ(3);
t204 = pkin(10) + t215;
t201 = cos(t204);
t220 = cos(qJ(5));
t222 = cos(qJ(1));
t274 = t220 * t222;
t217 = sin(qJ(5));
t219 = sin(qJ(1));
t277 = t217 * t219;
t165 = -t201 * t277 - t274;
t275 = t219 * t220;
t276 = t217 * t222;
t166 = t201 * t275 - t276;
t200 = sin(t204);
t281 = t200 * t219;
t93 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t281;
t95 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t281;
t341 = t93 + t95;
t167 = -t201 * t276 + t275;
t168 = t201 * t274 + t277;
t279 = t200 * t222;
t94 = Icges(7,5) * t168 + Icges(7,6) * t167 + Icges(7,3) * t279;
t96 = Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t279;
t340 = t94 + t96;
t97 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t281;
t99 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t281;
t339 = t97 + t99;
t100 = Icges(6,4) * t168 + Icges(6,2) * t167 + Icges(6,6) * t279;
t98 = Icges(7,4) * t168 + Icges(7,2) * t167 + Icges(7,6) * t279;
t338 = t100 + t98;
t101 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t281;
t103 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t281;
t337 = t101 + t103;
t102 = Icges(7,1) * t168 + Icges(7,4) * t167 + Icges(7,5) * t279;
t104 = Icges(6,1) * t168 + Icges(6,4) * t167 + Icges(6,5) * t279;
t336 = t102 + t104;
t216 = -qJ(6) - pkin(9);
t335 = rSges(7,3) - t216;
t334 = Icges(4,3) + Icges(5,3);
t205 = sin(t215);
t206 = cos(t215);
t333 = Icges(4,5) * t206 + Icges(5,5) * t201 - Icges(4,6) * t205 - Icges(5,6) * t200;
t332 = t339 * t165 + t337 * t166 + t281 * t341;
t331 = t165 * t338 + t166 * t336 + t281 * t340;
t330 = t339 * t167 + t337 * t168 + t279 * t341;
t329 = t167 * t338 + t168 * t336 + t279 * t340;
t121 = -Icges(7,3) * t201 + (Icges(7,5) * t220 - Icges(7,6) * t217) * t200;
t123 = -Icges(7,6) * t201 + (Icges(7,4) * t220 - Icges(7,2) * t217) * t200;
t125 = -Icges(7,5) * t201 + (Icges(7,1) * t220 - Icges(7,4) * t217) * t200;
t54 = t121 * t281 + t123 * t165 + t125 * t166;
t122 = -Icges(6,3) * t201 + (Icges(6,5) * t220 - Icges(6,6) * t217) * t200;
t124 = -Icges(6,6) * t201 + (Icges(6,4) * t220 - Icges(6,2) * t217) * t200;
t126 = -Icges(6,5) * t201 + (Icges(6,1) * t220 - Icges(6,4) * t217) * t200;
t55 = t122 * t281 + t124 * t165 + t126 * t166;
t328 = -t54 - t55;
t56 = t121 * t279 + t123 * t167 + t125 * t168;
t57 = t122 * t279 + t124 * t167 + t126 * t168;
t327 = -t56 - t57;
t202 = pkin(5) * t220 + pkin(4);
t278 = t201 * t222;
t326 = t168 * rSges(7,1) + t167 * rSges(7,2) + pkin(5) * t277 + t202 * t278 + t335 * t279;
t325 = t334 * t219 + t333 * t222;
t324 = -t333 * t219 + t334 * t222;
t283 = Icges(5,4) * t201;
t246 = -Icges(5,2) * t200 + t283;
t140 = Icges(5,6) * t219 + t222 * t246;
t284 = Icges(5,4) * t200;
t249 = Icges(5,1) * t201 - t284;
t142 = Icges(5,5) * t219 + t222 * t249;
t285 = Icges(4,4) * t206;
t247 = -Icges(4,2) * t205 + t285;
t152 = Icges(4,6) * t219 + t222 * t247;
t286 = Icges(4,4) * t205;
t250 = Icges(4,1) * t206 - t286;
t154 = Icges(4,5) * t219 + t222 * t250;
t323 = t140 * t200 - t142 * t201 + t152 * t205 - t154 * t206;
t139 = -Icges(5,6) * t222 + t219 * t246;
t141 = -Icges(5,5) * t222 + t219 * t249;
t151 = -Icges(4,6) * t222 + t219 * t247;
t153 = -Icges(4,5) * t222 + t219 * t250;
t322 = t139 * t200 - t141 * t201 + t151 * t205 - t153 * t206;
t213 = t219 ^ 2;
t321 = t328 * t201 + (t219 * t332 + t331 * t222) * t200;
t320 = t327 * t201 + (t219 * t330 + t222 * t329) * t200;
t319 = t219 * pkin(7);
t318 = t331 * t219 - t222 * t332;
t317 = t219 * t329 - t222 * t330;
t252 = -rSges(7,1) * t166 - rSges(7,2) * t165;
t301 = pkin(9) + t216;
t302 = -pkin(4) + t202;
t290 = rSges(7,3) * t281 - t252 - pkin(5) * t276 + (-t200 * t301 + t201 * t302) * t219;
t270 = pkin(4) * t278 + pkin(9) * t279;
t316 = -t270 + t326;
t315 = -t121 - t122;
t314 = Icges(4,5) * t205 + Icges(5,5) * t200 + Icges(4,6) * t206 + Icges(5,6) * t201;
t170 = Icges(5,2) * t201 + t284;
t171 = Icges(5,1) * t200 + t283;
t178 = Icges(4,2) * t206 + t286;
t179 = Icges(4,1) * t205 + t285;
t313 = -t170 * t200 + t171 * t201 - t178 * t205 + t179 * t206;
t255 = rSges(4,1) * t206 - rSges(4,2) * t205;
t312 = (-t123 - t124) * t217;
t311 = (t125 + t126) * t200 * t220;
t214 = t222 ^ 2;
t310 = -t318 + t324 * t214 + (t323 * t219 + (-t322 + t325) * t222) * t219;
t309 = t201 ^ 2;
t223 = -pkin(8) - pkin(7);
t307 = t219 / 0.2e1;
t306 = -t222 / 0.2e1;
t218 = sin(qJ(2));
t305 = pkin(2) * t218;
t304 = pkin(3) * t205;
t303 = pkin(4) * t201;
t221 = cos(qJ(2));
t203 = t221 * pkin(2) + pkin(1);
t300 = t200 * t312 + t201 * t315 + t311;
t299 = rSges(3,1) * t221;
t297 = rSges(3,2) * t218;
t295 = t222 * rSges(3,3);
t45 = -t201 * t93 + (t101 * t220 - t217 * t97) * t200;
t294 = t45 * t222;
t46 = -t201 * t94 + (t102 * t220 - t217 * t98) * t200;
t293 = t46 * t219;
t47 = -t201 * t95 + (t103 * t220 - t217 * t99) * t200;
t292 = t47 * t222;
t48 = -t201 * t96 + (-t100 * t217 + t104 * t220) * t200;
t291 = t48 * t219;
t288 = Icges(3,4) * t218;
t287 = Icges(3,4) * t221;
t181 = pkin(3) * t206 + t203;
t175 = t222 * t181;
t194 = t222 * t203;
t273 = t222 * (t175 - t194) + (t181 - t203) * t213;
t272 = (t301 - rSges(7,3)) * t201 + (rSges(7,1) * t220 - rSges(7,2) * t217 + t302) * t200;
t211 = t222 * pkin(7);
t271 = t219 * (t211 + (-pkin(1) + t203) * t219) + t222 * (-t222 * pkin(1) + t194 - t319);
t233 = t219 * rSges(4,3) + t222 * t255;
t92 = t219 * (-t222 * rSges(4,3) + t219 * t255) + t222 * t233;
t269 = t219 * rSges(3,3) + t222 * t299;
t267 = t213 + t214;
t108 = t168 * rSges(6,1) + t167 * rSges(6,2) + rSges(6,3) * t279;
t180 = rSges(4,1) * t205 + rSges(4,2) * t206;
t264 = -t180 - t305;
t172 = rSges(5,1) * t200 + rSges(5,2) * t201;
t263 = -t172 - t304;
t173 = pkin(4) * t200 - pkin(9) * t201;
t262 = -t173 - t304;
t261 = (t325 * t213 + ((-t323 + t324) * t219 + t322 * t222) * t222 + t317) * t219;
t212 = -qJ(4) + t223;
t260 = -t219 * t212 + t175;
t232 = rSges(5,1) * t278 - rSges(5,2) * t279 + t219 * rSges(5,3);
t254 = rSges(5,1) * t201 - rSges(5,2) * t200;
t53 = t219 * (-t222 * rSges(5,3) + t219 * t254) + t222 * t232 + t273;
t259 = t213 * (pkin(9) * t200 + t303) + t222 * t270 + t273;
t128 = -t201 * rSges(6,3) + (rSges(6,1) * t220 - rSges(6,2) * t217) * t200;
t258 = -t128 + t262;
t257 = -t304 - t305;
t256 = -t297 + t299;
t253 = -rSges(6,1) * t166 - rSges(6,2) * t165;
t251 = Icges(3,1) * t221 - t288;
t248 = -Icges(3,2) * t218 + t287;
t245 = Icges(3,5) * t221 - Icges(3,6) * t218;
t234 = t262 - t272;
t231 = -t172 + t257;
t230 = -t173 + t257;
t106 = rSges(6,3) * t281 - t253;
t25 = t219 * t106 + t222 * t108 + t259;
t228 = -t128 + t230;
t22 = t219 * t290 + t222 * t316 + t259;
t227 = t230 - t272;
t226 = -(t293 - t294 + t291 - t292) * t201 / 0.2e1 + t320 * t307 + t321 * t306 + t318 * t281 / 0.2e1 + t317 * t279 / 0.2e1;
t225 = t310 * t222 + t261;
t224 = -t294 / 0.2e1 + t293 / 0.2e1 - t292 / 0.2e1 + t291 / 0.2e1 + (t140 * t201 + t142 * t200 + t152 * t206 + t154 * t205 + t219 * t314 + t222 * t313 - t327) * t307 + (t139 * t201 + t141 * t200 + t151 * t206 + t153 * t205 + t219 * t313 - t222 * t314 - t328) * t306;
t193 = rSges(2,1) * t222 - rSges(2,2) * t219;
t192 = -rSges(2,1) * t219 - rSges(2,2) * t222;
t191 = rSges(3,1) * t218 + rSges(3,2) * t221;
t160 = Icges(3,3) * t219 + t222 * t245;
t159 = -Icges(3,3) * t222 + t219 * t245;
t148 = t264 * t222;
t147 = t264 * t219;
t134 = t319 + (pkin(1) - t297) * t222 + t269;
t133 = t295 + t211 + (-pkin(1) - t256) * t219;
t132 = t263 * t222;
t131 = t263 * t219;
t120 = t231 * t222;
t119 = t231 * t219;
t118 = -t219 * t223 + t194 + t233;
t117 = (rSges(4,3) - t223) * t222 + (-t203 - t255) * t219;
t111 = t222 * (-t222 * t297 + t269) + (t219 * t256 - t295) * t219;
t110 = t232 + t260;
t109 = (rSges(5,3) - t212) * t222 + (-t181 - t254) * t219;
t83 = t258 * t222;
t82 = t258 * t219;
t77 = t228 * t222;
t76 = t228 * t219;
t71 = t260 + t108 + t270;
t70 = -t222 * t212 + (-t303 - t181 + (-rSges(6,3) - pkin(9)) * t200) * t219 + t253;
t69 = -t108 * t201 - t128 * t279;
t68 = t106 * t201 + t128 * t281;
t67 = t234 * t222;
t66 = t234 * t219;
t65 = t92 + t271;
t64 = t260 + t326;
t63 = (pkin(5) * t217 - t212) * t222 + (-t335 * t200 - t201 * t202 - t181) * t219 + t252;
t62 = t227 * t222;
t61 = t227 * t219;
t58 = (t106 * t222 - t108 * t219) * t200;
t38 = t53 + t271;
t29 = -t201 * t316 - t272 * t279;
t28 = t201 * t290 + t272 * t281;
t24 = (-t219 * t316 + t222 * t290) * t200;
t23 = t25 + t271;
t19 = t22 + t271;
t1 = [t206 * t178 + t205 * t179 + t221 * (Icges(3,2) * t221 + t288) + t218 * (Icges(3,1) * t218 + t287) + Icges(2,3) + (t170 + t315) * t201 + (t171 + t312) * t200 + m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t70 ^ 2 + t71 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2) + m(4) * (t117 ^ 2 + t118 ^ 2) + m(3) * (t133 ^ 2 + t134 ^ 2) + m(2) * (t192 ^ 2 + t193 ^ 2) + t311; t224 + m(7) * (t61 * t64 + t62 * t63) + m(6) * (t70 * t77 + t71 * t76) + m(5) * (t109 * t120 + t110 * t119) + m(4) * (t117 * t148 + t118 * t147) + m(3) * (-t133 * t222 - t134 * t219) * t191 + (t221 * (-Icges(3,6) * t222 + t219 * t248) + t218 * (-Icges(3,5) * t222 + t219 * t251)) * t306 + (t221 * (Icges(3,6) * t219 + t222 * t248) + t218 * (Icges(3,5) * t219 + t222 * t251)) * t307 + (t213 / 0.2e1 + t214 / 0.2e1) * (Icges(3,5) * t218 + Icges(3,6) * t221); m(7) * (t19 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t23 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t119 ^ 2 + t120 ^ 2 + t38 ^ 2) + m(4) * (t147 ^ 2 + t148 ^ 2 + t65 ^ 2) + m(3) * (t191 ^ 2 * t267 + t111 ^ 2) + t219 * t213 * t160 + t261 + (-t214 * t159 + (-t219 * t159 + t222 * t160) * t219 + t310) * t222; t224 + m(4) * (-t117 * t222 - t118 * t219) * t180 + m(7) * (t63 * t67 + t64 * t66) + m(6) * (t70 * t83 + t71 * t82) + m(5) * (t109 * t132 + t110 * t131); m(7) * (t19 * t22 + t61 * t66 + t62 * t67) + m(6) * (t25 * t23 + t76 * t82 + t77 * t83) + m(5) * (t119 * t131 + t120 * t132 + t53 * t38) + m(4) * (t65 * t92 + (-t147 * t219 - t148 * t222) * t180) + t225; m(7) * (t22 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t25 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(5) * (t131 ^ 2 + t132 ^ 2 + t53 ^ 2) + m(4) * (t180 ^ 2 * t267 + t92 ^ 2) + t225; m(7) * (t219 * t63 - t222 * t64) + m(6) * (t219 * t70 - t222 * t71) + m(5) * (t109 * t219 - t110 * t222); m(7) * (t219 * t62 - t222 * t61) + m(6) * (t219 * t77 - t222 * t76) + m(5) * (-t119 * t222 + t120 * t219); m(7) * (t219 * t67 - t222 * t66) + m(6) * (t219 * t83 - t222 * t82) + m(5) * (-t131 * t222 + t132 * t219); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t267; -t300 * t201 + m(7) * (t28 * t63 + t29 * t64) + m(6) * (t68 * t70 + t69 * t71) + ((t46 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1 + t48 / 0.2e1) * t222 + (t55 / 0.2e1 + t54 / 0.2e1 + t47 / 0.2e1 + t45 / 0.2e1) * t219) * t200; m(7) * (t19 * t24 + t28 * t62 + t29 * t61) + m(6) * (t58 * t23 + t68 * t77 + t69 * t76) + t226; m(7) * (t22 * t24 + t28 * t67 + t29 * t66) + m(6) * (t58 * t25 + t68 * t83 + t69 * t82) + t226; m(6) * (t219 * t68 - t222 * t69) + m(7) * (t219 * t28 - t222 * t29); m(7) * (t24 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t58 ^ 2 + t68 ^ 2 + t69 ^ 2) + t300 * t309 + (t320 * t222 + t321 * t219 + ((-t46 - t48) * t222 + (-t45 - t47) * t219) * t201) * t200; m(7) * (t219 * t64 + t222 * t63) * t200; m(7) * (-t19 * t201 + (t219 * t61 + t222 * t62) * t200); m(7) * (-t201 * t22 + (t219 * t66 + t222 * t67) * t200); 0; m(7) * (-t201 * t24 + (t219 * t29 + t222 * t28) * t200); m(7) * (t200 ^ 2 * t267 + t309);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
