% Calculate joint inertia matrix for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:24
% EndTime: 2019-03-09 22:06:37
% DurationCPUTime: 6.51s
% Computational Cost: add. (17396->538), mult. (15104->753), div. (0->0), fcn. (16265->12), ass. (0->270)
t270 = qJ(4) + pkin(11);
t260 = sin(t270);
t261 = cos(t270);
t280 = cos(qJ(1));
t273 = qJ(2) + qJ(3);
t264 = cos(t273);
t277 = sin(qJ(1));
t341 = t264 * t277;
t208 = -t260 * t341 - t261 * t280;
t209 = -t260 * t280 + t261 * t341;
t263 = sin(t273);
t343 = t263 * t277;
t135 = Icges(6,5) * t209 + Icges(6,6) * t208 + Icges(6,3) * t343;
t278 = cos(qJ(4));
t336 = t278 * t280;
t275 = sin(qJ(4));
t339 = t275 * t277;
t223 = -t264 * t339 - t336;
t337 = t277 * t278;
t338 = t275 * t280;
t224 = t264 * t337 - t338;
t153 = Icges(5,5) * t224 + Icges(5,6) * t223 + Icges(5,3) * t343;
t390 = t135 + t153;
t340 = t264 * t280;
t210 = -t260 * t340 + t261 * t277;
t211 = t260 * t277 + t261 * t340;
t342 = t263 * t280;
t136 = Icges(6,5) * t211 + Icges(6,6) * t210 + Icges(6,3) * t342;
t225 = -t264 * t338 + t337;
t226 = t264 * t336 + t339;
t154 = Icges(5,5) * t226 + Icges(5,6) * t225 + Icges(5,3) * t342;
t389 = t136 + t154;
t137 = Icges(6,4) * t209 + Icges(6,2) * t208 + Icges(6,6) * t343;
t139 = Icges(6,1) * t209 + Icges(6,4) * t208 + Icges(6,5) * t343;
t155 = Icges(5,4) * t224 + Icges(5,2) * t223 + Icges(5,6) * t343;
t157 = Icges(5,1) * t224 + Icges(5,4) * t223 + Icges(5,5) * t343;
t388 = t137 * t208 + t139 * t209 + t155 * t223 + t157 * t224 + t390 * t343;
t138 = Icges(6,4) * t211 + Icges(6,2) * t210 + Icges(6,6) * t342;
t140 = Icges(6,1) * t211 + Icges(6,4) * t210 + Icges(6,5) * t342;
t156 = Icges(5,4) * t226 + Icges(5,2) * t225 + Icges(5,6) * t342;
t158 = Icges(5,1) * t226 + Icges(5,4) * t225 + Icges(5,5) * t342;
t387 = t138 * t208 + t140 * t209 + t156 * t223 + t158 * t224 + t389 * t343;
t386 = t137 * t210 + t139 * t211 + t155 * t225 + t157 * t226 + t390 * t342;
t385 = t138 * t210 + t140 * t211 + t156 * t225 + t158 * t226 + t389 * t342;
t177 = -Icges(6,3) * t264 + (Icges(6,5) * t261 - Icges(6,6) * t260) * t263;
t178 = -Icges(6,6) * t264 + (Icges(6,4) * t261 - Icges(6,2) * t260) * t263;
t179 = -Icges(6,5) * t264 + (Icges(6,1) * t261 - Icges(6,4) * t260) * t263;
t85 = t177 * t342 + t178 * t210 + t179 * t211;
t189 = -Icges(5,3) * t264 + (Icges(5,5) * t278 - Icges(5,6) * t275) * t263;
t190 = -Icges(5,6) * t264 + (Icges(5,4) * t278 - Icges(5,2) * t275) * t263;
t191 = -Icges(5,5) * t264 + (Icges(5,1) * t278 - Icges(5,4) * t275) * t263;
t95 = t189 * t342 + t190 * t225 + t191 * t226;
t384 = -t85 - t95;
t84 = t177 * t343 + t178 * t208 + t179 * t209;
t94 = t189 * t343 + t190 * t223 + t191 * t224;
t383 = -t94 - t84;
t262 = qJ(6) + t270;
t255 = sin(t262);
t256 = cos(t262);
t197 = -t255 * t340 + t256 * t277;
t198 = t255 * t277 + t256 * t340;
t131 = t198 * rSges(7,1) + t197 * rSges(7,2) + rSges(7,3) * t342;
t257 = t278 * pkin(4) + pkin(3);
t235 = pkin(5) * t261 + t257;
t236 = pkin(4) * t275 + pkin(5) * t260;
t382 = t235 * t340 + t277 * t236 + t131;
t381 = t383 * t264 + (t277 * t388 + t387 * t280) * t263;
t380 = t384 * t264 + (t277 * t386 + t280 * t385) * t263;
t379 = t387 * t277 - t280 * t388;
t378 = t277 * t385 - t280 * t386;
t274 = -qJ(5) - pkin(9);
t269 = -pkin(10) + t274;
t322 = t269 - t274;
t327 = -pkin(4) * t339 - t257 * t340;
t377 = -t322 * t342 + t327 + t382;
t376 = -t177 - t189;
t375 = -t178 * t260 - t190 * t275;
t326 = t235 - t257;
t161 = t263 * t326 + t264 * t322;
t176 = -rSges(7,3) * t264 + (rSges(7,1) * t256 - rSges(7,2) * t255) * t263;
t374 = -t161 - t176;
t373 = (t179 * t261 + t191 * t278) * t263;
t294 = Icges(4,5) * t264 - Icges(4,6) * t263;
t199 = -Icges(4,3) * t280 + t277 * t294;
t194 = -t255 * t341 - t256 * t280;
t195 = -t255 * t280 + t256 * t341;
t124 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t343;
t126 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t343;
t128 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t343;
t41 = t124 * t343 + t126 * t194 + t128 * t195;
t125 = Icges(7,5) * t198 + Icges(7,6) * t197 + Icges(7,3) * t342;
t127 = Icges(7,4) * t198 + Icges(7,2) * t197 + Icges(7,6) * t342;
t129 = Icges(7,1) * t198 + Icges(7,4) * t197 + Icges(7,5) * t342;
t42 = t125 * t343 + t127 * t194 + t129 * t195;
t20 = t277 * t42 - t280 * t41;
t200 = Icges(4,3) * t277 + t280 * t294;
t272 = t280 ^ 2;
t347 = Icges(4,4) * t264;
t296 = -Icges(4,2) * t263 + t347;
t202 = Icges(4,6) * t277 + t280 * t296;
t348 = Icges(4,4) * t263;
t298 = Icges(4,1) * t264 - t348;
t204 = Icges(4,5) * t277 + t280 * t298;
t292 = -t202 * t263 + t204 * t264;
t201 = -Icges(4,6) * t280 + t277 * t296;
t203 = -Icges(4,5) * t280 + t277 * t298;
t293 = t201 * t263 - t203 * t264;
t372 = -t20 - t272 * t199 - (t292 * t277 + (-t200 + t293) * t280) * t277 - t379;
t271 = t277 ^ 2;
t371 = m(6) / 0.2e1;
t370 = m(7) / 0.2e1;
t172 = -Icges(7,3) * t264 + (Icges(7,5) * t256 - Icges(7,6) * t255) * t263;
t173 = -Icges(7,6) * t264 + (Icges(7,4) * t256 - Icges(7,2) * t255) * t263;
t174 = -Icges(7,5) * t264 + (Icges(7,1) * t256 - Icges(7,4) * t255) * t263;
t77 = t172 * t343 + t173 * t194 + t174 * t195;
t5 = -t264 * t77 + (t277 * t41 + t280 * t42) * t263;
t43 = t124 * t342 + t126 * t197 + t128 * t198;
t44 = t125 * t342 + t127 * t197 + t129 * t198;
t78 = t172 * t342 + t173 * t197 + t174 * t198;
t6 = -t264 * t78 + (t277 * t43 + t280 * t44) * t263;
t369 = t6 * t342 + t5 * t343;
t368 = -t264 / 0.2e1;
t367 = t277 / 0.2e1;
t366 = -t280 / 0.2e1;
t276 = sin(qJ(2));
t365 = pkin(2) * t276;
t364 = pkin(3) * t264;
t363 = pkin(9) * t263;
t362 = -pkin(3) + t257;
t279 = cos(qJ(2));
t361 = rSges(3,1) * t279;
t360 = rSges(3,2) * t276;
t359 = t280 * rSges(3,3);
t54 = -t124 * t264 + (-t126 * t255 + t128 * t256) * t263;
t358 = t54 * t280;
t55 = -t125 * t264 + (-t127 * t255 + t129 * t256) * t263;
t357 = t55 * t277;
t59 = -t135 * t264 + (-t137 * t260 + t139 * t261) * t263;
t356 = t59 * t280;
t60 = -t136 * t264 + (-t138 * t260 + t140 * t261) * t263;
t355 = t60 * t277;
t71 = -t153 * t264 + (-t155 * t275 + t157 * t278) * t263;
t354 = t71 * t280;
t72 = -t154 * t264 + (-t156 * t275 + t158 * t278) * t263;
t353 = t72 * t277;
t165 = t263 * t256 * t174;
t346 = t173 * t255;
t89 = -t264 * t172 - t263 * t346 + t165;
t352 = t89 * t264;
t351 = t263 * t375 + t264 * t376 + t373;
t350 = Icges(3,4) * t276;
t349 = Icges(3,4) * t279;
t281 = -pkin(8) - pkin(7);
t335 = t280 * t281;
t142 = t211 * rSges(6,1) + t210 * rSges(6,2) + rSges(6,3) * t342;
t286 = -t274 * t342 - t327;
t324 = pkin(3) * t340 + pkin(9) * t342;
t152 = t286 - t324;
t334 = -t142 - t152;
t325 = pkin(4) * t338 + t274 * t343;
t151 = (t362 * t264 - t363) * t277 - t325;
t175 = (pkin(9) + t274) * t264 + t362 * t263;
t333 = t264 * t151 + t175 * t343;
t300 = -t195 * rSges(7,1) - t194 * rSges(7,2);
t130 = rSges(7,3) * t343 - t300;
t99 = t264 * t130 + t176 * t343;
t180 = -rSges(6,3) * t264 + (rSges(6,1) * t261 - rSges(6,2) * t260) * t263;
t332 = -t175 - t180;
t258 = pkin(2) * t279 + pkin(1);
t245 = t280 * t258;
t268 = t280 * pkin(7);
t331 = t277 * (t335 + t268 + (-pkin(1) + t258) * t277) + t280 * (-t280 * pkin(1) + t245 + (-pkin(7) - t281) * t277);
t287 = rSges(4,1) * t340 - rSges(4,2) * t342 + t277 * rSges(4,3);
t303 = rSges(4,1) * t264 - rSges(4,2) * t263;
t144 = t277 * (-t280 * rSges(4,3) + t277 * t303) + t280 * t287;
t192 = -rSges(5,3) * t264 + (rSges(5,1) * t278 - rSges(5,2) * t275) * t263;
t234 = pkin(3) * t263 - pkin(9) * t264;
t330 = -t192 - t234;
t329 = t271 * (t363 + t364) + t280 * t324;
t323 = t277 * rSges(3,3) + t280 * t361;
t321 = t271 + t272;
t320 = -t152 - t377;
t319 = -t175 + t374;
t318 = -t234 + t332;
t160 = t226 * rSges(5,1) + t225 * rSges(5,2) + rSges(5,3) * t342;
t317 = t343 / 0.2e1;
t316 = t342 / 0.2e1;
t233 = rSges(4,1) * t263 + rSges(4,2) * t264;
t315 = -t233 - t365;
t314 = -t234 - t365;
t21 = t277 * t44 - t280 * t43;
t313 = (t271 * t200 + t21 + (t293 * t280 + (-t199 + t292) * t277) * t280 + t378) * t277;
t312 = (t54 + t77) * t317 + (t55 + t78) * t316;
t311 = -t277 * t281 + t245;
t11 = -t352 + (t277 * t54 + t280 * t55) * t263;
t310 = -t264 * t11 + t369;
t309 = t277 * t151 + t280 * t152 + t329;
t308 = -t234 + t319;
t302 = -t224 * rSges(5,1) - t223 * rSges(5,2);
t159 = rSges(5,3) * t343 - t302;
t81 = t277 * t159 + t280 * t160 + t329;
t307 = t20 * t317 + t21 * t316 + t5 * t366 + t6 * t367 + (t357 - t358) * t368;
t306 = -t175 + t314;
t305 = -t192 + t314;
t304 = -t360 + t361;
t301 = -t209 * rSges(6,1) - t208 * rSges(6,2);
t299 = Icges(3,1) * t279 - t350;
t297 = -Icges(3,2) * t276 + t349;
t295 = Icges(3,5) * t279 - Icges(3,6) * t276;
t230 = Icges(4,2) * t264 + t348;
t231 = Icges(4,1) * t263 + t347;
t289 = -t230 * t263 + t231 * t264;
t288 = -t180 + t306;
t141 = rSges(6,3) * t343 - t301;
t40 = t277 * t141 + t280 * t142 + t309;
t285 = t306 + t374;
t114 = -t280 * t236 + (-t263 * t269 + t264 * t326) * t277 + t325;
t28 = t309 + t377 * t280 + (t114 + t130) * t277;
t284 = t280 * t372 + t313;
t283 = t307 + (t355 - t356 + t353 - t354) * t368 + t380 * t367 + t381 * t366 + t379 * t317 + t378 * t316;
t229 = Icges(4,5) * t263 + Icges(4,6) * t264;
t282 = -t358 / 0.2e1 + t357 / 0.2e1 - t356 / 0.2e1 + t355 / 0.2e1 - t354 / 0.2e1 + t353 / 0.2e1 + (t202 * t264 + t204 * t263 + t229 * t277 + t280 * t289 - t384 + t78) * t367 + (t201 * t264 + t203 * t263 - t229 * t280 + t277 * t289 - t383 + t77) * t366;
t244 = rSges(2,1) * t280 - rSges(2,2) * t277;
t243 = -rSges(2,1) * t277 - rSges(2,2) * t280;
t242 = rSges(3,1) * t276 + rSges(3,2) * t279;
t216 = Icges(3,3) * t277 + t280 * t295;
t215 = -Icges(3,3) * t280 + t277 * t295;
t196 = t315 * t280;
t193 = t315 * t277;
t182 = t277 * pkin(7) + (pkin(1) - t360) * t280 + t323;
t181 = t359 + t268 + (-pkin(1) - t304) * t277;
t170 = t287 + t311;
t169 = (rSges(4,3) - t281) * t280 + (-t258 - t303) * t277;
t164 = t330 * t280;
t163 = t330 * t277;
t162 = t280 * (-t280 * t360 + t323) + (t277 * t304 - t359) * t277;
t148 = t305 * t280;
t147 = t305 * t277;
t134 = t151 * t342;
t116 = t130 * t342;
t113 = t311 + t160 + t324;
t112 = -t335 + (-t364 - t258 + (-rSges(5,3) - pkin(9)) * t263) * t277 + t302;
t109 = t318 * t280;
t108 = t318 * t277;
t107 = -t160 * t264 - t192 * t342;
t106 = t159 * t264 + t192 * t343;
t105 = t288 * t280;
t104 = t288 * t277;
t102 = t286 + t311 + t142;
t101 = -t335 + (-rSges(6,3) * t263 - t257 * t264 - t258) * t277 + t301 + t325;
t100 = -t131 * t264 - t176 * t342;
t98 = t144 + t331;
t97 = (t159 * t280 - t160 * t277) * t263;
t93 = -t269 * t342 + t311 + t382;
t92 = (t236 - t281) * t280 + (-t235 * t264 - t258 + (-rSges(7,3) + t269) * t263) * t277 + t300;
t88 = t308 * t280;
t87 = t308 * t277;
t86 = -t131 * t343 + t116;
t83 = t285 * t280;
t82 = t285 * t277;
t62 = t264 * t334 + t332 * t342;
t61 = t141 * t264 + t180 * t343 + t333;
t58 = t81 + t331;
t45 = t134 + (t141 * t280 + t277 * t334) * t263;
t39 = t40 + t331;
t38 = t264 * t320 + t319 * t342;
t37 = t114 * t264 + t161 * t343 + t333 + t99;
t29 = t116 + t134 + (t114 * t280 + t277 * t320) * t263;
t19 = t28 + t331;
t1 = [t279 * (Icges(3,2) * t279 + t350) + t276 * (Icges(3,1) * t276 + t349) + Icges(2,3) + t165 + (-t172 + t230 + t376) * t264 + (t231 - t346 + t375) * t263 + m(7) * (t92 ^ 2 + t93 ^ 2) + m(6) * (t101 ^ 2 + t102 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2) + m(4) * (t169 ^ 2 + t170 ^ 2) + m(3) * (t181 ^ 2 + t182 ^ 2) + m(2) * (t243 ^ 2 + t244 ^ 2) + t373; m(7) * (t82 * t93 + t83 * t92) + m(6) * (t101 * t105 + t102 * t104) + m(5) * (t112 * t148 + t113 * t147) + m(4) * (t169 * t196 + t170 * t193) + t282 + m(3) * (-t181 * t280 - t182 * t277) * t242 + (t271 / 0.2e1 + t272 / 0.2e1) * (Icges(3,5) * t276 + Icges(3,6) * t279) + (t279 * (Icges(3,6) * t277 + t280 * t297) + t276 * (Icges(3,5) * t277 + t280 * t299)) * t367 + (t279 * (-Icges(3,6) * t280 + t277 * t297) + t276 * (-Icges(3,5) * t280 + t277 * t299)) * t366; m(7) * (t19 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(6) * (t104 ^ 2 + t105 ^ 2 + t39 ^ 2) + m(5) * (t147 ^ 2 + t148 ^ 2 + t58 ^ 2) + m(4) * (t193 ^ 2 + t196 ^ 2 + t98 ^ 2) + m(3) * (t242 ^ 2 * t321 + t162 ^ 2) + t277 * t271 * t216 + t313 + (-t272 * t215 + (-t277 * t215 + t280 * t216) * t277 + t372) * t280; t282 + m(5) * (t112 * t164 + t113 * t163) + m(7) * (t87 * t93 + t88 * t92) + m(6) * (t101 * t109 + t102 * t108) + m(4) * (-t169 * t280 - t170 * t277) * t233; m(7) * (t19 * t28 + t82 * t87 + t83 * t88) + m(6) * (t104 * t108 + t105 * t109 + t39 * t40) + m(5) * (t147 * t163 + t148 * t164 + t58 * t81) + m(4) * (t144 * t98 + (-t193 * t277 - t196 * t280) * t233) + t284; m(7) * (t28 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(6) * (t108 ^ 2 + t109 ^ 2 + t40 ^ 2) + m(5) * (t163 ^ 2 + t164 ^ 2 + t81 ^ 2) + m(4) * (t233 ^ 2 * t321 + t144 ^ 2) + t284; (-t89 - t351) * t264 + m(7) * (t37 * t92 + t38 * t93) + m(6) * (t101 * t61 + t102 * t62) + m(5) * (t106 * t112 + t107 * t113) + ((t95 / 0.2e1 + t85 / 0.2e1 + t72 / 0.2e1 + t60 / 0.2e1) * t280 + (t94 / 0.2e1 + t84 / 0.2e1 + t71 / 0.2e1 + t59 / 0.2e1) * t277) * t263 + t312; t283 + m(7) * (t19 * t29 + t37 * t83 + t38 * t82) + m(6) * (t104 * t62 + t105 * t61 + t39 * t45) + m(5) * (t106 * t148 + t107 * t147 + t58 * t97); t283 + m(7) * (t28 * t29 + t37 * t88 + t38 * t87) + m(6) * (t108 * t62 + t109 * t61 + t40 * t45) + m(5) * (t106 * t164 + t107 * t163 + t81 * t97); m(7) * (t29 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(6) * (t45 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) + (t264 * t351 - t11) * t264 + (t380 * t280 + t381 * t277 + ((-t60 - t72) * t280 + (-t59 - t71) * t277) * t264) * t263 + t369; 0.2e1 * ((t277 * t93 + t280 * t92) * t370 + (t101 * t280 + t102 * t277) * t371) * t263; m(7) * (-t19 * t264 + (t277 * t82 + t280 * t83) * t263) + m(6) * (-t264 * t39 + (t104 * t277 + t105 * t280) * t263); m(7) * (-t264 * t28 + (t277 * t87 + t280 * t88) * t263) + m(6) * (-t264 * t40 + (t108 * t277 + t109 * t280) * t263); m(7) * (-t264 * t29 + (t277 * t38 + t280 * t37) * t263) + m(6) * (-t264 * t45 + (t277 * t62 + t280 * t61) * t263); 0.2e1 * (t371 + t370) * (t263 ^ 2 * t321 + t264 ^ 2); m(7) * (t100 * t93 + t92 * t99) - t352 + t312; m(7) * (t100 * t82 + t19 * t86 + t83 * t99) + t307; m(7) * (t100 * t87 + t28 * t86 + t88 * t99) + t307; m(7) * (t100 * t38 + t29 * t86 + t37 * t99) + t310; m(7) * (-t264 * t86 + (t100 * t277 + t280 * t99) * t263); m(7) * (t100 ^ 2 + t86 ^ 2 + t99 ^ 2) + t310;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
