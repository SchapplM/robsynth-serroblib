% Calculate joint inertia matrix for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:36
% EndTime: 2019-03-10 01:26:51
% DurationCPUTime: 7.10s
% Computational Cost: add. (20236->567), mult. (21233->776), div. (0->0), fcn. (23196->10), ass. (0->286)
t290 = qJ(3) + qJ(4);
t283 = qJ(5) + t290;
t278 = sin(t283);
t279 = cos(t283);
t296 = cos(qJ(1));
t293 = sin(qJ(1));
t295 = cos(qJ(2));
t361 = t293 * t295;
t223 = t278 * t361 + t279 * t296;
t224 = -t278 * t296 + t279 * t361;
t292 = sin(qJ(2));
t364 = t292 * t293;
t145 = Icges(7,5) * t224 + Icges(7,6) * t364 + Icges(7,3) * t223;
t151 = Icges(6,4) * t224 - Icges(6,2) * t223 + Icges(6,6) * t364;
t412 = t145 - t151;
t360 = t295 * t296;
t225 = t278 * t360 - t279 * t293;
t226 = t278 * t293 + t279 * t360;
t363 = t292 * t296;
t146 = Icges(7,5) * t226 + Icges(7,6) * t363 + Icges(7,3) * t225;
t152 = Icges(6,4) * t226 - Icges(6,2) * t225 + Icges(6,6) * t363;
t411 = t146 - t152;
t147 = Icges(6,5) * t224 - Icges(6,6) * t223 + Icges(6,3) * t364;
t149 = Icges(7,4) * t224 + Icges(7,2) * t364 + Icges(7,6) * t223;
t410 = t147 + t149;
t148 = Icges(6,5) * t226 - Icges(6,6) * t225 + Icges(6,3) * t363;
t150 = Icges(7,4) * t226 + Icges(7,2) * t363 + Icges(7,6) * t225;
t409 = t148 + t150;
t153 = Icges(7,1) * t224 + Icges(7,4) * t364 + Icges(7,5) * t223;
t155 = Icges(6,1) * t224 - Icges(6,4) * t223 + Icges(6,5) * t364;
t408 = t153 + t155;
t154 = Icges(7,1) * t226 + Icges(7,4) * t363 + Icges(7,5) * t225;
t156 = Icges(6,1) * t226 - Icges(6,4) * t225 + Icges(6,5) * t363;
t407 = t154 + t156;
t387 = rSges(7,3) + qJ(6);
t394 = rSges(7,1) + pkin(5);
t406 = -t223 * t387 - t224 * t394;
t405 = t223 * t412 + t224 * t408 + t364 * t410;
t404 = t223 * t411 + t224 * t407 + t364 * t409;
t403 = t225 * t412 + t226 * t408 + t363 * t410;
t402 = t225 * t411 + t226 * t407 + t363 * t409;
t207 = -Icges(6,6) * t295 + (Icges(6,4) * t279 - Icges(6,2) * t278) * t292;
t368 = t278 * t292;
t204 = -Icges(7,6) * t295 + (Icges(7,5) * t279 + Icges(7,3) * t278) * t292;
t208 = -Icges(7,4) * t295 + (Icges(7,1) * t279 + Icges(7,5) * t278) * t292;
t209 = -Icges(6,5) * t295 + (Icges(6,1) * t279 - Icges(6,4) * t278) * t292;
t396 = t204 * t368 + (t208 + t209) * t279 * t292;
t205 = -Icges(6,3) * t295 + (Icges(6,5) * t279 - Icges(6,6) * t278) * t292;
t206 = -Icges(7,2) * t295 + (Icges(7,4) * t279 + Icges(7,6) * t278) * t292;
t398 = -t205 - t206;
t386 = -t207 * t368 + t295 * t398 + t396;
t401 = t386 * t295;
t100 = t204 * t223 + t206 * t364 + t208 * t224;
t101 = t205 * t364 - t207 * t223 + t209 * t224;
t400 = -t101 - t100;
t102 = t204 * t225 + t206 * t363 + t208 * t226;
t103 = t205 * t363 - t207 * t225 + t209 * t226;
t399 = -t102 - t103;
t397 = Icges(3,5) * t292;
t395 = t397 / 0.2e1;
t393 = t400 * t295 + (t293 * t405 + t296 * t404) * t292;
t392 = t399 * t295 + (t293 * t403 + t296 * t402) * t292;
t391 = t293 * t404 - t296 * t405;
t390 = t293 * t402 - t296 * t403;
t74 = -t295 * t149 + (t145 * t278 + t153 * t279) * t292;
t76 = -t295 * t147 + (-t151 * t278 + t155 * t279) * t292;
t389 = -t74 - t76;
t75 = -t295 * t150 + (t146 * t278 + t154 * t279) * t292;
t77 = -t295 * t148 + (-t152 * t278 + t156 * t279) * t292;
t388 = t75 + t77;
t385 = rSges(7,2) * t364 - t406;
t384 = -t295 * rSges(7,2) + (t278 * t387 + t279 * t394) * t292;
t383 = t293 ^ 2;
t382 = t296 ^ 2;
t297 = -pkin(9) - pkin(8);
t381 = t293 / 0.2e1;
t380 = -t295 / 0.2e1;
t379 = -t296 / 0.2e1;
t378 = t296 / 0.2e1;
t377 = pkin(2) * t295;
t376 = pkin(8) * t292;
t294 = cos(qJ(3));
t280 = pkin(3) * t294 + pkin(2);
t375 = -pkin(2) + t280;
t374 = t401 + (t293 * t389 - t296 * t388) * t292;
t373 = t296 * rSges(3,3);
t371 = Icges(3,4) * t295;
t281 = sin(t290);
t282 = cos(t290);
t215 = -Icges(5,6) * t295 + (Icges(5,4) * t282 - Icges(5,2) * t281) * t292;
t370 = t215 * t281;
t291 = sin(qJ(3));
t233 = -Icges(4,6) * t295 + (Icges(4,4) * t294 - Icges(4,2) * t291) * t292;
t369 = t233 * t291;
t366 = t291 * t293;
t365 = t291 * t296;
t362 = t292 * t297;
t341 = pkin(3) * t365 + t293 * t362;
t255 = pkin(4) * t282 + t280;
t343 = t255 - t280;
t257 = pkin(3) * t291 + pkin(4) * t281;
t289 = -pkin(10) + t297;
t344 = -t257 * t296 - t289 * t364;
t134 = t343 * t361 + t341 + t344;
t127 = t134 * t363;
t314 = -t224 * rSges(6,1) + t223 * rSges(6,2);
t158 = rSges(6,3) * t364 - t314;
t141 = t158 * t363;
t358 = t127 + t141;
t338 = t289 - t297;
t190 = t292 * t343 + t295 * t338;
t357 = t134 * t295 + t190 * t364;
t342 = -pkin(3) * t366 - t280 * t360;
t346 = t255 * t360 + t257 * t293;
t135 = -t338 * t363 + t342 + t346;
t160 = rSges(6,1) * t226 - rSges(6,2) * t225 + rSges(6,3) * t363;
t356 = -t135 - t160;
t355 = t385 * t363;
t354 = rSges(7,2) * t363 + t225 * t387 + t226 * t394;
t242 = -t281 * t360 + t282 * t293;
t243 = t281 * t293 + t282 * t360;
t171 = rSges(5,1) * t243 + rSges(5,2) * t242 + rSges(5,3) * t363;
t306 = -t296 * t362 - t342;
t340 = pkin(2) * t360 + pkin(8) * t363;
t189 = t306 - t340;
t353 = -t171 - t189;
t188 = (t295 * t375 - t376) * t293 - t341;
t213 = (pkin(8) + t297) * t295 + t375 * t292;
t352 = t188 * t295 + t213 * t364;
t212 = -t295 * rSges(6,3) + (rSges(6,1) * t279 - rSges(6,2) * t278) * t292;
t351 = -t190 - t212;
t121 = t158 * t295 + t212 * t364;
t240 = -t281 * t361 - t282 * t296;
t241 = -t281 * t296 + t282 * t361;
t315 = -t241 * rSges(5,1) - t240 * rSges(5,2);
t170 = rSges(5,3) * t364 - t315;
t222 = -t295 * rSges(5,3) + (rSges(5,1) * t282 - rSges(5,2) * t281) * t292;
t125 = t170 * t295 + t222 * t364;
t348 = -t213 - t222;
t239 = -t295 * rSges(4,3) + (rSges(4,1) * t294 - rSges(4,2) * t291) * t292;
t266 = t292 * pkin(2) - t295 * pkin(8);
t347 = -t239 - t266;
t345 = t383 * (t376 + t377) + t296 * t340;
t339 = pkin(1) * t296 + pkin(7) * t293;
t216 = -Icges(5,5) * t295 + (Icges(5,1) * t282 - Icges(5,4) * t281) * t292;
t196 = t292 * t282 * t216;
t214 = -Icges(5,3) * t295 + (Icges(5,5) * t282 - Icges(5,6) * t281) * t292;
t120 = -t295 * t214 - t292 * t370 + t196;
t164 = Icges(5,5) * t241 + Icges(5,6) * t240 + Icges(5,3) * t364;
t166 = Icges(5,4) * t241 + Icges(5,2) * t240 + Icges(5,6) * t364;
t168 = Icges(5,1) * t241 + Icges(5,4) * t240 + Icges(5,5) * t364;
t80 = -t295 * t164 + (-t166 * t281 + t168 * t282) * t292;
t165 = Icges(5,5) * t243 + Icges(5,6) * t242 + Icges(5,3) * t363;
t167 = Icges(5,4) * t243 + Icges(5,2) * t242 + Icges(5,6) * t363;
t169 = Icges(5,1) * t243 + Icges(5,4) * t242 + Icges(5,5) * t363;
t81 = -t295 * t165 + (-t167 * t281 + t169 * t282) * t292;
t337 = t120 * t295 - (t293 * t80 + t296 * t81) * t292 + t374;
t336 = t363 * t392 + t364 * t393;
t230 = -Icges(4,3) * t295 + (Icges(4,5) * t294 - Icges(4,6) * t291) * t292;
t236 = -Icges(4,5) * t295 + (Icges(4,1) * t294 - Icges(4,4) * t291) * t292;
t251 = -t291 * t360 + t293 * t294;
t252 = t294 * t360 + t366;
t118 = t230 * t363 + t233 * t251 + t236 * t252;
t181 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t363;
t183 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t363;
t185 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t363;
t93 = -t295 * t181 + (-t183 * t291 + t185 * t294) * t292;
t335 = t118 / 0.2e1 + t93 / 0.2e1;
t249 = -t291 * t361 - t294 * t296;
t250 = t294 * t361 - t365;
t117 = t230 * t364 + t233 * t249 + t236 * t250;
t180 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t364;
t182 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t364;
t184 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t364;
t92 = -t295 * t180 + (-t182 * t291 + t184 * t294) * t292;
t334 = t92 / 0.2e1 + t117 / 0.2e1;
t333 = -t120 - t386;
t332 = t127 + t355;
t331 = -t135 - t354;
t330 = -t189 + t356;
t329 = -t190 - t384;
t328 = -t213 + t351;
t327 = -t266 + t348;
t187 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t363;
t287 = t296 * pkin(7);
t326 = t287 - t344;
t325 = t364 / 0.2e1;
t324 = t363 / 0.2e1;
t323 = -t255 * t295 - pkin(1);
t322 = -t189 + t331;
t321 = t188 * t293 + t189 * t296 + t345;
t64 = t121 + t357;
t320 = -t213 + t329;
t319 = -t266 + t328;
t82 = t295 * t385 + t364 * t384;
t108 = t214 * t364 + t215 * t240 + t216 * t241;
t66 = t164 * t364 + t166 * t240 + t168 * t241;
t67 = t165 * t364 + t167 * t240 + t169 * t241;
t17 = -t108 * t295 + (t293 * t66 + t296 * t67) * t292;
t109 = t214 * t363 + t215 * t242 + t216 * t243;
t68 = t164 * t363 + t166 * t242 + t168 * t243;
t69 = t165 * t363 + t167 * t242 + t169 * t243;
t18 = -t109 * t295 + (t293 * t68 + t296 * t69) * t292;
t318 = t17 * t364 + t18 * t363 + t336;
t317 = rSges(3,1) * t295 - rSges(3,2) * t292;
t316 = -t250 * rSges(4,1) - t249 * rSges(4,2);
t313 = -t266 + t320;
t311 = -Icges(3,2) * t292 + t371;
t310 = Icges(3,5) * t295 - Icges(3,6) * t292;
t307 = rSges(3,1) * t360 - rSges(3,2) * t363 + rSges(3,3) * t293;
t305 = t134 * t293 + t135 * t296 + t321;
t48 = t82 + t357;
t304 = t295 * t374 + t336;
t303 = (-t389 - t400) * t325 + (t388 - t399) * t324;
t302 = -t289 * t363 + t339 + t346;
t301 = t392 * t381 + (t293 * t388 + t296 * t389) * t380 + t393 * t379 + t391 * t325 + t390 * t324;
t300 = t295 * t337 + t318;
t299 = t303 + (t108 + t80) * t325 + (t109 + t81) * t324;
t36 = t293 * t67 - t296 * t66;
t37 = t293 * t69 - t296 * t68;
t298 = t17 * t379 + t18 * t381 + t37 * t324 + t36 * t325 + t301 + (t81 * t293 - t80 * t296) * t380;
t264 = rSges(2,1) * t296 - rSges(2,2) * t293;
t263 = -rSges(2,1) * t293 - rSges(2,2) * t296;
t262 = rSges(3,1) * t292 + rSges(3,2) * t295;
t259 = Icges(3,6) * t295 + t397;
t232 = Icges(3,3) * t293 + t296 * t310;
t231 = -Icges(3,3) * t296 + t293 * t310;
t210 = t292 * t294 * t236;
t202 = t307 + t339;
t201 = t373 + t287 + (-pkin(1) - t317) * t293;
t192 = t347 * t296;
t191 = t347 * t293;
t186 = rSges(4,3) * t364 - t316;
t175 = t296 * t307 + (t293 * t317 - t373) * t293;
t174 = t188 * t363;
t144 = t170 * t363;
t139 = t187 + t339 + t340;
t138 = t287 + (-t377 - pkin(1) + (-rSges(4,3) - pkin(8)) * t292) * t293 + t316;
t137 = t327 * t296;
t136 = t327 * t293;
t133 = -t187 * t295 - t239 * t363;
t132 = t186 * t295 + t239 * t364;
t128 = -t295 * t230 - t292 * t369 + t210;
t126 = -t295 * t171 - t222 * t363;
t124 = t306 + t171 + t339;
t123 = t287 + (-rSges(5,3) * t292 - t280 * t295 - pkin(1)) * t293 + t315 + t341;
t122 = -t295 * t160 - t212 * t363;
t119 = (t186 * t296 - t187 * t293) * t292;
t114 = t302 + t160;
t113 = (-rSges(6,3) * t292 + t323) * t293 + t314 + t326;
t112 = t319 * t296;
t111 = t319 * t293;
t110 = -t171 * t364 + t144;
t107 = -t160 * t364 + t141;
t104 = t186 * t293 + t187 * t296 + t345;
t99 = t313 * t296;
t98 = t313 * t293;
t91 = t302 + t354;
t90 = (-rSges(7,2) * t292 + t323) * t293 + t326 + t406;
t89 = t295 * t353 + t348 * t363;
t88 = t125 + t352;
t87 = t181 * t363 + t183 * t251 + t185 * t252;
t86 = t180 * t363 + t182 * t251 + t184 * t252;
t85 = t181 * t364 + t183 * t249 + t185 * t250;
t84 = t180 * t364 + t182 * t249 + t184 * t250;
t83 = -t295 * t354 - t363 * t384;
t65 = t295 * t356 + t351 * t363;
t63 = t353 * t364 + t144 + t174;
t54 = -t354 * t364 + t355;
t53 = t170 * t293 + t171 * t296 + t321;
t52 = t356 * t364 + t358;
t51 = t295 * t330 + t328 * t363;
t50 = t64 + t352;
t49 = t295 * t331 + t329 * t363;
t47 = t293 * t87 - t296 * t86;
t46 = t293 * t85 - t296 * t84;
t44 = t295 * t322 + t320 * t363;
t43 = t48 + t352;
t42 = t330 * t364 + t174 + t358;
t41 = t331 * t364 + t332;
t38 = t158 * t293 + t160 * t296 + t305;
t29 = -t118 * t295 + (t293 * t86 + t296 * t87) * t292;
t28 = -t117 * t295 + (t293 * t84 + t296 * t85) * t292;
t22 = t322 * t364 + t174 + t332;
t21 = t293 * t385 + t296 * t354 + t305;
t1 = [Icges(2,3) + t196 + t210 + (Icges(3,1) * t292 - t207 * t278 - t369 - t370 + t371) * t292 + (Icges(3,4) * t292 + Icges(3,2) * t295 - t214 - t230 + t398) * t295 + m(7) * (t90 ^ 2 + t91 ^ 2) + m(6) * (t113 ^ 2 + t114 ^ 2) + m(5) * (t123 ^ 2 + t124 ^ 2) + m(4) * (t138 ^ 2 + t139 ^ 2) + m(3) * (t201 ^ 2 + t202 ^ 2) + m(2) * (t263 ^ 2 + t264 ^ 2) + t396; (-t80 / 0.2e1 - t76 / 0.2e1 - t74 / 0.2e1 - t108 / 0.2e1 - t100 / 0.2e1 - t101 / 0.2e1 + (-Icges(3,6) * t296 + t293 * t311) * t380 + t296 * t395 + t259 * t378 - t334) * t296 + (t81 / 0.2e1 + t77 / 0.2e1 + t75 / 0.2e1 + t109 / 0.2e1 + t102 / 0.2e1 + t103 / 0.2e1 + t295 * (Icges(3,6) * t293 + t296 * t311) / 0.2e1 + t293 * t395 + t259 * t381 + t335) * t293 + m(7) * (t90 * t99 + t91 * t98) + m(6) * (t111 * t114 + t112 * t113) + m(5) * (t123 * t137 + t124 * t136) + m(4) * (t138 * t192 + t139 * t191) + m(3) * (-t201 * t296 - t202 * t293) * t262; m(7) * (t21 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(6) * (t111 ^ 2 + t112 ^ 2 + t38 ^ 2) + m(5) * (t136 ^ 2 + t137 ^ 2 + t53 ^ 2) + m(4) * (t104 ^ 2 + t191 ^ 2 + t192 ^ 2) + m(3) * (t175 ^ 2 + (t382 + t383) * t262 ^ 2) + (-t382 * t231 - t36 - t391 - t46) * t296 + (t383 * t232 + t37 + t47 + (-t293 * t231 + t296 * t232) * t296 + t390) * t293; (t293 * t334 + t296 * t335) * t292 + m(7) * (t43 * t90 + t44 * t91) + m(6) * (t113 * t50 + t114 * t51) + m(5) * (t123 * t88 + t124 * t89) + m(4) * (t132 * t138 + t133 * t139) + (-t128 + t333) * t295 + t299; m(7) * (t21 * t22 + t43 * t99 + t44 * t98) + m(6) * (t111 * t51 + t112 * t50 + t38 * t42) + m(5) * (t136 * t89 + t137 * t88 + t53 * t63) + m(4) * (t104 * t119 + t132 * t192 + t133 * t191) + (t378 * t47 + t381 * t46) * t292 + t298 + t29 * t381 + (t93 * t293 - t92 * t296) * t380 + t28 * t379; (t293 * t28 + t296 * t29) * t292 + (t128 * t295 + (-t293 * t92 - t296 * t93) * t292 + t337) * t295 + m(5) * (t63 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t119 ^ 2 + t132 ^ 2 + t133 ^ 2) + m(7) * (t22 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t42 ^ 2 + t50 ^ 2 + t51 ^ 2) + t318; m(7) * (t48 * t90 + t49 * t91) + m(6) * (t113 * t64 + t114 * t65) + m(5) * (t123 * t125 + t124 * t126) + t333 * t295 + t299; m(7) * (t21 * t41 + t48 * t99 + t49 * t98) + m(6) * (t111 * t65 + t112 * t64 + t38 * t52) + m(5) * (t110 * t53 + t125 * t137 + t126 * t136) + t298; m(7) * (t22 * t41 + t43 * t48 + t44 * t49) + m(6) * (t42 * t52 + t50 * t64 + t51 * t65) + m(5) * (t110 * t63 + t125 * t88 + t126 * t89) + t300; m(7) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t52 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (t110 ^ 2 + t125 ^ 2 + t126 ^ 2) + t300; -t401 + m(7) * (t82 * t90 + t83 * t91) + m(6) * (t113 * t121 + t114 * t122) + t303; m(7) * (t21 * t54 + t82 * t99 + t83 * t98) + m(6) * (t107 * t38 + t111 * t122 + t112 * t121) + t301; m(7) * (t22 * t54 + t43 * t82 + t44 * t83) + m(6) * (t107 * t42 + t121 * t50 + t122 * t51) + t304; m(7) * (t41 * t54 + t48 * t82 + t49 * t83) + m(6) * (t107 * t52 + t121 * t64 + t122 * t65) + t304; m(7) * (t54 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(6) * (t107 ^ 2 + t121 ^ 2 + t122 ^ 2) + t304; m(7) * (t223 * t91 + t225 * t90); m(7) * (t21 * t368 + t223 * t98 + t225 * t99); m(7) * (t22 * t368 + t223 * t44 + t225 * t43); m(7) * (t223 * t49 + t225 * t48 + t368 * t41); m(7) * (t223 * t83 + t225 * t82 + t368 * t54); m(7) * (t278 ^ 2 * t292 ^ 2 + t223 ^ 2 + t225 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
