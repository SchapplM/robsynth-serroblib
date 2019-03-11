% Calculate joint inertia matrix for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:20
% EndTime: 2019-03-10 01:47:41
% DurationCPUTime: 8.81s
% Computational Cost: add. (38491->673), mult. (65165->912), div. (0->0), fcn. (82872->12), ass. (0->322)
t349 = cos(pkin(6));
t353 = sin(qJ(1));
t355 = cos(qJ(2));
t418 = t353 * t355;
t352 = sin(qJ(2));
t356 = cos(qJ(1));
t419 = t352 * t356;
t331 = t349 * t419 + t418;
t416 = qJ(3) + qJ(4);
t346 = sin(t416);
t380 = cos(t416);
t348 = sin(pkin(6));
t422 = t348 * t356;
t300 = t331 * t380 - t346 * t422;
t417 = t355 * t356;
t420 = t352 * t353;
t330 = -t349 * t417 + t420;
t350 = sin(qJ(5));
t437 = cos(qJ(5));
t262 = t300 * t350 - t330 * t437;
t263 = t300 * t437 + t330 * t350;
t373 = t348 * t380;
t299 = t331 * t346 + t356 * t373;
t451 = rSges(7,3) + qJ(6);
t456 = rSges(7,1) + pkin(5);
t414 = t299 * rSges(7,2) + t262 * t451 + t263 * t456;
t318 = t349 * t346 + t352 * t373;
t423 = t348 * t355;
t297 = t318 * t350 + t423 * t437;
t298 = t318 * t437 - t350 * t423;
t426 = t348 * t352;
t317 = t346 * t426 - t349 * t380;
t211 = Icges(7,5) * t298 + Icges(7,6) * t317 + Icges(7,3) * t297;
t213 = Icges(7,4) * t298 + Icges(7,2) * t317 + Icges(7,6) * t297;
t215 = Icges(7,1) * t298 + Icges(7,4) * t317 + Icges(7,5) * t297;
t107 = t211 * t262 + t213 * t299 + t215 * t263;
t212 = Icges(6,5) * t298 - Icges(6,6) * t297 + Icges(6,3) * t317;
t214 = Icges(6,4) * t298 - Icges(6,2) * t297 + Icges(6,6) * t317;
t216 = Icges(6,1) * t298 - Icges(6,4) * t297 + Icges(6,5) * t317;
t108 = t212 * t299 - t214 * t262 + t216 * t263;
t460 = -t107 - t108;
t333 = -t349 * t420 + t417;
t425 = t348 * t353;
t302 = t333 * t380 + t346 * t425;
t332 = t349 * t418 + t419;
t264 = t302 * t350 - t332 * t437;
t265 = t302 * t437 + t332 * t350;
t301 = t333 * t346 - t353 * t373;
t109 = t211 * t264 + t213 * t301 + t215 * t265;
t110 = t212 * t301 - t214 * t264 + t216 * t265;
t459 = -t109 - t110;
t121 = t297 * t211 + t317 * t213 + t298 * t215;
t122 = t317 * t212 - t297 * t214 + t298 * t216;
t444 = -t121 - t122;
t166 = Icges(7,5) * t263 + Icges(7,6) * t299 + Icges(7,3) * t262;
t170 = Icges(7,4) * t263 + Icges(7,2) * t299 + Icges(7,6) * t262;
t174 = Icges(7,1) * t263 + Icges(7,4) * t299 + Icges(7,5) * t262;
t79 = t166 * t262 + t170 * t299 + t174 * t263;
t167 = Icges(7,5) * t265 + Icges(7,6) * t301 + Icges(7,3) * t264;
t171 = Icges(7,4) * t265 + Icges(7,2) * t301 + Icges(7,6) * t264;
t175 = Icges(7,1) * t265 + Icges(7,4) * t301 + Icges(7,5) * t264;
t80 = t167 * t262 + t171 * t299 + t175 * t263;
t5 = t107 * t317 + t299 * t79 + t301 * t80;
t168 = Icges(6,5) * t263 - Icges(6,6) * t262 + Icges(6,3) * t299;
t172 = Icges(6,4) * t263 - Icges(6,2) * t262 + Icges(6,6) * t299;
t176 = Icges(6,1) * t263 - Icges(6,4) * t262 + Icges(6,5) * t299;
t81 = t168 * t299 - t172 * t262 + t176 * t263;
t169 = Icges(6,5) * t265 - Icges(6,6) * t264 + Icges(6,3) * t301;
t173 = Icges(6,4) * t265 - Icges(6,2) * t264 + Icges(6,6) * t301;
t177 = Icges(6,1) * t265 - Icges(6,4) * t264 + Icges(6,5) * t301;
t82 = t169 * t299 - t173 * t262 + t177 * t263;
t6 = t108 * t317 + t299 * t81 + t301 * t82;
t458 = t6 + t5;
t83 = t166 * t264 + t170 * t301 + t174 * t265;
t84 = t167 * t264 + t171 * t301 + t175 * t265;
t7 = t109 * t317 + t299 * t83 + t301 * t84;
t85 = t168 * t301 - t172 * t264 + t176 * t265;
t86 = t169 * t301 - t173 * t264 + t177 * t265;
t8 = t110 * t317 + t299 * t85 + t301 * t86;
t457 = t7 + t8;
t455 = t460 * t423 + (t80 + t82) * t332 + (t79 + t81) * t330;
t454 = t459 * t423 + (t84 + t86) * t332 + (t83 + t85) * t330;
t116 = t121 * t317;
t95 = t166 * t297 + t170 * t317 + t174 * t298;
t96 = t167 * t297 + t171 * t317 + t175 * t298;
t35 = t95 * t299 + t96 * t301 + t116;
t117 = t122 * t317;
t97 = t168 * t317 - t172 * t297 + t176 * t298;
t98 = t169 * t317 - t173 * t297 + t177 * t298;
t36 = t97 * t299 + t98 * t301 + t117;
t453 = t36 + t35;
t431 = t98 * t332;
t432 = t97 * t330;
t433 = t96 * t332;
t434 = t95 * t330;
t452 = t423 * t444 + t431 + t432 + t433 + t434;
t410 = rSges(7,2) * t317 + t451 * t297 + t298 * t456;
t220 = Icges(5,5) * t300 - Icges(5,6) * t299 + Icges(5,3) * t330;
t222 = Icges(5,4) * t300 - Icges(5,2) * t299 + Icges(5,6) * t330;
t224 = Icges(5,1) * t300 - Icges(5,4) * t299 + Icges(5,5) * t330;
t125 = t220 * t330 - t222 * t299 + t224 * t300;
t221 = Icges(5,5) * t302 - Icges(5,6) * t301 + Icges(5,3) * t332;
t223 = Icges(5,4) * t302 - Icges(5,2) * t301 + Icges(5,6) * t332;
t225 = Icges(5,1) * t302 - Icges(5,4) * t301 + Icges(5,5) * t332;
t126 = t221 * t330 - t223 * t299 + t225 * t300;
t268 = Icges(5,5) * t318 - Icges(5,6) * t317 - Icges(5,3) * t423;
t269 = Icges(5,4) * t318 - Icges(5,2) * t317 - Icges(5,6) * t423;
t270 = Icges(5,1) * t318 - Icges(5,4) * t317 - Icges(5,5) * t423;
t148 = t268 * t330 - t269 * t299 + t270 * t300;
t450 = t125 * t330 + t126 * t332 - t148 * t423 + t455;
t127 = t220 * t332 - t222 * t301 + t224 * t302;
t128 = t221 * t332 - t223 * t301 + t225 * t302;
t149 = t268 * t332 - t269 * t301 + t270 * t302;
t449 = t127 * t330 + t128 * t332 - t149 * t423 + t454;
t29 = t107 * t349 + (t353 * t80 - t356 * t79) * t348;
t30 = t108 * t349 + (t353 * t82 - t356 * t81) * t348;
t448 = t29 + t30 + t148 * t349 + (-t125 * t356 + t126 * t353) * t348;
t31 = t109 * t349 + (t353 * t84 - t356 * t83) * t348;
t32 = t110 * t349 + (t353 * t86 - t356 * t85) * t348;
t447 = t31 + t32 + t149 * t349 + (-t127 * t356 + t128 * t353) * t348;
t403 = -t317 * t269 + t318 * t270;
t154 = -t268 * t423 + t403;
t138 = -t221 * t423 - t223 * t317 + t225 * t318;
t429 = t138 * t332;
t137 = -t220 * t423 - t222 * t317 + t224 * t318;
t430 = t137 * t330;
t446 = -t154 * t423 + t429 + t430 + t452;
t152 = t154 * t349;
t118 = t121 * t349;
t45 = t118 + (t96 * t353 - t95 * t356) * t348;
t119 = t122 * t349;
t46 = t119 + (t98 * t353 - t97 * t356) * t348;
t445 = t45 + t46 + t152 + (-t137 * t356 + t138 * t353) * t348;
t440 = t330 / 0.2e1;
t439 = t332 / 0.2e1;
t438 = t349 / 0.2e1;
t354 = cos(qJ(3));
t345 = pkin(3) * t354 + pkin(2);
t436 = -pkin(2) + t345;
t277 = Icges(3,5) * t331 - Icges(3,6) * t330 - Icges(3,3) * t422;
t428 = t277 * t356;
t357 = -pkin(10) - pkin(9);
t427 = t330 * t357;
t424 = t348 * t354;
t351 = sin(qJ(3));
t421 = t349 * t351;
t179 = t263 * rSges(6,1) - t262 * rSges(6,2) + t299 * rSges(6,3);
t250 = t300 * pkin(4) + t299 * pkin(11);
t231 = t332 * t250;
t415 = t332 * t179 + t231;
t413 = t301 * rSges(7,2) + t264 * t451 + t265 * t456;
t181 = t265 * rSges(6,1) - t264 * rSges(6,2) + t301 * rSges(6,3);
t251 = t302 * pkin(4) + pkin(11) * t301;
t412 = -t181 - t251;
t326 = t330 * pkin(9);
t396 = t351 * t422;
t337 = pkin(3) * t396;
t229 = t331 * t436 - t326 - t337 - t427;
t283 = pkin(3) * t421 + ((pkin(9) + t357) * t355 + t436 * t352) * t348;
t411 = t229 * t423 + t330 * t283;
t218 = rSges(6,1) * t298 - rSges(6,2) * t297 + rSges(6,3) * t317;
t276 = pkin(4) * t318 + pkin(11) * t317;
t409 = -t218 - t276;
t290 = t333 * pkin(2) + pkin(9) * t332;
t397 = t351 * t425;
t386 = pkin(3) * t397 - t332 * t357 + t333 * t345;
t230 = -t290 + t386;
t288 = t349 * t290;
t408 = t349 * t230 + t288;
t227 = t302 * rSges(5,1) - t301 * rSges(5,2) + t332 * rSges(5,3);
t407 = -t227 - t230;
t289 = t331 * pkin(2) + t326;
t406 = -t229 - t289;
t305 = -t331 * t351 - t354 * t422;
t306 = t331 * t354 - t396;
t239 = rSges(4,1) * t306 + rSges(4,2) * t305 + rSges(4,3) * t330;
t405 = -t239 - t289;
t404 = t250 * t423 + t330 * t276;
t372 = -t300 * rSges(5,1) + t299 * rSges(5,2);
t226 = t330 * rSges(5,3) - t372;
t271 = rSges(5,1) * t318 - rSges(5,2) * t317 - rSges(5,3) * t423;
t162 = t226 * t423 + t330 * t271;
t328 = t349 * t354 - t351 * t426;
t329 = t352 * t424 + t421;
t273 = Icges(4,4) * t329 + Icges(4,2) * t328 - Icges(4,6) * t423;
t274 = Icges(4,1) * t329 + Icges(4,4) * t328 - Icges(4,5) * t423;
t402 = t328 * t273 + t329 * t274;
t401 = -t271 - t283;
t400 = t289 * t425 + t290 * t422;
t399 = t356 * pkin(1) + pkin(8) * t425;
t395 = -t154 + t444;
t394 = t332 * t414 + t231;
t393 = -t251 - t413;
t392 = -t230 + t412;
t391 = -t276 - t410;
t390 = -t283 + t409;
t389 = t349 * t251 + t408;
t388 = -t250 + t406;
t307 = -t333 * t351 + t353 * t424;
t308 = t333 * t354 + t397;
t240 = t308 * rSges(4,1) + t307 * rSges(4,2) + t332 * rSges(4,3);
t313 = Icges(3,3) * t349 + (Icges(3,5) * t352 + Icges(3,6) * t355) * t348;
t314 = Icges(3,6) * t349 + (Icges(3,4) * t352 + Icges(3,2) * t355) * t348;
t315 = Icges(3,5) * t349 + (Icges(3,1) * t352 + Icges(3,4) * t355) * t348;
t387 = t349 * t313 + t314 * t423 + t315 * t426;
t285 = t333 * rSges(3,1) - t332 * rSges(3,2) + rSges(3,3) * t425;
t384 = -t423 / 0.2e1;
t234 = Icges(4,5) * t308 + Icges(4,6) * t307 + Icges(4,3) * t332;
t236 = Icges(4,4) * t308 + Icges(4,2) * t307 + Icges(4,6) * t332;
t238 = Icges(4,1) * t308 + Icges(4,4) * t307 + Icges(4,5) * t332;
t144 = -t234 * t423 + t236 * t328 + t238 * t329;
t272 = Icges(4,5) * t329 + Icges(4,6) * t328 - Icges(4,3) * t423;
t151 = t272 * t332 + t273 * t307 + t274 * t308;
t382 = t144 / 0.2e1 + t151 / 0.2e1;
t233 = Icges(4,5) * t306 + Icges(4,6) * t305 + Icges(4,3) * t330;
t235 = Icges(4,4) * t306 + Icges(4,2) * t305 + Icges(4,6) * t330;
t237 = Icges(4,1) * t306 + Icges(4,4) * t305 + Icges(4,5) * t330;
t143 = -t233 * t423 + t235 * t328 + t237 * t329;
t150 = t272 * t330 + t273 * t305 + t274 * t306;
t381 = t150 / 0.2e1 + t143 / 0.2e1;
t379 = -t353 * pkin(1) + pkin(8) * t422;
t275 = rSges(4,1) * t329 + rSges(4,2) * t328 - rSges(4,3) * t423;
t334 = (pkin(2) * t352 - pkin(9) * t355) * t348;
t378 = t348 * (-t275 - t334);
t377 = -t230 + t393;
t376 = t229 * t425 + t230 * t422 + t400;
t375 = -t283 + t391;
t112 = t179 * t423 + t330 * t218 + t404;
t374 = t348 * (-t334 + t401);
t371 = t386 + t399;
t370 = t330 * t450 + t332 * t449;
t369 = t348 * (-t334 + t390);
t368 = t97 / 0.2e1 + t95 / 0.2e1 + t108 / 0.2e1 + t107 / 0.2e1;
t367 = t98 / 0.2e1 + t109 / 0.2e1 + t96 / 0.2e1 + t110 / 0.2e1;
t366 = t250 * t425 + t251 * t422 + t376;
t77 = t330 * t410 + t414 * t423 + t404;
t365 = t348 * (-t334 + t375);
t364 = -t331 * t345 + t337 + t379;
t284 = rSges(3,1) * t331 - rSges(3,2) * t330 - rSges(3,3) * t422;
t363 = t251 + t371;
t362 = t455 * t299 / 0.2e1 + t454 * t301 / 0.2e1 + t452 * t317 / 0.2e1 + t458 * t440 + t457 * t439 + t453 * t384;
t361 = -t423 * t446 + t370;
t360 = t430 / 0.2e1 + t429 / 0.2e1 + t434 / 0.2e1 + t433 / 0.2e1 + t432 / 0.2e1 + t431 / 0.2e1 + (t148 - t460) * t440 + (t149 - t459) * t439;
t359 = -t250 + t364 + t427;
t358 = t448 * t440 + t447 * t439 + t446 * t438 + t449 * t425 / 0.2e1 + t445 * t384 - t450 * t422 / 0.2e1;
t339 = rSges(2,1) * t356 - rSges(2,2) * t353;
t338 = -rSges(2,1) * t353 - rSges(2,2) * t356;
t316 = t349 * rSges(3,3) + (rSges(3,1) * t352 + rSges(3,2) * t355) * t348;
t282 = Icges(3,1) * t333 - Icges(3,4) * t332 + Icges(3,5) * t425;
t281 = Icges(3,1) * t331 - Icges(3,4) * t330 - Icges(3,5) * t422;
t280 = Icges(3,4) * t333 - Icges(3,2) * t332 + Icges(3,6) * t425;
t279 = Icges(3,4) * t331 - Icges(3,2) * t330 - Icges(3,6) * t422;
t278 = Icges(3,5) * t333 - Icges(3,6) * t332 + Icges(3,3) * t425;
t267 = t285 + t399;
t266 = -t284 + t379;
t245 = -t284 * t349 - t316 * t422;
t244 = t285 * t349 - t316 * t425;
t232 = t387 * t349;
t206 = (t284 * t353 + t285 * t356) * t348;
t205 = t313 * t425 - t314 * t332 + t315 * t333;
t204 = -t313 * t422 - t314 * t330 + t315 * t331;
t203 = t332 * t229;
t202 = t332 * t226;
t194 = t290 + t240 + t399;
t193 = t379 + t405;
t187 = -t240 * t423 - t275 * t332;
t186 = t239 * t423 + t275 * t330;
t185 = t349 * t278 + (t280 * t355 + t282 * t352) * t348;
t184 = t349 * t277 + (t279 * t355 + t281 * t352) * t348;
t183 = t371 + t227;
t182 = (-rSges(5,3) + t357) * t330 + t364 + t372;
t163 = -t227 * t423 - t271 * t332;
t159 = -t272 * t423 + t402;
t158 = t159 * t349;
t157 = t239 * t332 - t240 * t330;
t156 = t349 * t405 + t356 * t378;
t155 = t240 * t349 + t353 * t378 + t288;
t153 = -t227 * t330 + t202;
t147 = (t239 * t353 + t240 * t356) * t348 + t400;
t142 = t363 + t181;
t141 = -t179 + t359;
t140 = t181 * t317 - t218 * t301;
t139 = -t179 * t317 + t218 * t299;
t136 = t332 * t401 + t407 * t423;
t135 = t162 + t411;
t134 = t234 * t332 + t236 * t307 + t238 * t308;
t133 = t233 * t332 + t235 * t307 + t237 * t308;
t132 = t234 * t330 + t236 * t305 + t238 * t306;
t131 = t233 * t330 + t235 * t305 + t237 * t306;
t124 = (-t226 + t406) * t349 + t356 * t374;
t123 = t227 * t349 + t353 * t374 + t408;
t120 = t179 * t301 - t181 * t299;
t115 = t363 + t413;
t114 = t359 - t414;
t113 = t332 * t409 + t412 * t423;
t111 = t330 * t407 + t202 + t203;
t102 = (t226 * t353 + t227 * t356) * t348 + t376;
t101 = t330 * t412 + t415;
t100 = -t301 * t410 + t317 * t413;
t99 = t299 * t410 - t317 * t414;
t94 = t332 * t390 + t392 * t423;
t93 = t112 + t411;
t88 = (-t179 + t388) * t349 + t356 * t369;
t87 = t181 * t349 + t353 * t369 + t389;
t78 = t332 * t391 + t393 * t423;
t76 = -t299 * t413 + t301 * t414;
t75 = t330 * t392 + t203 + t415;
t74 = (t179 * t353 + t181 * t356) * t348 + t366;
t73 = t158 + (-t143 * t356 + t144 * t353) * t348;
t72 = t332 * t375 + t377 * t423;
t71 = t77 + t411;
t70 = t330 * t393 + t394;
t69 = (t388 - t414) * t349 + t356 * t365;
t68 = t349 * t413 + t353 * t365 + t389;
t67 = t143 * t330 + t144 * t332 - t159 * t423;
t64 = t151 * t349 + (-t133 * t356 + t134 * t353) * t348;
t63 = t150 * t349 + (-t131 * t356 + t132 * t353) * t348;
t60 = t133 * t330 + t134 * t332 - t151 * t423;
t59 = t131 * t330 + t132 * t332 - t150 * t423;
t52 = t330 * t377 + t203 + t394;
t51 = (t353 * t414 + t356 * t413) * t348 + t366;
t1 = [t387 + m(7) * (t114 ^ 2 + t115 ^ 2) + m(6) * (t141 ^ 2 + t142 ^ 2) + m(5) * (t182 ^ 2 + t183 ^ 2) + m(4) * (t193 ^ 2 + t194 ^ 2) + m(3) * (t266 ^ 2 + t267 ^ 2) + m(2) * (t338 ^ 2 + t339 ^ 2) + (-t268 - t272) * t423 + Icges(2,3) + t402 + t403 - t444; t152 + t158 + t118 + t232 + t119 + m(7) * (t114 * t69 + t115 * t68) + m(6) * (t141 * t88 + t142 * t87) + m(5) * (t123 * t183 + t124 * t182) + m(4) * (t155 * t194 + t156 * t193) + m(3) * (t244 * t267 + t245 * t266) + ((-t137 / 0.2e1 - t184 / 0.2e1 - t148 / 0.2e1 - t204 / 0.2e1 - t368 - t381) * t356 + (t185 / 0.2e1 + t149 / 0.2e1 + t205 / 0.2e1 + t138 / 0.2e1 + t367 + t382) * t353) * t348; (t73 + t232 + t445) * t349 + m(7) * (t51 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t74 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t102 ^ 2 + t123 ^ 2 + t124 ^ 2) + m(4) * (t147 ^ 2 + t155 ^ 2 + t156 ^ 2) + m(3) * (t206 ^ 2 + t244 ^ 2 + t245 ^ 2) + ((-t63 + (-t279 * t330 + t281 * t331 - t348 * t428) * t422 - t448) * t356 + (t64 + ((-t280 * t332 + t282 * t333 + (t278 * t353 - t428) * t348) * t353 + (t278 * t422 + t279 * t332 + t280 * t330 - t281 * t333 - t282 * t331) * t356) * t348 + t447) * t353 + ((-t184 - t204) * t356 + (t185 + t205) * t353) * t349) * t348; t360 + t382 * t332 + t381 * t330 + m(7) * (t114 * t71 + t115 * t72) + m(6) * (t141 * t93 + t142 * t94) + m(5) * (t135 * t182 + t136 * t183) + m(4) * (t186 * t193 + t187 * t194) + (-t159 + t395) * t423; t358 + m(7) * (t52 * t51 + t68 * t72 + t69 * t71) + m(6) * (t74 * t75 + t87 * t94 + t88 * t93) + m(5) * (t102 * t111 + t123 * t136 + t124 * t135) + m(4) * (t147 * t157 + t155 * t187 + t156 * t186) + (t353 * t60 / 0.2e1 - t355 * t73 / 0.2e1 - t356 * t59 / 0.2e1) * t348 + t63 * t440 + t64 * t439 + t67 * t438; t330 * t59 + t332 * t60 + (-t67 - t446) * t423 + m(7) * (t52 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t75 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(5) * (t111 ^ 2 + t135 ^ 2 + t136 ^ 2) + m(4) * (t157 ^ 2 + t186 ^ 2 + t187 ^ 2) + t370; t360 + t395 * t423 + m(7) * (t114 * t77 + t115 * t78) + m(6) * (t112 * t141 + t113 * t142) + m(5) * (t162 * t182 + t163 * t183); t358 + m(7) * (t51 * t70 + t68 * t78 + t69 * t77) + m(6) * (t101 * t74 + t112 * t88 + t113 * t87) + m(5) * (t102 * t153 + t123 * t163 + t124 * t162); m(7) * (t52 * t70 + t71 * t77 + t72 * t78) + m(6) * (t101 * t75 + t112 * t93 + t113 * t94) + m(5) * (t111 * t153 + t135 * t162 + t136 * t163) + t361; m(7) * (t70 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(6) * (t101 ^ 2 + t112 ^ 2 + t113 ^ 2) + m(5) * (t153 ^ 2 + t162 ^ 2 + t163 ^ 2) + t361; t117 + t116 + m(7) * (t100 * t115 + t114 * t99) + m(6) * (t139 * t141 + t140 * t142) + t367 * t301 + t368 * t299; (t35 / 0.2e1 + t36 / 0.2e1) * t349 + (t45 / 0.2e1 + t46 / 0.2e1) * t317 + (t31 / 0.2e1 + t32 / 0.2e1) * t301 + (t29 / 0.2e1 + t30 / 0.2e1) * t299 + m(7) * (t100 * t68 + t51 * t76 + t69 * t99) + m(6) * (t120 * t74 + t139 * t88 + t140 * t87) + ((-t5 / 0.2e1 - t6 / 0.2e1) * t356 + (t8 / 0.2e1 + t7 / 0.2e1) * t353) * t348; m(7) * (t100 * t72 + t52 * t76 + t71 * t99) + m(6) * (t120 * t75 + t139 * t93 + t140 * t94) + t362; m(7) * (t100 * t78 + t70 * t76 + t77 * t99) + m(6) * (t101 * t120 + t112 * t139 + t113 * t140) + t362; t453 * t317 + t457 * t301 + t458 * t299 + m(7) * (t100 ^ 2 + t76 ^ 2 + t99 ^ 2) + m(6) * (t120 ^ 2 + t139 ^ 2 + t140 ^ 2); m(7) * (t114 * t264 + t115 * t262); m(7) * (t262 * t68 + t264 * t69 + t297 * t51); m(7) * (t262 * t72 + t264 * t71 + t297 * t52); m(7) * (t262 * t78 + t264 * t77 + t297 * t70); m(7) * (t100 * t262 + t264 * t99 + t297 * t76); m(7) * (t262 ^ 2 + t264 ^ 2 + t297 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
