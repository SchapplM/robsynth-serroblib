% Calculate joint inertia matrix for
% S6RRRRRP7
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:33
% EndTime: 2019-03-10 01:34:52
% DurationCPUTime: 8.94s
% Computational Cost: add. (39654->684), mult. (66042->922), div. (0->0), fcn. (83724->12), ass. (0->327)
t350 = cos(pkin(6));
t355 = sin(qJ(1));
t358 = cos(qJ(2));
t422 = t355 * t358;
t354 = sin(qJ(2));
t359 = cos(qJ(1));
t423 = t354 * t359;
t331 = t350 * t423 + t422;
t420 = qJ(3) + qJ(4);
t347 = sin(t420);
t382 = cos(t420);
t349 = sin(pkin(6));
t426 = t349 * t359;
t301 = t331 * t382 - t347 * t426;
t421 = t358 * t359;
t424 = t354 * t355;
t330 = -t350 * t421 + t424;
t352 = sin(qJ(5));
t356 = cos(qJ(5));
t261 = -t301 * t352 + t330 * t356;
t433 = t330 * t352;
t262 = t301 * t356 + t433;
t375 = t349 * t382;
t300 = t331 * t347 + t359 * t375;
t463 = rSges(7,3) + qJ(6) + pkin(11);
t466 = -t262 * rSges(7,1) - t261 * rSges(7,2) - t300 * t463;
t319 = t347 * t350 + t354 * t375;
t427 = t349 * t358;
t298 = -t319 * t352 - t356 * t427;
t401 = t352 * t427;
t299 = t319 * t356 - t401;
t430 = t349 * t354;
t318 = t347 * t430 - t350 * t382;
t213 = Icges(7,5) * t299 + Icges(7,6) * t298 + Icges(7,3) * t318;
t215 = Icges(7,4) * t299 + Icges(7,2) * t298 + Icges(7,6) * t318;
t217 = Icges(7,1) * t299 + Icges(7,4) * t298 + Icges(7,5) * t318;
t107 = t213 * t300 + t215 * t261 + t217 * t262;
t214 = Icges(6,5) * t299 + Icges(6,6) * t298 + Icges(6,3) * t318;
t216 = Icges(6,4) * t299 + Icges(6,2) * t298 + Icges(6,6) * t318;
t218 = Icges(6,1) * t299 + Icges(6,4) * t298 + Icges(6,5) * t318;
t108 = t214 * t300 + t216 * t261 + t218 * t262;
t465 = -t107 - t108;
t333 = -t350 * t424 + t421;
t429 = t349 * t355;
t303 = t333 * t382 + t347 * t429;
t332 = t350 * t422 + t423;
t263 = -t303 * t352 + t332 * t356;
t431 = t332 * t352;
t264 = t303 * t356 + t431;
t302 = t333 * t347 - t355 * t375;
t109 = t213 * t302 + t215 * t263 + t217 * t264;
t110 = t214 * t302 + t216 * t263 + t218 * t264;
t464 = -t109 - t110;
t119 = t213 * t318 + t215 * t298 + t217 * t299;
t120 = t214 * t318 + t216 * t298 + t218 * t299;
t450 = -t119 - t120;
t170 = Icges(7,5) * t262 + Icges(7,6) * t261 + Icges(7,3) * t300;
t174 = Icges(7,4) * t262 + Icges(7,2) * t261 + Icges(7,6) * t300;
t178 = Icges(7,1) * t262 + Icges(7,4) * t261 + Icges(7,5) * t300;
t81 = t170 * t300 + t174 * t261 + t178 * t262;
t171 = Icges(7,5) * t264 + Icges(7,6) * t263 + Icges(7,3) * t302;
t175 = Icges(7,4) * t264 + Icges(7,2) * t263 + Icges(7,6) * t302;
t179 = Icges(7,1) * t264 + Icges(7,4) * t263 + Icges(7,5) * t302;
t82 = t171 * t300 + t175 * t261 + t179 * t262;
t5 = t107 * t318 + t300 * t81 + t302 * t82;
t172 = Icges(6,5) * t262 + Icges(6,6) * t261 + Icges(6,3) * t300;
t176 = Icges(6,4) * t262 + Icges(6,2) * t261 + Icges(6,6) * t300;
t180 = Icges(6,1) * t262 + Icges(6,4) * t261 + Icges(6,5) * t300;
t83 = t172 * t300 + t176 * t261 + t180 * t262;
t173 = Icges(6,5) * t264 + Icges(6,6) * t263 + Icges(6,3) * t302;
t177 = Icges(6,4) * t264 + Icges(6,2) * t263 + Icges(6,6) * t302;
t181 = Icges(6,1) * t264 + Icges(6,4) * t263 + Icges(6,5) * t302;
t84 = t173 * t300 + t177 * t261 + t181 * t262;
t6 = t108 * t318 + t300 * t83 + t302 * t84;
t462 = t5 + t6;
t85 = t170 * t302 + t174 * t263 + t178 * t264;
t86 = t171 * t302 + t175 * t263 + t179 * t264;
t7 = t109 * t318 + t300 * t85 + t302 * t86;
t87 = t172 * t302 + t176 * t263 + t180 * t264;
t88 = t173 * t302 + t177 * t263 + t181 * t264;
t8 = t110 * t318 + t300 * t87 + t302 * t88;
t461 = t7 + t8;
t460 = t465 * t427 + (t82 + t84) * t332 + (t81 + t83) * t330;
t459 = t464 * t427 + (t86 + t88) * t332 + (t85 + t87) * t330;
t114 = t119 * t318;
t97 = t170 * t318 + t174 * t298 + t178 * t299;
t98 = t171 * t318 + t175 * t298 + t179 * t299;
t35 = t97 * t300 + t98 * t302 + t114;
t100 = t173 * t318 + t177 * t298 + t181 * t299;
t115 = t120 * t318;
t99 = t172 * t318 + t176 * t298 + t180 * t299;
t36 = t100 * t302 + t99 * t300 + t115;
t458 = t36 + t35;
t437 = t100 * t332;
t438 = t99 * t330;
t439 = t98 * t332;
t440 = t97 * t330;
t457 = t427 * t450 + t437 + t438 + t439 + t440;
t296 = t300 * pkin(11);
t345 = pkin(5) * t356 + pkin(4);
t441 = -pkin(4) + t345;
t418 = pkin(5) * t433 + t301 * t441 - t296 - t466;
t415 = rSges(7,1) * t299 + rSges(7,2) * t298 - pkin(5) * t401 + t441 * t319 + (-pkin(11) + t463) * t318;
t222 = Icges(5,5) * t301 - Icges(5,6) * t300 + Icges(5,3) * t330;
t224 = Icges(5,4) * t301 - Icges(5,2) * t300 + Icges(5,6) * t330;
t226 = Icges(5,1) * t301 - Icges(5,4) * t300 + Icges(5,5) * t330;
t123 = t222 * t330 - t224 * t300 + t226 * t301;
t223 = Icges(5,5) * t303 - Icges(5,6) * t302 + Icges(5,3) * t332;
t225 = Icges(5,4) * t303 - Icges(5,2) * t302 + Icges(5,6) * t332;
t227 = Icges(5,1) * t303 - Icges(5,4) * t302 + Icges(5,5) * t332;
t124 = t223 * t330 - t225 * t300 + t227 * t301;
t267 = Icges(5,5) * t319 - Icges(5,6) * t318 - Icges(5,3) * t427;
t268 = Icges(5,4) * t319 - Icges(5,2) * t318 - Icges(5,6) * t427;
t269 = Icges(5,1) * t319 - Icges(5,4) * t318 - Icges(5,5) * t427;
t148 = t267 * t330 - t268 * t300 + t269 * t301;
t456 = t123 * t330 + t124 * t332 - t148 * t427 + t460;
t125 = t222 * t332 - t224 * t302 + t226 * t303;
t126 = t223 * t332 - t225 * t302 + t227 * t303;
t149 = t267 * t332 - t268 * t302 + t269 * t303;
t455 = t125 * t330 + t126 * t332 - t149 * t427 + t459;
t29 = t107 * t350 + (t355 * t82 - t359 * t81) * t349;
t30 = t108 * t350 + (t355 * t84 - t359 * t83) * t349;
t454 = t29 + t30 + t148 * t350 + (-t123 * t359 + t124 * t355) * t349;
t31 = t109 * t350 + (t355 * t86 - t359 * t85) * t349;
t32 = t110 * t350 + (t355 * t88 - t359 * t87) * t349;
t453 = t31 + t32 + t149 * t350 + (-t125 * t359 + t126 * t355) * t349;
t407 = -t268 * t318 + t269 * t319;
t154 = -t267 * t427 + t407;
t136 = -t223 * t427 - t225 * t318 + t227 * t319;
t435 = t136 * t332;
t135 = -t222 * t427 - t224 * t318 + t226 * t319;
t436 = t135 * t330;
t452 = -t154 * t427 + t435 + t436 + t457;
t152 = t154 * t350;
t116 = t119 * t350;
t45 = t116 + (t98 * t355 - t97 * t359) * t349;
t117 = t120 * t350;
t46 = t117 + (t100 * t355 - t99 * t359) * t349;
t451 = t45 + t46 + t152 + (-t135 * t359 + t136 * t355) * t349;
t449 = t264 * rSges(7,1) + t263 * rSges(7,2) + pkin(5) * t431 + t302 * t463 + t303 * t345;
t445 = t330 / 0.2e1;
t444 = t332 / 0.2e1;
t443 = t350 / 0.2e1;
t357 = cos(qJ(3));
t346 = pkin(3) * t357 + pkin(2);
t442 = -pkin(2) + t346;
t276 = Icges(3,5) * t331 - Icges(3,6) * t330 - Icges(3,3) * t426;
t434 = t276 * t359;
t360 = -pkin(10) - pkin(9);
t432 = t330 * t360;
t428 = t349 * t357;
t353 = sin(qJ(3));
t425 = t350 * t353;
t183 = rSges(6,1) * t262 + rSges(6,2) * t261 + rSges(6,3) * t300;
t250 = pkin(4) * t301 + t296;
t232 = t332 * t250;
t419 = t183 * t332 + t232;
t251 = pkin(4) * t303 + pkin(11) * t302;
t417 = -t251 + t449;
t185 = rSges(6,1) * t264 + rSges(6,2) * t263 + rSges(6,3) * t302;
t416 = -t185 - t251;
t326 = t330 * pkin(9);
t399 = t353 * t426;
t337 = pkin(3) * t399;
t230 = t331 * t442 - t326 - t337 - t432;
t282 = pkin(3) * t425 + ((pkin(9) + t360) * t358 + t442 * t354) * t349;
t414 = t230 * t427 + t282 * t330;
t220 = rSges(6,1) * t299 + rSges(6,2) * t298 + rSges(6,3) * t318;
t275 = pkin(4) * t319 + pkin(11) * t318;
t413 = -t220 - t275;
t291 = pkin(2) * t333 + pkin(9) * t332;
t400 = t353 * t429;
t388 = pkin(3) * t400 - t332 * t360 + t333 * t346;
t231 = -t291 + t388;
t288 = t350 * t291;
t412 = t231 * t350 + t288;
t229 = rSges(5,1) * t303 - rSges(5,2) * t302 + rSges(5,3) * t332;
t411 = -t229 - t231;
t290 = t331 * pkin(2) + t326;
t410 = -t230 - t290;
t306 = -t331 * t353 - t357 * t426;
t307 = t331 * t357 - t399;
t240 = rSges(4,1) * t307 + rSges(4,2) * t306 + rSges(4,3) * t330;
t409 = -t240 - t290;
t408 = t250 * t427 + t275 * t330;
t374 = -t301 * rSges(5,1) + t300 * rSges(5,2);
t228 = rSges(5,3) * t330 - t374;
t270 = rSges(5,1) * t319 - rSges(5,2) * t318 - rSges(5,3) * t427;
t166 = t228 * t427 + t270 * t330;
t328 = t350 * t357 - t353 * t430;
t329 = t354 * t428 + t425;
t272 = Icges(4,4) * t329 + Icges(4,2) * t328 - Icges(4,6) * t427;
t273 = Icges(4,1) * t329 + Icges(4,4) * t328 - Icges(4,5) * t427;
t406 = t272 * t328 + t273 * t329;
t405 = -t270 - t282;
t404 = t290 * t429 + t291 * t426;
t403 = pkin(1) * t359 + pkin(8) * t429;
t398 = -t154 + t450;
t397 = t332 * t418 + t232;
t396 = -t251 - t417;
t395 = -t231 + t416;
t394 = -t275 - t415;
t393 = -t282 + t413;
t392 = t251 * t350 + t412;
t391 = -t250 + t410;
t308 = -t333 * t353 + t355 * t428;
t309 = t333 * t357 + t400;
t241 = rSges(4,1) * t309 + rSges(4,2) * t308 + rSges(4,3) * t332;
t314 = Icges(3,3) * t350 + (Icges(3,5) * t354 + Icges(3,6) * t358) * t349;
t315 = Icges(3,6) * t350 + (Icges(3,4) * t354 + Icges(3,2) * t358) * t349;
t316 = Icges(3,5) * t350 + (Icges(3,1) * t354 + Icges(3,4) * t358) * t349;
t389 = t314 * t350 + t315 * t427 + t316 * t430;
t284 = rSges(3,1) * t333 - rSges(3,2) * t332 + rSges(3,3) * t429;
t386 = -t427 / 0.2e1;
t235 = Icges(4,5) * t309 + Icges(4,6) * t308 + Icges(4,3) * t332;
t237 = Icges(4,4) * t309 + Icges(4,2) * t308 + Icges(4,6) * t332;
t239 = Icges(4,1) * t309 + Icges(4,4) * t308 + Icges(4,5) * t332;
t144 = -t235 * t427 + t237 * t328 + t239 * t329;
t271 = Icges(4,5) * t329 + Icges(4,6) * t328 - Icges(4,3) * t427;
t151 = t271 * t332 + t272 * t308 + t273 * t309;
t384 = t144 / 0.2e1 + t151 / 0.2e1;
t234 = Icges(4,5) * t307 + Icges(4,6) * t306 + Icges(4,3) * t330;
t236 = Icges(4,4) * t307 + Icges(4,2) * t306 + Icges(4,6) * t330;
t238 = Icges(4,1) * t307 + Icges(4,4) * t306 + Icges(4,5) * t330;
t143 = -t234 * t427 + t236 * t328 + t238 * t329;
t150 = t271 * t330 + t272 * t306 + t273 * t307;
t383 = t150 / 0.2e1 + t143 / 0.2e1;
t381 = -t355 * pkin(1) + pkin(8) * t426;
t274 = rSges(4,1) * t329 + rSges(4,2) * t328 - rSges(4,3) * t427;
t334 = (pkin(2) * t354 - pkin(9) * t358) * t349;
t380 = t349 * (-t274 - t334);
t379 = -t231 + t396;
t378 = -t282 + t394;
t377 = t230 * t429 + t231 * t426 + t404;
t112 = t183 * t427 + t220 * t330 + t408;
t376 = t349 * (-t334 + t405);
t372 = t388 + t403;
t371 = t330 * t456 + t332 * t455;
t370 = t349 * (-t334 + t393);
t369 = t107 / 0.2e1 + t108 / 0.2e1 + t99 / 0.2e1 + t97 / 0.2e1;
t368 = t250 * t429 + t251 * t426 + t377;
t77 = t330 * t415 + t418 * t427 + t408;
t367 = t110 / 0.2e1 + t109 / 0.2e1 + t100 / 0.2e1 + t98 / 0.2e1;
t366 = t349 * (-t334 + t378);
t365 = -t331 * t346 + t337 + t381;
t283 = rSges(3,1) * t331 - rSges(3,2) * t330 - rSges(3,3) * t426;
t364 = t460 * t300 / 0.2e1 + t459 * t302 / 0.2e1 + t457 * t318 / 0.2e1 + t462 * t445 + t461 * t444 + t458 * t386;
t363 = -t427 * t452 + t371;
t362 = t436 / 0.2e1 + t435 / 0.2e1 + t440 / 0.2e1 + t439 / 0.2e1 + t438 / 0.2e1 + t437 / 0.2e1 + (t148 - t465) * t445 + (t149 - t464) * t444;
t361 = t454 * t445 + t453 * t444 + t452 * t443 + t455 * t429 / 0.2e1 + t451 * t386 - t456 * t426 / 0.2e1;
t339 = rSges(2,1) * t359 - rSges(2,2) * t355;
t338 = -rSges(2,1) * t355 - rSges(2,2) * t359;
t317 = rSges(3,3) * t350 + (rSges(3,1) * t354 + rSges(3,2) * t358) * t349;
t281 = Icges(3,1) * t333 - Icges(3,4) * t332 + Icges(3,5) * t429;
t280 = Icges(3,1) * t331 - Icges(3,4) * t330 - Icges(3,5) * t426;
t279 = Icges(3,4) * t333 - Icges(3,2) * t332 + Icges(3,6) * t429;
t278 = Icges(3,4) * t331 - Icges(3,2) * t330 - Icges(3,6) * t426;
t277 = Icges(3,5) * t333 - Icges(3,6) * t332 + Icges(3,3) * t429;
t266 = t284 + t403;
t265 = -t283 + t381;
t246 = -t283 * t350 - t317 * t426;
t245 = t284 * t350 - t317 * t429;
t233 = t389 * t350;
t208 = (t283 * t355 + t284 * t359) * t349;
t207 = t314 * t429 - t315 * t332 + t316 * t333;
t206 = -t314 * t426 - t315 * t330 + t316 * t331;
t205 = t332 * t230;
t204 = t332 * t228;
t197 = t291 + t241 + t403;
t196 = t381 + t409;
t191 = -t241 * t427 - t274 * t332;
t190 = t240 * t427 + t274 * t330;
t189 = t277 * t350 + (t279 * t358 + t281 * t354) * t349;
t188 = t276 * t350 + (t278 * t358 + t280 * t354) * t349;
t187 = t372 + t229;
t186 = (-rSges(5,3) + t360) * t330 + t365 + t374;
t167 = -t229 * t427 - t270 * t332;
t160 = -t271 * t427 + t406;
t159 = t160 * t350;
t157 = t240 * t332 - t241 * t330;
t156 = t350 * t409 + t359 * t380;
t155 = t241 * t350 + t355 * t380 + t288;
t153 = -t229 * t330 + t204;
t147 = (t240 * t355 + t241 * t359) * t349 + t404;
t142 = t372 - t416;
t141 = -t183 - t250 + t365 + t432;
t140 = t185 * t318 - t220 * t302;
t139 = -t183 * t318 + t220 * t300;
t138 = t372 + t449;
t137 = -t301 * t345 + (-pkin(5) * t352 + t360) * t330 + t365 + t466;
t134 = t332 * t405 + t411 * t427;
t133 = t166 + t414;
t132 = t235 * t332 + t237 * t308 + t239 * t309;
t131 = t234 * t332 + t236 * t308 + t238 * t309;
t130 = t235 * t330 + t237 * t306 + t239 * t307;
t129 = t234 * t330 + t236 * t306 + t238 * t307;
t122 = (-t228 + t410) * t350 + t359 * t376;
t121 = t229 * t350 + t355 * t376 + t412;
t118 = t183 * t302 - t185 * t300;
t113 = t332 * t413 + t416 * t427;
t111 = t330 * t411 + t204 + t205;
t102 = (t228 * t355 + t229 * t359) * t349 + t377;
t101 = t330 * t416 + t419;
t96 = t332 * t393 + t395 * t427;
t95 = t112 + t414;
t90 = (-t183 + t391) * t350 + t359 * t370;
t89 = t185 * t350 + t355 * t370 + t392;
t80 = -t302 * t415 + t318 * t417;
t79 = t300 * t415 - t318 * t418;
t78 = t332 * t394 + t396 * t427;
t76 = t330 * t395 + t205 + t419;
t75 = (t183 * t355 + t185 * t359) * t349 + t368;
t74 = t159 + (-t143 * t359 + t144 * t355) * t349;
t73 = -t300 * t417 + t302 * t418;
t72 = t143 * t330 + t144 * t332 - t160 * t427;
t70 = t332 * t378 + t379 * t427;
t69 = t77 + t414;
t67 = (t391 - t418) * t350 + t359 * t366;
t66 = t350 * t417 + t355 * t366 + t392;
t65 = t151 * t350 + (-t131 * t359 + t132 * t355) * t349;
t64 = t150 * t350 + (-t129 * t359 + t130 * t355) * t349;
t62 = t330 * t396 + t397;
t60 = t131 * t330 + t132 * t332 - t151 * t427;
t59 = t129 * t330 + t130 * t332 - t150 * t427;
t48 = t330 * t379 + t205 + t397;
t47 = (t355 * t418 + t359 * t417) * t349 + t368;
t1 = [m(7) * (t137 ^ 2 + t138 ^ 2) + m(6) * (t141 ^ 2 + t142 ^ 2) + m(5) * (t186 ^ 2 + t187 ^ 2) + m(4) * (t196 ^ 2 + t197 ^ 2) + m(3) * (t265 ^ 2 + t266 ^ 2) + m(2) * (t338 ^ 2 + t339 ^ 2) + (-t267 - t271) * t427 + t389 + Icges(2,3) + t406 + t407 - t450; t233 + t116 + t152 + t159 + t117 + m(5) * (t121 * t187 + t122 * t186) + m(6) * (t141 * t90 + t142 * t89) + m(7) * (t137 * t67 + t138 * t66) + m(4) * (t155 * t197 + t156 * t196) + m(3) * (t245 * t266 + t246 * t265) + ((-t135 / 0.2e1 - t188 / 0.2e1 - t148 / 0.2e1 - t206 / 0.2e1 - t369 - t383) * t359 + (t136 / 0.2e1 + t189 / 0.2e1 + t149 / 0.2e1 + t207 / 0.2e1 + t367 + t384) * t355) * t349; (t74 + t233 + t451) * t350 + m(7) * (t47 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(5) * (t102 ^ 2 + t121 ^ 2 + t122 ^ 2) + m(4) * (t147 ^ 2 + t155 ^ 2 + t156 ^ 2) + m(3) * (t208 ^ 2 + t245 ^ 2 + t246 ^ 2) + ((-t64 + (-t278 * t330 + t280 * t331 - t349 * t434) * t426 - t454) * t359 + (t65 + ((-t279 * t332 + t281 * t333 + (t277 * t355 - t434) * t349) * t355 + (t277 * t426 + t278 * t332 + t279 * t330 - t280 * t333 - t281 * t331) * t359) * t349 + t453) * t355 + ((-t188 - t206) * t359 + (t189 + t207) * t355) * t350) * t349; m(7) * (t137 * t69 + t138 * t70) + m(6) * (t141 * t95 + t142 * t96) + m(5) * (t133 * t186 + t134 * t187) + m(4) * (t190 * t196 + t191 * t197) + t384 * t332 + t383 * t330 + (-t160 + t398) * t427 + t362; t72 * t443 + t65 * t444 + t361 + t64 * t445 + (t355 * t60 / 0.2e1 - t358 * t74 / 0.2e1 - t359 * t59 / 0.2e1) * t349 + m(4) * (t147 * t157 + t155 * t191 + t156 * t190) + m(5) * (t102 * t111 + t121 * t134 + t122 * t133) + m(7) * (t47 * t48 + t66 * t70 + t67 * t69) + m(6) * (t75 * t76 + t89 * t96 + t90 * t95); t330 * t59 + t332 * t60 + (-t72 - t452) * t427 + m(7) * (t48 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t76 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(5) * (t111 ^ 2 + t133 ^ 2 + t134 ^ 2) + m(4) * (t157 ^ 2 + t190 ^ 2 + t191 ^ 2) + t371; t398 * t427 + m(5) * (t166 * t186 + t167 * t187) + m(7) * (t137 * t77 + t138 * t78) + m(6) * (t112 * t141 + t113 * t142) + t362; m(7) * (t47 * t62 + t66 * t78 + t67 * t77) + m(6) * (t101 * t75 + t112 * t90 + t113 * t89) + m(5) * (t102 * t153 + t121 * t167 + t122 * t166) + t361; m(7) * (t48 * t62 + t69 * t77 + t70 * t78) + m(6) * (t101 * t76 + t112 * t95 + t113 * t96) + m(5) * (t111 * t153 + t133 * t166 + t134 * t167) + t363; m(7) * (t62 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(6) * (t101 ^ 2 + t112 ^ 2 + t113 ^ 2) + m(5) * (t153 ^ 2 + t166 ^ 2 + t167 ^ 2) + t363; t115 + t114 + m(7) * (t137 * t79 + t138 * t80) + m(6) * (t139 * t141 + t140 * t142) + t367 * t302 + t369 * t300; (t36 / 0.2e1 + t35 / 0.2e1) * t350 + (t46 / 0.2e1 + t45 / 0.2e1) * t318 + (t32 / 0.2e1 + t31 / 0.2e1) * t302 + (t29 / 0.2e1 + t30 / 0.2e1) * t300 + m(7) * (t47 * t73 + t66 * t80 + t67 * t79) + m(6) * (t118 * t75 + t139 * t90 + t140 * t89) + ((-t5 / 0.2e1 - t6 / 0.2e1) * t359 + (t7 / 0.2e1 + t8 / 0.2e1) * t355) * t349; m(7) * (t48 * t73 + t69 * t79 + t70 * t80) + m(6) * (t118 * t76 + t139 * t95 + t140 * t96) + t364; m(7) * (t62 * t73 + t77 * t79 + t78 * t80) + m(6) * (t101 * t118 + t112 * t139 + t113 * t140) + t364; t458 * t318 + t461 * t302 + t462 * t300 + m(7) * (t73 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t118 ^ 2 + t139 ^ 2 + t140 ^ 2); m(7) * (t137 * t302 + t138 * t300); m(7) * (t300 * t66 + t302 * t67 + t318 * t47); m(7) * (t300 * t70 + t302 * t69 + t318 * t48); m(7) * (t300 * t78 + t302 * t77 + t318 * t62); m(7) * (t300 * t80 + t302 * t79 + t318 * t73); m(7) * (t300 ^ 2 + t302 ^ 2 + t318 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
