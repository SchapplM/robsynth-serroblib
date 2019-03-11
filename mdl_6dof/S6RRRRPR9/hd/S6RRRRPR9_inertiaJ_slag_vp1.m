% Calculate joint inertia matrix for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:48:01
% EndTime: 2019-03-09 22:48:20
% DurationCPUTime: 9.25s
% Computational Cost: add. (36649->666), mult. (56394->894), div. (0->0), fcn. (71213->14), ass. (0->329)
t349 = cos(pkin(6));
t353 = sin(qJ(1));
t355 = cos(qJ(2));
t418 = t353 * t355;
t352 = sin(qJ(2));
t356 = cos(qJ(1));
t419 = t352 * t356;
t325 = t349 * t419 + t418;
t416 = qJ(3) + qJ(4);
t343 = sin(t416);
t377 = cos(t416);
t347 = sin(pkin(6));
t422 = t347 * t356;
t295 = t325 * t377 - t343 * t422;
t417 = t355 * t356;
t420 = t352 * t353;
t324 = -t349 * t417 + t420;
t346 = sin(pkin(12));
t348 = cos(pkin(12));
t253 = -t295 * t346 + t324 * t348;
t429 = t324 * t346;
t254 = t295 * t348 + t429;
t369 = t347 * t377;
t294 = t325 * t343 + t356 * t369;
t165 = Icges(6,5) * t254 + Icges(6,6) * t253 + Icges(6,3) * t294;
t212 = Icges(5,4) * t295 - Icges(5,2) * t294 + Icges(5,6) * t324;
t465 = t165 - t212;
t327 = -t349 * t420 + t417;
t425 = t347 * t353;
t297 = t327 * t377 + t343 * t425;
t326 = t349 * t418 + t419;
t255 = -t297 * t346 + t326 * t348;
t427 = t326 * t346;
t256 = t297 * t348 + t427;
t296 = t327 * t343 - t353 * t369;
t166 = Icges(6,5) * t256 + Icges(6,6) * t255 + Icges(6,3) * t296;
t213 = Icges(5,4) * t297 - Icges(5,2) * t296 + Icges(5,6) * t326;
t464 = t166 - t213;
t314 = t349 * t343 + t352 * t369;
t423 = t347 * t355;
t292 = -t314 * t346 - t348 * t423;
t396 = t346 * t423;
t293 = t314 * t348 - t396;
t426 = t347 * t352;
t313 = t343 * t426 - t349 * t377;
t200 = Icges(6,5) * t293 + Icges(6,6) * t292 + Icges(6,3) * t313;
t201 = Icges(6,4) * t293 + Icges(6,2) * t292 + Icges(6,6) * t313;
t202 = Icges(6,1) * t293 + Icges(6,4) * t292 + Icges(6,5) * t313;
t107 = t313 * t200 + t292 * t201 + t293 * t202;
t259 = Icges(5,5) * t314 - Icges(5,6) * t313 - Icges(5,3) * t423;
t260 = Icges(5,4) * t314 - Icges(5,2) * t313 - Icges(5,6) * t423;
t261 = Icges(5,1) * t314 - Icges(5,4) * t313 - Icges(5,5) * t423;
t404 = -t313 * t260 + t314 * t261;
t141 = -t259 * t423 + t404;
t463 = t107 + t141;
t167 = Icges(6,4) * t254 + Icges(6,2) * t253 + Icges(6,6) * t294;
t169 = Icges(6,1) * t254 + Icges(6,4) * t253 + Icges(6,5) * t294;
t210 = Icges(5,5) * t295 - Icges(5,6) * t294 + Icges(5,3) * t324;
t214 = Icges(5,1) * t295 - Icges(5,4) * t294 + Icges(5,5) * t324;
t462 = t167 * t253 + t169 * t254 + t210 * t324 + t214 * t295 + t465 * t294;
t168 = Icges(6,4) * t256 + Icges(6,2) * t255 + Icges(6,6) * t296;
t170 = Icges(6,1) * t256 + Icges(6,4) * t255 + Icges(6,5) * t296;
t211 = Icges(5,5) * t297 - Icges(5,6) * t296 + Icges(5,3) * t326;
t215 = Icges(5,1) * t297 - Icges(5,4) * t296 + Icges(5,5) * t326;
t461 = t168 * t253 + t170 * t254 + t211 * t324 + t215 * t295 + t464 * t294;
t460 = t167 * t255 + t169 * t256 + t210 * t326 + t214 * t297 + t465 * t296;
t459 = t168 * t255 + t170 * t256 + t211 * t326 + t215 * t297 + t464 * t296;
t135 = t259 * t324 - t260 * t294 + t261 * t295;
t97 = t200 * t294 + t201 * t253 + t202 * t254;
t458 = t97 + t135;
t136 = t259 * t326 - t260 * t296 + t261 * t297;
t98 = t200 * t296 + t201 * t255 + t202 * t256;
t457 = t98 + t136;
t456 = t463 * t349;
t345 = pkin(12) + qJ(6);
t341 = sin(t345);
t342 = cos(t345);
t246 = -t295 * t341 + t324 * t342;
t247 = t295 * t342 + t324 * t341;
t366 = -t247 * rSges(7,1) - t246 * rSges(7,2);
t162 = t294 * rSges(7,3) - t366;
t286 = t294 * qJ(5);
t350 = -pkin(11) - qJ(5);
t339 = pkin(5) * t348 + pkin(4);
t437 = -pkin(4) + t339;
t414 = pkin(5) * t429 - t294 * t350 + t295 * t437 + t162 - t286;
t281 = -t314 * t341 - t342 * t423;
t282 = t314 * t342 - t341 * t423;
t198 = rSges(7,1) * t282 + rSges(7,2) * t281 + rSges(7,3) * t313;
t455 = -pkin(5) * t396 + t437 * t314 + (-qJ(5) - t350) * t313 + t198;
t156 = Icges(7,5) * t247 + Icges(7,6) * t246 + Icges(7,3) * t294;
t158 = Icges(7,4) * t247 + Icges(7,2) * t246 + Icges(7,6) * t294;
t160 = Icges(7,1) * t247 + Icges(7,4) * t246 + Icges(7,5) * t294;
t69 = t156 * t294 + t158 * t246 + t160 * t247;
t248 = -t297 * t341 + t326 * t342;
t249 = t297 * t342 + t326 * t341;
t157 = Icges(7,5) * t249 + Icges(7,6) * t248 + Icges(7,3) * t296;
t159 = Icges(7,4) * t249 + Icges(7,2) * t248 + Icges(7,6) * t296;
t161 = Icges(7,1) * t249 + Icges(7,4) * t248 + Icges(7,5) * t296;
t70 = t157 * t294 + t159 * t246 + t161 * t247;
t195 = Icges(7,5) * t282 + Icges(7,6) * t281 + Icges(7,3) * t313;
t196 = Icges(7,4) * t282 + Icges(7,2) * t281 + Icges(7,6) * t313;
t197 = Icges(7,1) * t282 + Icges(7,4) * t281 + Icges(7,5) * t313;
t92 = t195 * t294 + t196 * t246 + t197 * t247;
t11 = t324 * t69 + t326 * t70 - t423 * t92;
t454 = t324 * t462 + t461 * t326 - t458 * t423 + t11;
t71 = t156 * t296 + t158 * t248 + t160 * t249;
t72 = t157 * t296 + t159 * t248 + t161 * t249;
t93 = t195 * t296 + t196 * t248 + t197 * t249;
t12 = t324 * t71 + t326 * t72 - t423 * t93;
t453 = t324 * t460 + t326 * t459 - t423 * t457 + t12;
t15 = t92 * t349 + (t353 * t70 - t356 * t69) * t347;
t452 = t15 + t458 * t349 + (t461 * t353 - t356 * t462) * t347;
t16 = t93 * t349 + (t353 * t72 - t356 * t71) * t347;
t451 = t16 + t457 * t349 + (t353 * t459 - t356 * t460) * t347;
t104 = t313 * t195 + t281 * t196 + t282 * t197;
t82 = t157 * t313 + t159 * t281 + t161 * t282;
t435 = t82 * t326;
t81 = t156 * t313 + t158 * t281 + t160 * t282;
t436 = t81 * t324;
t31 = -t104 * t423 + t435 + t436;
t127 = -t211 * t423 - t213 * t313 + t215 * t314;
t431 = t127 * t326;
t126 = -t210 * t423 - t212 * t313 + t214 * t314;
t432 = t126 * t324;
t88 = t166 * t313 + t168 * t292 + t170 * t293;
t433 = t88 * t326;
t87 = t165 * t313 + t167 * t292 + t169 * t293;
t434 = t87 * t324;
t450 = -t463 * t423 + t31 + t431 + t432 + t433 + t434;
t101 = t104 * t349;
t33 = t101 + (t82 * t353 - t81 * t356) * t347;
t449 = t33 + t456 + ((-t126 - t87) * t356 + (t127 + t88) * t353) * t347;
t448 = -t104 - t107;
t163 = t249 * rSges(7,1) + t248 * rSges(7,2) + t296 * rSges(7,3);
t447 = pkin(5) * t427 - t296 * t350 + t297 * t339 + t163;
t446 = t294 / 0.2e1;
t445 = t296 / 0.2e1;
t444 = t313 / 0.2e1;
t443 = t324 / 0.2e1;
t442 = t326 / 0.2e1;
t441 = t349 / 0.2e1;
t440 = t353 / 0.2e1;
t439 = -t356 / 0.2e1;
t354 = cos(qJ(3));
t340 = pkin(3) * t354 + pkin(2);
t438 = -pkin(2) + t340;
t268 = Icges(3,5) * t325 - Icges(3,6) * t324 - Icges(3,3) * t422;
t430 = t268 * t356;
t357 = -pkin(10) - pkin(9);
t428 = t324 * t357;
t424 = t347 * t354;
t351 = sin(qJ(3));
t421 = t349 * t351;
t171 = t254 * rSges(6,1) + t253 * rSges(6,2) + t294 * rSges(6,3);
t238 = t295 * pkin(4) + t286;
t220 = t326 * t238;
t415 = t326 * t171 + t220;
t239 = t297 * pkin(4) + qJ(5) * t296;
t413 = -t239 + t447;
t172 = t256 * rSges(6,1) + t255 * rSges(6,2) + t296 * rSges(6,3);
t412 = -t172 - t239;
t203 = rSges(6,1) * t293 + rSges(6,2) * t292 + rSges(6,3) * t313;
t266 = t314 * pkin(4) + t313 * qJ(5);
t411 = -t203 - t266;
t320 = t324 * pkin(9);
t394 = t351 * t422;
t331 = pkin(3) * t394;
t218 = t325 * t438 - t320 - t331 - t428;
t274 = pkin(3) * t421 + ((pkin(9) + t357) * t355 + t438 * t352) * t347;
t410 = t218 * t423 + t324 * t274;
t285 = t327 * pkin(2) + pkin(9) * t326;
t395 = t351 * t425;
t383 = pkin(3) * t395 - t326 * t357 + t327 * t340;
t219 = -t285 + t383;
t280 = t349 * t285;
t409 = t349 * t219 + t280;
t217 = t297 * rSges(5,1) - t296 * rSges(5,2) + t326 * rSges(5,3);
t408 = -t217 - t219;
t284 = t325 * pkin(2) + t320;
t407 = -t218 - t284;
t300 = -t325 * t351 - t354 * t422;
t301 = t325 * t354 - t394;
t228 = rSges(4,1) * t301 + rSges(4,2) * t300 + rSges(4,3) * t324;
t406 = -t228 - t284;
t405 = t238 * t423 + t324 * t266;
t367 = -t295 * rSges(5,1) + t294 * rSges(5,2);
t216 = t324 * rSges(5,3) - t367;
t262 = rSges(5,1) * t314 - rSges(5,2) * t313 - rSges(5,3) * t423;
t154 = t216 * t423 + t324 * t262;
t322 = t349 * t354 - t351 * t426;
t323 = t352 * t424 + t421;
t264 = Icges(4,4) * t323 + Icges(4,2) * t322 - Icges(4,6) * t423;
t265 = Icges(4,1) * t323 + Icges(4,4) * t322 - Icges(4,5) * t423;
t403 = t322 * t264 + t323 * t265;
t402 = -t262 - t274;
t401 = t284 * t425 + t285 * t422;
t400 = t356 * pkin(1) + pkin(8) * t425;
t398 = t92 / 0.2e1 + t81 / 0.2e1;
t397 = t93 / 0.2e1 + t82 / 0.2e1;
t393 = -t141 + t448;
t392 = t414 * t326 + t220;
t391 = -t239 - t413;
t390 = -t219 + t412;
t389 = -t266 - t455;
t388 = -t274 + t411;
t387 = t349 * t239 + t409;
t386 = -t238 + t407;
t302 = -t327 * t351 + t353 * t424;
t303 = t327 * t354 + t395;
t229 = t303 * rSges(4,1) + t302 * rSges(4,2) + t326 * rSges(4,3);
t308 = Icges(3,3) * t349 + (Icges(3,5) * t352 + Icges(3,6) * t355) * t347;
t309 = Icges(3,6) * t349 + (Icges(3,4) * t352 + Icges(3,2) * t355) * t347;
t310 = Icges(3,5) * t349 + (Icges(3,1) * t352 + Icges(3,4) * t355) * t347;
t384 = t349 * t308 + t309 * t423 + t310 * t426;
t276 = t327 * rSges(3,1) - t326 * rSges(3,2) + rSges(3,3) * t425;
t381 = -t423 / 0.2e1;
t222 = Icges(4,5) * t301 + Icges(4,6) * t300 + Icges(4,3) * t324;
t224 = Icges(4,4) * t301 + Icges(4,2) * t300 + Icges(4,6) * t324;
t226 = Icges(4,1) * t301 + Icges(4,4) * t300 + Icges(4,5) * t324;
t130 = -t222 * t423 + t224 * t322 + t226 * t323;
t263 = Icges(4,5) * t323 + Icges(4,6) * t322 - Icges(4,3) * t423;
t137 = t263 * t324 + t264 * t300 + t265 * t301;
t379 = t137 / 0.2e1 + t130 / 0.2e1;
t223 = Icges(4,5) * t303 + Icges(4,6) * t302 + Icges(4,3) * t326;
t225 = Icges(4,4) * t303 + Icges(4,2) * t302 + Icges(4,6) * t326;
t227 = Icges(4,1) * t303 + Icges(4,4) * t302 + Icges(4,5) * t326;
t131 = -t223 * t423 + t225 * t322 + t227 * t323;
t138 = t263 * t326 + t264 * t302 + t265 * t303;
t378 = t138 / 0.2e1 + t131 / 0.2e1;
t376 = -t353 * pkin(1) + pkin(8) * t422;
t267 = rSges(4,1) * t323 + rSges(4,2) * t322 - rSges(4,3) * t423;
t328 = (pkin(2) * t352 - pkin(9) * t355) * t347;
t375 = t347 * (-t267 - t328);
t374 = -t219 + t391;
t373 = -t274 + t389;
t372 = t218 * t425 + t219 * t422 + t401;
t102 = t171 * t423 + t324 * t203 + t405;
t371 = t347 * (-t328 + t402);
t99 = t104 * t313;
t28 = t81 * t294 + t82 * t296 + t99;
t3 = t294 * t69 + t296 * t70 + t313 * t92;
t4 = t294 * t71 + t296 * t72 + t313 * t93;
t370 = t11 * t446 + t12 * t445 + t28 * t381 + t3 * t443 + t31 * t444 + t4 * t442;
t368 = t454 * t324 + t453 * t326;
t365 = t383 + t400;
t364 = t347 * (-t328 + t388);
t363 = t238 * t425 + t239 * t422 + t372;
t67 = t455 * t324 + t414 * t423 + t405;
t362 = t347 * (-t328 + t373);
t361 = -t325 * t340 + t331 + t376;
t275 = rSges(3,1) * t325 - rSges(3,2) * t324 - rSges(3,3) * t422;
t360 = -t423 * t450 + t368;
t359 = t432 / 0.2e1 + t431 / 0.2e1 + t436 / 0.2e1 + t435 / 0.2e1 + t434 / 0.2e1 + t433 / 0.2e1 + (t92 + t458) * t443 + (t93 + t457) * t442;
t358 = t452 * t443 + t451 * t442 + t450 * t441 + t453 * t425 / 0.2e1 + t449 * t381 - t454 * t422 / 0.2e1;
t333 = rSges(2,1) * t356 - rSges(2,2) * t353;
t332 = -rSges(2,1) * t353 - rSges(2,2) * t356;
t311 = t349 * rSges(3,3) + (rSges(3,1) * t352 + rSges(3,2) * t355) * t347;
t273 = Icges(3,1) * t327 - Icges(3,4) * t326 + Icges(3,5) * t425;
t272 = Icges(3,1) * t325 - Icges(3,4) * t324 - Icges(3,5) * t422;
t271 = Icges(3,4) * t327 - Icges(3,2) * t326 + Icges(3,6) * t425;
t270 = Icges(3,4) * t325 - Icges(3,2) * t324 - Icges(3,6) * t422;
t269 = Icges(3,5) * t327 - Icges(3,6) * t326 + Icges(3,3) * t425;
t258 = t276 + t400;
t257 = -t275 + t376;
t234 = -t275 * t349 - t311 * t422;
t233 = t276 * t349 - t311 * t425;
t221 = t384 * t349;
t199 = (t275 * t353 + t276 * t356) * t347;
t194 = t308 * t425 - t309 * t326 + t310 * t327;
t193 = -t308 * t422 - t309 * t324 + t310 * t325;
t192 = t326 * t218;
t191 = t326 * t216;
t184 = t285 + t229 + t400;
t183 = t376 + t406;
t178 = -t229 * t423 - t267 * t326;
t177 = t228 * t423 + t267 * t324;
t176 = t349 * t269 + (t271 * t355 + t273 * t352) * t347;
t175 = t349 * t268 + (t270 * t355 + t272 * t352) * t347;
t174 = t365 + t217;
t173 = (-rSges(5,3) + t357) * t324 + t361 + t367;
t155 = -t217 * t423 - t326 * t262;
t148 = -t263 * t423 + t403;
t146 = t148 * t349;
t144 = t228 * t326 - t229 * t324;
t143 = t349 * t406 + t356 * t375;
t142 = t349 * t229 + t353 * t375 + t280;
t140 = -t217 * t324 + t191;
t134 = (t228 * t353 + t229 * t356) * t347 + t401;
t129 = t365 - t412;
t128 = -t171 - t238 + t361 + t428;
t125 = t326 * t402 + t408 * t423;
t124 = t154 + t410;
t123 = t163 * t313 - t198 * t296;
t122 = -t162 * t313 + t198 * t294;
t121 = t223 * t326 + t225 * t302 + t227 * t303;
t120 = t222 * t326 + t224 * t302 + t226 * t303;
t119 = t223 * t324 + t225 * t300 + t227 * t301;
t118 = t222 * t324 + t224 * t300 + t226 * t301;
t117 = t365 + t447;
t116 = -t295 * t339 + (-pkin(5) * t346 + t357) * t324 + (-rSges(7,3) + t350) * t294 + t361 + t366;
t109 = (-t216 + t407) * t349 + t356 * t371;
t108 = t217 * t349 + t353 * t371 + t409;
t105 = t162 * t296 - t163 * t294;
t103 = t326 * t411 + t412 * t423;
t100 = t324 * t408 + t191 + t192;
t94 = (t216 * t353 + t217 * t356) * t347 + t372;
t89 = t324 * t412 + t415;
t86 = t326 * t388 + t390 * t423;
t85 = t102 + t410;
t80 = (-t171 + t386) * t349 + t356 * t364;
t79 = t172 * t349 + t353 * t364 + t387;
t68 = t326 * t389 + t391 * t423;
t66 = t324 * t390 + t192 + t415;
t65 = (t171 * t353 + t172 * t356) * t347 + t363;
t64 = t146 + (-t130 * t356 + t131 * t353) * t347;
t63 = t130 * t324 + t131 * t326 - t148 * t423;
t60 = t138 * t349 + (-t120 * t356 + t121 * t353) * t347;
t59 = t137 * t349 + (-t118 * t356 + t119 * t353) * t347;
t58 = t326 * t373 + t374 * t423;
t57 = t67 + t410;
t54 = (t386 - t414) * t349 + t356 * t362;
t53 = t349 * t413 + t353 * t362 + t387;
t52 = t120 * t324 + t121 * t326 - t138 * t423;
t51 = t118 * t324 + t119 * t326 - t137 * t423;
t50 = t324 * t391 + t392;
t39 = t324 * t374 + t192 + t392;
t38 = (t353 * t414 + t356 * t413) * t347 + t363;
t1 = [(-t259 - t263) * t423 + t384 + m(7) * (t116 ^ 2 + t117 ^ 2) + m(6) * (t128 ^ 2 + t129 ^ 2) + m(5) * (t173 ^ 2 + t174 ^ 2) + m(4) * (t183 ^ 2 + t184 ^ 2) + m(3) * (t257 ^ 2 + t258 ^ 2) + m(2) * (t332 ^ 2 + t333 ^ 2) + Icges(2,3) + t403 + t404 - t448; t221 + t146 + t101 + m(7) * (t116 * t54 + t117 * t53) + m(6) * (t128 * t80 + t129 * t79) + m(5) * (t108 * t174 + t109 * t173) + m(4) * (t142 * t184 + t143 * t183) + m(3) * (t233 * t258 + t234 * t257) + ((-t175 / 0.2e1 - t97 / 0.2e1 - t135 / 0.2e1 - t193 / 0.2e1 - t126 / 0.2e1 - t87 / 0.2e1 - t379 - t398) * t356 + (t127 / 0.2e1 + t88 / 0.2e1 + t176 / 0.2e1 + t98 / 0.2e1 + t136 / 0.2e1 + t194 / 0.2e1 + t378 + t397) * t353) * t347 + t456; (t64 + t221 + t449) * t349 + m(7) * (t38 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t65 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t108 ^ 2 + t109 ^ 2 + t94 ^ 2) + m(4) * (t134 ^ 2 + t142 ^ 2 + t143 ^ 2) + m(3) * (t199 ^ 2 + t233 ^ 2 + t234 ^ 2) + ((-t59 + (-t270 * t324 + t272 * t325 - t347 * t430) * t422 - t452) * t356 + (t60 + ((-t271 * t326 + t273 * t327 + (t269 * t353 - t430) * t347) * t353 + (t269 * t422 + t270 * t326 + t271 * t324 - t272 * t327 - t273 * t325) * t356) * t347 + t451) * t353 + ((-t175 - t193) * t356 + (t176 + t194) * t353) * t349) * t347; (-t148 + t393) * t423 + t378 * t326 + t379 * t324 + m(7) * (t116 * t57 + t117 * t58) + m(6) * (t128 * t85 + t129 * t86) + m(5) * (t124 * t173 + t125 * t174) + m(4) * (t177 * t183 + t178 * t184) + t359; m(7) * (t38 * t39 + t53 * t58 + t54 * t57) + m(6) * (t65 * t66 + t79 * t86 + t80 * t85) + m(5) * (t100 * t94 + t108 * t125 + t109 * t124) + m(4) * (t134 * t144 + t142 * t178 + t143 * t177) + (t51 * t439 - t355 * t64 / 0.2e1 + t52 * t440) * t347 + t358 + t59 * t443 + t60 * t442 + t63 * t441; t324 * t51 + t326 * t52 + (-t63 - t450) * t423 + m(7) * (t39 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t66 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(5) * (t100 ^ 2 + t124 ^ 2 + t125 ^ 2) + m(4) * (t144 ^ 2 + t177 ^ 2 + t178 ^ 2) + t368; t393 * t423 + m(7) * (t116 * t67 + t117 * t68) + m(6) * (t102 * t128 + t103 * t129) + m(5) * (t154 * t173 + t155 * t174) + t359; t358 + m(7) * (t38 * t50 + t53 * t68 + t54 * t67) + m(6) * (t102 * t80 + t103 * t79 + t65 * t89) + m(5) * (t108 * t155 + t109 * t154 + t140 * t94); m(7) * (t39 * t50 + t57 * t67 + t58 * t68) + m(6) * (t102 * t85 + t103 * t86 + t66 * t89) + m(5) * (t100 * t140 + t124 * t154 + t125 * t155) + t360; m(7) * (t50 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t102 ^ 2 + t103 ^ 2 + t89 ^ 2) + m(5) * (t140 ^ 2 + t154 ^ 2 + t155 ^ 2) + t360; m(7) * (t116 * t296 + t117 * t294) + m(6) * (t128 * t296 + t129 * t294); m(7) * (t294 * t53 + t296 * t54 + t313 * t38) + m(6) * (t294 * t79 + t296 * t80 + t313 * t65); m(7) * (t294 * t58 + t296 * t57 + t313 * t39) + m(6) * (t294 * t86 + t296 * t85 + t313 * t66); m(7) * (t294 * t68 + t296 * t67 + t313 * t50) + m(6) * (t102 * t296 + t103 * t294 + t313 * t89); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t294 ^ 2 + t296 ^ 2 + t313 ^ 2); m(7) * (t116 * t122 + t117 * t123) + t99 + t397 * t296 + t398 * t294; m(7) * (t105 * t38 + t122 * t54 + t123 * t53) + t16 * t445 + t15 * t446 + t33 * t444 + t28 * t441 + (t3 * t439 + t4 * t440) * t347; m(7) * (t105 * t39 + t122 * t57 + t123 * t58) + t370; m(7) * (t105 * t50 + t122 * t67 + t123 * t68) + t370; m(7) * (t105 * t313 + t122 * t296 + t123 * t294); t296 * t4 + t294 * t3 + t313 * t28 + m(7) * (t105 ^ 2 + t122 ^ 2 + t123 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
