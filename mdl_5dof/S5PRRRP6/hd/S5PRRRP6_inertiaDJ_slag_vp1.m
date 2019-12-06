% Calculate time derivative of joint inertia matrix for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:48
% EndTime: 2019-12-05 16:51:16
% DurationCPUTime: 17.14s
% Computational Cost: add. (23438->575), mult. (34801->861), div. (0->0), fcn. (33507->8), ass. (0->300)
t524 = Icges(5,4) - Icges(6,5);
t525 = Icges(5,1) + Icges(6,1);
t514 = Icges(6,4) + Icges(5,5);
t531 = Icges(5,2) + Icges(6,3);
t513 = Icges(5,6) - Icges(6,6);
t530 = Icges(6,2) + Icges(5,3);
t301 = qJ(3) + qJ(4);
t298 = sin(t301);
t529 = t524 * t298;
t299 = cos(t301);
t528 = t524 * t299;
t303 = cos(pkin(8));
t302 = sin(pkin(8));
t307 = cos(qJ(2));
t427 = t302 * t307;
t268 = t298 * t427 + t299 * t303;
t389 = t299 * t427;
t269 = -t298 * t303 + t389;
t305 = sin(qJ(2));
t428 = t302 * t305;
t508 = t531 * t268 - t524 * t269 - t513 * t428;
t424 = t303 * t307;
t270 = t298 * t424 - t302 * t299;
t271 = t298 * t302 + t299 * t424;
t425 = t303 * t305;
t507 = t531 * t270 - t524 * t271 - t513 * t425;
t504 = -t524 * t268 + t525 * t269 + t514 * t428;
t503 = -t524 * t270 + t525 * t271 + t514 * t425;
t527 = t531 * t298 - t528;
t526 = t525 * t299 - t529;
t300 = qJD(3) + qJD(4);
t397 = qJD(2) * t305;
t319 = t300 * t303 + t302 * t397;
t234 = t298 * t319 - t300 * t389;
t431 = t300 * t302;
t387 = t298 * t431;
t235 = -t299 * t319 - t307 * t387;
t396 = qJD(2) * t307;
t378 = t302 * t396;
t521 = t513 * t234 + t514 * t235 + t530 * t378;
t377 = t303 * t397;
t388 = t300 * t424;
t236 = t298 * t377 - t299 * t388 - t387;
t237 = -t298 * t388 + (-t377 + t431) * t299;
t376 = t303 * t396;
t520 = t513 * t236 + t514 * t237 + t530 * t376;
t506 = t513 * t268 - t514 * t269 - t530 * t428;
t505 = -t513 * t270 + t514 * t271 + t530 * t425;
t519 = -t513 * t298 + t514 * t299;
t502 = t519 * t305 - t307 * t530;
t483 = t502 * t307;
t485 = t527 * t305 + t513 * t307;
t501 = t526 * t305 - t514 * t307;
t518 = t507 * t298 + t503 * t299;
t517 = t508 * t298 + t504 * t299;
t516 = t268 * t507 + t269 * t503 + t428 * t505;
t515 = t270 * t508 + t271 * t504 - t425 * t506;
t512 = -t234 * t531 - t235 * t524 - t378 * t513;
t511 = -t236 * t531 - t237 * t524 - t376 * t513;
t510 = t524 * t234 + t235 * t525 + t514 * t378;
t509 = t524 * t236 + t237 * t525 + t514 * t376;
t500 = t520 * t305 + t396 * t505;
t499 = t521 * t305 - t396 * t506;
t498 = t517 * t305 + t506 * t307;
t497 = t518 * t305 - t505 * t307;
t496 = -t485 * t298 - t501 * t299;
t495 = rSges(6,1) + pkin(4);
t494 = t508 * t234 - t504 * t235 - t512 * t268 - t510 * t269 - t499 * t302;
t493 = -t507 * t234 + t503 * t235 + t511 * t268 + t509 * t269 + t500 * t302;
t492 = -t508 * t236 + t504 * t237 + t512 * t270 + t510 * t271 + t499 * t303;
t491 = -t507 * t236 + t503 * t237 + t511 * t270 + t509 * t271 + t500 * t303;
t471 = t508 * t268 + t504 * t269 - t506 * t428;
t470 = t507 * t270 + t503 * t271 + t505 * t425;
t490 = rSges(6,3) + qJ(5);
t489 = t485 * t268 + t501 * t269 + t502 * t428;
t488 = t485 * t270 + t501 * t271 + t502 * t425;
t430 = t300 * t305;
t487 = (-t299 * t531 - t529) * t430 + (t513 * t305 - t307 * t527) * qJD(2);
t486 = (t298 * t525 + t528) * t430 + (-t514 * t305 - t307 * t526) * qJD(2);
t482 = ((t514 * t298 + t513 * t299) * t430 + (-t305 * t530 - t307 * t519) * qJD(2)) * t307;
t481 = t516 * t303;
t480 = t515 * t302;
t479 = (-t517 * qJD(2) + t521) * t307 + ((-t508 * t300 - t510) * t299 + (t504 * t300 - t512) * t298 + t506 * qJD(2)) * t305;
t478 = (t518 * qJD(2) - t520) * t307 + ((t507 * t300 + t509) * t299 + (-t503 * t300 + t511) * t298 + t505 * qJD(2)) * t305;
t477 = t496 * t305 + t483;
t476 = t498 * t302 + t497 * t303;
t475 = (t485 * t234 - t235 * t501 + t487 * t268 + t486 * t269) * t307 + (t493 * t303 + (t482 - t494) * t302) * t305 + (((-t483 + t471) * t302 + t481) * t307 + t489 * t305) * qJD(2);
t474 = (t485 * t236 - t237 * t501 + t487 * t270 + t486 * t271) * t307 + ((t482 + t491) * t303 + t492 * t302) * t305 + (((-t483 + t470) * t303 + t480) * t307 + t488 * t305) * qJD(2);
t473 = t493 * t302 + t494 * t303;
t472 = t491 * t302 - t492 * t303;
t469 = rSges(6,2) * t378 + qJD(5) * t268 - t490 * t234 + t495 * t235;
t355 = rSges(6,1) * t299 + rSges(6,3) * t298;
t468 = (pkin(4) * t396 + qJ(5) * t430) * t299 + (qJ(5) * t396 + (-pkin(4) * t300 + qJD(5)) * t305) * t298 + (-rSges(6,1) * t298 + rSges(6,3) * t299) * t430 + (rSges(6,2) * t305 + t307 * t355) * qJD(2);
t407 = rSges(6,2) * t428 + t490 * t268 + t495 * t269;
t406 = rSges(6,2) * t425 + t490 * t270 + t495 * t271;
t467 = -t307 * rSges(6,2) + (pkin(4) * t299 + qJ(5) * t298 + t355) * t305;
t455 = t303 ^ 2;
t456 = t302 ^ 2;
t463 = t455 + t456;
t466 = ((-t496 * t307 - t476) * qJD(2) + t482) * t307 + (((t485 * t300 - t486) * t299 + (-t300 * t501 - t487) * t298) * t307 - t478 * t303 + t479 * t302 + (t477 + t483) * qJD(2)) * t305;
t462 = qJD(2) * (rSges(3,1) * t305 + rSges(3,2) * t307);
t306 = cos(qJ(3));
t450 = pkin(3) * t306;
t253 = -pkin(7) * t307 + t450 * t305;
t459 = 2 * m(4);
t458 = 2 * m(5);
t457 = 2 * m(6);
t454 = t302 / 0.2e1;
t453 = -t303 / 0.2e1;
t452 = t303 / 0.2e1;
t451 = -t307 / 0.2e1;
t304 = sin(qJ(3));
t423 = t304 * t307;
t287 = t302 * t306 - t303 * t423;
t422 = t306 * t307;
t429 = t302 * t304;
t288 = t303 * t422 + t429;
t219 = Icges(4,5) * t288 + Icges(4,6) * t287 + Icges(4,3) * t425;
t221 = Icges(4,4) * t288 + Icges(4,2) * t287 + Icges(4,6) * t425;
t223 = Icges(4,1) * t288 + Icges(4,4) * t287 + Icges(4,5) * t425;
t285 = -t302 * t423 - t303 * t306;
t426 = t303 * t304;
t286 = t302 * t422 - t426;
t99 = t219 * t428 + t221 * t285 + t223 * t286;
t444 = t303 * t99;
t441 = Icges(4,4) * t304;
t440 = Icges(4,4) * t306;
t218 = Icges(4,5) * t286 + Icges(4,6) * t285 + Icges(4,3) * t428;
t220 = Icges(4,4) * t286 + Icges(4,2) * t285 + Icges(4,6) * t428;
t222 = Icges(4,1) * t286 + Icges(4,4) * t285 + Icges(4,5) * t428;
t100 = t218 * t425 + t220 * t287 + t222 * t288;
t435 = t100 * t302;
t344 = Icges(4,5) * t306 - Icges(4,6) * t304;
t394 = qJD(3) * t305;
t419 = t307 * ((-Icges(4,5) * t304 - Icges(4,6) * t306) * t394 + (Icges(4,3) * t305 + t307 * t344) * qJD(2));
t272 = -Icges(4,3) * t307 + t305 * t344;
t416 = t307 * t272;
t415 = rSges(6,2) * t376 + qJD(5) * t270 - t490 * t236 + t495 * t237;
t140 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t378;
t204 = rSges(5,1) * t269 - rSges(5,2) * t268 + rSges(5,3) * t428;
t414 = t140 * t425 + t204 * t376;
t142 = rSges(5,1) * t237 + rSges(5,2) * t236 + rSges(5,3) * t376;
t395 = qJD(3) * t304;
t392 = pkin(3) * t395;
t309 = -t253 * qJD(2) - t307 * t392;
t393 = qJD(3) * t306;
t391 = pkin(3) * t393;
t168 = t302 * t391 + t303 * t309;
t413 = -t142 - t168;
t167 = t302 * t309 - t303 * t391;
t313 = pkin(7) * t305 + t450 * t307;
t228 = -pkin(3) * t426 + t302 * t313;
t412 = t167 * t425 + t228 * t376;
t411 = t407 * t425;
t410 = t406 * t397;
t356 = rSges(5,1) * t299 - rSges(5,2) * t298;
t189 = (-rSges(5,1) * t298 - rSges(5,2) * t299) * t430 + (rSges(5,3) * t305 + t307 * t356) * qJD(2);
t238 = qJD(2) * t313 - t305 * t392;
t408 = -t189 - t238;
t206 = rSges(5,1) * t271 - rSges(5,2) * t270 + rSges(5,3) * t425;
t229 = pkin(3) * t429 + t303 * t313;
t405 = -t206 - t229;
t404 = t307 * t228 + t253 * t428;
t357 = rSges(4,1) * t306 - rSges(4,2) * t304;
t233 = (-rSges(4,1) * t304 - rSges(4,2) * t306) * t394 + (rSges(4,3) * t305 + t307 * t357) * qJD(2);
t359 = pkin(2) * t307 + pkin(6) * t305;
t292 = t359 * qJD(2);
t403 = -t233 - t292;
t261 = -t307 * rSges(5,3) + t305 * t356;
t147 = t307 * t204 + t261 * t428;
t402 = -t253 - t261;
t295 = t305 * pkin(2) - t307 * pkin(6);
t400 = t463 * qJD(2) * t295;
t275 = -t307 * rSges(4,3) + t305 * t357;
t399 = -t275 - t295;
t398 = t463 * t359;
t390 = t299 * t430;
t386 = -t168 - t415;
t385 = t307 * t140 + t189 * t428 + t261 * t378;
t384 = t307 * t167 + t238 * t428 + t253 * t378;
t383 = -t238 - t468;
t382 = -t292 + t408;
t381 = -t229 - t406;
t380 = -t253 - t467;
t379 = -t295 + t402;
t375 = t304 * t397;
t374 = t306 * t397;
t369 = t406 * t307;
t368 = t405 * t307;
t367 = t407 * t376 + t469 * t425;
t366 = t302 * t167 + t303 * t168 - t400;
t365 = -t292 + t383;
t364 = t302 * t228 + t303 * t229 + t398;
t104 = t407 * t307 + t467 * t428;
t363 = -t295 + t380;
t360 = t381 * t307;
t351 = Icges(4,1) * t306 - t441;
t347 = -Icges(4,2) * t304 + t440;
t336 = -t220 * t304 + t222 * t306;
t106 = -t307 * t218 + t305 * t336;
t335 = -t221 * t304 + t223 * t306;
t107 = -t307 * t219 + t305 * t335;
t341 = t106 * t302 + t107 * t303;
t226 = rSges(4,1) * t286 + rSges(4,2) * t285 + rSges(4,3) * t428;
t227 = rSges(4,1) * t288 + rSges(4,2) * t287 + rSges(4,3) * t425;
t334 = t226 * t303 - t227 * t302;
t273 = -Icges(4,6) * t307 + t305 * t347;
t274 = -Icges(4,5) * t307 + t305 * t351;
t330 = t273 * t304 - t274 * t306;
t327 = t469 * t307 + t467 * t378 + t468 * t428;
t247 = -qJD(3) * t286 + t302 * t375;
t248 = qJD(3) * t285 - t302 * t374;
t152 = Icges(4,5) * t248 + Icges(4,6) * t247 + Icges(4,3) * t378;
t321 = t152 * t305 + t218 * t396;
t249 = -qJD(3) * t288 + t303 * t375;
t250 = qJD(3) * t287 - t303 * t374;
t153 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t376;
t320 = t153 * t305 + t219 * t396;
t315 = qJD(2) * (-Icges(3,5) * t305 - Icges(3,6) * t307);
t314 = t298 * t396 + t390;
t312 = t475 * t428 + t474 * t425 + (t476 * t305 + t477 * t307) * t397 + (-t489 * t307 + (t471 * t302 + t481) * t305) * t378 + (-t488 * t307 + (t470 * t303 + t480) * t305) * t376;
t311 = t466 * t307 + t312;
t310 = t474 * t454 + t475 * t453 + (t478 * t302 + t479 * t303) * t451 + t473 * t428 / 0.2e1 + t472 * t425 / 0.2e1 + (t497 * t302 - t498 * t303) * t397 / 0.2e1 + ((t516 * t302 - t471 * t303) * t302 + (t470 * t302 - t515 * t303) * t303) * t396 / 0.2e1;
t279 = t303 * t315;
t278 = t302 * t315;
t240 = t399 * t303;
t239 = t399 * t302;
t232 = (-Icges(4,1) * t304 - t440) * t394 + (Icges(4,5) * t305 + t307 * t351) * qJD(2);
t231 = (-Icges(4,2) * t306 - t441) * t394 + (Icges(4,6) * t305 + t307 * t347) * qJD(2);
t216 = t463 * t462;
t211 = t229 * t397;
t209 = t228 * t425;
t178 = t206 * t397;
t176 = t204 * t425;
t170 = t403 * t303;
t169 = t403 * t302;
t161 = t379 * t303;
t160 = t379 * t302;
t159 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t376;
t158 = rSges(4,1) * t248 + rSges(4,2) * t247 + rSges(4,3) * t378;
t157 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t376;
t156 = Icges(4,1) * t248 + Icges(4,4) * t247 + Icges(4,5) * t378;
t155 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t376;
t154 = Icges(4,4) * t248 + Icges(4,2) * t247 + Icges(4,6) * t378;
t151 = -t227 * t307 - t275 * t425;
t150 = t226 * t307 + t275 * t428;
t149 = -t305 * t330 - t416;
t148 = -t307 * t206 - t261 * t425;
t146 = t363 * t303;
t145 = t363 * t302;
t122 = t272 * t425 + t273 * t287 + t274 * t288;
t121 = t272 * t428 + t273 * t285 + t274 * t286;
t120 = t334 * t305;
t115 = t382 * t303;
t114 = t382 * t302;
t113 = -t206 * t428 + t176;
t108 = t226 * t302 + t227 * t303 + t398;
t105 = -t425 * t467 - t369;
t103 = t402 * t425 + t368;
t102 = t147 + t404;
t101 = t219 * t425 + t221 * t287 + t223 * t288;
t98 = t218 * t428 + t220 * t285 + t222 * t286;
t97 = t365 * t303;
t96 = t365 * t302;
t91 = t158 * t302 + t159 * t303 - t400;
t90 = -t233 * t425 - t307 * t159 + (t227 * t305 - t275 * t424) * qJD(2);
t89 = t233 * t428 + t307 * t158 + (-t226 * t305 + t275 * t427) * qJD(2);
t80 = t405 * t428 + t176 + t209;
t79 = -t406 * t428 + t411;
t78 = t380 * t425 + t360;
t77 = t104 + t404;
t76 = -t307 * t142 + t178 + (-t189 * t305 - t261 * t396) * t303;
t75 = -t204 * t397 + t385;
t74 = t204 * t302 + t206 * t303 + t364;
t73 = (t158 * t303 - t159 * t302) * t305 + t334 * t396;
t72 = (-t142 * t305 - t206 * t396) * t302 + t414;
t71 = t381 * t428 + t209 + t411;
t70 = t302 * t407 + t406 * t303 + t364;
t69 = t140 * t302 + t142 * t303 + t366;
t62 = t178 + t211 + t413 * t307 + (t305 * t408 + t396 * t402) * t303;
t61 = (-t204 - t228) * t397 + t384 + t385;
t60 = t155 * t287 + t157 * t288 + t221 * t249 + t223 * t250 + t320 * t303;
t59 = t154 * t287 + t156 * t288 + t220 * t249 + t222 * t250 + t321 * t303;
t58 = t155 * t285 + t157 * t286 + t221 * t247 + t223 * t248 + t302 * t320;
t57 = t154 * t285 + t156 * t286 + t220 * t247 + t222 * t248 + t302 * t321;
t54 = -t415 * t307 + (-t305 * t468 - t396 * t467) * t303 + t410;
t53 = -t397 * t407 + t327;
t52 = (qJD(2) * t335 - t153) * t307 + (qJD(2) * t219 - t155 * t304 + t157 * t306 + (-t221 * t306 - t223 * t304) * qJD(3)) * t305;
t51 = (qJD(2) * t336 - t152) * t307 + (qJD(2) * t218 - t154 * t304 + t156 * t306 + (-t220 * t306 - t222 * t304) * qJD(3)) * t305;
t46 = (qJD(2) * t368 + t305 * t413) * t302 + t412 + t414;
t37 = t469 * t302 + t415 * t303 + t366;
t32 = (-qJD(2) * t369 - t305 * t415) * t302 + t367;
t31 = t211 + t386 * t307 + (t305 * t383 + t380 * t396) * t303 + t410;
t30 = (-t228 - t407) * t397 + t327 + t384;
t29 = (qJD(2) * t360 + t305 * t386) * t302 + t367 + t412;
t28 = t302 * t60 - t303 * t59;
t27 = t302 * t58 - t303 * t57;
t16 = -(t231 * t287 + t232 * t288 + t249 * t273 + t250 * t274) * t307 + (t59 * t302 + (t60 - t419) * t303) * t305 + (t122 * t305 + (t435 + (t101 - t416) * t303) * t307) * qJD(2);
t15 = -(t231 * t285 + t232 * t286 + t247 * t273 + t248 * t274) * t307 + (t58 * t303 + (t57 - t419) * t302) * t305 + (t121 * t305 + (t444 + (t98 - t416) * t302) * t307) * qJD(2);
t1 = [0; -m(3) * t216 + m(4) * t91 + m(5) * t69 + m(6) * t37; (t114 * t160 + t115 * t161 + t69 * t74) * t458 + (t145 * t96 + t146 * t97 + t37 * t70) * t457 + (t108 * t91 + t169 * t239 + t170 * t240) * t459 + 0.2e1 * m(3) * (-t216 + t462) * t463 * (rSges(3,1) * t307 - rSges(3,2) * t305) + (-t455 * t278 - t27 - t473) * t303 + (t456 * t279 + t28 + (-t302 * t278 + t303 * t279) * t303 + t472) * t302; m(4) * t73 + m(5) * t46 + m(6) * t29; (t305 * (-t106 * t303 + t107 * t302) / 0.2e1 + ((t99 * t302 - t98 * t303) * t454 + (-t100 * t303 + t101 * t302) * t452) * t307) * qJD(2) + (t27 * t454 + t28 * t452) * t305 + t310 + m(4) * (t108 * t73 + t120 * t91 + t150 * t170 + t151 * t169 + t239 * t90 + t240 * t89) + m(5) * (t102 * t115 + t103 * t114 + t160 * t62 + t161 * t61 + t46 * t74 + t69 * t80) + m(6) * (t145 * t31 + t146 * t30 + t29 * t70 + t37 * t71 + t77 * t97 + t78 * t96) + (t52 * t302 - t51 * t303) * t451 + t15 * t453 + t16 * t454; t312 + (t29 * t71 + t30 * t77 + t31 * t78) * t457 + (t102 * t61 + t103 * t62 + t46 * t80) * t458 + (t120 * t73 + t150 * t89 + t151 * t90) * t459 + t15 * t428 + t16 * t425 + (t341 * t397 + (t101 * t303 + t435) * t376 + (t302 * t98 + t444) * t378) * t305 + (-t149 * t397 - t122 * t376 - t121 * t378 - (t149 * qJD(2) + t51 * t302 + t52 * t303) * t305 + (-t341 * qJD(2) - t419 - (-qJD(2) * t272 + t231 * t304 - t232 * t306 + t273 * t393 + t274 * t395) * t305 - t330 * t396) * t307 + t466) * t307; m(5) * t72 + m(6) * t32; t310 + m(5) * (t113 * t69 + t114 * t148 + t115 * t147 + t160 * t76 + t161 * t75 + t72 * t74) + m(6) * (t104 * t97 + t105 * t96 + t145 * t54 + t146 * t53 + t32 * t70 + t37 * t79); m(6) * (t104 * t30 + t105 * t31 + t29 * t79 + t32 * t71 + t53 * t77 + t54 * t78) + m(5) * (t102 * t75 + t103 * t76 + t113 * t46 + t147 * t61 + t148 * t62 + t72 * t80) + t311; t311 + (t104 * t53 + t105 * t54 + t32 * t79) * t457 + (t113 * t72 + t147 * t75 + t148 * t76) * t458; t314 * m(6); m(6) * (t70 * t390 - t145 * t234 - t146 * t236 + t268 * t96 + t270 * t97 + (t305 * t37 + t396 * t70) * t298); m(6) * (t71 * t390 - t234 * t78 - t236 * t77 + t268 * t31 + t270 * t30 + (t29 * t305 + t396 * t71) * t298); m(6) * (t79 * t390 - t104 * t236 - t105 * t234 + t268 * t54 + t270 * t53 + (t305 * t32 + t396 * t79) * t298); (t298 * t305 * t314 - t268 * t234 - t270 * t236) * t457;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
