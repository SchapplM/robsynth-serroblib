% Calculate time derivative of joint inertia matrix for
% S5PRRRP5
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:47
% EndTime: 2019-12-05 16:48:18
% DurationCPUTime: 17.83s
% Computational Cost: add. (24227->574), mult. (35990->872), div. (0->0), fcn. (34034->8), ass. (0->304)
t530 = Icges(5,4) + Icges(6,4);
t531 = Icges(5,1) + Icges(6,1);
t520 = Icges(5,5) + Icges(6,5);
t529 = Icges(5,2) + Icges(6,2);
t519 = Icges(5,6) + Icges(6,6);
t536 = Icges(5,3) + Icges(6,3);
t305 = qJ(3) + qJ(4);
t299 = cos(t305);
t535 = t530 * t299;
t298 = sin(t305);
t534 = t530 * t298;
t307 = cos(pkin(8));
t306 = sin(pkin(8));
t311 = cos(qJ(2));
t436 = t306 * t311;
t269 = -t298 * t436 - t299 * t307;
t270 = -t298 * t307 + t299 * t436;
t309 = sin(qJ(2));
t437 = t306 * t309;
t512 = t529 * t269 + t530 * t270 + t519 * t437;
t433 = t307 * t311;
t271 = -t298 * t433 + t299 * t306;
t272 = t298 * t306 + t299 * t433;
t434 = t307 * t309;
t511 = t529 * t271 + t530 * t272 + t519 * t434;
t510 = t530 * t269 + t531 * t270 + t520 * t437;
t509 = t530 * t271 + t531 * t272 + t520 * t434;
t533 = -t529 * t298 + t535;
t532 = t531 * t299 - t534;
t304 = qJD(3) + qJD(4);
t404 = qJD(2) * t309;
t324 = t304 * t307 + t306 * t404;
t396 = t304 * t436;
t236 = t298 * t324 - t299 * t396;
t237 = -t298 * t396 - t299 * t324;
t403 = qJD(2) * t311;
t386 = t306 * t403;
t527 = t519 * t236 + t520 * t237 + t536 * t386;
t325 = -t304 * t306 + t307 * t404;
t395 = t304 * t433;
t238 = t298 * t325 - t299 * t395;
t239 = -t298 * t395 - t299 * t325;
t385 = t307 * t403;
t526 = t519 * t238 + t520 * t239 + t536 * t385;
t514 = t519 * t269 + t520 * t270 + t536 * t437;
t513 = t519 * t271 + t520 * t272 + t536 * t434;
t525 = -t519 * t298 + t520 * t299;
t508 = t525 * t309 - t311 * t536;
t490 = t508 * t311;
t507 = t533 * t309 - t519 * t311;
t506 = t532 * t309 - t520 * t311;
t524 = -t298 * t511 + t299 * t509;
t523 = -t298 * t512 + t510 * t299;
t522 = t511 * t269 + t509 * t270 + t513 * t437;
t521 = t512 * t271 + t510 * t272 + t514 * t434;
t518 = t529 * t236 + t530 * t237 + t519 * t386;
t517 = t529 * t238 + t530 * t239 + t519 * t385;
t516 = t530 * t236 + t531 * t237 + t520 * t386;
t515 = t530 * t238 + t531 * t239 + t520 * t385;
t505 = t526 * t309 + t513 * t403;
t504 = t527 * t309 + t514 * t403;
t503 = t523 * t309 - t514 * t311;
t502 = t524 * t309 - t513 * t311;
t501 = t507 * t298 - t506 * t299;
t500 = -t512 * t236 - t510 * t237 - t518 * t269 - t516 * t270 - t504 * t306;
t499 = t511 * t236 + t509 * t237 + t517 * t269 + t515 * t270 + t505 * t306;
t498 = t512 * t238 + t510 * t239 + t518 * t271 + t516 * t272 + t504 * t307;
t497 = t511 * t238 + t509 * t239 + t517 * t271 + t515 * t272 + t505 * t307;
t478 = t512 * t269 + t510 * t270 + t514 * t437;
t477 = t511 * t271 + t509 * t272 + t513 * t434;
t496 = t507 * t269 + t506 * t270 + t508 * t437;
t495 = t507 * t271 + t506 * t272 + t508 * t434;
t439 = t304 * t309;
t494 = (t529 * t299 + t534) * t439 + (-t519 * t309 - t311 * t533) * qJD(2);
t493 = (-t531 * t298 - t535) * t439 + (t520 * t309 + t311 * t532) * qJD(2);
t489 = ((t520 * t298 + t519 * t299) * t439 + (-t309 * t536 - t525 * t311) * qJD(2)) * t311;
t488 = t522 * t307;
t487 = t521 * t306;
t486 = (-t523 * qJD(2) + t527) * t311 + ((t512 * t304 - t516) * t299 + (t510 * t304 + t518) * t298 - t514 * qJD(2)) * t309;
t485 = (t524 * qJD(2) - t526) * t311 + ((-t511 * t304 + t515) * t299 + (-t509 * t304 - t517) * t298 + t513 * qJD(2)) * t309;
t484 = t501 * t309 + t490;
t483 = t503 * t306 + t502 * t307;
t482 = (-t236 * t507 - t237 * t506 + t494 * t269 - t493 * t270) * t311 + (t499 * t307 + (t489 - t500) * t306) * t309 + (((-t490 + t478) * t306 + t488) * t311 + t496 * t309) * qJD(2);
t481 = (-t238 * t507 - t239 * t506 + t494 * t271 - t493 * t272) * t311 + ((t489 + t497) * t307 + t498 * t306) * t309 + (((-t490 + t477) * t307 + t487) * t311 + t495 * t309) * qJD(2);
t480 = t499 * t306 + t500 * t307;
t479 = t497 * t306 - t498 * t307;
t406 = pkin(4) * t299;
t241 = -qJ(5) * t311 + t309 * t406;
t441 = t298 * t304;
t368 = pkin(4) * t441;
t313 = -qJD(2) * t241 + qJD(5) * t309 - t368 * t311;
t440 = t299 * t304;
t367 = pkin(4) * t440;
t476 = rSges(6,1) * t237 + rSges(6,2) * t236 + rSges(6,3) * t386 + t306 * t313 - t307 * t367;
t317 = qJ(5) * t309 + t311 * t406;
t361 = rSges(6,1) * t299 - rSges(6,2) * t298;
t475 = -qJD(5) * t311 - t309 * t368 + (-rSges(6,1) * t298 - rSges(6,2) * t299) * t439 + (rSges(6,3) * t309 + t311 * t361 + t317) * qJD(2);
t378 = pkin(4) * t298;
t418 = rSges(6,1) * t270 + rSges(6,2) * t269 + rSges(6,3) * t437 + t306 * t317 - t307 * t378;
t417 = rSges(6,1) * t272 + rSges(6,2) * t271 + rSges(6,3) * t434 + t306 * t378 + t307 * t317;
t474 = -t311 * rSges(6,3) + t309 * t361 + t241;
t302 = t306 ^ 2;
t303 = t307 ^ 2;
t470 = t302 + t303;
t473 = ((-t501 * t311 - t483) * qJD(2) + t489) * t311 + ((t494 * t298 + t493 * t299 - t440 * t507 - t441 * t506) * t311 - t485 * t307 + t486 * t306 + (t484 + t490) * qJD(2)) * t309;
t469 = qJD(2) * (rSges(3,1) * t309 + rSges(3,2) * t311);
t310 = cos(qJ(3));
t458 = t310 * pkin(3);
t254 = -pkin(7) * t311 + t309 * t458;
t466 = 2 * m(4);
t465 = 2 * m(5);
t464 = 2 * m(6);
t463 = t306 / 0.2e1;
t462 = -t307 / 0.2e1;
t461 = t307 / 0.2e1;
t460 = -t311 / 0.2e1;
t308 = sin(qJ(3));
t450 = Icges(4,4) * t308;
t449 = Icges(4,4) * t310;
t432 = t308 * t311;
t287 = t306 * t310 - t307 * t432;
t431 = t310 * t311;
t438 = t306 * t308;
t288 = t307 * t431 + t438;
t223 = Icges(4,5) * t288 + Icges(4,6) * t287 + Icges(4,3) * t434;
t225 = Icges(4,4) * t288 + Icges(4,2) * t287 + Icges(4,6) * t434;
t227 = Icges(4,1) * t288 + Icges(4,4) * t287 + Icges(4,5) * t434;
t285 = -t306 * t432 - t307 * t310;
t435 = t307 * t308;
t286 = t306 * t431 - t435;
t101 = t223 * t437 + t225 * t285 + t227 * t286;
t444 = t101 * t307;
t222 = Icges(4,5) * t286 + Icges(4,6) * t285 + Icges(4,3) * t437;
t224 = Icges(4,4) * t286 + Icges(4,2) * t285 + Icges(4,6) * t437;
t226 = Icges(4,1) * t286 + Icges(4,4) * t285 + Icges(4,5) * t437;
t102 = t222 * t434 + t224 * t287 + t226 * t288;
t443 = t102 * t306;
t350 = Icges(4,5) * t310 - Icges(4,6) * t308;
t401 = qJD(3) * t309;
t428 = t311 * ((-Icges(4,5) * t308 - Icges(4,6) * t310) * t401 + (Icges(4,3) * t309 + t311 * t350) * qJD(2));
t273 = -Icges(4,3) * t311 + t309 * t350;
t425 = t311 * t273;
t424 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t385 + t306 * t367 + t307 * t313;
t142 = rSges(5,1) * t237 + rSges(5,2) * t236 + rSges(5,3) * t386;
t209 = rSges(5,1) * t270 + rSges(5,2) * t269 + rSges(5,3) * t437;
t423 = t142 * t434 + t209 * t385;
t144 = rSges(5,1) * t239 + rSges(5,2) * t238 + rSges(5,3) * t385;
t402 = qJD(3) * t308;
t398 = pkin(3) * t402;
t314 = -qJD(2) * t254 - t311 * t398;
t400 = qJD(3) * t310;
t397 = pkin(3) * t400;
t175 = t306 * t397 + t307 * t314;
t422 = -t144 - t175;
t421 = t418 * t434;
t420 = t417 * t404;
t174 = t306 * t314 - t307 * t397;
t319 = pkin(7) * t309 + t311 * t458;
t230 = -pkin(3) * t435 + t306 * t319;
t416 = t174 * t434 + t230 * t385;
t362 = rSges(5,1) * t299 - rSges(5,2) * t298;
t195 = (-rSges(5,1) * t298 - rSges(5,2) * t299) * t439 + (rSges(5,3) * t309 + t311 * t362) * qJD(2);
t240 = qJD(2) * t319 - t309 * t398;
t415 = -t195 - t240;
t211 = rSges(5,1) * t272 + rSges(5,2) * t271 + rSges(5,3) * t434;
t231 = pkin(3) * t438 + t307 * t319;
t414 = -t211 - t231;
t413 = t311 * t230 + t254 * t437;
t363 = rSges(4,1) * t310 - rSges(4,2) * t308;
t235 = (-rSges(4,1) * t308 - rSges(4,2) * t310) * t401 + (rSges(4,3) * t309 + t311 * t363) * qJD(2);
t365 = pkin(2) * t311 + pkin(6) * t309;
t293 = t365 * qJD(2);
t412 = -t235 - t293;
t262 = -t311 * rSges(5,3) + t309 * t362;
t147 = t311 * t209 + t262 * t437;
t410 = -t254 - t262;
t296 = t309 * pkin(2) - t311 * pkin(6);
t409 = t470 * qJD(2) * t296;
t276 = -t311 * rSges(4,3) + t309 * t363;
t408 = -t276 - t296;
t407 = t470 * t365;
t399 = m(6) * t404;
t394 = -t175 - t424;
t393 = t311 * t142 + t195 * t437 + t262 * t386;
t392 = -t240 - t475;
t391 = -t231 - t417;
t390 = t311 * t174 + t240 * t437 + t254 * t386;
t389 = -t293 + t415;
t388 = -t254 - t474;
t387 = -t296 + t410;
t384 = t308 * t404;
t383 = t310 * t404;
t377 = t417 * t311;
t376 = t414 * t311;
t375 = t418 * t385 + t476 * t434;
t374 = -t293 + t392;
t373 = t306 * t174 + t307 * t175 - t409;
t372 = t306 * t230 + t307 * t231 + t407;
t81 = t418 * t311 + t474 * t437;
t371 = -t296 + t388;
t366 = t391 * t311;
t357 = Icges(4,1) * t310 - t450;
t353 = -Icges(4,2) * t308 + t449;
t342 = -t224 * t308 + t226 * t310;
t106 = -t311 * t222 + t309 * t342;
t341 = -t225 * t308 + t227 * t310;
t107 = -t311 * t223 + t309 * t341;
t347 = t106 * t306 + t107 * t307;
t228 = rSges(4,1) * t286 + rSges(4,2) * t285 + rSges(4,3) * t437;
t229 = rSges(4,1) * t288 + rSges(4,2) * t287 + rSges(4,3) * t434;
t340 = t228 * t307 - t229 * t306;
t274 = -Icges(4,6) * t311 + t309 * t353;
t275 = -Icges(4,5) * t311 + t309 * t357;
t336 = t274 * t308 - t275 * t310;
t333 = t476 * t311 + t474 * t386 + t475 * t437;
t250 = -qJD(3) * t286 + t306 * t384;
t251 = qJD(3) * t285 - t306 * t383;
t160 = Icges(4,5) * t251 + Icges(4,6) * t250 + Icges(4,3) * t386;
t327 = t160 * t309 + t222 * t403;
t252 = -qJD(3) * t288 + t307 * t384;
t253 = qJD(3) * t287 - t307 * t383;
t161 = Icges(4,5) * t253 + Icges(4,6) * t252 + Icges(4,3) * t385;
t326 = t161 * t309 + t223 * t403;
t320 = qJD(2) * (-Icges(3,5) * t309 - Icges(3,6) * t311);
t318 = t482 * t437 + t481 * t434 + (t483 * t309 + t484 * t311) * t404 + (-t496 * t311 + (t478 * t306 + t488) * t309) * t386 + (-t495 * t311 + (t477 * t307 + t487) * t309) * t385;
t316 = t473 * t311 + t318;
t315 = t481 * t463 + t482 * t462 + (t485 * t306 + t486 * t307) * t460 + t480 * t437 / 0.2e1 + t479 * t434 / 0.2e1 + (t502 * t306 - t503 * t307) * t404 / 0.2e1 + ((t522 * t306 - t478 * t307) * t306 + (t477 * t306 - t521 * t307) * t307) * t403 / 0.2e1;
t280 = t307 * t320;
t279 = t306 * t320;
t243 = t408 * t307;
t242 = t408 * t306;
t234 = (-Icges(4,1) * t308 - t449) * t401 + (Icges(4,5) * t309 + t311 * t357) * qJD(2);
t233 = (-Icges(4,2) * t310 - t450) * t401 + (Icges(4,6) * t309 + t311 * t353) * qJD(2);
t219 = t470 * t469;
t214 = t231 * t404;
t213 = t230 * t434;
t185 = t211 * t404;
t183 = t209 * t434;
t177 = t412 * t307;
t176 = t412 * t306;
t169 = t387 * t307;
t168 = t387 * t306;
t167 = rSges(4,1) * t253 + rSges(4,2) * t252 + rSges(4,3) * t385;
t166 = rSges(4,1) * t251 + rSges(4,2) * t250 + rSges(4,3) * t386;
t165 = Icges(4,1) * t253 + Icges(4,4) * t252 + Icges(4,5) * t385;
t164 = Icges(4,1) * t251 + Icges(4,4) * t250 + Icges(4,5) * t386;
t163 = Icges(4,4) * t253 + Icges(4,2) * t252 + Icges(4,6) * t385;
t162 = Icges(4,4) * t251 + Icges(4,2) * t250 + Icges(4,6) * t386;
t157 = -t229 * t311 - t276 * t434;
t156 = t228 * t311 + t276 * t437;
t155 = -t309 * t336 - t425;
t148 = -t211 * t311 - t262 * t434;
t124 = t273 * t434 + t274 * t287 + t275 * t288;
t123 = t273 * t437 + t274 * t285 + t275 * t286;
t122 = t340 * t309;
t121 = t371 * t307;
t120 = t371 * t306;
t119 = t389 * t307;
t118 = t389 * t306;
t117 = -t211 * t437 + t183;
t108 = t228 * t306 + t229 * t307 + t407;
t105 = t410 * t434 + t376;
t104 = t147 + t413;
t103 = t223 * t434 + t225 * t287 + t227 * t288;
t100 = t222 * t437 + t224 * t285 + t226 * t286;
t95 = t374 * t307;
t94 = t374 * t306;
t93 = t166 * t306 + t167 * t307 - t409;
t92 = -t235 * t434 - t167 * t311 + (t229 * t309 - t276 * t433) * qJD(2);
t91 = t235 * t437 + t166 * t311 + (-t228 * t309 + t276 * t436) * qJD(2);
t82 = -t434 * t474 - t377;
t80 = t414 * t437 + t183 + t213;
t79 = -t144 * t311 + t185 + (-t195 * t309 - t262 * t403) * t307;
t78 = -t209 * t404 + t393;
t77 = t209 * t306 + t211 * t307 + t372;
t76 = (t166 * t307 - t167 * t306) * t309 + t340 * t403;
t75 = -t417 * t437 + t421;
t74 = t388 * t434 + t366;
t73 = t81 + t413;
t72 = (-t144 * t309 - t211 * t403) * t306 + t423;
t71 = t142 * t306 + t144 * t307 + t373;
t68 = t391 * t437 + t213 + t421;
t67 = t306 * t418 + t307 * t417 + t372;
t62 = t185 + t214 + t422 * t311 + (t309 * t415 + t403 * t410) * t307;
t61 = (-t209 - t230) * t404 + t390 + t393;
t60 = t163 * t287 + t165 * t288 + t225 * t252 + t227 * t253 + t307 * t326;
t59 = t162 * t287 + t164 * t288 + t224 * t252 + t226 * t253 + t307 * t327;
t58 = t163 * t285 + t165 * t286 + t225 * t250 + t227 * t251 + t306 * t326;
t57 = t162 * t285 + t164 * t286 + t224 * t250 + t226 * t251 + t306 * t327;
t54 = (qJD(2) * t341 - t161) * t311 + (qJD(2) * t223 - t163 * t308 + t165 * t310 + (-t225 * t310 - t227 * t308) * qJD(3)) * t309;
t53 = (qJD(2) * t342 - t160) * t311 + (qJD(2) * t222 - t162 * t308 + t164 * t310 + (-t224 * t310 - t226 * t308) * qJD(3)) * t309;
t48 = (qJD(2) * t376 + t309 * t422) * t306 + t416 + t423;
t39 = -t424 * t311 + (-t309 * t475 - t403 * t474) * t307 + t420;
t38 = -t404 * t418 + t333;
t33 = t476 * t306 + t424 * t307 + t373;
t32 = (-qJD(2) * t377 - t309 * t424) * t306 + t375;
t31 = t214 + t394 * t311 + (t309 * t392 + t388 * t403) * t307 + t420;
t30 = (-t230 - t418) * t404 + t333 + t390;
t29 = t306 * t60 - t307 * t59;
t28 = t306 * t58 - t307 * t57;
t27 = (qJD(2) * t366 + t309 * t394) * t306 + t375 + t416;
t16 = -(t233 * t287 + t234 * t288 + t252 * t274 + t253 * t275) * t311 + (t59 * t306 + (t60 - t428) * t307) * t309 + (t124 * t309 + (t443 + (t103 - t425) * t307) * t311) * qJD(2);
t15 = -(t233 * t285 + t234 * t286 + t250 * t274 + t251 * t275) * t311 + (t58 * t307 + (t57 - t428) * t306) * t309 + (t123 * t309 + (t444 + (t100 - t425) * t306) * t311) * qJD(2);
t1 = [0; -m(3) * t219 + m(4) * t93 + m(5) * t71 + m(6) * t33; (t118 * t168 + t119 * t169 + t71 * t77) * t465 + (t120 * t94 + t121 * t95 + t33 * t67) * t464 + (t108 * t93 + t176 * t242 + t177 * t243) * t466 + 0.2e1 * m(3) * (-t219 + t469) * t470 * (rSges(3,1) * t311 - rSges(3,2) * t309) + (-t303 * t279 - t28 - t480) * t307 + (t302 * t280 + t29 + (-t306 * t279 + t307 * t280) * t307 + t479) * t306; m(4) * t76 + m(5) * t48 + m(6) * t27; (t309 * (-t106 * t307 + t107 * t306) / 0.2e1 + ((-t100 * t307 + t101 * t306) * t463 + (-t102 * t307 + t103 * t306) * t461) * t311) * qJD(2) + (t28 * t463 + t29 * t461) * t309 + m(6) * (t120 * t31 + t121 * t30 + t27 * t67 + t33 * t68 + t73 * t95 + t74 * t94) + m(5) * (t104 * t119 + t105 * t118 + t168 * t62 + t169 * t61 + t48 * t77 + t71 * t80) + m(4) * (t108 * t76 + t122 * t93 + t156 * t177 + t157 * t176 + t242 * t92 + t243 * t91) + (t54 * t306 - t53 * t307) * t460 + t15 * t462 + t16 * t463 + t315; t318 + (t27 * t68 + t30 * t73 + t31 * t74) * t464 + (t104 * t61 + t105 * t62 + t48 * t80) * t465 + (t122 * t76 + t156 * t91 + t157 * t92) * t466 + t16 * t434 + t15 * t437 + ((t100 * t306 + t444) * t386 + (t103 * t307 + t443) * t385 + t347 * t404) * t309 + (-t123 * t386 - t124 * t385 - t155 * t404 - (t155 * qJD(2) + t53 * t306 + t54 * t307) * t309 + (-t347 * qJD(2) - t428 - (-qJD(2) * t273 + t233 * t308 - t234 * t310 + t274 * t400 + t275 * t402) * t309 - t336 * t403) * t311 + t473) * t311; m(5) * t72 + m(6) * t32; m(5) * (t117 * t71 + t118 * t148 + t119 * t147 + t168 * t79 + t169 * t78 + t72 * t77) + m(6) * (t120 * t39 + t121 * t38 + t32 * t67 + t33 * t75 + t81 * t95 + t82 * t94) + t315; m(6) * (t27 * t75 + t30 * t81 + t31 * t82 + t32 * t68 + t38 * t73 + t39 * t74) + m(5) * (t104 * t78 + t105 * t79 + t117 * t48 + t147 * t61 + t148 * t62 + t72 * t80) + t316; (t32 * t75 + t38 * t81 + t39 * t82) * t464 + (t117 * t72 + t147 * t78 + t148 * t79) * t465 + t316; t399; m(6) * (-t311 * t33 + (t306 * t94 + t307 * t95) * t309 + (t309 * t67 + (t120 * t306 + t121 * t307) * t311) * qJD(2)); m(6) * (-t27 * t311 + (t30 * t307 + t306 * t31) * t309 + (t309 * t68 + (t306 * t74 + t307 * t73) * t311) * qJD(2)); m(6) * (-t311 * t32 + (t306 * t39 + t307 * t38) * t309 + (t309 * t75 + (t306 * t82 + t307 * t81) * t311) * qJD(2)); 0.2e1 * (-0.1e1 + t470) * t311 * t399;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
