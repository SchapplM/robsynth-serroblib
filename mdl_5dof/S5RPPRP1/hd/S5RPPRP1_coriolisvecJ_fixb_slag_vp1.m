% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:08
% EndTime: 2020-01-03 11:25:46
% DurationCPUTime: 22.59s
% Computational Cost: add. (11133->452), mult. (15000->605), div. (0->0), fcn. (14318->8), ass. (0->223)
t510 = Icges(5,4) + Icges(6,4);
t489 = Icges(5,1) + Icges(6,1);
t467 = Icges(5,5) + Icges(6,5);
t509 = Icges(5,2) + Icges(6,2);
t465 = Icges(5,6) + Icges(6,6);
t508 = Icges(5,3) + Icges(6,3);
t246 = qJ(1) + pkin(7);
t240 = sin(t246);
t241 = cos(t246);
t252 = cos(qJ(4));
t248 = cos(pkin(8));
t250 = sin(qJ(4));
t378 = t248 * t250;
t175 = -t240 * t378 - t241 * t252;
t518 = t510 * t175;
t177 = -t240 * t252 + t241 * t378;
t517 = t510 * t177;
t377 = t248 * t252;
t382 = t240 * t250;
t178 = t241 * t377 + t382;
t516 = t510 * t178;
t380 = t241 * t250;
t176 = t240 * t377 - t380;
t515 = t510 * t176;
t247 = sin(pkin(8));
t384 = t240 * t247;
t464 = t509 * t175 + t465 * t384 + t515;
t381 = t241 * t247;
t463 = -t509 * t177 + t465 * t381 + t516;
t462 = t489 * t178 + t467 * t381 - t517;
t461 = t489 * t176 + t467 * t384 + t518;
t514 = t510 * t252;
t513 = t510 * t250;
t469 = t465 * t175 + t467 * t176 + t508 * t384;
t468 = -t465 * t177 + t467 * t178 + t508 * t381;
t500 = -t508 * t248 + (-t465 * t250 + t467 * t252) * t247;
t499 = -t465 * t248 + (-t509 * t250 + t514) * t247;
t498 = -t467 * t248 + (t489 * t252 - t513) * t247;
t512 = -t463 * t177 + t462 * t178;
t455 = -t464 * t177 + t461 * t178;
t511 = t463 * t175 + t462 * t176;
t504 = t464 * t175 + t461 * t176 + t469 * t384;
t503 = -t468 * t384 - t511;
t481 = -t469 * t381 - t455;
t480 = t468 * t381 + t512;
t433 = (-t467 * t250 - t465 * t252) * t247;
t507 = (-t509 * t252 - t513) * t247;
t506 = (-t489 * t250 - t514) * t247;
t505 = t499 * t175 + t498 * t176 + t500 * t384;
t460 = -t499 * t177 + t498 * t178 + t500 * t381;
t478 = -(-t463 * t250 + t462 * t252) * t247 + t468 * t248;
t126 = qJD(1) * t175 + qJD(4) * t178;
t268 = t177 * qJD(4);
t127 = qJD(1) * t176 + t268;
t346 = qJD(1) * t247;
t330 = t240 * t346;
t502 = t465 * t126 + t467 * t127 + t508 * t330;
t128 = -qJD(1) * t177 - qJD(4) * t176;
t267 = t175 * qJD(4);
t129 = qJD(1) * t178 + t267;
t329 = t241 * t346;
t474 = t465 * t128 + t467 * t129 + t508 * t329;
t501 = t509 * t126 + t510 * t127 + t465 * t330;
t472 = t509 * t128 + t510 * t129 + t465 * t329;
t471 = t510 * t126 + t489 * t127 + t467 * t330;
t470 = t510 * t128 + t489 * t129 + t467 * t329;
t497 = t433 * qJD(4);
t496 = t507 * qJD(4);
t495 = t506 * qJD(4);
t494 = (t481 * t240 - t480 * t241) * t247;
t493 = (t504 * t240 - t503 * t241) * t247;
t232 = -qJD(4) * t248 + qJD(1);
t491 = t505 * t232;
t490 = t502 * t241;
t388 = qJ(3) * t240;
t203 = pkin(2) * t241 + t388;
t244 = cos(qJ(1)) * pkin(1);
t239 = qJD(1) * t244;
t347 = qJD(1) * t241;
t348 = qJD(1) * t240;
t352 = pkin(2) * t347 + qJ(3) * t348;
t488 = -qJD(1) * t203 + t239 + t352;
t487 = t460 * t232;
t485 = t493 * qJD(4) + t491;
t484 = t494 * qJD(4) - t487;
t483 = -t474 * t248 + (t470 * t252 - t472 * t250 + (-t461 * t250 - t464 * t252) * qJD(4)) * t247;
t482 = t502 * t248 + (-t471 * t252 + t501 * t250 + (-t462 * t250 - t463 * t252) * qJD(4)) * t247;
t443 = (-t497 * t241 + t500 * t348) * t247 - t495 * t178 + t496 * t177 + t498 * t127 + t499 * t126;
t442 = (t497 * t240 + t500 * t347) * t247 + t495 * t176 + t496 * t175 + t498 * t129 + t499 * t128;
t479 = -t469 * t248 + (-t464 * t250 + t461 * t252) * t247;
t477 = -t497 * t248 + (t495 * t252 - t496 * t250 + (-t498 * t250 - t499 * t252) * qJD(4)) * t247;
t459 = t503 * t240 + t504 * t241;
t410 = pkin(3) * t248;
t304 = pkin(6) * t247 + t410;
t180 = t304 * t241;
t458 = -qJD(1) * t180 - t239 + t488;
t457 = t468 * t240;
t456 = t469 * t241;
t454 = -t504 + t512;
t345 = qJD(3) * t241;
t453 = t247 * t468;
t223 = pkin(4) * t382;
t237 = pkin(4) * t252 + pkin(3);
t385 = t237 * t248;
t450 = rSges(6,1) * t178 - rSges(6,2) * t177 + rSges(6,3) * t381 + t241 * t385 + t223;
t342 = qJD(5) * t247;
t212 = t240 * t342;
t249 = -qJ(5) - pkin(6);
t407 = pkin(6) + t249;
t421 = t247 * t407;
t366 = (t410 + t421) * t241 - t450;
t449 = -t366 * t232 + t212;
t427 = (-t489 * t177 - t463 - t516) * t241 + (t489 * t175 - t464 - t515) * t240;
t444 = rSges(6,1) + pkin(4);
t441 = t477 * t232;
t119 = rSges(5,1) * t178 - rSges(5,2) * t177 + rSges(5,3) * t381;
t439 = t119 * t232;
t174 = -rSges(5,3) * t248 + (rSges(5,1) * t252 - rSges(5,2) * t250) * t247;
t344 = qJD(4) * t247;
t438 = t174 * t344;
t383 = t240 * t248;
t432 = rSges(4,1) * t383 - rSges(4,2) * t384;
t431 = (t482 * t241 + t483 * t240 + (t240 * t478 + t241 * t479) * qJD(1)) * t247;
t430 = (t459 * qJD(1) + (t128 * t463 + t129 * t462 - t175 * t501 - t176 * t471 + t347 * t453) * t241 + ((t240 * t474 + t347 * t469 - t490) * t247 + t470 * t176 + t472 * t175 + t461 * t129 + t464 * t128) * t240) * t247;
t429 = (((t348 * t468 + t490) * t247 + t471 * t178 - t501 * t177 + t462 * t127 + t463 * t126 + t481 * qJD(1)) * t241 + ((-t241 * t474 + t348 * t469) * t247 - t470 * t178 + t472 * t177 + t461 * t127 + t464 * t126 + t480 * qJD(1)) * t240) * t247;
t428 = (-t178 * t509 + t462 - t517) * t241 + (-t176 * t509 + t461 + t518) * t240;
t426 = (t177 * t467 + t178 * t465) * t241 + (-t175 * t467 + t176 * t465) * t240;
t424 = t498 + t507;
t423 = t499 - t506;
t328 = t248 * t347;
t422 = t129 * rSges(6,1) + t128 * rSges(6,2) + rSges(6,3) * t329 + qJD(1) * t223 + t237 * t328 + t212;
t301 = -t240 * rSges(3,1) - t241 * rSges(3,2);
t243 = sin(qJ(1)) * pkin(1);
t349 = t240 * pkin(2) + t243;
t401 = rSges(4,3) * t241;
t419 = t401 - t432;
t213 = t241 * t342;
t418 = rSges(6,1) * t127 + rSges(6,2) * t126 - t213;
t379 = t247 * t249;
t417 = t176 * rSges(6,1) + t175 * rSges(6,2) + rSges(6,3) * t384 + t237 * t383 - t240 * t379;
t254 = qJD(1) ^ 2;
t409 = pkin(4) * t250;
t408 = -pkin(3) + t237;
t339 = pkin(4) * t380;
t406 = rSges(6,3) * t330 + pkin(4) * t268 + (-t339 + (t248 * t408 - t421) * t240) * qJD(1) + t418;
t353 = pkin(3) * t328 + pkin(6) * t329;
t405 = pkin(4) * t267 - t249 * t329 - t353 + t422;
t402 = rSges(4,3) * t240;
t400 = pkin(4) * qJD(4);
t179 = pkin(3) * t383 + pkin(6) * t384;
t367 = -t339 - t179 + t417;
t365 = -rSges(6,2) * t176 + t175 * t444;
t364 = rSges(6,2) * t178 + t177 * t444;
t238 = t254 * t244;
t363 = qJD(1) * (-t345 + t352) + t238;
t362 = (t407 - rSges(6,3)) * t248 + (rSges(6,1) * t252 - rSges(6,2) * t250 + t408) * t247;
t233 = qJD(3) * t240;
t351 = qJ(3) * t347 + t233;
t150 = pkin(2) * t348 - t351;
t361 = -t304 * t348 - t150;
t275 = (-rSges(6,1) * t250 - rSges(6,2) * t252) * t247;
t336 = t250 * t400;
t341 = qJD(5) * t248;
t356 = -qJD(4) * t275 + t247 * t336 + t341;
t354 = rSges(4,1) * t328 + rSges(4,3) * t348;
t340 = qJD(1) * qJD(3);
t338 = t254 * t243;
t335 = t252 * t400;
t72 = t129 * rSges(5,1) + t128 * rSges(5,2) + rSges(5,3) * t329;
t334 = qJD(1) * t353 + t363;
t117 = t176 * rSges(5,1) + t175 * rSges(5,2) + rSges(5,3) * t384;
t327 = t240 * t344;
t326 = t241 * t344;
t322 = -t344 / 0.2e1;
t321 = t344 / 0.2e1;
t320 = -qJ(3) * t241 + t349;
t269 = (t179 + t320) * qJD(1) - t233;
t302 = t362 * t344;
t27 = t232 * t367 - t240 * t302 - t213 + t269;
t319 = t27 * t362;
t318 = t244 + t388;
t317 = t239 - t345;
t315 = t248 * t336;
t314 = t174 * t327;
t311 = t240 * t322;
t310 = t240 * t321;
t309 = t241 * t322;
t308 = t241 * t321;
t307 = t240 * t340 - t338;
t306 = qJD(1) * t321;
t204 = rSges(3,1) * t241 - rSges(3,2) * t240;
t300 = rSges(4,1) * t248 - rSges(4,2) * t247;
t299 = rSges(5,1) * t127 + rSges(5,2) * t126;
t276 = (-rSges(5,1) * t250 - rSges(5,2) * t252) * t247;
t188 = qJD(4) * t276;
t70 = rSges(5,3) * t330 + t299;
t29 = -t188 * t326 - t232 * t70 + (t314 + t361) * qJD(1) + t307;
t274 = (-qJD(3) - t438) * t241;
t30 = qJD(1) * t274 - t188 * t327 + t232 * t72 + t334;
t296 = -t30 * t240 - t29 * t241;
t44 = t117 * t232 + t269 - t314;
t285 = t239 + (t180 + t203) * qJD(1);
t45 = t274 + t285 + t439;
t293 = -t240 * t44 - t241 * t45;
t292 = t240 * t70 + t241 * t72;
t289 = t117 * t241 - t119 * t240;
t284 = t240 * t306;
t283 = t241 * t306;
t282 = pkin(2) + t300;
t146 = t241 * t300 + t402;
t145 = rSges(5,1) * t177 + rSges(5,2) * t178;
t143 = rSges(5,1) * t175 - rSges(5,2) * t176;
t95 = (t146 + t203) * qJD(1) + t317;
t88 = ((-rSges(4,2) * t346 - qJD(3)) * t241 + t354) * qJD(1) + t363;
t87 = (qJD(1) * t419 - t150) * qJD(1) + t307;
t48 = t289 * t344 + qJD(2);
t28 = (-qJD(3) - t302) * t241 + t285 + t449;
t25 = -t341 + qJD(2) + (t240 * t366 + t241 * t367) * t344;
t20 = ((-t117 * t240 - t119 * t241) * qJD(1) + t292) * t344;
t19 = -t241 * t340 + t405 * t232 + (qJD(5) * t348 + (t240 * t356 - t347 * t362) * qJD(4)) * t247 + t334;
t18 = -t406 * t232 + t356 * t326 + ((qJD(4) * t240 * t362 + qJD(5) * t241) * t247 + t361) * qJD(1) + t307;
t1 = (t405 * t241 + t406 * t240 + (-t240 * t367 + t241 * t366) * qJD(1)) * t344;
t2 = [m(3) * ((t254 * t301 - t338) * (t204 + t244) + (-t238 + (-0.2e1 * rSges(3,1) * t347 + 0.2e1 * rSges(3,2) * t348 + qJD(1) * t204) * qJD(1)) * (t301 - t243)) + t441 + ((((t456 - t457) * t247 + t455 + t481 - t503) * t241 + (-t454 + t480) * t240) * t344 + t491) * t308 + (t18 * (t318 + t450) + t19 * (t349 + t417) + (t240 * t335 + t351 - t418 + (t339 - t243 + (-t385 - pkin(2) + (-rSges(6,3) + t249) * t247) * t240) * qJD(1)) * t28 + (-t240 * t315 + t28 + t422 - t449 + t458) * t27 + (t18 * (pkin(2) - t379) - t28 * t315 + t19 * (-qJ(3) - t409) + t319 * t344 + (-t379 * qJD(1) - t335) * t27) * t241) * m(6) + (t29 * (t119 + t318) + t30 * (t117 + t349 + t179) + (-t299 + t351 + (-t243 + (-t410 - pkin(2) + (-rSges(5,3) - pkin(6)) * t247) * t240) * qJD(1)) * t45 + (t72 + t353 + t45 - t439 + t458) * t44 + (t29 * (pkin(2) + t304) - t30 * qJ(3) + t438 * t44) * t241) * m(5) + (t87 * (t318 + t402) + t88 * (t349 + t432) + (t87 * t282 + t88 * (-rSges(4,3) - qJ(3))) * t241 + (t351 + (-t282 * t240 - t243 + t401) * qJD(1)) * t95 + (-t345 - t317 + t354 + t95 + (-rSges(4,2) * t381 - t146) * qJD(1) + t488) * (-t233 + (t320 - t419) * qJD(1))) * m(4) + (t442 + t483) * t310 + (-t460 + t478) * t284 + (t505 + t479) * t283 + (t443 - t482 + t485) * t309 + ((((t456 + t457) * t247 + t455 + t511) * t240 + t459 + (t241 * t453 + t454) * t241) * t344 + t484 + t487) * t311; m(5) * t20 + m(6) * t1; m(4) * (-t240 * t88 - t241 * t87) + m(5) * t296 + m(6) * (-t18 * t241 - t19 * t240); -(((-t250 * t424 - t252 * t423) * t232 + ((-t250 * t428 + t252 * t427) * t247 + t426 * t248) * qJD(4)) * t247 - t433 * t232 * t248) * t232 / 0.2e1 + (-t477 * t248 + t431) * t232 / 0.2e1 - (t431 * qJD(4) + t441) * t248 / 0.2e1 + (qJD(4) * t430 + t232 * t442) * t384 / 0.2e1 - (t429 * qJD(4) + t443 * t232) * t381 / 0.2e1 + ((t175 * t428 + t176 * t427 - t384 * t426) * t344 + (t175 * t424 - t176 * t423 + t384 * t433) * t232) * t311 + (-t248 * t442 + t430) * t310 + (-t248 * t443 + t429) * t309 + ((t428 * t177 - t178 * t427 + t426 * t381) * t344 + (t177 * t424 + t178 * t423 - t381 * t433) * t232) * t308 + (t248 * t460 + t494) * t284 + (-t248 * t505 + t493) * t283 + (-(t27 * t365 - t28 * t364) * t232 - (t25 * (t240 * t364 + t241 * t365) + (t247 * t409 - t275) * (t240 * t27 + t241 * t28)) * t344 + (t18 * t366 - t19 * t367 - t27 * t405 + t28 * t406) * t248 + ((t1 * t367 + t25 * t405 - t18 * t362 + t28 * t356 + (t25 * t366 - t319) * qJD(1)) * t241 + (t1 * t366 + t25 * t406 - t19 * t362 + t27 * t356 + (-t25 * t367 + t28 * t362) * qJD(1)) * t240) * t247) * m(6) + (-(t143 * t44 - t145 * t45) * t232 - (t48 * (t143 * t241 + t145 * t240) + t293 * t276) * t344 + (-t117 * t30 - t29 * t119 - t44 * t72 + t45 * t70) * t248 + (t20 * t289 + t48 * (-t117 * t348 - t119 * t347 + t292) + t293 * t188 + ((t240 * t45 - t241 * t44) * qJD(1) + t296) * t174) * t247) * m(5) + (t484 * t240 + t241 * t485) * t346 / 0.2e1; (-t1 * t248 + (t18 * t240 - t19 * t241) * t247) * m(6);];
tauc = t2(:);
