% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:32
% EndTime: 2020-01-03 12:03:42
% DurationCPUTime: 5.84s
% Computational Cost: add. (43773->411), mult. (27225->524), div. (0->0), fcn. (24452->10), ass. (0->279)
t413 = qJD(1) + qJD(2);
t369 = pkin(9) + qJ(4);
t361 = sin(t369);
t362 = cos(t369);
t326 = rSges(5,1) * t361 + rSges(5,2) * t362;
t370 = qJ(1) + qJ(2);
t364 = sin(t370);
t359 = t364 ^ 2;
t365 = cos(t370);
t360 = t365 ^ 2;
t414 = t359 + t360;
t566 = t414 * t326;
t579 = -m(5) / 0.2e1;
t549 = m(6) / 0.2e1;
t363 = qJ(5) + t369;
t355 = sin(t363);
t356 = cos(t363);
t319 = rSges(6,1) * t355 + rSges(6,2) * t356;
t499 = pkin(4) * t361;
t408 = t319 + t499;
t572 = t408 * t365;
t573 = t408 * t364;
t587 = (t364 * t573 + t365 * t572) * t549;
t434 = t566 * t579 - t587;
t301 = t326 * t364;
t302 = t326 * t365;
t388 = -t301 * t364 - t302 * t365;
t435 = t388 * t579 + t587;
t57 = t435 - t434;
t589 = t413 * t57;
t371 = cos(pkin(9));
t357 = t371 * pkin(3) + pkin(2);
t498 = pkin(4) * t362;
t330 = t357 + t498;
t311 = t365 * t330;
t372 = -pkin(7) - qJ(3);
t368 = -pkin(8) + t372;
t449 = t356 * t365;
t453 = t355 * t365;
t400 = rSges(6,1) * t449 - rSges(6,2) * t453;
t220 = t311 + (rSges(6,3) - t368) * t364 + t400;
t367 = cos(qJ(1)) * pkin(1);
t214 = t367 + t220;
t209 = t214 * t364;
t211 = t220 * t364;
t454 = t355 * t364;
t416 = -rSges(6,2) * t454 - t365 * rSges(6,3);
t493 = rSges(6,1) * t356;
t219 = t365 * t368 + (t330 + t493) * t364 + t416;
t500 = sin(qJ(1)) * pkin(1);
t213 = t219 + t500;
t339 = t365 * t357;
t441 = t362 * t365;
t445 = t361 * t365;
t401 = rSges(5,1) * t441 - rSges(5,2) * t445;
t231 = t339 + (rSges(5,3) - t372) * t364 + t401;
t228 = t367 + t231;
t222 = t228 * t364;
t226 = t231 * t364;
t446 = t361 * t364;
t415 = -rSges(5,2) * t446 - t365 * rSges(5,3);
t494 = rSges(5,1) * t362;
t230 = t365 * t372 + (t357 + t494) * t364 + t415;
t227 = t230 + t500;
t570 = -rSges(4,2) * sin(pkin(9)) + pkin(2) + rSges(4,1) * t371;
t575 = -qJ(3) - rSges(4,3);
t248 = -t575 * t364 + t570 * t365;
t241 = t367 + t248;
t232 = t241 * t364;
t239 = t248 * t364;
t247 = t570 * t364 + t575 * t365;
t240 = t247 + t500;
t550 = m(5) / 0.2e1;
t580 = m(4) / 0.2e1;
t411 = (t232 + t239 + (-t240 - t247) * t365) * t580 + (t209 + t211 + (-t213 - t219) * t365) * t549 + (t222 + t226 + (-t227 - t230) * t365) * t550;
t430 = -t227 + t230;
t432 = -t213 + t219;
t412 = (t232 - t239 + (-t240 + t247) * t365) * t580 + (t365 * t432 + t209 - t211) * t549 + (t365 * t430 + t222 - t226) * t550;
t14 = t412 - t411;
t588 = t14 * qJD(1);
t286 = t319 * t365;
t221 = t573 * t286;
t481 = Icges(6,4) * t355;
t318 = Icges(6,1) * t356 - t481;
t264 = -Icges(6,5) * t365 + t318 * t364;
t450 = t356 * t364;
t236 = t264 * t450;
t261 = Icges(6,5) * t449 - Icges(6,6) * t453 + Icges(6,3) * t364;
t333 = Icges(6,4) * t453;
t265 = Icges(6,1) * t449 + Icges(6,5) * t364 - t333;
t263 = Icges(6,4) * t449 - Icges(6,2) * t453 + Icges(6,6) * t364;
t468 = t263 * t355;
t391 = t265 * t356 - t468;
t586 = -t261 * t364 - t365 * t391 - t236;
t585 = -t214 + t220;
t314 = Icges(6,5) * t356 - Icges(6,6) * t355;
t462 = t314 * t364;
t260 = -Icges(6,3) * t365 + t462;
t584 = t260 * t364 + t264 * t449;
t346 = Icges(6,4) * t356;
t316 = -Icges(6,2) * t355 + t346;
t563 = Icges(6,1) * t355 + t346;
t583 = t316 + t563;
t347 = Icges(5,4) * t362;
t323 = -Icges(5,2) * t361 + t347;
t562 = Icges(5,1) * t361 + t347;
t582 = t323 + t562;
t528 = -t364 / 0.2e1;
t576 = t364 / 0.2e1;
t527 = -t365 / 0.2e1;
t524 = m(3) * (-t367 * (rSges(3,1) * t364 + rSges(3,2) * t365) + t500 * (t365 * rSges(3,1) - rSges(3,2) * t364));
t520 = m(4) * (t240 * t248 - t247 * t241);
t327 = -rSges(5,2) * t361 + t494;
t574 = t327 * t550;
t526 = t365 / 0.2e1;
t571 = t526 + t527;
t569 = m(6) * t319;
t251 = t364 * (rSges(6,1) * t450 + t416);
t266 = rSges(6,3) * t364 + t400;
t110 = t251 + (t330 - t357) * t359 + (t311 - t339 + t266) * t365;
t285 = t319 * t364;
t404 = t364 * t285 + t365 * t286;
t320 = -rSges(6,2) * t355 + t493;
t459 = t320 * t365;
t460 = t320 * t364;
t393 = Icges(6,5) * t355 + Icges(6,6) * t356;
t279 = t364 * t393;
t280 = t393 * t365;
t425 = Icges(6,2) * t449 - t265 + t333;
t315 = Icges(6,2) * t356 + t481;
t426 = -t315 * t364 + t264;
t427 = t365 * t563 + t263;
t262 = -Icges(6,6) * t365 + t316 * t364;
t428 = -t364 * t563 - t262;
t560 = t355 * (-t364 * t425 - t365 * t426) + t356 * (t364 * t427 + t365 * t428);
t497 = (-t360 * t279 + (t365 * t280 - t560) * t364) * t527 + (t359 * t280 + (-t364 * t279 + t560) * t365) * t528;
t18 = t497 + m(6) * (-t110 * t404 + t459 * t572 + t460 * t573);
t568 = t18 * qJD(5);
t387 = -t414 * t569 / 0.2e1;
t398 = m(6) * t404;
t151 = -t398 / 0.2e1 + t387;
t567 = t413 * t151;
t91 = t213 * t220 - t219 * t214;
t101 = t227 * t231 - t230 * t228;
t565 = -t227 * t365 + t222;
t564 = -t230 * t365 + t226;
t496 = (t432 * t573 + t572 * t585) * t549 + ((-t228 + t231) * t365 + t430 * t364) * t326 * t550;
t102 = -t213 * t573 - t214 * t572;
t104 = -t219 * t573 - t220 * t572;
t117 = -t227 * t301 - t228 * t302;
t120 = -t230 * t301 - t231 * t302;
t561 = (t104 + t102) * t549 + (t120 + t117) * t550;
t559 = t583 * t355 + (t315 - t318) * t356;
t338 = Icges(5,4) * t445;
t275 = Icges(5,1) * t441 + Icges(5,5) * t364 - t338;
t421 = Icges(5,2) * t441 - t275 + t338;
t482 = Icges(5,4) * t361;
t325 = Icges(5,1) * t362 - t482;
t274 = -Icges(5,5) * t365 + t325 * t364;
t322 = Icges(5,2) * t362 + t482;
t422 = -t322 * t364 + t274;
t273 = Icges(5,4) * t441 - Icges(5,2) * t445 + Icges(5,6) * t364;
t423 = t365 * t562 + t273;
t272 = -Icges(5,6) * t365 + t323 * t364;
t424 = -t364 * t562 - t272;
t558 = t361 * (-t364 * t421 - t365 * t422) + t362 * (t364 * t423 + t365 * t424);
t557 = t582 * t361 + (t322 - t325) * t362;
t556 = t582 * t362 / 0.2e1 + (t325 / 0.2e1 - t322 / 0.2e1) * t361;
t403 = t583 * t356 / 0.2e1 + (-t315 / 0.2e1 + t318 / 0.2e1) * t355;
t132 = -t260 * t365 - t262 * t454 + t236;
t237 = t265 * t450;
t133 = t261 * t365 + t263 * t454 - t237;
t467 = t264 * t356;
t469 = t262 * t355;
t406 = ((t237 + (t260 - t468) * t364 - t584) * t364 + ((t260 + t391) * t365 + (t467 + t469) * t364 + t586) * t365) * t528 + (-t132 * t365 - t133 * t364) * t576 + ((-t133 + (t261 - t467) * t365 + t584) * t365 + (t132 + (t261 + t469) * t364 + t586) * t364) * t527;
t555 = 0.4e1 * qJD(1);
t554 = 4 * qJD(2);
t553 = 2 * qJD(4);
t111 = -t213 * t285 - t214 * t286;
t112 = -t219 * t285 - t220 * t286;
t538 = m(6) * (t112 + t111);
t537 = (t432 * t364 + t365 * t585) * t569;
t181 = t213 * t459;
t431 = -t285 * t572 + t221;
t536 = m(6) * (-t214 * t460 + t181 + t431);
t470 = t572 * t319;
t477 = t214 * t320;
t535 = m(6) * (t181 - t221 + (t470 - t477) * t364);
t188 = t219 * t459;
t534 = m(6) * (-t220 * t460 + t188 + t431);
t475 = t220 * t320;
t533 = m(6) * (t188 - t221 + (t470 - t475) * t364);
t530 = m(6) * t91;
t245 = t272 * t445;
t321 = Icges(5,5) * t362 - Icges(5,6) * t361;
t458 = t321 * t364;
t270 = -Icges(5,3) * t365 + t458;
t149 = -t270 * t364 - t274 * t441 + t245;
t271 = Icges(5,5) * t441 - Icges(5,6) * t445 + Icges(5,3) * t364;
t465 = t273 * t361;
t389 = t275 * t362 - t465;
t150 = t271 * t364 + t365 * t389;
t442 = t362 * t364;
t243 = t274 * t442;
t244 = t275 * t442;
t464 = t274 * t362;
t466 = t272 * t361;
t31 = (t149 + t244 - t245 + (t270 - t465) * t364) * t364 + (-t243 - t150 + (t270 + t389) * t365 + (t464 + t466) * t364) * t365;
t147 = -t270 * t365 - t272 * t446 + t243;
t148 = t271 * t365 + t273 * t446 - t244;
t32 = (-t148 + t245 + (t271 - t464) * t365) * t365 + (t147 - t243 + (t271 + t466) * t364) * t364;
t92 = -t147 * t365 - t148 * t364;
t93 = -t149 * t365 - t150 * t364;
t2 = (-t32 / 0.2e1 - t93 / 0.2e1) * t365 + (t92 / 0.2e1 - t31 / 0.2e1) * t364 + t406;
t525 = -qJD(3) * t57 + t2 * qJD(4);
t518 = m(4) * (-t240 * t365 + t232);
t517 = m(4) * (-t247 * t365 + t239);
t516 = m(5) * t101;
t514 = m(5) * t117;
t513 = m(5) * t120;
t512 = m(5) * t565;
t511 = m(5) * t564;
t508 = m(6) * t102;
t507 = m(6) * t104;
t506 = m(6) * t111;
t505 = m(6) * t112;
t504 = m(6) * (-t213 * t365 + t209);
t503 = m(6) * (-t219 * t365 + t211);
t488 = t151 * qJD(3) + qJD(5) * t406;
t152 = t398 / 0.2e1 + t387;
t58 = t434 + t435;
t486 = t58 * qJD(4) + t152 * qJD(5);
t485 = t57 * qJD(4) - t151 * qJD(5);
t463 = t286 * t319;
t407 = t320 + t498;
t405 = t285 * t286;
t397 = t538 / 0.2e1 + t403;
t394 = Icges(5,5) * t361 + Icges(5,6) * t362;
t385 = -t406 + (t355 * t427 + t356 * t425 + t559 * t365 - t462) * t528 + (-t314 * t365 + t355 * t428 + t356 * t426 - t559 * t364) * t527;
t384 = -t403 + t571 * (t263 * t356 + t265 * t355);
t383 = t403 + t556;
t377 = t383 + t561;
t376 = t384 - t556 + t571 * (t273 * t362 + t275 * t361);
t375 = t58 * qJD(3) + (t31 * t576 + t385 + (t361 * t423 + t362 * t421 + t557 * t365 - t458 + t92) * t528 + (-t321 * t365 + t361 * t424 + t362 * t422 - t557 * t364) * t527 + (t32 + t93) * t526) * qJD(4);
t296 = t394 * t365;
t295 = t364 * t394;
t259 = t407 * t365;
t257 = t407 * t364;
t201 = t365 * t266 + t251;
t182 = -t414 * t499 - t404;
t144 = t152 * qJD(3);
t78 = t403 + t505;
t75 = t403 + t506;
t73 = t533 / 0.2e1;
t72 = t534 / 0.2e1;
t68 = t535 / 0.2e1;
t67 = t536 / 0.2e1;
t63 = t537 / 0.2e1;
t62 = t503 + t511 + t517;
t51 = t504 + t512 + t518;
t40 = t383 + t507 + t513;
t39 = t383 + t508 + t514;
t33 = t516 + t520 + t524 + t530;
t22 = -t537 / 0.2e1 + t397;
t21 = t63 + t397;
t20 = m(6) * (t414 * t319 * t320 - t201 * t404) + t497;
t19 = t20 * qJD(5);
t17 = t63 - t538 / 0.2e1 + t384;
t16 = t411 + t412;
t13 = t377 + t496;
t12 = t377 - t496;
t9 = t376 + t496 - t561;
t8 = t72 - t533 / 0.2e1 + t406;
t7 = t73 - t534 / 0.2e1 + t406;
t6 = t67 - t535 / 0.2e1 + t406;
t5 = t68 - t536 / 0.2e1 + t406;
t4 = t72 + t73 + t385;
t3 = t67 + t68 + t385;
t1 = [qJD(2) * t33 + qJD(3) * t51 + qJD(4) * t39 + qJD(5) * t75, t33 * qJD(1) + t16 * qJD(3) + t13 * qJD(4) + t21 * qJD(5) + 0.2e1 * (t520 / 0.2e1 + t101 * t550 + t524 / 0.2e1 + t91 * t549) * qJD(2), qJD(1) * t51 + qJD(2) * t16 + t486, t39 * qJD(1) + t13 * qJD(2) + t3 * qJD(5) + ((t213 * t259 - t214 * t257) * t549 - t565 * t574) * t553 + t375, t75 * qJD(1) + t21 * qJD(2) + t144 + t3 * qJD(4) + ((t181 + (t463 - t477) * t364 - t405) * m(6) + t385) * qJD(5); -t14 * qJD(3) + t12 * qJD(4) + t22 * qJD(5) + (-t524 / 0.4e1 - t520 / 0.4e1 - t516 / 0.4e1 - t530 / 0.4e1) * t555, qJD(3) * t62 + qJD(4) * t40 + qJD(5) * t78, qJD(2) * t62 + t486 - t588, t12 * qJD(1) + t40 * qJD(2) + t4 * qJD(5) + ((t219 * t259 - t220 * t257) * t549 - t564 * t574) * t553 + t375, t22 * qJD(1) + t78 * qJD(2) + t144 + t4 * qJD(4) + ((t188 + (t463 - t475) * t364 - t405) * m(6) + t385) * qJD(5); (-t518 / 0.4e1 - t512 / 0.4e1 - t504 / 0.4e1) * t555 + t14 * qJD(2) + t485, t588 + (-t503 / 0.4e1 - t511 / 0.4e1 - t517 / 0.4e1) * t554 + t485, 0, m(6) * (t257 * t365 - t259 * t364) * qJD(4) + t589, -t567; t376 * qJD(1) + t9 * qJD(2) + t6 * qJD(5) + (-t508 / 0.4e1 - t514 / 0.4e1) * t555 + t525, t9 * qJD(1) + t8 * qJD(5) + (-t507 / 0.4e1 - t513 / 0.4e1) * t554 + t525 + t376 * qJD(2), -t589, (m(5) * (t327 * t566 + (t365 * (rSges(5,3) * t364 + t401) + t364 * (rSges(5,1) * t442 + t415)) * t388) + (-t360 * t295 + (t365 * t296 - t558) * t364) * t527 + (t359 * t296 + (-t364 * t295 + t558) * t365) * t528 + m(6) * (t110 * t182 + t257 * t573 + t259 * t572) + t497) * qJD(4) + t568 + t413 * t2, t6 * qJD(1) + t8 * qJD(2) + t18 * qJD(4) + t568; (t384 - t506) * qJD(1) + t17 * qJD(2) + t5 * qJD(4) + t488, t17 * qJD(1) + (t384 - t505) * qJD(2) + t7 * qJD(4) + t488, t567, t5 * qJD(1) + t7 * qJD(2) + ((t182 * t201 + (t257 * t364 + t259 * t365) * t319) * m(6) + t497) * qJD(4) + t19, qJD(4) * t20 + t406 * t413 + t19;];
Cq = t1;
