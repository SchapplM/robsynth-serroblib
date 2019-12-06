% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:22
% EndTime: 2019-12-05 15:40:43
% DurationCPUTime: 13.87s
% Computational Cost: add. (12736->620), mult. (33991->902), div. (0->0), fcn. (36567->6), ass. (0->363)
t574 = rSges(6,1) + pkin(4);
t569 = rSges(6,3) + qJ(5);
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t335 = cos(qJ(4));
t333 = sin(qJ(4));
t483 = Icges(6,5) * t333;
t382 = -Icges(6,3) * t335 + t483;
t344 = -Icges(6,6) * t334 + t382 * t336;
t489 = Icges(5,4) * t333;
t386 = Icges(5,2) * t335 + t489;
t345 = -Icges(5,6) * t334 + t386 * t336;
t586 = t344 - t345;
t488 = Icges(5,4) * t335;
t391 = Icges(5,1) * t333 + t488;
t347 = -Icges(5,5) * t334 + t391 * t336;
t482 = Icges(6,5) * t335;
t390 = Icges(6,1) * t333 - t482;
t348 = -Icges(6,4) * t334 + t390 * t336;
t585 = t347 + t348;
t535 = t334 ^ 2;
t330 = t336 ^ 2;
t456 = t334 * t336;
t331 = sin(pkin(7));
t332 = cos(pkin(7));
t457 = t334 * t335;
t283 = t331 * t333 - t332 * t457;
t468 = t333 * t334;
t472 = t331 * t335;
t284 = t332 * t468 + t472;
t198 = -Icges(5,5) * t283 - Icges(5,6) * t284;
t200 = -Icges(6,4) * t283 + Icges(6,6) * t284;
t601 = t198 + t200;
t285 = t331 * t457 + t332 * t333;
t470 = t332 * t335;
t286 = -t331 * t468 + t470;
t199 = Icges(5,5) * t285 + Icges(5,6) * t286;
t201 = Icges(6,4) * t285 - Icges(6,6) * t286;
t600 = t199 + t201;
t383 = Icges(5,5) * t333 + Icges(5,6) * t335;
t343 = -Icges(5,3) * t334 + t383 * t336;
t385 = Icges(6,4) * t333 - Icges(6,6) * t335;
t346 = -Icges(6,2) * t334 + t385 * t336;
t599 = t343 + t346;
t280 = Icges(5,4) * t285;
t471 = t331 * t336;
t183 = -Icges(5,1) * t286 + Icges(5,5) * t471 + t280;
t443 = Icges(5,2) * t286 + t183 + t280;
t484 = Icges(6,5) * t285;
t181 = -Icges(6,1) * t286 + Icges(6,4) * t471 - t484;
t445 = Icges(6,3) * t286 + t181 - t484;
t598 = t443 + t445;
t279 = Icges(5,4) * t283;
t469 = t332 * t336;
t182 = Icges(5,1) * t284 + Icges(5,5) * t469 - t279;
t444 = -Icges(5,2) * t284 + t182 - t279;
t485 = Icges(6,5) * t283;
t180 = Icges(6,1) * t284 + Icges(6,4) * t469 + t485;
t446 = -Icges(6,3) * t284 + t180 + t485;
t597 = t444 + t446;
t490 = Icges(5,4) * t286;
t179 = Icges(5,2) * t285 + Icges(5,6) * t471 - t490;
t447 = -Icges(5,1) * t285 + t179 - t490;
t278 = Icges(6,5) * t286;
t173 = Icges(6,6) * t471 - Icges(6,3) * t285 - t278;
t449 = Icges(6,1) * t285 + t173 - t278;
t596 = t447 - t449;
t491 = Icges(5,4) * t284;
t178 = -Icges(5,2) * t283 + Icges(5,6) * t469 + t491;
t448 = Icges(5,1) * t283 + t178 + t491;
t277 = Icges(6,5) * t284;
t172 = Icges(6,6) * t469 + Icges(6,3) * t283 + t277;
t450 = -Icges(6,1) * t283 + t172 + t277;
t595 = t448 - t450;
t568 = Icges(4,2) - Icges(4,3);
t593 = t172 - t178;
t592 = t173 - t179;
t174 = Icges(5,5) * t284 - Icges(5,6) * t283 + Icges(5,3) * t469;
t176 = Icges(6,4) * t284 + Icges(6,2) * t469 + Icges(6,6) * t283;
t591 = t174 + t176;
t175 = -Icges(5,5) * t286 + Icges(5,6) * t285 + Icges(5,3) * t471;
t177 = -Icges(6,4) * t286 + Icges(6,2) * t471 - Icges(6,6) * t285;
t590 = t175 + t177;
t589 = t180 + t182;
t588 = t181 + t183;
t584 = t585 * t333 - t586 * t335;
t583 = t574 * t333 - t569 * t335;
t582 = -t597 * t283 - t595 * t284 + t601 * t469;
t581 = -t598 * t283 - t596 * t284 + t600 * t469;
t580 = t597 * t285 + t595 * t286 + t601 * t471;
t579 = t598 * t285 + t596 * t286 + t600 * t471;
t578 = -t585 + (t482 - t488 + (Icges(5,2) + Icges(6,3)) * t333) * t336;
t577 = t586 + (t483 - t489 + (Icges(5,1) + Icges(6,1)) * t335) * t336;
t300 = (-Icges(5,5) * t335 + Icges(5,6) * t333) * t336;
t301 = (-Icges(6,4) * t335 - Icges(6,6) * t333) * t336;
t576 = (t300 + t301) * t334;
t550 = t599 * t334;
t575 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t456 + (0.2e1 * t330 - 0.2e1 * t535) * Icges(3,4);
t529 = m(5) / 0.2e1;
t527 = m(6) / 0.2e1;
t221 = t344 * t331;
t227 = t345 * t331;
t229 = t348 * t331;
t231 = t347 * t331;
t372 = -t179 * t335 - t183 * t333;
t359 = t343 * t331 - t372;
t374 = t173 * t335 - t181 * t333;
t361 = t346 * t331 - t374;
t573 = ((t221 - t227) * t335 + (-t229 - t231) * t333 + t590) * t336 + (t361 + t359) * t334;
t222 = t344 * t332;
t228 = t345 * t332;
t230 = t348 * t332;
t232 = t347 * t332;
t373 = -t178 * t335 - t182 * t333;
t360 = t343 * t332 - t373;
t375 = t172 * t335 - t180 * t333;
t362 = t346 * t332 - t375;
t572 = ((t222 - t228) * t335 + (-t230 - t232) * t333 + t591) * t336 + (t360 + t362) * t334;
t571 = t593 * t283 + t589 * t284 + t591 * t469;
t570 = t592 * t283 + t588 * t284 + t590 * t469;
t384 = -Icges(3,5) * t334 - Icges(3,6) * t336;
t291 = t384 * t331;
t292 = t384 * t332;
t567 = -t593 * t285 - t589 * t286 + t591 * t471;
t566 = -t592 * t285 - t588 * t286 + t590 * t471;
t565 = -t586 * t283 - t585 * t284 - t469 * t599;
t564 = t586 * t285 + t585 * t286 - t471 * t599;
t563 = t584 * t336 - t550;
t562 = (-Icges(5,6) + Icges(6,6)) * t336 + (t382 - t386) * t334;
t561 = (Icges(6,4) + Icges(5,5)) * t336 + (t390 + t391) * t334;
t387 = Icges(4,4) * t334 + Icges(4,5) * t336;
t560 = t387 * t331 + t291;
t559 = t387 * t332 + t292;
t328 = t331 ^ 2;
t329 = t332 ^ 2;
t422 = t328 + t329;
t558 = -t334 * rSges(6,2) + t583 * t336;
t480 = Icges(4,6) * t334;
t553 = -0.2e1 * Icges(4,6) * t336 - t568 * t334;
t557 = t291 - t575 * t332 + (Icges(4,5) * t331 + t553 * t332) * t336 + (Icges(4,4) * t331 + 0.2e1 * t332 * t480 - t568 * t469) * t334;
t556 = t292 + t575 * t331 + (Icges(4,5) * t332 - t553 * t331) * t336 + (Icges(4,4) * t332 - 0.2e1 * t331 * t480 + t568 * t471) * t334;
t555 = ((t383 + t385) * t334 + (Icges(6,2) + Icges(5,3)) * t336 - t584) * t334;
t465 = t334 * t176;
t111 = t375 * t336 + t465;
t464 = t334 * t177;
t112 = t374 * t336 + t464;
t467 = t334 * t174;
t113 = t373 * t336 + t467;
t466 = t334 * t175;
t114 = t372 * t336 + t466;
t554 = (t111 + t113) * t332 + (t112 + t114) * t331;
t189 = t284 * rSges(5,1) - t283 * rSges(5,2) + rSges(5,3) * t469;
t191 = -t286 * rSges(5,1) + t285 * rSges(5,2) + rSges(5,3) * t471;
t127 = (-t189 * t331 + t191 * t332) * t336;
t397 = rSges(5,1) * t333 + rSges(5,2) * t335;
t349 = -t334 * rSges(5,3) + t397 * t336;
t462 = t334 * t191;
t145 = -t349 * t471 - t462;
t463 = t334 * t189;
t146 = t349 * t469 + t463;
t453 = t145 * t469 + t146 * t471;
t441 = rSges(6,2) * t471 - t569 * t285 - t574 * t286;
t404 = t441 * t334;
t117 = -t471 * t558 - t404;
t442 = rSges(6,2) * t469 + t569 * t283 + t574 * t284;
t405 = t442 * t334;
t118 = t469 * t558 + t405;
t454 = t117 * t469 + t118 * t471;
t540 = t422 * t334;
t92 = (-t442 * t331 + t441 * t332) * t336;
t500 = (t540 * t92 + t454) * t527 + (t127 * t540 + t453) * t529;
t234 = t349 * t331;
t236 = t349 * t332;
t105 = (t234 * t336 - t462) * t332 + (-t236 * t336 + t463) * t331;
t267 = rSges(5,3) * t336 + t397 * t334;
t369 = t267 * t336 + t334 * t349;
t120 = -t191 * t336 - t334 * t234 + t369 * t331;
t121 = t189 * t336 + t334 * t236 - t369 * t332;
t436 = t558 * t332;
t437 = t558 * t331;
t56 = (t437 * t336 - t404) * t332 + (-t436 * t336 + t405) * t331;
t430 = rSges(6,2) * t336 + t583 * t334;
t537 = -t334 * t558 - t430 * t336;
t83 = -t537 * t331 - t437 * t334 - t441 * t336;
t84 = t537 * t332 + t436 * t334 + t442 * t336;
t501 = (-t336 * t56 + (t331 * t84 + t332 * t83 + t92) * t334 + t454) * t527 + (-t105 * t336 + (t120 * t332 + t121 * t331 + t127) * t334 + t453) * t529;
t549 = t500 - t501;
t548 = -t331 / 0.2e1;
t510 = t331 / 0.2e1;
t509 = -t332 / 0.2e1;
t547 = t332 / 0.2e1;
t546 = -t334 / 0.2e1;
t507 = t334 / 0.2e1;
t545 = -t336 / 0.2e1;
t544 = t336 / 0.2e1;
t543 = ((t576 + t582) * t332 + t581 * t331) * t336 + (-t578 * t283 - t577 * t284) * t334;
t542 = (t580 * t332 + (t576 + t579) * t331) * t336 + (t578 * t285 + t577 * t286) * t334;
t526 = m(6) / 0.4e1;
t528 = m(5) / 0.4e1;
t402 = m(4) / 0.4e1 + t528 + t526;
t423 = t422 * t456;
t541 = t402 * (t423 - t456);
t192 = t283 * t332 - t285 * t331;
t455 = t335 * t336;
t161 = (-t192 - t457) * t455;
t409 = -t457 / 0.2e1;
t165 = (t409 - t192 / 0.2e1) * m(6);
t497 = m(6) * qJD(5);
t539 = t165 * qJD(1) + t161 * t497;
t439 = t574 * t285 - t569 * t286;
t440 = -t574 * t283 + t569 * t284;
t110 = t439 * t331 + t440 * t332;
t210 = -rSges(5,1) * t283 - rSges(5,2) * t284;
t214 = rSges(5,1) * t285 + rSges(5,2) * t286;
t137 = t210 * t332 + t214 * t331;
t426 = (-t569 * t333 - t574 * t335) * t336;
t194 = t426 * t331;
t195 = t426 * t332;
t309 = (-rSges(5,1) * t335 + rSges(5,2) * t333) * t336;
t494 = (-t137 * t336 - t309 * t540) * t529 + (-t110 * t336 + (-t194 * t331 - t195 * t332) * t334) * t527;
t534 = 2 * qJD(2);
t533 = 4 * qJD(2);
t532 = 2 * qJD(4);
t531 = 4 * qJD(4);
t530 = m(4) / 0.2e1;
t524 = m(5) * (t105 * t127 + t120 * t145 + t121 * t146);
t314 = pkin(2) * t336 + t334 * qJ(3);
t427 = t422 * t314;
t400 = t427 + (t331 * t471 + t332 * t469) * pkin(6);
t108 = t189 * t332 + t191 * t331 + t400;
t311 = t334 * pkin(2) - qJ(3) * t336;
t407 = -pkin(6) * t334 - t311;
t398 = t349 + t407;
t184 = t398 * t331;
t186 = t398 * t332;
t523 = m(5) * (t108 * t137 + (-t184 * t331 - t186 * t332) * t309);
t451 = t184 * t471 + t186 * t469;
t522 = m(5) * (t108 * t540 + t451);
t377 = -t117 * t332 - t118 * t331;
t520 = m(6) * (t283 * t83 - t285 * t84 + (-t92 * t334 + (t377 + t56) * t336) * t335);
t518 = m(6) * (t117 * t83 + t118 * t84 + t56 * t92);
t366 = t407 + t558;
t155 = t366 * t331;
t157 = t366 * t332;
t170 = (-t283 * t331 - t285 * t332) * t336;
t219 = t285 * t334 + t330 * t472;
t220 = t283 * t334 - t330 * t470;
t82 = t441 * t331 + t442 * t332 + t400;
t517 = m(6) * (t220 * t155 + t219 * t157 + t170 * t82 + t92 * t192 + t377 * t455);
t516 = m(6) * (-t286 * t155 + t284 * t157 + t285 * t194 - t283 * t195 + (t110 * t335 - t333 * t82) * t336);
t515 = m(6) * (t110 * t82 - t155 * t194 - t157 * t195);
t452 = t155 * t471 + t157 * t469;
t513 = m(6) * (t540 * t82 + t452);
t315 = -rSges(4,2) * t336 + t334 * rSges(4,3);
t133 = t422 * t315 + t427;
t395 = t334 * rSges(4,2) + rSges(4,3) * t336;
t425 = -t311 + t395;
t237 = t425 * t331;
t239 = t425 * t332;
t438 = t237 * t471 + t239 * t469;
t506 = m(4) * (t133 * t540 + t438);
t435 = t283 * t469 - t285 * t471;
t504 = m(6) * ((0.2e1 - t422) * t334 * t455 + t435);
t503 = m(6) * (-t335 * t540 - t192) * t336;
t502 = m(6) * (t455 * t540 + t435);
t499 = m(6) * qJD(2);
t498 = m(6) * qJD(4);
t23 = 0.2e1 * (t56 / 0.4e1 - t110 / 0.4e1) * m(6) + 0.2e1 * (t105 / 0.4e1 - t137 / 0.4e1) * m(5);
t475 = t23 * qJD(1);
t473 = t330 * t333;
t428 = t422 * t311;
t424 = -t314 - t315;
t421 = qJD(4) * t336;
t363 = 0.2e1 * t402 * t540;
t401 = t527 + t529 + t530;
t152 = -t401 * t334 + t363;
t420 = t152 * qJD(1);
t340 = t362 * t336 - t465;
t57 = -t285 * t222 - t286 * t230 + t340 * t331;
t339 = t361 * t336 - t464;
t58 = -t285 * t221 - t286 * t229 + t339 * t331;
t342 = t360 * t336 - t467;
t59 = t285 * t228 - t286 * t232 + t342 * t331;
t341 = t359 * t336 - t466;
t60 = t285 * t227 - t286 * t231 + t341 * t331;
t417 = ((t57 + t59) * t336 - t567 * t334) * t547 + t564 * t544 + ((t60 + t58 + t555) * t336 + (t550 - t566) * t334) * t510 + (-t562 * t285 - t561 * t286) * t507;
t61 = t283 * t222 + t284 * t230 + t340 * t332;
t62 = t283 * t221 + t284 * t229 + t339 * t332;
t63 = -t283 * t228 + t284 * t232 + t342 * t332;
t64 = -t283 * t227 + t284 * t231 + t341 * t332;
t416 = ((t63 + t61 + t555) * t336 + (t550 - t571) * t334) * t547 + t565 * t544 + ((t62 + t64) * t336 - t570 * t334) * t510 + (t562 * t283 + t561 * t284) * t507;
t415 = ((-t561 * t333 + t562 * t335 - t599) * t336 - t554 + t555) * t546 + (t573 * t331 + t572 * t332 + t563) * t545;
t414 = t581 * t509 + t582 * t510;
t413 = t579 * t547 + t580 * t548;
t412 = t565 * t546 + (t570 * t331 + t571 * t332) * t545;
t411 = t564 * t546 + (t566 * t331 + t567 * t332) * t545;
t410 = t554 * t545 + t563 * t546;
t408 = t333 * t545;
t406 = -pkin(6) * t336 - t314;
t399 = -t267 + t406;
t313 = t334 * rSges(3,1) + rSges(3,2) * t336;
t376 = -t155 * t331 - t157 * t332;
t162 = (t408 - t170 / 0.2e1) * m(6);
t337 = (-t170 * t336 + (t219 * t332 + t220 * t331) * t334) * t527;
t338 = m(6) * (t473 + (t284 * t332 - t286 * t331) * t334);
t94 = t337 - t338 / 0.2e1;
t368 = -t162 * qJD(1) + t94 * qJD(3);
t367 = t406 - t430;
t365 = -t413 - t417;
t364 = t414 - t416;
t357 = -pkin(6) * t540 - t428;
t240 = t424 * t332;
t238 = t424 * t331;
t193 = t422 * t313;
t187 = t399 * t332;
t185 = t399 * t331;
t166 = m(6) * t409 + t192 * t527;
t163 = m(6) * t408 + t170 * t527;
t160 = t334 * t210 - t309 * t469;
t159 = -t334 * t214 + t309 * t471;
t158 = t367 * t332;
t156 = t367 * t331;
t151 = t363 + (m(4) + m(5) + m(6)) * t507;
t148 = t502 / 0.2e1;
t147 = t503 / 0.2e1;
t144 = t283 * t284 + t285 * t286 - t335 * t473;
t140 = t422 * t395 - t428;
t136 = t504 / 0.2e1;
t134 = (-t210 * t331 + t214 * t332) * t336;
t132 = 0.4e1 * t541;
t126 = t234 * t331 + t236 * t332 + t357;
t125 = -t336 * t195 + t440 * t334;
t124 = -t439 * t334 + t426 * t471;
t107 = t437 * t331 + t436 * t332 + t357;
t104 = (-t440 * t331 + t439 * t332) * t336;
t93 = t337 + t338 / 0.2e1;
t91 = t148 + t136 - t503 / 0.2e1;
t90 = t147 + t148 - t504 / 0.2e1;
t89 = t147 + t136 - t502 / 0.2e1;
t88 = t334 * t199 + (t447 * t333 - t443 * t335) * t336;
t87 = t334 * t198 + (t448 * t333 - t444 * t335) * t336;
t86 = t334 * t201 + (-t449 * t333 - t445 * t335) * t336;
t85 = t334 * t200 + (-t450 * t333 - t446 * t335) * t336;
t50 = t82 * t192 + t376 * t455;
t39 = t117 * t219 + t118 * t220 + t170 * t92;
t38 = t331 * t63 - t332 * t64;
t37 = t331 * t61 - t332 * t62;
t36 = t331 * t59 - t332 * t60;
t35 = t331 * t57 - t332 * t58;
t31 = t516 / 0.2e1;
t24 = (t105 + t137) * t529 + (t110 + t56) * t527;
t22 = t506 + t513 + t522;
t20 = t517 / 0.2e1;
t9 = t520 / 0.2e1;
t8 = t31 + t9 - t517 / 0.2e1;
t7 = t20 + t31 - t520 / 0.2e1;
t6 = t20 + t9 - t516 / 0.2e1;
t5 = t500 + t501 - t494;
t4 = t494 - t549;
t3 = t494 + t549;
t2 = t414 * t331 + t413 * t332 + t515 + t523;
t1 = t524 + t518 + (t417 * t331 + t416 * t332 - t410) * t336 + (t411 * t331 + t412 * t332 - t415) * t334;
t10 = [0, t151 * qJD(3) + t24 * qJD(4) + t166 * qJD(5) + (-m(3) * t193 / 0.2e1 + t140 * t530 + t126 * t529 + t107 * t527) * t534, t151 * qJD(2), t24 * qJD(2) + (t104 * t527 + t134 * t529) * t532 + t163 * qJD(5), qJD(2) * t166 + qJD(4) * t163; qJD(3) * t152 - qJD(4) * t23 - qJD(5) * t165, t22 * qJD(3) + t2 * qJD(4) + t50 * t497 + (m(5) * (t108 * t126 + t184 * t185 + t186 * t187) + m(4) * (t133 * t140 + t237 * t238 + t239 * t240) + m(6) * (t107 * t82 + t155 * t156 + t157 * t158) + m(3) * (-t193 + t313) * t422 * (rSges(3,1) * t336 - t334 * rSges(3,2)) + (t37 + t38 + t559 * t328 + (t556 * t332 + (t557 - t560) * t331) * t332) * t510 + (t35 + t36 + t560 * t329 + (t557 * t331 + (t556 - t559) * t332) * t331) * t509) * qJD(2), t22 * qJD(2) + t3 * qJD(4) + t90 * qJD(5) + t420 + (-0.4e1 * t541 + 0.2e1 * t401 * (-t336 * t540 + t423)) * qJD(3), -t475 + t2 * qJD(2) + t3 * qJD(3) + t7 * qJD(5) + (-t524 / 0.4e1 - t518 / 0.4e1) * t531 + ((t108 * t134 + t127 * t137 + t159 * t186 + t160 * t184 + (-t145 * t332 - t146 * t331) * t309) * t529 + (t104 * t82 + t110 * t92 - t117 * t195 - t118 * t194 + t124 * t157 + t125 * t155) * t527) * t532 + (t331 * t365 + t332 * t364 + t410) * t421 + (((-t88 / 0.2e1 - t86 / 0.2e1 - t412) * t332 + (t87 / 0.2e1 + t85 / 0.2e1 - t411) * t331 + t415) * t334 + t543 * t510 + t542 * t509) * qJD(4), t90 * qJD(3) + t7 * qJD(4) + t499 * t50 - t539; -t152 * qJD(2), -t420 + t132 * qJD(3) + t4 * qJD(4) + t89 * qJD(5) + (-t513 / 0.4e1 - t522 / 0.4e1 - t506 / 0.4e1) * t533 + ((-t107 * t336 + t452) * t527 + (-t126 * t336 + t451) * t529 + (-t140 * t336 + t438) * t530 + ((t156 * t331 + t158 * t332 + t82) * t527 + (t185 * t331 + t187 * t332 + t108) * t529 + (t238 * t331 + t240 * t332 + t133) * t530) * t334) * t534, t132 * qJD(2), t4 * qJD(2) + ((-t134 * t336 + (t159 * t332 + t160 * t331) * t334) * t529 + (-t104 * t336 + (t124 * t332 + t125 * t331) * t334) * t527) * t532 + t93 * qJD(5), qJD(2) * t89 + qJD(4) * t93; qJD(2) * t23 - qJD(5) * t162, t475 + t5 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + (-t515 / 0.4e1 - t523 / 0.4e1) * t533 + ((t107 * t92 + t117 * t158 + t118 * t156 + t155 * t84 + t157 * t83 + t56 * t82) * t527 + (t105 * t108 + t120 * t186 + t121 * t184 + t126 * t127 + t145 * t187 + t146 * t185) * t529) * t534 + (((-t114 / 0.2e1 - t112 / 0.2e1 + t38 / 0.2e1 + t37 / 0.2e1) * t336 + t365) * t332 + ((t113 / 0.2e1 + t111 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1) * t336 - t364) * t331 + ((t567 * t331 - t566 * t332) * t548 + t572 * t510 + (t571 * t331 - t570 * t332 + t573) * t509) * t334) * qJD(2), qJD(2) * t5 + qJD(5) * t94, t1 * qJD(2) + ((t104 * t92 + t117 * t124 + t118 * t125) * t526 + (t127 * t134 + t145 * t159 + t146 * t160) * t528) * t531 + t39 * t497 + (t300 / 0.2e1 + t301 / 0.2e1) * qJD(4) * t334 * t535 + (t542 * t510 + t543 * t547 + ((t577 * t333 - t578 * t335) * t334 + (t87 + t85) * t332 + (t88 + t86) * t331) * t507) * t421, t6 * qJD(2) + t39 * t498 + (t170 * t455 + t283 * t219 - t285 * t220 - t144) * t497 + t368; qJD(2) * t165 + qJD(4) * t162, (-t285 * t156 + t283 * t158 - t50 + (-t334 * t82 + (t107 + t376) * t336) * t335) * t499 + t91 * qJD(3) + t8 * qJD(4) + t539, qJD(2) * t91 - qJD(4) * t94, t8 * qJD(2) + (t284 * t117 - t286 * t118 + t283 * t124 - t285 * t125 + (t104 * t335 - t333 * t92) * t336 - t39) * t498 + t144 * t497 - t368, 0.4e1 * (t161 * qJD(2) / 0.4e1 + t144 * qJD(4) / 0.4e1) * m(6);];
Cq = t10;
