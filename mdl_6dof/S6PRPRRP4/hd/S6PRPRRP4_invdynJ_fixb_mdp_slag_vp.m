% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:31
% EndTime: 2019-03-08 20:12:41
% DurationCPUTime: 7.37s
% Computational Cost: add. (5185->509), mult. (11975->652), div. (0->0), fcn. (9816->14), ass. (0->219)
t478 = sin(pkin(11));
t484 = sin(qJ(4));
t481 = cos(pkin(11));
t613 = cos(qJ(4));
t554 = t613 * t481;
t517 = -t478 * t484 + t554;
t605 = pkin(8) + qJ(3);
t453 = t605 * t478;
t454 = t605 * t481;
t518 = -t453 * t613 - t454 * t484;
t367 = qJD(3) * t517 + qJD(4) * t518;
t480 = sin(pkin(6));
t487 = cos(qJ(2));
t582 = t480 * t487;
t498 = t517 * t582;
t414 = qJD(1) * t498;
t631 = -t367 + t414;
t485 = sin(qJ(2));
t567 = qJD(1) * t480;
t553 = t485 * t567;
t450 = qJD(2) * qJ(3) + t553;
t602 = cos(pkin(6));
t546 = qJD(1) * t602;
t463 = t481 * t546;
t604 = pkin(8) * qJD(2);
t406 = t463 + (-t450 - t604) * t478;
t421 = t450 * t481 + t478 * t546;
t407 = t481 * t604 + t421;
t353 = t406 * t484 + t407 * t613;
t630 = t353 * qJD(4);
t629 = t517 * qJD(2);
t444 = t517 * qJD(4);
t448 = t478 * t613 + t481 * t484;
t501 = t448 * qJDD(2);
t489 = qJD(2) * t444 + t501;
t628 = -qJD(4) * qJD(5) - t489;
t627 = t481 * MDP(5) - t478 * MDP(6);
t620 = t448 * qJD(2);
t626 = qJD(4) * t620;
t433 = qJD(5) - t629;
t568 = t478 ^ 2 + t481 ^ 2;
t625 = MDP(7) * t568;
t477 = pkin(11) + qJ(4);
t470 = sin(t477);
t479 = sin(pkin(10));
t601 = cos(pkin(10));
t538 = t602 * t601;
t437 = t479 * t485 - t487 * t538;
t548 = t479 * t602;
t439 = t485 * t601 + t487 * t548;
t541 = g(1) * t439 + g(2) * t437;
t504 = -g(3) * t582 + t541;
t624 = t504 * t470;
t552 = t487 * t567;
t559 = qJDD(2) * qJ(3);
t560 = qJDD(1) * t480;
t419 = t485 * t560 + t559 + (qJD(3) + t552) * qJD(2);
t545 = qJDD(1) * t602;
t461 = t481 * t545;
t373 = t461 + (-pkin(8) * qJDD(2) - t419) * t478;
t387 = t419 * t481 + t478 * t545;
t557 = qJDD(2) * t481;
t374 = pkin(8) * t557 + t387;
t521 = t373 * t484 + t374 * t613;
t621 = t406 * t613 - t407 * t484;
t322 = qJDD(4) * pkin(9) + qJD(4) * t621 + t521;
t445 = t448 * qJD(4);
t558 = qJDD(2) * t478;
t534 = -qJDD(2) * t554 + t484 * t558;
t395 = qJD(2) * t445 + t534;
t549 = t487 * t560;
t566 = qJD(2) * t485;
t550 = qJD(1) * t566;
t619 = t480 * t550 + qJDD(3);
t522 = -t549 + t619;
t599 = qJDD(2) * pkin(2);
t424 = t522 - t599;
t339 = -pkin(3) * t557 + t395 * pkin(4) - pkin(9) * t489 + t424;
t349 = qJD(4) * pkin(9) + t353;
t468 = t481 * pkin(3) + pkin(2);
t535 = qJD(3) - t552;
t432 = -qJD(2) * t468 + t535;
t361 = -pkin(4) * t629 - pkin(9) * t620 + t432;
t483 = sin(qJ(5));
t486 = cos(qJ(5));
t563 = qJD(5) * t486;
t564 = qJD(5) * t483;
t542 = t322 * t483 - t339 * t486 + t349 * t563 + t361 * t564;
t389 = qJDD(5) + t395;
t612 = pkin(5) * t389;
t316 = qJDD(6) + t542 - t612;
t330 = t349 * t486 + t361 * t483;
t325 = qJ(6) * t433 + t330;
t598 = t325 * t433;
t623 = -t316 + t598;
t412 = -t453 * t484 + t454 * t613;
t499 = t448 * t582;
t571 = -qJD(1) * t499 + qJD(3) * t448 + qJD(4) * t412;
t391 = pkin(4) * t445 - pkin(9) * t444;
t392 = -pkin(4) * t517 - pkin(9) * t448 - t468;
t622 = -t392 * t563 + t412 * t564 + t631 * t486 + (-t391 + t553) * t483;
t570 = t392 * t483 + t412 * t486;
t583 = t480 * t485;
t435 = -t478 * t583 + t481 * t602;
t436 = t478 * t602 + t481 * t583;
t378 = t435 * t484 + t436 * t613;
t459 = t483 * t582;
t364 = t378 * t486 - t459;
t438 = t479 * t487 + t485 * t538;
t440 = -t485 * t548 + t487 * t601;
t471 = cos(t477);
t547 = t480 * t601;
t584 = t479 * t480;
t508 = -g(3) * (-t470 * t583 + t471 * t602) - g(2) * (-t438 * t470 - t471 * t547) - g(1) * (-t440 * t470 + t471 * t584);
t348 = -qJD(4) * pkin(4) - t621;
t415 = -qJD(4) * t486 + t483 * t620;
t565 = qJD(4) * t483;
t417 = t486 * t620 + t565;
t333 = pkin(5) * t415 - qJ(6) * t417 + t348;
t611 = pkin(9) * t389;
t618 = t333 * t433 - t611;
t532 = (-t450 * t478 + t463) * t478 - t421 * t481;
t616 = t487 * t532 - (-qJD(2) * pkin(2) + t535) * t485;
t615 = t417 ^ 2;
t614 = t433 ^ 2;
t603 = pkin(9) * qJD(5);
t600 = qJ(6) * t389;
t597 = t330 * t433;
t346 = -qJDD(4) * t483 + t486 * t628 + t564 * t620;
t596 = t346 * t483;
t594 = t415 * t629;
t593 = t415 * t620;
t592 = t417 * t415;
t544 = t417 * t433;
t591 = t417 * t620;
t590 = t433 * t483;
t589 = t448 * t486;
t587 = t471 * t483;
t586 = t471 * t486;
t380 = t483 * t389;
t381 = t486 * t389;
t579 = t486 * t487;
t488 = qJD(2) ^ 2;
t578 = t487 * t488;
t577 = qJDD(1) - g(3);
t576 = qJ(6) * t445 - qJD(6) * t517 - t622;
t375 = t414 * t483 - t486 * t553;
t575 = -pkin(5) * t445 + qJD(5) * t570 + t367 * t483 - t391 * t486 - t375;
t536 = pkin(5) * t483 - qJ(6) * t486;
t537 = pkin(5) * t486 + qJ(6) * t483;
t574 = t536 * t444 + (qJD(5) * t537 - qJD(6) * t486) * t448 + t571;
t539 = -qJDD(4) * t486 + t563 * t620;
t347 = t483 * t501 + (qJD(5) + t629) * t565 + t539;
t573 = -t347 * t483 - t415 * t563;
t390 = pkin(4) * t620 - pkin(9) * t629;
t572 = t390 * t483 + t486 * t621;
t569 = -qJD(6) * t483 + t433 * t536 - t353;
t329 = -t349 * t483 + t361 * t486;
t562 = qJD(6) - t329;
t556 = t480 * t579;
t551 = t480 * t566;
t540 = g(1) * t440 + g(2) * t438;
t324 = -pkin(5) * t433 + t562;
t533 = t324 * t486 - t325 * t483;
t512 = t322 * t486 + t339 * t483 - t349 * t564 + t361 * t563;
t315 = qJD(6) * t433 + t512 + t600;
t531 = t324 * t433 + t315;
t529 = pkin(4) + t537;
t528 = pkin(4) * t471 + pkin(9) * t470 + t468;
t527 = -t424 + t541;
t525 = t380 + (-t486 * t629 + t563) * t433;
t524 = -t433 * t564 + t590 * t629 + t381;
t363 = t378 * t483 + t556;
t519 = t435 * t613 - t436 * t484;
t516 = t444 * t483 + t448 * t563;
t515 = -t444 * t486 + t448 * t564;
t514 = t348 * t433 - t611;
t513 = -t373 * t613 + t374 * t484 + t630;
t369 = -t437 * t587 - t438 * t486;
t371 = -t439 * t587 - t440 * t486;
t422 = t459 * t471 - t486 * t583;
t510 = g(1) * t371 + g(2) * t369 + g(3) * t422;
t370 = -t437 * t586 + t438 * t483;
t372 = -t439 * t586 + t440 * t483;
t423 = (t471 * t579 + t483 * t485) * t480;
t509 = -g(1) * t372 - g(2) * t370 - g(3) * t423;
t397 = t438 * t471 - t470 * t547;
t399 = t440 * t471 + t470 * t584;
t427 = t470 * t602 + t471 * t583;
t507 = -g(1) * t399 - g(2) * t397 - g(3) * t427;
t323 = -qJDD(4) * pkin(4) + t513;
t497 = t504 + t549;
t386 = -t419 * t478 + t461;
t496 = -t386 * t478 + t387 * t481 - t540;
t494 = -t433 * t603 + t508;
t357 = t397 * t483 - t437 * t486;
t359 = t399 * t483 - t439 * t486;
t400 = t427 * t483 + t556;
t493 = g(1) * t359 + g(2) * t357 + g(3) * t400 - t542;
t317 = pkin(5) * t347 + qJ(6) * t346 - qJD(6) * t417 + t323;
t492 = -t317 + t494;
t358 = t397 * t486 + t437 * t483;
t360 = t399 * t486 + t439 * t483;
t401 = t427 * t486 - t459;
t491 = -g(1) * t360 - g(2) * t358 - g(3) * t401 + t512;
t490 = t333 * t417 + qJDD(6) - t493;
t409 = -qJDD(2) * t468 + t522;
t362 = pkin(5) * t417 + qJ(6) * t415;
t354 = t448 * t536 - t518;
t351 = qJD(2) * t499 + qJD(4) * t378;
t350 = qJD(2) * t498 + qJD(4) * t519;
t342 = pkin(5) * t517 - t392 * t486 + t412 * t483;
t341 = -qJ(6) * t517 + t570;
t334 = t415 * t433 - t346;
t332 = -pkin(5) * t620 - t390 * t486 + t483 * t621;
t331 = qJ(6) * t620 + t572;
t328 = qJD(5) * t364 + t350 * t483 - t486 * t551;
t327 = -qJD(5) * t363 + t350 * t486 + t483 * t551;
t1 = [t577 * MDP(1) + (-t435 * t478 + t436 * t481) * qJDD(2) * MDP(7) + (t386 * t435 + t387 * t436 - g(3)) * MDP(8) + (-qJD(4) * t351 + qJDD(4) * t519) * MDP(14) + (-t350 * qJD(4) - t378 * qJDD(4)) * MDP(15) + (-t327 * t415 + t328 * t417 - t346 * t363 - t347 * t364) * MDP(24) + (t315 * t364 + t316 * t363 - t317 * t519 + t324 * t328 + t325 * t327 + t333 * t351 - g(3)) * MDP(26) + (MDP(21) + MDP(23)) * (-t328 * t433 - t347 * t519 + t351 * t415 - t363 * t389) + (-MDP(22) + MDP(25)) * (t327 * t433 - t346 * t519 - t351 * t417 + t364 * t389) + ((-qJDD(2) * t485 - t578) * MDP(4) + t578 * t625 + (-qJD(2) * t616 - t424 * t487) * MDP(8) + (-t395 * t487 - t566 * t629) * MDP(14) + (-t487 * t501 + (-t444 * t487 + t485 * t620) * qJD(2)) * MDP(15) + (MDP(3) + t627) * (qJDD(2) * t487 - t485 * t488)) * t480; qJDD(2) * MDP(2) + t497 * MDP(3) + (-t577 * t583 + t540) * MDP(4) + (-g(3) * t583 + t496 + (qJD(2) * t535 + t559) * t568) * MDP(7) + (-t532 * qJD(3) + t527 * pkin(2) + t496 * qJ(3) + (-g(3) * (pkin(2) * t487 + qJ(3) * t485) + t616 * qJD(1)) * t480) * MDP(8) + (t444 * t620 + t448 * t489) * MDP(9) + (-t448 * t395 + t444 * t629 - t445 * t620 + t489 * t517) * MDP(10) + (qJD(4) * t444 + qJDD(4) * t448) * MDP(11) + (-qJD(4) * t445 + qJDD(4) * t517) * MDP(12) + (-qJD(4) * t571 + qJDD(4) * t518 - t395 * t468 - t409 * t517 + t432 * t445 + t471 * t504 + t553 * t629) * MDP(14) + (qJD(4) * t631 - t412 * qJDD(4) + t409 * t448 + t432 * t444 - t468 * t489 - t620 * t553 - t624) * MDP(15) + (-t346 * t589 - t417 * t515) * MDP(16) + ((-t415 * t486 - t417 * t483) * t444 + (t596 - t347 * t486 + (t415 * t483 - t417 * t486) * qJD(5)) * t448) * MDP(17) + (t346 * t517 + t381 * t448 + t417 * t445 - t433 * t515) * MDP(18) + (t347 * t517 - t380 * t448 - t415 * t445 - t433 * t516) * MDP(19) + (-t389 * t517 + t433 * t445) * MDP(20) + (t542 * t517 + t329 * t445 - t518 * t347 + t375 * t433 + t571 * t415 + ((-qJD(5) * t412 + t391) * t433 + t392 * t389 + t348 * qJD(5) * t448) * t486 + ((-qJD(5) * t392 - t367) * t433 - t412 * t389 + t323 * t448 + t348 * t444) * t483 + t509) * MDP(21) + (t323 * t589 - t330 * t445 + t346 * t518 - t515 * t348 - t570 * t389 + t571 * t417 + t433 * t622 + t512 * t517 + t510) * MDP(22) + (t317 * t448 * t483 + t316 * t517 - t324 * t445 + t333 * t516 - t342 * t389 + t347 * t354 + t415 * t574 - t433 * t575 + t509) * MDP(23) + (-t341 * t347 - t342 * t346 + t533 * t444 + t575 * t417 - t576 * t415 + t624 + (-t315 * t483 + t316 * t486 + (-t324 * t483 - t325 * t486) * qJD(5)) * t448) * MDP(24) + (-t315 * t517 - t317 * t589 + t325 * t445 + t333 * t515 + t341 * t389 + t346 * t354 - t417 * t574 + t433 * t576 - t510) * MDP(25) + (t315 * t341 + t317 * t354 + t316 * t342 - g(1) * (pkin(5) * t372 + qJ(6) * t371 + t440 * t605) - g(2) * (pkin(5) * t370 + qJ(6) * t369 + t438 * t605) + t574 * t333 + t576 * t325 + t575 * t324 + t541 * t528 + (-pkin(5) * t423 - qJ(6) * t422 - (t485 * t605 + t487 * t528) * t480) * g(3)) * MDP(26) + t627 * (t480 * (-g(3) * t487 + t550) + t527 + t599); -MDP(5) * t557 + MDP(6) * t558 - t488 * t625 + (qJD(2) * t532 - t497 - t599 + t619) * MDP(8) + (t534 + 0.2e1 * t626) * MDP(14) + (0.2e1 * qJD(4) * t629 + t501) * MDP(15) + (t524 - t593) * MDP(21) + (-t486 * t614 - t380 - t591) * MDP(22) + (-t433 * t590 + t381 - t593) * MDP(23) + ((t346 + t594) * t486 + t483 * t544 + t573) * MDP(24) + (t525 + t591) * MDP(25) + (-t333 * t620 + t531 * t483 + t486 * t623 - t504) * MDP(26); -t629 ^ 2 * MDP(10) + t501 * MDP(11) - t534 * MDP(12) + qJDD(4) * MDP(13) + (t508 - t513 + t630) * MDP(14) + (-t432 * t629 - t507 - t521) * MDP(15) + (t486 * t544 - t596) * MDP(16) + ((-t346 + t594) * t486 - t417 * t590 + t573) * MDP(17) + (t525 - t591) * MDP(18) + (t524 + t593) * MDP(19) + (-pkin(4) * t347 - t353 * t415 + (t433 * t621 + t514) * t483 + (-t323 + (-t390 - t603) * t433 + t508) * t486) * MDP(21) + (pkin(4) * t346 + t572 * t433 - t353 * t417 + t514 * t486 + (t323 - t494) * t483) * MDP(22) + (t332 * t433 - t347 * t529 + t569 * t415 + t483 * t618 + t492 * t486) * MDP(23) + (t331 * t415 - t332 * t417 + ((qJD(5) * t417 - t347) * pkin(9) + t531) * t486 + ((qJD(5) * t415 - t346) * pkin(9) - t623) * t483 + t507) * MDP(24) + (-t331 * t433 - t346 * t529 - t569 * t417 + t492 * t483 - t486 * t618) * MDP(25) + (-t324 * t332 - t325 * t331 + t569 * t333 + (qJD(5) * t533 + t315 * t486 + t316 * t483 + t507) * pkin(9) + (-t317 + t508) * t529) * MDP(26) + (MDP(10) * t620 - MDP(14) * t432 - MDP(20) * t433 - MDP(21) * t329 + MDP(22) * t330 + MDP(23) * t324 - MDP(25) * t325 - MDP(9) * t629) * t620; MDP(16) * t592 + (-t415 ^ 2 + t615) * MDP(17) + t334 * MDP(18) + (t483 * t628 - t539 + t544) * MDP(19) + t389 * MDP(20) + (-t348 * t417 + t493 + t597) * MDP(21) + (t329 * t433 + t348 * t415 - t491) * MDP(22) + (-t362 * t415 - t490 + t597 + 0.2e1 * t612) * MDP(23) + (pkin(5) * t346 - qJ(6) * t347 + (t325 - t330) * t417 + (t324 - t562) * t415) * MDP(24) + (0.2e1 * t600 - t333 * t415 + t362 * t417 + (0.2e1 * qJD(6) - t329) * t433 + t491) * MDP(25) + (t315 * qJ(6) - t316 * pkin(5) - t333 * t362 - t324 * t330 - g(1) * (-pkin(5) * t359 + qJ(6) * t360) - g(2) * (-pkin(5) * t357 + qJ(6) * t358) - g(3) * (-pkin(5) * t400 + qJ(6) * t401) + t562 * t325) * MDP(26); (-qJDD(5) - t534 + t592 - t626) * MDP(23) + t334 * MDP(24) + (-t614 - t615) * MDP(25) + (t490 - t598 - t612) * MDP(26);];
tau  = t1;
