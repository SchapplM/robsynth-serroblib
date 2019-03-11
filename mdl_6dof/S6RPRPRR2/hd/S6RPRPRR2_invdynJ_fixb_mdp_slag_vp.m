% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:39:10
% EndTime: 2019-03-09 03:39:18
% DurationCPUTime: 6.42s
% Computational Cost: add. (4853->472), mult. (10626->616), div. (0->0), fcn. (7924->18), ass. (0->229)
t487 = sin(pkin(11));
t494 = sin(qJ(3));
t560 = qJD(1) * t494;
t489 = cos(pkin(11));
t498 = cos(qJ(3));
t573 = t489 * t498;
t437 = qJD(1) * t573 - t487 * t560;
t607 = qJD(5) + qJD(6);
t617 = t437 - t607;
t477 = t498 * qJDD(2);
t488 = sin(pkin(10));
t465 = pkin(1) * t488 + pkin(7);
t455 = t465 * qJDD(1);
t504 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t455;
t571 = qJ(4) + t465;
t535 = t571 * qJD(1);
t518 = t535 * qJD(3);
t356 = qJDD(3) * pkin(3) - t494 * t504 - t498 * t518 + t477;
t360 = (qJDD(2) - t518) * t494 + t504 * t498;
t328 = t356 * t489 - t487 * t360;
t326 = -qJDD(3) * pkin(4) - t328;
t433 = qJD(5) - t437;
t464 = pkin(3) * t487 + pkin(8);
t482 = qJ(3) + pkin(11);
t471 = sin(t482);
t473 = cos(t482);
t483 = qJ(1) + pkin(10);
t472 = sin(t483);
t474 = cos(t483);
t531 = g(1) * t474 + g(2) * t472;
t507 = -g(3) * t473 + t471 * t531;
t616 = -qJD(5) * t464 * t433 - t326 + t507;
t447 = t487 * t498 + t489 * t494;
t439 = t447 * qJD(1);
t493 = sin(qJ(5));
t497 = cos(qJ(5));
t555 = t497 * qJD(3);
t418 = t439 * t493 - t555;
t496 = cos(qJ(6));
t420 = qJD(3) * t493 + t439 * t497;
t492 = sin(qJ(6));
t588 = t420 * t492;
t362 = t496 * t418 + t588;
t431 = qJD(6) + t433;
t615 = t362 * t431;
t520 = t418 * t492 - t496 * t420;
t614 = t431 * t520;
t572 = t492 * t497;
t449 = t493 * t496 + t572;
t563 = t617 * t449;
t559 = qJD(5) * t493;
t587 = t437 * t493;
t613 = t559 - t587;
t421 = t498 * qJD(2) - t494 * t535;
t600 = qJD(3) * pkin(3);
t412 = t421 + t600;
t422 = qJD(2) * t494 + t498 * t535;
t574 = t489 * t422;
t359 = t487 * t412 + t574;
t353 = qJD(3) * pkin(8) + t359;
t470 = pkin(3) * t498 + pkin(2);
t490 = cos(pkin(10));
t605 = pkin(1) * t490;
t519 = -t470 - t605;
t436 = qJD(1) * t519 + qJD(4);
t374 = -pkin(4) * t437 - pkin(8) * t439 + t436;
t334 = t353 * t497 + t374 * t493;
t322 = -pkin(9) * t418 + t334;
t557 = qJD(6) * t492;
t320 = t322 * t557;
t409 = t487 * t422;
t358 = t412 * t489 - t409;
t352 = -qJD(3) * pkin(4) - t358;
t341 = pkin(5) * t418 + t352;
t486 = qJ(5) + qJ(6);
t480 = sin(t486);
t578 = t474 * t480;
t481 = cos(t486);
t581 = t472 * t481;
t414 = -t473 * t581 + t578;
t577 = t474 * t481;
t582 = t472 * t480;
t416 = t473 * t577 + t582;
t604 = g(3) * t471;
t612 = g(1) * t416 - g(2) * t414 + t341 * t362 + t481 * t604 + t320;
t413 = t473 * t582 + t577;
t415 = -t473 * t578 + t581;
t329 = t487 * t356 + t489 * t360;
t327 = qJDD(3) * pkin(8) + t329;
t552 = qJDD(1) * t498;
t553 = qJDD(1) * t494;
t524 = -t487 * t553 + t489 * t552;
t401 = -qJD(3) * t439 + t524;
t554 = qJD(1) * qJD(3);
t542 = t498 * t554;
t543 = t494 * t554;
t402 = qJDD(1) * t447 - t487 * t543 + t489 * t542;
t608 = pkin(3) * t543 + t519 * qJDD(1) + qJDD(4);
t339 = -pkin(4) * t401 - pkin(8) * t402 + t608;
t338 = t497 * t339;
t349 = qJD(5) * t555 + t493 * qJDD(3) + t497 * t402 - t439 * t559;
t438 = t447 * qJD(3);
t400 = qJD(1) * t438 + qJDD(5) - t524;
t309 = pkin(5) * t400 - pkin(9) * t349 - qJD(5) * t334 - t327 * t493 + t338;
t350 = t420 * qJD(5) - t497 * qJDD(3) + t402 * t493;
t558 = qJD(5) * t497;
t513 = t497 * t327 + t493 * t339 - t353 * t559 + t374 * t558;
t310 = -pkin(9) * t350 + t513;
t540 = t496 * t309 - t492 * t310;
t611 = -g(1) * t415 + g(2) * t413 + t341 * t520 + t480 * t604 + t540;
t396 = qJDD(6) + t400;
t610 = t396 * MDP(25) + (-t362 ^ 2 + t520 ^ 2) * MDP(22) - t362 * MDP(21) * t520;
t393 = t449 * t447;
t448 = t492 * t493 - t496 * t497;
t564 = t617 * t448;
t606 = -t396 * t449 - t431 * t564;
t539 = t349 * t492 + t496 * t350;
t317 = -qJD(6) * t520 + t539;
t602 = g(3) * t498;
t601 = pkin(9) + t464;
t333 = -t353 * t493 + t497 * t374;
t321 = -pkin(9) * t420 + t333;
t319 = pkin(5) * t433 + t321;
t599 = t319 * t496;
t598 = t322 * t496;
t597 = t349 * t493;
t596 = t362 * t439;
t595 = t520 * t439;
t593 = t400 * t493;
t592 = t418 * t433;
t591 = t418 * t439;
t590 = t420 * t433;
t589 = t420 * t439;
t446 = t487 * t494 - t573;
t441 = t446 * qJD(3);
t586 = t441 * t493;
t585 = t441 * t497;
t584 = t447 * t493;
t583 = t447 * t497;
t580 = t472 * t493;
t579 = t472 * t497;
t576 = t474 * t493;
t575 = t474 * t497;
t442 = t571 * t494;
t443 = t571 * t498;
t398 = -t442 * t487 + t443 * t489;
t389 = t497 * t398;
t392 = t497 * t400;
t570 = qJDD(2) - g(3);
t556 = qJD(6) * t496;
t547 = t496 * t349 - t492 * t350 - t418 * t556;
t316 = -t420 * t557 + t547;
t569 = t316 * t446 - t438 * t520;
t546 = t447 * t559;
t331 = -t441 * t572 - t492 * t546 - t557 * t584 + (t583 * t607 - t586) * t496;
t568 = -t331 * t431 - t393 * t396;
t567 = t349 * t446 + t420 * t438;
t367 = t421 * t489 - t409;
t386 = pkin(3) * t560 + pkin(4) * t439 - pkin(8) * t437;
t566 = t497 * t367 + t493 * t386;
t390 = pkin(4) * t446 - pkin(8) * t447 + t519;
t565 = t493 * t390 + t389;
t484 = t494 ^ 2;
t562 = -t498 ^ 2 + t484;
t467 = -pkin(2) - t605;
t458 = qJD(1) * t467;
t550 = t494 * t600;
t549 = t400 * t584;
t548 = t447 * t392;
t466 = -pkin(3) * t489 - pkin(4);
t545 = t447 * t558;
t541 = qJD(5) * t601;
t366 = t421 * t487 + t574;
t537 = qJD(3) * t571;
t423 = qJD(4) * t498 - t494 * t537;
t424 = -qJD(4) * t494 - t498 * t537;
t372 = t423 * t487 - t489 * t424;
t397 = t489 * t442 + t443 * t487;
t536 = t433 * t497;
t534 = -qJD(5) * t374 - t327;
t533 = qJD(6) * t319 + t310;
t532 = pkin(5) * t613 - t366;
t530 = g(1) * t472 - g(2) * t474;
t495 = sin(qJ(1));
t499 = cos(qJ(1));
t529 = g(1) * t495 - g(2) * t499;
t528 = -t448 * t396 + t431 * t563;
t527 = -t353 * t558 + t338;
t379 = t497 * t386;
t445 = t601 * t497;
t526 = pkin(5) * t439 + qJD(6) * t445 - t367 * t493 + t379 + (-pkin(9) * t437 + t541) * t497;
t444 = t601 * t493;
t525 = -pkin(9) * t587 + qJD(6) * t444 + t493 * t541 + t566;
t523 = -t317 * t446 - t362 * t438;
t312 = t319 * t492 + t598;
t330 = -t393 * t607 + t448 * t441;
t394 = t448 * t447;
t522 = -t330 * t431 + t394 * t396;
t521 = -t350 * t446 - t418 * t438;
t517 = -t433 * t613 + t392;
t516 = t545 - t586;
t515 = t546 + t585;
t373 = t423 * t489 + t424 * t487;
t387 = pkin(4) * t438 + pkin(8) * t441 + t550;
t512 = t497 * t373 + t493 * t387 + t390 * t558 - t398 * t559;
t511 = t352 * t433 - t464 * t400;
t509 = -qJD(1) * t458 - t455 + t531;
t508 = 0.2e1 * qJD(3) * t458 - qJDD(3) * t465;
t500 = qJD(3) ^ 2;
t503 = -0.2e1 * qJDD(1) * t467 - t465 * t500 + t530;
t491 = -qJ(4) - pkin(7);
t454 = qJDD(3) * t498 - t494 * t500;
t453 = qJDD(3) * t494 + t498 * t500;
t452 = -pkin(5) * t497 + t466;
t428 = t473 * t575 + t580;
t427 = -t473 * t576 + t579;
t426 = -t473 * t579 + t576;
t425 = t473 * t580 + t575;
t385 = t497 * t390;
t380 = t497 * t387;
t371 = pkin(5) * t584 + t397;
t340 = pkin(5) * t516 + t372;
t336 = -pkin(9) * t584 + t565;
t335 = pkin(5) * t446 - pkin(9) * t583 - t398 * t493 + t385;
t318 = pkin(5) * t350 + t326;
t314 = -pkin(9) * t516 + t512;
t313 = pkin(9) * t585 + pkin(5) * t438 - t373 * t493 + t380 + (-t389 + (pkin(9) * t447 - t390) * t493) * qJD(5);
t311 = -t322 * t492 + t599;
t1 = [(t494 * t508 + t498 * t503) * MDP(10) + (-t494 * t503 + t498 * t508) * MDP(11) + (-g(1) * t413 - g(2) * t415 - t312 * t438 + t371 * t316 - t318 * t394 + t320 * t446 + t341 * t330 - t340 * t520 + (-(-qJD(6) * t336 + t313) * t431 - t335 * t396 - t309 * t446) * t492 + (-(qJD(6) * t335 + t314) * t431 - t336 * t396 - t533 * t446) * t496) * MDP(27) + (-t316 * t394 - t330 * t520) * MDP(21) + (-t316 * t393 + t317 * t394 - t330 * t362 + t331 * t520) * MDP(22) + (t329 * t398 + t359 * t373 - t328 * t397 - t358 * t372 + t436 * t550 - g(1) * (-pkin(1) * t495 - t470 * t472 - t474 * t491) - g(2) * (pkin(1) * t499 + t470 * t474 - t472 * t491) + t608 * t519) * MDP(13) + ((-t398 * t558 + t380) * t433 + t385 * t400 + t527 * t446 + t333 * t438 + t372 * t418 + t397 * t350 + t352 * t545 - g(1) * t426 - g(2) * t428 + ((-qJD(5) * t390 - t373) * t433 - t398 * t400 + t534 * t446 + t326 * t447 - t352 * t441) * t493) * MDP(19) + (-(-t418 * t497 - t420 * t493) * t441 + (-t597 - t350 * t497 + (t418 * t493 - t420 * t497) * qJD(5)) * t447) * MDP(15) + (t529 + (t488 ^ 2 + t490 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (-t328 * t447 - t329 * t446 + t358 * t441 - t359 * t438 + t372 * t439 + t373 * t437 + t397 * t402 + t398 * t401 - t531) * MDP(12) + t529 * MDP(2) + ((t313 * t496 - t314 * t492) * t431 + (t335 * t496 - t336 * t492) * t396 + t540 * t446 + t311 * t438 + t340 * t362 + t371 * t317 + t318 * t393 + t341 * t331 - g(1) * t414 - g(2) * t416 + ((-t335 * t492 - t336 * t496) * t431 - t312 * t446) * qJD(6)) * MDP(26) + (qJDD(1) * t484 + 0.2e1 * t494 * t542) * MDP(5) + (-t433 * t516 + t521 - t549) * MDP(17) + 0.2e1 * (t494 * t552 - t554 * t562) * MDP(6) + (-t433 * t515 + t548 + t567) * MDP(16) + (t523 + t568) * MDP(24) + (-t522 + t569) * MDP(23) + (t349 * t583 - t420 * t515) * MDP(14) + (-g(1) * t425 - g(2) * t427 + t326 * t583 - t334 * t438 + t397 * t349 - t352 * t515 + t372 * t420 - t400 * t565 - t433 * t512 - t446 * t513) * MDP(20) + qJDD(1) * MDP(1) + (t396 * t446 + t431 * t438) * MDP(25) + (t400 * t446 + t433 * t438) * MDP(18) + t453 * MDP(7) + t454 * MDP(8) + (g(1) * t499 + g(2) * t495) * MDP(3); t570 * MDP(4) + t454 * MDP(10) - t453 * MDP(11) + (t401 * t447 + t402 * t446 - t437 * t441 + t438 * t439) * MDP(12) + (-t328 * t446 + t329 * t447 - t358 * t438 - t359 * t441 - g(3)) * MDP(13) + (-t521 - t549) * MDP(19) + (-t548 + t567) * MDP(20) + (-t523 + t568) * MDP(26) + (t522 + t569) * MDP(27) + (-MDP(19) * t516 + t515 * MDP(20)) * t433; MDP(7) * t553 + MDP(8) * t552 + qJDD(3) * MDP(9) + (t494 * t509 + t477 - t602) * MDP(10) + (-t494 * t570 + t498 * t509) * MDP(11) + ((t358 - t367) * t437 + (t401 * t487 - t402 * t489) * pkin(3)) * MDP(12) + (t358 * t366 - t359 * t367 + (-t602 + t328 * t489 + t329 * t487 + (-qJD(1) * t436 + t531) * t494) * pkin(3)) * MDP(13) + (t420 * t536 + t597) * MDP(14) + ((t349 - t592) * t497 + (-t350 - t590) * t493) * MDP(15) + (t433 * t536 - t589 + t593) * MDP(16) + (t517 + t591) * MDP(17) + (t466 * t350 - t366 * t418 - t379 * t433 + (t367 * t433 + t511) * t493 + t616 * t497) * MDP(19) + (t466 * t349 - t366 * t420 + t566 * t433 - t616 * t493 + t511 * t497) * MDP(20) + (t316 * t449 - t520 * t564) * MDP(21) + (-t316 * t448 - t317 * t449 - t362 * t564 - t520 * t563) * MDP(22) + (t595 - t606) * MDP(23) + (t528 + t596) * MDP(24) + ((-t444 * t496 - t445 * t492) * t396 + t452 * t317 + t318 * t448 + (t492 * t525 - t496 * t526) * t431 + t532 * t362 - t563 * t341 + t507 * t481) * MDP(26) + (-(-t444 * t492 + t445 * t496) * t396 + t452 * t316 + t318 * t449 + (t492 * t526 + t496 * t525) * t431 - t532 * t520 + t564 * t341 - t507 * t480) * MDP(27) + (-MDP(5) * t494 * t498 + MDP(6) * t562) * qJD(1) ^ 2 - ((-t359 + t366) * MDP(12) + t433 * MDP(18) + t333 * MDP(19) - t334 * MDP(20) + t431 * MDP(25) + t311 * MDP(26) - t312 * MDP(27)) * t439; (-t437 ^ 2 - t439 ^ 2) * MDP(12) + (t517 - t591) * MDP(19) + (-t433 ^ 2 * t497 - t589 - t593) * MDP(20) + (t528 - t596) * MDP(26) + (t595 + t606) * MDP(27) + (t358 * t439 - t359 * t437 - t530 + t608) * MDP(13); t420 * t418 * MDP(14) + (-t418 ^ 2 + t420 ^ 2) * MDP(15) + (t349 + t592) * MDP(16) + (-t350 + t590) * MDP(17) + t400 * MDP(18) + (-g(1) * t427 + g(2) * t425 + t334 * t433 - t352 * t420 + (t534 + t604) * t493 + t527) * MDP(19) + (g(1) * t428 - g(2) * t426 + t333 * t433 + t352 * t418 + t497 * t604 - t513) * MDP(20) + (t316 + t615) * MDP(23) + (-t317 - t614) * MDP(24) + (-(-t321 * t492 - t598) * t431 - t312 * qJD(6) + (-t362 * t420 + t396 * t496 - t431 * t557) * pkin(5) + t611) * MDP(26) + ((-t322 * t431 - t309) * t492 + (t321 * t431 - t533) * t496 + (-t396 * t492 + t420 * t520 - t431 * t556) * pkin(5) + t612) * MDP(27) + t610; (t547 + t615) * MDP(23) + (-t539 - t614) * MDP(24) + (t312 * t431 + t611) * MDP(26) + (-t492 * t309 - t496 * t310 + t311 * t431 + t612) * MDP(27) + (-MDP(23) * t588 + MDP(24) * t520 - MDP(26) * t312 - MDP(27) * t599) * qJD(6) + t610;];
tau  = t1;
