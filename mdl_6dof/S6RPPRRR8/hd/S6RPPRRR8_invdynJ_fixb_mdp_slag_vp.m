% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:31
% EndTime: 2019-03-09 02:36:40
% DurationCPUTime: 6.29s
% Computational Cost: add. (4418->492), mult. (9160->624), div. (0->0), fcn. (6976->14), ass. (0->220)
t470 = sin(pkin(10));
t475 = sin(qJ(4));
t541 = qJD(1) * t475;
t520 = t470 * t541;
t471 = cos(pkin(10));
t479 = cos(qJ(4));
t540 = qJD(1) * t479;
t525 = t471 * t540;
t415 = -t520 + t525;
t474 = sin(qJ(5));
t478 = cos(qJ(5));
t532 = t478 * qJD(4);
t388 = t415 * t474 - t532;
t421 = t470 * t479 + t471 * t475;
t413 = t421 * qJD(1);
t603 = qJD(5) + t413;
t609 = t388 * t603;
t472 = -pkin(1) - qJ(3);
t586 = -qJD(1) * qJD(3) + qJDD(1) * t472;
t424 = qJDD(2) + t586;
t516 = -pkin(7) * qJDD(1) + t424;
t393 = t516 * t470;
t394 = t516 * t471;
t436 = qJD(1) * t472 + qJD(2);
t519 = -pkin(7) * qJD(1) + t436;
t405 = t519 * t470;
t406 = t519 * t471;
t538 = qJD(4) * t479;
t539 = qJD(4) * t475;
t495 = -t475 * t393 + t394 * t479 - t405 * t538 - t406 * t539;
t323 = -qJDD(4) * pkin(4) - t495;
t465 = pkin(10) + qJ(4);
t450 = sin(t465);
t451 = cos(t465);
t476 = sin(qJ(1));
t480 = cos(qJ(1));
t594 = g(1) * t476 - g(2) * t480;
t485 = -g(3) * t450 + t451 * t594;
t608 = pkin(8) * qJD(5) * t603 + t323 + t485;
t477 = cos(qJ(6));
t390 = qJD(4) * t474 + t415 * t478;
t473 = sin(qJ(6));
t571 = t390 * t473;
t339 = t477 * t388 + t571;
t404 = qJD(6) + t603;
t607 = t339 * t404;
t499 = t388 * t473 - t477 * t390;
t606 = t404 * t499;
t557 = t473 * t478;
t423 = t474 * t477 + t557;
t589 = qJD(5) + qJD(6);
t383 = t589 * t423;
t547 = t423 * t413 + t383;
t537 = qJD(5) * t474;
t568 = t413 * t474;
t605 = t537 + t568;
t512 = t478 * t603;
t524 = t471 * t538;
t435 = qJD(4) * t520;
t587 = t421 * qJDD(1) - t435;
t378 = qJD(1) * t524 + t587;
t374 = qJDD(5) + t378;
t556 = t474 * t374;
t604 = -t512 * t603 - t556;
t358 = t479 * t405 + t475 * t406;
t353 = qJD(4) * pkin(8) + t358;
t449 = qJD(1) * qJ(2) + qJD(3);
t457 = t470 * pkin(3);
t429 = qJD(1) * t457 + t449;
t354 = pkin(4) * t413 - pkin(8) * t415 + t429;
t325 = t353 * t478 + t354 * t474;
t315 = -pkin(9) * t388 + t325;
t535 = qJD(6) * t473;
t313 = t315 * t535;
t596 = -t475 * t405 + t406 * t479;
t352 = -qJD(4) * pkin(4) - t596;
t331 = pkin(5) * t388 + t352;
t469 = qJ(5) + qJ(6);
t456 = cos(t469);
t561 = t456 * t476;
t455 = sin(t469);
t562 = t455 * t480;
t399 = t450 * t561 + t562;
t560 = t456 * t480;
t563 = t455 * t476;
t401 = t450 * t560 - t563;
t582 = g(3) * t451;
t602 = g(1) * t399 - g(2) * t401 + t331 * t339 + t456 * t582 + t313;
t398 = -t450 * t563 + t560;
t400 = t450 * t562 + t561;
t498 = t393 * t479 + t394 * t475;
t322 = qJDD(4) * pkin(8) + qJD(4) * t596 + t498;
t551 = t479 * t471;
t437 = qJDD(1) * t551;
t529 = qJDD(1) * t470;
t500 = -t475 * t529 + t437;
t377 = -qJD(4) * t413 + t500;
t466 = qJDD(1) * qJ(2);
t467 = qJD(1) * qJD(2);
t593 = t466 + t467;
t432 = qJDD(3) + t593;
t419 = pkin(3) * t529 + t432;
t330 = pkin(4) * t378 - pkin(8) * t377 + t419;
t328 = t478 * t330;
t335 = qJD(5) * t532 + t474 * qJDD(4) + t478 * t377 - t415 * t537;
t303 = pkin(5) * t374 - pkin(9) * t335 - qJD(5) * t325 - t322 * t474 + t328;
t513 = -t478 * qJDD(4) + t377 * t474;
t336 = qJD(5) * t390 + t513;
t536 = qJD(5) * t478;
t487 = t478 * t322 + t474 * t330 - t353 * t537 + t354 * t536;
t304 = -pkin(9) * t336 + t487;
t515 = t477 * t303 - t473 * t304;
t601 = -g(1) * t398 - g(2) * t400 + t331 * t499 + t455 * t582 + t515;
t371 = qJDD(6) + t374;
t600 = t371 * MDP(29) + (-t339 ^ 2 + t499 ^ 2) * MDP(26) - t339 * MDP(25) * t499;
t543 = t470 ^ 2 + t471 ^ 2;
t598 = t436 * t543;
t506 = g(1) * t480 + g(2) * t476;
t595 = t432 - t506;
t592 = t473 * t537 + t474 * t535;
t591 = -qJD(6) * t478 - t536;
t590 = t470 * MDP(7) + t471 * MDP(8);
t422 = t473 * t474 - t477 * t478;
t548 = (-t413 - t589) * t422;
t574 = t371 * t423;
t588 = -t404 * t548 - t574;
t514 = t335 * t473 + t477 * t336;
t310 = -qJD(6) * t499 + t514;
t585 = 0.2e1 * t467;
t584 = pkin(8) + pkin(9);
t581 = -pkin(7) + t472;
t580 = pkin(1) * qJDD(1);
t324 = -t353 * t474 + t478 * t354;
t314 = -pkin(9) * t390 + t324;
t312 = pkin(5) * t603 + t314;
t579 = t312 * t477;
t578 = t315 * t477;
t577 = t335 * t474;
t576 = t339 * t415;
t575 = t499 * t415;
t573 = t388 * t415;
t572 = t390 * t415;
t417 = -t470 * t539 + t524;
t570 = t404 * t417;
t416 = -t470 * t538 - t471 * t539;
t567 = t416 * t474;
t566 = t416 * t478;
t420 = t470 * t475 - t551;
t565 = t420 * t474;
t564 = t420 * t478;
t345 = t422 * t371;
t555 = t474 * t476;
t554 = t474 * t480;
t553 = t476 * t478;
t364 = t478 * t374;
t427 = t581 * t470;
t428 = t581 * t471;
t381 = t427 * t479 + t428 * t475;
t379 = t478 * t381;
t552 = t478 * t480;
t442 = qJ(2) + t457;
t376 = pkin(4) * t415 + pkin(8) * t413;
t550 = t474 * t376 + t478 * t596;
t375 = pkin(4) * t421 + pkin(8) * t420 + t442;
t549 = t474 * t375 + t379;
t546 = t416 * qJD(4) - t420 * qJDD(4);
t545 = t480 * pkin(1) + t476 * qJ(2);
t534 = qJD(6) * t477;
t527 = t477 * t335 - t473 * t336 - t388 * t534;
t526 = qJD(5) * t584;
t523 = t420 * t536;
t518 = t543 * MDP(9);
t517 = t543 * t424;
t511 = -qJD(5) * t354 - t322;
t510 = qJD(6) * t312 + t304;
t509 = qJDD(2) - t580;
t508 = qJD(5) * t421 + qJD(1);
t507 = pkin(5) * t605 - t358;
t504 = -t404 * t547 - t345;
t503 = -t353 * t536 + t328;
t367 = t478 * t376;
t434 = t584 * t478;
t502 = pkin(5) * t415 + qJD(6) * t434 - t596 * t474 + t367 + (pkin(9) * t413 + t526) * t478;
t433 = t584 * t474;
t501 = pkin(9) * t568 + qJD(6) * t433 + t474 * t526 + t550;
t306 = t312 * t473 + t578;
t380 = t427 * t475 - t428 * t479;
t496 = -qJD(4) * t417 - qJDD(4) * t421;
t494 = -t603 * t605 + t364;
t493 = t523 - t567;
t492 = t420 * t537 + t566;
t491 = -pkin(8) * t374 + t352 * t603;
t489 = t404 * t422;
t348 = -qJD(3) * t421 - qJD(4) * t380;
t372 = pkin(4) * t417 - pkin(8) * t416 + qJD(2);
t486 = t478 * t348 + t474 * t372 + t375 * t536 - t381 * t537;
t309 = -t390 * t535 + t527;
t349 = -qJD(3) * t420 + qJD(4) * t381;
t481 = qJD(1) ^ 2;
t459 = t480 * qJ(2);
t448 = -pkin(5) * t478 - pkin(4);
t412 = t450 * t552 - t555;
t411 = t450 * t554 + t553;
t410 = t450 * t553 + t554;
t409 = -t450 * t555 + t552;
t369 = t422 * t420;
t368 = t423 * t420;
t365 = t478 * t375;
t361 = t478 * t372;
t350 = -pkin(5) * t565 + t380;
t329 = -pkin(5) * t493 + t349;
t326 = pkin(9) * t565 + t549;
t320 = pkin(5) * t421 + pkin(9) * t564 - t381 * t474 + t365;
t317 = t416 * t557 + (-t564 * t589 + t567) * t477 + t592 * t420;
t316 = t383 * t420 - t416 * t422;
t311 = pkin(5) * t336 + t323;
t308 = pkin(9) * t493 + t486;
t307 = -pkin(9) * t566 + pkin(5) * t417 - t348 * t474 + t361 + (-t379 + (-pkin(9) * t420 - t375) * t474) * qJD(5);
t305 = -t315 * t473 + t579;
t1 = [t496 * MDP(14) + (t432 * qJ(2) + t449 * qJD(2) - g(1) * (t472 * t476 + t459) - g(2) * (qJ(3) * t480 + t545) + t472 * t517 - qJD(3) * t598) * MDP(10) + t590 * (t593 + t595) + (-t377 * t420 + t415 * t416) * MDP(11) + (-t377 * t421 + t378 * t420 - t413 * t416 - t415 * t417) * MDP(12) + (qJD(2) * t413 - qJD(4) * t349 - qJDD(4) * t380 + t378 * t442 + t417 * t429 + t419 * t421 - t450 * t506) * MDP(16) + (-t310 * t421 - t317 * t404 - t339 * t417 + t368 * t371) * MDP(28) + ((t307 * t477 - t308 * t473) * t404 + (t320 * t477 - t326 * t473) * t371 + t515 * t421 + t305 * t417 + t329 * t339 + t350 * t310 - t311 * t368 + t331 * t317 - g(1) * t401 - g(2) * t399 + ((-t320 * t473 - t326 * t477) * t404 - t306 * t421) * qJD(6)) * MDP(30) + (t374 * t421 + t417 * t603) * MDP(22) + ((-t381 * t536 + t361) * t603 + t365 * t374 + t503 * t421 + t324 * t417 + t349 * t388 + t380 * t336 - t352 * t523 - g(1) * t412 - g(2) * t410 + ((-qJD(5) * t375 - t348) * t603 - t381 * t374 + t511 * t421 - t323 * t420 + t352 * t416) * t474) * MDP(23) + (t335 * t421 - t364 * t420 + t390 * t417 + t492 * t603) * MDP(20) + (-t336 * t421 - t388 * t417 + t420 * t556 + t493 * t603) * MDP(21) + (g(1) * t411 - g(2) * t409 - t323 * t564 - t325 * t417 + t380 * t335 + t349 * t390 + t352 * t492 - t374 * t549 - t421 * t487 - t486 * t603) * MDP(24) + (t309 * t421 + t316 * t404 + t369 * t371 - t417 * t499) * MDP(27) + (t309 * t369 - t316 * t499) * MDP(25) + (g(1) * t400 - g(2) * t398 - t306 * t417 + t350 * t309 + t311 * t369 + t313 * t421 + t331 * t316 - t329 * t499 + (-(-qJD(6) * t326 + t307) * t404 - t320 * t371 - t303 * t421) * t473 + (-(qJD(6) * t320 + t308) * t404 - t326 * t371 - t510 * t421) * t477) * MDP(31) + (t309 * t368 - t310 * t369 - t316 * t339 + t317 * t499) * MDP(26) + (t594 + t543 * (-t424 - t586)) * MDP(9) + t594 * MDP(2) + (qJDD(2) - t594 - 0.2e1 * t580) * MDP(4) + qJDD(1) * MDP(1) + t506 * MDP(3) + (qJD(2) * t415 - qJD(4) * t348 - qJDD(4) * t381 + t377 * t442 + t416 * t429 - t419 * t420 - t451 * t506) * MDP(17) + t546 * MDP(13) + (-t335 * t564 + t390 * t492) * MDP(18) + (t371 * t421 + t570) * MDP(29) + ((-t388 * t478 - t390 * t474) * t416 + (t577 + t336 * t478 + (-t388 * t474 + t390 * t478) * qJD(5)) * t420) * MDP(19) + (0.2e1 * t466 + t585 - t506) * MDP(5) + (-t509 * pkin(1) - g(1) * (-pkin(1) * t476 + t459) - g(2) * t545 + (t466 + t585) * qJ(2)) * MDP(6); (t509 - t594) * MDP(6) + (-qJD(1) * t449 + t517 - t594) * MDP(10) + (-qJD(1) * t413 + t546) * MDP(16) + (-qJD(1) * t415 + t496) * MDP(17) + (-t421 * t556 + t336 * t420 - t388 * t416 + (-t417 * t474 - t478 * t508) * t603) * MDP(23) + (-t421 * t364 + t335 * t420 - t390 * t416 + (-t417 * t478 + t474 * t508) * t603) * MDP(24) + (t420 * t310 - t416 * t339 - t423 * t570 + qJD(1) * t489 + ((t477 * t591 + t592) * t404 - t574) * t421) * MDP(30) + (t420 * t309 + t416 * t499 + t417 * t489 + t423 * t404 * qJD(1) + (-(t473 * t591 - t474 * t534 - t477 * t537) * t404 + t345) * t421) * MDP(31) + (MDP(4) - t518) * qJDD(1) + (-MDP(6) * qJ(2) - MDP(5) - t590) * t481; (qJD(1) * t598 + t595) * MDP(10) - t435 * MDP(16) + t437 * MDP(17) + (t494 - t573) * MDP(23) + (-t572 + t604) * MDP(24) + (t504 - t576) * MDP(30) + (t575 + t588) * MDP(31) - t481 * t518 + ((MDP(16) * t475 + MDP(8)) * t471 + (MDP(16) * t479 - MDP(17) * t475 + MDP(7)) * t470) * qJDD(1) + ((t415 + t525) * MDP(16) + (-t470 * t540 - t471 * t541 - t413) * MDP(17)) * qJD(4); -t413 ^ 2 * MDP(12) + t500 * MDP(13) + (-qJD(4) * t525 - t587) * MDP(14) + qJDD(4) * MDP(15) + (qJD(4) * t358 - t485 + t495) * MDP(16) + (t413 * t429 + t450 * t594 - t498 + t582) * MDP(17) + (t390 * t512 + t577) * MDP(18) + ((t335 - t609) * t478 + (-t390 * t603 - t336) * t474) * MDP(19) + (-t572 - t604) * MDP(20) + (t494 + t573) * MDP(21) + (-pkin(4) * t336 - t358 * t388 - t367 * t603 + (t596 * t603 + t491) * t474 - t608 * t478) * MDP(23) + (-pkin(4) * t335 - t358 * t390 + t474 * t608 + t491 * t478 + t550 * t603) * MDP(24) + (t309 * t423 - t499 * t548) * MDP(25) + (-t309 * t422 - t310 * t423 - t339 * t548 + t499 * t547) * MDP(26) + (t575 - t588) * MDP(27) + (t504 + t576) * MDP(28) + ((-t433 * t477 - t434 * t473) * t371 + t448 * t310 + t311 * t422 + (t473 * t501 - t477 * t502) * t404 + t507 * t339 + t547 * t331 - t485 * t456) * MDP(30) + (-(-t433 * t473 + t434 * t477) * t371 + t448 * t309 + t311 * t423 + (t473 * t502 + t477 * t501) * t404 - t507 * t499 + t548 * t331 + t485 * t455) * MDP(31) + (MDP(11) * t413 + MDP(12) * t415 + MDP(14) * qJD(4) - MDP(16) * t429 - MDP(22) * t603 - MDP(23) * t324 + MDP(24) * t325 - MDP(29) * t404 - MDP(30) * t305 + MDP(31) * t306) * t415; t390 * t388 * MDP(18) + (-t388 ^ 2 + t390 ^ 2) * MDP(19) + (t335 + t609) * MDP(20) + (-t513 + (-qJD(5) + t603) * t390) * MDP(21) + t374 * MDP(22) + (-g(1) * t409 - g(2) * t411 + t325 * t603 - t352 * t390 + (t511 + t582) * t474 + t503) * MDP(23) + (g(1) * t410 - g(2) * t412 + t324 * t603 + t352 * t388 + t478 * t582 - t487) * MDP(24) + (t309 + t607) * MDP(27) + (-t310 - t606) * MDP(28) + (-(-t314 * t473 - t578) * t404 - t306 * qJD(6) + (-t339 * t390 + t371 * t477 - t404 * t535) * pkin(5) + t601) * MDP(30) + ((-t315 * t404 - t303) * t473 + (t314 * t404 - t510) * t477 + (-t371 * t473 + t390 * t499 - t404 * t534) * pkin(5) + t602) * MDP(31) + t600; (t527 + t607) * MDP(27) + (-t514 - t606) * MDP(28) + (t306 * t404 + t601) * MDP(30) + (-t473 * t303 - t477 * t304 + t305 * t404 + t602) * MDP(31) + (-MDP(27) * t571 + MDP(28) * t499 - MDP(30) * t306 - MDP(31) * t579) * qJD(6) + t600;];
tau  = t1;
