% Calculate vector of inverse dynamics joint torques for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR11_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:28
% EndTime: 2019-12-31 21:35:40
% DurationCPUTime: 7.99s
% Computational Cost: add. (2964->525), mult. (6490->671), div. (0->0), fcn. (4342->8), ass. (0->217)
t464 = sin(qJ(3));
t468 = cos(qJ(3));
t544 = t468 * qJD(2);
t465 = sin(qJ(2));
t558 = qJD(1) * t465;
t408 = t464 * t558 - t544;
t532 = t468 * t558;
t555 = qJD(2) * t464;
t410 = t532 + t555;
t428 = -qJD(2) * pkin(2) + pkin(6) * t558;
t489 = qJ(4) * t410 - t428;
t602 = pkin(3) + pkin(4);
t336 = -t408 * t602 + t489;
t469 = cos(qJ(2));
t470 = cos(qJ(1));
t569 = t470 * t464;
t466 = sin(qJ(1));
t572 = t466 * t468;
t395 = t469 * t569 - t572;
t570 = t468 * t470;
t396 = t464 * t466 + t469 * t570;
t463 = sin(qJ(5));
t467 = cos(qJ(5));
t349 = t395 * t467 - t396 * t463;
t358 = t408 * t463 + t410 * t467;
t573 = t465 * t468;
t574 = t465 * t467;
t386 = t463 * t573 - t464 * t574;
t540 = qJD(1) * qJD(2);
t525 = t469 * t540;
t539 = qJDD(1) * t465;
t551 = qJD(3) * t465;
t620 = qJD(1) * t551 - qJDD(2);
t351 = -qJD(3) * t544 + (-t525 - t539) * t468 + t620 * t464;
t458 = t469 * qJDD(1);
t613 = -t465 * t540 + t458;
t406 = qJDD(3) - t613;
t516 = pkin(2) * t465 - pkin(7) * t469;
t419 = t516 * qJD(2);
t495 = pkin(2) * t469 + pkin(7) * t465 + pkin(1);
t370 = qJD(1) * t419 - qJDD(1) * t495;
t390 = pkin(6) * t613 + qJDD(2) * pkin(7);
t401 = t495 * qJD(1);
t557 = qJD(1) * t469;
t455 = pkin(6) * t557;
t429 = qJD(2) * pkin(7) + t455;
t550 = qJD(3) * t468;
t552 = qJD(3) * t464;
t518 = -t468 * t370 + t464 * t390 - t401 * t552 + t429 * t550;
t496 = qJDD(4) + t518;
t314 = pkin(8) * t351 - t406 * t602 + t496;
t399 = t406 * qJ(4);
t621 = qJD(3) - t557;
t432 = t621 * qJD(4);
t482 = t464 * t370 + t468 * t390 - t401 * t550 - t429 * t552;
t322 = t399 + t432 + t482;
t352 = t464 * (qJD(2) * (qJD(3) + t557) + t539) + t620 * t468;
t316 = pkin(8) * t352 + t322;
t523 = t467 * t314 - t463 * t316;
t575 = t464 * t469;
t393 = t466 * t575 + t570;
t571 = t468 * t469;
t394 = t466 * t571 - t569;
t611 = t393 * t467 - t394 * t463;
t625 = g(1) * t349 + g(2) * t611 - g(3) * t386 + t336 * t358 - t523;
t400 = -qJDD(5) + t406;
t498 = -t467 * t408 + t410 * t463;
t624 = MDP(22) * t358 * t498 + (t358 ^ 2 - t498 ^ 2) * MDP(23) - t400 * MDP(26);
t541 = -qJD(5) + t621;
t623 = t358 * t541;
t592 = g(3) * t469;
t513 = g(1) * t470 + g(2) * t466;
t609 = t465 * t513;
t622 = t592 - t609;
t610 = qJD(3) - qJD(5);
t616 = t541 * t498;
t615 = qJD(4) * t464 + t455;
t527 = t464 * t551;
t612 = t469 * t544 - t527;
t362 = -t468 * t401 - t464 * t429;
t542 = qJD(4) - t362;
t348 = pkin(3) * t408 - t489;
t598 = pkin(7) * t406;
t608 = -t348 * t621 + t598;
t589 = qJ(4) * t464;
t607 = -t468 * t602 - t589;
t350 = t395 * t463 + t396 * t467;
t412 = t463 * t464 + t467 * t468;
t387 = t412 * t465;
t499 = t393 * t463 + t394 * t467;
t543 = pkin(8) * t410 - t542;
t328 = -t602 * t621 - t543;
t546 = qJD(5) * t467;
t536 = t463 * t314 + t467 * t316 + t328 * t546;
t604 = -g(1) * t350 - g(2) * t499 - g(3) * t387 - t336 * t498 + t536;
t603 = t410 ^ 2;
t601 = pkin(7) - pkin(8);
t600 = pkin(3) * t406;
t599 = pkin(3) * t464;
t593 = g(3) * t465;
t591 = pkin(7) * qJD(3);
t590 = qJ(4) * t408;
t588 = qJ(4) * t468;
t363 = -t464 * t401 + t468 * t429;
t434 = t621 * qJ(4);
t343 = t434 + t363;
t587 = t343 * t621;
t586 = t351 * t464;
t585 = t363 * t621;
t583 = t408 * t410;
t582 = t408 * t621;
t581 = t410 * t621;
t580 = t410 * t468;
t416 = t516 * qJD(1);
t579 = t416 * t468;
t339 = pkin(8) * t408 + t363;
t334 = t339 + t434;
t578 = t463 * t334;
t577 = t463 * t468;
t576 = t464 * t465;
t484 = t412 * t469;
t568 = -qJD(1) * t484 + t412 * t610;
t537 = t467 * t575;
t547 = qJD(5) * t463;
t567 = qJD(1) * t537 + t463 * t550 + t464 * t546 - t467 * t552 - t468 * t547 - t557 * t577;
t538 = t602 * t464;
t488 = -t538 + t588;
t566 = t488 * t621 + t615;
t506 = -t588 + t599;
t565 = t506 * t621 - t615;
t397 = t464 * t416;
t564 = qJ(4) * t558 + t397;
t563 = t464 * t419 - t495 * t550;
t562 = (g(1) * t570 + g(2) * t572) * t465;
t446 = pkin(6) * t571;
t561 = -t464 * t495 + t446;
t461 = t465 ^ 2;
t560 = -t469 ^ 2 + t461;
t556 = qJD(2) * t410;
t554 = qJD(2) * t465;
t553 = qJD(2) * t469;
t548 = qJD(4) * t468;
t431 = t601 * t468;
t535 = -t467 * t351 + t463 * t352 + t408 * t546;
t534 = qJ(4) * t554 + t563;
t533 = -pkin(6) * t464 - pkin(3);
t531 = t621 * t555;
t530 = t621 * t544;
t528 = t621 * t552;
t522 = -t351 * t463 - t467 * t352;
t445 = pkin(6) * t575;
t520 = -t468 * t495 - t445;
t517 = t533 * t465;
t515 = -g(1) * t393 + g(2) * t395;
t514 = g(1) * t394 - g(2) * t396;
t512 = g(1) * t466 - g(2) * t470;
t453 = pkin(6) * t539;
t511 = qJDD(2) * pkin(2) - pkin(6) * t525 - t453;
t372 = -qJ(4) * t469 + t561;
t510 = -qJD(3) * t446 + t419 * t468 + t495 * t552;
t479 = -pkin(8) * t571 + (-pkin(4) + t533) * t465;
t509 = qJD(1) * t479 - t431 * t610 - t579;
t430 = t601 * t464;
t508 = -qJD(5) * t430 + (-t573 * pkin(6) + pkin(8) * t575) * qJD(1) + t564 + t601 * t552;
t507 = pkin(3) * t468 + t589;
t505 = qJD(3) * t428 - t598;
t504 = qJ(4) * t467 - t463 * t602;
t503 = qJ(4) * t463 + t467 * t602;
t319 = t463 * t328 + t467 * t334;
t342 = -pkin(3) * t621 + t542;
t502 = t342 * t468 - t343 * t464;
t460 = t469 * pkin(3);
t353 = pkin(4) * t469 + t445 + t460 + (-pkin(8) * t465 + t495) * t468;
t361 = pkin(8) * t576 + t372;
t501 = t353 * t467 - t361 * t463;
t500 = t353 * t463 + t361 * t467;
t497 = -t464 * t467 + t577;
t493 = pkin(2) + t507;
t492 = t591 * t621 + t592;
t487 = -0.2e1 * pkin(1) * t540 - pkin(6) * qJDD(2);
t486 = t406 * t464 + t550 * t621;
t485 = t406 * t468 - t528;
t476 = -qJ(4) * t351 + qJD(4) * t410 + t511;
t323 = pkin(3) * t352 - t476;
t483 = -t323 - t492;
t320 = -t410 * t547 + t535;
t473 = qJD(1) ^ 2;
t481 = pkin(1) * t473 + t513;
t480 = -t469 * t513 - t593;
t472 = qJD(2) ^ 2;
t478 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t472 + t512;
t477 = g(1) * t395 + g(2) * t393 + g(3) * t576 - t518;
t321 = qJD(5) * t358 + t522;
t475 = t348 * t410 + qJDD(4) - t477;
t474 = g(1) * t396 + g(2) * t394 + g(3) * t573 + t362 * t621 - t482;
t441 = qJ(4) * t573;
t404 = pkin(2) - t607;
t382 = -t441 + (pkin(6) + t599) * t465;
t373 = t460 - t520;
t371 = t441 + (-pkin(6) - t538) * t465;
t369 = pkin(3) * t410 + t590;
t368 = qJD(1) * t517 - t579;
t367 = -pkin(6) * t532 + t564;
t340 = -t410 * t602 - t590;
t337 = (qJD(3) * t507 - t548) * t465 + (pkin(6) + t506) * t553;
t335 = qJD(2) * t517 - t510;
t333 = -t351 + t582;
t332 = -qJD(4) * t469 + (-t465 * t544 - t469 * t552) * pkin(6) + t534;
t331 = t465 * t497 * t610 + qJD(2) * t484;
t330 = -qJD(2) * t537 + qJD(5) * t387 + t463 * t612 - t550 * t574;
t329 = (qJD(3) * t607 + t548) * t465 + (-pkin(6) + t488) * t553;
t326 = (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t573 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t464) * t469 + t534;
t325 = pkin(8) * t527 + qJD(2) * t479 - t510;
t324 = t496 - t600;
t318 = t328 * t467 - t578;
t317 = -t352 * t602 + t476;
t1 = [(-t351 * t573 + t410 * t612) * MDP(11) + (t465 * t487 + t469 * t478) * MDP(9) + (-t465 * t478 + t469 * t487) * MDP(10) + ((t351 + t530) * t469 + (t485 + t556) * t465) * MDP(13) + ((t352 - t531) * t469 + (-qJD(2) * t408 - t486) * t465) * MDP(14) + ((qJD(5) * t501 + t325 * t463 + t326 * t467) * t541 + t500 * t400 - (-t334 * t547 + t536) * t469 + t319 * t554 + t329 * t358 + t371 * t320 + t317 * t387 + t336 * t331 + g(1) * t611 - g(2) * t349) * MDP(28) + t513 * MDP(3) + t512 * MDP(2) + (-t335 * t621 + t337 * t408 + t352 * t382 - t373 * t406 + (t348 * t555 + t324) * t469 + (-qJD(2) * t342 + t323 * t464 + t348 * t550) * t465 + t514) * MDP(18) + (-t563 * t621 - t561 * t406 + (t428 * t544 + (t528 + t556) * pkin(6) + t482) * t469 + (-t428 * t552 - t363 * qJD(2) - t511 * t468 + (-t351 + t530) * pkin(6)) * t465 + t515) * MDP(17) + (t510 * t621 + t520 * t406 + ((pkin(6) * t408 + t428 * t464) * qJD(2) + t518) * t469 + (t428 * t550 + t362 * qJD(2) - t511 * t464 + (t352 + t531) * pkin(6)) * t465 + t514) * MDP(16) + (t332 * t621 - t337 * t410 + t351 * t382 + t372 * t406 + (-t348 * t544 - t322) * t469 + (qJD(2) * t343 - t323 * t468 + t348 * t552) * t465 - t515) * MDP(20) + (-t406 * t469 + t554 * t621) * MDP(15) + (t322 * t372 + t343 * t332 + t323 * t382 + t348 * t337 + t324 * t373 + t342 * t335 - g(1) * (-pkin(3) * t394 - qJ(4) * t393) - g(2) * (pkin(3) * t396 + qJ(4) * t395) + (-g(1) * pkin(6) - g(2) * t495) * t470 + (-g(2) * pkin(6) + g(1) * t495) * t466) * MDP(21) + (-t320 * t386 - t321 * t387 - t330 * t358 - t331 * t498) * MDP(23) + (qJDD(2) * t465 + t469 * t472) * MDP(6) + (qJDD(2) * t469 - t465 * t472) * MDP(7) + (t320 * t387 + t331 * t358) * MDP(22) + (-(t325 * t467 - t326 * t463) * t541 - t501 * t400 + t523 * t469 - t318 * t554 + t329 * t498 + t371 * t321 + t317 * t386 + t336 * t330 + g(1) * t499 - g(2) * t350 + (-t319 * t469 + t500 * t541) * qJD(5)) * MDP(27) + (-t321 * t469 + t330 * t541 + t386 * t400 + t498 * t554) * MDP(25) + (-t400 * t469 + t541 * t554) * MDP(26) + (t320 * t469 - t331 * t541 - t358 * t554 - t387 * t400) * MDP(24) + qJDD(1) * MDP(1) + (qJDD(1) * t461 + 0.2e1 * t465 * t525) * MDP(4) + (-t332 * t408 + t335 * t410 - t351 * t373 - t352 * t372 + t502 * t553 + (-t322 * t464 + t324 * t468 + (-t342 * t464 - t343 * t468) * qJD(3) + t512) * t465) * MDP(19) + 0.2e1 * (t458 * t465 - t540 * t560) * MDP(5) + ((-t408 * t468 - t410 * t464) * t553 + (t586 - t352 * t468 + (t408 * t464 - t580) * qJD(3)) * t465) * MDP(12); MDP(7) * t458 + MDP(6) * t539 + (-pkin(2) * t352 + t505 * t464 + (-t592 + t511 - (t416 + t591) * t621) * t468 + (-t428 * t575 - t362 * t465 + (-t408 * t469 - t576 * t621) * pkin(6)) * qJD(1) + t562) * MDP(16) + (-t320 * t412 + t321 * t497 - t358 * t567 - t498 * t568) * MDP(23) + (-t320 * t497 + t358 * t568) * MDP(22) + (t400 * t497 - t541 * t568) * MDP(24) + (t400 * t412 + t541 * t567) * MDP(25) + (-t352 * t493 + t368 * t621 + t565 * t408 - t464 * t608 + t483 * t468 + t562) * MDP(18) + (-(t430 * t467 - t431 * t463) * t400 + t404 * t321 - g(3) * t484 - (t463 * t508 - t467 * t509) * t541 + t566 * t498 + t567 * t336 + (t317 + t609) * t412) * MDP(27) + qJDD(2) * MDP(8) + ((t430 * t463 + t431 * t467) * t400 + t404 * t320 - (t463 * t509 + t467 * t508) * t541 + t566 * t358 + t568 * t336 + (-t317 + t622) * t497) * MDP(28) + (pkin(2) * t351 + t397 * t621 + t505 * t468 + (-t428 * t571 + t363 * t465 + (-t410 * t469 - t573 * t621) * pkin(6)) * qJD(1) + (-t511 + t492 - t609) * t464) * MDP(17) + (-t351 * t493 - t367 * t621 - t565 * t410 + t608 * t468 + (t483 + t609) * t464) * MDP(20) + ((-t410 * t465 - t571 * t621) * qJD(1) + t486) * MDP(13) + ((t408 * t465 + t575 * t621) * qJD(1) + t485) * MDP(14) + (-t342 * t368 - t343 * t367 + t565 * t348 + (qJD(3) * t502 + t322 * t468 + t324 * t464 + t480) * pkin(7) + (-t323 - t622) * t493) * MDP(21) + ((-t351 - t582) * t468 + (-t352 - t581) * t464) * MDP(12) + (t580 * t621 - t586) * MDP(11) + (t367 * t408 - t368 * t410 + (t322 + t621 * t342 + (qJD(3) * t410 - t352) * pkin(7)) * t468 + (t324 - t587 + (qJD(3) * t408 - t351) * pkin(7)) * t464 + t480) * MDP(19) + (t465 * t481 - t453 - t592) * MDP(9) + (t593 + (-pkin(6) * qJDD(1) + t481) * t469) * MDP(10) + (-MDP(15) * t621 + t342 * MDP(18) - t343 * MDP(20) + t358 * MDP(24) - MDP(25) * t498 - MDP(26) * t541 + t318 * MDP(27) - t319 * MDP(28)) * t558 + (-MDP(4) * t465 * t469 + MDP(5) * t560) * t473; MDP(11) * t583 + (-t408 ^ 2 + t603) * MDP(12) + t333 * MDP(13) + (-t352 + t581) * MDP(14) + t406 * MDP(15) + (-t410 * t428 + t477 + t585) * MDP(16) + (t408 * t428 + t474) * MDP(17) + (-t369 * t408 - t475 + t585 + 0.2e1 * t600) * MDP(18) + (pkin(3) * t351 - qJ(4) * t352 + (t343 - t363) * t410 + (t342 - t542) * t408) * MDP(19) + (-t348 * t408 + t369 * t410 + 0.2e1 * t399 + 0.2e1 * t432 - t474) * MDP(20) + (t322 * qJ(4) - t324 * pkin(3) - t348 * t369 - t342 * t363 - g(1) * (-pkin(3) * t395 + qJ(4) * t396) - g(2) * (-pkin(3) * t393 + qJ(4) * t394) - g(3) * (-pkin(3) * t576 + t441) + t542 * t343) * MDP(21) + (-t320 + t616) * MDP(24) + (t321 + t623) * MDP(25) + (t503 * t400 - t340 * t498 - (-t467 * t339 + t463 * t543) * t541 + (t504 * t541 + t319) * qJD(5) + t625) * MDP(27) + (t504 * t400 - t340 * t358 - (t463 * t339 + t467 * t543) * t541 + (-t503 * t541 - t578) * qJD(5) + t604) * MDP(28) - t624; (-t406 + t583) * MDP(18) + t333 * MDP(19) + (-t621 ^ 2 - t603) * MDP(20) + (t475 - t587 - t600) * MDP(21) + (-t467 * t400 - t410 * t498) * MDP(27) + (-t358 * t410 + t463 * t400) * MDP(28) - (MDP(27) * t463 + MDP(28) * t467) * t541 ^ 2; (t535 - t616) * MDP(24) + (-t522 - t623) * MDP(25) + (-t319 * t541 - t625) * MDP(27) + (-t318 * t541 - t604) * MDP(28) + ((-MDP(25) * t410 - MDP(27) * t334) * t467 + (-MDP(24) * t410 - MDP(25) * t408 - MDP(27) * t328 + MDP(28) * t334) * t463) * qJD(5) + t624;];
tau = t1;
