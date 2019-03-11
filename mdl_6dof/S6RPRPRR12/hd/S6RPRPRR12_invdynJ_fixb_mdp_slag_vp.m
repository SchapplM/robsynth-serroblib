% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:20:15
% EndTime: 2019-03-09 04:20:21
% DurationCPUTime: 6.36s
% Computational Cost: add. (2912->516), mult. (5429->648), div. (0->0), fcn. (3395->10), ass. (0->239)
t494 = -pkin(1) - pkin(7);
t620 = pkin(4) - t494;
t490 = cos(qJ(5));
t487 = sin(qJ(3));
t584 = qJD(1) * t487;
t551 = t490 * t584;
t486 = sin(qJ(5));
t578 = qJD(3) * t486;
t419 = -t551 + t578;
t489 = cos(qJ(6));
t576 = qJD(3) * t490;
t421 = t486 * t584 + t576;
t485 = sin(qJ(6));
t608 = t421 * t485;
t357 = t489 * t419 + t608;
t491 = cos(qJ(3));
t583 = qJD(1) * t491;
t456 = qJD(5) + t583;
t445 = qJD(6) + t456;
t642 = t357 * t445;
t518 = t419 * t485 - t489 * t421;
t641 = t445 * t518;
t467 = t491 * qJDD(1);
t565 = qJD(1) * qJD(3);
t549 = t487 * t565;
t633 = -t549 + t467;
t418 = -qJDD(5) - t633;
t410 = -qJDD(6) + t418;
t600 = t489 * t490;
t607 = t485 * t486;
t517 = -t600 + t607;
t640 = t517 * t410;
t424 = t485 * t490 + t486 * t489;
t626 = qJD(5) + qJD(6);
t596 = (-t583 - t626) * t424;
t548 = t491 * t565;
t564 = qJDD(1) * t487;
t639 = t548 + t564;
t591 = pkin(3) * t584 + qJD(1) * qJ(2);
t401 = -qJ(4) * t583 + t591;
t492 = cos(qJ(1));
t474 = g(2) * t492;
t488 = sin(qJ(1));
t475 = g(1) * t488;
t634 = t475 - t474;
t638 = -qJD(1) * t401 - t634;
t493 = -pkin(3) - pkin(8);
t446 = qJD(1) * t494 + qJD(2);
t429 = t491 * t446;
t627 = qJD(4) - t429;
t567 = pkin(4) * t583 + t627;
t371 = qJD(3) * t493 + t567;
t616 = qJ(4) * t491;
t525 = pkin(8) * t487 - t616;
t373 = qJD(1) * t525 + t591;
t339 = t371 * t486 + t373 * t490;
t333 = -pkin(9) * t419 + t339;
t570 = qJD(6) * t485;
t331 = t333 * t570;
t428 = t487 * t446;
t390 = -pkin(4) * t584 + t428;
t480 = qJD(3) * qJ(4);
t375 = t390 + t480;
t355 = pkin(5) * t419 + t375;
t484 = qJ(5) + qJ(6);
t468 = sin(t484);
t469 = cos(t484);
t598 = t491 * t492;
t382 = -t468 * t598 - t469 * t488;
t601 = t488 * t491;
t384 = -t468 * t601 + t469 * t492;
t622 = g(3) * t487;
t637 = g(1) * t384 - g(2) * t382 + t355 * t357 + t468 * t622 + t331;
t381 = t468 * t488 - t469 * t598;
t383 = -t468 * t492 - t469 * t601;
t573 = qJD(5) * t486;
t353 = -qJD(3) * t573 + qJD(5) * t551 + t490 * qJDD(3) + t486 * t639;
t477 = qJDD(1) * qJ(2);
t479 = qJD(1) * qJD(2);
t521 = pkin(3) * t639 + qJ(4) * t549 + t477 + t479;
t542 = qJD(3) * pkin(8) - qJD(4);
t566 = qJ(4) * qJDD(1);
t340 = pkin(8) * t564 + (qJD(1) * t542 - t566) * t491 + t521;
t443 = qJDD(1) * t494 + qJDD(2);
t577 = qJD(3) * t487;
t560 = t446 * t577 + qJDD(4);
t522 = -t443 * t491 + t560;
t344 = pkin(4) * t633 + t493 * qJDD(3) + t522;
t539 = -t340 * t486 + t490 * t344;
t501 = -t339 * qJD(5) + t539;
t321 = -pkin(5) * t418 - pkin(9) * t353 + t501;
t354 = t421 * qJD(5) + qJDD(3) * t486 - t490 * t639;
t572 = qJD(5) * t490;
t558 = -t490 * t340 - t486 * t344 - t371 * t572;
t322 = -pkin(9) * t354 - t373 * t573 - t558;
t540 = t489 * t321 - t485 * t322;
t636 = -g(1) * t383 + g(2) * t381 + t355 * t518 - t469 * t622 + t540;
t635 = (-t357 ^ 2 + t518 ^ 2) * MDP(26) - t410 * MDP(29) - t357 * MDP(25) * t518;
t389 = t424 * t487;
t472 = t487 * pkin(3);
t411 = qJ(2) + t472 + t525;
t434 = t620 * t491;
t412 = t486 * t434;
t593 = t490 * t411 + t412;
t397 = t490 * t418;
t632 = -t456 * t573 - t397;
t631 = qJD(6) * t490 + t572;
t580 = qJD(3) * t419;
t630 = t580 + t397;
t629 = -t616 + t472;
t559 = 0.2e1 * t479;
t628 = 0.2e1 * t477 + t559;
t625 = t375 * t456 - t493 * t418;
t427 = t487 * t443;
t476 = qJDD(3) * qJ(4);
t478 = qJD(3) * qJD(4);
t575 = qJD(3) * t491;
t361 = -t446 * t575 - t427 - t476 - t478;
t615 = qJDD(3) * pkin(3);
t362 = t522 - t615;
t543 = qJD(3) * pkin(3) - qJD(4);
t398 = -t429 - t543;
t402 = -t428 - t480;
t499 = qJD(3) * (t398 * t487 - t402 * t491) - t361 * t487 - t362 * t491;
t430 = qJ(2) + t629;
t561 = qJDD(3) * t494;
t624 = (qJD(1) * t430 + t401) * qJD(3) + t561;
t538 = t353 * t485 + t489 * t354;
t326 = -qJD(6) * t518 + t538;
t621 = g(3) * t491;
t619 = pkin(9) - t493;
t618 = pkin(1) * qJDD(1);
t496 = qJD(1) ^ 2;
t617 = qJ(2) * t496;
t338 = t490 * t371 - t373 * t486;
t332 = -pkin(9) * t421 + t338;
t330 = pkin(5) * t456 + t332;
t614 = t330 * t489;
t613 = t333 * t489;
t612 = t353 * t490;
t611 = t410 * t424;
t610 = t419 * t456;
t609 = t421 * t456;
t606 = t486 * t418;
t605 = t486 * t487;
t604 = t486 * t491;
t603 = t487 * t490;
t602 = t487 * t496;
t599 = t490 * t491;
t505 = -t485 * t573 - t486 * t570;
t555 = t490 * t583;
t595 = t489 * t555 - t583 * t607 + t600 * t626 + t505;
t426 = pkin(3) * t583 + qJ(4) * t584;
t387 = pkin(8) * t583 + t426;
t594 = t490 * t387 + t486 * t390;
t556 = -pkin(5) * t490 - pkin(4);
t592 = pkin(5) * t572 - t556 * t583 + t627;
t590 = t492 * pkin(1) + t488 * qJ(2);
t482 = t487 ^ 2;
t483 = t491 ^ 2;
t588 = t482 - t483;
t495 = qJD(3) ^ 2;
t587 = t495 + t496;
t582 = qJD(3) * t357;
t581 = qJD(3) * t518;
t579 = qJD(3) * t421;
t574 = qJD(5) * t373;
t571 = qJD(5) * t493;
t569 = qJD(6) * t489;
t563 = qJDD(3) * t487;
t562 = qJDD(3) * t491;
t557 = t489 * t353 - t485 * t354 - t419 * t569;
t554 = t486 * t575;
t553 = t490 * t575;
t550 = pkin(3) * t575 + qJ(4) * t577 + qJD(2);
t433 = t619 * t490;
t547 = -pkin(9) * t487 - t411;
t546 = -t443 - t474;
t457 = g(1) * t601;
t545 = -t457 + t622;
t544 = qJD(1) * t426 - g(3);
t372 = t491 * t542 + t550;
t414 = t620 * t577;
t537 = -t372 * t486 - t490 * t414;
t535 = qJD(6) * t330 + t322;
t534 = (-t482 - t483) * qJDD(1);
t533 = qJDD(2) - t618;
t532 = qJD(5) * t491 + qJD(1);
t415 = t620 * t575;
t531 = g(1) * t492 + g(2) * t488;
t529 = t445 * t596 + t640;
t377 = t490 * t390;
t431 = t619 * t486;
t514 = -pkin(5) * t487 - pkin(9) * t604;
t528 = qJD(1) * t514 - qJD(6) * t431 - t387 * t486 - t619 * t573 + t377;
t527 = pkin(9) * t555 + t433 * t626 + t594;
t526 = pkin(3) * t491 + qJ(4) * t487;
t524 = -t617 + t474;
t324 = t330 * t485 + t613;
t516 = t579 - t606;
t515 = t456 * t486;
t513 = -t456 * t572 + t606;
t511 = t445 * t517;
t510 = -t445 * t595 + t611;
t509 = 0.2e1 * qJ(2) * t565 + t561;
t508 = t490 * t372 - t411 * t573 - t486 * t414 + t434 * t572;
t325 = -t421 * t570 + t557;
t507 = qJD(3) * t424;
t506 = -t494 * t495 - t531;
t504 = -t487 * t573 + t553;
t503 = -t487 * t634 - t621;
t345 = -pkin(4) * t639 - t361;
t502 = t345 + t503;
t498 = t506 + t628;
t347 = (-qJD(1) * qJD(4) - t566) * t491 + t521;
t385 = -qJD(4) * t491 + t550;
t497 = -qJD(1) * t385 - qJDD(1) * t430 - t347 - t506;
t471 = t492 * qJ(2);
t462 = t487 * t494;
t458 = pkin(5) * t486 + qJ(4);
t432 = -pkin(4) * t487 + t462;
t413 = t490 * t434;
t407 = -t486 * t601 + t490 * t492;
t406 = -t486 * t492 - t488 * t599;
t405 = -t486 * t598 - t488 * t490;
t404 = t486 * t488 - t490 * t598;
t392 = t487 * t556 + t462;
t388 = t485 * t605 - t487 * t600;
t380 = t401 * t583;
t365 = -pkin(5) * t504 - t415;
t351 = pkin(9) * t603 + t593;
t346 = pkin(5) * t491 + t486 * t547 + t413;
t335 = t389 * t626 + t485 * t554 - t489 * t553;
t334 = -t487 * t517 * t626 + t491 * t507;
t329 = pkin(5) * t354 + t345;
t328 = pkin(9) * t504 + t508;
t327 = t514 * qJD(3) + (t490 * t547 - t412) * qJD(5) + t537;
t323 = -t333 * t485 + t614;
t1 = [(t487 * t624 + t497 * t491) * MDP(16) + (t497 * t487 - t491 * t624) * MDP(15) + (-t531 + t628) * MDP(5) + (t347 * t430 + t401 * t385 - g(1) * (-qJ(4) * t598 + t472 * t492 + t471) - g(2) * (pkin(7) * t492 + t590) + (-g(1) * t494 - g(2) * t629) * t488 + t499 * t494) * MDP(17) + ((t456 * t576 - t354) * t491 + (t580 + t632) * t487) * MDP(21) + (t537 * t456 - (-t411 * t486 + t413) * t418 + t539 * t491 - t415 * t419 + t432 * t354 - t345 * t603 - g(1) * t405 - g(2) * t407 + (-t338 * t487 - t375 * t599) * qJD(3) + (-t339 * t491 + t375 * t605 - t456 * t593) * qJD(5)) * MDP(23) + (-t508 * t456 + t593 * t418 - t415 * t421 + t432 * t353 - g(1) * t404 - g(2) * t406 + ((qJD(3) * t375 + t574) * t486 + t558) * t491 + (qJD(3) * t339 + t345 * t486 + t375 * t572) * t487) * MDP(24) + t531 * MDP(3) + qJDD(1) * MDP(1) + ((t327 * t489 - t328 * t485) * t445 - (t346 * t489 - t351 * t485) * t410 + t540 * t491 - t323 * t577 + t365 * t357 + t392 * t326 + t329 * t388 + t355 * t335 - g(1) * t382 - g(2) * t384 + ((-t346 * t485 - t351 * t489) * t445 - t324 * t491) * qJD(6)) * MDP(30) + (-t418 * t491 - t456 * t577) * MDP(22) + (-t410 * t491 - t445 * t577) * MDP(29) + (-t326 * t491 - t335 * t445 + t357 * t577 + t388 * t410) * MDP(28) + (-t487 * t495 + t562) * MDP(9) + (-t491 * t495 - t563) * MDP(10) + (qJDD(1) * t483 - 0.2e1 * t487 * t548) * MDP(7) + ((t456 * t578 + t353) * t491 + (-t513 - t579) * t487) * MDP(20) + 0.2e1 * (-t467 * t487 + t565 * t588) * MDP(8) + (t487 * t498 + t491 * t509) * MDP(12) + (-t487 * t509 + t491 * t498) * MDP(13) + (t325 * t491 + t334 * t445 - t389 * t410 + t518 * t577) * MDP(27) + (t324 * t577 - g(1) * t381 - g(2) * t383 + t392 * t325 + t329 * t389 + t331 * t491 + t355 * t334 - t365 * t518 + (-(-qJD(6) * t351 + t327) * t445 + t346 * t410 - t321 * t491) * t485 + (-(qJD(6) * t346 + t328) * t445 + t351 * t410 - t535 * t491) * t489) * MDP(31) + (t325 * t389 - t334 * t518) * MDP(25) + (-t325 * t388 - t326 * t389 - t334 * t357 + t335 * t518) * MDP(26) + (t494 * t534 - t499 + t634) * MDP(14) + (qJDD(2) - t634 - 0.2e1 * t618) * MDP(4) + t634 * MDP(2) + (-t533 * pkin(1) - g(1) * (-pkin(1) * t488 + t471) - g(2) * t590 + (t559 + t477) * qJ(2)) * MDP(6) + ((-t419 * t486 + t421 * t490) * t575 + (t612 - t354 * t486 + (-t419 * t490 - t421 * t486) * qJD(5)) * t487) * MDP(19) + (t353 * t605 + (t487 * t572 + t554) * t421) * MDP(18); qJDD(1) * MDP(4) - t496 * MDP(5) + (t533 - t634 - t617) * MDP(6) + MDP(14) * t534 + (t499 + t638) * MDP(17) + (t354 * t487 + t630 * t491 + (t486 * t532 + t487 * t576) * t456) * MDP(23) + (t353 * t487 + t516 * t491 + (-t486 * t577 + t490 * t532) * t456) * MDP(24) + (t424 * t445 * qJD(1) + (-qJD(3) * t511 + t326) * t487 + ((t485 * t631 + t486 * t569 + t489 * t573) * t445 - t640 + t582) * t491) * MDP(30) + (-qJD(1) * t511 + (-t445 * t507 + t325) * t487 + (-(-t489 * t631 - t505) * t445 - t611 - t581) * t491) * MDP(31) + (MDP(12) - MDP(15)) * (-t487 * t587 + t562) + (-MDP(13) + MDP(16)) * (t491 * t587 + t563); (qJ(4) * t354 - t377 * t456 + t567 * t419 + t625 * t490 + ((t387 - t571) * t456 + t502) * t486) * MDP(23) + (qJ(4) * t353 + t594 * t456 + t567 * t421 - t625 * t486 + (-t456 * t571 + t502) * t490) * MDP(24) + ((t421 * t487 - t456 * t604) * qJD(1) + t632) * MDP(20) + (-t362 * pkin(3) + g(3) * t629 - t361 * qJ(4) - t398 * t428 - t401 * t426 - t402 * t627 - t526 * t634) * MDP(17) + t529 * MDP(27) + (-t325 * t517 - t518 * t596) * MDP(25) + (-t325 * t424 + t326 * t517 - t357 * t596 + t518 * t595) * MDP(26) + ((-t419 * t487 - t456 * t599) * qJD(1) + t513) * MDP(21) + ((t443 + t524) * t491 + t545) * MDP(12) + qJDD(3) * MDP(11) - MDP(10) * t564 + ((-t354 - t609) * t490 + (-t353 + t610) * t486) * MDP(19) + (-t526 * qJDD(1) + ((-t402 - t480) * t491 + (t398 + t543) * t487) * qJD(1)) * MDP(14) + t510 * MDP(28) + (t487 * t638 + t544 * t491 + t427 + 0.2e1 * t476 + 0.2e1 * t478) * MDP(16) + (t621 - t427 + (-t524 + t475) * t487) * MDP(13) - t588 * MDP(8) * t496 + ((-t431 * t489 - t433 * t485) * t410 + t458 * t325 - t329 * t517 + (t485 * t528 + t489 * t527) * t445 - t592 * t518 + t596 * t355 + t503 * t469) * MDP(31) + (-(t431 * t485 - t433 * t489) * t410 + t458 * t326 + t329 * t424 + (t485 * t527 - t489 * t528) * t445 + t592 * t357 + t595 * t355 + t503 * t468) * MDP(30) + (-t421 * t515 + t612) * MDP(18) + (t487 * t544 + t491 * t546 + qJDD(4) + t380 + t457 - 0.2e1 * t615) * MDP(15) + t491 * MDP(7) * t602 + MDP(9) * t467 + (t456 * MDP(22) + t338 * MDP(23) - t339 * MDP(24) - MDP(27) * t518 - t357 * MDP(28) + t445 * MDP(29) + t323 * MDP(30) - t324 * MDP(31)) * t584; qJDD(3) * MDP(15) + (-t483 * t496 - t495) * MDP(16) + (qJD(3) * t402 + t380 - t545 + t560 - t615) * MDP(17) - t630 * MDP(23) - t516 * MDP(24) + (t529 - t582) * MDP(30) + (t510 + t581) * MDP(31) + (qJDD(1) * MDP(14) - MDP(15) * t602 + MDP(17) * t546) * t491 + (-t456 * t490 * MDP(24) - MDP(23) * t515) * t456; t421 * t419 * MDP(18) + (-t419 ^ 2 + t421 ^ 2) * MDP(19) + (t353 + t610) * MDP(20) + (-t354 + t609) * MDP(21) - t418 * MDP(22) + (-g(1) * t406 + g(2) * t404 - g(3) * t603 + t339 * t456 - t375 * t421 + t501) * MDP(23) + (g(1) * t407 - g(2) * t405 + t338 * t456 + t375 * t419 + (t574 + t622) * t486 + t558) * MDP(24) + (t325 + t642) * MDP(27) + (-t326 - t641) * MDP(28) + (-(-t332 * t485 - t613) * t445 - t324 * qJD(6) + (-t357 * t421 - t410 * t489 - t445 * t570) * pkin(5) + t636) * MDP(30) + ((-t333 * t445 - t321) * t485 + (t332 * t445 - t535) * t489 + (t410 * t485 + t421 * t518 - t445 * t569) * pkin(5) + t637) * MDP(31) + t635; (t557 + t642) * MDP(27) + (-t538 - t641) * MDP(28) + (t324 * t445 + t636) * MDP(30) + (-t485 * t321 - t489 * t322 + t323 * t445 + t637) * MDP(31) + (-MDP(27) * t608 + MDP(28) * t518 - MDP(30) * t324 - MDP(31) * t614) * qJD(6) + t635;];
tau  = t1;
