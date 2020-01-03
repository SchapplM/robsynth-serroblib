% Calculate vector of inverse dynamics joint torques for
% S5RRRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP9_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:18
% EndTime: 2019-12-31 22:06:30
% DurationCPUTime: 7.26s
% Computational Cost: add. (5138->535), mult. (11175->679), div. (0->0), fcn. (7546->10), ass. (0->224)
t507 = sin(qJ(3));
t508 = sin(qJ(2));
t593 = qJD(1) * t508;
t567 = t507 * t593;
t510 = cos(qJ(3));
t588 = qJD(2) * t510;
t443 = -t567 + t588;
t590 = qJD(2) * t507;
t444 = t510 * t593 + t590;
t506 = sin(qJ(4));
t637 = cos(qJ(4));
t387 = -t443 * t637 + t444 * t506;
t537 = t443 * t506 + t444 * t637;
t665 = t387 * t537;
t511 = cos(qJ(2));
t610 = t510 * t511;
t545 = pkin(3) * t508 - pkin(8) * t610;
t638 = pkin(8) + pkin(7);
t573 = qJD(3) * t638;
t552 = pkin(2) * t508 - pkin(7) * t511;
t447 = t552 * qJD(1);
t598 = pkin(6) * t567 + t447 * t510;
t664 = qJD(1) * t545 + t510 * t573 + t598;
t432 = t507 * t447;
t615 = t508 * t510;
t616 = t507 * t511;
t663 = -t432 - (-pkin(6) * t615 - pkin(8) * t616) * qJD(1) - t507 * t573;
t580 = qJD(1) * qJD(2);
t563 = t511 * t580;
t577 = qJDD(1) * t508;
t656 = qJD(2) * qJD(3) + t563 + t577;
t557 = t507 * qJDD(2) + t510 * t656;
t579 = qJD(1) * qJD(3);
t562 = t508 * t579;
t527 = -t507 * t562 + t557;
t558 = t507 * t656 + t510 * t562;
t530 = qJDD(2) * t510 - t558;
t565 = t637 * qJD(4);
t584 = qJD(4) * t506;
t345 = -t443 * t565 + t444 * t584 - t506 * t530 - t527 * t637;
t592 = qJD(1) * t511;
t485 = -qJD(3) + t592;
t469 = -qJD(4) + t485;
t336 = -t387 * t469 - t345;
t346 = qJD(4) * t537 + t506 * t527 - t530 * t637;
t499 = t511 * qJDD(1);
t648 = -t508 * t580 + t499;
t440 = qJDD(3) - t648;
t436 = qJDD(4) + t440;
t639 = t537 ^ 2;
t662 = t436 * MDP(22) + (-t469 * t537 - t346) * MDP(21) + MDP(18) * t665 + (-t387 ^ 2 + t639) * MDP(19) + t336 * MDP(20);
t459 = -qJD(2) * pkin(2) + pkin(6) * t593;
t405 = -pkin(3) * t443 + t459;
t351 = pkin(4) * t387 - qJ(5) * t537 + t405;
t661 = t351 * t387;
t660 = t387 * t405;
t618 = t506 * t507;
t535 = t510 * t637 - t618;
t644 = qJD(3) + qJD(4);
t647 = qJD(3) * t637 + t565;
t602 = -t510 * t647 + t535 * t592 + t618 * t644;
t446 = t506 * t510 + t507 * t637;
t400 = t644 * t446;
t601 = -t446 * t592 + t400;
t585 = qJD(3) * t510;
t568 = t508 * t585;
t587 = qJD(2) * t511;
t571 = t507 * t587;
t659 = t568 + t571;
t461 = t638 * t507;
t462 = t638 * t510;
t404 = -t461 * t506 + t462 * t637;
t505 = qJ(3) + qJ(4);
t500 = sin(t505);
t509 = sin(qJ(1));
t512 = cos(qJ(1));
t549 = g(1) * t512 + g(2) * t509;
t629 = g(3) * t511;
t526 = t508 * t549 - t629;
t658 = t404 * t436 + t500 * t526;
t359 = pkin(4) * t537 + qJ(5) * t387;
t536 = -t461 * t637 - t462 * t506;
t652 = qJD(4) * t536 - t506 * t664 + t637 * t663;
t651 = qJD(4) * t404 + t506 * t663 + t637 * t664;
t455 = -pkin(2) * t511 - pkin(7) * t508 - pkin(1);
t442 = t510 * t455;
t635 = pkin(6) * t507;
t395 = -pkin(8) * t615 + t442 + (-pkin(3) - t635) * t511;
t487 = pkin(6) * t610;
t596 = t455 * t507 + t487;
t617 = t507 * t508;
t402 = -pkin(8) * t617 + t596;
t650 = t395 * t506 + t402 * t637;
t422 = t436 * qJ(5);
t454 = t469 * qJD(5);
t649 = t422 - t454;
t496 = pkin(6) * t592;
t586 = qJD(3) * t507;
t636 = pkin(3) * t507;
t554 = pkin(3) * t586 - t592 * t636 - t496;
t494 = pkin(6) * t577;
t628 = qJDD(2) * pkin(2);
t424 = pkin(6) * t563 + t494 - t628;
t646 = pkin(7) * qJD(3) * t485 - t424;
t613 = t509 * t511;
t428 = t507 * t613 + t510 * t512;
t609 = t511 * t512;
t430 = -t507 * t609 + t509 * t510;
t645 = -g(1) * t430 + g(2) * t428;
t425 = t436 * pkin(4);
t643 = t425 - qJDD(5);
t501 = cos(t505);
t620 = t501 * t512;
t413 = t500 * t613 + t620;
t608 = t512 * t500;
t614 = t509 * t501;
t415 = t511 * t608 - t614;
t450 = t552 * qJD(2);
t401 = qJD(1) * t450 + qJDD(1) * t455;
t393 = t510 * t401;
t423 = pkin(6) * t648 + qJDD(2) * pkin(7);
t437 = t455 * qJD(1);
t460 = qJD(2) * pkin(7) + t496;
t611 = t510 * t460;
t339 = -t507 * t423 + t393 - t557 * pkin(8) + t440 * pkin(3) + (-t611 + (pkin(8) * t593 - t437) * t507) * qJD(3);
t533 = t401 * t507 + t423 * t510 + t437 * t585 - t460 * t586;
t342 = pkin(8) * t530 + t533;
t396 = t437 * t510 - t460 * t507;
t375 = -pkin(8) * t444 + t396;
t371 = -pkin(3) * t485 + t375;
t397 = t437 * t507 + t611;
t376 = pkin(8) * t443 + t397;
t559 = -t339 * t637 + t342 * t506 + t371 * t584 + t376 * t565;
t622 = t500 * t508;
t523 = g(1) * t415 + g(2) * t413 + g(3) * t622 - t559;
t517 = t351 * t537 - t523 - t643;
t642 = -t405 * t537 + t523;
t589 = qJD(2) * t508;
t599 = t450 * t510 + t589 * t635;
t360 = t545 * qJD(2) + (-t487 + (pkin(8) * t508 - t455) * t507) * qJD(3) + t599;
t529 = -t508 * t588 - t511 * t586;
t600 = t450 * t507 + t455 * t585;
t362 = t529 * pkin(6) - pkin(8) * t659 + t600;
t641 = -qJD(4) * t650 + t360 * t637 - t362 * t506;
t630 = g(3) * t508;
t572 = t637 * t376;
t348 = t371 * t506 + t572;
t627 = t348 * t469;
t625 = t443 * t485;
t624 = t444 * t510;
t623 = t444 * t511;
t621 = t501 * t508;
t619 = t506 * t376;
t612 = t510 * t440;
t607 = pkin(4) * t601 + qJ(5) * t602 - qJD(5) * t446 + t554;
t606 = -qJ(5) * t593 + t652;
t605 = pkin(4) * t593 + t651;
t350 = t375 * t637 - t619;
t597 = pkin(3) * t565 + qJD(5) - t350;
t486 = pkin(3) * t617;
t451 = pkin(6) * t508 + t486;
t503 = t508 ^ 2;
t595 = -t511 ^ 2 + t503;
t591 = qJD(2) * t443;
t583 = t444 * qJD(2);
t582 = t459 * qJD(3);
t347 = t371 * t637 - t619;
t581 = qJD(5) - t347;
t406 = pkin(3) * t659 + pkin(6) * t587;
t493 = pkin(3) * t510 + pkin(2);
t574 = pkin(6) + t636;
t570 = t485 * t585;
t569 = t508 * t586;
t561 = -qJD(3) * t437 - t423;
t560 = t339 * t506 + t342 * t637 + t371 * t565 - t376 * t584;
t556 = t637 * t587;
t555 = -pkin(4) * t622 + qJ(5) * t621;
t349 = t375 * t506 + t572;
t553 = pkin(3) * t584 - t349;
t551 = -g(1) * t413 + g(2) * t415;
t414 = t501 * t613 - t608;
t416 = t500 * t509 + t501 * t609;
t550 = g(1) * t414 - g(2) * t416;
t548 = g(1) * t509 - g(2) * t512;
t547 = t460 * t585 - t393;
t546 = -pkin(7) * t440 + t582;
t544 = t493 * t511 + t508 * t638 + pkin(1);
t543 = pkin(4) * t501 + qJ(5) * t500 + t493;
t542 = -0.2e1 * pkin(1) * t580 - pkin(6) * qJDD(2);
t540 = t395 * t637 - t402 * t506;
t534 = t440 * t507 - t570;
t532 = t360 * t506 + t362 * t637 + t395 * t565 - t402 * t584;
t515 = qJD(1) ^ 2;
t531 = pkin(1) * t515 + t549;
t528 = t536 * t436 - t501 * t629 + (g(1) * t620 + g(2) * t614) * t508;
t514 = qJD(2) ^ 2;
t525 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t514 + t548;
t524 = g(1) * t416 + g(2) * t414 + g(3) * t621 - t560;
t522 = t527 * t508;
t521 = -g(1) * (-pkin(4) * t415 + qJ(5) * t416) - g(2) * (-pkin(4) * t413 + qJ(5) * t414);
t518 = -t347 * t469 + t524;
t370 = -pkin(3) * t530 + t424;
t492 = -pkin(3) * t637 - pkin(4);
t488 = pkin(3) * t506 + qJ(5);
t431 = t507 * t509 + t510 * t609;
t429 = t507 * t512 - t509 * t610;
t419 = t535 * t508;
t418 = t446 * t508;
t385 = -pkin(4) * t535 - qJ(5) * t446 - t493;
t372 = pkin(4) * t418 - qJ(5) * t419 + t451;
t364 = t507 * t556 - t506 * t569 - t584 * t617 + (t506 * t587 + t508 * t647) * t510;
t363 = t400 * t508 + t506 * t571 - t510 * t556;
t358 = pkin(4) * t511 - t540;
t357 = -qJ(5) * t511 + t650;
t353 = pkin(3) * t444 + t359;
t344 = -qJ(5) * t469 + t348;
t343 = pkin(4) * t469 + t581;
t335 = pkin(4) * t364 + qJ(5) * t363 - qJD(5) * t419 + t406;
t334 = -pkin(4) * t589 - t641;
t333 = qJ(5) * t589 - qJD(5) * t511 + t532;
t332 = pkin(4) * t346 + qJ(5) * t345 - qJD(5) * t537 + t370;
t331 = t559 - t643;
t330 = t560 + t649;
t1 = [(t345 * t511 + t363 * t469 + t419 * t436 + t537 * t589) * MDP(20) + (-t330 * t511 - t332 * t419 - t333 * t469 - t335 * t537 + t344 * t589 + t345 * t372 + t351 * t363 + t357 * t436 - t551) * MDP(27) + (-t345 * t419 - t363 * t537) * MDP(18) + (t345 * t418 - t346 * t419 + t363 * t387 - t364 * t537) * MDP(19) + (-t330 * t418 + t331 * t419 - t333 * t387 + t334 * t537 - t343 * t363 - t344 * t364 - t345 * t358 - t346 * t357 + t508 * t548) * MDP(26) + (-t451 * t345 - t348 * t589 - t405 * t363 + t370 * t419 + t406 * t537 - t436 * t650 + t469 * t532 + t511 * t560 + t551) * MDP(24) + ((t443 * t510 - t444 * t507) * t587 + (t510 * t530 - t557 * t507 + (-t624 + (-t443 + t567) * t507) * qJD(3)) * t508) * MDP(12) + (t451 * t346 + t347 * t589 + t405 * t364 + t370 * t418 + t406 * t387 + t540 * t436 - t469 * t641 + t559 * t511 + t550) * MDP(23) + ((-t485 * t588 - t557) * t511 + (t583 + t612 + (t485 + t592) * t586) * t508) * MDP(13) + (-(-t455 * t586 + t599) * t485 + t442 * t440 - g(1) * t429 - g(2) * t431 + ((t570 - t591) * pkin(6) + (-pkin(6) * t440 + qJD(2) * t459 - t561) * t507 + t547) * t511 + (-pkin(6) * t530 + t396 * qJD(2) + t424 * t507 + t510 * t582) * t508) * MDP(16) + (t600 * t485 - t596 * t440 - g(1) * t428 - g(2) * t430 + (t459 * t588 + t533) * t511 + (-t397 * qJD(2) + t424 * t510 - t507 * t582) * t508 + (t485 * t529 + t511 * t583 + t522) * pkin(6)) * MDP(17) + 0.2e1 * (t499 * t508 - t580 * t595) * MDP(5) + ((t485 * t590 - t530) * t511 + (-t534 + t591) * t508) * MDP(14) + (-t440 * t511 - t485 * t589) * MDP(15) + (-t436 * t511 - t469 * t589) * MDP(22) + (t346 * t511 + t364 * t469 - t387 * t589 - t418 * t436) * MDP(21) + (t331 * t511 + t332 * t418 + t334 * t469 + t335 * t387 - t343 * t589 + t346 * t372 + t351 * t364 - t358 * t436 + t550) * MDP(25) + (t510 * t522 + (t510 * t587 - t569) * t444) * MDP(11) + (t330 * t357 + t344 * t333 + t332 * t372 + t351 * t335 + t331 * t358 + t343 * t334 - g(1) * (-pkin(4) * t414 - qJ(5) * t413) - g(2) * (pkin(4) * t416 + qJ(5) * t415) + (-g(1) * t574 - g(2) * t544) * t512 + (g(1) * t544 - g(2) * t574) * t509) * MDP(28) + (qJDD(1) * t503 + 0.2e1 * t508 * t563) * MDP(4) + (qJDD(2) * t508 + t511 * t514) * MDP(6) + (qJDD(2) * t511 - t508 * t514) * MDP(7) + qJDD(1) * MDP(1) + (t508 * t542 + t511 * t525) * MDP(9) + (-t508 * t525 + t511 * t542) * MDP(10) + t548 * MDP(2) + t549 * MDP(3); MDP(7) * t499 + MDP(6) * t577 + (t330 * t404 - t331 * t536 + t332 * t385 + t607 * t351 + t606 * t344 + t605 * t343 + (-g(3) * t543 - t549 * t638) * t511 + (-g(3) * t638 + t543 * t549) * t508) * MDP(28) + (-pkin(2) * t558 + t598 * t485 + t546 * t507 + (-t396 * t508 + (pkin(6) * t443 - t459 * t507) * t511) * qJD(1) + (t526 + t628 + t646) * t510) * MDP(16) + (t508 * t531 - t494 - t629) * MDP(9) + (t493 * t345 + t370 * t446 - t602 * t405 + t469 * t652 + t537 * t554 - t658) * MDP(24) + (-pkin(2) * t557 - t432 * t485 + t546 * t510 + (-t459 * t610 + t397 * t508 + (t485 * t615 - t623) * pkin(6)) * qJD(1) + (t629 + (pkin(2) * t579 - t549) * t508 - t646) * t507) * MDP(17) + (t330 * t535 + t331 * t446 - t343 * t602 - t344 * t601 + t345 * t536 - t346 * t404 - t387 * t606 - t511 * t549 + t537 * t605 - t630) * MDP(26) + (t630 + (-pkin(6) * qJDD(1) + t531) * t511) * MDP(10) + (-t332 * t446 + t345 * t385 + t351 * t602 - t469 * t606 - t537 * t607 + t658) * MDP(27) + ((t557 - t625) * t510 + (-t444 * qJD(3) + (-t568 + t623) * qJD(1) + t530) * t507) * MDP(12) + (-t485 * t624 + t507 * t527) * MDP(11) + (t485 * t586 + t612 + (-t443 * t508 - t485 * t616) * qJD(1)) * MDP(14) + ((-t444 * t508 + t485 * t610) * qJD(1) + t534) * MDP(13) + (-t332 * t535 + t346 * t385 + t351 * t601 + t387 * t607 + t469 * t605 + t528) * MDP(25) + (-t493 * t346 - t370 * t535 + t554 * t387 + t601 * t405 + t469 * t651 + t528) * MDP(23) + (t436 * t535 + t469 * t601) * MDP(21) + (-t345 * t535 - t346 * t446 + t387 * t602 - t537 * t601) * MDP(19) + (t436 * t446 + t469 * t602) * MDP(20) + (-t345 * t446 - t537 * t602) * MDP(18) + qJDD(2) * MDP(8) + (t485 * MDP(15) - MDP(20) * t537 + t387 * MDP(21) + t469 * MDP(22) - t347 * MDP(23) + t348 * MDP(24) + t343 * MDP(25) - t344 * MDP(27)) * t593 + (-MDP(4) * t508 * t511 + MDP(5) * t595) * t515; -t444 * t443 * MDP(11) + (-t443 ^ 2 + t444 ^ 2) * MDP(12) + (t527 + t625) * MDP(13) + (-t444 * t485 + t530) * MDP(14) + t440 * MDP(15) + (-t397 * t485 - t444 * t459 + (t561 + t630) * t507 - t547 + t645) * MDP(16) + (g(1) * t431 - g(2) * t429 + g(3) * t615 - t396 * t485 - t443 * t459 - t533) * MDP(17) + (-t349 * t469 + (-t387 * t444 + t436 * t637 + t469 * t584) * pkin(3) + t642) * MDP(23) + (-t350 * t469 + t660 + (-t436 * t506 - t444 * t537 + t469 * t565) * pkin(3) + t524) * MDP(24) + (-t353 * t387 - t436 * t492 + t469 * t553 - t517) * MDP(25) + (-t345 * t492 - t346 * t488 + (t344 + t553) * t537 + (t343 - t597) * t387) * MDP(26) + (t353 * t537 + t436 * t488 - t469 * t597 - t524 + t649 - t661) * MDP(27) + (t330 * t488 + t331 * t492 - t351 * t353 - t343 * t349 - g(3) * (-t486 + t555) + t597 * t344 + (t343 * t584 + t645) * pkin(3) + t521) * MDP(28) + t662; (-t627 + t642) * MDP(23) + (t518 + t660) * MDP(24) + (-t359 * t387 + t425 - t517 - t627) * MDP(25) + (pkin(4) * t345 - qJ(5) * t346 + (t344 - t348) * t537 + (t343 - t581) * t387) * MDP(26) + (t359 * t537 + 0.2e1 * t422 - 0.2e1 * t454 - t518 - t661) * MDP(27) + (-t331 * pkin(4) - g(3) * t555 + t330 * qJ(5) - t343 * t348 + t344 * t581 - t351 * t359 + t521) * MDP(28) + t662; (-t436 + t665) * MDP(25) + t336 * MDP(26) + (-t469 ^ 2 - t639) * MDP(27) + (t344 * t469 + t517) * MDP(28);];
tau = t1;
