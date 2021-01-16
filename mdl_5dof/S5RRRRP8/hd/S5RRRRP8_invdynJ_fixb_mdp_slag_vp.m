% Calculate vector of inverse dynamics joint torques for
% S5RRRRP8
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
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:35
% EndTime: 2021-01-16 00:21:52
% DurationCPUTime: 7.31s
% Computational Cost: add. (5197->527), mult. (11589->666), div. (0->0), fcn. (7987->10), ass. (0->228)
t534 = cos(qJ(3));
t532 = sin(qJ(2));
t535 = cos(qJ(2));
t631 = t534 * t535;
t560 = pkin(3) * t532 - pkin(8) * t631;
t661 = pkin(7) + pkin(8);
t594 = qJD(3) * t661;
t570 = pkin(2) * t532 - pkin(7) * t535;
t472 = t570 * qJD(1);
t531 = sin(qJ(3));
t613 = qJD(1) * t532;
t587 = t531 * t613;
t617 = pkin(6) * t587 + t534 * t472;
t683 = t560 * qJD(1) + t534 * t594 + t617;
t449 = t531 * t472;
t633 = t532 * t534;
t634 = t531 * t535;
t682 = -t449 - (-pkin(6) * t633 - pkin(8) * t634) * qJD(1) - t531 * t594;
t605 = qJD(3) * t532;
t582 = qJD(1) * t605;
t599 = qJD(1) * qJD(2);
t583 = t535 * t599;
t597 = qJDD(1) * t532;
t678 = qJD(2) * qJD(3) + t583 + t597;
t397 = (qJDD(2) - t582) * t531 + t678 * t534;
t608 = qJD(2) * t534;
t464 = -t587 + t608;
t610 = qJD(2) * t531;
t465 = t534 * t613 + t610;
t530 = sin(qJ(4));
t572 = t678 * t531 + t534 * t582;
t548 = t534 * qJDD(2) - t572;
t660 = cos(qJ(4));
t585 = t660 * qJD(4);
t603 = qJD(4) * t530;
t358 = -t660 * t397 - t464 * t585 + t465 * t603 - t530 * t548;
t556 = t530 * t464 + t660 * t465;
t359 = t556 * qJD(4) + t530 * t397 - t660 * t548;
t405 = -t660 * t464 + t465 * t530;
t403 = t405 ^ 2;
t521 = t535 * qJDD(1);
t670 = -t532 * t599 + t521;
t461 = qJDD(3) - t670;
t455 = qJDD(4) + t461;
t612 = qJD(1) * t535;
t507 = -qJD(3) + t612;
t495 = -qJD(4) + t507;
t642 = t556 * t495;
t643 = t405 * t495;
t662 = t556 ^ 2;
t681 = t455 * MDP(22) + t405 * MDP(18) * t556 + (-t359 - t642) * MDP(21) + (-t358 - t643) * MDP(20) + (-t403 + t662) * MDP(19);
t680 = t405 * qJ(5);
t636 = t530 * t531;
t555 = t660 * t534 - t636;
t668 = qJD(3) + qJD(4);
t669 = t660 * qJD(3) + t585;
t622 = -t669 * t534 + t555 * t612 + t668 * t636;
t467 = t530 * t534 + t660 * t531;
t419 = t668 * t467;
t621 = -t467 * t612 + t419;
t518 = pkin(6) * t612;
t606 = qJD(3) * t531;
t657 = pkin(3) * t531;
t674 = pkin(3) * t606 - t612 * t657 - t518;
t607 = qJD(2) * t535;
t588 = t531 * t607;
t604 = qJD(3) * t534;
t679 = t532 * t604 + t588;
t485 = -qJD(2) * pkin(2) + pkin(6) * t613;
t422 = -pkin(3) * t464 + t485;
t529 = qJ(3) + qJ(4);
t522 = sin(t529);
t523 = cos(t529);
t536 = cos(qJ(1));
t533 = sin(qJ(1));
t632 = t533 * t535;
t430 = t522 * t536 - t523 * t632;
t630 = t535 * t536;
t432 = t522 * t533 + t523 * t630;
t475 = t570 * qJD(2);
t481 = -pkin(2) * t535 - pkin(7) * t532 - pkin(1);
t420 = qJD(1) * t475 + t481 * qJDD(1);
t411 = t534 * t420;
t456 = t481 * qJD(1);
t486 = qJD(2) * pkin(7) + t518;
t415 = t456 * t531 + t486 * t534;
t439 = t670 * pkin(6) + qJDD(2) * pkin(7);
t353 = pkin(3) * t461 - pkin(8) * t397 - t415 * qJD(3) - t439 * t531 + t411;
t552 = t531 * t420 + t534 * t439 + t456 * t604 - t486 * t606;
t357 = t548 * pkin(8) + t552;
t414 = t534 * t456 - t486 * t531;
t386 = -pkin(8) * t465 + t414;
t380 = -pkin(3) * t507 + t386;
t387 = pkin(8) * t464 + t415;
t573 = -t530 * t353 - t660 * t357 - t380 * t585 + t387 * t603;
t651 = g(3) * t532;
t545 = g(1) * t432 - g(2) * t430 + t523 * t651 + t573;
t677 = t405 * t422 + t545;
t675 = qJ(5) * t556;
t463 = t534 * t481;
t656 = pkin(6) * t531;
t413 = -pkin(8) * t633 + t463 + (-pkin(3) - t656) * t535;
t509 = pkin(6) * t631;
t616 = t531 * t481 + t509;
t635 = t531 * t532;
t421 = -pkin(8) * t635 + t616;
t623 = t530 * t413 + t660 * t421;
t658 = pkin(3) * t495;
t673 = -t530 * pkin(3) * t455 + t585 * t658;
t487 = t661 * t531;
t488 = t661 * t534;
t618 = -t530 * t487 + t660 * t488;
t672 = t618 * qJD(4) + t682 * t530 + t683 * t660;
t671 = -t487 * t585 - t488 * t603 - t683 * t530 + t682 * t660;
t567 = g(1) * t536 + g(2) * t533;
t516 = pkin(6) * t597;
t646 = qJDD(2) * pkin(2);
t440 = pkin(6) * t583 + t516 - t646;
t650 = g(3) * t535;
t547 = -t567 * t532 + t650;
t667 = -qJD(3) * pkin(7) * t507 + t440 + t547;
t441 = t455 * pkin(4);
t645 = t358 * qJ(5);
t666 = -t556 * qJD(5) + t441 + t645;
t637 = t523 * t536;
t429 = t522 * t632 + t637;
t638 = t523 * t533;
t431 = -t522 * t630 + t638;
t665 = -g(1) * t431 + g(2) * t429 + t522 * t651;
t385 = t660 * t387;
t361 = t530 * t380 + t385;
t580 = t660 * t353 - t530 * t357;
t543 = -t361 * qJD(4) + t580;
t539 = t543 + t665;
t664 = -t422 * t556 + t539;
t581 = -pkin(4) * t405 - qJD(5);
t376 = t422 - t581;
t663 = -t376 * t556 + t665;
t659 = pkin(3) * t465;
t525 = t534 * pkin(3);
t649 = pkin(2) + t525;
t479 = pkin(4) * t522 + t657;
t648 = pkin(6) + t479;
t647 = qJ(5) * t359;
t644 = t397 * t531;
t641 = t464 * t507;
t640 = t465 * t507;
t639 = t465 * t534;
t383 = t530 * t387;
t360 = t660 * t380 - t383;
t351 = t360 - t675;
t348 = -pkin(4) * t495 + t351;
t629 = -t351 + t348;
t628 = -t621 * qJ(5) + qJD(5) * t555 + t671;
t627 = -pkin(4) * t613 + t622 * qJ(5) - t467 * qJD(5) - t672;
t626 = t660 * t386 - t383;
t625 = t621 * pkin(4) + t674;
t620 = t531 * t475 + t481 * t604;
t609 = qJD(2) * t532;
t619 = t534 * t475 + t609 * t656;
t480 = pkin(4) * t523 + t525;
t476 = pkin(3) * t635 + t532 * pkin(6);
t527 = t532 ^ 2;
t615 = -t535 ^ 2 + t527;
t611 = qJD(2) * t465;
t602 = qJD(5) * t405;
t600 = t464 * qJD(2);
t423 = t679 * pkin(3) + pkin(6) * t607;
t593 = t507 * t608;
t592 = t507 * t606;
t591 = t531 * t605;
t590 = t507 * t604;
t579 = -t386 * t530 - t385;
t576 = t660 * t413 - t421 * t530;
t575 = -t660 * t487 - t488 * t530;
t574 = -qJD(3) * t456 - t439;
t571 = t660 * t607;
t569 = -g(1) * t429 - g(2) * t431;
t568 = -g(1) * t430 - g(2) * t432;
t566 = g(1) * t533 - g(2) * t536;
t565 = t486 * t604 - t411;
t564 = -pkin(7) * t461 + qJD(3) * t485;
t471 = pkin(2) + t480;
t526 = -qJ(5) - t661;
t563 = t471 * t535 - t526 * t532;
t561 = -t523 * t650 + (g(1) * t637 + g(2) * t638) * t532;
t559 = pkin(1) + t563;
t558 = t567 * t522;
t557 = -0.2e1 * pkin(1) * t599 - pkin(6) * qJDD(2);
t554 = t531 * t461 - t590;
t553 = t461 * t534 + t592;
t371 = t560 * qJD(2) + (-t509 + (pkin(8) * t532 - t481) * t531) * qJD(3) + t619;
t373 = -t679 * pkin(8) + (-t532 * t608 - t535 * t606) * pkin(6) + t620;
t551 = t530 * t371 + t660 * t373 + t413 * t585 - t421 * t603;
t538 = qJD(1) ^ 2;
t549 = pkin(1) * t538 + t567;
t537 = qJD(2) ^ 2;
t546 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t537 + t566;
t542 = -t623 * qJD(4) + t660 * t371 - t530 * t373;
t540 = t545 + t647;
t379 = -t548 * pkin(3) + t440;
t347 = t359 * pkin(4) + qJDD(5) + t379;
t514 = t660 * pkin(3) + pkin(4);
t494 = t522 * t650;
t447 = t531 * t533 + t534 * t630;
t446 = -t531 * t630 + t533 * t534;
t445 = t531 * t536 - t533 * t631;
t444 = t531 * t632 + t534 * t536;
t436 = t555 * t532;
t435 = t467 * t532;
t426 = -pkin(4) * t555 - t649;
t416 = pkin(4) * t435 + t476;
t391 = qJ(5) * t555 + t618;
t390 = -qJ(5) * t467 + t575;
t382 = pkin(4) * t556 + t659;
t375 = t531 * t571 - t530 * t591 - t603 * t635 + (t530 * t607 + t669 * t532) * t534;
t374 = t419 * t532 + t530 * t588 - t534 * t571;
t368 = pkin(4) * t375 + t423;
t367 = -qJ(5) * t435 + t623;
t366 = -pkin(4) * t535 - qJ(5) * t436 + t576;
t355 = t626 - t675;
t354 = t579 + t680;
t352 = t361 - t680;
t346 = -qJ(5) * t375 - qJD(5) * t435 + t551;
t345 = pkin(4) * t609 + t374 * qJ(5) - t436 * qJD(5) + t542;
t344 = -t573 - t602 - t647;
t343 = t543 + t666;
t1 = [(t343 * t366 + t344 * t367 + t348 * t345 + t352 * t346 + t347 * t416 + t376 * t368 + (-g(1) * t648 - g(2) * t559) * t536 + (g(1) * t559 - g(2) * t648) * t533) * MDP(28) + ((t464 * t534 - t465 * t531) * t607 + (t534 * t548 - t644 + (-t531 * t464 - t639) * qJD(3)) * t532) * MDP(12) + (t397 * t633 + (t534 * t607 - t591) * t465) * MDP(11) + 0.2e1 * (t532 * t521 - t615 * t599) * MDP(5) + (-(-t481 * t606 + t619) * t507 + t463 * t461 - g(1) * t445 - g(2) * t447 + ((t590 - t600) * pkin(6) + (-pkin(6) * t461 + qJD(2) * t485 - t574) * t531 + t565) * t535 + (-pkin(6) * t548 + t414 * qJD(2) + t440 * t531 + t485 * t604) * t532) * MDP(16) + (t620 * t507 - t616 * t461 - g(1) * t444 - g(2) * t446 + (t485 * t608 + (-t592 + t611) * pkin(6) + t552) * t535 + (-t485 * t606 - t415 * qJD(2) + t440 * t534 + (t397 - t593) * pkin(6)) * t532) * MDP(17) + (t476 * t359 + t360 * t609 + t422 * t375 + t379 * t435 + t423 * t405 + t576 * t455 - t542 * t495 - t543 * t535 + t568) * MDP(23) + (-t343 * t535 - t345 * t495 + t347 * t435 + t348 * t609 + t359 * t416 + t366 * t455 + t368 * t405 + t375 * t376 + t568) * MDP(25) + ((t507 * t610 - t548) * t535 + (-t554 + t600) * t532) * MDP(14) + ((-t397 - t593) * t535 + (t553 + t611) * t532) * MDP(13) + (-t461 * t535 - t507 * t609) * MDP(15) + (-t455 * t535 - t495 * t609) * MDP(22) + (t359 * t535 + t375 * t495 - t405 * t609 - t435 * t455) * MDP(21) + (qJDD(1) * t527 + 0.2e1 * t532 * t583) * MDP(4) + (t557 * t532 + t546 * t535) * MDP(9) + (-t546 * t532 + t557 * t535) * MDP(10) + (-t476 * t358 - t361 * t609 - t422 * t374 + t379 * t436 + t423 * t556 - t623 * t455 + t551 * t495 - t573 * t535 + t569) * MDP(24) + (t344 * t535 + t346 * t495 + t347 * t436 - t352 * t609 - t358 * t416 - t367 * t455 + t368 * t556 - t374 * t376 + t569) * MDP(26) + (t358 * t535 + t374 * t495 + t436 * t455 + t556 * t609) * MDP(20) + (-t343 * t436 - t344 * t435 - t345 * t556 - t346 * t405 + t348 * t374 - t352 * t375 + t358 * t366 - t359 * t367 + t566 * t532) * MDP(27) + (-t358 * t436 - t374 * t556) * MDP(18) + (t358 * t435 - t359 * t436 + t374 * t405 - t375 * t556) * MDP(19) + qJDD(1) * MDP(1) + t566 * MDP(2) + t567 * MDP(3) + (qJDD(2) * t532 + t535 * t537) * MDP(6) + (qJDD(2) * t535 - t532 * t537) * MDP(7); (-pkin(2) * t572 + t617 * t507 + t564 * t531 + (-t414 * t532 + (pkin(6) * t464 - t485 * t531) * t535) * qJD(1) + (t646 - t667) * t534) * MDP(16) + (t549 * t532 - t516 - t650) * MDP(9) + (-t343 * t467 + t344 * t555 + t622 * t348 - t621 * t352 + t358 * t390 - t359 * t391 - t628 * t405 - t567 * t535 - t556 * t627 - t651) * MDP(27) + (t651 + (-pkin(6) * qJDD(1) + t549) * t535) * MDP(10) + (-g(3) * t563 + t343 * t390 + t344 * t391 + t347 * t426 + t627 * t348 + t628 * t352 + t625 * t376 + t567 * (t471 * t532 + t526 * t535)) * MDP(28) + (-t507 * t639 + t644) * MDP(11) + ((t397 - t641) * t534 + (t548 + t640) * t531) * MDP(12) + ((-t464 * t532 - t507 * t634) * qJD(1) + t553) * MDP(14) + (-pkin(2) * t397 - t449 * t507 + t564 * t534 + (-t485 * t631 + t415 * t532 + (-t465 * t535 + t507 * t633) * pkin(6)) * qJD(1) + t667 * t531) * MDP(17) + ((-t465 * t532 + t507 * t631) * qJD(1) + t554) * MDP(13) + (t455 * t467 - t556 * t613) * MDP(20) + (-t358 * t467 - t556 * t622) * MDP(18) + (-t618 * t455 + t649 * t358 + t379 * t467 + t494 - t622 * t422 + t674 * t556 + (t361 * qJD(1) - t558) * t532) * MDP(24) + (-t347 * t555 - t348 * t613 + t359 * t426 + t621 * t376 + t390 * t455 + t625 * t405 + t561) * MDP(25) + (t347 * t467 - t358 * t426 - t391 * t455 + t494 + t625 * t556 - t622 * t376 + (qJD(1) * t352 - t558) * t532) * MDP(26) + (t405 * t613 + t455 * t555) * MDP(21) + (-t358 * t555 - t359 * t467 + t622 * t405 - t556 * t621) * MDP(19) + (-t649 * t359 - t360 * t613 - t379 * t555 + t674 * t405 + t621 * t422 + t575 * t455 + t561) * MDP(23) + t507 * MDP(15) * t613 + qJDD(2) * MDP(8) + MDP(7) * t521 + MDP(6) * t597 + (-t532 * t535 * MDP(4) + t615 * MDP(5)) * t538 + (t622 * MDP(20) + t621 * MDP(21) + MDP(22) * t613 + t672 * MDP(23) + t671 * MDP(24) - t627 * MDP(25) + t628 * MDP(26)) * t495; -t465 * t464 * MDP(11) + (-t464 ^ 2 + t465 ^ 2) * MDP(12) + (t397 + t641) * MDP(13) + (t548 - t640) * MDP(14) + t461 * MDP(15) + (-g(1) * t446 + g(2) * t444 - t415 * t507 - t465 * t485 + (t574 + t651) * t531 - t565) * MDP(16) + (g(1) * t447 - g(2) * t445 + g(3) * t633 - t414 * t507 - t464 * t485 - t552) * MDP(17) + (t579 * t495 + (-t465 * t405 + t660 * t455 + t495 * t603) * pkin(3) + t664) * MDP(23) + (-t626 * t495 - t556 * t659 + t673 + t677) * MDP(24) + (t354 * t495 - t382 * t405 + t514 * t455 + (-t385 + (-t380 + t658) * t530) * qJD(4) + t580 + t663 + t666) * MDP(25) + (-t355 * t495 + t376 * t405 - t382 * t556 + t540 + t602 + t673) * MDP(26) + (-t348 * t405 + t352 * t556 + t354 * t556 + t355 * t405 + t514 * t358 + (-t359 * t530 + (-t660 * t405 + t530 * t556) * qJD(4)) * pkin(3)) * MDP(27) + (t343 * t514 - t352 * t355 - t348 * t354 - t376 * t382 - g(1) * (-t479 * t630 + t480 * t533) - g(2) * (-t479 * t632 - t480 * t536) + t479 * t651 + (t344 * t530 + (-t348 * t530 + t660 * t352) * qJD(4)) * pkin(3)) * MDP(28) + t681; (-t361 * t495 + t664) * MDP(23) + (-t360 * t495 + t677) * MDP(24) + (t645 - t352 * t495 + 0.2e1 * t441 + (-t376 + t581) * t556 + t539) * MDP(25) + (-pkin(4) * t662 - t351 * t495 + (qJD(5) + t376) * t405 + t540) * MDP(26) + (pkin(4) * t358 - t629 * t405) * MDP(27) + (t629 * t352 + (t343 + t663) * pkin(4)) * MDP(28) + t681; (t359 - t642) * MDP(25) + (-t358 + t643) * MDP(26) + (-t403 - t662) * MDP(27) + (t348 * t556 + t352 * t405 + t347 + t547) * MDP(28);];
tau = t1;
