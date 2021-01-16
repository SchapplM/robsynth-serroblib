% Calculate vector of inverse dynamics joint torques for
% S5RRRPR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:24:21
% EndTime: 2021-01-15 23:24:48
% DurationCPUTime: 9.58s
% Computational Cost: add. (5227->544), mult. (11903->721), div. (0->0), fcn. (8459->14), ass. (0->224)
t566 = sin(qJ(3));
t567 = sin(qJ(2));
t631 = qJD(1) * t567;
t613 = t566 * t631;
t570 = cos(qJ(3));
t619 = t570 * qJD(2);
t513 = t613 - t619;
t627 = qJD(2) * t566;
t515 = t570 * t631 + t627;
t562 = sin(pkin(9));
t563 = cos(pkin(9));
t441 = t563 * t513 + t515 * t562;
t569 = cos(qJ(5));
t647 = t569 * t441;
t565 = sin(qJ(5));
t589 = -t513 * t562 + t563 * t515;
t657 = t589 * t565;
t392 = t647 + t657;
t571 = cos(qJ(2));
t630 = qJD(1) * t571;
t537 = -qJD(3) + t630;
t529 = -qJD(5) + t537;
t659 = t392 * t529;
t618 = qJD(1) * qJD(2);
t606 = t571 * t618;
t617 = qJDD(1) * t567;
t623 = qJD(3) * t567;
t675 = -qJD(1) * t623 + qJDD(2);
t433 = qJD(3) * t619 + (t606 + t617) * t570 + t675 * t566;
t597 = pkin(2) * t567 - pkin(7) * t571;
t517 = t597 * qJD(2);
t522 = -pkin(2) * t571 - pkin(7) * t567 - pkin(1);
t452 = qJD(1) * t517 + qJDD(1) * t522;
t446 = t570 * t452;
t502 = t522 * qJD(1);
t549 = pkin(6) * t630;
t526 = qJD(2) * pkin(7) + t549;
t451 = t502 * t566 + t526 * t570;
t555 = t571 * qJDD(1);
t671 = -t567 * t618 + t555;
t485 = pkin(6) * t671 + qJDD(2) * pkin(7);
t505 = qJDD(3) - t671;
t359 = pkin(3) * t505 - qJ(4) * t433 - qJD(3) * t451 - qJD(4) * t515 - t485 * t566 + t446;
t434 = t566 * (qJD(2) * (qJD(3) + t630) + t617) - t675 * t570;
t622 = qJD(3) * t570;
t624 = qJD(3) * t566;
t581 = t566 * t452 + t570 * t485 + t502 * t622 - t526 * t624;
t364 = -qJ(4) * t434 - qJD(4) * t513 + t581;
t347 = t563 * t359 - t364 * t562;
t382 = t433 * t563 - t434 * t562;
t345 = pkin(4) * t505 - pkin(8) * t382 + t347;
t348 = t562 * t359 + t563 * t364;
t381 = t433 * t562 + t563 * t434;
t346 = -pkin(8) * t381 + t348;
t450 = t570 * t502 - t526 * t566;
t414 = -qJ(4) * t515 + t450;
t406 = -pkin(3) * t537 + t414;
t415 = -qJ(4) * t513 + t451;
t653 = t563 * t415;
t371 = t562 * t406 + t653;
t677 = pkin(8) * t441;
t363 = t371 - t677;
t620 = qJD(5) * t565;
t361 = t363 * t620;
t525 = -qJD(2) * pkin(2) + pkin(6) * t631;
t456 = pkin(3) * t513 + qJD(4) + t525;
t400 = pkin(4) * t441 + t456;
t559 = qJ(3) + pkin(9);
t556 = qJ(5) + t559;
t541 = sin(t556);
t542 = cos(t556);
t572 = cos(qJ(1));
t568 = sin(qJ(1));
t649 = t568 * t571;
t461 = t541 * t572 - t542 * t649;
t645 = t571 * t572;
t463 = t541 * t568 + t542 * t645;
t662 = g(3) * t567;
t684 = g(1) * t463 - g(2) * t461 - t565 * t345 - t569 * t346 + t392 * t400 + t542 * t662 + t361;
t501 = qJDD(5) + t505;
t590 = t441 * t565 - t569 * t589;
t683 = t501 * MDP(26) + (-t392 ^ 2 + t590 ^ 2) * MDP(23) - t392 * MDP(22) * t590;
t660 = t529 * t590;
t460 = t541 * t649 + t542 * t572;
t462 = -t541 * t645 + t542 * t568;
t603 = t569 * t345 - t565 * t346;
t681 = -g(1) * t462 + g(2) * t460 + t400 * t590 + t541 * t662 + t603;
t646 = t570 * t571;
t666 = pkin(3) * t567;
t588 = -qJ(4) * t646 + t666;
t564 = qJ(4) + pkin(7);
t604 = qJD(3) * t564;
t516 = t597 * qJD(1);
t635 = pkin(6) * t613 + t570 * t516;
t680 = -qJD(1) * t588 - qJD(4) * t566 - t570 * t604 - t635;
t497 = t566 * t516;
t621 = qJD(4) * t570;
t650 = t567 * t570;
t651 = t566 * t571;
t679 = t497 + (-pkin(6) * t650 - qJ(4) * t651) * qJD(1) + t566 * t604 - t621;
t547 = pkin(6) * t617;
t486 = -qJDD(2) * pkin(2) + pkin(6) * t606 + t547;
t596 = g(1) * t572 + g(2) * t568;
t661 = g(3) * t571;
t579 = t567 * t596 - t661;
t678 = pkin(7) * qJD(3) * t537 - t486 + t579;
t507 = t562 * t570 + t563 * t566;
t490 = t507 * qJD(3);
t639 = t507 * t630 - t490;
t506 = t562 * t566 - t563 * t570;
t676 = t537 * t506;
t674 = pkin(8) * t589;
t642 = t679 * t562 + t680 * t563;
t641 = t680 * t562 - t679 * t563;
t667 = pkin(3) * t566;
t598 = pkin(3) * t624 - t630 * t667 - t549;
t493 = t566 * t649 + t570 * t572;
t495 = -t566 * t645 + t568 * t570;
t670 = -g(1) * t495 + g(2) * t493;
t602 = t569 * t381 + t382 * t565;
t352 = -qJD(5) * t590 + t602;
t668 = pkin(3) * t562;
t665 = pkin(6) * t566;
t658 = t433 * t566;
t656 = t513 * t537;
t655 = t515 * t537;
t654 = t515 * t570;
t410 = t562 * t415;
t652 = t566 * t567;
t370 = t563 * t406 - t410;
t360 = -pkin(4) * t537 + t370 - t674;
t648 = t569 * t360;
t438 = t569 * t506 + t507 * t565;
t644 = -qJD(5) * t438 + t565 * t639 + t569 * t676;
t439 = -t506 * t565 + t507 * t569;
t643 = qJD(5) * t439 + t565 * t676 - t569 * t639;
t539 = pkin(6) * t646;
t626 = qJD(2) * t567;
t636 = t570 * t517 + t626 * t665;
t387 = -t567 * t621 + t588 * qJD(2) + (-t539 + (qJ(4) * t567 - t522) * t566) * qJD(3) + t636;
t637 = t566 * t517 + t522 * t622;
t396 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t650 + (-qJD(4) * t567 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t571) * t566 + t637;
t358 = t562 * t387 + t563 * t396;
t375 = t563 * t414 - t410;
t509 = t570 * t522;
t447 = -qJ(4) * t650 + t509 + (-pkin(3) - t665) * t571;
t634 = t566 * t522 + t539;
t453 = -qJ(4) * t652 + t634;
t398 = t562 * t447 + t563 * t453;
t640 = -pkin(4) * t639 + t598;
t523 = t564 * t566;
t524 = t564 * t570;
t455 = -t562 * t523 + t563 * t524;
t518 = pkin(3) * t652 + t567 * pkin(6);
t560 = t567 ^ 2;
t633 = -t571 ^ 2 + t560;
t629 = qJD(2) * t513;
t628 = qJD(2) * t515;
t625 = qJD(2) * t571;
t614 = -qJD(5) * t647 - t565 * t381 + t569 * t382;
t611 = t566 * t625;
t459 = pkin(3) * t611 + pkin(6) * t625 + t622 * t666;
t546 = pkin(3) * t570 + pkin(2);
t612 = t537 * t619;
t610 = t571 * t619;
t609 = t537 * t624;
t608 = t537 * t622;
t357 = t563 * t387 - t396 * t562;
t374 = -t414 * t562 - t653;
t397 = t563 * t447 - t453 * t562;
t454 = -t563 * t523 - t524 * t562;
t601 = t546 * t571 + t564 * t567;
t600 = -qJD(3) * t502 - t485;
t595 = g(1) * t568 - g(2) * t572;
t594 = t526 * t622 - t446;
t426 = -pkin(8) * t506 + t455;
t593 = pkin(4) * t631 + pkin(8) * t676 + qJD(5) * t426 - t642;
t425 = -pkin(8) * t507 + t454;
t592 = pkin(8) * t639 + qJD(5) * t425 + t641;
t591 = -pkin(7) * t505 + qJD(3) * t525;
t350 = t565 * t360 + t569 * t363;
t478 = t507 * t567;
t479 = t506 * t567;
t417 = t569 * t478 - t479 * t565;
t418 = -t478 * t565 - t479 * t569;
t543 = pkin(3) * t563 + pkin(4);
t587 = t543 * t565 + t569 * t668;
t586 = t543 * t569 - t565 * t668;
t584 = -0.2e1 * pkin(1) * t618 - pkin(6) * qJDD(2);
t583 = t505 * t566 - t608;
t582 = t505 * t570 + t609;
t351 = -t589 * t620 + t614;
t574 = qJD(1) ^ 2;
t580 = pkin(1) * t574 + t596;
t405 = pkin(3) * t434 + qJDD(4) + t486;
t573 = qJD(2) ^ 2;
t576 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t573 + t595;
t552 = cos(t559);
t551 = sin(t559);
t545 = pkin(6) + t667;
t496 = t566 * t568 + t570 * t645;
t494 = t566 * t572 - t568 * t646;
t484 = pkin(1) + t601;
t475 = t551 * t568 + t552 * t645;
t474 = -t551 * t645 + t552 * t568;
t473 = t551 * t572 - t552 * t649;
t472 = t551 * t649 + t552 * t572;
t467 = pkin(4) * t506 - t546;
t448 = pkin(4) * t478 + t518;
t420 = t490 * t567 + t562 * t611 - t563 * t610;
t419 = t506 * t623 - t507 * t625;
t407 = pkin(3) * t515 + pkin(4) * t589;
t399 = -pkin(4) * t419 + t459;
t377 = -pkin(8) * t478 + t398;
t376 = -pkin(4) * t571 + pkin(8) * t479 + t397;
t369 = t375 - t674;
t368 = t374 + t677;
t367 = qJD(5) * t418 - t569 * t419 - t420 * t565;
t366 = -qJD(5) * t417 + t419 * t565 - t420 * t569;
t365 = pkin(4) * t381 + t405;
t354 = pkin(8) * t419 + t358;
t353 = pkin(4) * t626 + pkin(8) * t420 + t357;
t349 = -t363 * t565 + t648;
t1 = [(qJDD(1) * t560 + 0.2e1 * t567 * t606) * MDP(4) + ((-t513 * t570 - t515 * t566) * t625 + (-t658 - t434 * t570 + (t513 * t566 - t654) * qJD(3)) * t567) * MDP(12) + (t433 * t650 + (-t566 * t623 + t610) * t515) * MDP(11) + (-t350 * t626 - g(1) * t460 - g(2) * t462 + t448 * t351 - t361 * t571 + t365 * t418 + t400 * t366 - t399 * t590 + ((-qJD(5) * t377 + t353) * t529 - t376 * t501 + t345 * t571) * t565 + ((qJD(5) * t376 + t354) * t529 - t377 * t501 + (qJD(5) * t360 + t346) * t571) * t569) * MDP(28) + (-t351 * t571 - t366 * t529 + t418 * t501 - t590 * t626) * MDP(24) + (t351 * t418 - t366 * t590) * MDP(22) + (-t351 * t417 - t352 * t418 - t366 * t392 + t367 * t590) * MDP(23) + (-(-t522 * t624 + t636) * t537 + t509 * t505 - g(1) * t494 - g(2) * t496 + ((t608 + t629) * pkin(6) + (-pkin(6) * t505 + qJD(2) * t525 - t600) * t566 + t594) * t571 + (pkin(6) * t434 + qJD(2) * t450 + t486 * t566 + t525 * t622) * t567) * MDP(16) + (t637 * t537 - t634 * t505 - g(1) * t493 - g(2) * t495 + (t525 * t619 + (-t609 + t628) * pkin(6) + t581) * t571 + (-t525 * t624 - t451 * qJD(2) + t486 * t570 + (t433 - t612) * pkin(6)) * t567) * MDP(17) + 0.2e1 * (t555 * t567 - t618 * t633) * MDP(5) + (-t505 * t571 - t537 * t626) * MDP(15) + (-t501 * t571 - t529 * t626) * MDP(26) + (t352 * t571 + t367 * t529 - t392 * t626 - t417 * t501) * MDP(25) + (-(t353 * t569 - t354 * t565) * t529 + (t376 * t569 - t377 * t565) * t501 - t603 * t571 + t349 * t626 + t399 * t392 + t448 * t352 + t365 * t417 + t400 * t367 - g(1) * t461 - g(2) * t463 + (-(-t376 * t565 - t377 * t569) * t529 + t350 * t571) * qJD(5)) * MDP(27) + (-g(1) * t473 - g(2) * t475 - t347 * t571 - t357 * t537 + t370 * t626 + t381 * t518 + t397 * t505 + t405 * t478 - t419 * t456 + t441 * t459) * MDP(18) + ((-t433 - t612) * t571 + (t582 + t628) * t567) * MDP(13) + ((t537 * t627 + t434) * t571 + (-t583 - t629) * t567) * MDP(14) + t595 * MDP(2) + t596 * MDP(3) + (-g(1) * t472 - g(2) * t474 + t348 * t571 + t358 * t537 - t371 * t626 + t382 * t518 - t398 * t505 - t405 * t479 - t420 * t456 + t459 * t589) * MDP(19) + (t347 * t479 - t348 * t478 - t357 * t589 - t358 * t441 + t370 * t420 + t371 * t419 - t381 * t398 - t382 * t397 + t567 * t595) * MDP(20) + (t348 * t398 + t371 * t358 + t347 * t397 + t370 * t357 + t405 * t518 + t456 * t459 - g(1) * (-t484 * t568 + t545 * t572) - g(2) * (t484 * t572 + t545 * t568)) * MDP(21) + (qJDD(2) * t567 + t571 * t573) * MDP(6) + (qJDD(2) * t571 - t567 * t573) * MDP(7) + (t567 * t584 + t571 * t576) * MDP(9) + (-t567 * t576 + t571 * t584) * MDP(10) + qJDD(1) * MDP(1); (-t537 * t654 + t658) * MDP(11) + ((t433 + t656) * t570 + (-t434 + t655) * t566) * MDP(12) + (-pkin(2) * t433 - t497 * t537 + t591 * t570 + (-t525 * t646 + t451 * t567 + (-t515 * t571 + t537 * t650) * pkin(6)) * qJD(1) - t678 * t566) * MDP(17) + ((t513 * t567 - t537 * t651) * qJD(1) + t582) * MDP(14) + ((-t515 * t567 + t537 * t646) * qJD(1) + t583) * MDP(13) + (-t382 * t546 + t405 * t507 - t455 * t505 + t456 * t676 + t537 * t641 - t551 * t579 + t589 * t598) * MDP(19) + (t348 * t455 + t347 * t454 - t405 * t546 - g(3) * t601 - t596 * (-t546 * t567 + t564 * t571) + t598 * t456 + t641 * t371 + t642 * t370) * MDP(21) + (-t381 * t546 + t405 * t506 + t441 * t598 + t454 * t505 - t456 * t639 - t537 * t642 + t552 * t579) * MDP(18) + ((t425 * t569 - t426 * t565) * t501 + t467 * t352 + t365 * t438 + (t565 * t592 + t569 * t593) * t529 + t643 * t400 + t640 * t392 + t579 * t542) * MDP(27) + (-t438 * t501 + t529 * t643) * MDP(25) + (-(t425 * t565 + t426 * t569) * t501 + t467 * t351 + t365 * t439 + (-t565 * t593 + t569 * t592) * t529 + t644 * t400 - t640 * t590 - t579 * t541) * MDP(28) + (t351 * t439 - t590 * t644) * MDP(22) + (-t351 * t438 - t352 * t439 - t392 * t644 + t590 * t643) * MDP(23) + (t439 * t501 - t529 * t644) * MDP(24) + (-pkin(2) * t434 + t635 * t537 + t591 * t566 + (-t450 * t567 + (-pkin(6) * t513 - t525 * t566) * t571) * qJD(1) + t678 * t570) * MDP(16) + (t567 * t580 - t547 - t661) * MDP(9) + (-t347 * t507 - t348 * t506 - t370 * t676 + t371 * t639 - t381 * t455 - t382 * t454 - t441 * t641 - t571 * t596 - t589 * t642 - t662) * MDP(20) + (t662 + (-pkin(6) * qJDD(1) + t580) * t571) * MDP(10) + qJDD(2) * MDP(8) + MDP(7) * t555 + MDP(6) * t617 + (t537 * MDP(15) - t370 * MDP(18) + t371 * MDP(19) + MDP(24) * t590 + t392 * MDP(25) + t529 * MDP(26) - t349 * MDP(27) + t350 * MDP(28)) * t631 + (-MDP(4) * t567 * t571 + MDP(5) * t633) * t574; t515 * t513 * MDP(11) + (-t513 ^ 2 + t515 ^ 2) * MDP(12) + (t433 - t656) * MDP(13) + (-t434 - t655) * MDP(14) + t505 * MDP(15) + (-t451 * t537 - t515 * t525 + (t600 + t662) * t566 - t594 + t670) * MDP(16) + (g(1) * t496 - g(2) * t494 + g(3) * t650 - t450 * t537 + t513 * t525 - t581) * MDP(17) + (t551 * t662 - g(1) * t474 + g(2) * t472 + t374 * t537 - t589 * t456 + (-t441 * t515 + t505 * t563) * pkin(3) + t347) * MDP(18) + (t552 * t662 + g(1) * t475 - g(2) * t473 - t375 * t537 + t441 * t456 + (-t505 * t562 - t515 * t589) * pkin(3) - t348) * MDP(19) + ((-t381 * t562 - t382 * t563) * pkin(3) + (t371 + t374) * t589 + (-t370 + t375) * t441) * MDP(20) + (-t370 * t374 - t371 * t375 + (g(3) * t652 + t347 * t563 + t348 * t562 - t456 * t515 + t670) * pkin(3)) * MDP(21) + (t351 - t659) * MDP(24) + (-t352 + t660) * MDP(25) + (t586 * t501 + (t368 * t569 - t369 * t565) * t529 - t407 * t392 + (t529 * t587 - t350) * qJD(5) + t681) * MDP(27) + (-t587 * t501 - (t368 * t565 + t369 * t569) * t529 + t407 * t590 + (t529 * t586 - t648) * qJD(5) + t684) * MDP(28) + t683; (-t537 * t589 + t381) * MDP(18) + (t441 * t537 + t382) * MDP(19) + (-t441 ^ 2 - t589 ^ 2) * MDP(20) + (t370 * t589 + t371 * t441 + t405 - t579) * MDP(21) + (t352 + t660) * MDP(27) + (t351 + t659) * MDP(28); (t614 - t659) * MDP(24) + (-t602 + t660) * MDP(25) + (-t350 * t529 + t681) * MDP(27) + (-t349 * t529 + t684) * MDP(28) + (-MDP(24) * t657 + MDP(25) * t590 - MDP(27) * t350 - MDP(28) * t648) * qJD(5) + t683;];
tau = t1;
