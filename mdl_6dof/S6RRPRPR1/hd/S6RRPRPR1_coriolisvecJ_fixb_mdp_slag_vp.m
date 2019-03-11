% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:10:13
% EndTime: 2019-03-09 10:10:25
% DurationCPUTime: 8.48s
% Computational Cost: add. (9128->453), mult. (23723->615), div. (0->0), fcn. (18567->10), ass. (0->217)
t566 = sin(pkin(10));
t568 = cos(pkin(10));
t571 = sin(qJ(2));
t573 = cos(qJ(2));
t537 = -t566 * t571 + t568 * t573;
t525 = t537 * qJD(1);
t539 = t566 * t573 + t568 * t571;
t527 = t539 * qJD(1);
t570 = sin(qJ(4));
t654 = cos(qJ(4));
t482 = t525 * t654 - t527 * t570;
t478 = qJD(6) - t482;
t567 = cos(pkin(11));
t572 = cos(qJ(6));
t628 = t572 * t567;
t565 = sin(pkin(11));
t569 = sin(qJ(6));
t631 = t565 * t569;
t538 = -t628 + t631;
t685 = t478 * t538;
t540 = t565 * t572 + t567 * t569;
t608 = qJD(1) * qJD(2);
t603 = t573 * t608;
t604 = t571 * t608;
t517 = -t566 * t604 + t568 * t603;
t526 = t539 * qJD(2);
t580 = qJD(1) * t526;
t588 = -t525 * t570 - t527 * t654;
t577 = qJD(4) * t588 - t517 * t570 - t580 * t654;
t687 = -t478 * t685 - t540 * t577;
t686 = t478 * t540;
t562 = qJD(2) + qJD(4);
t466 = t562 * t565 - t567 * t588;
t605 = qJD(4) * t654;
t614 = qJD(4) * t570;
t439 = t517 * t654 + t525 * t605 - t527 * t614 - t570 * t580;
t644 = t439 * t565;
t590 = -qJD(6) * t466 - t644;
t464 = -t562 * t567 - t565 * t588;
t610 = qJD(6) * t572;
t619 = t439 * t628 - t464 * t610;
t366 = t569 * t590 + t619;
t662 = t464 * t569 - t466 * t572;
t367 = -qJD(6) * t662 + t439 * t540;
t670 = t572 * t464;
t418 = t466 * t569 + t670;
t596 = -t478 * t686 + t538 * t577;
t638 = t482 * t562;
t641 = t588 * t562;
t645 = t662 * t588;
t646 = t418 * t588;
t684 = (t596 - t646) * MDP(27) - t482 ^ 2 * MDP(14) + (MDP(13) * t482 + MDP(14) * t588 + MDP(28) * t478) * t588 + (t439 - t638) * MDP(15) + (t577 - t641) * MDP(16) + t366 * t540 * MDP(24) + (-t645 + t687) * MDP(26) + (-t366 * t538 - t540 * t367 + t418 * t685) * MDP(25) + (MDP(24) * t685 + MDP(25) * t686) * t662;
t674 = t482 * t565;
t683 = pkin(5) * t674;
t682 = pkin(9) * t674;
t673 = t482 * t567;
t680 = -pkin(5) * t588 - pkin(9) * t673;
t649 = -qJ(3) - pkin(7);
t602 = qJD(2) * t649;
t522 = qJD(3) * t573 + t571 * t602;
t506 = t522 * qJD(1);
t523 = -qJD(3) * t571 + t573 * t602;
t507 = t523 * qJD(1);
t468 = -t506 * t566 + t507 * t568;
t448 = -pkin(8) * t517 + t468;
t469 = t506 * t568 + t507 * t566;
t449 = -pkin(8) * t580 + t469;
t551 = t649 * t573;
t545 = qJD(1) * t551;
t532 = t566 * t545;
t550 = t649 * t571;
t544 = qJD(1) * t550;
t648 = qJD(2) * pkin(2);
t536 = t544 + t648;
t484 = t536 * t568 + t532;
t651 = pkin(8) * t527;
t460 = qJD(2) * pkin(3) + t484 - t651;
t630 = t568 * t545;
t485 = t536 * t566 - t630;
t652 = pkin(8) * t525;
t463 = t485 + t652;
t576 = t570 * t448 + t449 * t654 + t460 * t605 - t463 * t614;
t371 = qJD(5) * t562 + t576;
t555 = pkin(2) * t604;
t491 = pkin(3) * t580 + t555;
t381 = -pkin(4) * t577 - qJ(5) * t439 + qJD(5) * t588 + t491;
t353 = t371 * t567 + t381 * t565;
t351 = t353 * t567;
t411 = t460 * t570 + t463 * t654;
t408 = qJ(5) * t562 + t411;
t607 = -pkin(2) * t573 - pkin(1);
t595 = t607 * qJD(1);
t547 = qJD(3) + t595;
t490 = -pkin(3) * t525 + t547;
t415 = -pkin(4) * t482 + qJ(5) * t588 + t490;
t382 = -t408 * t565 + t415 * t567;
t679 = t382 * t673 + t351;
t360 = -pkin(5) * t482 - pkin(9) * t466 + t382;
t383 = t408 * t567 + t415 * t565;
t368 = -pkin(9) * t464 + t383;
t348 = t360 * t569 + t368 * t572;
t373 = -t448 * t654 + t449 * t570 + t460 * t614 + t463 * t605;
t359 = pkin(5) * t644 + t373;
t410 = t460 * t654 - t463 * t570;
t407 = -pkin(4) * t562 + qJD(5) - t410;
t399 = pkin(5) * t464 + t407;
t678 = -t348 * t588 + t359 * t540 - t399 * t685;
t347 = t360 * t572 - t368 * t569;
t677 = t347 * t588 + t359 * t538 + t399 * t686;
t352 = -t371 * t565 + t381 * t567;
t667 = t383 * t482 - t352;
t666 = t373 * t565 - t383 * t588;
t665 = t490 * t588 - t373;
t664 = -t373 * t567 + t382 * t588;
t663 = -t490 * t482 - t576;
t441 = -pkin(4) * t588 - qJ(5) * t482;
t660 = -0.2e1 * t608;
t659 = MDP(4) * t571;
t658 = MDP(5) * (t571 ^ 2 - t573 ^ 2);
t488 = -t544 * t566 + t630;
t470 = t488 - t652;
t489 = t544 * t568 + t532;
t471 = t489 - t651;
t422 = t470 * t570 + t471 * t654;
t615 = qJD(1) * t571;
t499 = pkin(2) * t615 + pkin(3) * t527;
t423 = t441 + t499;
t386 = -t422 * t565 + t423 * t567;
t653 = pkin(2) * t566;
t553 = t570 * t653;
t557 = pkin(2) * t568 + pkin(3);
t589 = -qJD(4) * t553 + t557 * t605;
t505 = qJD(5) + t589;
t657 = -t505 * t565 - t386;
t387 = t422 * t567 + t423 * t565;
t600 = t505 * t567 - t387;
t394 = t410 * t567 + t441 * t565;
t656 = -qJD(5) * t567 + t394;
t582 = t557 * t570 + t653 * t654;
t622 = -qJD(4) * t582 - t470 * t654 + t471 * t570;
t494 = t550 * t568 + t551 * t566;
t475 = -pkin(8) * t539 + t494;
t495 = t550 * t566 - t551 * t568;
t476 = pkin(8) * t537 + t495;
t655 = t475 * t654 - t476 * t570;
t650 = t567 * pkin(5);
t561 = t567 * pkin(9);
t643 = t439 * t567;
t529 = t537 * qJD(2);
t587 = t537 * t654 - t539 * t570;
t444 = qJD(4) * t587 - t526 * t570 + t529 * t654;
t642 = t444 * t565;
t487 = t537 * t570 + t539 * t654;
t635 = t487 * t565;
t634 = t487 * t567;
t574 = qJD(2) ^ 2;
t629 = t571 * t574;
t627 = t573 * t574;
t575 = qJD(1) ^ 2;
t626 = t573 * t575;
t473 = -t522 * t566 + t523 * t568;
t451 = -pkin(8) * t529 + t473;
t474 = t522 * t568 + t523 * t566;
t452 = -pkin(8) * t526 + t474;
t388 = qJD(4) * t655 + t451 * t570 + t452 * t654;
t445 = qJD(4) * t487 + t526 * t654 + t529 * t570;
t560 = t571 * t648;
t500 = pkin(3) * t526 + t560;
t392 = pkin(4) * t445 - qJ(5) * t444 - qJD(5) * t487 + t500;
t356 = t388 * t567 + t392 * t565;
t509 = -pkin(3) * t537 + t607;
t434 = -pkin(4) * t587 - qJ(5) * t487 + t509;
t436 = t475 * t570 + t476 * t654;
t398 = t434 * t565 + t436 * t567;
t618 = -t622 - t683;
t612 = qJD(6) * t368;
t611 = qJD(6) * t487;
t601 = pkin(1) * t660;
t355 = -t388 * t565 + t392 * t567;
t393 = -t410 * t565 + t441 * t567;
t397 = t434 * t567 - t436 * t565;
t349 = -pkin(9) * t644 + t353;
t598 = -qJD(6) * t360 - t349;
t594 = -t352 * t565 + t351;
t375 = -pkin(5) * t587 - pkin(9) * t634 + t397;
t384 = -pkin(9) * t635 + t398;
t593 = t375 * t572 - t384 * t569;
t592 = t375 * t569 + t384 * t572;
t591 = t382 * t565 - t383 * t567;
t521 = -t557 * t654 - pkin(4) + t553;
t520 = qJ(5) + t582;
t493 = t520 * t567 + t561;
t586 = qJD(6) * t493 - t657 + t680;
t492 = (-pkin(9) - t520) * t565;
t585 = -qJD(6) * t492 - t600 - t682;
t549 = qJ(5) * t567 + t561;
t584 = qJD(5) * t565 + qJD(6) * t549 + t393 + t680;
t548 = (-pkin(9) - qJ(5)) * t565;
t583 = -qJD(6) * t548 + t656 - t682;
t581 = t373 * t487 + t407 * t444 - t439 * t655;
t579 = -pkin(4) * t439 + qJ(5) * t577 - (-qJD(5) + t407) * t482;
t578 = t439 * t521 + t577 * t520 + (-t407 + t505) * t482;
t389 = qJD(4) * t436 - t451 * t654 + t452 * t570;
t558 = -pkin(4) - t650;
t501 = t521 - t650;
t443 = t538 * t487;
t442 = t540 * t487;
t402 = pkin(5) * t635 - t655;
t400 = t411 + t683;
t380 = t444 * t540 + t610 * t634 - t611 * t631;
t379 = -t444 * t538 - t540 * t611;
t365 = pkin(5) * t642 + t389;
t354 = -pkin(9) * t642 + t356;
t350 = pkin(5) * t445 - t444 * t561 + t355;
t346 = -pkin(5) * t577 - pkin(9) * t643 + t352;
t345 = t572 * t346;
t1 = [-MDP(7) * t629 + (pkin(7) * t629 + t573 * t601) * MDP(10) + (-pkin(7) * t627 + t571 * t601) * MDP(9) + (t468 * t494 + t469 * t495 + t484 * t473 + t485 * t474 + (t547 + t595) * t560) * MDP(12) + ((t350 * t572 - t354 * t569) * t478 - t593 * t577 - (-t349 * t569 + t345) * t587 + t347 * t445 + t365 * t418 + t402 * t367 + t359 * t442 + t399 * t380 + (t348 * t587 - t478 * t592) * qJD(6)) * MDP(29) + (-(t350 * t569 + t354 * t572) * t478 + t592 * t577 + (t346 * t569 + t349 * t572) * t587 - t348 * t445 - t365 * t662 + t402 * t366 - t359 * t443 + t399 * t379 + (t347 * t587 - t478 * t593) * qJD(6)) * MDP(30) + (-t355 * t466 - t356 * t464 + (-t352 * t487 - t382 * t444 - t397 * t439) * t567 + (-t353 * t487 - t383 * t444 - t398 * t439) * t565) * MDP(22) + (-t468 * t539 + t469 * t537 - t473 * t527 + t474 * t525 - t484 * t529 - t485 * t526 - t494 * t517 - t495 * t580) * MDP(11) + (t353 * t587 + t356 * t482 - t383 * t445 + t389 * t466 + t398 * t577 + t567 * t581) * MDP(21) + (-t352 * t587 - t355 * t482 + t382 * t445 + t389 * t464 - t397 * t577 + t565 * t581) * MDP(20) + 0.2e1 * t603 * t659 + t658 * t660 + (t445 * t490 - t482 * t500 - t491 * t587 - t509 * t577) * MDP(18) + (t439 * t509 + t444 * t490 + t487 * t491 - t500 * t588) * MDP(19) + (t367 * t587 - t380 * t478 - t418 * t445 + t442 * t577) * MDP(27) + (-t366 * t587 + t379 * t478 + t443 * t577 - t445 * t662) * MDP(26) + (t445 * t478 + t577 * t587) * MDP(28) + (t439 * t587 + t444 * t482 + t445 * t588 + t487 * t577) * MDP(14) + (t439 * t487 - t444 * t588) * MDP(13) + (-t366 * t442 + t367 * t443 - t379 * t418 + t380 * t662) * MDP(25) + (-t366 * t443 - t379 * t662) * MDP(24) + (t352 * t397 + t353 * t398 + t355 * t382 + t356 * t383 - t373 * t655 + t389 * t407) * MDP(23) + MDP(6) * t627 + (MDP(15) * t444 - MDP(16) * t445 - MDP(18) * t389 - MDP(19) * t388) * t562; (t373 * t521 + t382 * t657 + t600 * t383 - t622 * t407 + t594 * t520) * MDP(23) + t575 * t658 - t626 * t659 + (-t484 * t488 - t485 * t489 + (t468 * t568 + t469 * t566 - t547 * t615) * pkin(2)) * MDP(12) + (MDP(9) * t571 * t575 + MDP(10) * t626) * pkin(1) + (-t387 * t482 - t466 * t622 + t567 * t578 + t666) * MDP(21) + (t499 * t588 + (t422 - t589) * t562 + t663) * MDP(19) + (t386 * t482 - t464 * t622 + t565 * t578 + t664) * MDP(20) + (t482 * t499 + t562 * t622 + t665) * MDP(18) + ((t492 * t569 + t493 * t572) * t577 + t501 * t366 + (t569 * t586 + t572 * t585) * t478 - t618 * t662 + t678) * MDP(30) + (-(t492 * t572 - t493 * t569) * t577 + t501 * t367 + (t569 * t585 - t572 * t586) * t478 + t618 * t418 + t677) * MDP(29) + (t386 * t466 - t600 * t464 + (t466 * t505 + t667) * t565 + t679) * MDP(22) + ((t485 + t488) * t527 + (-t489 + t484) * t525 + (-t568 * t517 - t566 * t580) * pkin(2)) * MDP(11) + t684; (-t525 ^ 2 - t527 ^ 2) * MDP(11) + (t484 * t527 - t485 * t525 + t555) * MDP(12) + (-t577 - t641) * MDP(18) + (t439 + t638) * MDP(19) + (t464 * t588 - t482 * t674 - t567 * t577) * MDP(20) + (t466 * t588 - t482 * t673 + t565 * t577) * MDP(21) + ((t464 * t567 - t466 * t565) * t482 + (-t565 ^ 2 - t567 ^ 2) * t439) * MDP(22) + (t352 * t567 + t353 * t565 + t407 * t588 + t482 * t591) * MDP(23) + (t596 + t646) * MDP(29) + (-t645 - t687) * MDP(30); (t411 * t562 + t665) * MDP(18) + (t410 * t562 + t663) * MDP(19) + (t393 * t482 - t411 * t464 + t565 * t579 + t664) * MDP(20) + (-t394 * t482 - t411 * t466 + t567 * t579 + t666) * MDP(21) + (t393 * t466 + t656 * t464 + (qJD(5) * t466 + t667) * t565 + t679) * MDP(22) + (-pkin(4) * t373 + qJ(5) * t594 - qJD(5) * t591 - t382 * t393 - t383 * t394 - t407 * t411) * MDP(23) + (-(t548 * t572 - t549 * t569) * t577 + t558 * t367 - t400 * t418 + (t569 * t583 - t572 * t584) * t478 + t677) * MDP(29) + ((t548 * t569 + t549 * t572) * t577 + t558 * t366 + t400 * t662 + (t569 * t584 + t572 * t583) * t478 + t678) * MDP(30) + t684; (-t466 * t482 + t644) * MDP(20) + (t464 * t482 + t643) * MDP(21) + (-t464 ^ 2 - t466 ^ 2) * MDP(22) + (t382 * t466 + t383 * t464 + t373) * MDP(23) + (-t478 * t662 + t367) * MDP(29) + (-t478 * t670 + (-t466 * t478 + t590) * t569 + t619) * MDP(30); -t418 ^ 2 * MDP(25) + (t418 * t478 + t619) * MDP(26) - t577 * MDP(28) + (t348 * t478 + t345) * MDP(29) + (t347 * t478 + t399 * t418) * MDP(30) - (MDP(24) * t418 - MDP(25) * t662 + MDP(27) * t478 - MDP(29) * t399) * t662 + (MDP(27) * t590 - MDP(29) * t612 + MDP(30) * t598) * t572 + (t590 * MDP(26) + (qJD(6) * t464 - t643) * MDP(27) + t598 * MDP(29) + (-t346 + t612) * MDP(30)) * t569;];
tauc  = t1;
