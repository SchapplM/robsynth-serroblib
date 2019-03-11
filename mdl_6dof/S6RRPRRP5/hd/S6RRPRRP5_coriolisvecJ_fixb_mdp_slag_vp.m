% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:41
% EndTime: 2019-03-09 12:05:56
% DurationCPUTime: 8.84s
% Computational Cost: add. (9396->521), mult. (27937->703), div. (0->0), fcn. (22317->10), ass. (0->229)
t521 = sin(pkin(6));
t522 = cos(pkin(11));
t529 = cos(qJ(2));
t622 = t522 * t529;
t584 = t521 * t622;
t505 = qJD(1) * t584;
t520 = sin(pkin(11));
t526 = sin(qJ(2));
t601 = qJD(1) * t521;
t583 = t526 * t601;
t556 = t520 * t583 - t505;
t659 = qJD(4) + t556;
t528 = cos(qJ(4));
t649 = qJD(4) * t659;
t660 = t556 * t659 + t649;
t665 = t528 * t660;
t523 = cos(pkin(6));
t645 = pkin(1) * t529;
t589 = t523 * t645;
t512 = qJD(1) * t589;
t642 = pkin(8) + qJ(3);
t579 = t642 * t526;
t568 = t521 * t579;
t475 = -qJD(1) * t568 + t512;
t646 = pkin(1) * t526;
t587 = t523 * t646;
t624 = t521 * t529;
t483 = t624 * t642 + t587;
t476 = t483 * qJD(1);
t623 = t522 * t476;
t432 = t475 * t520 + t623;
t525 = sin(qJ(4));
t664 = t432 - t659 * (pkin(4) * t525 - pkin(10) * t528);
t559 = t520 * t529 + t522 * t526;
t488 = t559 * t601;
t600 = qJD(1) * t523;
t571 = qJD(2) + t600;
t458 = t525 * t488 - t528 * t571;
t599 = qJD(2) * t521;
t582 = t526 * t599;
t566 = qJD(1) * t582;
t482 = qJD(2) * t505 - t520 * t566;
t616 = t528 * t482;
t532 = -qJD(4) * t458 + t616;
t663 = qJD(5) * t659 + t532;
t524 = sin(qJ(5));
t527 = cos(qJ(5));
t617 = t527 * t528;
t449 = t488 * t524 - t556 * t617;
t596 = qJD(5) * t524;
t597 = qJD(4) * t528;
t662 = t525 * t596 - t527 * t597 + t449;
t641 = -qJ(6) - pkin(10);
t661 = -qJ(6) * t458 + qJD(5) * t641;
t658 = pkin(2) * t582;
t657 = MDP(5) * (t526 ^ 2 - t529 ^ 2);
t656 = MDP(6) * t529;
t460 = t528 * t488 + t525 * t571;
t416 = t460 * t524 - t527 * t659;
t655 = t416 * t659;
t456 = qJD(5) + t458;
t654 = t662 * t456;
t625 = t521 * t526;
t490 = t520 * t625 - t584;
t474 = (pkin(2) + t645) * t523 - t568;
t441 = t520 * t474 + t522 * t483;
t431 = pkin(9) * t523 + t441;
t491 = t559 * t521;
t565 = (-pkin(2) * t529 - pkin(1)) * t521;
t447 = pkin(3) * t490 - pkin(9) * t491 + t565;
t609 = t528 * t431 + t525 * t447;
t381 = pkin(10) * t490 + t609;
t440 = t474 * t522 - t520 * t483;
t430 = -pkin(3) * t523 - t440;
t468 = t491 * t525 - t523 * t528;
t469 = t491 * t528 + t523 * t525;
t385 = pkin(4) * t468 - pkin(10) * t469 + t430;
t612 = t527 * t381 + t524 * t385;
t515 = pkin(2) * t520 + pkin(9);
t593 = t515 * qJD(4);
t620 = t524 * t525;
t653 = -t527 * t664 + t593 * t620;
t517 = t521 ^ 2;
t552 = -pkin(8) * t624 - t587;
t652 = -t517 * t646 + t552 * t523;
t630 = t556 * t528;
t651 = t597 + t630;
t465 = t520 * t476;
t433 = t475 * t522 - t465;
t443 = pkin(2) * t583 + pkin(3) * t488 + pkin(9) * t556;
t608 = t528 * t433 + t525 * t443;
t376 = pkin(10) * t488 + t608;
t516 = -pkin(2) * t522 - pkin(3);
t500 = -pkin(4) * t528 - pkin(10) * t525 + t516;
t595 = qJD(5) * t527;
t650 = t527 * t376 - t500 * t595 + t524 * t664;
t548 = qJD(2) * t559;
t486 = t521 * t548;
t536 = qJD(1) * t486;
t366 = t460 * t596 - t524 * t536 - t527 * t663;
t418 = t527 * t460 + t524 * t659;
t448 = -t527 * t488 - t524 * t630;
t563 = -t524 * t597 + t448;
t544 = t525 * t595 - t563;
t648 = t366 * t620 - t418 * t544;
t647 = t418 ^ 2;
t530 = qJD(1) ^ 2;
t644 = pkin(4) * t488;
t643 = pkin(5) * t524;
t508 = t512 * qJD(2);
t535 = (-qJD(2) * t579 + qJD(3) * t529) * t521;
t453 = qJD(1) * t535 + t508;
t463 = -qJD(2) * t483 - qJD(3) * t625;
t531 = qJD(1) * t463;
t393 = t522 * t453 + t520 * t531;
t461 = qJD(2) * pkin(2) + t512 + (t523 * pkin(2) - t568) * qJD(1);
t411 = t520 * t461 + t623;
t403 = pkin(9) * t571 + t411;
t557 = qJD(1) * t565;
t497 = qJD(3) + t557;
t428 = pkin(3) * t556 - pkin(9) * t488 + t497;
t507 = pkin(2) * t566;
t537 = qJD(2) * t488;
t429 = pkin(3) * t537 - t482 * pkin(9) + t507;
t598 = qJD(4) * t525;
t569 = t525 * t393 + t403 * t597 + t428 * t598 - t528 * t429;
t348 = -pkin(4) * t536 + t569;
t639 = t348 * t524;
t638 = t348 * t527;
t367 = t460 * t595 + t524 * t663 - t527 * t536;
t637 = t367 * t528;
t619 = t525 * t482;
t415 = t460 * qJD(4) + t619;
t636 = t415 * t527;
t635 = t416 * t456;
t634 = t418 * t456;
t633 = t418 * t556;
t572 = t456 * t527;
t632 = t460 * t488;
t631 = t556 * t525;
t629 = t488 * t458;
t628 = t515 * t528;
t627 = t517 * t530;
t626 = t520 * t526;
t621 = t524 * t415;
t618 = t525 * t527;
t374 = t528 * t403 + t525 * t428;
t364 = pkin(10) * t659 + t374;
t410 = t522 * t461 - t465;
t402 = -pkin(3) * t571 - t410;
t370 = t458 * pkin(4) - t460 * pkin(10) + t402;
t345 = -t364 * t524 + t527 * t370;
t340 = -qJ(6) * t418 + t345;
t339 = pkin(5) * t456 + t340;
t615 = t339 - t340;
t373 = -t525 * t403 + t528 * t428;
t406 = pkin(4) * t460 + pkin(10) * t458;
t614 = t527 * t373 + t524 * t406;
t503 = t515 * t617;
t594 = qJD(6) * t527;
t611 = pkin(5) * t631 + qJ(6) * t449 + t376 * t524 - t525 * t594 + (pkin(5) * t525 - qJ(6) * t617) * qJD(4) + (-t503 + (qJ(6) * t525 - t500) * t524) * qJD(5) + t653;
t610 = qJ(6) * t448 + (-qJ(6) * qJD(5) - t593) * t618 + (-qJD(6) * t525 + (-qJ(6) * qJD(4) - qJD(5) * t515) * t528) * t524 - t650;
t606 = t524 * t661 + t594 - t614;
t400 = t527 * t406;
t605 = -pkin(5) * t460 - t400 + t661 * t527 + (-qJD(6) + t373) * t524;
t603 = t524 * t500 + t503;
t591 = qJD(4) - t505;
t590 = qJD(1) * qJD(2);
t585 = t529 * t627;
t580 = t521 * t523 * t530;
t578 = t517 * t590;
t576 = t366 * t528 + t418 * t598;
t575 = -t381 * t524 + t527 * t385;
t574 = -t525 * t431 + t447 * t528;
t426 = t525 * t433;
t573 = t443 * t528 - t426;
t392 = t520 * t453 - t522 * t531;
t513 = qJD(2) * t589;
t462 = t513 + t535;
t404 = t462 * t520 - t522 * t463;
t567 = t529 * t578;
t346 = t364 * t527 + t370 * t524;
t341 = -qJ(6) * t416 + t346;
t561 = -t339 * t527 - t341 * t524;
t560 = t339 * t524 - t341 * t527;
t446 = t469 * t527 + t490 * t524;
t445 = t469 * t524 - t490 * t527;
t405 = t462 * t522 + t463 * t520;
t487 = (t622 - t626) * t599;
t444 = pkin(3) * t486 - pkin(9) * t487 + t658;
t558 = -t525 * t405 - t431 * t597 + t444 * t528 - t447 * t598;
t380 = -pkin(4) * t490 - t574;
t555 = t563 * t456;
t554 = -t456 * t595 - t621;
t363 = -pkin(4) * t659 - t373;
t553 = -pkin(10) * t415 + t363 * t456;
t551 = -t528 * t393 + t403 * t598 - t428 * t597 - t525 * t429;
t550 = t528 * t405 - t431 * t598 + t525 * t444 + t447 * t597;
t347 = pkin(10) * t536 - t551;
t359 = t415 * pkin(4) - pkin(10) * t532 + t392;
t336 = t527 * t347 + t524 * t359 - t364 * t596 + t370 * t595;
t353 = pkin(10) * t486 + t550;
t438 = qJD(4) * t469 + t487 * t525;
t439 = -qJD(4) * t468 + t487 * t528;
t362 = pkin(4) * t438 - pkin(10) * t439 + t404;
t549 = t527 * t353 + t524 * t362 - t381 * t596 + t385 * t595;
t354 = -pkin(4) * t486 - t558;
t538 = -t367 * t618 + t416 * t662;
t337 = -qJD(5) * t346 - t347 * t524 + t527 * t359;
t534 = -qJD(5) * t612 - t353 * t524 + t527 * t362;
t533 = -t525 * t660 + t528 * t536;
t338 = t367 * pkin(5) + t348;
t510 = t641 * t527;
t509 = t641 * t524;
t494 = t527 * t500;
t457 = -qJ(6) * t620 + t603;
t452 = -qJ(6) * t618 + t494 + (-t515 * t524 - pkin(5)) * t528;
t414 = t416 ^ 2;
t379 = -qJD(5) * t445 + t439 * t527 + t486 * t524;
t378 = qJD(5) * t446 + t439 * t524 - t486 * t527;
t375 = -t573 - t644;
t356 = t416 * pkin(5) + qJD(6) + t363;
t349 = -qJ(6) * t445 + t612;
t342 = pkin(5) * t468 - qJ(6) * t446 + t575;
t335 = -qJ(6) * t378 - qJD(6) * t445 + t549;
t334 = pkin(5) * t438 - qJ(6) * t379 - qJD(6) * t446 + t534;
t333 = -qJ(6) * t367 - qJD(6) * t416 + t336;
t332 = pkin(5) * t415 + qJ(6) * t366 - qJD(6) * t418 + t337;
t1 = [-0.2e1 * t578 * t657 + (t552 * qJD(2) ^ 2 + 0.2e1 * t590 * t652) * MDP(9) + (-0.2e1 * pkin(1) * t567 - (-pkin(8) * t582 + t513) * t571 - (-pkin(8) * t566 + t508) * t523) * MDP(10) + (t392 * t491 - t393 * t490 + t404 * t488 - t405 * t556 - t410 * t487 - t411 * t486 - t440 * t482 - t441 * t537) * MDP(11) + (-t392 * t440 + t393 * t441 - t410 * t404 + t411 * t405 + (t497 + t557) * t658) * MDP(12) + (t460 * t439 + t469 * t532) * MDP(13) + (-t469 * t415 - t460 * t438 - t439 * t458 - t468 * t532) * MDP(14) + (t439 * t659 + t460 * t486 + t469 * t536 + t490 * t532) * MDP(15) + (-t438 * t591 - t415 * t490 - t458 * t486 + (-t438 * t626 - t468 * t548) * t601) * MDP(16) + (t591 * t486 + (t486 * t626 + t490 * t548) * t601) * MDP(17) + (t373 * t486 + t392 * t468 + t402 * t438 + t404 * t458 + t430 * t415 - t490 * t569 + t536 * t574 + t558 * t659) * MDP(18) + (-t374 * t486 + t392 * t469 + t402 * t439 + t404 * t460 + t430 * t532 + t490 * t551 - t537 * t609 - t550 * t659) * MDP(19) + (-t366 * t446 + t379 * t418) * MDP(20) + (t366 * t445 - t367 * t446 - t378 * t418 - t379 * t416) * MDP(21) + (-t366 * t468 + t379 * t456 + t415 * t446 + t418 * t438) * MDP(22) + (-t367 * t468 - t378 * t456 - t415 * t445 - t416 * t438) * MDP(23) + (t415 * t468 + t438 * t456) * MDP(24) + (t337 * t468 + t345 * t438 + t348 * t445 + t354 * t416 + t363 * t378 + t380 * t367 + t415 * t575 + t456 * t534) * MDP(25) + (-t336 * t468 - t346 * t438 + t348 * t446 + t354 * t418 + t363 * t379 - t380 * t366 - t415 * t612 - t456 * t549) * MDP(26) + (-t332 * t446 - t333 * t445 - t334 * t418 - t335 * t416 - t339 * t379 - t341 * t378 + t342 * t366 - t349 * t367) * MDP(27) + (t333 * t349 + t341 * t335 + t332 * t342 + t339 * t334 + t338 * (pkin(5) * t445 + t380) + t356 * (pkin(5) * t378 + t354)) * MDP(28) + 0.2e1 * t526 * MDP(4) * t567 + (-MDP(7) * t582 + t599 * t656) * (qJD(2) + 0.2e1 * t600); t627 * t657 - t580 * t656 + (pkin(1) * t585 + (-pkin(8) * t583 + t512) * t600) * MDP(10) + ((t411 - t432) * t488 - (-t433 + t410) * t556 + (-t522 * t482 - t520 * t537) * pkin(2)) * MDP(11) + (t410 * t432 - t411 * t433 + (-t392 * t522 + t393 * t520 - t497 * t583) * pkin(2)) * MDP(12) + (-qJD(4) * t525 ^ 2 * t488 + ((qJD(4) * t571 + t482) * t525 + t659 * t460) * t528) * MDP(13) + (-t525 * t415 + t528 * t532 + (-t598 - t631) * t460 - t651 * t458) * MDP(14) + (t525 * t536 - t632 + t665) * MDP(15) + (t533 + t629) * MDP(16) + (-t649 * t628 + t516 * t415 - t392 * t528 - t573 * t659 - t373 * t488 - t432 * t458 + (t402 * t659 - t515 * t536) * t525) * MDP(18) + (t374 * t488 - t432 * t460 + t516 * t532 - t537 * t628 + t659 * t608 + (t515 * t649 + t392) * t525 + t651 * t402) * MDP(19) + (-t366 * t618 - t418 * t662) * MDP(20) + (t538 + t648) * MDP(21) + ((t633 + t636) * t525 - t654 + t576) * MDP(22) + (t637 + t555 + (t554 - t655) * t525) * MDP(23) + (t456 * t525 * t659 - t415 * t528) * MDP(24) + (-t363 * t448 - t375 * t416 + t494 * t415 + ((-qJD(5) * t500 + t376) * t524 + t653) * t456 + (t363 * t524 * qJD(4) - t337 + (qJD(4) * t416 + t554) * t515) * t528 + (t345 * t659 + t363 * t595 + t515 * t367 + t639) * t525) * MDP(25) + (-t603 * t415 - t375 * t418 - t363 * t449 + t650 * t456 + (t515 * t456 * t596 + t336 + (t363 * t527 + t418 * t515) * qJD(4)) * t528 + (-t363 * t596 - t346 * t556 + t638 - t515 * t366 + (t515 * t572 - t346) * qJD(4)) * t525) * MDP(26) + (t339 * t449 + t341 * t448 + t366 * t452 - t367 * t457 - t611 * t418 - t610 * t416 + t561 * t597 + (qJD(5) * t560 - t332 * t527 - t333 * t524) * t525) * MDP(27) + (t332 * t452 + t333 * t457 + t338 * (t515 + t643) * t525 + (t644 - t426 + (t443 + t593) * t528 + t544 * pkin(5)) * t356 + t610 * t341 + t611 * t339) * MDP(28) - t659 * t488 * MDP(17) - t652 * MDP(9) * t530 + (-MDP(4) * t585 + MDP(7) * t580) * t526; (-t488 ^ 2 - t556 ^ 2) * MDP(11) + (t410 * t488 + t411 * t556 + t507) * MDP(12) + (t533 - t629) * MDP(18) + (-t525 * t537 - t632 - t665) * MDP(19) + (-t637 + t555 + (t554 + t655) * t525) * MDP(25) + ((t633 - t636) * t525 + t654 + t576) * MDP(26) + (t538 - t648) * MDP(27) + (t339 * t448 - t341 * t449 + (-qJD(4) * t560 - t338) * t528 + (qJD(5) * t561 - t332 * t524 + t333 * t527 + t356 * t659) * t525) * MDP(28); -t458 ^ 2 * MDP(14) + (t458 * t556 + t616) * MDP(15) - t619 * MDP(16) + MDP(17) * t536 + (t374 * t659 - t569) * MDP(18) + (t373 * t659 + t402 * t458 + t551) * MDP(19) + (-t366 * t524 + t418 * t572) * MDP(20) + ((-t366 - t635) * t527 + (-t367 - t634) * t524) * MDP(21) + (t456 * t572 + t621) * MDP(22) + (-t456 ^ 2 * t524 + t636) * MDP(23) + (-pkin(4) * t367 - t638 - t374 * t416 + (-pkin(10) * t595 - t400) * t456 + (t373 * t456 + t553) * t524) * MDP(25) + (pkin(4) * t366 + t639 - t374 * t418 + (pkin(10) * t596 + t614) * t456 + t553 * t527) * MDP(26) + (t366 * t509 + t367 * t510 - t605 * t418 - t606 * t416 + (-t339 * t456 + t333) * t527 + (-t341 * t456 - t332) * t524) * MDP(27) + (-t333 * t510 + t332 * t509 + t338 * (-pkin(5) * t527 - pkin(4)) + (t456 * t643 - t374) * t356 + t606 * t341 + t605 * t339) * MDP(28) + (MDP(13) * t458 + MDP(14) * t460 + MDP(16) * t556 - t402 * MDP(18) - t418 * MDP(22) + t416 * MDP(23) - t456 * MDP(24) - t345 * MDP(25) + t346 * MDP(26)) * t460; t418 * t416 * MDP(20) + (-t414 + t647) * MDP(21) + (-t366 + t635) * MDP(22) + (-t367 + t634) * MDP(23) + t415 * MDP(24) + (t346 * t456 - t363 * t418 + t337) * MDP(25) + (t345 * t456 + t363 * t416 - t336) * MDP(26) + (pkin(5) * t366 - t416 * t615) * MDP(27) + (t615 * t341 + (-t356 * t418 + t332) * pkin(5)) * MDP(28); (-t414 - t647) * MDP(27) + (t339 * t418 + t341 * t416 + t338) * MDP(28);];
tauc  = t1;
