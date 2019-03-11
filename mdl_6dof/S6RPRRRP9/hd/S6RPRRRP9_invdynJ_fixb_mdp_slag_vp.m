% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP9_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:45
% EndTime: 2019-03-09 06:28:55
% DurationCPUTime: 7.09s
% Computational Cost: add. (5151->536), mult. (10189->679), div. (0->0), fcn. (6728->10), ass. (0->232)
t526 = sin(qJ(3));
t530 = cos(qJ(3));
t563 = pkin(3) * t530 + pkin(8) * t526;
t473 = t563 * qJD(1);
t529 = cos(qJ(4));
t452 = t529 * t473;
t532 = -pkin(1) - pkin(7);
t496 = qJD(1) * t532 + qJD(2);
t662 = pkin(8) + pkin(9);
t587 = qJD(4) * t662;
t638 = t526 * t529;
t593 = pkin(9) * t638;
t525 = sin(qJ(4));
t641 = t525 * t530;
t681 = -t496 * t641 + t452 + (pkin(4) * t530 + t593) * qJD(1) + t529 * t587;
t615 = qJD(1) * t526;
t586 = t525 * t615;
t634 = t529 * t530;
t621 = t525 * t473 + t496 * t634;
t677 = pkin(9) * t586 + t525 * t587 + t621;
t611 = qJD(3) * t526;
t650 = qJDD(3) * pkin(3);
t575 = t496 * t611 - t650;
t492 = qJDD(1) * t532 + qJDD(2);
t632 = t530 * t492;
t419 = t575 - t632;
t503 = qJD(4) + t615;
t657 = g(3) * t526;
t527 = sin(qJ(1));
t531 = cos(qJ(1));
t668 = -g(1) * t527 + g(2) * t531;
t545 = -t530 * t668 - t657;
t680 = qJD(4) * pkin(8) * t503 + t419 + t545;
t610 = qJD(3) * t529;
t584 = t526 * t610;
t604 = qJD(4) * t530;
t547 = -t525 * t604 - t584;
t596 = qJDD(1) * t530;
t675 = qJD(3) * qJD(4) + t596;
t401 = qJD(1) * t547 + t525 * qJDD(3) + t529 * t675;
t614 = qJD(1) * t530;
t580 = t525 * t614;
t464 = -t580 + t610;
t612 = qJD(3) * t525;
t465 = t529 * t614 + t612;
t524 = sin(qJ(5));
t528 = cos(qJ(5));
t581 = t529 * t604;
t588 = qJD(1) * t581 + t525 * t675;
t554 = t529 * qJDD(3) - t588;
t585 = t525 * t611;
t541 = qJD(1) * t585 + t554;
t602 = qJD(5) * t528;
t603 = qJD(5) * t524;
t359 = -t528 * t401 - t464 * t602 + t465 * t603 - t524 * t541;
t558 = t464 * t524 + t528 * t465;
t360 = qJD(5) * t558 + t401 * t524 - t528 * t541;
t408 = -t528 * t464 + t465 * t524;
t405 = t408 ^ 2;
t600 = qJD(1) * qJD(3);
t578 = t530 * t600;
t597 = qJDD(1) * t526;
t461 = qJDD(4) + t578 + t597;
t456 = qJDD(5) + t461;
t494 = qJD(5) + t503;
t663 = t558 ^ 2;
t679 = t456 * MDP(25) + (t494 * t558 - t360) * MDP(24) + t408 * MDP(21) * t558 + (t408 * t494 - t359) * MDP(23) + (-t405 + t663) * MDP(22);
t678 = qJ(6) * t408;
t594 = qJD(4) + qJD(5);
t605 = qJD(4) * t529;
t635 = t528 * t529;
t643 = t524 * t525;
t623 = t524 * t586 - t528 * t605 - t529 * t602 + t594 * t643 - t615 * t635;
t468 = t524 * t529 + t525 * t528;
t416 = t594 * t468;
t442 = t468 * qJD(1);
t622 = t526 * t442 + t416;
t676 = t581 - t585;
t645 = t496 * t530;
t653 = qJD(3) * pkin(3);
t454 = -t645 - t653;
t418 = -pkin(4) * t464 + t454;
t523 = qJ(4) + qJ(5);
t513 = sin(t523);
t514 = cos(t523);
t639 = t526 * t527;
t427 = t513 * t531 + t514 * t639;
t637 = t526 * t531;
t429 = -t513 * t527 + t514 * t637;
t479 = pkin(3) * t526 - pkin(8) * t530 + qJ(2);
t443 = t479 * qJD(1);
t476 = t526 * t496;
t453 = qJD(3) * pkin(8) + t476;
t398 = t443 * t525 + t453 * t529;
t462 = qJD(3) * t563 + qJD(2);
t406 = qJD(1) * t462 + qJDD(1) * t479;
t403 = t529 * t406;
t609 = qJD(3) * t530;
t420 = qJDD(3) * pkin(8) + t492 * t526 + t496 * t609;
t349 = pkin(4) * t461 - pkin(9) * t401 - qJD(4) * t398 - t525 * t420 + t403;
t590 = -t525 * t406 - t529 * t420 - t443 * t605;
t607 = qJD(4) * t525;
t550 = -t453 * t607 - t590;
t353 = pkin(9) * t541 + t550;
t397 = t529 * t443 - t453 * t525;
t384 = -pkin(9) * t465 + t397;
t379 = pkin(4) * t503 + t384;
t385 = pkin(9) * t464 + t398;
t565 = -t524 * t349 - t528 * t353 - t379 * t602 + t385 * t603;
t656 = g(3) * t530;
t674 = g(1) * t427 - g(2) * t429 + t408 * t418 + t514 * t656 + t565;
t672 = qJ(6) * t558;
t372 = pkin(5) * t408 + qJD(6) + t418;
t671 = t372 * t558;
t467 = -t635 + t643;
t433 = t467 * t526;
t431 = t468 * t526;
t460 = t529 * t479;
t640 = t525 * t532;
t577 = pkin(4) - t640;
t404 = -pkin(9) * t634 + t526 * t577 + t460;
t493 = t532 * t638;
t619 = t525 * t479 + t493;
t417 = -pkin(9) * t641 + t619;
t624 = t524 * t404 + t528 * t417;
t669 = t681 * t528;
t487 = t662 * t525;
t488 = t662 * t529;
t620 = -t524 * t487 + t528 * t488;
t667 = t487 * t602 + t488 * t603 + t524 * t681 + t677 * t528;
t426 = -t513 * t639 + t514 * t531;
t428 = t513 * t637 + t514 * t527;
t666 = -g(1) * t426 - g(2) * t428 + t513 * t656;
t383 = t528 * t385;
t357 = t379 * t524 + t383;
t574 = t528 * t349 - t524 * t353;
t539 = -qJD(5) * t357 + t574;
t665 = -t418 * t558 + t539 + t666;
t661 = pkin(8) * t461;
t517 = t529 * pkin(4);
t655 = pkin(3) + t517;
t477 = pkin(4) * t525 + pkin(5) * t513;
t654 = pkin(7) + t477;
t652 = pkin(1) * qJDD(1);
t534 = qJD(1) ^ 2;
t651 = qJ(2) * t534;
t649 = t401 * t525;
t648 = t461 * t529;
t647 = t464 * t503;
t646 = t465 * t529;
t644 = t503 * t525;
t381 = t524 * t385;
t642 = t525 * t461;
t636 = t527 * t529;
t633 = t529 * t531;
t356 = t528 * t379 - t381;
t350 = t356 - t672;
t348 = pkin(5) * t494 + t350;
t631 = -t350 + t348;
t630 = -qJ(6) * t622 - qJD(6) * t467 - t667;
t629 = -pkin(5) * t614 + qJ(6) * t623 - t620 * qJD(5) - qJD(6) * t468 + t524 * t677 - t669;
t434 = t467 * t530;
t628 = -qJD(3) * t434 - t431 * t594 - t442;
t627 = t467 * qJD(1) + t433 * t594 - t468 * t609;
t626 = t528 * t384 - t381;
t478 = pkin(5) * t514 + t517;
t522 = t530 ^ 2;
t618 = t526 ^ 2 - t522;
t533 = qJD(3) ^ 2;
t617 = -t533 - t534;
t613 = qJD(3) * t465;
t608 = qJD(3) * t532;
t606 = qJD(4) * t526;
t601 = t454 * qJD(4);
t598 = qJDD(1) * qJ(2);
t595 = qJDD(3) * t526;
t591 = 0.2e1 * qJD(1) * qJD(2);
t583 = t529 * t609;
t589 = t525 * t462 + t479 * t605 + t532 * t583;
t579 = g(2) * (t531 * pkin(1) + t527 * qJ(2));
t438 = t529 * t462;
t368 = t438 + (-t493 + (pkin(9) * t530 - t479) * t525) * qJD(4) + (t530 * t577 + t593) * qJD(3);
t370 = -pkin(9) * t676 - t606 * t640 + t589;
t573 = t528 * t368 - t370 * t524;
t572 = -t384 * t524 - t383;
t570 = t528 * t404 - t417 * t524;
t569 = t503 * t532 + t453;
t568 = -t528 * t487 - t488 * t524;
t463 = pkin(4) * t641 - t530 * t532;
t567 = -qJD(4) * t443 - t420;
t566 = qJD(1) + t606;
t564 = -t476 + (t586 + t607) * pkin(4);
t562 = g(1) * t531 + g(2) * t527;
t560 = qJDD(2) + t668;
t472 = pkin(3) + t478;
t520 = -qJ(6) - t662;
t556 = t472 * t526 + t520 * t530;
t555 = -t492 - t668;
t553 = t503 * t605 + t642;
t552 = -t503 * t607 + t648;
t421 = pkin(4) * t676 + t526 * t608;
t551 = 0.2e1 * qJ(2) * t600 + qJDD(3) * t532;
t549 = t524 * t368 + t528 * t370 + t404 * t602 - t417 * t603;
t546 = t555 + t651;
t542 = -t562 + t591 + 0.2e1 * t598;
t538 = -t532 * t533 + t542;
t537 = -pkin(4) * t541 + t575;
t535 = t360 * pkin(5) + qJDD(6) + t537;
t516 = t531 * qJ(2);
t512 = qJDD(3) * t530;
t509 = pkin(4) * t528 + pkin(5);
t447 = -t525 * t527 + t526 * t633;
t446 = t525 * t637 + t636;
t445 = t525 * t531 + t526 * t636;
t444 = -t525 * t639 + t633;
t432 = t468 * t530;
t391 = -qJ(6) * t467 + t620;
t390 = -qJ(6) * t468 + t568;
t378 = -t603 * t641 + (t594 * t634 - t585) * t528 + t547 * t524;
t376 = t416 * t530 - t524 * t585 + t528 * t584;
t371 = t537 - t632;
t365 = -qJ(6) * t432 + t624;
t364 = pkin(5) * t526 + qJ(6) * t434 + t570;
t355 = t626 - t672;
t354 = t572 + t678;
t351 = t357 - t678;
t345 = t535 - t632;
t344 = -qJ(6) * t378 - qJD(6) * t432 + t549;
t343 = pkin(5) * t609 + qJ(6) * t376 - qJD(5) * t624 + qJD(6) * t434 + t573;
t342 = -qJ(6) * t360 - qJD(6) * t408 - t565;
t341 = pkin(5) * t456 + qJ(6) * t359 - qJD(6) * t558 + t539;
t1 = [t562 * MDP(3) + (-t360 * t526 - t378 * t494 - t408 * t609 - t432 * t456) * MDP(24) + (t461 * t526 + t503 * t609) * MDP(18) + (t456 * t526 + t494 * t609) * MDP(25) + ((-t503 * t610 + t401) * t526 + (t552 + t613) * t530) * MDP(16) + (((t503 + t615) * t612 + t554) * t526 + (qJD(3) * t464 - t553) * t530) * MDP(17) + 0.2e1 * (-t526 * t596 + t600 * t618) * MDP(8) + (-t589 * t503 - t619 * t461 + g(1) * t446 - g(2) * t444 + (t569 * t607 + (-t454 * t529 + t465 * t532) * qJD(3) + t590) * t526 + (-qJD(3) * t398 - t401 * t532 + t419 * t529 - t525 * t601) * t530) * MDP(20) + (t401 * t634 + t465 * t547) * MDP(14) + ((-t479 * t607 + t438) * t503 + t460 * t461 - g(1) * t447 - g(2) * t445 + (-t464 * t608 + t403 - t569 * t605 + (-qJD(3) * t454 - t461 * t532 + t567) * t525) * t526 + (t532 * t554 + t419 * t525 + t529 * t601 + (t397 + (-t503 + t615) * t640) * qJD(3)) * t530) * MDP(19) + (-(qJDD(2) - t652) * pkin(1) - g(1) * (-pkin(1) * t527 + t516) - t579 + (t591 + t598) * qJ(2)) * MDP(6) + (t560 - 0.2e1 * t652) * MDP(4) + (t342 * t365 + t351 * t344 + t341 * t364 + t348 * t343 + t345 * (pkin(5) * t432 + t463) + t372 * (pkin(5) * t378 + t421) - g(1) * t516 - t579 + (-g(1) * t556 - g(2) * t654) * t531 + (-g(1) * (-pkin(1) - t654) - g(2) * t556) * t527) * MDP(29) + (-t530 * t533 - t595) * MDP(10) + t542 * MDP(5) + ((-t464 * t529 + t465 * t525) * t611 + (t529 * t541 - t649 + (-t464 * t525 - t646) * qJD(4)) * t530) * MDP(15) + (-t526 * t533 + t512) * MDP(9) + (t341 * t434 - t342 * t432 - t343 * t558 - t344 * t408 + t348 * t376 - t351 * t378 + t359 * t364 - t360 * t365 + t530 * t562) * MDP(28) + (-t359 * t526 - t376 * t494 - t434 * t456 + t558 * t609) * MDP(23) + (g(1) * t428 - g(2) * t426 - t357 * t609 - t463 * t359 - t371 * t434 - t418 * t376 + t421 * t558 - t456 * t624 - t494 * t549 + t526 * t565) * MDP(27) + (t359 * t434 - t376 * t558) * MDP(21) + (t359 * t432 + t360 * t434 + t376 * t408 - t378 * t558) * MDP(22) + (qJDD(1) * t522 - 0.2e1 * t526 * t578) * MDP(7) + (t573 * t494 + t570 * t456 + t574 * t526 + t356 * t609 + t421 * t408 + t463 * t360 + t371 * t432 + t418 * t378 - g(1) * t429 - g(2) * t427 + (-t357 * t526 - t494 * t624) * qJD(5)) * MDP(26) - t668 * MDP(2) + (t526 * t538 + t530 * t551) * MDP(12) + (-t526 * t551 + t530 * t538) * MDP(13) + qJDD(1) * MDP(1); qJDD(1) * MDP(4) - t534 * MDP(5) + (t560 - t651 - t652) * MDP(6) + (t526 * t617 + t512) * MDP(12) + (t530 * t617 - t595) * MDP(13) + (t530 * t554 + (-t642 + (-t464 + t580) * qJD(3)) * t526 + (-t525 * t609 - t529 * t566) * t503) * MDP(19) + (-t401 * t530 + (t613 - t648) * t526 + (t525 * t566 - t583) * t503) * MDP(20) + (-t360 * t530 + t408 * t611 - t431 * t456 + t494 * t627) * MDP(26) + (t359 * t530 + t433 * t456 - t494 * t628 + t558 * t611) * MDP(27) + (-t359 * t431 + t360 * t433 - t408 * t628 - t558 * t627) * MDP(28) + (-t341 * t431 - t342 * t433 - t345 * t530 + t348 * t627 + t351 * t628 + t372 * t611 + t668) * MDP(29); MDP(9) * t596 - MDP(10) * t597 + qJDD(3) * MDP(11) + (-t530 * t546 + t657) * MDP(12) + (t526 * t546 + t656) * MDP(13) + (t503 * t646 + t649) * MDP(14) + ((t401 + t647) * t529 + (-t465 * qJD(4) + (-t465 + t612) * t615 + t554) * t525) * MDP(15) + ((-t465 * t530 + t503 * t638) * qJD(1) + t553) * MDP(16) + ((-t464 * t530 - t526 * t644) * qJD(1) + t552) * MDP(17) + (-pkin(3) * t588 - t452 * t503 + t464 * t476 + (t503 * t645 - t661 + t601 + (t454 + t653) * t615) * t525 + (t650 - t680) * t529) * MDP(19) + (-pkin(3) * t401 + t621 * t503 - t465 * t476 + (t454 * t503 - t661) * t529 + t680 * t525) * MDP(20) + (-t359 * t468 - t558 * t623) * MDP(21) + (t359 * t467 - t360 * t468 + t408 * t623 - t558 * t622) * MDP(22) + (t456 * t468 - t494 * t623) * MDP(23) + (-t456 * t467 - t494 * t622) * MDP(24) + (t568 * t456 - t655 * t360 + t371 * t467 + (-t488 * t602 + (qJD(5) * t487 + t677) * t524 - t669) * t494 + t622 * t418 + t564 * t408 - t545 * t514) * MDP(26) + (t359 * t655 + t371 * t468 - t623 * t418 - t620 * t456 + t494 * t667 + t545 * t513 + t564 * t558) * MDP(27) + (-t341 * t468 - t342 * t467 + t348 * t623 - t351 * t622 + t359 * t390 - t360 * t391 - t408 * t630 + t526 * t668 - t558 * t629 - t656) * MDP(28) + (t342 * t391 + t341 * t390 + t345 * (pkin(5) * t467 - t655) + g(3) * t556 + (pkin(4) * t644 + pkin(5) * t622 - t476) * t372 + t630 * t351 + t629 * t348 + t668 * (t472 * t530 - t520 * t526)) * MDP(29) + (-t503 * MDP(18) - MDP(19) * t397 + t398 * MDP(20) - MDP(23) * t558 + t408 * MDP(24) - t494 * MDP(25) - t356 * MDP(26) + t357 * MDP(27)) * t614 + (MDP(7) * t526 * t530 - MDP(8) * t618) * t534; -t465 * t464 * MDP(14) + (-t464 ^ 2 + t465 ^ 2) * MDP(15) + (t401 - t647) * MDP(16) + (t465 * t503 + t541) * MDP(17) + t461 * MDP(18) + (-t453 * t605 - g(1) * t444 - g(2) * t446 + t398 * t503 - t454 * t465 + t403 + (t567 + t656) * t525) * MDP(19) + (g(1) * t445 - g(2) * t447 + g(3) * t634 + t397 * t503 - t454 * t464 - t550) * MDP(20) + (-t572 * t494 + (-t408 * t465 + t456 * t528 - t494 * t603) * pkin(4) + t665) * MDP(26) + (t626 * t494 + (-t456 * t524 - t465 * t558 - t494 * t602) * pkin(4) + t674) * MDP(27) + (-t348 * t408 + t351 * t558 + t354 * t558 + t355 * t408 + t359 * t509 + (-t360 * t524 + (-t408 * t528 + t524 * t558) * qJD(5)) * pkin(4)) * MDP(28) + (t341 * t509 - t351 * t355 - t348 * t354 - pkin(5) * t671 - g(1) * (-t477 * t639 + t478 * t531) - g(2) * (t477 * t637 + t478 * t527) + t477 * t656 + (t342 * t524 - t372 * t465 + (-t348 * t524 + t351 * t528) * qJD(5)) * pkin(4)) * MDP(29) + t679; (t357 * t494 + t665) * MDP(26) + (t356 * t494 + t674) * MDP(27) + (pkin(5) * t359 - t408 * t631) * MDP(28) + (t631 * t351 + (t341 + t666 - t671) * pkin(5)) * MDP(29) + t679; (-t405 - t663) * MDP(28) + (t348 * t558 + t351 * t408 + t530 * t555 + t535 - t657) * MDP(29);];
tau  = t1;
