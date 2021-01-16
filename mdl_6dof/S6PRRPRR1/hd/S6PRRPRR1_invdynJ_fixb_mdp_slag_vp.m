% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:44
% EndTime: 2021-01-16 03:29:06
% DurationCPUTime: 10.32s
% Computational Cost: add. (5801->537), mult. (13818->714), div. (0->0), fcn. (11451->18), ass. (0->229)
t553 = cos(qJ(6));
t620 = qJD(6) * t553;
t543 = sin(pkin(12));
t546 = cos(pkin(12));
t555 = cos(qJ(3));
t623 = qJD(2) * t555;
t608 = t546 * t623;
t551 = sin(qJ(3));
t624 = qJD(2) * t551;
t494 = t543 * t624 - t608;
t554 = cos(qJ(5));
t481 = t554 * t494;
t503 = t543 * t555 + t546 * t551;
t497 = t503 * qJD(2);
t550 = sin(qJ(5));
t437 = -t497 * t550 - t481;
t683 = t437 * t553;
t687 = t620 - t683;
t549 = sin(qJ(6));
t619 = -qJD(6) + t437;
t595 = t619 * t549;
t548 = qJ(4) + pkin(8);
t603 = qJD(3) * t548;
t487 = qJD(4) * t555 - t551 * t603;
t488 = -qJD(4) * t551 - t555 * t603;
t545 = sin(pkin(6));
t556 = cos(qJ(2));
t641 = t545 * t556;
t610 = qJD(1) * t641;
t634 = -t487 * t543 + t546 * t488 + t503 * t610;
t647 = t543 * t551;
t502 = -t546 * t555 + t647;
t633 = t546 * t487 + t543 * t488 + t502 * t610;
t580 = t494 * t550 - t554 * t497;
t496 = t503 * qJD(3);
t615 = qJDD(2) * t555;
t517 = t546 * t615;
t616 = qJDD(2) * t551;
t450 = qJD(2) * t496 + t543 * t616 - t517;
t617 = qJD(2) * qJD(3);
t605 = t551 * t617;
t514 = t543 * t605;
t604 = t555 * t617;
t451 = qJDD(2) * t503 + t546 * t604 - t514;
t622 = qJD(5) * t550;
t386 = -qJD(5) * t481 - t550 * t450 + t554 * t451 - t497 * t622;
t538 = qJDD(3) + qJDD(5);
t539 = qJD(3) + qJD(5);
t612 = t553 * t386 + t549 * t538 + t539 * t620;
t621 = qJD(6) * t549;
t374 = t580 * t621 + t612;
t371 = t374 * t553;
t423 = t539 * t549 - t553 * t580;
t597 = t386 * t549 - t553 * t538;
t375 = t423 * qJD(6) + t597;
t648 = t580 * t549;
t421 = -t553 * t539 - t648;
t686 = -t549 * t375 - t687 * t421 + t371;
t370 = t374 * t549;
t563 = qJD(5) * t580 - t554 * t450 - t451 * t550;
t384 = qJDD(6) - t563;
t381 = t549 * t384;
t635 = -t619 * t620 + t381;
t650 = t437 * t539;
t652 = t580 * t539;
t654 = t423 * t580;
t685 = -t437 ^ 2 * MDP(17) + t538 * MDP(20) + (t437 * MDP(16) + MDP(17) * t580 - MDP(27) * t619) * t580 + (t386 - t650) * MDP(18) + (t563 - t652) * MDP(19) + (t687 * t423 + t370) * MDP(23) + (t619 * t683 + t635 + t654) * MDP(25);
t552 = sin(qJ(2));
t625 = qJD(1) * t552;
t611 = t545 * t625;
t594 = qJD(2) * t548 + t611;
t547 = cos(pkin(6));
t626 = qJD(1) * t547;
t462 = t551 * t626 + t555 * t594;
t454 = t543 * t462;
t461 = -t551 * t594 + t555 * t626;
t660 = qJD(3) * pkin(3);
t458 = t461 + t660;
t409 = t546 * t458 - t454;
t666 = pkin(9) * t497;
t393 = qJD(3) * pkin(4) + t409 - t666;
t640 = t546 * t462;
t410 = t543 * t458 + t640;
t667 = pkin(9) * t494;
t396 = t410 - t667;
t372 = t393 * t554 - t396 * t550;
t367 = -pkin(5) * t539 - t372;
t684 = t367 * t437;
t655 = t421 * t580;
t614 = t547 * qJDD(1);
t518 = t555 * t614;
t618 = qJD(1) * qJD(2);
t473 = qJDD(2) * pkin(8) + (qJDD(1) * t552 + t556 * t618) * t545;
t564 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t626 + t473;
t579 = t594 * qJD(3);
t404 = qJDD(3) * pkin(3) - t551 * t564 - t555 * t579 + t518;
t405 = (-t579 + t614) * t551 + t564 * t555;
t379 = t546 * t404 - t543 * t405;
t363 = qJDD(3) * pkin(4) - pkin(9) * t451 + t379;
t380 = t543 * t404 + t546 * t405;
t364 = -pkin(9) * t450 + t380;
t373 = t393 * t550 + t396 * t554;
t671 = -qJD(5) * t373 + t554 * t363 - t550 * t364;
t353 = -pkin(5) * t538 - t671;
t659 = cos(pkin(11));
t600 = t659 * t556;
t544 = sin(pkin(11));
t645 = t544 * t552;
t490 = t547 * t645 - t600;
t601 = t659 * t552;
t644 = t544 * t556;
t492 = t547 * t601 + t644;
t540 = qJ(3) + pkin(12);
t537 = qJ(5) + t540;
t528 = sin(t537);
t529 = cos(t537);
t602 = t545 * t659;
t643 = t545 * t552;
t646 = t544 * t545;
t572 = -g(3) * (-t528 * t643 + t529 * t547) - g(2) * (-t492 * t528 - t529 * t602) - g(1) * (t490 * t528 + t529 * t646);
t568 = -t353 + t572;
t499 = t502 * qJD(3);
t682 = -pkin(9) * t499 - t634;
t681 = -pkin(9) * t496 + t633;
t452 = t554 * t502 + t503 * t550;
t453 = -t502 * t550 + t503 * t554;
t532 = pkin(3) * t555 + pkin(2);
t475 = pkin(4) * t502 - t532;
t389 = pkin(5) * t452 - pkin(10) * t453 + t475;
t491 = -t547 * t600 + t645;
t493 = t547 * t644 + t601;
t589 = g(1) * t493 + g(2) * t491;
t569 = g(3) * t641 - t589;
t566 = t569 * t529;
t680 = t389 * t384 - t566;
t368 = pkin(10) * t539 + t373;
t484 = -qJD(2) * t532 + qJD(4) - t610;
t446 = pkin(4) * t494 + t484;
t385 = -pkin(5) * t437 + pkin(10) * t580 + t446;
t355 = t368 * t553 + t385 * t549;
t679 = -t355 * t580 + t367 * t620 - t549 * t568;
t583 = t368 * t549 - t385 * t553;
t678 = t367 * t621 - t583 * t580;
t443 = t492 * t529 - t528 * t602;
t445 = -t490 * t529 + t528 * t646;
t471 = t528 * t547 + t529 * t643;
t672 = t554 * (qJD(5) * t393 + t364) + t550 * t363 - t396 * t622;
t677 = g(1) * t445 + g(2) * t443 + g(3) * t471 - t446 * t437 - t672;
t676 = t446 * t580 + t572 + t671;
t402 = -pkin(5) * t580 - pkin(10) * t437;
t534 = t551 * t660;
t587 = pkin(4) * t496 + t534 - t611;
t669 = pkin(3) * t546;
t530 = pkin(4) + t669;
t670 = pkin(3) * t543;
t630 = t550 * t530 + t554 * t670;
t382 = t553 * t384;
t574 = -t619 * t621 - t382;
t352 = pkin(10) * t538 + t672;
t510 = t548 * t551;
t511 = t548 * t555;
t466 = -t546 * t510 - t511 * t543;
t427 = -pkin(9) * t503 + t466;
t467 = -t543 * t510 + t546 * t511;
t428 = -pkin(9) * t502 + t467;
t392 = t427 * t550 + t428 * t554;
t406 = -qJD(5) * t452 - t496 * t550 - t499 * t554;
t590 = g(1) * t490 - g(2) * t492;
t570 = -g(3) * t643 + t590;
t582 = t427 * t554 - t428 * t550;
t638 = -qJD(5) * t582 + t550 * t682 - t554 * t681;
t674 = -(qJD(6) * t385 + t352) * t452 + t353 * t453 + t367 * t406 - (-qJD(6) * t389 + t638) * t619 - t392 * t384 + t570;
t557 = qJD(3) ^ 2;
t606 = t552 * t618;
t586 = -qJDD(1) * t641 + t545 * t606;
t673 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t557 + t545 * (-g(3) * t556 + t606) - t586 + t589;
t668 = pkin(3) * t551;
t500 = t547 * t555 - t551 * t643;
t662 = g(3) * t500;
t661 = qJD(2) * pkin(2);
t657 = t367 * t453;
t653 = t423 * t549;
t642 = t545 * t555;
t639 = qJDD(1) - g(3);
t637 = qJD(5) * t392 + t550 * t681 + t554 * t682;
t414 = t546 * t461 - t454;
t412 = -t461 * t543 - t640;
t399 = t412 + t667;
t400 = t414 - t666;
t577 = t530 * t554 - t550 * t670;
t632 = -t577 * qJD(5) + t399 * t550 + t400 * t554;
t631 = t630 * qJD(5) + t399 * t554 - t400 * t550;
t541 = t551 ^ 2;
t629 = -t555 ^ 2 + t541;
t533 = pkin(3) * t624;
t609 = qJD(2) * t641;
t468 = pkin(4) * t497 + t533;
t486 = pkin(10) + t630;
t591 = qJD(6) * t486 + t402 + t468;
t407 = qJD(5) * t453 + t554 * t496 - t499 * t550;
t588 = pkin(5) * t407 - pkin(10) * t406 + t587;
t584 = -t486 * t384 - t684;
t501 = t547 * t551 + t552 * t642;
t432 = t500 * t546 - t501 * t543;
t433 = t500 * t543 + t501 * t546;
t581 = t432 * t554 - t433 * t550;
t395 = t432 * t550 + t433 * t554;
t578 = -t437 * t595 - t574;
t576 = -g(1) * t544 + g(2) * t659;
t573 = t406 * t553 - t453 * t621;
t571 = pkin(3) * t605 + qJDD(4) + t586;
t509 = -t610 - t661;
t567 = -qJD(2) * t509 - t473 - t590;
t562 = -pkin(8) * qJDD(3) + (t509 + t610 - t661) * qJD(3);
t449 = -qJDD(2) * t532 + t571;
t408 = pkin(4) * t450 + t449;
t558 = qJD(2) ^ 2;
t536 = cos(t540);
t535 = sin(t540);
t485 = -pkin(5) - t577;
t460 = -qJD(3) * t501 - t551 * t609;
t459 = qJD(3) * t500 + t555 * t609;
t413 = t459 * t546 + t460 * t543;
t411 = -t459 * t543 + t460 * t546;
t359 = qJD(5) * t395 - t411 * t554 + t413 * t550;
t358 = qJD(5) * t581 + t411 * t550 + t413 * t554;
t357 = -pkin(5) * t563 - pkin(10) * t386 + t408;
t356 = t553 * t357;
t1 = [t639 * MDP(1) + (qJD(3) * t460 + qJDD(3) * t500) * MDP(10) + (-qJD(3) * t459 - qJDD(3) * t501) * MDP(11) + (qJD(3) * t411 + qJDD(3) * t432) * MDP(12) + (-qJD(3) * t413 - qJDD(3) * t433) * MDP(13) + (-t411 * t497 - t413 * t494 - t432 * t451 - t433 * t450) * MDP(14) + (t379 * t432 + t380 * t433 + t409 * t411 + t410 * t413 - g(3)) * MDP(15) + (-t359 * t539 + t538 * t581) * MDP(21) + (-t358 * t539 - t395 * t538) * MDP(22) + (-(-t358 * t549 - t395 * t620) * t619 - t395 * t381 + t359 * t421 - t581 * t375) * MDP(28) + ((t358 * t553 - t395 * t621) * t619 - t395 * t382 + t359 * t423 - t581 * t374) * MDP(29) + ((-qJDD(2) * MDP(4) + (-t555 * MDP(10) + t551 * MDP(11) - MDP(3)) * t558 + (MDP(12) * t494 + MDP(13) * t497 + MDP(15) * t484 - MDP(21) * t437 - MDP(22) * t580 - (MDP(28) * t553 - MDP(29) * t549) * t619) * qJD(2)) * t552 + (qJDD(2) * MDP(3) - t558 * MDP(4) + (-t605 + t615) * MDP(10) + (-t604 - t616) * MDP(11) - t450 * MDP(12) - t451 * MDP(13) - t449 * MDP(15) + t563 * MDP(21) - t386 * MDP(22) + t574 * MDP(28) + t635 * MDP(29)) * t556) * t545; (-t355 * t407 - t582 * t374 + t637 * t423 + (-(-qJD(6) * t368 + t357) * t452 - qJD(6) * t657 - (qJD(6) * t392 - t588) * t619 - t680) * t549 + t674 * t553) * MDP(29) + (-t583 * t407 + t356 * t452 - t582 * t375 + t637 * t421 + (-t588 * t619 + (-t368 * t452 + t392 * t619 + t657) * qJD(6) + t680) * t553 + t674 * t549) * MDP(28) + (-t379 * t503 - t380 * t502 + t409 * t499 - t410 * t496 - t450 * t467 - t451 * t466 - t494 * t633 - t634 * t497 + t570) * MDP(14) + (-t497 * t611 - qJDD(3) * t467 + t449 * t503 - t451 * t532 - t484 * t499 + t569 * t535 + (t497 * t668 - t633) * qJD(3)) * MDP(13) + (t562 * t551 + t673 * t555) * MDP(10) + (-t673 * t551 + t562 * t555) * MDP(11) + (t384 * t452 - t407 * t619) * MDP(27) + (t374 * t452 + t382 * t453 + t407 * t423 - t573 * t619) * MDP(25) + (-t453 * t381 - t375 * t452 - t407 * t421 - (-t406 * t549 - t453 * t620) * t619) * MDP(26) + (t386 * t453 - t406 * t580) * MDP(16) + (t386 * t475 - t392 * t538 + t406 * t446 + t408 * t453 + t528 * t569 + t539 * t638 - t580 * t587) * MDP(22) + qJDD(2) * MDP(2) + (t407 * t446 + t408 * t452 - t437 * t587 - t475 * t563 + t538 * t582 - t539 * t637 - t566) * MDP(21) + (-t386 * t452 + t406 * t437 + t407 * t580 + t453 * t563) * MDP(17) + (qJDD(2) * t541 + 0.2e1 * t551 * t604) * MDP(5) + 0.2e1 * (t551 * t615 - t617 * t629) * MDP(6) + (t380 * t467 + t379 * t466 - t449 * t532 + t484 * t534 - g(1) * (-t490 * t548 - t493 * t532) - g(2) * (-t491 * t532 + t492 * t548) + t633 * t410 + t634 * t409 + (-t484 * t625 - g(3) * (t532 * t556 + t548 * t552)) * t545) * MDP(15) + (t639 * t641 + t589) * MDP(3) + (-t639 * t643 - t590) * MDP(4) + (-t407 * t539 - t452 * t538) * MDP(19) + (t406 * t539 + t453 * t538) * MDP(18) + (t371 * t453 + t423 * t573) * MDP(23) + ((-t421 * t553 - t653) * t406 + (-t370 - t375 * t553 + (t421 * t549 - t423 * t553) * qJD(6)) * t453) * MDP(24) + (-t494 * t611 + qJDD(3) * t466 + t449 * t502 - t450 * t532 + t484 * t496 - t569 * t536 + (t494 * t668 + t634) * qJD(3)) * MDP(12) + (qJDD(3) * t551 + t555 * t557) * MDP(7) + (qJDD(3) * t555 - t551 * t557) * MDP(8); (t437 * t468 + t538 * t577 - t539 * t631 + t676) * MDP(21) + (t468 * t580 - t538 * t630 + t539 * t632 + t677) * MDP(22) + (t485 * t375 + t631 * t421 + (-t619 * t632 + t584) * t549 + (t591 * t619 + t568) * t553 + t678) * MDP(28) + (t485 * t374 + t584 * t553 + t631 * t423 - (t549 * t591 + t553 * t632) * t619 + t679) * MDP(29) + (t414 * qJD(3) + t484 * t494 - g(1) * (t490 * t536 - t535 * t646) - g(2) * (-t492 * t536 + t535 * t602) - g(3) * (-t535 * t547 - t536 * t643) + (-qJDD(3) * t543 - t497 * t624) * pkin(3) - t380) * MDP(13) + (-t412 * qJD(3) - t484 * t497 - g(1) * (t490 * t535 + t536 * t646) - g(2) * (-t492 * t535 - t536 * t602) - g(3) * (-t535 * t643 + t536 * t547) + (qJDD(3) * t546 - t494 * t624) * pkin(3) + t379) * MDP(12) + (t578 - t655) * MDP(26) + (-t551 * t555 * MDP(5) + t629 * MDP(6)) * t558 + (t379 * t669 + t380 * t670 - t409 * t412 - t410 * t414 - t484 * t533 + (-g(1) * (t490 * t551 + t544 * t642) - g(2) * (-t492 * t551 - t555 * t602) - t662) * pkin(3)) * MDP(15) + ((t410 + t412) * t497 + (-t409 + t414) * t494 + (-t450 * t543 - t451 * t546) * pkin(3)) * MDP(14) + t685 + MDP(8) * t615 + MDP(7) * t616 + qJDD(3) * MDP(9) + (t619 * t653 + t686) * MDP(24) + (g(3) * t501 + (-t545 * t576 - t614) * t551 + t567 * t555) * MDP(11) + (t551 * t567 + t576 * t642 + t518 - t662) * MDP(10); -t517 * MDP(12) - t514 * MDP(13) + (-t494 ^ 2 - t497 ^ 2) * MDP(14) + (t409 * t497 + t410 * t494 + t569 + t571) * MDP(15) + (-t563 - t652) * MDP(21) + (t386 + t650) * MDP(22) + (t578 + t655) * MDP(28) + (-t553 * t619 ^ 2 - t381 + t654) * MDP(29) + (MDP(12) * t647 + t503 * MDP(13) - t532 * MDP(15)) * qJDD(2) + ((t543 * t623 + t546 * t624 + t497) * MDP(12) + (-t494 + t608) * MDP(13)) * qJD(3); (t373 * t539 + t676) * MDP(21) + (t372 * t539 + t677) * MDP(22) + (t423 * t595 + t686) * MDP(24) + (-t595 * t619 + t382 - t655) * MDP(26) + (-pkin(5) * t375 - t373 * t421 + (-pkin(10) * t384 - t372 * t619 - t684) * t549 + (-(-pkin(10) * qJD(6) - t402) * t619 + t568) * t553 + t678) * MDP(28) + (-pkin(5) * t374 - (t372 * t553 + t402 * t549) * t619 - t373 * t423 - t367 * t683 + t574 * pkin(10) + t679) * MDP(29) + t685; t423 * t421 * MDP(23) + (-t421 ^ 2 + t423 ^ 2) * MDP(24) + (-t421 * t619 + t612) * MDP(25) + (-t423 * t619 - t597) * MDP(26) + t384 * MDP(27) + (-t549 * t352 + t356 - t355 * t619 - t367 * t423 - g(1) * (-t445 * t549 + t493 * t553) - g(2) * (-t443 * t549 + t491 * t553) - g(3) * (-t471 * t549 - t553 * t641)) * MDP(28) + (-t553 * t352 - t549 * t357 + t583 * t619 + t367 * t421 - g(1) * (-t445 * t553 - t493 * t549) - g(2) * (-t443 * t553 - t491 * t549) - g(3) * (-t471 * t553 + t549 * t641)) * MDP(29) + (MDP(25) * t648 - MDP(26) * t423 - MDP(28) * t355 + MDP(29) * t583) * qJD(6);];
tau = t1;
