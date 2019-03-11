% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:21
% EndTime: 2019-03-09 18:04:33
% DurationCPUTime: 8.08s
% Computational Cost: add. (11267->437), mult. (29690->579), div. (0->0), fcn. (23376->10), ass. (0->224)
t554 = cos(qJ(6));
t619 = qJD(6) * t554;
t555 = cos(qJ(5));
t556 = cos(qJ(3));
t557 = cos(qJ(2));
t623 = qJD(1) * t557;
t610 = t556 * t623;
t552 = sin(qJ(3));
t553 = sin(qJ(2));
t624 = qJD(1) * t553;
t612 = t552 * t624;
t507 = -t610 + t612;
t509 = -t552 * t623 - t556 * t624;
t548 = sin(pkin(11));
t549 = cos(pkin(11));
t579 = t507 * t549 - t509 * t548;
t479 = t555 * t579;
t483 = t507 * t548 + t509 * t549;
t551 = sin(qJ(5));
t437 = -t483 * t551 + t479;
t683 = t437 * t554;
t686 = t619 + t683;
t550 = sin(qJ(6));
t620 = qJD(6) * t550;
t545 = qJD(2) + qJD(3);
t615 = qJD(1) * qJD(2);
t609 = t557 * t615;
t487 = qJD(3) * t610 - t545 * t612 + t556 * t609;
t521 = t552 * t557 + t553 * t556;
t492 = t545 * t521;
t488 = t492 * qJD(1);
t444 = -t487 * t548 - t488 * t549;
t445 = t487 * t549 - t488 * t548;
t621 = qJD(5) * t551;
t384 = -qJD(5) * t479 + t551 * t444 + t555 * t445 + t483 * t621;
t544 = qJD(5) + t545;
t631 = t554 * t384 + t544 * t619;
t672 = -t555 * t483 - t551 * t579;
t359 = -t620 * t672 + t631;
t358 = t359 * t554;
t429 = t544 * t550 + t554 * t672;
t654 = t384 * t550;
t360 = t429 * qJD(6) + t654;
t645 = t672 * t550;
t427 = -t554 * t544 + t645;
t685 = -t550 * t360 - t686 * t427 + t358;
t357 = t359 * t550;
t385 = qJD(5) * t672 - t555 * t444 + t551 * t445;
t381 = t550 * t385;
t682 = -qJD(6) - t437;
t632 = -t619 * t682 + t381;
t647 = t437 * t544;
t649 = t672 * t544;
t651 = t429 * t672;
t684 = (-t385 + t649) * MDP(23) - t437 ^ 2 * MDP(21) + (MDP(20) * t437 + MDP(21) * t672 + MDP(31) * t682) * t672 + (t384 + t647) * MDP(22) + (t686 * t429 + t357) * MDP(27) + (-t682 * t683 + t632 - t651) * MDP(29);
t505 = t509 * qJ(4);
t663 = pkin(7) + pkin(8);
t530 = t663 * t557;
t524 = qJD(1) * t530;
t510 = t552 * t524;
t529 = t663 * t553;
t522 = qJD(1) * t529;
t658 = qJD(2) * pkin(2);
t516 = -t522 + t658;
t600 = t556 * t516 - t510;
t467 = t505 + t600;
t458 = pkin(3) * t545 + t467;
t514 = t556 * t524;
t578 = -t516 * t552 - t514;
t657 = qJ(4) * t507;
t468 = -t578 - t657;
t459 = t548 * t468;
t413 = t549 * t458 - t459;
t480 = t483 * pkin(9);
t398 = pkin(4) * t545 + t413 + t480;
t641 = t549 * t468;
t414 = t548 * t458 + t641;
t660 = pkin(9) * t579;
t400 = t414 - t660;
t371 = t398 * t555 - t400 * t551;
t369 = -pkin(5) * t544 - t371;
t656 = t369 * t437;
t397 = pkin(5) * t672 + pkin(10) * t437;
t613 = qJD(2) * t663;
t589 = qJD(1) * t613;
t518 = t557 * t589;
t622 = qJD(3) * t552;
t598 = -t552 * t518 - t524 * t622;
t517 = t553 * t589;
t668 = (qJD(3) * t516 - t517) * t556;
t408 = -qJ(4) * t488 - qJD(4) * t507 + t598 + t668;
t599 = t552 * t517 - t556 * t518;
t565 = qJD(3) * t578 + t599;
t409 = -qJ(4) * t487 + qJD(4) * t509 + t565;
t379 = -t408 * t548 + t549 * t409;
t364 = -pkin(9) * t445 + t379;
t380 = t549 * t408 + t548 * t409;
t365 = pkin(9) * t444 + t380;
t344 = t555 * (qJD(5) * t398 + t365) + t364 * t551 - t400 * t621;
t541 = -pkin(2) * t557 - pkin(1);
t528 = t541 * qJD(1);
t493 = t507 * pkin(3) + qJD(4) + t528;
t451 = pkin(4) * t579 + t493;
t563 = t437 * t451 - t344;
t652 = t427 * t672;
t597 = t522 * t552 - t514;
t470 = t597 + t657;
t627 = -t556 * t522 - t510;
t471 = t505 + t627;
t640 = t549 * t552;
t659 = pkin(2) * qJD(3);
t628 = -t549 * t470 + t471 * t548 + (-t548 * t556 - t640) * t659;
t642 = t548 * t552;
t667 = -t548 * t470 - t549 * t471 + (t549 * t556 - t642) * t659;
t372 = t398 * t551 + t400 * t555;
t345 = qJD(5) * t372 - t555 * t364 + t365 * t551;
t562 = -t451 * t672 - t345;
t370 = pkin(10) * t544 + t372;
t388 = t437 * pkin(5) - pkin(10) * t672 + t451;
t354 = t370 * t554 + t388 * t550;
t588 = t345 * t550 + t354 * t672 + t369 * t619;
t585 = t370 * t550 - t388 * t554;
t572 = -t345 * t554 + t369 * t620 + t585 * t672;
t596 = t682 * t550;
t674 = t660 - t628;
t673 = t480 - t667;
t671 = -0.2e1 * t615;
t670 = MDP(4) * t553;
t669 = MDP(5) * (t553 ^ 2 - t557 ^ 2);
t383 = t554 * t385;
t666 = -t620 * t682 - t383;
t665 = qJD(1) * t521;
t520 = t552 * t553 - t556 * t557;
t523 = t553 * t613;
t525 = t557 * t613;
t643 = t529 * t556;
t569 = -qJD(3) * t643 - t556 * t523 - t552 * t525 - t530 * t622;
t421 = -qJ(4) * t492 - qJD(4) * t520 + t569;
t491 = t545 * t520;
t577 = t529 * t552 - t530 * t556;
t564 = qJD(3) * t577 + t523 * t552 - t556 * t525;
t422 = qJ(4) * t491 - qJD(4) * t521 + t564;
t393 = -t421 * t548 + t549 * t422;
t450 = -t491 * t549 - t492 * t548;
t377 = -pkin(9) * t450 + t393;
t394 = t549 * t421 + t548 * t422;
t449 = t491 * t548 - t492 * t549;
t378 = pkin(9) * t449 + t394;
t481 = -qJ(4) * t521 - t530 * t552 - t643;
t482 = -qJ(4) * t520 - t577;
t431 = t549 * t481 - t482 * t548;
t490 = -t520 * t548 + t521 * t549;
t411 = -pkin(9) * t490 + t431;
t432 = t548 * t481 + t549 * t482;
t489 = -t520 * t549 - t521 * t548;
t412 = pkin(9) * t489 + t432;
t584 = t411 * t555 - t412 * t551;
t346 = qJD(5) * t584 + t377 * t551 + t378 * t555;
t387 = t411 * t551 + t412 * t555;
t582 = t555 * t489 - t490 * t551;
t391 = qJD(5) * t582 + t449 * t551 + t450 * t555;
t447 = t489 * t551 + t490 * t555;
t576 = pkin(3) * t520 + t541;
t463 = -pkin(4) * t489 + t576;
t395 = -pkin(5) * t582 - pkin(10) * t447 + t463;
t664 = t345 * t447 + t369 * t391 - t387 * t385 + (qJD(6) * t395 + t346) * t682 + (qJD(6) * t388 + t344) * t582;
t662 = pkin(3) * t509;
t661 = pkin(3) * t548;
t655 = t369 * t447;
t653 = t395 * t385;
t650 = t429 * t550;
t644 = t528 * t509;
t558 = qJD(2) ^ 2;
t639 = t553 * t558;
t638 = t557 * t558;
t559 = qJD(1) ^ 2;
t637 = t557 * t559;
t415 = -t467 * t548 - t641;
t401 = t415 + t660;
t416 = t549 * t467 - t459;
t402 = t480 + t416;
t538 = pkin(3) * t549 + pkin(4);
t573 = t538 * t555 - t551 * t661;
t634 = -t573 * qJD(5) + t401 * t551 + t402 * t555;
t540 = pkin(2) * t556 + pkin(3);
t503 = -pkin(2) * t642 + t549 * t540;
t496 = pkin(4) + t503;
t504 = pkin(2) * t640 + t540 * t548;
t581 = t496 * t555 - t504 * t551;
t633 = -qJD(5) * t581 + t551 * t674 + t555 * t673;
t580 = t496 * t551 + t504 * t555;
t630 = qJD(5) * t580 - t551 * t673 + t555 * t674;
t574 = t538 * t551 + t555 * t661;
t629 = t574 * qJD(5) + t401 * t555 - t402 * t551;
t543 = t553 * t658;
t542 = pkin(2) * t624;
t608 = -pkin(2) * t545 - t516;
t607 = pkin(3) * t488 + qJD(2) * t542;
t606 = pkin(3) * t492 + t543;
t605 = pkin(1) * t671;
t455 = -pkin(4) * t483 - t662;
t390 = t397 + t455;
t473 = pkin(10) + t580;
t591 = qJD(6) * t473 + t390 + t542;
t498 = pkin(10) + t574;
t590 = qJD(6) * t498 + t390;
t587 = -t385 * t473 + t656;
t586 = -t385 * t498 + t656;
t583 = -t413 * t579 - t414 * t483;
t410 = -pkin(4) * t444 + t607;
t420 = -pkin(4) * t449 + t606;
t575 = t437 * t596 - t666;
t571 = t528 * t507 - t598;
t570 = t391 * t554 - t447 * t620;
t560 = -t509 * t507 * MDP(11) + (t650 * t682 + t685) * MDP(28) + (t575 + t652) * MDP(30) + t487 * MDP(13) + (-t507 ^ 2 + t509 ^ 2) * MDP(12) + (t507 * MDP(13) + (-t509 - t665) * MDP(14)) * t545 + t684;
t497 = -pkin(5) - t573;
t472 = -pkin(5) - t581;
t452 = t455 + t542;
t392 = qJD(5) * t447 - t555 * t449 + t450 * t551;
t355 = pkin(5) * t392 - pkin(10) * t391 + t420;
t352 = pkin(5) * t385 - pkin(10) * t384 + t410;
t351 = t554 * t352;
t347 = qJD(5) * t387 - t377 * t555 + t378 * t551;
t1 = [0.2e1 * t609 * t670 + (t541 * t488 + t528 * t492 + (qJD(1) * t520 + t507) * t543) * MDP(16) + (MDP(22) * t391 - MDP(23) * t392 - MDP(25) * t347 - MDP(26) * t346) * t544 + (t385 * t463 + t392 * t451 - t410 * t582 + t420 * t437) * MDP(25) + (-t379 * t490 + t380 * t489 + t393 * t483 - t394 * t579 - t413 * t450 + t414 * t449 - t431 * t445 + t432 * t444) * MDP(18) + (-t491 * MDP(13) - t492 * MDP(14) + t564 * MDP(16) - t569 * MDP(17)) * t545 + (t541 * t487 - t528 * t491 + (-t509 + t665) * t543) * MDP(17) + (-t487 * t520 - t488 * t521 + t491 * t507 + t492 * t509) * MDP(12) + (t487 * t521 + t491 * t509) * MDP(11) + MDP(6) * t638 + (t384 * t582 - t385 * t447 - t391 * t437 - t392 * t672) * MDP(21) + (t384 * t447 + t391 * t672) * MDP(20) + (t384 * t463 + t391 * t451 + t410 * t447 + t420 * t672) * MDP(26) + (t347 * t427 - t351 * t582 - t585 * t392 - t584 * t360 + (-t355 * t682 + t653 + (t370 * t582 + t387 * t682 + t655) * qJD(6)) * t554 + t664 * t550) * MDP(32) + (t347 * t429 - t354 * t392 - t584 * t359 + ((-qJD(6) * t387 + t355) * t682 - t653 + (-qJD(6) * t370 + t352) * t582 - qJD(6) * t655) * t550 + t664 * t554) * MDP(33) + (-t447 * t381 + t360 * t582 - t392 * t427 - (-t391 * t550 - t447 * t619) * t682) * MDP(30) + (-t385 * t582 - t392 * t682) * MDP(31) + (-t359 * t582 + t383 * t447 + t392 * t429 - t570 * t682) * MDP(29) + (t447 * t358 + t570 * t429) * MDP(27) + ((-t427 * t554 - t650) * t391 + (-t357 - t360 * t554 + (t427 * t550 - t429 * t554) * qJD(6)) * t447) * MDP(28) - MDP(7) * t639 + (pkin(7) * t639 + t557 * t605) * MDP(10) + (-pkin(7) * t638 + t553 * t605) * MDP(9) + (t379 * t431 + t380 * t432 + t413 * t393 + t414 * t394 + t493 * t606 + t576 * t607) * MDP(19) + t669 * t671; t560 + (t472 * t359 + t587 * t554 + t630 * t429 - (t550 * t591 + t554 * t633) * t682 + t588) * MDP(33) + (-t452 * t672 + t544 * t633 + t563) * MDP(26) + (t509 * t542 + t627 * t545 + (qJD(3) * t608 + t517) * t556 + t571) * MDP(17) + (-t507 * t542 + t644 - t597 * t545 + (t552 * t608 - t514) * qJD(3) + t599) * MDP(16) + (-t437 * t452 - t544 * t630 + t562) * MDP(25) + (t504 * t444 - t503 * t445 + t483 * t628 - t579 * t667 + t583) * MDP(18) - t637 * t670 + (t380 * t504 + t379 * t503 - t493 * (t542 - t662) + t667 * t414 + t628 * t413) * MDP(19) + t559 * t669 + (t472 * t360 + t587 * t550 + t630 * t427 - (t550 * t633 - t554 * t591) * t682 + t572) * MDP(32) + (MDP(9) * t553 * t559 + MDP(10) * t637) * pkin(1); t560 + (-t437 * t455 - t544 * t629 + t562) * MDP(25) + (-t413 * t415 - t414 * t416 + (t379 * t549 + t380 * t548 + t493 * t509) * pkin(3)) * MDP(19) + (t416 * t579 - t415 * t483 + (t444 * t548 - t445 * t549) * pkin(3) + t583) * MDP(18) + (-t545 * t578 + t565 + t644) * MDP(16) + (t545 * t600 + t571 - t668) * MDP(17) + (t497 * t359 + t586 * t554 + t629 * t429 - (t550 * t590 + t554 * t634) * t682 + t588) * MDP(33) + (-t455 * t672 + t544 * t634 + t563) * MDP(26) + (t497 * t360 + t586 * t550 + t629 * t427 - (t550 * t634 - t554 * t590) * t682 + t572) * MDP(32); (-t483 ^ 2 - t579 ^ 2) * MDP(18) + (-t413 * t483 + t414 * t579 + t607) * MDP(19) + (t385 + t649) * MDP(25) + (t384 - t647) * MDP(26) + (t575 - t652) * MDP(32) + (-t554 * t682 ^ 2 - t381 - t651) * MDP(33); (t372 * t544 + t562) * MDP(25) + (t371 * t544 + t563) * MDP(26) + (t429 * t596 + t685) * MDP(28) + (-t596 * t682 + t383 + t652) * MDP(30) + (-pkin(5) * t360 + (-t371 * t550 + t397 * t554) * t682 - t372 * t427 + t550 * t656 - t632 * pkin(10) + t572) * MDP(32) + (-pkin(5) * t359 - (t371 * t554 + t397 * t550) * t682 - t372 * t429 + t369 * t683 + t666 * pkin(10) + t588) * MDP(33) + t684; t429 * t427 * MDP(27) + (-t427 ^ 2 + t429 ^ 2) * MDP(28) + (-t427 * t682 + t631) * MDP(29) + (-t429 * t682 - t654) * MDP(30) + t385 * MDP(31) + (-t344 * t550 - t354 * t682 - t369 * t429 + t351) * MDP(32) + (-t344 * t554 - t352 * t550 + t369 * t427 + t585 * t682) * MDP(33) + (-MDP(29) * t645 - MDP(30) * t429 - MDP(32) * t354 + MDP(33) * t585) * qJD(6);];
tauc  = t1;
