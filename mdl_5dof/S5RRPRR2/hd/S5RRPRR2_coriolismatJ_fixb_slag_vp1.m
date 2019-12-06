% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:05
% EndTime: 2019-12-05 18:27:42
% DurationCPUTime: 14.89s
% Computational Cost: add. (56152->637), mult. (40761->815), div. (0->0), fcn. (37944->10), ass. (0->398)
t493 = sin(qJ(1));
t489 = t493 ^ 2;
t495 = cos(qJ(1));
t490 = t495 ^ 2;
t580 = t489 + t490;
t488 = qJ(2) + pkin(9);
t467 = qJ(4) + t488;
t460 = sin(t467);
t697 = pkin(4) * t460;
t463 = qJ(5) + t467;
t451 = cos(t463);
t647 = t451 * t493;
t450 = sin(t463);
t651 = t450 * t493;
t358 = rSges(6,1) * t651 + rSges(6,2) * t647;
t393 = rSges(6,1) * t450 + rSges(6,2) * t451;
t359 = t393 * t495;
t772 = t493 * t358 + t495 * t359;
t188 = -t580 * t697 - t772;
t465 = sin(t488);
t492 = sin(qJ(2));
t698 = pkin(2) * t492;
t536 = pkin(3) * t465 + t698;
t392 = -t536 - t697;
t505 = -t392 + t393;
t229 = t505 * t493;
t231 = t505 * t495;
t688 = rSges(6,1) * t451;
t394 = -rSges(6,2) * t450 + t688;
t461 = cos(t467);
t696 = pkin(4) * t461;
t568 = -t394 - t696;
t273 = t568 * t493;
t275 = t568 * t495;
t615 = -t229 * t273 - t231 * t275;
t305 = rSges(6,1) * t647 - rSges(6,2) * t651 - t495 * rSges(6,3);
t650 = t450 * t495;
t561 = -rSges(6,2) * t650 + t493 * rSges(6,3);
t646 = t451 * t495;
t202 = t493 * t305 + t495 * (rSges(6,1) * t646 + t561);
t466 = cos(t488);
t459 = pkin(3) * t466;
t494 = cos(qJ(2));
t485 = t494 * pkin(2);
t464 = t485 + pkin(1);
t418 = t459 + t464;
t373 = t418 + t696;
t491 = -qJ(3) - pkin(6);
t487 = -pkin(7) + t491;
t572 = -pkin(8) + t487;
t448 = t493 * t572;
t457 = t493 * t487;
t589 = -t493 * t373 - t495 * t572;
t633 = t487 * t495;
t121 = -t493 * (t493 * t418 + t589 + t633) + t202 + t495 * (-t448 + t457 + (t373 - t418) * t495);
t462 = t493 * t491;
t584 = t418 - t464;
t486 = t495 * pkin(6);
t632 = t491 * t495;
t695 = pkin(1) - t464;
t600 = -t493 * (t695 * t493 - t486 - t632) + t495 * (-t493 * pkin(6) - t695 * t495 - t462);
t543 = -t493 * ((-t487 + t491) * t495 - t584 * t493) + t495 * (t495 * t584 - t457 + t462) + t600;
t80 = t543 + t121;
t55 = t80 * t188 + t615;
t642 = t460 * t495;
t562 = -rSges(5,2) * t642 + t493 * rSges(5,3);
t643 = t460 * t493;
t583 = rSges(5,2) * t643 + t495 * rSges(5,3);
t638 = t461 * t495;
t639 = t461 * t493;
t213 = t493 * (rSges(5,1) * t639 - t583) + t495 * (rSges(5,1) * t638 + t562);
t114 = t543 + t213;
t689 = rSges(5,1) * t461;
t401 = -rSges(5,2) * t460 + t689;
t400 = rSges(5,1) * t460 + rSges(5,2) * t461;
t371 = t400 * t493;
t372 = t400 * t495;
t771 = t493 * t371 + t495 * t372;
t516 = t400 + t536;
t261 = t516 * t495;
t785 = t516 * t493;
t794 = t493 * t785;
t773 = t261 * t495 + t794;
t67 = -t114 * t771 + t773 * t401;
t799 = -m(5) * t67 - m(6) * t55;
t235 = -t392 * t493 + t358;
t364 = t495 * t392;
t236 = t364 - t359;
t410 = t495 * t536;
t264 = -t410 - t372;
t752 = m(6) / 0.2e1;
t753 = m(5) / 0.2e1;
t754 = m(4) / 0.2e1;
t416 = rSges(4,1) * t465 + rSges(4,2) * t466;
t507 = t416 + t698;
t777 = t507 * t495;
t778 = t507 * t493;
t781 = t493 * t778 + t495 * t777;
t573 = (t493 * t235 - t236 * t495) * t752 + (-t264 * t495 + t794) * t753 + t781 * t754;
t775 = t493 * t229 + t231 * t495;
t789 = -m(6) / 0.2e1;
t790 = -m(5) / 0.2e1;
t574 = t775 * t789 + t773 * t790 - m(4) * t781 / 0.2e1;
t33 = t574 - t573;
t798 = t33 * qJD(1);
t569 = t393 + t697;
t272 = t569 * t493;
t534 = t569 * t495;
t666 = t534 * t495;
t768 = t272 * t493 + t666;
t795 = t400 * t580;
t613 = t768 * t789 + t790 * t795;
t270 = pkin(4) * t643 + t358;
t614 = (t493 * t270 + t666) * t752 + t771 * t753;
t73 = t614 - t613;
t797 = t73 * qJD(1);
t796 = t393 * t580;
t446 = Icges(6,4) * t451;
t388 = -Icges(6,2) * t450 + t446;
t389 = Icges(6,1) * t450 + t446;
t793 = t388 + t389;
t447 = Icges(5,4) * t461;
t397 = -Icges(5,2) * t460 + t447;
t398 = Icges(5,1) * t460 + t447;
t792 = t398 + t397;
t791 = -Icges(3,5) * t492 - Icges(4,5) * t465 - Icges(3,6) * t494 - Icges(4,6) * t466;
t728 = t493 / 0.2e1;
t726 = -t495 / 0.2e1;
t786 = t495 / 0.2e1;
t637 = t465 * t493;
t635 = t466 * t493;
t631 = t492 * t493;
t625 = t493 * t494;
t564 = t418 + t689;
t233 = -t493 * t564 + t583 - t633;
t234 = t495 * t564 - t457 + t562;
t774 = t233 * t495 + t234 * t493;
t511 = t774 * t401;
t210 = -t305 + t589;
t211 = -t448 + (t373 + t688) * t495 + t561;
t616 = t275 * t210 + t273 * t211;
t692 = (-t235 * t534 - t236 * t272 + t616) * t752 + (-t511 + (-t264 * t493 - t495 * t785) * t400) * t753;
t693 = (t229 * t534 - t231 * t270 + t616) * t752 + (-t261 * t371 + t372 * t785 - t511) * t753;
t784 = t692 - t693;
t729 = -t493 / 0.2e1;
t783 = t728 + t729;
t678 = Icges(4,4) * t465;
t412 = Icges(4,2) * t466 + t678;
t415 = Icges(4,1) * t466 - t678;
t679 = Icges(3,4) * t492;
t433 = Icges(3,2) * t494 + t679;
t436 = Icges(3,1) * t494 - t679;
t782 = -(t415 / 0.2e1 - t412 / 0.2e1) * t465 - (t436 / 0.2e1 - t433 / 0.2e1) * t492;
t776 = t210 * t495 + t211 * t493;
t770 = t791 * t493;
t769 = t791 * t495;
t458 = Icges(4,4) * t466;
t413 = -Icges(4,2) * t465 + t458;
t414 = Icges(4,1) * t465 + t458;
t480 = Icges(3,4) * t494;
t434 = -Icges(3,2) * t492 + t480;
t435 = Icges(3,1) * t492 + t480;
t320 = Icges(5,6) * t493 + t397 * t495;
t677 = Icges(5,4) * t460;
t399 = Icges(5,1) * t461 - t677;
t322 = Icges(5,5) * t493 + t399 * t495;
t267 = t322 * t639;
t395 = Icges(5,5) * t461 - Icges(5,6) * t460;
t655 = t395 * t495;
t318 = Icges(5,3) * t493 + t655;
t551 = t318 * t495 - t267;
t157 = -t320 * t643 - t551;
t319 = Icges(5,4) * t639 - Icges(5,2) * t643 - Icges(5,6) * t495;
t317 = Icges(5,5) * t639 - Icges(5,6) * t643 - Icges(5,3) * t495;
t429 = Icges(5,4) * t643;
t321 = Icges(5,1) * t639 - Icges(5,5) * t495 - t429;
t609 = -t493 * t317 - t321 * t638;
t158 = -t319 * t642 - t609;
t608 = t493 * t318 + t322 * t638;
t159 = -t320 * t642 + t608;
t548 = t320 * t460 - t317;
t662 = t319 * t460;
t764 = (-t158 * t495 + t159 * t493) * t786 + ((t157 - t267 + (t318 + t662) * t495 + t609) * t495 + t608 * t493) * t726 + ((t493 * t548 + t157 + t158 + t551) * t493 + (t493 * (-t321 * t461 + t662) + t159 - t608 + (t317 + t548) * t495) * t495) * t728;
t396 = Icges(5,2) * t461 + t677;
t763 = t792 * t461 / 0.2e1 + (t399 / 0.2e1 - t396 / 0.2e1) * t460;
t676 = Icges(6,4) * t450;
t387 = Icges(6,2) * t451 + t676;
t390 = Icges(6,1) * t451 - t676;
t542 = t793 * t451 / 0.2e1 + (-t387 / 0.2e1 + t390 / 0.2e1) * t450;
t294 = Icges(6,6) * t493 + t388 * t495;
t296 = Icges(6,5) * t493 + t390 * t495;
t251 = t296 * t647;
t386 = Icges(6,5) * t451 - Icges(6,6) * t450;
t657 = t386 * t495;
t292 = Icges(6,3) * t493 + t657;
t552 = t292 * t495 - t251;
t143 = -t294 * t651 - t552;
t293 = Icges(6,4) * t647 - Icges(6,2) * t651 - Icges(6,6) * t495;
t291 = Icges(6,5) * t647 - Icges(6,6) * t651 - Icges(6,3) * t495;
t421 = Icges(6,4) * t651;
t295 = Icges(6,1) * t647 - Icges(6,5) * t495 - t421;
t611 = -t493 * t291 - t295 * t646;
t144 = -t293 * t650 - t611;
t610 = t493 * t292 + t296 * t646;
t145 = -t294 * t650 + t610;
t549 = t294 * t450 - t291;
t664 = t293 * t450;
t567 = ((t143 - t251 + (t292 + t664) * t495 + t611) * t495 + t610 * t493) * t726 + (-t144 * t495 + t145 * t493) * t786 + ((t493 * t549 + t143 + t144 + t552) * t493 + (t145 - t610 + t493 * (-t295 * t451 + t664) + (t549 + t291) * t495) * t495) * t728;
t379 = Icges(3,5) * t493 + t436 * t495;
t585 = -t433 * t495 + t379;
t377 = Icges(3,6) * t493 + t434 * t495;
t587 = -t435 * t495 - t377;
t357 = Icges(4,5) * t493 + t415 * t495;
t590 = -t412 * t495 + t357;
t355 = Icges(4,6) * t493 + t413 * t495;
t592 = -t414 * t495 - t355;
t762 = -t585 * t631 + t587 * t625 - t590 * t637 + t592 * t635;
t454 = Icges(3,4) * t631;
t378 = Icges(3,1) * t625 - Icges(3,5) * t495 - t454;
t586 = -Icges(3,2) * t625 + t378 - t454;
t376 = Icges(3,4) * t625 - Icges(3,2) * t631 - Icges(3,6) * t495;
t588 = t435 * t493 + t376;
t443 = Icges(4,4) * t637;
t356 = Icges(4,1) * t635 - Icges(4,5) * t495 - t443;
t591 = -Icges(4,2) * t635 + t356 - t443;
t354 = Icges(4,4) * t635 - Icges(4,2) * t637 - Icges(4,6) * t495;
t593 = t414 * t493 + t354;
t761 = t591 * t465 + t593 * t466 + t586 * t492 + t588 * t494;
t596 = -t396 * t495 + t322;
t597 = -Icges(5,2) * t639 + t321 - t429;
t598 = -t398 * t495 - t320;
t599 = t398 * t493 + t319;
t760 = (-t596 * t493 + t495 * t597) * t460 + (t598 * t493 + t495 * t599) * t461;
t602 = -t387 * t495 + t296;
t603 = -Icges(6,2) * t647 + t295 - t421;
t604 = -t389 * t495 - t294;
t605 = t389 * t493 + t293;
t759 = (-t602 * t493 + t495 * t603) * t450 + (t604 * t493 + t495 * t605) * t451;
t758 = 0.4e1 * qJD(1);
t757 = 2 * qJD(2);
t755 = 2 * qJD(4);
t101 = t121 * t772;
t64 = t80 * t772;
t743 = m(6) * (-t101 - t64 + ((t231 + t534) * t495 + (t229 + t272) * t493) * t394);
t132 = t202 * t188;
t510 = (-t273 * t493 - t275 * t495) * t393;
t58 = t775 * t394 - t64;
t741 = m(6) * (t132 + t510 + t58);
t133 = t489 * (t392 + t536) + t495 * (t410 + t364) - t772;
t72 = t768 * t394 - t101;
t497 = t510 + t72;
t740 = m(6) * (t202 * t133 + t497);
t612 = -t272 * t273 - t275 * t534;
t737 = m(6) * (t121 * t133 + t612);
t512 = t776 * t394;
t734 = m(6) * (t229 * t359 - t231 * t358 - t512);
t733 = m(6) * (-t512 + (-t235 * t495 - t236 * t493) * t393);
t732 = m(6) * (t272 * t359 - t358 * t534 - t512);
t731 = m(6) * (-t512 + (-t270 * t495 + t493 * t534) * t393);
t691 = rSges(3,1) * t494;
t571 = pkin(1) + t691;
t581 = rSges(3,2) * t631 + t495 * rSges(3,3);
t311 = -t493 * t571 + t486 + t581;
t630 = t492 * t495;
t456 = rSges(3,2) * t630;
t312 = -t456 + t571 * t495 + (rSges(3,3) + pkin(6)) * t493;
t437 = rSges(3,1) * t492 + rSges(3,2) * t494;
t408 = t437 * t493;
t409 = t437 * t495;
t724 = m(3) * (t311 * t408 - t312 * t409);
t690 = rSges(4,1) * t466;
t565 = t464 + t690;
t582 = rSges(4,2) * t637 + t495 * rSges(4,3);
t255 = -t493 * t565 + t582 - t632;
t636 = t465 * t495;
t563 = -rSges(4,2) * t636 + t493 * rSges(4,3);
t256 = t495 * t565 - t462 + t563;
t723 = m(4) * (t255 * t778 - t256 * t777);
t722 = m(4) * (t255 * t495 + t256 * t493);
t120 = -t213 * t771 + t401 * t795;
t119 = m(5) * t120;
t718 = m(5) * (t233 * t785 + t234 * t264);
t717 = m(5) * (t233 * t371 - t234 * t372);
t716 = m(5) * t774;
t710 = m(6) * (t210 * t235 + t211 * t236);
t112 = -t202 * t772 + t394 * t796;
t709 = m(6) * t112;
t708 = m(6) * (t210 * t270 - t211 * t534);
t707 = m(6) * (t210 * t358 - t211 * t359);
t706 = m(6) * t776;
t701 = m(6) * (-t273 * t495 + t493 * t275);
t524 = Icges(6,5) * t450 + Icges(6,6) * t451;
t345 = t524 * t493;
t346 = t495 * t524;
t694 = (-t489 * t346 + (t493 * t345 + t759) * t495) * t728 + (-t490 * t345 + (t495 * t346 + t759) * t493) * t726;
t660 = t354 * t465;
t658 = t376 * t492;
t634 = t466 * t495;
t624 = t494 * t495;
t352 = Icges(4,5) * t635 - Icges(4,6) * t637 - Icges(4,3) * t495;
t607 = -t493 * t352 - t356 * t634;
t527 = Icges(4,5) * t466 - Icges(4,6) * t465;
t353 = Icges(4,3) * t493 + t495 * t527;
t606 = t493 * t353 + t357 * t634;
t374 = Icges(3,5) * t625 - Icges(3,6) * t631 - Icges(3,3) * t495;
t595 = -t493 * t374 - t378 * t624;
t529 = Icges(3,5) * t494 - Icges(3,6) * t492;
t375 = Icges(3,3) * t493 + t495 * t529;
t594 = t493 * t375 + t379 * t624;
t513 = t772 * t752;
t541 = m(6) * t796;
t151 = t513 + t541 / 0.2e1;
t579 = t151 * qJD(1);
t578 = t740 / 0.2e1 + t694;
t577 = -t709 + t694;
t576 = t121 * t188 + t612;
t570 = rSges(4,2) * t465 - t485 - t690;
t525 = Icges(5,5) * t460 + Icges(5,6) * t461;
t365 = t525 * t493;
t366 = t495 * t525;
t566 = (-t489 * t366 + (t493 * t365 + t760) * t495) * t728 + (-t490 * t365 + (t495 * t366 + t760) * t493) * t726 + t694;
t276 = t357 * t635;
t550 = t353 * t495 - t276;
t326 = t379 * t625;
t547 = t375 * t495 - t326;
t546 = t355 * t465 - t352;
t545 = t377 * t492 - t374;
t540 = t580 * t698;
t537 = t119 + t566;
t535 = -t485 - t459;
t515 = -t401 + t535;
t509 = t489 * (-t536 + t698) + t495 * (pkin(2) * t630 - t410) - t540;
t508 = t567 + t764;
t501 = (-t387 + t390) * t451 - t793 * t450;
t506 = -t567 + (t493 * t386 + t450 * t604 + t451 * t602 + t495 * t501) * t728 + (-t450 * t605 + t451 * t603 + t493 * t501 - t657) * t726;
t504 = t535 + t568;
t503 = -t542 + t783 * (t293 * t451 + t295 * t450);
t502 = t542 + t763;
t500 = (-t396 + t399) * t461 - t792 * t460;
t498 = t506 - t764 + (t493 * t395 + t460 * t598 + t461 * t596 + t495 * t500) * t728 + (-t460 * t599 + t461 * t597 + t493 * t500 - t655) * t726;
t496 = t503 - t763 + t783 * (t319 * t461 + t321 * t460);
t439 = -rSges(3,2) * t492 + t691;
t344 = t570 * t495;
t342 = t570 * t493;
t262 = t515 * t495;
t260 = t515 * t493;
t232 = t504 * t495;
t230 = t504 * t493;
t193 = -t377 * t630 + t594;
t192 = -t376 * t630 - t595;
t191 = -t377 * t631 - t547;
t181 = qJD(4) * t701;
t166 = -t355 * t636 + t606;
t165 = -t354 * t636 - t607;
t164 = -t355 * t637 - t550;
t150 = t513 - t541 / 0.2e1;
t141 = t509 - t771;
t126 = -t192 * t495 + t193 * t493;
t125 = -(-t493 * (-t378 * t494 + t658) - t374 * t495) * t495 + t191 * t493;
t117 = -t165 * t495 + t166 * t493;
t116 = -(-t493 * (-t356 * t466 + t660) - t352 * t495) * t495 + t164 * t493;
t113 = t509 + t133;
t88 = t731 / 0.2e1;
t87 = t732 / 0.2e1;
t86 = t542 + t707;
t84 = t733 / 0.2e1;
t82 = t734 / 0.2e1;
t78 = t706 + t716 + t722;
t74 = t613 + t614;
t61 = (t191 - t326 + (t375 + t658) * t495 + t595) * t495 + t594 * t493;
t60 = (t495 * t545 + t193 - t594) * t495 + (t493 * t545 + t192 + t547) * t493;
t56 = t502 + t708 + t717;
t54 = (t164 - t276 + (t353 + t660) * t495 + t607) * t495 + t606 * t493;
t53 = (t495 * t546 + t166 - t606) * t495 + (t493 * t546 + t165 + t550) * t493;
t34 = t573 + t574;
t30 = t741 / 0.2e1;
t27 = t743 / 0.2e1;
t24 = t694 + t709;
t23 = t24 * qJD(5);
t22 = t502 + (t435 / 0.2e1 + t434 / 0.2e1) * t494 + (t414 / 0.2e1 + t413 / 0.2e1) * t466 + t724 + t718 + t723 + t710 - t782;
t21 = m(6) * t72 + t694;
t20 = m(6) * t58 + t694;
t16 = t537 + t737;
t15 = t566 - t799;
t14 = t87 - t731 / 0.2e1 + t567;
t13 = t88 - t732 / 0.2e1 + t567;
t12 = t27 - t741 / 0.2e1 + t578;
t11 = t27 + t30 - t740 / 0.2e1 + t694;
t10 = t30 - t743 / 0.2e1 + t578;
t9 = t82 - t733 / 0.2e1 + t567;
t8 = t84 - t734 / 0.2e1 + t567;
t7 = t87 + t88 + t506;
t6 = t82 + t84 + t506;
t4 = t508 + t784;
t3 = t508 - t784;
t2 = t498 + t692 + t693;
t1 = (t117 / 0.2e1 - t54 / 0.2e1 + t126 / 0.2e1 - t61 / 0.2e1) * t495 + (t53 / 0.2e1 + t116 / 0.2e1 + t60 / 0.2e1 + t125 / 0.2e1) * t493 + t508;
t5 = [qJD(2) * t22 + qJD(3) * t78 + qJD(4) * t56 + qJD(5) * t86, t22 * qJD(1) + t34 * qJD(3) + t2 * qJD(4) + t6 * qJD(5) + ((t233 * t262 + t234 * t260 + (-t261 - t264) * t785) * t753 + (t210 * t232 + t211 * t230 - t229 * t236 - t231 * t235) * t752 + (t255 * t344 + t256 * t342) * t754 + m(3) * ((-t311 * t495 - t312 * t493) * t439 + (-t408 * t495 + t409 * t493) * t437) / 0.2e1) * t757 + (t498 + (t465 * t592 + t466 * t590 + t492 * t587 + t494 * t585) * t728 + (t54 + t61) * t786 + (t529 + t527) * (t490 / 0.2e1 + t489 / 0.2e1) + (t53 + t116 + t60 + t125) * t729 + (-t465 * t593 + t466 * t591 - t492 * t588 + t494 * t586 + t117 + t126) * t726) * qJD(2), qJD(1) * t78 + qJD(2) * t34 + qJD(4) * t74 + qJD(5) * t150, t56 * qJD(1) + t2 * qJD(2) + t74 * qJD(3) + t498 * qJD(4) + t7 * qJD(5) + ((t616 + (-t270 + t272) * t534) * t752 + (-t511 + (-t371 * t495 + t372 * t493) * t400) * t753) * t755, t86 * qJD(1) + t6 * qJD(2) + t150 * qJD(3) + t7 * qJD(4) + ((-t512 + (-t358 * t495 + t359 * t493) * t393) * m(6) + t506) * qJD(5); (t496 - (t414 + t413) * t466 / 0.2e1 - (t435 + t434) * t494 / 0.2e1 + t782) * qJD(1) + t1 * qJD(2) + t33 * qJD(3) + t3 * qJD(4) + t9 * qJD(5) + (-t724 / 0.4e1 - t718 / 0.4e1 - t723 / 0.4e1 - t710 / 0.4e1) * t758, t1 * qJD(1) + (m(6) * (t113 * t80 - t229 * t230 - t231 * t232) + m(5) * (t114 * t141 - t260 * t785 - t261 * t262) + m(4) * (-t777 * t344 - t778 * t342 + (t493 * (rSges(4,1) * t635 - t582) + t495 * (rSges(4,1) * t634 + t563) + t600) * (-t580 * t416 - t540)) + m(3) * ((t493 * (rSges(3,1) * t625 - t581) + t495 * (rSges(3,1) * t624 + t493 * rSges(3,3) - t456)) * (-t493 * t408 - t409 * t495) + t580 * t439 * t437) + t566 + ((-t770 * t493 + t761 * t495 + t762) * t495 + t769 * t489) * t728 + (((t761 - t769) * t495 + t762) * t493 + t770 * t490) * t726) * qJD(2) + t15 * qJD(4) + t20 * qJD(5), t798, t3 * qJD(1) + t15 * qJD(2) + t11 * qJD(5) + ((t576 + t55) * t752 + (t67 + t120) * t753) * t755 + (t566 - t119 - t737) * qJD(4), t9 * qJD(1) + t20 * qJD(2) + t11 * qJD(4) + (m(6) * (t58 + t112) + t577) * qJD(5); -t33 * qJD(2) + t73 * qJD(4) + t151 * qJD(5) + (-t706 / 0.4e1 - t716 / 0.4e1 - t722 / 0.4e1) * t758, -t798 + t181 + ((-t230 * t495 + t493 * t232) * t752 + (-t260 * t495 + t493 * t262) * t753 + (-t342 * t495 + t493 * t344) * t754) * t757, 0, qJD(2) * t701 + t181 + t797, t579; t496 * qJD(1) + t4 * qJD(2) - t73 * qJD(3) + t508 * qJD(4) + t14 * qJD(5) + (-t708 / 0.4e1 - t717 / 0.4e1) * t758, t4 * qJD(1) + t16 * qJD(4) + t12 * qJD(5) + ((t113 * t121 + t133 * t80 - t230 * t272 - t232 * t534 + t615) * t752 + (t213 * t141 + (-t260 * t493 - t262 * t495) * t400 + t67) * t753) * t757 + (t566 + t799) * qJD(2), -t797, t508 * qJD(1) + t16 * qJD(2) + (m(6) * t576 + t537) * qJD(4) + t21 * qJD(5), t14 * qJD(1) + t12 * qJD(2) + t21 * qJD(4) + (m(6) * (t72 + t112) + t577) * qJD(5); (t503 - t707) * qJD(1) + t8 * qJD(2) - t151 * qJD(3) + t13 * qJD(4) + t567 * qJD(5), t8 * qJD(1) + ((t202 * t113 + (-t230 * t493 - t232 * t495) * t393) * m(6) + t694) * qJD(2) + t10 * qJD(4) + t23, -t579, t13 * qJD(1) + t10 * qJD(2) + ((t132 - t72 + t497) * m(6) + t694) * qJD(4) + t23, qJD(1) * t567 + t23 + (qJD(2) + qJD(4)) * t24;];
Cq = t5;
