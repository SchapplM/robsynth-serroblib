% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:37
% EndTime: 2019-12-31 19:32:19
% DurationCPUTime: 29.18s
% Computational Cost: add. (52576->806), mult. (52007->1134), div. (0->0), fcn. (55552->10), ass. (0->469)
t507 = sin(qJ(1));
t500 = qJ(2) + pkin(8);
t486 = sin(t500);
t488 = cos(t500);
t499 = pkin(9) + qJ(5);
t485 = sin(t499);
t487 = cos(t499);
t724 = rSges(6,1) * t487;
t569 = -rSges(6,2) * t485 + t724;
t349 = -rSges(6,3) * t488 + t486 * t569;
t693 = qJ(4) * t488;
t732 = pkin(3) * t486;
t506 = sin(qJ(2));
t734 = pkin(2) * t506;
t594 = -t693 + t732 + t734;
t505 = -pkin(7) - qJ(4);
t653 = qJ(4) + t505;
t504 = cos(pkin(9));
t481 = pkin(4) * t504 + pkin(3);
t730 = -pkin(3) + t481;
t543 = t730 * t486 + t653 * t488 + t349 + t594;
t220 = t543 * t507;
t509 = cos(qJ(1));
t222 = t543 * t509;
t503 = sin(pkin(9));
t725 = rSges(5,1) * t504;
t572 = -rSges(5,2) * t503 + t725;
t529 = t572 * t486;
t837 = rSges(5,3) * t488 - t529;
t577 = t594 - t837;
t292 = t577 * t507;
t294 = t577 * t509;
t729 = -qJ(3) - pkin(6);
t482 = t507 * t729;
t508 = cos(qJ(2));
t733 = pkin(2) * t508;
t484 = pkin(1) + t733;
t655 = t509 * t503;
t661 = t507 * t504;
t428 = -t488 * t655 + t661;
t662 = t507 * t503;
t667 = t504 * t509;
t429 = t488 * t667 + t662;
t573 = t429 * rSges(5,1) + t428 * rSges(5,2);
t717 = rSges(5,3) + qJ(4);
t731 = pkin(3) * t488;
t804 = -t717 * t486 - t731;
t242 = -t482 + (t484 - t804) * t509 + t573;
t668 = t488 * t509;
t669 = t488 * t507;
t617 = -t507 * t484 - t509 * t729;
t426 = t488 * t662 + t667;
t427 = t488 * t661 - t655;
t815 = -t427 * rSges(5,1) + t426 * rSges(5,2);
t839 = t507 * t804 + t617 + t815;
t646 = t242 * t669 + t668 * t839;
t656 = t509 * t485;
t663 = t507 * t487;
t404 = -t488 * t656 + t663;
t664 = t507 * t485;
t405 = t487 * t668 + t664;
t570 = t405 * rSges(6,1) + t404 * rSges(6,2);
t675 = t486 * t509;
t578 = pkin(4) * t662 - t505 * t675;
t680 = t481 * t488;
t722 = rSges(6,3) * t486;
t210 = -t482 + (t484 + t680 + t722) * t509 + t570 + t578;
t579 = -pkin(4) * t655 + t481 * t669;
t676 = t486 * t507;
t716 = -rSges(6,3) + t505;
t672 = t487 * t509;
t402 = t488 * t664 + t672;
t403 = t488 * t663 - t656;
t816 = -t403 * rSges(6,1) + t402 * rSges(6,2);
t840 = t716 * t676 - t579 + t617 + t816;
t651 = t210 * t669 + t668 * t840;
t793 = m(6) / 0.2e1;
t795 = m(5) / 0.2e1;
t713 = (-t292 * t675 + t294 * t676 + t646) * t795 + (-t220 * t675 + t222 * t676 + t651) * t793;
t253 = (t734 + t716 * t488 + (t481 + t569) * t486) * t507;
t620 = t486 * rSges(6,2) * t656 + rSges(6,3) * t668;
t670 = t488 * t505;
t254 = (-t734 - t670 + (-t481 - t724) * t486) * t509 + t620;
t471 = pkin(3) * t676;
t271 = t471 + (-t717 * t488 + t529 + t734) * t507;
t463 = qJ(4) * t668;
t618 = t486 * rSges(5,2) * t655 + rSges(5,3) * t668;
t272 = t463 + (-t734 + (-pkin(3) - t725) * t486) * t509 + t618;
t714 = ((t271 * t509 + t272 * t507) * t486 + t646) * t795 + ((t253 * t509 + t254 * t507) * t486 + t651) * t793;
t7 = t714 - t713;
t858 = t7 * qJD(1);
t556 = Icges(6,5) * t487 - Icges(6,6) * t485;
t343 = -Icges(6,3) * t488 + t486 * t556;
t702 = Icges(6,4) * t487;
t562 = -Icges(6,2) * t485 + t702;
t345 = -Icges(6,6) * t488 + t486 * t562;
t703 = Icges(6,4) * t485;
t564 = Icges(6,1) * t487 - t703;
t347 = -Icges(6,5) * t488 + t486 * t564;
t178 = t343 * t676 - t345 * t402 + t347 * t403;
t277 = Icges(6,5) * t403 - Icges(6,6) * t402 + Icges(6,3) * t676;
t388 = Icges(6,4) * t403;
t280 = -Icges(6,2) * t402 + Icges(6,6) * t676 + t388;
t387 = Icges(6,4) * t402;
t284 = -Icges(6,1) * t403 - Icges(6,5) * t676 + t387;
t137 = t277 * t676 - t280 * t402 - t284 * t403;
t279 = Icges(6,5) * t405 + Icges(6,6) * t404 + Icges(6,3) * t675;
t704 = Icges(6,4) * t405;
t282 = Icges(6,2) * t404 + Icges(6,6) * t675 + t704;
t389 = Icges(6,4) * t404;
t285 = Icges(6,1) * t405 + Icges(6,5) * t675 + t389;
t138 = t279 * t676 - t402 * t282 + t403 * t285;
t555 = t137 * t507 + t138 * t509;
t17 = t178 * t488 - t486 * t555;
t180 = t343 * t675 + t404 * t345 + t405 * t347;
t139 = t277 * t675 + t404 * t280 - t405 * t284;
t140 = t279 * t675 + t404 * t282 + t405 * t285;
t554 = t507 * t139 + t140 * t509;
t856 = -t180 * t488 + t486 * t554;
t297 = Icges(5,5) * t427 - Icges(5,6) * t426 + Icges(5,3) * t676;
t300 = Icges(5,4) * t427 - Icges(5,2) * t426 + Icges(5,6) * t676;
t303 = Icges(5,1) * t427 - Icges(5,4) * t426 + Icges(5,5) * t676;
t834 = t428 * t300 + t429 * t303;
t161 = t297 * t675 + t834;
t299 = Icges(5,5) * t429 + Icges(5,6) * t428 + Icges(5,3) * t675;
t302 = Icges(5,4) * t429 + Icges(5,2) * t428 + Icges(5,6) * t675;
t305 = Icges(5,1) * t429 + Icges(5,4) * t428 + Icges(5,5) * t675;
t852 = -t139 * t509 + t140 * t507;
t855 = (t299 * t675 + t428 * t302 + t429 * t305) * t507 - t161 * t509 + t852;
t73 = -t137 * t509 + t138 * t507;
t289 = rSges(6,3) * t675 + t570;
t218 = t488 * t289 + t349 * t675;
t287 = rSges(6,3) * t676 - t816;
t851 = t287 * t488 + t349 * t676;
t818 = t218 * t507 - t509 * t851;
t687 = t277 * t488;
t848 = t280 * t485 + t284 * t487;
t155 = t848 * t486 + t687;
t796 = m(4) / 0.2e1;
t447 = rSges(4,1) * t486 + rSges(4,2) * t488;
t526 = t447 + t734;
t822 = t526 * t509;
t823 = t526 * t507;
t833 = t507 * t823 + t509 * t822;
t602 = (t507 * t253 - t254 * t509) * t793 + (t507 * t271 - t272 * t509) * t795 + t833 * t796;
t286 = t507 * t292;
t219 = t507 * t220;
t817 = t222 * t509 + t219;
t831 = -m(6) / 0.2e1;
t603 = t817 * t831 + (-t294 * t509 - t286) * t795 - m(4) * t833 / 0.2e1;
t32 = t603 - t602;
t854 = t32 * qJD(1);
t501 = t507 ^ 2;
t502 = t509 ^ 2;
t615 = t501 + t502;
t705 = Icges(4,4) * t486;
t443 = Icges(4,1) * t488 - t705;
t383 = Icges(4,5) * t507 + t443 * t509;
t440 = Icges(4,2) * t488 + t705;
t695 = Icges(5,3) * t488;
t557 = Icges(5,5) * t504 - Icges(5,6) * t503;
t825 = t486 * t557;
t519 = t695 - t825;
t850 = t302 * t503 - t305 * t504 - t383 + (t440 + t519) * t509;
t849 = -Icges(3,5) * t506 - Icges(4,5) * t486 - Icges(3,6) * t508 - Icges(4,6) * t488;
t819 = -t210 * t507 - t509 * t840;
t773 = -t488 / 0.2e1;
t666 = t506 * t507;
t660 = t507 * t508;
t344 = Icges(6,3) * t486 + t488 * t556;
t706 = Icges(3,4) * t506;
t455 = Icges(3,2) * t508 + t706;
t458 = Icges(3,1) * t508 - t706;
t479 = Icges(4,4) * t488;
t674 = t487 * t347;
t679 = t485 * t345;
t698 = Icges(4,2) * t486;
t707 = Icges(4,1) * t486;
t771 = t504 / 0.2e1;
t772 = -t503 / 0.2e1;
t565 = Icges(5,1) * t504 - Icges(5,4) * t503;
t805 = -Icges(5,5) * t488 + t486 * t565;
t563 = Icges(5,4) * t504 - Icges(5,2) * t503;
t806 = -Icges(5,6) * t488 + t486 * t563;
t835 = -(t458 / 0.2e1 - t455 / 0.2e1) * t506 - (t805 * t771 + t806 * t772 + t674 / 0.2e1 - t679 / 0.2e1 + t479 + t707 / 0.2e1 - t698 / 0.2e1 - Icges(5,3) * t486 / 0.2e1 + t557 * t773 - t344 / 0.2e1) * t488;
t549 = -t300 * t426 + t303 * t427;
t430 = t615 * t486;
t832 = -0.2e1 * t430;
t826 = t349 * t507;
t314 = -rSges(6,1) * t402 - rSges(6,2) * t403;
t315 = rSges(6,1) * t404 - rSges(6,2) * t405;
t194 = (t314 * t509 - t315 * t507) * t486;
t611 = m(5) / 0.4e1 + m(6) / 0.4e1;
t677 = t486 * t488;
t619 = t615 * t677;
t821 = t611 * (t619 - t677);
t820 = t653 * t486;
t814 = t849 * t507;
t813 = t849 * t509;
t495 = Icges(3,4) * t508;
t456 = -Icges(3,2) * t506 + t495;
t457 = Icges(3,1) * t506 + t495;
t466 = Icges(4,4) * t676;
t382 = Icges(4,1) * t669 - Icges(4,5) * t509 - t466;
t812 = -Icges(4,2) * t669 - t300 * t503 + t303 * t504 - t519 * t507 + t382 - t466;
t766 = -t509 / 0.2e1;
t811 = qJD(2) * t766;
t769 = t507 / 0.2e1;
t810 = qJD(2) * t769;
t809 = t299 * t676 - t426 * t302 + t427 * t305;
t385 = (-Icges(6,2) * t487 - t703) * t486;
t386 = (-Icges(6,1) * t485 - t702) * t486;
t803 = -t485 * (t347 / 0.2e1 + t385 / 0.2e1) + t487 * (t386 / 0.2e1 - t345 / 0.2e1);
t476 = Icges(3,4) * t666;
t416 = Icges(3,1) * t660 - Icges(3,5) * t509 - t476;
t622 = -Icges(3,2) * t660 + t416 - t476;
t414 = Icges(3,4) * t660 - Icges(3,2) * t666 - Icges(3,6) * t509;
t624 = t457 * t507 + t414;
t380 = Icges(4,4) * t669 - Icges(4,2) * t676 - Icges(4,6) * t509;
t566 = -t479 - t707;
t630 = -t566 * t507 + t380;
t802 = t630 * t488 + t622 * t506 + t624 * t508;
t417 = Icges(3,5) * t507 + t458 * t509;
t621 = -t455 * t509 + t417;
t415 = Icges(3,6) * t507 + t456 * t509;
t623 = -t457 * t509 - t415;
t381 = Icges(4,6) * t507 + (t479 - t698) * t509;
t629 = t566 * t509 - t381;
t801 = (-t297 * t509 + t299 * t507) * t488 + t623 * t660 + t629 * t669 + (t850 * t507 + t509 * t812) * t486 - t621 * t666;
t799 = 0.4e1 * qJD(1);
t798 = 0.2e1 * qJD(2);
t791 = -t856 / 0.2e1;
t790 = t73 / 0.2e1;
t789 = t852 / 0.2e1;
t330 = -rSges(6,1) * t486 * t672 + t620;
t550 = t287 * t509 - t289 * t507;
t143 = t550 * t488 + (-t330 * t507 - t509 * t826) * t486;
t350 = t488 * t569 + t722;
t165 = (t350 * t507 - t287) * t486;
t166 = (-t349 * t509 - t330) * t488 + (-t350 * t509 + t289) * t486;
t183 = t550 * t486;
t650 = -t218 * t669 + t668 * t851;
t785 = m(6) * (-t143 * t488 + (t165 * t509 + t166 * t507 + t183) * t486 + t650);
t784 = m(6) * (t165 * t840 + t166 * t210 - t218 * t254 + t253 * t851);
t448 = qJ(4) * t486 + t731;
t498 = t509 * pkin(6);
t633 = -t507 * (pkin(1) * t507 - t498 + t617) + t509 * (-t507 * pkin(6) - t482 + (-pkin(1) + t484) * t509);
t581 = t615 * t448 + t633;
t113 = (t287 - (t731 + t820) * t507 + t579) * t507 + (t289 + (-t448 + t680) * t509 + t578) * t509 + t581;
t648 = -t488 * t219 - t222 * t668;
t782 = m(6) * (t113 * t430 + t648);
t390 = (-rSges(6,1) * t485 - rSges(6,2) * t487) * t486;
t781 = m(6) * (-t220 * t315 + t222 * t314 + t819 * t390);
t308 = -Icges(6,5) * t402 - Icges(6,6) * t403;
t641 = -Icges(6,2) * t403 - t284 - t387;
t643 = -Icges(6,1) * t402 - t280 - t388;
t121 = -t308 * t488 + (-t641 * t485 + t643 * t487) * t486;
t776 = t121 / 0.2e1;
t774 = t486 / 0.2e1;
t768 = t507 / 0.4e1;
t765 = -t509 / 0.4e1;
t764 = t509 / 0.2e1;
t727 = rSges(3,1) * t508;
t595 = pkin(1) + t727;
t616 = rSges(3,2) * t666 + t509 * rSges(3,3);
t356 = -t507 * t595 + t498 + t616;
t665 = t506 * t509;
t478 = rSges(3,2) * t665;
t357 = -t478 + t595 * t509 + (rSges(3,3) + pkin(6)) * t507;
t459 = rSges(3,1) * t506 + rSges(3,2) * t508;
t437 = t459 * t507;
t438 = t459 * t509;
t763 = m(3) * (t356 * t437 - t357 * t438);
t540 = rSges(4,1) * t669 - rSges(4,2) * t676 - t509 * rSges(4,3);
t320 = -t540 + t617;
t591 = -rSges(4,2) * t675 + t507 * rSges(4,3);
t726 = rSges(4,1) * t488;
t321 = -t482 + (t484 + t726) * t509 + t591;
t761 = m(4) * (t320 * t823 - t321 * t822);
t760 = m(4) * (t320 * t509 + t321 * t507);
t150 = t507 * (rSges(5,3) * t676 - t815) + t509 * (rSges(5,3) * t675 + t573) + t581;
t644 = -t488 * t286 - t294 * t668;
t756 = m(5) * (t150 * t430 + t644);
t752 = m(5) * (t242 * t272 + t271 * t839);
t230 = t242 * t675;
t751 = m(5) * (-t676 * t839 + t230);
t750 = m(5) * (t242 * t507 + t509 * t839);
t746 = m(6) * (t183 * t430 + t650);
t745 = m(6) * (t210 * t254 + t253 * t840);
t198 = t210 * t675;
t744 = m(6) * (-t676 * t840 + t198);
t743 = m(6) * (-t218 * t675 - t676 * t851);
t742 = m(6) * t819;
t741 = m(6) * t818;
t214 = t507 * t314 + t315 * t509;
t738 = m(6) * (-t214 * t488 - t390 * t430);
t736 = m(6) * t194;
t735 = m(6) * t214;
t728 = m(6) * qJD(5);
t719 = t507 * t17;
t718 = t509 * t856;
t324 = t345 * t507;
t326 = t347 * t507;
t535 = -t343 * t507 + t848;
t109 = -t535 * t488 + (t324 * t485 - t326 * t487 + t277) * t486;
t692 = t109 * t509;
t325 = t345 * t509;
t327 = t347 * t509;
t551 = -t282 * t485 + t285 * t487;
t534 = -t343 * t509 - t551;
t110 = -t534 * t488 + (t325 * t485 - t327 * t487 + t279) * t486;
t691 = t110 * t507;
t686 = t279 * t488;
t685 = t343 * t488;
t683 = t380 * t486;
t681 = t414 * t506;
t346 = Icges(6,6) * t486 + t488 * t562;
t678 = t485 * t346;
t348 = Icges(6,5) * t486 + t488 * t564;
t673 = t487 * t348;
t384 = (-Icges(6,5) * t485 - Icges(6,6) * t487) * t486;
t671 = t488 * t384;
t659 = t508 * t509;
t98 = t165 * t507 - t166 * t509;
t654 = t98 * qJD(3);
t642 = Icges(6,1) * t404 - t282 - t704;
t640 = -Icges(6,2) * t405 + t285 + t389;
t378 = Icges(4,5) * t669 - Icges(4,6) * t676 - Icges(4,3) * t509;
t638 = -t507 * t378 - t382 * t668;
t559 = Icges(4,5) * t488 - Icges(4,6) * t486;
t379 = Icges(4,3) * t507 + t509 * t559;
t637 = t507 * t379 + t383 * t668;
t636 = -t345 + t386;
t635 = t347 + t385;
t412 = Icges(3,5) * t660 - Icges(3,6) * t666 - Icges(3,3) * t509;
t632 = -t507 * t412 - t416 * t659;
t561 = Icges(3,5) * t508 - Icges(3,6) * t506;
t413 = Icges(3,3) * t507 + t509 * t561;
t631 = t507 * t413 + t417 * t659;
t614 = qJD(1) * t486;
t613 = qJD(5) * t486;
t610 = t793 + t795;
t237 = t610 * t832;
t612 = t237 * qJD(1);
t609 = t98 * t793;
t605 = t791 + t856 / 0.2e1;
t601 = t676 / 0.4e1;
t593 = -t448 - t733;
t592 = rSges(4,2) * t486 - t726 - t733;
t340 = t383 * t669;
t586 = t379 * t509 - t340;
t370 = t417 * t660;
t585 = t413 * t509 - t370;
t584 = t381 * t486 - t378;
t583 = t415 * t506 - t412;
t580 = t615 * t734;
t576 = -rSges(5,3) * t486 - t488 * t572 + t593;
t309 = Icges(6,5) * t404 - Icges(6,6) * t405;
t122 = -t309 * t488 + (-t640 * t485 + t642 * t487) * t486;
t135 = t384 * t676 - t635 * t402 + t636 * t403;
t136 = t384 * t675 + t635 * t404 + t636 * t405;
t575 = t781 / 0.2e1 + (t122 + t136) * t768 + (t121 + t135) * t765;
t156 = t486 * t551 - t686;
t553 = -t155 * t507 + t156 * t509;
t547 = t674 - t679;
t544 = -t481 * t486 - t670;
t542 = -t730 * t488 - t350 + t593 + t820;
t541 = m(6) * (t210 * t315 - t314 * t840) - t671 / 0.2e1;
t101 = t308 * t676 - t641 * t402 + t643 * t403;
t102 = t309 * t676 - t640 * t402 + t642 * t403;
t50 = -t101 * t509 + t102 * t507;
t103 = t308 * t675 + t641 * t404 + t643 * t405;
t104 = t309 * t675 + t640 * t404 + t642 * t405;
t51 = -t103 * t509 + t104 * t507;
t538 = t50 * t766 + t51 * t769;
t531 = t344 - t547;
t527 = t507 * (qJ(4) * t669 - t471) + t509 * (-pkin(3) * t675 + t463) - t580;
t524 = t17 * t768 + t856 * t765 - t719 / 0.4e1 + t718 / 0.4e1 + (t601 - t676 / 0.4e1) * t852;
t516 = t486 * t535 + t687;
t515 = t486 * t534 + t686;
t514 = t486 * t531 + t685;
t128 = -t346 * t402 + t348 * t403 + t507 * t514;
t129 = t404 * t346 + t405 * t348 + t509 * t514;
t131 = -t531 * t488 + (t343 + t673 - t678) * t486;
t190 = t486 * t547 - t685;
t513 = t131 * t773 + t190 * t774 + t784 / 0.2e1 + (t109 + t128) * t601 + (t110 + t129) * t675 / 0.4e1 + (-t155 + t178) * t669 / 0.4e1 + (t156 + t180) * t668 / 0.4e1;
t363 = Icges(5,6) * t486 + t488 * t563;
t365 = Icges(5,5) * t486 + t488 * t565;
t510 = t365 * t771 + t363 * t772 + t673 / 0.2e1 - t678 / 0.2e1 + t443 / 0.2e1 - t440 / 0.2e1 + t825 / 0.2e1 - t695 / 0.2e1 + t343 / 0.2e1;
t461 = -rSges(3,2) * t506 + t727;
t376 = t592 * t509;
t374 = t592 * t507;
t337 = t805 * t509;
t336 = t805 * t507;
t335 = t806 * t509;
t334 = t806 * t507;
t295 = t576 * t509;
t293 = t576 * t507;
t252 = -t415 * t665 + t631;
t251 = -t414 * t665 - t632;
t250 = -t415 * t666 - t585;
t239 = 0.4e1 * t821;
t236 = t611 * t832 + (m(5) + m(6)) * t430 / 0.2e1;
t234 = -t488 * t315 - t390 * t675;
t233 = t314 * t488 + t390 * t676;
t227 = -t381 * t675 + t637;
t226 = -t380 * t675 - t638;
t225 = -t381 * t676 - t586;
t223 = t542 * t509;
t221 = t542 * t507;
t203 = -t735 / 0.2e1;
t191 = -t736 / 0.2e1;
t181 = t509 * (-rSges(5,1) * t486 * t667 + t618) + t837 * t501 + t527;
t175 = -t251 * t509 + t252 * t507;
t174 = -(-t507 * (-t416 * t508 + t681) - t412 * t509) * t509 + t250 * t507;
t169 = t738 / 0.2e1;
t168 = -t671 + (-t635 * t485 + t636 * t487) * t486;
t164 = -t226 * t509 + t227 * t507;
t163 = -(-t507 * (-t382 * t488 + t683) - t378 * t509) * t509 + t225 * t507;
t146 = -t741 / 0.2e1;
t142 = t743 / 0.2e1;
t134 = (-t463 + t330 + (t544 + t732) * t509) * t509 + (t471 - t826 + (t544 - t693) * t507) * t507 + t527;
t99 = t746 / 0.2e1;
t97 = -t404 * t325 - t405 * t327 + t509 * t515;
t96 = -t404 * t324 - t405 * t326 + t509 * t516;
t95 = t325 * t402 - t327 * t403 + t507 * t515;
t94 = t324 * t402 - t326 * t403 + t507 * t516;
t93 = qJD(2) * t609;
t89 = -(t297 * t676 + t549) * t509 + t507 * t809;
t85 = t744 + t751;
t76 = (t250 - t370 + (t413 + t681) * t509 + t632) * t509 + t631 * t507;
t75 = (t509 * t583 + t252 - t631) * t509 + (t507 * t583 + t251 + t585) * t507;
t70 = -t742 + t750 + t760;
t69 = t486 * t803 + t541;
t68 = -t190 * t488 + t486 * t553;
t62 = (t225 - t340 + (t379 + t683) * t509 + t638) * t509 + t637 * t507;
t61 = (t509 * t584 + t227 - t637) * t509 + (t507 * t584 + t226 + t586) * t507;
t52 = t113 * t214 + t817 * t390;
t49 = t146 + t735 / 0.2e1;
t48 = t203 + t146;
t47 = t203 + t741 / 0.2e1;
t46 = t97 * t507 - t509 * t96;
t45 = t95 * t507 - t509 * t94;
t44 = t756 + t782;
t42 = t142 + t736 / 0.2e1;
t41 = t191 + t142;
t40 = t191 - t743 / 0.2e1;
t38 = t143 * t183 + t165 * t851 - t166 * t218;
t36 = -t136 * t488 + (t103 * t507 + t104 * t509) * t486;
t35 = -t135 * t488 + (t101 * t507 + t102 * t509) * t486;
t34 = t785 / 0.2e1;
t31 = t602 + t603;
t28 = t549 * t509 + (t161 - t809 - t834) * t507;
t23 = (t457 / 0.2e1 + t456 / 0.2e1) * t508 + t763 + t761 + t752 + t745 + t510 * t486 - t835;
t22 = (-t131 + t553) * t488 + (t109 * t507 + t110 * t509 + t190) * t486;
t21 = t99 + t34 - t738 / 0.2e1;
t20 = t169 + t99 - t785 / 0.2e1;
t19 = t169 + t34 - t746 / 0.2e1;
t14 = (-t129 + t554) * t488 + (t507 * t96 + t509 * t97 + t180) * t486;
t13 = (-t128 + t555) * t488 + (t507 * t94 + t509 * t95 + t178) * t486;
t10 = m(6) * t52 + t538;
t8 = t713 + t714;
t6 = t605 * t676;
t5 = m(6) * t38 + (t718 / 0.2e1 - t719 / 0.2e1 - t22 / 0.2e1) * t488 + (t14 * t764 + t13 * t769 + t68 / 0.2e1) * t486;
t4 = (t164 / 0.2e1 - t62 / 0.2e1 - t76 / 0.2e1 + t175 / 0.2e1 - t852 / 0.2e1 + t789) * t509 + (t28 / 0.2e1 + t61 / 0.2e1 + t163 / 0.2e1 + t174 / 0.2e1 + t75 / 0.2e1 + t89 / 0.2e1 + t790 - t73 / 0.2e1) * t507;
t3 = t524 + (-t190 / 0.2e1 + (-t129 / 0.4e1 - t110 / 0.4e1) * t509 + (-t128 / 0.4e1 - t109 / 0.4e1) * t507) * t486 + (t131 / 0.2e1 + (-t180 / 0.4e1 - t156 / 0.4e1) * t509 + (-t178 / 0.4e1 + t155 / 0.4e1) * t507) * t488 - t784 / 0.2e1 + t575;
t2 = t513 + (t135 / 0.4e1 + t121 / 0.4e1) * t509 + t524 + (-t136 / 0.4e1 - t122 / 0.4e1) * t507 - t781 / 0.2e1;
t1 = t513 + t575;
t9 = [t23 * qJD(2) + t70 * qJD(3) + t85 * qJD(4) + t69 * qJD(5), t23 * qJD(1) + t31 * qJD(3) + t8 * qJD(4) + t1 * qJD(5) + (m(3) * ((-t356 * t509 - t357 * t507) * t461 + (-t437 * t509 + t438 * t507) * t459) / 0.2e1 + (t320 * t376 + t321 * t374) * t796 + (t242 * t293 - t271 * t294 - t272 * t292 + t295 * t839) * t795 + (t210 * t221 - t220 * t254 - t222 * t253 + t223 * t840) * t793) * t798 + (t428 * t363 + t429 * t365 + t506 * t623 + t508 * t621 + t129 - t850 * t488 + (t335 * t503 - t337 * t504 + t299 + t629) * t486) * t810 + (-t363 * t426 + t365 * t427 - t506 * t624 + t508 * t622 + t128 + t164 + t175 + t812 * t488 + (t334 * t503 - t336 * t504 + t297 - t630) * t486 + t855) * t811 + (t691 / 0.2e1 - t692 / 0.2e1 + (t561 + t559) * (t502 / 0.2e1 + t501 / 0.2e1) + (t76 + t62 + t855) * t764 - (t75 + t61 + t28 + t174 + t163 + t89) * t507 / 0.2e1) * qJD(2), qJD(1) * t70 + qJD(2) * t31 + qJD(4) * t236 + qJD(5) * t48, qJD(1) * t85 + qJD(2) * t8 + qJD(3) * t236 + qJD(5) * t41, t69 * qJD(1) + t1 * qJD(2) + t48 * qJD(3) + t41 * qJD(4) + (m(6) * (t210 * t234 - t218 * t315 + t233 * t840 - t314 * t851) - t168 * t488) * qJD(5) + ((t136 / 0.2e1 + t122 / 0.2e1) * t509 + (t135 / 0.2e1 + t776 - t605) * t507) * t613; t4 * qJD(2) + t32 * qJD(3) - t7 * qJD(4) + t3 * qJD(5) + (-t745 / 0.4e1 - t752 / 0.4e1 - t761 / 0.4e1 - t763 / 0.4e1) * t799 - t510 * t614 + (-(t457 + t456) * t508 / 0.2e1 + t835) * qJD(1), t4 * qJD(1) + (m(5) * (t150 * t181 - t292 * t293 - t294 * t295) + m(4) * (-t822 * t376 - t823 * t374 + (t507 * t540 + t509 * (rSges(4,1) * t668 + t591) + t633) * (-t447 * t615 - t580)) + m(3) * ((t507 * (rSges(3,1) * t660 - t616) + t509 * (rSges(3,1) * t659 + t507 * rSges(3,3) - t478)) * (-t507 * t437 - t438 * t509) + t615 * t461 * t459) + m(6) * (t113 * t134 - t220 * t221 - t222 * t223)) * qJD(2) + t44 * qJD(4) + t10 * qJD(5) + ((-t428 * t335 - t429 * t337) * t507 + t46 + (t428 * t334 + t429 * t336 - t507 * t814 + t509 * t802 + t801) * t509 + t813 * t501) * t810 + (-(t426 * t334 - t427 * t336) * t509 + t45 + t814 * t502 + (t426 * t335 - t427 * t337 + (t802 - t813) * t509 + t801) * t507) * t811, t854 - t98 * t728 / 0.2e1, -t858 + t44 * qJD(2) + t20 * qJD(5) + (-0.4e1 * t821 + 0.2e1 * t610 * (-t430 * t488 + t619)) * qJD(4), t3 * qJD(1) + t10 * qJD(2) + t20 * qJD(4) + t654 * t831 + (-t68 / 0.2e1 + (t51 / 0.2e1 - t14 / 0.2e1) * t509 + (t50 / 0.2e1 - t13 / 0.2e1) * t507) * t613 + (t35 * t766 + t36 * t769 + (t22 / 0.2e1 + (t776 + t791) * t509 + (-t122 / 0.2e1 + t17 / 0.2e1) * t507) * t488 + (t194 * t113 + t183 * t214 - t234 * t220 - t233 * t222 + t818 * t390 - t38) * m(6)) * qJD(5); -t32 * qJD(2) + t237 * qJD(4) + t47 * qJD(5) + (t742 / 0.4e1 - t750 / 0.4e1 - t760 / 0.4e1) * t799, -t854 + ((-t221 * t509 + t507 * t223) * t793 + (-t293 * t509 + t507 * t295) * t795 + (-t374 * t509 + t507 * t376) * t796) * t798 + qJD(5) * t609, 0, t612, t47 * qJD(1) + t93 + (t233 * t507 - t234 * t509) * t728; t7 * qJD(2) - t237 * qJD(3) + t40 * qJD(5) + (-t744 / 0.4e1 - t751 / 0.4e1) * t799 + 0.2e1 * (t198 * t793 + t230 * t795 + (-t210 * t793 - t242 * t795) * t675) * qJD(1), t858 + t239 * qJD(4) + t19 * qJD(5) + 0.4e1 * (-t782 / 0.4e1 - t756 / 0.4e1) * qJD(2) + ((-t488 * t134 + t648) * t793 + (-t488 * t181 + t644) * t795 + ((t221 * t507 + t223 * t509 + t113) * t793 + (t293 * t507 + t295 * t509 + t150) * t795) * t486) * t798, -t612, t239 * qJD(2), t40 * qJD(1) + t19 * qJD(2) + (-t194 * t488 + (t233 * t509 + t234 * t507) * t486) * t728; -t541 * qJD(1) + t2 * qJD(2) + t49 * qJD(3) + t42 * qJD(4) + t6 * qJD(5) - t803 * t614, t2 * qJD(1) + (t13 * t766 + t668 * t789 + t46 * t675 / 0.2e1 + (t155 * t509 + t156 * t507) * t774 + (t691 - t692) * t773 + t669 * t790 + t45 * t676 / 0.2e1 + t14 * t769 - t538) * qJD(2) + t21 * qJD(4) + t5 * qJD(5) + ((t113 * t143 + t134 * t183 - t165 * t222 - t166 * t220 - t218 * t221 + t223 * t851 - t52) * qJD(2) + t654 / 0.2e1) * m(6), qJD(1) * t49 + t93, qJD(1) * t42 + qJD(2) * t21, t6 * qJD(1) + t5 * qJD(2) + (m(6) * (t183 * t194 - t218 * t234 + t233 * t851) + t488 ^ 2 * t168 / 0.2e1 + (t36 * t764 + t35 * t769 + (t121 * t507 + t122 * t509) * t773) * t486) * qJD(5);];
Cq = t9;
