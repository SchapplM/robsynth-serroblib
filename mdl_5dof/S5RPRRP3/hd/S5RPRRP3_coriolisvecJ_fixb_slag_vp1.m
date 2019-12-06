% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:50
% DurationCPUTime: 27.14s
% Computational Cost: add. (17917->655), mult. (15561->808), div. (0->0), fcn. (12116->8), ass. (0->372)
t751 = Icges(5,4) + Icges(6,4);
t725 = Icges(5,2) + Icges(6,2);
t752 = Icges(5,1) + Icges(6,1);
t750 = Icges(5,5) + Icges(6,5);
t749 = Icges(5,6) + Icges(6,6);
t349 = qJ(3) + qJ(4);
t342 = cos(t349);
t753 = t751 * t342;
t341 = sin(t349);
t741 = t725 * t341 - t753;
t748 = Icges(5,3) + Icges(6,3);
t348 = qJ(1) + pkin(8);
t340 = cos(t348);
t544 = t340 * t341;
t747 = t751 * t544;
t339 = sin(t348);
t746 = t749 * t339;
t745 = t750 * t339;
t744 = t751 * t341;
t743 = t752 * t341 + t753;
t728 = t741 * t339 + t749 * t340;
t543 = t340 * t342;
t727 = t751 * t543 - t725 * t544 + t746;
t721 = t752 * t543 + t745 - t747;
t719 = -t749 * t341 + t750 * t342;
t736 = t725 * t342 + t744;
t735 = t752 * t342 - t744;
t548 = t339 * t341;
t742 = t751 * t548;
t740 = t750 * t340;
t547 = t339 * t342;
t738 = -t752 * t547 + t740 + t742;
t737 = t750 * t341 + t749 * t342;
t734 = t748 * t339;
t733 = t743 - t741;
t732 = -t735 + t736;
t731 = t743 * t340 + t727;
t730 = -t743 * t339 + t728;
t729 = -t719 * t339 + t748 * t340;
t724 = t750 * t543 - t749 * t544 + t734;
t726 = t736 * t341 - t342 * t743;
t718 = t727 * t341 - t721 * t342;
t720 = t737 * t339;
t347 = qJD(3) + qJD(4);
t717 = t732 * t347;
t716 = t733 * t347;
t715 = t730 * t347 + (t735 * t340 + t745) * qJD(1);
t714 = t731 * t347 + (t735 * t339 - t740) * qJD(1);
t244 = t339 * t347;
t713 = t738 * t347 + t736 * t244 + (t741 * t340 - t746) * qJD(1);
t245 = t340 * t347;
t712 = t728 * qJD(1) - t736 * t245 + t721 * t347;
t711 = t737 * t340;
t710 = t738 * t342;
t709 = t728 * t341;
t678 = t724 * t340 - t721 * t547 + t727 * t548;
t708 = t726 * t339 + t711;
t707 = t726 * t340 - t720;
t706 = t729 * qJD(1);
t705 = -t729 * t339 - t543 * t738;
t704 = t729 * t340 + t728 * t548;
t703 = t718 * t340;
t679 = -t547 * t738 + t704;
t677 = -t544 * t728 - t705;
t676 = t724 * t339 - t703;
t701 = (t725 * t547 + t738 + t742) * t245 + (-t725 * t543 + t721 - t747) * t244 + t733 * qJD(1);
t700 = t726 * qJD(1) + t719 * t347;
t699 = t720 * t347 + (-t719 * t340 + t709 - t710 - t734) * qJD(1);
t698 = t718 * qJD(1) - t711 * t347 + t706;
t697 = -t724 * qJD(1) + t712 * t341 + t714 * t342;
t696 = t713 * t341 + t715 * t342 - t706;
t695 = -t737 * qJD(1) + t716 * t341 + t717 * t342;
t665 = -t732 * qJD(1) - t731 * t244 - t730 * t245;
t694 = -t701 * t341 + t342 * t665;
t352 = cos(qJ(3));
t345 = t352 * pkin(3);
t335 = t345 + pkin(2);
t601 = pkin(4) * t342;
t279 = t335 + t601;
t692 = rSges(6,2) * t544 - t340 * t279;
t691 = t707 * qJD(1);
t690 = t708 * qJD(1) + t678 * t244;
t689 = t699 * t339 - t696 * t340;
t688 = t698 * t339 - t697 * t340;
t687 = t696 * t339 + t699 * t340;
t686 = t697 * t339 + t698 * t340;
t685 = t679 * t245 + t690;
t684 = t676 * t244 + t677 * t245 - t691;
t683 = -t715 * t341 + t713 * t342;
t682 = -t714 * t341 + t712 * t342;
t681 = t700 * t339 - t695 * t340;
t680 = t695 * t339 + t700 * t340;
t675 = -t341 * t738 - t342 * t728;
t674 = t721 * t341 + t727 * t342;
t265 = rSges(6,1) * t341 + rSges(6,2) * t342;
t210 = t265 * t340;
t350 = sin(qJ(3));
t482 = qJD(3) * t350;
t475 = pkin(3) * t482;
t540 = t341 * t347;
t241 = -pkin(4) * t540 - t475;
t591 = rSges(6,1) * t342;
t267 = -rSges(6,2) * t341 + t591;
t327 = qJD(5) * t339;
t330 = t340 * rSges(6,3);
t354 = -pkin(7) - pkin(6);
t346 = -qJ(5) + t354;
t486 = qJD(1) * t346;
t488 = qJD(1) * t339;
t485 = qJD(1) * t354;
t509 = t335 * t488 + t340 * t485;
t585 = -t279 * t488 + t327 + (t241 + t475 - t486) * t340 + t509 - t347 * t210 + (-t267 * t339 + t330) * qJD(1);
t290 = t340 * t335;
t394 = -rSges(6,1) * t543 - rSges(6,3) * t339;
t497 = t346 - t354;
t535 = t339 * t497 + t290 + t394 + t692;
t667 = t267 + t601;
t666 = (-t724 - t710) * t339 + t703 + t704;
t664 = t719 * qJD(1) - t711 * t244 + t720 * t245;
t417 = rSges(3,1) * t339 + rSges(3,2) * t340;
t351 = sin(qJ(1));
t607 = pkin(1) * t351;
t663 = t417 + t607;
t487 = qJD(1) * t340;
t465 = t341 * t487;
t539 = t342 * t347;
t472 = t339 * t539;
t662 = t465 + t472;
t331 = t340 * rSges(5,3);
t506 = rSges(5,2) * t548 + t331;
t661 = t506 - t607;
t343 = Icges(4,4) * t352;
t412 = -Icges(4,2) * t350 + t343;
t640 = Icges(4,1) * t350 + t343;
t499 = t640 + t412;
t580 = Icges(4,4) * t350;
t308 = Icges(4,2) * t352 + t580;
t311 = Icges(4,1) * t352 - t580;
t500 = t308 - t311;
t657 = (t350 * t499 + t352 * t500) * qJD(1);
t656 = 0.2e1 * qJD(3);
t592 = rSges(5,1) * t342;
t268 = -rSges(5,2) * t341 + t592;
t266 = rSges(5,1) * t341 + rSges(5,2) * t342;
t552 = t266 * t340;
t114 = -t347 * t552 + (-t268 * t339 + t331) * qJD(1);
t209 = t266 * t339;
t653 = t340 * t114 + t209 * t244 + t245 * t552;
t596 = pkin(6) + t354;
t652 = t340 * t596;
t541 = t340 * t352;
t542 = t340 * t350;
t568 = Icges(4,6) * t339;
t187 = Icges(4,4) * t541 - Icges(4,2) * t542 + t568;
t321 = Icges(4,4) * t542;
t577 = Icges(4,5) * t339;
t189 = Icges(4,1) * t541 - t321 + t577;
t406 = t187 * t350 - t189 * t352;
t649 = t406 * t340;
t326 = pkin(2) * t488;
t134 = t326 + (-pkin(6) * qJD(1) - t475) * t340 - t509;
t338 = t340 ^ 2;
t646 = t340 * t134 - (-t339 ^ 2 - t338) * t475;
t604 = pkin(2) * t340;
t161 = t339 * t596 - t290 + t604;
t600 = pkin(6) * t339;
t252 = t600 + t604;
t243 = qJD(1) * t252;
t645 = -qJD(1) * t161 + t243;
t353 = cos(qJ(1));
t606 = pkin(1) * t353;
t644 = -t290 - t606;
t420 = -rSges(3,1) * t340 - t606;
t643 = t339 * rSges(3,2) + t420;
t296 = pkin(4) * t548;
t275 = qJD(1) * t296;
t549 = t339 * t265;
t639 = qJD(1) * t549 - t245 * t667 - t265 * t488 - t275;
t481 = qJD(3) * t352;
t474 = pkin(3) * t481;
t638 = t340 * t474;
t221 = t267 * t347;
t276 = pkin(4) * t465;
t637 = pkin(4) * t472 + t339 * t221 - t244 * t667 + t265 * t487 + t276;
t636 = qJD(1) * t209 - t245 * t268 - t266 * t488;
t328 = qJD(5) * t340;
t473 = t339 * t540;
t635 = rSges(6,1) * t473 + rSges(6,2) * t662 + t339 * t486 + t328;
t323 = rSges(4,2) * t542;
t396 = -rSges(4,1) * t541 - rSges(4,3) * t339;
t192 = -t323 - t396;
t312 = rSges(4,1) * t350 + rSges(4,2) * t352;
t484 = qJD(3) * t339;
t239 = t312 * t484;
t451 = -t252 - t606;
t100 = t239 + (-t192 + t451) * qJD(1);
t332 = t340 * rSges(4,3);
t546 = t339 * t350;
t498 = rSges(4,2) * t546 + t332;
t545 = t339 * t352;
t399 = rSges(4,1) * t545 - t498;
t476 = qJD(1) * t607;
t422 = -qJD(1) * (-pkin(2) * t339 + t340 * pkin(6)) + t476;
t483 = qJD(3) * t340;
t463 = t312 * t483;
t99 = qJD(1) * t399 + t422 + t463;
t634 = t100 * t339 + t340 * t99;
t633 = t245 * t210 + t340 * t585;
t462 = t339 * t482;
t301 = pkin(3) * t462;
t602 = pkin(4) * t341;
t450 = t265 + t602;
t632 = t450 * t244 + t301 + t328;
t507 = rSges(6,2) * t548 + t330;
t393 = t340 * t497 - t507;
t186 = Icges(4,6) * t340 - t339 * t412;
t320 = Icges(4,4) * t546;
t576 = Icges(4,5) * t340;
t188 = -Icges(4,1) * t545 + t320 + t576;
t117 = t186 * t352 + t188 * t350;
t126 = t308 * t484 + (-t340 * t412 - t568) * qJD(1);
t229 = t640 * t339;
t128 = qJD(3) * t229 + (-t311 * t340 - t577) * qJD(1);
t307 = Icges(4,5) * t352 - Icges(4,6) * t350;
t184 = Icges(4,3) * t340 - t307 * t339;
t492 = qJD(1) * t184;
t631 = qJD(3) * t117 + t126 * t350 - t128 * t352 - t492;
t118 = t187 * t352 + t189 * t350;
t125 = qJD(1) * t186 - t308 * t483;
t230 = t640 * t340;
t127 = -qJD(3) * t230 + (-t311 * t339 + t576) * qJD(1);
t565 = Icges(4,3) * t339;
t185 = Icges(4,5) * t541 - Icges(4,6) * t542 + t565;
t630 = -qJD(1) * t185 + qJD(3) * t118 + t125 * t350 - t127 * t352;
t271 = t412 * qJD(3);
t272 = t311 * qJD(3);
t306 = Icges(4,5) * t350 + Icges(4,6) * t352;
t629 = -qJD(1) * t306 + qJD(3) * (t308 * t352 + t350 * t640) + t271 * t350 - t272 * t352;
t222 = t268 * t347;
t294 = rSges(5,2) * t544;
t395 = rSges(5,1) * t543 + rSges(5,3) * t339;
t180 = -t294 + t395;
t425 = t161 + t451;
t432 = t244 * t266 + t301;
t68 = (-t180 + t425) * qJD(1) + t432;
t628 = (-qJD(1) * t552 + t222 * t339 - t244 * t268 + t266 * t487) * t68;
t515 = -Icges(4,2) * t541 + t189 - t321;
t517 = t187 + t230;
t626 = t350 * t515 + t352 * t517;
t516 = Icges(4,2) * t545 + t188 + t320;
t518 = t186 - t229;
t625 = -t350 * t516 - t352 * t518;
t356 = qJD(1) ^ 2;
t235 = t347 * t488;
t616 = -t235 / 0.2e1;
t236 = qJD(1) * t245;
t615 = t236 / 0.2e1;
t614 = -t244 / 0.2e1;
t613 = t244 / 0.2e1;
t612 = -t245 / 0.2e1;
t611 = t245 / 0.2e1;
t610 = t339 / 0.2e1;
t609 = t340 / 0.2e1;
t608 = -rSges(4,3) - pkin(6);
t605 = pkin(1) * t356;
t603 = pkin(3) * t350;
t599 = -qJD(1) / 0.2e1;
t598 = qJD(1) / 0.2e1;
t597 = pkin(2) - t335;
t593 = rSges(4,1) * t352;
t459 = t340 * t597;
t503 = t339 * t485 + t301;
t135 = (t459 + t600) * qJD(1) + t503;
t478 = t353 * t605;
t427 = -t356 * t252 - t478;
t430 = qJD(1) * t475;
t479 = qJD(3) ^ 2 * t345;
t371 = qJD(1) * t135 + t339 * t430 - t340 * t479 + t427;
t508 = -t279 + t335;
t584 = -qJD(1) * t394 + t241 * t339 - t487 * t508 + t503 - t635;
t26 = -t221 * t245 + t235 * t265 + (t235 * t341 - t245 * t539) * pkin(4) + (t328 - t584) * qJD(1) + t371;
t588 = t26 * t340;
t586 = -rSges(6,3) + t346;
t557 = t186 * t350;
t556 = t188 * t352;
t225 = t306 * t339;
t551 = t306 * t340;
t231 = t312 * t339;
t550 = t312 * t340;
t238 = pkin(6) * t487 - t326;
t536 = -t134 - t238;
t530 = t535 * t340;
t529 = t340 * t184 + t186 * t546;
t528 = t339 * t184 + t188 * t541;
t505 = t296 + t549;
t489 = qJD(1) * t307;
t402 = t308 * t350 - t352 * t640;
t130 = -t340 * t402 + t225;
t480 = t130 * qJD(1);
t324 = pkin(3) * t546;
t477 = rSges(5,1) * t547;
t470 = rSges(5,1) * t473 + rSges(5,2) * t662;
t336 = t351 * t605;
t469 = t339 * t479 + t340 * t430 + t336;
t468 = rSges(4,1) * t462 + (t339 * t481 + t350 * t487) * rSges(4,2);
t460 = -pkin(2) - t593;
t458 = -t488 / 0.2e1;
t457 = t487 / 0.2e1;
t456 = -t484 / 0.2e1;
t454 = -t483 / 0.2e1;
t449 = -t335 - t592;
t448 = t279 + t591;
t438 = -t185 - t556;
t428 = -pkin(2) - t449;
t426 = -pkin(4) * t539 - t221;
t421 = -t335 + t448;
t295 = -t602 - t603;
t416 = -rSges(4,2) * t350 + t593;
t407 = -t556 + t557;
t398 = t477 - t506;
t80 = t340 * t185 + t187 * t546 - t189 * t545;
t397 = -t506 + t652;
t79 = -t188 * t545 + t529;
t391 = (t339 * t80 + t340 * t79) * qJD(3);
t81 = -t186 * t542 + t528;
t82 = t185 * t339 - t649;
t390 = (t339 * t82 + t340 * t81) * qJD(3);
t388 = t339 * t597 - t652;
t378 = qJD(1) * t406 - qJD(3) * t551 + t492;
t377 = qJD(3) * t225 + (-t307 * t340 + t407 - t565) * qJD(1);
t374 = t402 * qJD(1) + t307 * qJD(3);
t373 = -qJD(1) * t388 + t340 * t475 + t422;
t372 = t340 * t192 + t339 * t399;
t370 = (-pkin(2) + t448) * t339 + (pkin(6) + t346) * t340 - t507;
t365 = t339 * t421 + t393;
t364 = -t161 * t483 - t388 * t484 + qJD(2);
t131 = -qJD(3) * t550 + (-t416 * t339 + t332) * qJD(1);
t132 = qJD(1) * t396 + t468;
t359 = t340 * t131 - t339 * t132 + (-t339 * t192 + t340 * t399) * qJD(1);
t358 = t134 * t483 + (-t339 * t135 + (t338 * t596 + (t161 - t459) * t339) * qJD(1)) * qJD(3);
t357 = (t339 * t678 + t340 * t679) * t616 + (t339 * t676 + t340 * t677) * t615 + (t664 * t339 + t694 * t340) * t614 + (t689 * t340 + t688 * t339 + (-t339 * t677 + t340 * t676) * qJD(1)) * t613 + (-t694 * t339 + t664 * t340) * t612 + (t687 * t340 + t686 * t339 + (-t339 * t679 + t340 * t678) * qJD(1)) * t611 + (t681 * qJD(1) - t677 * t235 + t676 * t236 + t688 * t244 + t245 * t689) * t610 + (qJD(1) * t680 - t235 * t679 + t236 * t678 + t244 * t686 + t245 * t687) * t609 + (t665 * t341 + t701 * t342) * t599 + (t683 * t340 + t682 * t339 + (t339 * t675 + t340 * t674) * qJD(1)) * t598 + t685 * t458 + t684 * t457;
t325 = rSges(3,2) * t488;
t277 = t416 * qJD(3);
t212 = qJD(1) * t417 + t476;
t191 = (t295 + t603) * t340;
t190 = -t295 * t339 - t324;
t155 = t340 * t180;
t145 = t340 * t161;
t129 = t339 * t402 + t551;
t121 = t129 * qJD(1);
t116 = -qJD(1) * t395 + t470;
t94 = qJD(3) * t372 + qJD(2);
t78 = -t277 * t483 + (t132 + t239) * qJD(1) + t427;
t77 = t277 * t484 + t336 + (-t131 - t238 + t463) * qJD(1);
t67 = qJD(1) * t398 + t245 * t266 + t373;
t66 = t339 * t629 + t374 * t340;
t65 = t374 * t339 - t340 * t629;
t64 = t245 * t180 + t244 * t398 + t364;
t63 = -t406 * qJD(3) + t125 * t352 + t127 * t350;
t62 = -qJD(3) * t407 + t126 * t352 + t128 * t350;
t61 = t359 * qJD(3);
t60 = (t425 + t535) * qJD(1) + t632;
t59 = t450 * t245 - t327 + t373 + (rSges(6,1) * t547 - t339 * t508 + t393) * qJD(1);
t49 = qJD(1) * t116 - t222 * t245 + t235 * t266 + t371;
t48 = t222 * t244 + t236 * t266 + (-t114 + t536) * qJD(1) + t469;
t43 = t390 + t480;
t42 = t121 + t391;
t41 = t244 * t365 - t245 * t535 + t364;
t25 = t221 * t244 + t236 * t265 + (t236 * t341 + t244 * t539) * pkin(4) + (-t327 + t536 - t585) * qJD(1) + t469;
t16 = t245 * t114 - t244 * t116 - t235 * t180 + t236 * t398 + t358;
t9 = t235 * t535 + t236 * t365 + t244 * t584 + t245 * t585 + t358;
t1 = [m(3) * ((t356 * t417 + t336) * t643 - t212 * t325 + t478 * t663) + (-qJD(3) * t402 + t271 * t352 + t272 * t350 + m(3) * (-t212 * t420 + (rSges(3,1) * t487 + qJD(1) * t643 - t325) * t663) + t716 * t342 - t717 * t341) * qJD(1) + ((t666 + t676) * t245 + t690) * t614 + (t43 - t480 + ((t80 + (-t185 + t557) * t340 - t528) * t340 + (t339 * t438 + t529 - t79) * t339) * qJD(3)) * t454 + (t25 * (-t606 + t692) - t60 * t327 + t26 * (t507 - t607) - t59 * t635 + (-t25 * t591 + t60 * (rSges(6,1) * t540 + rSges(6,2) * t539 - t241) - t26 * t346) * t340 + (t59 * t241 + t25 * t586 - t26 * t448) * t339 + ((t351 * t60 + t353 * t59) * pkin(1) + (t448 * t59 + t586 * t60) * t340 + (t60 * (t267 + t279) + t59 * rSges(6,3)) * t339) * qJD(1) - (t60 + (-t535 + t606) * qJD(1) - t632 + t645) * t59) * m(6) + (t48 * (t294 + t644) + t49 * t661 + (-t49 * t354 - t48 * t592) * t340 + (t48 * (-rSges(5,3) + t354) + t49 * t449) * t339 + (t509 + (rSges(5,1) * t540 + rSges(5,2) * t539 + t475) * t340 + (t477 - t661) * qJD(1)) * t68 + (-t470 - t503 - t68 + t432 - t645 + (t395 - t644 - t180 - t606) * qJD(1)) * t67) * m(5) + (t77 * (t323 - t606) + t100 * t326 + t78 * (t498 - t607) - t99 * t468 + (t100 * t312 * qJD(3) + t78 * pkin(6) + t77 * t460) * t340 + (t78 * t460 + t77 * t608) * t339 + ((t100 * t351 + t353 * t99) * pkin(1) + (t100 * t608 - t460 * t99) * t340 + (t100 * t416 - t608 * t99) * t339) * qJD(1) - (t100 - t239 + t243 + (t192 + t606) * qJD(1)) * t99) * m(4) + (t63 + t65 + t42) * t484 / 0.2e1 + (-t675 + t708) * t616 + (t674 - t707) * t615 + (((-t724 + t709) * t340 - t718 * t339 + t678 + t705) * t245 + (t666 - t679) * t244 + t684 + t691) * t612 + (t680 + t683) * t611 + (t681 + t682 + t685) * t613 + (t121 + ((t529 + t82 + t649) * t340 + (-t81 + (t438 - t557) * t340 + t80 + t528) * t339) * qJD(3) + (t117 + t129) * qJD(1)) * t456 + (t62 + t66 + (t118 + t130) * qJD(1)) * t483 / 0.2e1; m(4) * t61 + m(5) * t16 + m(6) * t9; ((t225 * t483 + t489) * t340 + (t657 + (t626 * t339 + (-t551 - t625) * t340) * qJD(3)) * t339) * t454 + t357 + ((-t350 * t500 + t352 * t499) * qJD(1) + ((t339 * t515 + t340 * t516) * t352 + (-t339 * t517 - t340 * t518) * t350) * qJD(3)) * t599 + ((-t484 * t551 + t489) * t339 + (-t657 + (t625 * t340 + (t225 - t626) * t339) * qJD(3)) * t340) * t456 + (t339 * t63 + t340 * t62 + (-t117 * t339 + t118 * t340) * qJD(1)) * t598 + (t65 * qJD(1) + ((t377 * t339 - t340 * t631) * t340 + (t378 * t339 - t340 * t630) * t339 + (-t81 * t339 + t82 * t340) * qJD(1)) * t656) * t610 + (t66 * qJD(1) + ((t339 * t631 + t377 * t340) * t340 + (t339 * t630 + t378 * t340) * t339 + (-t79 * t339 + t80 * t340) * qJD(1)) * t656) * t609 + (t42 + t391) * t458 + (t390 + t43) * t457 + (t9 * (t339 * t370 - t145 - t530) + t25 * (t324 + t505) + (-t265 + t295) * t588 + (-(-t191 + t210) * qJD(1) + t637) * t60 + (-(t426 - t474) * t340 + qJD(1) * t190 - t638 + t639) * t59 + ((-t135 + t584) * t339 + ((t161 + t535) * t339 + t370 * t340) * qJD(1) - t191 * t245 - (-t190 - t549) * t244 + t633 + t646) * t41) * m(6) + (t16 * (-t145 + t155 + (t339 * t428 + t397) * t339) + t48 * (t324 + t209) + t49 * (-t266 - t603) * t340 + t628 + (-(-t222 - t474) * t340 + t636 - t638) * t67 + ((-t135 - t116) * t339 + (t397 * t340 + (t340 * t428 + t161 - t180) * t339) * qJD(1) + t646 + t653) * t64) * m(5) + (-(t100 * t550 - t231 * t99) * qJD(1) - (t94 * (-t231 * t339 - t340 * t550) + t634 * t416) * qJD(3) + t77 * t231 + t359 * t94 + t372 * t61 - t78 * t550 + (t100 * t487 - t488 * t99) * t312 + t634 * t277) * m(4); t357 + (t9 * (t339 * t365 - t530) + t25 * t505 - t450 * t588 + (-qJD(1) * t210 - t276 + t637) * t60 + (-t340 * t426 + t275 + t639) * t59 + (t584 * t339 + (t393 * t340 + (t340 * t421 + t535) * t339) * qJD(1) + t549 * t244 - (-t244 * t339 - t245 * t340) * t602 + t633) * t41) * m(6) + (t16 * (t339 * t398 + t155) + t48 * t209 - t49 * t552 + (-t339 * t116 + (-t339 * t180 + t340 * t398) * qJD(1) + t653) * t64 + t628 + (t222 * t340 + t636) * t67) * m(5); 0.2e1 * (t25 * t609 + t26 * t610 - t41 * (-t244 * t340 + t245 * t339) / 0.2e1) * m(6);];
tauc = t1(:);
