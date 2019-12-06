% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP4
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:33
% EndTime: 2019-12-05 18:06:22
% DurationCPUTime: 26.39s
% Computational Cost: add. (56417->598), mult. (76013->828), div. (0->0), fcn. (82179->8), ass. (0->370)
t856 = -Icges(5,4) - Icges(6,4);
t508 = cos(pkin(8));
t510 = sin(qJ(1));
t661 = qJ(3) + qJ(4);
t566 = cos(t661);
t559 = t510 * t566;
t503 = sin(t661);
t511 = cos(qJ(1));
t663 = t511 * t503;
t459 = t508 * t663 - t559;
t434 = t459 * pkin(4);
t558 = t511 * t566;
t665 = t510 * t503;
t460 = t508 * t558 + t665;
t552 = t459 * rSges(6,1) + t460 * rSges(6,2);
t295 = t552 + t434;
t457 = t508 * t665 + t558;
t544 = t508 * t559;
t458 = t544 - t663;
t365 = t457 * rSges(5,1) + t458 * rSges(5,2);
t367 = -t459 * rSges(5,1) - t460 * rSges(5,2);
t364 = t457 * rSges(6,1) + t458 * rSges(6,2);
t620 = -pkin(4) * t457 - t364;
t744 = m(6) / 0.2e1;
t745 = m(5) / 0.2e1;
t653 = (t295 * t511 - t510 * t620) * t744 + (t510 * t365 - t367 * t511) * t745;
t507 = sin(pkin(8));
t560 = t507 * t566;
t787 = (pkin(4) * t560 + (rSges(6,1) * t566 - rSges(6,2) * t503) * t507 + (-rSges(6,3) - qJ(5)) * t508) * t507;
t626 = t787 * t510;
t509 = sin(qJ(3));
t662 = t511 * t509;
t499 = pkin(3) * t662;
t743 = -pkin(7) - pkin(6);
t589 = t507 * t743;
t602 = t510 * t589 + t499;
t673 = t507 * t510;
t693 = pkin(4) * t503;
t694 = pkin(3) * t509;
t561 = t693 + t694;
t600 = -qJ(5) + t743;
t564 = t507 * t600;
t763 = t458 * rSges(6,1) - t457 * rSges(6,2) - t510 * t564 - t511 * t561;
t644 = -rSges(6,3) * t673 - pkin(4) * t544 - t602 - t763;
t159 = t508 * t644 - t626;
t625 = t787 * t511;
t494 = t511 * t589;
t719 = cos(qJ(3));
t596 = t719 * pkin(3);
t545 = t508 * (t596 + pkin(2));
t664 = t510 * t509;
t598 = pkin(3) * t664;
t672 = t507 * t511;
t487 = pkin(4) * t566 + t596;
t522 = t508 * (pkin(2) + t487);
t839 = -t460 * rSges(6,1) + t459 * rSges(6,2) - t510 * t561 - t511 * t522;
t642 = t598 - t494 + (t545 + t564) * t511 - rSges(6,3) * t672 + t839;
t161 = -t642 * t508 + t625;
t611 = t458 * rSges(5,1) - t457 * rSges(5,2);
t330 = -rSges(5,3) * t673 - t611;
t516 = t507 * (-t508 * rSges(5,3) + (rSges(5,1) * t566 - rSges(5,2) * t503) * t507);
t395 = t510 * t516;
t261 = t330 * t508 - t395;
t765 = -t460 * rSges(5,1) + t459 * rSges(5,2);
t332 = rSges(5,3) * t672 - t765;
t309 = t508 * t332;
t397 = t511 * t516;
t262 = t309 + t397;
t660 = (-t159 * t511 - t161 * t510) * t744 + (-t261 * t511 - t262 * t510) * t745;
t18 = t660 - t653;
t855 = qJD(1) * t18;
t847 = Icges(5,1) + Icges(6,1);
t846 = Icges(5,5) + Icges(6,5);
t845 = Icges(5,2) + Icges(6,2);
t844 = Icges(5,6) + Icges(6,6);
t854 = t856 * t566;
t853 = t856 * t503;
t579 = t511 * t719;
t481 = t508 * t579 + t664;
t369 = -pkin(3) * t481 + pkin(6) * t672 + t494;
t340 = t508 * t369;
t530 = t507 * (-t508 * pkin(7) + t507 * t596);
t393 = t511 * t530;
t132 = -t340 + t393 + t161;
t843 = t369 + t642;
t659 = -t843 * t508 - t132 + t393 + t625;
t852 = m(6) * t659;
t851 = t659 * t744;
t618 = t393 + t397;
t200 = -t340 + t309 + t618;
t828 = t332 - t369;
t656 = t508 * t828 - t200 + t618;
t850 = m(5) * t656;
t849 = t656 * t745;
t841 = t844 * t508 + (t503 * t845 + t854) * t507;
t848 = -t846 * t508 + (t566 * t847 + t853) * t507;
t788 = (-t503 * t847 + t854) * t507;
t312 = Icges(6,5) * t460 - Icges(6,6) * t459 + Icges(6,3) * t672;
t315 = Icges(5,5) * t460 - Icges(5,6) * t459 + Icges(5,3) * t672;
t829 = -t312 - t315;
t824 = t848 + (-t566 * t845 + t853) * t507;
t823 = -t788 - t841;
t687 = Icges(4,4) * t509;
t450 = -Icges(4,5) * t508 + (Icges(4,1) * t719 - t687) * t507;
t476 = (-Icges(4,2) * t719 - t687) * t507;
t578 = t719 * Icges(4,4);
t477 = (-Icges(4,1) * t509 - t578) * t507;
t576 = t719 * t477;
t449 = -Icges(4,6) * t508 + (-Icges(4,2) * t509 + t578) * t507;
t577 = t719 * t449;
t475 = (-Icges(4,5) * t509 - Icges(4,6) * t719) * t507;
t667 = t508 * t475;
t842 = t667 / 0.2e1 - (-(t450 / 0.2e1 + t476 / 0.2e1) * t509 - t577 / 0.2e1 + t576 / 0.2e1) * t507;
t506 = t507 ^ 2;
t774 = m(6) * t507;
t570 = t672 / 0.2e1;
t571 = -t672 / 0.2e1;
t573 = -t673 / 0.2e1;
t314 = -Icges(5,5) * t458 + Icges(5,6) * t457 - Icges(5,3) * t673;
t676 = t314 * t510;
t311 = -Icges(6,5) * t458 + Icges(6,6) * t457 - Icges(6,3) * t673;
t677 = t311 * t510;
t686 = Icges(5,4) * t458;
t320 = Icges(5,2) * t457 - Icges(5,6) * t673 - t686;
t430 = Icges(5,4) * t457;
t326 = -Icges(5,1) * t458 - Icges(5,5) * t673 + t430;
t645 = -t457 * t320 + t458 * t326;
t830 = t315 * t672 - t645;
t684 = Icges(6,4) * t458;
t317 = Icges(6,2) * t457 - Icges(6,6) * t673 - t684;
t427 = Icges(6,4) * t457;
t323 = -Icges(6,1) * t458 - Icges(6,5) * t673 + t427;
t646 = -t457 * t317 + t458 * t323;
t831 = t312 * t672 - t646;
t432 = Icges(5,4) * t460;
t321 = -Icges(5,2) * t459 + Icges(5,6) * t672 + t432;
t431 = Icges(5,4) * t459;
t328 = -Icges(5,1) * t460 - Icges(5,5) * t672 + t431;
t647 = t457 * t321 + t328 * t458;
t429 = Icges(6,4) * t460;
t318 = -Icges(6,2) * t459 + Icges(6,6) * t672 + t429;
t428 = Icges(6,4) * t459;
t325 = -Icges(6,1) * t460 - Icges(6,5) * t672 + t428;
t648 = t457 * t318 + t325 * t458;
t838 = t647 + t648;
t832 = t829 * t673 + t838;
t833 = t645 + t646 + (t311 + t314) * t673;
t754 = t510 ^ 2;
t834 = t754 * t507;
t778 = (t829 * t834 + (-t832 + t838) * t510 + ((t511 * t829 - t676 - t677) * t507 + t830 + t831 + t833) * t511) * t507;
t451 = (-Icges(6,5) * t503 - Icges(6,6) * t566) * t507;
t452 = (-Icges(5,5) * t503 - Icges(5,6) * t566) * t507;
t825 = t451 + t452;
t795 = t459 * t824 + t460 * t823 - t672 * t825;
t796 = -t457 * t824 - t458 * t823 + t673 * t825;
t783 = t847 * t459 + t318 + t321 + t429 + t432;
t785 = -t845 * t460 - t325 - t328 - t428 - t431;
t827 = -t846 * t459 - t844 * t460;
t797 = -t827 * t508 + (-t503 * t785 - t566 * t783) * t507;
t782 = -t847 * t457 + t317 + t320 - t684 - t686;
t784 = t845 * t458 + t323 + t326 + t427 + t430;
t826 = t846 * t457 + t844 * t458;
t798 = t826 * t508 + (t503 * t784 + t566 * t782) * t507;
t816 = (t848 * t458 + t841 * t457 + ((-Icges(5,3) - Icges(6,3)) * t508 + (-t844 * t503 + t846 * t566) * t507) * t673) * t508;
t806 = t816 + (t510 * t833 + t511 * t832) * t507;
t807 = (t648 * t511 + (t677 * t507 - t831) * t510) * t507 + (t647 * t511 + (t676 * t507 - t830) * t510) * t507 + t816;
t514 = t807 * t571 + (-t796 - t798) * t573 + t778 * t673 / 0.2e1 + (-t795 + t797 + t806) * t570;
t668 = t508 * t452;
t669 = t508 * t451;
t654 = t668 + t669 + (t824 * t503 + t823 * t566) * t507;
t780 = t654 * t508;
t840 = t514 + t780;
t504 = t511 * qJ(2);
t247 = -t504 + (t507 * rSges(6,3) + pkin(1) + t522) * t510 + t763;
t666 = t510 * qJ(2);
t248 = (-pkin(1) + (-rSges(6,3) + t600) * t507) * t511 - t666 + t839;
t836 = (t247 * t510 - t248 * t511) * t774;
t580 = t510 * t719;
t588 = t508 * t662;
t480 = -t580 + t588;
t356 = Icges(4,5) * t481 - Icges(4,6) * t480 + Icges(4,3) * t672;
t478 = t508 * t664 + t579;
t479 = t508 * t580 - t662;
t688 = Icges(4,4) * t479;
t358 = Icges(4,2) * t478 - Icges(4,6) * t673 - t688;
t466 = Icges(4,4) * t478;
t361 = -Icges(4,1) * t479 - Icges(4,5) * t673 + t466;
t639 = -t478 * t358 + t479 * t361;
t822 = t356 * t672 - t639;
t764 = -t481 * rSges(4,1) + t480 * rSges(4,2);
t371 = rSges(4,3) * t672 - t764;
t519 = t507 * (-t508 * rSges(4,3) + (rSges(4,1) * t719 - rSges(4,2) * t509) * t507);
t271 = t371 * t508 + t511 * t519;
t821 = qJD(1) * t836;
t814 = t295 * t510 + t620 * t511;
t468 = Icges(4,4) * t481;
t359 = -Icges(4,2) * t480 + Icges(4,6) * t672 + t468;
t467 = Icges(4,4) * t480;
t363 = -Icges(4,1) * t481 - Icges(4,5) * t672 + t467;
t640 = t478 * t359 + t363 * t479;
t810 = -t672 * t827 + t673 * t826;
t746 = m(4) / 0.2e1;
t526 = t507 * rSges(5,3) + pkin(1) + t545;
t260 = t494 - t526 * t511 + (-qJ(2) - t694) * t510 + t765;
t756 = t508 * pkin(2) + pkin(1) + (rSges(4,3) + pkin(6)) * t507;
t292 = -t511 * t756 - t666 + t764;
t392 = t510 * t530;
t368 = (t507 * pkin(6) - t508 * t596) * t510 + t602;
t583 = t368 + t644;
t130 = t508 * t583 - t392 - t626;
t630 = t330 + t368;
t198 = t508 * t630 - t392 - t395;
t605 = t479 * rSges(4,1) - t478 * rSges(4,2);
t370 = -rSges(4,3) * t673 - t605;
t270 = t370 * t508 - t510 * t519;
t590 = (-t198 * t511 - t200 * t510) * t745 + (-t270 * t511 - t271 * t510) * t746 + (-t130 * t511 - t132 * t510) * t744;
t597 = (t851 + t850 / 0.2e1) * t510;
t786 = t590 - t597;
t461 = (-rSges(6,1) * t503 - rSges(6,2) * t566) * t507;
t412 = t461 * t673;
t674 = t503 * t507;
t426 = pkin(4) * t674;
t617 = -t426 * t673 + t412;
t540 = t508 * t561;
t267 = t511 * t487 + t510 * t540 + t364;
t469 = t478 * pkin(3);
t629 = t469 - t267;
t212 = t508 * t629 + t617;
t490 = pkin(3) * t588;
t532 = t511 * t540;
t337 = t490 - t532 + (-t596 + t487) * t510;
t627 = t461 * t672 - t508 * t552;
t213 = t508 * t337 - t426 * t672 + t627;
t268 = -t510 * t487 + t532 + t552;
t470 = pkin(3) * t580 - t490;
t299 = -t367 - t470;
t619 = -t365 - t469;
t259 = t510 * t526 - t504 - t602 + t611;
t462 = (-rSges(5,1) * t503 - rSges(5,2) * t566) * t507;
t413 = t462 * t673;
t273 = -t365 * t508 + t413;
t274 = t508 * t367 + t462 * t672;
t655 = -t273 * t259 + t274 * t260;
t691 = (-t159 * t267 + t161 * t268 - t212 * t247 + t213 * t248) * t744 + (t261 * t619 + t262 * t299 + t655) * t745;
t599 = t506 * t693;
t241 = t508 * t620 - t510 * t599 + t412;
t242 = -t508 * t434 - t511 * t599 + t627;
t658 = -t241 * t247 + t242 * t248;
t692 = (t130 * t620 + t132 * t295 + t658) * t744 + (-t198 * t365 - t200 * t367 + t655) * t745;
t779 = t691 - t692;
t777 = -t668 / 0.2e1 - t669 / 0.2e1 - t824 * t674 / 0.2e1;
t776 = -t507 / 0.2e1;
t721 = -t508 / 0.2e1;
t624 = -Icges(4,1) * t478 + t358 - t688;
t623 = Icges(4,1) * t480 + t359 + t468;
t757 = -t159 * t851 - t261 * t849;
t381 = rSges(4,1) * t478 + rSges(4,2) * t479;
t382 = -rSges(4,1) * t480 - rSges(4,2) * t481;
t584 = (t510 * t267 + t268 * t511) * t744 + (t299 * t511 - t510 * t619) * t745 + (t510 * t381 - t382 * t511) * t746;
t751 = 0.4e1 * qJD(1);
t750 = 2 * qJD(3);
t749 = 4 * qJD(3);
t748 = 2 * qJD(4);
t740 = t198 * t850;
t152 = (-t510 * t828 - t630 * t511) * t507;
t252 = (-t365 * t511 - t367 * t510) * t507;
t63 = t152 * t252 - t198 * t273 + t200 * t274;
t739 = m(5) * t63;
t239 = (-t330 * t511 - t332 * t510) * t507;
t94 = t239 * t252 - t261 * t273 + t262 * t274;
t93 = m(5) * t94;
t733 = t130 * t852;
t732 = t672 * t852;
t107 = (t843 * t510 - t583 * t511) * t507;
t214 = t814 * t507;
t40 = t107 * t214 - t130 * t241 + t132 * t242;
t730 = m(6) * t40;
t141 = (t510 * t642 - t511 * t644) * t507;
t628 = -t337 + t552;
t157 = (t510 * t628 + t511 * t629) * t507;
t729 = m(6) * (t141 * t157 - t159 * t212 + t161 * t213);
t723 = (t130 * t510 - t132 * t511) * t774;
t718 = m(3) * (-((-rSges(3,3) - qJ(2)) * t510 + rSges(3,2) * t672) * t510 - t511 * (-rSges(3,2) * t673 - t511 * rSges(3,3) - t504));
t291 = t510 * t756 - t504 + t605;
t716 = m(4) * (-t291 * t381 - t292 * t382);
t714 = m(4) * (-t511 * t291 - t292 * t510);
t710 = m(5) * (t259 * t619 + t260 * t299);
t709 = m(5) * (-t259 * t365 - t260 * t367);
t708 = m(5) * (-t511 * t259 - t260 * t510);
t704 = (t159 * t510 - t161 * t511) * t774;
t702 = m(6) * (-t247 * t267 + t248 * t268);
t701 = m(6) * (t247 * t620 + t248 * t295);
t700 = m(6) * (-t511 * t247 - t248 * t510);
t699 = (t267 * t511 - t268 * t510) * t774;
t697 = t814 * t774;
t355 = -Icges(4,5) * t479 + Icges(4,6) * t478 - Icges(4,3) * t673;
t675 = t355 * t510;
t622 = Icges(4,2) * t479 + t361 + t466;
t621 = -Icges(4,2) * t481 - t363 - t467;
t608 = t449 - t477;
t607 = t450 + t476;
t373 = (t511 ^ 2 + t754) * t774;
t601 = t373 * qJD(1);
t192 = -t355 * t673 - t639;
t193 = -t356 * t673 + t640;
t592 = (-t356 * t834 + (-t193 + t640) * t510 + (-t192 + (-t356 * t511 - t675) * t507 + t822) * t511) * t776 + (t721 + t508 / 0.2e1) * ((-Icges(4,3) * t508 + (Icges(4,5) * t719 - Icges(4,6) * t509) * t507) * t672 - t480 * t449 + t481 * t450);
t591 = (-t192 * t510 + t193 * t511) * t776 + (t640 * t511 + (t675 * t507 - t822) * t510) * t507 / 0.2e1;
t585 = t141 * t214 - t159 * t241 + t161 * t242;
t582 = -t469 + t629;
t565 = t506 * t598;
t551 = ((t798 * t510 + t797 * t511) * t507 + t780) * t721 + ((t457 * t785 + t458 * t783) * t672 + t796 * t508 + (-t784 * t457 - t782 * t458 + t810) * t673) * t573 + ((t459 * t784 + t460 * t782) * t673 + t795 * t508 + (-t785 * t459 - t783 * t460 - t810) * t672) * t570;
t518 = (-t508 * t157 + (t212 * t511 - t213 * t510) * t507) * t744;
t521 = m(6) * (-t508 * t214 + (t241 * t511 - t242 * t510) * t507);
t54 = t518 - t521 / 0.2e1;
t531 = (t212 * t510 + t213 * t511) * t744;
t533 = m(6) * (t241 * t510 + t242 * t511);
t79 = t531 - t533 / 0.2e1;
t550 = -t79 * qJD(2) - t54 * qJD(5);
t546 = t508 * t470 - t499 * t506;
t543 = -t560 / 0.2e1;
t542 = t560 / 0.2e1;
t541 = t93 + t551;
t527 = t807 * t570 + t806 * t571 + t778 * t573;
t525 = t788 * t542 - t543 * t841 + t777;
t523 = t527 + t757;
t512 = -t542 * t841 + t788 * t543 - t777;
t482 = (-rSges(4,1) * t509 - rSges(4,2) * t719) * t507;
t376 = -Icges(4,5) * t480 - Icges(4,6) * t481;
t375 = Icges(4,5) * t478 + Icges(4,6) * t479;
t297 = t508 * t382 + t482 * t672;
t296 = -t381 * t508 + t482 * t673;
t245 = t546 + t274;
t244 = t508 * t619 + t413 - t565;
t238 = -t667 + (-t509 * t607 + t576 - t577) * t507;
t219 = -t697 / 0.2e1;
t216 = (t299 * t510 + t511 * t619) * t507;
t210 = m(5) * (t273 * t510 + t274 * t511);
t204 = t475 * t672 - t480 * t607 - t481 * t608;
t203 = -t475 * t673 + t478 * t607 + t479 * t608;
t196 = t699 / 0.2e1;
t180 = t546 + t213;
t179 = t508 * t582 - t565 + t617;
t149 = -t508 * t376 + (-t621 * t509 - t623 * t719) * t507;
t148 = -t508 * t375 + (-t622 * t509 - t624 * t719) * t507;
t144 = (t582 * t511 + (-t470 + t628) * t510) * t507;
t99 = t704 / 0.2e1;
t74 = t723 / 0.2e1;
t72 = t700 + t708 + t714 + t718;
t66 = t210 + t533 / 0.2e1 + t531;
t56 = t525 + t701 + t709;
t55 = t521 / 0.2e1 + t518;
t41 = t716 + t710 + t702 + t525 - t842;
t25 = t732 / 0.2e1;
t22 = t99 + t219;
t21 = t99 + t697 / 0.2e1;
t20 = t219 - t704 / 0.2e1;
t16 = t653 + t660;
t15 = t74 + t196 - t732 / 0.2e1;
t14 = t74 + t25 - t699 / 0.2e1;
t13 = t196 + t25 - t723 / 0.2e1;
t10 = t590 + t597 - t584;
t9 = t584 - t786;
t8 = t584 + t786;
t7 = t541 + t729;
t6 = t551 + t730 + t739;
t4 = (t510 * t592 + t511 * t591) * t507 - t740 - t733 + t527;
t3 = t523 + t779;
t2 = t523 - t779;
t1 = t691 + t692 - t757 + t840;
t5 = [t72 * qJD(2) + t41 * qJD(3) + t56 * qJD(4) + qJD(5) * t836, qJD(1) * t72 + qJD(3) * t8 + qJD(4) * t16, t41 * qJD(1) + t8 * qJD(2) + t1 * qJD(4) + t15 * qJD(5) + (t740 / 0.4e1 + t733 / 0.4e1) * t749 + ((-t270 * t381 - t271 * t382 - t291 * t296 + t292 * t297) * t746 + (t198 * t619 + t200 * t299 - t244 * t259 + t245 * t260) * t745 + (-t130 * t267 + t132 * t268 - t179 * t247 + t180 * t248) * t744) * t750 + (t514 + (-t238 + t654) * t508 + ((t149 / 0.2e1 + t204 / 0.2e1 - t591) * t511 + (-t203 / 0.2e1 - t148 / 0.2e1 - t592) * t510) * t507) * qJD(3), t56 * qJD(1) + t16 * qJD(2) + t1 * qJD(3) + t840 * qJD(4) + t22 * qJD(5) + ((t159 * t620 + t161 * t295 + t658) * t744 + (-t261 * t365 - t262 * t367 + t655) * t745) * t748, t15 * qJD(3) + t22 * qJD(4) + t821; t9 * qJD(3) - t18 * qJD(4) - t373 * qJD(5) + (-t700 / 0.4e1 - t718 / 0.4e1 - t714 / 0.4e1 - t708 / 0.4e1) * t751, 0, t9 * qJD(1) + t66 * qJD(4) + ((t296 * t510 + t297 * t511) * t746 + (t244 * t510 + t245 * t511) * t745 + (t179 * t510 + t180 * t511) * t744) * t750, -t855 + t66 * qJD(3) + (t210 + t533) * qJD(4), -t601; t10 * qJD(2) + t4 * qJD(3) + t2 * qJD(4) + t14 * qJD(5) + (-t716 / 0.4e1 - t710 / 0.4e1 - t702 / 0.4e1) * t751 + 0.2e1 * (-t247 * t851 - t259 * t849) * qJD(1) + (t512 + t842) * qJD(1), qJD(1) * t10 - qJD(4) * t79, t4 * qJD(1) + (m(6) * (t107 * t144 - t130 * t179 + t132 * t180) + m(5) * (t152 * t216 - t198 * t244 + t200 * t245) + (-t238 * t508 + (-t148 * t510 + t149 * t511) * t507) * t721 + m(4) * (-t270 * t296 + t271 * t297 + (-t370 * t511 - t371 * t510) * t506 * (-t381 * t511 - t382 * t510)) + (-t203 * t508 - (-t375 * t673 + t478 * t622 + t479 * t624) * t673 + (-t376 * t673 + t478 * t621 + t479 * t623) * t672) * t573 + (-t204 * t508 - (t375 * t672 - t480 * t622 - t481 * t624) * t673 + (t376 * t672 - t480 * t621 - t481 * t623) * t672) * t570 + t551) * qJD(3) + t6 * qJD(4), t2 * qJD(1) + t6 * qJD(3) + ((t585 + t40) * t744 + (t94 + t63) * t745) * t748 + t550 + (t551 - t93 - t729) * qJD(4), qJD(1) * t14 - qJD(4) * t54; t512 * qJD(1) + t18 * qJD(2) + t3 * qJD(3) + t527 * qJD(4) + t21 * qJD(5) + (-t701 / 0.4e1 - t709 / 0.4e1) * t751, qJD(3) * t79 + t855, t3 * qJD(1) + t551 * qJD(3) + t7 * qJD(4) + (-t730 / 0.4e1 - t739 / 0.4e1) * t749 + ((t107 * t157 - t130 * t212 + t132 * t213 + t141 * t144 - t159 * t179 + t161 * t180) * t744 + (t216 * t239 - t244 * t261 + t245 * t262 + t63) * t745) * t750 - t550, t527 * qJD(1) + t7 * qJD(3) + (m(6) * t585 + t541) * qJD(4), qJD(1) * t21 + qJD(3) * t54; t373 * qJD(2) + t13 * qJD(3) + t20 * qJD(4) - t821, t601, t13 * qJD(1) + m(6) * (-t508 * t144 + (t179 * t511 - t180 * t510) * t507) * qJD(3) + t55 * qJD(4), t20 * qJD(1) + t55 * qJD(3) + qJD(4) * t521, 0;];
Cq = t5;
