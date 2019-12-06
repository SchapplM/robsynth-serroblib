% Calculate vector of inverse dynamics joint torques for
% S5RPRRP3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:52
% DurationCPUTime: 30.93s
% Computational Cost: add. (19387->738), mult. (16589->880), div. (0->0), fcn. (12988->8), ass. (0->397)
t782 = Icges(5,4) + Icges(6,4);
t759 = Icges(5,2) + Icges(6,2);
t783 = Icges(5,1) + Icges(6,1);
t781 = Icges(5,5) + Icges(6,5);
t780 = Icges(5,6) + Icges(6,6);
t378 = qJ(3) + qJ(4);
t371 = cos(t378);
t784 = t782 * t371;
t370 = sin(t378);
t773 = t759 * t370 - t784;
t779 = Icges(5,3) + Icges(6,3);
t377 = qJ(1) + pkin(8);
t369 = cos(t377);
t578 = t369 * t370;
t778 = t782 * t578;
t368 = sin(t377);
t777 = t780 * t368;
t776 = t781 * t368;
t775 = t782 * t370;
t771 = t783 * t370 + t784;
t757 = t773 * t368 + t780 * t369;
t577 = t369 * t371;
t761 = t782 * t577 - t759 * t578 + t777;
t755 = t783 * t577 + t776 - t778;
t753 = -t780 * t370 + t781 * t371;
t769 = t759 * t371 + t775;
t768 = t783 * t371 - t775;
t582 = t368 * t370;
t774 = t782 * t582;
t772 = t781 * t369;
t581 = t368 * t371;
t756 = -t783 * t581 + t772 + t774;
t770 = t781 * t370 + t780 * t371;
t767 = t779 * t368;
t766 = t771 - t773;
t765 = -t768 + t769;
t764 = t771 * t369 + t761;
t763 = -t771 * t368 + t757;
t762 = -t753 * t368 + t779 * t369;
t758 = t781 * t577 - t780 * t578 + t767;
t760 = t769 * t370 - t771 * t371;
t752 = -t761 * t370 + t755 * t371;
t754 = t770 * t368;
t376 = qJD(3) + qJD(4);
t751 = t765 * t376;
t750 = t766 * t376;
t749 = t763 * t376 + (t768 * t369 + t776) * qJD(1);
t748 = t764 * t376 + (t768 * t368 - t772) * qJD(1);
t265 = t368 * t376;
t747 = t756 * t376 + t769 * t265 + (t773 * t369 - t777) * qJD(1);
t266 = t369 * t376;
t746 = t757 * qJD(1) - t769 * t266 + t755 * t376;
t745 = t770 * t369;
t744 = t756 * t371;
t743 = t757 * t370;
t712 = t758 * t369 - t755 * t581 + t761 * t582;
t707 = t760 * t368 + t745;
t706 = -t760 * t369 + t754;
t742 = t762 * qJD(1);
t741 = -t762 * t368 - t756 * t577;
t740 = t762 * t369 + t757 * t582;
t739 = t752 * t369;
t713 = -t756 * t581 + t740;
t711 = -t757 * t578 - t741;
t710 = t758 * t368 + t739;
t737 = (t759 * t581 + t756 + t774) * t266 + (-t759 * t577 + t755 - t778) * t265 + t766 * qJD(1);
t736 = t760 * qJD(1) + t753 * t376;
t735 = t754 * t376 + (-t753 * t369 + t743 - t744 - t767) * qJD(1);
t734 = -t752 * qJD(1) - t745 * t376 + t742;
t733 = -t758 * qJD(1) + t746 * t370 + t748 * t371;
t732 = t747 * t370 + t749 * t371 - t742;
t731 = -t770 * qJD(1) + t750 * t370 + t751 * t371;
t697 = -t765 * qJD(1) - t764 * t265 - t763 * t266;
t730 = -t737 * t370 + t371 * t697;
t363 = pkin(4) * t371;
t381 = cos(qJ(3));
t374 = t381 * pkin(3);
t364 = t374 + pkin(2);
t509 = t363 + t364;
t728 = rSges(6,2) * t578 - t369 * t509;
t727 = t706 * qJD(1);
t726 = t707 * qJD(1) + t712 * t265;
t385 = qJD(1) ^ 2;
t352 = qJDD(3) * t369;
t521 = qJD(1) * t368;
t177 = qJDD(4) * t369 - t376 * t521 + t352;
t361 = t371 * rSges(6,1);
t623 = rSges(6,2) * t370;
t673 = t361 - t623;
t238 = t673 * t376;
t618 = t371 * rSges(6,2);
t286 = rSges(6,1) * t370 + t618;
t354 = qJD(5) * t369;
t629 = pkin(2) - t364;
t489 = t369 * t629;
t379 = sin(qJ(3));
t515 = qJD(3) * t379;
t494 = t368 * t515;
t324 = pkin(3) * t494;
t383 = -pkin(7) - pkin(6);
t518 = qJD(1) * t383;
t536 = t368 * t518 + t324;
t636 = pkin(6) * t368;
t142 = (t489 + t636) * qJD(1) + t536;
t511 = qJD(1) * qJD(3);
t259 = -t368 * t511 + t352;
t628 = pkin(6) + t383;
t679 = t628 * t369;
t412 = -t368 * t629 + t679;
t639 = pkin(2) * t369;
t273 = t636 + t639;
t635 = pkin(6) * t369;
t452 = -pkin(2) * t368 + t635;
t545 = qJDD(1) * t452 - t385 * t273;
t572 = t381 * qJD(3) ^ 2;
t395 = qJD(1) * t142 - qJDD(1) * t412 + (-t259 * t379 - t369 * t572) * pkin(3) + t545;
t382 = cos(qJ(1));
t571 = t382 * t385;
t510 = pkin(1) * t571;
t460 = t364 - t509;
t375 = -qJ(5) + t383;
t530 = t383 - t375;
t356 = t369 * rSges(6,3);
t539 = rSges(6,2) * t582 + t356;
t567 = rSges(6,1) * t581 - t368 * t460 - t369 * t530 - t539;
t573 = t371 * t376;
t503 = rSges(6,1) * t577;
t417 = -t368 * rSges(6,3) - t503;
t520 = qJD(1) * t369;
t574 = t370 * t376;
t502 = t368 * t574;
t496 = t370 * t520;
t501 = t368 * t573;
t695 = t496 + t501;
t506 = pkin(3) * t515;
t261 = -pkin(4) * t574 - t506;
t701 = qJD(1) * t375 - t261;
t668 = rSges(6,1) * t502 + rSges(6,2) * t695 + t368 * t701 + t354;
t616 = qJD(1) * t417 + t460 * t520 - t536 + t668;
t380 = sin(qJ(1));
t633 = t380 * pkin(1);
t12 = -t510 + qJDD(5) * t368 - t177 * t286 - t266 * t238 + (-t177 * t370 - t266 * t573) * pkin(4) + (-t567 - t633) * qJDD(1) + (t354 + t616) * qJD(1) + t395;
t725 = t12 - g(3);
t258 = qJDD(3) * t368 + t369 * t511;
t176 = qJD(4) * t520 + qJDD(4) * t368 + t258;
t353 = qJD(5) * t368;
t309 = t369 * t364;
t172 = t368 * t628 - t309 + t639;
t632 = t382 * pkin(1);
t479 = -t273 - t632;
t456 = t172 + t479;
t566 = -t368 * t530 + t309 + t417 + t728;
t416 = t456 + t566;
t365 = t385 * t633;
t638 = pkin(3) * t379;
t500 = t368 * pkin(3) * t572 + t258 * t638 + t365;
t350 = pkin(2) * t521;
t540 = t364 * t521 + t369 * t518;
t141 = t350 + (-pkin(6) * qJD(1) - t506) * t369 - t540;
t256 = pkin(6) * t520 - t350;
t568 = -t141 - t256;
t226 = t286 * t369;
t614 = -t509 * t521 + t353 + (t506 - t701) * t369 + t540 - t376 * t226 + (-t368 * t673 + t356) * qJD(1);
t13 = qJDD(5) * t369 + t176 * t286 + t238 * t265 + (t176 * t370 + t265 * t573) * pkin(4) + t416 * qJDD(1) + (-t353 + t568 - t614) * qJD(1) + t500;
t724 = t13 - g(2);
t723 = t735 * t368 - t732 * t369;
t722 = t734 * t368 - t733 * t369;
t721 = t732 * t368 + t735 * t369;
t720 = t733 * t368 + t734 * t369;
t719 = t713 * t266 + t726;
t718 = t710 * t265 + t711 * t266 + t727;
t717 = -t749 * t370 + t747 * t371;
t716 = -t748 * t370 + t746 * t371;
t715 = t736 * t368 - t731 * t369;
t714 = t731 * t368 + t736 * t369;
t709 = t756 * t370 + t757 * t371;
t708 = t755 * t370 + t761 * t371;
t698 = (-t758 - t744) * t368 - t739 + t740;
t696 = t753 * qJD(1) - t745 * t265 + t754 * t266;
t357 = t369 * rSges(5,3);
t538 = rSges(5,2) * t582 + t357;
t694 = t633 - t538;
t453 = t673 + t363;
t372 = Icges(4,4) * t381;
t437 = -Icges(4,2) * t379 + t372;
t671 = Icges(4,1) * t379 + t372;
t532 = t671 + t437;
t610 = Icges(4,4) * t379;
t331 = Icges(4,2) * t381 + t610;
t334 = Icges(4,1) * t381 - t610;
t533 = t331 - t334;
t690 = (t379 * t532 + t381 * t533) * qJD(1);
t689 = t12 * t369;
t287 = rSges(5,1) * t370 + rSges(5,2) * t371;
t227 = t287 * t369;
t362 = t371 * rSges(5,1);
t672 = -rSges(5,2) * t370 + t362;
t116 = -t376 * t227 + (-t368 * t672 + t357) * qJD(1);
t225 = rSges(5,1) * t582 + rSges(5,2) * t581;
t686 = t369 * t116 + t225 * t265 + t266 * t227;
t575 = t369 * t381;
t576 = t369 * t379;
t598 = Icges(4,6) * t368;
t202 = Icges(4,4) * t575 - Icges(4,2) * t576 + t598;
t344 = Icges(4,4) * t576;
t607 = Icges(4,5) * t368;
t204 = Icges(4,1) * t575 - t344 + t607;
t429 = t202 * t379 - t204 * t381;
t684 = t369 * t429;
t678 = t369 * t141 - (-t368 ^ 2 - t369 ^ 2) * t506;
t264 = qJD(1) * t273;
t677 = -qJD(1) * t172 + t264;
t676 = -t309 - t632;
t626 = rSges(3,1) * t369;
t450 = -t626 - t632;
t253 = t368 * rSges(3,2) + t450;
t224 = rSges(6,1) * t582 + rSges(6,2) * t581;
t319 = pkin(4) * t582;
t296 = qJD(1) * t319;
t670 = qJD(1) * t224 - t266 * t453 - t286 * t521 - t296;
t514 = qJD(3) * t381;
t491 = t369 * t514;
t669 = pkin(3) * t491;
t297 = pkin(4) * t496;
t667 = pkin(4) * t501 + t368 * t238 - t265 * t453 + t286 * t520 + t297;
t666 = t266 * t226 + t369 * t614;
t637 = pkin(4) * t370;
t477 = t286 + t637;
t665 = -t477 * t265 - t324 - t354;
t201 = Icges(4,6) * t369 - t368 * t437;
t580 = t368 * t379;
t343 = Icges(4,4) * t580;
t579 = t368 * t381;
t606 = Icges(4,5) * t369;
t203 = -Icges(4,1) * t579 + t343 + t606;
t121 = t201 * t381 + t203 * t379;
t517 = qJD(3) * t368;
t132 = t331 * t517 + (-t369 * t437 - t598) * qJD(1);
t246 = t671 * t368;
t134 = qJD(3) * t246 + (-t334 * t369 - t607) * qJD(1);
t330 = Icges(4,5) * t381 - Icges(4,6) * t379;
t199 = Icges(4,3) * t369 - t330 * t368;
t525 = qJD(1) * t199;
t664 = qJD(3) * t121 + t132 * t379 - t134 * t381 - t525;
t122 = t202 * t381 + t204 * t379;
t516 = qJD(3) * t369;
t131 = qJD(1) * t201 - t331 * t516;
t247 = t671 * t369;
t133 = -qJD(3) * t247 + (-t334 * t368 + t606) * qJD(1);
t595 = Icges(4,3) * t368;
t200 = Icges(4,5) * t575 - Icges(4,6) * t576 + t595;
t663 = -qJD(1) * t200 + qJD(3) * t122 + t131 * t379 - t133 * t381;
t292 = t437 * qJD(3);
t293 = t334 * qJD(3);
t329 = Icges(4,5) * t379 + Icges(4,6) * t381;
t425 = t331 * t381 + t379 * t671;
t662 = -qJD(1) * t329 + qJD(3) * t425 + t292 * t379 - t293 * t381;
t547 = -Icges(4,2) * t575 + t204 - t344;
t549 = t202 + t247;
t660 = t379 * t547 + t381 * t549;
t548 = Icges(4,2) * t579 + t203 + t343;
t550 = t201 - t246;
t659 = -t379 * t548 - t381 * t550;
t650 = t176 / 0.2e1;
t649 = t177 / 0.2e1;
t648 = t258 / 0.2e1;
t647 = t259 / 0.2e1;
t646 = -t265 / 0.2e1;
t645 = t265 / 0.2e1;
t644 = -t266 / 0.2e1;
t643 = t266 / 0.2e1;
t642 = t368 / 0.2e1;
t641 = t369 / 0.2e1;
t640 = -rSges(4,3) - pkin(6);
t634 = g(3) * t369;
t631 = -qJD(1) / 0.2e1;
t630 = qJD(1) / 0.2e1;
t625 = rSges(4,1) * t381;
t358 = t369 * rSges(4,3);
t317 = rSges(5,2) * t578;
t507 = rSges(5,1) * t577;
t418 = rSges(5,3) * t368 + t507;
t195 = -t317 + t418;
t422 = -t195 + t456;
t461 = t265 * t287 + t324;
t72 = qJD(1) * t422 + t461;
t620 = t369 * t72;
t498 = rSges(5,1) * t502 + rSges(5,2) * t695;
t118 = -qJD(1) * t418 + t498;
t508 = rSges(5,1) * t581;
t193 = -t508 + t538;
t239 = t672 * t376;
t414 = (-qJDD(1) * t380 - t571) * pkin(1);
t37 = qJD(1) * t118 + qJDD(1) * t193 - t177 * t287 - t239 * t266 + t395 + t414;
t619 = t37 * t369;
t617 = qJDD(1) / 0.2e1;
t615 = t375 - rSges(6,3);
t588 = t201 * t379;
t587 = t203 * t381;
t318 = -t637 - t638;
t584 = t318 * t368;
t242 = t329 * t368;
t583 = t329 * t369;
t569 = t566 * t369;
t561 = t369 * t199 + t201 * t580;
t560 = t368 * t199 + t203 * t575;
t537 = t368 * t286 + t319;
t531 = rSges(4,2) * t580 + t358;
t522 = qJD(1) * t330;
t424 = t331 * t379 - t381 * t671;
t136 = -t369 * t424 + t242;
t512 = t136 * qJD(1);
t347 = pkin(3) * t580;
t505 = pkin(3) * t514;
t504 = qJD(1) * t633;
t497 = rSges(4,1) * t494 + (t368 * t514 + t379 * t520) * rSges(4,2);
t492 = t369 * t515;
t490 = -pkin(2) - t625;
t486 = -t521 / 0.2e1;
t485 = t520 / 0.2e1;
t484 = -t517 / 0.2e1;
t483 = t517 / 0.2e1;
t482 = -t516 / 0.2e1;
t481 = t516 / 0.2e1;
t467 = -t200 - t587;
t466 = qJD(1) * t225 - t266 * t672;
t459 = qJD(1) * t227 + t265 * t672;
t457 = -pkin(4) * t573 - t238;
t346 = rSges(4,2) * t576;
t419 = -rSges(4,1) * t575 - rSges(4,3) * t368;
t208 = -t346 - t419;
t455 = -t208 + t479;
t454 = -qJD(1) * t452 + t504;
t449 = t172 - t489;
t338 = rSges(2,1) * t382 - t380 * rSges(2,2);
t448 = rSges(2,1) * t380 + rSges(2,2) * t382;
t447 = rSges(3,1) * t368 + rSges(3,2) * t369;
t337 = -rSges(4,2) * t379 + t625;
t335 = rSges(4,1) * t379 + rSges(4,2) * t381;
t430 = -t587 + t588;
t400 = qJD(3) * t242 + (-t330 * t369 + t430 - t595) * qJD(1);
t401 = qJD(1) * t429 - qJD(3) * t583 + t525;
t444 = (t400 * t368 - t369 * t664) * t369 + (t401 * t368 - t369 * t663) * t368;
t443 = (t368 * t664 + t400 * t369) * t369 + (t368 * t663 + t401 * t369) * t368;
t81 = -t203 * t579 + t561;
t82 = t369 * t200 + t202 * t580 - t204 * t579;
t442 = t368 * t82 + t369 * t81;
t83 = -t201 * t576 + t560;
t84 = t200 * t368 - t684;
t441 = t368 * t84 + t369 * t83;
t207 = -rSges(4,1) * t579 + t531;
t101 = -qJD(1) * t207 + t335 * t516 + t454;
t257 = t335 * t517;
t102 = qJD(1) * t455 + t257;
t434 = t101 * t369 + t102 * t368;
t249 = t335 * t369;
t137 = -qJD(3) * t249 + (-t337 * t368 + t358) * qJD(1);
t138 = qJD(1) * t419 + t497;
t433 = t137 * t369 - t138 * t368;
t428 = -t207 * t368 + t208 * t369;
t423 = -t509 - t361;
t413 = t447 + t633;
t397 = t424 * qJD(1) + t330 * qJD(3);
t396 = pkin(3) * t492 + qJD(1) * t412 + t454;
t390 = -t172 * t516 + t412 * t517 + qJD(2);
t387 = qJDD(2) + t141 * t516 + (-qJD(3) * t142 - t258 * t629) * t368 - t259 * t172 + t258 * t679;
t386 = (t368 * t710 + t369 * t711) * t650 + (t368 * t712 + t369 * t713) * t649 + (t696 * t368 + t730 * t369) * t646 + (t723 * t369 + t722 * t368 + (-t368 * t711 + t369 * t710) * qJD(1)) * t645 + (-t730 * t368 + t696 * t369) * t644 + (t721 * t369 + t720 * t368 + (-t368 * t713 + t369 * t712) * qJD(1)) * t643 + (qJD(1) * t715 + qJDD(1) * t706 + t176 * t710 + t177 * t711 + t265 * t722 + t266 * t723) * t642 + (qJD(1) * t714 + qJDD(1) * t707 + t176 * t712 + t177 * t713 + t265 * t720 + t266 * t721) * t641 + (t697 * t370 + t737 * t371) * t631 + (t717 * t369 + t716 * t368 + (-t368 * t709 + t369 * t708) * qJD(1)) * t630 + (t368 * t708 + t369 * t709) * t617 + t719 * t486 + t718 * t485;
t349 = rSges(3,2) * t521;
t298 = t337 * qJD(3);
t255 = t369 * t318;
t248 = t335 * t368;
t228 = t447 * qJD(1) + t504;
t206 = pkin(3) * t576 + t255;
t205 = -t347 - t584;
t165 = t369 * t195;
t155 = t369 * t172;
t135 = t368 * t424 + t583;
t125 = t135 * qJD(1);
t96 = qJD(3) * t428 + qJD(2);
t71 = -qJD(1) * t193 + t266 * t287 + t396;
t70 = t368 * t662 + t397 * t369;
t69 = t397 * t368 - t369 * t662;
t66 = -t265 * t193 + t266 * t195 + t390;
t65 = t298 * t517 + t258 * t335 + t365 + (-t137 - t256) * qJD(1) + t455 * qJDD(1);
t64 = qJD(1) * t138 + qJDD(1) * t207 - t259 * t335 - t298 * t516 + t414 + t545;
t63 = -qJD(3) * t429 + t131 * t381 + t133 * t379;
t62 = -qJD(3) * t430 + t132 * t381 + t134 * t379;
t61 = qJD(3) * t433 - t207 * t258 + t208 * t259 + qJDD(2);
t60 = qJD(1) * t416 - t665;
t59 = qJD(1) * t567 + t266 * t477 - t353 + t396;
t49 = qJD(3) * t441 + t512;
t48 = qJD(3) * t442 + t125;
t43 = t265 * t567 - t266 * t566 + t390;
t38 = t176 * t287 + t239 * t265 + (-t116 + t568) * qJD(1) + t422 * qJDD(1) + t500;
t18 = t266 * t116 - t265 * t118 - t176 * t193 + t177 * t195 + t387;
t9 = t176 * t567 - t177 * t566 - t265 * t616 + t266 * t614 + t387;
t1 = [(t125 + ((t561 + t84 + t684) * t369 + (-t83 + (t467 - t588) * t369 + t82 + t560) * t368) * qJD(3)) * t484 + (t136 + t122) * t648 + (t135 + t121) * t647 + ((t698 + t710) * t266 + t726) * t646 + (t49 - t512 + ((t82 + (-t200 + t588) * t369 - t560) * t369 + (t368 * t467 + t561 - t81) * t368) * qJD(3)) * t482 + (t70 + t62) * t481 + (-qJD(3) * t424 + t292 * t381 + t293 * t379 - t751 * t370 + t750 * t371) * qJD(1) + (-t228 * t349 + (qJD(1) * t253 * t413 - t228 * t450) * qJD(1) + (t385 * t447 - g(2) + t365) * t253 + (t510 + (t626 * qJD(1) - t349) * qJD(1) + g(3)) * t413) * m(3) + (g(2) * t338 + g(3) * t448) * m(2) + (t63 + t69 + t48) * t483 + (-t60 * t353 - t59 * t668 + t60 * (t286 * t376 - t261) * t369 + ((t380 * t60 + t382 * t59) * pkin(1) + (-t423 * t59 + t60 * t615) * t369 + (t60 * (-t423 - t623) + t59 * rSges(6,3)) * t368) * qJD(1) - (t60 + (-t566 + t632) * qJD(1) + t665 + t677) * t59 + t724 * (t368 * t615 - t503 - t632 + t728) + t725 * (t368 * t423 - t369 * t375 + t539 - t633)) * m(6) + ((t287 * t376 + t506) * t620 + (t540 + (t508 + t694) * qJD(1)) * t72 + (-t498 - t536 - t72 + t461 - t677 + (t418 - t676 - t195 - t632) * qJD(1)) * t71 + (t38 - g(2)) * (-t507 + t317 + (-rSges(5,3) + t383) * t368 + t676) + (t37 - g(3)) * (-t369 * t383 + (-t364 - t362) * t368 - t694)) * m(5) + (-(t102 - t257 + t264 + (t208 + t632) * qJD(1)) * t101 + t102 * (rSges(4,1) * t492 + rSges(4,2) * t491 + t350) - t101 * t497 + ((t101 * t382 + t102 * t380) * pkin(1) + (-t101 * t490 + t102 * t640) * t369 + (-t101 * t640 + t102 * t337) * t368) * qJD(1) + (t65 - g(2)) * (t640 * t368 + t490 * t369 + t346 - t632) + (t64 - g(3)) * (t368 * t490 + t531 - t633 + t635)) * m(4) + (t706 + t708) * t650 + (t707 + t709) * t649 + (((-t758 + t743) * t369 + t752 * t368 + t712 + t741) * t266 + (t698 - t713) * t265 + t718 - t727) * t644 + (t714 + t717) * t643 + (m(3) * (t253 ^ 2 + t413 ^ 2) + t425 + m(2) * (t338 ^ 2 + t448 ^ 2) + Icges(2,3) + Icges(3,3) + t769 * t371 + t771 * t370) * qJDD(1) + (t715 + t716 + t719) * t645; m(3) * qJDD(2) + m(4) * t61 + m(5) * t18 + m(6) * t9 + (-m(3) - m(4) - m(5) - m(6)) * g(1); (t368 * t63 + t369 * t62 + (-t121 * t368 + t122 * t369) * qJD(1)) * t630 + (qJD(1) * t70 + qJD(3) * t443 + qJDD(1) * t135 + t258 * t82 + t259 * t81) * t641 + t441 * t648 + t442 * t647 + (t121 * t369 + t122 * t368) * t617 + t386 + ((-t81 * t368 + t82 * t369) * qJD(1) + t443) * t481 + ((-t379 * t533 + t381 * t532) * qJD(1) + ((t368 * t547 + t369 * t548) * t381 + (-t368 * t549 - t369 * t550) * t379) * qJD(3)) * t631 + ((-t517 * t583 + t522) * t368 + (-t690 + (t659 * t369 + (t242 - t660) * t368) * qJD(3)) * t369) * t484 + ((-t83 * t368 + t84 * t369) * qJD(1) + t444) * t483 + (qJD(1) * t69 + qJD(3) * t444 + qJDD(1) * t136 + t258 * t84 + t259 * t83) * t642 + ((t242 * t516 + t522) * t369 + (t690 + (t660 * t368 + (-t583 - t659) * t369) * qJD(3)) * t368) * t482 + t48 * t486 + t49 * t485 + (-g(1) * (t374 + t453) - g(2) * (t224 - t584) - g(3) * (t255 - t226) + t9 * (-t155 + (t412 + t567) * t368 - t569) + t13 * (t347 + t537) + (-t286 + t318) * t689 + (qJD(1) * t205 - (t457 - t505) * t369 - t669 + t670) * t59 + (-t206 * t266 - (-t205 - t224) * t265 + (-t142 - t616) * t368 + ((t679 + t567) * t369 + (t449 + t566) * t368) * qJD(1) + t666 + t678) * t43 + (-(-t206 + t226) * qJD(1) + t667) * t60) * m(6) + (-g(1) * (t672 + t374) - g(2) * (t347 + t225) + t18 * (-t155 + t165 + (-t193 + t412) * t368) + t38 * (t287 * t368 + t347) + (t466 - t287 * t521 - (-t239 - t505) * t369 - t669) * t71 + ((-t142 - t118) * t368 + ((-t193 + t679) * t369 + (-t195 + t449) * t368) * qJD(1) + t678 + t686) * t66 + (-t634 + t619) * (-t287 - t638) + (t239 * t368 + t287 * t520 - t459) * t72) * m(5) + (-g(1) * t337 - g(2) * t248 + g(3) * t249 - (-t101 * t248 + t102 * t249) * qJD(1) - (t96 * (-t248 * t368 - t249 * t369) + t434 * t337) * qJD(3) + t61 * t428 + t96 * ((-t207 * t369 - t208 * t368) * qJD(1) + t433) + t434 * t298 + (t65 * t368 - t64 * t369 + (-t101 * t368 + t102 * t369) * qJD(1)) * t335) * m(4); t386 + (-g(1) * t453 - g(2) * (t319 + t224) - (-t618 + (-rSges(6,1) - pkin(4)) * t370) * t634 + t13 * t537 - t477 * t689 + (t368 * t567 - t569) * t9 + (-qJD(1) * t226 - t297 + t667) * t60 + (-t369 * t457 + t296 + t670) * t59 + (t224 * t265 - (-t265 * t368 - t266 * t369) * t637 - t616 * t368 + (t566 * t368 + t567 * t369) * qJD(1) + t666) * t43) * m(6) + (-g(1) * t672 - g(2) * t225 + g(3) * t227 - t459 * t72 + t466 * t71 + t18 * (-t368 * t193 + t165) + (t368 * t72 + t369 * t71) * t239 + (t38 * t368 - t619 + (-t368 * t71 + t620) * qJD(1)) * t287 + (-t368 * t118 + (-t369 * t193 - t368 * t195) * qJD(1) + t686) * t66) * m(5); ((t265 * t43 + t724) * t369 + (-t266 * t43 + t725) * t368) * m(6);];
tau = t1;
