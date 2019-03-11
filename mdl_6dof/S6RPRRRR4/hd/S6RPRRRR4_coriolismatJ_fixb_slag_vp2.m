% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:58
% EndTime: 2019-03-09 07:05:25
% DurationCPUTime: 18.04s
% Computational Cost: add. (58450->606), mult. (112889->797), div. (0->0), fcn. (140724->10), ass. (0->367)
t356 = sin(qJ(6));
t354 = sin(pkin(11));
t355 = cos(pkin(11));
t618 = sin(qJ(3));
t621 = cos(qJ(3));
t333 = -t618 * t354 + t621 * t355;
t413 = t354 * t621 + t355 * t618;
t620 = cos(qJ(4));
t397 = t620 * t413;
t617 = sin(qJ(4));
t299 = t617 * t333 + t397;
t357 = sin(qJ(5));
t619 = cos(qJ(5));
t396 = t617 * t413;
t701 = t620 * t333 - t396;
t742 = t299 * t619 + t357 * t701;
t690 = mrSges(7,1) * t742;
t264 = t299 * t357 - t619 * t701;
t358 = cos(qJ(6));
t771 = t264 * t358;
t783 = mrSges(7,3) * t771 + t690;
t798 = t356 * t783;
t803 = -t798 / 0.2e1;
t604 = pkin(7) + qJ(2);
t446 = t618 * t604;
t447 = t621 * t604;
t279 = (-pkin(8) * t618 - t446) * t355 + (-pkin(8) * t621 - t447) * t354;
t303 = t333 * t604;
t381 = t333 * pkin(8) + t303;
t670 = t620 * t279 - t617 * t381;
t708 = t670 * mrSges(5,2);
t208 = t617 * t279 + t381 * t620;
t731 = t208 * mrSges(5,1);
t732 = Ifges(5,5) * t701;
t585 = Ifges(7,6) * t358;
t593 = Ifges(7,5) * t356;
t337 = t585 + t593;
t755 = t742 * t337;
t757 = Ifges(6,6) * t742;
t442 = mrSges(7,1) * t358 - t356 * mrSges(7,2);
t364 = pkin(9) * t701 + t208;
t743 = -t299 * pkin(9) + t670;
t749 = t357 * t743 + t364 * t619;
t764 = t749 * t442;
t598 = Ifges(7,4) * t358;
t439 = -Ifges(7,2) * t356 + t598;
t587 = Ifges(7,6) * t742;
t696 = -t264 * t439 + t587;
t767 = t358 * t696;
t599 = Ifges(7,4) * t356;
t441 = Ifges(7,1) * t358 - t599;
t594 = Ifges(7,5) * t742;
t695 = -t264 * t441 + t594;
t768 = t356 * t695;
t774 = t749 * mrSges(6,1);
t89 = t357 * t364 - t619 * t743;
t791 = t89 * mrSges(6,2);
t738 = t791 + t755 / 0.2e1 - t757 - t764 + t767 / 0.2e1 + t768 / 0.2e1 - t774;
t748 = Ifges(5,6) * t299;
t775 = Ifges(6,5) * t264;
t802 = t738 + t732 - t708 - t731 - t748 - t775;
t339 = Ifges(7,1) * t356 + t598;
t515 = t358 * t339;
t338 = t358 * Ifges(7,2) + t599;
t523 = t356 * t338;
t801 = t767 / 0.4e1 + t768 / 0.4e1 + t755 / 0.4e1 - t757 / 0.2e1 - t264 * (t515 / 0.4e1 - t523 / 0.4e1);
t745 = t791 / 0.2e1 - t764 / 0.2e1 - t774 / 0.2e1;
t800 = t745 - t775 / 0.2e1 + t801;
t689 = mrSges(7,2) * t742;
t772 = t264 * t356;
t784 = mrSges(7,3) * t772 - t689;
t797 = t358 * t784;
t799 = t797 - t798;
t491 = -t355 * pkin(2) - pkin(1);
t309 = -t333 * pkin(3) + t491;
t273 = -pkin(4) * t701 + t309;
t796 = m(6) * t273 + mrSges(6,1) * t264 + mrSges(6,2) * t742;
t728 = t357 * t89;
t795 = -t619 * t749 - t728;
t113 = t264 * pkin(5) - pkin(10) * t742 + t273;
t53 = t113 * t358 - t356 * t749;
t54 = t356 * t113 + t358 * t749;
t566 = t358 * mrSges(7,2);
t568 = t356 * mrSges(7,1);
t336 = t566 + t568;
t168 = t336 * t742;
t770 = t336 * t264;
t715 = t749 * t168 - t89 * t770;
t793 = t53 * t783 + t54 * t784 + t715;
t348 = pkin(3) * t620 + pkin(4);
t449 = t619 * t617;
t323 = pkin(3) * t449 + t357 * t348;
t730 = t323 * t89;
t729 = t356 * t89;
t727 = t358 * t89;
t773 = t749 * t89;
t511 = t619 * pkin(4);
t347 = -t511 - pkin(5);
t626 = t347 / 0.2e1;
t790 = t770 * t626;
t661 = pkin(5) / 0.2e1;
t789 = t770 * t661;
t416 = t439 * t742;
t417 = t441 * t742;
t567 = t356 * mrSges(7,3);
t495 = t567 / 0.2e1;
t496 = -t567 / 0.2e1;
t112 = Ifges(7,5) * t264 + t417;
t520 = t358 * t112;
t109 = Ifges(7,6) * t264 + t416;
t531 = t356 * t109;
t537 = t742 * t358;
t538 = t742 * t356;
t553 = t89 * t336;
t629 = -t339 / 0.4e1;
t641 = t264 / 0.4e1;
t349 = Ifges(7,5) * t358;
t586 = Ifges(7,6) * t356;
t678 = t349 - t586;
t386 = t520 / 0.4e1 + 0.2e1 * t538 * t629 - t338 * t537 / 0.2e1 - t356 * t416 / 0.4e1 + t358 * t417 / 0.4e1 + t678 * t641 + t553 / 0.2e1 - t531 / 0.4e1 + (t495 + t496) * t54;
t680 = (t349 / 0.2e1 - t586 / 0.2e1) * t264;
t788 = t386 - t680;
t718 = -m(7) * t749 + t770;
t642 = t264 / 0.2e1;
t759 = -t742 / 0.2e1;
t692 = Ifges(6,4) * t759;
t758 = t742 / 0.2e1;
t428 = t273 * mrSges(6,1) + t678 * t758 + t692 + (Ifges(6,2) + Ifges(7,3)) * t642;
t506 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t778 = -t264 / 0.2e1;
t782 = Ifges(6,1) * t778 + t264 * t506 + t428;
t781 = Ifges(7,5) * t778 - t112 / 0.4e1;
t780 = t109 / 0.4e1 + (t642 + t641) * Ifges(7,6);
t776 = pkin(5) * t749;
t710 = Ifges(7,3) * t758;
t578 = t264 * mrSges(6,3);
t769 = t347 * t749;
t614 = pkin(5) * t742;
t180 = pkin(10) * t264 + t614;
t295 = Ifges(5,4) * t701;
t440 = -Ifges(5,2) * t299 + t295;
t760 = t440 / 0.2e1;
t756 = t742 * mrSges(6,3);
t726 = t701 * mrSges(5,2);
t747 = t299 * mrSges(5,1);
t463 = t726 + t747;
t753 = -t708 / 0.2e1 - t731 / 0.2e1 + t732 / 0.2e1 - t748 / 0.2e1;
t752 = -t726 / 0.2e1 - t747 / 0.2e1;
t750 = t701 * t295 / 0.2e1 + (-Ifges(5,4) * t299 + (-Ifges(5,2) / 0.2e1 + Ifges(5,1)) * t701) * t299;
t616 = pkin(4) * t299;
t419 = t358 * t439;
t420 = t356 * t441;
t422 = t523 / 0.2e1 - t515 / 0.2e1;
t452 = t419 / 0.2e1 + t420 / 0.2e1 - t422;
t622 = t358 / 0.2e1;
t625 = -t356 / 0.2e1;
t746 = t622 * t695 + t625 * t696 + t692;
t353 = t358 ^ 2;
t602 = mrSges(7,3) * t353;
t352 = t356 ^ 2;
t603 = mrSges(7,3) * t352;
t719 = -t603 / 0.2e1 - t602 / 0.2e1;
t401 = t413 * pkin(3);
t735 = m(5) * t401;
t482 = t617 * t357;
t322 = -pkin(3) * t482 + t348 * t619;
t316 = -pkin(5) - t322;
t634 = t316 / 0.2e1;
t720 = t770 * t634;
t711 = mrSges(6,1) * t758;
t513 = t352 + t353;
t671 = mrSges(7,3) * t513;
t707 = -mrSges(6,2) + t671;
t706 = t316 * t749;
t117 = t180 + t616;
t63 = t117 * t358 + t729;
t64 = t356 * t117 - t727;
t434 = -t63 * t356 + t64 * t358;
t317 = pkin(10) + t323;
t328 = (t619 * t620 - t482) * pkin(3);
t565 = t358 * mrSges(7,3);
t493 = t565 / 0.2e1;
t390 = t64 * t493 + t63 * t496 + t800;
t510 = mrSges(7,3) * t538;
t172 = -mrSges(7,2) * t264 - t510;
t516 = t358 * t172;
t176 = t264 * mrSges(7,1) - mrSges(7,3) * t537;
t525 = t356 * t176;
t421 = t516 / 0.2e1 - t525 / 0.2e1;
t564 = t358 * t54;
t436 = -t356 * t53 + t564;
t170 = t264 * t567 - t689;
t518 = t358 * t170;
t174 = t264 * t565 + t690;
t527 = t356 * t174;
t327 = (t357 * t620 + t449) * pkin(3);
t570 = t327 * t89;
t631 = t327 / 0.2e1;
t632 = -t322 / 0.2e1;
t662 = m(7) / 0.2e1;
t664 = m(6) / 0.2e1;
t699 = (t317 * t434 + t328 * t436 + t570 + t706) * t662 + (-t730 + t570 + (-t322 + t328) * t749) * t664 + t390 + t168 * t631 - t720 + (t518 / 0.2e1 - t527 / 0.2e1) * t317 + (-t264 * t632 + t323 * t759 + t328 * t778 + t631 * t742) * mrSges(6,3) + t421 * t328 + t753;
t687 = mrSges(7,3) * (t353 / 0.2e1 + t352 / 0.2e1);
t683 = -mrSges(6,1) - t442;
t676 = t797 / 0.2e1 + t803;
t623 = -t358 / 0.2e1;
t423 = t172 * t625 + t176 * t623;
t69 = t180 * t358 + t729;
t70 = t356 * t180 - t727;
t433 = -t69 * t356 + t70 * t358;
t275 = t401 + t616;
t114 = t180 + t275;
t57 = t114 * t358 + t729;
t58 = t356 * t114 - t727;
t435 = -t57 * t356 + t58 * t358;
t674 = Ifges(6,1) * t758 - Ifges(6,4) * t264 - t531 / 0.2e1 + t520 / 0.2e1;
t672 = (-t617 * mrSges(5,1) - t620 * mrSges(5,2)) * pkin(3);
t668 = t719 * t742 + t423;
t663 = -m(7) / 0.2e1;
t660 = -pkin(10) / 0.2e1;
t659 = pkin(10) / 0.2e1;
t658 = m(5) * pkin(3);
t657 = m(6) * pkin(4);
t656 = m(7) * pkin(4);
t655 = -mrSges(7,2) / 0.2e1;
t654 = t57 / 0.2e1;
t653 = -t58 / 0.2e1;
t652 = -t69 / 0.2e1;
t651 = t70 / 0.2e1;
t260 = t264 * mrSges(6,2);
t648 = t260 / 0.2e1;
t635 = -t316 / 0.2e1;
t633 = t317 / 0.2e1;
t630 = t338 / 0.4e1;
t615 = pkin(4) * t357;
t346 = pkin(10) + t615;
t628 = t346 / 0.2e1;
t627 = -t347 / 0.2e1;
t624 = t356 / 0.2e1;
t613 = pkin(5) * t336;
t580 = t742 * mrSges(6,1);
t577 = t273 * mrSges(6,2);
t375 = t642 * t678 + t674;
t377 = t413 ^ 2;
t395 = mrSges(4,1) * t413 + t333 * mrSges(4,2);
t504 = mrSges(5,3) * t620;
t3 = m(7) * (t53 * t57 + t54 * t58 + t773) + t491 * t395 + t750 - Ifges(4,4) * t377 - t273 * t260 + (Ifges(4,4) * t333 + t620 * t760 - pkin(3) * mrSges(5,1) * t397 + (-Ifges(4,2) + Ifges(4,1)) * t413) * t333 - t440 * t396 / 0.2e1 + (mrSges(5,1) * t396 + t299 * mrSges(5,2)) * t401 + (t463 + t735) * t309 + t58 * t172 + t57 * t176 + ((-t701 - t396) * mrSges(5,3) + t504 * t333) * t670 + ((-Ifges(7,1) * t771 + t594) * t622 + (Ifges(7,2) * t772 + t587) * t625 + t692 + t772 * t598 + t782) * t742 - t375 * t264 + t796 * t275 + t793;
t574 = t3 * qJD(1);
t573 = t322 * mrSges(6,2);
t572 = t323 * mrSges(6,1);
t571 = t327 * mrSges(6,1);
t569 = t328 * mrSges(6,2);
t4 = t701 * t760 + t309 * t463 + t54 * t170 + t64 * t172 + t53 * t174 + t63 * t176 + m(7) * (t53 * t63 + t54 * t64 + t773) + (t746 + t782) * t742 + (-t375 - t577) * t264 + t715 + t796 * t616 + t750;
t563 = t4 * qJD(1);
t554 = t89 * t742;
t9 = t69 * t176 + t70 * t172 + m(7) * (t53 * t69 + t54 * t70 + t773) + (t428 + t746) * t742 + (-t577 - t680 + (-Ifges(6,1) / 0.2e1 + t506) * t742 - t674) * t264 + t793;
t552 = t9 * qJD(1);
t302 = -t354 * t447 - t355 * t446;
t17 = (t168 + t756) * t742 + (t299 ^ 2 + t701 ^ 2) * mrSges(5,3) + (t333 ^ 2 + t377) * mrSges(4,3) - (t516 - t525 - t578) * t264 + m(7) * (-t264 * t436 + t554) + m(6) * (-t264 * t749 + t554) + m(5) * (t208 * t701 - t299 * t670) + m(4) * (-t302 * t413 + t303 * t333) + (m(3) * qJ(2) + mrSges(3,3)) * (t354 ^ 2 + t355 ^ 2);
t550 = qJD(1) * t17;
t418 = t442 * t742;
t10 = t54 * t176 - t89 * t418 + (mrSges(7,3) * t564 + t109 * t622 + t112 * t624 + t337 * t642 - t422 * t742) * t742 + (-t510 - t172) * t53;
t549 = t10 * qJD(1);
t481 = t442 * t758;
t389 = t264 * t719 - t481;
t461 = t513 * t264;
t366 = (-t264 * t323 - t322 * t742) * t664 + (t316 * t742 - t317 * t461) * t662 + (-t299 * t620 + t617 * t701) * t658 / 0.2e1 + t389;
t367 = t275 * t664 + (t356 * t58 + t358 * t57) * t662 + t784 * t624 + t783 * t622 + t735 / 0.2e1;
t20 = t260 + t366 - t367 - t395 - t463 - t580;
t547 = t20 * qJD(1);
t512 = t657 / 0.2e1;
t363 = (t356 * t64 + t358 * t63) * t662 + t711 + mrSges(6,2) * t778 + t170 * t624 + t174 * t622 + t299 * t512 - t752;
t458 = t513 * t346;
t497 = -t580 / 0.2e1;
t365 = t648 + (-t264 * t458 + t347 * t742) * t662 + t497 + (-t264 * t357 - t619 * t742) * t512 + t389 + t752;
t21 = -t363 + t365;
t546 = t21 * qJD(1);
t373 = -t264 * t687 + t648 + (-pkin(10) * t461 - t614) * t662 - t481;
t376 = (t356 * t70 + t358 * t69) * t663 + mrSges(6,2) * t642 + t784 * t625 + t783 * t623;
t25 = 0.2e1 * t759 * mrSges(6,1) + t373 + t376;
t545 = t25 * qJD(1);
t425 = -t566 / 0.2e1 - t568 / 0.2e1;
t411 = t425 * t264;
t28 = -t411 - t421;
t536 = t28 * qJD(1);
t535 = t316 * t336;
t534 = t323 * t442;
t533 = t327 * t442;
t532 = t347 * t336;
t509 = -t615 / 0.2e1;
t508 = t615 / 0.2e1;
t507 = Ifges(7,2) / 0.4e1 - Ifges(7,1) / 0.4e1;
t494 = -t565 / 0.2e1;
t492 = t337 / 0.4e1 - Ifges(6,6) / 0.2e1;
t490 = t356 * t619;
t489 = t358 * t619;
t460 = t513 * t322;
t459 = t513 * t328;
t456 = -t511 / 0.2e1;
t455 = t511 / 0.2e1;
t448 = -t490 / 0.2e1;
t437 = t706 + t730;
t432 = -t613 / 0.2e1 + t452;
t362 = t662 * t769 - t790 + t795 * t512 + t509 * t756 - t456 * t578 + (t435 * t662 + t676) * t346 + t753;
t2 = t58 * t494 + t57 * t495 - t362 + t699 - t800;
t76 = t683 * t327 + t672 + t707 * t328 + m(7) * (t316 * t327 + t317 * t459) + m(6) * (-t322 * t327 + t323 * t328);
t431 = t2 * qJD(1) + t76 * qJD(3);
t382 = (pkin(10) * t435 - t776) * t663 - t789;
t387 = t492 * t742 - t720 + t323 * t168 / 0.2e1;
t404 = t695 / 0.4e1 - t317 * t783 / 0.2e1 + t176 * t632;
t405 = t696 / 0.4e1 + t784 * t633 + t322 * t172 / 0.2e1;
t5 = (t660 * t784 + t405) * t358 + (t659 * t783 + t404) * t356 + (t317 * t433 + t322 * t436 + t437) * t662 + (-t585 / 0.4e1 - t593 / 0.4e1 - t492) * t742 + ((t651 + t653) * t358 + (t652 + t654) * t356) * mrSges(7,3) + t382 + t387 + (t419 / 0.4e1 + t420 / 0.4e1) * t264;
t383 = mrSges(7,3) * t460 - t534 - t572 - t573;
t74 = m(7) * (t316 * t323 + t317 * t460) + t383;
t430 = t5 * qJD(1) + t74 * qJD(3);
t398 = t598 - t507 * t356 + t339 / 0.4e1;
t399 = mrSges(7,1) * t654 + mrSges(7,2) * t653 + t710;
t415 = t358 * t507 + t630;
t424 = -t264 * t349 / 0.4e1 - t553 / 0.2e1;
t11 = (t172 * t633 + t780) * t356 + (t176 * t633 + t781) * t358 + (t317 * t687 + (mrSges(7,1) * t635 + t415) * t358 + (mrSges(7,2) * t634 + t398) * t356) * t742 + t399 + t424;
t226 = t452 + t535;
t429 = -t11 * qJD(1) + t226 * qJD(3);
t427 = mrSges(7,1) * t652 + mrSges(7,2) * t651;
t414 = -pkin(5) * t442 / 0.2e1;
t412 = t513 * t619;
t410 = t425 * t322;
t409 = t425 * t328;
t408 = t418 / 0.2e1;
t372 = t511 * t707 + t615 * t683;
t244 = (t346 * t412 + t347 * t357) * t656 + t372;
t361 = (t347 * t323 + (t316 * t357 + t317 * t412) * pkin(4)) * t662 - t573 / 0.2e1 - t572 / 0.2e1 - t534 / 0.2e1 + mrSges(6,1) * t509 - t442 * t508 + mrSges(6,2) * t456 + t455 * t671 + (t458 * t662 - t719) * t322;
t368 = (-pkin(5) * t327 + pkin(10) * t459) * t663 + t571 / 0.2e1 + t533 / 0.2e1 + t569 / 0.2e1 + t719 * t328;
t48 = t361 + t368;
t384 = Ifges(6,5) * t778 + t745;
t359 = (t769 + t433 * t346 + (t489 * t54 - t490 * t53 + t728) * pkin(4)) * t662 - t790 + t346 * t803 + t168 * t508 + t797 * t628 + t69 * t496 + t70 * t493 + pkin(4) * t176 * t448 + t455 * t516 + t384 + t801;
t370 = (pkin(10) * t434 - t776) * t663 - t789 + t527 * t659 + t518 * t660;
t7 = t64 * t494 + t63 * t495 + t359 + t370 - t800;
t406 = t7 * qJD(1) + t48 * qJD(3) + t244 * qJD(4);
t400 = t710 + t63 * mrSges(7,1) / 0.2e1 + t64 * t655;
t13 = (t172 * t628 + t780) * t356 + (t176 * t628 + t781) * t358 + (t346 * t687 + (mrSges(7,1) * t627 + t415) * t358 + (mrSges(7,2) * t626 + t398) * t356) * t742 + t400 + t424;
t182 = (t627 + t635) * t336 + t409 - t452;
t270 = t452 + t532;
t402 = -t13 * qJD(1) - t182 * qJD(3) + t270 * qJD(4);
t16 = t386 + (t414 - Ifges(7,3) / 0.2e1 - pkin(10) * t687) * t742 + t423 * pkin(10) + t680 + t427;
t184 = (t661 + t635) * t336 + t410 - t452;
t385 = (mrSges(7,1) * t448 + t489 * t655) * pkin(4);
t191 = (t661 + t627) * t336 + t385 - t452;
t271 = -t452 + t613;
t393 = t16 * qJD(1) - t184 * qJD(3) - t191 * qJD(4) - t271 * qJD(5);
t388 = t58 * t493 + t57 * t496 + t800;
t320 = t532 / 0.2e1;
t300 = t535 / 0.2e1;
t192 = t320 + t385 + t432;
t185 = t300 + t410 + t432;
t183 = t320 + t300 + t409 + t452;
t47 = t361 - t368;
t29 = -t411 + t421;
t26 = t373 - t376 + t497 + t711;
t24 = t366 + t367;
t22 = t363 + t365;
t15 = t386 + t742 * t414 + (-t687 * t742 + t423) * pkin(10) - Ifges(7,5) * t771 / 0.2e1 + Ifges(7,6) * t772 / 0.2e1 + t710 - t427;
t14 = t346 * t668 + t347 * t408 + t400 + t788;
t12 = t316 * t408 + t317 * t668 + t399 + t788;
t8 = t359 - t370 + t390;
t6 = t437 * t662 + ((t317 * t70 + t322 * t54) * t662 + mrSges(7,3) * t651 + t264 * t629 + t405) * t358 + ((-t317 * t69 - t322 * t53) * t662 + mrSges(7,3) * t652 + t264 * t630 + t404) * t356 + t387 + t384 - t382 + t388 + t676 * pkin(10);
t1 = t362 + t388 + t699;
t18 = [qJD(2) * t17 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t9 - qJD(6) * t10, qJD(3) * t24 + qJD(4) * t22 + qJD(5) * t26 + qJD(6) * t29 + t550, t574 + t24 * qJD(2) + ((-mrSges(5,3) * t299 * t617 - t504 * t701) * pkin(3) + t322 * t578 - t718 * t316 + m(6) * (-t322 * t749 - t730) + t422 * t264 - t323 * t756 + (m(7) * t435 + t799) * t317 - Ifges(4,6) * t413 + t435 * mrSges(7,3) + Ifges(4,5) * t333 - t302 * mrSges(4,2) - t303 * mrSges(4,1) + (-t208 * t620 + t617 * t670) * t658 + t802) * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + t12 * qJD(6), t22 * qJD(2) + t1 * qJD(3) + t8 * qJD(5) + t14 * qJD(6) + t563 + (t795 * t657 - t615 * t756 - (-mrSges(6,3) * t511 - t422) * t264 - t718 * t347 + (m(7) * t434 + t518 - t527) * t346 + t434 * mrSges(7,3) + t802) * qJD(4), t26 * qJD(2) + t6 * qJD(3) + t8 * qJD(4) + t15 * qJD(6) + t552 + ((-Ifges(6,5) + t422) * t264 + t718 * pkin(5) + (m(7) * t433 + t799) * pkin(10) + t433 * mrSges(7,3) + t738) * qJD(5), -t549 + t29 * qJD(2) + t12 * qJD(3) + t14 * qJD(4) + t15 * qJD(5) + (-t54 * mrSges(7,1) - t53 * mrSges(7,2) - t755) * qJD(6); -qJD(3) * t20 - qJD(4) * t21 - qJD(5) * t25 - qJD(6) * t28 - t550, 0, -t547, -t546, -t545, -qJD(6) * t336 - t536; qJD(2) * t20 + qJD(4) * t2 + qJD(5) * t5 - qJD(6) * t11 - t574, t547, qJD(4) * t76 + qJD(5) * t74 + qJD(6) * t226 (-t533 - t569 - t571 + (m(7) * t347 - t619 * t657) * t327 + (m(7) * t458 + t357 * t657 + t602 + t603) * t328 + t672) * qJD(4) + t47 * qJD(5) + t183 * qJD(6) + t431, t47 * qJD(4) + (m(7) * (-pkin(5) * t323 + pkin(10) * t460) + t383) * qJD(5) + t185 * qJD(6) + t430, t183 * qJD(4) + t185 * qJD(5) + (-t317 * t442 + t678) * qJD(6) + t429; qJD(2) * t21 - qJD(3) * t2 + qJD(5) * t7 - qJD(6) * t13 - t563, t546, qJD(5) * t48 - qJD(6) * t182 - t431, qJD(5) * t244 + qJD(6) * t270 ((-pkin(5) * t357 + pkin(10) * t412) * t656 + t372) * qJD(5) + t192 * qJD(6) + t406, t192 * qJD(5) + (-t346 * t442 + t678) * qJD(6) + t402; qJD(2) * t25 - qJD(3) * t5 - qJD(4) * t7 + qJD(6) * t16 - t552, t545, -qJD(4) * t48 - qJD(6) * t184 - t430, -qJD(6) * t191 - t406, -t271 * qJD(6) (-pkin(10) * t442 + t678) * qJD(6) + t393; qJD(2) * t28 + qJD(3) * t11 + qJD(4) * t13 - qJD(5) * t16 + t549, t536, qJD(4) * t182 + qJD(5) * t184 - t429, qJD(5) * t191 - t402, -t393, 0;];
Cq  = t18;
