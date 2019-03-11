% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:21
% EndTime: 2019-03-08 23:31:50
% DurationCPUTime: 16.31s
% Computational Cost: add. (15834->786), mult. (36806->1074), div. (0->0), fcn. (36808->10), ass. (0->401)
t508 = sin(qJ(4));
t798 = pkin(9) - pkin(10);
t465 = t798 * t508;
t512 = cos(qJ(4));
t466 = t798 * t512;
t507 = sin(qJ(6));
t511 = cos(qJ(6));
t300 = t465 * t511 - t466 * t507;
t561 = t465 * t507 + t466 * t511;
t559 = t507 * t508 + t511 * t512;
t648 = t508 * t511;
t560 = t507 * t512 - t648;
t632 = -Ifges(7,5) * t559 + Ifges(7,6) * t560;
t825 = -t561 * mrSges(7,1) - t300 * mrSges(7,2) + t632;
t834 = qJD(6) * t825;
t766 = -m(7) / 0.2e1;
t509 = sin(qJ(3));
t513 = cos(qJ(3));
t720 = pkin(9) * t513;
t467 = pkin(3) * t509 - t720;
t613 = -pkin(8) * t508 - pkin(4);
t227 = (-pkin(10) * t513 - t467) * t512 + (-pkin(5) + t613) * t509;
t646 = t509 * t512;
t353 = -pkin(8) * t646 + t508 * t467;
t310 = t509 * qJ(5) + t353;
t647 = t508 * t513;
t250 = pkin(10) * t647 + t310;
t112 = t227 * t511 - t250 * t507;
t113 = t227 * t507 + t250 * t511;
t380 = t560 * t513;
t307 = mrSges(7,2) * t509 - mrSges(7,3) * t380;
t382 = t559 * t513;
t309 = -mrSges(7,1) * t509 - mrSges(7,3) * t382;
t655 = t467 * t512;
t315 = t509 * t613 - t655;
t649 = t508 * t509;
t352 = pkin(8) * t649 + t655;
t438 = -mrSges(6,2) * t647 + mrSges(6,3) * t509;
t758 = pkin(4) + pkin(5);
t442 = -t507 * qJ(5) - t511 * t758;
t443 = t511 * qJ(5) - t507 * t758;
t728 = t509 / 0.2e1;
t745 = t382 / 0.2e1;
t747 = -t380 / 0.2e1;
t785 = Ifges(7,5) * t745 + Ifges(7,6) * t747;
t815 = mrSges(7,2) / 0.2e1;
t816 = mrSges(7,1) / 0.2e1;
t568 = -Ifges(7,3) * t728 + t112 * t816 - t113 * t815 + t785;
t734 = -t443 / 0.2e1;
t643 = t512 * t513;
t683 = t509 * mrSges(6,1);
t437 = mrSges(6,2) * t643 - t683;
t736 = t437 / 0.2e1;
t750 = t315 / 0.2e1;
t761 = -mrSges(6,3) / 0.2e1;
t762 = mrSges(5,2) / 0.2e1;
t764 = -mrSges(5,1) / 0.2e1;
t768 = -m(6) / 0.2e1;
t833 = (-pkin(4) * t315 + qJ(5) * t310) * t768 + (t112 * t442 + t113 * t443) * t766 + pkin(4) * t736 - qJ(5) * t438 / 0.2e1 + t310 * t761 + mrSges(6,1) * t750 + t352 * t764 + t353 * t762 - t442 * t309 / 0.2e1 + t307 * t734 + t568;
t255 = -mrSges(7,1) * t560 - mrSges(7,2) * t559;
t423 = Ifges(7,4) * t560;
t260 = -Ifges(7,2) * t559 - t423;
t422 = Ifges(7,4) * t559;
t261 = -Ifges(7,2) * t560 + t422;
t263 = -Ifges(7,1) * t560 - t422;
t264 = Ifges(7,1) * t559 - t423;
t677 = qJ(5) * t508;
t774 = -t512 * t758 - t677;
t418 = pkin(3) - t774;
t832 = -t418 * t255 + (t263 / 0.2e1 - t261 / 0.2e1) * t559 - (t264 / 0.2e1 + t260 / 0.2e1) * t560;
t506 = sin(pkin(6));
t510 = sin(qJ(2));
t654 = t506 * t510;
t678 = cos(pkin(6));
t399 = t509 * t678 + t513 * t654;
t514 = cos(qJ(2));
t642 = t512 * t514;
t270 = t399 * t508 + t506 * t642;
t653 = t506 * t514;
t616 = t508 * t653;
t271 = t399 * t512 - t616;
t128 = t270 * t511 - t271 * t507;
t831 = -t128 / 0.2e1;
t830 = t255 / 0.2e1;
t794 = Ifges(5,5) + Ifges(6,4);
t563 = t270 * t507 + t271 * t511;
t26 = t563 * mrSges(7,1) + t128 * mrSges(7,2);
t829 = t26 * qJD(6);
t445 = -pkin(3) * t513 - pkin(9) * t509 - pkin(2);
t629 = pkin(8) * t647 - t512 * t445;
t276 = pkin(10) * t646 - t629;
t501 = t513 * pkin(4);
t226 = pkin(5) * t513 - t276 + t501;
t656 = t445 * t508;
t722 = pkin(8) * t512;
t303 = t656 + (-qJ(5) + t722) * t513;
t481 = pkin(10) * t649;
t247 = t481 + t303;
t110 = t226 * t507 + t247 * t511;
t349 = pkin(8) * t643 + t656;
t277 = t349 + t481;
t133 = -t276 * t507 + t277 * t511;
t828 = -t110 + t133;
t398 = t509 * t654 - t513 * t678;
t193 = t560 * t398;
t194 = t559 * t398;
t640 = t193 * t816 + t194 * t815;
t826 = t640 + (t193 * t442 - t194 * t443) * t766;
t379 = -t507 * t646 + t509 * t648;
t358 = Ifges(7,4) * t379;
t381 = t559 * t509;
t197 = Ifges(7,1) * t381 + t513 * Ifges(7,5) + t358;
t213 = Ifges(7,2) * t381 - t358;
t824 = t197 / 0.4e1 - t213 / 0.4e1;
t823 = -t263 / 0.4e1 + t261 / 0.4e1;
t820 = t264 / 0.4e1 + t260 / 0.4e1;
t336 = -t512 * t654 + t513 * t616;
t337 = (t508 * t510 + t513 * t642) * t506;
t168 = t336 * t511 - t337 * t507;
t169 = t336 * t507 + t337 * t511;
t571 = pkin(4) * t512 + t677;
t444 = -pkin(3) - t571;
t665 = t337 * t512;
t666 = t336 * t508;
t546 = (t665 + t666) * pkin(9);
t617 = t509 * t653;
t795 = mrSges(6,2) + mrSges(5,3);
t819 = -t795 * (t665 / 0.2e1 + t666 / 0.2e1) - m(5) * (-pkin(3) * t617 + t546) / 0.2e1 + (t444 * t617 + t546) * t768 + (t168 * t300 + t169 * t561 - t418 * t617) * t766;
t208 = mrSges(7,1) * t381 + mrSges(7,2) * t379;
t676 = qJ(5) * t512;
t426 = -t508 * t758 + t676;
t545 = -pkin(8) + t426;
t296 = t545 * t509;
t724 = t513 / 0.4e1;
t687 = t379 * mrSges(7,3);
t306 = -mrSges(7,2) * t513 + t687;
t753 = t306 / 0.2e1;
t818 = t632 * t724 + t418 * t208 / 0.2e1 + t296 * t830 + t300 * t753;
t356 = Ifges(7,6) * t381;
t357 = Ifges(7,5) * t379;
t633 = t357 - t356;
t725 = t513 / 0.2e1;
t817 = t296 * t208 + t633 * t725;
t814 = -t509 / 0.2e1;
t813 = -t513 / 0.2e1;
t812 = m(7) * t418;
t809 = t306 * t831;
t109 = t226 * t511 - t247 * t507;
t134 = t276 * t511 + t277 * t507;
t806 = t109 + t134;
t502 = t508 ^ 2;
t504 = t512 ^ 2;
t805 = -t502 - t504;
t577 = t512 * mrSges(6,1) + t508 * mrSges(6,3);
t578 = t512 * mrSges(5,1) - t508 * mrSges(5,2);
t804 = -t578 - t577;
t495 = Ifges(6,5) * t508;
t784 = Ifges(6,1) * t512 + t495;
t369 = -t513 * Ifges(6,4) + t509 * t784;
t714 = Ifges(5,4) * t508;
t460 = Ifges(5,1) * t512 - t714;
t371 = -t513 * Ifges(5,5) + t460 * t509;
t620 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t581 = t620 * t513;
t803 = -t581 + t369 / 0.2e1 + t371 / 0.2e1;
t576 = t507 * mrSges(7,1) + t511 * mrSges(7,2);
t802 = qJD(6) * t576;
t564 = t128 * t507 - t511 * t563;
t562 = t300 * t507 - t511 * t561;
t760 = -mrSges(7,3) / 0.2e1;
t799 = m(6) + m(5);
t797 = mrSges(5,1) + mrSges(6,1);
t796 = mrSges(5,2) - mrSges(6,3);
t793 = Ifges(6,2) + Ifges(5,3);
t685 = t559 * mrSges(7,3);
t727 = t511 / 0.2e1;
t730 = t507 / 0.2e1;
t792 = (t379 * t727 + t381 * t730) * mrSges(7,3);
t746 = t381 / 0.2e1;
t748 = t379 / 0.2e1;
t791 = (t128 * t748 + t563 * t746) * mrSges(7,3);
t359 = Ifges(7,4) * t381;
t195 = Ifges(7,2) * t379 + t513 * Ifges(7,6) + t359;
t215 = -Ifges(7,1) * t379 + t359;
t786 = t215 + t195;
t783 = Ifges(6,6) * t508 + t512 * t794;
t498 = Ifges(5,4) * t512;
t782 = -Ifges(5,2) * t508 + t498;
t459 = Ifges(5,1) * t508 + t498;
t449 = pkin(4) * t508 - t676;
t553 = m(6) * t449;
t451 = t508 * mrSges(5,1) + t512 * mrSges(5,2);
t731 = t451 / 0.2e1;
t450 = t508 * mrSges(6,1) - t512 * mrSges(6,3);
t732 = t450 / 0.2e1;
t781 = t830 + t731 + t732 + t553 / 0.2e1;
t778 = -t352 * t508 + t353 * t512;
t777 = t310 * t512 + t315 * t508;
t776 = Ifges(6,3) * t512 - t495;
t432 = mrSges(5,2) * t513 - mrSges(5,3) * t649;
t493 = t513 * mrSges(6,3);
t625 = mrSges(6,2) * t649;
t439 = -t493 - t625;
t591 = -t439 / 0.2e1 - t432 / 0.2e1;
t555 = t357 / 0.2e1 - t356 / 0.2e1;
t713 = Ifges(6,5) * t512;
t454 = Ifges(6,3) * t508 + t713;
t775 = Ifges(5,6) * t814 + Ifges(6,6) * t728 + t454 * t725 + t782 * t813;
t434 = -mrSges(5,1) * t513 - mrSges(5,3) * t646;
t435 = mrSges(6,1) * t513 + mrSges(6,2) * t646;
t592 = t434 / 0.2e1 - t435 / 0.2e1;
t645 = t511 * t306;
t686 = t381 * mrSges(7,3);
t308 = mrSges(7,1) * t513 - t686;
t652 = t507 * t308;
t637 = -t652 / 0.2e1 + t645 / 0.2e1;
t641 = t168 * t816 - t169 * mrSges(7,2) / 0.2e1;
t757 = t128 / 0.2e1;
t598 = t757 + t831;
t773 = t598 * t685;
t771 = 2 * qJD(3);
t770 = 2 * qJD(4);
t769 = m(5) / 0.2e1;
t767 = m(6) / 0.2e1;
t765 = m(7) / 0.2e1;
t763 = -mrSges(6,1) / 0.2e1;
t759 = -Ifges(6,6) / 0.2e1;
t256 = mrSges(7,1) * t559 - mrSges(7,2) * t560;
t756 = t256 / 0.2e1;
t752 = -t308 / 0.2e1;
t751 = t308 / 0.2e1;
t749 = -t379 / 0.2e1;
t744 = -t398 / 0.2e1;
t743 = t398 / 0.2e1;
t408 = t577 * t509;
t742 = t408 / 0.2e1;
t409 = t578 * t509;
t741 = -t409 / 0.2e1;
t410 = t450 * t509;
t740 = t410 / 0.2e1;
t733 = -t577 / 0.2e1;
t729 = t508 / 0.2e1;
t723 = m(6) * t444;
t719 = mrSges(4,2) * t509;
t705 = t109 * mrSges(7,2);
t704 = t110 * mrSges(7,1);
t697 = t133 * mrSges(7,1);
t696 = t134 * mrSges(7,2);
t684 = t560 * mrSges(7,3);
t682 = t512 * mrSges(6,2);
t681 = t513 * mrSges(4,2);
t680 = t513 * Ifges(5,6);
t679 = -t578 - mrSges(4,1);
t675 = t168 * t560;
t674 = t169 * t559;
t673 = t193 * t511;
t672 = t194 * t507;
t659 = t398 * t512;
t661 = t398 * t508;
t662 = t398 * t399;
t22 = m(7) * (t128 * t193 - t194 * t563 + t662) + t799 * (-t270 * t661 - t271 * t659 + t662);
t671 = t22 * qJD(1);
t660 = t398 * t509;
t23 = 0.4e1 * (t399 * t724 + t660 / 0.4e1 - t654 / 0.4e1) * t653 * m(4) + t799 * (t270 * t336 + t271 * t337 + t398 * t617) + (t128 * t168 + t169 * t563 + t653 * t660) * m(7);
t670 = t23 * qJD(1);
t658 = t442 * t379;
t657 = t443 * t381;
t651 = t507 * t560;
t650 = t507 * t513;
t644 = t511 * t513;
t210 = -mrSges(7,1) * t379 + mrSges(7,2) * t381;
t638 = t210 - t410;
t636 = -t303 + t349;
t304 = t501 + t629;
t635 = t304 - t629;
t634 = t805 * pkin(9) * t398;
t631 = -t577 - t256;
t630 = t644 * t815 + t650 * t816;
t108 = t443 * mrSges(7,1) + mrSges(7,2) * t442;
t626 = qJD(6) * t108;
t622 = t764 + t763;
t621 = mrSges(6,3) / 0.2e1 - mrSges(5,2) / 0.2e1;
t619 = Ifges(5,6) / 0.2e1 + t759;
t618 = t398 * t646;
t615 = -t687 / 0.2e1;
t614 = -t686 / 0.2e1;
t610 = t255 * t744;
t609 = -t661 / 0.2e1;
t606 = -t650 / 0.2e1;
t605 = -t646 / 0.2e1;
t604 = t559 * t727;
t603 = -t644 / 0.2e1;
t601 = t109 / 0.2e1 + t134 / 0.2e1;
t600 = t110 / 0.2e1 - t133 / 0.2e1;
t597 = -t210 / 0.2e1 + t740;
t596 = -t256 / 0.2e1 + t733;
t480 = Ifges(6,5) * t646;
t365 = -Ifges(6,6) * t513 + Ifges(6,3) * t649 + t480;
t367 = t509 * t782 - t680;
t595 = t365 / 0.2e1 - t367 / 0.2e1;
t455 = Ifges(5,2) * t512 + t714;
t590 = -t776 / 0.2e1 - t455 / 0.2e1;
t457 = Ifges(6,1) * t508 - t713;
t589 = t457 / 0.2e1 + t459 / 0.2e1;
t582 = t509 * (mrSges(6,2) / 0.2e1 + mrSges(5,3) / 0.2e1);
t452 = t509 * mrSges(4,1) + t681;
t211 = mrSges(7,1) * t380 + mrSges(7,2) * t382;
t297 = t545 * t513;
t557 = pkin(8) + t449;
t373 = t557 * t509;
t374 = t557 * t513;
t411 = t451 * t509;
t412 = t450 * t513;
t413 = t451 * t513;
t433 = -mrSges(5,2) * t509 - mrSges(5,3) * t647;
t436 = mrSges(5,1) * t509 - mrSges(5,3) * t643;
t517 = (-t436 / 0.2e1 + t736) * t270 + (t438 / 0.2e1 + t433 / 0.2e1) * t271 + (t411 / 0.2e1 + t597) * t399 + (t412 / 0.2e1 + t413 / 0.2e1 - t211 / 0.2e1 + t591 * t512 + t592 * t508) * t398 + (pkin(8) * t399 * t509 - t270 * t352 + t271 * t353 + (pkin(8) * t513 - t349 * t512 - t508 * t629) * t398) * t769 + (t270 * t315 + t271 * t310 + t373 * t399 + (-t303 * t512 - t304 * t508 + t374) * t398) * t767 + (t109 * t193 - t110 * t194 + t112 * t128 + t113 * t563 - t296 * t399 - t297 * t398) * t765 + t309 * t757 + t563 * t307 / 0.2e1 + t193 * t751 - t194 * t753;
t3 = t517 + (-t452 / 0.2e1 + t681 / 0.2e1 + (mrSges(4,1) / 0.2e1 + t578 / 0.2e1 - t596) * t509) * t653 + (-t675 / 0.2e1 + t674 / 0.2e1) * mrSges(7,3) + t819;
t196 = Ifges(7,4) * t382 - Ifges(7,2) * t380 - t509 * Ifges(7,6);
t198 = Ifges(7,1) * t382 - Ifges(7,4) * t380 - t509 * Ifges(7,5);
t370 = t509 * Ifges(6,4) + t513 * t784;
t372 = t509 * Ifges(5,5) + t460 * t513;
t5 = m(6) * (t303 * t310 + t304 * t315 + t373 * t374) + m(7) * (t109 * t112 + t110 * t113 + t296 * t297) - pkin(2) * t452 + (Ifges(4,4) * t513 + pkin(8) * t411 + t803 * t512 + (t513 * t619 + t595) * t508 + t785) * t513 + t315 * t435 + t304 * t437 + t303 * t438 + t310 * t439 + t353 * t432 + t349 * t433 + t352 * t434 + t374 * t410 + t373 * t412 + t113 * t306 + t110 * t307 + t112 * t308 + t109 * t309 + t296 * t211 + t297 * t210 - t629 * t436 + m(5) * (t349 * t353 - t352 * t629) + t197 * t745 + t198 * t746 + t195 * t747 + t196 * t748 + (-Ifges(7,5) * t381 / 0.2e1 + Ifges(7,6) * t749 + pkin(8) * t413 - Ifges(4,4) * t509 + (t372 / 0.2e1 + t370 / 0.2e1 + t620 * t509) * t512 + (-t619 * t509 + t775) * t508 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(7,3) - t793) * t513) * t509;
t570 = t3 * qJD(1) + t5 * qJD(2);
t332 = t774 * t509;
t407 = t571 * t509;
t521 = (t270 * t636 + t271 * t635 + t398 * t407) * t768 + (t128 * t828 - t332 * t398) * t766 - t809 + (t766 * t806 + t752) * t563;
t525 = t622 * t336 + t621 * t337 + (-pkin(4) * t336 + qJ(5) * t337) * t767 + (t168 * t442 + t169 * t443) * t765 - t641;
t6 = -t791 + (-t408 / 0.2e1 + t741 - t208 / 0.2e1) * t398 + (t512 * t582 + t592) * t271 + (t508 * t582 - t591) * t270 + t521 + t525;
t414 = t776 * t509;
t415 = t509 * t455;
t416 = -Ifges(6,1) * t649 + t480;
t417 = t509 * t459;
t479 = Ifges(6,6) * t646;
t8 = t479 * t813 + t373 * t408 + t407 * t410 + t197 * t749 + t213 * t748 + t332 * t210 + t134 * t306 + t133 * t308 + (-t434 + t435) * t349 - (t439 + t432) * t629 + (t109 * t379 + t110 * t381) * mrSges(7,3) + m(6) * (-t303 * t629 + t304 * t349 + t373 * t407) + m(7) * (t109 * t133 + t110 * t134 + t296 * t332) + (pkin(8) * t409 + (t680 / 0.2e1 + t416 / 0.2e1 - t417 / 0.2e1 - t349 * mrSges(5,3) - t303 * mrSges(6,2) + t595) * t512 + (t414 / 0.2e1 + t415 / 0.2e1 - t629 * mrSges(5,3) - t304 * mrSges(6,2) - t803) * t508) * t509 + t786 * t746 - t817;
t569 = -t6 * qJD(1) + t8 * qJD(2);
t537 = t208 * t743 + t563 * t751 + t809;
t13 = t537 + t641 + t791;
t15 = t109 * t306 - t110 * t308 + (-t110 * mrSges(7,3) - t195 / 0.2e1 - t215 / 0.2e1) * t381 + (-t109 * mrSges(7,3) - t213 / 0.2e1 + t197 / 0.2e1) * t379 + t817;
t566 = -t13 * qJD(1) + t15 * qJD(2);
t36 = t638 * t646 + (-t439 - t645 + t652) * t513 + m(7) * (t296 * t646 + (t109 * t507 - t110 * t511) * t513) + m(6) * (-t303 * t513 - t373 * t646);
t526 = (-t271 * t513 - t618) * t768 + (t513 * t564 - t618) * t766;
t535 = t336 * t767 + (t168 * t511 + t169 * t507) * t765;
t41 = t526 + t535;
t565 = qJD(1) * t41 - qJD(2) * t36;
t47 = t630 - t637 + t792;
t558 = t47 * qJD(2) - qJD(4) * t576;
t549 = m(7) * t564;
t547 = m(7) * t562;
t544 = t564 * t765;
t543 = t562 * t765;
t11 = (t621 * t512 + t622 * t508 + t426 * t766 - t553 / 0.2e1 + t781) * t398 - t773 + t826;
t18 = -t449 * t577 - pkin(3) * t451 + (t450 + t553) * t444 + (t782 / 0.2e1 - t454 / 0.2e1 + t589) * t512 + (t460 / 0.2e1 + t784 / 0.2e1 + t590) * t508 + (t812 + t256) * t426 + t832;
t516 = t823 * t379 + t820 * t381 + (t300 * t748 - t560 * t600) * mrSges(7,3) + (t373 * t449 + t407 * t444) * t767 + (t296 * t426 + t300 * t828 + t332 * t418) * t765 + pkin(3) * t741 + t332 * t756 + t373 * t732 + t407 * t733 + t426 * t210 / 0.2e1 + t444 * t742 + t449 * t740 - t786 * t560 / 0.4e1 - t783 * t513 / 0.4e1 - (t601 * mrSges(7,3) - t824) * t559 + (t746 * mrSges(7,3) + t765 * t806 + t751) * t561 - t818;
t522 = pkin(8) * t731 + (t460 / 0.4e1 + t784 / 0.4e1 - t455 / 0.4e1 - t776 / 0.4e1) * t512 + (-t459 / 0.4e1 - t457 / 0.4e1 - t782 / 0.4e1 + t454 / 0.4e1) * t508 + t795 * pkin(9) * (-t504 / 0.2e1 - t502 / 0.2e1);
t527 = (t349 / 0.2e1 - t303 / 0.2e1) * mrSges(6,2) + (t636 * t767 + t591) * pkin(9) + t365 / 0.4e1 - t367 / 0.4e1 + t416 / 0.4e1 - t417 / 0.4e1;
t528 = -t415 / 0.4e1 - t414 / 0.4e1 + t371 / 0.4e1 + t369 / 0.4e1 + (-t629 / 0.2e1 + t304 / 0.2e1) * mrSges(6,2) + (t635 * t767 - t592) * pkin(9);
t2 = ((0.3e1 / 0.4e1 * Ifges(5,6) + t759) * t513 + t527) * t508 + (-t581 + t528) * t512 + t516 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + t522) * t509 + t833;
t542 = t11 * qJD(1) + t2 * qJD(2) + t18 * qJD(3);
t518 = -t824 * t559 - (-t215 / 0.4e1 - t195 / 0.4e1) * t560 + (t300 * t760 - t823) * t379 + (t561 * t760 - t820) * t381 + t561 * t752 + t818;
t10 = t518 - t568;
t20 = t610 - t640;
t541 = t20 * qJD(1) + t10 * qJD(2) - qJD(3) * t832;
t524 = (-t373 * t508 + (-t444 * t509 - t720) * t512) * t768 + (t296 * t508 + t418 * t646 + t513 * t562) * t766;
t530 = m(6) * t750 + (t112 * t511 + t113 * t507) * t765 + t307 * t730 + t309 * t727;
t27 = t597 * t508 + (t682 + (-t604 - t651 / 0.2e1) * mrSges(7,3)) * t513 + (t512 * t596 + t763) * t509 + t524 + t530;
t61 = (t609 - t673 / 0.2e1 + t672 / 0.2e1) * m(7);
t99 = (-t631 - t723 + t812) * t508;
t540 = qJD(1) * t61 - qJD(2) * t27 + qJD(3) * t99;
t539 = t128 * t615 + t563 * t614;
t538 = t507 * t614 + t511 * t615;
t534 = t349 * t768 + (t133 * t511 + t134 * t507) * t766;
t529 = (-t657 / 0.2e1 - t658 / 0.2e1) * mrSges(7,3) + t442 * t753 + t308 * t734 - t555;
t17 = mrSges(7,1) * t600 + mrSges(7,2) * t601 + t529 + t555;
t25 = mrSges(7,2) * t598;
t533 = t25 * qJD(1) + t17 * qJD(2) + t108 * qJD(4);
t185 = m(7) * (-t442 * t507 + t443 * t511) + mrSges(6,3) + m(6) * qJ(5) + t576;
t523 = mrSges(7,1) * t606 + mrSges(7,2) * t603 - t493 + (t656 + (-0.2e1 * qJ(5) + t722) * t513) * t767 + ((-t443 * t513 + t110) * t511 + (t442 * t513 - t109) * t507) * t765 + t637;
t31 = t523 + t534 - t792;
t46 = -t544 + t549 / 0.2e1;
t63 = -t543 + t547 / 0.2e1;
t532 = -t46 * qJD(1) + t31 * qJD(2) - t63 * qJD(3) + t185 * qJD(4);
t505 = t513 ^ 2;
t503 = t509 ^ 2;
t446 = t503 * pkin(8) * t653;
t60 = m(6) * t609 + (-t672 + t673) * t765 + (t768 + t766) * t661;
t55 = -t547 / 0.2e1 - t543 + (m(6) * pkin(9) + mrSges(6,2)) * t512 + (-0.2e1 * t604 - t651) * mrSges(7,3);
t48 = t538 + t630 + t637;
t44 = -t549 / 0.2e1 + m(6) * t271 - t544;
t42 = -t526 + t535;
t30 = t523 - t534 - t538 - t625;
t28 = t210 * t729 + t646 * t756 - t508 * t410 / 0.2e1 - t577 * t605 - t603 * t685 - t606 * t684 - t683 / 0.2e1 - t524 + t530;
t19 = t610 + t640;
t16 = t705 / 0.2e1 + t704 / 0.2e1 - t696 / 0.2e1 + t697 / 0.2e1 + t529 - t555;
t14 = -t537 + t539 + t641;
t12 = t553 * t743 + (t761 + t762) * t659 + t797 * t661 / 0.2e1 + (-t426 * t765 + t781) * t398 - t773 - t826;
t9 = t518 + t568;
t7 = -t208 * t744 + t591 * t270 - t592 * t271 + t398 * t742 + t409 * t743 - t521 + t525 - t539 + t795 * (t271 * t605 - t270 * t649 / 0.2e1);
t4 = t517 + t804 * t617 / 0.2e1 + (-t675 + t674) * t760 - (t452 + t681 + (mrSges(4,1) + t256) * t509) * t653 / 0.2e1 - t819;
t1 = t528 * t512 + (t680 / 0.4e1 + t527) * t508 + t516 + t522 * t509 + t793 * t728 + (-Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t647 + t794 * t643 / 0.2e1 - t833;
t21 = [t23 * qJD(2) + t22 * qJD(3), t4 * qJD(3) + t7 * qJD(4) + t42 * qJD(5) + t14 * qJD(6) + t670 + (t168 * t308 + t169 * t306 - t336 * t434 + t336 * t435 + t337 * t432 + t337 * t439 + ((-mrSges(4,1) * t513 - mrSges(3,1) + t719) * t510 + (-mrSges(3,2) + (t503 + t505) * mrSges(4,3) + (t411 - t638) * t509) * t514) * t506 + 0.2e1 * (t303 * t337 + t304 * t336 + t373 * t617) * t767 + 0.2e1 * (t109 * t168 + t110 * t169 - t296 * t617) * t765 + 0.2e1 * (t336 * t629 + t337 * t349 + t446) * t769 + m(4) * (t446 + (pkin(8) * t505 * t514 - pkin(2) * t510) * t506)) * qJD(2), t671 + t4 * qJD(2) + t12 * qJD(4) + t60 * qJD(5) + t19 * qJD(6) + ((t399 * t444 + t634) * t767 + (t193 * t300 - t194 * t561 - t399 * t418) * t765 + (-pkin(3) * t399 + t634) * t769) * t771 + ((t193 * t560 + t194 * t559) * mrSges(7,3) + (t631 + t679) * t399 + (t795 * t805 + mrSges(4,2)) * t398) * qJD(3), t7 * qJD(2) + t12 * qJD(3) + (t270 * t796 - t271 * t797 - t26) * qJD(4) + t44 * qJD(5) + t829 + ((-pkin(4) * t271 - qJ(5) * t270) * t767 + (-t128 * t443 + t442 * t563) * t765) * t770, qJD(2) * t42 + qJD(3) * t60 + qJD(4) * t44, t14 * qJD(2) + t19 * qJD(3) + t26 * qJD(4) - t829; qJD(3) * t3 - qJD(4) * t6 - qJD(5) * t41 - qJD(6) * t13 - t670, qJD(3) * t5 + qJD(4) * t8 + qJD(5) * t36 + qJD(6) * t15, t1 * qJD(4) + t28 * qJD(5) + t9 * qJD(6) + (t374 * t723 / 0.2e1 + (t112 * t300 + t113 * t561 + t297 * t418) * t765) * t771 + t570 + ((m(5) * t778 + m(6) * t777) * pkin(9) + (Ifges(4,5) + (-m(5) * pkin(3) + t679) * pkin(8)) * t513 - Ifges(4,6) * t509 + (-Ifges(7,5) * t560 - Ifges(7,6) * t559) * t814 + t444 * t412 - t374 * t577 - t560 * t198 / 0.2e1 - t559 * t196 / 0.2e1 + t418 * t211 - pkin(3) * t413 + t561 * t307 + t300 * t309 + t297 * t256 + t112 * t684 + t263 * t745 + t260 * t747 + pkin(8) * t719 - t113 * t685 + (t372 + t370) * t729 + ((-t436 + t437) * pkin(9) + t590 * t513 + t794 * t728) * t508 + ((t433 + t438) * pkin(9) + t589 * t513 + (Ifges(5,6) - Ifges(6,6)) * t728 - t775) * t512 + t778 * mrSges(5,3) + t777 * mrSges(6,2)) * qJD(3), t1 * qJD(3) + t30 * qJD(5) + t16 * qJD(6) + ((t133 * t442 + t134 * t443) * t765 + (-pkin(4) * t349 - qJ(5) * t629) * t767) * t770 + t569 + (t479 + t633 + t696 - t697 + ((-mrSges(6,2) * qJ(5) - Ifges(5,6)) * t512 + (mrSges(6,2) * pkin(4) - t794) * t508) * t509 - t797 * t349 + t796 * t629 + (t657 + t658) * mrSges(7,3)) * qJD(4), qJD(3) * t28 + qJD(4) * t30 + qJD(6) * t48 - t565, t9 * qJD(3) + t16 * qJD(4) + t48 * qJD(5) + (t633 - t704 - t705) * qJD(6) + t566; -qJD(2) * t3 + qJD(4) * t11 + qJD(5) * t61 + qJD(6) * t20 - t671, qJD(4) * t2 - qJD(5) * t27 + qJD(6) * t10 - t570, qJD(4) * t18 + qJD(5) * t99 - qJD(6) * t832, t55 * qJD(5) - t834 + t542 + (m(7) * (-t300 * t443 + t442 * t561) - Ifges(5,6) * t508 - t443 * t684 - t442 * t685 - pkin(4) * t682 - mrSges(6,2) * t677 + (-m(6) * t571 + t804) * pkin(9) + t783 + t825) * qJD(4), qJD(4) * t55 + t540, -qJD(4) * t825 + t541 + t834; qJD(2) * t6 - qJD(3) * t11 - qJD(5) * t46 + qJD(6) * t25, -qJD(3) * t2 + qJD(5) * t31 + qJD(6) * t17 - t569, -qJD(5) * t63 - t542, qJD(5) * t185 + t626, t532, t533 - t626; qJD(2) * t41 - qJD(3) * t61 + qJD(4) * t46, qJD(3) * t27 - qJD(4) * t31 - qJD(6) * t47 + t565, qJD(4) * t63 - t540, -t532 + t802, 0, -t558 - t802; t13 * qJD(2) - t20 * qJD(3) - t25 * qJD(4), -qJD(3) * t10 - qJD(4) * t17 + qJD(5) * t47 - t566, -t541, -qJD(5) * t576 - t533, t558, 0;];
Cq  = t21;
