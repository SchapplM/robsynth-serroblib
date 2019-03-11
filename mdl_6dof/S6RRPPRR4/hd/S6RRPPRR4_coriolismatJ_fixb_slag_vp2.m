% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:22
% EndTime: 2019-03-09 09:01:47
% DurationCPUTime: 13.64s
% Computational Cost: add. (29259->841), mult. (73917->1141), div. (0->0), fcn. (80066->10), ass. (0->453)
t793 = mrSges(5,2) - mrSges(4,1);
t479 = sin(qJ(6));
t472 = t479 ^ 2;
t482 = cos(qJ(6));
t474 = t482 ^ 2;
t610 = t472 + t474;
t792 = mrSges(7,3) * t610;
t480 = sin(qJ(5));
t473 = t480 ^ 2;
t483 = cos(qJ(5));
t475 = t483 ^ 2;
t608 = t473 + t475;
t476 = sin(pkin(11));
t477 = sin(pkin(6));
t481 = sin(qJ(2));
t654 = cos(pkin(11));
t718 = cos(qJ(2));
t401 = (t476 * t718 + t481 * t654) * t477;
t742 = -t401 / 0.2e1;
t791 = t608 * t742;
t661 = t482 * mrSges(7,1);
t662 = t479 * mrSges(7,2);
t435 = -t661 + t662;
t656 = -mrSges(6,1) + t435;
t768 = m(7) * pkin(5);
t790 = t656 - t768;
t597 = t477 * t718;
t636 = t477 * t481;
t400 = t476 * t636 - t597 * t654;
t478 = cos(pkin(6));
t342 = t400 * t480 + t478 * t483;
t259 = -t342 * t479 + t401 * t482;
t260 = t342 * t482 + t401 * t479;
t341 = -t400 * t483 + t478 * t480;
t683 = t259 * mrSges(7,3);
t174 = -mrSges(7,2) * t341 + t683;
t681 = t260 * mrSges(7,3);
t175 = mrSges(7,1) * t341 - t681;
t723 = -t482 / 0.2e1;
t729 = -t479 / 0.2e1;
t521 = t174 * t729 + t175 * t723;
t728 = t479 / 0.2e1;
t497 = (t259 * t728 + t260 * t723) * mrSges(7,3) + t521;
t459 = pkin(2) * t636;
t789 = m(4) * t459;
t788 = pkin(8) + qJ(3);
t561 = mrSges(7,3) * (t474 / 0.2e1 + t472 / 0.2e1);
t787 = Ifges(3,6) * t636;
t669 = t342 * mrSges(6,3);
t139 = -mrSges(7,1) * t259 + mrSges(7,2) * t260;
t271 = t401 * mrSges(6,1) - t669;
t786 = -t139 + t271;
t471 = Ifges(7,4) * t482;
t785 = -Ifges(7,2) * t479 + t471;
t441 = Ifges(7,1) * t479 + t471;
t630 = t479 * t480;
t285 = -t400 * t482 - t401 * t630;
t639 = t401 * t483;
t202 = mrSges(7,2) * t639 + t285 * mrSges(7,3);
t626 = t480 * t482;
t286 = -t400 * t479 + t401 * t626;
t203 = -mrSges(7,1) * t639 - t286 * mrSges(7,3);
t722 = t482 / 0.2e1;
t784 = t202 * t722 + t203 * t729;
t710 = t480 * pkin(10);
t714 = pkin(5) * t483;
t448 = t710 + t714;
t596 = t654 * pkin(2);
t464 = -t596 - pkin(3);
t461 = -pkin(9) + t464;
t629 = t479 * t483;
t350 = t482 * t448 - t461 * t629;
t618 = t482 * t483;
t351 = t479 * t448 + t461 * t618;
t539 = -t350 * t479 + t351 * t482;
t717 = pkin(1) * t478;
t460 = t718 * t717;
t366 = -t636 * t788 + t460;
t605 = t481 * t717;
t367 = t597 * t788 + t605;
t573 = t654 * t367;
t267 = t366 * t476 + t573;
t716 = pkin(4) * t400;
t204 = t267 - t716;
t567 = qJ(4) * t400 + t459;
t761 = pkin(3) + pkin(9);
t216 = t401 * t761 + t567;
t105 = t204 * t483 - t480 * t216;
t106 = t480 * t204 + t483 * t216;
t783 = -t105 * t483 - t480 * t106;
t782 = (Ifges(4,4) + Ifges(5,6)) * t401;
t432 = t480 * mrSges(7,1) - mrSges(7,3) * t618;
t619 = t482 * t432;
t582 = -t619 / 0.2e1;
t720 = -t483 / 0.2e1;
t781 = t720 * t792 + t582;
t780 = Ifges(3,5) * t597 - t787 + (Ifges(5,5) - Ifges(4,6)) * t401 + (Ifges(5,4) - Ifges(4,5)) * t400;
t779 = -mrSges(6,2) + t792;
t777 = -(pkin(8) * t597 + t605) * mrSges(3,1) - (-pkin(8) * t636 + t460) * mrSges(3,2);
t345 = pkin(2) * t478 + t366;
t254 = t476 * t345 + t573;
t242 = -t478 * qJ(4) - t254;
t470 = t478 * mrSges(5,3);
t667 = t400 * mrSges(5,1);
t346 = -t470 + t667;
t705 = mrSges(4,3) * t400;
t776 = m(4) * t254 - m(5) * t242 - mrSges(4,2) * t478 - t346 - t705;
t347 = t476 * t367;
t547 = -t345 * t654 + t347;
t245 = -t478 * pkin(3) + t547;
t704 = mrSges(4,3) * t401;
t706 = mrSges(5,1) * t401;
t775 = m(4) * t547 + m(5) * t245 + t478 * t793 + t704 + t706;
t773 = m(5) / 0.2e1;
t772 = m(6) / 0.2e1;
t771 = -m(7) / 0.2e1;
t770 = m(7) / 0.2e1;
t769 = m(4) * pkin(2);
t767 = mrSges(7,1) / 0.2e1;
t766 = -mrSges(7,2) / 0.2e1;
t765 = -mrSges(7,3) / 0.2e1;
t764 = mrSges(7,3) / 0.2e1;
t713 = t401 * pkin(4);
t178 = -t478 * t761 + t547 + t713;
t428 = (-pkin(2) * t718 - pkin(1)) * t477;
t504 = -t401 * qJ(4) + t428;
t201 = t400 * t761 + t504;
t87 = t178 * t483 - t480 * t201;
t78 = -t401 * pkin(5) - t87;
t763 = -t78 / 0.2e1;
t762 = t78 / 0.2e1;
t677 = t286 * mrSges(7,2);
t678 = t285 * mrSges(7,1);
t151 = t677 - t678;
t759 = t151 / 0.2e1;
t758 = t174 / 0.2e1;
t757 = t175 / 0.2e1;
t756 = t259 / 0.2e1;
t755 = t259 / 0.4e1;
t754 = t260 / 0.2e1;
t753 = t260 / 0.4e1;
t752 = t285 / 0.2e1;
t751 = t286 / 0.2e1;
t327 = -t401 * t618 + t479 * t478;
t750 = -t327 / 0.2e1;
t712 = t476 * pkin(2);
t463 = qJ(4) + t712;
t711 = t480 * pkin(5);
t417 = -pkin(10) * t483 + t463 + t711;
t328 = t417 * t482 - t461 * t630;
t749 = t328 / 0.2e1;
t609 = t473 - t475;
t333 = (0.1e1 - t610) * t609;
t748 = -t333 / 0.2e1;
t747 = t341 / 0.2e1;
t746 = t341 / 0.4e1;
t745 = t350 / 0.2e1;
t744 = -t400 / 0.2e1;
t743 = t400 / 0.2e1;
t741 = t401 / 0.2e1;
t660 = t482 * mrSges(7,2);
t663 = t479 * mrSges(7,1);
t557 = t660 + t663;
t420 = t483 * t557;
t740 = t420 / 0.2e1;
t430 = -t480 * mrSges(7,2) - mrSges(7,3) * t629;
t739 = t430 / 0.2e1;
t738 = -t432 / 0.2e1;
t737 = t432 / 0.2e1;
t736 = -t435 / 0.2e1;
t437 = Ifges(7,5) * t479 + Ifges(7,6) * t482;
t735 = t437 / 0.2e1;
t696 = Ifges(7,4) * t479;
t438 = Ifges(7,2) * t482 + t696;
t734 = -t438 / 0.4e1;
t733 = -t441 / 0.4e1;
t732 = t461 / 0.2e1;
t731 = -t478 / 0.2e1;
t730 = t478 / 0.2e1;
t727 = t479 / 0.4e1;
t726 = -t480 / 0.2e1;
t725 = t480 / 0.2e1;
t724 = t480 / 0.4e1;
t721 = t482 / 0.4e1;
t719 = t483 / 0.2e1;
t418 = -mrSges(7,1) * t618 + mrSges(7,2) * t629;
t715 = pkin(5) * t418;
t181 = -t242 - t716;
t110 = pkin(5) * t341 - pkin(10) * t342 + t181;
t88 = t480 * t178 + t201 * t483;
t79 = t401 * pkin(10) + t88;
t54 = t110 * t482 - t479 * t79;
t709 = t54 * mrSges(7,3);
t55 = t110 * t479 + t482 * t79;
t708 = t55 * mrSges(7,3);
t707 = m(7) * qJD(5);
t703 = mrSges(7,3) * t341;
t701 = Ifges(3,4) * t481;
t699 = Ifges(6,4) * t342;
t698 = Ifges(6,4) * t480;
t697 = Ifges(6,4) * t483;
t695 = Ifges(7,5) * t286;
t694 = Ifges(7,5) * t480;
t693 = Ifges(7,5) * t482;
t690 = Ifges(7,6) * t285;
t689 = Ifges(7,6) * t479;
t688 = Ifges(7,6) * t480;
t687 = Ifges(6,3) * t400;
t686 = Ifges(7,3) * t342;
t685 = Ifges(7,3) * t483;
t107 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t341;
t671 = t341 * Ifges(7,6);
t679 = t260 * Ifges(7,4);
t682 = t259 * Ifges(7,2);
t108 = t671 + t679 + t682;
t252 = Ifges(7,4) * t259;
t672 = t341 * Ifges(7,5);
t680 = t260 * Ifges(7,1);
t109 = t252 + t672 + t680;
t120 = -Ifges(7,3) * t639 + t690 + t695;
t121 = Ifges(7,4) * t286 + Ifges(7,2) * t285 - Ifges(7,6) * t639;
t122 = Ifges(7,1) * t286 + Ifges(7,4) * t285 - Ifges(7,5) * t639;
t664 = t401 * Ifges(6,6);
t179 = -Ifges(6,2) * t341 + t664 + t699;
t332 = Ifges(6,4) * t341;
t665 = t401 * Ifges(6,5);
t180 = Ifges(6,1) * t342 - t332 + t665;
t268 = t366 * t654 - t347;
t205 = t268 - t713;
t553 = Ifges(6,2) * t483 + t698;
t217 = -Ifges(6,6) * t400 + t401 * t553;
t555 = Ifges(6,1) * t480 + t697;
t218 = -Ifges(6,5) * t400 + t401 * t555;
t670 = t342 * mrSges(6,2);
t674 = t341 * mrSges(6,1);
t236 = t670 + t674;
t262 = t400 * pkin(3) + t504;
t666 = t401 * mrSges(6,2);
t673 = t341 * mrSges(6,3);
t270 = -t666 - t673;
t284 = pkin(3) * t401 + t567;
t436 = t483 * mrSges(6,1) - t480 * mrSges(6,2);
t287 = t436 * t401;
t640 = t401 * t480;
t295 = -mrSges(6,1) * t400 - mrSges(6,3) * t640;
t296 = t400 * mrSges(6,2) + mrSges(6,3) * t639;
t297 = -mrSges(5,2) * t400 - mrSges(5,3) * t401;
t563 = t597 / 0.2e1;
t571 = -t401 * mrSges(5,2) + t400 * mrSges(5,3);
t572 = t401 * mrSges(4,1) - t400 * mrSges(4,2);
t591 = t639 / 0.2e1;
t592 = -t639 / 0.2e1;
t593 = t640 / 0.2e1;
t134 = (-pkin(4) - t448) * t401 + t268;
t85 = -pkin(10) * t400 + t106;
t62 = t134 * t482 - t479 * t85;
t63 = t134 * t479 + t482 * t85;
t84 = t400 * pkin(5) - t105;
t1 = (t572 + t789) * t428 + t179 * t591 + t107 * t592 + t180 * t593 - t245 * t667 + (mrSges(4,1) * t400 + mrSges(4,2) * t401) * t459 + (0.2e1 * Ifges(5,6) * t400 + (-Ifges(5,2) + Ifges(5,3)) * t401) * t743 + (Ifges(6,5) * t342 - Ifges(6,6) * t341 - 0.2e1 * Ifges(4,4) * t400 + (Ifges(4,1) - Ifges(4,2) + Ifges(6,3)) * t401) * t744 + (0.2e1 * Ifges(3,4) * t597 + (Ifges(3,1) - Ifges(3,2)) * t636) * t563 + (-t787 / 0.2e1 + Ifges(3,5) * t563 + Ifges(5,5) * t741 + Ifges(4,6) * t742 + Ifges(5,4) * t743 + Ifges(4,5) * t744 + t777) * t478 + (-t687 + (Ifges(6,5) * t480 + Ifges(6,6) * t483) * t401 + (Ifges(5,3) - Ifges(4,1)) * t400 - t782) * t741 + ((-Ifges(4,2) + Ifges(5,2)) * t400 + t782) * t742 + t780 * t730 + t776 * t268 + t775 * t267 + (m(5) * t284 + t571) * t262 - t547 * t705 - t341 * t217 / 0.2e1 + t342 * t218 / 0.2e1 + t284 * t297 - t181 * t287 + t87 * t295 + t88 * t296 + t106 * t270 + t105 * t271 + t205 * t236 + t55 * t202 + t54 * t203 + t63 * t174 + t62 * t175 + t78 * t151 + t84 * t139 + m(6) * (t105 * t87 + t106 * t88 + t181 * t205) + m(7) * (t54 * t62 + t55 * t63 + t78 * t84) - t254 * t704 + (-(Ifges(3,2) * t718 + t701) * t636 / 0.2e1 + (-pkin(1) * (mrSges(3,1) * t481 + mrSges(3,2) * t718) + t481 * (Ifges(3,1) * t718 - t701) / 0.2e1) * t477) * t477 + t242 * t706 + t120 * t747 + t109 * t751 + t108 * t752 + t122 * t754 + t121 * t756;
t684 = t1 * qJD(1);
t676 = t328 * mrSges(7,3);
t329 = t417 * t479 + t461 * t626;
t675 = t329 * mrSges(7,3);
t551 = -t689 + t693;
t152 = -t341 * t551 + t686;
t153 = Ifges(7,6) * t342 - t341 * t785;
t554 = Ifges(7,1) * t482 - t696;
t154 = Ifges(7,5) * t342 - t341 * t554;
t219 = t557 * t341;
t226 = -mrSges(7,2) * t342 + t479 * t703;
t227 = mrSges(7,1) * t342 + t482 * t703;
t235 = mrSges(6,1) * t342 - mrSges(6,2) * t341;
t237 = -Ifges(6,2) * t342 - t332;
t238 = -Ifges(6,1) * t341 - t699;
t613 = -Ifges(6,5) * t341 - Ifges(6,6) * t342;
t624 = t482 * t109;
t635 = t479 * t108;
t239 = pkin(5) * t342 + pkin(10) * t341;
t72 = t239 * t482 - t479 * t87;
t73 = t239 * t479 + t482 * t87;
t4 = t87 * t270 - t78 * t219 + t73 * t174 + t55 * t226 + t72 * t175 + t54 * t227 + t154 * t754 + t153 * t756 + m(7) * (t54 * t72 + t55 * t73) + t181 * t235 + t613 * t741 + (t238 / 0.2e1 - t179 / 0.2e1 + t107 / 0.2e1) * t342 + (-t180 / 0.2e1 - t237 / 0.2e1 + t152 / 0.2e1 - t624 / 0.2e1 + t635 / 0.2e1 + t87 * mrSges(6,3)) * t341 + (m(7) * t78 - t669 - t786) * t88;
t668 = t4 * qJD(1);
t659 = t62 * t479;
t658 = t63 * t482;
t138 = mrSges(7,1) * t260 + mrSges(7,2) * t259;
t140 = Ifges(7,5) * t259 - Ifges(7,6) * t260;
t141 = -Ifges(7,2) * t260 + t252;
t142 = Ifges(7,1) * t259 - t679;
t7 = t78 * t138 + t140 * t747 + t54 * t174 - t55 * t175 + (-t708 + t142 / 0.2e1 - t108 / 0.2e1) * t260 + (-t709 + t109 / 0.2e1 + t141 / 0.2e1) * t259;
t657 = t7 * qJD(1);
t632 = t479 * t430;
t588 = -t632 / 0.2e1;
t519 = t588 + t582;
t101 = t418 * t726 - t475 * t561 + t483 * t519;
t429 = -mrSges(7,2) * t483 + mrSges(7,3) * t630;
t620 = t482 * t429;
t431 = mrSges(7,1) * t483 + mrSges(7,3) * t626;
t631 = t479 * t431;
t506 = t740 - t631 / 0.2e1 + t620 / 0.2e1;
t419 = t480 * t557;
t507 = t419 / 0.2e1 + t432 * t729 + t430 * t722;
t517 = m(7) * t539;
t540 = t328 * t479 - t329 * t482;
t65 = m(7) * t609 * t732 + (t517 / 0.2e1 + t506) * t483 + (t540 * t770 - t507) * t480;
t655 = t65 * qJD(5) + t101 * qJD(6);
t368 = t475 * t401;
t558 = t480 * mrSges(6,1) + t483 * mrSges(6,2);
t641 = t401 * t461;
t642 = t400 * t463;
t484 = (t401 * t464 - t642) * t773 + (t608 * t641 - t642) * t772 + (t285 * t328 + t286 * t329 + t368 * t461) * t770 + t285 * t737 + t286 * t739 + t558 * t744 + (-t400 * t476 - t401 * t654) * t769 / 0.2e1 + t420 * t592 + mrSges(6,3) * t791;
t550 = t658 - t659;
t490 = t284 * t773 + (-t480 * t105 + t106 * t483) * t772 + (t480 * t84 + t483 * t550) * t770 + t789 / 0.2e1;
t581 = t618 / 0.2e1;
t586 = -t629 / 0.2e1;
t18 = t151 * t725 + t202 * t581 + t203 * t586 + t295 * t726 + t296 * t719 - t484 + t490 + t571 + t572;
t653 = qJD(1) * t18;
t645 = t327 * t482;
t326 = t401 * t629 + t482 * t478;
t646 = t326 * t479;
t541 = t645 - t646;
t637 = t473 * t401;
t499 = m(7) * (t483 * t541 - t637);
t647 = t286 * t482;
t648 = t285 * t479;
t542 = t647 - t648;
t502 = (t480 * t542 + t368) * t770;
t41 = t502 + m(5) * t401 - t499 / 0.2e1 + 0.2e1 * (t368 / 0.4e1 + (t473 / 0.2e1 + t475 / 0.4e1) * t401) * m(6);
t652 = qJD(1) * t41;
t615 = t483 * t271;
t617 = t483 * t139;
t627 = t480 * t270;
t16 = -t285 * t175 - t286 * t174 - m(7) * (t54 * t285 + t55 * t286 - t639 * t78) + (m(6) * t181 + t236 + t776) * t400 + (t617 - t627 - t615 - m(6) * (t480 * t88 + t483 * t87) - t775) * t401;
t650 = t16 * qJD(1);
t616 = t483 * t270;
t17 = t327 * t174 + t326 * t175 + (t236 - t346) * t478 + (t786 * t480 - t297 - t616) * t401 + m(7) * (t326 * t54 + t327 * t55 - t640 * t78) + m(6) * (t181 * t478 + (t480 * t87 - t483 * t88) * t401) + m(5) * (-t242 * t478 - t262 * t401);
t649 = t17 * qJD(1);
t638 = t461 * t483;
t634 = t479 * t227;
t405 = t483 * t785 + t688;
t633 = t479 * t405;
t625 = t480 * t483;
t623 = t482 * t174;
t622 = t482 * t226;
t516 = t554 * t483;
t407 = t516 + t694;
t621 = t482 * t407;
t614 = t483 * t418;
t607 = qJD(5) * t480;
t606 = qJD(5) * t483;
t604 = t88 * t771;
t568 = t610 * t483;
t363 = t480 * t568 - t625;
t603 = t363 * t707;
t599 = -t669 / 0.2e1;
t598 = mrSges(6,3) * t719;
t594 = -t640 / 0.2e1;
t587 = -t630 / 0.2e1;
t580 = -t108 / 0.4e1 + t142 / 0.4e1;
t579 = -t109 / 0.4e1 - t141 / 0.4e1;
t578 = t759 - t295 / 0.2e1;
t577 = t219 / 0.2e1 + t270 / 0.2e1;
t423 = t441 * t483;
t576 = -t405 / 0.4e1 - t423 / 0.4e1;
t422 = t438 * t483;
t575 = t407 / 0.4e1 - t422 / 0.4e1;
t569 = qJD(5) * t656;
t549 = -t479 * t72 + t482 * t73;
t548 = t479 * t785;
t524 = t139 / 0.2e1 - t271 / 0.2e1 + t599;
t534 = t549 + t78;
t489 = t534 * t771 - t622 / 0.2e1 + t634 / 0.2e1 + (t768 / 0.2e1 + t736 + mrSges(6,1) / 0.2e1) * t401 - t524;
t523 = t673 / 0.2e1 + t577;
t496 = -t666 / 0.2e1 + t175 * t729 + t623 / 0.2e1 + t523;
t512 = (-t648 / 0.2e1 + t647 / 0.2e1) * mrSges(7,3);
t515 = t542 * pkin(10);
t535 = t479 * t54 - t482 * t55 + t88;
t12 = t515 * t770 + t512 + t489 * t483 + (t535 * t771 + t496) * t480;
t546 = t12 * qJD(1) - t65 * qJD(2);
t511 = (t645 / 0.2e1 - t646 / 0.2e1) * mrSges(7,3);
t514 = t541 * pkin(10);
t14 = t514 * t770 + t511 + t489 * t480 + (t535 * t770 - t496) * t483;
t67 = (t540 * t771 + t507) * t483 + ((t539 - 0.2e1 * t638) * t770 + t506) * t480;
t545 = t14 * qJD(1) - t67 * qJD(2);
t494 = t138 * t720 + t480 * t497;
t530 = mrSges(7,2) * t750 + t326 * t767;
t26 = t494 - t530;
t526 = -t662 / 0.2e1 + t661 / 0.2e1;
t89 = -t614 / 0.2e1 + (t483 * t561 - t519) * t480 + t526;
t544 = -t26 * qJD(1) + t89 * qJD(2);
t493 = t138 * t725 + t483 * t497;
t531 = -t678 / 0.2e1 + t677 / 0.2e1;
t25 = t493 + t531;
t538 = -qJD(1) * t25 - qJD(2) * t101;
t536 = t463 * t478 - t242;
t486 = -t470 - m(5) * t536 / 0.2e1 - m(6) * (t536 - t716) / 0.2e1 + (t328 * t326 + t329 * t327 + t479 * t55 + t482 * t54 + t625 * t641) * t771 + t326 * t738 + t430 * t750 - t674 / 0.2e1 - t670 / 0.2e1 + t521;
t491 = t267 * t773 - t783 * t772 + (t480 * t550 - t483 * t84) * t770;
t508 = t296 / 0.2e1 + t784;
t11 = t486 + (mrSges(6,2) * t731 - t578) * t483 + (mrSges(6,1) * t731 + t401 * t740 + t508) * t480 + t491;
t113 = t632 + t619 + mrSges(5,3) + m(7) * (t328 * t482 + t329 * t479) + 0.4e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * t463 + t558;
t537 = qJD(1) * t11 - qJD(2) * t113;
t533 = -t693 / 0.2e1 + t689 / 0.2e1;
t532 = pkin(5) * t138 / 0.2e1 + t438 * t753;
t529 = qJD(3) * t748 - t363 * qJD(4);
t528 = mrSges(7,1) * t745 + t351 * t766;
t527 = t363 * qJD(3) + qJD(4) * t748;
t525 = t660 / 0.2e1 + t663 / 0.2e1;
t522 = t635 / 0.4e1 - t624 / 0.4e1;
t520 = t633 / 0.4e1 - t621 / 0.4e1;
t518 = t438 * t728 + t441 * t723;
t421 = t483 * t437;
t513 = -t420 / 0.2e1;
t40 = -t421 * t725 + t328 * t430 - t329 * t432 + (t461 * t418 + (-t675 - t423 / 0.2e1 - t405 / 0.2e1) * t482 + (t676 - t407 / 0.2e1 + t422 / 0.2e1) * t479) * t483;
t487 = (-t676 / 0.2e1 + t575) * t259 + (-t675 / 0.2e1 + t576) * t260 + t174 * t749 - t329 * t175 / 0.2e1 - t421 * t746 + t140 * t724 + t54 * t739 + t55 * t738 + t418 * t763;
t495 = -t461 * t138 / 0.2e1 + (-t708 / 0.2e1 + t580) * t482 + (t709 / 0.2e1 + t579) * t479;
t498 = -t695 / 0.2e1 - t690 / 0.2e1 - t62 * mrSges(7,1) / 0.2e1 + t63 * mrSges(7,2) / 0.2e1;
t5 = (Ifges(7,3) * t741 + t495) * t483 + t487 + t498;
t510 = t5 * qJD(1) + t40 * qJD(2) + t101 * qJD(3);
t228 = -pkin(5) * t557 + t554 * t728 + t722 * t785 - t518;
t38 = -t715 / 0.2e1 + (-0.3e1 / 0.4e1 * t694 + pkin(10) * t737 - t575) * t482 + (0.3e1 / 0.4e1 * t688 + pkin(10) * t739 - t576) * t479 + (Ifges(7,3) / 0.2e1 + t557 * t732 + t441 * t727 - t482 * t554 / 0.4e1 + t438 * t721 + t548 / 0.4e1 + pkin(10) * t561) * t483 + t528;
t505 = t686 / 0.2e1 + t72 * t767 + t73 * t766;
t8 = (t733 - t471 / 0.4e1) * t259 + (0.3e1 / 0.4e1 * t671 + t679 / 0.4e1 + t682 / 0.4e1 + mrSges(7,1) * t763 + (t758 - t683 / 0.2e1) * pkin(10) - t580) * t479 + (-0.3e1 / 0.4e1 * t672 - t680 / 0.4e1 + mrSges(7,2) * t763 + (t681 / 0.2e1 + t757) * pkin(10) + t579) * t482 + t505 + t532;
t509 = t8 * qJD(1) + t38 * qJD(2) - t228 * qJD(5);
t403 = Ifges(7,3) * t480 + t483 * t551;
t404 = Ifges(7,6) * t483 - t480 * t785;
t406 = Ifges(7,5) * t483 - t480 * t554;
t440 = -Ifges(6,2) * t480 + t697;
t485 = (-t440 / 0.4e1 - t555 / 0.4e1 + t403 / 0.4e1) * t342 + (t328 * t72 + t329 * t73 + t350 * t54 + t351 * t55) * t770 + t181 * t436 / 0.2e1 + t404 * t755 + t406 * t753 + t227 * t749 + t329 * t226 / 0.2e1 + t175 * t745 + t351 * t758 + t463 * t235 / 0.2e1 + t54 * t431 / 0.2e1 + t55 * t429 / 0.2e1 + t72 * t737 + t73 * t739 - t419 * t762 + t88 * t740;
t488 = (-pkin(5) * t84 + pkin(10) * t550) * t771 + t687 / 0.2e1 + pkin(5) * t759 - t105 * mrSges(6,1) / 0.2e1 + t106 * mrSges(6,2) / 0.2e1 + t285 * t734 + t286 * t733 + t84 * t736;
t492 = (m(7) * t762 + t524) * t461 + t152 / 0.4e1 - t180 / 0.4e1 - t237 / 0.4e1 + t522;
t402 = -t480 * t551 + t685;
t443 = Ifges(6,1) * t483 - t698;
t500 = t553 / 0.4e1 - t443 / 0.4e1 + t402 / 0.4e1 + t520;
t501 = t107 / 0.4e1 - t179 / 0.4e1 + t238 / 0.4e1 - t479 * t153 / 0.4e1 + t154 * t721;
t2 = t485 + (t63 * t765 - pkin(10) * t202 / 0.2e1 - t121 / 0.4e1) * t482 + (t62 * t764 + pkin(10) * t203 / 0.2e1 - t122 / 0.4e1) * t479 + ((-0.3e1 / 0.4e1 * Ifges(6,6) + t437 / 0.4e1) * t401 + (t604 + t523) * t461 + t501) * t483 + (-0.3e1 / 0.4e1 * t665 + t492) * t480 + t500 * t341 + t488;
t36 = t351 * t430 + t329 * t429 + m(7) * (t328 * t350 + t329 * t351) + t350 * t432 + t328 * t431 + t463 * t436 + (t403 / 0.2e1 - t440 / 0.2e1 - t555 / 0.2e1 + t461 * t419 + t406 * t722 + t404 * t729) * t483 + (t402 / 0.2e1 + t553 / 0.2e1 - t443 / 0.2e1 - t621 / 0.2e1 + t633 / 0.2e1 + (-m(7) * t638 + t420) * t461) * t480;
t503 = t2 * qJD(1) + t36 * qJD(2) + t65 * qJD(3) + t67 * qJD(4);
t340 = -t483 * t525 + t513;
t339 = t480 * t525 + t557 * t725;
t311 = t333 * t707 / 0.2e1;
t90 = t614 / 0.2e1 + t430 * t587 + t526 + t781 * t480;
t66 = t67 * qJD(5);
t42 = t499 / 0.2e1 + m(6) * t791 + t502 + (t368 + t637) * t772;
t39 = t715 / 0.2e1 + t551 * t724 - t423 * t727 + t461 * t513 + t629 * t733 + t618 * t734 - t483 * t548 / 0.4e1 + t685 / 0.2e1 + t533 * t480 - t520 + t528 + (-t422 + t516) * t721 + (t588 + t781) * pkin(10);
t27 = t494 + t530;
t24 = t493 - t531;
t19 = t480 * t578 + t483 * t508 + t484 + t490;
t15 = -t219 * t720 + t174 * t581 + t175 * t586 + t227 * t587 + t271 * t726 + t341 * t598 + t616 / 0.2e1 + t480 * t599 + t435 * t594 + mrSges(6,1) * t593 + mrSges(6,2) * t591 + t511 + (pkin(5) * t640 + t480 * t534 - t483 * t535 + t514) * t770 + (t139 + t622) * t725;
t13 = t617 / 0.2e1 - t219 * t725 + t226 * t581 + t630 * t757 + t227 * t586 - t615 / 0.2e1 - t627 / 0.2e1 + t483 * t599 + mrSges(6,1) * t591 + t435 * t592 + mrSges(6,2) * t594 + t512 + (pkin(5) * t639 + t480 * t535 + t483 * t534 + t515) * t770 + (t623 + t673) * t726;
t10 = t420 * t594 + t480 * t508 - t483 * t578 + t558 * t730 - t486 + t491 - t667;
t9 = t141 * t721 + t142 * t727 + t341 * t533 + t551 * t746 + t554 * t753 + t557 * t762 + t505 - t522 - t532 + (t441 + t785) * t755 + t497 * pkin(10);
t6 = Ifges(7,3) * t592 + t483 * t495 + t487 - t498;
t3 = -t488 + t485 + (-t665 / 0.4e1 + t492) * t480 + (-t664 / 0.4e1 + (t604 + t577) * t461 + t501) * t483 + Ifges(6,6) * t591 - t437 * t639 / 0.4e1 + t658 * t764 + Ifges(6,5) * t593 + t659 * t765 + t121 * t721 + t122 * t727 + (t461 * t598 + t500) * t341 + t784 * pkin(10);
t20 = [qJD(2) * t1 - qJD(3) * t16 + qJD(4) * t17 + qJD(5) * t4 + qJD(6) * t7, t684 + (t440 * t591 + t403 * t592 + t443 * t593 + t122 * t581 + t121 * t586 - t464 * t667 - t704 * t712 + (m(6) * t205 - t287 - t706) * t463 + t783 * mrSges(6,3) + t62 * t432 + t63 * t430 + t84 * t420 + t329 * t202 + t328 * t203 + (-m(6) * t783 + t480 * t296) * t461 + (m(5) * t464 - t654 * t769 + t793) * t267 + (-m(7) * t84 - t151 + t295) * t638 + (m(5) * t463 + t476 * t769 - mrSges(4,2) + mrSges(5,3)) * t268 + t205 * t558 + t780 + t596 * t705 + t777 + t218 * t719 + t120 * t725 + t217 * t726 + (Ifges(6,5) * t483 - Ifges(6,6) * t480) * t744 + t407 * t751 + t405 * t752 + m(7) * (t328 * t62 + t329 * t63)) * qJD(2) + t19 * qJD(3) + t10 * qJD(4) + t3 * qJD(5) + t6 * qJD(6), -t650 + t19 * qJD(2) + t42 * qJD(4) + t13 * qJD(5) + t24 * qJD(6) + m(7) * (t542 - t640) * qJD(3) * t483, t649 + t10 * qJD(2) + t42 * qJD(3) + t15 * qJD(5) + t27 * qJD(6) + m(7) * (t541 + t639) * qJD(4) * t480, t668 + t3 * qJD(2) + t13 * qJD(3) + t15 * qJD(4) + (-t87 * mrSges(6,2) + t549 * mrSges(7,3) + pkin(5) * t219 + t153 * t722 + t154 * t728 + t518 * t341 + t342 * t735 + t613 + t790 * t88 + (m(7) * t549 + t622 - t634) * pkin(10)) * qJD(5) + t9 * qJD(6), t657 + t6 * qJD(2) + t24 * qJD(3) + t27 * qJD(4) + t9 * qJD(5) + (-mrSges(7,1) * t55 - mrSges(7,2) * t54 + t140) * qJD(6); -qJD(3) * t18 - qJD(4) * t11 + qJD(5) * t2 + qJD(6) * t5 - t684, qJD(4) * t113 + qJD(5) * t36 + qJD(6) * t40, -t653 + t655, qJD(6) * t90 - t537 + t66, t39 * qJD(6) + (t790 * t461 - Ifges(6,5) + t518) * t607 + (-t461 * mrSges(6,2) - Ifges(6,6) + t735) * t606 + t503 + (pkin(5) * t419 + t404 * t722 + t406 * t728 + (t517 + t620 - t631) * pkin(10) + t539 * mrSges(7,3)) * qJD(5), t90 * qJD(4) + t39 * qJD(5) + (-t329 * mrSges(7,1) - t328 * mrSges(7,2) - t421) * qJD(6) + t510; qJD(2) * t18 - qJD(4) * t41 - qJD(5) * t12 + qJD(6) * t25 + t650, t653 + t655, -t603, t311 - t652, t339 * qJD(6) + t483 * t569 - t779 * t607 + ((-t610 * t710 - t714) * qJD(5) - t527) * m(7) - t546, qJD(5) * t339 + qJD(6) * t418 - t538; qJD(2) * t11 + qJD(3) * t41 - qJD(5) * t14 + qJD(6) * t26 - t649, -qJD(6) * t89 + t537 + t66, t311 + t652, t603, t340 * qJD(6) + t480 * t569 + t779 * t606 + ((pkin(10) * t568 - t711) * qJD(5) - t529) * m(7) - t545, qJD(6) * t435 * t480 + t340 * qJD(5) - t544; -qJD(2) * t2 + qJD(3) * t12 + qJD(4) * t14 - qJD(6) * t8 - t668, -qJD(6) * t38 - t503, m(7) * t527 + t546, m(7) * t529 + t545, t228 * qJD(6) (pkin(10) * t435 + t551) * qJD(6) - t509; -qJD(2) * t5 - qJD(3) * t25 - qJD(4) * t26 + qJD(5) * t8 - t657, qJD(4) * t89 + qJD(5) * t38 - t510, t538, t544, t509, 0;];
Cq  = t20;
