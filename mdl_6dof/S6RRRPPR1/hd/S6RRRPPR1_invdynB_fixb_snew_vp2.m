% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 04:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:08:01
% EndTime: 2019-05-07 04:08:29
% DurationCPUTime: 27.30s
% Computational Cost: add. (415113->386), mult. (946952->490), div. (0->0), fcn. (706120->12), ass. (0->148)
t766 = -2 * qJD(4);
t738 = sin(qJ(2));
t742 = cos(qJ(2));
t758 = qJD(1) * qJD(2);
t721 = qJDD(1) * t738 + t742 * t758;
t739 = sin(qJ(1));
t743 = cos(qJ(1));
t727 = -g(1) * t743 - g(2) * t739;
t744 = qJD(1) ^ 2;
t716 = -pkin(1) * t744 + qJDD(1) * pkin(7) + t727;
t763 = t738 * t716;
t765 = pkin(2) * t744;
t677 = qJDD(2) * pkin(2) - t721 * pkin(8) - t763 + (pkin(8) * t758 + t738 * t765 - g(3)) * t742;
t703 = -g(3) * t738 + t742 * t716;
t722 = qJDD(1) * t742 - t738 * t758;
t761 = qJD(1) * t738;
t725 = qJD(2) * pkin(2) - pkin(8) * t761;
t732 = t742 ^ 2;
t678 = pkin(8) * t722 - qJD(2) * t725 - t732 * t765 + t703;
t737 = sin(qJ(3));
t741 = cos(qJ(3));
t654 = t741 * t677 - t737 * t678;
t713 = (-t737 * t738 + t741 * t742) * qJD(1);
t687 = qJD(3) * t713 + t721 * t741 + t722 * t737;
t714 = (t737 * t742 + t738 * t741) * qJD(1);
t730 = qJDD(2) + qJDD(3);
t731 = qJD(2) + qJD(3);
t637 = (t713 * t731 - t687) * qJ(4) + (t713 * t714 + t730) * pkin(3) + t654;
t655 = t737 * t677 + t741 * t678;
t686 = -qJD(3) * t714 - t721 * t737 + t722 * t741;
t705 = pkin(3) * t731 - qJ(4) * t714;
t709 = t713 ^ 2;
t641 = -pkin(3) * t709 + qJ(4) * t686 - t705 * t731 + t655;
t734 = sin(pkin(10));
t764 = cos(pkin(10));
t700 = t734 * t713 + t764 * t714;
t624 = t764 * t637 - t734 * t641 + t700 * t766;
t699 = -t764 * t713 + t714 * t734;
t625 = t734 * t637 + t764 * t641 + t699 * t766;
t662 = -t764 * t686 + t687 * t734;
t673 = mrSges(5,1) * t699 + mrSges(5,2) * t700;
t690 = mrSges(5,1) * t731 - mrSges(5,3) * t700;
t672 = pkin(4) * t699 - qJ(5) * t700;
t729 = t731 ^ 2;
t623 = -pkin(4) * t729 + qJ(5) * t730 - t672 * t699 + t625;
t726 = t739 * g(1) - t743 * g(2);
t750 = -qJDD(1) * pkin(1) - t726;
t688 = -t722 * pkin(2) + t725 * t761 + (-pkin(8) * t732 - pkin(7)) * t744 + t750;
t643 = -t686 * pkin(3) - t709 * qJ(4) + t714 * t705 + qJDD(4) + t688;
t663 = t734 * t686 + t764 * t687;
t628 = (t699 * t731 - t663) * qJ(5) + (t700 * t731 + t662) * pkin(4) + t643;
t733 = sin(pkin(11));
t735 = cos(pkin(11));
t684 = t700 * t735 + t731 * t733;
t618 = -0.2e1 * qJD(5) * t684 - t733 * t623 + t735 * t628;
t653 = t663 * t735 + t730 * t733;
t683 = -t700 * t733 + t731 * t735;
t616 = (t683 * t699 - t653) * pkin(9) + (t683 * t684 + t662) * pkin(5) + t618;
t619 = 0.2e1 * qJD(5) * t683 + t735 * t623 + t733 * t628;
t652 = -t663 * t733 + t730 * t735;
t666 = pkin(5) * t699 - pkin(9) * t684;
t682 = t683 ^ 2;
t617 = -pkin(5) * t682 + pkin(9) * t652 - t666 * t699 + t619;
t736 = sin(qJ(6));
t740 = cos(qJ(6));
t614 = t616 * t740 - t617 * t736;
t656 = t683 * t740 - t684 * t736;
t631 = qJD(6) * t656 + t652 * t736 + t653 * t740;
t657 = t683 * t736 + t684 * t740;
t638 = -mrSges(7,1) * t656 + mrSges(7,2) * t657;
t693 = qJD(6) + t699;
t644 = -mrSges(7,2) * t693 + mrSges(7,3) * t656;
t660 = qJDD(6) + t662;
t612 = m(7) * t614 + mrSges(7,1) * t660 - mrSges(7,3) * t631 - t638 * t657 + t644 * t693;
t615 = t616 * t736 + t617 * t740;
t630 = -qJD(6) * t657 + t652 * t740 - t653 * t736;
t645 = mrSges(7,1) * t693 - mrSges(7,3) * t657;
t613 = m(7) * t615 - mrSges(7,2) * t660 + mrSges(7,3) * t630 + t638 * t656 - t645 * t693;
t604 = t740 * t612 + t736 * t613;
t661 = -mrSges(6,1) * t683 + mrSges(6,2) * t684;
t664 = -mrSges(6,2) * t699 + mrSges(6,3) * t683;
t602 = m(6) * t618 + mrSges(6,1) * t662 - mrSges(6,3) * t653 - t661 * t684 + t664 * t699 + t604;
t665 = mrSges(6,1) * t699 - mrSges(6,3) * t684;
t752 = -t612 * t736 + t740 * t613;
t603 = m(6) * t619 - mrSges(6,2) * t662 + mrSges(6,3) * t652 + t661 * t683 - t665 * t699 + t752;
t753 = -t602 * t733 + t735 * t603;
t597 = m(5) * t625 - mrSges(5,2) * t730 - mrSges(5,3) * t662 - t673 * t699 - t690 * t731 + t753;
t689 = -mrSges(5,2) * t731 - mrSges(5,3) * t699;
t622 = -t730 * pkin(4) - t729 * qJ(5) + t700 * t672 + qJDD(5) - t624;
t620 = -t652 * pkin(5) - t682 * pkin(9) + t684 * t666 + t622;
t748 = m(7) * t620 - t630 * mrSges(7,1) + mrSges(7,2) * t631 - t656 * t644 + t645 * t657;
t746 = -m(6) * t622 + t652 * mrSges(6,1) - mrSges(6,2) * t653 + t683 * t664 - t665 * t684 - t748;
t608 = m(5) * t624 + mrSges(5,1) * t730 - mrSges(5,3) * t663 - t673 * t700 + t689 * t731 + t746;
t590 = t734 * t597 + t764 * t608;
t701 = -mrSges(4,1) * t713 + mrSges(4,2) * t714;
t704 = -mrSges(4,2) * t731 + mrSges(4,3) * t713;
t588 = m(4) * t654 + mrSges(4,1) * t730 - mrSges(4,3) * t687 - t701 * t714 + t704 * t731 + t590;
t706 = mrSges(4,1) * t731 - mrSges(4,3) * t714;
t754 = t764 * t597 - t608 * t734;
t589 = m(4) * t655 - mrSges(4,2) * t730 + mrSges(4,3) * t686 + t701 * t713 - t706 * t731 + t754;
t583 = t741 * t588 + t737 * t589;
t702 = -t742 * g(3) - t763;
t720 = (-mrSges(3,1) * t742 + mrSges(3,2) * t738) * qJD(1);
t760 = qJD(1) * t742;
t724 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t760;
t581 = m(3) * t702 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t721 + qJD(2) * t724 - t720 * t761 + t583;
t723 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t761;
t755 = -t588 * t737 + t741 * t589;
t582 = m(3) * t703 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t722 - qJD(2) * t723 + t720 * t760 + t755;
t756 = -t581 * t738 + t742 * t582;
t575 = m(2) * t727 - mrSges(2,1) * t744 - qJDD(1) * mrSges(2,2) + t756;
t715 = -t744 * pkin(7) + t750;
t598 = t735 * t602 + t733 * t603;
t749 = m(5) * t643 + t662 * mrSges(5,1) + t663 * mrSges(5,2) + t699 * t689 + t700 * t690 + t598;
t747 = m(4) * t688 - t686 * mrSges(4,1) + mrSges(4,2) * t687 - t713 * t704 + t706 * t714 + t749;
t745 = -m(3) * t715 + t722 * mrSges(3,1) - mrSges(3,2) * t721 - t723 * t761 + t724 * t760 - t747;
t594 = m(2) * t726 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t744 + t745;
t762 = t739 * t575 + t743 * t594;
t576 = t742 * t581 + t738 * t582;
t757 = t743 * t575 - t594 * t739;
t712 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t738 + Ifges(3,4) * t742) * qJD(1);
t711 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t738 + Ifges(3,2) * t742) * qJD(1);
t710 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t738 + Ifges(3,6) * t742) * qJD(1);
t696 = Ifges(4,1) * t714 + Ifges(4,4) * t713 + Ifges(4,5) * t731;
t695 = Ifges(4,4) * t714 + Ifges(4,2) * t713 + Ifges(4,6) * t731;
t694 = Ifges(4,5) * t714 + Ifges(4,6) * t713 + Ifges(4,3) * t731;
t669 = Ifges(5,1) * t700 - Ifges(5,4) * t699 + Ifges(5,5) * t731;
t668 = Ifges(5,4) * t700 - Ifges(5,2) * t699 + Ifges(5,6) * t731;
t667 = Ifges(5,5) * t700 - Ifges(5,6) * t699 + Ifges(5,3) * t731;
t648 = Ifges(6,1) * t684 + Ifges(6,4) * t683 + Ifges(6,5) * t699;
t647 = Ifges(6,4) * t684 + Ifges(6,2) * t683 + Ifges(6,6) * t699;
t646 = Ifges(6,5) * t684 + Ifges(6,6) * t683 + Ifges(6,3) * t699;
t634 = Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t693;
t633 = Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t693;
t632 = Ifges(7,5) * t657 + Ifges(7,6) * t656 + Ifges(7,3) * t693;
t606 = mrSges(7,2) * t620 - mrSges(7,3) * t614 + Ifges(7,1) * t631 + Ifges(7,4) * t630 + Ifges(7,5) * t660 + t632 * t656 - t633 * t693;
t605 = -mrSges(7,1) * t620 + mrSges(7,3) * t615 + Ifges(7,4) * t631 + Ifges(7,2) * t630 + Ifges(7,6) * t660 - t632 * t657 + t634 * t693;
t592 = mrSges(6,2) * t622 - mrSges(6,3) * t618 + Ifges(6,1) * t653 + Ifges(6,4) * t652 + Ifges(6,5) * t662 - pkin(9) * t604 - t605 * t736 + t606 * t740 + t646 * t683 - t647 * t699;
t591 = -mrSges(6,1) * t622 + mrSges(6,3) * t619 + Ifges(6,4) * t653 + Ifges(6,2) * t652 + Ifges(6,6) * t662 - pkin(5) * t748 + pkin(9) * t752 + t740 * t605 + t736 * t606 - t684 * t646 + t699 * t648;
t584 = Ifges(5,4) * t663 + Ifges(5,6) * t730 - t700 * t667 + t731 * t669 - mrSges(5,1) * t643 + mrSges(5,3) * t625 - Ifges(6,5) * t653 - Ifges(6,6) * t652 - t684 * t647 + t683 * t648 - mrSges(6,1) * t618 + mrSges(6,2) * t619 - Ifges(7,5) * t631 - Ifges(7,6) * t630 - Ifges(7,3) * t660 - t657 * t633 + t656 * t634 - mrSges(7,1) * t614 + mrSges(7,2) * t615 - pkin(5) * t604 - pkin(4) * t598 + (-Ifges(5,2) - Ifges(6,3)) * t662;
t577 = mrSges(5,2) * t643 - mrSges(5,3) * t624 + Ifges(5,1) * t663 - Ifges(5,4) * t662 + Ifges(5,5) * t730 - qJ(5) * t598 - t591 * t733 + t592 * t735 - t667 * t699 - t668 * t731;
t572 = mrSges(4,2) * t688 - mrSges(4,3) * t654 + Ifges(4,1) * t687 + Ifges(4,4) * t686 + Ifges(4,5) * t730 - qJ(4) * t590 + t764 * t577 - t734 * t584 + t713 * t694 - t731 * t695;
t571 = -mrSges(4,1) * t688 + mrSges(4,3) * t655 + Ifges(4,4) * t687 + Ifges(4,2) * t686 + Ifges(4,6) * t730 - pkin(3) * t749 + qJ(4) * t754 + t734 * t577 + t764 * t584 - t714 * t694 + t731 * t696;
t570 = mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) + (-Ifges(4,3) - Ifges(5,3)) * t730 + (-t711 * t738 + t712 * t742) * qJD(1) + t744 * Ifges(2,5) - t735 * t591 - t733 * t592 + mrSges(2,3) * t727 - Ifges(3,5) * t721 - Ifges(3,6) * t722 + t713 * t696 - t714 * t695 - t700 * t668 - mrSges(3,1) * t702 + mrSges(3,2) * t703 - t699 * t669 - Ifges(4,6) * t686 - Ifges(4,5) * t687 - Ifges(5,5) * t663 + Ifges(5,6) * t662 + mrSges(4,2) * t655 - mrSges(4,1) * t654 + mrSges(5,2) * t625 - mrSges(5,1) * t624 - qJ(5) * t753 - pkin(4) * t746 - pkin(1) * t576 - pkin(2) * t583 - pkin(3) * t590 - Ifges(3,3) * qJDD(2);
t569 = mrSges(3,2) * t715 - mrSges(3,3) * t702 + Ifges(3,1) * t721 + Ifges(3,4) * t722 + Ifges(3,5) * qJDD(2) - pkin(8) * t583 - qJD(2) * t711 - t571 * t737 + t572 * t741 + t710 * t760;
t568 = -mrSges(3,1) * t715 + mrSges(3,3) * t703 + Ifges(3,4) * t721 + Ifges(3,2) * t722 + Ifges(3,6) * qJDD(2) - pkin(2) * t747 + pkin(8) * t755 + qJD(2) * t712 + t741 * t571 + t737 * t572 - t710 * t761;
t567 = -mrSges(2,2) * g(3) - mrSges(2,3) * t726 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t744 - pkin(7) * t576 - t568 * t738 + t569 * t742;
t1 = [-m(1) * g(1) + t757; -m(1) * g(2) + t762; (-m(1) - m(2)) * g(3) + t576; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t762 + t743 * t567 - t739 * t570; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t757 + t739 * t567 + t743 * t570; -mrSges(1,1) * g(2) + mrSges(2,1) * t726 + mrSges(1,2) * g(1) - mrSges(2,2) * t727 + Ifges(2,3) * qJDD(1) + pkin(1) * t745 + pkin(7) * t756 + t742 * t568 + t738 * t569;];
tauB  = t1;
