% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:08:56
% EndTime: 2019-05-07 20:09:38
% DurationCPUTime: 30.00s
% Computational Cost: add. (507445->387), mult. (1025000->489), div. (0->0), fcn. (757576->12), ass. (0->149)
t744 = qJD(1) ^ 2;
t762 = pkin(2) * t744;
t738 = sin(qJ(1));
t743 = cos(qJ(1));
t725 = -g(1) * t743 - g(2) * t738;
t713 = -pkin(1) * t744 + qJDD(1) * pkin(7) + t725;
t737 = sin(qJ(2));
t761 = t737 * t713;
t742 = cos(qJ(2));
t757 = qJD(1) * qJD(2);
t719 = qJDD(1) * t737 + t742 * t757;
t678 = qJDD(2) * pkin(2) - t719 * pkin(8) - t761 + (pkin(8) * t757 + t737 * t762 - g(3)) * t742;
t700 = -g(3) * t737 + t742 * t713;
t720 = qJDD(1) * t742 - t737 * t757;
t759 = qJD(1) * t737;
t723 = qJD(2) * pkin(2) - pkin(8) * t759;
t731 = t742 ^ 2;
t679 = pkin(8) * t720 - qJD(2) * t723 - t731 * t762 + t700;
t736 = sin(qJ(3));
t741 = cos(qJ(3));
t657 = t736 * t678 + t741 * t679;
t711 = (t736 * t742 + t737 * t741) * qJD(1);
t684 = -t711 * qJD(3) - t736 * t719 + t720 * t741;
t758 = qJD(1) * t742;
t710 = -t736 * t759 + t741 * t758;
t694 = -mrSges(4,1) * t710 + mrSges(4,2) * t711;
t730 = qJD(2) + qJD(3);
t702 = mrSges(4,1) * t730 - mrSges(4,3) * t711;
t729 = qJDD(2) + qJDD(3);
t685 = qJD(3) * t710 + t719 * t741 + t720 * t736;
t724 = t738 * g(1) - t743 * g(2);
t749 = -qJDD(1) * pkin(1) - t724;
t686 = -t720 * pkin(2) + t723 * t759 + (-pkin(8) * t731 - pkin(7)) * t744 + t749;
t640 = (-t710 * t730 - t685) * pkin(9) + (t711 * t730 - t684) * pkin(3) + t686;
t695 = -pkin(3) * t710 - pkin(9) * t711;
t728 = t730 ^ 2;
t646 = -pkin(3) * t728 + pkin(9) * t729 + t695 * t710 + t657;
t735 = sin(qJ(4));
t740 = cos(qJ(4));
t629 = t740 * t640 - t735 * t646;
t697 = -t711 * t735 + t730 * t740;
t660 = qJD(4) * t697 + t685 * t740 + t729 * t735;
t683 = qJDD(4) - t684;
t698 = t711 * t740 + t730 * t735;
t706 = qJD(4) - t710;
t622 = (t697 * t706 - t660) * qJ(5) + (t697 * t698 + t683) * pkin(4) + t629;
t630 = t735 * t640 + t740 * t646;
t659 = -qJD(4) * t698 - t685 * t735 + t729 * t740;
t688 = pkin(4) * t706 - qJ(5) * t698;
t696 = t697 ^ 2;
t624 = -pkin(4) * t696 + qJ(5) * t659 - t688 * t706 + t630;
t732 = sin(pkin(11));
t733 = cos(pkin(11));
t673 = t697 * t732 + t698 * t733;
t616 = -0.2e1 * qJD(5) * t673 + t733 * t622 - t732 * t624;
t643 = t659 * t732 + t660 * t733;
t672 = t697 * t733 - t698 * t732;
t614 = (t672 * t706 - t643) * pkin(10) + (t672 * t673 + t683) * pkin(5) + t616;
t617 = 0.2e1 * qJD(5) * t672 + t732 * t622 + t733 * t624;
t642 = t659 * t733 - t660 * t732;
t663 = pkin(5) * t706 - pkin(10) * t673;
t671 = t672 ^ 2;
t615 = -pkin(5) * t671 + pkin(10) * t642 - t663 * t706 + t617;
t734 = sin(qJ(6));
t739 = cos(qJ(6));
t612 = t614 * t739 - t615 * t734;
t653 = t672 * t739 - t673 * t734;
t628 = qJD(6) * t653 + t642 * t734 + t643 * t739;
t654 = t672 * t734 + t673 * t739;
t637 = -mrSges(7,1) * t653 + mrSges(7,2) * t654;
t704 = qJD(6) + t706;
t647 = -mrSges(7,2) * t704 + mrSges(7,3) * t653;
t680 = qJDD(6) + t683;
t608 = m(7) * t612 + mrSges(7,1) * t680 - mrSges(7,3) * t628 - t637 * t654 + t647 * t704;
t613 = t614 * t734 + t615 * t739;
t627 = -qJD(6) * t654 + t642 * t739 - t643 * t734;
t648 = mrSges(7,1) * t704 - mrSges(7,3) * t654;
t609 = m(7) * t613 - mrSges(7,2) * t680 + mrSges(7,3) * t627 + t637 * t653 - t648 * t704;
t602 = t739 * t608 + t734 * t609;
t655 = -mrSges(6,1) * t672 + mrSges(6,2) * t673;
t661 = -mrSges(6,2) * t706 + mrSges(6,3) * t672;
t600 = m(6) * t616 + mrSges(6,1) * t683 - mrSges(6,3) * t643 - t655 * t673 + t661 * t706 + t602;
t662 = mrSges(6,1) * t706 - mrSges(6,3) * t673;
t751 = -t608 * t734 + t739 * t609;
t601 = m(6) * t617 - mrSges(6,2) * t683 + mrSges(6,3) * t642 + t655 * t672 - t662 * t706 + t751;
t596 = t733 * t600 + t732 * t601;
t677 = -mrSges(5,1) * t697 + mrSges(5,2) * t698;
t687 = -mrSges(5,2) * t706 + mrSges(5,3) * t697;
t594 = m(5) * t629 + mrSges(5,1) * t683 - mrSges(5,3) * t660 - t677 * t698 + t687 * t706 + t596;
t689 = mrSges(5,1) * t706 - mrSges(5,3) * t698;
t752 = -t600 * t732 + t733 * t601;
t595 = m(5) * t630 - mrSges(5,2) * t683 + mrSges(5,3) * t659 + t677 * t697 - t689 * t706 + t752;
t753 = -t594 * t735 + t740 * t595;
t587 = m(4) * t657 - mrSges(4,2) * t729 + mrSges(4,3) * t684 + t694 * t710 - t702 * t730 + t753;
t656 = t678 * t741 - t736 * t679;
t701 = -mrSges(4,2) * t730 + mrSges(4,3) * t710;
t645 = -pkin(3) * t729 - pkin(9) * t728 + t711 * t695 - t656;
t631 = -pkin(4) * t659 - qJ(5) * t696 + t698 * t688 + qJDD(5) + t645;
t619 = -pkin(5) * t642 - pkin(10) * t671 + t663 * t673 + t631;
t750 = m(7) * t619 - t627 * mrSges(7,1) + t628 * mrSges(7,2) - t653 * t647 + t654 * t648;
t747 = m(6) * t631 - t642 * mrSges(6,1) + mrSges(6,2) * t643 - t672 * t661 + t662 * t673 + t750;
t745 = -m(5) * t645 + t659 * mrSges(5,1) - mrSges(5,2) * t660 + t697 * t687 - t689 * t698 - t747;
t611 = m(4) * t656 + mrSges(4,1) * t729 - mrSges(4,3) * t685 - t694 * t711 + t701 * t730 + t745;
t582 = t736 * t587 + t741 * t611;
t699 = -t742 * g(3) - t761;
t718 = (-mrSges(3,1) * t742 + mrSges(3,2) * t737) * qJD(1);
t722 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t758;
t580 = m(3) * t699 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t719 + qJD(2) * t722 - t718 * t759 + t582;
t721 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t759;
t754 = t741 * t587 - t611 * t736;
t581 = m(3) * t700 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t720 - qJD(2) * t721 + t718 * t758 + t754;
t755 = -t580 * t737 + t742 * t581;
t572 = m(2) * t725 - mrSges(2,1) * t744 - qJDD(1) * mrSges(2,2) + t755;
t712 = -t744 * pkin(7) + t749;
t588 = t740 * t594 + t735 * t595;
t748 = m(4) * t686 - t684 * mrSges(4,1) + mrSges(4,2) * t685 - t710 * t701 + t702 * t711 + t588;
t746 = -m(3) * t712 + t720 * mrSges(3,1) - mrSges(3,2) * t719 - t721 * t759 + t722 * t758 - t748;
t584 = m(2) * t724 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t744 + t746;
t760 = t738 * t572 + t743 * t584;
t573 = t742 * t580 + t737 * t581;
t756 = t743 * t572 - t584 * t738;
t709 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t737 + Ifges(3,4) * t742) * qJD(1);
t708 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t737 + Ifges(3,2) * t742) * qJD(1);
t707 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t737 + Ifges(3,6) * t742) * qJD(1);
t692 = Ifges(4,1) * t711 + Ifges(4,4) * t710 + Ifges(4,5) * t730;
t691 = Ifges(4,4) * t711 + Ifges(4,2) * t710 + Ifges(4,6) * t730;
t690 = Ifges(4,5) * t711 + Ifges(4,6) * t710 + Ifges(4,3) * t730;
t666 = Ifges(5,1) * t698 + Ifges(5,4) * t697 + Ifges(5,5) * t706;
t665 = Ifges(5,4) * t698 + Ifges(5,2) * t697 + Ifges(5,6) * t706;
t664 = Ifges(5,5) * t698 + Ifges(5,6) * t697 + Ifges(5,3) * t706;
t651 = Ifges(6,1) * t673 + Ifges(6,4) * t672 + Ifges(6,5) * t706;
t650 = Ifges(6,4) * t673 + Ifges(6,2) * t672 + Ifges(6,6) * t706;
t649 = Ifges(6,5) * t673 + Ifges(6,6) * t672 + Ifges(6,3) * t706;
t634 = Ifges(7,1) * t654 + Ifges(7,4) * t653 + Ifges(7,5) * t704;
t633 = Ifges(7,4) * t654 + Ifges(7,2) * t653 + Ifges(7,6) * t704;
t632 = Ifges(7,5) * t654 + Ifges(7,6) * t653 + Ifges(7,3) * t704;
t604 = mrSges(7,2) * t619 - mrSges(7,3) * t612 + Ifges(7,1) * t628 + Ifges(7,4) * t627 + Ifges(7,5) * t680 + t632 * t653 - t633 * t704;
t603 = -mrSges(7,1) * t619 + mrSges(7,3) * t613 + Ifges(7,4) * t628 + Ifges(7,2) * t627 + Ifges(7,6) * t680 - t632 * t654 + t634 * t704;
t590 = mrSges(6,2) * t631 - mrSges(6,3) * t616 + Ifges(6,1) * t643 + Ifges(6,4) * t642 + Ifges(6,5) * t683 - pkin(10) * t602 - t603 * t734 + t604 * t739 + t649 * t672 - t650 * t706;
t589 = -mrSges(6,1) * t631 + mrSges(6,3) * t617 + Ifges(6,4) * t643 + Ifges(6,2) * t642 + Ifges(6,6) * t683 - pkin(5) * t750 + pkin(10) * t751 + t739 * t603 + t734 * t604 - t673 * t649 + t706 * t651;
t576 = mrSges(5,2) * t645 - mrSges(5,3) * t629 + Ifges(5,1) * t660 + Ifges(5,4) * t659 + Ifges(5,5) * t683 - qJ(5) * t596 - t589 * t732 + t590 * t733 + t664 * t697 - t665 * t706;
t575 = -mrSges(5,1) * t645 + mrSges(5,3) * t630 + Ifges(5,4) * t660 + Ifges(5,2) * t659 + Ifges(5,6) * t683 - pkin(4) * t747 + qJ(5) * t752 + t733 * t589 + t732 * t590 - t698 * t664 + t706 * t666;
t574 = (-Ifges(6,3) - Ifges(5,3)) * t683 + Ifges(4,6) * t729 + t730 * t692 - t711 * t690 + t697 * t666 - t698 * t665 + Ifges(4,2) * t684 + Ifges(4,4) * t685 - mrSges(4,1) * t686 - Ifges(7,3) * t680 + t672 * t651 - t673 * t650 - Ifges(5,6) * t659 - Ifges(5,5) * t660 - t654 * t633 + mrSges(4,3) * t657 + t653 * t634 - Ifges(6,6) * t642 - Ifges(6,5) * t643 + mrSges(5,2) * t630 - Ifges(7,5) * t628 - mrSges(5,1) * t629 - Ifges(7,6) * t627 + mrSges(6,2) * t617 - mrSges(6,1) * t616 + mrSges(7,2) * t613 - mrSges(7,1) * t612 - pkin(5) * t602 - pkin(4) * t596 - pkin(3) * t588;
t569 = mrSges(4,2) * t686 - mrSges(4,3) * t656 + Ifges(4,1) * t685 + Ifges(4,4) * t684 + Ifges(4,5) * t729 - pkin(9) * t588 - t575 * t735 + t576 * t740 + t690 * t710 - t691 * t730;
t568 = mrSges(3,2) * t712 - mrSges(3,3) * t699 + Ifges(3,1) * t719 + Ifges(3,4) * t720 + Ifges(3,5) * qJDD(2) - pkin(8) * t582 - qJD(2) * t708 + t569 * t741 - t574 * t736 + t707 * t758;
t567 = -pkin(1) * t573 + mrSges(2,3) * t725 - pkin(2) * t582 - Ifges(3,5) * t719 - Ifges(3,6) * t720 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t699 + mrSges(3,2) * t700 - t735 * t576 - t740 * t575 - pkin(3) * t745 - pkin(9) * t753 - Ifges(4,5) * t685 - Ifges(4,6) * t684 - Ifges(4,3) * t729 - mrSges(4,1) * t656 + mrSges(4,2) * t657 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + t744 * Ifges(2,5) - t711 * t691 + t710 * t692 + (-t708 * t737 + t709 * t742) * qJD(1);
t566 = -mrSges(3,1) * t712 + mrSges(3,3) * t700 + Ifges(3,4) * t719 + Ifges(3,2) * t720 + Ifges(3,6) * qJDD(2) - pkin(2) * t748 + pkin(8) * t754 + qJD(2) * t709 + t736 * t569 + t741 * t574 - t707 * t759;
t565 = -mrSges(2,2) * g(3) - mrSges(2,3) * t724 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t744 - pkin(7) * t573 - t566 * t737 + t568 * t742;
t1 = [-m(1) * g(1) + t756; -m(1) * g(2) + t760; (-m(1) - m(2)) * g(3) + t573; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t760 + t743 * t565 - t738 * t567; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t756 + t738 * t565 + t743 * t567; -mrSges(1,1) * g(2) + mrSges(2,1) * t724 + mrSges(1,2) * g(1) - mrSges(2,2) * t725 + Ifges(2,3) * qJDD(1) + pkin(1) * t746 + pkin(7) * t755 + t742 * t566 + t737 * t568;];
tauB  = t1;
