% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:40:10
% EndTime: 2019-05-08 05:40:52
% DurationCPUTime: 26.38s
% Computational Cost: add. (441580->375), mult. (938700->473), div. (0->0), fcn. (756505->12), ass. (0->154)
t797 = Ifges(6,1) + Ifges(7,1);
t790 = Ifges(6,4) - Ifges(7,5);
t796 = -Ifges(6,5) - Ifges(7,4);
t795 = Ifges(6,2) + Ifges(7,3);
t788 = Ifges(6,6) - Ifges(7,6);
t794 = -Ifges(6,3) - Ifges(7,2);
t793 = cos(qJ(5));
t752 = cos(pkin(6));
t792 = t752 * g(3);
t791 = -mrSges(6,3) - mrSges(7,2);
t751 = sin(pkin(6));
t756 = sin(qJ(2));
t787 = t751 * t756;
t760 = cos(qJ(2));
t786 = t751 * t760;
t785 = t752 * t756;
t784 = t752 * t760;
t757 = sin(qJ(1));
t761 = cos(qJ(1));
t743 = t757 * g(1) - t761 * g(2);
t762 = qJD(1) ^ 2;
t733 = t762 * t751 * pkin(8) + qJDD(1) * pkin(1) + t743;
t744 = -t761 * g(1) - t757 * g(2);
t775 = qJDD(1) * t751;
t734 = -t762 * pkin(1) + pkin(8) * t775 + t744;
t778 = t733 * t785 + t760 * t734;
t708 = -g(3) * t787 + t778;
t748 = t752 * qJD(1) + qJD(2);
t777 = qJD(1) * t751;
t773 = t756 * t777;
t731 = t748 * mrSges(3,1) - mrSges(3,3) * t773;
t735 = (-mrSges(3,1) * t760 + mrSges(3,2) * t756) * t777;
t738 = -qJD(2) * t773 + t760 * t775;
t747 = t752 * qJDD(1) + qJDD(2);
t736 = (-pkin(2) * t760 - pkin(9) * t756) * t777;
t746 = t748 ^ 2;
t776 = qJD(1) * t760;
t693 = -t746 * pkin(2) + t747 * pkin(9) + (-g(3) * t756 + t736 * t776) * t751 + t778;
t737 = (qJD(2) * t776 + qJDD(1) * t756) * t751;
t694 = -t738 * pkin(2) - t737 * pkin(9) - t792 + (-t733 + (pkin(2) * t756 - pkin(9) * t760) * t748 * qJD(1)) * t751;
t755 = sin(qJ(3));
t759 = cos(qJ(3));
t657 = -t755 * t693 + t759 * t694;
t725 = t759 * t748 - t755 * t773;
t706 = t725 * qJD(3) + t759 * t737 + t755 * t747;
t726 = t755 * t748 + t759 * t773;
t730 = qJDD(3) - t738;
t772 = t751 * t776;
t742 = qJD(3) - t772;
t648 = (t725 * t742 - t706) * pkin(10) + (t725 * t726 + t730) * pkin(3) + t657;
t658 = t759 * t693 + t755 * t694;
t705 = -t726 * qJD(3) - t755 * t737 + t759 * t747;
t716 = t742 * pkin(3) - t726 * pkin(10);
t724 = t725 ^ 2;
t651 = -t724 * pkin(3) + t705 * pkin(10) - t742 * t716 + t658;
t754 = sin(qJ(4));
t758 = cos(qJ(4));
t646 = t754 * t648 + t758 * t651;
t712 = t754 * t725 + t758 * t726;
t673 = -t712 * qJD(4) + t758 * t705 - t754 * t706;
t711 = t758 * t725 - t754 * t726;
t687 = -t711 * mrSges(5,1) + t712 * mrSges(5,2);
t741 = qJD(4) + t742;
t698 = t741 * mrSges(5,1) - t712 * mrSges(5,3);
t729 = qJDD(4) + t730;
t688 = -t711 * pkin(4) - t712 * pkin(11);
t740 = t741 ^ 2;
t642 = -t740 * pkin(4) + t729 * pkin(11) + t711 * t688 + t646;
t707 = -g(3) * t786 + t733 * t784 - t756 * t734;
t692 = -t747 * pkin(2) - t746 * pkin(9) + t736 * t773 - t707;
t656 = -t705 * pkin(3) - t724 * pkin(10) + t726 * t716 + t692;
t674 = t711 * qJD(4) + t754 * t705 + t758 * t706;
t644 = (-t711 * t741 - t674) * pkin(11) + (t712 * t741 - t673) * pkin(4) + t656;
t753 = sin(qJ(5));
t639 = t793 * t642 + t753 * t644;
t696 = t793 * t712 + t753 * t741;
t654 = t696 * qJD(5) + t753 * t674 - t793 * t729;
t672 = qJDD(5) - t673;
t710 = qJD(5) - t711;
t681 = t710 * mrSges(6,1) - t696 * mrSges(6,3);
t695 = t753 * t712 - t793 * t741;
t676 = t695 * pkin(5) - t696 * qJ(6);
t709 = t710 ^ 2;
t635 = -t709 * pkin(5) + t672 * qJ(6) + 0.2e1 * qJD(6) * t710 - t695 * t676 + t639;
t682 = -t710 * mrSges(7,1) + t696 * mrSges(7,2);
t774 = m(7) * t635 + t672 * mrSges(7,3) + t710 * t682;
t677 = t695 * mrSges(7,1) - t696 * mrSges(7,3);
t779 = -t695 * mrSges(6,1) - t696 * mrSges(6,2) - t677;
t630 = m(6) * t639 - t672 * mrSges(6,2) + t791 * t654 - t710 * t681 + t779 * t695 + t774;
t638 = -t753 * t642 + t793 * t644;
t655 = -t695 * qJD(5) + t793 * t674 + t753 * t729;
t680 = -t710 * mrSges(6,2) - t695 * mrSges(6,3);
t636 = -t672 * pkin(5) - t709 * qJ(6) + t696 * t676 + qJDD(6) - t638;
t679 = -t695 * mrSges(7,2) + t710 * mrSges(7,3);
t767 = -m(7) * t636 + t672 * mrSges(7,1) + t710 * t679;
t632 = m(6) * t638 + t672 * mrSges(6,1) + t791 * t655 + t710 * t680 + t779 * t696 + t767;
t768 = t793 * t630 - t753 * t632;
t622 = m(5) * t646 - t729 * mrSges(5,2) + t673 * mrSges(5,3) + t711 * t687 - t741 * t698 + t768;
t645 = t758 * t648 - t754 * t651;
t697 = -t741 * mrSges(5,2) + t711 * mrSges(5,3);
t641 = -t729 * pkin(4) - t740 * pkin(11) + t712 * t688 - t645;
t637 = -0.2e1 * qJD(6) * t696 + (t695 * t710 - t655) * qJ(6) + (t696 * t710 + t654) * pkin(5) + t641;
t633 = m(7) * t637 + t654 * mrSges(7,1) - t655 * mrSges(7,3) + t695 * t679 - t696 * t682;
t764 = -m(6) * t641 - t654 * mrSges(6,1) - t655 * mrSges(6,2) - t695 * t680 - t696 * t681 - t633;
t627 = m(5) * t645 + t729 * mrSges(5,1) - t674 * mrSges(5,3) - t712 * t687 + t741 * t697 + t764;
t616 = t754 * t622 + t758 * t627;
t713 = -t725 * mrSges(4,1) + t726 * mrSges(4,2);
t714 = -t742 * mrSges(4,2) + t725 * mrSges(4,3);
t614 = m(4) * t657 + t730 * mrSges(4,1) - t706 * mrSges(4,3) - t726 * t713 + t742 * t714 + t616;
t715 = t742 * mrSges(4,1) - t726 * mrSges(4,3);
t769 = t758 * t622 - t754 * t627;
t615 = m(4) * t658 - t730 * mrSges(4,2) + t705 * mrSges(4,3) + t725 * t713 - t742 * t715 + t769;
t770 = -t755 * t614 + t759 * t615;
t605 = m(3) * t708 - t747 * mrSges(3,2) + t738 * mrSges(3,3) - t748 * t731 + t735 * t772 + t770;
t608 = t759 * t614 + t755 * t615;
t720 = -t751 * t733 - t792;
t732 = -t748 * mrSges(3,2) + mrSges(3,3) * t772;
t607 = m(3) * t720 - t738 * mrSges(3,1) + t737 * mrSges(3,2) + (t731 * t756 - t732 * t760) * t777 + t608;
t625 = t753 * t630 + t793 * t632;
t765 = m(5) * t656 - t673 * mrSges(5,1) + t674 * mrSges(5,2) - t711 * t697 + t712 * t698 + t625;
t763 = -m(4) * t692 + t705 * mrSges(4,1) - t706 * mrSges(4,2) + t725 * t714 - t726 * t715 - t765;
t619 = m(3) * t707 + t747 * mrSges(3,1) - t737 * mrSges(3,3) + t748 * t732 - t735 * t773 + t763;
t596 = t605 * t785 - t751 * t607 + t619 * t784;
t594 = m(2) * t743 + qJDD(1) * mrSges(2,1) - t762 * mrSges(2,2) + t596;
t601 = t760 * t605 - t756 * t619;
t600 = m(2) * t744 - t762 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t601;
t783 = t761 * t594 + t757 * t600;
t782 = t695 * t795 - t696 * t790 - t710 * t788;
t781 = t695 * t788 + t696 * t796 + t710 * t794;
t780 = -t790 * t695 + t696 * t797 - t796 * t710;
t595 = t605 * t787 + t752 * t607 + t619 * t786;
t771 = -t757 * t594 + t761 * t600;
t623 = -mrSges(6,1) * t641 - mrSges(7,1) * t637 + mrSges(7,2) * t635 + mrSges(6,3) * t639 - pkin(5) * t633 - t654 * t795 + t790 * t655 + t788 * t672 + t781 * t696 + t780 * t710;
t624 = mrSges(6,2) * t641 + mrSges(7,2) * t636 - mrSges(6,3) * t638 - mrSges(7,3) * t637 - qJ(6) * t633 - t790 * t654 + t655 * t797 - t672 * t796 + t781 * t695 + t782 * t710;
t683 = Ifges(5,5) * t712 + Ifges(5,6) * t711 + Ifges(5,3) * t741;
t684 = Ifges(5,4) * t712 + Ifges(5,2) * t711 + Ifges(5,6) * t741;
t609 = mrSges(5,2) * t656 - mrSges(5,3) * t645 + Ifges(5,1) * t674 + Ifges(5,4) * t673 + Ifges(5,5) * t729 - pkin(11) * t625 - t753 * t623 + t793 * t624 + t711 * t683 - t741 * t684;
t685 = Ifges(5,1) * t712 + Ifges(5,4) * t711 + Ifges(5,5) * t741;
t610 = Ifges(5,4) * t674 + Ifges(5,2) * t673 + Ifges(5,6) * t729 - t712 * t683 + t741 * t685 - mrSges(5,1) * t656 + mrSges(5,3) * t646 - mrSges(6,1) * t638 + mrSges(6,2) * t639 + mrSges(7,1) * t636 - mrSges(7,3) * t635 - pkin(5) * t767 - qJ(6) * t774 - pkin(4) * t625 + (pkin(5) * t677 + t782) * t696 + (qJ(6) * t677 - t780) * t695 + t794 * t672 + (pkin(5) * mrSges(7,2) + t796) * t655 + (qJ(6) * mrSges(7,2) + t788) * t654;
t699 = Ifges(4,5) * t726 + Ifges(4,6) * t725 + Ifges(4,3) * t742;
t701 = Ifges(4,1) * t726 + Ifges(4,4) * t725 + Ifges(4,5) * t742;
t592 = -mrSges(4,1) * t692 + mrSges(4,3) * t658 + Ifges(4,4) * t706 + Ifges(4,2) * t705 + Ifges(4,6) * t730 - pkin(3) * t765 + pkin(10) * t769 + t754 * t609 + t758 * t610 - t726 * t699 + t742 * t701;
t700 = Ifges(4,4) * t726 + Ifges(4,2) * t725 + Ifges(4,6) * t742;
t597 = mrSges(4,2) * t692 - mrSges(4,3) * t657 + Ifges(4,1) * t706 + Ifges(4,4) * t705 + Ifges(4,5) * t730 - pkin(10) * t616 + t758 * t609 - t754 * t610 + t725 * t699 - t742 * t700;
t717 = Ifges(3,3) * t748 + (Ifges(3,5) * t756 + Ifges(3,6) * t760) * t777;
t718 = Ifges(3,6) * t748 + (Ifges(3,4) * t756 + Ifges(3,2) * t760) * t777;
t590 = mrSges(3,2) * t720 - mrSges(3,3) * t707 + Ifges(3,1) * t737 + Ifges(3,4) * t738 + Ifges(3,5) * t747 - pkin(9) * t608 - t755 * t592 + t759 * t597 + t717 * t772 - t748 * t718;
t719 = Ifges(3,5) * t748 + (Ifges(3,1) * t756 + Ifges(3,4) * t760) * t777;
t591 = -pkin(4) * t764 - pkin(3) * t616 - pkin(2) * t608 - t753 * t624 + Ifges(3,6) * t747 + t748 * t719 + Ifges(3,4) * t737 + Ifges(3,2) * t738 + t725 * t701 - t726 * t700 - Ifges(5,3) * t729 - Ifges(4,3) * t730 - t712 * t684 - mrSges(3,1) * t720 + t711 * t685 - Ifges(4,6) * t705 - Ifges(4,5) * t706 + mrSges(3,3) * t708 - Ifges(5,6) * t673 - Ifges(5,5) * t674 - mrSges(4,1) * t657 + mrSges(4,2) * t658 + mrSges(5,2) * t646 - mrSges(5,1) * t645 - pkin(11) * t768 - t717 * t773 - t793 * t623;
t766 = pkin(8) * t601 + t590 * t756 + t591 * t760;
t589 = Ifges(3,5) * t737 + Ifges(3,6) * t738 + Ifges(3,3) * t747 + mrSges(3,1) * t707 - mrSges(3,2) * t708 + t755 * t597 + t759 * t592 + pkin(2) * t763 + pkin(9) * t770 + (t718 * t756 - t719 * t760) * t777;
t588 = -mrSges(2,2) * g(3) - mrSges(2,3) * t743 + Ifges(2,5) * qJDD(1) - t762 * Ifges(2,6) + t760 * t590 - t756 * t591 + (-t595 * t751 - t596 * t752) * pkin(8);
t587 = mrSges(2,1) * g(3) + mrSges(2,3) * t744 + t762 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t595 - t751 * t589 + t766 * t752;
t1 = [-m(1) * g(1) + t771; -m(1) * g(2) + t783; (-m(1) - m(2)) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t783 - t757 * t587 + t761 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t771 + t761 * t587 + t757 * t588; -mrSges(1,1) * g(2) + mrSges(2,1) * t743 + mrSges(1,2) * g(1) - mrSges(2,2) * t744 + Ifges(2,3) * qJDD(1) + pkin(1) * t596 + t752 * t589 + t766 * t751;];
tauB  = t1;
