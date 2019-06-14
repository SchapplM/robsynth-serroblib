% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-05-06 12:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:42:05
% EndTime: 2019-05-06 12:42:11
% DurationCPUTime: 3.75s
% Computational Cost: add. (22582->329), mult. (45462->366), div. (0->0), fcn. (24482->6), ass. (0->126)
t760 = cos(qJ(2));
t757 = sin(qJ(2));
t793 = qJD(1) * qJD(2);
t782 = t757 * t793;
t730 = qJDD(1) * t760 - t782;
t795 = qJD(1) * t757;
t738 = pkin(3) * t795 - qJD(2) * pkin(8);
t754 = t760 ^ 2;
t763 = qJD(1) ^ 2;
t758 = sin(qJ(1));
t761 = cos(qJ(1));
t740 = -g(1) * t761 - g(2) * t758;
t712 = -pkin(1) * t763 + qJDD(1) * pkin(7) + t740;
t699 = -t757 * g(3) + t760 * t712;
t726 = (-pkin(2) * t760 - qJ(3) * t757) * qJD(1);
t762 = qJD(2) ^ 2;
t794 = qJD(1) * t760;
t775 = t762 * pkin(2) - qJDD(2) * qJ(3) - t726 * t794 - t699;
t768 = t754 * t763 * pkin(8) - t730 * pkin(3) - qJD(2) * t738 + t775;
t792 = qJD(3) * qJD(2);
t657 = 0.2e1 * t792 - t768;
t756 = sin(qJ(4));
t759 = cos(qJ(4));
t725 = qJD(2) * t759 - t756 * t794;
t683 = qJD(4) * t725 + qJDD(2) * t756 + t759 * t730;
t744 = qJD(4) + t795;
t812 = -0.2e1 * t725;
t724 = qJD(2) * t756 + t759 * t794;
t684 = -qJD(4) * t724 + qJDD(2) * t759 - t730 * t756;
t804 = t724 * t744;
t814 = (-t684 + t804) * qJ(5);
t649 = qJD(5) * t812 + t814 + (t725 * t744 + t683) * pkin(4) + t657;
t693 = -mrSges(6,2) * t724 + mrSges(6,3) * t744;
t695 = -mrSges(7,1) * t744 - mrSges(7,3) * t725;
t697 = -mrSges(6,1) * t744 + mrSges(6,2) * t725;
t694 = -pkin(5) * t744 - qJ(6) * t725;
t722 = t724 ^ 2;
t748 = -0.2e1 * t792;
t811 = 0.2e1 * qJD(5);
t644 = -t722 * qJ(6) + qJDD(6) + t748 + (-pkin(4) - pkin(5)) * t683 - t814 + (-pkin(4) * t744 + t694 + t811) * t725 + t768;
t691 = mrSges(7,2) * t744 + mrSges(7,3) * t724;
t776 = m(7) * t644 - t683 * mrSges(7,1) - t724 * t691;
t639 = m(6) * t649 + t683 * mrSges(6,1) + t724 * t693 - t776 - (t695 + t697) * t725 - (mrSges(7,2) + mrSges(6,3)) * t684;
t692 = -mrSges(5,2) * t744 - mrSges(5,3) * t724;
t696 = mrSges(5,1) * t744 - mrSges(5,3) * t725;
t822 = -m(5) * t657 - t683 * mrSges(5,1) - t684 * mrSges(5,2) - t724 * t692 - t725 * t696 - t639;
t660 = t748 + t775;
t727 = (mrSges(4,2) * t760 - mrSges(4,3) * t757) * qJD(1);
t736 = mrSges(4,1) * t795 + qJD(2) * mrSges(4,2);
t821 = -m(4) * t660 + qJDD(2) * mrSges(4,3) + qJD(2) * t736 + t727 * t794 - t822;
t820 = Ifges(3,1) + Ifges(4,2);
t807 = Ifges(3,4) + Ifges(4,6);
t806 = Ifges(3,5) - Ifges(4,4);
t819 = Ifges(3,2) + Ifges(4,3);
t805 = Ifges(3,6) - Ifges(4,5);
t818 = Ifges(3,3) + Ifges(4,1);
t817 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t791 = -Ifges(5,4) + Ifges(6,5) + Ifges(7,4);
t790 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t816 = Ifges(5,2) + Ifges(6,3) + Ifges(7,2);
t815 = -Ifges(6,2) - Ifges(5,3) - Ifges(7,3);
t789 = Ifges(6,6) - Ifges(5,6) - Ifges(7,6);
t810 = t763 * pkin(7);
t808 = -mrSges(5,3) - mrSges(6,2);
t803 = t744 * t691;
t698 = -t760 * g(3) - t757 * t712;
t728 = (-mrSges(3,1) * t760 + mrSges(3,2) * t757) * qJD(1);
t783 = t760 * t793;
t729 = qJDD(1) * t757 + t783;
t734 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t794;
t735 = -mrSges(4,1) * t794 - qJD(2) * mrSges(4,3);
t739 = t758 * g(1) - t761 * g(2);
t774 = -qJDD(1) * pkin(1) - t739;
t767 = pkin(2) * t782 - 0.2e1 * qJD(3) * t795 + (-t729 - t783) * qJ(3) + t774;
t654 = -t738 * t795 + (-pkin(3) * t754 - pkin(7)) * t763 + (-pkin(2) - pkin(8)) * t730 + t767;
t661 = -qJDD(2) * pkin(2) - t762 * qJ(3) + t726 * t795 + qJDD(3) - t698;
t658 = (-t757 * t760 * t763 - qJDD(2)) * pkin(8) + (t729 - t783) * pkin(3) + t661;
t650 = -t756 * t654 + t759 * t658;
t723 = qJDD(4) + t729;
t687 = pkin(4) * t724 - qJ(5) * t725;
t741 = t744 ^ 2;
t647 = -t723 * pkin(4) - t741 * qJ(5) + t725 * t687 + qJDD(5) - t650;
t641 = qJD(6) * t812 + (-t684 - t804) * qJ(6) + (t724 * t725 - t723) * pkin(5) + t647;
t689 = -mrSges(7,1) * t724 + mrSges(7,2) * t725;
t777 = -m(7) * t641 + t684 * mrSges(7,3) + t725 * t689;
t769 = -m(6) * t647 + t723 * mrSges(6,1) + t744 * t693 + t777;
t688 = mrSges(6,1) * t724 - mrSges(6,3) * t725;
t800 = -mrSges(5,1) * t724 - mrSges(5,2) * t725 - t688;
t637 = m(5) * t650 + (t691 + t692) * t744 + t800 * t725 + (mrSges(5,1) + mrSges(7,1)) * t723 + t808 * t684 + t769;
t651 = t759 * t654 + t756 * t658;
t646 = -pkin(4) * t741 + t723 * qJ(5) - t724 * t687 + t744 * t811 + t651;
t643 = -pkin(5) * t722 + qJ(6) * t683 + 0.2e1 * qJD(6) * t724 + t694 * t744 + t646;
t788 = m(7) * t643 + t683 * mrSges(7,3) + t724 * t689;
t773 = m(6) * t646 + t723 * mrSges(6,3) + t744 * t697 + t788;
t638 = m(5) * t651 + (t695 - t696) * t744 + t800 * t724 + (-mrSges(5,2) + mrSges(7,2)) * t723 + t808 * t683 + t773;
t631 = t759 * t637 + t756 * t638;
t770 = -m(4) * t661 - t729 * mrSges(4,1) - t631;
t629 = m(3) * t698 - t729 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t734 - t735) * qJD(2) + (-t727 - t728) * t795 + t770;
t733 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t795;
t634 = (mrSges(3,3) + mrSges(4,1)) * t730 + t728 * t794 - qJDD(2) * mrSges(3,2) - qJD(2) * t733 + m(3) * t699 + t821;
t779 = -t629 * t757 + t760 * t634;
t622 = m(2) * t740 - mrSges(2,1) * t763 - qJDD(1) * mrSges(2,2) + t779;
t711 = t774 - t810;
t659 = -t730 * pkin(2) + t767 - t810;
t801 = -t756 * t637 + t759 * t638;
t772 = -m(4) * t659 - t730 * mrSges(4,2) + t736 * t795 - t801;
t764 = -m(3) * t711 + t734 * t794 + t730 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t729 + (-t733 * t757 - t735 * t760) * qJD(1) + t772;
t627 = m(2) * t739 + qJDD(1) * mrSges(2,1) - t763 * mrSges(2,2) + t764;
t802 = t758 * t622 + t761 * t627;
t623 = t760 * t629 + t757 * t634;
t798 = t818 * qJD(2) + (t757 * t806 + t760 * t805) * qJD(1);
t797 = -t805 * qJD(2) + (-t757 * t807 - t760 * t819) * qJD(1);
t796 = t806 * qJD(2) + (t757 * t820 + t807 * t760) * qJD(1);
t787 = -t724 * t789 - t725 * t790 + t744 * t815;
t786 = t724 * t816 + t725 * t791 + t744 * t789;
t785 = t724 * t791 + t725 * t817 + t744 * t790;
t780 = t761 * t622 - t627 * t758;
t640 = -t723 * mrSges(7,1) - t777 - t803;
t630 = -t729 * mrSges(4,3) + t735 * t794 - t772;
t625 = mrSges(5,2) * t657 + mrSges(6,2) * t647 + mrSges(7,2) * t644 - mrSges(5,3) * t650 - mrSges(6,3) * t649 - mrSges(7,3) * t641 - qJ(5) * t639 - qJ(6) * t640 + t791 * t683 + t684 * t817 + t790 * t723 + t787 * t724 + t786 * t744;
t624 = -mrSges(5,1) * t657 + mrSges(5,3) * t651 - mrSges(6,1) * t649 + mrSges(6,2) * t646 + mrSges(7,1) * t644 - mrSges(7,3) * t643 + pkin(5) * t776 - qJ(6) * t788 - pkin(4) * t639 + (-qJ(6) * t695 + t785) * t744 + (pkin(5) * t695 + t787) * t725 + (-mrSges(7,2) * qJ(6) - t789) * t723 + (mrSges(7,2) * pkin(5) - t791) * t684 - t816 * t683;
t619 = qJ(5) * (t744 * t695 + t773) + (-qJ(5) * t688 + t785) * t724 + (-pkin(4) * t688 - t786) * t725 + (-mrSges(6,2) * qJ(5) + t789) * t683 + (-mrSges(6,2) * pkin(4) + t790) * t684 + (mrSges(7,1) * pkin(4) + mrSges(7,2) * qJ(5) - t815) * t723 + t820 * t729 + pkin(4) * (t769 + t803) + t806 * qJDD(2) + t807 * t730 + t797 * qJD(2) + t798 * t794 + mrSges(3,2) * t711 - mrSges(3,3) * t698 - mrSges(4,3) * t659 + mrSges(4,1) * t661 + mrSges(5,1) * t650 - mrSges(5,2) * t651 + mrSges(7,2) * t643 + mrSges(6,3) * t646 - mrSges(6,1) * t647 - pkin(5) * t640 - mrSges(7,1) * t641 + pkin(3) * t631 - qJ(3) * t630;
t618 = -mrSges(3,1) * t711 - mrSges(4,1) * t660 + mrSges(4,2) * t659 + mrSges(3,3) * t699 - pkin(2) * t630 - pkin(3) * t822 - pkin(8) * t801 + t796 * qJD(2) + t805 * qJDD(2) - t759 * t624 - t756 * t625 + t807 * t729 + t819 * t730 - t798 * t795;
t617 = -pkin(1) * t623 + mrSges(2,3) * t740 - pkin(2) * (-qJD(2) * t735 + t770) - qJ(3) * t821 + t756 * t624 + pkin(8) * t631 - mrSges(3,1) * t698 + mrSges(3,2) * t699 - mrSges(4,2) * t661 + mrSges(4,3) * t660 - t759 * t625 + mrSges(2,1) * g(3) + t763 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-mrSges(4,1) * qJ(3) - t805) * t730 - t806 * t729 + (mrSges(4,2) * pkin(2) - t818) * qJDD(2) + (t796 * t760 + (pkin(2) * t727 + t797) * t757) * qJD(1);
t616 = -mrSges(2,2) * g(3) - mrSges(2,3) * t739 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t763 - pkin(7) * t623 - t618 * t757 + t619 * t760;
t1 = [-m(1) * g(1) + t780; -m(1) * g(2) + t802; (-m(1) - m(2)) * g(3) + t623; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t802 + t761 * t616 - t758 * t617; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t780 + t758 * t616 + t761 * t617; -mrSges(1,1) * g(2) + mrSges(2,1) * t739 + mrSges(1,2) * g(1) - mrSges(2,2) * t740 + Ifges(2,3) * qJDD(1) + pkin(1) * t764 + pkin(7) * t779 + t760 * t618 + t757 * t619;];
tauB  = t1;
