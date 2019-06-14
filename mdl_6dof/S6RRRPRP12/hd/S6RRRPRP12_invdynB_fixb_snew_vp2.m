% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:24:29
% EndTime: 2019-05-07 09:24:55
% DurationCPUTime: 9.58s
% Computational Cost: add. (142451->356), mult. (304954->431), div. (0->0), fcn. (231355->10), ass. (0->147)
t824 = Ifges(4,1) + Ifges(5,2);
t823 = Ifges(6,1) + Ifges(7,1);
t813 = Ifges(6,4) - Ifges(7,5);
t812 = Ifges(7,4) + Ifges(6,5);
t811 = -Ifges(4,5) + Ifges(5,4);
t822 = Ifges(4,2) + Ifges(5,3);
t821 = Ifges(6,2) + Ifges(7,3);
t810 = Ifges(4,6) - Ifges(5,5);
t809 = -Ifges(5,6) - Ifges(4,4);
t808 = Ifges(6,6) - Ifges(7,6);
t820 = Ifges(4,3) + Ifges(5,1);
t819 = Ifges(6,3) + Ifges(7,2);
t765 = sin(pkin(6));
t769 = sin(qJ(2));
t771 = cos(qJ(2));
t789 = qJD(1) * qJD(2);
t751 = (-qJDD(1) * t771 + t769 * t789) * t765;
t792 = qJD(1) * t765;
t749 = (-pkin(2) * t771 - pkin(9) * t769) * t792;
t766 = cos(pkin(6));
t762 = t766 * qJD(1) + qJD(2);
t760 = t762 ^ 2;
t761 = t766 * qJDD(1) + qJDD(2);
t791 = qJD(1) * t771;
t770 = sin(qJ(1));
t772 = cos(qJ(1));
t758 = t770 * g(1) - t772 * g(2);
t773 = qJD(1) ^ 2;
t816 = pkin(8) * t765;
t746 = qJDD(1) * pkin(1) + t773 * t816 + t758;
t759 = -t772 * g(1) - t770 * g(2);
t747 = -t773 * pkin(1) + qJDD(1) * t816 + t759;
t804 = t766 * t769;
t793 = t746 * t804 + t771 * t747;
t680 = -t760 * pkin(2) + t761 * pkin(9) + (-g(3) * t769 + t749 * t791) * t765 + t793;
t750 = (qJDD(1) * t769 + t771 * t789) * t765;
t815 = t766 * g(3);
t681 = t751 * pkin(2) - t750 * pkin(9) - t815 + (-t746 + (pkin(2) * t769 - pkin(9) * t771) * t762 * qJD(1)) * t765;
t768 = sin(qJ(3));
t818 = cos(qJ(3));
t663 = t818 * t680 + t768 * t681;
t787 = t769 * t792;
t737 = -t818 * t762 + t768 * t787;
t738 = t768 * t762 + t818 * t787;
t713 = t737 * pkin(3) - t738 * qJ(4);
t743 = qJDD(3) + t751;
t786 = t765 * t791;
t756 = -qJD(3) + t786;
t755 = t756 ^ 2;
t659 = t755 * pkin(3) - t743 * qJ(4) + 0.2e1 * qJD(4) * t756 + t737 * t713 - t663;
t817 = cos(qJ(5));
t814 = -mrSges(6,3) - mrSges(7,2);
t807 = t737 * t756;
t806 = t765 * t769;
t805 = t765 * t771;
t803 = t766 * t771;
t712 = -g(3) * t806 + t793;
t744 = t762 * mrSges(3,1) - mrSges(3,3) * t787;
t748 = (-mrSges(3,1) * t771 + mrSges(3,2) * t769) * t792;
t662 = -t768 * t680 + t818 * t681;
t710 = -t737 * qJD(3) + t818 * t750 + t768 * t761;
t714 = t737 * mrSges(4,1) + t738 * mrSges(4,2);
t720 = t756 * mrSges(4,2) - t737 * mrSges(4,3);
t722 = t737 * mrSges(5,1) + t756 * mrSges(5,3);
t660 = -t743 * pkin(3) - t755 * qJ(4) + t738 * t713 + qJDD(4) - t662;
t654 = (t737 * t738 - t743) * pkin(10) + (t710 - t807) * pkin(4) + t660;
t709 = t738 * qJD(3) + t768 * t750 - t818 * t761;
t724 = t738 * pkin(4) + t756 * pkin(10);
t736 = t737 ^ 2;
t711 = -g(3) * t805 + t746 * t803 - t769 * t747;
t679 = -t761 * pkin(2) - t760 * pkin(9) + t749 * t787 - t711;
t775 = (-t710 - t807) * qJ(4) + t679 + (-t756 * pkin(3) - 0.2e1 * qJD(4)) * t738;
t658 = -t736 * pkin(4) - t738 * t724 + (pkin(3) + pkin(10)) * t709 + t775;
t767 = sin(qJ(5));
t650 = t767 * t654 + t817 * t658;
t719 = t767 * t737 - t817 * t756;
t666 = t719 * qJD(5) - t817 * t709 + t767 * t743;
t735 = qJD(5) + t738;
t692 = t735 * mrSges(6,1) - t719 * mrSges(6,3);
t706 = qJDD(5) + t710;
t718 = -t817 * t737 - t767 * t756;
t684 = t718 * pkin(5) - t719 * qJ(6);
t734 = t735 ^ 2;
t647 = -t734 * pkin(5) + t706 * qJ(6) + 0.2e1 * qJD(6) * t735 - t718 * t684 + t650;
t693 = -t735 * mrSges(7,1) + t719 * mrSges(7,2);
t788 = m(7) * t647 + t706 * mrSges(7,3) + t735 * t693;
t685 = t718 * mrSges(7,1) - t719 * mrSges(7,3);
t797 = -t718 * mrSges(6,1) - t719 * mrSges(6,2) - t685;
t642 = m(6) * t650 - t706 * mrSges(6,2) + t814 * t666 - t735 * t692 + t797 * t718 + t788;
t649 = t817 * t654 - t767 * t658;
t667 = -t718 * qJD(5) + t767 * t709 + t817 * t743;
t691 = -t735 * mrSges(6,2) - t718 * mrSges(6,3);
t648 = -t706 * pkin(5) - t734 * qJ(6) + t719 * t684 + qJDD(6) - t649;
t690 = -t718 * mrSges(7,2) + t735 * mrSges(7,3);
t783 = -m(7) * t648 + t706 * mrSges(7,1) + t735 * t690;
t644 = m(6) * t649 + t706 * mrSges(6,1) + t814 * t667 + t735 * t691 + t797 * t719 + t783;
t637 = t767 * t642 + t817 * t644;
t715 = -t737 * mrSges(5,2) - t738 * mrSges(5,3);
t778 = -m(5) * t660 - t710 * mrSges(5,1) - t738 * t715 - t637;
t633 = m(4) * t662 - t710 * mrSges(4,3) - t738 * t714 + (-t720 + t722) * t756 + (mrSges(4,1) - mrSges(5,2)) * t743 + t778;
t721 = -t756 * mrSges(4,1) - t738 * mrSges(4,3);
t723 = t738 * mrSges(5,1) - t756 * mrSges(5,2);
t657 = -t709 * pkin(4) - t736 * pkin(10) - t756 * t724 - t659;
t652 = -0.2e1 * qJD(6) * t719 + (t718 * t735 - t667) * qJ(6) + (t719 * t735 + t666) * pkin(5) + t657;
t645 = m(7) * t652 + t666 * mrSges(7,1) - t667 * mrSges(7,3) + t718 * t690 - t719 * t693;
t777 = m(6) * t657 + t666 * mrSges(6,1) + t667 * mrSges(6,2) + t718 * t691 + t719 * t692 + t645;
t774 = -m(5) * t659 + t743 * mrSges(5,3) - t756 * t723 + t777;
t640 = (-mrSges(4,3) - mrSges(5,1)) * t709 + t756 * t721 + (-t714 - t715) * t737 - t743 * mrSges(4,2) + t774 + m(4) * t663;
t784 = -t768 * t633 + t818 * t640;
t625 = m(3) * t712 - t761 * mrSges(3,2) - t751 * mrSges(3,3) - t762 * t744 + t748 * t786 + t784;
t628 = t818 * t633 + t768 * t640;
t729 = -t765 * t746 - t815;
t745 = -t762 * mrSges(3,2) + mrSges(3,3) * t786;
t627 = m(3) * t729 + t751 * mrSges(3,1) + t750 * mrSges(3,2) + (t744 * t769 - t745 * t771) * t792 + t628;
t661 = t709 * pkin(3) + t775;
t801 = t817 * t642 - t767 * t644;
t782 = -m(5) * t661 + t709 * mrSges(5,2) + t737 * t722 - t801;
t776 = -m(4) * t679 - t709 * mrSges(4,1) - t737 * t720 + (-t721 + t723) * t738 + (-mrSges(4,2) + mrSges(5,3)) * t710 + t782;
t631 = m(3) * t711 + t761 * mrSges(3,1) - t750 * mrSges(3,3) + t762 * t745 - t748 * t787 + t776;
t616 = t625 * t804 - t765 * t627 + t631 * t803;
t614 = m(2) * t758 + qJDD(1) * mrSges(2,1) - t773 * mrSges(2,2) + t616;
t621 = t771 * t625 - t769 * t631;
t620 = m(2) * t759 - t773 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t621;
t802 = t772 * t614 + t770 * t620;
t800 = t821 * t718 - t813 * t719 - t808 * t735;
t799 = t808 * t718 - t812 * t719 - t819 * t735;
t798 = -t813 * t718 + t823 * t719 + t812 * t735;
t796 = t810 * t737 + t811 * t738 + t820 * t756;
t795 = t822 * t737 + t809 * t738 + t810 * t756;
t794 = -t809 * t737 - t824 * t738 - t811 * t756;
t615 = t625 * t806 + t766 * t627 + t631 * t805;
t785 = -t770 * t614 + t772 * t620;
t634 = -t710 * mrSges(5,3) - t738 * t723 - t782;
t635 = -mrSges(6,1) * t657 - mrSges(7,1) * t652 + mrSges(7,2) * t647 + mrSges(6,3) * t650 - pkin(5) * t645 - t821 * t666 + t813 * t667 + t808 * t706 + t799 * t719 + t798 * t735;
t636 = mrSges(6,2) * t657 + mrSges(7,2) * t648 - mrSges(6,3) * t649 - mrSges(7,3) * t652 - qJ(6) * t645 - t813 * t666 + t823 * t667 + t812 * t706 + t799 * t718 + t800 * t735;
t612 = -mrSges(4,1) * t679 - mrSges(5,1) * t659 + mrSges(5,2) * t661 + mrSges(4,3) * t663 - pkin(3) * t634 + pkin(4) * t777 - pkin(10) * t801 - t817 * t635 - t767 * t636 - t822 * t709 - t809 * t710 + t796 * t738 + t810 * t743 + t794 * t756;
t617 = pkin(4) * t637 + qJ(6) * t788 + pkin(5) * t783 + mrSges(4,2) * t679 + mrSges(5,1) * t660 - mrSges(5,3) * t661 - mrSges(4,3) * t662 - qJ(4) * t634 + mrSges(7,3) * t647 - mrSges(7,1) * t648 + mrSges(6,1) * t649 - mrSges(6,2) * t650 - t795 * t756 - t811 * t743 + t796 * t737 + (-pkin(5) * t685 - t800) * t719 + (-qJ(6) * t685 + t798) * t718 + t824 * t710 + t809 * t709 + t819 * t706 + (-pkin(5) * mrSges(7,2) + t812) * t667 + (-qJ(6) * mrSges(7,2) - t808) * t666;
t726 = Ifges(3,3) * t762 + (Ifges(3,5) * t769 + Ifges(3,6) * t771) * t792;
t727 = Ifges(3,6) * t762 + (Ifges(3,4) * t769 + Ifges(3,2) * t771) * t792;
t610 = mrSges(3,2) * t729 - mrSges(3,3) * t711 + Ifges(3,1) * t750 - Ifges(3,4) * t751 + Ifges(3,5) * t761 - pkin(9) * t628 - t768 * t612 + t818 * t617 + t726 * t786 - t762 * t727;
t728 = Ifges(3,5) * t762 + (Ifges(3,1) * t769 + Ifges(3,4) * t771) * t792;
t611 = -t726 * t787 + Ifges(3,6) * t761 + t762 * t728 + t767 * t635 - pkin(2) * t628 + pkin(10) * t637 + Ifges(3,4) * t750 - Ifges(3,2) * t751 + mrSges(5,3) * t659 - mrSges(5,2) * t660 - mrSges(4,1) * t662 + mrSges(4,2) * t663 - qJ(4) * t774 - mrSges(3,1) * t729 - pkin(3) * (t756 * t722 + t778) - t817 * t636 + mrSges(3,3) * t712 + (pkin(3) * mrSges(5,2) - t820) * t743 + t795 * t738 + (qJ(4) * t715 + t794) * t737 + t811 * t710 + (qJ(4) * mrSges(5,1) + t810) * t709;
t779 = pkin(8) * t621 + t610 * t769 + t611 * t771;
t609 = Ifges(3,5) * t750 - Ifges(3,6) * t751 + Ifges(3,3) * t761 + mrSges(3,1) * t711 - mrSges(3,2) * t712 + t768 * t617 + t818 * t612 + pkin(2) * t776 + pkin(9) * t784 + (t727 * t769 - t728 * t771) * t792;
t608 = -mrSges(2,2) * g(3) - mrSges(2,3) * t758 + Ifges(2,5) * qJDD(1) - t773 * Ifges(2,6) + t771 * t610 - t769 * t611 + (-t615 * t765 - t616 * t766) * pkin(8);
t607 = mrSges(2,1) * g(3) + mrSges(2,3) * t759 + t773 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t615 - t765 * t609 + t779 * t766;
t1 = [-m(1) * g(1) + t785; -m(1) * g(2) + t802; (-m(1) - m(2)) * g(3) + t615; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t802 - t770 * t607 + t772 * t608; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t785 + t772 * t607 + t770 * t608; -mrSges(1,1) * g(2) + mrSges(2,1) * t758 + mrSges(1,2) * g(1) - mrSges(2,2) * t759 + Ifges(2,3) * qJDD(1) + pkin(1) * t616 + t766 * t609 + t779 * t765;];
tauB  = t1;
