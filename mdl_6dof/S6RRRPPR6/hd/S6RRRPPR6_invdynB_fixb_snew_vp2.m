% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 05:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:28:03
% EndTime: 2019-05-07 05:28:26
% DurationCPUTime: 20.52s
% Computational Cost: add. (320523->381), mult. (710665->475), div. (0->0), fcn. (562429->12), ass. (0->159)
t828 = -2 * qJD(4);
t827 = Ifges(5,1) + Ifges(6,2);
t826 = Ifges(6,1) + Ifges(5,3);
t821 = Ifges(5,4) + Ifges(6,6);
t820 = Ifges(5,5) - Ifges(6,4);
t825 = -Ifges(5,2) - Ifges(6,3);
t819 = Ifges(5,6) - Ifges(6,5);
t776 = sin(pkin(6));
t780 = sin(qJ(2));
t784 = cos(qJ(2));
t802 = qJD(1) * qJD(2);
t764 = (-qJDD(1) * t784 + t780 * t802) * t776;
t805 = qJD(1) * t776;
t762 = (-pkin(2) * t784 - pkin(9) * t780) * t805;
t777 = cos(pkin(6));
t772 = qJD(1) * t777 + qJD(2);
t770 = t772 ^ 2;
t771 = qJDD(1) * t777 + qJDD(2);
t804 = qJD(1) * t784;
t781 = sin(qJ(1));
t785 = cos(qJ(1));
t768 = t781 * g(1) - g(2) * t785;
t786 = qJD(1) ^ 2;
t823 = pkin(8) * t776;
t759 = qJDD(1) * pkin(1) + t786 * t823 + t768;
t769 = -g(1) * t785 - g(2) * t781;
t760 = -pkin(1) * t786 + qJDD(1) * t823 + t769;
t814 = t777 * t780;
t806 = t759 * t814 + t784 * t760;
t710 = -t770 * pkin(2) + t771 * pkin(9) + (-g(3) * t780 + t762 * t804) * t776 + t806;
t763 = (qJDD(1) * t780 + t784 * t802) * t776;
t822 = t777 * g(3);
t711 = t764 * pkin(2) - t763 * pkin(9) - t822 + (-t759 + (pkin(2) * t780 - pkin(9) * t784) * t772 * qJD(1)) * t776;
t779 = sin(qJ(3));
t783 = cos(qJ(3));
t676 = -t779 * t710 + t783 * t711;
t801 = t780 * t805;
t751 = t772 * t783 - t779 * t801;
t729 = qJD(3) * t751 + t763 * t783 + t771 * t779;
t752 = t772 * t779 + t783 * t801;
t756 = qJDD(3) + t764;
t800 = t776 * t804;
t766 = -qJD(3) + t800;
t668 = (-t751 * t766 - t729) * qJ(4) + (t751 * t752 + t756) * pkin(3) + t676;
t677 = t783 * t710 + t779 * t711;
t728 = -qJD(3) * t752 - t763 * t779 + t771 * t783;
t741 = -pkin(3) * t766 - qJ(4) * t752;
t750 = t751 ^ 2;
t671 = -pkin(3) * t750 + qJ(4) * t728 + t741 * t766 + t677;
t775 = sin(pkin(11));
t818 = cos(pkin(11));
t738 = t775 * t751 + t752 * t818;
t663 = t668 * t818 - t775 * t671 + t738 * t828;
t824 = -2 * qJD(5);
t737 = -t751 * t818 + t752 * t775;
t817 = t737 * t766;
t816 = t776 * t780;
t815 = t776 * t784;
t813 = t777 * t784;
t731 = -g(3) * t816 + t806;
t757 = mrSges(3,1) * t772 - mrSges(3,3) * t801;
t761 = (-mrSges(3,1) * t784 + mrSges(3,2) * t780) * t805;
t693 = t775 * t728 + t729 * t818;
t704 = mrSges(5,1) * t737 + mrSges(5,2) * t738;
t715 = mrSges(6,1) * t737 + mrSges(6,3) * t766;
t717 = mrSges(5,2) * t766 - mrSges(5,3) * t737;
t703 = pkin(4) * t737 - qJ(5) * t738;
t765 = t766 ^ 2;
t662 = -t756 * pkin(4) - t765 * qJ(5) + t738 * t703 + qJDD(5) - t663;
t657 = (t737 * t738 - t756) * pkin(10) + (t693 - t817) * pkin(5) + t662;
t692 = -t728 * t818 + t729 * t775;
t719 = pkin(5) * t738 + pkin(10) * t766;
t736 = t737 ^ 2;
t730 = -g(3) * t815 + t759 * t813 - t780 * t760;
t709 = -t771 * pkin(2) - t770 * pkin(9) + t762 * t801 - t730;
t672 = -t728 * pkin(3) - t750 * qJ(4) + t752 * t741 + qJDD(4) + t709;
t787 = (-t693 - t817) * qJ(5) + t672 + (-pkin(4) * t766 + t824) * t738;
t660 = t787 - t736 * pkin(5) - t738 * t719 + (pkin(4) + pkin(10)) * t692;
t778 = sin(qJ(6));
t782 = cos(qJ(6));
t655 = t657 * t782 - t660 * t778;
t713 = t737 * t782 + t766 * t778;
t675 = qJD(6) * t713 + t692 * t778 + t756 * t782;
t714 = t737 * t778 - t766 * t782;
t684 = -mrSges(7,1) * t713 + mrSges(7,2) * t714;
t735 = qJD(6) + t738;
t685 = -mrSges(7,2) * t735 + mrSges(7,3) * t713;
t691 = qJDD(6) + t693;
t653 = m(7) * t655 + mrSges(7,1) * t691 - mrSges(7,3) * t675 - t684 * t714 + t685 * t735;
t656 = t657 * t778 + t660 * t782;
t674 = -qJD(6) * t714 + t692 * t782 - t756 * t778;
t686 = mrSges(7,1) * t735 - mrSges(7,3) * t714;
t654 = m(7) * t656 - mrSges(7,2) * t691 + mrSges(7,3) * t674 + t684 * t713 - t686 * t735;
t645 = t782 * t653 + t778 * t654;
t705 = -mrSges(6,2) * t737 - mrSges(6,3) * t738;
t791 = -m(6) * t662 - t693 * mrSges(6,1) - t738 * t705 - t645;
t640 = m(5) * t663 - t693 * mrSges(5,3) - t738 * t704 + (t715 - t717) * t766 + (mrSges(5,1) - mrSges(6,2)) * t756 + t791;
t733 = t737 * t828;
t810 = t775 * t668 + t818 * t671;
t664 = t733 + t810;
t718 = -mrSges(5,1) * t766 - mrSges(5,3) * t738;
t794 = t765 * pkin(4) - t756 * qJ(5) - t810;
t661 = 0.2e1 * qJD(5) * t766 + ((2 * qJD(4)) + t703) * t737 + t794;
t716 = mrSges(6,1) * t738 - mrSges(6,2) * t766;
t659 = -t692 * pkin(5) - t736 * pkin(10) - t737 * t703 + t733 + (t824 - t719) * t766 - t794;
t792 = -m(7) * t659 + t674 * mrSges(7,1) - t675 * mrSges(7,2) + t713 * t685 - t714 * t686;
t790 = -m(6) * t661 + t756 * mrSges(6,3) - t766 * t716 - t792;
t650 = m(5) * t664 - t756 * mrSges(5,2) + t766 * t718 + (-t704 - t705) * t737 + (-mrSges(5,3) - mrSges(6,1)) * t692 + t790;
t638 = t818 * t640 + t775 * t650;
t739 = -mrSges(4,1) * t751 + mrSges(4,2) * t752;
t740 = mrSges(4,2) * t766 + mrSges(4,3) * t751;
t636 = m(4) * t676 + mrSges(4,1) * t756 - mrSges(4,3) * t729 - t739 * t752 - t740 * t766 + t638;
t742 = -mrSges(4,1) * t766 - mrSges(4,3) * t752;
t797 = -t640 * t775 + t818 * t650;
t637 = m(4) * t677 - mrSges(4,2) * t756 + mrSges(4,3) * t728 + t739 * t751 + t742 * t766 + t797;
t798 = -t636 * t779 + t783 * t637;
t628 = m(3) * t731 - mrSges(3,2) * t771 - mrSges(3,3) * t764 - t757 * t772 + t761 * t800 + t798;
t631 = t783 * t636 + t779 * t637;
t746 = -t776 * t759 - t822;
t758 = -mrSges(3,2) * t772 + mrSges(3,3) * t800;
t630 = m(3) * t746 + t764 * mrSges(3,1) + t763 * mrSges(3,2) + (t757 * t780 - t758 * t784) * t805 + t631;
t666 = t692 * pkin(4) + t787;
t811 = -t778 * t653 + t782 * t654;
t644 = m(6) * t666 - t692 * mrSges(6,2) - t693 * mrSges(6,3) - t737 * t715 - t738 * t716 + t811;
t789 = m(5) * t672 + t692 * mrSges(5,1) + t693 * mrSges(5,2) + t737 * t717 + t738 * t718 + t644;
t788 = -m(4) * t709 + t728 * mrSges(4,1) - t729 * mrSges(4,2) + t751 * t740 - t752 * t742 - t789;
t643 = m(3) * t730 + t771 * mrSges(3,1) - t763 * mrSges(3,3) + t772 * t758 - t761 * t801 + t788;
t619 = t628 * t814 - t630 * t776 + t643 * t813;
t617 = m(2) * t768 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t786 + t619;
t623 = t784 * t628 - t643 * t780;
t622 = m(2) * t769 - mrSges(2,1) * t786 - qJDD(1) * mrSges(2,2) + t623;
t812 = t785 * t617 + t781 * t622;
t809 = t737 * t819 - t738 * t820 + t766 * t826;
t808 = t737 * t825 + t738 * t821 - t766 * t819;
t807 = t737 * t821 - t738 * t827 + t766 * t820;
t618 = t628 * t816 + t777 * t630 + t643 * t815;
t799 = -t617 * t781 + t785 * t622;
t678 = Ifges(7,5) * t714 + Ifges(7,6) * t713 + Ifges(7,3) * t735;
t680 = Ifges(7,1) * t714 + Ifges(7,4) * t713 + Ifges(7,5) * t735;
t648 = -mrSges(7,1) * t659 + mrSges(7,3) * t656 + Ifges(7,4) * t675 + Ifges(7,2) * t674 + Ifges(7,6) * t691 - t678 * t714 + t680 * t735;
t679 = Ifges(7,4) * t714 + Ifges(7,2) * t713 + Ifges(7,6) * t735;
t649 = mrSges(7,2) * t659 - mrSges(7,3) * t655 + Ifges(7,1) * t675 + Ifges(7,4) * t674 + Ifges(7,5) * t691 + t678 * t713 - t679 * t735;
t624 = -mrSges(5,1) * t672 - mrSges(6,1) * t661 + mrSges(6,2) * t666 + mrSges(5,3) * t664 - pkin(4) * t644 - pkin(5) * t792 - pkin(10) * t811 - t782 * t648 - t778 * t649 + t692 * t825 + t821 * t693 + t809 * t738 + t819 * t756 + t807 * t766;
t632 = mrSges(6,1) * t662 + mrSges(7,1) * t655 + mrSges(5,2) * t672 - mrSges(7,2) * t656 - mrSges(5,3) * t663 - mrSges(6,3) * t666 + Ifges(7,5) * t675 + Ifges(7,6) * t674 + Ifges(7,3) * t691 + pkin(5) * t645 - qJ(5) * t644 + t714 * t679 - t713 * t680 + t808 * t766 + t820 * t756 + t809 * t737 + t827 * t693 - t821 * t692;
t722 = Ifges(4,5) * t752 + Ifges(4,6) * t751 - Ifges(4,3) * t766;
t724 = Ifges(4,1) * t752 + Ifges(4,4) * t751 - Ifges(4,5) * t766;
t613 = -mrSges(4,1) * t709 + mrSges(4,3) * t677 + Ifges(4,4) * t729 + Ifges(4,2) * t728 + Ifges(4,6) * t756 - pkin(3) * t789 + qJ(4) * t797 + t624 * t818 + t775 * t632 - t752 * t722 - t766 * t724;
t723 = Ifges(4,4) * t752 + Ifges(4,2) * t751 - Ifges(4,6) * t766;
t615 = mrSges(4,2) * t709 - mrSges(4,3) * t676 + Ifges(4,1) * t729 + Ifges(4,4) * t728 + Ifges(4,5) * t756 - qJ(4) * t638 - t775 * t624 + t632 * t818 + t751 * t722 + t766 * t723;
t743 = Ifges(3,3) * t772 + (Ifges(3,5) * t780 + Ifges(3,6) * t784) * t805;
t744 = Ifges(3,6) * t772 + (Ifges(3,4) * t780 + Ifges(3,2) * t784) * t805;
t612 = mrSges(3,2) * t746 - mrSges(3,3) * t730 + Ifges(3,1) * t763 - Ifges(3,4) * t764 + Ifges(3,5) * t771 - pkin(9) * t631 - t613 * t779 + t615 * t783 + t743 * t800 - t744 * t772;
t745 = Ifges(3,5) * t772 + (Ifges(3,1) * t780 + Ifges(3,4) * t784) * t805;
t614 = -t743 * t801 + (mrSges(6,1) * qJ(5) + t819) * t692 - t820 * t693 - t782 * t649 - Ifges(4,5) * t729 + mrSges(3,3) * t731 - pkin(2) * t631 - t752 * t723 - qJ(5) * t790 - pkin(4) * (t766 * t715 + t791) + t778 * t648 - mrSges(3,1) * t746 + t751 * t724 - pkin(3) * t638 + Ifges(3,4) * t763 - Ifges(3,2) * t764 - Ifges(4,6) * t728 + (mrSges(6,2) * pkin(4) - Ifges(4,3) - t826) * t756 - mrSges(4,1) * t676 + mrSges(4,2) * t677 - mrSges(5,1) * t663 + mrSges(5,2) * t664 + mrSges(6,3) * t661 - mrSges(6,2) * t662 + Ifges(3,6) * t771 + t772 * t745 + pkin(10) * t645 + (qJ(5) * t705 + t807) * t737 - t808 * t738;
t793 = pkin(8) * t623 + t612 * t780 + t614 * t784;
t611 = Ifges(3,5) * t763 - Ifges(3,6) * t764 + Ifges(3,3) * t771 + mrSges(3,1) * t730 - mrSges(3,2) * t731 + t779 * t615 + t783 * t613 + pkin(2) * t788 + pkin(9) * t798 + (t744 * t780 - t745 * t784) * t805;
t610 = -mrSges(2,2) * g(3) - mrSges(2,3) * t768 + Ifges(2,5) * qJDD(1) - t786 * Ifges(2,6) + t784 * t612 - t780 * t614 + (-t618 * t776 - t619 * t777) * pkin(8);
t609 = mrSges(2,1) * g(3) + mrSges(2,3) * t769 + t786 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t618 - t776 * t611 + t777 * t793;
t1 = [-m(1) * g(1) + t799; -m(1) * g(2) + t812; (-m(1) - m(2)) * g(3) + t618; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t812 - t781 * t609 + t785 * t610; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t799 + t785 * t609 + t781 * t610; -mrSges(1,1) * g(2) + mrSges(2,1) * t768 + mrSges(1,2) * g(1) - mrSges(2,2) * t769 + Ifges(2,3) * qJDD(1) + pkin(1) * t619 + t777 * t611 + t776 * t793;];
tauB  = t1;
