% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 10:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 10:08:12
% EndTime: 2019-05-08 10:09:42
% DurationCPUTime: 65.99s
% Computational Cost: add. (1092361->399), mult. (2366761->515), div. (0->0), fcn. (1953529->14), ass. (0->164)
t768 = cos(pkin(6));
t803 = t768 * g(3);
t767 = sin(pkin(6));
t773 = sin(qJ(2));
t802 = t767 * t773;
t779 = cos(qJ(2));
t801 = t767 * t779;
t800 = t768 * t773;
t799 = t768 * t779;
t774 = sin(qJ(1));
t780 = cos(qJ(1));
t759 = t774 * g(1) - g(2) * t780;
t781 = qJD(1) ^ 2;
t749 = pkin(8) * t767 * t781 + qJDD(1) * pkin(1) + t759;
t760 = -g(1) * t780 - g(2) * t774;
t794 = qJDD(1) * t767;
t750 = -pkin(1) * t781 + pkin(8) * t794 + t760;
t797 = t749 * t800 + t779 * t750;
t725 = -g(3) * t802 + t797;
t764 = qJD(1) * t768 + qJD(2);
t796 = qJD(1) * t767;
t793 = t773 * t796;
t747 = mrSges(3,1) * t764 - mrSges(3,3) * t793;
t751 = (-mrSges(3,1) * t779 + mrSges(3,2) * t773) * t796;
t754 = -qJD(2) * t793 + t779 * t794;
t763 = qJDD(1) * t768 + qJDD(2);
t752 = (-pkin(2) * t779 - pkin(9) * t773) * t796;
t762 = t764 ^ 2;
t795 = qJD(1) * t779;
t711 = -t762 * pkin(2) + t763 * pkin(9) + (-g(3) * t773 + t752 * t795) * t767 + t797;
t753 = (qJD(2) * t795 + qJDD(1) * t773) * t767;
t712 = -t754 * pkin(2) - t753 * pkin(9) - t803 + (-t749 + (pkin(2) * t773 - pkin(9) * t779) * t764 * qJD(1)) * t767;
t772 = sin(qJ(3));
t778 = cos(qJ(3));
t689 = -t772 * t711 + t778 * t712;
t741 = t764 * t778 - t772 * t793;
t723 = qJD(3) * t741 + t753 * t778 + t763 * t772;
t742 = t764 * t772 + t778 * t793;
t746 = qJDD(3) - t754;
t792 = t767 * t795;
t758 = qJD(3) - t792;
t675 = (t741 * t758 - t723) * pkin(10) + (t741 * t742 + t746) * pkin(3) + t689;
t690 = t778 * t711 + t772 * t712;
t722 = -qJD(3) * t742 - t753 * t772 + t763 * t778;
t732 = pkin(3) * t758 - pkin(10) * t742;
t740 = t741 ^ 2;
t677 = -pkin(3) * t740 + pkin(10) * t722 - t732 * t758 + t690;
t771 = sin(qJ(4));
t777 = cos(qJ(4));
t657 = t777 * t675 - t771 * t677;
t727 = t741 * t777 - t742 * t771;
t693 = qJD(4) * t727 + t722 * t771 + t723 * t777;
t728 = t741 * t771 + t742 * t777;
t745 = qJDD(4) + t746;
t757 = qJD(4) + t758;
t654 = (t727 * t757 - t693) * pkin(11) + (t727 * t728 + t745) * pkin(4) + t657;
t658 = t771 * t675 + t777 * t677;
t692 = -qJD(4) * t728 + t722 * t777 - t723 * t771;
t715 = pkin(4) * t757 - pkin(11) * t728;
t726 = t727 ^ 2;
t656 = -pkin(4) * t726 + pkin(11) * t692 - t715 * t757 + t658;
t770 = sin(qJ(5));
t776 = cos(qJ(5));
t651 = t770 * t654 + t776 * t656;
t705 = t727 * t770 + t728 * t776;
t667 = -qJD(5) * t705 + t692 * t776 - t693 * t770;
t704 = t727 * t776 - t728 * t770;
t687 = -mrSges(6,1) * t704 + mrSges(6,2) * t705;
t756 = qJD(5) + t757;
t697 = mrSges(6,1) * t756 - mrSges(6,3) * t705;
t739 = qJDD(5) + t745;
t688 = -pkin(5) * t704 - pkin(12) * t705;
t755 = t756 ^ 2;
t649 = -pkin(5) * t755 + pkin(12) * t739 + t688 * t704 + t651;
t724 = -g(3) * t801 + t749 * t799 - t773 * t750;
t710 = -pkin(2) * t763 - pkin(9) * t762 + t752 * t793 - t724;
t681 = -pkin(3) * t722 - pkin(10) * t740 + t742 * t732 + t710;
t663 = -t692 * pkin(4) - pkin(11) * t726 + t728 * t715 + t681;
t668 = qJD(5) * t704 + t692 * t770 + t693 * t776;
t652 = (-t704 * t756 - t668) * pkin(12) + (t705 * t756 - t667) * pkin(5) + t663;
t769 = sin(qJ(6));
t775 = cos(qJ(6));
t646 = -t649 * t769 + t652 * t775;
t694 = -t705 * t769 + t756 * t775;
t661 = qJD(6) * t694 + t668 * t775 + t739 * t769;
t666 = qJDD(6) - t667;
t695 = t705 * t775 + t756 * t769;
t678 = -mrSges(7,1) * t694 + mrSges(7,2) * t695;
t703 = qJD(6) - t704;
t679 = -mrSges(7,2) * t703 + mrSges(7,3) * t694;
t644 = m(7) * t646 + mrSges(7,1) * t666 - mrSges(7,3) * t661 - t678 * t695 + t679 * t703;
t647 = t649 * t775 + t652 * t769;
t660 = -qJD(6) * t695 - t668 * t769 + t739 * t775;
t680 = mrSges(7,1) * t703 - mrSges(7,3) * t695;
t645 = m(7) * t647 - mrSges(7,2) * t666 + mrSges(7,3) * t660 + t678 * t694 - t680 * t703;
t787 = -t644 * t769 + t775 * t645;
t635 = m(6) * t651 - mrSges(6,2) * t739 + mrSges(6,3) * t667 + t687 * t704 - t697 * t756 + t787;
t650 = t654 * t776 - t656 * t770;
t696 = -mrSges(6,2) * t756 + mrSges(6,3) * t704;
t648 = -pkin(5) * t739 - pkin(12) * t755 + t688 * t705 - t650;
t784 = -m(7) * t648 + t660 * mrSges(7,1) - mrSges(7,2) * t661 + t694 * t679 - t680 * t695;
t640 = m(6) * t650 + mrSges(6,1) * t739 - mrSges(6,3) * t668 - t687 * t705 + t696 * t756 + t784;
t629 = t770 * t635 + t776 * t640;
t706 = -mrSges(5,1) * t727 + mrSges(5,2) * t728;
t713 = -mrSges(5,2) * t757 + mrSges(5,3) * t727;
t627 = m(5) * t657 + mrSges(5,1) * t745 - mrSges(5,3) * t693 - t706 * t728 + t713 * t757 + t629;
t714 = mrSges(5,1) * t757 - mrSges(5,3) * t728;
t788 = t776 * t635 - t640 * t770;
t628 = m(5) * t658 - mrSges(5,2) * t745 + mrSges(5,3) * t692 + t706 * t727 - t714 * t757 + t788;
t621 = t777 * t627 + t771 * t628;
t729 = -mrSges(4,1) * t741 + mrSges(4,2) * t742;
t730 = -mrSges(4,2) * t758 + mrSges(4,3) * t741;
t619 = m(4) * t689 + mrSges(4,1) * t746 - mrSges(4,3) * t723 - t729 * t742 + t730 * t758 + t621;
t731 = mrSges(4,1) * t758 - mrSges(4,3) * t742;
t789 = -t627 * t771 + t777 * t628;
t620 = m(4) * t690 - mrSges(4,2) * t746 + mrSges(4,3) * t722 + t729 * t741 - t731 * t758 + t789;
t790 = -t619 * t772 + t778 * t620;
t611 = m(3) * t725 - mrSges(3,2) * t763 + mrSges(3,3) * t754 - t747 * t764 + t751 * t792 + t790;
t614 = t778 * t619 + t772 * t620;
t736 = -t767 * t749 - t803;
t748 = -mrSges(3,2) * t764 + mrSges(3,3) * t792;
t613 = m(3) * t736 - t754 * mrSges(3,1) + t753 * mrSges(3,2) + (t747 * t773 - t748 * t779) * t796 + t614;
t636 = t775 * t644 + t769 * t645;
t786 = m(6) * t663 - t667 * mrSges(6,1) + t668 * mrSges(6,2) - t704 * t696 + t705 * t697 + t636;
t783 = m(5) * t681 - t692 * mrSges(5,1) + t693 * mrSges(5,2) - t727 * t713 + t728 * t714 + t786;
t782 = -m(4) * t710 + t722 * mrSges(4,1) - t723 * mrSges(4,2) + t741 * t730 - t742 * t731 - t783;
t632 = m(3) * t724 + t763 * mrSges(3,1) - t753 * mrSges(3,3) + t764 * t748 - t751 * t793 + t782;
t602 = t611 * t800 - t613 * t767 + t632 * t799;
t600 = m(2) * t759 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t781 + t602;
t606 = t779 * t611 - t632 * t773;
t605 = m(2) * t760 - mrSges(2,1) * t781 - qJDD(1) * mrSges(2,2) + t606;
t798 = t780 * t600 + t774 * t605;
t601 = t611 * t802 + t768 * t613 + t632 * t801;
t791 = -t600 * t774 + t780 * t605;
t669 = Ifges(7,5) * t695 + Ifges(7,6) * t694 + Ifges(7,3) * t703;
t671 = Ifges(7,1) * t695 + Ifges(7,4) * t694 + Ifges(7,5) * t703;
t637 = -mrSges(7,1) * t648 + mrSges(7,3) * t647 + Ifges(7,4) * t661 + Ifges(7,2) * t660 + Ifges(7,6) * t666 - t669 * t695 + t671 * t703;
t670 = Ifges(7,4) * t695 + Ifges(7,2) * t694 + Ifges(7,6) * t703;
t638 = mrSges(7,2) * t648 - mrSges(7,3) * t646 + Ifges(7,1) * t661 + Ifges(7,4) * t660 + Ifges(7,5) * t666 + t669 * t694 - t670 * t703;
t682 = Ifges(6,5) * t705 + Ifges(6,6) * t704 + Ifges(6,3) * t756;
t683 = Ifges(6,4) * t705 + Ifges(6,2) * t704 + Ifges(6,6) * t756;
t622 = mrSges(6,2) * t663 - mrSges(6,3) * t650 + Ifges(6,1) * t668 + Ifges(6,4) * t667 + Ifges(6,5) * t739 - pkin(12) * t636 - t637 * t769 + t638 * t775 + t682 * t704 - t683 * t756;
t684 = Ifges(6,1) * t705 + Ifges(6,4) * t704 + Ifges(6,5) * t756;
t623 = -mrSges(6,1) * t663 - mrSges(7,1) * t646 + mrSges(7,2) * t647 + mrSges(6,3) * t651 + Ifges(6,4) * t668 - Ifges(7,5) * t661 + Ifges(6,2) * t667 + Ifges(6,6) * t739 - Ifges(7,6) * t660 - Ifges(7,3) * t666 - pkin(5) * t636 - t670 * t695 + t671 * t694 - t682 * t705 + t684 * t756;
t698 = Ifges(5,5) * t728 + Ifges(5,6) * t727 + Ifges(5,3) * t757;
t700 = Ifges(5,1) * t728 + Ifges(5,4) * t727 + Ifges(5,5) * t757;
t607 = -mrSges(5,1) * t681 + mrSges(5,3) * t658 + Ifges(5,4) * t693 + Ifges(5,2) * t692 + Ifges(5,6) * t745 - pkin(4) * t786 + pkin(11) * t788 + t770 * t622 + t776 * t623 - t728 * t698 + t757 * t700;
t699 = Ifges(5,4) * t728 + Ifges(5,2) * t727 + Ifges(5,6) * t757;
t615 = mrSges(5,2) * t681 - mrSges(5,3) * t657 + Ifges(5,1) * t693 + Ifges(5,4) * t692 + Ifges(5,5) * t745 - pkin(11) * t629 + t622 * t776 - t623 * t770 + t698 * t727 - t699 * t757;
t716 = Ifges(4,5) * t742 + Ifges(4,6) * t741 + Ifges(4,3) * t758;
t718 = Ifges(4,1) * t742 + Ifges(4,4) * t741 + Ifges(4,5) * t758;
t596 = -mrSges(4,1) * t710 + mrSges(4,3) * t690 + Ifges(4,4) * t723 + Ifges(4,2) * t722 + Ifges(4,6) * t746 - pkin(3) * t783 + pkin(10) * t789 + t777 * t607 + t771 * t615 - t742 * t716 + t758 * t718;
t717 = Ifges(4,4) * t742 + Ifges(4,2) * t741 + Ifges(4,6) * t758;
t597 = mrSges(4,2) * t710 - mrSges(4,3) * t689 + Ifges(4,1) * t723 + Ifges(4,4) * t722 + Ifges(4,5) * t746 - pkin(10) * t621 - t607 * t771 + t615 * t777 + t716 * t741 - t717 * t758;
t733 = Ifges(3,3) * t764 + (Ifges(3,5) * t773 + Ifges(3,6) * t779) * t796;
t734 = Ifges(3,6) * t764 + (Ifges(3,4) * t773 + Ifges(3,2) * t779) * t796;
t595 = mrSges(3,2) * t736 - mrSges(3,3) * t724 + Ifges(3,1) * t753 + Ifges(3,4) * t754 + Ifges(3,5) * t763 - pkin(9) * t614 - t596 * t772 + t597 * t778 + t733 * t792 - t734 * t764;
t735 = Ifges(3,5) * t764 + (Ifges(3,1) * t773 + Ifges(3,4) * t779) * t796;
t598 = -pkin(5) * t784 - pkin(2) * t614 - t775 * t637 + Ifges(3,6) * t763 + t764 * t735 - t769 * t638 - Ifges(5,3) * t745 - Ifges(4,3) * t746 + Ifges(3,4) * t753 + Ifges(3,2) * t754 - Ifges(6,3) * t739 + t741 * t718 - t742 * t717 + t727 * t700 - t728 * t699 - mrSges(3,1) * t736 - Ifges(4,6) * t722 - Ifges(4,5) * t723 + mrSges(3,3) * t725 - t705 * t683 + t704 * t684 - Ifges(5,6) * t692 - Ifges(5,5) * t693 - mrSges(4,1) * t689 + mrSges(4,2) * t690 - Ifges(6,6) * t667 - Ifges(6,5) * t668 - mrSges(5,1) * t657 + mrSges(5,2) * t658 + mrSges(6,2) * t651 - mrSges(6,1) * t650 - pkin(4) * t629 - pkin(12) * t787 - t733 * t793 - pkin(3) * t621;
t785 = pkin(8) * t606 + t595 * t773 + t598 * t779;
t594 = Ifges(3,5) * t753 + Ifges(3,6) * t754 + Ifges(3,3) * t763 + mrSges(3,1) * t724 - mrSges(3,2) * t725 + t772 * t597 + t778 * t596 + pkin(2) * t782 + pkin(9) * t790 + (t734 * t773 - t735 * t779) * t796;
t593 = -mrSges(2,2) * g(3) - mrSges(2,3) * t759 + Ifges(2,5) * qJDD(1) - t781 * Ifges(2,6) + t779 * t595 - t773 * t598 + (-t601 * t767 - t602 * t768) * pkin(8);
t592 = mrSges(2,1) * g(3) + mrSges(2,3) * t760 + t781 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t601 - t767 * t594 + t768 * t785;
t1 = [-m(1) * g(1) + t791; -m(1) * g(2) + t798; (-m(1) - m(2)) * g(3) + t601; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t798 - t774 * t592 + t780 * t593; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t791 + t780 * t592 + t774 * t593; -mrSges(1,1) * g(2) + mrSges(2,1) * t759 + mrSges(1,2) * g(1) - mrSges(2,2) * t760 + Ifges(2,3) * qJDD(1) + pkin(1) * t602 + t768 * t594 + t767 * t785;];
tauB  = t1;
