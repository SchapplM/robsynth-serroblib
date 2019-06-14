% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 07:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:59:45
% EndTime: 2019-05-07 07:00:08
% DurationCPUTime: 19.65s
% Computational Cost: add. (310741->379), mult. (678567->474), div. (0->0), fcn. (525926->12), ass. (0->155)
t816 = -Ifges(4,1) - Ifges(5,2);
t810 = Ifges(4,5) - Ifges(5,4);
t815 = -Ifges(4,2) - Ifges(5,3);
t809 = Ifges(4,6) - Ifges(5,5);
t808 = -Ifges(5,6) - Ifges(4,4);
t814 = Ifges(4,3) + Ifges(5,1);
t767 = sin(pkin(6));
t772 = sin(qJ(2));
t775 = cos(qJ(2));
t793 = qJD(1) * qJD(2);
t752 = (-qJDD(1) * t775 + t772 * t793) * t767;
t796 = qJD(1) * t767;
t750 = (-pkin(2) * t775 - pkin(9) * t772) * t796;
t769 = cos(pkin(6));
t763 = t769 * qJD(1) + qJD(2);
t761 = t763 ^ 2;
t762 = t769 * qJDD(1) + qJDD(2);
t795 = qJD(1) * t775;
t773 = sin(qJ(1));
t776 = cos(qJ(1));
t759 = t773 * g(1) - t776 * g(2);
t777 = qJD(1) ^ 2;
t812 = pkin(8) * t767;
t747 = qJDD(1) * pkin(1) + t777 * t812 + t759;
t760 = -t776 * g(1) - t773 * g(2);
t748 = -t777 * pkin(1) + qJDD(1) * t812 + t760;
t804 = t769 * t772;
t797 = t747 * t804 + t775 * t748;
t683 = -t761 * pkin(2) + t762 * pkin(9) + (-g(3) * t772 + t750 * t795) * t767 + t797;
t751 = (qJDD(1) * t772 + t775 * t793) * t767;
t811 = t769 * g(3);
t684 = t752 * pkin(2) - t751 * pkin(9) - t811 + (-t747 + (pkin(2) * t772 - pkin(9) * t775) * t763 * qJD(1)) * t767;
t771 = sin(qJ(3));
t813 = cos(qJ(3));
t667 = t813 * t683 + t771 * t684;
t792 = t772 * t796;
t738 = -t813 * t763 + t771 * t792;
t739 = t771 * t763 + t813 * t792;
t713 = t738 * pkin(3) - t739 * qJ(4);
t744 = qJDD(3) + t752;
t791 = t767 * t795;
t757 = -qJD(3) + t791;
t756 = t757 ^ 2;
t660 = t756 * pkin(3) - t744 * qJ(4) + 0.2e1 * qJD(4) * t757 + t738 * t713 - t667;
t807 = t738 * t757;
t806 = t767 * t772;
t805 = t767 * t775;
t803 = t769 * t775;
t712 = -g(3) * t806 + t797;
t745 = t763 * mrSges(3,1) - mrSges(3,3) * t792;
t749 = (-mrSges(3,1) * t775 + mrSges(3,2) * t772) * t796;
t666 = -t771 * t683 + t813 * t684;
t710 = -t738 * qJD(3) + t813 * t751 + t771 * t762;
t714 = t738 * mrSges(4,1) + t739 * mrSges(4,2);
t723 = t757 * mrSges(4,2) - t738 * mrSges(4,3);
t726 = t738 * mrSges(5,1) + t757 * mrSges(5,3);
t661 = -t744 * pkin(3) - t756 * qJ(4) + t739 * t713 + qJDD(4) - t666;
t655 = (t738 * t739 - t744) * qJ(5) + (t710 - t807) * pkin(4) + t661;
t709 = t739 * qJD(3) + t771 * t751 - t813 * t762;
t725 = t739 * pkin(4) + t757 * qJ(5);
t737 = t738 ^ 2;
t711 = -g(3) * t805 + t747 * t803 - t772 * t748;
t682 = -t762 * pkin(2) - t761 * pkin(9) + t750 * t792 - t711;
t779 = (-t710 - t807) * qJ(4) + t682 + (-t757 * pkin(3) - 0.2e1 * qJD(4)) * t739;
t659 = -t737 * pkin(4) - t739 * t725 + (pkin(3) + qJ(5)) * t709 + t779;
t766 = sin(pkin(11));
t768 = cos(pkin(11));
t722 = t766 * t738 - t768 * t757;
t649 = -0.2e1 * qJD(5) * t722 + t768 * t655 - t766 * t659;
t689 = t766 * t709 + t768 * t744;
t721 = t768 * t738 + t766 * t757;
t647 = (t721 * t739 - t689) * pkin(10) + (t721 * t722 + t710) * pkin(5) + t649;
t650 = 0.2e1 * qJD(5) * t721 + t766 * t655 + t768 * t659;
t688 = t768 * t709 - t766 * t744;
t696 = t739 * pkin(5) - t722 * pkin(10);
t720 = t721 ^ 2;
t648 = -t720 * pkin(5) + t688 * pkin(10) - t739 * t696 + t650;
t770 = sin(qJ(6));
t774 = cos(qJ(6));
t645 = t774 * t647 - t770 * t648;
t685 = t774 * t721 - t770 * t722;
t665 = t685 * qJD(6) + t770 * t688 + t774 * t689;
t686 = t770 * t721 + t774 * t722;
t672 = -t685 * mrSges(7,1) + t686 * mrSges(7,2);
t736 = qJD(6) + t739;
t673 = -t736 * mrSges(7,2) + t685 * mrSges(7,3);
t706 = qJDD(6) + t710;
t643 = m(7) * t645 + t706 * mrSges(7,1) - t665 * mrSges(7,3) - t686 * t672 + t736 * t673;
t646 = t770 * t647 + t774 * t648;
t664 = -t686 * qJD(6) + t774 * t688 - t770 * t689;
t674 = t736 * mrSges(7,1) - t686 * mrSges(7,3);
t644 = m(7) * t646 - t706 * mrSges(7,2) + t664 * mrSges(7,3) + t685 * t672 - t736 * t674;
t634 = t774 * t643 + t770 * t644;
t690 = -t721 * mrSges(6,1) + t722 * mrSges(6,2);
t694 = -t739 * mrSges(6,2) + t721 * mrSges(6,3);
t632 = m(6) * t649 + t710 * mrSges(6,1) - t689 * mrSges(6,3) - t722 * t690 + t739 * t694 + t634;
t695 = t739 * mrSges(6,1) - t722 * mrSges(6,3);
t788 = -t770 * t643 + t774 * t644;
t633 = m(6) * t650 - t710 * mrSges(6,2) + t688 * mrSges(6,3) + t721 * t690 - t739 * t695 + t788;
t629 = t768 * t632 + t766 * t633;
t715 = -t738 * mrSges(5,2) - t739 * mrSges(5,3);
t782 = -m(5) * t661 - t710 * mrSges(5,1) - t739 * t715 - t629;
t627 = m(4) * t666 - t710 * mrSges(4,3) - t739 * t714 + (-t723 + t726) * t757 + (mrSges(4,1) - mrSges(5,2)) * t744 + t782;
t724 = -t757 * mrSges(4,1) - t739 * mrSges(4,3);
t727 = t739 * mrSges(5,1) - t757 * mrSges(5,2);
t658 = -t709 * pkin(4) - t737 * qJ(5) - t757 * t725 + qJDD(5) - t660;
t652 = -t688 * pkin(5) - t720 * pkin(10) + t722 * t696 + t658;
t783 = m(7) * t652 - t664 * mrSges(7,1) + t665 * mrSges(7,2) - t685 * t673 + t686 * t674;
t781 = -m(6) * t658 + t688 * mrSges(6,1) - t689 * mrSges(6,2) + t721 * t694 - t722 * t695 - t783;
t778 = -m(5) * t660 + t744 * mrSges(5,3) - t757 * t727 - t781;
t639 = m(4) * t667 + t757 * t724 + t778 + (-t714 - t715) * t738 + (-mrSges(4,3) - mrSges(5,1)) * t709 - t744 * mrSges(4,2);
t789 = -t771 * t627 + t813 * t639;
t617 = m(3) * t712 - t762 * mrSges(3,2) - t752 * mrSges(3,3) - t763 * t745 + t749 * t791 + t789;
t620 = t813 * t627 + t771 * t639;
t732 = -t767 * t747 - t811;
t746 = -t763 * mrSges(3,2) + mrSges(3,3) * t791;
t619 = m(3) * t732 + t752 * mrSges(3,1) + t751 * mrSges(3,2) + (t745 * t772 - t746 * t775) * t796 + t620;
t662 = t709 * pkin(3) + t779;
t801 = -t766 * t632 + t768 * t633;
t787 = -m(5) * t662 + t709 * mrSges(5,2) + t738 * t726 - t801;
t780 = -m(4) * t682 - t709 * mrSges(4,1) - t738 * t723 + (-t724 + t727) * t739 + (-mrSges(4,2) + mrSges(5,3)) * t710 + t787;
t625 = m(3) * t711 + t762 * mrSges(3,1) - t751 * mrSges(3,3) + t763 * t746 - t749 * t792 + t780;
t608 = t617 * t804 - t767 * t619 + t625 * t803;
t606 = m(2) * t759 + qJDD(1) * mrSges(2,1) - t777 * mrSges(2,2) + t608;
t613 = t775 * t617 - t772 * t625;
t612 = m(2) * t760 - t777 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t613;
t802 = t776 * t606 + t773 * t612;
t800 = t738 * t809 - t739 * t810 + t757 * t814;
t799 = t738 * t815 - t739 * t808 - t757 * t809;
t798 = -t808 * t738 + t739 * t816 + t810 * t757;
t607 = t617 * t806 + t769 * t619 + t625 * t805;
t790 = -t773 * t606 + t776 * t612;
t668 = Ifges(7,5) * t686 + Ifges(7,6) * t685 + Ifges(7,3) * t736;
t670 = Ifges(7,1) * t686 + Ifges(7,4) * t685 + Ifges(7,5) * t736;
t635 = -mrSges(7,1) * t652 + mrSges(7,3) * t646 + Ifges(7,4) * t665 + Ifges(7,2) * t664 + Ifges(7,6) * t706 - t686 * t668 + t736 * t670;
t669 = Ifges(7,4) * t686 + Ifges(7,2) * t685 + Ifges(7,6) * t736;
t636 = mrSges(7,2) * t652 - mrSges(7,3) * t645 + Ifges(7,1) * t665 + Ifges(7,4) * t664 + Ifges(7,5) * t706 + t685 * t668 - t736 * t669;
t675 = Ifges(6,5) * t722 + Ifges(6,6) * t721 + Ifges(6,3) * t739;
t677 = Ifges(6,1) * t722 + Ifges(6,4) * t721 + Ifges(6,5) * t739;
t621 = -mrSges(6,1) * t658 + mrSges(6,3) * t650 + Ifges(6,4) * t689 + Ifges(6,2) * t688 + Ifges(6,6) * t710 - pkin(5) * t783 + pkin(10) * t788 + t774 * t635 + t770 * t636 - t722 * t675 + t739 * t677;
t676 = Ifges(6,4) * t722 + Ifges(6,2) * t721 + Ifges(6,6) * t739;
t622 = mrSges(6,2) * t658 - mrSges(6,3) * t649 + Ifges(6,1) * t689 + Ifges(6,4) * t688 + Ifges(6,5) * t710 - pkin(10) * t634 - t770 * t635 + t774 * t636 + t721 * t675 - t739 * t676;
t628 = -t710 * mrSges(5,3) - t739 * t727 - t787;
t604 = -mrSges(4,1) * t682 - mrSges(5,1) * t660 + mrSges(5,2) * t662 + mrSges(4,3) * t667 - pkin(3) * t628 - pkin(4) * t781 - qJ(5) * t801 - t768 * t621 - t766 * t622 + t709 * t815 - t808 * t710 + t800 * t739 + t809 * t744 + t798 * t757;
t609 = (Ifges(6,3) - t816) * t710 - t685 * t670 + t686 * t669 + mrSges(4,2) * t682 + pkin(5) * t634 + pkin(4) * t629 + Ifges(7,5) * t665 - mrSges(4,3) * t666 - mrSges(7,2) * t646 - qJ(4) * t628 + Ifges(7,3) * t706 + t808 * t709 + t810 * t744 + t799 * t757 + t800 * t738 + mrSges(7,1) * t645 + mrSges(5,1) * t661 - mrSges(5,3) * t662 + Ifges(7,6) * t664 + mrSges(6,1) * t649 - mrSges(6,2) * t650 - t721 * t677 + t722 * t676 + Ifges(6,6) * t688 + Ifges(6,5) * t689;
t729 = Ifges(3,3) * t763 + (Ifges(3,5) * t772 + Ifges(3,6) * t775) * t796;
t730 = Ifges(3,6) * t763 + (Ifges(3,4) * t772 + Ifges(3,2) * t775) * t796;
t602 = mrSges(3,2) * t732 - mrSges(3,3) * t711 + Ifges(3,1) * t751 - Ifges(3,4) * t752 + Ifges(3,5) * t762 - pkin(9) * t620 - t771 * t604 + t813 * t609 + t729 * t791 - t763 * t730;
t731 = Ifges(3,5) * t763 + (Ifges(3,1) * t772 + Ifges(3,4) * t775) * t796;
t603 = Ifges(3,6) * t762 + t763 * t731 + t766 * t621 - pkin(3) * (t757 * t726 + t782) - t768 * t622 + qJ(5) * t629 - t729 * t792 - mrSges(4,1) * t666 + mrSges(4,2) * t667 - mrSges(3,1) * t732 - pkin(2) * t620 + mrSges(5,3) * t660 - mrSges(5,2) * t661 + mrSges(3,3) * t712 - qJ(4) * t778 + Ifges(3,4) * t751 - Ifges(3,2) * t752 + (pkin(3) * mrSges(5,2) - t814) * t744 - t799 * t739 + (qJ(4) * t715 + t798) * t738 - t810 * t710 + (qJ(4) * mrSges(5,1) + t809) * t709;
t784 = pkin(8) * t613 + t602 * t772 + t603 * t775;
t601 = Ifges(3,5) * t751 - Ifges(3,6) * t752 + Ifges(3,3) * t762 + mrSges(3,1) * t711 - mrSges(3,2) * t712 + t771 * t609 + t813 * t604 + pkin(2) * t780 + pkin(9) * t789 + (t730 * t772 - t731 * t775) * t796;
t600 = -mrSges(2,2) * g(3) - mrSges(2,3) * t759 + Ifges(2,5) * qJDD(1) - t777 * Ifges(2,6) + t775 * t602 - t772 * t603 + (-t607 * t767 - t608 * t769) * pkin(8);
t599 = mrSges(2,1) * g(3) + mrSges(2,3) * t760 + t777 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t607 - t767 * t601 + t784 * t769;
t1 = [-m(1) * g(1) + t790; -m(1) * g(2) + t802; (-m(1) - m(2)) * g(3) + t607; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t802 - t773 * t599 + t776 * t600; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t790 + t776 * t599 + t773 * t600; -mrSges(1,1) * g(2) + mrSges(2,1) * t759 + mrSges(1,2) * g(1) - mrSges(2,2) * t760 + Ifges(2,3) * qJDD(1) + pkin(1) * t608 + t769 * t601 + t784 * t767;];
tauB  = t1;
