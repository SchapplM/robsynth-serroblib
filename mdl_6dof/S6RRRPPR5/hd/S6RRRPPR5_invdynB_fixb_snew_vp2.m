% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 05:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:02:16
% EndTime: 2019-05-07 05:03:10
% DurationCPUTime: 52.15s
% Computational Cost: add. (866311->397), mult. (1930398->517), div. (0->0), fcn. (1561982->14), ass. (0->162)
t801 = -2 * qJD(4);
t762 = sin(pkin(6));
t767 = sin(qJ(2));
t771 = cos(qJ(2));
t788 = qJD(1) * qJD(2);
t750 = (-qJDD(1) * t771 + t767 * t788) * t762;
t791 = qJD(1) * t762;
t748 = (-pkin(2) * t771 - pkin(9) * t767) * t791;
t764 = cos(pkin(6));
t757 = qJD(1) * t764 + qJD(2);
t755 = t757 ^ 2;
t756 = qJDD(1) * t764 + qJDD(2);
t790 = qJD(1) * t771;
t768 = sin(qJ(1));
t772 = cos(qJ(1));
t753 = t768 * g(1) - g(2) * t772;
t773 = qJD(1) ^ 2;
t800 = pkin(8) * t762;
t745 = qJDD(1) * pkin(1) + t773 * t800 + t753;
t754 = -g(1) * t772 - g(2) * t768;
t746 = -pkin(1) * t773 + qJDD(1) * t800 + t754;
t795 = t764 * t767;
t792 = t745 * t795 + t771 * t746;
t701 = -t755 * pkin(2) + t756 * pkin(9) + (-g(3) * t767 + t748 * t790) * t762 + t792;
t749 = (qJDD(1) * t767 + t771 * t788) * t762;
t799 = t764 * g(3);
t702 = t750 * pkin(2) - t749 * pkin(9) - t799 + (-t745 + (pkin(2) * t767 - pkin(9) * t771) * t757 * qJD(1)) * t762;
t766 = sin(qJ(3));
t770 = cos(qJ(3));
t671 = -t766 * t701 + t770 * t702;
t787 = t767 * t791;
t738 = t757 * t770 - t766 * t787;
t718 = qJD(3) * t738 + t749 * t770 + t756 * t766;
t739 = t757 * t766 + t770 * t787;
t742 = qJDD(3) + t750;
t786 = t762 * t790;
t752 = qJD(3) - t786;
t659 = (t738 * t752 - t718) * qJ(4) + (t738 * t739 + t742) * pkin(3) + t671;
t672 = t770 * t701 + t766 * t702;
t717 = -qJD(3) * t739 - t749 * t766 + t756 * t770;
t728 = pkin(3) * t752 - qJ(4) * t739;
t737 = t738 ^ 2;
t666 = -pkin(3) * t737 + qJ(4) * t717 - t728 * t752 + t672;
t761 = sin(pkin(11));
t798 = cos(pkin(11));
t725 = t761 * t738 + t739 * t798;
t650 = t659 * t798 - t761 * t666 + t725 * t801;
t797 = t762 * t767;
t796 = t762 * t771;
t794 = t764 * t771;
t720 = -g(3) * t797 + t792;
t743 = mrSges(3,1) * t757 - mrSges(3,3) * t787;
t747 = (-mrSges(3,1) * t771 + mrSges(3,2) * t767) * t791;
t724 = -t798 * t738 + t739 * t761;
t651 = t761 * t659 + t798 * t666 + t724 * t801;
t689 = -t798 * t717 + t718 * t761;
t696 = mrSges(5,1) * t724 + mrSges(5,2) * t725;
t709 = mrSges(5,1) * t752 - mrSges(5,3) * t725;
t695 = pkin(4) * t724 - qJ(5) * t725;
t751 = t752 ^ 2;
t649 = -pkin(4) * t751 + qJ(5) * t742 - t695 * t724 + t651;
t719 = -g(3) * t796 + t745 * t794 - t767 * t746;
t700 = -t756 * pkin(2) - t755 * pkin(9) + t748 * t787 - t719;
t668 = -t717 * pkin(3) - t737 * qJ(4) + t739 * t728 + qJDD(4) + t700;
t690 = t761 * t717 + t718 * t798;
t654 = (t724 * t752 - t690) * qJ(5) + (t725 * t752 + t689) * pkin(4) + t668;
t760 = sin(pkin(12));
t763 = cos(pkin(12));
t707 = t725 * t763 + t752 * t760;
t644 = -0.2e1 * qJD(5) * t707 - t760 * t649 + t763 * t654;
t680 = t690 * t763 + t742 * t760;
t706 = -t725 * t760 + t752 * t763;
t642 = (t706 * t724 - t680) * pkin(10) + (t706 * t707 + t689) * pkin(5) + t644;
t645 = 0.2e1 * qJD(5) * t706 + t763 * t649 + t760 * t654;
t679 = -t690 * t760 + t742 * t763;
t686 = pkin(5) * t724 - pkin(10) * t707;
t705 = t706 ^ 2;
t643 = -pkin(5) * t705 + pkin(10) * t679 - t686 * t724 + t645;
t765 = sin(qJ(6));
t769 = cos(qJ(6));
t640 = t642 * t769 - t643 * t765;
t681 = t706 * t769 - t707 * t765;
t657 = qJD(6) * t681 + t679 * t765 + t680 * t769;
t682 = t706 * t765 + t707 * t769;
t667 = -mrSges(7,1) * t681 + mrSges(7,2) * t682;
t723 = qJD(6) + t724;
t669 = -mrSges(7,2) * t723 + mrSges(7,3) * t681;
t688 = qJDD(6) + t689;
t638 = m(7) * t640 + mrSges(7,1) * t688 - mrSges(7,3) * t657 - t667 * t682 + t669 * t723;
t641 = t642 * t765 + t643 * t769;
t656 = -qJD(6) * t682 + t679 * t769 - t680 * t765;
t670 = mrSges(7,1) * t723 - mrSges(7,3) * t682;
t639 = m(7) * t641 - mrSges(7,2) * t688 + mrSges(7,3) * t656 + t667 * t681 - t670 * t723;
t630 = t769 * t638 + t765 * t639;
t683 = -mrSges(6,1) * t706 + mrSges(6,2) * t707;
t684 = -mrSges(6,2) * t724 + mrSges(6,3) * t706;
t628 = m(6) * t644 + mrSges(6,1) * t689 - mrSges(6,3) * t680 - t683 * t707 + t684 * t724 + t630;
t685 = mrSges(6,1) * t724 - mrSges(6,3) * t707;
t781 = -t638 * t765 + t769 * t639;
t629 = m(6) * t645 - mrSges(6,2) * t689 + mrSges(6,3) * t679 + t683 * t706 - t685 * t724 + t781;
t782 = -t628 * t760 + t763 * t629;
t623 = m(5) * t651 - mrSges(5,2) * t742 - mrSges(5,3) * t689 - t696 * t724 - t709 * t752 + t782;
t708 = -mrSges(5,2) * t752 - mrSges(5,3) * t724;
t648 = -t742 * pkin(4) - t751 * qJ(5) + t725 * t695 + qJDD(5) - t650;
t646 = -t679 * pkin(5) - t705 * pkin(10) + t707 * t686 + t648;
t777 = m(7) * t646 - t656 * mrSges(7,1) + mrSges(7,2) * t657 - t681 * t669 + t670 * t682;
t775 = -m(6) * t648 + t679 * mrSges(6,1) - mrSges(6,2) * t680 + t706 * t684 - t685 * t707 - t777;
t634 = m(5) * t650 + mrSges(5,1) * t742 - mrSges(5,3) * t690 - t696 * t725 + t708 * t752 + t775;
t615 = t761 * t623 + t798 * t634;
t726 = -mrSges(4,1) * t738 + mrSges(4,2) * t739;
t727 = -mrSges(4,2) * t752 + mrSges(4,3) * t738;
t613 = m(4) * t671 + mrSges(4,1) * t742 - mrSges(4,3) * t718 - t726 * t739 + t727 * t752 + t615;
t729 = mrSges(4,1) * t752 - mrSges(4,3) * t739;
t783 = t798 * t623 - t634 * t761;
t614 = m(4) * t672 - mrSges(4,2) * t742 + mrSges(4,3) * t717 + t726 * t738 - t729 * t752 + t783;
t784 = -t613 * t766 + t770 * t614;
t605 = m(3) * t720 - mrSges(3,2) * t756 - mrSges(3,3) * t750 - t743 * t757 + t747 * t786 + t784;
t608 = t770 * t613 + t766 * t614;
t733 = -t762 * t745 - t799;
t744 = -mrSges(3,2) * t757 + mrSges(3,3) * t786;
t607 = m(3) * t733 + t750 * mrSges(3,1) + t749 * mrSges(3,2) + (t743 * t767 - t744 * t771) * t791 + t608;
t624 = t763 * t628 + t760 * t629;
t776 = m(5) * t668 + t689 * mrSges(5,1) + mrSges(5,2) * t690 + t724 * t708 + t709 * t725 + t624;
t774 = -m(4) * t700 + t717 * mrSges(4,1) - mrSges(4,2) * t718 + t738 * t727 - t729 * t739 - t776;
t620 = m(3) * t719 + mrSges(3,1) * t756 - mrSges(3,3) * t749 + t744 * t757 - t747 * t787 + t774;
t596 = t605 * t795 - t607 * t762 + t620 * t794;
t594 = m(2) * t753 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t773 + t596;
t600 = t771 * t605 - t620 * t767;
t599 = m(2) * t754 - mrSges(2,1) * t773 - qJDD(1) * mrSges(2,2) + t600;
t793 = t772 * t594 + t768 * t599;
t595 = t605 * t797 + t764 * t607 + t620 * t796;
t785 = -t594 * t768 + t772 * t599;
t660 = Ifges(7,5) * t682 + Ifges(7,6) * t681 + Ifges(7,3) * t723;
t662 = Ifges(7,1) * t682 + Ifges(7,4) * t681 + Ifges(7,5) * t723;
t631 = -mrSges(7,1) * t646 + mrSges(7,3) * t641 + Ifges(7,4) * t657 + Ifges(7,2) * t656 + Ifges(7,6) * t688 - t660 * t682 + t662 * t723;
t661 = Ifges(7,4) * t682 + Ifges(7,2) * t681 + Ifges(7,6) * t723;
t632 = mrSges(7,2) * t646 - mrSges(7,3) * t640 + Ifges(7,1) * t657 + Ifges(7,4) * t656 + Ifges(7,5) * t688 + t660 * t681 - t661 * t723;
t673 = Ifges(6,5) * t707 + Ifges(6,6) * t706 + Ifges(6,3) * t724;
t675 = Ifges(6,1) * t707 + Ifges(6,4) * t706 + Ifges(6,5) * t724;
t616 = -mrSges(6,1) * t648 + mrSges(6,3) * t645 + Ifges(6,4) * t680 + Ifges(6,2) * t679 + Ifges(6,6) * t689 - pkin(5) * t777 + pkin(10) * t781 + t769 * t631 + t765 * t632 - t707 * t673 + t724 * t675;
t674 = Ifges(6,4) * t707 + Ifges(6,2) * t706 + Ifges(6,6) * t724;
t617 = mrSges(6,2) * t648 - mrSges(6,3) * t644 + Ifges(6,1) * t680 + Ifges(6,4) * t679 + Ifges(6,5) * t689 - pkin(10) * t630 - t631 * t765 + t632 * t769 + t673 * t706 - t674 * t724;
t691 = Ifges(5,5) * t725 - Ifges(5,6) * t724 + Ifges(5,3) * t752;
t692 = Ifges(5,4) * t725 - Ifges(5,2) * t724 + Ifges(5,6) * t752;
t601 = mrSges(5,2) * t668 - mrSges(5,3) * t650 + Ifges(5,1) * t690 - Ifges(5,4) * t689 + Ifges(5,5) * t742 - qJ(5) * t624 - t616 * t760 + t617 * t763 - t691 * t724 - t692 * t752;
t693 = Ifges(5,1) * t725 - Ifges(5,4) * t724 + Ifges(5,5) * t752;
t609 = Ifges(5,4) * t690 + Ifges(5,6) * t742 - t725 * t691 + t752 * t693 - mrSges(5,1) * t668 + mrSges(5,3) * t651 - Ifges(6,5) * t680 - Ifges(6,6) * t679 - t707 * t674 + t706 * t675 - mrSges(6,1) * t644 + mrSges(6,2) * t645 - Ifges(7,5) * t657 - Ifges(7,6) * t656 - Ifges(7,3) * t688 - t682 * t661 + t681 * t662 - mrSges(7,1) * t640 + mrSges(7,2) * t641 - pkin(5) * t630 - pkin(4) * t624 + (-Ifges(5,2) - Ifges(6,3)) * t689;
t711 = Ifges(4,5) * t739 + Ifges(4,6) * t738 + Ifges(4,3) * t752;
t713 = Ifges(4,1) * t739 + Ifges(4,4) * t738 + Ifges(4,5) * t752;
t590 = -mrSges(4,1) * t700 + mrSges(4,3) * t672 + Ifges(4,4) * t718 + Ifges(4,2) * t717 + Ifges(4,6) * t742 - pkin(3) * t776 + qJ(4) * t783 + t761 * t601 + t609 * t798 - t739 * t711 + t752 * t713;
t712 = Ifges(4,4) * t739 + Ifges(4,2) * t738 + Ifges(4,6) * t752;
t592 = mrSges(4,2) * t700 - mrSges(4,3) * t671 + Ifges(4,1) * t718 + Ifges(4,4) * t717 + Ifges(4,5) * t742 - qJ(4) * t615 + t601 * t798 - t761 * t609 + t738 * t711 - t752 * t712;
t730 = Ifges(3,3) * t757 + (Ifges(3,5) * t767 + Ifges(3,6) * t771) * t791;
t731 = Ifges(3,6) * t757 + (Ifges(3,4) * t767 + Ifges(3,2) * t771) * t791;
t589 = mrSges(3,2) * t733 - mrSges(3,3) * t719 + Ifges(3,1) * t749 - Ifges(3,4) * t750 + Ifges(3,5) * t756 - pkin(9) * t608 - t590 * t766 + t592 * t770 + t730 * t786 - t731 * t757;
t732 = Ifges(3,5) * t757 + (Ifges(3,1) * t767 + Ifges(3,4) * t771) * t791;
t591 = -Ifges(4,6) * t717 - Ifges(4,5) * t718 + Ifges(3,4) * t749 - Ifges(3,2) * t750 - t760 * t617 - t763 * t616 - mrSges(4,1) * t671 + mrSges(4,2) * t672 - pkin(4) * t775 + t738 * t713 - t739 * t712 - t724 * t693 - t725 * t692 - mrSges(3,1) * t733 + mrSges(5,2) * t651 - mrSges(5,1) * t650 + Ifges(3,6) * t756 + t757 * t732 - qJ(5) * t782 - t730 * t787 + Ifges(5,6) * t689 - Ifges(5,5) * t690 - pkin(2) * t608 - pkin(3) * t615 + mrSges(3,3) * t720 + (-Ifges(4,3) - Ifges(5,3)) * t742;
t778 = pkin(8) * t600 + t589 * t767 + t591 * t771;
t588 = Ifges(3,5) * t749 - Ifges(3,6) * t750 + Ifges(3,3) * t756 + mrSges(3,1) * t719 - mrSges(3,2) * t720 + t766 * t592 + t770 * t590 + pkin(2) * t774 + pkin(9) * t784 + (t731 * t767 - t732 * t771) * t791;
t587 = -mrSges(2,2) * g(3) - mrSges(2,3) * t753 + Ifges(2,5) * qJDD(1) - t773 * Ifges(2,6) + t771 * t589 - t767 * t591 + (-t595 * t762 - t596 * t764) * pkin(8);
t586 = mrSges(2,1) * g(3) + mrSges(2,3) * t754 + t773 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t595 - t762 * t588 + t764 * t778;
t1 = [-m(1) * g(1) + t785; -m(1) * g(2) + t793; (-m(1) - m(2)) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t793 - t768 * t586 + t772 * t587; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t785 + t772 * t586 + t768 * t587; -mrSges(1,1) * g(2) + mrSges(2,1) * t753 + mrSges(1,2) * g(1) - mrSges(2,2) * t754 + Ifges(2,3) * qJDD(1) + pkin(1) * t596 + t764 * t588 + t762 * t778;];
tauB  = t1;
