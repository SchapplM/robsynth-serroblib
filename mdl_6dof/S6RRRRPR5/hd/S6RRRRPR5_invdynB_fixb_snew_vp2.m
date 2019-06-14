% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 20:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:26:38
% EndTime: 2019-05-07 20:26:57
% DurationCPUTime: 9.86s
% Computational Cost: add. (148057->367), mult. (295666->446), div. (0->0), fcn. (210220->10), ass. (0->144)
t788 = Ifges(5,1) + Ifges(6,1);
t780 = Ifges(5,4) - Ifges(6,5);
t779 = Ifges(5,5) + Ifges(6,4);
t787 = Ifges(5,2) + Ifges(6,3);
t778 = Ifges(5,6) - Ifges(6,6);
t786 = -Ifges(5,3) - Ifges(6,2);
t746 = sin(qJ(3));
t750 = cos(qJ(3));
t751 = cos(qJ(2));
t769 = qJD(1) * t751;
t747 = sin(qJ(2));
t770 = qJD(1) * t747;
t720 = -t746 * t770 + t750 * t769;
t768 = qJD(1) * qJD(2);
t729 = qJDD(1) * t747 + t751 * t768;
t730 = qJDD(1) * t751 - t747 * t768;
t692 = qJD(3) * t720 + t729 * t750 + t730 * t746;
t721 = (t746 * t751 + t747 * t750) * qJD(1);
t742 = qJD(2) + qJD(3);
t745 = sin(qJ(4));
t783 = cos(qJ(4));
t706 = t721 * t745 - t783 * t742;
t741 = qJDD(2) + qJDD(3);
t658 = -t706 * qJD(4) + t783 * t692 + t745 * t741;
t748 = sin(qJ(1));
t752 = cos(qJ(1));
t735 = -g(1) * t752 - g(2) * t748;
t753 = qJD(1) ^ 2;
t723 = -pkin(1) * t753 + qJDD(1) * pkin(7) + t735;
t776 = t723 * t747;
t782 = pkin(2) * t753;
t682 = qJDD(2) * pkin(2) - pkin(8) * t729 - t776 + (pkin(8) * t768 + t747 * t782 - g(3)) * t751;
t709 = -g(3) * t747 + t751 * t723;
t733 = qJD(2) * pkin(2) - pkin(8) * t770;
t743 = t751 ^ 2;
t683 = pkin(8) * t730 - qJD(2) * t733 - t743 * t782 + t709;
t653 = t750 * t682 - t746 * t683;
t704 = -pkin(3) * t720 - pkin(9) * t721;
t740 = t742 ^ 2;
t758 = pkin(3) * t741 + pkin(9) * t740 - t721 * t704 + t653;
t716 = qJD(4) - t720;
t777 = t706 * t716;
t785 = (-t658 + t777) * qJ(5) - t758;
t784 = 2 * qJD(5);
t781 = -mrSges(5,3) - mrSges(6,2);
t654 = t746 * t682 + t750 * t683;
t691 = -qJD(3) * t721 - t729 * t746 + t750 * t730;
t703 = -mrSges(4,1) * t720 + mrSges(4,2) * t721;
t711 = mrSges(4,1) * t742 - mrSges(4,3) * t721;
t734 = g(1) * t748 - t752 * g(2);
t760 = -qJDD(1) * pkin(1) - t734;
t693 = -pkin(2) * t730 + t733 * t770 + (-pkin(8) * t743 - pkin(7)) * t753 + t760;
t642 = (-t720 * t742 - t692) * pkin(9) + (t721 * t742 - t691) * pkin(3) + t693;
t646 = -pkin(3) * t740 + pkin(9) * t741 + t704 * t720 + t654;
t636 = t745 * t642 + t783 * t646;
t707 = t783 * t721 + t745 * t742;
t657 = qJD(4) * t707 + t692 * t745 - t783 * t741;
t690 = qJDD(4) - t691;
t696 = mrSges(5,1) * t716 - mrSges(5,3) * t707;
t679 = pkin(4) * t706 - qJ(5) * t707;
t715 = t716 ^ 2;
t632 = -pkin(4) * t715 + t690 * qJ(5) - t706 * t679 + t716 * t784 + t636;
t697 = -mrSges(6,1) * t716 + mrSges(6,2) * t707;
t635 = t783 * t642 - t745 * t646;
t633 = -t690 * pkin(4) - t715 * qJ(5) + t707 * t679 + qJDD(5) - t635;
t627 = (-t658 - t777) * pkin(10) + (t706 * t707 - t690) * pkin(5) + t633;
t698 = -pkin(5) * t716 - pkin(10) * t707;
t705 = t706 ^ 2;
t628 = -pkin(5) * t705 + pkin(10) * t657 + t698 * t716 + t632;
t744 = sin(qJ(6));
t749 = cos(qJ(6));
t625 = t627 * t749 - t628 * t744;
t673 = t706 * t749 - t707 * t744;
t640 = qJD(6) * t673 + t657 * t744 + t658 * t749;
t674 = t706 * t744 + t707 * t749;
t652 = -mrSges(7,1) * t673 + mrSges(7,2) * t674;
t714 = qJD(6) - t716;
t661 = -mrSges(7,2) * t714 + mrSges(7,3) * t673;
t687 = qJDD(6) - t690;
t623 = m(7) * t625 + mrSges(7,1) * t687 - mrSges(7,3) * t640 - t652 * t674 + t661 * t714;
t626 = t627 * t744 + t628 * t749;
t639 = -qJD(6) * t674 + t657 * t749 - t658 * t744;
t662 = mrSges(7,1) * t714 - mrSges(7,3) * t674;
t624 = m(7) * t626 - mrSges(7,2) * t687 + mrSges(7,3) * t639 + t652 * t673 - t662 * t714;
t763 = -t623 * t744 + t749 * t624;
t759 = m(6) * t632 + t690 * mrSges(6,3) + t716 * t697 + t763;
t680 = mrSges(6,1) * t706 - mrSges(6,3) * t707;
t771 = -mrSges(5,1) * t706 - mrSges(5,2) * t707 - t680;
t613 = m(5) * t636 - mrSges(5,2) * t690 + t781 * t657 - t696 * t716 + t771 * t706 + t759;
t695 = -mrSges(5,2) * t716 - mrSges(5,3) * t706;
t616 = t623 * t749 + t624 * t744;
t694 = -mrSges(6,2) * t706 + mrSges(6,3) * t716;
t757 = -m(6) * t633 + t690 * mrSges(6,1) + t716 * t694 - t616;
t615 = m(5) * t635 + mrSges(5,1) * t690 + t781 * t658 + t695 * t716 + t771 * t707 + t757;
t764 = t783 * t613 - t615 * t745;
t609 = m(4) * t654 - mrSges(4,2) * t741 + mrSges(4,3) * t691 + t703 * t720 - t711 * t742 + t764;
t710 = -mrSges(4,2) * t742 + mrSges(4,3) * t720;
t634 = -0.2e1 * qJD(5) * t707 + (t707 * t716 + t657) * pkin(4) + t785;
t630 = -pkin(10) * t705 + (-pkin(4) - pkin(5)) * t657 + (-pkin(4) * t716 + t698 + t784) * t707 - t785;
t761 = -m(7) * t630 + t639 * mrSges(7,1) - t640 * mrSges(7,2) + t673 * t661 - t674 * t662;
t621 = m(6) * t634 + mrSges(6,1) * t657 - t658 * mrSges(6,3) + t694 * t706 - t707 * t697 + t761;
t754 = m(5) * t758 - t657 * mrSges(5,1) - mrSges(5,2) * t658 - t706 * t695 - t696 * t707 - t621;
t620 = m(4) * t653 + mrSges(4,1) * t741 - mrSges(4,3) * t692 - t703 * t721 + t710 * t742 + t754;
t604 = t746 * t609 + t750 * t620;
t708 = -g(3) * t751 - t776;
t728 = (-mrSges(3,1) * t751 + mrSges(3,2) * t747) * qJD(1);
t732 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t769;
t602 = m(3) * t708 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t729 + qJD(2) * t732 - t728 * t770 + t604;
t731 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t770;
t765 = t750 * t609 - t620 * t746;
t603 = m(3) * t709 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t730 - qJD(2) * t731 + t728 * t769 + t765;
t766 = -t602 * t747 + t751 * t603;
t595 = m(2) * t735 - mrSges(2,1) * t753 - qJDD(1) * mrSges(2,2) + t766;
t722 = -pkin(7) * t753 + t760;
t610 = t745 * t613 + t783 * t615;
t756 = m(4) * t693 - t691 * mrSges(4,1) + mrSges(4,2) * t692 - t720 * t710 + t711 * t721 + t610;
t755 = -m(3) * t722 + t730 * mrSges(3,1) - mrSges(3,2) * t729 - t731 * t770 + t732 * t769 - t756;
t606 = m(2) * t734 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t753 + t755;
t775 = t748 * t595 + t752 * t606;
t596 = t751 * t602 + t747 * t603;
t774 = t706 * t787 - t707 * t780 - t716 * t778;
t773 = t706 * t778 - t707 * t779 + t716 * t786;
t772 = -t780 * t706 + t707 * t788 + t779 * t716;
t767 = t752 * t595 - t606 * t748;
t719 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t747 + Ifges(3,4) * t751) * qJD(1);
t718 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t747 + Ifges(3,2) * t751) * qJD(1);
t717 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t747 + Ifges(3,6) * t751) * qJD(1);
t701 = Ifges(4,1) * t721 + Ifges(4,4) * t720 + Ifges(4,5) * t742;
t700 = Ifges(4,4) * t721 + Ifges(4,2) * t720 + Ifges(4,6) * t742;
t699 = Ifges(4,5) * t721 + Ifges(4,6) * t720 + Ifges(4,3) * t742;
t649 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t714;
t648 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t714;
t647 = Ifges(7,5) * t674 + Ifges(7,6) * t673 + Ifges(7,3) * t714;
t618 = mrSges(7,2) * t630 - mrSges(7,3) * t625 + Ifges(7,1) * t640 + Ifges(7,4) * t639 + Ifges(7,5) * t687 + t647 * t673 - t648 * t714;
t617 = -mrSges(7,1) * t630 + mrSges(7,3) * t626 + Ifges(7,4) * t640 + Ifges(7,2) * t639 + Ifges(7,6) * t687 - t647 * t674 + t649 * t714;
t598 = -mrSges(5,2) * t758 + mrSges(6,2) * t633 - mrSges(5,3) * t635 - mrSges(6,3) * t634 - pkin(10) * t616 - qJ(5) * t621 - t617 * t744 + t618 * t749 - t780 * t657 + t658 * t788 + t779 * t690 + t773 * t706 + t774 * t716;
t597 = mrSges(5,1) * t758 - mrSges(6,1) * t634 + mrSges(6,2) * t632 + mrSges(5,3) * t636 - pkin(4) * t621 - pkin(5) * t761 - pkin(10) * t763 - t749 * t617 - t744 * t618 - t657 * t787 + t780 * t658 + t778 * t690 + t773 * t707 + t772 * t716;
t592 = (qJ(5) * mrSges(6,2) + t778) * t657 + (pkin(4) * mrSges(6,2) - t779) * t658 - qJ(5) * t759 - pkin(4) * t757 + (qJ(5) * t680 - t772) * t706 + (pkin(4) * t680 + t774) * t707 + pkin(5) * t616 - pkin(3) * t610 - mrSges(7,2) * t626 + mrSges(7,1) * t625 + t742 * t701 - mrSges(6,3) * t632 + mrSges(6,1) * t633 - mrSges(5,1) * t635 + mrSges(5,2) * t636 + Ifges(7,6) * t639 + Ifges(4,6) * t741 + Ifges(7,5) * t640 + Ifges(7,3) * t687 + Ifges(4,2) * t691 + Ifges(4,4) * t692 - mrSges(4,1) * t693 - t673 * t649 + t674 * t648 + t786 * t690 + mrSges(4,3) * t654 - t721 * t699;
t591 = mrSges(4,2) * t693 - mrSges(4,3) * t653 + Ifges(4,1) * t692 + Ifges(4,4) * t691 + Ifges(4,5) * t741 - pkin(9) * t610 - t745 * t597 + t783 * t598 + t720 * t699 - t742 * t700;
t590 = mrSges(2,1) * g(3) + t753 * Ifges(2,5) - t721 * t700 + t720 * t701 + Ifges(2,6) * qJDD(1) - pkin(1) * t596 + mrSges(2,3) * t735 - pkin(2) * t604 - Ifges(3,5) * t729 - Ifges(3,6) * t730 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t708 + mrSges(3,2) * t709 - pkin(9) * t764 - t745 * t598 - t783 * t597 - pkin(3) * t754 - Ifges(4,5) * t692 - Ifges(4,6) * t691 - Ifges(4,3) * t741 - mrSges(4,1) * t653 + mrSges(4,2) * t654 + (-t718 * t747 + t719 * t751) * qJD(1);
t589 = mrSges(3,2) * t722 - mrSges(3,3) * t708 + Ifges(3,1) * t729 + Ifges(3,4) * t730 + Ifges(3,5) * qJDD(2) - pkin(8) * t604 - qJD(2) * t718 + t591 * t750 - t592 * t746 + t717 * t769;
t588 = -mrSges(3,1) * t722 + mrSges(3,3) * t709 + Ifges(3,4) * t729 + Ifges(3,2) * t730 + Ifges(3,6) * qJDD(2) - pkin(2) * t756 + pkin(8) * t765 + qJD(2) * t719 + t746 * t591 + t750 * t592 - t717 * t770;
t587 = -mrSges(2,2) * g(3) - mrSges(2,3) * t734 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t753 - pkin(7) * t596 - t588 * t747 + t589 * t751;
t1 = [-m(1) * g(1) + t767; -m(1) * g(2) + t775; (-m(1) - m(2)) * g(3) + t596; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t775 + t752 * t587 - t748 * t590; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t767 + t748 * t587 + t752 * t590; -mrSges(1,1) * g(2) + mrSges(2,1) * t734 + mrSges(1,2) * g(1) - mrSges(2,2) * t735 + Ifges(2,3) * qJDD(1) + pkin(1) * t755 + pkin(7) * t766 + t751 * t588 + t747 * t589;];
tauB  = t1;
