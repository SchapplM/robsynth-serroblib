% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:01:45
% EndTime: 2019-05-07 18:01:57
% DurationCPUTime: 6.56s
% Computational Cost: add. (75417->343), mult. (151128->406), div. (0->0), fcn. (104024->8), ass. (0->130)
t780 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t762 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t761 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t779 = Ifges(5,2) + Ifges(6,3) + Ifges(7,2);
t760 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t778 = -Ifges(5,3) - Ifges(6,2) - Ifges(7,3);
t736 = sin(qJ(3));
t739 = cos(qJ(2));
t764 = qJD(1) * t739;
t737 = sin(qJ(2));
t765 = qJD(1) * t737;
t774 = cos(qJ(3));
t711 = -t736 * t765 + t774 * t764;
t763 = qJD(1) * qJD(2);
t720 = t737 * qJDD(1) + t739 * t763;
t721 = t739 * qJDD(1) - t737 * t763;
t682 = t711 * qJD(3) + t774 * t720 + t736 * t721;
t712 = (t736 * t739 + t774 * t737) * qJD(1);
t733 = qJD(2) + qJD(3);
t735 = sin(qJ(4));
t773 = cos(qJ(4));
t698 = t735 * t712 - t773 * t733;
t732 = qJDD(2) + qJDD(3);
t645 = -t698 * qJD(4) + t773 * t682 + t735 * t732;
t738 = sin(qJ(1));
t740 = cos(qJ(1));
t726 = -t740 * g(1) - t738 * g(2);
t741 = qJD(1) ^ 2;
t714 = -t741 * pkin(1) + qJDD(1) * pkin(7) + t726;
t768 = t737 * t714;
t772 = pkin(2) * t741;
t673 = qJDD(2) * pkin(2) - t720 * pkin(8) - t768 + (pkin(8) * t763 + t737 * t772 - g(3)) * t739;
t701 = -t737 * g(3) + t739 * t714;
t724 = qJD(2) * pkin(2) - pkin(8) * t765;
t734 = t739 ^ 2;
t674 = t721 * pkin(8) - qJD(2) * t724 - t734 * t772 + t701;
t636 = t774 * t673 - t736 * t674;
t696 = -t711 * pkin(3) - t712 * pkin(9);
t731 = t733 ^ 2;
t746 = t732 * pkin(3) + t731 * pkin(9) - t712 * t696 + t636;
t707 = qJD(4) - t711;
t770 = t698 * t707;
t777 = (-t645 + t770) * qJ(5) - t746;
t699 = t773 * t712 + t735 * t733;
t776 = -0.2e1 * t699;
t775 = 2 * qJD(5);
t771 = -mrSges(5,3) - mrSges(6,2);
t685 = t707 * mrSges(7,2) + t698 * mrSges(7,3);
t769 = t707 * t685;
t637 = t736 * t673 + t774 * t674;
t681 = -t712 * qJD(3) - t736 * t720 + t774 * t721;
t695 = -t711 * mrSges(4,1) + t712 * mrSges(4,2);
t703 = t733 * mrSges(4,1) - t712 * mrSges(4,3);
t725 = t738 * g(1) - t740 * g(2);
t748 = -qJDD(1) * pkin(1) - t725;
t683 = -t721 * pkin(2) + t724 * t765 + (-pkin(8) * t734 - pkin(7)) * t741 + t748;
t631 = (-t711 * t733 - t682) * pkin(9) + (t712 * t733 - t681) * pkin(3) + t683;
t635 = -t731 * pkin(3) + t732 * pkin(9) + t711 * t696 + t637;
t628 = t773 * t631 - t735 * t635;
t680 = qJDD(4) - t681;
t686 = -t707 * mrSges(5,2) - t698 * mrSges(5,3);
t669 = t698 * pkin(4) - t699 * qJ(5);
t706 = t707 ^ 2;
t626 = -t680 * pkin(4) - t706 * qJ(5) + t699 * t669 + qJDD(5) - t628;
t684 = -t698 * mrSges(6,2) + t707 * mrSges(6,3);
t619 = qJD(6) * t776 + (-t645 - t770) * qJ(6) + (t698 * t699 - t680) * pkin(5) + t626;
t671 = -t698 * mrSges(7,1) + t699 * mrSges(7,2);
t750 = -m(7) * t619 + t645 * mrSges(7,3) + t699 * t671;
t745 = -m(6) * t626 + t680 * mrSges(6,1) + t707 * t684 + t750;
t670 = t698 * mrSges(6,1) - t699 * mrSges(6,3);
t766 = -t698 * mrSges(5,1) - t699 * mrSges(5,2) - t670;
t611 = m(5) * t628 + (t685 + t686) * t707 + t766 * t699 + (mrSges(5,1) + mrSges(7,1)) * t680 + t771 * t645 + t745;
t629 = t735 * t631 + t773 * t635;
t644 = t699 * qJD(4) + t735 * t682 - t773 * t732;
t688 = -t707 * mrSges(7,1) - t699 * mrSges(7,3);
t689 = t707 * mrSges(5,1) - t699 * mrSges(5,3);
t625 = -t706 * pkin(4) + t680 * qJ(5) - t698 * t669 + t707 * t775 + t629;
t690 = -t707 * mrSges(6,1) + t699 * mrSges(6,2);
t687 = -t707 * pkin(5) - t699 * qJ(6);
t697 = t698 ^ 2;
t621 = -t697 * pkin(5) + t644 * qJ(6) + 0.2e1 * qJD(6) * t698 + t707 * t687 + t625;
t759 = m(7) * t621 + t644 * mrSges(7,3) + t698 * t671;
t747 = m(6) * t625 + t680 * mrSges(6,3) + t707 * t690 + t759;
t616 = m(5) * t629 + (t688 - t689) * t707 + t766 * t698 + (-mrSges(5,2) + mrSges(7,2)) * t680 + t771 * t644 + t747;
t752 = -t735 * t611 + t773 * t616;
t608 = m(4) * t637 - t732 * mrSges(4,2) + t681 * mrSges(4,3) + t711 * t695 - t733 * t703 + t752;
t702 = -t733 * mrSges(4,2) + t711 * mrSges(4,3);
t627 = qJD(5) * t776 + (t699 * t707 + t644) * pkin(4) + t777;
t623 = -t697 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t644 + (-pkin(4) * t707 + t687 + t775) * t699 - t777;
t749 = -m(7) * t623 + t644 * mrSges(7,1) - t645 * mrSges(7,2) + t698 * t685 - t699 * t688;
t617 = m(6) * t627 + t644 * mrSges(6,1) - t645 * mrSges(6,3) + t698 * t684 - t699 * t690 + t749;
t742 = m(5) * t746 - t644 * mrSges(5,1) - t645 * mrSges(5,2) - t698 * t686 - t699 * t689 - t617;
t613 = m(4) * t636 + t732 * mrSges(4,1) - t682 * mrSges(4,3) - t712 * t695 + t733 * t702 + t742;
t602 = t736 * t608 + t774 * t613;
t700 = -t739 * g(3) - t768;
t719 = (-mrSges(3,1) * t739 + mrSges(3,2) * t737) * qJD(1);
t723 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t764;
t600 = m(3) * t700 + qJDD(2) * mrSges(3,1) - t720 * mrSges(3,3) + qJD(2) * t723 - t719 * t765 + t602;
t722 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t765;
t753 = t774 * t608 - t736 * t613;
t601 = m(3) * t701 - qJDD(2) * mrSges(3,2) + t721 * mrSges(3,3) - qJD(2) * t722 + t719 * t764 + t753;
t754 = -t737 * t600 + t739 * t601;
t593 = m(2) * t726 - t741 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t754;
t713 = -t741 * pkin(7) + t748;
t609 = t773 * t611 + t735 * t616;
t744 = m(4) * t683 - t681 * mrSges(4,1) + t682 * mrSges(4,2) - t711 * t702 + t712 * t703 + t609;
t743 = -m(3) * t713 + t721 * mrSges(3,1) - t720 * mrSges(3,2) - t722 * t765 + t723 * t764 - t744;
t605 = m(2) * t725 + qJDD(1) * mrSges(2,1) - t741 * mrSges(2,2) + t743;
t767 = t738 * t593 + t740 * t605;
t594 = t739 * t600 + t737 * t601;
t758 = t760 * t698 - t761 * t699 + t778 * t707;
t757 = t779 * t698 - t762 * t699 - t760 * t707;
t756 = t762 * t698 - t780 * t699 - t761 * t707;
t755 = t740 * t593 - t738 * t605;
t710 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t737 + Ifges(3,4) * t739) * qJD(1);
t709 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t737 + Ifges(3,2) * t739) * qJD(1);
t708 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t737 + Ifges(3,6) * t739) * qJD(1);
t693 = Ifges(4,1) * t712 + Ifges(4,4) * t711 + Ifges(4,5) * t733;
t692 = Ifges(4,4) * t712 + Ifges(4,2) * t711 + Ifges(4,6) * t733;
t691 = Ifges(4,5) * t712 + Ifges(4,6) * t711 + Ifges(4,3) * t733;
t618 = -t680 * mrSges(7,1) - t750 - t769;
t603 = -mrSges(5,2) * t746 + mrSges(6,2) * t626 + mrSges(7,2) * t623 - mrSges(5,3) * t628 - mrSges(6,3) * t627 - mrSges(7,3) * t619 - qJ(5) * t617 - qJ(6) * t618 - t762 * t644 + t780 * t645 + t761 * t680 + t758 * t698 + t757 * t707;
t596 = mrSges(5,1) * t746 + mrSges(5,3) * t629 - mrSges(6,1) * t627 + mrSges(6,2) * t625 + mrSges(7,1) * t623 - mrSges(7,3) * t621 - pkin(5) * t749 - qJ(6) * t759 - pkin(4) * t617 + (-qJ(6) * t688 - t756) * t707 + t758 * t699 + (-qJ(6) * mrSges(7,2) + t760) * t680 + t762 * t645 - t779 * t644;
t595 = -pkin(4) * (t745 + t769) - qJ(5) * (t707 * t688 + t747) + Ifges(4,6) * t732 + t733 * t693 - t712 * t691 + Ifges(4,2) * t681 + Ifges(4,4) * t682 - mrSges(4,1) * t683 + mrSges(4,3) * t637 + mrSges(6,1) * t626 - mrSges(5,1) * t628 + mrSges(5,2) * t629 - mrSges(6,3) * t625 + mrSges(7,1) * t619 - mrSges(7,2) * t621 + pkin(5) * t618 - pkin(3) * t609 + (pkin(4) * t670 + t757) * t699 + (qJ(5) * t670 + t756) * t698 + (-pkin(4) * mrSges(7,1) - qJ(5) * mrSges(7,2) + t778) * t680 + (pkin(4) * mrSges(6,2) - t761) * t645 + (qJ(5) * mrSges(6,2) + t760) * t644;
t590 = mrSges(4,2) * t683 - mrSges(4,3) * t636 + Ifges(4,1) * t682 + Ifges(4,4) * t681 + Ifges(4,5) * t732 - pkin(9) * t609 - t735 * t596 + t773 * t603 + t711 * t691 - t733 * t692;
t589 = mrSges(3,2) * t713 - mrSges(3,3) * t700 + Ifges(3,1) * t720 + Ifges(3,4) * t721 + Ifges(3,5) * qJDD(2) - pkin(8) * t602 - qJD(2) * t709 + t774 * t590 - t736 * t595 + t708 * t764;
t588 = -pkin(1) * t594 + mrSges(2,3) * t726 - pkin(2) * t602 - Ifges(3,5) * t720 - Ifges(3,6) * t721 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t700 + mrSges(3,2) * t701 - pkin(9) * t752 - mrSges(4,1) * t636 + mrSges(4,2) * t637 - t735 * t603 - t773 * t596 - pkin(3) * t742 - Ifges(4,5) * t682 - Ifges(4,6) * t681 - Ifges(4,3) * t732 + mrSges(2,1) * g(3) - t712 * t692 + t711 * t693 + t741 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-t737 * t709 + t739 * t710) * qJD(1);
t587 = -mrSges(3,1) * t713 + mrSges(3,3) * t701 + Ifges(3,4) * t720 + Ifges(3,2) * t721 + Ifges(3,6) * qJDD(2) - pkin(2) * t744 + pkin(8) * t753 + qJD(2) * t710 + t736 * t590 + t774 * t595 - t708 * t765;
t586 = -mrSges(2,2) * g(3) - mrSges(2,3) * t725 + Ifges(2,5) * qJDD(1) - t741 * Ifges(2,6) - pkin(7) * t594 - t737 * t587 + t739 * t589;
t1 = [-m(1) * g(1) + t755; -m(1) * g(2) + t767; (-m(1) - m(2)) * g(3) + t594; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t767 + t740 * t586 - t738 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t755 + t738 * t586 + t740 * t588; -mrSges(1,1) * g(2) + mrSges(2,1) * t725 + mrSges(1,2) * g(1) - mrSges(2,2) * t726 + Ifges(2,3) * qJDD(1) + pkin(1) * t743 + pkin(7) * t754 + t739 * t587 + t737 * t589;];
tauB  = t1;
