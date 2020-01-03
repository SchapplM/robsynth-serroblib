% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:34
% EndTime: 2019-12-31 21:56:41
% DurationCPUTime: 5.01s
% Computational Cost: add. (55163->291), mult. (111369->357), div. (0->0), fcn. (74604->8), ass. (0->118)
t781 = Ifges(5,1) + Ifges(6,1);
t773 = Ifges(5,4) - Ifges(6,5);
t772 = -Ifges(5,5) - Ifges(6,4);
t780 = Ifges(5,2) + Ifges(6,3);
t771 = Ifges(5,6) - Ifges(6,6);
t779 = -Ifges(5,3) - Ifges(6,2);
t738 = sin(qJ(3));
t739 = sin(qJ(2));
t741 = cos(qJ(3));
t742 = cos(qJ(2));
t712 = (t738 * t739 - t741 * t742) * qJD(1);
t761 = qJD(1) * qJD(2);
t720 = t739 * qJDD(1) + t742 * t761;
t740 = sin(qJ(1));
t743 = cos(qJ(1));
t727 = -t743 * g(1) - t740 * g(2);
t744 = qJD(1) ^ 2;
t715 = -t744 * pkin(1) + qJDD(1) * pkin(6) + t727;
t769 = t739 * t715;
t775 = pkin(2) * t744;
t678 = qJDD(2) * pkin(2) - t720 * pkin(7) - t769 + (pkin(7) * t761 + t739 * t775 - g(3)) * t742;
t703 = -t739 * g(3) + t742 * t715;
t721 = t742 * qJDD(1) - t739 * t761;
t763 = qJD(1) * t739;
t725 = qJD(2) * pkin(2) - pkin(7) * t763;
t736 = t742 ^ 2;
t679 = t721 * pkin(7) - qJD(2) * t725 - t736 * t775 + t703;
t656 = t738 * t678 + t741 * t679;
t713 = (t738 * t742 + t739 * t741) * qJD(1);
t687 = -t713 * qJD(3) - t738 * t720 + t741 * t721;
t698 = t712 * mrSges(4,1) + t713 * mrSges(4,2);
t734 = qJD(2) + qJD(3);
t705 = t734 * mrSges(4,1) - t713 * mrSges(4,3);
t733 = qJDD(2) + qJDD(3);
t688 = -t712 * qJD(3) + t741 * t720 + t738 * t721;
t726 = t740 * g(1) - t743 * g(2);
t752 = -qJDD(1) * pkin(1) - t726;
t689 = -t721 * pkin(2) + t725 * t763 + (-pkin(7) * t736 - pkin(6)) * t744 + t752;
t650 = (t712 * t734 - t688) * pkin(8) + (t713 * t734 - t687) * pkin(3) + t689;
t699 = t712 * pkin(3) - t713 * pkin(8);
t732 = t734 ^ 2;
t653 = -t732 * pkin(3) + t733 * pkin(8) - t712 * t699 + t656;
t737 = sin(qJ(4));
t776 = cos(qJ(4));
t648 = t737 * t650 + t776 * t653;
t701 = t776 * t713 + t737 * t734;
t659 = t701 * qJD(4) + t737 * t688 - t776 * t733;
t686 = qJDD(4) - t687;
t708 = qJD(4) + t712;
t692 = t708 * mrSges(5,1) - t701 * mrSges(5,3);
t700 = t737 * t713 - t776 * t734;
t674 = t700 * pkin(4) - t701 * qJ(5);
t707 = t708 ^ 2;
t644 = -t707 * pkin(4) + t686 * qJ(5) + 0.2e1 * qJD(5) * t708 - t700 * t674 + t648;
t693 = -t708 * mrSges(6,1) + t701 * mrSges(6,2);
t760 = m(6) * t644 + t686 * mrSges(6,3) + t708 * t693;
t675 = t700 * mrSges(6,1) - t701 * mrSges(6,3);
t764 = -t700 * mrSges(5,1) - t701 * mrSges(5,2) - t675;
t774 = -mrSges(5,3) - mrSges(6,2);
t635 = m(5) * t648 - t686 * mrSges(5,2) + t774 * t659 - t708 * t692 + t764 * t700 + t760;
t647 = t776 * t650 - t737 * t653;
t660 = -t700 * qJD(4) + t776 * t688 + t737 * t733;
t691 = -t708 * mrSges(5,2) - t700 * mrSges(5,3);
t645 = -t686 * pkin(4) - t707 * qJ(5) + t701 * t674 + qJDD(5) - t647;
t690 = -t700 * mrSges(6,2) + t708 * mrSges(6,3);
t755 = -m(6) * t645 + t686 * mrSges(6,1) + t708 * t690;
t637 = m(5) * t647 + t686 * mrSges(5,1) + t774 * t660 + t708 * t691 + t764 * t701 + t755;
t756 = t776 * t635 - t737 * t637;
t623 = m(4) * t656 - t733 * mrSges(4,2) + t687 * mrSges(4,3) - t712 * t698 - t734 * t705 + t756;
t655 = t741 * t678 - t738 * t679;
t704 = -t734 * mrSges(4,2) - t712 * mrSges(4,3);
t652 = -t733 * pkin(3) - t732 * pkin(8) + t713 * t699 - t655;
t646 = -0.2e1 * qJD(5) * t701 + (t700 * t708 - t660) * qJ(5) + (t701 * t708 + t659) * pkin(4) + t652;
t642 = m(6) * t646 + t659 * mrSges(6,1) - t660 * mrSges(6,3) + t700 * t690 - t701 * t693;
t747 = -m(5) * t652 - t659 * mrSges(5,1) - t660 * mrSges(5,2) - t700 * t691 - t701 * t692 - t642;
t632 = m(4) * t655 + t733 * mrSges(4,1) - t688 * mrSges(4,3) - t713 * t698 + t734 * t704 + t747;
t617 = t738 * t623 + t741 * t632;
t702 = -t742 * g(3) - t769;
t710 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t739 + Ifges(3,2) * t742) * qJD(1);
t711 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t739 + Ifges(3,4) * t742) * qJD(1);
t765 = t773 * t700 - t781 * t701 + t772 * t708;
t767 = t771 * t700 + t772 * t701 + t779 * t708;
t625 = -mrSges(5,1) * t652 - mrSges(6,1) * t646 + mrSges(6,2) * t644 + mrSges(5,3) * t648 - pkin(4) * t642 - t780 * t659 + t773 * t660 + t771 * t686 + t767 * t701 - t765 * t708;
t766 = t780 * t700 - t773 * t701 - t771 * t708;
t627 = mrSges(5,2) * t652 + mrSges(6,2) * t645 - mrSges(5,3) * t647 - mrSges(6,3) * t646 - qJ(5) * t642 - t773 * t659 + t781 * t660 - t772 * t686 + t767 * t700 + t766 * t708;
t695 = Ifges(4,4) * t713 - Ifges(4,2) * t712 + Ifges(4,6) * t734;
t696 = Ifges(4,1) * t713 - Ifges(4,4) * t712 + Ifges(4,5) * t734;
t749 = -mrSges(4,1) * t655 + mrSges(4,2) * t656 - Ifges(4,5) * t688 - Ifges(4,6) * t687 - Ifges(4,3) * t733 - pkin(3) * t747 - pkin(8) * t756 - t776 * t625 - t737 * t627 - t713 * t695 - t712 * t696;
t778 = mrSges(3,1) * t702 - mrSges(3,2) * t703 + Ifges(3,5) * t720 + Ifges(3,6) * t721 + Ifges(3,3) * qJDD(2) + pkin(2) * t617 + (t739 * t710 - t742 * t711) * qJD(1) - t749;
t641 = t660 * mrSges(6,2) + t701 * t675 - t755;
t777 = -t771 * t659 - t772 * t660 - t779 * t686 - t765 * t700 - t766 * t701 + mrSges(5,1) * t647 - mrSges(6,1) * t645 - mrSges(5,2) * t648 + mrSges(6,3) * t644 - pkin(4) * t641 + qJ(5) * (-t659 * mrSges(6,2) - t700 * t675 + t760);
t719 = (-mrSges(3,1) * t742 + mrSges(3,2) * t739) * qJD(1);
t762 = qJD(1) * t742;
t724 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t762;
t615 = m(3) * t702 + qJDD(2) * mrSges(3,1) - t720 * mrSges(3,3) + qJD(2) * t724 - t719 * t763 + t617;
t723 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t763;
t757 = t741 * t623 - t738 * t632;
t616 = m(3) * t703 - qJDD(2) * mrSges(3,2) + t721 * mrSges(3,3) - qJD(2) * t723 + t719 * t762 + t757;
t758 = -t739 * t615 + t742 * t616;
t608 = m(2) * t727 - t744 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t758;
t714 = -t744 * pkin(6) + t752;
t629 = t737 * t635 + t776 * t637;
t750 = m(4) * t689 - t687 * mrSges(4,1) + t688 * mrSges(4,2) + t712 * t704 + t713 * t705 + t629;
t746 = -m(3) * t714 + t721 * mrSges(3,1) - t720 * mrSges(3,2) - t723 * t763 + t724 * t762 - t750;
t619 = m(2) * t726 + qJDD(1) * mrSges(2,1) - t744 * mrSges(2,2) + t746;
t768 = t740 * t608 + t743 * t619;
t610 = t742 * t615 + t739 * t616;
t759 = t743 * t608 - t740 * t619;
t694 = Ifges(4,5) * t713 - Ifges(4,6) * t712 + Ifges(4,3) * t734;
t605 = mrSges(4,2) * t689 - mrSges(4,3) * t655 + Ifges(4,1) * t688 + Ifges(4,4) * t687 + Ifges(4,5) * t733 - pkin(8) * t629 - t737 * t625 + t776 * t627 - t712 * t694 - t734 * t695;
t611 = -mrSges(4,1) * t689 + mrSges(4,3) * t656 + Ifges(4,4) * t688 + Ifges(4,2) * t687 + Ifges(4,6) * t733 - pkin(3) * t629 - t713 * t694 + t734 * t696 - t777;
t709 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t739 + Ifges(3,6) * t742) * qJD(1);
t602 = -mrSges(3,1) * t714 + mrSges(3,3) * t703 + Ifges(3,4) * t720 + Ifges(3,2) * t721 + Ifges(3,6) * qJDD(2) - pkin(2) * t750 + pkin(7) * t757 + qJD(2) * t711 + t738 * t605 + t741 * t611 - t709 * t763;
t604 = mrSges(3,2) * t714 - mrSges(3,3) * t702 + Ifges(3,1) * t720 + Ifges(3,4) * t721 + Ifges(3,5) * qJDD(2) - pkin(7) * t617 - qJD(2) * t710 + t741 * t605 - t738 * t611 + t709 * t762;
t751 = mrSges(2,1) * t726 - mrSges(2,2) * t727 + Ifges(2,3) * qJDD(1) + pkin(1) * t746 + pkin(6) * t758 + t742 * t602 + t739 * t604;
t600 = mrSges(2,1) * g(3) + mrSges(2,3) * t727 + t744 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t610 - t778;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t726 + Ifges(2,5) * qJDD(1) - t744 * Ifges(2,6) - pkin(6) * t610 - t739 * t602 + t742 * t604;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t768 + t743 * t599 - t740 * t600; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t759 + t740 * t599 + t743 * t600; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t751; t751; t778; -t749; t777; t641;];
tauJB = t1;
