% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:57:11
% EndTime: 2019-12-31 19:57:16
% DurationCPUTime: 4.55s
% Computational Cost: add. (46175->291), mult. (104451->357), div. (0->0), fcn. (68791->8), ass. (0->117)
t778 = -2 * qJD(3);
t777 = Ifges(5,1) + Ifges(6,1);
t770 = Ifges(5,4) + Ifges(6,4);
t769 = Ifges(5,5) + Ifges(6,5);
t776 = Ifges(5,2) + Ifges(6,2);
t768 = Ifges(5,6) + Ifges(6,6);
t775 = Ifges(5,3) + Ifges(6,3);
t735 = sin(qJ(2));
t738 = cos(qJ(2));
t756 = qJD(1) * qJD(2);
t721 = t735 * qJDD(1) + t738 * t756;
t736 = sin(qJ(1));
t739 = cos(qJ(1));
t728 = -t739 * g(1) - t736 * g(2);
t741 = qJD(1) ^ 2;
t716 = -t741 * pkin(1) + qJDD(1) * pkin(6) + t728;
t765 = t735 * t716;
t772 = pkin(2) * t741;
t677 = qJDD(2) * pkin(2) - t721 * qJ(3) - t765 + (qJ(3) * t756 + t735 * t772 - g(3)) * t738;
t703 = -t735 * g(3) + t738 * t716;
t722 = t738 * qJDD(1) - t735 * t756;
t759 = qJD(1) * t735;
t724 = qJD(2) * pkin(2) - qJ(3) * t759;
t731 = t738 ^ 2;
t678 = t722 * qJ(3) - qJD(2) * t724 - t731 * t772 + t703;
t732 = sin(pkin(8));
t733 = cos(pkin(8));
t711 = (t732 * t738 + t733 * t735) * qJD(1);
t657 = t733 * t677 - t732 * t678 + t711 * t778;
t710 = (t732 * t735 - t733 * t738) * qJD(1);
t658 = t732 * t677 + t733 * t678 + t710 * t778;
t693 = t710 * pkin(3) - t711 * pkin(7);
t740 = qJD(2) ^ 2;
t653 = -t740 * pkin(3) + qJDD(2) * pkin(7) - t710 * t693 + t658;
t727 = t736 * g(1) - t739 * g(2);
t746 = -qJDD(1) * pkin(1) - t727;
t682 = -t722 * pkin(2) + qJDD(3) + t724 * t759 + (-qJ(3) * t731 - pkin(6)) * t741 + t746;
t697 = -t732 * t721 + t733 * t722;
t698 = t733 * t721 + t732 * t722;
t656 = (qJD(2) * t710 - t698) * pkin(7) + (qJD(2) * t711 - t697) * pkin(3) + t682;
t734 = sin(qJ(4));
t737 = cos(qJ(4));
t649 = -t734 * t653 + t737 * t656;
t700 = t737 * qJD(2) - t734 * t711;
t672 = t700 * qJD(4) + t734 * qJDD(2) + t737 * t698;
t701 = t734 * qJD(2) + t737 * t711;
t679 = -t700 * mrSges(6,1) + t701 * mrSges(6,2);
t680 = -t700 * mrSges(5,1) + t701 * mrSges(5,2);
t709 = qJD(4) + t710;
t684 = -t709 * mrSges(5,2) + t700 * mrSges(5,3);
t696 = qJDD(4) - t697;
t645 = -0.2e1 * qJD(5) * t701 + (t700 * t709 - t672) * qJ(5) + (t700 * t701 + t696) * pkin(4) + t649;
t683 = -t709 * mrSges(6,2) + t700 * mrSges(6,3);
t755 = m(6) * t645 + t696 * mrSges(6,1) + t709 * t683;
t635 = m(5) * t649 + t696 * mrSges(5,1) + t709 * t684 + (-t679 - t680) * t701 + (-mrSges(5,3) - mrSges(6,3)) * t672 + t755;
t650 = t737 * t653 + t734 * t656;
t671 = -t701 * qJD(4) + t737 * qJDD(2) - t734 * t698;
t685 = t709 * pkin(4) - t701 * qJ(5);
t699 = t700 ^ 2;
t647 = -t699 * pkin(4) + t671 * qJ(5) + 0.2e1 * qJD(5) * t700 - t709 * t685 + t650;
t754 = m(6) * t647 + t671 * mrSges(6,3) + t700 * t679;
t686 = t709 * mrSges(6,1) - t701 * mrSges(6,3);
t760 = -t709 * mrSges(5,1) + t701 * mrSges(5,3) - t686;
t771 = -mrSges(5,2) - mrSges(6,2);
t638 = m(5) * t650 + t671 * mrSges(5,3) + t700 * t680 + t696 * t771 + t709 * t760 + t754;
t633 = -t734 * t635 + t737 * t638;
t692 = t710 * mrSges(4,1) + t711 * mrSges(4,2);
t705 = qJD(2) * mrSges(4,1) - t711 * mrSges(4,3);
t629 = m(4) * t658 - qJDD(2) * mrSges(4,2) + t697 * mrSges(4,3) - qJD(2) * t705 - t710 * t692 + t633;
t652 = -qJDD(2) * pkin(3) - t740 * pkin(7) + t711 * t693 - t657;
t648 = -t671 * pkin(4) - t699 * qJ(5) + t701 * t685 + qJDD(5) + t652;
t749 = -m(6) * t648 + t671 * mrSges(6,1) + t700 * t683;
t641 = -m(5) * t652 + t671 * mrSges(5,1) + t771 * t672 + t700 * t684 + t760 * t701 + t749;
t704 = -qJD(2) * mrSges(4,2) - t710 * mrSges(4,3);
t640 = m(4) * t657 + qJDD(2) * mrSges(4,1) - t698 * mrSges(4,3) + qJD(2) * t704 - t711 * t692 + t641;
t622 = t732 * t629 + t733 * t640;
t643 = t672 * mrSges(6,2) + t701 * t686 - t749;
t761 = t770 * t700 + t777 * t701 + t769 * t709;
t763 = -t768 * t700 - t769 * t701 - t775 * t709;
t623 = -mrSges(5,1) * t652 + mrSges(5,3) * t650 - mrSges(6,1) * t648 + mrSges(6,3) * t647 - pkin(4) * t643 + qJ(5) * t754 + (-qJ(5) * t686 + t761) * t709 + t763 * t701 + (-qJ(5) * mrSges(6,2) + t768) * t696 + t770 * t672 + t776 * t671;
t642 = -t672 * mrSges(6,3) - t701 * t679 + t755;
t762 = -t776 * t700 - t770 * t701 - t768 * t709;
t630 = mrSges(5,2) * t652 + mrSges(6,2) * t648 - mrSges(5,3) * t649 - mrSges(6,3) * t645 - qJ(5) * t642 + t770 * t671 + t777 * t672 + t769 * t696 - t763 * t700 + t762 * t709;
t689 = Ifges(4,4) * t711 - Ifges(4,2) * t710 + Ifges(4,6) * qJD(2);
t690 = Ifges(4,1) * t711 - Ifges(4,4) * t710 + Ifges(4,5) * qJD(2);
t702 = -t738 * g(3) - t765;
t713 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t735 + Ifges(3,2) * t738) * qJD(1);
t714 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t735 + Ifges(3,4) * t738) * qJD(1);
t774 = (t735 * t713 - t738 * t714) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t702 + mrSges(4,1) * t657 - mrSges(3,2) * t703 - mrSges(4,2) * t658 + Ifges(3,5) * t721 + Ifges(4,5) * t698 + Ifges(3,6) * t722 + Ifges(4,6) * t697 + pkin(2) * t622 + pkin(3) * t641 + pkin(7) * t633 + t737 * t623 + t734 * t630 + t711 * t689 + t710 * t690;
t773 = mrSges(5,1) * t649 + mrSges(6,1) * t645 - mrSges(5,2) * t650 - mrSges(6,2) * t647 + pkin(4) * t642 + t671 * t768 + t672 * t769 + t775 * t696 - t700 * t761 - t762 * t701;
t720 = (-mrSges(3,1) * t738 + mrSges(3,2) * t735) * qJD(1);
t758 = qJD(1) * t738;
t726 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t758;
t620 = m(3) * t702 + qJDD(2) * mrSges(3,1) - t721 * mrSges(3,3) + qJD(2) * t726 - t720 * t759 + t622;
t725 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t759;
t751 = t733 * t629 - t732 * t640;
t621 = m(3) * t703 - qJDD(2) * mrSges(3,2) + t722 * mrSges(3,3) - qJD(2) * t725 + t720 * t758 + t751;
t752 = -t735 * t620 + t738 * t621;
t613 = m(2) * t728 - t741 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t752;
t632 = t737 * t635 + t734 * t638;
t631 = m(4) * t682 - t697 * mrSges(4,1) + t698 * mrSges(4,2) + t710 * t704 + t711 * t705 + t632;
t715 = -t741 * pkin(6) + t746;
t743 = -m(3) * t715 + t722 * mrSges(3,1) - t721 * mrSges(3,2) - t725 * t759 + t726 * t758 - t631;
t625 = m(2) * t727 + qJDD(1) * mrSges(2,1) - t741 * mrSges(2,2) + t743;
t764 = t736 * t613 + t739 * t625;
t615 = t738 * t620 + t735 * t621;
t753 = t739 * t613 - t736 * t625;
t688 = Ifges(4,5) * t711 - Ifges(4,6) * t710 + Ifges(4,3) * qJD(2);
t610 = mrSges(4,2) * t682 - mrSges(4,3) * t657 + Ifges(4,1) * t698 + Ifges(4,4) * t697 + Ifges(4,5) * qJDD(2) - pkin(7) * t632 - qJD(2) * t689 - t734 * t623 + t737 * t630 - t710 * t688;
t616 = -mrSges(4,1) * t682 + mrSges(4,3) * t658 + Ifges(4,4) * t698 + Ifges(4,2) * t697 + Ifges(4,6) * qJDD(2) - pkin(3) * t632 + qJD(2) * t690 - t711 * t688 - t773;
t712 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t735 + Ifges(3,6) * t738) * qJD(1);
t607 = -mrSges(3,1) * t715 + mrSges(3,3) * t703 + Ifges(3,4) * t721 + Ifges(3,2) * t722 + Ifges(3,6) * qJDD(2) - pkin(2) * t631 + qJ(3) * t751 + qJD(2) * t714 + t732 * t610 + t733 * t616 - t712 * t759;
t609 = mrSges(3,2) * t715 - mrSges(3,3) * t702 + Ifges(3,1) * t721 + Ifges(3,4) * t722 + Ifges(3,5) * qJDD(2) - qJ(3) * t622 - qJD(2) * t713 + t733 * t610 - t732 * t616 + t712 * t758;
t745 = mrSges(2,1) * t727 - mrSges(2,2) * t728 + Ifges(2,3) * qJDD(1) + pkin(1) * t743 + pkin(6) * t752 + t738 * t607 + t735 * t609;
t605 = mrSges(2,1) * g(3) + mrSges(2,3) * t728 + t741 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t615 - t774;
t604 = -mrSges(2,2) * g(3) - mrSges(2,3) * t727 + Ifges(2,5) * qJDD(1) - t741 * Ifges(2,6) - pkin(6) * t615 - t735 * t607 + t738 * t609;
t1 = [-m(1) * g(1) + t753; -m(1) * g(2) + t764; (-m(1) - m(2)) * g(3) + t615; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t764 + t739 * t604 - t736 * t605; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t753 + t736 * t604 + t739 * t605; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t745; t745; t774; t631; t773; t643;];
tauJB = t1;
