% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:48
% EndTime: 2019-12-31 18:53:53
% DurationCPUTime: 4.19s
% Computational Cost: add. (40447->267), mult. (94691->323), div. (0->0), fcn. (66013->8), ass. (0->118)
t768 = Ifges(5,1) + Ifges(6,1);
t761 = Ifges(5,4) - Ifges(6,5);
t760 = -Ifges(5,5) - Ifges(6,4);
t767 = Ifges(5,2) + Ifges(6,3);
t759 = Ifges(5,6) - Ifges(6,6);
t766 = -Ifges(5,3) - Ifges(6,2);
t727 = qJD(1) ^ 2;
t719 = sin(pkin(8));
t720 = cos(pkin(8));
t722 = sin(qJ(3));
t724 = cos(qJ(3));
t734 = t719 * t722 - t720 * t724;
t700 = t734 * qJD(1);
t735 = t719 * t724 + t720 * t722;
t701 = t735 * qJD(1);
t748 = t701 * qJD(3);
t688 = -t734 * qJDD(1) - t748;
t749 = t700 * qJD(3);
t689 = t735 * qJDD(1) - t749;
t721 = sin(qJ(4));
t764 = cos(qJ(4));
t692 = -t764 * qJD(3) + t721 * t701;
t661 = -t692 * qJD(4) + t721 * qJDD(3) + t764 * t689;
t693 = t721 * qJD(3) + t764 * t701;
t665 = t692 * mrSges(6,1) - t693 * mrSges(6,3);
t723 = sin(qJ(1));
t725 = cos(qJ(1));
t707 = -t725 * g(1) - t723 * g(2);
t702 = -t727 * pkin(1) + qJDD(1) * qJ(2) + t707;
t747 = qJD(1) * qJD(2);
t744 = -t720 * g(3) - 0.2e1 * t719 * t747;
t763 = pkin(2) * t720;
t675 = (-pkin(6) * qJDD(1) + t727 * t763 - t702) * t719 + t744;
t691 = -t719 * g(3) + (t702 + 0.2e1 * t747) * t720;
t746 = qJDD(1) * t720;
t717 = t720 ^ 2;
t756 = t717 * t727;
t676 = -pkin(2) * t756 + pkin(6) * t746 + t691;
t649 = t722 * t675 + t724 * t676;
t686 = t700 * pkin(3) - t701 * pkin(7);
t726 = qJD(3) ^ 2;
t645 = -t726 * pkin(3) + qJDD(3) * pkin(7) - t700 * t686 + t649;
t716 = t719 ^ 2;
t706 = t723 * g(1) - t725 * g(2);
t739 = qJDD(2) - t706;
t687 = (-pkin(1) - t763) * qJDD(1) + (-qJ(2) + (-t716 - t717) * pkin(6)) * t727 + t739;
t647 = (-t689 + t749) * pkin(7) + (-t688 + t748) * pkin(3) + t687;
t641 = -t721 * t645 + t764 * t647;
t664 = t692 * pkin(4) - t693 * qJ(5);
t685 = qJDD(4) - t688;
t698 = qJD(4) + t700;
t697 = t698 ^ 2;
t639 = -t685 * pkin(4) - t697 * qJ(5) + t693 * t664 + qJDD(5) - t641;
t669 = -t692 * mrSges(6,2) + t698 * mrSges(6,3);
t740 = -m(6) * t639 + t685 * mrSges(6,1) + t698 * t669;
t636 = t661 * mrSges(6,2) + t693 * t665 - t740;
t642 = t764 * t645 + t721 * t647;
t638 = -t697 * pkin(4) + t685 * qJ(5) + 0.2e1 * qJD(5) * t698 - t692 * t664 + t642;
t660 = t693 * qJD(4) - t764 * qJDD(3) + t721 * t689;
t672 = -t698 * mrSges(6,1) + t693 * mrSges(6,2);
t745 = m(6) * t638 + t685 * mrSges(6,3) + t698 * t672;
t752 = t761 * t692 - t768 * t693 + t760 * t698;
t753 = t767 * t692 - t761 * t693 - t759 * t698;
t765 = -t759 * t660 - t760 * t661 - t766 * t685 - t752 * t692 - t753 * t693 + mrSges(5,1) * t641 - mrSges(6,1) * t639 - mrSges(5,2) * t642 + mrSges(6,3) * t638 - pkin(4) * t636 + qJ(5) * (-t660 * mrSges(6,2) - t692 * t665 + t745);
t762 = -mrSges(5,3) - mrSges(6,2);
t757 = mrSges(3,2) * t719;
t671 = t698 * mrSges(5,1) - t693 * mrSges(5,3);
t751 = -t692 * mrSges(5,1) - t693 * mrSges(5,2) - t665;
t630 = m(5) * t642 - t685 * mrSges(5,2) + t762 * t660 - t698 * t671 + t751 * t692 + t745;
t670 = -t698 * mrSges(5,2) - t692 * mrSges(5,3);
t632 = m(5) * t641 + t685 * mrSges(5,1) + t762 * t661 + t698 * t670 + t751 * t693 + t740;
t625 = t764 * t630 - t721 * t632;
t683 = t700 * mrSges(4,1) + t701 * mrSges(4,2);
t695 = qJD(3) * mrSges(4,1) - t701 * mrSges(4,3);
t623 = m(4) * t649 - qJDD(3) * mrSges(4,2) + t688 * mrSges(4,3) - qJD(3) * t695 - t700 * t683 + t625;
t648 = t724 * t675 - t722 * t676;
t644 = -qJDD(3) * pkin(3) - t726 * pkin(7) + t701 * t686 - t648;
t640 = -0.2e1 * qJD(5) * t693 + (t692 * t698 - t661) * qJ(5) + (t693 * t698 + t660) * pkin(4) + t644;
t635 = m(6) * t640 + t660 * mrSges(6,1) - t661 * mrSges(6,3) + t692 * t669 - t693 * t672;
t633 = -m(5) * t644 - t660 * mrSges(5,1) - t661 * mrSges(5,2) - t692 * t670 - t693 * t671 - t635;
t694 = -qJD(3) * mrSges(4,2) - t700 * mrSges(4,3);
t627 = m(4) * t648 + qJDD(3) * mrSges(4,1) - t689 * mrSges(4,3) + qJD(3) * t694 - t701 * t683 + t633;
t614 = t722 * t623 + t724 * t627;
t690 = -t719 * t702 + t744;
t733 = mrSges(3,3) * qJDD(1) + t727 * (-mrSges(3,1) * t720 + t757);
t612 = m(3) * t690 - t733 * t719 + t614;
t741 = t724 * t623 - t722 * t627;
t613 = m(3) * t691 + t733 * t720 + t741;
t742 = -t719 * t612 + t720 * t613;
t605 = m(2) * t707 - t727 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t742;
t699 = -qJDD(1) * pkin(1) - t727 * qJ(2) + t739;
t624 = t721 * t630 + t764 * t632;
t731 = m(4) * t687 - t688 * mrSges(4,1) + t689 * mrSges(4,2) + t700 * t694 + t701 * t695 + t624;
t729 = -m(3) * t699 + mrSges(3,1) * t746 - t731 + (t716 * t727 + t756) * mrSges(3,3);
t616 = (mrSges(2,1) - t757) * qJDD(1) + t729 - t727 * mrSges(2,2) + m(2) * t706;
t755 = t723 * t605 + t725 * t616;
t607 = t720 * t612 + t719 * t613;
t754 = t759 * t692 + t760 * t693 + t766 * t698;
t736 = Ifges(3,5) * t719 + Ifges(3,6) * t720;
t750 = t727 * t736;
t743 = t725 * t605 - t723 * t616;
t738 = Ifges(3,1) * t719 + Ifges(3,4) * t720;
t737 = Ifges(3,4) * t719 + Ifges(3,2) * t720;
t617 = -mrSges(5,1) * t644 - mrSges(6,1) * t640 + mrSges(6,2) * t638 + mrSges(5,3) * t642 - pkin(4) * t635 - t767 * t660 + t761 * t661 + t759 * t685 + t754 * t693 - t752 * t698;
t620 = mrSges(5,2) * t644 + mrSges(6,2) * t639 - mrSges(5,3) * t641 - mrSges(6,3) * t640 - qJ(5) * t635 - t761 * t660 + t768 * t661 - t760 * t685 + t754 * t692 + t753 * t698;
t677 = Ifges(4,5) * t701 - Ifges(4,6) * t700 + Ifges(4,3) * qJD(3);
t678 = Ifges(4,4) * t701 - Ifges(4,2) * t700 + Ifges(4,6) * qJD(3);
t602 = mrSges(4,2) * t687 - mrSges(4,3) * t648 + Ifges(4,1) * t689 + Ifges(4,4) * t688 + Ifges(4,5) * qJDD(3) - pkin(7) * t624 - qJD(3) * t678 - t721 * t617 + t764 * t620 - t700 * t677;
t679 = Ifges(4,1) * t701 - Ifges(4,4) * t700 + Ifges(4,5) * qJD(3);
t608 = -mrSges(4,1) * t687 + mrSges(4,3) * t649 + Ifges(4,4) * t689 + Ifges(4,2) * t688 + Ifges(4,6) * qJDD(3) - pkin(3) * t624 + qJD(3) * t679 - t701 * t677 - t765;
t599 = -mrSges(3,1) * t699 + mrSges(3,3) * t691 - pkin(2) * t731 + pkin(6) * t741 + t737 * qJDD(1) + t722 * t602 + t724 * t608 - t719 * t750;
t601 = mrSges(3,2) * t699 - mrSges(3,3) * t690 - pkin(6) * t614 + t738 * qJDD(1) + t724 * t602 - t722 * t608 + t720 * t750;
t619 = qJDD(1) * t757 - t729;
t732 = mrSges(2,1) * t706 - mrSges(2,2) * t707 + Ifges(2,3) * qJDD(1) - pkin(1) * t619 + qJ(2) * t742 + t720 * t599 + t719 * t601;
t728 = mrSges(4,1) * t648 - mrSges(4,2) * t649 + Ifges(4,5) * t689 + Ifges(4,6) * t688 + Ifges(4,3) * qJDD(3) + pkin(3) * t633 + pkin(7) * t625 + t764 * t617 + t721 * t620 + t701 * t678 + t700 * t679;
t597 = -pkin(1) * t607 - pkin(2) * t614 + mrSges(2,1) * g(3) + (Ifges(2,6) - t736) * qJDD(1) + mrSges(2,3) * t707 - mrSges(3,1) * t690 + mrSges(3,2) * t691 - t728 + (-t719 * t737 + t720 * t738 + Ifges(2,5)) * t727;
t596 = -mrSges(2,2) * g(3) - mrSges(2,3) * t706 + Ifges(2,5) * qJDD(1) - t727 * Ifges(2,6) - qJ(2) * t607 - t719 * t599 + t720 * t601;
t1 = [-m(1) * g(1) + t743; -m(1) * g(2) + t755; (-m(1) - m(2)) * g(3) + t607; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t755 + t725 * t596 - t723 * t597; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t743 + t723 * t596 + t725 * t597; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t732; t732; t619; t728; t765; t636;];
tauJB = t1;
