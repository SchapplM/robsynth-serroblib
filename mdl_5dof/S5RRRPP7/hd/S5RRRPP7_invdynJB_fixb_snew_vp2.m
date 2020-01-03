% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:04:12
% EndTime: 2019-12-31 21:04:15
% DurationCPUTime: 2.28s
% Computational Cost: add. (20170->266), mult. (39299->311), div. (0->0), fcn. (23045->6), ass. (0->103)
t767 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t750 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t749 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t766 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t748 = -Ifges(5,6) + Ifges(6,6) + Ifges(4,6);
t765 = Ifges(6,3) + Ifges(4,3) + Ifges(5,2);
t726 = sin(qJ(1));
t728 = cos(qJ(1));
t713 = t726 * g(1) - t728 * g(2);
t730 = qJD(1) ^ 2;
t694 = -qJDD(1) * pkin(1) - t730 * pkin(6) - t713;
t725 = sin(qJ(2));
t727 = cos(qJ(2));
t751 = qJD(1) * qJD(2);
t742 = t727 * t751;
t707 = t725 * qJDD(1) + t742;
t718 = t725 * t751;
t708 = t727 * qJDD(1) - t718;
t643 = (-t707 - t742) * pkin(7) + (-t708 + t718) * pkin(2) + t694;
t714 = -t728 * g(1) - t726 * g(2);
t695 = -t730 * pkin(1) + qJDD(1) * pkin(6) + t714;
t687 = -t725 * g(3) + t727 * t695;
t706 = (-pkin(2) * t727 - pkin(7) * t725) * qJD(1);
t729 = qJD(2) ^ 2;
t752 = t727 * qJD(1);
t647 = -t729 * pkin(2) + qJDD(2) * pkin(7) + t706 * t752 + t687;
t724 = sin(qJ(3));
t759 = cos(qJ(3));
t640 = t759 * t643 - t724 * t647;
t753 = qJD(1) * t725;
t703 = -t759 * qJD(2) + t724 * t753;
t704 = t724 * qJD(2) + t759 * t753;
t675 = t703 * pkin(3) - t704 * qJ(4);
t702 = -qJDD(3) + t708;
t716 = -qJD(3) + t752;
t715 = t716 ^ 2;
t639 = t702 * pkin(3) - t715 * qJ(4) + t704 * t675 + qJDD(4) - t640;
t685 = -t703 * mrSges(5,2) - t716 * mrSges(5,3);
t764 = -m(5) * t639 - t702 * mrSges(5,1) - t716 * t685;
t672 = -t703 * qJD(3) + t724 * qJDD(2) + t759 * t707;
t686 = -t727 * g(3) - t725 * t695;
t735 = qJDD(2) * pkin(2) + t729 * pkin(7) - t706 * t753 + t686;
t756 = t703 * t716;
t763 = -(t672 + t756) * qJ(4) - t735;
t671 = t704 * qJD(3) - t759 * qJDD(2) + t724 * t707;
t681 = t716 * pkin(4) - t704 * qJ(5);
t701 = t703 ^ 2;
t635 = -t701 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t671 + (pkin(3) * t716 + (2 * qJD(4)) + t681) * t704 - t763;
t679 = -t716 * mrSges(6,2) + t703 * mrSges(6,3);
t682 = t716 * mrSges(6,1) - t704 * mrSges(6,3);
t630 = m(6) * t635 - t671 * mrSges(6,1) + t672 * mrSges(6,2) - t703 * t679 + t704 * t682;
t760 = -2 * qJD(4);
t638 = t704 * t760 + (-t704 * t716 + t671) * pkin(3) + t763;
t684 = t716 * mrSges(5,1) + t704 * mrSges(5,2);
t626 = m(5) * t638 + t671 * mrSges(5,1) - t672 * mrSges(5,3) - t704 * t684 + t703 * t685 - t630;
t641 = t724 * t643 + t759 * t647;
t637 = -t715 * pkin(3) - t702 * qJ(4) - t703 * t675 + t716 * t760 + t641;
t633 = -t701 * pkin(4) + t671 * qJ(5) + 0.2e1 * qJD(5) * t703 - t716 * t681 + t637;
t743 = -t750 * t703 + t767 * t704 - t749 * t716;
t745 = t748 * t703 - t749 * t704 + t765 * t716;
t677 = -t703 * mrSges(6,1) + t704 * mrSges(6,2);
t746 = m(6) * t633 + t671 * mrSges(6,3) + t703 * t677;
t605 = mrSges(4,1) * t735 + mrSges(4,3) * t641 - mrSges(5,1) * t638 + mrSges(5,2) * t637 + mrSges(6,1) * t635 - mrSges(6,3) * t633 + pkin(4) * t630 - qJ(5) * t746 - pkin(3) * t626 + (qJ(5) * t682 - t743) * t716 + t745 * t704 + (qJ(5) * mrSges(6,2) - t748) * t702 + t750 * t672 + t766 * t671;
t631 = -0.2e1 * qJD(5) * t704 + (-t672 + t756) * qJ(5) + (t703 * t704 + t702) * pkin(4) + t639;
t738 = -m(6) * t631 + t672 * mrSges(6,3) + t704 * t677;
t629 = t702 * mrSges(6,1) + t716 * t679 - t738;
t744 = t766 * t703 + t750 * t704 - t748 * t716;
t611 = -mrSges(4,2) * t735 + mrSges(5,2) * t639 + mrSges(6,2) * t635 - mrSges(4,3) * t640 - mrSges(5,3) * t638 - mrSges(6,3) * t631 - qJ(4) * t626 - qJ(5) * t629 - t750 * t671 + t767 * t672 - t749 * t702 + t745 * t703 + t744 * t716;
t680 = t716 * mrSges(4,2) - t703 * mrSges(4,3);
t676 = t703 * mrSges(5,1) - t704 * mrSges(5,3);
t754 = -t703 * mrSges(4,1) - t704 * mrSges(4,2) - t676;
t757 = -mrSges(4,3) - mrSges(5,2);
t623 = m(4) * t640 + (-t679 - t680) * t716 + t754 * t704 + (-mrSges(4,1) - mrSges(6,1)) * t702 + t757 * t672 + t738 + t764;
t683 = -t716 * mrSges(4,1) - t704 * mrSges(4,3);
t736 = m(5) * t637 - t702 * mrSges(5,3) - t716 * t684 + t746;
t624 = m(4) * t641 + (-t682 + t683) * t716 + t754 * t703 + (mrSges(4,2) - mrSges(6,2)) * t702 + t757 * t671 + t736;
t619 = -t724 * t623 + t759 * t624;
t625 = m(4) * t735 - t671 * mrSges(4,1) - t672 * mrSges(4,2) - t703 * t680 - t704 * t683 - t626;
t692 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t725 + Ifges(3,2) * t727) * qJD(1);
t693 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t725 + Ifges(3,4) * t727) * qJD(1);
t762 = mrSges(3,1) * t686 - mrSges(3,2) * t687 + Ifges(3,5) * t707 + Ifges(3,6) * t708 + Ifges(3,3) * qJDD(2) + pkin(2) * t625 + pkin(7) * t619 + (t692 * t725 - t693 * t727) * qJD(1) + t759 * t605 + t724 * t611;
t627 = t672 * mrSges(5,2) + t704 * t676 + t629 - t764;
t761 = -t748 * t671 + t749 * t672 - t765 * t702 + t743 * t703 + t744 * t704 + mrSges(4,1) * t640 - mrSges(5,1) * t639 - mrSges(6,1) * t631 - mrSges(4,2) * t641 + mrSges(6,2) * t633 + mrSges(5,3) * t637 - pkin(3) * t627 - pkin(4) * t629 + qJ(4) * (-t671 * mrSges(5,2) - t702 * mrSges(6,2) - t703 * t676 - t716 * t682 + t736);
t705 = (-mrSges(3,1) * t727 + mrSges(3,2) * t725) * qJD(1);
t710 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t753;
t617 = m(3) * t687 - qJDD(2) * mrSges(3,2) + t708 * mrSges(3,3) - qJD(2) * t710 + t705 * t752 + t619;
t711 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t752;
t621 = m(3) * t686 + qJDD(2) * mrSges(3,1) - t707 * mrSges(3,3) + qJD(2) * t711 - t705 * t753 + t625;
t740 = t727 * t617 - t725 * t621;
t608 = m(2) * t714 - t730 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t740;
t618 = t759 * t623 + t724 * t624;
t733 = -m(3) * t694 + t708 * mrSges(3,1) - t707 * mrSges(3,2) - t710 * t753 + t711 * t752 - t618;
t613 = m(2) * t713 + qJDD(1) * mrSges(2,1) - t730 * mrSges(2,2) + t733;
t755 = t726 * t608 + t728 * t613;
t610 = t725 * t617 + t727 * t621;
t741 = t728 * t608 - t726 * t613;
t691 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t725 + Ifges(3,6) * t727) * qJD(1);
t602 = mrSges(3,2) * t694 - mrSges(3,3) * t686 + Ifges(3,1) * t707 + Ifges(3,4) * t708 + Ifges(3,5) * qJDD(2) - pkin(7) * t618 - qJD(2) * t692 - t724 * t605 + t759 * t611 + t691 * t752;
t604 = -mrSges(3,1) * t694 + mrSges(3,3) * t687 + Ifges(3,4) * t707 + Ifges(3,2) * t708 + Ifges(3,6) * qJDD(2) - pkin(2) * t618 + qJD(2) * t693 - t691 * t753 - t761;
t734 = mrSges(2,1) * t713 - mrSges(2,2) * t714 + Ifges(2,3) * qJDD(1) + pkin(1) * t733 + pkin(6) * t740 + t725 * t602 + t727 * t604;
t600 = mrSges(2,1) * g(3) + mrSges(2,3) * t714 + t730 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t610 - t762;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t713 + Ifges(2,5) * qJDD(1) - t730 * Ifges(2,6) - pkin(6) * t610 + t727 * t602 - t725 * t604;
t1 = [-m(1) * g(1) + t741; -m(1) * g(2) + t755; (-m(1) - m(2)) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t755 + t728 * t599 - t726 * t600; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t741 + t726 * t599 + t728 * t600; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t734; t734; t762; t761; t627; t630;];
tauJB = t1;
