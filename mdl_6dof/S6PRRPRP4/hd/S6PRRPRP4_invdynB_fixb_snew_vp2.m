% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:00:58
% EndTime: 2019-05-05 04:01:05
% DurationCPUTime: 4.72s
% Computational Cost: add. (52824->310), mult. (103521->360), div. (0->0), fcn. (61536->10), ass. (0->135)
t770 = -2 * qJD(4);
t769 = Ifges(4,1) + Ifges(5,2);
t768 = Ifges(6,1) + Ifges(7,1);
t757 = Ifges(4,4) + Ifges(5,6);
t756 = Ifges(6,4) + Ifges(7,4);
t755 = Ifges(4,5) - Ifges(5,4);
t754 = Ifges(6,5) + Ifges(7,5);
t767 = Ifges(4,2) + Ifges(5,3);
t766 = Ifges(6,2) + Ifges(7,2);
t753 = Ifges(4,6) - Ifges(5,5);
t752 = Ifges(6,6) + Ifges(7,6);
t765 = Ifges(4,3) + Ifges(5,1);
t764 = Ifges(6,3) + Ifges(7,3);
t708 = sin(qJ(5));
t711 = cos(qJ(5));
t712 = cos(qJ(3));
t737 = qJD(2) * t712;
t679 = t711 * qJD(3) - t708 * t737;
t709 = sin(qJ(3));
t735 = qJD(2) * qJD(3);
t730 = t709 * t735;
t684 = t712 * qJDD(2) - t730;
t646 = -t679 * qJD(5) - t708 * qJDD(3) - t711 * t684;
t763 = (mrSges(6,1) + mrSges(7,1)) * t646;
t704 = sin(pkin(10));
t706 = cos(pkin(10));
t687 = t704 * g(1) - t706 * g(2);
t688 = -t706 * g(1) - t704 * g(2);
t703 = -g(3) + qJDD(1);
t713 = cos(qJ(2));
t707 = cos(pkin(6));
t710 = sin(qJ(2));
t748 = t707 * t710;
t705 = sin(pkin(6));
t749 = t705 * t710;
t634 = t687 * t748 + t713 * t688 + t703 * t749;
t715 = qJD(2) ^ 2;
t631 = -t715 * pkin(2) + qJDD(2) * pkin(8) + t634;
t657 = -t705 * t687 + t707 * t703;
t627 = t712 * t631 + t709 * t657;
t680 = (-pkin(3) * t712 - qJ(4) * t709) * qJD(2);
t714 = qJD(3) ^ 2;
t623 = t714 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t770 - t680 * t737 - t627;
t633 = -t710 * t688 + (t687 * t707 + t703 * t705) * t713;
t678 = -t708 * qJD(3) - t711 * t737;
t736 = t709 * qJD(2);
t696 = qJD(5) + t736;
t652 = -t696 * mrSges(7,2) + t678 * mrSges(7,3);
t653 = -t696 * mrSges(6,2) + t678 * mrSges(6,3);
t762 = -t763 - (t652 + t653) * t678;
t761 = -pkin(3) - pkin(9);
t760 = pkin(9) * t715;
t759 = t715 * pkin(8);
t719 = -qJDD(2) * pkin(2) - t633;
t630 = t719 - t759;
t729 = t712 * t735;
t683 = t709 * qJDD(2) + t729;
t689 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t736;
t690 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t737;
t691 = -mrSges(5,1) * t737 - qJD(3) * mrSges(5,3);
t717 = pkin(3) * t730 + t736 * t770 + (-t683 - t729) * qJ(4) + t719;
t625 = -t684 * pkin(3) + t717 - t759;
t692 = mrSges(5,1) * t736 + qJD(3) * mrSges(5,2);
t628 = t709 * t631;
t725 = -t714 * qJ(4) + t680 * t736 + qJDD(4) + t628;
t620 = t683 * pkin(4) + t761 * qJDD(3) + (-pkin(4) * t735 - t709 * t760 - t657) * t712 + t725;
t694 = pkin(4) * t736 - qJD(3) * pkin(9);
t702 = t712 ^ 2;
t622 = -t694 * t736 + (-pkin(4) * t702 - pkin(8)) * t715 + t761 * t684 + t717;
t612 = t711 * t620 - t708 * t622;
t647 = t678 * qJD(5) + t711 * qJDD(3) - t708 * t684;
t649 = -t678 * mrSges(7,1) + t679 * mrSges(7,2);
t650 = -t678 * mrSges(6,1) + t679 * mrSges(6,2);
t675 = qJDD(5) + t683;
t609 = -0.2e1 * qJD(6) * t679 + (t678 * t696 - t647) * qJ(6) + (t678 * t679 + t675) * pkin(5) + t612;
t733 = m(7) * t609 + t675 * mrSges(7,1) + t696 * t652;
t604 = m(6) * t612 + t675 * mrSges(6,1) + t696 * t653 + (-t649 - t650) * t679 + (-mrSges(6,3) - mrSges(7,3)) * t647 + t733;
t613 = t708 * t620 + t711 * t622;
t655 = t696 * mrSges(7,1) - t679 * mrSges(7,3);
t656 = t696 * mrSges(6,1) - t679 * mrSges(6,3);
t654 = t696 * pkin(5) - t679 * qJ(6);
t674 = t678 ^ 2;
t611 = -t674 * pkin(5) + t646 * qJ(6) + 0.2e1 * qJD(6) * t678 - t696 * t654 + t613;
t732 = m(7) * t611 + t646 * mrSges(7,3) + t678 * t649;
t606 = m(6) * t613 + t646 * mrSges(6,3) + t678 * t650 + (-t655 - t656) * t696 + (-mrSges(6,2) - mrSges(7,2)) * t675 + t732;
t745 = -t708 * t604 + t711 * t606;
t723 = -m(5) * t625 - t684 * mrSges(5,2) + t692 * t736 - t745;
t716 = -m(4) * t630 + t690 * t737 + t684 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t683 + (-t689 * t709 - t691 * t712) * qJD(2) + t723;
t594 = m(3) * t633 + qJDD(2) * mrSges(3,1) - t715 * mrSges(3,2) + t716;
t751 = t594 * t713;
t750 = t678 * t652;
t747 = t712 * t657;
t626 = -t628 + t747;
t681 = (mrSges(5,2) * t712 - mrSges(5,3) * t709) * qJD(2);
t682 = (-mrSges(4,1) * t712 + mrSges(4,2) * t709) * qJD(2);
t599 = t711 * t604 + t708 * t606;
t624 = -qJDD(3) * pkin(3) + t725 - t747;
t720 = -m(5) * t624 - t683 * mrSges(5,1) - t599;
t596 = m(4) * t626 - t683 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t690 - t691) * qJD(3) + (-t681 - t682) * t736 + t720;
t619 = t684 * pkin(4) + qJD(3) * t694 - t702 * t760 - t623;
t615 = -t646 * pkin(5) - t674 * qJ(6) + t679 * t654 + qJDD(6) + t619;
t731 = m(7) * t615 + t647 * mrSges(7,2) + t679 * t655;
t724 = -m(6) * t619 - t647 * mrSges(6,2) - t679 * t656 - t731;
t718 = -m(5) * t623 + qJDD(3) * mrSges(5,3) + qJD(3) * t692 + t681 * t737 - t724;
t602 = (mrSges(4,3) + mrSges(5,1)) * t684 + t682 * t737 + t718 - qJDD(3) * mrSges(4,2) + m(4) * t627 - qJD(3) * t689 + t762;
t727 = -t709 * t596 + t712 * t602;
t587 = m(3) * t634 - t715 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t727;
t590 = t712 * t596 + t709 * t602;
t589 = m(3) * t657 + t590;
t578 = t587 * t748 - t705 * t589 + t707 * t751;
t576 = m(2) * t687 + t578;
t583 = t713 * t587 - t710 * t594;
t582 = m(2) * t688 + t583;
t746 = t706 * t576 + t704 * t582;
t744 = -t752 * t678 - t754 * t679 - t764 * t696;
t743 = t766 * t678 + t756 * t679 + t752 * t696;
t742 = -t756 * t678 - t768 * t679 - t754 * t696;
t740 = t765 * qJD(3) + (t755 * t709 + t753 * t712) * qJD(2);
t739 = -t753 * qJD(3) + (-t757 * t709 - t767 * t712) * qJD(2);
t738 = t755 * qJD(3) + (t769 * t709 + t757 * t712) * qJD(2);
t577 = t587 * t749 + t707 * t589 + t705 * t751;
t728 = -t704 * t576 + t706 * t582;
t591 = -mrSges(6,1) * t619 + mrSges(6,3) * t613 - mrSges(7,1) * t615 + mrSges(7,3) * t611 - pkin(5) * (t731 - t750) + qJ(6) * t732 + (-qJ(6) * t655 - t742) * t696 + t744 * t679 + (-qJ(6) * mrSges(7,2) + t752) * t675 + t756 * t647 + (pkin(5) * mrSges(7,1) + t766) * t646;
t597 = -t683 * mrSges(5,3) + t691 * t737 - t723;
t607 = -t647 * mrSges(7,3) - t679 * t649 + t733;
t598 = mrSges(6,2) * t619 + mrSges(7,2) * t615 - mrSges(6,3) * t612 - mrSges(7,3) * t609 - qJ(6) * t607 + t756 * t646 + t768 * t647 + t754 * t675 - t744 * t678 - t743 * t696;
t574 = -mrSges(4,1) * t630 + mrSges(4,3) * t627 - mrSges(5,1) * t623 + mrSges(5,2) * t625 - t708 * t598 - t711 * t591 - pkin(4) * (t724 - t762) - pkin(9) * t745 - pkin(3) * t597 + t767 * t684 + t757 * t683 + t753 * qJDD(3) + t738 * qJD(3) - t740 * t736;
t579 = mrSges(5,1) * t624 + mrSges(6,1) * t612 + mrSges(7,1) * t609 + mrSges(4,2) * t630 - mrSges(6,2) * t613 - mrSges(7,2) * t611 - mrSges(4,3) * t626 - mrSges(5,3) * t625 + pkin(4) * t599 + pkin(5) * t607 - qJ(4) * t597 + t757 * t684 + t769 * t683 + t743 * t679 + t742 * t678 + t764 * t675 + t754 * t647 + t752 * t646 + t755 * qJDD(3) + t739 * qJD(3) + t740 * t737;
t572 = mrSges(3,2) * t657 - mrSges(3,3) * t633 + Ifges(3,5) * qJDD(2) - t715 * Ifges(3,6) - pkin(8) * t590 - t709 * t574 + t712 * t579;
t573 = t715 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t590 + mrSges(3,3) * t634 - mrSges(3,1) * t657 - pkin(3) * (-qJD(3) * t691 + t720) - qJ(4) * (-t678 * t653 + t718 - t750 - t763) - mrSges(5,2) * t624 + mrSges(5,3) * t623 - t711 * t598 + t708 * t591 + pkin(9) * t599 - mrSges(4,1) * t626 + mrSges(4,2) * t627 + (-qJ(4) * mrSges(5,1) - t753) * t684 - t755 * t683 + (pkin(3) * mrSges(5,2) - t765) * qJDD(3) + (t738 * t712 + (pkin(3) * t681 + t739) * t709) * qJD(2);
t721 = pkin(7) * t583 + t572 * t710 + t573 * t713;
t571 = mrSges(3,1) * t633 - mrSges(3,2) * t634 + Ifges(3,3) * qJDD(2) + pkin(2) * t716 + pkin(8) * t727 + t712 * t574 + t709 * t579;
t570 = mrSges(2,2) * t703 - mrSges(2,3) * t687 + t713 * t572 - t710 * t573 + (-t577 * t705 - t578 * t707) * pkin(7);
t569 = -mrSges(2,1) * t703 + mrSges(2,3) * t688 - pkin(1) * t577 - t705 * t571 + t721 * t707;
t1 = [-m(1) * g(1) + t728; -m(1) * g(2) + t746; -m(1) * g(3) + m(2) * t703 + t577; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t746 - t704 * t569 + t706 * t570; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t728 + t706 * t569 + t704 * t570; -mrSges(1,1) * g(2) + mrSges(2,1) * t687 + mrSges(1,2) * g(1) - mrSges(2,2) * t688 + pkin(1) * t578 + t707 * t571 + t721 * t705;];
tauB  = t1;
