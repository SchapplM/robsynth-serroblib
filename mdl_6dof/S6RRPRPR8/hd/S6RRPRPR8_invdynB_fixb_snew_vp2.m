% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:55:30
% EndTime: 2019-05-06 14:55:45
% DurationCPUTime: 9.88s
% Computational Cost: add. (142662->364), mult. (308217->442), div. (0->0), fcn. (216731->10), ass. (0->141)
t769 = Ifges(5,1) + Ifges(6,1);
t762 = Ifges(5,4) - Ifges(6,5);
t761 = Ifges(5,5) + Ifges(6,4);
t768 = Ifges(5,2) + Ifges(6,3);
t760 = Ifges(5,6) - Ifges(6,6);
t767 = -Ifges(5,3) - Ifges(6,2);
t725 = sin(pkin(10));
t726 = cos(pkin(10));
t729 = sin(qJ(2));
t753 = qJD(1) * t729;
t706 = t726 * qJD(2) - t725 * t753;
t707 = t725 * qJD(2) + t726 * t753;
t728 = sin(qJ(4));
t764 = cos(qJ(4));
t677 = -t706 * t764 + t728 * t707;
t732 = cos(qJ(2));
t751 = qJD(1) * qJD(2);
t750 = t732 * t751;
t712 = t729 * qJDD(1) + t750;
t686 = t726 * qJDD(2) - t725 * t712;
t687 = t725 * qJDD(2) + t726 * t712;
t643 = -t677 * qJD(4) + t728 * t686 + t687 * t764;
t752 = t732 * qJD(1);
t688 = -pkin(3) * t752 - t707 * pkin(8);
t705 = t706 ^ 2;
t730 = sin(qJ(1));
t733 = cos(qJ(1));
t718 = -t733 * g(1) - t730 * g(2);
t735 = qJD(1) ^ 2;
t699 = -t735 * pkin(1) + qJDD(1) * pkin(7) + t718;
t681 = -t732 * g(3) - t729 * t699;
t710 = (-pkin(2) * t732 - qJ(3) * t729) * qJD(1);
t734 = qJD(2) ^ 2;
t740 = qJDD(2) * pkin(2) + t734 * qJ(3) - t710 * t753 - qJDD(3) + t681;
t739 = t686 * pkin(3) + t705 * pkin(8) - t707 * t688 + t740;
t721 = qJD(4) - t752;
t759 = t677 * t721;
t766 = (-t643 + t759) * qJ(5) - t739;
t765 = 2 * qJD(5);
t763 = -mrSges(5,3) - mrSges(6,2);
t682 = -t729 * g(3) + t732 * t699;
t711 = (-mrSges(3,1) * t732 + mrSges(3,2) * t729) * qJD(1);
t722 = t729 * t751;
t713 = t732 * qJDD(1) - t722;
t714 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t753;
t717 = t730 * g(1) - t733 * g(2);
t698 = -qJDD(1) * pkin(1) - t735 * pkin(7) - t717;
t660 = (-t712 - t750) * qJ(3) + (-t713 + t722) * pkin(2) + t698;
t665 = -t734 * pkin(2) + qJDD(2) * qJ(3) + t710 * t752 + t682;
t632 = -0.2e1 * qJD(3) * t707 + t726 * t660 - t725 * t665;
t622 = (-t706 * t752 - t687) * pkin(8) + (t706 * t707 - t713) * pkin(3) + t632;
t633 = 0.2e1 * qJD(3) * t706 + t725 * t660 + t726 * t665;
t625 = -t705 * pkin(3) + t686 * pkin(8) + t688 * t752 + t633;
t616 = t728 * t622 + t764 * t625;
t678 = t728 * t706 + t707 * t764;
t642 = t678 * qJD(4) - t686 * t764 + t728 * t687;
t667 = t721 * mrSges(5,1) - t678 * mrSges(5,3);
t709 = qJDD(4) - t713;
t655 = t677 * pkin(4) - t678 * qJ(5);
t720 = t721 ^ 2;
t611 = -t720 * pkin(4) + t709 * qJ(5) - t677 * t655 + t721 * t765 + t616;
t668 = -t721 * mrSges(6,1) + t678 * mrSges(6,2);
t615 = t622 * t764 - t728 * t625;
t612 = -t709 * pkin(4) - t720 * qJ(5) + t678 * t655 + qJDD(5) - t615;
t606 = (-t643 - t759) * pkin(9) + (t677 * t678 - t709) * pkin(5) + t612;
t670 = -t721 * pkin(5) - t678 * pkin(9);
t676 = t677 ^ 2;
t607 = -t676 * pkin(5) + t642 * pkin(9) + t721 * t670 + t611;
t727 = sin(qJ(6));
t731 = cos(qJ(6));
t604 = t731 * t606 - t727 * t607;
t653 = t731 * t677 - t727 * t678;
t620 = t653 * qJD(6) + t727 * t642 + t731 * t643;
t654 = t727 * t677 + t731 * t678;
t631 = -t653 * mrSges(7,1) + t654 * mrSges(7,2);
t719 = qJD(6) - t721;
t637 = -t719 * mrSges(7,2) + t653 * mrSges(7,3);
t703 = qJDD(6) - t709;
t602 = m(7) * t604 + t703 * mrSges(7,1) - t620 * mrSges(7,3) - t654 * t631 + t719 * t637;
t605 = t727 * t606 + t731 * t607;
t619 = -t654 * qJD(6) + t731 * t642 - t727 * t643;
t638 = t719 * mrSges(7,1) - t654 * mrSges(7,3);
t603 = m(7) * t605 - t703 * mrSges(7,2) + t619 * mrSges(7,3) + t653 * t631 - t719 * t638;
t745 = -t727 * t602 + t731 * t603;
t742 = m(6) * t611 + t709 * mrSges(6,3) + t721 * t668 + t745;
t656 = t677 * mrSges(6,1) - t678 * mrSges(6,3);
t754 = -t677 * mrSges(5,1) - t678 * mrSges(5,2) - t656;
t592 = m(5) * t616 - t709 * mrSges(5,2) + t642 * t763 - t721 * t667 + t677 * t754 + t742;
t666 = -t721 * mrSges(5,2) - t677 * mrSges(5,3);
t595 = t731 * t602 + t727 * t603;
t669 = -t677 * mrSges(6,2) + t721 * mrSges(6,3);
t741 = -m(6) * t612 + t709 * mrSges(6,1) + t721 * t669 - t595;
t594 = m(5) * t615 + t709 * mrSges(5,1) + t643 * t763 + t721 * t666 + t678 * t754 + t741;
t589 = t728 * t592 + t764 * t594;
t679 = -t706 * mrSges(4,1) + t707 * mrSges(4,2);
t684 = mrSges(4,2) * t752 + t706 * mrSges(4,3);
t587 = m(4) * t632 - t713 * mrSges(4,1) - t687 * mrSges(4,3) - t707 * t679 - t684 * t752 + t589;
t685 = -mrSges(4,1) * t752 - t707 * mrSges(4,3);
t746 = t764 * t592 - t728 * t594;
t588 = m(4) * t633 + t713 * mrSges(4,2) + t686 * mrSges(4,3) + t706 * t679 + t685 * t752 + t746;
t747 = -t725 * t587 + t726 * t588;
t582 = m(3) * t682 - qJDD(2) * mrSges(3,2) + t713 * mrSges(3,3) - qJD(2) * t714 + t711 * t752 + t747;
t715 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t752;
t614 = -0.2e1 * qJD(5) * t678 + (t678 * t721 + t642) * pkin(4) + t766;
t609 = -t676 * pkin(9) + (-pkin(4) - pkin(5)) * t642 + (-pkin(4) * t721 + t670 + t765) * t678 - t766;
t743 = -m(7) * t609 + t619 * mrSges(7,1) - t620 * mrSges(7,2) + t653 * t637 - t654 * t638;
t600 = m(6) * t614 + t642 * mrSges(6,1) - t643 * mrSges(6,3) - t678 * t668 + t677 * t669 + t743;
t737 = -m(5) * t739 + t642 * mrSges(5,1) + t643 * mrSges(5,2) + t677 * t666 + t678 * t667 + t600;
t736 = m(4) * t740 + t686 * mrSges(4,1) - t687 * mrSges(4,2) + t706 * t684 - t707 * t685 - t737;
t599 = m(3) * t681 + qJDD(2) * mrSges(3,1) - t712 * mrSges(3,3) + qJD(2) * t715 - t711 * t753 + t736;
t748 = t732 * t582 - t729 * t599;
t576 = m(2) * t718 - t735 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t748;
t583 = t726 * t587 + t725 * t588;
t738 = -m(3) * t698 + t713 * mrSges(3,1) - t712 * mrSges(3,2) - t714 * t753 + t715 * t752 - t583;
t579 = m(2) * t717 + qJDD(1) * mrSges(2,1) - t735 * mrSges(2,2) + t738;
t758 = t730 * t576 + t733 * t579;
t577 = t729 * t582 + t732 * t599;
t757 = t768 * t677 - t762 * t678 - t760 * t721;
t756 = t760 * t677 - t761 * t678 + t767 * t721;
t755 = -t762 * t677 + t769 * t678 + t761 * t721;
t749 = t733 * t576 - t730 * t579;
t697 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t729 + Ifges(3,4) * t732) * qJD(1);
t696 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t729 + Ifges(3,2) * t732) * qJD(1);
t695 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t729 + Ifges(3,6) * t732) * qJD(1);
t673 = Ifges(4,1) * t707 + Ifges(4,4) * t706 - Ifges(4,5) * t752;
t672 = Ifges(4,4) * t707 + Ifges(4,2) * t706 - Ifges(4,6) * t752;
t671 = Ifges(4,5) * t707 + Ifges(4,6) * t706 - Ifges(4,3) * t752;
t628 = Ifges(7,1) * t654 + Ifges(7,4) * t653 + Ifges(7,5) * t719;
t627 = Ifges(7,4) * t654 + Ifges(7,2) * t653 + Ifges(7,6) * t719;
t626 = Ifges(7,5) * t654 + Ifges(7,6) * t653 + Ifges(7,3) * t719;
t597 = mrSges(7,2) * t609 - mrSges(7,3) * t604 + Ifges(7,1) * t620 + Ifges(7,4) * t619 + Ifges(7,5) * t703 + t653 * t626 - t719 * t627;
t596 = -mrSges(7,1) * t609 + mrSges(7,3) * t605 + Ifges(7,4) * t620 + Ifges(7,2) * t619 + Ifges(7,6) * t703 - t654 * t626 + t719 * t628;
t585 = -mrSges(5,2) * t739 + mrSges(6,2) * t612 - mrSges(5,3) * t615 - mrSges(6,3) * t614 - pkin(9) * t595 - qJ(5) * t600 - t727 * t596 + t731 * t597 - t762 * t642 + t769 * t643 + t756 * t677 + t761 * t709 + t757 * t721;
t584 = mrSges(5,1) * t739 - mrSges(6,1) * t614 + mrSges(6,2) * t611 + mrSges(5,3) * t616 - pkin(4) * t600 - pkin(5) * t743 - pkin(9) * t745 - t731 * t596 - t727 * t597 - t768 * t642 + t762 * t643 + t756 * t678 + t760 * t709 + t755 * t721;
t573 = -mrSges(4,2) * t740 - mrSges(4,3) * t632 + Ifges(4,1) * t687 + Ifges(4,4) * t686 - Ifges(4,5) * t713 - pkin(8) * t589 - t728 * t584 + t585 * t764 + t706 * t671 + t672 * t752;
t572 = mrSges(4,1) * t740 + mrSges(4,3) * t633 + Ifges(4,4) * t687 + Ifges(4,2) * t686 - Ifges(4,6) * t713 - pkin(3) * t737 + pkin(8) * t746 + t584 * t764 + t728 * t585 - t707 * t671 - t673 * t752;
t571 = -qJ(5) * t742 + (pkin(4) * mrSges(6,2) - t761) * t643 + (qJ(5) * mrSges(6,2) + t760) * t642 - t695 * t753 + (qJ(5) * t656 - t755) * t677 + (pkin(4) * t656 + t757) * t678 + t767 * t709 - pkin(4) * t741 + Ifges(3,6) * qJDD(2) - pkin(2) * t583 + (Ifges(3,2) + Ifges(4,3)) * t713 + t706 * t673 - t707 * t672 + Ifges(3,4) * t712 + qJD(2) * t697 - mrSges(3,1) * t698 + Ifges(7,3) * t703 + mrSges(3,3) * t682 - Ifges(4,6) * t686 - Ifges(4,5) * t687 - t653 * t628 + t654 * t627 - mrSges(4,1) * t632 + mrSges(4,2) * t633 + Ifges(7,6) * t619 + Ifges(7,5) * t620 - mrSges(5,1) * t615 + mrSges(5,2) * t616 + mrSges(6,1) * t612 - mrSges(6,3) * t611 - mrSges(7,2) * t605 + mrSges(7,1) * t604 - pkin(3) * t589 + pkin(5) * t595;
t570 = mrSges(3,2) * t698 - mrSges(3,3) * t681 + Ifges(3,1) * t712 + Ifges(3,4) * t713 + Ifges(3,5) * qJDD(2) - qJ(3) * t583 - qJD(2) * t696 - t725 * t572 + t726 * t573 + t695 * t752;
t569 = Ifges(2,6) * qJDD(1) + t735 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t718 - Ifges(3,5) * t712 - Ifges(3,6) * t713 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t681 + mrSges(3,2) * t682 - t725 * t573 - t726 * t572 - pkin(2) * t736 - qJ(3) * t747 - pkin(1) * t577 + (-t729 * t696 + t732 * t697) * qJD(1);
t568 = -mrSges(2,2) * g(3) - mrSges(2,3) * t717 + Ifges(2,5) * qJDD(1) - t735 * Ifges(2,6) - pkin(7) * t577 + t732 * t570 - t729 * t571;
t1 = [-m(1) * g(1) + t749; -m(1) * g(2) + t758; (-m(1) - m(2)) * g(3) + t577; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t758 + t733 * t568 - t730 * t569; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t749 + t730 * t568 + t733 * t569; -mrSges(1,1) * g(2) + mrSges(2,1) * t717 + mrSges(1,2) * g(1) - mrSges(2,2) * t718 + Ifges(2,3) * qJDD(1) + pkin(1) * t738 + pkin(7) * t748 + t729 * t570 + t732 * t571;];
tauB  = t1;
