% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR13
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:51
% EndTime: 2019-12-31 20:33:00
% DurationCPUTime: 8.86s
% Computational Cost: add. (132763->312), mult. (288848->394), div. (0->0), fcn. (200421->10), ass. (0->124)
t731 = sin(qJ(1));
t735 = cos(qJ(1));
t718 = g(1) * t731 - t735 * g(2);
t737 = qJD(1) ^ 2;
t702 = -qJDD(1) * pkin(1) - pkin(6) * t737 - t718;
t730 = sin(qJ(2));
t734 = cos(qJ(2));
t752 = qJD(1) * qJD(2);
t751 = t734 * t752;
t713 = qJDD(1) * t730 + t751;
t722 = t730 * t752;
t714 = qJDD(1) * t734 - t722;
t670 = (-t713 - t751) * qJ(3) + (-t714 + t722) * pkin(2) + t702;
t719 = -g(1) * t735 - g(2) * t731;
t703 = -pkin(1) * t737 + qJDD(1) * pkin(6) + t719;
t687 = -g(3) * t730 + t734 * t703;
t711 = (-pkin(2) * t734 - qJ(3) * t730) * qJD(1);
t736 = qJD(2) ^ 2;
t753 = qJD(1) * t734;
t673 = -pkin(2) * t736 + qJDD(2) * qJ(3) + t711 * t753 + t687;
t726 = sin(pkin(9));
t727 = cos(pkin(9));
t754 = qJD(1) * t730;
t708 = qJD(2) * t726 + t727 * t754;
t653 = -0.2e1 * qJD(3) * t708 + t727 * t670 - t673 * t726;
t691 = qJDD(2) * t726 + t713 * t727;
t707 = qJD(2) * t727 - t726 * t754;
t645 = (-t707 * t753 - t691) * pkin(7) + (t707 * t708 - t714) * pkin(3) + t653;
t654 = 0.2e1 * qJD(3) * t707 + t726 * t670 + t727 * t673;
t690 = qJDD(2) * t727 - t713 * t726;
t692 = -pkin(3) * t753 - pkin(7) * t708;
t706 = t707 ^ 2;
t647 = -pkin(3) * t706 + pkin(7) * t690 + t692 * t753 + t654;
t729 = sin(qJ(4));
t733 = cos(qJ(4));
t634 = t733 * t645 - t647 * t729;
t683 = t707 * t733 - t708 * t729;
t660 = qJD(4) * t683 + t690 * t729 + t691 * t733;
t684 = t707 * t729 + t708 * t733;
t710 = qJDD(4) - t714;
t721 = qJD(4) - t753;
t632 = (t683 * t721 - t660) * pkin(8) + (t683 * t684 + t710) * pkin(4) + t634;
t635 = t729 * t645 + t733 * t647;
t659 = -qJD(4) * t684 + t690 * t733 - t691 * t729;
t676 = pkin(4) * t721 - pkin(8) * t684;
t682 = t683 ^ 2;
t633 = -pkin(4) * t682 + pkin(8) * t659 - t676 * t721 + t635;
t728 = sin(qJ(5));
t732 = cos(qJ(5));
t631 = t632 * t728 + t633 * t732;
t686 = -t734 * g(3) - t730 * t703;
t672 = -qJDD(2) * pkin(2) - qJ(3) * t736 + t711 * t754 + qJDD(3) - t686;
t655 = -pkin(3) * t690 - pkin(7) * t706 + t708 * t692 + t672;
t636 = -pkin(4) * t659 - pkin(8) * t682 + t676 * t684 + t655;
t666 = t683 * t728 + t684 * t732;
t640 = -qJD(5) * t666 + t659 * t732 - t660 * t728;
t665 = t683 * t732 - t684 * t728;
t641 = qJD(5) * t665 + t659 * t728 + t660 * t732;
t720 = qJD(5) + t721;
t648 = Ifges(6,5) * t666 + Ifges(6,6) * t665 + Ifges(6,3) * t720;
t650 = Ifges(6,1) * t666 + Ifges(6,4) * t665 + Ifges(6,5) * t720;
t704 = qJDD(5) + t710;
t619 = -mrSges(6,1) * t636 + mrSges(6,3) * t631 + Ifges(6,4) * t641 + Ifges(6,2) * t640 + Ifges(6,6) * t704 - t648 * t666 + t650 * t720;
t630 = t632 * t732 - t633 * t728;
t649 = Ifges(6,4) * t666 + Ifges(6,2) * t665 + Ifges(6,6) * t720;
t620 = mrSges(6,2) * t636 - mrSges(6,3) * t630 + Ifges(6,1) * t641 + Ifges(6,4) * t640 + Ifges(6,5) * t704 + t648 * t665 - t649 * t720;
t661 = Ifges(5,5) * t684 + Ifges(5,6) * t683 + Ifges(5,3) * t721;
t663 = Ifges(5,1) * t684 + Ifges(5,4) * t683 + Ifges(5,5) * t721;
t656 = -mrSges(6,2) * t720 + mrSges(6,3) * t665;
t657 = mrSges(6,1) * t720 - mrSges(6,3) * t666;
t743 = m(6) * t636 - t640 * mrSges(6,1) + t641 * mrSges(6,2) - t665 * t656 + t666 * t657;
t652 = -mrSges(6,1) * t665 + mrSges(6,2) * t666;
t626 = m(6) * t630 + mrSges(6,1) * t704 - mrSges(6,3) * t641 - t652 * t666 + t656 * t720;
t627 = m(6) * t631 - mrSges(6,2) * t704 + mrSges(6,3) * t640 + t652 * t665 - t657 * t720;
t747 = -t626 * t728 + t732 * t627;
t607 = -mrSges(5,1) * t655 + mrSges(5,3) * t635 + Ifges(5,4) * t660 + Ifges(5,2) * t659 + Ifges(5,6) * t710 - pkin(4) * t743 + pkin(8) * t747 + t732 * t619 + t728 * t620 - t684 * t661 + t721 * t663;
t618 = t732 * t626 + t728 * t627;
t662 = Ifges(5,4) * t684 + Ifges(5,2) * t683 + Ifges(5,6) * t721;
t608 = mrSges(5,2) * t655 - mrSges(5,3) * t634 + Ifges(5,1) * t660 + Ifges(5,4) * t659 + Ifges(5,5) * t710 - pkin(8) * t618 - t619 * t728 + t620 * t732 + t661 * t683 - t662 * t721;
t677 = Ifges(4,5) * t708 + Ifges(4,6) * t707 - Ifges(4,3) * t753;
t679 = Ifges(4,1) * t708 + Ifges(4,4) * t707 - Ifges(4,5) * t753;
t674 = -mrSges(5,2) * t721 + mrSges(5,3) * t683;
t675 = mrSges(5,1) * t721 - mrSges(5,3) * t684;
t740 = m(5) * t655 - t659 * mrSges(5,1) + t660 * mrSges(5,2) - t683 * t674 + t684 * t675 + t743;
t667 = -mrSges(5,1) * t683 + mrSges(5,2) * t684;
t616 = m(5) * t634 + mrSges(5,1) * t710 - mrSges(5,3) * t660 - t667 * t684 + t674 * t721 + t618;
t617 = m(5) * t635 - mrSges(5,2) * t710 + mrSges(5,3) * t659 + t667 * t683 - t675 * t721 + t747;
t748 = -t729 * t616 + t733 * t617;
t592 = -mrSges(4,1) * t672 + mrSges(4,3) * t654 + Ifges(4,4) * t691 + Ifges(4,2) * t690 - Ifges(4,6) * t714 - pkin(3) * t740 + pkin(7) * t748 + t733 * t607 + t729 * t608 - t708 * t677 - t679 * t753;
t612 = t733 * t616 + t729 * t617;
t678 = Ifges(4,4) * t708 + Ifges(4,2) * t707 - Ifges(4,6) * t753;
t593 = mrSges(4,2) * t672 - mrSges(4,3) * t653 + Ifges(4,1) * t691 + Ifges(4,4) * t690 - Ifges(4,5) * t714 - pkin(7) * t612 - t607 * t729 + t608 * t733 + t677 * t707 + t678 * t753;
t685 = -mrSges(4,1) * t707 + mrSges(4,2) * t708;
t745 = mrSges(4,2) * t753 + mrSges(4,3) * t707;
t610 = m(4) * t653 - t714 * mrSges(4,1) - t691 * mrSges(4,3) - t708 * t685 - t745 * t753 + t612;
t689 = -mrSges(4,1) * t753 - mrSges(4,3) * t708;
t611 = m(4) * t654 + mrSges(4,2) * t714 + mrSges(4,3) * t690 + t685 * t707 + t689 * t753 + t748;
t606 = -t610 * t726 + t727 * t611;
t628 = m(4) * t672 - t690 * mrSges(4,1) + t691 * mrSges(4,2) + t708 * t689 - t707 * t745 + t740;
t700 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t730 + Ifges(3,2) * t734) * qJD(1);
t701 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t730 + Ifges(3,4) * t734) * qJD(1);
t756 = mrSges(3,1) * t686 - mrSges(3,2) * t687 + Ifges(3,5) * t713 + Ifges(3,6) * t714 + Ifges(3,3) * qJDD(2) - pkin(2) * t628 + qJ(3) * t606 + t727 * t592 + t726 * t593 + (t730 * t700 - t734 * t701) * qJD(1);
t712 = (-mrSges(3,1) * t734 + mrSges(3,2) * t730) * qJD(1);
t716 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t754;
t604 = m(3) * t687 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t714 - qJD(2) * t716 + t712 * t753 + t606;
t717 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t753;
t622 = m(3) * t686 + qJDD(2) * mrSges(3,1) - t713 * mrSges(3,3) + qJD(2) * t717 - t712 * t754 - t628;
t749 = t734 * t604 - t622 * t730;
t596 = m(2) * t719 - mrSges(2,1) * t737 - qJDD(1) * mrSges(2,2) + t749;
t605 = t610 * t727 + t611 * t726;
t741 = -m(3) * t702 + t714 * mrSges(3,1) - mrSges(3,2) * t713 - t716 * t754 + t717 * t753 - t605;
t600 = m(2) * t718 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t737 + t741;
t755 = t731 * t596 + t735 * t600;
t598 = t730 * t604 + t734 * t622;
t750 = t735 * t596 - t600 * t731;
t699 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t730 + Ifges(3,6) * t734) * qJD(1);
t589 = mrSges(3,2) * t702 - mrSges(3,3) * t686 + Ifges(3,1) * t713 + Ifges(3,4) * t714 + Ifges(3,5) * qJDD(2) - qJ(3) * t605 - qJD(2) * t700 - t592 * t726 + t593 * t727 + t699 * t753;
t742 = -mrSges(6,1) * t630 + mrSges(6,2) * t631 - Ifges(6,5) * t641 - Ifges(6,6) * t640 - Ifges(6,3) * t704 - t666 * t649 + t665 * t650;
t738 = mrSges(5,1) * t634 - mrSges(5,2) * t635 + Ifges(5,5) * t660 + Ifges(5,6) * t659 + Ifges(5,3) * t710 + pkin(4) * t618 + t684 * t662 - t683 * t663 - t742;
t591 = (Ifges(3,2) + Ifges(4,3)) * t714 + Ifges(3,4) * t713 + t707 * t679 - t708 * t678 + mrSges(3,3) * t687 - Ifges(4,6) * t690 - Ifges(4,5) * t691 + qJD(2) * t701 - mrSges(3,1) * t702 - mrSges(4,1) * t653 + mrSges(4,2) * t654 + Ifges(3,6) * qJDD(2) - pkin(3) * t612 - pkin(2) * t605 - t738 - t699 * t754;
t744 = mrSges(2,1) * t718 - mrSges(2,2) * t719 + Ifges(2,3) * qJDD(1) + pkin(1) * t741 + pkin(6) * t749 + t730 * t589 + t734 * t591;
t587 = mrSges(2,1) * g(3) + mrSges(2,3) * t719 + t737 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t598 - t756;
t586 = -mrSges(2,2) * g(3) - mrSges(2,3) * t718 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t737 - pkin(6) * t598 + t589 * t734 - t591 * t730;
t1 = [-m(1) * g(1) + t750; -m(1) * g(2) + t755; (-m(1) - m(2)) * g(3) + t598; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t755 + t735 * t586 - t731 * t587; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t750 + t731 * t586 + t735 * t587; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t744; t744; t756; t628; t738; -t742;];
tauJB = t1;
