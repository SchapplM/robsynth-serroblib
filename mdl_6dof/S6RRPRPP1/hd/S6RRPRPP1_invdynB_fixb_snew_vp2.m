% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-05-06 12:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:16:54
% EndTime: 2019-05-06 12:17:06
% DurationCPUTime: 11.28s
% Computational Cost: add. (162989->363), mult. (376232->448), div. (0->0), fcn. (266848->10), ass. (0->140)
t754 = -2 * qJD(3);
t753 = Ifges(6,1) + Ifges(7,1);
t747 = Ifges(6,4) - Ifges(7,5);
t746 = Ifges(6,5) + Ifges(7,4);
t752 = Ifges(6,2) + Ifges(7,3);
t751 = -Ifges(7,2) - Ifges(6,3);
t745 = Ifges(6,6) - Ifges(7,6);
t712 = sin(qJ(2));
t715 = cos(qJ(2));
t734 = qJD(1) * qJD(2);
t699 = t712 * qJDD(1) + t715 * t734;
t713 = sin(qJ(1));
t716 = cos(qJ(1));
t705 = -t716 * g(1) - t713 * g(2);
t718 = qJD(1) ^ 2;
t694 = -t718 * pkin(1) + qJDD(1) * pkin(7) + t705;
t743 = t712 * t694;
t749 = pkin(2) * t718;
t655 = qJDD(2) * pkin(2) - t699 * qJ(3) - t743 + (qJ(3) * t734 + t712 * t749 - g(3)) * t715;
t680 = -t712 * g(3) + t715 * t694;
t700 = t715 * qJDD(1) - t712 * t734;
t737 = qJD(1) * t712;
t701 = qJD(2) * pkin(2) - qJ(3) * t737;
t707 = t715 ^ 2;
t656 = t700 * qJ(3) - qJD(2) * t701 - t707 * t749 + t680;
t709 = sin(pkin(9));
t710 = cos(pkin(9));
t689 = (t709 * t715 + t710 * t712) * qJD(1);
t629 = t710 * t655 - t709 * t656 + t689 * t754;
t688 = (t709 * t712 - t710 * t715) * qJD(1);
t750 = -2 * qJD(5);
t748 = -mrSges(6,3) - mrSges(7,2);
t744 = cos(pkin(10));
t630 = t709 * t655 + t710 * t656 + t688 * t754;
t667 = t688 * mrSges(4,1) + t689 * mrSges(4,2);
t673 = -t709 * t699 + t710 * t700;
t682 = qJD(2) * mrSges(4,1) - t689 * mrSges(4,3);
t668 = t688 * pkin(3) - t689 * pkin(8);
t717 = qJD(2) ^ 2;
t611 = -t717 * pkin(3) + qJDD(2) * pkin(8) - t688 * t668 + t630;
t704 = t713 * g(1) - t716 * g(2);
t724 = -qJDD(1) * pkin(1) - t704;
t659 = -t700 * pkin(2) + qJDD(3) + t701 * t737 + (-qJ(3) * t707 - pkin(7)) * t718 + t724;
t674 = t710 * t699 + t709 * t700;
t614 = (qJD(2) * t688 - t674) * pkin(8) + (qJD(2) * t689 - t673) * pkin(3) + t659;
t711 = sin(qJ(4));
t714 = cos(qJ(4));
t607 = -t711 * t611 + t714 * t614;
t677 = t714 * qJD(2) - t711 * t689;
t648 = t677 * qJD(4) + t711 * qJDD(2) + t714 * t674;
t672 = qJDD(4) - t673;
t678 = t711 * qJD(2) + t714 * t689;
t687 = qJD(4) + t688;
t603 = (t677 * t687 - t648) * qJ(5) + (t677 * t678 + t672) * pkin(4) + t607;
t608 = t714 * t611 + t711 * t614;
t647 = -t678 * qJD(4) + t714 * qJDD(2) - t711 * t674;
t661 = t687 * pkin(4) - t678 * qJ(5);
t676 = t677 ^ 2;
t605 = -t676 * pkin(4) + t647 * qJ(5) - t687 * t661 + t608;
t708 = sin(pkin(10));
t653 = -t744 * t677 + t708 * t678;
t599 = t708 * t603 + t744 * t605 + t653 * t750;
t618 = -t744 * t647 + t708 * t648;
t654 = t708 * t677 + t744 * t678;
t638 = t687 * mrSges(6,1) - t654 * mrSges(6,3);
t631 = t653 * pkin(5) - t654 * qJ(6);
t686 = t687 ^ 2;
t596 = -t686 * pkin(5) + t672 * qJ(6) + 0.2e1 * qJD(6) * t687 - t653 * t631 + t599;
t639 = -t687 * mrSges(7,1) + t654 * mrSges(7,2);
t733 = m(7) * t596 + t672 * mrSges(7,3) + t687 * t639;
t632 = t653 * mrSges(7,1) - t654 * mrSges(7,3);
t738 = -t653 * mrSges(6,1) - t654 * mrSges(6,2) - t632;
t589 = m(6) * t599 - t672 * mrSges(6,2) + t748 * t618 - t687 * t638 + t738 * t653 + t733;
t723 = t744 * t603 - t708 * t605;
t598 = t654 * t750 + t723;
t619 = t708 * t647 + t744 * t648;
t637 = -t687 * mrSges(6,2) - t653 * mrSges(6,3);
t597 = -t672 * pkin(5) - t686 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t631) * t654 - t723;
t636 = -t653 * mrSges(7,2) + t687 * mrSges(7,3);
t726 = -m(7) * t597 + t672 * mrSges(7,1) + t687 * t636;
t591 = m(6) * t598 + t672 * mrSges(6,1) + t748 * t619 + t687 * t637 + t738 * t654 + t726;
t584 = t708 * t589 + t744 * t591;
t657 = -t677 * mrSges(5,1) + t678 * mrSges(5,2);
t660 = -t687 * mrSges(5,2) + t677 * mrSges(5,3);
t582 = m(5) * t607 + t672 * mrSges(5,1) - t648 * mrSges(5,3) - t678 * t657 + t687 * t660 + t584;
t662 = t687 * mrSges(5,1) - t678 * mrSges(5,3);
t728 = t744 * t589 - t708 * t591;
t583 = m(5) * t608 - t672 * mrSges(5,2) + t647 * mrSges(5,3) + t677 * t657 - t687 * t662 + t728;
t729 = -t711 * t582 + t714 * t583;
t577 = m(4) * t630 - qJDD(2) * mrSges(4,2) + t673 * mrSges(4,3) - qJD(2) * t682 - t688 * t667 + t729;
t681 = -qJD(2) * mrSges(4,2) - t688 * mrSges(4,3);
t610 = -qJDD(2) * pkin(3) - t717 * pkin(8) + t689 * t668 - t629;
t606 = -t647 * pkin(4) - t676 * qJ(5) + t678 * t661 + qJDD(5) + t610;
t601 = -0.2e1 * qJD(6) * t654 + (t653 * t687 - t619) * qJ(6) + (t654 * t687 + t618) * pkin(5) + t606;
t594 = m(7) * t601 + t618 * mrSges(7,1) - t619 * mrSges(7,3) + t653 * t636 - t654 * t639;
t721 = m(6) * t606 + t618 * mrSges(6,1) + t619 * mrSges(6,2) + t653 * t637 + t654 * t638 + t594;
t719 = -m(5) * t610 + t647 * mrSges(5,1) - t648 * mrSges(5,2) + t677 * t660 - t678 * t662 - t721;
t593 = m(4) * t629 + qJDD(2) * mrSges(4,1) - t674 * mrSges(4,3) + qJD(2) * t681 - t689 * t667 + t719;
t572 = t709 * t577 + t710 * t593;
t679 = -t715 * g(3) - t743;
t698 = (-mrSges(3,1) * t715 + mrSges(3,2) * t712) * qJD(1);
t736 = qJD(1) * t715;
t703 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t736;
t570 = m(3) * t679 + qJDD(2) * mrSges(3,1) - t699 * mrSges(3,3) + qJD(2) * t703 - t698 * t737 + t572;
t702 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t737;
t730 = t710 * t577 - t709 * t593;
t571 = m(3) * t680 - qJDD(2) * mrSges(3,2) + t700 * mrSges(3,3) - qJD(2) * t702 + t698 * t736 + t730;
t731 = -t712 * t570 + t715 * t571;
t562 = m(2) * t705 - t718 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t731;
t693 = -t718 * pkin(7) + t724;
t578 = t714 * t582 + t711 * t583;
t722 = m(4) * t659 - t673 * mrSges(4,1) + t674 * mrSges(4,2) + t688 * t681 + t689 * t682 + t578;
t720 = -m(3) * t693 + t700 * mrSges(3,1) - t699 * mrSges(3,2) - t702 * t737 + t703 * t736 - t722;
t574 = m(2) * t704 + qJDD(1) * mrSges(2,1) - t718 * mrSges(2,2) + t720;
t742 = t713 * t562 + t716 * t574;
t563 = t715 * t570 + t712 * t571;
t741 = t752 * t653 - t747 * t654 - t745 * t687;
t740 = t745 * t653 - t746 * t654 + t751 * t687;
t739 = -t747 * t653 + t753 * t654 + t746 * t687;
t732 = t716 * t562 - t713 * t574;
t692 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t712 + Ifges(3,4) * t715) * qJD(1);
t691 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t712 + Ifges(3,2) * t715) * qJD(1);
t690 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t712 + Ifges(3,6) * t715) * qJD(1);
t665 = Ifges(4,1) * t689 - Ifges(4,4) * t688 + Ifges(4,5) * qJD(2);
t664 = Ifges(4,4) * t689 - Ifges(4,2) * t688 + Ifges(4,6) * qJD(2);
t663 = Ifges(4,5) * t689 - Ifges(4,6) * t688 + Ifges(4,3) * qJD(2);
t642 = Ifges(5,1) * t678 + Ifges(5,4) * t677 + Ifges(5,5) * t687;
t641 = Ifges(5,4) * t678 + Ifges(5,2) * t677 + Ifges(5,6) * t687;
t640 = Ifges(5,5) * t678 + Ifges(5,6) * t677 + Ifges(5,3) * t687;
t586 = mrSges(6,2) * t606 + mrSges(7,2) * t597 - mrSges(6,3) * t598 - mrSges(7,3) * t601 - qJ(6) * t594 - t747 * t618 + t753 * t619 + t740 * t653 + t746 * t672 + t741 * t687;
t585 = -mrSges(6,1) * t606 - mrSges(7,1) * t601 + mrSges(7,2) * t596 + mrSges(6,3) * t599 - pkin(5) * t594 - t752 * t618 + t747 * t619 + t740 * t654 + t745 * t672 + t739 * t687;
t566 = mrSges(5,2) * t610 - mrSges(5,3) * t607 + Ifges(5,1) * t648 + Ifges(5,4) * t647 + Ifges(5,5) * t672 - qJ(5) * t584 - t708 * t585 + t744 * t586 + t677 * t640 - t687 * t641;
t565 = -mrSges(5,1) * t610 + mrSges(5,3) * t608 + Ifges(5,4) * t648 + Ifges(5,2) * t647 + Ifges(5,6) * t672 - pkin(4) * t721 + qJ(5) * t728 + t744 * t585 + t708 * t586 - t678 * t640 + t687 * t642;
t564 = (qJ(6) * t632 - t739) * t653 + (pkin(5) * t632 + t741) * t654 + (-Ifges(5,3) + t751) * t672 - qJ(6) * t733 + Ifges(4,6) * qJDD(2) - t689 * t663 + Ifges(4,4) * t674 + t677 * t642 - t678 * t641 + qJD(2) * t665 + Ifges(4,2) * t673 - mrSges(4,1) * t659 - Ifges(5,6) * t647 - Ifges(5,5) * t648 + mrSges(4,3) * t630 - mrSges(5,1) * t607 + mrSges(5,2) * t608 + mrSges(6,2) * t599 - mrSges(6,1) * t598 + mrSges(7,1) * t597 - mrSges(7,3) * t596 + (qJ(6) * mrSges(7,2) + t745) * t618 + (pkin(5) * mrSges(7,2) - t746) * t619 - pkin(4) * t584 - pkin(3) * t578 - pkin(5) * t726;
t559 = mrSges(4,2) * t659 - mrSges(4,3) * t629 + Ifges(4,1) * t674 + Ifges(4,4) * t673 + Ifges(4,5) * qJDD(2) - pkin(8) * t578 - qJD(2) * t664 - t711 * t565 + t714 * t566 - t688 * t663;
t558 = mrSges(3,2) * t693 - mrSges(3,3) * t679 + Ifges(3,1) * t699 + Ifges(3,4) * t700 + Ifges(3,5) * qJDD(2) - qJ(3) * t572 - qJD(2) * t691 + t710 * t559 - t709 * t564 + t690 * t736;
t557 = -pkin(1) * t563 + mrSges(2,3) * t705 - pkin(2) * t572 - Ifges(3,5) * t699 - Ifges(3,6) * t700 - mrSges(3,1) * t679 + mrSges(3,2) * t680 - t714 * t565 - pkin(3) * t719 - pkin(8) * t729 - mrSges(4,1) * t629 + mrSges(4,2) * t630 - t711 * t566 - Ifges(4,5) * t674 - Ifges(4,6) * t673 + mrSges(2,1) * g(3) + t718 * Ifges(2,5) - t689 * t664 - t688 * t665 + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t712 * t691 + t715 * t692) * qJD(1);
t556 = -mrSges(3,1) * t693 + mrSges(3,3) * t680 + Ifges(3,4) * t699 + Ifges(3,2) * t700 + Ifges(3,6) * qJDD(2) - pkin(2) * t722 + qJ(3) * t730 + qJD(2) * t692 + t709 * t559 + t710 * t564 - t690 * t737;
t555 = -mrSges(2,2) * g(3) - mrSges(2,3) * t704 + Ifges(2,5) * qJDD(1) - t718 * Ifges(2,6) - pkin(7) * t563 - t712 * t556 + t715 * t558;
t1 = [-m(1) * g(1) + t732; -m(1) * g(2) + t742; (-m(1) - m(2)) * g(3) + t563; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t742 + t716 * t555 - t713 * t557; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t732 + t713 * t555 + t716 * t557; -mrSges(1,1) * g(2) + mrSges(2,1) * t704 + mrSges(1,2) * g(1) - mrSges(2,2) * t705 + Ifges(2,3) * qJDD(1) + pkin(1) * t720 + pkin(7) * t731 + t715 * t556 + t712 * t558;];
tauB  = t1;
