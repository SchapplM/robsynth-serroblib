% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:17:35
% EndTime: 2019-05-06 17:17:52
% DurationCPUTime: 13.71s
% Computational Cost: add. (186893->365), mult. (432604->448), div. (0->0), fcn. (318556->10), ass. (0->140)
t751 = Ifges(6,1) + Ifges(7,1);
t745 = Ifges(6,4) + Ifges(7,4);
t744 = Ifges(6,5) + Ifges(7,5);
t750 = Ifges(6,2) + Ifges(7,2);
t749 = Ifges(6,6) + Ifges(7,6);
t748 = Ifges(6,3) + Ifges(7,3);
t720 = qJD(1) ^ 2;
t747 = pkin(2) * t720;
t746 = -mrSges(6,2) - mrSges(7,2);
t715 = sin(qJ(1));
t719 = cos(qJ(1));
t704 = -g(1) * t719 - g(2) * t715;
t693 = -pkin(1) * t720 + qJDD(1) * pkin(7) + t704;
t714 = sin(qJ(2));
t742 = t714 * t693;
t718 = cos(qJ(2));
t734 = qJD(1) * qJD(2);
t698 = qJDD(1) * t714 + t718 * t734;
t658 = qJDD(2) * pkin(2) - t698 * qJ(3) - t742 + (qJ(3) * t734 + t714 * t747 - g(3)) * t718;
t679 = -g(3) * t714 + t718 * t693;
t699 = qJDD(1) * t718 - t714 * t734;
t736 = qJD(1) * t714;
t700 = qJD(2) * pkin(2) - qJ(3) * t736;
t709 = t718 ^ 2;
t659 = qJ(3) * t699 - qJD(2) * t700 - t709 * t747 + t679;
t710 = sin(pkin(10));
t711 = cos(pkin(10));
t688 = (t710 * t718 + t711 * t714) * qJD(1);
t632 = -0.2e1 * qJD(3) * t688 + t711 * t658 - t710 * t659;
t677 = t698 * t711 + t699 * t710;
t687 = (-t710 * t714 + t711 * t718) * qJD(1);
t611 = (qJD(2) * t687 - t677) * pkin(8) + (t687 * t688 + qJDD(2)) * pkin(3) + t632;
t633 = 0.2e1 * qJD(3) * t687 + t710 * t658 + t711 * t659;
t676 = -t698 * t710 + t699 * t711;
t682 = qJD(2) * pkin(3) - pkin(8) * t688;
t686 = t687 ^ 2;
t614 = -pkin(3) * t686 + pkin(8) * t676 - qJD(2) * t682 + t633;
t713 = sin(qJ(4));
t717 = cos(qJ(4));
t609 = t713 * t611 + t717 * t614;
t671 = t687 * t713 + t688 * t717;
t639 = -qJD(4) * t671 + t676 * t717 - t677 * t713;
t670 = t687 * t717 - t688 * t713;
t653 = -mrSges(5,1) * t670 + mrSges(5,2) * t671;
t708 = qJD(2) + qJD(4);
t665 = mrSges(5,1) * t708 - mrSges(5,3) * t671;
t707 = qJDD(2) + qJDD(4);
t654 = -pkin(4) * t670 - pkin(9) * t671;
t706 = t708 ^ 2;
t604 = -pkin(4) * t706 + pkin(9) * t707 + t654 * t670 + t609;
t703 = t715 * g(1) - t719 * g(2);
t725 = -qJDD(1) * pkin(1) - t703;
t661 = -t699 * pkin(2) + qJDD(3) + t700 * t736 + (-qJ(3) * t709 - pkin(7)) * t720 + t725;
t630 = -pkin(3) * t676 - pkin(8) * t686 + t688 * t682 + t661;
t640 = qJD(4) * t670 + t676 * t713 + t677 * t717;
t607 = (-t670 * t708 - t640) * pkin(9) + (t671 * t708 - t639) * pkin(4) + t630;
t712 = sin(qJ(5));
t716 = cos(qJ(5));
t599 = -t712 * t604 + t716 * t607;
t662 = -t671 * t712 + t708 * t716;
t619 = qJD(5) * t662 + t640 * t716 + t707 * t712;
t638 = qJDD(5) - t639;
t663 = t671 * t716 + t708 * t712;
t641 = -mrSges(7,1) * t662 + mrSges(7,2) * t663;
t642 = -mrSges(6,1) * t662 + mrSges(6,2) * t663;
t666 = qJD(5) - t670;
t644 = -mrSges(6,2) * t666 + mrSges(6,3) * t662;
t596 = -0.2e1 * qJD(6) * t663 + (t662 * t666 - t619) * qJ(6) + (t662 * t663 + t638) * pkin(5) + t599;
t643 = -mrSges(7,2) * t666 + mrSges(7,3) * t662;
t733 = m(7) * t596 + t638 * mrSges(7,1) + t666 * t643;
t588 = m(6) * t599 + t638 * mrSges(6,1) + t666 * t644 + (-t641 - t642) * t663 + (-mrSges(6,3) - mrSges(7,3)) * t619 + t733;
t600 = t716 * t604 + t712 * t607;
t618 = -qJD(5) * t663 - t640 * t712 + t707 * t716;
t645 = pkin(5) * t666 - qJ(6) * t663;
t660 = t662 ^ 2;
t598 = -pkin(5) * t660 + qJ(6) * t618 + 0.2e1 * qJD(6) * t662 - t645 * t666 + t600;
t732 = m(7) * t598 + t618 * mrSges(7,3) + t662 * t641;
t646 = mrSges(7,1) * t666 - mrSges(7,3) * t663;
t737 = -mrSges(6,1) * t666 + mrSges(6,3) * t663 - t646;
t591 = m(6) * t600 + t618 * mrSges(6,3) + t638 * t746 + t662 * t642 + t666 * t737 + t732;
t727 = -t712 * t588 + t716 * t591;
t584 = m(5) * t609 - mrSges(5,2) * t707 + mrSges(5,3) * t639 + t653 * t670 - t665 * t708 + t727;
t608 = t611 * t717 - t713 * t614;
t664 = -mrSges(5,2) * t708 + mrSges(5,3) * t670;
t603 = -pkin(4) * t707 - pkin(9) * t706 + t671 * t654 - t608;
t601 = -pkin(5) * t618 - qJ(6) * t660 + t645 * t663 + qJDD(6) + t603;
t726 = m(7) * t601 - t618 * mrSges(7,1) - t662 * t643;
t722 = -m(6) * t603 + t618 * mrSges(6,1) + t619 * t746 + t662 * t644 + t663 * t737 - t726;
t593 = m(5) * t608 + mrSges(5,1) * t707 - mrSges(5,3) * t640 - t653 * t671 + t664 * t708 + t722;
t578 = t713 * t584 + t717 * t593;
t674 = -mrSges(4,1) * t687 + mrSges(4,2) * t688;
t680 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t687;
t576 = m(4) * t632 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t677 + qJD(2) * t680 - t674 * t688 + t578;
t681 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t688;
t728 = t717 * t584 - t593 * t713;
t577 = m(4) * t633 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t676 - qJD(2) * t681 + t674 * t687 + t728;
t571 = t711 * t576 + t710 * t577;
t678 = -t718 * g(3) - t742;
t697 = (-mrSges(3,1) * t718 + mrSges(3,2) * t714) * qJD(1);
t735 = qJD(1) * t718;
t702 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t735;
t569 = m(3) * t678 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t698 + qJD(2) * t702 - t697 * t736 + t571;
t701 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t736;
t729 = -t576 * t710 + t711 * t577;
t570 = m(3) * t679 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t699 - qJD(2) * t701 + t697 * t735 + t729;
t730 = -t569 * t714 + t718 * t570;
t563 = m(2) * t704 - mrSges(2,1) * t720 - qJDD(1) * mrSges(2,2) + t730;
t692 = -pkin(7) * t720 + t725;
t586 = t716 * t588 + t712 * t591;
t724 = m(5) * t630 - t639 * mrSges(5,1) + t640 * mrSges(5,2) - t670 * t664 + t671 * t665 + t586;
t723 = m(4) * t661 - t676 * mrSges(4,1) + mrSges(4,2) * t677 - t687 * t680 + t681 * t688 + t724;
t721 = -m(3) * t692 + t699 * mrSges(3,1) - mrSges(3,2) * t698 - t701 * t736 + t702 * t735 - t723;
t581 = m(2) * t703 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t720 + t721;
t741 = t715 * t563 + t719 * t581;
t564 = t718 * t569 + t714 * t570;
t740 = t749 * t662 + t744 * t663 + t748 * t666;
t739 = -t750 * t662 - t745 * t663 - t749 * t666;
t738 = t745 * t662 + t751 * t663 + t744 * t666;
t731 = t719 * t563 - t581 * t715;
t691 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t714 + Ifges(3,4) * t718) * qJD(1);
t690 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t714 + Ifges(3,2) * t718) * qJD(1);
t689 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t714 + Ifges(3,6) * t718) * qJD(1);
t669 = Ifges(4,1) * t688 + Ifges(4,4) * t687 + Ifges(4,5) * qJD(2);
t668 = Ifges(4,4) * t688 + Ifges(4,2) * t687 + Ifges(4,6) * qJD(2);
t667 = Ifges(4,5) * t688 + Ifges(4,6) * t687 + Ifges(4,3) * qJD(2);
t650 = Ifges(5,1) * t671 + Ifges(5,4) * t670 + Ifges(5,5) * t708;
t649 = Ifges(5,4) * t671 + Ifges(5,2) * t670 + Ifges(5,6) * t708;
t648 = Ifges(5,5) * t671 + Ifges(5,6) * t670 + Ifges(5,3) * t708;
t594 = -t619 * mrSges(7,3) - t663 * t641 + t733;
t585 = mrSges(6,2) * t603 + mrSges(7,2) * t601 - mrSges(6,3) * t599 - mrSges(7,3) * t596 - qJ(6) * t594 + t745 * t618 + t751 * t619 + t744 * t638 + t740 * t662 + t739 * t666;
t579 = -mrSges(6,1) * t603 + mrSges(6,3) * t600 - mrSges(7,1) * t601 + mrSges(7,3) * t598 - pkin(5) * t726 + qJ(6) * t732 + (-qJ(6) * t646 + t738) * t666 + (-pkin(5) * t646 - t740) * t663 + (-mrSges(7,2) * qJ(6) + t749) * t638 + (-mrSges(7,2) * pkin(5) + t745) * t619 + t750 * t618;
t572 = -mrSges(5,1) * t630 - mrSges(6,1) * t599 - mrSges(7,1) * t596 + mrSges(6,2) * t600 + mrSges(7,2) * t598 + mrSges(5,3) * t609 + Ifges(5,4) * t640 + Ifges(5,2) * t639 + Ifges(5,6) * t707 - pkin(4) * t586 - pkin(5) * t594 - t671 * t648 + t708 * t650 + t739 * t663 + t738 * t662 - t748 * t638 - t744 * t619 - t749 * t618;
t565 = mrSges(5,2) * t630 - mrSges(5,3) * t608 + Ifges(5,1) * t640 + Ifges(5,4) * t639 + Ifges(5,5) * t707 - pkin(9) * t586 - t579 * t712 + t585 * t716 + t648 * t670 - t649 * t708;
t560 = mrSges(4,2) * t661 - mrSges(4,3) * t632 + Ifges(4,1) * t677 + Ifges(4,4) * t676 + Ifges(4,5) * qJDD(2) - pkin(8) * t578 - qJD(2) * t668 + t565 * t717 - t572 * t713 + t667 * t687;
t559 = -mrSges(4,1) * t661 + mrSges(4,3) * t633 + Ifges(4,4) * t677 + Ifges(4,2) * t676 + Ifges(4,6) * qJDD(2) - pkin(3) * t724 + pkin(8) * t728 + qJD(2) * t669 + t713 * t565 + t717 * t572 - t688 * t667;
t558 = -pkin(9) * t727 - pkin(4) * t722 - mrSges(4,1) * t632 + mrSges(2,1) * g(3) - pkin(3) * t578 - Ifges(4,6) * t676 - mrSges(5,1) * t608 + mrSges(5,2) * t609 - pkin(2) * t571 - Ifges(5,5) * t640 - Ifges(4,5) * t677 - mrSges(3,1) * t678 + mrSges(3,2) * t679 + t670 * t650 - t671 * t649 + mrSges(4,2) * t633 - Ifges(5,6) * t639 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t690 * t714 + t691 * t718) * qJD(1) - pkin(1) * t564 + Ifges(2,6) * qJDD(1) + t687 * t669 - t688 * t668 - Ifges(3,5) * t698 - Ifges(3,6) * t699 + mrSges(2,3) * t704 - Ifges(5,3) * t707 - t712 * t585 - t716 * t579 + t720 * Ifges(2,5);
t557 = mrSges(3,2) * t692 - mrSges(3,3) * t678 + Ifges(3,1) * t698 + Ifges(3,4) * t699 + Ifges(3,5) * qJDD(2) - qJ(3) * t571 - qJD(2) * t690 - t559 * t710 + t560 * t711 + t689 * t735;
t556 = -mrSges(3,1) * t692 + mrSges(3,3) * t679 + Ifges(3,4) * t698 + Ifges(3,2) * t699 + Ifges(3,6) * qJDD(2) - pkin(2) * t723 + qJ(3) * t729 + qJD(2) * t691 + t711 * t559 + t710 * t560 - t689 * t736;
t555 = -mrSges(2,2) * g(3) - mrSges(2,3) * t703 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t720 - pkin(7) * t564 - t556 * t714 + t557 * t718;
t1 = [-m(1) * g(1) + t731; -m(1) * g(2) + t741; (-m(1) - m(2)) * g(3) + t564; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t741 + t719 * t555 - t715 * t558; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t731 + t715 * t555 + t719 * t558; -mrSges(1,1) * g(2) + mrSges(2,1) * t703 + mrSges(1,2) * g(1) - mrSges(2,2) * t704 + Ifges(2,3) * qJDD(1) + pkin(1) * t721 + pkin(7) * t730 + t718 * t556 + t714 * t557;];
tauB  = t1;
