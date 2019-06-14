% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:54:53
% EndTime: 2019-05-07 17:55:08
% DurationCPUTime: 12.45s
% Computational Cost: add. (201067->364), mult. (406618->448), div. (0->0), fcn. (292432->10), ass. (0->141)
t752 = Ifges(6,1) + Ifges(7,1);
t746 = Ifges(6,4) - Ifges(7,5);
t745 = Ifges(6,5) + Ifges(7,4);
t751 = Ifges(6,2) + Ifges(7,3);
t750 = -Ifges(7,2) - Ifges(6,3);
t744 = Ifges(6,6) - Ifges(7,6);
t712 = sin(qJ(3));
t713 = sin(qJ(2));
t716 = cos(qJ(3));
t717 = cos(qJ(2));
t690 = (t712 * t713 - t716 * t717) * qJD(1);
t749 = -2 * qJD(5);
t719 = qJD(1) ^ 2;
t748 = pkin(2) * t719;
t747 = -mrSges(6,3) - mrSges(7,2);
t743 = cos(pkin(10));
t714 = sin(qJ(1));
t718 = cos(qJ(1));
t704 = -t718 * g(1) - t714 * g(2);
t693 = -t719 * pkin(1) + qJDD(1) * pkin(7) + t704;
t742 = t713 * t693;
t734 = qJD(1) * qJD(2);
t698 = t713 * qJDD(1) + t717 * t734;
t657 = qJDD(2) * pkin(2) - t698 * pkin(8) - t742 + (pkin(8) * t734 + t713 * t748 - g(3)) * t717;
t681 = -t713 * g(3) + t717 * t693;
t699 = t717 * qJDD(1) - t713 * t734;
t736 = qJD(1) * t713;
t702 = qJD(2) * pkin(2) - pkin(8) * t736;
t709 = t717 ^ 2;
t658 = t699 * pkin(8) - qJD(2) * t702 - t709 * t748 + t681;
t634 = t712 * t657 + t716 * t658;
t691 = (t712 * t717 + t713 * t716) * qJD(1);
t664 = -t691 * qJD(3) - t712 * t698 + t716 * t699;
t674 = t690 * mrSges(4,1) + t691 * mrSges(4,2);
t708 = qJD(2) + qJD(3);
t683 = t708 * mrSges(4,1) - t691 * mrSges(4,3);
t707 = qJDD(2) + qJDD(3);
t665 = -t690 * qJD(3) + t716 * t698 + t712 * t699;
t703 = t714 * g(1) - t718 * g(2);
t725 = -qJDD(1) * pkin(1) - t703;
t666 = -t699 * pkin(2) + t702 * t736 + (-pkin(8) * t709 - pkin(7)) * t719 + t725;
t612 = (t690 * t708 - t665) * pkin(9) + (t691 * t708 - t664) * pkin(3) + t666;
t675 = t690 * pkin(3) - t691 * pkin(9);
t706 = t708 ^ 2;
t620 = -t706 * pkin(3) + t707 * pkin(9) - t690 * t675 + t634;
t711 = sin(qJ(4));
t715 = cos(qJ(4));
t607 = t715 * t612 - t711 * t620;
t678 = -t711 * t691 + t715 * t708;
t640 = t678 * qJD(4) + t715 * t665 + t711 * t707;
t663 = qJDD(4) - t664;
t679 = t715 * t691 + t711 * t708;
t686 = qJD(4) + t690;
t604 = (t678 * t686 - t640) * qJ(5) + (t678 * t679 + t663) * pkin(4) + t607;
t608 = t711 * t612 + t715 * t620;
t639 = -t679 * qJD(4) - t711 * t665 + t715 * t707;
t668 = t686 * pkin(4) - t679 * qJ(5);
t677 = t678 ^ 2;
t606 = -t677 * pkin(4) + t639 * qJ(5) - t686 * t668 + t608;
t710 = sin(pkin(10));
t651 = -t743 * t678 + t710 * t679;
t600 = t710 * t604 + t743 * t606 + t651 * t749;
t616 = -t743 * t639 + t710 * t640;
t652 = t710 * t678 + t743 * t679;
t643 = t686 * mrSges(6,1) - t652 * mrSges(6,3);
t630 = t651 * pkin(5) - t652 * qJ(6);
t685 = t686 ^ 2;
t597 = -t685 * pkin(5) + t663 * qJ(6) + 0.2e1 * qJD(6) * t686 - t651 * t630 + t600;
t644 = -t686 * mrSges(7,1) + t652 * mrSges(7,2);
t733 = m(7) * t597 + t663 * mrSges(7,3) + t686 * t644;
t631 = t651 * mrSges(7,1) - t652 * mrSges(7,3);
t737 = -t651 * mrSges(6,1) - t652 * mrSges(6,2) - t631;
t590 = m(6) * t600 - t663 * mrSges(6,2) + t747 * t616 - t686 * t643 + t737 * t651 + t733;
t724 = t743 * t604 - t710 * t606;
t599 = t652 * t749 + t724;
t617 = t710 * t639 + t743 * t640;
t642 = -t686 * mrSges(6,2) - t651 * mrSges(6,3);
t598 = -t663 * pkin(5) - t685 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t630) * t652 - t724;
t641 = -t651 * mrSges(7,2) + t686 * mrSges(7,3);
t727 = -m(7) * t598 + t663 * mrSges(7,1) + t686 * t641;
t592 = m(6) * t599 + t663 * mrSges(6,1) + t747 * t617 + t686 * t642 + t737 * t652 + t727;
t585 = t710 * t590 + t743 * t592;
t656 = -t678 * mrSges(5,1) + t679 * mrSges(5,2);
t667 = -t686 * mrSges(5,2) + t678 * mrSges(5,3);
t583 = m(5) * t607 + t663 * mrSges(5,1) - t640 * mrSges(5,3) - t679 * t656 + t686 * t667 + t585;
t669 = t686 * mrSges(5,1) - t679 * mrSges(5,3);
t728 = t743 * t590 - t710 * t592;
t584 = m(5) * t608 - t663 * mrSges(5,2) + t639 * mrSges(5,3) + t678 * t656 - t686 * t669 + t728;
t729 = -t711 * t583 + t715 * t584;
t578 = m(4) * t634 - t707 * mrSges(4,2) + t664 * mrSges(4,3) - t690 * t674 - t708 * t683 + t729;
t633 = t716 * t657 - t712 * t658;
t682 = -t708 * mrSges(4,2) - t690 * mrSges(4,3);
t619 = -t707 * pkin(3) - t706 * pkin(9) + t691 * t675 - t633;
t609 = -t639 * pkin(4) - t677 * qJ(5) + t679 * t668 + qJDD(5) + t619;
t602 = -0.2e1 * qJD(6) * t652 + (t651 * t686 - t617) * qJ(6) + (t652 * t686 + t616) * pkin(5) + t609;
t595 = m(7) * t602 + t616 * mrSges(7,1) - t617 * mrSges(7,3) + t651 * t641 - t652 * t644;
t722 = m(6) * t609 + t616 * mrSges(6,1) + t617 * mrSges(6,2) + t651 * t642 + t652 * t643 + t595;
t720 = -m(5) * t619 + t639 * mrSges(5,1) - t640 * mrSges(5,2) + t678 * t667 - t679 * t669 - t722;
t594 = m(4) * t633 + t707 * mrSges(4,1) - t665 * mrSges(4,3) - t691 * t674 + t708 * t682 + t720;
t573 = t712 * t578 + t716 * t594;
t680 = -t717 * g(3) - t742;
t697 = (-mrSges(3,1) * t717 + mrSges(3,2) * t713) * qJD(1);
t735 = qJD(1) * t717;
t701 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t735;
t571 = m(3) * t680 + qJDD(2) * mrSges(3,1) - t698 * mrSges(3,3) + qJD(2) * t701 - t697 * t736 + t573;
t700 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t736;
t730 = t716 * t578 - t712 * t594;
t572 = m(3) * t681 - qJDD(2) * mrSges(3,2) + t699 * mrSges(3,3) - qJD(2) * t700 + t697 * t735 + t730;
t731 = -t713 * t571 + t717 * t572;
t563 = m(2) * t704 - t719 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t731;
t692 = -t719 * pkin(7) + t725;
t579 = t715 * t583 + t711 * t584;
t723 = m(4) * t666 - t664 * mrSges(4,1) + t665 * mrSges(4,2) + t690 * t682 + t691 * t683 + t579;
t721 = -m(3) * t692 + t699 * mrSges(3,1) - t698 * mrSges(3,2) - t700 * t736 + t701 * t735 - t723;
t575 = m(2) * t703 + qJDD(1) * mrSges(2,1) - t719 * mrSges(2,2) + t721;
t741 = t714 * t563 + t718 * t575;
t564 = t717 * t571 + t713 * t572;
t740 = t651 * t751 - t652 * t746 - t686 * t744;
t739 = t651 * t744 - t652 * t745 + t686 * t750;
t738 = -t746 * t651 + t652 * t752 + t745 * t686;
t732 = t718 * t563 - t714 * t575;
t689 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t713 + Ifges(3,4) * t717) * qJD(1);
t688 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t713 + Ifges(3,2) * t717) * qJD(1);
t687 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t713 + Ifges(3,6) * t717) * qJD(1);
t672 = Ifges(4,1) * t691 - Ifges(4,4) * t690 + Ifges(4,5) * t708;
t671 = Ifges(4,4) * t691 - Ifges(4,2) * t690 + Ifges(4,6) * t708;
t670 = Ifges(4,5) * t691 - Ifges(4,6) * t690 + Ifges(4,3) * t708;
t647 = Ifges(5,1) * t679 + Ifges(5,4) * t678 + Ifges(5,5) * t686;
t646 = Ifges(5,4) * t679 + Ifges(5,2) * t678 + Ifges(5,6) * t686;
t645 = Ifges(5,5) * t679 + Ifges(5,6) * t678 + Ifges(5,3) * t686;
t587 = mrSges(6,2) * t609 + mrSges(7,2) * t598 - mrSges(6,3) * t599 - mrSges(7,3) * t602 - qJ(6) * t595 - t746 * t616 + t617 * t752 + t739 * t651 + t745 * t663 + t740 * t686;
t586 = -mrSges(6,1) * t609 - mrSges(7,1) * t602 + mrSges(7,2) * t597 + mrSges(6,3) * t600 - pkin(5) * t595 - t616 * t751 + t746 * t617 + t739 * t652 + t744 * t663 + t738 * t686;
t567 = mrSges(5,2) * t619 - mrSges(5,3) * t607 + Ifges(5,1) * t640 + Ifges(5,4) * t639 + Ifges(5,5) * t663 - qJ(5) * t585 - t710 * t586 + t743 * t587 + t678 * t645 - t686 * t646;
t566 = -mrSges(5,1) * t619 + mrSges(5,3) * t608 + Ifges(5,4) * t640 + Ifges(5,2) * t639 + Ifges(5,6) * t663 - pkin(4) * t722 + qJ(5) * t728 + t743 * t586 + t710 * t587 - t679 * t645 + t686 * t647;
t565 = (qJ(6) * mrSges(7,2) + t744) * t616 + (pkin(5) * mrSges(7,2) - t745) * t617 + (-Ifges(5,3) + t750) * t663 + Ifges(4,6) * t707 + t708 * t672 - t691 * t670 + t678 * t647 - t679 * t646 + Ifges(4,2) * t664 + Ifges(4,4) * t665 - mrSges(4,1) * t666 - Ifges(5,6) * t639 - Ifges(5,5) * t640 + mrSges(4,3) * t634 - qJ(6) * t733 - pkin(4) * t585 + mrSges(7,1) * t598 - mrSges(6,1) * t599 - pkin(3) * t579 - mrSges(7,3) * t597 - mrSges(5,1) * t607 + mrSges(5,2) * t608 + (qJ(6) * t631 - t738) * t651 + (pkin(5) * t631 + t740) * t652 + mrSges(6,2) * t600 - pkin(5) * t727;
t560 = mrSges(4,2) * t666 - mrSges(4,3) * t633 + Ifges(4,1) * t665 + Ifges(4,4) * t664 + Ifges(4,5) * t707 - pkin(9) * t579 - t711 * t566 + t715 * t567 - t690 * t670 - t708 * t671;
t559 = mrSges(3,2) * t692 - mrSges(3,3) * t680 + Ifges(3,1) * t698 + Ifges(3,4) * t699 + Ifges(3,5) * qJDD(2) - pkin(8) * t573 - qJD(2) * t688 + t716 * t560 - t712 * t565 + t687 * t735;
t558 = -pkin(1) * t564 + mrSges(2,3) * t704 - pkin(2) * t573 - Ifges(3,5) * t698 - Ifges(3,6) * t699 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t680 + mrSges(3,2) * t681 - t711 * t567 - t715 * t566 - pkin(3) * t720 - pkin(9) * t729 - Ifges(4,5) * t665 - Ifges(4,6) * t664 - Ifges(4,3) * t707 - mrSges(4,1) * t633 + mrSges(4,2) * t634 + mrSges(2,1) * g(3) + t719 * Ifges(2,5) - t691 * t671 - t690 * t672 + Ifges(2,6) * qJDD(1) + (-t713 * t688 + t717 * t689) * qJD(1);
t557 = -mrSges(3,1) * t692 + mrSges(3,3) * t681 + Ifges(3,4) * t698 + Ifges(3,2) * t699 + Ifges(3,6) * qJDD(2) - pkin(2) * t723 + pkin(8) * t730 + qJD(2) * t689 + t712 * t560 + t716 * t565 - t687 * t736;
t556 = -mrSges(2,2) * g(3) - mrSges(2,3) * t703 + Ifges(2,5) * qJDD(1) - t719 * Ifges(2,6) - pkin(7) * t564 - t713 * t557 + t717 * t559;
t1 = [-m(1) * g(1) + t732; -m(1) * g(2) + t741; (-m(1) - m(2)) * g(3) + t564; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t741 + t718 * t556 - t714 * t558; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t732 + t714 * t556 + t718 * t558; -mrSges(1,1) * g(2) + mrSges(2,1) * t703 + mrSges(1,2) * g(1) - mrSges(2,2) * t704 + Ifges(2,3) * qJDD(1) + pkin(1) * t721 + pkin(7) * t731 + t717 * t557 + t713 * t559;];
tauB  = t1;
