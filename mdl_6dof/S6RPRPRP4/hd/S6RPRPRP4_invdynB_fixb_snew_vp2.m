% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:42:50
% EndTime: 2019-05-05 17:42:55
% DurationCPUTime: 3.24s
% Computational Cost: add. (31262->304), mult. (60773->350), div. (0->0), fcn. (31820->8), ass. (0->125)
t718 = -2 * qJD(4);
t717 = Ifges(4,1) + Ifges(5,2);
t716 = Ifges(6,1) + Ifges(7,1);
t707 = Ifges(4,4) + Ifges(5,6);
t706 = Ifges(6,4) - Ifges(7,5);
t705 = Ifges(7,4) + Ifges(6,5);
t704 = Ifges(4,5) - Ifges(5,4);
t715 = Ifges(4,2) + Ifges(5,3);
t714 = Ifges(6,2) + Ifges(7,3);
t703 = Ifges(4,6) - Ifges(5,5);
t702 = Ifges(6,6) - Ifges(7,6);
t713 = (Ifges(4,3) + Ifges(5,1));
t712 = Ifges(6,3) + Ifges(7,2);
t665 = sin(qJ(1));
t668 = cos(qJ(1));
t646 = t665 * g(1) - g(2) * t668;
t632 = qJDD(1) * pkin(1) + t646;
t647 = -g(1) * t668 - g(2) * t665;
t670 = qJD(1) ^ 2;
t636 = -pkin(1) * t670 + t647;
t661 = sin(pkin(9));
t662 = cos(pkin(9));
t599 = t661 * t632 + t662 * t636;
t585 = -pkin(2) * t670 + qJDD(1) * pkin(7) + t599;
t660 = -g(3) + qJDD(2);
t664 = sin(qJ(3));
t667 = cos(qJ(3));
t581 = t667 * t585 + t664 * t660;
t633 = (-pkin(3) * t667 - qJ(4) * t664) * qJD(1);
t669 = qJD(3) ^ 2;
t690 = qJD(1) * t667;
t578 = t669 * pkin(3) - qJDD(3) * qJ(4) + (qJD(3) * t718) - t633 * t690 - t581;
t711 = -pkin(3) - pkin(8);
t710 = pkin(8) * t670;
t709 = t670 * pkin(7);
t708 = -mrSges(6,3) - mrSges(7,2);
t701 = t667 * t660;
t582 = t664 * t585;
t580 = -t582 + t701;
t634 = (mrSges(5,2) * t667 - mrSges(5,3) * t664) * qJD(1);
t635 = (-mrSges(4,1) * t667 + mrSges(4,2) * t664) * qJD(1);
t689 = qJD(1) * qJD(3);
t685 = t667 * t689;
t637 = qJDD(1) * t664 + t685;
t642 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t690;
t643 = -mrSges(5,1) * t690 - qJD(3) * mrSges(5,3);
t684 = t664 * t689;
t638 = qJDD(1) * t667 - t684;
t691 = qJD(1) * t664;
t645 = pkin(4) * t691 - qJD(3) * pkin(8);
t659 = t667 ^ 2;
t598 = t662 * t632 - t661 * t636;
t678 = -qJDD(1) * pkin(2) - t598;
t674 = pkin(3) * t684 + t691 * t718 + (-t637 - t685) * qJ(4) + t678;
t572 = -t645 * t691 + (-pkin(4) * t659 - pkin(7)) * t670 + t711 * t638 + t674;
t679 = -t669 * qJ(4) + t633 * t691 + qJDD(4) + t582;
t576 = t637 * pkin(4) + t711 * qJDD(3) + (-pkin(4) * t689 - t664 * t710 - t660) * t667 + t679;
t663 = sin(qJ(5));
t666 = cos(qJ(5));
t570 = t666 * t572 + t663 * t576;
t631 = qJD(3) * t666 - t663 * t690;
t596 = qJD(5) * t631 + qJDD(3) * t663 + t666 * t638;
t650 = qJD(5) + t691;
t607 = mrSges(6,1) * t650 - mrSges(6,3) * t631;
t629 = qJDD(5) + t637;
t630 = qJD(3) * t663 + t666 * t690;
t602 = pkin(5) * t630 - qJ(6) * t631;
t648 = t650 ^ 2;
t565 = -pkin(5) * t648 + qJ(6) * t629 + 0.2e1 * qJD(6) * t650 - t602 * t630 + t570;
t608 = -mrSges(7,1) * t650 + mrSges(7,2) * t631;
t687 = m(7) * t565 + t629 * mrSges(7,3) + t650 * t608;
t603 = mrSges(7,1) * t630 - mrSges(7,3) * t631;
t695 = -mrSges(6,1) * t630 - mrSges(6,2) * t631 - t603;
t560 = m(6) * t570 - t629 * mrSges(6,2) + t708 * t596 - t650 * t607 + t695 * t630 + t687;
t569 = -t572 * t663 + t576 * t666;
t597 = -qJD(5) * t630 + qJDD(3) * t666 - t638 * t663;
t605 = -mrSges(6,2) * t650 - mrSges(6,3) * t630;
t566 = -pkin(5) * t629 - qJ(6) * t648 + t602 * t631 + qJDD(6) - t569;
t606 = -mrSges(7,2) * t630 + mrSges(7,3) * t650;
t680 = -m(7) * t566 + t629 * mrSges(7,1) + t650 * t606;
t562 = m(6) * t569 + t629 * mrSges(6,1) + t708 * t597 + t650 * t605 + t695 * t631 + t680;
t555 = t663 * t560 + t666 * t562;
t579 = -qJDD(3) * pkin(3) + t679 - t701;
t675 = -m(5) * t579 - t637 * mrSges(5,1) - t555;
t551 = m(4) * t580 - t637 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t642 - t643) * qJD(3) + (-t634 - t635) * t691 + t675;
t641 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t691;
t644 = mrSges(5,1) * t691 + qJD(3) * mrSges(5,2);
t575 = t638 * pkin(4) + qJD(3) * t645 - t659 * t710 - t578;
t568 = -0.2e1 * qJD(6) * t631 + (t630 * t650 - t597) * qJ(6) + (t631 * t650 + t596) * pkin(5) + t575;
t563 = m(7) * t568 + t596 * mrSges(7,1) - t597 * mrSges(7,3) + t630 * t606 - t631 * t608;
t673 = m(6) * t575 + t596 * mrSges(6,1) + t597 * mrSges(6,2) + t630 * t605 + t631 * t607 + t563;
t672 = -m(5) * t578 + qJDD(3) * mrSges(5,3) + qJD(3) * t644 + t634 * t690 + t673;
t558 = t672 - qJD(3) * t641 - qJDD(3) * mrSges(4,2) + m(4) * t581 + t635 * t690 + (mrSges(4,3) + mrSges(5,1)) * t638;
t681 = -t551 * t664 + t667 * t558;
t546 = m(3) * t599 - mrSges(3,1) * t670 - qJDD(1) * mrSges(3,2) + t681;
t584 = t678 - t709;
t577 = -t638 * pkin(3) + t674 - t709;
t699 = t666 * t560 - t663 * t562;
t677 = -m(5) * t577 - t638 * mrSges(5,2) + t644 * t691 - t699;
t671 = -m(4) * t584 + t642 * t690 + t638 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t637 + (-t641 * t664 - t643 * t667) * qJD(1) + t677;
t549 = m(3) * t598 + qJDD(1) * mrSges(3,1) - t670 * mrSges(3,2) + t671;
t543 = t661 * t546 + t662 * t549;
t541 = m(2) * t646 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t670 + t543;
t682 = t662 * t546 - t549 * t661;
t542 = m(2) * t647 - mrSges(2,1) * t670 - qJDD(1) * mrSges(2,2) + t682;
t700 = t668 * t541 + t665 * t542;
t547 = t667 * t551 + t664 * t558;
t698 = t714 * t630 - t706 * t631 - t702 * t650;
t697 = t702 * t630 - t705 * t631 - t712 * t650;
t696 = -t706 * t630 + t716 * t631 + t705 * t650;
t694 = (t713 * qJD(3)) + (t704 * t664 + t703 * t667) * qJD(1);
t693 = -t703 * qJD(3) + (-t707 * t664 - t715 * t667) * qJD(1);
t692 = t704 * qJD(3) + (t717 * t664 + t707 * t667) * qJD(1);
t686 = m(3) * t660 + t547;
t683 = -t541 * t665 + t668 * t542;
t554 = mrSges(6,2) * t575 + mrSges(7,2) * t566 - mrSges(6,3) * t569 - mrSges(7,3) * t568 - qJ(6) * t563 - t706 * t596 + t716 * t597 + t705 * t629 + t697 * t630 + t698 * t650;
t553 = -mrSges(6,1) * t575 - mrSges(7,1) * t568 + mrSges(7,2) * t565 + mrSges(6,3) * t570 - pkin(5) * t563 - t714 * t596 + t706 * t597 + t702 * t629 + t697 * t631 + t696 * t650;
t552 = -t637 * mrSges(5,3) + t643 * t690 - t677;
t537 = pkin(5) * t680 + qJ(6) * t687 + mrSges(4,2) * t584 - mrSges(5,3) * t577 + mrSges(5,1) * t579 - mrSges(4,3) * t580 - mrSges(7,1) * t566 + mrSges(6,1) * t569 - mrSges(6,2) * t570 + mrSges(7,3) * t565 + pkin(4) * t555 - qJ(4) * t552 + t707 * t638 + t717 * t637 + (-pkin(5) * t603 - t698) * t631 + (-qJ(6) * t603 + t696) * t630 + t712 * t629 + (-mrSges(7,2) * pkin(5) + t705) * t597 + (-mrSges(7,2) * qJ(6) - t702) * t596 + t704 * qJDD(3) + t693 * qJD(3) + t694 * t690;
t536 = -mrSges(4,1) * t584 - mrSges(5,1) * t578 + mrSges(5,2) * t577 + mrSges(4,3) * t581 - pkin(3) * t552 + pkin(4) * t673 - pkin(8) * t699 + t692 * qJD(3) + t703 * qJDD(3) - t666 * t553 - t663 * t554 + t707 * t637 + t715 * t638 - t694 * t691;
t535 = -pkin(2) * t547 - mrSges(3,1) * t660 + mrSges(3,3) * t599 - pkin(3) * (-qJD(3) * t643 + t675) - qJ(4) * t672 - t666 * t554 + t663 * t553 + pkin(8) * t555 - mrSges(4,1) * t580 + mrSges(4,2) * t581 - mrSges(5,2) * t579 + mrSges(5,3) * t578 + t670 * Ifges(3,5) + Ifges(3,6) * qJDD(1) + (-mrSges(5,1) * qJ(4) - t703) * t638 - t704 * t637 + (mrSges(5,2) * pkin(3) - t713) * qJDD(3) + (t692 * t667 + (pkin(3) * t634 + t693) * t664) * qJD(1);
t534 = mrSges(3,2) * t660 - mrSges(3,3) * t598 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t670 - pkin(7) * t547 - t536 * t664 + t537 * t667;
t533 = -mrSges(2,2) * g(3) - mrSges(2,3) * t646 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t670 - qJ(2) * t543 + t534 * t662 - t535 * t661;
t532 = mrSges(2,1) * g(3) + mrSges(2,3) * t647 + t670 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t686 + qJ(2) * t682 + t661 * t534 + t662 * t535;
t1 = [-m(1) * g(1) + t683; -m(1) * g(2) + t700; (-m(1) - m(2)) * g(3) + t686; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t700 - t665 * t532 + t668 * t533; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t683 + t668 * t532 + t665 * t533; pkin(1) * t543 + mrSges(2,1) * t646 - mrSges(2,2) * t647 + t667 * t536 + pkin(2) * t671 + pkin(7) * t681 + t664 * t537 + mrSges(3,1) * t598 - mrSges(3,2) * t599 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
