% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:15:35
% EndTime: 2019-05-05 22:15:45
% DurationCPUTime: 6.26s
% Computational Cost: add. (74536->322), mult. (141925->386), div. (0->0), fcn. (88403->10), ass. (0->129)
t698 = Ifges(5,1) + Ifges(6,1);
t691 = Ifges(5,4) - Ifges(6,5);
t690 = Ifges(5,5) + Ifges(6,4);
t697 = Ifges(5,2) + Ifges(6,3);
t689 = Ifges(5,6) - Ifges(6,6);
t696 = -Ifges(5,3) - Ifges(6,2);
t657 = sin(qJ(4));
t658 = sin(qJ(3));
t682 = qJD(1) * t658;
t693 = cos(qJ(4));
t630 = -qJD(3) * t693 + t657 * t682;
t661 = cos(qJ(3));
t680 = qJD(1) * qJD(3);
t677 = t661 * t680;
t636 = qJDD(1) * t658 + t677;
t600 = -t630 * qJD(4) + t657 * qJDD(3) + t636 * t693;
t659 = sin(qJ(1));
t662 = cos(qJ(1));
t641 = t659 * g(1) - g(2) * t662;
t632 = qJDD(1) * pkin(1) + t641;
t642 = -g(1) * t662 - g(2) * t659;
t664 = qJD(1) ^ 2;
t634 = -pkin(1) * t664 + t642;
t654 = sin(pkin(10));
t655 = cos(pkin(10));
t604 = t654 * t632 + t655 * t634;
t588 = -pkin(2) * t664 + qJDD(1) * pkin(7) + t604;
t653 = -g(3) + qJDD(2);
t579 = -t658 * t588 + t661 * t653;
t635 = (-pkin(3) * t661 - pkin(8) * t658) * qJD(1);
t663 = qJD(3) ^ 2;
t668 = qJDD(3) * pkin(3) + pkin(8) * t663 - t635 * t682 + t579;
t681 = qJD(1) * t661;
t645 = qJD(4) - t681;
t688 = t630 * t645;
t695 = (-t600 + t688) * qJ(5) - t668;
t694 = 2 * qJD(5);
t692 = -mrSges(5,3) - mrSges(6,2);
t580 = t661 * t588 + t658 * t653;
t633 = (-mrSges(4,1) * t661 + mrSges(4,2) * t658) * qJD(1);
t678 = t658 * t680;
t637 = t661 * qJDD(1) - t678;
t638 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t682;
t603 = t632 * t655 - t654 * t634;
t587 = -qJDD(1) * pkin(2) - pkin(7) * t664 - t603;
t571 = (-t636 - t677) * pkin(8) + (-t637 + t678) * pkin(3) + t587;
t577 = -pkin(3) * t663 + qJDD(3) * pkin(8) + t635 * t681 + t580;
t562 = t657 * t571 + t693 * t577;
t631 = t657 * qJD(3) + t682 * t693;
t599 = qJD(4) * t631 - qJDD(3) * t693 + t636 * t657;
t611 = mrSges(5,1) * t645 - mrSges(5,3) * t631;
t629 = qJDD(4) - t637;
t607 = pkin(4) * t630 - qJ(5) * t631;
t644 = t645 ^ 2;
t558 = -pkin(4) * t644 + t629 * qJ(5) - t630 * t607 + t645 * t694 + t562;
t612 = -mrSges(6,1) * t645 + mrSges(6,2) * t631;
t561 = t571 * t693 - t657 * t577;
t559 = -t629 * pkin(4) - t644 * qJ(5) + t631 * t607 + qJDD(5) - t561;
t553 = (-t600 - t688) * pkin(9) + (t630 * t631 - t629) * pkin(5) + t559;
t614 = -pkin(5) * t645 - pkin(9) * t631;
t628 = t630 ^ 2;
t554 = -pkin(5) * t628 + pkin(9) * t599 + t614 * t645 + t558;
t656 = sin(qJ(6));
t660 = cos(qJ(6));
t551 = t553 * t660 - t554 * t656;
t601 = t630 * t660 - t631 * t656;
t566 = qJD(6) * t601 + t599 * t656 + t600 * t660;
t602 = t630 * t656 + t631 * t660;
t578 = -mrSges(7,1) * t601 + mrSges(7,2) * t602;
t643 = qJD(6) - t645;
t581 = -mrSges(7,2) * t643 + mrSges(7,3) * t601;
t625 = qJDD(6) - t629;
t549 = m(7) * t551 + mrSges(7,1) * t625 - mrSges(7,3) * t566 - t578 * t602 + t581 * t643;
t552 = t553 * t656 + t554 * t660;
t565 = -qJD(6) * t602 + t599 * t660 - t600 * t656;
t582 = mrSges(7,1) * t643 - mrSges(7,3) * t602;
t550 = m(7) * t552 - mrSges(7,2) * t625 + mrSges(7,3) * t565 + t578 * t601 - t582 * t643;
t672 = -t549 * t656 + t660 * t550;
t669 = m(6) * t558 + t629 * mrSges(6,3) + t645 * t612 + t672;
t608 = mrSges(6,1) * t630 - mrSges(6,3) * t631;
t683 = -mrSges(5,1) * t630 - mrSges(5,2) * t631 - t608;
t540 = m(5) * t562 - mrSges(5,2) * t629 + t599 * t692 - t611 * t645 + t630 * t683 + t669;
t610 = -mrSges(5,2) * t645 - mrSges(5,3) * t630;
t542 = t549 * t660 + t550 * t656;
t613 = -mrSges(6,2) * t630 + mrSges(6,3) * t645;
t667 = -m(6) * t559 + t629 * mrSges(6,1) + t645 * t613 - t542;
t541 = m(5) * t561 + mrSges(5,1) * t629 + t600 * t692 + t610 * t645 + t631 * t683 + t667;
t673 = t693 * t540 - t541 * t657;
t537 = m(4) * t580 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t637 - qJD(3) * t638 + t633 * t681 + t673;
t639 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t681;
t560 = -0.2e1 * qJD(5) * t631 + (t631 * t645 + t599) * pkin(4) + t695;
t556 = -pkin(9) * t628 + (-pkin(4) - pkin(5)) * t599 + (-pkin(4) * t645 + t614 + t694) * t631 - t695;
t670 = -m(7) * t556 + t565 * mrSges(7,1) - t566 * mrSges(7,2) + t601 * t581 - t602 * t582;
t547 = m(6) * t560 + t599 * mrSges(6,1) - t600 * mrSges(6,3) - t631 * t612 + t630 * t613 + t670;
t665 = m(5) * t668 - t599 * mrSges(5,1) - t600 * mrSges(5,2) - t630 * t610 - t631 * t611 - t547;
t546 = m(4) * t579 + qJDD(3) * mrSges(4,1) - t636 * mrSges(4,3) + qJD(3) * t639 - t633 * t682 + t665;
t674 = t661 * t537 - t546 * t658;
t531 = m(3) * t604 - mrSges(3,1) * t664 - qJDD(1) * mrSges(3,2) + t674;
t538 = t657 * t540 + t541 * t693;
t666 = -m(4) * t587 + t637 * mrSges(4,1) - t636 * mrSges(4,2) - t638 * t682 + t639 * t681 - t538;
t534 = m(3) * t603 + qJDD(1) * mrSges(3,1) - t664 * mrSges(3,2) + t666;
t526 = t654 * t531 + t655 * t534;
t524 = m(2) * t641 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t664 + t526;
t675 = t655 * t531 - t534 * t654;
t525 = m(2) * t642 - mrSges(2,1) * t664 - qJDD(1) * mrSges(2,2) + t675;
t687 = t662 * t524 + t659 * t525;
t532 = t658 * t537 + t661 * t546;
t686 = t697 * t630 - t691 * t631 - t689 * t645;
t685 = t689 * t630 - t690 * t631 + t696 * t645;
t684 = -t691 * t630 + t698 * t631 + t690 * t645;
t679 = m(3) * t653 + t532;
t676 = -t524 * t659 + t662 * t525;
t621 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t658 + Ifges(4,4) * t661) * qJD(1);
t620 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t658 + Ifges(4,2) * t661) * qJD(1);
t619 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t658 + Ifges(4,6) * t661) * qJD(1);
t570 = Ifges(7,1) * t602 + Ifges(7,4) * t601 + Ifges(7,5) * t643;
t569 = Ifges(7,4) * t602 + Ifges(7,2) * t601 + Ifges(7,6) * t643;
t568 = Ifges(7,5) * t602 + Ifges(7,6) * t601 + Ifges(7,3) * t643;
t544 = mrSges(7,2) * t556 - mrSges(7,3) * t551 + Ifges(7,1) * t566 + Ifges(7,4) * t565 + Ifges(7,5) * t625 + t568 * t601 - t569 * t643;
t543 = -mrSges(7,1) * t556 + mrSges(7,3) * t552 + Ifges(7,4) * t566 + Ifges(7,2) * t565 + Ifges(7,6) * t625 - t568 * t602 + t570 * t643;
t528 = -mrSges(5,2) * t668 + mrSges(6,2) * t559 - mrSges(5,3) * t561 - mrSges(6,3) * t560 - pkin(9) * t542 - qJ(5) * t547 - t543 * t656 + t544 * t660 - t691 * t599 + t698 * t600 + t690 * t629 + t685 * t630 + t686 * t645;
t527 = mrSges(5,1) * t668 - mrSges(6,1) * t560 + mrSges(6,2) * t558 + mrSges(5,3) * t562 - pkin(4) * t547 - pkin(5) * t670 - pkin(9) * t672 - t660 * t543 - t656 * t544 - t697 * t599 + t691 * t600 + t689 * t629 + t685 * t631 + t684 * t645;
t520 = (mrSges(6,2) * qJ(5) + t689) * t599 + (mrSges(6,2) * pkin(4) - t690) * t600 - pkin(4) * t667 - qJ(5) * t669 - t619 * t682 + Ifges(4,6) * qJDD(3) - mrSges(7,2) * t552 - t601 * t570 + t602 * t569 - mrSges(4,1) * t587 + qJD(3) * t621 + Ifges(7,3) * t625 + (qJ(5) * t608 - t684) * t630 + (pkin(4) * t608 + t686) * t631 + mrSges(4,3) * t580 + t696 * t629 - pkin(3) * t538 + mrSges(7,1) * t551 + pkin(5) * t542 + Ifges(4,4) * t636 + Ifges(4,2) * t637 - mrSges(6,3) * t558 + mrSges(6,1) * t559 - mrSges(5,1) * t561 + mrSges(5,2) * t562 + Ifges(7,6) * t565 + Ifges(7,5) * t566;
t519 = mrSges(4,2) * t587 - mrSges(4,3) * t579 + Ifges(4,1) * t636 + Ifges(4,4) * t637 + Ifges(4,5) * qJDD(3) - pkin(8) * t538 - qJD(3) * t620 - t657 * t527 + t528 * t693 + t619 * t681;
t518 = Ifges(3,6) * qJDD(1) + t664 * Ifges(3,5) - mrSges(3,1) * t653 + mrSges(3,3) * t604 - Ifges(4,5) * t636 - Ifges(4,6) * t637 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t579 + mrSges(4,2) * t580 - t657 * t528 - t693 * t527 - pkin(3) * t665 - pkin(8) * t673 - pkin(2) * t532 + (-t620 * t658 + t621 * t661) * qJD(1);
t517 = mrSges(3,2) * t653 - mrSges(3,3) * t603 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t664 - pkin(7) * t532 + t519 * t661 - t520 * t658;
t516 = -mrSges(2,2) * g(3) - mrSges(2,3) * t641 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t664 - qJ(2) * t526 + t517 * t655 - t518 * t654;
t515 = mrSges(2,1) * g(3) + mrSges(2,3) * t642 + t664 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t679 + qJ(2) * t675 + t654 * t517 + t655 * t518;
t1 = [-m(1) * g(1) + t676; -m(1) * g(2) + t687; (-m(1) - m(2)) * g(3) + t679; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t687 - t659 * t515 + t662 * t516; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t676 + t662 * t515 + t659 * t516; pkin(1) * t526 + mrSges(2,1) * t641 - mrSges(2,2) * t642 + pkin(7) * t674 + t658 * t519 + t661 * t520 + pkin(2) * t666 + mrSges(3,1) * t603 - mrSges(3,2) * t604 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
