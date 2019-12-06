% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:48:04
% EndTime: 2019-12-05 16:48:08
% DurationCPUTime: 2.70s
% Computational Cost: add. (23993->238), mult. (48343->292), div. (0->0), fcn. (30154->8), ass. (0->99)
t668 = Ifges(5,4) + Ifges(6,4);
t677 = Ifges(5,2) + Ifges(6,2);
t673 = Ifges(5,6) + Ifges(6,6);
t640 = sin(qJ(4));
t641 = sin(qJ(3));
t643 = cos(qJ(4));
t644 = cos(qJ(3));
t612 = (-t640 * t641 + t643 * t644) * qJD(2);
t661 = qJD(2) * qJD(3);
t657 = t644 * t661;
t622 = t641 * qJDD(2) + t657;
t623 = t644 * qJDD(2) - t641 * t661;
t583 = t612 * qJD(4) + t643 * t622 + t640 * t623;
t613 = (t640 * t644 + t641 * t643) * qJD(2);
t595 = -t612 * mrSges(6,1) + t613 * mrSges(6,2);
t639 = sin(pkin(8));
t667 = cos(pkin(8));
t625 = -t667 * g(1) - t639 * g(2);
t638 = -g(3) + qJDD(1);
t642 = sin(qJ(2));
t645 = cos(qJ(2));
t607 = t645 * t625 + t642 * t638;
t646 = qJD(2) ^ 2;
t600 = -t646 * pkin(2) + qJDD(2) * pkin(6) + t607;
t624 = t639 * g(1) - t667 * g(2);
t584 = -t641 * t600 - t644 * t624;
t568 = (-t622 + t657) * pkin(7) + (t641 * t644 * t646 + qJDD(3)) * pkin(3) + t584;
t585 = t644 * t600 - t641 * t624;
t663 = qJD(2) * t641;
t628 = qJD(3) * pkin(3) - pkin(7) * t663;
t637 = t644 ^ 2;
t569 = -t637 * t646 * pkin(3) + t623 * pkin(7) - qJD(3) * t628 + t585;
t563 = t643 * t568 - t640 * t569;
t635 = qJDD(3) + qJDD(4);
t636 = qJD(3) + qJD(4);
t557 = -0.2e1 * qJD(5) * t613 + (t612 * t636 - t583) * qJ(5) + (t612 * t613 + t635) * pkin(4) + t563;
t601 = -t636 * mrSges(6,2) + t612 * mrSges(6,3);
t660 = m(6) * t557 + t635 * mrSges(6,1) + t636 * t601;
t553 = -t583 * mrSges(6,3) - t613 * t595 + t660;
t564 = t640 * t568 + t643 * t569;
t582 = -t613 * qJD(4) - t640 * t622 + t643 * t623;
t603 = t636 * pkin(4) - t613 * qJ(5);
t608 = t612 ^ 2;
t559 = -t608 * pkin(4) + t582 * qJ(5) + 0.2e1 * qJD(5) * t612 - t636 * t603 + t564;
t674 = Ifges(5,5) + Ifges(6,5);
t675 = Ifges(5,1) + Ifges(6,1);
t664 = -t668 * t612 - t675 * t613 - t674 * t636;
t671 = t677 * t612 + t668 * t613 + t673 * t636;
t672 = Ifges(5,3) + Ifges(6,3);
t676 = mrSges(5,1) * t563 + mrSges(6,1) * t557 - mrSges(5,2) * t564 - mrSges(6,2) * t559 + pkin(4) * t553 + t673 * t582 + t674 * t583 + t664 * t612 + t671 * t613 + t672 * t635;
t596 = -t612 * mrSges(5,1) + t613 * mrSges(5,2);
t602 = -t636 * mrSges(5,2) + t612 * mrSges(5,3);
t545 = m(5) * t563 + t635 * mrSges(5,1) + t636 * t602 + (-t595 - t596) * t613 + (-mrSges(5,3) - mrSges(6,3)) * t583 + t660;
t604 = t636 * mrSges(6,1) - t613 * mrSges(6,3);
t605 = t636 * mrSges(5,1) - t613 * mrSges(5,3);
t659 = m(6) * t559 + t582 * mrSges(6,3) + t612 * t595;
t550 = m(5) * t564 + t582 * mrSges(5,3) + t612 * t596 + (-t604 - t605) * t636 + (-mrSges(5,2) - mrSges(6,2)) * t635 + t659;
t543 = t643 * t545 + t640 * t550;
t610 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t641 + Ifges(4,2) * t644) * qJD(2);
t611 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t641 + Ifges(4,4) * t644) * qJD(2);
t670 = mrSges(4,1) * t584 - mrSges(4,2) * t585 + Ifges(4,5) * t622 + Ifges(4,6) * t623 + Ifges(4,3) * qJDD(3) + pkin(3) * t543 + (t641 * t610 - t644 * t611) * qJD(2) + t676;
t621 = (-mrSges(4,1) * t644 + mrSges(4,2) * t641) * qJD(2);
t662 = qJD(2) * t644;
t627 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t662;
t541 = m(4) * t584 + qJDD(3) * mrSges(4,1) - t622 * mrSges(4,3) + qJD(3) * t627 - t621 * t663 + t543;
t626 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t663;
t653 = -t640 * t545 + t643 * t550;
t542 = m(4) * t585 - qJDD(3) * mrSges(4,2) + t623 * mrSges(4,3) - qJD(3) * t626 + t621 * t662 + t653;
t537 = -t641 * t541 + t644 * t542;
t533 = m(3) * t607 - t646 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t537;
t606 = -t642 * t625 + t645 * t638;
t651 = -qJDD(2) * pkin(2) - t606;
t599 = -t646 * pkin(6) + t651;
t570 = -t623 * pkin(3) + t628 * t663 + (-pkin(7) * t637 - pkin(6)) * t646 + t651;
t561 = -t582 * pkin(4) - t608 * qJ(5) + t613 * t603 + qJDD(5) + t570;
t554 = m(6) * t561 - t582 * mrSges(6,1) + t583 * mrSges(6,2) - t612 * t601 + t613 * t604;
t650 = m(5) * t570 - t582 * mrSges(5,1) + t583 * mrSges(5,2) - t612 * t602 + t613 * t605 + t554;
t551 = -m(4) * t599 + t623 * mrSges(4,1) - t622 * mrSges(4,2) - t626 * t663 + t627 * t662 - t650;
t549 = m(3) * t606 + qJDD(2) * mrSges(3,1) - t646 * mrSges(3,2) + t551;
t654 = t645 * t533 - t642 * t549;
t529 = m(2) * t625 + t654;
t536 = t644 * t541 + t641 * t542;
t535 = (m(2) + m(3)) * t624 - t536;
t666 = t639 * t529 + t667 * t535;
t530 = t642 * t533 + t645 * t549;
t665 = -t673 * t612 - t674 * t613 - t672 * t636;
t658 = m(2) * t638 + t530;
t655 = t667 * t529 - t639 * t535;
t538 = -mrSges(5,1) * t570 + mrSges(5,3) * t564 - mrSges(6,1) * t561 + mrSges(6,3) * t559 - pkin(4) * t554 + qJ(5) * t659 + (-qJ(5) * t604 - t664) * t636 + (-qJ(5) * mrSges(6,2) + t673) * t635 + t665 * t613 + t668 * t583 + t677 * t582;
t540 = mrSges(5,2) * t570 + mrSges(6,2) * t561 - mrSges(5,3) * t563 - mrSges(6,3) * t557 - qJ(5) * t553 + t668 * t582 + t675 * t583 - t665 * t612 + t674 * t635 - t671 * t636;
t609 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t641 + Ifges(4,6) * t644) * qJD(2);
t525 = -mrSges(4,1) * t599 + mrSges(4,3) * t585 + Ifges(4,4) * t622 + Ifges(4,2) * t623 + Ifges(4,6) * qJDD(3) - pkin(3) * t650 + pkin(7) * t653 + qJD(3) * t611 + t643 * t538 + t640 * t540 - t609 * t663;
t526 = mrSges(4,2) * t599 - mrSges(4,3) * t584 + Ifges(4,1) * t622 + Ifges(4,4) * t623 + Ifges(4,5) * qJDD(3) - pkin(7) * t543 - qJD(3) * t610 - t640 * t538 + t643 * t540 + t609 * t662;
t648 = mrSges(3,1) * t606 - mrSges(3,2) * t607 + Ifges(3,3) * qJDD(2) + pkin(2) * t551 + pkin(6) * t537 + t644 * t525 + t641 * t526;
t524 = mrSges(3,1) * t624 + mrSges(3,3) * t607 + t646 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t536 - t670;
t523 = -mrSges(3,2) * t624 - mrSges(3,3) * t606 + Ifges(3,5) * qJDD(2) - t646 * Ifges(3,6) - pkin(6) * t536 - t641 * t525 + t644 * t526;
t522 = -mrSges(2,1) * t638 + mrSges(2,3) * t625 - pkin(1) * t530 - t648;
t521 = mrSges(2,2) * t638 - mrSges(2,3) * t624 - pkin(5) * t530 + t645 * t523 - t642 * t524;
t1 = [-m(1) * g(1) + t655; -m(1) * g(2) + t666; -m(1) * g(3) + t658; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t666 + t667 * t521 - t639 * t522; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t655 + t639 * t521 + t667 * t522; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t624 - mrSges(2,2) * t625 + t642 * t523 + t645 * t524 + pkin(1) * (m(3) * t624 - t536) + pkin(5) * t654; t658; t648; t670; t676; t554;];
tauJB = t1;
