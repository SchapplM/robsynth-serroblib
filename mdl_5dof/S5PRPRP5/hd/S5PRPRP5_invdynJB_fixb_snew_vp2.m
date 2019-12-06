% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:51
% EndTime: 2019-12-05 15:37:54
% DurationCPUTime: 2.04s
% Computational Cost: add. (17813->215), mult. (39326->263), div. (0->0), fcn. (25229->8), ass. (0->103)
t673 = Ifges(5,1) + Ifges(6,1);
t666 = Ifges(5,4) - Ifges(6,5);
t665 = Ifges(5,5) + Ifges(6,4);
t672 = Ifges(5,2) + Ifges(6,3);
t664 = Ifges(5,6) - Ifges(6,6);
t671 = Ifges(5,3) + Ifges(6,2);
t630 = qJD(2) ^ 2;
t626 = sin(qJ(4));
t625 = cos(pkin(8));
t669 = cos(qJ(4));
t645 = t625 * t669;
t623 = sin(pkin(8));
t651 = t623 * qJD(2);
t598 = -qJD(2) * t645 + t626 * t651;
t635 = t669 * t623 + t625 * t626;
t599 = t635 * qJD(2);
t581 = t598 * mrSges(6,1) - t599 * mrSges(6,3);
t653 = t598 * qJD(4);
t586 = t635 * qJDD(2) - t653;
t624 = sin(pkin(7));
t661 = cos(pkin(7));
t607 = -t661 * g(1) - t624 * g(2);
t622 = -g(3) + qJDD(1);
t627 = sin(qJ(2));
t628 = cos(qJ(2));
t597 = t628 * t607 + t627 * t622;
t589 = -t630 * pkin(2) + qJDD(2) * qJ(3) + t597;
t606 = t624 * g(1) - t661 * g(2);
t650 = qJD(2) * qJD(3);
t654 = -t625 * t606 - 0.2e1 * t623 * t650;
t668 = pkin(3) * t625;
t567 = (-pkin(6) * qJDD(2) + t630 * t668 - t589) * t623 + t654;
t570 = -t623 * t606 + (t589 + 0.2e1 * t650) * t625;
t648 = qJDD(2) * t625;
t619 = t625 ^ 2;
t660 = t619 * t630;
t568 = -pkin(3) * t660 + pkin(6) * t648 + t570;
t563 = t669 * t567 - t626 * t568;
t580 = t598 * pkin(4) - t599 * qJ(5);
t629 = qJD(4) ^ 2;
t561 = -qJDD(4) * pkin(4) - t629 * qJ(5) + t599 * t580 + qJDD(5) - t563;
t595 = -t598 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t641 = -m(6) * t561 + qJDD(4) * mrSges(6,1) + qJD(4) * t595;
t558 = t586 * mrSges(6,2) + t599 * t581 - t641;
t564 = t626 * t567 + t669 * t568;
t560 = -t629 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t598 * t580 + t564;
t649 = qJDD(2) * t623;
t652 = t599 * qJD(4);
t585 = -qJDD(2) * t645 + t626 * t649 + t652;
t594 = -qJD(4) * mrSges(6,1) + t599 * mrSges(6,2);
t647 = m(6) * t560 + qJDD(4) * mrSges(6,3) + qJD(4) * t594;
t656 = t665 * qJD(4) - t666 * t598 + t673 * t599;
t657 = -t664 * qJD(4) + t672 * t598 - t666 * t599;
t670 = t671 * qJDD(4) - t664 * t585 + t665 * t586 + t656 * t598 - t657 * t599 + mrSges(5,1) * t563 - mrSges(6,1) * t561 - mrSges(5,2) * t564 + mrSges(6,3) * t560 - pkin(4) * t558 + qJ(5) * (-t585 * mrSges(6,2) - t598 * t581 + t647);
t667 = -mrSges(5,3) - mrSges(6,2);
t662 = mrSges(4,2) * t623;
t593 = qJD(4) * mrSges(5,1) - t599 * mrSges(5,3);
t655 = -t598 * mrSges(5,1) - t599 * mrSges(5,2) - t581;
t554 = m(5) * t564 - qJDD(4) * mrSges(5,2) - qJD(4) * t593 + t667 * t585 + t655 * t598 + t647;
t592 = -qJD(4) * mrSges(5,2) - t598 * mrSges(5,3);
t555 = m(5) * t563 + qJDD(4) * mrSges(5,1) + qJD(4) * t592 + t667 * t586 + t655 * t599 + t641;
t547 = t626 * t554 + t669 * t555;
t569 = -t623 * t589 + t654;
t636 = mrSges(4,3) * qJDD(2) + t630 * (-mrSges(4,1) * t625 + t662);
t545 = m(4) * t569 - t636 * t623 + t547;
t642 = t669 * t554 - t626 * t555;
t546 = m(4) * t570 + t636 * t625 + t642;
t541 = -t623 * t545 + t625 * t546;
t537 = m(3) * t597 - t630 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t541;
t596 = -t627 * t607 + t628 * t622;
t637 = qJDD(3) - t596;
t588 = -qJDD(2) * pkin(2) - t630 * qJ(3) + t637;
t618 = t623 ^ 2;
t571 = (-pkin(2) - t668) * qJDD(2) + (-qJ(3) + (-t618 - t619) * pkin(6)) * t630 + t637;
t562 = -0.2e1 * qJD(5) * t599 + (-t586 + t653) * qJ(5) + (t585 + t652) * pkin(4) + t571;
t556 = m(6) * t562 + t585 * mrSges(6,1) - t586 * mrSges(6,3) - t599 * t594 + t598 * t595;
t632 = m(5) * t571 + t585 * mrSges(5,1) + t586 * mrSges(5,2) + t598 * t592 + t599 * t593 + t556;
t631 = -m(4) * t588 + mrSges(4,1) * t648 - t632 + (t618 * t630 + t660) * mrSges(4,3);
t549 = t631 - t630 * mrSges(3,2) + m(3) * t596 + (mrSges(3,1) - t662) * qJDD(2);
t643 = t628 * t537 - t627 * t549;
t533 = m(2) * t607 + t643;
t540 = t625 * t545 + t623 * t546;
t539 = (m(2) + m(3)) * t606 - t540;
t659 = t624 * t533 + t661 * t539;
t534 = t627 * t537 + t628 * t549;
t658 = -t671 * qJD(4) + t664 * t598 - t665 * t599;
t646 = m(2) * t622 + t534;
t644 = t661 * t533 - t624 * t539;
t640 = Ifges(4,1) * t623 + Ifges(4,4) * t625;
t639 = Ifges(4,4) * t623 + Ifges(4,2) * t625;
t638 = Ifges(4,5) * t623 + Ifges(4,6) * t625;
t542 = -mrSges(5,1) * t571 - mrSges(6,1) * t562 + mrSges(6,2) * t560 + mrSges(5,3) * t564 - pkin(4) * t556 + t656 * qJD(4) + t664 * qJDD(4) - t672 * t585 + t666 * t586 + t658 * t599;
t543 = mrSges(5,2) * t571 + mrSges(6,2) * t561 - mrSges(5,3) * t563 - mrSges(6,3) * t562 - qJ(5) * t556 + t657 * qJD(4) + t665 * qJDD(4) - t666 * t585 + t673 * t586 + t658 * t598;
t605 = t638 * qJD(2);
t529 = -mrSges(4,1) * t588 + mrSges(4,3) * t570 - pkin(3) * t632 + pkin(6) * t642 + t639 * qJDD(2) + t669 * t542 + t626 * t543 - t605 * t651;
t530 = t625 * qJD(2) * t605 + mrSges(4,2) * t588 - mrSges(4,3) * t569 - pkin(6) * t547 + t640 * qJDD(2) - t626 * t542 + t669 * t543;
t550 = mrSges(4,2) * t649 - t631;
t633 = mrSges(3,1) * t596 - mrSges(3,2) * t597 + Ifges(3,3) * qJDD(2) - pkin(2) * t550 + qJ(3) * t541 + t625 * t529 + t623 * t530;
t528 = (Ifges(3,6) - t638) * qJDD(2) + mrSges(3,1) * t606 + mrSges(3,3) * t597 - pkin(3) * t547 - mrSges(4,1) * t569 + mrSges(4,2) * t570 - pkin(2) * t540 + (-t623 * t639 + t625 * t640 + Ifges(3,5)) * t630 - t670;
t527 = -mrSges(3,2) * t606 - mrSges(3,3) * t596 + Ifges(3,5) * qJDD(2) - t630 * Ifges(3,6) - qJ(3) * t540 - t623 * t529 + t625 * t530;
t526 = -mrSges(2,1) * t622 + mrSges(2,3) * t607 - pkin(1) * t534 - t633;
t525 = mrSges(2,2) * t622 - mrSges(2,3) * t606 - pkin(5) * t534 + t628 * t527 - t627 * t528;
t1 = [-m(1) * g(1) + t644; -m(1) * g(2) + t659; -m(1) * g(3) + t646; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t659 + t661 * t525 - t624 * t526; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t644 + t624 * t525 + t661 * t526; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t606 - mrSges(2,2) * t607 + t627 * t527 + t628 * t528 + pkin(1) * (m(3) * t606 - t540) + pkin(5) * t643; t646; t633; t550; t670; t558;];
tauJB = t1;
