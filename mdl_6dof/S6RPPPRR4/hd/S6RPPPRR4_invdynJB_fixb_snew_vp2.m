% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:45:59
% EndTime: 2019-05-05 13:46:01
% DurationCPUTime: 2.05s
% Computational Cost: add. (21925->250), mult. (37159->290), div. (0->0), fcn. (16021->8), ass. (0->104)
t632 = qJD(1) ^ 2;
t627 = sin(qJ(1));
t630 = cos(qJ(1));
t603 = -t630 * g(1) - t627 * g(2);
t639 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t603;
t663 = -pkin(1) - pkin(2);
t580 = t632 * t663 + t639;
t602 = t627 * g(1) - t630 * g(2);
t638 = -t632 * qJ(2) + qJDD(2) - t602;
t583 = qJDD(1) * t663 + t638;
t623 = sin(pkin(9));
t624 = cos(pkin(9));
t566 = t624 * t580 + t623 * t583;
t666 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t566;
t618 = g(3) + qJDD(3);
t565 = -t623 * t580 + t583 * t624;
t563 = qJDD(1) * pkin(3) - qJ(4) * t632 + qJDD(4) - t565;
t561 = qJDD(1) * pkin(7) + t563;
t626 = sin(qJ(5));
t629 = cos(qJ(5));
t556 = t626 * t561 + t629 * t618;
t595 = (-mrSges(6,1) * t626 - mrSges(6,2) * t629) * qJD(1);
t651 = qJD(1) * qJD(5);
t648 = t629 * t651;
t597 = qJDD(1) * t626 + t648;
t652 = qJD(1) * t629;
t601 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t652;
t559 = (-pkin(3) - pkin(7)) * t632 + t666;
t649 = t626 * t651;
t598 = -qJDD(1) * t629 + t649;
t551 = (-t598 - t649) * pkin(8) + (-t597 - t648) * pkin(5) + t559;
t596 = (-pkin(5) * t626 + pkin(8) * t629) * qJD(1);
t631 = qJD(5) ^ 2;
t653 = qJD(1) * t626;
t553 = -pkin(5) * t631 + qJDD(5) * pkin(8) + t596 * t653 + t556;
t625 = sin(qJ(6));
t628 = cos(qJ(6));
t549 = t551 * t628 - t553 * t625;
t593 = qJD(5) * t628 + t625 * t652;
t573 = qJD(6) * t593 + qJDD(5) * t625 + t598 * t628;
t594 = qJD(5) * t625 - t628 * t652;
t574 = -mrSges(7,1) * t593 + mrSges(7,2) * t594;
t604 = qJD(6) - t653;
t578 = -mrSges(7,2) * t604 + mrSges(7,3) * t593;
t592 = qJDD(6) - t597;
t546 = m(7) * t549 + mrSges(7,1) * t592 - t573 * mrSges(7,3) - t574 * t594 + t578 * t604;
t550 = t551 * t625 + t553 * t628;
t572 = -qJD(6) * t594 + qJDD(5) * t628 - t598 * t625;
t579 = mrSges(7,1) * t604 - mrSges(7,3) * t594;
t547 = m(7) * t550 - mrSges(7,2) * t592 + t572 * mrSges(7,3) + t574 * t593 - t579 * t604;
t645 = -t546 * t625 + t628 * t547;
t534 = m(6) * t556 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t597 - qJD(5) * t601 + t595 * t653 + t645;
t656 = t626 * t618;
t555 = t561 * t629 - t656;
t600 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t653;
t552 = -qJDD(5) * pkin(5) - t631 * pkin(8) + t656 + (-qJD(1) * t596 - t561) * t629;
t637 = -m(7) * t552 + t572 * mrSges(7,1) - t573 * mrSges(7,2) + t593 * t578 - t579 * t594;
t542 = m(6) * t555 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t598 + qJD(5) * t600 + t595 * t652 + t637;
t654 = -t629 * t534 + t626 * t542;
t526 = (m(4) + m(5)) * t618 - t654;
t536 = t628 * t546 + t625 * t547;
t665 = -m(6) * t559 + t597 * mrSges(6,1) - t598 * mrSges(6,2) + (t600 * t626 + t601 * t629) * qJD(1) - t536;
t662 = mrSges(2,1) + mrSges(3,1);
t661 = mrSges(4,2) - mrSges(5,3);
t660 = Ifges(3,4) + Ifges(2,5);
t659 = Ifges(5,4) - Ifges(4,5);
t658 = Ifges(5,5) - Ifges(4,6);
t657 = Ifges(2,6) - Ifges(3,6);
t584 = -pkin(1) * t632 + t639;
t529 = t534 * t626 + t542 * t629;
t524 = m(5) * t563 - qJDD(1) * mrSges(5,2) - t632 * mrSges(5,3) + t529;
t523 = m(4) * t565 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t632 - t524;
t562 = t632 * pkin(3) - t666;
t636 = -m(5) * t562 + t632 * mrSges(5,2) - t665;
t531 = m(4) * t566 - t632 * mrSges(4,1) + qJDD(1) * t661 + t636;
t646 = -t623 * t523 + t624 * t531;
t641 = m(3) * t584 + qJDD(1) * mrSges(3,3) + t646;
t514 = m(2) * t603 - qJDD(1) * mrSges(2,2) - t632 * t662 + t641;
t519 = t523 * t624 + t531 * t623;
t585 = -qJDD(1) * pkin(1) + t638;
t518 = m(3) * t585 - qJDD(1) * mrSges(3,1) - t632 * mrSges(3,3) + t519;
t515 = m(2) * t602 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t632 - t518;
t655 = t627 * t514 + t630 * t515;
t647 = t630 * t514 - t515 * t627;
t568 = Ifges(7,4) * t594 + Ifges(7,2) * t593 + Ifges(7,6) * t604;
t569 = Ifges(7,1) * t594 + Ifges(7,4) * t593 + Ifges(7,5) * t604;
t635 = mrSges(7,1) * t549 - mrSges(7,2) * t550 + Ifges(7,5) * t573 + Ifges(7,6) * t572 + Ifges(7,3) * t592 + t568 * t594 - t569 * t593;
t567 = Ifges(7,5) * t594 + Ifges(7,6) * t593 + Ifges(7,3) * t604;
t539 = -mrSges(7,1) * t552 + mrSges(7,3) * t550 + Ifges(7,4) * t573 + Ifges(7,2) * t572 + Ifges(7,6) * t592 - t567 * t594 + t569 * t604;
t540 = mrSges(7,2) * t552 - mrSges(7,3) * t549 + Ifges(7,1) * t573 + Ifges(7,4) * t572 + Ifges(7,5) * t592 + t567 * t593 - t568 * t604;
t587 = (Ifges(6,6) * qJD(5)) + (-Ifges(6,4) * t629 + Ifges(6,2) * t626) * qJD(1);
t588 = (Ifges(6,5) * qJD(5)) + (-Ifges(6,1) * t629 + Ifges(6,4) * t626) * qJD(1);
t634 = -mrSges(6,2) * t556 + pkin(8) * t645 + t625 * t540 + t628 * t539 + pkin(5) * t637 + mrSges(6,1) * t555 + Ifges(6,6) * t597 + Ifges(6,5) * t598 + Ifges(6,3) * qJDD(5) + (-t587 * t629 - t588 * t626) * qJD(1);
t586 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t629 + Ifges(6,6) * t626) * qJD(1);
t520 = mrSges(6,2) * t559 - mrSges(6,3) * t555 + Ifges(6,1) * t598 + Ifges(6,4) * t597 + Ifges(6,5) * qJDD(5) - pkin(8) * t536 - qJD(5) * t587 - t539 * t625 + t540 * t628 + t586 * t653;
t522 = -mrSges(6,1) * t559 + mrSges(6,3) * t556 + Ifges(6,4) * t598 + Ifges(6,2) * t597 + Ifges(6,6) * qJDD(5) - pkin(5) * t536 + qJD(5) * t588 + t586 * t652 - t635;
t633 = -mrSges(3,1) * t585 - mrSges(4,1) * t565 - mrSges(2,2) * t603 - mrSges(5,2) * t563 - pkin(2) * t519 + pkin(3) * t524 - qJ(4) * t636 - t629 * t520 + qJ(2) * (-mrSges(3,1) * t632 + t641) - pkin(1) * t518 + t626 * t522 + pkin(7) * t529 + mrSges(5,3) * t562 + mrSges(4,2) * t566 + mrSges(3,3) * t584 + mrSges(2,1) * t602 + (mrSges(5,3) * qJ(4) + Ifges(5,1) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t527 = m(5) * t618 - t654;
t525 = -m(3) * g(3) - t526;
t510 = mrSges(5,1) * t563 - mrSges(4,3) * t565 + pkin(4) * t529 - qJ(4) * t527 + qJDD(1) * t659 + t618 * t661 + t632 * t658 + t634;
t509 = mrSges(4,3) * t566 - mrSges(5,1) * t562 - t626 * t520 - t629 * t522 - pkin(4) * t665 + pkin(7) * t654 - pkin(3) * t527 - t659 * t632 + (-mrSges(4,1) + mrSges(5,2)) * t618 + t658 * qJDD(1);
t508 = mrSges(3,2) * t585 - mrSges(2,3) * t602 - qJ(2) * t525 - qJ(3) * t519 - t623 * t509 + t624 * t510 - t657 * t632 + t660 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t507 = mrSges(3,2) * t584 + mrSges(2,3) * t603 - pkin(1) * t525 + pkin(2) * t526 + g(3) * t662 - qJ(3) * t646 + qJDD(1) * t657 - t624 * t509 - t623 * t510 + t632 * t660;
t1 = [-m(1) * g(1) + t647; -m(1) * g(2) + t655; (-m(1) - m(2) - m(3)) * g(3) - t526; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t655 - t627 * t507 + t630 * t508; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t647 + t630 * t507 + t627 * t508; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t633; t633; t518; t526; t524; t634; t635;];
tauJB  = t1;
