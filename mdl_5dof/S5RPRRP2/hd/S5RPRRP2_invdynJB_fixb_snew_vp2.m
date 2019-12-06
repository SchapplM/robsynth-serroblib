% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:35
% EndTime: 2019-12-05 18:01:37
% DurationCPUTime: 1.75s
% Computational Cost: add. (21678->207), mult. (29571->247), div. (0->0), fcn. (14363->8), ass. (0->94)
t630 = Ifges(5,1) + Ifges(6,1);
t622 = Ifges(5,4) + Ifges(6,4);
t621 = Ifges(5,5) + Ifges(6,5);
t629 = Ifges(5,2) + Ifges(6,2);
t620 = Ifges(5,6) + Ifges(6,6);
t628 = Ifges(5,3) + Ifges(6,3);
t581 = qJD(1) + qJD(3);
t589 = sin(qJ(4));
t592 = cos(qJ(4));
t553 = (-mrSges(6,1) * t592 + mrSges(6,2) * t589) * t581;
t580 = qJDD(1) + qJDD(3);
t612 = qJD(4) * t581;
t607 = t592 * t612;
t555 = t589 * t580 + t607;
t591 = sin(qJ(1));
t594 = cos(qJ(1));
t570 = t594 * g(2) + t591 * g(3);
t561 = qJDD(1) * pkin(1) + t570;
t569 = t591 * g(2) - t594 * g(3);
t595 = qJD(1) ^ 2;
t562 = -t595 * pkin(1) + t569;
t587 = sin(pkin(8));
t588 = cos(pkin(8));
t539 = t588 * t561 - t587 * t562;
t536 = qJDD(1) * pkin(2) + t539;
t540 = t587 * t561 + t588 * t562;
t537 = -t595 * pkin(2) + t540;
t590 = sin(qJ(3));
t593 = cos(qJ(3));
t532 = t590 * t536 + t593 * t537;
t579 = t581 ^ 2;
t529 = -t579 * pkin(3) + t580 * pkin(7) + t532;
t586 = -g(1) + qJDD(2);
t572 = t592 * t586;
t611 = qJD(5) * t581;
t624 = pkin(4) * t579;
t522 = qJDD(4) * pkin(4) + t572 + (-t555 + t607) * qJ(5) + (t592 * t624 - t529 - 0.2e1 * t611) * t589;
t617 = t581 * t592;
t566 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t617;
t610 = m(6) * t522 + qJDD(4) * mrSges(6,1) + qJD(4) * t566;
t618 = t581 * t589;
t519 = -t555 * mrSges(6,3) - t553 * t618 + t610;
t526 = t592 * t529 + t589 * t586;
t556 = t592 * t580 - t589 * t612;
t563 = qJD(4) * pkin(4) - qJ(5) * t618;
t585 = t592 ^ 2;
t523 = t556 * qJ(5) - qJD(4) * t563 - t585 * t624 + 0.2e1 * t592 * t611 + t526;
t525 = -t589 * t529 + t572;
t614 = (t630 * t589 + t622 * t592) * t581 + t621 * qJD(4);
t615 = (t622 * t589 + t629 * t592) * t581 + t620 * qJD(4);
t627 = mrSges(5,1) * t525 + mrSges(6,1) * t522 - mrSges(5,2) * t526 - mrSges(6,2) * t523 + pkin(4) * t519 + t628 * qJDD(4) + t621 * t555 + t620 * t556 + (t615 * t589 - t614 * t592) * t581;
t623 = -mrSges(5,2) - mrSges(6,2);
t554 = (-mrSges(5,1) * t592 + mrSges(5,2) * t589) * t581;
t567 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t617;
t516 = m(5) * t525 + qJDD(4) * mrSges(5,1) + qJD(4) * t567 + (-t553 - t554) * t618 + (-mrSges(5,3) - mrSges(6,3)) * t555 + t610;
t609 = m(6) * t523 + t556 * mrSges(6,3) + t553 * t617;
t564 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t618;
t613 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t618 - t564;
t517 = m(5) * t526 + t556 * mrSges(5,3) + t613 * qJD(4) + t623 * qJDD(4) + t554 * t617 + t609;
t603 = -t589 * t516 + t592 * t517;
t506 = m(4) * t532 - t579 * mrSges(4,1) - t580 * mrSges(4,2) + t603;
t531 = t593 * t536 - t590 * t537;
t600 = -t580 * pkin(3) - t531;
t528 = -t579 * pkin(7) + t600;
t524 = t563 * t618 - t556 * pkin(4) + qJDD(5) + (-qJ(5) * t585 - pkin(7)) * t579 + t600;
t602 = -m(6) * t524 + t556 * mrSges(6,1) + t566 * t617;
t597 = -m(5) * t528 + t556 * mrSges(5,1) + t623 * t555 + t567 * t617 + t613 * t618 + t602;
t511 = m(4) * t531 + t580 * mrSges(4,1) - t579 * mrSges(4,2) + t597;
t499 = t590 * t506 + t593 * t511;
t496 = m(3) * t539 + qJDD(1) * mrSges(3,1) - t595 * mrSges(3,2) + t499;
t604 = t593 * t506 - t590 * t511;
t497 = m(3) * t540 - t595 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t604;
t491 = t588 * t496 + t587 * t497;
t509 = t592 * t516 + t589 * t517;
t616 = (t621 * t589 + t620 * t592) * t581 + t628 * qJD(4);
t608 = m(4) * t586 + t509;
t605 = -t587 * t496 + t588 * t497;
t488 = m(2) * t569 - t595 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t605;
t489 = m(2) * t570 + qJDD(1) * mrSges(2,1) - t595 * mrSges(2,2) + t491;
t606 = t594 * t488 - t591 * t489;
t507 = m(3) * t586 + t608;
t601 = -t591 * t488 - t594 * t489;
t518 = t555 * mrSges(6,2) + t564 * t618 - t602;
t501 = -mrSges(5,1) * t528 + mrSges(5,3) * t526 - mrSges(6,1) * t524 + mrSges(6,3) * t523 - pkin(4) * t518 + qJ(5) * t609 - t616 * t618 + t629 * t556 + t622 * t555 + (-qJ(5) * mrSges(6,2) + t620) * qJDD(4) + (-qJ(5) * t564 + t614) * qJD(4);
t503 = mrSges(5,2) * t528 + mrSges(6,2) * t524 - mrSges(5,3) * t525 - mrSges(6,3) * t522 - qJ(5) * t519 - t615 * qJD(4) + t621 * qJDD(4) + t630 * t555 + t622 * t556 + t616 * t617;
t599 = mrSges(4,1) * t531 - mrSges(4,2) * t532 + Ifges(4,3) * t580 + pkin(3) * t597 + pkin(7) * t603 + t592 * t501 + t589 * t503;
t596 = mrSges(2,1) * t570 + mrSges(3,1) * t539 - mrSges(2,2) * t569 - mrSges(3,2) * t540 + pkin(1) * t491 + pkin(2) * t499 + t599 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t492 = -mrSges(4,1) * t586 + mrSges(4,3) * t532 + t579 * Ifges(4,5) + Ifges(4,6) * t580 - pkin(3) * t509 - t627;
t486 = mrSges(4,2) * t586 - mrSges(4,3) * t531 + Ifges(4,5) * t580 - t579 * Ifges(4,6) - pkin(7) * t509 - t589 * t501 + t592 * t503;
t485 = mrSges(3,2) * t586 - mrSges(3,3) * t539 + Ifges(3,5) * qJDD(1) - t595 * Ifges(3,6) - pkin(6) * t499 + t593 * t486 - t590 * t492;
t484 = -mrSges(3,1) * t586 + mrSges(3,3) * t540 + t595 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t608 + pkin(6) * t604 + t590 * t486 + t593 * t492;
t483 = -mrSges(2,2) * g(1) - mrSges(2,3) * t570 + Ifges(2,5) * qJDD(1) - t595 * Ifges(2,6) - qJ(2) * t491 - t587 * t484 + t588 * t485;
t482 = mrSges(2,1) * g(1) + mrSges(2,3) * t569 + t595 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t507 + qJ(2) * t605 + t588 * t484 + t587 * t485;
t1 = [(-m(1) - m(2)) * g(1) + t507; -m(1) * g(2) + t601; -m(1) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t596; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t606 - t594 * t482 - t591 * t483; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t601 - t591 * t482 + t594 * t483; t596; t507; t599; t627; t518;];
tauJB = t1;
