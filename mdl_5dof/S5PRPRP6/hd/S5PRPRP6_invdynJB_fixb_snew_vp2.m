% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:29
% EndTime: 2019-12-05 15:40:31
% DurationCPUTime: 1.13s
% Computational Cost: add. (5644->188), mult. (9806->223), div. (0->0), fcn. (4440->6), ass. (0->83)
t567 = Ifges(5,1) + Ifges(6,1);
t554 = Ifges(5,4) - Ifges(6,5);
t565 = Ifges(6,4) + Ifges(5,5);
t566 = Ifges(5,2) + Ifges(6,3);
t564 = Ifges(5,6) - Ifges(6,6);
t563 = Ifges(5,3) + Ifges(6,2);
t524 = sin(qJ(4));
t526 = cos(qJ(4));
t562 = -t564 * qJD(4) + (t566 * t524 - t554 * t526) * qJD(2);
t561 = t565 * qJD(4) + (-t554 * t524 + t567 * t526) * qJD(2);
t523 = sin(pkin(7));
t550 = cos(pkin(7));
t505 = -t550 * g(1) - t523 * g(2);
t520 = -g(3) + qJDD(1);
t525 = sin(qJ(2));
t527 = cos(qJ(2));
t475 = t527 * t505 + t525 * t520;
t560 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t475;
t529 = qJD(2) ^ 2;
t557 = -pkin(2) - pkin(6);
t470 = t557 * t529 - t560;
t544 = qJD(2) * qJD(4);
t500 = t524 * qJDD(2) + t526 * t544;
t501 = t526 * qJDD(2) - t524 * t544;
t460 = t500 * pkin(4) - t501 * qJ(5) + (-0.2e1 * qJD(5) * t526 + (pkin(4) * t526 + qJ(5) * t524) * qJD(4)) * qJD(2) + t470;
t545 = qJD(2) * t526;
t508 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t545;
t546 = qJD(2) * t524;
t509 = -mrSges(6,2) * t546 + qJD(4) * mrSges(6,3);
t455 = m(6) * t460 + t500 * mrSges(6,1) - t501 * mrSges(6,3) - t508 * t545 + t509 * t546;
t474 = -t525 * t505 + t527 * t520;
t534 = -t529 * qJ(3) + qJDD(3) - t474;
t471 = t557 * qJDD(2) + t534;
t504 = t523 * g(1) - t550 * g(2);
t467 = t524 * t471 - t526 * t504;
t497 = (pkin(4) * t524 - qJ(5) * t526) * qJD(2);
t528 = qJD(4) ^ 2;
t463 = -t528 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t497 * t546 + t467;
t537 = qJD(2) * (-t563 * qJD(4) + (t564 * t524 - t565 * t526) * qJD(2));
t438 = -mrSges(5,1) * t470 - mrSges(6,1) * t460 + mrSges(6,2) * t463 + mrSges(5,3) * t467 - pkin(4) * t455 + t561 * qJD(4) + t564 * qJDD(4) - t566 * t500 + t554 * t501 + t526 * t537;
t549 = t524 * t504;
t464 = -qJDD(4) * pkin(4) - t528 * qJ(5) - t549 + qJDD(5) + (qJD(2) * t497 - t471) * t526;
t466 = t526 * t471 + t549;
t439 = mrSges(5,2) * t470 + mrSges(6,2) * t464 - mrSges(5,3) * t466 - mrSges(6,3) * t460 - qJ(5) * t455 + t562 * qJD(4) + t565 * qJDD(4) - t554 * t500 + t567 * t501 + t524 * t537;
t507 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t545;
t498 = (mrSges(6,1) * t524 - mrSges(6,3) * t526) * qJD(2);
t536 = qJD(2) * (-t498 - (mrSges(5,1) * t524 + mrSges(5,2) * t526) * qJD(2));
t542 = m(6) * t463 + qJDD(4) * mrSges(6,3) + qJD(4) * t508;
t555 = -mrSges(5,3) - mrSges(6,2);
t453 = m(5) * t467 - qJDD(4) * mrSges(5,2) - qJD(4) * t507 + t555 * t500 + t524 * t536 + t542;
t506 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t546;
t535 = -m(6) * t464 + qJDD(4) * mrSges(6,1) + qJD(4) * t509;
t454 = m(5) * t466 + qJDD(4) * mrSges(5,1) + qJD(4) * t506 + t555 * t501 + t526 * t536 + t535;
t446 = t524 * t453 + t526 * t454;
t473 = -qJDD(2) * pkin(2) + t534;
t533 = -m(4) * t473 + t529 * mrSges(4,3) - t446;
t443 = qJDD(2) * mrSges(4,2) - t533;
t472 = t529 * pkin(2) + t560;
t532 = -m(5) * t470 - t500 * mrSges(5,1) - t501 * mrSges(5,2) - t506 * t546 - t507 * t545 - t455;
t450 = -m(4) * t472 + t529 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t532;
t559 = mrSges(3,1) * t474 - mrSges(3,2) * t475 + mrSges(4,2) * t473 - mrSges(4,3) * t472 - pkin(2) * t443 - pkin(6) * t446 + qJ(3) * t450 - t524 * t438 + t526 * t439 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t558 = m(3) + m(4);
t556 = mrSges(3,1) - mrSges(4,2);
t553 = Ifges(3,5) - Ifges(4,4);
t552 = -Ifges(3,6) + Ifges(4,5);
t441 = m(3) * t474 - t529 * mrSges(3,2) + t556 * qJDD(2) + t533;
t449 = m(3) * t475 - t529 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t450;
t538 = -t525 * t441 + t527 * t449;
t436 = m(2) * t505 + t538;
t547 = t526 * t453 - t524 * t454;
t444 = (m(2) + t558) * t504 - t547;
t548 = t523 * t436 + t550 * t444;
t437 = t527 * t441 + t525 * t449;
t540 = m(2) * t520 + t437;
t539 = t550 * t436 - t523 * t444;
t458 = t501 * mrSges(6,2) + t498 * t545 - t535;
t531 = mrSges(5,1) * t466 - mrSges(6,1) * t464 - mrSges(5,2) * t467 + mrSges(6,3) * t463 - pkin(4) * t458 + qJ(5) * t542 + (-qJ(5) * t498 + t561) * t546 - t562 * t545 + t565 * t501 + (-qJ(5) * mrSges(6,2) - t564) * t500 + t563 * qJDD(4);
t445 = -m(4) * t504 + t547;
t433 = mrSges(4,1) * t473 - mrSges(3,3) * t474 - qJ(3) * t445 + pkin(3) * t446 + (-mrSges(3,2) + mrSges(4,3)) * t504 + t552 * t529 + t531 + t553 * qJDD(2);
t432 = -mrSges(4,1) * t472 + mrSges(3,3) * t475 - pkin(2) * t445 - pkin(3) * t532 - pkin(6) * t547 - t552 * qJDD(2) - t526 * t438 - t524 * t439 + t556 * t504 + t553 * t529;
t431 = -mrSges(2,1) * t520 + mrSges(2,3) * t505 - pkin(1) * t437 - t559;
t430 = mrSges(2,2) * t520 - mrSges(2,3) * t504 - pkin(5) * t437 - t525 * t432 + t527 * t433;
t1 = [-m(1) * g(1) + t539; -m(1) * g(2) + t548; -m(1) * g(3) + t540; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t548 + t550 * t430 - t523 * t431; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t539 + t523 * t430 + t550 * t431; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t505 + t525 * t433 + t527 * t432 - pkin(1) * t547 + pkin(5) * t538 + (pkin(1) * t558 + mrSges(2,1)) * t504; t540; t559; t443; t531; t458;];
tauJB = t1;
