% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:46
% EndTime: 2019-12-05 16:15:49
% DurationCPUTime: 2.60s
% Computational Cost: add. (36948->195), mult. (52874->250), div. (0->0), fcn. (34717->10), ass. (0->98)
t552 = qJD(2) + qJD(3);
t548 = t552 ^ 2;
t557 = cos(pkin(9));
t593 = pkin(4) * t557;
t555 = sin(pkin(9));
t592 = mrSges(5,2) * t555;
t551 = t557 ^ 2;
t591 = t548 * t551;
t549 = qJDD(2) + qJDD(3);
t590 = t549 * t557;
t576 = Ifges(5,5) * t555 + Ifges(5,6) * t557;
t589 = t548 * t576;
t556 = sin(pkin(8));
t558 = cos(pkin(8));
t535 = t556 * g(1) - t558 * g(2);
t536 = -t558 * g(1) - t556 * g(2);
t561 = sin(qJ(2));
t564 = cos(qJ(2));
t524 = t564 * t535 - t561 * t536;
t521 = qJDD(2) * pkin(2) + t524;
t525 = t561 * t535 + t564 * t536;
t565 = qJD(2) ^ 2;
t522 = -t565 * pkin(2) + t525;
t560 = sin(qJ(3));
t563 = cos(qJ(3));
t509 = t560 * t521 + t563 * t522;
t506 = -t548 * pkin(3) + t549 * qJ(4) + t509;
t554 = -g(3) + qJDD(1);
t586 = qJD(4) * t552;
t587 = t557 * t554 - 0.2e1 * t555 * t586;
t499 = (-pkin(7) * t549 + t548 * t593 - t506) * t555 + t587;
t503 = t555 * t554 + (t506 + 0.2e1 * t586) * t557;
t500 = -pkin(4) * t591 + pkin(7) * t590 + t503;
t559 = sin(qJ(5));
t562 = cos(qJ(5));
t497 = t562 * t499 - t559 * t500;
t571 = -t555 * t559 + t557 * t562;
t528 = t571 * t552;
t572 = t555 * t562 + t557 * t559;
t529 = t572 * t552;
t515 = -t528 * mrSges(6,1) + t529 * mrSges(6,2);
t517 = t528 * qJD(5) + t572 * t549;
t526 = -qJD(5) * mrSges(6,2) + t528 * mrSges(6,3);
t495 = m(6) * t497 + qJDD(5) * mrSges(6,1) - t517 * mrSges(6,3) + qJD(5) * t526 - t529 * t515;
t498 = t559 * t499 + t562 * t500;
t516 = -t529 * qJD(5) + t571 * t549;
t527 = qJD(5) * mrSges(6,1) - t529 * mrSges(6,3);
t496 = m(6) * t498 - qJDD(5) * mrSges(6,2) + t516 * mrSges(6,3) - qJD(5) * t527 + t528 * t515;
t485 = t562 * t495 + t559 * t496;
t502 = -t555 * t506 + t587;
t573 = mrSges(5,3) * t549 + (-mrSges(5,1) * t557 + t592) * t548;
t483 = m(5) * t502 - t573 * t555 + t485;
t580 = -t559 * t495 + t562 * t496;
t484 = m(5) * t503 + t573 * t557 + t580;
t581 = -t555 * t483 + t557 * t484;
t477 = m(4) * t509 - t548 * mrSges(4,1) - t549 * mrSges(4,2) + t581;
t508 = t563 * t521 - t560 * t522;
t575 = qJDD(4) - t508;
t505 = -t549 * pkin(3) - t548 * qJ(4) + t575;
t550 = t555 ^ 2;
t501 = (-pkin(3) - t593) * t549 + (-qJ(4) + (-t550 - t551) * pkin(7)) * t548 + t575;
t569 = m(6) * t501 - t516 * mrSges(6,1) + t517 * mrSges(6,2) - t528 * t526 + t529 * t527;
t567 = -m(5) * t505 + mrSges(5,1) * t590 - t569 + (t548 * t550 + t591) * mrSges(5,3);
t489 = m(4) * t508 - t548 * mrSges(4,2) + (mrSges(4,1) - t592) * t549 + t567;
t472 = t560 * t477 + t563 * t489;
t467 = m(3) * t524 + qJDD(2) * mrSges(3,1) - t565 * mrSges(3,2) + t472;
t582 = t563 * t477 - t560 * t489;
t468 = m(3) * t525 - t565 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t582;
t462 = t564 * t467 + t561 * t468;
t460 = m(2) * t535 + t462;
t583 = -t561 * t467 + t564 * t468;
t461 = m(2) * t536 + t583;
t588 = t558 * t460 + t556 * t461;
t479 = t557 * t483 + t555 * t484;
t585 = m(4) * t554 + t479;
t584 = -t556 * t460 + t558 * t461;
t579 = m(3) * t554 + t585;
t578 = Ifges(5,1) * t555 + Ifges(5,4) * t557;
t577 = Ifges(5,4) * t555 + Ifges(5,2) * t557;
t574 = m(2) * t554 + t579;
t510 = Ifges(6,5) * t529 + Ifges(6,6) * t528 + Ifges(6,3) * qJD(5);
t512 = Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * qJD(5);
t486 = -mrSges(6,1) * t501 + mrSges(6,3) * t498 + Ifges(6,4) * t517 + Ifges(6,2) * t516 + Ifges(6,6) * qJDD(5) + qJD(5) * t512 - t529 * t510;
t511 = Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * qJD(5);
t487 = mrSges(6,2) * t501 - mrSges(6,3) * t497 + Ifges(6,1) * t517 + Ifges(6,4) * t516 + Ifges(6,5) * qJDD(5) - qJD(5) * t511 + t528 * t510;
t470 = -mrSges(5,1) * t505 + mrSges(5,3) * t503 - pkin(4) * t569 + pkin(7) * t580 + t562 * t486 + t559 * t487 + t577 * t549 - t555 * t589;
t474 = mrSges(5,2) * t505 - mrSges(5,3) * t502 - pkin(7) * t485 - t559 * t486 + t562 * t487 + t578 * t549 + t557 * t589;
t491 = t549 * t592 - t567;
t570 = mrSges(4,1) * t508 - mrSges(4,2) * t509 + Ifges(4,3) * t549 - pkin(3) * t491 + qJ(4) * t581 + t557 * t470 + t555 * t474;
t568 = mrSges(3,1) * t524 - mrSges(3,2) * t525 + Ifges(3,3) * qJDD(2) + pkin(2) * t472 + t570;
t566 = mrSges(6,1) * t497 - mrSges(6,2) * t498 + Ifges(6,5) * t517 + Ifges(6,6) * t516 + Ifges(6,3) * qJDD(5) + t529 * t511 - t528 * t512;
t463 = -mrSges(4,1) * t554 - mrSges(5,1) * t502 + mrSges(5,2) * t503 + mrSges(4,3) * t509 - pkin(3) * t479 - pkin(4) * t485 + (Ifges(4,6) - t576) * t549 - t566 + (-t555 * t577 + t557 * t578 + Ifges(4,5)) * t548;
t456 = mrSges(4,2) * t554 - mrSges(4,3) * t508 + Ifges(4,5) * t549 - t548 * Ifges(4,6) - qJ(4) * t479 - t555 * t470 + t557 * t474;
t455 = mrSges(3,2) * t554 - mrSges(3,3) * t524 + Ifges(3,5) * qJDD(2) - t565 * Ifges(3,6) - pkin(6) * t472 + t563 * t456 - t560 * t463;
t454 = -mrSges(3,1) * t554 + mrSges(3,3) * t525 + t565 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t585 + pkin(6) * t582 + t560 * t456 + t563 * t463;
t453 = mrSges(2,2) * t554 - mrSges(2,3) * t535 - pkin(5) * t462 - t561 * t454 + t564 * t455;
t452 = -mrSges(2,1) * t554 + mrSges(2,3) * t536 - pkin(1) * t579 + pkin(5) * t583 + t564 * t454 + t561 * t455;
t1 = [-m(1) * g(1) + t584; -m(1) * g(2) + t588; -m(1) * g(3) + t574; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t588 - t556 * t452 + t558 * t453; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t584 + t558 * t452 + t556 * t453; -mrSges(1,1) * g(2) + mrSges(2,1) * t535 + mrSges(1,2) * g(1) - mrSges(2,2) * t536 + pkin(1) * t462 + t568; t574; t568; t570; t491; t566;];
tauJB = t1;
