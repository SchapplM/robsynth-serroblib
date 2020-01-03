% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:52
% EndTime: 2019-12-31 16:44:54
% DurationCPUTime: 1.22s
% Computational Cost: add. (8378->194), mult. (19826->237), div. (0->0), fcn. (12087->6), ass. (0->91)
t594 = Ifges(4,1) + Ifges(5,1);
t587 = Ifges(4,4) - Ifges(5,5);
t586 = Ifges(4,5) + Ifges(5,4);
t593 = Ifges(4,2) + Ifges(5,3);
t585 = Ifges(4,6) - Ifges(5,6);
t592 = Ifges(4,3) + Ifges(5,2);
t553 = qJD(1) ^ 2;
t549 = sin(qJ(3));
t548 = cos(pkin(6));
t590 = cos(qJ(3));
t569 = t548 * t590;
t547 = sin(pkin(6));
t574 = t547 * qJD(1);
t526 = -qJD(1) * t569 + t549 * t574;
t558 = t590 * t547 + t548 * t549;
t527 = t558 * qJD(1);
t508 = t526 * mrSges(5,1) - t527 * mrSges(5,3);
t576 = t526 * qJD(3);
t515 = t558 * qJDD(1) - t576;
t550 = sin(qJ(1));
t551 = cos(qJ(1));
t533 = -t551 * g(1) - t550 * g(2);
t528 = -t553 * pkin(1) + qJDD(1) * qJ(2) + t533;
t573 = qJD(1) * qJD(2);
t568 = -t548 * g(3) - 0.2e1 * t547 * t573;
t589 = pkin(2) * t548;
t496 = (-pkin(5) * qJDD(1) + t553 * t589 - t528) * t547 + t568;
t517 = -t547 * g(3) + (t528 + 0.2e1 * t573) * t548;
t571 = qJDD(1) * t548;
t543 = t548 ^ 2;
t582 = t543 * t553;
t497 = -pkin(2) * t582 + pkin(5) * t571 + t517;
t492 = t590 * t496 - t549 * t497;
t507 = t526 * pkin(3) - t527 * qJ(4);
t552 = qJD(3) ^ 2;
t491 = -qJDD(3) * pkin(3) - t552 * qJ(4) + t527 * t507 + qJDD(4) - t492;
t523 = -t526 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t564 = -m(5) * t491 + qJDD(3) * mrSges(5,1) + qJD(3) * t523;
t486 = t515 * mrSges(5,2) + t527 * t508 - t564;
t493 = t549 * t496 + t590 * t497;
t490 = -t552 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t526 * t507 + t493;
t572 = qJDD(1) * t547;
t575 = t527 * qJD(3);
t514 = -qJDD(1) * t569 + t549 * t572 + t575;
t522 = -qJD(3) * mrSges(5,1) + t527 * mrSges(5,2);
t570 = m(5) * t490 + qJDD(3) * mrSges(5,3) + qJD(3) * t522;
t578 = t586 * qJD(3) - t587 * t526 + t594 * t527;
t580 = -t585 * qJD(3) + t593 * t526 - t587 * t527;
t591 = t592 * qJDD(3) - t585 * t514 + t586 * t515 + t578 * t526 - t580 * t527 + mrSges(4,1) * t492 - mrSges(5,1) * t491 - mrSges(4,2) * t493 + mrSges(5,3) * t490 - pkin(3) * t486 + qJ(4) * (-t514 * mrSges(5,2) - t526 * t508 + t570);
t588 = -mrSges(4,3) - mrSges(5,2);
t583 = mrSges(3,2) * t547;
t521 = qJD(3) * mrSges(4,1) - t527 * mrSges(4,3);
t577 = -t526 * mrSges(4,1) - t527 * mrSges(4,2) - t508;
t482 = m(4) * t493 - qJDD(3) * mrSges(4,2) - qJD(3) * t521 + t588 * t514 + t577 * t526 + t570;
t520 = -qJD(3) * mrSges(4,2) - t526 * mrSges(4,3);
t483 = m(4) * t492 + qJDD(3) * mrSges(4,1) + qJD(3) * t520 + t588 * t515 + t577 * t527 + t564;
t474 = t549 * t482 + t590 * t483;
t516 = -t547 * t528 + t568;
t559 = mrSges(3,3) * qJDD(1) + t553 * (-mrSges(3,1) * t548 + t583);
t472 = m(3) * t516 - t559 * t547 + t474;
t565 = t590 * t482 - t549 * t483;
t473 = m(3) * t517 + t559 * t548 + t565;
t566 = -t547 * t472 + t548 * t473;
t464 = m(2) * t533 - t553 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t566;
t532 = t550 * g(1) - t551 * g(2);
t563 = qJDD(2) - t532;
t525 = -qJDD(1) * pkin(1) - t553 * qJ(2) + t563;
t542 = t547 ^ 2;
t513 = (-pkin(1) - t589) * qJDD(1) + (-qJ(2) + (-t542 - t543) * pkin(5)) * t553 + t563;
t488 = -0.2e1 * qJD(4) * t527 + (-t515 + t576) * qJ(4) + (t514 + t575) * pkin(3) + t513;
t484 = m(5) * t488 + t514 * mrSges(5,1) - t515 * mrSges(5,3) - t527 * t522 + t526 * t523;
t556 = m(4) * t513 + t514 * mrSges(4,1) + t515 * mrSges(4,2) + t526 * t520 + t527 * t521 + t484;
t554 = -m(3) * t525 + mrSges(3,1) * t571 - t556 + (t542 * t553 + t582) * mrSges(3,3);
t476 = (mrSges(2,1) - t583) * qJDD(1) + t554 + m(2) * t532 - t553 * mrSges(2,2);
t581 = t550 * t464 + t551 * t476;
t466 = t548 * t472 + t547 * t473;
t579 = -t592 * qJD(3) + t585 * t526 - t586 * t527;
t567 = t551 * t464 - t550 * t476;
t562 = Ifges(3,1) * t547 + Ifges(3,4) * t548;
t561 = Ifges(3,4) * t547 + Ifges(3,2) * t548;
t560 = Ifges(3,5) * t547 + Ifges(3,6) * t548;
t467 = -mrSges(4,1) * t513 - mrSges(5,1) * t488 + mrSges(5,2) * t490 + mrSges(4,3) * t493 - pkin(3) * t484 + t578 * qJD(3) + t585 * qJDD(3) - t593 * t514 + t587 * t515 + t579 * t527;
t468 = mrSges(4,2) * t513 + mrSges(5,2) * t491 - mrSges(4,3) * t492 - mrSges(5,3) * t488 - qJ(4) * t484 + t580 * qJD(3) + t586 * qJDD(3) - t587 * t514 + t594 * t515 + t579 * t526;
t530 = t560 * qJD(1);
t459 = -mrSges(3,1) * t525 + mrSges(3,3) * t517 - pkin(2) * t556 + pkin(5) * t565 + t561 * qJDD(1) + t590 * t467 + t549 * t468 - t530 * t574;
t461 = t548 * qJD(1) * t530 + mrSges(3,2) * t525 - mrSges(3,3) * t516 - pkin(5) * t474 + t562 * qJDD(1) - t549 * t467 + t590 * t468;
t478 = mrSges(3,2) * t572 - t554;
t557 = mrSges(2,1) * t532 - mrSges(2,2) * t533 + Ifges(2,3) * qJDD(1) - pkin(1) * t478 + qJ(2) * t566 + t548 * t459 + t547 * t461;
t457 = -pkin(1) * t466 - pkin(2) * t474 - mrSges(3,1) * t516 + mrSges(3,2) * t517 + mrSges(2,1) * g(3) + mrSges(2,3) * t533 + (Ifges(2,6) - t560) * qJDD(1) + (-t547 * t561 + t548 * t562 + Ifges(2,5)) * t553 - t591;
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t532 + Ifges(2,5) * qJDD(1) - t553 * Ifges(2,6) - qJ(2) * t466 - t547 * t459 + t548 * t461;
t1 = [-m(1) * g(1) + t567; -m(1) * g(2) + t581; (-m(1) - m(2)) * g(3) + t466; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t581 + t551 * t456 - t550 * t457; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t567 + t550 * t456 + t551 * t457; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t557; t557; t478; t591; t486;];
tauJB = t1;
