% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPP3
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:30
% EndTime: 2019-12-31 16:57:32
% DurationCPUTime: 1.37s
% Computational Cost: add. (9890->216), mult. (22387->267), div. (0->0), fcn. (12843->6), ass. (0->88)
t580 = Ifges(4,1) + Ifges(5,1);
t573 = Ifges(4,4) - Ifges(5,5);
t572 = Ifges(4,5) + Ifges(5,4);
t579 = -Ifges(4,2) - Ifges(5,3);
t578 = -Ifges(5,2) - Ifges(4,3);
t571 = Ifges(4,6) - Ifges(5,6);
t543 = sin(qJ(2));
t545 = cos(qJ(2));
t561 = qJD(1) * qJD(2);
t528 = t543 * qJDD(1) + t545 * t561;
t544 = sin(qJ(1));
t546 = cos(qJ(1));
t535 = -t546 * g(1) - t544 * g(2);
t548 = qJD(1) ^ 2;
t523 = -t548 * pkin(1) + qJDD(1) * pkin(5) + t535;
t569 = t543 * t523;
t575 = pkin(2) * t548;
t486 = qJDD(2) * pkin(2) - t528 * qJ(3) - t569 + (qJ(3) * t561 + t543 * t575 - g(3)) * t545;
t509 = -t543 * g(3) + t545 * t523;
t529 = t545 * qJDD(1) - t543 * t561;
t563 = qJD(1) * t543;
t531 = qJD(2) * pkin(2) - qJ(3) * t563;
t541 = t545 ^ 2;
t487 = t529 * qJ(3) - qJD(2) * t531 - t541 * t575 + t509;
t542 = sin(pkin(6));
t562 = qJD(1) * t545;
t570 = cos(pkin(6));
t516 = t542 * t563 - t570 * t562;
t576 = -2 * qJD(3);
t483 = t542 * t486 + t570 * t487 + t516 * t576;
t504 = t542 * t528 - t570 * t529;
t517 = (t542 * t545 + t570 * t543) * qJD(1);
t511 = qJD(2) * mrSges(4,1) - t517 * mrSges(4,3);
t498 = t516 * pkin(3) - t517 * qJ(4);
t547 = qJD(2) ^ 2;
t478 = -t547 * pkin(3) + qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) - t516 * t498 + t483;
t512 = -qJD(2) * mrSges(5,1) + t517 * mrSges(5,2);
t559 = m(5) * t478 + qJDD(2) * mrSges(5,3) + qJD(2) * t512;
t499 = t516 * mrSges(5,1) - t517 * mrSges(5,3);
t564 = -t516 * mrSges(4,1) - t517 * mrSges(4,2) - t499;
t574 = -mrSges(4,3) - mrSges(5,2);
t471 = m(4) * t483 - qJDD(2) * mrSges(4,2) - qJD(2) * t511 + t574 * t504 + t564 * t516 + t559;
t552 = t570 * t486 - t542 * t487;
t482 = t517 * t576 + t552;
t505 = t570 * t528 + t542 * t529;
t510 = -qJD(2) * mrSges(4,2) - t516 * mrSges(4,3);
t479 = -qJDD(2) * pkin(3) - t547 * qJ(4) + qJDD(4) + ((2 * qJD(3)) + t498) * t517 - t552;
t513 = -t516 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t555 = -m(5) * t479 + qJDD(2) * mrSges(5,1) + qJD(2) * t513;
t472 = m(4) * t482 + qJDD(2) * mrSges(4,1) + qJD(2) * t510 + t574 * t505 + t564 * t517 + t555;
t464 = t542 * t471 + t570 * t472;
t476 = t505 * mrSges(5,2) + t517 * t499 - t555;
t508 = -t545 * g(3) - t569;
t519 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t543 + Ifges(3,2) * t545) * qJD(1);
t520 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t543 + Ifges(3,4) * t545) * qJD(1);
t565 = t572 * qJD(2) - t573 * t516 + t580 * t517;
t566 = t571 * qJD(2) + t579 * t516 + t573 * t517;
t577 = (t543 * t519 - t545 * t520) * qJD(1) + (Ifges(3,3) - t578) * qJDD(2) - t571 * t504 + t572 * t505 + t565 * t516 + t566 * t517 + mrSges(3,1) * t508 + mrSges(4,1) * t482 - mrSges(5,1) * t479 - mrSges(3,2) * t509 - mrSges(4,2) * t483 + mrSges(5,3) * t478 + Ifges(3,5) * t528 + Ifges(3,6) * t529 + pkin(2) * t464 - pkin(3) * t476 + qJ(4) * (-t504 * mrSges(5,2) - t516 * t499 + t559);
t527 = (-mrSges(3,1) * t545 + mrSges(3,2) * t543) * qJD(1);
t533 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t562;
t462 = m(3) * t508 + qJDD(2) * mrSges(3,1) - t528 * mrSges(3,3) + qJD(2) * t533 - t527 * t563 + t464;
t532 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t563;
t556 = t570 * t471 - t542 * t472;
t463 = m(3) * t509 - qJDD(2) * mrSges(3,2) + t529 * mrSges(3,3) - qJD(2) * t532 + t527 * t562 + t556;
t557 = -t543 * t462 + t545 * t463;
t454 = m(2) * t535 - t548 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t557;
t534 = t544 * g(1) - t546 * g(2);
t553 = -qJDD(1) * pkin(1) - t534;
t488 = -t529 * pkin(2) + qJDD(3) + t531 * t563 + (-qJ(3) * t541 - pkin(5)) * t548 + t553;
t481 = -0.2e1 * qJD(4) * t517 + (qJD(2) * t516 - t505) * qJ(4) + (qJD(2) * t517 + t504) * pkin(3) + t488;
t474 = m(5) * t481 + t504 * mrSges(5,1) - t505 * mrSges(5,3) - t517 * t512 + t516 * t513;
t473 = m(4) * t488 + t504 * mrSges(4,1) + t505 * mrSges(4,2) + t516 * t510 + t517 * t511 + t474;
t522 = -t548 * pkin(5) + t553;
t550 = -m(3) * t522 + t529 * mrSges(3,1) - t528 * mrSges(3,2) - t532 * t563 + t533 * t562 - t473;
t466 = m(2) * t534 + qJDD(1) * mrSges(2,1) - t548 * mrSges(2,2) + t550;
t568 = t544 * t454 + t546 * t466;
t456 = t545 * t462 + t543 * t463;
t567 = t578 * qJD(2) + t571 * t516 - t572 * t517;
t558 = t546 * t454 - t544 * t466;
t457 = -mrSges(4,1) * t488 - mrSges(5,1) * t481 + mrSges(5,2) * t478 + mrSges(4,3) * t483 - pkin(3) * t474 + t565 * qJD(2) + t571 * qJDD(2) + t579 * t504 + t573 * t505 + t567 * t517;
t458 = mrSges(4,2) * t488 + mrSges(5,2) * t479 - mrSges(4,3) * t482 - mrSges(5,3) * t481 - qJ(4) * t474 - t566 * qJD(2) + t572 * qJDD(2) - t573 * t504 + t580 * t505 + t567 * t516;
t518 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t543 + Ifges(3,6) * t545) * qJD(1);
t449 = -mrSges(3,1) * t522 + mrSges(3,3) * t509 + Ifges(3,4) * t528 + Ifges(3,2) * t529 + Ifges(3,6) * qJDD(2) - pkin(2) * t473 + qJ(3) * t556 + qJD(2) * t520 + t570 * t457 + t542 * t458 - t518 * t563;
t451 = mrSges(3,2) * t522 - mrSges(3,3) * t508 + Ifges(3,1) * t528 + Ifges(3,4) * t529 + Ifges(3,5) * qJDD(2) - qJ(3) * t464 - qJD(2) * t519 - t542 * t457 + t570 * t458 + t518 * t562;
t551 = mrSges(2,1) * t534 - mrSges(2,2) * t535 + Ifges(2,3) * qJDD(1) + pkin(1) * t550 + pkin(5) * t557 + t545 * t449 + t543 * t451;
t447 = mrSges(2,1) * g(3) + mrSges(2,3) * t535 + t548 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t456 - t577;
t446 = -mrSges(2,2) * g(3) - mrSges(2,3) * t534 + Ifges(2,5) * qJDD(1) - t548 * Ifges(2,6) - pkin(5) * t456 - t543 * t449 + t545 * t451;
t1 = [-m(1) * g(1) + t558; -m(1) * g(2) + t568; (-m(1) - m(2)) * g(3) + t456; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t568 + t546 * t446 - t544 * t447; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t558 + t544 * t446 + t546 * t447; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t551; t551; t577; t473; t476;];
tauJB = t1;
