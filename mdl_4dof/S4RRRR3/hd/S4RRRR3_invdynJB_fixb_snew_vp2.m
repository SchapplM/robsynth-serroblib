% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:30
% EndTime: 2019-12-31 17:24:32
% DurationCPUTime: 2.83s
% Computational Cost: add. (30004->242), mult. (64917->309), div. (0->0), fcn. (42692->8), ass. (0->99)
t569 = sin(qJ(2));
t573 = cos(qJ(2));
t589 = qJD(1) * qJD(2);
t550 = t569 * qJDD(1) + t573 * t589;
t570 = sin(qJ(1));
t574 = cos(qJ(1));
t557 = -t574 * g(1) - t570 * g(2);
t575 = qJD(1) ^ 2;
t545 = -t575 * pkin(1) + qJDD(1) * pkin(5) + t557;
t593 = t569 * t545;
t594 = pkin(2) * t575;
t515 = qJDD(2) * pkin(2) - t550 * pkin(6) - t593 + (pkin(6) * t589 + t569 * t594 - g(3)) * t573;
t533 = -t569 * g(3) + t573 * t545;
t551 = t573 * qJDD(1) - t569 * t589;
t591 = qJD(1) * t569;
t555 = qJD(2) * pkin(2) - pkin(6) * t591;
t566 = t573 ^ 2;
t516 = t551 * pkin(6) - qJD(2) * t555 - t566 * t594 + t533;
t568 = sin(qJ(3));
t572 = cos(qJ(3));
t503 = t572 * t515 - t568 * t516;
t542 = (-t568 * t569 + t572 * t573) * qJD(1);
t521 = t542 * qJD(3) + t572 * t550 + t568 * t551;
t543 = (t568 * t573 + t569 * t572) * qJD(1);
t563 = qJDD(2) + qJDD(3);
t564 = qJD(2) + qJD(3);
t491 = (t542 * t564 - t521) * pkin(7) + (t542 * t543 + t563) * pkin(3) + t503;
t504 = t568 * t515 + t572 * t516;
t520 = -t543 * qJD(3) - t568 * t550 + t572 * t551;
t536 = t564 * pkin(3) - t543 * pkin(7);
t538 = t542 ^ 2;
t492 = -t538 * pkin(3) + t520 * pkin(7) - t564 * t536 + t504;
t567 = sin(qJ(4));
t571 = cos(qJ(4));
t489 = t571 * t491 - t567 * t492;
t529 = t571 * t542 - t567 * t543;
t500 = t529 * qJD(4) + t567 * t520 + t571 * t521;
t530 = t567 * t542 + t571 * t543;
t510 = -t529 * mrSges(5,1) + t530 * mrSges(5,2);
t561 = qJD(4) + t564;
t523 = -t561 * mrSges(5,2) + t529 * mrSges(5,3);
t560 = qJDD(4) + t563;
t486 = m(5) * t489 + t560 * mrSges(5,1) - t500 * mrSges(5,3) - t530 * t510 + t561 * t523;
t490 = t567 * t491 + t571 * t492;
t499 = -t530 * qJD(4) + t571 * t520 - t567 * t521;
t524 = t561 * mrSges(5,1) - t530 * mrSges(5,3);
t487 = m(5) * t490 - t560 * mrSges(5,2) + t499 * mrSges(5,3) + t529 * t510 - t561 * t524;
t477 = t571 * t486 + t567 * t487;
t531 = -t542 * mrSges(4,1) + t543 * mrSges(4,2);
t534 = -t564 * mrSges(4,2) + t542 * mrSges(4,3);
t474 = m(4) * t503 + t563 * mrSges(4,1) - t521 * mrSges(4,3) - t543 * t531 + t564 * t534 + t477;
t535 = t564 * mrSges(4,1) - t543 * mrSges(4,3);
t585 = -t567 * t486 + t571 * t487;
t475 = m(4) * t504 - t563 * mrSges(4,2) + t520 * mrSges(4,3) + t542 * t531 - t564 * t535 + t585;
t470 = t572 * t474 + t568 * t475;
t532 = -t573 * g(3) - t593;
t540 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t569 + Ifges(3,2) * t573) * qJD(1);
t541 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t569 + Ifges(3,4) * t573) * qJD(1);
t526 = Ifges(4,4) * t543 + Ifges(4,2) * t542 + Ifges(4,6) * t564;
t527 = Ifges(4,1) * t543 + Ifges(4,4) * t542 + Ifges(4,5) * t564;
t506 = Ifges(5,4) * t530 + Ifges(5,2) * t529 + Ifges(5,6) * t561;
t507 = Ifges(5,1) * t530 + Ifges(5,4) * t529 + Ifges(5,5) * t561;
t580 = -mrSges(5,1) * t489 + mrSges(5,2) * t490 - Ifges(5,5) * t500 - Ifges(5,6) * t499 - Ifges(5,3) * t560 - t530 * t506 + t529 * t507;
t578 = -mrSges(4,1) * t503 + mrSges(4,2) * t504 - Ifges(4,5) * t521 - Ifges(4,6) * t520 - Ifges(4,3) * t563 - pkin(3) * t477 - t543 * t526 + t542 * t527 + t580;
t595 = mrSges(3,1) * t532 - mrSges(3,2) * t533 + Ifges(3,5) * t550 + Ifges(3,6) * t551 + Ifges(3,3) * qJDD(2) + pkin(2) * t470 + (t569 * t540 - t573 * t541) * qJD(1) - t578;
t549 = (-mrSges(3,1) * t573 + mrSges(3,2) * t569) * qJD(1);
t590 = qJD(1) * t573;
t554 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t590;
t468 = m(3) * t532 + qJDD(2) * mrSges(3,1) - t550 * mrSges(3,3) + qJD(2) * t554 - t549 * t591 + t470;
t553 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t591;
t586 = -t568 * t474 + t572 * t475;
t469 = m(3) * t533 - qJDD(2) * mrSges(3,2) + t551 * mrSges(3,3) - qJD(2) * t553 + t549 * t590 + t586;
t587 = -t569 * t468 + t573 * t469;
t460 = m(2) * t557 - t575 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t587;
t556 = t570 * g(1) - t574 * g(2);
t582 = -qJDD(1) * pkin(1) - t556;
t544 = -t575 * pkin(5) + t582;
t522 = -t551 * pkin(2) + t555 * t591 + (-pkin(6) * t566 - pkin(5)) * t575 + t582;
t494 = -t520 * pkin(3) - t538 * pkin(7) + t543 * t536 + t522;
t584 = m(5) * t494 - t499 * mrSges(5,1) + t500 * mrSges(5,2) - t529 * t523 + t530 * t524;
t579 = m(4) * t522 - t520 * mrSges(4,1) + t521 * mrSges(4,2) - t542 * t534 + t543 * t535 + t584;
t577 = -m(3) * t544 + t551 * mrSges(3,1) - t550 * mrSges(3,2) - t553 * t591 + t554 * t590 - t579;
t481 = m(2) * t556 + qJDD(1) * mrSges(2,1) - t575 * mrSges(2,2) + t577;
t592 = t570 * t460 + t574 * t481;
t462 = t573 * t468 + t569 * t469;
t588 = t574 * t460 - t570 * t481;
t505 = Ifges(5,5) * t530 + Ifges(5,6) * t529 + Ifges(5,3) * t561;
t478 = -mrSges(5,1) * t494 + mrSges(5,3) * t490 + Ifges(5,4) * t500 + Ifges(5,2) * t499 + Ifges(5,6) * t560 - t530 * t505 + t561 * t507;
t479 = mrSges(5,2) * t494 - mrSges(5,3) * t489 + Ifges(5,1) * t500 + Ifges(5,4) * t499 + Ifges(5,5) * t560 + t529 * t505 - t561 * t506;
t525 = Ifges(4,5) * t543 + Ifges(4,6) * t542 + Ifges(4,3) * t564;
t463 = -mrSges(4,1) * t522 + mrSges(4,3) * t504 + Ifges(4,4) * t521 + Ifges(4,2) * t520 + Ifges(4,6) * t563 - pkin(3) * t584 + pkin(7) * t585 + t571 * t478 + t567 * t479 - t543 * t525 + t564 * t527;
t464 = mrSges(4,2) * t522 - mrSges(4,3) * t503 + Ifges(4,1) * t521 + Ifges(4,4) * t520 + Ifges(4,5) * t563 - pkin(7) * t477 - t567 * t478 + t571 * t479 + t542 * t525 - t564 * t526;
t539 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t569 + Ifges(3,6) * t573) * qJD(1);
t455 = -mrSges(3,1) * t544 + mrSges(3,3) * t533 + Ifges(3,4) * t550 + Ifges(3,2) * t551 + Ifges(3,6) * qJDD(2) - pkin(2) * t579 + pkin(6) * t586 + qJD(2) * t541 + t572 * t463 + t568 * t464 - t539 * t591;
t457 = mrSges(3,2) * t544 - mrSges(3,3) * t532 + Ifges(3,1) * t550 + Ifges(3,4) * t551 + Ifges(3,5) * qJDD(2) - pkin(6) * t470 - qJD(2) * t540 - t568 * t463 + t572 * t464 + t539 * t590;
t581 = mrSges(2,1) * t556 - mrSges(2,2) * t557 + Ifges(2,3) * qJDD(1) + pkin(1) * t577 + pkin(5) * t587 + t573 * t455 + t569 * t457;
t453 = mrSges(2,1) * g(3) + mrSges(2,3) * t557 + t575 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t462 - t595;
t452 = -mrSges(2,2) * g(3) - mrSges(2,3) * t556 + Ifges(2,5) * qJDD(1) - t575 * Ifges(2,6) - pkin(5) * t462 - t569 * t455 + t573 * t457;
t1 = [-m(1) * g(1) + t588; -m(1) * g(2) + t592; (-m(1) - m(2)) * g(3) + t462; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t592 + t574 * t452 - t570 * t453; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t588 + t570 * t452 + t574 * t453; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t581; t581; t595; -t578; -t580;];
tauJB = t1;
