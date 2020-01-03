% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR7
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:51
% EndTime: 2019-12-31 16:53:53
% DurationCPUTime: 1.93s
% Computational Cost: add. (17874->217), mult. (41982->275), div. (0->0), fcn. (28090->8), ass. (0->99)
t544 = qJD(1) ^ 2;
t535 = sin(pkin(7));
t536 = cos(pkin(7));
t538 = sin(qJ(3));
t541 = cos(qJ(3));
t551 = t535 * t538 - t536 * t541;
t518 = t551 * qJD(1);
t552 = t535 * t541 + t536 * t538;
t519 = t552 * qJD(1);
t563 = t519 * qJD(3);
t507 = -t551 * qJDD(1) - t563;
t569 = pkin(2) * t536;
t568 = mrSges(3,2) * t535;
t533 = t536 ^ 2;
t567 = t533 * t544;
t539 = sin(qJ(1));
t542 = cos(qJ(1));
t525 = -t542 * g(1) - t539 * g(2);
t520 = -t544 * pkin(1) + qJDD(1) * qJ(2) + t525;
t562 = qJD(1) * qJD(2);
t560 = -t536 * g(3) - 0.2e1 * t535 * t562;
t496 = (-pkin(5) * qJDD(1) + t544 * t569 - t520) * t535 + t560;
t510 = -t535 * g(3) + (t520 + 0.2e1 * t562) * t536;
t561 = qJDD(1) * t536;
t497 = -pkin(2) * t567 + pkin(5) * t561 + t510;
t483 = t538 * t496 + t541 * t497;
t505 = t518 * pkin(3) - t519 * pkin(6);
t543 = qJD(3) ^ 2;
t480 = -t543 * pkin(3) + qJDD(3) * pkin(6) - t518 * t505 + t483;
t532 = t535 ^ 2;
t524 = t539 * g(1) - t542 * g(2);
t556 = qJDD(2) - t524;
t506 = (-pkin(1) - t569) * qJDD(1) + (-qJ(2) + (-t532 - t533) * pkin(5)) * t544 + t556;
t564 = t518 * qJD(3);
t508 = t552 * qJDD(1) - t564;
t481 = (-t508 + t564) * pkin(6) + (-t507 + t563) * pkin(3) + t506;
t537 = sin(qJ(4));
t540 = cos(qJ(4));
t477 = -t537 * t480 + t540 * t481;
t511 = t540 * qJD(3) - t537 * t519;
t490 = t511 * qJD(4) + t537 * qJDD(3) + t540 * t508;
t512 = t537 * qJD(3) + t540 * t519;
t491 = -t511 * mrSges(5,1) + t512 * mrSges(5,2);
t516 = qJD(4) + t518;
t493 = -t516 * mrSges(5,2) + t511 * mrSges(5,3);
t504 = qJDD(4) - t507;
t474 = m(5) * t477 + t504 * mrSges(5,1) - t490 * mrSges(5,3) - t512 * t491 + t516 * t493;
t478 = t540 * t480 + t537 * t481;
t489 = -t512 * qJD(4) + t540 * qJDD(3) - t537 * t508;
t494 = t516 * mrSges(5,1) - t512 * mrSges(5,3);
t475 = m(5) * t478 - t504 * mrSges(5,2) + t489 * mrSges(5,3) + t511 * t491 - t516 * t494;
t466 = -t537 * t474 + t540 * t475;
t502 = t518 * mrSges(4,1) + t519 * mrSges(4,2);
t514 = qJD(3) * mrSges(4,1) - t519 * mrSges(4,3);
t464 = m(4) * t483 - qJDD(3) * mrSges(4,2) + t507 * mrSges(4,3) - qJD(3) * t514 - t518 * t502 + t466;
t482 = t541 * t496 - t538 * t497;
t479 = -qJDD(3) * pkin(3) - t543 * pkin(6) + t519 * t505 - t482;
t476 = -m(5) * t479 + t489 * mrSges(5,1) - t490 * mrSges(5,2) + t511 * t493 - t512 * t494;
t513 = -qJD(3) * mrSges(4,2) - t518 * mrSges(4,3);
t470 = m(4) * t482 + qJDD(3) * mrSges(4,1) - t508 * mrSges(4,3) + qJD(3) * t513 - t519 * t502 + t476;
t457 = t538 * t464 + t541 * t470;
t509 = -t535 * t520 + t560;
t550 = mrSges(3,3) * qJDD(1) + t544 * (-mrSges(3,1) * t536 + t568);
t455 = m(3) * t509 - t550 * t535 + t457;
t557 = t541 * t464 - t538 * t470;
t456 = m(3) * t510 + t550 * t536 + t557;
t558 = -t535 * t455 + t536 * t456;
t447 = m(2) * t525 - t544 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t558;
t517 = -qJDD(1) * pkin(1) - t544 * qJ(2) + t556;
t465 = t540 * t474 + t537 * t475;
t548 = m(4) * t506 - t507 * mrSges(4,1) + t508 * mrSges(4,2) + t518 * t513 + t519 * t514 + t465;
t547 = -m(3) * t517 + mrSges(3,1) * t561 - t548 + (t532 * t544 + t567) * mrSges(3,3);
t459 = m(2) * t524 - t544 * mrSges(2,2) + t547 + (mrSges(2,1) - t568) * qJDD(1);
t566 = t539 * t447 + t542 * t459;
t449 = t536 * t455 + t535 * t456;
t553 = Ifges(3,5) * t535 + Ifges(3,6) * t536;
t565 = t544 * t553;
t559 = t542 * t447 - t539 * t459;
t555 = Ifges(3,1) * t535 + Ifges(3,4) * t536;
t554 = Ifges(3,4) * t535 + Ifges(3,2) * t536;
t484 = Ifges(5,5) * t512 + Ifges(5,6) * t511 + Ifges(5,3) * t516;
t486 = Ifges(5,1) * t512 + Ifges(5,4) * t511 + Ifges(5,5) * t516;
t467 = -mrSges(5,1) * t479 + mrSges(5,3) * t478 + Ifges(5,4) * t490 + Ifges(5,2) * t489 + Ifges(5,6) * t504 - t512 * t484 + t516 * t486;
t485 = Ifges(5,4) * t512 + Ifges(5,2) * t511 + Ifges(5,6) * t516;
t468 = mrSges(5,2) * t479 - mrSges(5,3) * t477 + Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * t504 + t511 * t484 - t516 * t485;
t498 = Ifges(4,5) * t519 - Ifges(4,6) * t518 + Ifges(4,3) * qJD(3);
t499 = Ifges(4,4) * t519 - Ifges(4,2) * t518 + Ifges(4,6) * qJD(3);
t450 = mrSges(4,2) * t506 - mrSges(4,3) * t482 + Ifges(4,1) * t508 + Ifges(4,4) * t507 + Ifges(4,5) * qJDD(3) - pkin(6) * t465 - qJD(3) * t499 - t537 * t467 + t540 * t468 - t518 * t498;
t500 = Ifges(4,1) * t519 - Ifges(4,4) * t518 + Ifges(4,5) * qJD(3);
t546 = mrSges(5,1) * t477 - mrSges(5,2) * t478 + Ifges(5,5) * t490 + Ifges(5,6) * t489 + Ifges(5,3) * t504 + t512 * t485 - t511 * t486;
t451 = -mrSges(4,1) * t506 + mrSges(4,3) * t483 + Ifges(4,4) * t508 + Ifges(4,2) * t507 + Ifges(4,6) * qJDD(3) - pkin(3) * t465 + qJD(3) * t500 - t519 * t498 - t546;
t442 = -mrSges(3,1) * t517 + mrSges(3,3) * t510 - pkin(2) * t548 + pkin(5) * t557 + t554 * qJDD(1) + t538 * t450 + t541 * t451 - t535 * t565;
t444 = mrSges(3,2) * t517 - mrSges(3,3) * t509 - pkin(5) * t457 + t555 * qJDD(1) + t541 * t450 - t538 * t451 + t536 * t565;
t461 = qJDD(1) * t568 - t547;
t549 = mrSges(2,1) * t524 - mrSges(2,2) * t525 + Ifges(2,3) * qJDD(1) - pkin(1) * t461 + qJ(2) * t558 + t536 * t442 + t535 * t444;
t545 = mrSges(4,1) * t482 - mrSges(4,2) * t483 + Ifges(4,5) * t508 + Ifges(4,6) * t507 + Ifges(4,3) * qJDD(3) + pkin(3) * t476 + pkin(6) * t466 + t540 * t467 + t537 * t468 + t519 * t499 + t518 * t500;
t440 = mrSges(2,3) * t525 + (Ifges(2,6) - t553) * qJDD(1) - mrSges(3,1) * t509 + mrSges(3,2) * t510 + mrSges(2,1) * g(3) - pkin(1) * t449 - pkin(2) * t457 - t545 + (-t535 * t554 + t536 * t555 + Ifges(2,5)) * t544;
t439 = -mrSges(2,2) * g(3) - mrSges(2,3) * t524 + Ifges(2,5) * qJDD(1) - t544 * Ifges(2,6) - qJ(2) * t449 - t535 * t442 + t536 * t444;
t1 = [-m(1) * g(1) + t559; -m(1) * g(2) + t566; (-m(1) - m(2)) * g(3) + t449; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t566 + t542 * t439 - t539 * t440; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t559 + t539 * t439 + t542 * t440; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t549; t549; t461; t545; t546;];
tauJB = t1;
