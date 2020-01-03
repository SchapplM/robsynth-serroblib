% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR7
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:08
% EndTime: 2019-12-31 17:06:10
% DurationCPUTime: 2.15s
% Computational Cost: add. (20531->240), mult. (46257->310), div. (0->0), fcn. (29091->8), ass. (0->98)
t536 = sin(pkin(7));
t537 = cos(pkin(7));
t539 = sin(qJ(2));
t542 = cos(qJ(2));
t514 = (t536 * t539 - t537 * t542) * qJD(1);
t557 = qJD(1) * qJD(2);
t525 = t539 * qJDD(1) + t542 * t557;
t540 = sin(qJ(1));
t543 = cos(qJ(1));
t532 = -t543 * g(1) - t540 * g(2);
t545 = qJD(1) ^ 2;
t520 = -t545 * pkin(1) + qJDD(1) * pkin(5) + t532;
t561 = t539 * t520;
t563 = pkin(2) * t545;
t489 = qJDD(2) * pkin(2) - t525 * qJ(3) - t561 + (qJ(3) * t557 + t539 * t563 - g(3)) * t542;
t508 = -t539 * g(3) + t542 * t520;
t526 = t542 * qJDD(1) - t539 * t557;
t559 = qJD(1) * t539;
t528 = qJD(2) * pkin(2) - qJ(3) * t559;
t535 = t542 ^ 2;
t490 = t526 * qJ(3) - qJD(2) * t528 - t535 * t563 + t508;
t564 = 2 * qJD(3);
t479 = t536 * t489 + t537 * t490 - t514 * t564;
t515 = (t536 * t542 + t537 * t539) * qJD(1);
t500 = t514 * pkin(3) - t515 * pkin(6);
t544 = qJD(2) ^ 2;
t476 = -t544 * pkin(3) + qJDD(2) * pkin(6) - t514 * t500 + t479;
t531 = t540 * g(1) - t543 * g(2);
t550 = -qJDD(1) * pkin(1) - t531;
t492 = -t526 * pkin(2) + qJDD(3) + t528 * t559 + (-qJ(3) * t535 - pkin(5)) * t545 + t550;
t503 = -t536 * t525 + t537 * t526;
t504 = t537 * t525 + t536 * t526;
t477 = (qJD(2) * t514 - t504) * pkin(6) + (qJD(2) * t515 - t503) * pkin(3) + t492;
t538 = sin(qJ(4));
t541 = cos(qJ(4));
t473 = -t538 * t476 + t541 * t477;
t505 = t541 * qJD(2) - t538 * t515;
t486 = t505 * qJD(4) + t538 * qJDD(2) + t541 * t504;
t506 = t538 * qJD(2) + t541 * t515;
t491 = -t505 * mrSges(5,1) + t506 * mrSges(5,2);
t513 = qJD(4) + t514;
t493 = -t513 * mrSges(5,2) + t505 * mrSges(5,3);
t502 = qJDD(4) - t503;
t470 = m(5) * t473 + t502 * mrSges(5,1) - t486 * mrSges(5,3) - t506 * t491 + t513 * t493;
t474 = t541 * t476 + t538 * t477;
t485 = -t506 * qJD(4) + t541 * qJDD(2) - t538 * t504;
t494 = t513 * mrSges(5,1) - t506 * mrSges(5,3);
t471 = m(5) * t474 - t502 * mrSges(5,2) + t485 * mrSges(5,3) + t505 * t491 - t513 * t494;
t462 = -t538 * t470 + t541 * t471;
t499 = t514 * mrSges(4,1) + t515 * mrSges(4,2);
t510 = qJD(2) * mrSges(4,1) - t515 * mrSges(4,3);
t459 = m(4) * t479 - qJDD(2) * mrSges(4,2) + t503 * mrSges(4,3) - qJD(2) * t510 - t514 * t499 + t462;
t553 = -t537 * t489 + t536 * t490;
t475 = -qJDD(2) * pkin(3) - t544 * pkin(6) + (t564 + t500) * t515 + t553;
t472 = -m(5) * t475 + t485 * mrSges(5,1) - t486 * mrSges(5,2) + t505 * t493 - t506 * t494;
t478 = -0.2e1 * qJD(3) * t515 - t553;
t509 = -qJD(2) * mrSges(4,2) - t514 * mrSges(4,3);
t466 = m(4) * t478 + qJDD(2) * mrSges(4,1) - t504 * mrSges(4,3) + qJD(2) * t509 - t515 * t499 + t472;
t453 = t536 * t459 + t537 * t466;
t480 = Ifges(5,5) * t506 + Ifges(5,6) * t505 + Ifges(5,3) * t513;
t482 = Ifges(5,1) * t506 + Ifges(5,4) * t505 + Ifges(5,5) * t513;
t463 = -mrSges(5,1) * t475 + mrSges(5,3) * t474 + Ifges(5,4) * t486 + Ifges(5,2) * t485 + Ifges(5,6) * t502 - t506 * t480 + t513 * t482;
t481 = Ifges(5,4) * t506 + Ifges(5,2) * t505 + Ifges(5,6) * t513;
t464 = mrSges(5,2) * t475 - mrSges(5,3) * t473 + Ifges(5,1) * t486 + Ifges(5,4) * t485 + Ifges(5,5) * t502 + t505 * t480 - t513 * t481;
t496 = Ifges(4,4) * t515 - Ifges(4,2) * t514 + Ifges(4,6) * qJD(2);
t497 = Ifges(4,1) * t515 - Ifges(4,4) * t514 + Ifges(4,5) * qJD(2);
t507 = -t542 * g(3) - t561;
t517 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t539 + Ifges(3,2) * t542) * qJD(1);
t518 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t539 + Ifges(3,4) * t542) * qJD(1);
t565 = (t539 * t517 - t542 * t518) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(3,1) * t507 + mrSges(4,1) * t478 - mrSges(3,2) * t508 - mrSges(4,2) * t479 + Ifges(3,5) * t525 + Ifges(4,5) * t504 + Ifges(3,6) * t526 + Ifges(4,6) * t503 + pkin(2) * t453 + pkin(3) * t472 + pkin(6) * t462 + t541 * t463 + t538 * t464 + t515 * t496 + t514 * t497;
t524 = (-mrSges(3,1) * t542 + mrSges(3,2) * t539) * qJD(1);
t558 = qJD(1) * t542;
t530 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t558;
t451 = m(3) * t507 + qJDD(2) * mrSges(3,1) - t525 * mrSges(3,3) + qJD(2) * t530 - t524 * t559 + t453;
t529 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t559;
t554 = t537 * t459 - t536 * t466;
t452 = m(3) * t508 - qJDD(2) * mrSges(3,2) + t526 * mrSges(3,3) - qJD(2) * t529 + t524 * t558 + t554;
t555 = -t539 * t451 + t542 * t452;
t443 = m(2) * t532 - t545 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t555;
t461 = t541 * t470 + t538 * t471;
t460 = m(4) * t492 - t503 * mrSges(4,1) + t504 * mrSges(4,2) + t514 * t509 + t515 * t510 + t461;
t519 = -t545 * pkin(5) + t550;
t547 = -m(3) * t519 + t526 * mrSges(3,1) - t525 * mrSges(3,2) - t529 * t559 + t530 * t558 - t460;
t455 = m(2) * t531 + qJDD(1) * mrSges(2,1) - t545 * mrSges(2,2) + t547;
t560 = t540 * t443 + t543 * t455;
t445 = t542 * t451 + t539 * t452;
t556 = t543 * t443 - t540 * t455;
t495 = Ifges(4,5) * t515 - Ifges(4,6) * t514 + Ifges(4,3) * qJD(2);
t446 = mrSges(4,2) * t492 - mrSges(4,3) * t478 + Ifges(4,1) * t504 + Ifges(4,4) * t503 + Ifges(4,5) * qJDD(2) - pkin(6) * t461 - qJD(2) * t496 - t538 * t463 + t541 * t464 - t514 * t495;
t548 = mrSges(5,1) * t473 - mrSges(5,2) * t474 + Ifges(5,5) * t486 + Ifges(5,6) * t485 + Ifges(5,3) * t502 + t506 * t481 - t505 * t482;
t447 = -mrSges(4,1) * t492 + mrSges(4,3) * t479 + Ifges(4,4) * t504 + Ifges(4,2) * t503 + Ifges(4,6) * qJDD(2) - pkin(3) * t461 + qJD(2) * t497 - t515 * t495 - t548;
t516 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t539 + Ifges(3,6) * t542) * qJD(1);
t438 = -mrSges(3,1) * t519 + mrSges(3,3) * t508 + Ifges(3,4) * t525 + Ifges(3,2) * t526 + Ifges(3,6) * qJDD(2) - pkin(2) * t460 + qJ(3) * t554 + qJD(2) * t518 + t536 * t446 + t537 * t447 - t516 * t559;
t440 = mrSges(3,2) * t519 - mrSges(3,3) * t507 + Ifges(3,1) * t525 + Ifges(3,4) * t526 + Ifges(3,5) * qJDD(2) - qJ(3) * t453 - qJD(2) * t517 + t537 * t446 - t536 * t447 + t516 * t558;
t549 = mrSges(2,1) * t531 - mrSges(2,2) * t532 + Ifges(2,3) * qJDD(1) + pkin(1) * t547 + pkin(5) * t555 + t542 * t438 + t539 * t440;
t436 = mrSges(2,1) * g(3) + mrSges(2,3) * t532 + t545 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t445 - t565;
t435 = -mrSges(2,2) * g(3) - mrSges(2,3) * t531 + Ifges(2,5) * qJDD(1) - t545 * Ifges(2,6) - pkin(5) * t445 - t539 * t438 + t542 * t440;
t1 = [-m(1) * g(1) + t556; -m(1) * g(2) + t560; (-m(1) - m(2)) * g(3) + t445; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t560 + t543 * t435 - t540 * t436; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t556 + t540 * t435 + t543 * t436; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t549; t549; t565; t460; t548;];
tauJB = t1;
