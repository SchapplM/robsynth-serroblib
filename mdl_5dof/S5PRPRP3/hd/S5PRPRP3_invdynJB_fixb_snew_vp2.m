% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRP3
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:54
% EndTime: 2019-12-05 15:32:56
% DurationCPUTime: 1.37s
% Computational Cost: add. (11535->196), mult. (20000->236), div. (0->0), fcn. (10526->8), ass. (0->89)
t579 = Ifges(5,1) + Ifges(6,1);
t570 = Ifges(5,4) + Ifges(6,4);
t569 = Ifges(5,5) + Ifges(6,5);
t578 = Ifges(5,2) + Ifges(6,2);
t568 = Ifges(5,6) + Ifges(6,6);
t577 = Ifges(5,3) + Ifges(6,3);
t541 = sin(qJ(4));
t543 = cos(qJ(4));
t517 = (-mrSges(6,1) * t543 + mrSges(6,2) * t541) * qJD(2);
t558 = qJD(2) * qJD(4);
t553 = t543 * t558;
t519 = t541 * qJDD(2) + t553;
t538 = sin(pkin(7));
t540 = cos(pkin(7));
t524 = -t540 * g(1) - t538 * g(2);
t536 = -g(3) + qJDD(1);
t542 = sin(qJ(2));
t544 = cos(qJ(2));
t500 = -t542 * t524 + t544 * t536;
t498 = qJDD(2) * pkin(2) + t500;
t501 = t544 * t524 + t542 * t536;
t545 = qJD(2) ^ 2;
t499 = -t545 * pkin(2) + t501;
t537 = sin(pkin(8));
t539 = cos(pkin(8));
t494 = t537 * t498 + t539 * t499;
t492 = -t545 * pkin(3) + qJDD(2) * pkin(6) + t494;
t523 = t538 * g(1) - t540 * g(2);
t522 = qJDD(3) - t523;
t510 = t543 * t522;
t557 = qJD(2) * qJD(5);
t572 = pkin(4) * t545;
t485 = qJDD(4) * pkin(4) + t510 + (-t519 + t553) * qJ(5) + (t543 * t572 - t492 - 0.2e1 * t557) * t541;
t559 = qJD(2) * t543;
t528 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t559;
t556 = m(6) * t485 + qJDD(4) * mrSges(6,1) + qJD(4) * t528;
t560 = qJD(2) * t541;
t481 = -t519 * mrSges(6,3) - t517 * t560 + t556;
t489 = t543 * t492 + t541 * t522;
t520 = t543 * qJDD(2) - t541 * t558;
t525 = qJD(4) * pkin(4) - qJ(5) * t560;
t535 = t543 ^ 2;
t486 = t520 * qJ(5) - qJD(4) * t525 - t535 * t572 + 0.2e1 * t543 * t557 + t489;
t488 = -t541 * t492 + t510;
t562 = t569 * qJD(4) + (t579 * t541 + t570 * t543) * qJD(2);
t563 = t568 * qJD(4) + (t570 * t541 + t578 * t543) * qJD(2);
t576 = mrSges(5,1) * t488 + mrSges(6,1) * t485 - mrSges(5,2) * t489 - mrSges(6,2) * t486 + pkin(4) * t481 + (t563 * t541 - t562 * t543) * qJD(2) + t577 * qJDD(4) + t569 * t519 + t568 * t520;
t518 = (-mrSges(5,1) * t543 + mrSges(5,2) * t541) * qJD(2);
t529 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t559;
t478 = m(5) * t488 + qJDD(4) * mrSges(5,1) + qJD(4) * t529 + (-mrSges(5,3) - mrSges(6,3)) * t519 + (-t517 - t518) * t560 + t556;
t555 = m(6) * t486 + t520 * mrSges(6,3) + t517 * t559;
t526 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t560;
t561 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t560 - t526;
t571 = -mrSges(5,2) - mrSges(6,2);
t479 = m(5) * t489 + t520 * mrSges(5,3) + t561 * qJD(4) + t571 * qJDD(4) + t518 * t559 + t555;
t472 = -t541 * t478 + t543 * t479;
t467 = m(4) * t494 - t545 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t472;
t493 = t539 * t498 - t537 * t499;
t548 = -qJDD(2) * pkin(3) - t493;
t491 = -t545 * pkin(6) + t548;
t487 = t525 * t560 - t520 * pkin(4) + qJDD(5) + (-qJ(5) * t535 - pkin(6)) * t545 + t548;
t549 = -m(6) * t487 + t520 * mrSges(6,1) + t528 * t559;
t480 = -m(5) * t491 + t520 * mrSges(5,1) + t571 * t519 + t529 * t559 + t561 * t560 + t549;
t474 = m(4) * t493 + qJDD(2) * mrSges(4,1) - t545 * mrSges(4,2) + t480;
t462 = t537 * t467 + t539 * t474;
t482 = t519 * mrSges(6,2) + t526 * t560 - t549;
t564 = t577 * qJD(4) + (t569 * t541 + t568 * t543) * qJD(2);
t463 = -mrSges(5,1) * t491 + mrSges(5,3) * t489 - mrSges(6,1) * t487 + mrSges(6,3) * t486 - pkin(4) * t482 + qJ(5) * t555 + t578 * t520 + t570 * t519 + (-qJ(5) * mrSges(6,2) + t568) * qJDD(4) + (-qJ(5) * t526 + t562) * qJD(4) - t564 * t560;
t464 = mrSges(5,2) * t491 + mrSges(6,2) * t487 - mrSges(5,3) * t488 - mrSges(6,3) * t485 - qJ(5) * t481 - t563 * qJD(4) + t569 * qJDD(4) + t579 * t519 + t570 * t520 + t564 * t559;
t575 = mrSges(3,1) * t500 + mrSges(4,1) * t493 - mrSges(3,2) * t501 - mrSges(4,2) * t494 + pkin(2) * t462 + pkin(3) * t480 + pkin(6) * t472 + t543 * t463 + t541 * t464 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t460 = m(3) * t500 + qJDD(2) * mrSges(3,1) - t545 * mrSges(3,2) + t462;
t550 = t539 * t467 - t537 * t474;
t461 = m(3) * t501 - t545 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t550;
t551 = -t542 * t460 + t544 * t461;
t454 = m(2) * t524 + t551;
t471 = t543 * t478 + t541 * t479;
t470 = m(4) * t522 + t471;
t469 = (m(2) + m(3)) * t523 - t470;
t565 = t538 * t454 + t540 * t469;
t455 = t544 * t460 + t542 * t461;
t554 = m(2) * t536 + t455;
t552 = t540 * t454 - t538 * t469;
t456 = -mrSges(4,1) * t522 + mrSges(4,3) * t494 + t545 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t471 - t576;
t451 = mrSges(4,2) * t522 - mrSges(4,3) * t493 + Ifges(4,5) * qJDD(2) - t545 * Ifges(4,6) - pkin(6) * t471 - t541 * t463 + t543 * t464;
t450 = -mrSges(3,2) * t523 - mrSges(3,3) * t500 + Ifges(3,5) * qJDD(2) - t545 * Ifges(3,6) - qJ(3) * t462 + t539 * t451 - t537 * t456;
t449 = mrSges(3,1) * t523 + mrSges(3,3) * t501 + t545 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t470 + qJ(3) * t550 + t537 * t451 + t539 * t456;
t448 = -mrSges(2,1) * t536 + mrSges(2,3) * t524 - pkin(1) * t455 - t575;
t447 = mrSges(2,2) * t536 - mrSges(2,3) * t523 - pkin(5) * t455 - t542 * t449 + t544 * t450;
t1 = [-m(1) * g(1) + t552; -m(1) * g(2) + t565; -m(1) * g(3) + t554; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t565 + t540 * t447 - t538 * t448; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t552 + t538 * t447 + t540 * t448; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t523 - mrSges(2,2) * t524 + t542 * t450 + t544 * t449 + pkin(1) * (m(3) * t523 - t470) + pkin(5) * t551; t554; t575; t470; t576; t482;];
tauJB = t1;
