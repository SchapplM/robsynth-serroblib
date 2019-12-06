% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:28
% EndTime: 2019-12-05 17:51:31
% DurationCPUTime: 2.75s
% Computational Cost: add. (36512->207), mult. (51081->271), div. (0->0), fcn. (28386->10), ass. (0->104)
t607 = 2 * qJD(4);
t571 = sin(qJ(1));
t574 = cos(qJ(1));
t548 = t574 * g(2) + t571 * g(3);
t541 = qJDD(1) * pkin(1) + t548;
t547 = t571 * g(2) - t574 * g(3);
t575 = qJD(1) ^ 2;
t542 = -t575 * pkin(1) + t547;
t566 = sin(pkin(8));
t568 = cos(pkin(8));
t527 = t568 * t541 - t566 * t542;
t521 = qJDD(1) * pkin(2) + t527;
t528 = t566 * t541 + t568 * t542;
t522 = -t575 * pkin(2) + t528;
t570 = sin(qJ(3));
t573 = cos(qJ(3));
t517 = t570 * t521 + t573 * t522;
t561 = (qJD(1) + qJD(3));
t559 = t561 ^ 2;
t560 = qJDD(1) + qJDD(3);
t514 = -t559 * pkin(3) + t560 * qJ(4) + t517;
t606 = (t561 * t607) + t514;
t565 = sin(pkin(9));
t605 = mrSges(5,2) * t565;
t603 = mrSges(5,3) * t560;
t602 = t565 * t561;
t569 = sin(qJ(5));
t601 = t565 * t569;
t572 = cos(qJ(5));
t600 = t565 * t572;
t567 = cos(pkin(9));
t599 = t567 * t560;
t598 = t567 * t561;
t564 = -g(1) + qJDD(2);
t597 = t567 * t564;
t510 = t565 * t564 + t606 * t567;
t535 = (-mrSges(5,1) * t567 + t605) * t561;
t586 = -pkin(4) * t567 - pkin(7) * t565;
t537 = t586 * t561;
t508 = t537 * t598 + t510;
t516 = t573 * t521 - t570 * t522;
t580 = -t559 * qJ(4) + qJDD(4) - t516;
t511 = (-pkin(3) + t586) * t560 + t580;
t505 = -t569 * t508 + t572 * t511;
t545 = qJD(5) - t598;
t594 = t561 * t601;
t530 = -t545 * mrSges(6,2) - mrSges(6,3) * t594;
t532 = (mrSges(6,1) * t569 + mrSges(6,2) * t572) * t602;
t595 = qJD(5) * t561;
t534 = (t560 * t572 - t569 * t595) * t565;
t544 = qJDD(5) - t599;
t593 = t561 * t600;
t503 = m(6) * t505 + t544 * mrSges(6,1) - t534 * mrSges(6,3) + t545 * t530 - t532 * t593;
t506 = t572 * t508 + t569 * t511;
t531 = t545 * mrSges(6,1) - mrSges(6,3) * t593;
t533 = (-t560 * t569 - t572 * t595) * t565;
t504 = m(6) * t506 - t544 * mrSges(6,2) + t533 * mrSges(6,3) - t545 * t531 - t532 * t594;
t587 = -t569 * t503 + t572 * t504;
t494 = m(5) * t510 + (t535 * t561 + t603) * t567 + t587;
t509 = -t606 * t565 + t597;
t507 = -t597 + (t514 + (t607 + t537) * t561) * t565;
t581 = -m(6) * t507 + t533 * mrSges(6,1) - t534 * mrSges(6,2);
t501 = m(5) * t509 + (-t603 + (-t530 * t569 - t531 * t572 - t535) * t561) * t565 + t581;
t588 = t567 * t494 - t565 * t501;
t486 = m(4) * t517 - t559 * mrSges(4,1) - t560 * mrSges(4,2) + t588;
t497 = t572 * t503 + t569 * t504;
t513 = -t560 * pkin(3) + t580;
t578 = -m(5) * t513 + mrSges(5,1) * t599 - t497 + (t565 ^ 2 + t567 ^ 2) * mrSges(5,3) * t559;
t491 = m(4) * t516 - t559 * mrSges(4,2) + (mrSges(4,1) - t605) * t560 + t578;
t479 = t570 * t486 + t573 * t491;
t476 = m(3) * t527 + qJDD(1) * mrSges(3,1) - t575 * mrSges(3,2) + t479;
t589 = t573 * t486 - t570 * t491;
t477 = m(3) * t528 - t575 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t589;
t471 = t568 * t476 + t566 * t477;
t489 = t565 * t494 + t567 * t501;
t592 = m(4) * t564 + t489;
t590 = -t566 * t476 + t568 * t477;
t468 = m(2) * t547 - t575 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t590;
t469 = m(2) * t548 + qJDD(1) * mrSges(2,1) - t575 * mrSges(2,2) + t471;
t591 = t574 * t468 - t571 * t469;
t487 = m(3) * t564 + t592;
t585 = Ifges(5,1) * t565 + Ifges(5,4) * t567;
t584 = Ifges(5,5) * t565 + Ifges(5,6) * t567;
t583 = -t571 * t468 - t574 * t469;
t525 = Ifges(6,6) * t545 + (Ifges(6,4) * t572 - Ifges(6,2) * t569) * t602;
t526 = Ifges(6,5) * t545 + (Ifges(6,1) * t572 - Ifges(6,4) * t569) * t602;
t582 = t525 * t572 + t526 * t569;
t524 = Ifges(6,3) * t545 + (Ifges(6,5) * t572 - Ifges(6,6) * t569) * t602;
t498 = -mrSges(6,1) * t507 + mrSges(6,3) * t506 + Ifges(6,4) * t534 + Ifges(6,2) * t533 + Ifges(6,6) * t544 - t524 * t593 + t545 * t526;
t499 = mrSges(6,2) * t507 - mrSges(6,3) * t505 + Ifges(6,1) * t534 + Ifges(6,4) * t533 + Ifges(6,5) * t544 - t524 * t594 - t545 * t525;
t536 = t584 * t561;
t481 = mrSges(5,2) * t513 - mrSges(5,3) * t509 - pkin(7) * t497 - t569 * t498 + t572 * t499 + t536 * t598 + t585 * t560;
t577 = mrSges(6,1) * t505 - mrSges(6,2) * t506 + Ifges(6,5) * t534 + Ifges(6,6) * t533 + Ifges(6,3) * t544;
t483 = Ifges(5,2) * t599 - mrSges(5,1) * t513 + mrSges(5,3) * t510 - pkin(4) * t497 + (Ifges(5,4) * t560 + (-t536 - t582) * t561) * t565 - t577;
t496 = t560 * t605 - t578;
t579 = mrSges(4,1) * t516 - mrSges(4,2) * t517 + Ifges(4,3) * t560 - pkin(3) * t496 + qJ(4) * t588 + t565 * t481 + t567 * t483;
t576 = mrSges(2,1) * t548 + mrSges(3,1) * t527 - mrSges(2,2) * t547 - mrSges(3,2) * t528 + pkin(1) * t471 + pkin(2) * t479 + t579 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t472 = t559 * Ifges(4,5) - mrSges(4,1) * t564 + mrSges(4,3) * t517 - mrSges(5,1) * t509 + mrSges(5,2) * t510 - t569 * t499 - t572 * t498 - pkin(4) * t581 - pkin(7) * t587 - pkin(3) * t489 + (Ifges(4,6) - t584) * t560 + (-pkin(4) * (-t530 * t601 - t531 * t600) + (-t565 * (Ifges(5,4) * t565 + Ifges(5,2) * t567) + t567 * t585) * t561) * t561;
t466 = mrSges(4,2) * t564 - mrSges(4,3) * t516 + Ifges(4,5) * t560 - t559 * Ifges(4,6) - qJ(4) * t489 + t567 * t481 - t565 * t483;
t465 = mrSges(3,2) * t564 - mrSges(3,3) * t527 + Ifges(3,5) * qJDD(1) - t575 * Ifges(3,6) - pkin(6) * t479 + t573 * t466 - t570 * t472;
t464 = -mrSges(3,1) * t564 + mrSges(3,3) * t528 + t575 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t592 + pkin(6) * t589 + t570 * t466 + t573 * t472;
t463 = -mrSges(2,2) * g(1) - mrSges(2,3) * t548 + Ifges(2,5) * qJDD(1) - t575 * Ifges(2,6) - qJ(2) * t471 - t566 * t464 + t568 * t465;
t462 = mrSges(2,1) * g(1) + mrSges(2,3) * t547 + t575 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t487 + qJ(2) * t590 + t568 * t464 + t566 * t465;
t1 = [(-m(1) - m(2)) * g(1) + t487; -m(1) * g(2) + t583; -m(1) * g(3) + t591; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t576; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t591 - t574 * t462 - t571 * t463; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t583 - t571 * t462 + t574 * t463; t576; t487; t579; t496; t582 * t602 + t577;];
tauJB = t1;
