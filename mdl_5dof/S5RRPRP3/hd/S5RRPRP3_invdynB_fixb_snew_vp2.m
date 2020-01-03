% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:05
% EndTime: 2019-12-31 19:51:07
% DurationCPUTime: 2.25s
% Computational Cost: add. (28980->226), mult. (39979->270), div. (0->0), fcn. (24046->8), ass. (0->100)
t549 = Ifges(5,1) + Ifges(6,1);
t542 = Ifges(5,4) - Ifges(6,5);
t541 = Ifges(5,5) + Ifges(6,4);
t548 = Ifges(5,2) + Ifges(6,3);
t540 = Ifges(5,6) - Ifges(6,6);
t547 = -Ifges(5,3) - Ifges(6,2);
t499 = qJD(1) + qJD(2);
t495 = t499 ^ 2;
t502 = sin(pkin(8));
t503 = cos(pkin(8));
t504 = sin(qJ(4));
t545 = cos(qJ(4));
t546 = t502 * t504 - t503 * t545;
t544 = pkin(3) * t503;
t543 = -mrSges(5,3) - mrSges(6,2);
t539 = mrSges(4,2) * t502;
t516 = Ifges(4,5) * t502 + Ifges(4,6) * t503;
t538 = t516 * t495;
t498 = t503 ^ 2;
t537 = t495 * t498;
t496 = qJDD(1) + qJDD(2);
t536 = t496 * t503;
t506 = sin(qJ(1));
t508 = cos(qJ(1));
t489 = t506 * g(1) - t508 * g(2);
t483 = qJDD(1) * pkin(1) + t489;
t490 = -t508 * g(1) - t506 * g(2);
t510 = qJD(1) ^ 2;
t484 = -t510 * pkin(1) + t490;
t505 = sin(qJ(2));
t507 = cos(qJ(2));
t469 = t505 * t483 + t507 * t484;
t467 = -t495 * pkin(2) + t496 * qJ(3) + t469;
t529 = qJD(3) * t499;
t524 = -t503 * g(3) - 0.2e1 * t502 * t529;
t443 = (-pkin(7) * t496 + t495 * t544 - t467) * t502 + t524;
t447 = -t502 * g(3) + (t467 + 0.2e1 * t529) * t503;
t444 = -pkin(3) * t537 + pkin(7) * t536 + t447;
t440 = t504 * t443 + t545 * t444;
t513 = t545 * t502 + t503 * t504;
t477 = t513 * t499;
t527 = t477 * qJD(4);
t465 = t546 * t496 + t527;
t473 = qJD(4) * mrSges(5,1) - t477 * mrSges(5,3);
t476 = t546 * t499;
t461 = t476 * pkin(4) - t477 * qJ(5);
t509 = qJD(4) ^ 2;
t435 = -t509 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t476 * t461 + t440;
t474 = -qJD(4) * mrSges(6,1) + t477 * mrSges(6,2);
t526 = m(6) * t435 + qJDD(4) * mrSges(6,3) + qJD(4) * t474;
t462 = t476 * mrSges(6,1) - t477 * mrSges(6,3);
t530 = -t476 * mrSges(5,1) - t477 * mrSges(5,2) - t462;
t431 = m(5) * t440 - qJDD(4) * mrSges(5,2) - qJD(4) * t473 + t543 * t465 + t530 * t476 + t526;
t439 = t545 * t443 - t504 * t444;
t528 = t476 * qJD(4);
t466 = t513 * t496 - t528;
t472 = -qJD(4) * mrSges(5,2) - t476 * mrSges(5,3);
t436 = -qJDD(4) * pkin(4) - t509 * qJ(5) + t477 * t461 + qJDD(5) - t439;
t475 = -t476 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t519 = -m(6) * t436 + qJDD(4) * mrSges(6,1) + qJD(4) * t475;
t432 = m(5) * t439 + qJDD(4) * mrSges(5,1) + qJD(4) * t472 + t543 * t466 + t530 * t477 + t519;
t425 = t504 * t431 + t545 * t432;
t446 = -t502 * t467 + t524;
t514 = mrSges(4,3) * t496 + (-mrSges(4,1) * t503 + t539) * t495;
t423 = m(4) * t446 - t514 * t502 + t425;
t520 = t545 * t431 - t504 * t432;
t424 = m(4) * t447 + t514 * t503 + t520;
t521 = -t502 * t423 + t503 * t424;
t416 = m(3) * t469 - t495 * mrSges(3,1) - t496 * mrSges(3,2) + t521;
t468 = t507 * t483 - t505 * t484;
t515 = qJDD(3) - t468;
t464 = -t496 * pkin(2) - t495 * qJ(3) + t515;
t497 = t502 ^ 2;
t445 = (-pkin(2) - t544) * t496 + (-qJ(3) + (-t497 - t498) * pkin(7)) * t495 + t515;
t438 = -0.2e1 * qJD(5) * t477 + (-t466 + t528) * qJ(5) + (t465 + t527) * pkin(4) + t445;
t433 = m(6) * t438 + t465 * mrSges(6,1) - t466 * mrSges(6,3) - t477 * t474 + t476 * t475;
t512 = m(5) * t445 + t465 * mrSges(5,1) + t466 * mrSges(5,2) + t476 * t472 + t477 * t473 + t433;
t511 = -m(4) * t464 + mrSges(4,1) * t536 - t512 + (t495 * t497 + t537) * mrSges(4,3);
t427 = t511 - t495 * mrSges(3,2) + m(3) * t468 + (mrSges(3,1) - t539) * t496;
t413 = t505 * t416 + t507 * t427;
t411 = m(2) * t489 + qJDD(1) * mrSges(2,1) - t510 * mrSges(2,2) + t413;
t522 = t507 * t416 - t505 * t427;
t412 = m(2) * t490 - t510 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t522;
t534 = t508 * t411 + t506 * t412;
t417 = t503 * t423 + t502 * t424;
t533 = -t540 * qJD(4) + t548 * t476 - t542 * t477;
t532 = t547 * qJD(4) + t540 * t476 - t541 * t477;
t531 = t541 * qJD(4) - t542 * t476 + t549 * t477;
t523 = -t506 * t411 + t508 * t412;
t518 = Ifges(4,1) * t502 + Ifges(4,4) * t503;
t517 = Ifges(4,4) * t502 + Ifges(4,2) * t503;
t419 = mrSges(5,2) * t445 + mrSges(6,2) * t436 - mrSges(5,3) * t439 - mrSges(6,3) * t438 - qJ(5) * t433 + t533 * qJD(4) + t541 * qJDD(4) - t542 * t465 + t549 * t466 + t532 * t476;
t418 = -mrSges(5,1) * t445 - mrSges(6,1) * t438 + mrSges(6,2) * t435 + mrSges(5,3) * t440 - pkin(4) * t433 + t531 * qJD(4) + t540 * qJDD(4) - t548 * t465 + t542 * t466 + t532 * t477;
t407 = mrSges(4,2) * t464 - mrSges(4,3) * t446 - pkin(7) * t425 - t504 * t418 + t545 * t419 + t518 * t496 + t503 * t538;
t406 = -mrSges(4,1) * t464 + mrSges(4,3) * t447 - pkin(3) * t512 + pkin(7) * t520 + t545 * t418 + t504 * t419 + t517 * t496 - t502 * t538;
t405 = mrSges(3,1) * g(3) - qJ(5) * t526 - pkin(4) * t519 + mrSges(3,3) * t469 + mrSges(5,2) * t440 - mrSges(4,1) * t446 + mrSges(4,2) * t447 - mrSges(6,3) * t435 + mrSges(6,1) * t436 - mrSges(5,1) * t439 - pkin(3) * t425 - pkin(2) * t417 + (pkin(4) * t462 + t533) * t477 + (qJ(5) * t462 - t531) * t476 + (pkin(4) * mrSges(6,2) - t541) * t466 + (qJ(5) * mrSges(6,2) + t540) * t465 + t547 * qJDD(4) + (Ifges(3,6) - t516) * t496 + (-t502 * t517 + t503 * t518 + Ifges(3,5)) * t495;
t404 = -mrSges(3,2) * g(3) - mrSges(3,3) * t468 + Ifges(3,5) * t496 - t495 * Ifges(3,6) - qJ(3) * t417 - t502 * t406 + t503 * t407;
t403 = -mrSges(2,2) * g(3) - mrSges(2,3) * t489 + Ifges(2,5) * qJDD(1) - t510 * Ifges(2,6) - pkin(6) * t413 + t507 * t404 - t505 * t405;
t402 = Ifges(2,6) * qJDD(1) + t510 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t490 + t505 * t404 + t507 * t405 - pkin(1) * (-m(3) * g(3) + t417) + pkin(6) * t522;
t1 = [-m(1) * g(1) + t523; -m(1) * g(2) + t534; (-m(1) - m(2) - m(3)) * g(3) + t417; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t534 - t506 * t402 + t508 * t403; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t523 + t508 * t402 + t506 * t403; pkin(1) * t413 + mrSges(2,1) * t489 - mrSges(2,2) * t490 + qJ(3) * t521 + t502 * t407 + t503 * t406 + pkin(2) * (-t496 * t539 + t511) + mrSges(3,1) * t468 - mrSges(3,2) * t469 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * t496 + Ifges(2,3) * qJDD(1);];
tauB = t1;
