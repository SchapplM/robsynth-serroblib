% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:13
% EndTime: 2019-12-31 17:49:15
% DurationCPUTime: 1.98s
% Computational Cost: add. (18564->225), mult. (39979->271), div. (0->0), fcn. (24046->8), ass. (0->100)
t543 = Ifges(5,1) + Ifges(6,1);
t537 = Ifges(5,4) - Ifges(6,5);
t536 = Ifges(5,5) + Ifges(6,4);
t542 = Ifges(5,2) + Ifges(6,3);
t535 = Ifges(5,6) - Ifges(6,6);
t541 = -Ifges(5,3) - Ifges(6,2);
t505 = qJD(1) ^ 2;
t540 = cos(qJ(4));
t499 = cos(pkin(8));
t539 = pkin(3) * t499;
t538 = -mrSges(5,3) - mrSges(6,2);
t497 = sin(pkin(8));
t534 = mrSges(4,2) * t497;
t493 = t499 ^ 2;
t533 = t493 * t505;
t502 = sin(qJ(1));
t503 = cos(qJ(1));
t479 = t502 * g(1) - t503 * g(2);
t477 = qJDD(1) * pkin(1) + t479;
t480 = -t503 * g(1) - t502 * g(2);
t478 = -t505 * pkin(1) + t480;
t498 = sin(pkin(7));
t500 = cos(pkin(7));
t463 = t498 * t477 + t500 * t478;
t450 = -t505 * pkin(2) + qJDD(1) * qJ(3) + t463;
t496 = -g(3) + qJDD(2);
t523 = qJD(1) * qJD(3);
t527 = t499 * t496 - 0.2e1 * t497 * t523;
t437 = (-pkin(6) * qJDD(1) + t505 * t539 - t450) * t497 + t527;
t441 = t497 * t496 + (t450 + 0.2e1 * t523) * t499;
t522 = qJDD(1) * t499;
t438 = -pkin(3) * t533 + pkin(6) * t522 + t441;
t501 = sin(qJ(4));
t434 = t501 * t437 + t540 * t438;
t519 = t499 * t540;
t508 = t540 * t497 + t499 * t501;
t471 = t508 * qJD(1);
t525 = t471 * qJD(4);
t460 = t525 + (t497 * t501 - t519) * qJDD(1);
t467 = qJD(4) * mrSges(5,1) - t471 * mrSges(5,3);
t524 = t497 * qJD(1);
t470 = -qJD(1) * t519 + t501 * t524;
t454 = t470 * pkin(4) - t471 * qJ(5);
t504 = qJD(4) ^ 2;
t429 = -t504 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t470 * t454 + t434;
t468 = -qJD(4) * mrSges(6,1) + t471 * mrSges(6,2);
t521 = m(6) * t429 + qJDD(4) * mrSges(6,3) + qJD(4) * t468;
t455 = t470 * mrSges(6,1) - t471 * mrSges(6,3);
t528 = -t470 * mrSges(5,1) - t471 * mrSges(5,2) - t455;
t425 = m(5) * t434 - qJDD(4) * mrSges(5,2) - qJD(4) * t467 + t538 * t460 + t528 * t470 + t521;
t433 = t540 * t437 - t501 * t438;
t526 = t470 * qJD(4);
t461 = t508 * qJDD(1) - t526;
t466 = -qJD(4) * mrSges(5,2) - t470 * mrSges(5,3);
t430 = -qJDD(4) * pkin(4) - t504 * qJ(5) + t471 * t454 + qJDD(5) - t433;
t469 = -t470 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t514 = -m(6) * t430 + qJDD(4) * mrSges(6,1) + qJD(4) * t469;
t426 = m(5) * t433 + qJDD(4) * mrSges(5,1) + qJD(4) * t466 + t538 * t461 + t528 * t471 + t514;
t419 = t501 * t425 + t540 * t426;
t440 = -t497 * t450 + t527;
t509 = mrSges(4,3) * qJDD(1) + t505 * (-mrSges(4,1) * t499 + t534);
t417 = m(4) * t440 - t509 * t497 + t419;
t515 = t540 * t425 - t501 * t426;
t418 = m(4) * t441 + t509 * t499 + t515;
t516 = -t497 * t417 + t499 * t418;
t410 = m(3) * t463 - t505 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t516;
t462 = t500 * t477 - t498 * t478;
t510 = qJDD(3) - t462;
t443 = -qJDD(1) * pkin(2) - t505 * qJ(3) + t510;
t492 = t497 ^ 2;
t439 = (-pkin(2) - t539) * qJDD(1) + (-qJ(3) + (-t492 - t493) * pkin(6)) * t505 + t510;
t432 = -0.2e1 * qJD(5) * t471 + (-t461 + t526) * qJ(5) + (t460 + t525) * pkin(4) + t439;
t427 = m(6) * t432 + t460 * mrSges(6,1) - t461 * mrSges(6,3) - t471 * t468 + t470 * t469;
t507 = m(5) * t439 + t460 * mrSges(5,1) + t461 * mrSges(5,2) + t470 * t466 + t471 * t467 + t427;
t506 = -m(4) * t443 + mrSges(4,1) * t522 - t507 + (t492 * t505 + t533) * mrSges(4,3);
t421 = (mrSges(3,1) - t534) * qJDD(1) + t506 - t505 * mrSges(3,2) + m(3) * t462;
t407 = t498 * t410 + t500 * t421;
t405 = m(2) * t479 + qJDD(1) * mrSges(2,1) - t505 * mrSges(2,2) + t407;
t517 = t500 * t410 - t498 * t421;
t406 = m(2) * t480 - t505 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t517;
t532 = t503 * t405 + t502 * t406;
t411 = t499 * t417 + t497 * t418;
t531 = -t535 * qJD(4) + t542 * t470 - t537 * t471;
t530 = t541 * qJD(4) + t535 * t470 - t536 * t471;
t529 = t536 * qJD(4) - t537 * t470 + t543 * t471;
t520 = m(3) * t496 + t411;
t518 = -t502 * t405 + t503 * t406;
t513 = Ifges(4,1) * t497 + Ifges(4,4) * t499;
t512 = Ifges(4,4) * t497 + Ifges(4,2) * t499;
t511 = Ifges(4,5) * t497 + Ifges(4,6) * t499;
t476 = t511 * qJD(1);
t413 = mrSges(5,2) * t439 + mrSges(6,2) * t430 - mrSges(5,3) * t433 - mrSges(6,3) * t432 - qJ(5) * t427 + t531 * qJD(4) + t536 * qJDD(4) - t537 * t460 + t543 * t461 + t530 * t470;
t412 = -mrSges(5,1) * t439 - mrSges(6,1) * t432 + mrSges(6,2) * t429 + mrSges(5,3) * t434 - pkin(4) * t427 + t529 * qJD(4) + t535 * qJDD(4) - t542 * t460 + t537 * t461 + t530 * t471;
t401 = t499 * qJD(1) * t476 + mrSges(4,2) * t443 - mrSges(4,3) * t440 - pkin(6) * t419 + t513 * qJDD(1) - t501 * t412 + t540 * t413;
t400 = -mrSges(4,1) * t443 + mrSges(4,3) * t441 - pkin(3) * t507 + pkin(6) * t515 + t512 * qJDD(1) + t540 * t412 + t501 * t413 - t476 * t524;
t399 = -qJ(5) * t521 - pkin(4) * t514 - mrSges(3,1) * t496 + mrSges(3,3) * t463 - mrSges(4,1) * t440 + mrSges(4,2) * t441 + mrSges(6,1) * t430 - mrSges(5,1) * t433 + mrSges(5,2) * t434 - mrSges(6,3) * t429 - pkin(3) * t419 - pkin(2) * t411 + (pkin(4) * t455 + t531) * t471 + (qJ(5) * t455 - t529) * t470 + (pkin(4) * mrSges(6,2) - t536) * t461 + (qJ(5) * mrSges(6,2) + t535) * t460 + t541 * qJDD(4) + (Ifges(3,6) - t511) * qJDD(1) + (-t497 * t512 + t499 * t513 + Ifges(3,5)) * t505;
t398 = mrSges(3,2) * t496 - mrSges(3,3) * t462 + Ifges(3,5) * qJDD(1) - t505 * Ifges(3,6) - qJ(3) * t411 - t497 * t400 + t499 * t401;
t397 = -mrSges(2,2) * g(3) - mrSges(2,3) * t479 + Ifges(2,5) * qJDD(1) - t505 * Ifges(2,6) - qJ(2) * t407 + t500 * t398 - t498 * t399;
t396 = mrSges(2,1) * g(3) + mrSges(2,3) * t480 + t505 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t520 + qJ(2) * t517 + t498 * t398 + t500 * t399;
t1 = [-m(1) * g(1) + t518; -m(1) * g(2) + t532; (-m(1) - m(2)) * g(3) + t520; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t532 - t502 * t396 + t503 * t397; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t518 + t503 * t396 + t502 * t397; pkin(1) * t407 + mrSges(2,1) * t479 - mrSges(2,2) * t480 + t497 * t401 + t499 * t400 + pkin(2) * t506 + qJ(3) * t516 + mrSges(3,1) * t462 - mrSges(3,2) * t463 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t534 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
