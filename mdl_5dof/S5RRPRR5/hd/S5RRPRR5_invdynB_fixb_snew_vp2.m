% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:46
% EndTime: 2020-01-03 12:03:51
% DurationCPUTime: 4.53s
% Computational Cost: add. (77172->250), mult. (106816->313), div. (0->0), fcn. (71884->10), ass. (0->109)
t508 = qJD(1) + qJD(2);
t502 = t508 ^ 2;
t510 = cos(pkin(9));
t543 = pkin(3) * t510;
t509 = sin(pkin(9));
t542 = mrSges(4,2) * t509;
t506 = t510 ^ 2;
t541 = t502 * t506;
t504 = qJDD(1) + qJDD(2);
t540 = t504 * t510;
t527 = Ifges(4,5) * t509 + Ifges(4,6) * t510;
t539 = t502 * t527;
t514 = sin(qJ(1));
t518 = cos(qJ(1));
t496 = -g(2) * t514 + t518 * g(3);
t519 = qJD(1) ^ 2;
t497 = -g(2) * t518 - g(3) * t514;
t492 = qJDD(1) * pkin(1) + t497;
t493 = -pkin(1) * t519 + t496;
t513 = sin(qJ(2));
t517 = cos(qJ(2));
t480 = t513 * t492 + t517 * t493;
t478 = -pkin(2) * t502 + qJ(3) * t504 + t480;
t537 = qJD(3) * t508;
t535 = -t510 * g(1) - 0.2e1 * t509 * t537;
t459 = (-pkin(7) * t504 + t502 * t543 - t478) * t509 + t535;
t463 = -t509 * g(1) + (t478 + 0.2e1 * t537) * t510;
t460 = -pkin(3) * t541 + pkin(7) * t540 + t463;
t512 = sin(qJ(4));
t516 = cos(qJ(4));
t444 = t516 * t459 - t512 * t460;
t523 = t509 * t516 + t510 * t512;
t522 = -t509 * t512 + t510 * t516;
t485 = t522 * t508;
t536 = t485 * qJD(4);
t477 = t504 * t523 + t536;
t486 = t523 * t508;
t440 = (-t477 + t536) * pkin(8) + (t485 * t486 + qJDD(4)) * pkin(4) + t444;
t445 = t512 * t459 + t516 * t460;
t476 = -t486 * qJD(4) + t504 * t522;
t483 = qJD(4) * pkin(4) - pkin(8) * t486;
t484 = t485 ^ 2;
t441 = -pkin(4) * t484 + pkin(8) * t476 - qJD(4) * t483 + t445;
t511 = sin(qJ(5));
t515 = cos(qJ(5));
t438 = t440 * t515 - t441 * t511;
t469 = t485 * t515 - t486 * t511;
t449 = qJD(5) * t469 + t476 * t511 + t477 * t515;
t470 = t485 * t511 + t486 * t515;
t455 = -mrSges(6,1) * t469 + mrSges(6,2) * t470;
t507 = qJD(4) + qJD(5);
t464 = -mrSges(6,2) * t507 + mrSges(6,3) * t469;
t503 = qJDD(4) + qJDD(5);
t436 = m(6) * t438 + mrSges(6,1) * t503 - mrSges(6,3) * t449 - t455 * t470 + t464 * t507;
t439 = t440 * t511 + t441 * t515;
t448 = -qJD(5) * t470 + t476 * t515 - t477 * t511;
t465 = mrSges(6,1) * t507 - mrSges(6,3) * t470;
t437 = m(6) * t439 - mrSges(6,2) * t503 + mrSges(6,3) * t448 + t455 * t469 - t465 * t507;
t428 = t515 * t436 + t511 * t437;
t474 = -mrSges(5,1) * t485 + mrSges(5,2) * t486;
t481 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t485;
t426 = m(5) * t444 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t477 + qJD(4) * t481 - t474 * t486 + t428;
t482 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t486;
t530 = -t436 * t511 + t515 * t437;
t427 = m(5) * t445 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t476 - qJD(4) * t482 + t474 * t485 + t530;
t422 = t516 * t426 + t512 * t427;
t462 = -t509 * t478 + t535;
t524 = mrSges(4,3) * t504 + (-mrSges(4,1) * t510 + t542) * t502;
t420 = m(4) * t462 - t509 * t524 + t422;
t531 = -t512 * t426 + t516 * t427;
t421 = m(4) * t463 + t510 * t524 + t531;
t532 = -t420 * t509 + t510 * t421;
t413 = m(3) * t480 - mrSges(3,1) * t502 - mrSges(3,2) * t504 + t532;
t479 = t517 * t492 - t513 * t493;
t526 = qJDD(3) - t479;
t475 = -t504 * pkin(2) - t502 * qJ(3) + t526;
t505 = t509 ^ 2;
t461 = (-pkin(2) - t543) * t504 + (-qJ(3) + (-t505 - t506) * pkin(7)) * t502 + t526;
t443 = -t476 * pkin(4) - t484 * pkin(8) + t486 * t483 + t461;
t525 = m(6) * t443 - t448 * mrSges(6,1) + t449 * mrSges(6,2) - t469 * t464 + t470 * t465;
t521 = m(5) * t461 - t476 * mrSges(5,1) + t477 * mrSges(5,2) - t485 * t481 + t486 * t482 + t525;
t520 = -m(4) * t475 + mrSges(4,1) * t540 - t521 + (t502 * t505 + t541) * mrSges(4,3);
t432 = t520 + (mrSges(3,1) - t542) * t504 - t502 * mrSges(3,2) + m(3) * t479;
t533 = t517 * t413 - t432 * t513;
t408 = m(2) * t496 - mrSges(2,1) * t519 - qJDD(1) * mrSges(2,2) + t533;
t410 = t513 * t413 + t517 * t432;
t409 = m(2) * t497 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t519 + t410;
t538 = t514 * t408 + t518 * t409;
t414 = t510 * t420 + t509 * t421;
t534 = -t408 * t518 + t514 * t409;
t529 = Ifges(4,1) * t509 + Ifges(4,4) * t510;
t528 = Ifges(4,4) * t509 + Ifges(4,2) * t510;
t468 = Ifges(5,1) * t486 + Ifges(5,4) * t485 + Ifges(5,5) * qJD(4);
t467 = Ifges(5,4) * t486 + Ifges(5,2) * t485 + Ifges(5,6) * qJD(4);
t466 = Ifges(5,5) * t486 + Ifges(5,6) * t485 + Ifges(5,3) * qJD(4);
t452 = Ifges(6,1) * t470 + Ifges(6,4) * t469 + Ifges(6,5) * t507;
t451 = Ifges(6,4) * t470 + Ifges(6,2) * t469 + Ifges(6,6) * t507;
t450 = Ifges(6,5) * t470 + Ifges(6,6) * t469 + Ifges(6,3) * t507;
t430 = mrSges(6,2) * t443 - mrSges(6,3) * t438 + Ifges(6,1) * t449 + Ifges(6,4) * t448 + Ifges(6,5) * t503 + t450 * t469 - t451 * t507;
t429 = -mrSges(6,1) * t443 + mrSges(6,3) * t439 + Ifges(6,4) * t449 + Ifges(6,2) * t448 + Ifges(6,6) * t503 - t450 * t470 + t452 * t507;
t416 = mrSges(5,2) * t461 - mrSges(5,3) * t444 + Ifges(5,1) * t477 + Ifges(5,4) * t476 + Ifges(5,5) * qJDD(4) - pkin(8) * t428 - qJD(4) * t467 - t429 * t511 + t430 * t515 + t466 * t485;
t415 = -mrSges(5,1) * t461 + mrSges(5,3) * t445 + Ifges(5,4) * t477 + Ifges(5,2) * t476 + Ifges(5,6) * qJDD(4) - pkin(4) * t525 + pkin(8) * t530 + qJD(4) * t468 + t515 * t429 + t511 * t430 - t486 * t466;
t404 = mrSges(4,2) * t475 - mrSges(4,3) * t462 - pkin(7) * t422 - t512 * t415 + t516 * t416 + t504 * t529 + t510 * t539;
t403 = -mrSges(4,1) * t475 + mrSges(4,3) * t463 - pkin(3) * t521 + pkin(7) * t531 + t516 * t415 + t512 * t416 + t504 * t528 - t509 * t539;
t402 = mrSges(3,1) * g(1) - Ifges(5,3) * qJDD(4) - Ifges(6,3) * t503 + t485 * t468 - t486 * t467 - Ifges(5,6) * t476 - Ifges(5,5) * t477 + mrSges(3,3) * t480 + t469 * t452 - t470 * t451 - mrSges(4,1) * t462 + mrSges(4,2) * t463 - mrSges(5,1) * t444 + mrSges(5,2) * t445 - Ifges(6,6) * t448 - Ifges(6,5) * t449 - mrSges(6,1) * t438 + mrSges(6,2) * t439 - pkin(4) * t428 - pkin(3) * t422 - pkin(2) * t414 + (Ifges(3,6) - t527) * t504 + (-t509 * t528 + t510 * t529 + Ifges(3,5)) * t502;
t401 = -mrSges(3,2) * g(1) - mrSges(3,3) * t479 + Ifges(3,5) * t504 - Ifges(3,6) * t502 - qJ(3) * t414 - t403 * t509 + t404 * t510;
t400 = -mrSges(2,2) * g(1) - mrSges(2,3) * t497 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t519 - pkin(6) * t410 + t401 * t517 - t402 * t513;
t399 = Ifges(2,6) * qJDD(1) + t519 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t496 + t513 * t401 + t517 * t402 - pkin(1) * (-m(3) * g(1) + t414) + pkin(6) * t533;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t414; -m(1) * g(2) + t538; -m(1) * g(3) + t534; pkin(1) * t410 + t509 * t404 + t510 * t403 + pkin(2) * (-t504 * t542 + t520) + qJ(3) * t532 + mrSges(3,1) * t479 - mrSges(3,2) * t480 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t504 + mrSges(2,1) * t497 - mrSges(2,2) * t496 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t534 + t518 * t399 + t514 * t400; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t538 + t399 * t514 - t400 * t518;];
tauB = t1;
