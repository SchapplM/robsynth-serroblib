% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:47
% EndTime: 2019-12-31 19:14:52
% DurationCPUTime: 2.49s
% Computational Cost: add. (26021->264), mult. (50170->322), div. (0->0), fcn. (30397->8), ass. (0->105)
t487 = sin(qJ(1));
t491 = cos(qJ(1));
t472 = -t491 * g(1) - t487 * g(2);
t514 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t472;
t513 = -pkin(1) - pkin(6);
t512 = mrSges(2,1) - mrSges(3,2);
t511 = -Ifges(3,4) + Ifges(2,5);
t510 = (Ifges(3,5) - Ifges(2,6));
t471 = t487 * g(1) - t491 * g(2);
t493 = qJD(1) ^ 2;
t499 = -t493 * qJ(2) + qJDD(2) - t471;
t450 = t513 * qJDD(1) + t499;
t486 = sin(qJ(3));
t490 = cos(qJ(3));
t443 = -t490 * g(3) + t486 * t450;
t465 = (mrSges(4,1) * t486 + mrSges(4,2) * t490) * qJD(1);
t507 = qJD(1) * qJD(3);
t475 = t490 * t507;
t467 = -t486 * qJDD(1) - t475;
t508 = qJD(1) * t490;
t470 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t508;
t477 = t486 * qJD(1);
t449 = t513 * t493 - t514;
t505 = t486 * t507;
t468 = t490 * qJDD(1) - t505;
t426 = (-t468 + t505) * pkin(7) + (-t467 + t475) * pkin(3) + t449;
t466 = (pkin(3) * t486 - pkin(7) * t490) * qJD(1);
t492 = qJD(3) ^ 2;
t429 = -t492 * pkin(3) + qJDD(3) * pkin(7) - t466 * t477 + t443;
t485 = sin(qJ(4));
t489 = cos(qJ(4));
t414 = t489 * t426 - t485 * t429;
t463 = t489 * qJD(3) - t485 * t508;
t438 = t463 * qJD(4) + t485 * qJDD(3) + t489 * t468;
t462 = qJDD(4) - t467;
t464 = t485 * qJD(3) + t489 * t508;
t474 = t477 + qJD(4);
t411 = (t463 * t474 - t438) * pkin(8) + (t463 * t464 + t462) * pkin(4) + t414;
t415 = t485 * t426 + t489 * t429;
t437 = -t464 * qJD(4) + t489 * qJDD(3) - t485 * t468;
t448 = t474 * pkin(4) - t464 * pkin(8);
t461 = t463 ^ 2;
t412 = -t461 * pkin(4) + t437 * pkin(8) - t474 * t448 + t415;
t484 = sin(qJ(5));
t488 = cos(qJ(5));
t409 = t488 * t411 - t484 * t412;
t439 = t488 * t463 - t484 * t464;
t418 = t439 * qJD(5) + t484 * t437 + t488 * t438;
t440 = t484 * t463 + t488 * t464;
t423 = -t439 * mrSges(6,1) + t440 * mrSges(6,2);
t473 = qJD(5) + t474;
t430 = -t473 * mrSges(6,2) + t439 * mrSges(6,3);
t457 = qJDD(5) + t462;
t407 = m(6) * t409 + t457 * mrSges(6,1) - t418 * mrSges(6,3) - t440 * t423 + t473 * t430;
t410 = t484 * t411 + t488 * t412;
t417 = -t440 * qJD(5) + t488 * t437 - t484 * t438;
t431 = t473 * mrSges(6,1) - t440 * mrSges(6,3);
t408 = m(6) * t410 - t457 * mrSges(6,2) + t417 * mrSges(6,3) + t439 * t423 - t473 * t431;
t400 = t488 * t407 + t484 * t408;
t441 = -t463 * mrSges(5,1) + t464 * mrSges(5,2);
t444 = -t474 * mrSges(5,2) + t463 * mrSges(5,3);
t398 = m(5) * t414 + t462 * mrSges(5,1) - t438 * mrSges(5,3) - t464 * t441 + t474 * t444 + t400;
t445 = t474 * mrSges(5,1) - t464 * mrSges(5,3);
t501 = -t484 * t407 + t488 * t408;
t399 = m(5) * t415 - t462 * mrSges(5,2) + t437 * mrSges(5,3) + t463 * t441 - t474 * t445 + t501;
t502 = -t485 * t398 + t489 * t399;
t393 = m(4) * t443 - qJDD(3) * mrSges(4,2) + t467 * mrSges(4,3) - qJD(3) * t470 - t465 * t477 + t502;
t442 = t486 * g(3) + t490 * t450;
t469 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t477;
t428 = -qJDD(3) * pkin(3) - t492 * pkin(7) + t466 * t508 - t442;
t413 = -t437 * pkin(4) - t461 * pkin(8) + t464 * t448 + t428;
t496 = m(6) * t413 - t417 * mrSges(6,1) + t418 * mrSges(6,2) - t439 * t430 + t440 * t431;
t494 = -m(5) * t428 + t437 * mrSges(5,1) - t438 * mrSges(5,2) + t463 * t444 - t464 * t445 - t496;
t403 = m(4) * t442 + qJDD(3) * mrSges(4,1) - t468 * mrSges(4,3) + qJD(3) * t469 - t465 * t508 + t494;
t387 = t486 * t393 + t490 * t403;
t452 = -qJDD(1) * pkin(1) + t499;
t498 = -m(3) * t452 + (t493 * mrSges(3,3)) - t387;
t385 = m(2) * t471 - (t493 * mrSges(2,2)) + t512 * qJDD(1) + t498;
t451 = t493 * pkin(1) + t514;
t394 = t489 * t398 + t485 * t399;
t497 = -m(4) * t449 + t467 * mrSges(4,1) - t468 * mrSges(4,2) - t469 * t477 - t470 * t508 - t394;
t495 = -m(3) * t451 + (t493 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t497;
t391 = m(2) * t472 - (t493 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t495;
t509 = t491 * t385 + t487 * t391;
t504 = -t487 * t385 + t491 * t391;
t503 = t490 * t393 - t486 * t403;
t456 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t490 - Ifges(4,4) * t486) * qJD(1);
t455 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t490 - Ifges(4,2) * t486) * qJD(1);
t454 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t490 - Ifges(4,6) * t486) * qJD(1);
t434 = Ifges(5,1) * t464 + Ifges(5,4) * t463 + Ifges(5,5) * t474;
t433 = Ifges(5,4) * t464 + Ifges(5,2) * t463 + Ifges(5,6) * t474;
t432 = Ifges(5,5) * t464 + Ifges(5,6) * t463 + Ifges(5,3) * t474;
t421 = Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * t473;
t420 = Ifges(6,4) * t440 + Ifges(6,2) * t439 + Ifges(6,6) * t473;
t419 = Ifges(6,5) * t440 + Ifges(6,6) * t439 + Ifges(6,3) * t473;
t402 = mrSges(6,2) * t413 - mrSges(6,3) * t409 + Ifges(6,1) * t418 + Ifges(6,4) * t417 + Ifges(6,5) * t457 + t439 * t419 - t473 * t420;
t401 = -mrSges(6,1) * t413 + mrSges(6,3) * t410 + Ifges(6,4) * t418 + Ifges(6,2) * t417 + Ifges(6,6) * t457 - t440 * t419 + t473 * t421;
t388 = mrSges(5,2) * t428 - mrSges(5,3) * t414 + Ifges(5,1) * t438 + Ifges(5,4) * t437 + Ifges(5,5) * t462 - pkin(8) * t400 - t484 * t401 + t488 * t402 + t463 * t432 - t474 * t433;
t386 = -m(3) * g(3) + t503;
t383 = -mrSges(5,1) * t428 + mrSges(5,3) * t415 + Ifges(5,4) * t438 + Ifges(5,2) * t437 + Ifges(5,6) * t462 - pkin(4) * t496 + pkin(8) * t501 + t488 * t401 + t484 * t402 - t464 * t432 + t474 * t434;
t382 = Ifges(4,4) * t468 + Ifges(4,2) * t467 + Ifges(4,6) * qJDD(3) - t454 * t508 + qJD(3) * t456 - mrSges(4,1) * t449 + mrSges(4,3) * t443 - Ifges(5,5) * t438 - Ifges(5,6) * t437 - Ifges(5,3) * t462 - t464 * t433 + t463 * t434 - mrSges(5,1) * t414 + mrSges(5,2) * t415 - Ifges(6,5) * t418 - Ifges(6,6) * t417 - Ifges(6,3) * t457 - t440 * t420 + t439 * t421 - mrSges(6,1) * t409 + mrSges(6,2) * t410 - pkin(4) * t400 - pkin(3) * t394;
t381 = mrSges(4,2) * t449 - mrSges(4,3) * t442 + Ifges(4,1) * t468 + Ifges(4,4) * t467 + Ifges(4,5) * qJDD(3) - pkin(7) * t394 - qJD(3) * t455 - t485 * t383 + t489 * t388 - t454 * t477;
t380 = -qJ(2) * t386 - mrSges(2,3) * t471 + pkin(2) * t387 + mrSges(3,1) * t452 + t485 * t388 + t489 * t383 + pkin(3) * t494 + pkin(7) * t502 + mrSges(4,1) * t442 - mrSges(4,2) * t443 + Ifges(4,5) * t468 + Ifges(4,6) * t467 + Ifges(4,3) * qJDD(3) + (t510 * t493) + t511 * qJDD(1) + (t490 * t455 + t486 * t456) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t379 = -mrSges(3,1) * t451 + mrSges(2,3) * t472 - pkin(1) * t386 - pkin(2) * t497 - pkin(6) * t503 + t512 * g(3) - t510 * qJDD(1) - t486 * t381 - t490 * t382 + t511 * t493;
t1 = [-m(1) * g(1) + t504; -m(1) * g(2) + t509; (-m(1) - m(2) - m(3)) * g(3) + t503; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t509 - t487 * t379 + t491 * t380; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t504 + t491 * t379 + t487 * t380; pkin(1) * t498 + qJ(2) * t495 + t490 * t381 - t486 * t382 - pkin(6) * t387 + mrSges(2,1) * t471 - mrSges(2,2) * t472 + mrSges(3,2) * t452 - mrSges(3,3) * t451 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
