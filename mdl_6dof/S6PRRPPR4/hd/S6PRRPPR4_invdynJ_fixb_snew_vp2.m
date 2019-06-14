% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:08:50
% EndTime: 2019-05-05 03:08:54
% DurationCPUTime: 2.23s
% Computational Cost: add. (12476->270), mult. (26032->328), div. (0->0), fcn. (17443->12), ass. (0->116)
t532 = Ifges(5,1) + Ifges(6,1);
t518 = Ifges(5,4) - Ifges(6,5);
t531 = Ifges(5,5) + Ifges(6,4);
t530 = Ifges(5,2) + Ifges(6,3);
t527 = Ifges(5,6) - Ifges(6,6);
t529 = -2 * qJD(4);
t528 = Ifges(6,2) + Ifges(5,3);
t481 = sin(pkin(10));
t483 = cos(pkin(10));
t468 = g(1) * t481 - g(2) * t483;
t479 = -g(3) + qJDD(1);
t482 = sin(pkin(6));
t484 = cos(pkin(6));
t526 = t468 * t484 + t479 * t482;
t480 = sin(pkin(11));
t486 = sin(qJ(3));
t510 = qJD(2) * t486;
t516 = cos(pkin(11));
t456 = -qJD(3) * t516 + t480 * t510;
t457 = t480 * qJD(3) + t510 * t516;
t489 = cos(qJ(3));
t509 = qJD(2) * t489;
t523 = t530 * t456 - t518 * t457 + t527 * t509;
t522 = t518 * t456 - t532 * t457 + t531 * t509;
t469 = -g(1) * t483 - g(2) * t481;
t487 = sin(qJ(2));
t490 = cos(qJ(2));
t419 = t490 * t469 + t487 * t526;
t492 = qJD(2) ^ 2;
t413 = -pkin(2) * t492 + qJDD(2) * pkin(8) + t419;
t444 = -t468 * t482 + t479 * t484;
t408 = t489 * t413 + t486 * t444;
t464 = (-pkin(3) * t489 - qJ(4) * t486) * qJD(2);
t491 = qJD(3) ^ 2;
t393 = -pkin(3) * t491 + qJDD(3) * qJ(4) + t464 * t509 + t408;
t418 = -t487 * t469 + t490 * t526;
t412 = -qJDD(2) * pkin(2) - t492 * pkin(8) - t418;
t507 = qJD(2) * qJD(3);
t505 = t489 * t507;
t466 = qJDD(2) * t486 + t505;
t475 = t486 * t507;
t467 = qJDD(2) * t489 - t475;
t399 = (-t466 - t505) * qJ(4) + (-t467 + t475) * pkin(3) + t412;
t386 = -t480 * t393 + t399 * t516 + t457 * t529;
t442 = t480 * qJDD(3) + t466 * t516;
t407 = -t486 * t413 + t489 * t444;
t495 = qJDD(3) * pkin(3) + t491 * qJ(4) - t464 * t510 - qJDD(4) + t407;
t506 = t456 * t509;
t521 = -(t442 + t506) * qJ(5) - t495;
t520 = -2 * qJD(5);
t519 = -mrSges(5,3) - mrSges(6,2);
t517 = t486 * Ifges(4,4);
t514 = t489 ^ 2 * t492;
t512 = t456 * t527 - t457 * t531 + t509 * t528;
t431 = mrSges(6,1) * t456 - mrSges(6,3) * t457;
t511 = -mrSges(5,1) * t456 - mrSges(5,2) * t457 - t431;
t387 = t516 * t393 + t480 * t399 + t456 * t529;
t465 = (-mrSges(4,1) * t489 + mrSges(4,2) * t486) * qJD(2);
t470 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t510;
t439 = -mrSges(5,1) * t509 - mrSges(5,3) * t457;
t440 = mrSges(6,1) * t509 + mrSges(6,2) * t457;
t441 = -qJDD(3) * t516 + t480 * t466;
t430 = pkin(4) * t456 - qJ(5) * t457;
t382 = -pkin(4) * t514 - t467 * qJ(5) - t456 * t430 + t509 * t520 + t387;
t384 = t467 * pkin(4) - qJ(5) * t514 + t457 * t430 + qJDD(5) - t386;
t379 = (-t442 + t506) * pkin(9) + (t456 * t457 + t467) * pkin(5) + t384;
t443 = pkin(5) * t509 - pkin(9) * t457;
t454 = t456 ^ 2;
t380 = -pkin(5) * t454 + pkin(9) * t441 - t443 * t509 + t382;
t485 = sin(qJ(6));
t488 = cos(qJ(6));
t377 = t379 * t488 - t380 * t485;
t428 = t456 * t488 - t457 * t485;
t401 = qJD(6) * t428 + t441 * t485 + t442 * t488;
t429 = t456 * t485 + t457 * t488;
t409 = -mrSges(7,1) * t428 + mrSges(7,2) * t429;
t473 = qJD(6) + t509;
t414 = -mrSges(7,2) * t473 + mrSges(7,3) * t428;
t461 = qJDD(6) + t467;
t373 = m(7) * t377 + mrSges(7,1) * t461 - mrSges(7,3) * t401 - t409 * t429 + t414 * t473;
t378 = t379 * t485 + t380 * t488;
t400 = -qJD(6) * t429 + t441 * t488 - t442 * t485;
t415 = mrSges(7,1) * t473 - mrSges(7,3) * t429;
t374 = m(7) * t378 - mrSges(7,2) * t461 + mrSges(7,3) * t400 + t409 * t428 - t415 * t473;
t502 = -t485 * t373 + t488 * t374;
t497 = m(6) * t382 - t467 * mrSges(6,3) + t502;
t364 = m(5) * t387 + t467 * mrSges(5,2) + t511 * t456 + t519 * t441 + (t439 - t440) * t509 + t497;
t437 = -mrSges(6,2) * t456 - mrSges(6,3) * t509;
t438 = mrSges(5,2) * t509 - mrSges(5,3) * t456;
t367 = t488 * t373 + t485 * t374;
t496 = -m(6) * t384 - t467 * mrSges(6,1) - t367;
t365 = m(5) * t386 - t467 * mrSges(5,1) + t511 * t457 + t519 * t442 + (-t437 - t438) * t509 + t496;
t503 = t516 * t364 - t365 * t480;
t361 = m(4) * t408 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t467 - qJD(3) * t470 + t465 * t509 + t503;
t388 = t457 * t520 + (-t457 * t509 + t441) * pkin(4) + t521;
t385 = -t454 * pkin(9) + (-pkin(4) - pkin(5)) * t441 + (pkin(4) * t509 + (2 * qJD(5)) + t443) * t457 - t521;
t499 = -m(7) * t385 + t400 * mrSges(7,1) - t401 * mrSges(7,2) + t428 * t414 - t429 * t415;
t375 = m(6) * t388 + t441 * mrSges(6,1) - t442 * mrSges(6,3) + t456 * t437 - t457 * t440 + t499;
t371 = -m(5) * t495 + t441 * mrSges(5,1) + t442 * mrSges(5,2) + t456 * t438 + t457 * t439 + t375;
t471 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t509;
t370 = m(4) * t407 + qJDD(3) * mrSges(4,1) - t466 * mrSges(4,3) + qJD(3) * t471 - t465 * t510 - t371;
t504 = t489 * t361 - t370 * t486;
t362 = t480 * t364 + t365 * t516;
t403 = Ifges(7,4) * t429 + Ifges(7,2) * t428 + Ifges(7,6) * t473;
t404 = Ifges(7,1) * t429 + Ifges(7,4) * t428 + Ifges(7,5) * t473;
t494 = mrSges(7,1) * t377 - mrSges(7,2) * t378 + Ifges(7,5) * t401 + Ifges(7,6) * t400 + Ifges(7,3) * t461 + t429 * t403 - t428 * t404;
t493 = -m(4) * t412 + t467 * mrSges(4,1) - t466 * mrSges(4,2) - t470 * t510 + t471 * t509 - t362;
t451 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t486 + Ifges(4,4) * t489) * qJD(2);
t450 = Ifges(4,6) * qJD(3) + (Ifges(4,2) * t489 + t517) * qJD(2);
t402 = Ifges(7,5) * t429 + Ifges(7,6) * t428 + Ifges(7,3) * t473;
t369 = mrSges(7,2) * t385 - mrSges(7,3) * t377 + Ifges(7,1) * t401 + Ifges(7,4) * t400 + Ifges(7,5) * t461 + t402 * t428 - t403 * t473;
t368 = -mrSges(7,1) * t385 + mrSges(7,3) * t378 + Ifges(7,4) * t401 + Ifges(7,2) * t400 + Ifges(7,6) * t461 - t402 * t429 + t404 * t473;
t366 = t442 * mrSges(6,2) + t457 * t431 + t437 * t509 - t496;
t359 = -mrSges(5,2) * t495 + mrSges(6,2) * t384 - mrSges(5,3) * t386 - mrSges(6,3) * t388 - pkin(9) * t367 - qJ(5) * t375 - t485 * t368 + t488 * t369 - t518 * t441 + t532 * t442 + t512 * t456 - t467 * t531 - t523 * t509;
t358 = mrSges(5,1) * t495 - mrSges(6,1) * t388 + mrSges(6,2) * t382 + mrSges(5,3) * t387 - pkin(4) * t375 - pkin(5) * t499 - pkin(9) * t502 - t488 * t368 - t485 * t369 - t530 * t441 + t518 * t442 + t512 * t457 - t467 * t527 + t522 * t509;
t1 = [m(2) * t479 + t484 * (m(3) * t444 + t361 * t486 + t370 * t489) + (t487 * (m(3) * t419 - mrSges(3,1) * t492 - qJDD(2) * mrSges(3,2) + t504) + t490 * (m(3) * t418 + qJDD(2) * mrSges(3,1) - t492 * mrSges(3,2) + t493)) * t482; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t418 - mrSges(3,2) * t419 + t486 * (mrSges(4,2) * t412 - mrSges(4,3) * t407 + Ifges(4,1) * t466 + Ifges(4,5) * qJDD(3) - qJ(4) * t362 - qJD(3) * t450 - t480 * t358 + t359 * t516) + t489 * (Ifges(4,6) * qJDD(3) + t494 + Ifges(4,4) * t466 + qJD(3) * t451 - mrSges(4,1) * t412 + mrSges(4,3) * t408 - mrSges(5,1) * t386 + mrSges(5,2) * t387 - mrSges(6,3) * t382 + mrSges(6,1) * t384 - pkin(3) * t362 + pkin(4) * t366 + pkin(5) * t367 - qJ(5) * (-t440 * t509 + t497) + t523 * t457 + (qJ(5) * t431 + t522) * t456 - t531 * t442 + (qJ(5) * mrSges(6,2) + t527) * t441) + pkin(2) * t493 + pkin(8) * t504 + (t517 + t489 * (Ifges(4,2) + t528)) * t467; Ifges(4,5) * t466 + Ifges(4,6) * t467 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t407 - mrSges(4,2) * t408 + t480 * t359 + t516 * t358 - pkin(3) * t371 + qJ(4) * t503 + (t450 * t486 - t451 * t489) * qJD(2); t371; t366; t494;];
tauJ  = t1;
