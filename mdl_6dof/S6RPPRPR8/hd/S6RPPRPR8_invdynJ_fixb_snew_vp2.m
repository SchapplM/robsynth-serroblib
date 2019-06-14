% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-05-05 14:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:37:48
% EndTime: 2019-05-05 14:37:51
% DurationCPUTime: 1.60s
% Computational Cost: add. (9124->242), mult. (20880->292), div. (0->0), fcn. (13708->8), ass. (0->111)
t524 = Ifges(5,1) + Ifges(6,2);
t514 = Ifges(5,4) + Ifges(6,6);
t513 = Ifges(5,5) - Ifges(6,4);
t523 = -Ifges(5,2) - Ifges(6,3);
t512 = Ifges(5,6) - Ifges(6,5);
t522 = Ifges(5,3) + Ifges(6,1);
t476 = qJD(1) ^ 2;
t471 = sin(qJ(1));
t474 = cos(qJ(1));
t496 = t471 * g(1) - t474 * g(2);
t486 = -t476 * qJ(2) + qJDD(2) - t496;
t511 = -pkin(1) - qJ(3);
t521 = -(2 * qJD(1) * qJD(3)) + t511 * qJDD(1) + t486;
t469 = cos(pkin(9));
t520 = t469 ^ 2;
t493 = -t474 * g(1) - t471 * g(2);
t519 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t493;
t518 = -2 * qJD(5);
t517 = sin(qJ(4));
t516 = pkin(3) * t476;
t515 = mrSges(5,1) - mrSges(6,2);
t468 = sin(pkin(9));
t427 = t468 * g(3) + t521 * t469;
t405 = (-pkin(7) * qJDD(1) - t468 * t516) * t469 + t427;
t428 = -g(3) * t469 + t521 * t468;
t464 = t468 ^ 2;
t500 = qJDD(1) * t468;
t406 = -pkin(7) * t500 - t464 * t516 + t428;
t473 = cos(qJ(4));
t387 = t473 * t405 - t517 * t406;
t488 = t468 * t473 + t517 * t469;
t450 = t488 * qJD(1);
t497 = t468 * t517;
t504 = qJD(1) * t469;
t451 = -qJD(1) * t497 + t473 * t504;
t421 = mrSges(5,1) * t450 + mrSges(5,2) * t451;
t499 = qJDD(1) * t469;
t503 = qJD(4) * t450;
t430 = -qJDD(1) * t497 + t473 * t499 - t503;
t502 = t451 * qJD(4);
t429 = t488 * qJDD(1) + t502;
t442 = pkin(5) * t451 - qJD(4) * pkin(8);
t449 = t450 ^ 2;
t485 = qJDD(3) + t519;
t505 = -t464 - t520;
t410 = pkin(3) * t500 + (t505 * pkin(7) + t511) * t476 + t485;
t477 = pkin(4) * t502 + t451 * t518 + (-t430 + t503) * qJ(5) + t410;
t379 = -t449 * pkin(5) - t451 * t442 + (pkin(4) + pkin(8)) * t429 + t477;
t420 = pkin(4) * t450 - qJ(5) * t451;
t475 = qJD(4) ^ 2;
t386 = -qJDD(4) * pkin(4) - t475 * qJ(5) + t451 * t420 + qJDD(5) - t387;
t380 = (t450 * t451 - qJDD(4)) * pkin(8) + (t430 + t503) * pkin(5) + t386;
t470 = sin(qJ(6));
t472 = cos(qJ(6));
t377 = -t379 * t470 + t380 * t472;
t432 = -qJD(4) * t470 + t450 * t472;
t398 = qJD(6) * t432 + qJDD(4) * t472 + t429 * t470;
t433 = qJD(4) * t472 + t450 * t470;
t400 = -mrSges(7,1) * t432 + mrSges(7,2) * t433;
t447 = qJD(6) + t451;
t407 = -mrSges(7,2) * t447 + mrSges(7,3) * t432;
t426 = qJDD(6) + t430;
t374 = m(7) * t377 + mrSges(7,1) * t426 - mrSges(7,3) * t398 - t400 * t433 + t407 * t447;
t378 = t379 * t472 + t380 * t470;
t397 = -qJD(6) * t433 - qJDD(4) * t470 + t429 * t472;
t408 = mrSges(7,1) * t447 - mrSges(7,3) * t433;
t375 = m(7) * t378 - mrSges(7,2) * t426 + mrSges(7,3) * t397 + t400 * t432 - t408 * t447;
t367 = t472 * t374 + t470 * t375;
t422 = -mrSges(6,2) * t450 - mrSges(6,3) * t451;
t483 = -m(6) * t386 - t430 * mrSges(6,1) - t451 * t422 - t367;
t440 = mrSges(6,1) * t450 - qJD(4) * mrSges(6,3);
t506 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t450 - t440;
t364 = m(5) * t387 - t430 * mrSges(5,3) + t506 * qJD(4) + t515 * qJDD(4) - t451 * t421 + t483;
t388 = t517 * t405 + t473 * t406;
t439 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t451;
t482 = -t475 * pkin(4) + qJDD(4) * qJ(5) - t450 * t420 + t388;
t385 = qJD(4) * t518 - t482;
t441 = mrSges(6,1) * t451 + qJD(4) * mrSges(6,2);
t382 = -t429 * pkin(5) - t449 * pkin(8) + ((2 * qJD(5)) + t442) * qJD(4) + t482;
t484 = -m(7) * t382 + t397 * mrSges(7,1) - t398 * mrSges(7,2) + t432 * t407 - t433 * t408;
t480 = -m(6) * t385 + qJDD(4) * mrSges(6,3) + qJD(4) * t441 - t484;
t372 = m(5) * t388 - qJDD(4) * mrSges(5,2) - qJD(4) * t439 + (-t421 - t422) * t450 + (-mrSges(5,3) - mrSges(6,1)) * t429 + t480;
t510 = t473 * t364 + t517 * t372;
t509 = -t522 * qJD(4) + t512 * t450 - t513 * t451;
t508 = t512 * qJD(4) + t523 * t450 + t514 * t451;
t507 = t513 * qJD(4) - t514 * t450 + t524 * t451;
t495 = t505 * mrSges(4,3);
t494 = -t470 * t374 + t472 * t375;
t492 = -t517 * t364 + t473 * t372;
t490 = -qJDD(1) * mrSges(4,3) - t476 * (mrSges(4,1) * t468 + mrSges(4,2) * t469);
t491 = (m(4) * t427 + t490 * t469 + t510) * t469 + (m(4) * t428 + t490 * t468 + t492) * t468;
t384 = t429 * pkin(4) + t477;
t487 = m(6) * t384 - t430 * mrSges(6,3) - t451 * t441 + t494;
t391 = Ifges(7,4) * t433 + Ifges(7,2) * t432 + Ifges(7,6) * t447;
t392 = Ifges(7,1) * t433 + Ifges(7,4) * t432 + Ifges(7,5) * t447;
t481 = mrSges(7,1) * t377 - mrSges(7,2) * t378 + Ifges(7,5) * t398 + Ifges(7,6) * t397 + Ifges(7,3) * t426 + t433 * t391 - t432 * t392;
t479 = m(5) * t410 + t430 * mrSges(5,2) + t515 * t429 + t451 * t439 + t506 * t450 + t487;
t436 = t511 * t476 + t485;
t478 = m(4) * t436 + mrSges(4,1) * t500 + mrSges(4,2) * t499 + t479;
t453 = (Ifges(4,5) * t469 - Ifges(4,6) * t468) * qJD(1);
t448 = -qJDD(1) * pkin(1) + t486;
t444 = t476 * pkin(1) - t519;
t390 = Ifges(7,5) * t433 + Ifges(7,6) * t432 + Ifges(7,3) * t447;
t369 = mrSges(7,2) * t382 - mrSges(7,3) * t377 + Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t426 + t390 * t432 - t391 * t447;
t368 = -mrSges(7,1) * t382 + mrSges(7,3) * t378 + Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t426 - t390 * t433 + t392 * t447;
t366 = qJDD(4) * mrSges(6,2) + qJD(4) * t440 - t483;
t365 = -t429 * mrSges(6,2) - t450 * t440 + t487;
t360 = mrSges(6,1) * t386 + mrSges(5,2) * t410 - mrSges(5,3) * t387 - mrSges(6,3) * t384 + pkin(5) * t367 - qJ(5) * t365 - t508 * qJD(4) + t513 * qJDD(4) - t514 * t429 + t524 * t430 + t509 * t450 + t481;
t359 = m(3) * t448 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t476 + t491;
t358 = -mrSges(5,1) * t410 - mrSges(6,1) * t385 + mrSges(6,2) * t384 + mrSges(5,3) * t388 - pkin(4) * t365 - pkin(5) * t484 - pkin(8) * t494 + t507 * qJD(4) + t512 * qJDD(4) - t472 * t368 - t470 * t369 + t523 * t429 + t514 * t430 + t509 * t451;
t1 = [mrSges(2,1) * t496 - mrSges(2,2) * t493 + mrSges(3,2) * t448 - mrSges(3,3) * t444 + t469 * (-t468 * qJD(1) * t453 + mrSges(4,2) * t436 - mrSges(4,3) * t427 - pkin(7) * t510 - t517 * t358 + t473 * t360) - t468 * (-mrSges(4,1) * t436 + mrSges(4,3) * t428 - pkin(3) * t479 + pkin(7) * t492 + t473 * t358 + t517 * t360 - t453 * t504) - qJ(3) * t491 - pkin(1) * t359 + qJ(2) * (t478 - m(3) * t444 + (mrSges(3,2) + t495) * t476) + (Ifges(4,1) * t520 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t469 + Ifges(4,2) * t468) * t468) * qJDD(1); t359; t476 * t495 + t478; mrSges(5,1) * t387 - mrSges(5,2) * t388 + mrSges(6,2) * t386 - mrSges(6,3) * t385 + t472 * t369 - t470 * t368 - pkin(8) * t367 - pkin(4) * t366 + qJ(5) * t480 + t508 * t451 + (-qJ(5) * t422 + t507) * t450 + t513 * t430 + (-mrSges(6,1) * qJ(5) - t512) * t429 + t522 * qJDD(4); t366; t481;];
tauJ  = t1;
