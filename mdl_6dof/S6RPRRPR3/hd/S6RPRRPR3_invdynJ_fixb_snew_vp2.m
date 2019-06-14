% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-05 22:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:15:31
% EndTime: 2019-05-05 22:15:35
% DurationCPUTime: 2.17s
% Computational Cost: add. (14912->272), mult. (28410->327), div. (0->0), fcn. (17720->10), ass. (0->115)
t514 = Ifges(5,1) + Ifges(6,1);
t504 = Ifges(5,4) - Ifges(6,5);
t503 = Ifges(5,5) + Ifges(6,4);
t513 = -Ifges(5,2) - Ifges(6,3);
t502 = Ifges(5,6) - Ifges(6,6);
t511 = Ifges(5,3) + Ifges(6,2);
t470 = sin(qJ(4));
t471 = sin(qJ(3));
t496 = qJD(1) * t471;
t506 = cos(qJ(4));
t446 = -t506 * qJD(3) + t470 * t496;
t474 = cos(qJ(3));
t494 = qJD(1) * qJD(3);
t492 = t474 * t494;
t452 = qJDD(1) * t471 + t492;
t415 = -t446 * qJD(4) + t470 * qJDD(3) + t506 * t452;
t447 = t470 * qJD(3) + t506 * t496;
t423 = mrSges(6,1) * t446 - mrSges(6,3) * t447;
t472 = sin(qJ(1));
t475 = cos(qJ(1));
t491 = t472 * g(1) - g(2) * t475;
t448 = qJDD(1) * pkin(1) + t491;
t477 = qJD(1) ^ 2;
t486 = -g(1) * t475 - g(2) * t472;
t450 = -pkin(1) * t477 + t486;
t467 = sin(pkin(10));
t468 = cos(pkin(10));
t418 = t468 * t448 - t467 * t450;
t402 = -qJDD(1) * pkin(2) - t477 * pkin(7) - t418;
t493 = t471 * t494;
t453 = t474 * qJDD(1) - t493;
t386 = (-t452 - t492) * pkin(8) + (-t453 + t493) * pkin(3) + t402;
t419 = t467 * t448 + t468 * t450;
t403 = -pkin(2) * t477 + qJDD(1) * pkin(7) + t419;
t466 = -g(3) + qJDD(2);
t395 = t474 * t403 + t471 * t466;
t451 = (-pkin(3) * t474 - pkin(8) * t471) * qJD(1);
t476 = qJD(3) ^ 2;
t495 = qJD(1) * t474;
t392 = -pkin(3) * t476 + qJDD(3) * pkin(8) + t451 * t495 + t395;
t373 = t506 * t386 - t470 * t392;
t422 = pkin(4) * t446 - qJ(5) * t447;
t445 = qJDD(4) - t453;
t459 = qJD(4) - t495;
t458 = t459 ^ 2;
t371 = -t445 * pkin(4) - t458 * qJ(5) + t447 * t422 + qJDD(5) - t373;
t501 = t446 * t459;
t365 = (-t415 - t501) * pkin(9) + (t446 * t447 - t445) * pkin(5) + t371;
t374 = t470 * t386 + t506 * t392;
t507 = 2 * qJD(5);
t370 = -pkin(4) * t458 + t445 * qJ(5) - t446 * t422 + t459 * t507 + t374;
t414 = qJD(4) * t447 - t506 * qJDD(3) + t452 * t470;
t429 = -pkin(5) * t459 - pkin(9) * t447;
t444 = t446 ^ 2;
t366 = -pkin(5) * t444 + pkin(9) * t414 + t429 * t459 + t370;
t469 = sin(qJ(6));
t473 = cos(qJ(6));
t363 = t365 * t473 - t366 * t469;
t416 = t446 * t473 - t447 * t469;
t381 = qJD(6) * t416 + t414 * t469 + t415 * t473;
t417 = t446 * t469 + t447 * t473;
t393 = -mrSges(7,1) * t416 + mrSges(7,2) * t417;
t457 = qJD(6) - t459;
t396 = -mrSges(7,2) * t457 + mrSges(7,3) * t416;
t441 = qJDD(6) - t445;
t360 = m(7) * t363 + mrSges(7,1) * t441 - mrSges(7,3) * t381 - t393 * t417 + t396 * t457;
t364 = t365 * t469 + t366 * t473;
t380 = -qJD(6) * t417 + t414 * t473 - t415 * t469;
t397 = mrSges(7,1) * t457 - mrSges(7,3) * t417;
t361 = m(7) * t364 - mrSges(7,2) * t441 + mrSges(7,3) * t380 + t393 * t416 - t397 * t457;
t354 = t473 * t360 + t469 * t361;
t428 = -mrSges(6,2) * t446 + mrSges(6,3) * t459;
t482 = -m(6) * t371 + t445 * mrSges(6,1) + t459 * t428 - t354;
t353 = t415 * mrSges(6,2) + t447 * t423 - t482;
t384 = Ifges(7,4) * t417 + Ifges(7,2) * t416 + Ifges(7,6) * t457;
t385 = Ifges(7,1) * t417 + Ifges(7,4) * t416 + Ifges(7,5) * t457;
t481 = -mrSges(7,1) * t363 + mrSges(7,2) * t364 - Ifges(7,5) * t381 - Ifges(7,6) * t380 - Ifges(7,3) * t441 - t417 * t384 + t416 * t385;
t427 = -mrSges(6,1) * t459 + mrSges(6,2) * t447;
t488 = -t469 * t360 + t473 * t361;
t484 = m(6) * t370 + t445 * mrSges(6,3) + t459 * t427 + t488;
t498 = -t504 * t446 + t514 * t447 + t503 * t459;
t499 = t513 * t446 + t504 * t447 + t502 * t459;
t512 = -t502 * t414 + t503 * t415 + t511 * t445 + mrSges(5,1) * t373 - mrSges(6,1) * t371 - mrSges(5,2) * t374 + mrSges(6,3) * t370 - pkin(4) * t353 - pkin(5) * t354 + qJ(5) * (-t414 * mrSges(6,2) - t446 * t423 + t484) + t481 + t499 * t447 + t498 * t446;
t394 = -t471 * t403 + t474 * t466;
t483 = qJDD(3) * pkin(3) + t476 * pkin(8) - t451 * t496 + t394;
t508 = (-t415 + t501) * qJ(5) - t483;
t505 = -mrSges(5,3) - mrSges(6,2);
t500 = t502 * t446 - t503 * t447 - t511 * t459;
t497 = -mrSges(5,1) * t446 - mrSges(5,2) * t447 - t423;
t449 = (-mrSges(4,1) * t474 + mrSges(4,2) * t471) * qJD(1);
t454 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t496;
t426 = mrSges(5,1) * t459 - mrSges(5,3) * t447;
t350 = m(5) * t374 - t445 * mrSges(5,2) + t505 * t414 - t459 * t426 + t497 * t446 + t484;
t425 = -mrSges(5,2) * t459 - mrSges(5,3) * t446;
t351 = m(5) * t373 + t445 * mrSges(5,1) + t505 * t415 + t459 * t425 + t497 * t447 + t482;
t489 = t506 * t350 - t351 * t470;
t347 = m(4) * t395 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t453 - qJD(3) * t454 + t449 * t495 + t489;
t455 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t495;
t372 = -0.2e1 * qJD(5) * t447 + (t447 * t459 + t414) * pkin(4) + t508;
t368 = -t444 * pkin(9) + (-pkin(4) - pkin(5)) * t414 + (-pkin(4) * t459 + t429 + t507) * t447 - t508;
t485 = -m(7) * t368 + t380 * mrSges(7,1) - t381 * mrSges(7,2) + t416 * t396 - t417 * t397;
t358 = m(6) * t372 + t414 * mrSges(6,1) - t415 * mrSges(6,3) - t447 * t427 + t446 * t428 + t485;
t479 = m(5) * t483 - t414 * mrSges(5,1) - t415 * mrSges(5,2) - t446 * t425 - t447 * t426 - t358;
t357 = m(4) * t394 + qJDD(3) * mrSges(4,1) - t452 * mrSges(4,3) + qJD(3) * t455 - t449 * t496 + t479;
t490 = t474 * t347 - t357 * t471;
t348 = t470 * t350 + t506 * t351;
t480 = -m(4) * t402 + t453 * mrSges(4,1) - t452 * mrSges(4,2) - t454 * t496 + t455 * t495 - t348;
t437 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t471 + Ifges(4,4) * t474) * qJD(1);
t436 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t471 + Ifges(4,2) * t474) * qJD(1);
t383 = Ifges(7,5) * t417 + Ifges(7,6) * t416 + Ifges(7,3) * t457;
t356 = mrSges(7,2) * t368 - mrSges(7,3) * t363 + Ifges(7,1) * t381 + Ifges(7,4) * t380 + Ifges(7,5) * t441 + t383 * t416 - t384 * t457;
t355 = -mrSges(7,1) * t368 + mrSges(7,3) * t364 + Ifges(7,4) * t381 + Ifges(7,2) * t380 + Ifges(7,6) * t441 - t383 * t417 + t385 * t457;
t345 = -mrSges(5,2) * t483 + mrSges(6,2) * t371 - mrSges(5,3) * t373 - mrSges(6,3) * t372 - pkin(9) * t354 - qJ(5) * t358 - t469 * t355 + t473 * t356 - t504 * t414 + t514 * t415 + t503 * t445 + t500 * t446 - t499 * t459;
t344 = mrSges(5,1) * t483 - mrSges(6,1) * t372 + mrSges(6,2) * t370 + mrSges(5,3) * t374 - pkin(4) * t358 - pkin(5) * t485 - pkin(9) * t488 - t473 * t355 - t469 * t356 + t513 * t414 + t504 * t415 + t502 * t445 + t500 * t447 + t498 * t459;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t491 - mrSges(2,2) * t486 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t418 - mrSges(3,2) * t419 + t471 * (mrSges(4,2) * t402 - mrSges(4,3) * t394 + Ifges(4,1) * t452 + Ifges(4,4) * t453 + Ifges(4,5) * qJDD(3) - pkin(8) * t348 - qJD(3) * t436 - t470 * t344 + t506 * t345) + t474 * (-mrSges(4,1) * t402 + mrSges(4,3) * t395 + Ifges(4,4) * t452 + Ifges(4,2) * t453 + Ifges(4,6) * qJDD(3) - pkin(3) * t348 + qJD(3) * t437 - t512) + pkin(2) * t480 + pkin(7) * t490 + pkin(1) * (t467 * (m(3) * t419 - mrSges(3,1) * t477 - qJDD(1) * mrSges(3,2) + t490) + t468 * (m(3) * t418 + qJDD(1) * mrSges(3,1) - t477 * mrSges(3,2) + t480)); m(3) * t466 + t347 * t471 + t357 * t474; Ifges(4,5) * t452 + Ifges(4,6) * t453 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t394 - mrSges(4,2) * t395 + t470 * t345 + t506 * t344 + pkin(3) * t479 + pkin(8) * t489 + (t436 * t471 - t437 * t474) * qJD(1); t512; t353; -t481;];
tauJ  = t1;
