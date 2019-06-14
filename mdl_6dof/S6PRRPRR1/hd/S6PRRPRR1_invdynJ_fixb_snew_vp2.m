% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:15:12
% EndTime: 2019-05-05 04:15:15
% DurationCPUTime: 3.81s
% Computational Cost: add. (38318->290), mult. (84182->376), div. (0->0), fcn. (62739->14), ass. (0->121)
t485 = sin(pkin(11));
t488 = cos(pkin(11));
t472 = t485 * g(1) - t488 * g(2);
t483 = -g(3) + qJDD(1);
t486 = sin(pkin(6));
t489 = cos(pkin(6));
t516 = t472 * t489 + t483 * t486;
t473 = -t488 * g(1) - t485 * g(2);
t493 = sin(qJ(2));
t497 = cos(qJ(2));
t441 = -t493 * t473 + t516 * t497;
t442 = t497 * t473 + t516 * t493;
t498 = qJD(2) ^ 2;
t434 = -t498 * pkin(2) + qJDD(2) * pkin(8) + t442;
t454 = -t486 * t472 + t489 * t483;
t492 = sin(qJ(3));
t496 = cos(qJ(3));
t426 = -t492 * t434 + t496 * t454;
t511 = qJD(2) * qJD(3);
t510 = t496 * t511;
t470 = t492 * qJDD(2) + t510;
t413 = (-t470 + t510) * qJ(4) + (t492 * t496 * t498 + qJDD(3)) * pkin(3) + t426;
t427 = t496 * t434 + t492 * t454;
t471 = t496 * qJDD(2) - t492 * t511;
t512 = t492 * qJD(2);
t474 = qJD(3) * pkin(3) - qJ(4) * t512;
t482 = t496 ^ 2;
t415 = -t482 * t498 * pkin(3) + t471 * qJ(4) - qJD(3) * t474 + t427;
t484 = sin(pkin(12));
t487 = cos(pkin(12));
t459 = (t496 * t484 + t492 * t487) * qJD(2);
t391 = -0.2e1 * qJD(4) * t459 + t487 * t413 - t484 * t415;
t448 = t487 * t470 + t484 * t471;
t458 = (-t492 * t484 + t496 * t487) * qJD(2);
t388 = (qJD(3) * t458 - t448) * pkin(9) + (t458 * t459 + qJDD(3)) * pkin(4) + t391;
t392 = 0.2e1 * qJD(4) * t458 + t484 * t413 + t487 * t415;
t447 = -t484 * t470 + t487 * t471;
t453 = qJD(3) * pkin(4) - pkin(9) * t459;
t457 = t458 ^ 2;
t390 = -pkin(4) * t457 + pkin(9) * t447 - qJD(3) * t453 + t392;
t491 = sin(qJ(5));
t495 = cos(qJ(5));
t385 = t491 * t388 + t495 * t390;
t440 = t491 * t458 + t495 * t459;
t408 = -t440 * qJD(5) + t495 * t447 - t491 * t448;
t439 = t495 * t458 - t491 * t459;
t423 = -mrSges(6,1) * t439 + mrSges(6,2) * t440;
t481 = qJD(3) + qJD(5);
t432 = t481 * mrSges(6,1) - t440 * mrSges(6,3);
t480 = qJDD(3) + qJDD(5);
t424 = -pkin(5) * t439 - pkin(10) * t440;
t479 = t481 ^ 2;
t382 = -t479 * pkin(5) + t480 * pkin(10) + t439 * t424 + t385;
t502 = -qJDD(2) * pkin(2) - t441;
t425 = -t471 * pkin(3) + qJDD(4) + t474 * t512 + (-qJ(4) * t482 - pkin(8)) * t498 + t502;
t397 = -t447 * pkin(4) - t457 * pkin(9) + t459 * t453 + t425;
t409 = t439 * qJD(5) + t491 * t447 + t495 * t448;
t386 = (-t439 * t481 - t409) * pkin(10) + (t440 * t481 - t408) * pkin(5) + t397;
t490 = sin(qJ(6));
t494 = cos(qJ(6));
t379 = -t490 * t382 + t494 * t386;
t428 = -t490 * t440 + t494 * t481;
t395 = t428 * qJD(6) + t494 * t409 + t490 * t480;
t407 = qJDD(6) - t408;
t429 = t494 * t440 + t490 * t481;
t414 = -mrSges(7,1) * t428 + mrSges(7,2) * t429;
t435 = qJD(6) - t439;
t416 = -mrSges(7,2) * t435 + mrSges(7,3) * t428;
t376 = m(7) * t379 + mrSges(7,1) * t407 - mrSges(7,3) * t395 - t414 * t429 + t416 * t435;
t380 = t494 * t382 + t490 * t386;
t394 = -t429 * qJD(6) - t490 * t409 + t494 * t480;
t417 = mrSges(7,1) * t435 - mrSges(7,3) * t429;
t377 = m(7) * t380 - mrSges(7,2) * t407 + mrSges(7,3) * t394 + t414 * t428 - t417 * t435;
t506 = -t490 * t376 + t494 * t377;
t363 = m(6) * t385 - t480 * mrSges(6,2) + t408 * mrSges(6,3) + t439 * t423 - t481 * t432 + t506;
t384 = t495 * t388 - t491 * t390;
t431 = -t481 * mrSges(6,2) + t439 * mrSges(6,3);
t381 = -t480 * pkin(5) - t479 * pkin(10) + t440 * t424 - t384;
t503 = -m(7) * t381 + t394 * mrSges(7,1) - t395 * mrSges(7,2) + t428 * t416 - t429 * t417;
t372 = m(6) * t384 + t480 * mrSges(6,1) - t409 * mrSges(6,3) - t440 * t423 + t481 * t431 + t503;
t360 = t491 * t363 + t495 * t372;
t445 = -mrSges(5,1) * t458 + mrSges(5,2) * t459;
t451 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t458;
t358 = m(5) * t391 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t448 + qJD(3) * t451 - t445 * t459 + t360;
t452 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t459;
t507 = t495 * t363 - t491 * t372;
t359 = m(5) * t392 - qJDD(3) * mrSges(5,2) + t447 * mrSges(5,3) - qJD(3) * t452 + t458 * t445 + t507;
t352 = t487 * t358 + t484 * t359;
t366 = t494 * t376 + t490 * t377;
t513 = qJD(2) * t496;
t469 = (-t496 * mrSges(4,1) + t492 * mrSges(4,2)) * qJD(2);
t476 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t513;
t350 = m(4) * t426 + qJDD(3) * mrSges(4,1) - t470 * mrSges(4,3) + qJD(3) * t476 - t469 * t512 + t352;
t475 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t512;
t508 = -t484 * t358 + t487 * t359;
t351 = m(4) * t427 - qJDD(3) * mrSges(4,2) + t471 * mrSges(4,3) - qJD(3) * t475 + t469 * t513 + t508;
t509 = -t492 * t350 + t496 * t351;
t504 = m(6) * t397 - t408 * mrSges(6,1) + t409 * mrSges(6,2) - t439 * t431 + t440 * t432 + t366;
t398 = Ifges(7,5) * t429 + Ifges(7,6) * t428 + Ifges(7,3) * t435;
t400 = Ifges(7,1) * t429 + Ifges(7,4) * t428 + Ifges(7,5) * t435;
t369 = -mrSges(7,1) * t381 + mrSges(7,3) * t380 + Ifges(7,4) * t395 + Ifges(7,2) * t394 + Ifges(7,6) * t407 - t398 * t429 + t400 * t435;
t399 = Ifges(7,4) * t429 + Ifges(7,2) * t428 + Ifges(7,6) * t435;
t370 = mrSges(7,2) * t381 - mrSges(7,3) * t379 + Ifges(7,1) * t395 + Ifges(7,4) * t394 + Ifges(7,5) * t407 + t398 * t428 - t399 * t435;
t419 = Ifges(6,4) * t440 + Ifges(6,2) * t439 + Ifges(6,6) * t481;
t420 = Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * t481;
t501 = mrSges(6,1) * t384 - mrSges(6,2) * t385 + Ifges(6,5) * t409 + Ifges(6,6) * t408 + Ifges(6,3) * t480 + pkin(5) * t503 + pkin(10) * t506 + t494 * t369 + t490 * t370 + t440 * t419 - t439 * t420;
t364 = m(5) * t425 - t447 * mrSges(5,1) + t448 * mrSges(5,2) - t458 * t451 + t459 * t452 + t504;
t500 = mrSges(7,1) * t379 - mrSges(7,2) * t380 + Ifges(7,5) * t395 + Ifges(7,6) * t394 + Ifges(7,3) * t407 + t429 * t399 - t428 * t400;
t433 = -t498 * pkin(8) + t502;
t499 = -m(4) * t433 + t471 * mrSges(4,1) - t470 * mrSges(4,2) - t475 * t512 + t476 * t513 - t364;
t463 = Ifges(4,5) * qJD(3) + (t492 * Ifges(4,1) + t496 * Ifges(4,4)) * qJD(2);
t462 = Ifges(4,6) * qJD(3) + (t492 * Ifges(4,4) + t496 * Ifges(4,2)) * qJD(2);
t438 = Ifges(5,1) * t459 + Ifges(5,4) * t458 + Ifges(5,5) * qJD(3);
t437 = Ifges(5,4) * t459 + Ifges(5,2) * t458 + Ifges(5,6) * qJD(3);
t436 = Ifges(5,5) * t459 + Ifges(5,6) * t458 + Ifges(5,3) * qJD(3);
t418 = Ifges(6,5) * t440 + Ifges(6,6) * t439 + Ifges(6,3) * t481;
t354 = -mrSges(6,1) * t397 + mrSges(6,3) * t385 + Ifges(6,4) * t409 + Ifges(6,2) * t408 + Ifges(6,6) * t480 - pkin(5) * t366 - t440 * t418 + t481 * t420 - t500;
t353 = mrSges(6,2) * t397 - mrSges(6,3) * t384 + Ifges(6,1) * t409 + Ifges(6,4) * t408 + Ifges(6,5) * t480 - pkin(10) * t366 - t490 * t369 + t494 * t370 + t439 * t418 - t481 * t419;
t348 = mrSges(5,2) * t425 - mrSges(5,3) * t391 + Ifges(5,1) * t448 + Ifges(5,4) * t447 + Ifges(5,5) * qJDD(3) - pkin(9) * t360 - qJD(3) * t437 + t495 * t353 - t491 * t354 + t458 * t436;
t347 = -mrSges(5,1) * t425 + mrSges(5,3) * t392 + Ifges(5,4) * t448 + Ifges(5,2) * t447 + Ifges(5,6) * qJDD(3) - pkin(4) * t504 + pkin(9) * t507 + qJD(3) * t438 + t491 * t353 + t495 * t354 - t459 * t436;
t1 = [m(2) * t483 + t489 * (m(3) * t454 + t496 * t350 + t492 * t351) + (t493 * (m(3) * t442 - t498 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t509) + t497 * (m(3) * t441 + qJDD(2) * mrSges(3,1) - t498 * mrSges(3,2) + t499)) * t486; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t441 - mrSges(3,2) * t442 + t492 * (mrSges(4,2) * t433 - mrSges(4,3) * t426 + Ifges(4,1) * t470 + Ifges(4,4) * t471 + Ifges(4,5) * qJDD(3) - qJ(4) * t352 - qJD(3) * t462 - t484 * t347 + t487 * t348) + t496 * (-mrSges(4,1) * t433 + mrSges(4,3) * t427 + Ifges(4,4) * t470 + Ifges(4,2) * t471 + Ifges(4,6) * qJDD(3) - pkin(3) * t364 + qJ(4) * t508 + qJD(3) * t463 + t487 * t347 + t484 * t348) + pkin(2) * t499 + pkin(8) * t509; mrSges(5,1) * t391 - mrSges(5,2) * t392 + mrSges(4,1) * t426 - mrSges(4,2) * t427 + pkin(4) * t360 + pkin(3) * t352 + Ifges(4,5) * t470 + Ifges(4,6) * t471 - t458 * t438 + t459 * t437 + Ifges(5,6) * t447 + Ifges(5,5) * t448 + t501 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t492 * t462 - t496 * t463) * qJD(2); t364; t501; t500;];
tauJ  = t1;
