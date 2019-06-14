% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR2
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
% Datum: 2019-05-05 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:29:08
% EndTime: 2019-05-05 04:29:12
% DurationCPUTime: 3.72s
% Computational Cost: add. (36378->291), mult. (78765->375), div. (0->0), fcn. (57344->14), ass. (0->121)
t485 = sin(pkin(11));
t488 = cos(pkin(11));
t475 = t485 * g(1) - t488 * g(2);
t483 = -g(3) + qJDD(1);
t486 = sin(pkin(6));
t489 = cos(pkin(6));
t519 = t475 * t489 + t483 * t486;
t476 = -t488 * g(1) - t485 * g(2);
t493 = sin(qJ(2));
t497 = cos(qJ(2));
t438 = t497 * t476 + t519 * t493;
t499 = qJD(2) ^ 2;
t433 = -t499 * pkin(2) + qJDD(2) * pkin(8) + t438;
t455 = -t486 * t475 + t489 * t483;
t492 = sin(qJ(3));
t496 = cos(qJ(3));
t415 = -t492 * t433 + t496 * t455;
t513 = qJD(2) * qJD(3);
t512 = t496 * t513;
t473 = t492 * qJDD(2) + t512;
t410 = (-t473 + t512) * qJ(4) + (t492 * t496 * t499 + qJDD(3)) * pkin(3) + t415;
t416 = t496 * t433 + t492 * t455;
t474 = t496 * qJDD(2) - t492 * t513;
t514 = t492 * qJD(2);
t477 = qJD(3) * pkin(3) - qJ(4) * t514;
t482 = t496 ^ 2;
t411 = -t482 * t499 * pkin(3) + t474 * qJ(4) - qJD(3) * t477 + t416;
t484 = sin(pkin(12));
t487 = cos(pkin(12));
t461 = (t496 * t484 + t492 * t487) * qJD(2);
t391 = -0.2e1 * qJD(4) * t461 + t487 * t410 - t484 * t411;
t437 = -t493 * t476 + t497 * t519;
t516 = qJD(2) * t496;
t460 = -t484 * t514 + t487 * t516;
t392 = 0.2e1 * qJD(4) * t460 + t484 * t410 + t487 * t411;
t440 = -t460 * mrSges(5,1) + t461 * mrSges(5,2);
t446 = -t484 * t473 + t487 * t474;
t454 = qJD(3) * mrSges(5,1) - t461 * mrSges(5,3);
t442 = -t460 * pkin(4) - t461 * pkin(9);
t498 = qJD(3) ^ 2;
t390 = -t498 * pkin(4) + qJDD(3) * pkin(9) + t460 * t442 + t392;
t503 = -qJDD(2) * pkin(2) - t437;
t412 = -t474 * pkin(3) + qJDD(4) + t477 * t514 + (-qJ(4) * t482 - pkin(8)) * t499 + t503;
t447 = t487 * t473 + t484 * t474;
t401 = (-qJD(3) * t460 - t447) * pkin(9) + (qJD(3) * t461 - t446) * pkin(4) + t412;
t491 = sin(qJ(5));
t495 = cos(qJ(5));
t385 = -t491 * t390 + t495 * t401;
t449 = t495 * qJD(3) - t491 * t461;
t423 = t449 * qJD(5) + t491 * qJDD(3) + t495 * t447;
t445 = qJDD(5) - t446;
t450 = t491 * qJD(3) + t495 * t461;
t459 = qJD(5) - t460;
t383 = (t449 * t459 - t423) * pkin(10) + (t449 * t450 + t445) * pkin(5) + t385;
t386 = t495 * t390 + t491 * t401;
t422 = -t450 * qJD(5) + t495 * qJDD(3) - t491 * t447;
t431 = t459 * pkin(5) - t450 * pkin(10);
t448 = t449 ^ 2;
t384 = -t448 * pkin(5) + t422 * pkin(10) - t459 * t431 + t386;
t490 = sin(qJ(6));
t494 = cos(qJ(6));
t381 = t494 * t383 - t490 * t384;
t424 = t494 * t449 - t490 * t450;
t397 = t424 * qJD(6) + t490 * t422 + t494 * t423;
t425 = t490 * t449 + t494 * t450;
t406 = -t424 * mrSges(7,1) + t425 * mrSges(7,2);
t456 = qJD(6) + t459;
t413 = -t456 * mrSges(7,2) + t424 * mrSges(7,3);
t443 = qJDD(6) + t445;
t378 = m(7) * t381 + t443 * mrSges(7,1) - t397 * mrSges(7,3) - t425 * t406 + t456 * t413;
t382 = t490 * t383 + t494 * t384;
t396 = -t425 * qJD(6) + t494 * t422 - t490 * t423;
t414 = t456 * mrSges(7,1) - t425 * mrSges(7,3);
t379 = m(7) * t382 - t443 * mrSges(7,2) + t396 * mrSges(7,3) + t424 * t406 - t456 * t414;
t370 = t494 * t378 + t490 * t379;
t426 = -t449 * mrSges(6,1) + t450 * mrSges(6,2);
t429 = -t459 * mrSges(6,2) + t449 * mrSges(6,3);
t368 = m(6) * t385 + t445 * mrSges(6,1) - t423 * mrSges(6,3) - t450 * t426 + t459 * t429 + t370;
t430 = t459 * mrSges(6,1) - t450 * mrSges(6,3);
t508 = -t490 * t378 + t494 * t379;
t369 = m(6) * t386 - t445 * mrSges(6,2) + t422 * mrSges(6,3) + t449 * t426 - t459 * t430 + t508;
t509 = -t491 * t368 + t495 * t369;
t362 = m(5) * t392 - qJDD(3) * mrSges(5,2) + t446 * mrSges(5,3) - qJD(3) * t454 + t460 * t440 + t509;
t453 = -qJD(3) * mrSges(5,2) + t460 * mrSges(5,3);
t389 = -qJDD(3) * pkin(4) - t498 * pkin(9) + t461 * t442 - t391;
t387 = -t422 * pkin(5) - t448 * pkin(10) + t450 * t431 + t389;
t505 = m(7) * t387 - t396 * mrSges(7,1) + t397 * mrSges(7,2) - t424 * t413 + t425 * t414;
t502 = -m(6) * t389 + t422 * mrSges(6,1) - t423 * mrSges(6,2) + t449 * t429 - t450 * t430 - t505;
t374 = m(5) * t391 + qJDD(3) * mrSges(5,1) - t447 * mrSges(5,3) + qJD(3) * t453 - t461 * t440 + t502;
t357 = t484 * t362 + t487 * t374;
t364 = t495 * t368 + t491 * t369;
t472 = (-t496 * mrSges(4,1) + t492 * mrSges(4,2)) * qJD(2);
t479 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t516;
t355 = m(4) * t415 + qJDD(3) * mrSges(4,1) - t473 * mrSges(4,3) + qJD(3) * t479 - t472 * t514 + t357;
t478 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t514;
t510 = t487 * t362 - t484 * t374;
t356 = m(4) * t416 - qJDD(3) * mrSges(4,2) + t474 * mrSges(4,3) - qJD(3) * t478 + t472 * t516 + t510;
t511 = -t492 * t355 + t496 * t356;
t403 = Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t456;
t404 = Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t456;
t504 = -mrSges(7,1) * t381 + mrSges(7,2) * t382 - Ifges(7,5) * t397 - Ifges(7,6) * t396 - Ifges(7,3) * t443 - t425 * t403 + t424 * t404;
t363 = m(5) * t412 - t446 * mrSges(5,1) + t447 * mrSges(5,2) - t460 * t453 + t461 * t454 + t364;
t432 = -t499 * pkin(8) + t503;
t501 = -m(4) * t432 + t474 * mrSges(4,1) - t473 * mrSges(4,2) - t478 * t514 + t479 * t516 - t363;
t418 = Ifges(6,4) * t450 + Ifges(6,2) * t449 + Ifges(6,6) * t459;
t419 = Ifges(6,1) * t450 + Ifges(6,4) * t449 + Ifges(6,5) * t459;
t500 = mrSges(6,1) * t385 - mrSges(6,2) * t386 + Ifges(6,5) * t423 + Ifges(6,6) * t422 + Ifges(6,3) * t445 + pkin(5) * t370 + t450 * t418 - t449 * t419 - t504;
t465 = Ifges(4,5) * qJD(3) + (t492 * Ifges(4,1) + t496 * Ifges(4,4)) * qJD(2);
t464 = Ifges(4,6) * qJD(3) + (t492 * Ifges(4,4) + t496 * Ifges(4,2)) * qJD(2);
t436 = Ifges(5,1) * t461 + Ifges(5,4) * t460 + Ifges(5,5) * qJD(3);
t435 = Ifges(5,4) * t461 + Ifges(5,2) * t460 + Ifges(5,6) * qJD(3);
t434 = Ifges(5,5) * t461 + Ifges(5,6) * t460 + Ifges(5,3) * qJD(3);
t417 = Ifges(6,5) * t450 + Ifges(6,6) * t449 + Ifges(6,3) * t459;
t402 = Ifges(7,5) * t425 + Ifges(7,6) * t424 + Ifges(7,3) * t456;
t372 = mrSges(7,2) * t387 - mrSges(7,3) * t381 + Ifges(7,1) * t397 + Ifges(7,4) * t396 + Ifges(7,5) * t443 + t424 * t402 - t456 * t403;
t371 = -mrSges(7,1) * t387 + mrSges(7,3) * t382 + Ifges(7,4) * t397 + Ifges(7,2) * t396 + Ifges(7,6) * t443 - t425 * t402 + t456 * t404;
t359 = mrSges(6,2) * t389 - mrSges(6,3) * t385 + Ifges(6,1) * t423 + Ifges(6,4) * t422 + Ifges(6,5) * t445 - pkin(10) * t370 - t490 * t371 + t494 * t372 + t449 * t417 - t459 * t418;
t358 = -mrSges(6,1) * t389 + mrSges(6,3) * t386 + Ifges(6,4) * t423 + Ifges(6,2) * t422 + Ifges(6,6) * t445 - pkin(5) * t505 + pkin(10) * t508 + t494 * t371 + t490 * t372 - t450 * t417 + t459 * t419;
t353 = -mrSges(5,1) * t412 + mrSges(5,3) * t392 + Ifges(5,4) * t447 + Ifges(5,2) * t446 + Ifges(5,6) * qJDD(3) - pkin(4) * t364 + qJD(3) * t436 - t461 * t434 - t500;
t352 = mrSges(5,2) * t412 - mrSges(5,3) * t391 + Ifges(5,1) * t447 + Ifges(5,4) * t446 + Ifges(5,5) * qJDD(3) - pkin(9) * t364 - qJD(3) * t435 - t491 * t358 + t495 * t359 + t460 * t434;
t1 = [m(2) * t483 + t489 * (m(3) * t455 + t496 * t355 + t492 * t356) + (t493 * (m(3) * t438 - t499 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t511) + t497 * (m(3) * t437 + qJDD(2) * mrSges(3,1) - t499 * mrSges(3,2) + t501)) * t486; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t437 - mrSges(3,2) * t438 + t492 * (mrSges(4,2) * t432 - mrSges(4,3) * t415 + Ifges(4,1) * t473 + Ifges(4,4) * t474 + Ifges(4,5) * qJDD(3) - qJ(4) * t357 - qJD(3) * t464 + t487 * t352 - t484 * t353) + t496 * (-mrSges(4,1) * t432 + mrSges(4,3) * t416 + Ifges(4,4) * t473 + Ifges(4,2) * t474 + Ifges(4,6) * qJDD(3) - pkin(3) * t363 + qJ(4) * t510 + qJD(3) * t465 + t484 * t352 + t487 * t353) + pkin(2) * t501 + pkin(8) * t511; Ifges(4,5) * t473 + Ifges(4,6) * t474 + mrSges(4,1) * t415 - mrSges(4,2) * t416 + Ifges(5,5) * t447 + Ifges(5,6) * t446 + t461 * t435 - t460 * t436 + mrSges(5,1) * t391 - mrSges(5,2) * t392 + t491 * t359 + t495 * t358 + pkin(4) * t502 + pkin(9) * t509 + pkin(3) * t357 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t492 * t464 - t496 * t465) * qJD(2); t363; t500; -t504;];
tauJ  = t1;
