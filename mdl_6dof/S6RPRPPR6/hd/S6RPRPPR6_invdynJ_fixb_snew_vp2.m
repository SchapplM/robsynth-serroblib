% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-05 17:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:07:50
% EndTime: 2019-05-05 17:07:54
% DurationCPUTime: 3.13s
% Computational Cost: add. (28676->288), mult. (64946->367), div. (0->0), fcn. (43198->10), ass. (0->112)
t498 = -2 * qJD(4);
t478 = qJD(1) ^ 2;
t473 = sin(qJ(1));
t476 = cos(qJ(1));
t490 = t473 * g(1) - t476 * g(2);
t481 = -t478 * qJ(2) + qJDD(2) - t490;
t497 = -pkin(1) - pkin(7);
t442 = t497 * qJDD(1) + t481;
t472 = sin(qJ(3));
t475 = cos(qJ(3));
t430 = t472 * g(3) + t475 * t442;
t493 = qJD(1) * qJD(3);
t491 = t472 * t493;
t457 = t475 * qJDD(1) - t491;
t406 = (-t457 - t491) * qJ(4) + (-t472 * t475 * t478 + qJDD(3)) * pkin(3) + t430;
t431 = -t475 * g(3) + t472 * t442;
t456 = -t472 * qJDD(1) - t475 * t493;
t494 = t475 * qJD(1);
t459 = qJD(3) * pkin(3) - qJ(4) * t494;
t466 = t472 ^ 2;
t407 = -t466 * t478 * pkin(3) + t456 * qJ(4) - qJD(3) * t459 + t431;
t468 = sin(pkin(9));
t470 = cos(pkin(9));
t496 = qJD(1) * t472;
t449 = -t468 * t496 + t470 * t494;
t390 = t470 * t406 - t468 * t407 + t449 * t498;
t485 = -t476 * g(1) - t473 * g(2);
t482 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t485;
t448 = (t475 * t468 + t472 * t470) * qJD(1);
t391 = t468 * t406 + t470 * t407 + t448 * t498;
t425 = t448 * mrSges(5,1) + t449 * mrSges(5,2);
t428 = -t470 * t456 + t468 * t457;
t441 = qJD(3) * mrSges(5,1) - t449 * mrSges(5,3);
t424 = t448 * pkin(4) - t449 * qJ(5);
t477 = qJD(3) ^ 2;
t383 = -t477 * pkin(4) + qJDD(3) * qJ(5) - t448 * t424 + t391;
t411 = -t456 * pkin(3) + qJDD(4) + t459 * t494 + (-qJ(4) * t466 + t497) * t478 + t482;
t429 = t468 * t456 + t470 * t457;
t386 = (qJD(3) * t448 - t429) * qJ(5) + (qJD(3) * t449 + t428) * pkin(4) + t411;
t467 = sin(pkin(10));
t469 = cos(pkin(10));
t436 = t467 * qJD(3) + t469 * t449;
t378 = -0.2e1 * qJD(5) * t436 - t467 * t383 + t469 * t386;
t419 = t467 * qJDD(3) + t469 * t429;
t435 = t469 * qJD(3) - t467 * t449;
t376 = (t435 * t448 - t419) * pkin(8) + (t435 * t436 + t428) * pkin(5) + t378;
t379 = 0.2e1 * qJD(5) * t435 + t469 * t383 + t467 * t386;
t416 = t448 * pkin(5) - t436 * pkin(8);
t418 = t469 * qJDD(3) - t467 * t429;
t434 = t435 ^ 2;
t377 = -t434 * pkin(5) + t418 * pkin(8) - t448 * t416 + t379;
t471 = sin(qJ(6));
t474 = cos(qJ(6));
t374 = t474 * t376 - t471 * t377;
t409 = t474 * t435 - t471 * t436;
t389 = t409 * qJD(6) + t471 * t418 + t474 * t419;
t410 = t471 * t435 + t474 * t436;
t396 = -t409 * mrSges(7,1) + t410 * mrSges(7,2);
t446 = qJD(6) + t448;
t397 = -t446 * mrSges(7,2) + t409 * mrSges(7,3);
t427 = qJDD(6) + t428;
t371 = m(7) * t374 + t427 * mrSges(7,1) - t389 * mrSges(7,3) - t410 * t396 + t446 * t397;
t375 = t471 * t376 + t474 * t377;
t388 = -t410 * qJD(6) + t474 * t418 - t471 * t419;
t398 = t446 * mrSges(7,1) - t410 * mrSges(7,3);
t372 = m(7) * t375 - t427 * mrSges(7,2) + t388 * mrSges(7,3) + t409 * t396 - t446 * t398;
t363 = t474 * t371 + t471 * t372;
t412 = -t435 * mrSges(6,1) + t436 * mrSges(6,2);
t414 = -t448 * mrSges(6,2) + t435 * mrSges(6,3);
t361 = m(6) * t378 + t428 * mrSges(6,1) - t419 * mrSges(6,3) - t436 * t412 + t448 * t414 + t363;
t415 = t448 * mrSges(6,1) - t436 * mrSges(6,3);
t487 = -t471 * t371 + t474 * t372;
t362 = m(6) * t379 - t428 * mrSges(6,2) + t418 * mrSges(6,3) + t435 * t412 - t448 * t415 + t487;
t488 = -t467 * t361 + t469 * t362;
t355 = m(5) * t391 - qJDD(3) * mrSges(5,2) - t428 * mrSges(5,3) - qJD(3) * t441 - t448 * t425 + t488;
t382 = -qJDD(3) * pkin(4) - t477 * qJ(5) + t449 * t424 + qJDD(5) - t390;
t380 = -t418 * pkin(5) - t434 * pkin(8) + t436 * t416 + t382;
t480 = m(7) * t380 - t388 * mrSges(7,1) + t389 * mrSges(7,2) - t409 * t397 + t410 * t398;
t373 = m(6) * t382 - t418 * mrSges(6,1) + t419 * mrSges(6,2) - t435 * t414 + t436 * t415 + t480;
t440 = -qJD(3) * mrSges(5,2) - t448 * mrSges(5,3);
t367 = m(5) * t390 + qJDD(3) * mrSges(5,1) - t429 * mrSges(5,3) + qJD(3) * t440 - t449 * t425 - t373;
t350 = t468 * t355 + t470 * t367;
t357 = t469 * t361 + t467 * t362;
t489 = t470 * t355 - t468 * t367;
t455 = (t472 * mrSges(4,1) + t475 * mrSges(4,2)) * qJD(1);
t458 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t496;
t460 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t494;
t484 = t475 * (m(4) * t430 + qJDD(3) * mrSges(4,1) - t457 * mrSges(4,3) + qJD(3) * t458 - t455 * t494 + t350) + t472 * (m(4) * t431 - qJDD(3) * mrSges(4,2) + t456 * mrSges(4,3) - qJD(3) * t460 - t455 * t496 + t489);
t356 = m(5) * t411 + t428 * mrSges(5,1) + t429 * mrSges(5,2) + t448 * t440 + t449 * t441 + t357;
t393 = Ifges(7,4) * t410 + Ifges(7,2) * t409 + Ifges(7,6) * t446;
t394 = Ifges(7,1) * t410 + Ifges(7,4) * t409 + Ifges(7,5) * t446;
t479 = mrSges(7,1) * t374 - mrSges(7,2) * t375 + Ifges(7,5) * t389 + Ifges(7,6) * t388 + Ifges(7,3) * t427 + t410 * t393 - t409 * t394;
t452 = Ifges(4,5) * qJD(3) + (t475 * Ifges(4,1) - t472 * Ifges(4,4)) * qJD(1);
t451 = Ifges(4,6) * qJD(3) + (t475 * Ifges(4,4) - t472 * Ifges(4,2)) * qJD(1);
t447 = -qJDD(1) * pkin(1) + t481;
t443 = t478 * pkin(1) - t482;
t439 = t497 * t478 + t482;
t422 = Ifges(5,1) * t449 - Ifges(5,4) * t448 + Ifges(5,5) * qJD(3);
t421 = Ifges(5,4) * t449 - Ifges(5,2) * t448 + Ifges(5,6) * qJD(3);
t420 = Ifges(5,5) * t449 - Ifges(5,6) * t448 + Ifges(5,3) * qJD(3);
t401 = Ifges(6,1) * t436 + Ifges(6,4) * t435 + Ifges(6,5) * t448;
t400 = Ifges(6,4) * t436 + Ifges(6,2) * t435 + Ifges(6,6) * t448;
t399 = Ifges(6,5) * t436 + Ifges(6,6) * t435 + Ifges(6,3) * t448;
t392 = Ifges(7,5) * t410 + Ifges(7,6) * t409 + Ifges(7,3) * t446;
t365 = mrSges(7,2) * t380 - mrSges(7,3) * t374 + Ifges(7,1) * t389 + Ifges(7,4) * t388 + Ifges(7,5) * t427 + t409 * t392 - t446 * t393;
t364 = -mrSges(7,1) * t380 + mrSges(7,3) * t375 + Ifges(7,4) * t389 + Ifges(7,2) * t388 + Ifges(7,6) * t427 - t410 * t392 + t446 * t394;
t352 = mrSges(6,2) * t382 - mrSges(6,3) * t378 + Ifges(6,1) * t419 + Ifges(6,4) * t418 + Ifges(6,5) * t428 - pkin(8) * t363 - t471 * t364 + t474 * t365 + t435 * t399 - t448 * t400;
t351 = -mrSges(6,1) * t382 + mrSges(6,3) * t379 + Ifges(6,4) * t419 + Ifges(6,2) * t418 + Ifges(6,6) * t428 - pkin(5) * t480 + pkin(8) * t487 + t474 * t364 + t471 * t365 - t436 * t399 + t448 * t401;
t347 = -t479 + (-Ifges(5,2) - Ifges(6,3)) * t428 + Ifges(5,6) * qJDD(3) - t449 * t420 + t435 * t401 - t436 * t400 + Ifges(5,4) * t429 - Ifges(6,6) * t418 - Ifges(6,5) * t419 + qJD(3) * t422 - mrSges(5,1) * t411 + mrSges(5,3) * t391 + mrSges(6,2) * t379 - mrSges(6,1) * t378 - pkin(5) * t363 - pkin(4) * t357;
t346 = m(3) * t447 + qJDD(1) * mrSges(3,2) - t478 * mrSges(3,3) + t484;
t345 = mrSges(5,2) * t411 - mrSges(5,3) * t390 + Ifges(5,1) * t429 - Ifges(5,4) * t428 + Ifges(5,5) * qJDD(3) - qJ(5) * t357 - qJD(3) * t421 - t467 * t351 + t469 * t352 - t448 * t420;
t1 = [mrSges(2,1) * t490 - mrSges(2,2) * t485 + mrSges(3,2) * t447 - mrSges(3,3) * t443 + t475 * (mrSges(4,2) * t439 - mrSges(4,3) * t430 + Ifges(4,1) * t457 + Ifges(4,4) * t456 + Ifges(4,5) * qJDD(3) - qJ(4) * t350 - qJD(3) * t451 + t470 * t345 - t468 * t347) - t472 * (-mrSges(4,1) * t439 + mrSges(4,3) * t431 + Ifges(4,4) * t457 + Ifges(4,2) * t456 + Ifges(4,6) * qJDD(3) - pkin(3) * t356 + qJ(4) * t489 + qJD(3) * t452 + t468 * t345 + t470 * t347) - pkin(7) * t484 - pkin(1) * t346 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t443 + m(4) * t439 - t456 * mrSges(4,1) + t478 * mrSges(3,2) + t457 * mrSges(4,2) + t356 + qJDD(1) * mrSges(3,3) + (t458 * t472 + t460 * t475) * qJD(1)) * qJ(2); t346; Ifges(4,5) * t457 + Ifges(4,6) * t456 + mrSges(4,1) * t430 - mrSges(4,2) * t431 + Ifges(5,5) * t429 - Ifges(5,6) * t428 + t449 * t421 + t448 * t422 + mrSges(5,1) * t390 - mrSges(5,2) * t391 + t467 * t352 + t469 * t351 - pkin(4) * t373 + qJ(5) * t488 + pkin(3) * t350 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t475 * t451 + t472 * t452) * qJD(1); t356; t373; t479;];
tauJ  = t1;
