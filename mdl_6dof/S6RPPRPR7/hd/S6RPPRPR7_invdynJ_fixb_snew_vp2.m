% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-05 14:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:31:37
% EndTime: 2019-05-05 14:31:40
% DurationCPUTime: 2.82s
% Computational Cost: add. (25120->257), mult. (59079->326), div. (0->0), fcn. (41640->10), ass. (0->112)
t473 = qJD(1) ^ 2;
t468 = sin(qJ(1));
t471 = cos(qJ(1));
t490 = t468 * g(1) - t471 * g(2);
t479 = -t473 * qJ(2) + qJDD(2) - t490;
t500 = -pkin(1) - qJ(3);
t504 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t500 + t479;
t465 = cos(pkin(9));
t503 = t465 ^ 2;
t485 = -t471 * g(1) - t468 * g(2);
t502 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t485;
t501 = pkin(3) * t473;
t499 = t465 * mrSges(4,2);
t463 = sin(pkin(9));
t428 = t463 * g(3) + t465 * t504;
t410 = (-pkin(7) * qJDD(1) - t463 * t501) * t465 + t428;
t429 = -t465 * g(3) + t463 * t504;
t459 = t463 ^ 2;
t492 = t463 * qJDD(1);
t411 = -pkin(7) * t492 - t459 * t501 + t429;
t467 = sin(qJ(4));
t470 = cos(qJ(4));
t396 = t467 * t410 + t470 * t411;
t483 = t463 * t470 + t465 * t467;
t447 = t483 * qJD(1);
t482 = -t463 * t467 + t465 * t470;
t448 = t482 * qJD(1);
t425 = t447 * mrSges(5,1) + t448 * mrSges(5,2);
t494 = t448 * qJD(4);
t430 = qJDD(1) * t483 + t494;
t442 = qJD(4) * mrSges(5,1) - t448 * mrSges(5,3);
t424 = t447 * pkin(4) - t448 * qJ(5);
t472 = qJD(4) ^ 2;
t383 = -t472 * pkin(4) + qJDD(4) * qJ(5) - t447 * t424 + t396;
t478 = qJDD(3) + t502;
t497 = -t459 - t503;
t419 = pkin(3) * t492 + (pkin(7) * t497 + t500) * t473 + t478;
t495 = t447 * qJD(4);
t431 = qJDD(1) * t482 - t495;
t389 = (-t431 + t495) * qJ(5) + (t430 + t494) * pkin(4) + t419;
t462 = sin(pkin(10));
t464 = cos(pkin(10));
t436 = t462 * qJD(4) + t464 * t448;
t378 = -0.2e1 * qJD(5) * t436 - t462 * t383 + t464 * t389;
t417 = t462 * qJDD(4) + t464 * t431;
t435 = t464 * qJD(4) - t462 * t448;
t376 = (t435 * t447 - t417) * pkin(8) + (t435 * t436 + t430) * pkin(5) + t378;
t379 = 0.2e1 * qJD(5) * t435 + t464 * t383 + t462 * t389;
t415 = t447 * pkin(5) - t436 * pkin(8);
t416 = t464 * qJDD(4) - t462 * t431;
t434 = t435 ^ 2;
t377 = -t434 * pkin(5) + t416 * pkin(8) - t447 * t415 + t379;
t466 = sin(qJ(6));
t469 = cos(qJ(6));
t374 = t469 * t376 - t466 * t377;
t403 = t469 * t435 - t466 * t436;
t388 = t403 * qJD(6) + t466 * t416 + t469 * t417;
t404 = t466 * t435 + t469 * t436;
t394 = -t403 * mrSges(7,1) + t404 * mrSges(7,2);
t445 = qJD(6) + t447;
t397 = -t445 * mrSges(7,2) + t403 * mrSges(7,3);
t427 = qJDD(6) + t430;
t371 = m(7) * t374 + t427 * mrSges(7,1) - t388 * mrSges(7,3) - t404 * t394 + t445 * t397;
t375 = t466 * t376 + t469 * t377;
t387 = -t404 * qJD(6) + t469 * t416 - t466 * t417;
t398 = t445 * mrSges(7,1) - t404 * mrSges(7,3);
t372 = m(7) * t375 - t427 * mrSges(7,2) + t387 * mrSges(7,3) + t403 * t394 - t445 * t398;
t363 = t469 * t371 + t466 * t372;
t405 = -t435 * mrSges(6,1) + t436 * mrSges(6,2);
t413 = -t447 * mrSges(6,2) + t435 * mrSges(6,3);
t361 = m(6) * t378 + t430 * mrSges(6,1) - t417 * mrSges(6,3) - t436 * t405 + t447 * t413 + t363;
t414 = t447 * mrSges(6,1) - t436 * mrSges(6,3);
t486 = -t466 * t371 + t469 * t372;
t362 = m(6) * t379 - t430 * mrSges(6,2) + t416 * mrSges(6,3) + t435 * t405 - t447 * t414 + t486;
t487 = -t462 * t361 + t464 * t362;
t356 = m(5) * t396 - qJDD(4) * mrSges(5,2) - t430 * mrSges(5,3) - qJD(4) * t442 - t447 * t425 + t487;
t395 = t470 * t410 - t467 * t411;
t382 = -qJDD(4) * pkin(4) - t472 * qJ(5) + t448 * t424 + qJDD(5) - t395;
t380 = -t416 * pkin(5) - t434 * pkin(8) + t436 * t415 + t382;
t477 = m(7) * t380 - t387 * mrSges(7,1) + t388 * mrSges(7,2) - t403 * t397 + t404 * t398;
t373 = m(6) * t382 - t416 * mrSges(6,1) + t417 * mrSges(6,2) - t435 * t413 + t436 * t414 + t477;
t441 = -qJD(4) * mrSges(5,2) - t447 * mrSges(5,3);
t367 = m(5) * t395 + qJDD(4) * mrSges(5,1) - t431 * mrSges(5,3) + qJD(4) * t441 - t448 * t425 - t373;
t498 = t467 * t356 + t470 * t367;
t357 = t464 * t361 + t462 * t362;
t489 = t497 * mrSges(4,3);
t488 = t470 * t356 - t467 * t367;
t481 = -qJDD(1) * mrSges(4,3) - t473 * (t463 * mrSges(4,1) + t499);
t484 = t465 * (m(4) * t428 + t465 * t481 + t498) + t463 * (m(4) * t429 + t463 * t481 + t488);
t476 = m(5) * t419 + t430 * mrSges(5,1) + t431 * mrSges(5,2) + t447 * t441 + t448 * t442 + t357;
t440 = t473 * t500 + t478;
t475 = m(4) * t440 + mrSges(4,1) * t492 + qJDD(1) * t499 + t476;
t391 = Ifges(7,4) * t404 + Ifges(7,2) * t403 + Ifges(7,6) * t445;
t392 = Ifges(7,1) * t404 + Ifges(7,4) * t403 + Ifges(7,5) * t445;
t474 = mrSges(7,1) * t374 - mrSges(7,2) * t375 + Ifges(7,5) * t388 + Ifges(7,6) * t387 + Ifges(7,3) * t427 + t404 * t391 - t403 * t392;
t446 = -qJDD(1) * pkin(1) + t479;
t444 = t473 * pkin(1) - t502;
t422 = Ifges(5,1) * t448 - Ifges(5,4) * t447 + Ifges(5,5) * qJD(4);
t421 = Ifges(5,4) * t448 - Ifges(5,2) * t447 + Ifges(5,6) * qJD(4);
t420 = Ifges(5,5) * t448 - Ifges(5,6) * t447 + Ifges(5,3) * qJD(4);
t401 = Ifges(6,1) * t436 + Ifges(6,4) * t435 + Ifges(6,5) * t447;
t400 = Ifges(6,4) * t436 + Ifges(6,2) * t435 + Ifges(6,6) * t447;
t399 = Ifges(6,5) * t436 + Ifges(6,6) * t435 + Ifges(6,3) * t447;
t390 = Ifges(7,5) * t404 + Ifges(7,6) * t403 + Ifges(7,3) * t445;
t365 = mrSges(7,2) * t380 - mrSges(7,3) * t374 + Ifges(7,1) * t388 + Ifges(7,4) * t387 + Ifges(7,5) * t427 + t403 * t390 - t445 * t391;
t364 = -mrSges(7,1) * t380 + mrSges(7,3) * t375 + Ifges(7,4) * t388 + Ifges(7,2) * t387 + Ifges(7,6) * t427 - t404 * t390 + t445 * t392;
t353 = mrSges(6,2) * t382 - mrSges(6,3) * t378 + Ifges(6,1) * t417 + Ifges(6,4) * t416 + Ifges(6,5) * t430 - pkin(8) * t363 - t466 * t364 + t469 * t365 + t435 * t399 - t447 * t400;
t352 = -mrSges(6,1) * t382 + mrSges(6,3) * t379 + Ifges(6,4) * t417 + Ifges(6,2) * t416 + Ifges(6,6) * t430 - pkin(5) * t477 + pkin(8) * t486 + t469 * t364 + t466 * t365 - t436 * t399 + t447 * t401;
t349 = -t474 - pkin(4) * t357 + (-Ifges(5,2) - Ifges(6,3)) * t430 + Ifges(5,6) * qJDD(4) - t448 * t420 + t435 * t401 - t436 * t400 + Ifges(5,4) * t431 + qJD(4) * t422 - Ifges(6,6) * t416 - Ifges(6,5) * t417 - mrSges(5,1) * t419 - pkin(5) * t363 - mrSges(6,1) * t378 + mrSges(6,2) * t379 + mrSges(5,3) * t396;
t348 = m(3) * t446 + qJDD(1) * mrSges(3,2) - t473 * mrSges(3,3) + t484;
t347 = mrSges(5,2) * t419 - mrSges(5,3) * t395 + Ifges(5,1) * t431 - Ifges(5,4) * t430 + Ifges(5,5) * qJDD(4) - qJ(5) * t357 - qJD(4) * t421 - t462 * t352 + t464 * t353 - t447 * t420;
t1 = [mrSges(2,1) * t490 - mrSges(2,2) * t485 + mrSges(3,2) * t446 - mrSges(3,3) * t444 + t465 * (mrSges(4,2) * t440 - mrSges(4,3) * t428 - pkin(7) * t498 + t470 * t347 - t467 * t349) - t463 * (-mrSges(4,1) * t440 + mrSges(4,3) * t429 - pkin(3) * t476 + pkin(7) * t488 + t467 * t347 + t470 * t349) - qJ(3) * t484 - pkin(1) * t348 + qJ(2) * (-m(3) * t444 + (mrSges(3,2) + t489) * t473 + t475) + (Ifges(4,1) * t503 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t465 + Ifges(4,2) * t463) * t463) * qJDD(1); t348; t473 * t489 + t475; mrSges(5,1) * t395 - mrSges(5,2) * t396 + Ifges(5,5) * t431 - Ifges(5,6) * t430 + Ifges(5,3) * qJDD(4) - pkin(4) * t373 + qJ(5) * t487 + t464 * t352 + t462 * t353 + t448 * t421 + t447 * t422; t373; t474;];
tauJ  = t1;
