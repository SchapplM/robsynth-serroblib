% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:02:22
% EndTime: 2019-05-05 11:02:27
% DurationCPUTime: 4.59s
% Computational Cost: add. (54264->292), mult. (107076->370), div. (0->0), fcn. (78953->14), ass. (0->124)
t475 = sin(pkin(12));
t477 = cos(pkin(12));
t465 = g(1) * t475 - g(2) * t477;
t474 = -g(3) + qJDD(1);
t476 = sin(pkin(6));
t478 = cos(pkin(6));
t509 = t465 * t478 + t474 * t476;
t466 = -g(1) * t477 - g(2) * t475;
t483 = sin(qJ(2));
t488 = cos(qJ(2));
t425 = -t483 * t466 + t488 * t509;
t426 = t488 * t466 + t483 * t509;
t490 = qJD(2) ^ 2;
t421 = -pkin(2) * t490 + qJDD(2) * pkin(8) + t426;
t443 = -t465 * t476 + t474 * t478;
t482 = sin(qJ(3));
t487 = cos(qJ(3));
t414 = t487 * t421 + t482 * t443;
t462 = (-pkin(3) * t487 - pkin(9) * t482) * qJD(2);
t489 = qJD(3) ^ 2;
t505 = qJD(2) * t487;
t403 = -pkin(3) * t489 + qJDD(3) * pkin(9) + t462 * t505 + t414;
t420 = -qJDD(2) * pkin(2) - t490 * pkin(8) - t425;
t504 = qJD(2) * qJD(3);
t503 = t487 * t504;
t463 = qJDD(2) * t482 + t503;
t473 = t482 * t504;
t464 = qJDD(2) * t487 - t473;
t408 = (-t463 - t503) * pkin(9) + (-t464 + t473) * pkin(3) + t420;
t481 = sin(qJ(4));
t486 = cos(qJ(4));
t386 = -t481 * t403 + t486 * t408;
t506 = qJD(2) * t482;
t459 = qJD(3) * t486 - t481 * t506;
t434 = qJD(4) * t459 + qJDD(3) * t481 + t463 * t486;
t456 = qJDD(4) - t464;
t460 = qJD(3) * t481 + t486 * t506;
t472 = qJD(4) - t505;
t382 = (t459 * t472 - t434) * pkin(10) + (t459 * t460 + t456) * pkin(4) + t386;
t387 = t486 * t403 + t481 * t408;
t433 = -qJD(4) * t460 + qJDD(3) * t486 - t463 * t481;
t442 = pkin(4) * t472 - pkin(10) * t460;
t455 = t459 ^ 2;
t384 = -pkin(4) * t455 + pkin(10) * t433 - t442 * t472 + t387;
t480 = sin(qJ(5));
t485 = cos(qJ(5));
t370 = t485 * t382 - t480 * t384;
t436 = t459 * t485 - t460 * t480;
t400 = qJD(5) * t436 + t433 * t480 + t434 * t485;
t437 = t459 * t480 + t460 * t485;
t452 = qJDD(5) + t456;
t471 = qJD(5) + t472;
t367 = (t436 * t471 - t400) * pkin(11) + (t436 * t437 + t452) * pkin(5) + t370;
t371 = t480 * t382 + t485 * t384;
t399 = -qJD(5) * t437 + t433 * t485 - t434 * t480;
t424 = pkin(5) * t471 - pkin(11) * t437;
t435 = t436 ^ 2;
t368 = -pkin(5) * t435 + pkin(11) * t399 - t424 * t471 + t371;
t479 = sin(qJ(6));
t484 = cos(qJ(6));
t365 = t367 * t484 - t368 * t479;
t415 = t436 * t484 - t437 * t479;
t379 = qJD(6) * t415 + t399 * t479 + t400 * t484;
t416 = t436 * t479 + t437 * t484;
t394 = -mrSges(7,1) * t415 + mrSges(7,2) * t416;
t467 = qJD(6) + t471;
t406 = -mrSges(7,2) * t467 + mrSges(7,3) * t415;
t447 = qJDD(6) + t452;
t362 = m(7) * t365 + mrSges(7,1) * t447 - mrSges(7,3) * t379 - t394 * t416 + t406 * t467;
t366 = t367 * t479 + t368 * t484;
t378 = -qJD(6) * t416 + t399 * t484 - t400 * t479;
t407 = mrSges(7,1) * t467 - mrSges(7,3) * t416;
t363 = m(7) * t366 - mrSges(7,2) * t447 + mrSges(7,3) * t378 + t394 * t415 - t407 * t467;
t355 = t484 * t362 + t479 * t363;
t417 = -mrSges(6,1) * t436 + mrSges(6,2) * t437;
t422 = -mrSges(6,2) * t471 + mrSges(6,3) * t436;
t352 = m(6) * t370 + mrSges(6,1) * t452 - mrSges(6,3) * t400 - t417 * t437 + t422 * t471 + t355;
t423 = mrSges(6,1) * t471 - mrSges(6,3) * t437;
t499 = -t362 * t479 + t484 * t363;
t353 = m(6) * t371 - mrSges(6,2) * t452 + mrSges(6,3) * t399 + t417 * t436 - t423 * t471 + t499;
t348 = t485 * t352 + t480 * t353;
t461 = (-mrSges(4,1) * t487 + mrSges(4,2) * t482) * qJD(2);
t468 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t506;
t438 = -mrSges(5,1) * t459 + mrSges(5,2) * t460;
t440 = -mrSges(5,2) * t472 + mrSges(5,3) * t459;
t346 = m(5) * t386 + mrSges(5,1) * t456 - mrSges(5,3) * t434 - t438 * t460 + t440 * t472 + t348;
t441 = mrSges(5,1) * t472 - mrSges(5,3) * t460;
t500 = -t352 * t480 + t485 * t353;
t347 = m(5) * t387 - mrSges(5,2) * t456 + mrSges(5,3) * t433 + t438 * t459 - t441 * t472 + t500;
t501 = -t346 * t481 + t486 * t347;
t341 = m(4) * t414 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t464 - qJD(3) * t468 + t461 * t505 + t501;
t413 = -t482 * t421 + t443 * t487;
t469 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t505;
t402 = -qJDD(3) * pkin(3) - pkin(9) * t489 + t462 * t506 - t413;
t388 = -pkin(4) * t433 - pkin(10) * t455 + t460 * t442 + t402;
t373 = -pkin(5) * t399 - pkin(11) * t435 + t424 * t437 + t388;
t498 = m(7) * t373 - t378 * mrSges(7,1) + t379 * mrSges(7,2) - t415 * t406 + t416 * t407;
t495 = m(6) * t388 - t399 * mrSges(6,1) + t400 * mrSges(6,2) - t436 * t422 + t437 * t423 + t498;
t492 = -m(5) * t402 + t433 * mrSges(5,1) - t434 * mrSges(5,2) + t459 * t440 - t460 * t441 - t495;
t358 = m(4) * t413 + qJDD(3) * mrSges(4,1) - t463 * mrSges(4,3) + qJD(3) * t469 - t461 * t506 + t492;
t502 = t487 * t341 - t358 * t482;
t342 = t346 * t486 + t347 * t481;
t390 = Ifges(7,4) * t416 + Ifges(7,2) * t415 + Ifges(7,6) * t467;
t391 = Ifges(7,1) * t416 + Ifges(7,4) * t415 + Ifges(7,5) * t467;
t496 = -mrSges(7,1) * t365 + mrSges(7,2) * t366 - Ifges(7,5) * t379 - Ifges(7,6) * t378 - Ifges(7,3) * t447 - t416 * t390 + t415 * t391;
t494 = -m(4) * t420 + t464 * mrSges(4,1) - mrSges(4,2) * t463 - t468 * t506 + t469 * t505 - t342;
t410 = Ifges(6,4) * t437 + Ifges(6,2) * t436 + Ifges(6,6) * t471;
t411 = Ifges(6,1) * t437 + Ifges(6,4) * t436 + Ifges(6,5) * t471;
t493 = -mrSges(6,1) * t370 + mrSges(6,2) * t371 - Ifges(6,5) * t400 - Ifges(6,6) * t399 - Ifges(6,3) * t452 - pkin(5) * t355 - t437 * t410 + t436 * t411 + t496;
t428 = Ifges(5,4) * t460 + Ifges(5,2) * t459 + Ifges(5,6) * t472;
t429 = Ifges(5,1) * t460 + Ifges(5,4) * t459 + Ifges(5,5) * t472;
t491 = mrSges(5,1) * t386 - mrSges(5,2) * t387 + Ifges(5,5) * t434 + Ifges(5,6) * t433 + Ifges(5,3) * t456 + pkin(4) * t348 + t460 * t428 - t459 * t429 - t493;
t451 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t482 + Ifges(4,4) * t487) * qJD(2);
t450 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t482 + Ifges(4,2) * t487) * qJD(2);
t427 = Ifges(5,5) * t460 + Ifges(5,6) * t459 + Ifges(5,3) * t472;
t409 = Ifges(6,5) * t437 + Ifges(6,6) * t436 + Ifges(6,3) * t471;
t389 = Ifges(7,5) * t416 + Ifges(7,6) * t415 + Ifges(7,3) * t467;
t357 = mrSges(7,2) * t373 - mrSges(7,3) * t365 + Ifges(7,1) * t379 + Ifges(7,4) * t378 + Ifges(7,5) * t447 + t389 * t415 - t390 * t467;
t356 = -mrSges(7,1) * t373 + mrSges(7,3) * t366 + Ifges(7,4) * t379 + Ifges(7,2) * t378 + Ifges(7,6) * t447 - t389 * t416 + t391 * t467;
t344 = mrSges(6,2) * t388 - mrSges(6,3) * t370 + Ifges(6,1) * t400 + Ifges(6,4) * t399 + Ifges(6,5) * t452 - pkin(11) * t355 - t356 * t479 + t357 * t484 + t409 * t436 - t410 * t471;
t343 = -mrSges(6,1) * t388 + mrSges(6,3) * t371 + Ifges(6,4) * t400 + Ifges(6,2) * t399 + Ifges(6,6) * t452 - pkin(5) * t498 + pkin(11) * t499 + t484 * t356 + t479 * t357 - t437 * t409 + t471 * t411;
t339 = mrSges(5,2) * t402 - mrSges(5,3) * t386 + Ifges(5,1) * t434 + Ifges(5,4) * t433 + Ifges(5,5) * t456 - pkin(10) * t348 - t343 * t480 + t344 * t485 + t427 * t459 - t428 * t472;
t338 = -mrSges(5,1) * t402 + mrSges(5,3) * t387 + Ifges(5,4) * t434 + Ifges(5,2) * t433 + Ifges(5,6) * t456 - pkin(4) * t495 + pkin(10) * t500 + t485 * t343 + t480 * t344 - t460 * t427 + t472 * t429;
t1 = [m(2) * t474 + t478 * (m(3) * t443 + t341 * t482 + t358 * t487) + (t483 * (m(3) * t426 - mrSges(3,1) * t490 - qJDD(2) * mrSges(3,2) + t502) + t488 * (m(3) * t425 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t490 + t494)) * t476; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t425 - mrSges(3,2) * t426 + t482 * (mrSges(4,2) * t420 - mrSges(4,3) * t413 + Ifges(4,1) * t463 + Ifges(4,4) * t464 + Ifges(4,5) * qJDD(3) - pkin(9) * t342 - qJD(3) * t450 - t338 * t481 + t339 * t486) + t487 * (-mrSges(4,1) * t420 + mrSges(4,3) * t414 + Ifges(4,4) * t463 + Ifges(4,2) * t464 + Ifges(4,6) * qJDD(3) - pkin(3) * t342 + qJD(3) * t451 - t491) + pkin(2) * t494 + pkin(8) * t502; Ifges(4,5) * t463 + Ifges(4,6) * t464 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t413 - mrSges(4,2) * t414 + t481 * t339 + t486 * t338 + pkin(3) * t492 + pkin(9) * t501 + (t450 * t482 - t451 * t487) * qJD(2); t491; -t493; -t496;];
tauJ  = t1;
