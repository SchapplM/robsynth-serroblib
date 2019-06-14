% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:07:57
% EndTime: 2019-05-06 01:08:01
% DurationCPUTime: 2.26s
% Computational Cost: add. (19039->269), mult. (37159->332), div. (0->0), fcn. (24583->10), ass. (0->112)
t502 = Ifges(6,1) + Ifges(7,1);
t496 = Ifges(6,4) - Ifges(7,5);
t495 = -Ifges(6,5) - Ifges(7,4);
t501 = Ifges(6,2) + Ifges(7,3);
t494 = Ifges(6,6) - Ifges(7,6);
t500 = -Ifges(6,3) - Ifges(7,2);
t464 = sin(qJ(4));
t465 = sin(qJ(3));
t467 = cos(qJ(4));
t468 = cos(qJ(3));
t437 = (t464 * t465 - t467 * t468) * qJD(1);
t486 = qJD(1) * qJD(3);
t484 = t468 * t486;
t445 = qJDD(1) * t465 + t484;
t446 = qJDD(1) * t468 - t465 * t486;
t409 = -qJD(4) * t437 + t445 * t467 + t446 * t464;
t438 = (t464 * t468 + t465 * t467) * qJD(1);
t458 = qJD(3) + qJD(4);
t463 = sin(qJ(5));
t498 = cos(qJ(5));
t425 = t463 * t438 - t458 * t498;
t457 = qJDD(3) + qJDD(4);
t378 = -t425 * qJD(5) + t409 * t498 + t463 * t457;
t426 = t438 * t498 + t463 * t458;
t396 = mrSges(7,1) * t425 - mrSges(7,3) * t426;
t466 = sin(qJ(1));
t469 = cos(qJ(1));
t483 = t466 * g(1) - g(2) * t469;
t442 = qJDD(1) * pkin(1) + t483;
t470 = qJD(1) ^ 2;
t478 = -g(1) * t469 - g(2) * t466;
t444 = -pkin(1) * t470 + t478;
t461 = sin(pkin(10));
t462 = cos(pkin(10));
t424 = t461 * t442 + t462 * t444;
t420 = -pkin(2) * t470 + qJDD(1) * pkin(7) + t424;
t460 = -g(3) + qJDD(2);
t406 = -t465 * t420 + t468 * t460;
t382 = (-t445 + t484) * pkin(8) + (t465 * t468 * t470 + qJDD(3)) * pkin(3) + t406;
t407 = t468 * t420 + t465 * t460;
t488 = qJD(1) * t465;
t449 = qJD(3) * pkin(3) - pkin(8) * t488;
t459 = t468 ^ 2;
t383 = -pkin(3) * t459 * t470 + pkin(8) * t446 - qJD(3) * t449 + t407;
t374 = t464 * t382 + t467 * t383;
t422 = pkin(4) * t437 - pkin(9) * t438;
t456 = t458 ^ 2;
t369 = -pkin(4) * t456 + pkin(9) * t457 - t422 * t437 + t374;
t423 = t462 * t442 - t461 * t444;
t476 = -qJDD(1) * pkin(2) - t423;
t392 = -t446 * pkin(3) + t449 * t488 + (-pkin(8) * t459 - pkin(7)) * t470 + t476;
t408 = -qJD(4) * t438 - t445 * t464 + t446 * t467;
t371 = (t437 * t458 - t409) * pkin(9) + (t438 * t458 - t408) * pkin(4) + t392;
t365 = -t463 * t369 + t371 * t498;
t395 = pkin(5) * t425 - qJ(6) * t426;
t405 = qJDD(5) - t408;
t430 = qJD(5) + t437;
t429 = t430 ^ 2;
t363 = -t405 * pkin(5) - t429 * qJ(6) + t426 * t395 + qJDD(6) - t365;
t410 = -mrSges(7,2) * t425 + mrSges(7,3) * t430;
t479 = -m(7) * t363 + t405 * mrSges(7,1) + t430 * t410;
t359 = t378 * mrSges(7,2) + t426 * t396 - t479;
t366 = t369 * t498 + t463 * t371;
t362 = -pkin(5) * t429 + qJ(6) * t405 + 0.2e1 * qJD(6) * t430 - t395 * t425 + t366;
t377 = t426 * qJD(5) + t463 * t409 - t457 * t498;
t413 = -mrSges(7,1) * t430 + mrSges(7,2) * t426;
t485 = m(7) * t362 + t405 * mrSges(7,3) + t430 * t413;
t490 = t496 * t425 - t502 * t426 + t495 * t430;
t491 = t501 * t425 - t496 * t426 - t494 * t430;
t499 = -t377 * t494 - t378 * t495 - t500 * t405 - t490 * t425 - t491 * t426 + mrSges(6,1) * t365 - mrSges(7,1) * t363 - mrSges(6,2) * t366 + mrSges(7,3) * t362 - pkin(5) * t359 + qJ(6) * (-t377 * mrSges(7,2) - t425 * t396 + t485);
t497 = -mrSges(6,3) - mrSges(7,2);
t421 = mrSges(5,1) * t437 + mrSges(5,2) * t438;
t428 = mrSges(5,1) * t458 - mrSges(5,3) * t438;
t412 = mrSges(6,1) * t430 - mrSges(6,3) * t426;
t489 = -mrSges(6,1) * t425 - mrSges(6,2) * t426 - t396;
t354 = m(6) * t366 - t405 * mrSges(6,2) + t377 * t497 - t430 * t412 + t489 * t425 + t485;
t411 = -mrSges(6,2) * t430 - mrSges(6,3) * t425;
t356 = m(6) * t365 + t405 * mrSges(6,1) + t378 * t497 + t430 * t411 + t489 * t426 + t479;
t480 = t354 * t498 - t356 * t463;
t343 = m(5) * t374 - mrSges(5,2) * t457 + mrSges(5,3) * t408 - t421 * t437 - t428 * t458 + t480;
t373 = t467 * t382 - t464 * t383;
t427 = -mrSges(5,2) * t458 - mrSges(5,3) * t437;
t368 = -t457 * pkin(4) - t456 * pkin(9) + t438 * t422 - t373;
t364 = -0.2e1 * qJD(6) * t426 + (t425 * t430 - t378) * qJ(6) + (t426 * t430 + t377) * pkin(5) + t368;
t360 = m(7) * t364 + mrSges(7,1) * t377 - t378 * mrSges(7,3) + t410 * t425 - t426 * t413;
t472 = -m(6) * t368 - t377 * mrSges(6,1) - mrSges(6,2) * t378 - t425 * t411 - t412 * t426 - t360;
t351 = m(5) * t373 + mrSges(5,1) * t457 - mrSges(5,3) * t409 - t421 * t438 + t427 * t458 + t472;
t340 = t464 * t343 + t467 * t351;
t349 = t463 * t354 + t356 * t498;
t492 = t494 * t425 + t495 * t426 + t500 * t430;
t487 = qJD(1) * t468;
t443 = (-mrSges(4,1) * t468 + mrSges(4,2) * t465) * qJD(1);
t448 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t487;
t338 = m(4) * t406 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t445 + qJD(3) * t448 - t443 * t488 + t340;
t447 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t488;
t481 = t467 * t343 - t351 * t464;
t339 = m(4) * t407 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t446 - qJD(3) * t447 + t443 * t487 + t481;
t482 = -t338 * t465 + t468 * t339;
t475 = m(5) * t392 - t408 * mrSges(5,1) + mrSges(5,2) * t409 + t437 * t427 + t428 * t438 + t349;
t345 = -mrSges(6,1) * t368 - mrSges(7,1) * t364 + mrSges(7,2) * t362 + mrSges(6,3) * t366 - pkin(5) * t360 - t501 * t377 + t496 * t378 + t494 * t405 + t492 * t426 - t490 * t430;
t347 = mrSges(6,2) * t368 + mrSges(7,2) * t363 - mrSges(6,3) * t365 - mrSges(7,3) * t364 - qJ(6) * t360 - t496 * t377 + t502 * t378 - t495 * t405 + t492 * t425 + t491 * t430;
t416 = Ifges(5,4) * t438 - Ifges(5,2) * t437 + Ifges(5,6) * t458;
t417 = Ifges(5,1) * t438 - Ifges(5,4) * t437 + Ifges(5,5) * t458;
t474 = mrSges(5,1) * t373 - mrSges(5,2) * t374 + Ifges(5,5) * t409 + Ifges(5,6) * t408 + Ifges(5,3) * t457 + pkin(4) * t472 + pkin(9) * t480 + t345 * t498 + t463 * t347 + t438 * t416 + t437 * t417;
t419 = -t470 * pkin(7) + t476;
t471 = -m(4) * t419 + t446 * mrSges(4,1) - mrSges(4,2) * t445 - t447 * t488 + t448 * t487 - t475;
t436 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t465 + Ifges(4,4) * t468) * qJD(1);
t435 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t465 + Ifges(4,2) * t468) * qJD(1);
t415 = Ifges(5,5) * t438 - Ifges(5,6) * t437 + Ifges(5,3) * t458;
t336 = -mrSges(5,1) * t392 + mrSges(5,3) * t374 + Ifges(5,4) * t409 + Ifges(5,2) * t408 + Ifges(5,6) * t457 - pkin(4) * t349 - t438 * t415 + t458 * t417 - t499;
t335 = mrSges(5,2) * t392 - mrSges(5,3) * t373 + Ifges(5,1) * t409 + Ifges(5,4) * t408 + Ifges(5,5) * t457 - pkin(9) * t349 - t463 * t345 + t347 * t498 - t437 * t415 - t458 * t416;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t483 - mrSges(2,2) * t478 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t423 - mrSges(3,2) * t424 + t465 * (mrSges(4,2) * t419 - mrSges(4,3) * t406 + Ifges(4,1) * t445 + Ifges(4,4) * t446 + Ifges(4,5) * qJDD(3) - pkin(8) * t340 - qJD(3) * t435 + t467 * t335 - t464 * t336) + t468 * (-mrSges(4,1) * t419 + mrSges(4,3) * t407 + Ifges(4,4) * t445 + Ifges(4,2) * t446 + Ifges(4,6) * qJDD(3) - pkin(3) * t475 + pkin(8) * t481 + qJD(3) * t436 + t464 * t335 + t467 * t336) + pkin(2) * t471 + pkin(7) * t482 + pkin(1) * (t461 * (m(3) * t424 - mrSges(3,1) * t470 - qJDD(1) * mrSges(3,2) + t482) + t462 * (m(3) * t423 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t470 + t471)); m(3) * t460 + t338 * t468 + t339 * t465; Ifges(4,3) * qJDD(3) + t474 + (t435 * t465 - t436 * t468) * qJD(1) + Ifges(4,6) * t446 + Ifges(4,5) * t445 + mrSges(4,1) * t406 - mrSges(4,2) * t407 + pkin(3) * t340; t474; t499; t359;];
tauJ  = t1;
