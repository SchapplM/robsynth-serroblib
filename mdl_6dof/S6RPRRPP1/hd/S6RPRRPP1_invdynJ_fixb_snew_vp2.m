% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:15:42
% EndTime: 2019-05-05 21:15:47
% DurationCPUTime: 2.43s
% Computational Cost: add. (19171->268), mult. (37974->328), div. (0->0), fcn. (24599->10), ass. (0->111)
t485 = Ifges(6,1) + Ifges(7,1);
t476 = Ifges(6,4) - Ifges(7,5);
t475 = Ifges(6,5) + Ifges(7,4);
t484 = -Ifges(6,2) - Ifges(7,3);
t474 = Ifges(6,6) - Ifges(7,6);
t483 = -Ifges(7,2) - Ifges(6,3);
t446 = sin(qJ(1));
t449 = cos(qJ(1));
t461 = t446 * g(1) - g(2) * t449;
t428 = qJDD(1) * pkin(1) + t461;
t451 = qJD(1) ^ 2;
t456 = -g(1) * t449 - g(2) * t446;
t430 = -pkin(1) * t451 + t456;
t442 = sin(pkin(9));
t443 = cos(pkin(9));
t406 = t428 * t443 - t442 * t430;
t393 = -qJDD(1) * pkin(2) - pkin(7) * t451 - t406;
t445 = sin(qJ(3));
t448 = cos(qJ(3));
t466 = qJD(1) * qJD(3);
t462 = t448 * t466;
t432 = qJDD(1) * t445 + t462;
t463 = t445 * t466;
t433 = qJDD(1) * t448 - t463;
t370 = (-t432 - t462) * pkin(8) + (-t433 + t463) * pkin(3) + t393;
t407 = t442 * t428 + t443 * t430;
t394 = -pkin(2) * t451 + qJDD(1) * pkin(7) + t407;
t440 = -g(3) + qJDD(2);
t385 = t448 * t394 + t445 * t440;
t431 = (-t448 * pkin(3) - t445 * pkin(8)) * qJD(1);
t450 = qJD(3) ^ 2;
t467 = qJD(1) * t448;
t378 = -pkin(3) * t450 + qJDD(3) * pkin(8) + t431 * t467 + t385;
t444 = sin(qJ(4));
t447 = cos(qJ(4));
t357 = t447 * t370 - t378 * t444;
t468 = qJD(1) * t445;
t426 = qJD(3) * t447 - t444 * t468;
t403 = qJD(4) * t426 + qJDD(3) * t444 + t432 * t447;
t425 = qJDD(4) - t433;
t427 = qJD(3) * t444 + t447 * t468;
t437 = qJD(4) - t467;
t353 = (t426 * t437 - t403) * qJ(5) + (t426 * t427 + t425) * pkin(4) + t357;
t358 = t444 * t370 + t447 * t378;
t402 = -qJD(4) * t427 + qJDD(3) * t447 - t432 * t444;
t410 = pkin(4) * t437 - qJ(5) * t427;
t424 = t426 ^ 2;
t355 = -pkin(4) * t424 + qJ(5) * t402 - t410 * t437 + t358;
t441 = sin(pkin(10));
t473 = cos(pkin(10));
t404 = -t426 * t473 + t427 * t441;
t478 = -2 * qJD(5);
t349 = t441 * t353 + t473 * t355 + t404 * t478;
t371 = -t402 * t473 + t403 * t441;
t405 = t441 * t426 + t427 * t473;
t387 = mrSges(6,1) * t437 - mrSges(6,3) * t405;
t379 = pkin(5) * t404 - qJ(6) * t405;
t436 = t437 ^ 2;
t346 = -pkin(5) * t436 + qJ(6) * t425 + 0.2e1 * qJD(6) * t437 - t379 * t404 + t349;
t388 = -mrSges(7,1) * t437 + mrSges(7,2) * t405;
t464 = m(7) * t346 + t425 * mrSges(7,3) + t437 * t388;
t380 = mrSges(7,1) * t404 - mrSges(7,3) * t405;
t469 = -mrSges(6,1) * t404 - mrSges(6,2) * t405 - t380;
t477 = -mrSges(6,3) - mrSges(7,2);
t337 = m(6) * t349 - mrSges(6,2) * t425 + t371 * t477 - t387 * t437 + t404 * t469 + t464;
t455 = t353 * t473 - t441 * t355;
t348 = t405 * t478 + t455;
t372 = t441 * t402 + t403 * t473;
t386 = -mrSges(6,2) * t437 - mrSges(6,3) * t404;
t347 = -t425 * pkin(5) - t436 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t379) * t405 - t455;
t389 = -mrSges(7,2) * t404 + mrSges(7,3) * t437;
t457 = -m(7) * t347 + t425 * mrSges(7,1) + t437 * t389;
t339 = m(6) * t348 + mrSges(6,1) * t425 + t372 * t477 + t386 * t437 + t405 * t469 + t457;
t332 = t441 * t337 + t473 * t339;
t343 = mrSges(7,2) * t372 + t380 * t405 - t457;
t396 = Ifges(5,4) * t427 + Ifges(5,2) * t426 + Ifges(5,6) * t437;
t397 = Ifges(5,1) * t427 + Ifges(5,4) * t426 + Ifges(5,5) * t437;
t470 = -t476 * t404 + t485 * t405 + t475 * t437;
t471 = t484 * t404 + t476 * t405 + t474 * t437;
t482 = -t371 * t474 + t372 * t475 + (Ifges(5,3) - t483) * t425 + mrSges(5,1) * t357 + mrSges(6,1) * t348 - mrSges(7,1) * t347 - mrSges(5,2) * t358 - mrSges(6,2) * t349 + mrSges(7,3) * t346 + Ifges(5,5) * t403 + Ifges(5,6) * t402 + pkin(4) * t332 - pkin(5) * t343 + qJ(6) * (-mrSges(7,2) * t371 - t380 * t404 + t464) + t427 * t396 - t426 * t397 + t405 * t471 + t404 * t470;
t472 = t474 * t404 - t475 * t405 + t483 * t437;
t429 = (-mrSges(4,1) * t448 + mrSges(4,2) * t445) * qJD(1);
t434 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t468;
t408 = -mrSges(5,1) * t426 + mrSges(5,2) * t427;
t409 = -mrSges(5,2) * t437 + mrSges(5,3) * t426;
t330 = m(5) * t357 + mrSges(5,1) * t425 - mrSges(5,3) * t403 - t408 * t427 + t409 * t437 + t332;
t411 = mrSges(5,1) * t437 - mrSges(5,3) * t427;
t458 = t473 * t337 - t339 * t441;
t331 = m(5) * t358 - mrSges(5,2) * t425 + mrSges(5,3) * t402 + t408 * t426 - t411 * t437 + t458;
t459 = -t330 * t444 + t447 * t331;
t327 = m(4) * t385 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t433 - qJD(3) * t434 + t429 * t467 + t459;
t384 = -t445 * t394 + t440 * t448;
t435 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t467;
t377 = -qJDD(3) * pkin(3) - pkin(8) * t450 + t431 * t468 - t384;
t356 = -pkin(4) * t402 - qJ(5) * t424 + t427 * t410 + qJDD(5) + t377;
t351 = -0.2e1 * qJD(6) * t405 + (t404 * t437 - t372) * qJ(6) + (t405 * t437 + t371) * pkin(5) + t356;
t344 = m(7) * t351 + t371 * mrSges(7,1) - t372 * mrSges(7,3) - t405 * t388 + t404 * t389;
t341 = m(6) * t356 + t371 * mrSges(6,1) + t372 * mrSges(6,2) + t404 * t386 + t405 * t387 + t344;
t453 = -m(5) * t377 + t402 * mrSges(5,1) - t403 * mrSges(5,2) + t426 * t409 - t427 * t411 - t341;
t340 = m(4) * t384 + qJDD(3) * mrSges(4,1) - t432 * mrSges(4,3) + qJD(3) * t435 - t429 * t468 + t453;
t460 = t448 * t327 - t340 * t445;
t328 = t447 * t330 + t444 * t331;
t454 = -m(4) * t393 + t433 * mrSges(4,1) - t432 * mrSges(4,2) - t434 * t468 + t435 * t467 - t328;
t419 = Ifges(4,5) * qJD(3) + (t445 * Ifges(4,1) + t448 * Ifges(4,4)) * qJD(1);
t418 = Ifges(4,6) * qJD(3) + (t445 * Ifges(4,4) + t448 * Ifges(4,2)) * qJD(1);
t395 = Ifges(5,5) * t427 + Ifges(5,6) * t426 + Ifges(5,3) * t437;
t334 = mrSges(6,2) * t356 + mrSges(7,2) * t347 - mrSges(6,3) * t348 - mrSges(7,3) * t351 - qJ(6) * t344 - t476 * t371 + t485 * t372 + t472 * t404 + t475 * t425 - t471 * t437;
t333 = -mrSges(6,1) * t356 - mrSges(7,1) * t351 + mrSges(7,2) * t346 + mrSges(6,3) * t349 - pkin(5) * t344 + t484 * t371 + t476 * t372 + t472 * t405 + t474 * t425 + t470 * t437;
t325 = mrSges(5,2) * t377 - mrSges(5,3) * t357 + Ifges(5,1) * t403 + Ifges(5,4) * t402 + Ifges(5,5) * t425 - qJ(5) * t332 - t441 * t333 + t334 * t473 + t426 * t395 - t437 * t396;
t324 = -mrSges(5,1) * t377 + mrSges(5,3) * t358 + Ifges(5,4) * t403 + Ifges(5,2) * t402 + Ifges(5,6) * t425 - pkin(4) * t341 + qJ(5) * t458 + t333 * t473 + t441 * t334 - t427 * t395 + t437 * t397;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t461 - mrSges(2,2) * t456 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t406 - mrSges(3,2) * t407 + t445 * (mrSges(4,2) * t393 - mrSges(4,3) * t384 + Ifges(4,1) * t432 + Ifges(4,4) * t433 + Ifges(4,5) * qJDD(3) - pkin(8) * t328 - qJD(3) * t418 - t444 * t324 + t447 * t325) + pkin(2) * t454 + pkin(7) * t460 + pkin(1) * (t442 * (m(3) * t407 - mrSges(3,1) * t451 - qJDD(1) * mrSges(3,2) + t460) + t443 * (m(3) * t406 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t451 + t454)) + (-mrSges(4,1) * t393 + mrSges(4,3) * t385 + Ifges(4,4) * t432 + Ifges(4,2) * t433 + Ifges(4,6) * qJDD(3) - pkin(3) * t328 + qJD(3) * t419 - t482) * t448; m(3) * t440 + t327 * t445 + t340 * t448; Ifges(4,5) * t432 + Ifges(4,6) * t433 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t384 - mrSges(4,2) * t385 + t444 * t325 + t447 * t324 + pkin(3) * t453 + pkin(8) * t459 + (t445 * t418 - t448 * t419) * qJD(1); t482; t341; t343;];
tauJ  = t1;
