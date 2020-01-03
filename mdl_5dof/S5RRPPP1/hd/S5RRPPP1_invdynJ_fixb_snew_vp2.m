% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:40
% EndTime: 2019-12-31 19:23:43
% DurationCPUTime: 2.06s
% Computational Cost: add. (11769->254), mult. (29613->311), div. (0->0), fcn. (20624->8), ass. (0->109)
t490 = -2 * qJD(3);
t489 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t488 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t465 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t487 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t464 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t486 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t439 = cos(qJ(2));
t478 = cos(pkin(5));
t452 = qJD(1) * t478;
t436 = sin(pkin(5));
t469 = qJD(2) * t436;
t414 = (t439 * t452 + t469) * qJ(3);
t441 = qJD(1) ^ 2;
t438 = sin(qJ(1));
t440 = cos(qJ(1));
t457 = t438 * g(1) - t440 * g(2);
t418 = -qJDD(1) * pkin(1) - t441 * pkin(7) - t457;
t437 = sin(qJ(2));
t475 = qJ(3) * t437;
t424 = qJD(2) * pkin(2) - t452 * t475;
t467 = qJD(1) * qJD(2);
t426 = t437 * qJDD(1) + t439 * t467;
t427 = t439 * qJDD(1) - t437 * t467;
t476 = qJ(3) * t436;
t366 = -t426 * t476 - t427 * pkin(2) + (-t414 * t439 + t424 * t437) * qJD(1) + t418;
t450 = -t440 * g(1) - t438 * g(2);
t419 = -t441 * pkin(1) + qJDD(1) * pkin(7) + t450;
t420 = (-t439 * pkin(2) - t436 * t475) * qJD(1);
t456 = qJ(3) * t478;
t481 = t439 * g(3);
t367 = -t426 * t456 + qJDD(2) * pkin(2) - t481 + qJD(2) * t414 + (-qJD(1) * t420 - t419) * t437;
t407 = -t437 * g(3) + t439 * t419;
t446 = qJDD(2) * t436 + t478 * t427;
t470 = qJD(1) * t439;
t368 = t446 * qJ(3) - qJD(2) * t424 + t420 * t470 + t407;
t435 = sin(pkin(8));
t455 = t435 * t478;
t477 = cos(pkin(8));
t405 = t435 * t469 + (t477 * t437 + t439 * t455) * qJD(1);
t449 = t478 * t477;
t454 = t436 * t477;
t356 = t366 * t454 + t367 * t449 - t435 * t368 + t405 * t490;
t471 = qJD(1) * t437;
t404 = -qJD(2) * t454 + t435 * t471 - t449 * t470;
t384 = t404 * pkin(3) - t405 * qJ(4);
t411 = t478 * qJDD(2) - t436 * t427;
t422 = -t478 * qJD(2) + t436 * t470;
t421 = t422 ^ 2;
t353 = -t411 * pkin(3) - t421 * qJ(4) + t405 * t384 + qJDD(4) - t356;
t386 = -t404 * mrSges(5,2) - t405 * mrSges(5,3);
t398 = t477 * t426 + t446 * t435;
t485 = -m(5) * t353 - t398 * mrSges(5,1) - t405 * t386;
t383 = -t405 * mrSges(6,2) + t404 * mrSges(6,3);
t385 = t404 * mrSges(4,1) + t405 * mrSges(4,2);
t474 = t404 * t422;
t483 = 2 * qJD(5);
t347 = t422 * t483 + (t404 * t405 - t411) * qJ(5) + (t398 - t474) * pkin(4) + t353;
t392 = -t404 * mrSges(6,1) - t422 * mrSges(6,2);
t451 = -m(6) * t347 + t411 * mrSges(6,3) - t422 * t392;
t391 = t404 * mrSges(5,1) + t422 * mrSges(5,3);
t472 = t422 * mrSges(4,2) - t404 * mrSges(4,3) - t391;
t479 = -mrSges(6,1) - mrSges(4,3);
t480 = mrSges(4,1) - mrSges(5,2);
t338 = m(4) * t356 - t472 * t422 + t480 * t411 + (-t383 - t385) * t405 + t479 * t398 + t451 + t485;
t401 = t404 * t490;
t461 = t436 * t435 * t366 + t367 * t455 + t477 * t368;
t357 = t401 + t461;
t388 = -t422 * mrSges(4,1) - t405 * mrSges(4,3);
t397 = -qJDD(2) * t454 + t435 * t426 - t427 * t449;
t444 = t421 * pkin(3) - t411 * qJ(4) - t461;
t352 = 0.2e1 * qJD(4) * t422 + ((2 * qJD(3)) + t384) * t404 + t444;
t393 = t405 * mrSges(5,1) - t422 * mrSges(5,2);
t389 = t405 * pkin(4) + t422 * qJ(5);
t403 = t404 ^ 2;
t484 = -0.2e1 * qJD(4);
t349 = -t397 * pkin(4) - t403 * qJ(5) - t404 * t384 + qJDD(5) + t401 + (t484 - t389) * t422 - t444;
t390 = t405 * mrSges(6,1) + t422 * mrSges(6,3);
t463 = m(6) * t349 + t411 * mrSges(6,2) - t422 * t390;
t445 = -m(5) * t352 + t411 * mrSges(5,3) - t422 * t393 + t463;
t473 = -t383 - t386;
t341 = m(4) * t357 - t411 * mrSges(4,2) + t422 * t388 + (-t385 + t473) * t404 + (-mrSges(5,1) + t479) * t397 + t445;
t358 = t478 * t366 - t436 * t367 + qJDD(3);
t442 = (-t398 - t474) * qJ(4) + t358 + (-pkin(3) * t422 + t484) * t405;
t355 = t397 * pkin(3) + t442;
t351 = -t403 * pkin(4) + t404 * t483 - t405 * t389 + (pkin(3) + qJ(5)) * t397 + t442;
t462 = m(6) * t351 + t397 * mrSges(6,3) + t404 * t392;
t448 = m(5) * t355 - t398 * mrSges(5,3) - t405 * t393 + t462;
t342 = m(4) * t358 + (t388 - t390) * t405 + t472 * t404 + (mrSges(4,2) - mrSges(6,2)) * t398 + t480 * t397 + t448;
t333 = (t477 * t338 + t341 * t435) * t436 + t478 * t342;
t460 = t464 * t404 - t465 * t405 + t486 * t422;
t459 = t487 * t404 + t488 * t405 - t464 * t422;
t458 = t488 * t404 - t489 * t405 + t465 * t422;
t336 = -t435 * t338 + t477 * t341;
t334 = t338 * t449 + t341 * t455 - t436 * t342;
t345 = t398 * mrSges(6,1) + t405 * t383 - t451;
t430 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t470;
t429 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t471;
t425 = (-t439 * mrSges(3,1) + t437 * mrSges(3,2)) * qJD(1);
t417 = Ifges(3,5) * qJD(2) + (t437 * Ifges(3,1) + t439 * Ifges(3,4)) * qJD(1);
t416 = Ifges(3,6) * qJD(2) + (t437 * Ifges(3,4) + t439 * Ifges(3,2)) * qJD(1);
t406 = -t437 * t419 - t481;
t346 = -t397 * mrSges(6,1) - t404 * t383 + t463;
t344 = -t397 * mrSges(5,2) - t398 * mrSges(6,2) - t405 * t390 - t404 * t391 + t448;
t343 = t411 * mrSges(5,2) - t422 * t391 + t345 - t485;
t335 = mrSges(5,1) * t353 + mrSges(6,1) * t347 + mrSges(4,2) * t358 - mrSges(6,2) * t351 - mrSges(4,3) * t356 - mrSges(5,3) * t355 + pkin(4) * t345 - qJ(4) * t344 - t397 * t488 + t489 * t398 + t460 * t404 + t465 * t411 + t459 * t422;
t332 = -mrSges(4,1) * t358 + mrSges(4,3) * t357 - mrSges(5,1) * t352 + mrSges(5,2) * t355 + mrSges(6,1) * t349 - mrSges(6,3) * t351 + pkin(4) * t346 - qJ(5) * t462 - pkin(3) * t344 + t458 * t422 + t464 * t411 + (qJ(5) * t390 + t460) * t405 + (qJ(5) * mrSges(6,2) + t488) * t398 + t487 * t397;
t331 = mrSges(4,1) * t356 - mrSges(4,2) * t357 + mrSges(5,2) * t353 - mrSges(5,3) * t352 + mrSges(6,2) * t349 - mrSges(6,3) * t347 - qJ(5) * t345 - pkin(3) * t343 + qJ(4) * t445 + t486 * t411 + t459 * t405 + t465 * t398 + (qJ(4) * t473 - t458) * t404 + (qJ(4) * (-mrSges(5,1) - mrSges(6,1)) - t464) * t397;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t457 - mrSges(2,2) * t450 + t437 * (mrSges(3,2) * t418 - mrSges(3,3) * t406 + Ifges(3,1) * t426 + Ifges(3,4) * t427 + Ifges(3,5) * qJDD(2) - qJD(2) * t416 - t435 * t332 - t333 * t476 - t334 * t456 + t477 * t335) + t439 * (-mrSges(3,1) * t418 + mrSges(3,3) * t407 + Ifges(3,4) * t426 + Ifges(3,2) * t427 + Ifges(3,6) * qJDD(2) - pkin(2) * t333 + qJD(2) * t417 - t436 * t331 + t332 * t449 + t335 * t455 + t336 * t456) + pkin(1) * (-m(3) * t418 + t427 * mrSges(3,1) - t426 * mrSges(3,2) + (-t429 * t437 + t430 * t439) * qJD(1) - t333) + pkin(7) * (t439 * (m(3) * t407 - qJDD(2) * mrSges(3,2) + t427 * mrSges(3,3) - qJD(2) * t429 + t425 * t470 + t336) - t437 * (m(3) * t406 + qJDD(2) * mrSges(3,1) - t426 * mrSges(3,3) + qJD(2) * t430 - t425 * t471 + t334)); mrSges(3,1) * t406 - mrSges(3,2) * t407 + t478 * t331 + Ifges(3,5) * t426 + Ifges(3,6) * t427 + Ifges(3,3) * qJDD(2) + pkin(2) * t334 + (t437 * t416 - t439 * t417) * qJD(1) + (qJ(3) * t336 + t477 * t332 + t335 * t435) * t436; t342; t343; t346;];
tauJ = t1;
