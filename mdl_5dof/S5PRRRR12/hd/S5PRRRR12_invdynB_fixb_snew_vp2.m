% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:51
% EndTime: 2024-09-28 18:07:58
% DurationCPUTime: 3.94s
% Computational Cost: add. (126981->193), mult. (172245->266), div. (0->0), fcn. (131563->14), ass. (0->103)
t443 = sin(pkin(6));
t476 = pkin(10) * t443;
t440 = qJD(2) + qJD(3);
t437 = qJD(4) + t440;
t475 = t437 * t443;
t448 = sin(qJ(5));
t474 = t443 * t448;
t452 = cos(qJ(5));
t473 = t443 * t452;
t444 = sin(pkin(5));
t451 = sin(qJ(2));
t472 = t444 * t451;
t455 = cos(qJ(2));
t471 = t444 * t455;
t447 = cos(pkin(5));
t470 = t447 * t451;
t469 = t447 * t455;
t442 = sin(pkin(11));
t445 = cos(pkin(11));
t431 = t442 * g(1) - t445 * g(2);
t432 = -t445 * g(1) - t442 * g(2);
t441 = -g(3) + qJDD(1);
t413 = t431 * t469 - t451 * t432 + t441 * t471;
t411 = qJDD(2) * pkin(2) + t413;
t414 = t431 * t470 + t455 * t432 + t441 * t472;
t456 = qJD(2) ^ 2;
t412 = -t456 * pkin(2) + t414;
t450 = sin(qJ(3));
t454 = cos(qJ(3));
t406 = t454 * t411 - t450 * t412;
t439 = qJDD(2) + qJDD(3);
t404 = t439 * pkin(3) + t406;
t407 = t450 * t411 + t454 * t412;
t438 = t440 ^ 2;
t405 = -t438 * pkin(3) + t407;
t449 = sin(qJ(4));
t453 = cos(qJ(4));
t400 = t449 * t404 + t453 * t405;
t435 = t437 ^ 2;
t436 = qJDD(4) + t439;
t398 = -t435 * pkin(4) + t436 * t476 + t400;
t399 = t453 * t404 - t449 * t405;
t397 = t436 * pkin(4) + t435 * t476 + t399;
t425 = -t444 * t431 + t447 * t441;
t446 = cos(pkin(6));
t459 = t397 * t446 + t425 * t443;
t392 = -t448 * t398 + t459 * t452;
t430 = t446 * t437 + qJD(5);
t465 = t437 * t473;
t419 = -t430 * mrSges(6,2) + mrSges(6,3) * t465;
t420 = (-mrSges(6,1) * t452 + mrSges(6,2) * t448) * t475;
t467 = qJD(5) * t437;
t421 = (t436 * t448 + t452 * t467) * t443;
t429 = t446 * t436 + qJDD(5);
t466 = t437 * t474;
t390 = m(6) * t392 + t429 * mrSges(6,1) - t421 * mrSges(6,3) + t430 * t419 - t420 * t466;
t393 = t452 * t398 + t459 * t448;
t418 = t430 * mrSges(6,1) - mrSges(6,3) * t466;
t422 = (t436 * t452 - t448 * t467) * t443;
t391 = m(6) * t393 - t429 * mrSges(6,2) + t422 * mrSges(6,3) - t430 * t418 + t420 * t465;
t396 = -t443 * t397 + t446 * t425;
t395 = m(6) * t396 - t422 * mrSges(6,1) + t421 * mrSges(6,2) + (t418 * t448 - t419 * t452) * t475;
t377 = -t443 * t395 + (t390 * t452 + t391 * t448) * t446;
t373 = m(5) * t399 + t436 * mrSges(5,1) - t435 * mrSges(5,2) + t377;
t381 = -t448 * t390 + t452 * t391;
t380 = m(5) * t400 - t435 * mrSges(5,1) - t436 * mrSges(5,2) + t381;
t371 = t453 * t373 + t449 * t380;
t369 = m(4) * t406 + t439 * mrSges(4,1) - t438 * mrSges(4,2) + t371;
t462 = -t449 * t373 + t453 * t380;
t370 = m(4) * t407 - t438 * mrSges(4,1) - t439 * mrSges(4,2) + t462;
t363 = t454 * t369 + t450 * t370;
t361 = m(3) * t413 + qJDD(2) * mrSges(3,1) - t456 * mrSges(3,2) + t363;
t463 = -t450 * t369 + t454 * t370;
t362 = m(3) * t414 - t456 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t463;
t376 = t390 * t473 + t391 * t474 + t446 * t395;
t461 = m(5) * t425 + t376;
t460 = m(4) * t425 + t461;
t375 = m(3) * t425 + t460;
t350 = t361 * t469 + t362 * t470 - t444 * t375;
t348 = m(2) * t431 + t350;
t355 = -t451 * t361 + t455 * t362;
t354 = m(2) * t432 + t355;
t468 = t445 * t348 + t442 * t354;
t349 = t361 * t471 + t362 * t472 + t447 * t375;
t464 = -t442 * t348 + t445 * t354;
t416 = Ifges(6,6) * t430 + (Ifges(6,4) * t448 + Ifges(6,2) * t452) * t475;
t417 = Ifges(6,5) * t430 + (Ifges(6,1) * t448 + Ifges(6,4) * t452) * t475;
t382 = mrSges(6,1) * t392 - mrSges(6,2) * t393 + Ifges(6,5) * t421 + Ifges(6,6) * t422 + Ifges(6,3) * t429 + (t416 * t448 - t417 * t452) * t475;
t415 = Ifges(6,3) * t430 + (Ifges(6,5) * t448 + Ifges(6,6) * t452) * t475;
t383 = -mrSges(6,1) * t396 + mrSges(6,3) * t393 + Ifges(6,4) * t421 + Ifges(6,2) * t422 + Ifges(6,6) * t429 - t415 * t466 + t430 * t417;
t384 = mrSges(6,2) * t396 - mrSges(6,3) * t392 + Ifges(6,1) * t421 + Ifges(6,4) * t422 + Ifges(6,5) * t429 + t415 * t465 - t430 * t416;
t457 = pkin(10) * t381 + t452 * t383 + t448 * t384;
t364 = -mrSges(5,1) * t425 + mrSges(5,3) * t400 + t435 * Ifges(5,5) + Ifges(5,6) * t436 - pkin(4) * t376 - t443 * t382 + t457 * t446;
t365 = mrSges(5,2) * t425 - mrSges(5,3) * t399 + Ifges(5,5) * t436 - t435 * Ifges(5,6) - t448 * t383 + t452 * t384 + (-t376 * t443 - t377 * t446) * pkin(10);
t346 = -mrSges(4,1) * t425 + mrSges(4,3) * t407 + t438 * Ifges(4,5) + Ifges(4,6) * t439 - pkin(3) * t461 + pkin(9) * t462 + t453 * t364 + t449 * t365;
t351 = mrSges(4,2) * t425 - mrSges(4,3) * t406 + Ifges(4,5) * t439 - t438 * Ifges(4,6) - pkin(9) * t371 - t449 * t364 + t453 * t365;
t343 = -mrSges(3,1) * t425 + mrSges(3,3) * t414 + t456 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t460 + pkin(8) * t463 + t454 * t346 + t450 * t351;
t344 = mrSges(3,2) * t425 - mrSges(3,3) * t413 + Ifges(3,5) * qJDD(2) - t456 * Ifges(3,6) - pkin(8) * t363 - t450 * t346 + t454 * t351;
t458 = pkin(7) * t355 + t343 * t455 + t344 * t451;
t345 = mrSges(3,1) * t413 + mrSges(4,1) * t406 + mrSges(5,1) * t399 - mrSges(3,2) * t414 - mrSges(4,2) * t407 - mrSges(5,2) * t400 + Ifges(3,3) * qJDD(2) + Ifges(4,3) * t439 + Ifges(5,3) * t436 + pkin(2) * t363 + pkin(3) * t371 + pkin(4) * t377 + t446 * t382 + t457 * t443;
t342 = mrSges(2,2) * t441 - mrSges(2,3) * t431 - t451 * t343 + t455 * t344 + (-t349 * t444 - t350 * t447) * pkin(7);
t341 = -mrSges(2,1) * t441 + mrSges(2,3) * t432 - pkin(1) * t349 - t444 * t345 + t458 * t447;
t1 = [-m(1) * g(1) + t464; -m(1) * g(2) + t468; -m(1) * g(3) + m(2) * t441 + t349; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t468 - t442 * t341 + t445 * t342; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t464 + t445 * t341 + t442 * t342; -mrSges(1,1) * g(2) + mrSges(2,1) * t431 + mrSges(1,2) * g(1) - mrSges(2,2) * t432 + pkin(1) * t350 + t447 * t345 + t458 * t444;];
tauB = t1;
