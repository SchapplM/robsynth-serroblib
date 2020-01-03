% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:39:46
% EndTime: 2019-12-31 22:39:54
% DurationCPUTime: 5.28s
% Computational Cost: add. (59085->289), mult. (126323->381), div. (0->0), fcn. (99013->12), ass. (0->126)
t429 = sin(pkin(5));
t434 = sin(qJ(2));
t439 = cos(qJ(2));
t454 = qJD(1) * qJD(2);
t420 = (-qJDD(1) * t439 + t434 * t454) * t429;
t462 = t429 * pkin(7);
t430 = cos(pkin(5));
t461 = t430 * g(3);
t441 = qJD(1) ^ 2;
t435 = sin(qJ(1));
t440 = cos(qJ(1));
t451 = t435 * g(1) - t440 * g(2);
t415 = qJDD(1) * pkin(1) + t441 * t462 + t451;
t460 = t415 * t430;
t459 = t429 * t434;
t458 = t429 * t439;
t456 = qJD(1) * t429;
t418 = (-t439 * pkin(2) - t434 * pkin(8)) * t456;
t426 = t430 * qJD(1) + qJD(2);
t424 = t426 ^ 2;
t425 = t430 * qJDD(1) + qJDD(2);
t455 = qJD(1) * t439;
t448 = -t440 * g(1) - t435 * g(2);
t416 = -t441 * pkin(1) + qJDD(1) * t462 + t448;
t457 = t439 * t416 + t434 * t460;
t370 = -t424 * pkin(2) + t425 * pkin(8) + (-g(3) * t434 + t418 * t455) * t429 + t457;
t419 = (qJDD(1) * t434 + t439 * t454) * t429;
t371 = t420 * pkin(2) - t419 * pkin(8) - t461 + (-t415 + (pkin(2) * t434 - pkin(8) * t439) * t426 * qJD(1)) * t429;
t433 = sin(qJ(3));
t438 = cos(qJ(3));
t351 = t438 * t370 + t433 * t371;
t453 = t434 * t456;
t407 = t438 * t426 - t433 * t453;
t408 = t433 * t426 + t438 * t453;
t392 = -t407 * pkin(3) - t408 * pkin(9);
t412 = qJDD(3) + t420;
t452 = t429 * t455;
t423 = qJD(3) - t452;
t421 = t423 ^ 2;
t345 = -t421 * pkin(3) + t412 * pkin(9) + t407 * t392 + t351;
t389 = -g(3) * t458 - t434 * t416 + t439 * t460;
t369 = -t425 * pkin(2) - t424 * pkin(8) + t418 * t453 - t389;
t387 = -t408 * qJD(3) - t433 * t419 + t438 * t425;
t388 = t407 * qJD(3) + t438 * t419 + t433 * t425;
t349 = (-t407 * t423 - t388) * pkin(9) + (t408 * t423 - t387) * pkin(3) + t369;
t432 = sin(qJ(4));
t437 = cos(qJ(4));
t335 = -t432 * t345 + t437 * t349;
t394 = -t432 * t408 + t437 * t423;
t359 = t394 * qJD(4) + t437 * t388 + t432 * t412;
t385 = qJDD(4) - t387;
t395 = t437 * t408 + t432 * t423;
t406 = qJD(4) - t407;
t333 = (t394 * t406 - t359) * pkin(10) + (t394 * t395 + t385) * pkin(4) + t335;
t336 = t437 * t345 + t432 * t349;
t358 = -t395 * qJD(4) - t432 * t388 + t437 * t412;
t378 = t406 * pkin(4) - t395 * pkin(10);
t393 = t394 ^ 2;
t334 = -t393 * pkin(4) + t358 * pkin(10) - t406 * t378 + t336;
t431 = sin(qJ(5));
t436 = cos(qJ(5));
t331 = t436 * t333 - t431 * t334;
t372 = t436 * t394 - t431 * t395;
t342 = t372 * qJD(5) + t431 * t358 + t436 * t359;
t373 = t431 * t394 + t436 * t395;
t356 = -t372 * mrSges(6,1) + t373 * mrSges(6,2);
t404 = qJD(5) + t406;
t360 = -t404 * mrSges(6,2) + t372 * mrSges(6,3);
t380 = qJDD(5) + t385;
t327 = m(6) * t331 + t380 * mrSges(6,1) - t342 * mrSges(6,3) - t373 * t356 + t404 * t360;
t332 = t431 * t333 + t436 * t334;
t341 = -t373 * qJD(5) + t436 * t358 - t431 * t359;
t361 = t404 * mrSges(6,1) - t373 * mrSges(6,3);
t328 = m(6) * t332 - t380 * mrSges(6,2) + t341 * mrSges(6,3) + t372 * t356 - t404 * t361;
t319 = t436 * t327 + t431 * t328;
t374 = -t394 * mrSges(5,1) + t395 * mrSges(5,2);
t376 = -t406 * mrSges(5,2) + t394 * mrSges(5,3);
t317 = m(5) * t335 + t385 * mrSges(5,1) - t359 * mrSges(5,3) - t395 * t374 + t406 * t376 + t319;
t377 = t406 * mrSges(5,1) - t395 * mrSges(5,3);
t449 = -t431 * t327 + t436 * t328;
t318 = m(5) * t336 - t385 * mrSges(5,2) + t358 * mrSges(5,3) + t394 * t374 - t406 * t377 + t449;
t315 = -t432 * t317 + t437 * t318;
t391 = -t407 * mrSges(4,1) + t408 * mrSges(4,2);
t397 = t423 * mrSges(4,1) - t408 * mrSges(4,3);
t313 = m(4) * t351 - t412 * mrSges(4,2) + t387 * mrSges(4,3) + t407 * t391 - t423 * t397 + t315;
t350 = -t433 * t370 + t438 * t371;
t344 = -t412 * pkin(3) - t421 * pkin(9) + t408 * t392 - t350;
t337 = -t358 * pkin(4) - t393 * pkin(10) + t395 * t378 + t344;
t446 = m(6) * t337 - t341 * mrSges(6,1) + t342 * mrSges(6,2) - t372 * t360 + t373 * t361;
t329 = -m(5) * t344 + t358 * mrSges(5,1) - t359 * mrSges(5,2) + t394 * t376 - t395 * t377 - t446;
t396 = -t423 * mrSges(4,2) + t407 * mrSges(4,3);
t323 = m(4) * t350 + t412 * mrSges(4,1) - t388 * mrSges(4,3) - t408 * t391 + t423 * t396 + t329;
t307 = t433 * t313 + t438 * t323;
t450 = t438 * t313 - t433 * t323;
t314 = t437 * t317 + t432 * t318;
t353 = Ifges(6,4) * t373 + Ifges(6,2) * t372 + Ifges(6,6) * t404;
t354 = Ifges(6,1) * t373 + Ifges(6,4) * t372 + Ifges(6,5) * t404;
t445 = -mrSges(6,1) * t331 + mrSges(6,2) * t332 - Ifges(6,5) * t342 - Ifges(6,6) * t341 - Ifges(6,3) * t380 - t373 * t353 + t372 * t354;
t444 = -m(4) * t369 + t387 * mrSges(4,1) - t388 * mrSges(4,2) + t407 * t396 - t408 * t397 - t314;
t352 = Ifges(6,5) * t373 + Ifges(6,6) * t372 + Ifges(6,3) * t404;
t320 = -mrSges(6,1) * t337 + mrSges(6,3) * t332 + Ifges(6,4) * t342 + Ifges(6,2) * t341 + Ifges(6,6) * t380 - t373 * t352 + t404 * t354;
t321 = mrSges(6,2) * t337 - mrSges(6,3) * t331 + Ifges(6,1) * t342 + Ifges(6,4) * t341 + Ifges(6,5) * t380 + t372 * t352 - t404 * t353;
t362 = Ifges(5,5) * t395 + Ifges(5,6) * t394 + Ifges(5,3) * t406;
t364 = Ifges(5,1) * t395 + Ifges(5,4) * t394 + Ifges(5,5) * t406;
t308 = -mrSges(5,1) * t344 + mrSges(5,3) * t336 + Ifges(5,4) * t359 + Ifges(5,2) * t358 + Ifges(5,6) * t385 - pkin(4) * t446 + pkin(10) * t449 + t436 * t320 + t431 * t321 - t395 * t362 + t406 * t364;
t363 = Ifges(5,4) * t395 + Ifges(5,2) * t394 + Ifges(5,6) * t406;
t309 = mrSges(5,2) * t344 - mrSges(5,3) * t335 + Ifges(5,1) * t359 + Ifges(5,4) * t358 + Ifges(5,5) * t385 - pkin(10) * t319 - t431 * t320 + t436 * t321 + t394 * t362 - t406 * t363;
t382 = Ifges(4,4) * t408 + Ifges(4,2) * t407 + Ifges(4,6) * t423;
t383 = Ifges(4,1) * t408 + Ifges(4,4) * t407 + Ifges(4,5) * t423;
t443 = mrSges(4,1) * t350 - mrSges(4,2) * t351 + Ifges(4,5) * t388 + Ifges(4,6) * t387 + Ifges(4,3) * t412 + pkin(3) * t329 + pkin(9) * t315 + t437 * t308 + t432 * t309 + t408 * t382 - t407 * t383;
t442 = mrSges(5,1) * t335 - mrSges(5,2) * t336 + Ifges(5,5) * t359 + Ifges(5,6) * t358 + Ifges(5,3) * t385 + pkin(4) * t319 + t395 * t363 - t394 * t364 - t445;
t417 = (-t439 * mrSges(3,1) + t434 * mrSges(3,2)) * t456;
t414 = -t426 * mrSges(3,2) + mrSges(3,3) * t452;
t413 = t426 * mrSges(3,1) - mrSges(3,3) * t453;
t401 = -t429 * t415 - t461;
t400 = Ifges(3,5) * t426 + (t434 * Ifges(3,1) + t439 * Ifges(3,4)) * t456;
t399 = Ifges(3,6) * t426 + (t434 * Ifges(3,4) + t439 * Ifges(3,2)) * t456;
t398 = Ifges(3,3) * t426 + (t434 * Ifges(3,5) + t439 * Ifges(3,6)) * t456;
t390 = -g(3) * t459 + t457;
t381 = Ifges(4,5) * t408 + Ifges(4,6) * t407 + Ifges(4,3) * t423;
t310 = m(3) * t389 + t425 * mrSges(3,1) - t419 * mrSges(3,3) + t426 * t414 - t417 * t453 + t444;
t306 = m(3) * t390 - t425 * mrSges(3,2) - t420 * mrSges(3,3) - t426 * t413 + t417 * t452 + t450;
t305 = -mrSges(4,1) * t369 + mrSges(4,3) * t351 + Ifges(4,4) * t388 + Ifges(4,2) * t387 + Ifges(4,6) * t412 - pkin(3) * t314 - t408 * t381 + t423 * t383 - t442;
t304 = mrSges(4,2) * t369 - mrSges(4,3) * t350 + Ifges(4,1) * t388 + Ifges(4,4) * t387 + Ifges(4,5) * t412 - pkin(9) * t314 - t432 * t308 + t437 * t309 + t407 * t381 - t423 * t382;
t303 = Ifges(3,5) * t419 - Ifges(3,6) * t420 + Ifges(3,3) * t425 + mrSges(3,1) * t389 - mrSges(3,2) * t390 + t433 * t304 + t438 * t305 + pkin(2) * t444 + pkin(8) * t450 + (t434 * t399 - t439 * t400) * t456;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t451 - mrSges(2,2) * t448 + (mrSges(3,2) * t401 - mrSges(3,3) * t389 + Ifges(3,1) * t419 - Ifges(3,4) * t420 + Ifges(3,5) * t425 - pkin(8) * t307 + t438 * t304 - t433 * t305 + t398 * t452 - t426 * t399) * t459 + (-mrSges(3,1) * t401 + mrSges(3,3) * t390 + Ifges(3,4) * t419 - Ifges(3,2) * t420 + Ifges(3,6) * t425 - pkin(2) * t307 - t398 * t453 + t426 * t400 - t443) * t458 + t430 * t303 + pkin(1) * ((t434 * t306 + t439 * t310) * t430 + (-m(3) * t401 - t420 * mrSges(3,1) - t419 * mrSges(3,2) + (-t413 * t434 + t414 * t439) * t456 - t307) * t429) + (t439 * t306 - t434 * t310) * t462; t303; t443; t442; -t445;];
tauJ = t1;
