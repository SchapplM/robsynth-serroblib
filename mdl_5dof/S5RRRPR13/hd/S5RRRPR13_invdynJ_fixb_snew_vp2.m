% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR13
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:43
% EndTime: 2019-12-31 21:43:49
% DurationCPUTime: 2.53s
% Computational Cost: add. (19143->270), mult. (41351->338), div. (0->0), fcn. (30819->10), ass. (0->120)
t476 = Ifges(4,1) + Ifges(5,2);
t468 = Ifges(4,4) + Ifges(5,6);
t467 = Ifges(4,5) - Ifges(5,4);
t475 = -Ifges(4,2) - Ifges(5,3);
t466 = Ifges(4,6) - Ifges(5,5);
t474 = Ifges(4,3) + Ifges(5,1);
t428 = sin(pkin(5));
t432 = sin(qJ(2));
t435 = cos(qJ(2));
t453 = qJD(1) * qJD(2);
t417 = (-qJDD(1) * t435 + t432 * t453) * t428;
t429 = cos(pkin(5));
t425 = qJD(1) * t429 + qJD(2);
t431 = sin(qJ(3));
t455 = qJD(1) * t428;
t452 = t432 * t455;
t471 = cos(qJ(3));
t404 = -t471 * t425 + t431 * t452;
t454 = qJD(1) * t435;
t451 = t428 * t454;
t421 = -qJD(3) + t451;
t391 = mrSges(5,1) * t404 + mrSges(5,3) * t421;
t409 = qJDD(3) + t417;
t415 = (-t435 * pkin(2) - t432 * pkin(8)) * t455;
t423 = t425 ^ 2;
t424 = qJDD(1) * t429 + qJDD(2);
t437 = qJD(1) ^ 2;
t433 = sin(qJ(1));
t436 = cos(qJ(1));
t448 = -g(1) * t436 - g(2) * t433;
t470 = pkin(7) * t428;
t413 = -pkin(1) * t437 + qJDD(1) * t470 + t448;
t450 = t433 * g(1) - g(2) * t436;
t412 = qJDD(1) * pkin(1) + t437 * t470 + t450;
t463 = t412 * t429;
t456 = t435 * t413 + t432 * t463;
t359 = -pkin(2) * t423 + pkin(8) * t424 + (-g(3) * t432 + t415 * t454) * t428 + t456;
t416 = (qJDD(1) * t432 + t435 * t453) * t428;
t469 = g(3) * t429;
t360 = pkin(2) * t417 - pkin(8) * t416 - t469 + (-t412 + (pkin(2) * t432 - pkin(8) * t435) * t425 * qJD(1)) * t428;
t343 = -t431 * t359 + t471 * t360;
t405 = t431 * t425 + t471 * t452;
t383 = pkin(3) * t404 - qJ(4) * t405;
t420 = t421 ^ 2;
t341 = -t409 * pkin(3) - t420 * qJ(4) + t405 * t383 + qJDD(4) - t343;
t380 = -t404 * qJD(3) + t471 * t416 + t431 * t424;
t464 = t404 * t421;
t336 = (t404 * t405 - t409) * pkin(9) + (t380 - t464) * pkin(4) + t341;
t379 = qJD(3) * t405 + t416 * t431 - t471 * t424;
t393 = pkin(4) * t405 + pkin(9) * t421;
t403 = t404 ^ 2;
t461 = t428 * t435;
t381 = -g(3) * t461 - t432 * t413 + t435 * t463;
t358 = -pkin(2) * t424 - pkin(8) * t423 + t415 * t452 - t381;
t472 = -2 * qJD(4);
t439 = (-t380 - t464) * qJ(4) + t358 + (-t421 * pkin(3) + t472) * t405;
t339 = -pkin(4) * t403 - t393 * t405 + (pkin(3) + pkin(9)) * t379 + t439;
t430 = sin(qJ(5));
t434 = cos(qJ(5));
t334 = t336 * t434 - t339 * t430;
t387 = t404 * t434 + t421 * t430;
t350 = qJD(5) * t387 + t379 * t430 + t409 * t434;
t388 = t404 * t430 - t421 * t434;
t361 = -mrSges(6,1) * t387 + mrSges(6,2) * t388;
t402 = qJD(5) + t405;
t364 = -mrSges(6,2) * t402 + mrSges(6,3) * t387;
t376 = qJDD(5) + t380;
t331 = m(6) * t334 + mrSges(6,1) * t376 - mrSges(6,3) * t350 - t361 * t388 + t364 * t402;
t335 = t336 * t430 + t339 * t434;
t349 = -qJD(5) * t388 + t379 * t434 - t409 * t430;
t365 = mrSges(6,1) * t402 - mrSges(6,3) * t388;
t332 = m(6) * t335 - mrSges(6,2) * t376 + mrSges(6,3) * t349 + t361 * t387 - t365 * t402;
t322 = t331 * t434 + t332 * t430;
t385 = -mrSges(5,2) * t404 - mrSges(5,3) * t405;
t444 = -m(5) * t341 - t380 * mrSges(5,1) - t405 * t385 - t322;
t321 = mrSges(5,2) * t409 - t391 * t421 - t444;
t344 = t471 * t359 + t431 * t360;
t443 = -pkin(3) * t420 + qJ(4) * t409 - t383 * t404 + t344;
t338 = -pkin(4) * t379 - pkin(9) * t403 + (t472 - t393) * t421 + t443;
t351 = Ifges(6,5) * t388 + Ifges(6,6) * t387 + Ifges(6,3) * t402;
t353 = Ifges(6,1) * t388 + Ifges(6,4) * t387 + Ifges(6,5) * t402;
t323 = -mrSges(6,1) * t338 + mrSges(6,3) * t335 + Ifges(6,4) * t350 + Ifges(6,2) * t349 + Ifges(6,6) * t376 - t351 * t388 + t353 * t402;
t352 = Ifges(6,4) * t388 + Ifges(6,2) * t387 + Ifges(6,6) * t402;
t324 = mrSges(6,2) * t338 - mrSges(6,3) * t334 + Ifges(6,1) * t350 + Ifges(6,4) * t349 + Ifges(6,5) * t376 + t351 * t387 - t352 * t402;
t340 = 0.2e1 * qJD(4) * t421 - t443;
t392 = mrSges(5,1) * t405 - mrSges(5,2) * t421;
t445 = -m(6) * t338 + mrSges(6,1) * t349 - t350 * mrSges(6,2) + t364 * t387 - t388 * t365;
t441 = -m(5) * t340 + t409 * mrSges(5,3) - t421 * t392 - t445;
t457 = -t468 * t404 + t476 * t405 - t467 * t421;
t458 = t475 * t404 + t468 * t405 - t466 * t421;
t473 = -t466 * t379 + t467 * t380 + t457 * t404 + t458 * t405 + t474 * t409 + mrSges(4,1) * t343 - mrSges(4,2) * t344 + mrSges(5,2) * t341 - mrSges(5,3) * t340 - pkin(3) * t321 - pkin(9) * t322 + qJ(4) * (-mrSges(5,1) * t379 - t385 * t404 + t441) - t430 * t323 + t434 * t324;
t462 = t428 * t432;
t384 = mrSges(4,1) * t404 + mrSges(4,2) * t405;
t389 = mrSges(4,2) * t421 - mrSges(4,3) * t404;
t319 = m(4) * t343 - mrSges(4,3) * t380 - t384 * t405 + (-t389 + t391) * t421 + (mrSges(4,1) - mrSges(5,2)) * t409 + t444;
t390 = -mrSges(4,1) * t421 - mrSges(4,3) * t405;
t327 = m(4) * t344 - mrSges(4,2) * t409 + t390 * t421 + (-t384 - t385) * t404 + (-mrSges(4,3) - mrSges(5,1)) * t379 + t441;
t316 = t471 * t319 + t431 * t327;
t460 = -t430 * t331 + t434 * t332;
t459 = t466 * t404 - t467 * t405 + t474 * t421;
t449 = -t319 * t431 + t471 * t327;
t342 = pkin(3) * t379 + t439;
t447 = -m(5) * t342 + t379 * mrSges(5,2) + t404 * t391 - t460;
t442 = mrSges(6,1) * t334 - mrSges(6,2) * t335 + Ifges(6,5) * t350 + Ifges(6,6) * t349 + Ifges(6,3) * t376 + t388 * t352 - t387 * t353;
t440 = -m(4) * t358 - t379 * mrSges(4,1) - t404 * t389 + (-t390 + t392) * t405 + (-mrSges(4,2) + mrSges(5,3)) * t380 + t447;
t414 = (-t435 * mrSges(3,1) + t432 * mrSges(3,2)) * t455;
t411 = -mrSges(3,2) * t425 + mrSges(3,3) * t451;
t410 = mrSges(3,1) * t425 - mrSges(3,3) * t452;
t398 = -t412 * t428 - t469;
t397 = Ifges(3,5) * t425 + (t432 * Ifges(3,1) + t435 * Ifges(3,4)) * t455;
t396 = Ifges(3,6) * t425 + (t432 * Ifges(3,4) + t435 * Ifges(3,2)) * t455;
t395 = Ifges(3,3) * t425 + (t432 * Ifges(3,5) + t435 * Ifges(3,6)) * t455;
t382 = -g(3) * t462 + t456;
t320 = -mrSges(5,3) * t380 - t392 * t405 - t447;
t317 = m(3) * t381 + mrSges(3,1) * t424 - mrSges(3,3) * t416 + t411 * t425 - t414 * t452 + t440;
t315 = m(3) * t382 - mrSges(3,2) * t424 - mrSges(3,3) * t417 - t410 * t425 + t414 * t451 + t449;
t314 = mrSges(5,1) * t341 + mrSges(4,2) * t358 - mrSges(4,3) * t343 - mrSges(5,3) * t342 + pkin(4) * t322 - qJ(4) * t320 - t468 * t379 + t476 * t380 + t459 * t404 + t467 * t409 + t458 * t421 + t442;
t313 = -mrSges(4,1) * t358 - mrSges(5,1) * t340 + mrSges(5,2) * t342 + mrSges(4,3) * t344 - pkin(3) * t320 - pkin(4) * t445 - pkin(9) * t460 - t434 * t323 - t430 * t324 + t475 * t379 + t468 * t380 + t459 * t405 + t466 * t409 - t457 * t421;
t312 = Ifges(3,5) * t416 - Ifges(3,6) * t417 + Ifges(3,3) * t424 + mrSges(3,1) * t381 - mrSges(3,2) * t382 + t431 * t314 + t471 * t313 + pkin(2) * t440 + pkin(8) * t449 + (t432 * t396 - t435 * t397) * t455;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t450 - mrSges(2,2) * t448 + (mrSges(3,2) * t398 - mrSges(3,3) * t381 + Ifges(3,1) * t416 - Ifges(3,4) * t417 + Ifges(3,5) * t424 - pkin(8) * t316 - t431 * t313 + t471 * t314 + t395 * t451 - t425 * t396) * t462 + (-mrSges(3,1) * t398 + mrSges(3,3) * t382 + Ifges(3,4) * t416 - Ifges(3,2) * t417 + Ifges(3,6) * t424 - pkin(2) * t316 - t395 * t452 + t425 * t397 - t473) * t461 + t429 * t312 + pkin(1) * ((t315 * t432 + t317 * t435) * t429 + (-m(3) * t398 - t417 * mrSges(3,1) - t416 * mrSges(3,2) + (-t410 * t432 + t411 * t435) * t455 - t316) * t428) + (t315 * t435 - t317 * t432) * t470; t312; t473; t321; t442;];
tauJ = t1;
