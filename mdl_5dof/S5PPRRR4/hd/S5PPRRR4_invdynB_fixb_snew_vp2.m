% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRR4
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:09
% EndTime: 2019-12-05 15:19:20
% DurationCPUTime: 7.16s
% Computational Cost: add. (88824->221), mult. (157504->297), div. (0->0), fcn. (122706->14), ass. (0->107)
t429 = sin(pkin(10));
t433 = cos(pkin(10));
t423 = -t433 * g(1) - t429 * g(2);
t428 = sin(pkin(11));
t432 = cos(pkin(11));
t422 = t429 * g(1) - t433 * g(2);
t427 = -g(3) + qJDD(1);
t431 = sin(pkin(5));
t435 = cos(pkin(5));
t448 = t422 * t435 + t427 * t431;
t396 = -t428 * t423 + t448 * t432;
t397 = t432 * t423 + t448 * t428;
t409 = -t431 * t422 + t435 * t427 + qJDD(2);
t441 = cos(qJ(3));
t434 = cos(pkin(6));
t438 = sin(qJ(3));
t460 = t434 * t438;
t430 = sin(pkin(6));
t461 = t430 * t438;
t390 = t396 * t460 + t441 * t397 + t409 * t461;
t443 = qJD(3) ^ 2;
t388 = -t443 * pkin(3) + qJDD(3) * pkin(8) + t390;
t392 = -t430 * t396 + t434 * t409;
t437 = sin(qJ(4));
t440 = cos(qJ(4));
t384 = t440 * t388 + t437 * t392;
t418 = (-mrSges(5,1) * t440 + mrSges(5,2) * t437) * qJD(3);
t455 = qJD(3) * qJD(4);
t454 = t437 * t455;
t421 = t440 * qJDD(3) - t454;
t457 = qJD(3) * t437;
t424 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t457;
t419 = (-pkin(4) * t440 - pkin(9) * t437) * qJD(3);
t442 = qJD(4) ^ 2;
t456 = t440 * qJD(3);
t382 = -t442 * pkin(4) + qJDD(4) * pkin(9) + t419 * t456 + t384;
t389 = -t438 * t397 + (t396 * t434 + t409 * t430) * t441;
t387 = -qJDD(3) * pkin(3) - t443 * pkin(8) - t389;
t453 = t440 * t455;
t420 = t437 * qJDD(3) + t453;
t385 = (-t420 - t453) * pkin(9) + (-t421 + t454) * pkin(4) + t387;
t436 = sin(qJ(5));
t439 = cos(qJ(5));
t379 = -t436 * t382 + t439 * t385;
t416 = t439 * qJD(4) - t436 * t457;
t404 = t416 * qJD(5) + t436 * qJDD(4) + t439 * t420;
t417 = t436 * qJD(4) + t439 * t457;
t405 = -t416 * mrSges(6,1) + t417 * mrSges(6,2);
t426 = qJD(5) - t456;
t407 = -t426 * mrSges(6,2) + t416 * mrSges(6,3);
t415 = qJDD(5) - t421;
t377 = m(6) * t379 + t415 * mrSges(6,1) - t404 * mrSges(6,3) - t417 * t405 + t426 * t407;
t380 = t439 * t382 + t436 * t385;
t403 = -t417 * qJD(5) + t439 * qJDD(4) - t436 * t420;
t408 = t426 * mrSges(6,1) - t417 * mrSges(6,3);
t378 = m(6) * t380 - t415 * mrSges(6,2) + t403 * mrSges(6,3) + t416 * t405 - t426 * t408;
t450 = -t436 * t377 + t439 * t378;
t370 = m(5) * t384 - qJDD(4) * mrSges(5,2) + t421 * mrSges(5,3) - qJD(4) * t424 + t418 * t456 + t450;
t459 = t440 * t392;
t383 = -t437 * t388 + t459;
t425 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t456;
t381 = -qJDD(4) * pkin(4) - t442 * pkin(9) - t459 + (qJD(3) * t419 + t388) * t437;
t445 = -m(6) * t381 + t403 * mrSges(6,1) - t404 * mrSges(6,2) + t416 * t407 - t417 * t408;
t375 = m(5) * t383 + qJDD(4) * mrSges(5,1) - t420 * mrSges(5,3) + qJD(4) * t425 - t418 * t457 + t445;
t451 = t440 * t370 - t437 * t375;
t361 = m(4) * t390 - t443 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t451;
t364 = t437 * t370 + t440 * t375;
t363 = m(4) * t392 + t364;
t371 = t439 * t377 + t436 * t378;
t444 = -m(5) * t387 + t421 * mrSges(5,1) - t420 * mrSges(5,2) - t424 * t457 + t425 * t456 - t371;
t367 = m(4) * t389 + qJDD(3) * mrSges(4,1) - t443 * mrSges(4,2) + t444;
t462 = t367 * t441;
t350 = t361 * t460 - t430 * t363 + t434 * t462;
t346 = m(3) * t396 + t350;
t355 = t441 * t361 - t438 * t367;
t354 = m(3) * t397 + t355;
t465 = t346 * t432 + t354 * t428;
t349 = t361 * t461 + t434 * t363 + t430 * t462;
t348 = m(3) * t409 + t349;
t336 = -t431 * t348 + t465 * t435;
t334 = m(2) * t422 + t336;
t340 = -t428 * t346 + t432 * t354;
t339 = m(2) * t423 + t340;
t458 = t433 * t334 + t429 * t339;
t335 = t435 * t348 + t465 * t431;
t452 = -t429 * t334 + t433 * t339;
t398 = Ifges(6,5) * t417 + Ifges(6,6) * t416 + Ifges(6,3) * t426;
t400 = Ifges(6,1) * t417 + Ifges(6,4) * t416 + Ifges(6,5) * t426;
t372 = -mrSges(6,1) * t381 + mrSges(6,3) * t380 + Ifges(6,4) * t404 + Ifges(6,2) * t403 + Ifges(6,6) * t415 - t417 * t398 + t426 * t400;
t399 = Ifges(6,4) * t417 + Ifges(6,2) * t416 + Ifges(6,6) * t426;
t373 = mrSges(6,2) * t381 - mrSges(6,3) * t379 + Ifges(6,1) * t404 + Ifges(6,4) * t403 + Ifges(6,5) * t415 + t416 * t398 - t426 * t399;
t410 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t437 + Ifges(5,6) * t440) * qJD(3);
t411 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t437 + Ifges(5,2) * t440) * qJD(3);
t356 = mrSges(5,2) * t387 - mrSges(5,3) * t383 + Ifges(5,1) * t420 + Ifges(5,4) * t421 + Ifges(5,5) * qJDD(4) - pkin(9) * t371 - qJD(4) * t411 - t436 * t372 + t439 * t373 + t410 * t456;
t412 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t437 + Ifges(5,4) * t440) * qJD(3);
t357 = -mrSges(5,1) * t387 - mrSges(6,1) * t379 + mrSges(6,2) * t380 + mrSges(5,3) * t384 + Ifges(5,4) * t420 - Ifges(6,5) * t404 + Ifges(5,2) * t421 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t403 - Ifges(6,3) * t415 - pkin(4) * t371 + qJD(4) * t412 - t417 * t399 + t416 * t400 - t410 * t457;
t342 = mrSges(4,2) * t392 - mrSges(4,3) * t389 + Ifges(4,5) * qJDD(3) - t443 * Ifges(4,6) - pkin(8) * t364 + t440 * t356 - t437 * t357;
t343 = Ifges(4,6) * qJDD(3) + t443 * Ifges(4,5) - mrSges(4,1) * t392 + mrSges(4,3) * t390 - Ifges(5,5) * t420 - Ifges(5,6) * t421 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t383 + mrSges(5,2) * t384 - t436 * t373 - t439 * t372 - pkin(4) * t445 - pkin(9) * t450 - pkin(3) * t364 + (-t437 * t411 + t440 * t412) * qJD(3);
t447 = pkin(7) * t355 + t342 * t438 + t343 * t441;
t341 = mrSges(4,1) * t389 - mrSges(4,2) * t390 + Ifges(4,3) * qJDD(3) + pkin(3) * t444 + pkin(8) * t451 + t437 * t356 + t440 * t357;
t331 = -mrSges(3,1) * t409 + mrSges(3,3) * t397 - pkin(2) * t349 - t430 * t341 + t447 * t434;
t332 = mrSges(3,2) * t409 - mrSges(3,3) * t396 + t441 * t342 - t438 * t343 + (-t349 * t430 - t350 * t434) * pkin(7);
t446 = qJ(2) * t340 + t331 * t432 + t332 * t428;
t330 = mrSges(3,1) * t396 - mrSges(3,2) * t397 + pkin(2) * t350 + t434 * t341 + t447 * t430;
t329 = mrSges(2,2) * t427 - mrSges(2,3) * t422 - t428 * t331 + t432 * t332 + (-t335 * t431 - t336 * t435) * qJ(2);
t328 = -mrSges(2,1) * t427 + mrSges(2,3) * t423 - pkin(1) * t335 - t431 * t330 + t446 * t435;
t1 = [-m(1) * g(1) + t452; -m(1) * g(2) + t458; -m(1) * g(3) + m(2) * t427 + t335; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t458 - t429 * t328 + t433 * t329; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t452 + t433 * t328 + t429 * t329; -mrSges(1,1) * g(2) + mrSges(2,1) * t422 + mrSges(1,2) * g(1) - mrSges(2,2) * t423 + pkin(1) * t336 + t435 * t330 + t446 * t431;];
tauB = t1;
