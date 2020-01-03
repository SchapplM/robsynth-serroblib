% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:17
% EndTime: 2019-12-31 20:29:20
% DurationCPUTime: 1.69s
% Computational Cost: add. (8563->260), mult. (17618->319), div. (0->0), fcn. (10390->8), ass. (0->105)
t465 = Ifges(3,1) + Ifges(4,1);
t460 = Ifges(3,4) - Ifges(4,5);
t459 = Ifges(3,5) + Ifges(4,4);
t464 = Ifges(3,2) + Ifges(4,3);
t458 = Ifges(3,6) - Ifges(4,6);
t463 = Ifges(3,3) + Ifges(4,2);
t429 = sin(qJ(4));
t430 = sin(qJ(2));
t433 = cos(qJ(4));
t434 = cos(qJ(2));
t393 = (t429 * t430 + t433 * t434) * qJD(1);
t462 = 2 * qJD(3);
t461 = mrSges(3,3) + mrSges(4,2);
t437 = qJD(1) ^ 2;
t457 = t434 ^ 2 * t437;
t450 = qJD(1) * qJD(2);
t448 = t430 * t450;
t404 = qJDD(1) * t434 - t448;
t452 = t430 * qJD(1);
t411 = -qJD(2) * pkin(3) - pkin(7) * t452;
t431 = sin(qJ(1));
t435 = cos(qJ(1));
t453 = t431 * g(1) - t435 * g(2);
t395 = -qJDD(1) * pkin(1) - t437 * pkin(6) - t453;
t449 = t434 * t450;
t403 = qJDD(1) * t430 + t449;
t444 = -t404 * pkin(2) + t395 + (-t403 - t449) * qJ(3);
t346 = -pkin(2) * t448 + pkin(3) * t404 - pkin(7) * t457 - t444 + (t411 + t462) * t452;
t394 = (-t429 * t434 + t430 * t433) * qJD(1);
t365 = -qJD(4) * t394 - t403 * t429 - t404 * t433;
t366 = -qJD(4) * t393 + t403 * t433 - t404 * t429;
t421 = -qJD(2) + qJD(4);
t341 = (t394 * t421 - t365) * pkin(4) + (t393 * t421 - t366) * pkin(8) + t346;
t446 = -g(1) * t435 - g(2) * t431;
t396 = -pkin(1) * t437 + qJDD(1) * pkin(6) + t446;
t378 = -g(3) * t430 + t434 * t396;
t400 = (-pkin(2) * t434 - qJ(3) * t430) * qJD(1);
t436 = qJD(2) ^ 2;
t451 = t434 * qJD(1);
t361 = -pkin(2) * t436 + qJDD(2) * qJ(3) + qJD(2) * t462 + t400 * t451 + t378;
t352 = -pkin(3) * t457 - pkin(7) * t404 + qJD(2) * t411 + t361;
t377 = -t434 * g(3) - t430 * t396;
t364 = -qJDD(2) * pkin(2) - qJ(3) * t436 + t400 * t452 + qJDD(3) - t377;
t353 = (-t403 + t449) * pkin(7) + (-t430 * t434 * t437 - qJDD(2)) * pkin(3) + t364;
t345 = t352 * t433 + t353 * t429;
t374 = pkin(4) * t393 - pkin(8) * t394;
t419 = t421 ^ 2;
t420 = -qJDD(2) + qJDD(4);
t343 = -pkin(4) * t419 + pkin(8) * t420 - t374 * t393 + t345;
t428 = sin(qJ(5));
t432 = cos(qJ(5));
t339 = t341 * t432 - t343 * t428;
t375 = -t394 * t428 + t421 * t432;
t349 = qJD(5) * t375 + t366 * t432 + t420 * t428;
t376 = t394 * t432 + t421 * t428;
t359 = -mrSges(6,1) * t375 + mrSges(6,2) * t376;
t363 = qJDD(5) - t365;
t386 = qJD(5) + t393;
t367 = -mrSges(6,2) * t386 + mrSges(6,3) * t375;
t336 = m(6) * t339 + mrSges(6,1) * t363 - mrSges(6,3) * t349 - t359 * t376 + t367 * t386;
t340 = t341 * t428 + t343 * t432;
t348 = -qJD(5) * t376 - t366 * t428 + t420 * t432;
t368 = mrSges(6,1) * t386 - mrSges(6,3) * t376;
t337 = m(6) * t340 - mrSges(6,2) * t363 + mrSges(6,3) * t348 + t359 * t375 - t368 * t386;
t328 = t336 * t432 + t337 * t428;
t456 = t463 * qJD(2) + (t430 * t459 + t434 * t458) * qJD(1);
t455 = -t458 * qJD(2) + (-t430 * t460 - t464 * t434) * qJD(1);
t454 = t459 * qJD(2) + (t430 * t465 + t434 * t460) * qJD(1);
t329 = -t336 * t428 + t337 * t432;
t373 = mrSges(5,1) * t393 + mrSges(5,2) * t394;
t380 = mrSges(5,1) * t421 - mrSges(5,3) * t394;
t327 = m(5) * t345 - mrSges(5,2) * t420 + mrSges(5,3) * t365 - t373 * t393 - t380 * t421 + t329;
t344 = -t352 * t429 + t353 * t433;
t342 = -pkin(4) * t420 - pkin(8) * t419 + t374 * t394 - t344;
t338 = -m(6) * t342 + mrSges(6,1) * t348 - mrSges(6,2) * t349 + t375 * t367 - t368 * t376;
t379 = -mrSges(5,2) * t421 - mrSges(5,3) * t393;
t332 = m(5) * t344 + mrSges(5,1) * t420 - mrSges(5,3) * t366 - t373 * t394 + t379 * t421 + t338;
t447 = t327 * t433 - t429 * t332;
t324 = t429 * t327 + t433 * t332;
t401 = (-mrSges(4,1) * t434 - mrSges(4,3) * t430) * qJD(1);
t408 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t452;
t443 = m(4) * t361 + qJDD(2) * mrSges(4,3) + qJD(2) * t408 + t401 * t451 + t447;
t410 = mrSges(4,2) * t451 + qJD(2) * mrSges(4,3);
t442 = m(4) * t364 - qJDD(2) * mrSges(4,1) - qJD(2) * t410 + t324;
t441 = m(5) * t346 - t365 * mrSges(5,1) + mrSges(5,2) * t366 + t393 * t379 + t380 * t394 + t328;
t358 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t452 + t444;
t440 = m(4) * t358 - t441;
t355 = Ifges(6,4) * t376 + Ifges(6,2) * t375 + Ifges(6,6) * t386;
t356 = Ifges(6,1) * t376 + Ifges(6,4) * t375 + Ifges(6,5) * t386;
t439 = mrSges(6,1) * t339 - mrSges(6,2) * t340 + Ifges(6,5) * t349 + Ifges(6,6) * t348 + Ifges(6,3) * t363 + t355 * t376 - t356 * t375;
t354 = Ifges(6,5) * t376 + Ifges(6,6) * t375 + Ifges(6,3) * t386;
t330 = -mrSges(6,1) * t342 + mrSges(6,3) * t340 + Ifges(6,4) * t349 + Ifges(6,2) * t348 + Ifges(6,6) * t363 - t354 * t376 + t356 * t386;
t331 = mrSges(6,2) * t342 - mrSges(6,3) * t339 + Ifges(6,1) * t349 + Ifges(6,4) * t348 + Ifges(6,5) * t363 + t354 * t375 - t355 * t386;
t370 = Ifges(5,4) * t394 - Ifges(5,2) * t393 + Ifges(5,6) * t421;
t371 = Ifges(5,1) * t394 - Ifges(5,4) * t393 + Ifges(5,5) * t421;
t438 = mrSges(5,1) * t344 - mrSges(5,2) * t345 + Ifges(5,5) * t366 + Ifges(5,6) * t365 + Ifges(5,3) * t420 + pkin(4) * t338 + pkin(8) * t329 + t432 * t330 + t428 * t331 + t394 * t370 + t393 * t371;
t409 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t451;
t407 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t452;
t402 = (-mrSges(3,1) * t434 + mrSges(3,2) * t430) * qJD(1);
t369 = Ifges(5,5) * t394 - Ifges(5,6) * t393 + Ifges(5,3) * t421;
t325 = -mrSges(4,1) * t404 - mrSges(4,3) * t403 + (-t408 * t430 - t410 * t434) * qJD(1) + t440;
t323 = t403 * mrSges(4,2) + t401 * t452 + t442;
t322 = -mrSges(5,1) * t346 + mrSges(5,3) * t345 + Ifges(5,4) * t366 + Ifges(5,2) * t365 + Ifges(5,6) * t420 - pkin(4) * t328 - t369 * t394 + t371 * t421 - t439;
t321 = mrSges(5,2) * t346 - mrSges(5,3) * t344 + Ifges(5,1) * t366 + Ifges(5,4) * t365 + Ifges(5,5) * t420 - pkin(8) * t328 - t330 * t428 + t331 * t432 - t369 * t393 - t370 * t421;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t453 - mrSges(2,2) * t446 + t430 * (mrSges(3,2) * t395 + mrSges(4,2) * t364 - mrSges(3,3) * t377 - mrSges(4,3) * t358 - pkin(7) * t324 - qJ(3) * t325 + t455 * qJD(2) + t459 * qJDD(2) + t433 * t321 - t429 * t322 + t465 * t403 + t460 * t404 + t456 * t451) + t434 * (-mrSges(3,1) * t395 - mrSges(4,1) * t358 + mrSges(4,2) * t361 + mrSges(3,3) * t378 - pkin(2) * t325 + pkin(3) * t441 - pkin(7) * t447 + t454 * qJD(2) + t458 * qJDD(2) - t429 * t321 - t433 * t322 + t460 * t403 + t464 * t404 - t456 * t452) + pkin(1) * (-m(3) * t395 + (mrSges(3,1) + mrSges(4,1)) * t404 + (-mrSges(3,2) + mrSges(4,3)) * t403 + ((t409 + t410) * t434 + (-t407 + t408) * t430) * qJD(1) - t440) + pkin(6) * (t434 * (m(3) * t378 - qJDD(2) * mrSges(3,2) - qJD(2) * t407 + t402 * t451 + t404 * t461 + t443) + (-m(3) * t377 - qJDD(2) * mrSges(3,1) - qJD(2) * t409 + t461 * t403 + (t401 + t402) * t452 + t442) * t430); -t438 + qJ(3) * t443 + (mrSges(4,2) * qJ(3) + t458) * t404 + t459 * t403 + t463 * qJDD(2) - pkin(3) * t324 - pkin(2) * t323 + mrSges(3,1) * t377 - mrSges(3,2) * t378 + mrSges(4,3) * t361 - mrSges(4,1) * t364 + (-t455 * t430 - t454 * t434) * qJD(1); t323; t438; t439;];
tauJ = t1;
