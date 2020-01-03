% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:24:34
% EndTime: 2019-12-31 20:24:38
% DurationCPUTime: 4.01s
% Computational Cost: add. (35921->287), mult. (94202->381), div. (0->0), fcn. (72292->12), ass. (0->123)
t463 = -2 * qJD(3);
t428 = sin(pkin(10));
t430 = cos(pkin(10));
t434 = sin(qJ(2));
t438 = cos(qJ(2));
t429 = sin(pkin(5));
t457 = qJD(1) * t429;
t411 = (t428 * t434 - t430 * t438) * t457;
t455 = qJD(1) * qJD(2);
t420 = (qJDD(1) * t434 + t438 * t455) * t429;
t431 = cos(pkin(5));
t423 = qJDD(1) * t431 + qJDD(2);
t424 = qJD(1) * t431 + qJD(2);
t440 = qJD(1) ^ 2;
t435 = sin(qJ(1));
t439 = cos(qJ(1));
t446 = -g(1) * t439 - g(2) * t435;
t462 = t429 * pkin(7);
t418 = -pkin(1) * t440 + qJDD(1) * t462 + t446;
t451 = t435 * g(1) - g(2) * t439;
t417 = qJDD(1) * pkin(1) + t440 * t462 + t451;
t461 = t417 * t431;
t447 = -t434 * t418 + t438 * t461;
t460 = t429 ^ 2 * t440;
t362 = t423 * pkin(2) - t420 * qJ(3) + (pkin(2) * t434 * t460 + (qJ(3) * qJD(1) * t424 - g(3)) * t429) * t438 + t447;
t459 = t429 * t434;
t388 = -g(3) * t459 + t438 * t418 + t434 * t461;
t452 = t434 * t457;
t414 = pkin(2) * t424 - qJ(3) * t452;
t421 = (qJDD(1) * t438 - t434 * t455) * t429;
t454 = t438 ^ 2 * t460;
t365 = -pkin(2) * t454 + qJ(3) * t421 - t414 * t424 + t388;
t412 = (t428 * t438 + t430 * t434) * t457;
t347 = t430 * t362 - t365 * t428 + t412 * t463;
t458 = t429 * t438;
t348 = t362 * t428 + t365 * t430 + t411 * t463;
t389 = mrSges(4,1) * t411 + mrSges(4,2) * t412;
t393 = -t420 * t428 + t421 * t430;
t399 = mrSges(4,1) * t424 - mrSges(4,3) * t412;
t390 = pkin(3) * t411 - pkin(8) * t412;
t422 = t424 ^ 2;
t346 = -pkin(3) * t422 + pkin(8) * t423 - t390 * t411 + t348;
t403 = -t431 * g(3) - t429 * t417;
t375 = -t421 * pkin(2) - qJ(3) * t454 + t414 * t452 + qJDD(3) + t403;
t394 = t420 * t430 + t421 * t428;
t350 = (t411 * t424 - t394) * pkin(8) + (t412 * t424 - t393) * pkin(3) + t375;
t433 = sin(qJ(4));
t437 = cos(qJ(4));
t343 = t346 * t437 + t350 * t433;
t396 = -t412 * t433 + t424 * t437;
t397 = t412 * t437 + t424 * t433;
t377 = -pkin(4) * t396 - pkin(9) * t397;
t392 = qJDD(4) - t393;
t410 = qJD(4) + t411;
t409 = t410 ^ 2;
t340 = -pkin(4) * t409 + pkin(9) * t392 + t377 * t396 + t343;
t345 = -t423 * pkin(3) - t422 * pkin(8) + t412 * t390 - t347;
t372 = -qJD(4) * t397 - t394 * t433 + t423 * t437;
t373 = qJD(4) * t396 + t394 * t437 + t423 * t433;
t341 = (-t396 * t410 - t373) * pkin(9) + (t397 * t410 - t372) * pkin(4) + t345;
t432 = sin(qJ(5));
t436 = cos(qJ(5));
t337 = -t340 * t432 + t341 * t436;
t379 = -t397 * t432 + t410 * t436;
t353 = qJD(5) * t379 + t373 * t436 + t392 * t432;
t380 = t397 * t436 + t410 * t432;
t358 = -mrSges(6,1) * t379 + mrSges(6,2) * t380;
t395 = qJD(5) - t396;
t363 = -mrSges(6,2) * t395 + mrSges(6,3) * t379;
t371 = qJDD(5) - t372;
t334 = m(6) * t337 + mrSges(6,1) * t371 - mrSges(6,3) * t353 - t358 * t380 + t363 * t395;
t338 = t340 * t436 + t341 * t432;
t352 = -qJD(5) * t380 - t373 * t432 + t392 * t436;
t364 = mrSges(6,1) * t395 - mrSges(6,3) * t380;
t335 = m(6) * t338 - mrSges(6,2) * t371 + mrSges(6,3) * t352 + t358 * t379 - t364 * t395;
t328 = -t334 * t432 + t335 * t436;
t376 = -mrSges(5,1) * t396 + mrSges(5,2) * t397;
t382 = mrSges(5,1) * t410 - mrSges(5,3) * t397;
t326 = m(5) * t343 - mrSges(5,2) * t392 + mrSges(5,3) * t372 + t376 * t396 - t382 * t410 + t328;
t342 = -t346 * t433 + t350 * t437;
t339 = -pkin(4) * t392 - pkin(9) * t409 + t377 * t397 - t342;
t336 = -m(6) * t339 + mrSges(6,1) * t352 - mrSges(6,2) * t353 + t363 * t379 - t364 * t380;
t381 = -mrSges(5,2) * t410 + mrSges(5,3) * t396;
t332 = m(5) * t342 + mrSges(5,1) * t392 - mrSges(5,3) * t373 - t376 * t397 + t381 * t410 + t336;
t449 = t326 * t437 - t332 * t433;
t319 = m(4) * t348 - mrSges(4,2) * t423 + mrSges(4,3) * t393 - t389 * t411 - t399 * t424 + t449;
t398 = -mrSges(4,2) * t424 - mrSges(4,3) * t411;
t327 = t334 * t436 + t335 * t432;
t443 = -m(5) * t345 + mrSges(5,1) * t372 - mrSges(5,2) * t373 + t381 * t396 - t382 * t397 - t327;
t323 = m(4) * t347 + mrSges(4,1) * t423 - mrSges(4,3) * t394 - t389 * t412 + t398 * t424 + t443;
t314 = t319 * t428 + t323 * t430;
t321 = t326 * t433 + t332 * t437;
t453 = t438 * t457;
t450 = t319 * t430 - t323 * t428;
t444 = -m(4) * t375 + t393 * mrSges(4,1) - t394 * mrSges(4,2) - t411 * t398 - t412 * t399 - t321;
t355 = Ifges(6,4) * t380 + Ifges(6,2) * t379 + Ifges(6,6) * t395;
t356 = Ifges(6,1) * t380 + Ifges(6,4) * t379 + Ifges(6,5) * t395;
t442 = mrSges(6,1) * t337 - mrSges(6,2) * t338 + Ifges(6,5) * t353 + Ifges(6,6) * t352 + Ifges(6,3) * t371 + t355 * t380 - t356 * t379;
t354 = Ifges(6,5) * t380 + Ifges(6,6) * t379 + Ifges(6,3) * t395;
t329 = -mrSges(6,1) * t339 + mrSges(6,3) * t338 + Ifges(6,4) * t353 + Ifges(6,2) * t352 + Ifges(6,6) * t371 - t354 * t380 + t356 * t395;
t330 = mrSges(6,2) * t339 - mrSges(6,3) * t337 + Ifges(6,1) * t353 + Ifges(6,4) * t352 + Ifges(6,5) * t371 + t354 * t379 - t355 * t395;
t367 = Ifges(5,4) * t397 + Ifges(5,2) * t396 + Ifges(5,6) * t410;
t368 = Ifges(5,1) * t397 + Ifges(5,4) * t396 + Ifges(5,5) * t410;
t441 = mrSges(5,1) * t342 - mrSges(5,2) * t343 + Ifges(5,5) * t373 + Ifges(5,6) * t372 + Ifges(5,3) * t392 + pkin(4) * t336 + pkin(9) * t328 + t329 * t436 + t330 * t432 + t367 * t397 - t368 * t396;
t419 = (-mrSges(3,1) * t438 + mrSges(3,2) * t434) * t457;
t416 = -mrSges(3,2) * t424 + mrSges(3,3) * t453;
t415 = mrSges(3,1) * t424 - mrSges(3,3) * t452;
t402 = Ifges(3,5) * t424 + (Ifges(3,1) * t434 + Ifges(3,4) * t438) * t457;
t401 = Ifges(3,6) * t424 + (Ifges(3,4) * t434 + Ifges(3,2) * t438) * t457;
t400 = Ifges(3,3) * t424 + (Ifges(3,5) * t434 + Ifges(3,6) * t438) * t457;
t387 = -g(3) * t458 + t447;
t385 = Ifges(4,1) * t412 - Ifges(4,4) * t411 + Ifges(4,5) * t424;
t384 = Ifges(4,4) * t412 - Ifges(4,2) * t411 + Ifges(4,6) * t424;
t383 = Ifges(4,5) * t412 - Ifges(4,6) * t411 + Ifges(4,3) * t424;
t366 = Ifges(5,5) * t397 + Ifges(5,6) * t396 + Ifges(5,3) * t410;
t316 = -mrSges(5,1) * t345 + mrSges(5,3) * t343 + Ifges(5,4) * t373 + Ifges(5,2) * t372 + Ifges(5,6) * t392 - pkin(4) * t327 - t366 * t397 + t368 * t410 - t442;
t315 = mrSges(5,2) * t345 - mrSges(5,3) * t342 + Ifges(5,1) * t373 + Ifges(5,4) * t372 + Ifges(5,5) * t392 - pkin(9) * t327 - t329 * t432 + t330 * t436 + t366 * t396 - t367 * t410;
t313 = m(3) * t388 - mrSges(3,2) * t423 + mrSges(3,3) * t421 - t415 * t424 + t419 * t453 + t450;
t312 = m(3) * t387 + mrSges(3,1) * t423 - mrSges(3,3) * t420 + t416 * t424 - t419 * t452 + t314;
t311 = -mrSges(4,1) * t375 + mrSges(4,3) * t348 + Ifges(4,4) * t394 + Ifges(4,2) * t393 + Ifges(4,6) * t423 - pkin(3) * t321 - t383 * t412 + t385 * t424 - t441;
t310 = mrSges(4,2) * t375 - mrSges(4,3) * t347 + Ifges(4,1) * t394 + Ifges(4,4) * t393 + Ifges(4,5) * t423 - pkin(8) * t321 + t315 * t437 - t316 * t433 - t383 * t411 - t384 * t424;
t309 = Ifges(3,5) * t420 + Ifges(3,6) * t421 + mrSges(3,1) * t387 - mrSges(3,2) * t388 + Ifges(4,5) * t394 + Ifges(4,6) * t393 + t412 * t384 + t411 * t385 + mrSges(4,1) * t347 - mrSges(4,2) * t348 + t433 * t315 + t437 * t316 + pkin(3) * t443 + pkin(8) * t449 + pkin(2) * t314 + (Ifges(3,3) + Ifges(4,3)) * t423 + (t401 * t434 - t402 * t438) * t457;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t451 - mrSges(2,2) * t446 + (mrSges(3,2) * t403 - mrSges(3,3) * t387 + Ifges(3,1) * t420 + Ifges(3,4) * t421 + Ifges(3,5) * t423 - qJ(3) * t314 + t310 * t430 - t311 * t428 + t400 * t453 - t401 * t424) * t459 + (-mrSges(3,1) * t403 + mrSges(3,3) * t388 + Ifges(3,4) * t420 + Ifges(3,2) * t421 + Ifges(3,6) * t423 + pkin(2) * t444 + qJ(3) * t450 + t428 * t310 + t430 * t311 - t400 * t452 + t424 * t402) * t458 + t431 * t309 + pkin(1) * ((t312 * t438 + t313 * t434) * t431 + (-m(3) * t403 + t421 * mrSges(3,1) - t420 * mrSges(3,2) + (-t415 * t434 + t416 * t438) * t457 + t444) * t429) + (-t312 * t434 + t313 * t438) * t462; t309; -t444; t441; t442;];
tauJ = t1;
