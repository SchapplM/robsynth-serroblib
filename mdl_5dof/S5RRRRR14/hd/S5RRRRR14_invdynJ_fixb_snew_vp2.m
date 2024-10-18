% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR14
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:56
% EndTime: 2024-09-27 18:42:59
% DurationCPUTime: 1.75s
% Computational Cost: add. (44190->230), mult. (67806->310), div. (0->0), fcn. (48435->12), ass. (0->108)
t413 = sin(pkin(5));
t447 = pkin(8) * t413;
t411 = qJD(1) + qJD(2);
t408 = t411 ^ 2;
t446 = t408 * t413 ^ 2;
t445 = t411 * t413;
t417 = sin(qJ(3));
t444 = t413 * t417;
t422 = cos(qJ(3));
t443 = t413 * t422;
t414 = cos(pkin(5));
t442 = t414 * t417;
t441 = t414 * t422;
t409 = qJDD(1) + qJDD(2);
t440 = qJD(3) * t411;
t389 = (t409 * t417 + t422 * t440) * t413;
t402 = t414 * t409 + qJDD(3);
t403 = t414 * t411 + qJD(3);
t419 = sin(qJ(1));
t424 = cos(qJ(1));
t436 = t419 * g(1) - t424 * g(2);
t397 = qJDD(1) * pkin(1) + t436;
t431 = -t424 * g(1) - t419 * g(2);
t399 = -qJD(1) ^ 2 * pkin(1) + t431;
t418 = sin(qJ(2));
t423 = cos(qJ(2));
t380 = t423 * t397 - t418 * t399;
t370 = t409 * pkin(2) + t408 * t447 + t380;
t381 = t418 * t397 + t423 * t399;
t371 = -t408 * pkin(2) + t409 * t447 + t381;
t432 = t370 * t441 - t417 * t371;
t340 = t402 * pkin(3) - t389 * pkin(9) + (pkin(3) * t417 * t446 + (pkin(9) * t403 * t411 - g(3)) * t413) * t422 + t432;
t350 = -g(3) * t444 + t370 * t442 + t422 * t371;
t438 = t411 * t444;
t388 = t403 * pkin(3) - pkin(9) * t438;
t390 = (t409 * t422 - t417 * t440) * t413;
t439 = t422 ^ 2 * t446;
t341 = -pkin(3) * t439 + t390 * pkin(9) - t403 * t388 + t350;
t416 = sin(qJ(4));
t421 = cos(qJ(4));
t327 = t421 * t340 - t416 * t341;
t383 = (-t416 * t417 + t421 * t422) * t445;
t356 = t383 * qJD(4) + t421 * t389 + t416 * t390;
t384 = (t416 * t422 + t417 * t421) * t445;
t400 = qJDD(4) + t402;
t401 = qJD(4) + t403;
t324 = (t383 * t401 - t356) * pkin(10) + (t383 * t384 + t400) * pkin(4) + t327;
t328 = t416 * t340 + t421 * t341;
t355 = -t384 * qJD(4) - t416 * t389 + t421 * t390;
t374 = t401 * pkin(4) - t384 * pkin(10);
t382 = t383 ^ 2;
t325 = -t382 * pkin(4) + t355 * pkin(10) - t401 * t374 + t328;
t415 = sin(qJ(5));
t420 = cos(qJ(5));
t322 = t420 * t324 - t415 * t325;
t363 = t420 * t383 - t415 * t384;
t335 = t363 * qJD(5) + t415 * t355 + t420 * t356;
t364 = t415 * t383 + t420 * t384;
t347 = -t363 * mrSges(6,1) + t364 * mrSges(6,2);
t398 = qJD(5) + t401;
t357 = -t398 * mrSges(6,2) + t363 * mrSges(6,3);
t396 = qJDD(5) + t400;
t319 = m(6) * t322 + t396 * mrSges(6,1) - t335 * mrSges(6,3) - t364 * t347 + t398 * t357;
t323 = t415 * t324 + t420 * t325;
t334 = -t364 * qJD(5) + t420 * t355 - t415 * t356;
t358 = t398 * mrSges(6,1) - t364 * mrSges(6,3);
t320 = m(6) * t323 - t396 * mrSges(6,2) + t334 * mrSges(6,3) + t363 * t347 - t398 * t358;
t312 = t420 * t319 + t415 * t320;
t365 = -t383 * mrSges(5,1) + t384 * mrSges(5,2);
t372 = -t401 * mrSges(5,2) + t383 * mrSges(5,3);
t309 = m(5) * t327 + t400 * mrSges(5,1) - t356 * mrSges(5,3) - t384 * t365 + t401 * t372 + t312;
t373 = t401 * mrSges(5,1) - t384 * mrSges(5,3);
t433 = -t415 * t319 + t420 * t320;
t310 = m(5) * t328 - t400 * mrSges(5,2) + t355 * mrSges(5,3) + t383 * t365 - t401 * t373 + t433;
t305 = t421 * t309 + t416 * t310;
t437 = t411 * t443;
t349 = -g(3) * t443 + t432;
t386 = -t403 * mrSges(4,2) + mrSges(4,3) * t437;
t387 = (-mrSges(4,1) * t422 + mrSges(4,2) * t417) * t445;
t303 = m(4) * t349 + t402 * mrSges(4,1) - t389 * mrSges(4,3) + t403 * t386 - t387 * t438 + t305;
t385 = t403 * mrSges(4,1) - mrSges(4,3) * t438;
t434 = -t416 * t309 + t421 * t310;
t304 = m(4) * t350 - t402 * mrSges(4,2) + t390 * mrSges(4,3) - t403 * t385 + t387 * t437 + t434;
t435 = -t417 * t303 + t422 * t304;
t366 = -t414 * g(3) - t413 * t370;
t348 = -t390 * pkin(3) - pkin(9) * t439 + t388 * t438 + t366;
t330 = -t355 * pkin(4) - t382 * pkin(10) + t384 * t374 + t348;
t429 = m(6) * t330 - t334 * mrSges(6,1) + t335 * mrSges(6,2) - t363 * t357 + t364 * t358;
t426 = m(5) * t348 - t355 * mrSges(5,1) + t356 * mrSges(5,2) - t383 * t372 + t384 * t373 + t429;
t430 = t303 * t441 + t304 * t442 - t413 * ((t385 * t417 - t386 * t422) * t445 + t426 + t389 * mrSges(4,2) - t390 * mrSges(4,1) + m(4) * t366);
t376 = Ifges(4,6) * t403 + (Ifges(4,4) * t417 + Ifges(4,2) * t422) * t445;
t377 = Ifges(4,5) * t403 + (Ifges(4,1) * t417 + Ifges(4,4) * t422) * t445;
t360 = Ifges(5,4) * t384 + Ifges(5,2) * t383 + Ifges(5,6) * t401;
t361 = Ifges(5,1) * t384 + Ifges(5,4) * t383 + Ifges(5,5) * t401;
t343 = Ifges(6,4) * t364 + Ifges(6,2) * t363 + Ifges(6,6) * t398;
t344 = Ifges(6,1) * t364 + Ifges(6,4) * t363 + Ifges(6,5) * t398;
t427 = mrSges(6,1) * t322 - mrSges(6,2) * t323 + Ifges(6,5) * t335 + Ifges(6,6) * t334 + Ifges(6,3) * t396 + t364 * t343 - t363 * t344;
t425 = mrSges(5,1) * t327 - mrSges(5,2) * t328 + Ifges(5,5) * t356 + Ifges(5,6) * t355 + Ifges(5,3) * t400 + pkin(4) * t312 + t384 * t360 - t383 * t361 + t427;
t297 = Ifges(4,3) * t402 + Ifges(4,5) * t389 + Ifges(4,6) * t390 + mrSges(4,1) * t349 - mrSges(4,2) * t350 + t425 + (t376 * t417 - t377 * t422) * t445 + pkin(3) * t305;
t342 = Ifges(6,5) * t364 + Ifges(6,6) * t363 + Ifges(6,3) * t398;
t313 = -mrSges(6,1) * t330 + mrSges(6,3) * t323 + Ifges(6,4) * t335 + Ifges(6,2) * t334 + Ifges(6,6) * t396 - t364 * t342 + t398 * t344;
t314 = mrSges(6,2) * t330 - mrSges(6,3) * t322 + Ifges(6,1) * t335 + Ifges(6,4) * t334 + Ifges(6,5) * t396 + t363 * t342 - t398 * t343;
t359 = Ifges(5,5) * t384 + Ifges(5,6) * t383 + Ifges(5,3) * t401;
t298 = -mrSges(5,1) * t348 + mrSges(5,3) * t328 + Ifges(5,4) * t356 + Ifges(5,2) * t355 + Ifges(5,6) * t400 - pkin(4) * t429 + pkin(10) * t433 + t420 * t313 + t415 * t314 - t384 * t359 + t401 * t361;
t299 = mrSges(5,2) * t348 - mrSges(5,3) * t327 + Ifges(5,1) * t356 + Ifges(5,4) * t355 + Ifges(5,5) * t400 - pkin(10) * t312 - t415 * t313 + t420 * t314 + t383 * t359 - t401 * t360;
t375 = Ifges(4,3) * t403 + (Ifges(4,5) * t417 + Ifges(4,6) * t422) * t445;
t428 = -mrSges(3,2) * t381 + (-mrSges(4,1) * t366 + mrSges(4,3) * t350 + Ifges(4,4) * t389 + Ifges(4,2) * t390 + Ifges(4,6) * t402 - pkin(3) * t426 + pkin(9) * t434 + t421 * t298 + t416 * t299 - t375 * t438 + t403 * t377) * t443 + (mrSges(4,2) * t366 - mrSges(4,3) * t349 + Ifges(4,1) * t389 + Ifges(4,4) * t390 + Ifges(4,5) * t402 - pkin(9) * t305 - t416 * t298 + t421 * t299 + t375 * t437 - t403 * t376) * t444 + pkin(2) * t430 + t435 * t447 + t414 * t297 + mrSges(3,1) * t380 + Ifges(3,3) * t409;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t436 - mrSges(2,2) * t431 + pkin(1) * (t418 * (m(3) * t381 - t408 * mrSges(3,1) - t409 * mrSges(3,2) + t435) + t423 * (m(3) * t380 + t409 * mrSges(3,1) - t408 * mrSges(3,2) + t430)) + t428; t428; t297; t425; t427;];
tauJ = t1;
