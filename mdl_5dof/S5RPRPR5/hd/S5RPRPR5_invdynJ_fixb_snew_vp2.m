% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:42:10
% EndTime: 2020-01-03 11:42:15
% DurationCPUTime: 2.61s
% Computational Cost: add. (15870->240), mult. (40702->322), div. (0->0), fcn. (27700->10), ass. (0->105)
t423 = qJD(1) ^ 2;
t419 = sin(qJ(1));
t422 = cos(qJ(1));
t433 = -t419 * g(2) + t422 * g(3);
t449 = -t423 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t433;
t414 = sin(pkin(8));
t416 = cos(pkin(8));
t374 = -t416 * g(1) - t449 * t414;
t375 = -t414 * g(1) + t449 * t416;
t430 = -t416 * pkin(2) - t414 * pkin(6);
t399 = t430 * qJD(1);
t441 = t416 * qJD(1);
t365 = t399 * t441 + t375;
t444 = -t422 * g(2) - t419 * g(3);
t427 = -t423 * qJ(2) + qJDD(2) - t444;
t376 = (-pkin(1) + t430) * qJDD(1) + t427;
t421 = cos(qJ(3));
t373 = t421 * t376;
t418 = sin(qJ(3));
t439 = qJD(1) * qJD(3);
t393 = (qJDD(1) * t421 - t418 * t439) * t414;
t438 = t416 * qJDD(1);
t403 = qJDD(3) - t438;
t404 = qJD(3) - t441;
t443 = qJD(1) * t414;
t411 = t414 ^ 2;
t446 = t411 * t423;
t346 = t403 * pkin(3) - t393 * qJ(4) + t373 + (-pkin(3) * t421 * t446 - qJ(4) * t404 * t443 - t365) * t418;
t355 = t421 * t365 + t418 * t376;
t435 = t421 * t443;
t389 = t404 * pkin(3) - qJ(4) * t435;
t392 = (-qJDD(1) * t418 - t421 * t439) * t414;
t437 = t418 ^ 2 * t446;
t347 = -pkin(3) * t437 + t392 * qJ(4) - t404 * t389 + t355;
t413 = sin(pkin(9));
t415 = cos(pkin(9));
t385 = (-t418 * t413 + t421 * t415) * t443;
t333 = -0.2e1 * qJD(4) * t385 + t415 * t346 - t413 * t347;
t368 = t413 * t392 + t415 * t393;
t384 = (-t421 * t413 - t418 * t415) * t443;
t331 = (t384 * t404 - t368) * pkin(7) + (t384 * t385 + t403) * pkin(4) + t333;
t334 = 0.2e1 * qJD(4) * t384 + t413 * t346 + t415 * t347;
t367 = t415 * t392 - t413 * t393;
t371 = t404 * pkin(4) - t385 * pkin(7);
t383 = t384 ^ 2;
t332 = -t383 * pkin(4) + t367 * pkin(7) - t404 * t371 + t334;
t417 = sin(qJ(5));
t420 = cos(qJ(5));
t329 = t420 * t331 - t417 * t332;
t362 = t420 * t384 - t417 * t385;
t342 = t362 * qJD(5) + t417 * t367 + t420 * t368;
t363 = t417 * t384 + t420 * t385;
t352 = -t362 * mrSges(6,1) + t363 * mrSges(6,2);
t402 = qJD(5) + t404;
t356 = -t402 * mrSges(6,2) + t362 * mrSges(6,3);
t401 = qJDD(5) + t403;
t325 = m(6) * t329 + t401 * mrSges(6,1) - t342 * mrSges(6,3) - t363 * t352 + t402 * t356;
t330 = t417 * t331 + t420 * t332;
t341 = -t363 * qJD(5) + t420 * t367 - t417 * t368;
t357 = t402 * mrSges(6,1) - t363 * mrSges(6,3);
t326 = m(6) * t330 - t401 * mrSges(6,2) + t341 * mrSges(6,3) + t362 * t352 - t402 * t357;
t319 = t420 * t325 + t417 * t326;
t366 = -t384 * mrSges(5,1) + t385 * mrSges(5,2);
t369 = -t404 * mrSges(5,2) + t384 * mrSges(5,3);
t317 = m(5) * t333 + t403 * mrSges(5,1) - t368 * mrSges(5,3) - t385 * t366 + t404 * t369 + t319;
t370 = t404 * mrSges(5,1) - t385 * mrSges(5,3);
t431 = -t417 * t325 + t420 * t326;
t318 = m(5) * t334 - t403 * mrSges(5,2) + t367 * mrSges(5,3) + t384 * t366 - t404 * t370 + t431;
t313 = t415 * t317 + t413 * t318;
t354 = -t418 * t365 + t373;
t359 = Ifges(5,4) * t385 + Ifges(5,2) * t384 + Ifges(5,6) * t404;
t360 = Ifges(5,1) * t385 + Ifges(5,4) * t384 + Ifges(5,5) * t404;
t349 = Ifges(6,4) * t363 + Ifges(6,2) * t362 + Ifges(6,6) * t402;
t350 = Ifges(6,1) * t363 + Ifges(6,4) * t362 + Ifges(6,5) * t402;
t425 = -mrSges(6,1) * t329 + mrSges(6,2) * t330 - Ifges(6,5) * t342 - Ifges(6,6) * t341 - Ifges(6,3) * t401 - t363 * t349 + t362 * t350;
t448 = -mrSges(4,1) * t354 - mrSges(5,1) * t333 + mrSges(4,2) * t355 + mrSges(5,2) * t334 - Ifges(4,5) * t393 - Ifges(5,5) * t368 - Ifges(4,6) * t392 - Ifges(5,6) * t367 - pkin(3) * t313 - pkin(4) * t319 - t385 * t359 + t384 * t360 - (Ifges(4,3) + Ifges(5,3)) * t403 + t425;
t442 = qJDD(1) * mrSges(3,3);
t436 = t418 * t443;
t432 = -t413 * t317 + t415 * t318;
t364 = t399 * t443 - t374;
t429 = -t416 * mrSges(3,1) + t414 * mrSges(3,2);
t388 = -t404 * mrSges(4,2) - mrSges(4,3) * t436;
t391 = (t418 * mrSges(4,1) + t421 * mrSges(4,2)) * t443;
t311 = m(4) * t354 + t403 * mrSges(4,1) - t393 * mrSges(4,3) + t404 * t388 - t391 * t435 + t313;
t390 = t404 * mrSges(4,1) - mrSges(4,3) * t435;
t312 = m(4) * t355 - t403 * mrSges(4,2) + t392 * mrSges(4,3) - t404 * t390 - t391 * t436 + t432;
t308 = t421 * t311 + t418 * t312;
t378 = Ifges(4,6) * t404 + (t421 * Ifges(4,4) - t418 * Ifges(4,2)) * t443;
t379 = Ifges(4,5) * t404 + (t421 * Ifges(4,1) - t418 * Ifges(4,4)) * t443;
t428 = t421 * t378 + t418 * t379;
t353 = -t392 * pkin(3) - qJ(4) * t437 + t389 * t435 + qJDD(4) + t364;
t336 = -t367 * pkin(4) - t383 * pkin(7) + t385 * t371 + t353;
t426 = m(6) * t336 - t341 * mrSges(6,1) + t342 * mrSges(6,2) - t362 * t356 + t363 * t357;
t327 = m(5) * t353 - t367 * mrSges(5,1) + t368 * mrSges(5,2) - t384 * t369 + t385 * t370 + t426;
t398 = (Ifges(3,5) * t414 + Ifges(3,6) * t416) * qJD(1);
t397 = t429 * qJD(1);
t395 = -qJDD(1) * pkin(1) + t427;
t358 = Ifges(5,5) * t385 + Ifges(5,6) * t384 + Ifges(5,3) * t404;
t348 = Ifges(6,5) * t363 + Ifges(6,6) * t362 + Ifges(6,3) * t402;
t321 = mrSges(6,2) * t336 - mrSges(6,3) * t329 + Ifges(6,1) * t342 + Ifges(6,4) * t341 + Ifges(6,5) * t401 + t362 * t348 - t402 * t349;
t320 = -mrSges(6,1) * t336 + mrSges(6,3) * t330 + Ifges(6,4) * t342 + Ifges(6,2) * t341 + Ifges(6,6) * t401 - t363 * t348 + t402 * t350;
t310 = mrSges(5,2) * t353 - mrSges(5,3) * t333 + Ifges(5,1) * t368 + Ifges(5,4) * t367 + Ifges(5,5) * t403 - pkin(7) * t319 - t417 * t320 + t420 * t321 + t384 * t358 - t404 * t359;
t309 = -mrSges(5,1) * t353 + mrSges(5,3) * t334 + Ifges(5,4) * t368 + Ifges(5,2) * t367 + Ifges(5,6) * t403 - pkin(4) * t426 + pkin(7) * t431 + t420 * t320 + t417 * t321 - t385 * t358 + t404 * t360;
t307 = m(3) * t395 + t429 * qJDD(1) + (-t416 ^ 2 - t411) * t423 * mrSges(3,3) + t308;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t444 - mrSges(2,2) * t433 + t414 * (mrSges(3,2) * t395 - mrSges(3,3) * t374 + t421 * (mrSges(4,2) * t364 - mrSges(4,3) * t354 + Ifges(4,1) * t393 + Ifges(4,4) * t392 + Ifges(4,5) * t403 - qJ(4) * t313 - t413 * t309 + t415 * t310 - t404 * t378) - t418 * (-mrSges(4,1) * t364 + mrSges(4,3) * t355 + Ifges(4,4) * t393 + Ifges(4,2) * t392 + Ifges(4,6) * t403 - pkin(3) * t327 + qJ(4) * t432 + t415 * t309 + t413 * t310 + t404 * t379) - pkin(6) * t308 + (Ifges(3,1) * t414 + Ifges(3,4) * t416) * qJDD(1) + t398 * t441) + t416 * (Ifges(3,2) * t438 - mrSges(3,1) * t395 + mrSges(3,3) * t375 + (Ifges(3,4) * qJDD(1) + (-t398 - t428) * qJD(1)) * t414 - pkin(2) * t308 + t448) - pkin(1) * t307 + qJ(2) * ((m(3) * t375 - t418 * t311 + t421 * t312 + (qJD(1) * t397 + t442) * t416) * t416 + (t414 * t442 + (t388 * t418 + t390 * t421 + t397) * t443 + t393 * mrSges(4,2) - t392 * mrSges(4,1) - m(3) * t374 + m(4) * t364 + t327) * t414); t307; t428 * t443 - t448; t327; -t425;];
tauJ = t1;
