% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:27
% EndTime: 2019-12-31 21:13:30
% DurationCPUTime: 2.90s
% Computational Cost: add. (27619->277), mult. (61955->356), div. (0->0), fcn. (44009->10), ass. (0->111)
t443 = 2 * qJD(4);
t425 = qJD(1) ^ 2;
t442 = pkin(2) * t425;
t420 = sin(qJ(1));
t424 = cos(qJ(1));
t433 = -t424 * g(1) - t420 * g(2);
t400 = -t425 * pkin(1) + qJDD(1) * pkin(6) + t433;
t419 = sin(qJ(2));
t441 = t419 * t400;
t423 = cos(qJ(2));
t438 = qJD(1) * qJD(2);
t403 = t419 * qJDD(1) + t423 * t438;
t367 = qJDD(2) * pkin(2) - t403 * pkin(7) - t441 + (pkin(7) * t438 + t419 * t442 - g(3)) * t423;
t388 = -t419 * g(3) + t423 * t400;
t404 = t423 * qJDD(1) - t419 * t438;
t440 = qJD(1) * t419;
t407 = qJD(2) * pkin(2) - pkin(7) * t440;
t414 = t423 ^ 2;
t368 = t404 * pkin(7) - qJD(2) * t407 - t414 * t442 + t388;
t418 = sin(qJ(3));
t422 = cos(qJ(3));
t346 = t422 * t367 - t418 * t368;
t397 = (-t419 * t418 + t423 * t422) * qJD(1);
t374 = t397 * qJD(3) + t422 * t403 + t418 * t404;
t398 = (t423 * t418 + t419 * t422) * qJD(1);
t412 = qJDD(2) + qJDD(3);
t413 = qJD(2) + qJD(3);
t332 = (t397 * t413 - t374) * qJ(4) + (t397 * t398 + t412) * pkin(3) + t346;
t347 = t418 * t367 + t422 * t368;
t373 = -t398 * qJD(3) - t418 * t403 + t422 * t404;
t390 = t413 * pkin(3) - t398 * qJ(4);
t393 = t397 ^ 2;
t334 = -t393 * pkin(3) + t373 * qJ(4) - t413 * t390 + t347;
t415 = sin(pkin(9));
t416 = cos(pkin(9));
t384 = t416 * t397 - t415 * t398;
t329 = t415 * t332 + t416 * t334 + t384 * t443;
t352 = t416 * t373 - t415 * t374;
t385 = t415 * t397 + t416 * t398;
t361 = -t384 * mrSges(5,1) + t385 * mrSges(5,2);
t377 = t413 * mrSges(5,1) - t385 * mrSges(5,3);
t362 = -t384 * pkin(4) - t385 * pkin(8);
t411 = t413 ^ 2;
t326 = -t411 * pkin(4) + t412 * pkin(8) + t384 * t362 + t329;
t437 = t420 * g(1) - t424 * g(2);
t431 = -qJDD(1) * pkin(1) - t437;
t375 = -t404 * pkin(2) + t407 * t440 + (-pkin(7) * t414 - pkin(6)) * t425 + t431;
t339 = -t373 * pkin(3) - t393 * qJ(4) + t398 * t390 + qJDD(4) + t375;
t353 = t415 * t373 + t416 * t374;
t330 = (-t384 * t413 - t353) * pkin(8) + (t385 * t413 - t352) * pkin(4) + t339;
t417 = sin(qJ(5));
t421 = cos(qJ(5));
t323 = -t417 * t326 + t421 * t330;
t371 = -t417 * t385 + t421 * t413;
t337 = t371 * qJD(5) + t421 * t353 + t417 * t412;
t351 = qJDD(5) - t352;
t372 = t421 * t385 + t417 * t413;
t354 = -t371 * mrSges(6,1) + t372 * mrSges(6,2);
t379 = qJD(5) - t384;
t355 = -t379 * mrSges(6,2) + t371 * mrSges(6,3);
t320 = m(6) * t323 + t351 * mrSges(6,1) - t337 * mrSges(6,3) - t372 * t354 + t379 * t355;
t324 = t421 * t326 + t417 * t330;
t336 = -t372 * qJD(5) - t417 * t353 + t421 * t412;
t356 = t379 * mrSges(6,1) - t372 * mrSges(6,3);
t321 = m(6) * t324 - t351 * mrSges(6,2) + t336 * mrSges(6,3) + t371 * t354 - t379 * t356;
t434 = -t417 * t320 + t421 * t321;
t307 = m(5) * t329 - t412 * mrSges(5,2) + t352 * mrSges(5,3) + t384 * t361 - t413 * t377 + t434;
t432 = -t416 * t332 + t415 * t334;
t328 = -0.2e1 * qJD(4) * t385 - t432;
t376 = -t413 * mrSges(5,2) + t384 * mrSges(5,3);
t325 = -t412 * pkin(4) - t411 * pkin(8) + (t443 + t362) * t385 + t432;
t430 = -m(6) * t325 + t336 * mrSges(6,1) - t337 * mrSges(6,2) + t371 * t355 - t372 * t356;
t316 = m(5) * t328 + t412 * mrSges(5,1) - t353 * mrSges(5,3) - t385 * t361 + t413 * t376 + t430;
t304 = t415 * t307 + t416 * t316;
t386 = -t397 * mrSges(4,1) + t398 * mrSges(4,2);
t389 = -t413 * mrSges(4,2) + t397 * mrSges(4,3);
t301 = m(4) * t346 + t412 * mrSges(4,1) - t374 * mrSges(4,3) - t398 * t386 + t413 * t389 + t304;
t391 = t413 * mrSges(4,1) - t398 * mrSges(4,3);
t435 = t416 * t307 - t415 * t316;
t302 = m(4) * t347 - t412 * mrSges(4,2) + t373 * mrSges(4,3) + t397 * t386 - t413 * t391 + t435;
t295 = t422 * t301 + t418 * t302;
t310 = t421 * t320 + t417 * t321;
t439 = qJD(1) * t423;
t436 = -t418 * t301 + t422 * t302;
t429 = -m(5) * t339 + t352 * mrSges(5,1) - t353 * mrSges(5,2) + t384 * t376 - t385 * t377 - t310;
t341 = Ifges(6,4) * t372 + Ifges(6,2) * t371 + Ifges(6,6) * t379;
t342 = Ifges(6,1) * t372 + Ifges(6,4) * t371 + Ifges(6,5) * t379;
t428 = mrSges(6,1) * t323 - mrSges(6,2) * t324 + Ifges(6,5) * t337 + Ifges(6,6) * t336 + Ifges(6,3) * t351 + t372 * t341 - t371 * t342;
t427 = m(4) * t375 - t373 * mrSges(4,1) + t374 * mrSges(4,2) - t397 * t389 + t398 * t391 - t429;
t340 = Ifges(6,5) * t372 + Ifges(6,6) * t371 + Ifges(6,3) * t379;
t313 = -mrSges(6,1) * t325 + mrSges(6,3) * t324 + Ifges(6,4) * t337 + Ifges(6,2) * t336 + Ifges(6,6) * t351 - t372 * t340 + t379 * t342;
t314 = mrSges(6,2) * t325 - mrSges(6,3) * t323 + Ifges(6,1) * t337 + Ifges(6,4) * t336 + Ifges(6,5) * t351 + t371 * t340 - t379 * t341;
t358 = Ifges(5,4) * t385 + Ifges(5,2) * t384 + Ifges(5,6) * t413;
t359 = Ifges(5,1) * t385 + Ifges(5,4) * t384 + Ifges(5,5) * t413;
t381 = Ifges(4,4) * t398 + Ifges(4,2) * t397 + Ifges(4,6) * t413;
t382 = Ifges(4,1) * t398 + Ifges(4,4) * t397 + Ifges(4,5) * t413;
t426 = mrSges(4,1) * t346 + mrSges(5,1) * t328 - mrSges(4,2) * t347 - mrSges(5,2) * t329 + pkin(3) * t304 + pkin(4) * t430 + pkin(8) * t434 + t421 * t313 + t417 * t314 + t385 * t358 - t384 * t359 - t397 * t382 + Ifges(5,6) * t352 + Ifges(5,5) * t353 + t398 * t381 + Ifges(4,6) * t373 + Ifges(4,5) * t374 + (Ifges(5,3) + Ifges(4,3)) * t412;
t406 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t439;
t405 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t440;
t402 = (-t423 * mrSges(3,1) + t419 * mrSges(3,2)) * qJD(1);
t399 = -t425 * pkin(6) + t431;
t396 = Ifges(3,5) * qJD(2) + (t419 * Ifges(3,1) + t423 * Ifges(3,4)) * qJD(1);
t395 = Ifges(3,6) * qJD(2) + (t419 * Ifges(3,4) + t423 * Ifges(3,2)) * qJD(1);
t387 = -t423 * g(3) - t441;
t380 = Ifges(4,5) * t398 + Ifges(4,6) * t397 + Ifges(4,3) * t413;
t357 = Ifges(5,5) * t385 + Ifges(5,6) * t384 + Ifges(5,3) * t413;
t297 = -mrSges(5,1) * t339 + mrSges(5,3) * t329 + Ifges(5,4) * t353 + Ifges(5,2) * t352 + Ifges(5,6) * t412 - pkin(4) * t310 - t385 * t357 + t413 * t359 - t428;
t296 = mrSges(5,2) * t339 - mrSges(5,3) * t328 + Ifges(5,1) * t353 + Ifges(5,4) * t352 + Ifges(5,5) * t412 - pkin(8) * t310 - t417 * t313 + t421 * t314 + t384 * t357 - t413 * t358;
t294 = mrSges(4,2) * t375 - mrSges(4,3) * t346 + Ifges(4,1) * t374 + Ifges(4,4) * t373 + Ifges(4,5) * t412 - qJ(4) * t304 + t416 * t296 - t415 * t297 + t397 * t380 - t413 * t381;
t293 = -mrSges(4,1) * t375 + mrSges(4,3) * t347 + Ifges(4,4) * t374 + Ifges(4,2) * t373 + Ifges(4,6) * t412 + pkin(3) * t429 + qJ(4) * t435 + t415 * t296 + t416 * t297 - t398 * t380 + t413 * t382;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t437 - mrSges(2,2) * t433 + t419 * (mrSges(3,2) * t399 - mrSges(3,3) * t387 + Ifges(3,1) * t403 + Ifges(3,4) * t404 + Ifges(3,5) * qJDD(2) - pkin(7) * t295 - qJD(2) * t395 - t418 * t293 + t422 * t294) + t423 * (-mrSges(3,1) * t399 + mrSges(3,3) * t388 + Ifges(3,4) * t403 + Ifges(3,2) * t404 + Ifges(3,6) * qJDD(2) - pkin(2) * t427 + pkin(7) * t436 + qJD(2) * t396 + t422 * t293 + t418 * t294) + pkin(1) * ((-t405 * t419 + t406 * t423) * qJD(1) - t427 - m(3) * t399 - t403 * mrSges(3,2) + t404 * mrSges(3,1)) + pkin(6) * (t423 * (m(3) * t388 - qJDD(2) * mrSges(3,2) + t404 * mrSges(3,3) - qJD(2) * t405 + t402 * t439 + t436) - t419 * (m(3) * t387 + qJDD(2) * mrSges(3,1) - t403 * mrSges(3,3) + qJD(2) * t406 - t402 * t440 + t295)); t426 + Ifges(3,3) * qJDD(2) + Ifges(3,5) * t403 + Ifges(3,6) * t404 + mrSges(3,1) * t387 - mrSges(3,2) * t388 + pkin(2) * t295 + (t419 * t395 - t423 * t396) * qJD(1); t426; -t429; t428;];
tauJ = t1;
