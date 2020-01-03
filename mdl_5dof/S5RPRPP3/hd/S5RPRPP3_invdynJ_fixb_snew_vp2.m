% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:27
% EndTime: 2019-12-31 18:12:29
% DurationCPUTime: 1.02s
% Computational Cost: add. (3307->206), mult. (7988->240), div. (0->0), fcn. (5051->6), ass. (0->90)
t445 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t444 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t424 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t443 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t423 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t442 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t402 = qJD(1) ^ 2;
t399 = sin(qJ(1));
t400 = cos(qJ(1));
t412 = -t400 * g(1) - t399 * g(2);
t383 = -t402 * pkin(1) + qJDD(1) * qJ(2) + t412;
t396 = sin(pkin(7));
t397 = cos(pkin(7));
t426 = qJD(1) * qJD(2);
t415 = -t397 * g(3) - 0.2e1 * t396 * t426;
t434 = pkin(6) * qJDD(1);
t437 = pkin(2) * t402;
t341 = (t397 * t437 - t383 - t434) * t396 + t415;
t365 = -t396 * g(3) + (t383 + 0.2e1 * t426) * t397;
t392 = t397 ^ 2;
t342 = -t392 * t437 + t397 * t434 + t365;
t398 = sin(qJ(3));
t439 = cos(qJ(3));
t334 = t439 * t341 - t398 * t342;
t417 = t397 * t439;
t429 = qJD(1) * t396;
t381 = -qJD(1) * t417 + t398 * t429;
t407 = t439 * t396 + t397 * t398;
t382 = t407 * qJD(1);
t355 = t381 * pkin(3) - t382 * qJ(4);
t401 = qJD(3) ^ 2;
t333 = -qJDD(3) * pkin(3) - t401 * qJ(4) + t382 * t355 + qJDD(4) - t334;
t357 = -t381 * mrSges(5,2) - t382 * mrSges(5,3);
t428 = t381 * qJD(3);
t363 = t407 * qJDD(1) - t428;
t441 = -m(5) * t333 - t363 * mrSges(5,1) - t382 * t357;
t440 = -2 * qJD(4);
t436 = mrSges(4,1) - mrSges(5,2);
t435 = -mrSges(6,1) - mrSges(4,3);
t354 = -t382 * mrSges(6,2) + t381 * mrSges(6,3);
t356 = t381 * mrSges(4,1) + t382 * mrSges(4,2);
t327 = -0.2e1 * qJD(5) * qJD(3) + (t381 * t382 - qJDD(3)) * qJ(5) + (t363 + t428) * pkin(4) + t333;
t374 = -t381 * mrSges(6,1) + qJD(3) * mrSges(6,2);
t413 = -m(6) * t327 + qJDD(3) * mrSges(6,3) + qJD(3) * t374;
t373 = t381 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t431 = -qJD(3) * mrSges(4,2) - t381 * mrSges(4,3) - t373;
t317 = m(4) * t334 + (-t354 - t356) * t382 + t435 * t363 + t436 * qJDD(3) + t431 * qJD(3) + t413 + t441;
t335 = t398 * t341 + t439 * t342;
t427 = t382 * qJD(3);
t362 = t427 + (t396 * t398 - t417) * qJDD(1);
t370 = qJD(3) * mrSges(4,1) - t382 * mrSges(4,3);
t405 = -t401 * pkin(3) + qJDD(3) * qJ(4) - t381 * t355 + t335;
t332 = qJD(3) * t440 - t405;
t375 = t382 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t371 = t382 * pkin(4) - qJD(3) * qJ(5);
t380 = t381 ^ 2;
t329 = -t362 * pkin(4) - t380 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t371) * qJD(3) + t405;
t372 = t382 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t421 = m(6) * t329 + qJDD(3) * mrSges(6,2) + qJD(3) * t372;
t406 = -m(5) * t332 + qJDD(3) * mrSges(5,3) + qJD(3) * t375 + t421;
t432 = -t354 - t357;
t320 = m(4) * t335 - qJDD(3) * mrSges(4,2) - qJD(3) * t370 + (-t356 + t432) * t381 + (-mrSges(5,1) + t435) * t362 + t406;
t433 = t439 * t317 + t398 * t320;
t430 = -t396 ^ 2 - t392;
t416 = t399 * g(1) - t400 * g(2);
t411 = qJDD(2) - t416;
t361 = (-pkin(2) * t397 - pkin(1)) * qJDD(1) + (t430 * pkin(6) - qJ(2)) * t402 + t411;
t403 = pkin(3) * t427 + t382 * t440 + (-t363 + t428) * qJ(4) + t361;
t326 = -t380 * pkin(4) + 0.2e1 * qJD(5) * t381 - t382 * t371 + (pkin(3) + qJ(5)) * t362 + t403;
t422 = m(6) * t326 + t362 * mrSges(6,3) + t381 * t374;
t420 = -t442 * qJD(3) + t423 * t381 - t424 * t382;
t419 = t423 * qJD(3) + t443 * t381 + t444 * t382;
t418 = t424 * qJD(3) - t444 * t381 + t445 * t382;
t414 = -t398 * t317 + t439 * t320;
t410 = -t397 * mrSges(3,1) + t396 * mrSges(3,2);
t409 = mrSges(3,3) * qJDD(1) + t402 * t410;
t331 = t362 * pkin(3) + t403;
t408 = m(5) * t331 - t363 * mrSges(5,3) - t382 * t375 + t422;
t323 = t363 * mrSges(6,1) + t382 * t354 - t413;
t404 = m(4) * t361 + (t370 - t372) * t382 + t431 * t381 + (mrSges(4,2) - mrSges(6,2)) * t363 + t436 * t362 + t408;
t385 = (Ifges(3,5) * t396 + Ifges(3,6) * t397) * qJD(1);
t379 = -qJDD(1) * pkin(1) - t402 * qJ(2) + t411;
t364 = -t396 * t383 + t415;
t324 = -t362 * mrSges(6,1) - t381 * t354 + t421;
t322 = qJDD(3) * mrSges(5,2) + qJD(3) * t373 + t323 - t441;
t321 = -t362 * mrSges(5,2) - t363 * mrSges(6,2) - t382 * t372 - t381 * t373 + t408;
t315 = t430 * t402 * mrSges(3,3) + m(3) * t379 + t410 * qJDD(1) + t404;
t314 = mrSges(5,1) * t333 + mrSges(6,1) * t327 + mrSges(4,2) * t361 - mrSges(6,2) * t326 - mrSges(4,3) * t334 - mrSges(5,3) * t331 + pkin(4) * t323 - qJ(4) * t321 - t419 * qJD(3) + t424 * qJDD(3) - t362 * t444 + t445 * t363 + t420 * t381;
t313 = -mrSges(4,1) * t361 + mrSges(4,3) * t335 - mrSges(5,1) * t332 + mrSges(5,2) * t331 + mrSges(6,1) * t329 - mrSges(6,3) * t326 + pkin(4) * t324 - qJ(5) * t422 - pkin(3) * t321 + (qJ(5) * t372 + t420) * t382 + (qJ(5) * mrSges(6,2) + t444) * t363 + t443 * t362 + t423 * qJDD(3) + t418 * qJD(3);
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t416 - mrSges(2,2) * t412 + t396 * (t397 * qJD(1) * t385 + mrSges(3,2) * t379 - mrSges(3,3) * t364 + t439 * t314 - t398 * t313 - pkin(6) * t433 + (Ifges(3,1) * t396 + Ifges(3,4) * t397) * qJDD(1)) + t397 * (-t385 * t429 - mrSges(3,1) * t379 + mrSges(3,3) * t365 + t398 * t314 + t439 * t313 - pkin(2) * t404 + pkin(6) * t414 + (Ifges(3,4) * t396 + Ifges(3,2) * t397) * qJDD(1)) - pkin(1) * t315 + qJ(2) * ((m(3) * t365 + t409 * t397 + t414) * t397 + (-m(3) * t364 + t409 * t396 - t433) * t396); t315; mrSges(4,1) * t334 - mrSges(4,2) * t335 + mrSges(5,2) * t333 - mrSges(5,3) * t332 + mrSges(6,2) * t329 - mrSges(6,3) * t327 - qJ(5) * t323 - pkin(3) * t322 + qJ(4) * t406 + t419 * t382 + t424 * t363 + t442 * qJDD(3) + (qJ(4) * t432 + t418) * t381 + (qJ(4) * (-mrSges(5,1) - mrSges(6,1)) - t423) * t362; t322; t324;];
tauJ = t1;
