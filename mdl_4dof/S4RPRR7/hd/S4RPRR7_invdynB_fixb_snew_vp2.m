% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:49
% EndTime: 2019-12-31 16:53:51
% DurationCPUTime: 1.66s
% Computational Cost: add. (14230->217), mult. (33418->275), div. (0->0), fcn. (22354->8), ass. (0->95)
t414 = qJD(1) ^ 2;
t405 = sin(pkin(7));
t406 = cos(pkin(7));
t408 = sin(qJ(3));
t411 = cos(qJ(3));
t419 = t405 * t408 - t406 * t411;
t390 = t419 * qJD(1);
t420 = t405 * t411 + t406 * t408;
t391 = t420 * qJD(1);
t432 = t391 * qJD(3);
t379 = -t419 * qJDD(1) - t432;
t438 = pkin(2) * t406;
t437 = mrSges(3,2) * t405;
t404 = t406 ^ 2;
t436 = t404 * t414;
t409 = sin(qJ(1));
t412 = cos(qJ(1));
t396 = -t412 * g(1) - t409 * g(2);
t392 = -t414 * pkin(1) + qJDD(1) * qJ(2) + t396;
t431 = qJD(1) * qJD(2);
t429 = -t406 * g(3) - 0.2e1 * t405 * t431;
t368 = (-pkin(5) * qJDD(1) + t414 * t438 - t392) * t405 + t429;
t382 = -t405 * g(3) + (t392 + 0.2e1 * t431) * t406;
t430 = qJDD(1) * t406;
t369 = -pkin(2) * t436 + pkin(5) * t430 + t382;
t355 = t408 * t368 + t411 * t369;
t374 = t390 * mrSges(4,1) + t391 * mrSges(4,2);
t386 = qJD(3) * mrSges(4,1) - t391 * mrSges(4,3);
t377 = t390 * pkin(3) - t391 * pkin(6);
t413 = qJD(3) ^ 2;
t352 = -t413 * pkin(3) + qJDD(3) * pkin(6) - t390 * t377 + t355;
t403 = t405 ^ 2;
t395 = t409 * g(1) - t412 * g(2);
t424 = qJDD(2) - t395;
t378 = (-pkin(1) - t438) * qJDD(1) + (-qJ(2) + (-t403 - t404) * pkin(5)) * t414 + t424;
t433 = t390 * qJD(3);
t380 = t420 * qJDD(1) - t433;
t353 = (-t380 + t433) * pkin(6) + (-t379 + t432) * pkin(3) + t378;
t407 = sin(qJ(4));
t410 = cos(qJ(4));
t349 = -t407 * t352 + t410 * t353;
t383 = t410 * qJD(3) - t407 * t391;
t362 = t383 * qJD(4) + t407 * qJDD(3) + t410 * t380;
t384 = t407 * qJD(3) + t410 * t391;
t363 = -t383 * mrSges(5,1) + t384 * mrSges(5,2);
t388 = qJD(4) + t390;
t365 = -t388 * mrSges(5,2) + t383 * mrSges(5,3);
t376 = qJDD(4) - t379;
t347 = m(5) * t349 + t376 * mrSges(5,1) - t362 * mrSges(5,3) - t384 * t363 + t388 * t365;
t350 = t410 * t352 + t407 * t353;
t361 = -t384 * qJD(4) + t410 * qJDD(3) - t407 * t380;
t366 = t388 * mrSges(5,1) - t384 * mrSges(5,3);
t348 = m(5) * t350 - t376 * mrSges(5,2) + t361 * mrSges(5,3) + t383 * t363 - t388 * t366;
t425 = -t407 * t347 + t410 * t348;
t338 = m(4) * t355 - qJDD(3) * mrSges(4,2) + t379 * mrSges(4,3) - qJD(3) * t386 - t390 * t374 + t425;
t354 = t411 * t368 - t408 * t369;
t385 = -qJD(3) * mrSges(4,2) - t390 * mrSges(4,3);
t351 = -qJDD(3) * pkin(3) - t413 * pkin(6) + t391 * t377 - t354;
t417 = -m(5) * t351 + t361 * mrSges(5,1) - t362 * mrSges(5,2) + t383 * t365 - t384 * t366;
t343 = m(4) * t354 + qJDD(3) * mrSges(4,1) - t380 * mrSges(4,3) + qJD(3) * t385 - t391 * t374 + t417;
t333 = t408 * t338 + t411 * t343;
t381 = -t405 * t392 + t429;
t418 = mrSges(3,3) * qJDD(1) + t414 * (-mrSges(3,1) * t406 + t437);
t331 = m(3) * t381 - t418 * t405 + t333;
t426 = t411 * t338 - t408 * t343;
t332 = m(3) * t382 + t418 * t406 + t426;
t427 = -t405 * t331 + t406 * t332;
t324 = m(2) * t396 - t414 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t427;
t389 = -qJDD(1) * pkin(1) - t414 * qJ(2) + t424;
t339 = t410 * t347 + t407 * t348;
t416 = m(4) * t378 - t379 * mrSges(4,1) + t380 * mrSges(4,2) + t390 * t385 + t391 * t386 + t339;
t415 = -m(3) * t389 + mrSges(3,1) * t430 - t416 + (t403 * t414 + t436) * mrSges(3,3);
t335 = (mrSges(2,1) - t437) * qJDD(1) + t415 - t414 * mrSges(2,2) + m(2) * t395;
t435 = t409 * t324 + t412 * t335;
t325 = t406 * t331 + t405 * t332;
t421 = Ifges(3,5) * t405 + Ifges(3,6) * t406;
t434 = t414 * t421;
t428 = t412 * t324 - t409 * t335;
t423 = Ifges(3,1) * t405 + Ifges(3,4) * t406;
t422 = Ifges(3,4) * t405 + Ifges(3,2) * t406;
t372 = Ifges(4,1) * t391 - Ifges(4,4) * t390 + Ifges(4,5) * qJD(3);
t371 = Ifges(4,4) * t391 - Ifges(4,2) * t390 + Ifges(4,6) * qJD(3);
t370 = Ifges(4,5) * t391 - Ifges(4,6) * t390 + Ifges(4,3) * qJD(3);
t358 = Ifges(5,1) * t384 + Ifges(5,4) * t383 + Ifges(5,5) * t388;
t357 = Ifges(5,4) * t384 + Ifges(5,2) * t383 + Ifges(5,6) * t388;
t356 = Ifges(5,5) * t384 + Ifges(5,6) * t383 + Ifges(5,3) * t388;
t341 = mrSges(5,2) * t351 - mrSges(5,3) * t349 + Ifges(5,1) * t362 + Ifges(5,4) * t361 + Ifges(5,5) * t376 + t383 * t356 - t388 * t357;
t340 = -mrSges(5,1) * t351 + mrSges(5,3) * t350 + Ifges(5,4) * t362 + Ifges(5,2) * t361 + Ifges(5,6) * t376 - t384 * t356 + t388 * t358;
t327 = -mrSges(4,1) * t378 - mrSges(5,1) * t349 + mrSges(5,2) * t350 + mrSges(4,3) * t355 + Ifges(4,4) * t380 - Ifges(5,5) * t362 + Ifges(4,2) * t379 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t361 - Ifges(5,3) * t376 - pkin(3) * t339 + qJD(3) * t372 - t384 * t357 + t383 * t358 - t391 * t370;
t326 = mrSges(4,2) * t378 - mrSges(4,3) * t354 + Ifges(4,1) * t380 + Ifges(4,4) * t379 + Ifges(4,5) * qJDD(3) - pkin(6) * t339 - qJD(3) * t371 - t407 * t340 + t410 * t341 - t390 * t370;
t321 = mrSges(3,2) * t389 - mrSges(3,3) * t381 - pkin(5) * t333 + t423 * qJDD(1) + t411 * t326 - t408 * t327 + t406 * t434;
t320 = -mrSges(3,1) * t389 + mrSges(3,3) * t382 - pkin(2) * t416 + pkin(5) * t426 + t422 * qJDD(1) + t408 * t326 + t411 * t327 - t405 * t434;
t319 = mrSges(2,1) * g(3) - pkin(1) * t325 + mrSges(2,3) * t396 - pkin(2) * t333 - mrSges(3,1) * t381 + mrSges(3,2) * t382 - mrSges(4,1) * t354 + mrSges(4,2) * t355 - t407 * t341 - t410 * t340 - pkin(3) * t417 - pkin(6) * t425 - Ifges(4,5) * t380 - Ifges(4,6) * t379 - Ifges(4,3) * qJDD(3) - t391 * t371 - t390 * t372 + (Ifges(2,6) - t421) * qJDD(1) + (-t405 * t422 + t406 * t423 + Ifges(2,5)) * t414;
t318 = -mrSges(2,2) * g(3) - mrSges(2,3) * t395 + Ifges(2,5) * qJDD(1) - t414 * Ifges(2,6) - qJ(2) * t325 - t405 * t320 + t406 * t321;
t1 = [-m(1) * g(1) + t428; -m(1) * g(2) + t435; (-m(1) - m(2)) * g(3) + t325; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t435 + t412 * t318 - t409 * t319; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t428 + t409 * t318 + t412 * t319; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t395 - mrSges(2,2) * t396 + t405 * t321 + t406 * t320 + pkin(1) * (-qJDD(1) * t437 + t415) + qJ(2) * t427;];
tauB = t1;
