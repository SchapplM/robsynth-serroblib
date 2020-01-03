% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:40
% EndTime: 2019-12-31 21:07:43
% DurationCPUTime: 1.20s
% Computational Cost: add. (4960->233), mult. (9626->269), div. (0->0), fcn. (5637->6), ass. (0->93)
t436 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t421 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t420 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t435 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t419 = -Ifges(4,6) + Ifges(5,5) + Ifges(6,4);
t434 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t400 = qJD(1) ^ 2;
t396 = sin(qJ(1));
t398 = cos(qJ(1));
t411 = g(1) * t396 - t398 * g(2);
t372 = -qJDD(1) * pkin(1) - pkin(6) * t400 - t411;
t395 = sin(qJ(2));
t397 = cos(qJ(2));
t422 = qJD(1) * qJD(2);
t412 = t397 * t422;
t383 = qJDD(1) * t395 + t412;
t413 = t395 * t422;
t384 = qJDD(1) * t397 - t413;
t321 = (-t383 - t412) * pkin(7) + (-t384 + t413) * pkin(2) + t372;
t408 = -g(1) * t398 - t396 * g(2);
t373 = -pkin(1) * t400 + qJDD(1) * pkin(6) + t408;
t363 = -g(3) * t395 + t397 * t373;
t382 = (-t397 * pkin(2) - t395 * pkin(7)) * qJD(1);
t399 = qJD(2) ^ 2;
t423 = qJD(1) * t397;
t325 = -pkin(2) * t399 + qJDD(2) * pkin(7) + t382 * t423 + t363;
t394 = sin(qJ(3));
t429 = cos(qJ(3));
t318 = t429 * t321 - t394 * t325;
t424 = qJD(1) * t395;
t379 = -t429 * qJD(2) + t394 * t424;
t380 = t394 * qJD(2) + t429 * t424;
t352 = pkin(3) * t379 - qJ(4) * t380;
t378 = qJDD(3) - t384;
t388 = -qJD(3) + t423;
t387 = t388 ^ 2;
t317 = -t378 * pkin(3) - t387 * qJ(4) + t380 * t352 + qJDD(4) - t318;
t347 = -t379 * qJD(3) + t394 * qJDD(2) + t429 * t383;
t354 = -mrSges(5,2) * t379 - mrSges(5,3) * t380;
t433 = -m(5) * t317 - t347 * mrSges(5,1) - t380 * t354;
t351 = -mrSges(6,2) * t380 + mrSges(6,3) * t379;
t426 = t379 * t388;
t430 = 2 * qJD(5);
t310 = t388 * t430 + (t379 * t380 - t378) * qJ(5) + (t347 - t426) * pkin(4) + t317;
t360 = -mrSges(6,1) * t379 - mrSges(6,2) * t388;
t409 = -m(6) * t310 + t378 * mrSges(6,3) - t388 * t360;
t308 = mrSges(6,1) * t347 + t351 * t380 - t409;
t359 = mrSges(5,1) * t379 + mrSges(5,3) * t388;
t306 = mrSges(5,2) * t378 - t359 * t388 + t308 - t433;
t346 = qJD(3) * t380 - t429 * qJDD(2) + t383 * t394;
t357 = pkin(4) * t380 + qJ(5) * t388;
t377 = t379 ^ 2;
t319 = t394 * t321 + t429 * t325;
t404 = -pkin(3) * t387 + qJ(4) * t378 - t352 * t379 + t319;
t431 = -2 * qJD(4);
t314 = -pkin(4) * t346 - qJ(5) * t377 + qJDD(5) + (t431 - t357) * t388 + t404;
t315 = 0.2e1 * qJD(4) * t388 - t404;
t361 = mrSges(5,1) * t380 - mrSges(5,2) * t388;
t358 = mrSges(6,1) * t380 + mrSges(6,3) * t388;
t417 = m(6) * t314 + t378 * mrSges(6,2) - t388 * t358;
t406 = -m(5) * t315 + t378 * mrSges(5,3) - t388 * t361 + t417;
t414 = -t421 * t379 + t436 * t380 - t420 * t388;
t415 = t435 * t379 + t421 * t380 + t419 * t388;
t425 = -t351 - t354;
t432 = t419 * t346 + t420 * t347 + t434 * t378 + t414 * t379 + t415 * t380 + mrSges(4,1) * t318 - mrSges(4,2) * t319 + mrSges(5,2) * t317 + mrSges(6,2) * t314 - mrSges(5,3) * t315 - mrSges(6,3) * t310 - pkin(3) * t306 + qJ(4) * (t425 * t379 + (-mrSges(5,1) - mrSges(6,1)) * t346 + t406) - qJ(5) * t308;
t427 = -mrSges(6,1) - mrSges(4,3);
t362 = -t397 * g(3) - t395 * t373;
t416 = -t419 * t379 - t420 * t380 + t434 * t388;
t353 = mrSges(4,1) * t379 + mrSges(4,2) * t380;
t355 = mrSges(4,2) * t388 - mrSges(4,3) * t379;
t302 = m(4) * t318 + (-t355 + t359) * t388 + (-t351 - t353) * t380 + (mrSges(4,1) - mrSges(5,2)) * t378 + t427 * t347 + t409 + t433;
t356 = -mrSges(4,1) * t388 - mrSges(4,3) * t380;
t304 = m(4) * t319 - mrSges(4,2) * t378 + t356 * t388 + (-t353 + t425) * t379 + (-mrSges(5,1) + t427) * t346 + t406;
t410 = -t302 * t394 + t429 * t304;
t324 = -qJDD(2) * pkin(2) - pkin(7) * t399 + t382 * t424 - t362;
t403 = (-t347 - t426) * qJ(4) + t324 + (-pkin(3) * t388 + t431) * t380;
t312 = -pkin(4) * t377 + t379 * t430 - t357 * t380 + (pkin(3) + qJ(5)) * t346 + t403;
t407 = m(6) * t312 - t347 * mrSges(6,2) + t346 * mrSges(6,3) - t380 * t358 + t379 * t360;
t301 = t429 * t302 + t394 * t304;
t316 = pkin(3) * t346 + t403;
t405 = -m(5) * t316 + t346 * mrSges(5,2) + t379 * t359 - t407;
t402 = -m(4) * t324 - t346 * mrSges(4,1) - t379 * t355 + (-t356 + t361) * t380 + (-mrSges(4,2) + mrSges(5,3)) * t347 + t405;
t386 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t423;
t385 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t424;
t381 = (-t397 * mrSges(3,1) + t395 * mrSges(3,2)) * qJD(1);
t371 = Ifges(3,5) * qJD(2) + (t395 * Ifges(3,1) + t397 * Ifges(3,4)) * qJD(1);
t370 = Ifges(3,6) * qJD(2) + (t395 * Ifges(3,4) + t397 * Ifges(3,2)) * qJD(1);
t369 = Ifges(3,3) * qJD(2) + (t395 * Ifges(3,5) + t397 * Ifges(3,6)) * qJD(1);
t309 = -mrSges(6,1) * t346 - t351 * t379 + t417;
t305 = -mrSges(5,3) * t347 - t361 * t380 - t405;
t300 = mrSges(5,1) * t317 + mrSges(6,1) * t310 + mrSges(4,2) * t324 - mrSges(6,2) * t312 - mrSges(4,3) * t318 - mrSges(5,3) * t316 + pkin(4) * t308 - qJ(4) * t305 - t421 * t346 + t436 * t347 + t420 * t378 + t416 * t379 + t415 * t388;
t299 = -mrSges(4,1) * t324 - mrSges(5,1) * t315 + mrSges(6,1) * t314 + mrSges(5,2) * t316 + mrSges(4,3) * t319 - mrSges(6,3) * t312 - pkin(3) * t305 + pkin(4) * t309 - qJ(5) * t407 + t435 * t346 + t421 * t347 - t419 * t378 + t416 * t380 - t414 * t388;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t411 - mrSges(2,2) * t408 + t395 * (mrSges(3,2) * t372 - mrSges(3,3) * t362 + Ifges(3,1) * t383 + Ifges(3,4) * t384 + Ifges(3,5) * qJDD(2) - pkin(7) * t301 - qJD(2) * t370 - t394 * t299 + t429 * t300 + t369 * t423) + t397 * (-mrSges(3,1) * t372 + mrSges(3,3) * t363 + Ifges(3,4) * t383 + Ifges(3,2) * t384 + Ifges(3,6) * qJDD(2) - pkin(2) * t301 + qJD(2) * t371 - t369 * t424 - t432) + pkin(1) * (-m(3) * t372 + t384 * mrSges(3,1) - t383 * mrSges(3,2) + (-t385 * t395 + t386 * t397) * qJD(1) - t301) + pkin(6) * (t397 * (m(3) * t363 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t384 - qJD(2) * t385 + t381 * t423 + t410) - t395 * (m(3) * t362 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t383 + qJD(2) * t386 - t381 * t424 + t402)); Ifges(3,5) * t383 + Ifges(3,6) * t384 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t362 - mrSges(3,2) * t363 + t394 * t300 + t429 * t299 + pkin(2) * t402 + pkin(7) * t410 + (t395 * t370 - t397 * t371) * qJD(1); t432; t306; t309;];
tauJ = t1;
