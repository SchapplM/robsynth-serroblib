% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:33:35
% EndTime: 2019-05-06 08:33:48
% DurationCPUTime: 6.97s
% Computational Cost: add. (84669->392), mult. (184847->462), div. (0->0), fcn. (98735->8), ass. (0->148)
t380 = sin(qJ(1));
t383 = cos(qJ(1));
t351 = -g(1) * t383 - g(2) * t380;
t385 = qJD(1) ^ 2;
t302 = -pkin(1) * t385 + qJDD(1) * pkin(7) + t351;
t379 = sin(qJ(2));
t382 = cos(qJ(2));
t268 = -t382 * g(3) - t379 * t302;
t269 = -g(3) * t379 + t382 * t302;
t293 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t379 - Ifges(4,3) * t382) * qJD(1);
t297 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t379 + Ifges(3,2) * t382) * qJD(1);
t332 = (-mrSges(4,1) * t382 - mrSges(4,3) * t379) * qJD(1);
t416 = qJD(1) * qJD(2);
t411 = t382 * t416;
t336 = qJDD(1) * t379 + t411;
t409 = t379 * t416;
t337 = qJDD(1) * t382 - t409;
t414 = qJD(3) * qJD(2);
t359 = 0.2e1 * t414;
t384 = qJD(2) ^ 2;
t432 = pkin(2) * t384;
t331 = (-pkin(2) * t382 - qJ(3) * t379) * qJD(1);
t417 = qJD(1) * t382;
t439 = qJDD(2) * qJ(3) + t331 * t417 + t269;
t256 = t359 - t432 + t439;
t418 = qJD(1) * t379;
t346 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t418;
t343 = -qJD(2) * pkin(3) - qJ(4) * t418;
t415 = qJD(1) * qJD(4);
t424 = t382 ^ 2 * t385;
t440 = pkin(3) * t424 + t337 * qJ(4) - qJD(2) * t343 + 0.2e1 * t382 * t415 - t439;
t243 = -0.2e1 * t414 + t432 + t440;
t344 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t418;
t334 = (pkin(4) * t379 + qJ(5) * t382) * qJD(1);
t429 = -pkin(2) - qJ(5);
t232 = qJDD(2) * pkin(4) - t334 * t417 + t384 * t429 + qJDD(5) + t359 - t440;
t375 = sin(pkin(9));
t376 = cos(pkin(9));
t313 = -qJD(2) * t376 + t375 * t417;
t273 = -mrSges(6,2) * t418 + mrSges(6,3) * t313;
t314 = qJD(2) * t375 + t376 * t417;
t274 = mrSges(6,1) * t418 + mrSges(6,3) * t314;
t275 = -qJDD(2) * t376 + t337 * t375;
t276 = -qJDD(2) * t375 - t337 * t376;
t277 = pkin(5) * t418 + pkin(8) * t314;
t310 = t313 ^ 2;
t225 = -pkin(5) * t275 - pkin(8) * t310 - t277 * t314 + t232;
t378 = sin(qJ(6));
t381 = cos(qJ(6));
t266 = t313 * t378 - t314 * t381;
t245 = -qJD(6) * t266 + t275 * t381 - t276 * t378;
t265 = t313 * t381 + t314 * t378;
t246 = qJD(6) * t265 + t275 * t378 + t276 * t381;
t354 = qJD(6) + t418;
t259 = -mrSges(7,2) * t354 + mrSges(7,3) * t265;
t260 = mrSges(7,1) * t354 - mrSges(7,3) * t266;
t403 = m(7) * t225 - t245 * mrSges(7,1) + t246 * mrSges(7,2) - t265 * t259 + t266 * t260;
t395 = -m(6) * t232 + t275 * mrSges(6,1) - t276 * mrSges(6,2) + t313 * t273 + t314 * t274 - t403;
t392 = -m(5) * t243 + qJDD(2) * mrSges(5,1) - t337 * mrSges(5,3) + qJD(2) * t344 - t395;
t389 = m(4) * t256 + qJDD(2) * mrSges(4,3) + qJD(2) * t346 + t332 * t417 + t392;
t347 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t417;
t406 = t331 * t418 + qJDD(3) - t268;
t258 = -qJDD(2) * pkin(2) - qJ(3) * t384 + t406;
t423 = t382 * t385;
t437 = -0.2e1 * t379 * t415 + (-t336 + t411) * qJ(4);
t244 = (-t379 * t423 - qJDD(2)) * pkin(3) + t258 + t437;
t335 = (mrSges(5,1) * t379 - mrSges(5,2) * t382) * qJD(1);
t350 = t380 * g(1) - t383 * g(2);
t301 = -qJDD(1) * pkin(1) - t385 * pkin(7) - t350;
t404 = -t337 * pkin(2) + t301 + (-t336 - t411) * qJ(3);
t396 = -qJ(4) * t424 + qJDD(4) - t404 + (0.2e1 * qJD(3) + t343) * t418;
t428 = pkin(3) + qJ(5);
t228 = t396 + t428 * t337 + pkin(4) * t336 + (pkin(4) * t382 + t379 * t429) * t416;
t233 = (-pkin(4) - qJ(3)) * t384 + (-pkin(3) * t423 - qJD(1) * t334) * t379 + (-pkin(2) - t428) * qJDD(2) + t406 + t437;
t434 = 2 * qJD(5);
t222 = t376 * t228 - t233 * t375 + t314 * t434;
t219 = (t313 * t418 - t276) * pkin(8) + (-t313 * t314 + t336) * pkin(5) + t222;
t223 = t375 * t228 + t376 * t233 + t313 * t434;
t220 = -pkin(5) * t310 + pkin(8) * t275 - t277 * t418 + t223;
t217 = t219 * t381 - t220 * t378;
t251 = -mrSges(7,1) * t265 + mrSges(7,2) * t266;
t329 = qJDD(6) + t336;
t211 = m(7) * t217 + mrSges(7,1) * t329 - mrSges(7,3) * t246 - t251 * t266 + t259 * t354;
t218 = t219 * t378 + t220 * t381;
t212 = m(7) * t218 - mrSges(7,2) * t329 + mrSges(7,3) * t245 + t251 * t265 - t260 * t354;
t202 = t381 * t211 + t378 * t212;
t267 = -mrSges(6,1) * t313 - mrSges(6,2) * t314;
t199 = m(6) * t222 + mrSges(6,1) * t336 - mrSges(6,3) * t276 + t267 * t314 + t273 * t418 + t202;
t407 = -t211 * t378 + t381 * t212;
t200 = m(6) * t223 - mrSges(6,2) * t336 + mrSges(6,3) * t275 + t267 * t313 - t274 * t418 + t407;
t422 = -t375 * t199 + t376 * t200;
t405 = -m(5) * t244 + t335 * t418 - t422;
t399 = -qJDD(2) * mrSges(5,2) - qJD(2) * t347 + t405;
t190 = -mrSges(5,3) * t336 - t399;
t247 = Ifges(7,5) * t266 + Ifges(7,6) * t265 + Ifges(7,3) * t354;
t249 = Ifges(7,1) * t266 + Ifges(7,4) * t265 + Ifges(7,5) * t354;
t203 = -mrSges(7,1) * t225 + mrSges(7,3) * t218 + Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t329 - t247 * t266 + t249 * t354;
t248 = Ifges(7,4) * t266 + Ifges(7,2) * t265 + Ifges(7,6) * t354;
t204 = mrSges(7,2) * t225 - mrSges(7,3) * t217 + Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t329 + t247 * t265 - t248 * t354;
t261 = -Ifges(6,5) * t314 + Ifges(6,6) * t313 + Ifges(6,3) * t418;
t263 = -Ifges(6,1) * t314 + Ifges(6,4) * t313 + Ifges(6,5) * t418;
t180 = -mrSges(6,1) * t232 + mrSges(6,3) * t223 + Ifges(6,4) * t276 + Ifges(6,2) * t275 + Ifges(6,6) * t336 - pkin(5) * t403 + pkin(8) * t407 + t381 * t203 + t378 * t204 + t314 * t261 + t263 * t418;
t262 = -Ifges(6,4) * t314 + Ifges(6,2) * t313 + Ifges(6,6) * t418;
t183 = mrSges(6,2) * t232 - mrSges(6,3) * t222 + Ifges(6,1) * t276 + Ifges(6,4) * t275 + Ifges(6,5) * t336 - pkin(8) * t202 - t203 * t378 + t204 * t381 + t261 * t313 - t262 * t418;
t295 = -Ifges(5,6) * qJD(2) + (-Ifges(5,4) * t382 - Ifges(5,2) * t379) * qJD(1);
t298 = -Ifges(5,5) * qJD(2) + (-Ifges(5,1) * t382 - Ifges(5,4) * t379) * qJD(1);
t398 = mrSges(5,1) * t243 - mrSges(5,2) * t244 - Ifges(5,5) * t337 - Ifges(5,6) * t336 - Ifges(5,3) * qJDD(2) + pkin(4) * t395 + qJ(5) * t422 + t376 * t180 + t375 * t183 - t295 * t417 + t298 * t418;
t394 = -mrSges(4,1) * t258 + mrSges(4,3) * t256 + Ifges(4,4) * t336 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t337 - pkin(3) * t190 - t398;
t412 = t335 * t417;
t299 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t379 - Ifges(4,5) * t382) * qJD(1);
t419 = t299 + Ifges(3,5) * qJD(2) + (Ifges(3,1) * t379 + Ifges(3,4) * t382) * qJD(1);
t349 = mrSges(4,2) * t417 + qJD(2) * mrSges(4,3);
t438 = -m(4) * t258 + qJDD(2) * mrSges(4,1) + qJD(2) * t349;
t441 = -qJD(1) * ((t293 - t297) * t379 + t382 * t419) + mrSges(3,1) * t268 - mrSges(3,2) * t269 + Ifges(3,5) * t336 + Ifges(3,6) * t337 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t332 * t418 + (-mrSges(4,2) + mrSges(5,3)) * t336 + t399 + t438) + qJ(3) * (t337 * mrSges(4,2) + t389 - t412) + t394;
t292 = -Ifges(5,3) * qJD(2) + (-Ifges(5,5) * t382 - Ifges(5,6) * t379) * qJD(1);
t436 = Ifges(5,4) * t337 + Ifges(5,2) * t336 - t292 * t417;
t431 = mrSges(3,3) + mrSges(4,2);
t430 = Ifges(4,6) - Ifges(5,5);
t193 = t376 * t199 + t375 * t200;
t296 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t379 - Ifges(4,6) * t382) * qJD(1);
t421 = -t292 + t296;
t333 = (-mrSges(3,1) * t382 + mrSges(3,2) * t379) * qJD(1);
t348 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t417;
t186 = m(3) * t268 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t347 + t348) * qJD(2) + (-t332 - t333) * t418 + (mrSges(5,3) - t431) * t336 + t405 + t438;
t345 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t418;
t206 = -qJD(2) * t345 + m(3) * t269 + (t333 - t335) * t417 - qJDD(2) * mrSges(3,2) + t389 + t431 * t337;
t408 = -t186 * t379 + t382 * t206;
t237 = -pkin(2) * t409 + pkin(3) * t337 + t396;
t189 = -m(5) * t237 - t336 * mrSges(5,1) + t337 * mrSges(5,2) - t344 * t418 + t347 * t417 - t193;
t252 = (pkin(2) * qJD(2) - 0.2e1 * qJD(3)) * t418 + t404;
t187 = m(4) * t252 - mrSges(4,1) * t337 - t336 * mrSges(4,3) - t346 * t418 - t349 * t417 + t189;
t294 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t379 + Ifges(3,6) * t382) * qJD(1);
t397 = -mrSges(5,2) * t237 + mrSges(5,3) * t243 + Ifges(5,1) * t337 + Ifges(5,4) * t336 + qJ(5) * t193 - qJD(2) * t295 + t375 * t180 - t376 * t183;
t391 = mrSges(4,1) * t252 - mrSges(4,2) * t256 + pkin(3) * t189 + qJ(4) * (t392 - t412) - t397;
t173 = (-t294 - t421) * t418 - mrSges(3,1) * t301 + mrSges(3,3) * t269 - pkin(2) * t187 - t391 + t419 * qJD(2) + (Ifges(3,6) - t430) * qJDD(2) + (-Ifges(4,5) + Ifges(3,4)) * t336 + (Ifges(3,2) + Ifges(4,3)) * t337;
t401 = -mrSges(7,1) * t217 + mrSges(7,2) * t218 - Ifges(7,5) * t246 - Ifges(7,6) * t245 - Ifges(7,3) * t329 - t266 * t248 + t265 * t249;
t393 = -mrSges(6,1) * t222 + mrSges(6,2) * t223 - Ifges(6,5) * t276 - Ifges(6,6) * t275 - Ifges(6,3) * t336 - pkin(5) * t202 + t314 * t262 + t313 * t263 + t401;
t388 = -mrSges(5,1) * t237 + mrSges(5,3) * t244 - Ifges(5,6) * qJDD(2) - pkin(4) * t193 - qJD(2) * t298 + t393;
t386 = mrSges(4,2) * t258 - mrSges(4,3) * t252 + Ifges(4,1) * t336 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t337 - qJ(4) * t190 + qJD(2) * t293 + t296 * t417 - t388;
t175 = (-t292 + t294) * t417 - qJD(2) * t297 + mrSges(3,2) * t301 - mrSges(3,3) * t268 - qJ(3) * t187 + Ifges(3,5) * qJDD(2) + t386 + (Ifges(3,4) + Ifges(5,4)) * t337 + (Ifges(5,2) + Ifges(3,1)) * t336;
t390 = -m(3) * t301 + t337 * mrSges(3,1) - mrSges(3,2) * t336 - t345 * t418 + t348 * t417 - t187;
t402 = mrSges(2,1) * t350 - mrSges(2,2) * t351 + Ifges(2,3) * qJDD(1) + pkin(1) * t390 + pkin(7) * t408 + t382 * t173 + t379 * t175;
t184 = m(2) * t350 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t385 + t390;
t178 = t186 * t382 + t206 * t379;
t176 = m(2) * t351 - mrSges(2,1) * t385 - qJDD(1) * mrSges(2,2) + t408;
t171 = mrSges(2,1) * g(3) + mrSges(2,3) * t351 + t385 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t178 - t441;
t170 = -mrSges(2,2) * g(3) - mrSges(2,3) * t350 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t385 - pkin(7) * t178 - t173 * t379 + t175 * t382;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t383 * t170 - t380 * t171 - pkin(6) * (t176 * t380 + t184 * t383), t170, t175, t386 + t436, -Ifges(5,5) * qJDD(2) - t292 * t418 - t397, t183, t204; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t380 * t170 + t383 * t171 + pkin(6) * (t176 * t383 - t184 * t380), t171, t173, t394 + (-t379 * t293 - t382 * t299) * qJD(1), t388 - t436, t180, t203; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t402, t402, t441, Ifges(4,5) * t336 - Ifges(4,3) * t337 - qJD(2) * t299 + qJDD(2) * t430 + t418 * t421 + t391, t398, -t393, -t401;];
m_new  = t1;
