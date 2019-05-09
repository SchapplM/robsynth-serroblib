% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-05-06 09:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:14:27
% EndTime: 2019-05-06 09:14:37
% DurationCPUTime: 4.34s
% Computational Cost: add. (39310->389), mult. (81176->442), div. (0->0), fcn. (40118->6), ass. (0->139)
t378 = sin(qJ(5));
t381 = cos(qJ(5));
t382 = cos(qJ(2));
t421 = qJD(1) * t382;
t331 = -qJD(2) * t381 + t378 * t421;
t379 = sin(qJ(2));
t420 = qJD(1) * qJD(2);
t410 = t379 * t420;
t339 = qJDD(1) * t382 - t410;
t268 = qJD(5) * t331 - qJDD(2) * t378 - t339 * t381;
t332 = qJD(2) * t378 + t381 * t421;
t270 = -mrSges(7,1) * t331 - mrSges(7,2) * t332;
t412 = t382 * t420;
t338 = qJDD(1) * t379 + t412;
t422 = qJD(1) * t379;
t345 = -qJD(2) * pkin(3) - qJ(4) * t422;
t380 = sin(qJ(1));
t383 = cos(qJ(1));
t352 = t380 * g(1) - t383 * g(2);
t385 = qJD(1) ^ 2;
t300 = -qJDD(1) * pkin(1) - t385 * pkin(7) - t352;
t404 = -t339 * pkin(2) + t300 + (-t338 - t412) * qJ(3);
t430 = t382 ^ 2 * t385;
t395 = -qJ(4) * t430 + qJDD(4) - t404 + ((2 * qJD(3)) + t345) * t422;
t439 = pkin(3) + pkin(8);
t440 = -pkin(2) - pkin(8);
t224 = (pkin(4) * t382 + t379 * t440) * t420 + t395 + t439 * t339 + pkin(4) * t338;
t337 = (pkin(4) * t379 + pkin(8) * t382) * qJD(1);
t384 = qJD(2) ^ 2;
t353 = -g(1) * t383 - g(2) * t380;
t301 = -pkin(1) * t385 + qJDD(1) * pkin(7) + t353;
t277 = -t382 * g(3) - t379 * t301;
t333 = (-pkin(2) * t382 - qJ(3) * t379) * qJD(1);
t407 = t333 * t422 + qJDD(3) - t277;
t429 = t382 * t385;
t419 = qJD(1) * qJD(4);
t445 = -0.2e1 * t379 * t419 + (-t338 + t412) * qJ(4);
t228 = (-pkin(4) - qJ(3)) * t384 + (-pkin(3) * t429 - qJD(1) * t337) * t379 + (-pkin(2) - t439) * qJDD(2) + t407 + t445;
t217 = t381 * t224 - t228 * t378;
t329 = qJDD(5) + t338;
t356 = qJD(5) + t422;
t441 = 2 * qJD(6);
t212 = t332 * t441 + (t331 * t356 - t268) * qJ(6) + (-t331 * t332 + t329) * pkin(5) + t217;
t272 = -mrSges(7,2) * t356 + mrSges(7,3) * t331;
t417 = m(7) * t212 + t329 * mrSges(7,1) + t356 * t272;
t209 = -mrSges(7,3) * t268 + t270 * t332 + t417;
t218 = t378 * t224 + t381 * t228;
t251 = -Ifges(6,4) * t332 + Ifges(6,2) * t331 + Ifges(6,6) * t356;
t252 = -Ifges(7,1) * t332 + Ifges(7,4) * t331 + Ifges(7,5) * t356;
t253 = -Ifges(6,1) * t332 + Ifges(6,4) * t331 + Ifges(6,5) * t356;
t267 = qJD(5) * t332 - qJDD(2) * t381 + t339 * t378;
t274 = pkin(5) * t356 + qJ(6) * t332;
t328 = t331 ^ 2;
t215 = -pkin(5) * t328 + qJ(6) * t267 - t274 * t356 + t331 * t441 + t218;
t250 = -Ifges(7,4) * t332 + Ifges(7,2) * t331 + Ifges(7,6) * t356;
t402 = -mrSges(7,1) * t212 + mrSges(7,2) * t215 - Ifges(7,5) * t268 - Ifges(7,6) * t267 - Ifges(7,3) * t329 + t332 * t250;
t451 = mrSges(6,1) * t217 - mrSges(6,2) * t218 + Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * t329 + pkin(5) * t209 - t332 * t251 - t402 - (t252 + t253) * t331;
t271 = -mrSges(6,1) * t331 - mrSges(6,2) * t332;
t273 = -mrSges(6,2) * t356 + mrSges(6,3) * t331;
t200 = m(6) * t217 + mrSges(6,1) * t329 + t273 * t356 + (t270 + t271) * t332 + (-mrSges(6,3) - mrSges(7,3)) * t268 + t417;
t275 = mrSges(7,1) * t356 + mrSges(7,3) * t332;
t276 = mrSges(6,1) * t356 + mrSges(6,3) * t332;
t416 = m(7) * t215 + t267 * mrSges(7,3) + t331 * t270;
t203 = m(6) * t218 + mrSges(6,3) * t267 + t271 * t331 + (-t275 - t276) * t356 + (-mrSges(6,2) - mrSges(7,2)) * t329 + t416;
t194 = t381 * t200 + t378 * t203;
t232 = -pkin(2) * t410 + pkin(3) * t339 + t395;
t245 = -qJDD(2) * pkin(2) - qJ(3) * t384 + t407;
t236 = (-t379 * t429 - qJDD(2)) * pkin(3) + t245 + t445;
t297 = -Ifges(5,5) * qJD(2) + (-Ifges(5,1) * t382 - Ifges(5,4) * t379) * qJD(1);
t450 = -mrSges(5,1) * t232 + mrSges(5,3) * t236 - Ifges(5,6) * qJDD(2) - pkin(4) * t194 - qJD(2) * t297 - t451;
t278 = -g(3) * t379 + t382 * t301;
t292 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t379 - Ifges(4,3) * t382) * qJD(1);
t296 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t379 + Ifges(3,2) * t382) * qJD(1);
t334 = (-mrSges(4,1) * t382 - mrSges(4,3) * t379) * qJD(1);
t418 = qJD(3) * qJD(2);
t361 = 0.2e1 * t418;
t437 = pkin(2) * t384;
t447 = qJDD(2) * qJ(3) + t333 * t421 + t278;
t243 = t361 - t437 + t447;
t348 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t422;
t448 = pkin(3) * t430 + t339 * qJ(4) - qJD(2) * t345 + 0.2e1 * t382 * t419 - t447;
t235 = -0.2e1 * t418 + t437 + t448;
t346 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t422;
t227 = qJDD(2) * pkin(4) - t337 * t421 + t384 * t440 + t361 - t448;
t221 = -pkin(5) * t267 - qJ(6) * t328 - t274 * t332 + qJDD(6) + t227;
t415 = -m(7) * t221 - t268 * mrSges(7,2) + t332 * t275;
t442 = -m(6) * t227 - t268 * mrSges(6,2) + (mrSges(6,1) + mrSges(7,1)) * t267 + t332 * t276 - (-t273 - t272) * t331 + t415;
t393 = -m(5) * t235 + qJDD(2) * mrSges(5,1) - t339 * mrSges(5,3) + qJD(2) * t346 - t442;
t389 = m(4) * t243 + qJDD(2) * mrSges(4,3) + qJD(2) * t348 + t334 * t421 + t393;
t349 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t421;
t336 = (mrSges(5,1) * t379 - mrSges(5,2) * t382) * qJD(1);
t428 = -t378 * t200 + t381 * t203;
t406 = -m(5) * t236 + t336 * t422 - t428;
t398 = -qJDD(2) * mrSges(5,2) - qJD(2) * t349 + t406;
t188 = -mrSges(5,3) * t338 - t398;
t248 = -Ifges(7,5) * t332 + Ifges(7,6) * t331 + Ifges(7,3) * t356;
t249 = -Ifges(6,5) * t332 + Ifges(6,6) * t331 + Ifges(6,3) * t356;
t403 = -mrSges(7,1) * t221 + mrSges(7,3) * t215 + Ifges(7,4) * t268 + Ifges(7,2) * t267 + Ifges(7,6) * t329 + t356 * t252;
t181 = Ifges(6,4) * t268 + Ifges(6,2) * t267 + Ifges(6,6) * t329 + t356 * t253 - mrSges(6,1) * t227 + mrSges(6,3) * t218 - pkin(5) * (-t267 * mrSges(7,1) - t331 * t272 - t415) + qJ(6) * (-mrSges(7,2) * t329 - t275 * t356 + t416) + (t249 + t248) * t332 + t403;
t401 = mrSges(7,2) * t221 - mrSges(7,3) * t212 + Ifges(7,1) * t268 + Ifges(7,4) * t267 + Ifges(7,5) * t329 + t331 * t248;
t191 = mrSges(6,2) * t227 - mrSges(6,3) * t217 + Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * t329 - qJ(6) * t209 + t249 * t331 + (-t250 - t251) * t356 + t401;
t294 = -Ifges(5,6) * qJD(2) + (-Ifges(5,4) * t382 - Ifges(5,2) * t379) * qJD(1);
t397 = mrSges(5,1) * t235 - mrSges(5,2) * t236 - Ifges(5,5) * t339 - Ifges(5,6) * t338 - Ifges(5,3) * qJDD(2) + pkin(4) * t442 + pkin(8) * t428 + t381 * t181 + t378 * t191 - t294 * t421 + t297 * t422;
t392 = -mrSges(4,1) * t245 + mrSges(4,3) * t243 + Ifges(4,4) * t338 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t339 - pkin(3) * t188 - t397;
t413 = t336 * t421;
t298 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t379 - Ifges(4,5) * t382) * qJD(1);
t423 = t298 + Ifges(3,5) * qJD(2) + (Ifges(3,1) * t379 + Ifges(3,4) * t382) * qJD(1);
t351 = mrSges(4,2) * t421 + qJD(2) * mrSges(4,3);
t446 = -m(4) * t245 + qJDD(2) * mrSges(4,1) + qJD(2) * t351;
t449 = -((t292 - t296) * t379 + t382 * t423) * qJD(1) + mrSges(3,1) * t277 - mrSges(3,2) * t278 + Ifges(3,5) * t338 + Ifges(3,6) * t339 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t334 * t422 + (-mrSges(4,2) + mrSges(5,3)) * t338 + t398 + t446) + qJ(3) * (t339 * mrSges(4,2) + t389 - t413) + t392;
t291 = -Ifges(5,3) * qJD(2) + (-Ifges(5,5) * t382 - Ifges(5,6) * t379) * qJD(1);
t444 = Ifges(5,4) * t339 + Ifges(5,2) * t338 - t291 * t421;
t435 = mrSges(3,3) + mrSges(4,2);
t434 = Ifges(4,6) - Ifges(5,5);
t295 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t379 - Ifges(4,6) * t382) * qJD(1);
t425 = -t291 + t295;
t335 = (-mrSges(3,1) * t382 + mrSges(3,2) * t379) * qJD(1);
t350 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t421;
t184 = m(3) * t277 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t349 + t350) * qJD(2) + (-t334 - t335) * t422 + (mrSges(5,3) - t435) * t338 + t406 + t446;
t347 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t422;
t197 = -qJDD(2) * mrSges(3,2) + t435 * t339 - qJD(2) * t347 + m(3) * t278 + t389 + (t335 - t336) * t421;
t408 = -t184 * t379 + t382 * t197;
t187 = -m(5) * t232 - t338 * mrSges(5,1) + t339 * mrSges(5,2) - t346 * t422 + t349 * t421 - t194;
t237 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t422 + t404;
t185 = m(4) * t237 - mrSges(4,1) * t339 - t338 * mrSges(4,3) - t348 * t422 - t351 * t421 + t187;
t293 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t379 + Ifges(3,6) * t382) * qJD(1);
t396 = -mrSges(5,2) * t232 + mrSges(5,3) * t235 + Ifges(5,1) * t339 + Ifges(5,4) * t338 + pkin(8) * t194 - qJD(2) * t294 + t378 * t181 - t381 * t191;
t391 = mrSges(4,1) * t237 - mrSges(4,2) * t243 + pkin(3) * t187 + qJ(4) * (t393 - t413) - t396;
t174 = (Ifges(4,3) + Ifges(3,2)) * t339 + (Ifges(3,4) - Ifges(4,5)) * t338 + (Ifges(3,6) - t434) * qJDD(2) + t423 * qJD(2) - mrSges(3,1) * t300 + mrSges(3,3) * t278 - t391 - pkin(2) * t185 + (-t293 - t425) * t422;
t386 = mrSges(4,2) * t245 - mrSges(4,3) * t237 + Ifges(4,1) * t338 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t339 - qJ(4) * t188 + qJD(2) * t292 + t295 * t421 - t450;
t176 = (Ifges(5,4) + Ifges(3,4)) * t339 + (Ifges(3,1) + Ifges(5,2)) * t338 + Ifges(3,5) * qJDD(2) - qJD(2) * t296 + mrSges(3,2) * t300 - mrSges(3,3) * t277 + t386 - qJ(3) * t185 + (-t291 + t293) * t421;
t388 = -m(3) * t300 + t339 * mrSges(3,1) - mrSges(3,2) * t338 - t347 * t422 + t350 * t421 - t185;
t400 = mrSges(2,1) * t352 - mrSges(2,2) * t353 + Ifges(2,3) * qJDD(1) + pkin(1) * t388 + pkin(7) * t408 + t382 * t174 + t379 * t176;
t182 = m(2) * t352 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t385 + t388;
t179 = t184 * t382 + t197 * t379;
t177 = m(2) * t353 - mrSges(2,1) * t385 - qJDD(1) * mrSges(2,2) + t408;
t172 = mrSges(2,1) * g(3) + mrSges(2,3) * t353 + t385 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t179 - t449;
t171 = -mrSges(2,2) * g(3) - mrSges(2,3) * t352 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t385 - pkin(7) * t179 - t174 * t379 + t176 * t382;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t383 * t171 - t380 * t172 - pkin(6) * (t177 * t380 + t182 * t383), t171, t176, t386 + t444, -Ifges(5,5) * qJDD(2) - t291 * t422 - t396, t191, -t250 * t356 + t401; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t380 * t171 + t383 * t172 + pkin(6) * (t177 * t383 - t182 * t380), t172, t174, t392 + (-t379 * t292 - t382 * t298) * qJD(1), -t444 + t450, t181, t332 * t248 + t403; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t400, t400, t449, Ifges(4,5) * t338 - Ifges(4,3) * t339 - qJD(2) * t298 + qJDD(2) * t434 + t422 * t425 + t391, t397, t451, -t331 * t252 - t402;];
m_new  = t1;
