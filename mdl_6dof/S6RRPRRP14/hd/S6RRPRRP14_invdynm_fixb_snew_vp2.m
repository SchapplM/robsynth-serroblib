% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:17:29
% EndTime: 2019-05-06 19:18:01
% DurationCPUTime: 13.52s
% Computational Cost: add. (221466->402), mult. (492940->487), div. (0->0), fcn. (354899->10), ass. (0->153)
t403 = -2 * qJD(3);
t348 = cos(pkin(6));
t341 = qJD(1) * t348 + qJD(2);
t351 = sin(qJ(2));
t347 = sin(pkin(6));
t382 = qJD(1) * t347;
t375 = t351 * t382;
t402 = (pkin(2) * t341 + t403) * t375;
t352 = sin(qJ(1));
t355 = cos(qJ(1));
t336 = t352 * g(1) - g(2) * t355;
t356 = qJD(1) ^ 2;
t398 = pkin(8) * t347;
t318 = qJDD(1) * pkin(1) + t356 * t398 + t336;
t337 = -g(1) * t355 - g(2) * t352;
t379 = qJDD(1) * t347;
t319 = -pkin(1) * t356 + pkin(8) * t379 + t337;
t354 = cos(qJ(2));
t390 = t348 * t351;
t392 = t347 * t351;
t279 = -g(3) * t392 + t318 * t390 + t354 * t319;
t320 = (-pkin(2) * t354 - qJ(3) * t351) * t382;
t339 = t341 ^ 2;
t340 = qJDD(1) * t348 + qJDD(2);
t381 = qJD(1) * t354;
t374 = t347 * t381;
t241 = t339 * pkin(2) - t340 * qJ(3) - t320 * t374 + t341 * t403 - t279;
t323 = pkin(3) * t375 - pkin(9) * t341;
t324 = (qJD(2) * t381 + qJDD(1) * t351) * t347;
t325 = -qJD(2) * t375 + t354 * t379;
t393 = t347 ^ 2 * t356;
t378 = t354 ^ 2 * t393;
t397 = t348 * g(3);
t400 = -pkin(2) - pkin(9);
t227 = -pkin(3) * t378 - t397 - t324 * qJ(3) + t400 * t325 + (-t318 + (-qJ(3) * t341 * t354 - t323 * t351) * qJD(1)) * t347 + t402;
t391 = t347 * t354;
t383 = g(3) * t391 + t351 * t319;
t369 = -t339 * qJ(3) + t320 * t375 + qJDD(3) + t383;
t229 = t324 * pkin(3) + t400 * t340 + (-pkin(3) * t341 * t382 - pkin(9) * t351 * t393 - t318 * t348) * t354 + t369;
t350 = sin(qJ(4));
t353 = cos(qJ(4));
t222 = t353 * t227 + t350 * t229;
t304 = -t341 * t350 - t353 * t374;
t305 = t341 * t353 - t350 * t374;
t281 = -pkin(4) * t304 - pkin(10) * t305;
t313 = qJDD(4) + t324;
t331 = qJD(4) + t375;
t329 = t331 ^ 2;
t217 = -pkin(4) * t329 + pkin(10) * t313 + t281 * t304 + t222;
t226 = t325 * pkin(3) - pkin(9) * t378 + t341 * t323 - t241;
t276 = -qJD(4) * t305 - t325 * t353 - t340 * t350;
t277 = qJD(4) * t304 - t325 * t350 + t340 * t353;
t219 = (-t304 * t331 - t277) * pkin(10) + (t305 * t331 - t276) * pkin(4) + t226;
t349 = sin(qJ(5));
t399 = cos(qJ(5));
t213 = -t349 * t217 + t399 * t219;
t214 = t399 * t217 + t349 * t219;
t283 = t399 * t305 + t349 * t331;
t238 = qJD(5) * t283 + t277 * t349 - t399 * t313;
t282 = t305 * t349 - t399 * t331;
t239 = -t282 * qJD(5) + t399 * t277 + t349 * t313;
t302 = qJD(5) - t304;
t243 = Ifges(7,5) * t283 + Ifges(7,6) * t302 + Ifges(7,3) * t282;
t246 = Ifges(6,4) * t283 - Ifges(6,2) * t282 + Ifges(6,6) * t302;
t248 = Ifges(6,1) * t283 - Ifges(6,4) * t282 + Ifges(6,5) * t302;
t256 = mrSges(7,1) * t282 - mrSges(7,3) * t283;
t274 = qJDD(5) - t276;
t255 = pkin(5) * t282 - qJ(6) * t283;
t301 = t302 ^ 2;
t209 = -pkin(5) * t301 + qJ(6) * t274 + 0.2e1 * qJD(6) * t302 - t255 * t282 + t214;
t211 = -t274 * pkin(5) - t301 * qJ(6) + t283 * t255 + qJDD(6) - t213;
t247 = Ifges(7,1) * t283 + Ifges(7,4) * t302 + Ifges(7,5) * t282;
t368 = mrSges(7,1) * t211 - mrSges(7,3) * t209 - Ifges(7,4) * t239 - Ifges(7,2) * t274 - Ifges(7,6) * t238 - t282 * t247;
t260 = -mrSges(7,2) * t282 + mrSges(7,3) * t302;
t372 = -m(7) * t211 + t274 * mrSges(7,1) + t302 * t260;
t263 = -mrSges(7,1) * t302 + mrSges(7,2) * t283;
t376 = m(7) * t209 + t274 * mrSges(7,3) + t302 * t263;
t401 = -(-t246 + t243) * t283 + mrSges(6,1) * t213 - mrSges(6,2) * t214 + Ifges(6,5) * t239 - Ifges(6,6) * t238 + Ifges(6,3) * t274 + pkin(5) * (-t239 * mrSges(7,2) - t283 * t256 + t372) + qJ(6) * (-t238 * mrSges(7,2) - t282 * t256 + t376) + t282 * t248 - t368;
t396 = mrSges(3,1) - mrSges(4,2);
t395 = -mrSges(6,3) - mrSges(7,2);
t394 = Ifges(3,4) + Ifges(4,6);
t389 = t348 * t354;
t262 = mrSges(6,1) * t302 - mrSges(6,3) * t283;
t386 = -mrSges(6,1) * t282 - mrSges(6,2) * t283 - t256;
t199 = m(6) * t214 - t274 * mrSges(6,2) + t395 * t238 - t302 * t262 + t386 * t282 + t376;
t261 = -mrSges(6,2) * t302 - mrSges(6,3) * t282;
t201 = m(6) * t213 + t274 * mrSges(6,1) + t395 * t239 + t302 * t261 + t386 * t283 + t372;
t194 = t349 * t199 + t399 * t201;
t245 = Ifges(7,4) * t283 + Ifges(7,2) * t302 + Ifges(7,6) * t282;
t388 = -Ifges(6,5) * t283 + Ifges(6,6) * t282 - Ifges(6,3) * t302 - t245;
t292 = Ifges(4,1) * t341 + (-Ifges(4,4) * t351 - Ifges(4,5) * t354) * t382;
t385 = Ifges(3,3) * t341 + (Ifges(3,5) * t351 + Ifges(3,6) * t354) * t382 + t292;
t290 = Ifges(4,5) * t341 + (-Ifges(4,6) * t351 - Ifges(4,3) * t354) * t382;
t384 = -Ifges(3,6) * t341 - (Ifges(3,4) * t351 + Ifges(3,2) * t354) * t382 + t290;
t377 = t318 * t389;
t278 = t377 - t383;
t315 = -mrSges(3,2) * t341 + mrSges(3,3) * t374;
t316 = -mrSges(4,1) * t374 - mrSges(4,3) * t341;
t321 = (mrSges(4,2) * t354 - mrSges(4,3) * t351) * t382;
t322 = (-mrSges(3,1) * t354 + mrSges(3,2) * t351) * t382;
t280 = -mrSges(5,1) * t304 + mrSges(5,2) * t305;
t285 = mrSges(5,1) * t331 - mrSges(5,3) * t305;
t373 = t399 * t199 - t201 * t349;
t187 = m(5) * t222 - mrSges(5,2) * t313 + mrSges(5,3) * t276 + t280 * t304 - t285 * t331 + t373;
t221 = -t350 * t227 + t353 * t229;
t284 = -mrSges(5,2) * t331 + mrSges(5,3) * t304;
t216 = -t313 * pkin(4) - t329 * pkin(10) + t305 * t281 - t221;
t212 = -0.2e1 * qJD(6) * t283 + (t282 * t302 - t239) * qJ(6) + (t283 * t302 + t238) * pkin(5) + t216;
t206 = m(7) * t212 + mrSges(7,1) * t238 - t239 * mrSges(7,3) + t260 * t282 - t283 * t263;
t359 = -m(6) * t216 - t238 * mrSges(6,1) - mrSges(6,2) * t239 - t282 * t261 - t262 * t283 - t206;
t196 = m(5) * t221 + mrSges(5,1) * t313 - mrSges(5,3) * t277 - t280 * t305 + t284 * t331 + t359;
t181 = t350 * t187 + t353 * t196;
t252 = -t340 * pkin(2) + t369 - t377;
t366 = -m(4) * t252 - t324 * mrSges(4,1) - t181;
t179 = m(3) * t278 - t324 * mrSges(3,3) + (t315 - t316) * t341 + t396 * t340 + (-t321 - t322) * t375 + t366;
t314 = mrSges(3,1) * t341 - mrSges(3,3) * t375;
t190 = -m(5) * t226 + t276 * mrSges(5,1) - t277 * mrSges(5,2) + t304 * t284 - t305 * t285 - t194;
t317 = mrSges(4,1) * t375 + mrSges(4,2) * t341;
t360 = -m(4) * t241 + t340 * mrSges(4,3) + t341 * t317 + t321 * t374 - t190;
t185 = t322 * t374 + t360 - t340 * mrSges(3,2) - t341 * t314 + m(3) * t279 + (mrSges(3,3) + mrSges(4,1)) * t325;
t175 = -t179 * t351 + t354 * t185;
t182 = t353 * t187 - t350 * t196;
t293 = -t347 * t318 - t397;
t242 = -t325 * pkin(2) + (-t341 * t374 - t324) * qJ(3) + t293 + t402;
t370 = m(4) * t242 - t324 * mrSges(4,3) + t316 * t374 + t182;
t178 = m(3) * t293 + t324 * mrSges(3,2) - t396 * t325 + (-t315 * t354 + (t314 - t317) * t351) * t382 + t370;
t170 = -t178 * t347 + t179 * t389 + t185 * t390;
t371 = -mrSges(7,1) * t212 + mrSges(7,2) * t209;
t367 = mrSges(7,2) * t211 - mrSges(7,3) * t212 + Ifges(7,1) * t239 + Ifges(7,4) * t274 + Ifges(7,5) * t238 + t302 * t243;
t289 = Ifges(3,5) * t341 + (Ifges(3,1) * t351 + Ifges(3,4) * t354) * t382;
t189 = -mrSges(6,1) * t216 + mrSges(6,3) * t214 - pkin(5) * t206 + (t247 + t248) * t302 + t388 * t283 + (Ifges(6,6) - Ifges(7,6)) * t274 + (Ifges(6,4) - Ifges(7,5)) * t239 + (-Ifges(6,2) - Ifges(7,3)) * t238 + t371;
t192 = mrSges(6,2) * t216 - mrSges(6,3) * t213 + Ifges(6,1) * t239 - Ifges(6,4) * t238 + Ifges(6,5) * t274 - qJ(6) * t206 - t302 * t246 + t388 * t282 + t367;
t268 = Ifges(5,5) * t305 + Ifges(5,6) * t304 + Ifges(5,3) * t331;
t269 = Ifges(5,4) * t305 + Ifges(5,2) * t304 + Ifges(5,6) * t331;
t172 = mrSges(5,2) * t226 - mrSges(5,3) * t221 + Ifges(5,1) * t277 + Ifges(5,4) * t276 + Ifges(5,5) * t313 - pkin(10) * t194 - t349 * t189 + t399 * t192 + t304 * t268 - t331 * t269;
t270 = Ifges(5,1) * t305 + Ifges(5,4) * t304 + Ifges(5,5) * t331;
t176 = -mrSges(5,1) * t226 + mrSges(5,3) * t222 + Ifges(5,4) * t277 + Ifges(5,2) * t276 + Ifges(5,6) * t313 - pkin(4) * t194 - t305 * t268 + t331 * t270 - t401;
t291 = Ifges(4,4) * t341 + (-Ifges(4,2) * t351 - Ifges(4,6) * t354) * t382;
t363 = mrSges(4,2) * t252 - mrSges(4,3) * t241 + Ifges(4,1) * t340 - Ifges(4,4) * t324 - Ifges(4,5) * t325 - pkin(9) * t181 + t353 * t172 - t350 * t176 + t291 * t374;
t162 = t363 + qJ(3) * (mrSges(4,1) * t325 + t360) + (-t289 * t354 + (-pkin(2) * t321 - t384) * t351) * t382 + Ifges(3,3) * t340 + Ifges(3,5) * t324 + Ifges(3,6) * t325 + mrSges(3,1) * t278 - mrSges(3,2) * t279 + pkin(2) * (-t340 * mrSges(4,2) - t341 * t316 + t366);
t180 = t325 * mrSges(4,2) - t317 * t375 + t370;
t361 = -mrSges(4,1) * t241 + mrSges(4,2) * t242 - pkin(3) * t190 - pkin(9) * t182 - t350 * t172 - t353 * t176;
t164 = -mrSges(3,1) * t293 + mrSges(3,3) * t279 - pkin(2) * t180 + (t289 - t291) * t341 + (Ifges(3,6) - Ifges(4,5)) * t340 + (Ifges(3,2) + Ifges(4,3)) * t325 + t394 * t324 - t385 * t375 + t361;
t362 = -mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t277 - Ifges(5,6) * t276 - Ifges(5,3) * t313 - pkin(4) * t359 - pkin(10) * t373 - t399 * t189 - t349 * t192 - t305 * t269 + t304 * t270;
t358 = -mrSges(4,1) * t252 + mrSges(4,3) * t242 - pkin(3) * t181 + t362;
t166 = -t358 + t384 * t341 + (Ifges(3,5) - Ifges(4,4)) * t340 + t394 * t325 + (Ifges(3,1) + Ifges(4,2)) * t324 + mrSges(3,2) * t293 - mrSges(3,3) * t278 - qJ(3) * t180 + t385 * t374;
t365 = mrSges(2,1) * t336 - mrSges(2,2) * t337 + Ifges(2,3) * qJDD(1) + pkin(1) * t170 + t348 * t162 + t164 * t391 + t166 * t392 + t175 * t398;
t173 = m(2) * t337 - mrSges(2,1) * t356 - qJDD(1) * mrSges(2,2) + t175;
t169 = t348 * t178 + (t179 * t354 + t185 * t351) * t347;
t167 = m(2) * t336 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t356 + t170;
t160 = -mrSges(2,2) * g(3) - mrSges(2,3) * t336 + Ifges(2,5) * qJDD(1) - t356 * Ifges(2,6) - t351 * t164 + t354 * t166 + (-t169 * t347 - t170 * t348) * pkin(8);
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t337 + t356 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t169 - t347 * t162 + (pkin(8) * t175 + t164 * t354 + t166 * t351) * t348;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t160 - t352 * t159 - pkin(7) * (t167 * t355 + t173 * t352), t160, t166, -t290 * t375 + t363, t172, t192, -t245 * t282 + t367; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t160 + t355 * t159 + pkin(7) * (-t167 * t352 + t173 * t355), t159, t164, Ifges(4,4) * t340 - Ifges(4,2) * t324 - Ifges(4,6) * t325 - t341 * t290 - t292 * t374 + t358, t176, t189, -t283 * t243 - t368; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t365, t365, t162, Ifges(4,5) * t340 - Ifges(4,6) * t324 - Ifges(4,3) * t325 + t341 * t291 + t292 * t375 - t361, -t362, t401, Ifges(7,5) * t239 + Ifges(7,6) * t274 + Ifges(7,3) * t238 + t283 * t245 - t302 * t247 - t371;];
m_new  = t1;
