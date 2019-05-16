% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:12:38
% EndTime: 2019-05-06 18:12:55
% DurationCPUTime: 6.87s
% Computational Cost: add. (108929->378), mult. (221723->453), div. (0->0), fcn. (137702->8), ass. (0->134)
t353 = sin(qJ(1));
t356 = cos(qJ(1));
t327 = -g(1) * t356 - g(2) * t353;
t358 = qJD(1) ^ 2;
t300 = -pkin(1) * t358 + qJDD(1) * pkin(7) + t327;
t352 = sin(qJ(2));
t355 = cos(qJ(2));
t278 = -t355 * g(3) - t352 * t300;
t279 = -g(3) * t352 + t355 * t300;
t291 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t352 - Ifges(4,3) * t355) * qJD(1);
t294 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t352 + Ifges(3,2) * t355) * qJD(1);
t314 = (-mrSges(4,1) * t355 - mrSges(4,3) * t352) * qJD(1);
t381 = qJD(1) * qJD(2);
t379 = t355 * t381;
t316 = qJDD(1) * t352 + t379;
t378 = t352 * t381;
t317 = qJDD(1) * t355 - t378;
t313 = (-pkin(2) * t355 - qJ(3) * t352) * qJD(1);
t357 = qJD(2) ^ 2;
t382 = qJD(1) * t355;
t393 = 2 * qJD(3);
t255 = -pkin(2) * t357 + qJDD(2) * qJ(3) + qJD(2) * t393 + t313 * t382 + t279;
t383 = qJD(1) * t352;
t325 = -qJD(2) * pkin(3) - pkin(8) * t383;
t389 = t355 ^ 2 * t358;
t230 = -pkin(3) * t389 - pkin(8) * t317 + qJD(2) * t325 + t255;
t262 = -qJDD(2) * pkin(2) - qJ(3) * t357 + t313 * t383 + qJDD(3) - t278;
t231 = (-t316 + t379) * pkin(8) + (-t352 * t355 * t358 - qJDD(2)) * pkin(3) + t262;
t351 = sin(qJ(4));
t354 = cos(qJ(4));
t214 = t354 * t230 + t351 * t231;
t298 = (-t351 * t355 + t352 * t354) * qJD(1);
t263 = -qJD(4) * t298 - t316 * t351 - t317 * t354;
t297 = (t351 * t352 + t354 * t355) * qJD(1);
t274 = mrSges(5,1) * t297 + mrSges(5,2) * t298;
t340 = -qJD(2) + qJD(4);
t281 = mrSges(5,1) * t340 - mrSges(5,3) * t298;
t339 = -qJDD(2) + qJDD(4);
t326 = t353 * g(1) - t356 * g(2);
t299 = -qJDD(1) * pkin(1) - t358 * pkin(7) - t326;
t372 = -t317 * pkin(2) + t299 + (-t316 - t379) * qJ(3);
t217 = -pkin(2) * t378 + pkin(3) * t317 - pkin(8) * t389 - t372 + (t325 + t393) * t383;
t264 = -qJD(4) * t297 + t316 * t354 - t317 * t351;
t207 = t217 + (t298 * t340 - t263) * pkin(4) + (t297 * t340 - t264) * pkin(9);
t275 = pkin(4) * t297 - pkin(9) * t298;
t338 = t340 ^ 2;
t210 = -pkin(4) * t338 + pkin(9) * t339 - t275 * t297 + t214;
t350 = sin(qJ(5));
t392 = cos(qJ(5));
t205 = t350 * t207 + t392 * t210;
t277 = t392 * t298 + t350 * t340;
t225 = qJD(5) * t277 + t264 * t350 - t392 * t339;
t261 = qJDD(5) - t263;
t290 = qJD(5) + t297;
t267 = mrSges(6,1) * t290 - mrSges(6,3) * t277;
t276 = t298 * t350 - t392 * t340;
t243 = pkin(5) * t276 - qJ(6) * t277;
t286 = t290 ^ 2;
t200 = -pkin(5) * t286 + qJ(6) * t261 + 0.2e1 * qJD(6) * t290 - t243 * t276 + t205;
t268 = -mrSges(7,1) * t290 + mrSges(7,2) * t277;
t380 = m(7) * t200 + t261 * mrSges(7,3) + t290 * t268;
t244 = mrSges(7,1) * t276 - mrSges(7,3) * t277;
t386 = -mrSges(6,1) * t276 - mrSges(6,2) * t277 - t244;
t390 = -mrSges(6,3) - mrSges(7,2);
t190 = m(6) * t205 - mrSges(6,2) * t261 + t390 * t225 - t267 * t290 + t386 * t276 + t380;
t204 = t392 * t207 - t350 * t210;
t226 = -t276 * qJD(5) + t392 * t264 + t350 * t339;
t266 = -mrSges(6,2) * t290 - mrSges(6,3) * t276;
t202 = -t261 * pkin(5) - t286 * qJ(6) + t277 * t243 + qJDD(6) - t204;
t265 = -mrSges(7,2) * t276 + mrSges(7,3) * t290;
t375 = -m(7) * t202 + t261 * mrSges(7,1) + t290 * t265;
t192 = m(6) * t204 + mrSges(6,1) * t261 + t390 * t226 + t266 * t290 + t386 * t277 + t375;
t376 = t392 * t190 - t192 * t350;
t178 = m(5) * t214 - mrSges(5,2) * t339 + mrSges(5,3) * t263 - t274 * t297 - t281 * t340 + t376;
t213 = -t351 * t230 + t231 * t354;
t280 = -mrSges(5,2) * t340 - mrSges(5,3) * t297;
t209 = -pkin(4) * t339 - pkin(9) * t338 + t298 * t275 - t213;
t203 = -0.2e1 * qJD(6) * t277 + (t276 * t290 - t226) * qJ(6) + (t277 * t290 + t225) * pkin(5) + t209;
t197 = m(7) * t203 + t225 * mrSges(7,1) - t226 * mrSges(7,3) + t265 * t276 - t277 * t268;
t363 = -m(6) * t209 - t225 * mrSges(6,1) - t226 * mrSges(6,2) - t276 * t266 - t267 * t277 - t197;
t187 = m(5) * t213 + mrSges(5,1) * t339 - mrSges(5,3) * t264 - t274 * t298 + t280 * t340 + t363;
t172 = t178 * t351 + t187 * t354;
t236 = Ifges(7,1) * t277 + Ifges(7,4) * t290 + Ifges(7,5) * t276;
t237 = Ifges(6,1) * t277 - Ifges(6,4) * t276 + Ifges(6,5) * t290;
t374 = -mrSges(7,1) * t203 + mrSges(7,2) * t200;
t234 = Ifges(7,4) * t277 + Ifges(7,2) * t290 + Ifges(7,6) * t276;
t388 = -Ifges(6,5) * t277 + Ifges(6,6) * t276 - Ifges(6,3) * t290 - t234;
t181 = -mrSges(6,1) * t209 + mrSges(6,3) * t205 - pkin(5) * t197 + (t236 + t237) * t290 + t388 * t277 + (Ifges(6,6) - Ifges(7,6)) * t261 + (Ifges(6,4) - Ifges(7,5)) * t226 + (-Ifges(6,2) - Ifges(7,3)) * t225 + t374;
t235 = Ifges(6,4) * t277 - Ifges(6,2) * t276 + Ifges(6,6) * t290;
t232 = Ifges(7,5) * t277 + Ifges(7,6) * t290 + Ifges(7,3) * t276;
t370 = mrSges(7,2) * t202 - mrSges(7,3) * t203 + Ifges(7,1) * t226 + Ifges(7,4) * t261 + Ifges(7,5) * t225 + t290 * t232;
t183 = mrSges(6,2) * t209 - mrSges(6,3) * t204 + Ifges(6,1) * t226 - Ifges(6,4) * t225 + Ifges(6,5) * t261 - qJ(6) * t197 - t235 * t290 + t388 * t276 + t370;
t270 = Ifges(5,4) * t298 - Ifges(5,2) * t297 + Ifges(5,6) * t340;
t271 = Ifges(5,1) * t298 - Ifges(5,4) * t297 + Ifges(5,5) * t340;
t366 = mrSges(5,1) * t213 - mrSges(5,2) * t214 + Ifges(5,5) * t264 + Ifges(5,6) * t263 + Ifges(5,3) * t339 + pkin(4) * t363 + pkin(9) * t376 + t392 * t181 + t350 * t183 + t298 * t270 + t297 * t271;
t362 = -mrSges(4,1) * t262 + mrSges(4,3) * t255 + Ifges(4,4) * t316 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t317 - pkin(3) * t172 - t366;
t324 = mrSges(4,2) * t382 + qJD(2) * mrSges(4,3);
t367 = -m(4) * t262 + qJDD(2) * mrSges(4,1) + qJD(2) * t324 - t172;
t173 = t354 * t178 - t187 * t351;
t322 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t383;
t369 = m(4) * t255 + qJDD(2) * mrSges(4,3) + qJD(2) * t322 + t314 * t382 + t173;
t295 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t352 - Ifges(4,5) * t355) * qJD(1);
t384 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t352 + Ifges(3,4) * t355) * qJD(1) + t295;
t396 = -((t291 - t294) * t352 + t384 * t355) * qJD(1) + mrSges(3,1) * t278 - mrSges(3,2) * t279 + Ifges(3,5) * t316 + Ifges(3,6) * t317 + Ifges(3,3) * qJDD(2) + pkin(2) * (-mrSges(4,2) * t316 - t314 * t383 + t367) + qJ(3) * (mrSges(4,2) * t317 + t369) + t362;
t371 = mrSges(7,1) * t202 - mrSges(7,3) * t200 - Ifges(7,4) * t226 - Ifges(7,2) * t261 - Ifges(7,6) * t225 - t276 * t236;
t395 = -(-t235 + t232) * t277 + mrSges(6,1) * t204 - mrSges(6,2) * t205 + Ifges(6,5) * t226 - Ifges(6,6) * t225 + Ifges(6,3) * t261 + pkin(5) * (-t226 * mrSges(7,2) - t244 * t277 + t375) + qJ(6) * (-t225 * mrSges(7,2) - t244 * t276 + t380) + t276 * t237 - t371;
t391 = mrSges(3,3) + mrSges(4,2);
t185 = t350 * t190 + t392 * t192;
t315 = (-mrSges(3,1) * t355 + mrSges(3,2) * t352) * qJD(1);
t321 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t383;
t168 = m(3) * t279 - qJDD(2) * mrSges(3,2) - qJD(2) * t321 + t315 * t382 + t391 * t317 + t369;
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t382;
t169 = m(3) * t278 + qJDD(2) * mrSges(3,1) + qJD(2) * t323 - t391 * t316 + (-t314 - t315) * t383 + t367;
t377 = t355 * t168 - t169 * t352;
t179 = -m(5) * t217 + t263 * mrSges(5,1) - t264 * mrSges(5,2) - t297 * t280 - t298 * t281 - t185;
t242 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t383 + t372;
t176 = m(4) * t242 - mrSges(4,1) * t317 - t316 * mrSges(4,3) - t322 * t383 - t324 * t382 + t179;
t292 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t352 + Ifges(3,6) * t355) * qJD(1);
t293 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t352 - Ifges(4,6) * t355) * qJD(1);
t269 = Ifges(5,5) * t298 - Ifges(5,6) * t297 + Ifges(5,3) * t340;
t165 = mrSges(5,2) * t217 - mrSges(5,3) * t213 + Ifges(5,1) * t264 + Ifges(5,4) * t263 + Ifges(5,5) * t339 - pkin(9) * t185 - t350 * t181 + t392 * t183 - t297 * t269 - t340 * t270;
t166 = -mrSges(5,1) * t217 + mrSges(5,3) * t214 + Ifges(5,4) * t264 + Ifges(5,2) * t263 + Ifges(5,6) * t339 - pkin(4) * t185 - t298 * t269 + t340 * t271 - t395;
t364 = -mrSges(4,1) * t242 + mrSges(4,2) * t255 - pkin(3) * t179 - pkin(8) * t173 - t351 * t165 - t354 * t166;
t158 = -mrSges(3,1) * t299 + mrSges(3,3) * t279 - pkin(2) * t176 + (Ifges(3,2) + Ifges(4,3)) * t317 + (Ifges(3,4) - Ifges(4,5)) * t316 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t384 * qJD(2) + (-t292 - t293) * t383 + t364;
t365 = mrSges(4,2) * t262 - mrSges(4,3) * t242 + Ifges(4,1) * t316 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t317 - pkin(8) * t172 + qJD(2) * t291 + t354 * t165 - t166 * t351 + t293 * t382;
t160 = mrSges(3,2) * t299 - mrSges(3,3) * t278 + Ifges(3,1) * t316 + Ifges(3,4) * t317 + Ifges(3,5) * qJDD(2) - qJ(3) * t176 - qJD(2) * t294 + t292 * t382 + t365;
t361 = -m(3) * t299 + t317 * mrSges(3,1) - mrSges(3,2) * t316 - t321 * t383 + t323 * t382 - t176;
t368 = mrSges(2,1) * t326 - mrSges(2,2) * t327 + Ifges(2,3) * qJDD(1) + pkin(1) * t361 + pkin(7) * t377 + t355 * t158 + t352 * t160;
t174 = m(2) * t326 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t358 + t361;
t163 = t168 * t352 + t169 * t355;
t161 = m(2) * t327 - mrSges(2,1) * t358 - qJDD(1) * mrSges(2,2) + t377;
t156 = mrSges(2,1) * g(3) + mrSges(2,3) * t327 + t358 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t163 - t396;
t155 = -mrSges(2,2) * g(3) - mrSges(2,3) * t326 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t358 - pkin(7) * t163 - t158 * t352 + t160 * t355;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t356 * t155 - t353 * t156 - pkin(6) * (t161 * t353 + t174 * t356), t155, t160, t365, t165, t183, -t234 * t276 + t370; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t353 * t155 + t356 * t156 + pkin(6) * (t161 * t356 - t174 * t353), t156, t158, t362 + (-t352 * t291 - t355 * t295) * qJD(1), t166, t181, -t277 * t232 - t371; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t368, t368, t396, Ifges(4,5) * t316 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t317 - qJD(2) * t295 + t293 * t383 - t364, t366, t395, Ifges(7,5) * t226 + Ifges(7,6) * t261 + Ifges(7,3) * t225 + t234 * t277 - t290 * t236 - t374;];
m_new  = t1;
