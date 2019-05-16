% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:30:00
% EndTime: 2019-05-06 03:30:55
% DurationCPUTime: 41.81s
% Computational Cost: add. (718398->365), mult. (1717970->455), div. (0->0), fcn. (1373355->12), ass. (0->154)
t341 = qJD(1) ^ 2;
t335 = sin(qJ(1));
t340 = cos(qJ(1));
t314 = -t340 * g(1) - t335 * g(2);
t307 = -t341 * pkin(1) + qJDD(1) * qJ(2) + t314;
t329 = sin(pkin(11));
t330 = cos(pkin(11));
t368 = qJD(1) * qJD(2);
t366 = -t330 * g(3) - 0.2e1 * t329 * t368;
t373 = pkin(2) * t330;
t276 = (-pkin(7) * qJDD(1) + t341 * t373 - t307) * t329 + t366;
t297 = -t329 * g(3) + (t307 + 0.2e1 * t368) * t330;
t367 = qJDD(1) * t330;
t325 = t330 ^ 2;
t371 = t325 * t341;
t277 = -pkin(2) * t371 + pkin(7) * t367 + t297;
t334 = sin(qJ(3));
t339 = cos(qJ(3));
t255 = t339 * t276 - t334 * t277;
t355 = t329 * t339 + t330 * t334;
t354 = -t329 * t334 + t330 * t339;
t305 = t354 * qJD(1);
t369 = t305 * qJD(3);
t295 = t355 * qJDD(1) + t369;
t306 = t355 * qJD(1);
t225 = (-t295 + t369) * pkin(8) + (t305 * t306 + qJDD(3)) * pkin(3) + t255;
t256 = t334 * t276 + t339 * t277;
t294 = -t306 * qJD(3) + t354 * qJDD(1);
t300 = qJD(3) * pkin(3) - t306 * pkin(8);
t304 = t305 ^ 2;
t231 = -t304 * pkin(3) + t294 * pkin(8) - qJD(3) * t300 + t256;
t333 = sin(qJ(4));
t338 = cos(qJ(4));
t218 = t333 * t225 + t338 * t231;
t285 = t333 * t305 + t338 * t306;
t249 = -t285 * qJD(4) + t338 * t294 - t333 * t295;
t284 = t338 * t305 - t333 * t306;
t265 = -t284 * mrSges(5,1) + t285 * mrSges(5,2);
t326 = qJD(3) + qJD(4);
t274 = t326 * mrSges(5,1) - t285 * mrSges(5,3);
t323 = qJDD(3) + qJDD(4);
t266 = -t284 * pkin(4) - t285 * pkin(9);
t322 = t326 ^ 2;
t206 = -t322 * pkin(4) + t323 * pkin(9) + t284 * t266 + t218;
t324 = t329 ^ 2;
t313 = t335 * g(1) - t340 * g(2);
t360 = qJDD(2) - t313;
t293 = (-pkin(1) - t373) * qJDD(1) + (-qJ(2) + (-t324 - t325) * pkin(7)) * t341 + t360;
t240 = -t294 * pkin(3) - t304 * pkin(8) + t306 * t300 + t293;
t250 = t284 * qJD(4) + t333 * t294 + t338 * t295;
t214 = (-t284 * t326 - t250) * pkin(9) + (t285 * t326 - t249) * pkin(4) + t240;
t332 = sin(qJ(5));
t337 = cos(qJ(5));
t201 = -t332 * t206 + t337 * t214;
t269 = -t332 * t285 + t337 * t326;
t228 = t269 * qJD(5) + t337 * t250 + t332 * t323;
t248 = qJDD(5) - t249;
t270 = t337 * t285 + t332 * t326;
t280 = qJD(5) - t284;
t199 = (t269 * t280 - t228) * pkin(10) + (t269 * t270 + t248) * pkin(5) + t201;
t202 = t337 * t206 + t332 * t214;
t227 = -t270 * qJD(5) - t332 * t250 + t337 * t323;
t259 = t280 * pkin(5) - t270 * pkin(10);
t268 = t269 ^ 2;
t200 = -t268 * pkin(5) + t227 * pkin(10) - t280 * t259 + t202;
t331 = sin(qJ(6));
t336 = cos(qJ(6));
t197 = t336 * t199 - t331 * t200;
t251 = t336 * t269 - t331 * t270;
t211 = t251 * qJD(6) + t331 * t227 + t336 * t228;
t252 = t331 * t269 + t336 * t270;
t223 = -t251 * mrSges(7,1) + t252 * mrSges(7,2);
t278 = qJD(6) + t280;
t232 = -t278 * mrSges(7,2) + t251 * mrSges(7,3);
t243 = qJDD(6) + t248;
t192 = m(7) * t197 + t243 * mrSges(7,1) - t211 * mrSges(7,3) - t252 * t223 + t278 * t232;
t198 = t331 * t199 + t336 * t200;
t210 = -t252 * qJD(6) + t336 * t227 - t331 * t228;
t233 = t278 * mrSges(7,1) - t252 * mrSges(7,3);
t193 = m(7) * t198 - t243 * mrSges(7,2) + t210 * mrSges(7,3) + t251 * t223 - t278 * t233;
t184 = t336 * t192 + t331 * t193;
t253 = -t269 * mrSges(6,1) + t270 * mrSges(6,2);
t257 = -t280 * mrSges(6,2) + t269 * mrSges(6,3);
t182 = m(6) * t201 + t248 * mrSges(6,1) - t228 * mrSges(6,3) - t270 * t253 + t280 * t257 + t184;
t258 = t280 * mrSges(6,1) - t270 * mrSges(6,3);
t361 = -t331 * t192 + t336 * t193;
t183 = m(6) * t202 - t248 * mrSges(6,2) + t227 * mrSges(6,3) + t269 * t253 - t280 * t258 + t361;
t362 = -t332 * t182 + t337 * t183;
t175 = m(5) * t218 - t323 * mrSges(5,2) + t249 * mrSges(5,3) + t284 * t265 - t326 * t274 + t362;
t217 = t338 * t225 - t333 * t231;
t273 = -t326 * mrSges(5,2) + t284 * mrSges(5,3);
t205 = -t323 * pkin(4) - t322 * pkin(9) + t285 * t266 - t217;
t203 = -t227 * pkin(5) - t268 * pkin(10) + t270 * t259 + t205;
t350 = m(7) * t203 - t210 * mrSges(7,1) + t211 * mrSges(7,2) - t251 * t232 + t252 * t233;
t346 = -m(6) * t205 + t227 * mrSges(6,1) - t228 * mrSges(6,2) + t269 * t257 - t270 * t258 - t350;
t188 = m(5) * t217 + t323 * mrSges(5,1) - t250 * mrSges(5,3) - t285 * t265 + t326 * t273 + t346;
t166 = t333 * t175 + t338 * t188;
t289 = -t305 * mrSges(4,1) + t306 * mrSges(4,2);
t298 = -qJD(3) * mrSges(4,2) + t305 * mrSges(4,3);
t163 = m(4) * t255 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t295 + qJD(3) * t298 - t289 * t306 + t166;
t299 = qJD(3) * mrSges(4,1) - t306 * mrSges(4,3);
t363 = t338 * t175 - t188 * t333;
t164 = m(4) * t256 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t294 - qJD(3) * t299 + t289 * t305 + t363;
t158 = t339 * t163 + t334 * t164;
t296 = -t329 * t307 + t366;
t282 = Ifges(4,4) * t306 + Ifges(4,2) * t305 + Ifges(4,6) * qJD(3);
t283 = Ifges(4,1) * t306 + Ifges(4,4) * t305 + Ifges(4,5) * qJD(3);
t219 = Ifges(7,5) * t252 + Ifges(7,6) * t251 + Ifges(7,3) * t278;
t221 = Ifges(7,1) * t252 + Ifges(7,4) * t251 + Ifges(7,5) * t278;
t185 = -mrSges(7,1) * t203 + mrSges(7,3) * t198 + Ifges(7,4) * t211 + Ifges(7,2) * t210 + Ifges(7,6) * t243 - t252 * t219 + t278 * t221;
t220 = Ifges(7,4) * t252 + Ifges(7,2) * t251 + Ifges(7,6) * t278;
t186 = mrSges(7,2) * t203 - mrSges(7,3) * t197 + Ifges(7,1) * t211 + Ifges(7,4) * t210 + Ifges(7,5) * t243 + t251 * t219 - t278 * t220;
t234 = Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t280;
t236 = Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t280;
t168 = -mrSges(6,1) * t205 + mrSges(6,3) * t202 + Ifges(6,4) * t228 + Ifges(6,2) * t227 + Ifges(6,6) * t248 - pkin(5) * t350 + pkin(10) * t361 + t336 * t185 + t331 * t186 - t270 * t234 + t280 * t236;
t235 = Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t280;
t170 = mrSges(6,2) * t205 - mrSges(6,3) * t201 + Ifges(6,1) * t228 + Ifges(6,4) * t227 + Ifges(6,5) * t248 - pkin(10) * t184 - t185 * t331 + t186 * t336 + t234 * t269 - t235 * t280;
t261 = Ifges(5,4) * t285 + Ifges(5,2) * t284 + Ifges(5,6) * t326;
t262 = Ifges(5,1) * t285 + Ifges(5,4) * t284 + Ifges(5,5) * t326;
t348 = -mrSges(5,1) * t217 + mrSges(5,2) * t218 - Ifges(5,5) * t250 - Ifges(5,6) * t249 - Ifges(5,3) * t323 - pkin(4) * t346 - pkin(9) * t362 - t337 * t168 - t332 * t170 - t285 * t261 + t284 * t262;
t344 = -mrSges(4,1) * t255 + mrSges(4,2) * t256 - Ifges(4,5) * t295 - Ifges(4,6) * t294 - Ifges(4,3) * qJDD(3) - pkin(3) * t166 - t306 * t282 + t305 * t283 + t348;
t358 = Ifges(3,4) * t329 + Ifges(3,2) * t330;
t359 = Ifges(3,1) * t329 + Ifges(3,4) * t330;
t374 = -mrSges(3,1) * t296 + mrSges(3,2) * t297 - pkin(2) * t158 - (t329 * t358 - t330 * t359) * t341 + t344;
t372 = mrSges(3,2) * t329;
t177 = t337 * t182 + t332 * t183;
t357 = Ifges(3,5) * t329 + Ifges(3,6) * t330;
t370 = t341 * t357;
t353 = mrSges(3,3) * qJDD(1) + t341 * (-mrSges(3,1) * t330 + t372);
t156 = m(3) * t296 - t353 * t329 + t158;
t364 = -t334 * t163 + t339 * t164;
t157 = m(3) * t297 + t353 * t330 + t364;
t365 = -t156 * t329 + t330 * t157;
t352 = m(5) * t240 - t249 * mrSges(5,1) + t250 * mrSges(5,2) - t284 * t273 + t285 * t274 + t177;
t260 = Ifges(5,5) * t285 + Ifges(5,6) * t284 + Ifges(5,3) * t326;
t154 = mrSges(5,2) * t240 - mrSges(5,3) * t217 + Ifges(5,1) * t250 + Ifges(5,4) * t249 + Ifges(5,5) * t323 - pkin(9) * t177 - t168 * t332 + t170 * t337 + t260 * t284 - t261 * t326;
t349 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t211 - Ifges(7,6) * t210 - Ifges(7,3) * t243 - t252 * t220 + t251 * t221;
t343 = mrSges(6,1) * t201 - mrSges(6,2) * t202 + Ifges(6,5) * t228 + Ifges(6,6) * t227 + Ifges(6,3) * t248 + pkin(5) * t184 + t270 * t235 - t269 * t236 - t349;
t159 = -mrSges(5,1) * t240 + mrSges(5,3) * t218 + Ifges(5,4) * t250 + Ifges(5,2) * t249 + Ifges(5,6) * t323 - pkin(4) * t177 - t285 * t260 + t326 * t262 - t343;
t281 = Ifges(4,5) * t306 + Ifges(4,6) * t305 + Ifges(4,3) * qJD(3);
t149 = -mrSges(4,1) * t293 + mrSges(4,3) * t256 + Ifges(4,4) * t295 + Ifges(4,2) * t294 + Ifges(4,6) * qJDD(3) - pkin(3) * t352 + pkin(8) * t363 + qJD(3) * t283 + t333 * t154 + t338 * t159 - t306 * t281;
t150 = mrSges(4,2) * t293 - mrSges(4,3) * t255 + Ifges(4,1) * t295 + Ifges(4,4) * t294 + Ifges(4,5) * qJDD(3) - pkin(8) * t166 - qJD(3) * t282 + t154 * t338 - t159 * t333 + t281 * t305;
t303 = -qJDD(1) * pkin(1) - t341 * qJ(2) + t360;
t347 = m(4) * t293 - t294 * mrSges(4,1) + t295 * mrSges(4,2) - t305 * t298 + t306 * t299 + t352;
t145 = -mrSges(3,1) * t303 + mrSges(3,3) * t297 - pkin(2) * t347 + pkin(7) * t364 + t358 * qJDD(1) + t339 * t149 + t334 * t150 - t329 * t370;
t147 = mrSges(3,2) * t303 - mrSges(3,3) * t296 - pkin(7) * t158 + t359 * qJDD(1) - t334 * t149 + t339 * t150 + t330 * t370;
t345 = -m(3) * t303 + mrSges(3,1) * t367 - t347 + (t324 * t341 + t371) * mrSges(3,3);
t351 = -mrSges(2,2) * t314 + qJ(2) * t365 + t330 * t145 + t329 * t147 + pkin(1) * (-qJDD(1) * t372 + t345) + mrSges(2,1) * t313 + Ifges(2,3) * qJDD(1);
t171 = (mrSges(2,1) - t372) * qJDD(1) - t341 * mrSges(2,2) + m(2) * t313 + t345;
t153 = t156 * t330 + t157 * t329;
t151 = m(2) * t314 - mrSges(2,1) * t341 - qJDD(1) * mrSges(2,2) + t365;
t148 = mrSges(2,1) * g(3) + (Ifges(2,6) - t357) * qJDD(1) + t341 * Ifges(2,5) + mrSges(2,3) * t314 - pkin(1) * t153 + t374;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t313 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t341 - qJ(2) * t153 - t145 * t329 + t147 * t330;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t340 * t143 - t335 * t148 - pkin(6) * (t151 * t335 + t171 * t340), t143, t147, t150, t154, t170, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t143 + t340 * t148 + pkin(6) * (t151 * t340 - t171 * t335), t148, t145, t149, t159, t168, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t351, t351, t357 * qJDD(1) - t374, -t344, -t348, t343, -t349;];
m_new  = t1;
