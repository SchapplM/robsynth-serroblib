% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 10:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:24:48
% EndTime: 2019-05-07 10:26:02
% DurationCPUTime: 48.14s
% Computational Cost: add. (910571->387), mult. (1877920->487), div. (0->0), fcn. (1386534->12), ass. (0->153)
t337 = sin(qJ(2));
t341 = cos(qJ(2));
t362 = qJD(1) * qJD(2);
t316 = qJDD(1) * t337 + t341 * t362;
t338 = sin(qJ(1));
t342 = cos(qJ(1));
t323 = -g(1) * t342 - g(2) * t338;
t343 = qJD(1) ^ 2;
t310 = -pkin(1) * t343 + qJDD(1) * pkin(7) + t323;
t365 = t337 * t310;
t366 = pkin(2) * t343;
t270 = qJDD(2) * pkin(2) - t316 * pkin(8) - t365 + (pkin(8) * t362 + t337 * t366 - g(3)) * t341;
t298 = -g(3) * t337 + t341 * t310;
t317 = qJDD(1) * t341 - t337 * t362;
t364 = qJD(1) * t337;
t321 = qJD(2) * pkin(2) - pkin(8) * t364;
t331 = t341 ^ 2;
t271 = pkin(8) * t317 - qJD(2) * t321 - t331 * t366 + t298;
t336 = sin(qJ(3));
t367 = cos(qJ(3));
t250 = t336 * t270 + t367 * t271;
t308 = (t336 * t341 + t367 * t337) * qJD(1);
t280 = qJD(3) * t308 + t316 * t336 - t367 * t317;
t363 = qJD(1) * t341;
t307 = t336 * t364 - t367 * t363;
t291 = mrSges(4,1) * t307 + mrSges(4,2) * t308;
t329 = qJD(2) + qJD(3);
t300 = mrSges(4,1) * t329 - mrSges(4,3) * t308;
t328 = qJDD(2) + qJDD(3);
t281 = -t307 * qJD(3) + t367 * t316 + t336 * t317;
t322 = t338 * g(1) - t342 * g(2);
t354 = -qJDD(1) * pkin(1) - t322;
t282 = -t317 * pkin(2) + t321 * t364 + (-pkin(8) * t331 - pkin(7)) * t343 + t354;
t235 = (t307 * t329 - t281) * qJ(4) + (t308 * t329 + t280) * pkin(3) + t282;
t290 = pkin(3) * t307 - qJ(4) * t308;
t327 = t329 ^ 2;
t238 = -pkin(3) * t327 + qJ(4) * t328 - t290 * t307 + t250;
t332 = sin(pkin(11));
t333 = cos(pkin(11));
t296 = t308 * t333 + t329 * t332;
t218 = -0.2e1 * qJD(4) * t296 + t333 * t235 - t332 * t238;
t264 = t281 * t333 + t328 * t332;
t295 = -t308 * t332 + t329 * t333;
t208 = (t295 * t307 - t264) * pkin(9) + (t295 * t296 + t280) * pkin(4) + t218;
t219 = 0.2e1 * qJD(4) * t295 + t332 * t235 + t333 * t238;
t263 = -t281 * t332 + t328 * t333;
t285 = pkin(4) * t307 - pkin(9) * t296;
t294 = t295 ^ 2;
t210 = -pkin(4) * t294 + pkin(9) * t263 - t285 * t307 + t219;
t335 = sin(qJ(5));
t340 = cos(qJ(5));
t202 = t340 * t208 - t335 * t210;
t261 = t295 * t340 - t296 * t335;
t234 = qJD(5) * t261 + t263 * t335 + t264 * t340;
t262 = t295 * t335 + t296 * t340;
t279 = qJDD(5) + t280;
t303 = qJD(5) + t307;
t199 = (t261 * t303 - t234) * pkin(10) + (t261 * t262 + t279) * pkin(5) + t202;
t203 = t335 * t208 + t340 * t210;
t233 = -qJD(5) * t262 + t263 * t340 - t264 * t335;
t253 = pkin(5) * t303 - pkin(10) * t262;
t260 = t261 ^ 2;
t200 = -pkin(5) * t260 + pkin(10) * t233 - t253 * t303 + t203;
t334 = sin(qJ(6));
t339 = cos(qJ(6));
t197 = t199 * t339 - t200 * t334;
t245 = t261 * t339 - t262 * t334;
t216 = qJD(6) * t245 + t233 * t334 + t234 * t339;
t246 = t261 * t334 + t262 * t339;
t226 = -mrSges(7,1) * t245 + mrSges(7,2) * t246;
t302 = qJD(6) + t303;
t239 = -mrSges(7,2) * t302 + mrSges(7,3) * t245;
t274 = qJDD(6) + t279;
t190 = m(7) * t197 + mrSges(7,1) * t274 - mrSges(7,3) * t216 - t226 * t246 + t239 * t302;
t198 = t199 * t334 + t200 * t339;
t215 = -qJD(6) * t246 + t233 * t339 - t234 * t334;
t240 = mrSges(7,1) * t302 - mrSges(7,3) * t246;
t191 = m(7) * t198 - mrSges(7,2) * t274 + mrSges(7,3) * t215 + t226 * t245 - t240 * t302;
t184 = t339 * t190 + t334 * t191;
t247 = -mrSges(6,1) * t261 + mrSges(6,2) * t262;
t251 = -mrSges(6,2) * t303 + mrSges(6,3) * t261;
t181 = m(6) * t202 + mrSges(6,1) * t279 - mrSges(6,3) * t234 - t247 * t262 + t251 * t303 + t184;
t252 = mrSges(6,1) * t303 - mrSges(6,3) * t262;
t357 = -t190 * t334 + t339 * t191;
t182 = m(6) * t203 - mrSges(6,2) * t279 + mrSges(6,3) * t233 + t247 * t261 - t252 * t303 + t357;
t177 = t340 * t181 + t335 * t182;
t266 = -mrSges(5,1) * t295 + mrSges(5,2) * t296;
t283 = -mrSges(5,2) * t307 + mrSges(5,3) * t295;
t175 = m(5) * t218 + mrSges(5,1) * t280 - mrSges(5,3) * t264 - t266 * t296 + t283 * t307 + t177;
t284 = mrSges(5,1) * t307 - mrSges(5,3) * t296;
t358 = -t181 * t335 + t340 * t182;
t176 = m(5) * t219 - mrSges(5,2) * t280 + mrSges(5,3) * t263 + t266 * t295 - t284 * t307 + t358;
t359 = -t175 * t332 + t333 * t176;
t166 = m(4) * t250 - mrSges(4,2) * t328 - mrSges(4,3) * t280 - t291 * t307 - t300 * t329 + t359;
t249 = t367 * t270 - t336 * t271;
t299 = -mrSges(4,2) * t329 - mrSges(4,3) * t307;
t237 = -t328 * pkin(3) - t327 * qJ(4) + t308 * t290 + qJDD(4) - t249;
t220 = -t263 * pkin(4) - t294 * pkin(9) + t296 * t285 + t237;
t205 = -t233 * pkin(5) - t260 * pkin(10) + t262 * t253 + t220;
t356 = m(7) * t205 - t215 * mrSges(7,1) + t216 * mrSges(7,2) - t245 * t239 + t246 * t240;
t349 = m(6) * t220 - t233 * mrSges(6,1) + mrSges(6,2) * t234 - t261 * t251 + t252 * t262 + t356;
t346 = -m(5) * t237 + t263 * mrSges(5,1) - mrSges(5,2) * t264 + t295 * t283 - t284 * t296 - t349;
t193 = m(4) * t249 + mrSges(4,1) * t328 - mrSges(4,3) * t281 - t291 * t308 + t299 * t329 + t346;
t161 = t336 * t166 + t367 * t193;
t297 = -t341 * g(3) - t365;
t305 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t337 + Ifges(3,2) * t341) * qJD(1);
t306 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t337 + Ifges(3,4) * t341) * qJD(1);
t221 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t302;
t223 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t302;
t185 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t216 + Ifges(7,2) * t215 + Ifges(7,6) * t274 - t221 * t246 + t223 * t302;
t222 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t302;
t186 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t216 + Ifges(7,4) * t215 + Ifges(7,5) * t274 + t221 * t245 - t222 * t302;
t241 = Ifges(6,5) * t262 + Ifges(6,6) * t261 + Ifges(6,3) * t303;
t243 = Ifges(6,1) * t262 + Ifges(6,4) * t261 + Ifges(6,5) * t303;
t170 = -mrSges(6,1) * t220 + mrSges(6,3) * t203 + Ifges(6,4) * t234 + Ifges(6,2) * t233 + Ifges(6,6) * t279 - pkin(5) * t356 + pkin(10) * t357 + t339 * t185 + t334 * t186 - t262 * t241 + t303 * t243;
t242 = Ifges(6,4) * t262 + Ifges(6,2) * t261 + Ifges(6,6) * t303;
t171 = mrSges(6,2) * t220 - mrSges(6,3) * t202 + Ifges(6,1) * t234 + Ifges(6,4) * t233 + Ifges(6,5) * t279 - pkin(10) * t184 - t185 * t334 + t186 * t339 + t241 * t261 - t242 * t303;
t254 = Ifges(5,5) * t296 + Ifges(5,6) * t295 + Ifges(5,3) * t307;
t256 = Ifges(5,1) * t296 + Ifges(5,4) * t295 + Ifges(5,5) * t307;
t155 = -mrSges(5,1) * t237 + mrSges(5,3) * t219 + Ifges(5,4) * t264 + Ifges(5,2) * t263 + Ifges(5,6) * t280 - pkin(4) * t349 + pkin(9) * t358 + t340 * t170 + t335 * t171 - t296 * t254 + t307 * t256;
t255 = Ifges(5,4) * t296 + Ifges(5,2) * t295 + Ifges(5,6) * t307;
t157 = mrSges(5,2) * t237 - mrSges(5,3) * t218 + Ifges(5,1) * t264 + Ifges(5,4) * t263 + Ifges(5,5) * t280 - pkin(9) * t177 - t170 * t335 + t171 * t340 + t254 * t295 - t255 * t307;
t287 = Ifges(4,4) * t308 - Ifges(4,2) * t307 + Ifges(4,6) * t329;
t288 = Ifges(4,1) * t308 - Ifges(4,4) * t307 + Ifges(4,5) * t329;
t350 = -mrSges(4,1) * t249 + mrSges(4,2) * t250 - Ifges(4,5) * t281 + Ifges(4,6) * t280 - Ifges(4,3) * t328 - pkin(3) * t346 - qJ(4) * t359 - t333 * t155 - t332 * t157 - t308 * t287 - t307 * t288;
t368 = mrSges(3,1) * t297 - mrSges(3,2) * t298 + Ifges(3,5) * t316 + Ifges(3,6) * t317 + Ifges(3,3) * qJDD(2) + pkin(2) * t161 + (t305 * t337 - t306 * t341) * qJD(1) - t350;
t168 = t333 * t175 + t332 * t176;
t315 = (-mrSges(3,1) * t341 + mrSges(3,2) * t337) * qJD(1);
t320 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t159 = m(3) * t297 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t316 + qJD(2) * t320 - t315 * t364 + t161;
t319 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t360 = t367 * t166 - t193 * t336;
t160 = m(3) * t298 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t317 - qJD(2) * t319 + t315 * t363 + t360;
t361 = -t159 * t337 + t341 * t160;
t286 = Ifges(4,5) * t308 - Ifges(4,6) * t307 + Ifges(4,3) * t329;
t149 = mrSges(4,2) * t282 - mrSges(4,3) * t249 + Ifges(4,1) * t281 - Ifges(4,4) * t280 + Ifges(4,5) * t328 - qJ(4) * t168 - t155 * t332 + t157 * t333 - t286 * t307 - t287 * t329;
t352 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t216 - Ifges(7,6) * t215 - Ifges(7,3) * t274 - t246 * t222 + t245 * t223;
t348 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t234 - Ifges(6,6) * t233 - Ifges(6,3) * t279 - pkin(5) * t184 - t262 * t242 + t261 * t243 + t352;
t344 = mrSges(5,1) * t218 - mrSges(5,2) * t219 + Ifges(5,5) * t264 + Ifges(5,6) * t263 + pkin(4) * t177 + t296 * t255 - t295 * t256 - t348;
t153 = -pkin(3) * t168 - t344 + (-Ifges(5,3) - Ifges(4,2)) * t280 + mrSges(4,3) * t250 + Ifges(4,4) * t281 - mrSges(4,1) * t282 - t308 * t286 + Ifges(4,6) * t328 + t329 * t288;
t304 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t337 + Ifges(3,6) * t341) * qJD(1);
t309 = -t343 * pkin(7) + t354;
t351 = m(4) * t282 + t280 * mrSges(4,1) + mrSges(4,2) * t281 + t307 * t299 + t300 * t308 + t168;
t145 = -mrSges(3,1) * t309 + mrSges(3,3) * t298 + Ifges(3,4) * t316 + Ifges(3,2) * t317 + Ifges(3,6) * qJDD(2) - pkin(2) * t351 + pkin(8) * t360 + qJD(2) * t306 + t336 * t149 + t367 * t153 - t304 * t364;
t148 = mrSges(3,2) * t309 - mrSges(3,3) * t297 + Ifges(3,1) * t316 + Ifges(3,4) * t317 + Ifges(3,5) * qJDD(2) - pkin(8) * t161 - qJD(2) * t305 + t367 * t149 - t336 * t153 + t304 * t363;
t347 = -m(3) * t309 + t317 * mrSges(3,1) - mrSges(3,2) * t316 - t319 * t364 + t320 * t363 - t351;
t353 = mrSges(2,1) * t322 - mrSges(2,2) * t323 + Ifges(2,3) * qJDD(1) + pkin(1) * t347 + pkin(7) * t361 + t341 * t145 + t337 * t148;
t162 = m(2) * t322 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t343 + t347;
t152 = t159 * t341 + t160 * t337;
t150 = m(2) * t323 - mrSges(2,1) * t343 - qJDD(1) * mrSges(2,2) + t361;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t323 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t152 - t368;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t322 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t343 - pkin(7) * t152 - t145 * t337 + t148 * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t143 - t338 * t146 - pkin(6) * (t150 * t338 + t162 * t342), t143, t148, t149, t157, t171, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t338 * t143 + t342 * t146 + pkin(6) * (t150 * t342 - t162 * t338), t146, t145, t153, t155, t170, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, t368, -t350, Ifges(5,3) * t280 + t344, -t348, -t352;];
m_new  = t1;
