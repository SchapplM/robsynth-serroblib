% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 19:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:41:19
% EndTime: 2019-05-07 19:42:20
% DurationCPUTime: 46.26s
% Computational Cost: add. (859609->387), mult. (1884422->488), div. (0->0), fcn. (1429279->12), ass. (0->153)
t337 = sin(qJ(2));
t341 = cos(qJ(2));
t362 = qJD(1) * qJD(2);
t314 = qJDD(1) * t337 + t341 * t362;
t338 = sin(qJ(1));
t342 = cos(qJ(1));
t321 = -g(1) * t342 - g(2) * t338;
t343 = qJD(1) ^ 2;
t309 = -pkin(1) * t343 + qJDD(1) * pkin(7) + t321;
t365 = t309 * t337;
t366 = pkin(2) * t343;
t270 = qJDD(2) * pkin(2) - pkin(8) * t314 - t365 + (pkin(8) * t362 + t337 * t366 - g(3)) * t341;
t296 = -g(3) * t337 + t341 * t309;
t315 = qJDD(1) * t341 - t337 * t362;
t364 = qJD(1) * t337;
t319 = qJD(2) * pkin(2) - pkin(8) * t364;
t331 = t341 ^ 2;
t273 = pkin(8) * t315 - qJD(2) * t319 - t331 * t366 + t296;
t336 = sin(qJ(3));
t340 = cos(qJ(3));
t251 = t340 * t270 - t273 * t336;
t306 = (-t336 * t337 + t340 * t341) * qJD(1);
t282 = qJD(3) * t306 + t314 * t340 + t315 * t336;
t307 = (t336 * t341 + t337 * t340) * qJD(1);
t328 = qJDD(2) + qJDD(3);
t329 = qJD(2) + qJD(3);
t224 = (t306 * t329 - t282) * pkin(9) + (t306 * t307 + t328) * pkin(3) + t251;
t252 = t336 * t270 + t340 * t273;
t281 = -qJD(3) * t307 - t314 * t336 + t315 * t340;
t299 = pkin(3) * t329 - pkin(9) * t307;
t302 = t306 ^ 2;
t228 = -pkin(3) * t302 + pkin(9) * t281 - t299 * t329 + t252;
t335 = sin(qJ(4));
t367 = cos(qJ(4));
t212 = t335 * t224 + t367 * t228;
t293 = t335 * t306 + t367 * t307;
t246 = qJD(4) * t293 - t367 * t281 + t282 * t335;
t292 = -t367 * t306 + t307 * t335;
t265 = mrSges(5,1) * t292 + mrSges(5,2) * t293;
t326 = qJD(4) + t329;
t285 = mrSges(5,1) * t326 - mrSges(5,3) * t293;
t325 = qJDD(4) + t328;
t264 = pkin(4) * t292 - qJ(5) * t293;
t324 = t326 ^ 2;
t206 = -pkin(4) * t324 + qJ(5) * t325 - t264 * t292 + t212;
t320 = g(1) * t338 - t342 * g(2);
t355 = -qJDD(1) * pkin(1) - t320;
t283 = -pkin(2) * t315 + t319 * t364 + (-pkin(8) * t331 - pkin(7)) * t343 + t355;
t233 = -pkin(3) * t281 - pkin(9) * t302 + t307 * t299 + t283;
t247 = -t292 * qJD(4) + t335 * t281 + t367 * t282;
t209 = (t292 * t326 - t247) * qJ(5) + (t293 * t326 + t246) * pkin(4) + t233;
t332 = sin(pkin(11));
t333 = cos(pkin(11));
t277 = t293 * t333 + t326 * t332;
t201 = -0.2e1 * qJD(5) * t277 - t206 * t332 + t333 * t209;
t231 = t247 * t333 + t325 * t332;
t276 = -t293 * t332 + t326 * t333;
t199 = (t276 * t292 - t231) * pkin(10) + (t276 * t277 + t246) * pkin(5) + t201;
t202 = 0.2e1 * qJD(5) * t276 + t333 * t206 + t332 * t209;
t230 = -t247 * t332 + t325 * t333;
t258 = pkin(5) * t292 - pkin(10) * t277;
t274 = t276 ^ 2;
t200 = -pkin(5) * t274 + pkin(10) * t230 - t258 * t292 + t202;
t334 = sin(qJ(6));
t339 = cos(qJ(6));
t197 = t199 * t339 - t200 * t334;
t253 = t276 * t339 - t277 * t334;
t217 = qJD(6) * t253 + t230 * t334 + t231 * t339;
t254 = t276 * t334 + t277 * t339;
t225 = -mrSges(7,1) * t253 + mrSges(7,2) * t254;
t289 = qJD(6) + t292;
t234 = -mrSges(7,2) * t289 + mrSges(7,3) * t253;
t244 = qJDD(6) + t246;
t190 = m(7) * t197 + mrSges(7,1) * t244 - mrSges(7,3) * t217 - t225 * t254 + t234 * t289;
t198 = t199 * t334 + t200 * t339;
t216 = -qJD(6) * t254 + t230 * t339 - t231 * t334;
t235 = mrSges(7,1) * t289 - mrSges(7,3) * t254;
t191 = m(7) * t198 - mrSges(7,2) * t244 + mrSges(7,3) * t216 + t225 * t253 - t235 * t289;
t184 = t339 * t190 + t334 * t191;
t255 = -mrSges(6,1) * t276 + mrSges(6,2) * t277;
t256 = -mrSges(6,2) * t292 + mrSges(6,3) * t276;
t182 = m(6) * t201 + mrSges(6,1) * t246 - mrSges(6,3) * t231 - t255 * t277 + t256 * t292 + t184;
t257 = mrSges(6,1) * t292 - mrSges(6,3) * t277;
t357 = -t190 * t334 + t339 * t191;
t183 = m(6) * t202 - mrSges(6,2) * t246 + mrSges(6,3) * t230 + t255 * t276 - t257 * t292 + t357;
t358 = -t182 * t332 + t333 * t183;
t175 = m(5) * t212 - mrSges(5,2) * t325 - mrSges(5,3) * t246 - t265 * t292 - t285 * t326 + t358;
t211 = t367 * t224 - t335 * t228;
t284 = -mrSges(5,2) * t326 - mrSges(5,3) * t292;
t205 = -t325 * pkin(4) - t324 * qJ(5) + t293 * t264 + qJDD(5) - t211;
t203 = -t230 * pkin(5) - t274 * pkin(10) + t277 * t258 + t205;
t352 = m(7) * t203 - t216 * mrSges(7,1) + mrSges(7,2) * t217 - t253 * t234 + t235 * t254;
t348 = -m(6) * t205 + t230 * mrSges(6,1) - mrSges(6,2) * t231 + t276 * t256 - t257 * t277 - t352;
t193 = m(5) * t211 + mrSges(5,1) * t325 - mrSges(5,3) * t247 - t265 * t293 + t284 * t326 + t348;
t166 = t335 * t175 + t367 * t193;
t294 = -mrSges(4,1) * t306 + mrSges(4,2) * t307;
t297 = -mrSges(4,2) * t329 + mrSges(4,3) * t306;
t163 = m(4) * t251 + mrSges(4,1) * t328 - mrSges(4,3) * t282 - t294 * t307 + t297 * t329 + t166;
t298 = mrSges(4,1) * t329 - mrSges(4,3) * t307;
t359 = t367 * t175 - t193 * t335;
t164 = m(4) * t252 - mrSges(4,2) * t328 + mrSges(4,3) * t281 + t294 * t306 - t298 * t329 + t359;
t158 = t340 * t163 + t336 * t164;
t295 = -g(3) * t341 - t365;
t304 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t337 + Ifges(3,2) * t341) * qJD(1);
t305 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t337 + Ifges(3,4) * t341) * qJD(1);
t287 = Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t329;
t288 = Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t329;
t219 = Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t289;
t221 = Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t289;
t185 = -mrSges(7,1) * t203 + mrSges(7,3) * t198 + Ifges(7,4) * t217 + Ifges(7,2) * t216 + Ifges(7,6) * t244 - t219 * t254 + t221 * t289;
t220 = Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t289;
t186 = mrSges(7,2) * t203 - mrSges(7,3) * t197 + Ifges(7,1) * t217 + Ifges(7,4) * t216 + Ifges(7,5) * t244 + t219 * t253 - t220 * t289;
t237 = Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * t292;
t239 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t292;
t168 = -mrSges(6,1) * t205 + mrSges(6,3) * t202 + Ifges(6,4) * t231 + Ifges(6,2) * t230 + Ifges(6,6) * t246 - pkin(5) * t352 + pkin(10) * t357 + t339 * t185 + t334 * t186 - t277 * t237 + t292 * t239;
t238 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t292;
t170 = mrSges(6,2) * t205 - mrSges(6,3) * t201 + Ifges(6,1) * t231 + Ifges(6,4) * t230 + Ifges(6,5) * t246 - pkin(10) * t184 - t185 * t334 + t186 * t339 + t237 * t276 - t238 * t292;
t260 = Ifges(5,4) * t293 - Ifges(5,2) * t292 + Ifges(5,6) * t326;
t261 = Ifges(5,1) * t293 - Ifges(5,4) * t292 + Ifges(5,5) * t326;
t350 = -mrSges(5,1) * t211 + mrSges(5,2) * t212 - Ifges(5,5) * t247 + Ifges(5,6) * t246 - Ifges(5,3) * t325 - pkin(4) * t348 - qJ(5) * t358 - t333 * t168 - t332 * t170 - t293 * t260 - t292 * t261;
t347 = -mrSges(4,1) * t251 + mrSges(4,2) * t252 - Ifges(4,5) * t282 - Ifges(4,6) * t281 - Ifges(4,3) * t328 - pkin(3) * t166 - t307 * t287 + t306 * t288 + t350;
t368 = mrSges(3,1) * t295 - mrSges(3,2) * t296 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t158 + (t304 * t337 - t305 * t341) * qJD(1) - t347;
t177 = t333 * t182 + t332 * t183;
t363 = qJD(1) * t341;
t313 = (-mrSges(3,1) * t341 + mrSges(3,2) * t337) * qJD(1);
t318 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t156 = m(3) * t295 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t314 + qJD(2) * t318 - t313 * t364 + t158;
t317 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t360 = -t163 * t336 + t340 * t164;
t157 = m(3) * t296 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t315 - qJD(2) * t317 + t313 * t363 + t360;
t361 = -t156 * t337 + t341 * t157;
t354 = m(5) * t233 + t246 * mrSges(5,1) + t247 * mrSges(5,2) + t292 * t284 + t293 * t285 + t177;
t259 = Ifges(5,5) * t293 - Ifges(5,6) * t292 + Ifges(5,3) * t326;
t154 = mrSges(5,2) * t233 - mrSges(5,3) * t211 + Ifges(5,1) * t247 - Ifges(5,4) * t246 + Ifges(5,5) * t325 - qJ(5) * t177 - t168 * t332 + t170 * t333 - t259 * t292 - t260 * t326;
t351 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t217 - Ifges(7,6) * t216 - Ifges(7,3) * t244 - t254 * t220 + t253 * t221;
t345 = -mrSges(6,1) * t201 + mrSges(6,2) * t202 - Ifges(6,5) * t231 - Ifges(6,6) * t230 - pkin(5) * t184 - t277 * t238 + t276 * t239 + t351;
t159 = (-Ifges(5,2) - Ifges(6,3)) * t246 + Ifges(5,6) * t325 + t326 * t261 - t293 * t259 + Ifges(5,4) * t247 - mrSges(5,1) * t233 + mrSges(5,3) * t212 + t345 - pkin(4) * t177;
t286 = Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t329;
t149 = -mrSges(4,1) * t283 + mrSges(4,3) * t252 + Ifges(4,4) * t282 + Ifges(4,2) * t281 + Ifges(4,6) * t328 - pkin(3) * t354 + pkin(9) * t359 + t335 * t154 + t367 * t159 - t307 * t286 + t329 * t288;
t150 = mrSges(4,2) * t283 - mrSges(4,3) * t251 + Ifges(4,1) * t282 + Ifges(4,4) * t281 + Ifges(4,5) * t328 - pkin(9) * t166 + t367 * t154 - t335 * t159 + t306 * t286 - t329 * t287;
t303 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t337 + Ifges(3,6) * t341) * qJD(1);
t308 = -pkin(7) * t343 + t355;
t349 = m(4) * t283 - t281 * mrSges(4,1) + mrSges(4,2) * t282 - t306 * t297 + t298 * t307 + t354;
t145 = -mrSges(3,1) * t308 + mrSges(3,3) * t296 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t349 + pkin(8) * t360 + qJD(2) * t305 + t340 * t149 + t336 * t150 - t303 * t364;
t147 = mrSges(3,2) * t308 - mrSges(3,3) * t295 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - pkin(8) * t158 - qJD(2) * t304 - t149 * t336 + t150 * t340 + t303 * t363;
t346 = -m(3) * t308 + t315 * mrSges(3,1) - mrSges(3,2) * t314 - t317 * t364 + t318 * t363 - t349;
t353 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(7) * t361 + t341 * t145 + t337 * t147;
t171 = m(2) * t320 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t343 + t346;
t153 = t156 * t341 + t157 * t337;
t151 = m(2) * t321 - mrSges(2,1) * t343 - qJDD(1) * mrSges(2,2) + t361;
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t153 - t368;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t343 - pkin(7) * t153 - t145 * t337 + t147 * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t143 - t338 * t148 - pkin(6) * (t151 * t338 + t171 * t342), t143, t147, t150, t154, t170, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t338 * t143 + t342 * t148 + pkin(6) * (t151 * t342 - t171 * t338), t148, t145, t149, t159, t168, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, t368, -t347, -t350, Ifges(6,3) * t246 - t345, -t351;];
m_new  = t1;
