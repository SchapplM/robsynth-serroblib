% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 09:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:50:28
% EndTime: 2019-05-07 08:52:23
% DurationCPUTime: 39.25s
% Computational Cost: add. (702831->395), mult. (1532336->500), div. (0->0), fcn. (1216007->12), ass. (0->154)
t335 = sin(pkin(6));
t340 = sin(qJ(2));
t342 = cos(qJ(2));
t361 = qJD(1) * qJD(2);
t320 = (-qJDD(1) * t342 + t340 * t361) * t335;
t375 = cos(qJ(3));
t374 = cos(qJ(5));
t373 = pkin(8) * t335;
t337 = cos(pkin(6));
t372 = t337 * g(3);
t371 = -mrSges(6,3) - mrSges(7,2);
t370 = t335 * t340;
t369 = t335 * t342;
t368 = t337 * t340;
t367 = t337 * t342;
t363 = qJD(1) * t335;
t318 = (-pkin(2) * t342 - pkin(9) * t340) * t363;
t330 = qJD(1) * t337 + qJD(2);
t328 = t330 ^ 2;
t329 = qJDD(1) * t337 + qJDD(2);
t362 = qJD(1) * t342;
t341 = sin(qJ(1));
t343 = cos(qJ(1));
t326 = t341 * g(1) - g(2) * t343;
t344 = qJD(1) ^ 2;
t315 = qJDD(1) * pkin(1) + t344 * t373 + t326;
t327 = -g(1) * t343 - g(2) * t341;
t316 = -pkin(1) * t344 + qJDD(1) * t373 + t327;
t364 = t315 * t368 + t342 * t316;
t262 = -t328 * pkin(2) + t329 * pkin(9) + (-g(3) * t340 + t318 * t362) * t335 + t364;
t319 = (qJDD(1) * t340 + t342 * t361) * t335;
t263 = t320 * pkin(2) - t319 * pkin(9) - t372 + (-t315 + (pkin(2) * t340 - pkin(9) * t342) * t330 * qJD(1)) * t335;
t339 = sin(qJ(3));
t234 = t375 * t262 + t339 * t263;
t359 = t340 * t363;
t308 = -t375 * t330 + t339 * t359;
t309 = t339 * t330 + t375 * t359;
t290 = pkin(3) * t308 - qJ(4) * t309;
t312 = qJDD(3) + t320;
t358 = t335 * t362;
t325 = qJD(3) - t358;
t324 = t325 ^ 2;
t215 = -pkin(3) * t324 + qJ(4) * t312 - t290 * t308 + t234;
t288 = -g(3) * t369 + t315 * t367 - t340 * t316;
t261 = -t329 * pkin(2) - t328 * pkin(9) + t318 * t359 - t288;
t286 = qJD(3) * t309 + t319 * t339 - t375 * t329;
t287 = -t308 * qJD(3) + t375 * t319 + t339 * t329;
t220 = (t308 * t325 - t287) * qJ(4) + (t309 * t325 + t286) * pkin(3) + t261;
t334 = sin(pkin(11));
t336 = cos(pkin(11));
t297 = t309 * t336 + t325 * t334;
t210 = -0.2e1 * qJD(4) * t297 - t334 * t215 + t336 * t220;
t269 = t287 * t336 + t312 * t334;
t296 = -t309 * t334 + t325 * t336;
t207 = (t296 * t308 - t269) * pkin(10) + (t296 * t297 + t286) * pkin(4) + t210;
t211 = 0.2e1 * qJD(4) * t296 + t336 * t215 + t334 * t220;
t268 = -t287 * t334 + t312 * t336;
t274 = pkin(4) * t308 - pkin(10) * t297;
t295 = t296 ^ 2;
t209 = -pkin(4) * t295 + pkin(10) * t268 - t274 * t308 + t211;
t338 = sin(qJ(5));
t203 = t338 * t207 + t374 * t209;
t266 = t338 * t296 + t374 * t297;
t230 = t266 * qJD(5) - t374 * t268 + t338 * t269;
t307 = qJD(5) + t308;
t251 = mrSges(6,1) * t307 - mrSges(6,3) * t266;
t265 = -t374 * t296 + t338 * t297;
t284 = qJDD(5) + t286;
t244 = pkin(5) * t265 - qJ(6) * t266;
t306 = t307 ^ 2;
t198 = -pkin(5) * t306 + qJ(6) * t284 + 0.2e1 * qJD(6) * t307 - t244 * t265 + t203;
t252 = -mrSges(7,1) * t307 + mrSges(7,2) * t266;
t360 = m(7) * t198 + t284 * mrSges(7,3) + t307 * t252;
t245 = mrSges(7,1) * t265 - mrSges(7,3) * t266;
t365 = -mrSges(6,1) * t265 - mrSges(6,2) * t266 - t245;
t185 = m(6) * t203 - t284 * mrSges(6,2) + t371 * t230 - t307 * t251 + t365 * t265 + t360;
t202 = t374 * t207 - t338 * t209;
t231 = -t265 * qJD(5) + t338 * t268 + t374 * t269;
t250 = -mrSges(6,2) * t307 - mrSges(6,3) * t265;
t200 = -t284 * pkin(5) - t306 * qJ(6) + t266 * t244 + qJDD(6) - t202;
t249 = -mrSges(7,2) * t265 + mrSges(7,3) * t307;
t355 = -m(7) * t200 + t284 * mrSges(7,1) + t307 * t249;
t187 = m(6) * t202 + t284 * mrSges(6,1) + t371 * t231 + t307 * t250 + t365 * t266 + t355;
t180 = t338 * t185 + t374 * t187;
t270 = -mrSges(5,1) * t296 + mrSges(5,2) * t297;
t272 = -mrSges(5,2) * t308 + mrSges(5,3) * t296;
t178 = m(5) * t210 + mrSges(5,1) * t286 - mrSges(5,3) * t269 - t270 * t297 + t272 * t308 + t180;
t273 = mrSges(5,1) * t308 - mrSges(5,3) * t297;
t356 = t374 * t185 - t187 * t338;
t179 = m(5) * t211 - mrSges(5,2) * t286 + mrSges(5,3) * t268 + t270 * t296 - t273 * t308 + t356;
t176 = -t178 * t334 + t336 * t179;
t291 = mrSges(4,1) * t308 + mrSges(4,2) * t309;
t299 = mrSges(4,1) * t325 - mrSges(4,3) * t309;
t174 = m(4) * t234 - mrSges(4,2) * t312 - mrSges(4,3) * t286 - t291 * t308 - t299 * t325 + t176;
t233 = -t339 * t262 + t375 * t263;
t214 = -t312 * pkin(3) - t324 * qJ(4) + t309 * t290 + qJDD(4) - t233;
t212 = -t268 * pkin(4) - t295 * pkin(10) + t297 * t274 + t214;
t205 = -0.2e1 * qJD(6) * t266 + (t265 * t307 - t231) * qJ(6) + (t266 * t307 + t230) * pkin(5) + t212;
t195 = m(7) * t205 + t230 * mrSges(7,1) - t231 * mrSges(7,3) + t265 * t249 - t266 * t252;
t349 = m(6) * t212 + t230 * mrSges(6,1) + mrSges(6,2) * t231 + t265 * t250 + t251 * t266 + t195;
t190 = -m(5) * t214 + t268 * mrSges(5,1) - mrSges(5,2) * t269 + t296 * t272 - t273 * t297 - t349;
t298 = -mrSges(4,2) * t325 - mrSges(4,3) * t308;
t189 = m(4) * t233 + mrSges(4,1) * t312 - mrSges(4,3) * t287 - t291 * t309 + t298 * t325 + t190;
t169 = t339 * t174 + t375 * t189;
t237 = Ifges(7,4) * t266 + Ifges(7,2) * t307 + Ifges(7,6) * t265;
t366 = -Ifges(6,5) * t266 + Ifges(6,6) * t265 - Ifges(6,3) * t307 - t237;
t289 = -g(3) * t370 + t364;
t313 = mrSges(3,1) * t330 - mrSges(3,3) * t359;
t317 = (-mrSges(3,1) * t342 + mrSges(3,2) * t340) * t363;
t357 = t375 * t174 - t189 * t339;
t167 = m(3) * t289 - mrSges(3,2) * t329 - mrSges(3,3) * t320 - t313 * t330 + t317 * t358 + t357;
t314 = -mrSges(3,2) * t330 + mrSges(3,3) * t358;
t175 = t178 * t336 + t179 * t334;
t348 = -m(4) * t261 - t286 * mrSges(4,1) - mrSges(4,2) * t287 - t308 * t298 - t299 * t309 - t175;
t171 = m(3) * t288 + mrSges(3,1) * t329 - mrSges(3,3) * t319 + t314 * t330 - t317 * t359 + t348;
t161 = t342 * t167 - t171 * t340;
t303 = -t335 * t315 - t372;
t168 = m(3) * t303 + t320 * mrSges(3,1) + t319 * mrSges(3,2) + (t313 * t340 - t314 * t342) * t363 + t169;
t158 = t167 * t368 - t168 * t335 + t171 * t367;
t354 = -mrSges(7,1) * t205 + mrSges(7,2) * t198;
t235 = Ifges(7,5) * t266 + Ifges(7,6) * t307 + Ifges(7,3) * t265;
t352 = mrSges(7,2) * t200 - mrSges(7,3) * t205 + Ifges(7,1) * t231 + Ifges(7,4) * t284 + Ifges(7,5) * t230 + t307 * t235;
t239 = Ifges(7,1) * t266 + Ifges(7,4) * t307 + Ifges(7,5) * t265;
t240 = Ifges(6,1) * t266 - Ifges(6,4) * t265 + Ifges(6,5) * t307;
t181 = -mrSges(6,1) * t212 + mrSges(6,3) * t203 - pkin(5) * t195 + (t239 + t240) * t307 + (Ifges(6,6) - Ifges(7,6)) * t284 + t366 * t266 + (Ifges(6,4) - Ifges(7,5)) * t231 + (-Ifges(6,2) - Ifges(7,3)) * t230 + t354;
t238 = Ifges(6,4) * t266 - Ifges(6,2) * t265 + Ifges(6,6) * t307;
t182 = mrSges(6,2) * t212 - mrSges(6,3) * t202 + Ifges(6,1) * t231 - Ifges(6,4) * t230 + Ifges(6,5) * t284 - qJ(6) * t195 - t307 * t238 + t366 * t265 + t352;
t253 = Ifges(5,5) * t297 + Ifges(5,6) * t296 + Ifges(5,3) * t308;
t255 = Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t308;
t163 = -mrSges(5,1) * t214 + mrSges(5,3) * t211 + Ifges(5,4) * t269 + Ifges(5,2) * t268 + Ifges(5,6) * t286 - pkin(4) * t349 + pkin(10) * t356 + t374 * t181 + t338 * t182 - t297 * t253 + t308 * t255;
t254 = Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t308;
t164 = mrSges(5,2) * t214 - mrSges(5,3) * t210 + Ifges(5,1) * t269 + Ifges(5,4) * t268 + Ifges(5,5) * t286 - pkin(10) * t180 - t338 * t181 + t374 * t182 + t296 * t253 - t308 * t254;
t280 = Ifges(4,5) * t309 - Ifges(4,6) * t308 + Ifges(4,3) * t325;
t281 = Ifges(4,4) * t309 - Ifges(4,2) * t308 + Ifges(4,6) * t325;
t154 = mrSges(4,2) * t261 - mrSges(4,3) * t233 + Ifges(4,1) * t287 - Ifges(4,4) * t286 + Ifges(4,5) * t312 - qJ(4) * t175 - t163 * t334 + t164 * t336 - t280 * t308 - t281 * t325;
t282 = Ifges(4,1) * t309 - Ifges(4,4) * t308 + Ifges(4,5) * t325;
t350 = mrSges(7,1) * t200 - mrSges(7,3) * t198 - Ifges(7,4) * t231 - Ifges(7,2) * t284 - Ifges(7,6) * t230 + t266 * t235 - t265 * t239;
t347 = mrSges(6,2) * t203 - t265 * t240 - qJ(6) * (-t230 * mrSges(7,2) - t265 * t245 + t360) - pkin(5) * (-t231 * mrSges(7,2) - t266 * t245 + t355) - mrSges(6,1) * t202 - t266 * t238 + Ifges(6,6) * t230 - Ifges(6,5) * t231 - Ifges(6,3) * t284 + t350;
t345 = mrSges(5,1) * t210 - mrSges(5,2) * t211 + Ifges(5,5) * t269 + Ifges(5,6) * t268 + pkin(4) * t180 + t297 * t254 - t296 * t255 - t347;
t162 = -t345 - pkin(3) * t175 + t325 * t282 - t309 * t280 + Ifges(4,6) * t312 + Ifges(4,4) * t287 - mrSges(4,1) * t261 + mrSges(4,3) * t234 + (-Ifges(5,3) - Ifges(4,2)) * t286;
t301 = Ifges(3,6) * t330 + (Ifges(3,4) * t340 + Ifges(3,2) * t342) * t363;
t302 = Ifges(3,5) * t330 + (Ifges(3,1) * t340 + Ifges(3,4) * t342) * t363;
t149 = Ifges(3,5) * t319 - Ifges(3,6) * t320 + Ifges(3,3) * t329 + mrSges(3,1) * t288 - mrSges(3,2) * t289 + t339 * t154 + t375 * t162 + pkin(2) * t348 + pkin(9) * t357 + (t301 * t340 - t302 * t342) * t363;
t300 = Ifges(3,3) * t330 + (Ifges(3,5) * t340 + Ifges(3,6) * t342) * t363;
t151 = mrSges(3,2) * t303 - mrSges(3,3) * t288 + Ifges(3,1) * t319 - Ifges(3,4) * t320 + Ifges(3,5) * t329 - pkin(9) * t169 + t375 * t154 - t339 * t162 + t300 * t358 - t330 * t301;
t346 = mrSges(4,1) * t233 - mrSges(4,2) * t234 + Ifges(4,5) * t287 - Ifges(4,6) * t286 + Ifges(4,3) * t312 + pkin(3) * t190 + qJ(4) * t176 + t336 * t163 + t334 * t164 + t309 * t281 + t308 * t282;
t153 = -mrSges(3,1) * t303 + mrSges(3,3) * t289 + Ifges(3,4) * t319 - Ifges(3,2) * t320 + Ifges(3,6) * t329 - pkin(2) * t169 - t300 * t359 + t330 * t302 - t346;
t351 = mrSges(2,1) * t326 - mrSges(2,2) * t327 + Ifges(2,3) * qJDD(1) + pkin(1) * t158 + t337 * t149 + t151 * t370 + t153 * t369 + t161 * t373;
t159 = m(2) * t327 - mrSges(2,1) * t344 - qJDD(1) * mrSges(2,2) + t161;
t157 = t337 * t168 + (t167 * t340 + t171 * t342) * t335;
t155 = m(2) * t326 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t344 + t158;
t147 = -mrSges(2,2) * g(3) - mrSges(2,3) * t326 + Ifges(2,5) * qJDD(1) - t344 * Ifges(2,6) + t342 * t151 - t340 * t153 + (-t157 * t335 - t158 * t337) * pkin(8);
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t327 + t344 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t157 - t335 * t149 + (pkin(8) * t161 + t151 * t340 + t153 * t342) * t337;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t343 * t147 - t341 * t146 - pkin(7) * (t155 * t343 + t159 * t341), t147, t151, t154, t164, t182, -t237 * t265 + t352; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t341 * t147 + t343 * t146 + pkin(7) * (-t155 * t341 + t159 * t343), t146, t153, t162, t163, t181, -t350; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t351, t351, t149, t346, Ifges(5,3) * t286 + t345, -t347, Ifges(7,5) * t231 + Ifges(7,6) * t284 + Ifges(7,3) * t230 + t266 * t237 - t307 * t239 - t354;];
m_new  = t1;
