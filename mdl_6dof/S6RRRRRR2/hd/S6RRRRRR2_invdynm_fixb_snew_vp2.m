% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 08:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:24:04
% EndTime: 2019-05-08 08:26:05
% DurationCPUTime: 49.76s
% Computational Cost: add. (895136->387), mult. (1929489->486), div. (0->0), fcn. (1473427->12), ass. (0->155)
t338 = sin(qJ(2));
t344 = cos(qJ(2));
t365 = qJD(1) * qJD(2);
t316 = qJDD(1) * t338 + t344 * t365;
t339 = sin(qJ(1));
t345 = cos(qJ(1));
t323 = -g(1) * t345 - g(2) * t339;
t346 = qJD(1) ^ 2;
t311 = -pkin(1) * t346 + qJDD(1) * pkin(7) + t323;
t368 = t338 * t311;
t369 = pkin(2) * t346;
t272 = qJDD(2) * pkin(2) - t316 * pkin(8) - t368 + (pkin(8) * t365 + t338 * t369 - g(3)) * t344;
t298 = -g(3) * t338 + t344 * t311;
t317 = qJDD(1) * t344 - t338 * t365;
t367 = qJD(1) * t338;
t321 = qJD(2) * pkin(2) - pkin(8) * t367;
t333 = t344 ^ 2;
t273 = pkin(8) * t317 - qJD(2) * t321 - t333 * t369 + t298;
t337 = sin(qJ(3));
t343 = cos(qJ(3));
t253 = t343 * t272 - t337 * t273;
t308 = (-t337 * t338 + t343 * t344) * qJD(1);
t282 = qJD(3) * t308 + t316 * t343 + t317 * t337;
t309 = (t337 * t344 + t338 * t343) * qJD(1);
t330 = qJDD(2) + qJDD(3);
t331 = qJD(2) + qJD(3);
t225 = (t308 * t331 - t282) * pkin(9) + (t308 * t309 + t330) * pkin(3) + t253;
t254 = t337 * t272 + t343 * t273;
t281 = -qJD(3) * t309 - t316 * t337 + t317 * t343;
t301 = pkin(3) * t331 - pkin(9) * t309;
t304 = t308 ^ 2;
t230 = -pkin(3) * t304 + pkin(9) * t281 - t301 * t331 + t254;
t336 = sin(qJ(4));
t342 = cos(qJ(4));
t218 = t336 * t225 + t342 * t230;
t295 = t308 * t336 + t309 * t342;
t248 = -t295 * qJD(4) + t281 * t342 - t336 * t282;
t294 = t308 * t342 - t336 * t309;
t266 = -mrSges(5,1) * t294 + mrSges(5,2) * t295;
t328 = qJD(4) + t331;
t285 = mrSges(5,1) * t328 - mrSges(5,3) * t295;
t327 = qJDD(4) + t330;
t267 = -pkin(4) * t294 - pkin(10) * t295;
t326 = t328 ^ 2;
t207 = -pkin(4) * t326 + pkin(10) * t327 + t267 * t294 + t218;
t322 = t339 * g(1) - t345 * g(2);
t358 = -qJDD(1) * pkin(1) - t322;
t283 = -t317 * pkin(2) + t321 * t367 + (-pkin(8) * t333 - pkin(7)) * t346 + t358;
t234 = -t281 * pkin(3) - t304 * pkin(9) + t309 * t301 + t283;
t249 = qJD(4) * t294 + t281 * t336 + t282 * t342;
t215 = (-t294 * t328 - t249) * pkin(10) + (t295 * t328 - t248) * pkin(4) + t234;
t335 = sin(qJ(5));
t341 = cos(qJ(5));
t202 = -t335 * t207 + t341 * t215;
t276 = -t295 * t335 + t328 * t341;
t232 = qJD(5) * t276 + t249 * t341 + t327 * t335;
t246 = qJDD(5) - t248;
t277 = t295 * t341 + t328 * t335;
t291 = qJD(5) - t294;
t200 = (t276 * t291 - t232) * pkin(11) + (t276 * t277 + t246) * pkin(5) + t202;
t203 = t341 * t207 + t335 * t215;
t231 = -qJD(5) * t277 - t249 * t335 + t327 * t341;
t260 = pkin(5) * t291 - pkin(11) * t277;
t275 = t276 ^ 2;
t201 = -pkin(5) * t275 + pkin(11) * t231 - t260 * t291 + t203;
t334 = sin(qJ(6));
t340 = cos(qJ(6));
t198 = t200 * t340 - t201 * t334;
t255 = t276 * t340 - t277 * t334;
t214 = qJD(6) * t255 + t231 * t334 + t232 * t340;
t256 = t276 * t334 + t277 * t340;
t226 = -mrSges(7,1) * t255 + mrSges(7,2) * t256;
t286 = qJD(6) + t291;
t235 = -mrSges(7,2) * t286 + mrSges(7,3) * t255;
t238 = qJDD(6) + t246;
t193 = m(7) * t198 + mrSges(7,1) * t238 - mrSges(7,3) * t214 - t226 * t256 + t235 * t286;
t199 = t200 * t334 + t201 * t340;
t213 = -qJD(6) * t256 + t231 * t340 - t232 * t334;
t236 = mrSges(7,1) * t286 - mrSges(7,3) * t256;
t194 = m(7) * t199 - mrSges(7,2) * t238 + mrSges(7,3) * t213 + t226 * t255 - t236 * t286;
t185 = t340 * t193 + t334 * t194;
t257 = -mrSges(6,1) * t276 + mrSges(6,2) * t277;
t258 = -mrSges(6,2) * t291 + mrSges(6,3) * t276;
t183 = m(6) * t202 + mrSges(6,1) * t246 - mrSges(6,3) * t232 - t257 * t277 + t258 * t291 + t185;
t259 = mrSges(6,1) * t291 - mrSges(6,3) * t277;
t360 = -t193 * t334 + t340 * t194;
t184 = m(6) * t203 - mrSges(6,2) * t246 + mrSges(6,3) * t231 + t257 * t276 - t259 * t291 + t360;
t361 = -t183 * t335 + t341 * t184;
t176 = m(5) * t218 - mrSges(5,2) * t327 + mrSges(5,3) * t248 + t266 * t294 - t285 * t328 + t361;
t217 = t225 * t342 - t336 * t230;
t284 = -mrSges(5,2) * t328 + mrSges(5,3) * t294;
t206 = -pkin(4) * t327 - pkin(10) * t326 + t295 * t267 - t217;
t204 = -pkin(5) * t231 - pkin(11) * t275 + t260 * t277 + t206;
t355 = m(7) * t204 - t213 * mrSges(7,1) + mrSges(7,2) * t214 - t255 * t235 + t236 * t256;
t351 = -m(6) * t206 + t231 * mrSges(6,1) - mrSges(6,2) * t232 + t276 * t258 - t259 * t277 - t355;
t189 = m(5) * t217 + mrSges(5,1) * t327 - mrSges(5,3) * t249 - t266 * t295 + t284 * t328 + t351;
t167 = t336 * t176 + t342 * t189;
t296 = -mrSges(4,1) * t308 + mrSges(4,2) * t309;
t299 = -mrSges(4,2) * t331 + mrSges(4,3) * t308;
t164 = m(4) * t253 + mrSges(4,1) * t330 - mrSges(4,3) * t282 - t296 * t309 + t299 * t331 + t167;
t300 = mrSges(4,1) * t331 - mrSges(4,3) * t309;
t362 = t342 * t176 - t189 * t336;
t165 = m(4) * t254 - mrSges(4,2) * t330 + mrSges(4,3) * t281 + t296 * t308 - t300 * t331 + t362;
t159 = t343 * t164 + t337 * t165;
t297 = -t344 * g(3) - t368;
t306 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t338 + Ifges(3,2) * t344) * qJD(1);
t307 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t338 + Ifges(3,4) * t344) * qJD(1);
t288 = Ifges(4,4) * t309 + Ifges(4,2) * t308 + Ifges(4,6) * t331;
t289 = Ifges(4,1) * t309 + Ifges(4,4) * t308 + Ifges(4,5) * t331;
t220 = Ifges(7,5) * t256 + Ifges(7,6) * t255 + Ifges(7,3) * t286;
t222 = Ifges(7,1) * t256 + Ifges(7,4) * t255 + Ifges(7,5) * t286;
t186 = -mrSges(7,1) * t204 + mrSges(7,3) * t199 + Ifges(7,4) * t214 + Ifges(7,2) * t213 + Ifges(7,6) * t238 - t220 * t256 + t222 * t286;
t221 = Ifges(7,4) * t256 + Ifges(7,2) * t255 + Ifges(7,6) * t286;
t187 = mrSges(7,2) * t204 - mrSges(7,3) * t198 + Ifges(7,1) * t214 + Ifges(7,4) * t213 + Ifges(7,5) * t238 + t220 * t255 - t221 * t286;
t239 = Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * t291;
t241 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t291;
t169 = -mrSges(6,1) * t206 + mrSges(6,3) * t203 + Ifges(6,4) * t232 + Ifges(6,2) * t231 + Ifges(6,6) * t246 - pkin(5) * t355 + pkin(11) * t360 + t340 * t186 + t334 * t187 - t277 * t239 + t291 * t241;
t240 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t291;
t171 = mrSges(6,2) * t206 - mrSges(6,3) * t202 + Ifges(6,1) * t232 + Ifges(6,4) * t231 + Ifges(6,5) * t246 - pkin(11) * t185 - t186 * t334 + t187 * t340 + t239 * t276 - t240 * t291;
t262 = Ifges(5,4) * t295 + Ifges(5,2) * t294 + Ifges(5,6) * t328;
t263 = Ifges(5,1) * t295 + Ifges(5,4) * t294 + Ifges(5,5) * t328;
t353 = -mrSges(5,1) * t217 + mrSges(5,2) * t218 - Ifges(5,5) * t249 - Ifges(5,6) * t248 - Ifges(5,3) * t327 - pkin(4) * t351 - pkin(10) * t361 - t341 * t169 - t335 * t171 - t295 * t262 + t294 * t263;
t350 = -mrSges(4,1) * t253 + mrSges(4,2) * t254 - Ifges(4,5) * t282 - Ifges(4,6) * t281 - Ifges(4,3) * t330 - pkin(3) * t167 - t309 * t288 + t308 * t289 + t353;
t370 = mrSges(3,1) * t297 - mrSges(3,2) * t298 + Ifges(3,5) * t316 + Ifges(3,6) * t317 + Ifges(3,3) * qJDD(2) + pkin(2) * t159 + (t306 * t338 - t307 * t344) * qJD(1) - t350;
t178 = t341 * t183 + t335 * t184;
t366 = qJD(1) * t344;
t315 = (-mrSges(3,1) * t344 + mrSges(3,2) * t338) * qJD(1);
t320 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t157 = m(3) * t297 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t316 + qJD(2) * t320 - t315 * t367 + t159;
t319 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t363 = -t164 * t337 + t343 * t165;
t158 = m(3) * t298 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t317 - qJD(2) * t319 + t315 * t366 + t363;
t364 = -t157 * t338 + t344 * t158;
t357 = m(5) * t234 - t248 * mrSges(5,1) + t249 * mrSges(5,2) - t294 * t284 + t295 * t285 + t178;
t261 = Ifges(5,5) * t295 + Ifges(5,6) * t294 + Ifges(5,3) * t328;
t155 = mrSges(5,2) * t234 - mrSges(5,3) * t217 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t327 - pkin(10) * t178 - t169 * t335 + t171 * t341 + t261 * t294 - t262 * t328;
t354 = -mrSges(7,1) * t198 + mrSges(7,2) * t199 - Ifges(7,5) * t214 - Ifges(7,6) * t213 - Ifges(7,3) * t238 - t256 * t221 + t255 * t222;
t348 = mrSges(6,1) * t202 - mrSges(6,2) * t203 + Ifges(6,5) * t232 + Ifges(6,6) * t231 + Ifges(6,3) * t246 + pkin(5) * t185 + t277 * t240 - t276 * t241 - t354;
t160 = -mrSges(5,1) * t234 + mrSges(5,3) * t218 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t327 - pkin(4) * t178 - t295 * t261 + t328 * t263 - t348;
t287 = Ifges(4,5) * t309 + Ifges(4,6) * t308 + Ifges(4,3) * t331;
t150 = -mrSges(4,1) * t283 + mrSges(4,3) * t254 + Ifges(4,4) * t282 + Ifges(4,2) * t281 + Ifges(4,6) * t330 - pkin(3) * t357 + pkin(9) * t362 + t336 * t155 + t342 * t160 - t309 * t287 + t331 * t289;
t151 = mrSges(4,2) * t283 - mrSges(4,3) * t253 + Ifges(4,1) * t282 + Ifges(4,4) * t281 + Ifges(4,5) * t330 - pkin(9) * t167 + t155 * t342 - t160 * t336 + t287 * t308 - t288 * t331;
t305 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t338 + Ifges(3,6) * t344) * qJD(1);
t310 = -t346 * pkin(7) + t358;
t352 = m(4) * t283 - t281 * mrSges(4,1) + mrSges(4,2) * t282 - t308 * t299 + t300 * t309 + t357;
t146 = -mrSges(3,1) * t310 + mrSges(3,3) * t298 + Ifges(3,4) * t316 + Ifges(3,2) * t317 + Ifges(3,6) * qJDD(2) - pkin(2) * t352 + pkin(8) * t363 + qJD(2) * t307 + t343 * t150 + t337 * t151 - t305 * t367;
t148 = mrSges(3,2) * t310 - mrSges(3,3) * t297 + Ifges(3,1) * t316 + Ifges(3,4) * t317 + Ifges(3,5) * qJDD(2) - pkin(8) * t159 - qJD(2) * t306 - t150 * t337 + t151 * t343 + t305 * t366;
t349 = -m(3) * t310 + t317 * mrSges(3,1) - mrSges(3,2) * t316 - t319 * t367 + t320 * t366 - t352;
t356 = mrSges(2,1) * t322 - mrSges(2,2) * t323 + Ifges(2,3) * qJDD(1) + pkin(1) * t349 + pkin(7) * t364 + t344 * t146 + t338 * t148;
t172 = m(2) * t322 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t346 + t349;
t154 = t157 * t344 + t158 * t338;
t152 = m(2) * t323 - mrSges(2,1) * t346 - qJDD(1) * mrSges(2,2) + t364;
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t323 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t154 - t370;
t144 = -mrSges(2,2) * g(3) - mrSges(2,3) * t322 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t346 - pkin(7) * t154 - t146 * t338 + t148 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t144 - t339 * t149 - pkin(6) * (t152 * t339 + t172 * t345), t144, t148, t151, t155, t171, t187; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t339 * t144 + t345 * t149 + pkin(6) * (t152 * t345 - t172 * t339), t149, t146, t150, t160, t169, t186; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t370, -t350, -t353, t348, -t354;];
m_new  = t1;
