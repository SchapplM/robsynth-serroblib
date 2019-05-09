% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR4
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
% Datum: 2019-05-07 20:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:07:05
% EndTime: 2019-05-07 20:08:41
% DurationCPUTime: 50.61s
% Computational Cost: add. (948725->386), mult. (1915698->487), div. (0->0), fcn. (1417530->12), ass. (0->153)
t339 = sin(qJ(2));
t344 = cos(qJ(2));
t365 = qJD(1) * qJD(2);
t318 = qJDD(1) * t339 + t344 * t365;
t340 = sin(qJ(1));
t345 = cos(qJ(1));
t325 = -g(1) * t345 - g(2) * t340;
t346 = qJD(1) ^ 2;
t312 = -pkin(1) * t346 + qJDD(1) * pkin(7) + t325;
t368 = t339 * t312;
t369 = pkin(2) * t346;
t273 = qJDD(2) * pkin(2) - t318 * pkin(8) - t368 + (pkin(8) * t365 + t339 * t369 - g(3)) * t344;
t299 = -g(3) * t339 + t344 * t312;
t319 = qJDD(1) * t344 - t339 * t365;
t367 = qJD(1) * t339;
t323 = qJD(2) * pkin(2) - pkin(8) * t367;
t333 = t344 ^ 2;
t274 = pkin(8) * t319 - qJD(2) * t323 - t333 * t369 + t299;
t338 = sin(qJ(3));
t343 = cos(qJ(3));
t251 = t338 * t273 + t343 * t274;
t310 = (t338 * t344 + t339 * t343) * qJD(1);
t283 = -t310 * qJD(3) - t338 * t318 + t319 * t343;
t366 = qJD(1) * t344;
t309 = -t338 * t367 + t343 * t366;
t293 = -mrSges(4,1) * t309 + mrSges(4,2) * t310;
t331 = qJD(2) + qJD(3);
t301 = mrSges(4,1) * t331 - mrSges(4,3) * t310;
t330 = qJDD(2) + qJDD(3);
t284 = qJD(3) * t309 + t318 * t343 + t319 * t338;
t324 = t340 * g(1) - t345 * g(2);
t357 = -qJDD(1) * pkin(1) - t324;
t285 = -t319 * pkin(2) + t323 * t367 + (-pkin(8) * t333 - pkin(7)) * t346 + t357;
t231 = (-t309 * t331 - t284) * pkin(9) + (t310 * t331 - t283) * pkin(3) + t285;
t294 = -pkin(3) * t309 - pkin(9) * t310;
t329 = t331 ^ 2;
t239 = -pkin(3) * t329 + pkin(9) * t330 + t294 * t309 + t251;
t337 = sin(qJ(4));
t342 = cos(qJ(4));
t219 = t342 * t231 - t337 * t239;
t296 = -t310 * t337 + t331 * t342;
t254 = qJD(4) * t296 + t284 * t342 + t330 * t337;
t282 = qJDD(4) - t283;
t297 = t310 * t342 + t331 * t337;
t305 = qJD(4) - t309;
t209 = (t296 * t305 - t254) * qJ(5) + (t296 * t297 + t282) * pkin(4) + t219;
t220 = t337 * t231 + t342 * t239;
t253 = -qJD(4) * t297 - t284 * t337 + t330 * t342;
t287 = pkin(4) * t305 - qJ(5) * t297;
t295 = t296 ^ 2;
t211 = -pkin(4) * t295 + qJ(5) * t253 - t287 * t305 + t220;
t334 = sin(pkin(11));
t335 = cos(pkin(11));
t267 = t296 * t334 + t297 * t335;
t203 = -0.2e1 * qJD(5) * t267 + t335 * t209 - t334 * t211;
t236 = t253 * t334 + t254 * t335;
t266 = t296 * t335 - t297 * t334;
t200 = (t266 * t305 - t236) * pkin(10) + (t266 * t267 + t282) * pkin(5) + t203;
t204 = 0.2e1 * qJD(5) * t266 + t334 * t209 + t335 * t211;
t235 = t253 * t335 - t254 * t334;
t257 = pkin(5) * t305 - pkin(10) * t267;
t265 = t266 ^ 2;
t201 = -pkin(5) * t265 + pkin(10) * t235 - t257 * t305 + t204;
t336 = sin(qJ(6));
t341 = cos(qJ(6));
t198 = t200 * t341 - t201 * t336;
t246 = t266 * t341 - t267 * t336;
t217 = qJD(6) * t246 + t235 * t336 + t236 * t341;
t247 = t266 * t336 + t267 * t341;
t227 = -mrSges(7,1) * t246 + mrSges(7,2) * t247;
t303 = qJD(6) + t305;
t240 = -mrSges(7,2) * t303 + mrSges(7,3) * t246;
t277 = qJDD(6) + t282;
t191 = m(7) * t198 + mrSges(7,1) * t277 - mrSges(7,3) * t217 - t227 * t247 + t240 * t303;
t199 = t200 * t336 + t201 * t341;
t216 = -qJD(6) * t247 + t235 * t341 - t236 * t336;
t241 = mrSges(7,1) * t303 - mrSges(7,3) * t247;
t192 = m(7) * t199 - mrSges(7,2) * t277 + mrSges(7,3) * t216 + t227 * t246 - t241 * t303;
t185 = t341 * t191 + t336 * t192;
t248 = -mrSges(6,1) * t266 + mrSges(6,2) * t267;
t255 = -mrSges(6,2) * t305 + mrSges(6,3) * t266;
t182 = m(6) * t203 + mrSges(6,1) * t282 - mrSges(6,3) * t236 - t248 * t267 + t255 * t305 + t185;
t256 = mrSges(6,1) * t305 - mrSges(6,3) * t267;
t360 = -t191 * t336 + t341 * t192;
t183 = m(6) * t204 - mrSges(6,2) * t282 + mrSges(6,3) * t235 + t248 * t266 - t256 * t305 + t360;
t178 = t335 * t182 + t334 * t183;
t271 = -mrSges(5,1) * t296 + mrSges(5,2) * t297;
t286 = -mrSges(5,2) * t305 + mrSges(5,3) * t296;
t176 = m(5) * t219 + mrSges(5,1) * t282 - mrSges(5,3) * t254 - t271 * t297 + t286 * t305 + t178;
t288 = mrSges(5,1) * t305 - mrSges(5,3) * t297;
t361 = -t182 * t334 + t335 * t183;
t177 = m(5) * t220 - mrSges(5,2) * t282 + mrSges(5,3) * t253 + t271 * t296 - t288 * t305 + t361;
t362 = -t176 * t337 + t342 * t177;
t167 = m(4) * t251 - mrSges(4,2) * t330 + mrSges(4,3) * t283 + t293 * t309 - t301 * t331 + t362;
t250 = t273 * t343 - t338 * t274;
t300 = -mrSges(4,2) * t331 + mrSges(4,3) * t309;
t238 = -pkin(3) * t330 - pkin(9) * t329 + t310 * t294 - t250;
t221 = -pkin(4) * t253 - qJ(5) * t295 + t297 * t287 + qJDD(5) + t238;
t206 = -pkin(5) * t235 - pkin(10) * t265 + t257 * t267 + t221;
t359 = m(7) * t206 - t216 * mrSges(7,1) + t217 * mrSges(7,2) - t246 * t240 + t247 * t241;
t352 = m(6) * t221 - t235 * mrSges(6,1) + mrSges(6,2) * t236 - t266 * t255 + t256 * t267 + t359;
t349 = -m(5) * t238 + t253 * mrSges(5,1) - mrSges(5,2) * t254 + t296 * t286 - t288 * t297 - t352;
t194 = m(4) * t250 + mrSges(4,1) * t330 - mrSges(4,3) * t284 - t293 * t310 + t300 * t331 + t349;
t162 = t338 * t167 + t343 * t194;
t298 = -t344 * g(3) - t368;
t307 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t344) * qJD(1);
t308 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t344) * qJD(1);
t222 = Ifges(7,5) * t247 + Ifges(7,6) * t246 + Ifges(7,3) * t303;
t224 = Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t303;
t186 = -mrSges(7,1) * t206 + mrSges(7,3) * t199 + Ifges(7,4) * t217 + Ifges(7,2) * t216 + Ifges(7,6) * t277 - t222 * t247 + t224 * t303;
t223 = Ifges(7,4) * t247 + Ifges(7,2) * t246 + Ifges(7,6) * t303;
t187 = mrSges(7,2) * t206 - mrSges(7,3) * t198 + Ifges(7,1) * t217 + Ifges(7,4) * t216 + Ifges(7,5) * t277 + t222 * t246 - t223 * t303;
t242 = Ifges(6,5) * t267 + Ifges(6,6) * t266 + Ifges(6,3) * t305;
t244 = Ifges(6,1) * t267 + Ifges(6,4) * t266 + Ifges(6,5) * t305;
t171 = -mrSges(6,1) * t221 + mrSges(6,3) * t204 + Ifges(6,4) * t236 + Ifges(6,2) * t235 + Ifges(6,6) * t282 - pkin(5) * t359 + pkin(10) * t360 + t341 * t186 + t336 * t187 - t267 * t242 + t305 * t244;
t243 = Ifges(6,4) * t267 + Ifges(6,2) * t266 + Ifges(6,6) * t305;
t172 = mrSges(6,2) * t221 - mrSges(6,3) * t203 + Ifges(6,1) * t236 + Ifges(6,4) * t235 + Ifges(6,5) * t282 - pkin(10) * t185 - t186 * t336 + t187 * t341 + t242 * t266 - t243 * t305;
t258 = Ifges(5,5) * t297 + Ifges(5,6) * t296 + Ifges(5,3) * t305;
t260 = Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t305;
t156 = -mrSges(5,1) * t238 + mrSges(5,3) * t220 + Ifges(5,4) * t254 + Ifges(5,2) * t253 + Ifges(5,6) * t282 - pkin(4) * t352 + qJ(5) * t361 + t335 * t171 + t334 * t172 - t297 * t258 + t305 * t260;
t259 = Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t305;
t158 = mrSges(5,2) * t238 - mrSges(5,3) * t219 + Ifges(5,1) * t254 + Ifges(5,4) * t253 + Ifges(5,5) * t282 - qJ(5) * t178 - t171 * t334 + t172 * t335 + t258 * t296 - t259 * t305;
t290 = Ifges(4,4) * t310 + Ifges(4,2) * t309 + Ifges(4,6) * t331;
t291 = Ifges(4,1) * t310 + Ifges(4,4) * t309 + Ifges(4,5) * t331;
t353 = -mrSges(4,1) * t250 + mrSges(4,2) * t251 - Ifges(4,5) * t284 - Ifges(4,6) * t283 - Ifges(4,3) * t330 - pkin(3) * t349 - pkin(9) * t362 - t342 * t156 - t337 * t158 - t310 * t290 + t309 * t291;
t370 = mrSges(3,1) * t298 - mrSges(3,2) * t299 + Ifges(3,5) * t318 + Ifges(3,6) * t319 + Ifges(3,3) * qJDD(2) + pkin(2) * t162 + (t307 * t339 - t308 * t344) * qJD(1) - t353;
t169 = t342 * t176 + t337 * t177;
t317 = (-mrSges(3,1) * t344 + mrSges(3,2) * t339) * qJD(1);
t322 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t160 = m(3) * t298 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t318 + qJD(2) * t322 - t317 * t367 + t162;
t321 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t363 = t343 * t167 - t194 * t338;
t161 = m(3) * t299 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t319 - qJD(2) * t321 + t317 * t366 + t363;
t364 = -t160 * t339 + t344 * t161;
t289 = Ifges(4,5) * t310 + Ifges(4,6) * t309 + Ifges(4,3) * t331;
t150 = mrSges(4,2) * t285 - mrSges(4,3) * t250 + Ifges(4,1) * t284 + Ifges(4,4) * t283 + Ifges(4,5) * t330 - pkin(9) * t169 - t156 * t337 + t158 * t342 + t289 * t309 - t290 * t331;
t355 = -mrSges(7,1) * t198 + mrSges(7,2) * t199 - Ifges(7,5) * t217 - Ifges(7,6) * t216 - Ifges(7,3) * t277 - t247 * t223 + t246 * t224;
t351 = -mrSges(6,1) * t203 + mrSges(6,2) * t204 - Ifges(6,5) * t236 - Ifges(6,6) * t235 - Ifges(6,3) * t282 - pkin(5) * t185 - t267 * t243 + t266 * t244 + t355;
t347 = mrSges(5,1) * t219 - mrSges(5,2) * t220 + Ifges(5,5) * t254 + Ifges(5,6) * t253 + Ifges(5,3) * t282 + pkin(4) * t178 + t297 * t259 - t296 * t260 - t351;
t154 = -mrSges(4,1) * t285 + mrSges(4,3) * t251 + Ifges(4,4) * t284 + Ifges(4,2) * t283 + Ifges(4,6) * t330 - pkin(3) * t169 - t310 * t289 + t331 * t291 - t347;
t306 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t344) * qJD(1);
t311 = -t346 * pkin(7) + t357;
t354 = m(4) * t285 - t283 * mrSges(4,1) + mrSges(4,2) * t284 - t309 * t300 + t301 * t310 + t169;
t146 = -mrSges(3,1) * t311 + mrSges(3,3) * t299 + Ifges(3,4) * t318 + Ifges(3,2) * t319 + Ifges(3,6) * qJDD(2) - pkin(2) * t354 + pkin(8) * t363 + qJD(2) * t308 + t338 * t150 + t343 * t154 - t306 * t367;
t149 = mrSges(3,2) * t311 - mrSges(3,3) * t298 + Ifges(3,1) * t318 + Ifges(3,4) * t319 + Ifges(3,5) * qJDD(2) - pkin(8) * t162 - qJD(2) * t307 + t150 * t343 - t154 * t338 + t306 * t366;
t350 = -m(3) * t311 + t319 * mrSges(3,1) - mrSges(3,2) * t318 - t321 * t367 + t322 * t366 - t354;
t356 = mrSges(2,1) * t324 - mrSges(2,2) * t325 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t364 + t344 * t146 + t339 * t149;
t163 = m(2) * t324 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t346 + t350;
t153 = t160 * t344 + t161 * t339;
t151 = m(2) * t325 - mrSges(2,1) * t346 - qJDD(1) * mrSges(2,2) + t364;
t147 = mrSges(2,1) * g(3) + mrSges(2,3) * t325 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t153 - t370;
t144 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t346 - pkin(7) * t153 - t146 * t339 + t149 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t144 - t340 * t147 - pkin(6) * (t151 * t340 + t163 * t345), t144, t149, t150, t158, t172, t187; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t144 + t345 * t147 + pkin(6) * (t151 * t345 - t163 * t340), t147, t146, t154, t156, t171, t186; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t370, -t353, t347, -t351, -t355;];
m_new  = t1;
