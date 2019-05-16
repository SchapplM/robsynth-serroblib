% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR1
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
% Datum: 2019-05-08 08:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:57:53
% EndTime: 2019-05-08 07:59:31
% DurationCPUTime: 56.13s
% Computational Cost: add. (980517->386), mult. (2250970->486), div. (0->0), fcn. (1740209->12), ass. (0->155)
t335 = sin(qJ(2));
t341 = cos(qJ(2));
t362 = qJD(1) * qJD(2);
t310 = qJDD(1) * t335 + t341 * t362;
t336 = sin(qJ(1));
t342 = cos(qJ(1));
t317 = -g(1) * t342 - g(2) * t336;
t343 = qJD(1) ^ 2;
t305 = -pkin(1) * t343 + qJDD(1) * pkin(7) + t317;
t365 = t335 * t305;
t366 = pkin(2) * t343;
t272 = qJDD(2) * pkin(2) - t310 * pkin(8) - t365 + (pkin(8) * t362 + t335 * t366 - g(3)) * t341;
t293 = -g(3) * t335 + t341 * t305;
t311 = qJDD(1) * t341 - t335 * t362;
t364 = qJD(1) * t335;
t315 = qJD(2) * pkin(2) - pkin(8) * t364;
t330 = t341 ^ 2;
t273 = pkin(8) * t311 - qJD(2) * t315 - t330 * t366 + t293;
t334 = sin(qJ(3));
t340 = cos(qJ(3));
t252 = t340 * t272 - t334 * t273;
t302 = (-t334 * t335 + t340 * t341) * qJD(1);
t278 = qJD(3) * t302 + t310 * t340 + t311 * t334;
t303 = (t334 * t341 + t335 * t340) * qJD(1);
t327 = qJDD(2) + qJDD(3);
t328 = qJD(2) + qJD(3);
t229 = (t302 * t328 - t278) * pkin(9) + (t302 * t303 + t327) * pkin(3) + t252;
t253 = t334 * t272 + t340 * t273;
t277 = -qJD(3) * t303 - t310 * t334 + t311 * t340;
t296 = pkin(3) * t328 - pkin(9) * t303;
t298 = t302 ^ 2;
t231 = -pkin(3) * t298 + pkin(9) * t277 - t296 * t328 + t253;
t333 = sin(qJ(4));
t339 = cos(qJ(4));
t211 = t339 * t229 - t333 * t231;
t289 = t302 * t339 - t303 * t333;
t249 = qJD(4) * t289 + t277 * t333 + t278 * t339;
t290 = t302 * t333 + t303 * t339;
t324 = qJDD(4) + t327;
t325 = qJD(4) + t328;
t204 = (t289 * t325 - t249) * pkin(10) + (t289 * t290 + t324) * pkin(4) + t211;
t212 = t333 * t229 + t339 * t231;
t248 = -qJD(4) * t290 + t277 * t339 - t278 * t333;
t282 = pkin(4) * t325 - pkin(10) * t290;
t288 = t289 ^ 2;
t206 = -pkin(4) * t288 + pkin(10) * t248 - t282 * t325 + t212;
t332 = sin(qJ(5));
t338 = cos(qJ(5));
t201 = t332 * t204 + t338 * t206;
t266 = t289 * t332 + t290 * t338;
t220 = -qJD(5) * t266 + t248 * t338 - t249 * t332;
t265 = t289 * t338 - t290 * t332;
t240 = -mrSges(6,1) * t265 + mrSges(6,2) * t266;
t322 = qJD(5) + t325;
t257 = mrSges(6,1) * t322 - mrSges(6,3) * t266;
t320 = qJDD(5) + t324;
t241 = -pkin(5) * t265 - pkin(11) * t266;
t321 = t322 ^ 2;
t198 = -pkin(5) * t321 + pkin(11) * t320 + t241 * t265 + t201;
t316 = t336 * g(1) - t342 * g(2);
t355 = -qJDD(1) * pkin(1) - t316;
t279 = -t311 * pkin(2) + t315 * t364 + (-pkin(8) * t330 - pkin(7)) * t343 + t355;
t243 = -t277 * pkin(3) - t298 * pkin(9) + t303 * t296 + t279;
t214 = -t248 * pkin(4) - t288 * pkin(10) + t290 * t282 + t243;
t221 = qJD(5) * t265 + t248 * t332 + t249 * t338;
t202 = t214 + (-t265 * t322 - t221) * pkin(11) + (t266 * t322 - t220) * pkin(5);
t331 = sin(qJ(6));
t337 = cos(qJ(6));
t195 = -t198 * t331 + t202 * t337;
t254 = -t266 * t331 + t322 * t337;
t209 = qJD(6) * t254 + t221 * t337 + t320 * t331;
t219 = qJDD(6) - t220;
t255 = t266 * t337 + t322 * t331;
t232 = -mrSges(7,1) * t254 + mrSges(7,2) * t255;
t264 = qJD(6) - t265;
t233 = -mrSges(7,2) * t264 + mrSges(7,3) * t254;
t191 = m(7) * t195 + mrSges(7,1) * t219 - mrSges(7,3) * t209 - t232 * t255 + t233 * t264;
t196 = t198 * t337 + t202 * t331;
t208 = -qJD(6) * t255 - t221 * t331 + t320 * t337;
t234 = mrSges(7,1) * t264 - mrSges(7,3) * t255;
t192 = m(7) * t196 - mrSges(7,2) * t219 + mrSges(7,3) * t208 + t232 * t254 - t234 * t264;
t357 = -t191 * t331 + t337 * t192;
t178 = m(6) * t201 - mrSges(6,2) * t320 + mrSges(6,3) * t220 + t240 * t265 - t257 * t322 + t357;
t200 = t204 * t338 - t206 * t332;
t256 = -mrSges(6,2) * t322 + mrSges(6,3) * t265;
t197 = -pkin(5) * t320 - pkin(11) * t321 + t241 * t266 - t200;
t352 = -m(7) * t197 + t208 * mrSges(7,1) - mrSges(7,2) * t209 + t254 * t233 - t234 * t255;
t187 = m(6) * t200 + mrSges(6,1) * t320 - mrSges(6,3) * t221 - t240 * t266 + t256 * t322 + t352;
t173 = t332 * t178 + t338 * t187;
t267 = -mrSges(5,1) * t289 + mrSges(5,2) * t290;
t280 = -mrSges(5,2) * t325 + mrSges(5,3) * t289;
t170 = m(5) * t211 + mrSges(5,1) * t324 - mrSges(5,3) * t249 - t267 * t290 + t280 * t325 + t173;
t281 = mrSges(5,1) * t325 - mrSges(5,3) * t290;
t358 = t338 * t178 - t187 * t332;
t171 = m(5) * t212 - mrSges(5,2) * t324 + mrSges(5,3) * t248 + t267 * t289 - t281 * t325 + t358;
t164 = t339 * t170 + t333 * t171;
t291 = -mrSges(4,1) * t302 + mrSges(4,2) * t303;
t294 = -mrSges(4,2) * t328 + mrSges(4,3) * t302;
t161 = m(4) * t252 + mrSges(4,1) * t327 - mrSges(4,3) * t278 - t291 * t303 + t294 * t328 + t164;
t295 = mrSges(4,1) * t328 - mrSges(4,3) * t303;
t359 = -t170 * t333 + t339 * t171;
t162 = m(4) * t253 - mrSges(4,2) * t327 + mrSges(4,3) * t277 + t291 * t302 - t295 * t328 + t359;
t156 = t340 * t161 + t334 * t162;
t292 = -t341 * g(3) - t365;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t335 + Ifges(3,2) * t341) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t335 + Ifges(3,4) * t341) * qJD(1);
t284 = Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * t328;
t285 = Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * t328;
t259 = Ifges(5,4) * t290 + Ifges(5,2) * t289 + Ifges(5,6) * t325;
t260 = Ifges(5,1) * t290 + Ifges(5,4) * t289 + Ifges(5,5) * t325;
t222 = Ifges(7,5) * t255 + Ifges(7,6) * t254 + Ifges(7,3) * t264;
t224 = Ifges(7,1) * t255 + Ifges(7,4) * t254 + Ifges(7,5) * t264;
t184 = -mrSges(7,1) * t197 + mrSges(7,3) * t196 + Ifges(7,4) * t209 + Ifges(7,2) * t208 + Ifges(7,6) * t219 - t222 * t255 + t224 * t264;
t223 = Ifges(7,4) * t255 + Ifges(7,2) * t254 + Ifges(7,6) * t264;
t185 = mrSges(7,2) * t197 - mrSges(7,3) * t195 + Ifges(7,1) * t209 + Ifges(7,4) * t208 + Ifges(7,5) * t219 + t222 * t254 - t223 * t264;
t236 = Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t322;
t237 = Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t322;
t350 = -mrSges(6,1) * t200 + mrSges(6,2) * t201 - Ifges(6,5) * t221 - Ifges(6,6) * t220 - Ifges(6,3) * t320 - pkin(5) * t352 - pkin(11) * t357 - t337 * t184 - t331 * t185 - t266 * t236 + t265 * t237;
t347 = -mrSges(5,1) * t211 + mrSges(5,2) * t212 - Ifges(5,5) * t249 - Ifges(5,6) * t248 - Ifges(5,3) * t324 - pkin(4) * t173 - t290 * t259 + t289 * t260 + t350;
t345 = mrSges(4,1) * t252 - mrSges(4,2) * t253 + Ifges(4,5) * t278 + Ifges(4,6) * t277 + Ifges(4,3) * t327 + pkin(3) * t164 + t303 * t284 - t302 * t285 - t347;
t367 = mrSges(3,1) * t292 - mrSges(3,2) * t293 + Ifges(3,5) * t310 + Ifges(3,6) * t311 + Ifges(3,3) * qJDD(2) + pkin(2) * t156 + (t300 * t335 - t301 * t341) * qJD(1) + t345;
t180 = t337 * t191 + t331 * t192;
t363 = qJD(1) * t341;
t309 = (-mrSges(3,1) * t341 + mrSges(3,2) * t335) * qJD(1);
t314 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t154 = m(3) * t292 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t310 + qJD(2) * t314 - t309 * t364 + t156;
t313 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t360 = -t161 * t334 + t340 * t162;
t155 = m(3) * t293 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t311 - qJD(2) * t313 + t309 * t363 + t360;
t361 = -t154 * t335 + t341 * t155;
t354 = m(6) * t214 - t220 * mrSges(6,1) + t221 * mrSges(6,2) - t265 * t256 + t266 * t257 + t180;
t235 = Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t322;
t165 = mrSges(6,2) * t214 - mrSges(6,3) * t200 + Ifges(6,1) * t221 + Ifges(6,4) * t220 + Ifges(6,5) * t320 - pkin(11) * t180 - t184 * t331 + t185 * t337 + t235 * t265 - t236 * t322;
t349 = mrSges(7,1) * t195 - mrSges(7,2) * t196 + Ifges(7,5) * t209 + Ifges(7,6) * t208 + Ifges(7,3) * t219 + t223 * t255 - t224 * t254;
t166 = -mrSges(6,1) * t214 + mrSges(6,3) * t201 + Ifges(6,4) * t221 + Ifges(6,2) * t220 + Ifges(6,6) * t320 - pkin(5) * t180 - t235 * t266 + t237 * t322 - t349;
t258 = Ifges(5,5) * t290 + Ifges(5,6) * t289 + Ifges(5,3) * t325;
t152 = -mrSges(5,1) * t243 + mrSges(5,3) * t212 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t324 - pkin(4) * t354 + pkin(10) * t358 + t332 * t165 + t338 * t166 - t290 * t258 + t325 * t260;
t157 = mrSges(5,2) * t243 - mrSges(5,3) * t211 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t324 - pkin(10) * t173 + t165 * t338 - t166 * t332 + t258 * t289 - t259 * t325;
t283 = Ifges(4,5) * t303 + Ifges(4,6) * t302 + Ifges(4,3) * t328;
t351 = m(5) * t243 - t248 * mrSges(5,1) + t249 * mrSges(5,2) - t289 * t280 + t290 * t281 + t354;
t147 = -mrSges(4,1) * t279 + mrSges(4,3) * t253 + Ifges(4,4) * t278 + Ifges(4,2) * t277 + Ifges(4,6) * t327 - pkin(3) * t351 + pkin(9) * t359 + t339 * t152 + t333 * t157 - t303 * t283 + t328 * t285;
t148 = mrSges(4,2) * t279 - mrSges(4,3) * t252 + Ifges(4,1) * t278 + Ifges(4,4) * t277 + Ifges(4,5) * t327 - pkin(9) * t164 - t152 * t333 + t157 * t339 + t283 * t302 - t284 * t328;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t335 + Ifges(3,6) * t341) * qJD(1);
t304 = -t343 * pkin(7) + t355;
t348 = m(4) * t279 - t277 * mrSges(4,1) + t278 * mrSges(4,2) - t302 * t294 + t303 * t295 + t351;
t143 = -mrSges(3,1) * t304 + mrSges(3,3) * t293 + Ifges(3,4) * t310 + Ifges(3,2) * t311 + Ifges(3,6) * qJDD(2) - pkin(2) * t348 + pkin(8) * t360 + qJD(2) * t301 + t340 * t147 + t334 * t148 - t299 * t364;
t145 = mrSges(3,2) * t304 - mrSges(3,3) * t292 + Ifges(3,1) * t310 + Ifges(3,4) * t311 + Ifges(3,5) * qJDD(2) - pkin(8) * t156 - qJD(2) * t300 - t147 * t334 + t148 * t340 + t299 * t363;
t346 = -m(3) * t304 + t311 * mrSges(3,1) - t310 * mrSges(3,2) - t313 * t364 + t314 * t363 - t348;
t353 = mrSges(2,1) * t316 - mrSges(2,2) * t317 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(7) * t361 + t341 * t143 + t335 * t145;
t174 = m(2) * t316 + qJDD(1) * mrSges(2,1) - t343 * mrSges(2,2) + t346;
t151 = t154 * t341 + t155 * t335;
t149 = m(2) * t317 - mrSges(2,1) * t343 - qJDD(1) * mrSges(2,2) + t361;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t317 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t367;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t316 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t343 - pkin(7) * t151 - t143 * t335 + t145 * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t141 - t336 * t146 - pkin(6) * (t149 * t336 + t174 * t342), t141, t145, t148, t157, t165, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t336 * t141 + t342 * t146 + pkin(6) * (t149 * t342 - t174 * t336), t146, t143, t147, t152, t166, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, t367, t345, -t347, -t350, t349;];
m_new  = t1;
