% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR1
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
% Datum: 2019-05-07 09:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:40:14
% EndTime: 2019-05-07 09:42:11
% DurationCPUTime: 53.52s
% Computational Cost: add. (901975->385), mult. (2110800->488), div. (0->0), fcn. (1608609->12), ass. (0->153)
t336 = sin(qJ(2));
t341 = cos(qJ(2));
t362 = qJD(1) * qJD(2);
t312 = qJDD(1) * t336 + t341 * t362;
t337 = sin(qJ(1));
t342 = cos(qJ(1));
t319 = -g(1) * t342 - g(2) * t337;
t343 = qJD(1) ^ 2;
t307 = -pkin(1) * t343 + qJDD(1) * pkin(7) + t319;
t365 = t307 * t336;
t366 = pkin(2) * t343;
t272 = qJDD(2) * pkin(2) - pkin(8) * t312 - t365 + (pkin(8) * t362 + t336 * t366 - g(3)) * t341;
t295 = -g(3) * t336 + t341 * t307;
t313 = qJDD(1) * t341 - t336 * t362;
t364 = qJD(1) * t336;
t317 = qJD(2) * pkin(2) - pkin(8) * t364;
t330 = t341 ^ 2;
t273 = pkin(8) * t313 - qJD(2) * t317 - t330 * t366 + t295;
t335 = sin(qJ(3));
t340 = cos(qJ(3));
t246 = t340 * t272 - t273 * t335;
t304 = (-t335 * t336 + t340 * t341) * qJD(1);
t278 = qJD(3) * t304 + t312 * t340 + t313 * t335;
t305 = (t335 * t341 + t336 * t340) * qJD(1);
t327 = qJDD(2) + qJDD(3);
t328 = qJD(2) + qJD(3);
t229 = (t304 * t328 - t278) * qJ(4) + (t304 * t305 + t327) * pkin(3) + t246;
t247 = t335 * t272 + t340 * t273;
t277 = -qJD(3) * t305 - t312 * t335 + t313 * t340;
t297 = pkin(3) * t328 - qJ(4) * t305;
t300 = t304 ^ 2;
t231 = -pkin(3) * t300 + qJ(4) * t277 - t297 * t328 + t247;
t331 = sin(pkin(11));
t332 = cos(pkin(11));
t292 = t304 * t331 + t305 * t332;
t208 = -0.2e1 * qJD(4) * t292 + t332 * t229 - t231 * t331;
t253 = t277 * t331 + t278 * t332;
t291 = t304 * t332 - t305 * t331;
t204 = (t291 * t328 - t253) * pkin(9) + (t291 * t292 + t327) * pkin(4) + t208;
t209 = 0.2e1 * qJD(4) * t291 + t331 * t229 + t332 * t231;
t252 = t277 * t332 - t278 * t331;
t282 = pkin(4) * t328 - pkin(9) * t292;
t288 = t291 ^ 2;
t206 = -pkin(4) * t288 + pkin(9) * t252 - t282 * t328 + t209;
t334 = sin(qJ(5));
t339 = cos(qJ(5));
t201 = t334 * t204 + t339 * t206;
t266 = t291 * t334 + t292 * t339;
t220 = -qJD(5) * t266 + t252 * t339 - t253 * t334;
t265 = t291 * t339 - t292 * t334;
t240 = -mrSges(6,1) * t265 + mrSges(6,2) * t266;
t325 = qJD(5) + t328;
t257 = mrSges(6,1) * t325 - mrSges(6,3) * t266;
t324 = qJDD(5) + t327;
t241 = -pkin(5) * t265 - pkin(10) * t266;
t323 = t325 ^ 2;
t198 = -pkin(5) * t323 + pkin(10) * t324 + t241 * t265 + t201;
t318 = g(1) * t337 - t342 * g(2);
t355 = -qJDD(1) * pkin(1) - t318;
t279 = -pkin(2) * t313 + t317 * t364 + (-pkin(8) * t330 - pkin(7)) * t343 + t355;
t243 = -pkin(3) * t277 - qJ(4) * t300 + t305 * t297 + qJDD(4) + t279;
t214 = -pkin(4) * t252 - pkin(9) * t288 + t292 * t282 + t243;
t221 = qJD(5) * t265 + t252 * t334 + t253 * t339;
t202 = t214 + (-t265 * t325 - t221) * pkin(10) + (t266 * t325 - t220) * pkin(5);
t333 = sin(qJ(6));
t338 = cos(qJ(6));
t195 = -t198 * t333 + t202 * t338;
t254 = -t266 * t333 + t325 * t338;
t212 = qJD(6) * t254 + t221 * t338 + t324 * t333;
t219 = qJDD(6) - t220;
t255 = t266 * t338 + t325 * t333;
t232 = -mrSges(7,1) * t254 + mrSges(7,2) * t255;
t261 = qJD(6) - t265;
t233 = -mrSges(7,2) * t261 + mrSges(7,3) * t254;
t191 = m(7) * t195 + mrSges(7,1) * t219 - mrSges(7,3) * t212 - t232 * t255 + t233 * t261;
t196 = t198 * t338 + t202 * t333;
t211 = -qJD(6) * t255 - t221 * t333 + t324 * t338;
t234 = mrSges(7,1) * t261 - mrSges(7,3) * t255;
t192 = m(7) * t196 - mrSges(7,2) * t219 + mrSges(7,3) * t211 + t232 * t254 - t234 * t261;
t357 = -t191 * t333 + t338 * t192;
t178 = m(6) * t201 - mrSges(6,2) * t324 + mrSges(6,3) * t220 + t240 * t265 - t257 * t325 + t357;
t200 = t204 * t339 - t206 * t334;
t256 = -mrSges(6,2) * t325 + mrSges(6,3) * t265;
t197 = -pkin(5) * t324 - pkin(10) * t323 + t241 * t266 - t200;
t352 = -m(7) * t197 + t211 * mrSges(7,1) - mrSges(7,2) * t212 + t254 * t233 - t234 * t255;
t187 = m(6) * t200 + mrSges(6,1) * t324 - mrSges(6,3) * t221 - t240 * t266 + t256 * t325 + t352;
t173 = t334 * t178 + t339 * t187;
t267 = -mrSges(5,1) * t291 + mrSges(5,2) * t292;
t280 = -mrSges(5,2) * t328 + mrSges(5,3) * t291;
t170 = m(5) * t208 + mrSges(5,1) * t327 - mrSges(5,3) * t253 - t267 * t292 + t280 * t328 + t173;
t281 = mrSges(5,1) * t328 - mrSges(5,3) * t292;
t358 = t339 * t178 - t187 * t334;
t171 = m(5) * t209 - mrSges(5,2) * t327 + mrSges(5,3) * t252 + t267 * t291 - t281 * t328 + t358;
t164 = t332 * t170 + t331 * t171;
t293 = -mrSges(4,1) * t304 + mrSges(4,2) * t305;
t296 = -mrSges(4,2) * t328 + mrSges(4,3) * t304;
t161 = m(4) * t246 + mrSges(4,1) * t327 - mrSges(4,3) * t278 - t293 * t305 + t296 * t328 + t164;
t298 = mrSges(4,1) * t328 - mrSges(4,3) * t305;
t359 = -t170 * t331 + t332 * t171;
t162 = m(4) * t247 - mrSges(4,2) * t327 + mrSges(4,3) * t277 + t293 * t304 - t298 * t328 + t359;
t156 = t340 * t161 + t335 * t162;
t294 = -g(3) * t341 - t365;
t302 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t336 + Ifges(3,2) * t341) * qJD(1);
t303 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t336 + Ifges(3,4) * t341) * qJD(1);
t286 = Ifges(4,4) * t305 + Ifges(4,2) * t304 + Ifges(4,6) * t328;
t287 = Ifges(4,1) * t305 + Ifges(4,4) * t304 + Ifges(4,5) * t328;
t259 = Ifges(5,4) * t292 + Ifges(5,2) * t291 + Ifges(5,6) * t328;
t260 = Ifges(5,1) * t292 + Ifges(5,4) * t291 + Ifges(5,5) * t328;
t222 = Ifges(7,5) * t255 + Ifges(7,6) * t254 + Ifges(7,3) * t261;
t224 = Ifges(7,1) * t255 + Ifges(7,4) * t254 + Ifges(7,5) * t261;
t184 = -mrSges(7,1) * t197 + mrSges(7,3) * t196 + Ifges(7,4) * t212 + Ifges(7,2) * t211 + Ifges(7,6) * t219 - t222 * t255 + t224 * t261;
t223 = Ifges(7,4) * t255 + Ifges(7,2) * t254 + Ifges(7,6) * t261;
t185 = mrSges(7,2) * t197 - mrSges(7,3) * t195 + Ifges(7,1) * t212 + Ifges(7,4) * t211 + Ifges(7,5) * t219 + t222 * t254 - t223 * t261;
t236 = Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t325;
t237 = Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t325;
t350 = -mrSges(6,1) * t200 + mrSges(6,2) * t201 - Ifges(6,5) * t221 - Ifges(6,6) * t220 - Ifges(6,3) * t324 - pkin(5) * t352 - pkin(10) * t357 - t338 * t184 - t333 * t185 - t266 * t236 + t265 * t237;
t347 = -mrSges(5,1) * t208 + mrSges(5,2) * t209 - Ifges(5,5) * t253 - Ifges(5,6) * t252 - Ifges(5,3) * t327 - pkin(4) * t173 - t292 * t259 + t291 * t260 + t350;
t345 = mrSges(4,1) * t246 - mrSges(4,2) * t247 + Ifges(4,5) * t278 + Ifges(4,6) * t277 + Ifges(4,3) * t327 + pkin(3) * t164 + t305 * t286 - t304 * t287 - t347;
t367 = mrSges(3,1) * t294 - mrSges(3,2) * t295 + Ifges(3,5) * t312 + Ifges(3,6) * t313 + Ifges(3,3) * qJDD(2) + pkin(2) * t156 + (t336 * t302 - t341 * t303) * qJD(1) + t345;
t180 = t338 * t191 + t333 * t192;
t363 = qJD(1) * t341;
t311 = (-mrSges(3,1) * t341 + mrSges(3,2) * t336) * qJD(1);
t316 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t154 = m(3) * t294 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t312 + qJD(2) * t316 - t311 * t364 + t156;
t315 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t360 = -t335 * t161 + t340 * t162;
t155 = m(3) * t295 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t313 - qJD(2) * t315 + t311 * t363 + t360;
t361 = -t154 * t336 + t341 * t155;
t354 = m(6) * t214 - t220 * mrSges(6,1) + t221 * mrSges(6,2) - t265 * t256 + t266 * t257 + t180;
t235 = Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t325;
t165 = mrSges(6,2) * t214 - mrSges(6,3) * t200 + Ifges(6,1) * t221 + Ifges(6,4) * t220 + Ifges(6,5) * t324 - pkin(10) * t180 - t184 * t333 + t185 * t338 + t235 * t265 - t236 * t325;
t349 = mrSges(7,1) * t195 - mrSges(7,2) * t196 + Ifges(7,5) * t212 + Ifges(7,6) * t211 + Ifges(7,3) * t219 + t223 * t255 - t224 * t254;
t166 = -mrSges(6,1) * t214 + mrSges(6,3) * t201 + Ifges(6,4) * t221 + Ifges(6,2) * t220 + Ifges(6,6) * t324 - pkin(5) * t180 - t235 * t266 + t237 * t325 - t349;
t258 = Ifges(5,5) * t292 + Ifges(5,6) * t291 + Ifges(5,3) * t328;
t152 = -mrSges(5,1) * t243 + mrSges(5,3) * t209 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * t327 - pkin(4) * t354 + pkin(9) * t358 + t334 * t165 + t339 * t166 - t292 * t258 + t328 * t260;
t157 = mrSges(5,2) * t243 - mrSges(5,3) * t208 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * t327 - pkin(9) * t173 + t165 * t339 - t166 * t334 + t258 * t291 - t259 * t328;
t285 = Ifges(4,5) * t305 + Ifges(4,6) * t304 + Ifges(4,3) * t328;
t351 = m(5) * t243 - t252 * mrSges(5,1) + t253 * mrSges(5,2) - t291 * t280 + t292 * t281 + t354;
t147 = -mrSges(4,1) * t279 + mrSges(4,3) * t247 + Ifges(4,4) * t278 + Ifges(4,2) * t277 + Ifges(4,6) * t327 - pkin(3) * t351 + qJ(4) * t359 + t332 * t152 + t331 * t157 - t305 * t285 + t328 * t287;
t148 = mrSges(4,2) * t279 - mrSges(4,3) * t246 + Ifges(4,1) * t278 + Ifges(4,4) * t277 + Ifges(4,5) * t327 - qJ(4) * t164 - t152 * t331 + t157 * t332 + t285 * t304 - t286 * t328;
t301 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t336 + Ifges(3,6) * t341) * qJD(1);
t306 = -pkin(7) * t343 + t355;
t348 = m(4) * t279 - t277 * mrSges(4,1) + t278 * mrSges(4,2) - t304 * t296 + t305 * t298 + t351;
t143 = -mrSges(3,1) * t306 + mrSges(3,3) * t295 + Ifges(3,4) * t312 + Ifges(3,2) * t313 + Ifges(3,6) * qJDD(2) - pkin(2) * t348 + pkin(8) * t360 + qJD(2) * t303 + t340 * t147 + t335 * t148 - t301 * t364;
t145 = mrSges(3,2) * t306 - mrSges(3,3) * t294 + Ifges(3,1) * t312 + Ifges(3,4) * t313 + Ifges(3,5) * qJDD(2) - pkin(8) * t156 - qJD(2) * t302 - t147 * t335 + t148 * t340 + t301 * t363;
t346 = -m(3) * t306 + t313 * mrSges(3,1) - t312 * mrSges(3,2) - t315 * t364 + t316 * t363 - t348;
t353 = mrSges(2,1) * t318 - mrSges(2,2) * t319 + Ifges(2,3) * qJDD(1) + pkin(1) * t346 + pkin(7) * t361 + t341 * t143 + t336 * t145;
t174 = m(2) * t318 + qJDD(1) * mrSges(2,1) - t343 * mrSges(2,2) + t346;
t151 = t154 * t341 + t155 * t336;
t149 = m(2) * t319 - mrSges(2,1) * t343 - qJDD(1) * mrSges(2,2) + t361;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t319 + t343 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t367;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t318 + Ifges(2,5) * qJDD(1) - t343 * Ifges(2,6) - pkin(7) * t151 - t143 * t336 + t341 * t145;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t342 * t141 - t337 * t146 - pkin(6) * (t337 * t149 + t342 * t174), t141, t145, t148, t157, t165, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t337 * t141 + t342 * t146 + pkin(6) * (t149 * t342 - t337 * t174), t146, t143, t147, t152, t166, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, t367, t345, -t347, -t350, t349;];
m_new  = t1;
