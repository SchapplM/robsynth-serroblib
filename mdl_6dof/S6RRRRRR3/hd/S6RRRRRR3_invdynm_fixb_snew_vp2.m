% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR3
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
% Datum: 2019-05-08 09:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:51:55
% EndTime: 2019-05-08 08:54:01
% DurationCPUTime: 51.96s
% Computational Cost: add. (987544->387), mult. (1966283->485), div. (0->0), fcn. (1467468->12), ass. (0->155)
t338 = sin(qJ(2));
t344 = cos(qJ(2));
t365 = qJD(1) * qJD(2);
t318 = t338 * qJDD(1) + t344 * t365;
t339 = sin(qJ(1));
t345 = cos(qJ(1));
t325 = -t345 * g(1) - t339 * g(2);
t346 = qJD(1) ^ 2;
t312 = -t346 * pkin(1) + qJDD(1) * pkin(7) + t325;
t368 = t338 * t312;
t369 = pkin(2) * t346;
t271 = qJDD(2) * pkin(2) - t318 * pkin(8) - t368 + (pkin(8) * t365 + t338 * t369 - g(3)) * t344;
t298 = -t338 * g(3) + t344 * t312;
t319 = t344 * qJDD(1) - t338 * t365;
t367 = qJD(1) * t338;
t323 = qJD(2) * pkin(2) - pkin(8) * t367;
t333 = t344 ^ 2;
t273 = t319 * pkin(8) - qJD(2) * t323 - t333 * t369 + t298;
t337 = sin(qJ(3));
t343 = cos(qJ(3));
t251 = t337 * t271 + t343 * t273;
t310 = (t337 * t344 + t338 * t343) * qJD(1);
t282 = -t310 * qJD(3) - t337 * t318 + t343 * t319;
t366 = qJD(1) * t344;
t309 = -t337 * t367 + t343 * t366;
t292 = -t309 * mrSges(4,1) + t310 * mrSges(4,2);
t331 = qJD(2) + qJD(3);
t300 = t331 * mrSges(4,1) - t310 * mrSges(4,3);
t330 = qJDD(2) + qJDD(3);
t283 = t309 * qJD(3) + t343 * t318 + t337 * t319;
t324 = t339 * g(1) - t345 * g(2);
t357 = -qJDD(1) * pkin(1) - t324;
t284 = -t319 * pkin(2) + t323 * t367 + (-pkin(8) * t333 - pkin(7)) * t346 + t357;
t236 = (-t309 * t331 - t283) * pkin(9) + (t310 * t331 - t282) * pkin(3) + t284;
t293 = -t309 * pkin(3) - t310 * pkin(9);
t329 = t331 ^ 2;
t239 = -t329 * pkin(3) + t330 * pkin(9) + t309 * t293 + t251;
t336 = sin(qJ(4));
t342 = cos(qJ(4));
t219 = t342 * t236 - t336 * t239;
t295 = -t336 * t310 + t342 * t331;
t254 = t295 * qJD(4) + t342 * t283 + t336 * t330;
t281 = qJDD(4) - t282;
t296 = t342 * t310 + t336 * t331;
t305 = qJD(4) - t309;
t215 = (t295 * t305 - t254) * pkin(10) + (t295 * t296 + t281) * pkin(4) + t219;
t220 = t336 * t236 + t342 * t239;
t253 = -t296 * qJD(4) - t336 * t283 + t342 * t330;
t287 = t305 * pkin(4) - t296 * pkin(10);
t294 = t295 ^ 2;
t217 = -t294 * pkin(4) + t253 * pkin(10) - t305 * t287 + t220;
t335 = sin(qJ(5));
t341 = cos(qJ(5));
t203 = t341 * t215 - t335 * t217;
t264 = t341 * t295 - t335 * t296;
t232 = t264 * qJD(5) + t335 * t253 + t341 * t254;
t265 = t335 * t295 + t341 * t296;
t276 = qJDD(5) + t281;
t303 = qJD(5) + t305;
t200 = (t264 * t303 - t232) * pkin(11) + (t264 * t265 + t276) * pkin(5) + t203;
t204 = t335 * t215 + t341 * t217;
t231 = -t265 * qJD(5) + t341 * t253 - t335 * t254;
t257 = t303 * pkin(5) - t265 * pkin(11);
t263 = t264 ^ 2;
t201 = -t263 * pkin(5) + t231 * pkin(11) - t303 * t257 + t204;
t334 = sin(qJ(6));
t340 = cos(qJ(6));
t198 = t340 * t200 - t334 * t201;
t246 = t340 * t264 - t334 * t265;
t212 = t246 * qJD(6) + t334 * t231 + t340 * t232;
t247 = t334 * t264 + t340 * t265;
t227 = -t246 * mrSges(7,1) + t247 * mrSges(7,2);
t301 = qJD(6) + t303;
t240 = -t301 * mrSges(7,2) + t246 * mrSges(7,3);
t275 = qJDD(6) + t276;
t193 = m(7) * t198 + t275 * mrSges(7,1) - t212 * mrSges(7,3) - t247 * t227 + t301 * t240;
t199 = t334 * t200 + t340 * t201;
t211 = -t247 * qJD(6) + t340 * t231 - t334 * t232;
t241 = t301 * mrSges(7,1) - t247 * mrSges(7,3);
t194 = m(7) * t199 - t275 * mrSges(7,2) + t211 * mrSges(7,3) + t246 * t227 - t301 * t241;
t185 = t340 * t193 + t334 * t194;
t248 = -t264 * mrSges(6,1) + t265 * mrSges(6,2);
t255 = -t303 * mrSges(6,2) + t264 * mrSges(6,3);
t182 = m(6) * t203 + t276 * mrSges(6,1) - t232 * mrSges(6,3) - t265 * t248 + t303 * t255 + t185;
t256 = t303 * mrSges(6,1) - t265 * mrSges(6,3);
t360 = -t334 * t193 + t340 * t194;
t183 = m(6) * t204 - t276 * mrSges(6,2) + t231 * mrSges(6,3) + t264 * t248 - t303 * t256 + t360;
t178 = t341 * t182 + t335 * t183;
t269 = -t295 * mrSges(5,1) + t296 * mrSges(5,2);
t285 = -t305 * mrSges(5,2) + t295 * mrSges(5,3);
t176 = m(5) * t219 + t281 * mrSges(5,1) - t254 * mrSges(5,3) - t296 * t269 + t305 * t285 + t178;
t286 = t305 * mrSges(5,1) - t296 * mrSges(5,3);
t361 = -t335 * t182 + t341 * t183;
t177 = m(5) * t220 - t281 * mrSges(5,2) + t253 * mrSges(5,3) + t295 * t269 - t305 * t286 + t361;
t362 = -t336 * t176 + t342 * t177;
t167 = m(4) * t251 - t330 * mrSges(4,2) + t282 * mrSges(4,3) + t309 * t292 - t331 * t300 + t362;
t250 = t343 * t271 - t337 * t273;
t299 = -t331 * mrSges(4,2) + t309 * mrSges(4,3);
t238 = -t330 * pkin(3) - t329 * pkin(9) + t310 * t293 - t250;
t221 = -t253 * pkin(4) - t294 * pkin(10) + t296 * t287 + t238;
t206 = -t231 * pkin(5) - t263 * pkin(11) + t265 * t257 + t221;
t359 = m(7) * t206 - t211 * mrSges(7,1) + t212 * mrSges(7,2) - t246 * t240 + t247 * t241;
t352 = m(6) * t221 - t231 * mrSges(6,1) + t232 * mrSges(6,2) - t264 * t255 + t265 * t256 + t359;
t349 = -m(5) * t238 + t253 * mrSges(5,1) - t254 * mrSges(5,2) + t295 * t285 - t296 * t286 - t352;
t189 = m(4) * t250 + t330 * mrSges(4,1) - t283 * mrSges(4,3) - t310 * t292 + t331 * t299 + t349;
t162 = t337 * t167 + t343 * t189;
t297 = -t344 * g(3) - t368;
t307 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t338 + Ifges(3,2) * t344) * qJD(1);
t308 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t338 + Ifges(3,4) * t344) * qJD(1);
t222 = Ifges(7,5) * t247 + Ifges(7,6) * t246 + Ifges(7,3) * t301;
t224 = Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t301;
t186 = -mrSges(7,1) * t206 + mrSges(7,3) * t199 + Ifges(7,4) * t212 + Ifges(7,2) * t211 + Ifges(7,6) * t275 - t247 * t222 + t301 * t224;
t223 = Ifges(7,4) * t247 + Ifges(7,2) * t246 + Ifges(7,6) * t301;
t187 = mrSges(7,2) * t206 - mrSges(7,3) * t198 + Ifges(7,1) * t212 + Ifges(7,4) * t211 + Ifges(7,5) * t275 + t246 * t222 - t301 * t223;
t242 = Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * t303;
t244 = Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t303;
t171 = -mrSges(6,1) * t221 + mrSges(6,3) * t204 + Ifges(6,4) * t232 + Ifges(6,2) * t231 + Ifges(6,6) * t276 - pkin(5) * t359 + pkin(11) * t360 + t340 * t186 + t334 * t187 - t265 * t242 + t303 * t244;
t243 = Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t303;
t172 = mrSges(6,2) * t221 - mrSges(6,3) * t203 + Ifges(6,1) * t232 + Ifges(6,4) * t231 + Ifges(6,5) * t276 - pkin(11) * t185 - t334 * t186 + t340 * t187 + t264 * t242 - t303 * t243;
t258 = Ifges(5,5) * t296 + Ifges(5,6) * t295 + Ifges(5,3) * t305;
t260 = Ifges(5,1) * t296 + Ifges(5,4) * t295 + Ifges(5,5) * t305;
t156 = -mrSges(5,1) * t238 + mrSges(5,3) * t220 + Ifges(5,4) * t254 + Ifges(5,2) * t253 + Ifges(5,6) * t281 - pkin(4) * t352 + pkin(10) * t361 + t341 * t171 + t335 * t172 - t296 * t258 + t305 * t260;
t259 = Ifges(5,4) * t296 + Ifges(5,2) * t295 + Ifges(5,6) * t305;
t158 = mrSges(5,2) * t238 - mrSges(5,3) * t219 + Ifges(5,1) * t254 + Ifges(5,4) * t253 + Ifges(5,5) * t281 - pkin(10) * t178 - t335 * t171 + t341 * t172 + t295 * t258 - t305 * t259;
t289 = Ifges(4,4) * t310 + Ifges(4,2) * t309 + Ifges(4,6) * t331;
t290 = Ifges(4,1) * t310 + Ifges(4,4) * t309 + Ifges(4,5) * t331;
t353 = -mrSges(4,1) * t250 + mrSges(4,2) * t251 - Ifges(4,5) * t283 - Ifges(4,6) * t282 - Ifges(4,3) * t330 - pkin(3) * t349 - pkin(9) * t362 - t342 * t156 - t336 * t158 - t310 * t289 + t309 * t290;
t370 = mrSges(3,1) * t297 - mrSges(3,2) * t298 + Ifges(3,5) * t318 + Ifges(3,6) * t319 + Ifges(3,3) * qJDD(2) + pkin(2) * t162 + (t338 * t307 - t344 * t308) * qJD(1) - t353;
t169 = t342 * t176 + t336 * t177;
t317 = (-mrSges(3,1) * t344 + mrSges(3,2) * t338) * qJD(1);
t322 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t366;
t160 = m(3) * t297 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,3) + qJD(2) * t322 - t317 * t367 + t162;
t321 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t367;
t363 = t343 * t167 - t337 * t189;
t161 = m(3) * t298 - qJDD(2) * mrSges(3,2) + t319 * mrSges(3,3) - qJD(2) * t321 + t317 * t366 + t363;
t364 = -t338 * t160 + t344 * t161;
t288 = Ifges(4,5) * t310 + Ifges(4,6) * t309 + Ifges(4,3) * t331;
t150 = mrSges(4,2) * t284 - mrSges(4,3) * t250 + Ifges(4,1) * t283 + Ifges(4,4) * t282 + Ifges(4,5) * t330 - pkin(9) * t169 - t336 * t156 + t342 * t158 + t309 * t288 - t331 * t289;
t355 = -mrSges(7,1) * t198 + mrSges(7,2) * t199 - Ifges(7,5) * t212 - Ifges(7,6) * t211 - Ifges(7,3) * t275 - t247 * t223 + t246 * t224;
t351 = -mrSges(6,1) * t203 + mrSges(6,2) * t204 - Ifges(6,5) * t232 - Ifges(6,6) * t231 - Ifges(6,3) * t276 - pkin(5) * t185 - t265 * t243 + t264 * t244 + t355;
t347 = mrSges(5,1) * t219 - mrSges(5,2) * t220 + Ifges(5,5) * t254 + Ifges(5,6) * t253 + Ifges(5,3) * t281 + pkin(4) * t178 + t296 * t259 - t295 * t260 - t351;
t154 = -mrSges(4,1) * t284 + mrSges(4,3) * t251 + Ifges(4,4) * t283 + Ifges(4,2) * t282 + Ifges(4,6) * t330 - pkin(3) * t169 - t310 * t288 + t331 * t290 - t347;
t306 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t338 + Ifges(3,6) * t344) * qJD(1);
t311 = -t346 * pkin(7) + t357;
t354 = m(4) * t284 - t282 * mrSges(4,1) + t283 * mrSges(4,2) - t309 * t299 + t310 * t300 + t169;
t146 = -mrSges(3,1) * t311 + mrSges(3,3) * t298 + Ifges(3,4) * t318 + Ifges(3,2) * t319 + Ifges(3,6) * qJDD(2) - pkin(2) * t354 + pkin(8) * t363 + qJD(2) * t308 + t337 * t150 + t343 * t154 - t306 * t367;
t149 = mrSges(3,2) * t311 - mrSges(3,3) * t297 + Ifges(3,1) * t318 + Ifges(3,4) * t319 + Ifges(3,5) * qJDD(2) - pkin(8) * t162 - qJD(2) * t307 + t343 * t150 - t337 * t154 + t306 * t366;
t350 = -m(3) * t311 + t319 * mrSges(3,1) - t318 * mrSges(3,2) - t321 * t367 + t322 * t366 - t354;
t356 = mrSges(2,1) * t324 - mrSges(2,2) * t325 + Ifges(2,3) * qJDD(1) + pkin(1) * t350 + pkin(7) * t364 + t344 * t146 + t338 * t149;
t163 = m(2) * t324 + qJDD(1) * mrSges(2,1) - t346 * mrSges(2,2) + t350;
t153 = t344 * t160 + t338 * t161;
t151 = m(2) * t325 - t346 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t364;
t147 = mrSges(2,1) * g(3) + mrSges(2,3) * t325 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t153 - t370;
t144 = -mrSges(2,2) * g(3) - mrSges(2,3) * t324 + Ifges(2,5) * qJDD(1) - t346 * Ifges(2,6) - pkin(7) * t153 - t338 * t146 + t344 * t149;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t144 - t339 * t147 - pkin(6) * (t339 * t151 + t345 * t163), t144, t149, t150, t158, t172, t187; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t339 * t144 + t345 * t147 + pkin(6) * (t345 * t151 - t339 * t163), t147, t146, t154, t156, t171, t186; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t370, -t353, t347, -t351, -t355;];
m_new  = t1;
