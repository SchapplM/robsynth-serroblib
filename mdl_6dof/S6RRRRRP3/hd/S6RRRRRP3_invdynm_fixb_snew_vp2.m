% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 04:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:37:58
% EndTime: 2019-05-08 04:38:54
% DurationCPUTime: 21.68s
% Computational Cost: add. (404573->381), mult. (805512->467), div. (0->0), fcn. (587042->10), ass. (0->144)
t340 = sin(qJ(3));
t345 = cos(qJ(3));
t346 = cos(qJ(2));
t372 = qJD(1) * t346;
t341 = sin(qJ(2));
t373 = qJD(1) * t341;
t313 = -t340 * t373 + t345 * t372;
t371 = qJD(1) * qJD(2);
t322 = qJDD(1) * t341 + t346 * t371;
t323 = qJDD(1) * t346 - t341 * t371;
t288 = qJD(3) * t313 + t322 * t345 + t323 * t340;
t314 = (t340 * t346 + t341 * t345) * qJD(1);
t335 = qJD(2) + qJD(3);
t339 = sin(qJ(4));
t344 = cos(qJ(4));
t301 = t314 * t344 + t335 * t339;
t334 = qJDD(2) + qJDD(3);
t254 = -qJD(4) * t301 - t288 * t339 + t334 * t344;
t300 = -t314 * t339 + t335 * t344;
t255 = qJD(4) * t300 + t288 * t344 + t334 * t339;
t338 = sin(qJ(5));
t343 = cos(qJ(5));
t267 = t300 * t343 - t301 * t338;
t226 = qJD(5) * t267 + t254 * t338 + t255 * t343;
t268 = t300 * t338 + t301 * t343;
t247 = -mrSges(7,1) * t267 + mrSges(7,2) * t268;
t287 = -t314 * qJD(3) - t340 * t322 + t323 * t345;
t327 = qJD(2) * pkin(2) - pkin(8) * t373;
t337 = t346 ^ 2;
t348 = qJD(1) ^ 2;
t342 = sin(qJ(1));
t347 = cos(qJ(1));
t328 = g(1) * t342 - t347 * g(2);
t361 = -qJDD(1) * pkin(1) - t328;
t289 = -pkin(2) * t323 + t327 * t373 + (-pkin(8) * t337 - pkin(7)) * t348 + t361;
t232 = (-t313 * t335 - t288) * pkin(9) + (t314 * t335 - t287) * pkin(3) + t289;
t329 = -g(1) * t347 - g(2) * t342;
t316 = -pkin(1) * t348 + qJDD(1) * pkin(7) + t329;
t375 = t316 * t341;
t376 = pkin(2) * t348;
t274 = qJDD(2) * pkin(2) - pkin(8) * t322 - t375 + (pkin(8) * t371 + t341 * t376 - g(3)) * t346;
t303 = -g(3) * t341 + t346 * t316;
t275 = pkin(8) * t323 - qJD(2) * t327 - t337 * t376 + t303;
t251 = t340 * t274 + t345 * t275;
t298 = -pkin(3) * t313 - pkin(9) * t314;
t333 = t335 ^ 2;
t236 = -pkin(3) * t333 + pkin(9) * t334 + t298 * t313 + t251;
t210 = t344 * t232 - t236 * t339;
t286 = qJDD(4) - t287;
t309 = qJD(4) - t313;
t207 = (t300 * t309 - t255) * pkin(10) + (t300 * t301 + t286) * pkin(4) + t210;
t211 = t339 * t232 + t344 * t236;
t292 = pkin(4) * t309 - pkin(10) * t301;
t299 = t300 ^ 2;
t209 = -pkin(4) * t299 + pkin(10) * t254 - t292 * t309 + t211;
t200 = t343 * t207 - t209 * t338;
t281 = qJDD(5) + t286;
t307 = qJD(5) + t309;
t195 = -0.2e1 * qJD(6) * t268 + (t267 * t307 - t226) * qJ(6) + (t267 * t268 + t281) * pkin(5) + t200;
t256 = -mrSges(7,2) * t307 + mrSges(7,3) * t267;
t370 = m(7) * t195 + t281 * mrSges(7,1) + t307 * t256;
t192 = -mrSges(7,3) * t226 - t247 * t268 + t370;
t201 = t338 * t207 + t343 * t209;
t225 = -qJD(5) * t268 + t254 * t343 - t255 * t338;
t241 = Ifges(6,4) * t268 + Ifges(6,2) * t267 + Ifges(6,6) * t307;
t242 = Ifges(7,1) * t268 + Ifges(7,4) * t267 + Ifges(7,5) * t307;
t243 = Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * t307;
t258 = pkin(5) * t307 - qJ(6) * t268;
t266 = t267 ^ 2;
t198 = -pkin(5) * t266 + qJ(6) * t225 + 0.2e1 * qJD(6) * t267 - t258 * t307 + t201;
t240 = Ifges(7,4) * t268 + Ifges(7,2) * t267 + Ifges(7,6) * t307;
t359 = -mrSges(7,1) * t195 + mrSges(7,2) * t198 - Ifges(7,5) * t226 - Ifges(7,6) * t225 - Ifges(7,3) * t281 - t268 * t240;
t379 = mrSges(6,1) * t200 - mrSges(6,2) * t201 + Ifges(6,5) * t226 + Ifges(6,6) * t225 + Ifges(6,3) * t281 + pkin(5) * t192 + t268 * t241 - t359 - (t243 + t242) * t267;
t248 = -mrSges(6,1) * t267 + mrSges(6,2) * t268;
t257 = -mrSges(6,2) * t307 + mrSges(6,3) * t267;
t183 = m(6) * t200 + mrSges(6,1) * t281 + t257 * t307 + (-t247 - t248) * t268 + (-mrSges(6,3) - mrSges(7,3)) * t226 + t370;
t259 = mrSges(7,1) * t307 - mrSges(7,3) * t268;
t260 = mrSges(6,1) * t307 - mrSges(6,3) * t268;
t369 = m(7) * t198 + t225 * mrSges(7,3) + t267 * t247;
t186 = m(6) * t201 + mrSges(6,3) * t225 + t248 * t267 + (-t259 - t260) * t307 + (-mrSges(6,2) - mrSges(7,2)) * t281 + t369;
t181 = t343 * t183 + t338 * t186;
t262 = Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * t309;
t263 = Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * t309;
t378 = mrSges(5,1) * t210 - mrSges(5,2) * t211 + Ifges(5,5) * t255 + Ifges(5,6) * t254 + Ifges(5,3) * t286 + pkin(4) * t181 + t301 * t262 - t300 * t263 + t379;
t297 = -mrSges(4,1) * t313 + mrSges(4,2) * t314;
t305 = mrSges(4,1) * t335 - mrSges(4,3) * t314;
t272 = -mrSges(5,1) * t300 + mrSges(5,2) * t301;
t290 = -mrSges(5,2) * t309 + mrSges(5,3) * t300;
t178 = m(5) * t210 + mrSges(5,1) * t286 - mrSges(5,3) * t255 - t272 * t301 + t290 * t309 + t181;
t291 = mrSges(5,1) * t309 - mrSges(5,3) * t301;
t364 = -t183 * t338 + t343 * t186;
t179 = m(5) * t211 - mrSges(5,2) * t286 + mrSges(5,3) * t254 + t272 * t300 - t291 * t309 + t364;
t365 = -t178 * t339 + t344 * t179;
t170 = m(4) * t251 - mrSges(4,2) * t334 + mrSges(4,3) * t287 + t297 * t313 - t305 * t335 + t365;
t250 = t274 * t345 - t340 * t275;
t304 = -mrSges(4,2) * t335 + mrSges(4,3) * t313;
t235 = -pkin(3) * t334 - pkin(9) * t333 + t314 * t298 - t250;
t212 = -pkin(4) * t254 - pkin(10) * t299 + t301 * t292 + t235;
t204 = -pkin(5) * t225 - qJ(6) * t266 + t258 * t268 + qJDD(6) + t212;
t363 = m(7) * t204 - t225 * mrSges(7,1) + t226 * mrSges(7,2) - t267 * t256 + t268 * t259;
t354 = m(6) * t212 - t225 * mrSges(6,1) + mrSges(6,2) * t226 - t267 * t257 + t260 * t268 + t363;
t351 = -m(5) * t235 + t254 * mrSges(5,1) - mrSges(5,2) * t255 + t300 * t290 - t291 * t301 - t354;
t188 = m(4) * t250 + mrSges(4,1) * t334 - mrSges(4,3) * t288 - t297 * t314 + t304 * t335 + t351;
t165 = t340 * t170 + t345 * t188;
t302 = -g(3) * t346 - t375;
t311 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t341 + Ifges(3,2) * t346) * qJD(1);
t312 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t341 + Ifges(3,4) * t346) * qJD(1);
t238 = Ifges(7,5) * t268 + Ifges(7,6) * t267 + Ifges(7,3) * t307;
t239 = Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * t307;
t360 = -mrSges(7,1) * t204 + mrSges(7,3) * t198 + Ifges(7,4) * t226 + Ifges(7,2) * t225 + Ifges(7,6) * t281 + t307 * t242;
t174 = Ifges(6,4) * t226 + Ifges(6,2) * t225 + Ifges(6,6) * t281 + t307 * t243 - mrSges(6,1) * t212 + mrSges(6,3) * t201 - pkin(5) * t363 + qJ(6) * (-mrSges(7,2) * t281 - t259 * t307 + t369) + (-t239 - t238) * t268 + t360;
t358 = mrSges(7,2) * t204 - mrSges(7,3) * t195 + Ifges(7,1) * t226 + Ifges(7,4) * t225 + Ifges(7,5) * t281 + t267 * t238;
t180 = mrSges(6,2) * t212 - mrSges(6,3) * t200 + Ifges(6,1) * t226 + Ifges(6,4) * t225 + Ifges(6,5) * t281 - qJ(6) * t192 + t239 * t267 + (-t240 - t241) * t307 + t358;
t261 = Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * t309;
t159 = -mrSges(5,1) * t235 + mrSges(5,3) * t211 + Ifges(5,4) * t255 + Ifges(5,2) * t254 + Ifges(5,6) * t286 - pkin(4) * t354 + pkin(10) * t364 + t343 * t174 + t338 * t180 - t301 * t261 + t309 * t263;
t161 = mrSges(5,2) * t235 - mrSges(5,3) * t210 + Ifges(5,1) * t255 + Ifges(5,4) * t254 + Ifges(5,5) * t286 - pkin(10) * t181 - t174 * t338 + t180 * t343 + t261 * t300 - t262 * t309;
t294 = Ifges(4,4) * t314 + Ifges(4,2) * t313 + Ifges(4,6) * t335;
t295 = Ifges(4,1) * t314 + Ifges(4,4) * t313 + Ifges(4,5) * t335;
t355 = -mrSges(4,1) * t250 + mrSges(4,2) * t251 - Ifges(4,5) * t288 - Ifges(4,6) * t287 - Ifges(4,3) * t334 - pkin(3) * t351 - pkin(9) * t365 - t344 * t159 - t339 * t161 - t314 * t294 + t313 * t295;
t377 = mrSges(3,1) * t302 - mrSges(3,2) * t303 + Ifges(3,5) * t322 + Ifges(3,6) * t323 + Ifges(3,3) * qJDD(2) + pkin(2) * t165 + (t341 * t311 - t346 * t312) * qJD(1) - t355;
t172 = t344 * t178 + t339 * t179;
t321 = (-mrSges(3,1) * t346 + mrSges(3,2) * t341) * qJD(1);
t326 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t372;
t163 = m(3) * t302 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t322 + qJD(2) * t326 - t321 * t373 + t165;
t325 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t373;
t366 = t345 * t170 - t340 * t188;
t164 = m(3) * t303 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t323 - qJD(2) * t325 + t321 * t372 + t366;
t367 = -t163 * t341 + t346 * t164;
t293 = Ifges(4,5) * t314 + Ifges(4,6) * t313 + Ifges(4,3) * t335;
t153 = mrSges(4,2) * t289 - mrSges(4,3) * t250 + Ifges(4,1) * t288 + Ifges(4,4) * t287 + Ifges(4,5) * t334 - pkin(9) * t172 - t159 * t339 + t161 * t344 + t293 * t313 - t294 * t335;
t157 = -mrSges(4,1) * t289 + mrSges(4,3) * t251 + Ifges(4,4) * t288 + Ifges(4,2) * t287 + Ifges(4,6) * t334 - pkin(3) * t172 - t314 * t293 + t335 * t295 - t378;
t310 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t341 + Ifges(3,6) * t346) * qJD(1);
t315 = -pkin(7) * t348 + t361;
t356 = m(4) * t289 - t287 * mrSges(4,1) + mrSges(4,2) * t288 - t313 * t304 + t305 * t314 + t172;
t149 = -mrSges(3,1) * t315 + mrSges(3,3) * t303 + Ifges(3,4) * t322 + Ifges(3,2) * t323 + Ifges(3,6) * qJDD(2) - pkin(2) * t356 + pkin(8) * t366 + qJD(2) * t312 + t340 * t153 + t345 * t157 - t310 * t373;
t152 = mrSges(3,2) * t315 - mrSges(3,3) * t302 + Ifges(3,1) * t322 + Ifges(3,4) * t323 + Ifges(3,5) * qJDD(2) - pkin(8) * t165 - qJD(2) * t311 + t153 * t345 - t157 * t340 + t310 * t372;
t352 = -m(3) * t315 + t323 * mrSges(3,1) - mrSges(3,2) * t322 - t325 * t373 + t326 * t372 - t356;
t357 = mrSges(2,1) * t328 - mrSges(2,2) * t329 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t367 + t346 * t149 + t341 * t152;
t166 = m(2) * t328 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t348 + t352;
t156 = t163 * t346 + t164 * t341;
t154 = m(2) * t329 - mrSges(2,1) * t348 - qJDD(1) * mrSges(2,2) + t367;
t150 = mrSges(2,1) * g(3) + mrSges(2,3) * t329 + t348 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t156 - t377;
t147 = -mrSges(2,2) * g(3) - mrSges(2,3) * t328 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t348 - pkin(7) * t156 - t149 * t341 + t346 * t152;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t347 * t147 - t342 * t150 - pkin(6) * (t342 * t154 + t347 * t166), t147, t152, t153, t161, t180, -t240 * t307 + t358; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t147 + t347 * t150 + pkin(6) * (t154 * t347 - t342 * t166), t150, t149, t157, t159, t174, -t268 * t238 + t360; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t357, t357, t377, -t355, t378, t379, -t267 * t242 - t359;];
m_new  = t1;
