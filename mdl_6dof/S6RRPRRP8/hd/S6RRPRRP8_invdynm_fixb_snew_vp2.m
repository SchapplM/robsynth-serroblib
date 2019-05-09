% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:19:19
% EndTime: 2019-05-06 18:20:03
% DurationCPUTime: 21.79s
% Computational Cost: add. (383572->381), mult. (834706->464), div. (0->0), fcn. (604294->10), ass. (0->141)
t329 = sin(qJ(1));
t332 = cos(qJ(1));
t316 = -g(1) * t332 - g(2) * t329;
t334 = qJD(1) ^ 2;
t299 = -pkin(1) * t334 + qJDD(1) * pkin(7) + t316;
t328 = sin(qJ(2));
t331 = cos(qJ(2));
t277 = -t331 * g(3) - t328 * t299;
t308 = (-pkin(2) * t331 - qJ(3) * t328) * qJD(1);
t333 = qJD(2) ^ 2;
t355 = qJD(1) * t328;
t262 = -qJDD(2) * pkin(2) - t333 * qJ(3) + t308 * t355 + qJDD(3) - t277;
t353 = qJD(1) * qJD(2);
t351 = t331 * t353;
t310 = qJDD(1) * t328 + t351;
t324 = sin(pkin(10));
t325 = cos(pkin(10));
t282 = qJDD(2) * t325 - t310 * t324;
t305 = qJD(2) * t324 + t325 * t355;
t354 = qJD(1) * t331;
t284 = -pkin(3) * t354 - pkin(8) * t305;
t304 = qJD(2) * t325 - t324 * t355;
t303 = t304 ^ 2;
t236 = -t282 * pkin(3) - t303 * pkin(8) + t305 * t284 + t262;
t327 = sin(qJ(4));
t330 = cos(qJ(4));
t275 = t304 * t327 + t305 * t330;
t283 = qJDD(2) * t324 + t310 * t325;
t248 = -qJD(4) * t275 + t282 * t330 - t283 * t327;
t319 = qJD(4) - t354;
t266 = pkin(4) * t319 - pkin(9) * t275;
t274 = t304 * t330 - t305 * t327;
t273 = t274 ^ 2;
t200 = -t248 * pkin(4) - t273 * pkin(9) + t275 * t266 + t236;
t249 = qJD(4) * t274 + t282 * t327 + t283 * t330;
t326 = sin(qJ(5));
t359 = cos(qJ(5));
t256 = t326 * t274 + t359 * t275;
t211 = t256 * qJD(5) - t359 * t248 + t326 * t249;
t255 = -t359 * t274 + t326 * t275;
t212 = -t255 * qJD(5) + t326 * t248 + t359 * t249;
t318 = qJD(5) + t319;
t191 = t200 + (t255 * t318 - t212) * qJ(6) + (t256 * t318 + t211) * pkin(5) - 0.2e1 * qJD(6) * t256;
t241 = -mrSges(7,2) * t255 + mrSges(7,3) * t318;
t244 = -mrSges(7,1) * t318 + mrSges(7,2) * t256;
t181 = m(7) * t191 + t211 * mrSges(7,1) - t212 * mrSges(7,3) + t255 * t241 - t256 * t244;
t315 = t329 * g(1) - t332 * g(2);
t298 = -qJDD(1) * pkin(1) - t334 * pkin(7) - t315;
t320 = t328 * t353;
t311 = qJDD(1) * t331 - t320;
t260 = (-t310 - t351) * qJ(3) + (-t311 + t320) * pkin(2) + t298;
t278 = -g(3) * t328 + t331 * t299;
t263 = -pkin(2) * t333 + qJDD(2) * qJ(3) + t308 * t354 + t278;
t234 = -0.2e1 * qJD(3) * t305 + t325 * t260 - t324 * t263;
t217 = (-t304 * t354 - t283) * pkin(8) + (t304 * t305 - t311) * pkin(3) + t234;
t235 = 0.2e1 * qJD(3) * t304 + t324 * t260 + t325 * t263;
t219 = -pkin(3) * t303 + pkin(8) * t282 + t284 * t354 + t235;
t197 = t330 * t217 - t327 * t219;
t307 = qJDD(4) - t311;
t193 = (t274 * t319 - t249) * pkin(9) + (t274 * t275 + t307) * pkin(4) + t197;
t198 = t327 * t217 + t330 * t219;
t195 = -pkin(4) * t273 + pkin(9) * t248 - t266 * t319 + t198;
t189 = t326 * t193 + t359 * t195;
t225 = Ifges(7,1) * t256 + Ifges(7,4) * t318 + Ifges(7,5) * t255;
t226 = Ifges(6,1) * t256 - Ifges(6,4) * t255 + Ifges(6,5) * t318;
t301 = qJDD(5) + t307;
t231 = pkin(5) * t255 - qJ(6) * t256;
t317 = t318 ^ 2;
t184 = -pkin(5) * t317 + qJ(6) * t301 + 0.2e1 * qJD(6) * t318 - t231 * t255 + t189;
t346 = -mrSges(7,1) * t191 + mrSges(7,2) * t184;
t223 = Ifges(7,4) * t256 + Ifges(7,2) * t318 + Ifges(7,6) * t255;
t357 = -Ifges(6,5) * t256 + Ifges(6,6) * t255 - Ifges(6,3) * t318 - t223;
t166 = -mrSges(6,1) * t200 + mrSges(6,3) * t189 - pkin(5) * t181 + (t225 + t226) * t318 + (Ifges(6,6) - Ifges(7,6)) * t301 + t357 * t256 + (Ifges(6,4) - Ifges(7,5)) * t212 + (-Ifges(6,2) - Ifges(7,3)) * t211 + t346;
t188 = t359 * t193 - t326 * t195;
t224 = Ifges(6,4) * t256 - Ifges(6,2) * t255 + Ifges(6,6) * t318;
t186 = -t301 * pkin(5) - t317 * qJ(6) + t256 * t231 + qJDD(6) - t188;
t221 = Ifges(7,5) * t256 + Ifges(7,6) * t318 + Ifges(7,3) * t255;
t344 = mrSges(7,2) * t186 - mrSges(7,3) * t191 + Ifges(7,1) * t212 + Ifges(7,4) * t301 + Ifges(7,5) * t211 + t318 * t221;
t167 = mrSges(6,2) * t200 - mrSges(6,3) * t188 + Ifges(6,1) * t212 - Ifges(6,4) * t211 + Ifges(6,5) * t301 - qJ(6) * t181 - t318 * t224 + t357 * t255 + t344;
t250 = Ifges(5,5) * t275 + Ifges(5,6) * t274 + Ifges(5,3) * t319;
t252 = Ifges(5,1) * t275 + Ifges(5,4) * t274 + Ifges(5,5) * t319;
t242 = -mrSges(6,2) * t318 - mrSges(6,3) * t255;
t243 = mrSges(6,1) * t318 - mrSges(6,3) * t256;
t342 = m(6) * t200 + t211 * mrSges(6,1) + t212 * mrSges(6,2) + t255 * t242 + t256 * t243 + t181;
t352 = m(7) * t184 + t301 * mrSges(7,3) + t318 * t244;
t232 = mrSges(7,1) * t255 - mrSges(7,3) * t256;
t356 = -mrSges(6,1) * t255 - mrSges(6,2) * t256 - t232;
t358 = -mrSges(6,3) - mrSges(7,2);
t172 = m(6) * t189 - t301 * mrSges(6,2) + t358 * t211 - t318 * t243 + t356 * t255 + t352;
t347 = -m(7) * t186 + t301 * mrSges(7,1) + t318 * t241;
t174 = m(6) * t188 + t301 * mrSges(6,1) + t358 * t212 + t318 * t242 + t356 * t256 + t347;
t348 = t359 * t172 - t174 * t326;
t155 = -mrSges(5,1) * t236 + mrSges(5,3) * t198 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t307 - pkin(4) * t342 + pkin(9) * t348 + t359 * t166 + t326 * t167 - t275 * t250 + t319 * t252;
t169 = t326 * t172 + t359 * t174;
t251 = Ifges(5,4) * t275 + Ifges(5,2) * t274 + Ifges(5,6) * t319;
t156 = mrSges(5,2) * t236 - mrSges(5,3) * t197 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t307 - pkin(9) * t169 - t326 * t166 + t359 * t167 + t274 * t250 - t319 * t251;
t268 = Ifges(4,5) * t305 + Ifges(4,6) * t304 - Ifges(4,3) * t354;
t270 = Ifges(4,1) * t305 + Ifges(4,4) * t304 - Ifges(4,5) * t354;
t264 = -mrSges(5,2) * t319 + mrSges(5,3) * t274;
t265 = mrSges(5,1) * t319 - mrSges(5,3) * t275;
t339 = m(5) * t236 - t248 * mrSges(5,1) + t249 * mrSges(5,2) - t274 * t264 + t275 * t265 + t342;
t257 = -mrSges(5,1) * t274 + mrSges(5,2) * t275;
t164 = m(5) * t197 + mrSges(5,1) * t307 - mrSges(5,3) * t249 - t257 * t275 + t264 * t319 + t169;
t165 = m(5) * t198 - mrSges(5,2) * t307 + mrSges(5,3) * t248 + t257 * t274 - t265 * t319 + t348;
t349 = -t164 * t327 + t330 * t165;
t142 = -mrSges(4,1) * t262 + mrSges(4,3) * t235 + Ifges(4,4) * t283 + Ifges(4,2) * t282 - Ifges(4,6) * t311 - pkin(3) * t339 + pkin(8) * t349 + t330 * t155 + t327 * t156 - t305 * t268 - t270 * t354;
t160 = t330 * t164 + t327 * t165;
t269 = Ifges(4,4) * t305 + Ifges(4,2) * t304 - Ifges(4,6) * t354;
t145 = mrSges(4,2) * t262 - mrSges(4,3) * t234 + Ifges(4,1) * t283 + Ifges(4,4) * t282 - Ifges(4,5) * t311 - pkin(8) * t160 - t155 * t327 + t156 * t330 + t268 * t304 + t269 * t354;
t276 = -mrSges(4,1) * t304 + mrSges(4,2) * t305;
t280 = mrSges(4,2) * t354 + mrSges(4,3) * t304;
t158 = m(4) * t234 - mrSges(4,1) * t311 - mrSges(4,3) * t283 - t276 * t305 - t280 * t354 + t160;
t281 = -mrSges(4,1) * t354 - mrSges(4,3) * t305;
t159 = m(4) * t235 + mrSges(4,2) * t311 + mrSges(4,3) * t282 + t276 * t304 + t281 * t354 + t349;
t154 = -t158 * t324 + t325 * t159;
t176 = -m(4) * t262 + t282 * mrSges(4,1) - t283 * mrSges(4,2) + t304 * t280 - t305 * t281 - t339;
t296 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t328 + Ifges(3,2) * t331) * qJD(1);
t297 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t328 + Ifges(3,4) * t331) * qJD(1);
t360 = mrSges(3,1) * t277 - mrSges(3,2) * t278 + Ifges(3,5) * t310 + Ifges(3,6) * t311 + Ifges(3,3) * qJDD(2) + pkin(2) * t176 + qJ(3) * t154 + t325 * t142 + t324 * t145 + (t296 * t328 - t297 * t331) * qJD(1);
t309 = (-mrSges(3,1) * t331 + mrSges(3,2) * t328) * qJD(1);
t313 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t355;
t152 = m(3) * t278 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t311 - qJD(2) * t313 + t309 * t354 + t154;
t314 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t354;
t175 = m(3) * t277 + qJDD(2) * mrSges(3,1) - t310 * mrSges(3,3) + qJD(2) * t314 - t309 * t355 + t176;
t350 = t331 * t152 - t175 * t328;
t153 = t158 * t325 + t159 * t324;
t295 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t328 + Ifges(3,6) * t331) * qJD(1);
t141 = mrSges(3,2) * t298 - mrSges(3,3) * t277 + Ifges(3,1) * t310 + Ifges(3,4) * t311 + Ifges(3,5) * qJDD(2) - qJ(3) * t153 - qJD(2) * t296 - t142 * t324 + t145 * t325 + t295 * t354;
t341 = mrSges(7,1) * t186 - mrSges(7,3) * t184 - Ifges(7,4) * t212 - Ifges(7,2) * t301 - Ifges(7,6) * t211 + t256 * t221 - t255 * t225;
t338 = mrSges(6,2) * t189 - t255 * t226 - qJ(6) * (-t211 * mrSges(7,2) - t255 * t232 + t352) - pkin(5) * (-t212 * mrSges(7,2) - t256 * t232 + t347) - mrSges(6,1) * t188 + Ifges(6,6) * t211 - Ifges(6,5) * t212 - t256 * t224 - Ifges(6,3) * t301 + t341;
t336 = -mrSges(5,1) * t197 + mrSges(5,2) * t198 - Ifges(5,5) * t249 - Ifges(5,6) * t248 - Ifges(5,3) * t307 - pkin(4) * t169 - t275 * t251 + t274 * t252 + t338;
t335 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t283 + Ifges(4,6) * t282 + pkin(3) * t160 + t305 * t269 - t304 * t270 - t336;
t144 = Ifges(3,6) * qJDD(2) - pkin(2) * t153 - t335 + (Ifges(3,2) + Ifges(4,3)) * t311 + Ifges(3,4) * t310 + qJD(2) * t297 - mrSges(3,1) * t298 + mrSges(3,3) * t278 - t295 * t355;
t340 = -m(3) * t298 + t311 * mrSges(3,1) - mrSges(3,2) * t310 - t313 * t355 + t314 * t354 - t153;
t343 = mrSges(2,1) * t315 - mrSges(2,2) * t316 + Ifges(2,3) * qJDD(1) + pkin(1) * t340 + pkin(7) * t350 + t328 * t141 + t331 * t144;
t149 = m(2) * t315 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t334 + t340;
t148 = t152 * t328 + t175 * t331;
t146 = m(2) * t316 - mrSges(2,1) * t334 - qJDD(1) * mrSges(2,2) + t350;
t139 = mrSges(2,1) * g(3) + mrSges(2,3) * t316 + t334 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t148 - t360;
t138 = -mrSges(2,2) * g(3) - mrSges(2,3) * t315 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t334 - pkin(7) * t148 + t141 * t331 - t144 * t328;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t332 * t138 - t329 * t139 - pkin(6) * (t146 * t329 + t149 * t332), t138, t141, t145, t156, t167, -t223 * t255 + t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t329 * t138 + t332 * t139 + pkin(6) * (t146 * t332 - t149 * t329), t139, t144, t142, t155, t166, -t341; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t343, t343, t360, -Ifges(4,3) * t311 + t335, -t336, -t338, Ifges(7,5) * t212 + Ifges(7,6) * t301 + Ifges(7,3) * t211 + t256 * t223 - t318 * t225 - t346;];
m_new  = t1;
