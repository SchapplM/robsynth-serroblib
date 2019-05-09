% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP6
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
% Datum: 2019-05-08 05:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:08:23
% EndTime: 2019-05-08 05:09:30
% DurationCPUTime: 23.74s
% Computational Cost: add. (437097->381), mult. (891470->462), div. (0->0), fcn. (650821->10), ass. (0->143)
t328 = sin(qJ(1));
t332 = cos(qJ(1));
t317 = -g(1) * t332 - g(2) * t328;
t334 = qJD(1) ^ 2;
t299 = -pkin(1) * t334 + qJDD(1) * pkin(7) + t317;
t327 = sin(qJ(2));
t331 = cos(qJ(2));
t282 = -t331 * g(3) - t327 * t299;
t308 = (-pkin(2) * t331 - pkin(8) * t327) * qJD(1);
t333 = qJD(2) ^ 2;
t355 = qJD(1) * t327;
t262 = -qJDD(2) * pkin(2) - t333 * pkin(8) + t308 * t355 - t282;
t326 = sin(qJ(3));
t330 = cos(qJ(3));
t306 = qJD(2) * t326 + t330 * t355;
t353 = qJD(1) * qJD(2);
t351 = t331 * t353;
t309 = qJDD(1) * t327 + t351;
t274 = -qJD(3) * t306 + qJDD(2) * t330 - t309 * t326;
t354 = qJD(1) * t331;
t319 = qJD(3) - t354;
t284 = pkin(3) * t319 - pkin(9) * t306;
t305 = qJD(2) * t330 - t326 * t355;
t303 = t305 ^ 2;
t234 = -t274 * pkin(3) - t303 * pkin(9) + t306 * t284 + t262;
t275 = qJD(3) * t305 + qJDD(2) * t326 + t309 * t330;
t325 = sin(qJ(4));
t329 = cos(qJ(4));
t278 = t305 * t325 + t306 * t329;
t240 = -qJD(4) * t278 + t274 * t329 - t275 * t325;
t318 = qJD(4) + t319;
t266 = pkin(4) * t318 - pkin(10) * t278;
t277 = t305 * t329 - t306 * t325;
t276 = t277 ^ 2;
t200 = -t240 * pkin(4) - t276 * pkin(10) + t278 * t266 + t234;
t241 = qJD(4) * t277 + t274 * t325 + t275 * t329;
t324 = sin(qJ(5));
t359 = cos(qJ(5));
t256 = t324 * t277 + t359 * t278;
t211 = t256 * qJD(5) - t359 * t240 + t324 * t241;
t255 = -t359 * t277 + t324 * t278;
t212 = -t255 * qJD(5) + t324 * t240 + t359 * t241;
t313 = qJD(5) + t318;
t191 = t200 - 0.2e1 * qJD(6) * t256 + (t255 * t313 - t212) * qJ(6) + (t256 * t313 + t211) * pkin(5);
t246 = -mrSges(7,2) * t255 + mrSges(7,3) * t313;
t249 = -mrSges(7,1) * t313 + mrSges(7,2) * t256;
t181 = m(7) * t191 + t211 * mrSges(7,1) - t212 * mrSges(7,3) + t255 * t246 - t256 * t249;
t316 = t328 * g(1) - t332 * g(2);
t298 = -qJDD(1) * pkin(1) - t334 * pkin(7) - t316;
t320 = t327 * t353;
t310 = qJDD(1) * t331 - t320;
t260 = (-t309 - t351) * pkin(8) + (-t310 + t320) * pkin(2) + t298;
t283 = -g(3) * t327 + t331 * t299;
t263 = -pkin(2) * t333 + qJDD(2) * pkin(8) + t308 * t354 + t283;
t242 = t330 * t260 - t326 * t263;
t304 = qJDD(3) - t310;
t217 = (t305 * t319 - t275) * pkin(9) + (t305 * t306 + t304) * pkin(3) + t242;
t243 = t326 * t260 + t330 * t263;
t219 = -pkin(3) * t303 + pkin(9) * t274 - t284 * t319 + t243;
t197 = t329 * t217 - t325 * t219;
t300 = qJDD(4) + t304;
t193 = (t277 * t318 - t241) * pkin(10) + (t277 * t278 + t300) * pkin(4) + t197;
t198 = t325 * t217 + t329 * t219;
t195 = -pkin(4) * t276 + pkin(10) * t240 - t266 * t318 + t198;
t189 = t324 * t193 + t359 * t195;
t225 = Ifges(7,1) * t256 + Ifges(7,4) * t313 + Ifges(7,5) * t255;
t226 = Ifges(6,1) * t256 - Ifges(6,4) * t255 + Ifges(6,5) * t313;
t294 = qJDD(5) + t300;
t231 = pkin(5) * t255 - qJ(6) * t256;
t312 = t313 ^ 2;
t184 = -pkin(5) * t312 + qJ(6) * t294 + 0.2e1 * qJD(6) * t313 - t231 * t255 + t189;
t346 = -mrSges(7,1) * t191 + mrSges(7,2) * t184;
t223 = Ifges(7,4) * t256 + Ifges(7,2) * t313 + Ifges(7,6) * t255;
t357 = -Ifges(6,5) * t256 + Ifges(6,6) * t255 - Ifges(6,3) * t313 - t223;
t166 = -mrSges(6,1) * t200 + mrSges(6,3) * t189 - pkin(5) * t181 + (t225 + t226) * t313 + (Ifges(6,6) - Ifges(7,6)) * t294 + t357 * t256 + (Ifges(6,4) - Ifges(7,5)) * t212 + (-Ifges(6,2) - Ifges(7,3)) * t211 + t346;
t188 = t359 * t193 - t324 * t195;
t224 = Ifges(6,4) * t256 - Ifges(6,2) * t255 + Ifges(6,6) * t313;
t186 = -t294 * pkin(5) - t312 * qJ(6) + t256 * t231 + qJDD(6) - t188;
t221 = Ifges(7,5) * t256 + Ifges(7,6) * t313 + Ifges(7,3) * t255;
t344 = mrSges(7,2) * t186 - mrSges(7,3) * t191 + Ifges(7,1) * t212 + Ifges(7,4) * t294 + Ifges(7,5) * t211 + t313 * t221;
t167 = mrSges(6,2) * t200 - mrSges(6,3) * t188 + Ifges(6,1) * t212 - Ifges(6,4) * t211 + Ifges(6,5) * t294 - qJ(6) * t181 - t313 * t224 + t357 * t255 + t344;
t250 = Ifges(5,5) * t278 + Ifges(5,6) * t277 + Ifges(5,3) * t318;
t252 = Ifges(5,1) * t278 + Ifges(5,4) * t277 + Ifges(5,5) * t318;
t247 = -mrSges(6,2) * t313 - mrSges(6,3) * t255;
t248 = mrSges(6,1) * t313 - mrSges(6,3) * t256;
t342 = m(6) * t200 + t211 * mrSges(6,1) + t212 * mrSges(6,2) + t255 * t247 + t256 * t248 + t181;
t352 = m(7) * t184 + t294 * mrSges(7,3) + t313 * t249;
t232 = mrSges(7,1) * t255 - mrSges(7,3) * t256;
t356 = -mrSges(6,1) * t255 - mrSges(6,2) * t256 - t232;
t358 = -mrSges(6,3) - mrSges(7,2);
t174 = m(6) * t189 - t294 * mrSges(6,2) + t358 * t211 - t313 * t248 + t356 * t255 + t352;
t347 = -m(7) * t186 + t294 * mrSges(7,1) + t313 * t246;
t176 = m(6) * t188 + t294 * mrSges(6,1) + t358 * t212 + t313 * t247 + t356 * t256 + t347;
t348 = t359 * t174 - t176 * t324;
t155 = -mrSges(5,1) * t234 + mrSges(5,3) * t198 + Ifges(5,4) * t241 + Ifges(5,2) * t240 + Ifges(5,6) * t300 - pkin(4) * t342 + pkin(10) * t348 + t359 * t166 + t324 * t167 - t278 * t250 + t318 * t252;
t169 = t324 * t174 + t359 * t176;
t251 = Ifges(5,4) * t278 + Ifges(5,2) * t277 + Ifges(5,6) * t318;
t156 = mrSges(5,2) * t234 - mrSges(5,3) * t197 + Ifges(5,1) * t241 + Ifges(5,4) * t240 + Ifges(5,5) * t300 - pkin(10) * t169 - t324 * t166 + t359 * t167 + t277 * t250 - t318 * t251;
t267 = Ifges(4,5) * t306 + Ifges(4,6) * t305 + Ifges(4,3) * t319;
t269 = Ifges(4,1) * t306 + Ifges(4,4) * t305 + Ifges(4,5) * t319;
t264 = -mrSges(5,2) * t318 + mrSges(5,3) * t277;
t265 = mrSges(5,1) * t318 - mrSges(5,3) * t278;
t339 = m(5) * t234 - t240 * mrSges(5,1) + t241 * mrSges(5,2) - t277 * t264 + t278 * t265 + t342;
t257 = -mrSges(5,1) * t277 + mrSges(5,2) * t278;
t164 = m(5) * t197 + mrSges(5,1) * t300 - mrSges(5,3) * t241 - t257 * t278 + t264 * t318 + t169;
t165 = m(5) * t198 - mrSges(5,2) * t300 + mrSges(5,3) * t240 + t257 * t277 - t265 * t318 + t348;
t349 = -t164 * t325 + t329 * t165;
t142 = -mrSges(4,1) * t262 + mrSges(4,3) * t243 + Ifges(4,4) * t275 + Ifges(4,2) * t274 + Ifges(4,6) * t304 - pkin(3) * t339 + pkin(9) * t349 + t329 * t155 + t325 * t156 - t306 * t267 + t319 * t269;
t160 = t329 * t164 + t325 * t165;
t268 = Ifges(4,4) * t306 + Ifges(4,2) * t305 + Ifges(4,6) * t319;
t145 = mrSges(4,2) * t262 - mrSges(4,3) * t242 + Ifges(4,1) * t275 + Ifges(4,4) * t274 + Ifges(4,5) * t304 - pkin(9) * t160 - t155 * t325 + t156 * t329 + t267 * t305 - t268 * t319;
t279 = -mrSges(4,1) * t305 + mrSges(4,2) * t306;
t280 = -mrSges(4,2) * t319 + mrSges(4,3) * t305;
t158 = m(4) * t242 + mrSges(4,1) * t304 - mrSges(4,3) * t275 - t279 * t306 + t280 * t319 + t160;
t281 = mrSges(4,1) * t319 - mrSges(4,3) * t306;
t159 = m(4) * t243 - mrSges(4,2) * t304 + mrSges(4,3) * t274 + t279 * t305 - t281 * t319 + t349;
t154 = -t158 * t326 + t330 * t159;
t171 = -m(4) * t262 + t274 * mrSges(4,1) - t275 * mrSges(4,2) + t305 * t280 - t306 * t281 - t339;
t296 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t327 + Ifges(3,2) * t331) * qJD(1);
t297 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t327 + Ifges(3,4) * t331) * qJD(1);
t360 = mrSges(3,1) * t282 - mrSges(3,2) * t283 + Ifges(3,5) * t309 + Ifges(3,6) * t310 + Ifges(3,3) * qJDD(2) + pkin(2) * t171 + pkin(8) * t154 + t330 * t142 + t326 * t145 + (t296 * t327 - t297 * t331) * qJD(1);
t307 = (-mrSges(3,1) * t331 + mrSges(3,2) * t327) * qJD(1);
t314 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t355;
t152 = m(3) * t283 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t310 - qJD(2) * t314 + t307 * t354 + t154;
t315 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t354;
t170 = m(3) * t282 + qJDD(2) * mrSges(3,1) - t309 * mrSges(3,3) + qJD(2) * t315 - t307 * t355 + t171;
t350 = t331 * t152 - t170 * t327;
t153 = t158 * t330 + t159 * t326;
t295 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t327 + Ifges(3,6) * t331) * qJD(1);
t141 = mrSges(3,2) * t298 - mrSges(3,3) * t282 + Ifges(3,1) * t309 + Ifges(3,4) * t310 + Ifges(3,5) * qJDD(2) - pkin(8) * t153 - qJD(2) * t296 - t142 * t326 + t145 * t330 + t295 * t354;
t341 = mrSges(7,1) * t186 - mrSges(7,3) * t184 - Ifges(7,4) * t212 - Ifges(7,2) * t294 - Ifges(7,6) * t211 + t256 * t221 - t255 * t225;
t338 = mrSges(6,2) * t189 - t255 * t226 - qJ(6) * (-t211 * mrSges(7,2) - t255 * t232 + t352) - pkin(5) * (-t212 * mrSges(7,2) - t256 * t232 + t347) - mrSges(6,1) * t188 + Ifges(6,6) * t211 - Ifges(6,5) * t212 - t256 * t224 - Ifges(6,3) * t294 + t341;
t336 = -mrSges(5,1) * t197 + mrSges(5,2) * t198 - Ifges(5,5) * t241 - Ifges(5,6) * t240 - Ifges(5,3) * t300 - pkin(4) * t169 - t278 * t251 + t277 * t252 + t338;
t335 = mrSges(4,1) * t242 - mrSges(4,2) * t243 + Ifges(4,5) * t275 + Ifges(4,6) * t274 + Ifges(4,3) * t304 + pkin(3) * t160 + t306 * t268 - t305 * t269 - t336;
t144 = -mrSges(3,1) * t298 + mrSges(3,3) * t283 + Ifges(3,4) * t309 + Ifges(3,2) * t310 + Ifges(3,6) * qJDD(2) - pkin(2) * t153 + qJD(2) * t297 - t295 * t355 - t335;
t340 = -m(3) * t298 + t310 * mrSges(3,1) - mrSges(3,2) * t309 - t314 * t355 + t315 * t354 - t153;
t343 = mrSges(2,1) * t316 - mrSges(2,2) * t317 + Ifges(2,3) * qJDD(1) + pkin(1) * t340 + pkin(7) * t350 + t327 * t141 + t331 * t144;
t149 = m(2) * t316 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t334 + t340;
t148 = t152 * t327 + t170 * t331;
t146 = m(2) * t317 - mrSges(2,1) * t334 - qJDD(1) * mrSges(2,2) + t350;
t139 = mrSges(2,1) * g(3) + mrSges(2,3) * t317 + t334 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t148 - t360;
t138 = -mrSges(2,2) * g(3) - mrSges(2,3) * t316 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t334 - pkin(7) * t148 + t141 * t331 - t144 * t327;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t332 * t138 - t328 * t139 - pkin(6) * (t146 * t328 + t149 * t332), t138, t141, t145, t156, t167, -t223 * t255 + t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t328 * t138 + t332 * t139 + pkin(6) * (t146 * t332 - t149 * t328), t139, t144, t142, t155, t166, -t341; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t343, t343, t360, t335, -t336, -t338, Ifges(7,5) * t212 + Ifges(7,6) * t294 + Ifges(7,3) * t211 + t256 * t223 - t313 * t225 - t346;];
m_new  = t1;
