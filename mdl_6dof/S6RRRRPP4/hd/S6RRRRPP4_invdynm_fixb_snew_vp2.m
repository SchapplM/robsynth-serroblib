% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:14:56
% EndTime: 2019-05-07 18:15:46
% DurationCPUTime: 23.69s
% Computational Cost: add. (423103->381), mult. (872770->464), div. (0->0), fcn. (631993->10), ass. (0->143)
t327 = sin(qJ(1));
t331 = cos(qJ(1));
t315 = -g(1) * t331 - g(2) * t327;
t333 = qJD(1) ^ 2;
t299 = -pkin(1) * t333 + qJDD(1) * pkin(7) + t315;
t326 = sin(qJ(2));
t330 = cos(qJ(2));
t283 = -t330 * g(3) - t326 * t299;
t308 = (-pkin(2) * t330 - pkin(8) * t326) * qJD(1);
t332 = qJD(2) ^ 2;
t355 = qJD(1) * t326;
t263 = -qJDD(2) * pkin(2) - pkin(8) * t332 + t308 * t355 - t283;
t325 = sin(qJ(3));
t329 = cos(qJ(3));
t306 = qJD(2) * t325 + t329 * t355;
t353 = qJD(1) * qJD(2);
t351 = t330 * t353;
t309 = qJDD(1) * t326 + t351;
t275 = -qJD(3) * t306 + qJDD(2) * t329 - t309 * t325;
t354 = qJD(1) * t330;
t318 = qJD(3) - t354;
t285 = pkin(3) * t318 - pkin(9) * t306;
t305 = qJD(2) * t329 - t325 * t355;
t303 = t305 ^ 2;
t234 = -pkin(3) * t275 - pkin(9) * t303 + t306 * t285 + t263;
t276 = qJD(3) * t305 + qJDD(2) * t325 + t309 * t329;
t324 = sin(qJ(4));
t328 = cos(qJ(4));
t279 = t305 * t324 + t306 * t328;
t240 = -qJD(4) * t279 + t275 * t328 - t276 * t324;
t317 = qJD(4) + t318;
t266 = pkin(4) * t317 - qJ(5) * t279;
t278 = t305 * t328 - t306 * t324;
t277 = t278 ^ 2;
t200 = -pkin(4) * t240 - qJ(5) * t277 + t279 * t266 + qJDD(5) + t234;
t241 = qJD(4) * t278 + t275 * t324 + t276 * t328;
t323 = sin(pkin(10));
t358 = cos(pkin(10));
t213 = -t358 * t240 + t241 * t323;
t214 = t323 * t240 + t358 * t241;
t256 = -t358 * t278 + t279 * t323;
t257 = t323 * t278 + t358 * t279;
t191 = t200 + (t256 * t317 - t214) * qJ(6) + (t257 * t317 + t213) * pkin(5) - 0.2e1 * qJD(6) * t257;
t246 = -mrSges(7,2) * t256 + mrSges(7,3) * t317;
t249 = -mrSges(7,1) * t317 + mrSges(7,2) * t257;
t181 = m(7) * t191 + t213 * mrSges(7,1) - t214 * mrSges(7,3) + t256 * t246 - t257 * t249;
t314 = t327 * g(1) - t331 * g(2);
t298 = -qJDD(1) * pkin(1) - t333 * pkin(7) - t314;
t319 = t326 * t353;
t310 = qJDD(1) * t330 - t319;
t261 = (-t309 - t351) * pkin(8) + (-t310 + t319) * pkin(2) + t298;
t284 = -g(3) * t326 + t330 * t299;
t264 = -pkin(2) * t332 + qJDD(2) * pkin(8) + t308 * t354 + t284;
t242 = t329 * t261 - t325 * t264;
t304 = qJDD(3) - t310;
t217 = (t305 * t318 - t276) * pkin(9) + (t305 * t306 + t304) * pkin(3) + t242;
t243 = t325 * t261 + t329 * t264;
t219 = -pkin(3) * t303 + pkin(9) * t275 - t285 * t318 + t243;
t197 = t328 * t217 - t324 * t219;
t300 = qJDD(4) + t304;
t193 = (t278 * t317 - t241) * qJ(5) + (t278 * t279 + t300) * pkin(4) + t197;
t198 = t324 * t217 + t328 * t219;
t195 = -pkin(4) * t277 + qJ(5) * t240 - t266 * t317 + t198;
t360 = -2 * qJD(5);
t189 = t323 * t193 + t358 * t195 + t256 * t360;
t225 = Ifges(7,1) * t257 + Ifges(7,4) * t317 + Ifges(7,5) * t256;
t226 = Ifges(6,1) * t257 - Ifges(6,4) * t256 + Ifges(6,5) * t317;
t231 = pkin(5) * t256 - qJ(6) * t257;
t316 = t317 ^ 2;
t184 = -pkin(5) * t316 + qJ(6) * t300 + 0.2e1 * qJD(6) * t317 - t231 * t256 + t189;
t346 = -mrSges(7,1) * t191 + mrSges(7,2) * t184;
t223 = Ifges(7,4) * t257 + Ifges(7,2) * t317 + Ifges(7,6) * t256;
t357 = -Ifges(6,5) * t257 + Ifges(6,6) * t256 - Ifges(6,3) * t317 - t223;
t168 = -mrSges(6,1) * t200 + mrSges(6,3) * t189 - pkin(5) * t181 + (t225 + t226) * t317 + (Ifges(6,6) - Ifges(7,6)) * t300 + t357 * t257 + (Ifges(6,4) - Ifges(7,5)) * t214 + (-Ifges(6,2) - Ifges(7,3)) * t213 + t346;
t344 = t358 * t193 - t323 * t195;
t188 = t257 * t360 + t344;
t224 = Ifges(6,4) * t257 - Ifges(6,2) * t256 + Ifges(6,6) * t317;
t186 = -t300 * pkin(5) - t316 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t231) * t257 - t344;
t221 = Ifges(7,5) * t257 + Ifges(7,6) * t317 + Ifges(7,3) * t256;
t343 = mrSges(7,2) * t186 - mrSges(7,3) * t191 + Ifges(7,1) * t214 + Ifges(7,4) * t300 + Ifges(7,5) * t213 + t317 * t221;
t169 = mrSges(6,2) * t200 - mrSges(6,3) * t188 + Ifges(6,1) * t214 - Ifges(6,4) * t213 + Ifges(6,5) * t300 - qJ(6) * t181 - t317 * t224 + t357 * t256 + t343;
t250 = Ifges(5,5) * t279 + Ifges(5,6) * t278 + Ifges(5,3) * t317;
t252 = Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t317;
t247 = -mrSges(6,2) * t317 - mrSges(6,3) * t256;
t248 = mrSges(6,1) * t317 - mrSges(6,3) * t257;
t341 = m(6) * t200 + t213 * mrSges(6,1) + t214 * mrSges(6,2) + t256 * t247 + t257 * t248 + t181;
t352 = m(7) * t184 + t300 * mrSges(7,3) + t317 * t249;
t232 = mrSges(7,1) * t256 - mrSges(7,3) * t257;
t356 = -mrSges(6,1) * t256 - mrSges(6,2) * t257 - t232;
t359 = -mrSges(6,3) - mrSges(7,2);
t172 = m(6) * t189 - t300 * mrSges(6,2) + t359 * t213 - t317 * t248 + t356 * t256 + t352;
t347 = -m(7) * t186 + t300 * mrSges(7,1) + t317 * t246;
t174 = m(6) * t188 + t300 * mrSges(6,1) + t359 * t214 + t317 * t247 + t356 * t257 + t347;
t348 = t358 * t172 - t174 * t323;
t155 = -mrSges(5,1) * t234 + mrSges(5,3) * t198 + Ifges(5,4) * t241 + Ifges(5,2) * t240 + Ifges(5,6) * t300 - pkin(4) * t341 + qJ(5) * t348 + t358 * t168 + t323 * t169 - t279 * t250 + t317 * t252;
t167 = t323 * t172 + t358 * t174;
t251 = Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t317;
t156 = mrSges(5,2) * t234 - mrSges(5,3) * t197 + Ifges(5,1) * t241 + Ifges(5,4) * t240 + Ifges(5,5) * t300 - qJ(5) * t167 - t323 * t168 + t358 * t169 + t278 * t250 - t317 * t251;
t268 = Ifges(4,5) * t306 + Ifges(4,6) * t305 + Ifges(4,3) * t318;
t270 = Ifges(4,1) * t306 + Ifges(4,4) * t305 + Ifges(4,5) * t318;
t265 = -mrSges(5,2) * t317 + mrSges(5,3) * t278;
t267 = mrSges(5,1) * t317 - mrSges(5,3) * t279;
t338 = m(5) * t234 - t240 * mrSges(5,1) + t241 * mrSges(5,2) - t278 * t265 + t279 * t267 + t341;
t258 = -mrSges(5,1) * t278 + mrSges(5,2) * t279;
t164 = m(5) * t197 + mrSges(5,1) * t300 - mrSges(5,3) * t241 - t258 * t279 + t265 * t317 + t167;
t165 = m(5) * t198 - mrSges(5,2) * t300 + mrSges(5,3) * t240 + t258 * t278 - t267 * t317 + t348;
t349 = -t164 * t324 + t328 * t165;
t144 = -mrSges(4,1) * t263 + mrSges(4,3) * t243 + Ifges(4,4) * t276 + Ifges(4,2) * t275 + Ifges(4,6) * t304 - pkin(3) * t338 + pkin(9) * t349 + t328 * t155 + t324 * t156 - t306 * t268 + t318 * t270;
t160 = t328 * t164 + t324 * t165;
t269 = Ifges(4,4) * t306 + Ifges(4,2) * t305 + Ifges(4,6) * t318;
t145 = mrSges(4,2) * t263 - mrSges(4,3) * t242 + Ifges(4,1) * t276 + Ifges(4,4) * t275 + Ifges(4,5) * t304 - pkin(9) * t160 - t155 * t324 + t156 * t328 + t268 * t305 - t269 * t318;
t280 = -mrSges(4,1) * t305 + mrSges(4,2) * t306;
t281 = -mrSges(4,2) * t318 + mrSges(4,3) * t305;
t158 = m(4) * t242 + mrSges(4,1) * t304 - mrSges(4,3) * t276 - t280 * t306 + t281 * t318 + t160;
t282 = mrSges(4,1) * t318 - mrSges(4,3) * t306;
t159 = m(4) * t243 - mrSges(4,2) * t304 + mrSges(4,3) * t275 + t280 * t305 - t282 * t318 + t349;
t154 = -t158 * t325 + t329 * t159;
t176 = -m(4) * t263 + t275 * mrSges(4,1) - t276 * mrSges(4,2) + t305 * t281 - t306 * t282 - t338;
t296 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t326 + Ifges(3,2) * t330) * qJD(1);
t297 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t326 + Ifges(3,4) * t330) * qJD(1);
t361 = mrSges(3,1) * t283 - mrSges(3,2) * t284 + Ifges(3,5) * t309 + Ifges(3,6) * t310 + Ifges(3,3) * qJDD(2) + pkin(2) * t176 + pkin(8) * t154 + t329 * t144 + t325 * t145 + (t296 * t326 - t297 * t330) * qJD(1);
t307 = (-mrSges(3,1) * t330 + mrSges(3,2) * t326) * qJD(1);
t312 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t355;
t152 = m(3) * t284 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t310 - qJD(2) * t312 + t307 * t354 + t154;
t313 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t354;
t175 = m(3) * t283 + qJDD(2) * mrSges(3,1) - t309 * mrSges(3,3) + qJD(2) * t313 - t307 * t355 + t176;
t350 = t330 * t152 - t175 * t326;
t153 = t158 * t329 + t159 * t325;
t295 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t326 + Ifges(3,6) * t330) * qJD(1);
t141 = mrSges(3,2) * t298 - mrSges(3,3) * t283 + Ifges(3,1) * t309 + Ifges(3,4) * t310 + Ifges(3,5) * qJDD(2) - pkin(8) * t153 - qJD(2) * t296 - t144 * t325 + t145 * t329 + t295 * t354;
t340 = mrSges(7,1) * t186 - mrSges(7,3) * t184 - Ifges(7,4) * t214 - Ifges(7,2) * t300 - Ifges(7,6) * t213 + t257 * t221 - t256 * t225;
t337 = mrSges(6,2) * t189 - t256 * t226 - qJ(6) * (-t213 * mrSges(7,2) - t256 * t232 + t352) - pkin(5) * (-t214 * mrSges(7,2) - t257 * t232 + t347) - mrSges(6,1) * t188 + Ifges(6,6) * t213 - Ifges(6,5) * t214 - t257 * t224 - Ifges(6,3) * t300 + t340;
t335 = -mrSges(5,1) * t197 + mrSges(5,2) * t198 - Ifges(5,5) * t241 - Ifges(5,6) * t240 - Ifges(5,3) * t300 - pkin(4) * t167 - t279 * t251 + t278 * t252 + t337;
t334 = mrSges(4,1) * t242 - mrSges(4,2) * t243 + Ifges(4,5) * t276 + Ifges(4,6) * t275 + Ifges(4,3) * t304 + pkin(3) * t160 + t306 * t269 - t305 * t270 - t335;
t143 = -mrSges(3,1) * t298 + mrSges(3,3) * t284 + Ifges(3,4) * t309 + Ifges(3,2) * t310 + Ifges(3,6) * qJDD(2) - pkin(2) * t153 + qJD(2) * t297 - t295 * t355 - t334;
t339 = -m(3) * t298 + t310 * mrSges(3,1) - mrSges(3,2) * t309 - t312 * t355 + t313 * t354 - t153;
t342 = mrSges(2,1) * t314 - mrSges(2,2) * t315 + Ifges(2,3) * qJDD(1) + pkin(1) * t339 + pkin(7) * t350 + t326 * t141 + t330 * t143;
t149 = m(2) * t314 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t333 + t339;
t148 = t152 * t326 + t175 * t330;
t146 = m(2) * t315 - mrSges(2,1) * t333 - qJDD(1) * mrSges(2,2) + t350;
t139 = mrSges(2,1) * g(3) + mrSges(2,3) * t315 + t333 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t148 - t361;
t138 = -mrSges(2,2) * g(3) - mrSges(2,3) * t314 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t333 - pkin(7) * t148 + t141 * t330 - t143 * t326;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t331 * t138 - t327 * t139 - pkin(6) * (t146 * t327 + t149 * t331), t138, t141, t145, t156, t169, -t223 * t256 + t343; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t327 * t138 + t331 * t139 + pkin(6) * (t146 * t331 - t149 * t327), t139, t143, t144, t155, t168, -t340; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t342, t342, t361, t334, -t335, -t337, Ifges(7,5) * t214 + Ifges(7,6) * t300 + Ifges(7,3) * t213 + t257 * t223 - t317 * t225 - t346;];
m_new  = t1;
