% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:51:11
% EndTime: 2019-05-07 07:52:00
% DurationCPUTime: 22.91s
% Computational Cost: add. (406010->380), mult. (849814->464), div. (0->0), fcn. (611493->10), ass. (0->141)
t329 = sin(qJ(1));
t332 = cos(qJ(1));
t316 = -t332 * g(1) - t329 * g(2);
t334 = qJD(1) ^ 2;
t299 = -t334 * pkin(1) + qJDD(1) * pkin(7) + t316;
t328 = sin(qJ(2));
t331 = cos(qJ(2));
t285 = -t331 * g(3) - t328 * t299;
t309 = (-pkin(2) * t331 - pkin(8) * t328) * qJD(1);
t333 = qJD(2) ^ 2;
t355 = qJD(1) * t328;
t262 = -qJDD(2) * pkin(2) - t333 * pkin(8) + t309 * t355 - t285;
t327 = sin(qJ(3));
t330 = cos(qJ(3));
t307 = t327 * qJD(2) + t330 * t355;
t353 = qJD(1) * qJD(2);
t351 = t331 * t353;
t310 = t328 * qJDD(1) + t351;
t276 = -t307 * qJD(3) + t330 * qJDD(2) - t327 * t310;
t354 = t331 * qJD(1);
t319 = qJD(3) - t354;
t283 = t319 * pkin(3) - t307 * qJ(4);
t306 = t330 * qJD(2) - t327 * t355;
t304 = t306 ^ 2;
t234 = -t276 * pkin(3) - t304 * qJ(4) + t307 * t283 + qJDD(4) + t262;
t277 = t306 * qJD(3) + t327 * qJDD(2) + t330 * t310;
t324 = sin(pkin(10));
t325 = cos(pkin(10));
t251 = t325 * t276 - t324 * t277;
t280 = t324 * t306 + t325 * t307;
t266 = t319 * pkin(4) - t280 * pkin(9);
t279 = t325 * t306 - t324 * t307;
t278 = t279 ^ 2;
t200 = -t251 * pkin(4) - t278 * pkin(9) + t280 * t266 + t234;
t252 = t324 * t276 + t325 * t277;
t326 = sin(qJ(5));
t359 = cos(qJ(5));
t256 = t326 * t279 + t359 * t280;
t211 = t256 * qJD(5) - t359 * t251 + t326 * t252;
t255 = -t359 * t279 + t326 * t280;
t212 = -t255 * qJD(5) + t326 * t251 + t359 * t252;
t318 = qJD(5) + t319;
t191 = (t255 * t318 - t212) * qJ(6) + (t256 * t318 + t211) * pkin(5) + t200 - 0.2e1 * qJD(6) * t256;
t240 = -t255 * mrSges(7,2) + t318 * mrSges(7,3);
t243 = -t318 * mrSges(7,1) + t256 * mrSges(7,2);
t181 = m(7) * t191 + t211 * mrSges(7,1) - t212 * mrSges(7,3) + t255 * t240 - t256 * t243;
t315 = t329 * g(1) - t332 * g(2);
t298 = -qJDD(1) * pkin(1) - t334 * pkin(7) - t315;
t320 = t328 * t353;
t311 = t331 * qJDD(1) - t320;
t260 = (-t310 - t351) * pkin(8) + (-t311 + t320) * pkin(2) + t298;
t286 = -t328 * g(3) + t331 * t299;
t263 = -t333 * pkin(2) + qJDD(2) * pkin(8) + t309 * t354 + t286;
t236 = t330 * t260 - t327 * t263;
t305 = qJDD(3) - t311;
t217 = (t306 * t319 - t277) * qJ(4) + (t306 * t307 + t305) * pkin(3) + t236;
t237 = t327 * t260 + t330 * t263;
t219 = -t304 * pkin(3) + t276 * qJ(4) - t319 * t283 + t237;
t197 = -0.2e1 * qJD(4) * t280 + t325 * t217 - t324 * t219;
t193 = (t279 * t319 - t252) * pkin(9) + (t279 * t280 + t305) * pkin(4) + t197;
t198 = 0.2e1 * qJD(4) * t279 + t324 * t217 + t325 * t219;
t195 = -t278 * pkin(4) + t251 * pkin(9) - t319 * t266 + t198;
t189 = t326 * t193 + t359 * t195;
t225 = Ifges(7,1) * t256 + Ifges(7,4) * t318 + Ifges(7,5) * t255;
t226 = Ifges(6,1) * t256 - Ifges(6,4) * t255 + Ifges(6,5) * t318;
t301 = qJDD(5) + t305;
t231 = t255 * pkin(5) - t256 * qJ(6);
t317 = t318 ^ 2;
t184 = -t317 * pkin(5) + t301 * qJ(6) + 0.2e1 * qJD(6) * t318 - t255 * t231 + t189;
t346 = -mrSges(7,1) * t191 + mrSges(7,2) * t184;
t223 = Ifges(7,4) * t256 + Ifges(7,2) * t318 + Ifges(7,6) * t255;
t357 = -Ifges(6,5) * t256 + Ifges(6,6) * t255 - Ifges(6,3) * t318 - t223;
t168 = -mrSges(6,1) * t200 + mrSges(6,3) * t189 - pkin(5) * t181 + (t225 + t226) * t318 + (Ifges(6,6) - Ifges(7,6)) * t301 + t357 * t256 + (Ifges(6,4) - Ifges(7,5)) * t212 + (-Ifges(6,2) - Ifges(7,3)) * t211 + t346;
t188 = t359 * t193 - t326 * t195;
t224 = Ifges(6,4) * t256 - Ifges(6,2) * t255 + Ifges(6,6) * t318;
t186 = -t301 * pkin(5) - t317 * qJ(6) + t256 * t231 + qJDD(6) - t188;
t221 = Ifges(7,5) * t256 + Ifges(7,6) * t318 + Ifges(7,3) * t255;
t344 = mrSges(7,2) * t186 - mrSges(7,3) * t191 + Ifges(7,1) * t212 + Ifges(7,4) * t301 + Ifges(7,5) * t211 + t318 * t221;
t169 = mrSges(6,2) * t200 - mrSges(6,3) * t188 + Ifges(6,1) * t212 - Ifges(6,4) * t211 + Ifges(6,5) * t301 - qJ(6) * t181 - t318 * t224 + t357 * t255 + t344;
t248 = Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t319;
t250 = Ifges(5,1) * t280 + Ifges(5,4) * t279 + Ifges(5,5) * t319;
t241 = -t318 * mrSges(6,2) - t255 * mrSges(6,3);
t242 = t318 * mrSges(6,1) - t256 * mrSges(6,3);
t342 = m(6) * t200 + t211 * mrSges(6,1) + t212 * mrSges(6,2) + t255 * t241 + t256 * t242 + t181;
t352 = m(7) * t184 + t301 * mrSges(7,3) + t318 * t243;
t232 = t255 * mrSges(7,1) - t256 * mrSges(7,3);
t356 = -t255 * mrSges(6,1) - t256 * mrSges(6,2) - t232;
t358 = -mrSges(6,3) - mrSges(7,2);
t172 = m(6) * t189 - t301 * mrSges(6,2) + t358 * t211 - t318 * t242 + t356 * t255 + t352;
t347 = -m(7) * t186 + t301 * mrSges(7,1) + t318 * t240;
t174 = m(6) * t188 + t301 * mrSges(6,1) + t358 * t212 + t318 * t241 + t356 * t256 + t347;
t348 = t359 * t172 - t326 * t174;
t155 = -mrSges(5,1) * t234 + mrSges(5,3) * t198 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * t305 - pkin(4) * t342 + pkin(9) * t348 + t359 * t168 + t326 * t169 - t280 * t248 + t319 * t250;
t167 = t326 * t172 + t359 * t174;
t249 = Ifges(5,4) * t280 + Ifges(5,2) * t279 + Ifges(5,6) * t319;
t156 = mrSges(5,2) * t234 - mrSges(5,3) * t197 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * t305 - pkin(9) * t167 - t326 * t168 + t359 * t169 + t279 * t248 - t319 * t249;
t267 = Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t319;
t269 = Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t319;
t264 = -t319 * mrSges(5,2) + t279 * mrSges(5,3);
t265 = t319 * mrSges(5,1) - t280 * mrSges(5,3);
t339 = m(5) * t234 - t251 * mrSges(5,1) + t252 * mrSges(5,2) - t279 * t264 + t280 * t265 + t342;
t257 = -t279 * mrSges(5,1) + t280 * mrSges(5,2);
t164 = m(5) * t197 + t305 * mrSges(5,1) - t252 * mrSges(5,3) - t280 * t257 + t319 * t264 + t167;
t165 = m(5) * t198 - t305 * mrSges(5,2) + t251 * mrSges(5,3) + t279 * t257 - t319 * t265 + t348;
t349 = -t324 * t164 + t325 * t165;
t144 = -mrSges(4,1) * t262 + mrSges(4,3) * t237 + Ifges(4,4) * t277 + Ifges(4,2) * t276 + Ifges(4,6) * t305 - pkin(3) * t339 + qJ(4) * t349 + t325 * t155 + t324 * t156 - t307 * t267 + t319 * t269;
t160 = t325 * t164 + t324 * t165;
t268 = Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t319;
t145 = mrSges(4,2) * t262 - mrSges(4,3) * t236 + Ifges(4,1) * t277 + Ifges(4,4) * t276 + Ifges(4,5) * t305 - qJ(4) * t160 - t324 * t155 + t325 * t156 + t306 * t267 - t319 * t268;
t281 = -t306 * mrSges(4,1) + t307 * mrSges(4,2);
t282 = -t319 * mrSges(4,2) + t306 * mrSges(4,3);
t158 = m(4) * t236 + t305 * mrSges(4,1) - t277 * mrSges(4,3) - t307 * t281 + t319 * t282 + t160;
t284 = t319 * mrSges(4,1) - t307 * mrSges(4,3);
t159 = m(4) * t237 - t305 * mrSges(4,2) + t276 * mrSges(4,3) + t306 * t281 - t319 * t284 + t349;
t154 = -t327 * t158 + t330 * t159;
t176 = -m(4) * t262 + t276 * mrSges(4,1) - t277 * mrSges(4,2) + t306 * t282 - t307 * t284 - t339;
t296 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t328 + Ifges(3,2) * t331) * qJD(1);
t297 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t328 + Ifges(3,4) * t331) * qJD(1);
t360 = mrSges(3,1) * t285 - mrSges(3,2) * t286 + Ifges(3,5) * t310 + Ifges(3,6) * t311 + Ifges(3,3) * qJDD(2) + pkin(2) * t176 + pkin(8) * t154 + t330 * t144 + t327 * t145 + (t328 * t296 - t331 * t297) * qJD(1);
t308 = (-mrSges(3,1) * t331 + mrSges(3,2) * t328) * qJD(1);
t313 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t355;
t152 = m(3) * t286 - qJDD(2) * mrSges(3,2) + t311 * mrSges(3,3) - qJD(2) * t313 + t308 * t354 + t154;
t314 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t354;
t175 = m(3) * t285 + qJDD(2) * mrSges(3,1) - t310 * mrSges(3,3) + qJD(2) * t314 - t308 * t355 + t176;
t350 = t331 * t152 - t328 * t175;
t153 = t330 * t158 + t327 * t159;
t295 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t328 + Ifges(3,6) * t331) * qJD(1);
t141 = mrSges(3,2) * t298 - mrSges(3,3) * t285 + Ifges(3,1) * t310 + Ifges(3,4) * t311 + Ifges(3,5) * qJDD(2) - pkin(8) * t153 - qJD(2) * t296 - t327 * t144 + t330 * t145 + t295 * t354;
t341 = mrSges(7,1) * t186 - mrSges(7,3) * t184 - Ifges(7,4) * t212 - Ifges(7,2) * t301 - Ifges(7,6) * t211 + t256 * t221 - t255 * t225;
t338 = mrSges(6,2) * t189 - t255 * t226 - qJ(6) * (-t211 * mrSges(7,2) - t255 * t232 + t352) - pkin(5) * (-t212 * mrSges(7,2) - t256 * t232 + t347) - mrSges(6,1) * t188 + Ifges(6,6) * t211 - Ifges(6,5) * t212 - t256 * t224 - Ifges(6,3) * t301 + t341;
t336 = -mrSges(5,1) * t197 + mrSges(5,2) * t198 - Ifges(5,5) * t252 - Ifges(5,6) * t251 - Ifges(5,3) * t305 - pkin(4) * t167 - t280 * t249 + t279 * t250 + t338;
t335 = mrSges(4,1) * t236 - mrSges(4,2) * t237 + Ifges(4,5) * t277 + Ifges(4,6) * t276 + Ifges(4,3) * t305 + pkin(3) * t160 + t307 * t268 - t306 * t269 - t336;
t143 = -mrSges(3,1) * t298 + mrSges(3,3) * t286 + Ifges(3,4) * t310 + Ifges(3,2) * t311 + Ifges(3,6) * qJDD(2) - pkin(2) * t153 + qJD(2) * t297 - t295 * t355 - t335;
t340 = -m(3) * t298 + t311 * mrSges(3,1) - t310 * mrSges(3,2) - t313 * t355 + t314 * t354 - t153;
t343 = mrSges(2,1) * t315 - mrSges(2,2) * t316 + Ifges(2,3) * qJDD(1) + pkin(1) * t340 + pkin(7) * t350 + t328 * t141 + t331 * t143;
t149 = m(2) * t315 + qJDD(1) * mrSges(2,1) - t334 * mrSges(2,2) + t340;
t148 = t328 * t152 + t331 * t175;
t146 = m(2) * t316 - t334 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t350;
t139 = mrSges(2,1) * g(3) + mrSges(2,3) * t316 + t334 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t148 - t360;
t138 = -mrSges(2,2) * g(3) - mrSges(2,3) * t315 + Ifges(2,5) * qJDD(1) - t334 * Ifges(2,6) - pkin(7) * t148 + t331 * t141 - t328 * t143;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t332 * t138 - t329 * t139 - pkin(6) * (t329 * t146 + t332 * t149), t138, t141, t145, t156, t169, -t255 * t223 + t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t329 * t138 + t332 * t139 + pkin(6) * (t332 * t146 - t329 * t149), t139, t143, t144, t155, t168, -t341; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t343, t343, t360, t335, -t336, -t338, Ifges(7,5) * t212 + Ifges(7,6) * t301 + Ifges(7,3) * t211 + t256 * t223 - t318 * t225 - t346;];
m_new  = t1;
