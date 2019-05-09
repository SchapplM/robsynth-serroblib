% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:47:06
% EndTime: 2019-05-05 17:47:33
% DurationCPUTime: 16.40s
% Computational Cost: add. (264113->362), mult. (648341->442), div. (0->0), fcn. (485291->10), ass. (0->145)
t334 = qJD(1) ^ 2;
t326 = sin(pkin(9));
t367 = qJD(1) * t326;
t328 = cos(pkin(9));
t377 = qJD(1) * t328;
t331 = sin(qJ(1));
t332 = cos(qJ(1));
t312 = -g(1) * t332 - g(2) * t331;
t305 = -pkin(1) * t334 + qJDD(1) * qJ(2) + t312;
t364 = qJD(1) * qJD(2);
t359 = -t328 * g(3) - 0.2e1 * t326 * t364;
t373 = pkin(2) * t328;
t263 = (-pkin(7) * qJDD(1) + t334 * t373 - t305) * t326 + t359;
t290 = -t326 * g(3) + (t305 + 0.2e1 * t364) * t328;
t362 = qJDD(1) * t328;
t322 = t328 ^ 2;
t370 = t322 * t334;
t270 = -pkin(2) * t370 + pkin(7) * t362 + t290;
t330 = sin(qJ(3));
t375 = cos(qJ(3));
t244 = t263 * t330 + t270 * t375;
t360 = t328 * t375;
t303 = -qJD(1) * t360 + t330 * t367;
t346 = t326 * t375 + t328 * t330;
t304 = t346 * qJD(1);
t281 = mrSges(4,1) * t303 + mrSges(4,2) * t304;
t363 = qJDD(1) * t326;
t365 = t304 * qJD(3);
t287 = -qJDD(1) * t360 + t330 * t363 + t365;
t298 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t304;
t280 = pkin(3) * t303 - qJ(4) * t304;
t333 = qJD(3) ^ 2;
t224 = -pkin(3) * t333 + qJDD(3) * qJ(4) - t280 * t303 + t244;
t321 = t326 ^ 2;
t311 = t331 * g(1) - g(2) * t332;
t353 = qJDD(2) - t311;
t286 = (-pkin(1) - t373) * qJDD(1) + (-qJ(2) + (-t321 - t322) * pkin(7)) * t334 + t353;
t366 = t303 * qJD(3);
t288 = qJDD(1) * t346 - t366;
t227 = (-t288 + t366) * qJ(4) + (t287 + t365) * pkin(3) + t286;
t325 = sin(pkin(10));
t327 = cos(pkin(10));
t296 = qJD(3) * t325 + t304 * t327;
t205 = -0.2e1 * qJD(4) * t296 - t325 * t224 + t227 * t327;
t269 = qJDD(3) * t325 + t288 * t327;
t295 = qJD(3) * t327 - t304 * t325;
t202 = (t295 * t303 - t269) * pkin(8) + (t295 * t296 + t287) * pkin(4) + t205;
t206 = 0.2e1 * qJD(4) * t295 + t224 * t327 + t227 * t325;
t267 = pkin(4) * t303 - pkin(8) * t296;
t268 = qJDD(3) * t327 - t288 * t325;
t294 = t295 ^ 2;
t204 = -pkin(4) * t294 + pkin(8) * t268 - t267 * t303 + t206;
t329 = sin(qJ(5));
t374 = cos(qJ(5));
t198 = t202 * t329 + t204 * t374;
t255 = t295 * t329 + t296 * t374;
t222 = qJD(5) * t255 - t268 * t374 + t269 * t329;
t301 = qJD(5) + t303;
t247 = mrSges(6,1) * t301 - mrSges(6,3) * t255;
t254 = -t295 * t374 + t296 * t329;
t285 = qJDD(5) + t287;
t237 = pkin(5) * t254 - qJ(6) * t255;
t300 = t301 ^ 2;
t193 = -pkin(5) * t300 + qJ(6) * t285 + 0.2e1 * qJD(6) * t301 - t237 * t254 + t198;
t248 = -mrSges(7,1) * t301 + mrSges(7,2) * t255;
t361 = m(7) * t193 + mrSges(7,3) * t285 + t248 * t301;
t238 = mrSges(7,1) * t254 - mrSges(7,3) * t255;
t368 = -mrSges(6,1) * t254 - mrSges(6,2) * t255 - t238;
t372 = -mrSges(6,3) - mrSges(7,2);
t179 = m(6) * t198 - t285 * mrSges(6,2) + t222 * t372 - t301 * t247 + t254 * t368 + t361;
t197 = t202 * t374 - t204 * t329;
t223 = -qJD(5) * t254 + t268 * t329 + t269 * t374;
t246 = -mrSges(6,2) * t301 - mrSges(6,3) * t254;
t195 = -pkin(5) * t285 - qJ(6) * t300 + t237 * t255 + qJDD(6) - t197;
t245 = -mrSges(7,2) * t254 + mrSges(7,3) * t301;
t354 = -m(7) * t195 + mrSges(7,1) * t285 + t245 * t301;
t181 = m(6) * t197 + t285 * mrSges(6,1) + t223 * t372 + t301 * t246 + t255 * t368 + t354;
t174 = t179 * t329 + t181 * t374;
t257 = -mrSges(5,1) * t295 + mrSges(5,2) * t296;
t265 = -mrSges(5,2) * t303 + mrSges(5,3) * t295;
t172 = m(5) * t205 + mrSges(5,1) * t287 - mrSges(5,3) * t269 - t257 * t296 + t265 * t303 + t174;
t266 = mrSges(5,1) * t303 - mrSges(5,3) * t296;
t355 = t179 * t374 - t181 * t329;
t173 = m(5) * t206 - mrSges(5,2) * t287 + mrSges(5,3) * t268 + t257 * t295 - t266 * t303 + t355;
t356 = -t172 * t325 + t173 * t327;
t165 = m(4) * t244 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t287 - qJD(3) * t298 - t281 * t303 + t356;
t243 = t263 * t375 - t270 * t330;
t297 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t303;
t221 = -qJDD(3) * pkin(3) - t333 * qJ(4) + t280 * t304 + qJDD(4) - t243;
t207 = -t268 * pkin(4) - t294 * pkin(8) + t267 * t296 + t221;
t200 = -0.2e1 * qJD(6) * t255 + (t254 * t301 - t223) * qJ(6) + (t255 * t301 + t222) * pkin(5) + t207;
t190 = m(7) * t200 + mrSges(7,1) * t222 - mrSges(7,3) * t223 + t245 * t254 - t248 * t255;
t340 = m(6) * t207 + mrSges(6,1) * t222 + mrSges(6,2) * t223 + t246 * t254 + t247 * t255 + t190;
t336 = -m(5) * t221 + mrSges(5,1) * t268 - mrSges(5,2) * t269 + t265 * t295 - t266 * t296 - t340;
t183 = m(4) * t243 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t288 + qJD(3) * t297 - t281 * t304 + t336;
t160 = t165 * t330 + t183 * t375;
t289 = -t326 * t305 + t359;
t232 = Ifges(7,1) * t255 + Ifges(7,4) * t301 + Ifges(7,5) * t254;
t233 = Ifges(6,1) * t255 - Ifges(6,4) * t254 + Ifges(6,5) * t301;
t352 = -mrSges(7,1) * t200 + mrSges(7,2) * t193;
t230 = Ifges(7,4) * t255 + Ifges(7,2) * t301 + Ifges(7,6) * t254;
t369 = -Ifges(6,5) * t255 + Ifges(6,6) * t254 - Ifges(6,3) * t301 - t230;
t175 = -mrSges(6,1) * t207 + mrSges(6,3) * t198 - pkin(5) * t190 + (t232 + t233) * t301 + (Ifges(6,6) - Ifges(7,6)) * t285 + t369 * t255 + (Ifges(6,4) - Ifges(7,5)) * t223 + (-Ifges(6,2) - Ifges(7,3)) * t222 + t352;
t231 = Ifges(6,4) * t255 - Ifges(6,2) * t254 + Ifges(6,6) * t301;
t228 = Ifges(7,5) * t255 + Ifges(7,6) * t301 + Ifges(7,3) * t254;
t345 = mrSges(7,2) * t195 - mrSges(7,3) * t200 + Ifges(7,1) * t223 + Ifges(7,4) * t285 + Ifges(7,5) * t222 + t228 * t301;
t176 = mrSges(6,2) * t207 - mrSges(6,3) * t197 + Ifges(6,1) * t223 - Ifges(6,4) * t222 + Ifges(6,5) * t285 - qJ(6) * t190 - t301 * t231 + t254 * t369 + t345;
t249 = Ifges(5,5) * t296 + Ifges(5,6) * t295 + Ifges(5,3) * t303;
t251 = Ifges(5,1) * t296 + Ifges(5,4) * t295 + Ifges(5,5) * t303;
t154 = -mrSges(5,1) * t221 + mrSges(5,3) * t206 + Ifges(5,4) * t269 + Ifges(5,2) * t268 + Ifges(5,6) * t287 - pkin(4) * t340 + pkin(8) * t355 + t175 * t374 + t329 * t176 - t296 * t249 + t303 * t251;
t250 = Ifges(5,4) * t296 + Ifges(5,2) * t295 + Ifges(5,6) * t303;
t156 = mrSges(5,2) * t221 - mrSges(5,3) * t205 + Ifges(5,1) * t269 + Ifges(5,4) * t268 + Ifges(5,5) * t287 - pkin(8) * t174 - t175 * t329 + t176 * t374 + t249 * t295 - t250 * t303;
t272 = Ifges(4,4) * t304 - Ifges(4,2) * t303 + Ifges(4,6) * qJD(3);
t273 = Ifges(4,1) * t304 - Ifges(4,4) * t303 + Ifges(4,5) * qJD(3);
t341 = -mrSges(4,1) * t243 + mrSges(4,2) * t244 - Ifges(4,5) * t288 + Ifges(4,6) * t287 - Ifges(4,3) * qJDD(3) - pkin(3) * t336 - qJ(4) * t356 - t154 * t327 - t156 * t325 - t272 * t304 - t303 * t273;
t350 = Ifges(3,4) * t326 + Ifges(3,2) * t328;
t351 = Ifges(3,1) * t326 + Ifges(3,4) * t328;
t376 = -mrSges(3,1) * t289 + mrSges(3,2) * t290 - pkin(2) * t160 - (t350 * t367 - t351 * t377) * qJD(1) + t341;
t371 = mrSges(3,2) * t326;
t167 = t172 * t327 + t173 * t325;
t347 = mrSges(3,3) * qJDD(1) + t334 * (-mrSges(3,1) * t328 + t371);
t158 = m(3) * t289 - t326 * t347 + t160;
t357 = t165 * t375 - t330 * t183;
t159 = m(3) * t290 + t328 * t347 + t357;
t358 = -t158 * t326 + t159 * t328;
t349 = Ifges(3,5) * t326 + Ifges(3,6) * t328;
t271 = Ifges(4,5) * t304 - Ifges(4,6) * t303 + Ifges(4,3) * qJD(3);
t148 = mrSges(4,2) * t286 - mrSges(4,3) * t243 + Ifges(4,1) * t288 - Ifges(4,4) * t287 + Ifges(4,5) * qJDD(3) - qJ(4) * t167 - qJD(3) * t272 - t154 * t325 + t156 * t327 - t271 * t303;
t343 = mrSges(7,1) * t195 - mrSges(7,3) * t193 - Ifges(7,4) * t223 - Ifges(7,2) * t285 - Ifges(7,6) * t222 + t255 * t228 - t232 * t254;
t337 = mrSges(6,2) * t198 - t254 * t233 - qJ(6) * (-t222 * mrSges(7,2) - t254 * t238 + t361) - pkin(5) * (-t223 * mrSges(7,2) - t255 * t238 + t354) - mrSges(6,1) * t197 - t255 * t231 + Ifges(6,6) * t222 - Ifges(6,5) * t223 - Ifges(6,3) * t285 + t343;
t335 = mrSges(5,1) * t205 - mrSges(5,2) * t206 + Ifges(5,5) * t269 + Ifges(5,6) * t268 + pkin(4) * t174 + t296 * t250 - t295 * t251 - t337;
t152 = -t304 * t271 - mrSges(4,1) * t286 + Ifges(4,4) * t288 + qJD(3) * t273 + mrSges(4,3) * t244 + Ifges(4,6) * qJDD(3) - pkin(3) * t167 - t335 + (-Ifges(5,3) - Ifges(4,2)) * t287;
t302 = -qJDD(1) * pkin(1) - t334 * qJ(2) + t353;
t307 = t349 * qJD(1);
t342 = m(4) * t286 + mrSges(4,1) * t287 + t288 * mrSges(4,2) + t297 * t303 + t304 * t298 + t167;
t144 = -mrSges(3,1) * t302 + mrSges(3,3) * t290 - pkin(2) * t342 + pkin(7) * t357 + qJDD(1) * t350 + t330 * t148 + t152 * t375 - t307 * t367;
t147 = mrSges(3,2) * t302 - mrSges(3,3) * t289 - pkin(7) * t160 + qJDD(1) * t351 + t148 * t375 - t330 * t152 + t307 * t377;
t339 = -m(3) * t302 + mrSges(3,1) * t362 - t342 + (t321 * t334 + t370) * mrSges(3,3);
t344 = -mrSges(2,2) * t312 + qJ(2) * t358 + t328 * t144 + t326 * t147 + pkin(1) * (-mrSges(3,2) * t363 + t339) + mrSges(2,1) * t311 + Ifges(2,3) * qJDD(1);
t161 = m(2) * t311 + (mrSges(2,1) - t371) * qJDD(1) - t334 * mrSges(2,2) + t339;
t151 = t158 * t328 + t159 * t326;
t149 = m(2) * t312 - mrSges(2,1) * t334 - qJDD(1) * mrSges(2,2) + t358;
t145 = mrSges(2,3) * t312 + mrSges(2,1) * g(3) + t334 * Ifges(2,5) - pkin(1) * t151 + (Ifges(2,6) - t349) * qJDD(1) + t376;
t142 = -mrSges(2,2) * g(3) - mrSges(2,3) * t311 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t334 - qJ(2) * t151 - t144 * t326 + t147 * t328;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t332 * t142 - t331 * t145 - pkin(6) * (t149 * t331 + t161 * t332), t142, t147, t148, t156, t176, -t230 * t254 + t345; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t331 * t142 + t332 * t145 + pkin(6) * (t149 * t332 - t161 * t331), t145, t144, t152, t154, t175, -t343; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t344, t344, qJDD(1) * t349 - t376, -t341, Ifges(5,3) * t287 + t335, -t337, Ifges(7,5) * t223 + Ifges(7,6) * t285 + Ifges(7,3) * t222 + t255 * t230 - t301 * t232 - t352;];
m_new  = t1;
