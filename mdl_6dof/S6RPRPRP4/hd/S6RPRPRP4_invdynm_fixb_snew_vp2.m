% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:42:37
% EndTime: 2019-05-05 17:42:47
% DurationCPUTime: 4.59s
% Computational Cost: add. (56499->345), mult. (109952->399), div. (0->0), fcn. (57693->8), ass. (0->127)
t349 = -2 * qJD(4);
t301 = sin(qJ(1));
t304 = cos(qJ(1));
t280 = t301 * g(1) - g(2) * t304;
t265 = qJDD(1) * pkin(1) + t280;
t281 = -g(1) * t304 - g(2) * t301;
t306 = qJD(1) ^ 2;
t269 = -pkin(1) * t306 + t281;
t297 = sin(pkin(9));
t298 = cos(pkin(9));
t226 = t297 * t265 + t298 * t269;
t205 = -pkin(2) * t306 + qJDD(1) * pkin(7) + t226;
t300 = sin(qJ(3));
t202 = t300 * t205;
t296 = -g(3) + qJDD(2);
t303 = cos(qJ(3));
t340 = t303 * t296;
t199 = -t202 + t340;
t200 = t303 * t205 + t300 * t296;
t245 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t300 + Ifges(4,4) * t303) * qJD(1);
t267 = (mrSges(5,2) * t303 - mrSges(5,3) * t300) * qJD(1);
t333 = qJD(1) * qJD(3);
t330 = t303 * t333;
t270 = qJDD(1) * t300 + t330;
t329 = t300 * t333;
t271 = qJDD(1) * t303 - t329;
t334 = qJD(1) * t303;
t277 = -mrSges(5,1) * t334 - qJD(3) * mrSges(5,3);
t266 = (-pkin(3) * t303 - qJ(4) * t300) * qJD(1);
t305 = qJD(3) ^ 2;
t194 = t305 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t349 - t266 * t334 - t200;
t335 = qJD(1) * t300;
t278 = mrSges(5,1) * t335 + qJD(3) * mrSges(5,2);
t279 = pkin(4) * t335 - qJD(3) * pkin(8);
t295 = t303 ^ 2;
t345 = pkin(8) * t306;
t191 = t271 * pkin(4) + qJD(3) * t279 - t295 * t345 - t194;
t299 = sin(qJ(5));
t302 = cos(qJ(5));
t264 = qJD(3) * t302 - t299 * t334;
t222 = qJD(5) * t264 + qJDD(3) * t299 + t302 * t271;
t263 = qJD(3) * t299 + t302 * t334;
t223 = -qJD(5) * t263 + qJDD(3) * t302 - t271 * t299;
t284 = qJD(5) + t335;
t183 = -0.2e1 * qJD(6) * t264 + (t263 * t284 - t223) * qJ(6) + (t264 * t284 + t222) * pkin(5) + t191;
t233 = -mrSges(7,2) * t263 + mrSges(7,3) * t284;
t235 = -mrSges(7,1) * t284 + mrSges(7,2) * t264;
t176 = m(7) * t183 + t222 * mrSges(7,1) - t223 * mrSges(7,3) + t263 * t233 - t264 * t235;
t232 = -mrSges(6,2) * t284 - mrSges(6,3) * t263;
t234 = mrSges(6,1) * t284 - mrSges(6,3) * t264;
t312 = m(6) * t191 + t222 * mrSges(6,1) + t223 * mrSges(6,2) + t263 * t232 + t264 * t234 + t176;
t311 = -m(5) * t194 + qJDD(3) * mrSges(5,3) + qJD(3) * t278 + t267 * t334 + t312;
t225 = t298 * t265 - t297 * t269;
t323 = -qJDD(1) * pkin(2) - t225;
t313 = pkin(3) * t329 + t335 * t349 + (-t270 - t330) * qJ(4) + t323;
t346 = -pkin(3) - pkin(8);
t188 = -t279 * t335 + (-pkin(4) * t295 - pkin(7)) * t306 + t346 * t271 + t313;
t324 = -t305 * qJ(4) + t266 * t335 + qJDD(4) + t202;
t192 = t270 * pkin(4) + t346 * qJDD(3) + (-pkin(4) * t333 - t300 * t345 - t296) * t303 + t324;
t186 = t302 * t188 + t299 * t192;
t210 = Ifges(7,1) * t264 + Ifges(7,4) * t284 + Ifges(7,5) * t263;
t211 = Ifges(6,1) * t264 - Ifges(6,4) * t263 + Ifges(6,5) * t284;
t262 = qJDD(5) + t270;
t229 = pkin(5) * t263 - qJ(6) * t264;
t282 = t284 ^ 2;
t179 = -pkin(5) * t282 + qJ(6) * t262 + 0.2e1 * qJD(6) * t284 - t229 * t263 + t186;
t325 = -mrSges(7,1) * t183 + mrSges(7,2) * t179;
t208 = Ifges(7,4) * t264 + Ifges(7,2) * t284 + Ifges(7,6) * t263;
t339 = -Ifges(6,5) * t264 + Ifges(6,6) * t263 - Ifges(6,3) * t284 - t208;
t158 = -mrSges(6,1) * t191 + mrSges(6,3) * t186 - pkin(5) * t176 + (t210 + t211) * t284 + t339 * t264 + (Ifges(6,6) - Ifges(7,6)) * t262 + (Ifges(6,4) - Ifges(7,5)) * t223 + (-Ifges(6,2) - Ifges(7,3)) * t222 + t325;
t185 = -t188 * t299 + t192 * t302;
t209 = Ifges(6,4) * t264 - Ifges(6,2) * t263 + Ifges(6,6) * t284;
t181 = -pkin(5) * t262 - qJ(6) * t282 + t229 * t264 + qJDD(6) - t185;
t206 = Ifges(7,5) * t264 + Ifges(7,6) * t284 + Ifges(7,3) * t263;
t320 = mrSges(7,2) * t181 - mrSges(7,3) * t183 + Ifges(7,1) * t223 + Ifges(7,4) * t262 + Ifges(7,5) * t222 + t284 * t206;
t160 = mrSges(6,2) * t191 - mrSges(6,3) * t185 + Ifges(6,1) * t223 - Ifges(6,4) * t222 + Ifges(6,5) * t262 - qJ(6) * t176 - t284 * t209 + t339 * t263 + t320;
t331 = m(7) * t179 + t262 * mrSges(7,3) + t284 * t235;
t230 = mrSges(7,1) * t263 - mrSges(7,3) * t264;
t338 = -mrSges(6,1) * t263 - mrSges(6,2) * t264 - t230;
t343 = -mrSges(6,3) - mrSges(7,2);
t167 = m(6) * t186 - t262 * mrSges(6,2) + t343 * t222 - t284 * t234 + t338 * t263 + t331;
t326 = -m(7) * t181 + t262 * mrSges(7,1) + t284 * t233;
t169 = m(6) * t185 + t262 * mrSges(6,1) + t343 * t223 + t284 * t232 + t338 * t264 + t326;
t161 = t299 * t167 + t302 * t169;
t196 = -qJDD(3) * pkin(3) + t324 - t340;
t247 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t300 - Ifges(5,6) * t303) * qJD(1);
t316 = -mrSges(5,2) * t196 + mrSges(5,3) * t194 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t270 + Ifges(5,5) * t271 + pkin(8) * t161 + t299 * t158 - t302 * t160 - t247 * t334;
t319 = -m(5) * t196 - t270 * mrSges(5,1) - t161;
t246 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t300 - Ifges(5,3) * t303) * qJD(1);
t336 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t300 + Ifges(4,2) * t303) * qJD(1) - t246;
t348 = (-t303 * t245 + t300 * t336) * qJD(1) + mrSges(4,1) * t199 - mrSges(4,2) * t200 + Ifges(4,5) * t270 + Ifges(4,6) * t271 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t277 - t267 * t335 + t319) + qJ(4) * (mrSges(5,1) * t271 + t311) - t316;
t344 = t306 * pkin(7);
t342 = Ifges(4,4) + Ifges(5,6);
t268 = (-mrSges(4,1) * t303 + mrSges(4,2) * t300) * qJD(1);
t276 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t334;
t155 = m(4) * t199 - t270 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t276 - t277) * qJD(3) + (-t267 - t268) * t335 + t319;
t275 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t335;
t165 = -qJDD(3) * mrSges(4,2) - qJD(3) * t275 + m(4) * t200 + t311 + t268 * t334 + (mrSges(4,3) + mrSges(5,1)) * t271;
t327 = -t155 * t300 + t303 * t165;
t148 = m(3) * t226 - mrSges(3,1) * t306 - qJDD(1) * mrSges(3,2) + t327;
t204 = t323 - t344;
t162 = t302 * t167 - t299 * t169;
t193 = -t271 * pkin(3) + t313 - t344;
t322 = -m(5) * t193 - t271 * mrSges(5,2) + t278 * t335 - t162;
t310 = -m(4) * t204 + t276 * t334 + t271 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t270 + (-t275 * t300 - t277 * t303) * qJD(1) + t322;
t152 = m(3) * t225 + qJDD(1) * mrSges(3,1) - t306 * mrSges(3,2) + t310;
t145 = t297 * t148 + t298 * t152;
t150 = t303 * t155 + t300 * t165;
t248 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t300 - Ifges(5,5) * t303) * qJD(1);
t337 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t300 + Ifges(4,6) * t303) * qJD(1) + t248;
t328 = t298 * t148 - t152 * t297;
t156 = -t270 * mrSges(5,3) + t277 * t334 - t322;
t314 = -mrSges(5,1) * t194 + mrSges(5,2) * t193 + pkin(4) * t312 - pkin(8) * t162 - t302 * t158 - t299 * t160;
t139 = -mrSges(4,1) * t204 + mrSges(4,3) * t200 - pkin(3) * t156 + (Ifges(4,2) + Ifges(5,3)) * t271 + t342 * t270 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t245 - t247) * qJD(3) - t337 * t335 + t314;
t317 = mrSges(7,1) * t181 - mrSges(7,3) * t179 - Ifges(7,4) * t223 - Ifges(7,2) * t262 - Ifges(7,6) * t222 + t264 * t206 - t263 * t210;
t309 = mrSges(6,2) * t186 - t263 * t211 - qJ(6) * (-t222 * mrSges(7,2) - t263 * t230 + t331) - pkin(5) * (-t223 * mrSges(7,2) - t264 * t230 + t326) - mrSges(6,1) * t185 - t264 * t209 + Ifges(6,6) * t222 - Ifges(6,5) * t223 - Ifges(6,3) * t262 + t317;
t308 = -mrSges(5,1) * t196 + mrSges(5,3) * t193 - pkin(4) * t161 + t309;
t141 = t342 * t271 + (Ifges(4,1) + Ifges(5,2)) * t270 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t336 * qJD(3) - t308 + mrSges(4,2) * t204 - mrSges(4,3) * t199 - qJ(4) * t156 + t337 * t334;
t318 = mrSges(3,1) * t225 - mrSges(3,2) * t226 + Ifges(3,3) * qJDD(1) + pkin(2) * t310 + pkin(7) * t327 + t303 * t139 + t300 * t141;
t315 = mrSges(2,1) * t280 - mrSges(2,2) * t281 + Ifges(2,3) * qJDD(1) + pkin(1) * t145 + t318;
t143 = m(2) * t281 - mrSges(2,1) * t306 - qJDD(1) * mrSges(2,2) + t328;
t142 = m(2) * t280 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t306 + t145;
t137 = -mrSges(3,1) * t296 + mrSges(3,3) * t226 + t306 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t150 - t348;
t136 = mrSges(3,2) * t296 - mrSges(3,3) * t225 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t306 - pkin(7) * t150 - t139 * t300 + t141 * t303;
t135 = -mrSges(2,2) * g(3) - mrSges(2,3) * t280 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t306 - qJ(2) * t145 + t136 * t298 - t137 * t297;
t134 = Ifges(2,6) * qJDD(1) + t306 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t281 + t297 * t136 + t298 * t137 - pkin(1) * (m(3) * t296 + t150) + qJ(2) * t328;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t135 - t301 * t134 - pkin(6) * (t142 * t304 + t143 * t301), t135, t136, t141, -t246 * t335 - t316, t160, -t208 * t263 + t320; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t135 + t304 * t134 + pkin(6) * (-t142 * t301 + t143 * t304), t134, t137, t139, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t270 - Ifges(5,6) * t271 - qJD(3) * t246 - t248 * t334 + t308, t158, -t317; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t315, t315, t318, t348, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t270 - Ifges(5,3) * t271 + qJD(3) * t247 + t248 * t335 - t314, -t309, Ifges(7,5) * t223 + Ifges(7,6) * t262 + Ifges(7,3) * t222 + t264 * t208 - t284 * t210 - t325;];
m_new  = t1;
