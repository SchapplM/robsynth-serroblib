% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:06:16
% EndTime: 2019-05-05 15:06:30
% DurationCPUTime: 6.98s
% Computational Cost: add. (86952->316), mult. (194309->371), div. (0->0), fcn. (131726->8), ass. (0->128)
t298 = qJD(1) ^ 2;
t290 = sin(pkin(9));
t281 = t290 ^ 2;
t291 = cos(pkin(9));
t335 = t291 ^ 2 + t281;
t325 = t335 * mrSges(4,3);
t349 = t298 * t325;
t294 = sin(qJ(1));
t296 = cos(qJ(1));
t268 = g(1) * t294 - t296 * g(2);
t314 = -qJ(2) * t298 + qJDD(2) - t268;
t340 = -pkin(1) - qJ(3);
t348 = -(2 * qJD(1) * qJD(3)) + t340 * qJDD(1) + t314;
t269 = -g(1) * t296 - g(2) * t294;
t347 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t269;
t293 = sin(qJ(4));
t295 = cos(qJ(4));
t318 = t290 * t295 + t291 * t293;
t261 = t318 * qJD(1);
t240 = t290 * g(3) + t348 * t291;
t344 = pkin(3) * t298;
t218 = (-pkin(7) * qJDD(1) - t290 * t344) * t291 + t240;
t241 = -g(3) * t291 + t348 * t290;
t330 = qJDD(1) * t290;
t219 = -pkin(7) * t330 - t281 * t344 + t241;
t189 = t293 * t218 + t295 * t219;
t317 = -t290 * t293 + t291 * t295;
t262 = t317 * qJD(1);
t242 = pkin(4) * t261 - pkin(8) * t262;
t297 = qJD(4) ^ 2;
t184 = -pkin(4) * t297 + qJDD(4) * pkin(8) - t242 * t261 + t189;
t311 = qJDD(3) + t347;
t225 = pkin(3) * t330 + (-t335 * pkin(7) + t340) * t298 + t311;
t332 = qJD(4) * t262;
t243 = -t318 * qJDD(1) - t332;
t333 = qJD(4) * t261;
t244 = t317 * qJDD(1) - t333;
t186 = (-t244 + t333) * pkin(8) + (-t243 + t332) * pkin(4) + t225;
t292 = sin(qJ(5));
t345 = cos(qJ(5));
t180 = -t292 * t184 + t345 * t186;
t181 = t345 * t184 + t292 * t186;
t246 = -t345 * qJD(4) + t262 * t292;
t247 = t292 * qJD(4) + t345 * t262;
t259 = qJD(5) + t261;
t192 = Ifges(7,5) * t247 + Ifges(7,6) * t259 + Ifges(7,3) * t246;
t195 = Ifges(6,4) * t247 - Ifges(6,2) * t246 + Ifges(6,6) * t259;
t197 = Ifges(6,1) * t247 - Ifges(6,4) * t246 + Ifges(6,5) * t259;
t206 = qJD(5) * t247 - t345 * qJDD(4) + t244 * t292;
t207 = -t246 * qJD(5) + t292 * qJDD(4) + t345 * t244;
t213 = mrSges(7,1) * t246 - mrSges(7,3) * t247;
t239 = qJDD(5) - t243;
t212 = pkin(5) * t246 - qJ(6) * t247;
t256 = t259 ^ 2;
t176 = -pkin(5) * t256 + qJ(6) * t239 + 0.2e1 * qJD(6) * t259 - t212 * t246 + t181;
t178 = -t239 * pkin(5) - t256 * qJ(6) + t247 * t212 + qJDD(6) - t180;
t196 = Ifges(7,1) * t247 + Ifges(7,4) * t259 + Ifges(7,5) * t246;
t313 = mrSges(7,1) * t178 - mrSges(7,3) * t176 - Ifges(7,4) * t207 - Ifges(7,2) * t239 - Ifges(7,6) * t206 - t246 * t196;
t220 = -mrSges(7,2) * t246 + mrSges(7,3) * t259;
t322 = -m(7) * t178 + t239 * mrSges(7,1) + t259 * t220;
t223 = -mrSges(7,1) * t259 + mrSges(7,2) * t247;
t327 = m(7) * t176 + t239 * mrSges(7,3) + t259 * t223;
t346 = -(-t195 + t192) * t247 + mrSges(6,1) * t180 - mrSges(6,2) * t181 + Ifges(6,5) * t207 - Ifges(6,6) * t206 + Ifges(6,3) * t239 + pkin(5) * (-mrSges(7,2) * t207 - t213 * t247 + t322) + qJ(6) * (-mrSges(7,2) * t206 - t213 * t246 + t327) + t246 * t197 - t313;
t343 = mrSges(2,1) - mrSges(3,2);
t342 = -mrSges(6,3) - mrSges(7,2);
t341 = -Ifges(2,6) + Ifges(3,5);
t339 = Ifges(4,6) * t290;
t235 = mrSges(5,1) * t261 + mrSges(5,2) * t262;
t252 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t262;
t222 = mrSges(6,1) * t259 - mrSges(6,3) * t247;
t336 = -mrSges(6,1) * t246 - mrSges(6,2) * t247 - t213;
t166 = m(6) * t181 - mrSges(6,2) * t239 + t342 * t206 - t222 * t259 + t336 * t246 + t327;
t221 = -mrSges(6,2) * t259 - mrSges(6,3) * t246;
t168 = m(6) * t180 + mrSges(6,1) * t239 + t342 * t207 + t221 * t259 + t336 * t247 + t322;
t323 = t345 * t166 - t168 * t292;
t157 = m(5) * t189 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t243 - qJD(4) * t252 - t235 * t261 + t323;
t188 = t218 * t295 - t293 * t219;
t251 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t261;
t183 = -qJDD(4) * pkin(4) - pkin(8) * t297 + t262 * t242 - t188;
t179 = -0.2e1 * qJD(6) * t247 + (t246 * t259 - t207) * qJ(6) + (t247 * t259 + t206) * pkin(5) + t183;
t173 = m(7) * t179 + mrSges(7,1) * t206 - t207 * mrSges(7,3) + t220 * t246 - t247 * t223;
t302 = -m(6) * t183 - t206 * mrSges(6,1) - mrSges(6,2) * t207 - t246 * t221 - t222 * t247 - t173;
t163 = m(5) * t188 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t244 + qJD(4) * t251 - t235 * t262 + t302;
t148 = t293 * t157 + t295 * t163;
t160 = t292 * t166 + t345 * t168;
t194 = Ifges(7,4) * t247 + Ifges(7,2) * t259 + Ifges(7,6) * t246;
t338 = -Ifges(6,5) * t247 + Ifges(6,6) * t246 - Ifges(6,3) * t259 - t194;
t334 = t298 * (Ifges(4,5) * t291 - t339);
t329 = qJDD(1) * t291;
t326 = Ifges(3,4) + t339;
t316 = -qJDD(1) * mrSges(4,3) - t298 * (mrSges(4,1) * t290 + mrSges(4,2) * t291);
t145 = m(4) * t240 + t316 * t291 + t148;
t324 = t295 * t157 - t163 * t293;
t146 = m(4) * t241 + t316 * t290 + t324;
t142 = -t145 * t290 + t291 * t146;
t321 = -mrSges(7,1) * t179 + mrSges(7,2) * t176;
t320 = Ifges(4,1) * t291 - Ifges(4,4) * t290;
t319 = Ifges(4,4) * t291 - Ifges(4,2) * t290;
t141 = t145 * t291 + t146 * t290;
t312 = mrSges(7,2) * t178 - mrSges(7,3) * t179 + Ifges(7,1) * t207 + Ifges(7,4) * t239 + Ifges(7,5) * t206 + t259 * t192;
t260 = -qJDD(1) * pkin(1) + t314;
t310 = -m(3) * t260 + t298 * mrSges(3,3) - t141;
t309 = -m(5) * t225 + mrSges(5,1) * t243 - t244 * mrSges(5,2) - t251 * t261 - t262 * t252 - t160;
t153 = -mrSges(6,1) * t183 + mrSges(6,3) * t181 - pkin(5) * t173 + (t196 + t197) * t259 + t338 * t247 + (Ifges(6,6) - Ifges(7,6)) * t239 + (Ifges(6,4) - Ifges(7,5)) * t207 + (-Ifges(6,2) - Ifges(7,3)) * t206 + t321;
t158 = mrSges(6,2) * t183 - mrSges(6,3) * t180 + Ifges(6,1) * t207 - Ifges(6,4) * t206 + Ifges(6,5) * t239 - qJ(6) * t173 - t195 * t259 + t338 * t246 + t312;
t226 = Ifges(5,5) * t262 - Ifges(5,6) * t261 + Ifges(5,3) * qJD(4);
t227 = Ifges(5,4) * t262 - Ifges(5,2) * t261 + Ifges(5,6) * qJD(4);
t138 = mrSges(5,2) * t225 - mrSges(5,3) * t188 + Ifges(5,1) * t244 + Ifges(5,4) * t243 + Ifges(5,5) * qJDD(4) - pkin(8) * t160 - qJD(4) * t227 - t292 * t153 + t345 * t158 - t261 * t226;
t228 = Ifges(5,1) * t262 - Ifges(5,4) * t261 + Ifges(5,5) * qJD(4);
t143 = -mrSges(5,1) * t225 + mrSges(5,3) * t189 + Ifges(5,4) * t244 + Ifges(5,2) * t243 + Ifges(5,6) * qJDD(4) - pkin(4) * t160 + qJD(4) * t228 - t262 * t226 - t346;
t250 = t340 * t298 + t311;
t134 = -mrSges(4,1) * t250 + mrSges(4,3) * t241 + pkin(3) * t309 + pkin(7) * t324 + t319 * qJDD(1) + t293 * t138 + t295 * t143 - t291 * t334;
t136 = mrSges(4,2) * t250 - mrSges(4,3) * t240 - pkin(7) * t148 + t320 * qJDD(1) + t138 * t295 - t143 * t293 - t290 * t334;
t254 = pkin(1) * t298 - t347;
t308 = mrSges(3,2) * t260 - mrSges(3,3) * t254 + Ifges(3,1) * qJDD(1) - qJ(3) * t141 - t134 * t290 + t291 * t136;
t306 = -m(4) * t250 - mrSges(4,1) * t330 - mrSges(4,2) * t329 + t309;
t307 = -mrSges(3,1) * t254 - pkin(2) * (t306 + t349) - qJ(3) * t142 - t134 * t291 - t136 * t290;
t305 = -mrSges(5,1) * t188 + mrSges(5,2) * t189 - Ifges(5,5) * t244 - Ifges(5,6) * t243 - Ifges(5,3) * qJDD(4) - pkin(4) * t302 - pkin(8) * t323 - t345 * t153 - t292 * t158 - t262 * t227 - t261 * t228;
t303 = -m(3) * t254 + t298 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t306;
t304 = -mrSges(2,2) * t269 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t310) + qJ(2) * (t303 - t349) + mrSges(2,1) * t268 + Ifges(2,3) * qJDD(1) + t308;
t301 = -mrSges(4,1) * t240 + mrSges(4,2) * t241 - Ifges(4,5) * t329 - pkin(3) * t148 + t305 + (-t290 * t320 - t291 * t319) * t298;
t300 = -mrSges(3,1) * t260 - pkin(2) * t141 + t301;
t149 = t303 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t325) * t298 + m(2) * t269;
t140 = -m(3) * g(3) + t142;
t137 = m(2) * t268 - mrSges(2,2) * t298 + t343 * qJDD(1) + t310;
t133 = -t300 - qJ(2) * t140 + t341 * t298 + (Ifges(2,5) - t326) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t268;
t132 = mrSges(2,3) * t269 - pkin(1) * t140 + (-Ifges(3,4) + Ifges(2,5)) * t298 - t341 * qJDD(1) + t343 * g(3) + t307;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t296 * t133 - t294 * t132 - pkin(6) * (t137 * t296 + t149 * t294), t133, t308, t136, t138, t158, -t194 * t246 + t312; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t294 * t133 + t296 * t132 + pkin(6) * (-t137 * t294 + t149 * t296), t132, -mrSges(3,3) * g(3) - t298 * Ifges(3,5) + t326 * qJDD(1) + t300, t134, t143, t153, -t247 * t192 - t313; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t304, t304, mrSges(3,2) * g(3) + Ifges(3,4) * t298 + Ifges(3,5) * qJDD(1) - t307, -Ifges(4,6) * t330 - t301, -t305, t346, Ifges(7,5) * t207 + Ifges(7,6) * t239 + Ifges(7,3) * t206 + t194 * t247 - t196 * t259 - t321;];
m_new  = t1;
