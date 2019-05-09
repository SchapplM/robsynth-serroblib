% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP1
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
% Datum: 2019-05-05 17:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:28:28
% EndTime: 2019-05-05 17:28:45
% DurationCPUTime: 8.92s
% Computational Cost: add. (152597->336), mult. (328051->414), div. (0->0), fcn. (213279->10), ass. (0->128)
t342 = -2 * qJD(4);
t304 = sin(qJ(1));
t307 = cos(qJ(1));
t288 = t304 * g(1) - g(2) * t307;
t279 = qJDD(1) * pkin(1) + t288;
t289 = -g(1) * t307 - g(2) * t304;
t309 = qJD(1) ^ 2;
t281 = -pkin(1) * t309 + t289;
t299 = sin(pkin(9));
t301 = cos(pkin(9));
t257 = t299 * t279 + t301 * t281;
t243 = -pkin(2) * t309 + qJDD(1) * pkin(7) + t257;
t297 = -g(3) + qJDD(2);
t303 = sin(qJ(3));
t306 = cos(qJ(3));
t231 = -t303 * t243 + t306 * t297;
t333 = qJD(1) * qJD(3);
t330 = t306 * t333;
t282 = qJDD(1) * t303 + t330;
t203 = (-t282 + t330) * qJ(4) + (t303 * t306 * t309 + qJDD(3)) * pkin(3) + t231;
t232 = t306 * t243 + t303 * t297;
t283 = qJDD(1) * t306 - t303 * t333;
t336 = qJD(1) * t303;
t285 = qJD(3) * pkin(3) - qJ(4) * t336;
t296 = t306 ^ 2;
t204 = -pkin(3) * t296 * t309 + qJ(4) * t283 - qJD(3) * t285 + t232;
t298 = sin(pkin(10));
t300 = cos(pkin(10));
t269 = (t298 * t306 + t300 * t303) * qJD(1);
t192 = t203 * t300 - t298 * t204 + t269 * t342;
t268 = (t298 * t303 - t300 * t306) * qJD(1);
t259 = t282 * t300 + t283 * t298;
t302 = sin(qJ(5));
t305 = cos(qJ(5));
t261 = qJD(3) * t305 - t269 * t302;
t225 = qJD(5) * t261 + qJDD(3) * t302 + t259 * t305;
t262 = qJD(3) * t302 + t269 * t305;
t228 = -mrSges(7,1) * t261 + mrSges(7,2) * t262;
t193 = t298 * t203 + t300 * t204 + t268 * t342;
t246 = pkin(4) * t268 - pkin(8) * t269;
t308 = qJD(3) ^ 2;
t190 = -pkin(4) * t308 + qJDD(3) * pkin(8) - t246 * t268 + t193;
t256 = t301 * t279 - t299 * t281;
t321 = -qJDD(1) * pkin(2) - t256;
t206 = -t283 * pkin(3) + qJDD(4) + t285 * t336 + (-qJ(4) * t296 - pkin(7)) * t309 + t321;
t258 = -t282 * t298 + t283 * t300;
t196 = (qJD(3) * t268 - t259) * pkin(8) + (qJD(3) * t269 - t258) * pkin(4) + t206;
t184 = -t302 * t190 + t305 * t196;
t255 = qJDD(5) - t258;
t267 = qJD(5) + t268;
t180 = -0.2e1 * qJD(6) * t262 + (t261 * t267 - t225) * qJ(6) + (t261 * t262 + t255) * pkin(5) + t184;
t233 = -mrSges(7,2) * t267 + mrSges(7,3) * t261;
t332 = m(7) * t180 + t255 * mrSges(7,1) + t267 * t233;
t177 = -t225 * mrSges(7,3) - t262 * t228 + t332;
t185 = t305 * t190 + t302 * t196;
t210 = Ifges(6,4) * t262 + Ifges(6,2) * t261 + Ifges(6,6) * t267;
t211 = Ifges(7,1) * t262 + Ifges(7,4) * t261 + Ifges(7,5) * t267;
t212 = Ifges(6,1) * t262 + Ifges(6,4) * t261 + Ifges(6,5) * t267;
t224 = -qJD(5) * t262 + qJDD(3) * t305 - t259 * t302;
t235 = pkin(5) * t267 - qJ(6) * t262;
t260 = t261 ^ 2;
t183 = -pkin(5) * t260 + qJ(6) * t224 + 0.2e1 * qJD(6) * t261 - t235 * t267 + t185;
t209 = Ifges(7,4) * t262 + Ifges(7,2) * t261 + Ifges(7,6) * t267;
t319 = -mrSges(7,1) * t180 + mrSges(7,2) * t183 - Ifges(7,5) * t225 - Ifges(7,6) * t224 - Ifges(7,3) * t255 - t262 * t209;
t341 = mrSges(6,1) * t184 - mrSges(6,2) * t185 + Ifges(6,5) * t225 + Ifges(6,6) * t224 + Ifges(6,3) * t255 + pkin(5) * t177 + t262 * t210 - (t212 + t211) * t261 - t319;
t245 = mrSges(5,1) * t268 + mrSges(5,2) * t269;
t264 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t269;
t229 = -mrSges(6,1) * t261 + mrSges(6,2) * t262;
t234 = -mrSges(6,2) * t267 + mrSges(6,3) * t261;
t169 = m(6) * t184 + t255 * mrSges(6,1) + t267 * t234 + (-t228 - t229) * t262 + (-mrSges(6,3) - mrSges(7,3)) * t225 + t332;
t331 = m(7) * t183 + t224 * mrSges(7,3) + t261 * t228;
t236 = mrSges(7,1) * t267 - mrSges(7,3) * t262;
t337 = -mrSges(6,1) * t267 + mrSges(6,3) * t262 - t236;
t339 = -mrSges(6,2) - mrSges(7,2);
t172 = m(6) * t185 + t224 * mrSges(6,3) + t261 * t229 + t339 * t255 + t337 * t267 + t331;
t326 = -t169 * t302 + t305 * t172;
t162 = m(5) * t193 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t258 - qJD(3) * t264 - t245 * t268 + t326;
t263 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t268;
t189 = -qJDD(3) * pkin(4) - pkin(8) * t308 + t269 * t246 - t192;
t187 = -pkin(5) * t224 - qJ(6) * t260 + t235 * t262 + qJDD(6) + t189;
t324 = -m(7) * t187 + t224 * mrSges(7,1) + t261 * t233;
t313 = -m(6) * t189 + t224 * mrSges(6,1) + t339 * t225 + t261 * t234 + t337 * t262 + t324;
t174 = m(5) * t192 + qJDD(3) * mrSges(5,1) - t259 * mrSges(5,3) + qJD(3) * t263 - t269 * t245 + t313;
t154 = t298 * t162 + t300 * t174;
t274 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t303 + Ifges(4,2) * t306) * qJD(1);
t275 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t303 + Ifges(4,4) * t306) * qJD(1);
t207 = Ifges(7,5) * t262 + Ifges(7,6) * t261 + Ifges(7,3) * t267;
t208 = Ifges(6,5) * t262 + Ifges(6,6) * t261 + Ifges(6,3) * t267;
t320 = -mrSges(7,1) * t187 + mrSges(7,3) * t183 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t255 + t267 * t211;
t156 = Ifges(6,4) * t225 + Ifges(6,2) * t224 + Ifges(6,6) * t255 + t267 * t212 - mrSges(6,1) * t189 + mrSges(6,3) * t185 - pkin(5) * (t225 * mrSges(7,2) - t324) + qJ(6) * (-t255 * mrSges(7,2) - t267 * t236 + t331) + (-pkin(5) * t236 - t207 - t208) * t262 + t320;
t318 = mrSges(7,2) * t187 - mrSges(7,3) * t180 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t255 + t261 * t207;
t164 = mrSges(6,2) * t189 - mrSges(6,3) * t184 + Ifges(6,1) * t225 + Ifges(6,4) * t224 + Ifges(6,5) * t255 - qJ(6) * t177 + t261 * t208 + (-t209 - t210) * t267 + t318;
t240 = Ifges(5,4) * t269 - Ifges(5,2) * t268 + Ifges(5,6) * qJD(3);
t241 = Ifges(5,1) * t269 - Ifges(5,4) * t268 + Ifges(5,5) * qJD(3);
t314 = -mrSges(5,1) * t192 + mrSges(5,2) * t193 - Ifges(5,5) * t259 - Ifges(5,6) * t258 - Ifges(5,3) * qJDD(3) - pkin(4) * t313 - pkin(8) * t326 - t305 * t156 - t302 * t164 - t269 * t240 - t268 * t241;
t340 = mrSges(4,1) * t231 - mrSges(4,2) * t232 + Ifges(4,5) * t282 + Ifges(4,6) * t283 + Ifges(4,3) * qJDD(3) + pkin(3) * t154 + (t274 * t303 - t275 * t306) * qJD(1) - t314;
t280 = (-mrSges(4,1) * t306 + mrSges(4,2) * t303) * qJD(1);
t335 = qJD(1) * t306;
t287 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t335;
t152 = m(4) * t231 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t282 + qJD(3) * t287 - t280 * t336 + t154;
t286 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t336;
t327 = t300 * t162 - t174 * t298;
t153 = m(4) * t232 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t283 - qJD(3) * t286 + t280 * t335 + t327;
t328 = -t152 * t303 + t306 * t153;
t145 = m(3) * t257 - mrSges(3,1) * t309 - qJDD(1) * mrSges(3,2) + t328;
t242 = -t309 * pkin(7) + t321;
t166 = t305 * t169 + t302 * t172;
t316 = m(5) * t206 - t258 * mrSges(5,1) + mrSges(5,2) * t259 + t268 * t263 + t264 * t269 + t166;
t312 = -m(4) * t242 + t283 * mrSges(4,1) - mrSges(4,2) * t282 - t286 * t336 + t287 * t335 - t316;
t158 = m(3) * t256 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t309 + t312;
t141 = t299 * t145 + t301 * t158;
t147 = t306 * t152 + t303 * t153;
t329 = t301 * t145 - t158 * t299;
t239 = Ifges(5,5) * t269 - Ifges(5,6) * t268 + Ifges(5,3) * qJD(3);
t142 = mrSges(5,2) * t206 - mrSges(5,3) * t192 + Ifges(5,1) * t259 + Ifges(5,4) * t258 + Ifges(5,5) * qJDD(3) - pkin(8) * t166 - qJD(3) * t240 - t156 * t302 + t164 * t305 - t239 * t268;
t148 = -mrSges(5,1) * t206 + mrSges(5,3) * t193 + Ifges(5,4) * t259 + Ifges(5,2) * t258 + Ifges(5,6) * qJDD(3) - pkin(4) * t166 + qJD(3) * t241 - t269 * t239 - t341;
t273 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t303 + Ifges(4,6) * t306) * qJD(1);
t135 = -mrSges(4,1) * t242 + mrSges(4,3) * t232 + Ifges(4,4) * t282 + Ifges(4,2) * t283 + Ifges(4,6) * qJDD(3) - pkin(3) * t316 + qJ(4) * t327 + qJD(3) * t275 + t298 * t142 + t300 * t148 - t273 * t336;
t137 = mrSges(4,2) * t242 - mrSges(4,3) * t231 + Ifges(4,1) * t282 + Ifges(4,4) * t283 + Ifges(4,5) * qJDD(3) - qJ(4) * t154 - qJD(3) * t274 + t142 * t300 - t148 * t298 + t273 * t335;
t317 = mrSges(3,1) * t256 - mrSges(3,2) * t257 + Ifges(3,3) * qJDD(1) + pkin(2) * t312 + pkin(7) * t328 + t306 * t135 + t303 * t137;
t315 = mrSges(2,1) * t288 - mrSges(2,2) * t289 + Ifges(2,3) * qJDD(1) + pkin(1) * t141 + t317;
t139 = m(2) * t289 - mrSges(2,1) * t309 - qJDD(1) * mrSges(2,2) + t329;
t138 = m(2) * t288 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t309 + t141;
t133 = -mrSges(3,1) * t297 + mrSges(3,3) * t257 + t309 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t147 - t340;
t132 = mrSges(3,2) * t297 - mrSges(3,3) * t256 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t309 - pkin(7) * t147 - t135 * t303 + t137 * t306;
t131 = -mrSges(2,2) * g(3) - mrSges(2,3) * t288 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t309 - qJ(2) * t141 + t132 * t301 - t133 * t299;
t130 = Ifges(2,6) * qJDD(1) + t309 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t289 + t299 * t132 + t301 * t133 - pkin(1) * (m(3) * t297 + t147) + qJ(2) * t329;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t307 * t131 - t304 * t130 - pkin(6) * (t138 * t307 + t139 * t304), t131, t132, t137, t142, t164, -t209 * t267 + t318; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t131 + t307 * t130 + pkin(6) * (-t138 * t304 + t139 * t307), t130, t133, t135, t148, t156, -t262 * t207 + t320; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t315, t315, t317, t340, -t314, t341, -t261 * t211 - t319;];
m_new  = t1;
