% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-05-05 13:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:54:37
% EndTime: 2019-05-05 13:54:56
% DurationCPUTime: 16.36s
% Computational Cost: add. (298827->322), mult. (690958->406), div. (0->0), fcn. (490431->12), ass. (0->141)
t308 = qJD(1) ^ 2;
t297 = sin(pkin(10));
t337 = qJD(1) * t297;
t300 = cos(pkin(10));
t344 = qJD(1) * t300;
t304 = sin(qJ(1));
t306 = cos(qJ(1));
t278 = t304 * g(1) - g(2) * t306;
t275 = qJDD(1) * pkin(1) + t278;
t279 = -g(1) * t306 - g(2) * t304;
t276 = -pkin(1) * t308 + t279;
t298 = sin(pkin(9));
t301 = cos(pkin(9));
t257 = t298 * t275 + t301 * t276;
t244 = -pkin(2) * t308 + qJDD(1) * qJ(3) + t257;
t295 = -g(3) + qJDD(2);
t334 = qJD(1) * qJD(3);
t338 = t300 * t295 - 0.2e1 * t297 * t334;
t341 = pkin(3) * t300;
t223 = (-pkin(7) * qJDD(1) + t308 * t341 - t244) * t297 + t338;
t229 = t297 * t295 + (t244 + 0.2e1 * t334) * t300;
t332 = qJDD(1) * t300;
t291 = t300 ^ 2;
t339 = t291 * t308;
t224 = -pkin(3) * t339 + pkin(7) * t332 + t229;
t303 = sin(qJ(4));
t342 = cos(qJ(4));
t208 = t303 * t223 + t342 * t224;
t331 = t300 * t342;
t266 = -qJD(1) * t331 + t303 * t337;
t319 = t342 * t297 + t300 * t303;
t267 = t319 * qJD(1);
t248 = mrSges(5,1) * t266 + mrSges(5,2) * t267;
t333 = qJDD(1) * t297;
t336 = qJD(4) * t267;
t253 = -qJDD(1) * t331 + t303 * t333 + t336;
t264 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t267;
t247 = pkin(4) * t266 - qJ(5) * t267;
t307 = qJD(4) ^ 2;
t196 = -pkin(4) * t307 + qJDD(4) * qJ(5) - t247 * t266 + t208;
t290 = t297 ^ 2;
t256 = t301 * t275 - t298 * t276;
t322 = qJDD(3) - t256;
t225 = (-pkin(2) - t341) * qJDD(1) + (-qJ(3) + (-t290 - t291) * pkin(7)) * t308 + t322;
t335 = t266 * qJD(4);
t254 = t319 * qJDD(1) - t335;
t200 = (-t254 + t335) * qJ(5) + (t253 + t336) * pkin(4) + t225;
t296 = sin(pkin(11));
t299 = cos(pkin(11));
t262 = qJD(4) * t296 + t267 * t299;
t191 = -0.2e1 * qJD(5) * t262 - t296 * t196 + t299 * t200;
t238 = qJDD(4) * t296 + t254 * t299;
t261 = qJD(4) * t299 - t267 * t296;
t189 = (t261 * t266 - t238) * pkin(8) + (t261 * t262 + t253) * pkin(5) + t191;
t192 = 0.2e1 * qJD(5) * t261 + t299 * t196 + t296 * t200;
t236 = pkin(5) * t266 - pkin(8) * t262;
t237 = qJDD(4) * t299 - t254 * t296;
t260 = t261 ^ 2;
t190 = -pkin(5) * t260 + pkin(8) * t237 - t236 * t266 + t192;
t302 = sin(qJ(6));
t305 = cos(qJ(6));
t187 = t189 * t305 - t190 * t302;
t226 = t261 * t305 - t262 * t302;
t205 = qJD(6) * t226 + t237 * t302 + t238 * t305;
t227 = t261 * t302 + t262 * t305;
t213 = -mrSges(7,1) * t226 + mrSges(7,2) * t227;
t265 = qJD(6) + t266;
t214 = -mrSges(7,2) * t265 + mrSges(7,3) * t226;
t252 = qJDD(6) + t253;
t182 = m(7) * t187 + mrSges(7,1) * t252 - mrSges(7,3) * t205 - t213 * t227 + t214 * t265;
t188 = t189 * t302 + t190 * t305;
t204 = -qJD(6) * t227 + t237 * t305 - t238 * t302;
t215 = mrSges(7,1) * t265 - mrSges(7,3) * t227;
t183 = m(7) * t188 - mrSges(7,2) * t252 + mrSges(7,3) * t204 + t213 * t226 - t215 * t265;
t174 = t305 * t182 + t302 * t183;
t231 = -mrSges(6,1) * t261 + mrSges(6,2) * t262;
t234 = -mrSges(6,2) * t266 + mrSges(6,3) * t261;
t172 = m(6) * t191 + mrSges(6,1) * t253 - mrSges(6,3) * t238 - t231 * t262 + t234 * t266 + t174;
t235 = mrSges(6,1) * t266 - mrSges(6,3) * t262;
t326 = -t182 * t302 + t305 * t183;
t173 = m(6) * t192 - mrSges(6,2) * t253 + mrSges(6,3) * t237 + t231 * t261 - t235 * t266 + t326;
t327 = -t172 * t296 + t299 * t173;
t165 = m(5) * t208 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t253 - qJD(4) * t264 - t248 * t266 + t327;
t207 = t342 * t223 - t303 * t224;
t263 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t266;
t195 = -qJDD(4) * pkin(4) - t307 * qJ(5) + t267 * t247 + qJDD(5) - t207;
t193 = -t237 * pkin(5) - t260 * pkin(8) + t262 * t236 + t195;
t317 = m(7) * t193 - t204 * mrSges(7,1) + mrSges(7,2) * t205 - t226 * t214 + t215 * t227;
t311 = -m(6) * t195 + t237 * mrSges(6,1) - mrSges(6,2) * t238 + t261 * t234 - t235 * t262 - t317;
t178 = m(5) * t207 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t254 + qJD(4) * t263 - t248 * t267 + t311;
t155 = t303 * t165 + t342 * t178;
t228 = -t244 * t297 + t338;
t209 = Ifges(7,5) * t227 + Ifges(7,6) * t226 + Ifges(7,3) * t265;
t211 = Ifges(7,1) * t227 + Ifges(7,4) * t226 + Ifges(7,5) * t265;
t175 = -mrSges(7,1) * t193 + mrSges(7,3) * t188 + Ifges(7,4) * t205 + Ifges(7,2) * t204 + Ifges(7,6) * t252 - t209 * t227 + t211 * t265;
t210 = Ifges(7,4) * t227 + Ifges(7,2) * t226 + Ifges(7,6) * t265;
t176 = mrSges(7,2) * t193 - mrSges(7,3) * t187 + Ifges(7,1) * t205 + Ifges(7,4) * t204 + Ifges(7,5) * t252 + t209 * t226 - t210 * t265;
t217 = Ifges(6,5) * t262 + Ifges(6,6) * t261 + Ifges(6,3) * t266;
t219 = Ifges(6,1) * t262 + Ifges(6,4) * t261 + Ifges(6,5) * t266;
t157 = -mrSges(6,1) * t195 + mrSges(6,3) * t192 + Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t253 - pkin(5) * t317 + pkin(8) * t326 + t305 * t175 + t302 * t176 - t262 * t217 + t266 * t219;
t218 = Ifges(6,4) * t262 + Ifges(6,2) * t261 + Ifges(6,6) * t266;
t159 = mrSges(6,2) * t195 - mrSges(6,3) * t191 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t253 - pkin(8) * t174 - t175 * t302 + t176 * t305 + t217 * t261 - t218 * t266;
t242 = Ifges(5,4) * t267 - Ifges(5,2) * t266 + Ifges(5,6) * qJD(4);
t243 = Ifges(5,1) * t267 - Ifges(5,4) * t266 + Ifges(5,5) * qJD(4);
t313 = -mrSges(5,1) * t207 + mrSges(5,2) * t208 - Ifges(5,5) * t254 + Ifges(5,6) * t253 - Ifges(5,3) * qJDD(4) - pkin(4) * t311 - qJ(5) * t327 - t299 * t157 - t296 * t159 - t267 * t242 - t266 * t243;
t324 = Ifges(4,4) * t297 + Ifges(4,2) * t300;
t325 = Ifges(4,1) * t297 + Ifges(4,4) * t300;
t343 = -mrSges(4,1) * t228 + mrSges(4,2) * t229 - pkin(3) * t155 - (t324 * t337 - t325 * t344) * qJD(1) + t313;
t340 = mrSges(4,2) * t297;
t320 = mrSges(4,3) * qJDD(1) + t308 * (-mrSges(4,1) * t300 + t340);
t153 = m(4) * t228 - t320 * t297 + t155;
t328 = t342 * t165 - t303 * t178;
t154 = m(4) * t229 + t320 * t300 + t328;
t329 = -t153 * t297 + t300 * t154;
t146 = m(3) * t257 - mrSges(3,1) * t308 - qJDD(1) * mrSges(3,2) + t329;
t240 = -qJDD(1) * pkin(2) - t308 * qJ(3) + t322;
t167 = t299 * t172 + t296 * t173;
t315 = m(5) * t225 + t253 * mrSges(5,1) + t254 * mrSges(5,2) + t266 * t263 + t267 * t264 + t167;
t312 = -m(4) * t240 + mrSges(4,1) * t332 - t315 + (t290 * t308 + t339) * mrSges(4,3);
t161 = t312 + (mrSges(3,1) - t340) * qJDD(1) + m(3) * t256 - t308 * mrSges(3,2);
t142 = t298 * t146 + t301 * t161;
t148 = t300 * t153 + t297 * t154;
t330 = t301 * t146 - t161 * t298;
t323 = Ifges(4,5) * t297 + Ifges(4,6) * t300;
t241 = Ifges(5,5) * t267 - Ifges(5,6) * t266 + Ifges(5,3) * qJD(4);
t143 = mrSges(5,2) * t225 - mrSges(5,3) * t207 + Ifges(5,1) * t254 - Ifges(5,4) * t253 + Ifges(5,5) * qJDD(4) - qJ(5) * t167 - qJD(4) * t242 - t157 * t296 + t159 * t299 - t241 * t266;
t316 = -mrSges(7,1) * t187 + mrSges(7,2) * t188 - Ifges(7,5) * t205 - Ifges(7,6) * t204 - Ifges(7,3) * t252 - t227 * t210 + t226 * t211;
t309 = -mrSges(6,1) * t191 + mrSges(6,2) * t192 - Ifges(6,5) * t238 - Ifges(6,6) * t237 - pkin(5) * t174 - t262 * t218 + t261 * t219 + t316;
t149 = mrSges(5,3) * t208 - pkin(4) * t167 + Ifges(5,6) * qJDD(4) - mrSges(5,1) * t225 + t309 + qJD(4) * t243 + Ifges(5,4) * t254 - t267 * t241 + (-Ifges(5,2) - Ifges(6,3)) * t253;
t272 = t323 * qJD(1);
t135 = -mrSges(4,1) * t240 + mrSges(4,3) * t229 - pkin(3) * t315 + pkin(7) * t328 + t324 * qJDD(1) + t303 * t143 + t342 * t149 - t272 * t337;
t138 = mrSges(4,2) * t240 - mrSges(4,3) * t228 - pkin(7) * t155 + t325 * qJDD(1) + t342 * t143 - t303 * t149 + t272 * t344;
t318 = -mrSges(3,2) * t257 + qJ(3) * t329 + t300 * t135 + t297 * t138 + pkin(2) * (-mrSges(4,2) * t333 + t312) + mrSges(3,1) * t256 + Ifges(3,3) * qJDD(1);
t314 = mrSges(2,1) * t278 - mrSges(2,2) * t279 + Ifges(2,3) * qJDD(1) + pkin(1) * t142 + t318;
t140 = m(2) * t279 - mrSges(2,1) * t308 - qJDD(1) * mrSges(2,2) + t330;
t139 = m(2) * t278 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t308 + t142;
t136 = -pkin(2) * t148 + (Ifges(3,6) - t323) * qJDD(1) + mrSges(3,3) * t257 - mrSges(3,1) * t295 + t308 * Ifges(3,5) + t343;
t133 = mrSges(3,2) * t295 - mrSges(3,3) * t256 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t308 - qJ(3) * t148 - t135 * t297 + t138 * t300;
t132 = -mrSges(2,2) * g(3) - mrSges(2,3) * t278 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t308 - qJ(2) * t142 + t133 * t301 - t136 * t298;
t131 = Ifges(2,6) * qJDD(1) + t308 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t279 + t298 * t133 + t301 * t136 - pkin(1) * (m(3) * t295 + t148) + qJ(2) * t330;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t306 * t132 - t304 * t131 - pkin(6) * (t139 * t306 + t140 * t304), t132, t133, t138, t143, t159, t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t304 * t132 + t306 * t131 + pkin(6) * (-t139 * t304 + t140 * t306), t131, t136, t135, t149, t157, t175; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t314, t314, t318, t323 * qJDD(1) - t343, -t313, Ifges(6,3) * t253 - t309, -t316;];
m_new  = t1;
