% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-05-05 21:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:25:27
% EndTime: 2019-05-05 21:25:37
% DurationCPUTime: 4.97s
% Computational Cost: add. (65675->344), mult. (123586->388), div. (0->0), fcn. (72171->8), ass. (0->127)
t292 = sin(qJ(1));
t294 = cos(qJ(1));
t276 = t292 * g(1) - g(2) * t294;
t267 = qJDD(1) * pkin(1) + t276;
t277 = -g(1) * t294 - g(2) * t292;
t296 = qJD(1) ^ 2;
t269 = -pkin(1) * t296 + t277;
t288 = sin(pkin(9));
t289 = cos(pkin(9));
t229 = t267 * t289 - t288 * t269;
t198 = -qJDD(1) * pkin(2) - pkin(7) * t296 - t229;
t291 = sin(qJ(3));
t293 = cos(qJ(3));
t319 = qJD(1) * qJD(3);
t316 = t293 * t319;
t271 = qJDD(1) * t291 + t316;
t317 = t291 * t319;
t272 = qJDD(1) * t293 - t317;
t183 = (-t271 - t316) * pkin(8) + (-t272 + t317) * pkin(3) + t198;
t230 = t288 * t267 + t289 * t269;
t199 = -pkin(2) * t296 + qJDD(1) * pkin(7) + t230;
t287 = -g(3) + qJDD(2);
t192 = t293 * t199 + t291 * t287;
t270 = (-pkin(3) * t293 - pkin(8) * t291) * qJD(1);
t295 = qJD(3) ^ 2;
t320 = qJD(1) * t293;
t187 = -pkin(3) * t295 + qJDD(3) * pkin(8) + t270 * t320 + t192;
t290 = sin(qJ(4));
t331 = cos(qJ(4));
t180 = t331 * t183 - t290 * t187;
t321 = qJD(1) * t291;
t265 = -t331 * qJD(3) + t290 * t321;
t266 = t290 * qJD(3) + t331 * t321;
t235 = pkin(4) * t265 - qJ(5) * t266;
t264 = qJDD(4) - t272;
t279 = -qJD(4) + t320;
t278 = t279 ^ 2;
t178 = -t264 * pkin(4) - t278 * qJ(5) + t266 * t235 + qJDD(5) - t180;
t227 = -t265 * qJD(4) + t290 * qJDD(3) + t331 * t271;
t237 = -mrSges(6,2) * t265 - mrSges(6,3) * t266;
t336 = -m(6) * t178 - t227 * mrSges(6,1) - t266 * t237;
t244 = mrSges(6,1) * t266 - mrSges(6,2) * t279;
t226 = qJD(4) * t266 - t331 * qJDD(3) + t271 * t290;
t240 = pkin(5) * t266 + qJ(6) * t279;
t263 = t265 ^ 2;
t191 = -t291 * t199 + t287 * t293;
t186 = -qJDD(3) * pkin(3) - pkin(8) * t295 + t270 * t321 - t191;
t327 = t265 * t279;
t333 = -2 * qJD(5);
t301 = (-t227 - t327) * qJ(5) + t186 + (-t279 * pkin(4) + t333) * t266;
t332 = 2 * qJD(6);
t175 = -pkin(5) * t263 + t265 * t332 - t240 * t266 + (pkin(4) + qJ(6)) * t226 + t301;
t241 = mrSges(7,1) * t266 + mrSges(7,3) * t279;
t243 = -mrSges(7,1) * t265 - mrSges(7,2) * t279;
t166 = m(7) * t175 - t227 * mrSges(7,2) + t226 * mrSges(7,3) - t266 * t241 + t265 * t243;
t179 = pkin(4) * t226 + t301;
t242 = mrSges(6,1) * t265 + mrSges(6,3) * t279;
t306 = -m(6) * t179 + t226 * mrSges(6,2) + t265 * t242 - t166;
t163 = -mrSges(6,3) * t227 - t244 * t266 - t306;
t181 = t290 * t183 + t331 * t187;
t200 = -Ifges(7,5) * t279 + Ifges(7,6) * t265 + Ifges(7,3) * t266;
t202 = Ifges(5,5) * t266 - Ifges(5,6) * t265 - Ifges(5,3) * t279;
t207 = -Ifges(6,1) * t279 - Ifges(6,4) * t266 + Ifges(6,5) * t265;
t305 = -pkin(4) * t278 + qJ(5) * t264 - t235 * t265 + t181;
t176 = 0.2e1 * qJD(5) * t279 - t305;
t234 = -mrSges(7,2) * t266 + mrSges(7,3) * t265;
t172 = -pkin(5) * t226 - qJ(6) * t263 + qJDD(6) + (t333 - t240) * t279 + t305;
t206 = -Ifges(7,1) * t279 + Ifges(7,4) * t265 + Ifges(7,5) * t266;
t310 = mrSges(7,1) * t172 - mrSges(7,3) * t175 - Ifges(7,4) * t264 - Ifges(7,2) * t226 - Ifges(7,6) * t227 - t266 * t206;
t318 = m(7) * t172 + t264 * mrSges(7,2) - t279 * t241;
t300 = mrSges(6,1) * t176 - mrSges(6,2) * t179 + pkin(5) * (mrSges(7,1) * t226 + t234 * t265 - t318) + qJ(6) * t166 - t310;
t204 = -Ifges(6,4) * t279 - Ifges(6,2) * t266 + Ifges(6,6) * t265;
t323 = Ifges(5,1) * t266 - Ifges(5,4) * t265 - Ifges(5,5) * t279 - t204;
t328 = Ifges(5,4) + Ifges(6,6);
t141 = -t300 + mrSges(5,3) * t181 - mrSges(5,1) * t186 - pkin(4) * t163 + (-t200 - t323) * t279 + (-t202 - t207) * t266 + (Ifges(5,6) - Ifges(6,5)) * t264 + t328 * t227 + (-Ifges(5,2) - Ifges(6,3)) * t226;
t201 = -Ifges(6,5) * t279 - Ifges(6,6) * t266 + Ifges(6,3) * t265;
t205 = Ifges(5,4) * t266 - Ifges(5,2) * t265 - Ifges(5,6) * t279;
t169 = t279 * t332 + (t265 * t266 - t264) * qJ(6) + (t227 - t327) * pkin(5) + t178;
t313 = -m(7) * t169 + t264 * mrSges(7,3) - t279 * t243;
t165 = mrSges(7,1) * t227 + t234 * t266 - t313;
t203 = -Ifges(7,4) * t279 + Ifges(7,2) * t265 + Ifges(7,6) * t266;
t309 = -mrSges(7,1) * t169 + mrSges(7,2) * t175 - Ifges(7,5) * t264 - Ifges(7,6) * t226 - Ifges(7,3) * t227 + t279 * t203;
t303 = -mrSges(6,1) * t178 + mrSges(6,3) * t179 - pkin(5) * t165 + t309;
t324 = t206 + t207;
t147 = -t303 + mrSges(5,2) * t186 - mrSges(5,3) * t180 - qJ(5) * t163 - t328 * t226 + (Ifges(5,1) + Ifges(6,2)) * t227 + (Ifges(5,5) - Ifges(6,4)) * t264 + (-t202 - t324) * t265 + (t205 - t201) * t279;
t236 = mrSges(5,1) * t265 + mrSges(5,2) * t266;
t238 = mrSges(5,2) * t279 - mrSges(5,3) * t265;
t329 = -mrSges(7,1) - mrSges(5,3);
t158 = m(5) * t180 + (-t238 + t242) * t279 + (-t234 - t236) * t266 + (mrSges(5,1) - mrSges(6,2)) * t264 + t329 * t227 + t313 + t336;
t239 = -mrSges(5,1) * t279 - mrSges(5,3) * t266;
t311 = -m(6) * t176 + t264 * mrSges(6,3) - t279 * t244 + t318;
t322 = -t234 - t237;
t160 = m(5) * t181 - mrSges(5,2) * t264 + t239 * t279 + (-t236 + t322) * t265 + (-mrSges(6,1) + t329) * t226 + t311;
t155 = -t158 * t290 + t331 * t160;
t161 = -m(5) * t186 - t226 * mrSges(5,1) - t265 * t238 + (-t239 + t244) * t266 + (-mrSges(5,2) + mrSges(6,3)) * t227 + t306;
t252 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t291 + Ifges(4,2) * t293) * qJD(1);
t253 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t291 + Ifges(4,4) * t293) * qJD(1);
t335 = mrSges(4,1) * t191 - mrSges(4,2) * t192 + Ifges(4,5) * t271 + Ifges(4,6) * t272 + Ifges(4,3) * qJDD(3) + pkin(3) * t161 + pkin(8) * t155 + (t252 * t291 - t253 * t293) * qJD(1) + t331 * t141 + t290 * t147;
t308 = -mrSges(7,2) * t172 + mrSges(7,3) * t169 - Ifges(7,1) * t264 - Ifges(7,4) * t226 - Ifges(7,5) * t227 - t265 * t200;
t299 = -mrSges(6,2) * t178 + mrSges(6,3) * t176 - Ifges(6,1) * t264 + Ifges(6,4) * t227 - Ifges(6,5) * t226 + qJ(6) * t165 + t266 * t201 + t308;
t334 = t323 * t265 + (t205 - t203) * t266 + mrSges(5,1) * t180 - mrSges(5,2) * t181 + Ifges(5,5) * t227 - Ifges(5,6) * t226 + Ifges(5,3) * t264 + pkin(4) * (-mrSges(6,2) * t264 + t242 * t279 - t165 + t336) + qJ(5) * (t322 * t265 + (-mrSges(6,1) - mrSges(7,1)) * t226 + t311) - t299;
t326 = t266 * t203;
t268 = (-mrSges(4,1) * t293 + mrSges(4,2) * t291) * qJD(1);
t274 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t321;
t153 = m(4) * t192 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t272 - qJD(3) * t274 + t268 * t320 + t155;
t275 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t320;
t157 = m(4) * t191 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t271 + qJD(3) * t275 - t268 * t321 + t161;
t314 = t293 * t153 - t157 * t291;
t144 = m(3) * t230 - mrSges(3,1) * t296 - qJDD(1) * mrSges(3,2) + t314;
t154 = t331 * t158 + t290 * t160;
t302 = -m(4) * t198 + t272 * mrSges(4,1) - t271 * mrSges(4,2) - t274 * t321 + t275 * t320 - t154;
t149 = m(3) * t229 + qJDD(1) * mrSges(3,1) - t296 * mrSges(3,2) + t302;
t140 = t288 * t144 + t289 * t149;
t146 = t291 * t153 + t293 * t157;
t315 = t289 * t144 - t149 * t288;
t251 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t291 + Ifges(4,6) * t293) * qJD(1);
t134 = mrSges(4,2) * t198 - mrSges(4,3) * t191 + Ifges(4,1) * t271 + Ifges(4,4) * t272 + Ifges(4,5) * qJDD(3) - pkin(8) * t154 - qJD(3) * t252 - t290 * t141 + t331 * t147 + t251 * t320;
t136 = -mrSges(4,1) * t198 + mrSges(4,3) * t192 + Ifges(4,4) * t271 + Ifges(4,2) * t272 + Ifges(4,6) * qJDD(3) - pkin(3) * t154 + qJD(3) * t253 - t251 * t321 - t334;
t307 = mrSges(3,1) * t229 - mrSges(3,2) * t230 + Ifges(3,3) * qJDD(1) + pkin(2) * t302 + pkin(7) * t314 + t291 * t134 + t293 * t136;
t304 = mrSges(2,1) * t276 - mrSges(2,2) * t277 + Ifges(2,3) * qJDD(1) + pkin(1) * t140 + t307;
t138 = m(2) * t277 - mrSges(2,1) * t296 - qJDD(1) * mrSges(2,2) + t315;
t137 = m(2) * t276 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t296 + t140;
t132 = -mrSges(3,1) * t287 + mrSges(3,3) * t230 + t296 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t146 - t335;
t131 = mrSges(3,2) * t287 - mrSges(3,3) * t229 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t296 - pkin(7) * t146 + t134 * t293 - t136 * t291;
t130 = -mrSges(2,2) * g(3) - mrSges(2,3) * t276 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t296 - qJ(2) * t140 + t131 * t289 - t132 * t288;
t129 = Ifges(2,6) * qJDD(1) + t296 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t277 + t288 * t131 + t289 * t132 - pkin(1) * (m(3) * t287 + t146) + qJ(2) * t315;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t294 * t130 - t292 * t129 - pkin(6) * (t137 * t294 + t138 * t292), t130, t131, t134, t147, -t265 * t204 - t299 - t326, -t308 - t326; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t292 * t130 + t294 * t129 + pkin(6) * (-t137 * t292 + t138 * t294), t129, t132, t136, t141, Ifges(6,4) * t264 - Ifges(6,2) * t227 + Ifges(6,6) * t226 + t279 * t201 + t324 * t265 + t303, t279 * t200 - t310; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t304, t304, t307, t335, t334, t300 + t266 * t207 + Ifges(6,5) * t264 + Ifges(6,3) * t226 - Ifges(6,6) * t227 + (t200 - t204) * t279, -t265 * t206 - t309;];
m_new  = t1;
