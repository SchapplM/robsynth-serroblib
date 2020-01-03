% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:07
% EndTime: 2019-12-31 19:54:19
% DurationCPUTime: 6.68s
% Computational Cost: add. (78343->308), mult. (182676->380), div. (0->0), fcn. (123708->8), ass. (0->113)
t265 = sin(qJ(2));
t267 = cos(qJ(2));
t287 = qJD(1) * qJD(2);
t242 = t265 * qJDD(1) + t267 * t287;
t266 = sin(qJ(1));
t268 = cos(qJ(1));
t249 = -t268 * g(1) - t266 * g(2);
t269 = qJD(1) ^ 2;
t237 = -t269 * pkin(1) + qJDD(1) * pkin(6) + t249;
t292 = t265 * t237;
t294 = pkin(2) * t269;
t198 = qJDD(2) * pkin(2) - t242 * qJ(3) - t292 + (qJ(3) * t287 + t265 * t294 - g(3)) * t267;
t222 = -t265 * g(3) + t267 * t237;
t243 = t267 * qJDD(1) - t265 * t287;
t289 = qJD(1) * t265;
t245 = qJD(2) * pkin(2) - qJ(3) * t289;
t261 = t267 ^ 2;
t199 = t243 * qJ(3) - qJD(2) * t245 - t261 * t294 + t222;
t262 = sin(pkin(8));
t263 = cos(pkin(8));
t232 = (t262 * t267 + t263 * t265) * qJD(1);
t164 = -0.2e1 * qJD(3) * t232 + t263 * t198 - t262 * t199;
t220 = t263 * t242 + t262 * t243;
t231 = (-t262 * t265 + t263 * t267) * qJD(1);
t158 = (qJD(2) * t231 - t220) * pkin(7) + (t231 * t232 + qJDD(2)) * pkin(3) + t164;
t165 = 0.2e1 * qJD(3) * t231 + t262 * t198 + t263 * t199;
t219 = -t262 * t242 + t263 * t243;
t225 = qJD(2) * pkin(3) - t232 * pkin(7);
t230 = t231 ^ 2;
t160 = -t230 * pkin(3) + t219 * pkin(7) - qJD(2) * t225 + t165;
t264 = sin(qJ(4));
t295 = cos(qJ(4));
t156 = t264 * t158 + t295 * t160;
t211 = t264 * t231 + t295 * t232;
t178 = t211 * qJD(4) - t295 * t219 + t264 * t220;
t258 = qJD(2) + qJD(4);
t204 = t258 * mrSges(5,1) - t211 * mrSges(5,3);
t210 = -t295 * t231 + t264 * t232;
t257 = qJDD(2) + qJDD(4);
t191 = t210 * pkin(4) - t211 * qJ(5);
t256 = t258 ^ 2;
t149 = -t256 * pkin(4) + t257 * qJ(5) + 0.2e1 * qJD(5) * t258 - t210 * t191 + t156;
t205 = -t258 * mrSges(6,1) + t211 * mrSges(6,2);
t286 = m(6) * t149 + t257 * mrSges(6,3) + t258 * t205;
t192 = t210 * mrSges(6,1) - t211 * mrSges(6,3);
t290 = -t210 * mrSges(5,1) - t211 * mrSges(5,2) - t192;
t293 = -mrSges(5,3) - mrSges(6,2);
t139 = m(5) * t156 - t257 * mrSges(5,2) + t293 * t178 - t258 * t204 + t290 * t210 + t286;
t155 = t295 * t158 - t264 * t160;
t179 = -t210 * qJD(4) + t264 * t219 + t295 * t220;
t203 = -t258 * mrSges(5,2) - t210 * mrSges(5,3);
t151 = -t257 * pkin(4) - t256 * qJ(5) + t211 * t191 + qJDD(5) - t155;
t206 = -t210 * mrSges(6,2) + t258 * mrSges(6,3);
t282 = -m(6) * t151 + t257 * mrSges(6,1) + t258 * t206;
t141 = m(5) * t155 + t257 * mrSges(5,1) + t293 * t179 + t258 * t203 + t290 * t211 + t282;
t134 = t264 * t139 + t295 * t141;
t214 = -t231 * mrSges(4,1) + t232 * mrSges(4,2);
t223 = -qJD(2) * mrSges(4,2) + t231 * mrSges(4,3);
t129 = m(4) * t164 + qJDD(2) * mrSges(4,1) - t220 * mrSges(4,3) + qJD(2) * t223 - t232 * t214 + t134;
t224 = qJD(2) * mrSges(4,1) - t232 * mrSges(4,3);
t283 = t295 * t139 - t264 * t141;
t130 = m(4) * t165 - qJDD(2) * mrSges(4,2) + t219 * mrSges(4,3) - qJD(2) * t224 + t231 * t214 + t283;
t125 = t263 * t129 + t262 * t130;
t221 = -t267 * g(3) - t292;
t234 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t265 + Ifges(3,2) * t267) * qJD(1);
t235 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t265 + Ifges(3,4) * t267) * qJD(1);
t208 = Ifges(4,4) * t232 + Ifges(4,2) * t231 + Ifges(4,6) * qJD(2);
t209 = Ifges(4,1) * t232 + Ifges(4,4) * t231 + Ifges(4,5) * qJD(2);
t184 = Ifges(5,4) * t211 - Ifges(5,2) * t210 + Ifges(5,6) * t258;
t186 = Ifges(5,1) * t211 - Ifges(5,4) * t210 + Ifges(5,5) * t258;
t181 = Ifges(6,5) * t211 + Ifges(6,6) * t258 + Ifges(6,3) * t210;
t185 = Ifges(6,1) * t211 + Ifges(6,4) * t258 + Ifges(6,5) * t210;
t275 = mrSges(6,1) * t151 - mrSges(6,3) * t149 - Ifges(6,4) * t179 - Ifges(6,2) * t257 - Ifges(6,6) * t178 + t211 * t181 - t210 * t185;
t273 = mrSges(5,2) * t156 - t210 * t186 - qJ(5) * (-t178 * mrSges(6,2) - t210 * t192 + t286) - pkin(4) * (-t179 * mrSges(6,2) - t211 * t192 + t282) - mrSges(5,1) * t155 - t211 * t184 + Ifges(5,6) * t178 - Ifges(5,5) * t179 - Ifges(5,3) * t257 + t275;
t271 = -mrSges(4,1) * t164 + mrSges(4,2) * t165 - Ifges(4,5) * t220 - Ifges(4,6) * t219 - Ifges(4,3) * qJDD(2) - pkin(3) * t134 - t232 * t208 + t231 * t209 + t273;
t296 = mrSges(3,1) * t221 - mrSges(3,2) * t222 + Ifges(3,5) * t242 + Ifges(3,6) * t243 + Ifges(3,3) * qJDD(2) + pkin(2) * t125 + (t265 * t234 - t267 * t235) * qJD(1) - t271;
t183 = Ifges(6,4) * t211 + Ifges(6,2) * t258 + Ifges(6,6) * t210;
t291 = -Ifges(5,5) * t211 + Ifges(5,6) * t210 - Ifges(5,3) * t258 - t183;
t288 = qJD(1) * t267;
t248 = t266 * g(1) - t268 * g(2);
t241 = (-mrSges(3,1) * t267 + mrSges(3,2) * t265) * qJD(1);
t247 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t288;
t123 = m(3) * t221 + qJDD(2) * mrSges(3,1) - t242 * mrSges(3,3) + qJD(2) * t247 - t241 * t289 + t125;
t246 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t289;
t284 = -t262 * t129 + t263 * t130;
t124 = m(3) * t222 - qJDD(2) * mrSges(3,2) + t243 * mrSges(3,3) - qJD(2) * t246 + t241 * t288 + t284;
t285 = -t265 * t123 + t267 * t124;
t279 = -qJDD(1) * pkin(1) - t248;
t202 = -t243 * pkin(2) + qJDD(3) + t245 * t289 + (-qJ(3) * t261 - pkin(6)) * t269 + t279;
t162 = -t219 * pkin(3) - t230 * pkin(7) + t232 * t225 + t202;
t153 = -0.2e1 * qJD(5) * t211 + (t210 * t258 - t179) * qJ(5) + (t211 * t258 + t178) * pkin(4) + t162;
t281 = -mrSges(6,1) * t153 + mrSges(6,2) * t149;
t142 = m(6) * t153 + t178 * mrSges(6,1) - t179 * mrSges(6,3) - t211 * t205 + t210 * t206;
t278 = mrSges(6,2) * t151 - mrSges(6,3) * t153 + Ifges(6,1) * t179 + Ifges(6,4) * t257 + Ifges(6,5) * t178 + t258 * t181;
t131 = -mrSges(5,1) * t162 + mrSges(5,3) * t156 - pkin(4) * t142 + (t185 + t186) * t258 + (Ifges(5,6) - Ifges(6,6)) * t257 + t291 * t211 + (Ifges(5,4) - Ifges(6,5)) * t179 + (-Ifges(5,2) - Ifges(6,3)) * t178 + t281;
t132 = mrSges(5,2) * t162 - mrSges(5,3) * t155 + Ifges(5,1) * t179 - Ifges(5,4) * t178 + Ifges(5,5) * t257 - qJ(5) * t142 - t258 * t184 + t291 * t210 + t278;
t207 = Ifges(4,5) * t232 + Ifges(4,6) * t231 + Ifges(4,3) * qJD(2);
t276 = m(5) * t162 + t178 * mrSges(5,1) + t179 * mrSges(5,2) + t210 * t203 + t211 * t204 + t142;
t120 = -mrSges(4,1) * t202 + mrSges(4,3) * t165 + Ifges(4,4) * t220 + Ifges(4,2) * t219 + Ifges(4,6) * qJDD(2) - pkin(3) * t276 + pkin(7) * t283 + qJD(2) * t209 + t295 * t131 + t264 * t132 - t232 * t207;
t121 = mrSges(4,2) * t202 - mrSges(4,3) * t164 + Ifges(4,1) * t220 + Ifges(4,4) * t219 + Ifges(4,5) * qJDD(2) - pkin(7) * t134 - qJD(2) * t208 - t264 * t131 + t295 * t132 + t231 * t207;
t233 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t265 + Ifges(3,6) * t267) * qJD(1);
t236 = -t269 * pkin(6) + t279;
t274 = m(4) * t202 - t219 * mrSges(4,1) + t220 * mrSges(4,2) - t231 * t223 + t232 * t224 + t276;
t113 = -mrSges(3,1) * t236 + mrSges(3,3) * t222 + Ifges(3,4) * t242 + Ifges(3,2) * t243 + Ifges(3,6) * qJDD(2) - pkin(2) * t274 + qJ(3) * t284 + qJD(2) * t235 + t263 * t120 + t262 * t121 - t233 * t289;
t116 = mrSges(3,2) * t236 - mrSges(3,3) * t221 + Ifges(3,1) * t242 + Ifges(3,4) * t243 + Ifges(3,5) * qJDD(2) - qJ(3) * t125 - qJD(2) * t234 - t262 * t120 + t263 * t121 + t233 * t288;
t272 = -m(3) * t236 + t243 * mrSges(3,1) - t242 * mrSges(3,2) - t246 * t289 + t247 * t288 - t274;
t277 = mrSges(2,1) * t248 - mrSges(2,2) * t249 + Ifges(2,3) * qJDD(1) + pkin(1) * t272 + pkin(6) * t285 + t267 * t113 + t265 * t116;
t135 = m(2) * t248 + qJDD(1) * mrSges(2,1) - t269 * mrSges(2,2) + t272;
t119 = t267 * t123 + t265 * t124;
t117 = m(2) * t249 - t269 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t285;
t114 = mrSges(2,1) * g(3) + mrSges(2,3) * t249 + t269 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t119 - t296;
t111 = -mrSges(2,2) * g(3) - mrSges(2,3) * t248 + Ifges(2,5) * qJDD(1) - t269 * Ifges(2,6) - pkin(6) * t119 - t265 * t113 + t267 * t116;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t268 * t111 - t266 * t114 - pkin(5) * (t266 * t117 + t268 * t135), t111, t116, t121, t132, -t210 * t183 + t278; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t266 * t111 + t268 * t114 + pkin(5) * (t268 * t117 - t266 * t135), t114, t113, t120, t131, -t275; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t277, t277, t296, -t271, -t273, Ifges(6,5) * t179 + Ifges(6,6) * t257 + Ifges(6,3) * t178 + t211 * t183 - t258 * t185 - t281;];
m_new = t1;
