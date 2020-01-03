% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP9
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:47
% EndTime: 2019-12-31 20:05:59
% DurationCPUTime: 5.60s
% Computational Cost: add. (67203->306), mult. (144142->374), div. (0->0), fcn. (94802->8), ass. (0->112)
t257 = sin(qJ(1));
t259 = cos(qJ(1));
t246 = -t259 * g(1) - t257 * g(2);
t261 = qJD(1) ^ 2;
t226 = -t261 * pkin(1) + qJDD(1) * pkin(6) + t246;
t256 = sin(qJ(2));
t258 = cos(qJ(2));
t208 = -t258 * g(3) - t256 * t226;
t238 = (-pkin(2) * t258 - qJ(3) * t256) * qJD(1);
t260 = qJD(2) ^ 2;
t279 = qJD(1) * t256;
t193 = -qJDD(2) * pkin(2) - t260 * qJ(3) + t238 * t279 + qJDD(3) - t208;
t277 = qJD(1) * qJD(2);
t275 = t258 * t277;
t240 = t256 * qJDD(1) + t275;
t253 = sin(pkin(8));
t254 = cos(pkin(8));
t213 = t254 * qJDD(2) - t253 * t240;
t235 = t253 * qJD(2) + t254 * t279;
t278 = t258 * qJD(1);
t215 = -pkin(3) * t278 - t235 * pkin(7);
t234 = t254 * qJD(2) - t253 * t279;
t233 = t234 ^ 2;
t160 = -t213 * pkin(3) - t233 * pkin(7) + t235 * t215 + t193;
t255 = sin(qJ(4));
t283 = cos(qJ(4));
t205 = t255 * t234 + t283 * t235;
t214 = t253 * qJDD(2) + t254 * t240;
t173 = t205 * qJD(4) - t283 * t213 + t255 * t214;
t204 = -t283 * t234 + t255 * t235;
t174 = -t204 * qJD(4) + t255 * t213 + t283 * t214;
t248 = qJD(4) - t278;
t150 = -0.2e1 * qJD(5) * t205 + (t204 * t248 - t174) * qJ(5) + (t205 * t248 + t173) * pkin(4) + t160;
t197 = -t248 * mrSges(6,1) + t205 * mrSges(6,2);
t198 = -t204 * mrSges(6,2) + t248 * mrSges(6,3);
t143 = m(6) * t150 + t173 * mrSges(6,1) - t174 * mrSges(6,3) - t205 * t197 + t204 * t198;
t245 = t257 * g(1) - t259 * g(2);
t225 = -qJDD(1) * pkin(1) - t261 * pkin(6) - t245;
t249 = t256 * t277;
t241 = t258 * qJDD(1) - t249;
t189 = (-t240 - t275) * qJ(3) + (-t241 + t249) * pkin(2) + t225;
t209 = -t256 * g(3) + t258 * t226;
t194 = -t260 * pkin(2) + qJDD(2) * qJ(3) + t238 * t278 + t209;
t158 = -0.2e1 * qJD(3) * t235 + t254 * t189 - t253 * t194;
t155 = (-t234 * t278 - t214) * pkin(7) + (t234 * t235 - t241) * pkin(3) + t158;
t159 = 0.2e1 * qJD(3) * t234 + t253 * t189 + t254 * t194;
t157 = -t233 * pkin(3) + t213 * pkin(7) + t215 * t278 + t159;
t153 = t255 * t155 + t283 * t157;
t179 = Ifges(6,1) * t205 + Ifges(6,4) * t248 + Ifges(6,5) * t204;
t180 = Ifges(5,1) * t205 - Ifges(5,4) * t204 + Ifges(5,5) * t248;
t237 = qJDD(4) - t241;
t184 = t204 * pkin(4) - t205 * qJ(5);
t247 = t248 ^ 2;
t146 = -t247 * pkin(4) + t237 * qJ(5) + 0.2e1 * qJD(5) * t248 - t204 * t184 + t153;
t271 = -mrSges(6,1) * t150 + mrSges(6,2) * t146;
t177 = Ifges(6,4) * t205 + Ifges(6,2) * t248 + Ifges(6,6) * t204;
t281 = -Ifges(5,5) * t205 + Ifges(5,6) * t204 - Ifges(5,3) * t248 - t177;
t129 = -mrSges(5,1) * t160 + mrSges(5,3) * t153 - pkin(4) * t143 + (t179 + t180) * t248 + (Ifges(5,6) - Ifges(6,6)) * t237 + t281 * t205 + (Ifges(5,4) - Ifges(6,5)) * t174 + (-Ifges(5,2) - Ifges(6,3)) * t173 + t271;
t152 = t283 * t155 - t255 * t157;
t178 = Ifges(5,4) * t205 - Ifges(5,2) * t204 + Ifges(5,6) * t248;
t148 = -t237 * pkin(4) - t247 * qJ(5) + t205 * t184 + qJDD(5) - t152;
t175 = Ifges(6,5) * t205 + Ifges(6,6) * t248 + Ifges(6,3) * t204;
t269 = mrSges(6,2) * t148 - mrSges(6,3) * t150 + Ifges(6,1) * t174 + Ifges(6,4) * t237 + Ifges(6,5) * t173 + t248 * t175;
t130 = mrSges(5,2) * t160 - mrSges(5,3) * t152 + Ifges(5,1) * t174 - Ifges(5,4) * t173 + Ifges(5,5) * t237 - qJ(5) * t143 - t248 * t178 + t281 * t204 + t269;
t199 = Ifges(4,5) * t235 + Ifges(4,6) * t234 - Ifges(4,3) * t278;
t201 = Ifges(4,1) * t235 + Ifges(4,4) * t234 - Ifges(4,5) * t278;
t195 = -t248 * mrSges(5,2) - t204 * mrSges(5,3);
t196 = t248 * mrSges(5,1) - t205 * mrSges(5,3);
t266 = m(5) * t160 + t173 * mrSges(5,1) + t174 * mrSges(5,2) + t204 * t195 + t205 * t196 + t143;
t276 = m(6) * t146 + t237 * mrSges(6,3) + t248 * t197;
t185 = t204 * mrSges(6,1) - t205 * mrSges(6,3);
t280 = -t204 * mrSges(5,1) - t205 * mrSges(5,2) - t185;
t282 = -mrSges(5,3) - mrSges(6,2);
t134 = m(5) * t153 - t237 * mrSges(5,2) + t282 * t173 - t248 * t196 + t280 * t204 + t276;
t272 = -m(6) * t148 + t237 * mrSges(6,1) + t248 * t198;
t136 = m(5) * t152 + t237 * mrSges(5,1) + t282 * t174 + t248 * t195 + t280 * t205 + t272;
t273 = t283 * t134 - t255 * t136;
t115 = -mrSges(4,1) * t193 + mrSges(4,3) * t159 + Ifges(4,4) * t214 + Ifges(4,2) * t213 - Ifges(4,6) * t241 - pkin(3) * t266 + pkin(7) * t273 + t283 * t129 + t255 * t130 - t235 * t199 - t201 * t278;
t131 = t255 * t134 + t283 * t136;
t200 = Ifges(4,4) * t235 + Ifges(4,2) * t234 - Ifges(4,6) * t278;
t116 = mrSges(4,2) * t193 - mrSges(4,3) * t158 + Ifges(4,1) * t214 + Ifges(4,4) * t213 - Ifges(4,5) * t241 - pkin(7) * t131 - t255 * t129 + t283 * t130 + t234 * t199 + t200 * t278;
t206 = -t234 * mrSges(4,1) + t235 * mrSges(4,2);
t211 = mrSges(4,2) * t278 + t234 * mrSges(4,3);
t127 = m(4) * t158 - t241 * mrSges(4,1) - t214 * mrSges(4,3) - t235 * t206 - t211 * t278 + t131;
t212 = -mrSges(4,1) * t278 - t235 * mrSges(4,3);
t128 = m(4) * t159 + t241 * mrSges(4,2) + t213 * mrSges(4,3) + t234 * t206 + t212 * t278 + t273;
t125 = -t253 * t127 + t254 * t128;
t138 = -m(4) * t193 + t213 * mrSges(4,1) - t214 * mrSges(4,2) + t234 * t211 - t235 * t212 - t266;
t223 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t256 + Ifges(3,2) * t258) * qJD(1);
t224 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t256 + Ifges(3,4) * t258) * qJD(1);
t284 = mrSges(3,1) * t208 - mrSges(3,2) * t209 + Ifges(3,5) * t240 + Ifges(3,6) * t241 + Ifges(3,3) * qJDD(2) + pkin(2) * t138 + qJ(3) * t125 + t254 * t115 + t253 * t116 + (t256 * t223 - t258 * t224) * qJD(1);
t239 = (-mrSges(3,1) * t258 + mrSges(3,2) * t256) * qJD(1);
t243 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t279;
t123 = m(3) * t209 - qJDD(2) * mrSges(3,2) + t241 * mrSges(3,3) - qJD(2) * t243 + t239 * t278 + t125;
t244 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t278;
t137 = m(3) * t208 + qJDD(2) * mrSges(3,1) - t240 * mrSges(3,3) + qJD(2) * t244 - t239 * t279 + t138;
t274 = t258 * t123 - t256 * t137;
t124 = t254 * t127 + t253 * t128;
t222 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t256 + Ifges(3,6) * t258) * qJD(1);
t112 = mrSges(3,2) * t225 - mrSges(3,3) * t208 + Ifges(3,1) * t240 + Ifges(3,4) * t241 + Ifges(3,5) * qJDD(2) - qJ(3) * t124 - qJD(2) * t223 - t253 * t115 + t254 * t116 + t222 * t278;
t267 = mrSges(6,1) * t148 - mrSges(6,3) * t146 - Ifges(6,4) * t174 - Ifges(6,2) * t237 - Ifges(6,6) * t173 + t205 * t175 - t204 * t179;
t264 = mrSges(5,2) * t153 - t204 * t180 - qJ(5) * (-t173 * mrSges(6,2) - t204 * t185 + t276) - pkin(4) * (-t174 * mrSges(6,2) - t205 * t185 + t272) - mrSges(5,1) * t152 - t205 * t178 + Ifges(5,6) * t173 - Ifges(5,5) * t174 - Ifges(5,3) * t237 + t267;
t262 = mrSges(4,1) * t158 - mrSges(4,2) * t159 + Ifges(4,5) * t214 + Ifges(4,6) * t213 + pkin(3) * t131 + t235 * t200 - t234 * t201 - t264;
t114 = -t222 * t279 + Ifges(3,6) * qJDD(2) + (Ifges(4,3) + Ifges(3,2)) * t241 + Ifges(3,4) * t240 + qJD(2) * t224 - mrSges(3,1) * t225 + mrSges(3,3) * t209 - pkin(2) * t124 - t262;
t265 = -m(3) * t225 + t241 * mrSges(3,1) - t240 * mrSges(3,2) - t243 * t279 + t244 * t278 - t124;
t268 = mrSges(2,1) * t245 - mrSges(2,2) * t246 + Ifges(2,3) * qJDD(1) + pkin(1) * t265 + pkin(6) * t274 + t256 * t112 + t258 * t114;
t120 = m(2) * t245 + qJDD(1) * mrSges(2,1) - t261 * mrSges(2,2) + t265;
t119 = t256 * t123 + t258 * t137;
t117 = m(2) * t246 - t261 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t274;
t110 = mrSges(2,1) * g(3) + mrSges(2,3) * t246 + t261 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t119 - t284;
t109 = -mrSges(2,2) * g(3) - mrSges(2,3) * t245 + Ifges(2,5) * qJDD(1) - t261 * Ifges(2,6) - pkin(6) * t119 + t258 * t112 - t256 * t114;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t259 * t109 - t257 * t110 - pkin(5) * (t257 * t117 + t259 * t120), t109, t112, t116, t130, -t204 * t177 + t269; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t257 * t109 + t259 * t110 + pkin(5) * (t259 * t117 - t257 * t120), t110, t114, t115, t129, -t267; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t268, t268, t284, -Ifges(4,3) * t241 + t262, -t264, Ifges(6,5) * t174 + Ifges(6,6) * t237 + Ifges(6,3) * t173 + t205 * t177 - t248 * t179 - t271;];
m_new = t1;
