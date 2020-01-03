% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR13
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:20
% EndTime: 2019-12-31 20:32:38
% DurationCPUTime: 10.67s
% Computational Cost: add. (177951->310), mult. (387167->392), div. (0->0), fcn. (269191->10), ass. (0->123)
t262 = sin(qJ(1));
t266 = cos(qJ(1));
t249 = t262 * g(1) - t266 * g(2);
t268 = qJD(1) ^ 2;
t232 = -qJDD(1) * pkin(1) - t268 * pkin(6) - t249;
t261 = sin(qJ(2));
t265 = cos(qJ(2));
t282 = qJD(1) * qJD(2);
t281 = t265 * t282;
t244 = qJDD(1) * t261 + t281;
t253 = t261 * t282;
t245 = qJDD(1) * t265 - t253;
t199 = (-t244 - t281) * qJ(3) + (-t245 + t253) * pkin(2) + t232;
t250 = -g(1) * t266 - g(2) * t262;
t233 = -pkin(1) * t268 + qJDD(1) * pkin(6) + t250;
t216 = -g(3) * t261 + t265 * t233;
t242 = (-pkin(2) * t265 - qJ(3) * t261) * qJD(1);
t267 = qJD(2) ^ 2;
t283 = t265 * qJD(1);
t202 = -pkin(2) * t267 + qJDD(2) * qJ(3) + t242 * t283 + t216;
t257 = sin(pkin(9));
t258 = cos(pkin(9));
t284 = qJD(1) * t261;
t239 = qJD(2) * t257 + t258 * t284;
t179 = -0.2e1 * qJD(3) * t239 + t258 * t199 - t257 * t202;
t221 = qJDD(2) * t257 + t244 * t258;
t238 = qJD(2) * t258 - t257 * t284;
t170 = (-t238 * t283 - t221) * pkin(7) + (t238 * t239 - t245) * pkin(3) + t179;
t180 = 0.2e1 * qJD(3) * t238 + t257 * t199 + t258 * t202;
t220 = qJDD(2) * t258 - t244 * t257;
t222 = -pkin(3) * t283 - pkin(7) * t239;
t237 = t238 ^ 2;
t172 = -pkin(3) * t237 + pkin(7) * t220 + t222 * t283 + t180;
t260 = sin(qJ(4));
t264 = cos(qJ(4));
t157 = t170 * t264 - t260 * t172;
t212 = t238 * t264 - t239 * t260;
t189 = qJD(4) * t212 + t220 * t260 + t221 * t264;
t213 = t238 * t260 + t239 * t264;
t241 = qJDD(4) - t245;
t252 = qJD(4) - t283;
t154 = (t212 * t252 - t189) * pkin(8) + (t212 * t213 + t241) * pkin(4) + t157;
t158 = t170 * t260 + t172 * t264;
t188 = -qJD(4) * t213 + t220 * t264 - t221 * t260;
t205 = pkin(4) * t252 - pkin(8) * t213;
t211 = t212 ^ 2;
t155 = -pkin(4) * t211 + pkin(8) * t188 - t205 * t252 + t158;
t259 = sin(qJ(5));
t263 = cos(qJ(5));
t153 = t154 * t259 + t155 * t263;
t215 = -t265 * g(3) - t261 * t233;
t201 = -qJDD(2) * pkin(2) - qJ(3) * t267 + t242 * t284 + qJDD(3) - t215;
t181 = -pkin(3) * t220 - pkin(7) * t237 + t239 * t222 + t201;
t160 = -pkin(4) * t188 - pkin(8) * t211 + t205 * t213 + t181;
t195 = t212 * t259 + t213 * t263;
t165 = -qJD(5) * t195 + t188 * t263 - t189 * t259;
t194 = t212 * t263 - t213 * t259;
t166 = qJD(5) * t194 + t188 * t259 + t189 * t263;
t251 = qJD(5) + t252;
t173 = Ifges(6,5) * t195 + Ifges(6,6) * t194 + Ifges(6,3) * t251;
t175 = Ifges(6,1) * t195 + Ifges(6,4) * t194 + Ifges(6,5) * t251;
t235 = qJDD(5) + t241;
t142 = -mrSges(6,1) * t160 + mrSges(6,3) * t153 + Ifges(6,4) * t166 + Ifges(6,2) * t165 + Ifges(6,6) * t235 - t173 * t195 + t175 * t251;
t152 = t154 * t263 - t155 * t259;
t174 = Ifges(6,4) * t195 + Ifges(6,2) * t194 + Ifges(6,6) * t251;
t143 = mrSges(6,2) * t160 - mrSges(6,3) * t152 + Ifges(6,1) * t166 + Ifges(6,4) * t165 + Ifges(6,5) * t235 + t173 * t194 - t174 * t251;
t190 = Ifges(5,5) * t213 + Ifges(5,6) * t212 + Ifges(5,3) * t252;
t192 = Ifges(5,1) * t213 + Ifges(5,4) * t212 + Ifges(5,5) * t252;
t183 = -mrSges(6,2) * t251 + mrSges(6,3) * t194;
t184 = mrSges(6,1) * t251 - mrSges(6,3) * t195;
t277 = m(6) * t160 - mrSges(6,1) * t165 + mrSges(6,2) * t166 - t183 * t194 + t184 * t195;
t178 = -mrSges(6,1) * t194 + mrSges(6,2) * t195;
t148 = m(6) * t152 + mrSges(6,1) * t235 - mrSges(6,3) * t166 - t178 * t195 + t183 * t251;
t149 = m(6) * t153 - mrSges(6,2) * t235 + mrSges(6,3) * t165 + t178 * t194 - t184 * t251;
t278 = -t148 * t259 + t149 * t263;
t129 = -mrSges(5,1) * t181 + mrSges(5,3) * t158 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * t241 - pkin(4) * t277 + pkin(8) * t278 + t263 * t142 + t259 * t143 - t213 * t190 + t252 * t192;
t141 = t148 * t263 + t149 * t259;
t191 = Ifges(5,4) * t213 + Ifges(5,2) * t212 + Ifges(5,6) * t252;
t130 = mrSges(5,2) * t181 - mrSges(5,3) * t157 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * t241 - pkin(8) * t141 - t142 * t259 + t143 * t263 + t190 * t212 - t191 * t252;
t206 = Ifges(4,5) * t239 + Ifges(4,6) * t238 - Ifges(4,3) * t283;
t208 = Ifges(4,1) * t239 + Ifges(4,4) * t238 - Ifges(4,5) * t283;
t203 = -mrSges(5,2) * t252 + mrSges(5,3) * t212;
t204 = mrSges(5,1) * t252 - mrSges(5,3) * t213;
t273 = m(5) * t181 - t188 * mrSges(5,1) + mrSges(5,2) * t189 - t212 * t203 + t204 * t213 + t277;
t196 = -mrSges(5,1) * t212 + mrSges(5,2) * t213;
t138 = m(5) * t157 + mrSges(5,1) * t241 - mrSges(5,3) * t189 - t196 * t213 + t203 * t252 + t141;
t139 = m(5) * t158 - mrSges(5,2) * t241 + mrSges(5,3) * t188 + t196 * t212 - t204 * t252 + t278;
t279 = -t138 * t260 + t139 * t264;
t118 = -mrSges(4,1) * t201 + mrSges(4,3) * t180 + Ifges(4,4) * t221 + Ifges(4,2) * t220 - Ifges(4,6) * t245 - pkin(3) * t273 + pkin(7) * t279 + t264 * t129 + t260 * t130 - t239 * t206 - t208 * t283;
t134 = t138 * t264 + t139 * t260;
t207 = Ifges(4,4) * t239 + Ifges(4,2) * t238 - Ifges(4,6) * t283;
t119 = mrSges(4,2) * t201 - mrSges(4,3) * t179 + Ifges(4,1) * t221 + Ifges(4,4) * t220 - Ifges(4,5) * t245 - pkin(7) * t134 - t129 * t260 + t130 * t264 + t206 * t238 + t207 * t283;
t214 = -mrSges(4,1) * t238 + mrSges(4,2) * t239;
t218 = mrSges(4,2) * t283 + mrSges(4,3) * t238;
t132 = m(4) * t179 - mrSges(4,1) * t245 - mrSges(4,3) * t221 - t214 * t239 - t218 * t283 + t134;
t219 = -mrSges(4,1) * t283 - mrSges(4,3) * t239;
t133 = m(4) * t180 + mrSges(4,2) * t245 + mrSges(4,3) * t220 + t214 * t238 + t219 * t283 + t279;
t128 = -t132 * t257 + t133 * t258;
t150 = -m(4) * t201 + t220 * mrSges(4,1) - mrSges(4,2) * t221 + t238 * t218 - t219 * t239 - t273;
t230 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t261 + Ifges(3,2) * t265) * qJD(1);
t231 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t261 + Ifges(3,4) * t265) * qJD(1);
t285 = mrSges(3,1) * t215 - mrSges(3,2) * t216 + Ifges(3,5) * t244 + Ifges(3,6) * t245 + Ifges(3,3) * qJDD(2) + pkin(2) * t150 + qJ(3) * t128 + t258 * t118 + t257 * t119 + (t230 * t261 - t231 * t265) * qJD(1);
t243 = (-mrSges(3,1) * t265 + mrSges(3,2) * t261) * qJD(1);
t247 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t284;
t126 = m(3) * t216 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t245 - qJD(2) * t247 + t243 * t283 + t128;
t248 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t283;
t144 = m(3) * t215 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t244 + qJD(2) * t248 - t243 * t284 + t150;
t280 = t126 * t265 - t144 * t261;
t127 = t132 * t258 + t133 * t257;
t229 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t261 + Ifges(3,6) * t265) * qJD(1);
t115 = mrSges(3,2) * t232 - mrSges(3,3) * t215 + Ifges(3,1) * t244 + Ifges(3,4) * t245 + Ifges(3,5) * qJDD(2) - qJ(3) * t127 - qJD(2) * t230 - t118 * t257 + t119 * t258 + t229 * t283;
t274 = -mrSges(6,1) * t152 + mrSges(6,2) * t153 - Ifges(6,5) * t166 - Ifges(6,6) * t165 - Ifges(6,3) * t235 - t174 * t195 + t194 * t175;
t271 = -mrSges(5,1) * t157 + mrSges(5,2) * t158 - Ifges(5,5) * t189 - Ifges(5,6) * t188 - Ifges(5,3) * t241 - pkin(4) * t141 - t191 * t213 + t212 * t192 + t274;
t269 = mrSges(4,1) * t179 - mrSges(4,2) * t180 + Ifges(4,5) * t221 + Ifges(4,6) * t220 + pkin(3) * t134 + t239 * t207 - t238 * t208 - t271;
t117 = -t269 + Ifges(3,6) * qJDD(2) - t229 * t284 + (Ifges(3,2) + Ifges(4,3)) * t245 + Ifges(3,4) * t244 + qJD(2) * t231 - mrSges(3,1) * t232 + mrSges(3,3) * t216 - pkin(2) * t127;
t272 = -m(3) * t232 + t245 * mrSges(3,1) - mrSges(3,2) * t244 - t247 * t284 + t248 * t283 - t127;
t275 = mrSges(2,1) * t249 - mrSges(2,2) * t250 + Ifges(2,3) * qJDD(1) + pkin(1) * t272 + pkin(6) * t280 + t115 * t261 + t117 * t265;
t123 = m(2) * t249 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t268 + t272;
t122 = t126 * t261 + t144 * t265;
t120 = m(2) * t250 - mrSges(2,1) * t268 - qJDD(1) * mrSges(2,2) + t280;
t113 = mrSges(2,1) * g(3) + mrSges(2,3) * t250 + t268 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t122 - t285;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t249 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t268 - pkin(6) * t122 + t115 * t265 - t117 * t261;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t266 * t112 - t262 * t113 - pkin(5) * (t120 * t262 + t123 * t266), t112, t115, t119, t130, t143; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t262 * t112 + t266 * t113 + pkin(5) * (t120 * t266 - t123 * t262), t113, t117, t118, t129, t142; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t275, t275, t285, -Ifges(4,3) * t245 + t269, -t271, -t274;];
m_new = t1;
