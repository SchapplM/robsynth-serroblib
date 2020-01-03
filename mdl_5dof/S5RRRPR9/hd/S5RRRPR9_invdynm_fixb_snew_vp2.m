% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:27
% EndTime: 2019-12-31 21:22:47
% DurationCPUTime: 11.10s
% Computational Cost: add. (190107->309), mult. (396097->392), div. (0->0), fcn. (274792->10), ass. (0->123)
t262 = sin(qJ(1));
t266 = cos(qJ(1));
t249 = t262 * g(1) - t266 * g(2);
t268 = qJD(1) ^ 2;
t232 = -qJDD(1) * pkin(1) - t268 * pkin(6) - t249;
t261 = sin(qJ(2));
t265 = cos(qJ(2));
t282 = qJD(1) * qJD(2);
t281 = t265 * t282;
t244 = t261 * qJDD(1) + t281;
t253 = t261 * t282;
t245 = t265 * qJDD(1) - t253;
t199 = (-t244 - t281) * pkin(7) + (-t245 + t253) * pkin(2) + t232;
t250 = -t266 * g(1) - t262 * g(2);
t233 = -t268 * pkin(1) + qJDD(1) * pkin(6) + t250;
t224 = -t261 * g(3) + t265 * t233;
t243 = (-pkin(2) * t265 - pkin(7) * t261) * qJD(1);
t267 = qJD(2) ^ 2;
t283 = t265 * qJD(1);
t202 = -t267 * pkin(2) + qJDD(2) * pkin(7) + t243 * t283 + t224;
t260 = sin(qJ(3));
t264 = cos(qJ(3));
t181 = t264 * t199 - t202 * t260;
t284 = qJD(1) * t261;
t240 = t264 * qJD(2) - t260 * t284;
t215 = t240 * qJD(3) + t260 * qJDD(2) + t264 * t244;
t239 = qJDD(3) - t245;
t241 = t260 * qJD(2) + t264 * t284;
t252 = qJD(3) - t283;
t170 = (t240 * t252 - t215) * qJ(4) + (t240 * t241 + t239) * pkin(3) + t181;
t182 = t260 * t199 + t264 * t202;
t214 = -t241 * qJD(3) + t264 * qJDD(2) - t260 * t244;
t221 = pkin(3) * t252 - qJ(4) * t241;
t238 = t240 ^ 2;
t172 = -pkin(3) * t238 + qJ(4) * t214 - t221 * t252 + t182;
t257 = sin(pkin(9));
t258 = cos(pkin(9));
t218 = t240 * t257 + t241 * t258;
t157 = -0.2e1 * qJD(4) * t218 + t258 * t170 - t172 * t257;
t192 = t214 * t257 + t215 * t258;
t217 = t240 * t258 - t241 * t257;
t154 = (t217 * t252 - t192) * pkin(8) + (t217 * t218 + t239) * pkin(4) + t157;
t158 = 0.2e1 * qJD(4) * t217 + t257 * t170 + t258 * t172;
t191 = t214 * t258 - t215 * t257;
t205 = pkin(4) * t252 - pkin(8) * t218;
t216 = t217 ^ 2;
t155 = -pkin(4) * t216 + pkin(8) * t191 - t205 * t252 + t158;
t259 = sin(qJ(5));
t263 = cos(qJ(5));
t153 = t259 * t154 + t155 * t263;
t223 = -t265 * g(3) - t261 * t233;
t201 = -qJDD(2) * pkin(2) - t267 * pkin(7) + t243 * t284 - t223;
t179 = -pkin(3) * t214 - qJ(4) * t238 + t241 * t221 + qJDD(4) + t201;
t160 = -pkin(4) * t191 - pkin(8) * t216 + t205 * t218 + t179;
t195 = t259 * t217 + t218 * t263;
t165 = -t195 * qJD(5) + t191 * t263 - t259 * t192;
t194 = t217 * t263 - t259 * t218;
t166 = t194 * qJD(5) + t259 * t191 + t192 * t263;
t251 = qJD(5) + t252;
t173 = Ifges(6,5) * t195 + Ifges(6,6) * t194 + Ifges(6,3) * t251;
t175 = Ifges(6,1) * t195 + Ifges(6,4) * t194 + Ifges(6,5) * t251;
t235 = qJDD(5) + t239;
t142 = -mrSges(6,1) * t160 + mrSges(6,3) * t153 + Ifges(6,4) * t166 + Ifges(6,2) * t165 + Ifges(6,6) * t235 - t173 * t195 + t175 * t251;
t152 = t154 * t263 - t259 * t155;
t174 = Ifges(6,4) * t195 + Ifges(6,2) * t194 + Ifges(6,6) * t251;
t143 = mrSges(6,2) * t160 - mrSges(6,3) * t152 + Ifges(6,1) * t166 + Ifges(6,4) * t165 + Ifges(6,5) * t235 + t173 * t194 - t174 * t251;
t188 = Ifges(5,5) * t218 + Ifges(5,6) * t217 + Ifges(5,3) * t252;
t190 = Ifges(5,1) * t218 + Ifges(5,4) * t217 + Ifges(5,5) * t252;
t183 = -mrSges(6,2) * t251 + mrSges(6,3) * t194;
t184 = mrSges(6,1) * t251 - mrSges(6,3) * t195;
t277 = m(6) * t160 - t165 * mrSges(6,1) + t166 * mrSges(6,2) - t194 * t183 + t195 * t184;
t178 = -mrSges(6,1) * t194 + mrSges(6,2) * t195;
t148 = m(6) * t152 + mrSges(6,1) * t235 - mrSges(6,3) * t166 - t178 * t195 + t183 * t251;
t149 = m(6) * t153 - mrSges(6,2) * t235 + mrSges(6,3) * t165 + t178 * t194 - t184 * t251;
t278 = -t148 * t259 + t263 * t149;
t129 = -mrSges(5,1) * t179 + mrSges(5,3) * t158 + Ifges(5,4) * t192 + Ifges(5,2) * t191 + Ifges(5,6) * t239 - pkin(4) * t277 + pkin(8) * t278 + t263 * t142 + t259 * t143 - t218 * t188 + t252 * t190;
t141 = t263 * t148 + t259 * t149;
t189 = Ifges(5,4) * t218 + Ifges(5,2) * t217 + Ifges(5,6) * t252;
t130 = mrSges(5,2) * t179 - mrSges(5,3) * t157 + Ifges(5,1) * t192 + Ifges(5,4) * t191 + Ifges(5,5) * t239 - pkin(8) * t141 - t259 * t142 + t143 * t263 + t217 * t188 - t252 * t189;
t206 = Ifges(4,5) * t241 + Ifges(4,6) * t240 + Ifges(4,3) * t252;
t208 = Ifges(4,1) * t241 + Ifges(4,4) * t240 + Ifges(4,5) * t252;
t203 = -mrSges(5,2) * t252 + mrSges(5,3) * t217;
t204 = mrSges(5,1) * t252 - mrSges(5,3) * t218;
t273 = m(5) * t179 - t191 * mrSges(5,1) + t192 * mrSges(5,2) - t217 * t203 + t218 * t204 + t277;
t196 = -mrSges(5,1) * t217 + mrSges(5,2) * t218;
t138 = m(5) * t157 + mrSges(5,1) * t239 - mrSges(5,3) * t192 - t196 * t218 + t203 * t252 + t141;
t139 = m(5) * t158 - mrSges(5,2) * t239 + mrSges(5,3) * t191 + t196 * t217 - t204 * t252 + t278;
t279 = -t138 * t257 + t258 * t139;
t118 = -mrSges(4,1) * t201 + mrSges(4,3) * t182 + Ifges(4,4) * t215 + Ifges(4,2) * t214 + Ifges(4,6) * t239 - pkin(3) * t273 + qJ(4) * t279 + t258 * t129 + t257 * t130 - t241 * t206 + t252 * t208;
t134 = t258 * t138 + t257 * t139;
t207 = Ifges(4,4) * t241 + Ifges(4,2) * t240 + Ifges(4,6) * t252;
t119 = mrSges(4,2) * t201 - mrSges(4,3) * t181 + Ifges(4,1) * t215 + Ifges(4,4) * t214 + Ifges(4,5) * t239 - qJ(4) * t134 - t129 * t257 + t130 * t258 + t206 * t240 - t207 * t252;
t219 = -mrSges(4,1) * t240 + mrSges(4,2) * t241;
t220 = -mrSges(4,2) * t252 + mrSges(4,3) * t240;
t132 = m(4) * t181 + mrSges(4,1) * t239 - mrSges(4,3) * t215 - t219 * t241 + t220 * t252 + t134;
t222 = mrSges(4,1) * t252 - mrSges(4,3) * t241;
t133 = m(4) * t182 - mrSges(4,2) * t239 + mrSges(4,3) * t214 + t219 * t240 - t222 * t252 + t279;
t128 = -t260 * t132 + t264 * t133;
t150 = -m(4) * t201 + t214 * mrSges(4,1) - t215 * mrSges(4,2) + t240 * t220 - t241 * t222 - t273;
t230 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t261 + Ifges(3,2) * t265) * qJD(1);
t231 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t261 + Ifges(3,4) * t265) * qJD(1);
t285 = mrSges(3,1) * t223 - mrSges(3,2) * t224 + Ifges(3,5) * t244 + Ifges(3,6) * t245 + Ifges(3,3) * qJDD(2) + pkin(2) * t150 + pkin(7) * t128 + t264 * t118 + t260 * t119 + (t261 * t230 - t265 * t231) * qJD(1);
t242 = (-mrSges(3,1) * t265 + mrSges(3,2) * t261) * qJD(1);
t247 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t284;
t126 = m(3) * t224 - qJDD(2) * mrSges(3,2) + t245 * mrSges(3,3) - qJD(2) * t247 + t242 * t283 + t128;
t248 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t283;
t144 = m(3) * t223 + qJDD(2) * mrSges(3,1) - t244 * mrSges(3,3) + qJD(2) * t248 - t242 * t284 + t150;
t280 = t265 * t126 - t144 * t261;
t127 = t264 * t132 + t260 * t133;
t229 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t261 + Ifges(3,6) * t265) * qJD(1);
t115 = mrSges(3,2) * t232 - mrSges(3,3) * t223 + Ifges(3,1) * t244 + Ifges(3,4) * t245 + Ifges(3,5) * qJDD(2) - pkin(7) * t127 - qJD(2) * t230 - t260 * t118 + t264 * t119 + t229 * t283;
t274 = -mrSges(6,1) * t152 + mrSges(6,2) * t153 - Ifges(6,5) * t166 - Ifges(6,6) * t165 - Ifges(6,3) * t235 - t195 * t174 + t194 * t175;
t271 = -mrSges(5,1) * t157 + mrSges(5,2) * t158 - Ifges(5,5) * t192 - Ifges(5,6) * t191 - Ifges(5,3) * t239 - pkin(4) * t141 - t218 * t189 + t217 * t190 + t274;
t269 = mrSges(4,1) * t181 - mrSges(4,2) * t182 + Ifges(4,5) * t215 + Ifges(4,6) * t214 + Ifges(4,3) * t239 + pkin(3) * t134 + t241 * t207 - t240 * t208 - t271;
t117 = -mrSges(3,1) * t232 + mrSges(3,3) * t224 + Ifges(3,4) * t244 + Ifges(3,2) * t245 + Ifges(3,6) * qJDD(2) - pkin(2) * t127 + qJD(2) * t231 - t229 * t284 - t269;
t272 = -m(3) * t232 + t245 * mrSges(3,1) - t244 * mrSges(3,2) - t247 * t284 + t248 * t283 - t127;
t275 = mrSges(2,1) * t249 - mrSges(2,2) * t250 + Ifges(2,3) * qJDD(1) + pkin(1) * t272 + pkin(6) * t280 + t261 * t115 + t265 * t117;
t123 = m(2) * t249 + qJDD(1) * mrSges(2,1) - t268 * mrSges(2,2) + t272;
t122 = t261 * t126 + t265 * t144;
t120 = m(2) * t250 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t280;
t113 = mrSges(2,1) * g(3) + mrSges(2,3) * t250 + t268 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t122 - t285;
t112 = -mrSges(2,2) * g(3) - mrSges(2,3) * t249 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - pkin(6) * t122 + t265 * t115 - t261 * t117;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t266 * t112 - t262 * t113 - pkin(5) * (t262 * t120 + t266 * t123), t112, t115, t119, t130, t143; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t262 * t112 + t266 * t113 + pkin(5) * (t266 * t120 - t262 * t123), t113, t117, t118, t129, t142; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t275, t275, t285, t269, -t271, -t274;];
m_new = t1;
