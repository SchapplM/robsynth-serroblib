% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:38
% EndTime: 2024-09-27 18:42:56
% DurationCPUTime: 12.05s
% Computational Cost: add. (367435->280), mult. (563944->372), div. (0->0), fcn. (403604->12), ass. (0->124)
t258 = sin(pkin(5));
t291 = pkin(8) * t258;
t255 = qJD(1) + qJD(2);
t252 = t255 ^ 2;
t290 = t252 * t258 ^ 2;
t289 = t255 * t258;
t262 = sin(qJ(3));
t288 = t258 * t262;
t267 = cos(qJ(3));
t287 = t258 * t267;
t259 = cos(pkin(5));
t286 = t259 * t262;
t285 = t259 * t267;
t253 = qJDD(1) + qJDD(2);
t284 = qJD(3) * t255;
t230 = (t253 * t262 + t267 * t284) * t258;
t244 = t259 * t253 + qJDD(3);
t245 = t259 * t255 + qJD(3);
t264 = sin(qJ(1));
t269 = cos(qJ(1));
t246 = t264 * g(1) - g(2) * t269;
t238 = qJDD(1) * pkin(1) + t246;
t247 = -g(1) * t269 - g(2) * t264;
t270 = qJD(1) ^ 2;
t240 = -pkin(1) * t270 + t247;
t263 = sin(qJ(2));
t268 = cos(qJ(2));
t221 = t268 * t238 - t240 * t263;
t211 = pkin(2) * t253 + t252 * t291 + t221;
t222 = t263 * t238 + t268 * t240;
t212 = -pkin(2) * t252 + t253 * t291 + t222;
t277 = t211 * t285 - t212 * t262;
t181 = pkin(3) * t244 - pkin(9) * t230 + (pkin(3) * t262 * t290 + (pkin(9) * t245 * t255 - g(3)) * t258) * t267 + t277;
t191 = -g(3) * t288 + t211 * t286 + t267 * t212;
t282 = t255 * t288;
t229 = pkin(3) * t245 - pkin(9) * t282;
t231 = (t253 * t267 - t262 * t284) * t258;
t283 = t267 ^ 2 * t290;
t182 = -pkin(3) * t283 + pkin(9) * t231 - t229 * t245 + t191;
t261 = sin(qJ(4));
t266 = cos(qJ(4));
t168 = t266 * t181 - t182 * t261;
t224 = (-t261 * t262 + t266 * t267) * t289;
t197 = qJD(4) * t224 + t230 * t266 + t231 * t261;
t225 = (t261 * t267 + t262 * t266) * t289;
t241 = qJDD(4) + t244;
t243 = qJD(4) + t245;
t165 = (t224 * t243 - t197) * pkin(10) + (t224 * t225 + t241) * pkin(4) + t168;
t169 = t261 * t181 + t266 * t182;
t196 = -qJD(4) * t225 - t230 * t261 + t231 * t266;
t215 = pkin(4) * t243 - pkin(10) * t225;
t223 = t224 ^ 2;
t166 = -pkin(4) * t223 + pkin(10) * t196 - t215 * t243 + t169;
t260 = sin(qJ(5));
t265 = cos(qJ(5));
t163 = t165 * t265 - t166 * t260;
t204 = t224 * t265 - t225 * t260;
t176 = qJD(5) * t204 + t196 * t260 + t197 * t265;
t205 = t224 * t260 + t225 * t265;
t188 = -mrSges(6,1) * t204 + mrSges(6,2) * t205;
t239 = qJD(5) + t243;
t198 = -mrSges(6,2) * t239 + mrSges(6,3) * t204;
t237 = qJDD(5) + t241;
t160 = m(6) * t163 + mrSges(6,1) * t237 - mrSges(6,3) * t176 - t188 * t205 + t198 * t239;
t164 = t165 * t260 + t166 * t265;
t175 = -qJD(5) * t205 + t196 * t265 - t197 * t260;
t199 = mrSges(6,1) * t239 - mrSges(6,3) * t205;
t161 = m(6) * t164 - mrSges(6,2) * t237 + mrSges(6,3) * t175 + t188 * t204 - t199 * t239;
t152 = t265 * t160 + t260 * t161;
t206 = -mrSges(5,1) * t224 + mrSges(5,2) * t225;
t213 = -mrSges(5,2) * t243 + mrSges(5,3) * t224;
t149 = m(5) * t168 + mrSges(5,1) * t241 - mrSges(5,3) * t197 - t206 * t225 + t213 * t243 + t152;
t214 = mrSges(5,1) * t243 - mrSges(5,3) * t225;
t278 = -t160 * t260 + t265 * t161;
t150 = m(5) * t169 - mrSges(5,2) * t241 + mrSges(5,3) * t196 + t206 * t224 - t214 * t243 + t278;
t145 = t266 * t149 + t261 * t150;
t190 = -g(3) * t287 + t277;
t281 = t255 * t287;
t227 = -mrSges(4,2) * t245 + mrSges(4,3) * t281;
t228 = (-mrSges(4,1) * t267 + mrSges(4,2) * t262) * t289;
t143 = m(4) * t190 + mrSges(4,1) * t244 - mrSges(4,3) * t230 + t227 * t245 - t228 * t282 + t145;
t226 = mrSges(4,1) * t245 - mrSges(4,3) * t282;
t279 = -t149 * t261 + t266 * t150;
t144 = m(4) * t191 - mrSges(4,2) * t244 + mrSges(4,3) * t231 - t226 * t245 + t228 * t281 + t279;
t207 = -g(3) * t259 - t258 * t211;
t189 = -pkin(3) * t231 - pkin(9) * t283 + t229 * t282 + t207;
t171 = -pkin(4) * t196 - pkin(10) * t223 + t215 * t225 + t189;
t276 = m(6) * t171 - mrSges(6,1) * t175 + t176 * mrSges(6,2) - t198 * t204 + t205 * t199;
t272 = m(5) * t189 - mrSges(5,1) * t196 + t197 * mrSges(5,2) - t213 * t224 + t225 * t214 + t276;
t156 = m(4) * t207 - mrSges(4,1) * t231 + mrSges(4,2) * t230 + (t226 * t262 - t227 * t267) * t289 + t272;
t128 = t143 * t285 + t144 * t286 - t258 * t156;
t125 = m(3) * t221 + mrSges(3,1) * t253 - mrSges(3,2) * t252 + t128;
t133 = -t143 * t262 + t267 * t144;
t131 = m(3) * t222 - mrSges(3,1) * t252 - mrSges(3,2) * t253 + t133;
t119 = t268 * t125 + t263 * t131;
t127 = t143 * t287 + t144 * t288 + t259 * t156;
t280 = -t125 * t263 + t268 * t131;
t183 = Ifges(6,5) * t205 + Ifges(6,6) * t204 + Ifges(6,3) * t239;
t185 = Ifges(6,1) * t205 + Ifges(6,4) * t204 + Ifges(6,5) * t239;
t153 = -mrSges(6,1) * t171 + mrSges(6,3) * t164 + Ifges(6,4) * t176 + Ifges(6,2) * t175 + Ifges(6,6) * t237 - t183 * t205 + t185 * t239;
t184 = Ifges(6,4) * t205 + Ifges(6,2) * t204 + Ifges(6,6) * t239;
t154 = mrSges(6,2) * t171 - mrSges(6,3) * t163 + Ifges(6,1) * t176 + Ifges(6,4) * t175 + Ifges(6,5) * t237 + t183 * t204 - t184 * t239;
t200 = Ifges(5,5) * t225 + Ifges(5,6) * t224 + Ifges(5,3) * t243;
t202 = Ifges(5,1) * t225 + Ifges(5,4) * t224 + Ifges(5,5) * t243;
t136 = -mrSges(5,1) * t189 + mrSges(5,3) * t169 + Ifges(5,4) * t197 + Ifges(5,2) * t196 + Ifges(5,6) * t241 - pkin(4) * t276 + pkin(10) * t278 + t265 * t153 + t260 * t154 - t225 * t200 + t243 * t202;
t201 = Ifges(5,4) * t225 + Ifges(5,2) * t224 + Ifges(5,6) * t243;
t137 = mrSges(5,2) * t189 - mrSges(5,3) * t168 + Ifges(5,1) * t197 + Ifges(5,4) * t196 + Ifges(5,5) * t241 - pkin(10) * t152 - t153 * t260 + t154 * t265 + t200 * t224 - t201 * t243;
t216 = Ifges(4,3) * t245 + (Ifges(4,5) * t262 + Ifges(4,6) * t267) * t289;
t218 = Ifges(4,5) * t245 + (Ifges(4,1) * t262 + Ifges(4,4) * t267) * t289;
t121 = -mrSges(4,1) * t207 + mrSges(4,3) * t191 + Ifges(4,4) * t230 + Ifges(4,2) * t231 + Ifges(4,6) * t244 - pkin(3) * t272 + pkin(9) * t279 + t266 * t136 + t261 * t137 - t216 * t282 + t245 * t218;
t217 = Ifges(4,6) * t245 + (Ifges(4,4) * t262 + Ifges(4,2) * t267) * t289;
t123 = mrSges(4,2) * t207 - mrSges(4,3) * t190 + Ifges(4,1) * t230 + Ifges(4,4) * t231 + Ifges(4,5) * t244 - pkin(9) * t145 - t136 * t261 + t137 * t266 + t216 * t281 - t217 * t245;
t274 = mrSges(6,1) * t163 - mrSges(6,2) * t164 + Ifges(6,5) * t176 + Ifges(6,6) * t175 + Ifges(6,3) * t237 + t205 * t184 - t204 * t185;
t271 = mrSges(5,1) * t168 - mrSges(5,2) * t169 + Ifges(5,5) * t197 + Ifges(5,6) * t196 + Ifges(5,3) * t241 + pkin(4) * t152 + t225 * t201 - t224 * t202 + t274;
t135 = (t217 * t262 - t218 * t267) * t289 + Ifges(4,3) * t244 + Ifges(4,5) * t230 + Ifges(4,6) * t231 + mrSges(4,1) * t190 - mrSges(4,2) * t191 + pkin(3) * t145 + t271;
t275 = mrSges(3,1) * t221 - mrSges(3,2) * t222 + Ifges(3,3) * t253 + pkin(2) * t128 + t121 * t287 + t123 * t288 + t133 * t291 + t259 * t135;
t273 = mrSges(2,1) * t246 - mrSges(2,2) * t247 + Ifges(2,3) * qJDD(1) + pkin(1) * t119 + t275;
t117 = m(2) * t247 - mrSges(2,1) * t270 - qJDD(1) * mrSges(2,2) + t280;
t116 = m(2) * t246 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t270 + t119;
t115 = -mrSges(3,2) * g(3) - mrSges(3,3) * t221 + Ifges(3,5) * t253 - Ifges(3,6) * t252 - t121 * t262 + t123 * t267 + (-t127 * t258 - t128 * t259) * pkin(8);
t114 = mrSges(3,1) * g(3) + mrSges(3,3) * t222 + Ifges(3,5) * t252 + Ifges(3,6) * t253 - pkin(2) * t127 - t135 * t258 + (pkin(8) * t133 + t121 * t267 + t123 * t262) * t259;
t113 = -mrSges(2,2) * g(3) - mrSges(2,3) * t246 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t270 - pkin(7) * t119 - t114 * t263 + t115 * t268;
t112 = Ifges(2,6) * qJDD(1) + t270 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t247 + t263 * t115 + t268 * t114 - pkin(1) * (-m(3) * g(3) + t127) + pkin(7) * t280;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t113 - t264 * t112 - pkin(6) * (t116 * t269 + t117 * t264), t113, t115, t123, t137, t154; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t113 + t269 * t112 + pkin(6) * (-t116 * t264 + t117 * t269), t112, t114, t121, t136, t153; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t273, t273, t275, t135, t271, t274;];
m_new = t1;
