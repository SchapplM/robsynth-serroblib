% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:24
% EndTime: 2024-09-27 21:45:34
% DurationCPUTime: 4.77s
% Computational Cost: add. (213135->267), mult. (531542->357), div. (0->0), fcn. (403604->12), ass. (0->119)
t257 = sin(pkin(5));
t288 = pkin(7) * t257;
t268 = qJD(2) ^ 2;
t287 = t257 ^ 2 * t268;
t262 = sin(qJ(3));
t286 = t257 * t262;
t266 = cos(qJ(3));
t285 = t257 * t266;
t259 = cos(pkin(5));
t284 = t259 * t262;
t283 = t259 * t266;
t256 = sin(pkin(10));
t258 = cos(pkin(10));
t241 = t256 * g(1) - t258 * g(2);
t242 = -t258 * g(1) - t256 * g(2);
t263 = sin(qJ(2));
t267 = cos(qJ(2));
t221 = t267 * t241 - t263 * t242;
t212 = qJDD(2) * pkin(2) + t268 * t288 + t221;
t222 = t263 * t241 + t267 * t242;
t213 = -t268 * pkin(2) + qJDD(2) * t288 + t222;
t255 = -g(3) + qJDD(1);
t191 = t212 * t283 - t262 * t213 + t255 * t285;
t281 = qJD(2) * qJD(3);
t231 = (qJDD(2) * t262 + t266 * t281) * t257;
t247 = t259 * qJDD(2) + qJDD(3);
t248 = t259 * qJD(2) + qJD(3);
t282 = qJD(2) * t257;
t278 = t266 * t282;
t182 = (t248 * t278 - t231) * pkin(8) + (t262 * t266 * t287 + t247) * pkin(3) + t191;
t192 = t212 * t284 + t266 * t213 + t255 * t286;
t279 = t262 * t282;
t230 = t248 * pkin(3) - pkin(8) * t279;
t232 = (qJDD(2) * t266 - t262 * t281) * t257;
t280 = t266 ^ 2 * t287;
t183 = -pkin(3) * t280 + t232 * pkin(8) - t248 * t230 + t192;
t261 = sin(qJ(4));
t265 = cos(qJ(4));
t169 = t265 * t182 - t261 * t183;
t225 = (-t261 * t262 + t265 * t266) * t282;
t198 = t225 * qJD(4) + t265 * t231 + t261 * t232;
t226 = (t261 * t266 + t262 * t265) * t282;
t245 = qJDD(4) + t247;
t246 = qJD(4) + t248;
t166 = (t225 * t246 - t198) * pkin(9) + (t225 * t226 + t245) * pkin(4) + t169;
t170 = t261 * t182 + t265 * t183;
t197 = -t226 * qJD(4) - t261 * t231 + t265 * t232;
t216 = t246 * pkin(4) - t226 * pkin(9);
t224 = t225 ^ 2;
t167 = -t224 * pkin(4) + t197 * pkin(9) - t246 * t216 + t170;
t260 = sin(qJ(5));
t264 = cos(qJ(5));
t164 = t264 * t166 - t260 * t167;
t205 = t264 * t225 - t260 * t226;
t177 = t205 * qJD(5) + t260 * t197 + t264 * t198;
t206 = t260 * t225 + t264 * t226;
t188 = -t205 * mrSges(6,1) + t206 * mrSges(6,2);
t240 = qJD(5) + t246;
t199 = -t240 * mrSges(6,2) + t205 * mrSges(6,3);
t239 = qJDD(5) + t245;
t161 = m(6) * t164 + t239 * mrSges(6,1) - t177 * mrSges(6,3) - t206 * t188 + t240 * t199;
t165 = t260 * t166 + t264 * t167;
t176 = -t206 * qJD(5) + t264 * t197 - t260 * t198;
t200 = t240 * mrSges(6,1) - t206 * mrSges(6,3);
t162 = m(6) * t165 - t239 * mrSges(6,2) + t176 * mrSges(6,3) + t205 * t188 - t240 * t200;
t153 = t264 * t161 + t260 * t162;
t208 = -t225 * mrSges(5,1) + t226 * mrSges(5,2);
t214 = -t246 * mrSges(5,2) + t225 * mrSges(5,3);
t150 = m(5) * t169 + t245 * mrSges(5,1) - t198 * mrSges(5,3) - t226 * t208 + t246 * t214 + t153;
t215 = t246 * mrSges(5,1) - t226 * mrSges(5,3);
t275 = -t260 * t161 + t264 * t162;
t151 = m(5) * t170 - t245 * mrSges(5,2) + t197 * mrSges(5,3) + t225 * t208 - t246 * t215 + t275;
t146 = t265 * t150 + t261 * t151;
t228 = -t248 * mrSges(4,2) + mrSges(4,3) * t278;
t229 = (-mrSges(4,1) * t266 + mrSges(4,2) * t262) * t282;
t144 = m(4) * t191 + t247 * mrSges(4,1) - t231 * mrSges(4,3) + t248 * t228 - t229 * t279 + t146;
t227 = t248 * mrSges(4,1) - mrSges(4,3) * t279;
t276 = -t261 * t150 + t265 * t151;
t145 = m(4) * t192 - t247 * mrSges(4,2) + t232 * mrSges(4,3) - t248 * t227 + t229 * t278 + t276;
t207 = -t257 * t212 + t259 * t255;
t190 = -t232 * pkin(3) - pkin(8) * t280 + t230 * t279 + t207;
t172 = -t197 * pkin(4) - t224 * pkin(9) + t226 * t216 + t190;
t274 = m(6) * t172 - t176 * mrSges(6,1) + t177 * mrSges(6,2) - t205 * t199 + t206 * t200;
t270 = m(5) * t190 - t197 * mrSges(5,1) + t198 * mrSges(5,2) - t225 * t214 + t226 * t215 + t274;
t157 = t231 * mrSges(4,2) - t232 * mrSges(4,1) + m(4) * t207 + t270 + (t227 * t262 - t228 * t266) * t282;
t129 = t144 * t283 + t145 * t284 - t257 * t157;
t126 = m(3) * t221 + qJDD(2) * mrSges(3,1) - t268 * mrSges(3,2) + t129;
t134 = -t262 * t144 + t266 * t145;
t132 = m(3) * t222 - t268 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t134;
t120 = t267 * t126 + t263 * t132;
t128 = t144 * t285 + t145 * t286 + t259 * t157;
t277 = -t263 * t126 + t267 * t132;
t184 = Ifges(6,5) * t206 + Ifges(6,6) * t205 + Ifges(6,3) * t240;
t186 = Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t240;
t154 = -mrSges(6,1) * t172 + mrSges(6,3) * t165 + Ifges(6,4) * t177 + Ifges(6,2) * t176 + Ifges(6,6) * t239 - t206 * t184 + t240 * t186;
t185 = Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t240;
t155 = mrSges(6,2) * t172 - mrSges(6,3) * t164 + Ifges(6,1) * t177 + Ifges(6,4) * t176 + Ifges(6,5) * t239 + t205 * t184 - t240 * t185;
t201 = Ifges(5,5) * t226 + Ifges(5,6) * t225 + Ifges(5,3) * t246;
t203 = Ifges(5,1) * t226 + Ifges(5,4) * t225 + Ifges(5,5) * t246;
t137 = -mrSges(5,1) * t190 + mrSges(5,3) * t170 + Ifges(5,4) * t198 + Ifges(5,2) * t197 + Ifges(5,6) * t245 - pkin(4) * t274 + pkin(9) * t275 + t264 * t154 + t260 * t155 - t226 * t201 + t246 * t203;
t202 = Ifges(5,4) * t226 + Ifges(5,2) * t225 + Ifges(5,6) * t246;
t138 = mrSges(5,2) * t190 - mrSges(5,3) * t169 + Ifges(5,1) * t198 + Ifges(5,4) * t197 + Ifges(5,5) * t245 - pkin(9) * t153 - t260 * t154 + t264 * t155 + t225 * t201 - t246 * t202;
t218 = Ifges(4,3) * t248 + (Ifges(4,5) * t262 + Ifges(4,6) * t266) * t282;
t220 = Ifges(4,5) * t248 + (Ifges(4,1) * t262 + Ifges(4,4) * t266) * t282;
t122 = -mrSges(4,1) * t207 + mrSges(4,3) * t192 + Ifges(4,4) * t231 + Ifges(4,2) * t232 + Ifges(4,6) * t247 - pkin(3) * t270 + pkin(8) * t276 + t265 * t137 + t261 * t138 - t218 * t279 + t248 * t220;
t219 = Ifges(4,6) * t248 + (Ifges(4,4) * t262 + Ifges(4,2) * t266) * t282;
t124 = mrSges(4,2) * t207 - mrSges(4,3) * t191 + Ifges(4,1) * t231 + Ifges(4,4) * t232 + Ifges(4,5) * t247 - pkin(8) * t146 - t261 * t137 + t265 * t138 + t218 * t278 - t248 * t219;
t272 = mrSges(6,1) * t164 - mrSges(6,2) * t165 + Ifges(6,5) * t177 + Ifges(6,6) * t176 + Ifges(6,3) * t239 + t206 * t185 - t205 * t186;
t269 = mrSges(5,1) * t169 - mrSges(5,2) * t170 + Ifges(5,5) * t198 + Ifges(5,6) * t197 + Ifges(5,3) * t245 + pkin(4) * t153 + t226 * t202 - t225 * t203 + t272;
t136 = (t219 * t262 - t220 * t266) * t282 + Ifges(4,3) * t247 + Ifges(4,5) * t231 + Ifges(4,6) * t232 + mrSges(4,1) * t191 - mrSges(4,2) * t192 + pkin(3) * t146 + t269;
t273 = mrSges(3,1) * t221 - mrSges(3,2) * t222 + Ifges(3,3) * qJDD(2) + pkin(2) * t129 + t122 * t285 + t124 * t286 + t134 * t288 + t259 * t136;
t271 = mrSges(2,1) * t241 - mrSges(2,2) * t242 + pkin(1) * t120 + t273;
t118 = m(2) * t242 + t277;
t117 = m(2) * t241 + t120;
t116 = mrSges(3,2) * t255 - mrSges(3,3) * t221 + Ifges(3,5) * qJDD(2) - t268 * Ifges(3,6) - t262 * t122 + t266 * t124 + (-t128 * t257 - t129 * t259) * pkin(7);
t115 = -mrSges(3,1) * t255 + mrSges(3,3) * t222 + t268 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t128 - t257 * t136 + (pkin(7) * t134 + t122 * t266 + t124 * t262) * t259;
t114 = mrSges(2,2) * t255 - mrSges(2,3) * t241 - pkin(6) * t120 - t263 * t115 + t267 * t116;
t113 = -mrSges(2,1) * t255 + mrSges(2,3) * t242 + t263 * t116 + t267 * t115 - pkin(1) * (m(3) * t255 + t128) + pkin(6) * t277;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t258 * t114 - t256 * t113 - qJ(1) * (t258 * t117 + t256 * t118), t114, t116, t124, t138, t155; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t256 * t114 + t258 * t113 + qJ(1) * (-t256 * t117 + t258 * t118), t113, t115, t122, t137, t154; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t271, t271, t273, t136, t269, t272;];
m_new = t1;
