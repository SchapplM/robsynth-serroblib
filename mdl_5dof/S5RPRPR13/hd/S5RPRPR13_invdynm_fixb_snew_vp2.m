% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:02
% EndTime: 2019-12-31 18:32:09
% DurationCPUTime: 4.39s
% Computational Cost: add. (44167->292), mult. (106035->350), div. (0->0), fcn. (71659->8), ass. (0->121)
t259 = sin(qJ(1));
t261 = cos(qJ(1));
t240 = -g(1) * t261 - g(2) * t259;
t263 = qJD(1) ^ 2;
t233 = -pkin(1) * t263 + qJDD(1) * qJ(2) + t240;
t255 = sin(pkin(8));
t256 = cos(pkin(8));
t291 = qJD(1) * qJD(2);
t287 = -t256 * g(3) - 0.2e1 * t255 * t291;
t300 = pkin(2) * t256;
t188 = (-pkin(6) * qJDD(1) + t263 * t300 - t233) * t255 + t287;
t216 = -t255 * g(3) + (t233 + 0.2e1 * t291) * t256;
t289 = qJDD(1) * t256;
t250 = t256 ^ 2;
t297 = t250 * t263;
t189 = -pkin(2) * t297 + pkin(6) * t289 + t216;
t258 = sin(qJ(3));
t301 = cos(qJ(3));
t167 = t188 * t301 - t189 * t258;
t168 = t188 * t258 + t189 * t301;
t288 = t256 * t301;
t294 = qJD(1) * t255;
t231 = -qJD(1) * t288 + t258 * t294;
t277 = t255 * t301 + t256 * t258;
t232 = t277 * qJD(1);
t193 = Ifges(4,4) * t232 - Ifges(4,2) * t231 + Ifges(4,6) * qJD(3);
t202 = -mrSges(5,2) * t231 - mrSges(5,3) * t232;
t290 = qJDD(1) * t255;
t292 = t232 * qJD(3);
t213 = -qJDD(1) * t288 + t258 * t290 + t292;
t293 = t231 * qJD(3);
t214 = qJDD(1) * t277 - t293;
t222 = mrSges(5,1) * t231 - qJD(3) * mrSges(5,3);
t224 = pkin(4) * t232 - qJD(3) * pkin(7);
t230 = t231 ^ 2;
t249 = t255 ^ 2;
t239 = t259 * g(1) - t261 * g(2);
t283 = qJDD(2) - t239;
t212 = (-pkin(1) - t300) * qJDD(1) + (-qJ(2) + (-t249 - t250) * pkin(6)) * t263 + t283;
t302 = -2 * qJD(4);
t265 = pkin(3) * t292 + t232 * t302 + (-t214 + t293) * qJ(4) + t212;
t156 = -t230 * pkin(4) - t232 * t224 + (pkin(3) + pkin(7)) * t213 + t265;
t200 = pkin(3) * t231 - qJ(4) * t232;
t262 = qJD(3) ^ 2;
t164 = -qJDD(3) * pkin(3) - t262 * qJ(4) + t200 * t232 + qJDD(4) - t167;
t157 = (t231 * t232 - qJDD(3)) * pkin(7) + (t214 + t293) * pkin(4) + t164;
t257 = sin(qJ(5));
t260 = cos(qJ(5));
t154 = -t156 * t257 + t157 * t260;
t217 = -qJD(3) * t257 + t231 * t260;
t177 = qJD(5) * t217 + qJDD(3) * t260 + t213 * t257;
t218 = qJD(3) * t260 + t231 * t257;
t180 = -mrSges(6,1) * t217 + mrSges(6,2) * t218;
t228 = qJD(5) + t232;
t184 = -mrSges(6,2) * t228 + mrSges(6,3) * t217;
t211 = qJDD(5) + t214;
t150 = m(6) * t154 + mrSges(6,1) * t211 - t177 * mrSges(6,3) - t180 * t218 + t184 * t228;
t155 = t156 * t260 + t157 * t257;
t176 = -qJD(5) * t218 - qJDD(3) * t257 + t213 * t260;
t185 = mrSges(6,1) * t228 - mrSges(6,3) * t218;
t151 = m(6) * t155 - mrSges(6,2) * t211 + t176 * mrSges(6,3) + t180 * t217 - t185 * t228;
t138 = t260 * t150 + t257 * t151;
t274 = -t262 * pkin(3) + qJDD(3) * qJ(4) - t231 * t200 + t168;
t159 = -t213 * pkin(4) - t230 * pkin(7) + ((2 * qJD(4)) + t224) * qJD(3) + t274;
t169 = Ifges(6,5) * t218 + Ifges(6,6) * t217 + Ifges(6,3) * t228;
t171 = Ifges(6,1) * t218 + Ifges(6,4) * t217 + Ifges(6,5) * t228;
t141 = -mrSges(6,1) * t159 + mrSges(6,3) * t155 + Ifges(6,4) * t177 + Ifges(6,2) * t176 + Ifges(6,6) * t211 - t169 * t218 + t171 * t228;
t170 = Ifges(6,4) * t218 + Ifges(6,2) * t217 + Ifges(6,6) * t228;
t142 = mrSges(6,2) * t159 - mrSges(6,3) * t154 + Ifges(6,1) * t177 + Ifges(6,4) * t176 + Ifges(6,5) * t211 + t169 * t217 - t170 * t228;
t162 = qJD(3) * t302 - t274;
t190 = Ifges(5,5) * qJD(3) - Ifges(5,6) * t232 + Ifges(5,3) * t231;
t271 = -mrSges(5,2) * t164 + mrSges(5,3) * t162 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t214 - Ifges(5,5) * t213 + pkin(7) * t138 + t257 * t141 - t142 * t260 + t190 * t232;
t152 = -m(6) * t159 + t176 * mrSges(6,1) - t177 * mrSges(6,2) + t217 * t184 - t185 * t218;
t223 = mrSges(5,1) * t232 + qJD(3) * mrSges(5,2);
t272 = -m(5) * t162 + qJDD(3) * mrSges(5,3) + qJD(3) * t223 - t152;
t275 = -m(5) * t164 - t214 * mrSges(5,1) - t202 * t232 - t138;
t192 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t232 + Ifges(5,6) * t231;
t295 = Ifges(4,1) * t232 - Ifges(4,4) * t231 + Ifges(4,5) * qJD(3) - t192;
t306 = -mrSges(4,2) * t168 + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t222 + t275) + qJ(4) * (-t213 * mrSges(5,1) - t231 * t202 + t272) + mrSges(4,1) * t167 + t232 * t193 - Ifges(4,6) * t213 + Ifges(4,5) * t214 + Ifges(4,3) * qJDD(3) - t271 + t295 * t231;
t201 = mrSges(4,1) * t231 + mrSges(4,2) * t232;
t220 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t231;
t134 = m(4) * t167 - t214 * mrSges(4,3) - t232 * t201 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t220 - t222) * qJD(3) + t275;
t221 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t232;
t145 = m(4) * t168 - qJDD(3) * mrSges(4,2) - qJD(3) * t221 + (-t201 - t202) * t231 + (-mrSges(4,3) - mrSges(5,1)) * t213 + t272;
t130 = t134 * t301 + t145 * t258;
t215 = -t255 * t233 + t287;
t281 = Ifges(3,4) * t255 + Ifges(3,2) * t256;
t282 = Ifges(3,1) * t255 + Ifges(3,4) * t256;
t304 = qJD(1) * t256;
t305 = -mrSges(3,1) * t215 + mrSges(3,2) * t216 - pkin(2) * t130 - (t281 * t294 - t282 * t304) * qJD(1) - t306;
t299 = Ifges(4,4) + Ifges(5,6);
t298 = mrSges(3,2) * t255;
t139 = -t150 * t257 + t151 * t260;
t194 = Ifges(5,1) * qJD(3) - Ifges(5,4) * t232 + Ifges(5,5) * t231;
t296 = -Ifges(4,5) * t232 + Ifges(4,6) * t231 - Ifges(4,3) * qJD(3) - t194;
t278 = mrSges(3,3) * qJDD(1) + t263 * (-mrSges(3,1) * t256 + t298);
t128 = m(3) * t215 - t255 * t278 + t130;
t284 = -t258 * t134 + t145 * t301;
t129 = m(3) * t216 + t256 * t278 + t284;
t285 = -t128 * t255 + t129 * t256;
t280 = Ifges(3,5) * t255 + Ifges(3,6) * t256;
t161 = t213 * pkin(3) + t265;
t135 = m(5) * t161 - t213 * mrSges(5,2) - t214 * mrSges(5,3) - t222 * t231 - t223 * t232 + t139;
t270 = -mrSges(5,1) * t162 + mrSges(5,2) * t161 - pkin(4) * t152 - pkin(7) * t139 - t260 * t141 - t257 * t142;
t122 = -mrSges(4,1) * t212 + mrSges(4,3) * t168 - pkin(3) * t135 + t296 * t232 + t299 * t214 + (-Ifges(4,2) - Ifges(5,3)) * t213 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + t295 * qJD(3) + t270;
t273 = mrSges(6,1) * t154 - mrSges(6,2) * t155 + Ifges(6,5) * t177 + Ifges(6,6) * t176 + Ifges(6,3) * t211 + t170 * t218 - t217 * t171;
t268 = mrSges(5,1) * t164 - mrSges(5,3) * t161 + pkin(4) * t138 + t273;
t126 = mrSges(4,2) * t212 - mrSges(4,3) * t167 - qJ(4) * t135 + (-t193 + t190) * qJD(3) + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t299 * t213 + (Ifges(4,1) + Ifges(5,2)) * t214 + t268 + t296 * t231;
t229 = -qJDD(1) * pkin(1) - t263 * qJ(2) + t283;
t235 = t280 * qJD(1);
t269 = m(4) * t212 + t213 * mrSges(4,1) + t214 * mrSges(4,2) + t220 * t231 + t232 * t221 + t135;
t118 = -mrSges(3,1) * t229 + mrSges(3,3) * t216 - pkin(2) * t269 + pkin(6) * t284 + qJDD(1) * t281 + t122 * t301 + t258 * t126 - t235 * t294;
t121 = mrSges(3,2) * t229 - mrSges(3,3) * t215 - pkin(6) * t130 + qJDD(1) * t282 - t258 * t122 + t126 * t301 + t235 * t304;
t267 = -m(3) * t229 + mrSges(3,1) * t289 - t269 + (t249 * t263 + t297) * mrSges(3,3);
t276 = -mrSges(2,2) * t240 + qJ(2) * t285 + t256 * t118 + t255 * t121 + pkin(1) * (-mrSges(3,2) * t290 + t267) + mrSges(2,1) * t239 + Ifges(2,3) * qJDD(1);
t131 = t267 - t263 * mrSges(2,2) + m(2) * t239 + (mrSges(2,1) - t298) * qJDD(1);
t125 = t128 * t256 + t129 * t255;
t123 = m(2) * t240 - mrSges(2,1) * t263 - qJDD(1) * mrSges(2,2) + t285;
t119 = t263 * Ifges(2,5) + mrSges(2,3) * t240 + mrSges(2,1) * g(3) - pkin(1) * t125 + (Ifges(2,6) - t280) * qJDD(1) + t305;
t116 = -mrSges(2,2) * g(3) - mrSges(2,3) * t239 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t263 - qJ(2) * t125 - t118 * t255 + t121 * t256;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t261 * t116 - t259 * t119 - pkin(5) * (t123 * t259 + t131 * t261), t116, t121, t126, -t231 * t192 - t271, t142; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t259 * t116 + t261 * t119 + pkin(5) * (t123 * t261 - t131 * t259), t119, t118, t122, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t214 + Ifges(5,6) * t213 - qJD(3) * t190 + t231 * t194 - t268, t141; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t276, t276, qJDD(1) * t280 - t305, t306, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t214 + Ifges(5,3) * t213 + qJD(3) * t192 + t232 * t194 - t270, t273;];
m_new = t1;
