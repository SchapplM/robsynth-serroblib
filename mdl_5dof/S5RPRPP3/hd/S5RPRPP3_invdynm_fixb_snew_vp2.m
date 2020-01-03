% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:23
% EndTime: 2019-12-31 18:12:27
% DurationCPUTime: 2.73s
% Computational Cost: add. (22211->291), mult. (53450->332), div. (0->0), fcn. (33636->6), ass. (0->113)
t258 = sin(qJ(3));
t257 = cos(pkin(7));
t306 = cos(qJ(3));
t289 = t257 * t306;
t256 = sin(pkin(7));
t296 = qJD(1) * t256;
t227 = -qJD(1) * t289 + t258 * t296;
t277 = t306 * t256 + t257 * t258;
t228 = t277 * qJD(1);
t187 = -mrSges(6,2) * t228 + mrSges(6,3) * t227;
t294 = t227 * qJD(3);
t209 = t277 * qJDD(1) - t294;
t259 = sin(qJ(1));
t260 = cos(qJ(1));
t236 = -g(1) * t260 - g(2) * t259;
t262 = qJD(1) ^ 2;
t229 = -pkin(1) * t262 + qJDD(1) * qJ(2) + t236;
t293 = qJD(1) * qJD(2);
t288 = -g(3) * t257 - 0.2e1 * t256 * t293;
t304 = pkin(2) * t257;
t170 = (-pkin(6) * qJDD(1) + t262 * t304 - t229) * t256 + t288;
t211 = -g(3) * t256 + (t229 + 0.2e1 * t293) * t257;
t291 = qJDD(1) * t257;
t246 = t257 ^ 2;
t300 = t246 * t262;
t171 = -pkin(2) * t300 + pkin(6) * t291 + t211;
t159 = t306 * t170 - t258 * t171;
t188 = pkin(3) * t227 - qJ(4) * t228;
t261 = qJD(3) ^ 2;
t157 = -qJDD(3) * pkin(3) - t261 * qJ(4) + t228 * t188 + qJDD(4) - t159;
t149 = -0.2e1 * qJD(5) * qJD(3) + (t227 * t228 - qJDD(3)) * qJ(5) + (t209 + t294) * pkin(4) + t157;
t220 = -mrSges(6,1) * t227 + qJD(3) * mrSges(6,2);
t284 = -m(6) * t149 + qJDD(3) * mrSges(6,3) + qJD(3) * t220;
t143 = mrSges(6,1) * t209 + t187 * t228 - t284;
t160 = t258 * t170 + t306 * t171;
t178 = Ifges(4,4) * t228 - Ifges(4,2) * t227 + Ifges(4,6) * qJD(3);
t292 = qJDD(1) * t256;
t295 = qJD(3) * t228;
t208 = -qJDD(1) * t289 + t258 * t292 + t295;
t219 = mrSges(5,1) * t227 - qJD(3) * mrSges(5,3);
t272 = -pkin(3) * t261 + qJDD(3) * qJ(4) - t188 * t227 + t160;
t307 = -2 * qJD(4);
t155 = qJD(3) * t307 - t272;
t174 = Ifges(5,5) * qJD(3) - Ifges(5,6) * t228 + Ifges(5,3) * t227;
t217 = pkin(4) * t228 - qJD(3) * qJ(5);
t226 = t227 ^ 2;
t152 = -pkin(4) * t208 - qJ(5) * t226 + qJDD(5) + ((2 * qJD(4)) + t217) * qJD(3) + t272;
t173 = Ifges(6,5) * qJD(3) + Ifges(6,6) * t227 + Ifges(6,3) * t228;
t176 = Ifges(6,4) * qJD(3) + Ifges(6,2) * t227 + Ifges(6,6) * t228;
t271 = -mrSges(6,2) * t152 + mrSges(6,3) * t149 - Ifges(6,1) * qJDD(3) - Ifges(6,4) * t208 - Ifges(6,5) * t209 - t227 * t173 + t228 * t176;
t267 = -mrSges(5,2) * t157 + mrSges(5,3) * t155 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t209 - Ifges(5,5) * t208 + qJ(5) * t143 + t228 * t174 + t271;
t221 = mrSges(5,1) * t228 + qJD(3) * mrSges(5,2);
t218 = mrSges(6,1) * t228 - qJD(3) * mrSges(6,3);
t290 = m(6) * t152 + qJDD(3) * mrSges(6,2) + qJD(3) * t218;
t276 = -m(5) * t155 + qJDD(3) * mrSges(5,3) + qJD(3) * t221 + t290;
t190 = -mrSges(5,2) * t227 - mrSges(5,3) * t228;
t297 = -t187 - t190;
t177 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t228 + Ifges(5,6) * t227;
t298 = Ifges(4,1) * t228 - Ifges(4,4) * t227 + Ifges(4,5) * qJD(3) - t177;
t309 = -m(5) * t157 - t209 * mrSges(5,1) - t228 * t190;
t312 = -mrSges(4,2) * t160 + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t219 - t143 + t309) + qJ(4) * (t297 * t227 + (-mrSges(5,1) - mrSges(6,1)) * t208 + t276) + mrSges(4,1) * t159 + t228 * t178 - Ifges(4,6) * t208 + Ifges(4,5) * t209 + Ifges(4,3) * qJDD(3) - t267 + t298 * t227;
t189 = mrSges(4,1) * t227 + mrSges(4,2) * t228;
t215 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t227;
t303 = -mrSges(6,1) - mrSges(4,3);
t133 = m(4) * t159 + (-t187 - t189) * t228 + t303 * t209 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t215 - t219) * qJD(3) + t284 + t309;
t216 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t228;
t136 = m(4) * t160 - qJDD(3) * mrSges(4,2) - qJD(3) * t216 + (-t189 + t297) * t227 + (-mrSges(5,1) + t303) * t208 + t276;
t129 = t306 * t133 + t258 * t136;
t210 = -t229 * t256 + t288;
t281 = Ifges(3,4) * t256 + Ifges(3,2) * t257;
t282 = Ifges(3,1) * t256 + Ifges(3,4) * t257;
t310 = qJD(1) * t257;
t311 = -mrSges(3,1) * t210 + mrSges(3,2) * t211 - pkin(2) * t129 - (t281 * t296 - t282 * t310) * qJD(1) - t312;
t302 = Ifges(4,4) + Ifges(5,6);
t301 = mrSges(3,2) * t256;
t179 = Ifges(6,1) * qJD(3) + Ifges(6,4) * t227 + Ifges(6,5) * t228;
t180 = Ifges(5,1) * qJD(3) - Ifges(5,4) * t228 + Ifges(5,5) * t227;
t299 = t179 + t180;
t235 = g(1) * t259 - t260 * g(2);
t278 = mrSges(3,3) * qJDD(1) + t262 * (-mrSges(3,1) * t257 + t301);
t127 = m(3) * t210 - t278 * t256 + t129;
t285 = -t258 * t133 + t306 * t136;
t128 = m(3) * t211 + t278 * t257 + t285;
t286 = -t127 * t256 + t257 * t128;
t283 = qJDD(2) - t235;
t280 = Ifges(3,5) * t256 + Ifges(3,6) * t257;
t245 = t256 ^ 2;
t207 = (-pkin(1) - t304) * qJDD(1) + (-qJ(2) + (-t245 - t246) * pkin(6)) * t262 + t283;
t265 = pkin(3) * t295 + t228 * t307 + (-t209 + t294) * qJ(4) + t207;
t147 = -pkin(4) * t226 + 0.2e1 * qJD(5) * t227 - t217 * t228 + (pkin(3) + qJ(5)) * t208 + t265;
t142 = m(6) * t147 - t209 * mrSges(6,2) + t208 * mrSges(6,3) - t228 * t218 + t227 * t220;
t275 = mrSges(6,1) * t152 - mrSges(6,3) * t147 - Ifges(6,4) * qJDD(3) - Ifges(6,2) * t208 - Ifges(6,6) * t209 - t228 * t179;
t274 = -mrSges(6,1) * t149 + mrSges(6,2) * t147 - Ifges(6,5) * qJDD(3) - Ifges(6,6) * t208 - Ifges(6,3) * t209 - qJD(3) * t176;
t154 = pkin(3) * t208 + t265;
t137 = m(5) * t154 - t208 * mrSges(5,2) - t209 * mrSges(5,3) - t227 * t219 - t228 * t221 + t142;
t175 = Ifges(4,5) * t228 - Ifges(4,6) * t227 + Ifges(4,3) * qJD(3);
t269 = mrSges(5,1) * t155 - mrSges(5,2) * t154 + pkin(4) * (mrSges(6,1) * t208 + t187 * t227 - t290) + qJ(5) * t142 - t275;
t124 = -t269 - mrSges(4,1) * t207 + mrSges(4,3) * t160 - pkin(3) * t137 + (-Ifges(4,2) - Ifges(5,3)) * t208 + (-t175 - t180) * t228 + t302 * t209 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t173 + t298) * qJD(3);
t270 = -mrSges(5,1) * t157 + mrSges(5,3) * t154 - pkin(4) * t143 + t274;
t125 = -t270 + (-t175 - t299) * t227 + (Ifges(4,1) + Ifges(5,2)) * t209 - t302 * t208 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) + (-t178 + t174) * qJD(3) + mrSges(4,2) * t207 - mrSges(4,3) * t159 - qJ(4) * t137;
t225 = -qJDD(1) * pkin(1) - qJ(2) * t262 + t283;
t231 = t280 * qJD(1);
t268 = m(4) * t207 + t208 * mrSges(4,1) + mrSges(4,2) * t209 + t227 * t215 + t216 * t228 + t137;
t118 = -mrSges(3,1) * t225 + mrSges(3,3) * t211 - pkin(2) * t268 + pkin(6) * t285 + t281 * qJDD(1) + t306 * t124 + t258 * t125 - t231 * t296;
t120 = mrSges(3,2) * t225 - mrSges(3,3) * t210 - pkin(6) * t129 + t282 * qJDD(1) - t258 * t124 + t306 * t125 + t231 * t310;
t266 = -m(3) * t225 + mrSges(3,1) * t291 - t268 + (t245 * t262 + t300) * mrSges(3,3);
t273 = -mrSges(2,2) * t236 + qJ(2) * t286 + t257 * t118 + t256 * t120 + pkin(1) * (-mrSges(3,2) * t292 + t266) + mrSges(2,1) * t235 + Ifges(2,3) * qJDD(1);
t130 = t266 + (mrSges(2,1) - t301) * qJDD(1) + m(2) * t235 - mrSges(2,2) * t262;
t123 = t127 * t257 + t128 * t256;
t121 = m(2) * t236 - mrSges(2,1) * t262 - qJDD(1) * mrSges(2,2) + t286;
t116 = t262 * Ifges(2,5) + mrSges(2,3) * t236 + mrSges(2,1) * g(3) - pkin(1) * t123 + (Ifges(2,6) - t280) * qJDD(1) + t311;
t115 = -mrSges(2,2) * g(3) - mrSges(2,3) * t235 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t262 - qJ(2) * t123 - t118 * t256 + t120 * t257;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t260 * t115 - t259 * t116 - pkin(5) * (t121 * t259 + t130 * t260), t115, t120, t125, -t227 * t177 - t267, -t271; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t259 * t115 + t260 * t116 + pkin(5) * (t121 * t260 - t130 * t259), t116, t118, t124, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t209 + Ifges(5,6) * t208 - qJD(3) * t174 + t299 * t227 + t270, -qJD(3) * t173 - t275; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t273, t273, t280 * qJDD(1) - t311, t312, t269 + t228 * t180 + Ifges(5,5) * qJDD(3) - Ifges(5,6) * t209 + Ifges(5,3) * t208 + (-t173 + t177) * qJD(3), -t227 * t179 - t274;];
m_new = t1;
