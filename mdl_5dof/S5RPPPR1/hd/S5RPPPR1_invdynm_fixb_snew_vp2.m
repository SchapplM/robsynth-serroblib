% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:51
% EndTime: 2019-12-05 17:28:56
% DurationCPUTime: 4.65s
% Computational Cost: add. (53134->235), mult. (123562->328), div. (0->0), fcn. (74901->10), ass. (0->117)
t234 = sin(qJ(1));
t236 = cos(qJ(1));
t210 = t236 * g(2) + t234 * g(3);
t203 = qJDD(1) * pkin(1) + t210;
t209 = t234 * g(2) - t236 * g(3);
t237 = qJD(1) ^ 2;
t204 = -t237 * pkin(1) + t209;
t229 = sin(pkin(7));
t232 = cos(pkin(7));
t184 = t229 * t203 + t232 * t204;
t281 = -t237 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t184;
t183 = t232 * t203 - t229 * t204;
t248 = -t237 * qJ(3) + qJDD(3) - t183;
t228 = sin(pkin(8));
t231 = cos(pkin(8));
t257 = -pkin(3) * t231 - qJ(4) * t228;
t273 = qJD(1) * t228;
t280 = (-pkin(2) + t257) * qJDD(1) + t248 - 0.2e1 * qJD(4) * t273;
t226 = -g(1) + qJDD(2);
t166 = t231 * t226 - t281 * t228;
t272 = t231 * qJD(1);
t167 = t228 * t226 + t281 * t231;
t198 = t257 * qJD(1);
t160 = t198 * t272 + t167;
t223 = t228 ^ 2;
t227 = sin(pkin(9));
t230 = cos(pkin(9));
t275 = t228 * t230;
t253 = -pkin(4) * t231 - pkin(6) * t275;
t274 = t280 * t230;
t152 = t253 * qJDD(1) + (-t160 + (-pkin(4) * t223 * t230 + pkin(6) * t228 * t231) * t237) * t227 + t274;
t155 = t230 * t160 + t280 * t227;
t197 = t253 * qJD(1);
t277 = t223 * t237;
t268 = t227 ^ 2 * t277;
t270 = qJDD(1) * t228;
t153 = -t227 * pkin(6) * t270 - pkin(4) * t268 + t197 * t272 + t155;
t233 = sin(qJ(5));
t235 = cos(qJ(5));
t151 = t233 * t152 + t235 * t153;
t159 = t198 * t273 + qJDD(4) - t166;
t156 = -pkin(6) * t268 + (pkin(4) * qJDD(1) * t227 + qJD(1) * t197 * t230) * t228 + t159;
t250 = (-t227 * t235 - t230 * t233) * t228;
t188 = qJD(1) * t250;
t249 = (-t227 * t233 + t230 * t235) * t228;
t189 = qJD(1) * t249;
t212 = qJD(5) - t272;
t161 = Ifges(6,5) * t189 + Ifges(6,6) * t188 + Ifges(6,3) * t212;
t163 = Ifges(6,1) * t189 + Ifges(6,4) * t188 + Ifges(6,5) * t212;
t174 = -t189 * qJD(5) + qJDD(1) * t250;
t175 = t188 * qJD(5) + qJDD(1) * t249;
t269 = t231 * qJDD(1);
t211 = qJDD(5) - t269;
t139 = -mrSges(6,1) * t156 + mrSges(6,3) * t151 + Ifges(6,4) * t175 + Ifges(6,2) * t174 + Ifges(6,6) * t211 - t189 * t161 + t212 * t163;
t150 = t235 * t152 - t233 * t153;
t162 = Ifges(6,4) * t189 + Ifges(6,2) * t188 + Ifges(6,6) * t212;
t140 = mrSges(6,2) * t156 - mrSges(6,3) * t150 + Ifges(6,1) * t175 + Ifges(6,4) * t174 + Ifges(6,5) * t211 + t188 * t161 - t212 * t162;
t258 = Ifges(5,5) * t230 - Ifges(5,6) * t227;
t185 = (-Ifges(5,3) * t231 + t258 * t228) * qJD(1);
t245 = -Ifges(5,5) * t231 + (Ifges(5,1) * t230 - Ifges(5,4) * t227) * t228;
t187 = t245 * qJD(1);
t244 = -Ifges(5,6) * t231 + (Ifges(5,4) * t230 - Ifges(5,2) * t227) * t228;
t181 = -t212 * mrSges(6,2) + t188 * mrSges(6,3);
t182 = t212 * mrSges(6,1) - t189 * mrSges(6,3);
t246 = m(6) * t156 - t174 * mrSges(6,1) + t175 * mrSges(6,2) - t188 * t181 + t189 * t182;
t170 = -t188 * mrSges(6,1) + t189 * mrSges(6,2);
t146 = m(6) * t150 + t211 * mrSges(6,1) - t175 * mrSges(6,3) - t189 * t170 + t212 * t181;
t147 = m(6) * t151 - t211 * mrSges(6,2) + t174 * mrSges(6,3) + t188 * t170 - t212 * t182;
t264 = -t233 * t146 + t235 * t147;
t125 = -mrSges(5,1) * t159 + mrSges(5,3) * t155 + t233 * t140 + t235 * t139 - pkin(4) * t246 + pkin(6) * t264 + (-t185 * t275 - t231 * t187) * qJD(1) + t244 * qJDD(1);
t138 = t235 * t146 + t233 * t147;
t154 = -t227 * t160 + t274;
t186 = t244 * qJD(1);
t276 = t227 * t228;
t126 = mrSges(5,2) * t159 - mrSges(5,3) * t154 - pkin(6) * t138 - t233 * t139 + t235 * t140 + (-t185 * t276 + t186 * t231) * qJD(1) + t245 * qJDD(1);
t261 = mrSges(5,1) * t227 + mrSges(5,2) * t230;
t190 = t261 * t273;
t251 = mrSges(5,2) * t231 - mrSges(5,3) * t276;
t192 = t251 * qJD(1);
t252 = -mrSges(5,1) * t231 - mrSges(5,3) * t275;
t136 = m(5) * t154 + t252 * qJDD(1) + (-t190 * t275 - t192 * t231) * qJD(1) + t138;
t193 = t252 * qJD(1);
t137 = m(5) * t155 + t251 * qJDD(1) + (-t190 * t276 + t193 * t231) * qJD(1) + t264;
t134 = -t227 * t136 + t230 * t137;
t241 = -m(5) * t159 - t246;
t255 = -t192 * t227 - t193 * t230;
t260 = Ifges(4,1) * t228 + Ifges(4,4) * t231;
t279 = -((Ifges(4,4) * t228 + Ifges(4,2) * t231) * t273 - t260 * t272) * qJD(1) - mrSges(4,1) * t166 + mrSges(4,2) * t167 - pkin(3) * ((t255 * qJD(1) - t261 * qJDD(1)) * t228 + t241) - qJ(4) * t134 - t230 * t125 - t227 * t126;
t278 = mrSges(4,2) * t228;
t199 = (-mrSges(4,1) * t231 + t278) * qJD(1);
t131 = m(4) * t167 + (qJDD(1) * mrSges(4,3) + qJD(1) * t199) * t231 + t134;
t142 = m(4) * t166 + ((-mrSges(4,3) - t261) * qJDD(1) + (-t199 + t255) * qJD(1)) * t228 + t241;
t265 = t231 * t131 - t228 * t142;
t122 = m(3) * t184 - t237 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t265;
t133 = t230 * t136 + t227 * t137;
t178 = -qJDD(1) * pkin(2) + t248;
t242 = -m(4) * t178 + mrSges(4,1) * t269 - t133 + (t231 ^ 2 * t237 + t277) * mrSges(4,3);
t128 = m(3) * t183 - t237 * mrSges(3,2) + (mrSges(3,1) - t278) * qJDD(1) + t242;
t117 = t229 * t122 + t232 * t128;
t124 = t228 * t131 + t231 * t142;
t266 = t232 * t122 - t229 * t128;
t259 = Ifges(4,5) * t228 + Ifges(4,6) * t231;
t256 = t186 * t230 + t187 * t227;
t200 = t259 * qJD(1);
t113 = mrSges(4,2) * t178 - mrSges(4,3) * t166 - qJ(4) * t133 + t260 * qJDD(1) - t227 * t125 + t230 * t126 + t200 * t272;
t243 = -mrSges(6,1) * t150 + mrSges(6,2) * t151 - Ifges(6,5) * t175 - Ifges(6,6) * t174 - Ifges(6,3) * t211 - t189 * t162 + t188 * t163;
t238 = mrSges(5,1) * t154 - mrSges(5,2) * t155 + pkin(4) * t138 - t243;
t119 = -mrSges(4,1) * t178 + mrSges(4,3) * t167 - pkin(3) * t133 + (Ifges(4,2) + Ifges(5,3)) * t269 + ((Ifges(4,4) - t258) * qJDD(1) + (-t200 - t256) * qJD(1)) * t228 - t238;
t247 = -mrSges(3,2) * t184 + qJ(3) * t265 + t228 * t113 + t231 * t119 + pkin(2) * (-mrSges(4,2) * t270 + t242) + mrSges(3,1) * t183 + Ifges(3,3) * qJDD(1);
t240 = mrSges(2,1) * t210 - mrSges(2,2) * t209 + Ifges(2,3) * qJDD(1) + pkin(1) * t117 + t247;
t115 = m(2) * t210 + qJDD(1) * mrSges(2,1) - t237 * mrSges(2,2) + t117;
t114 = m(2) * t209 - t237 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t266;
t111 = -mrSges(3,1) * t226 + mrSges(3,3) * t184 + t237 * Ifges(3,5) - pkin(2) * t124 + (Ifges(3,6) - t259) * qJDD(1) + t279;
t110 = mrSges(3,2) * t226 - mrSges(3,3) * t183 + Ifges(3,5) * qJDD(1) - t237 * Ifges(3,6) - qJ(3) * t124 + t231 * t113 - t228 * t119;
t109 = -mrSges(2,2) * g(1) - mrSges(2,3) * t210 + Ifges(2,5) * qJDD(1) - t237 * Ifges(2,6) - qJ(2) * t117 + t232 * t110 - t229 * t111;
t108 = Ifges(2,6) * qJDD(1) + t237 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t209 + t229 * t110 + t232 * t111 - pkin(1) * (m(3) * t226 + t124) + qJ(2) * t266;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t240, t109, t110, t113, t126, t140; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t234 * t109 - t236 * t108 - pkin(5) * (t236 * t114 - t234 * t115), t108, t111, t119, t125, t139; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t236 * t109 - t234 * t108 + pkin(5) * (-t234 * t114 - t236 * t115), t240, t247, t259 * qJDD(1) - t279, -Ifges(5,3) * t269 + (t256 * qJD(1) + t258 * qJDD(1)) * t228 + t238, -t243;];
m_new = t1;
