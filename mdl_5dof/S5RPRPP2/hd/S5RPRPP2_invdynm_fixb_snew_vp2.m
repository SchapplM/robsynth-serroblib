% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPP2
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:43
% EndTime: 2019-12-31 18:10:45
% DurationCPUTime: 1.69s
% Computational Cost: add. (13092->268), mult. (25249->318), div. (0->0), fcn. (10839->6), ass. (0->100)
t245 = sin(qJ(1));
t247 = cos(qJ(1));
t224 = t245 * g(1) - t247 * g(2);
t204 = qJDD(1) * pkin(1) + t224;
t225 = -t247 * g(1) - t245 * g(2);
t249 = qJD(1) ^ 2;
t209 = -t249 * pkin(1) + t225;
t242 = sin(pkin(7));
t243 = cos(pkin(7));
t160 = t242 * t204 + t243 * t209;
t157 = -t249 * pkin(2) + qJDD(1) * pkin(6) + t160;
t244 = sin(qJ(3));
t154 = t244 * t157;
t241 = -g(3) + qJDD(2);
t246 = cos(qJ(3));
t278 = t246 * t241;
t152 = -t154 + t278;
t153 = t246 * t157 + t244 * t241;
t175 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t244 - Ifges(5,3) * t246) * qJD(1);
t179 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t244 + Ifges(4,2) * t246) * qJD(1);
t206 = (-mrSges(5,1) * t246 - mrSges(5,3) * t244) * qJD(1);
t271 = qJD(1) * qJD(3);
t210 = t244 * qJDD(1) + t246 * t271;
t211 = t246 * qJDD(1) - t244 * t271;
t205 = (-pkin(3) * t246 - qJ(4) * t244) * qJD(1);
t248 = qJD(3) ^ 2;
t273 = qJD(1) * t244;
t262 = -t248 * qJ(4) + t205 * t273 + qJDD(4) + t154;
t270 = -0.2e1 * qJD(1) * qJD(5);
t283 = pkin(4) * t249;
t284 = pkin(3) + pkin(4);
t145 = t244 * t270 - t210 * qJ(5) - t284 * qJDD(3) + (qJ(5) * t271 - t244 * t283 - t241) * t246 + t262;
t207 = (mrSges(6,1) * t246 + mrSges(6,2) * t244) * qJD(1);
t272 = qJD(1) * t246;
t221 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t272;
t136 = m(6) * t145 - qJDD(3) * mrSges(6,1) - t210 * mrSges(6,3) - qJD(3) * t221 - t207 * t273;
t285 = 2 * qJD(4);
t149 = -t248 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t285 + t205 * t272 + t153;
t151 = -qJDD(3) * pkin(3) + t262 - t278;
t217 = -qJD(3) * pkin(4) - qJ(5) * t273;
t240 = t246 ^ 2;
t144 = -t211 * qJ(5) + qJD(3) * t217 - t240 * t283 + t246 * t270 + t149;
t177 = -Ifges(6,6) * qJD(3) + (Ifges(6,4) * t244 - Ifges(6,2) * t246) * qJD(1);
t180 = -Ifges(6,5) * qJD(3) + (Ifges(6,1) * t244 - Ifges(6,4) * t246) * qJD(1);
t260 = mrSges(6,1) * t145 - mrSges(6,2) * t144 + Ifges(6,5) * t210 - Ifges(6,6) * t211 - Ifges(6,3) * qJDD(3) + t177 * t273 + t180 * t272;
t253 = -mrSges(5,1) * t151 + mrSges(5,3) * t149 + Ifges(5,4) * t210 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t211 - pkin(4) * t136 - t260;
t223 = mrSges(5,2) * t272 + qJD(3) * mrSges(5,3);
t256 = -m(5) * t151 + qJDD(3) * mrSges(5,1) + qJD(3) * t223 - t136;
t220 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t273;
t218 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t273;
t265 = m(6) * t144 + qJDD(3) * mrSges(6,2) - t211 * mrSges(6,3) + qJD(3) * t218;
t258 = m(5) * t149 + qJDD(3) * mrSges(5,3) + qJD(3) * t220 + t206 * t272 + t265;
t268 = t207 * t272;
t181 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t244 - Ifges(5,5) * t246) * qJD(1);
t275 = t181 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t244 + Ifges(4,4) * t246) * qJD(1);
t288 = -(t275 * t246 + (t175 - t179) * t244) * qJD(1) + mrSges(4,1) * t152 - mrSges(4,2) * t153 + Ifges(4,5) * t210 + Ifges(4,6) * t211 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t210 * mrSges(5,2) - t206 * t273 + t256) + qJ(4) * (t211 * mrSges(5,2) + t258 - t268) + t253;
t174 = -Ifges(6,3) * qJD(3) + (Ifges(6,5) * t244 - Ifges(6,6) * t246) * qJD(1);
t287 = -Ifges(6,5) * qJDD(3) - t174 * t272;
t282 = t249 * pkin(6);
t281 = mrSges(4,3) + mrSges(5,2);
t280 = Ifges(5,6) - Ifges(6,6);
t279 = qJ(4) * t246;
t208 = (-mrSges(4,1) * t246 + mrSges(4,2) * t244) * qJD(1);
t219 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t273;
t129 = m(4) * t153 - qJDD(3) * mrSges(4,2) - qJD(3) * t219 + t281 * t211 + (-t207 + t208) * t272 + t258;
t222 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t272;
t130 = m(4) * t152 + qJDD(3) * mrSges(4,1) + qJD(3) * t222 - t281 * t210 + (-t206 - t208) * t273 + t256;
t266 = t246 * t129 - t244 * t130;
t120 = m(3) * t160 - t249 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t266;
t159 = t243 * t204 - t242 * t209;
t264 = qJDD(1) * pkin(2) + t159;
t259 = -t210 * qJ(4) - t264;
t139 = qJDD(5) + (-qJ(5) * t240 + pkin(6)) * t249 + t284 * t211 + (qJD(3) * t279 + (-pkin(3) * qJD(3) + t217 + t285) * t244) * qJD(1) - t259;
t134 = -m(6) * t139 - t211 * mrSges(6,1) - t210 * mrSges(6,2) - t218 * t273 - t221 * t272;
t146 = -t211 * pkin(3) - t282 + (-0.2e1 * qJD(4) * t244 + (pkin(3) * t244 - t279) * qJD(3)) * qJD(1) + t259;
t131 = m(5) * t146 - t211 * mrSges(5,1) - t210 * mrSges(5,3) - t220 * t273 - t223 * t272 + t134;
t156 = -t264 - t282;
t251 = -m(4) * t156 + t211 * mrSges(4,1) - t210 * mrSges(4,2) - t219 * t273 + t222 * t272 - t131;
t124 = m(3) * t159 + qJDD(1) * mrSges(3,1) - t249 * mrSges(3,2) + t251;
t115 = t242 * t120 + t243 * t124;
t122 = t244 * t129 + t246 * t130;
t178 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t244 - Ifges(5,6) * t246) * qJD(1);
t277 = -t174 + t178;
t267 = t243 * t120 - t242 * t124;
t263 = mrSges(6,1) * t139 - mrSges(6,3) * t144 - Ifges(6,4) * t210 + Ifges(6,2) * t211;
t261 = mrSges(6,2) * t139 - mrSges(6,3) * t145 + Ifges(6,1) * t210 - Ifges(6,4) * t211 + qJD(3) * t177;
t176 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t244 + Ifges(4,6) * t246) * qJD(1);
t254 = mrSges(5,1) * t146 - mrSges(5,2) * t149 + pkin(4) * t134 + qJ(5) * (t265 - t268) - t263;
t111 = -t254 + (Ifges(4,2) + Ifges(5,3)) * t211 + (Ifges(4,4) - Ifges(5,5)) * t210 + (Ifges(4,6) - t280) * qJDD(3) + (t180 + t275) * qJD(3) + mrSges(4,3) * t153 - mrSges(4,1) * t156 - pkin(3) * t131 + (-t176 - t277) * t273;
t252 = mrSges(5,2) * t151 - mrSges(5,3) * t146 + Ifges(5,1) * t210 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t211 - qJ(5) * t136 + qJD(3) * t175 + t178 * t272 + t261;
t117 = t252 + (Ifges(4,5) - Ifges(6,5)) * qJDD(3) + Ifges(4,1) * t210 + Ifges(4,4) * t211 - mrSges(4,3) * t152 + mrSges(4,2) * t156 - qJD(3) * t179 - qJ(4) * t131 + (-t174 + t176) * t272;
t257 = mrSges(3,1) * t159 - mrSges(3,2) * t160 + Ifges(3,3) * qJDD(1) + pkin(2) * t251 + pkin(6) * t266 + t246 * t111 + t244 * t117;
t255 = mrSges(2,1) * t224 - mrSges(2,2) * t225 + Ifges(2,3) * qJDD(1) + pkin(1) * t115 + t257;
t113 = m(2) * t225 - t249 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t267;
t112 = m(2) * t224 + qJDD(1) * mrSges(2,1) - t249 * mrSges(2,2) + t115;
t109 = -mrSges(3,1) * t241 + mrSges(3,3) * t160 + t249 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t122 - t288;
t108 = mrSges(3,2) * t241 - mrSges(3,3) * t159 + Ifges(3,5) * qJDD(1) - t249 * Ifges(3,6) - pkin(6) * t122 - t244 * t111 + t246 * t117;
t107 = -mrSges(2,2) * g(3) - mrSges(2,3) * t224 + Ifges(2,5) * qJDD(1) - t249 * Ifges(2,6) - qJ(2) * t115 + t243 * t108 - t242 * t109;
t106 = Ifges(2,6) * qJDD(1) + t249 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t225 + t242 * t108 + t243 * t109 - pkin(1) * (m(3) * t241 + t122) + qJ(2) * t267;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t247 * t107 - t245 * t106 - pkin(5) * (t247 * t112 + t245 * t113), t107, t108, t117, t252 + t287, t261 + t287; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t245 * t107 + t247 * t106 + pkin(5) * (-t245 * t112 + t247 * t113), t106, t109, t111, t253 + (-t244 * t175 - t246 * t181) * qJD(1), -Ifges(6,6) * qJDD(3) - qJD(3) * t180 - t174 * t273 - t263; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t255, t255, t257, t288, Ifges(5,5) * t210 - Ifges(5,3) * t211 + t280 * qJDD(3) + (-t180 - t181) * qJD(3) + t277 * t273 + t254, t260;];
m_new = t1;
