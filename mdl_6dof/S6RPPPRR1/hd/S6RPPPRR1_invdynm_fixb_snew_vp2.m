% Calculate vector of cutting torques with Newton-Euler for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:31:58
% EndTime: 2019-05-05 13:32:02
% DurationCPUTime: 2.34s
% Computational Cost: add. (32103->254), mult. (56486->298), div. (0->0), fcn. (27371->8), ass. (0->104)
t248 = sin(qJ(1));
t251 = cos(qJ(1));
t215 = t248 * g(1) - g(2) * t251;
t206 = qJDD(1) * pkin(1) + t215;
t216 = -g(1) * t251 - g(2) * t248;
t253 = qJD(1) ^ 2;
t208 = -pkin(1) * t253 + t216;
t244 = sin(pkin(9));
t245 = cos(pkin(9));
t183 = t245 * t206 - t244 * t208;
t174 = -qJDD(1) * pkin(2) - t253 * qJ(3) + qJDD(3) - t183;
t170 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t174;
t184 = t244 * t206 + t245 * t208;
t287 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t184;
t286 = mrSges(3,1) - mrSges(4,2);
t285 = mrSges(5,2) + mrSges(4,3);
t284 = Ifges(4,4) - Ifges(5,5);
t283 = Ifges(3,6) - Ifges(4,5);
t282 = mrSges(5,3) * t253;
t238 = -g(3) + qJDD(2);
t247 = sin(qJ(5));
t281 = t247 * t238;
t172 = pkin(2) * t253 - t287;
t171 = qJDD(4) + (-pkin(2) - qJ(4)) * t253 + t287;
t166 = -qJDD(1) * pkin(7) + t171;
t250 = cos(qJ(5));
t162 = t247 * t166 + t250 * t238;
t207 = (mrSges(6,1) * t247 + mrSges(6,2) * t250) * qJD(1);
t278 = qJD(1) * qJD(5);
t273 = t250 * t278;
t210 = -qJDD(1) * t247 - t273;
t279 = qJD(1) * t250;
t214 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t279;
t165 = -t253 * pkin(7) - t170;
t274 = t247 * t278;
t211 = qJDD(1) * t250 - t274;
t156 = (-t211 + t274) * pkin(8) + (-t210 + t273) * pkin(5) + t165;
t209 = (pkin(5) * t247 - pkin(8) * t250) * qJD(1);
t252 = qJD(5) ^ 2;
t280 = qJD(1) * t247;
t158 = -pkin(5) * t252 + qJDD(5) * pkin(8) - t209 * t280 + t162;
t246 = sin(qJ(6));
t249 = cos(qJ(6));
t154 = t156 * t249 - t158 * t246;
t204 = qJD(5) * t249 - t246 * t279;
t181 = qJD(6) * t204 + qJDD(5) * t246 + t211 * t249;
t205 = qJD(5) * t246 + t249 * t279;
t185 = -mrSges(7,1) * t204 + mrSges(7,2) * t205;
t217 = qJD(6) + t280;
t186 = -mrSges(7,2) * t217 + mrSges(7,3) * t204;
t203 = qJDD(6) - t210;
t150 = m(7) * t154 + mrSges(7,1) * t203 - mrSges(7,3) * t181 - t185 * t205 + t186 * t217;
t155 = t156 * t246 + t158 * t249;
t180 = -qJD(6) * t205 + qJDD(5) * t249 - t211 * t246;
t187 = mrSges(7,1) * t217 - mrSges(7,3) * t205;
t151 = m(7) * t155 - mrSges(7,2) * t203 + mrSges(7,3) * t180 + t185 * t204 - t187 * t217;
t270 = -t150 * t246 + t249 * t151;
t136 = m(6) * t162 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t210 - qJD(5) * t214 - t207 * t280 + t270;
t161 = t166 * t250 - t281;
t213 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t280;
t157 = -qJDD(5) * pkin(5) - t252 * pkin(8) + t281 + (qJD(1) * t209 - t166) * t250;
t264 = -m(7) * t157 + t180 * mrSges(7,1) - mrSges(7,2) * t181 + t204 * t186 - t187 * t205;
t146 = m(6) * t161 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t211 + qJD(5) * t213 - t207 * t279 + t264;
t127 = t247 * t136 + t250 * t146;
t269 = m(5) * t171 + qJDD(1) * mrSges(5,2) + t127;
t265 = -m(4) * t172 + t253 * mrSges(4,2) + qJDD(1) * mrSges(4,3) + t269;
t121 = m(3) * t184 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) - mrSges(5,3)) * t253 + t265;
t139 = t249 * t150 + t246 * t151;
t266 = -m(6) * t165 + t210 * mrSges(6,1) - t211 * mrSges(6,2) - t213 * t280 - t214 * t279 - t139;
t132 = m(5) * t170 - t253 * mrSges(5,2) - qJDD(1) * mrSges(5,3) + t266;
t260 = -m(4) * t174 + t253 * mrSges(4,3) - t132;
t130 = m(3) * t183 - t253 * mrSges(3,2) + qJDD(1) * t286 + t260;
t113 = t244 * t121 + t245 * t130;
t272 = t245 * t121 - t130 * t244;
t271 = t250 * t136 - t146 * t247;
t125 = m(5) * t238 + t271;
t124 = m(4) * t238 + t125;
t175 = Ifges(7,5) * t205 + Ifges(7,6) * t204 + Ifges(7,3) * t217;
t177 = Ifges(7,1) * t205 + Ifges(7,4) * t204 + Ifges(7,5) * t217;
t143 = -mrSges(7,1) * t157 + mrSges(7,3) * t155 + Ifges(7,4) * t181 + Ifges(7,2) * t180 + Ifges(7,6) * t203 - t175 * t205 + t177 * t217;
t176 = Ifges(7,4) * t205 + Ifges(7,2) * t204 + Ifges(7,6) * t217;
t144 = mrSges(7,2) * t157 - mrSges(7,3) * t154 + Ifges(7,1) * t181 + Ifges(7,4) * t180 + Ifges(7,5) * t203 + t175 * t204 - t176 * t217;
t194 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t250 - Ifges(6,6) * t247) * qJD(1);
t195 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t250 - Ifges(6,2) * t247) * qJD(1);
t116 = mrSges(6,2) * t165 - mrSges(6,3) * t161 + Ifges(6,1) * t211 + Ifges(6,4) * t210 + Ifges(6,5) * qJDD(5) - pkin(8) * t139 - qJD(5) * t195 - t143 * t246 + t144 * t249 - t194 * t280;
t196 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t250 - Ifges(6,4) * t247) * qJD(1);
t259 = mrSges(7,1) * t154 - mrSges(7,2) * t155 + Ifges(7,5) * t181 + Ifges(7,6) * t180 + Ifges(7,3) * t203 + t176 * t205 - t177 * t204;
t118 = -mrSges(6,1) * t165 + mrSges(6,3) * t162 + Ifges(6,4) * t211 + Ifges(6,2) * t210 + Ifges(6,6) * qJDD(5) - pkin(5) * t139 + qJD(5) * t196 - t194 * t279 - t259;
t267 = mrSges(5,1) * t170 + t253 * Ifges(5,4) + Ifges(5,5) * qJDD(1) + pkin(4) * t266 + pkin(7) * t271 + t247 * t116 + t250 * t118;
t263 = mrSges(5,2) * t171 - mrSges(5,3) * t170 + Ifges(5,1) * qJDD(1) - pkin(7) * t127 + t250 * t116 - t118 * t247;
t262 = mrSges(6,1) * t161 - mrSges(6,2) * t162 + Ifges(6,5) * t211 + Ifges(6,6) * t210 + Ifges(6,3) * qJDD(5) + pkin(5) * t264 + pkin(8) * t270 + t249 * t143 + t246 * t144 + t195 * t279 + t196 * t280;
t261 = -mrSges(4,1) * t174 - pkin(3) * t132 - t267;
t258 = mrSges(4,2) * t174 - mrSges(4,3) * t172 + Ifges(4,1) * qJDD(1) - qJ(4) * t132 + t263;
t257 = mrSges(5,1) * t171 - mrSges(5,3) * t238 - Ifges(5,4) * qJDD(1) + pkin(4) * t127 + t262;
t256 = -mrSges(3,2) * t184 + qJ(3) * (t265 - t282) + pkin(2) * (-qJDD(1) * mrSges(4,2) + t260) + mrSges(3,1) * t183 + Ifges(3,3) * qJDD(1) + t258;
t255 = mrSges(4,1) * t172 + pkin(3) * (-t269 + t282) + qJ(4) * t125 - t257;
t254 = mrSges(2,1) * t215 - mrSges(2,2) * t216 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t256;
t111 = m(2) * t216 - mrSges(2,1) * t253 - qJDD(1) * mrSges(2,2) + t272;
t110 = m(2) * t215 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t253 + t113;
t109 = -t255 + (Ifges(3,5) - t284) * t253 - t286 * t238 + t283 * qJDD(1) + mrSges(3,3) * t184 - pkin(2) * t124;
t108 = -mrSges(3,3) * t183 - qJ(3) * t124 - t283 * t253 + (-Ifges(4,4) + Ifges(3,5)) * qJDD(1) + (mrSges(3,2) - t285) * t238 - t261;
t107 = -mrSges(2,2) * g(3) - mrSges(2,3) * t215 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t253 - qJ(2) * t113 + t108 * t245 - t109 * t244;
t106 = Ifges(2,6) * qJDD(1) + t253 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t216 + t244 * t108 + t245 * t109 - pkin(1) * (m(3) * t238 + t124) + qJ(2) * t272;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t251 * t107 - t248 * t106 - pkin(6) * (t110 * t251 + t111 * t248), t107, t108, t258, t263, t116, t144; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t248 * t107 + t251 * t106 + pkin(6) * (-t110 * t248 + t111 * t251), t106, t109, Ifges(4,4) * qJDD(1) - t253 * Ifges(4,5) + t238 * t285 + t261, -t253 * Ifges(5,5) - t257, t118, t143; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t254, t254, t256, -mrSges(4,2) * t238 + Ifges(4,5) * qJDD(1) + t253 * t284 + t255, -mrSges(5,2) * t238 + t267, t262, t259;];
m_new  = t1;
