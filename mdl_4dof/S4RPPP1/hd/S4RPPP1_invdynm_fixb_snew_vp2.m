% Calculate vector of cutting torques with Newton-Euler for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:07:14
% EndTime: 2019-05-04 19:07:16
% DurationCPUTime: 1.65s
% Computational Cost: add. (8397->265), mult. (24442->371), div. (0->0), fcn. (14166->6), ass. (0->116)
t209 = sin(qJ(1));
t210 = cos(qJ(1));
t191 = -g(1) * t210 - g(2) * t209;
t211 = qJD(1) ^ 2;
t206 = sin(pkin(4));
t238 = qJDD(1) * t206;
t171 = -pkin(1) * t211 + qJ(2) * t238 + t191;
t205 = sin(pkin(6));
t240 = qJD(1) * t206;
t233 = qJD(2) * t240;
t207 = cos(pkin(6));
t251 = t206 * t207;
t269 = -g(3) * t251 + (-t171 - 0.2e1 * t233) * t205;
t268 = -2 * qJD(3);
t267 = -mrSges(4,1) - mrSges(5,1);
t266 = -mrSges(5,2) - mrSges(4,3);
t265 = -Ifges(4,4) + Ifges(3,5);
t264 = -Ifges(4,5) + Ifges(3,6);
t263 = -Ifges(5,2) - Ifges(4,3);
t262 = -Ifges(4,6) + Ifges(5,6);
t261 = Ifges(4,4) * t205;
t208 = cos(pkin(4));
t260 = Ifges(4,4) * t208;
t259 = Ifges(4,5) * t208;
t258 = Ifges(4,2) * t205;
t257 = Ifges(4,6) * t205;
t190 = t209 * g(1) - g(2) * t210;
t250 = t206 * t211;
t170 = qJDD(1) * pkin(1) + qJ(2) * t250 + t190;
t252 = t205 * t208;
t253 = t205 * t206;
t225 = g(3) * t253 - t170 * t252 - t207 * t171;
t141 = 0.2e1 * t207 * t233 - t225;
t168 = (-mrSges(3,1) * t207 + mrSges(3,2) * t205) * t240;
t172 = (mrSges(3,1) * t208 - mrSges(3,3) * t253) * qJD(1);
t166 = (-pkin(2) * t207 - qJ(3) * t205) * t240;
t174 = (pkin(3) * t253 - qJ(4) * t208) * qJD(1);
t201 = t206 ^ 2;
t202 = t207 ^ 2;
t203 = t208 ^ 2;
t131 = qJDD(4) + (-qJ(4) * t201 * t202 - pkin(2) * t203) * t211 + (pkin(3) * t251 + qJ(3) * t208) * qJDD(1) + (t166 * t251 + ((2 * qJD(3)) + t174) * t208) * qJD(1) + t141;
t226 = -mrSges(5,2) * t205 - mrSges(5,3) * t207;
t169 = t226 * t240;
t175 = (mrSges(5,1) * t253 - mrSges(5,3) * t208) * qJD(1);
t237 = qJDD(1) * t207;
t231 = t206 * t237;
t234 = t207 * t240;
t236 = qJDD(1) * t208;
t239 = qJD(1) * t208;
t125 = -m(5) * t131 - mrSges(5,1) * t231 - mrSges(5,2) * t236 - t169 * t234 - t175 * t239;
t254 = t203 * t211;
t135 = pkin(2) * t254 - qJ(3) * t236 + (t208 * t268 + (-0.2e1 * qJD(2) - t166) * t251) * qJD(1) + t225;
t167 = (mrSges(4,2) * t207 - mrSges(4,3) * t205) * t240;
t178 = (mrSges(4,1) * t253 + mrSges(4,2) * t208) * qJD(1);
t215 = -m(4) * t135 + mrSges(4,1) * t231 + mrSges(4,3) * t236 + t167 * t234 + t178 * t239 - t125;
t222 = -mrSges(3,2) * t208 + mrSges(3,3) * t251;
t120 = m(3) * t141 + t222 * qJDD(1) + (t168 * t251 - t172 * t208) * qJD(1) + t215;
t248 = t207 * t208;
t140 = t170 * t248 + t269;
t173 = t222 * qJD(1);
t176 = (-mrSges(4,1) * t251 - mrSges(4,3) * t208) * qJD(1);
t235 = t205 * t240;
t217 = -qJ(3) * t254 + t166 * t235 + qJDD(3) - t269;
t137 = (-pkin(2) * qJDD(1) - t170 * t207) * t208 + t217;
t219 = (-pkin(2) - qJ(4)) * qJDD(1) - 0.2e1 * qJD(1) * qJD(4);
t255 = t201 * t211;
t130 = (-qJ(4) * t207 * t255 + pkin(3) * t238) * t205 + ((-pkin(3) * t250 - t170) * t207 + t219) * t208 + t217;
t177 = (mrSges(5,1) * t251 + mrSges(5,2) * t208) * qJD(1);
t227 = m(5) * t130 - mrSges(5,3) * t236 - t177 * t239;
t220 = -m(4) * t137 - t227;
t242 = -t167 - t169;
t121 = m(3) * t140 + ((mrSges(3,1) - mrSges(4,2)) * qJDD(1) + (t173 - t176) * qJD(1)) * t208 + ((-mrSges(3,3) + t267) * qJDD(1) + (-t168 + t242) * qJD(1)) * t253 + t220;
t116 = t207 * t120 - t121 * t205;
t256 = qJ(2) * t116;
t154 = (Ifges(5,5) * t208 + (-Ifges(5,6) * t207 + Ifges(5,3) * t205) * t206) * qJD(1);
t249 = t207 * t154;
t247 = t208 * t211;
t159 = (Ifges(4,1) * t208 + (-Ifges(4,5) * t207 - t261) * t206) * qJD(1);
t246 = (Ifges(3,3) * t208 + (Ifges(3,5) * t205 + Ifges(3,6) * t207) * t206) * qJD(1) + t159;
t157 = (t260 + (-Ifges(4,6) * t207 - t258) * t206) * qJD(1);
t245 = -t154 + t157;
t155 = (t259 + (-Ifges(4,3) * t207 - t257) * t206) * qJD(1);
t156 = (Ifges(5,4) * t208 + (-Ifges(5,2) * t207 + Ifges(5,6) * t205) * t206) * qJD(1);
t244 = -t155 - t156;
t241 = -t175 - t178;
t232 = t205 * t238;
t230 = -t208 * g(3) + qJDD(2);
t223 = pkin(2) * t247 * t253 + t235 * t268 + t230;
t134 = -pkin(3) * t202 * t255 + (-t170 + (-qJ(3) * qJDD(1) - qJD(1) * t174) * t205 + (-qJ(3) * t247 + t219) * t207) * t206 + t223;
t133 = m(5) * t134;
t139 = (-pkin(2) * t237 - t170 + (-qJDD(1) * t205 - t207 * t247) * qJ(3)) * t206 + t223;
t229 = m(4) * t139 + mrSges(4,2) * t231 + t176 * t234 + t133;
t146 = -t206 * t170 + t230;
t122 = m(3) * t146 + (((-mrSges(3,1) - mrSges(5,3)) * qJDD(1) + (-t173 - t177) * qJD(1)) * t207 + ((mrSges(3,2) + t266) * qJDD(1) + (t172 + t241) * qJD(1)) * t205) * t206 + t229;
t111 = t120 * t252 + t121 * t248 - t122 * t206;
t224 = mrSges(5,2) * t131 - mrSges(5,3) * t130 + Ifges(5,1) * t236 + Ifges(5,5) * t232;
t158 = (Ifges(5,1) * t208 + (-Ifges(5,4) * t207 + Ifges(5,5) * t205) * t206) * qJD(1);
t221 = mrSges(5,1) * t131 - mrSges(5,3) * t134 - Ifges(5,4) * t236 - Ifges(5,6) * t232 - t158 * t235;
t218 = -mrSges(5,1) * t130 + mrSges(5,2) * t134 - Ifges(5,5) * t236 - Ifges(5,3) * t232 - t156 * t239 - t158 * t234;
t152 = (Ifges(3,6) * t208 + (Ifges(3,4) * t205 + Ifges(3,2) * t207) * t206) * qJD(1);
t153 = (Ifges(3,5) * t208 + (Ifges(3,1) * t205 + Ifges(3,4) * t207) * t206) * qJD(1);
t124 = (mrSges(5,1) * qJDD(1) + qJD(1) * t169) * t253 + t227;
t214 = mrSges(4,2) * t137 - mrSges(4,3) * t135 + Ifges(4,1) * t236 - qJ(4) * t124 + t157 * t234 + t224;
t105 = (((-t153 - t154) * qJD(1) + (-Ifges(5,4) + t264) * qJDD(1)) * t207 + ((pkin(2) * t267 + t265) * qJDD(1) + (pkin(2) * t242 + t152 + t244) * qJD(1)) * t205) * t206 + t214 + pkin(2) * t220 + (Ifges(3,3) * qJDD(1) + pkin(2) * (-qJDD(1) * mrSges(4,2) - qJD(1) * t176)) * t208 + qJ(3) * t215 + mrSges(3,1) * t140 - mrSges(3,2) * t141;
t123 = ((-mrSges(5,3) * qJDD(1) - qJD(1) * t177) * t207 + (t241 * qJD(1) + t266 * qJDD(1)) * t205) * t206 + t229;
t212 = mrSges(4,1) * t135 - mrSges(4,2) * t139 + pkin(3) * t125 + qJ(4) * (t133 + (t226 * qJDD(1) + (-t175 * t205 - t177 * t207) * qJD(1)) * t206) - t221;
t107 = -mrSges(3,1) * t146 + mrSges(3,3) * t141 - pkin(2) * t123 + (-t246 * t253 + (t153 - t245) * t208) * qJD(1) + (t264 * t208 + ((Ifges(3,4) + Ifges(4,6)) * t205 + (Ifges(3,2) - t263) * t207) * t206) * qJDD(1) - t212;
t213 = -mrSges(4,1) * t137 + mrSges(4,3) * t139 - pkin(3) * t124 + t218;
t113 = mrSges(3,2) * t146 - mrSges(3,3) * t140 - qJ(3) * t123 + ((-t152 + t155) * t208 + t246 * t251) * qJD(1) + (t265 * t208 + ((Ifges(3,1) + Ifges(4,2)) * t205 + (Ifges(3,4) - t262) * t207) * t206) * qJDD(1) - t213;
t216 = mrSges(2,1) * t190 - mrSges(2,2) * t191 + Ifges(2,3) * qJDD(1) + pkin(1) * t111 + t208 * t105 + t107 * t251 + t113 * t253 + t206 * t256;
t114 = m(2) * t191 - mrSges(2,1) * t211 - qJDD(1) * mrSges(2,2) + t116;
t110 = t208 * t122 + (t120 * t205 + t121 * t207) * t206;
t108 = m(2) * t190 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t211 + t111;
t103 = -mrSges(2,2) * g(3) - mrSges(2,3) * t190 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t211 - t107 * t205 + t113 * t207 + (-t110 * t206 - t111 * t208) * qJ(2);
t102 = mrSges(2,1) * g(3) + mrSges(2,3) * t191 + Ifges(2,5) * t211 + Ifges(2,6) * qJDD(1) - pkin(1) * t110 - t105 * t206 + (t107 * t207 + t113 * t205 + t256) * t208;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t210 * t103 - t209 * t102 - pkin(5) * (t108 * t210 + t114 * t209), t103, t113, ((-t261 + (-Ifges(5,4) - Ifges(4,5)) * t207) * qJDD(1) + (t244 * t205 - t249) * qJD(1)) * t206 + t214, (-Ifges(5,4) * t237 + (-t205 * t156 - t249) * qJD(1)) * t206 + t224; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t209 * t103 + t210 * t102 + pkin(5) * (-t108 * t209 + t114 * t210), t102, t107, (-t208 * t155 - t159 * t251) * qJD(1) + (t260 + (t262 * t207 - t258) * t206) * qJDD(1) + t213, -Ifges(5,2) * t231 - t154 * t239 - t221; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t216, t216, t105, (t159 * t253 + t245 * t208) * qJD(1) + (t259 + (t263 * t207 - t257) * t206) * qJDD(1) + t212, -Ifges(5,6) * t231 - t218;];
m_new  = t1;
