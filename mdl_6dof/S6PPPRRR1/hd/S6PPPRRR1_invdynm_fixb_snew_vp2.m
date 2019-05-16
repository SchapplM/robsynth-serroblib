% Calculate vector of cutting torques with Newton-Euler for
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-05-04 19:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPPRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:32:07
% EndTime: 2019-05-04 19:32:52
% DurationCPUTime: 45.20s
% Computational Cost: add. (962577->244), mult. (1629952->332), div. (0->0), fcn. (1447275->18), ass. (0->127)
t223 = sin(pkin(12));
t229 = cos(pkin(12));
t216 = -t229 * g(1) - t223 * g(2);
t222 = sin(pkin(13));
t228 = cos(pkin(13));
t215 = t223 * g(1) - t229 * g(2);
t220 = -g(3) + qJDD(1);
t226 = sin(pkin(6));
t232 = cos(pkin(6));
t247 = t215 * t232 + t220 * t226;
t190 = t228 * t216 + t247 * t222;
t221 = sin(pkin(14));
t227 = cos(pkin(14));
t189 = -t222 * t216 + t247 * t228;
t201 = -t226 * t215 + t232 * t220 + qJDD(2);
t225 = sin(pkin(7));
t231 = cos(pkin(7));
t249 = t189 * t231 + t201 * t225;
t185 = -t221 * t190 + t249 * t227;
t186 = t227 * t190 + t249 * t221;
t188 = -t225 * t189 + t231 * t201 + qJDD(3);
t238 = cos(qJ(4));
t230 = cos(pkin(8));
t235 = sin(qJ(4));
t258 = t230 * t235;
t224 = sin(pkin(8));
t259 = t224 * t235;
t179 = t185 * t258 + t238 * t186 + t188 * t259;
t240 = qJD(4) ^ 2;
t177 = -t240 * pkin(4) + qJDD(4) * pkin(10) + t179;
t181 = -t224 * t185 + t230 * t188;
t234 = sin(qJ(5));
t237 = cos(qJ(5));
t173 = t237 * t177 + t234 * t181;
t211 = (-pkin(5) * t237 - pkin(11) * t234) * qJD(4);
t239 = qJD(5) ^ 2;
t255 = t237 * qJD(4);
t171 = -t239 * pkin(5) + qJDD(5) * pkin(11) + t211 * t255 + t173;
t178 = -t235 * t186 + (t185 * t230 + t188 * t224) * t238;
t176 = -qJDD(4) * pkin(4) - t240 * pkin(10) - t178;
t254 = qJD(4) * qJD(5);
t252 = t237 * t254;
t212 = t234 * qJDD(4) + t252;
t253 = t234 * t254;
t213 = t237 * qJDD(4) - t253;
t174 = (-t212 - t252) * pkin(11) + (-t213 + t253) * pkin(5) + t176;
t233 = sin(qJ(6));
t236 = cos(qJ(6));
t167 = -t233 * t171 + t236 * t174;
t256 = qJD(4) * t234;
t208 = t236 * qJD(5) - t233 * t256;
t197 = t208 * qJD(6) + t233 * qJDD(5) + t236 * t212;
t209 = t233 * qJD(5) + t236 * t256;
t198 = -t208 * mrSges(7,1) + t209 * mrSges(7,2);
t219 = qJD(6) - t255;
t199 = -t219 * mrSges(7,2) + t208 * mrSges(7,3);
t207 = qJDD(6) - t213;
t165 = m(7) * t167 + t207 * mrSges(7,1) - t197 * mrSges(7,3) - t209 * t198 + t219 * t199;
t168 = t236 * t171 + t233 * t174;
t196 = -t209 * qJD(6) + t236 * qJDD(5) - t233 * t212;
t200 = t219 * mrSges(7,1) - t209 * mrSges(7,3);
t166 = m(7) * t168 - t207 * mrSges(7,2) + t196 * mrSges(7,3) + t208 * t198 - t219 * t200;
t159 = -t233 * t165 + t236 * t166;
t210 = (-mrSges(6,1) * t237 + mrSges(6,2) * t234) * qJD(4);
t217 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t256;
t157 = m(6) * t173 - qJDD(5) * mrSges(6,2) + t213 * mrSges(6,3) - qJD(5) * t217 + t210 * t255 + t159;
t257 = t237 * t181;
t170 = -qJDD(5) * pkin(5) - t239 * pkin(11) - t257 + (qJD(4) * t211 + t177) * t234;
t169 = -m(7) * t170 + t196 * mrSges(7,1) - t197 * mrSges(7,2) + t208 * t199 - t209 * t200;
t172 = -t234 * t177 + t257;
t218 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t255;
t163 = m(6) * t172 + qJDD(5) * mrSges(6,1) - t212 * mrSges(6,3) + qJD(5) * t218 - t210 * t256 + t169;
t251 = t237 * t157 - t234 * t163;
t148 = m(5) * t179 - t240 * mrSges(5,1) - qJDD(4) * mrSges(5,2) + t251;
t151 = t234 * t157 + t237 * t163;
t150 = m(5) * t181 + t151;
t158 = t236 * t165 + t233 * t166;
t243 = -m(6) * t176 + t213 * mrSges(6,1) - t212 * mrSges(6,2) - t217 * t256 + t218 * t255 - t158;
t154 = m(5) * t178 + qJDD(4) * mrSges(5,1) - t240 * mrSges(5,2) + t243;
t260 = t154 * t238;
t137 = t148 * t258 - t224 * t150 + t230 * t260;
t133 = m(4) * t185 + t137;
t142 = t238 * t148 - t235 * t154;
t141 = m(4) * t186 + t142;
t271 = t133 * t227 + t141 * t221;
t136 = t148 * t259 + t230 * t150 + t224 * t260;
t135 = m(4) * t188 + t136;
t123 = -t225 * t135 + t271 * t231;
t120 = m(3) * t189 + t123;
t127 = -t221 * t133 + t227 * t141;
t126 = m(3) * t190 + t127;
t270 = t120 * t228 + t126 * t222;
t191 = Ifges(7,5) * t209 + Ifges(7,6) * t208 + Ifges(7,3) * t219;
t193 = Ifges(7,1) * t209 + Ifges(7,4) * t208 + Ifges(7,5) * t219;
t160 = -mrSges(7,1) * t170 + mrSges(7,3) * t168 + Ifges(7,4) * t197 + Ifges(7,2) * t196 + Ifges(7,6) * t207 - t209 * t191 + t219 * t193;
t192 = Ifges(7,4) * t209 + Ifges(7,2) * t208 + Ifges(7,6) * t219;
t161 = mrSges(7,2) * t170 - mrSges(7,3) * t167 + Ifges(7,1) * t197 + Ifges(7,4) * t196 + Ifges(7,5) * t207 + t208 * t191 - t219 * t192;
t202 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t234 + Ifges(6,6) * t237) * qJD(4);
t203 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t234 + Ifges(6,2) * t237) * qJD(4);
t143 = mrSges(6,2) * t176 - mrSges(6,3) * t172 + Ifges(6,1) * t212 + Ifges(6,4) * t213 + Ifges(6,5) * qJDD(5) - pkin(11) * t158 - qJD(5) * t203 - t233 * t160 + t236 * t161 + t202 * t255;
t204 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t234 + Ifges(6,4) * t237) * qJD(4);
t242 = mrSges(7,1) * t167 - mrSges(7,2) * t168 + Ifges(7,5) * t197 + Ifges(7,6) * t196 + Ifges(7,3) * t207 + t209 * t192 - t208 * t193;
t144 = -mrSges(6,1) * t176 + mrSges(6,3) * t173 + Ifges(6,4) * t212 + Ifges(6,2) * t213 + Ifges(6,6) * qJDD(5) - pkin(5) * t158 + qJD(5) * t204 - t202 * t256 - t242;
t128 = mrSges(5,1) * t178 - mrSges(5,2) * t179 + Ifges(5,3) * qJDD(4) + pkin(4) * t243 + pkin(10) * t251 + t234 * t143 + t237 * t144;
t129 = mrSges(5,2) * t181 - mrSges(5,3) * t178 + Ifges(5,5) * qJDD(4) - t240 * Ifges(5,6) - pkin(10) * t151 + t237 * t143 - t234 * t144;
t268 = mrSges(6,1) * t172 - mrSges(6,2) * t173 + Ifges(6,5) * t212 + Ifges(6,6) * t213 + Ifges(6,3) * qJDD(5) + pkin(5) * t169 + pkin(11) * t159 + t236 * t160 + t233 * t161 + (t234 * t203 - t237 * t204) * qJD(4);
t130 = -mrSges(5,1) * t181 + mrSges(5,3) * t179 + t240 * Ifges(5,5) + Ifges(5,6) * qJDD(4) - pkin(4) * t151 - t268;
t246 = pkin(9) * t142 + t129 * t235 + t130 * t238;
t113 = mrSges(4,1) * t185 - mrSges(4,2) * t186 + pkin(3) * t137 + t230 * t128 + t246 * t224;
t122 = t231 * t135 + t271 * t225;
t114 = -mrSges(4,1) * t188 + mrSges(4,3) * t186 - pkin(3) * t136 - t224 * t128 + t246 * t230;
t115 = mrSges(4,2) * t188 - mrSges(4,3) * t185 + t238 * t129 - t235 * t130 + (-t136 * t224 - t137 * t230) * pkin(9);
t245 = qJ(3) * t127 + t114 * t227 + t115 * t221;
t106 = -mrSges(3,1) * t201 + mrSges(3,3) * t190 - pkin(2) * t122 - t225 * t113 + t245 * t231;
t108 = mrSges(3,2) * t201 - mrSges(3,3) * t189 - t221 * t114 + t227 * t115 + (-t122 * t225 - t123 * t231) * qJ(3);
t118 = -t222 * t120 + t228 * t126;
t269 = qJ(2) * t118 + t106 * t228 + t108 * t222;
t121 = m(3) * t201 + t122;
t112 = -t226 * t121 + t270 * t232;
t104 = mrSges(3,1) * t189 - mrSges(3,2) * t190 + pkin(2) * t123 + t231 * t113 + t245 * t225;
t244 = mrSges(2,1) * t215 - mrSges(2,2) * t216 + pkin(1) * t112 + t232 * t104 + t269 * t226;
t116 = m(2) * t216 + t118;
t111 = t232 * t121 + t270 * t226;
t109 = m(2) * t215 + t112;
t102 = mrSges(2,2) * t220 - mrSges(2,3) * t215 - t222 * t106 + t228 * t108 + (-t111 * t226 - t112 * t232) * qJ(2);
t101 = -mrSges(2,1) * t220 + mrSges(2,3) * t216 - pkin(1) * t111 - t226 * t104 + t269 * t232;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t229 * t102 - t223 * t101 - qJ(1) * (t229 * t109 + t223 * t116), t102, t108, t115, t129, t143, t161; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t223 * t102 + t229 * t101 + qJ(1) * (-t223 * t109 + t229 * t116), t101, t106, t114, t130, t144, t160; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t244, t244, t104, t113, t128, t268, t242;];
m_new  = t1;
