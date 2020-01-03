% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:46
% EndTime: 2019-12-31 17:57:53
% DurationCPUTime: 5.03s
% Computational Cost: add. (62990->246), mult. (137501->310), div. (0->0), fcn. (90285->10), ass. (0->111)
t233 = qJD(1) ^ 2;
t222 = sin(pkin(9));
t224 = cos(pkin(9));
t227 = sin(qJ(4));
t230 = cos(qJ(4));
t243 = t222 * t227 - t224 * t230;
t194 = t243 * qJD(1);
t228 = sin(qJ(1));
t231 = cos(qJ(1));
t206 = t228 * g(1) - t231 * g(2);
t203 = qJDD(1) * pkin(1) + t206;
t207 = -t231 * g(1) - t228 * g(2);
t204 = -t233 * pkin(1) + t207;
t223 = sin(pkin(8));
t225 = cos(pkin(8));
t188 = t223 * t203 + t225 * t204;
t176 = -t233 * pkin(2) + qJDD(1) * qJ(3) + t188;
t221 = -g(3) + qJDD(2);
t255 = qJD(1) * qJD(3);
t259 = t224 * t221 - 0.2e1 * t222 * t255;
t262 = pkin(3) * t224;
t160 = (-pkin(6) * qJDD(1) + t233 * t262 - t176) * t222 + t259;
t166 = t222 * t221 + (t176 + 0.2e1 * t255) * t224;
t254 = qJDD(1) * t224;
t217 = t224 ^ 2;
t260 = t217 * t233;
t163 = -pkin(3) * t260 + pkin(6) * t254 + t166;
t152 = t227 * t160 + t230 * t163;
t244 = t222 * t230 + t224 * t227;
t195 = t244 * qJD(1);
t178 = t194 * mrSges(5,1) + t195 * mrSges(5,2);
t256 = t195 * qJD(4);
t184 = -t243 * qJDD(1) - t256;
t192 = qJD(4) * mrSges(5,1) - t195 * mrSges(5,3);
t183 = t194 * pkin(4) - t195 * pkin(7);
t232 = qJD(4) ^ 2;
t148 = -t232 * pkin(4) + qJDD(4) * pkin(7) - t194 * t183 + t152;
t216 = t222 ^ 2;
t187 = t225 * t203 - t223 * t204;
t246 = qJDD(3) - t187;
t164 = (-pkin(2) - t262) * qJDD(1) + (-qJ(3) + (-t216 - t217) * pkin(6)) * t233 + t246;
t257 = t194 * qJD(4);
t185 = t244 * qJDD(1) - t257;
t149 = (-t185 + t257) * pkin(7) + (-t184 + t256) * pkin(4) + t164;
t226 = sin(qJ(5));
t229 = cos(qJ(5));
t145 = -t226 * t148 + t229 * t149;
t189 = t229 * qJD(4) - t226 * t195;
t162 = t189 * qJD(5) + t226 * qJDD(4) + t229 * t185;
t190 = t226 * qJD(4) + t229 * t195;
t168 = -t189 * mrSges(6,1) + t190 * mrSges(6,2);
t193 = qJD(5) + t194;
t169 = -t193 * mrSges(6,2) + t189 * mrSges(6,3);
t182 = qJDD(5) - t184;
t141 = m(6) * t145 + t182 * mrSges(6,1) - t162 * mrSges(6,3) - t190 * t168 + t193 * t169;
t146 = t229 * t148 + t226 * t149;
t161 = -t190 * qJD(5) + t229 * qJDD(4) - t226 * t185;
t170 = t193 * mrSges(6,1) - t190 * mrSges(6,3);
t142 = m(6) * t146 - t182 * mrSges(6,2) + t161 * mrSges(6,3) + t189 * t168 - t193 * t170;
t250 = -t226 * t141 + t229 * t142;
t128 = m(5) * t152 - qJDD(4) * mrSges(5,2) + t184 * mrSges(5,3) - qJD(4) * t192 - t194 * t178 + t250;
t151 = t230 * t160 - t227 * t163;
t191 = -qJD(4) * mrSges(5,2) - t194 * mrSges(5,3);
t147 = -qJDD(4) * pkin(4) - t232 * pkin(7) + t195 * t183 - t151;
t240 = -m(6) * t147 + t161 * mrSges(6,1) - t162 * mrSges(6,2) + t189 * t169 - t190 * t170;
t137 = m(5) * t151 + qJDD(4) * mrSges(5,1) - t185 * mrSges(5,3) + qJD(4) * t191 - t195 * t178 + t240;
t122 = t227 * t128 + t230 * t137;
t165 = -t222 * t176 + t259;
t154 = Ifges(6,5) * t190 + Ifges(6,6) * t189 + Ifges(6,3) * t193;
t156 = Ifges(6,1) * t190 + Ifges(6,4) * t189 + Ifges(6,5) * t193;
t134 = -mrSges(6,1) * t147 + mrSges(6,3) * t146 + Ifges(6,4) * t162 + Ifges(6,2) * t161 + Ifges(6,6) * t182 - t190 * t154 + t193 * t156;
t155 = Ifges(6,4) * t190 + Ifges(6,2) * t189 + Ifges(6,6) * t193;
t135 = mrSges(6,2) * t147 - mrSges(6,3) * t145 + Ifges(6,1) * t162 + Ifges(6,4) * t161 + Ifges(6,5) * t182 + t189 * t154 - t193 * t155;
t174 = Ifges(5,4) * t195 - Ifges(5,2) * t194 + Ifges(5,6) * qJD(4);
t175 = Ifges(5,1) * t195 - Ifges(5,4) * t194 + Ifges(5,5) * qJD(4);
t237 = -mrSges(5,1) * t151 + mrSges(5,2) * t152 - Ifges(5,5) * t185 - Ifges(5,6) * t184 - Ifges(5,3) * qJDD(4) - pkin(4) * t240 - pkin(7) * t250 - t229 * t134 - t226 * t135 - t195 * t174 - t194 * t175;
t248 = Ifges(4,4) * t222 + Ifges(4,2) * t224;
t249 = Ifges(4,1) * t222 + Ifges(4,4) * t224;
t263 = -mrSges(4,1) * t165 + mrSges(4,2) * t166 - pkin(3) * t122 - (t222 * t248 - t224 * t249) * t233 + t237;
t261 = mrSges(4,2) * t222;
t242 = mrSges(4,3) * qJDD(1) + t233 * (-mrSges(4,1) * t224 + t261);
t120 = m(4) * t165 - t242 * t222 + t122;
t251 = t230 * t128 - t227 * t137;
t121 = m(4) * t166 + t242 * t224 + t251;
t252 = -t222 * t120 + t224 * t121;
t112 = m(3) * t188 - t233 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t252;
t172 = -qJDD(1) * pkin(2) - t233 * qJ(3) + t246;
t130 = t229 * t141 + t226 * t142;
t239 = m(5) * t164 - t184 * mrSges(5,1) + t185 * mrSges(5,2) + t194 * t191 + t195 * t192 + t130;
t236 = -m(4) * t172 + mrSges(4,1) * t254 - t239 + (t216 * t233 + t260) * mrSges(4,3);
t124 = t236 + (mrSges(3,1) - t261) * qJDD(1) - t233 * mrSges(3,2) + m(3) * t187;
t109 = t223 * t112 + t225 * t124;
t114 = t224 * t120 + t222 * t121;
t247 = Ifges(4,5) * t222 + Ifges(4,6) * t224;
t258 = t233 * t247;
t253 = t225 * t112 - t223 * t124;
t173 = Ifges(5,5) * t195 - Ifges(5,6) * t194 + Ifges(5,3) * qJD(4);
t115 = mrSges(5,2) * t164 - mrSges(5,3) * t151 + Ifges(5,1) * t185 + Ifges(5,4) * t184 + Ifges(5,5) * qJDD(4) - pkin(7) * t130 - qJD(4) * t174 - t226 * t134 + t229 * t135 - t194 * t173;
t235 = mrSges(6,1) * t145 - mrSges(6,2) * t146 + Ifges(6,5) * t162 + Ifges(6,6) * t161 + Ifges(6,3) * t182 + t190 * t155 - t189 * t156;
t116 = -mrSges(5,1) * t164 + mrSges(5,3) * t152 + Ifges(5,4) * t185 + Ifges(5,2) * t184 + Ifges(5,6) * qJDD(4) - pkin(4) * t130 + qJD(4) * t175 - t195 * t173 - t235;
t103 = -mrSges(4,1) * t172 + mrSges(4,3) * t166 - pkin(3) * t239 + pkin(6) * t251 + t248 * qJDD(1) + t227 * t115 + t230 * t116 - t222 * t258;
t105 = mrSges(4,2) * t172 - mrSges(4,3) * t165 - pkin(6) * t122 + t249 * qJDD(1) + t230 * t115 - t227 * t116 + t224 * t258;
t241 = -mrSges(3,2) * t188 + qJ(3) * t252 + t224 * t103 + t222 * t105 + pkin(2) * (-qJDD(1) * t261 + t236) + mrSges(3,1) * t187 + Ifges(3,3) * qJDD(1);
t238 = mrSges(2,1) * t206 - mrSges(2,2) * t207 + Ifges(2,3) * qJDD(1) + pkin(1) * t109 + t241;
t107 = m(2) * t207 - t233 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t253;
t106 = m(2) * t206 + qJDD(1) * mrSges(2,1) - t233 * mrSges(2,2) + t109;
t101 = (Ifges(3,6) - t247) * qJDD(1) + t233 * Ifges(3,5) - mrSges(3,1) * t221 + mrSges(3,3) * t188 - pkin(2) * t114 + t263;
t100 = mrSges(3,2) * t221 - mrSges(3,3) * t187 + Ifges(3,5) * qJDD(1) - t233 * Ifges(3,6) - qJ(3) * t114 - t222 * t103 + t224 * t105;
t99 = -mrSges(2,2) * g(3) - mrSges(2,3) * t206 + Ifges(2,5) * qJDD(1) - t233 * Ifges(2,6) - qJ(2) * t109 + t225 * t100 - t223 * t101;
t98 = Ifges(2,6) * qJDD(1) + t233 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t207 + t223 * t100 + t225 * t101 - pkin(1) * (m(3) * t221 + t114) + qJ(2) * t253;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t231 * t99 - t228 * t98 - pkin(5) * (t231 * t106 + t228 * t107), t99, t100, t105, t115, t135; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t228 * t99 + t231 * t98 + pkin(5) * (-t228 * t106 + t231 * t107), t98, t101, t103, t116, t134; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t238, t238, t241, t247 * qJDD(1) - t263, -t237, t235;];
m_new = t1;
