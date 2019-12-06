% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:46
% EndTime: 2019-12-05 15:42:53
% DurationCPUTime: 6.22s
% Computational Cost: add. (75605->236), mult. (172426->300), div. (0->0), fcn. (122991->10), ass. (0->108)
t233 = qJD(2) ^ 2;
t224 = sin(pkin(8));
t226 = cos(pkin(8));
t206 = t224 * g(1) - t226 * g(2);
t207 = -t226 * g(1) - t224 * g(2);
t229 = sin(qJ(2));
t232 = cos(qJ(2));
t191 = t229 * t206 + t232 * t207;
t188 = -t233 * pkin(2) + qJDD(2) * qJ(3) + t191;
t223 = sin(pkin(9));
t222 = -g(3) + qJDD(1);
t225 = cos(pkin(9));
t255 = qJD(2) * qJD(3);
t258 = t225 * t222 - 0.2e1 * t223 * t255;
t261 = pkin(3) * t225;
t166 = (-pkin(6) * qJDD(2) + t233 * t261 - t188) * t223 + t258;
t171 = t223 * t222 + (t188 + 0.2e1 * t255) * t225;
t254 = qJDD(2) * t225;
t218 = t225 ^ 2;
t259 = t218 * t233;
t167 = -pkin(3) * t259 + pkin(6) * t254 + t171;
t228 = sin(qJ(4));
t231 = cos(qJ(4));
t148 = t231 * t166 - t228 * t167;
t243 = t223 * t231 + t225 * t228;
t242 = -t223 * t228 + t225 * t231;
t196 = t242 * qJD(2);
t256 = t196 * qJD(4);
t187 = t243 * qJDD(2) + t256;
t197 = t243 * qJD(2);
t143 = (-t187 + t256) * pkin(7) + (t196 * t197 + qJDD(4)) * pkin(4) + t148;
t149 = t228 * t166 + t231 * t167;
t186 = -t197 * qJD(4) + t242 * qJDD(2);
t194 = qJD(4) * pkin(4) - t197 * pkin(7);
t195 = t196 ^ 2;
t144 = -t195 * pkin(4) + t186 * pkin(7) - qJD(4) * t194 + t149;
t227 = sin(qJ(5));
t230 = cos(qJ(5));
t141 = t230 * t143 - t227 * t144;
t177 = t230 * t196 - t227 * t197;
t156 = t177 * qJD(5) + t227 * t186 + t230 * t187;
t178 = t227 * t196 + t230 * t197;
t162 = -t177 * mrSges(6,1) + t178 * mrSges(6,2);
t219 = qJD(4) + qJD(5);
t172 = -t219 * mrSges(6,2) + t177 * mrSges(6,3);
t216 = qJDD(4) + qJDD(5);
t138 = m(6) * t141 + t216 * mrSges(6,1) - t156 * mrSges(6,3) - t178 * t162 + t219 * t172;
t142 = t227 * t143 + t230 * t144;
t155 = -t178 * qJD(5) + t230 * t186 - t227 * t187;
t173 = t219 * mrSges(6,1) - t178 * mrSges(6,3);
t139 = m(6) * t142 - t216 * mrSges(6,2) + t155 * mrSges(6,3) + t177 * t162 - t219 * t173;
t129 = t230 * t138 + t227 * t139;
t181 = -t196 * mrSges(5,1) + t197 * mrSges(5,2);
t192 = -qJD(4) * mrSges(5,2) + t196 * mrSges(5,3);
t126 = m(5) * t148 + qJDD(4) * mrSges(5,1) - t187 * mrSges(5,3) + qJD(4) * t192 - t197 * t181 + t129;
t193 = qJD(4) * mrSges(5,1) - t197 * mrSges(5,3);
t250 = -t227 * t138 + t230 * t139;
t127 = m(5) * t149 - qJDD(4) * mrSges(5,2) + t186 * mrSges(5,3) - qJD(4) * t193 + t196 * t181 + t250;
t122 = t231 * t126 + t228 * t127;
t170 = -t223 * t188 + t258;
t175 = Ifges(5,4) * t197 + Ifges(5,2) * t196 + Ifges(5,6) * qJD(4);
t176 = Ifges(5,1) * t197 + Ifges(5,4) * t196 + Ifges(5,5) * qJD(4);
t158 = Ifges(6,4) * t178 + Ifges(6,2) * t177 + Ifges(6,6) * t219;
t159 = Ifges(6,1) * t178 + Ifges(6,4) * t177 + Ifges(6,5) * t219;
t239 = -mrSges(6,1) * t141 + mrSges(6,2) * t142 - Ifges(6,5) * t156 - Ifges(6,6) * t155 - Ifges(6,3) * t216 - t178 * t158 + t177 * t159;
t235 = -mrSges(5,1) * t148 + mrSges(5,2) * t149 - Ifges(5,5) * t187 - Ifges(5,6) * t186 - Ifges(5,3) * qJDD(4) - pkin(4) * t129 - t197 * t175 + t196 * t176 + t239;
t248 = Ifges(4,4) * t223 + Ifges(4,2) * t225;
t249 = Ifges(4,1) * t223 + Ifges(4,4) * t225;
t262 = -mrSges(4,1) * t170 + mrSges(4,2) * t171 - pkin(3) * t122 - (t223 * t248 - t225 * t249) * t233 + t235;
t260 = mrSges(4,2) * t223;
t241 = mrSges(4,3) * qJDD(2) + t233 * (-mrSges(4,1) * t225 + t260);
t120 = m(4) * t170 - t241 * t223 + t122;
t251 = -t228 * t126 + t231 * t127;
t121 = m(4) * t171 + t241 * t225 + t251;
t252 = -t223 * t120 + t225 * t121;
t112 = m(3) * t191 - t233 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t252;
t190 = t232 * t206 - t229 * t207;
t246 = qJDD(3) - t190;
t185 = -qJDD(2) * pkin(2) - t233 * qJ(3) + t246;
t217 = t223 ^ 2;
t169 = (-pkin(2) - t261) * qJDD(2) + (-qJ(3) + (-t217 - t218) * pkin(6)) * t233 + t246;
t146 = -t186 * pkin(4) - t195 * pkin(7) + t197 * t194 + t169;
t245 = m(6) * t146 - t155 * mrSges(6,1) + t156 * mrSges(6,2) - t177 * t172 + t178 * t173;
t237 = m(5) * t169 - t186 * mrSges(5,1) + t187 * mrSges(5,2) - t196 * t192 + t197 * t193 + t245;
t236 = -m(4) * t185 + mrSges(4,1) * t254 - t237 + (t217 * t233 + t259) * mrSges(4,3);
t133 = t236 - t233 * mrSges(3,2) + m(3) * t190 + (mrSges(3,1) - t260) * qJDD(2);
t109 = t229 * t112 + t232 * t133;
t114 = t225 * t120 + t223 * t121;
t247 = Ifges(4,5) * t223 + Ifges(4,6) * t225;
t257 = t233 * t247;
t253 = t232 * t112 - t229 * t133;
t157 = Ifges(6,5) * t178 + Ifges(6,6) * t177 + Ifges(6,3) * t219;
t130 = -mrSges(6,1) * t146 + mrSges(6,3) * t142 + Ifges(6,4) * t156 + Ifges(6,2) * t155 + Ifges(6,6) * t216 - t178 * t157 + t219 * t159;
t131 = mrSges(6,2) * t146 - mrSges(6,3) * t141 + Ifges(6,1) * t156 + Ifges(6,4) * t155 + Ifges(6,5) * t216 + t177 * t157 - t219 * t158;
t174 = Ifges(5,5) * t197 + Ifges(5,6) * t196 + Ifges(5,3) * qJD(4);
t115 = -mrSges(5,1) * t169 + mrSges(5,3) * t149 + Ifges(5,4) * t187 + Ifges(5,2) * t186 + Ifges(5,6) * qJDD(4) - pkin(4) * t245 + pkin(7) * t250 + qJD(4) * t176 + t230 * t130 + t227 * t131 - t197 * t174;
t116 = mrSges(5,2) * t169 - mrSges(5,3) * t148 + Ifges(5,1) * t187 + Ifges(5,4) * t186 + Ifges(5,5) * qJDD(4) - pkin(7) * t129 - qJD(4) * t175 - t227 * t130 + t230 * t131 + t196 * t174;
t103 = -mrSges(4,1) * t185 + mrSges(4,3) * t171 - pkin(3) * t237 + pkin(6) * t251 + t248 * qJDD(2) + t231 * t115 + t228 * t116 - t223 * t257;
t105 = mrSges(4,2) * t185 - mrSges(4,3) * t170 - pkin(6) * t122 + t249 * qJDD(2) - t228 * t115 + t231 * t116 + t225 * t257;
t240 = -mrSges(3,2) * t191 + qJ(3) * t252 + t225 * t103 + t223 * t105 + pkin(2) * (-qJDD(2) * t260 + t236) + mrSges(3,1) * t190 + Ifges(3,3) * qJDD(2);
t238 = mrSges(2,1) * t206 - mrSges(2,2) * t207 + pkin(1) * t109 + t240;
t107 = m(2) * t207 + t253;
t106 = m(2) * t206 + t109;
t101 = t233 * Ifges(3,5) - mrSges(3,1) * t222 + mrSges(3,3) * t191 - pkin(2) * t114 + (Ifges(3,6) - t247) * qJDD(2) + t262;
t100 = mrSges(3,2) * t222 - mrSges(3,3) * t190 + Ifges(3,5) * qJDD(2) - t233 * Ifges(3,6) - qJ(3) * t114 - t223 * t103 + t225 * t105;
t99 = mrSges(2,2) * t222 - mrSges(2,3) * t206 - pkin(5) * t109 + t232 * t100 - t229 * t101;
t98 = -mrSges(2,1) * t222 + mrSges(2,3) * t207 + t229 * t100 + t232 * t101 - pkin(1) * (m(3) * t222 + t114) + pkin(5) * t253;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t226 * t99 - t224 * t98 - qJ(1) * (t226 * t106 + t224 * t107), t99, t100, t105, t116, t131; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t224 * t99 + t226 * t98 + qJ(1) * (-t224 * t106 + t226 * t107), t98, t101, t103, t115, t130; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t238, t238, t240, t247 * qJDD(2) - t262, -t235, -t239;];
m_new = t1;
