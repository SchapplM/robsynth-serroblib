% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:25
% EndTime: 2019-12-31 18:21:34
% DurationCPUTime: 5.88s
% Computational Cost: add. (77448->266), mult. (158256->338), div. (0->0), fcn. (97852->10), ass. (0->108)
t226 = sin(qJ(1));
t229 = cos(qJ(1));
t211 = t226 * g(1) - t229 * g(2);
t202 = qJDD(1) * pkin(1) + t211;
t212 = -t229 * g(1) - t226 * g(2);
t231 = qJD(1) ^ 2;
t205 = -t231 * pkin(1) + t212;
t221 = sin(pkin(8));
t223 = cos(pkin(8));
t177 = t223 * t202 - t221 * t205;
t168 = -qJDD(1) * pkin(2) - t231 * pkin(6) - t177;
t225 = sin(qJ(3));
t228 = cos(qJ(3));
t244 = qJD(1) * qJD(3);
t243 = t228 * t244;
t206 = t225 * qJDD(1) + t243;
t215 = t225 * t244;
t207 = t228 * qJDD(1) - t215;
t157 = (-t206 - t243) * qJ(4) + (-t207 + t215) * pkin(3) + t168;
t178 = t221 * t202 + t223 * t205;
t169 = -t231 * pkin(2) + qJDD(1) * pkin(6) + t178;
t219 = -g(3) + qJDD(2);
t163 = t228 * t169 + t225 * t219;
t203 = (-pkin(3) * t228 - qJ(4) * t225) * qJD(1);
t230 = qJD(3) ^ 2;
t245 = t228 * qJD(1);
t161 = -t230 * pkin(3) + qJDD(3) * qJ(4) + t203 * t245 + t163;
t220 = sin(pkin(9));
t222 = cos(pkin(9));
t246 = qJD(1) * t225;
t199 = t220 * qJD(3) + t222 * t246;
t142 = -0.2e1 * qJD(4) * t199 + t222 * t157 - t220 * t161;
t184 = t220 * qJDD(3) + t222 * t206;
t198 = t222 * qJD(3) - t220 * t246;
t140 = (-t198 * t245 - t184) * pkin(7) + (t198 * t199 - t207) * pkin(4) + t142;
t143 = 0.2e1 * qJD(4) * t198 + t220 * t157 + t222 * t161;
t183 = t222 * qJDD(3) - t220 * t206;
t185 = -pkin(4) * t245 - t199 * pkin(7);
t197 = t198 ^ 2;
t141 = -t197 * pkin(4) + t183 * pkin(7) + t185 * t245 + t143;
t224 = sin(qJ(5));
t227 = cos(qJ(5));
t139 = t224 * t140 + t227 * t141;
t162 = -t225 * t169 + t228 * t219;
t159 = -qJDD(3) * pkin(3) - t230 * qJ(4) + t203 * t246 + qJDD(4) - t162;
t144 = -t183 * pkin(4) - t197 * pkin(7) + t199 * t185 + t159;
t175 = t224 * t198 + t227 * t199;
t149 = -t175 * qJD(5) + t227 * t183 - t224 * t184;
t174 = t227 * t198 - t224 * t199;
t150 = t174 * qJD(5) + t224 * t183 + t227 * t184;
t213 = qJD(5) - t245;
t153 = Ifges(6,5) * t175 + Ifges(6,6) * t174 + Ifges(6,3) * t213;
t155 = Ifges(6,1) * t175 + Ifges(6,4) * t174 + Ifges(6,5) * t213;
t201 = qJDD(5) - t207;
t127 = -mrSges(6,1) * t144 + mrSges(6,3) * t139 + Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t201 - t175 * t153 + t213 * t155;
t138 = t227 * t140 - t224 * t141;
t154 = Ifges(6,4) * t175 + Ifges(6,2) * t174 + Ifges(6,6) * t213;
t128 = mrSges(6,2) * t144 - mrSges(6,3) * t138 + Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t201 + t174 * t153 - t213 * t154;
t170 = Ifges(5,5) * t199 + Ifges(5,6) * t198 - Ifges(5,3) * t245;
t172 = Ifges(5,1) * t199 + Ifges(5,4) * t198 - Ifges(5,5) * t245;
t164 = -t213 * mrSges(6,2) + t174 * mrSges(6,3);
t165 = t213 * mrSges(6,1) - t175 * mrSges(6,3);
t237 = m(6) * t144 - t149 * mrSges(6,1) + t150 * mrSges(6,2) - t174 * t164 + t175 * t165;
t160 = -t174 * mrSges(6,1) + t175 * mrSges(6,2);
t134 = m(6) * t138 + t201 * mrSges(6,1) - t150 * mrSges(6,3) - t175 * t160 + t213 * t164;
t135 = m(6) * t139 - t201 * mrSges(6,2) + t149 * mrSges(6,3) + t174 * t160 - t213 * t165;
t240 = -t224 * t134 + t227 * t135;
t112 = -mrSges(5,1) * t159 + mrSges(5,3) * t143 + Ifges(5,4) * t184 + Ifges(5,2) * t183 - Ifges(5,6) * t207 - pkin(4) * t237 + pkin(7) * t240 + t227 * t127 + t224 * t128 - t199 * t170 - t172 * t245;
t126 = t227 * t134 + t224 * t135;
t171 = Ifges(5,4) * t199 + Ifges(5,2) * t198 - Ifges(5,6) * t245;
t114 = mrSges(5,2) * t159 - mrSges(5,3) * t142 + Ifges(5,1) * t184 + Ifges(5,4) * t183 - Ifges(5,5) * t207 - pkin(7) * t126 - t224 * t127 + t227 * t128 + t198 * t170 + t171 * t245;
t179 = -t198 * mrSges(5,1) + t199 * mrSges(5,2);
t181 = mrSges(5,2) * t245 + t198 * mrSges(5,3);
t124 = m(5) * t142 - t207 * mrSges(5,1) - t184 * mrSges(5,3) - t199 * t179 - t181 * t245 + t126;
t182 = -mrSges(5,1) * t245 - t199 * mrSges(5,3);
t125 = m(5) * t143 + t207 * mrSges(5,2) + t183 * mrSges(5,3) + t198 * t179 + t182 * t245 + t240;
t122 = -t220 * t124 + t222 * t125;
t136 = -m(5) * t159 + t183 * mrSges(5,1) - t184 * mrSges(5,2) + t198 * t181 - t199 * t182 - t237;
t193 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t225 + Ifges(4,2) * t228) * qJD(1);
t194 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t225 + Ifges(4,4) * t228) * qJD(1);
t247 = mrSges(4,1) * t162 - mrSges(4,2) * t163 + Ifges(4,5) * t206 + Ifges(4,6) * t207 + Ifges(4,3) * qJDD(3) + pkin(3) * t136 + qJ(4) * t122 + t222 * t112 + t220 * t114 + (t225 * t193 - t228 * t194) * qJD(1);
t204 = (-mrSges(4,1) * t228 + mrSges(4,2) * t225) * qJD(1);
t209 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t246;
t120 = m(4) * t163 - qJDD(3) * mrSges(4,2) + t207 * mrSges(4,3) - qJD(3) * t209 + t204 * t245 + t122;
t210 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t245;
t130 = m(4) * t162 + qJDD(3) * mrSges(4,1) - t206 * mrSges(4,3) + qJD(3) * t210 - t204 * t246 + t136;
t241 = t228 * t120 - t225 * t130;
t110 = m(3) * t178 - t231 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t241;
t121 = t222 * t124 + t220 * t125;
t234 = -m(4) * t168 + t207 * mrSges(4,1) - t206 * mrSges(4,2) - t209 * t246 + t210 * t245 - t121;
t116 = m(3) * t177 + qJDD(1) * mrSges(3,1) - t231 * mrSges(3,2) + t234;
t105 = t221 * t110 + t223 * t116;
t113 = t225 * t120 + t228 * t130;
t242 = t223 * t110 - t221 * t116;
t192 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t225 + Ifges(4,6) * t228) * qJD(1);
t101 = mrSges(4,2) * t168 - mrSges(4,3) * t162 + Ifges(4,1) * t206 + Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3) - qJ(4) * t121 - qJD(3) * t193 - t220 * t112 + t222 * t114 + t192 * t245;
t236 = -mrSges(6,1) * t138 + mrSges(6,2) * t139 - Ifges(6,5) * t150 - Ifges(6,6) * t149 - Ifges(6,3) * t201 - t175 * t154 + t174 * t155;
t232 = -mrSges(5,1) * t142 + mrSges(5,2) * t143 - Ifges(5,5) * t184 - Ifges(5,6) * t183 - pkin(4) * t126 - t199 * t171 + t198 * t172 + t236;
t107 = t232 + Ifges(4,6) * qJDD(3) - t192 * t246 + (Ifges(4,2) + Ifges(5,3)) * t207 + Ifges(4,4) * t206 + qJD(3) * t194 + mrSges(4,3) * t163 - mrSges(4,1) * t168 - pkin(3) * t121;
t238 = mrSges(3,1) * t177 - mrSges(3,2) * t178 + Ifges(3,3) * qJDD(1) + pkin(2) * t234 + pkin(6) * t241 + t225 * t101 + t228 * t107;
t235 = mrSges(2,1) * t211 - mrSges(2,2) * t212 + Ifges(2,3) * qJDD(1) + pkin(1) * t105 + t238;
t103 = m(2) * t212 - t231 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t242;
t102 = m(2) * t211 + qJDD(1) * mrSges(2,1) - t231 * mrSges(2,2) + t105;
t99 = -mrSges(3,1) * t219 + mrSges(3,3) * t178 + t231 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t113 - t247;
t98 = mrSges(3,2) * t219 - mrSges(3,3) * t177 + Ifges(3,5) * qJDD(1) - t231 * Ifges(3,6) - pkin(6) * t113 + t228 * t101 - t225 * t107;
t97 = -mrSges(2,2) * g(3) - mrSges(2,3) * t211 + Ifges(2,5) * qJDD(1) - t231 * Ifges(2,6) - qJ(2) * t105 - t221 * t99 + t223 * t98;
t96 = Ifges(2,6) * qJDD(1) + t231 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t212 + t221 * t98 + t223 * t99 - pkin(1) * (m(3) * t219 + t113) + qJ(2) * t242;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t229 * t97 - t226 * t96 - pkin(5) * (t229 * t102 + t226 * t103), t97, t98, t101, t114, t128; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t226 * t97 + t229 * t96 + pkin(5) * (-t226 * t102 + t229 * t103), t96, t99, t107, t112, t127; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t235, t235, t238, t247, -Ifges(5,3) * t207 - t232, -t236;];
m_new = t1;
