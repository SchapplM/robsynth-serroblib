% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:37
% EndTime: 2019-12-31 19:14:46
% DurationCPUTime: 3.94s
% Computational Cost: add. (52635->266), mult. (101536->325), div. (0->0), fcn. (61692->8), ass. (0->110)
t231 = sin(qJ(1));
t235 = cos(qJ(1));
t213 = -t235 * g(1) - t231 * g(2);
t261 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t213;
t260 = (-pkin(1) - pkin(6));
t259 = mrSges(2,1) - mrSges(3,2);
t258 = -Ifges(3,4) + Ifges(2,5);
t257 = (Ifges(3,5) - Ifges(2,6));
t237 = qJD(1) ^ 2;
t183 = (t260 * t237) - t261;
t234 = cos(qJ(3));
t255 = qJD(1) * qJD(3);
t216 = t234 * t255;
t230 = sin(qJ(3));
t207 = -t230 * qJDD(1) - t216;
t253 = t230 * t255;
t208 = qJDD(1) * t234 - t253;
t159 = (-t208 + t253) * pkin(7) + (-t207 + t216) * pkin(3) + t183;
t212 = g(1) * t231 - t235 * g(2);
t249 = -qJ(2) * t237 + qJDD(2) - t212;
t184 = t260 * qJDD(1) + t249;
t177 = -g(3) * t234 + t230 * t184;
t206 = (pkin(3) * t230 - pkin(7) * t234) * qJD(1);
t218 = t230 * qJD(1);
t236 = qJD(3) ^ 2;
t162 = -pkin(3) * t236 + qJDD(3) * pkin(7) - t206 * t218 + t177;
t229 = sin(qJ(4));
t233 = cos(qJ(4));
t144 = t233 * t159 - t162 * t229;
t256 = qJD(1) * t234;
t203 = qJD(3) * t233 - t229 * t256;
t171 = qJD(4) * t203 + qJDD(3) * t229 + t208 * t233;
t202 = qJDD(4) - t207;
t204 = qJD(3) * t229 + t233 * t256;
t215 = t218 + qJD(4);
t141 = (t203 * t215 - t171) * pkin(8) + (t203 * t204 + t202) * pkin(4) + t144;
t145 = t229 * t159 + t233 * t162;
t170 = -qJD(4) * t204 + qJDD(3) * t233 - t208 * t229;
t182 = pkin(4) * t215 - pkin(8) * t204;
t201 = t203 ^ 2;
t142 = -pkin(4) * t201 + pkin(8) * t170 - t182 * t215 + t145;
t228 = sin(qJ(5));
t232 = cos(qJ(5));
t139 = t141 * t232 - t142 * t228;
t172 = t203 * t232 - t204 * t228;
t151 = qJD(5) * t172 + t170 * t228 + t171 * t232;
t173 = t203 * t228 + t204 * t232;
t156 = -mrSges(6,1) * t172 + mrSges(6,2) * t173;
t214 = qJD(5) + t215;
t163 = -mrSges(6,2) * t214 + mrSges(6,3) * t172;
t195 = qJDD(5) + t202;
t136 = m(6) * t139 + mrSges(6,1) * t195 - mrSges(6,3) * t151 - t156 * t173 + t163 * t214;
t140 = t141 * t228 + t142 * t232;
t150 = -qJD(5) * t173 + t170 * t232 - t171 * t228;
t164 = mrSges(6,1) * t214 - mrSges(6,3) * t173;
t137 = m(6) * t140 - mrSges(6,2) * t195 + mrSges(6,3) * t150 + t156 * t172 - t164 * t214;
t127 = t232 * t136 + t228 * t137;
t175 = -mrSges(5,1) * t203 + mrSges(5,2) * t204;
t178 = -mrSges(5,2) * t215 + mrSges(5,3) * t203;
t125 = m(5) * t144 + mrSges(5,1) * t202 - mrSges(5,3) * t171 - t175 * t204 + t178 * t215 + t127;
t179 = mrSges(5,1) * t215 - mrSges(5,3) * t204;
t251 = -t228 * t136 + t232 * t137;
t126 = m(5) * t145 - mrSges(5,2) * t202 + mrSges(5,3) * t170 + t175 * t203 - t179 * t215 + t251;
t120 = t233 * t125 + t229 * t126;
t205 = (mrSges(4,1) * t230 + mrSges(4,2) * t234) * qJD(1);
t211 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t256;
t252 = -t125 * t229 + t233 * t126;
t118 = m(4) * t177 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t207 - qJD(3) * t211 - t205 * t218 + t252;
t176 = g(3) * t230 + t184 * t234;
t210 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t218;
t161 = -qJDD(3) * pkin(3) - pkin(7) * t236 + t206 * t256 - t176;
t143 = -pkin(4) * t170 - pkin(8) * t201 + t182 * t204 + t161;
t247 = m(6) * t143 - t150 * mrSges(6,1) + mrSges(6,2) * t151 - t172 * t163 + t164 * t173;
t239 = -m(5) * t161 + t170 * mrSges(5,1) - mrSges(5,2) * t171 + t203 * t178 - t179 * t204 - t247;
t130 = m(4) * t176 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t208 + qJD(3) * t210 - t205 * t256 + t239;
t111 = t234 * t118 - t130 * t230;
t110 = t118 * t230 + t130 * t234;
t189 = -qJDD(1) * pkin(1) + t249;
t248 = -m(3) * t189 + (t237 * mrSges(3,3)) - t110;
t116 = -m(4) * t183 + mrSges(4,1) * t207 - t208 * mrSges(4,2) - t210 * t218 - t211 * t256 - t120;
t153 = Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t214;
t154 = Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t214;
t246 = -mrSges(6,1) * t139 + mrSges(6,2) * t140 - Ifges(6,5) * t151 - Ifges(6,6) * t150 - Ifges(6,3) * t195 - t173 * t153 + t172 * t154;
t152 = Ifges(6,5) * t173 + Ifges(6,6) * t172 + Ifges(6,3) * t214;
t128 = -mrSges(6,1) * t143 + mrSges(6,3) * t140 + Ifges(6,4) * t151 + Ifges(6,2) * t150 + Ifges(6,6) * t195 - t152 * t173 + t154 * t214;
t129 = mrSges(6,2) * t143 - mrSges(6,3) * t139 + Ifges(6,1) * t151 + Ifges(6,4) * t150 + Ifges(6,5) * t195 + t152 * t172 - t153 * t214;
t165 = Ifges(5,5) * t204 + Ifges(5,6) * t203 + Ifges(5,3) * t215;
t167 = Ifges(5,1) * t204 + Ifges(5,4) * t203 + Ifges(5,5) * t215;
t106 = -mrSges(5,1) * t161 + mrSges(5,3) * t145 + Ifges(5,4) * t171 + Ifges(5,2) * t170 + Ifges(5,6) * t202 - pkin(4) * t247 + pkin(8) * t251 + t232 * t128 + t228 * t129 - t204 * t165 + t215 * t167;
t166 = Ifges(5,4) * t204 + Ifges(5,2) * t203 + Ifges(5,6) * t215;
t113 = mrSges(5,2) * t161 - mrSges(5,3) * t144 + Ifges(5,1) * t171 + Ifges(5,4) * t170 + Ifges(5,5) * t202 - pkin(8) * t127 - t128 * t228 + t129 * t232 + t165 * t203 - t166 * t215;
t192 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t234 - Ifges(4,6) * t230) * qJD(1);
t193 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t234 - Ifges(4,2) * t230) * qJD(1);
t103 = mrSges(4,2) * t183 - mrSges(4,3) * t176 + Ifges(4,1) * t208 + Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3) - pkin(7) * t120 - qJD(3) * t193 - t106 * t229 + t113 * t233 - t192 * t218;
t194 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t234 - Ifges(4,4) * t230) * qJD(1);
t238 = mrSges(5,1) * t144 - mrSges(5,2) * t145 + Ifges(5,5) * t171 + Ifges(5,6) * t170 + Ifges(5,3) * t202 + pkin(4) * t127 + t204 * t166 - t203 * t167 - t246;
t104 = -mrSges(4,1) * t183 + mrSges(4,3) * t177 + Ifges(4,4) * t208 + Ifges(4,2) * t207 + Ifges(4,6) * qJDD(3) - pkin(3) * t120 + qJD(3) * t194 - t192 * t256 - t238;
t187 = pkin(1) * t237 + t261;
t245 = mrSges(3,2) * t189 - mrSges(3,3) * t187 + Ifges(3,1) * qJDD(1) - pkin(6) * t110 + t234 * t103 - t104 * t230;
t244 = -mrSges(3,1) * t187 - pkin(2) * t116 - pkin(6) * t111 - t103 * t230 - t104 * t234;
t243 = mrSges(4,1) * t176 - mrSges(4,2) * t177 + Ifges(4,5) * t208 + Ifges(4,6) * t207 + Ifges(4,3) * qJDD(3) + pkin(3) * t239 + pkin(7) * t252 + t233 * t106 + t229 * t113 + t193 * t256 + t194 * t218;
t242 = -m(3) * t187 + t237 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t116;
t241 = -mrSges(2,2) * t213 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t248) + qJ(2) * t242 + mrSges(2,1) * t212 + Ifges(2,3) * qJDD(1) + t245;
t240 = mrSges(3,1) * t189 + pkin(2) * t110 + t243;
t114 = m(2) * t213 - mrSges(2,1) * t237 - qJDD(1) * mrSges(2,2) + t242;
t109 = -m(3) * g(3) + t111;
t107 = m(2) * t212 - mrSges(2,2) * t237 + t259 * qJDD(1) + t248;
t101 = (t257 * t237) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t258 * qJDD(1) - mrSges(2,3) * t212 - qJ(2) * t109 + t240;
t100 = mrSges(2,3) * t213 - pkin(1) * t109 + t259 * g(3) - t257 * qJDD(1) + t258 * t237 + t244;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t235 * t101 - t231 * t100 - pkin(5) * (t107 * t235 + t114 * t231), t101, t245, t103, t113, t129; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t231 * t101 + t235 * t100 + pkin(5) * (-t107 * t231 + t114 * t235), t100, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (Ifges(3,5) * t237) - t240, t104, t106, t128; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t241, t241, mrSges(3,2) * g(3) + Ifges(3,4) * t237 + Ifges(3,5) * qJDD(1) - t244, t243, t238, -t246;];
m_new = t1;
