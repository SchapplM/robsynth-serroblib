% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:34
% EndTime: 2019-12-31 18:58:39
% DurationCPUTime: 1.92s
% Computational Cost: add. (21199->262), mult. (39725->307), div. (0->0), fcn. (21506->6), ass. (0->99)
t224 = sin(qJ(1));
t226 = cos(qJ(1));
t208 = -t226 * g(1) - t224 * g(2);
t262 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t208;
t228 = qJD(1) ^ 2;
t260 = (-pkin(1) - pkin(6));
t177 = (t260 * t228) - t262;
t223 = sin(qJ(3));
t225 = cos(qJ(3));
t249 = qJD(1) * qJD(3);
t245 = t225 * t249;
t202 = -t223 * qJDD(1) - t245;
t246 = t223 * t249;
t203 = t225 * qJDD(1) - t246;
t141 = (-t203 + t246) * pkin(7) + (-t202 + t245) * pkin(3) + t177;
t207 = t224 * g(1) - t226 * g(2);
t240 = -t228 * qJ(2) + qJDD(2) - t207;
t178 = t260 * qJDD(1) + t240;
t170 = -t225 * g(3) + t223 * t178;
t201 = (pkin(3) * t223 - pkin(7) * t225) * qJD(1);
t227 = qJD(3) ^ 2;
t250 = t223 * qJD(1);
t145 = -t227 * pkin(3) + qJDD(3) * pkin(7) - t201 * t250 + t170;
t222 = sin(qJ(4));
t259 = cos(qJ(4));
t138 = t259 * t141 - t222 * t145;
t139 = t222 * t141 + t259 * t145;
t251 = qJD(1) * t225;
t198 = -t259 * qJD(3) + t222 * t251;
t199 = t222 * qJD(3) + t259 * t251;
t210 = qJD(4) + t250;
t147 = Ifges(6,5) * t199 + Ifges(6,6) * t210 + Ifges(6,3) * t198;
t150 = Ifges(5,4) * t199 - Ifges(5,2) * t198 + Ifges(5,6) * t210;
t152 = Ifges(5,1) * t199 - Ifges(5,4) * t198 + Ifges(5,5) * t210;
t161 = t199 * qJD(4) - t259 * qJDD(3) + t222 * t203;
t162 = -t198 * qJD(4) + t222 * qJDD(3) + t259 * t203;
t167 = t198 * mrSges(6,1) - t199 * mrSges(6,3);
t197 = qJDD(4) - t202;
t166 = t198 * pkin(4) - t199 * qJ(5);
t209 = t210 ^ 2;
t134 = -t209 * pkin(4) + t197 * qJ(5) + 0.2e1 * qJD(5) * t210 - t198 * t166 + t139;
t136 = -t197 * pkin(4) - t209 * qJ(5) + t199 * t166 + qJDD(5) - t138;
t151 = Ifges(6,1) * t199 + Ifges(6,4) * t210 + Ifges(6,5) * t198;
t239 = mrSges(6,1) * t136 - mrSges(6,3) * t134 - Ifges(6,4) * t162 - Ifges(6,2) * t197 - Ifges(6,6) * t161 - t198 * t151;
t174 = -t198 * mrSges(6,2) + t210 * mrSges(6,3);
t243 = -m(6) * t136 + t197 * mrSges(6,1) + t210 * t174;
t173 = -t210 * mrSges(6,1) + t199 * mrSges(6,2);
t247 = m(6) * t134 + t197 * mrSges(6,3) + t210 * t173;
t261 = -(-t150 + t147) * t199 + mrSges(5,1) * t138 - mrSges(5,2) * t139 + Ifges(5,5) * t162 - Ifges(5,6) * t161 + Ifges(5,3) * t197 + pkin(4) * (-t162 * mrSges(6,2) - t199 * t167 + t243) + qJ(5) * (-t161 * mrSges(6,2) - t198 * t167 + t247) + t198 * t152 - t239;
t258 = mrSges(2,1) - mrSges(3,2);
t257 = -mrSges(5,3) - mrSges(6,2);
t256 = -Ifges(3,4) + Ifges(2,5);
t255 = (Ifges(3,5) - Ifges(2,6));
t172 = t210 * mrSges(5,1) - t199 * mrSges(5,3);
t252 = -t198 * mrSges(5,1) - t199 * mrSges(5,2) - t167;
t124 = m(5) * t139 - t197 * mrSges(5,2) + t257 * t161 - t210 * t172 + t252 * t198 + t247;
t171 = -t210 * mrSges(5,2) - t198 * mrSges(5,3);
t126 = m(5) * t138 + t197 * mrSges(5,1) + t257 * t162 + t210 * t171 + t252 * t199 + t243;
t119 = t222 * t124 + t259 * t126;
t149 = Ifges(6,4) * t199 + Ifges(6,2) * t210 + Ifges(6,6) * t198;
t254 = -Ifges(5,5) * t199 + Ifges(5,6) * t198 - Ifges(5,3) * t210 - t149;
t200 = (mrSges(4,1) * t223 + mrSges(4,2) * t225) * qJD(1);
t206 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t251;
t244 = t259 * t124 - t222 * t126;
t117 = m(4) * t170 - qJDD(3) * mrSges(4,2) + t202 * mrSges(4,3) - qJD(3) * t206 - t200 * t250 + t244;
t169 = t223 * g(3) + t225 * t178;
t205 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t250;
t144 = -qJDD(3) * pkin(3) - t227 * pkin(7) + t201 * t251 - t169;
t137 = -0.2e1 * qJD(5) * t199 + (t198 * t210 - t162) * qJ(5) + (t199 * t210 + t161) * pkin(4) + t144;
t129 = m(6) * t137 + t161 * mrSges(6,1) - t162 * mrSges(6,3) - t199 * t173 + t198 * t174;
t230 = -m(5) * t144 - t161 * mrSges(5,1) - t162 * mrSges(5,2) - t198 * t171 - t199 * t172 - t129;
t121 = m(4) * t169 + qJDD(3) * mrSges(4,1) - t203 * mrSges(4,3) + qJD(3) * t205 - t200 * t251 + t230;
t108 = t225 * t117 - t223 * t121;
t242 = -mrSges(6,1) * t137 + mrSges(6,2) * t134;
t107 = t223 * t117 + t225 * t121;
t238 = mrSges(6,2) * t136 - mrSges(6,3) * t137 + Ifges(6,1) * t162 + Ifges(6,4) * t197 + Ifges(6,5) * t161 + t210 * t147;
t183 = -qJDD(1) * pkin(1) + t240;
t237 = -m(3) * t183 + (t228 * mrSges(3,3)) - t107;
t115 = -m(4) * t177 + t202 * mrSges(4,1) - t203 * mrSges(4,2) - t205 * t250 - t206 * t251 - t119;
t111 = -mrSges(5,1) * t144 + mrSges(5,3) * t139 - pkin(4) * t129 + (t151 + t152) * t210 + t254 * t199 + (Ifges(5,6) - Ifges(6,6)) * t197 + (Ifges(5,4) - Ifges(6,5)) * t162 + (-Ifges(5,2) - Ifges(6,3)) * t161 + t242;
t113 = mrSges(5,2) * t144 - mrSges(5,3) * t138 + Ifges(5,1) * t162 - Ifges(5,4) * t161 + Ifges(5,5) * t197 - qJ(5) * t129 - t210 * t150 + t254 * t198 + t238;
t185 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t225 - Ifges(4,6) * t223) * qJD(1);
t186 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t225 - Ifges(4,2) * t223) * qJD(1);
t102 = mrSges(4,2) * t177 - mrSges(4,3) * t169 + Ifges(4,1) * t203 + Ifges(4,4) * t202 + Ifges(4,5) * qJDD(3) - pkin(7) * t119 - qJD(3) * t186 - t222 * t111 + t259 * t113 - t185 * t250;
t187 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t225 - Ifges(4,4) * t223) * qJD(1);
t103 = -mrSges(4,1) * t177 + mrSges(4,3) * t170 + Ifges(4,4) * t203 + Ifges(4,2) * t202 + Ifges(4,6) * qJDD(3) - pkin(3) * t119 + qJD(3) * t187 - t185 * t251 - t261;
t181 = t228 * pkin(1) + t262;
t236 = mrSges(3,2) * t183 - mrSges(3,3) * t181 + Ifges(3,1) * qJDD(1) - pkin(6) * t107 + t225 * t102 - t223 * t103;
t235 = -mrSges(3,1) * t181 - pkin(2) * t115 - pkin(6) * t108 - t223 * t102 - t225 * t103;
t234 = mrSges(4,1) * t169 - mrSges(4,2) * t170 + Ifges(4,5) * t203 + Ifges(4,6) * t202 + Ifges(4,3) * qJDD(3) + pkin(3) * t230 + pkin(7) * t244 + t259 * t111 + t222 * t113 + t186 * t251 + t187 * t250;
t233 = -m(3) * t181 + t228 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t115;
t232 = -mrSges(2,2) * t208 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t237) + qJ(2) * t233 + mrSges(2,1) * t207 + Ifges(2,3) * qJDD(1) + t236;
t231 = mrSges(3,1) * t183 + pkin(2) * t107 + t234;
t109 = m(2) * t208 - t228 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t233;
t106 = -m(3) * g(3) + t108;
t104 = m(2) * t207 - t228 * mrSges(2,2) + t258 * qJDD(1) + t237;
t100 = t231 + (t255 * t228) + t256 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t207 - qJ(2) * t106;
t99 = mrSges(2,3) * t208 - pkin(1) * t106 + t258 * g(3) - t255 * qJDD(1) + t256 * t228 + t235;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t226 * t100 - t224 * t99 - pkin(5) * (t226 * t104 + t224 * t109), t100, t236, t102, t113, -t198 * t149 + t238; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t224 * t100 + t226 * t99 + pkin(5) * (-t224 * t104 + t226 * t109), t99, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t228 * Ifges(3,5)) - t231, t103, t111, -t199 * t147 - t239; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t232, t232, mrSges(3,2) * g(3) + t228 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t235, t234, t261, Ifges(6,5) * t162 + Ifges(6,6) * t197 + Ifges(6,3) * t161 + t199 * t149 - t210 * t151 - t242;];
m_new = t1;
