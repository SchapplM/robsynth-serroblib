% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP12
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:19
% EndTime: 2019-12-31 18:56:24
% DurationCPUTime: 1.94s
% Computational Cost: add. (21824->261), mult. (41305->307), div. (0->0), fcn. (22524->6), ass. (0->100)
t226 = sin(qJ(1));
t229 = cos(qJ(1));
t213 = -t229 * g(1) - t226 * g(2);
t264 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t213;
t224 = sin(qJ(4));
t227 = cos(qJ(4));
t228 = cos(qJ(3));
t255 = qJD(1) * t228;
t203 = qJD(3) * t227 - t224 * t255;
t225 = sin(qJ(3));
t253 = qJD(1) * qJD(3);
t249 = t225 * t253;
t208 = qJDD(1) * t228 - t249;
t167 = qJD(4) * t203 + qJDD(3) * t224 + t208 * t227;
t204 = qJD(3) * t224 + t227 * t255;
t170 = -t203 * mrSges(6,1) + t204 * mrSges(6,2);
t231 = qJD(1) ^ 2;
t262 = (-pkin(1) - pkin(6));
t181 = (t231 * t262) - t264;
t248 = t228 * t253;
t207 = -qJDD(1) * t225 - t248;
t141 = (-t208 + t249) * pkin(7) + (-t207 + t248) * pkin(3) + t181;
t212 = t226 * g(1) - t229 * g(2);
t244 = -t231 * qJ(2) + qJDD(2) - t212;
t182 = qJDD(1) * t262 + t244;
t173 = -g(3) * t228 + t225 * t182;
t206 = (pkin(3) * t225 - pkin(7) * t228) * qJD(1);
t230 = qJD(3) ^ 2;
t254 = t225 * qJD(1);
t146 = -pkin(3) * t230 + qJDD(3) * pkin(7) - t206 * t254 + t173;
t137 = t141 * t227 - t224 * t146;
t202 = qJDD(4) - t207;
t214 = qJD(4) + t254;
t131 = -0.2e1 * qJD(5) * t204 + (t203 * t214 - t167) * qJ(5) + (t203 * t204 + t202) * pkin(4) + t137;
t174 = -mrSges(6,2) * t214 + mrSges(6,3) * t203;
t251 = m(6) * t131 + t202 * mrSges(6,1) + t214 * t174;
t128 = -t167 * mrSges(6,3) - t204 * t170 + t251;
t138 = t141 * t224 + t146 * t227;
t152 = Ifges(5,4) * t204 + Ifges(5,2) * t203 + Ifges(5,6) * t214;
t153 = Ifges(6,1) * t204 + Ifges(6,4) * t203 + Ifges(6,5) * t214;
t154 = Ifges(5,1) * t204 + Ifges(5,4) * t203 + Ifges(5,5) * t214;
t166 = -qJD(4) * t204 + qJDD(3) * t227 - t208 * t224;
t176 = pkin(4) * t214 - qJ(5) * t204;
t201 = t203 ^ 2;
t134 = -pkin(4) * t201 + t166 * qJ(5) + 0.2e1 * qJD(5) * t203 - t176 * t214 + t138;
t151 = Ifges(6,4) * t204 + Ifges(6,2) * t203 + Ifges(6,6) * t214;
t242 = -mrSges(6,1) * t131 + mrSges(6,2) * t134 - Ifges(6,5) * t167 - Ifges(6,6) * t166 - Ifges(6,3) * t202 - t151 * t204;
t263 = mrSges(5,1) * t137 - mrSges(5,2) * t138 + Ifges(5,5) * t167 + Ifges(5,6) * t166 + Ifges(5,3) * t202 + pkin(4) * t128 + t204 * t152 - (t154 + t153) * t203 - t242;
t261 = mrSges(2,1) - mrSges(3,2);
t260 = -mrSges(5,2) - mrSges(6,2);
t259 = -Ifges(3,4) + Ifges(2,5);
t258 = (Ifges(3,5) - Ifges(2,6));
t171 = -mrSges(5,1) * t203 + mrSges(5,2) * t204;
t175 = -mrSges(5,2) * t214 + mrSges(5,3) * t203;
t121 = m(5) * t137 + t202 * mrSges(5,1) + t214 * t175 + (-t170 - t171) * t204 + (-mrSges(5,3) - mrSges(6,3)) * t167 + t251;
t250 = m(6) * t134 + t166 * mrSges(6,3) + t203 * t170;
t177 = mrSges(6,1) * t214 - mrSges(6,3) * t204;
t256 = -mrSges(5,1) * t214 + mrSges(5,3) * t204 - t177;
t124 = m(5) * t138 + t166 * mrSges(5,3) + t203 * t171 + t202 * t260 + t214 * t256 + t250;
t118 = t121 * t227 + t124 * t224;
t205 = (mrSges(4,1) * t225 + mrSges(4,2) * t228) * qJD(1);
t211 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t255;
t247 = -t121 * t224 + t124 * t227;
t116 = m(4) * t173 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t207 - qJD(3) * t211 - t205 * t254 + t247;
t172 = g(3) * t225 + t182 * t228;
t210 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t254;
t145 = -qJDD(3) * pkin(3) - pkin(7) * t230 + t206 * t255 - t172;
t136 = -t166 * pkin(4) - qJ(5) * t201 + t176 * t204 + qJDD(5) + t145;
t246 = -m(6) * t136 + t166 * mrSges(6,1) + t203 * t174;
t235 = -m(5) * t145 + t166 * mrSges(5,1) + t167 * t260 + t203 * t175 + t204 * t256 + t246;
t125 = m(4) * t172 + qJDD(3) * mrSges(4,1) - t208 * mrSges(4,3) + qJD(3) * t210 - t205 * t255 + t235;
t107 = t116 * t228 - t125 * t225;
t106 = t225 * t116 + t228 * t125;
t243 = -mrSges(6,1) * t136 + mrSges(6,3) * t134 + Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * t202 + t214 * t153;
t149 = Ifges(6,5) * t204 + Ifges(6,6) * t203 + Ifges(6,3) * t214;
t241 = mrSges(6,2) * t136 - mrSges(6,3) * t131 + Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * t202 + t149 * t203;
t187 = -qJDD(1) * pkin(1) + t244;
t240 = -m(3) * t187 + (t231 * mrSges(3,3)) - t106;
t114 = -m(4) * t181 + mrSges(4,1) * t207 - t208 * mrSges(4,2) - t210 * t254 - t211 * t255 - t118;
t150 = Ifges(5,5) * t204 + Ifges(5,6) * t203 + Ifges(5,3) * t214;
t109 = Ifges(5,4) * t167 + Ifges(5,2) * t166 + Ifges(5,6) * t202 + t214 * t154 - mrSges(5,1) * t145 + mrSges(5,3) * t138 - pkin(4) * (t167 * mrSges(6,2) - t246) + qJ(5) * (-t202 * mrSges(6,2) - t214 * t177 + t250) + (-pkin(4) * t177 - t149 - t150) * t204 + t243;
t113 = mrSges(5,2) * t145 - mrSges(5,3) * t137 + Ifges(5,1) * t167 + Ifges(5,4) * t166 + Ifges(5,5) * t202 - qJ(5) * t128 + t203 * t150 + (-t151 - t152) * t214 + t241;
t189 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t228 - Ifges(4,6) * t225) * qJD(1);
t190 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t228 - Ifges(4,2) * t225) * qJD(1);
t101 = mrSges(4,2) * t181 - mrSges(4,3) * t172 + Ifges(4,1) * t208 + Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3) - pkin(7) * t118 - qJD(3) * t190 - t109 * t224 + t113 * t227 - t189 * t254;
t191 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t228 - Ifges(4,4) * t225) * qJD(1);
t102 = -mrSges(4,1) * t181 + mrSges(4,3) * t173 + Ifges(4,4) * t208 + Ifges(4,2) * t207 + Ifges(4,6) * qJDD(3) - pkin(3) * t118 + qJD(3) * t191 - t189 * t255 - t263;
t185 = t231 * pkin(1) + t264;
t239 = mrSges(3,2) * t187 - mrSges(3,3) * t185 + Ifges(3,1) * qJDD(1) - pkin(6) * t106 + t101 * t228 - t102 * t225;
t238 = -mrSges(3,1) * t185 - pkin(2) * t114 - pkin(6) * t107 - t225 * t101 - t228 * t102;
t237 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,5) * t208 + Ifges(4,6) * t207 + Ifges(4,3) * qJDD(3) + pkin(3) * t235 + pkin(7) * t247 + t109 * t227 + t113 * t224 + t190 * t255 + t191 * t254;
t236 = -m(3) * t185 + t231 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t114;
t234 = -mrSges(2,2) * t213 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t240) + qJ(2) * t236 + mrSges(2,1) * t212 + Ifges(2,3) * qJDD(1) + t239;
t233 = mrSges(3,1) * t187 + pkin(2) * t106 + t237;
t110 = m(2) * t213 - mrSges(2,1) * t231 - qJDD(1) * mrSges(2,2) + t236;
t105 = -m(3) * g(3) + t107;
t103 = m(2) * t212 - t231 * mrSges(2,2) + qJDD(1) * t261 + t240;
t99 = -mrSges(2,3) * t212 + (t258 * t231) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t259 * qJDD(1) + t233 - qJ(2) * t105;
t98 = mrSges(2,3) * t213 - pkin(1) * t105 + g(3) * t261 - qJDD(1) * t258 + t231 * t259 + t238;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t229 * t99 - t226 * t98 - pkin(5) * (t103 * t229 + t110 * t226), t99, t239, t101, t113, -t151 * t214 + t241; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t226 * t99 + t229 * t98 + pkin(5) * (-t103 * t226 + t110 * t229), t98, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t231 * Ifges(3,5)) - t233, t102, t109, -t204 * t149 + t243; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t234, t234, mrSges(3,2) * g(3) + t231 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t238, t237, t263, -t203 * t153 - t242;];
m_new = t1;
