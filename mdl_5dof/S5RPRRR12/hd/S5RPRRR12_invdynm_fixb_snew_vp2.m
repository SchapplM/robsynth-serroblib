% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR12
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:34
% EndTime: 2019-12-31 19:12:40
% DurationCPUTime: 3.91s
% Computational Cost: add. (49819->268), mult. (97647->330), div. (0->0), fcn. (60841->8), ass. (0->109)
t232 = sin(qJ(1));
t236 = cos(qJ(1));
t213 = -t236 * g(1) - t232 * g(2);
t250 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t213;
t230 = sin(qJ(4));
t231 = sin(qJ(3));
t234 = cos(qJ(4));
t235 = cos(qJ(3));
t197 = (t230 * t235 + t231 * t234) * qJD(1);
t263 = -pkin(1) - pkin(6);
t262 = mrSges(2,1) - mrSges(3,2);
t261 = Ifges(2,5) - Ifges(3,4);
t260 = (-Ifges(2,6) + Ifges(3,5));
t212 = t232 * g(1) - t236 * g(2);
t237 = qJD(1) ^ 2;
t249 = -t237 * qJ(2) + qJDD(2) - t212;
t187 = t263 * qJDD(1) + t249;
t177 = t231 * g(3) + t235 * t187;
t257 = qJD(1) * qJD(3);
t255 = t231 * t257;
t207 = t235 * qJDD(1) - t255;
t156 = (-t207 - t255) * pkin(7) + (-t231 * t235 * t237 + qJDD(3)) * pkin(3) + t177;
t178 = -t235 * g(3) + t231 * t187;
t206 = -t231 * qJDD(1) - t235 * t257;
t258 = qJD(1) * t235;
t211 = qJD(3) * pkin(3) - pkin(7) * t258;
t226 = t231 ^ 2;
t157 = -t226 * t237 * pkin(3) + t206 * pkin(7) - qJD(3) * t211 + t178;
t146 = t230 * t156 + t234 * t157;
t198 = (-t230 * t231 + t234 * t235) * qJD(1);
t166 = -t198 * qJD(4) + t234 * t206 - t230 * t207;
t174 = t197 * mrSges(5,1) + t198 * mrSges(5,2);
t220 = qJD(3) + qJD(4);
t185 = t220 * mrSges(5,1) - t198 * mrSges(5,3);
t219 = qJDD(3) + qJDD(4);
t161 = -t206 * pkin(3) + t211 * t258 + (-pkin(7) * t226 + t263) * t237 + t250;
t167 = -t197 * qJD(4) + t230 * t206 + t234 * t207;
t141 = (t197 * t220 - t167) * pkin(8) + (t198 * t220 - t166) * pkin(4) + t161;
t175 = t197 * pkin(4) - t198 * pkin(8);
t218 = t220 ^ 2;
t143 = -t218 * pkin(4) + t219 * pkin(8) - t197 * t175 + t146;
t229 = sin(qJ(5));
t233 = cos(qJ(5));
t139 = t233 * t141 - t229 * t143;
t179 = -t229 * t198 + t233 * t220;
t149 = t179 * qJD(5) + t233 * t167 + t229 * t219;
t180 = t233 * t198 + t229 * t220;
t158 = -t179 * mrSges(6,1) + t180 * mrSges(6,2);
t165 = qJDD(5) - t166;
t193 = qJD(5) + t197;
t168 = -t193 * mrSges(6,2) + t179 * mrSges(6,3);
t135 = m(6) * t139 + t165 * mrSges(6,1) - t149 * mrSges(6,3) - t180 * t158 + t193 * t168;
t140 = t229 * t141 + t233 * t143;
t148 = -t180 * qJD(5) - t229 * t167 + t233 * t219;
t169 = t193 * mrSges(6,1) - t180 * mrSges(6,3);
t136 = m(6) * t140 - t165 * mrSges(6,2) + t148 * mrSges(6,3) + t179 * t158 - t193 * t169;
t253 = -t229 * t135 + t233 * t136;
t122 = m(5) * t146 - t219 * mrSges(5,2) + t166 * mrSges(5,3) - t197 * t174 - t220 * t185 + t253;
t145 = t234 * t156 - t230 * t157;
t184 = -t220 * mrSges(5,2) - t197 * mrSges(5,3);
t142 = -t219 * pkin(4) - t218 * pkin(8) + t198 * t175 - t145;
t247 = -m(6) * t142 + t148 * mrSges(6,1) - t149 * mrSges(6,2) + t179 * t168 - t180 * t169;
t131 = m(5) * t145 + t219 * mrSges(5,1) - t167 * mrSges(5,3) - t198 * t174 + t220 * t184 + t247;
t116 = t230 * t122 + t234 * t131;
t124 = t233 * t135 + t229 * t136;
t259 = qJD(1) * t231;
t205 = (mrSges(4,1) * t231 + mrSges(4,2) * t235) * qJD(1);
t209 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t259;
t113 = m(4) * t177 + qJDD(3) * mrSges(4,1) - t207 * mrSges(4,3) + qJD(3) * t209 - t205 * t258 + t116;
t210 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t258;
t254 = t234 * t122 - t230 * t131;
t114 = m(4) * t178 - qJDD(3) * mrSges(4,2) + t206 * mrSges(4,3) - qJD(3) * t210 - t205 * t259 + t254;
t109 = -t231 * t113 + t235 * t114;
t108 = t235 * t113 + t231 * t114;
t192 = -qJDD(1) * pkin(1) + t249;
t248 = -m(3) * t192 + (t237 * mrSges(3,3)) - t108;
t246 = m(5) * t161 - t166 * mrSges(5,1) + t167 * mrSges(5,2) + t197 * t184 + t198 * t185 + t124;
t150 = Ifges(6,5) * t180 + Ifges(6,6) * t179 + Ifges(6,3) * t193;
t152 = Ifges(6,1) * t180 + Ifges(6,4) * t179 + Ifges(6,5) * t193;
t128 = -mrSges(6,1) * t142 + mrSges(6,3) * t140 + Ifges(6,4) * t149 + Ifges(6,2) * t148 + Ifges(6,6) * t165 - t180 * t150 + t193 * t152;
t151 = Ifges(6,4) * t180 + Ifges(6,2) * t179 + Ifges(6,6) * t193;
t129 = mrSges(6,2) * t142 - mrSges(6,3) * t139 + Ifges(6,1) * t149 + Ifges(6,4) * t148 + Ifges(6,5) * t165 + t179 * t150 - t193 * t151;
t170 = Ifges(5,5) * t198 - Ifges(5,6) * t197 + Ifges(5,3) * t220;
t171 = Ifges(5,4) * t198 - Ifges(5,2) * t197 + Ifges(5,6) * t220;
t110 = mrSges(5,2) * t161 - mrSges(5,3) * t145 + Ifges(5,1) * t167 + Ifges(5,4) * t166 + Ifges(5,5) * t219 - pkin(8) * t124 - t229 * t128 + t233 * t129 - t197 * t170 - t220 * t171;
t172 = Ifges(5,1) * t198 - Ifges(5,4) * t197 + Ifges(5,5) * t220;
t241 = mrSges(6,1) * t139 - mrSges(6,2) * t140 + Ifges(6,5) * t149 + Ifges(6,6) * t148 + Ifges(6,3) * t165 + t180 * t151 - t179 * t152;
t111 = -mrSges(5,1) * t161 + mrSges(5,3) * t146 + Ifges(5,4) * t167 + Ifges(5,2) * t166 + Ifges(5,6) * t219 - pkin(4) * t124 - t198 * t170 + t220 * t172 - t241;
t186 = t263 * t237 + t250;
t194 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t235 - Ifges(4,6) * t231) * qJD(1);
t196 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t235 - Ifges(4,4) * t231) * qJD(1);
t102 = -mrSges(4,1) * t186 + mrSges(4,3) * t178 + Ifges(4,4) * t207 + Ifges(4,2) * t206 + Ifges(4,6) * qJDD(3) - pkin(3) * t246 + pkin(7) * t254 + qJD(3) * t196 + t230 * t110 + t234 * t111 - t194 * t258;
t195 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t235 - Ifges(4,2) * t231) * qJD(1);
t104 = mrSges(4,2) * t186 - mrSges(4,3) * t177 + Ifges(4,1) * t207 + Ifges(4,4) * t206 + Ifges(4,5) * qJDD(3) - pkin(7) * t116 - qJD(3) * t195 + t234 * t110 - t230 * t111 - t194 * t259;
t190 = t237 * pkin(1) - t250;
t245 = mrSges(3,2) * t192 - mrSges(3,3) * t190 + Ifges(3,1) * qJDD(1) - pkin(6) * t108 - t231 * t102 + t235 * t104;
t119 = -m(4) * t186 + t206 * mrSges(4,1) - t207 * mrSges(4,2) - t209 * t259 - t210 * t258 - t246;
t244 = -mrSges(3,1) * t190 - pkin(2) * t119 - pkin(6) * t109 - t235 * t102 - t231 * t104;
t243 = -mrSges(5,1) * t145 + mrSges(5,2) * t146 - Ifges(5,5) * t167 - Ifges(5,6) * t166 - Ifges(5,3) * t219 - pkin(4) * t247 - pkin(8) * t253 - t233 * t128 - t229 * t129 - t198 * t171 - t197 * t172;
t240 = -m(3) * t190 + t237 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t119;
t242 = -mrSges(2,2) * t213 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t248) + qJ(2) * t240 + mrSges(2,1) * t212 + Ifges(2,3) * qJDD(1) + t245;
t239 = -mrSges(4,1) * t177 + mrSges(4,2) * t178 - Ifges(4,5) * t207 - Ifges(4,6) * t206 - Ifges(4,3) * qJDD(3) - pkin(3) * t116 - t195 * t258 - t196 * t259 + t243;
t238 = -mrSges(3,1) * t192 - pkin(2) * t108 + t239;
t117 = m(2) * t213 - t237 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t240;
t107 = -m(3) * g(3) + t109;
t105 = m(2) * t212 - t237 * mrSges(2,2) + t262 * qJDD(1) + t248;
t101 = -t238 + (t260 * t237) + t261 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t212 - qJ(2) * t107;
t100 = mrSges(2,3) * t213 - pkin(1) * t107 + t262 * g(3) - t260 * qJDD(1) + t261 * t237 + t244;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t236 * t101 - t232 * t100 - pkin(5) * (t236 * t105 + t232 * t117), t101, t245, t104, t110, t129; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t232 * t101 + t236 * t100 + pkin(5) * (-t232 * t105 + t236 * t117), t100, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t237 * Ifges(3,5)) + t238, t102, t111, t128; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t242, t242, mrSges(3,2) * g(3) + t237 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t244, -t239, -t243, t241;];
m_new = t1;
