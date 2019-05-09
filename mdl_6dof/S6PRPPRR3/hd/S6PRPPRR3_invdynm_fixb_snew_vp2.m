% Calculate vector of cutting torques with Newton-Euler for
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-04 22:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:59:59
% EndTime: 2019-05-04 22:00:12
% DurationCPUTime: 6.88s
% Computational Cost: add. (105058->255), mult. (181445->319), div. (0->0), fcn. (108453->12), ass. (0->116)
t234 = qJD(2) ^ 2;
t222 = sin(pkin(10));
t225 = cos(pkin(10));
t204 = t222 * g(1) - t225 * g(2);
t205 = -t225 * g(1) - t222 * g(2);
t217 = -g(3) + qJDD(1);
t232 = cos(qJ(2));
t226 = cos(pkin(6));
t229 = sin(qJ(2));
t254 = t226 * t229;
t223 = sin(pkin(6));
t255 = t223 * t229;
t173 = t204 * t254 + t232 * t205 + t217 * t255;
t247 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t173;
t260 = -pkin(2) - pkin(3);
t165 = t260 * t234 + t247;
t172 = -t229 * t205 + (t204 * t226 + t217 * t223) * t232;
t240 = -t234 * qJ(3) + qJDD(3) - t172;
t168 = t260 * qJDD(2) + t240;
t221 = sin(pkin(11));
t224 = cos(pkin(11));
t161 = t224 * t165 + t221 * t168;
t158 = -t234 * pkin(4) - qJDD(2) * pkin(8) + t161;
t187 = -t223 * t204 + t226 * t217;
t185 = qJDD(4) - t187;
t228 = sin(qJ(5));
t231 = cos(qJ(5));
t155 = t231 * t158 + t228 * t185;
t200 = (pkin(5) * t231 + pkin(9) * t228) * qJD(2);
t233 = qJD(5) ^ 2;
t251 = t231 * qJD(2);
t152 = -t233 * pkin(5) + qJDD(5) * pkin(9) - t200 * t251 + t155;
t160 = -t221 * t165 + t224 * t168;
t157 = qJDD(2) * pkin(4) - t234 * pkin(8) - t160;
t250 = qJD(2) * qJD(5);
t248 = t231 * t250;
t201 = -t228 * qJDD(2) - t248;
t249 = t228 * t250;
t202 = -t231 * qJDD(2) + t249;
t153 = (-t201 + t248) * pkin(9) + (-t202 - t249) * pkin(5) + t157;
t227 = sin(qJ(6));
t230 = cos(qJ(6));
t148 = -t227 * t152 + t230 * t153;
t252 = qJD(2) * t228;
t197 = t230 * qJD(5) + t227 * t252;
t180 = t197 * qJD(6) + t227 * qJDD(5) + t230 * t201;
t198 = t227 * qJD(5) - t230 * t252;
t181 = -t197 * mrSges(7,1) + t198 * mrSges(7,2);
t209 = qJD(6) + t251;
t183 = -t209 * mrSges(7,2) + t197 * mrSges(7,3);
t194 = qJDD(6) - t202;
t146 = m(7) * t148 + t194 * mrSges(7,1) - t180 * mrSges(7,3) - t198 * t181 + t209 * t183;
t149 = t230 * t152 + t227 * t153;
t179 = -t198 * qJD(6) + t230 * qJDD(5) - t227 * t201;
t184 = t209 * mrSges(7,1) - t198 * mrSges(7,3);
t147 = m(7) * t149 - t194 * mrSges(7,2) + t179 * mrSges(7,3) + t197 * t181 - t209 * t184;
t141 = -t227 * t146 + t230 * t147;
t199 = (mrSges(6,1) * t231 - mrSges(6,2) * t228) * qJD(2);
t206 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t252;
t138 = m(6) * t155 - qJDD(5) * mrSges(6,2) + t202 * mrSges(6,3) - qJD(5) * t206 - t199 * t251 + t141;
t253 = t231 * t185;
t151 = -qJDD(5) * pkin(5) - t233 * pkin(9) - t253 + (-qJD(2) * t200 + t158) * t228;
t150 = -m(7) * t151 + t179 * mrSges(7,1) - t180 * mrSges(7,2) + t197 * t183 - t198 * t184;
t154 = -t228 * t158 + t253;
t207 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t251;
t144 = m(6) * t154 + qJDD(5) * mrSges(6,1) - t201 * mrSges(6,3) + qJD(5) * t207 + t199 * t252 + t150;
t134 = t228 * t138 + t231 * t144;
t133 = -m(5) * t185 - t134;
t132 = m(4) * t187 + t133;
t140 = t230 * t146 + t227 * t147;
t174 = Ifges(7,5) * t198 + Ifges(7,6) * t197 + Ifges(7,3) * t209;
t176 = Ifges(7,1) * t198 + Ifges(7,4) * t197 + Ifges(7,5) * t209;
t142 = -mrSges(7,1) * t151 + mrSges(7,3) * t149 + Ifges(7,4) * t180 + Ifges(7,2) * t179 + Ifges(7,6) * t194 - t198 * t174 + t209 * t176;
t175 = Ifges(7,4) * t198 + Ifges(7,2) * t197 + Ifges(7,6) * t209;
t143 = mrSges(7,2) * t151 - mrSges(7,3) * t148 + Ifges(7,1) * t180 + Ifges(7,4) * t179 + Ifges(7,5) * t194 + t197 * t174 - t209 * t175;
t189 = (Ifges(6,3) * qJD(5)) + (-Ifges(6,5) * t228 - Ifges(6,6) * t231) * qJD(2);
t190 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t228 - Ifges(6,2) * t231) * qJD(2);
t127 = mrSges(6,2) * t157 - mrSges(6,3) * t154 + Ifges(6,1) * t201 + Ifges(6,4) * t202 + Ifges(6,5) * qJDD(5) - pkin(9) * t140 - qJD(5) * t190 - t227 * t142 + t230 * t143 - t189 * t251;
t191 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t228 - Ifges(6,4) * t231) * qJD(2);
t237 = mrSges(7,1) * t148 - mrSges(7,2) * t149 + Ifges(7,5) * t180 + Ifges(7,6) * t179 + Ifges(7,3) * t194 + t198 * t175 - t197 * t176;
t128 = -mrSges(6,1) * t157 + mrSges(6,3) * t155 + Ifges(6,4) * t201 + Ifges(6,2) * t202 + Ifges(6,6) * qJDD(5) - pkin(5) * t140 + qJD(5) * t191 + t189 * t252 - t237;
t115 = mrSges(5,2) * t185 - mrSges(5,3) * t160 - Ifges(5,5) * qJDD(2) - t234 * Ifges(5,6) - pkin(8) * t134 + t231 * t127 - t228 * t128;
t261 = mrSges(6,1) * t154 - mrSges(6,2) * t155 + Ifges(6,5) * t201 + Ifges(6,6) * t202 + Ifges(6,3) * qJDD(5) + pkin(5) * t150 + pkin(9) * t141 - (t228 * t190 - t231 * t191) * qJD(2) + t230 * t142 + t227 * t143;
t119 = -mrSges(5,1) * t185 + mrSges(5,3) * t161 + t234 * Ifges(5,5) - Ifges(5,6) * qJDD(2) - pkin(4) * t134 - t261;
t135 = t231 * t138 - t228 * t144;
t130 = m(5) * t161 - t234 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + t135;
t139 = -m(6) * t157 + t202 * mrSges(6,1) - t201 * mrSges(6,2) + t206 * t252 - t207 * t251 - t140;
t136 = m(5) * t160 - qJDD(2) * mrSges(5,1) - t234 * mrSges(5,2) + t139;
t126 = t224 * t130 - t221 * t136;
t169 = -t234 * pkin(2) + t247;
t239 = mrSges(4,2) * t169 - pkin(3) * t133 - qJ(4) * t126 - t221 * t115 - t224 * t119;
t258 = -mrSges(3,1) - mrSges(4,1);
t105 = mrSges(3,3) * t173 - pkin(2) * t132 + (Ifges(4,4) + Ifges(3,5)) * t234 + t258 * t187 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t239;
t244 = m(4) * t169 + qJDD(2) * mrSges(4,3) + t126;
t123 = m(3) * t173 - qJDD(2) * mrSges(3,2) + t258 * t234 + t244;
t125 = t221 * t130 + t224 * t136;
t171 = -qJDD(2) * pkin(2) + t240;
t241 = -m(4) * t171 + qJDD(2) * mrSges(4,1) + t234 * mrSges(4,3) - t125;
t124 = m(3) * t172 + qJDD(2) * mrSges(3,1) - t234 * mrSges(3,2) + t241;
t118 = t232 * t123 - t229 * t124;
t262 = pkin(7) * t118 + t105 * t232;
t256 = t124 * t232;
t131 = m(3) * t187 + t132;
t113 = t123 * t254 - t223 * t131 + t226 * t256;
t238 = mrSges(5,1) * t160 - mrSges(5,2) * t161 - Ifges(5,3) * qJDD(2) + pkin(4) * t139 + pkin(8) * t135 + t228 * t127 + t231 * t128;
t236 = -mrSges(4,1) * t171 + mrSges(4,3) * t169 + Ifges(4,2) * qJDD(2) - pkin(3) * t125 - t238;
t107 = t236 + qJ(3) * (-t234 * mrSges(4,1) + t244) + pkin(2) * t241 + mrSges(3,1) * t172 - mrSges(3,2) * t173 + Ifges(3,3) * qJDD(2);
t242 = mrSges(4,2) * t171 + Ifges(4,4) * qJDD(2) + t234 * Ifges(4,6) - qJ(4) * t125 + t224 * t115 - t221 * t119;
t109 = -mrSges(3,3) * t172 + Ifges(3,5) * qJDD(2) - t234 * Ifges(3,6) - qJ(3) * t132 + (mrSges(3,2) - mrSges(4,3)) * t187 + t242;
t243 = mrSges(2,1) * t204 - mrSges(2,2) * t205 + pkin(1) * t113 + t226 * t107 + t109 * t255 + t262 * t223;
t116 = m(2) * t205 + t118;
t112 = t226 * t131 + (t123 * t229 + t256) * t223;
t110 = m(2) * t204 + t113;
t103 = mrSges(2,2) * t217 - mrSges(2,3) * t204 - t229 * t105 + t232 * t109 + (-t112 * t223 - t113 * t226) * pkin(7);
t102 = -mrSges(2,1) * t217 + mrSges(2,3) * t205 - pkin(1) * t112 - t223 * t107 + (t109 * t229 + t262) * t226;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t225 * t103 - t222 * t102 - qJ(1) * (t225 * t110 + t222 * t116), t103, t109, -mrSges(4,3) * t187 + t242, t115, t127, t143; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t222 * t103 + t225 * t102 + qJ(1) * (-t222 * t110 + t225 * t116), t102, t105, t236, t119, t128, t142; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t243, t243, t107, mrSges(4,1) * t187 - t234 * Ifges(4,4) + Ifges(4,6) * qJDD(2) - t239, t238, t261, t237;];
m_new  = t1;
