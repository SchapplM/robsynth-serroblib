% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:25
% EndTime: 2019-12-05 16:30:42
% DurationCPUTime: 7.58s
% Computational Cost: add. (122491->265), mult. (250168->346), div. (0->0), fcn. (171715->12), ass. (0->115)
t228 = sin(pkin(9));
t231 = cos(pkin(9));
t219 = t228 * g(1) - t231 * g(2);
t220 = -t231 * g(1) - t228 * g(2);
t226 = -g(3) + qJDD(1);
t238 = cos(qJ(2));
t232 = cos(pkin(5));
t235 = sin(qJ(2));
t255 = t232 * t235;
t229 = sin(pkin(5));
t256 = t229 * t235;
t182 = t219 * t255 + t238 * t220 + t226 * t256;
t240 = qJD(2) ^ 2;
t178 = -t240 * pkin(2) + qJDD(2) * pkin(7) + t182;
t197 = -t229 * t219 + t232 * t226;
t234 = sin(qJ(3));
t237 = cos(qJ(3));
t173 = t237 * t178 + t234 * t197;
t214 = (-pkin(3) * t237 - qJ(4) * t234) * qJD(2);
t239 = qJD(3) ^ 2;
t253 = t237 * qJD(2);
t159 = -t239 * pkin(3) + qJDD(3) * qJ(4) + t214 * t253 + t173;
t181 = -t235 * t220 + (t219 * t232 + t226 * t229) * t238;
t177 = -qJDD(2) * pkin(2) - t240 * pkin(7) - t181;
t252 = qJD(2) * qJD(3);
t251 = t237 * t252;
t216 = t234 * qJDD(2) + t251;
t225 = t234 * t252;
t217 = t237 * qJDD(2) - t225;
t165 = (-t216 - t251) * qJ(4) + (-t217 + t225) * pkin(3) + t177;
t227 = sin(pkin(10));
t230 = cos(pkin(10));
t254 = qJD(2) * t234;
t209 = t227 * qJD(3) + t230 * t254;
t153 = -0.2e1 * qJD(4) * t209 - t227 * t159 + t230 * t165;
t195 = t227 * qJDD(3) + t230 * t216;
t208 = t230 * qJD(3) - t227 * t254;
t151 = (-t208 * t253 - t195) * pkin(8) + (t208 * t209 - t217) * pkin(4) + t153;
t154 = 0.2e1 * qJD(4) * t208 + t230 * t159 + t227 * t165;
t194 = t230 * qJDD(3) - t227 * t216;
t196 = -pkin(4) * t253 - t209 * pkin(8);
t207 = t208 ^ 2;
t152 = -t207 * pkin(4) + t194 * pkin(8) + t196 * t253 + t154;
t233 = sin(qJ(5));
t236 = cos(qJ(5));
t149 = t236 * t151 - t233 * t152;
t187 = t236 * t208 - t233 * t209;
t167 = t187 * qJD(5) + t233 * t194 + t236 * t195;
t188 = t233 * t208 + t236 * t209;
t174 = -t187 * mrSges(6,1) + t188 * mrSges(6,2);
t224 = qJD(5) - t253;
t179 = -t224 * mrSges(6,2) + t187 * mrSges(6,3);
t211 = qJDD(5) - t217;
t145 = m(6) * t149 + t211 * mrSges(6,1) - t167 * mrSges(6,3) - t188 * t174 + t224 * t179;
t150 = t233 * t151 + t236 * t152;
t166 = -t188 * qJD(5) + t236 * t194 - t233 * t195;
t180 = t224 * mrSges(6,1) - t188 * mrSges(6,3);
t146 = m(6) * t150 - t211 * mrSges(6,2) + t166 * mrSges(6,3) + t187 * t174 - t224 * t180;
t137 = t236 * t145 + t233 * t146;
t189 = -t208 * mrSges(5,1) + t209 * mrSges(5,2);
t192 = mrSges(5,2) * t253 + t208 * mrSges(5,3);
t135 = m(5) * t153 - t217 * mrSges(5,1) - t195 * mrSges(5,3) - t209 * t189 - t192 * t253 + t137;
t193 = -mrSges(5,1) * t253 - t209 * mrSges(5,3);
t249 = -t233 * t145 + t236 * t146;
t136 = m(5) * t154 + t217 * mrSges(5,2) + t194 * mrSges(5,3) + t208 * t189 + t193 * t253 + t249;
t133 = -t227 * t135 + t230 * t136;
t215 = (-mrSges(4,1) * t237 + mrSges(4,2) * t234) * qJD(2);
t221 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t254;
t131 = m(4) * t173 - qJDD(3) * mrSges(4,2) + t217 * mrSges(4,3) - qJD(3) * t221 + t215 * t253 + t133;
t172 = -t234 * t178 + t237 * t197;
t158 = -qJDD(3) * pkin(3) - t239 * qJ(4) + t214 * t254 + qJDD(4) - t172;
t155 = -t194 * pkin(4) - t207 * pkin(8) + t209 * t196 + t158;
t245 = m(6) * t155 - t166 * mrSges(6,1) + t167 * mrSges(6,2) - t187 * t179 + t188 * t180;
t147 = -m(5) * t158 + t194 * mrSges(5,1) - t195 * mrSges(5,2) + t208 * t192 - t209 * t193 - t245;
t222 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t253;
t141 = m(4) * t172 + qJDD(3) * mrSges(4,1) - t216 * mrSges(4,3) + qJD(3) * t222 - t215 * t254 + t147;
t124 = t234 * t131 + t237 * t141;
t168 = Ifges(6,5) * t188 + Ifges(6,6) * t187 + Ifges(6,3) * t224;
t170 = Ifges(6,1) * t188 + Ifges(6,4) * t187 + Ifges(6,5) * t224;
t138 = -mrSges(6,1) * t155 + mrSges(6,3) * t150 + Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * t211 - t188 * t168 + t224 * t170;
t169 = Ifges(6,4) * t188 + Ifges(6,2) * t187 + Ifges(6,6) * t224;
t139 = mrSges(6,2) * t155 - mrSges(6,3) * t149 + Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * t211 + t187 * t168 - t224 * t169;
t183 = Ifges(5,5) * t209 + Ifges(5,6) * t208 - Ifges(5,3) * t253;
t185 = Ifges(5,1) * t209 + Ifges(5,4) * t208 - Ifges(5,5) * t253;
t125 = -mrSges(5,1) * t158 + mrSges(5,3) * t154 + Ifges(5,4) * t195 + Ifges(5,2) * t194 - Ifges(5,6) * t217 - pkin(4) * t245 + pkin(8) * t249 + t236 * t138 + t233 * t139 - t209 * t183 - t185 * t253;
t184 = Ifges(5,4) * t209 + Ifges(5,2) * t208 - Ifges(5,6) * t253;
t126 = mrSges(5,2) * t158 - mrSges(5,3) * t153 + Ifges(5,1) * t195 + Ifges(5,4) * t194 - Ifges(5,5) * t217 - pkin(8) * t137 - t233 * t138 + t236 * t139 + t208 * t183 + t184 * t253;
t203 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t234 + Ifges(4,2) * t237) * qJD(2);
t204 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t234 + Ifges(4,4) * t237) * qJD(2);
t260 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,5) * t216 + Ifges(4,6) * t217 + Ifges(4,3) * qJDD(3) + pkin(3) * t147 + qJ(4) * t133 + t230 * t125 + t227 * t126 + (t234 * t203 - t237 * t204) * qJD(2);
t110 = -mrSges(3,1) * t197 + mrSges(3,3) * t182 + t240 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t124 - t260;
t250 = t237 * t131 - t234 * t141;
t122 = m(3) * t182 - t240 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t250;
t132 = t230 * t135 + t227 * t136;
t243 = -m(4) * t177 + t217 * mrSges(4,1) - t216 * mrSges(4,2) - t221 * t254 + t222 * t253 - t132;
t128 = m(3) * t181 + qJDD(2) * mrSges(3,1) - t240 * mrSges(3,2) + t243;
t118 = t238 * t122 - t235 * t128;
t261 = pkin(6) * t118 + t110 * t238;
t257 = t128 * t238;
t123 = m(3) * t197 + t124;
t114 = t122 * t255 - t229 * t123 + t232 * t257;
t202 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t234 + Ifges(4,6) * t237) * qJD(2);
t115 = mrSges(4,2) * t177 - mrSges(4,3) * t172 + Ifges(4,1) * t216 + Ifges(4,4) * t217 + Ifges(4,5) * qJDD(3) - qJ(4) * t132 - qJD(3) * t203 - t227 * t125 + t230 * t126 + t202 * t253;
t244 = -mrSges(6,1) * t149 + mrSges(6,2) * t150 - Ifges(6,5) * t167 - Ifges(6,6) * t166 - Ifges(6,3) * t211 - t188 * t169 + t187 * t170;
t241 = -mrSges(5,1) * t153 + mrSges(5,2) * t154 - Ifges(5,5) * t195 - Ifges(5,6) * t194 - pkin(4) * t137 - t209 * t184 + t208 * t185 + t244;
t119 = Ifges(4,6) * qJDD(3) + t241 + (Ifges(4,2) + Ifges(5,3)) * t217 + Ifges(4,4) * t216 + qJD(3) * t204 - mrSges(4,1) * t177 + mrSges(4,3) * t173 - pkin(3) * t132 - t202 * t254;
t106 = mrSges(3,1) * t181 - mrSges(3,2) * t182 + Ifges(3,3) * qJDD(2) + pkin(2) * t243 + pkin(7) * t250 + t234 * t115 + t237 * t119;
t108 = mrSges(3,2) * t197 - mrSges(3,3) * t181 + Ifges(3,5) * qJDD(2) - t240 * Ifges(3,6) - pkin(7) * t124 + t237 * t115 - t234 * t119;
t246 = mrSges(2,1) * t219 - mrSges(2,2) * t220 + pkin(1) * t114 + t232 * t106 + t108 * t256 + t261 * t229;
t116 = m(2) * t220 + t118;
t113 = t232 * t123 + (t122 * t235 + t257) * t229;
t111 = m(2) * t219 + t114;
t104 = mrSges(2,2) * t226 - mrSges(2,3) * t219 + t238 * t108 - t235 * t110 + (-t113 * t229 - t114 * t232) * pkin(6);
t103 = -mrSges(2,1) * t226 + mrSges(2,3) * t220 - pkin(1) * t113 - t229 * t106 + (t108 * t235 + t261) * t232;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t231 * t104 - t228 * t103 - qJ(1) * (t231 * t111 + t228 * t116), t104, t108, t115, t126, t139; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t228 * t104 + t231 * t103 + qJ(1) * (-t228 * t111 + t231 * t116), t103, t110, t119, t125, t138; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t246, t246, t106, t260, -Ifges(5,3) * t217 - t241, -t244;];
m_new = t1;
