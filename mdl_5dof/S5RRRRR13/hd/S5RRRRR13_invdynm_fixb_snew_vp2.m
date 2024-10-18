% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR13
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
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:38
% EndTime: 2024-09-27 17:30:49
% DurationCPUTime: 5.86s
% Computational Cost: add. (249879->240), mult. (291728->321), div. (0->0), fcn. (185501->12), ass. (0->112)
t227 = sin(pkin(5));
t260 = pkin(9) * t227;
t228 = cos(pkin(5));
t259 = t228 * g(3);
t224 = qJD(1) + qJD(2);
t219 = qJD(3) + t224;
t217 = t219 ^ 2;
t258 = t217 * t227 ^ 2;
t257 = t219 * t227;
t230 = sin(qJ(4));
t256 = t227 * t230;
t235 = cos(qJ(4));
t255 = t227 * t235;
t254 = t228 * t230;
t253 = t228 * t235;
t222 = qJDD(1) + qJDD(2);
t218 = qJDD(3) + t222;
t252 = qJD(4) * t219;
t198 = (t218 * t230 + t235 * t252) * t227;
t208 = t228 * t218 + qJDD(4);
t210 = t228 * t219 + qJD(4);
t233 = sin(qJ(1));
t238 = cos(qJ(1));
t211 = t233 * g(1) - t238 * g(2);
t205 = qJDD(1) * pkin(1) + t211;
t212 = -t238 * g(1) - t233 * g(2);
t239 = qJD(1) ^ 2;
t206 = -t239 * pkin(1) + t212;
t232 = sin(qJ(2));
t237 = cos(qJ(2));
t190 = t237 * t205 - t232 * t206;
t187 = t222 * pkin(2) + t190;
t191 = t232 * t205 + t237 * t206;
t221 = t224 ^ 2;
t188 = -t221 * pkin(2) + t191;
t231 = sin(qJ(3));
t236 = cos(qJ(3));
t172 = t236 * t187 - t231 * t188;
t164 = t218 * pkin(3) + t217 * t260 + t172;
t173 = t231 * t187 + t236 * t188;
t165 = -t217 * pkin(3) + t218 * t260 + t173;
t245 = t164 * t253 - t230 * t165;
t153 = t208 * pkin(4) - t198 * pkin(10) + (pkin(4) * t230 * t258 + (pkin(10) * t210 * t219 - g(3)) * t227) * t235 + t245;
t156 = -g(3) * t256 + t164 * t254 + t235 * t165;
t250 = t219 * t256;
t197 = t210 * pkin(4) - pkin(10) * t250;
t199 = (t218 * t235 - t230 * t252) * t227;
t251 = t235 ^ 2 * t258;
t154 = -pkin(4) * t251 + t199 * pkin(10) - t210 * t197 + t156;
t229 = sin(qJ(5));
t234 = cos(qJ(5));
t151 = t234 * t153 - t229 * t154;
t192 = (-t229 * t230 + t234 * t235) * t257;
t170 = t192 * qJD(5) + t234 * t198 + t229 * t199;
t193 = (t229 * t235 + t230 * t234) * t257;
t178 = -t192 * mrSges(6,1) + t193 * mrSges(6,2);
t207 = qJD(5) + t210;
t179 = -t207 * mrSges(6,2) + t192 * mrSges(6,3);
t204 = qJDD(5) + t208;
t146 = m(6) * t151 + t204 * mrSges(6,1) - t170 * mrSges(6,3) - t193 * t178 + t207 * t179;
t152 = t229 * t153 + t234 * t154;
t169 = -t193 * qJD(5) - t229 * t198 + t234 * t199;
t180 = t207 * mrSges(6,1) - t193 * mrSges(6,3);
t147 = m(6) * t152 - t204 * mrSges(6,2) + t169 * mrSges(6,3) + t192 * t178 - t207 * t180;
t140 = t234 * t146 + t229 * t147;
t155 = -g(3) * t255 + t245;
t249 = t219 * t255;
t195 = -t210 * mrSges(5,2) + mrSges(5,3) * t249;
t196 = (-mrSges(5,1) * t235 + mrSges(5,2) * t230) * t257;
t138 = m(5) * t155 + t208 * mrSges(5,1) - t198 * mrSges(5,3) + t210 * t195 - t196 * t250 + t140;
t194 = t210 * mrSges(5,1) - mrSges(5,3) * t250;
t246 = -t229 * t146 + t234 * t147;
t139 = m(5) * t156 - t208 * mrSges(5,2) + t199 * mrSges(5,3) - t210 * t194 + t196 * t249 + t246;
t159 = -t227 * t164 - t259;
t158 = -pkin(10) * t251 - t199 * pkin(4) - t259 + (t197 * t219 * t230 - t164) * t227;
t244 = m(6) * t158 - t169 * mrSges(6,1) + t170 * mrSges(6,2) - t192 * t179 + t193 * t180;
t149 = m(5) * t159 - t199 * mrSges(5,1) + t198 * mrSges(5,2) + (t194 * t230 - t195 * t235) * t257 + t244;
t123 = t138 * t253 + t139 * t254 - t227 * t149;
t120 = m(4) * t172 + t218 * mrSges(4,1) - t217 * mrSges(4,2) + t123;
t130 = -t230 * t138 + t235 * t139;
t128 = m(4) * t173 - t217 * mrSges(4,1) - t218 * mrSges(4,2) + t130;
t116 = t236 * t120 + t231 * t128;
t113 = m(3) * t190 + t222 * mrSges(3,1) - t221 * mrSges(3,2) + t116;
t247 = -t231 * t120 + t236 * t128;
t114 = m(3) * t191 - t221 * mrSges(3,1) - t222 * mrSges(3,2) + t247;
t109 = t237 * t113 + t232 * t114;
t122 = t138 * t255 + t139 * t256 + t228 * t149;
t248 = -t232 * t113 + t237 * t114;
t174 = Ifges(6,5) * t193 + Ifges(6,6) * t192 + Ifges(6,3) * t207;
t176 = Ifges(6,1) * t193 + Ifges(6,4) * t192 + Ifges(6,5) * t207;
t141 = -mrSges(6,1) * t158 + mrSges(6,3) * t152 + Ifges(6,4) * t170 + Ifges(6,2) * t169 + Ifges(6,6) * t204 - t193 * t174 + t207 * t176;
t175 = Ifges(6,4) * t193 + Ifges(6,2) * t192 + Ifges(6,6) * t207;
t142 = mrSges(6,2) * t158 - mrSges(6,3) * t151 + Ifges(6,1) * t170 + Ifges(6,4) * t169 + Ifges(6,5) * t204 + t192 * t174 - t207 * t175;
t184 = Ifges(5,3) * t210 + (Ifges(5,5) * t230 + Ifges(5,6) * t235) * t257;
t186 = Ifges(5,5) * t210 + (Ifges(5,1) * t230 + Ifges(5,4) * t235) * t257;
t118 = -mrSges(5,1) * t159 + mrSges(5,3) * t156 + Ifges(5,4) * t198 + Ifges(5,2) * t199 + Ifges(5,6) * t208 - pkin(4) * t244 + pkin(10) * t246 + t234 * t141 + t229 * t142 - t184 * t250 + t210 * t186;
t185 = Ifges(5,6) * t210 + (Ifges(5,4) * t230 + Ifges(5,2) * t235) * t257;
t125 = mrSges(5,2) * t159 - mrSges(5,3) * t155 + Ifges(5,1) * t198 + Ifges(5,4) * t199 + Ifges(5,5) * t208 - pkin(10) * t140 - t229 * t141 + t234 * t142 + t184 * t249 - t210 * t185;
t242 = mrSges(6,1) * t151 - mrSges(6,2) * t152 + Ifges(6,5) * t170 + Ifges(6,6) * t169 + Ifges(6,3) * t204 + t193 * t175 - t192 * t176;
t132 = mrSges(5,1) * t155 - mrSges(5,2) * t156 + Ifges(5,5) * t198 + Ifges(5,6) * t199 + Ifges(5,3) * t208 + pkin(4) * t140 + (t185 * t230 - t186 * t235) * t257 + t242;
t243 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,3) * t218 + pkin(3) * t123 + t118 * t255 + t125 * t256 + t130 * t260 + t228 * t132;
t241 = mrSges(3,1) * t190 - mrSges(3,2) * t191 + Ifges(3,3) * t222 + pkin(2) * t116 + t243;
t240 = mrSges(2,1) * t211 - mrSges(2,2) * t212 + Ifges(2,3) * qJDD(1) + pkin(1) * t109 + t241;
t107 = m(2) * t212 - t239 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t248;
t106 = m(2) * t211 + qJDD(1) * mrSges(2,1) - t239 * mrSges(2,2) + t109;
t105 = -mrSges(4,2) * g(3) - mrSges(4,3) * t172 + Ifges(4,5) * t218 - t217 * Ifges(4,6) - t230 * t118 + t235 * t125 + (-t122 * t227 - t123 * t228) * pkin(9);
t104 = mrSges(4,1) * g(3) + mrSges(4,3) * t173 + t217 * Ifges(4,5) + Ifges(4,6) * t218 - pkin(3) * t122 - t227 * t132 + (pkin(9) * t130 + t118 * t235 + t125 * t230) * t228;
t103 = -mrSges(3,2) * g(3) - mrSges(3,3) * t190 + Ifges(3,5) * t222 - t221 * Ifges(3,6) - pkin(8) * t116 - t231 * t104 + t236 * t105;
t102 = Ifges(3,6) * t222 + t221 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t191 + t231 * t105 + t236 * t104 - pkin(2) * (-m(4) * g(3) + t122) + pkin(8) * t247;
t101 = -mrSges(2,2) * g(3) - mrSges(2,3) * t211 + Ifges(2,5) * qJDD(1) - t239 * Ifges(2,6) - pkin(7) * t109 - t232 * t102 + t237 * t103;
t100 = Ifges(2,6) * qJDD(1) + t239 * Ifges(2,5) + mrSges(2,3) * t212 + t232 * t103 + t237 * t102 - pkin(1) * t122 + pkin(7) * t248 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t238 * t101 - t233 * t100 - pkin(6) * (t238 * t106 + t233 * t107), t101, t103, t105, t125, t142; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t233 * t101 + t238 * t100 + pkin(6) * (-t233 * t106 + t238 * t107), t100, t102, t104, t118, t141; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t240, t240, t241, t243, t132, t242;];
m_new = t1;
