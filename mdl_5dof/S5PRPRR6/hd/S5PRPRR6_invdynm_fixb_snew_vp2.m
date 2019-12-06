% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:39
% EndTime: 2019-12-05 15:56:48
% DurationCPUTime: 6.40s
% Computational Cost: add. (102906->245), mult. (223626->318), div. (0->0), fcn. (162434->12), ass. (0->118)
t242 = qJD(2) ^ 2;
t230 = sin(pkin(9));
t233 = cos(pkin(9));
t217 = t230 * g(1) - t233 * g(2);
t218 = -t233 * g(1) - t230 * g(2);
t228 = -g(3) + qJDD(1);
t240 = cos(qJ(2));
t234 = cos(pkin(5));
t237 = sin(qJ(2));
t268 = t234 * t237;
t231 = sin(pkin(5));
t269 = t231 * t237;
t189 = t217 * t268 + t240 * t218 + t228 * t269;
t184 = -t242 * pkin(2) + qJDD(2) * qJ(3) + t189;
t229 = sin(pkin(10));
t205 = -t231 * t217 + t234 * t228;
t232 = cos(pkin(10));
t263 = qJD(2) * qJD(3);
t267 = t232 * t205 - 0.2e1 * t229 * t263;
t275 = pkin(3) * t232;
t166 = (-pkin(7) * qJDD(2) + t242 * t275 - t184) * t229 + t267;
t169 = t229 * t205 + (t184 + 0.2e1 * t263) * t232;
t262 = qJDD(2) * t232;
t226 = t232 ^ 2;
t270 = t226 * t242;
t167 = -pkin(3) * t270 + pkin(7) * t262 + t169;
t236 = sin(qJ(4));
t239 = cos(qJ(4));
t162 = t236 * t166 + t239 * t167;
t252 = t229 * t236 - t232 * t239;
t207 = t252 * qJD(2);
t253 = t229 * t239 + t232 * t236;
t208 = t253 * qJD(2);
t191 = t207 * mrSges(5,1) + t208 * mrSges(5,2);
t264 = t208 * qJD(4);
t197 = -t252 * qJDD(2) - t264;
t204 = qJD(4) * mrSges(5,1) - t208 * mrSges(5,3);
t196 = t207 * pkin(4) - t208 * pkin(8);
t241 = qJD(4) ^ 2;
t159 = -pkin(4) * t241 + qJDD(4) * pkin(8) - t196 * t207 + t162;
t225 = t229 ^ 2;
t188 = -t237 * t218 + (t217 * t234 + t228 * t231) * t240;
t248 = qJDD(3) - t188;
t177 = (-pkin(2) - t275) * qJDD(2) + (-qJ(3) + (-t225 - t226) * pkin(7)) * t242 + t248;
t265 = t207 * qJD(4);
t198 = t253 * qJDD(2) - t265;
t163 = (-t198 + t265) * pkin(8) + (-t197 + t264) * pkin(4) + t177;
t235 = sin(qJ(5));
t238 = cos(qJ(5));
t156 = -t159 * t235 + t163 * t238;
t199 = t238 * qJD(4) - t235 * t208;
t176 = t199 * qJD(5) + t235 * qJDD(4) + t238 * t198;
t200 = t235 * qJD(4) + t238 * t208;
t179 = -t199 * mrSges(6,1) + t200 * mrSges(6,2);
t206 = qJD(5) + t207;
t182 = -t206 * mrSges(6,2) + t199 * mrSges(6,3);
t195 = qJDD(5) - t197;
t152 = m(6) * t156 + mrSges(6,1) * t195 - t176 * mrSges(6,3) - t179 * t200 + t182 * t206;
t157 = t159 * t238 + t163 * t235;
t175 = -t200 * qJD(5) + t238 * qJDD(4) - t235 * t198;
t183 = t206 * mrSges(6,1) - t200 * mrSges(6,3);
t153 = m(6) * t157 - mrSges(6,2) * t195 + t175 * mrSges(6,3) + t179 * t199 - t183 * t206;
t259 = -t152 * t235 + t238 * t153;
t139 = m(5) * t162 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t197 - qJD(4) * t204 - t191 * t207 + t259;
t161 = t166 * t239 - t167 * t236;
t203 = -qJD(4) * mrSges(5,2) - t207 * mrSges(5,3);
t158 = -qJDD(4) * pkin(4) - pkin(8) * t241 + t196 * t208 - t161;
t249 = -m(6) * t158 + t175 * mrSges(6,1) - t176 * mrSges(6,2) + t199 * t182 - t183 * t200;
t148 = m(5) * t161 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t198 + qJD(4) * t203 - t191 * t208 + t249;
t134 = t236 * t139 + t239 * t148;
t168 = -t229 * t184 + t267;
t273 = mrSges(4,2) * t229;
t251 = mrSges(4,3) * qJDD(2) + t242 * (-mrSges(4,1) * t232 + t273);
t132 = m(4) * t168 - t251 * t229 + t134;
t260 = t239 * t139 - t236 * t148;
t133 = m(4) * t169 + t251 * t232 + t260;
t126 = t232 * t132 + t229 * t133;
t256 = Ifges(4,5) * t229 + Ifges(4,6) * t232;
t170 = Ifges(6,5) * t200 + Ifges(6,6) * t199 + Ifges(6,3) * t206;
t172 = Ifges(6,1) * t200 + Ifges(6,4) * t199 + Ifges(6,5) * t206;
t145 = -mrSges(6,1) * t158 + mrSges(6,3) * t157 + Ifges(6,4) * t176 + Ifges(6,2) * t175 + Ifges(6,6) * t195 - t170 * t200 + t172 * t206;
t171 = Ifges(6,4) * t200 + Ifges(6,2) * t199 + Ifges(6,6) * t206;
t146 = mrSges(6,2) * t158 - mrSges(6,3) * t156 + Ifges(6,1) * t176 + Ifges(6,4) * t175 + Ifges(6,5) * t195 + t170 * t199 - t171 * t206;
t186 = Ifges(5,4) * t208 - Ifges(5,2) * t207 + Ifges(5,6) * qJD(4);
t187 = Ifges(5,1) * t208 - Ifges(5,4) * t207 + Ifges(5,5) * qJD(4);
t246 = -mrSges(5,1) * t161 + mrSges(5,2) * t162 - Ifges(5,5) * t198 - Ifges(5,6) * t197 - Ifges(5,3) * qJDD(4) - pkin(4) * t249 - pkin(8) * t259 - t238 * t145 - t235 * t146 - t208 * t186 - t207 * t187;
t257 = Ifges(4,4) * t229 + Ifges(4,2) * t232;
t258 = Ifges(4,1) * t229 + Ifges(4,4) * t232;
t276 = -mrSges(4,1) * t168 + mrSges(4,2) * t169 - pkin(3) * t134 - (t229 * t257 - t232 * t258) * t242 + t246;
t112 = (Ifges(3,6) - t256) * qJDD(2) + t242 * Ifges(3,5) - mrSges(3,1) * t205 + mrSges(3,3) * t189 - pkin(2) * t126 + t276;
t261 = -t132 * t229 + t232 * t133;
t124 = m(3) * t189 - mrSges(3,1) * t242 - qJDD(2) * mrSges(3,2) + t261;
t181 = -qJDD(2) * pkin(2) - t242 * qJ(3) + t248;
t141 = t238 * t152 + t235 * t153;
t247 = m(5) * t177 - t197 * mrSges(5,1) + t198 * mrSges(5,2) + t207 * t203 + t208 * t204 + t141;
t245 = -m(4) * t181 + mrSges(4,1) * t262 - t247 + (t225 * t242 + t270) * mrSges(4,3);
t136 = t245 + (mrSges(3,1) - t273) * qJDD(2) - t242 * mrSges(3,2) + m(3) * t188;
t121 = t240 * t124 - t136 * t237;
t277 = pkin(6) * t121 + t112 * t240;
t271 = t136 * t240;
t266 = t242 * t256;
t125 = m(3) * t205 + t126;
t116 = t124 * t268 - t125 * t231 + t234 * t271;
t185 = Ifges(5,5) * t208 - Ifges(5,6) * t207 + Ifges(5,3) * qJD(4);
t127 = mrSges(5,2) * t177 - mrSges(5,3) * t161 + Ifges(5,1) * t198 + Ifges(5,4) * t197 + Ifges(5,5) * qJDD(4) - pkin(8) * t141 - qJD(4) * t186 - t145 * t235 + t146 * t238 - t185 * t207;
t244 = mrSges(6,1) * t156 - mrSges(6,2) * t157 + Ifges(6,5) * t176 + Ifges(6,6) * t175 + Ifges(6,3) * t195 + t171 * t200 - t172 * t199;
t128 = -mrSges(5,1) * t177 + mrSges(5,3) * t162 + Ifges(5,4) * t198 + Ifges(5,2) * t197 + Ifges(5,6) * qJDD(4) - pkin(4) * t141 + qJD(4) * t187 - t185 * t208 - t244;
t117 = -mrSges(4,1) * t181 + mrSges(4,3) * t169 - pkin(3) * t247 + pkin(7) * t260 + t257 * qJDD(2) + t236 * t127 + t239 * t128 - t229 * t266;
t118 = mrSges(4,2) * t181 - mrSges(4,3) * t168 - pkin(7) * t134 + t258 * qJDD(2) + t239 * t127 - t236 * t128 + t232 * t266;
t108 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t188 - mrSges(3,2) * t189 + t229 * t118 + t232 * t117 + pkin(2) * (-qJDD(2) * t273 + t245) + qJ(3) * t261;
t110 = mrSges(3,2) * t205 - mrSges(3,3) * t188 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t242 - qJ(3) * t126 - t117 * t229 + t118 * t232;
t250 = mrSges(2,1) * t217 - mrSges(2,2) * t218 + pkin(1) * t116 + t234 * t108 + t110 * t269 + t277 * t231;
t119 = m(2) * t218 + t121;
t115 = t234 * t125 + (t124 * t237 + t271) * t231;
t113 = m(2) * t217 + t116;
t106 = mrSges(2,2) * t228 - mrSges(2,3) * t217 + t240 * t110 - t237 * t112 + (-t115 * t231 - t116 * t234) * pkin(6);
t105 = -mrSges(2,1) * t228 + mrSges(2,3) * t218 - pkin(1) * t115 - t231 * t108 + (t110 * t237 + t277) * t234;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t233 * t106 - t230 * t105 - qJ(1) * (t113 * t233 + t119 * t230), t106, t110, t118, t127, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t230 * t106 + t233 * t105 + qJ(1) * (-t113 * t230 + t119 * t233), t105, t112, t117, t128, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t250, t250, t108, t256 * qJDD(2) - t276, -t246, t244;];
m_new = t1;
