% Calculate vector of inverse dynamics joint torques for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:41
% EndTime: 2019-12-31 17:25:51
% DurationCPUTime: 5.70s
% Computational Cost: add. (3397->414), mult. (7722->590), div. (0->0), fcn. (5075->10), ass. (0->197)
t154 = qJ(2) + qJ(3);
t150 = sin(t154);
t151 = cos(t154);
t155 = sin(qJ(4));
t242 = t155 * mrSges(5,2);
t303 = -t150 * t242 + t151 * (-m(5) * pkin(7) - mrSges(5,3));
t291 = t151 * pkin(3) + t150 * pkin(7);
t302 = m(5) * t291;
t301 = -t151 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t150;
t157 = sin(qJ(2));
t161 = cos(qJ(2));
t215 = qJD(1) * qJD(2);
t123 = qJDD(1) * t161 - t157 * t215;
t124 = qJDD(1) * t157 + t161 * t215;
t156 = sin(qJ(3));
t160 = cos(qJ(3));
t118 = t156 * t161 + t157 * t160;
t171 = t118 * qJD(3);
t63 = -qJD(1) * t171 + t123 * t160 - t124 * t156;
t159 = cos(qJ(4));
t252 = mrSges(5,1) * t159;
t300 = t242 - t252;
t152 = qJDD(2) + qJDD(3);
t117 = t156 * t157 - t160 * t161;
t170 = t117 * qJD(3);
t62 = -qJD(1) * t170 + t123 * t156 + t124 * t160;
t111 = t118 * qJD(1);
t153 = qJD(2) + qJD(3);
t94 = -t111 * t155 + t153 * t159;
t29 = qJD(4) * t94 + t152 * t155 + t159 * t62;
t277 = t29 / 0.2e1;
t95 = t111 * t159 + t153 * t155;
t30 = -qJD(4) * t95 + t152 * t159 - t155 * t62;
t276 = t30 / 0.2e1;
t110 = t117 * qJD(1);
t107 = qJD(4) + t110;
t264 = Ifges(5,4) * t95;
t41 = Ifges(5,2) * t94 + Ifges(5,6) * t107 + t264;
t297 = -t41 / 0.2e1;
t61 = qJDD(4) - t63;
t274 = t61 / 0.2e1;
t296 = m(4) + m(5);
t295 = t123 / 0.2e1;
t265 = t161 / 0.2e1;
t244 = t111 * mrSges(4,3);
t294 = mrSges(4,1) * t153 + mrSges(5,1) * t94 - mrSges(5,2) * t95 - t244;
t293 = t161 * Ifges(3,2);
t185 = mrSges(5,1) * t155 + mrSges(5,2) * t159;
t163 = -pkin(6) - pkin(5);
t131 = t163 * t157;
t120 = qJD(1) * t131;
t114 = qJD(2) * pkin(2) + t120;
t132 = t163 * t161;
t121 = qJD(1) * t132;
t231 = t121 * t156;
t82 = t114 * t160 + t231;
t74 = -pkin(3) * t153 - t82;
t292 = t185 * t74;
t290 = t151 * t300 + t301;
t216 = qJD(4) * t159;
t88 = -qJD(2) * t117 - t170;
t174 = t118 * t216 + t155 * t88;
t221 = qJD(1) * t161;
t222 = qJD(1) * t157;
t258 = pkin(5) * t161;
t259 = pkin(5) * t157;
t289 = (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t221) * t259 + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t222) * t258;
t115 = t123 * pkin(5);
t116 = t124 * pkin(5);
t288 = t115 * t161 + t116 * t157;
t285 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t217 = qJD(4) * t155;
t260 = pkin(2) * t161;
t146 = pkin(1) + t260;
t130 = t146 * qJD(1);
t66 = pkin(3) * t110 - pkin(7) * t111 - t130;
t224 = t160 * t121;
t83 = t114 * t156 - t224;
t75 = pkin(7) * t153 + t83;
t25 = -t155 * t75 + t159 * t66;
t238 = qJDD(1) * pkin(1);
t105 = -pkin(2) * t123 - t238;
t13 = -pkin(3) * t63 - pkin(7) * t62 + t105;
t218 = qJD(3) * t160;
t219 = qJD(3) * t156;
t91 = qJDD(2) * pkin(2) - pkin(6) * t124 - t116;
t93 = pkin(6) * t123 + t115;
t18 = t114 * t218 + t121 * t219 + t156 * t91 + t160 * t93;
t15 = pkin(7) * t152 + t18;
t26 = t155 * t66 + t159 * t75;
t3 = -qJD(4) * t26 + t13 * t159 - t15 * t155;
t255 = t155 * t3;
t284 = -t25 * t216 - t26 * t217 - t255;
t19 = -qJD(3) * t83 - t156 * t93 + t160 * t91;
t129 = -mrSges(3,1) * t161 + mrSges(3,2) * t157;
t283 = m(3) * pkin(1) + mrSges(2,1) - t129 - t301;
t11 = mrSges(5,1) * t61 - mrSges(5,3) * t29;
t12 = -mrSges(5,2) * t61 + mrSges(5,3) * t30;
t64 = -mrSges(5,2) * t107 + mrSges(5,3) * t94;
t65 = mrSges(5,1) * t107 - mrSges(5,3) * t95;
t282 = -t155 * t11 + t159 * t12 - t65 * t216 - t64 * t217;
t2 = qJD(4) * t25 + t13 * t155 + t15 * t159;
t281 = t3 * mrSges(5,1) - t2 * mrSges(5,2);
t279 = Ifges(5,1) * t277 + Ifges(5,4) * t276 + Ifges(5,5) * t274;
t278 = m(5) * pkin(3);
t273 = -t94 / 0.2e1;
t272 = -t95 / 0.2e1;
t271 = t95 / 0.2e1;
t270 = -t107 / 0.2e1;
t267 = t111 / 0.2e1;
t266 = t159 / 0.2e1;
t263 = pkin(2) * t156;
t262 = pkin(2) * t157;
t261 = pkin(2) * t160;
t254 = t159 * t2;
t250 = mrSges(4,2) * t151;
t249 = mrSges(4,3) * t110;
t248 = Ifges(3,4) * t157;
t247 = Ifges(3,4) * t161;
t246 = Ifges(5,4) * t155;
t245 = Ifges(5,4) * t159;
t243 = t111 * Ifges(4,4);
t237 = t110 * t155;
t236 = t110 * t159;
t233 = t118 * t155;
t232 = t118 * t159;
t158 = sin(qJ(1));
t228 = t155 * t158;
t162 = cos(qJ(1));
t227 = t155 * t162;
t226 = t158 * t159;
t225 = t159 * t162;
t220 = qJD(2) * t157;
t213 = Ifges(5,5) * t29 + Ifges(5,6) * t30 + Ifges(5,3) * t61;
t209 = pkin(2) * t222;
t208 = pkin(2) * t220;
t92 = Ifges(5,4) * t94;
t42 = t95 * Ifges(5,1) + t107 * Ifges(5,5) + t92;
t200 = t42 * t266;
t199 = qJD(2) * t163;
t195 = -t217 / 0.2e1;
t194 = t215 / 0.2e1;
t191 = t161 * t199;
t190 = t303 * t158;
t189 = t303 * t162;
t80 = pkin(3) * t111 + pkin(7) * t110;
t187 = mrSges(3,1) * t157 + mrSges(3,2) * t161;
t184 = Ifges(5,1) * t159 - t246;
t183 = t248 + t293;
t182 = -Ifges(5,2) * t155 + t245;
t181 = Ifges(3,5) * t161 - Ifges(3,6) * t157;
t180 = Ifges(5,5) * t159 - Ifges(5,6) * t155;
t81 = pkin(3) * t117 - pkin(7) * t118 - t146;
t97 = t131 * t156 - t132 * t160;
t44 = -t155 * t97 + t159 * t81;
t45 = t155 * t81 + t159 * t97;
t177 = t160 * t131 + t132 * t156;
t175 = pkin(1) * t187;
t173 = t118 * t217 - t159 * t88;
t172 = t157 * (Ifges(3,1) * t161 - t248);
t167 = -t255 + (-t155 * t26 - t159 * t25) * qJD(4);
t166 = m(5) * (-pkin(3) * t150 - t262) - t150 * t252;
t165 = t250 + (mrSges(4,1) + t252 + t278) * t150;
t106 = Ifges(4,4) * t110;
t16 = -pkin(3) * t152 - t19;
t40 = t95 * Ifges(5,5) + t94 * Ifges(5,6) + t107 * Ifges(5,3);
t6 = t29 * Ifges(5,4) + t30 * Ifges(5,2) + t61 * Ifges(5,6);
t71 = -t110 * Ifges(4,2) + t153 * Ifges(4,6) + t243;
t72 = t111 * Ifges(4,1) + t153 * Ifges(4,5) - t106;
t164 = Ifges(4,3) * t152 - t25 * (mrSges(5,1) * t111 + mrSges(5,3) * t236) - t26 * (-mrSges(5,2) * t111 + mrSges(5,3) * t237) - t153 * (-Ifges(4,5) * t110 - Ifges(4,6) * t111) / 0.2e1 + t130 * (mrSges(4,1) * t111 - mrSges(4,2) * t110) + (Ifges(5,5) * t111 - t110 * t184) * t272 + (Ifges(5,6) * t111 - t110 * t182) * t273 + (Ifges(5,3) * t111 - t110 * t180) * t270 + Ifges(4,5) * t62 + Ifges(4,6) * t63 - t18 * mrSges(4,2) + t19 * mrSges(4,1) + t16 * t300 + (Ifges(5,5) * t155 + Ifges(5,6) * t159) * t274 + (Ifges(5,2) * t159 + t246) * t276 + (Ifges(5,1) * t155 + t245) * t277 + t155 * t279 + t110 * t292 + t237 * t297 + mrSges(5,3) * t254 + t6 * t266 + t71 * t267 - t82 * t249 + (t292 + t200) * qJD(4) + (t107 * t180 + t182 * t94 + t184 * t95) * qJD(4) / 0.2e1 - (-Ifges(4,1) * t110 - t243 + t40) * t111 / 0.2e1 + (-Ifges(4,2) * t111 - t106 + t72) * t110 / 0.2e1 + t41 * t195 + t42 * t236 / 0.2e1;
t148 = Ifges(3,4) * t221;
t145 = -pkin(3) - t261;
t122 = t157 * t199;
t109 = Ifges(3,1) * t222 + Ifges(3,5) * qJD(2) + t148;
t108 = Ifges(3,6) * qJD(2) + qJD(1) * t183;
t104 = t151 * t225 + t228;
t103 = -t151 * t227 + t226;
t102 = -t151 * t226 + t227;
t101 = t151 * t228 + t225;
t98 = -mrSges(4,2) * t153 - t249;
t89 = qJD(2) * t118 + t171;
t86 = t120 * t160 + t231;
t85 = t120 * t156 - t224;
t79 = mrSges(4,1) * t110 + mrSges(4,2) * t111;
t70 = t80 + t209;
t53 = -mrSges(4,2) * t152 + mrSges(4,3) * t63;
t52 = mrSges(4,1) * t152 - mrSges(4,3) * t62;
t47 = qJD(3) * t177 + t160 * t122 + t156 * t191;
t43 = pkin(3) * t89 - pkin(7) * t88 + t208;
t39 = t155 * t80 + t159 * t82;
t38 = -t155 * t82 + t159 * t80;
t36 = t155 * t70 + t159 * t86;
t35 = -t155 * t86 + t159 * t70;
t10 = -mrSges(5,1) * t30 + mrSges(5,2) * t29;
t9 = -qJD(4) * t45 - t155 * t47 + t159 * t43;
t8 = qJD(4) * t44 + t155 * t43 + t159 * t47;
t1 = [t124 * t247 / 0.2e1 + t153 * (Ifges(4,5) * t88 - Ifges(4,6) * t89) / 0.2e1 - t146 * (-mrSges(4,1) * t63 + mrSges(4,2) * t62) - t130 * (mrSges(4,1) * t89 + mrSges(4,2) * t88) - (-m(4) * t19 + m(5) * t16 + t10 - t52) * t177 + m(4) * (-t105 * t146 - t130 * t208 + t18 * t97 + t47 * t83) - pkin(1) * (-mrSges(3,1) * t123 + mrSges(3,2) * t124) - t110 * (Ifges(4,4) * t88 - Ifges(4,2) * t89) / 0.2e1 + t97 * t53 + t47 * t98 + t88 * t72 / 0.2e1 + t89 * t40 / 0.2e1 - t89 * t71 / 0.2e1 + t8 * t64 + t9 * t65 + (-t104 * mrSges(5,1) - t103 * mrSges(5,2) - t296 * (t162 * t146 - t158 * t163) + t285 * t158 + (-t283 - t302) * t162) * g(2) + t44 * t11 + t45 * t12 + t25 * mrSges(5,1) * t89 - t26 * mrSges(5,2) * t89 + (Ifges(3,4) * t124 + Ifges(3,2) * t123) * t265 + (t173 * t25 - t174 * t26 - t2 * t233 - t232 * t3) * mrSges(5,3) + (t161 * t247 + t172) * t194 + (t123 * t258 + t124 * t259 + t288) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t288) + (t109 * t265 + t181 * qJD(2) / 0.2e1 - t289) * qJD(2) + t79 * t208 + m(5) * (t2 * t45 + t25 * t9 + t26 * t8 + t3 * t44) + (-t18 * mrSges(4,3) - Ifges(4,4) * t62 - Ifges(4,2) * t63 - Ifges(4,6) * t152 + t105 * mrSges(4,1) + Ifges(5,3) * t274 + Ifges(5,6) * t276 + Ifges(5,5) * t277 + t213 / 0.2e1 + t281) * t117 + (-Ifges(5,1) * t173 - Ifges(5,4) * t174 + Ifges(5,5) * t89) * t271 + t232 * t279 + t183 * t295 + t174 * t297 + (Ifges(4,1) * t88 - Ifges(4,4) * t89) * t267 + (-t82 * t88 - t83 * t89) * mrSges(4,3) + t74 * (mrSges(5,1) * t174 - mrSges(5,2) * t173) + t107 * (-Ifges(5,5) * t173 - Ifges(5,6) * t174 + Ifges(5,3) * t89) / 0.2e1 + t94 * (-Ifges(5,4) * t173 - Ifges(5,2) * t174 + Ifges(5,6) * t89) / 0.2e1 + (t105 * mrSges(4,2) - t19 * mrSges(4,3) + Ifges(4,1) * t62 + Ifges(4,4) * t63 + Ifges(4,5) * t152 + t16 * t185 + t180 * t274 + t182 * t276 + t184 * t277 + t195 * t42) * t118 + (-t102 * mrSges(5,1) - t101 * mrSges(5,2) + (t163 * t296 + t285) * t162 + (m(4) * t146 - m(5) * (-t146 - t291) + t283) * t158) * g(1) + (-m(4) * t82 + m(5) * t74 - t294) * (qJD(3) * t97 + t122 * t156 - t160 * t191) + (-mrSges(3,1) * t259 - mrSges(3,2) * t258 + 0.2e1 * Ifges(3,6) * t265) * qJDD(2) + (Ifges(3,1) * t124 + Ifges(3,4) * t295 + Ifges(3,5) * qJDD(2) - t194 * t293) * t157 + t88 * t200 - t175 * t215 - t108 * t220 / 0.2e1 - t6 * t233 / 0.2e1 - t129 * t238 + Ifges(2,3) * qJDD(1); t145 * t10 - (-Ifges(3,2) * t222 + t109 + t148) * t221 / 0.2e1 + Ifges(3,6) * t123 + Ifges(3,5) * t124 - t115 * mrSges(3,2) - t116 * mrSges(3,1) - t86 * t98 - t36 * t64 - t35 * t65 + t284 * mrSges(5,3) + (m(5) * (t167 + t254) + t282) * (pkin(7) + t263) - g(1) * (t162 * t166 - t189) - g(2) * (t158 * t166 - t190) + t52 * t261 + t53 * t263 + t83 * t244 + (m(4) * t262 + mrSges(4,1) * t150 + t187 + t250) * (g(1) * t162 + g(2) * t158) + (t129 - m(4) * t260 - m(5) * (t291 + t260) + t290) * g(3) + t294 * t85 + t164 - t79 * t209 - t181 * t215 / 0.2e1 + t108 * t222 / 0.2e1 + (t130 * t209 + t82 * t85 - t83 * t86) * m(4) + (t145 * t16 - t25 * t35 - t26 * t36 - t74 * t85) * m(5) + ((-t155 * t65 + t159 * t64 + t98) * t218 + (t156 * t18 + t160 * t19 + (-t156 * t82 + t160 * t83) * qJD(3)) * m(4) - t294 * t219 + (t156 * t74 + (-t155 * t25 + t159 * t26) * t160) * qJD(3) * m(5)) * pkin(2) + Ifges(3,3) * qJDD(2) + (t289 + (t175 - t172 / 0.2e1) * qJD(1)) * qJD(1); -t82 * t98 - t39 * t64 - t38 * t65 - pkin(3) * t10 + (m(5) * (t254 + t284) + t282) * pkin(7) + t167 * mrSges(5,3) + (t162 * t165 + t189) * g(1) + (t158 * t165 + t190) * g(2) + (t244 + t294) * t83 + t164 + (t290 - t302) * g(3) - t16 * t278 - m(5) * (t25 * t38 + t26 * t39 + t74 * t83); -t74 * (mrSges(5,1) * t95 + mrSges(5,2) * t94) + (Ifges(5,1) * t94 - t264) * t272 + t41 * t271 + (Ifges(5,5) * t94 - Ifges(5,6) * t95) * t270 - t25 * t64 + t26 * t65 - g(1) * (mrSges(5,1) * t103 - mrSges(5,2) * t104) - g(2) * (-mrSges(5,1) * t101 + mrSges(5,2) * t102) + g(3) * t185 * t150 + (t25 * t94 + t26 * t95) * mrSges(5,3) + t213 + (-Ifges(5,2) * t95 + t42 + t92) * t273 + t281;];
tau = t1;
