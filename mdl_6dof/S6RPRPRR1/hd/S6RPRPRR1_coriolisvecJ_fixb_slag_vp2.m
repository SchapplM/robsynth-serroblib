% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:00
% EndTime: 2018-11-23 16:03:05
% DurationCPUTime: 4.63s
% Computational Cost: add. (9463->455), mult. (23122->622), div. (0->0), fcn. (16769->10), ass. (0->213)
t189 = cos(qJ(3));
t235 = -cos(pkin(10)) * pkin(1) - pkin(2);
t164 = -pkin(3) * t189 + t235;
t245 = qJD(1) * t164;
t152 = qJD(4) + t245;
t181 = sin(pkin(11));
t183 = cos(pkin(11));
t187 = sin(qJ(3));
t162 = -t181 * t187 + t183 * t189;
t154 = t162 * qJD(1);
t108 = -t154 * pkin(4) + t152;
t180 = qJD(3) + qJD(5);
t185 = sin(qJ(6));
t188 = cos(qJ(6));
t186 = sin(qJ(5));
t280 = cos(qJ(5));
t243 = qJD(1) * t189;
t244 = qJD(1) * t187;
t155 = -t181 * t243 - t183 * t244;
t273 = pkin(8) * t155;
t174 = sin(pkin(10)) * pkin(1) + pkin(7);
t165 = t174 * qJD(1);
t227 = qJ(4) * qJD(1) + t165;
t240 = t187 * qJD(2);
t124 = t189 * t227 + t240;
t114 = t181 * t124;
t179 = t189 * qJD(2);
t123 = -t227 * t187 + t179;
t260 = qJD(3) * pkin(3);
t118 = t123 + t260;
t80 = t183 * t118 - t114;
t68 = qJD(3) * pkin(4) + t273 + t80;
t274 = pkin(8) * t154;
t248 = t183 * t124;
t81 = t181 * t118 + t248;
t70 = t81 + t274;
t33 = t186 * t68 + t280 * t70;
t31 = t180 * pkin(9) + t33;
t202 = t186 * t154 - t155 * t280;
t229 = t280 * t154 + t186 * t155;
t47 = -pkin(5) * t229 - pkin(9) * t202 + t108;
t14 = -t185 * t31 + t188 * t47;
t15 = t185 * t47 + t188 * t31;
t215 = t14 * t188 + t15 * t185;
t203 = t215 * mrSges(7,3);
t268 = Ifges(6,1) * t202;
t95 = Ifges(6,4) * t229;
t295 = t268 / 0.2e1 + t95 / 0.2e1;
t217 = Ifges(7,5) * t188 - Ifges(7,6) * t185;
t264 = Ifges(7,4) * t188;
t219 = -Ifges(7,2) * t185 + t264;
t265 = Ifges(7,4) * t185;
t221 = Ifges(7,1) * t188 - t265;
t222 = mrSges(7,1) * t185 + mrSges(7,2) * t188;
t281 = t188 / 0.2e1;
t282 = -t185 / 0.2e1;
t92 = t180 * t185 + t188 * t202;
t289 = t92 / 0.2e1;
t32 = -t186 * t70 + t280 * t68;
t30 = -t180 * pkin(5) - t32;
t279 = Ifges(7,4) * t92;
t91 = t180 * t188 - t185 * t202;
t96 = qJD(6) - t229;
t41 = Ifges(7,2) * t91 + Ifges(7,6) * t96 + t279;
t90 = Ifges(7,4) * t91;
t42 = Ifges(7,1) * t92 + Ifges(7,5) * t96 + t90;
t301 = t30 * t222 + t41 * t282 + t42 * t281 + t91 * t219 / 0.2e1 + t221 * t289 + t96 * t217 / 0.2e1;
t302 = -t108 * mrSges(6,2) - Ifges(6,5) * t180 + t203 - t295 - t301;
t175 = pkin(3) * t183 + pkin(4);
t275 = pkin(3) * t181;
t151 = t186 * t175 + t280 * t275;
t82 = -t123 * t181 - t248;
t207 = t82 - t274;
t83 = t183 * t123 - t114;
t75 = t83 + t273;
t300 = -t151 * qJD(5) + t186 * t75 - t207 * t280;
t299 = -t165 * t187 + t179;
t298 = -t95 / 0.2e1 + t302;
t297 = -t14 * t185 + t15 * t188;
t157 = t162 * qJD(3);
t141 = qJD(1) * t157;
t143 = t165 * t189 + t240;
t163 = t181 * t189 + t183 * t187;
t239 = qJD(1) * qJD(4);
t71 = -t163 * t239 + (-t181 * (-qJ(4) * t244 + t299) + t183 * (-qJ(4) * t243 - t143)) * qJD(3);
t192 = -t141 * pkin(8) + t71;
t156 = t163 * qJD(3);
t140 = qJD(1) * t156;
t72 = t183 * (qJD(3) * t123 + t189 * t239) + t181 * (-qJD(3) * t124 - t187 * t239);
t54 = -pkin(8) * t140 + t72;
t12 = qJD(5) * t32 + t186 * t192 + t280 * t54;
t177 = pkin(3) * t244;
t173 = qJD(3) * t177;
t113 = pkin(4) * t140 + t173;
t64 = qJD(5) * t229 - t186 * t140 + t280 * t141;
t65 = qJD(5) * t202 + t140 * t280 + t186 * t141;
t27 = pkin(5) * t65 - pkin(9) * t64 + t113;
t2 = qJD(6) * t14 + t12 * t188 + t185 * t27;
t3 = -qJD(6) * t15 - t12 * t185 + t188 * t27;
t45 = qJD(6) * t91 + t188 * t64;
t46 = -qJD(6) * t92 - t185 * t64;
t296 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t45 + Ifges(7,6) * t46;
t67 = pkin(5) * t202 - pkin(9) * t229;
t294 = t45 / 0.2e1;
t293 = t46 / 0.2e1;
t292 = t65 / 0.2e1;
t291 = -t91 / 0.2e1;
t290 = -t92 / 0.2e1;
t287 = -t96 / 0.2e1;
t285 = -t155 / 0.2e1;
t284 = -t156 / 0.2e1;
t283 = t157 / 0.2e1;
t13 = qJD(5) * t33 + t186 * t54 - t192 * t280;
t246 = qJ(4) + t174;
t160 = t246 * t187;
t161 = t246 * t189;
t105 = -t183 * t160 - t161 * t181;
t88 = -pkin(8) * t163 + t105;
t106 = -t181 * t160 + t183 * t161;
t89 = pkin(8) * t162 + t106;
t205 = -t186 * t89 + t280 * t88;
t272 = t13 * t205;
t271 = t188 * t2;
t270 = t3 * t185;
t269 = -mrSges(6,1) * t180 - mrSges(7,1) * t91 + mrSges(7,2) * t92 + mrSges(6,3) * t202;
t267 = Ifges(4,4) * t187;
t262 = Ifges(6,2) * t229;
t201 = t162 * t280 - t186 * t163;
t259 = t201 * t13;
t256 = t155 * Ifges(5,4);
t167 = t235 * qJD(1);
t255 = t167 * mrSges(4,2);
t253 = Ifges(4,5) * qJD(3);
t252 = Ifges(4,6) * qJD(3);
t251 = t162 * t141;
t250 = t163 * t140;
t228 = qJD(3) * t246;
t127 = qJD(4) * t189 - t187 * t228;
t128 = -t187 * qJD(4) - t189 * t228;
t87 = t183 * t127 + t181 * t128;
t242 = qJD(6) * t185;
t241 = qJD(6) * t188;
t176 = Ifges(4,4) * t243;
t234 = t253 / 0.2e1;
t233 = -t252 / 0.2e1;
t232 = m(4) * t174 + mrSges(4,3);
t231 = t65 * mrSges(6,1) + t64 * mrSges(6,2);
t125 = -pkin(4) * t155 + t177;
t126 = pkin(4) * t156 + t187 * t260;
t230 = t140 * mrSges(5,1) + t141 * mrSges(5,2);
t86 = -t127 * t181 + t183 * t128;
t226 = t252 / 0.2e1 + (Ifges(4,2) * t189 + t267) * qJD(1) / 0.2e1 - t167 * mrSges(4,1);
t224 = -t185 * t2 - t188 * t3;
t223 = mrSges(7,1) * t188 - mrSges(7,2) * t185;
t220 = Ifges(7,1) * t185 + t264;
t218 = Ifges(7,2) * t188 + t265;
t216 = Ifges(7,5) * t185 + Ifges(7,6) * t188;
t23 = mrSges(7,1) * t65 - mrSges(7,3) * t45;
t24 = -mrSges(7,2) * t65 + mrSges(7,3) * t46;
t213 = -t185 * t23 + t188 * t24;
t49 = t186 * t88 + t280 * t89;
t110 = t186 * t162 + t163 * t280;
t122 = -t162 * pkin(4) + t164;
t55 = -pkin(5) * t201 - t110 * pkin(9) + t122;
t26 = t185 * t55 + t188 * t49;
t25 = -t185 * t49 + t188 * t55;
t56 = -mrSges(7,2) * t96 + mrSges(7,3) * t91;
t57 = mrSges(7,1) * t96 - mrSges(7,3) * t92;
t212 = -t185 * t57 + t188 * t56;
t211 = -t185 * t56 - t188 * t57;
t93 = -mrSges(6,2) * t180 + mrSges(6,3) * t229;
t208 = -t212 - t93;
t206 = -pkin(8) * t157 + t86;
t150 = t175 * t280 - t186 * t275;
t197 = -qJD(6) * t215 - t270;
t196 = t197 + t271;
t10 = t45 * Ifges(7,1) + t46 * Ifges(7,4) + t65 * Ifges(7,5);
t9 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t65 * Ifges(7,6);
t195 = -t12 * mrSges(6,2) + mrSges(7,3) * t271 + t220 * t294 + t218 * t293 + t216 * t292 + t185 * t10 / 0.2e1 - Ifges(6,6) * t65 + Ifges(6,5) * t64 + t9 * t281 + (-mrSges(6,1) - t223) * t13 + t301 * qJD(6);
t194 = t15 * mrSges(7,2) - Ifges(7,3) * t96 - Ifges(7,6) * t91 - Ifges(7,5) * t92 + Ifges(6,6) * t180 + t262 / 0.2e1 + Ifges(6,4) * t202 - t108 * mrSges(6,1) - t14 * mrSges(7,1);
t193 = -t262 / 0.2e1 - t194;
t168 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t243;
t166 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t244;
t159 = Ifges(4,1) * t244 + t176 + t253;
t149 = Ifges(5,4) * t154;
t148 = pkin(9) + t151;
t147 = -pkin(5) - t150;
t135 = t150 * qJD(5);
t133 = t143 * qJD(3);
t132 = t299 * qJD(3);
t131 = qJD(3) * mrSges(5,1) + t155 * mrSges(5,3);
t130 = -qJD(3) * mrSges(5,2) + t154 * mrSges(5,3);
t107 = -mrSges(5,1) * t154 - mrSges(5,2) * t155;
t98 = -t155 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t149;
t97 = t154 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t256;
t77 = -pkin(8) * t156 + t87;
t74 = qJD(5) * t110 + t156 * t280 + t186 * t157;
t73 = qJD(5) * t201 - t186 * t156 + t157 * t280;
t66 = -mrSges(6,1) * t229 + mrSges(6,2) * t202;
t61 = Ifges(7,3) * t65;
t50 = t125 + t67;
t35 = t186 * t207 + t280 * t75;
t29 = pkin(5) * t74 - pkin(9) * t73 + t126;
t22 = t185 * t67 + t188 * t32;
t21 = -t185 * t32 + t188 * t67;
t20 = qJD(5) * t49 + t186 * t77 - t206 * t280;
t19 = qJD(5) * t205 + t186 * t206 + t280 * t77;
t18 = t185 * t50 + t188 * t35;
t17 = -t185 * t35 + t188 * t50;
t16 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t5 = -qJD(6) * t26 - t185 * t19 + t188 * t29;
t4 = qJD(6) * t25 + t185 * t29 + t188 * t19;
t1 = [(t295 - t302) * t73 + m(6) * (t108 * t126 + t113 * t122 + t12 * t49 + t19 * t33 - t20 * t32 - t272) + m(7) * (t14 * t5 + t15 * t4 + t2 * t26 + t20 * t30 + t25 * t3 - t272) + (-t205 * t64 - t32 * t73 - t33 * t74 - t49 * t65) * mrSges(6,3) - t205 * t16 - (-t12 * mrSges(6,3) + t61 / 0.2e1 - Ifges(6,4) * t64 + t113 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t65 + t296) * t201 + (t232 * t132 + (-t174 * t166 + t159 / 0.2e1 + 0.3e1 / 0.2e1 * t176 + t234 - t232 * t299 + 0.2e1 * t255) * qJD(3)) * t189 + (-t105 * t141 - t106 * t140 - t156 * t81 - t157 * t80 + t162 * t72 - t163 * t71) * mrSges(5,3) + (-t162 * t140 + t154 * t284) * Ifges(5,2) + (t154 * t283 - t156 * t285 - t250 + t251) * Ifges(5,4) + m(5) * (t105 * t71 + t106 * t72 + t80 * t86 + t81 * t87) + (t232 * t133 + (-t174 * t168 + t233 - t232 * t143 + (t235 * mrSges(4,1) - 0.3e1 / 0.2e1 * t267 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t189) * qJD(1) + (t107 + qJD(1) * (-mrSges(5,1) * t162 + mrSges(5,2) * t163) + m(5) * (t152 + t245)) * pkin(3) - t226) * qJD(3)) * t187 + qJD(3) * (Ifges(5,5) * t157 - Ifges(5,6) * t156) / 0.2e1 + t152 * (mrSges(5,1) * t156 + mrSges(5,2) * t157) + (t221 * t294 + t219 * t293 + t217 * t292 + t9 * t282 + t10 * t281 + Ifges(6,1) * t64 - Ifges(6,4) * t65 + t113 * mrSges(6,2) + (mrSges(6,3) + t222) * t13 + t224 * mrSges(7,3) + (t42 * t282 - t188 * t41 / 0.2e1 + t30 * t223 + t218 * t291 + t220 * t290 + t216 * t287 - t297 * mrSges(7,3)) * qJD(6)) * t110 + t193 * t74 + t98 * t283 + t97 * t284 + t164 * t230 + t122 * t231 + t25 * t23 + t26 * t24 + t269 * t20 + t4 * t56 + t5 * t57 + (t163 * t141 + t157 * t285) * Ifges(5,1) + t19 * t93 + t126 * t66 + t87 * t130 + t86 * t131; t157 * t130 - t156 * t131 + t269 * t74 - (t64 * mrSges(6,3) + t16) * t201 + (-t250 - t251) * mrSges(5,3) - t208 * t73 + (-t187 * t166 + t189 * t168 + (-t187 ^ 2 - t189 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t65 * mrSges(6,3) + t211 * qJD(6) + t213) * t110 + m(7) * (t110 * t196 + t297 * t73 + t30 * t74 - t259) + m(4) * (t132 * t187 - t133 * t189 + (t143 * t189 - t187 * t299) * qJD(3)) + m(6) * (t110 * t12 - t32 * t74 + t33 * t73 - t259) + m(5) * (-t156 * t80 + t157 * t81 + t162 * t71 + t163 * t72); (-t268 / 0.2e1 + t298) * t229 + ((t299 * mrSges(4,3) + t234 - t255 - t159 / 0.2e1 - t176 / 0.2e1) * t189 + (t143 * mrSges(4,3) + t233 + (t267 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t189) * qJD(1) + (-m(5) * t152 - t107) * pkin(3) + t226) * t187) * qJD(1) - t299 * t168 + (t154 * t80 - t155 * t81 + (-t140 * t181 - t141 * t183) * pkin(3)) * mrSges(5,3) + (-t150 * t64 - t151 * t65 + t202 * t33 + t229 * t32) * mrSges(6,3) + (t148 * t211 - t203) * qJD(6) - (Ifges(5,2) * t155 + t149 + t98) * t154 / 0.2e1 - t193 * t202 + (-t108 * t125 + t12 * t151 - t13 * t150 + (-t35 + t135) * t33 + t300 * t32) * m(6) + (t13 * t147 + t135 * t297 - t14 * t17 + t148 * t196 - t15 * t18 - t30 * t300) * m(7) - t300 * t269 + t195 - m(5) * (t80 * t82 + t81 * t83) + m(5) * (t181 * t72 + t183 * t71) * pkin(3) - t208 * t135 + t213 * t148 + t97 * t285 + t155 * (Ifges(5,1) * t154 + t256) / 0.2e1 - t18 * t56 - t17 * t57 - mrSges(7,3) * t270 + t71 * mrSges(5,1) - t72 * mrSges(5,2) - t35 * t93 - t125 * t66 - t83 * t130 - t82 * t131 - t132 * mrSges(4,2) - t133 * mrSges(4,1) - Ifges(5,6) * t140 + Ifges(5,5) * t141 + t147 * t16 - qJD(3) * (Ifges(5,5) * t154 + Ifges(5,6) * t155) / 0.2e1 - t152 * (-mrSges(5,1) * t155 + mrSges(5,2) * t154) + t143 * t166; t212 * qJD(6) + t208 * t229 - t154 * t130 - t155 * t131 + t185 * t24 + t188 * t23 - t269 * t202 + t230 + t231 + (-t30 * t202 + t297 * t96 - t224) * m(7) + (t202 * t32 - t229 * t33 + t113) * m(6) + (-t154 * t81 - t155 * t80 + t173) * m(5); (t32 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t202 + t298) * t229 + (t33 * mrSges(6,3) + t194) * t202 + (-t57 * t241 - t56 * t242 + t213) * pkin(9) + t195 + t197 * mrSges(7,3) - t269 * t33 - pkin(5) * t16 - t22 * t56 - t21 * t57 - t32 * t93 + (-t14 * t21 - t15 * t22 - t30 * t33 + (-t14 * t241 - t15 * t242 - t270 + t271) * pkin(9) - pkin(5) * t13) * m(7); t61 - t30 * (mrSges(7,1) * t92 + mrSges(7,2) * t91) + (Ifges(7,1) * t91 - t279) * t290 + t41 * t289 + (Ifges(7,5) * t91 - Ifges(7,6) * t92) * t287 - t14 * t56 + t15 * t57 + (t14 * t91 + t15 * t92) * mrSges(7,3) + (-Ifges(7,2) * t92 + t42 + t90) * t291 + t296;];
tauc  = t1(:);
