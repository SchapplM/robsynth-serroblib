% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2018-11-23 16:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:46:08
% EndTime: 2018-11-23 16:46:14
% DurationCPUTime: 5.86s
% Computational Cost: add. (3142->476), mult. (7079->598), div. (0->0), fcn. (3466->4), ass. (0->210)
t301 = Ifges(6,4) + Ifges(7,4);
t302 = Ifges(6,1) + Ifges(7,1);
t291 = Ifges(6,5) + Ifges(7,5);
t300 = Ifges(6,2) + Ifges(7,2);
t290 = Ifges(6,6) + Ifges(7,6);
t155 = sin(qJ(5));
t157 = cos(qJ(5));
t156 = sin(qJ(2));
t158 = cos(qJ(2));
t195 = pkin(4) * t156 + pkin(8) * t158;
t226 = qJD(1) * t158;
t227 = qJD(1) * t156;
t97 = -qJD(1) * pkin(1) - pkin(2) * t226 - qJ(3) * t227;
t75 = pkin(3) * t226 + qJD(4) - t97;
t48 = qJD(1) * t195 + t75;
t159 = -pkin(2) - pkin(3);
t150 = -pkin(8) + t159;
t133 = qJ(4) * t227;
t140 = pkin(7) * t227;
t218 = qJD(3) + t140;
t206 = -t133 + t218;
t73 = qJD(2) * t150 + t206;
t12 = -t155 * t73 + t157 * t48;
t13 = t155 * t48 + t157 * t73;
t174 = t12 * t157 + t13 * t155;
t185 = mrSges(7,1) * t155 + mrSges(7,2) * t157;
t187 = mrSges(6,1) * t155 + mrSges(6,2) * t157;
t224 = qJD(2) * t157;
t101 = t155 * t226 - t224;
t11 = qJ(6) * t101 + t13;
t102 = qJD(2) * t155 + t157 * t226;
t10 = qJ(6) * t102 + t12;
t131 = qJD(5) + t227;
t9 = pkin(5) * t131 + t10;
t189 = t11 * t155 + t157 * t9;
t256 = -t157 / 0.2e1;
t257 = t155 / 0.2e1;
t298 = t301 * t101;
t275 = -t102 * t302 + t131 * t291 + t298;
t297 = t301 * t102;
t276 = t101 * t300 + t131 * t290 - t297;
t141 = pkin(7) * t226;
t107 = -qJ(4) * t226 + t141;
t151 = qJD(2) * qJ(3);
t89 = -t107 - t151;
t83 = qJD(2) * pkin(4) - t89;
t42 = -pkin(5) * t101 + qJD(6) + t83;
t299 = t174 * mrSges(6,3) + t189 * mrSges(7,3) - t42 * t185 - t83 * t187 + t256 * t275 + t257 * t276;
t296 = t301 * t157;
t295 = t301 * t155;
t117 = t141 + t151;
t120 = mrSges(4,2) * t226 + qJD(2) * mrSges(4,3);
t138 = Ifges(4,5) * t227;
t259 = t131 / 0.2e1;
t262 = -t102 / 0.2e1;
t263 = t101 / 0.2e1;
t271 = t157 * t302 - t295;
t272 = -t155 * t300 + t296;
t240 = Ifges(7,6) * t155;
t241 = Ifges(6,6) * t155;
t242 = Ifges(7,5) * t157;
t243 = Ifges(6,5) * t157;
t273 = -t241 + t243 - t240 + t242;
t277 = qJD(2) / 0.2e1;
t279 = qJD(1) / 0.2e1;
t280 = -qJD(1) / 0.2e1;
t294 = -(m(4) * t117 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t226 + t120) * pkin(7) + t75 * mrSges(5,2) + t97 * mrSges(4,1) + Ifges(4,6) * t277 - Ifges(4,3) * t226 / 0.2e1 + t138 / 0.2e1 + (Ifges(3,4) * t156 + t158 * Ifges(3,2)) * t280 + (-t158 * Ifges(5,1) - Ifges(5,4) * t156) * t279 - t117 * mrSges(4,2) - t89 * mrSges(5,3) + t272 * t263 + t271 * t262 + t273 * t259 - (Ifges(3,6) + Ifges(5,5)) * qJD(2) / 0.2e1 - t299;
t104 = t140 - t133;
t293 = -qJD(3) - t104;
t292 = -mrSges(3,1) - mrSges(4,1);
t217 = qJD(1) * qJD(2);
t204 = t158 * t217;
t216 = qJD(2) * qJD(5);
t220 = qJD(5) * t158;
t66 = -t157 * t216 + (t155 * t220 + t156 * t224) * qJD(1);
t207 = t157 * t220;
t225 = qJD(2) * t156;
t167 = t155 * t225 - t207;
t67 = -qJD(1) * t167 + t155 * t216;
t289 = t204 * t290 + t300 * t67 + t301 * t66;
t288 = t204 * t291 + t301 * t67 + t302 * t66;
t232 = qJ(6) - t150;
t202 = qJD(5) * t232;
t208 = t155 * t227;
t135 = qJ(3) * t226;
t168 = pkin(4) * t158 + t150 * t156;
t54 = qJD(1) * t168 + t135;
t25 = t107 * t157 + t155 * t54;
t287 = qJ(6) * t208 - qJD(6) * t157 + t155 * t202 - t25;
t233 = t156 * t157;
t172 = pkin(5) * t158 - qJ(6) * t233;
t24 = -t107 * t155 + t157 * t54;
t286 = -qJD(1) * t172 + qJD(6) * t155 + t157 * t202 - t24;
t222 = qJD(5) * t155;
t285 = -t293 + (-t208 - t222) * pkin(5);
t284 = t157 * t300 + t295;
t283 = t155 * t302 + t296;
t282 = -t66 / 0.2e1;
t281 = -t67 / 0.2e1;
t252 = pkin(7) - qJ(4);
t121 = t252 * t156;
t103 = t157 * t121;
t113 = -pkin(2) * t158 - qJ(3) * t156 - pkin(1);
t98 = pkin(3) * t158 - t113;
t72 = t195 + t98;
t31 = t155 * t72 + t103;
t221 = qJD(5) * t157;
t164 = t168 * qJD(2);
t144 = t156 * qJD(3);
t230 = qJ(3) * t204 + qJD(1) * t144;
t29 = qJD(1) * t164 + t230;
t223 = qJD(2) * t158;
t87 = -qJD(4) * t156 + t223 * t252;
t77 = t87 * qJD(1);
t3 = t155 * t29 + t157 * t77 + t221 * t48 - t222 * t73;
t4 = -qJD(5) * t13 - t155 * t77 + t157 * t29;
t270 = -t155 * t4 + t157 * t3;
t149 = qJD(2) * qJD(3);
t170 = pkin(7) * t225 + qJD(4) * t158;
t205 = t156 * t217;
t64 = -qJ(4) * t205 + qJD(1) * t170 - t149;
t112 = -qJD(2) * pkin(2) + t218;
t139 = Ifges(3,4) * t226;
t212 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t213 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t214 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t209 = t159 * qJD(2);
t80 = t209 + t206;
t267 = (m(4) * t112 + (mrSges(4,2) + mrSges(3,3)) * t227 + t292 * qJD(2)) * pkin(7) - t213 * t101 - t214 * t102 + t212 * t131 + t112 * mrSges(4,2) + t12 * mrSges(6,1) + t75 * mrSges(5,1) + t9 * mrSges(7,1) + (-Ifges(5,4) * t158 - t156 * Ifges(5,2)) * t280 + (t156 * Ifges(4,1) - Ifges(4,5) * t158) * t279 + Ifges(3,1) * t227 / 0.2e1 + t139 / 0.2e1 - t11 * mrSges(7,2) - t13 * mrSges(6,2) - t80 * mrSges(5,3) - t97 * mrSges(4,3) + t290 * t263 + t291 * t262 + (Ifges(7,3) + Ifges(6,3)) * t259 + (Ifges(5,6) + Ifges(4,4) + Ifges(3,5)) * t277;
t266 = pkin(1) * mrSges(3,1);
t265 = pkin(1) * mrSges(3,2);
t264 = -t101 / 0.2e1;
t261 = t102 / 0.2e1;
t260 = -t131 / 0.2e1;
t154 = qJ(3) + pkin(4);
t43 = mrSges(7,1) * t204 - mrSges(7,3) * t66;
t44 = mrSges(6,1) * t204 - mrSges(6,3) * t66;
t251 = t43 + t44;
t45 = -mrSges(7,2) * t204 + mrSges(7,3) * t67;
t46 = -mrSges(6,2) * t204 + mrSges(6,3) * t67;
t250 = t45 + t46;
t68 = -mrSges(7,2) * t131 + mrSges(7,3) * t101;
t69 = -mrSges(6,2) * t131 + mrSges(6,3) * t101;
t249 = t68 + t69;
t70 = mrSges(7,1) * t131 + t102 * mrSges(7,3);
t71 = mrSges(6,1) * t131 + mrSges(6,3) * t102;
t248 = -t70 - t71;
t146 = t158 * qJ(4);
t122 = pkin(7) * t158 - t146;
t237 = t122 * t64;
t236 = -qJD(2) * mrSges(5,1) + mrSges(6,1) * t101 + mrSges(6,2) * t102 + mrSges(5,3) * t226;
t235 = qJ(6) * t158;
t234 = t155 * t156;
t229 = mrSges(5,1) * t204 + mrSges(5,2) * t205;
t228 = qJ(3) * t223 + t144;
t219 = qJD(6) * t158;
t47 = t164 + t228;
t215 = t155 * t47 + t157 * t87 + t221 * t72;
t211 = t120 - t236;
t210 = m(4) * pkin(7) + mrSges(4,2);
t20 = -t67 * mrSges(7,1) + mrSges(7,2) * t66;
t203 = -t155 * t87 + t157 * t47;
t30 = -t121 * t155 + t157 * t72;
t201 = t156 * t209;
t200 = -0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(4,5);
t199 = -Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1 - Ifges(3,6) / 0.2e1;
t198 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t1 = pkin(5) * t204 - qJ(6) * t66 + qJD(6) * t102 + t4;
t2 = qJ(6) * t67 + qJD(6) * t101 + t3;
t194 = t1 * t157 + t155 * t2;
t193 = t1 * t155 - t157 * t2;
t191 = t155 * t3 + t157 * t4;
t190 = t11 * t157 - t155 * t9;
t188 = -mrSges(6,1) * t157 + mrSges(6,2) * t155;
t186 = mrSges(7,1) * t157 - mrSges(7,2) * t155;
t173 = -t12 * t155 + t13 * t157;
t85 = -qJ(4) * t225 + t170;
t166 = t155 * t248 + t157 * t249;
t165 = -t155 * t249 + t157 * t248;
t163 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t130 = pkin(5) * t157 + t154;
t129 = Ifges(6,3) * t204;
t128 = Ifges(7,3) * t204;
t114 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t227;
t111 = -pkin(7) * t205 + t149;
t110 = t232 * t157;
t109 = t232 * t155;
t108 = (mrSges(5,1) * t156 - mrSges(5,2) * t158) * qJD(1);
t105 = (-mrSges(4,1) * t158 - mrSges(4,3) * t156) * qJD(1);
t88 = -t146 + (-pkin(5) * t155 + pkin(7)) * t158;
t86 = pkin(2) * t225 - t228;
t84 = t159 * t227 + t135;
t76 = pkin(2) * t205 - t230;
t74 = t201 + t228;
t63 = Ifges(6,5) * t66;
t62 = Ifges(7,5) * t66;
t61 = Ifges(6,6) * t67;
t60 = Ifges(7,6) * t67;
t55 = qJD(1) * t201 + t230;
t52 = -mrSges(7,1) * t101 - mrSges(7,2) * t102;
t49 = pkin(5) * t167 - t85;
t26 = t155 * t235 + t31;
t23 = -pkin(5) * t67 - t64;
t22 = pkin(5) * t156 + t157 * t235 + t30;
t21 = -mrSges(6,1) * t67 + mrSges(6,2) * t66;
t8 = -qJD(5) * t31 + t203;
t7 = -t121 * t222 + t215;
t6 = qJ(6) * t207 + (-qJ(6) * t225 - qJD(5) * t121 + t219) * t155 + t215;
t5 = t157 * t219 + t172 * qJD(2) + (-t103 + (-t72 - t235) * t155) * qJD(5) + t203;
t14 = [t98 * t229 + t122 * t21 + t74 * t108 + t87 * t114 + t86 * t105 + t88 * t20 + t6 * t68 + t7 * t69 + t5 * t70 + t8 * t71 + t22 * t43 + t30 * t44 + t26 * t45 + t31 * t46 + t49 * t52 + t236 * t85 + m(6) * (t12 * t8 + t13 * t7 + t3 * t31 + t30 * t4 - t83 * t85 - t237) + m(5) * (t121 * t77 + t55 * t98 + t74 * t75 + t80 * t87 + t85 * t89 - t237) + m(7) * (t1 * t22 + t11 * t6 + t2 * t26 + t23 * t88 + t42 * t49 + t5 * t9) + m(4) * (t113 * t76 + t86 * t97) + (t128 / 0.2e1 + t129 / 0.2e1 + t62 / 0.2e1 + t63 / 0.2e1 + t60 / 0.2e1 + t61 / 0.2e1 - t77 * mrSges(5,3) - t76 * mrSges(4,3) + t55 * mrSges(5,1) - t213 * t67 + t214 * t66 + ((mrSges(4,1) * t113 + mrSges(5,3) * t122 + t156 * t200 - 0.2e1 * t266) * qJD(1) + t199 * qJD(2) + t294) * qJD(2) + t163) * t156 + (-t23 * t185 - t76 * mrSges(4,1) - t55 * mrSges(5,2) + (mrSges(5,3) + t187) * t64 + t210 * t111 + t194 * mrSges(7,3) + t191 * mrSges(6,3) + (t173 * mrSges(6,3) + t190 * mrSges(7,3) - t42 * t186 + t83 * t188 + t284 * t263 + t283 * t262 + (t155 * t291 + t157 * t290) * t259 + t276 * t157 / 0.2e1) * qJD(5) + ((-t113 * mrSges(4,3) - t121 * mrSges(5,3) - 0.2e1 * t265 + (-t242 / 0.2e1 + t240 / 0.2e1 - t243 / 0.2e1 + t241 / 0.2e1 - t200) * t158 + (0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(5,1) + t210 * pkin(7) + t212) * t156) * qJD(1) + t198 * qJD(2) + t267) * qJD(2) + t271 * t282 + t272 * t281 + (qJD(5) * t275 + t289) * t257 + t288 * t256) * t158; ((-m(6) * t174 - t155 * t69 - t157 * t71) * t150 + t272 * t264 + t271 * t261 + t273 * t260 + t299) * qJD(5) - t270 * mrSges(6,3) + (-t155 * t44 + t157 * t46) * t150 + (-t64 * qJ(3) - t107 * t80 + t159 * t77 + t293 * t89 - t75 * t84) * m(5) + ((-t139 / 0.2e1 + (-pkin(2) * mrSges(4,2) - t159 * mrSges(5,3) + t213 * t157 - t214 * t155 + (-m(4) * pkin(2) + t292) * pkin(7) + t198) * qJD(2) + ((-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t158 + t265) * qJD(1) - t267) * t158 + (-t138 / 0.2e1 + ((Ifges(3,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t156 + t266) * qJD(1) + (pkin(7) * mrSges(3,2) + (mrSges(5,3) - mrSges(4,2)) * qJ(3) + t199) * qJD(2) + (Ifges(3,2) / 0.2e1 - Ifges(4,1) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t226 - t294) * t156) * qJD(1) - t236 * t104 + t211 * qJD(3) + t154 * t21 + t130 * t20 - t84 * t108 + t109 * t43 - t110 * t45 + t111 * mrSges(4,3) - t107 * t114 + t77 * mrSges(5,2) - t25 * t69 - t24 * t71 + t283 * t282 + t284 * t281 + (-m(4) * t97 - t105) * (pkin(2) * t227 - t135) + m(4) * (qJ(3) * t111 + qJD(3) * t117) + (-mrSges(5,1) + t188) * t64 + t193 * mrSges(7,3) + t23 * t186 - m(6) * (-t104 * t83 + t12 * t24 + t13 * t25) + m(6) * (qJD(3) * t83 + t150 * t270 - t154 * t64) + t285 * t52 + t286 * t70 + (t1 * t109 + t11 * t287 - t110 * t2 + t130 * t23 + t285 * t42 + t286 * t9) * m(7) + t287 * t68 - t288 * t155 / 0.2e1 + t289 * t256; t250 * t157 - t251 * t155 + t165 * qJD(5) + (-t52 - t211) * qJD(2) + ((-mrSges(5,3) + t210) * t223 + (t105 - t108 + t165) * t156) * qJD(1) - m(4) * (qJD(2) * t117 - t227 * t97) + (-qJD(2) * t42 - t131 * t189 - t193) * m(7) + (-qJD(2) * t83 - t131 * t174 + t270) * m(6) + (qJD(2) * t89 - t227 * t75 + t77) * m(5); t251 * t157 + t250 * t155 + t166 * qJD(5) + m(5) * t55 + m(6) * (qJD(5) * t173 + t191) + m(7) * (qJD(5) * t190 + t194) + ((t52 - t236) * t158 + (t114 + t166) * t156 - m(5) * (-t156 * t80 + t158 * t89) - m(6) * (t12 * t234 - t13 * t233 - t158 * t83) - m(7) * (-t11 * t233 - t158 * t42 + t234 * t9)) * qJD(1) + t229; (t102 * t52 + t43) * pkin(5) + t128 + t129 + t62 + t63 + t60 + t61 + (t101 * t9 - t102 * t11) * mrSges(7,3) + (t101 * t12 - t102 * t13) * mrSges(6,3) - t83 * (-mrSges(6,1) * t102 + mrSges(6,2) * t101) - t42 * (-mrSges(7,1) * t102 + mrSges(7,2) * t101) - t10 * t68 - t12 * t69 + t11 * t70 + t13 * t71 + t163 + (-(t10 - t9) * t11 + (t102 * t42 + t1) * pkin(5)) * m(7) + t276 * t262 + (t101 * t302 + t297) * t261 + (t101 * t291 + t102 * t290) * t260 + (t102 * t300 + t275 + t298) * t264; -t101 * t68 - t102 * t70 + 0.2e1 * (t23 / 0.2e1 + t11 * t264 + t9 * t262) * m(7) + t20;];
tauc  = t14(:);
