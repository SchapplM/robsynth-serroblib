% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2018-11-23 16:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:14:41
% EndTime: 2018-11-23 16:14:47
% DurationCPUTime: 6.49s
% Computational Cost: add. (3347->494), mult. (7424->598), div. (0->0), fcn. (3880->4), ass. (0->229)
t327 = Ifges(7,4) + Ifges(6,5);
t313 = Ifges(6,4) + Ifges(5,5);
t328 = -Ifges(7,5) + t313;
t323 = Ifges(7,2) + Ifges(6,3);
t326 = Ifges(5,6) - Ifges(6,6);
t325 = Ifges(6,6) - Ifges(7,6);
t152 = -pkin(1) - pkin(7);
t137 = qJD(1) * t152 + qJD(2);
t150 = cos(qJ(3));
t107 = -qJD(3) * pkin(3) - t137 * t150;
t147 = sin(qJ(4));
t149 = cos(qJ(4));
t148 = sin(qJ(3));
t123 = pkin(3) * t148 - pkin(8) * t150 + qJ(2);
t104 = t123 * qJD(1);
t122 = t148 * t137;
t106 = qJD(3) * pkin(8) + t122;
t47 = t149 * t104 - t147 * t106;
t48 = t147 * t104 + t149 * t106;
t168 = t147 * t48 + t149 * t47;
t233 = qJD(1) * t148;
t141 = qJD(4) + t233;
t298 = qJD(5) - t47;
t32 = -pkin(4) * t141 + t298;
t134 = t141 * qJ(5);
t33 = t134 + t48;
t170 = t147 * t33 - t149 * t32;
t283 = pkin(4) + pkin(5);
t232 = qJD(1) * t150;
t119 = qJD(3) * t147 + t149 * t232;
t26 = qJ(6) * t119 + t47;
t299 = qJD(5) - t26;
t11 = -t141 * t283 + t299;
t222 = t149 * qJD(3);
t118 = t147 * t232 - t222;
t27 = qJ(6) * t118 + t48;
t14 = t134 + t27;
t172 = t11 * t149 - t14 * t147;
t254 = Ifges(5,4) * t149;
t186 = -Ifges(5,2) * t147 + t254;
t193 = mrSges(7,1) * t147 - mrSges(7,2) * t149;
t195 = mrSges(6,1) * t147 - mrSges(6,3) * t149;
t197 = mrSges(5,1) * t147 + mrSges(5,2) * t149;
t160 = qJ(5) * t119 - t107;
t21 = -t118 * t283 + qJD(6) + t160;
t273 = t149 / 0.2e1;
t275 = t147 / 0.2e1;
t276 = -t147 / 0.2e1;
t277 = t141 / 0.2e1;
t278 = -t141 / 0.2e1;
t279 = t119 / 0.2e1;
t281 = t118 / 0.2e1;
t282 = -t118 / 0.2e1;
t255 = Ifges(5,4) * t147;
t316 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t317 = t327 * t147;
t293 = t149 * t316 - t255 + t317;
t113 = Ifges(5,4) * t118;
t304 = t327 * t118;
t295 = t316 * t119 + t141 * t328 - t113 + t304;
t318 = t327 * t149;
t296 = t147 * t323 + t318;
t307 = t327 * t119;
t300 = t323 * t118 + t141 * t325 + t307;
t35 = pkin(4) * t118 - t160;
t256 = Ifges(5,4) * t119;
t43 = -t118 * Ifges(5,2) + t141 * Ifges(5,6) + t256;
t324 = t273 * t295 + t275 * t300 + t276 * t43 + t107 * t197 + (Ifges(7,5) * t149 + Ifges(7,6) * t147) * t278 + t186 * t282 + t35 * t195 + t296 * t281 + (-t147 * t326 + t149 * t313) * t277 + t293 * t279 - t21 * t193 - t170 * mrSges(6,2) - t168 * mrSges(5,3) - t172 * mrSges(7,3);
t321 = qJD(5) * t147 + t122;
t239 = qJ(5) * t149;
t165 = t147 * t283 - t239;
t320 = t152 - t165;
t243 = Ifges(4,5) * qJD(3);
t258 = Ifges(4,4) * t148;
t315 = qJD(1) / 0.2e1;
t319 = t243 / 0.2e1 + (t150 * Ifges(4,1) - t258) * t315 + t324;
t314 = -qJD(3) / 0.2e1;
t230 = qJD(3) * t150;
t206 = qJD(1) * t230;
t227 = qJD(4) * t150;
t75 = qJD(4) * t222 + (-t147 * t227 - t148 * t222) * qJD(1);
t213 = t147 * t233;
t76 = -qJD(3) * t213 + qJD(4) * t119;
t312 = t206 * t325 + t323 * t76 + t327 * t75;
t311 = -t141 * t165 + t321;
t223 = qJD(6) * t149;
t229 = qJD(4) * t147;
t267 = pkin(8) - qJ(6);
t202 = pkin(3) * t150 + pkin(8) * t148;
t121 = t202 * qJD(1);
t234 = t149 * t150;
t60 = t147 * t121 + t137 * t234;
t50 = qJ(5) * t232 + t60;
t310 = qJ(6) * t213 - t229 * t267 - t223 - t50;
t173 = pkin(4) * t147 - t239;
t309 = t141 * t173 - t321;
t131 = t267 * t149;
t238 = qJ(6) * t149;
t236 = t147 * t150;
t59 = t121 * t149 - t137 * t236;
t308 = qJD(4) * t131 - qJD(6) * t147 - (t148 * t238 - t150 * t283) * qJD(1) + t59;
t306 = -t149 * t323 + t317;
t305 = t147 * t313 + t326 * t149;
t303 = qJD(1) * (qJ(2) * (m(3) + m(4)) + mrSges(3,3));
t302 = (-Ifges(5,4) + t327) * t76 + t316 * t75 + t328 * t206;
t301 = t147 * t316 + t254 - t318;
t280 = -t119 / 0.2e1;
t207 = Ifges(4,6) * t314;
t297 = qJ(5) * t230 + t148 * qJD(5);
t211 = t137 * t230;
t228 = qJD(4) * t149;
t117 = qJD(3) * t202 + qJD(2);
t95 = t117 * qJD(1);
t10 = -t104 * t229 - t106 * t228 - t147 * t211 + t149 * t95;
t9 = t104 * t228 - t106 * t229 + t147 * t95 + t149 * t211;
t199 = -t10 * t147 + t149 * t9;
t4 = qJ(5) * t206 + t141 * qJD(5) + t9;
t7 = -pkin(4) * t206 - t10;
t200 = t147 * t7 + t149 * t4;
t240 = qJ(5) * t147;
t291 = -t149 * t283 - t240;
t290 = -m(5) * t168 - m(6) * t170;
t216 = -Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t203 = Ifges(7,6) / 0.2e1 - t216;
t218 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t205 = Ifges(7,5) / 0.2e1 - t218;
t217 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t257 = Ifges(4,4) * t150;
t289 = t203 * t118 + t205 * t119 - (Ifges(7,3) / 0.2e1 + t217) * t141 - t14 * mrSges(7,2) - t33 * mrSges(6,3) - t47 * mrSges(5,1) - t207 + (-Ifges(4,2) * t148 + t257) * t315 - Ifges(7,5) * t280 - Ifges(6,6) * t281 + t11 * mrSges(7,1) + t32 * mrSges(6,1) + t48 * mrSges(5,2) - (Ifges(7,6) + Ifges(5,6)) * t282 - t313 * t279 - (Ifges(7,3) + Ifges(5,3) + Ifges(6,2)) * t277;
t287 = t75 / 0.2e1;
t286 = -t76 / 0.2e1;
t285 = t76 / 0.2e1;
t274 = -t149 / 0.2e1;
t269 = t75 * mrSges(7,3);
t53 = -mrSges(6,2) * t76 + mrSges(6,3) * t206;
t58 = -mrSges(5,2) * t206 - mrSges(5,3) * t76;
t266 = t53 + t58;
t55 = mrSges(5,1) * t206 - mrSges(5,3) * t75;
t70 = t75 * mrSges(6,2);
t56 = -mrSges(6,1) * t206 + t70;
t265 = -t55 + t56;
t62 = mrSges(6,1) * t118 - mrSges(6,3) * t119;
t63 = -mrSges(7,1) * t118 + mrSges(7,2) * t119;
t264 = t62 - t63;
t260 = mrSges(5,3) * t118;
t78 = -mrSges(5,2) * t141 - t260;
t82 = -mrSges(6,2) * t118 + mrSges(6,3) * t141;
t263 = t78 + t82;
t259 = mrSges(5,3) * t119;
t80 = mrSges(5,1) * t141 - t259;
t81 = -mrSges(6,1) * t141 + mrSges(6,2) * t119;
t262 = -t80 + t81;
t77 = mrSges(7,2) * t141 + mrSges(7,3) * t118;
t261 = t82 + t77;
t249 = qJ(2) * mrSges(4,1);
t248 = qJ(2) * mrSges(4,2);
t244 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t118 + mrSges(5,2) * t119 + mrSges(4,3) * t232;
t241 = qJ(5) * t118;
t237 = qJD(3) * mrSges(4,2);
t235 = t148 * t152;
t84 = t147 * t123 + t149 * t235;
t231 = qJD(3) * t148;
t226 = qJD(4) * t152;
t224 = qJD(5) * t149;
t221 = qJD(1) * qJD(2);
t220 = t78 + t261;
t79 = -mrSges(7,1) * t141 - mrSges(7,3) * t119;
t219 = t79 + t262;
t210 = t152 * t230;
t215 = t147 * t117 + t123 * t228 + t149 * t210;
t209 = t148 * t226;
t214 = t123 * t229 + t147 * t210 + t149 * t209;
t67 = t148 * qJ(5) + t84;
t212 = t137 * t231;
t208 = -t243 / 0.2e1;
t24 = -t76 * mrSges(7,1) + t75 * mrSges(7,2);
t135 = t147 * t235;
t83 = t123 * t149 - t135;
t1 = qJ(6) * t76 + qJD(6) * t118 + t4;
t2 = -qJ(6) * t75 - qJD(6) * t119 - t206 * t283 - t10;
t201 = -t1 * t149 - t147 * t2;
t198 = mrSges(5,1) * t149 - mrSges(5,2) * t147;
t196 = mrSges(6,1) * t149 + mrSges(6,3) * t147;
t194 = mrSges(7,1) * t149 + mrSges(7,2) * t147;
t185 = Ifges(5,2) * t149 + t255;
t175 = Ifges(7,5) * t147 - Ifges(7,6) * t149;
t174 = pkin(4) * t149 + t240;
t171 = t11 * t147 + t14 * t149;
t169 = t147 * t32 + t149 * t33;
t167 = t147 * t47 - t149 * t48;
t29 = t117 * t149 - t214;
t166 = -t152 + t173;
t28 = -t147 * t209 + t215;
t159 = qJ(5) * t75 + qJD(5) * t119 - t212;
t158 = -t147 * t220 + t149 * t219;
t157 = t10 * mrSges(5,1) - t7 * mrSges(6,1) - t2 * mrSges(7,1) - t9 * mrSges(5,2) + t1 * mrSges(7,2) + t4 * mrSges(6,3);
t140 = Ifges(6,2) * t206;
t139 = Ifges(5,3) * t206;
t130 = t267 * t147;
t128 = -mrSges(4,3) * t233 - t237;
t124 = -pkin(3) - t174;
t120 = (mrSges(4,1) * t148 + mrSges(4,2) * t150) * qJD(1);
t110 = pkin(3) - t291;
t87 = t166 * t150;
t74 = Ifges(6,4) * t75;
t73 = Ifges(5,5) * t75;
t72 = Ifges(5,6) * t76;
t71 = Ifges(6,6) * t76;
t68 = -pkin(4) * t148 - t83;
t66 = t320 * t150;
t61 = pkin(4) * t119 + t241;
t57 = mrSges(7,2) * t206 + mrSges(7,3) * t76;
t54 = -mrSges(7,1) * t206 - t269;
t52 = qJ(6) * t236 + t67;
t51 = -pkin(4) * t232 - t59;
t37 = t135 + (-qJ(6) * t150 - t123) * t149 - t283 * t148;
t36 = -t119 * t283 - t241;
t31 = (qJD(4) * t174 - t224) * t150 - t166 * t231;
t25 = mrSges(5,1) * t76 + mrSges(5,2) * t75;
t23 = mrSges(6,1) * t76 - mrSges(6,3) * t75;
t22 = -pkin(4) * t230 - t29;
t17 = t75 * Ifges(5,4) - t76 * Ifges(5,2) + Ifges(5,6) * t206;
t13 = t28 + t297;
t12 = (qJD(4) * t291 + t224) * t150 - t320 * t231;
t8 = pkin(4) * t76 - t159;
t6 = t227 * t238 + (qJD(6) * t150 + (-qJ(6) * qJD(3) - t226) * t148) * t147 + t215 + t297;
t5 = (qJ(6) * t231 - t117) * t149 + (qJ(6) * t229 - qJD(3) * t283 - t223) * t150 + t214;
t3 = -t283 * t76 + t159;
t15 = [t29 * t80 + t22 * t81 + t13 * t82 + t83 * t55 + t84 * t58 + t87 * t23 + t6 * t77 + t28 * t78 + t5 * t79 + t31 * t62 + t12 * t63 + t66 * t24 + t67 * t53 + t68 * t56 + t37 * t54 + t52 * t57 + m(6) * (t13 * t33 + t22 * t32 + t31 * t35 + t4 * t67 + t68 * t7 + t8 * t87) + m(5) * (t10 * t83 + t28 * t48 + t29 * t47 + t84 * t9) + m(7) * (t1 * t52 + t11 * t5 + t12 * t21 + t14 * t6 + t2 * t37 + t3 * t66) + (t120 + 0.2e1 * t303) * qJD(2) + (t139 / 0.2e1 + t140 / 0.2e1 + t74 / 0.2e1 + t71 / 0.2e1 - t72 / 0.2e1 + t73 / 0.2e1 + mrSges(4,1) * t221 + (-Ifges(7,6) + t216) * t76 + (-Ifges(7,5) + t218) * t75 + ((-0.2e1 * t248 + 0.3e1 / 0.2e1 * t258) * qJD(1) + (m(5) * t107 + t244) * t152 + t208 - t319) * qJD(3) + t157) * t148 + (-t152 * t25 + t186 * t286 + t8 * t195 - t3 * t193 + mrSges(4,2) * t221 + (t1 * t147 - t149 * t2) * mrSges(7,3) + (-t10 * t149 - t147 * t9) * mrSges(5,3) + (-t147 * t4 + t149 * t7) * mrSges(6,2) + ((-m(5) * t152 + t197) * t122 + (-0.3e1 / 0.2e1 * t257 + 0.2e1 * t249 - t205 * t234 - t203 * t236 + (Ifges(7,3) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1) + t217) * t148) * qJD(1) + t152 * t128 + t207 - t289) * qJD(3) + t296 * t285 + t312 * t275 + t293 * t287 + t17 * t276 + t302 * t273 + (-t169 * mrSges(6,2) + t167 * mrSges(5,3) + t171 * mrSges(7,3) + t107 * t198 + t175 * t277 + t185 * t281 - t21 * t194 + t35 * t196 + t300 * t273 + t43 * t274 + t295 * t276 + t278 * t305 + t280 * t301 + t282 * t306) * qJD(4)) * t150; (-t23 + t24 - t25 - m(6) * t8 + m(7) * t3 + (-m(5) * t167 + m(6) * t169 + m(7) * t171 + t147 * t219 + t149 * t220 + t128) * qJD(3)) * t150 + ((t57 + t266) * t149 + (t54 + t265) * t147 + (t244 + t264) * qJD(3) + t158 * qJD(4) + m(5) * (qJD(3) * t107 - t228 * t47 - t229 * t48 + t199 - t211) + m(6) * (qJD(3) * t35 + t228 * t32 - t229 * t33 + t200) + m(7) * (-qJD(3) * t21 + t11 * t228 - t14 * t229 - t201)) * t148 + (m(7) * t172 - t120 + t158 + t290 - t303) * qJD(1); ((t175 * t314 + (t257 / 0.2e1 - t249) * qJD(1) + t207 + t305 * qJD(3) / 0.2e1 + t289) * t150 + ((t248 - t258 / 0.2e1 + (-Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t150) * qJD(1) + t208 + t319) * t148) * qJD(1) + ((-t147 * t263 + t149 * t262 + t290) * pkin(8) + t324) * qJD(4) + (-pkin(3) * t212 + pkin(8) * t199 - t107 * t122 - t47 * t59 - t48 * t60) * m(5) + t131 * t57 + t124 * t23 + t130 * t54 + t110 * t24 - t59 * t80 - t51 * t81 - t50 * t82 - t60 * t78 + (t147 * t265 + t149 * t266) * pkin(8) + ((-t128 - t237) * t150 + ((-mrSges(4,1) - t198) * qJD(3) - t244) * t148) * t137 + t3 * t194 - pkin(3) * t25 - t8 * t196 + t199 * mrSges(5,3) + t200 * mrSges(6,2) + t201 * mrSges(7,3) + t308 * t79 + (pkin(8) * t200 + t124 * t8 + t309 * t35 - t32 * t51 - t33 * t50) * m(6) + t309 * t62 + t310 * t77 + t311 * t63 + (t1 * t131 + t11 * t308 + t110 * t3 + t130 * t2 + t14 * t310 + t21 * t311) * m(7) + t312 * t274 + t302 * t275 + t306 * t285 + t301 * t287 + t185 * t286 + t17 * t273; (-t118 * t316 - t256 + t300 + t307) * t280 + (t119 * t323 - t304) * t282 + (-t313 * t118 - t119 * t326) * t278 - t283 * t54 + (qJ(5) * t1 - t11 * t27 + t14 * t299 - t2 * t283 - t21 * t36) * m(7) + t139 + t140 - t107 * (mrSges(5,1) * t119 - mrSges(5,2) * t118) - t35 * (mrSges(6,1) * t119 + mrSges(6,3) * t118) - t21 * (-mrSges(7,1) * t119 - mrSges(7,2) * t118) - Ifges(7,5) * t75 - Ifges(7,6) * t76 - t26 * t77 - t27 * t79 + (-pkin(4) * t7 + qJ(5) * t4 + t298 * t33 - t32 * t48 - t35 * t61) * m(6) + (-Ifges(5,2) * t119 - t113 + t295) * t281 + t74 + t71 - t72 + t73 + (-t260 - t263) * t47 + t261 * qJD(5) + (t259 - t262) * t48 + (t53 + t57) * qJ(5) - t36 * t63 - pkin(4) * t56 - t61 * t62 + (-t11 * t118 - t119 * t14) * mrSges(7,3) + t157 + (t118 * t32 + t119 * t33) * mrSges(6,2) + Ifges(7,3) * t206 + (-Ifges(7,5) * t118 + Ifges(7,6) * t119) * t277 + t43 * t279; -t269 + t70 - t261 * t141 + t264 * t119 + (-mrSges(6,1) - mrSges(7,1)) * t206 + (-t119 * t21 - t14 * t141 + t2) * m(7) + (t119 * t35 - t141 * t33 + t7) * m(6); -t118 * t77 + t119 * t79 + 0.2e1 * (t3 / 0.2e1 + t11 * t279 + t14 * t282) * m(7) + t24;];
tauc  = t15(:);
