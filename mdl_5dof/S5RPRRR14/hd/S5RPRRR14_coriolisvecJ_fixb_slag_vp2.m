% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:54
% EndTime: 2019-12-31 19:17:18
% DurationCPUTime: 10.37s
% Computational Cost: add. (12802->547), mult. (42768->809), div. (0->0), fcn. (35902->12), ass. (0->257)
t186 = sin(pkin(11));
t188 = sin(pkin(5));
t189 = cos(pkin(11));
t191 = cos(pkin(5));
t194 = sin(qJ(3));
t190 = cos(pkin(6));
t197 = cos(qJ(3));
t257 = t190 * t197;
t187 = sin(pkin(6));
t261 = t187 * t197;
t202 = t191 * t261 + t188 * (-t186 * t194 + t189 * t257);
t144 = t202 * qJD(1);
t142 = qJD(4) - t144;
t258 = t190 * t194;
t206 = (-t186 * t258 + t189 * t197) * t188;
t161 = qJD(1) * t206;
t249 = qJD(3) * t197;
t334 = t187 * t249 - t161;
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t260 = t188 * t189;
t182 = qJ(2) * t260;
t298 = pkin(1) * t191;
t243 = qJD(1) * t298;
t168 = qJD(1) * t182 + t186 * t243;
t208 = (t187 * t191 + t190 * t260) * pkin(8);
t137 = qJD(1) * t208 + t168;
t181 = t189 * t243;
t264 = t186 * t188;
t205 = pkin(2) * t191 + (-pkin(8) * t190 - qJ(2)) * t264;
t143 = qJD(1) * t205 + t181;
t162 = (-pkin(8) * t186 * t187 - pkin(2) * t189 - pkin(1)) * t188;
t156 = qJD(1) * t162 + qJD(2);
t218 = t143 * t190 + t156 * t187;
t84 = t137 * t197 + t194 * t218;
t333 = -pkin(9) * qJD(5) * t196 - t84 + t142 * (pkin(4) * t193 - pkin(10) * t196);
t262 = t187 * t194;
t151 = t191 * t262 + (t186 * t197 + t189 * t258) * t188;
t147 = t151 * qJD(1);
t177 = -pkin(4) * t196 - pkin(10) * t193 - pkin(3);
t248 = qJD(4) * t193;
t112 = pkin(3) * t147 - pkin(9) * t144;
t83 = -t194 * t137 + t218 * t197;
t55 = t193 * t112 + t196 * t83;
t332 = pkin(9) * t248 + pkin(10) * t147 - qJD(5) * t177 + t55;
t172 = -t196 * t190 + t193 * t262;
t252 = qJD(1) * t188;
t240 = t186 * t252;
t236 = t187 * t240;
t318 = -qJD(4) * t172 - t193 * t236 + t334 * t196;
t207 = (t186 * t257 + t189 * t194) * t188;
t160 = qJD(1) * t207;
t250 = qJD(3) * t194;
t331 = t187 * t250 - t160;
t192 = sin(qJ(5));
t195 = cos(qJ(5));
t224 = Ifges(6,5) * t195 - Ifges(6,6) * t192;
t283 = Ifges(6,4) * t195;
t226 = -Ifges(6,2) * t192 + t283;
t284 = Ifges(6,4) * t192;
t228 = Ifges(6,1) * t195 - t284;
t229 = mrSges(6,1) * t192 + mrSges(6,2) * t195;
t117 = -t143 * t187 + t190 * t156;
t68 = -pkin(3) * t144 - pkin(9) * t147 + t117;
t169 = -t187 * t260 + t190 * t191;
t163 = qJD(1) * t169 + qJD(3);
t70 = t163 * pkin(9) + t84;
t32 = t193 * t68 + t196 * t70;
t29 = pkin(10) * t142 + t32;
t120 = -t147 * t193 + t163 * t196;
t121 = t147 * t196 + t163 * t193;
t69 = -t163 * pkin(3) - t83;
t38 = -t120 * pkin(4) - t121 * pkin(10) + t69;
t8 = -t192 * t29 + t195 * t38;
t9 = t192 * t38 + t195 * t29;
t231 = t192 * t9 + t195 * t8;
t31 = -t193 * t70 + t196 * t68;
t28 = -pkin(4) * t142 - t31;
t300 = t195 / 0.2e1;
t301 = -t192 / 0.2e1;
t119 = qJD(5) - t120;
t302 = t119 / 0.2e1;
t92 = t121 * t195 + t142 * t192;
t305 = t92 / 0.2e1;
t91 = -t121 * t192 + t142 * t195;
t307 = t91 / 0.2e1;
t299 = Ifges(6,4) * t92;
t36 = Ifges(6,2) * t91 + Ifges(6,6) * t119 + t299;
t88 = Ifges(6,4) * t91;
t37 = Ifges(6,1) * t92 + Ifges(6,5) * t119 + t88;
t330 = -t231 * mrSges(6,3) + t224 * t302 + t226 * t307 + t228 * t305 + t229 * t28 + t300 * t37 + t301 * t36;
t329 = Ifges(5,5) / 0.2e1;
t145 = t202 * qJD(3);
t135 = qJD(1) * t145;
t89 = qJD(4) * t120 + t135 * t196;
t328 = t89 / 0.2e1;
t90 = qJD(4) * t121 + t135 * t193;
t327 = -t90 / 0.2e1;
t326 = t333 * t192 - t332 * t195;
t325 = t332 * t192 + t333 * t195;
t253 = t186 * t298 + t182;
t148 = t208 + t253;
t184 = t189 * t298;
t152 = t184 + t205;
t217 = t152 * t190 + t162 * t187;
t95 = -t194 * t148 + t217 * t197;
t273 = t142 * Ifges(5,5);
t118 = Ifges(5,4) * t120;
t278 = t121 * Ifges(5,1);
t61 = t118 + t273 + t278;
t209 = -t61 / 0.2e1 - t273 / 0.2e1 - t69 * mrSges(5,2);
t324 = t31 * mrSges(5,3) + t209 - t330;
t125 = t151 * t193 - t169 * t196;
t323 = -t125 / 0.2e1;
t322 = t169 / 0.2e1;
t270 = t147 * mrSges(4,3);
t265 = mrSges(4,1) * t163 + mrSges(5,1) * t120 - mrSges(5,2) * t121 - t270;
t173 = t190 * t193 + t196 * t262;
t213 = -t195 * t173 + t192 * t261;
t320 = qJD(5) * t213 - t192 * t318 + t331 * t195;
t157 = -t192 * t173 - t195 * t261;
t319 = qJD(5) * t157 + t331 * t192 + t195 * t318;
t317 = qJD(4) * t173 + t334 * t193 + t196 * t236;
t146 = t151 * qJD(3);
t136 = qJD(1) * t146;
t251 = qJD(2) * t188;
t239 = t186 * t251;
t234 = qJD(1) * t239;
t215 = t187 * t234;
t102 = pkin(3) * t136 - pkin(9) * t135 + t215;
t247 = qJD(4) * t196;
t203 = qJD(2) * t206;
t65 = qJD(1) * t203 + qJD(3) * t83;
t12 = t193 * t102 + t196 * t65 + t68 * t247 - t248 * t70;
t13 = -qJD(4) * t32 + t102 * t196 - t193 * t65;
t316 = t12 * t196 - t13 * t193;
t10 = pkin(10) * t136 + t12;
t204 = qJD(2) * t207;
t66 = qJD(1) * t204 + qJD(3) * t84;
t27 = t90 * pkin(4) - t89 * pkin(10) + t66;
t1 = qJD(5) * t8 + t10 * t195 + t192 * t27;
t2 = -qJD(5) * t9 - t10 * t192 + t195 * t27;
t315 = t1 * t195 - t192 * t2;
t235 = t187 * t239;
t106 = pkin(3) * t146 - pkin(9) * t145 + t235;
t122 = -t152 * t187 + t190 * t162;
t78 = -pkin(3) * t202 - pkin(9) * t151 + t122;
t140 = t197 * t148;
t96 = t152 * t258 + t162 * t262 + t140;
t82 = pkin(9) * t169 + t96;
t286 = t193 * t78 + t196 * t82;
t72 = qJD(3) * t95 + t203;
t20 = -qJD(4) * t286 + t106 * t196 - t193 * t72;
t310 = Ifges(5,1) * t328 + Ifges(5,4) * t327 + t329 * t136;
t41 = qJD(5) * t91 + t136 * t192 + t195 * t89;
t42 = -qJD(5) * t92 + t136 * t195 - t192 * t89;
t7 = t41 * Ifges(6,1) + t42 * Ifges(6,4) + t90 * Ifges(6,5);
t314 = t7 / 0.2e1;
t313 = -t36 / 0.2e1;
t312 = t41 / 0.2e1;
t311 = t42 / 0.2e1;
t309 = t90 / 0.2e1;
t308 = -t91 / 0.2e1;
t306 = -t92 / 0.2e1;
t304 = -t118 / 0.2e1;
t303 = -t119 / 0.2e1;
t40 = Ifges(6,5) * t41;
t39 = Ifges(6,6) * t42;
t295 = t32 * mrSges(5,3);
t294 = t83 * mrSges(4,3);
t292 = t89 * Ifges(5,4);
t290 = t91 * Ifges(6,6);
t289 = t92 * Ifges(6,5);
t14 = -mrSges(6,1) * t42 + mrSges(6,2) * t41;
t63 = mrSges(5,1) * t136 - mrSges(5,3) * t89;
t288 = t14 - t63;
t53 = -mrSges(6,1) * t91 + mrSges(6,2) * t92;
t94 = mrSges(5,1) * t142 - mrSges(5,3) * t121;
t287 = t53 - t94;
t285 = Ifges(5,4) * t121;
t282 = t119 * Ifges(6,3);
t280 = t120 * Ifges(5,2);
t279 = t120 * Ifges(5,6);
t277 = t121 * Ifges(5,5);
t274 = t136 * Ifges(5,6);
t272 = t142 * Ifges(5,6);
t271 = t142 * Ifges(5,3);
t269 = t147 * Ifges(4,4);
t268 = t163 * Ifges(4,5);
t267 = t163 * Ifges(4,6);
t266 = t197 * t66;
t259 = t189 * (-mrSges(3,2) * t191 + mrSges(3,3) * t260) * qJD(1);
t256 = t192 * t196;
t255 = t195 * t196;
t254 = Ifges(4,5) * t135 - Ifges(4,6) * t136;
t5 = Ifges(6,3) * t90 + t39 + t40;
t244 = Ifges(5,5) * t89 - Ifges(5,6) * t90 + Ifges(5,3) * t136;
t233 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t230 = mrSges(6,1) * t195 - mrSges(6,2) * t192;
t227 = Ifges(6,1) * t192 + t283;
t225 = Ifges(6,2) * t195 + t284;
t223 = Ifges(6,5) * t192 + Ifges(6,6) * t195;
t34 = -pkin(10) * t202 + t286;
t126 = t151 * t196 + t169 * t193;
t81 = -t169 * pkin(3) - t95;
t49 = t125 * pkin(4) - t126 * pkin(10) + t81;
t16 = t192 * t49 + t195 * t34;
t15 = -t192 * t34 + t195 * t49;
t222 = t193 * t32 + t196 * t31;
t45 = -t193 * t82 + t196 * t78;
t54 = t112 * t196 - t193 * t83;
t104 = t126 * t195 - t192 * t202;
t103 = -t126 * t192 - t195 * t202;
t216 = -(-qJ(2) * t240 + t181) * t186 + t168 * t189;
t19 = t193 * t106 + t196 * t72 + t78 * t247 - t248 * t82;
t170 = (mrSges(3,1) * t191 - mrSges(3,3) * t264) * qJD(1);
t73 = t204 + (t194 * t217 + t140) * qJD(3);
t35 = t282 + t289 + t290;
t60 = t272 + t280 + t285;
t201 = t9 * mrSges(6,2) - t35 / 0.2e1 + t60 / 0.2e1 + t285 / 0.2e1 - t282 / 0.2e1 + t272 / 0.2e1 - t69 * mrSges(5,1) - t8 * mrSges(6,1) - t290 / 0.2e1 - t289 / 0.2e1;
t200 = t280 / 0.2e1 + t201;
t165 = pkin(9) * t255 + t177 * t192;
t164 = -pkin(9) * t256 + t177 * t195;
t141 = Ifges(4,4) * t144;
t123 = -mrSges(4,2) * t163 + mrSges(4,3) * t144;
t111 = -mrSges(4,1) * t144 + mrSges(4,2) * t147;
t109 = t144 * t255 + t147 * t192;
t108 = -t144 * t256 + t147 * t195;
t107 = mrSges(4,1) * t136 + mrSges(4,2) * t135;
t101 = -qJD(4) * t125 + t145 * t196;
t100 = qJD(4) * t126 + t145 * t193;
t99 = t147 * Ifges(4,1) + t141 + t268;
t98 = t144 * Ifges(4,2) + t267 + t269;
t93 = -mrSges(5,2) * t142 + mrSges(5,3) * t120;
t76 = pkin(4) * t121 - pkin(10) * t120;
t64 = -mrSges(5,2) * t136 - mrSges(5,3) * t90;
t59 = t271 + t277 + t279;
t58 = mrSges(6,1) * t119 - mrSges(6,3) * t92;
t57 = -mrSges(6,2) * t119 + mrSges(6,3) * t91;
t52 = mrSges(5,1) * t90 + mrSges(5,2) * t89;
t51 = qJD(5) * t103 + t101 * t195 + t146 * t192;
t50 = -qJD(5) * t104 - t101 * t192 + t146 * t195;
t47 = -pkin(4) * t147 - t54;
t43 = -t90 * Ifges(5,2) + t274 + t292;
t33 = pkin(4) * t202 - t45;
t30 = t100 * pkin(4) - t101 * pkin(10) + t73;
t26 = -mrSges(6,2) * t90 + mrSges(6,3) * t42;
t25 = mrSges(6,1) * t90 - mrSges(6,3) * t41;
t24 = t192 * t76 + t195 * t31;
t23 = -t192 * t31 + t195 * t76;
t18 = -pkin(4) * t146 - t20;
t17 = pkin(10) * t146 + t19;
t11 = -pkin(4) * t136 - t13;
t6 = t41 * Ifges(6,4) + t42 * Ifges(6,2) + t90 * Ifges(6,6);
t4 = -qJD(5) * t16 - t17 * t192 + t195 * t30;
t3 = qJD(5) * t15 + t17 * t195 + t192 * t30;
t21 = [m(5) * (t12 * t286 + t13 * t45 + t19 * t32 + t20 * t31) + t286 * t64 + t104 * t314 + (Ifges(6,5) * t51 + Ifges(6,6) * t50 + Ifges(6,3) * t100) * t302 + (Ifges(6,1) * t51 + Ifges(6,4) * t50 + Ifges(6,5) * t100) * t305 + (Ifges(6,4) * t51 + Ifges(6,2) * t50 + Ifges(6,6) * t100) * t307 + (Ifges(6,5) * t104 + Ifges(6,6) * t103 + Ifges(6,3) * t125) * t309 + t126 * t310 + (Ifges(6,4) * t104 + Ifges(6,2) * t103 + Ifges(6,6) * t125) * t311 + (Ifges(6,1) * t104 + Ifges(6,4) * t103 + Ifges(6,5) * t125) * t312 + (-m(4) * t95 + m(5) * t81 - mrSges(4,1) * t169 + mrSges(5,1) * t125 + mrSges(5,2) * t126 + mrSges(4,3) * t151) * t66 + m(6) * (t1 * t16 + t11 * t33 + t15 * t2 + t18 * t28 + t3 * t9 + t4 * t8) + (m(3) * ((t189 * t253 + (qJ(2) * t264 - t184) * t186) * qJD(1) + t216) + 0.2e1 * t259) * t251 + m(4) * (t65 * t96 + t72 * t84 + (qJD(1) * t122 + t117) * t235) + (-m(4) * t83 + m(5) * t69 - t265) * t73 + t111 * t235 + (Ifges(5,4) * t126 - Ifges(5,2) * t125 - Ifges(5,6) * t202) * t327 + (Ifges(5,1) * t126 - Ifges(5,4) * t125 - Ifges(5,5) * t202) * t328 + (-Ifges(4,4) * t151 - Ifges(4,6) * t169 / 0.2e1 + t126 * t329 + Ifges(5,6) * t323 - t96 * mrSges(4,3) - (Ifges(4,2) + Ifges(5,3) / 0.2e1) * t202) * t136 + (-t95 * mrSges(4,3) + Ifges(4,1) * t151 + Ifges(4,4) * t202 + Ifges(4,5) * t322) * t135 + t65 * (-mrSges(4,2) * t169 + mrSges(4,3) * t202) + t12 * (mrSges(5,2) * t202 - mrSges(5,3) * t125) + t13 * (-mrSges(5,1) * t202 - mrSges(5,3) * t126) + (-mrSges(4,1) * t202 + mrSges(4,2) * t151) * t215 - t202 * t244 / 0.2e1 + t163 * (Ifges(4,5) * t145 - Ifges(4,6) * t146) / 0.2e1 + t147 * (Ifges(4,1) * t145 - Ifges(4,4) * t146) / 0.2e1 + t145 * t99 / 0.2e1 + t142 * (Ifges(5,5) * t101 - Ifges(5,6) * t100 + Ifges(5,3) * t146) / 0.2e1 + t121 * (Ifges(5,1) * t101 - Ifges(5,4) * t100 + Ifges(5,5) * t146) / 0.2e1 + t120 * (Ifges(5,4) * t101 - Ifges(5,2) * t100 + Ifges(5,6) * t146) / 0.2e1 + t32 * (-mrSges(5,2) * t146 - mrSges(5,3) * t100) + t31 * (mrSges(5,1) * t146 - mrSges(5,3) * t101) + t146 * t59 / 0.2e1 + t144 * (Ifges(4,4) * t145 - Ifges(4,2) * t146) / 0.2e1 + t117 * (mrSges(4,1) * t146 + mrSges(4,2) * t145) - t146 * t98 / 0.2e1 + t1 * (-mrSges(6,2) * t125 + mrSges(6,3) * t103) + t2 * (mrSges(6,1) * t125 - mrSges(6,3) * t104) + t125 * t5 / 0.2e1 + t122 * t107 + t72 * t123 + t9 * (-mrSges(6,2) * t100 + mrSges(6,3) * t50) + t8 * (mrSges(6,1) * t100 - mrSges(6,3) * t51) + t100 * t35 / 0.2e1 - t100 * t60 / 0.2e1 + t69 * (mrSges(5,1) * t100 + mrSges(5,2) * t101) + t101 * t61 / 0.2e1 + t103 * t6 / 0.2e1 + t11 * (-mrSges(6,1) * t103 + mrSges(6,2) * t104) + t19 * t93 + t20 * t94 + t81 * t52 + t45 * t63 + t3 * t57 + t4 * t58 + t18 * t53 + t50 * t36 / 0.2e1 + t28 * (-mrSges(6,1) * t50 + mrSges(6,2) * t51) + t51 * t37 / 0.2e1 + t254 * t322 + t43 * t323 - 0.2e1 * t170 * t239 - t84 * t146 * mrSges(4,3) - t145 * t294 + t15 * t25 + t16 * t26 + t33 * t14; t190 * t107 - t161 * t123 + t157 * t25 - t213 * t26 + t173 * t64 + t318 * t93 + t320 * t58 + t319 * t57 + t288 * t172 + t265 * t160 - m(4) * (-t160 * t83 + t161 * t84) + t317 * t287 + (-m(3) * t216 + t186 * t170 - t259) * t252 + (-t1 * t213 + t11 * t172 + t157 * t2 + t28 * t317 + t319 * t9 + t320 * t8) * m(6) + (t12 * t173 - t13 * t172 - t160 * t69 - t317 * t31 + t318 * t32) * m(5) + (-t111 * t240 - t197 * t52 + (-t135 * t197 - t136 * t194) * mrSges(4,3) + (t123 * t197 - t194 * t265) * qJD(3) + m(5) * (t250 * t69 - t266) + (-t117 * t240 + t190 * t234 + t194 * t65 + t249 * t84 - t250 * t83 - t266) * m(4)) * t187; (Ifges(6,5) * t109 + Ifges(6,6) * t108) * t303 + (Ifges(6,1) * t109 + Ifges(6,4) * t108) * t306 + (Ifges(6,4) * t109 + Ifges(6,2) * t108) * t308 + t108 * t313 + t325 * t58 + t326 * t57 + (t1 * t165 + t164 * t2 - t28 * t47 + t325 * t8 + t326 * t9) * m(6) + (t144 * t222 + t316) * mrSges(5,3) + (-t268 / 0.2e1 - t141 / 0.2e1 - t117 * mrSges(4,2) - t99 / 0.2e1 + t294 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t147) * t144 + ((t278 / 0.2e1 + t118 / 0.2e1 + (m(6) * t28 + t287) * pkin(9) - t324) * qJD(4) + t64 * pkin(9) - t66 * mrSges(5,1) + t43 / 0.2e1 - t5 / 0.2e1 - t40 / 0.2e1 - t39 / 0.2e1 + t274 / 0.2e1 + t292 / 0.2e1 + (-Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t90 + t233 + (-t278 / 0.2e1 + t304 + t209) * t144) * t196 + (t66 * mrSges(5,2) + t11 * t229 + t224 * t309 + t226 * t311 + t228 * t312 + t7 * t300 + t6 * t301 + 0.2e1 * t310 + (-t200 - t295) * qJD(4) + (-t1 * t192 - t195 * t2) * mrSges(6,3) + (t223 * t303 + t225 * t308 + t227 * t306 + t28 * t230 + t195 * t313 + t37 * t301 + (t192 * t8 - t195 * t9) * mrSges(6,3)) * qJD(5) + t200 * t144 + (m(6) * t11 - qJD(4) * t93 + t288) * pkin(9)) * t193 + (-pkin(3) * t66 - t31 * t54 - t32 * t55 - t69 * t84 + (-qJD(4) * t222 + t316) * pkin(9)) * m(5) + (-t9 * t108 + t8 * t109) * mrSges(6,3) + t164 * t25 + t165 * t26 - t83 * t123 - t28 * (-mrSges(6,1) * t108 + mrSges(6,2) * t109) - t109 * t37 / 0.2e1 - t55 * t93 - t54 * t94 - t65 * mrSges(4,2) - t66 * mrSges(4,1) - pkin(3) * t52 - t47 * t53 + (t265 + t270) * t84 + (-t117 * mrSges(4,1) - t271 / 0.2e1 - t277 / 0.2e1 - t279 / 0.2e1 + t267 / 0.2e1 - t59 / 0.2e1 + t98 / 0.2e1 - t31 * mrSges(5,1) + t32 * mrSges(5,2) + t269 / 0.2e1) * t147 + t254; t244 + t315 * mrSges(6,3) - t287 * t32 + t6 * t300 + t192 * t314 - t31 * t93 - t24 * t57 - t23 * t58 + t223 * t309 + t225 * t311 + t227 * t312 - t11 * t230 + (t304 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t121 + t324) * t120 + t330 * qJD(5) + (t201 + t295) * t121 - t12 * mrSges(5,2) + t13 * mrSges(5,1) - pkin(4) * t14 + (-pkin(4) * t11 - t23 * t8 - t24 * t9 - t28 * t32) * m(6) + (m(6) * t315 - t192 * t25 + t195 * t26 + (-m(6) * t231 - t192 * t57 - t195 * t58) * qJD(5)) * pkin(10); -t28 * (mrSges(6,1) * t92 + mrSges(6,2) * t91) + (Ifges(6,1) * t91 - t299) * t306 + t36 * t305 + (Ifges(6,5) * t91 - Ifges(6,6) * t92) * t303 - t8 * t57 + t9 * t58 + (t8 * t91 + t9 * t92) * mrSges(6,3) - t233 + t5 + (-Ifges(6,2) * t92 + t37 + t88) * t308;];
tauc = t21(:);
