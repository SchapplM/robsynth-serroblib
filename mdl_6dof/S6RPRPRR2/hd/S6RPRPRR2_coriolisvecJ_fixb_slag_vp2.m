% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRR2
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

function tauc = S6RPRPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:03:18
% EndTime: 2018-11-23 16:03:26
% DurationCPUTime: 7.74s
% Computational Cost: add. (8931->549), mult. (21511->758), div. (0->0), fcn. (15107->10), ass. (0->267)
t189 = sin(qJ(5));
t258 = sin(pkin(11));
t229 = t258 * pkin(3);
t180 = t229 + pkin(8);
t285 = pkin(9) + t180;
t225 = qJD(5) * t285;
t190 = sin(qJ(3));
t193 = cos(qJ(3));
t220 = qJD(1) * t258;
t259 = cos(pkin(11));
t221 = qJD(1) * t259;
t156 = -t190 * t220 + t193 * t221;
t252 = t156 * t189;
t157 = -t190 * t221 - t193 * t220;
t247 = qJD(1) * t190;
t237 = pkin(3) * t247;
t107 = -pkin(4) * t157 - pkin(8) * t156 + t237;
t192 = cos(qJ(5));
t181 = sin(pkin(10)) * pkin(1) + pkin(7);
t172 = t181 * qJD(1);
t217 = qJ(4) * qJD(1) + t172;
t241 = t190 * qJD(2);
t135 = t193 * t217 + t241;
t126 = t258 * t135;
t185 = t193 * qJD(2);
t134 = -t190 * t217 + t185;
t80 = t134 * t259 - t126;
t49 = t189 * t107 + t192 * t80;
t343 = -pkin(9) * t252 + t189 * t225 + t49;
t251 = t156 * t192;
t48 = t192 * t107 - t189 * t80;
t342 = pkin(5) * t157 + pkin(9) * t251 - t192 * t225 - t48;
t243 = qJD(5) * t189;
t341 = t243 - t252;
t293 = -t156 / 0.2e1;
t340 = Ifges(5,2) * t293 - Ifges(5,6) * qJD(3) / 0.2e1;
t191 = cos(qJ(6));
t188 = sin(qJ(6));
t132 = qJD(3) * t192 + t157 * t189;
t129 = qJD(3) * pkin(3) + t134;
t224 = t259 * t135;
t73 = t258 * t129 + t224;
t70 = qJD(3) * pkin(8) + t73;
t231 = -cos(pkin(10)) * pkin(1) - pkin(2);
t171 = -pkin(3) * t193 + t231;
t248 = qJD(1) * t171;
t155 = qJD(4) + t248;
t88 = -t156 * pkin(4) + t157 * pkin(8) + t155;
t45 = t189 * t88 + t192 * t70;
t33 = pkin(9) * t132 + t45;
t264 = t188 * t33;
t152 = qJD(5) - t156;
t133 = qJD(3) * t189 - t157 * t192;
t44 = -t189 * t70 + t192 * t88;
t32 = -pkin(9) * t133 + t44;
t29 = pkin(5) * t152 + t32;
t10 = t191 * t29 - t264;
t262 = t191 * t33;
t11 = t188 * t29 + t262;
t167 = t190 * t259 + t193 * t258;
t158 = t167 * qJD(3);
t146 = qJD(1) * t158;
t143 = Ifges(7,3) * t146;
t219 = t191 * t132 - t133 * t188;
t78 = t132 * t188 + t133 * t191;
t289 = Ifges(7,4) * t78;
t145 = qJD(6) + t152;
t297 = -t145 / 0.2e1;
t306 = -t78 / 0.2e1;
t72 = t129 * t259 - t126;
t69 = -qJD(3) * pkin(4) - t72;
t53 = -t132 * pkin(5) + t69;
t339 = t143 + (Ifges(7,5) * t219 - Ifges(7,6) * t78) * t297 + (t10 * t219 + t11 * t78) * mrSges(7,3) - t53 * (mrSges(7,1) * t78 + mrSges(7,2) * t219) + (Ifges(7,1) * t219 - t289) * t306;
t242 = qJD(5) * t192;
t245 = qJD(3) * t190;
t140 = qJD(3) * t185 - t172 * t245;
t244 = qJD(4) * t193;
t121 = (-qJ(4) * t245 + t244) * qJD(1) + t140;
t240 = t190 * qJD(4);
t194 = -qJD(1) * t240 - t135 * qJD(3);
t62 = t121 * t259 + t194 * t258;
t166 = t190 * t258 - t193 * t259;
t159 = t166 * qJD(3);
t147 = qJD(1) * t159;
t239 = qJD(1) * qJD(3);
t227 = t190 * t239;
t216 = pkin(3) * t227;
t95 = pkin(4) * t146 + pkin(8) * t147 + t216;
t17 = t189 * t95 + t192 * t62 + t88 * t242 - t243 * t70;
t94 = -qJD(5) * t133 + t147 * t189;
t12 = pkin(9) * t94 + t17;
t18 = -qJD(5) * t45 - t189 * t62 + t192 * t95;
t93 = qJD(5) * t132 - t147 * t192;
t9 = pkin(5) * t146 - pkin(9) * t93 + t18;
t2 = qJD(6) * t10 + t12 * t191 + t188 * t9;
t3 = -qJD(6) * t11 - t12 * t188 + t191 * t9;
t30 = t219 * qJD(6) + t188 * t94 + t191 * t93;
t31 = -qJD(6) * t78 - t188 * t93 + t191 * t94;
t338 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t256 = Ifges(5,5) * qJD(3);
t337 = t155 * mrSges(5,2) + t256 / 0.2e1;
t74 = Ifges(7,4) * t219;
t336 = -Ifges(7,2) * t78 + t74;
t314 = t30 / 0.2e1;
t313 = t31 / 0.2e1;
t295 = t146 / 0.2e1;
t164 = t285 * t189;
t165 = t285 * t192;
t120 = -t164 * t188 + t165 * t191;
t333 = -qJD(6) * t120 + t188 * t343 + t342 * t191;
t119 = -t164 * t191 - t165 * t188;
t332 = qJD(6) * t119 + t342 * t188 - t191 * t343;
t79 = t134 * t258 + t224;
t326 = pkin(5) * t341 - t79;
t270 = t157 * mrSges(5,3);
t260 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t132 + mrSges(6,2) * t133 - t270;
t203 = t188 * t189 - t191 * t192;
t114 = t203 * t167;
t249 = qJ(4) + t181;
t162 = t249 * t190;
t163 = t249 * t193;
t117 = -t162 * t258 + t163 * t259;
t109 = t192 * t117;
t110 = t166 * pkin(4) - t167 * pkin(8) + t171;
t55 = t189 * t110 + t109;
t103 = t203 * t156;
t321 = qJD(5) + qJD(6);
t123 = t321 * t203;
t325 = t103 - t123;
t169 = t188 * t192 + t189 * t191;
t102 = t169 * t156;
t124 = t321 * t169;
t324 = t102 - t124;
t67 = mrSges(6,1) * t146 - mrSges(6,3) * t93;
t68 = -mrSges(6,2) * t146 + mrSges(6,3) * t94;
t323 = -t189 * t67 + t192 * t68;
t322 = t17 * t192 - t18 * t189;
t268 = t157 * Ifges(5,4);
t301 = t268 / 0.2e1 + t340;
t319 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(6,6) * t94 + t338;
t214 = mrSges(6,1) * t189 + mrSges(6,2) * t192;
t201 = t69 * t214;
t209 = Ifges(6,5) * t192 - Ifges(6,6) * t189;
t280 = Ifges(6,4) * t192;
t211 = -Ifges(6,2) * t189 + t280;
t281 = Ifges(6,4) * t189;
t213 = Ifges(6,1) * t192 - t281;
t290 = t192 / 0.2e1;
t291 = -t189 / 0.2e1;
t298 = t133 / 0.2e1;
t282 = Ifges(6,4) * t133;
t65 = Ifges(6,2) * t132 + Ifges(6,6) * t152 + t282;
t130 = Ifges(6,4) * t132;
t66 = t133 * Ifges(6,1) + t152 * Ifges(6,5) + t130;
t318 = t66 * t290 + t65 * t291 + t152 * t209 / 0.2e1 + t213 * t298 + t132 * t211 / 0.2e1 + t201;
t317 = t155 * mrSges(5,1) + t44 * mrSges(6,1) + t10 * mrSges(7,1) - t45 * mrSges(6,2) - t11 * mrSges(7,2) + t301;
t316 = Ifges(7,4) * t314 + Ifges(7,2) * t313 + Ifges(7,6) * t295;
t315 = Ifges(7,1) * t314 + Ifges(7,4) * t313 + Ifges(7,5) * t295;
t36 = Ifges(7,2) * t219 + Ifges(7,6) * t145 + t289;
t312 = -t36 / 0.2e1;
t311 = t36 / 0.2e1;
t37 = Ifges(7,1) * t78 + Ifges(7,5) * t145 + t74;
t310 = -t37 / 0.2e1;
t309 = t37 / 0.2e1;
t308 = -t219 / 0.2e1;
t307 = t219 / 0.2e1;
t305 = t78 / 0.2e1;
t304 = t93 / 0.2e1;
t303 = t94 / 0.2e1;
t300 = -t132 / 0.2e1;
t299 = -t133 / 0.2e1;
t296 = t145 / 0.2e1;
t294 = -t152 / 0.2e1;
t292 = t157 / 0.2e1;
t288 = pkin(9) * t192;
t287 = t219 * Ifges(7,6);
t286 = t78 * Ifges(7,5);
t284 = mrSges(5,3) * t147;
t283 = Ifges(4,4) * t190;
t151 = Ifges(5,4) * t156;
t116 = t259 * t162 + t163 * t258;
t61 = t121 * t258 - t259 * t194;
t278 = t116 * t61;
t277 = t132 * Ifges(6,6);
t276 = t133 * Ifges(6,5);
t275 = t145 * Ifges(7,3);
t274 = t146 * mrSges(5,3);
t273 = t152 * Ifges(6,3);
t272 = t156 * mrSges(5,3);
t269 = t157 * Ifges(5,1);
t267 = t166 * t61;
t257 = Ifges(4,5) * qJD(3);
t255 = Ifges(4,6) * qJD(3);
t250 = t167 * t189;
t246 = qJD(1) * t193;
t46 = -mrSges(7,1) * t219 + mrSges(7,2) * t78;
t238 = t46 + t260;
t235 = mrSges(4,3) * t247;
t234 = mrSges(4,3) * t246;
t184 = Ifges(4,4) * t246;
t230 = t259 * pkin(3);
t228 = m(4) * t181 + mrSges(4,3);
t223 = t146 * mrSges(5,1) - t147 * mrSges(5,2);
t108 = pkin(3) * t245 + pkin(4) * t158 + pkin(8) * t159;
t218 = qJD(3) * t249;
t136 = -t190 * t218 + t244;
t137 = -t193 * t218 - t240;
t87 = t136 * t259 + t137 * t258;
t222 = t192 * t108 - t189 * t87;
t54 = t192 * t110 - t117 * t189;
t182 = -t230 - pkin(4);
t174 = t231 * qJD(1);
t215 = mrSges(6,1) * t192 - mrSges(6,2) * t189;
t212 = Ifges(6,1) * t189 + t280;
t210 = Ifges(6,2) * t192 + t281;
t208 = Ifges(6,5) * t189 + Ifges(6,6) * t192;
t207 = -t17 * t189 - t18 * t192;
t47 = pkin(5) * t166 - t167 * t288 + t54;
t51 = -pkin(9) * t250 + t55;
t21 = -t188 * t51 + t191 * t47;
t22 = t188 * t47 + t191 * t51;
t206 = -t189 * t45 - t192 * t44;
t205 = t189 * t44 - t192 * t45;
t96 = -mrSges(6,2) * t152 + mrSges(6,3) * t132;
t97 = mrSges(6,1) * t152 - mrSges(6,3) * t133;
t204 = -t189 * t97 + t192 * t96;
t86 = t136 * t258 - t259 * t137;
t149 = t172 * t193 + t241;
t138 = -qJD(3) * mrSges(5,2) + t272;
t202 = -t138 - t204;
t200 = -t159 * t189 + t167 * t242;
t25 = t189 * t108 + t110 * t242 - t117 * t243 + t192 * t87;
t195 = qJD(5) * t206 + t322;
t175 = -qJD(3) * mrSges(4,2) + t234;
t173 = qJD(3) * mrSges(4,1) - t235;
t170 = -t192 * pkin(5) + t182;
t161 = Ifges(4,1) * t247 + t184 + t257;
t160 = t255 + (t193 * Ifges(4,2) + t283) * qJD(1);
t148 = -t172 * t190 + t185;
t144 = Ifges(6,3) * t146;
t141 = t149 * qJD(3);
t118 = -mrSges(5,1) * t156 - mrSges(5,2) * t157;
t113 = t169 * t167;
t112 = t151 + t256 - t269;
t84 = pkin(5) * t250 + t116;
t64 = t273 + t276 + t277;
t58 = mrSges(7,1) * t145 - mrSges(7,3) * t78;
t57 = -mrSges(7,2) * t145 + mrSges(7,3) * t219;
t52 = pkin(5) * t200 + t86;
t50 = -mrSges(6,1) * t94 + mrSges(6,2) * t93;
t43 = t93 * Ifges(6,1) + t94 * Ifges(6,4) + t146 * Ifges(6,5);
t42 = t93 * Ifges(6,4) + t94 * Ifges(6,2) + t146 * Ifges(6,6);
t40 = t321 * t114 + t169 * t159;
t39 = -t124 * t167 + t159 * t203;
t38 = -t94 * pkin(5) + t61;
t35 = t275 + t286 + t287;
t26 = -t55 * qJD(5) + t222;
t24 = -mrSges(7,2) * t146 + mrSges(7,3) * t31;
t23 = mrSges(7,1) * t146 - mrSges(7,3) * t30;
t20 = -pkin(9) * t200 + t25;
t19 = t159 * t288 + pkin(5) * t158 + (-t109 + (pkin(9) * t167 - t110) * t189) * qJD(5) + t222;
t14 = t191 * t32 - t264;
t13 = -t188 * t32 - t262;
t8 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t5 = -qJD(6) * t22 - t188 * t20 + t19 * t191;
t4 = qJD(6) * t21 + t188 * t19 + t191 * t20;
t1 = [(t277 / 0.2e1 + t276 / 0.2e1 + t273 / 0.2e1 + t35 / 0.2e1 + t64 / 0.2e1 + t287 / 0.2e1 + t286 / 0.2e1 + t275 / 0.2e1 + t317 + t301) * t158 + m(5) * (t117 * t62 - t72 * t86 + t73 * t87 + t278) + m(6) * (t17 * t55 + t18 * t54 + t25 * t45 + t26 * t44 + t69 * t86 + t278) + (t228 * t140 + (t161 / 0.2e1 - t181 * t173 + 0.3e1 / 0.2e1 * t184 + t257 / 0.2e1 - t228 * t148 + 0.2e1 * t174 * mrSges(4,2)) * qJD(3)) * t193 + t260 * t86 + t171 * t223 + (t228 * t141 + (-t160 / 0.2e1 - t181 * t175 + t174 * mrSges(4,1) - t255 / 0.2e1 - t228 * t149 + (t231 * mrSges(4,1) - 0.3e1 / 0.2e1 * t283 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t193) * qJD(1) + (qJD(1) * (mrSges(5,1) * t166 + mrSges(5,2) * t167) + m(5) * (t155 + t248) + t118) * pkin(3)) * qJD(3)) * t190 - (t112 / 0.2e1 + t151 / 0.2e1 - t269 / 0.2e1 + t206 * mrSges(6,3) + t318 + t337) * t159 + m(7) * (t10 * t5 + t11 * t4 + t2 * t22 + t21 * t3 + t38 * t84 + t52 * t53) + (-t10 * t39 + t11 * t40 - t113 * t2 + t114 * t3) * mrSges(7,3) + (-Ifges(7,5) * t114 - Ifges(7,6) * t113) * t295 + (-Ifges(7,4) * t114 - Ifges(7,2) * t113) * t313 + (-Ifges(7,1) * t114 - Ifges(7,4) * t113) * t314 + t38 * (mrSges(7,1) * t113 - mrSges(7,2) * t114) + t21 * t23 + t22 * t24 + t52 * t46 + t53 * (-mrSges(7,1) * t40 + mrSges(7,2) * t39) + (Ifges(7,5) * t39 + Ifges(7,6) * t40) * t296 + (-t62 * mrSges(5,3) + t143 / 0.2e1 + t144 / 0.2e1 + Ifges(5,4) * t147 + (Ifges(7,3) / 0.2e1 + Ifges(5,2) + Ifges(6,3) / 0.2e1) * t146 + t319) * t166 + (-t116 * t147 - t117 * t146 - t158 * t73 + t159 * t72) * mrSges(5,3) + (t42 * t291 + t43 * t290 - Ifges(5,1) * t147 - Ifges(5,4) * t146 + t213 * t304 + t211 * t303 + t209 * t295 + (mrSges(5,3) + t214) * t61 + t207 * mrSges(6,3) + (t66 * t291 - t192 * t65 / 0.2e1 + t69 * t215 + t210 * t300 + t212 * t299 + t208 * t294 + t205 * mrSges(6,3)) * qJD(5)) * t167 + (Ifges(7,1) * t39 + Ifges(7,4) * t40) * t305 + (Ifges(7,4) * t39 + Ifges(7,2) * t40) * t307 + t39 * t309 + t40 * t311 - t114 * t315 - t113 * t316 + t4 * t57 + t5 * t58 + t54 * t67 + t55 * t68 + t84 * t8 + t25 * t96 + t26 * t97 + t116 * t50 + t87 * t138; -t113 * t23 - t114 * t24 + t39 * t57 + t40 * t58 + (t50 + t8 - t284) * t166 + t202 * t159 + t238 * t158 + (-t190 * t173 + t193 * t175 + (-t190 ^ 2 - t193 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t274 + (-t189 * t96 - t192 * t97) * qJD(5) + t323) * t167 + m(7) * (t10 * t40 + t11 * t39 - t113 * t3 - t114 * t2 + t158 * t53 + t166 * t38) + m(5) * (-t158 * t72 - t159 * t73 + t167 * t62 + t267) + m(6) * (t158 * t69 + t159 * t205 + t167 * t195 + t267) + m(4) * (t140 * t190 - t141 * t193 + (-t148 * t190 + t149 * t193) * qJD(3)); (-m(6) * t69 - t260) * t79 + t326 * t46 + (-t324 * mrSges(7,1) + t325 * mrSges(7,2)) * t53 + (-t174 * (mrSges(4,1) * t190 + mrSges(4,2) * t193) - (Ifges(4,1) * t193 - t283) * t247 / 0.2e1) * qJD(1) + (m(6) * t195 - t242 * t97 - t243 * t96 + t323) * t180 + (-t155 * t237 + t72 * t79 - t73 * t80 + (t258 * t62 - t259 * t61) * pkin(3)) * m(5) + t318 * qJD(5) + t230 * t284 - t73 * t270 - t229 * t274 - t66 * t251 / 0.2e1 + t65 * t252 / 0.2e1 + t160 * t247 / 0.2e1 - t118 * t237 - Ifges(4,6) * t227 / 0.2e1 + t239 * Ifges(4,5) * t193 / 0.2e1 + (-t341 * t45 + (-t242 + t251) * t44 + t322) * mrSges(6,3) + (-Ifges(6,5) * t299 - Ifges(7,5) * t306 - Ifges(6,6) * t300 - Ifges(7,6) * t308 - Ifges(6,3) * t294 - Ifges(7,3) * t297 + t317 + t340) * t157 - (-Ifges(4,2) * t247 + t161 + t184) * t246 / 0.2e1 + (t235 + t173) * t149 + (t234 - t175) * t148 + (t268 + t64 + t35) * t292 + (t151 + t112) * t293 + (m(6) * t182 - mrSges(5,1) - t215) * t61 - m(6) * (t44 * t48 + t45 * t49) + t72 * t272 + (Ifges(5,1) * t292 + t209 * t294 + t211 * t300 + t213 * t299 - t201 - t337) * t156 + (-Ifges(7,1) * t103 - Ifges(7,4) * t102) * t306 + (-Ifges(7,4) * t103 - Ifges(7,2) * t102) * t308 + (-Ifges(7,5) * t103 - Ifges(7,6) * t102) * t297 + t332 * t57 + t333 * t58 + (t333 * t10 + t332 * t11 + t119 * t3 + t120 * t2 + t170 * t38 + t326 * t53) * m(7) + t42 * t290 + (-t325 * t10 + t324 * t11 - t169 * t3 - t2 * t203) * mrSges(7,3) + (Ifges(7,5) * t169 - Ifges(7,6) * t203 + t208) * t295 + (Ifges(7,4) * t169 - Ifges(7,2) * t203) * t313 + (Ifges(7,1) * t169 - Ifges(7,4) * t203) * t314 + t38 * (mrSges(7,1) * t203 + mrSges(7,2) * t169) - t203 * t316 + (-Ifges(7,5) * t123 - Ifges(7,6) * t124) * t296 + (-Ifges(7,1) * t123 - Ifges(7,4) * t124) * t305 + (-Ifges(7,4) * t123 - Ifges(7,2) * t124) * t307 + t210 * t303 + t212 * t304 - t123 * t309 - t103 * t310 - t124 * t311 - t102 * t312 + t169 * t315 - t62 * mrSges(5,2) - t49 * t96 - t48 * t97 + t119 * t23 + t120 * t24 - t80 * t138 - t140 * mrSges(4,2) - t141 * mrSges(4,1) - Ifges(5,6) * t146 - Ifges(5,5) * t147 + t170 * t8 + t182 * t50 + t189 * t43 / 0.2e1; -t203 * t23 + t169 * t24 + t189 * t68 + t192 * t67 + t324 * t58 + t325 * t57 + t204 * qJD(5) + t238 * t157 + t202 * t156 + t223 + (t324 * t10 + t325 * t11 + t157 * t53 + t169 * t2 - t203 * t3) * m(7) + (-t152 * t205 + t157 * t69 - t207) * m(6) + (-t156 * t73 - t157 * t72 + t216) * m(5); t219 * t310 + t319 + t144 - m(7) * (t10 * t13 + t11 * t14) + (t132 * t44 + t133 * t45) * mrSges(6,3) + (-Ifges(6,2) * t133 + t130 + t66) * t300 + (Ifges(6,5) * t132 - Ifges(6,6) * t133) * t294 + t65 * t298 + (Ifges(6,1) * t132 - t282) * t299 - t78 * t312 + t336 * t308 - t14 * t57 - t13 * t58 - t44 * t96 + t45 * t97 - t69 * (mrSges(6,1) * t133 + mrSges(6,2) * t132) + (-t133 * t46 + t188 * t24 + t191 * t23 + (-t188 * t58 + t191 * t57) * qJD(6) + (-t133 * t53 + t188 * t2 + t191 * t3 + (-t10 * t188 + t11 * t191) * qJD(6)) * m(7)) * pkin(5) + t339; t36 * t305 - t10 * t57 + t11 * t58 + (t336 + t37) * t308 + t338 + t339;];
tauc  = t1(:);
