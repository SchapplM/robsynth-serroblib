% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:30
% EndTime: 2019-03-09 02:25:46
% DurationCPUTime: 8.93s
% Computational Cost: add. (6236->526), mult. (12637->759), div. (0->0), fcn. (7487->8), ass. (0->242)
t173 = cos(qJ(4));
t164 = t173 * qJD(1);
t169 = sin(qJ(5));
t293 = -pkin(9) - pkin(8);
t219 = qJD(5) * t293;
t170 = sin(qJ(4));
t199 = -pkin(4) * t170 + pkin(8) * t173;
t141 = t199 * qJD(1);
t172 = cos(qJ(5));
t174 = -pkin(1) - pkin(2);
t155 = qJD(1) * t174 + qJD(2);
t165 = sin(pkin(10));
t166 = cos(pkin(10));
t243 = qJD(1) * t166;
t120 = qJ(2) * t243 + t165 * t155;
t109 = -qJD(1) * pkin(7) + t120;
t86 = qJD(3) * t173 - t170 * t109;
t55 = t169 * t141 + t172 * t86;
t327 = -t55 + (-pkin(9) * t164 + t219) * t169;
t245 = t172 * t173;
t185 = -pkin(5) * t170 + pkin(9) * t245;
t54 = t172 * t141 - t169 * t86;
t326 = -qJD(1) * t185 + t172 * t219 - t54;
t248 = t169 * t173;
t126 = -t165 * t248 - t166 * t172;
t236 = qJD(4) * t172;
t212 = t170 * t236;
t319 = qJD(5) * t126 - t165 * t212 - (t165 * t169 + t166 * t245) * qJD(1);
t127 = t165 * t245 - t166 * t169;
t231 = t169 * qJD(4);
t214 = t170 * t231;
t318 = -qJD(5) * t127 + t165 * t214 - (t165 * t172 - t166 * t248) * qJD(1);
t235 = qJD(4) * t173;
t211 = t172 * t235;
t233 = qJD(5) * t170;
t178 = -t169 * t233 + t211;
t228 = qJD(4) * qJD(5);
t101 = -qJD(1) * t178 + t172 * t228;
t241 = qJD(2) * t165;
t124 = qJD(4) * t199 + t241;
t110 = t124 * qJD(1);
t239 = qJD(3) * t170;
t87 = t109 * t173 + t239;
t80 = qJD(4) * pkin(8) + t87;
t161 = t165 * qJ(2);
t119 = -qJD(1) * t161 + t155 * t166;
t108 = qJD(1) * pkin(3) - t119;
t200 = pkin(4) * t173 + pkin(8) * t170;
t81 = qJD(1) * t200 + t108;
t40 = t169 * t81 + t172 * t80;
t230 = qJD(1) * qJD(2);
t209 = t166 * t230;
t64 = qJD(4) * t86 + t173 * t209;
t16 = -qJD(5) * t40 + t172 * t110 - t169 * t64;
t229 = qJD(1) * qJD(4);
t208 = t170 * t229;
t13 = -pkin(5) * t208 - pkin(9) * t101 + t16;
t213 = t173 * t231;
t232 = qJD(5) * t172;
t177 = t170 * t232 + t213;
t102 = qJD(1) * t177 - t169 * t228;
t234 = qJD(5) * t169;
t15 = t169 * t110 + t172 * t64 + t81 * t232 - t234 * t80;
t14 = pkin(9) * t102 + t15;
t168 = sin(qJ(6));
t171 = cos(qJ(6));
t242 = qJD(1) * t170;
t134 = t169 * t242 + t236;
t35 = pkin(9) * t134 + t40;
t261 = t168 * t35;
t156 = t164 + qJD(5);
t135 = t172 * t242 - t231;
t39 = -t169 * t80 + t172 * t81;
t34 = pkin(9) * t135 + t39;
t29 = pkin(5) * t156 + t34;
t8 = t171 * t29 - t261;
t2 = qJD(6) * t8 + t13 * t168 + t14 * t171;
t203 = t171 * t134 + t135 * t168;
t85 = t134 * t168 - t135 * t171;
t280 = Ifges(7,4) * t85;
t227 = qJD(5) + qJD(6);
t154 = t164 + t227;
t285 = -t154 / 0.2e1;
t295 = -t85 / 0.2e1;
t297 = -t203 / 0.2e1;
t259 = t171 * t35;
t9 = t168 * t29 + t259;
t3 = -qJD(6) * t9 + t13 * t171 - t14 * t168;
t75 = Ifges(7,4) * t203;
t38 = Ifges(7,1) * t85 + Ifges(7,5) * t154 + t75;
t79 = -qJD(4) * pkin(4) - t86;
t56 = -pkin(5) * t134 + t79;
t325 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + (Ifges(7,5) * t203 - Ifges(7,6) * t85) * t285 + (t203 * t8 + t85 * t9) * mrSges(7,3) + (-Ifges(7,2) * t85 + t38 + t75) * t297 - t56 * (mrSges(7,1) * t85 + mrSges(7,2) * t203) + (Ifges(7,1) * t203 - t280) * t295;
t151 = t293 * t169;
t152 = t293 * t172;
t97 = t151 * t168 - t152 * t171;
t324 = -qJD(6) * t97 - t327 * t168 + t326 * t171;
t96 = t151 * t171 + t152 * t168;
t323 = qJD(6) * t96 + t326 * t168 + t327 * t171;
t67 = t126 * t171 - t127 * t168;
t322 = qJD(6) * t67 + t318 * t168 + t319 * t171;
t68 = t126 * t168 + t127 * t171;
t321 = -qJD(6) * t68 - t319 * t168 + t318 * t171;
t224 = mrSges(5,3) * t242;
t320 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t134 - mrSges(6,2) * t135 - t224;
t251 = qJD(4) * t87;
t65 = t170 * t209 + t251;
t260 = t170 * t65;
t188 = t173 * t64 + t260;
t237 = qJD(4) * t170;
t317 = t86 * t235 + t87 * t237 - t188;
t32 = qJD(6) * t203 + t101 * t171 + t102 * t168;
t300 = t32 / 0.2e1;
t33 = -qJD(6) * t85 - t101 * t168 + t102 * t171;
t299 = t33 / 0.2e1;
t37 = Ifges(7,2) * t203 + Ifges(7,6) * t154 + t280;
t315 = t37 / 0.2e1;
t290 = t101 / 0.2e1;
t289 = t102 / 0.2e1;
t314 = -t208 / 0.2e1;
t313 = m(6) * (t235 * t79 + t260);
t311 = Ifges(6,3) + Ifges(7,3);
t306 = -m(3) * qJ(2) - mrSges(3,3);
t292 = m(7) * t56;
t41 = -mrSges(7,1) * t203 + mrSges(7,2) * t85;
t305 = t41 + t292;
t139 = t168 * t172 + t169 * t171;
t122 = t139 * t170;
t202 = t166 * t174 - t161;
t136 = pkin(3) - t202;
t111 = t136 + t200;
t244 = t166 * qJ(2) + t165 * t174;
t137 = -pkin(7) + t244;
t121 = t137 * t245;
t63 = t169 * t111 + t121;
t240 = qJD(2) * t166;
t217 = t170 * t240;
t304 = t137 * t235 + t217;
t191 = t15 * t172 - t16 * t169;
t190 = -t169 * t40 - t172 * t39;
t193 = Ifges(6,5) * t172 - Ifges(6,6) * t169;
t268 = Ifges(6,4) * t172;
t195 = -Ifges(6,2) * t169 + t268;
t269 = Ifges(6,4) * t169;
t197 = Ifges(6,1) * t172 - t269;
t198 = mrSges(6,1) * t169 + mrSges(6,2) * t172;
t133 = Ifges(6,4) * t134;
t72 = -t135 * Ifges(6,1) + t156 * Ifges(6,5) + t133;
t258 = t172 * t72;
t283 = t156 / 0.2e1;
t286 = -t135 / 0.2e1;
t287 = t134 / 0.2e1;
t270 = Ifges(6,4) * t135;
t71 = t134 * Ifges(6,2) + t156 * Ifges(6,6) - t270;
t302 = t193 * t283 + t195 * t287 + t197 * t286 + t79 * t198 + t190 * mrSges(6,3) - t169 * t71 / 0.2e1 + t258 / 0.2e1;
t301 = Ifges(7,4) * t300 + Ifges(7,2) * t299 + Ifges(7,6) * t314;
t298 = Ifges(6,4) * t290 + Ifges(6,2) * t289 + Ifges(6,6) * t314;
t296 = t203 / 0.2e1;
t294 = t85 / 0.2e1;
t284 = t154 / 0.2e1;
t279 = pkin(5) * t169;
t278 = t203 * Ifges(7,6);
t277 = t85 * Ifges(7,5);
t10 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t52 = -mrSges(6,1) * t102 + mrSges(6,2) * t101;
t275 = t10 + t52;
t274 = Ifges(7,5) * t32 + Ifges(7,6) * t33;
t273 = Ifges(6,5) * t101 + Ifges(6,6) * t102;
t272 = Ifges(5,4) * t170;
t271 = Ifges(5,4) * t173;
t267 = t134 * Ifges(6,6);
t266 = t135 * Ifges(6,5);
t264 = t154 * Ifges(7,3);
t263 = t156 * Ifges(6,3);
t112 = t139 * t164;
t89 = t227 * t139;
t256 = t112 + t89;
t246 = t171 * t172;
t186 = t168 * t169 - t246;
t113 = t186 * t164;
t88 = t227 * t186;
t255 = t113 + t88;
t253 = Ifges(5,5) * qJD(4);
t252 = Ifges(5,6) * qJD(4);
t249 = t169 * t170;
t247 = t170 * t172;
t226 = -t41 - t320;
t223 = mrSges(5,3) * t164;
t220 = t168 * t249;
t218 = t166 * t242;
t216 = t173 * t240;
t196 = Ifges(6,1) * t169 + t268;
t194 = Ifges(6,2) * t172 + t269;
t192 = Ifges(6,5) * t169 + Ifges(6,6) * t172;
t104 = t172 * t111;
t51 = pkin(9) * t247 + t104 + (-t137 * t169 + pkin(5)) * t173;
t53 = pkin(9) * t249 + t63;
t21 = -t168 * t53 + t171 * t51;
t22 = t168 * t51 + t171 * t53;
t76 = -mrSges(6,1) * t208 - mrSges(6,3) * t101;
t77 = mrSges(6,2) * t208 + mrSges(6,3) * t102;
t189 = -t169 * t76 + t172 * t77;
t105 = -mrSges(6,2) * t156 + mrSges(6,3) * t134;
t106 = mrSges(6,1) * t156 + mrSges(6,3) * t135;
t187 = -t169 * t105 - t172 * t106;
t183 = -Ifges(7,3) * t208 + t274;
t140 = (mrSges(5,1) * t173 - mrSges(5,2) * t170) * qJD(1);
t182 = (-mrSges(5,1) * t170 - mrSges(5,2) * t173) * qJD(4);
t181 = (-t170 * Ifges(5,1) - t271) * qJD(1);
t180 = (-Ifges(5,2) * t173 - t272) * qJD(1);
t179 = t172 * t124 + t137 * t214 - t169 * t216;
t27 = t169 * t124 + t172 * t216 + t111 * t232 + (-t173 * t234 - t212) * t137;
t175 = qJD(1) ^ 2;
t160 = -pkin(5) * t172 - pkin(4);
t150 = -qJD(4) * mrSges(5,2) - t223;
t132 = qJD(1) * t182;
t129 = t181 + t253;
t128 = t180 + t252;
t123 = t170 * t246 - t220;
t107 = (t137 - t279) * t170;
t70 = t263 - t266 + t267;
t69 = t239 + (-qJD(1) * t279 + t109) * t173;
t66 = -pkin(5) * t177 + t304;
t62 = -t137 * t248 + t104;
t59 = mrSges(7,1) * t154 - mrSges(7,3) * t85;
t58 = -mrSges(7,2) * t154 + mrSges(7,3) * t203;
t47 = t101 * Ifges(6,1) + t102 * Ifges(6,4) - Ifges(6,5) * t208;
t45 = -qJD(6) * t220 + (t227 * t247 + t213) * t171 + t178 * t168;
t44 = t122 * t227 + t168 * t213 - t171 * t211;
t42 = -pkin(5) * t102 + t65;
t36 = t264 + t277 + t278;
t28 = -qJD(5) * t63 + t179;
t26 = mrSges(7,2) * t208 + mrSges(7,3) * t33;
t25 = -mrSges(7,1) * t208 - mrSges(7,3) * t32;
t20 = pkin(9) * t177 + t27;
t19 = t185 * qJD(4) + (-t121 + (-pkin(9) * t170 - t111) * t169) * qJD(5) + t179;
t12 = t171 * t34 - t261;
t11 = -t168 * t34 - t259;
t7 = t32 * Ifges(7,1) + t33 * Ifges(7,4) - Ifges(7,5) * t208;
t5 = -qJD(6) * t22 - t168 * t20 + t171 * t19;
t4 = qJD(6) * t21 + t168 * t19 + t171 * t20;
t1 = [-(t70 + t36 + (-Ifges(7,5) * t123 + Ifges(7,6) * t122 - t170 * t193 + t173 * t311) * qJD(1)) * t237 / 0.2e1 + m(6) * (t15 * t63 + t16 * t62 + t217 * t79 + t27 * t40 + t28 * t39) + (m(5) * ((-t170 * t86 + t173 * t87) * t166 + (qJD(1) * t136 + t108) * t165) + m(4) * (-t119 * t165 + t120 * t166 + (-t165 * t202 + t166 * t244) * qJD(1))) * qJD(2) + (t169 * t72 + t172 * t71) * t233 / 0.2e1 - (t258 + t181 + t129) * t235 / 0.2e1 + (t128 + t180) * t237 / 0.2e1 + (-Ifges(6,3) * t208 + t183 + t273) * t173 / 0.2e1 + (m(5) * ((-t170 * t87 - t173 * t86) * qJD(4) + t188) - t150 * t237 + t170 * t52 + t313) * t137 + 0.2e1 * (mrSges(4,1) * t165 - t306) * t230 + t320 * t304 + t317 * mrSges(5,3) + (t194 * t233 + (-Ifges(6,6) * t170 - t173 * t195) * qJD(4)) * t287 + (Ifges(6,6) * t173 - t170 * t195) * t289 + (Ifges(6,5) * t173 - t170 * t197) * t290 + (Ifges(7,1) * t44 + Ifges(7,4) * t45 - Ifges(7,5) * t237) * t294 + (Ifges(7,4) * t44 + Ifges(7,2) * t45 - Ifges(7,6) * t237) * t296 + t249 * t298 + (-Ifges(7,4) * t123 + Ifges(7,2) * t122 + Ifges(7,6) * t173) * t299 + (-Ifges(7,1) * t123 + Ifges(7,4) * t122 + Ifges(7,5) * t173) * t300 + t122 * t301 + (t192 * t233 + (-Ifges(6,3) * t170 - t173 * t193) * qJD(4)) * t283 + (Ifges(7,5) * t44 + Ifges(7,6) * t45 - Ifges(7,3) * t237) * t284 + (t196 * t233 + (-Ifges(6,5) * t170 - t173 * t197) * qJD(4)) * t286 + t79 * (-mrSges(6,1) * t177 - mrSges(6,2) * t178) + qJD(4) ^ 2 * (-Ifges(5,5) * t173 + Ifges(5,6) * t170) / 0.2e1 + t108 * t182 + m(7) * (t107 * t42 + t2 * t22 + t21 * t3 + t4 * t9 + t5 * t8 + t56 * t66) + t45 * t315 - t173 * (Ifges(5,2) * t170 - t271) * t229 + t150 * t216 + t21 * t25 + t22 * t26 + t44 * t38 / 0.2e1 + t56 * (-mrSges(7,1) * t45 + mrSges(7,2) * t44) + t4 * t58 + t5 * t59 + t66 * t41 + 0.2e1 * mrSges(4,2) * t209 + t62 * t76 + t63 * t77 + t71 * t213 / 0.2e1 + t27 * t105 + t28 * t106 + t107 * t10 - t123 * t7 / 0.2e1 + t42 * (-mrSges(7,1) * t122 - mrSges(7,2) * t123) + t136 * t132 + t8 * (-mrSges(7,1) * t237 - mrSges(7,3) * t44) + t39 * (-mrSges(6,1) * t237 + mrSges(6,3) * t178) + t9 * (mrSges(7,2) * t237 + mrSges(7,3) * t45) + t40 * (mrSges(6,2) * t237 + mrSges(6,3) * t177) + 0.2e1 * t140 * t241 + t2 * (-mrSges(7,2) * t173 + mrSges(7,3) * t122) + t3 * (mrSges(7,1) * t173 + mrSges(7,3) * t123) - t47 * t247 / 0.2e1 + t16 * (mrSges(6,1) * t173 + mrSges(6,3) * t247) + t15 * (-mrSges(6,2) * t173 + mrSges(6,3) * t249) - t198 * t260 - (-Ifges(5,1) * t173 + t272) * t208; t126 * t76 + t127 * t77 + t67 * t25 + t68 * t26 + t321 * t59 + t322 * t58 + t306 * t175 + t318 * t106 + t319 * t105 + (-t175 * mrSges(4,2) - t132 + (-t173 * t150 + t170 * t226) * qJD(1)) * t166 - m(4) * t120 * t243 - m(5) * (t164 * t166 * t87 - t218 * t86) + (t2 * t68 - t218 * t56 + t3 * t67 + t321 * t8 + t322 * t9) * m(7) + (t126 * t16 + t127 * t15 - t218 * t79 + t318 * t39 + t319 * t40) * m(6) + (-t175 * mrSges(4,1) + t275 * t170 + (-t150 * t170 - t173 * t226) * qJD(4) + t313 + m(7) * (t170 * t42 + t235 * t56) + m(5) * (-t209 - t317) + (m(4) * t119 - m(5) * t108 - t140) * qJD(1)) * t165; t123 * t26 - t122 * t25 - t44 * t58 - t45 * t59 + m(7) * (-t122 * t3 + t123 * t2 - t44 * t9 - t45 * t8) + ((t105 * t172 - t106 * t169 + t150 + t223) * qJD(4) + m(6) * (-t231 * t39 + t236 * t40 - t65) - m(7) * t42 + m(5) * (-t65 + t251) - t275) * t173 + (t187 * qJD(5) + m(6) * (-t232 * t39 - t234 * t40 + t191) + m(5) * t64 + t189 + (-m(5) * t86 + m(6) * t79 + t224 - t226 + t292) * qJD(4)) * t170; (-Ifges(7,1) * t88 - Ifges(7,4) * t89) * t294 + (-Ifges(7,4) * t88 - Ifges(7,2) * t89) * t296 + (-Ifges(7,5) * t88 - Ifges(7,6) * t89) * t284 + (-t89 / 0.2e1 - t112 / 0.2e1) * t37 + t42 * (mrSges(7,1) * t186 + mrSges(7,2) * t139) + (-t139 * t3 - t186 * t2 + t255 * t8 - t256 * t9) * mrSges(7,3) + (Ifges(7,4) * t139 - Ifges(7,2) * t186) * t299 + (Ifges(7,1) * t139 - Ifges(7,4) * t186) * t300 - t186 * t301 + ((-t86 * mrSges(5,3) - t253 / 0.2e1 + t108 * mrSges(5,2) + t129 / 0.2e1 - Ifges(5,4) * t164 / 0.2e1 + t302) * t173 + (-t87 * mrSges(5,3) + t8 * mrSges(7,1) - t9 * mrSges(7,2) + t278 / 0.2e1 + t277 / 0.2e1 + t264 / 0.2e1 + t252 / 0.2e1 + t108 * mrSges(5,1) + t39 * mrSges(6,1) - t40 * mrSges(6,2) + t267 / 0.2e1 - t266 / 0.2e1 + t263 / 0.2e1 - t128 / 0.2e1 + t36 / 0.2e1 + t70 / 0.2e1 + (t272 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t173) * qJD(1) - (Ifges(7,5) * t139 - Ifges(7,6) * t186 + t192) * qJD(4) / 0.2e1) * t170) * qJD(1) + (-t88 / 0.2e1 - t113 / 0.2e1) * t38 + (t305 * t279 + (m(6) * t190 + t187) * pkin(8) + t302) * qJD(5) + m(6) * (-pkin(4) * t65 + pkin(8) * t191) + t323 * t58 + (t160 * t42 + t2 * t97 + t3 * t96 + t323 * t9 + t324 * t8 - t56 * t69) * m(7) + t324 * t59 + t194 * t289 + t196 * t290 + (Ifges(7,1) * t113 + Ifges(7,4) * t112) * t295 + (Ifges(7,4) * t113 + Ifges(7,2) * t112) * t297 + t172 * t298 + (Ifges(7,5) * t113 + Ifges(7,6) * t112) * t285 + (-mrSges(6,1) * t172 + mrSges(6,2) * t169 - mrSges(5,1)) * t65 + t189 * pkin(8) + t191 * mrSges(6,3) - t320 * t87 - m(6) * (t39 * t54 + t40 * t55 + t79 * t87) - pkin(4) * t52 - t64 * mrSges(5,2) - t69 * t41 + t96 * t25 + t97 * t26 - t55 * t105 - t54 * t106 + t139 * t7 / 0.2e1 - t86 * t150 + t160 * t10 + t169 * t47 / 0.2e1 + (mrSges(7,1) * t256 - mrSges(7,2) * t255) * t56; -(Ifges(6,2) * t135 + t133 + t72) * t134 / 0.2e1 + t85 * t315 - t311 * t208 + (t168 * t26 + t171 * t25 + m(7) * (t168 * t2 + t171 * t3) + t305 * t135 + (-t168 * t59 + t171 * t58 + m(7) * (-t168 * t8 + t171 * t9)) * qJD(6)) * pkin(5) + t71 * t286 + (t134 * t39 - t135 * t40) * mrSges(6,3) + t273 + t274 - m(7) * (t11 * t8 + t12 * t9) - t15 * mrSges(6,2) + t16 * mrSges(6,1) - t12 * t58 - t11 * t59 - t39 * t105 + t40 * t106 - t79 * (-mrSges(6,1) * t135 + mrSges(6,2) * t134) - t156 * (Ifges(6,5) * t134 + Ifges(6,6) * t135) / 0.2e1 + t135 * (Ifges(6,1) * t134 + t270) / 0.2e1 + t325; t37 * t294 - t8 * t58 + t9 * t59 + t183 + t325;];
tauc  = t1(:);
