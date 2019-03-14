% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:16
% EndTime: 2019-03-09 09:15:31
% DurationCPUTime: 5.68s
% Computational Cost: add. (10105->463), mult. (23670->611), div. (0->0), fcn. (16830->8), ass. (0->242)
t188 = sin(qJ(6));
t191 = cos(qJ(6));
t186 = cos(pkin(10));
t190 = sin(qJ(2));
t247 = qJD(1) * qJD(2);
t240 = t190 * t247;
t157 = t186 * t240;
t185 = sin(pkin(10));
t192 = cos(qJ(2));
t239 = t192 * t247;
t113 = t185 * t239 - t157;
t212 = t185 * t190 + t186 * t192;
t204 = t212 * qJD(2);
t114 = qJD(1) * t204;
t254 = qJD(1) * t192;
t255 = qJD(1) * t190;
t121 = -t185 * t255 - t186 * t254;
t242 = t185 * t254;
t123 = t186 * t255 - t242;
t189 = sin(qJ(5));
t290 = cos(qJ(5));
t241 = qJD(5) * t290;
t250 = qJD(5) * t189;
t199 = t189 * t113 - t290 * t114 - t121 * t241 + t123 * t250;
t207 = t189 * t121 + t290 * t123;
t246 = qJD(2) - qJD(5);
t225 = t191 * t246;
t249 = qJD(6) * t188;
t21 = qJD(6) * t225 + t191 * t199 + t207 * t249;
t277 = t188 * t21;
t300 = -t290 * t121 + t189 * t123;
t309 = qJD(6) + t300;
t322 = t191 * t309;
t55 = -t188 * t246 + t191 * t207;
t330 = t322 * t55 - t277;
t269 = qJD(6) * t55;
t22 = -t188 * t199 + t269;
t248 = qJD(6) * t191;
t53 = t188 * t207 + t225;
t281 = -t188 * t22 - t53 * t248;
t314 = t300 * t53;
t323 = t188 * t309;
t329 = t323 * t55 + t191 * (t21 + t314) - t281;
t227 = t290 * t113 + t189 * t114;
t298 = qJD(5) * t207;
t36 = t227 + t298;
t276 = t188 * t36;
t301 = t207 * t55;
t328 = t309 * t322 + t276 - t301;
t302 = t207 * t53;
t35 = t191 * t36;
t327 = t309 * t323 - t302 - t35;
t288 = pkin(8) * t123;
t172 = pkin(7) * t255;
t143 = qJ(4) * t255 - t172;
t193 = -pkin(2) - pkin(3);
t243 = t193 * qJD(2);
t106 = qJD(3) + t243 - t143;
t173 = pkin(7) * t254;
t145 = -qJ(4) * t254 + t173;
t182 = qJD(2) * qJ(3);
t131 = t145 + t182;
t61 = t186 * t106 - t131 * t185;
t45 = -qJD(2) * pkin(4) - t288 + t61;
t289 = pkin(8) * t121;
t62 = t185 * t106 + t186 * t131;
t47 = t62 + t289;
t25 = -t189 * t47 + t290 * t45;
t23 = t246 * pkin(5) - t25;
t326 = t23 * t309;
t26 = t189 * t45 + t290 * t47;
t137 = t290 * t185 + t189 * t186;
t82 = -t143 * t185 + t145 * t186;
t203 = t82 + t289;
t83 = t186 * t143 + t185 * t145;
t56 = t83 + t288;
t146 = -qJ(3) * t185 + t186 * t193;
t142 = -pkin(4) + t146;
t147 = qJ(3) * t186 + t185 * t193;
t85 = t189 * t142 + t290 * t147;
t279 = t137 * qJD(3) + t85 * qJD(5) - t189 * t56 + t290 * t203;
t325 = t26 - t279;
t304 = t207 ^ 2;
t317 = t300 ^ 2;
t320 = t317 - t304;
t223 = t190 * t243;
t176 = t190 * qJD(3);
t258 = qJ(3) * t239 + qJD(1) * t176;
t88 = qJD(1) * t223 + t258;
t59 = t113 * pkin(4) + t88;
t12 = t36 * pkin(5) + pkin(9) * t199 + t59;
t164 = pkin(7) * t239;
t251 = qJD(4) * t190;
t252 = qJD(2) * t192;
t101 = t164 + (-qJ(4) * t252 - t251) * qJD(1);
t253 = qJD(2) * t190;
t282 = pkin(7) - qJ(4);
t117 = -qJD(4) * t192 - t282 * t253;
t181 = qJD(2) * qJD(3);
t95 = qJD(1) * t117 + t181;
t238 = -t186 * t101 + t185 * t95;
t211 = -pkin(8) * t114 - t238;
t51 = t185 * t101 + t186 * t95;
t43 = -pkin(8) * t113 + t51;
t208 = -t189 * t211 - t45 * t241 + t47 * t250 - t290 * t43;
t24 = -t246 * pkin(9) + t26;
t132 = -qJD(1) * pkin(1) - pkin(2) * t254 - qJ(3) * t255;
t99 = pkin(3) * t254 + qJD(4) - t132;
t69 = -pkin(4) * t121 + t99;
t27 = pkin(5) * t300 - pkin(9) * t207 + t69;
t8 = t188 * t27 + t191 * t24;
t2 = -qJD(6) * t8 + t191 * t12 + t188 * t208;
t319 = t309 * t8 + t2;
t214 = t188 * t24 - t191 * t27;
t1 = -qJD(6) * t214 + t12 * t188 - t191 * t208;
t315 = t214 * t309 + t1;
t316 = t207 * t8;
t313 = t309 * t207;
t312 = t207 * t214;
t284 = t207 * t300;
t311 = t246 * t300;
t297 = t207 * qJD(2);
t310 = -t227 - t297;
t237 = t189 * t43 - t290 * t211;
t307 = t207 * t69 + t237;
t306 = t300 * t69 + t208;
t305 = t199 + t311;
t37 = pkin(5) * t207 + pkin(9) * t300;
t30 = t189 * t203 + t290 * t56;
t206 = -t189 * t185 + t290 * t186;
t84 = t290 * t142 - t189 * t147;
t57 = qJD(3) * t206 + qJD(5) * t84;
t280 = t57 - t30;
t125 = t206 * qJD(2);
t126 = t206 * qJD(5);
t260 = t126 - t125;
t259 = t246 * t137;
t299 = qJD(1) * t212;
t296 = t121 ^ 2;
t295 = t123 ^ 2;
t194 = qJD(2) ^ 2;
t266 = t185 * t192;
t213 = -t186 * t190 + t266;
t153 = t282 * t190;
t154 = t282 * t192;
t91 = t186 * t153 - t154 * t185;
t65 = pkin(8) * t213 + t91;
t92 = t185 * t153 + t186 * t154;
t66 = -pkin(8) * t212 + t92;
t32 = t189 * t66 - t290 * t65;
t6 = qJD(5) * t26 + t237;
t294 = t32 * t6;
t79 = -t189 * t212 - t213 * t290;
t293 = t6 * t79;
t287 = t206 * t6;
t202 = t290 * t212;
t78 = -t189 * t213 + t202;
t286 = t36 * t78;
t285 = t55 * t53;
t283 = t79 * t36;
t278 = qJD(2) * pkin(2);
t275 = t188 * t53;
t274 = t188 * t55;
t272 = t191 * t22;
t271 = t191 * t53;
t270 = t191 * t55;
t268 = qJD(6) * t309;
t267 = t123 * t121;
t195 = qJD(1) ^ 2;
t264 = t192 * t195;
t263 = t194 * t190;
t177 = t194 * t192;
t118 = qJD(2) * t154 - t251;
t64 = t186 * t117 + t185 * t118;
t257 = qJ(3) * t252 + t176;
t183 = t190 ^ 2;
t256 = t192 ^ 2 - t183;
t150 = -t192 * pkin(2) - t190 * qJ(3) - pkin(1);
t245 = t79 * t249;
t244 = t79 * t248;
t236 = t309 ^ 2;
t231 = -0.2e1 * pkin(1) * t247;
t230 = qJD(3) - t278;
t229 = qJD(3) * t185 + t82;
t228 = qJD(3) * t186 - t83;
t134 = t192 * pkin(3) - t150;
t226 = qJD(1) * t150 + t132;
t224 = pkin(9) * t268 + t6;
t81 = -pkin(9) + t85;
t222 = t81 * t268 - t6;
t221 = t190 * t239;
t205 = t213 * qJD(2);
t41 = qJD(5) * t202 + t189 * t205 - t290 * t204 - t213 * t250;
t220 = -t23 * t41 + t293;
t219 = t188 * t8 - t191 * t214;
t218 = -t188 * t214 - t191 * t8;
t217 = -t309 * t41 + t283;
t215 = t185 * t61 - t186 * t62;
t86 = pkin(4) * t212 + t134;
t31 = t78 * pkin(5) - t79 * pkin(9) + t86;
t33 = t189 * t65 + t290 * t66;
t17 = -t188 * t33 + t191 * t31;
t18 = t188 * t31 + t191 * t33;
t63 = -t185 * t117 + t186 * t118;
t170 = qJ(3) * t254;
t116 = t193 * t255 + t170;
t102 = pkin(2) * t240 - t258;
t119 = pkin(2) * t253 - t257;
t210 = -pkin(7) * t194 - qJD(1) * t119 - t102;
t209 = -pkin(9) * t36 + t326;
t77 = -pkin(4) * t123 + t116;
t201 = -t309 * t57 - t36 * t81 - t326;
t200 = -qJD(6) * t219 + t1 * t191 - t188 * t2;
t198 = -pkin(8) * t204 + t63;
t68 = (pkin(4) * t266 + (-pkin(4) * t186 + t193) * t190) * qJD(2) + t257;
t148 = -pkin(7) * t240 + t181;
t149 = t172 + t230;
t152 = t173 + t182;
t196 = t148 * t192 + (t149 * t192 + (-t152 + t173) * t190) * qJD(2);
t163 = t190 * t264;
t156 = -0.2e1 * t221;
t155 = 0.2e1 * t221;
t151 = t256 * t195;
t144 = pkin(2) * t255 - t170;
t141 = t256 * t247;
t98 = t223 + t257;
t90 = t125 * t191 + t188 * t255;
t89 = -t125 * t188 + t191 * t255;
t80 = pkin(5) - t84;
t48 = -pkin(8) * t205 + t64;
t42 = t79 * qJD(5) + t189 * t204 + t290 * t205;
t28 = -t37 + t77;
t19 = t42 * pkin(5) + t41 * pkin(9) + t68;
t16 = t33 * qJD(5) + t189 * t48 - t290 * t198;
t15 = t189 * t198 + t65 * t241 - t66 * t250 + t290 * t48;
t14 = t188 * t37 + t191 * t25;
t13 = -t188 * t25 + t191 * t37;
t10 = t188 * t28 + t191 * t30;
t9 = -t188 * t30 + t191 * t28;
t4 = -qJD(6) * t18 - t15 * t188 + t19 * t191;
t3 = qJD(6) * t17 + t15 * t191 + t188 * t19;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, 0.2e1 * t141, t177, t156, -t263, 0, -pkin(7) * t177 + t190 * t231, pkin(7) * t263 + t192 * t231, 0, 0, t155, t177, -0.2e1 * t141, 0, t263, t156, t192 * t210 + t226 * t253, t196, t190 * t210 - t226 * t252, pkin(7) * t196 + t102 * t150 + t119 * t132, -t114 * t213 + t123 * t204, t213 * t113 - t114 * t212 + (t121 * t212 - t123 * t213) * qJD(2), -t212 * t194, t113 * t212 - t121 * t205, t213 * t194, 0, -t98 * t121 + t134 * t113 + t88 * t212 + (t213 * t99 - t63) * qJD(2), t134 * t114 + t98 * t123 - t88 * t213 + (t212 * t99 + t64) * qJD(2), t64 * t121 - t92 * t113 - t51 * t212 - t63 * t123 - t91 * t114 - t238 * t213 + (-t212 * t61 - t213 * t62) * qJD(2), t134 * t88 - t238 * t91 + t51 * t92 + t61 * t63 + t62 * t64 + t98 * t99, -t199 * t79 - t207 * t41, t199 * t78 - t207 * t42 + t300 * t41 - t283, t41 * t246, t300 * t42 + t286, t42 * t246, 0, t16 * t246 + t300 * t68 + t86 * t36 + t69 * t42 + t59 * t78, t15 * t246 - t86 * t199 + t207 * t68 - t69 * t41 + t59 * t79, -t15 * t300 + t16 * t207 - t199 * t32 + t208 * t78 + t25 * t41 - t26 * t42 - t33 * t36 + t293, t15 * t26 - t16 * t25 - t208 * t33 + t59 * t86 + t68 * t69 + t294, -t55 * t245 + (-t21 * t79 - t41 * t55) * t191 (t271 + t274) * t41 + (t277 - t272 + (-t270 + t275) * qJD(6)) * t79, t191 * t217 - t21 * t78 - t245 * t309 + t42 * t55, t53 * t244 + (t22 * t79 - t41 * t53) * t188, -t188 * t217 - t22 * t78 - t244 * t309 - t42 * t53, t309 * t42 + t286, t16 * t53 + t17 * t36 + t188 * t220 + t2 * t78 - t214 * t42 + t22 * t32 + t23 * t244 + t309 * t4, -t1 * t78 + t16 * t55 - t18 * t36 + t191 * t220 - t21 * t32 - t23 * t245 - t3 * t309 - t42 * t8, t17 * t21 - t18 * t22 - t3 * t53 - t4 * t55 + t219 * t41 + (qJD(6) * t218 - t1 * t188 - t191 * t2) * t79, t1 * t18 + t16 * t23 + t17 * t2 - t214 * t4 + t3 * t8 + t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t151, 0, t163, 0, 0, t195 * pkin(1) * t190, pkin(1) * t264, 0, 0, -t163, 0, t151, 0, 0, t163 (-t132 * t190 + t144 * t192) * qJD(1) ((t152 - t182) * t190 + (-t149 + t230) * t192) * qJD(1), 0.2e1 * t181 + (t132 * t192 + t144 * t190) * qJD(1), qJ(3) * t148 + qJD(3) * t152 - t132 * t144 + (t152 * t190 + (-t149 - t278) * t192) * qJD(1) * pkin(7), t267, -t295 + t296 (-t121 - t299) * qJD(2), -t267, -t157 + (t123 + t242) * qJD(2), 0, qJD(2) * t229 + t116 * t121 + t123 * t99 + t238, qJD(2) * t228 - t116 * t123 + t121 * t99 + t51, -t113 * t147 - t114 * t146 + (t229 - t62) * t123 + (t228 - t61) * t121, -qJD(3) * t215 - t116 * t99 - t146 * t238 + t147 * t51 - t61 * t82 - t62 * t83, -t284, t320, t305, t284, -t310, 0, t279 * qJD(2) + qJD(5) * t325 - t77 * t300 + t307, -t77 * t207 + t280 * t246 - t306, t84 * t199 - t85 * t36 + (t25 - t280) * t300 - t325 * t207, -t208 * t85 - t279 * t25 + t280 * t26 - t6 * t84 - t69 * t77, -t330, t329, -t328, -t275 * t309 + t272, t327, t313, t188 * t201 - t191 * t222 + t22 * t80 + t279 * t53 - t309 * t9 - t312, t10 * t309 + t188 * t222 + t191 * t201 - t21 * t80 + t279 * t55 - t316, t10 * t53 + t55 * t9 + (-t22 * t81 - t53 * t57 - t214 * t300 - t1 + (t55 * t81 - t214) * qJD(6)) * t191 + (-t21 * t81 + t55 * t57 + t300 * t8 + t2 + (t53 * t81 + t8) * qJD(6)) * t188, -t10 * t8 + t200 * t81 + t214 * t9 - t218 * t57 + t23 * t279 + t6 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, 0, -t183 * t195 - t194, -qJD(2) * t152 + t132 * t255 + t164, 0, 0, 0, 0, 0, 0, t121 * t255 - t185 * t194, -t123 * t255 - t186 * t194, -t113 * t185 - t114 * t186 + (-t121 * t186 - t123 * t185) * qJD(2), qJD(2) * t215 + t185 * t51 - t186 * t238 - t99 * t255, 0, 0, 0, 0, 0, 0, -t259 * t246 - t255 * t300, -t207 * t255 + t260 * t246, -t137 * t36 + t199 * t206 - t207 * t259 - t260 * t300, -t137 * t208 + t259 * t25 - t69 * t255 + t260 * t26 - t287, 0, 0, 0, 0, 0, 0, -t137 * t276 - t206 * t22 - t259 * t53 + (-t126 * t188 - t137 * t248 - t89) * t309, -t137 * t35 + t206 * t21 - t259 * t55 + (-t126 * t191 + t137 * t249 + t90) * t309, t53 * t90 + t55 * t89 + (-t271 + t274) * t126 + (-t277 - t272 + (t270 + t275) * qJD(6)) * t137, -t126 * t218 + t137 * t200 + t214 * t89 - t23 * t259 - t8 * t90 - t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157 + (-t123 + t242) * qJD(2) (-t121 + t299) * qJD(2), -t295 - t296, -t121 * t62 + t123 * t61 + t88, 0, 0, 0, 0, 0, 0, t227 - t297 + 0.2e1 * t298, -t199 + t311, -t304 - t317, t207 * t25 + t26 * t300 + t59, 0, 0, 0, 0, 0, 0, -t188 * t236 - t302 + t35, -t191 * t236 - t276 - t301 (t21 - t314) * t191 + t309 * t274 + t281, t188 * t315 + t191 * t319 - t23 * t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, -t320, -t305, -t284, t310, 0, -t26 * qJD(2) - t307, -t25 * t246 + t306, 0, 0, t330, -t329, t328, t323 * t53 - t272, -t327, -t313, -pkin(5) * t22 - t13 * t309 + t188 * t209 - t191 * t224 - t26 * t53 + t312, pkin(5) * t21 + t14 * t309 + t188 * t224 + t191 * t209 - t26 * t55 + t316, t13 * t55 + t14 * t53 + ((-t22 + t269) * pkin(9) + t315) * t191 + ((qJD(6) * t53 - t21) * pkin(9) - t319) * t188, -pkin(5) * t6 + pkin(9) * t200 + t13 * t214 - t14 * t8 - t23 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t53 ^ 2 + t55 ^ 2, t309 * t53 - t21, -t285, t309 * t55 - t22, t36, -t23 * t55 + t319, t23 * t53 - t315, 0, 0;];
tauc_reg  = t5;