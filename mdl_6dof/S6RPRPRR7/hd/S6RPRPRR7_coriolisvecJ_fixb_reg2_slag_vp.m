% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:04
% EndTime: 2019-03-09 03:56:17
% DurationCPUTime: 5.02s
% Computational Cost: add. (10531->414), mult. (22804->536), div. (0->0), fcn. (16414->8), ass. (0->227)
t175 = sin(qJ(5));
t173 = sin(pkin(10));
t176 = sin(qJ(3));
t178 = cos(qJ(3));
t248 = cos(pkin(10));
t291 = t173 * t178 + t248 * t176;
t289 = qJD(3) * t291;
t121 = qJD(1) * t289;
t179 = -pkin(1) - pkin(7);
t150 = t179 * qJD(1) + qJD(2);
t230 = qJD(4) * t178;
t233 = qJD(3) * t176;
t95 = -t150 * t233 + (qJ(4) * t233 - t230) * qJD(1);
t231 = qJD(4) * t176;
t232 = qJD(3) * t178;
t96 = t150 * t232 + (-qJ(4) * t232 - t231) * qJD(1);
t60 = -t173 * t96 + t248 * t95;
t194 = pkin(8) * t121 + t60;
t272 = cos(qJ(5));
t216 = qJD(5) * t272;
t229 = qJD(5) * t175;
t212 = qJD(1) * t248;
t202 = t178 * t212;
t148 = qJD(3) * t202;
t235 = qJD(1) * t176;
t217 = t173 * t235;
t120 = qJD(3) * t217 - t148;
t61 = t173 * t95 + t248 * t96;
t49 = pkin(8) * t120 + t61;
t136 = t202 - t217;
t270 = pkin(8) * t136;
t127 = -qJ(4) * t235 + t150 * t176;
t106 = t173 * t127;
t234 = qJD(1) * t178;
t128 = -qJ(4) * t234 + t178 * t150;
t109 = qJD(3) * pkin(3) + t128;
t69 = t248 * t109 - t106;
t55 = qJD(3) * pkin(4) - t270 + t69;
t135 = -t173 * t234 - t176 * t212;
t268 = t135 * pkin(8);
t107 = t248 * t127;
t70 = t173 * t109 + t107;
t58 = t70 + t268;
t182 = -t175 * t194 - t55 * t216 + t58 * t229 - t272 * t49;
t213 = t248 * t178;
t141 = -t173 * t176 + t213;
t199 = -qJD(3) * t213 + t173 * t233;
t292 = t175 * t291;
t298 = -qJD(5) * t292 + t141 * t216 - t175 * t289 - t272 * t199;
t200 = t272 * t291;
t303 = t175 * t141 + t200;
t214 = t175 * t49 - t272 * t194;
t33 = t175 * t55 + t272 * t58;
t8 = qJD(5) * t33 + t214;
t84 = -t272 * t141 + t292;
t329 = t8 * t84;
t330 = -t182 * t303 + t298 * t33 + t329;
t177 = cos(qJ(6));
t193 = t175 * t135 + t272 * t136;
t224 = qJD(3) + qJD(5);
t205 = t177 * t224;
t174 = sin(qJ(6));
t228 = qJD(6) * t174;
t44 = -t175 * t120 + t272 * t121 - t135 * t216 + t136 * t229;
t30 = -qJD(6) * t205 + t177 * t44 + t193 * t228;
t328 = t30 * t84;
t67 = t174 * t224 + t177 * t193;
t247 = qJD(6) * t67;
t31 = -t174 * t44 + t247;
t327 = t31 * t84;
t326 = t44 * t84;
t207 = t272 * t120 + t175 * t121;
t288 = qJD(5) * t193;
t45 = -t207 + t288;
t169 = qJD(1) * qJD(2);
t225 = qJD(1) * qJD(3);
t215 = t178 * t225;
t143 = pkin(3) * t215 + t169;
t88 = -pkin(4) * t120 + t143;
t18 = pkin(5) * t45 + pkin(9) * t44 + t88;
t29 = pkin(9) * t224 + t33;
t80 = -t272 * t135 + t136 * t175;
t144 = pkin(3) * t235 + qJD(1) * qJ(2) + qJD(4);
t94 = -pkin(4) * t135 + t144;
t36 = pkin(5) * t80 - pkin(9) * t193 + t94;
t196 = t174 * t29 - t177 * t36;
t301 = qJD(6) * t196;
t2 = t174 * t18 - t177 * t182 - t301;
t304 = qJD(6) + t80;
t276 = t304 * t196;
t325 = t2 + t276;
t318 = t174 * t304;
t324 = t67 * t318;
t267 = t45 * t303;
t323 = t298 * t80 + t267;
t27 = t31 * t177;
t65 = t174 * t193 - t205;
t322 = t318 * t65 - t27;
t10 = t174 * t36 + t177 * t29;
t321 = t10 * t298;
t295 = t224 * t80;
t319 = -t44 + t295;
t317 = t298 * t224;
t261 = t80 ^ 2;
t262 = t193 ^ 2;
t316 = -t261 + t262;
t3 = -qJD(6) * t10 + t174 * t182 + t177 * t18;
t315 = t10 * t304 + t3;
t227 = qJD(6) * t177;
t259 = -t174 * t31 - t65 * t227;
t306 = t177 * t80;
t314 = -t177 * t30 - t306 * t65 + t259;
t25 = t30 * t174;
t313 = -t25 + (t227 + t306) * t67;
t41 = t174 * t45;
t258 = t227 * t304 + t41;
t263 = t67 * t193;
t312 = t304 * t306 + t258 - t263;
t311 = t196 * t298 - t3 * t303;
t310 = t298 * t67 - t30 * t303;
t309 = -t298 * t65 - t303 * t31;
t32 = -t175 * t58 + t272 * t55;
t28 = -pkin(5) * t224 - t32;
t307 = t28 * t80;
t260 = t80 * t193;
t161 = t248 * pkin(3) + pkin(4);
t271 = pkin(3) * t173;
t131 = t175 * t161 + t272 * t271;
t76 = -t173 * t128 - t107;
t187 = t76 - t268;
t77 = t248 * t128 - t106;
t62 = t77 - t270;
t249 = qJD(5) * t131 - t175 * t62 + t272 * t187;
t305 = -t33 - t249;
t48 = pkin(5) * t193 + pkin(9) * t80;
t300 = t94 * t80 + t182;
t299 = t60 * t141 - t199 * t70 - t289 * t69 + t291 * t61;
t265 = t65 * t193;
t294 = t304 * t193;
t130 = t272 * t161 - t175 * t271;
t113 = t130 * qJD(5);
t35 = t175 * t187 + t272 * t62;
t250 = t113 - t35;
t43 = t177 * t45;
t293 = t228 * t304 - t43;
t241 = t193 * qJD(3);
t290 = t241 + t207;
t287 = t199 * qJD(3);
t22 = t28 * t228;
t286 = t193 * t196 + t22;
t23 = t28 * t227;
t285 = t10 * t193 + t8 * t174 + t23;
t284 = -t94 * t193 - t214;
t283 = -t175 * t199 + t272 * t289;
t281 = -t120 * t291 + t135 * t199;
t279 = t136 ^ 2;
t278 = 0.2e1 * t169;
t180 = qJD(3) ^ 2;
t240 = qJ(4) - t179;
t145 = t240 * t176;
t146 = t240 * t178;
t92 = t145 * t173 - t248 * t146;
t74 = -pkin(8) * t141 + t92;
t93 = -t248 * t145 - t173 * t146;
t75 = -pkin(8) * t291 + t93;
t38 = t175 * t75 - t272 * t74;
t277 = t38 * t8;
t1 = t2 * t177;
t266 = t45 * t84;
t264 = t67 * t65;
t255 = t174 * t67;
t246 = t136 * t135;
t244 = t180 * t176;
t243 = t180 * t178;
t181 = qJD(1) ^ 2;
t242 = t181 * qJ(2);
t162 = t176 * pkin(3) + qJ(2);
t122 = t233 * t240 - t230;
t123 = -qJD(3) * t146 - t231;
t73 = t173 * t122 + t248 * t123;
t237 = t176 ^ 2 - t178 ^ 2;
t236 = -t180 - t181;
t226 = t144 * qJD(1);
t151 = pkin(3) * t232 + qJD(2);
t222 = 0.2e1 * qJD(1);
t165 = pkin(3) * t234;
t220 = t84 * t228;
t219 = t84 * t227;
t218 = t178 * t181 * t176;
t100 = pkin(4) * t136 + t165;
t206 = qJD(6) * t303 + qJD(1);
t204 = t176 * t215;
t51 = qJD(5) * t200 + t141 * t229 + t283;
t203 = -t28 * t51 - t329;
t201 = -t304 * t51 - t266;
t198 = t10 * t174 - t177 * t196;
t39 = t175 * t74 + t272 * t75;
t110 = pkin(4) * t291 + t162;
t40 = pkin(5) * t303 + pkin(9) * t84 + t110;
t20 = t174 * t40 + t177 * t39;
t19 = -t174 * t39 + t177 * t40;
t195 = -t318 * t80 - t293;
t72 = t248 * t122 - t173 * t123;
t125 = pkin(9) + t131;
t191 = -t113 * t304 - t125 * t45 + t307;
t188 = t291 * qJD(1);
t97 = -pkin(4) * t199 + t151;
t186 = -qJD(6) * t198 - t174 * t3 + t1;
t185 = -t141 * t121 - t136 * t289;
t183 = pkin(8) * t289 + t72;
t167 = qJ(2) * t278;
t132 = t135 ^ 2;
t126 = t291 * t180;
t124 = -pkin(5) - t130;
t59 = pkin(8) * t199 + t73;
t52 = qJD(5) * t303 + t283;
t37 = t100 + t48;
t21 = pkin(5) * t298 + t51 * pkin(9) + t97;
t16 = t174 * t48 + t177 * t32;
t15 = -t174 * t32 + t177 * t48;
t14 = t39 * qJD(5) + t175 * t59 - t272 * t183;
t13 = t175 * t183 + t74 * t216 - t75 * t229 + t272 * t59;
t12 = t174 * t37 + t177 * t35;
t11 = -t174 * t35 + t177 * t37;
t5 = -qJD(6) * t20 - t13 * t174 + t177 * t21;
t4 = qJD(6) * t19 + t13 * t177 + t174 * t21;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, t167, -0.2e1 * t204, 0.2e1 * t237 * t225, -t244, 0.2e1 * t204, -t243, 0, -t179 * t244 + (qJ(2) * t232 + qJD(2) * t176) * t222, -t179 * t243 + (-qJ(2) * t233 + qJD(2) * t178) * t222, 0, t167, t185, t141 * t120 + t121 * t291 - t135 * t289 + t136 * t199, -t126, t281, t287, 0, t72 * qJD(3) - t162 * t120 - t151 * t135 + t143 * t291 - t144 * t199, -t73 * qJD(3) - t162 * t121 + t151 * t136 + t143 * t141 - t144 * t289, t93 * t120 + t92 * t121 + t73 * t135 - t72 * t136 - t299, t143 * t162 + t144 * t151 + t60 * t92 + t61 * t93 + t69 * t72 + t70 * t73, -t193 * t51 + t326, -t193 * t298 + t303 * t44 + t51 * t80 + t266, -t51 * t224, t323, -t317, 0, t110 * t45 - t14 * t224 + t298 * t94 + t303 * t88 + t97 * t80, -t110 * t44 - t13 * t224 + t193 * t97 - t94 * t51 - t84 * t88, -t13 * t80 + t14 * t193 + t32 * t51 - t38 * t44 - t39 * t45 - t330, t110 * t88 + t13 * t33 - t14 * t32 - t182 * t39 + t94 * t97 + t277, t67 * t220 + (-t51 * t67 + t328) * t177 (t177 * t65 + t255) * t51 - (t25 - t27 + (t174 * t65 - t177 * t67) * qJD(6)) * t84, t177 * t201 + t220 * t304 + t310, -t65 * t219 + (-t51 * t65 - t327) * t174, -t174 * t201 + t219 * t304 + t309, t298 * t304 + t267, t14 * t65 + t174 * t203 + t19 * t45 - t23 * t84 + t304 * t5 + t31 * t38 - t311, t14 * t67 + t177 * t203 - t2 * t303 - t20 * t45 + t22 * t84 - t30 * t38 - t304 * t4 - t321, t19 * t30 - t20 * t31 - t4 * t65 - t5 * t67 + t198 * t51 - (-t174 * t2 - t177 * t3 + (-t10 * t177 - t174 * t196) * qJD(6)) * t84, t10 * t4 + t14 * t28 + t19 * t3 - t196 * t5 + t2 * t20 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t242, 0, 0, 0, 0, 0, 0, t236 * t176, t236 * t178, 0, -t242, 0, 0, 0, 0, 0, 0, qJD(1) * t135 - t126, -qJD(1) * t136 + t287, -t185 - t281, -t226 + t299, 0, 0, 0, 0, 0, 0, -qJD(1) * t80 - t224 * t52, -qJD(1) * t193 - t317, t193 * t52 - t323 - t326, -qJD(1) * t94 - t32 * t52 + t330, 0, 0, 0, 0, 0, 0, -t303 * t41 + t327 + t52 * t65 + (-t174 * t298 - t177 * t206) * t304, -t303 * t43 - t328 + t52 * t67 + (t174 * t206 - t177 * t298) * t304 (t206 * t67 + t309) * t177 + (t206 * t65 + t310) * t174, t28 * t52 + t329 + (qJD(1) * t196 + t321 + (t2 + t301) * t303) * t177 + (-t10 * t206 + t311) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t237 * t181, 0, -t218, 0, 0, -t178 * t242, t176 * t242, 0, 0, -t246, -t132 + t279 (-t135 - t188) * qJD(3), t246, -t148 + (t136 + t217) * qJD(3), 0, -qJD(3) * t76 + t135 * t165 - t136 * t144 + t60, qJD(3) * t77 - t135 * t144 - t136 * t165 - t61 (t70 + t76) * t136 + (t69 - t77) * t135 + (t120 * t173 + t121 * t248) * pkin(3), -t69 * t76 - t70 * t77 + (t173 * t61 - t178 * t226 + t248 * t60) * pkin(3), t260, t316, t319, -t260, t290, 0, -t249 * qJD(3) + t305 * qJD(5) - t100 * t80 + t284, -t100 * t193 - t224 * t250 + t300, t130 * t44 - t131 * t45 - t305 * t193 + (-t250 - t32) * t80, -t100 * t94 - t130 * t8 - t131 * t182 - t249 * t32 + t250 * t33, t313, -t255 * t304 + t314, t312, t322, t195 + t265, -t294, -t11 * t304 + t124 * t31 + t249 * t65 + (-qJD(6) * t125 * t304 - t8) * t177 + t191 * t174 + t286, -t124 * t30 + (t125 * t228 + t12) * t304 + t249 * t67 + t191 * t177 + t285, t11 * t67 + t12 * t65 + t1 + (-t113 * t65 - t125 * t31 + t80 * t196 + (t125 * t67 + t196) * qJD(6)) * t177 + (-t10 * t80 + t113 * t67 - t125 * t30 - t3 + (t125 * t65 - t10) * qJD(6)) * t174, t124 * t8 - (-t113 * t174 - t11) * t196 + t249 * t28 + (t113 * t177 - t12) * t10 + t186 * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148 + (t136 - t217) * qJD(3) (t135 - t188) * qJD(3), -t132 - t279, -t135 * t70 + t136 * t69 + t143, 0, 0, 0, 0, 0, 0, -t207 + t241 + 0.2e1 * t288, -t44 - t295, -t261 - t262, t193 * t32 + t33 * t80 + t88, 0, 0, 0, 0, 0, 0, t195 - t265, -t177 * t304 ^ 2 - t263 - t41 (-t65 * t80 + t30) * t177 + t324 + t259, t174 * t325 + t315 * t177 - t28 * t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t316, t319, -t260, t290, 0, t33 * qJD(3) + t284, t224 * t32 + t300, 0, 0, t313, t314 - t324, t312, t322, -t304 * t318 + t265 + t43, -t294, -pkin(5) * t31 - pkin(9) * t258 - t15 * t304 + t174 * t307 - t177 * t8 - t33 * t65 + t286, pkin(5) * t30 + pkin(9) * t293 + t16 * t304 + t28 * t306 - t33 * t67 + t285, t15 * t67 + t16 * t65 + t1 + (t276 + (-t31 + t247) * pkin(9)) * t177 + ((qJD(6) * t65 - t30) * pkin(9) - t315) * t174, -pkin(5) * t8 + pkin(9) * t186 - t10 * t16 + t15 * t196 - t28 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, -t65 ^ 2 + t67 ^ 2, t304 * t65 - t30, -t264, t304 * t67 - t31, t45, -t28 * t67 + t315, t28 * t65 - t325, 0, 0;];
tauc_reg  = t6;
