% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:49
% EndTime: 2019-12-05 17:17:03
% DurationCPUTime: 5.41s
% Computational Cost: add. (5652->488), mult. (13241->680), div. (0->0), fcn. (10297->14), ass. (0->248)
t166 = sin(qJ(2));
t161 = sin(pkin(5));
t258 = qJD(1) * t161;
t226 = t166 * t258;
t165 = sin(qJ(3));
t254 = qJD(3) * t165;
t182 = pkin(3) * t254 - t226;
t312 = pkin(8) + pkin(7);
t231 = qJD(3) * t312;
t126 = t165 * t231;
t164 = sin(qJ(4));
t168 = cos(qJ(3));
t309 = cos(qJ(4));
t230 = t309 * t168;
t186 = -t164 * t165 + t230;
t133 = t312 * t165;
t134 = t312 * t168;
t187 = -t309 * t133 - t164 * t134;
t169 = cos(qJ(2));
t225 = t169 * t258;
t267 = t164 * t168;
t300 = t187 * qJD(4) - t309 * t126 - t186 * t225 - t231 * t267;
t207 = qJD(3) * t230;
t243 = qJD(3) + qJD(4);
t210 = t164 * t243;
t223 = qJD(4) * t309;
t83 = t165 * t210 - t168 * t223 - t207;
t123 = t309 * t165 + t267;
t84 = t243 * t123;
t321 = t84 * pkin(4) + t83 * pkin(9) + t182;
t162 = cos(pkin(5));
t272 = t161 * t166;
t116 = t162 * t168 - t165 * t272;
t159 = qJ(3) + qJ(4);
t154 = sin(t159);
t155 = cos(t159);
t105 = -t154 * t272 + t162 * t155;
t283 = cos(pkin(10));
t215 = t283 * t166;
t160 = sin(pkin(10));
t273 = t160 * t169;
t113 = t162 * t215 + t273;
t216 = t161 * t283;
t75 = -t113 * t154 - t155 * t216;
t214 = t283 * t169;
t274 = t160 * t166;
t115 = -t162 * t274 + t214;
t275 = t160 * t161;
t77 = -t115 * t154 + t155 * t275;
t320 = g(1) * t77 + g(2) * t75 + g(3) * t105;
t163 = sin(qJ(5));
t167 = cos(qJ(5));
t250 = qJD(5) * t167;
t127 = qJD(2) * pkin(7) + t226;
t213 = pkin(8) * qJD(2) + t127;
t257 = qJD(1) * t162;
t224 = t165 * t257;
t88 = t213 * t168 + t224;
t291 = t164 * t88;
t140 = t168 * t257;
t87 = -t213 * t165 + t140;
t82 = qJD(3) * pkin(3) + t87;
t44 = t309 * t82 - t291;
t38 = -t243 * pkin(4) - t44;
t242 = qJDD(3) + qJDD(4);
t248 = qJD(1) * qJD(2);
t221 = t169 * t248;
t281 = qJDD(2) * pkin(7);
t104 = t281 + (qJDD(1) * t166 + t221) * t161;
t246 = qJDD(1) * t162;
t139 = t168 * t246;
t37 = qJDD(3) * pkin(3) + t139 + (-pkin(8) * qJDD(2) - t104) * t165 - t88 * qJD(3);
t247 = qJD(2) * qJD(3);
t220 = t165 * t247;
t244 = t168 * qJDD(2);
t237 = qJD(3) * t140 + t168 * t104 + t165 * t246;
t52 = -t127 * t254 + t237;
t40 = (-t220 + t244) * pkin(8) + t52;
t217 = t164 * t40 - t309 * t37;
t236 = t309 * t88;
t45 = t164 * t82 + t236;
t9 = -t45 * qJD(4) - t217;
t7 = -t242 * pkin(4) - t9;
t319 = t7 * t163 + t38 * t250;
t95 = -t164 * t133 + t309 * t134;
t299 = t95 * qJD(4) - t123 * t225 - t164 * t126 + t312 * t207;
t253 = qJD(4) * t164;
t48 = t164 * t87 + t236;
t205 = pkin(3) * t253 - t48;
t112 = -t162 * t214 + t274;
t114 = t162 * t273 + t215;
t202 = g(1) * t114 + g(2) * t112;
t270 = t161 * t169;
t180 = g(3) * t270 - t202;
t318 = t180 * t154;
t317 = t123 * qJDD(2);
t39 = t243 * pkin(9) + t45;
t153 = t168 * pkin(3) + pkin(2);
t110 = -t153 * qJD(2) - t225;
t256 = qJD(2) * t165;
t118 = -qJD(2) * t230 + t164 * t256;
t120 = t123 * qJD(2);
t58 = t118 * pkin(4) - t120 * pkin(9) + t110;
t195 = t163 * t39 - t167 * t58;
t313 = t243 * qJD(2);
t59 = -t186 * t313 - t317;
t245 = t165 * qJDD(2);
t198 = -qJDD(2) * t230 + t164 * t245;
t60 = t84 * qJD(2) + t198;
t222 = t166 * t248;
t199 = -qJDD(1) * t270 + t161 * t222;
t74 = pkin(3) * t220 - t153 * qJDD(2) + t199;
t21 = t60 * pkin(4) + t59 * pkin(9) + t74;
t8 = t164 * t37 + t82 * t223 - t88 * t253 + t309 * t40;
t6 = t242 * pkin(9) + t8;
t2 = -t195 * qJD(5) + t163 * t21 + t167 * t6;
t1 = t2 * t167;
t279 = t118 * t167;
t23 = t163 * t58 + t167 * t39;
t295 = t163 * t23;
t316 = -t118 * t295 + t195 * t279 + t1;
t306 = t164 * pkin(3);
t151 = pkin(9) + t306;
t208 = pkin(3) * t223;
t49 = t309 * t87 - t291;
t239 = pkin(3) * t256;
t72 = t120 * pkin(4) + t118 * pkin(9);
t64 = t72 + t239;
t25 = t163 * t64 + t167 * t49;
t251 = qJD(5) * t163;
t315 = -t151 * t251 + t167 * t208 - t25;
t93 = t167 * t120 + t163 * t243;
t31 = t93 * qJD(5) - t163 * t59 - t167 * t242;
t171 = qJD(3) ^ 2;
t282 = qJDD(2) * pkin(2);
t103 = t199 - t282;
t191 = -t103 + t202;
t314 = -pkin(7) * t171 + t161 * (-g(3) * t169 + t222) + t191 + t282;
t307 = g(3) * t161;
t3 = -qJD(5) * t23 - t163 * t6 + t167 * t21;
t305 = t3 * t163;
t209 = t167 * t243;
t91 = t163 * t120 - t209;
t304 = t93 * t91;
t73 = -pkin(4) * t186 - t123 * pkin(9) - t153;
t46 = -t163 * t95 + t167 * t73;
t303 = t46 * qJD(5) + t321 * t163 + t300 * t167;
t47 = t163 * t73 + t167 * t95;
t302 = -t47 * qJD(5) - t300 * t163 + t321 * t167;
t301 = t45 - t48;
t298 = qJD(2) * pkin(2);
t297 = t118 * t38;
t296 = t163 * t195;
t55 = qJDD(5) + t60;
t294 = t163 * t55;
t293 = t163 * t83;
t292 = t163 * t91;
t290 = t167 * t55;
t289 = t167 * t83;
t288 = t167 * t93;
t30 = -qJD(5) * t209 + t120 * t251 - t163 * t242 + t167 * t59;
t287 = t30 * t163;
t286 = t31 * t167;
t111 = qJD(5) + t118;
t285 = t91 * t111;
t284 = t93 * t111;
t280 = t111 * t120;
t278 = t120 * t118;
t277 = t155 * t163;
t276 = t155 * t167;
t271 = t161 * t168;
t268 = t163 * t169;
t266 = t165 * t127;
t265 = t166 * t312;
t264 = t167 * t169;
t263 = qJDD(1) - g(3);
t262 = -t112 * t153 + t113 * t312;
t261 = -t114 * t153 + t115 * t312;
t157 = t165 ^ 2;
t158 = t168 ^ 2;
t260 = t157 - t158;
t259 = t157 + t158;
t255 = qJD(2) * t166;
t252 = qJD(5) * t111;
t128 = -t225 - t298;
t249 = t128 * qJD(2);
t241 = t309 * pkin(3);
t234 = t161 * t268;
t233 = t161 * t264;
t172 = qJD(2) ^ 2;
t232 = t165 * t172 * t168;
t32 = t38 * t251;
t229 = t161 * t255;
t228 = qJD(2) * t270;
t219 = t169 * t247;
t218 = t120 * t195 + t32;
t211 = t111 * t167;
t206 = t168 * t220;
t204 = (-t115 * t165 + t160 * t271) * pkin(3);
t203 = pkin(4) * t155 + pkin(9) * t154;
t201 = g(1) * t115 + g(2) * t113;
t197 = -t151 * t55 + t297;
t196 = -t167 * t195 + t295;
t193 = t116 * pkin(3);
t192 = qJDD(2) * t169 - t166 * t172;
t117 = t162 * t165 + t166 * t271;
t66 = t164 * t116 + t309 * t117;
t56 = -t163 * t66 - t233;
t190 = -t167 * t66 + t234;
t188 = t309 * t116 - t164 * t117;
t99 = t168 * t127 + t224;
t106 = t162 * t154 + t155 * t272;
t76 = t113 * t155 - t154 * t216;
t78 = t115 * t155 + t154 * t275;
t184 = g(1) * t78 + g(2) * t76 + g(3) * t106;
t183 = t23 * t120 + t320 * t163 + t319;
t181 = -t320 - t7;
t179 = (-t113 * t165 - t168 * t216) * pkin(3);
t177 = -pkin(7) * qJDD(3) + (t128 + t225 - t298) * qJD(3);
t176 = -t110 * t120 - t217 - t320;
t175 = t110 * t118 + t184 - t8;
t174 = -t196 * qJD(5) - t184 - t305;
t53 = -t99 * qJD(3) - t165 * t104 + t139;
t98 = t140 - t266;
t173 = -t53 * t165 + t52 * t168 + (-t165 * t99 - t168 * t98) * qJD(3) - t201;
t152 = -t241 - pkin(4);
t125 = t153 * t270;
t102 = t105 * pkin(4);
t86 = -t117 * qJD(3) - t165 * t228;
t85 = t116 * qJD(3) + t168 * t228;
t71 = t77 * pkin(4);
t70 = t75 * pkin(4);
t61 = -t118 ^ 2 + t120 ^ 2;
t42 = t120 * t243 - t123 * t313 - t198;
t41 = t317 + (qJD(2) * t186 + t118) * t243;
t29 = t66 * qJD(4) + t164 * t85 - t309 * t86;
t28 = t188 * qJD(4) + t164 * t86 + t309 * t85;
t27 = t163 * t72 + t167 * t44;
t26 = -t163 * t44 + t167 * t72;
t24 = -t163 * t49 + t167 * t64;
t17 = t111 * t211 - t93 * t120 + t294;
t16 = -t111 ^ 2 * t163 + t91 * t120 + t290;
t15 = t163 * t285 - t286;
t14 = t93 * t211 - t287;
t13 = t190 * qJD(5) - t163 * t28 + t167 * t229;
t12 = t56 * qJD(5) + t163 * t229 + t167 * t28;
t4 = (-t30 - t285) * t167 + (-t31 - t284) * t163;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t263, 0, 0, 0, 0, 0, 0, t192 * t161, (-qJDD(2) * t166 - t169 * t172) * t161, 0, -g(3) + (t162 ^ 2 + (t166 ^ 2 + t169 ^ 2) * t161 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t86 * qJD(3) + t116 * qJDD(3) + (-t165 * t219 + t168 * t192) * t161, -t85 * qJD(3) - t117 * qJDD(3) + (-t165 * t192 - t168 * t219) * t161, (-t116 * t165 + t117 * t168) * qJDD(2) + (-t165 * t86 + t168 * t85 + (-t116 * t168 - t117 * t165) * qJD(3)) * qJD(2), t53 * t116 + t52 * t117 + t99 * t85 + t98 * t86 - g(3) + (-t103 * t169 + t166 * t249) * t161, 0, 0, 0, 0, 0, 0, -t29 * t243 + t188 * t242 + (t118 * t255 - t169 * t60) * t161, -t28 * t243 - t66 * t242 + (t120 * t255 + t169 * t59) * t161, -t118 * t28 + t120 * t29 + t188 * t59 - t60 * t66, t45 * t28 - t44 * t29 + t9 * t188 + t8 * t66 - g(3) + (t110 * t255 - t169 * t74) * t161, 0, 0, 0, 0, 0, 0, t111 * t13 - t188 * t31 + t29 * t91 + t56 * t55, -t111 * t12 + t188 * t30 + t190 * t55 + t29 * t93, -t12 * t91 - t13 * t93 + t190 * t31 + t30 * t56, t12 * t23 - t13 * t195 - t188 * t7 - t190 * t2 + t29 * t38 + t3 * t56 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t263 * t270 + t202, -t263 * t272 + t201, 0, 0, t157 * qJDD(2) + 0.2e1 * t206, 0.2e1 * t165 * t244 - 0.2e1 * t260 * t247, t165 * qJDD(3) + t171 * t168, t158 * qJDD(2) - 0.2e1 * t206, qJDD(3) * t168 - t171 * t165, 0, t177 * t165 + t314 * t168, -t314 * t165 + t177 * t168, t259 * t281 + (-g(3) * t166 - t259 * t221) * t161 + t173, t191 * pkin(2) + t173 * pkin(7) + (-g(3) * (pkin(2) * t169 + pkin(7) * t166) + (-t128 * t166 + (t165 * t98 - t168 * t99) * t169) * qJD(1)) * t161, -t120 * t83 - t123 * t59, t118 * t83 - t120 * t84 - t123 * t60 - t186 * t59, t123 * t242 - t83 * t243, t118 * t84 - t186 * t60, t186 * t242 - t84 * t243, 0, t110 * t84 + t118 * t182 - t153 * t60 - t155 * t180 - t186 * t74 + t187 * t242 - t299 * t243, -t110 * t83 + t120 * t182 + t74 * t123 + t153 * t59 - t242 * t95 - t300 * t243 + t318, -g(3) * t272 - t300 * t118 + t299 * t120 - t9 * t123 + t186 * t8 + t187 * t59 + t44 * t83 - t45 * t84 - t95 * t60 - t201, t8 * t95 + t9 * t187 - t74 * t153 - g(1) * t261 - g(2) * t262 - g(3) * (t161 * t265 + t125) + t300 * t45 - t299 * t44 + t182 * t110, -t83 * t288 + (-t30 * t167 - t251 * t93) * t123, (t163 * t93 + t167 * t91) * t83 + (t287 - t286 + (-t288 + t292) * qJD(5)) * t123, t123 * t290 + t30 * t186 + t93 * t84 + (-t123 * t251 - t289) * t111, -t83 * t292 + (t163 * t31 + t250 * t91) * t123, -t123 * t294 + t31 * t186 - t91 * t84 + (-t123 * t250 + t293) * t111, t111 * t84 - t186 * t55, t46 * t55 - t3 * t186 - t195 * t84 - t187 * t31 - t38 * t293 - g(1) * (-t114 * t276 + t115 * t163) - g(2) * (-t112 * t276 + t113 * t163) + t299 * t91 - (t155 * t264 + t163 * t166) * t307 + t319 * t123 + t302 * t111, -t47 * t55 + t2 * t186 - t23 * t84 + t187 * t30 - t38 * t289 - g(1) * (t114 * t277 + t115 * t167) - g(2) * (t112 * t277 + t113 * t167) + t299 * t93 - (-t155 * t268 + t166 * t167) * t307 + (t7 * t167 - t32) * t123 - t303 * t111, t46 * t30 - t47 * t31 - t302 * t93 - t303 * t91 + t196 * t83 - t318 + (-t163 * t2 - t167 * t3 + (-t167 * t23 - t296) * qJD(5)) * t123, t2 * t47 + t3 * t46 - t7 * t187 - g(1) * (-t114 * t203 + t261) - g(2) * (-t112 * t203 + t262) - g(3) * t125 + t299 * t38 + t303 * t23 - t302 * t195 - (t169 * t203 + t265) * t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, t260 * t172, t245, t232, t244, qJDD(3), -g(3) * t116 + t139 + (-g(1) * t160 + g(2) * t283) * t271 + (-t104 + t201 - t249) * t165, -t168 * t249 - g(1) * (-t115 * t168 - t165 * t275) - g(2) * (-t113 * t168 + t165 * t216) + g(3) * t117 + (t98 + t266) * qJD(3) - t237, 0, 0, t278, t61, t41, -t278, t42, t242, t48 * qJD(3) - t301 * qJD(4) + (-qJD(4) * t210 - t118 * t256 + t309 * t242) * pkin(3) + t176, t49 * t243 + (-t120 * t256 - t164 * t242 - t223 * t243) * pkin(3) + t175, t301 * t120 + (-t44 + t49) * t118 + (t309 * t59 - t164 * t60 + (-t309 * t118 + t120 * t164) * qJD(4)) * pkin(3), -g(1) * t204 - g(2) * t179 - g(3) * t193 - t110 * t239 + t241 * t9 + t8 * t306 + (t208 - t49) * t45 - t205 * t44, t14, t4, t17, t15, t16, -t280, -t24 * t111 + t152 * t31 + t205 * t91 + (-t111 * t208 + t197) * t163 + (-t151 * t252 + t181) * t167 + t218, -t315 * t111 - t152 * t30 + t197 * t167 + t205 * t93 + t183, t24 * t93 + t25 * t91 + (-t91 * t208 - t151 * t31 + (t151 * t93 + t195) * qJD(5)) * t167 + (t93 * t208 - t151 * t30 - t3 + (t151 * t91 - t23) * qJD(5)) * t163 - t184 + t316, t208 * t296 + t7 * t152 + t195 * t24 - g(1) * (t78 * pkin(9) + t204 + t71) - g(2) * (t76 * pkin(9) + t179 + t70) - g(3) * (t106 * pkin(9) + t102 + t193) + t205 * t38 + t315 * t23 + (t195 * t250 + t1 - t305) * t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, t61, t41, -t278, t42, t242, t45 * qJD(3) + t176, t243 * t44 + t175, 0, 0, t14, t4, t17, t15, t16, -t280, -pkin(4) * t31 - t26 * t111 - t45 * t91 + (-pkin(9) * t55 + t297) * t163 + (-pkin(9) * t252 + t181) * t167 + t218, t38 * t279 + pkin(4) * t30 + t27 * t111 - t45 * t93 + (t111 * t251 - t290) * pkin(9) + t183, t26 * t93 + t27 * t91 + (-t287 - t286 + (t288 + t292) * qJD(5)) * pkin(9) + t174 + t316, -t7 * pkin(4) - g(1) * t71 - g(2) * t70 - g(3) * t102 + t195 * t26 - t23 * t27 - t38 * t45 + (t174 + t1) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t304, -t91 ^ 2 + t93 ^ 2, -t30 + t285, -t304, t284 - t31, t55, t23 * t111 - t38 * t93 - g(1) * (t114 * t167 - t78 * t163) - g(2) * (t112 * t167 - t76 * t163) - g(3) * (-t106 * t163 - t233) + t3, -t195 * t111 + t38 * t91 - g(1) * (-t114 * t163 - t78 * t167) - g(2) * (-t112 * t163 - t76 * t167) - g(3) * (-t106 * t167 + t234) - t2, 0, 0;];
tau_reg = t5;
