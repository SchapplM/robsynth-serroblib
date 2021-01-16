% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:43
% EndTime: 2021-01-15 22:03:08
% DurationCPUTime: 5.19s
% Computational Cost: add. (5414->509), mult. (15415->722), div. (0->0), fcn. (12610->18), ass. (0->277)
t194 = sin(pkin(5));
t193 = sin(pkin(10));
t198 = sin(qJ(2));
t202 = cos(qJ(2));
t327 = cos(pkin(10));
t216 = t193 * t202 + t327 * t198;
t132 = t216 * t194;
t124 = qJD(1) * t132;
t201 = cos(qJ(4));
t195 = cos(pkin(5));
t300 = qJD(1) * t195;
t246 = qJD(2) + t300;
t160 = t201 * t246;
t197 = sin(qJ(4));
t93 = t124 * t197 - t160;
t92 = qJD(5) + t93;
t253 = t327 * t202;
t239 = qJD(1) * t253;
t299 = qJD(1) * t198;
t267 = t194 * t299;
t121 = t193 * t267 - t194 * t239;
t118 = qJD(4) + t121;
t199 = sin(qJ(1));
t309 = t199 * t202;
t203 = cos(qJ(1));
t313 = t198 * t203;
t146 = -t195 * t309 - t313;
t304 = t202 * t203;
t314 = t198 * t199;
t220 = t195 * t304 - t314;
t360 = -g(1) * t146 - g(2) * t220;
t190 = qJ(2) + pkin(10);
t187 = sin(t190);
t303 = t203 * t187;
t188 = cos(t190);
t311 = t199 * t188;
t130 = t195 * t303 + t311;
t209 = -t201 * t124 - t197 * t246;
t236 = g(1) * t199 - g(2) * t203;
t320 = t194 * t201;
t179 = t195 * t201;
t322 = t194 * t197;
t345 = g(3) * (-t187 * t322 + t179);
t170 = t203 * t188;
t312 = t199 * t187;
t356 = t195 * t312 - t170;
t317 = t195 * t198;
t173 = pkin(1) * t317;
t339 = pkin(7) + qJ(3);
t263 = t194 * t339;
t110 = (t202 * t263 + t173) * qJD(1);
t254 = t327 * t110;
t347 = pkin(1) * t195;
t174 = t202 * t347;
t168 = qJD(1) * t174;
t262 = t198 * t339;
t242 = t194 * t262;
t214 = t195 * pkin(2) - t242;
t96 = qJD(2) * pkin(2) + t214 * qJD(1) + t168;
t50 = t193 * t96 + t254;
t43 = t246 * pkin(8) + t50;
t184 = pkin(2) * t202 + pkin(1);
t152 = t184 * t194;
t143 = -qJD(1) * t152 + qJD(3);
t60 = pkin(3) * t121 - pkin(8) * t124 + t143;
t23 = t197 * t60 + t201 * t43;
t288 = qJDD(1) * t195;
t171 = qJDD(2) + t288;
t286 = t202 * qJDD(1);
t281 = pkin(1) * t286;
t167 = t195 * t281;
t282 = qJD(2) * t347;
t243 = qJD(1) * t282;
t255 = qJD(2) * t339;
t297 = qJD(3) * t198;
t47 = -t198 * t243 + pkin(2) * t171 + t167 + (-qJDD(1) * t262 + (-t202 * t255 - t297) * qJD(1)) * t194;
t213 = qJD(3) * t202 - t198 * t255;
t259 = t194 * t286;
t269 = pkin(7) * t259 + qJDD(1) * t173 + t202 * t243;
t57 = (qJ(3) * t286 + t213 * qJD(1)) * t194 + t269;
t21 = t193 * t47 + t327 * t57;
t18 = pkin(8) * t171 + t21;
t234 = t184 * qJDD(1);
t298 = qJD(2) * t198;
t266 = t194 * t298;
t241 = qJD(1) * t266;
t285 = pkin(2) * t241 + qJDD(3);
t108 = -t194 * t234 + t285;
t211 = qJD(2) * t216;
t240 = t194 * t253;
t287 = t198 * qJDD(1);
t88 = -qJDD(1) * t240 + (qJD(1) * t211 + t193 * t287) * t194;
t89 = -t193 * t241 + (qJD(2) * t239 + t216 * qJDD(1)) * t194;
t33 = pkin(3) * t88 - pkin(8) * t89 + t108;
t257 = t18 * t197 - t201 * t33;
t87 = qJDD(4) + t88;
t4 = -pkin(4) * t87 + t23 * qJD(4) + t257;
t359 = (g(1) * t356 - g(2) * t130) * t197 + t236 * t320 + (-pkin(4) * t209 + t92 * pkin(9)) * t92 + t4 + t345;
t186 = pkin(5) - t190;
t358 = sin(t186) / 0.2e1;
t196 = sin(qJ(5));
t200 = cos(qJ(5));
t54 = -t200 * t118 - t196 * t209;
t357 = t118 * t54;
t355 = -t132 * t197 + t179;
t319 = t194 * t202;
t354 = (-t339 * t319 - t173) * qJD(2) - t194 * t297;
t302 = pkin(7) * t319 + t173;
t189 = t194 ^ 2;
t348 = pkin(1) * t189;
t353 = t195 * t302 + t198 * t348;
t294 = qJD(4) * t197;
t37 = qJD(4) * t160 - t124 * t294 + t197 * t171 + t201 * t89;
t56 = t118 * t196 - t200 * t209;
t12 = t56 * qJD(5) + t196 * t37 - t200 * t87;
t38 = -t209 * qJD(4) - t201 * t171 + t197 * t89;
t15 = pkin(9) * t118 + t23;
t101 = t193 * t110;
t49 = t327 * t96 - t101;
t42 = -t246 * pkin(3) - t49;
t19 = t93 * pkin(4) + pkin(9) * t209 + t42;
t233 = t15 * t196 - t19 * t200;
t292 = qJD(4) * t201;
t219 = -t201 * t18 - t197 * t33 - t60 * t292 + t43 * t294;
t3 = pkin(9) * t87 - t219;
t20 = -t193 * t57 + t327 * t47;
t17 = -pkin(3) * t171 - t20;
t6 = pkin(4) * t38 - pkin(9) * t37 + t17;
t1 = -t233 * qJD(5) + t196 * t6 + t200 * t3;
t350 = 0.2e1 * t189;
t204 = qJD(1) ^ 2;
t185 = pkin(5) + t190;
t175 = sin(t185);
t349 = -t175 / 0.2e1 + t358;
t344 = g(3) * t202;
t343 = t54 * t92;
t342 = t56 * t92;
t80 = -t121 * t196 * t201 - t200 * t124;
t341 = t80 * t92;
t307 = t200 * t201;
t81 = -t121 * t307 + t124 * t196;
t340 = t81 * t92;
t107 = t174 + t214;
t119 = qJ(3) * t319 + t302;
t73 = t193 * t107 + t327 * t119;
t62 = pkin(8) * t195 + t73;
t321 = t194 * t198;
t131 = t193 * t321 - t240;
t79 = pkin(3) * t131 - pkin(8) * t132 - t152;
t229 = t197 * t79 + t201 * t62;
t109 = -qJD(1) * t242 + t168;
t64 = t327 * t109 - t101;
t75 = pkin(2) * t267 + pkin(3) * t124 + pkin(8) * t121;
t338 = t197 * t75 + t201 * t64;
t290 = qJD(5) * t200;
t291 = qJD(5) * t196;
t11 = t118 * t290 + t196 * t87 + t200 * t37 + t209 * t291;
t337 = t11 * t196;
t336 = t118 * t93;
t335 = t121 * t56;
t334 = t124 * t93;
t333 = t124 * t209;
t36 = qJDD(5) + t38;
t332 = t196 * t36;
t331 = t197 * t87;
t330 = t200 * t36;
t251 = t200 * t92;
t329 = t209 * t118;
t63 = t109 * t193 + t254;
t328 = -t63 + t118 * (pkin(4) * t197 - pkin(9) * t201);
t326 = t118 * t197;
t324 = t188 * t194;
t323 = t189 * t204;
t318 = t195 * t197;
t316 = t196 * t199;
t315 = t196 * t203;
t310 = t199 * t201;
t308 = t200 * t199;
t306 = t200 * t203;
t305 = t201 * t203;
t191 = t198 ^ 2;
t301 = -t202 ^ 2 + t191;
t180 = pkin(2) * t193 + pkin(8);
t296 = qJD(4) * t180;
t295 = qJD(4) * t196;
t293 = qJD(4) * t200;
t289 = qJD(1) * qJD(2);
t284 = g(3) * t187 * t194;
t280 = t92 * t295;
t279 = t92 * t293;
t277 = t196 * t324;
t276 = t200 * t324;
t275 = t202 * t323;
t274 = t199 * t322;
t273 = t196 * t310;
t272 = t196 * t305;
t271 = t199 * t307;
t270 = t200 * t305;
t164 = t203 * t322;
t264 = t194 * t195 * t204;
t261 = t202 * t289;
t260 = t194 * t287;
t169 = t202 * t282;
t97 = t213 * t194 + t169;
t44 = t193 * t97 - t327 * t354;
t181 = -t327 * pkin(2) - pkin(3);
t148 = -t201 * pkin(4) - t197 * pkin(9) + t181;
t250 = pkin(9) * t124 - qJD(5) * t148 + t338;
t249 = t130 * t201 - t164;
t248 = t118 * t201;
t247 = pkin(2) * t266;
t245 = qJD(2) + 0.2e1 * t300;
t244 = t171 + t288;
t237 = g(1) * t203 + g(2) * t199;
t8 = t15 * t200 + t19 * t196;
t29 = pkin(9) * t131 + t229;
t105 = t132 * t201 + t318;
t72 = t327 * t107 - t193 * t119;
t61 = -t195 * pkin(3) - t72;
t32 = -pkin(4) * t355 - t105 * pkin(9) + t61;
t232 = t196 * t32 + t200 * t29;
t231 = -t196 * t29 + t200 * t32;
t22 = -t197 * t43 + t201 * t60;
t45 = t354 * t193 + t327 * t97;
t123 = t194 * t211;
t126 = (-t193 * t198 + t253) * t194 * qJD(2);
t76 = pkin(3) * t123 - pkin(8) * t126 + t247;
t230 = -t197 * t45 + t201 * t76;
t228 = -t197 * t62 + t201 * t79;
t78 = t105 * t200 + t131 * t196;
t77 = t105 * t196 - t200 * t131;
t227 = -t118 * t294 - t121 * t326 + t201 * t87;
t177 = cos(t185);
t178 = cos(t186);
t151 = t177 + t178;
t226 = t312 - t203 * t151 / 0.2e1;
t225 = t303 + t199 * t151 / 0.2e1;
t224 = -t92 * t290 - t332;
t223 = t92 * t291 - t330;
t222 = -t201 * t356 + t274;
t218 = t197 * t76 + t201 * t45 + t79 * t292 - t62 * t294;
t217 = t118 * t42 - t180 * t87;
t14 = -pkin(4) * t118 - t22;
t210 = -pkin(9) * t36 + (t14 + t22) * t92;
t208 = g(1) * (t195 * t311 + t303) - g(2) * (t195 * t170 - t312) - g(3) * t324 - t17;
t134 = t195 * t271 - t315;
t142 = t195 * t316 + t270;
t207 = -t134 * t187 + t142 * t188 + t200 * t274;
t135 = t195 * t272 - t308;
t140 = -t195 * t306 - t273;
t206 = -t135 * t187 + t140 * t188 + t196 * t164;
t2 = -t8 * qJD(5) - t196 * t3 + t200 * t6;
t149 = pkin(2) * t317 - t263;
t147 = -t195 * t314 + t304;
t145 = -t195 * t313 - t309;
t141 = t195 * t308 - t272;
t139 = t195 * t315 - t271;
t138 = t187 * t320 + t318;
t136 = t195 * t270 + t316;
t133 = t195 * t273 + t306;
t113 = t199 * t349 + t170;
t112 = t203 * t349 - t311;
t99 = t194 * t310 + t197 * t356;
t98 = t130 * t197 + t194 * t305;
t71 = t105 * qJD(4) + t126 * t197;
t70 = t355 * qJD(4) + t126 * t201;
t66 = t133 * t187 + t141 * t188 - t196 * t274;
t65 = -t136 * t187 + t139 * t188 + t200 * t164;
t48 = t56 * t294;
t28 = -pkin(4) * t131 - t228;
t27 = t78 * qJD(5) - t200 * t123 + t196 * t70;
t26 = -t77 * qJD(5) + t123 * t196 + t200 * t70;
t24 = -pkin(4) * t124 + t197 * t64 - t201 * t75;
t13 = pkin(4) * t71 - pkin(9) * t70 + t44;
t10 = -pkin(4) * t123 + t229 * qJD(4) - t230;
t9 = pkin(9) * t123 + t218;
t5 = [qJDD(1), t236, t237, (qJDD(1) * t191 + 0.2e1 * t198 * t261) * t189, (t198 * t286 - t289 * t301) * t350, (qJD(2) * t202 * t245 + t198 * t244) * t194, (t202 * t244 - t245 * t298) * t194, t171 * t195, t281 * t350 + (-pkin(7) * t321 + t174) * t171 + (-pkin(7) * t260 + t167) * t195 - g(1) * t145 - g(2) * t147 - t302 * qJD(2) ^ 2 - 0.2e1 * t353 * t289, -(-pkin(7) * t266 + t169) * t246 - t302 * t171 - (-pkin(7) * t241 + t269) * t195 + g(1) * t220 - g(2) * t146 + 0.2e1 * (-t261 - t287) * t348, -g(1) * t112 - g(2) * t113 + t108 * t131 + t121 * t247 + t143 * t123 - t152 * t88 + t72 * t171 + t20 * t195 - t246 * t44, -g(1) * t226 + g(2) * t225 + t108 * t132 + t124 * t247 + t143 * t126 - t152 * t89 - t73 * t171 - t21 * t195 - t246 * t45, -t121 * t45 - t123 * t50 + t124 * t44 - t126 * t49 - t131 * t21 - t132 * t20 - t194 * t237 - t72 * t89 - t73 * t88, t21 * t73 + t50 * t45 + t20 * t72 - t49 * t44 - t108 * t152 + t143 * t247 - g(1) * (-t149 * t203 - t184 * t199) - g(2) * (-t149 * t199 + t184 * t203), t105 * t37 - t209 * t70, -t105 * t38 + t209 * t71 + t355 * t37 - t70 * t93, t105 * t87 + t118 * t70 - t123 * t209 + t131 * t37, -t118 * t71 - t123 * t93 - t131 * t38 + t355 * t87, t118 * t123 + t131 * t87, t230 * t118 + t228 * t87 - t257 * t131 + t22 * t123 + t44 * t93 + t61 * t38 - t17 * t355 + t42 * t71 + g(1) * t249 - g(2) * t222 + (-t118 * t229 - t131 * t23) * qJD(4), -g(1) * t98 - g(2) * t99 + t17 * t105 - t118 * t218 - t23 * t123 + t131 * t219 - t209 * t44 - t229 * t87 + t61 * t37 + t42 * t70, t11 * t78 + t26 * t56, -t11 * t77 - t12 * t78 - t26 * t54 - t27 * t56, -t11 * t355 + t26 * t92 + t36 * t78 + t56 * t71, t12 * t355 - t27 * t92 - t36 * t77 - t54 * t71, -t355 * t36 + t71 * t92, (-qJD(5) * t232 + t13 * t200 - t196 * t9) * t92 + t231 * t36 - t2 * t355 - t233 * t71 + t10 * t54 + t28 * t12 + t4 * t77 + t14 * t27 - g(1) * t65 - g(2) * t207, -(qJD(5) * t231 + t13 * t196 + t200 * t9) * t92 - t232 * t36 + t1 * t355 - t8 * t71 + t10 * t56 + t28 * t11 + t4 * t78 + t14 * t26 + g(1) * t206 - g(2) * t66; 0, 0, 0, -t198 * t275, t301 * t323, -t202 * t264 + t260, t198 * t264 + t259, t171, t167 + (-pkin(7) * t287 - t344) * t194 + t353 * t204 + t360, pkin(1) * t275 + (-pkin(7) * t267 + t168) * t300 + g(1) * t147 - g(2) * t145 + g(3) * t321 + t168 * qJD(2) - t269, t63 * t246 - t143 * t124 + g(1) * t225 + g(2) * t226 - g(3) * (t358 + t175 / 0.2e1) + (-t121 * t267 + t171 * t327) * pkin(2) + t20, t64 * t246 + t143 * t121 + g(1) * t113 - g(2) * t112 - g(3) * (-t178 / 0.2e1 + t177 / 0.2e1) + (-t124 * t267 - t171 * t193) * pkin(2) - t21, (t50 - t63) * t124 + (-t49 + t64) * t121 + (-t193 * t88 - t327 * t89) * pkin(2), t49 * t63 - t50 * t64 + (t21 * t193 + t20 * t327 + (-t143 * t299 - t344) * t194 + t360) * pkin(2), t197 * t37 - t209 * t248, (t37 - t336) * t201 + (-t38 + t329) * t197, t118 * t248 + t331 + t333, t227 + t334, -t118 * t124, -t22 * t124 + t181 * t38 - t63 * t93 + (t64 * t118 + t217) * t197 + ((-t75 - t296) * t118 + t208) * t201, t181 * t37 + t338 * t118 + t23 * t124 + t63 * t209 + t217 * t201 + (t118 * t296 - t208) * t197, t11 * t197 * t200 + (-t197 * t291 + t200 * t292 - t81) * t56, t54 * t81 + t56 * t80 + (-t196 * t56 - t200 * t54) * t292 + (-t337 - t12 * t200 + (t196 * t54 - t200 * t56) * qJD(5)) * t197, -t340 + t48 + (-t11 + t279) * t201 + (-t223 + t335) * t197, t341 + (t12 - t280) * t201 + (t224 - t357) * t197, -t201 * t36 + t326 * t92, t148 * t330 - t24 * t54 - t14 * t80 - g(1) * (-t134 * t188 - t142 * t187) - g(2) * (t136 * t188 + t139 * t187) - t196 * t284 + (t196 * t250 + t200 * t328) * t92 + (-g(3) * t276 + t14 * t295 - t2 + (qJD(4) * t54 + t224) * t180) * t201 + (t14 * t290 + t180 * t12 - t233 * t121 + t4 * t196 + (t180 * t196 * t92 - t233) * qJD(4)) * t197, -t148 * t332 - t24 * t56 - t14 * t81 - g(1) * (t133 * t188 - t141 * t187) - g(2) * (-t135 * t188 - t140 * t187) - t200 * t284 + (-t196 * t328 + t200 * t250) * t92 + (g(3) * t277 + t14 * t293 + t1 + (qJD(4) * t56 + t223) * t180) * t201 + (-t14 * t291 + t180 * t11 - t8 * t121 + t4 * t200 + (t180 * t251 - t8) * qJD(4)) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124 * t246 + t88, -t121 * t246 + t89, -t121 ^ 2 - t124 ^ 2, -g(3) * t195 + t121 * t50 + t124 * t49 + (-t234 - t236) * t194 + t285, 0, 0, 0, 0, 0, t227 - t334, -t118 ^ 2 * t201 - t331 + t333, 0, 0, 0, 0, 0, t341 + (-t12 - t280) * t201 + (t224 + t357) * t197, t340 + t48 + (-t11 - t279) * t201 + (t223 + t335) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209 * t93, t209 ^ 2 - t93 ^ 2, t37 + t336, -t38 - t329, t87, -g(1) * t99 + g(2) * t98 + t42 * t209 - t257 - t345 + (-qJD(4) + t118) * t23, g(1) * t222 + g(2) * t249 + g(3) * t138 + t22 * t118 + t42 * t93 + t219, t251 * t56 + t337, (t11 - t343) * t200 + (-t12 - t342) * t196, t209 * t56 + t251 * t92 + t332, -t196 * t92 ^ 2 - t209 * t54 + t330, t92 * t209, -pkin(4) * t12 + t210 * t196 - t359 * t200 - t233 * t209 - t23 * t54, -pkin(4) * t11 + t359 * t196 + t210 * t200 - t8 * t209 - t23 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, t11 + t343, -t12 + t342, t36, t8 * t92 - t14 * t56 - g(1) * t66 - g(2) * t206 - g(3) * (-t138 * t196 - t276) + t2, -t233 * t92 + t14 * t54 + g(1) * t207 - g(2) * t65 - g(3) * (-t138 * t200 + t277) - t1;];
tau_reg = t5;
