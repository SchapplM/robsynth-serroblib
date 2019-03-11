% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:31
% EndTime: 2019-03-09 04:48:42
% DurationCPUTime: 4.33s
% Computational Cost: add. (5613->490), mult. (10948->620), div. (0->0), fcn. (6973->10), ass. (0->244)
t192 = sin(qJ(3));
t278 = qJD(1) * t192;
t160 = qJD(4) + t278;
t191 = sin(qJ(4));
t194 = cos(qJ(4));
t266 = t194 * qJD(3);
t195 = cos(qJ(3));
t277 = qJD(1) * t195;
t128 = t191 * t277 - t266;
t274 = qJD(3) * t191;
t130 = t194 * t277 + t274;
t188 = sin(pkin(9));
t189 = cos(pkin(9));
t74 = t189 * t128 + t130 * t188;
t339 = t160 * t74;
t126 = t188 * t194 + t189 * t191;
t106 = t126 * qJD(4);
t107 = t126 * qJD(1);
t315 = t192 * t107 + t106;
t303 = t189 * t194;
t304 = t188 * t191;
t125 = -t303 + t304;
t269 = qJD(4) * t192;
t272 = qJD(3) * t195;
t338 = t125 * t272 + t126 * t269 + t107;
t241 = qJD(1) * t303;
t250 = t191 * t278;
t268 = qJD(4) * t194;
t270 = qJD(4) * t191;
t334 = t188 * t270 - t189 * t268;
t337 = t188 * t250 - t192 * t241 + t334;
t193 = sin(qJ(1));
t296 = t193 * t194;
t196 = cos(qJ(1));
t297 = t192 * t196;
t112 = t191 * t297 + t296;
t221 = -t128 * t188 + t189 * t130;
t336 = t221 ^ 2;
t324 = qJ(5) + pkin(8);
t234 = qJD(4) * t324;
t265 = t194 * qJD(5);
t100 = -t191 * t234 + t265;
t207 = -t191 * qJD(5) - t194 * t234;
t227 = pkin(3) * t195 + pkin(8) * t192;
t132 = t227 * qJD(1);
t115 = t194 * t132;
t197 = -pkin(1) - pkin(7);
t155 = t197 * qJD(1) + qJD(2);
t298 = t192 * t194;
t301 = t191 * t195;
t54 = -t155 * t301 + t115 + (pkin(4) * t195 + qJ(5) * t298) * qJD(1);
t293 = t194 * t195;
t285 = t191 * t132 + t155 * t293;
t63 = qJ(5) * t250 + t285;
t322 = (t207 - t54) * t189 + (-t100 + t63) * t188;
t137 = pkin(3) * t192 - pkin(8) * t195 + qJ(2);
t283 = t191 * t137 + t197 * t298;
t184 = g(2) * t196;
t328 = g(1) * t193;
t335 = -t184 + t328;
t327 = g(3) * t192;
t205 = -t195 * t335 + t327;
t273 = qJD(3) * t192;
t233 = -qJDD(3) * pkin(3) + t155 * t273;
t151 = t197 * qJDD(1) + qJDD(2);
t307 = t151 * t195;
t83 = t233 - t307;
t333 = -qJD(4) * pkin(8) * t160 + t205 - t83;
t247 = t192 * t266;
t267 = qJD(4) * t195;
t208 = -t191 * t267 - t247;
t260 = t195 * qJDD(1);
t67 = qJD(1) * t208 + qJD(4) * t266 + t191 * qJDD(3) + t194 * t260;
t249 = t191 * t273;
t68 = -qJD(1) * t249 + qJD(4) * t130 - t194 * qJDD(3) + t191 * t260;
t33 = t188 * t67 + t189 * t68;
t34 = -t188 * t68 + t189 * t67;
t332 = t33 * pkin(5) - t34 * qJ(6) - t221 * qJD(6);
t185 = qJ(4) + pkin(9);
t174 = sin(t185);
t306 = t155 * t195;
t118 = -qJD(3) * pkin(3) - t306;
t81 = pkin(4) * t128 + qJD(5) + t118;
t24 = pkin(5) * t74 - qJ(6) * t221 + t81;
t136 = t192 * t155;
t117 = qJD(3) * pkin(8) + t136;
t109 = t137 * qJD(1);
t124 = qJD(3) * t227 + qJD(2);
t79 = qJD(1) * t124 + qJDD(1) * t137;
t84 = qJDD(3) * pkin(8) + t151 * t192 + t155 * t272;
t258 = -t109 * t268 - t191 * t79 - t194 * t84;
t212 = -t117 * t270 - t258;
t15 = -qJ(5) * t68 - qJD(5) * t128 + t212;
t264 = qJD(1) * qJD(3);
t239 = t195 * t264;
t261 = t192 * qJDD(1);
t123 = qJDD(4) + t239 + t261;
t66 = t109 * t191 + t117 * t194;
t70 = t194 * t79;
t8 = t123 * pkin(4) - t67 * qJ(5) - qJD(4) * t66 - t130 * qJD(5) - t191 * t84 + t70;
t3 = -t188 * t15 + t189 * t8;
t255 = -qJDD(6) + t3;
t326 = g(3) * t195;
t175 = cos(t185);
t289 = t196 * t175;
t299 = t192 * t193;
t91 = t174 * t299 - t289;
t93 = t174 * t297 + t175 * t193;
t331 = g(1) * t91 - g(2) * t93 + t174 * t326 - t24 * t221 + t255;
t4 = t189 * t15 + t188 * t8;
t329 = pkin(5) * t123;
t52 = -qJ(5) * t128 + t66;
t48 = t189 * t52;
t65 = t194 * t109 - t117 * t191;
t51 = -qJ(5) * t130 + t65;
t21 = t188 * t51 + t48;
t325 = t21 * t221;
t102 = t194 * t124;
t238 = -t191 * t197 + pkin(4);
t32 = qJ(5) * t247 + t102 - t283 * qJD(4) + (qJ(5) * t270 + qJD(3) * t238 - t265) * t195;
t243 = t194 * t267;
t271 = qJD(3) * t197;
t245 = t195 * t271;
t252 = t191 * t124 + t137 * t268 + t194 * t245;
t37 = -qJ(5) * t243 + (-qJD(5) * t195 + (qJ(5) * qJD(3) - qJD(4) * t197) * t192) * t191 + t252;
t12 = t188 * t32 + t189 * t37;
t323 = -qJD(6) * t126 - t136 + t337 * qJ(6) + t315 * pkin(5) + (t250 + t270) * pkin(4);
t46 = pkin(4) * t160 + t51;
t20 = t188 * t46 + t48;
t29 = t188 * t54 + t189 * t63;
t321 = pkin(5) * t277 - t322;
t26 = qJ(6) * t277 + t29;
t61 = t189 * t100 + t188 * t207;
t320 = t61 - t26;
t122 = t194 * t137;
t72 = -qJ(5) * t293 + t192 * t238 + t122;
t82 = -qJ(5) * t301 + t283;
t42 = t188 * t72 + t189 * t82;
t319 = t188 * t52;
t317 = t67 * t191;
t246 = t195 * t266;
t248 = t191 * t272;
t316 = qJD(1) * t304 - t188 * t246 - t189 * t248 + t334 * t192 - t241;
t312 = pkin(1) * qJDD(1);
t199 = qJD(1) ^ 2;
t311 = qJ(2) * t199;
t310 = t128 * t160;
t309 = t130 * t160;
t308 = t130 * t194;
t305 = t160 * t191;
t302 = t191 * t123;
t300 = t191 * t196;
t295 = t193 * t195;
t294 = t194 * t123;
t292 = t194 * t196;
t291 = t195 * t324;
t290 = t195 * t196;
t198 = qJD(3) ^ 2;
t288 = t197 * t198;
t22 = t189 * t51 - t319;
t286 = qJD(6) - t22;
t284 = t112 * pkin(4);
t282 = g(1) * t290 + g(2) * t295;
t281 = t196 * pkin(1) + t193 * qJ(2);
t187 = t195 ^ 2;
t280 = t192 ^ 2 - t187;
t279 = -t198 - t199;
t276 = qJD(3) * t128;
t275 = qJD(3) * t130;
t263 = qJDD(1) * qJ(2);
t262 = qJDD(3) * t192;
t259 = t123 * qJ(6) + t4;
t256 = 0.2e1 * qJD(1) * qJD(2);
t254 = t191 * t299;
t173 = pkin(4) * t194 + pkin(3);
t181 = t196 * qJ(2);
t251 = t173 * t297 - t290 * t324 + t181;
t240 = t324 * t191;
t236 = -t151 - t184;
t235 = -g(1) * t295 + t327;
t232 = -qJD(4) * t109 - t84;
t231 = t160 * t197 + t117;
t161 = pkin(4) * t301;
t230 = -t195 * t197 + t161;
t229 = qJDD(2) - t312;
t228 = qJD(1) + t269;
t226 = g(1) * t196 + g(2) * t193;
t224 = -t74 ^ 2 - t336;
t223 = pkin(5) * t175 + qJ(6) * t174;
t222 = t311 + t328;
t11 = -t188 * t37 + t189 * t32;
t19 = t189 * t46 - t319;
t41 = -t188 * t82 + t189 * t72;
t220 = g(1) * t299 + t326;
t219 = t256 + 0.2e1 * t263;
t218 = t173 + t223;
t216 = (-pkin(4) * t191 + t197) * t328;
t215 = pkin(4) * t300 + t196 * pkin(7) + t173 * t299 - t193 * t291 + t281;
t214 = t160 * t268 + t302;
t213 = -t160 * t270 + t294;
t211 = t192 * t271 + (t243 - t249) * pkin(4);
t95 = t126 * t192;
t210 = 0.2e1 * qJ(2) * t264 + qJDD(3) * t197;
t209 = t68 * pkin(4) + qJDD(5) + t233;
t206 = -pkin(8) * t123 + t118 * t160;
t43 = t209 - t307;
t97 = t125 * t192;
t204 = -t221 * t316 + t97 * t33 + t338 * t74 + t95 * t34;
t203 = t219 - t226;
t166 = g(2) * t297;
t144 = t324 * t194;
t85 = t188 * t144 + t189 * t240;
t86 = t189 * t144 - t188 * t240;
t201 = -t86 * t33 + t85 * t34 - t61 * t74 + t166 - t220;
t200 = t195 * t236 + t209 - t235;
t178 = qJDD(3) * t195;
t171 = -pkin(4) * t189 - pkin(5);
t169 = pkin(4) * t188 + qJ(6);
t164 = pkin(4) * t292;
t149 = t324 * t297;
t134 = t173 * t295;
t113 = -t191 * t193 + t192 * t292;
t111 = t192 * t296 + t300;
t110 = -t254 + t292;
t98 = -t188 * t301 + t189 * t293;
t96 = t126 * t195;
t94 = -t174 * t193 + t192 * t289;
t92 = t174 * t196 + t175 * t299;
t71 = pkin(5) * t125 - qJ(6) * t126 - t173;
t59 = t106 * t195 - t188 * t249 + t189 * t247;
t57 = qJD(3) * t95 + t125 * t267;
t50 = pkin(5) * t96 - qJ(6) * t98 + t230;
t38 = -pkin(5) * t192 - t41;
t36 = qJ(6) * t192 + t42;
t30 = pkin(4) * t130 + pkin(5) * t221 + qJ(6) * t74;
t18 = qJ(6) * t160 + t20;
t17 = -pkin(5) * t160 + qJD(6) - t19;
t16 = -pkin(5) * t57 + qJ(6) * t59 - qJD(6) * t98 + t211;
t10 = -pkin(5) * t272 - t11;
t9 = qJ(6) * t272 + qJD(6) * t192 + t12;
t5 = t43 + t332;
t2 = -t255 - t329;
t1 = qJD(6) * t160 + t259;
t6 = [qJDD(1), t335, t226, qJDD(2) - t335 - 0.2e1 * t312, t203, -t229 * pkin(1) - g(1) * (-pkin(1) * t193 + t181) - g(2) * t281 + (t256 + t263) * qJ(2), qJDD(1) * t187 - 0.2e1 * t192 * t239, -0.2e1 * t192 * t260 + 0.2e1 * t264 * t280, -t192 * t198 + t178, -t195 * t198 - t262, 0, t210 * t195 + (t203 - t288) * t192, -t210 * t192 + (t219 - t288) * t195 - t282, t130 * t208 + t293 * t67 (t128 * t194 + t130 * t191) * t273 + (-t317 - t194 * t68 + (t128 * t191 - t308) * qJD(4)) * t195 (-t160 * t266 + t67) * t192 + (t213 + t275) * t195 (t160 * t274 - t68) * t192 + (-t214 - t276) * t195, t123 * t192 + t160 * t272, -g(1) * t113 - g(2) * t111 + t102 * t160 + t122 * t123 + (t128 * t271 - t231 * t268 + t70) * t192 + (qJD(3) * t65 + t118 * t268 - t197 * t68) * t195 + ((-qJD(4) * t137 - t245) * t160 + t83 * t195 + (-qJD(3) * t118 - t123 * t197 + t232) * t192) * t191, -t252 * t160 - t283 * t123 + g(1) * t112 - g(2) * t110 + (t231 * t270 + (-t118 * t194 + t130 * t197) * qJD(3) + t258) * t192 + (-qJD(3) * t66 - t118 * t270 + t83 * t194 - t197 * t67) * t195, -t11 * t221 - t12 * t74 + t19 * t59 + t20 * t57 - t3 * t98 - t33 * t42 - t34 * t41 - t4 * t96 + t282, -g(1) * t251 - g(2) * t215 + t19 * t11 + t20 * t12 + t211 * t81 + t230 * t43 + t3 * t41 + t4 * t42 - t216, -g(1) * t94 - g(2) * t92 - t10 * t160 - t123 * t38 + t16 * t74 - t17 * t272 - t192 * t2 - t24 * t57 + t33 * t50 + t5 * t96, -t1 * t96 + t10 * t221 - t17 * t59 + t18 * t57 + t2 * t98 - t33 * t36 + t34 * t38 - t74 * t9 + t282, -g(1) * t93 - g(2) * t91 + t1 * t192 + t123 * t36 - t16 * t221 + t160 * t9 + t18 * t272 + t24 * t59 - t34 * t50 - t5 * t98, t1 * t36 + t18 * t9 + t5 * t50 + t24 * t16 + t2 * t38 + t17 * t10 - g(1) * (t94 * pkin(5) + t93 * qJ(6) + t251) - g(2) * (pkin(5) * t92 + qJ(6) * t91 + t215) - t216; 0, 0, 0, qJDD(1), -t199, t184 - t222 + t229, 0, 0, 0, 0, 0, t192 * t279 + t178, t195 * t279 - t262, 0, 0, 0, 0, 0, -t195 * t68 + (t276 - t302) * t192 + (-t194 * t228 - t248) * t160, -t195 * t67 + (t275 - t294) * t192 + (t191 * t228 - t246) * t160, t204, t19 * t316 - t43 * t195 - t20 * t338 + t273 * t81 - t3 * t95 - t4 * t97 - t335, -t95 * t123 + t160 * t316 - t195 * t33 + t273 * t74, t204, -t97 * t123 - t160 * t338 + t195 * t34 - t221 * t273, -t1 * t97 - t17 * t316 - t18 * t338 - t5 * t195 + t2 * t95 + t24 * t273 - t335; 0, 0, 0, 0, 0, 0, t195 * t199 * t192, -t280 * t199, t260, -t261, qJDD(3) (-t236 - t311) * t195 + t235, t326 - t166 + (-t151 + t222) * t192, t160 * t308 + t317 (t67 - t310) * t194 + (-t68 - t309) * t191 (-t130 * t195 + t160 * t298) * qJD(1) + t214 (t128 * t195 - t192 * t305) * qJD(1) + t213, -t160 * t277, -t65 * t277 - t128 * t136 - pkin(3) * t68 - t115 * t160 + (t160 * t306 + t206) * t191 + t333 * t194, -pkin(3) * t67 - t130 * t136 + t285 * t160 - t191 * t333 + t206 * t194 + t66 * t277, -t4 * t125 - t3 * t126 + t19 * t337 - t20 * t315 - t221 * t322 + t29 * t74 + t201, t4 * t86 - t3 * t85 - t43 * t173 - g(1) * (t299 * t324 + t134) - g(2) * (-t173 * t290 - t149) - g(3) * (-t173 * t192 + t291) + (pkin(4) * t305 - t136) * t81 + (t61 - t29) * t20 + t322 * t19, -t85 * t123 + t5 * t125 - t321 * t160 + t17 * t277 + t205 * t175 + t315 * t24 + t323 * t74 + t71 * t33, -t1 * t125 + t2 * t126 - t17 * t337 - t18 * t315 + t221 * t321 + t26 * t74 + t201, t86 * t123 - t5 * t126 + t320 * t160 + t205 * t174 - t18 * t277 - t221 * t323 + t24 * t337 - t71 * t34, -g(1) * t134 + g(2) * t149 + t1 * t86 + t2 * t85 + t5 * t71 + t323 * t24 + t320 * t18 + t321 * t17 + (g(3) * t218 - t324 * t328) * t192 + (-g(3) * t324 + t218 * t184 - t223 * t328) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t128, -t128 ^ 2 + t130 ^ 2, t67 + t310, t309 - t68, t123, -t117 * t268 - g(1) * t110 - g(2) * t112 - t118 * t130 + t66 * t160 + t70 + (t232 + t326) * t191, g(1) * t111 - g(2) * t113 + g(3) * t293 + t118 * t128 + t160 * t65 - t212, t20 * t221 - t325 + (-t188 * t33 - t189 * t34) * pkin(4) + (-t19 + t22) * t74, -t20 * t22 + t19 * t21 - g(1) * t164 - g(2) * t284 + (-t81 * t130 + t4 * t188 + t3 * t189 + t191 * t220) * pkin(4), t21 * t160 - t30 * t74 + (pkin(5) - t171) * t123 + t331, -t169 * t33 + t171 * t34 + t18 * t221 - t325 + (t17 - t286) * t74, -t175 * t326 - g(1) * t92 + g(2) * t94 + t169 * t123 - t24 * t74 + t30 * t221 + (0.2e1 * qJD(6) - t22) * t160 + t259, t1 * t169 + t2 * t171 - t24 * t30 - t17 * t21 - g(1) * (-pkin(4) * t254 - pkin(5) * t91 + qJ(6) * t92 + t164) - g(2) * (pkin(5) * t93 - qJ(6) * t94 + t284) - g(3) * (-t161 + (-pkin(5) * t174 + qJ(6) * t175) * t195) + t286 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t19 * t221 + t20 * t74 + t200, t160 * t221 + t33, t224, -t34 + t339, -t17 * t221 + t18 * t74 + t200 + t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221 * t74 - t123, t34 + t339, -t160 ^ 2 - t336, -t160 * t18 - t329 - t331;];
tau_reg  = t6;
