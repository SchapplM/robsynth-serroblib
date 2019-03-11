% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:07
% EndTime: 2019-03-09 09:05:21
% DurationCPUTime: 5.37s
% Computational Cost: add. (7058->554), mult. (20175->743), div. (0->0), fcn. (16702->12), ass. (0->285)
t214 = sin(pkin(6));
t223 = cos(qJ(2));
t207 = pkin(2) * t223 + pkin(1);
t264 = t207 * qJDD(1);
t219 = sin(qJ(2));
t320 = qJD(2) * t219;
t295 = t214 * t320;
t272 = qJD(1) * t295;
t307 = pkin(2) * t272 + qJDD(3);
t389 = t214 * t264 - t307;
t213 = sin(pkin(11));
t321 = qJD(1) * t219;
t296 = t214 * t321;
t215 = cos(pkin(11));
t330 = t223 * t215;
t299 = t214 * t330;
t156 = -qJD(1) * t299 + t213 * t296;
t216 = cos(pkin(6));
t322 = qJD(1) * t216;
t197 = qJD(2) + t322;
t218 = sin(qJ(5));
t222 = cos(qJ(5));
t125 = -t222 * t156 + t197 * t218;
t124 = qJD(6) + t125;
t255 = t213 * t223 + t215 * t219;
t311 = qJD(1) * qJD(2);
t292 = t223 * t311;
t111 = -t213 * t272 + t214 * (qJDD(1) * t255 + t215 * t292);
t245 = qJD(1) * t255;
t160 = t214 * t245;
t312 = qJD(5) + t160;
t127 = t156 * t218 + t197 * t222;
t163 = t255 * t214;
t149 = qJD(1) * t163 + qJD(5);
t217 = sin(qJ(6));
t221 = cos(qJ(6));
t80 = t127 * t217 - t221 * t149;
t388 = t312 * t80;
t281 = t222 * t312;
t337 = t216 * t219;
t200 = pkin(1) * t337;
t340 = t214 * t223;
t365 = pkin(8) + qJ(3);
t140 = (t340 * t365 + t200) * qJD(1);
t132 = t213 * t140;
t336 = t216 * t223;
t201 = pkin(1) * t336;
t194 = qJD(1) * t201;
t293 = t365 * t219;
t273 = t214 * t293;
t139 = -qJD(1) * t273 + t194;
t88 = t139 * t215 - t132;
t326 = qJD(4) - t88;
t150 = t160 ^ 2;
t386 = -t156 ^ 2 - t150;
t385 = -qJD(5) + t149;
t109 = qJDD(5) + t111;
t206 = -pkin(2) * t215 - pkin(3);
t203 = -pkin(9) + t206;
t370 = pkin(4) * t156;
t128 = pkin(2) * t197 + t139;
t338 = t215 * t140;
t77 = t213 * t128 + t338;
t66 = -t197 * qJ(4) - t77;
t45 = -t66 - t370;
t384 = -t203 * t109 - t312 * t45;
t308 = qJDD(1) * t223;
t290 = t214 * t308;
t175 = t215 * t290;
t244 = t255 * qJD(2);
t309 = qJDD(1) * t219;
t291 = t213 * t309;
t110 = -t175 + (qJD(1) * t244 + t291) * t214;
t310 = qJDD(1) * t216;
t196 = qJDD(2) + t310;
t316 = qJD(5) * t222;
t317 = qJD(5) * t218;
t49 = t218 * t110 + t156 * t316 + t222 * t196 - t197 * t317;
t82 = t127 * t221 + t149 * t217;
t15 = qJD(6) * t82 - t221 * t109 + t217 * t49;
t250 = -t312 * t149 * t218 + t222 * t109;
t315 = qJD(6) * t217;
t347 = t160 * t218;
t98 = -t156 * t217 + t221 * t347;
t235 = t221 * t317 + t222 * t315 + t98;
t331 = t221 * t222;
t283 = -t222 * t110 + t196 * t218;
t50 = qJD(5) * t127 + t283;
t47 = qJDD(6) + t50;
t383 = -t124 * t235 + t47 * t331;
t220 = sin(qJ(1));
t224 = cos(qJ(1));
t172 = t213 * t219 - t330;
t247 = t172 * t216;
t115 = -t220 * t255 - t224 * t247;
t339 = t214 * t224;
t101 = -t115 * t222 + t218 * t339;
t342 = t214 * t219;
t162 = t213 * t342 - t299;
t257 = t162 * t222 - t216 * t218;
t118 = t220 * t247 - t224 * t255;
t341 = t214 * t220;
t99 = -t118 * t222 - t218 * t341;
t246 = g(1) * t99 + g(2) * t101 + g(3) * t257;
t76 = t128 * t215 - t132;
t265 = qJD(4) - t76;
t369 = pkin(4) * t160;
t375 = pkin(3) + pkin(9);
t40 = -t197 * t375 + t265 + t369;
t270 = t207 * t214;
t165 = -qJD(1) * t270 + qJD(3);
t232 = -qJ(4) * t160 + t165;
t57 = t156 * t375 + t232;
t22 = t218 * t40 + t222 * t57;
t302 = pkin(1) * t308;
t191 = t216 * t302;
t304 = pkin(1) * qJD(2) * t216;
t274 = qJD(1) * t304;
t287 = qJD(2) * t365;
t319 = qJD(3) * t219;
t73 = -t219 * t274 + pkin(2) * t196 + t191 + (-qJDD(1) * t293 + (-t223 * t287 - t319) * qJD(1)) * t214;
t239 = qJD(3) * t223 - t219 * t287;
t297 = pkin(8) * t290 + qJDD(1) * t200 + t223 * t274;
t83 = (qJ(3) * t308 + qJD(1) * t239) * t214 + t297;
t32 = -t213 * t83 + t215 * t73;
t266 = qJDD(4) - t32;
t19 = pkin(4) * t111 - t196 * t375 + t266;
t355 = qJ(4) * t111;
t227 = -qJD(4) * t160 - t355 - t389;
t26 = t110 * t375 + t227;
t288 = -t222 * t19 + t218 * t26;
t4 = -pkin(5) * t109 + qJD(5) * t22 + t288;
t382 = t124 * (pkin(5) * t127 + t124 * pkin(10)) + t246 + t4;
t204 = pkin(2) * t213 + qJ(4);
t170 = pkin(5) * t218 - pkin(10) * t222 + t204;
t269 = pkin(5) * t222 + pkin(10) * t218;
t381 = t124 * (t269 * qJD(5) - (-pkin(4) - t269) * t160 + t326) + t170 * t47;
t325 = -t213 * t336 - t215 * t337;
t114 = t220 * t172 + t224 * t325;
t253 = t115 * t218 + t222 * t339;
t380 = -t114 * t221 + t217 * t253;
t119 = -t224 * t172 + t220 * t325;
t379 = t114 * t217 + t221 * t253;
t12 = pkin(10) * t149 + t22;
t29 = pkin(5) * t125 - pkin(10) * t127 + t45;
t263 = t12 * t217 - t221 * t29;
t252 = t218 * t19 + t222 * t26 + t40 * t316 - t317 * t57;
t3 = pkin(10) * t109 + t252;
t33 = t213 * t73 + t215 * t83;
t278 = t196 * qJ(4) + t197 * qJD(4) + t33;
t20 = -pkin(4) * t110 + t278;
t6 = pkin(5) * t50 - pkin(10) * t49 + t20;
t1 = -t263 * qJD(6) + t217 * t6 + t221 * t3;
t210 = t214 ^ 2;
t376 = 0.2e1 * t210;
t374 = pkin(1) * t210;
t373 = pkin(3) * t110;
t372 = pkin(3) * t162;
t371 = pkin(3) * t196;
t332 = t220 * t223;
t334 = t219 * t224;
t168 = -t216 * t332 - t334;
t368 = g(1) * t168;
t367 = g(1) * t220;
t366 = g(3) * t223;
t87 = t139 * t213 + t338;
t59 = t87 - t370;
t284 = pkin(2) * t296 + qJ(4) * t156;
t64 = t160 * t375 + t284;
t364 = t218 * t59 + t222 * t64;
t138 = pkin(2) * t216 + t201 - t273;
t324 = pkin(8) * t340 + t200;
t151 = qJ(3) * t340 + t324;
t91 = t138 * t215 - t213 * t151;
t58 = pkin(4) * t163 - t216 * t375 - t91;
t354 = qJ(4) * t163;
t240 = -t270 - t354;
t74 = t162 * t375 + t240;
t258 = t218 * t58 + t222 * t74;
t363 = t124 * t80;
t362 = t124 * t82;
t314 = qJD(6) * t221;
t14 = t217 * t109 - t127 * t315 + t149 * t314 + t221 * t49;
t361 = t14 * t217;
t360 = t15 * t218;
t358 = t217 * t47;
t357 = t221 * t47;
t351 = t125 * t156;
t350 = t127 * t156;
t348 = t156 * t197;
t345 = t203 * t217;
t344 = t203 * t221;
t225 = qJD(1) ^ 2;
t343 = t210 * t225;
t335 = t219 * t220;
t329 = t223 * t224;
t327 = t369 + t326;
t195 = t223 * t304;
t129 = t214 * t239 + t195;
t294 = t365 * t214;
t130 = -t214 * t319 + (-t223 * t294 - t200) * qJD(2);
t71 = t215 * t129 + t213 * t130;
t92 = t213 * t138 + t215 * t151;
t211 = t219 ^ 2;
t323 = -t223 ^ 2 + t211;
t318 = qJD(5) * t203;
t313 = qJD(2) - t197;
t306 = t14 * t218 + (t160 * t222 + t316) * t82;
t305 = t219 * t374;
t301 = pkin(8) * t309;
t300 = t223 * t343;
t298 = t216 * t329;
t63 = -t216 * qJD(4) - t71;
t85 = -t216 * qJ(4) - t92;
t97 = t221 * t156 + t217 * t347;
t285 = (t217 * t317 + t97) * t124;
t70 = t129 * t213 - t215 * t130;
t164 = pkin(2) * t337 - t294;
t282 = -t164 * t220 + t224 * t207;
t280 = t312 * t127;
t279 = t124 * t221;
t277 = t197 + t322;
t276 = t313 * qJD(1);
t275 = t196 + t310;
t268 = g(1) * t224 + g(2) * t220;
t8 = t12 * t221 + t217 * t29;
t28 = pkin(10) * t163 + t258;
t136 = t162 * t218 + t216 * t222;
t62 = -pkin(4) * t162 - t85;
t35 = -pkin(5) * t257 - pkin(10) * t136 + t62;
t262 = t217 * t35 + t221 * t28;
t261 = -t217 * t28 + t221 * t35;
t21 = -t218 * t57 + t222 * t40;
t158 = t214 * t244;
t159 = qJD(2) * t299 - t213 * t295;
t249 = pkin(2) * t295 - qJ(4) * t159 - qJD(4) * t163;
t46 = t158 * t375 + t249;
t48 = pkin(4) * t159 + t70;
t260 = -t218 * t46 + t222 * t48;
t259 = -t218 * t74 + t222 * t58;
t95 = t136 * t221 + t163 * t217;
t94 = t136 * t217 - t163 * t221;
t256 = -t164 * t224 - t207 * t220;
t41 = -pkin(4) * t158 - t63;
t254 = -t124 * t314 - t358;
t251 = t218 * t48 + t222 * t46 + t58 * t316 - t317 * t74;
t243 = -g(1) * t119 + g(2) * t114 - g(3) * t163;
t242 = -g(1) * t118 - g(2) * t115 + g(3) * t162;
t238 = t20 + t243;
t237 = -t109 * t218 - t149 * t281;
t236 = t324 * t197;
t11 = -pkin(5) * t149 - t21;
t234 = -pkin(10) * t47 + (t11 + t21) * t124;
t233 = t160 * t70 - t214 * t268;
t2 = -qJD(6) * t8 - t217 * t3 + t221 * t6;
t231 = t254 + t388;
t230 = qJD(6) * t124 * t203 - t243;
t229 = (-pkin(10) * t156 - qJD(6) * t170 + t364) * t124 + t242;
t84 = pkin(3) * t156 + t232;
t228 = t160 * t84 - t242 + t266;
t226 = g(2) * t339 + (-t264 - t367) * t214 - g(3) * t216 + t307;
t181 = pkin(2) * t298;
t169 = -t216 * t335 + t329;
t167 = -t216 * t334 - t332;
t166 = -t298 + t335;
t100 = -t118 * t218 + t222 * t341;
t96 = t240 + t372;
t93 = pkin(3) * t160 + t284;
t90 = qJD(5) * t136 - t158 * t222;
t89 = qJD(5) * t257 + t158 * t218;
t86 = -pkin(3) * t216 - t91;
t69 = pkin(3) * t158 + t249;
t65 = -pkin(3) * t197 + t265;
t52 = t100 * t221 + t119 * t217;
t51 = -t100 * t217 + t119 * t221;
t37 = -qJD(6) * t94 + t159 * t217 + t221 * t89;
t36 = qJD(6) * t95 - t159 * t221 + t217 * t89;
t34 = t227 + t373;
t31 = t266 - t371;
t27 = -pkin(5) * t163 - t259;
t24 = pkin(5) * t156 + t218 * t64 - t222 * t59;
t16 = pkin(5) * t90 - pkin(10) * t89 + t41;
t10 = -pkin(5) * t159 + qJD(5) * t258 - t260;
t9 = pkin(10) * t159 + t251;
t5 = [qJDD(1), -g(2) * t224 + t367, t268 (qJDD(1) * t211 + 0.2e1 * t219 * t292) * t210 (t219 * t308 - t311 * t323) * t376 (qJD(2) * t223 * t277 + t219 * t275) * t214 (t223 * t275 - t277 * t320) * t214, t196 * t216, t302 * t376 + (-pkin(8) * t342 + t201) * t196 + (-t214 * t301 + t191) * t216 - g(1) * t167 - g(2) * t169 + (-t236 + (-t216 * t324 - 0.2e1 * t305) * qJD(1)) * qJD(2) -(-pkin(8) * t295 + t195) * t197 - t324 * t196 - (-pkin(8) * t272 + t297) * t216 - g(1) * t166 - g(2) * t168 + 0.2e1 * (-t292 - t309) * t374, -t110 * t92 - t111 * t91 - t156 * t71 - t158 * t77 - t159 * t76 - t162 * t33 - t163 * t32 + t233, t33 * t92 + t77 * t71 + t32 * t91 - t76 * t70 - g(1) * t256 - g(2) * t282 + (t165 * pkin(2) * t320 + t389 * t207) * t214, t110 * t85 + t111 * t86 + t156 * t63 + t158 * t66 + t159 * t65 - t162 * t278 + t163 * t31 + t233, g(1) * t114 + g(2) * t119 - t110 * t96 - t156 * t69 - t158 * t84 - t162 * t34 + t196 * t86 + t197 * t70 + t216 * t31, -g(1) * t115 + g(2) * t118 - t111 * t96 - t159 * t84 - t160 * t69 - t163 * t34 - t196 * t85 - t197 * t63 + t216 * t278, t34 * t96 + t84 * t69 - t278 * t85 + t66 * t63 + t31 * t86 + t65 * t70 - g(1) * (pkin(3) * t114 + qJ(4) * t115 + t256) - g(2) * (pkin(3) * t119 - qJ(4) * t118 + t282) t127 * t89 + t136 * t49, -t125 * t89 - t127 * t90 - t136 * t50 + t257 * t49, t109 * t136 + t127 * t159 + t149 * t89 + t163 * t49, t109 * t257 - t125 * t159 - t149 * t90 - t163 * t50, t109 * t163 + t149 * t159, t260 * t149 + t259 * t109 - t288 * t163 + t21 * t159 + t41 * t125 + t62 * t50 - t20 * t257 + t45 * t90 - g(1) * t253 - g(2) * t100 + (-t149 * t258 - t163 * t22) * qJD(5), g(1) * t101 - g(2) * t99 - t109 * t258 + t41 * t127 + t20 * t136 - t149 * t251 - t22 * t159 - t163 * t252 + t45 * t89 + t62 * t49, t14 * t95 + t37 * t82, -t14 * t94 - t15 * t95 - t36 * t82 - t37 * t80, t124 * t37 - t14 * t257 + t47 * t95 + t82 * t90, -t124 * t36 + t15 * t257 - t47 * t94 - t80 * t90, t124 * t90 - t257 * t47 (-qJD(6) * t262 + t16 * t221 - t217 * t9) * t124 + t261 * t47 - t2 * t257 - t263 * t90 + t10 * t80 + t27 * t15 + t4 * t94 + t11 * t36 - g(1) * t379 - g(2) * t52 -(qJD(6) * t261 + t16 * t217 + t221 * t9) * t124 - t262 * t47 + t1 * t257 - t8 * t90 + t10 * t82 + t27 * t14 + t4 * t95 + t11 * t37 + g(1) * t380 - g(2) * t51; 0, 0, 0, -t219 * t300, t323 * t343 (t223 * t276 + t309) * t214 (-t313 * t321 + t308) * t214, t196, t225 * t305 - t368 + g(2) * t166 + t191 + (-t301 - t366) * t214 + (-qJD(2) * t324 + t236) * qJD(1), pkin(1) * t300 + g(1) * t169 - g(2) * t167 + t194 * t197 + (pkin(8) * t276 + g(3)) * t342 - t297 (t77 - t87) * t160 + (-t76 + t88) * t156 + (-t110 * t213 - t111 * t215) * pkin(2), -g(2) * t181 + t76 * t87 - t77 * t88 + (t33 * t213 + t32 * t215 - t368 + g(2) * t335 + (-t165 * t321 - t366) * t214) * pkin(2), -t110 * t204 + t111 * t206 + (-t66 - t87) * t160 + (t65 - t326) * t156, t156 * t93 - t197 * t87 + (-pkin(3) + t206) * t196 + t228, -t156 * t84 + t160 * t93 + t196 * t204 + t197 * t326 + t243 + t278, t278 * t204 + t31 * t206 - t84 * t93 - t65 * t87 - g(1) * (pkin(2) * t168 + pkin(3) * t118 + qJ(4) * t119) - g(2) * (-pkin(2) * t335 + pkin(3) * t115 - qJ(4) * t114 + t181) - g(3) * (pkin(2) * t340 + t354 - t372) - t326 * t66, -t218 * t280 + t222 * t49 (-t50 - t280) * t222 + (t125 * t312 - t49) * t218, t250 + t350, t237 - t351, t149 * t156, t21 * t156 + t204 * t50 + t327 * t125 + (-t59 * t149 - t384) * t222 + ((t64 - t318) * t149 + t238) * t218, t204 * t49 + t364 * t149 - t22 * t156 + t327 * t127 + t384 * t218 + (-t149 * t318 + t238) * t222, t14 * t331 - t235 * t82, t80 * t98 + t82 * t97 + (t217 * t82 + t221 * t80) * t317 + (-t361 - t15 * t221 + (t217 * t80 - t221 * t82) * qJD(6)) * t222, t306 + t383, -t360 + (t254 - t388) * t222 + t285, t124 * t281 + t218 * t47, -t11 * t97 - t24 * t80 + t381 * t221 + t229 * t217 + (-t47 * t345 + t2 + (-t11 * t217 + t203 * t80) * qJD(5) - t230 * t221) * t218 + (t11 * t314 - t203 * t15 - t263 * t160 + t4 * t217 + (-t124 * t345 - t263) * qJD(5)) * t222, -t11 * t98 - t24 * t82 - t381 * t217 + t229 * t221 + (-t47 * t344 - t1 + (-t11 * t221 + t203 * t82) * qJD(5) + t230 * t217) * t218 + (-t11 * t315 - t203 * t14 - t8 * t160 + t4 * t221 + (-t124 * t344 - t8) * qJD(5)) * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386, t156 * t77 + t160 * t76 + t226, t386, -t160 * t197 + t175 + (-qJD(2) * t245 - t291) * t214, -t111 + t348, t373 - t355 - t156 * t66 + (-qJD(4) - t65) * t160 + t226, 0, 0, 0, 0, 0, t237 + t351, -t250 + t350, 0, 0, 0, 0, 0, t222 * t231 + t285 + t360, t306 - t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111 + t348, -t156 * t160 + t196, -t197 ^ 2 - t150, t197 * t66 + t228 - t371, 0, 0, 0, 0, 0, -t125 * t197 + t250, -t127 * t197 + t237, 0, 0, 0, 0, 0, -t222 * t15 + (-t221 * t197 - t217 * t281) * t124 + t231 * t218, -t222 * t14 + (t217 * t197 - t221 * t281) * t124 + (t124 * t315 + t312 * t82 - t357) * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 * t125, -t125 ^ 2 + t127 ^ 2, t125 * t149 + t49, t127 * t385 - t283, t109, -t127 * t45 + t22 * t385 - t246 - t288, g(1) * t100 - g(2) * t253 + g(3) * t136 + t125 * t45 + t149 * t21 - t252, t279 * t82 + t361 (t14 - t363) * t221 + (-t15 - t362) * t217, t124 * t279 - t127 * t82 + t358, -t124 ^ 2 * t217 + t127 * t80 + t357, -t124 * t127, -pkin(5) * t15 + t127 * t263 + t234 * t217 - t22 * t80 - t221 * t382, -pkin(5) * t14 + t8 * t127 + t217 * t382 - t22 * t82 + t234 * t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t80, -t80 ^ 2 + t82 ^ 2, t14 + t363, -t15 + t362, t47, -g(1) * t51 - g(2) * t380 + g(3) * t94 - t11 * t82 + t8 * t124 + t2, g(1) * t52 - g(2) * t379 + g(3) * t95 + t11 * t80 - t124 * t263 - t1;];
tau_reg  = t5;
