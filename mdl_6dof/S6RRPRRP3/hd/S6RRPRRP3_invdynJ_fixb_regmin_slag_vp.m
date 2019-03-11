% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:38
% EndTime: 2019-03-09 11:50:51
% DurationCPUTime: 5.30s
% Computational Cost: add. (8765->477), mult. (20461->622), div. (0->0), fcn. (15484->14), ass. (0->253)
t226 = sin(qJ(4));
t222 = sin(pkin(10));
t202 = pkin(2) * t222 + pkin(8);
t341 = pkin(9) + t202;
t270 = qJD(4) * t341;
t227 = sin(qJ(2));
t298 = qJD(1) * t227;
t223 = cos(pkin(10));
t230 = cos(qJ(2));
t308 = t223 * t230;
t167 = qJD(1) * t308 - t222 * t298;
t319 = t167 * t226;
t340 = qJ(3) + pkin(7);
t193 = t340 * t230;
t184 = qJD(1) * t193;
t172 = t222 * t184;
t192 = t340 * t227;
t183 = qJD(1) * t192;
t124 = -t183 * t223 - t172;
t229 = cos(qJ(4));
t179 = t222 * t230 + t223 * t227;
t169 = t179 * qJD(1);
t98 = pkin(2) * t298 + pkin(3) * t169 - pkin(8) * t167;
t328 = t229 * t124 + t226 * t98;
t378 = -pkin(9) * t319 + t226 * t270 + t328;
t91 = t229 * t98;
t377 = pkin(4) * t169 - t124 * t226 + t91 + (-pkin(9) * t167 + t270) * t229;
t133 = qJD(2) * t229 - t226 * t169;
t134 = qJD(2) * t226 + t169 * t229;
t225 = sin(qJ(5));
t353 = cos(qJ(5));
t253 = t225 * t133 + t134 * t353;
t354 = t253 ^ 2;
t75 = -t353 * t133 + t134 * t225;
t73 = t75 ^ 2;
t376 = -t73 + t354;
t297 = qJD(4) * t226;
t375 = t297 - t319;
t374 = t253 * t75;
t373 = t75 * qJ(6);
t307 = t225 * t226;
t252 = t353 * t229 - t307;
t361 = qJD(4) + qJD(5);
t279 = t353 * qJD(5);
t362 = t353 * qJD(4) + t279;
t327 = t252 * t167 - t229 * t362 + t307 * t361;
t286 = t353 * t226;
t182 = t225 * t229 + t286;
t127 = t361 * t182;
t326 = -t182 * t167 + t127;
t294 = qJD(1) * qJD(2);
t278 = t227 * t294;
t194 = t222 * t278;
t277 = t230 * t294;
t257 = t223 * t277 - t194;
t235 = qJDD(1) * t179 + t257;
t372 = qJD(2) * qJD(4) + t235;
t178 = t222 * t227 - t308;
t171 = t178 * qJD(2);
t296 = qJD(4) * t229;
t283 = t179 * t296;
t371 = -t171 * t226 + t283;
t157 = qJD(4) - t167;
t151 = qJD(5) + t157;
t287 = t169 * t296 + t226 * t372;
t256 = t229 * qJDD(2) - t287;
t295 = qJD(5) * t225;
t71 = t226 * qJDD(2) - t169 * t297 + t229 * t372;
t26 = -t133 * t279 + t134 * t295 - t225 * t256 - t353 * t71;
t370 = t151 * t75 - t26;
t218 = qJ(2) + pkin(10);
t211 = cos(t218);
t221 = qJ(4) + qJ(5);
t214 = cos(t221);
t228 = sin(qJ(1));
t311 = t214 * t228;
t213 = sin(t221);
t231 = cos(qJ(1));
t312 = t213 * t231;
t137 = -t211 * t311 + t312;
t310 = t214 * t231;
t313 = t213 * t228;
t139 = t211 * t310 + t313;
t210 = sin(t218);
t333 = qJD(2) * pkin(2);
t175 = -t183 + t333;
t309 = t223 * t184;
t120 = t222 * t175 + t309;
t108 = qJD(2) * pkin(8) + t120;
t291 = t230 * qJDD(1);
t292 = t227 * qJDD(1);
t260 = -t222 * t292 + t223 * t291;
t121 = -qJD(2) * t169 + t260;
t290 = pkin(2) * t278 + qJDD(3);
t342 = t230 * pkin(2);
t209 = pkin(1) + t342;
t366 = -pkin(8) * t179 - t209;
t58 = -t121 * pkin(3) - t257 * pkin(8) + qJDD(1) * t366 + t290;
t271 = qJD(2) * t340;
t166 = -t227 * qJD(3) - t230 * t271;
t115 = qJDD(2) * pkin(2) + qJD(1) * t166 - qJDD(1) * t192;
t165 = t230 * qJD(3) - t227 * t271;
t122 = qJD(1) * t165 + qJDD(1) * t193;
t67 = t222 * t115 + t223 * t122;
t65 = qJDD(2) * pkin(8) + t67;
t190 = -qJD(1) * t209 + qJD(3);
t87 = -pkin(3) * t167 - pkin(8) * t169 + t190;
t249 = -t108 * t297 + t226 * t58 + t229 * t65 + t87 * t296;
t10 = pkin(9) * t256 + t249;
t56 = -t108 * t226 + t229 * t87;
t46 = -pkin(9) * t134 + t56;
t33 = pkin(4) * t157 + t46;
t57 = t108 * t229 + t226 * t87;
t47 = pkin(9) * t133 + t57;
t168 = t179 * qJD(2);
t116 = qJD(1) * t168 + qJDD(4) - t260;
t55 = t229 * t58;
t7 = pkin(4) * t116 - pkin(9) * t71 - qJD(4) * t57 - t226 * t65 + t55;
t281 = -t353 * t10 - t225 * t7 - t33 * t279 + t47 * t295;
t344 = g(3) * t214;
t119 = t175 * t223 - t172;
t107 = -qJD(2) * pkin(3) - t119;
t72 = -t133 * pkin(4) + t107;
t369 = g(1) * t139 - g(2) * t137 + t210 * t344 + t72 * t75 + t281;
t38 = pkin(5) * t75 + qJD(6) + t72;
t368 = t253 * t38;
t123 = -t222 * t183 + t309;
t264 = pkin(4) * t375 - t123;
t367 = qJ(6) * t253;
t176 = t341 * t226;
t177 = t341 * t229;
t301 = -t225 * t176 + t353 * t177;
t365 = g(1) * t228 - g(2) * t231;
t364 = -qJD(5) * t301 + t225 * t378 - t377 * t353;
t363 = t176 * t279 + t177 * t295 + t377 * t225 + t353 * t378;
t263 = g(1) * t231 + g(2) * t228;
t242 = g(3) * t211 - t210 * t263;
t66 = t223 * t115 - t222 * t122;
t64 = -qJDD(2) * pkin(3) - t66;
t360 = qJD(4) * t202 * t157 + t242 + t64;
t136 = t211 * t313 + t310;
t138 = -t211 * t312 + t311;
t346 = g(3) * t210;
t359 = -g(1) * t138 + g(2) * t136 + t213 * t346;
t44 = t353 * t47;
t19 = t225 * t33 + t44;
t239 = -qJD(5) * t19 - t225 * t10 + t353 * t7;
t358 = -t72 * t253 + t239 + t359;
t27 = qJD(5) * t253 + t225 * t71 - t353 * t256;
t357 = t151 * t253 - t27;
t111 = qJDD(5) + t116;
t356 = t111 * t182 - t151 * t327;
t355 = -t252 * t26 - t253 * t326;
t343 = g(3) * t230;
t42 = t225 * t47;
t18 = t353 * t33 - t42;
t12 = t18 - t367;
t11 = pkin(5) * t151 + t12;
t339 = -t12 + t11;
t338 = -qJ(6) * t326 + qJD(6) * t252 - t363;
t337 = -pkin(5) * t169 + qJ(6) * t327 - t182 * qJD(6) + t364;
t336 = t353 * t46 - t42;
t118 = pkin(3) * t178 + t366;
t106 = t229 * t118;
t132 = -t192 * t222 + t193 * t223;
t315 = t179 * t229;
t52 = pkin(4) * t178 - pkin(9) * t315 - t132 * t226 + t106;
t125 = t229 * t132;
t302 = t226 * t118 + t125;
t316 = t179 * t226;
t61 = -pkin(9) * t316 + t302;
t334 = t225 * t52 + t353 * t61;
t332 = t169 * t75;
t331 = t169 * t253;
t329 = t71 * t226;
t324 = t116 * t226;
t323 = t133 * t157;
t322 = t133 * t169;
t321 = t134 * t157;
t320 = t134 * t169;
t317 = t171 * t229;
t186 = pkin(4) * t226 + pkin(5) * t213;
t314 = t186 * t211;
t306 = t226 * t228;
t305 = t226 * t231;
t304 = t228 * t229;
t102 = t229 * t116;
t303 = t229 * t231;
t300 = t186 + t340;
t215 = t229 * pkin(4);
t187 = pkin(5) * t214 + t215;
t219 = t227 ^ 2;
t299 = -t230 ^ 2 + t219;
t289 = t227 * t333;
t203 = -pkin(2) * t223 - pkin(3);
t284 = t179 * t297;
t274 = -t225 * t46 - t44;
t273 = -t225 * t61 + t353 * t52;
t269 = -qJD(4) * t87 - t65;
t96 = t165 * t222 - t223 * t166;
t268 = -t353 * t176 - t177 * t225;
t131 = t223 * t192 + t193 * t222;
t267 = t157 * t229;
t191 = t203 - t215;
t266 = -t182 * t27 + t327 * t75;
t265 = t252 * t111 - t151 * t326;
t261 = -t108 * t296 + t55;
t95 = pkin(4) * t316 + t131;
t185 = pkin(3) + t187;
t217 = -qJ(6) - pkin(9) - pkin(8);
t259 = t185 * t211 - t210 * t217;
t63 = pkin(4) * t371 + t96;
t258 = -t157 * t375 + t102;
t255 = -0.2e1 * pkin(1) * t294 - pkin(7) * qJDD(2);
t99 = pkin(3) * t168 + pkin(8) * t171 + t289;
t92 = t229 * t99;
t97 = t165 * t223 + t166 * t222;
t25 = pkin(9) * t317 + pkin(4) * t168 - t226 * t97 + t92 + (-t125 + (pkin(9) * t179 - t118) * t226) * qJD(4);
t248 = t118 * t296 - t132 * t297 + t226 * t99 + t229 * t97;
t29 = -pkin(9) * t371 + t248;
t254 = t225 * t25 + t52 * t279 + t353 * t29 - t295 * t61;
t250 = -t284 - t317;
t246 = t107 * t157 - t202 * t116;
t244 = -qJDD(1) * t209 + t290;
t232 = qJD(2) ^ 2;
t241 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t232 + t365;
t233 = qJD(1) ^ 2;
t240 = pkin(1) * t233 - pkin(7) * qJDD(1) + t263;
t238 = -qJD(5) * t334 - t225 * t29 + t353 * t25;
t30 = -pkin(4) * t256 + t64;
t8 = t27 * pkin(5) + qJDD(6) + t30;
t208 = pkin(4) * t353 + pkin(5);
t195 = t231 * t209;
t161 = t211 * t303 + t306;
t160 = -t211 * t305 + t304;
t159 = -t211 * t304 + t305;
t158 = t211 * t306 + t303;
t104 = t252 * t179;
t103 = t182 * t179;
t84 = qJ(6) * t252 + t301;
t83 = -qJ(6) * t182 + t268;
t35 = -t171 * t286 - t225 * t284 - t295 * t316 + (-t171 * t225 + t179 * t362) * t229;
t34 = t127 * t179 + t171 * t252;
t21 = -qJ(6) * t103 + t334;
t20 = pkin(5) * t178 - qJ(6) * t104 + t273;
t15 = t336 - t367;
t14 = t274 + t373;
t13 = t19 - t373;
t4 = -qJ(6) * t35 - qJD(6) * t103 + t254;
t3 = t168 * pkin(5) + t34 * qJ(6) - t104 * qJD(6) + t238;
t2 = -qJ(6) * t27 - qJD(6) * t75 - t281;
t1 = t111 * pkin(5) + t26 * qJ(6) - qJD(6) * t253 + t239;
t5 = [qJDD(1), t365, t263, qJDD(1) * t219 + 0.2e1 * t227 * t277, 0.2e1 * t227 * t291 - 0.2e1 * t294 * t299, qJDD(2) * t227 + t230 * t232, qJDD(2) * t230 - t227 * t232, 0, t227 * t255 + t230 * t241, -t227 * t241 + t230 * t255, t119 * t171 - t120 * t168 + t132 * t121 + t131 * t235 + t97 * t167 + t96 * t169 - t67 * t178 - t66 * t179 - t263, t67 * t132 + t120 * t97 - t66 * t131 - t119 * t96 - t244 * t209 + t190 * t289 - g(1) * (-t209 * t228 + t231 * t340) - g(2) * (t228 * t340 + t195) t134 * t250 + t315 * t71 -(t133 * t229 - t134 * t226) * t171 + (t229 * t256 - t329 + (-t133 * t226 - t134 * t229) * qJD(4)) * t179, t102 * t179 + t134 * t168 + t157 * t250 + t178 * t71, -t116 * t316 + t133 * t168 - t157 * t371 + t178 * t256, t116 * t178 + t157 * t168 (-t132 * t296 + t92) * t157 + t106 * t116 + t261 * t178 + t56 * t168 - t96 * t133 - t131 * t256 + t107 * t283 - g(1) * t159 - g(2) * t161 + ((-qJD(4) * t118 - t97) * t157 - t132 * t116 + t269 * t178 + t64 * t179 - t107 * t171) * t226, -g(1) * t158 - g(2) * t160 + t107 * t250 - t116 * t302 + t131 * t71 + t96 * t134 - t157 * t248 - t57 * t168 - t178 * t249 + t315 * t64, -t104 * t26 - t253 * t34, t103 * t26 - t104 * t27 - t253 * t35 + t34 * t75, t104 * t111 - t151 * t34 + t168 * t253 - t178 * t26, -t103 * t111 - t151 * t35 - t168 * t75 - t178 * t27, t111 * t178 + t151 * t168, -g(1) * t137 - g(2) * t139 + t30 * t103 + t111 * t273 + t151 * t238 + t18 * t168 + t178 * t239 + t95 * t27 + t72 * t35 + t63 * t75, -g(1) * t136 - g(2) * t138 + t30 * t104 - t111 * t334 - t151 * t254 - t19 * t168 + t178 * t281 + t253 * t63 - t95 * t26 - t72 * t34, -t1 * t104 - t103 * t2 + t11 * t34 - t13 * t35 + t20 * t26 - t21 * t27 + t210 * t365 - t253 * t3 - t4 * t75, t2 * t21 + t13 * t4 + t1 * t20 + t11 * t3 + t8 * (pkin(5) * t103 + t95) + t38 * (pkin(5) * t35 + t63) - g(2) * t195 + (-g(1) * t300 - g(2) * t259) * t231 + (-g(1) * (-t209 - t259) - g(2) * t300) * t228; 0, 0, 0, -t227 * t233 * t230, t299 * t233, t292, t291, qJDD(2), t227 * t240 - t343, g(3) * t227 + t230 * t240 (t120 - t123) * t169 + (-t124 + t119) * t167 + (t222 * t121 + (-t222 * t291 + t194 + (-t277 - t292) * t223) * t223) * pkin(2), t119 * t123 - t120 * t124 + (-t343 + t222 * t67 + t223 * t66 + (-qJD(1) * t190 + t263) * t227) * pkin(2), t134 * t267 + t329 (t71 + t323) * t229 + (t256 - t321) * t226, t157 * t267 - t320 + t324, t258 - t322, -t157 * t169, t203 * t287 - t91 * t157 - t56 * t169 + t123 * t133 + (t124 * t157 + t246) * t226 + (-t203 * qJDD(2) - t360) * t229, -t123 * t134 + t328 * t157 + t57 * t169 + t203 * t71 + t226 * t360 + t246 * t229, -t26 * t182 - t253 * t327, t266 + t355, -t331 + t356, t265 + t332, -t151 * t169, t111 * t268 - t18 * t169 - t30 * t252 + t191 * t27 - t211 * t344 + t264 * t75 + t326 * t72 + (g(1) * t310 + g(2) * t311) * t210 + t364 * t151, -t301 * t111 + t151 * t363 + t19 * t169 + t30 * t182 - t191 * t26 + t242 * t213 + t264 * t253 - t327 * t72, -t1 * t182 + t11 * t327 - t13 * t326 + t2 * t252 - t211 * t263 - t253 * t337 + t26 * t83 - t27 * t84 - t338 * t75 - t346, t2 * t84 + t1 * t83 + t8 * (-pkin(5) * t252 + t191) - g(3) * (t259 + t342) + (pkin(5) * t326 + t264) * t38 + t338 * t13 + t337 * t11 + t263 * (pkin(2) * t227 + t185 * t210 + t211 * t217); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167 ^ 2 - t169 ^ 2, t119 * t169 - t120 * t167 + t244 - t365, 0, 0, 0, 0, 0, t258 + t322, -t157 ^ 2 * t229 - t320 - t324, 0, 0, 0, 0, 0, t265 - t332, -t331 - t356, t266 - t355, t1 * t252 - t11 * t326 - t13 * t327 - t169 * t38 + t182 * t2 - t365; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t133, -t133 ^ 2 + t134 ^ 2, t71 - t323, t256 + t321, t116, -g(1) * t160 + g(2) * t158 - t107 * t134 + t157 * t57 + (t269 + t346) * t226 + t261, g(1) * t161 - g(2) * t159 - t107 * t133 + t157 * t56 + t229 * t346 - t249, t374, t376, t370, t357, t111, -t274 * t151 + (t111 * t353 - t134 * t75 - t151 * t295) * pkin(4) + t358, t336 * t151 + (-t225 * t111 - t134 * t253 - t151 * t279) * pkin(4) + t369, -t11 * t75 + t13 * t253 + t14 * t253 + t15 * t75 + t208 * t26 + (-t225 * t27 + (t225 * t253 - t353 * t75) * qJD(5)) * pkin(4), t1 * t208 - t13 * t15 - t11 * t14 - pkin(5) * t368 - g(1) * (t187 * t228 - t231 * t314) - g(2) * (-t187 * t231 - t228 * t314) + t186 * t346 + (-t38 * t134 + t2 * t225 + (-t11 * t225 + t13 * t353) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t374, t376, t370, t357, t111, t19 * t151 + t358, t151 * t18 + t369, pkin(5) * t26 - t339 * t75, t339 * t13 + (t1 + t359 - t368) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 - t354, t11 * t253 + t13 * t75 + t242 + t8;];
tau_reg  = t5;
