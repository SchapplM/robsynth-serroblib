% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:11:03
% EndTime: 2019-03-08 21:11:12
% DurationCPUTime: 5.10s
% Computational Cost: add. (3651->518), mult. (7995->652), div. (0->0), fcn. (5617->10), ass. (0->264)
t176 = sin(qJ(2));
t169 = sin(pkin(6));
t288 = qJD(1) * t169;
t100 = qJD(2) * pkin(8) + t176 * t288;
t171 = cos(pkin(6));
t178 = cos(qJ(3));
t287 = qJD(1) * t178;
t124 = t171 * t287;
t175 = sin(qJ(3));
t241 = t175 * t100 - t124;
t340 = qJD(4) + t241;
t52 = -qJD(3) * pkin(3) + t340;
t284 = qJD(2) * t178;
t151 = t178 * qJDD(2);
t174 = sin(qJ(6));
t177 = cos(qJ(6));
t277 = qJD(6) * t178;
t259 = t174 * t277;
t276 = t177 * qJD(3);
t29 = -t177 * t151 - qJD(6) * t276 - t174 * qJDD(3) + (t175 * t276 + t259) * qJD(2);
t286 = qJD(2) * t175;
t133 = qJD(6) + t286;
t94 = t174 * t284 - t276;
t313 = t94 * t133;
t336 = t29 - t313;
t275 = qJD(2) * qJD(3);
t254 = t175 * t275;
t282 = qJD(3) * t174;
t95 = t177 * t284 + t282;
t30 = qJD(6) * t95 - t177 * qJDD(3) + (t151 - t254) * t174;
t312 = t95 * t133;
t211 = t30 - t312;
t342 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t253 = t178 * t275;
t269 = t175 * qJDD(2);
t341 = t253 + t269;
t327 = pkin(8) - qJ(5);
t232 = pkin(5) * t175 + pkin(9) * t178;
t179 = cos(qJ(2));
t123 = t179 * t288;
t172 = qJD(2) * pkin(2);
t101 = -t123 - t172;
t62 = -pkin(3) * t284 - qJ(4) * t286 + t101;
t40 = pkin(4) * t284 + qJD(5) - t62;
t25 = qJD(2) * t232 + t40;
t180 = -pkin(3) - pkin(4);
t164 = -pkin(9) + t180;
t138 = qJ(5) * t286;
t217 = -t138 + t340;
t28 = qJD(3) * t164 + t217;
t227 = t174 * t28 - t177 * t25;
t153 = t175 * qJD(4);
t285 = qJD(2) * t176;
t261 = t169 * t285;
t116 = qJD(1) * t261;
t301 = t169 * t179;
t121 = qJDD(1) * t301;
t168 = qJDD(2) * pkin(2);
t64 = t116 - t121 - t168;
t207 = pkin(3) * t151 + qJ(4) * t341 + qJD(2) * t153 - t64;
t199 = pkin(4) * t151 + qJDD(5) + t207;
t212 = pkin(5) * t178 + t164 * t175;
t201 = t212 * qJD(3);
t6 = qJD(2) * t201 + qJDD(2) * t232 + t199;
t281 = qJD(3) * t175;
t255 = qJD(1) * t281;
t273 = qJDD(1) * t171;
t280 = qJD(3) * t178;
t283 = qJD(2) * t179;
t256 = qJD(1) * t283;
t272 = qJDD(1) * t176;
t65 = qJDD(2) * pkin(8) + (t256 + t272) * t169;
t239 = t100 * t280 + t171 * t255 + t175 * t65 - t178 * t273;
t224 = -qJDD(4) - t239;
t274 = qJD(2) * qJD(5);
t184 = -qJ(5) * t341 - t175 * t274 - t224;
t7 = qJDD(3) * t164 + t184;
t1 = -t227 * qJD(6) + t174 * t6 + t177 * t7;
t339 = t133 * t227 + t1;
t9 = t174 * t25 + t177 * t28;
t2 = -qJD(6) * t9 - t174 * t7 + t177 * t6;
t338 = -t9 * t133 - t2;
t165 = qJD(3) * qJ(4);
t298 = t171 * t175;
t258 = qJD(1) * t298;
t49 = -t178 * (qJ(5) * qJD(2) - t100) + t258;
t337 = t165 + t49;
t59 = t100 * t178 + t258;
t53 = t165 + t59;
t335 = 0.2e1 * t342;
t262 = qJD(3) * t180;
t257 = qJDD(3) * t180;
t127 = qJ(5) * t254;
t240 = -qJD(3) * t124 + t100 * t281 - t175 * t273 - t178 * t65;
t20 = -t240 + t342;
t12 = t178 * (qJ(5) * qJDD(2) + t274) - t127 - t20;
t260 = t169 * t283;
t302 = t169 * t178;
t89 = t176 * t302 + t298;
t47 = qJD(3) * t89 + t175 * t260;
t303 = t169 * t176;
t264 = t175 * t303;
t88 = -t171 * t178 + t264;
t334 = -(t175 * t89 - t178 * t88) * qJD(3) + t175 * t47;
t132 = g(3) * t303;
t167 = t178 ^ 2;
t270 = qJDD(2) * t167;
t166 = t175 ^ 2;
t271 = qJDD(2) * t166;
t289 = t166 + t167;
t333 = -t289 * t169 * t256 - t132 + (t270 + t271) * pkin(8);
t296 = t177 * t179;
t203 = t169 * (-t174 * t176 + t175 * t296);
t111 = t327 * t175;
t155 = t175 * qJ(4);
t102 = -pkin(3) * t178 - pkin(2) - t155;
t90 = pkin(4) * t178 - t102;
t61 = t232 + t90;
t26 = -t111 * t174 + t177 * t61;
t290 = qJ(4) * t280 + t153;
t33 = t201 + t290;
t80 = -t175 * qJD(5) + t280 * t327;
t332 = -qJD(1) * t203 + qJD(6) * t26 + t174 * t33 + t177 * t80;
t297 = t175 * t179;
t202 = t169 * (-t174 * t297 - t176 * t177);
t27 = t111 * t177 + t174 * t61;
t331 = -qJD(1) * t202 - qJD(6) * t27 - t174 * t80 + t177 * t33;
t328 = t94 * t95;
t173 = qJ(4) + pkin(5);
t326 = t116 * t178 + t255 * t301;
t93 = qJDD(6) + t341;
t324 = t174 * t93;
t322 = t177 * t93;
t321 = t177 * t94;
t320 = t178 * t337;
t311 = sin(pkin(10));
t247 = t311 * t176;
t170 = cos(pkin(10));
t299 = t170 * t179;
t82 = t171 * t299 - t247;
t319 = t178 * t82;
t246 = t311 * t179;
t300 = t170 * t176;
t84 = -t171 * t246 - t300;
t318 = t178 * t84;
t317 = t178 * t94;
t316 = t178 * t95;
t315 = t29 * t174;
t314 = t30 * t177;
t310 = pkin(8) * qJDD(3);
t309 = qJ(4) * t178;
t308 = qJD(3) * t94;
t307 = qJD(3) * t95;
t306 = qJDD(3) * pkin(3);
t220 = t133 * t174;
t305 = t133 * t177;
t182 = qJD(2) ^ 2;
t304 = t167 * t182;
t295 = t337 * qJD(3);
t48 = -t138 + t241;
t294 = qJD(4) + t48;
t292 = -qJD(5) - t40;
t291 = pkin(2) * t301 + pkin(8) * t303;
t279 = qJD(6) * t174;
t278 = qJD(6) * t177;
t75 = t82 * pkin(2);
t268 = pkin(3) * t319 + t155 * t82 + t75;
t76 = t84 * pkin(2);
t267 = pkin(3) * t318 + t155 * t84 + t76;
t266 = t175 * t220;
t265 = t175 * t305;
t263 = t178 * t301;
t252 = t179 * t275;
t83 = t171 * t300 + t246;
t41 = t170 * t302 + t175 * t83;
t37 = t41 * pkin(3);
t42 = -t169 * t170 * t175 + t178 * t83;
t251 = qJ(4) * t42 - t37;
t248 = t169 * t311;
t85 = -t171 * t247 + t299;
t43 = t175 * t85 - t178 * t248;
t39 = t43 * pkin(3);
t44 = t175 * t248 + t178 * t85;
t250 = qJ(4) * t44 - t39;
t78 = t88 * pkin(3);
t249 = qJ(4) * t89 - t78;
t245 = t101 - t172;
t244 = qJD(2) * t90 + t40;
t243 = t292 * t175;
t242 = qJD(2) * t102 + t62;
t237 = qJ(4) * t169 * t297 + pkin(3) * t263 + t291;
t236 = t175 * t262;
t235 = t175 * t253;
t234 = g(1) * t84 + g(2) * t82;
t233 = g(1) * t85 + g(2) * t83;
t230 = -t174 * t9 + t177 * t227;
t229 = t174 * t227 + t177 * t9;
t225 = pkin(4) * t263 + t237;
t221 = t132 + t233;
t219 = qJDD(2) * t179 - t176 * t182;
t50 = t169 * t296 - t174 * t88;
t51 = t174 * t301 + t177 * t88;
t216 = -t133 * t278 - t324;
t215 = t133 * t279 - t322;
t214 = g(1) * t43 + g(2) * t41 + g(3) * t88;
t213 = g(1) * t44 + g(2) * t42 + g(3) * t89;
t206 = pkin(4) * t319 + t327 * t83 + t268;
t205 = pkin(4) * t318 + t327 * t85 + t267;
t204 = g(3) * t301 + t234;
t32 = qJD(3) * pkin(5) + t337;
t200 = -t133 * t32 - t164 * t93;
t181 = qJD(3) ^ 2;
t198 = pkin(8) * t181 + t204;
t197 = t214 - t239;
t196 = -t213 - t240;
t10 = qJDD(3) * pkin(5) - t12;
t195 = qJD(6) * t133 * t164 - t10 + t213;
t194 = -qJDD(4) + t197;
t193 = qJD(3) * t59 + t197;
t192 = qJD(3) * t241 + t196;
t191 = t198 + t64 - t168;
t17 = qJD(2) * t236 + t199;
t63 = t236 + t290;
t190 = -qJD(2) * t63 - qJDD(2) * t90 - t17 + t204;
t189 = -qJ(5) * t269 - t194;
t24 = pkin(3) * t254 - t207;
t81 = pkin(3) * t281 - t290;
t188 = -qJD(2) * t81 - qJDD(2) * t102 - t198 - t24;
t21 = -t224 - t306;
t187 = t21 * t175 + t20 * t178 + (-t175 * t53 + t178 * t52) * qJD(3) - t233;
t186 = t239 * t175 - t240 * t178 + (-t175 * t59 + t178 * t241) * qJD(3) - t233;
t46 = -qJD(3) * t264 + (qJD(3) * t171 + t260) * t178;
t185 = t334 * qJD(2) + t151 * t89 + t269 * t88 + t284 * t46;
t183 = qJD(6) * t229 + t1 * t174 + t2 * t177 - t204;
t14 = t46 * qJD(3) + t89 * qJDD(3) + (t175 * t219 + t178 * t252) * t169;
t13 = -t47 * qJD(3) - t88 * qJDD(3) + (-t175 * t252 + t178 * t219) * t169;
t152 = t166 * t182;
t140 = qJ(4) * t284;
t129 = t175 * t182 * t178;
t125 = -t152 - t181;
t112 = t327 * t178;
t110 = qJDD(3) + t129;
t106 = -t152 + t304;
t104 = qJDD(3) * t178 - t175 * t181;
t103 = qJDD(3) * t175 + t178 * t181;
t98 = t175 * t116;
t96 = pkin(3) * t286 - t140;
t92 = -0.2e1 * t235 + t270;
t91 = 0.2e1 * t235 + t271;
t79 = qJD(5) * t178 + t281 * t327;
t77 = t88 * pkin(4);
t74 = t180 * t286 + t140;
t69 = t175 * t151 + (-t166 + t167) * t275;
t60 = 0.2e1 * t69;
t45 = qJD(2) * t212 + t140;
t38 = t43 * pkin(4);
t36 = t41 * pkin(4);
t31 = t262 + t217;
t19 = t174 * t45 + t177 * t49;
t18 = -t174 * t49 + t177 * t45;
t16 = qJD(6) * t50 - t174 * t261 + t47 * t177;
t15 = -qJD(6) * t51 - t47 * t174 - t177 * t261;
t11 = t257 + t184;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, t219 * t169 (-qJDD(2) * t176 - t179 * t182) * t169, 0, -g(3) + (t171 ^ 2 + (t176 ^ 2 + t179 ^ 2) * t169 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t13, -t14, t185, -t240 * t89 + t239 * t88 + t59 * t46 + t241 * t47 - g(3) + (t101 * t285 - t179 * t64) * t169, 0, 0, 0, 0, 0, 0, t13, t185, t14, t20 * t89 + t21 * t88 + t53 * t46 + t52 * t47 - g(3) + (-t179 * t24 + t285 * t62) * t169, 0, 0, 0, 0, 0, 0, t14, -t13 (-t175 * t88 - t178 * t89) * qJDD(2) + (-t178 * t46 - t334) * qJD(2), t11 * t88 - t12 * t89 + t31 * t47 + t337 * t46 - g(3) + (t17 * t179 - t285 * t40) * t169, 0, 0, 0, 0, 0, 0, t133 * t15 - t30 * t89 - t46 * t94 + t50 * t93, -t133 * t16 + t29 * t89 - t46 * t95 - t51 * t93, t15 * t95 + t16 * t94 - t29 * t50 + t30 * t51, t1 * t51 + t10 * t89 - t15 * t227 + t16 * t9 + t2 * t50 + t32 * t46 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t121 - t204, -t169 * t272 + t221, 0, 0, t91, t60, t103, t92, t104, 0 (qJD(3) * t245 - t310) * t175 - t191 * t178 + t326, -t98 + (-t310 + (t245 + t123) * qJD(3)) * t178 + t191 * t175, t186 + t333, -t64 * pkin(2) - g(1) * t76 - g(2) * t75 - g(3) * t291 + (-t101 * t176 + (-t175 * t241 - t178 * t59) * t179) * t288 + t186 * pkin(8), t91, t103, -0.2e1 * t69, 0, -t104, t92 (qJD(3) * t242 - t310) * t175 + t188 * t178 + t326, t187 + t333, t98 + (t310 + (-t242 - t123) * qJD(3)) * t178 + t188 * t175, t24 * t102 + t62 * t81 - g(1) * t267 - g(2) * t268 - g(3) * t237 + (-t176 * t62 + (-t175 * t52 - t178 * t53) * t179) * t288 + t187 * pkin(8), t92, t60, t104, t91, t103, 0, t112 * qJDD(3) + t98 + (-t79 + (t244 - t123) * t178) * qJD(3) - t190 * t175, t111 * qJDD(3) + (t175 * t244 + t80) * qJD(3) + t190 * t178 - t326 (-qJD(3) * t31 - qJDD(2) * t112 + t12) * t178 + (-qJDD(2) * t111 - t11 + t295) * t175 + (-t175 * t80 + t178 * t79 + (-t111 * t178 + t112 * t175) * qJD(3) + t289 * t123) * qJD(2) + t221, t11 * t111 + t31 * t80 - t12 * t112 - t337 * t79 + t17 * t90 + t40 * t63 - g(1) * t205 - g(2) * t206 - g(3) * (-qJ(5) * t303 + t225) + (t176 * t40 + (-t175 * t31 - t320) * t179) * t288, -t95 * t259 + (-t178 * t29 - t281 * t95) * t177 (t174 * t95 + t321) * t281 + (t315 - t314 + (t174 * t94 - t177 * t95) * qJD(6)) * t178 (t133 * t276 + t29) * t175 + (t215 - t307) * t178, t277 * t321 + (t178 * t30 - t281 * t94) * t174 (-t133 * t282 + t30) * t175 + (-t216 + t308) * t178, t133 * t280 + t175 * t93, -t112 * t30 + t26 * t93 + t79 * t94 + t233 * t174 + (-t177 * t234 + t282 * t32 + t2) * t175 - g(3) * t203 + t331 * t133 + (-qJD(3) * t227 - t10 * t174 + t123 * t94 - t278 * t32) * t178, t112 * t29 - t27 * t93 + t79 * t95 + t233 * t177 + (t174 * t234 + t276 * t32 - t1) * t175 - g(3) * t202 - t332 * t133 + (-qJD(3) * t9 - t10 * t177 + t123 * t95 + t279 * t32) * t178, t178 * t183 + t230 * t281 - t26 * t29 + t27 * t30 + t331 * t95 + t332 * t94, t1 * t27 + t2 * t26 + t10 * t112 - t32 * t79 - g(1) * (t232 * t84 + t205) - g(2) * (t232 * t82 + t206) - g(3) * t225 + t332 * t9 - t331 * t227 + (g(3) * qJ(5) * t176 + (-g(3) * t232 - t287 * t32) * t179) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t106, t269, t129, t151, qJDD(3), -t101 * t286 + t193, -t101 * t284 - t192, 0, 0, -t129, t269, t106, qJDD(3), -t151, t129, 0.2e1 * t306 - qJDD(4) + (-t175 * t62 + t178 * t96) * qJD(2) + t193 (-pkin(3) * t175 + t309) * qJDD(2) (t175 * t96 + t178 * t62) * qJD(2) + t192 + t335, -t21 * pkin(3) - g(1) * t250 - g(2) * t251 - g(3) * t249 + t20 * qJ(4) + t340 * t53 - t52 * t59 - t62 * t96, t129, -t106, t151, -t129, t269, qJDD(3), -qJ(5) * t151 + t48 * qJD(3) + t127 + (-t175 * t74 + t178 * t292) * qJD(2) + t196 + t335, -qJD(3) * t49 + 0.2e1 * t257 + ((-qJ(5) * qJD(3) + t74) * t178 + t243) * qJD(2) + t189 (-t175 * t180 - t309) * qJDD(2) + (-t294 + t31 - t262) * t284, t11 * t180 - t12 * qJ(4) - t31 * t49 - t40 * t74 - g(1) * (t250 - t38) - g(2) * (t251 - t36) - g(3) * (t249 - t77) + t294 * t337, t305 * t95 - t315 (-t29 - t313) * t177 + (-t30 - t312) * t174 (-t265 + t316) * qJD(2) + t216, t220 * t94 - t314 (t266 - t317) * qJD(2) + t215, -t133 * t284, -t18 * t133 - t173 * t30 + t174 * t200 - t177 * t195 + t227 * t284 - t294 * t94, t19 * t133 + t173 * t29 + t174 * t195 + t177 * t200 + t284 * t9 - t294 * t95, -t18 * t95 - t19 * t94 + (-t227 * t286 + t164 * t30 - t1 + (-t164 * t95 - t227) * qJD(6)) * t177 + (t9 * t286 + t164 * t29 + t2 + (-t164 * t94 + t9) * qJD(6)) * t174 + t214, t10 * t173 - t9 * t19 + t227 * t18 - g(1) * (-pkin(9) * t43 + t173 * t44 - t38 - t39) - g(2) * (-pkin(9) * t41 + t173 * t42 - t36 - t37) - g(3) * (-t88 * pkin(9) + t173 * t89 - t77 - t78) + t294 * t32 + (qJD(6) * t230 + t1 * t177 - t2 * t174) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t269, t125, -qJD(3) * t53 + t286 * t62 - t194 - t306, 0, 0, 0, 0, 0, 0, t125, t110, -t269, -t295 + t257 + (-qJ(5) * t280 + t243) * qJD(2) + t189, 0, 0, 0, 0, 0, 0, -t133 ^ 2 * t177 + t308 - t324, t133 * t220 + t307 - t322, t174 * t336 + t177 * t211, -t32 * qJD(3) + t174 * t338 + t177 * t339 - t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t253 + t269, -t151 + 0.2e1 * t254, -t152 - t304 (t320 + (t31 + t262) * t175) * qJD(2) + t199 - t204, 0, 0, 0, 0, 0, 0 (-t266 - t317) * qJD(2) - t215 (-t265 - t316) * qJD(2) + t216, t174 * t211 - t177 * t336 (t175 * t229 + t178 * t32) * qJD(2) + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, -t94 ^ 2 + t95 ^ 2, t336, -t328, t211, t93, t32 * t95 - g(1) * (-t174 * t43 + t177 * t84) - g(2) * (-t174 * t41 + t177 * t82) - g(3) * t50 - t338, -t32 * t94 - g(1) * (-t174 * t84 - t177 * t43) - g(2) * (-t174 * t82 - t177 * t41) + g(3) * t51 - t339, 0, 0;];
tau_reg  = t3;
