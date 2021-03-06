% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:45
% EndTime: 2019-03-09 02:51:56
% DurationCPUTime: 6.86s
% Computational Cost: add. (8950->522), mult. (20929->623), div. (0->0), fcn. (16062->14), ass. (0->268)
t214 = sin(pkin(9));
t221 = sin(qJ(3));
t300 = qJD(3) * t221;
t286 = t214 * t300;
t216 = cos(pkin(9));
t364 = cos(qJ(3));
t290 = t364 * t216;
t272 = qJD(1) * t290;
t281 = qJDD(1) * t364;
t294 = t216 * qJDD(1);
t292 = qJD(3) * t272 + t214 * t281 + t221 * t294;
t111 = qJD(1) * t286 - t292;
t109 = -qJDD(6) + t111;
t167 = t364 * t214 + t221 * t216;
t373 = t167 * qJD(1);
t148 = qJD(6) + t373;
t213 = sin(pkin(10));
t215 = cos(pkin(10));
t220 = sin(qJ(6));
t363 = cos(qJ(6));
t166 = -t220 * t213 + t363 * t215;
t245 = -t363 * t213 - t220 * t215;
t158 = t245 * qJD(6);
t393 = t245 * t373 + t158;
t381 = -t166 * t109 + t148 * t393;
t312 = t214 * t221;
t287 = qJD(1) * t312;
t154 = -t272 + t287;
t125 = qJD(3) * t215 + t154 * t213;
t198 = t216 * pkin(2) + pkin(1);
t173 = -t198 * qJD(1) + qJD(2);
t238 = -qJ(4) * t373 + t173;
t353 = pkin(3) + qJ(5);
t59 = t353 * t154 + t238;
t350 = pkin(7) + qJ(2);
t174 = t350 * t214;
t168 = qJD(1) * t174;
t175 = t350 * t216;
t169 = qJD(1) * t175;
t115 = t364 * t168 + t221 * t169;
t372 = qJD(4) + t115;
t307 = pkin(4) * t373 + t372;
t65 = -t353 * qJD(3) + t307;
t33 = -t213 * t59 + t215 * t65;
t18 = pkin(5) * t373 - pkin(8) * t125 + t33;
t123 = t213 * qJD(3) - t154 * t215;
t34 = t213 * t65 + t215 * t59;
t22 = -pkin(8) * t123 + t34;
t247 = -t363 * t18 + t220 * t22;
t161 = t167 * qJD(3);
t295 = t214 * qJDD(1);
t267 = -t216 * t281 + t221 * t295;
t112 = qJD(1) * t161 + t267;
t172 = -t198 * qJDD(1) + qJDD(2);
t330 = t111 * qJ(4);
t229 = -qJD(4) * t373 + t172 + t330;
t23 = t154 * qJD(5) + t353 * t112 + t229;
t296 = qJD(1) * qJD(2);
t368 = t350 * qJDD(1) + t296;
t140 = t368 * t214;
t141 = t368 * t216;
t284 = qJD(3) * t364;
t277 = t364 * t140 + t221 * t141 - t168 * t300 + t169 * t284;
t258 = qJDD(4) + t277;
t35 = -t111 * pkin(4) - qJD(3) * qJD(5) - t353 * qJDD(3) + t258;
t10 = -t213 * t23 + t215 * t35;
t100 = qJDD(3) * t215 + t112 * t213;
t6 = -pkin(5) * t111 - pkin(8) * t100 + t10;
t11 = t213 * t35 + t215 * t23;
t99 = qJDD(3) * t213 - t112 * t215;
t9 = -pkin(8) * t99 + t11;
t1 = -qJD(6) * t247 + t220 * t6 + t363 * t9;
t8 = t220 * t18 + t363 * t22;
t2 = -qJD(6) * t8 - t220 * t9 + t363 * t6;
t209 = pkin(9) + qJ(3);
t204 = cos(t209);
t195 = g(3) * t204;
t202 = sin(t209);
t223 = cos(qJ(1));
t318 = t202 * t223;
t222 = sin(qJ(1));
t319 = t202 * t222;
t291 = -g(1) * t318 - g(2) * t319 + t195;
t283 = qJD(6) * t363;
t299 = qJD(6) * t220;
t159 = -t213 * t299 + t215 * t283;
t390 = t166 * t373 + t159;
t392 = -t1 * t245 + t2 * t166 - t247 * t393 + t390 * t8 + t291;
t270 = -t245 * t109 - t148 * t390;
t280 = t220 * t100 + t363 * t99;
t70 = -t220 * t123 + t363 * t125;
t27 = qJD(6) * t70 + t280;
t67 = t363 * t123 + t125 * t220;
t391 = -t245 * t27 + t390 * t67;
t26 = -t363 * t100 + t123 * t283 + t125 * t299 + t220 * t99;
t389 = -t166 * t26 + t393 * t70;
t387 = t67 ^ 2;
t386 = t70 ^ 2;
t385 = t148 * t67;
t105 = t215 * t111;
t366 = t373 ^ 2;
t378 = t366 * t213 + t105;
t334 = qJDD(1) * pkin(1);
t377 = g(1) * t222 - g(2) * t223;
t251 = -qJDD(2) + t334 + t377;
t269 = g(1) * t223 + g(2) * t222;
t239 = -g(3) * t202 - t204 * t269;
t278 = t221 * t140 - t364 * t141 + t168 * t284 + t169 * t300;
t228 = t239 - t278;
t383 = -t115 * qJD(3) - t228;
t122 = -t221 * t174 + t364 * t175;
t246 = t377 * t202;
t93 = (qJD(2) * t214 + qJD(3) * t175) * t221 - qJD(2) * t290 + t174 * t284;
t382 = t93 * qJD(3) - t122 * qJDD(3) - t246;
t380 = 0.2e1 * qJD(3) * t373 + t267;
t304 = t204 * pkin(3) + t202 * qJ(4);
t210 = qJDD(3) * qJ(4);
t211 = qJD(3) * qJD(4);
t376 = -t210 - t211;
t367 = t154 ^ 2;
t375 = -t367 - t366;
t374 = -t367 + t366;
t371 = qJ(2) * qJDD(1);
t370 = -t112 * pkin(4) + qJDD(5);
t369 = -t245 * t26 - t390 * t70;
t329 = t111 * t213;
t262 = -t215 * t366 + t329;
t365 = t99 * pkin(5);
t362 = pkin(5) * t213;
t361 = pkin(8) * t215;
t356 = t112 * pkin(3);
t354 = t70 * t67;
t351 = pkin(4) + t350;
t349 = -pkin(8) - t353;
t165 = -t290 + t312;
t160 = -t216 * t284 + t286;
t257 = t160 * qJ(4) - t167 * qJD(4);
t48 = t165 * qJD(5) + t353 * t161 + t257;
t94 = qJD(2) * t167 + qJD(3) * t122;
t61 = -t160 * pkin(4) + t94;
t25 = t213 * t61 + t215 * t48;
t323 = t154 * qJ(4);
t76 = t353 * t373 + t323;
t116 = -t221 * t168 + t364 * t169;
t81 = -t154 * pkin(4) + t116;
t39 = t213 * t81 + t215 * t76;
t121 = t364 * t174 + t175 * t221;
t101 = pkin(4) * t167 + t121;
t252 = -t167 * qJ(4) - t198;
t79 = t353 * t165 + t252;
t42 = t213 * t101 + t215 * t79;
t170 = t349 * t213;
t171 = t349 * t215;
t117 = -t220 * t170 + t363 * t171;
t75 = t215 * t81;
t28 = -t154 * pkin(5) + t75 + (-pkin(8) * t373 - t76) * t213;
t321 = t373 * t215;
t30 = pkin(8) * t321 + t39;
t348 = qJD(5) * t245 + qJD(6) * t117 - t220 * t28 - t363 * t30;
t118 = t363 * t170 + t220 * t171;
t347 = -qJD(5) * t166 - qJD(6) * t118 + t220 * t30 - t363 * t28;
t346 = t154 * t67;
t342 = t70 * t154;
t341 = t99 * t213;
t335 = qJ(5) * t204;
t333 = qJDD(3) * pkin(3);
t332 = t100 * t213;
t331 = t100 * t215;
t328 = t123 * t154;
t327 = t123 * t213;
t326 = t125 * t154;
t325 = t125 * t213;
t324 = t125 * t215;
t322 = t154 * t373;
t208 = pkin(10) + qJ(6);
t203 = cos(t208);
t317 = t203 * t222;
t316 = t203 * t223;
t217 = -pkin(8) - qJ(5);
t315 = t204 * t217;
t314 = t204 * t222;
t313 = t204 * t223;
t309 = t223 * t350;
t197 = pkin(5) * t215 + pkin(4);
t308 = t197 * t373 + t372;
t212 = qJD(3) * qJ(4);
t72 = t212 + qJD(5) + t81;
t306 = qJD(5) - t72;
t303 = t197 + t350;
t206 = t214 ^ 2;
t207 = t216 ^ 2;
t302 = t206 + t207;
t301 = qJD(3) * t116;
t293 = t204 * t362;
t279 = t302 * qJD(1) ^ 2;
t276 = -t166 * t27 - t393 * t67;
t275 = 0.2e1 * t302;
t179 = t223 * t198;
t274 = g(2) * (pkin(3) * t313 + qJ(4) * t318 + t179);
t271 = -g(1) * t314 + g(2) * t313;
t266 = t10 * t215 + t11 * t213;
t265 = -t10 * t213 + t11 * t215;
t264 = -t213 * t34 - t215 * t33;
t263 = -t213 * t33 + t215 * t34;
t53 = -t111 * t167 - t160 * t373;
t261 = -t111 * t165 + t161 * t373;
t260 = t112 * t165 + t154 * t161;
t259 = -t324 + t327;
t254 = qJD(3) * t160 - qJDD(3) * t167;
t253 = qJD(3) * t161 + qJDD(3) * t165;
t249 = t202 * t362 - t315;
t248 = -t198 - t304;
t92 = t215 * t101;
t29 = t167 * pkin(5) + t92 + (-pkin(8) * t165 - t79) * t213;
t37 = t165 * t361 + t42;
t14 = -t220 * t37 + t363 * t29;
t15 = t220 * t29 + t363 * t37;
t43 = t278 + t376;
t244 = -t277 - t291;
t240 = t251 + t334;
t36 = -t43 + t370;
t237 = t161 * t72 + t165 * t36 + t269;
t235 = -qJD(3) * t94 - qJDD(3) * t121 - t271;
t234 = t36 + t239;
t233 = t172 - t377;
t232 = -t112 * t167 + t154 * t160 - t261;
t90 = t154 * pkin(3) + t238;
t231 = t373 * t90 + qJDD(4) - t244;
t230 = t275 * t296 - t269;
t227 = -t111 * t121 - t112 * t122 + t154 * t93 + t373 * t94 - t269;
t226 = t228 + t370 - t376;
t201 = sin(t208);
t193 = qJ(4) + t362;
t178 = qJ(4) * t313;
t176 = qJ(4) * t314;
t143 = qJD(3) * t154;
t138 = -t201 * t319 + t316;
t137 = t201 * t223 + t202 * t317;
t136 = t201 * t318 + t317;
t135 = -t201 * t222 + t202 * t316;
t110 = pkin(3) * t165 + t252;
t108 = pkin(3) * t373 + t323;
t107 = -t212 - t116;
t106 = -qJD(3) * pkin(3) + t372;
t104 = t245 * t165;
t103 = t166 * t165;
t102 = -t165 * pkin(4) + t122;
t86 = t215 * t99;
t78 = pkin(3) * t161 + t257;
t77 = t111 - t143;
t64 = -t197 * t165 + t122;
t60 = -t161 * pkin(4) - t93;
t58 = t215 * t61;
t50 = t123 * pkin(5) + t72;
t49 = -t197 * t161 - t93;
t47 = t258 - t333;
t46 = -t165 * t158 - t161 * t166;
t45 = -t159 * t165 + t161 * t245;
t41 = -t213 * t79 + t92;
t40 = t229 + t356;
t38 = -t213 * t76 + t75;
t24 = -t213 * t48 + t58;
t19 = t36 + t365;
t17 = t161 * t361 + t25;
t16 = -t160 * pkin(5) + t58 + (-pkin(8) * t161 - t48) * t213;
t4 = -t15 * qJD(6) + t363 * t16 - t220 * t17;
t3 = t14 * qJD(6) + t220 * t16 + t363 * t17;
t5 = [0, 0, 0, 0, 0, qJDD(1), t377, t269, 0, 0, t206 * qJDD(1), 0.2e1 * t214 * t294, 0, t207 * qJDD(1), 0, 0, t240 * t216, -t240 * t214, t275 * t371 + t230, pkin(1) * t251 + (t302 * t371 + t230) * qJ(2), t53, t232, -t254, t260, -t253, 0, -t112 * t198 + t161 * t173 + t165 * t172 + t235, t198 * t111 - t173 * t160 + t172 * t167 + t382, -t115 * t160 - t116 * t161 + t165 * t278 + t167 * t277 + t227, -t278 * t122 - t116 * t93 + t277 * t121 + t115 * t94 - t172 * t198 - g(1) * (-t222 * t198 + t309) - g(2) * (t222 * t350 + t179) 0, t254, t253, t53, t232, t260, -t106 * t160 + t107 * t161 + t165 * t43 + t167 * t47 + t227, -t110 * t112 - t154 * t78 - t161 * t90 - t165 * t40 - t235, t110 * t111 + t90 * t160 - t40 * t167 - t373 * t78 - t382, t40 * t110 + t90 * t78 - t43 * t122 + t107 * t93 + t47 * t121 + t106 * t94 - g(1) * t309 - t274 + (-g(1) * t248 - g(2) * t350) * t222 (t100 * t165 + t125 * t161) * t213 (t331 - t341) * t165 - t259 * t161, t100 * t167 - t125 * t160 + t213 * t261 (-t123 * t161 - t165 * t99) * t215, t123 * t160 - t99 * t167 + t215 * t261, t53, t10 * t167 + t102 * t99 - t41 * t111 + t60 * t123 - t33 * t160 + t213 * t246 - t215 * t237 + t24 * t373, t102 * t100 - t11 * t167 + t42 * t111 + t60 * t125 + t34 * t160 + t213 * t237 + t215 * t246 - t25 * t373, -t41 * t100 - t25 * t123 - t24 * t125 + t161 * t263 + t165 * t265 - t42 * t99 - t271, t11 * t42 + t34 * t25 + t10 * t41 + t33 * t24 + t36 * t102 + t72 * t60 - t274 + (-g(1) * t351 - g(2) * t335) * t223 + (-g(1) * (t248 - t335) - g(2) * t351) * t222, t104 * t26 - t45 * t70, -t103 * t26 + t104 * t27 + t45 * t67 - t46 * t70, t104 * t109 - t148 * t45 - t160 * t70 - t167 * t26, -t103 * t27 + t46 * t67, -t103 * t109 - t148 * t46 + t160 * t67 - t167 * t27, -t109 * t167 - t148 * t160, -g(1) * t138 - g(2) * t136 - t103 * t19 - t109 * t14 + t148 * t4 + t160 * t247 + t167 * t2 + t27 * t64 + t46 * t50 + t49 * t67, g(1) * t137 - g(2) * t135 - t1 * t167 - t104 * t19 + t109 * t15 - t148 * t3 + t160 * t8 - t26 * t64 - t45 * t50 + t49 * t70, t1 * t103 + t104 * t2 + t14 * t26 - t15 * t27 - t247 * t45 - t3 * t67 - t4 * t70 - t46 * t8 - t271, t1 * t15 + t8 * t3 + t2 * t14 - t247 * t4 + t19 * t64 + t50 * t49 - t274 + (-g(1) * t303 - g(2) * t249) * t223 + (-g(1) * (t248 - t249) - g(2) * t303) * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t294, t295, -t279, -qJ(2) * t279 - t251, 0, 0, 0, 0, 0, 0, t380 (-t154 - t287) * qJD(3) + t292, t375, -t115 * t373 + t116 * t154 + t233, 0, 0, 0, 0, 0, 0, t375, -t380, t111 + t143, t356 + t330 - t107 * t154 + (-qJD(4) - t106) * t373 + t233, 0, 0, 0, 0, 0, 0, t262 + t328, t326 + t378, t332 - t86 + (t324 + t327) * t373, t72 * t154 + t264 * t373 + t265 - t377, 0, 0, 0, 0, 0, 0, t270 + t346, t342 - t381, t276 - t369, t1 * t166 + t50 * t154 + t2 * t245 + t247 * t390 + t393 * t8 - t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, t374 (t154 - t287) * qJD(3) + t292, -t322, -t267, qJDD(3), -t173 * t373 + t244 + t301, t173 * t154 + t383, 0, 0, qJDD(3), t77, t267, t322, t374, -t322, pkin(3) * t111 - qJ(4) * t112 + (-t107 - t116) * t373 + (t106 - t372) * t154, t108 * t154 + t231 - t301 - 0.2e1 * t333, t108 * t373 - t90 * t154 + 0.2e1 * t210 + 0.2e1 * t211 - t383, -t43 * qJ(4) - t47 * pkin(3) - t90 * t108 - t106 * t116 - g(1) * (-pkin(3) * t318 + t178) - g(2) * (-pkin(3) * t319 + t176) - g(3) * t304 - t372 * t107, -t325 * t373 + t331, t259 * t373 - t332 - t86, t326 - t378, t123 * t321 + t341, t262 - t328, t322, t353 * t105 + qJ(4) * t99 + t33 * t154 + t307 * t123 + (-t215 * t306 - t38) * t373 + t234 * t213, -t353 * t329 + qJ(4) * t100 - t34 * t154 + t307 * t125 + (t213 * t306 + t39) * t373 + t234 * t215, t39 * t123 + t38 * t125 + (qJD(5) * t125 + t100 * t353 - t34 * t373 - t10) * t215 + (qJD(5) * t123 + t33 * t373 + t353 * t99 - t11) * t213 - t291, t36 * qJ(4) - t34 * t39 - t33 * t38 - g(1) * t178 - g(2) * t176 - g(3) * (t304 + t335) + t307 * t72 + t264 * qJD(5) + (t202 * t269 - t266) * t353, t389, t276 + t369, t342 + t381, t391, t270 - t346, t148 * t154, -t117 * t109 + t148 * t347 - t154 * t247 - t19 * t245 + t193 * t27 + t201 * t239 + t308 * t67 + t390 * t50, t118 * t109 - t148 * t348 - t8 * t154 + t19 * t166 - t193 * t26 + t203 * t239 + t308 * t70 + t393 * t50, t117 * t26 - t118 * t27 - t347 * t70 - t348 * t67 - t392, t1 * t118 + t2 * t117 + t19 * t193 - g(1) * (t223 * t293 + t178) - g(2) * (t222 * t293 + t176) - g(3) * (t304 - t315) + t348 * t8 - t347 * t247 + t308 * t50 + (-g(3) * t362 + t269 * (pkin(3) - t217)) * t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, qJDD(3) - t322, -qJD(3) ^ 2 - t366, qJD(3) * t107 + t231 - t333, 0, 0, 0, 0, 0, 0, -qJD(3) * t123 - t378, -qJD(3) * t125 + t262, -t331 - t341 + (-t123 * t215 + t325) * t373, -t72 * qJD(3) + t263 * t373 + t266 + t291, 0, 0, 0, 0, 0, 0, -qJD(3) * t67 + t381, -qJD(3) * t70 + t270, -t389 - t391, -t50 * qJD(3) + t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125 * t373 + t99, -t123 * t373 + t100, -t123 ^ 2 - t125 ^ 2, t34 * t123 + t33 * t125 + t226, 0, 0, 0, 0, 0, 0, t70 * t148 + t27, -t26 - t385, -t386 - t387, -t247 * t70 + t67 * t8 + t226 + t365; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, t386 - t387, -t26 + t385, -t354, -t280 + (-qJD(6) + t148) * t70, -t109, -g(1) * t135 - g(2) * t137 + t8 * t148 + t195 * t203 - t50 * t70 + t2, g(1) * t136 - g(2) * t138 - t148 * t247 - t195 * t201 + t50 * t67 - t1, 0, 0;];
tau_reg  = t5;
