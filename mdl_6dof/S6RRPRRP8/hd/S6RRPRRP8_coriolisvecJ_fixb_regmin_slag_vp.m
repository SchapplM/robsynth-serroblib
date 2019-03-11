% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:16
% EndTime: 2019-03-09 12:25:31
% DurationCPUTime: 5.89s
% Computational Cost: add. (10087->461), mult. (24923->610), div. (0->0), fcn. (18497->8), ass. (0->230)
t220 = sin(pkin(10));
t223 = sin(qJ(4));
t221 = cos(pkin(10));
t225 = cos(qJ(4));
t309 = t221 * t225;
t180 = t220 * t223 - t309;
t226 = cos(qJ(2));
t241 = t180 * t226;
t360 = -qJD(1) * t241 + qJD(4) * t180;
t224 = sin(qJ(2));
t297 = qJD(1) * t224;
t281 = t220 * t297;
t288 = t221 * qJD(2);
t173 = t281 - t288;
t280 = t221 * t297;
t295 = qJD(2) * t220;
t175 = t280 + t295;
t117 = t173 * t223 - t175 * t225;
t222 = sin(qJ(5));
t257 = -t173 * t225 - t175 * t223;
t331 = cos(qJ(5));
t67 = t117 * t222 + t257 * t331;
t362 = t67 ^ 2;
t296 = qJD(1) * t226;
t208 = -qJD(4) + t296;
t200 = -qJD(5) + t208;
t361 = t200 * t67;
t311 = t221 * t223;
t181 = t220 * t225 + t311;
t340 = t226 * t181;
t142 = qJD(1) * t340;
t238 = t181 * qJD(4);
t337 = t142 - t238;
t259 = pkin(2) * t224 - qJ(3) * t226;
t183 = t259 * qJD(1);
t133 = pkin(7) * t281 + t183 * t221;
t308 = t221 * t226;
t255 = pkin(3) * t224 - pkin(8) * t308;
t107 = qJD(1) * t255 + t133;
t164 = t220 * t183;
t310 = t221 * t224;
t312 = t220 * t226;
t244 = -pkin(7) * t310 - pkin(8) * t312;
t123 = qJD(1) * t244 + t164;
t327 = pkin(8) + qJ(3);
t192 = t327 * t220;
t193 = t327 * t221;
t299 = -t192 * t223 + t193 * t225;
t359 = qJD(3) * t181 + qJD(4) * t299 + t107 * t225 - t123 * t223;
t290 = qJD(4) * t225;
t358 = qJD(3) * t309 - t107 * t223 - t123 * t225 - t192 * t290;
t349 = -t117 * t331 + t222 * t257;
t332 = t349 ^ 2;
t357 = pkin(4) * t297 - pkin(9) * t360 + t359;
t292 = qJD(3) * t220;
t306 = t223 * t193;
t330 = pkin(9) * t181;
t356 = pkin(9) * t142 - t223 * t292 + (-t306 - t330) * qJD(4) + t358;
t355 = t200 * t349;
t343 = t180 * t224;
t354 = pkin(9) * t343;
t213 = pkin(7) * t297;
t320 = qJD(2) * pkin(2);
t271 = qJD(3) - t320;
t189 = t213 + t271;
t132 = pkin(3) * t173 + t189;
t77 = -pkin(4) * t257 + t132;
t24 = -pkin(5) * t67 - qJ(6) * t349 + t77;
t353 = t24 * t67;
t352 = t67 * t77;
t279 = qJD(5) * t331;
t289 = qJD(5) * t222;
t322 = t180 * t279 + t181 * t289 - t222 * t337 + t331 * t360;
t122 = -t180 * t222 + t181 * t331;
t321 = t122 * qJD(5) - t222 * t360 - t331 * t337;
t328 = t349 * t67;
t351 = t117 * t208;
t350 = t332 - t362;
t287 = qJD(1) * qJD(2);
t278 = t226 * t287;
t262 = t225 * t278;
t291 = qJD(4) * t223;
t264 = -t173 * t291 + t175 * t290 + t220 * t262 + t278 * t311;
t263 = t220 * t278;
t80 = -t173 * t290 + t221 * t262 + (-qJD(4) * t175 - t263) * t223;
t29 = -t117 * t289 + t222 * t264 - t257 * t279 - t331 * t80;
t18 = -t29 + t361;
t37 = pkin(5) * t349 - qJ(6) * t67;
t348 = -0.2e1 * t287;
t329 = t24 * t349;
t346 = t77 * t349;
t268 = -t192 * t225 - t306;
t101 = t268 - t330;
t102 = -pkin(9) * t180 + t299;
t246 = t101 * t331 - t102 * t222;
t345 = t246 * qJD(5) - t222 * t357 + t331 * t356;
t59 = t101 * t222 + t102 * t331;
t344 = t59 * qJD(5) + t222 * t356 + t331 * t357;
t342 = t208 * t257;
t214 = pkin(7) * t296;
t285 = pkin(3) * t296;
t167 = t220 * t285 + t214;
t338 = -pkin(4) * t337 - t167;
t190 = -pkin(2) * t226 - qJ(3) * t224 - pkin(1);
t172 = t221 * t190;
t124 = -pkin(8) * t310 + t172 + (-pkin(7) * t220 - pkin(3)) * t226;
t206 = pkin(7) * t308;
t140 = t190 * t220 + t206;
t313 = t220 * t224;
t131 = -pkin(8) * t313 + t140;
t302 = t124 * t223 + t131 * t225;
t30 = qJD(5) * t349 + t222 * t80 + t264 * t331;
t336 = -t30 - t355;
t307 = t223 * t131;
t269 = t124 * t225 - t307;
t55 = -pkin(4) * t226 + t269 + t354;
t153 = t181 * t224;
t60 = -pkin(9) * t153 + t302;
t249 = t222 * t55 + t331 * t60;
t103 = -qJD(2) * t241 - t224 * t238;
t160 = qJD(2) * t259 - qJD(3) * t224;
t146 = t220 * t160;
t109 = qJD(2) * t244 + t146;
t294 = qJD(2) * t224;
t284 = pkin(7) * t294;
t129 = t160 * t221 + t220 * t284;
t237 = t255 * qJD(2);
t98 = t237 + t129;
t273 = -t109 * t223 + t225 * t98;
t26 = pkin(4) * t294 - pkin(9) * t103 - qJD(4) * t302 + t273;
t235 = qJD(2) * t340;
t283 = t109 * t225 + t124 * t290 + t223 * t98;
t31 = -pkin(9) * t235 + (-t307 + t354) * qJD(4) + t283;
t333 = -qJD(5) * t249 - t222 * t31 + t26 * t331;
t326 = -qJ(6) * t297 + t345;
t325 = pkin(5) * t297 + t344;
t324 = pkin(5) * t321 + qJ(6) * t322 - t122 * qJD(6) + t338;
t166 = t190 * qJD(1);
t195 = qJD(2) * qJ(3) + t214;
t125 = t166 * t221 - t195 * t220;
t86 = -pkin(8) * t175 + t125 - t285;
t126 = t166 * t220 + t195 * t221;
t90 = -pkin(8) * t173 + t126;
t50 = t223 * t86 + t225 * t90;
t42 = pkin(9) * t257 + t50;
t319 = t222 * t42;
t49 = -t223 * t90 + t225 * t86;
t41 = pkin(9) * t117 + t49;
t38 = -pkin(4) * t208 + t41;
t9 = t331 * t38 - t319;
t318 = qJD(6) - t9;
t14 = t331 * t41 - t319;
t316 = -pkin(4) * t279 - qJD(6) + t14;
t315 = qJD(2) * t246;
t314 = qJD(2) * t59;
t228 = qJD(1) ^ 2;
t305 = t226 * t228;
t227 = qJD(2) ^ 2;
t304 = t227 * t224;
t303 = t227 * t226;
t144 = t160 * qJD(1);
t188 = (qJD(3) - t213) * qJD(2);
t106 = t144 * t220 + t188 * t221;
t207 = pkin(7) * t278;
t159 = pkin(3) * t263 + t207;
t293 = qJD(2) * t226;
t215 = pkin(7) * t293;
t168 = pkin(3) * t220 * t293 + t215;
t184 = pkin(3) * t313 + pkin(7) * t224;
t298 = t224 ^ 2 - t226 ^ 2;
t286 = pkin(7) * t312;
t282 = t331 * t42;
t210 = -pkin(3) * t221 - pkin(2);
t211 = t224 * t287;
t105 = t144 * t221 - t188 * t220;
t83 = qJD(1) * t237 + t105;
t87 = -pkin(8) * t263 + t106;
t276 = -t223 * t87 + t225 * t83;
t232 = -qJD(4) * t50 + t276;
t15 = pkin(4) * t211 - pkin(9) * t80 + t232;
t248 = t223 * t83 + t225 * t87 + t290 * t86 - t291 * t90;
t19 = -pkin(9) * t264 + t248;
t275 = -t15 * t222 - t19 * t331 - t279 * t38 + t289 * t42;
t274 = t15 * t331 - t19 * t222 - t279 * t42 - t289 * t38;
t272 = pkin(1) * t348;
t267 = t173 + t288;
t266 = -t175 + t295;
t265 = pkin(5) * t211;
t13 = t222 * t41 + t282;
t261 = pkin(4) * t289 - t13;
t128 = pkin(4) * t153 + t184;
t260 = -t189 + t271;
t145 = pkin(4) * t180 + t210;
t191 = t200 * qJD(6);
t204 = qJ(6) * t211;
t1 = t204 - t191 - t275;
t254 = -t200 * t9 + t275;
t10 = t222 * t38 + t282;
t253 = -t10 * t200 + t274;
t251 = -t222 * t60 + t331 * t55;
t247 = t222 * t26 + t279 * t55 - t289 * t60 + t31 * t331;
t97 = -t153 * t222 - t331 * t343;
t234 = qJD(4) * t343;
t61 = pkin(4) * t264 + t159;
t2 = -t265 - t274;
t233 = t29 + t361;
t230 = t30 - t355;
t5 = pkin(5) * t30 + qJ(6) * t29 - qJD(6) * t349 + t61;
t229 = t234 - t235;
t81 = -pkin(4) * t229 + t168;
t212 = -pkin(4) * t331 - pkin(5);
t209 = pkin(4) * t222 + qJ(6);
t139 = t172 - t286;
t134 = -pkin(7) * t280 + t164;
t130 = -t221 * t284 + t146;
t121 = t180 * t331 + t181 * t222;
t96 = t153 * t331 - t222 * t343;
t57 = pkin(5) * t121 - qJ(6) * t122 + t145;
t46 = pkin(5) * t96 - qJ(6) * t97 + t128;
t45 = qJD(5) * t97 + t103 * t222 - t229 * t331;
t44 = -t103 * t331 + t153 * t279 - t222 * t229 - t289 * t343;
t34 = -pkin(4) * t117 + t37;
t33 = pkin(5) * t226 - t251;
t32 = -qJ(6) * t226 + t249;
t8 = -qJ(6) * t200 + t10;
t7 = pkin(5) * t200 + t318;
t6 = t45 * pkin(5) + t44 * qJ(6) - t97 * qJD(6) + t81;
t4 = -pkin(5) * t294 - t333;
t3 = qJ(6) * t294 - qJD(6) * t226 + t247;
t11 = [0, 0, 0, 0.2e1 * t226 * t211, t298 * t348, t303, -t304, 0, -pkin(7) * t303 + t224 * t272, pkin(7) * t304 + t226 * t272 (-qJD(1) * t129 - t105) * t226 + ((pkin(7) * t173 + t189 * t220) * t226 + (t125 + (t139 + 0.2e1 * t286) * qJD(1)) * t224) * qJD(2) (qJD(1) * t130 + t106) * t226 + ((pkin(7) * t175 + t189 * t221) * t226 + (-t126 + (-t140 + 0.2e1 * t206) * qJD(1)) * t224) * qJD(2), -t129 * t175 - t130 * t173 + (-t105 * t221 - t106 * t220) * t224 + (-t125 * t221 - t126 * t220 + (-t139 * t221 - t140 * t220) * qJD(1)) * t293, t105 * t139 + t106 * t140 + t125 * t129 + t126 * t130 + (t189 + t213) * t215, -t103 * t117 - t343 * t80, t103 * t257 - t117 * t229 - t80 * t153 + t264 * t343, -t103 * t208 - t226 * t80 + (-qJD(1) * t343 - t117) * t294, t264 * t226 - t208 * t234 + (t208 * t340 + (-qJD(1) * t153 + t257) * t224) * qJD(2) (-t208 - t296) * t294, -t273 * t208 - t276 * t226 - t168 * t257 + t184 * t264 + t159 * t153 + (-t132 * t343 + t208 * t302 + t226 * t50) * qJD(4) + (t132 * t340 + (qJD(1) * t269 + t49) * t224) * qJD(2) (-t131 * t291 + t283) * t208 + t248 * t226 - t168 * t117 + t184 * t80 - t159 * t343 + t132 * t103 + (-qJD(1) * t302 - t50) * t294, -t29 * t97 - t349 * t44, t29 * t96 - t30 * t97 - t349 * t45 - t44 * t67, t200 * t44 + t226 * t29 + (qJD(1) * t97 + t349) * t294, t200 * t45 + t226 * t30 + (-qJD(1) * t96 + t67) * t294 (-t200 - t296) * t294, t128 * t30 - t200 * t333 + t211 * t251 - t226 * t274 + t294 * t9 + t77 * t45 + t61 * t96 - t67 * t81, t247 * t200 - t275 * t226 + t81 * t349 - t128 * t29 + t61 * t97 - t77 * t44 + (-qJD(1) * t249 - t10) * t294, t2 * t226 + t200 * t4 + t24 * t45 + t30 * t46 + t5 * t96 - t6 * t67 + (-qJD(1) * t33 - t7) * t294, -t1 * t96 + t2 * t97 - t29 * t33 + t3 * t67 - t30 * t32 + t349 * t4 - t44 * t7 - t45 * t8, -t1 * t226 - t200 * t3 + t24 * t44 + t29 * t46 - t5 * t97 - t6 * t349 + (qJD(1) * t32 + t8) * t294, t1 * t32 + t2 * t33 + t24 * t6 + t3 * t8 + t4 * t7 + t46 * t5; 0, 0, 0, -t224 * t305, t298 * t228, 0, 0, 0, t228 * pkin(1) * t224, pkin(1) * t305 ((-qJ(3) * t295 - t125) * t224 + (-pkin(7) * t267 + t220 * t260 + t133) * t226) * qJD(1) ((-qJ(3) * t288 + t126) * t224 + (pkin(7) * t266 + t221 * t260 - t134) * t226) * qJD(1), t133 * t175 + t134 * t173 + (-qJD(3) * t173 + t125 * t296 + t106) * t221 + (qJD(3) * t175 + t126 * t296 - t105) * t220, -t125 * t133 - t126 * t134 + (-t125 * t220 + t126 * t221) * qJD(3) + (-t105 * t220 + t106 * t221) * qJ(3) + (-t189 - t320) * t214, t117 * t360 + t181 * t80, -t117 * t337 - t80 * t180 - t181 * t264 - t257 * t360, t360 * t208 + (qJD(2) * t181 + t117) * t297, -t337 * t208 + (-qJD(2) * t180 - t257) * t297, t208 * t297, t210 * t264 + t159 * t180 + t167 * t257 + (qJD(2) * t268 - t49) * t297 - t337 * t132 + t359 * t208, t167 * t117 + t159 * t181 + t210 * t80 + ((-qJD(4) * t193 - t292) * t223 + t358) * t208 - t360 * t132 + (-qJD(2) * t299 + t50) * t297, -t122 * t29 - t322 * t349, t121 * t29 - t122 * t30 - t321 * t349 - t322 * t67, t322 * t200 + (qJD(2) * t122 - t349) * t297, t321 * t200 + (-qJD(2) * t121 - t67) * t297, t200 * t297, t61 * t121 + t145 * t30 + t321 * t77 + t344 * t200 + (-t9 + t315) * t297 - t338 * t67, t61 * t122 - t145 * t29 - t322 * t77 + t345 * t200 + (t10 - t314) * t297 + t338 * t349, t121 * t5 + t30 * t57 - t324 * t67 + t321 * t24 + t325 * t200 + (t7 + t315) * t297, -t1 * t121 + t122 * t2 + t246 * t29 - t30 * t59 - t321 * t8 - t322 * t7 + t325 * t349 + t326 * t67, -t122 * t5 + t29 * t57 - t324 * t349 + t322 * t24 - t326 * t200 + (-t8 + t314) * t297, t1 * t59 - t2 * t246 + t24 * t324 + t325 * t7 + t326 * t8 + t5 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266 * t296, t267 * t296, -t173 ^ 2 - t175 ^ 2, t125 * t175 + t126 * t173 + t207, 0, 0, 0, 0, 0, t264 + t351, t80 - t342, 0, 0, 0, 0, 0, t230, -t233, t230, -t332 - t362, t233, -t349 * t7 - t67 * t8 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117 * t257, t117 ^ 2 - t257 ^ 2, t80 + t342, -t264 + t351, t211, t117 * t132 - t208 * t50 + t232, -t132 * t257 - t208 * t49 - t248, -t328, t350, t18, t336, t211, -t13 * t200 - t346 + (-t117 * t67 + t200 * t289 + t211 * t331) * pkin(4) + t274, -t14 * t200 - t352 + (t117 * t349 + t200 * t279 - t211 * t222) * pkin(4) + t275, -t329 + t34 * t67 + t261 * t200 + (pkin(5) - t212) * t211 + t274, -t209 * t30 - t212 * t29 + (t261 + t8) * t349 - (t316 + t7) * t67, t200 * t316 + t209 * t211 + t34 * t349 + t1 + t353, t1 * t209 + t2 * t212 - t24 * t34 + t261 * t7 - t316 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328, t350, t18, t336, t211, t253 - t346, t254 - t352, t37 * t67 + t253 + 0.2e1 * t265 - t329, pkin(5) * t29 - qJ(6) * t30 + (-t10 + t8) * t349 - (t7 - t318) * t67, t349 * t37 - 0.2e1 * t191 + 0.2e1 * t204 - t254 + t353, -pkin(5) * t2 + qJ(6) * t1 - t10 * t7 - t24 * t37 + t318 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211 - t328, t18, -t200 ^ 2 - t332, t200 * t8 + t2 + t329;];
tauc_reg  = t11;
