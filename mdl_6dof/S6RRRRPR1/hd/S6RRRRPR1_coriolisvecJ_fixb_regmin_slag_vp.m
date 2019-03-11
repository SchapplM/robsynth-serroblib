% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x33]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:54:55
% EndTime: 2019-03-09 21:55:10
% DurationCPUTime: 5.42s
% Computational Cost: add. (12671->398), mult. (33265->532), div. (0->0), fcn. (25733->10), ass. (0->247)
t213 = cos(qJ(6));
t284 = qJD(6) * t213;
t215 = cos(qJ(3));
t216 = cos(qJ(2));
t289 = qJD(1) * t216;
t277 = t215 * t289;
t211 = sin(qJ(3));
t212 = sin(qJ(2));
t290 = qJD(1) * t212;
t278 = t211 * t290;
t164 = -t277 + t278;
t166 = -t211 * t289 - t215 * t290;
t210 = sin(qJ(4));
t214 = cos(qJ(4));
t130 = -t214 * t164 + t166 * t210;
t208 = sin(pkin(11));
t242 = t164 * t210 + t214 * t166;
t310 = cos(pkin(11));
t342 = t310 * t130 + t208 * t242;
t352 = t213 * t342;
t355 = t284 - t352;
t201 = -pkin(2) * t216 - pkin(1);
t185 = t201 * qJD(1);
t143 = t164 * pkin(3) + t185;
t205 = qJD(2) + qJD(3);
t283 = qJD(1) * qJD(2);
t274 = t216 * t283;
t136 = qJD(3) * t277 - t205 * t278 + t215 * t274;
t335 = pkin(7) + pkin(8);
t187 = t335 * t216;
t181 = qJD(1) * t187;
t171 = t215 * t181;
t186 = t335 * t212;
t179 = qJD(1) * t186;
t320 = qJD(2) * pkin(2);
t173 = -t179 + t320;
t241 = -t173 * t211 - t171;
t279 = qJD(2) * t335;
t252 = qJD(1) * t279;
t174 = t212 * t252;
t175 = t216 * t252;
t257 = t211 * t174 - t215 * t175;
t225 = qJD(3) * t241 + t257;
t69 = -pkin(9) * t136 + t225;
t332 = pkin(9) * t164;
t108 = -t241 - t332;
t287 = qJD(4) * t210;
t98 = t108 * t287;
t354 = -t143 * t130 - t210 * t69 + t98;
t209 = sin(qJ(6));
t204 = qJD(4) + t205;
t285 = qJD(6) * t209;
t341 = t130 * t208 - t242 * t310;
t178 = t211 * t216 + t212 * t215;
t142 = t205 * t178;
t137 = t142 * qJD(1);
t243 = t210 * t136 + t214 * t137;
t222 = -qJD(4) * t242 + t243;
t286 = qJD(4) * t214;
t67 = t214 * t136 - t210 * t137 - t164 * t286 + t166 * t287;
t42 = -t208 * t222 + t310 * t67;
t22 = t204 * t284 + t213 * t42 - t285 * t341;
t77 = t204 * t209 + t213 * t341;
t23 = qJD(6) * t77 + t209 * t42;
t75 = -t213 * t204 + t209 * t341;
t353 = -t209 * t23 + t22 * t213 - t355 * t75;
t20 = t22 * t209;
t11 = t355 * t77 + t20;
t41 = t208 * t67 + t310 * t222;
t38 = t209 * t41;
t84 = qJD(6) - t342;
t83 = t84 * t284;
t324 = t38 + t83;
t327 = t77 * t341;
t8 = -t352 * t84 + t324 - t327;
t316 = t209 * t84;
t328 = t75 * t341;
t326 = t84 * t341;
t351 = qJD(6) - t84;
t255 = t179 * t211 - t171;
t110 = t255 + t332;
t160 = t166 * pkin(9);
t167 = t211 * t181;
t292 = -t215 * t179 - t167;
t111 = t160 + t292;
t200 = pkin(2) * t215 + pkin(3);
t300 = t210 * t211;
t350 = -t200 * t286 - (-t211 * t287 + (t214 * t215 - t300) * qJD(3)) * pkin(2) + t210 * t110 + t214 * t111;
t299 = t211 * t214;
t349 = -t200 * t287 + (-t211 * t286 + (-t210 * t215 - t299) * qJD(3)) * pkin(2) - t214 * t110 + t111 * t210;
t103 = t214 * t108;
t258 = t215 * t173 - t167;
t107 = t160 + t258;
t97 = pkin(3) * t205 + t107;
t245 = -t210 * t97 - t103;
t288 = qJD(3) * t211;
t256 = -t211 * t175 - t181 * t288;
t338 = (qJD(3) * t173 - t174) * t215;
t68 = -pkin(9) * t137 + t256 + t338;
t246 = -t210 * t68 + t214 * t69;
t223 = qJD(4) * t245 + t246;
t306 = t143 * t242;
t348 = t223 + t306;
t267 = -qJD(4) * t97 - t68;
t347 = t214 * t267 + t354;
t121 = t242 * qJ(5);
t101 = t210 * t108;
t266 = t214 * t97 - t101;
t57 = t266 + t121;
t54 = pkin(4) * t204 + t57;
t309 = qJ(5) * t130;
t58 = -t245 + t309;
t55 = t310 * t58;
t29 = t208 * t54 + t55;
t27 = pkin(10) * t204 + t29;
t95 = -pkin(4) * t130 + qJD(5) + t143;
t43 = -pkin(5) * t342 - pkin(10) * t341 + t95;
t13 = t209 * t43 + t213 * t27;
t319 = t208 * t58;
t28 = t310 * t54 - t319;
t26 = -t204 * pkin(5) - t28;
t16 = t130 * qJD(5) - t98 + (t69 + (qJD(4) * t164 - t136) * qJ(5)) * t210 + ((qJD(4) * t166 - t137) * qJ(5) - t267) * t214;
t220 = -qJ(5) * t67 + qJD(5) * t242 + t223;
t3 = t16 * t208 - t310 * t220;
t253 = t13 * t341 + t3 * t209 + t26 * t284;
t247 = t209 * t27 - t213 * t43;
t275 = t247 * t341 + t26 * t285;
t346 = pkin(5) * t341 - pkin(10) * t342;
t345 = pkin(4) * t242;
t344 = -t121 - t350;
t307 = t242 * t130;
t343 = t309 + t349;
t59 = -t130 ^ 2 + t242 ^ 2;
t52 = -t130 * t204 + t67;
t53 = -t204 * t242 - t222;
t340 = -0.2e1 * t283;
t39 = t213 * t41;
t339 = t84 * t285 - t39;
t337 = qJD(1) * t178;
t302 = t186 * t215;
t122 = -pkin(9) * t178 - t187 * t211 - t302;
t177 = t211 * t212 - t215 * t216;
t240 = t186 * t211 - t187 * t215;
t123 = -pkin(9) * t177 - t240;
t140 = -t177 * t210 + t178 * t214;
t229 = -qJ(5) * t140 + t122 * t214 - t123 * t210;
t139 = t214 * t177 + t178 * t210;
t244 = -t122 * t210 - t123 * t214;
t66 = -qJ(5) * t139 - t244;
t37 = t208 * t229 + t310 * t66;
t94 = -t208 * t139 + t310 * t140;
t251 = t3 * t94 - t37 * t41;
t4 = t310 * t16 + t208 * t220;
t141 = t205 * t177;
t71 = -qJD(4) * t139 - t214 * t141 - t210 * t142;
t72 = qJD(4) * t140 - t210 * t141 + t214 * t142;
t47 = -t208 * t72 + t310 * t71;
t148 = pkin(3) * t177 + t201;
t230 = pkin(4) * t139 + t148;
t93 = t310 * t139 + t140 * t208;
t49 = pkin(5) * t93 - pkin(10) * t94 + t230;
t180 = t212 * t279;
t182 = t216 * t279;
t231 = -qJD(3) * t302 - t215 * t180 - t211 * t182 - t187 * t288;
t81 = -pkin(9) * t142 + t231;
t224 = qJD(3) * t240 + t211 * t180 - t215 * t182;
t82 = pkin(9) * t141 + t224;
t234 = t122 * t286 - t123 * t287 + t210 * t82 + t214 * t81;
t18 = -qJ(5) * t72 - qJD(5) * t139 + t234;
t269 = -t210 * t81 + t214 * t82;
t221 = -qJ(5) * t71 + qJD(4) * t244 - qJD(5) * t140 + t269;
t6 = t310 * t18 + t208 * t221;
t336 = t26 * t47 - (qJD(6) * t49 + t6) * t84 - (qJD(6) * t43 + t4) * t93 + t251;
t333 = pkin(3) * t166;
t331 = t26 * t342;
t330 = t26 * t94;
t329 = t49 * t41;
t323 = t208 * t344 - t310 * t343;
t322 = t208 * t343 + t310 * t344;
t321 = pkin(3) * qJD(4);
t317 = t209 * t77;
t315 = t209 * t342;
t261 = -t107 * t210 - t103;
t237 = t261 - t309;
t268 = t310 * t210;
t294 = t214 * t107 - t101;
t61 = t121 + t294;
t312 = -t208 * t61 + t310 * t237 + (t208 * t214 + t268) * t321;
t301 = t208 * t210;
t311 = -t208 * t237 - t310 * t61 + (t310 * t214 - t301) * t321;
t304 = t166 * t164;
t303 = t185 * t166;
t218 = qJD(1) ^ 2;
t298 = t216 * t218;
t217 = qJD(2) ^ 2;
t297 = t217 * t212;
t296 = t217 * t216;
t157 = -pkin(2) * t300 + t200 * t214 + pkin(4);
t161 = pkin(2) * t299 + t200 * t210;
t117 = t208 * t157 + t310 * t161;
t199 = pkin(3) * t214 + pkin(4);
t159 = pkin(3) * t268 + t208 * t199;
t291 = t212 ^ 2 - t216 ^ 2;
t203 = t212 * t320;
t202 = pkin(2) * t290;
t280 = t94 * t285;
t276 = -pkin(3) * t204 - t97;
t273 = -pkin(2) * t205 - t173;
t112 = t137 * pkin(3) + qJD(2) * t202;
t132 = pkin(3) * t142 + t203;
t264 = pkin(1) * t340;
t115 = pkin(10) + t117;
t250 = -t333 - t345;
t45 = t250 + t346;
t263 = qJD(6) * t115 + t202 + t45;
t152 = pkin(10) + t159;
t262 = qJD(6) * t152 + t45;
t249 = t28 * t342 + t29 * t341;
t248 = t41 * t94 + t47 * t84;
t239 = t84 * t315 - t339;
t238 = pkin(4) * t72 + t132;
t235 = t185 * t164 - t256;
t116 = t310 * t157 - t208 * t161;
t158 = -pkin(3) * t301 + t310 * t199;
t228 = -t115 * t41 - t322 * t84 - t331;
t227 = -t152 * t41 - t311 * t84 - t331;
t219 = pkin(4) * t222 + t112;
t197 = -t310 * pkin(4) - pkin(5);
t196 = pkin(4) * t208 + pkin(10);
t151 = -pkin(5) - t158;
t144 = t202 - t333;
t114 = -pkin(5) - t116;
t109 = -t164 ^ 2 + t166 ^ 2;
t100 = (-t166 - t337) * t205;
t99 = t164 * t205 + t136;
t48 = -t345 + t346;
t46 = t208 * t71 + t310 * t72;
t36 = t208 * t66 - t310 * t229;
t31 = t310 * t57 - t319;
t30 = t208 * t57 + t55;
t14 = pkin(5) * t46 - pkin(10) * t47 + t238;
t10 = t41 * pkin(5) - t42 * pkin(10) + t219;
t9 = t213 * t10;
t7 = -t316 * t84 + t328 + t39;
t5 = t18 * t208 - t310 * t221;
t1 = -t316 * t77 + t353;
t2 = [0, 0, 0, 0.2e1 * t212 * t274, t291 * t340, t296, -t297, 0, -pkin(7) * t296 + t212 * t264, pkin(7) * t297 + t216 * t264, t136 * t178 + t141 * t166, -t136 * t177 - t137 * t178 + t141 * t164 + t142 * t166, -t141 * t205, -t142 * t205, 0, t201 * t137 + t185 * t142 + t224 * t205 + (qJD(1) * t177 + t164) * t203, t201 * t136 - t185 * t141 - t231 * t205 + (-t166 + t337) * t203, t140 * t67 - t242 * t71, t130 * t71 - t67 * t139 - t140 * t222 + t242 * t72, t71 * t204, -t72 * t204, 0, -t132 * t130 + t148 * t243 + t112 * t139 + t143 * t72 + t269 * t204 + (-t148 * t242 + t204 * t244) * qJD(4), t112 * t140 - t132 * t242 + t143 * t71 + t148 * t67 - t204 * t234, -t28 * t47 - t29 * t46 + t341 * t5 + t342 * t6 + t36 * t42 - t4 * t93 + t251, t219 * t230 + t238 * t95 - t28 * t5 + t29 * t6 + t3 * t36 + t4 * t37, -t77 * t280 + (t22 * t94 + t47 * t77) * t213 (-t213 * t75 - t317) * t47 + (-t20 - t213 * t23 + (t209 * t75 - t213 * t77) * qJD(6)) * t94, t213 * t248 + t22 * t93 - t280 * t84 + t46 * t77, -t209 * t248 - t23 * t93 - t46 * t75 - t83 * t94, t41 * t93 + t46 * t84, -t247 * t46 + t36 * t23 + t5 * t75 + t9 * t93 + (t14 * t84 + t329 + (-t27 * t93 - t37 * t84 + t330) * qJD(6)) * t213 + t336 * t209, -t13 * t46 + t36 * t22 + t5 * t77 + (-(-qJD(6) * t37 + t14) * t84 - t329 - (-qJD(6) * t27 + t10) * t93 - qJD(6) * t330) * t209 + t336 * t213; 0, 0, 0, -t212 * t298, t291 * t218, 0, 0, 0, t218 * pkin(1) * t212, pkin(1) * t298, -t304, t109, t99, t100, 0, -t164 * t202 + t303 - t255 * t205 + (t211 * t273 - t171) * qJD(3) + t257, t166 * t202 + t292 * t205 + (qJD(3) * t273 + t174) * t215 + t235, t307, t59, t52, t53, 0, t144 * t130 + t349 * t204 + t348, t144 * t242 + t350 * t204 + t347, -t116 * t42 - t117 * t41 + t322 * t342 + t323 * t341 + t249, t4 * t117 - t3 * t116 - t95 * (t202 + t250) + t322 * t29 - t323 * t28, t11, t1, t8, t7, -t326, t114 * t23 + t323 * t75 + (-t263 * t84 - t3) * t213 + t228 * t209 + t275, t114 * t22 + t228 * t213 + t263 * t316 + t323 * t77 + t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t304, t109, t99, t100, 0, -t205 * t241 + t225 + t303, t205 * t258 + t235 - t338, t307, t59, t52, t53, 0, -t130 * t333 + t306 - t261 * t204 + (t210 * t276 - t103) * qJD(4) + t246, -t242 * t333 + t294 * t204 + (qJD(4) * t276 - t68) * t214 + t354, -t158 * t42 - t159 * t41 + t311 * t342 + t312 * t341 + t249, -t3 * t158 + t4 * t159 - t95 * t250 - t312 * t28 + t311 * t29, t11, t1, t8, t7, -t326, t151 * t23 + t312 * t75 + (-t262 * t84 - t3) * t213 + t227 * t209 + t275, t151 * t22 + t213 * t227 + t262 * t316 + t312 * t77 + t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, t59, t52, t53, 0, -t204 * t245 + t348, t204 * t266 + t347 (-t208 * t41 - t310 * t42) * pkin(4) + (t28 - t31) * t342 + (t29 - t30) * t341, t28 * t30 - t29 * t31 + (t208 * t4 + t242 * t95 - t310 * t3) * pkin(4), t11, -t317 * t84 + t353, t8, t239 + t328, -t326, t197 * t23 - t3 * t213 - (-t209 * t31 + t213 * t48) * t84 - t30 * t75 - t26 * t315 - t324 * t196 + t275, t197 * t22 + (t209 * t48 + t213 * t31) * t84 - t30 * t77 - t26 * t352 + t339 * t196 + t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341 ^ 2 - t342 ^ 2, t28 * t341 - t29 * t342 + t219, 0, 0, 0, 0, 0, t239 - t328, -t213 * t84 ^ 2 - t327 - t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75, -t75 ^ 2 + t77 ^ 2, t75 * t84 + t22, t77 * t84 - t23, t41, -t351 * t13 - t209 * t4 - t26 * t77 + t9, -t209 * t10 - t213 * t4 + t351 * t247 + t26 * t75;];
tauc_reg  = t2;
