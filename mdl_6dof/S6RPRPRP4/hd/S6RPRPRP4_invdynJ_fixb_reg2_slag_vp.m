% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:08
% EndTime: 2019-03-09 03:13:16
% DurationCPUTime: 4.32s
% Computational Cost: add. (4460->486), mult. (8613->559), div. (0->0), fcn. (5015->10), ass. (0->266)
t179 = cos(qJ(3));
t282 = qJD(1) * qJD(3);
t260 = t179 * t282;
t176 = sin(qJ(3));
t277 = t176 * qJDD(1);
t211 = t260 + t277;
t104 = qJDD(5) + t211;
t295 = qJD(1) * t176;
t139 = qJD(5) + t295;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t289 = qJD(3) * t178;
t266 = t176 * t289;
t284 = qJD(5) * t179;
t207 = t175 * t284 + t266;
t298 = t139 * t207;
t303 = t178 * t179;
t220 = -t104 * t303 + t298;
t291 = qJD(3) * t175;
t294 = qJD(1) * t179;
t107 = t178 * t294 + t291;
t288 = qJD(3) * t179;
t263 = t175 * t294;
t261 = t176 * t282;
t276 = t179 * qJDD(1);
t353 = -t261 + t276;
t48 = -qJD(5) * t263 + t175 * qJDD(3) + (qJD(3) * qJD(5) + t353) * t178;
t42 = t176 * t48;
t328 = t107 * t288 + t42;
t360 = t220 + t328;
t109 = -t263 + t289;
t224 = t109 * t139;
t225 = t107 * t139;
t47 = qJD(5) * t107 - t178 * qJDD(3) + t175 * t353;
t359 = t178 * (t48 + t224) - (t47 + t225) * t175;
t173 = sin(pkin(9));
t146 = pkin(1) * t173 + pkin(7);
t329 = pkin(4) + t146;
t285 = qJD(5) * t178;
t44 = t48 * t175;
t358 = t107 * t285 + t44;
t122 = t146 * qJD(1);
t159 = t179 * qJD(2);
t81 = t176 * t122 - t159;
t355 = -qJD(4) - t81;
t354 = qJD(2) * qJD(3) + qJDD(1) * t146;
t321 = t178 * t47;
t12 = -t175 * t224 - t321;
t253 = -t179 * qJDD(2) + t122 * t288 + t176 * t354;
t228 = qJDD(4) + t253;
t337 = pkin(3) + pkin(8);
t20 = t211 * pkin(4) - t337 * qJDD(3) + t228;
t138 = pkin(3) * t261;
t318 = qJ(4) * t179;
t239 = pkin(8) * t176 - t318;
t283 = t176 * qJD(4);
t195 = qJD(3) * t239 - t283;
t161 = t176 * qJ(4);
t258 = -pkin(2) - t161;
t210 = -t337 * t179 + t258;
t174 = cos(pkin(9));
t330 = t174 * pkin(1);
t198 = t210 - t330;
t27 = qJD(1) * t195 + qJDD(1) * t198 + t138;
t286 = qJD(5) * t175;
t344 = qJD(4) - t159;
t300 = (pkin(4) * qJD(1) + t122) * t176 + t344;
t52 = -t337 * qJD(3) + t300;
t58 = t198 * qJD(1);
t257 = t175 * t27 - t178 * t20 + t58 * t285 + t52 * t286;
t336 = pkin(5) * t104;
t2 = qJDD(6) + t257 - t336;
t16 = t175 * t52 + t178 * t58;
t14 = qJ(6) * t139 + t16;
t324 = t139 * t14;
t352 = -t2 + t324;
t323 = t139 * t16;
t351 = -t257 + t323;
t331 = g(3) * t176;
t167 = qJ(1) + pkin(9);
t154 = sin(t167);
t155 = cos(t167);
t345 = -g(1) * t155 - g(2) * t154;
t187 = qJD(5) * t139 * t337 + t179 * t345 - t331;
t168 = qJDD(3) * qJ(4);
t169 = qJD(3) * qJD(4);
t290 = qJD(3) * t176;
t254 = -t176 * qJDD(2) + t122 * t290 - t179 * t354;
t32 = -t168 - t169 + t254;
t21 = pkin(4) * t353 - t32;
t5 = pkin(5) * t48 + qJ(6) * t47 - qJD(6) * t109 + t21;
t350 = t5 + t187;
t91 = t175 * t104;
t216 = -t139 * t285 - t91;
t82 = qJD(2) * t176 + t122 * t179;
t67 = -qJD(3) * qJ(4) - t82;
t64 = -qJD(3) * pkin(3) - t355;
t348 = t48 - t224;
t293 = qJD(3) * t107;
t338 = t139 ^ 2;
t92 = t178 * t104;
t347 = -t175 * t338 - t293 + t92;
t319 = pkin(1) * qJDD(1);
t346 = t176 * t337;
t315 = t104 * t337;
t151 = pkin(4) * t294;
t56 = t151 - t67;
t343 = t139 * t56 - t315;
t166 = g(3) * t179;
t275 = t175 * t20 + t178 * t27 + t52 * t285;
t306 = t175 * t176;
t75 = t154 * t178 + t155 * t306;
t77 = -t154 * t306 + t155 * t178;
t342 = -g(1) * t75 + g(2) * t77 + (-qJD(5) * t58 + t166) * t175 + t275;
t162 = t179 * pkin(8);
t147 = -pkin(2) - t330;
t164 = t179 * pkin(3);
t297 = t164 + t161;
t93 = t147 - t297;
t78 = -t162 + t93;
t96 = t329 * t176;
t327 = t175 * t96 + t178 * t78;
t152 = pkin(3) * t290;
t65 = t152 + t195;
t89 = t329 * t288;
t10 = -qJD(5) * t327 - t175 * t65 + t178 * t89;
t15 = -t175 * t58 + t178 * t52;
t233 = t15 * t175 - t16 * t178;
t3 = -t286 * t58 + t275;
t341 = -qJD(5) * t233 + t3 * t175 - t178 * t257;
t316 = t104 * qJ(6);
t1 = t139 * qJD(6) + t3 + t316;
t299 = qJD(6) - t15;
t13 = -pkin(5) * t139 + t299;
t235 = t13 * t175 + t14 * t178;
t340 = qJD(5) * t235 + t1 * t175 - t2 * t178;
t339 = t109 ^ 2;
t335 = g(1) * t154;
t332 = g(2) * t155;
t62 = t151 + t82;
t153 = pkin(3) * t295;
t83 = qJD(1) * t239 + t153;
t31 = t175 * t62 + t178 * t83;
t241 = pkin(5) * t178 + qJ(6) * t175;
t221 = -pkin(4) - t241;
t326 = qJD(5) * t241 - t178 * qJD(6) - (qJD(1) * t221 - t122) * t176 + t344;
t43 = t176 * t47;
t325 = t109 * t288 - t43;
t322 = t15 * t139;
t320 = t337 * t47;
t317 = qJDD(3) * pkin(3);
t314 = t107 * t178;
t313 = t109 * t107;
t312 = t109 * t175;
t311 = t109 * t337;
t310 = t154 * t176;
t309 = t154 * t179;
t308 = t155 * t176;
t307 = t155 * t179;
t305 = t175 * t179;
t304 = t176 * t178;
t26 = pkin(5) * t107 - qJ(6) * t109 + t56;
t302 = t26 * qJD(3);
t301 = t56 * qJD(3);
t97 = t329 * t179;
t171 = t176 ^ 2;
t172 = t179 ^ 2;
t296 = t171 - t172;
t123 = qJD(1) * t147;
t292 = qJD(3) * t109;
t121 = qJDD(1) * t147;
t280 = qJDD(1) * t171;
t279 = qJDD(1) * t172;
t278 = qJDD(3) * t146;
t268 = t178 * t295;
t274 = t107 * t268 + t358;
t273 = t139 * t268 - t216;
t180 = cos(qJ(1));
t271 = t180 * pkin(1) + t155 * pkin(2) + t154 * pkin(7);
t270 = t162 + t297;
t269 = -g(1) * t308 - g(2) * t310 + t166;
t267 = t175 * t290;
t264 = t139 * t294;
t262 = t107 ^ 2 - t339;
t177 = sin(qJ(1));
t259 = -t177 * pkin(1) + t155 * pkin(7);
t252 = t109 * t266;
t251 = t176 * t260;
t74 = t154 * t175 - t155 * t304;
t76 = t154 * t304 + t155 * t175;
t250 = -g(1) * t76 - g(2) * t74;
t249 = -g(1) * t77 - g(2) * t75;
t248 = t155 * pkin(4) + t259;
t124 = qJ(4) * t309;
t126 = qJ(4) * t307;
t247 = -g(1) * t126 - g(2) * t124;
t245 = g(1) * t177 - g(2) * t180;
t242 = t337 * t358 - t269;
t240 = pkin(5) * t175 - qJ(6) * t178;
t182 = qJD(3) ^ 2;
t238 = t146 * t182 + t332;
t236 = -t13 * t178 + t14 * t175;
t234 = t15 * t178 + t16 * t175;
t30 = -t175 * t83 + t178 * t62;
t40 = -t175 * t78 + t178 * t96;
t229 = pkin(3) * t307 + t155 * t161 + t271;
t71 = t107 * t267;
t223 = -t71 - t252;
t222 = t258 - t164;
t9 = t175 * t89 + t178 * t65 + t96 * t285 - t286 * t78;
t215 = -qJ(4) * t288 - t283;
t214 = -t81 * qJD(3) + t254;
t213 = t345 + (t279 + t280) * t146;
t212 = t107 * t294 - t273;
t208 = t154 * pkin(4) + pkin(8) * t307 + t229;
t206 = -t178 * t284 + t267;
t35 = t228 - t317;
t205 = t222 - t330;
t204 = t139 * t26 - t315;
t68 = t205 * qJD(1);
t203 = t278 + (-qJD(1) * t93 - t68) * qJD(3);
t202 = -0.2e1 * t121 - t238;
t201 = g(1) * t74 - g(2) * t76 + g(3) * t303 - t257;
t200 = 0.2e1 * qJD(3) * t123 - t278;
t199 = t210 * t335;
t197 = qJD(3) * t82 - t253 - t269;
t194 = t44 + (-t312 + t314) * qJD(5);
t38 = qJD(1) * t215 + qJDD(1) * t205 + t138;
t90 = t152 + t215;
t192 = qJD(1) * t90 + qJDD(1) * t93 + t238 + t38;
t191 = t35 * t176 - t32 * t179 + (t176 * t67 + t179 * t64) * qJD(3);
t190 = t253 * t176 - t254 * t179 + (-t176 * t82 + t179 * t81) * qJD(3);
t189 = -t12 - t274;
t188 = t109 * t26 + qJDD(6) - t201;
t186 = -t107 * t207 + t303 * t48;
t185 = t21 + t187;
t183 = qJD(1) ^ 2;
t137 = t176 * t183 * t179;
t131 = g(1) * t309;
t119 = t296 * t183;
t118 = qJDD(3) * t179 - t182 * t176;
t117 = qJDD(3) * t176 + t179 * t182;
t116 = qJ(4) + t240;
t111 = -qJ(4) * t294 + t153;
t100 = t139 * t267;
t95 = -0.2e1 * t251 + t279;
t94 = 0.2e1 * t251 + t280;
t88 = t329 * t290;
t63 = 0.2e1 * t176 * t276 - 0.2e1 * t282 * t296;
t59 = t68 * t295;
t55 = t179 * t241 + t97;
t54 = t104 * t176 + t139 * t288;
t53 = pkin(5) * t109 + qJ(6) * t107;
t36 = t47 * t303;
t34 = -pkin(5) * t176 - t40;
t33 = qJ(6) * t176 + t327;
t29 = -pkin(5) * t294 - t30;
t28 = qJ(6) * t294 + t31;
t25 = t225 - t47;
t24 = (-qJD(5) * t240 + qJD(6) * t175) * t179 + (-t146 + t221) * t290;
t19 = -t139 * t286 + t92 + (-t109 * t179 - t139 * t306) * qJD(1);
t11 = t109 * t206 + t305 * t47;
t8 = -pkin(5) * t288 - t10;
t7 = qJ(6) * t288 + qJD(6) * t176 + t9;
t6 = t179 * t216 + t100 + t325;
t4 = [0, 0, 0, 0, 0, qJDD(1), t245, g(1) * t180 + g(2) * t177, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t174 * t319 - t332 + t335, -0.2e1 * t173 * t319 - t345, 0 (t245 + (t173 ^ 2 + t174 ^ 2) * t319) * pkin(1), t94, t63, t117, t95, t118, 0, t176 * t200 + t179 * t202 + t131, t200 * t179 + (-t202 - t335) * t176, t190 + t213, t121 * t147 - g(1) * (-t154 * pkin(2) + t259) - g(2) * t271 + t190 * t146, 0, -t117, -t118, t94, t63, t95, t191 + t213, t176 * t203 + t179 * t192 - t131, t203 * t179 + (-t192 + t335) * t176, -g(1) * t259 - g(2) * t229 + t146 * t191 - t222 * t335 + t38 * t93 + t68 * t90, t11, t252 + t36 - t71 + (t44 + (t312 + t314) * qJD(5)) * t179, t6, t186, -t42 + (-t293 - t92) * t179 + t298, t54, t10 * t139 + t40 * t104 - t88 * t107 + t97 * t48 + (-t289 * t56 - t257) * t176 + (qJD(3) * t15 + t21 * t178 - t286 * t56) * t179 + t249, -t327 * t104 - t88 * t109 - t9 * t139 - t97 * t47 + (t291 * t56 - t3) * t176 + (-qJD(3) * t16 - t21 * t175 - t285 * t56) * t179 - t250, -t10 * t109 - t9 * t107 + t40 * t47 - t327 * t48 + t131 - t233 * t290 + (qJD(5) * t234 - t175 * t257 - t178 * t3 - t332) * t179, -g(1) * t248 - g(2) * t208 + t15 * t10 + t16 * t9 + t21 * t97 - t257 * t40 + t3 * t327 - t56 * t88 - t199, t11, t6, -t36 + (-t107 * t284 - t109 * t290) * t178 + (t107 * t290 + (-qJD(5) * t109 - t48) * t179) * t175, t54, -t220 + t328, t186, -t34 * t104 + t24 * t107 - t8 * t139 + t55 * t48 + (-t26 * t289 - t2) * t176 + (-qJD(3) * t13 + t5 * t178 - t26 * t286) * t179 + t249, -t7 * t107 + t8 * t109 - t33 * t48 - t34 * t47 + t131 + t235 * t290 + (qJD(5) * t236 - t1 * t178 - t175 * t2 - t332) * t179, t33 * t104 - t24 * t109 + t7 * t139 + t55 * t47 + (-t26 * t291 + t1) * t176 + (qJD(3) * t14 + t5 * t175 + t26 * t285) * t179 + t250, t1 * t33 + t14 * t7 + t5 * t55 + t26 * t24 + t2 * t34 + t13 * t8 - g(1) * (t77 * pkin(5) + t76 * qJ(6) + t248) - g(2) * (pkin(5) * t75 + qJ(6) * t74 + t208) - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t118, -t117, 0, -t254 * t176 - t253 * t179 - g(3) + (t176 * t81 + t179 * t82) * qJD(3), 0, 0, 0, 0, 0, 0, 0, -t118, t117, -t32 * t176 - t35 * t179 - g(3) + (t176 * t64 - t179 * t67) * qJD(3), 0, 0, 0, 0, 0, 0, t360, t104 * t305 - t139 * t206 + t325 (t194 - t321) * t179 + t223, -g(3) + (qJD(3) * t234 + t21) * t176 + (t301 - t341) * t179, 0, 0, 0, 0, 0, 0, t360, t179 * t194 + t223 - t36, t43 + t100 + (t216 - t292) * t179, -g(3) + (qJD(3) * t236 + t5) * t176 + (t302 - t340) * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t119, t277, t137, t276, qJDD(3), -t123 * t295 + t197, t331 + (-qJD(1) * t123 - t345) * t179 + t214, 0, 0, qJDD(3), -t277, -t276, -t137, t119, t137 (-pkin(3) * t176 + t318) * qJDD(1), -t111 * t294 + qJDD(4) - t197 - 0.2e1 * t317 + t59, 0.2e1 * t168 + 0.2e1 * t169 + (qJD(1) * t111 - g(3)) * t176 + (qJD(1) * t68 + t345) * t179 - t214, -t32 * qJ(4) - t35 * pkin(3) - t68 * t111 - t64 * t82 - g(1) * (-pkin(3) * t308 + t126) - g(2) * (-pkin(3) * t310 + t124) - g(3) * t297 + t355 * t67, t12, -t359, t19, t274, t212, -t264, qJ(4) * t48 + t300 * t107 - t30 * t139 - t15 * t294 + t185 * t175 + t178 * t343, -qJ(4) * t47 + t300 * t109 + t31 * t139 + t16 * t294 - t175 * t343 + t185 * t178, t31 * t107 + t30 * t109 + (-t320 - t351) * t178 + (t15 * t295 - t3 + (t15 - t311) * qJD(5)) * t175 + t242, -g(3) * t270 + t21 * qJ(4) - t15 * t30 - t16 * t31 + t300 * t56 - t337 * t341 - t345 * t346 + t247, t12, t19, t359, -t264, -t212, t178 * t225 + t44, t326 * t107 + t116 * t48 + t13 * t294 + t29 * t139 + t175 * t350 + t204 * t178, t28 * t107 - t29 * t109 + (-t320 - t352) * t178 + (-t13 * t295 - t1 + (-t13 - t311) * qJD(5)) * t175 + t242, -t326 * t109 + t116 * t47 - t28 * t139 - t14 * t294 + t204 * t175 - t178 * t350, t5 * t116 - t14 * t28 - t13 * t29 - g(3) * (t176 * t240 + t270) + t326 * t26 - t340 * t337 + t247 + t345 * (t179 * t240 - t346); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, qJDD(3) + t137, -t171 * t183 - t182, qJD(3) * t67 + t269 + t35 + t59, 0, 0, 0, 0, 0, 0, t347, -t178 * t338 - t292 - t91, t189, -t301 + t351 * t178 + (t3 - t322) * t175 + t269, 0, 0, 0, 0, 0, 0, t347, t189, t273 + t292, -t302 + t352 * t178 + (t13 * t139 + t1) * t175 + t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, -t262, t25, -t313, -t348, t104, -t109 * t56 + t201 + t323, t56 * t107 + t322 - t342, 0, 0, t313, t25, t262, t104, t348, -t313, -t107 * t53 - t188 + t323 + 0.2e1 * t336, pkin(5) * t47 - t48 * qJ(6) + (t14 - t16) * t109 + (t13 - t299) * t107, 0.2e1 * t316 - t26 * t107 + t53 * t109 + (0.2e1 * qJD(6) - t15) * t139 + t342, t1 * qJ(6) - t2 * pkin(5) - t26 * t53 - t13 * t16 - g(1) * (-pkin(5) * t74 + qJ(6) * t75) - g(2) * (pkin(5) * t76 - qJ(6) * t77) + t241 * t166 + t299 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 + t313, t25, -t338 - t339, t188 - t324 - t336;];
tau_reg  = t4;
