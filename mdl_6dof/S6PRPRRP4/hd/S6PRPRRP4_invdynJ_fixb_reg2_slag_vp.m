% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:26
% EndTime: 2019-03-08 20:12:35
% DurationCPUTime: 5.75s
% Computational Cost: add. (7517->552), mult. (17703->697), div. (0->0), fcn. (14511->14), ass. (0->278)
t212 = cos(pkin(11));
t370 = cos(qJ(4));
t296 = qJD(4) * t370;
t209 = sin(pkin(11));
t216 = sin(qJ(4));
t315 = qJD(4) * t216;
t299 = t209 * t315;
t165 = -t212 * t296 + t299;
t170 = t209 * t370 + t216 * t212;
t166 = t170 * qJD(4);
t217 = sin(qJ(2));
t211 = sin(pkin(6));
t318 = qJD(1) * t211;
t301 = t217 * t318;
t394 = pkin(4) * t166 + pkin(9) * t165 - t301;
t302 = t370 * t212;
t326 = t216 * t209;
t255 = t302 - t326;
t219 = cos(qJ(2));
t328 = t211 * t219;
t237 = t255 * t328;
t366 = pkin(8) + qJ(3);
t181 = t366 * t209;
t182 = t366 * t212;
t386 = -t370 * t181 - t216 * t182;
t388 = -qJD(1) * t237 + t255 * qJD(3) + qJD(4) * t386;
t194 = qJD(2) * t302;
t232 = -qJD(2) * t299 + qJDD(2) * t170;
t385 = qJD(4) * t194 + t232;
t393 = qJD(4) * qJD(5) + t385;
t161 = qJD(2) * t326 - t194;
t152 = qJD(5) + t161;
t163 = t170 * qJD(2);
t215 = sin(qJ(5));
t218 = cos(qJ(5));
t124 = -t218 * qJD(4) + t163 * t215;
t126 = qJD(4) * t215 + t163 * t218;
t346 = t126 * t215;
t349 = t124 * t218;
t268 = t346 + t349;
t313 = qJD(5) * t218;
t51 = -t218 * qJDD(4) + t163 * t313 + t393 * t215;
t358 = t51 * t218;
t314 = qJD(5) * t215;
t50 = -t215 * qJDD(4) + t163 * t314 - t393 * t218;
t359 = t50 * t215;
t392 = t170 * ((t124 * t215 - t126 * t218) * qJD(5) - t358 + t359) + t268 * t165;
t121 = -t216 * t181 + t182 * t370;
t198 = t212 * pkin(3) + pkin(2);
t99 = -pkin(4) * t255 - pkin(9) * t170 - t198;
t365 = -t121 * t314 + t215 * t394 + t388 * t218 + t99 * t313;
t175 = qJD(2) * qJ(3) + t301;
t213 = cos(pkin(6));
t317 = qJD(1) * t213;
t130 = t212 * t175 + t209 * t317;
t363 = pkin(8) * qJD(2);
t115 = t212 * t363 + t130;
t192 = t212 * t317;
t114 = t192 + (-t175 - t363) * t209;
t298 = t219 * t318;
t309 = qJDD(2) * qJ(3);
t311 = qJDD(1) * t211;
t128 = t217 * t311 + t309 + (qJD(3) + t298) * qJD(2);
t310 = qJDD(1) * t213;
t190 = t212 * t310;
t80 = t190 + (-pkin(8) * qJDD(2) - t128) * t209;
t307 = t212 * qJDD(2);
t94 = t212 * t128 + t209 * t310;
t81 = pkin(8) * t307 + t94;
t305 = -t114 * t296 - t216 * t80 - t370 * t81;
t15 = -t115 * t315 - t305;
t12 = qJDD(4) * pkin(9) + t15;
t308 = t209 * qJDD(2);
t274 = -qJDD(2) * t302 + t216 * t308;
t105 = qJD(2) * t166 + t274;
t294 = t219 * t311;
t316 = qJD(2) * t217;
t295 = qJD(1) * t316;
t383 = t211 * t295 + qJDD(3);
t258 = -t294 + t383;
t350 = qJDD(2) * pkin(2);
t139 = t258 - t350;
t39 = -pkin(3) * t307 + t105 * pkin(4) - pkin(9) * t385 + t139;
t57 = t216 * t114 + t115 * t370;
t53 = qJD(4) * pkin(9) + t57;
t275 = qJD(3) - t298;
t151 = -qJD(2) * t198 + t275;
t65 = t161 * pkin(4) - t163 * pkin(9) + t151;
t290 = t215 * t12 - t218 * t39 + t53 * t313 + t65 * t314;
t96 = qJDD(5) + t105;
t374 = pkin(5) * t96;
t2 = qJDD(6) + t290 - t374;
t27 = t215 * t65 + t218 * t53;
t20 = qJ(6) * t152 + t27;
t362 = t152 * t20;
t391 = -t2 + t362;
t360 = t152 * t27;
t390 = -t290 + t360;
t389 = t218 * t121 + t215 * t99;
t238 = t170 * t328;
t353 = -qJD(1) * t238 + qJD(3) * t170 + qJD(4) * t121;
t220 = qJD(2) ^ 2;
t245 = (qJDD(2) * t219 - t217 * t220) * t211;
t352 = cos(pkin(10));
t291 = t352 * t219;
t210 = sin(pkin(10));
t331 = t210 * t217;
t157 = -t213 * t291 + t331;
t292 = t352 * t217;
t330 = t210 * t219;
t159 = t213 * t330 + t292;
t282 = g(1) * t159 + g(2) * t157;
t242 = -g(3) * t328 + t282;
t277 = pkin(5) * t218 + qJ(6) * t215;
t384 = t163 * qJD(4);
t112 = t216 * t115;
t56 = t114 * t370 - t112;
t52 = -qJD(4) * pkin(4) - t56;
t31 = t124 * pkin(5) - t126 * qJ(6) + t52;
t373 = pkin(9) * t96;
t382 = t152 * t31 - t373;
t339 = t165 * t215;
t254 = t170 * t313 - t339;
t337 = t170 * t215;
t381 = t166 * t124 + t152 * t254 - t255 * t51 + t337 * t96;
t261 = -t139 + t282;
t380 = t211 * (-g(3) * t219 + t295) + t261 + t350;
t266 = (-t175 * t209 + t192) * t209 - t130 * t212;
t379 = t219 * t266 - (-qJD(2) * pkin(2) + t275) * t217;
t364 = -qJD(5) * t389 - t388 * t215 + t218 * t394;
t377 = t126 ^ 2;
t376 = t152 ^ 2;
t375 = t163 ^ 2;
t372 = qJ(6) * t166 - qJD(6) * t255 + t365;
t371 = -t166 * pkin(5) - t364;
t208 = pkin(11) + qJ(4);
t201 = cos(t208);
t369 = pkin(4) * t201;
t367 = pkin(9) * t126;
t97 = pkin(4) * t163 + pkin(9) * t161;
t35 = t215 * t97 + t218 * t56;
t26 = -t215 * t53 + t218 * t65;
t361 = t152 * t26;
t87 = t215 * t96;
t88 = t218 * t96;
t357 = t96 * qJ(6);
t356 = -t124 * t313 - t215 * t51;
t276 = pkin(5) * t215 - qJ(6) * t218;
t355 = -t215 * qJD(6) + t152 * t276 - t57;
t354 = -t276 * t165 + (qJD(5) * t277 - qJD(6) * t218) * t170 + t353;
t286 = t124 * t152;
t348 = t126 * t124;
t288 = t126 * t152;
t347 = t126 * t163;
t345 = t152 * t161;
t344 = t152 * t163;
t200 = sin(t208);
t343 = t157 * t200;
t342 = t159 * t200;
t341 = t163 * t124;
t340 = t163 * t161;
t338 = t165 * t218;
t336 = t170 * t218;
t334 = t201 * t215;
t333 = t201 * t218;
t332 = t210 * t211;
t329 = t211 * t217;
t325 = t218 * t219;
t324 = t219 * t220;
t323 = qJD(6) - t26;
t322 = qJDD(1) - g(3);
t158 = t213 * t292 + t330;
t321 = -t157 * t198 + t158 * t366;
t160 = -t213 * t331 + t291;
t320 = -t159 * t198 + t160 * t366;
t205 = t209 ^ 2;
t207 = t212 ^ 2;
t319 = t205 + t207;
t306 = g(3) * t329;
t304 = t200 * t328;
t303 = t211 * t325;
t188 = t215 * t328;
t300 = t211 * t316;
t297 = t124 ^ 2 - t377;
t293 = t211 * t352;
t116 = t126 * t314;
t289 = -t218 * t50 - t116;
t285 = t114 * t315 + t115 * t296 + t216 * t81 - t370 * t80;
t284 = -pkin(9) * t343 - t157 * t369 + t321;
t283 = -pkin(9) * t342 - t159 * t369 + t320;
t281 = g(1) * t160 + g(2) * t158;
t280 = t198 * t328 + t329 * t366;
t279 = t198 * qJDD(2);
t278 = (qJD(5) * t124 - t50) * pkin(9);
t19 = -pkin(5) * t152 + t323;
t273 = t19 * t218 - t20 * t215;
t272 = t215 * t27 + t218 * t26;
t34 = -t215 * t56 + t218 * t97;
t47 = -t121 * t215 + t218 * t99;
t263 = -t152 * t314 - t215 * t345 + t88;
t262 = t152 * t313 + t218 * t345 + t87;
t155 = -t209 * t329 + t212 * t213;
t156 = t209 * t213 + t212 * t329;
t85 = t216 * t155 + t156 * t370;
t67 = t215 * t85 + t303;
t3 = t218 * t12 + t215 * t39 + t65 * t313 - t314 * t53;
t257 = t152 * t52 - t373;
t256 = t155 * t370 - t216 * t156;
t253 = -t170 * t314 - t338;
t13 = -qJDD(4) * pkin(4) + t285;
t54 = qJD(2) * t237 + qJD(4) * t256;
t24 = -qJD(5) * t67 + t215 * t300 + t218 * t54;
t25 = -qJD(5) * t188 + t215 * t54 - t218 * t300 + t313 * t85;
t68 = t218 * t85 - t188;
t252 = -t24 * t124 + t126 * t25 - t50 * t67 - t68 * t51;
t55 = qJD(2) * t238 + qJD(4) * t85;
t251 = t55 * t124 - t152 * t25 - t256 * t51 - t67 * t96;
t131 = t188 * t201 - t218 * t329;
t76 = -t157 * t334 - t158 * t218;
t78 = -t159 * t334 - t160 * t218;
t250 = -g(1) * t78 - g(2) * t76 - g(3) * t131;
t132 = (t201 * t325 + t215 * t217) * t211;
t77 = -t157 * t333 + t158 * t215;
t79 = -t159 * t333 + t160 * t215;
t249 = -g(1) * t79 - g(2) * t77 - g(3) * t132;
t248 = pkin(9) * t304 + t328 * t369 + t280;
t106 = -t158 * t200 - t201 * t293;
t108 = -t160 * t200 + t201 * t332;
t142 = -t200 * t329 + t201 * t213;
t247 = g(1) * t108 + g(2) * t106 + g(3) * t142;
t107 = t158 * t201 - t200 * t293;
t109 = t160 * t201 + t200 * t332;
t143 = t200 * t213 + t201 * t329;
t246 = -g(1) * t109 - g(2) * t107 - g(3) * t143;
t244 = t215 * t286 - t358;
t243 = -pkin(9) * t358 + t246;
t241 = -t263 - t341;
t236 = t242 * t200;
t235 = t242 + t294;
t93 = -t128 * t209 + t190;
t234 = -t93 * t209 + t94 * t212 - t281;
t233 = t126 * t55 - t152 * t24 + t256 * t50 - t68 * t96;
t110 = t143 * t215 + t303;
t61 = t107 * t215 - t157 * t218;
t63 = t109 * t215 - t159 * t218;
t231 = g(1) * t63 + g(2) * t61 + g(3) * t110 - t290;
t230 = pkin(9) * qJD(5) * t152 + t247;
t5 = pkin(5) * t51 + qJ(6) * t50 - qJD(6) * t126 + t13;
t229 = -t230 - t5;
t228 = t13 + t230;
t227 = -t235 + t383;
t226 = t124 * t254 + t337 * t51;
t111 = t143 * t218 - t188;
t62 = t107 * t218 + t157 * t215;
t64 = t109 * t218 + t159 * t215;
t225 = -g(1) * t64 - g(2) * t62 - g(3) * t111 + t3;
t224 = t126 * t31 + qJDD(6) - t231;
t221 = -t51 + t288;
t177 = -pkin(4) - t277;
t154 = t161 ^ 2;
t141 = t142 * pkin(4);
t118 = -t279 + t258;
t101 = t108 * pkin(4);
t100 = t106 * pkin(4);
t66 = pkin(5) * t126 + qJ(6) * t124;
t58 = t170 * t276 - t386;
t43 = t152 * t166 - t255 * t96;
t42 = pkin(5) * t255 - t47;
t41 = -qJ(6) * t255 + t389;
t32 = -t50 + t286;
t29 = -pkin(5) * t163 - t34;
t28 = qJ(6) * t163 + t35;
t23 = t262 - t347;
t21 = t218 * t288 - t359;
t9 = t126 * t253 - t336 * t50;
t6 = t126 * t166 + t152 * t253 + t255 * t50 + t336 * t96;
t1 = qJD(6) * t152 + t3 + t357;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t322, 0, 0, 0, 0, 0, 0, t245 (-qJDD(2) * t217 - t324) * t211, 0, -g(3) + (t213 ^ 2 + (t217 ^ 2 + t219 ^ 2) * t211 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t212 * t245, -t209 * t245, t319 * t211 * t324 + (-t155 * t209 + t156 * t212) * qJDD(2), t93 * t155 + t94 * t156 - g(3) + (-qJD(2) * t379 - t139 * t219) * t211, 0, 0, 0, 0, 0, 0, -t55 * qJD(4) + t256 * qJDD(4) + (-t105 * t219 + t161 * t316) * t211, -t54 * qJD(4) - t85 * qJDD(4) + t163 * t300 - t328 * t385, -t85 * t105 - t54 * t161 + t55 * t163 - t256 * t385, t15 * t85 - t285 * t256 + t57 * t54 - t56 * t55 - g(3) + (-t118 * t219 + t151 * t316) * t211, 0, 0, 0, 0, 0, 0, t251, t233, t252, -t13 * t256 + t24 * t27 - t25 * t26 + t290 * t67 + t3 * t68 + t52 * t55 - g(3), 0, 0, 0, 0, 0, 0, t251, t252, -t233, t1 * t68 + t19 * t25 + t2 * t67 + t20 * t24 - t256 * t5 + t31 * t55 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t235, -t322 * t329 + t281, 0, 0, t205 * qJDD(2), 0.2e1 * t209 * t307, 0, t207 * qJDD(2), 0, 0, t380 * t212, -t380 * t209, -t306 + t234 + (qJD(2) * t275 + t309) * t319, -t266 * qJD(3) + t261 * pkin(2) + t234 * qJ(3) + (-g(3) * (pkin(2) * t219 + qJ(3) * t217) + t379 * qJD(1)) * t211, -t163 * t165 + t170 * t385, -t170 * t105 + t165 * t161 - t163 * t166 + t255 * t385, -qJD(4) * t165 + qJDD(4) * t170, -t105 * t255 + t161 * t166, -qJD(4) * t166 + qJDD(4) * t255, 0, -qJD(4) * t353 + qJDD(4) * t386 - t198 * t105 - t118 * t255 + t151 * t166 - t161 * t301 + t201 * t242, -g(1) * t342 - g(2) * t343 + g(3) * t304 - qJD(4) * t388 - t121 * qJDD(4) + t118 * t170 - t151 * t165 - t163 * t301 - t198 * t385, -t121 * t105 + t15 * t255 - t161 * t388 + t163 * t353 + t56 * t165 - t57 * t166 + t170 * t285 - t385 * t386 - t281 - t306, -g(1) * t320 - g(2) * t321 - g(3) * t280 - t118 * t198 + t15 * t121 - t151 * t301 - t285 * t386 - t353 * t56 + t388 * t57, t9, t392, t6, t226, -t381, t43, -t52 * t339 - t386 * t51 + t26 * t166 + t290 * t255 + t47 * t96 + (t13 * t215 + t313 * t52) * t170 + t364 * t152 + t353 * t124 + t249, -t52 * t338 + t386 * t50 - t27 * t166 + t3 * t255 - t389 * t96 + (t13 * t218 - t314 * t52) * t170 - t365 * t152 + t353 * t126 - t250, t47 * t50 - t389 * t51 + t272 * t165 - t364 * t126 - t365 * t124 + t236 + (-t3 * t215 + t290 * t218 + (t215 * t26 - t218 * t27) * qJD(5)) * t170, -g(1) * t283 - g(2) * t284 - g(3) * t248 - t13 * t386 + t26 * t364 + t27 * t365 - t290 * t47 + t3 * t389 + t353 * t52, t9, t6, -t392, t43, t381, t226, -t31 * t339 - t19 * t166 + t2 * t255 - t42 * t96 + t58 * t51 + (t5 * t215 + t31 * t313) * t170 - t371 * t152 + t354 * t124 + t249, -t41 * t51 - t42 * t50 - t273 * t165 + t371 * t126 - t372 * t124 + t236 + (-t1 * t215 + t2 * t218 + (-t19 * t215 - t20 * t218) * qJD(5)) * t170, t31 * t338 - t1 * t255 + t20 * t166 + t41 * t96 + t58 * t50 + (-t5 * t218 + t31 * t314) * t170 + t372 * t152 - t354 * t126 + t250, t1 * t41 + t5 * t58 + t2 * t42 - g(1) * (pkin(5) * t79 + qJ(6) * t78 + t283) - g(2) * (pkin(5) * t77 + qJ(6) * t76 + t284) - g(3) * (pkin(5) * t132 + qJ(6) * t131 + t248) + t354 * t31 + t372 * t20 + t371 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, t308, -t319 * t220, qJD(2) * t266 + t227 - t350, 0, 0, 0, 0, 0, 0, t274 + 0.2e1 * t384 (t194 - t161) * qJD(4) + t232, -t154 - t375, t57 * t161 + t56 * t163 + t227 - t279, 0, 0, 0, 0, 0, 0, t263 - t341, -t218 * t376 - t347 - t87 (-t124 * t161 + t50) * t218 + t215 * t288 + t356, -t52 * t163 + t390 * t218 + (t3 - t361) * t215 - t242, 0, 0, 0, 0, 0, 0, -t215 * t376 - t341 + t88 (t346 - t349) * t161 - t289 + t356, t262 + t347, -t31 * t163 + t391 * t218 + (t152 * t19 + t1) * t215 - t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, -t154 + t375 (t194 + t161) * qJD(4) + t232, -t340, -t274, qJDD(4), qJD(4) * t57 - t151 * t163 - t247 - t285, t151 * t161 + (t56 + t112) * qJD(4) - t246 + t305, 0, 0, t21, -t161 * t268 + t289 + t356, t23, t244, -t241, -t344, -pkin(4) * t51 - t57 * t124 - t34 * t152 - t26 * t163 + t215 * t257 - t218 * t228, pkin(4) * t50 - t57 * t126 + t35 * t152 + t27 * t163 + t215 * t228 + t218 * t257, t124 * t35 + t126 * t34 + (-t161 * t26 + t3 + (-t26 + t367) * qJD(5)) * t218 + (t278 - t390) * t215 + t243, -t13 * pkin(4) - g(1) * t101 - g(2) * t100 - g(3) * t141 - t26 * t34 - t27 * t35 - t52 * t57 + (-qJD(5) * t272 + t215 * t290 + t3 * t218 + t246) * pkin(9), t21, t23, t116 + (t126 * t161 + t51) * t215 + (t50 + t286) * t218, -t344, t241, t244, t355 * t124 + t29 * t152 + t19 * t163 + t177 * t51 + t215 * t382 + t229 * t218, t124 * t28 - t126 * t29 + (t161 * t19 + t1 + (t19 + t367) * qJD(5)) * t218 + (t278 - t391) * t215 + t243, -t355 * t126 - t28 * t152 - t20 * t163 + t177 * t50 + t229 * t215 - t218 * t382, t5 * t177 - t20 * t28 - t19 * t29 - g(1) * (t108 * t277 + t101) - g(2) * (t106 * t277 + t100) - g(3) * (t142 * t277 + t141) + t355 * t31 + (qJD(5) * t273 + t1 * t218 + t2 * t215 + t246) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, -t297, t32, -t348, t221, t96, -t126 * t52 + t231 + t360, t124 * t52 - t225 + t361, 0, 0, t348, t32, t297, t96, -t221, -t348, -t124 * t66 - t224 + t360 + 0.2e1 * t374, pkin(5) * t50 - t51 * qJ(6) + (t20 - t27) * t126 + (t19 - t323) * t124, 0.2e1 * t357 - t31 * t124 + t66 * t126 + (0.2e1 * qJD(6) - t26) * t152 + t225, t1 * qJ(6) - t2 * pkin(5) - t31 * t66 - t19 * t27 - g(1) * (-pkin(5) * t63 + qJ(6) * t64) - g(2) * (-pkin(5) * t61 + qJ(6) * t62) - g(3) * (-pkin(5) * t110 + qJ(6) * t111) + t323 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) - t274 + t348 - t384, t32, -t376 - t377, t224 - t362 - t374;];
tau_reg  = t4;