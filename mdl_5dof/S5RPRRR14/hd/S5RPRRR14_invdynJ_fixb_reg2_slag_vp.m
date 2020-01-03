% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR14_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:49
% EndTime: 2019-12-31 19:20:15
% DurationCPUTime: 13.87s
% Computational Cost: add. (18708->658), mult. (60180->912), div. (0->0), fcn. (51430->14), ass. (0->292)
t247 = cos(pkin(5));
t246 = cos(pkin(11));
t413 = sin(qJ(1));
t350 = t413 * t246;
t244 = sin(pkin(11));
t253 = cos(qJ(1));
t378 = t253 * t244;
t196 = t247 * t378 + t350;
t250 = sin(qJ(3));
t394 = sin(pkin(6));
t414 = cos(qJ(3));
t320 = t394 * t414;
t245 = sin(pkin(5));
t382 = t245 * t253;
t238 = t413 * t244;
t377 = t253 * t246;
t313 = t247 * t377 - t238;
t395 = cos(pkin(6));
t418 = t313 * t395;
t119 = -t196 * t250 - t320 * t382 + t414 * t418;
t248 = sin(qJ(5));
t251 = cos(qJ(5));
t338 = t250 * t394;
t118 = t196 * t414 + t250 * t418 - t338 * t382;
t342 = t245 * t395;
t165 = t253 * t342 + t313 * t394;
t249 = sin(qJ(4));
t252 = cos(qJ(4));
t84 = t118 * t252 - t165 * t249;
t437 = t119 * t251 + t248 * t84;
t436 = -t119 * t248 + t251 * t84;
t339 = t250 * t395;
t272 = t245 * (-t244 * t339 + t414 * t246);
t179 = qJD(1) * t272;
t299 = qJD(3) * t320;
t432 = -t179 + t299;
t323 = t246 * t339;
t334 = t394 * t247;
t325 = t250 * t334;
t163 = (t414 * t244 + t323) * t245 + t325;
t431 = t163 * qJD(3);
t259 = qJD(1) * t431;
t359 = qJDD(1) * t245;
t346 = t244 * t359;
t321 = t395 * t414;
t383 = t245 * t246;
t294 = t321 * t383;
t427 = -t247 * t320 - t294;
t289 = t427 * qJDD(1) + t250 * t346;
t280 = qJDD(4) + t289;
t417 = t259 + t280;
t384 = t244 * t245;
t349 = qJD(1) * t384;
t372 = t427 * qJD(1);
t150 = t250 * t349 + t372;
t281 = qJD(4) + t150;
t85 = -t118 * t249 - t165 * t252;
t411 = pkin(9) * t252;
t367 = qJD(4) * t249;
t154 = t163 * qJD(1);
t100 = pkin(3) * t154 + pkin(9) * t150;
t231 = qJ(2) * t383;
t412 = pkin(1) * t247;
t355 = qJD(1) * t412;
t192 = qJD(1) * t231 + t244 * t355;
t268 = (t246 * t342 + t334) * pkin(8);
t139 = qJD(1) * t268 + t192;
t183 = (-t394 * pkin(8) * t244 - pkin(2) * t246 - pkin(1)) * t245;
t171 = qJD(1) * t183 + qJD(2);
t229 = t246 * t355;
t273 = t247 * pkin(2) + (-t395 * pkin(8) - qJ(2)) * t384;
t145 = qJD(1) * t273 + t229;
t304 = t145 * t321;
t72 = -t250 * t139 + t171 * t320 + t304;
t49 = t249 * t100 + t252 * t72;
t433 = pkin(9) * t367 + t49;
t340 = t247 * t395;
t341 = t245 * t394;
t195 = t246 * t341 - t340;
t187 = qJDD(1) * t195 - qJDD(3);
t361 = qJD(1) * qJD(2);
t347 = t245 * t361;
t357 = qJDD(1) * t247;
t354 = pkin(1) * t357;
t175 = qJDD(1) * t231 + t244 * t354 + t246 * t347;
t127 = qJDD(1) * t268 + t175;
t227 = t246 * t354;
t327 = t244 * t347;
t130 = qJDD(1) * t273 + t227 - t327;
t167 = qJDD(1) * t183 + qJDD(2);
t318 = qJD(3) * t338;
t337 = t395 * t145;
t351 = t414 * t139;
t369 = qJD(3) * t250;
t37 = -qJD(3) * t351 - t250 * t127 + t130 * t321 + t167 * t320 - t171 * t318 - t337 * t369;
t34 = pkin(3) * t187 - t37;
t332 = qJD(1) * t394;
t319 = t245 * t332;
t362 = -t246 * t319 + qJD(3);
t287 = -qJD(1) * t340 - t362;
t182 = t252 * t287;
t353 = t250 * t384;
t328 = qJD(3) * t353;
t98 = qJD(1) * t328 + t372 * qJD(3) - qJDD(1) * t325 - t323 * t359 - t414 * t346;
t53 = qJD(4) * t182 + t154 * t367 + t249 * t187 + t252 * t98;
t333 = t252 * t187 - t249 * t98;
t109 = t252 * t154 - t249 * t287;
t368 = qJD(4) * t109;
t54 = t333 + t368;
t10 = pkin(4) * t54 + pkin(10) * t53 + t34;
t105 = -t394 * t145 + t395 * t171;
t58 = t150 * pkin(3) - t154 * pkin(9) + t105;
t73 = t351 + (t394 * t171 + t337) * t250;
t60 = -pkin(9) * t287 + t73;
t29 = t249 * t58 + t252 * t60;
t26 = pkin(10) * t281 + t29;
t107 = t154 * t249 + t182;
t59 = pkin(3) * t287 - t72;
t35 = t107 * pkin(4) - t109 * pkin(10) + t59;
t309 = t248 * t26 - t251 * t35;
t277 = -qJD(3) * t304 - t414 * t127 - t130 * t339 + t139 * t369 - t167 * t338 - t171 * t299;
t33 = -pkin(9) * t187 - t277;
t366 = qJD(4) * t252;
t91 = -t394 * t130 + t395 * t167;
t99 = t259 + t289;
t45 = t99 * pkin(3) + t98 * pkin(9) + t91;
t296 = -t249 * t45 - t252 * t33 - t58 * t366 + t60 * t367;
t5 = t417 * pkin(10) - t296;
t1 = -t309 * qJD(5) + t248 * t10 + t251 * t5;
t106 = qJD(5) + t107;
t429 = t309 * t106 + t1;
t200 = t249 * t338 - t252 * t395;
t301 = t244 * t319;
t374 = -qJD(4) * t200 - t249 * t301 + t432 * t252;
t288 = t247 * t350 + t378;
t428 = t288 * t394 + t413 * t342;
t28 = -t249 * t60 + t252 * t58;
t426 = -t28 * t281 - t296;
t271 = t245 * (t244 * t321 + t246 * t250);
t178 = qJD(1) * t271;
t415 = t318 - t178;
t425 = pkin(10) * t154 + t433;
t424 = -qJD(5) * t411 - t73 + t281 * (pkin(4) * t249 - pkin(10) * t252);
t199 = t244 * t412 + t231;
t156 = t268 + t199;
t235 = t246 * t412;
t164 = t235 + t273;
t78 = t414 * t156 + (t395 * t164 + t394 * t183) * t250;
t12 = t248 * t35 + t251 * t26;
t2 = -qJD(5) * t12 + t251 * t10 - t248 * t5;
t423 = -t12 * t106 - t2;
t421 = t107 * t281;
t420 = t109 * t281;
t241 = t244 ^ 2;
t242 = t245 ^ 2;
t243 = t246 ^ 2;
t419 = t242 * (t241 + t243);
t77 = -t250 * t156 + t164 * t321 + t183 * t320;
t201 = t249 * t395 + t252 * t338;
t373 = qJD(4) * t201 + t432 * t249 + t252 * t301;
t416 = t288 * t395 - t413 * t341;
t255 = qJDD(4) + t99;
t76 = t251 * t109 + t248 * t281;
t393 = qJD(5) * t76;
t24 = -t248 * t53 - t251 * t255 + t393;
t343 = t249 * t33 - t252 * t45;
t8 = -qJD(4) * t29 - t343;
t110 = -t394 * t164 + t395 * t183;
t162 = t353 + t427;
t160 = t162 * pkin(3);
t67 = -t163 * pkin(9) + t110 + t160;
t71 = -pkin(9) * t195 + t78;
t408 = t249 * t67 + t252 * t71;
t62 = qJD(2) * t272 + t77 * qJD(3);
t152 = -qJD(3) * t294 - t247 * t299 + t328;
t300 = t244 * qJD(2) * t341;
t90 = pkin(3) * t431 + t152 * pkin(9) + t300;
t18 = -qJD(4) * t408 - t249 * t62 + t252 * t90;
t410 = g(2) * t253;
t276 = t251 * t281;
t74 = t109 * t248 - t276;
t409 = t76 * t74;
t407 = t106 * t74;
t365 = qJD(5) * t248;
t23 = -qJD(5) * t276 + t109 * t365 - t248 * t255 + t251 * t53;
t404 = t23 * t248;
t403 = t24 * t251;
t52 = qJDD(5) + t54;
t402 = t248 * t52;
t401 = t248 * t74;
t400 = t251 * t52;
t399 = t59 * t150;
t398 = t76 * t106;
t221 = -pkin(4) * t252 - pkin(10) * t249 - pkin(3);
t364 = qJD(5) * t251;
t397 = t221 * t364 + t248 * t424 - t251 * t425;
t396 = -t221 * t365 + t248 * t425 + t251 * t424;
t392 = t109 * t107;
t391 = t119 * t249;
t390 = t119 * t252;
t197 = -t247 * t238 + t377;
t121 = t197 * t250 + t416 * t414;
t389 = t121 * t249;
t388 = t121 * t252;
t387 = t154 * t150;
t386 = t162 * t249;
t385 = t162 * t252;
t254 = qJD(1) ^ 2;
t381 = t247 * t254;
t380 = t248 * t252;
t379 = t251 * t252;
t279 = -t251 * t201 + t248 * t320;
t376 = -t279 * qJD(5) + t374 * t248 - t415 * t251;
t172 = -t248 * t201 - t251 * t320;
t375 = -t172 * qJD(5) - t415 * t248 - t374 * t251;
t352 = t245 * t413;
t371 = t253 * pkin(1) + qJ(2) * t352;
t360 = qJDD(1) * t242;
t358 = qJDD(1) * t246;
t356 = g(1) * t413;
t348 = qJD(4) + t372;
t345 = g(2) * t382 - g(3) * t247;
t331 = t281 * t249;
t330 = t106 * t251;
t329 = 0.2e1 * t245 * t357;
t122 = t197 * t414 - t416 * t250;
t87 = t122 * t249 - t252 * t428;
t326 = g(1) * t85 + g(2) * t87;
t322 = -t413 * pkin(1) + qJ(2) * t382;
t316 = g(1) * t119 + g(2) * t121;
t92 = -t150 * t380 - t251 * t154;
t315 = -t248 * t366 + t92;
t93 = -t150 * t379 + t154 * t248;
t314 = t251 * t366 - t93;
t310 = -t12 * t248 + t251 * t309;
t32 = pkin(10) * t162 + t408;
t115 = t163 * t249 + t195 * t252;
t116 = t163 * t252 - t195 * t249;
t70 = t195 * pkin(3) - t77;
t44 = t115 * pkin(4) - t116 * pkin(10) + t70;
t14 = t248 * t44 + t251 * t32;
t13 = -t248 * t32 + t251 * t44;
t38 = -t249 * t71 + t252 * t67;
t48 = t100 * t252 - t249 * t72;
t82 = t116 * t251 + t162 * t248;
t81 = t116 * t248 - t162 * t251;
t306 = (-qJ(2) * t349 + t229) * t244 - t192 * t246;
t298 = g(1) * t253 + g(2) * t413;
t230 = -pkin(1) * t359 + qJDD(2);
t297 = pkin(1) * t360 - t230 * t245;
t17 = t249 * t90 + t252 * t62 + t67 * t366 - t71 * t367;
t25 = -pkin(4) * t281 - t28;
t295 = -pkin(10) * t52 + t106 * t25;
t293 = g(1) * t87 - g(2) * t85 + g(3) * t115;
t88 = t122 * t252 + t249 * t428;
t292 = -g(1) * t88 - g(2) * t84 - g(3) * t116;
t291 = g(1) * t121 - g(2) * t119 + g(3) * t162;
t290 = g(1) * t122 + g(2) * t118 + g(3) * t163;
t6 = -t417 * pkin(4) - t8;
t286 = t293 - t6;
t282 = -g(1) * t352 + t345;
t275 = qJD(4) * t281;
t267 = pkin(10) * qJD(5) * t106 - t286;
t265 = -t196 * pkin(2) + t165 * pkin(8) + t322;
t264 = t197 * pkin(2) + t428 * pkin(8) + t371;
t260 = -pkin(3) * t118 + t119 * pkin(9) + t265;
t258 = t122 * pkin(3) + t121 * pkin(9) + t264;
t257 = t150 * t281 + t275;
t63 = qJD(2) * t271 + qJD(3) * t78;
t198 = -qJ(2) * t384 + t235;
t189 = pkin(9) * t379 + t221 * t248;
t188 = -pkin(9) * t380 + t221 * t251;
t174 = t227 + (-qJ(2) * qJDD(1) - t361) * t384;
t113 = t121 * pkin(3);
t111 = t119 * pkin(3);
t80 = -qJD(4) * t115 - t152 * t252;
t79 = qJD(4) * t116 - t152 * t249;
t65 = pkin(4) * t109 + pkin(10) * t107;
t56 = t121 * t248 + t251 * t88;
t55 = t121 * t251 - t248 * t88;
t47 = -qJD(5) * t81 + t248 * t431 + t80 * t251;
t46 = qJD(5) * t82 + t80 * t248 - t251 * t431;
t40 = -pkin(4) * t154 - t48;
t31 = -pkin(4) * t162 - t38;
t27 = t79 * pkin(4) - t80 * pkin(10) + t63;
t22 = t248 * t65 + t251 * t28;
t21 = -t248 * t28 + t251 * t65;
t16 = -pkin(4) * t431 - t18;
t15 = pkin(10) * t431 + t17;
t4 = -qJD(5) * t14 - t248 * t15 + t251 * t27;
t3 = qJD(5) * t13 + t251 * t15 + t248 * t27;
t7 = [0, 0, 0, 0, 0, qJDD(1), t356 - t410, t298, 0, 0, t241 * t360, 0.2e1 * t242 * t244 * t358, t244 * t329, t243 * t360, t246 * t329, t247 ^ 2 * qJDD(1), g(1) * t196 - g(2) * t197 + t297 * t246 + (qJDD(1) * t198 + t174 - t327) * t247, -g(1) * t238 + (-t297 + t410) * t244 + (-t199 * qJDD(1) - t175 + (t298 - t347) * t246) * t247, t361 * t419 + (-t174 * t244 + t175 * t246 + (-t198 * t244 + t199 * t246) * qJDD(1) - t298) * t245, t175 * t199 + t174 * t198 - g(1) * t322 - g(2) * t371 + (-t230 * pkin(1) - qJD(2) * t306) * t245, -t152 * t154 - t163 * t98, t150 * t152 - t154 * t431 + t162 * t98 - t163 * t99, t152 * t287 - t163 * t187 + t98 * t195, t150 * t431 + t162 * t99, t162 * t187 + t99 * t195 + t287 * t431, t187 * t195, g(1) * t118 - g(2) * t122 + t105 * t431 + t110 * t99 + t150 * t300 + t91 * t162 - t77 * t187 - t37 * t195 + t287 * t63, -t105 * t152 - t110 * t98 + t154 * t300 + t91 * t163 + t78 * t187 - t195 * t277 + t287 * t62 + t316, -g(1) * t165 - g(2) * t428 - t62 * t150 + t72 * t152 + t63 * t154 + t162 * t277 - t37 * t163 - t431 * t73 + t77 * t98 - t78 * t99, -g(1) * t265 - g(2) * t264 + t105 * t300 + t91 * t110 - t277 * t78 + t37 * t77 + t73 * t62 - t72 * t63, t109 * t80 - t116 * t53, -t107 * t80 - t109 * t79 + t115 * t53 - t116 * t54, t80 * t348 + t116 * t280 - t53 * t162 + t109 * t431 + (t116 * t431 + t353 * t80) * qJD(1), t107 * t79 + t115 * t54, -t79 * t348 - t115 * t280 - t54 * t162 - t107 * t431 + (-t115 * t431 - t353 * t79) * qJD(1), t280 * t162 + (t348 + (t162 + t353) * qJD(1)) * t431, g(1) * t84 - g(2) * t88 + t63 * t107 + t34 * t115 + t8 * t162 + t18 * t281 + t255 * t38 + t28 * t431 + t70 * t54 + t59 * t79, t63 * t109 + t34 * t116 + t162 * t296 - t17 * t281 - t29 * t431 - t408 * t417 - t70 * t53 + t59 * t80 + t326, -t107 * t17 - t109 * t18 + t115 * t296 - t116 * t8 - t28 * t80 - t29 * t79 + t38 * t53 - t408 * t54 - t316, -g(1) * t260 - g(2) * t258 + t29 * t17 + t28 * t18 - t296 * t408 + t34 * t70 + t8 * t38 + t59 * t63, -t23 * t82 + t47 * t76, t23 * t81 - t24 * t82 - t46 * t76 - t47 * t74, t106 * t47 - t115 * t23 + t52 * t82 + t76 * t79, t24 * t81 + t46 * t74, -t106 * t46 - t115 * t24 - t52 * t81 - t74 * t79, t106 * t79 + t115 * t52, g(1) * t436 - g(2) * t56 + t4 * t106 + t2 * t115 + t13 * t52 + t16 * t74 + t31 * t24 + t25 * t46 - t309 * t79 + t6 * t81, -g(1) * t437 - g(2) * t55 - t1 * t115 - t3 * t106 - t12 * t79 - t14 * t52 + t16 * t76 - t31 * t23 + t25 * t47 + t6 * t82, -t1 * t81 - t12 * t46 + t13 * t23 - t14 * t24 - t2 * t82 - t3 * t74 + t309 * t47 - t4 * t76 - t326, t1 * t14 + t12 * t3 + t2 * t13 - t309 * t4 + t6 * t31 + t25 * t16 - g(1) * (-pkin(4) * t84 + t85 * pkin(10) + t260) - g(2) * (t88 * pkin(4) + t87 * pkin(10) + t258); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t244 * t381 - t358) * t245, (qJDD(1) * t244 + t246 * t381) * t245, -t254 * t419, qJDD(2) + (-pkin(1) * qJDD(1) + t306 * qJD(1) - t356) * t245 + t345, 0, 0, 0, 0, 0, 0, -t150 * t301 - t187 * t320 + t415 * t287 + t395 * t99, -t154 * t301 + t187 * t338 + t432 * t287 - t395 * t98, t98 * t320 - t99 * t338 + t179 * t150 - t178 * t154 + (-t150 * t320 + t154 * t338) * qJD(3), t37 * t320 - t277 * t338 + t91 * t395 + t72 * t178 - t73 * t179 + (-t105 * t244 * t332 - t356) * t245 + (t320 * t73 - t338 * t72) * qJD(3) + t345, 0, 0, 0, 0, 0, 0, -t200 * t280 - t54 * t320 - t178 * t107 + (t107 * t338 - t154 * t200) * qJD(3) - t373 * t281, -t201 * t280 + t53 * t320 - t178 * t109 + (t109 * t338 - t201 * t154) * qJD(3) - t374 * t281, -t107 * t374 + t109 * t373 - t200 * t53 - t201 * t54, -t8 * t200 - t201 * t296 - t373 * t28 + t374 * t29 - t34 * t320 + t415 * t59 + t282, 0, 0, 0, 0, 0, 0, -t106 * t376 + t172 * t52 + t200 * t24 + t373 * t74, t106 * t375 - t200 * t23 + t279 * t52 + t373 * t76, t172 * t23 + t24 * t279 + t375 * t74 + t376 * t76, -t1 * t279 - t12 * t375 + t2 * t172 + t6 * t200 + t25 * t373 + t309 * t376 + t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t387, -t150 ^ 2 + t154 ^ 2, -t150 * t287 - t98, -t387, t154 * t362 + (t154 * t340 - t431) * qJD(1) - t289, -t187, -t105 * t154 - t287 * t73 + t291 + t37, t105 * t150 - t287 * t72 + t277 + t290, 0, 0, -t53 * t249 + t252 * t420, (-t53 - t421) * t252 + (-t54 - t420) * t249, -t109 * t154 + t417 * t249 + t252 * t257, t107 * t331 - t54 * t252, t107 * t154 - t249 * t257 + t417 * t252, -t281 * t154, -pkin(3) * t54 + g(1) * t388 - g(2) * t390 + g(3) * t385 - t73 * t107 - t28 * t154 - t34 * t252 - t275 * t411 - t281 * t48 + t367 * t59 + (-pkin(9) * t255 + t399) * t249, pkin(3) * t53 - g(1) * t389 + g(2) * t391 - g(3) * t386 - t73 * t109 + t29 * t154 + t34 * t249 + t252 * t399 + t433 * t281 + t366 * t59 - t411 * t417, t107 * t49 + t109 * t48 + ((-t54 + t368) * pkin(9) + t426) * t252 + (-t8 - t281 * t29 + (qJD(4) * t107 - t53) * pkin(9)) * t249 - t290, -t34 * pkin(3) + g(1) * t113 - g(2) * t111 + g(3) * t160 - t28 * t48 - t29 * t49 - t59 * t73 + (-t8 * t249 - t296 * t252 + (-t249 * t29 - t252 * t28) * qJD(4) - t290) * pkin(9), -t23 * t249 * t251 + (-t249 * t365 + t314) * t76, t74 * t93 + t76 * t92 + (-t248 * t76 - t251 * t74) * t366 + (t404 - t403 + (-t251 * t76 + t401) * qJD(5)) * t249, t23 * t252 + t314 * t106 + (-t106 * t365 + t281 * t76 + t400) * t249, t24 * t248 * t249 + (t249 * t364 - t315) * t74, t24 * t252 + t315 * t106 + (-t106 * t364 - t281 * t74 - t402) * t249, t106 * t331 - t52 * t252, t188 * t52 - t25 * t92 - t40 * t74 + t396 * t106 - t290 * t248 + (-t2 + (pkin(9) * t74 + t248 * t25) * qJD(4) + t291 * t251) * t252 + (pkin(9) * t24 + t6 * t248 + t25 * t364 - t281 * t309) * t249, -t189 * t52 - t25 * t93 - t40 * t76 - t397 * t106 - t290 * t251 + (t1 + (pkin(9) * t76 + t25 * t251) * qJD(4) - t291 * t248) * t252 + (-pkin(9) * t23 - t12 * t281 - t25 * t365 + t6 * t251) * t249, -t309 * t93 + t12 * t92 + t188 * t23 - t189 * t24 - t396 * t76 - t397 * t74 + t310 * t366 + (-t1 * t248 - t2 * t251 + (-t12 * t251 - t248 * t309) * qJD(5) + t291) * t249, t1 * t189 + t2 * t188 - t25 * t40 - g(1) * (-pkin(4) * t388 - pkin(10) * t389 - t113) - g(2) * (pkin(4) * t390 + pkin(10) * t391 + t111) - g(3) * (-pkin(4) * t385 - pkin(10) * t386 - t160) + t397 * t12 - t396 * t309 + (t249 * t6 + t25 * t366 - t290) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392, -t107 ^ 2 + t109 ^ 2, -t53 + t421, -t392, t109 * t150 - t333, t255, -t59 * t109 + t150 * t29 + t293 - t343, t59 * t107 - t292 - t426, 0, 0, t330 * t76 - t404, (-t23 - t407) * t251 + (-t24 - t398) * t248, t106 * t330 - t76 * t109 + t402, t106 * t401 - t403, -t106 ^ 2 * t248 + t74 * t109 + t400, -t106 * t109, -pkin(4) * t24 - t106 * t21 + t109 * t309 + t248 * t295 - t251 * t267 - t29 * t74, pkin(4) * t23 + t106 * t22 + t109 * t12 + t248 * t267 + t251 * t295 - t29 * t76, t21 * t76 + t22 * t74 + ((-t24 + t393) * pkin(10) + t429) * t251 + ((qJD(5) * t74 - t23) * pkin(10) + t423) * t248 + t292, t309 * t21 - t12 * t22 - t25 * t29 + t286 * pkin(4) + (qJD(5) * t310 + t1 * t251 - t2 * t248 + t292) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t409, -t74 ^ 2 + t76 ^ 2, -t23 + t407, -t409, t398 - t24, t52, -g(1) * t55 + g(2) * t437 + g(3) * t81 - t25 * t76 - t423, g(1) * t56 + g(2) * t436 + g(3) * t82 + t25 * t74 - t429, 0, 0;];
tau_reg = t7;
