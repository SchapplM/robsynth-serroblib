% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPPRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:32:00
% EndTime: 2019-03-09 08:32:10
% DurationCPUTime: 4.81s
% Computational Cost: add. (6138->352), mult. (11172->474), div. (0->0), fcn. (12241->6), ass. (0->285)
t214 = cos(qJ(5));
t215 = cos(qJ(2));
t354 = cos(pkin(9));
t202 = t354 * t215;
t211 = sin(pkin(9));
t213 = sin(qJ(2));
t325 = t211 * t213;
t184 = -t202 + t325;
t383 = t184 * pkin(4);
t378 = -qJ(3) - pkin(7);
t194 = t378 * t215;
t189 = t354 * t194;
t193 = t378 * t213;
t327 = t211 * t193;
t401 = -t189 + t327;
t403 = t401 - t383;
t409 = t214 * t403;
t270 = t409 / 0.2e1;
t212 = sin(qJ(5));
t410 = t212 * t403;
t411 = -t410 / 0.2e1;
t181 = t184 ^ 2;
t186 = t211 * t215 + t354 * t213;
t392 = t186 ^ 2;
t399 = -t392 - t181;
t408 = qJD(3) * t399;
t407 = t399 * qJD(1);
t188 = t354 * t193;
t326 = t211 * t194;
t130 = -t188 - t326;
t380 = t186 * pkin(4);
t94 = t130 + t380;
t406 = t212 * t94;
t88 = t214 * t94;
t304 = qJD(5) * t212;
t278 = pkin(5) * t304;
t379 = t186 * pkin(5);
t269 = -t215 * pkin(2) - pkin(1);
t236 = -t186 * qJ(4) + t269;
t389 = pkin(3) + pkin(8);
t75 = t389 * t184 + t236;
t253 = qJ(6) * t184 + t75;
t53 = -t253 * t212 + t88;
t42 = t53 + t379;
t273 = -t53 / 0.2e1 + t42 / 0.2e1;
t274 = t379 / 0.2e1;
t235 = t274 + t273;
t225 = t235 * t212;
t393 = t225 * qJD(1);
t405 = t278 + t393;
t108 = t212 * t186;
t302 = t108 * qJD(1);
t404 = t302 + t304;
t111 = t214 * t184;
t252 = 0.2e1 * t212 * t111;
t366 = t42 * t212;
t54 = t253 * t214 + t406;
t245 = -t54 * t214 + t366;
t23 = t245 * t186;
t203 = -t354 * pkin(2) - pkin(3);
t200 = -pkin(8) + t203;
t402 = -qJ(6) + t200;
t390 = t53 / 0.2e1;
t400 = t390 + t274;
t242 = t130 * t186 - t184 * t401;
t398 = qJD(1) * t242;
t396 = qJD(3) * t242;
t395 = qJD(5) * t225;
t175 = t325 / 0.2e1 - t202 / 0.2e1;
t353 = qJ(6) * t186;
t382 = t184 * pkin(5);
t208 = t213 * pkin(2);
t334 = t184 * qJ(4);
t251 = t208 + t334;
t76 = t389 * t186 + t251;
t43 = -t382 + t409 + (-t76 - t353) * t212;
t391 = t43 / 0.2e1;
t388 = -t184 / 0.2e1;
t387 = -t186 / 0.2e1;
t386 = t186 / 0.2e1;
t248 = -t189 / 0.2e1;
t210 = t214 ^ 2;
t385 = t210 / 0.2e1;
t384 = pkin(5) * t214;
t381 = t186 * pkin(3);
t362 = t54 * t212;
t365 = t42 * t214;
t377 = -t365 / 0.2e1 - t362 / 0.2e1;
t376 = -t42 + t53;
t74 = t214 * t76;
t375 = -t74 - t410;
t374 = pkin(5) * qJD(5);
t333 = t184 * t212;
t281 = pkin(5) * t333;
t268 = -pkin(4) - t384;
t71 = t268 * t184 + t401;
t7 = t71 * t281 + t376 * t54;
t373 = qJD(1) * t7;
t8 = t376 * t111;
t372 = qJD(1) * t8;
t371 = qJD(2) * pkin(2);
t370 = t212 * t76;
t55 = t214 * t353 - t375;
t70 = t268 * t186 - t130;
t3 = t42 * t43 + t54 * t55 + t70 * t71;
t368 = t3 * qJD(1);
t364 = t43 * t212;
t363 = t43 * t214;
t361 = t55 * t212;
t360 = t55 * t214;
t6 = -t23 + (t360 - t364) * t184;
t359 = t6 * qJD(1);
t17 = -t71 * t184 + (t362 + t365) * t186;
t352 = qJD(1) * t17;
t201 = pkin(2) * t211 + qJ(4);
t190 = pkin(5) * t212 + t201;
t256 = t190 * t388;
t171 = t402 * t214;
t335 = t171 * t214;
t170 = t402 * t212;
t338 = t170 * t212;
t216 = (t338 / 0.2e1 + t335 / 0.2e1) * t186 + t256;
t234 = t364 / 0.2e1 - t360 / 0.2e1;
t19 = t216 + t234;
t351 = qJD(1) * t19;
t22 = t245 * t184;
t350 = qJD(1) * t22;
t349 = qJD(1) * t23;
t56 = t212 * t75 - t88;
t29 = t111 * t403 + t186 * t56;
t348 = qJD(1) * t29;
t57 = t214 * t75 + t406;
t30 = -t186 * t57 + t333 * t403;
t347 = qJD(1) * t30;
t107 = t184 * pkin(3) + t236;
t113 = t251 + t381;
t58 = -t107 * t186 - t113 * t184;
t346 = qJD(1) * t58;
t59 = t107 * t184 - t113 * t186;
t345 = qJD(1) * t59;
t68 = t399 * t212;
t342 = qJD(1) * t68;
t282 = t392 - t181;
t79 = t282 * t212;
t341 = qJD(1) * t79;
t80 = t282 * t214;
t340 = qJD(1) * t80;
t10 = t235 * t214;
t339 = t10 * qJD(1);
t337 = t170 * t214;
t336 = t171 * t212;
t332 = t186 * t201;
t20 = -t370 * t186 + (t56 + t88) * t184;
t331 = t20 * qJD(1);
t330 = t201 * t184;
t209 = t212 ^ 2;
t329 = t209 * t184;
t165 = t209 * t186;
t21 = (t410 + t375) * t186 + (t57 - t406) * t184;
t328 = t21 * qJD(1);
t166 = t210 * t186;
t28 = t107 * t113;
t324 = t28 * qJD(1);
t52 = t269 * t208;
t321 = t52 * qJD(1);
t204 = t208 / 0.2e1;
t63 = t204 + (pkin(3) / 0.2e1 - t203 / 0.2e1) * t186 + (qJ(4) / 0.2e1 + t201 / 0.2e1) * t184;
t320 = t63 * qJD(1);
t254 = t385 + t209 / 0.2e1;
t65 = t254 * t186 + t165 / 0.2e1 + t166 / 0.2e1;
t319 = t65 * qJD(1);
t313 = t209 + t210;
t69 = t313 * t186 * t184;
t318 = t69 * qJD(1);
t77 = -t254 * t184 - t175;
t317 = t77 * qJD(1);
t81 = t399 * t214;
t316 = t81 * qJD(1);
t219 = (t211 * t388 + t354 * t387) * pkin(2);
t90 = -t208 / 0.2e1 + t219;
t315 = t90 * qJD(1);
t197 = t209 - t210;
t312 = qJD(1) * t212;
t311 = qJD(1) * t215;
t310 = qJD(2) * t111;
t309 = qJD(3) * t184;
t308 = qJD(4) * t186;
t307 = qJD(4) * t212;
t306 = qJD(4) * t214;
t305 = qJD(5) * t108;
t303 = qJD(5) * t214;
t301 = t111 * qJD(1);
t114 = -t165 - t166;
t300 = t114 * qJD(1);
t115 = t313 * t181;
t299 = t115 * qJD(1);
t296 = t175 * qJD(1);
t295 = t392 * qJD(1);
t294 = t184 * qJD(1);
t293 = t184 * qJD(2);
t292 = t184 * qJD(4);
t291 = t186 * qJD(1);
t290 = t186 * qJD(2);
t172 = t186 * qJD(3);
t192 = -0.1e1 / 0.2e1 - t254;
t289 = t192 * qJD(2);
t288 = t313 * qJD(2);
t199 = -t213 ^ 2 + t215 ^ 2;
t287 = t199 * qJD(1);
t286 = t212 * qJD(2);
t285 = t213 * qJD(2);
t284 = t214 * qJD(2);
t283 = t215 * qJD(2);
t280 = pkin(1) * t213 * qJD(1);
t279 = pkin(1) * t311;
t277 = pkin(5) * t303;
t276 = t384 / 0.2e1;
t275 = -t379 / 0.2e1;
t272 = -t74 / 0.2e1 + t411;
t267 = t107 * t291;
t266 = t392 * t312;
t265 = t184 * t286;
t264 = t212 * t284;
t263 = t186 * t304;
t262 = t186 * t303;
t261 = t184 * t291;
t260 = t184 * t290;
t259 = t212 * t294;
t258 = t213 * t311;
t257 = t212 * t303;
t162 = t214 * t291;
t250 = qJD(5) + t291;
t249 = pkin(4) / 0.2e1 + t276;
t247 = qJD(2) * t252;
t246 = -t188 / 0.2e1 - t326 / 0.2e1;
t1 = t273 * t170 + (t391 + t212 * t256 - t71 * t214 / 0.2e1) * pkin(5);
t62 = t190 * t384;
t244 = -qJD(1) * t1 + qJD(2) * t62;
t241 = t184 * t200 + t332;
t100 = t335 + t338;
t233 = t337 / 0.2e1 - t336 / 0.2e1;
t218 = t184 * t233 + t377;
t14 = t249 * t186 + t218 + t246;
t240 = qJD(1) * t14 - qJD(2) * t100;
t221 = t327 / 0.2e1 + t248;
t16 = -t383 / 0.2e1 + (t171 * t386 - t55 / 0.2e1) * t212 + (t170 * t387 - t382 / 0.2e1 - t43 / 0.2e1) * t214 + t221;
t239 = qJD(1) * t16 + qJD(2) * t190;
t238 = t250 * t212;
t128 = t248 + t189 / 0.2e1;
t237 = qJD(1) * t128 + qJD(2) * t201;
t232 = t200 * t387 + t330 / 0.2e1;
t217 = t232 * t214 + t411;
t26 = t217 - t272;
t231 = -qJD(1) * t26 + t201 * t286;
t24 = t270 - t409 / 0.2e1 + (t76 / 0.2e1 + t232) * t212;
t230 = -qJD(1) * t24 - t201 * t284;
t229 = t184 * t238;
t157 = t184 * t385;
t110 = t157 - t329 / 0.2e1;
t228 = -qJD(1) * t110 + t264;
t227 = t259 + t284;
t226 = qJD(5) * t175 + t261;
t224 = t181 * t214 * t312 + qJD(2) * t110;
t116 = t197 * t181;
t223 = qJD(1) * t116 + t247;
t222 = qJD(1) * t252 - qJD(2) * t197;
t191 = 0.1e1 / 0.2e1 - t254;
t163 = t175 * qJD(2);
t139 = -t162 - t303;
t134 = t227 * pkin(5);
t105 = t110 * qJD(5);
t97 = 0.2e1 * t248 + t327;
t89 = t204 + t219;
t78 = t157 + t329 / 0.2e1 - t175;
t64 = -t330 / 0.2e1 + t203 * t386 + t204 + t334 / 0.2e1 + t381 / 0.2e1;
t27 = t217 + t272;
t25 = 0.2e1 * t270 - t370 / 0.2e1 + t232 * t212;
t18 = t216 - t234;
t15 = t361 / 0.2e1 + t363 / 0.2e1 - t233 * t186 - t249 * t184 + t221;
t13 = -t380 / 0.2e1 + t214 * t275 + t218 - t246;
t12 = t362 / 0.2e1 + t377 + t400 * t214;
t11 = -t366 / 0.2e1 + t400 * t212;
t5 = (t273 + t275) * t212;
t2 = t190 * t281 / 0.2e1 + t71 * t276 + pkin(5) * t391 + (t390 - t42 / 0.2e1) * t170;
t4 = [0, 0, 0, t213 * t283, t199 * qJD(2), 0, 0, 0, -pkin(1) * t285, -pkin(1) * t283, -t408, qJD(2) * t52 + t396, -t408, qJD(2) * t58 + t186 * t292, qJD(2) * t59 + qJD(4) * t392, qJD(2) * t28 - t107 * t308 + t396, t181 * t257 + t209 * t260, -qJD(5) * t116 + t186 * t247, qJD(2) * t79 + t184 * t262, qJD(2) * t80 - t184 * t263, -t260, qJD(2) * t20 - qJD(3) * t81 + qJD(5) * t30 + t392 * t307, qJD(2) * t21 + qJD(3) * t68 + qJD(5) * t29 + t306 * t392, qJD(2) * t6 - qJD(4) * t69 + qJD(5) * t8 + qJD(6) * t115, qJD(2) * t3 + qJD(3) * t17 + qJD(4) * t23 + qJD(5) * t7 - qJD(6) * t22; 0, 0, 0, t258, t287, t283, -t285, 0, -pkin(7) * t283 - t280, pkin(7) * t285 - t279 (t354 * t184 - t186 * t211) * t371, t321 + (-t130 * t211 - t354 * t401) * t371 + t89 * qJD(3) (-t184 * t203 - t332) * qJD(2) - t292, qJD(2) * t401 + t346, -qJD(2) * t130 + t345, t324 + (-t130 * t201 + t203 * t401) * qJD(2) + t64 * qJD(3) + t97 * qJD(4), t105 + (t209 * t294 + t264) * t186 (t166 - t165) * qJD(2) + (-qJD(5) + t291) * t252, -t184 * t284 + t341, t265 + t340, -t226, t331 + (-t214 * t241 - t406) * qJD(2) - t111 * qJD(4) + t25 * qJD(5), -t94 * t284 + t328 + t27 * qJD(5) + (qJD(2) * t241 + t292) * t212, t359 + (-t361 - t363 + (-t336 + t337) * t186) * qJD(2) + t5 * qJD(5), t368 + (t170 * t55 + t171 * t43 + t190 * t70) * qJD(2) + t18 * qJD(3) + t15 * qJD(4) + t2 * qJD(5) + t13 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t407, qJD(2) * t89 + t398, -t407, 0, 0, qJD(2) * t64 + t398, 0, 0, 0, 0, 0, -t316, t342, 0, qJD(2) * t18 + qJD(5) * t12 + qJD(6) * t78 + t352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, t261, t295, qJD(2) * t97 - t267, 0, 0, 0, 0, 0, t266 - t310, t214 * t295 + t265, -t318, qJD(2) * t15 + qJD(5) * t11 + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, -t223, t250 * t111, -t229, -t163, qJD(2) * t25 - qJD(5) * t57 + t347, qJD(2) * t27 + qJD(5) * t56 + t348, qJD(2) * t5 - t184 * t277 + t372, qJD(2) * t2 + qJD(3) * t12 + qJD(4) * t11 - t374 * t54 + t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, qJD(2) * t13 + qJD(3) * t78 - t350; 0, 0, 0, -t258, -t287, 0, 0, 0, t280, t279, 0, qJD(3) * t90 - t321, 0, t172 - t346, -t309 - t345, -qJD(3) * t63 + qJD(4) * t128 - t324, -t209 * t261 + t105, -0.2e1 * t214 * t229, -t263 - t341, -t262 - t340, t226, qJD(5) * t24 - t212 * t309 - t331, -qJD(3) * t111 + qJD(5) * t26 - t328, qJD(3) * t114 - t359 + t395, qJD(3) * t19 + qJD(4) * t16 - qJD(5) * t1 + qJD(6) * t14 - t368; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t201 * qJD(4), -t257, t197 * qJD(5), 0, 0, 0, t201 * t303 + t307, -t201 * t304 + t306, t313 * qJD(6), qJD(4) * t190 + qJD(5) * t62 - qJD(6) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, 0, t291, -t294, -t320, 0, 0, 0, 0, 0, -t259, -t301, t300, t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t237, 0, 0, 0, 0, 0, t286, t284, 0, qJD(6) * t191 + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, -t222, -t238, t139, t296, -t200 * t304 - t230, -t200 * t303 - t231, t405, -t170 * t374 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, qJD(4) * t191 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t407, -qJD(2) * t90 - t398, t407, -t290, t293, qJD(2) * t63 - t308 - t398, 0, 0, 0, 0, 0, -t262 + t265 + t316, t305 + t310 - t342, -t114 * qJD(2), -qJD(2) * t19 - qJD(4) * t65 - qJD(5) * t10 - qJD(6) * t77 - t352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t315, 0, -t291, t294, t320, 0, 0, 0, 0, 0, t259, t301, -t300, -t351; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, 0, 0, 0, 0, 0, 0, 0, 0, -t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t404, 0, -t277 - t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, -t295, -qJD(2) * t128 + t172 + t267, 0, 0, 0, 0, 0, -t266 - t305 (-qJD(5) * t186 - t295) * t214, t318, -qJD(2) * t16 + qJD(3) * t65 - t349 - t395; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t237, 0, 0, 0, 0, 0, -t286, -t284, 0, qJD(6) * t192 - t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, 0, 0, 0, 0, 0, 0, 0, 0, t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t404, -t250 * t214, 0, -t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, t223 (-t214 * t294 + t286) * t186, t227 * t186, -t163, -qJD(2) * t24 + qJD(4) * t108 + t172 * t214 - t347, -qJD(2) * t26 - qJD(3) * t108 + t186 * t306 - t348, -qJD(2) * t225 - t372, qJD(2) * t1 + qJD(3) * t10 + qJD(4) * t225 - qJD(6) * t281 - t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t222, t212 * t291, t162, -t296, t230, t231, -t393, -qJD(6) * t384 - t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, -t302, 0, t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t162, 0, t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, -qJD(2) * t14 + qJD(3) * t77 + t184 * t278 + t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, -qJD(4) * t192 - t240 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t4;
