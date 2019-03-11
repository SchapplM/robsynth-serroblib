% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:11
% EndTime: 2019-03-09 09:32:31
% DurationCPUTime: 8.94s
% Computational Cost: add. (4513->567), mult. (11087->720), div. (0->0), fcn. (8317->10), ass. (0->290)
t203 = sin(pkin(6));
t213 = cos(qJ(2));
t311 = qJD(1) * qJD(2);
t295 = t213 * t311;
t209 = sin(qJ(2));
t309 = qJDD(1) * t209;
t393 = t295 + t309;
t391 = t393 * t203;
t204 = cos(pkin(6));
t321 = qJD(1) * t204;
t183 = qJD(2) + t321;
t208 = sin(qJ(5));
t212 = cos(qJ(5));
t322 = qJD(1) * t203;
t301 = t209 * t322;
t106 = t183 * t212 + t208 * t301;
t300 = t213 * t322;
t146 = qJD(5) + t300;
t356 = t106 * t146;
t310 = qJDD(1) * t204;
t181 = qJDD(2) + t310;
t40 = qJD(5) * t106 + t208 * t181 - t212 * t391;
t392 = -t40 - t356;
t308 = qJDD(1) * t213;
t179 = t203 * t308;
t296 = t209 * t311;
t276 = t203 * t296;
t390 = -t179 + t276;
t306 = pkin(1) * t321;
t329 = -pkin(8) * t301 + t213 * t306;
t312 = qJD(3) - t329;
t139 = t212 * t301;
t104 = t183 * t208 - t139;
t100 = qJD(6) + t104;
t200 = t203 ^ 2;
t389 = 0.2e1 * t200;
t207 = sin(qJ(6));
t211 = cos(qJ(6));
t63 = t106 * t207 - t211 * t146;
t388 = t146 * t63;
t357 = t104 * t146;
t316 = qJD(5) * t208;
t39 = qJD(5) * t139 + t212 * t181 - t183 * t316 + t208 * t391;
t387 = -t39 + t357;
t115 = pkin(8) * t300 + t209 * t306;
t94 = pkin(3) * t300 + t115;
t333 = qJD(4) + t94;
t148 = t181 * qJ(3);
t150 = t183 * qJD(3);
t386 = -t148 - t150;
t169 = t183 * pkin(2);
t385 = -t183 * qJ(4) - t169;
t180 = t183 ^ 2;
t201 = t209 ^ 2;
t351 = t200 * qJD(1) ^ 2;
t384 = -t201 * t351 - t180;
t214 = cos(qJ(1));
t337 = t213 * t214;
t210 = sin(qJ(1));
t341 = t209 * t210;
t123 = -t204 * t337 + t341;
t339 = t210 * t213;
t340 = t209 * t214;
t125 = t204 * t339 + t340;
t347 = t203 * t213;
t232 = -g(1) * t125 - g(2) * t123 + g(3) * t347;
t383 = pkin(3) * t179 + qJDD(4);
t154 = t181 * pkin(2);
t382 = t154 - qJDD(3);
t372 = pkin(1) * t204;
t305 = qJD(2) * t372;
t281 = qJD(1) * t305;
t304 = pkin(1) * t310;
t279 = pkin(8) * t390 - t209 * t304 - t213 * t281;
t41 = t279 + t386;
t231 = -t41 + t383;
t240 = -t296 + t308;
t17 = -pkin(9) * t181 + (-pkin(3) * t296 + pkin(4) * t240) * t203 + t231;
t156 = t183 * qJ(3);
t334 = pkin(4) * t300 + t333;
t44 = -pkin(9) * t183 + t156 + t334;
t205 = qJ(3) - pkin(9);
t206 = pkin(2) + qJ(4);
t275 = -t206 * t213 - pkin(1);
t230 = -t205 * t209 + t275;
t56 = t230 * t322;
t21 = t208 * t44 + t212 * t56;
t318 = qJD(3) * t209;
t224 = -t318 + (-qJD(2) * t205 - qJD(4)) * t213;
t140 = pkin(2) * t276;
t332 = qJ(4) * t276 + t140;
t22 = (qJD(1) * t224 + qJDD(1) * t230) * t203 + t332;
t223 = -qJD(5) * t21 + t212 * t17 - t208 * t22;
t110 = -qJDD(5) + t390;
t65 = t106 * t211 + t146 * t207;
t12 = qJD(6) * t65 + t211 * t110 + t207 * t39;
t36 = qJDD(6) + t40;
t365 = t211 * t36;
t381 = -t100 ^ 2 * t207 + t365;
t348 = t203 * t212;
t251 = -t204 * t208 + t209 * t348;
t342 = t208 * t214;
t124 = t204 * t340 + t339;
t353 = t124 * t212;
t126 = -t204 * t341 + t337;
t349 = t203 * t210;
t83 = t126 * t212 - t208 * t349;
t242 = g(1) * t83 + g(2) * (t203 * t342 + t353) + g(3) * t251;
t4 = t110 * pkin(5) - t223;
t380 = t100 * (pkin(5) * t106 + pkin(10) * t100) + t242 + t4;
t133 = pkin(5) * t208 - pkin(10) * t212 + t206;
t274 = pkin(5) * t212 + pkin(10) * t208;
t379 = t100 * (qJD(5) * t274 - (-pkin(4) - t274) * t300 + t333) + t133 * t36;
t346 = t203 * t214;
t252 = -t124 * t208 + t212 * t346;
t378 = -t123 * t211 + t207 * t252;
t377 = t123 * t207 + t211 * t252;
t371 = pkin(1) * t209;
t191 = t204 * t371;
t325 = pkin(8) * t347 + t191;
t107 = -t204 * qJ(3) - t325;
t90 = pkin(3) * t347 - t107;
t62 = pkin(4) * t347 - pkin(9) * t204 + t90;
t294 = -qJ(4) * t213 - pkin(1);
t350 = t203 * t209;
t326 = pkin(2) * t347 + qJ(3) * t350;
t72 = (pkin(9) * t209 + t294) * t203 - t326;
t260 = t208 * t62 + t212 * t72;
t320 = qJD(2) * t209;
t299 = t203 * t320;
t172 = pkin(2) * t299;
t330 = qJ(4) * t299 + t172;
t48 = t203 * t224 + t330;
t373 = pkin(3) + pkin(8);
t287 = t203 * (-pkin(4) - t373);
t176 = t213 * t305;
t195 = t204 * qJD(3);
t328 = t176 + t195;
t58 = t287 * t320 + t328;
t376 = -qJD(5) * t260 - t208 * t48 + t212 * t58;
t14 = pkin(10) * t146 + t21;
t282 = (-pkin(3) - pkin(4)) * t350;
t70 = qJD(1) * t282 + t329;
t336 = qJD(3) - t70;
t43 = -t336 - t385;
t24 = t104 * pkin(5) - t106 * pkin(10) + t43;
t265 = t14 * t207 - t211 * t24;
t315 = qJD(5) * t212;
t249 = t208 * t17 + t212 * t22 + t44 * t315 - t316 * t56;
t3 = -pkin(10) * t110 + t249;
t375 = -pkin(3) * t391 + t181 * qJ(4) + t183 * qJD(4);
t278 = pkin(8) * t391 + t209 * t281 - t213 * t304;
t49 = t278 - t382;
t29 = t49 - t375;
t18 = -pkin(4) * t391 - t29;
t6 = pkin(5) * t40 - pkin(10) * t39 + t18;
t1 = -t265 * qJD(6) + t207 * t6 + t211 * t3;
t374 = 0.2e1 * t148;
t157 = qJ(4) * t301;
t177 = pkin(2) * t301;
t69 = -t205 * t300 + t157 + t177;
t370 = t208 * t70 + t212 * t69;
t369 = t100 * t63;
t313 = qJD(6) * t211;
t314 = qJD(6) * t207;
t11 = -t106 * t314 - t207 * t110 + t146 * t313 + t211 * t39;
t368 = t11 * t207;
t366 = t207 * t36;
t364 = t212 * t11;
t363 = t63 * t106;
t362 = t65 * t100;
t361 = t65 * t106;
t359 = qJ(3) * t209;
t358 = t100 * t205;
t286 = t100 * t211;
t352 = t146 * t212;
t345 = t205 * t110;
t344 = t208 * t110;
t343 = t208 * t213;
t338 = t211 * t213;
t93 = -pkin(3) * t301 + t329;
t335 = qJD(3) - t93;
t331 = qJ(3) * t179 + qJD(3) * t300;
t327 = -pkin(8) * t350 + t213 * t372;
t202 = t213 ^ 2;
t324 = t201 - t202;
t323 = qJ(3) * qJD(2);
t319 = qJD(2) * t213;
t317 = qJD(5) * t205;
t307 = g(3) * t350;
t302 = t213 * t351;
t109 = -t204 * pkin(2) - t327;
t298 = t203 * t319;
t297 = pkin(1) * t389;
t238 = t275 - t359;
t73 = t238 * t322;
t89 = t203 * t294 - t326;
t291 = qJD(1) * t89 + t73;
t289 = -t123 * pkin(2) + qJ(3) * t124;
t288 = -t125 * pkin(2) + qJ(3) * t126;
t285 = t183 + t321;
t284 = -qJD(4) - t323;
t283 = t181 + t310;
t280 = t209 * t302;
t277 = t204 * qJ(4) - t109;
t272 = g(1) * t123 - g(2) * t125;
t271 = g(1) * t126 + g(2) * t124;
t270 = g(1) * t124 - g(2) * t126;
t269 = g(1) * t214 + g(2) * t210;
t98 = (-t207 * t209 + t208 * t338) * t322;
t268 = -t211 * t316 - t98;
t8 = t14 * t211 + t207 * t24;
t26 = pkin(10) * t347 + t260;
t122 = t204 * t212 + t208 * t350;
t61 = t282 + t277;
t34 = -pkin(5) * t251 - t122 * pkin(10) + t61;
t264 = t207 * t34 + t211 * t26;
t263 = -t207 * t26 + t211 * t34;
t20 = -t208 * t56 + t212 * t44;
t261 = -t208 * t72 + t212 * t62;
t108 = -t203 * pkin(1) - t326;
t257 = -pkin(2) * t213 - pkin(1) - t359;
t96 = t257 * t322;
t259 = qJD(2) * (-qJD(1) * t108 - t96);
t258 = -pkin(8) * t299 + t176;
t256 = t214 * pkin(1) + t126 * pkin(2) + pkin(8) * t349 + qJ(3) * t125;
t112 = -qJ(3) * t300 + t177;
t255 = -t100 * t313 - t366;
t254 = -t100 * t314 + t365;
t253 = -t122 * t207 + t203 * t338;
t82 = t122 * t211 + t207 * t347;
t228 = t213 * t284 - t318;
t30 = (qJD(1) * t228 + qJDD(1) * t238) * t203 + t332;
t54 = t203 * t228 + t330;
t250 = -qJD(1) * t54 - qJDD(1) * t89 - t30;
t248 = t208 * t58 + t212 * t48 + t62 * t315 - t316 * t72;
t247 = t146 * t65;
t243 = -qJ(3) * t319 - t318;
t45 = t140 + (qJD(1) * t243 + qJDD(1) * t257) * t203;
t95 = t203 * t243 + t172;
t246 = qJD(1) * t95 + qJDD(1) * t108 + t45;
t244 = t146 * t208;
t74 = -t183 * t300 + t391;
t239 = -pkin(1) * t210 - t124 * pkin(2) + pkin(8) * t346 - qJ(3) * t123;
t116 = t325 * qJD(2);
t237 = t116 * t183 - t270;
t235 = t179 + (-qJD(2) + t183) * t301;
t234 = t100 * t286 + t366;
t13 = -pkin(5) * t146 - t20;
t233 = -qJD(3) * t100 - qJD(5) * t13 - t205 * t36;
t229 = -t271 - t279;
t227 = -pkin(10) * t36 + (t13 + t20) * t100;
t226 = t18 - t232;
t225 = -t278 - t232;
t2 = -qJD(6) * t8 - t207 * t3 + t211 * t6;
t222 = qJD(6) * t358 + t232;
t221 = t73 * t300 + t229 + t383;
t194 = t204 * qJD(4);
t59 = t194 + (t213 * t287 - t191) * qJD(2);
t220 = t225 + t382;
t219 = t115 * t183 + t225;
t218 = t307 + (-pkin(10) * t301 - qJD(6) * t133 + t370) * t100 + t271;
t216 = t220 + t375;
t113 = t181 + t280;
t103 = t212 * t110;
t99 = -t195 - t258;
t97 = (t207 * t343 + t209 * t211) * t322;
t92 = -t156 - t115;
t91 = t112 + t157;
t88 = -t169 + t312;
t84 = t126 * t208 + t210 * t348;
t80 = qJD(5) * t122 - t212 * t298;
t79 = qJD(5) * t251 + t208 * t298;
t78 = -t194 + (t347 * t373 + t191) * qJD(2);
t77 = -t299 * t373 + t328;
t76 = pkin(3) * t350 - t277;
t75 = t96 * t301;
t68 = t156 + t333;
t51 = t335 + t385;
t47 = -t125 * t207 + t211 * t84;
t46 = -t125 * t211 - t207 * t84;
t33 = qJD(6) * t253 - t207 * t299 + t79 * t211;
t32 = qJD(6) * t82 + t79 * t207 + t211 * t299;
t31 = -pkin(3) * t276 + t231;
t27 = pkin(5) * t301 + t208 * t69 - t212 * t70;
t25 = -pkin(5) * t347 - t261;
t23 = t80 * pkin(5) - t79 * pkin(10) + t59;
t10 = pkin(5) * t299 - t376;
t9 = -pkin(10) * t299 + t248;
t5 = [qJDD(1), g(1) * t210 - g(2) * t214, t269 (qJDD(1) * t201 + 0.2e1 * t209 * t295) * t200 (t209 * t308 - t311 * t324) * t389 (t209 * t283 + t285 * t319) * t203 (t213 * t283 - t285 * t320) * t203, t181 * t204, t181 * t327 - t204 * t278 + t240 * t297 - t237, -t181 * t325 - t183 * t258 + t204 * t279 - t297 * t393 - t272 ((qJD(2) * t88 - qJDD(1) * t107 - t41 + (qJD(2) * t109 - t99) * qJD(1)) * t213 + (qJD(2) * t92 + qJDD(1) * t109 + t49 + (qJD(2) * t107 + t116) * qJD(1)) * t209 - t269) * t203, t109 * t181 + t49 * t204 + (t209 * t259 + t213 * t246) * t203 + t237, -t107 * t181 - t99 * t183 - t41 * t204 + (-t209 * t246 + t213 * t259) * t203 + t272, -g(1) * t239 - g(2) * t256 + t41 * t107 + t45 * t108 + t49 * t109 + t88 * t116 + t92 * t99 + t96 * t95 ((qJD(2) * t51 + qJDD(1) * t90 + t31 + (qJD(2) * t76 + t77) * qJD(1)) * t213 + (-qJD(2) * t68 + qJDD(1) * t76 + t29 + (-qJD(2) * t90 + t78) * qJD(1)) * t209 - t269) * t203, t90 * t181 + t77 * t183 + t31 * t204 + (t209 * t250 - t291 * t319) * t203 + t272, -t76 * t181 - t78 * t183 - t29 * t204 + (t213 * t250 + t291 * t320) * t203 + t270, t30 * t89 + t73 * t54 + t29 * t76 + t51 * t78 + t31 * t90 + t68 * t77 - g(1) * (pkin(3) * t346 - qJ(4) * t124 + t239) - g(2) * (pkin(3) * t349 + qJ(4) * t126 + t256) t106 * t79 + t122 * t39, -t104 * t79 - t106 * t80 - t122 * t40 + t251 * t39, -t122 * t110 + t79 * t146 + (-t106 * t320 + t213 * t39) * t203, -t251 * t110 - t80 * t146 + (t104 * t320 - t213 * t40) * t203 (-t110 * t213 - t146 * t320) * t203, t376 * t146 - t261 * t110 + t59 * t104 + t61 * t40 - t18 * t251 + t43 * t80 - g(1) * t252 - g(2) * t84 + (-t20 * t320 + t213 * t223) * t203, -t248 * t146 + t260 * t110 + t59 * t106 + t61 * t39 + t18 * t122 + t43 * t79 + g(1) * t353 - g(2) * t83 + (g(1) * t342 + t21 * t320 - t213 * t249) * t203, t11 * t82 + t33 * t65, t11 * t253 - t12 * t82 - t32 * t65 - t33 * t63, t100 * t33 - t11 * t251 + t36 * t82 + t65 * t80, -t100 * t32 + t12 * t251 + t253 * t36 - t63 * t80, t100 * t80 - t251 * t36 (-qJD(6) * t264 - t207 * t9 + t211 * t23) * t100 + t263 * t36 - t2 * t251 - t265 * t80 + t10 * t63 + t25 * t12 - t4 * t253 + t13 * t32 - g(1) * t377 - g(2) * t47 -(qJD(6) * t263 + t207 * t23 + t211 * t9) * t100 - t264 * t36 + t1 * t251 - t8 * t80 + t10 * t65 + t25 * t11 + t4 * t82 + t13 * t33 + g(1) * t378 - g(2) * t46; 0, 0, 0, -t280, t324 * t351, t74, t235, t181, t351 * t371 + t219, pkin(1) * t302 + t183 * t329 - t229 + t307 (-pkin(2) * t309 + ((-pkin(2) * qJD(2) - t329 - t88) * t213 + (-t115 - t92 - t323) * t209) * qJD(1)) * t203 + t331, -t112 * t300 + qJDD(3) - 0.2e1 * t154 - t219 + t75, t374 + t150 + t312 * t183 + (-g(3) * t209 + (t112 * t209 + t213 * t96) * qJD(1)) * t203 + t229, -t49 * pkin(2) - g(1) * t288 - g(2) * t289 - g(3) * t326 - t41 * qJ(3) - t96 * t112 - t88 * t115 - t312 * t92 (-t206 * t309 + ((-qJD(2) * t206 - t51 - t93) * t213 + (t284 + t68 - t94) * t209) * qJD(1)) * t203 + t331, -t183 * t93 + t374 + 0.2e1 * t150 + (-g(3) + (-pkin(3) * qJD(2) + t91) * qJD(1)) * t350 + t221, t216 + t333 * t183 + t181 * t206 + (-t209 * t73 + t213 * t91) * t322, -t29 * t206 + t31 * qJ(3) - t73 * t91 - g(1) * (-qJ(4) * t125 + t288) - g(2) * (-qJ(4) * t123 + t289) - g(3) * (qJ(4) * t347 + t326) + t335 * t68 - t333 * t51, -t106 * t244 + t39 * t212, t387 * t208 + t212 * t392, -t146 * t316 - t103 + (t106 * t209 - t146 * t343) * t322, -t146 * t315 + t344 + (-t104 * t209 - t213 * t352) * t322, t146 * t301, t20 * t301 + t206 * t40 + t334 * t104 + (-t345 + (t336 + t43) * t146) * t212 + ((t69 - t317) * t146 + t226) * t208, t334 * t106 + t206 * t39 + t345 * t208 - t21 * t301 + t226 * t212 + (t370 + (-qJD(3) - t43) * t208 - t317 * t212) * t146, t211 * t364 + (-t212 * t314 + t268) * t65, t98 * t63 + t65 * t97 + (t207 * t65 + t211 * t63) * t316 + (-t368 - t12 * t211 + (t207 * t63 - t211 * t65) * qJD(6)) * t212, t11 * t208 + t268 * t100 + (t247 + t254) * t212, -t12 * t208 + (t207 * t316 + t97) * t100 + (t255 - t388) * t212, t100 * t352 + t36 * t208, -t13 * t97 - t27 * t63 + t379 * t211 + t218 * t207 + (t207 * t233 - t211 * t222 + t317 * t63 + t2) * t208 + (-t265 * t300 + t13 * t313 - qJD(3) * t63 - t205 * t12 + t4 * t207 + (-t207 * t358 - t265) * qJD(5)) * t212, -t13 * t98 - t27 * t65 - t379 * t207 + t218 * t211 + (t207 * t222 + t211 * t233 + t317 * t65 - t1) * t208 + (-t8 * t300 - t13 * t314 - qJD(3) * t65 - t205 * t11 + t4 * t211 + (-t205 * t286 - t8) * qJD(5)) * t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t113, t384, t183 * t92 - t220 + t75, t74, t384, -t113, -t183 * t68 + t301 * t73 - t216, 0, 0, 0, 0, 0, t392, t387, 0, 0, 0, 0, 0, t363 - t381, t234 + t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t181 - t280, -t202 * t351 - t180, t183 * t51 + (-pkin(3) * t311 - g(3)) * t350 + t221 - t386, 0, 0, 0, 0, 0, -t183 * t104 - t146 * t244 - t103, -t183 * t106 - t146 * t352 + t344, 0, 0, 0, 0, 0, -t212 * t12 + (-t211 * t183 - t207 * t352) * t100 + (t255 + t388) * t208, -t364 + (t207 * t183 - t211 * t352) * t100 + (t247 - t254) * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t104, -t104 ^ 2 + t106 ^ 2, t39 + t357, -t40 + t356, -t110, -t43 * t106 + t21 * t146 + t223 - t242, g(1) * t84 - g(2) * t252 + g(3) * t122 + t104 * t43 + t146 * t20 - t249, t286 * t65 + t368 (t11 - t369) * t211 + (-t12 - t362) * t207, t234 - t361, t363 + t381, -t100 * t106, -pkin(5) * t12 + t106 * t265 + t227 * t207 - t21 * t63 - t211 * t380, -pkin(5) * t11 + t8 * t106 + t207 * t380 - t21 * t65 + t227 * t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t63, -t63 ^ 2 + t65 ^ 2, t11 + t369, -t12 + t362, t36, -g(1) * t46 - g(2) * t378 - g(3) * t253 + t8 * t100 - t13 * t65 + t2, g(1) * t47 - g(2) * t377 + g(3) * t82 - t100 * t265 + t13 * t63 - t1;];
tau_reg  = t5;
