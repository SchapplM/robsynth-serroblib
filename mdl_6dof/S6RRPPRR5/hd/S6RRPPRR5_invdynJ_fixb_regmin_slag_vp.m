% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:11:03
% EndTime: 2019-03-09 09:11:19
% DurationCPUTime: 7.40s
% Computational Cost: add. (4491->560), mult. (10987->715), div. (0->0), fcn. (8249->10), ass. (0->278)
t206 = sin(pkin(6));
t215 = cos(qJ(2));
t318 = qJD(1) * qJD(2);
t298 = t215 * t318;
t211 = sin(qJ(2));
t316 = qJDD(1) * t211;
t243 = t298 + t316;
t395 = t243 * t206;
t329 = qJD(1) * t206;
t304 = t215 * t329;
t388 = -t304 - qJD(5);
t207 = cos(pkin(6));
t320 = t207 * qJD(1);
t182 = qJD(2) + t320;
t210 = sin(qJ(5));
t214 = cos(qJ(5));
t305 = t211 * t329;
t99 = -t182 * t210 + t214 * t305;
t356 = t99 * t388;
t317 = qJDD(1) * t207;
t179 = qJDD(2) + t317;
t323 = qJD(5) * t214;
t327 = qJD(2) * t215;
t238 = t210 * t327 + t211 * t323;
t324 = qJD(5) * t210;
t38 = t214 * t179 - t182 * t324 + t206 * (qJD(1) * t238 + t210 * t316);
t394 = t38 - t356;
t280 = t210 * t305;
t227 = qJD(5) * t280 + t210 * t179 + t182 * t323 - t214 * t395;
t383 = t214 * t182 + t280;
t362 = t388 * t383;
t393 = t227 - t362;
t209 = sin(qJ(6));
t93 = qJD(6) + t383;
t290 = t93 ^ 2;
t213 = cos(qJ(6));
t36 = qJDD(6) + t38;
t357 = t213 * t36;
t392 = -t209 * t290 + t357;
t60 = t209 * t99 + t213 * t388;
t391 = t60 * t99;
t62 = -t209 * t388 + t213 * t99;
t390 = t62 * t99;
t381 = t179 * qJ(3) + t182 * qJD(3);
t216 = cos(qJ(1));
t337 = t215 * t216;
t212 = sin(qJ(1));
t342 = t211 * t212;
t118 = -t207 * t337 + t342;
t340 = t212 * t215;
t341 = t211 * t216;
t119 = t207 * t341 + t340;
t343 = t206 * t216;
t78 = t119 * t214 + t210 * t343;
t387 = t118 * t213 + t209 * t78;
t386 = -t118 * t209 + t213 * t78;
t202 = t206 ^ 2;
t385 = 0.2e1 * t202;
t384 = t388 * t62;
t344 = t206 * t215;
t372 = pkin(1) * t211;
t245 = pkin(8) * t344 + t207 * t372;
t110 = t245 * qJD(2);
t347 = t206 * t211;
t176 = qJD(4) * t347;
t302 = t206 * t327;
t67 = -qJ(4) * t302 + t110 - t176;
t382 = 0.2e1 * t381;
t120 = t207 * t340 + t341;
t235 = -g(1) * t120 - g(2) * t118 + g(3) * t344;
t217 = -pkin(2) - pkin(3);
t380 = qJD(2) * t217;
t314 = pkin(1) * t320;
t108 = -pkin(8) * t305 + t215 * t314;
t319 = qJD(3) - t108;
t149 = t179 * pkin(2);
t379 = t149 - qJDD(3);
t116 = t207 * t214 + t210 * t347;
t338 = t214 * t216;
t350 = t119 * t210;
t121 = -t207 * t342 + t337;
t345 = t206 * t214;
t80 = -t121 * t210 - t212 * t345;
t244 = g(1) * t80 + g(2) * (t206 * t338 - t350) - g(3) * t116;
t315 = qJDD(1) * t215;
t175 = t206 * t315;
t299 = t211 * t318;
t279 = t206 * t299;
t105 = -qJDD(5) - t175 + t279;
t203 = pkin(4) - t217;
t246 = -pkin(9) * t215 - t203 * t211;
t234 = qJD(2) * t246;
t177 = qJD(3) * t347;
t197 = t206 * pkin(1);
t261 = pkin(2) * t175 + qJ(3) * t395 + qJD(1) * t177 + qJDD(1) * t197;
t240 = pkin(3) * t175 + qJDD(4) + t261;
t277 = pkin(4) * t215 - pkin(9) * t211;
t17 = (qJD(1) * t234 + t277 * qJDD(1)) * t206 + t240;
t255 = t277 * t206;
t90 = -pkin(1) * t329 - pkin(2) * t304 - qJ(3) * t305;
t68 = pkin(3) * t304 + qJD(4) - t90;
t47 = qJD(1) * t255 + t68;
t150 = t182 * qJ(3);
t109 = pkin(8) * t304 + t211 * t314;
t88 = -qJ(4) * t304 + t109;
t64 = t150 + t88;
t55 = -pkin(9) * t182 + t64;
t19 = t210 * t47 + t214 * t55;
t133 = qJ(4) * t279;
t300 = qJ(4) * t315;
t373 = pkin(1) * t207;
t313 = qJD(2) * t373;
t283 = qJD(1) * t313;
t312 = pkin(1) * t317;
t306 = pkin(8) * t175 + t211 * t312 + t215 * t283;
t326 = qJD(4) * t215;
t328 = qJD(2) * t211;
t30 = t133 + (-t300 + (-pkin(8) * t328 - t326) * qJD(1)) * t206 + t306 + t381;
t24 = -t179 * pkin(9) + t30;
t225 = -qJD(5) * t19 + t214 * t17 - t210 * t24;
t4 = t105 * pkin(5) - t225;
t378 = (pkin(5) * t99 + pkin(10) * t93) * t93 + t244 + t4;
t127 = pkin(5) * t214 + pkin(10) * t210 + t203;
t377 = (t88 + t388 * (pkin(5) * t210 - pkin(10) * t214)) * t93 + t127 * t36;
t12 = qJD(6) * t62 + t213 * t105 - t209 * t227;
t332 = pkin(2) * t344 + qJ(3) * t347;
t103 = -t197 - t332;
t185 = pkin(3) * t344;
t84 = t185 - t103;
t58 = t255 + t84;
t102 = t207 * qJ(3) + t245;
t83 = -qJ(4) * t344 + t102;
t70 = -pkin(9) * t207 + t83;
t262 = t210 * t58 + t214 * t70;
t334 = qJ(3) * t302 + t177;
t44 = t206 * t234 + t334;
t174 = t215 * t313;
t194 = t207 * qJD(3);
t59 = t174 + t194 + (-t326 + (-pkin(8) + qJ(4)) * t328) * t206;
t375 = -qJD(5) * t262 - t210 * t59 + t214 * t44;
t16 = -pkin(10) * t388 + t19;
t151 = qJ(4) * t305;
t82 = -t182 * pkin(2) + t319;
t53 = -t182 * pkin(3) - t151 + t82;
t42 = pkin(4) * t182 - t53;
t20 = pkin(5) * t383 - pkin(10) * t99 + t42;
t267 = t16 * t209 - t20 * t213;
t252 = t210 * t17 + t214 * t24 + t47 * t323 - t55 * t324;
t3 = -pkin(10) * t105 + t252;
t374 = -t179 * pkin(3) - qJ(4) * t395 - qJD(1) * t176;
t281 = pkin(8) * t395 + t211 * t283 - t215 * t312;
t49 = t281 - t379;
t29 = t49 + t374;
t23 = pkin(4) * t179 - t29;
t8 = pkin(5) * t38 + pkin(10) * t227 + t23;
t1 = -t267 * qJD(6) + t209 * t8 + t213 * t3;
t371 = g(3) * t211;
t370 = t60 * t93;
t369 = t62 * t93;
t339 = t214 * t215;
t91 = (t209 * t339 + t211 * t213) * t329;
t368 = t91 * t93;
t92 = (-t209 * t211 + t213 * t339) * t329;
t367 = t92 * t93;
t154 = qJ(3) * t304;
t56 = t246 * t329 + t154;
t86 = t151 + t108;
t366 = t210 * t56 + t214 * t86;
t365 = pkin(8) * qJD(2);
t321 = qJD(6) * t213;
t322 = qJD(6) * t209;
t11 = -t209 * t105 - t213 * t227 - t321 * t388 - t99 * t322;
t364 = t11 * t209;
t208 = qJ(3) - pkin(9);
t361 = t208 * t93;
t360 = t209 * t36;
t359 = t211 * t383;
t358 = t211 * t99;
t291 = t213 * t93;
t354 = qJ(3) * t215;
t353 = qJD(5) * t93;
t349 = t388 * t210;
t348 = t202 * qJD(1) ^ 2;
t346 = t206 * t212;
t336 = qJD(3) - t86;
t333 = -pkin(8) * t347 + t215 * t373;
t204 = t211 ^ 2;
t205 = t215 ^ 2;
t331 = t204 - t205;
t330 = qJ(3) * qJD(2);
t325 = qJD(5) * t208;
t311 = t209 * t353;
t310 = t213 * t353;
t309 = t215 * t349;
t308 = t388 * t339;
t307 = t215 * t348;
t104 = -t207 * pkin(2) - t333;
t303 = t206 * t328;
t301 = pkin(1) * t385;
t292 = qJD(1) * t84 + t68;
t289 = -t118 * pkin(2) + qJ(3) * t119;
t288 = -t120 * pkin(2) + qJ(3) * t121;
t287 = qJD(1) * t103 + t90;
t286 = t182 + t320;
t285 = t179 + t317;
t284 = t217 * t347;
t282 = t211 * t307;
t274 = g(1) * t118 - g(2) * t120;
t273 = g(1) * t121 + g(2) * t119;
t272 = g(1) * t119 - g(2) * t121;
t271 = g(1) * t216 + g(2) * t212;
t270 = g(1) * t212 - g(2) * t216;
t6 = t16 * t213 + t20 * t209;
t26 = pkin(10) * t344 + t262;
t117 = -t207 * t210 + t211 * t345;
t71 = -t207 * pkin(3) - qJ(4) * t347 + t104;
t63 = t207 * pkin(4) - t71;
t34 = pkin(5) * t116 - pkin(10) * t117 + t63;
t266 = t209 * t34 + t213 * t26;
t265 = -t209 * t26 + t213 * t34;
t18 = -t210 * t55 + t214 * t47;
t263 = -t210 * t70 + t214 * t58;
t260 = qJD(2) * t284;
t259 = -pkin(8) * t303 + t174;
t258 = t216 * pkin(1) + t121 * pkin(2) + pkin(8) * t346 + qJ(3) * t120;
t257 = -t93 * t321 - t360;
t256 = t93 * t322 - t357;
t254 = -t117 * t209 + t213 * t344;
t76 = t117 * t213 + t209 * t344;
t31 = qJD(1) * t260 + t240;
t66 = t260 + t334;
t253 = qJD(1) * t66 + qJDD(1) * t84 + t31;
t251 = t210 * t44 + t214 * t59 + t58 * t323 - t70 * t324;
t250 = t388 * t60;
t249 = t210 * t105 + t323 * t388;
t248 = -t214 * t105 + t324 * t388;
t43 = pkin(2) * t279 - t261;
t89 = pkin(2) * t303 - t334;
t247 = -qJD(1) * t89 - qJDD(1) * t103 - t43;
t242 = -pkin(1) * t212 - t119 * pkin(2) + pkin(8) * t343 - qJ(3) * t118;
t241 = -t110 * t182 + t272;
t15 = pkin(5) * t388 - t18;
t237 = -qJD(3) * t93 - qJD(5) * t15 - t208 * t36;
t236 = t273 - t306;
t232 = -pkin(10) * t36 + (t15 + t18) * t93;
t231 = -pkin(8) * t279 + t306;
t230 = t23 - t235;
t229 = t208 * t105 + t388 * t42;
t228 = -t281 - t235;
t2 = -t6 * qJD(6) - t209 * t3 + t213 * t8;
t226 = t108 * t182 + t236;
t224 = qJD(6) * t361 + t235;
t223 = t228 + t379;
t222 = t109 * t182 + t228;
t221 = g(3) * t347 + (-pkin(10) * t305 - qJD(6) * t127 + t366) * t93 + t273;
t220 = t223 - t374;
t122 = t182 * t304;
t115 = -t182 ^ 2 - t204 * t348;
t107 = -t179 - t282;
t106 = pkin(2) * t305 - t154;
t94 = t194 + t259;
t87 = qJD(1) * t284 + t154;
t85 = t150 + t109;
t81 = t121 * t214 - t210 * t346;
t74 = -qJD(5) * t116 + t214 * t302;
t73 = t206 * t238 - t207 * t324;
t69 = -t122 + t395;
t46 = -t120 * t209 + t213 * t81;
t45 = -t120 * t213 - t209 * t81;
t39 = t231 + t381;
t33 = qJD(6) * t254 - t209 * t303 + t74 * t213;
t32 = qJD(6) * t76 + t74 * t209 + t213 * t303;
t27 = pkin(5) * t305 + t210 * t86 - t214 * t56;
t25 = -pkin(5) * t344 - t263;
t22 = t73 * pkin(5) - t74 * pkin(10) - t67;
t10 = pkin(5) * t303 - t375;
t9 = -pkin(10) * t303 + t251;
t5 = [qJDD(1), t270, t271 (qJDD(1) * t204 + 0.2e1 * t211 * t298) * t202 (t211 * t315 - t331 * t318) * t385 (t285 * t211 + t286 * t327) * t206 (t215 * t285 - t286 * t328) * t206, t179 * t207, t333 * t179 - t281 * t207 + (-t299 + t315) * t301 + t241, -t245 * t179 - t259 * t182 - t231 * t207 - t243 * t301 - t274, -t104 * t179 - t49 * t207 + (t215 * t247 + t287 * t328) * t206 + t241 ((qJD(2) * t82 + qJDD(1) * t102 + t39 + (qJD(2) * t104 + t94) * qJD(1)) * t215 + (-qJD(2) * t85 + qJDD(1) * t104 + t49 + (-qJD(2) * t102 + t110) * qJD(1)) * t211 - t271) * t206, t102 * t179 + t94 * t182 + t39 * t207 + (t211 * t247 - t287 * t327) * t206 + t274, -g(1) * t242 - g(2) * t258 + t39 * t102 + t43 * t103 + t49 * t104 + t82 * t110 + t85 * t94 + t90 * t89, -t71 * t179 - t67 * t182 - t29 * t207 + (t215 * t253 - t292 * t328) * t206 + t272, t83 * t179 + t59 * t182 + t30 * t207 + (t211 * t253 + t292 * t327) * t206 + t274 ((-qJD(2) * t53 - qJDD(1) * t83 - t30 + (-qJD(2) * t71 - t59) * qJD(1)) * t215 + (qJD(2) * t64 - qJDD(1) * t71 - t29 + (qJD(2) * t83 - t67) * qJD(1)) * t211 + t271) * t206, t30 * t83 + t64 * t59 + t29 * t71 + t53 * t67 + t31 * t84 + t68 * t66 - g(1) * (-pkin(3) * t119 - qJ(4) * t343 + t242) - g(2) * (pkin(3) * t121 - qJ(4) * t346 + t258) -t117 * t227 + t74 * t99, t116 * t227 - t117 * t38 - t383 * t74 - t73 * t99, -t117 * t105 - t74 * t388 + (-t215 * t227 - t99 * t328) * t206, t116 * t105 + t73 * t388 + (-t215 * t38 + t328 * t383) * t206 (-t105 * t215 + t328 * t388) * t206, -t375 * t388 - t263 * t105 - t67 * t383 + t63 * t38 + t23 * t116 + t42 * t73 + g(1) * t78 - g(2) * t81 + (-t18 * t328 + t215 * t225) * t206, t251 * t388 + t262 * t105 - t67 * t99 - t63 * t227 + t23 * t117 + t42 * t74 - g(1) * t350 - g(2) * t80 + (g(1) * t338 + t19 * t328 - t215 * t252) * t206, t11 * t76 + t33 * t62, t11 * t254 - t12 * t76 - t32 * t62 - t33 * t60, t11 * t116 + t33 * t93 + t36 * t76 + t62 * t73, -t116 * t12 + t254 * t36 - t32 * t93 - t60 * t73, t116 * t36 + t73 * t93 (-qJD(6) * t266 - t209 * t9 + t213 * t22) * t93 + t265 * t36 + t2 * t116 - t267 * t73 + t10 * t60 + t25 * t12 - t4 * t254 + t15 * t32 + g(1) * t386 - g(2) * t46 -(qJD(6) * t265 + t209 * t22 + t213 * t9) * t93 - t266 * t36 - t1 * t116 - t6 * t73 + t10 * t62 + t25 * t11 + t4 * t76 + t15 * t33 - g(1) * t387 - g(2) * t45; 0, 0, 0, -t282, t331 * t348, t69, t175 + (-qJD(2) + t182) * t305, t179, t348 * t372 + t222, pkin(1) * t307 + (pkin(8) * t318 + g(3)) * t347 + t226, 0.2e1 * t149 - qJDD(3) + (t106 * t215 - t211 * t90) * t329 + t222 ((-pkin(2) * t211 + t354) * qJDD(1) + ((-t109 + t85 - t330) * t211 + (-pkin(2) * qJD(2) + t319 - t82) * t215) * qJD(1)) * t206 (-t371 + (t215 * t90 + (t106 - t365) * t211) * qJD(1)) * t206 - t226 + t382, -t49 * pkin(2) - g(1) * t288 - g(2) * t289 - g(3) * t332 + t39 * qJ(3) - t90 * t106 - t82 * t109 + t319 * t85, t220 - t179 * t217 + t182 * t88 + (t211 * t68 - t215 * t87) * t329, -t182 * t86 + t133 + (-t300 - t371 + ((-qJD(4) - t68) * t215 + (-t87 - t365) * t211) * qJD(1)) * t206 - t236 + t382 ((-t211 * t217 - t354) * qJDD(1) + ((-t64 + t88 + t330) * t211 + (-t336 + t53 - t380) * t215) * qJD(1)) * t206, t30 * qJ(3) + t29 * t217 - t53 * t88 - t68 * t87 - g(1) * (-pkin(3) * t120 + t288) - g(2) * (-pkin(3) * t118 + t289) - g(3) * (t185 + t332) + t336 * t64, t210 * t227 + t214 * t356, t394 * t210 + t393 * t214 (t308 + t358) * t329 + t249 (-t309 - t359) * t329 - t248, -t388 * t305, t18 * t305 + t203 * t38 + t88 * t383 + (t336 * t388 + t229) * t210 + (-(-t56 - t325) * t388 + t230) * t214, -t203 * t227 - t366 * t388 - t19 * t305 + t88 * t99 + (qJD(3) * t388 + t229) * t214 + (-t325 * t388 - t230) * t210, -t11 * t213 * t210 + (t210 * t322 - t213 * t323 - t92) * t62, t92 * t60 + t62 * t91 + (t209 * t62 + t213 * t60) * t323 + (t364 + t12 * t213 + (-t209 * t60 + t213 * t62) * qJD(6)) * t210, -t367 + (t11 - t310) * t214 + (t256 + t384) * t210, t368 + (-t12 + t311) * t214 + (-t250 - t257) * t210, t36 * t214 + t349 * t93, -t15 * t91 - t27 * t60 + t377 * t213 + t221 * t209 + (t209 * t237 - t213 * t224 + t325 * t60 + t2) * t214 + (t267 * t304 - t15 * t321 + qJD(3) * t60 + t208 * t12 - t4 * t209 + (t209 * t361 + t267) * qJD(5)) * t210, -t15 * t92 - t27 * t62 - t377 * t209 + t221 * t213 + (t209 * t224 + t213 * t237 + t325 * t62 - t1) * t214 + (t6 * t304 + t15 * t322 + qJD(3) * t62 + t208 * t11 - t4 * t213 + (t208 * t291 + t6) * qJD(5)) * t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t69, t115, -t182 * t85 + t305 * t90 - t223, t107, t115, -t69, -t182 * t64 - t305 * t68 - t220, 0, 0, 0, 0, 0, -t394, t393, 0, 0, 0, 0, 0, t391 - t392, t213 * t290 + t360 + t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175 + (-qJD(2) - t182) * t305, t122 + t395 (-t204 - t205) * t348, g(3) * t207 + ((t215 * t64 + (t53 + t380) * t211) * qJD(1) + t270) * t206 + t240, 0, 0, 0, 0, 0 (t309 - t359) * t329 + t248 (t308 - t358) * t329 + t249, 0, 0, 0, 0, 0, -t368 + (-t12 - t311) * t214 + (-t250 + t257) * t210, -t367 + (-t11 - t310) * t214 + (t256 - t384) * t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 * t383, -t383 ^ 2 + t99 ^ 2, -t227 - t362, -t38 - t356, -t105, -t19 * t388 - t42 * t99 + t225 - t244, g(1) * t81 + g(2) * t78 + g(3) * t117 - t18 * t388 + t383 * t42 - t252, t291 * t62 + t364 (t11 - t370) * t213 + (-t12 - t369) * t209, t291 * t93 + t360 - t390, t391 + t392, -t93 * t99, -pkin(5) * t12 - t19 * t60 + t232 * t209 - t213 * t378 + t267 * t99, -pkin(5) * t11 - t19 * t62 + t209 * t378 + t232 * t213 + t6 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, t11 + t370, -t12 + t369, t36, -g(1) * t45 + g(2) * t387 - g(3) * t254 - t15 * t62 + t6 * t93 + t2, g(1) * t46 + g(2) * t386 + g(3) * t76 + t15 * t60 - t267 * t93 - t1;];
tau_reg  = t5;
