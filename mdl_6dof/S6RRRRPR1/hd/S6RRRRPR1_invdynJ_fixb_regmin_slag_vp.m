% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:54:59
% EndTime: 2019-03-09 21:55:16
% DurationCPUTime: 6.94s
% Computational Cost: add. (14860->472), mult. (36593->615), div. (0->0), fcn. (28426->18), ass. (0->281)
t271 = cos(qJ(6));
t361 = qJD(6) * t271;
t269 = sin(qJ(2));
t418 = cos(qJ(3));
t348 = qJD(1) * t418;
t268 = sin(qJ(3));
t273 = cos(qJ(2));
t377 = t268 * t273;
t186 = -qJD(1) * t377 - t269 * t348;
t267 = sin(qJ(4));
t272 = cos(qJ(4));
t366 = qJD(1) * t269;
t435 = -t268 * t366 + t273 * t348;
t143 = t186 * t267 + t272 * t435;
t264 = sin(pkin(11));
t265 = cos(pkin(11));
t313 = t272 * t186 - t267 * t435;
t427 = t143 * t265 + t264 * t313;
t438 = t271 * t427;
t441 = t361 - t438;
t260 = qJD(2) + qJD(3);
t440 = t260 * t435;
t259 = qJDD(2) + qJDD(3);
t252 = qJDD(4) + t259;
t253 = qJD(4) + t260;
t266 = sin(qJ(6));
t362 = qJD(6) * t266;
t426 = t143 * t264 - t265 * t313;
t342 = qJDD(1) * t418;
t358 = t273 * qJDD(1);
t120 = t268 * t358 + t269 * t342 + t440;
t199 = t269 * t418 + t377;
t154 = t260 * t199;
t359 = t269 * qJDD(1);
t317 = t268 * t359 - t273 * t342;
t121 = qJD(1) * t154 + t317;
t288 = qJD(4) * t313 - t120 * t267 - t272 * t121;
t363 = qJD(4) * t272;
t364 = qJD(4) * t267;
t71 = t272 * t120 - t267 * t121 + t186 * t364 + t363 * t435;
t43 = t264 * t288 + t265 * t71;
t25 = t266 * t252 + t253 * t361 + t271 * t43 - t362 * t426;
t23 = t25 * t271;
t84 = t253 * t266 + t271 * t426;
t26 = qJD(6) * t84 - t271 * t252 + t266 * t43;
t82 = -t271 * t253 + t266 * t426;
t439 = -t266 * t26 - t441 * t82 + t23;
t22 = t25 * t266;
t13 = t441 * t84 + t22;
t42 = -t264 * t71 + t265 * t288;
t39 = qJDD(6) - t42;
t36 = t266 * t39;
t91 = qJD(6) - t427;
t404 = t91 * t361 + t36;
t408 = t84 * t426;
t12 = -t438 * t91 + t404 - t408;
t397 = t266 * t91;
t179 = t186 * pkin(9);
t419 = pkin(7) + pkin(8);
t216 = t419 * t273;
t205 = qJD(1) * t216;
t187 = t268 * t205;
t215 = t419 * t269;
t203 = qJD(1) * t215;
t400 = qJD(2) * pkin(2);
t193 = -t203 + t400;
t327 = t418 * t193 - t187;
t118 = t179 + t327;
t107 = pkin(3) * t260 + t118;
t191 = t418 * t205;
t303 = -t268 * t193 - t191;
t412 = t435 * pkin(9);
t119 = -t303 + t412;
t111 = t272 * t119;
t315 = -t107 * t267 - t111;
t389 = qJ(5) * t143;
t65 = -t315 + t389;
t399 = t264 * t65;
t136 = t313 * qJ(5);
t109 = t267 * t119;
t333 = t272 * t107 - t109;
t64 = t333 + t136;
t59 = pkin(4) * t253 + t64;
t32 = t265 * t59 - t399;
t30 = -pkin(5) * t253 - t32;
t411 = t30 * t427;
t409 = t82 * t426;
t407 = t91 * t426;
t326 = t203 * t268 - t191;
t123 = t326 - t412;
t370 = -t418 * t203 - t187;
t124 = t179 + t370;
t248 = pkin(2) * t418 + pkin(3);
t379 = t267 * t268;
t437 = -t248 * t363 - (-t268 * t364 + (t272 * t418 - t379) * qJD(3)) * pkin(2) + t267 * t123 + t272 * t124;
t378 = t268 * t272;
t436 = -t248 * t364 + (-t268 * t363 + (-t267 * t418 - t378) * qJD(3)) * pkin(2) - t272 * t123 + t124 * t267;
t263 = qJ(2) + qJ(3);
t257 = qJ(4) + t263;
t245 = sin(t257);
t246 = cos(t257);
t270 = sin(qJ(1));
t274 = cos(qJ(1));
t320 = g(1) * t274 + g(2) * t270;
t434 = -g(3) * t246 + t245 * t320;
t61 = t265 * t65;
t33 = t264 * t59 + t61;
t31 = pkin(10) * t253 + t33;
t258 = t273 * pkin(2);
t406 = pkin(1) + t258;
t214 = t406 * qJD(1);
t157 = -pkin(3) * t435 - t214;
t105 = -pkin(4) * t143 + qJD(5) + t157;
t48 = -pkin(5) * t427 - pkin(10) * t426 + t105;
t18 = -t266 * t31 + t271 * t48;
t243 = pkin(11) + t257;
t230 = sin(t243);
t375 = t271 * t274;
t376 = t270 * t271;
t312 = -t18 * t426 + t30 * t362 + (g(1) * t375 + g(2) * t376) * t230;
t19 = t266 * t48 + t271 * t31;
t231 = cos(t243);
t360 = qJD(1) * qJD(2);
t345 = t273 * t360;
t155 = qJDD(2) * pkin(2) + t419 * (-t345 - t359);
t346 = t269 * t360;
t156 = t419 * (-t346 + t358);
t287 = qJD(3) * t303 + t418 * t155 - t268 * t156;
t63 = t259 * pkin(3) - t120 * pkin(9) + t287;
t347 = qJD(3) * t418;
t365 = qJD(3) * t268;
t285 = t268 * t155 + t156 * t418 + t193 * t347 - t205 * t365;
t69 = -t121 * pkin(9) + t285;
t290 = qJD(4) * t315 - t267 * t69 + t272 * t63;
t15 = pkin(4) * t252 - qJ(5) * t71 + qJD(5) * t313 + t290;
t421 = (qJD(4) * t107 + t69) * t272 - t119 * t364 + t267 * t63;
t17 = qJ(5) * t288 + qJD(5) * t143 + t421;
t5 = t15 * t265 - t17 * t264;
t3 = -pkin(5) * t252 - t5;
t431 = g(3) * t231 + t3;
t316 = t19 * t426 + t266 * t431 + t30 * t361;
t433 = pkin(5) * t426 - pkin(10) * t427;
t432 = pkin(4) * t313;
t430 = -t136 - t437;
t387 = t313 * t143;
t428 = t389 + t436;
t66 = -t143 ^ 2 + t313 ^ 2;
t55 = -t143 * t253 + t71;
t282 = g(3) * t245 - t157 * t143 + t246 * t320 - t421;
t279 = t157 * t313 + t290 + t434;
t56 = -t253 * t313 + t288;
t239 = pkin(4) * t264 + pkin(10);
t425 = (qJD(6) * t239 - t432 + t433) * t91;
t247 = pkin(3) * t272 + pkin(4);
t382 = t265 * t267;
t178 = pkin(3) * t382 + t264 * t247;
t172 = pkin(10) + t178;
t417 = pkin(3) * t186;
t322 = -t417 - t432;
t50 = t322 + t433;
t424 = (qJD(6) * t172 + t50) * t91;
t323 = -pkin(2) * t379 + t272 * t248;
t176 = pkin(4) + t323;
t181 = pkin(2) * t378 + t248 * t267;
t132 = t264 * t176 + t265 * t181;
t127 = pkin(10) + t132;
t250 = pkin(2) * t366;
t423 = (qJD(6) * t127 + t250 + t50) * t91;
t37 = t271 * t39;
t422 = -t91 * t362 + t37;
t325 = -t418 * t215 - t216 * t268;
t137 = -pkin(9) * t199 + t325;
t302 = -t268 * t269 + t273 * t418;
t369 = -t268 * t215 + t418 * t216;
t138 = pkin(9) * t302 + t369;
t371 = t267 * t137 + t272 * t138;
t151 = t199 * t267 - t272 * t302;
t152 = t199 * t272 + t267 * t302;
t100 = t265 * t151 + t152 * t264;
t101 = -t151 * t264 + t152 * t265;
t6 = t264 * t15 + t265 * t17;
t344 = pkin(10) * t252 + qJD(6) * t48 + t6;
t329 = t272 * t137 - t138 * t267;
t304 = -qJ(5) * t152 + t329;
t76 = -qJ(5) * t151 + t371;
t47 = t264 * t304 + t265 * t76;
t153 = t260 * t302;
t78 = -qJD(4) * t151 + t153 * t272 - t154 * t267;
t79 = qJD(4) * t152 + t153 * t267 + t272 * t154;
t52 = -t264 * t79 + t265 * t78;
t162 = -pkin(3) * t302 - t406;
t301 = pkin(4) * t151 + t162;
t54 = pkin(5) * t100 - pkin(10) * t101 + t301;
t350 = qJD(2) * t419;
t204 = t269 * t350;
t206 = t273 * t350;
t297 = -t418 * t204 - t268 * t206 - t215 * t347 - t216 * t365;
t88 = -pkin(9) * t154 + t297;
t286 = -qJD(3) * t369 + t268 * t204 - t418 * t206;
t89 = -t153 * pkin(9) + t286;
t300 = t137 * t363 - t138 * t364 + t267 * t89 + t272 * t88;
t27 = -qJ(5) * t79 - qJD(5) * t151 + t300;
t289 = -qJD(4) * t371 - t267 * t88 + t272 * t89;
t280 = -qJ(5) * t78 - qJD(5) * t152 + t289;
t9 = t264 * t280 + t265 * t27;
t420 = -t100 * t344 + t3 * t101 + t30 * t52 - (qJD(6) * t54 + t9) * t91 - t47 * t39;
t410 = t54 * t39;
t403 = t264 * t430 - t265 * t428;
t402 = t264 * t428 + t265 * t430;
t401 = pkin(3) * qJD(4);
t398 = t266 * t84;
t396 = t271 * t84;
t395 = t271 * t91;
t392 = t30 * t101;
t332 = -t118 * t267 - t111;
t306 = t332 - t389;
t373 = t272 * t118 - t109;
t70 = t136 + t373;
t391 = -t264 * t70 + t265 * t306 + (t264 * t272 + t382) * t401;
t383 = t264 * t267;
t390 = -t264 * t306 - t265 * t70 + (t265 * t272 - t383) * t401;
t386 = t186 * t435;
t255 = cos(t263);
t385 = t255 * t270;
t384 = t255 * t274;
t381 = t266 * t270;
t380 = t266 * t274;
t368 = pkin(3) * t255 + pkin(4) * t246;
t261 = t269 ^ 2;
t367 = -t273 ^ 2 + t261;
t251 = t269 * t400;
t353 = t258 + t368;
t145 = pkin(3) * t154 + t251;
t180 = pkin(2) * t346 - qJDD(1) * t406;
t104 = pkin(3) * t121 + t180;
t284 = -pkin(4) * t288 + qJDD(5) + t104;
t10 = -pkin(5) * t42 - pkin(10) * t43 + t284;
t341 = qJD(6) * t31 - t10;
t254 = sin(t263);
t321 = -pkin(3) * t254 - pkin(4) * t245;
t319 = g(1) * t270 - g(2) * t274;
t318 = t32 * t427 + t33 * t426;
t131 = t176 * t265 - t181 * t264;
t311 = t397 * t427 + t422;
t310 = pkin(4) * t79 + t145;
t309 = g(3) * t230 - t344;
t177 = -pkin(3) * t383 + t247 * t265;
t308 = t320 * t230;
t307 = -0.2e1 * pkin(1) * t360 - pkin(7) * qJDD(2);
t35 = t265 * t64 - t399;
t295 = -t239 * t39 + t35 * t91 - t411;
t294 = -t127 * t39 - t402 * t91 - t411;
t293 = -t172 * t39 - t390 * t91 - t411;
t275 = qJD(2) ^ 2;
t292 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t275 + t319;
t276 = qJD(1) ^ 2;
t291 = pkin(1) * t276 - pkin(7) * qJDD(1) + t320;
t278 = g(1) * t384 + g(2) * t385 + g(3) * t254 + t214 * t435 - t285;
t277 = -g(3) * t255 - t214 * t186 + t254 * t320 + t287;
t256 = -qJ(5) - pkin(9) - t419;
t240 = -pkin(4) * t265 - pkin(5);
t171 = -pkin(5) - t177;
t167 = pkin(1) + t353;
t166 = t231 * t375 + t381;
t165 = -t231 * t380 + t376;
t164 = -t231 * t376 + t380;
t163 = t231 * t381 + t375;
t158 = t250 - t417;
t126 = -pkin(5) - t131;
t122 = t186 ^ 2 - t435 ^ 2;
t103 = -t317 + (-qJD(1) * t199 - t186) * t260;
t102 = t120 - t440;
t51 = t264 * t78 + t265 * t79;
t46 = t264 * t76 - t265 * t304;
t34 = t264 * t64 + t61;
t20 = pkin(5) * t51 - pkin(10) * t52 + t310;
t11 = -t397 * t91 + t37 + t409;
t8 = t264 * t27 - t265 * t280;
t7 = t271 * t10;
t1 = -t397 * t84 + t439;
t2 = [qJDD(1), t319, t320, qJDD(1) * t261 + 0.2e1 * t269 * t345, 0.2e1 * t269 * t358 - 0.2e1 * t360 * t367, qJDD(2) * t269 + t273 * t275, qJDD(2) * t273 - t269 * t275, 0, t269 * t307 + t273 * t292, -t269 * t292 + t273 * t307, t120 * t199 - t153 * t186, t120 * t302 - t121 * t199 + t153 * t435 + t154 * t186, t153 * t260 + t199 * t259, -t154 * t260 + t259 * t302, 0, g(1) * t385 - g(2) * t384 - t121 * t406 - t214 * t154 - t180 * t302 - t251 * t435 + t259 * t325 + t260 * t286, -t120 * t406 - t214 * t153 + t180 * t199 - t186 * t251 - t254 * t319 - t259 * t369 - t260 * t297, t152 * t71 - t313 * t78, t143 * t78 - t151 * t71 + t152 * t288 + t313 * t79, t152 * t252 + t253 * t78, -t151 * t252 - t253 * t79, 0, t104 * t151 - t143 * t145 + t157 * t79 - t162 * t288 + t246 * t319 + t252 * t329 + t253 * t289, t104 * t152 - t145 * t313 + t157 * t78 + t162 * t71 - t245 * t319 - t252 * t371 - t253 * t300, -t100 * t6 - t101 * t5 - t32 * t52 - t33 * t51 + t42 * t47 + t426 * t8 + t427 * t9 + t43 * t46 - t320, t6 * t47 + t33 * t9 - t5 * t46 - t32 * t8 + t284 * t301 + t105 * t310 - g(1) * (-t167 * t270 - t256 * t274) - g(2) * (t167 * t274 - t256 * t270) t52 * t396 + (-t362 * t84 + t23) * t101 (-t271 * t82 - t398) * t52 + (-t22 - t26 * t271 + (t266 * t82 - t396) * qJD(6)) * t101, t100 * t25 + t101 * t422 + t52 * t395 + t51 * t84, -t100 * t26 - t101 * t404 - t397 * t52 - t51 * t82, t100 * t39 + t51 * t91, -g(1) * t164 - g(2) * t166 + t7 * t100 + t18 * t51 + t46 * t26 + t8 * t82 + (t20 * t91 + t410 + (-t31 * t100 - t47 * t91 + t392) * qJD(6)) * t271 + t420 * t266, -g(1) * t163 - g(2) * t165 - t19 * t51 + t46 * t25 + t8 * t84 + (-(-qJD(6) * t47 + t20) * t91 - t410 + t341 * t100 - qJD(6) * t392) * t266 + t420 * t271; 0, 0, 0, -t269 * t276 * t273, t367 * t276, t359, t358, qJDD(2), -g(3) * t273 + t269 * t291, g(3) * t269 + t273 * t291, t386, t122, t102, t103, t259, -t326 * t260 + (t259 * t418 - t260 * t365 + t366 * t435) * pkin(2) + t277, t370 * t260 + (t186 * t366 - t268 * t259 - t260 * t347) * pkin(2) + t278, t387, t66, t55, t56, t252, t158 * t143 + t323 * t252 + t253 * t436 + t279, t158 * t313 - t181 * t252 + t253 * t437 + t282, -t131 * t43 + t132 * t42 + t402 * t427 + t403 * t426 + t318, t6 * t132 + t5 * t131 - t105 * (t250 + t322) - g(3) * t353 + t402 * t33 - t403 * t32 - t320 * (-pkin(2) * t269 + t321) t13, t1, t12, t11, -t407, t126 * t26 + t403 * t82 + (-t431 - t423) * t271 + t294 * t266 + t312, t126 * t25 + t403 * t84 + t294 * t271 + (-t308 + t423) * t266 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386, t122, t102, t103, t259, -t260 * t303 + t277, t260 * t327 + t278, t387, t66, t55, t56, t252, -t332 * t253 + (-t143 * t186 + t272 * t252 - t253 * t364) * pkin(3) + t279, t373 * t253 + (-t186 * t313 - t267 * t252 - t253 * t363) * pkin(3) + t282, -t177 * t43 + t178 * t42 + t390 * t427 + t391 * t426 + t318, -g(3) * t368 - t105 * t322 + t5 * t177 + t6 * t178 - t32 * t391 - t320 * t321 + t33 * t390, t13, t1, t12, t11, -t407, t171 * t26 + t391 * t82 + (-t431 - t424) * t271 + t293 * t266 + t312, t171 * t25 + t391 * t84 + t293 * t271 + (-t308 + t424) * t266 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t387, t66, t55, t56, t252, -t253 * t315 + t279, t253 * t333 + t282 (t264 * t42 - t265 * t43) * pkin(4) + (t32 - t35) * t427 + (t33 - t34) * t426, t32 * t34 - t33 * t35 + (t105 * t313 + t264 * t6 + t265 * t5 + t434) * pkin(4), t13, -t398 * t91 + t439, t12, t311 + t409, -t407, t240 * t26 - t34 * t82 + t295 * t266 + (-t431 - t425) * t271 + t312, t240 * t25 - t34 * t84 + t295 * t271 + (-t308 + t425) * t266 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426 ^ 2 - t427 ^ 2, t32 * t426 - t33 * t427 + t284 - t319, 0, 0, 0, 0, 0, t311 - t409, -t395 * t91 - t36 - t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t82, -t82 ^ 2 + t84 ^ 2, t82 * t91 + t25, t84 * t91 - t26, t39, -g(1) * t165 + g(2) * t163 + t19 * t91 + t266 * t309 - t30 * t84 - t31 * t361 + t7, g(1) * t166 - g(2) * t164 + t18 * t91 + t266 * t341 + t271 * t309 + t30 * t82;];
tau_reg  = t2;
