% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:17:55
% EndTime: 2019-03-09 16:18:23
% DurationCPUTime: 11.38s
% Computational Cost: add. (15775->704), mult. (41924->947), div. (0->0), fcn. (32537->10), ass. (0->309)
t287 = sin(qJ(3));
t289 = cos(qJ(3));
t414 = cos(pkin(6));
t354 = t414 * qJD(1);
t324 = t354 + qJD(2);
t284 = sin(pkin(6));
t288 = sin(qJ(2));
t383 = qJD(1) * t288;
t364 = t284 * t383;
t212 = t287 * t324 + t289 * t364;
t283 = sin(pkin(11));
t285 = cos(pkin(11));
t290 = cos(qJ(2));
t382 = qJD(1) * t290;
t269 = t284 * t382;
t333 = t269 - qJD(3);
t154 = t212 * t283 + t285 * t333;
t286 = sin(qJ(6));
t309 = t285 * t212 - t283 * t333;
t432 = cos(qJ(6));
t458 = -t154 * t432 + t286 * t309;
t472 = t458 ^ 2;
t459 = t286 * t154 + t309 * t432;
t471 = t459 ^ 2;
t371 = qJD(1) * qJD(2);
t356 = t284 * t371;
t341 = t290 * t356;
t210 = t287 * t364 - t289 * t324;
t467 = qJD(3) * t210;
t162 = -t289 * t341 + t467;
t342 = t288 * t356;
t128 = -t162 * t283 - t285 * t342;
t129 = -t162 * t285 + t283 * t342;
t35 = qJD(6) * t459 - t432 * t128 + t286 * t129;
t372 = qJD(6) - t210;
t466 = t372 * t459;
t470 = -t35 + t466;
t367 = pkin(1) * t414;
t347 = t290 * t367;
t232 = -pkin(8) * t364 + qJD(1) * t347;
t317 = t284 * (pkin(2) * t288 - pkin(9) * t290);
t233 = qJD(1) * t317;
t147 = t289 * t232 + t287 * t233;
t130 = qJ(4) * t364 + t147;
t276 = t288 * t367;
t335 = pkin(3) * t287 - qJ(4) * t289;
t397 = t284 * t290;
t152 = (t276 + (pkin(8) + t335) * t397) * qJD(1);
t239 = qJD(3) * t335 - qJD(4) * t287;
t452 = t283 * t130 + (-t152 + t239) * t285;
t345 = t287 * t269;
t378 = qJD(3) * t287;
t86 = t285 * t130 + t283 * t152;
t469 = -qJD(5) * t289 - t86 + (-t345 + t378) * qJ(5);
t468 = t372 * t458;
t384 = qJD(1) * t284;
t393 = t289 * t290;
t191 = (t283 * t288 + t285 * t393) * t384;
t368 = -pkin(9) * t283 - pkin(4);
t395 = t285 * t289;
t433 = pkin(4) + pkin(5);
t465 = -pkin(10) * t191 - t345 * t433 - (-pkin(10) * t395 + (-pkin(5) + t368) * t287) * qJD(3) + t452;
t399 = t283 * t289;
t190 = t269 * t399 - t285 * t364;
t224 = t283 * t239;
t396 = t285 * t287;
t464 = -pkin(10) * t190 + t224 + (-pkin(9) * t396 + pkin(10) * t399) * qJD(3) + t469;
t463 = t154 * t210;
t386 = pkin(8) * t397 + t276;
t237 = t386 * qJD(2);
t223 = qJD(1) * t237;
t349 = qJD(3) * t333;
t462 = -pkin(9) * t349 + t223;
t226 = t414 * pkin(9) + t386;
t186 = qJD(2) * pkin(9) + qJD(1) * t226;
t228 = (-pkin(2) * t290 - pkin(9) * t288 - pkin(1)) * t284;
t199 = qJD(1) * t228;
t110 = -t287 * t186 + t289 * t199;
t234 = qJD(2) * t317;
t221 = qJD(1) * t234;
t398 = t284 * t288;
t245 = -pkin(8) * t398 + t347;
t236 = t245 * qJD(2);
t222 = qJD(1) * t236;
t377 = qJD(3) * t289;
t310 = t186 * t378 - t199 * t377 - t287 * t221 - t289 * t222;
t461 = t333 * t110 - t310;
t185 = -pkin(2) * t324 - t232;
t100 = t210 * pkin(3) - t212 * qJ(4) + t185;
t111 = t289 * t186 + t287 * t199;
t104 = -qJ(4) * t333 + t111;
t50 = t100 * t285 - t283 * t104;
t336 = qJD(5) - t50;
t29 = -pkin(10) * t309 - t210 * t433 + t336;
t51 = t283 * t100 + t285 * t104;
t42 = t210 * qJ(5) + t51;
t31 = pkin(10) * t154 + t42;
t318 = t286 * t31 - t29 * t432;
t325 = t287 * t341;
t379 = qJD(3) * t212;
t163 = t325 + t379;
t63 = qJ(4) * t342 - qJD(4) * t333 - t310;
t67 = t163 * pkin(3) + t162 * qJ(4) - t212 * qJD(4) + t223;
t26 = -t283 * t63 + t285 * t67;
t8 = -pkin(10) * t129 - t163 * t433 - t26;
t27 = t283 * t67 + t285 * t63;
t16 = t163 * qJ(5) + t210 * qJD(5) + t27;
t9 = pkin(10) * t128 + t16;
t1 = -qJD(6) * t318 + t286 * t8 + t432 * t9;
t460 = t318 * t372 + t1;
t338 = t283 * t377 - t190;
t407 = t163 * t283;
t457 = t128 * t289 - t338 * t210 + t287 * (t154 * t333 - t407);
t456 = t309 ^ 2;
t280 = t284 ^ 2;
t357 = t280 * t371;
t455 = -0.2e1 * t357;
t18 = -pkin(4) * t163 - t26;
t453 = -t210 * t42 + t18;
t451 = t210 * t309;
t450 = t210 * t333;
t449 = t212 * t333;
t447 = t287 * t333;
t446 = t288 * t290;
t242 = t287 * t398 - t289 * t414;
t127 = t163 * t242;
t243 = t287 * t414 + t289 * t398;
t362 = qJD(2) * t397;
t175 = qJD(3) * t243 + t287 * t362;
t82 = t175 * t210 + t127;
t159 = t163 * t289;
t88 = -t210 * t447 - t159;
t146 = -t287 * t232 + t289 * t233;
t311 = pkin(3) * t364 + t146;
t442 = qJ(5) * t191 - qJD(5) * t396 + t311;
t358 = qJD(6) * t432;
t374 = qJD(6) * t286;
t441 = -t283 * t358 + t285 * t374;
t250 = t432 * t283 - t286 * t285;
t80 = -t186 * t377 - t199 * t378 + t289 * t221 - t287 * t222;
t440 = t333 * t111 - t80;
t300 = pkin(3) * t342 + t80;
t292 = qJ(5) * t129 + qJD(5) * t309 + t300;
t14 = -t128 * t433 + t292;
t408 = t154 * t285;
t327 = t283 * t309 + t408;
t438 = t128 * t283 - t129 * t285 + t210 * t327;
t413 = qJ(5) * t285;
t436 = t283 * t433 - t413;
t409 = t129 * t283;
t411 = t128 * t285;
t435 = -t154 * t191 - t190 * t309 + t287 * (t409 + t411) + t327 * t377;
t434 = t210 ^ 2;
t291 = qJD(1) ^ 2;
t431 = t459 * t458;
t430 = -pkin(10) + qJ(4);
t259 = -pkin(3) * t289 - qJ(4) * t287 - pkin(2);
t271 = pkin(9) * t399;
t279 = t289 * pkin(4);
t161 = pkin(5) * t289 + t271 + t279 + (-pkin(10) * t287 - t259) * t285;
t209 = pkin(9) * t395 + t283 * t259;
t187 = -qJ(5) * t289 + t209;
t400 = t283 * t287;
t165 = pkin(10) * t400 + t187;
t98 = t286 * t161 + t165 * t432;
t429 = qJD(6) * t98 + t464 * t286 + t465 * t432;
t97 = t161 * t432 - t286 * t165;
t428 = -qJD(6) * t97 + t465 * t286 - t464 * t432;
t381 = qJD(2) * t288;
t90 = -t226 * t378 + t228 * t377 + t287 * t234 + t289 * t236;
t81 = (qJ(4) * t381 - qJD(4) * t290) * t284 + t90;
t176 = -qJD(3) * t242 + t289 * t362;
t89 = t175 * pkin(3) - t176 * qJ(4) - t243 * qJD(4) + t237;
t33 = t283 * t89 + t285 * t81;
t24 = pkin(4) * t128 - t292;
t426 = t24 * t283;
t425 = t24 * t285;
t424 = t283 * t300;
t423 = t285 * t300;
t264 = t430 * t283;
t265 = t430 * t285;
t182 = t286 * t264 + t265 * t432;
t107 = t283 * t110;
t136 = pkin(3) * t212 + qJ(4) * t210;
t39 = t107 + (pkin(10) * t210 - t136) * t285 - t433 * t212;
t403 = t210 * t283;
t73 = t285 * t110 + t283 * t136;
t56 = t212 * qJ(5) + t73;
t48 = -pkin(10) * t403 + t56;
t422 = -qJD(4) * t250 + qJD(6) * t182 - t286 * t48 + t39 * t432;
t181 = t264 * t432 - t286 * t265;
t314 = -t286 * t283 - t285 * t432;
t421 = qJD(4) * t314 - qJD(6) * t181 + t286 * t39 + t432 * t48;
t370 = pkin(9) * t378;
t171 = -t285 * t370 + t224;
t420 = t171 + t469;
t308 = -pkin(9) - t436;
t419 = t190 * t433 + t308 * t377 - t442;
t418 = pkin(4) * t345 + t368 * t378 - t452;
t334 = -pkin(4) * t283 + t413;
t321 = pkin(9) - t334;
t417 = -pkin(4) * t190 + t321 * t377 + t442;
t416 = t283 * t370 + t452;
t415 = t171 - t86;
t103 = pkin(3) * t333 + qJD(4) - t110;
t299 = qJ(5) * t309 - t103;
t49 = t154 * pkin(4) - t299;
t412 = qJD(3) * t49;
t406 = t163 * t285;
t404 = t210 * t212;
t401 = t280 * t291;
t392 = qJD(4) - t49;
t225 = -pkin(2) * t414 - t245;
t124 = t242 * pkin(3) - t243 * qJ(4) + t225;
t138 = t289 * t226 + t287 * t228;
t125 = -qJ(4) * t397 + t138;
t70 = t283 * t124 + t285 * t125;
t391 = -qJ(4) * t411 - qJD(4) * t408;
t240 = t314 * qJD(6);
t344 = t432 * t377;
t361 = t286 * t377;
t390 = -t190 * t432 + t191 * t286 + t240 * t287 + t283 * t344 - t285 * t361;
t389 = t286 * t190 + t432 * t191 - t283 * t361 - t285 * t344 + t287 * t441;
t388 = t210 * t250 + t441;
t387 = t314 * t210 - t240;
t137 = -t287 * t226 + t289 * t228;
t385 = t288 ^ 2 - t290 ^ 2;
t380 = qJD(2) * t289;
t376 = qJD(5) * t283;
t373 = -qJD(4) + t103;
t369 = qJ(4) * t406;
t59 = t242 * qJ(5) + t70;
t126 = pkin(3) * t397 - t137;
t363 = t284 * t381;
t355 = qJ(5) * t283 + pkin(3);
t32 = -t283 * t81 + t285 * t89;
t353 = -t210 * t436 + t111 + t376;
t72 = t136 * t285 - t107;
t69 = t124 * t285 - t283 * t125;
t208 = t259 * t285 - t271;
t351 = t290 * t333;
t350 = t333 * t284;
t25 = t175 * qJ(5) + t242 * qJD(5) + t33;
t348 = t401 * t446;
t91 = -t226 * t377 - t228 * t378 + t289 * t234 - t287 * t236;
t340 = t284 * t291 * t414;
t339 = pkin(1) * t455;
t337 = t285 * t377 - t191;
t143 = t176 * t283 - t285 * t363;
t173 = t243 * t283 + t285 * t397;
t330 = t128 * t173 + t143 * t154;
t329 = -t154 ^ 2 - t456;
t326 = t357 * t446;
t323 = 0.2e1 * t354 + qJD(2);
t322 = qJ(4) * t129 + qJD(4) * t309;
t174 = t243 * t285 - t283 * t397;
t320 = qJ(5) * t174 - t126;
t40 = -pkin(10) * t174 - t242 * t433 - t69;
t45 = pkin(10) * t173 + t59;
t10 = -t286 * t45 + t40 * t432;
t6 = t286 * t29 + t31 * t432;
t11 = t286 * t40 + t432 * t45;
t316 = t154 * t403 - t411;
t315 = t128 + t451;
t105 = -t173 * t432 + t174 * t286;
t106 = t173 * t286 + t174 * t432;
t34 = -t286 * t128 - t432 * t129 - t154 * t358 + t309 * t374;
t306 = -t210 * t432 + t358;
t305 = -t154 * t212 + t283 * t434 - t406;
t303 = pkin(3) * t363 + t91;
t301 = pkin(1) * (-qJD(2) * t354 + t401);
t144 = t176 * t285 + t283 * t363;
t298 = t128 * t174 + t129 * t173 + t143 * t309 + t144 * t154;
t297 = t128 * t242 + t143 * t210 + t154 * t175 + t163 * t173;
t296 = t128 * t400 + t154 * t338;
t295 = t129 - t463;
t2 = -qJD(6) * t6 - t286 * t9 + t432 * t8;
t293 = qJ(5) * t144 + qJD(5) * t174 + t303;
t255 = -pkin(4) * t285 - t355;
t244 = t285 * t433 + t355;
t235 = t386 * qJD(1);
t230 = t314 * t287;
t229 = t250 * t287;
t227 = t321 * t287;
t188 = -t208 + t279;
t184 = t308 * t287;
t78 = t210 * t334 + t111;
t71 = pkin(4) * t173 - t320;
t66 = t285 * t451 + t409;
t60 = -pkin(4) * t242 - t69;
t58 = -pkin(4) * t212 - t72;
t54 = -t173 * t433 + t320;
t53 = -t212 * t309 + t285 * t434 + t407;
t47 = -qJD(6) * t105 + t143 * t286 + t144 * t432;
t46 = qJD(6) * t106 - t143 * t432 + t144 * t286;
t44 = t129 * t174 + t144 * t309;
t43 = t129 * t396 + t309 * t337;
t41 = -pkin(4) * t210 + t336;
t36 = -t154 * t433 + t299;
t30 = pkin(4) * t143 - t293;
t28 = -pkin(4) * t175 - t32;
t23 = -t129 * t289 + t337 * t210 + (-t309 * t333 + t406) * t287;
t22 = -t143 * t433 + t293;
t21 = t129 * t242 + t144 * t210 + t163 * t174 + t175 * t309;
t17 = pkin(10) * t143 + t25;
t15 = -pkin(10) * t144 - t175 * t433 - t32;
t4 = -qJD(6) * t11 + t15 * t432 - t286 * t17;
t3 = qJD(6) * t10 + t286 * t15 + t17 * t432;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t326, t385 * t455, t323 * t362, -0.2e1 * t326, -t323 * t363, 0, -t223 * t414 - t237 * t324 + t288 * t339, -t222 * t414 - t236 * t324 + t290 * t339 (t222 * t290 + t223 * t288 + (-t232 * t290 - t235 * t288) * qJD(2) + (t236 * t290 + t237 * t288 + (-t245 * t290 - t288 * t386) * qJD(2)) * qJD(1)) * t284, t222 * t386 - t223 * t245 - t232 * t237 + t235 * t236, -t162 * t243 + t176 * t212, t162 * t242 - t163 * t243 - t175 * t212 - t176 * t210, -t176 * t333 + (t162 * t290 + (qJD(1) * t243 + t212) * t381) * t284, t82, t175 * t333 + (t163 * t290 + (-qJD(1) * t242 - t210) * t381) * t284 (-t280 * t382 - t350) * t381, -t91 * t333 + t237 * t210 + t225 * t163 + t223 * t242 + t185 * t175 + (-t80 * t290 + (qJD(1) * t137 + t110) * t381) * t284, t90 * t333 + t237 * t212 - t225 * t162 + t223 * t243 + t185 * t176 + (-t310 * t290 + (-qJD(1) * t138 - t111) * t381) * t284, -t110 * t176 - t111 * t175 + t137 * t162 - t138 * t163 - t210 * t90 - t212 * t91 + t242 * t310 - t243 * t80, t110 * t91 + t111 * t90 + t137 * t80 - t138 * t310 + t185 * t237 + t223 * t225, t44, -t298, t21, t330, -t297, t82, t103 * t143 + t126 * t128 - t154 * t303 + t163 * t69 - t173 * t300 + t175 * t50 + t210 * t32 + t242 * t26, t103 * t144 + t126 * t129 - t163 * t70 - t174 * t300 - t175 * t51 - t210 * t33 - t242 * t27 - t303 * t309, -t128 * t70 - t129 * t69 - t143 * t51 - t144 * t50 - t154 * t33 - t173 * t27 - t174 * t26 - t309 * t32, -t103 * t303 - t126 * t300 + t26 * t69 + t27 * t70 + t32 * t50 + t33 * t51, t44, t21, t298, t82, t297, t330, t128 * t71 + t143 * t49 + t154 * t30 - t163 * t60 + t173 * t24 - t175 * t41 - t18 * t242 - t210 * t28, -t128 * t59 + t129 * t60 - t143 * t42 + t144 * t41 - t154 * t25 - t16 * t173 + t174 * t18 + t28 * t309, -t129 * t71 - t144 * t49 + t16 * t242 + t163 * t59 - t174 * t24 + t175 * t42 + t210 * t25 - t30 * t309, t16 * t59 + t18 * t60 + t24 * t71 + t25 * t42 + t28 * t41 + t30 * t49, -t106 * t34 + t459 * t47, t105 * t34 - t106 * t35 - t458 * t47 - t459 * t46, -t106 * t163 - t175 * t459 + t242 * t34 + t372 * t47, t105 * t35 + t458 * t46, t105 * t163 + t175 * t458 + t242 * t35 - t372 * t46, -t175 * t372 + t127, -t10 * t163 + t105 * t14 + t175 * t318 - t2 * t242 + t22 * t458 + t35 * t54 + t36 * t46 + t372 * t4, t1 * t242 + t106 * t14 + t11 * t163 + t175 * t6 + t22 * t459 - t3 * t372 - t34 * t54 + t36 * t47, -t1 * t105 + t10 * t34 - t106 * t2 - t11 * t35 - t3 * t458 + t318 * t47 - t4 * t459 - t46 * t6, t1 * t11 + t10 * t2 + t14 * t54 + t22 * t36 + t3 * t6 - t318 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t348, t385 * t401, -t290 * t340, t348, t288 * t340, 0, -pkin(8) * t341 + t235 * t324 + t288 * t301, pkin(8) * t342 + t232 * t324 + t290 * t301, 0, 0, -t162 * t287 - t289 * t449 (-t162 + t450) * t289 + (-t163 + t449) * t287, -t289 * t349 + (t289 * t351 + (t287 * qJD(2) - t212) * t288) * t384, t88, t287 * t349 + (-t287 * t351 + (t210 + t380) * t288) * t384, t350 * t383, -pkin(2) * t163 + t185 * t378 + t146 * t333 - t235 * t210 - t462 * t289 + (-t110 * t288 + (-pkin(9) * t381 - t185 * t290) * t287) * t384, pkin(2) * t162 + t185 * t377 - t147 * t333 - t235 * t212 + t462 * t287 + (-t185 * t393 + (-pkin(9) * t380 + t111) * t288) * t384, t146 * t212 + t147 * t210 + ((-t163 + t379) * pkin(9) + t461) * t289 + ((-t162 + t467) * pkin(9) + t440) * t287, -pkin(2) * t223 - t110 * t146 - t111 * t147 - t185 * t235 + (-t287 * t80 - t289 * t310 + (-t110 * t289 - t111 * t287) * qJD(3)) * pkin(9), t43, -t435, t23, t296, t457, t88, -t103 * t190 + t311 * t154 + t163 * t208 + t416 * t210 + (-t26 + (pkin(9) * t154 + t103 * t283) * qJD(3)) * t289 + (pkin(9) * t128 - t333 * t50 - t424) * t287, -t103 * t191 + t311 * t309 - t163 * t209 - t415 * t210 + (t27 + (pkin(9) * t309 + t103 * t285) * qJD(3)) * t289 + (pkin(9) * t129 + t333 * t51 - t423) * t287, -t128 * t209 - t129 * t208 + t190 * t51 + t191 * t50 + (-t26 * t285 - t27 * t283) * t287 - t416 * t309 - t415 * t154 + (-t283 * t51 - t285 * t50) * t377, t103 * t311 + t208 * t26 + t209 * t27 + t415 * t51 + t416 * t50 + (t103 * t377 - t287 * t300) * pkin(9), t43, t23, t435, t88, -t457, t296, t128 * t227 - t163 * t188 - t190 * t49 + (t283 * t412 + t18) * t289 - t418 * t210 + t417 * t154 + (t333 * t41 + t426) * t287, -t128 * t187 + t129 * t188 + t190 * t42 - t191 * t41 + (-t16 * t283 + t18 * t285) * t287 + t418 * t309 - t420 * t154 + (-t283 * t42 + t285 * t41) * t377, -t129 * t227 + t163 * t187 + t191 * t49 + (-t285 * t412 - t16) * t289 + t420 * t210 - t417 * t309 + (-t333 * t42 - t425) * t287, t16 * t187 + t18 * t188 + t227 * t24 + t41 * t418 + t417 * t49 + t42 * t420, t230 * t34 - t389 * t459, -t229 * t34 + t230 * t35 + t389 * t458 + t390 * t459, t163 * t230 - t289 * t34 - t372 * t389 + t447 * t459, -t229 * t35 - t390 * t458, -t163 * t229 - t289 * t35 + t372 * t390 - t447 * t458, t372 * t447 - t159, -t14 * t229 - t163 * t97 + t184 * t35 + t2 * t289 - t318 * t447 - t36 * t390 - t372 * t429 + t419 * t458, -t1 * t289 - t14 * t230 + t163 * t98 - t184 * t34 - t36 * t389 + t372 * t428 + t419 * t459 - t447 * t6, t1 * t229 + t2 * t230 - t318 * t389 + t34 * t97 - t35 * t98 + t390 * t6 + t428 * t458 + t429 * t459, t1 * t98 + t14 * t184 + t2 * t97 + t318 * t429 + t36 * t419 - t428 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t404, t212 ^ 2 - t434, -t162 - t450, -t404, -t212 * t269 - t325, t342, -t185 * t212 - t440, t185 * t210 - t461, 0, 0, t66, -t438, t53, t316, -t305, -t404, -qJ(4) * t407 - pkin(3) * t128 - t111 * t154 - t212 * t50 + t423 + (t283 * t373 - t72) * t210, -t369 - pkin(3) * t129 - t111 * t309 + t212 * t51 - t424 + (t285 * t373 + t73) * t210, t154 * t73 + t309 * t72 + (-t210 * t50 + t27) * t285 + (-t210 * t51 - t26 + t322) * t283 + t391, pkin(3) * t300 - t103 * t111 - t50 * t72 - t51 * t73 + (-t283 * t50 + t285 * t51) * qJD(4) + (-t26 * t283 + t27 * t285) * qJ(4), t66, t53, t438, -t404, t305, t316, t128 * t255 - t154 * t78 + t210 * t58 + t212 * t41 - t425 + (-qJ(4) * t163 - qJD(5) * t154 - t210 * t392) * t283, t154 * t56 - t309 * t58 + (t210 * t41 + t16) * t285 + (t322 + t453) * t283 + t391, t369 - t129 * t255 - t212 * t42 - t426 + (t78 + t376) * t309 + (t285 * t392 - t56) * t210, t24 * t255 - t41 * t58 - t42 * t56 - t49 * t78 + (qJ(4) * t16 + qJD(4) * t42) * t285 + (qJ(4) * t18 + qJD(4) * t41 - qJD(5) * t49) * t283, -t250 * t34 - t387 * t459, -t250 * t35 - t314 * t34 + t387 * t458 + t388 * t459, -t163 * t250 + t212 * t459 - t372 * t387, -t314 * t35 - t388 * t458, -t163 * t314 - t212 * t458 + t372 * t388, t372 * t212, -t14 * t314 - t163 * t181 - t212 * t318 + t244 * t35 + t353 * t458 - t36 * t388 - t372 * t422, t14 * t250 + t163 * t182 - t212 * t6 - t244 * t34 + t353 * t459 - t36 * t387 + t372 * t421, t1 * t314 + t181 * t34 - t182 * t35 - t2 * t250 - t318 * t387 + t388 * t6 + t421 * t458 + t422 * t459, t1 * t182 + t14 * t244 + t181 * t2 + t318 * t422 + t353 * t36 - t421 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, t295, t329, t154 * t51 + t309 * t50 - t300, 0, 0, 0, 0, 0, 0, t315, t329, -t295, t154 * t42 - t309 * t41 + t24, 0, 0, 0, 0, 0, 0, -t35 - t466, t34 + t468, t471 + t472, t318 * t459 - t458 * t6 - t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154 * t309 - t163, t129 + t463, -t434 - t456, t309 * t49 + t453, 0, 0, 0, 0, 0, 0, -t286 * t372 ^ 2 - t163 * t432 - t309 * t458, t286 * t163 - t306 * t372 - t309 * t459, t470 * t286 - t306 * t458 + t432 * t34, t2 * t432 + t286 * t460 + t306 * t6 - t36 * t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t431, t471 - t472, -t34 + t468, -t431, t470, -t163, -t36 * t459 + t372 * t6 + t2, t36 * t458 - t460, 0, 0;];
tauc_reg  = t5;