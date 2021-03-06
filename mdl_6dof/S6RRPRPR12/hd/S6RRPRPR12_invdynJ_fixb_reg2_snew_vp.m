% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 16:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRPR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:19:06
% EndTime: 2019-05-06 16:19:31
% DurationCPUTime: 11.22s
% Computational Cost: add. (75515->587), mult. (174185->842), div. (0->0), fcn. (131713->12), ass. (0->360)
t466 = 2 * qJD(5);
t331 = sin(pkin(11));
t334 = cos(pkin(6));
t327 = qJD(1) * t334 + qJD(2);
t336 = sin(qJ(4));
t340 = cos(qJ(4));
t332 = sin(pkin(6));
t341 = cos(qJ(2));
t401 = qJD(1) * t341;
t385 = t332 * t401;
t282 = -t336 * t327 - t340 * t385;
t283 = t327 * t340 - t336 * t385;
t333 = cos(pkin(11));
t257 = -t333 * t282 + t283 * t331;
t259 = t331 * t282 + t333 * t283;
t212 = t259 * t257;
t337 = sin(qJ(2));
t296 = (qJD(2) * t401 + qJDD(1) * t337) * t332;
t288 = qJDD(4) + t296;
t453 = -t212 + t288;
t465 = t331 * t453;
t464 = t333 * t453;
t335 = sin(qJ(6));
t402 = qJD(1) * t337;
t386 = t332 * t402;
t316 = qJD(4) + t386;
t339 = cos(qJ(6));
t230 = t259 * t335 - t339 * t316;
t232 = t259 * t339 + t316 * t335;
t189 = t232 * t230;
t318 = qJD(2) * t386;
t396 = qJDD(1) * t341;
t297 = t332 * t396 - t318;
t326 = qJDD(1) * t334 + qJDD(2);
t378 = t340 * t297 + t326 * t336;
t243 = -qJD(4) * t283 - t378;
t244 = t282 * qJD(4) - t336 * t297 + t340 * t326;
t380 = -t333 * t243 + t244 * t331;
t199 = qJDD(6) + t380;
t454 = -t189 + t199;
t463 = t335 * t454;
t261 = t282 * t283;
t452 = t261 + t288;
t462 = t336 * t452;
t461 = t339 * t454;
t460 = t340 * t452;
t422 = t259 * t316;
t174 = t380 + t422;
t323 = t327 ^ 2;
t330 = t332 ^ 2;
t343 = qJD(1) ^ 2;
t409 = t330 * t343;
t446 = t341 ^ 2;
t325 = t446 * t409;
t300 = -t325 - t323;
t313 = t341 * t337 * t409;
t450 = -t326 - t313;
t459 = pkin(1) * (t300 * t337 - t341 * t450);
t324 = t337 ^ 2 * t409;
t277 = -t324 - t323;
t295 = -t313 + t326;
t458 = pkin(1) * (t277 * t341 - t295 * t337);
t417 = t450 * t337;
t457 = pkin(8) * (-t300 * t341 - t417);
t416 = t295 * t341;
t456 = pkin(8) * (t277 * t337 + t416);
t208 = pkin(5) * t257 - pkin(10) * t259;
t447 = t316 ^ 2;
t338 = sin(qJ(1));
t342 = cos(qJ(1));
t383 = t338 * g(1) - g(2) * t342;
t442 = pkin(8) * t332;
t290 = qJDD(1) * pkin(1) + t343 * t442 + t383;
t293 = pkin(3) * t386 - pkin(9) * t327;
t384 = qJD(3) * t402;
t317 = -0.2e1 * t332 * t384;
t303 = t327 * pkin(2) * t386;
t439 = t334 * g(3);
t354 = -t296 * qJ(3) + t303 - t439;
t410 = t327 * t341;
t393 = qJ(3) * t410;
t444 = -pkin(2) - pkin(9);
t182 = -pkin(3) * t325 + t317 + t444 * t297 + (-t290 + (-t293 * t337 - t393) * qJD(1)) * t332 + t354;
t406 = t337 * qJ(3);
t370 = -pkin(2) * t341 - t406;
t403 = qJD(1) * t332;
t292 = t370 * t403;
t373 = g(1) * t342 + g(2) * t338;
t291 = -pkin(1) * t343 + qJDD(1) * t442 - t373;
t418 = t290 * t334;
t379 = t291 * t337 - t341 * t418;
t349 = -t326 * pkin(2) - t323 * qJ(3) + qJDD(3) + t379;
t441 = g(3) * t341;
t185 = pkin(3) * t296 + t450 * pkin(9) + (t441 + (-pkin(3) * t410 + t292 * t337) * qJD(1)) * t332 + t349;
t143 = t182 * t340 + t185 * t336;
t263 = pkin(4) * t316 - qJ(5) * t283;
t280 = t282 ^ 2;
t109 = -pkin(4) * t280 + qJ(5) * t243 - t263 * t316 + t143;
t142 = t182 * t336 - t340 * t185;
t419 = t282 * t316;
t222 = -t244 + t419;
t345 = pkin(4) * t452 + qJ(5) * t222 - t142;
t60 = -0.2e1 * qJD(5) * t257 + t333 * t109 + t331 * t345;
t56 = -pkin(5) * t447 + pkin(10) * t288 - t208 * t257 + t60;
t321 = t323 * pkin(2);
t408 = t332 * t337;
t372 = -g(3) * t408 + t337 * t418;
t346 = (t292 * t403 + t291) * t341 - t321 + t372;
t411 = t326 * qJ(3);
t445 = 0.2e1 * qJD(3);
t181 = t411 + t297 * pkin(3) - pkin(9) * t325 + (t445 + t293) * t327 + t346;
t144 = -t243 * pkin(4) - t280 * qJ(5) + t263 * t283 + qJDD(5) + t181;
t202 = t331 * t243 + t333 * t244;
t423 = t257 * t316;
t371 = t202 - t423;
t90 = pkin(5) * t174 - t371 * pkin(10) + t144;
t38 = t335 * t56 - t339 * t90;
t39 = t335 * t90 + t339 * t56;
t25 = t335 * t38 + t339 * t39;
t455 = t330 * (qJD(1) * t327 - t334 * t343);
t305 = t327 * t385;
t451 = t296 + t305;
t382 = t109 * t331 - t333 * t345;
t59 = t259 * t466 + t382;
t254 = qJD(6) + t257;
t381 = t202 * t335 - t339 * t288;
t136 = (qJD(6) - t254) * t232 + t381;
t387 = t327 * t402;
t268 = -t318 + (t387 + t396) * t332;
t449 = t332 * ((t325 - t323) * t337 + t416) + t268 * t334;
t267 = -t305 + t296;
t448 = t332 * ((-t324 + t323) * t341 - t417) + t267 * t334;
t228 = t230 ^ 2;
t229 = t232 ^ 2;
t253 = t254 ^ 2;
t255 = t257 ^ 2;
t256 = t259 ^ 2;
t281 = t283 ^ 2;
t443 = pkin(5) * t331;
t440 = t297 * pkin(2);
t55 = -t288 * pkin(5) - t447 * pkin(10) + (t466 + t208) * t259 + t382;
t52 = t335 * t55;
t34 = t331 * t60 - t333 * t59;
t438 = t336 * t34;
t53 = t339 * t55;
t437 = t34 * t340;
t436 = t144 * t331;
t435 = t144 * t333;
t147 = t189 + t199;
t434 = t147 * t335;
t433 = t147 * t339;
t348 = (-qJD(4) + t316) * t283 - t378;
t171 = t222 * t340 + t336 * t348;
t432 = t171 * t337;
t431 = t181 * t336;
t430 = t181 * t340;
t197 = t212 + t288;
t429 = t197 * t331;
t428 = t197 * t333;
t235 = -t261 + t288;
t427 = t235 * t336;
t426 = t235 * t340;
t425 = t254 * t335;
t424 = t254 * t339;
t415 = t316 * t331;
t414 = t316 * t333;
t413 = t316 * t336;
t412 = t316 * t340;
t407 = t332 * t341;
t302 = -t324 - t325;
t404 = pkin(1) * (-t302 * t332 + (-t267 * t341 + t268 * t337) * t334) + (t267 * t337 + t268 * t341) * t442;
t397 = qJD(6) + t254;
t16 = t25 * t331 - t333 * t55;
t395 = pkin(4) * t16 - pkin(5) * t55 + pkin(10) * t25;
t394 = t327 * t445;
t392 = t331 * t189;
t391 = t333 * t189;
t390 = t337 * t212;
t389 = t337 * t261;
t388 = -pkin(5) * t333 - pkin(4);
t35 = t331 * t59 + t333 * t60;
t180 = -t229 - t253;
t112 = -t180 * t335 - t433;
t367 = -t202 * t339 - t288 * t335;
t141 = t230 * t397 + t367;
t72 = t112 * t331 + t141 * t333;
t377 = pkin(4) * t72 + pkin(5) * t141 + pkin(10) * t112 + t52;
t166 = -t253 - t228;
t107 = t166 * t339 - t463;
t137 = -t232 * t397 - t381;
t69 = t107 * t331 + t137 * t333;
t376 = pkin(4) * t69 + pkin(5) * t137 + pkin(10) * t107 - t53;
t225 = -t256 - t447;
t160 = t225 * t333 - t429;
t375 = pkin(4) * t160 - t60;
t158 = t228 + t229;
t155 = -qJD(6) * t230 - t367;
t200 = t254 * t230;
t140 = t155 + t200;
t88 = -t136 * t339 + t140 * t335;
t64 = t158 * t333 + t331 * t88;
t374 = pkin(4) * t64 + pkin(5) * t158 + pkin(10) * t88 + t25;
t270 = t332 * t290 + t439;
t24 = t335 * t39 - t339 * t38;
t91 = -t142 * t340 + t143 * t336;
t368 = t142 * t336 + t143 * t340;
t359 = pkin(1) - t370;
t203 = -t447 - t255;
t152 = t203 * t331 + t464;
t358 = pkin(4) * t152 - t59;
t357 = t244 + t419;
t356 = -qJD(1) * t393 - t290;
t353 = t292 * t402 + t441;
t352 = -t380 + t422;
t248 = g(3) * t407 + t379;
t249 = t291 * t341 + t372;
t351 = (t248 * t337 + t249 * t341) * t332;
t350 = t341 * t444 - pkin(1) - t406;
t347 = -t354 + t440;
t214 = t332 * t353 + t349;
t344 = t346 + t394;
t308 = t334 * t326;
t301 = t324 - t325;
t271 = t288 * t408;
t269 = -t318 + (-t387 + t396) * t332;
t265 = -t281 + t447;
t264 = t280 - t447;
t262 = -t281 - t447;
t260 = t281 - t280;
t247 = -t447 - t280;
t246 = (t296 * t332 + t341 * t455) * t337;
t245 = (t297 * t332 - t337 * t455) * t341;
t238 = -t256 + t447;
t237 = t255 - t447;
t233 = -t280 - t281;
t227 = -pkin(2) * t267 + qJ(3) * t268;
t226 = (t282 * t340 + t283 * t336) * t316;
t218 = (qJD(4) + t316) * t283 + t378;
t217 = t244 * t340 - t283 * t413;
t216 = -t243 * t336 - t282 * t412;
t215 = t334 * t301 + (t337 * t269 + t341 * t451) * t332;
t211 = t264 * t340 - t427;
t210 = -t265 * t336 + t460;
t209 = t256 - t255;
t206 = t344 + t411;
t204 = t262 * t340 - t427;
t195 = -t229 + t253;
t194 = t228 - t253;
t192 = t247 * t336 + t460;
t191 = (-t257 * t333 + t259 * t331) * t316;
t190 = (-t257 * t331 - t259 * t333) * t316;
t188 = pkin(2) * t450 - qJ(3) * t300 + t214;
t187 = t229 - t228;
t186 = -pkin(2) * t277 + (t295 + t326) * qJ(3) + t344;
t184 = -t255 - t256;
t178 = -t202 - t423;
t172 = -t218 * t340 - t336 * t357;
t170 = t202 * t333 - t259 * t415;
t169 = t202 * t331 + t259 * t414;
t168 = t257 * t414 + t331 * t380;
t167 = t257 * t415 - t333 * t380;
t165 = t237 * t333 - t429;
t164 = -t238 * t331 + t464;
t163 = t237 * t331 + t428;
t162 = t238 * t333 + t465;
t161 = -t225 * t331 - t428;
t156 = -pkin(2) * t214 + qJ(3) * t206;
t154 = -qJD(6) * t232 - t381;
t153 = t203 * t333 - t465;
t150 = (-t230 * t339 + t232 * t335) * t254;
t149 = (-t230 * t335 - t232 * t339) * t254;
t145 = -t190 * t336 + t191 * t340;
t139 = t155 - t200;
t133 = t155 * t339 - t232 * t425;
t132 = t155 * t335 + t232 * t424;
t131 = -t154 * t335 + t230 * t424;
t130 = t154 * t339 + t230 * t425;
t129 = -t178 * t331 + t333 * t352;
t128 = -t174 * t333 - t331 * t371;
t127 = t178 * t333 + t331 * t352;
t126 = -t174 * t331 + t333 * t371;
t125 = pkin(4) * t127;
t124 = t150 * t333 + t199 * t331;
t123 = t150 * t331 - t199 * t333;
t122 = -t169 * t336 + t170 * t340;
t121 = -t167 * t336 + t168 * t340;
t120 = t194 * t339 - t434;
t119 = -t195 * t335 + t461;
t118 = t194 * t335 + t433;
t117 = t195 * t339 + t463;
t116 = -t163 * t336 + t165 * t340;
t115 = -t162 * t336 + t164 * t340;
t113 = t160 * t340 + t161 * t336;
t111 = t180 * t339 - t434;
t106 = t166 * t335 + t461;
t102 = qJ(3) * t357 + t204 * t444 + t430;
t101 = t133 * t333 + t392;
t100 = t131 * t333 - t392;
t99 = t133 * t331 - t391;
t98 = t131 * t331 + t391;
t96 = t152 * t340 + t153 * t336;
t95 = qJ(3) * t218 + t192 * t444 + t431;
t94 = -qJ(5) * t160 + t435;
t93 = -qJ(5) * t152 + t436;
t87 = t137 * t339 - t139 * t335;
t86 = -t136 * t335 - t140 * t339;
t85 = t137 * t335 + t139 * t339;
t83 = -pkin(4) * t371 + qJ(5) * t161 + t436;
t82 = -pkin(4) * t174 + qJ(5) * t153 - t435;
t81 = t120 * t333 - t136 * t331;
t80 = t119 * t333 + t140 * t331;
t79 = t120 * t331 + t136 * t333;
t78 = t119 * t331 - t140 * t333;
t76 = -t126 * t336 + t128 * t340;
t75 = t127 * t340 + t129 * t336;
t74 = -t123 * t336 + t124 * t340;
t73 = t112 * t333 - t141 * t331;
t70 = t107 * t333 - t137 * t331;
t67 = t187 * t331 + t333 * t87;
t66 = -t187 * t333 + t331 * t87;
t65 = -t158 * t331 + t333 * t88;
t62 = t101 * t340 - t336 * t99;
t61 = t100 * t340 - t336 * t98;
t57 = qJ(3) * t233 + t171 * t444 - t91;
t51 = qJ(3) * t181 + t444 * t91;
t50 = -t336 * t79 + t340 * t81;
t49 = -t336 * t78 + t340 * t80;
t47 = t336 * t73 + t340 * t72;
t45 = t336 * t70 + t340 * t69;
t44 = -pkin(10) * t111 + t53;
t43 = -pkin(10) * t106 + t52;
t42 = -t336 * t66 + t340 * t67;
t40 = t336 * t65 + t340 * t64;
t33 = pkin(4) * t34;
t32 = qJ(3) * t371 + t113 * t444 - t336 * t83 + t340 * t94;
t31 = -pkin(4) * t144 + qJ(5) * t35;
t30 = -pkin(5) * t111 + t39;
t29 = -pkin(5) * t106 + t38;
t28 = -qJ(5) * t127 - t34;
t27 = qJ(3) * t174 - t336 * t82 + t340 * t93 + t444 * t96;
t26 = -pkin(4) * t184 + qJ(5) * t129 + t35;
t21 = t336 * t35 + t437;
t20 = -pkin(10) * t86 - t24;
t19 = -qJ(5) * t72 - t30 * t331 + t333 * t44;
t18 = -qJ(5) * t69 - t29 * t331 + t333 * t43;
t17 = t25 * t333 + t331 * t55;
t14 = -pkin(4) * t111 + qJ(5) * t73 + t30 * t333 + t331 * t44;
t13 = -pkin(4) * t106 + qJ(5) * t70 + t29 * t333 + t331 * t43;
t12 = -qJ(5) * t64 + t20 * t333 + t443 * t86;
t11 = qJ(3) * t184 - t26 * t336 + t28 * t340 + t444 * t75;
t10 = qJ(5) * t65 + t20 * t331 + t388 * t86;
t8 = t16 * t340 + t17 * t336;
t7 = -qJ(5) * t16 + (-pkin(10) * t333 + t443) * t24;
t6 = qJ(3) * t144 - qJ(5) * t437 + t21 * t444 - t31 * t336;
t5 = qJ(3) * t111 - t14 * t336 + t19 * t340 + t444 * t47;
t4 = qJ(3) * t106 - t13 * t336 + t18 * t340 + t444 * t45;
t3 = qJ(5) * t17 + (-pkin(10) * t331 + t388) * t24;
t2 = qJ(3) * t86 - t10 * t336 + t12 * t340 + t40 * t444;
t1 = qJ(3) * t24 - t3 * t336 + t340 * t7 + t444 * t8;
t9 = [0, 0, 0, 0, 0, qJDD(1), t383, t373, 0, 0, t246, t215, t448, t245, t449, t308, (-t248 + t459) * t334 + (pkin(1) * t269 + t341 * t270 - t457) * t332, (-t249 + t458) * t334 + (-pkin(1) * t451 - t337 * t270 - t456) * t332, t351 + t404, pkin(1) * (t270 * t332 + (-t248 * t341 + t249 * t337) * t334) + pkin(8) * t351, t308, -t448, -t449, t246, t215, t245, t334 * t227 + (t337 * (-qJ(3) * t302 + t349) + t341 * (-pkin(2) * t302 + t249 - t321 + t394 + t411) + (qJD(1) * t292 * t446 + t337 * t353) * t332) * t332 + t404, (t188 - t459) * t334 + (t341 * (t317 - t347) + t457 + t356 * t407 - t359 * t269) * t332, (t186 - t458) * t334 + (t337 * t347 + t456 + (-t356 + 0.2e1 * t384) * t408 + t359 * t451) * t332, (t156 + pkin(1) * (t206 * t337 - t214 * t341)) * t334 + (pkin(8) * (t206 * t341 + t214 * t337) - t359 * (-qJ(3) * t451 - t270 + t303 + t317 - t440)) * t332, t334 * t217 + (-t389 + t341 * (-t244 * t336 - t283 * t412)) * t332, t334 * t172 + (t337 * t260 + t341 * (t218 * t336 - t340 * t357)) * t332, t334 * t210 + (-t337 * t222 + t341 * (-t265 * t340 - t462)) * t332, t334 * t216 + (t389 + t341 * (-t243 * t340 + t282 * t413)) * t332, t334 * t211 + (t337 * t348 + t341 * (-t264 * t336 - t426)) * t332, t334 * t226 + t271 + (-t282 * t336 + t283 * t340) * t316 * t407, (t95 + pkin(1) * (-t192 * t341 + t218 * t337)) * t334 + (t337 * (pkin(3) * t192 - t142) + t341 * (pkin(3) * t218 + t430) + pkin(8) * (t192 * t337 + t218 * t341) + t350 * (t247 * t340 - t462)) * t332, (t102 + pkin(1) * (-t204 * t341 + t337 * t357)) * t334 + (t337 * (pkin(3) * t204 - t143) + t341 * (pkin(3) * t357 - t431) + pkin(8) * (t204 * t337 + t341 * t357) + t350 * (-t262 * t336 - t426)) * t332, (t57 + pkin(1) * (-t171 * t341 + t233 * t337)) * t334 + (pkin(3) * t432 + t341 * (pkin(3) * t233 - t368) + pkin(8) * (t233 * t341 + t432) + t350 * (-t222 * t336 + t340 * t348)) * t332, (t51 + pkin(1) * (t181 * t337 - t341 * t91)) * t334 + (t350 * t368 + (pkin(3) + pkin(8)) * (t181 * t341 + t337 * t91)) * t332, t334 * t122 + (t390 + t341 * (-t169 * t340 - t170 * t336)) * t332, t334 * t76 + (t337 * t209 + t341 * (-t126 * t340 - t128 * t336)) * t332, t334 * t115 + (-t337 * t178 + t341 * (-t162 * t340 - t164 * t336)) * t332, t334 * t121 + (-t390 + t341 * (-t167 * t340 - t168 * t336)) * t332, t334 * t116 + (t337 * t352 + t341 * (-t163 * t340 - t165 * t336)) * t332, t271 + (-t190 * t340 - t191 * t336) * t407 + t334 * t145, (t27 + pkin(1) * (t174 * t337 - t341 * t96)) * t334 + (t337 * (pkin(3) * t96 + t358) + t341 * (pkin(3) * t174 - t336 * t93 - t340 * t82) + pkin(8) * (t174 * t341 + t337 * t96) + t350 * (-t152 * t336 + t153 * t340)) * t332, (t32 + pkin(1) * (-t113 * t341 + t337 * t371)) * t334 + (t337 * (pkin(3) * t113 + t375) + t341 * (pkin(3) * t371 - t336 * t94 - t340 * t83) + pkin(8) * (t113 * t337 + t341 * t371) + t350 * (-t160 * t336 + t161 * t340)) * t332, (t11 + pkin(1) * (t184 * t337 - t341 * t75)) * t334 + (t337 * (pkin(3) * t75 + t125) + t341 * (pkin(3) * t184 - t26 * t340 - t28 * t336) + pkin(8) * (t184 * t341 + t337 * t75) + t350 * (-t127 * t336 + t129 * t340)) * t332, (t6 + pkin(1) * (t144 * t337 - t21 * t341)) * t334 + (t337 * (pkin(3) * t21 + t33) + t341 * (pkin(3) * t144 + qJ(5) * t438 - t31 * t340) + pkin(8) * (t144 * t341 + t21 * t337) + t350 * (t340 * t35 - t438)) * t332, t334 * t62 + (t337 * t132 + t341 * (-t101 * t336 - t340 * t99)) * t332, t334 * t42 + (t337 * t85 + t341 * (-t336 * t67 - t340 * t66)) * t332, t334 * t49 + (t337 * t117 + t341 * (-t336 * t80 - t340 * t78)) * t332, t334 * t61 + (t337 * t130 + t341 * (-t100 * t336 - t340 * t98)) * t332, t334 * t50 + (t337 * t118 + t341 * (-t336 * t81 - t340 * t79)) * t332, t334 * t74 + (t337 * t149 + t341 * (-t123 * t340 - t124 * t336)) * t332, (t4 + pkin(1) * (t106 * t337 - t341 * t45)) * t334 + (t337 * (pkin(3) * t45 + t376) + t341 * (pkin(3) * t106 - t13 * t340 - t18 * t336) + pkin(8) * (t106 * t341 + t337 * t45) + t350 * (-t336 * t69 + t340 * t70)) * t332, (t5 + pkin(1) * (t111 * t337 - t341 * t47)) * t334 + (t337 * (pkin(3) * t47 + t377) + t341 * (pkin(3) * t111 - t14 * t340 - t19 * t336) + pkin(8) * (t111 * t341 + t337 * t47) + t350 * (-t336 * t72 + t340 * t73)) * t332, (t2 + pkin(1) * (t337 * t86 - t341 * t40)) * t334 + (t337 * (pkin(3) * t40 + t374) + t341 * (pkin(3) * t86 - t10 * t340 - t12 * t336) + pkin(8) * (t337 * t40 + t341 * t86) + t350 * (-t336 * t64 + t340 * t65)) * t332, (t1 + pkin(1) * (t24 * t337 - t341 * t8)) * t334 + (t337 * (pkin(3) * t8 + t395) + t341 * (pkin(3) * t24 - t3 * t340 - t336 * t7) + pkin(8) * (t24 * t341 + t337 * t8) + t350 * (-t16 * t336 + t17 * t340)) * t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, t301, t267, t313, t268, t326, -t248, -t249, 0, 0, t326, -t267, -t268, -t313, t301, t313, t227, t188, t186, t156, t217, t172, t210, t216, t211, t226, t95, t102, t57, t51, t122, t76, t115, t121, t116, t145, t27, t32, t11, t6, t62, t42, t49, t61, t50, t74, t4, t5, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, -t450, t277, t214, 0, 0, 0, 0, 0, 0, t192, t204, t171, t91, 0, 0, 0, 0, 0, 0, t96, t113, t75, t21, 0, 0, 0, 0, 0, 0, t45, t47, t40, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t260, -t222, t261, t348, t288, -t142, -t143, 0, 0, t212, t209, -t178, -t212, t352, t288, t358, t375, t125, t33, t132, t85, t117, t130, t118, t149, t376, t377, t374, t395; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t371, t184, t144, 0, 0, 0, 0, 0, 0, t106, t111, t86, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t187, t140, -t189, -t136, t199, -t38, -t39, 0, 0;];
tauJ_reg  = t9;
