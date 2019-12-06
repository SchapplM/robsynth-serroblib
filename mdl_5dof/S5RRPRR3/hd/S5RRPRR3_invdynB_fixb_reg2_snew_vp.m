% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:35
% EndTime: 2019-12-05 18:30:46
% DurationCPUTime: 8.86s
% Computational Cost: add. (36284->383), mult. (48058->557), div. (0->0), fcn. (28562->10), ass. (0->257)
t422 = qJD(1) + qJD(2);
t419 = qJD(4) + t422;
t417 = t419 ^ 2;
t433 = cos(qJ(4));
t421 = qJDD(1) + qJDD(2);
t418 = qJDD(4) + t421;
t429 = sin(qJ(4));
t461 = t429 * t418;
t378 = t433 * t417 + t461;
t455 = t433 * t418;
t381 = t429 * t417 - t455;
t426 = sin(pkin(9));
t427 = cos(pkin(9));
t321 = t427 * t378 - t426 * t381;
t425 = g(1) - qJDD(3);
t359 = pkin(7) * t378 - t433 * t425;
t488 = pkin(7) * t381 - t429 * t425;
t271 = qJ(3) * t321 + t427 * t359 - t426 * t488;
t325 = t426 * t378 + t427 * t381;
t430 = sin(qJ(2));
t434 = cos(qJ(2));
t280 = t434 * t321 - t430 * t325;
t501 = qJ(3) * t325 + t426 * t359 + t427 * t488;
t215 = pkin(6) * t280 + t434 * t271 - t430 * t501;
t431 = sin(qJ(1));
t435 = cos(qJ(1));
t284 = t430 * t321 + t434 * t325;
t515 = pkin(6) * t284 + t430 * t271 + t434 * t501;
t517 = t431 * t280 + t435 * t284;
t524 = pkin(5) * t517 + t431 * t215 + t435 * t515;
t238 = t435 * t280 - t431 * t284;
t523 = pkin(5) * t238 + t435 * t215 - t431 * t515;
t410 = t435 * g(2) + t431 * g(3);
t399 = qJDD(1) * pkin(1) + t410;
t409 = t431 * g(2) - t435 * g(3);
t437 = qJD(1) ^ 2;
t400 = -t437 * pkin(1) + t409;
t346 = -t434 * t399 + t430 * t400;
t337 = t421 * pkin(2) - t346;
t347 = t430 * t399 + t434 * t400;
t420 = t422 ^ 2;
t338 = -t420 * pkin(2) + t347;
t298 = t426 * t337 + t427 * t338;
t291 = -t420 * pkin(3) + t298;
t439 = t427 * t337 - t426 * t338;
t438 = t421 * pkin(3) + t439;
t247 = t429 * t291 - t433 * t438;
t248 = t433 * t291 + t429 * t438;
t447 = t429 * t247 + t433 * t248;
t206 = t433 * t247 - t429 * t248;
t467 = t427 * t206;
t180 = -t426 * t447 + t467;
t469 = t426 * t206;
t492 = t427 * t447 + t469;
t167 = t434 * t180 - t430 * t492;
t511 = t430 * t180 + t434 * t492;
t153 = t431 * t167 + t435 * t511;
t152 = t435 * t167 - t431 * t511;
t388 = t427 * t420 + t426 * t421;
t391 = t426 * t420 - t427 * t421;
t331 = t434 * t388 - t430 * t391;
t365 = qJ(3) * t388 - t427 * t425;
t487 = qJ(3) * t391 - t426 * t425;
t279 = pkin(6) * t331 + t434 * t365 - t430 * t487;
t335 = t430 * t388 + t434 * t391;
t502 = pkin(6) * t335 + t430 * t365 + t434 * t487;
t504 = t431 * t331 + t435 * t335;
t518 = pkin(5) * t504 + t431 * t279 + t435 * t502;
t293 = t435 * t331 - t431 * t335;
t516 = pkin(5) * t293 + t435 * t279 - t431 * t502;
t446 = t427 * t298 - t426 * t439;
t252 = -t426 * t298 - t427 * t439;
t453 = t434 * t252;
t210 = -t430 * t446 + t453;
t459 = t430 * t252;
t493 = t434 * t446 + t459;
t185 = t431 * t210 + t435 * t493;
t184 = t435 * t210 - t431 * t493;
t393 = t434 * t420 + t430 * t421;
t370 = pkin(6) * t393 - t434 * g(1);
t396 = t430 * t420 - t434 * t421;
t485 = t431 * t393 + t435 * t396;
t489 = pkin(6) * t396 - t430 * g(1);
t505 = pkin(5) * t485 + t431 * t370 + t435 * t489;
t339 = t435 * t393 - t431 * t396;
t503 = pkin(5) * t339 + t435 * t370 - t431 * t489;
t445 = t430 * t346 + t434 * t347;
t302 = t434 * t346 - t430 * t347;
t452 = t435 * t302;
t254 = -t431 * t445 + t452;
t458 = t431 * t302;
t491 = -t435 * t445 - t458;
t474 = pkin(1) * t425;
t473 = pkin(2) * t425;
t428 = sin(qJ(5));
t423 = t428 ^ 2;
t470 = t423 * t417;
t242 = -t418 * pkin(4) - t417 * pkin(8) + t247;
t465 = t428 * t242;
t432 = cos(qJ(5));
t405 = t432 * t417 * t428;
t397 = qJDD(5) + t405;
t464 = t428 * t397;
t398 = qJDD(5) - t405;
t463 = t428 * t398;
t462 = t428 * t418;
t457 = t432 * t242;
t456 = t432 * t398;
t411 = t432 * t418;
t243 = -t417 * pkin(4) + t418 * pkin(8) + t248;
t236 = t432 * t243 - t428 * t425;
t424 = t432 ^ 2;
t451 = t423 + t424;
t450 = qJD(5) * t419;
t449 = t428 * t450;
t448 = t432 * t450;
t235 = t428 * t243 + t432 * t425;
t200 = t428 * t235 + t432 * t236;
t442 = t429 * t405;
t441 = t433 * t405;
t406 = t431 * qJDD(1) + t435 * t437;
t440 = pkin(5) * t406 - t435 * g(1);
t199 = t432 * t235 - t428 * t236;
t361 = t435 * t409 - t431 * t410;
t360 = -t431 * t409 - t435 * t410;
t436 = qJD(5) ^ 2;
t412 = t424 * t417;
t407 = -t435 * qJDD(1) + t431 * t437;
t404 = -t412 - t436;
t403 = t412 - t436;
t402 = -t436 - t470;
t401 = t436 - t470;
t385 = t432 * t397;
t384 = -pkin(5) * t407 + t431 * g(1);
t383 = t412 - t470;
t382 = t412 + t470;
t377 = t451 * t418;
t375 = t411 - 0.2e1 * t449;
t374 = t411 - t449;
t373 = t448 + t462;
t372 = 0.2e1 * t448 + t462;
t371 = t451 * t450;
t355 = t429 * qJDD(5) + t433 * t371;
t354 = -t433 * qJDD(5) + t429 * t371;
t353 = -t428 * t402 - t456;
t352 = -t428 * t401 + t385;
t351 = t432 * t404 - t464;
t350 = t432 * t403 - t463;
t349 = t432 * t402 - t463;
t348 = t428 * t404 + t385;
t345 = t432 * t373 - t423 * t450;
t344 = -t428 * t374 - t424 * t450;
t328 = t433 * t377 - t429 * t382;
t327 = t429 * t377 + t433 * t382;
t320 = -t428 * t372 + t432 * t375;
t319 = t433 * t352 + t428 * t461;
t318 = t433 * t350 + t429 * t411;
t317 = t429 * t352 - t428 * t455;
t316 = t429 * t350 - t432 * t455;
t315 = t433 * t345 - t442;
t314 = t433 * t344 + t442;
t313 = t429 * t345 + t441;
t312 = t429 * t344 - t441;
t311 = t433 * t353 + t429 * t372;
t310 = t433 * t351 - t429 * t375;
t309 = t429 * t353 - t433 * t372;
t308 = t429 * t351 + t433 * t375;
t307 = -t426 * t354 + t427 * t355;
t306 = t427 * t354 + t426 * t355;
t305 = t433 * t320 - t429 * t383;
t304 = t429 * t320 + t433 * t383;
t299 = pkin(1) * g(1) + pkin(6) * t445;
t287 = -t426 * t327 + t427 * t328;
t286 = t427 * t327 + t426 * t328;
t275 = -t426 * t317 + t427 * t319;
t274 = -t426 * t316 + t427 * t318;
t273 = t427 * t317 + t426 * t319;
t272 = t427 * t316 + t426 * t318;
t267 = -t426 * t313 + t427 * t315;
t266 = -t426 * t312 + t427 * t314;
t265 = t427 * t313 + t426 * t315;
t264 = t427 * t312 + t426 * t314;
t263 = -t426 * t309 + t427 * t311;
t262 = -t426 * t308 + t427 * t310;
t261 = t427 * t309 + t426 * t311;
t260 = t427 * t308 + t426 * t310;
t259 = -t430 * t306 + t434 * t307;
t258 = t434 * t306 + t430 * t307;
t257 = -t426 * t304 + t427 * t305;
t256 = t427 * t304 + t426 * t305;
t249 = qJ(3) * t446 + t473;
t245 = -t430 * t286 + t434 * t287;
t244 = t434 * t286 + t430 * t287;
t233 = -t430 * t273 + t434 * t275;
t232 = -t430 * t272 + t434 * t274;
t231 = t434 * t273 + t430 * t275;
t230 = t434 * t272 + t430 * t274;
t229 = -t430 * t265 + t434 * t267;
t228 = -t430 * t264 + t434 * t266;
t227 = t434 * t265 + t430 * t267;
t226 = t434 * t264 + t430 * t266;
t225 = -pkin(8) * t349 + t457;
t224 = -pkin(8) * t348 + t465;
t223 = -pkin(4) * t349 + t236;
t222 = -pkin(4) * t348 + t235;
t221 = -t430 * t261 + t434 * t263;
t220 = -t430 * t260 + t434 * t262;
t219 = t434 * t261 + t430 * t263;
t218 = t434 * t260 + t430 * t262;
t217 = -t430 * t256 + t434 * t257;
t216 = t434 * t256 + t430 * t257;
t203 = pkin(3) * t425 + pkin(7) * t447;
t202 = -t431 * t244 + t435 * t245;
t201 = -t435 * t244 - t431 * t245;
t197 = -pkin(7) * t327 + t433 * t199;
t196 = pkin(7) * t328 + t429 * t199;
t195 = -t431 * t219 + t435 * t221;
t194 = -t431 * t218 + t435 * t220;
t193 = -t435 * t219 - t431 * t221;
t192 = -t435 * t218 - t431 * t220;
t191 = -pkin(7) * t309 - t429 * t223 + t433 * t225;
t190 = -pkin(7) * t308 - t429 * t222 + t433 * t224;
t189 = -pkin(3) * t349 + pkin(7) * t311 + t433 * t223 + t429 * t225;
t188 = -pkin(3) * t348 + pkin(7) * t310 + t433 * t222 + t429 * t224;
t187 = t433 * t200 + t429 * t242;
t186 = t429 * t200 - t433 * t242;
t183 = pkin(6) * t210 + qJ(3) * t453 - t430 * t249;
t182 = pkin(6) * t493 + qJ(3) * t459 + t434 * t249 + t474;
t177 = -qJ(3) * t286 - t426 * t196 + t427 * t197;
t176 = qJ(3) * t287 + t427 * t196 + t426 * t197;
t175 = -qJ(3) * t261 - t426 * t189 + t427 * t191;
t174 = -qJ(3) * t260 - t426 * t188 + t427 * t190;
t173 = -t426 * t186 + t427 * t187;
t172 = t427 * t186 + t426 * t187;
t171 = -pkin(2) * t349 + qJ(3) * t263 + t427 * t189 + t426 * t191;
t170 = -pkin(2) * t348 + qJ(3) * t262 + t427 * t188 + t426 * t190;
t169 = -pkin(7) * t186 - (pkin(4) * t429 - pkin(8) * t433) * t199;
t164 = pkin(7) * t467 + qJ(3) * t180 - t426 * t203;
t163 = pkin(7) * t469 + qJ(3) * t492 + t427 * t203 + t473;
t162 = pkin(7) * t187 - (-pkin(4) * t433 - pkin(8) * t429 - pkin(3)) * t199;
t161 = -pkin(6) * t244 - t430 * t176 + t434 * t177;
t160 = pkin(6) * t245 + t434 * t176 + t430 * t177;
t159 = -t430 * t172 + t434 * t173;
t158 = t434 * t172 + t430 * t173;
t157 = -pkin(6) * t219 - t430 * t171 + t434 * t175;
t156 = -pkin(6) * t218 - t430 * t170 + t434 * t174;
t155 = -pkin(1) * t349 + pkin(6) * t221 + t434 * t171 + t430 * t175;
t154 = -pkin(1) * t348 + pkin(6) * t220 + t434 * t170 + t430 * t174;
t151 = pkin(6) * t167 - t430 * t163 + t434 * t164;
t150 = pkin(6) * t511 + t434 * t163 + t430 * t164 + t474;
t149 = -qJ(3) * t172 - t426 * t162 + t427 * t169;
t148 = -t431 * t158 + t435 * t159;
t147 = -t435 * t158 - t431 * t159;
t146 = pkin(2) * t199 + qJ(3) * t173 + t427 * t162 + t426 * t169;
t145 = -pkin(6) * t158 - t430 * t146 + t434 * t149;
t144 = pkin(1) * t199 + pkin(6) * t159 + t434 * t146 + t430 * t149;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, 0, 0, 0, 0, 0, 0, t348, t349, 0, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t407, t406, 0, t360, 0, 0, 0, 0, 0, 0, t485, t339, 0, t254, 0, 0, 0, 0, 0, 0, t504, t293, 0, t184, 0, 0, 0, 0, 0, 0, t517, t238, 0, t152, 0, 0, 0, 0, 0, 0, t192, t193, t201, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t406, t407, 0, t361, 0, 0, 0, 0, 0, 0, -t339, t485, 0, -t491, 0, 0, 0, 0, 0, 0, -t293, t504, 0, t185, 0, 0, 0, 0, 0, 0, -t238, t517, 0, t153, 0, 0, 0, 0, 0, 0, t194, t195, t202, t148; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), t410, -t409, 0, 0, 0, 0, 0, 0, 0, t421, -pkin(1) * t396 - t346, -pkin(1) * t393 - t347, 0, -pkin(1) * t302, 0, 0, 0, 0, 0, t421, -pkin(1) * t335 - pkin(2) * t391 + t439, -pkin(1) * t331 - pkin(2) * t388 - t298, 0, -pkin(1) * t210 - pkin(2) * t252, 0, 0, 0, 0, 0, t418, -pkin(1) * t284 - pkin(2) * t325 - pkin(3) * t381 - t247, -pkin(1) * t280 - pkin(2) * t321 - pkin(3) * t378 - t248, 0, -pkin(1) * t167 - pkin(2) * t180 - pkin(3) * t206, (t373 + t448) * t428, t432 * t372 + t428 * t375, t432 * t401 + t464, (t374 - t449) * t432, t428 * t403 + t456, 0, pkin(1) * t218 + pkin(2) * t260 + pkin(3) * t308 + pkin(4) * t375 + pkin(8) * t351 - t457, pkin(1) * t219 + pkin(2) * t261 + pkin(3) * t309 - pkin(4) * t372 + pkin(8) * t353 + t465, pkin(1) * t244 + pkin(2) * t286 + pkin(3) * t327 + pkin(4) * t382 + pkin(8) * t377 + t200, pkin(1) * t158 + pkin(2) * t172 + pkin(3) * t186 - pkin(4) * t242 + pkin(8) * t200; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t406, 0, t407, 0, t440, t384, -t361, -pkin(5) * t361, 0, 0, -t339, 0, t485, 0, t503, -t505, t491, pkin(5) * t491 - pkin(6) * t458 - t435 * t299, 0, 0, -t293, 0, t504, 0, t516, -t518, -t185, -pkin(5) * t185 - t435 * t182 - t431 * t183, 0, 0, -t238, 0, t517, 0, t523, -t524, -t153, -pkin(5) * t153 - t435 * t150 - t431 * t151, -t435 * t227 - t431 * t229, -t435 * t216 - t431 * t217, -t435 * t231 - t431 * t233, -t435 * t226 - t431 * t228, -t435 * t230 - t431 * t232, -t435 * t258 - t431 * t259, -pkin(5) * t194 - t435 * t154 - t431 * t156, -pkin(5) * t195 - t435 * t155 - t431 * t157, -pkin(5) * t202 - t435 * t160 - t431 * t161, -pkin(5) * t148 - t435 * t144 - t431 * t145; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t407, 0, -t406, 0, -t384, t440, t360, pkin(5) * t360, 0, 0, -t485, 0, -t339, 0, t505, t503, t254, pkin(5) * t254 + pkin(6) * t452 - t431 * t299, 0, 0, -t504, 0, -t293, 0, t518, t516, t184, pkin(5) * t184 - t431 * t182 + t435 * t183, 0, 0, -t517, 0, -t238, 0, t524, t523, t152, pkin(5) * t152 - t431 * t150 + t435 * t151, -t431 * t227 + t435 * t229, -t431 * t216 + t435 * t217, -t431 * t231 + t435 * t233, -t431 * t226 + t435 * t228, -t431 * t230 + t435 * t232, -t431 * t258 + t435 * t259, pkin(5) * t192 - t431 * t154 + t435 * t156, pkin(5) * t193 - t431 * t155 + t435 * t157, pkin(5) * t201 - t431 * t160 + t435 * t161, pkin(5) * t147 - t431 * t144 + t435 * t145;];
tauB_reg = t1;
