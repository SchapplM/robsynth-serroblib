% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:41:02
% EndTime: 2019-12-05 18:41:12
% DurationCPUTime: 8.84s
% Computational Cost: add. (37767->383), mult. (48058->558), div. (0->0), fcn. (28562->10), ass. (0->254)
t423 = qJD(1) + qJD(2);
t420 = qJD(3) + t423;
t418 = t420 ^ 2;
t422 = qJDD(1) + qJDD(2);
t419 = qJDD(3) + t422;
t427 = sin(pkin(9));
t428 = cos(pkin(9));
t378 = t428 * t418 + t427 * t419;
t381 = t427 * t418 - t428 * t419;
t430 = sin(qJ(3));
t434 = cos(qJ(3));
t322 = t434 * t378 - t430 * t381;
t426 = g(1) - qJDD(4);
t360 = qJ(4) * t378 - t428 * t426;
t485 = qJ(4) * t381 - t427 * t426;
t272 = pkin(7) * t322 + t434 * t360 - t430 * t485;
t326 = t430 * t378 + t434 * t381;
t431 = sin(qJ(2));
t435 = cos(qJ(2));
t281 = t435 * t322 - t431 * t326;
t498 = pkin(7) * t326 + t430 * t360 + t434 * t485;
t216 = pkin(6) * t281 + t435 * t272 - t431 * t498;
t432 = sin(qJ(1));
t436 = cos(qJ(1));
t285 = t431 * t322 + t435 * t326;
t512 = pkin(6) * t285 + t431 * t272 + t435 * t498;
t514 = t432 * t281 + t436 * t285;
t521 = pkin(5) * t514 + t432 * t216 + t436 * t512;
t239 = t436 * t281 - t432 * t285;
t520 = pkin(5) * t239 + t436 * t216 - t432 * t512;
t411 = t436 * g(2) + t432 * g(3);
t400 = qJDD(1) * pkin(1) + t411;
t410 = t432 * g(2) - t436 * g(3);
t439 = qJD(1) ^ 2;
t401 = -t439 * pkin(1) + t410;
t347 = -t435 * t400 + t431 * t401;
t338 = t422 * pkin(2) - t347;
t348 = t431 * t400 + t435 * t401;
t421 = t423 ^ 2;
t339 = -t421 * pkin(2) + t348;
t299 = t430 * t338 + t434 * t339;
t296 = -t418 * pkin(3) + t299;
t441 = t434 * t338 - t430 * t339;
t440 = t419 * pkin(3) + t441;
t248 = t427 * t296 - t428 * t440;
t249 = t428 * t296 + t427 * t440;
t449 = t427 * t248 + t428 * t249;
t207 = t428 * t248 - t427 * t249;
t456 = t434 * t207;
t181 = -t430 * t449 + t456;
t461 = t430 * t207;
t489 = t434 * t449 + t461;
t168 = t435 * t181 - t431 * t489;
t508 = t431 * t181 + t435 * t489;
t154 = t432 * t168 + t436 * t508;
t153 = t436 * t168 - t432 * t508;
t383 = t434 * t418 + t430 * t419;
t386 = t430 * t418 - t434 * t419;
t330 = t435 * t383 - t431 * t386;
t366 = pkin(7) * t383 - t434 * g(1);
t486 = pkin(7) * t386 - t430 * g(1);
t280 = pkin(6) * t330 + t435 * t366 - t431 * t486;
t334 = t431 * t383 + t435 * t386;
t499 = pkin(6) * t334 + t431 * t366 + t435 * t486;
t501 = t432 * t330 + t436 * t334;
t515 = pkin(5) * t501 + t432 * t280 + t436 * t499;
t289 = t436 * t330 - t432 * t334;
t513 = pkin(5) * t289 + t436 * t280 - t432 * t499;
t448 = t434 * t299 - t430 * t441;
t253 = -t430 * t299 - t434 * t441;
t455 = t435 * t253;
t211 = -t431 * t448 + t455;
t460 = t431 * t253;
t490 = t435 * t448 + t460;
t185 = t432 * t211 + t436 * t490;
t184 = t436 * t211 - t432 * t490;
t394 = t435 * t421 + t431 * t422;
t371 = pkin(6) * t394 - t435 * g(1);
t397 = t431 * t421 - t435 * t422;
t483 = t432 * t394 + t436 * t397;
t487 = pkin(6) * t397 - t431 * g(1);
t502 = pkin(5) * t483 + t432 * t371 + t436 * t487;
t340 = t436 * t394 - t432 * t397;
t500 = pkin(5) * t340 + t436 * t371 - t432 * t487;
t447 = t431 * t347 + t435 * t348;
t305 = t435 * t347 - t431 * t348;
t454 = t436 * t305;
t256 = -t432 * t447 + t454;
t459 = t432 * t305;
t258 = t436 * t447 + t459;
t429 = sin(qJ(5));
t424 = t429 ^ 2;
t468 = t424 * t418;
t243 = -t419 * pkin(4) - t418 * pkin(8) + t248;
t465 = t429 * t243;
t433 = cos(qJ(5));
t406 = t433 * t418 * t429;
t398 = qJDD(5) + t406;
t464 = t429 * t398;
t399 = qJDD(5) - t406;
t463 = t429 * t399;
t462 = t429 * t419;
t458 = t433 * t243;
t457 = t433 * t399;
t412 = t433 * t419;
t244 = -t418 * pkin(4) + t419 * pkin(8) + t249;
t237 = t433 * t244 - t429 * t426;
t425 = t433 ^ 2;
t453 = t424 + t425;
t452 = qJD(5) * t420;
t451 = t429 * t452;
t450 = t433 * t452;
t236 = t429 * t244 + t433 * t426;
t201 = t429 * t236 + t433 * t237;
t444 = t427 * t406;
t443 = t428 * t406;
t407 = t432 * qJDD(1) + t436 * t439;
t442 = pkin(5) * t407 - t436 * g(1);
t200 = t433 * t236 - t429 * t237;
t362 = t436 * t410 - t432 * t411;
t361 = -t432 * t410 - t436 * t411;
t438 = qJD(5) ^ 2;
t437 = pkin(1) * g(1);
t413 = t425 * t418;
t408 = -t436 * qJDD(1) + t432 * t439;
t405 = -t413 - t438;
t404 = t413 - t438;
t403 = -t438 - t468;
t402 = t438 - t468;
t390 = t433 * t398;
t389 = -pkin(5) * t408 + t432 * g(1);
t388 = t413 - t468;
t387 = t413 + t468;
t382 = t453 * t419;
t377 = t412 - 0.2e1 * t451;
t376 = t412 - t451;
t375 = t450 + t462;
t374 = 0.2e1 * t450 + t462;
t373 = t453 * t452;
t356 = t427 * qJDD(5) + t428 * t373;
t355 = -t428 * qJDD(5) + t427 * t373;
t354 = -t429 * t403 - t457;
t353 = -t429 * t402 + t390;
t352 = t433 * t405 - t464;
t351 = t433 * t404 - t463;
t350 = t433 * t403 - t463;
t349 = t429 * t405 + t390;
t346 = t433 * t375 - t424 * t452;
t345 = -t429 * t376 - t425 * t452;
t329 = t428 * t382 - t427 * t387;
t328 = t427 * t382 + t428 * t387;
t321 = -t429 * t374 + t433 * t377;
t320 = t428 * t353 + t427 * t462;
t319 = t428 * t351 + t427 * t412;
t318 = t427 * t353 - t428 * t462;
t317 = t427 * t351 - t428 * t412;
t316 = t428 * t346 - t444;
t315 = t428 * t345 + t444;
t314 = t427 * t346 + t443;
t313 = t427 * t345 - t443;
t312 = t428 * t354 + t427 * t374;
t311 = t428 * t352 - t427 * t377;
t310 = t427 * t354 - t428 * t374;
t309 = t427 * t352 + t428 * t377;
t308 = -t430 * t355 + t434 * t356;
t307 = t434 * t355 + t430 * t356;
t304 = t428 * t321 - t427 * t388;
t301 = t427 * t321 + t428 * t388;
t300 = pkin(6) * t447 + t437;
t288 = -t430 * t328 + t434 * t329;
t287 = t434 * t328 + t430 * t329;
t276 = -t430 * t318 + t434 * t320;
t275 = -t430 * t317 + t434 * t319;
t274 = t434 * t318 + t430 * t320;
t273 = t434 * t317 + t430 * t319;
t268 = -t430 * t314 + t434 * t316;
t267 = -t430 * t313 + t434 * t315;
t266 = t434 * t314 + t430 * t316;
t265 = t434 * t313 + t430 * t315;
t264 = -t430 * t310 + t434 * t312;
t263 = -t430 * t309 + t434 * t311;
t262 = t434 * t310 + t430 * t312;
t261 = t434 * t309 + t430 * t311;
t260 = -t431 * t307 + t435 * t308;
t259 = t435 * t307 + t431 * t308;
t257 = -t430 * t301 + t434 * t304;
t255 = t434 * t301 + t430 * t304;
t250 = pkin(2) * g(1) + pkin(7) * t448;
t246 = -t431 * t287 + t435 * t288;
t245 = t435 * t287 + t431 * t288;
t234 = -t431 * t274 + t435 * t276;
t233 = -t431 * t273 + t435 * t275;
t232 = t435 * t274 + t431 * t276;
t231 = t435 * t273 + t431 * t275;
t230 = -t431 * t266 + t435 * t268;
t229 = -t431 * t265 + t435 * t267;
t228 = t435 * t266 + t431 * t268;
t227 = t435 * t265 + t431 * t267;
t226 = -pkin(8) * t350 + t458;
t225 = -pkin(8) * t349 + t465;
t224 = -pkin(4) * t350 + t237;
t223 = -pkin(4) * t349 + t236;
t222 = -t431 * t262 + t435 * t264;
t221 = -t431 * t261 + t435 * t263;
t220 = t435 * t262 + t431 * t264;
t219 = t435 * t261 + t431 * t263;
t218 = -t431 * t255 + t435 * t257;
t217 = t435 * t255 + t431 * t257;
t204 = pkin(3) * t426 + qJ(4) * t449;
t203 = -t432 * t245 + t436 * t246;
t202 = -t436 * t245 - t432 * t246;
t198 = -qJ(4) * t328 + t428 * t200;
t197 = qJ(4) * t329 + t427 * t200;
t196 = -t432 * t220 + t436 * t222;
t195 = -t432 * t219 + t436 * t221;
t194 = -t436 * t220 - t432 * t222;
t193 = -t436 * t219 - t432 * t221;
t192 = -qJ(4) * t310 - t427 * t224 + t428 * t226;
t191 = -qJ(4) * t309 - t427 * t223 + t428 * t225;
t190 = -pkin(3) * t350 + qJ(4) * t312 + t428 * t224 + t427 * t226;
t189 = -pkin(3) * t349 + qJ(4) * t311 + t428 * t223 + t427 * t225;
t188 = t428 * t201 + t427 * t243;
t187 = t427 * t201 - t428 * t243;
t186 = pkin(6) * t211 + pkin(7) * t455 - t431 * t250;
t183 = pkin(6) * t490 + pkin(7) * t460 + t435 * t250 + t437;
t178 = -pkin(7) * t287 - t430 * t197 + t434 * t198;
t177 = pkin(7) * t288 + t434 * t197 + t430 * t198;
t176 = -pkin(7) * t262 - t430 * t190 + t434 * t192;
t175 = -pkin(7) * t261 - t430 * t189 + t434 * t191;
t174 = -t430 * t187 + t434 * t188;
t173 = t434 * t187 + t430 * t188;
t172 = -pkin(2) * t350 + pkin(7) * t264 + t434 * t190 + t430 * t192;
t171 = -pkin(2) * t349 + pkin(7) * t263 + t434 * t189 + t430 * t191;
t170 = -qJ(4) * t187 - (pkin(4) * t427 - pkin(8) * t428) * t200;
t165 = pkin(7) * t181 + qJ(4) * t456 - t430 * t204;
t164 = pkin(2) * t426 + pkin(7) * t489 + qJ(4) * t461 + t434 * t204;
t163 = qJ(4) * t188 - (-pkin(4) * t428 - pkin(8) * t427 - pkin(3)) * t200;
t162 = -pkin(6) * t245 - t431 * t177 + t435 * t178;
t161 = pkin(6) * t246 + t435 * t177 + t431 * t178;
t160 = -t431 * t173 + t435 * t174;
t159 = t435 * t173 + t431 * t174;
t158 = -pkin(6) * t220 - t431 * t172 + t435 * t176;
t157 = -pkin(6) * t219 - t431 * t171 + t435 * t175;
t156 = -pkin(1) * t350 + pkin(6) * t222 + t435 * t172 + t431 * t176;
t155 = -pkin(1) * t349 + pkin(6) * t221 + t435 * t171 + t431 * t175;
t152 = pkin(6) * t168 - t431 * t164 + t435 * t165;
t151 = pkin(1) * t426 + pkin(6) * t508 + t435 * t164 + t431 * t165;
t150 = -pkin(7) * t173 - t430 * t163 + t434 * t170;
t149 = -t432 * t159 + t436 * t160;
t148 = -t436 * t159 - t432 * t160;
t147 = pkin(2) * t200 + pkin(7) * t174 + t434 * t163 + t430 * t170;
t146 = -pkin(6) * t159 - t431 * t147 + t435 * t150;
t145 = pkin(1) * t200 + pkin(6) * t160 + t435 * t147 + t431 * t150;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, 0, 0, 0, 0, 0, 0, t349, t350, 0, -t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t408, t407, 0, t361, 0, 0, 0, 0, 0, 0, t483, t340, 0, t256, 0, 0, 0, 0, 0, 0, t501, t289, 0, t184, 0, 0, 0, 0, 0, 0, t514, t239, 0, t153, 0, 0, 0, 0, 0, 0, t193, t194, t202, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t407, t408, 0, t362, 0, 0, 0, 0, 0, 0, -t340, t483, 0, t258, 0, 0, 0, 0, 0, 0, -t289, t501, 0, t185, 0, 0, 0, 0, 0, 0, -t239, t514, 0, t154, 0, 0, 0, 0, 0, 0, t195, t196, t203, t149; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), t411, -t410, 0, 0, 0, 0, 0, 0, 0, t422, -pkin(1) * t397 - t347, -pkin(1) * t394 - t348, 0, -pkin(1) * t305, 0, 0, 0, 0, 0, t419, -pkin(1) * t334 - pkin(2) * t386 + t441, -pkin(1) * t330 - pkin(2) * t383 - t299, 0, -pkin(1) * t211 - pkin(2) * t253, 0, 0, 0, 0, 0, t419, -pkin(1) * t285 - pkin(2) * t326 - pkin(3) * t381 - t248, -pkin(1) * t281 - pkin(2) * t322 - pkin(3) * t378 - t249, 0, -pkin(1) * t168 - pkin(2) * t181 - pkin(3) * t207, (t375 + t450) * t429, t433 * t374 + t429 * t377, t433 * t402 + t464, (t376 - t451) * t433, t429 * t404 + t457, 0, pkin(1) * t219 + pkin(2) * t261 + pkin(3) * t309 + pkin(4) * t377 + pkin(8) * t352 - t458, pkin(1) * t220 + pkin(2) * t262 + pkin(3) * t310 - pkin(4) * t374 + pkin(8) * t354 + t465, pkin(1) * t245 + pkin(2) * t287 + pkin(3) * t328 + pkin(4) * t387 + pkin(8) * t382 + t201, pkin(1) * t159 + pkin(2) * t173 + pkin(3) * t187 - pkin(4) * t243 + pkin(8) * t201; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t407, 0, t408, 0, t442, t389, -t362, -pkin(5) * t362, 0, 0, -t340, 0, t483, 0, t500, -t502, -t258, -pkin(5) * t258 - pkin(6) * t459 - t436 * t300, 0, 0, -t289, 0, t501, 0, t513, -t515, -t185, -pkin(5) * t185 - t436 * t183 - t432 * t186, 0, 0, -t239, 0, t514, 0, t520, -t521, -t154, -pkin(5) * t154 - t436 * t151 - t432 * t152, -t436 * t228 - t432 * t230, -t436 * t217 - t432 * t218, -t436 * t232 - t432 * t234, -t436 * t227 - t432 * t229, -t436 * t231 - t432 * t233, -t436 * t259 - t432 * t260, -pkin(5) * t195 - t436 * t155 - t432 * t157, -pkin(5) * t196 - t436 * t156 - t432 * t158, -pkin(5) * t203 - t436 * t161 - t432 * t162, -pkin(5) * t149 - t436 * t145 - t432 * t146; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t408, 0, -t407, 0, -t389, t442, t361, pkin(5) * t361, 0, 0, -t483, 0, -t340, 0, t502, t500, t256, pkin(5) * t256 + pkin(6) * t454 - t432 * t300, 0, 0, -t501, 0, -t289, 0, t515, t513, t184, pkin(5) * t184 - t432 * t183 + t436 * t186, 0, 0, -t514, 0, -t239, 0, t521, t520, t153, pkin(5) * t153 - t432 * t151 + t436 * t152, -t432 * t228 + t436 * t230, -t432 * t217 + t436 * t218, -t432 * t232 + t436 * t234, -t432 * t227 + t436 * t229, -t432 * t231 + t436 * t233, -t432 * t259 + t436 * t260, pkin(5) * t193 - t432 * t155 + t436 * t157, pkin(5) * t194 - t432 * t156 + t436 * t158, pkin(5) * t202 - t432 * t161 + t436 * t162, pkin(5) * t148 - t432 * t145 + t436 * t146;];
tauB_reg = t1;
