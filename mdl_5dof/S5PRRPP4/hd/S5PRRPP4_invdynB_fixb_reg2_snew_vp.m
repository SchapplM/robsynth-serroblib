% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRRPP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:14
% EndTime: 2019-12-31 17:41:21
% DurationCPUTime: 4.63s
% Computational Cost: add. (6308->379), mult. (13676->464), div. (0->0), fcn. (7429->6), ass. (0->253)
t396 = qJD(3) ^ 2;
t392 = sin(qJ(3));
t382 = t392 ^ 2;
t397 = qJD(2) ^ 2;
t467 = t382 * t397;
t362 = t396 + t467;
t394 = cos(qJ(3));
t430 = t392 * t397 * t394;
t360 = qJDD(3) - t430;
t445 = t394 * t360;
t316 = -t392 * t362 + t445;
t437 = qJD(2) * qJD(3);
t423 = t394 * t437;
t435 = t392 * qJDD(2);
t348 = 0.2e1 * t423 + t435;
t393 = sin(qJ(2));
t395 = cos(qJ(2));
t262 = t393 * t316 + t395 * t348;
t266 = t395 * t316 - t393 * t348;
t388 = sin(pkin(7));
t389 = cos(pkin(7));
t212 = t389 * t262 + t388 * t266;
t473 = qJ(1) * t212;
t216 = t388 * t262 - t389 * t266;
t210 = qJ(1) * t216;
t479 = pkin(5) * t262;
t428 = pkin(1) * t262 + pkin(2) * t348 + pkin(6) * t316;
t454 = t392 * t360;
t310 = t394 * t362 + t454;
t419 = -pkin(1) * t310 + pkin(5) * t266;
t359 = qJDD(3) + t430;
t342 = t394 * t359;
t363 = -t396 + t467;
t315 = t392 * t363 + t342;
t432 = t395 * qJDD(2);
t272 = t393 * t315 - t392 * t432;
t434 = t393 * qJDD(2);
t276 = t395 * t315 + t392 * t434;
t226 = t389 * t272 + t388 * t276;
t503 = t388 * t272 - t389 * t276;
t383 = t394 ^ 2;
t440 = t382 + t383;
t352 = t440 * qJDD(2);
t355 = t440 * t397;
t294 = t393 * t352 + t395 * t355;
t296 = t395 * t352 - t393 * t355;
t239 = t389 * t294 + t388 * t296;
t472 = qJ(1) * t239;
t242 = t388 * t294 - t389 * t296;
t238 = qJ(1) * t242;
t306 = t392 * t348;
t374 = t392 * t437;
t433 = t394 * qJDD(2);
t350 = -0.2e1 * t374 + t433;
t446 = t394 * t350;
t286 = -t446 + t306;
t356 = (t382 - t383) * t397;
t252 = t393 * t286 + t395 * t356;
t254 = t395 * t286 - t393 * t356;
t207 = t389 * t252 + t388 * t254;
t502 = t388 * t252 - t389 * t254;
t466 = t383 * t397;
t364 = -t396 + t466;
t314 = -t394 * t364 + t454;
t271 = t393 * t314 + t394 * t432;
t275 = t395 * t314 - t393 * t433;
t501 = t389 * t271 + t388 * t275;
t227 = t388 * t271 - t389 * t275;
t357 = t389 * g(1) + t388 * g(2);
t418 = t388 * g(1) - t389 * g(2);
t415 = t393 * t357 + t395 * t418;
t442 = t395 * t357 - t393 * t418;
t417 = -t393 * t415 - t395 * t442;
t247 = t393 * t442 - t395 * t415;
t460 = t389 * t247;
t500 = -t388 * t417 + t460;
t465 = t388 * t247;
t198 = t389 * t417 + t465;
t353 = t395 * t397 + t434;
t354 = -t393 * t397 + t432;
t291 = -t388 * t353 + t389 * t354;
t384 = g(3) - qJDD(1);
t326 = pkin(5) * t353 - t395 * t384;
t408 = -pkin(5) * t354 - t393 * t384;
t498 = -qJ(1) * t291 + t388 * t326 + t389 * t408;
t496 = pkin(2) * t310;
t478 = pkin(5) * t294;
t288 = pkin(5) * t296;
t476 = pkin(6) * t310;
t427 = pkin(1) * t294 + pkin(2) * t355 + pkin(6) * t352;
t490 = t389 * t353 + t388 * t354;
t494 = qJ(1) * t490 + t389 * t326 - t388 * t408;
t307 = t392 * t364 + t445;
t365 = -t396 - t466;
t308 = t392 * t365 + t342;
t482 = pkin(2) * t308;
t488 = -qJ(4) * t365 - t482;
t439 = qJD(2) * t392;
t358 = -qJD(3) * pkin(4) - qJ(5) * t439;
t487 = t358 * t439 + qJDD(5);
t436 = qJD(4) * qJD(3);
t379 = -0.2e1 * t436;
t486 = -qJ(4) * t360 + t379 - t496;
t349 = -t374 + t433;
t485 = t349 * pkin(4) + t487;
t299 = -t389 * t357 - t388 * t418;
t298 = -t388 * t357 + t389 * t418;
t373 = t394 * t384;
t411 = qJDD(3) * pkin(3) + t396 * qJ(4) - qJDD(4) - t373;
t283 = -t397 * pkin(2) + qJDD(2) * pkin(6) - t442;
t410 = -pkin(3) * t394 - qJ(4) * t392;
t346 = t410 * qJD(2);
t414 = qJD(2) * t346 + t283;
t399 = t414 * t392 - t411;
t483 = pkin(3) + pkin(4);
t455 = t392 * t359;
t313 = t394 * t365 - t455;
t261 = t393 * t313 + t395 * t350;
t480 = pkin(5) * t261;
t477 = pkin(6) * t308;
t265 = t395 * t313 - t393 * t350;
t211 = t389 * t261 + t388 * t265;
t474 = qJ(1) * t211;
t471 = qJ(4) * t355;
t468 = qJ(4) * t394;
t461 = t388 * t384;
t459 = t389 * t384;
t405 = qJDD(2) * pkin(2) + t397 * pkin(6) + t415;
t458 = t392 * t405;
t457 = t392 * t283;
t456 = t392 * t350;
t448 = t394 * t405;
t447 = t394 * t348;
t441 = -t394 * t283 + t392 * t384;
t438 = qJD(2) * t394;
t431 = 0.2e1 * t439;
t429 = pkin(1) * t261 + pkin(2) * t350 + pkin(6) * t313;
t426 = qJD(5) * t438;
t420 = -pkin(1) * t308 + pkin(5) * t265;
t260 = t373 + t457;
t219 = t392 * t260 - t394 * t441;
t413 = t393 * t430;
t412 = t395 * t430;
t409 = pkin(3) * t392 - t468;
t407 = t396 * pkin(3) - qJDD(3) * qJ(4) - t346 * t438 + t441;
t218 = t394 * t260 + t392 * t441;
t406 = t447 + t456;
t309 = -t394 * t363 + t455;
t378 = 0.2e1 * t436;
t233 = t378 - t407;
t404 = t423 + t435;
t403 = pkin(4) * t359 + t404 * qJ(5) + t411;
t402 = t349 * pkin(3) + t405 + (t404 + t423) * qJ(4);
t401 = pkin(4) * t466 - qJD(3) * t358 + t407;
t400 = qJD(4) * t431 + t402;
t398 = t349 * qJ(5) + t401;
t224 = (pkin(3) * qJD(3) - 0.2e1 * qJD(4)) * t439 - t402;
t221 = (t350 - t374) * pkin(3) + t400;
t220 = -pkin(3) * t374 + qJ(4) * t348 + t400;
t206 = t457 + (qJ(5) * qJD(3) * t394 + (-0.2e1 * qJD(5) + t346) * t392) * qJD(2) - t403;
t370 = 0.2e1 * t426;
t345 = t409 * qJDD(2);
t341 = t440 * t437;
t327 = (-t483 * t392 + t468) * qJDD(2);
t322 = t393 * qJDD(3) + t395 * t341;
t321 = -t382 * t437 + t394 * t404;
t320 = -t395 * qJDD(3) + t393 * t341;
t319 = -t392 * t349 - t383 * t437;
t305 = (t349 - t374) * t394;
t304 = qJ(4) * t350 + qJ(5) * t359;
t281 = t395 * t321 - t413;
t280 = t395 * t319 + t413;
t279 = t393 * t321 + t412;
t278 = t393 * t319 - t412;
t268 = -qJ(5) * t360 + t483 * t348;
t250 = -t388 * t320 + t389 * t322;
t249 = t389 * t320 + t388 * t322;
t244 = -t448 + t476;
t243 = -t458 - t477;
t237 = pkin(1) * t384 + pkin(5) * t417;
t235 = -t441 + t496;
t234 = t260 - t482;
t232 = -t388 * t279 + t389 * t281;
t231 = -t388 * t278 + t389 * t280;
t230 = t389 * t279 + t388 * t281;
t229 = t389 * t278 + t388 * t280;
t223 = t399 + t471;
t222 = pkin(3) * t355 + t233;
t214 = -t388 * t261 + t389 * t265;
t209 = qJ(1) * t214;
t205 = t378 - t398 - 0.2e1 * t426;
t204 = -pkin(3) * t359 + t399 + t488;
t203 = qJ(5) * t466 + t224 - t485;
t202 = -pkin(3) * t362 + t407 + t486;
t201 = -qJ(5) * t423 - t471 + qJD(5) * t431 + (qJ(5) * qJDD(2) - t414) * t392 + t403;
t200 = t395 * t218 - t478;
t199 = t393 * t218 + t288;
t196 = t395 * t219 - t393 * t405;
t195 = t393 * t219 + t395 * t405;
t194 = t370 + t379 - t483 * t355 + (t349 + t433) * qJ(5) + t401;
t193 = t220 + (t362 - t466) * qJ(5) + t485;
t192 = -pkin(3) * t306 + t394 * t220 - t476;
t191 = qJ(4) * t446 - t392 * t221 - t477;
t190 = t394 * t233 + t392 * t399;
t189 = t392 * t233 - t394 * t399;
t188 = (-t365 - t466) * qJ(5) + (t349 + t350) * pkin(4) + t221 + t487;
t187 = -t483 * t362 + t370 + t398 + t486;
t186 = -t483 * t359 + t206 + t488;
t185 = -t392 * t222 + t394 * t223;
t184 = -t393 * t235 + t395 * t244 + t479;
t183 = -t393 * t234 + t395 * t243 - t480;
t182 = -t392 * t188 + t394 * t304 - t477;
t181 = t394 * t193 - t392 * t268 - t476;
t180 = t395 * t235 + t393 * t244 - t419;
t179 = t395 * t234 + t393 * t243 + t420;
t178 = t395 * t185 - t393 * t345 - t478;
t177 = t393 * t185 + t395 * t345 + t288;
t176 = t394 * t205 + t392 * t206;
t175 = t392 * t205 - t394 * t206;
t174 = -qJ(4) * t203 - qJ(5) * t206;
t173 = t395 * t190 + t393 * t224;
t172 = t393 * t190 - t395 * t224;
t171 = -t392 * t194 + t394 * t201;
t170 = -t388 * t195 + t389 * t196;
t169 = t389 * t195 + t388 * t196;
t168 = -pkin(2) * t189 + pkin(3) * t399 - qJ(4) * t233;
t167 = t395 * t191 - t393 * t204 - t480;
t166 = t395 * t192 - t393 * t202 - t479;
t165 = -pkin(6) * t189 + t409 * t224;
t164 = t395 * t171 - t393 * t327 + t478;
t163 = t393 * t171 + t395 * t327 - t288;
t162 = t393 * t191 + t395 * t204 + t420;
t161 = -pkin(5) * t195 - (pkin(2) * t393 - pkin(6) * t395) * t218;
t160 = t393 * t192 + t395 * t202 + t419;
t159 = -qJ(5) * t205 - t483 * t203;
t158 = t395 * t176 + t393 * t203;
t157 = t393 * t176 - t395 * t203;
t156 = t395 * t181 - t393 * t187 - t479;
t155 = t395 * t182 - t393 * t186 - t480;
t154 = pkin(5) * t196 - (-pkin(2) * t395 - pkin(6) * t393 - pkin(1)) * t218;
t153 = t393 * t181 + t395 * t187 + t419;
t152 = t393 * t182 + t395 * t186 + t420;
t151 = -t388 * t172 + t389 * t173;
t150 = t389 * t172 + t388 * t173;
t149 = -pkin(2) * t175 - qJ(4) * t205 + t483 * t206;
t148 = -t388 * t157 + t389 * t158;
t147 = t389 * t157 + t388 * t158;
t146 = -pkin(6) * t175 - t392 * t159 + t394 * t174;
t145 = -pkin(5) * t172 + t395 * t165 - t393 * t168;
t144 = -pkin(1) * t189 + pkin(5) * t173 + t393 * t165 + t395 * t168;
t143 = -pkin(5) * t157 + t395 * t146 - t393 * t149;
t142 = -pkin(1) * t175 + pkin(5) * t158 + t393 * t146 + t395 * t149;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, 0, 0, 0, 0, 0, 0, -t490, -t291, 0, t198, 0, 0, 0, 0, 0, 0, t214, t216, -t242, t170, 0, 0, 0, 0, 0, 0, t214, -t242, -t216, t151, 0, 0, 0, 0, 0, 0, t214, -t216, t242, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t298, 0, 0, 0, 0, 0, 0, t291, -t490, 0, -t500, 0, 0, 0, 0, 0, 0, t211, -t212, t239, t169, 0, 0, 0, 0, 0, 0, t211, t239, t212, t150, 0, 0, 0, 0, 0, 0, t211, t212, -t239, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t384, 0, 0, 0, 0, 0, 0, t308, -t310, 0, -t218, 0, 0, 0, 0, 0, 0, t308, 0, t310, t189, 0, 0, 0, 0, 0, 0, t308, t310, 0, t175; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t461, -t459, -t298, -qJ(1) * t298, 0, 0, t291, 0, -t490, 0, t498, t494, t500, pkin(5) * t460 + qJ(1) * t500 - t388 * t237, t232, t502, -t503, t231, t227, t250, -t388 * t179 + t389 * t183 - t474, -t388 * t180 + t389 * t184 + t473, -t388 * t199 + t389 * t200 - t472, -qJ(1) * t169 - t388 * t154 + t389 * t161, t232, -t503, -t502, t250, -t227, t231, -t388 * t162 + t389 * t167 - t474, -t388 * t177 + t389 * t178 - t472, -t388 * t160 + t389 * t166 - t473, -qJ(1) * t150 - t388 * t144 + t389 * t145, t232, -t502, t503, t231, t227, t250, -t388 * t152 + t389 * t155 - t474, -t388 * t153 + t389 * t156 - t473, -t388 * t163 + t389 * t164 + t472, -qJ(1) * t147 - t388 * t142 + t389 * t143; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t459, -t461, t299, qJ(1) * t299, 0, 0, t490, 0, t291, 0, -t494, t498, t198, pkin(5) * t465 + qJ(1) * t198 + t389 * t237, t230, -t207, t226, t229, -t501, t249, t389 * t179 + t388 * t183 + t209, t389 * t180 + t388 * t184 + t210, t389 * t199 + t388 * t200 - t238, qJ(1) * t170 + t389 * t154 + t388 * t161, t230, t226, t207, t249, t501, t229, t389 * t162 + t388 * t167 + t209, t389 * t177 + t388 * t178 - t238, t389 * t160 + t388 * t166 - t210, qJ(1) * t151 + t389 * t144 + t388 * t145, t230, t207, -t226, t229, -t501, t249, t389 * t152 + t388 * t155 + t209, t389 * t153 + t388 * t156 - t210, t389 * t163 + t388 * t164 + t238, qJ(1) * t148 + t389 * t142 + t388 * t143; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t418, t357, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t354 + t415, -pkin(1) * t353 + t442, 0, -pkin(1) * t247, t306, t406, t309, t305, t307, 0, t429 + t448, -t428 - t458, t219 + t427, pkin(1) * t195 + pkin(2) * t405 + pkin(6) * t219, t306, t309, -t406, 0, -t307, t305, qJ(4) * t456 + t394 * t221 + t429, t394 * t222 + t392 * t223 + t427, pkin(3) * t447 + t392 * t220 + t428, pkin(1) * t172 + pkin(6) * t190 + (-pkin(2) + t410) * t224, t306, -t406, -t309, t305, t307, 0, t394 * t188 + t392 * t304 + t429, t392 * t193 + t394 * t268 + t428, t394 * t194 + t392 * t201 - t427, pkin(1) * t157 - pkin(2) * t203 + pkin(6) * t176 + t394 * t159 + t392 * t174;];
tauB_reg = t1;
