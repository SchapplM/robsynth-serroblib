% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRPRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRPRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:10
% EndTime: 2019-12-05 15:36:20
% DurationCPUTime: 6.14s
% Computational Cost: add. (9289->380), mult. (17201->520), div. (0->0), fcn. (10626->8), ass. (0->277)
t392 = qJD(4) ^ 2;
t388 = sin(qJ(4));
t379 = t388 ^ 2;
t393 = qJD(2) ^ 2;
t441 = t379 * t393;
t363 = t392 + t441;
t390 = cos(qJ(4));
t368 = t388 * t393 * t390;
t362 = qJDD(4) - t368;
t432 = t390 * t362;
t321 = -t388 * t363 + t432;
t426 = qJD(2) * qJD(4);
t418 = t390 * t426;
t423 = t388 * qJDD(2);
t349 = 0.2e1 * t418 + t423;
t383 = sin(pkin(8));
t385 = cos(pkin(8));
t268 = t321 * t383 + t349 * t385;
t271 = t321 * t385 - t349 * t383;
t389 = sin(qJ(2));
t391 = cos(qJ(2));
t225 = t268 * t389 - t271 * t391;
t434 = t388 * t362;
t313 = t363 * t390 + t434;
t384 = sin(pkin(7));
t386 = cos(pkin(7));
t192 = t225 * t384 + t313 * t386;
t500 = qJ(1) * t192;
t195 = t225 * t386 - t313 * t384;
t499 = qJ(1) * t195;
t221 = t268 * t391 + t271 * t389;
t498 = pkin(5) * t221;
t497 = -pkin(1) * t221 - pkin(2) * t268 - pkin(6) * t321;
t496 = -pkin(1) * t313 - pkin(5) * t225;
t449 = t349 * t388;
t419 = t388 * t426;
t422 = t390 * qJDD(2);
t400 = 0.2e1 * t419 - t422;
t472 = t400 * t390;
t290 = t472 + t449;
t380 = t390 ^ 2;
t358 = (t379 - t380) * t393;
t260 = t290 * t383 + t358 * t385;
t262 = t290 * t385 - t358 * t383;
t214 = t260 * t389 - t262 * t391;
t288 = t349 * t390 - t388 * t400;
t495 = t214 * t384 - t288 * t386;
t494 = t214 * t386 + t288 * t384;
t440 = t380 * t393;
t365 = -t392 + t440;
t319 = -t390 * t365 + t434;
t274 = t319 * t383 + t385 * t422;
t277 = t319 * t385 - t383 * t422;
t231 = t274 * t389 - t277 * t391;
t324 = t388 * t365 + t432;
t493 = t231 * t384 - t324 * t386;
t492 = t231 * t386 + t324 * t384;
t491 = qJ(3) * t268;
t489 = -pkin(2) * t313 + qJ(3) * t271;
t488 = t260 * t391 + t262 * t389;
t487 = t274 * t391 + t277 * t389;
t352 = t383 * qJDD(2) + t385 * t393;
t353 = t385 * qJDD(2) - t383 * t393;
t411 = -t352 * t389 + t391 * t353;
t486 = t384 * t411;
t485 = t386 * t411;
t360 = t386 * g(1) + t384 * g(2);
t381 = g(3) - qJDD(1);
t332 = -t391 * t360 - t389 * t381;
t326 = -t393 * pkin(2) + t332;
t331 = -t389 * t360 + t391 * t381;
t397 = qJDD(2) * pkin(2) - t331;
t255 = t383 * t326 - t385 * t397;
t256 = t385 * t326 + t383 * t397;
t412 = t255 * t383 + t385 * t256;
t201 = t255 * t385 - t256 * t383;
t455 = t201 * t389;
t173 = t391 * t412 + t455;
t454 = t201 * t391;
t172 = -t389 * t412 + t454;
t359 = t384 * g(1) - t386 * g(2);
t351 = -qJDD(3) + t359;
t301 = qJ(3) * t352 - t351 * t385;
t403 = -qJ(3) * t353 - t351 * t383;
t226 = -pkin(5) * t411 + t301 * t389 + t391 * t403;
t484 = 2 * qJD(5);
t481 = pkin(3) * t313;
t480 = pkin(6) * t313;
t467 = t391 * t352 + t353 * t389;
t227 = pkin(5) * t467 + t301 * t391 - t389 * t403;
t342 = t386 * t359;
t303 = -t360 * t384 + t342;
t456 = qJ(5) * t388;
t464 = pkin(4) * t390;
t406 = -t456 - t464;
t348 = t406 * qJD(2);
t251 = -pkin(3) * t393 + qJDD(2) * pkin(6) + t256;
t430 = -t390 * t251 + t388 * t351;
t404 = t390 * qJD(2) * t348 + qJDD(4) * qJ(5) + (qJD(4) * t484) - t430;
t427 = qJD(2) * t388;
t466 = t348 * t427 + qJDD(5);
t366 = -t392 - t440;
t361 = qJDD(4) + t368;
t435 = t388 * t361;
t317 = t390 * t366 - t435;
t267 = t317 * t383 - t385 * t400;
t270 = t317 * t385 + t383 * t400;
t220 = t267 * t391 + t270 * t389;
t463 = pkin(5) * t220;
t428 = t379 + t380;
t354 = t428 * qJDD(2);
t357 = t428 * t393;
t299 = t354 * t383 + t357 * t385;
t302 = t354 * t385 - t357 * t383;
t246 = t299 * t391 + t302 * t389;
t462 = pkin(5) * t246;
t433 = t390 * t361;
t312 = t366 * t388 + t433;
t461 = pkin(6) * t312;
t223 = -t267 * t389 + t270 * t391;
t190 = t223 * t384 - t312 * t386;
t460 = qJ(1) * t190;
t247 = -t299 * t389 + t302 * t391;
t459 = qJ(1) * t247;
t458 = qJ(3) * t267;
t457 = qJ(3) * t299;
t250 = -qJDD(2) * pkin(3) - t393 * pkin(6) + t255;
t453 = t250 * t388;
t452 = t250 * t390;
t443 = t359 * t384;
t356 = t391 * qJDD(2) - t389 * t393;
t439 = t384 * t356;
t438 = t384 * t381;
t242 = t386 * t247;
t437 = t386 * t356;
t436 = t386 * t381;
t431 = -pkin(1) * t312 + pkin(5) * t223;
t238 = t388 * t251 + t390 * t351;
t429 = t357 - t392;
t424 = t386 * qJDD(2);
t417 = -pkin(2) * t312 + qJ(3) * t270;
t416 = pkin(1) * t467 + pkin(2) * t352 - qJ(1) * t411 + t256;
t415 = -pkin(1) * t411 - pkin(2) * t353 - qJ(1) * t467 + t255;
t355 = t389 * qJDD(2) + t391 * t393;
t414 = -pkin(1) * t355 + qJ(1) * t356 - t332;
t413 = pkin(1) * t356 + qJ(1) * t355 - t331;
t265 = t331 * t389 + t391 * t332;
t304 = -t386 * t360 - t443;
t409 = t383 * t368;
t408 = t385 * t368;
t212 = -pkin(3) * t312 + t238;
t405 = pkin(4) * t388 - qJ(5) * t390;
t310 = pkin(5) * t355 - t359 * t391;
t309 = -pkin(5) * t356 - t359 * t389;
t184 = t238 * t390 + t388 * t430;
t185 = t238 * t388 - t390 * t430;
t264 = t331 * t391 - t332 * t389;
t402 = t418 + t423;
t401 = -t419 + t422;
t399 = -pkin(1) * t220 - pkin(2) * t267 - pkin(6) * t317;
t398 = -qJDD(4) * pkin(4) + t238 + t466;
t396 = -t401 * pkin(4) + t250 + (-t402 - t418) * qJ(5);
t395 = -pkin(1) * t246 - pkin(2) * t299 - pkin(3) * t357 - pkin(6) * t354;
t394 = t427 * t484 - t396;
t373 = t384 * qJDD(2);
t364 = t392 - t441;
t347 = t405 * qJDD(2);
t346 = t428 * t426;
t340 = t386 * t355;
t339 = t384 * t355;
t330 = -t379 * t426 + t390 * t402;
t329 = -t380 * t426 - t388 * t401;
t328 = t383 * qJDD(4) + t385 * t346;
t327 = -t385 * qJDD(4) + t383 * t346;
t320 = -t388 * t364 + t433;
t315 = -t364 * t390 - t435;
t291 = qJ(3) * t302;
t284 = t386 * t467;
t283 = t384 * t467;
t282 = t330 * t385 - t409;
t281 = t329 * t385 + t409;
t280 = t330 * t383 + t408;
t279 = t329 * t383 - t408;
t278 = t320 * t385 + t383 * t423;
t275 = t320 * t383 - t385 * t423;
t258 = -t327 * t389 + t328 * t391;
t257 = t327 * t391 + t328 * t389;
t254 = t386 * t258;
t253 = t384 * t258;
t245 = t265 * t386 - t443;
t244 = t265 * t384 + t342;
t243 = pkin(5) * t247;
t241 = t384 * t247;
t240 = qJ(1) * t242;
t237 = -t280 * t389 + t282 * t391;
t236 = -t279 * t389 + t281 * t391;
t235 = t280 * t391 + t282 * t389;
t234 = t279 * t391 + t281 * t389;
t233 = -t275 * t389 + t278 * t391;
t230 = t275 * t391 + t278 * t389;
t229 = t452 + t480;
t228 = t453 - t461;
t216 = t392 * qJ(5) - t398;
t213 = -t430 + t481;
t211 = -pkin(4) * t392 + t404;
t210 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t427 + t396;
t209 = t237 * t386 + t384 * t449;
t208 = t236 * t386 - t384 * t472;
t207 = t237 * t384 - t386 * t449;
t206 = t236 * t384 + t386 * t472;
t205 = t233 * t386 - t315 * t384;
t204 = t233 * t384 + t315 * t386;
t203 = t429 * qJ(5) + t398;
t198 = t429 * pkin(4) + t404;
t197 = (-t400 - t419) * pkin(4) + t394;
t196 = -pkin(4) * t419 + qJ(5) * t349 + t394;
t193 = t223 * t386 + t312 * t384;
t189 = qJ(1) * t193;
t188 = pkin(2) * t351 + qJ(3) * t412;
t187 = (-t366 - t392) * qJ(5) + (-qJDD(4) - t361) * pkin(4) + t212 + t466;
t186 = -t481 - qJ(5) * t362 + (-t363 + t392) * pkin(4) - t404;
t182 = -pkin(4) * t449 + t196 * t390 - t480;
t181 = -qJ(5) * t472 - t197 * t388 - t461;
t180 = t184 * t385 - t457;
t179 = t184 * t383 + t291;
t178 = t211 * t390 - t216 * t388;
t177 = t211 * t388 + t216 * t390;
t176 = t185 * t385 + t250 * t383;
t175 = t185 * t383 - t250 * t385;
t174 = -t198 * t388 + t203 * t390;
t170 = pkin(3) * t349 - t453 - t497;
t169 = pkin(3) * t400 + t399 + t452;
t168 = -t213 * t383 + t229 * t385 + t491;
t167 = -t212 * t383 + t228 * t385 - t458;
t166 = t173 * t386 - t351 * t384;
t165 = t173 * t384 + t351 * t386;
t164 = t213 * t385 + t229 * t383 - t489;
t163 = t212 * t385 + t228 * t383 + t417;
t162 = -t185 + t395;
t161 = t174 * t385 - t347 * t383 - t457;
t160 = t174 * t383 + t347 * t385 + t291;
t159 = -t388 * t196 + (-pkin(3) - t464) * t349 + t497;
t158 = -t390 * t197 - (-pkin(3) - t456) * t400 + t399;
t157 = t178 * t385 + t210 * t383;
t156 = t178 * t383 - t210 * t385;
t155 = -t198 * t390 - t203 * t388 + t395;
t154 = pkin(1) * t172 + pkin(2) * t201;
t153 = t181 * t385 - t187 * t383 - t458;
t152 = t182 * t385 - t186 * t383 - t491;
t151 = t181 * t383 + t187 * t385 + t417;
t150 = t182 * t383 + t186 * t385 + t489;
t149 = -pkin(3) * t177 - pkin(4) * t216 - qJ(5) * t211;
t148 = -pkin(6) * t177 + t405 * t210;
t147 = -t179 * t389 + t180 * t391 - t462;
t146 = -t175 * t389 + t176 * t391;
t145 = t175 * t391 + t176 * t389;
t144 = pkin(5) * t172 + qJ(3) * t454 - t188 * t389;
t143 = -qJ(3) * t175 - (pkin(3) * t383 - pkin(6) * t385) * t184;
t142 = -t164 * t389 + t168 * t391 + t498;
t141 = -t163 * t389 + t167 * t391 - t463;
t140 = t146 * t386 - t184 * t384;
t139 = t146 * t384 + t184 * t386;
t138 = -t160 * t389 + t161 * t391 - t462;
t137 = -t156 * t389 + t157 * t391;
t136 = t156 * t391 + t157 * t389;
t135 = qJ(3) * t176 - (-pkin(3) * t385 - pkin(6) * t383 - pkin(2)) * t184;
t134 = -t151 * t389 + t153 * t391 - t463;
t133 = -t150 * t389 + t152 * t391 - t498;
t132 = t137 * t386 + t177 * t384;
t131 = t137 * t384 - t177 * t386;
t130 = -pkin(1) * t145 - pkin(2) * t175 + pkin(3) * t250 - pkin(6) * t185;
t129 = -qJ(3) * t156 + t148 * t385 - t149 * t383;
t128 = -pkin(2) * t177 + qJ(3) * t157 + t148 * t383 + t149 * t385;
t127 = -pkin(1) * t136 - pkin(2) * t156 - pkin(6) * t178 + (pkin(3) - t406) * t210;
t126 = -pkin(5) * t145 - t135 * t389 + t143 * t391;
t125 = -pkin(5) * t136 - t128 * t389 + t129 * t391;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t304, 0, 0, 0, 0, 0, 0, -t340, -t437, 0, t245, 0, 0, 0, 0, 0, 0, -t284, -t485, 0, t166, 0, 0, 0, 0, 0, 0, t193, t195, t242, t140, 0, 0, 0, 0, 0, 0, t193, t242, -t195, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, 0, 0, 0, 0, 0, 0, -t339, -t439, 0, t244, 0, 0, 0, 0, 0, 0, -t283, -t486, 0, t165, 0, 0, 0, 0, 0, 0, t190, t192, t241, t139, 0, 0, 0, 0, 0, 0, t190, t241, -t192, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t381, 0, 0, 0, 0, 0, 0, t356, -t355, 0, -t264, 0, 0, 0, 0, 0, 0, t411, -t467, 0, -t172, 0, 0, 0, 0, 0, 0, t220, -t221, t246, t145, 0, 0, 0, 0, 0, 0, t220, t246, t221, t136; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t438, -t436, -t303, -qJ(1) * t303, 0, 0, t437, 0, -t340, t373, t386 * t309 + t413 * t384, t386 * t310 + t414 * t384, t386 * t264, -qJ(1) * t244 - (pkin(1) * t384 - pkin(5) * t386) * t264, 0, 0, t485, 0, -t284, t373, t386 * t226 - t384 * t415, t386 * t227 - t384 * t416, t386 * t172, -qJ(1) * t165 + t144 * t386 - t154 * t384, t209, t494, t205, t208, t492, t254, t141 * t386 - t169 * t384 - t460, t142 * t386 - t170 * t384 - t500, t386 * t147 + (-t162 - t459) * t384, -qJ(1) * t139 + t126 * t386 - t130 * t384, t209, t205, -t494, t254, -t492, t208, t134 * t386 - t158 * t384 - t460, t386 * t138 + (-t155 - t459) * t384, t133 * t386 - t159 * t384 + t500, -qJ(1) * t131 + t125 * t386 - t127 * t384; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t436, -t438, t304, qJ(1) * t304, 0, 0, t439, 0, -t339, -t424, t384 * t309 - t413 * t386, t384 * t310 - t414 * t386, t384 * t264, qJ(1) * t245 - (-pkin(1) * t386 - pkin(5) * t384) * t264, 0, 0, t486, 0, -t283, -t424, t384 * t226 + t386 * t415, t384 * t227 + t386 * t416, t384 * t172, qJ(1) * t166 + t144 * t384 + t154 * t386, t207, t495, t204, t206, t493, t253, t141 * t384 + t169 * t386 + t189, t142 * t384 + t170 * t386 + t499, t147 * t384 + t162 * t386 + t240, qJ(1) * t140 + t126 * t384 + t130 * t386, t207, t204, -t495, t253, -t493, t206, t134 * t384 + t158 * t386 + t189, t138 * t384 + t155 * t386 + t240, t133 * t384 + t159 * t386 - t499, qJ(1) * t132 + t125 * t384 + t127 * t386; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t359, t360, 0, 0, 0, 0, t355, 0, t356, 0, -t310, t309, t265, pkin(1) * t359 + pkin(5) * t265, 0, 0, t467, 0, t411, 0, -t227, t226, t173, pkin(1) * t351 + pkin(5) * t173 + qJ(3) * t455 + t188 * t391, t235, -t488, t230, t234, -t487, t257, t163 * t391 + t167 * t389 + t431, t164 * t391 + t168 * t389 - t496, t179 * t391 + t180 * t389 + t243, pkin(1) * t184 + pkin(5) * t146 + t135 * t391 + t143 * t389, t235, t230, t488, t257, t487, t234, t151 * t391 + t153 * t389 + t431, t160 * t391 + t161 * t389 + t243, t150 * t391 + t152 * t389 + t496, -pkin(1) * t177 + pkin(5) * t137 + t128 * t391 + t129 * t389;];
tauB_reg = t1;
