% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:58
% EndTime: 2019-12-31 19:53:04
% DurationCPUTime: 3.59s
% Computational Cost: add. (7645->326), mult. (10364->405), div. (0->0), fcn. (5316->6), ass. (0->230)
t373 = qJD(4) ^ 2;
t360 = qJD(1) + qJD(2);
t358 = t360 ^ 2;
t369 = cos(qJ(4));
t364 = t369 ^ 2;
t432 = t364 * t358;
t341 = t373 + t432;
t366 = sin(qJ(4));
t402 = t366 * t358 * t369;
t335 = qJDD(4) + t402;
t430 = t366 * t335;
t290 = t369 * t341 + t430;
t406 = qJD(4) * t360;
t399 = t366 * t406;
t359 = qJDD(1) + qJDD(2);
t416 = t369 * t359;
t323 = -0.2e1 * t399 + t416;
t367 = sin(qJ(2));
t370 = cos(qJ(2));
t246 = t370 * t290 + t367 * t323;
t250 = t367 * t290 - t370 * t323;
t368 = sin(qJ(1));
t371 = cos(qJ(1));
t205 = t371 * t246 - t368 * t250;
t476 = pkin(5) * t205;
t209 = t368 * t246 + t371 * t250;
t475 = pkin(5) * t209;
t474 = pkin(6) * t246;
t446 = -pkin(7) - pkin(2);
t473 = -pkin(1) * t246 + t446 * t290;
t418 = t369 * t335;
t294 = -t366 * t341 + t418;
t472 = -pkin(1) * t294 + pkin(6) * t250;
t398 = t369 * t406;
t428 = t366 * t359;
t320 = 0.2e1 * t398 + t428;
t266 = t366 * t320 - t369 * t323;
t363 = t366 ^ 2;
t331 = (-t363 + t364) * t358;
t232 = t370 * t266 + t367 * t331;
t235 = t367 * t266 - t370 * t331;
t471 = t371 * t232 - t368 * t235;
t470 = t368 * t232 + t371 * t235;
t433 = t363 * t358;
t339 = -t373 + t433;
t289 = t366 * t339 + t418;
t424 = t367 * t359;
t259 = t370 * t289 + t366 * t424;
t412 = t370 * t359;
t263 = t367 * t289 - t366 * t412;
t469 = t371 * t259 - t368 * t263;
t468 = t368 * t259 + t371 * t263;
t326 = t370 * t358 + t424;
t328 = t367 * t358 - t412;
t271 = t371 * t326 - t368 * t328;
t395 = pkin(6) * t326 - t370 * g(3);
t457 = -pkin(6) * t328 + t367 * g(3);
t459 = pkin(5) * t271 + t368 * t457 + t371 * t395;
t391 = t368 * t326 + t371 * t328;
t447 = pkin(5) * t391 + t368 * t395 - t371 * t457;
t343 = t368 * g(1) - t371 * g(2);
t333 = qJDD(1) * pkin(1) + t343;
t344 = t371 * g(1) + t368 * g(2);
t374 = qJD(1) ^ 2;
t334 = -t374 * pkin(1) - t344;
t281 = -t370 * t333 + t367 * t334;
t282 = t367 * t333 + t370 * t334;
t392 = t367 * t281 + t370 * t282;
t230 = t370 * t281 - t367 * t282;
t411 = t371 * t230;
t465 = -t368 * t392 + t411;
t423 = t368 * t230;
t190 = t371 * t392 + t423;
t463 = pkin(1) * t328;
t461 = t446 * t294;
t460 = -pkin(1) * t326 - t282;
t434 = t359 * qJ(3);
t379 = -t358 * pkin(2) + t282 + t434;
t407 = (qJD(3) * t360);
t404 = 2 * t407;
t253 = t379 + t404;
t349 = t359 * pkin(2);
t380 = qJDD(3) + t281 - t349;
t258 = -t358 * qJ(3) + t380;
t213 = t367 * t253 - t370 * t258;
t393 = t370 * t253 + t367 * t258;
t180 = -t368 * t213 + t371 * t393;
t179 = t371 * t213 + t368 * t393;
t456 = pkin(3) * t290 - qJ(3) * t294;
t340 = -t373 - t433;
t336 = qJDD(4) - t402;
t417 = t369 * t336;
t288 = t366 * t340 + t417;
t429 = t366 * t336;
t293 = t369 * t340 - t429;
t451 = pkin(3) * t288 - qJ(3) * t293;
t408 = t363 + t364;
t330 = t408 * t358;
t445 = pkin(3) * t330;
t444 = pkin(4) * t366;
t443 = pkin(4) * t369;
t245 = -t370 * t288 + t367 * t320;
t248 = t367 * t288 + t370 * t320;
t204 = t371 * t245 + t368 * t248;
t442 = pkin(5) * t204;
t325 = t408 * t359;
t414 = t370 * t325;
t269 = -t367 * t330 + t414;
t426 = t367 * t325;
t274 = -t370 * t330 - t426;
t225 = t371 * t269 + t368 * t274;
t441 = pkin(5) * t225;
t440 = pkin(6) * t245;
t439 = pkin(6) * t269;
t436 = qJ(5) * t366;
t385 = -qJ(5) * t369 + t444;
t435 = t385 * t358;
t351 = t358 * pkin(7);
t242 = t253 - t351;
t431 = t366 * t242;
t420 = t369 * t242;
t419 = t369 * t320;
t252 = -t359 * pkin(7) + t258;
t238 = t366 * g(3) + t369 * t252;
t405 = qJD(5) * t369;
t403 = t446 * t293;
t397 = -pkin(1) * t293 + pkin(6) * t248;
t394 = t369 * g(3) - t366 * t252;
t302 = -t368 * t343 - t371 * t344;
t389 = t367 * t402;
t388 = t370 * t402;
t338 = t371 * qJDD(1) - t368 * t374;
t387 = -pkin(5) * t338 - t368 * g(3);
t386 = t436 + t443;
t196 = t369 * t238 - t366 * t394;
t197 = -t366 * t238 - t369 * t394;
t384 = t366 * t323 + t419;
t383 = -t369 * t339 + t430;
t301 = t371 * t343 - t368 * t344;
t382 = pkin(1) * t245 + qJ(3) * t320 + t446 * t288;
t381 = pkin(1) * t269 - qJ(3) * t330 - t446 * t325;
t378 = qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t366 * t435 - t394;
t377 = -qJDD(4) * pkin(4) - t373 * qJ(5) + t369 * t435 + qJDD(5) - t238;
t321 = -t398 - t428;
t322 = -t399 + t416;
t376 = -t321 * pkin(4) - t322 * qJ(5) - t351 + t379;
t375 = 0.2e1 * t360 * t405 - t376 - (2 * t407);
t372 = pkin(1) * g(3);
t342 = t373 - t432;
t337 = t368 * qJDD(1) + t371 * t374;
t315 = -pkin(5) * t337 + t371 * g(3);
t314 = t408 * t406;
t300 = t370 * qJDD(4) - t367 * t314;
t299 = t367 * qJDD(4) + t370 * t314;
t298 = -t366 * t322 - t364 * t406;
t297 = -t369 * t321 - t363 * t406;
t295 = -t366 * t342 + t417;
t291 = -t369 * t342 - t429;
t284 = (t322 - t399) * t369;
t283 = (-t321 + t398) * t366;
t280 = -pkin(3) * t325 - t386 * t359;
t267 = pkin(6) * t274;
t264 = -t367 * t291 + t369 * t412;
t261 = t370 * t291 + t367 * t416;
t257 = -t367 * t297 - t388;
t256 = -t367 * t298 + t388;
t255 = t370 * t297 - t389;
t254 = t370 * t298 + t389;
t237 = -t368 * t299 + t371 * t300;
t236 = t371 * t299 + t368 * t300;
t227 = pkin(6) * t392 + t372;
t226 = -t368 * t269 + t371 * t274;
t224 = pkin(5) * t226;
t222 = -t373 * pkin(4) + t378;
t221 = -t368 * t261 + t371 * t264;
t220 = t371 * t261 + t368 * t264;
t219 = -t368 * t255 + t371 * t257;
t218 = -t368 * t254 + t371 * t256;
t217 = t371 * t255 + t368 * t257;
t216 = t371 * t254 + t368 * t256;
t211 = qJ(5) * t330 + t377;
t210 = (t330 - t373) * pkin(4) + t378;
t207 = -t368 * t245 + t371 * t248;
t203 = -pkin(6) * t213 + (-pkin(2) * t367 + qJ(3) * t370) * g(3);
t202 = pkin(5) * t207;
t201 = pkin(6) * t393 + t372 + (pkin(2) * t370 + qJ(3) * t367) * g(3);
t200 = (t386 * qJD(4) + (2 * qJD(3)) - 0.2e1 * t405) * t360 + t376;
t199 = t394 - t456;
t198 = t238 + t451;
t195 = -qJ(5) * t399 + (-t320 - t398) * pkin(4) + t375;
t194 = -pkin(4) * t398 + (t323 - t399) * qJ(5) + t375;
t193 = pkin(3) * t320 + t403 + t420;
t192 = pkin(3) * t323 - t431 - t461;
t191 = -t197 - t445;
t188 = pkin(4) * t336 + qJ(5) * t340 - t377 + t451;
t187 = qJ(5) * t335 + (t341 - t373) * pkin(4) + t378 + t456;
t186 = t367 * t196 + t370 * t242;
t185 = -t370 * t196 + t367 * t242;
t184 = t369 * t222 + t366 * t377;
t183 = t366 * t222 - t369 * t377;
t182 = -pkin(3) * t414 - t367 * t191 - t439;
t181 = -pkin(3) * t426 + t370 * t191 + t267;
t178 = -t369 * t210 - t366 * t211 - t445;
t177 = -t366 * t194 + (-pkin(3) - t443) * t323 + t461;
t176 = -t369 * t195 + (pkin(3) + t436) * t320 + t403;
t175 = pkin(3) * t196 - qJ(3) * t197;
t174 = pkin(3) * t242 + t446 * t197;
t173 = t367 * t183 + t370 * t200;
t172 = -t370 * t183 + t367 * t200;
t171 = -t367 * t178 + t370 * t280 - t439;
t170 = t370 * t178 + t367 * t280 + t267;
t169 = -t367 * t192 + t370 * t199 - t474;
t168 = -t367 * t193 + t370 * t198 - t440;
t167 = t370 * t192 + t367 * t199 - t472;
t166 = t370 * t193 + t367 * t198 + t397;
t165 = -t368 * t185 + t371 * t186;
t164 = t371 * t185 + t368 * t186;
t163 = -t367 * t176 + t370 * t188 - t440;
t162 = -t367 * t177 + t370 * t187 + t474;
t161 = t370 * t176 + t367 * t188 + t397;
t160 = t370 * t177 + t367 * t187 + t472;
t159 = pkin(3) * t183 - pkin(4) * t377 - qJ(3) * t184 + qJ(5) * t222;
t158 = -t368 * t172 + t371 * t173;
t157 = t371 * t172 + t368 * t173;
t156 = t446 * t184 + (pkin(3) + t386) * t200;
t155 = -pkin(6) * t185 - t367 * t174 + t370 * t175;
t154 = -pkin(1) * t197 + pkin(6) * t186 + t370 * t174 + t367 * t175;
t153 = -pkin(6) * t172 - t367 * t156 + t370 * t159;
t152 = -pkin(1) * t184 + pkin(6) * t173 + t370 * t156 + t367 * t159;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t337, -t338, 0, t302, 0, 0, 0, 0, 0, 0, -t271, t391, 0, t190, 0, 0, 0, 0, 0, 0, 0, t271, -t391, t180, 0, 0, 0, 0, 0, 0, t207, -t209, t226, t165, 0, 0, 0, 0, 0, 0, t207, t226, t209, t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t338, -t337, 0, t301, 0, 0, 0, 0, 0, 0, -t391, -t271, 0, -t465, 0, 0, 0, 0, 0, 0, 0, t391, t271, t179, 0, 0, 0, 0, 0, 0, t204, t205, t225, t164, 0, 0, 0, 0, 0, 0, t204, t225, -t205, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t293, -t294, 0, t197, 0, 0, 0, 0, 0, 0, t293, 0, t294, t184; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t338, 0, -t337, 0, t387, -t315, -t301, -pkin(5) * t301, 0, 0, -t391, 0, -t271, 0, t447, t459, t465, pkin(5) * t465 + pkin(6) * t411 - t368 * t227, 0, t391, t271, 0, 0, 0, -t179, -t447, -t459, -pkin(5) * t179 - t368 * t201 + t371 * t203, t218, -t470, t221, t219, t468, t237, -t368 * t166 + t371 * t168 - t442, -t368 * t167 + t371 * t169 - t476, -t368 * t181 + t371 * t182 - t441, -pkin(5) * t164 - t368 * t154 + t371 * t155, t218, t221, t470, t237, -t468, t219, -t368 * t161 + t371 * t163 - t442, -t368 * t170 + t371 * t171 - t441, -t368 * t160 + t371 * t162 + t476, -pkin(5) * t157 - t368 * t152 + t371 * t153; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t337, 0, t338, 0, t315, t387, t302, pkin(5) * t302, 0, 0, t271, 0, -t391, 0, -t459, t447, t190, pkin(5) * t190 + pkin(6) * t423 + t371 * t227, 0, -t271, t391, 0, 0, 0, t180, t459, -t447, pkin(5) * t180 + t371 * t201 + t368 * t203, t216, t471, t220, t217, -t469, t236, t371 * t166 + t368 * t168 + t202, t371 * t167 + t368 * t169 - t475, t371 * t181 + t368 * t182 + t224, pkin(5) * t165 + t371 * t154 + t368 * t155, t216, t220, -t471, t236, t469, t217, t371 * t161 + t368 * t163 + t202, t371 * t170 + t368 * t171 + t224, t371 * t160 + t368 * t162 + t475, pkin(5) * t158 + t371 * t152 + t368 * t153; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t343, t344, 0, 0, 0, 0, 0, 0, 0, t359, -t281 - t463, t460, 0, -pkin(1) * t230, t359, 0, 0, 0, 0, 0, 0, -t349 + t380 + t463, t404 + 0.2e1 * t434 - t460, pkin(1) * t213 - pkin(2) * t258 + qJ(3) * t253, t284, -t384, t295, t283, -t383, 0, t382 + t431, qJ(3) * t323 + t420 - t473, -t196 + t381, pkin(1) * t185 + qJ(3) * t242 + t446 * t196, t284, t295, t384, 0, t383, t283, -qJ(5) * t419 - t366 * t195 + t382, -t366 * t210 + t369 * t211 + t381, t369 * t194 + (-qJ(3) - t444) * t323 + t473, pkin(1) * t172 + t446 * t183 + (qJ(3) + t385) * t200;];
tauB_reg = t1;
