% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:47
% EndTime: 2022-01-23 09:34:55
% DurationCPUTime: 7.27s
% Computational Cost: add. (33385->386), mult. (48058->556), div. (0->0), fcn. (28562->10), ass. (0->257)
t414 = qJD(1) + qJD(3);
t409 = qJD(4) + t414;
t406 = t409 ^ 2;
t425 = cos(qJ(4));
t413 = qJDD(1) + qJDD(3);
t407 = qJDD(4) + t413;
t421 = sin(qJ(4));
t459 = t421 * t407;
t366 = t425 * t406 + t459;
t450 = t425 * t407;
t369 = t421 * t406 - t450;
t422 = sin(qJ(3));
t426 = cos(qJ(3));
t310 = t426 * t366 - t422 * t369;
t417 = g(3) - qJDD(2);
t347 = pkin(7) * t366 - t425 * t417;
t486 = pkin(7) * t369 - t421 * t417;
t259 = pkin(6) * t310 + t426 * t347 - t422 * t486;
t315 = t422 * t366 + t426 * t369;
t418 = sin(pkin(9));
t419 = cos(pkin(9));
t268 = t419 * t310 - t418 * t315;
t499 = pkin(6) * t315 + t422 * t347 + t426 * t486;
t203 = qJ(2) * t268 + t419 * t259 - t418 * t499;
t272 = t418 * t310 + t419 * t315;
t423 = sin(qJ(1));
t427 = cos(qJ(1));
t227 = t423 * t268 + t427 * t272;
t511 = qJ(2) * t272 + t418 * t259 + t419 * t499;
t520 = pkin(5) * t227 + t423 * t203 + t427 * t511;
t500 = t427 * t268 - t423 * t272;
t519 = pkin(5) * t500 + t427 * t203 - t423 * t511;
t398 = t423 * g(1) - t427 * g(2);
t386 = qJDD(1) * pkin(1) + t398;
t399 = t427 * g(1) + t423 * g(2);
t429 = qJD(1) ^ 2;
t387 = -t429 * pkin(1) - t399;
t330 = -t419 * t386 + t418 * t387;
t326 = qJDD(1) * pkin(2) - t330;
t331 = t418 * t386 + t419 * t387;
t327 = -t429 * pkin(2) + t331;
t287 = t422 * t326 + t426 * t327;
t412 = t414 ^ 2;
t279 = -t412 * pkin(3) + t287;
t431 = t426 * t326 - t422 * t327;
t430 = t413 * pkin(3) + t431;
t235 = t421 * t279 - t425 * t430;
t236 = t425 * t279 + t421 * t430;
t441 = t421 * t235 + t425 * t236;
t194 = t425 * t235 - t421 * t236;
t448 = t426 * t194;
t168 = -t422 * t441 + t448;
t457 = t422 * t194;
t490 = t426 * t441 + t457;
t155 = t419 * t168 - t418 * t490;
t506 = t418 * t168 + t419 * t490;
t141 = t423 * t155 + t427 * t506;
t516 = t427 * t155 - t423 * t506;
t377 = t426 * t412 + t422 * t413;
t380 = t422 * t412 - t426 * t413;
t319 = t419 * t377 - t418 * t380;
t353 = pkin(6) * t377 - t426 * t417;
t487 = pkin(6) * t380 - t422 * t417;
t267 = qJ(2) * t319 + t419 * t353 - t418 * t487;
t323 = t418 * t377 + t419 * t380;
t282 = t423 * t319 + t427 * t323;
t498 = qJ(2) * t323 + t418 * t353 + t419 * t487;
t513 = pkin(5) * t282 + t423 * t267 + t427 * t498;
t484 = t427 * t319 - t423 * t323;
t512 = pkin(5) * t484 + t427 * t267 - t423 * t498;
t440 = t426 * t287 - t422 * t431;
t240 = -t422 * t287 - t426 * t431;
t464 = t419 * t240;
t198 = -t418 * t440 + t464;
t465 = t418 * t240;
t489 = t419 * t440 + t465;
t173 = t423 * t198 + t427 * t489;
t505 = t427 * t198 - t423 * t489;
t439 = t418 * t330 + t419 * t331;
t290 = t419 * t330 - t418 * t331;
t446 = t427 * t290;
t491 = -t423 * t439 + t446;
t455 = t423 * t290;
t243 = t427 * t439 + t455;
t388 = t418 * qJDD(1) + t419 * t429;
t389 = t419 * qJDD(1) - t418 * t429;
t338 = -t423 * t388 + t427 * t389;
t358 = qJ(2) * t388 - t419 * t417;
t432 = -qJ(2) * t389 - t418 * t417;
t488 = -pkin(5) * t338 + t423 * t358 + t427 * t432;
t470 = t427 * t388 + t423 * t389;
t482 = pkin(5) * t470 + t427 * t358 - t423 * t432;
t468 = pkin(1) * t417;
t467 = pkin(2) * t417;
t420 = sin(qJ(5));
t415 = t420 ^ 2;
t466 = t415 * t406;
t232 = -t407 * pkin(4) - t406 * pkin(8) + t235;
t463 = t420 * t232;
t424 = cos(qJ(5));
t394 = t424 * t406 * t420;
t383 = qJDD(5) + t394;
t462 = t420 * t383;
t384 = qJDD(5) - t394;
t461 = t420 * t384;
t460 = t420 * t407;
t452 = t424 * t232;
t451 = t424 * t384;
t400 = t424 * t407;
t233 = -t406 * pkin(4) + t407 * pkin(8) + t236;
t224 = t424 * t233 - t420 * t417;
t416 = t424 ^ 2;
t445 = t415 + t416;
t444 = qJD(5) * t409;
t443 = t420 * t444;
t442 = t424 * t444;
t223 = t420 * t233 + t424 * t417;
t188 = t420 * t223 + t424 * t224;
t349 = -t423 * t398 - t427 * t399;
t435 = t421 * t394;
t434 = t425 * t394;
t396 = t427 * qJDD(1) - t423 * t429;
t433 = -pkin(5) * t396 - t423 * g(3);
t187 = t424 * t223 - t420 * t224;
t348 = t427 * t398 - t423 * t399;
t428 = qJD(5) ^ 2;
t401 = t416 * t406;
t395 = t423 * qJDD(1) + t427 * t429;
t393 = -t401 - t428;
t392 = t401 - t428;
t391 = -t428 - t466;
t390 = t428 - t466;
t373 = t424 * t383;
t372 = -pkin(5) * t395 + t427 * g(3);
t371 = t401 - t466;
t370 = t401 + t466;
t365 = t445 * t407;
t363 = t400 - 0.2e1 * t443;
t362 = t400 - t443;
t361 = t442 + t460;
t360 = 0.2e1 * t442 + t460;
t359 = t445 * t444;
t343 = t421 * qJDD(5) + t425 * t359;
t342 = -t425 * qJDD(5) + t421 * t359;
t337 = -t420 * t391 - t451;
t336 = -t420 * t390 + t373;
t335 = t424 * t393 - t462;
t334 = t424 * t392 - t461;
t333 = t424 * t391 - t461;
t332 = t420 * t393 + t373;
t329 = t424 * t361 - t415 * t444;
t328 = -t420 * t362 - t416 * t444;
t313 = t425 * t365 - t421 * t370;
t309 = t421 * t365 + t425 * t370;
t308 = -t420 * t360 + t424 * t363;
t307 = t425 * t336 + t420 * t459;
t306 = t425 * t334 + t421 * t400;
t305 = t421 * t336 - t420 * t450;
t304 = t421 * t334 - t424 * t450;
t303 = t425 * t329 - t435;
t302 = t425 * t328 + t435;
t301 = t421 * t329 + t434;
t300 = t421 * t328 - t434;
t299 = t425 * t337 + t421 * t360;
t298 = t425 * t335 - t421 * t363;
t297 = t421 * t337 - t425 * t360;
t296 = t421 * t335 + t425 * t363;
t295 = -t422 * t342 + t426 * t343;
t294 = t426 * t342 + t422 * t343;
t293 = t425 * t308 - t421 * t371;
t292 = t421 * t308 + t425 * t371;
t285 = qJ(2) * t439 + t468;
t275 = -t422 * t309 + t426 * t313;
t274 = t426 * t309 + t422 * t313;
t263 = -t422 * t305 + t426 * t307;
t262 = -t422 * t304 + t426 * t306;
t261 = t426 * t305 + t422 * t307;
t260 = t426 * t304 + t422 * t306;
t255 = -t422 * t301 + t426 * t303;
t254 = -t422 * t300 + t426 * t302;
t253 = t426 * t301 + t422 * t303;
t252 = t426 * t300 + t422 * t302;
t251 = -t422 * t297 + t426 * t299;
t250 = -t422 * t296 + t426 * t298;
t249 = t426 * t297 + t422 * t299;
t248 = t426 * t296 + t422 * t298;
t247 = -t418 * t294 + t419 * t295;
t246 = t419 * t294 + t418 * t295;
t245 = -t422 * t292 + t426 * t293;
t244 = t426 * t292 + t422 * t293;
t237 = pkin(6) * t440 + t467;
t230 = -t418 * t274 + t419 * t275;
t229 = t419 * t274 + t418 * t275;
t221 = -t418 * t261 + t419 * t263;
t220 = -t418 * t260 + t419 * t262;
t219 = t419 * t261 + t418 * t263;
t218 = t419 * t260 + t418 * t262;
t217 = -t418 * t253 + t419 * t255;
t216 = -t418 * t252 + t419 * t254;
t215 = t419 * t253 + t418 * t255;
t214 = t419 * t252 + t418 * t254;
t213 = -pkin(8) * t333 + t452;
t212 = -pkin(8) * t332 + t463;
t211 = -pkin(4) * t333 + t224;
t210 = -pkin(4) * t332 + t223;
t209 = -t418 * t249 + t419 * t251;
t208 = -t418 * t248 + t419 * t250;
t207 = t419 * t249 + t418 * t251;
t206 = t419 * t248 + t418 * t250;
t205 = -t418 * t244 + t419 * t245;
t204 = t419 * t244 + t418 * t245;
t191 = pkin(3) * t417 + pkin(7) * t441;
t190 = -t423 * t229 + t427 * t230;
t189 = t427 * t229 + t423 * t230;
t185 = -pkin(7) * t309 + t425 * t187;
t184 = pkin(7) * t313 + t421 * t187;
t183 = -t423 * t207 + t427 * t209;
t182 = -t423 * t206 + t427 * t208;
t181 = t427 * t207 + t423 * t209;
t180 = t427 * t206 + t423 * t208;
t179 = -pkin(7) * t297 - t421 * t211 + t425 * t213;
t178 = -pkin(7) * t296 - t421 * t210 + t425 * t212;
t177 = -pkin(3) * t333 + pkin(7) * t299 + t425 * t211 + t421 * t213;
t176 = -pkin(3) * t332 + pkin(7) * t298 + t425 * t210 + t421 * t212;
t175 = t425 * t188 + t421 * t232;
t174 = t421 * t188 - t425 * t232;
t171 = pkin(6) * t464 + qJ(2) * t198 - t418 * t237;
t170 = pkin(6) * t465 + qJ(2) * t489 + t419 * t237 + t468;
t165 = -pkin(6) * t274 - t422 * t184 + t426 * t185;
t164 = pkin(6) * t275 + t426 * t184 + t422 * t185;
t163 = -t422 * t174 + t426 * t175;
t162 = t426 * t174 + t422 * t175;
t161 = -pkin(6) * t249 - t422 * t177 + t426 * t179;
t160 = -pkin(6) * t248 - t422 * t176 + t426 * t178;
t159 = -pkin(2) * t333 + pkin(6) * t251 + t426 * t177 + t422 * t179;
t158 = -pkin(2) * t332 + pkin(6) * t250 + t426 * t176 + t422 * t178;
t157 = -pkin(7) * t174 - (pkin(4) * t421 - pkin(8) * t425) * t187;
t152 = pkin(6) * t168 + pkin(7) * t448 - t422 * t191;
t151 = pkin(6) * t490 + pkin(7) * t457 + t426 * t191 + t467;
t150 = pkin(7) * t175 - (-pkin(4) * t425 - pkin(8) * t421 - pkin(3)) * t187;
t149 = -qJ(2) * t229 - t418 * t164 + t419 * t165;
t148 = qJ(2) * t230 + t419 * t164 + t418 * t165;
t147 = -t418 * t162 + t419 * t163;
t146 = t419 * t162 + t418 * t163;
t145 = -qJ(2) * t207 - t418 * t159 + t419 * t161;
t144 = -qJ(2) * t206 - t418 * t158 + t419 * t160;
t143 = -pkin(1) * t333 + qJ(2) * t209 + t419 * t159 + t418 * t161;
t142 = -pkin(1) * t332 + qJ(2) * t208 + t419 * t158 + t418 * t160;
t139 = qJ(2) * t155 - t418 * t151 + t419 * t152;
t138 = qJ(2) * t506 + t419 * t151 + t418 * t152 + t468;
t137 = -pkin(6) * t162 - t422 * t150 + t426 * t157;
t136 = -t423 * t146 + t427 * t147;
t135 = t427 * t146 + t423 * t147;
t134 = pkin(2) * t187 + pkin(6) * t163 + t426 * t150 + t422 * t157;
t133 = -qJ(2) * t146 - t418 * t134 + t419 * t137;
t132 = pkin(1) * t187 + qJ(2) * t147 + t419 * t134 + t418 * t137;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t395, -t396, 0, t349, 0, 0, 0, 0, 0, 0, -t470, -t338, 0, t243, 0, 0, 0, 0, 0, 0, -t484, t282, 0, t173, 0, 0, 0, 0, 0, 0, -t500, t227, 0, t141, 0, 0, 0, 0, 0, 0, t182, t183, t190, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t396, -t395, 0, t348, 0, 0, 0, 0, 0, 0, t338, -t470, 0, -t491, 0, 0, 0, 0, 0, 0, -t282, -t484, 0, -t505, 0, 0, 0, 0, 0, 0, -t227, -t500, 0, -t516, 0, 0, 0, 0, 0, 0, t180, t181, t189, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, 0, 0, 0, 0, 0, 0, t332, t333, 0, -t187; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t396, 0, -t395, 0, t433, -t372, -t348, -pkin(5) * t348, 0, 0, t338, 0, -t470, 0, t488, t482, t491, pkin(5) * t491 + qJ(2) * t446 - t423 * t285, 0, 0, -t282, 0, -t484, 0, t513, t512, t505, pkin(5) * t505 - t423 * t170 + t427 * t171, 0, 0, -t227, 0, -t500, 0, t520, t519, t516, pkin(5) * t516 - t423 * t138 + t427 * t139, -t423 * t215 + t427 * t217, -t423 * t204 + t427 * t205, -t423 * t219 + t427 * t221, -t423 * t214 + t427 * t216, -t423 * t218 + t427 * t220, -t423 * t246 + t427 * t247, -pkin(5) * t180 - t423 * t142 + t427 * t144, -pkin(5) * t181 - t423 * t143 + t427 * t145, -pkin(5) * t189 - t423 * t148 + t427 * t149, -pkin(5) * t135 - t423 * t132 + t427 * t133; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t395, 0, t396, 0, t372, t433, t349, pkin(5) * t349, 0, 0, t470, 0, t338, 0, -t482, t488, t243, pkin(5) * t243 + qJ(2) * t455 + t427 * t285, 0, 0, t484, 0, -t282, 0, -t512, t513, t173, pkin(5) * t173 + t427 * t170 + t423 * t171, 0, 0, t500, 0, -t227, 0, -t519, t520, t141, pkin(5) * t141 + t427 * t138 + t423 * t139, t427 * t215 + t423 * t217, t427 * t204 + t423 * t205, t427 * t219 + t423 * t221, t427 * t214 + t423 * t216, t427 * t218 + t423 * t220, t427 * t246 + t423 * t247, pkin(5) * t182 + t427 * t142 + t423 * t144, pkin(5) * t183 + t427 * t143 + t423 * t145, pkin(5) * t190 + t427 * t148 + t423 * t149, pkin(5) * t136 + t427 * t132 + t423 * t133; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t398, t399, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t389 - t330, -pkin(1) * t388 - t331, 0, -pkin(1) * t290, 0, 0, 0, 0, 0, t413, -pkin(1) * t323 - pkin(2) * t380 + t431, -pkin(1) * t319 - pkin(2) * t377 - t287, 0, -pkin(1) * t198 - pkin(2) * t240, 0, 0, 0, 0, 0, t407, -pkin(1) * t272 - pkin(2) * t315 - pkin(3) * t369 - t235, -pkin(1) * t268 - pkin(2) * t310 - pkin(3) * t366 - t236, 0, -pkin(1) * t155 - pkin(2) * t168 - pkin(3) * t194, (t361 + t442) * t420, t424 * t360 + t420 * t363, t424 * t390 + t462, (t362 - t443) * t424, t420 * t392 + t451, 0, pkin(1) * t206 + pkin(2) * t248 + pkin(3) * t296 + pkin(4) * t363 + pkin(8) * t335 - t452, pkin(1) * t207 + pkin(2) * t249 + pkin(3) * t297 - pkin(4) * t360 + pkin(8) * t337 + t463, pkin(1) * t229 + pkin(2) * t274 + pkin(3) * t309 + pkin(4) * t370 + pkin(8) * t365 + t188, pkin(1) * t146 + pkin(2) * t162 + pkin(3) * t174 - pkin(4) * t232 + pkin(8) * t188;];
tauB_reg = t1;
