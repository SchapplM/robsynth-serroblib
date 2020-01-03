% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:18
% EndTime: 2019-12-31 17:45:25
% DurationCPUTime: 5.17s
% Computational Cost: add. (11334->372), mult. (22883->544), div. (0->0), fcn. (14149->8), ass. (0->246)
t365 = sin(pkin(7));
t374 = qJD(1) ^ 2;
t367 = cos(pkin(7));
t400 = t367 * qJDD(1);
t338 = t365 * t374 - t400;
t362 = g(3) - qJDD(2);
t307 = qJ(2) * t338 - t365 * t362;
t370 = sin(qJ(1));
t372 = cos(qJ(1));
t401 = t365 * qJDD(1);
t337 = t367 * t374 + t401;
t389 = t372 * t337 - t370 * t338;
t394 = -qJ(2) * t337 + t367 * t362;
t457 = -pkin(5) * t389 + t370 * t307 + t372 * t394;
t369 = sin(qJ(5));
t364 = sin(pkin(8));
t366 = cos(pkin(8));
t371 = cos(qJ(5));
t381 = t364 * t371 + t366 * t369;
t321 = t381 * qJD(1);
t323 = (-t364 * t369 + t366 * t371) * qJD(1);
t432 = t323 * t321;
t453 = qJDD(5) - t432;
t459 = t369 * t453;
t458 = t371 * t453;
t343 = t370 * g(1) - t372 * g(2);
t332 = qJDD(1) * pkin(1) + t343;
t344 = t372 * g(1) + t370 * g(2);
t333 = -t374 * pkin(1) - t344;
t282 = -t367 * t332 + t365 * t333;
t283 = t365 * t332 + t367 * t333;
t390 = t365 * t282 + t367 * t283;
t237 = t367 * t282 - t365 * t283;
t415 = t372 * t237;
t454 = -t370 * t390 + t415;
t422 = t370 * t237;
t202 = t372 * t390 + t422;
t387 = t370 * t337 + t372 * t338;
t436 = pkin(5) * t387 + t372 * t307 - t370 * t394;
t452 = -0.2e1 * t364;
t378 = 0.2e1 * qJD(3) * qJD(1) + t283;
t377 = -t374 * pkin(2) + t378;
t404 = qJDD(1) * qJ(3);
t267 = t377 + t404;
t363 = qJDD(1) * pkin(2);
t270 = -t374 * qJ(3) + qJDD(3) + t282 - t363;
t220 = t365 * t267 - t367 * t270;
t391 = t367 * t267 + t365 * t270;
t174 = -t370 * t220 + t372 * t391;
t173 = t372 * t220 + t370 * t391;
t358 = t364 ^ 2;
t359 = t366 ^ 2;
t408 = t358 + t359;
t442 = t408 * t374;
t330 = t364 * t442;
t296 = t367 * t330 + t364 * t401;
t299 = -t365 * t330 + t364 * t400;
t250 = t372 * t296 + t370 * t299;
t449 = t370 * t296 - t372 * t299;
t266 = -qJDD(1) * qJ(4) + t270;
t405 = qJD(1) * qJD(4);
t249 = t364 * t266 - t366 * t362 + t405 * t452;
t379 = t381 * qJDD(1);
t441 = -t374 * qJ(4) + qJDD(4);
t317 = t321 ^ 2;
t318 = t323 ^ 2;
t435 = pkin(1) * t337;
t434 = pkin(2) + qJ(4);
t431 = t358 * t374;
t263 = t366 * t266;
t396 = t366 * t405;
t348 = -0.2e1 * t396;
t402 = qJDD(1) * t366;
t427 = t366 * t374;
t229 = -pkin(6) * t402 + t263 + t348 + (-pkin(4) * t427 + t362) * t364;
t403 = qJDD(1) * t364;
t234 = -pkin(4) * t431 - pkin(6) * t403 + t249;
t195 = -t371 * t229 + t369 * t234;
t196 = t369 * t229 + t371 * t234;
t161 = -t371 * t195 + t369 * t196;
t430 = t364 * t161;
t334 = t408 * qJDD(1);
t429 = t365 * t334;
t428 = t366 * t161;
t426 = t367 * t334;
t241 = (pkin(4) * t364 + qJ(3)) * qJDD(1) + t377 + (-t359 * t374 - t431) * pkin(6) + t441;
t424 = t369 * t241;
t274 = qJDD(5) + t432;
t423 = t369 * t274;
t417 = t371 * t241;
t416 = t371 * t274;
t409 = t358 - t359;
t407 = t321 * qJD(5);
t406 = t323 * qJD(5);
t399 = t365 * t432;
t398 = t367 * t432;
t395 = t366 * t400;
t162 = t369 * t195 + t371 * t196;
t392 = -t364 * t362 - t263;
t303 = -t370 * t343 - t372 * t344;
t264 = t267 + t441;
t385 = t264 + t404;
t342 = t372 * qJDD(1) - t370 * t374;
t384 = -pkin(5) * t342 - t370 * g(3);
t320 = -t369 * t403 + t371 * t402;
t248 = t392 + 0.2e1 * t396;
t210 = -t366 * t248 + t364 * t249;
t211 = t364 * t248 + t366 * t249;
t304 = (t365 * t427 - t395) * t364;
t305 = t337 * t366 * t364;
t383 = t372 * t304 + t370 * t305;
t382 = t370 * t304 - t372 * t305;
t302 = t372 * t343 - t370 * t344;
t380 = -pkin(1) * t338 - t282;
t373 = qJD(5) ^ 2;
t341 = t370 * qJDD(1) + t372 * t374;
t340 = t409 * t374;
t335 = t409 * qJDD(1);
t329 = t366 * t442;
t316 = -pkin(5) * t341 + t372 * g(3);
t314 = -t318 - t373;
t313 = -t318 + t373;
t312 = t317 - t373;
t301 = -t365 * t329 + t395;
t298 = t367 * t329 + t366 * t401;
t287 = -t367 * t442 - t429;
t286 = -t365 * t335 - t367 * t340;
t285 = -t365 * t442 + t426;
t284 = t367 * t335 - t365 * t340;
t280 = t318 - t317;
t279 = t320 - t407;
t278 = t320 - 0.2e1 * t407;
t277 = -t379 - t406;
t276 = 0.2e1 * t406 + t379;
t272 = -t373 - t317;
t269 = (-t321 * t371 + t323 * t369) * qJD(5);
t268 = (-t321 * t369 - t323 * t371) * qJD(5);
t261 = -t317 - t318;
t259 = pkin(3) * t403 + t366 * t264;
t258 = pkin(3) * t402 - t364 * t264;
t257 = t371 * t279 - t369 * t406;
t256 = t369 * t279 + t371 * t406;
t255 = -t369 * t277 + t371 * t407;
t254 = t371 * t277 + t369 * t407;
t253 = -t370 * t298 + t372 * t301;
t251 = t372 * t298 + t370 * t301;
t247 = -t369 * t314 - t416;
t246 = -t369 * t313 + t458;
t245 = t371 * t312 - t423;
t244 = t371 * t314 - t423;
t243 = t371 * t313 + t459;
t242 = t369 * t312 + t416;
t240 = -t370 * t285 + t372 * t287;
t239 = t372 * t285 + t370 * t287;
t233 = -t371 * t276 - t369 * t278;
t232 = t369 * t320 - t371 * t379;
t231 = -t369 * t276 + t371 * t278;
t230 = -t371 * t320 - t369 * t379;
t228 = -pkin(3) * t330 + t348 - t392;
t227 = -pkin(3) * t329 - t249;
t226 = t371 * t272 - t459;
t225 = t369 * t272 + t458;
t223 = pkin(1) * t362 + qJ(2) * t390;
t218 = -t366 * t268 - t364 * t269;
t217 = t367 * qJDD(5) - t365 * t218;
t216 = t365 * qJDD(5) + t367 * t218;
t215 = -qJ(2) * t220 + (-pkin(2) * t365 + qJ(3) * t367) * t362;
t214 = qJ(2) * t391 + (pkin(2) * t367 + qJ(3) * t365 + pkin(1)) * t362;
t213 = -t366 * t256 - t364 * t257;
t212 = -t366 * t254 - t364 * t255;
t209 = -t364 * t244 + t366 * t247;
t208 = t366 * t244 + t364 * t247;
t207 = -t366 * t243 - t364 * t246;
t206 = -t366 * t242 - t364 * t245;
t205 = -pkin(6) * t244 + t417;
t204 = -pkin(3) * t442 - t211;
t203 = -pkin(6) * t225 + t424;
t200 = -t365 * t207 + t367 * t320;
t199 = -t365 * t206 - t367 * t379;
t198 = t367 * t207 + t365 * t320;
t197 = t367 * t206 - t365 * t379;
t193 = -t364 * t230 + t366 * t232;
t192 = -t366 * t231 - t364 * t233;
t191 = t366 * t230 + t364 * t232;
t190 = -t364 * t225 + t366 * t226;
t189 = t366 * t225 + t364 * t226;
t188 = -t365 * t213 + t398;
t187 = -t365 * t212 - t398;
t186 = t367 * t213 + t399;
t185 = t367 * t212 - t399;
t184 = -pkin(4) * t278 + pkin(6) * t247 + t424;
t183 = -qJ(2) * t296 + t367 * t228 - t365 * t259;
t182 = -qJ(2) * t298 + t367 * t227 - t365 * t258;
t181 = qJ(2) * t299 + t365 * t228 + t367 * t259;
t180 = qJ(2) * t301 + t365 * t227 + t367 * t258;
t179 = t365 * t208 + t367 * t278;
t178 = -t367 * t208 + t365 * t278;
t177 = -pkin(4) * t276 + pkin(6) * t226 - t417;
t176 = t365 * t210 + t367 * t264;
t175 = -t367 * t210 + t365 * t264;
t172 = -pkin(3) * t426 - qJ(2) * t285 - t365 * t204;
t171 = -pkin(3) * t429 + qJ(2) * t287 + t367 * t204;
t170 = -t365 * t192 + t367 * t280;
t169 = t367 * t192 + t365 * t280;
t168 = t365 * t189 + t367 * t276;
t167 = -t367 * t189 + t365 * t276;
t166 = t365 * t191 + t367 * t261;
t165 = -t367 * t191 + t365 * t261;
t164 = pkin(3) * t210 - qJ(3) * t211;
t163 = pkin(3) * t264 - t434 * t211;
t160 = -t370 * t178 + t372 * t179;
t159 = t372 * t178 + t370 * t179;
t158 = -t370 * t175 + t372 * t176;
t157 = t372 * t175 + t370 * t176;
t156 = -pkin(4) * t241 + pkin(6) * t162;
t155 = -pkin(6) * t230 - t161;
t154 = pkin(3) * t191 + pkin(4) * t230 - qJ(3) * t193;
t153 = -t370 * t167 + t372 * t168;
t152 = t372 * t167 + t370 * t168;
t151 = -t370 * t165 + t372 * t166;
t150 = t372 * t165 + t370 * t166;
t149 = -pkin(4) * t261 + pkin(6) * t232 + t162;
t148 = pkin(3) * t208 + pkin(4) * t244 - qJ(3) * t209 - t196;
t147 = pkin(3) * t189 + pkin(4) * t225 - qJ(3) * t190 - t195;
t146 = pkin(3) * t278 - t366 * t184 - t364 * t205 - t434 * t209;
t145 = t366 * t162 - t430;
t144 = t364 * t162 + t428;
t143 = pkin(3) * t276 - t366 * t177 - t434 * t190 - t364 * t203;
t142 = t365 * t144 + t367 * t241;
t141 = -t367 * t144 + t365 * t241;
t140 = -qJ(2) * t175 - t365 * t163 + t367 * t164;
t139 = -pkin(1) * t211 + qJ(2) * t176 + t367 * t163 + t365 * t164;
t138 = pkin(3) * t261 - t366 * t149 - t364 * t155 - t434 * t193;
t137 = -qJ(2) * t178 - t365 * t146 + t367 * t148;
t136 = -pkin(1) * t209 + qJ(2) * t179 + t367 * t146 + t365 * t148;
t135 = -t370 * t141 + t372 * t142;
t134 = t372 * t141 + t370 * t142;
t133 = -qJ(2) * t167 - t365 * t143 + t367 * t147;
t132 = pkin(3) * t144 + pkin(4) * t161 - qJ(3) * t145;
t131 = -pkin(1) * t190 + qJ(2) * t168 + t367 * t143 + t365 * t147;
t130 = -qJ(2) * t165 - t365 * t138 + t367 * t154;
t129 = -pkin(1) * t193 + qJ(2) * t166 + t367 * t138 + t365 * t154;
t128 = pkin(3) * t241 + pkin(6) * t430 - t434 * t145 - t366 * t156;
t127 = -qJ(2) * t141 - t365 * t128 + t367 * t132;
t126 = -pkin(1) * t145 + qJ(2) * t142 + t367 * t128 + t365 * t132;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t341, -t342, 0, t303, 0, 0, 0, 0, 0, 0, -t389, t387, 0, t202, 0, 0, 0, 0, 0, 0, 0, t389, -t387, t174, 0, 0, 0, 0, 0, 0, -t449, t253, t240, t158, 0, 0, 0, 0, 0, 0, t153, t160, t151, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t342, -t341, 0, t302, 0, 0, 0, 0, 0, 0, -t387, -t389, 0, -t454, 0, 0, 0, 0, 0, 0, 0, t387, t389, t173, 0, 0, 0, 0, 0, 0, t250, t251, t239, t157, 0, 0, 0, 0, 0, 0, t152, t159, t150, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t362, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t362, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, 0, 0, 0, 0, 0, t190, t209, t193, t145; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t342, 0, -t341, 0, t384, -t316, -t302, -pkin(5) * t302, 0, 0, -t387, 0, -t389, 0, t436, -t457, t454, pkin(5) * t454 + qJ(2) * t415 - t370 * t223, 0, t387, t389, 0, 0, 0, -t173, -t436, t457, -pkin(5) * t173 - t370 * t214 + t372 * t215, -t382, -t370 * t284 + t372 * t286, t253, t382, t449, 0, -pkin(5) * t250 - t370 * t181 + t372 * t183, -pkin(5) * t251 - t370 * t180 + t372 * t182, -pkin(5) * t239 - t370 * t171 + t372 * t172, -pkin(5) * t157 - t370 * t139 + t372 * t140, -t370 * t186 + t372 * t188, -t370 * t169 + t372 * t170, -t370 * t198 + t372 * t200, -t370 * t185 + t372 * t187, -t370 * t197 + t372 * t199, -t370 * t216 + t372 * t217, -pkin(5) * t152 - t370 * t131 + t372 * t133, -pkin(5) * t159 - t370 * t136 + t372 * t137, -pkin(5) * t150 - t370 * t129 + t372 * t130, -pkin(5) * t134 - t370 * t126 + t372 * t127; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t341, 0, t342, 0, t316, t384, t303, pkin(5) * t303, 0, 0, t389, 0, -t387, 0, t457, t436, t202, pkin(5) * t202 + qJ(2) * t422 + t372 * t223, 0, -t389, t387, 0, 0, 0, t174, -t457, -t436, pkin(5) * t174 + t372 * t214 + t370 * t215, t383, t372 * t284 + t370 * t286, t251, -t383, -t250, 0, -pkin(5) * t449 + t372 * t181 + t370 * t183, pkin(5) * t253 + t372 * t180 + t370 * t182, pkin(5) * t240 + t372 * t171 + t370 * t172, pkin(5) * t158 + t372 * t139 + t370 * t140, t372 * t186 + t370 * t188, t372 * t169 + t370 * t170, t372 * t198 + t370 * t200, t372 * t185 + t370 * t187, t372 * t197 + t370 * t199, t372 * t216 + t370 * t217, pkin(5) * t153 + t372 * t131 + t370 * t133, pkin(5) * t160 + t372 * t136 + t370 * t137, pkin(5) * t151 + t372 * t129 + t370 * t130, pkin(5) * t135 + t372 * t126 + t370 * t127; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t343, t344, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t380, -t283 - t435, 0, -pkin(1) * t237, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t363 - t380, t378 + 0.2e1 * t404 + t435, pkin(1) * t220 - pkin(2) * t270 + qJ(3) * t267, t359 * qJDD(1), t402 * t452, 0, t358 * qJDD(1), 0, 0, pkin(1) * t296 + t434 * t330 + t385 * t364, pkin(1) * t298 + t434 * t329 + t385 * t366, pkin(1) * t285 - qJ(3) * t442 + t434 * t334 - t210, pkin(1) * t175 + qJ(3) * t264 - t434 * t210, -t364 * t256 + t366 * t257, -t364 * t231 + t366 * t233, -t364 * t243 + t366 * t246, -t364 * t254 + t366 * t255, -t364 * t242 + t366 * t245, -t364 * t268 + t366 * t269, pkin(1) * t167 + qJ(3) * t276 - t364 * t177 - t434 * t189 + t366 * t203, pkin(1) * t178 + qJ(3) * t278 - t364 * t184 + t366 * t205 - t434 * t208, pkin(1) * t165 + qJ(3) * t261 - t364 * t149 + t366 * t155 - t434 * t191, pkin(1) * t141 - pkin(6) * t428 + qJ(3) * t241 - t434 * t144 - t364 * t156;];
tauB_reg = t1;
