% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRP3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:04
% EndTime: 2019-12-31 17:51:09
% DurationCPUTime: 3.38s
% Computational Cost: add. (5500->317), mult. (10699->412), div. (0->0), fcn. (5497->6), ass. (0->226)
t329 = sin(pkin(7));
t330 = cos(pkin(7));
t338 = qJD(1) ^ 2;
t296 = -t330 * qJDD(1) + t329 * t338;
t326 = g(3) - qJDD(2);
t267 = qJ(2) * t296 - t329 * t326;
t334 = sin(qJ(1));
t336 = cos(qJ(1));
t295 = t329 * qJDD(1) + t330 * t338;
t357 = t336 * t295 - t334 * t296;
t361 = -qJ(2) * t295 + t330 * t326;
t434 = -pkin(5) * t357 + t334 * t267 + t336 * t361;
t305 = t334 * g(1) - t336 * g(2);
t287 = qJDD(1) * pkin(1) + t305;
t306 = t336 * g(1) + t334 * g(2);
t288 = -t338 * pkin(1) - t306;
t229 = -t330 * t287 + t329 * t288;
t230 = t329 * t287 + t330 * t288;
t358 = t329 * t229 + t330 * t230;
t194 = t330 * t229 - t329 * t230;
t384 = t336 * t194;
t431 = -t334 * t358 + t384;
t394 = t334 * t194;
t153 = t336 * t358 + t394;
t355 = t334 * t295 + t336 * t296;
t415 = pkin(5) * t355 + t336 * t267 - t334 * t361;
t430 = pkin(2) + pkin(6);
t372 = qJDD(1) * qJ(3);
t341 = -t338 * pkin(2) + t230 + t372;
t373 = (qJD(3) * qJD(1));
t369 = 2 * t373;
t218 = t341 + t369;
t328 = qJDD(1) * pkin(2);
t219 = -t338 * qJ(3) + qJDD(3) + t229 - t328;
t178 = t329 * t218 - t330 * t219;
t359 = t330 * t218 + t329 * t219;
t143 = -t334 * t178 + t336 * t359;
t142 = t336 * t178 + t334 * t359;
t337 = qJD(4) ^ 2;
t333 = sin(qJ(4));
t324 = t333 ^ 2;
t403 = t324 * t338;
t308 = -t337 - t403;
t335 = cos(qJ(4));
t385 = t335 * t338;
t368 = t333 * t385;
t304 = qJDD(4) - t368;
t386 = t335 * t304;
t254 = t333 * t308 + t386;
t395 = t333 * t304;
t258 = t335 * t308 - t395;
t421 = pkin(3) * t254 - qJ(3) * t258;
t325 = t335 ^ 2;
t402 = t325 * t338;
t310 = -t337 - t402;
t303 = qJDD(4) + t368;
t396 = t333 * t303;
t256 = t335 * t310 - t396;
t387 = t335 * t303;
t261 = -t333 * t310 - t387;
t420 = pkin(3) * t256 - qJ(3) * t261;
t217 = -qJDD(1) * pkin(6) + t219;
t206 = t335 * t217;
t374 = qJD(1) * qJD(5);
t419 = -qJDD(4) * pkin(4) + 0.2e1 * t335 * t374 - t206;
t375 = qJD(1) * qJD(4);
t364 = t335 * t375;
t371 = t333 * qJDD(1);
t292 = -t364 - t371;
t376 = qJD(1) * t335;
t302 = qJD(4) * pkin(4) - qJ(5) * t376;
t378 = -t333 * t217 + t335 * t326;
t340 = t292 * qJ(5) - qJD(4) * t302 - 0.2e1 * t333 * t374 - t378;
t377 = t324 + t325;
t300 = t377 * t338;
t413 = pkin(3) * t300;
t291 = 0.2e1 * t364 + t371;
t213 = -t330 * t254 + t329 * t291;
t215 = t329 * t254 + t330 * t291;
t171 = t336 * t213 + t334 * t215;
t412 = pkin(5) * t171;
t365 = t333 * t375;
t370 = t335 * qJDD(1);
t294 = -0.2e1 * t365 + t370;
t214 = -t330 * t256 + t329 * t294;
t216 = t329 * t256 + t330 * t294;
t172 = t336 * t214 + t334 * t216;
t411 = pkin(5) * t172;
t297 = t377 * qJDD(1);
t400 = t330 * t297;
t243 = -t329 * t300 + t400;
t401 = t329 * t297;
t244 = -t330 * t300 - t401;
t197 = t336 * t243 + t334 * t244;
t410 = pkin(5) * t197;
t409 = qJ(2) * t213;
t408 = qJ(2) * t214;
t407 = qJ(2) * t243;
t293 = -t365 + t370;
t351 = pkin(4) * t385 - t326;
t175 = -t293 * qJ(5) + (-qJ(5) * t375 - t351) * t333 - t419;
t398 = t333 * t175;
t331 = t338 * pkin(6);
t212 = t218 - t331;
t397 = t333 * t212;
t389 = t335 * t175;
t388 = t335 * t212;
t363 = -pkin(1) * t258 + qJ(2) * t215;
t362 = -pkin(1) * t261 + qJ(2) * t216;
t199 = -t333 * t326 - t206;
t250 = -t334 * t305 - t336 * t306;
t353 = t329 * t368;
t352 = t330 * t368;
t350 = -pkin(1) * t295 - t230;
t299 = t336 * qJDD(1) - t334 * t338;
t349 = -pkin(5) * t299 - t334 * g(3);
t156 = -t335 * t199 - t333 * t378;
t157 = t333 * t199 - t335 * t378;
t249 = t336 * t305 - t334 * t306;
t348 = pkin(1) * t213 + qJ(3) * t291 - t430 * t254;
t347 = pkin(1) * t214 + qJ(3) * t294 - t430 * t256;
t346 = pkin(1) * t243 - qJ(3) * t300 + t430 * t297;
t344 = pkin(3) * t291 - t258 * t430;
t343 = pkin(3) * t294 - t261 * t430;
t342 = -pkin(1) * t296 - t229;
t339 = -t292 * pkin(4) - qJ(5) * t403 + qJDD(5) - t331 + t341;
t186 = (t302 * t335 + (2 * qJD(3))) * qJD(1) + t339;
t309 = t337 - t402;
t307 = -t337 + t403;
t301 = (-t324 + t325) * t338;
t298 = t334 * qJDD(1) + t336 * t338;
t285 = t377 * t375;
t273 = -pkin(5) * t298 + t336 * g(3);
t272 = -pkin(3) * t297 - pkin(4) * t370;
t265 = -t333 * t293 - t325 * t375;
t264 = -t335 * t292 - t324 * t375;
t263 = t330 * qJDD(4) - t329 * t285;
t262 = t329 * qJDD(4) + t330 * t285;
t260 = -t333 * t309 + t386;
t259 = (t293 - t365) * t335;
t257 = t335 * t307 - t396;
t255 = -t335 * t309 - t395;
t253 = -t333 * t307 - t387;
t252 = (-t292 + t364) * t333;
t251 = -pkin(4) * t294 - qJ(5) * t303;
t233 = qJ(2) * t244;
t232 = -t335 * t291 - t333 * t294;
t231 = t333 * t291 - t335 * t294;
t227 = -t329 * t264 - t352;
t226 = -t329 * t265 + t352;
t225 = t330 * t264 - t353;
t224 = t330 * t265 + t353;
t223 = -t329 * t255 + t330 * t370;
t222 = -t329 * t253 - t330 * t371;
t221 = t330 * t255 + t329 * t370;
t220 = t330 * t253 - t329 * t371;
t204 = -t329 * t231 + t330 * t301;
t203 = t330 * t231 + t329 * t301;
t202 = -t334 * t262 + t336 * t263;
t201 = t336 * t262 + t334 * t263;
t198 = -t334 * t243 + t336 * t244;
t196 = pkin(5) * t198;
t191 = pkin(1) * t326 + qJ(2) * t358;
t190 = -t334 * t225 + t336 * t227;
t189 = -t334 * t224 + t336 * t226;
t188 = t336 * t225 + t334 * t227;
t187 = t336 * t224 + t334 * t226;
t185 = -t334 * t221 + t336 * t223;
t184 = -t334 * t220 + t336 * t222;
t183 = t336 * t221 + t334 * t223;
t182 = t336 * t220 + t334 * t222;
t181 = -qJ(5) * t310 + t186;
t176 = -pkin(4) * t403 + t340;
t174 = -t334 * t214 + t336 * t216;
t173 = -t334 * t213 + t336 * t215;
t170 = pkin(5) * t174;
t169 = pkin(5) * t173;
t168 = t351 * t333 + (t293 + t365 + t370) * qJ(5) + t419;
t167 = -pkin(4) * t291 + qJ(5) * t308 - t302 * t376 - t339 - (2 * t373);
t166 = -t334 * t203 + t336 * t204;
t165 = t336 * t203 + t334 * t204;
t164 = t378 + t420;
t163 = -t199 + t421;
t162 = -qJ(5) * t371 + (t300 - t403) * pkin(4) + t340;
t161 = -qJ(2) * t178 + (-pkin(2) * t329 + qJ(3) * t330) * t326;
t160 = qJ(2) * t359 + (pkin(2) * t330 + qJ(3) * t329 + pkin(1)) * t326;
t159 = t344 + t388;
t158 = t343 - t397;
t155 = -t157 - t413;
t154 = (t310 + t403) * pkin(4) - t340 + t420;
t151 = pkin(4) * t304 + t175 + t421;
t150 = t329 * t156 + t330 * t212;
t149 = -t330 * t156 + t329 * t212;
t148 = -pkin(4) * t186 + qJ(5) * t176;
t147 = -pkin(3) * t400 - t329 * t155 - t407;
t146 = -pkin(3) * t401 + t330 * t155 + t233;
t145 = qJ(5) * t395 - t335 * t167 + t344;
t144 = -t333 * t181 - t335 * t251 + t343;
t141 = t335 * t176 - t398;
t140 = t333 * t176 + t389;
t139 = -t335 * t162 - t333 * t168 - t413;
t138 = pkin(3) * t156 - qJ(3) * t157;
t137 = -t329 * t158 + t330 * t164 - t408;
t136 = -t329 * t159 + t330 * t163 - t409;
t135 = t329 * t140 + t330 * t186;
t134 = -t330 * t140 + t329 * t186;
t133 = pkin(3) * t212 - t157 * t430;
t132 = -t329 * t139 + t330 * t272 - t407;
t131 = t330 * t139 + t329 * t272 + t233;
t130 = t330 * t158 + t329 * t164 + t362;
t129 = t330 * t159 + t329 * t163 + t363;
t128 = -t334 * t149 + t336 * t150;
t127 = t336 * t149 + t334 * t150;
t126 = -t329 * t144 + t330 * t154 - t408;
t125 = -t329 * t145 + t330 * t151 - t409;
t124 = t330 * t144 + t329 * t154 + t362;
t123 = t330 * t145 + t329 * t151 + t363;
t122 = pkin(3) * t140 + pkin(4) * t175 - qJ(3) * t141;
t121 = -t334 * t134 + t336 * t135;
t120 = t336 * t134 + t334 * t135;
t119 = -qJ(2) * t149 - t329 * t133 + t330 * t138;
t118 = pkin(3) * t186 + qJ(5) * t398 - t141 * t430 - t335 * t148;
t117 = -pkin(1) * t157 + qJ(2) * t150 + t330 * t133 + t329 * t138;
t116 = -qJ(2) * t134 - t329 * t118 + t330 * t122;
t115 = -pkin(1) * t141 + qJ(2) * t135 + t330 * t118 + t329 * t122;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t298, -t299, 0, t250, 0, 0, 0, 0, 0, 0, -t357, t355, 0, t153, 0, 0, 0, 0, 0, 0, 0, t357, -t355, t143, 0, 0, 0, 0, 0, 0, t173, t174, t198, t128, 0, 0, 0, 0, 0, 0, t173, t174, t198, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t299, -t298, 0, t249, 0, 0, 0, 0, 0, 0, -t355, -t357, 0, -t431, 0, 0, 0, 0, 0, 0, 0, t355, t357, t142, 0, 0, 0, 0, 0, 0, t171, t172, t197, t127, 0, 0, 0, 0, 0, 0, t171, t172, t197, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t326, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t326, 0, 0, 0, 0, 0, 0, t258, t261, 0, t157, 0, 0, 0, 0, 0, 0, t258, t261, 0, t141; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t299, 0, -t298, 0, t349, -t273, -t249, -pkin(5) * t249, 0, 0, -t355, 0, -t357, 0, t415, -t434, t431, pkin(5) * t431 + qJ(2) * t384 - t334 * t191, 0, t355, t357, 0, 0, 0, -t142, -t415, t434, -pkin(5) * t142 - t334 * t160 + t336 * t161, t189, t166, t185, t190, t184, t202, -t334 * t129 + t336 * t136 - t412, -t334 * t130 + t336 * t137 - t411, -t334 * t146 + t336 * t147 - t410, -pkin(5) * t127 - t334 * t117 + t336 * t119, t189, t166, t185, t190, t184, t202, -t334 * t123 + t336 * t125 - t412, -t334 * t124 + t336 * t126 - t411, -t334 * t131 + t336 * t132 - t410, -pkin(5) * t120 - t334 * t115 + t336 * t116; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t298, 0, t299, 0, t273, t349, t250, pkin(5) * t250, 0, 0, t357, 0, -t355, 0, t434, t415, t153, pkin(5) * t153 + qJ(2) * t394 + t336 * t191, 0, -t357, t355, 0, 0, 0, t143, -t434, -t415, pkin(5) * t143 + t336 * t160 + t334 * t161, t187, t165, t183, t188, t182, t201, t336 * t129 + t334 * t136 + t169, t336 * t130 + t334 * t137 + t170, t336 * t146 + t334 * t147 + t196, pkin(5) * t128 + t336 * t117 + t334 * t119, t187, t165, t183, t188, t182, t201, t336 * t123 + t334 * t125 + t169, t336 * t124 + t334 * t126 + t170, t336 * t131 + t334 * t132 + t196, pkin(5) * t121 + t336 * t115 + t334 * t116; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t305, t306, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t342, t350, 0, -pkin(1) * t194, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t328 - t342, -t350 + t369 + 0.2e1 * t372, pkin(1) * t178 - pkin(2) * t219 + qJ(3) * t218, t259, t232, t260, t252, t257, 0, t348 + t397, t347 + t388, -t156 + t346, pkin(1) * t149 + qJ(3) * t212 - t156 * t430, t259, t232, t260, t252, t257, 0, -qJ(5) * t386 - t333 * t167 + t348, t335 * t181 - t333 * t251 + t347, -t333 * t162 + t335 * t168 + t346, pkin(1) * t134 + qJ(3) * t186 - qJ(5) * t389 - t140 * t430 - t333 * t148;];
tauB_reg = t1;
