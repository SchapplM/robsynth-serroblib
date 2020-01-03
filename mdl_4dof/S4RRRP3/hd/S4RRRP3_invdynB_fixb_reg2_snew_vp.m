% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRP3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:17
% EndTime: 2019-12-31 17:14:21
% DurationCPUTime: 2.84s
% Computational Cost: add. (5243->278), mult. (7735->358), div. (0->0), fcn. (4216->6), ass. (0->206)
t329 = qJD(3) ^ 2;
t319 = qJD(1) + qJD(2);
t317 = t319 ^ 2;
t323 = sin(qJ(3));
t321 = t323 ^ 2;
t378 = t321 * t317;
t300 = t329 + t378;
t326 = cos(qJ(3));
t305 = t326 * t317 * t323;
t297 = qJDD(3) - t305;
t364 = t326 * t297;
t255 = -t323 * t300 + t364;
t356 = qJD(3) * t319;
t349 = t326 * t356;
t318 = qJDD(1) + qJDD(2);
t370 = t323 * t318;
t278 = 0.2e1 * t349 + t370;
t324 = sin(qJ(2));
t327 = cos(qJ(2));
t210 = t324 * t255 + t327 * t278;
t213 = t327 * t255 - t324 * t278;
t325 = sin(qJ(1));
t328 = cos(qJ(1));
t169 = t328 * t210 + t325 * t213;
t411 = pkin(4) * t169;
t173 = t325 * t210 - t328 * t213;
t410 = pkin(4) * t173;
t409 = pkin(5) * t210;
t408 = pkin(1) * t210 + pkin(6) * t255;
t373 = t323 * t297;
t249 = t326 * t300 + t373;
t407 = -pkin(1) * t249 + pkin(5) * t213;
t350 = t323 * t356;
t363 = t326 * t318;
t279 = -0.2e1 * t350 + t363;
t244 = t326 * t279;
t245 = t323 * t278;
t232 = -t244 + t245;
t322 = t326 ^ 2;
t293 = (t321 - t322) * t317;
t201 = t324 * t232 + t327 * t293;
t203 = t327 * t232 - t324 * t293;
t406 = t328 * t201 + t325 * t203;
t405 = t325 * t201 - t328 * t203;
t377 = t322 * t317;
t302 = -t329 + t377;
t253 = -t326 * t302 + t373;
t360 = t327 * t318;
t224 = t324 * t253 + t326 * t360;
t227 = t327 * t253 - t324 * t363;
t404 = t328 * t224 + t325 * t227;
t403 = t325 * t224 - t328 * t227;
t367 = t324 * t318;
t288 = t327 * t317 + t367;
t291 = t324 * t317 - t360;
t239 = t325 * t288 + t328 * t291;
t268 = pkin(5) * t288 - t327 * g(3);
t392 = pkin(5) * t291 - t324 * g(3);
t402 = pkin(4) * t239 + t325 * t268 + t328 * t392;
t339 = t328 * t288 - t325 * t291;
t401 = pkin(4) * t339 + t328 * t268 - t325 * t392;
t307 = t328 * g(1) + t325 * g(2);
t330 = qJD(1) ^ 2;
t295 = -t330 * pkin(1) - t307;
t306 = t325 * g(1) - t328 * g(2);
t337 = qJDD(1) * pkin(1) + t306;
t242 = t324 * t295 - t327 * t337;
t243 = t327 * t295 + t324 * t337;
t346 = t324 * t242 + t327 * t243;
t198 = t327 * t242 - t324 * t243;
t359 = t328 * t198;
t397 = -t325 * t346 + t359;
t366 = t325 * t198;
t162 = t328 * t346 + t366;
t396 = 2 * qJD(4);
t394 = pkin(2) * t249;
t393 = pkin(6) * t249;
t387 = pkin(3) * t326;
t342 = -qJ(4) * t323 - t387;
t277 = t342 * t319;
t230 = -t317 * pkin(2) + t318 * pkin(6) + t243;
t347 = t323 * g(3) - t326 * t230;
t334 = t326 * t319 * t277 + qJDD(3) * qJ(4) + (qJD(3) * t396) - t347;
t389 = t323 * t302 + t364;
t379 = t319 * t323;
t388 = t277 * t379 + qJDD(4);
t303 = -t329 - t377;
t296 = qJDD(3) + t305;
t374 = t323 * t296;
t252 = t326 * t303 - t374;
t209 = t324 * t252 + t327 * t279;
t212 = t327 * t252 - t324 * t279;
t168 = t328 * t209 + t325 * t212;
t386 = pkin(4) * t168;
t357 = t321 + t322;
t284 = t357 * t318;
t292 = t357 * t317;
t235 = t324 * t284 + t327 * t292;
t238 = t327 * t284 - t324 * t292;
t189 = t328 * t235 + t325 * t238;
t385 = pkin(4) * t189;
t384 = pkin(5) * t209;
t383 = pkin(5) * t235;
t285 = t326 * t296;
t247 = t323 * t303 + t285;
t382 = pkin(6) * t247;
t229 = -t318 * pkin(2) - t317 * pkin(6) + t242;
t376 = t323 * t229;
t375 = t323 * t279;
t365 = t326 * t229;
t215 = t326 * g(3) + t323 * t230;
t358 = t292 - t329;
t352 = pkin(1) * t209 + pkin(2) * t279 + pkin(6) * t252;
t351 = pkin(1) * t235 + pkin(2) * t292 + pkin(6) * t284;
t348 = -pkin(1) * t247 + pkin(5) * t212;
t176 = t323 * t215 - t326 * t347;
t262 = -t325 * t306 - t328 * t307;
t345 = t324 * t305;
t344 = t327 * t305;
t191 = -pkin(2) * t247 + t215;
t299 = t328 * qJDD(1) - t325 * t330;
t343 = -pkin(4) * t299 - t325 * g(3);
t341 = pkin(3) * t323 - qJ(4) * t326;
t175 = t326 * t215 + t323 * t347;
t340 = t326 * t278 + t375;
t261 = t328 * t306 - t325 * t307;
t336 = t349 + t370;
t335 = -t350 + t363;
t333 = -qJDD(3) * pkin(3) + t215 + t388;
t332 = -t335 * pkin(3) + t229 + (-t336 - t349) * qJ(4);
t331 = t379 * t396 - t332;
t301 = t329 - t378;
t298 = t325 * qJDD(1) + t328 * t330;
t275 = -pkin(4) * t298 + t328 * g(3);
t274 = t341 * t318;
t273 = t357 * t356;
t260 = t324 * qJDD(3) + t327 * t273;
t259 = -t327 * qJDD(3) + t324 * t273;
t258 = -t321 * t356 + t326 * t336;
t257 = -t322 * t356 - t323 * t335;
t254 = -t323 * t301 + t285;
t248 = t326 * t301 + t374;
t233 = pkin(5) * t238;
t228 = t327 * t254 + t323 * t367;
t225 = t324 * t254 - t323 * t360;
t222 = t327 * t258 - t345;
t221 = t327 * t257 + t345;
t220 = t324 * t258 + t344;
t219 = t324 * t257 - t344;
t205 = -t325 * t259 + t328 * t260;
t204 = t328 * t259 + t325 * t260;
t195 = pkin(1) * g(3) + pkin(5) * t346;
t194 = t365 + t393;
t193 = t376 - t382;
t192 = -t347 + t394;
t190 = -t325 * t235 + t328 * t238;
t188 = t329 * qJ(4) - t333;
t187 = pkin(4) * t190;
t186 = -t329 * pkin(3) + t334;
t185 = -t325 * t225 + t328 * t228;
t184 = t328 * t225 + t325 * t228;
t183 = -t325 * t220 + t328 * t222;
t182 = -t325 * t219 + t328 * t221;
t181 = t328 * t220 + t325 * t222;
t180 = t328 * t219 + t325 * t221;
t179 = t358 * qJ(4) + t333;
t178 = t358 * pkin(3) + t334;
t177 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t379 + t332;
t171 = -t325 * t209 + t328 * t212;
t167 = pkin(4) * t171;
t166 = (t279 - t350) * pkin(3) + t331;
t165 = -pkin(3) * t350 + qJ(4) * t278 + t331;
t164 = (-t303 - t329) * qJ(4) + (-qJDD(3) - t296) * pkin(3) + t191 + t388;
t163 = -t394 - qJ(4) * t297 + (-t300 + t329) * pkin(3) - t334;
t160 = t327 * t175 - t383;
t159 = t324 * t175 + t233;
t158 = t327 * t176 + t324 * t229;
t157 = t324 * t176 - t327 * t229;
t156 = t326 * t186 - t323 * t188;
t155 = t323 * t186 + t326 * t188;
t154 = -pkin(3) * t245 + t326 * t165 - t393;
t153 = qJ(4) * t244 - t323 * t166 - t382;
t152 = -t323 * t178 + t326 * t179;
t151 = -t324 * t192 + t327 * t194 + t409;
t150 = -t324 * t191 + t327 * t193 - t384;
t149 = t327 * t192 + t324 * t194 - t407;
t148 = t327 * t191 + t324 * t193 + t348;
t147 = t327 * t152 - t324 * t274 - t383;
t146 = t324 * t152 + t327 * t274 + t233;
t145 = t327 * t156 + t324 * t177;
t144 = t324 * t156 - t327 * t177;
t143 = -t325 * t157 + t328 * t158;
t142 = t328 * t157 + t325 * t158;
t141 = -pkin(2) * t155 - pkin(3) * t188 - qJ(4) * t186;
t140 = t327 * t153 - t324 * t164 - t384;
t139 = t327 * t154 - t324 * t163 - t409;
t138 = -pkin(5) * t157 - (pkin(2) * t324 - pkin(6) * t327) * t175;
t137 = -pkin(6) * t155 + t341 * t177;
t136 = t324 * t153 + t327 * t164 + t348;
t135 = t324 * t154 + t327 * t163 + t407;
t134 = pkin(5) * t158 - (-pkin(2) * t327 - pkin(6) * t324 - pkin(1)) * t175;
t133 = -t325 * t144 + t328 * t145;
t132 = t328 * t144 + t325 * t145;
t131 = -pkin(5) * t144 + t327 * t137 - t324 * t141;
t130 = -pkin(1) * t155 + pkin(5) * t145 + t324 * t137 + t327 * t141;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t298, -t299, 0, t262, 0, 0, 0, 0, 0, 0, -t339, t239, 0, t162, 0, 0, 0, 0, 0, 0, t171, t173, t190, t143, 0, 0, 0, 0, 0, 0, t171, t190, -t173, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t299, -t298, 0, t261, 0, 0, 0, 0, 0, 0, -t239, -t339, 0, -t397, 0, 0, 0, 0, 0, 0, t168, -t169, t189, t142, 0, 0, 0, 0, 0, 0, t168, t189, t169, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t247, -t249, 0, -t175, 0, 0, 0, 0, 0, 0, t247, 0, t249, t155; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t299, 0, -t298, 0, t343, -t275, -t261, -pkin(4) * t261, 0, 0, -t239, 0, -t339, 0, t402, t401, t397, pkin(4) * t397 + pkin(5) * t359 - t325 * t195, t183, t405, t185, t182, t403, t205, -t325 * t148 + t328 * t150 - t386, -t325 * t149 + t328 * t151 + t411, -t325 * t159 + t328 * t160 - t385, -pkin(4) * t142 - t325 * t134 + t328 * t138, t183, t185, -t405, t205, -t403, t182, -t325 * t136 + t328 * t140 - t386, -t325 * t146 + t328 * t147 - t385, -t325 * t135 + t328 * t139 - t411, -pkin(4) * t132 - t325 * t130 + t328 * t131; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t298, 0, t299, 0, t275, t343, t262, pkin(4) * t262, 0, 0, t339, 0, -t239, 0, -t401, t402, t162, pkin(4) * t162 + pkin(5) * t366 + t328 * t195, t181, -t406, t184, t180, -t404, t204, t328 * t148 + t325 * t150 + t167, t328 * t149 + t325 * t151 + t410, t328 * t159 + t325 * t160 + t187, pkin(4) * t143 + t328 * t134 + t325 * t138, t181, t184, t406, t204, t404, t180, t328 * t136 + t325 * t140 + t167, t328 * t146 + t325 * t147 + t187, t328 * t135 + t325 * t139 - t410, pkin(4) * t133 + t328 * t130 + t325 * t131; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t306, t307, 0, 0, 0, 0, 0, 0, 0, t318, -pkin(1) * t291 - t242, -pkin(1) * t288 - t243, 0, -pkin(1) * t198, t245, t340, t248, t244, t389, 0, t352 - t365, -pkin(2) * t278 + t376 - t408, t176 + t351, pkin(1) * t157 - pkin(2) * t229 + pkin(6) * t176, t245, t248, -t340, 0, -t389, t244, qJ(4) * t375 + t326 * t166 + t352, t326 * t178 + t323 * t179 + t351, t323 * t165 + (pkin(2) + t387) * t278 + t408, pkin(1) * t144 + pkin(6) * t156 + (-pkin(2) + t342) * t177;];
tauB_reg = t1;
