% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRRR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:53
% EndTime: 2019-12-05 15:13:01
% DurationCPUTime: 5.41s
% Computational Cost: add. (16769->359), mult. (24356->545), div. (0->0), fcn. (19030->10), ass. (0->242)
t344 = qJD(3) + qJD(4);
t342 = t344 ^ 2;
t356 = cos(qJ(4));
t343 = qJDD(3) + qJDD(4);
t353 = sin(qJ(4));
t384 = t353 * t343;
t310 = t356 * t342 + t384;
t378 = t356 * t343;
t313 = t353 * t342 - t378;
t354 = sin(qJ(3));
t357 = cos(qJ(3));
t259 = t357 * t310 - t354 * t313;
t264 = t354 * t310 + t357 * t313;
t348 = sin(pkin(9));
t350 = cos(pkin(9));
t224 = t348 * t259 + t350 * t264;
t349 = sin(pkin(8));
t425 = t349 * t224;
t351 = cos(pkin(8));
t424 = t351 * t224;
t327 = t349 * g(1) - t351 * g(2);
t324 = -qJDD(2) + t327;
t271 = pkin(6) * t310 - t356 * t324;
t413 = pkin(6) * t313 - t353 * t324;
t202 = pkin(5) * t259 + t357 * t271 - t354 * t413;
t419 = pkin(5) * t264 + t354 * t271 + t357 * t413;
t147 = qJ(2) * t224 + t348 * t202 + t350 * t419;
t411 = t350 * t259 - t348 * t264;
t148 = qJ(2) * t411 + t350 * t202 - t348 * t419;
t328 = t351 * g(1) + t349 * g(2);
t347 = g(3) - qJDD(1);
t300 = -t348 * t328 + t350 * t347;
t301 = -t350 * t328 - t348 * t347;
t247 = -t354 * t300 + t357 * t301;
t400 = qJD(3) ^ 2;
t238 = -t400 * pkin(3) + t247;
t246 = t357 * t300 + t354 * t301;
t359 = qJDD(3) * pkin(3) - t246;
t193 = t353 * t238 - t356 * t359;
t194 = t356 * t238 + t353 * t359;
t367 = t353 * t193 + t356 * t194;
t164 = t356 * t193 - t353 * t194;
t377 = t357 * t164;
t135 = -t354 * t367 + t377;
t383 = t354 * t164;
t414 = t357 * t367 + t383;
t123 = t348 * t135 + t350 * t414;
t122 = t350 * t135 - t348 * t414;
t325 = t354 * qJDD(3) + t357 * t400;
t326 = t357 * qJDD(3) - t354 * t400;
t364 = -t348 * t325 + t350 * t326;
t418 = t349 * t364;
t417 = t351 * t364;
t366 = t354 * t246 + t357 * t247;
t205 = t357 * t246 - t354 * t247;
t392 = t350 * t205;
t171 = -t348 * t366 + t392;
t397 = t348 * t205;
t172 = t350 * t366 + t397;
t294 = pkin(5) * t325 - t357 * t324;
t360 = -pkin(5) * t326 - t354 * t324;
t211 = -qJ(2) * t364 + t348 * t294 + t350 * t360;
t402 = t350 * t325 + t348 * t326;
t212 = qJ(2) * t402 + t350 * t294 - t348 * t360;
t399 = pkin(2) * t324;
t352 = sin(qJ(5));
t345 = t352 ^ 2;
t398 = t345 * t342;
t394 = t349 * t324;
t393 = t349 * t347;
t391 = t350 * t324;
t309 = t351 * t324;
t390 = t351 * t347;
t189 = -t343 * pkin(4) - t342 * pkin(7) + t193;
t389 = t352 * t189;
t355 = cos(qJ(5));
t333 = t352 * t342 * t355;
t322 = qJDD(5) + t333;
t388 = t352 * t322;
t323 = qJDD(5) - t333;
t387 = t352 * t323;
t386 = t352 * t343;
t382 = t355 * t189;
t381 = t355 * t322;
t380 = t355 * t323;
t335 = t355 * t343;
t190 = -t342 * pkin(4) + t343 * pkin(7) + t194;
t184 = t355 * t190 - t352 * t324;
t346 = t355 ^ 2;
t376 = t345 + t346;
t375 = qJD(5) * t344;
t374 = t352 * t375;
t373 = t355 * t375;
t183 = t352 * t190 + t355 * t324;
t151 = t352 * t183 + t355 * t184;
t308 = t376 * t343;
t337 = t346 * t342;
t316 = t337 + t398;
t258 = t353 * t308 + t356 * t316;
t262 = t356 * t308 - t353 * t316;
t226 = t357 * t258 + t354 * t262;
t227 = -t354 * t258 + t357 * t262;
t179 = t350 * t226 + t348 * t227;
t180 = -t348 * t226 + t350 * t227;
t372 = -pkin(1) * t179 - pkin(2) * t226 - pkin(3) * t258 - pkin(4) * t316 - pkin(7) * t308 + qJ(1) * t180 - t151;
t371 = pkin(1) * t411 + pkin(2) * t259 + pkin(3) * t310 + qJ(1) * t224 + t194;
t370 = pkin(1) * t224 + pkin(2) * t264 + pkin(3) * t313 - qJ(1) * t411 + t193;
t369 = pkin(1) * t402 + pkin(2) * t325 - qJ(1) * t364 + t247;
t368 = -pkin(1) * t364 - pkin(2) * t326 - qJ(1) * t402 + t246;
t245 = t348 * t300 + t350 * t301;
t290 = -t349 * t327 - t351 * t328;
t362 = t353 * t333;
t361 = t356 * t333;
t150 = t355 * t183 - t352 * t184;
t244 = t350 * t300 - t348 * t301;
t289 = t351 * t327 - t349 * t328;
t358 = qJD(5) ^ 2;
t332 = -t337 - t358;
t331 = t337 - t358;
t330 = -t358 - t398;
t329 = t358 - t398;
t321 = pkin(1) * t324;
t317 = t337 - t398;
t307 = t335 - 0.2e1 * t374;
t306 = t335 - t374;
t305 = t373 + t386;
t304 = 0.2e1 * t373 + t386;
t303 = t376 * t375;
t296 = t353 * qJDD(5) + t356 * t303;
t295 = -t356 * qJDD(5) + t353 * t303;
t288 = t355 * t305 - t345 * t375;
t287 = -t352 * t306 - t346 * t375;
t286 = -t352 * t330 - t380;
t285 = -t352 * t329 + t381;
t284 = t355 * t332 - t388;
t283 = t355 * t331 - t387;
t279 = t355 * t330 - t387;
t278 = -t355 * t329 - t388;
t277 = t352 * t332 + t381;
t276 = -t352 * t331 - t380;
t273 = (-t305 - t373) * t352;
t272 = (-t306 + t374) * t355;
t267 = t351 * t402;
t266 = t349 * t402;
t257 = -t352 * t304 + t355 * t307;
t256 = -t355 * t304 - t352 * t307;
t255 = t356 * t285 + t352 * t384;
t254 = t356 * t283 + t353 * t335;
t253 = t353 * t285 - t352 * t378;
t252 = t353 * t283 - t355 * t378;
t251 = t356 * t288 - t362;
t250 = t356 * t287 + t362;
t249 = t353 * t288 + t361;
t248 = t353 * t287 - t361;
t242 = t356 * t286 + t353 * t304;
t241 = t356 * t284 - t353 * t307;
t240 = t353 * t286 - t356 * t304;
t239 = t353 * t284 + t356 * t307;
t233 = -t354 * t295 + t357 * t296;
t232 = t357 * t295 + t354 * t296;
t231 = t356 * t257 - t353 * t317;
t230 = t353 * t257 + t356 * t317;
t229 = t351 * t245 - t394;
t228 = t349 * t245 + t309;
t220 = t351 * t411;
t219 = t349 * t411;
t216 = -t354 * t253 + t357 * t255;
t215 = -t354 * t252 + t357 * t254;
t214 = t357 * t253 + t354 * t255;
t213 = t357 * t252 + t354 * t254;
t210 = -t354 * t249 + t357 * t251;
t209 = -t354 * t248 + t357 * t250;
t208 = t357 * t249 + t354 * t251;
t207 = t357 * t248 + t354 * t250;
t198 = -t354 * t240 + t357 * t242;
t197 = -t354 * t239 + t357 * t241;
t196 = t357 * t240 + t354 * t242;
t195 = t357 * t239 + t354 * t241;
t191 = pkin(5) * t366 + t399;
t187 = -t348 * t232 + t350 * t233;
t186 = -t354 * t230 + t357 * t231;
t185 = t357 * t230 + t354 * t231;
t182 = -pkin(7) * t279 + t382;
t181 = -pkin(7) * t277 + t389;
t178 = -pkin(4) * t279 + t184;
t177 = -pkin(4) * t277 + t183;
t176 = -t348 * t214 + t350 * t216;
t175 = -t348 * t213 + t350 * t215;
t174 = -t348 * t208 + t350 * t210;
t173 = -t348 * t207 + t350 * t209;
t169 = -t348 * t196 + t350 * t198;
t168 = -t348 * t195 + t350 * t197;
t167 = t350 * t196 + t348 * t198;
t166 = t350 * t195 + t348 * t197;
t161 = t351 * t172 - t394;
t160 = t349 * t172 + t309;
t159 = pkin(3) * t324 + pkin(6) * t367;
t156 = -t348 * t185 + t350 * t186;
t155 = t351 * t169 + t349 * t279;
t154 = t351 * t168 + t349 * t277;
t153 = t349 * t169 - t351 * t279;
t152 = t349 * t168 - t351 * t277;
t146 = pkin(1) * t171 + pkin(2) * t205;
t145 = -pkin(6) * t258 + t356 * t150;
t144 = pkin(6) * t262 + t353 * t150;
t143 = -pkin(6) * t240 - t353 * t178 + t356 * t182;
t142 = -pkin(6) * t239 - t353 * t177 + t356 * t181;
t141 = -pkin(3) * t279 + pkin(6) * t242 + t356 * t178 + t353 * t182;
t140 = -pkin(3) * t277 + pkin(6) * t241 + t356 * t177 + t353 * t181;
t139 = t356 * t151 + t353 * t189;
t138 = t353 * t151 - t356 * t189;
t137 = pkin(5) * t392 + qJ(2) * t171 - t348 * t191;
t132 = -pkin(1) * t167 - pkin(2) * t196 - pkin(3) * t240 + pkin(4) * t304 - pkin(7) * t286 - t389;
t131 = -pkin(1) * t166 - pkin(2) * t195 - pkin(3) * t239 - pkin(4) * t307 - pkin(7) * t284 + t382;
t129 = -pkin(5) * t226 - t354 * t144 + t357 * t145;
t128 = pkin(5) * t227 + t357 * t144 + t354 * t145;
t127 = -t354 * t138 + t357 * t139;
t126 = t357 * t138 + t354 * t139;
t125 = -pkin(5) * t196 - t354 * t141 + t357 * t143;
t124 = -pkin(5) * t195 - t354 * t140 + t357 * t142;
t120 = pkin(5) * t135 + pkin(6) * t377 - t354 * t159;
t119 = -pkin(2) * t279 + pkin(5) * t198 + t357 * t141 + t354 * t143;
t118 = -pkin(2) * t277 + pkin(5) * t197 + t357 * t140 + t354 * t142;
t117 = t351 * t123 - t394;
t116 = t349 * t123 + t309;
t115 = pkin(5) * t414 + pkin(6) * t383 + t357 * t159 + t399;
t114 = -pkin(6) * t138 - (pkin(4) * t353 - pkin(7) * t356) * t150;
t113 = pkin(6) * t139 - (-pkin(4) * t356 - pkin(7) * t353 - pkin(3)) * t150;
t112 = pkin(1) * t122 + pkin(2) * t135 + pkin(3) * t164;
t111 = -qJ(2) * t179 - t348 * t128 + t350 * t129;
t110 = -t348 * t126 + t350 * t127;
t109 = t350 * t126 + t348 * t127;
t108 = t351 * t110 - t150 * t349;
t107 = t349 * t110 + t150 * t351;
t106 = -qJ(2) * t167 - t348 * t119 + t350 * t125;
t105 = -qJ(2) * t166 - t348 * t118 + t350 * t124;
t104 = qJ(2) * t122 - t348 * t115 + t350 * t120;
t103 = -pkin(5) * t126 - t354 * t113 + t357 * t114;
t102 = -pkin(1) * t109 - pkin(2) * t126 - pkin(3) * t138 + pkin(4) * t189 - pkin(7) * t151;
t101 = pkin(2) * t150 + pkin(5) * t127 + t357 * t113 + t354 * t114;
t100 = -qJ(2) * t109 - t348 * t101 + t350 * t103;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, 0, 0, 0, 0, 0, 0, -t267, -t417, 0, t161, 0, 0, 0, 0, 0, 0, -t220, t424, 0, t117, 0, 0, 0, 0, 0, 0, t154, t155, t351 * t180, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, 0, 0, 0, 0, 0, 0, -t266, -t418, 0, t160, 0, 0, 0, 0, 0, 0, -t219, t425, 0, t116, 0, 0, 0, 0, 0, 0, t152, t153, t349 * t180, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, 0, 0, 0, 0, 0, 0, t364, -t402, 0, -t171, 0, 0, 0, 0, 0, 0, -t224, -t411, 0, -t122, 0, 0, 0, 0, 0, 0, t166, t167, t179, t109; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t393, -t390, -t289, -qJ(1) * t289, 0, 0, 0, 0, 0, 0, -t349 * t300 - t348 * t309, -t349 * t301 - t350 * t309, t351 * t244, -qJ(1) * t228 - (pkin(1) * t349 - qJ(2) * t351) * t244, 0, 0, t417, 0, -t267, t349 * qJDD(3), t351 * t211 - t349 * t368, t351 * t212 - t349 * t369, t351 * t171, -qJ(1) * t160 + t351 * t137 - t349 * t146, 0, 0, -t424, 0, -t220, t349 * t343, t351 * t147 - t349 * t370, t351 * t148 - t349 * t371, t351 * t122, -qJ(1) * t116 + t351 * t104 - t349 * t112, t351 * t174 - t349 * t273, t351 * t156 - t349 * t256, t351 * t176 - t349 * t278, t351 * t173 - t349 * t272, t351 * t175 - t349 * t276, t351 * t187, -qJ(1) * t152 + t351 * t105 - t349 * t131, -qJ(1) * t153 + t351 * t106 - t349 * t132, t351 * t111 - t349 * t372, -qJ(1) * t107 + t351 * t100 - t349 * t102; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t390, -t393, t290, qJ(1) * t290, 0, 0, 0, 0, 0, 0, t351 * t300 - t348 * t394, t351 * t301 - t349 * t391, t349 * t244, qJ(1) * t229 - (-pkin(1) * t351 - qJ(2) * t349) * t244, 0, 0, t418, 0, -t266, -t351 * qJDD(3), t349 * t211 + t351 * t368, t349 * t212 + t351 * t369, t349 * t171, qJ(1) * t161 + t349 * t137 + t351 * t146, 0, 0, -t425, 0, -t219, -t351 * t343, t349 * t147 + t351 * t370, t349 * t148 + t351 * t371, t349 * t122, qJ(1) * t117 + t349 * t104 + t351 * t112, t349 * t174 + t351 * t273, t349 * t156 + t351 * t256, t349 * t176 + t351 * t278, t349 * t173 + t351 * t272, t349 * t175 + t351 * t276, t349 * t187, qJ(1) * t154 + t349 * t105 + t351 * t131, qJ(1) * t155 + t349 * t106 + t351 * t132, t349 * t111 + t351 * t372, qJ(1) * t108 + t349 * t100 + t351 * t102; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t327, t328, 0, 0, 0, 0, 0, 0, 0, 0, t391, -t348 * t324, t245, qJ(2) * t245 + t321, 0, 0, t402, 0, t364, 0, -t212, t211, t172, pkin(5) * t397 + qJ(2) * t172 + t350 * t191 + t321, 0, 0, t411, 0, -t224, 0, -t148, t147, t123, qJ(2) * t123 + t350 * t115 + t348 * t120 + t321, t350 * t208 + t348 * t210, t350 * t185 + t348 * t186, t350 * t214 + t348 * t216, t350 * t207 + t348 * t209, t350 * t213 + t348 * t215, t350 * t232 + t348 * t233, -pkin(1) * t277 + qJ(2) * t168 + t350 * t118 + t348 * t124, -pkin(1) * t279 + qJ(2) * t169 + t350 * t119 + t348 * t125, qJ(2) * t180 + t350 * t128 + t348 * t129, pkin(1) * t150 + qJ(2) * t110 + t350 * t101 + t348 * t103;];
tauB_reg = t1;
