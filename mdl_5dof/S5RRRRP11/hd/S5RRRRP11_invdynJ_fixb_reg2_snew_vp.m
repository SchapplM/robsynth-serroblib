% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:34
% EndTime: 2019-12-31 22:18:55
% DurationCPUTime: 9.26s
% Computational Cost: add. (24841->450), mult. (53626->611), div. (0->0), fcn. (41905->10), ass. (0->295)
t249 = cos(pkin(5));
t243 = qJD(1) * t249 + qJD(2);
t251 = sin(qJ(3));
t254 = cos(qJ(3));
t248 = sin(pkin(5));
t252 = sin(qJ(2));
t318 = qJD(1) * t252;
t302 = t248 * t318;
t220 = t243 * t251 + t254 * t302;
t255 = cos(qJ(2));
t317 = qJD(1) * t255;
t301 = t248 * t317;
t236 = -qJD(3) + t301;
t250 = sin(qJ(4));
t253 = cos(qJ(4));
t201 = t220 * t253 - t236 * t250;
t218 = -t243 * t254 + t251 * t302;
t214 = qJD(4) + t218;
t179 = t214 * t201;
t313 = qJDD(1) * t248;
t226 = qJD(2) * t301 + t252 * t313;
t242 = qJDD(1) * t249 + qJDD(2);
t284 = -t254 * t226 - t251 * t242;
t187 = -qJD(3) * t218 - t284;
t312 = qJDD(1) * t255;
t277 = qJD(2) * t318 - t312;
t271 = t277 * t248;
t267 = qJDD(3) + t271;
t299 = -t250 * t187 + t253 * t267;
t279 = qJD(4) * t201 - t299;
t113 = -t179 + t279;
t365 = t214 ^ 2;
t199 = t220 * t250 + t236 * t253;
t366 = t199 ^ 2;
t172 = t366 - t365;
t297 = t251 * t226 - t242 * t254;
t186 = -qJD(3) * t220 - t297;
t185 = qJDD(4) - t186;
t352 = t201 * t199;
t128 = -t352 - t185;
t339 = t250 * t128;
t96 = -t172 * t253 - t339;
t60 = -t254 * t113 + t251 * t96;
t328 = t253 * t128;
t92 = -t172 * t250 + t328;
t440 = t249 * t60 + (t252 * (t251 * t113 + t254 * t96) - t255 * t92) * t248;
t198 = t201 ^ 2;
t382 = t198 - t366;
t260 = -t253 * t187 - t250 * t267;
t258 = -t199 * qJD(4) - t260;
t353 = t199 * t214;
t378 = -t353 + t258;
t340 = t250 * t378;
t384 = t179 + t279;
t66 = t253 * t384 + t340;
t46 = t251 * t66 + t254 * t382;
t64 = -t250 * t384 + t253 * t378;
t438 = t248 * (t252 * (-t251 * t382 + t254 * t66) + t255 * t64) + t249 * t46;
t383 = -t198 - t365;
t86 = t253 * t383 + t339;
t437 = pkin(2) * t86;
t436 = pkin(3) * t86;
t435 = pkin(9) * t86;
t88 = -t250 * t383 + t328;
t434 = pkin(9) * t88;
t433 = t251 * t88;
t432 = t252 * t86;
t431 = t254 * t88;
t430 = t255 * t86;
t377 = t353 + t258;
t173 = -t198 + t365;
t379 = -t352 + t185;
t327 = t253 * t379;
t407 = -t173 * t250 + t327;
t338 = t250 * t379;
t408 = t173 * t253 + t338;
t420 = t251 * t407 - t254 * t377;
t425 = t249 * t420 + (t252 * (t251 * t377 + t254 * t407) - t255 * t408) * t248;
t376 = -t365 - t366;
t391 = t253 * t376 - t338;
t406 = t251 * t391 - t254 * t384;
t424 = pkin(2) * t406;
t423 = pkin(8) * t406;
t392 = t250 * t376 + t327;
t405 = t251 * t384 + t254 * t391;
t421 = -pkin(2) * t392 + pkin(8) * t405;
t419 = pkin(1) * (t252 * t405 - t255 * t392);
t418 = pkin(7) * (t252 * t392 + t255 * t405) - pkin(1) * t406;
t415 = pkin(3) * t392;
t414 = pkin(9) * t391;
t413 = pkin(9) * t392;
t381 = t198 + t366;
t404 = pkin(3) * t381;
t403 = qJ(5) * t378;
t400 = t251 * t381;
t395 = t254 * t381;
t351 = t214 * t250;
t107 = t199 * t351 - t253 * t279;
t350 = t214 * t253;
t308 = t199 * t350;
t280 = t250 * t279 + t308;
t309 = t251 * t352;
t307 = t254 * t352;
t369 = t251 * t280 + t307;
t390 = t249 * t369 + (-t255 * t107 + t252 * (t254 * t280 - t309)) * t248;
t170 = t201 * t351;
t289 = t170 - t308;
t371 = -t185 * t254 + t251 * t289;
t385 = (t199 * t250 + t201 * t253) * t214;
t389 = t249 * t371 + (t255 * t385 + t252 * (t185 * t251 + t254 * t289)) * t248;
t245 = t248 ^ 2;
t388 = t245 * t255;
t349 = t220 * t218;
t264 = t267 - t349;
t387 = t251 * t264;
t386 = t254 * t264;
t207 = t218 * t236;
t157 = t187 + t207;
t256 = qJD(1) ^ 2;
t361 = sin(qJ(1));
t362 = cos(qJ(1));
t276 = g(1) * t362 + g(2) * t361;
t222 = -pkin(1) * t256 + pkin(7) * t313 - t276;
t275 = g(1) * t361 - g(2) * t362;
t343 = t248 * t256;
t263 = qJDD(1) * pkin(1) + pkin(7) * t343 + t275;
t261 = t249 * t263;
t298 = t252 * t222 - t255 * t261;
t344 = t248 * t255;
t189 = g(3) * t344 + t298;
t345 = t248 * t252;
t259 = -g(3) * t345 + t252 * t261;
t190 = t255 * t222 + t259;
t380 = t252 * t189 + t255 * t190;
t162 = pkin(4) * t199 - qJ(5) * t201;
t294 = -pkin(2) * t255 - pkin(8) * t252;
t319 = qJD(1) * t248;
t225 = t294 * t319;
t241 = t243 ^ 2;
t151 = t242 * pkin(8) - t241 * pkin(2) + (t225 * t319 + t222) * t255 + t259;
t357 = t249 * g(3);
t257 = -t226 * pkin(8) - t357 + (-t243 * pkin(8) * t317 + (-t312 + (qJD(2) + t243) * t318) * pkin(2) - t263) * t248;
t100 = t151 * t254 + t251 * t257;
t191 = pkin(3) * t218 - pkin(9) * t220;
t364 = t236 ^ 2;
t74 = -pkin(3) * t364 + pkin(9) * t267 - t218 * t191 + t100;
t150 = -t242 * pkin(2) - t241 * pkin(8) + (g(3) * t255 + t225 * t318) * t248 + t298;
t77 = -t157 * pkin(9) + (-t220 * t236 - t186) * pkin(3) + t150;
t38 = t250 * t77 + t253 * t74;
t296 = qJ(5) * t185 - t162 * t199 + t38;
t375 = -(t383 + t365) * pkin(4) - qJ(5) * t128 + t296;
t316 = qJD(5) * t214;
t209 = 0.2e1 * t316;
t288 = t209 + t296;
t31 = -pkin(4) * t365 + t288;
t37 = t250 * t74 - t253 * t77;
t32 = -pkin(4) * t185 - qJ(5) * t365 + t162 * t201 + qJDD(5) + t37;
t15 = t250 * t32 + t253 * t31;
t300 = qJ(5) * t250 + pkin(3);
t99 = t251 * t151 - t254 * t257;
t73 = -pkin(3) * t267 - pkin(9) * t364 + t220 * t191 + t99;
t266 = pkin(4) * t279 - t403 + t73;
t34 = (pkin(4) * t214 - 0.2e1 * qJD(5)) * t201 + t266;
t359 = pkin(4) * t253;
t374 = -(t300 + t359) * t34 + pkin(9) * t15;
t265 = 0.2e1 * qJD(5) * t201 - t266;
t28 = (-t384 - t179) * pkin(4) + t265;
t373 = t253 * t28 - t300 * t384 + t414;
t27 = -pkin(4) * t179 + t265 + t403;
t372 = -t434 + t378 * (pkin(3) + t359) + t250 * t27;
t154 = (qJD(3) + t236) * t220 + t297;
t110 = t201 * t350 + t250 * t258;
t111 = t253 * t258 - t170;
t290 = t111 * t251 - t307;
t370 = t249 * t290 + (t111 * t254 + t309) * t345 - t110 * t344;
t216 = t218 ^ 2;
t217 = t220 ^ 2;
t360 = pkin(3) * t251;
t358 = t248 * pkin(7);
t356 = t250 * t73;
t355 = t253 * t73;
t354 = qJ(5) * t253;
t348 = t236 * t251;
t347 = t236 * t254;
t346 = t245 * t256;
t341 = t250 * t377;
t335 = t251 * t150;
t175 = -t267 - t349;
t334 = t251 * t175;
t235 = t255 * t252 * t346;
t223 = t235 + t242;
t331 = t252 * t223;
t329 = t253 * t377;
t324 = t254 * t150;
t323 = t254 * t175;
t224 = -t235 + t242;
t321 = t255 * t224;
t315 = qJD(3) - t236;
t246 = t252 ^ 2;
t311 = t246 * t346;
t247 = t255 ^ 2;
t310 = t247 * t346;
t306 = t255 * t349;
t230 = t243 * t301;
t305 = t230 + t226;
t304 = -pkin(3) * t254 - pkin(2);
t24 = t250 * t37 + t253 * t38;
t57 = t100 * t254 + t251 * t99;
t295 = -pkin(3) * t73 + pkin(9) * t24;
t293 = -pkin(4) * t32 + qJ(5) * t31;
t287 = -pkin(4) * t377 - qJ(5) * t113;
t23 = t250 * t38 - t253 * t37;
t286 = t100 * t251 - t254 * t99;
t283 = qJD(1) * t243 - t249 * t256;
t282 = -pkin(1) + t294;
t274 = -pkin(3) * t384 - t355 + t414;
t120 = (qJD(4) + t214) * t199 + t260;
t273 = pkin(3) * t120 + t356 + t434;
t115 = (-qJD(4) + t214) * t201 + t299;
t69 = t115 * t253 + t341;
t272 = pkin(9) * t69 + t24 + t404;
t25 = (t381 - t365) * pkin(4) + t288;
t26 = qJ(5) * t381 + t32;
t67 = -t113 * t253 + t341;
t269 = pkin(9) * t67 + t25 * t253 + t250 * t26 + t404;
t262 = pkin(4) * t379 + qJ(5) * t376 - t32;
t229 = t243 * t302;
t228 = (t246 - t247) * t346;
t227 = -t241 - t310;
t213 = -t311 - t241;
t208 = t248 * t263 + t357;
t206 = -t229 - t271;
t205 = t229 - t271;
t204 = -t230 + t226;
t203 = -t217 + t364;
t202 = t216 - t364;
t194 = -t217 - t364;
t193 = t217 - t216;
t188 = -t364 - t216;
t169 = t216 + t217;
t168 = (t218 * t251 + t220 * t254) * t236;
t159 = t218 * t315 + t284;
t158 = t187 - t207;
t155 = -t220 * t315 - t297;
t153 = t187 * t251 - t220 * t347;
t152 = t186 * t254 - t218 * t348;
t144 = t202 * t251 - t323;
t143 = t203 * t254 + t387;
t142 = -t194 * t251 + t323;
t141 = t194 * t254 + t334;
t135 = t188 * t254 - t387;
t134 = t188 * t251 + t386;
t123 = -t154 * t254 + t158 * t251;
t121 = t155 * t251 + t157 * t254;
t72 = pkin(2) * t159 + pkin(8) * t142 + t335;
t70 = pkin(2) * t155 + pkin(8) * t135 - t324;
t65 = t115 * t250 - t329;
t63 = -t113 * t250 - t329;
t54 = -t120 * t251 + t431;
t52 = t120 * t254 + t433;
t50 = -t251 * t378 - t431;
t48 = t254 * t378 - t433;
t45 = t254 * t69 - t400;
t44 = t254 * t67 - t400;
t43 = t251 * t69 + t395;
t42 = t251 * t67 + t395;
t41 = -pkin(2) * t150 + pkin(8) * t57;
t40 = t355 - t435;
t39 = t356 - t413;
t35 = pkin(2) * t169 + pkin(8) * t123 + t57;
t33 = -pkin(3) * t63 - t287;
t30 = t38 - t436;
t29 = t37 - t415;
t22 = -t262 - t415;
t21 = -0.2e1 * t316 - t375 + t436;
t20 = -t250 * t28 - t354 * t384 - t413;
t19 = -pkin(4) * t340 + t253 * t27 + t435;
t18 = t24 * t254 + t251 * t73;
t17 = t24 * t251 - t254 * t73;
t16 = -pkin(9) * t65 - t23;
t14 = t250 * t31 - t253 * t32;
t13 = pkin(8) * t54 + t251 * t40 + t254 * t30 - t437;
t12 = t251 * t39 + t254 * t29 + t421;
t11 = -pkin(9) * t63 - t25 * t250 + t253 * t26;
t10 = t15 * t254 + t251 * t34;
t9 = t15 * t251 - t254 * t34;
t8 = pkin(8) * t45 + t251 * t16 + t304 * t65;
t7 = -pkin(9) * t14 + (pkin(4) * t250 - t354) * t34;
t6 = -pkin(3) * t14 - t293;
t5 = t20 * t251 + t22 * t254 + t421;
t4 = pkin(8) * t50 + t19 * t251 + t21 * t254 + t437;
t3 = -pkin(2) * t63 + pkin(8) * t44 + t11 * t251 + t254 * t33;
t2 = pkin(8) * t18 + (-pkin(9) * t251 + t304) * t23;
t1 = -pkin(2) * t14 + pkin(8) * t10 + t251 * t7 + t254 * t6;
t36 = [0, 0, 0, 0, 0, qJDD(1), t275, t276, 0, 0, (t226 * t248 + t283 * t388) * t252, t249 * t228 + (t252 * t206 + t255 * t305) * t248, t249 * t204 + (t331 + t255 * (t241 - t311)) * t248, (-t252 * t283 - t277) * t388, t249 * t205 + (t252 * (-t241 + t310) + t321) * t248, t249 * t242, (-t189 + pkin(1) * (t223 * t255 + t227 * t252)) * t249 + (t255 * t208 + pkin(1) * t206 + pkin(7) * (t227 * t255 - t331)) * t248, -t208 * t345 - t249 * t190 + pkin(1) * (-t248 * t305 + (t213 * t255 - t224 * t252) * t249) + (-t213 * t252 - t321) * t358, pkin(1) * ((-t204 * t255 + t205 * t252) * t249 - (-t246 - t247) * t245 * t343) + (t204 * t252 + t205 * t255) * t358 + t380 * t248, pkin(1) * (t248 * t208 + (-t189 * t255 + t190 * t252) * t249) + t380 * t358, t249 * t153 + (t252 * (t187 * t254 + t220 * t348) - t306) * t248, t249 * t121 + (t252 * (t155 * t254 - t157 * t251) - t255 * t193) * t248, t249 * t143 + (t252 * (-t203 * t251 + t386) - t255 * t158) * t248, t249 * t152 + (t252 * (-t186 * t251 - t218 * t347) + t306) * t248, t249 * t144 + (t252 * (t202 * t254 + t334) + t255 * t154) * t248, -t267 * t344 + t249 * t168 + (t218 * t254 - t220 * t251) * t236 * t345, (t70 + pkin(1) * (t135 * t252 + t155 * t255)) * t249 + (t252 * (-pkin(8) * t134 + t335) + t255 * (-pkin(2) * t134 + t99) - pkin(1) * t134 + pkin(7) * (t135 * t255 - t155 * t252)) * t248, (t72 + pkin(1) * (t142 * t252 + t159 * t255)) * t249 + (t252 * (-pkin(8) * t141 + t324) + t255 * (-pkin(2) * t141 + t100) - pkin(1) * t141 + pkin(7) * (t142 * t255 - t159 * t252)) * t248, (t35 + pkin(1) * (t123 * t252 + t169 * t255)) * t249 + (-t252 * t286 + pkin(7) * (t123 * t255 - t169 * t252) + t282 * (-t154 * t251 - t158 * t254)) * t248, (t41 + pkin(1) * (-t150 * t255 + t252 * t57)) * t249 + (pkin(7) * (t150 * t252 + t255 * t57) + t282 * t286) * t248, t370, -t438, t425, t390, -t440, t389, (t12 + t419) * t249 + (t252 * (-t251 * t29 + t254 * t39 - t423) + t255 * (-t274 - t424) + t418) * t248, (t13 + pkin(1) * (t252 * t54 - t430)) * t249 + (t252 * (-pkin(8) * t52 - t251 * t30 + t254 * t40) + t255 * (-pkin(2) * t52 - t273) - pkin(1) * t52 + pkin(7) * (t255 * t54 + t432)) * t248, (t8 + pkin(1) * (t252 * t45 - t255 * t65)) * t249 + (t252 * (-pkin(8) * t43 + t16 * t254 + t360 * t65) + t255 * (-pkin(2) * t43 - t272) - pkin(1) * t43 + pkin(7) * (t252 * t65 + t255 * t45)) * t248, (t2 + pkin(1) * (t18 * t252 - t23 * t255)) * t249 + (t252 * (-pkin(8) * t17 + (-pkin(9) * t254 + t360) * t23) + t255 * (-pkin(2) * t17 - t295) - pkin(1) * t17 + pkin(7) * (t18 * t255 + t23 * t252)) * t248, t370, t425, t438, t389, t440, t390, (t5 + t419) * t249 + (t252 * (t20 * t254 - t22 * t251 - t423) + t255 * (-t373 - t424) + t418) * t248, (t3 + pkin(1) * (t252 * t44 - t255 * t63)) * t249 + (t252 * (-pkin(8) * t42 + t11 * t254 - t251 * t33) + t255 * (-pkin(2) * t42 - t269) - pkin(1) * t42 + pkin(7) * (t252 * t63 + t255 * t44)) * t248, (t4 + pkin(1) * (t252 * t50 + t430)) * t249 + (t252 * (-pkin(8) * t48 + t19 * t254 - t21 * t251) + t255 * (-pkin(2) * t48 - t372) - pkin(1) * t48 + pkin(7) * (t255 * t50 - t432)) * t248, (t1 + pkin(1) * (t10 * t252 - t14 * t255)) * t249 + (-pkin(1) * t9 + (pkin(7) * t14 - pkin(8) * t9 - t251 * t6 + t254 * t7) * t252 + (-pkin(2) * t9 + pkin(7) * t10 - t374) * t255) * t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, t228, t204, t235, t205, t242, -t189, -t190, 0, 0, t153, t121, t143, t152, t144, t168, t70, t72, t35, t41, t290, -t46, t420, t369, -t60, t371, t12, t13, t8, t2, t290, t420, t46, t371, t60, t369, t5, t3, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t349, t193, t158, -t349, -t154, t267, -t99, -t100, 0, 0, t110, t64, t408, t107, -t92, -t385, t274, t273, t272, t295, t110, t408, -t64, -t385, t92, t107, t373, t269, t372, t374; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t352, t382, t377, -t352, -t113, t185, -t37, -t38, 0, 0, t352, t377, -t382, t185, t113, -t352, t262, t287, t209 + t375, t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t379, t377, t383, t32;];
tauJ_reg = t36;
