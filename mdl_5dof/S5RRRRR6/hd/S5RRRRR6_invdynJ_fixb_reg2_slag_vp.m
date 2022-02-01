% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:23
% EndTime: 2022-01-20 12:08:35
% DurationCPUTime: 4.06s
% Computational Cost: add. (8658->432), mult. (13431->552), div. (0->0), fcn. (9125->16), ass. (0->259)
t256 = cos(qJ(4));
t257 = cos(qJ(3));
t340 = t256 * t257;
t251 = sin(qJ(4));
t252 = sin(qJ(3));
t344 = t251 * t252;
t169 = -t340 + t344;
t260 = -pkin(8) - pkin(7);
t316 = qJD(3) * t260;
t175 = t252 * t316;
t176 = t257 * t316;
t203 = t260 * t252;
t238 = t257 * pkin(8);
t204 = pkin(7) * t257 + t238;
t258 = cos(qJ(2));
t332 = qJD(1) * t258;
t321 = pkin(1) * t332;
t326 = qJD(4) * t256;
t327 = qJD(4) * t251;
t358 = t169 * t321 + t256 * t175 + t251 * t176 + t203 * t326 - t204 * t327;
t127 = t251 * t203 + t256 * t204;
t341 = t252 * t256;
t170 = t251 * t257 + t341;
t357 = -qJD(4) * t127 + t170 * t321 - t175 * t251 + t256 * t176;
t244 = qJD(1) + qJD(2);
t319 = t244 * t344;
t143 = -t244 * t340 + t319;
t145 = t170 * t244;
t250 = sin(qJ(5));
t255 = cos(qJ(5));
t283 = t143 * t250 - t255 * t145;
t92 = t255 * t143 + t145 * t250;
t390 = t92 * t283;
t243 = qJD(3) + qJD(4);
t285 = t243 * t344;
t328 = qJD(3) * t257;
t385 = -t256 * t328 - t257 * t326;
t113 = t285 + t385;
t376 = pkin(9) * t113;
t389 = t376 + t357;
t114 = t243 * t170;
t111 = t114 * pkin(9);
t388 = -t111 + t358;
t246 = t252 ^ 2;
t247 = t257 ^ 2;
t333 = t246 + t247;
t383 = t333 * t258;
t387 = t244 * t383;
t242 = qJDD(1) + qJDD(2);
t253 = sin(qJ(2));
t323 = qJDD(1) * t253;
t330 = qJD(2) * t258;
t148 = pkin(7) * t242 + (qJD(1) * t330 + t323) * pkin(1);
t302 = t333 * t148;
t331 = qJD(2) * t253;
t230 = pkin(1) * t331;
t380 = pkin(1) * t258;
t335 = -qJD(1) * t230 + qJDD(1) * t380;
t379 = pkin(2) * t242;
t147 = -t335 - t379;
t249 = qJ(1) + qJ(2);
t236 = cos(t249);
t372 = g(2) * t236;
t386 = t147 + t372;
t22 = t283 ^ 2 - t92 ^ 2;
t324 = qJD(5) * t255;
t325 = qJD(5) * t250;
t339 = t257 * t242;
t293 = -t242 * t341 + t385 * t244 - t251 * t339;
t67 = t244 * t285 + t293;
t342 = t252 * t242;
t287 = t251 * t342 - t256 * t339;
t68 = t114 * t244 + t287;
t23 = t143 * t324 + t145 * t325 + t250 * t68 + t255 * t67;
t232 = qJD(5) + t243;
t18 = t232 * t92 - t23;
t248 = qJ(3) + qJ(4);
t237 = qJ(5) + t248;
t222 = sin(t237);
t223 = cos(t237);
t351 = t223 * t236;
t234 = sin(t249);
t352 = t223 * t234;
t241 = qJDD(3) + qJDD(4);
t363 = pkin(1) * qJD(1);
t322 = t253 * t363;
t184 = pkin(7) * t244 + t322;
t313 = t244 * t328;
t73 = -t184 * t328 + qJDD(3) * pkin(3) - t148 * t252 + (-t313 - t342) * pkin(8);
t309 = pkin(8) * t244 + t184;
t133 = t309 * t257;
t124 = t256 * t133;
t132 = t309 * t252;
t125 = qJD(3) * pkin(3) - t132;
t77 = t125 * t251 + t124;
t329 = qJD(3) * t252;
t314 = t244 * t329;
t78 = -t184 * t329 + t148 * t257 + (-t314 + t339) * pkin(8);
t26 = -qJD(4) * t77 - t251 * t78 + t256 * t73;
t16 = pkin(4) * t241 + pkin(9) * t67 + t26;
t296 = -t125 * t326 + t133 * t327 - t251 * t73 - t256 * t78;
t17 = -pkin(9) * t68 - t296;
t141 = t145 * pkin(9);
t122 = t251 * t133;
t76 = t256 * t125 - t122;
t60 = -t141 + t76;
t58 = pkin(4) * t243 + t60;
t375 = pkin(9) * t143;
t61 = t77 - t375;
t4 = t255 * (qJD(5) * t58 + t17) + t16 * t250 - t61 * t325;
t377 = pkin(3) * t257;
t226 = pkin(2) + t377;
t146 = -t226 * t244 - t321;
t97 = pkin(4) * t143 + t146;
t266 = g(1) * t351 + g(2) * t352 + g(3) * t222 + t92 * t97 - t4;
t353 = t222 * t236;
t354 = t222 * t234;
t361 = t255 * t61;
t28 = t250 * t58 + t361;
t5 = -qJD(5) * t28 + t255 * t16 - t17 * t250;
t264 = g(1) * t353 + g(2) * t354 - g(3) * t223 + t283 * t97 + t5;
t378 = pkin(2) * t244;
t185 = -t321 - t378;
t384 = t184 * t383 + t185 * t253;
t269 = qJD(5) * t283 + t250 * t67 - t255 * t68;
t19 = -t232 * t283 + t269;
t224 = pkin(1) * t253 + pkin(7);
t366 = -pkin(8) - t224;
t166 = t366 * t252;
t167 = t224 * t257 + t238;
t105 = t251 * t166 + t256 * t167;
t382 = g(1) * t236 + g(2) * t234;
t254 = sin(qJ(1));
t381 = pkin(1) * t254;
t374 = pkin(9) * t170;
t218 = g(1) * t234;
t373 = g(1) * t254;
t235 = cos(t248);
t371 = g(3) * t235;
t370 = g(3) * t257;
t126 = t256 * t203 - t204 * t251;
t98 = t126 - t374;
t165 = t169 * pkin(9);
t99 = -t165 + t127;
t56 = -t250 * t99 + t255 * t98;
t365 = qJD(5) * t56 + t389 * t250 + t388 * t255;
t57 = t250 * t98 + t255 * t99;
t364 = -qJD(5) * t57 - t388 * t250 + t389 * t255;
t362 = t250 * t61;
t225 = pkin(3) * t256 + pkin(4);
t345 = t250 * t251;
t79 = t132 * t251 - t124;
t62 = t79 + t375;
t80 = -t256 * t132 - t122;
t63 = -t141 + t80;
t360 = -t250 * t62 - t255 * t63 + t225 * t324 + (-t251 * t325 + (t255 * t256 - t345) * qJD(4)) * pkin(3);
t343 = t251 * t255;
t359 = t250 * t63 - t255 * t62 - t225 * t325 + (-t251 * t324 + (-t250 * t256 - t343) * qJD(4)) * pkin(3);
t356 = t145 * t143;
t233 = sin(t248);
t350 = t233 * t234;
t349 = t233 * t236;
t348 = t234 * t235;
t347 = t235 * t236;
t346 = t244 * t252;
t338 = t185 * t329 + t257 * t218;
t337 = t236 * pkin(2) + t234 * pkin(7);
t334 = t246 - t247;
t320 = pkin(1) * t330;
t229 = pkin(3) * t329;
t240 = t244 ^ 2;
t318 = t252 * t240 * t257;
t317 = t185 * t328 + t386 * t252;
t315 = t244 * t331;
t310 = pkin(4) * t235 + t377;
t100 = pkin(4) * t114 + t229;
t305 = qJD(3) * t366;
t301 = t333 * t242;
t104 = t256 * t166 - t167 * t251;
t174 = pkin(2) + t310;
t245 = pkin(9) - t260;
t299 = t236 * t174 + t234 * t245;
t298 = -t174 * t234 + t245 * t236;
t297 = t236 * t226 - t234 * t260;
t295 = t244 * t322;
t294 = -t382 + t302;
t292 = t252 * t313;
t291 = t100 - t322;
t290 = g(1) * (-pkin(2) * t234 + t236 * pkin(7));
t259 = cos(qJ(1));
t288 = -g(2) * t259 + t373;
t27 = t255 * t58 - t362;
t286 = -t27 * t92 - t28 * t283;
t87 = t104 - t374;
t88 = -t165 + t105;
t43 = -t250 * t88 + t255 * t87;
t44 = t250 * t87 + t255 * t88;
t109 = t255 * t169 + t170 * t250;
t110 = -t169 * t250 + t170 * t255;
t45 = t255 * t113 + t250 * t114 + t169 * t324 + t170 * t325;
t46 = qJD(5) * t110 - t113 * t250 + t255 * t114;
t284 = -t4 * t109 - t5 * t110 + t27 * t45 - t28 * t46 - t382;
t282 = -t226 * t234 - t236 * t260;
t142 = pkin(4) * t169 - t226;
t281 = t218 + t335 - t372;
t280 = t76 * t113 - t77 * t114 + t169 * t296 - t26 * t170 - t382;
t101 = pkin(3) * t314 - t226 * t242 - t335;
t47 = pkin(4) * t68 + t101;
t279 = -g(1) * t354 + g(2) * t353 + t47 * t110 - t97 * t45;
t278 = -g(1) * t350 + g(2) * t349 + t101 * t170 - t146 * t113;
t277 = g(1) * t352 - g(2) * t351 + t47 * t109 + t97 * t46;
t276 = g(1) * t348 - g(2) * t347 + t101 * t169 + t146 * t114;
t128 = t252 * t305 + t257 * t320;
t129 = -t252 * t320 + t257 * t305;
t48 = t256 * t128 + t251 * t129 + t166 * t326 - t167 * t327;
t275 = t229 - t322;
t274 = -t185 * t244 - t148 + t382;
t261 = qJD(3) ^ 2;
t273 = pkin(7) * t261 - t295 - t379;
t272 = g(1) * t347 + g(2) * t348 + g(3) * t233 + t143 * t146 + t296;
t227 = -pkin(2) - t380;
t271 = pkin(1) * t315 + t224 * t261 + t227 * t242;
t270 = -pkin(7) * qJDD(3) + (t321 - t378) * qJD(3);
t49 = -qJD(4) * t105 - t128 * t251 + t256 * t129;
t268 = -qJDD(3) * t224 + (t227 * t244 - t320) * qJD(3);
t263 = g(1) * t349 + g(2) * t350 - t145 * t146 + t26 - t371;
t239 = t259 * pkin(1);
t231 = qJDD(5) + t241;
t198 = -t226 - t380;
t197 = qJDD(3) * t257 - t252 * t261;
t196 = qJDD(3) * t252 + t257 * t261;
t177 = t230 + t229;
t158 = pkin(3) * t343 + t225 * t250;
t157 = -pkin(3) * t345 + t225 * t255;
t150 = t242 * t247 - 0.2e1 * t292;
t149 = t242 * t246 + 0.2e1 * t292;
t131 = t142 - t380;
t116 = -0.2e1 * qJD(3) * t244 * t334 + 0.2e1 * t252 * t339;
t115 = pkin(3) * t346 + pkin(4) * t145;
t96 = t100 + t230;
t84 = -t114 * t243 - t169 * t241;
t83 = -t113 * t243 + t170 * t241;
t71 = -t143 ^ 2 + t145 ^ 2;
t52 = -t293 + (t143 - t319) * t243;
t42 = t49 + t376;
t41 = -t111 + t48;
t40 = t114 * t143 + t169 * t68;
t39 = -t113 * t145 - t170 * t67;
t34 = -t109 * t231 - t232 * t46;
t33 = t110 * t231 - t232 * t45;
t30 = t255 * t60 - t362;
t29 = -t250 * t60 - t361;
t14 = t113 * t143 - t114 * t145 + t169 * t67 - t170 * t68;
t9 = -t109 * t269 + t46 * t92;
t8 = -t110 * t23 + t283 * t45;
t7 = -qJD(5) * t44 - t250 * t41 + t255 * t42;
t6 = qJD(5) * t43 + t250 * t42 + t255 * t41;
t1 = t109 * t23 + t110 * t269 + t283 * t46 + t45 * t92;
t2 = [0, 0, 0, 0, 0, qJDD(1), t288, g(1) * t259 + g(2) * t254, 0, 0, 0, 0, 0, 0, 0, t242, (t242 * t258 - t315) * pkin(1) + t281, ((-qJDD(1) - t242) * t253 + (-qJD(1) - t244) * t330) * pkin(1) + t382, 0, (t288 + (t253 ^ 2 + t258 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t149, t116, t196, t150, t197, 0, t268 * t252 + (-t271 - t386) * t257 + t338, t268 * t257 + (t271 - t218) * t252 + t317, pkin(1) * qJD(2) * t387 + t224 * t301 + t294, t147 * t227 - t290 - g(2) * (t239 + t337) + t224 * t302 + (t384 * qJD(2) + t373) * pkin(1), t39, t14, t83, t40, t84, 0, t104 * t241 + t143 * t177 + t198 * t68 + t243 * t49 + t276, -t105 * t241 + t145 * t177 - t198 * t67 - t243 * t48 + t278, t104 * t67 - t105 * t68 - t143 * t48 - t145 * t49 + t280, -t296 * t105 + t77 * t48 + t26 * t104 + t76 * t49 + t101 * t198 + t146 * t177 - g(1) * (t282 - t381) - g(2) * (t239 + t297), t8, t1, t33, t9, t34, 0, -t131 * t269 + t231 * t43 + t232 * t7 + t92 * t96 + t277, -t131 * t23 - t231 * t44 - t232 * t6 - t283 * t96 + t279, t23 * t43 + t269 * t44 + t283 * t7 - t6 * t92 + t284, t4 * t44 + t28 * t6 + t5 * t43 + t27 * t7 + t47 * t131 + t97 * t96 - g(1) * (t298 - t381) - g(2) * (t239 + t299); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t281 + t295, (-t323 + (-qJD(2) + t244) * t332) * pkin(1) + t382, 0, 0, t149, t116, t196, t150, t197, 0, t270 * t252 + (-t273 - t386) * t257 + t338, t270 * t257 + (t273 - t218) * t252 + t317, pkin(7) * t301 - t363 * t387 + t294, -t147 * pkin(2) + pkin(7) * t302 - g(2) * t337 - t384 * t363 - t290, t39, t14, t83, t40, t84, 0, t126 * t241 + t275 * t143 - t226 * t68 + t357 * t243 + t276, -t127 * t241 + t275 * t145 + t226 * t67 - t358 * t243 + t278, t126 * t67 - t127 * t68 - t358 * t143 - t357 * t145 + t280, -g(1) * t282 - g(2) * t297 - t101 * t226 + t26 * t126 - t127 * t296 + t275 * t146 + t357 * t76 + t358 * t77, t8, t1, t33, t9, t34, 0, -t142 * t269 + t231 * t56 + t232 * t364 + t291 * t92 + t277, -t142 * t23 - t231 * t57 - t232 * t365 - t283 * t291 + t279, t23 * t56 + t269 * t57 + t283 * t364 - t365 * t92 + t284, -g(1) * t298 - g(2) * t299 + t47 * t142 + t27 * t364 + t28 * t365 + t291 * t97 + t4 * t57 + t5 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, t334 * t240, t342, t318, t339, qJDD(3), t252 * t274 - t370, g(3) * t252 + t257 * t274, 0, 0, t356, t71, t52, -t356, -t287, t241, -t243 * t79 + (-t143 * t346 + t241 * t256 - t243 * t327) * pkin(3) + t263, t243 * t80 + (-t145 * t346 - t241 * t251 - t243 * t326) * pkin(3) + t272, (t77 + t79) * t145 + (-t76 + t80) * t143 + (-t251 * t68 + t256 * t67 + (-t143 * t256 + t145 * t251) * qJD(4)) * pkin(3), -t76 * t79 - t77 * t80 + (-t370 - t296 * t251 + t256 * t26 + (-t251 * t76 + t256 * t77) * qJD(4) + (-t146 * t244 + t382) * t252) * pkin(3), -t390, t22, t18, t390, t19, t231, -t115 * t92 + t157 * t231 + t232 * t359 + t264, t115 * t283 - t158 * t231 - t232 * t360 + t266, t157 * t23 + t158 * t269 + t283 * t359 - t360 * t92 + t286, t4 * t158 + t5 * t157 - t97 * t115 - g(3) * t310 + t360 * t28 + t359 * t27 - t382 * (-pkin(3) * t252 - pkin(4) * t233); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t356, t71, t52, -t356, -t287, t241, t243 * t77 + t263, t243 * t76 + t272, 0, 0, -t390, t22, t18, t390, t19, t231, -t232 * t29 + (-t145 * t92 + t231 * t255 - t232 * t325) * pkin(4) + t264, t232 * t30 + (t145 * t283 - t231 * t250 - t232 * t324) * pkin(4) + t266, -t29 * t283 + t30 * t92 + (t23 * t255 + t269 * t250 + (-t250 * t283 - t255 * t92) * qJD(5)) * pkin(4) + t286, -t27 * t29 - t28 * t30 + (-t371 - t145 * t97 + t250 * t4 + t255 * t5 + t382 * t233 + (-t250 * t27 + t255 * t28) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t390, t22, t18, t390, t19, t231, t232 * t28 + t264, t232 * t27 + t266, 0, 0;];
tau_reg = t2;
