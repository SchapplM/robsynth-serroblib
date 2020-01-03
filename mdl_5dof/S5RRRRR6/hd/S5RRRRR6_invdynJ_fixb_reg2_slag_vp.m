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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:15:24
% EndTime: 2020-01-03 12:15:32
% DurationCPUTime: 3.95s
% Computational Cost: add. (8658->429), mult. (13431->543), div. (0->0), fcn. (9125->16), ass. (0->255)
t255 = cos(qJ(4));
t256 = cos(qJ(3));
t336 = t255 * t256;
t250 = sin(qJ(4));
t251 = sin(qJ(3));
t340 = t250 * t251;
t170 = -t336 + t340;
t259 = -pkin(8) - pkin(7);
t312 = qJD(3) * t259;
t176 = t251 * t312;
t177 = t256 * t312;
t201 = t259 * t251;
t237 = t256 * pkin(8);
t202 = pkin(7) * t256 + t237;
t257 = cos(qJ(2));
t327 = qJD(1) * t257;
t315 = pkin(1) * t327;
t321 = qJD(4) * t255;
t322 = qJD(4) * t250;
t350 = t170 * t315 + t176 * t255 + t177 * t250 + t201 * t321 - t202 * t322;
t127 = t201 * t250 + t202 * t255;
t337 = t255 * t251;
t171 = t250 * t256 + t337;
t349 = -qJD(4) * t127 + t171 * t315 - t250 * t176 + t177 * t255;
t252 = sin(qJ(2));
t356 = pkin(1) * qJD(1);
t316 = t252 * t356;
t363 = t257 * pkin(1);
t330 = -qJD(2) * t316 + qJDD(1) * t363;
t241 = qJDD(1) + qJDD(2);
t366 = t241 * pkin(2);
t147 = -t330 - t366;
t248 = qJ(1) + qJ(2);
t232 = sin(t248);
t234 = cos(t248);
t285 = -g(2) * t234 - g(3) * t232;
t383 = t147 - t285;
t243 = qJD(1) + qJD(2);
t314 = t243 * t340;
t143 = -t243 * t336 + t314;
t145 = t171 * t243;
t249 = sin(qJ(5));
t254 = cos(qJ(5));
t278 = t143 * t249 - t145 * t254;
t92 = t143 * t254 + t145 * t249;
t382 = t92 * t278;
t242 = qJD(3) + qJD(4);
t280 = t242 * t340;
t323 = qJD(3) * t256;
t378 = -t255 * t323 - t256 * t321;
t113 = t280 + t378;
t369 = t113 * pkin(9);
t381 = t369 + t349;
t114 = t242 * t171;
t111 = t114 * pkin(9);
t380 = -t111 + t350;
t245 = t251 ^ 2;
t246 = t256 ^ 2;
t328 = t245 + t246;
t376 = t328 * t257;
t379 = t243 * t376;
t318 = qJDD(1) * t252;
t325 = qJD(2) * t257;
t148 = t241 * pkin(7) + (qJD(1) * t325 + t318) * pkin(1);
t299 = t328 * t148;
t215 = g(3) * t234;
t331 = -g(2) * t232 + t215;
t22 = t278 ^ 2 - t92 ^ 2;
t319 = qJD(5) * t254;
t320 = qJD(5) * t249;
t335 = t256 * t241;
t288 = -t241 * t337 + t243 * t378 - t250 * t335;
t67 = t243 * t280 + t288;
t338 = t251 * t241;
t282 = t250 * t338 - t255 * t335;
t68 = t114 * t243 + t282;
t23 = t143 * t319 + t145 * t320 + t249 * t68 + t254 * t67;
t230 = qJD(5) + t242;
t18 = t230 * t92 - t23;
t247 = qJ(3) + qJ(4);
t235 = qJ(5) + t247;
t220 = sin(t235);
t221 = cos(t235);
t240 = qJDD(3) + qJDD(4);
t183 = pkin(7) * t243 + t316;
t309 = t243 * t323;
t73 = -t183 * t323 + qJDD(3) * pkin(3) - t251 * t148 + (-t309 - t338) * pkin(8);
t305 = pkin(8) * t243 + t183;
t133 = t305 * t256;
t124 = t255 * t133;
t132 = t305 * t251;
t125 = qJD(3) * pkin(3) - t132;
t77 = t125 * t250 + t124;
t324 = qJD(3) * t251;
t310 = t243 * t324;
t78 = -t183 * t324 + t256 * t148 + (-t310 + t335) * pkin(8);
t26 = -qJD(4) * t77 - t250 * t78 + t255 * t73;
t16 = t240 * pkin(4) + t67 * pkin(9) + t26;
t294 = -t125 * t321 + t133 * t322 - t250 * t73 - t255 * t78;
t17 = -pkin(9) * t68 - t294;
t141 = t145 * pkin(9);
t122 = t250 * t133;
t76 = t125 * t255 - t122;
t60 = -t141 + t76;
t58 = pkin(4) * t242 + t60;
t368 = t143 * pkin(9);
t61 = t77 - t368;
t4 = (qJD(5) * t58 + t17) * t254 + t249 * t16 - t61 * t320;
t364 = t256 * pkin(3);
t224 = pkin(2) + t364;
t146 = -t224 * t243 - t315;
t97 = t143 * pkin(4) + t146;
t265 = g(1) * t220 - t221 * t331 + t97 * t92 - t4;
t345 = t220 * t234;
t346 = t220 * t232;
t353 = t254 * t61;
t28 = t249 * t58 + t353;
t5 = -qJD(5) * t28 + t16 * t254 - t249 * t17;
t263 = -g(1) * t221 + g(2) * t346 - g(3) * t345 + t278 * t97 + t5;
t365 = t243 * pkin(2);
t184 = -t315 - t365;
t377 = t183 * t376 + t184 * t252;
t269 = qJD(5) * t278 + t249 * t67 - t254 * t68;
t19 = -t230 * t278 + t269;
t222 = pkin(1) * t252 + pkin(7);
t359 = -pkin(8) - t222;
t167 = t359 * t251;
t168 = t222 * t256 + t237;
t105 = t167 * t250 + t168 * t255;
t260 = qJD(3) ^ 2;
t375 = pkin(7) * t260 - t366;
t233 = cos(t247);
t373 = g(1) * t233;
t372 = g(1) * t256;
t367 = t171 * pkin(9);
t126 = t201 * t255 - t202 * t250;
t98 = t126 - t367;
t166 = t170 * pkin(9);
t99 = -t166 + t127;
t56 = -t249 * t99 + t254 * t98;
t358 = t56 * qJD(5) + t249 * t381 + t254 * t380;
t57 = t249 * t98 + t254 * t99;
t357 = -t57 * qJD(5) - t249 * t380 + t254 * t381;
t355 = pkin(1) * qJD(2);
t354 = t249 * t61;
t223 = pkin(3) * t255 + pkin(4);
t341 = t249 * t250;
t79 = t132 * t250 - t124;
t62 = t79 + t368;
t80 = -t132 * t255 - t122;
t63 = -t141 + t80;
t352 = -t249 * t62 - t254 * t63 + t223 * t319 + (-t250 * t320 + (t254 * t255 - t341) * qJD(4)) * pkin(3);
t339 = t250 * t254;
t351 = t249 * t63 - t254 * t62 - t223 * t320 + (-t250 * t319 + (-t249 * t255 - t339) * qJD(4)) * pkin(3);
t348 = t145 * t143;
t231 = sin(t247);
t344 = t231 * t232;
t343 = t231 * t234;
t342 = t243 * t251;
t306 = pkin(4) * t233 + t364;
t175 = pkin(2) + t306;
t244 = -pkin(9) + t259;
t334 = t175 * t232 + t234 * t244;
t333 = t224 * t232 + t234 * t259;
t332 = pkin(2) * t234 + pkin(7) * t232;
t329 = t245 - t246;
t326 = qJD(2) * t252;
t317 = pkin(1) * t325;
t227 = pkin(3) * t324;
t239 = t243 ^ 2;
t313 = t251 * t239 * t256;
t311 = t243 * t326;
t100 = pkin(4) * t114 + t227;
t302 = qJD(3) * t359;
t298 = t328 * t241;
t104 = t167 * t255 - t168 * t250;
t296 = t175 * t234 - t232 * t244;
t295 = t224 * t234 - t232 * t259;
t110 = -t170 * t249 + t171 * t254;
t45 = t113 * t254 + t114 * t249 + t170 * t319 + t171 * t320;
t101 = pkin(3) * t310 - t224 * t241 - t330;
t47 = t68 * pkin(4) + t101;
t293 = g(2) * t345 + g(3) * t346 + t110 * t47 - t45 * t97;
t292 = g(2) * t343 + g(3) * t344 + t101 * t171 - t113 * t146;
t291 = t243 * t316;
t290 = t331 + t299;
t289 = t184 * t323 + t251 * t383;
t287 = t251 * t309;
t286 = t100 - t316;
t253 = sin(qJ(1));
t258 = cos(qJ(1));
t283 = -g(2) * t258 - g(3) * t253;
t27 = t254 * t58 - t354;
t281 = -t27 * t92 - t278 * t28;
t87 = t104 - t367;
t88 = -t166 + t105;
t43 = -t249 * t88 + t254 * t87;
t44 = t249 * t87 + t254 * t88;
t109 = t170 * t254 + t171 * t249;
t46 = qJD(5) * t110 - t249 * t113 + t114 * t254;
t279 = -t109 * t4 - t110 * t5 + t27 * t45 - t28 * t46 + t331;
t142 = pkin(4) * t170 - t224;
t277 = t113 * t76 - t114 * t77 + t170 * t294 - t171 * t26 + t331;
t128 = t251 * t302 + t256 * t317;
t129 = -t251 * t317 + t256 * t302;
t48 = t128 * t255 + t129 * t250 + t167 * t321 - t168 * t322;
t276 = -t316 + t227;
t275 = -t184 * t243 - t148 - t331;
t274 = t109 * t47 + t221 * t285 + t46 * t97;
t273 = t101 * t170 + t114 * t146 + t233 * t285;
t272 = t285 + t291;
t225 = -pkin(2) - t363;
t271 = pkin(1) * t311 + t222 * t260 + t225 * t241;
t270 = -pkin(7) * qJDD(3) + (t315 - t365) * qJD(3);
t49 = -qJD(4) * t105 - t250 * t128 + t129 * t255;
t268 = -qJDD(3) * t222 + (t225 * t243 - t317) * qJD(3);
t267 = g(1) * t231 + t146 * t143 - t233 * t331 + t294;
t262 = g(2) * t344 - g(3) * t343 - t146 * t145 + t26 - t373;
t238 = t258 * pkin(1);
t236 = t253 * pkin(1);
t229 = qJDD(5) + t240;
t228 = pkin(1) * t326;
t213 = t232 * pkin(2);
t196 = -t224 - t363;
t195 = qJDD(3) * t256 - t251 * t260;
t194 = qJDD(3) * t251 + t256 * t260;
t178 = t228 + t227;
t160 = t184 * t324;
t159 = pkin(3) * t339 + t223 * t249;
t158 = -pkin(3) * t341 + t223 * t254;
t150 = t241 * t246 - 0.2e1 * t287;
t149 = t241 * t245 + 0.2e1 * t287;
t131 = t142 - t363;
t116 = -0.2e1 * qJD(3) * t243 * t329 + 0.2e1 * t251 * t335;
t115 = pkin(3) * t342 + pkin(4) * t145;
t96 = t100 + t228;
t84 = -t114 * t242 - t170 * t240;
t83 = -t113 * t242 + t171 * t240;
t71 = -t143 ^ 2 + t145 ^ 2;
t52 = -t288 + (t143 - t314) * t242;
t42 = t49 + t369;
t41 = -t111 + t48;
t40 = t114 * t143 + t170 * t68;
t39 = -t113 * t145 - t171 * t67;
t34 = -t109 * t229 - t230 * t46;
t33 = t110 * t229 - t230 * t45;
t30 = t254 * t60 - t354;
t29 = -t249 * t60 - t353;
t14 = t113 * t143 - t114 * t145 + t170 * t67 - t171 * t68;
t9 = -t109 * t269 + t46 * t92;
t8 = -t110 * t23 + t278 * t45;
t7 = -qJD(5) * t44 - t249 * t41 + t254 * t42;
t6 = qJD(5) * t43 + t249 * t42 + t254 * t41;
t1 = t109 * t23 + t110 * t269 + t278 * t46 + t45 * t92;
t2 = [0, 0, 0, 0, 0, qJDD(1), t283, g(2) * t253 - g(3) * t258, 0, 0, 0, 0, 0, 0, 0, t241, (t241 * t257 - t311) * pkin(1) + t285 + t330, ((-qJDD(1) - t241) * t252 + (-qJD(1) - t243) * t325) * pkin(1) - t331, 0, (t283 + (t252 ^ 2 + t257 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t149, t116, t194, t150, t195, 0, t160 + t268 * t251 + (-t271 - t383) * t256, t251 * t271 + t256 * t268 + t289, t222 * t298 + t355 * t379 + t290, t147 * t225 - g(2) * (t238 + t332) - g(3) * (-pkin(7) * t234 + t213 + t236) + t222 * t299 + t377 * t355, t39, t14, t83, t40, t84, 0, t104 * t240 + t178 * t143 + t196 * t68 + t49 * t242 + t273, -t105 * t240 + t145 * t178 - t196 * t67 - t242 * t48 + t292, t104 * t67 - t105 * t68 - t143 * t48 - t145 * t49 + t277, -t294 * t105 + t77 * t48 + t26 * t104 + t76 * t49 + t101 * t196 + t146 * t178 - g(2) * (t238 + t295) - g(3) * (t236 + t333), t8, t1, t33, t9, t34, 0, -t131 * t269 + t43 * t229 + t7 * t230 + t96 * t92 + t274, -t131 * t23 - t229 * t44 - t230 * t6 - t278 * t96 + t293, t23 * t43 + t269 * t44 + t278 * t7 - t6 * t92 + t279, t4 * t44 + t28 * t6 + t5 * t43 + t27 * t7 + t47 * t131 + t97 * t96 - g(2) * (t238 + t296) - g(3) * (t236 + t334); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t272 + t330, (-t318 + (-qJD(2) + t243) * t327) * pkin(1) - t331, 0, 0, t149, t116, t194, t150, t195, 0, t160 + t270 * t251 + (-t147 + t272 - t375) * t256, t270 * t256 + (-t291 + t375) * t251 + t289, pkin(7) * t298 - t356 * t379 + t290, -t147 * pkin(2) - g(2) * t332 - g(3) * t213 + (t299 + t215) * pkin(7) - t377 * t356, t39, t14, t83, t40, t84, 0, t126 * t240 + t143 * t276 - t224 * t68 + t242 * t349 + t273, -t127 * t240 + t145 * t276 + t224 * t67 - t242 * t350 + t292, t126 * t67 - t127 * t68 - t143 * t350 - t145 * t349 + t277, -g(2) * t295 - g(3) * t333 - t101 * t224 + t26 * t126 - t127 * t294 + t146 * t276 + t349 * t76 + t350 * t77, t8, t1, t33, t9, t34, 0, -t142 * t269 + t56 * t229 + t230 * t357 + t286 * t92 + t274, -t142 * t23 - t57 * t229 - t230 * t358 - t278 * t286 + t293, t56 * t23 + t269 * t57 + t278 * t357 - t358 * t92 + t279, -g(2) * t296 - g(3) * t334 + t47 * t142 + t27 * t357 + t28 * t358 + t286 * t97 + t4 * t57 + t5 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, t329 * t239, t338, t313, t335, qJDD(3), t251 * t275 - t372, g(1) * t251 + t256 * t275, 0, 0, t348, t71, t52, -t348, -t282, t240, -t79 * t242 + (-t143 * t342 + t240 * t255 - t242 * t322) * pkin(3) + t262, t80 * t242 + (-t145 * t342 - t240 * t250 - t242 * t321) * pkin(3) + t267, (t77 + t79) * t145 + (-t76 + t80) * t143 + (-t250 * t68 + t255 * t67 + (-t143 * t255 + t145 * t250) * qJD(4)) * pkin(3), -t76 * t79 - t77 * t80 + (-t372 - t294 * t250 + t255 * t26 + (-t250 * t76 + t255 * t77) * qJD(4) + (-t146 * t243 - t331) * t251) * pkin(3), -t382, t22, t18, t382, t19, t229, -t115 * t92 + t158 * t229 + t230 * t351 + t263, t115 * t278 - t159 * t229 - t230 * t352 + t265, t158 * t23 + t159 * t269 + t278 * t351 - t352 * t92 + t281, t4 * t159 + t5 * t158 - t97 * t115 - g(1) * t306 + t352 * t28 + t351 * t27 + t331 * (-pkin(3) * t251 - pkin(4) * t231); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, t71, t52, -t348, -t282, t240, t77 * t242 + t262, t242 * t76 + t267, 0, 0, -t382, t22, t18, t382, t19, t229, -t29 * t230 + (-t145 * t92 + t229 * t254 - t230 * t320) * pkin(4) + t263, t30 * t230 + (t145 * t278 - t229 * t249 - t230 * t319) * pkin(4) + t265, -t29 * t278 + t30 * t92 + (t23 * t254 + t269 * t249 + (-t249 * t278 - t254 * t92) * qJD(5)) * pkin(4) + t281, -t27 * t29 - t28 * t30 + (-t373 - t145 * t97 + t249 * t4 + t254 * t5 - t331 * t231 + (-t249 * t27 + t254 * t28) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t382, t22, t18, t382, t19, t229, t28 * t230 + t263, t27 * t230 + t265, 0, 0;];
tau_reg = t2;
