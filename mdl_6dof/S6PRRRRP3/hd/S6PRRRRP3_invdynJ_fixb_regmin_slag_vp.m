% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:11
% EndTime: 2019-03-09 00:11:26
% DurationCPUTime: 6.38s
% Computational Cost: add. (5643->495), mult. (12905->693), div. (0->0), fcn. (9970->14), ass. (0->245)
t213 = sin(qJ(3));
t217 = cos(qJ(3));
t252 = pkin(3) * t213 - pkin(9) * t217;
t156 = t252 * qJD(3);
t212 = sin(qJ(4));
t214 = sin(qJ(2));
t216 = cos(qJ(4));
t302 = qJD(3) * t213;
t209 = sin(pkin(6));
t310 = qJD(1) * t209;
t218 = cos(qJ(2));
t318 = t217 * t218;
t351 = pkin(8) * t212;
t388 = (-t212 * t318 + t214 * t216) * t310 - t216 * t156 - t302 * t351;
t162 = -pkin(3) * t217 - pkin(9) * t213 - pkin(2);
t296 = qJD(4) * t216;
t387 = -(t212 * t214 + t216 * t318) * t310 + t212 * t156 + t162 * t296;
t319 = t216 * t217;
t190 = pkin(8) * t319;
t245 = pkin(4) * t213 - pkin(10) * t319;
t386 = -t245 * qJD(3) - (-t190 + (pkin(10) * t213 - t162) * t212) * qJD(4) + t388;
t298 = qJD(4) * t212;
t301 = qJD(3) * t216;
t300 = qJD(3) * t217;
t275 = t212 * t300;
t376 = t213 * t296 + t275;
t385 = -t376 * pkin(10) + (-t213 * t301 - t217 * t298) * pkin(8) + t387;
t305 = qJD(2) * t217;
t279 = t212 * t305;
t357 = pkin(9) + pkin(10);
t285 = qJD(4) * t357;
t153 = t252 * qJD(2);
t158 = qJD(2) * pkin(8) + t214 * t310;
t210 = cos(pkin(6));
t323 = t210 * t217;
t366 = qJD(1) * t323 - t213 * t158;
t316 = t212 * t153 + t216 * t366;
t378 = -pkin(10) * t279 + t212 * t285 + t316;
t135 = t216 * t153;
t384 = qJD(2) * t245 - t212 * t366 + t216 * t285 + t135;
t306 = qJD(2) * t213;
t145 = -t212 * t306 + t301;
t303 = qJD(3) * t212;
t146 = t216 * t306 + t303;
t211 = sin(qJ(5));
t215 = cos(qJ(5));
t247 = t145 * t211 + t215 * t146;
t358 = t247 ^ 2;
t81 = -t215 * t145 + t146 * t211;
t79 = t81 ^ 2;
t383 = -t79 + t358;
t333 = sin(pkin(11));
t258 = t333 * t214;
t334 = cos(pkin(11));
t259 = t334 * t218;
t121 = -t210 * t259 + t258;
t257 = t333 * t218;
t260 = t334 * t214;
t123 = t210 * t257 + t260;
t251 = g(1) * t123 + g(2) * t121;
t324 = t209 * t218;
t382 = -g(3) * t324 + t251;
t380 = qJ(6) * t81;
t379 = t247 * t81;
t148 = t211 * t216 + t212 * t215;
t288 = qJD(4) + qJD(5);
t92 = t288 * t148;
t336 = t148 * t305 - t92;
t322 = t211 * t212;
t147 = -t215 * t216 + t322;
t294 = qJD(5) * t215;
t377 = -t147 * t305 - t215 * t296 - t216 * t294 + t288 * t322;
t122 = t210 * t260 + t257;
t124 = -t210 * t258 + t259;
t250 = g(1) * t124 + g(2) * t122;
t325 = t209 * t214;
t230 = -g(3) * t325 - t250;
t292 = qJD(2) * qJD(3);
t271 = t217 * t292;
t289 = t213 * qJDD(2);
t375 = qJD(3) * qJD(4) + t271 + t289;
t188 = -qJD(4) + t305;
t177 = -qJD(5) + t188;
t297 = qJD(4) * t213;
t269 = qJD(2) * t297;
t254 = t212 * t375 + t216 * t269;
t232 = t216 * qJDD(3) - t254;
t295 = qJD(5) * t211;
t75 = (qJDD(3) - t269) * t212 + t375 * t216;
t22 = -t145 * t294 + t146 * t295 - t211 * t232 - t215 * t75;
t374 = -t177 * t81 - t22;
t129 = t210 * t213 + t217 * t325;
t208 = qJ(4) + qJ(5);
t201 = sin(t208);
t202 = cos(t208);
t309 = qJD(1) * t213;
t182 = t210 * t309;
t111 = t217 * t158 + t182;
t105 = qJD(3) * pkin(9) + t111;
t308 = qJD(1) * t218;
t284 = t209 * t308;
t113 = qJD(2) * t162 - t284;
t293 = qJD(1) * qJD(2);
t115 = qJDD(2) * pkin(8) + (qJDD(1) * t214 + t218 * t293) * t209;
t290 = qJDD(1) * t210;
t268 = t213 * t290;
t49 = qJDD(3) * pkin(9) + qJD(3) * t366 + t115 * t217 + t268;
t272 = t214 * t293;
t249 = -qJDD(1) * t324 + t209 * t272;
t68 = qJD(2) * t156 + qJDD(2) * t162 + t249;
t287 = t113 * t296 + t212 * t68 + t216 * t49;
t237 = -t105 * t298 + t287;
t12 = pkin(10) * t232 + t237;
t60 = -t105 * t212 + t216 * t113;
t41 = -pkin(10) * t146 + t60;
t31 = -pkin(4) * t188 + t41;
t61 = t105 * t216 + t113 * t212;
t42 = pkin(10) * t145 + t61;
t200 = t217 * qJDD(2);
t142 = t213 * t292 + qJDD(4) - t200;
t67 = t216 * t68;
t224 = -qJD(4) * t61 - t212 * t49 + t67;
t7 = pkin(4) * t142 - pkin(10) * t75 + t224;
t267 = -t215 * t12 - t211 * t7 - t31 * t294 + t42 * t295;
t104 = -qJD(3) * pkin(3) - t366;
t69 = -pkin(4) * t145 + t104;
t262 = t209 * t334;
t94 = t122 * t217 - t213 * t262;
t261 = t209 * t333;
t96 = t124 * t217 + t213 * t261;
t373 = t69 * t81 - g(1) * (-t123 * t201 - t202 * t96) - g(2) * (-t121 * t201 - t202 * t94) - g(3) * (-t129 * t202 + t201 * t324) + t267;
t32 = pkin(5) * t81 + qJD(6) + t69;
t371 = t247 * t32;
t312 = t212 * t162 + t190;
t321 = t212 * t213;
t101 = -pkin(10) * t321 + t312;
t144 = t216 * t162;
t320 = t213 * t216;
t88 = -pkin(10) * t320 + t144 + (-pkin(4) - t351) * t217;
t370 = -t101 * t295 - t211 * t386 + t385 * t215 + t88 * t294;
t369 = t385 * t211 + t215 * t386;
t368 = t384 * t215;
t253 = -t111 + (-t279 + t298) * pkin(4);
t367 = qJ(6) * t247;
t170 = t357 * t212;
t171 = t357 * t216;
t313 = -t211 * t170 + t215 * t171;
t365 = -t170 * t294 - t171 * t295 - t211 * t384 - t215 * t378;
t274 = t216 * t300;
t364 = -t212 * t297 + t274;
t363 = -g(3) * (-t129 * t201 - t202 * t324) - g(2) * (t121 * t202 - t201 * t94) - g(1) * (t123 * t202 - t201 * t96);
t38 = t215 * t42;
t19 = t211 * t31 + t38;
t273 = -t211 * t12 + t215 * t7;
t227 = -qJD(5) * t19 + t273;
t362 = -t69 * t247 + t227 + t363;
t23 = qJD(5) * t247 + t211 * t75 - t215 * t232;
t361 = -t177 * t247 - t23;
t219 = qJD(3) ^ 2;
t360 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t219 + t209 * (-g(3) * t218 + t272) - t249 + t251;
t359 = qJD(4) * (pkin(8) * t188 + t105) + t382;
t119 = t147 * t213;
t339 = t215 * t101 + t211 * t88;
t51 = t211 * t275 + t213 * t92 - t215 * t274;
t353 = pkin(5) * t302 + qJ(6) * t51 - qJD(5) * t339 + qJD(6) * t119 - t369;
t118 = t148 * t213;
t52 = -t295 * t321 + (t288 * t320 + t275) * t215 + t364 * t211;
t352 = -qJ(6) * t52 - qJD(6) * t118 + t370;
t204 = t216 * pkin(4);
t346 = pkin(3) + t204;
t36 = t211 * t42;
t18 = t215 * t31 - t36;
t14 = t18 - t367;
t13 = -pkin(5) * t177 + t14;
t344 = -t14 + t13;
t343 = qJ(6) * t336 - qJD(6) * t147 + t365;
t342 = -pkin(5) * t306 + qJ(6) * t377 - t313 * qJD(5) - qJD(6) * t148 + t211 * t378 - t368;
t341 = t215 * t41 - t36;
t338 = qJD(2) * pkin(2);
t337 = t75 * t212;
t331 = qJDD(3) * pkin(3);
t330 = t145 * t188;
t329 = t146 * t188;
t328 = t146 * t216;
t327 = t201 * t217;
t326 = t202 * t217;
t317 = qJDD(1) - g(3);
t161 = pkin(5) * t202 + t204;
t157 = pkin(4) * t321 + t213 * pkin(8);
t206 = t213 ^ 2;
t311 = -t217 ^ 2 + t206;
t307 = qJD(2) * t209;
t304 = qJD(3) * t145;
t299 = qJD(4) * t188;
t112 = pkin(4) * t376 + pkin(8) * t300;
t282 = t213 * t308;
t281 = t214 * t307;
t280 = t218 * t307;
t278 = t188 * t301;
t270 = t218 * t292;
t265 = -t211 * t41 - t38;
t256 = -t101 * t211 + t215 * t88;
t255 = -t215 * t170 - t171 * t211;
t242 = -t129 * t216 + t212 * t324;
t99 = -t129 * t212 - t216 * t324;
t46 = t211 * t99 - t215 * t242;
t45 = t211 * t242 + t215 * t99;
t220 = qJD(2) ^ 2;
t246 = qJDD(2) * t218 - t214 * t220;
t128 = t213 * t325 - t323;
t240 = t142 * t212 - t188 * t296;
t239 = t142 * t216 + t188 * t298;
t93 = t122 * t213 + t217 * t262;
t95 = t124 * t213 - t217 * t261;
t236 = g(1) * t95 + g(2) * t93 + g(3) * t128;
t235 = g(1) * t96 + g(2) * t94 + g(3) * t129;
t233 = qJD(3) * t182 + t213 * t115 + t158 * t300 - t217 * t290;
t229 = -pkin(9) * t142 - t104 * t188;
t50 = t233 - t331;
t223 = -pkin(9) * t299 - t236 + t50;
t159 = -t284 - t338;
t222 = -pkin(8) * qJDD(3) + (t159 + t284 - t338) * qJD(3);
t26 = -pkin(4) * t232 + t50;
t10 = t23 * pkin(5) + qJDD(6) + t26;
t205 = -qJ(6) - t357;
t196 = pkin(4) * t215 + pkin(5);
t160 = pkin(4) * t212 + pkin(5) * t201;
t152 = pkin(3) + t161;
t137 = qJDD(5) + t142;
t98 = qJD(3) * t129 + t213 * t280;
t97 = -qJD(3) * t128 + t217 * t280;
t71 = -qJ(6) * t147 + t313;
t70 = -qJ(6) * t148 + t255;
t35 = qJD(4) * t99 + t212 * t281 + t216 * t97;
t34 = qJD(4) * t242 - t212 * t97 + t216 * t281;
t28 = -qJ(6) * t118 + t339;
t27 = -pkin(5) * t217 + qJ(6) * t119 + t256;
t17 = t341 - t367;
t16 = t265 + t380;
t15 = t19 - t380;
t9 = -qJD(5) * t46 - t211 * t35 + t215 * t34;
t8 = qJD(5) * t45 + t211 * t34 + t215 * t35;
t2 = -qJ(6) * t23 - qJD(6) * t81 - t267;
t1 = pkin(5) * t137 + qJ(6) * t22 - qJD(6) * t247 + t227;
t3 = [t317, 0, t246 * t209 (-qJDD(2) * t214 - t218 * t220) * t209, 0, 0, 0, 0, 0, -qJD(3) * t98 - qJDD(3) * t128 + (-t213 * t270 + t217 * t246) * t209, -qJD(3) * t97 - qJDD(3) * t129 + (-t213 * t246 - t217 * t270) * t209, 0, 0, 0, 0, 0, -t128 * t232 + t99 * t142 - t98 * t145 - t34 * t188, t128 * t75 + t142 * t242 + t146 * t98 + t188 * t35, 0, 0, 0, 0, 0, t128 * t23 + t137 * t45 - t177 * t9 + t81 * t98, -t128 * t22 - t137 * t46 + t177 * t8 + t247 * t98, t22 * t45 - t23 * t46 - t247 * t9 - t8 * t81, t1 * t45 + t10 * t128 + t13 * t9 + t15 * t8 + t2 * t46 + t32 * t98 - g(3); 0, qJDD(2), t317 * t324 + t251, -t317 * t325 + t250, qJDD(2) * t206 + 0.2e1 * t213 * t271, 0.2e1 * t200 * t213 - 0.2e1 * t292 * t311, qJDD(3) * t213 + t217 * t219, qJDD(3) * t217 - t213 * t219, 0, t222 * t213 + t360 * t217, -t360 * t213 + t222 * t217, t364 * t146 + t75 * t320 (t145 * t216 - t146 * t212) * t300 + (t216 * t232 - t337 + (-t145 * t212 - t328) * qJD(4)) * t213 (-t75 - t278) * t217 + (qJD(3) * t146 + t239) * t213 (t188 * t303 - t232) * t217 + (-t240 + t304) * t213, -t142 * t217 - t188 * t302, t144 * t142 + t388 * t188 + (t162 * t299 + t230) * t212 + (-pkin(8) * t304 - t67 + (-pkin(8) * t142 + qJD(3) * t104 + qJD(4) * t113 + t49) * t212 + t359 * t216) * t217 + (-pkin(8) * t232 + t60 * qJD(3) + t104 * t296 + t145 * t284 + t50 * t212) * t213, -t312 * t142 + t387 * t188 + t230 * t216 + ((pkin(8) * t146 + t104 * t216) * qJD(3) - t359 * t212 + t287) * t217 + (-t146 * t284 - t104 * t298 - t61 * qJD(3) + t50 * t216 + (t75 - t278) * pkin(8)) * t213, t119 * t22 - t247 * t51, t118 * t22 + t119 * t23 - t247 * t52 + t51 * t81, -t119 * t137 + t177 * t51 + t217 * t22 + t247 * t302, -t118 * t137 + t177 * t52 + t217 * t23 - t302 * t81, -t137 * t217 - t177 * t302, t256 * t137 - t273 * t217 + t18 * t302 + t112 * t81 + t157 * t23 + t26 * t118 + t69 * t52 - g(1) * (-t123 * t326 + t124 * t201) - g(2) * (-t121 * t326 + t122 * t201) + t369 * t177 + (t177 * t339 + t19 * t217) * qJD(5) + (-t81 * t282 - g(3) * (t201 * t214 + t202 * t318)) * t209, -t339 * t137 - t267 * t217 - t19 * t302 + t112 * t247 - t157 * t22 - t26 * t119 - t69 * t51 - g(1) * (t123 * t327 + t124 * t202) - g(2) * (t121 * t327 + t122 * t202) + t370 * t177 + (-t247 * t282 - g(3) * (-t201 * t318 + t202 * t214)) * t209, t1 * t119 - t118 * t2 + t13 * t51 - t15 * t52 + t213 * t382 + t22 * t27 - t23 * t28 - t247 * t353 - t352 * t81, t2 * t28 + t1 * t27 + t10 * (pkin(5) * t118 + t157) + t352 * t15 + t353 * t13 + t230 * (pkin(8) + t160) + (pkin(5) * t52 - t309 * t324 + t112) * t32 + t382 * (t152 * t217 - t205 * t213 + pkin(2)); 0, 0, 0, 0, -t213 * t220 * t217, t311 * t220, t289, t200, qJDD(3), qJD(3) * t111 - t159 * t306 - t233 + t236, -t268 + (-qJD(2) * t159 - t115) * t217 + t235, -t188 * t328 + t337 (t75 - t330) * t216 + (t232 + t329) * t212 (-t146 * t213 + t188 * t319) * qJD(2) + t240 (-t188 * t212 * t217 - t145 * t213) * qJD(2) + t239, t188 * t306, -pkin(3) * t254 + t135 * t188 - t60 * t306 + t111 * t145 + (-t188 * t366 + t229) * t212 + (-t223 + t331) * t216, -pkin(3) * t75 - t111 * t146 - t188 * t316 + t212 * t223 + t216 * t229 + t306 * t61, -t148 * t22 - t247 * t377, t147 * t22 - t148 * t23 + t247 * t336 + t377 * t81, t137 * t148 + t177 * t377 - t247 * t306, -t137 * t147 - t177 * t336 + t306 * t81, t177 * t306, t255 * t137 - t346 * t23 + t26 * t147 - t18 * t306 + t253 * t81 - t336 * t69 + (t171 * t294 + (-qJD(5) * t170 - t378) * t211 + t368) * t177 + t236 * t202, -t313 * t137 + t26 * t148 + t177 * t365 + t19 * t306 - t236 * t201 + t22 * t346 + t253 * t247 - t377 * t69, -t1 * t148 + t13 * t377 - t147 * t2 + t15 * t336 + t22 * t70 - t23 * t71 - t247 * t342 - t343 * t81 - t235, t2 * t71 + t1 * t70 + t10 * (pkin(5) * t147 - t346) - g(1) * (-t152 * t95 - t205 * t96) - g(2) * (-t152 * t93 - t205 * t94) - g(3) * (-t128 * t152 - t129 * t205) + (-pkin(5) * t336 + t253) * t32 + t343 * t15 + t342 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146 * t145, -t145 ^ 2 + t146 ^ 2, t75 + t330, t232 - t329, t142, -t61 * t188 - t104 * t146 - g(1) * (t123 * t216 - t212 * t96) - g(2) * (t121 * t216 - t212 * t94) - g(3) * t99 + t224, -t60 * t188 - t104 * t145 - g(1) * (-t123 * t212 - t216 * t96) - g(2) * (-t121 * t212 - t216 * t94) - g(3) * t242 - t237, t379, t383, t374, t361, t137, t265 * t177 + (t137 * t215 - t146 * t81 + t177 * t295) * pkin(4) + t362, -t341 * t177 + (-t211 * t137 - t146 * t247 + t177 * t294) * pkin(4) + t373, -t13 * t81 + t15 * t247 + t16 * t247 + t17 * t81 + t196 * t22 + (-t211 * t23 + (t211 * t247 - t215 * t81) * qJD(5)) * pkin(4), t1 * t196 - t15 * t17 - t13 * t16 - pkin(5) * t371 - g(1) * (t123 * t161 - t160 * t96) - g(2) * (t121 * t161 - t160 * t94) - g(3) * (-t129 * t160 - t161 * t324) + (-t32 * t146 + t2 * t211 + (-t13 * t211 + t15 * t215) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t379, t383, t374, t361, t137, -t177 * t19 + t362, -t177 * t18 + t373, pkin(5) * t22 - t344 * t81, t344 * t15 + (t1 + t363 - t371) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 - t358, t13 * t247 + t15 * t81 + t10 - t236;];
tau_reg  = t3;
