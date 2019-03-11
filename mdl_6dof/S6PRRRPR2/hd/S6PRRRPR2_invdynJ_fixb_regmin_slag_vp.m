% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:09:02
% EndTime: 2019-03-08 23:09:17
% DurationCPUTime: 6.83s
% Computational Cost: add. (6967->496), mult. (15882->705), div. (0->0), fcn. (12765->18), ass. (0->263)
t221 = qJD(3) + qJD(4);
t229 = sin(qJ(4));
t233 = cos(qJ(3));
t365 = cos(qJ(4));
t298 = qJD(2) * t365;
t230 = sin(qJ(3));
t318 = qJD(2) * t230;
t379 = -t229 * t318 + t233 * t298;
t383 = t221 * t379;
t161 = qJD(6) - t379;
t324 = t229 * t233;
t172 = -qJD(2) * t324 - t230 * t298;
t225 = sin(pkin(12));
t227 = cos(pkin(12));
t143 = t172 * t227 - t221 * t225;
t145 = t172 * t225 + t221 * t227;
t228 = sin(qJ(6));
t232 = cos(qJ(6));
t374 = t143 * t228 + t145 * t232;
t382 = t374 * t161;
t231 = sin(qJ(2));
t226 = sin(pkin(6));
t319 = qJD(1) * t226;
t303 = t231 * t319;
t355 = qJD(3) * pkin(3);
t253 = t230 * t355 - t303;
t366 = pkin(9) + pkin(8);
t304 = qJD(3) * t366;
t179 = t230 * t304;
t180 = t233 * t304;
t259 = -t229 * t230 + t233 * t365;
t190 = t366 * t230;
t191 = t366 * t233;
t260 = -t190 * t365 - t191 * t229;
t234 = cos(qJ(2));
t302 = t234 * t319;
t381 = -qJD(4) * t260 + t179 * t365 + t180 * t229 + t259 * t302;
t135 = t221 * t259;
t178 = t230 * t365 + t324;
t136 = t221 * t178;
t380 = pkin(4) * t136 - qJ(5) * t135 - qJD(5) * t178 + t253;
t378 = t161 - qJD(6);
t176 = t225 * t232 + t227 * t228;
t167 = t176 * qJD(6);
t377 = -t176 * t379 + t167;
t346 = cos(pkin(6));
t288 = qJD(1) * t346;
t289 = qJD(2) * t366 + t303;
t140 = -t230 * t289 + t233 * t288;
t132 = t140 + t355;
t141 = t230 * t288 + t233 * t289;
t297 = qJD(4) * t365;
t316 = qJD(4) * t229;
t286 = t346 * qJDD(1);
t196 = t233 * t286;
t312 = qJD(1) * qJD(2);
t154 = qJDD(2) * pkin(8) + (qJDD(1) * t231 + t234 * t312) * t226;
t287 = pkin(9) * qJDD(2) + t154;
t64 = qJDD(3) * pkin(3) - qJD(3) * t141 - t230 * t287 + t196;
t67 = qJD(3) * t140 + t230 * t286 + t233 * t287;
t285 = -t132 * t316 - t141 * t297 - t229 * t67 + t365 * t64;
t219 = qJDD(3) + qJDD(4);
t368 = -pkin(4) * t219 + qJDD(5);
t19 = -t285 + t368;
t345 = cos(pkin(11));
t277 = t346 * t345;
t344 = sin(pkin(11));
t163 = t231 * t277 + t234 * t344;
t224 = qJ(3) + qJ(4);
t216 = sin(t224);
t217 = cos(t224);
t291 = t226 * t345;
t124 = t163 * t216 + t217 * t291;
t276 = t346 * t344;
t165 = -t231 * t276 + t234 * t345;
t290 = t226 * t344;
t126 = t165 * t216 - t217 * t290;
t327 = t226 * t231;
t157 = t216 * t327 - t217 * t346;
t255 = g(1) * t126 + g(2) * t124 + g(3) * t157;
t376 = -t19 + t255;
t267 = t143 * t232 - t145 * t228;
t375 = t161 * t267;
t175 = t225 * t228 - t227 * t232;
t321 = t161 * t175;
t125 = t163 * t217 - t216 * t291;
t127 = t165 * t217 + t216 * t290;
t158 = t216 * t346 + t217 * t327;
t254 = g(1) * t127 + g(2) * t125 + g(3) * t158;
t243 = t132 * t297 - t141 * t316 + t229 * t64 + t365 * t67;
t17 = qJ(5) * t219 + qJD(5) * t221 + t243;
t213 = pkin(3) * t233 + pkin(2);
t296 = t231 * t312;
t325 = t226 * t234;
t275 = -qJDD(1) * t325 + t226 * t296;
t311 = qJD(2) * qJD(3);
t295 = t230 * t311;
t123 = pkin(3) * t295 - qJDD(2) * t213 + t275;
t292 = qJDD(2) * t365;
t309 = t233 * qJDD(2);
t94 = t229 * t309 + t230 * t292 + t383;
t310 = t230 * qJDD(2);
t274 = t229 * t310 - t233 * t292;
t95 = qJD(2) * t136 + t274;
t28 = pkin(4) * t95 - qJ(5) * t94 + qJD(5) * t172 + t123;
t7 = t17 * t227 + t225 * t28;
t373 = t227 * t7 - t254;
t350 = t225 * t381 + t227 * t380;
t349 = t225 * t380 - t227 * t381;
t197 = pkin(3) * t297 + qJD(5);
t118 = -pkin(4) * t172 - qJ(5) * t379;
t104 = pkin(3) * t318 + t118;
t130 = t229 * t141;
t79 = t140 * t365 - t130;
t45 = t104 * t227 - t225 * t79;
t372 = -t197 * t225 - t45;
t46 = t104 * t225 + t227 * t79;
t371 = -t197 * t227 + t46;
t73 = t132 * t365 - t130;
t48 = t118 * t225 + t227 * t73;
t370 = -qJD(5) * t227 + t48;
t131 = t365 * t141;
t78 = t140 * t229 + t131;
t283 = pkin(3) * t316 - t78;
t150 = -t190 * t229 + t191 * t365;
t347 = qJD(4) * t150 - t178 * t302 - t179 * t229 + t180 * t365;
t162 = t231 * t344 - t234 * t277;
t164 = t231 * t345 + t234 * t276;
t280 = g(1) * t164 + g(2) * t162;
t251 = g(3) * t325 - t280;
t369 = t251 * t216;
t236 = qJD(3) ^ 2;
t367 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t236 + t226 * (-g(3) * t234 + t296) - t275 + t280;
t80 = -t219 * t227 + t225 * t94;
t81 = t219 * t225 + t227 * t94;
t22 = -qJD(6) * t267 + t228 * t81 + t232 * t80;
t358 = g(3) * t226;
t357 = t227 * pkin(5);
t218 = t227 * pkin(10);
t74 = t132 * t229 + t131;
t65 = qJ(5) * t221 + t74;
t159 = -qJD(2) * t213 - t302;
t93 = -pkin(4) * t379 + qJ(5) * t172 + t159;
t37 = t225 * t93 + t227 * t65;
t356 = qJD(2) * pkin(2);
t354 = -pkin(5) * t136 + t135 * t218 - t350;
t341 = t135 * t225;
t353 = pkin(10) * t341 - t349;
t62 = -pkin(4) * t221 + qJD(5) - t73;
t351 = -t197 + t62;
t348 = pkin(5) * t341 + t347;
t343 = qJ(5) * t227;
t340 = t161 * t172;
t339 = t379 * t225;
t338 = t379 * t227;
t337 = t172 * t379;
t336 = t178 * t225;
t335 = t178 * t227;
t209 = pkin(3) * t229 + qJ(5);
t332 = t209 * t225;
t331 = t209 * t227;
t220 = pkin(12) + qJ(6);
t214 = sin(t220);
t330 = t214 * t217;
t215 = cos(t220);
t329 = t215 * t217;
t328 = t217 * t234;
t326 = t226 * t233;
t323 = -qJD(5) + t62;
t322 = qJDD(1) - g(3);
t122 = -pkin(4) * t259 - qJ(5) * t178 - t213;
t72 = t122 * t225 + t150 * t227;
t222 = t230 ^ 2;
t320 = -t233 ^ 2 + t222;
t317 = qJD(2) * t231;
t314 = qJD(6) * t228;
t313 = qJD(6) * t232;
t308 = pkin(10) * t339;
t6 = -t17 * t225 + t227 * t28;
t3 = pkin(5) * t95 - pkin(10) * t81 + t6;
t4 = -pkin(10) * t80 + t7;
t305 = -t228 * t4 + t232 * t3;
t301 = t226 * t317;
t300 = qJD(2) * t325;
t294 = t234 * t311;
t36 = -t225 * t65 + t227 * t93;
t47 = t118 * t227 - t225 * t73;
t71 = t122 * t227 - t150 * t225;
t212 = -pkin(3) * t365 - pkin(4);
t156 = pkin(5) * t339;
t284 = -t156 + t283;
t282 = -pkin(5) * t172 - pkin(10) * t338;
t279 = g(1) * t165 + g(2) * t163;
t278 = t228 * t3 + t232 * t4;
t24 = -pkin(5) * t379 + pkin(10) * t143 + t36;
t25 = pkin(10) * t145 + t37;
t273 = t228 * t25 - t232 * t24;
t10 = t228 * t24 + t232 * t25;
t52 = -pkin(5) * t259 - pkin(10) * t335 + t71;
t57 = -pkin(10) * t336 + t72;
t272 = -t228 * t57 + t232 * t52;
t271 = t228 * t52 + t232 * t57;
t168 = -t230 * t327 + t233 * t346;
t169 = t230 * t346 + t231 * t326;
t108 = t168 * t229 + t169 * t365;
t90 = -t108 * t225 - t227 * t325;
t91 = t108 * t227 - t225 * t325;
t270 = -t228 * t91 + t232 * t90;
t269 = t228 * t90 + t232 * t91;
t237 = qJD(2) ^ 2;
t266 = qJDD(2) * t234 - t231 * t237;
t173 = (-pkin(10) - t209) * t225;
t263 = -qJD(6) * t173 - t308 + t371;
t174 = t218 + t331;
t262 = qJD(6) * t174 + t282 - t372;
t261 = t168 * t365 - t169 * t229;
t186 = (-pkin(10) - qJ(5)) * t225;
t258 = -qJD(6) * t186 - t308 + t370;
t187 = t218 + t343;
t257 = qJD(5) * t225 + qJD(6) * t187 + t282 + t47;
t256 = -t172 * t37 - t225 * t376;
t21 = t143 * t314 + t145 * t313 - t228 * t80 + t232 * t81;
t252 = g(1) * t344 - g(2) * t345;
t182 = -t302 - t356;
t250 = -qJD(2) * t182 - t154 + t279;
t249 = t338 * t36 + t339 * t37 + t373;
t248 = t251 * t217;
t247 = t135 * t62 + t178 * t19 - t279;
t13 = pkin(5) * t80 + t19;
t51 = -pkin(5) * t145 + t62;
t246 = -t10 * t172 + t13 * t176 - t214 * t255 - t321 * t51;
t244 = t255 + t285;
t242 = -pkin(8) * qJDD(3) + (t182 + t302 - t356) * qJD(3);
t241 = t172 * t36 + t227 * t376;
t240 = t159 * t172 + t244;
t239 = t13 * t175 - t172 * t273 + t215 * t255 + t377 * t51;
t238 = -t159 * t379 - t243 + t254;
t210 = -pkin(4) - t357;
t185 = t212 - t357;
t152 = t157 * pkin(4);
t139 = -qJD(3) * t169 - t230 * t300;
t138 = qJD(3) * t168 + t233 * t300;
t121 = t126 * pkin(4);
t120 = t124 * pkin(4);
t112 = t175 * t178;
t111 = t176 * t178;
t101 = pkin(5) * t336 - t260;
t96 = t172 ^ 2 - t379 ^ 2;
t92 = qJDD(6) + t95;
t69 = -t274 + (-qJD(2) * t178 - t172) * t221;
t68 = t94 - t383;
t56 = t156 + t74;
t50 = qJD(4) * t108 + t138 * t229 - t139 * t365;
t49 = qJD(4) * t261 + t138 * t365 + t229 * t139;
t44 = t135 * t176 + t313 * t335 - t314 * t336;
t43 = -t135 * t175 - t167 * t178;
t40 = t225 * t301 + t227 * t49;
t39 = -t225 * t49 + t227 * t301;
t15 = -t161 * t377 + t172 * t374 - t175 * t92;
t14 = -t161 * t321 - t172 * t267 + t176 * t92;
t8 = t176 * t21 + t267 * t321;
t1 = -t175 * t21 - t176 * t22 + t267 * t377 - t321 * t374;
t2 = [t322, 0, t266 * t226 (-qJDD(2) * t231 - t234 * t237) * t226, 0, 0, 0, 0, 0, qJD(3) * t139 + qJDD(3) * t168 + (-t230 * t294 + t233 * t266) * t226, -qJD(3) * t138 - qJDD(3) * t169 + (-t230 * t266 - t233 * t294) * t226, 0, 0, 0, 0, 0, t261 * t219 - t221 * t50 + (-t234 * t95 - t317 * t379) * t226, -t108 * t219 - t221 * t49 + (-t172 * t317 - t234 * t94) * t226, -t145 * t50 - t261 * t80 - t379 * t39 + t90 * t95, -t143 * t50 - t261 * t81 + t379 * t40 - t91 * t95, t143 * t39 + t145 * t40 - t80 * t91 - t81 * t90, -t19 * t261 + t36 * t39 + t37 * t40 + t50 * t62 + t6 * t90 + t7 * t91 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t269 - t228 * t40 + t232 * t39) * t161 + t270 * t92 - t50 * t374 - t261 * t22 -(qJD(6) * t270 + t228 * t39 + t232 * t40) * t161 - t269 * t92 - t50 * t267 - t261 * t21; 0, qJDD(2), t322 * t325 + t280, -t322 * t327 + t279, qJDD(2) * t222 + 0.2e1 * t233 * t295, 0.2e1 * t230 * t309 - 0.2e1 * t311 * t320, qJDD(3) * t230 + t233 * t236, qJDD(3) * t233 - t230 * t236, 0, t230 * t242 + t233 * t367, -t230 * t367 + t233 * t242, -t135 * t172 + t178 * t94, t135 * t379 + t136 * t172 - t178 * t95 + t259 * t94, t135 * t221 + t178 * t219, -t136 * t221 + t219 * t259, 0, -t123 * t259 + t136 * t159 - t213 * t95 + t219 * t260 - t221 * t347 - t253 * t379 - t248, t123 * t178 + t135 * t159 - t150 * t219 - t253 * t172 - t213 * t94 + t221 * t381 + t369, t36 * t136 - t260 * t80 - t6 * t259 + t71 * t95 - t227 * t248 + (-g(3) * t327 + t247) * t225 - t350 * t379 - t347 * t145, -t37 * t136 - t260 * t81 + t7 * t259 - t72 * t95 - t280 * t225 * t217 + t247 * t227 - (-t225 * t328 + t227 * t231) * t358 + t349 * t379 - t347 * t143, -t71 * t81 - t72 * t80 + (-t225 * t7 - t227 * t6) * t178 + t350 * t143 + t349 * t145 + (-t225 * t37 - t227 * t36) * t135 - t369, -t19 * t260 + t347 * t62 + t349 * t37 + t350 * t36 + t6 * t71 + t7 * t72 - (t231 * t358 + t279) * t366 + (-t234 * t358 + t280) * (pkin(4) * t217 + qJ(5) * t216 + t213) -t112 * t21 - t267 * t43, -t111 * t21 + t112 * t22 + t267 * t44 + t374 * t43, -t112 * t92 - t136 * t267 + t161 * t43 - t21 * t259, -t111 * t92 + t136 * t374 - t161 * t44 + t22 * t259, t136 * t161 - t259 * t92, t272 * t92 - t305 * t259 - t273 * t136 + t101 * t22 + t13 * t111 + t51 * t44 - g(1) * (-t164 * t329 + t165 * t214) - g(2) * (-t162 * t329 + t163 * t214) - t348 * t374 - (t214 * t231 + t215 * t328) * t358 + (t228 * t353 - t232 * t354) * t161 + (t10 * t259 - t161 * t271) * qJD(6), -t271 * t92 + t278 * t259 - t10 * t136 + t101 * t21 - t13 * t112 + t51 * t43 - g(1) * (t164 * t330 + t165 * t215) - g(2) * (t162 * t330 + t163 * t215) - t348 * t267 - (-t214 * t328 + t215 * t231) * t358 + (t228 * t354 + t232 * t353) * t161 + (-t161 * t272 - t259 * t273) * qJD(6); 0, 0, 0, 0, -t230 * t237 * t233, t320 * t237, t310, t309, qJDD(3), -g(3) * t168 + t230 * t250 - t252 * t326 + t196, g(3) * t169 + (t226 * t252 - t286) * t230 + t250 * t233, t337, t96, t68, t69, t219, t78 * t221 + (t219 * t365 - t221 * t316 + t318 * t379) * pkin(3) + t240, t79 * t221 + (t172 * t318 - t219 * t229 - t221 * t297) * pkin(3) + t238, -t95 * t332 + t212 * t80 - t283 * t145 - (t225 * t351 - t45) * t379 + t241, -t95 * t331 + t212 * t81 - t283 * t143 - (t227 * t351 + t46) * t379 + t256, -t80 * t331 - t143 * t45 - t371 * t145 + (-t143 * t197 + t209 * t81 - t6) * t225 + t249, t7 * t331 - t6 * t332 + t19 * t212 - g(1) * (t127 * qJ(5) - t121 + (-t165 * t230 + t233 * t290) * pkin(3)) - g(2) * (t125 * qJ(5) - t120 + (-t163 * t230 - t233 * t291) * pkin(3)) - g(3) * (pkin(3) * t168 + qJ(5) * t158 - t152) + t283 * t62 - t371 * t37 + t372 * t36, t8, t1, t14, t15, t340 (t173 * t232 - t174 * t228) * t92 + t185 * t22 - t284 * t374 + (t228 * t263 - t232 * t262) * t161 + t239 -(t173 * t228 + t174 * t232) * t92 + t185 * t21 - t284 * t267 + (t228 * t262 + t232 * t263) * t161 + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, t96, t68, t69, t219, t221 * t74 + t240, t221 * t73 + t238, -qJ(5) * t225 * t95 - pkin(4) * t80 + t145 * t74 - (t225 * t323 - t47) * t379 + t241, -t95 * t343 - pkin(4) * t81 + t143 * t74 - (t227 * t323 + t48) * t379 + t256, -t80 * t343 - t143 * t47 - t370 * t145 + (qJ(5) * t81 - qJD(5) * t143 - t6) * t225 + t249, -t19 * pkin(4) + g(1) * t121 + g(2) * t120 + g(3) * t152 - t36 * t47 - t37 * t48 - t62 * t74 + (-t225 * t36 + t227 * t37) * qJD(5) + (-t225 * t6 + t373) * qJ(5), t8, t1, t14, t15, t340 (t186 * t232 - t187 * t228) * t92 + t210 * t22 + t56 * t374 + (t228 * t258 - t232 * t257) * t161 + t239 -(t186 * t228 + t187 * t232) * t92 + t210 * t21 + t56 * t267 + (t228 * t257 + t232 * t258) * t161 + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143 * t379 + t80, -t145 * t379 + t81, -t143 ^ 2 - t145 ^ 2, -t143 * t36 - t145 * t37 - t244 + t368, 0, 0, 0, 0, 0, t22 - t375, t21 + t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267 * t374, t267 ^ 2 - t374 ^ 2, t21 - t382, -t22 - t375, t92, t51 * t267 - g(1) * (-t127 * t214 + t164 * t215) - g(2) * (-t125 * t214 + t162 * t215) - g(3) * (-t158 * t214 - t215 * t325) + t305 + t378 * t10, -t51 * t374 - g(1) * (-t127 * t215 - t164 * t214) - g(2) * (-t125 * t215 - t162 * t214) - g(3) * (-t158 * t215 + t214 * t325) - t278 - t378 * t273;];
tau_reg  = t2;
