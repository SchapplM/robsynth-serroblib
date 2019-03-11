% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:56
% EndTime: 2019-03-09 02:34:07
% DurationCPUTime: 6.81s
% Computational Cost: add. (12178->522), mult. (24850->629), div. (0->0), fcn. (18471->14), ass. (0->271)
t204 = cos(qJ(6));
t283 = qJD(6) * t204;
t196 = sin(pkin(10));
t197 = cos(pkin(10));
t331 = cos(qJ(4));
t263 = qJD(1) * t331;
t202 = sin(qJ(4));
t287 = qJD(1) * t202;
t130 = -t196 * t263 - t197 * t287;
t242 = t197 * t263;
t265 = t196 * t287;
t131 = t242 - t265;
t201 = sin(qJ(5));
t330 = cos(qJ(5));
t86 = -t330 * t130 + t131 * t201;
t365 = t204 * t86;
t390 = t283 + t365;
t199 = -pkin(1) - qJ(3);
t338 = -qJD(1) * qJD(3) + qJDD(1) * t199;
t143 = qJDD(2) + t338;
t250 = -pkin(7) * qJDD(1) + t143;
t106 = t250 * t196;
t107 = t250 * t197;
t249 = -t202 * t106 + t331 * t107;
t156 = t199 * qJD(1) + qJD(2);
t257 = -pkin(7) * qJD(1) + t156;
t122 = t257 * t196;
t123 = t257 * t197;
t83 = t331 * t122 + t202 * t123;
t52 = -t83 * qJD(4) + t249;
t266 = t331 * t197;
t281 = t196 * qJDD(1);
t231 = qJDD(1) * t266 - t202 * t281;
t350 = t331 * t196 + t202 * t197;
t345 = qJD(4) * t350;
t93 = qJD(1) * t345 - t231;
t37 = qJDD(4) * pkin(4) + t93 * pkin(8) + t52;
t262 = qJD(4) * t331;
t268 = -t331 * t106 - t202 * t107 - t123 * t262;
t286 = qJD(4) * t202;
t51 = -t122 * t286 - t268;
t213 = -qJD(4) * t265 + qJDD(1) * t350;
t94 = qJD(4) * t242 + t213;
t38 = -pkin(8) * t94 + t51;
t259 = t201 * t38 - t330 * t37;
t69 = t130 * pkin(8) + t83;
t272 = t330 * t69;
t304 = t122 * t202;
t82 = t331 * t123 - t304;
t68 = -pkin(8) * t131 + t82;
t67 = qJD(4) * pkin(4) + t68;
t40 = t201 * t67 + t272;
t10 = -t40 * qJD(5) - t259;
t278 = qJDD(4) + qJDD(5);
t8 = -t278 * pkin(5) - t10;
t140 = -t202 * t196 + t266;
t353 = t201 * t350;
t95 = -t330 * t140 + t353;
t389 = t8 * t95;
t261 = qJD(5) * t330;
t285 = qJD(5) * t201;
t47 = -t130 * t261 + t131 * t285 + t201 * t94 + t330 * t93;
t388 = t47 * t95;
t200 = sin(qJ(6));
t223 = t201 * t130 + t330 * t131;
t279 = qJD(4) + qJD(5);
t246 = t204 * t279;
t284 = qJD(6) * t200;
t27 = -qJD(6) * t246 - t200 * t278 + t204 * t47 + t223 * t284;
t387 = t95 * t27;
t75 = t200 * t279 + t204 * t223;
t306 = qJD(6) * t75;
t28 = -t200 * t47 - t204 * t278 + t306;
t386 = t95 * t28;
t258 = t201 * t93 - t330 * t94;
t344 = t223 * qJD(5);
t48 = -t258 + t344;
t193 = qJDD(1) * qJ(2);
t194 = qJD(1) * qJD(2);
t348 = t193 + t194;
t154 = qJDD(3) + t348;
t175 = pkin(3) * t281;
t134 = t175 + t154;
t332 = pkin(4) * t94;
t72 = t134 + t332;
t12 = pkin(5) * t48 + pkin(9) * t47 + t72;
t11 = t204 * t12;
t35 = t279 * pkin(9) + t40;
t179 = qJD(1) * qJ(2) + qJD(3);
t183 = t196 * pkin(3);
t147 = qJD(1) * t183 + t179;
t101 = -pkin(4) * t130 + t147;
t49 = pkin(5) * t86 - pkin(9) * t223 + t101;
t14 = t200 * t49 + t204 * t35;
t9 = t201 * t37 + t67 * t261 - t69 * t285 + t330 * t38;
t7 = t278 * pkin(9) + t9;
t3 = -qJD(6) * t14 - t200 * t7 + t11;
t364 = qJD(6) + t86;
t325 = t14 * t364;
t385 = t3 + t325;
t384 = t278 * t95;
t379 = t200 * t364;
t383 = t75 * t379;
t205 = cos(qJ(1));
t187 = g(2) * t205;
t203 = sin(qJ(1));
t349 = g(1) * t203 - t187;
t232 = t196 * t286 - t197 * t262;
t358 = -qJD(5) * t353 + t140 * t261 - t201 * t345 - t330 * t232;
t233 = t330 * t350;
t363 = t201 * t140 + t233;
t382 = -t10 * t95 + t358 * t40 + t363 * t9 - t349;
t26 = t28 * t204;
t73 = t200 * t223 - t246;
t381 = t379 * t73 - t26;
t13 = -t200 * t35 + t204 * t49;
t326 = t13 * t364;
t356 = t279 * t86;
t380 = -t47 + t356;
t320 = t86 ^ 2;
t321 = t223 ^ 2;
t378 = -t320 + t321;
t317 = -t200 * t28 - t73 * t283;
t377 = -t204 * t27 - t365 * t73 + t317;
t24 = t27 * t200;
t376 = t390 * t75 - t24;
t322 = t75 * t223;
t46 = qJDD(6) + t48;
t43 = t200 * t46;
t375 = t364 * t390 - t322 + t43;
t374 = -t13 * t358 - t3 * t363;
t2 = qJD(6) * t13 + t200 * t12 + t204 * t7;
t373 = t14 * t358 + t2 * t363;
t372 = -t27 * t363 + t358 * t75;
t371 = -t28 * t363 - t358 * t73;
t370 = t358 * t86 + t363 * t48;
t369 = -t278 * t363 - t279 * t358;
t310 = t201 * t69;
t39 = t330 * t67 - t310;
t34 = -t279 * pkin(5) - t39;
t367 = t34 * t86;
t192 = pkin(10) + qJ(4);
t182 = qJ(5) + t192;
t174 = cos(t182);
t301 = t174 * t203;
t274 = g(1) * t301;
t366 = t8 + t274;
t319 = t86 * t223;
t173 = sin(t182);
t227 = g(3) * t173 + t174 * t187;
t361 = t140 * qJD(3);
t57 = pkin(5) * t223 + pkin(9) * t86;
t151 = t173 * t187;
t165 = g(3) * t174;
t302 = t173 * t203;
t267 = -g(1) * t302 + t151 - t165;
t360 = t101 * t86 - t267 - t9;
t359 = t52 * t140 - t232 * t83 - t345 * t82 + t350 * t51 - t349;
t324 = t73 * t223;
t355 = t364 * t223;
t190 = t196 ^ 2;
t191 = t197 ^ 2;
t288 = t190 + t191;
t354 = t156 * t288;
t236 = g(1) * t205 + g(2) * t203;
t352 = t236 * t174;
t318 = -pkin(7) + t199;
t144 = t318 * t196;
t145 = t318 * t197;
t100 = t331 * t144 + t202 * t145;
t163 = t174 * pkin(9);
t351 = -pkin(5) * t173 + t163;
t220 = -t236 + t154;
t294 = t223 * qJD(4);
t347 = t294 + t258;
t346 = -t13 * t200 + t14 * t204;
t30 = t34 * t284;
t343 = -t13 * t223 + t227 * t204 + t30;
t31 = t34 * t283;
t342 = t14 * t223 + t366 * t200 + t31;
t341 = -t101 * t223 + t227 - t259 - t274;
t340 = -t201 * t232 + t330 * t345;
t337 = t130 * t232 + t350 * t94;
t336 = qJD(4) * t232 - qJDD(4) * t350;
t334 = t131 ^ 2;
t333 = 0.2e1 * t194;
t1 = t2 * t204;
t323 = t75 * t73;
t198 = -pkin(7) - qJ(3);
t312 = t200 * t75;
t44 = t204 * t46;
t307 = pkin(1) * qJDD(1);
t305 = qJD(6) * t364;
t303 = t131 * t130;
t300 = t200 * t203;
t299 = t200 * t205;
t296 = t203 * t204;
t295 = t204 * t205;
t170 = qJ(2) + t183;
t292 = -qJD(4) ^ 2 * t350 + t140 * qJDD(4);
t290 = t205 * pkin(1) + t203 * qJ(2);
t280 = t197 * qJDD(1);
t275 = pkin(9) * t305;
t177 = pkin(4) * t201 + pkin(9);
t271 = t177 * t305;
t270 = t95 * t284;
t269 = t95 * t283;
t264 = g(1) * (pkin(5) * t301 + pkin(9) * t302);
t185 = t205 * qJ(2);
t260 = -t203 * pkin(1) + t185;
t251 = t288 * t143;
t99 = -t144 * t202 + t331 * t145;
t248 = t1 + t267;
t247 = qJD(6) * t363 + qJD(1);
t245 = qJDD(2) - t307;
t244 = pkin(4) * t261;
t41 = t201 * t68 + t272;
t241 = pkin(4) * t285 - t41;
t240 = -pkin(9) * t46 + t367;
t59 = qJD(5) * t233 + t140 * t285 + t340;
t239 = -t34 * t59 - t389;
t234 = -t364 * t59 - t46 * t95;
t230 = t13 * t204 + t14 * t200;
t79 = -pkin(8) * t140 + t99;
t80 = -pkin(8) * t350 + t100;
t54 = t201 * t79 + t330 * t80;
t109 = pkin(4) * t350 + t170;
t55 = pkin(5) * t363 + pkin(9) * t95 + t109;
t22 = t200 * t55 + t204 * t54;
t21 = -t200 * t54 + t204 * t55;
t229 = -t284 * t364 - t379 * t86 + t44;
t180 = sin(t192);
t146 = pkin(4) * t180 + t183;
t189 = -pkin(8) + t198;
t228 = t205 * t146 + t203 * t189 + t260;
t226 = -qJD(6) * t49 + t165 - t7;
t225 = t203 * t146 - t205 * t189 + t290;
t112 = -t232 * pkin(4) + qJD(2);
t217 = t350 * qJD(1);
t215 = t175 + t220;
t181 = cos(t192);
t214 = g(3) * t180 - t181 * t349;
t212 = -t177 * t46 - t244 * t364 + t367;
t211 = t220 + t348;
t210 = -t230 * qJD(6) - t3 * t200 + t1;
t209 = -t131 * t345 - t140 * t93;
t76 = -qJD(3) * t350 - t144 * t286 + t145 * t262;
t207 = pkin(8) * t345 - t144 * t262 - t145 * t286 - t361;
t206 = qJD(1) ^ 2;
t178 = -t330 * pkin(4) - pkin(5);
t126 = t130 ^ 2;
t121 = t173 * t295 - t300;
t120 = t173 * t299 + t296;
t119 = t173 * t296 + t299;
t118 = -t173 * t300 + t295;
t77 = -qJD(4) * t100 - t361;
t65 = t232 * pkin(8) + t76;
t60 = qJD(5) * t363 + t340;
t53 = t201 * t80 - t330 * t79;
t50 = pkin(4) * t131 + t57;
t42 = t330 * t68 - t310;
t29 = pkin(5) * t358 + t59 * pkin(9) + t112;
t20 = t200 * t57 + t204 * t39;
t19 = -t200 * t39 + t204 * t57;
t18 = t200 * t50 + t204 * t42;
t17 = -t200 * t42 + t204 * t50;
t16 = t54 * qJD(5) + t201 * t65 - t330 * t207;
t15 = t201 * t207 + t79 * t261 - t80 * t285 + t330 * t65;
t5 = -t22 * qJD(6) - t200 * t15 + t204 * t29;
t4 = t21 * qJD(6) + t204 * t15 + t200 * t29;
t6 = [0, 0, 0, 0, 0, qJDD(1), t349, t236, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t349 - 0.2e1 * t307, 0.2e1 * t193 + t333 - t236, -t245 * pkin(1) - g(1) * t260 - g(2) * t290 + (t193 + t333) * qJ(2), t191 * qJDD(1), -0.2e1 * t196 * t280, 0, t190 * qJDD(1), 0, 0, t211 * t196, t211 * t197, t349 + t288 * (-t143 - t338) t154 * qJ(2) + t179 * qJD(2) - g(1) * (t199 * t203 + t185) - g(2) * (qJ(3) * t205 + t290) + t199 * t251 - qJD(3) * t354, t209, -t130 * t345 + t131 * t232 - t140 * t94 + t350 * t93, t292, t337, t336, 0, -qJD(2) * t130 + t77 * qJD(4) + t99 * qJDD(4) + t134 * t350 - t147 * t232 + t170 * t94 - t180 * t236, qJD(2) * t131 - t76 * qJD(4) - t100 * qJDD(4) + t134 * t140 - t147 * t345 - t170 * t93 - t181 * t236, -t100 * t94 + t76 * t130 - t77 * t131 + t99 * t93 - t359, t51 * t100 + t83 * t76 + t52 * t99 + t82 * t77 + t134 * t170 + t147 * qJD(2) - g(1) * (t205 * t183 + t185 + (-pkin(1) + t198) * t203) - g(2) * (t203 * t183 - t198 * t205 + t290) -t223 * t59 + t388, -t223 * t358 + t363 * t47 + t48 * t95 + t59 * t86, -t279 * t59 - t384, t370, t369, 0, t101 * t358 + t109 * t48 + t112 * t86 - t16 * t279 - t173 * t236 - t278 * t53 + t363 * t72, -t101 * t59 - t109 * t47 + t112 * t223 - t15 * t279 - t278 * t54 - t72 * t95 - t352, -t15 * t86 + t16 * t223 + t39 * t59 - t47 * t53 - t48 * t54 - t382, -g(1) * t228 - g(2) * t225 - t10 * t53 + t101 * t112 + t72 * t109 + t40 * t15 - t39 * t16 + t9 * t54, t75 * t270 + (-t59 * t75 + t387) * t204 (t204 * t73 + t312) * t59 - (t24 - t26 + (t200 * t73 - t204 * t75) * qJD(6)) * t95, t204 * t234 + t270 * t364 + t372, -t73 * t269 + (-t59 * t73 - t386) * t200, -t200 * t234 + t269 * t364 + t371, t358 * t364 + t363 * t46, -g(1) * t121 - g(2) * t119 + t16 * t73 + t200 * t239 + t21 * t46 + t53 * t28 - t31 * t95 + t364 * t5 - t374, g(1) * t120 - g(2) * t118 + t16 * t75 + t204 * t239 - t22 * t46 - t53 * t27 + t30 * t95 - t364 * t4 - t373, t21 * t27 - t22 * t28 - t4 * t73 - t5 * t75 + t230 * t59 + t352 - (-qJD(6) * t346 - t2 * t200 - t3 * t204) * t95, t2 * t22 + t14 * t4 + t3 * t21 + t13 * t5 + t8 * t53 + t34 * t16 - g(1) * (-t205 * t351 + t228) - g(2) * (-t203 * t351 + t225); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t206, -qJ(2) * t206 + t245 - t349, 0, 0, 0, 0, 0, 0, -t206 * t196, -t206 * t197, -t288 * qJDD(1), -t179 * qJD(1) + t251 - t349, 0, 0, 0, 0, 0, 0, qJD(1) * t130 + t292, -qJD(1) * t131 + t336, -t209 - t337, -t147 * qJD(1) + t359, 0, 0, 0, 0, 0, 0, -qJD(1) * t86 - t279 * t60 - t384, -qJD(1) * t223 + t369, t223 * t60 - t370 - t388, -qJD(1) * t101 - t39 * t60 + t382, 0, 0, 0, 0, 0, 0, -t363 * t43 + t386 + t60 * t73 + (-t200 * t358 - t204 * t247) * t364, -t363 * t44 - t387 + t60 * t75 + (t200 * t247 - t204 * t358) * t364 (t247 * t75 + t371) * t204 + (t247 * t73 + t372) * t200, t34 * t60 + t389 + (-t13 * t247 + t373) * t204 + (-t14 * t247 + t374) * t200 - t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, t280, -t288 * t206, qJD(1) * t354 + t220, 0, 0, 0, 0, 0, 0 (t242 + t131) * qJD(4) + t213 (t130 - t217) * qJD(4) + t231, -t126 - t334, -t130 * t83 + t131 * t82 + t215, 0, 0, 0, 0, 0, 0, -t258 + t294 + 0.2e1 * t344, -t47 - t356, -t320 - t321, t223 * t39 + t40 * t86 + t215 + t332, 0, 0, 0, 0, 0, 0, t229 - t324, -t204 * t364 ^ 2 - t322 - t43 (-t73 * t86 + t27) * t204 + t383 + t317, -t34 * t223 + t385 * t204 + (t2 - t326) * t200 - t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t303, -t126 + t334 (-t130 - t217) * qJD(4) + t231, t303 (-t242 + t131) * qJD(4) - t213, qJDD(4), -t147 * t131 + t214 + t249, g(3) * t181 - t147 * t130 + t349 * t180 + (t82 + t304) * qJD(4) + t268, 0, 0, t319, t378, t380, -t319, t347, t278, t41 * qJD(4) + (t41 - t40) * qJD(5) + (-t131 * t86 + t330 * t278 - t279 * t285) * pkin(4) + t341, t42 * t279 + (-t131 * t223 - t201 * t278 - t261 * t279) * pkin(4) + t360, -t39 * t86 + t40 * t223 - t41 * t223 + t42 * t86 + (t330 * t47 - t201 * t48 + (t201 * t223 - t330 * t86) * qJD(5)) * pkin(4), t39 * t41 - t40 * t42 + (t330 * t10 - t101 * t131 + t201 * t9 + (-t201 * t39 + t330 * t40) * qJD(5) + t214) * pkin(4), t376, -t312 * t364 + t377, t375, t381, t229 + t324, -t355, -t17 * t364 + t178 * t28 + t241 * t73 + (-t366 - t271) * t204 + t212 * t200 + t343, -t178 * t27 + t18 * t364 + t241 * t75 + t212 * t204 + (-t227 + t271) * t200 + t342, t17 * t75 + t18 * t73 + (-t73 * t244 - t13 * t86 - t177 * t28 + (t177 * t75 - t13) * qJD(6)) * t204 + (t75 * t244 - t14 * t86 - t177 * t27 - t3 + (t177 * t73 - t14) * qJD(6)) * t200 + t248, t8 * t178 - t14 * t18 - t13 * t17 - t34 * t41 - t264 - g(3) * t351 - (-pkin(5) * t174 - pkin(9) * t173) * t187 + ((t201 * t34 + t330 * t346) * qJD(5) + t214) * pkin(4) + t210 * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t319, t378, t380, -t319, t347, t278, t40 * qJD(4) + t341, t279 * t39 + t360, 0, 0, t376, t377 - t383, t375, t381, -t364 * t379 + t324 + t44, -t355, -pkin(5) * t28 - t19 * t364 - t40 * t73 + t240 * t200 + (-t366 - t275) * t204 + t343, pkin(5) * t27 + t20 * t364 - t40 * t75 + t240 * t204 + (-t227 + t275) * t200 + t342, t19 * t75 + t20 * t73 + (-t326 + (-t28 + t306) * pkin(9)) * t204 + ((qJD(6) * t73 - t27) * pkin(9) - t385) * t200 + t248, -t14 * t20 - t13 * t19 - t34 * t40 - t264 - g(3) * t163 + (t227 - t8) * pkin(5) + (t210 + t151) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323, -t73 ^ 2 + t75 ^ 2, t364 * t73 - t27, -t323, t364 * t75 - t28, t46, -g(1) * t118 - g(2) * t120 + t200 * t226 - t283 * t35 - t34 * t75 + t11 + t325, g(1) * t119 - g(2) * t121 + t326 + t34 * t73 + (qJD(6) * t35 - t12) * t200 + t226 * t204, 0, 0;];
tau_reg  = t6;
