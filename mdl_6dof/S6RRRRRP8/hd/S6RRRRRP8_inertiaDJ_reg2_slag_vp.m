% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRRP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:50
% EndTime: 2019-03-10 01:56:20
% DurationCPUTime: 11.82s
% Computational Cost: add. (15968->587), mult. (41165->975), div. (0->0), fcn. (40686->10), ass. (0->276)
t194 = sin(qJ(5));
t191 = t194 ^ 2;
t198 = cos(qJ(5));
t192 = t198 ^ 2;
t389 = t191 + t192;
t195 = sin(qJ(4));
t196 = sin(qJ(3));
t199 = cos(qJ(3));
t193 = sin(pkin(6));
t197 = sin(qJ(2));
t331 = t193 * t197;
t344 = cos(pkin(6));
t236 = t196 * t344 + t199 * t331;
t278 = t344 * t199;
t237 = t196 * t331 - t278;
t368 = cos(qJ(4));
t103 = t195 * t236 + t237 * t368;
t281 = qJD(4) * t368;
t275 = pkin(3) * t281;
t388 = t389 * t275;
t276 = t344 * qJD(3);
t200 = cos(qJ(2));
t320 = qJD(2) * t200;
t285 = t193 * t320;
t387 = t276 + t285;
t380 = (t191 - t192) * qJD(5);
t290 = t368 * t199;
t271 = qJD(3) * t290;
t329 = t195 * t196;
t379 = qJD(3) + qJD(4);
t114 = -t199 * t281 + t329 * t379 - t271;
t291 = t368 * t196;
t328 = t195 * t199;
t155 = t291 + t328;
t104 = -t195 * t237 + t236 * t368;
t330 = t193 * t200;
t93 = t104 * t194 + t198 * t330;
t298 = t194 * t330;
t94 = t104 * t198 - t298;
t258 = t194 * t94 + t198 * t93;
t347 = t198 * t94;
t321 = qJD(2) * t197;
t176 = t193 * t321;
t189 = qJD(5) * t198;
t288 = qJD(3) * t331;
t296 = t196 * t387 + t199 * t288;
t162 = t196 * t288;
t377 = -t199 * t387 + t162;
t203 = qJD(4) * t103 + t195 * t296 + t368 * t377;
t43 = -qJD(5) * t298 + t104 * t189 - t176 * t198 - t194 * t203;
t348 = t198 * t43;
t349 = t194 * t93;
t42 = qJD(5) * t93 - t176 * t194 + t198 * t203;
t352 = t194 * t42;
t386 = t155 * (qJD(5) * (-t347 + t349) - t348 + t352) + t258 * t114;
t314 = qJD(6) * t103;
t64 = qJD(4) * t104 - t195 * t377 + t368 * t296;
t356 = qJ(6) * t64;
t188 = qJD(5) * t194;
t277 = qJD(2) * t344;
t269 = t197 * t277;
t272 = t296 * pkin(3);
t201 = pkin(1) * t269 + t64 * pkin(4) + pkin(8) * t285 + pkin(11) * t203 + t272;
t295 = pkin(1) * t344;
t229 = pkin(8) * t331 - t200 * t295;
t129 = -pkin(2) * t344 + t229;
t107 = pkin(3) * t237 + t129;
t206 = t103 * pkin(4) - t104 * pkin(11) + t107;
t205 = t198 * t206;
t145 = pkin(8) * t330 + t197 * t295;
t130 = pkin(9) * t344 + t145;
t135 = t229 * qJD(2);
t245 = qJD(2) * (pkin(2) * t197 - pkin(9) * t200);
t249 = -pkin(2) * t200 - pkin(9) * t197 - pkin(1);
t318 = qJD(3) * t199;
t319 = qJD(3) * t196;
t66 = (t199 * t245 - t249 * t319) * t193 - t130 * t318 + t196 * t135;
t202 = pkin(3) * t176 + pkin(10) * t377 + t66;
t242 = t249 * t193;
t122 = t199 * t242;
t65 = -t193 * t196 * t245 - qJD(3) * t122 + t130 * t319 + t135 * t199;
t213 = -pkin(10) * t296 - t65;
t316 = qJD(4) * t195;
t100 = -t196 * t130 + t122;
t77 = -pkin(3) * t330 - pkin(10) * t236 + t100;
t101 = t199 * t130 + t196 * t242;
t86 = -pkin(10) * t237 + t101;
t17 = -t195 * t202 - t213 * t368 - t281 * t77 + t316 * t86;
t217 = pkin(11) * t176 - t17;
t53 = t195 * t77 + t368 * t86;
t50 = -pkin(11) * t330 + t53;
t6 = -qJD(5) * t205 + t188 * t50 - t194 * t201 - t198 * t217;
t2 = t314 - t6 + t356;
t369 = t64 * pkin(5);
t358 = t194 * t206 + t198 * t50;
t7 = -qJD(5) * t358 - t194 * t217 + t198 * t201;
t4 = -t369 - t7;
t385 = t194 * t4 + t2 * t198;
t44 = t103 * t188 - t198 * t64;
t384 = pkin(11) * t44;
t52 = -t195 * t86 + t368 * t77;
t49 = pkin(4) * t330 - t52;
t29 = pkin(5) * t93 - qJ(6) * t94 + t49;
t18 = -t195 * t213 + t202 * t368 - t281 * t86 - t316 * t77;
t16 = -pkin(4) * t176 - t18;
t8 = pkin(5) * t43 + qJ(6) * t42 - qJD(6) * t94 + t16;
t282 = t188 * t29 - t8 * t198;
t280 = -t16 * t198 + t188 * t49;
t370 = -pkin(10) - pkin(9);
t169 = t370 * t199;
t118 = -t169 * t368 + t329 * t370;
t154 = -t290 + t329;
t186 = -pkin(3) * t199 - pkin(2);
t235 = -pkin(4) * t154 + pkin(11) * t155 - t186;
t382 = t198 * t118 - t194 * t235;
t222 = qJD(3) * t237;
t381 = t196 * t222 - t199 * t296;
t115 = t379 * t155;
t334 = t155 * t194;
t341 = t114 * t194;
t98 = t155 * t189 - t341;
t378 = t103 * t98 + t115 * t93 + t154 * t43 + t334 * t64;
t305 = pkin(3) * t319;
t228 = pkin(4) * t115 + pkin(11) * t114 + t305;
t156 = t370 * t291;
t294 = qJD(3) * t370;
t157 = t196 * t294;
t91 = -qJD(4) * t156 - t157 * t368 - t169 * t316 - t294 * t328;
t37 = -qJD(5) * t382 + t194 * t91 + t198 * t228;
t346 = t42 * t198;
t351 = t194 * t43;
t12 = qJD(5) * t258 + t346 + t351;
t266 = pkin(5) * t198 + qJ(6) * t194;
t375 = qJD(5) * t266 - t198 * qJD(6);
t374 = t199 ^ 2;
t373 = -0.2e1 * t380;
t372 = 0.2e1 * t193;
t371 = 0.2e1 * qJD(6);
t367 = pkin(3) * t193;
t366 = pkin(3) * t195;
t365 = pkin(9) * t193;
t364 = pkin(11) * t115;
t363 = pkin(11) * t154;
t362 = t115 * pkin(5);
t5 = t6 * t198;
t360 = t16 * t194 + t189 * t49;
t184 = pkin(11) + t366;
t251 = t198 * t275;
t359 = -t184 * t348 - t251 * t93;
t357 = pkin(3) * qJD(4);
t117 = -t195 * t169 - t156;
t92 = qJD(4) * t118 + t195 * t157 - t271 * t370;
t355 = t117 * t92;
t353 = t184 * t93;
t350 = t194 * t64;
t345 = t117 * t189 + t194 * t92;
t343 = qJ(6) * t115;
t265 = pkin(5) * t194 - qJ(6) * t198;
t90 = t155 * t265 + t117;
t342 = qJD(5) * t90;
t340 = t114 * t198;
t339 = t115 * t184;
t338 = t117 * t195;
t136 = t145 * qJD(2);
t337 = t136 * t196;
t336 = t154 * t184;
t335 = t155 * t114;
t333 = t155 * t198;
t332 = t192 * t114;
t327 = t198 * t115;
t311 = t194 * qJD(6);
t140 = pkin(5) * t188 - qJ(6) * t189 - t311;
t301 = pkin(3) * t316;
t123 = t140 + t301;
t326 = -t123 - t140;
t325 = t388 * t184;
t324 = t388 * pkin(11);
t306 = t368 * pkin(3);
t185 = -t306 - pkin(4);
t323 = t185 * t189 + t194 * t301;
t317 = qJD(3) * t200;
t313 = qJD(6) * t154;
t312 = t129 * qJD(3);
t309 = 0.2e1 * t93 * t43;
t24 = -t194 * t50 + t205;
t23 = -t103 * pkin(5) - t24;
t308 = t189 * t23 + t385;
t307 = -0.2e1 * pkin(2) * qJD(3);
t51 = 0.2e1 * t103 * t64;
t102 = 0.2e1 * t154 * t115;
t304 = pkin(4) * t188;
t303 = pkin(4) * t189;
t302 = pkin(9) * t318;
t300 = pkin(11) * t188;
t299 = pkin(11) * t189;
t297 = t194 * t340;
t85 = t90 * t188;
t293 = t194 * t368;
t292 = t198 * t368;
t190 = t193 ^ 2;
t289 = t190 * t320;
t287 = t196 * t317;
t286 = t184 * t189;
t284 = t194 * t189;
t283 = t196 * t318;
t151 = t155 ^ 2;
t274 = t151 * t284;
t273 = t197 * t289;
t268 = t42 * t93 - t43 * t94;
t267 = t103 * t43 + t64 * t93;
t22 = qJ(6) * t103 + t358;
t264 = t194 * t22 - t198 * t23;
t263 = t194 * t358 + t198 * t24;
t69 = qJ(6) * t154 + t382;
t223 = t198 * t235;
t83 = -t194 * t118 - t223;
t70 = -t154 * pkin(5) - t83;
t260 = t194 * t69 - t198 * t70;
t259 = t194 * t382 + t198 * t83;
t256 = t347 + t349;
t255 = -t66 * t196 - t65 * t199;
t254 = t103 * t275;
t38 = t103 * t115 + t154 * t64;
t250 = -t155 * t185 + t336;
t248 = t185 * t188 - t198 * t301;
t161 = -pkin(4) - t266;
t247 = t189 * t29 + t8 * t194;
t246 = t188 * t93 - t348;
t45 = t103 * t189 + t350;
t243 = t155 * t188 + t340;
t95 = t154 * t188 - t327;
t36 = qJD(5) * t223 + t118 * t188 - t194 * t228 + t198 * t91;
t240 = t45 * pkin(11);
t239 = t196 * t321 - t199 * t317;
t238 = -t114 * t161 + t140 * t155 - t364;
t233 = -t184 * t42 + t275 * t94;
t232 = (-t154 * t368 + t155 * t195) * qJD(4);
t40 = -t114 * t265 + t155 * t375 + t92;
t230 = -t40 + (t155 * t161 - t363) * qJD(5);
t147 = -t306 + t161;
t227 = -t40 + (t147 * t155 - t336) * qJD(5);
t221 = -t194 * t275 - t286;
t214 = -t93 * t341 + (t189 * t93 + t351) * t155;
t212 = -qJD(5) * t264 + t385;
t211 = -qJD(5) * t263 - t7 * t194 - t5;
t30 = t313 - t36 + t343;
t31 = -t362 - t37;
t9 = -qJD(5) * t260 + t31 * t194 + t30 * t198;
t10 = -qJD(5) * t259 - t37 * t194 - t36 * t198;
t209 = -t114 * t147 + t123 * t155 - t154 * t275 - t339;
t207 = pkin(3) * t232 - t114 * t185 - t339;
t175 = -0.2e1 * t284;
t174 = 0.2e1 * t284;
t158 = -0.2e1 * t273;
t146 = t161 * t188;
t138 = t389 * t368 * t357;
t133 = t147 * t188;
t131 = t184 * t188 - t251;
t128 = 0.2e1 * t138;
t112 = t117 * t188;
t99 = t272 + t136;
t97 = t115 * t194 + t154 * t189;
t79 = -0.2e1 * t155 * t332 - 0.2e1 * t274;
t78 = -0.2e1 * t191 * t335 + 0.2e1 * t274;
t76 = t155 * t380 + t297;
t68 = t151 * t380 + 0.2e1 * t155 * t297;
t67 = -t114 * t191 + 0.4e1 * t155 * t284 + t332;
t60 = t115 * t334 + t154 * t98;
t58 = -0.2e1 * t154 * t243 + 0.2e1 * t155 * t327;
t41 = pkin(11) * t348;
t35 = -0.2e1 * t94 * t42;
t32 = t189 * t94 - t352;
t20 = -t94 * t340 + (-t188 * t94 - t346) * t155;
t19 = -0.2e1 * t103 * t42 + 0.2e1 * t64 * t94;
t11 = -t103 * t243 + t115 * t94 - t154 * t42 + t333 * t64;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t273, 0.2e1 * (-t197 ^ 2 + t200 ^ 2) * t190 * qJD(2), 0.2e1 * t277 * t330, t158, -0.2e1 * t193 * t269, 0, -0.2e1 * pkin(1) * t190 * t321 - 0.2e1 * t136 * t344, -0.2e1 * pkin(1) * t289 + 0.2e1 * t135 * t344 (-t135 * t200 + t136 * t197 + (-t145 * t197 + t200 * t229) * qJD(2)) * t372, -0.2e1 * t135 * t145 + 0.2e1 * t136 * t229, -0.2e1 * t236 * t377, -0.2e1 * t236 * t296 + 0.2e1 * t237 * t377, 0.2e1 * t176 * t236 + 0.2e1 * t330 * t377, 0.2e1 * t237 * t296, -0.2e1 * t176 * t237 + 0.2e1 * t296 * t330, t158, -0.2e1 * t136 * t278 + 0.2e1 * t129 * t296 + 0.2e1 * (-t66 * t200 + (qJD(2) * t100 + t337) * t197) * t193, -0.2e1 * t101 * t176 - 0.2e1 * t129 * t377 + 0.2e1 * t136 * t236 - 0.2e1 * t330 * t65, 0.2e1 * t100 * t377 - 0.2e1 * t101 * t296 - 0.2e1 * t236 * t66 + 0.2e1 * t237 * t65, 0.2e1 * t100 * t66 - 0.2e1 * t101 * t65 + 0.2e1 * t129 * t136, -0.2e1 * t104 * t203, 0.2e1 * t103 * t203 - 0.2e1 * t104 * t64, 0.2e1 * t104 * t176 + 0.2e1 * t203 * t330, t51 (-t103 * t321 + t200 * t64) * t372, t158, 0.2e1 * t99 * t103 + 0.2e1 * t107 * t64 + 0.2e1 * (-t18 * t200 + t321 * t52) * t193, 0.2e1 * t99 * t104 - 0.2e1 * t107 * t203 - 0.2e1 * t17 * t330 - 0.2e1 * t176 * t53, 0.2e1 * t17 * t103 - 0.2e1 * t18 * t104 + 0.2e1 * t203 * t52 - 0.2e1 * t53 * t64, 0.2e1 * t107 * t99 - 0.2e1 * t17 * t53 + 0.2e1 * t18 * t52, t35, 0.2e1 * t268, t19, t309, -0.2e1 * t267, t51, 0.2e1 * t103 * t7 + 0.2e1 * t16 * t93 + 0.2e1 * t24 * t64 + 0.2e1 * t43 * t49, 0.2e1 * t103 * t6 + 0.2e1 * t16 * t94 - 0.2e1 * t358 * t64 - 0.2e1 * t42 * t49, 0.2e1 * t24 * t42 - 0.2e1 * t358 * t43 + 0.2e1 * t6 * t93 - 0.2e1 * t7 * t94, 0.2e1 * t16 * t49 + 0.2e1 * t24 * t7 - 0.2e1 * t358 * t6, t35, t19, -0.2e1 * t268, t51, 0.2e1 * t267, t309, -0.2e1 * t103 * t4 - 0.2e1 * t23 * t64 + 0.2e1 * t29 * t43 + 0.2e1 * t8 * t93, -0.2e1 * t2 * t93 - 0.2e1 * t22 * t43 - 0.2e1 * t23 * t42 + 0.2e1 * t4 * t94, 0.2e1 * t103 * t2 + 0.2e1 * t22 * t64 + 0.2e1 * t29 * t42 - 0.2e1 * t8 * t94, 0.2e1 * t2 * t22 + 0.2e1 * t23 * t4 + 0.2e1 * t29 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, 0, -t176, 0, -t136, t135, 0, 0, t374 * t288 + (-t162 + (0.2e1 * t276 + t285) * t199) * t196, -t196 * t296 - t236 * t319 + (-t377 - t222) * t199, t239 * t193, t381 (t199 * t321 + t287) * t193, 0, -pkin(2) * t296 - t136 * t199 + t196 * t312 - t239 * t365, pkin(2) * t377 - t287 * t365 + t337 + (-pkin(9) * t176 + t312) * t199, -t100 * t318 - t101 * t319 + t236 * t302 + t255 + (-t196 * t377 + t381) * pkin(9), -t136 * pkin(2) + ((-t100 * t199 - t101 * t196) * qJD(3) + t255) * pkin(9), -t104 * t114 - t155 * t203, t114 * t103 - t104 * t115 + t154 * t203 - t155 * t64 (t114 * t200 + t155 * t321) * t193, t38 (t115 * t200 - t154 * t321) * t193, 0, t103 * t305 + t107 * t115 + t99 * t154 + t186 * t64 + (-t117 * t321 + t200 * t92) * t193, t104 * t305 - t107 * t114 - t118 * t176 + t99 * t155 - t186 * t203 - t330 * t91, t91 * t103 + t92 * t104 + t52 * t114 - t53 * t115 - t117 * t203 - t118 * t64 + t17 * t154 - t18 * t155, t107 * t305 - t117 * t18 - t118 * t17 + t186 * t99 - t52 * t92 - t53 * t91, t20, t386, t11, t214, -t378, t38, t103 * t37 + t115 * t24 + t117 * t43 + t154 * t7 + t155 * t360 - t341 * t49 + t64 * t83 + t92 * t93, t103 * t36 - t115 * t358 - t117 * t42 + t154 * t6 - t155 * t280 - t340 * t49 - t382 * t64 + t92 * t94, t36 * t93 - t37 * t94 + t42 * t83 - t43 * t382 + t263 * t114 + (t194 * t6 - t198 * t7 + (t194 * t24 - t198 * t358) * qJD(5)) * t155, t117 * t16 + t24 * t37 - t358 * t36 - t382 * t6 + t49 * t92 + t7 * t83, t20, t11, -t386, t38, t378, t214, -t103 * t31 - t115 * t23 - t154 * t4 + t155 * t247 - t29 * t341 + t40 * t93 + t43 * t90 - t64 * t70, -t30 * t93 + t31 * t94 - t42 * t70 - t43 * t69 + t264 * t114 + (-t194 * t2 + t198 * t4 + (-t194 * t23 - t198 * t22) * qJD(5)) * t155, t103 * t30 + t115 * t22 + t154 * t2 + t155 * t282 + t29 * t340 - t40 * t94 + t42 * t90 + t64 * t69, t2 * t69 + t22 * t30 + t23 * t31 + t29 * t40 + t4 * t70 + t8 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t283, 0.2e1 * (-t196 ^ 2 + t374) * qJD(3), 0, -0.2e1 * t283, 0, 0, t196 * t307, t199 * t307, 0, 0, -0.2e1 * t335, 0.2e1 * t114 * t154 - 0.2e1 * t115 * t155, 0, t102, 0, 0, 0.2e1 * t115 * t186 + 0.2e1 * t154 * t305, -0.2e1 * t114 * t186 + 0.2e1 * t155 * t305, -0.2e1 * t114 * t117 - 0.2e1 * t115 * t118 + 0.2e1 * t154 * t91 + 0.2e1 * t155 * t92, -0.2e1 * t118 * t91 + 0.2e1 * t186 * t305 + 0.2e1 * t355, t79, 0.2e1 * t68, t58, t78, -0.2e1 * t60, t102, 0.2e1 * t115 * t83 + 0.2e1 * t117 * t98 + 0.2e1 * t154 * t37 + 0.2e1 * t334 * t92, -0.2e1 * t115 * t382 - 0.2e1 * t117 * t243 + 0.2e1 * t154 * t36 + 0.2e1 * t333 * t92, 0.2e1 * t259 * t114 + 0.2e1 * (t194 * t36 - t198 * t37 + (t194 * t83 - t198 * t382) * qJD(5)) * t155, -0.2e1 * t36 * t382 + 0.2e1 * t37 * t83 + 0.2e1 * t355, t79, t58, -0.2e1 * t68, t102, 0.2e1 * t60, t78, -0.2e1 * t90 * t341 - 0.2e1 * t115 * t70 - 0.2e1 * t154 * t31 + 0.2e1 * (t189 * t90 + t40 * t194) * t155, 0.2e1 * t260 * t114 + 0.2e1 * (-t194 * t30 + t198 * t31 + (-t194 * t70 - t198 * t69) * qJD(5)) * t155, 0.2e1 * t90 * t340 + 0.2e1 * t115 * t69 + 0.2e1 * t154 * t30 + 0.2e1 * (-t40 * t198 + t85) * t155, 0.2e1 * t30 * t69 + 0.2e1 * t31 * t70 + 0.2e1 * t40 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t377, 0, -t296, t176, t66, t65, 0, 0, 0, 0, -t203, 0, -t64, t176 (t200 * t316 + t321 * t368) * t367 + t18 (-t195 * t321 + t200 * t281) * t367 + t17, t104 * t301 + t203 * t306 - t366 * t64 - t254 (t368 * t18 - t17 * t195 + (-t195 * t52 + t368 * t53) * qJD(4)) * pkin(3), t32, -t12, t45, t246, -t44, 0, t185 * t43 - t45 * t184 + (-t103 * t293 + t195 * t93) * t357 + t280, -t185 * t42 + t44 * t184 + (-t103 * t292 + t195 * t94) * t357 + t360, -t5 + (t184 * t94 - t24) * t189 + (-t7 + (-t358 + t353) * qJD(5) + t233) * t194 + t359, t16 * t185 + (t195 * t49 - t24 * t293 + t292 * t358) * t357 + t211 * t184, t32, t45, t12, 0, t44, t246, t103 * t221 + t123 * t93 + t147 * t43 - t184 * t350 + t282, t94 * t286 + ((-t22 + t353) * qJD(5) + t233) * t194 + t308 + t359, -t123 * t94 + t147 * t42 + (-qJD(5) * t103 * t184 - t8) * t194 + (-qJD(5) * t29 + t184 * t64 + t254) * t198, t29 * t123 + t8 * t147 + (t22 * t292 + t23 * t293) * t357 + t212 * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, 0, -t319, 0, -t302, pkin(9) * t319, 0, 0, 0, 0, -t114, 0, -t115, 0, -t92, t91 (t114 * t368 - t115 * t195 + t232) * pkin(3) (-t368 * t92 - t195 * t91 + (t118 * t368 + t338) * qJD(4)) * pkin(3), -t76, -t67, t97, t76, -t95, 0, t112 + (-qJD(5) * t250 - t92) * t198 + t207 * t194, t188 * t250 + t198 * t207 + t345, t10, t92 * t185 + (t292 * t382 - t293 * t83 + t338) * t357 + t10 * t184, -t76, t97, t67, 0, t95, t76, t194 * t209 + t198 * t227 + t85, t9, t227 * t194 + (-t209 - t342) * t198, t90 * t123 + t40 * t147 + (t292 * t69 + t293 * t70) * t357 + t9 * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t301, -0.2e1 * t275, 0, 0, t174, t373, 0, t175, 0, 0, 0.2e1 * t248, 0.2e1 * t323, t128, 0.2e1 * t185 * t301 + 0.2e1 * t325, t174, 0, -t373, 0, 0, t175, -0.2e1 * t123 * t198 + 0.2e1 * t133, t128, -0.2e1 * t123 * t194 - 0.2e1 * t147 * t189, 0.2e1 * t123 * t147 + 0.2e1 * t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203, 0, -t64, t176, t18, t17, 0, 0, t32, -t12, t45, t246, -t44, 0, -pkin(4) * t43 - t240 + t280, pkin(4) * t42 + t360 + t384, -t41 - t5 + (-pkin(11) * t42 - t7) * t194 + (pkin(11) * t256 - t263) * qJD(5), -t16 * pkin(4) + pkin(11) * t211, t32, t45, t12, 0, t44, t246, t140 * t93 + t161 * t43 - t240 + t282, -t22 * t188 - t41 + (qJD(5) * t256 - t352) * pkin(11) + t308, -t140 * t94 + t161 * t42 - t247 - t384, pkin(11) * t212 + t140 * t29 + t161 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, 0, -t115, 0, -t92, t91, 0, 0, -t76, -t67, t97, t76, -t95, 0, t112 + (pkin(4) * t114 - t364) * t194 + (-t92 + (-pkin(4) * t155 - t363) * qJD(5)) * t198, pkin(4) * t243 + pkin(11) * t95 + t345, t10, -t92 * pkin(4) + pkin(11) * t10, -t76, t97, t67, 0, t95, t76, t194 * t238 + t198 * t230 + t85, t9, t230 * t194 + (-t238 - t342) * t198, pkin(11) * t9 + t140 * t90 + t161 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, -t275, 0, 0, t174, t373, 0, t175, 0, 0, t248 - t304, -t303 + t323, t138, -pkin(4) * t301 + t324, t174, 0, -t373, 0, 0, t175, t198 * t326 + t133 + t146, t138, t326 * t194 + (-t147 - t161) * t189, t123 * t161 + t140 * t147 + t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t373, 0, t175, 0, 0, -0.2e1 * t304, -0.2e1 * t303, 0, 0, t174, 0, -t373, 0, 0, t175, -0.2e1 * t140 * t198 + 0.2e1 * t146, 0, -0.2e1 * t140 * t194 - 0.2e1 * t161 * t189, 0.2e1 * t161 * t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, 0, -t43, t64, t7, t6, 0, 0, 0, -t42, 0, t64, t43, 0, t7 + 0.2e1 * t369, pkin(5) * t42 - qJ(6) * t43 - qJD(6) * t93, 0.2e1 * t314 - t6 + 0.2e1 * t356, -pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, 0, -t98, t115, t37, t36, 0, 0, 0, -t243, 0, t115, t98, 0, t37 + 0.2e1 * t362, t266 * t114 + (qJD(5) * t265 - t311) * t155, 0.2e1 * t313 - t36 + 0.2e1 * t343, -pkin(5) * t31 + qJ(6) * t30 + qJD(6) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, -t188, 0, t221, t131, 0, 0, 0, t189, 0, 0, t188, 0, t221, -t375, -t131 (-pkin(5) * t293 + qJ(6) * t292) * t357 - t375 * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, -t188, 0, -t299, t300, 0, 0, 0, t189, 0, 0, t188, 0, -t299, -t375, -t300, -t375 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t371, qJ(6) * t371; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t42, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t243, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, -t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
