% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRRP9
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRRP9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:40
% EndTime: 2019-03-10 02:12:10
% DurationCPUTime: 12.75s
% Computational Cost: add. (15394->648), mult. (40717->1106), div. (0->0), fcn. (39476->10), ass. (0->293)
t189 = sin(qJ(4));
t192 = cos(qJ(4));
t188 = sin(pkin(6));
t194 = cos(qJ(2));
t348 = t188 * t194;
t162 = t192 * t348;
t190 = sin(qJ(3));
t351 = cos(pkin(6));
t284 = t351 * t190;
t191 = sin(qJ(2));
t193 = cos(qJ(3));
t345 = t191 * t193;
t237 = t188 * t345 + t284;
t230 = qJD(4) * t237;
t335 = qJD(2) * t191;
t293 = t188 * t335;
t349 = t188 * t191;
t307 = qJD(3) * t349;
t154 = t190 * t307;
t282 = t351 * qJD(3);
t334 = qJD(2) * t194;
t292 = t188 * t334;
t367 = -t193 * (t282 + t292) + t154;
t212 = -qJD(4) * t162 - t192 * t367 + (-t230 + t293) * t189;
t114 = -t189 * t348 + t192 * t237;
t327 = t114 * qJD(4);
t374 = t212 * t189 + t192 * t327;
t311 = pkin(1) * t351;
t271 = t194 * t311;
t322 = pkin(8) * t349;
t234 = -t271 + t322;
t227 = t234 * qJD(2);
t314 = -pkin(9) * t191 - pkin(1);
t247 = -pkin(2) * t194 + t314;
t242 = t247 * t188;
t380 = -qJD(3) * t242 + t227;
t329 = qJD(4) * t192;
t332 = qJD(3) * t193;
t379 = t189 * t332 + t190 * t329;
t339 = t237 * t189 + t162;
t362 = sin(qJ(5));
t256 = t362 * t339;
t363 = cos(qJ(5));
t285 = t363 * qJD(5);
t330 = qJD(4) * t189;
t298 = t194 * t330;
t340 = t189 * t367 - t192 * t230;
t366 = -t188 * (t192 * t335 + t298) - t340;
t34 = -qJD(5) * t256 + t114 * t285 + t212 * t362 + t363 * t366;
t257 = t363 * t339;
t75 = t114 * t362 + t257;
t378 = t34 * qJ(6) + t75 * qJD(6);
t148 = t189 * t363 + t192 * t362;
t136 = t148 * t190;
t288 = qJD(3) * t362;
t264 = t193 * t288;
t289 = qJD(3) * t363;
t265 = t193 * t289;
t309 = t362 * t189;
t269 = t190 * t309;
t368 = qJD(4) + qJD(5);
t372 = t363 * qJD(4) + t285;
t80 = t189 * t265 + (t372 * t190 + t264) * t192 - t368 * t269;
t377 = t80 * qJ(6) + t136 * qJD(6);
t232 = -pkin(2) * t351 - t271;
t132 = t232 + t322;
t138 = t190 * t349 - t193 * t351;
t356 = t138 * pkin(3);
t209 = -pkin(10) * t237 + t132 + t356;
t272 = t191 * t311;
t231 = pkin(9) * t351 + t272;
t320 = pkin(8) * t348;
t219 = t231 + t320;
t216 = qJD(3) * t219;
t260 = pkin(2) * t191 - pkin(9) * t194;
t336 = qJD(2) * t188;
t241 = t260 * t336;
t52 = t380 * t193 + (t216 - t241) * t190;
t376 = -pkin(10) * t293 - qJD(4) * t209 + t52;
t229 = t237 * qJD(3);
t113 = t190 * t292 + t229;
t226 = t193 * t231;
t207 = t226 + (t190 * t314 + (-pkin(2) * t190 + pkin(8) * t193 - pkin(10)) * t194) * t188;
t283 = qJD(2) * t351;
t261 = t191 * t283;
t375 = -pkin(1) * t261 - t113 * pkin(3) - pkin(8) * t292 - pkin(10) * t367 + qJD(4) * t207;
t279 = t339 * qJD(4);
t373 = t189 * t279 - t192 * t366;
t300 = t190 * t330;
t303 = t192 * t332;
t370 = t300 - t303;
t369 = -t113 * t192 + t138 * t330;
t181 = qJD(3) * t190;
t86 = -t193 * t113 + t138 * t181;
t184 = t189 ^ 2;
t186 = t192 ^ 2;
t338 = t184 - t186;
t281 = qJD(4) * t338;
t354 = t193 * pkin(3);
t364 = -pkin(11) - pkin(10);
t365 = t190 * t364 - pkin(2) - t354;
t187 = t193 ^ 2;
t361 = pkin(3) * t190;
t360 = pkin(3) * t192;
t359 = pkin(9) * t188;
t358 = pkin(9) * t189;
t357 = pkin(9) * t193;
t355 = t190 * pkin(4);
t182 = t190 * pkin(9);
t276 = pkin(4) * t285;
t318 = t362 * pkin(4);
t353 = -t75 * t276 - t34 * t318;
t199 = t138 * pkin(4) - t114 * pkin(11) - t189 * t207 + t192 * t209;
t39 = t362 * t199;
t48 = t189 * t209 + t192 * t207;
t44 = -pkin(11) * t339 + t48;
t23 = t363 * t44 + t39;
t352 = -t136 * t276 - t80 * t318;
t347 = t189 * t190;
t346 = t190 * t192;
t344 = t192 * t193;
t342 = t193 * t194;
t112 = t368 * t148;
t310 = t363 * t192;
t147 = t309 - t310;
t341 = -t112 * t318 - t147 * t276;
t233 = t365 * t192;
t312 = -pkin(4) - t358;
t215 = t193 * t312 + t233;
t100 = t362 * t215;
t173 = pkin(9) * t344;
t259 = -pkin(10) * t190 - t354;
t246 = pkin(2) - t259;
t127 = -t189 * t246 + t173;
t115 = -pkin(11) * t347 + t127;
t73 = t363 * t115 + t100;
t159 = t364 * t192;
t273 = t364 * t362;
t117 = -t363 * t159 + t189 * t273;
t152 = pkin(4) * t347 + t182;
t185 = t190 ^ 2;
t337 = t185 - t187;
t333 = qJD(3) * t192;
t331 = qJD(3) * t194;
t328 = qJD(4) * t193;
t326 = t132 * qJD(3);
t325 = -0.2e1 * pkin(2) * qJD(3);
t324 = -0.2e1 * pkin(3) * qJD(4);
t93 = 0.2e1 * t138 * t113;
t323 = t189 * t357;
t321 = pkin(9) * t346;
t319 = t363 * pkin(4);
t317 = pkin(4) * t330;
t180 = pkin(9) * t332;
t316 = t364 * t193;
t315 = -t190 * t380 + t193 * t216;
t123 = pkin(4) * t379 + t180;
t178 = -pkin(4) * t192 - pkin(3);
t183 = t188 ^ 2;
t308 = t183 * t334;
t304 = t190 * t331;
t302 = t138 * t332;
t299 = t189 * t328;
t296 = t192 * t328;
t295 = t114 * t332;
t291 = t189 * t329;
t290 = t190 * t332;
t287 = qJD(5) * t362;
t280 = t339 * qJD(3);
t278 = t337 * qJD(3);
t277 = 0.2e1 * t290;
t275 = pkin(4) * t287;
t274 = t364 * t363;
t76 = t114 * t363 - t256;
t270 = t76 * t287;
t268 = t192 * t290;
t267 = t185 * t291;
t266 = t191 * t308;
t137 = t190 * t310 - t269;
t263 = t137 * t287;
t262 = t148 * t287;
t40 = t363 * t199;
t22 = -t362 * t44 + t40;
t258 = -pkin(10) * t193 + t361;
t255 = t193 * t280;
t101 = t363 * t215;
t72 = -t115 * t362 + t101;
t150 = t189 * t274;
t116 = t159 * t362 + t150;
t24 = t375 * t189 + t376 * t192;
t25 = t376 * t189 - t375 * t192;
t253 = -t25 * t189 - t24 * t192;
t53 = t193 * t241 - t315;
t252 = -t53 * t190 - t52 * t193;
t251 = t138 * t275;
t225 = t192 * t246 + t323;
t250 = -t127 * t189 + t192 * t225;
t249 = qJD(4) * t274;
t248 = qJD(4) * t273;
t51 = (-t191 * pkin(3) - t193 * t260) * t336 + t315;
t91 = -t190 * t219 + t193 * t242;
t85 = pkin(3) * t348 - t91;
t245 = t51 * t189 + t329 * t85;
t244 = t189 * t258;
t243 = t189 * t113 + t138 * t329;
t196 = t113 * pkin(4) - pkin(11) * t212 + t25;
t197 = -pkin(11) * t366 - t24;
t3 = -qJD(5) * t40 - t362 * t196 - t363 * t197 + t287 * t44;
t203 = (-t189 * t365 - t173) * qJD(4) + (t192 * t316 + (-t312 + t360) * t190) * qJD(3);
t204 = (t233 - t323) * qJD(4) + (-t321 + (t316 + t361) * t189) * qJD(3);
t35 = -qJD(5) * t101 + t115 * t287 - t362 * t203 - t363 * t204;
t240 = t190 * t335 - t193 * t331;
t81 = -qJD(5) * t150 - t159 * t287 - t189 * t249 - t192 * t248;
t143 = t272 + t320;
t228 = -t113 * t318 - t138 * t276 + t3;
t224 = pkin(8) * t342 + t190 * t247;
t89 = t225 * qJD(4) + (-t244 + t321) * qJD(3);
t90 = -t127 * qJD(4) + (pkin(9) * t347 + t192 * t258) * qJD(3);
t218 = qJD(4) * t250 - t189 * t90 - t192 * t89;
t214 = t193 * t276 - t288 * t355 + t35;
t62 = pkin(4) * t339 + t85;
t213 = t189 * t366 + t192 * t279;
t82 = t159 * t285 + t192 * t249 + (-qJD(5) * t273 - t248) * t189;
t210 = t212 * t192;
t37 = -t340 * pkin(4) + (-pkin(4) * t298 + (pkin(9) * t342 + (-pkin(2) * t193 + t178) * t191) * qJD(2)) * t188 + t315;
t202 = t363 * t203 - t204 * t362;
t200 = -qJD(5) * t100 - t115 * t285 + t202;
t79 = t112 * t190 + t189 * t264 - t192 * t265;
t198 = t79 * qJ(6) - t137 * qJD(6) + t200;
t179 = pkin(5) * t181;
t26 = t179 + t198;
t4 = -qJD(5) * t39 + t363 * t196 - t197 * t362 - t285 * t44;
t33 = qJD(5) * t257 + t114 * t287 - t212 * t363 + t362 * t366;
t195 = t33 * qJ(6) - t76 * qJD(6) + t4;
t110 = t113 * pkin(5);
t1 = t110 + t195;
t177 = t319 + pkin(5);
t175 = -0.2e1 * t276;
t174 = -0.2e1 * t275;
t168 = -0.2e1 * t290;
t165 = t193 * t275;
t153 = -0.2e1 * t266;
t135 = qJD(2) * t143;
t130 = pkin(5) * t147 + t178;
t118 = -t189 * t303 + t190 * t281;
t111 = (qJD(4) * t362 + t287) * t189 - t372 * t192;
t108 = pkin(5) * t136 + t152;
t98 = pkin(5) * t112 + t317;
t97 = -qJ(6) * t147 + t117;
t96 = -t148 * qJ(6) + t116;
t95 = -0.2e1 * t148 * t111;
t94 = 0.2e1 * t147 * t112;
t92 = t188 * t224 + t226;
t88 = t112 * t193 - t147 * t181;
t87 = t111 * t193 + t148 * t181;
t64 = -0.2e1 * t137 * t79;
t63 = 0.2e1 * t136 * t80;
t61 = pkin(5) * t80 + t123;
t60 = -qJ(6) * t136 + t73;
t59 = -0.2e1 * t136 * t181 + 0.2e1 * t193 * t80;
t58 = 0.2e1 * t137 * t181 + 0.2e1 * t193 * t79;
t57 = -t112 * t138 - t113 * t147;
t56 = -t111 * t138 + t113 * t148;
t55 = 0.2e1 * t111 * t147 - 0.2e1 * t112 * t148;
t54 = -t193 * pkin(5) - t137 * qJ(6) + t72;
t50 = t111 * qJ(6) - t148 * qJD(6) + t82;
t49 = t112 * qJ(6) + t147 * qJD(6) + t81;
t47 = -t189 * t226 + t192 * (-pkin(10) * t284 + t232 + t356) + (-t189 * (-t194 * pkin(10) + t224) + t192 * (pkin(8) * t191 - pkin(10) * t345)) * t188;
t46 = t112 * t136 + t147 * t80;
t45 = -t111 * t137 - t148 * t79;
t42 = t75 * pkin(5) + t62;
t41 = 0.2e1 * t136 * t79 - 0.2e1 * t137 * t80;
t36 = -qJD(5) * t73 + t202;
t30 = t111 * t136 - t112 * t137 + t147 * t79 - t148 * t80;
t29 = -0.2e1 * t76 * t33;
t28 = 0.2e1 * t75 * t34;
t27 = t35 + t377;
t21 = t112 * t75 + t147 * t34;
t20 = -t111 * t76 - t148 * t33;
t19 = -0.2e1 * t113 * t75 - 0.2e1 * t138 * t34;
t18 = 0.2e1 * t113 * t76 - 0.2e1 * t138 * t33;
t17 = t136 * t34 + t75 * t80;
t16 = -t137 * t33 - t76 * t79;
t14 = -qJ(6) * t75 + t23;
t13 = t34 * pkin(5) + t37;
t12 = t138 * pkin(5) - t76 * qJ(6) + t22;
t9 = -t113 * t136 - t138 * t80 - t181 * t75 + t193 * t34;
t8 = t113 * t137 - t138 * t79 + t181 * t76 + t193 * t33;
t7 = 0.2e1 * t33 * t75 - 0.2e1 * t34 * t76;
t6 = t111 * t75 - t112 * t76 + t147 * t33 - t148 * t34;
t5 = t136 * t33 - t137 * t34 + t75 * t79 - t76 * t80;
t2 = t3 + t378;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t266, 0.2e1 * (-t191 ^ 2 + t194 ^ 2) * t183 * qJD(2), 0.2e1 * t283 * t348, t153, -0.2e1 * t188 * t261, 0, -0.2e1 * pkin(1) * t183 * t335 - 0.2e1 * t135 * t351, -0.2e1 * pkin(1) * t308 + 0.2e1 * t227 * t351, 0.2e1 * t135 * t349 - 0.2e1 * t143 * t293 - 0.2e1 * t227 * t348 + 0.2e1 * t234 * t292, 0.2e1 * t135 * t234 - 0.2e1 * t143 * t227, -0.2e1 * t237 * t367, -0.2e1 * t113 * t237 + 0.2e1 * t138 * t367, 0.2e1 * t237 * t293 + 0.2e1 * t348 * t367, t93, 0.2e1 * (t113 * t194 - t138 * t335) * t188, t153, 0.2e1 * t132 * t113 + 0.2e1 * t135 * t138 + 0.2e1 * (-t194 * t53 + t335 * t91) * t188, -0.2e1 * t132 * t367 + 0.2e1 * t135 * t237 - 0.2e1 * t293 * t92 - 0.2e1 * t348 * t52, -0.2e1 * t92 * t113 + 0.2e1 * t52 * t138 - 0.2e1 * t237 * t53 + 0.2e1 * t367 * t91, 0.2e1 * t132 * t135 - 0.2e1 * t52 * t92 + 0.2e1 * t53 * t91, 0.2e1 * t114 * t212, -0.2e1 * t114 * t366 - 0.2e1 * t212 * t339, 0.2e1 * t114 * t113 + 0.2e1 * t138 * t212, 0.2e1 * t339 * t366, -0.2e1 * t113 * t339 - 0.2e1 * t138 * t366, t93, 0.2e1 * t47 * t113 + 0.2e1 * t25 * t138 + 0.2e1 * t339 * t51 + 0.2e1 * t366 * t85, -0.2e1 * t48 * t113 + 0.2e1 * t51 * t114 + 0.2e1 * t24 * t138 + 0.2e1 * t212 * t85, -0.2e1 * t25 * t114 - 0.2e1 * t212 * t47 + 0.2e1 * t24 * t339 - 0.2e1 * t366 * t48, -0.2e1 * t24 * t48 + 0.2e1 * t25 * t47 + 0.2e1 * t51 * t85, t29, t7, t18, t28, t19, t93, 0.2e1 * t113 * t22 + 0.2e1 * t138 * t4 + 0.2e1 * t34 * t62 + 0.2e1 * t37 * t75, -0.2e1 * t113 * t23 + 0.2e1 * t138 * t3 - 0.2e1 * t33 * t62 + 0.2e1 * t37 * t76, 0.2e1 * t22 * t33 - 0.2e1 * t23 * t34 + 0.2e1 * t3 * t75 - 0.2e1 * t4 * t76, 0.2e1 * t22 * t4 - 0.2e1 * t23 * t3 + 0.2e1 * t37 * t62, t29, t7, t18, t28, t19, t93, 0.2e1 * t1 * t138 + 0.2e1 * t113 * t12 + 0.2e1 * t13 * t75 + 0.2e1 * t34 * t42, -0.2e1 * t113 * t14 + 0.2e1 * t13 * t76 + 0.2e1 * t138 * t2 - 0.2e1 * t33 * t42, -0.2e1 * t1 * t76 + 0.2e1 * t12 * t33 - 0.2e1 * t14 * t34 + 0.2e1 * t2 * t75, 0.2e1 * t1 * t12 + 0.2e1 * t13 * t42 - 0.2e1 * t14 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, 0, -t293, 0, -t135, t227, 0, 0, t187 * t307 + (-t154 + (0.2e1 * t282 + t292) * t193) * t190, -t193 * t367 - t302 + (-t113 - t229) * t190, t240 * t188, t86 (t193 * t335 + t304) * t188, 0, -pkin(2) * t113 - t135 * t193 + t190 * t326 - t240 * t359, pkin(2) * t367 + t135 * t190 + t193 * t326 - t293 * t357 - t304 * t359, pkin(9) * t86 + t180 * t237 - t181 * t92 - t182 * t367 - t332 * t91 + t252, -pkin(2) * t135 + ((-t190 * t92 - t193 * t91) * qJD(3) + t252) * pkin(9), -t114 * t370 + t190 * t210, -t189 * t295 - t192 * t255 + (t373 - t374) * t190, t113 * t346 + t114 * t181 - t138 * t300 + t192 * t302 - t193 * t212, t189 * t255 + t190 * t213 (-t189 * qJD(3) * t138 + t366) * t193 + (-t280 - t243) * t190, t86, -t225 * t113 + t90 * t138 + (-t25 + (pkin(9) * t339 + t85 * t189) * qJD(3)) * t193 + (pkin(9) * t366 + t47 * qJD(3) + t245) * t190, pkin(9) * t295 - t127 * t113 + t89 * t138 - t181 * t48 + t212 * t182 - t24 * t193 + t51 * t346 - t370 * t85, -t90 * t114 - t127 * t366 + t212 * t225 + t24 * t347 - t25 * t346 + t339 * t89 + t370 * t47 - t379 * t48, -t225 * t25 - t127 * t24 + t47 * t90 - t48 * t89 + (t51 * t190 + t332 * t85) * pkin(9), t16, t5, t8, t17, t9, t86, t113 * t72 + t123 * t75 + t136 * t37 + t138 * t36 + t152 * t34 + t181 * t22 - t193 * t4 + t62 * t80, -t113 * t73 + t123 * t76 + t137 * t37 + t138 * t35 - t152 * t33 - t181 * t23 - t193 * t3 - t62 * t79, t136 * t3 - t137 * t4 + t22 * t79 - t23 * t80 + t33 * t72 - t34 * t73 + t35 * t75 - t36 * t76, t123 * t62 + t152 * t37 + t22 * t36 - t23 * t35 - t3 * t73 + t4 * t72, t16, t5, t8, t17, t9, t86, -t1 * t193 + t108 * t34 + t113 * t54 + t12 * t181 + t13 * t136 + t138 * t26 + t42 * t80 + t61 * t75, -t108 * t33 - t113 * t60 + t13 * t137 + t138 * t27 - t14 * t181 - t193 * t2 - t42 * t79 + t61 * t76, -t1 * t137 + t12 * t79 + t136 * t2 - t14 * t80 - t26 * t76 + t27 * t75 + t33 * t54 - t34 * t60, t1 * t54 + t108 * t13 + t12 * t26 - t14 * t27 - t2 * t60 + t42 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, -0.2e1 * t278, 0, t168, 0, 0, t190 * t325, t193 * t325, 0, 0, 0.2e1 * t186 * t290 - 0.2e1 * t267, 0.2e1 * t185 * t281 - 0.4e1 * t189 * t268, 0.2e1 * t190 * t299 + 0.2e1 * t333 * t337, 0.2e1 * t184 * t290 + 0.2e1 * t267, -0.2e1 * t189 * t278 + 0.2e1 * t190 * t296, t168, -0.2e1 * t225 * t181 - 0.2e1 * t193 * t90 + 0.2e1 * (t185 * t329 + t189 * t277) * pkin(9), -0.2e1 * t127 * t181 - 0.2e1 * t193 * t89 + 0.2e1 * (-t185 * t330 + 0.2e1 * t268) * pkin(9), 0.2e1 * t250 * t332 + 0.2e1 * (t189 * t89 - t192 * t90 + (-t127 * t192 - t189 * t225) * qJD(4)) * t190, 0.2e1 * pkin(9) ^ 2 * t290 - 0.2e1 * t127 * t89 - 0.2e1 * t225 * t90, t64, t41, t58, t63, t59, t168, 0.2e1 * t123 * t136 + 0.2e1 * t152 * t80 + 0.2e1 * t181 * t72 - 0.2e1 * t193 * t36, 0.2e1 * t123 * t137 - 0.2e1 * t152 * t79 - 0.2e1 * t181 * t73 - 0.2e1 * t193 * t35, 0.2e1 * t136 * t35 - 0.2e1 * t137 * t36 + 0.2e1 * t72 * t79 - 0.2e1 * t73 * t80, 0.2e1 * t123 * t152 - 0.2e1 * t35 * t73 + 0.2e1 * t36 * t72, t64, t41, t58, t63, t59, t168, 0.2e1 * t108 * t80 + 0.2e1 * t136 * t61 + 0.2e1 * t181 * t54 - 0.2e1 * t193 * t26, -0.2e1 * t108 * t79 + 0.2e1 * t137 * t61 - 0.2e1 * t181 * t60 - 0.2e1 * t193 * t27, 0.2e1 * t136 * t27 - 0.2e1 * t137 * t26 + 0.2e1 * t54 * t79 - 0.2e1 * t60 * t80, 0.2e1 * t108 * t61 + 0.2e1 * t26 * t54 - 0.2e1 * t27 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t367, 0, -t113, t293, t53, t52, 0, 0, t374, -t189 * t327 + t210 - t213, t243, t373, -t369, 0, -pkin(3) * t366 - pkin(10) * t243 - t51 * t192 + t330 * t85, -pkin(3) * t212 + t369 * pkin(10) + t245, -t329 * t47 - t330 * t48 + t253 + (t373 + t374) * pkin(10), -pkin(3) * t51 + ((-t48 * t189 - t47 * t192) * qJD(4) + t253) * pkin(10), t20, t6, t56, t21, t57, 0, t112 * t62 + t113 * t116 + t138 * t82 + t147 * t37 + t178 * t34 + t317 * t75, -t111 * t62 - t113 * t117 + t138 * t81 + t148 * t37 - t178 * t33 + t317 * t76, t111 * t22 - t112 * t23 + t116 * t33 - t117 * t34 + t147 * t3 - t148 * t4 + t75 * t81 - t76 * t82, t116 * t4 - t117 * t3 + t178 * t37 + t22 * t82 - t23 * t81 + t317 * t62, t20, t6, t56, t21, t57, 0, t112 * t42 + t113 * t96 + t13 * t147 + t130 * t34 + t138 * t50 + t75 * t98, -t111 * t42 - t113 * t97 + t13 * t148 - t130 * t33 + t138 * t49 + t76 * t98, -t1 * t148 + t111 * t12 - t112 * t14 + t147 * t2 + t33 * t96 - t34 * t97 + t49 * t75 - t50 * t76, t1 * t96 + t12 * t50 + t13 * t130 - t14 * t49 - t2 * t97 + t42 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, 0, -t181, 0, -t180, pkin(9) * t181, 0, 0, -t118, -0.4e1 * t190 * t291 - t332 * t338, t181 * t189 - t296, t118, t190 * t333 + t299, 0 (pkin(10) * t344 + (t358 - t360) * t190) * qJD(4) + (t189 * t259 - t173) * qJD(3) (t244 + t321) * qJD(4) + (t192 * t259 + t323) * qJD(3), t218, -pkin(3) * t180 + pkin(10) * t218, t45, t30, t87, t46, t88, 0, t112 * t152 + t116 * t181 + t123 * t147 + t136 * t317 + t178 * t80 - t193 * t82, -t111 * t152 - t117 * t181 + t123 * t148 + t137 * t317 - t178 * t79 - t193 * t81, t111 * t72 - t112 * t73 + t116 * t79 - t117 * t80 + t136 * t81 - t137 * t82 + t147 * t35 - t148 * t36, t116 * t36 - t117 * t35 + t123 * t178 + t152 * t317 + t72 * t82 - t73 * t81, t45, t30, t87, t46, t88, 0, t108 * t112 + t130 * t80 + t136 * t98 + t147 * t61 + t181 * t96 - t193 * t50, -t108 * t111 - t130 * t79 + t137 * t98 + t148 * t61 - t181 * t97 - t193 * t49, t111 * t54 - t112 * t60 + t136 * t49 - t137 * t50 + t147 * t27 - t148 * t26 + t79 * t96 - t80 * t97, t108 * t98 + t130 * t61 + t26 * t96 - t27 * t97 - t49 * t60 + t50 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t291, -0.2e1 * t281, 0, -0.2e1 * t291, 0, 0, t189 * t324, t192 * t324, 0, 0, t95, t55, 0, t94, 0, 0, 0.2e1 * t112 * t178 + 0.2e1 * t147 * t317, -0.2e1 * t111 * t178 + 0.2e1 * t148 * t317, 0.2e1 * t111 * t116 - 0.2e1 * t112 * t117 + 0.2e1 * t147 * t81 - 0.2e1 * t148 * t82, 0.2e1 * t116 * t82 - 0.2e1 * t117 * t81 + 0.2e1 * t178 * t317, t95, t55, 0, t94, 0, 0, 0.2e1 * t112 * t130 + 0.2e1 * t147 * t98, -0.2e1 * t111 * t130 + 0.2e1 * t148 * t98, 0.2e1 * t111 * t96 - 0.2e1 * t112 * t97 + 0.2e1 * t147 * t49 - 0.2e1 * t148 * t50, 0.2e1 * t130 * t98 - 0.2e1 * t49 * t97 + 0.2e1 * t50 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, 0, -t366, t113, t25, t24, 0, 0, 0, 0, -t33, 0, -t34, t113, t113 * t319 - t251 + t4, t228 (t33 * t363 + t270) * pkin(4) + t353 (-t362 * t3 + t363 * t4 + (-t22 * t362 + t23 * t363) * qJD(5)) * pkin(4), 0, 0, -t33, 0, -t34, t113, t177 * t113 + t1 - t251, t228 + t378, pkin(4) * t270 + t177 * t33 + t353, t1 * t177 + (-t362 * t2 + (-t12 * t362 + t14 * t363) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t370, 0, -t379, t181, t90, t89, 0, 0, 0, 0, -t79, 0, -t80, t181, t289 * t355 + t165 + t200, t214 (t363 * t79 + t263) * pkin(4) + t352 (-t362 * t35 + t363 * t36 + (-t362 * t72 + t363 * t73) * qJD(5)) * pkin(4), 0, 0, -t79, 0, -t80, t181, t177 * t181 + t165 + t26, t214 + t377, pkin(4) * t263 + t177 * t79 + t352, t26 * t177 + (-t362 * t27 + (-t362 * t54 + t363 * t60) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, 0, -t330, 0, -pkin(10) * t329, pkin(10) * t330, 0, 0, 0, 0, -t111, 0, -t112, 0, t82, t81 (t111 * t363 + t262) * pkin(4) + t341 (-t362 * t81 + t363 * t82 + (-t116 * t362 + t117 * t363) * qJD(5)) * pkin(4), 0, 0, -t111, 0, -t112, 0, t50, t49, pkin(4) * t262 + t177 * t111 + t341, t50 * t177 + (-t362 * t49 + (-t362 * t96 + t363 * t97) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t175, 0, 0, 0, 0, 0, 0, 0, 0, t174, t175, 0, 0.2e1 * (-t177 * t362 + t318 * t363) * qJD(5) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, -t34, t113, t4, t3, 0, 0, 0, 0, -t33, 0, -t34, t113, 0.2e1 * t110 + t195, t2, t33 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, -t80, t181, t36, t35, 0, 0, 0, 0, -t79, 0, -t80, t181, 0.2e1 * t179 + t198, t27, t79 * pkin(5), t26 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, 0, -t112, 0, t82, t81, 0, 0, 0, 0, -t111, 0, -t112, 0, t50, t49, t111 * pkin(5), t50 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, -t276, 0, 0, 0, 0, 0, 0, 0, 0, -t275, -t276, 0, -pkin(5) * t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t79, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t111, 0, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
