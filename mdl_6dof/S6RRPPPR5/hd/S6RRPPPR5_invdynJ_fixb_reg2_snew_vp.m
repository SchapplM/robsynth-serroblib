% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:54:18
% EndTime: 2019-05-06 08:54:37
% DurationCPUTime: 7.11s
% Computational Cost: add. (12982->434), mult. (29139->518), div. (0->0), fcn. (18764->8), ass. (0->255)
t235 = sin(qJ(2));
t238 = cos(qJ(2));
t232 = sin(pkin(9));
t233 = cos(pkin(9));
t288 = qJD(1) * t235;
t203 = qJD(2) * t232 + t233 * t288;
t347 = t203 ^ 2;
t201 = -t233 * qJD(2) + t232 * t288;
t348 = t201 ^ 2;
t355 = -t348 - t347;
t287 = qJD(1) * t238;
t180 = t201 * t287;
t225 = t235 * qJDD(1);
t280 = qJD(1) * qJD(2);
t275 = t238 * t280;
t209 = t225 + t275;
t254 = qJDD(2) * t232 + t209 * t233;
t132 = t254 - t180;
t181 = t203 * t287;
t263 = qJDD(2) * t233 - t232 * t209;
t368 = t263 - t181;
t391 = t132 * t232 + t368 * t233;
t396 = pkin(7) * (t355 * t235 + t238 * t391);
t77 = t132 * t233 - t368 * t232;
t416 = -pkin(1) * t77 - t396;
t363 = t180 + t254;
t230 = t238 ^ 2;
t240 = qJD(1) ^ 2;
t294 = t230 * t240;
t356 = t347 + t294;
t221 = t235 * t280;
t279 = t238 * qJDD(1);
t210 = -t221 + t279;
t297 = t203 * t201;
t357 = -t297 + t210;
t382 = t357 * t232;
t388 = t356 * t233 - t382;
t381 = t357 * t233;
t390 = -t356 * t232 - t381;
t394 = pkin(7) * (-t363 * t235 + t238 * t390) - pkin(1) * t388;
t361 = t181 + t263;
t147 = t294 + t348;
t358 = t210 + t297;
t370 = t358 * t233;
t89 = t147 * t232 + t370;
t371 = t358 * t232;
t92 = t233 * t147 - t371;
t395 = -pkin(7) * (t361 * t235 + t238 * t92) + pkin(1) * t89;
t413 = pkin(2) * t89;
t411 = pkin(2) * t388;
t410 = qJ(3) * t77;
t409 = qJ(3) * t89;
t408 = qJ(3) * t92;
t407 = qJ(3) * t388;
t177 = t347 - t294;
t406 = t238 * t132 - t235 * (t177 * t232 - t370);
t405 = -qJ(4) * t147 - t413;
t404 = pkin(2) * t77 - qJ(4) * t368;
t403 = qJ(4) * t357 - t411;
t402 = pkin(2) * t361 - t408;
t400 = pkin(2) * t363 + qJ(3) * t390;
t399 = -pkin(2) * t355 + qJ(3) * t391;
t176 = t348 - t294;
t385 = t235 * (t176 * t233 + t382) - t238 * t368;
t380 = t177 * t233 + t371;
t389 = t176 * t232 - t381;
t323 = t363 * t232;
t369 = t361 * t233;
t386 = t235 * (t323 - t369) + t238 * (-t348 + t347);
t375 = pkin(3) * t361;
t384 = pkin(4) * t132;
t322 = t363 * t233;
t378 = t361 * t232 + t322;
t234 = sin(qJ(6));
t237 = cos(qJ(6));
t144 = t201 * t234 - t237 * t203;
t216 = -qJD(6) + t287;
t123 = t144 * t216;
t84 = -t144 * qJD(6) + t234 * t254 - t237 * t263;
t367 = t123 + t84;
t236 = sin(qJ(1));
t239 = cos(qJ(1));
t261 = g(1) * t239 + g(2) * t236;
t331 = qJDD(1) * pkin(7);
t191 = -pkin(1) * t240 - t261 + t331;
t258 = -pkin(2) * t238 - qJ(3) * t235;
t207 = t258 * qJD(1);
t264 = qJD(1) * t207 + t191;
t364 = t264 * t235;
t286 = qJD(3) * t201;
t185 = -0.2e1 * t286;
t274 = qJD(4) * t287;
t360 = t185 - 0.2e1 * t274;
t359 = 0.2e1 * t286 + 0.2e1 * t274;
t353 = -t254 * qJ(4) - t375;
t166 = pkin(4) * t203 + qJ(5) * t287;
t193 = t348 * pkin(4);
t283 = qJD(5) * t201;
t352 = t203 * t166 + t193 - 0.2e1 * t283;
t305 = t263 * qJ(5);
t351 = (-0.2e1 * qJD(4) - t166) * t203 - t193 - t305;
t146 = t201 * t237 + t203 * t234;
t268 = -t234 * t263 - t237 * t254;
t66 = (qJD(6) + t216) * t146 + t268;
t171 = t233 * t180;
t276 = t238 * t297;
t260 = t235 * (-t232 * t263 - t171) + t276;
t342 = t238 * g(3);
t346 = qJD(2) ^ 2;
t252 = -qJDD(2) * pkin(2) - t346 * qJ(3) + qJDD(3) + t342;
t248 = t252 + t353;
t293 = t235 * t191;
t241 = -(qJ(4) * t201 * t238 - t207 * t235) * qJD(1) + t248 + t293;
t349 = qJ(5) * t358 + t384;
t142 = t144 ^ 2;
t143 = t146 ^ 2;
t214 = t216 ^ 2;
t345 = pkin(4) + pkin(8);
t343 = t235 * g(3);
t341 = pkin(3) + qJ(5);
t340 = -pkin(5) - qJ(4);
t110 = t252 + t364;
t282 = -pkin(5) * t287 - pkin(8) * t201 + 0.2e1 * qJD(5);
t31 = -t254 * pkin(5) - t347 * pkin(8) + (-qJ(4) * t287 + t282) * t201 + t110 + t351 + t353;
t339 = t234 * t31;
t205 = -qJDD(6) + t210;
t307 = t146 * t144;
t87 = t205 - t307;
t338 = t234 * t87;
t88 = -t205 - t307;
t337 = t234 * t88;
t336 = t237 * t31;
t335 = t237 * t87;
t334 = t237 * t88;
t333 = qJ(4) * t355;
t330 = t110 * t232;
t329 = t110 * t233;
t296 = t216 * t234;
t295 = t216 * t237;
t215 = t238 * t240 * t235;
t292 = t235 * (qJDD(2) + t215);
t291 = t238 * (-t215 + qJDD(2));
t271 = t236 * g(1) - t239 * g(2);
t190 = qJDD(1) * pkin(1) + t240 * pkin(7) + t271;
t257 = t209 + t275;
t105 = -t257 * qJ(3) + (-t210 + t221) * pkin(2) - t190;
t111 = -pkin(2) * t346 + qJDD(2) * qJ(3) + t238 * t264 - t343;
t290 = t232 * t105 + t233 * t111;
t170 = t233 * t181;
t113 = t232 * t254 - t170;
t169 = t232 * t180;
t112 = t233 * t263 - t169;
t289 = t169 + t170;
t285 = qJD(3) * t203;
t284 = qJD(4) * t203;
t278 = -t143 - t214;
t277 = t238 * t307;
t72 = t185 + t290;
t273 = qJD(5) * t287;
t272 = -qJ(4) * t232 - pkin(2);
t270 = -t233 * t105 + t232 * t111;
t71 = t270 + 0.2e1 * t285;
t42 = t232 * t71 + t233 * t72;
t159 = t293 + t342;
t160 = t191 * t238 - t343;
t269 = t235 * t159 + t238 * t160;
t262 = t232 * t181;
t259 = t235 * (t233 * t254 + t262) - t276;
t148 = pkin(3) * t201 - qJ(4) * t203;
t249 = t210 * pkin(3) - qJ(4) * t294 + qJDD(4) + t270;
t47 = (0.2e1 * qJD(3) + t148) * t203 + t249;
t242 = t47 + t349;
t27 = -pkin(5) * t347 + pkin(8) * t254 + t282 * t287 + t242;
t255 = pkin(3) * t294 + t210 * qJ(4) + t201 * t148 - t290;
t246 = pkin(4) * t263 - qJ(5) * t348 - t166 * t287 + qJDD(5) - t255;
t36 = t246 + t360;
t28 = -pkin(5) * t357 + pkin(8) * t368 + t36;
t15 = t234 * t27 - t237 * t28;
t16 = t234 * t28 + t237 * t27;
t7 = -t237 * t15 + t234 * t16;
t8 = t234 * t15 + t237 * t16;
t256 = t232 * t72 - t233 * t71;
t253 = -pkin(1) + t258;
t46 = -t255 + t360;
t192 = t238 * t210;
t247 = t235 * (t171 - t262) + t192;
t245 = -t203 * t148 - t249 - 0.2e1 * t285;
t243 = -t246 + t359;
t35 = t242 + 0.2e1 * t273;
t184 = 0.2e1 * t284;
t48 = t184 - t364 + (t363 + t180) * qJ(4) - t248;
t52 = t241 - 0.2e1 * t284;
t229 = t235 ^ 2;
t226 = t229 * t240;
t211 = -0.2e1 * t221 + t279;
t208 = t225 + 0.2e1 * t275;
t116 = -t143 + t214;
t115 = t142 - t214;
t96 = t143 - t142;
t95 = -t214 - t142;
t86 = pkin(4) * t358 + qJ(4) * t361;
t83 = -qJD(6) * t146 - t268;
t76 = (-t144 * t237 + t146 * t234) * t216;
t75 = (t144 * t234 + t146 * t237) * t216;
t74 = -t142 - t143;
t73 = -pkin(4) * t357 + t341 * t363;
t70 = -t123 + t84;
t65 = (qJD(6) - t216) * t146 + t268;
t62 = -t115 * t237 - t338;
t61 = t116 * t234 - t334;
t60 = t115 * t234 - t335;
t59 = t116 * t237 + t337;
t58 = -t146 * t296 - t237 * t84;
t57 = -t146 * t295 + t234 * t84;
t56 = t144 * t295 + t234 * t83;
t55 = -t144 * t296 + t237 * t83;
t54 = -t234 * t278 + t335;
t53 = t237 * t278 + t338;
t51 = t237 * t95 - t337;
t50 = t234 * t95 + t334;
t49 = t52 - t375;
t45 = t47 - t333;
t44 = t241 + 0.2e1 * t283 + t351;
t43 = -pkin(3) * t355 + t46;
t40 = t234 * t70 - t237 * t66;
t39 = t234 * t367 + t237 * t65;
t38 = -t234 * t66 - t237 * t70;
t37 = -t234 * t65 + t237 * t367;
t34 = t232 * t54 + t233 * t53;
t33 = t232 * t53 - t233 * t54;
t32 = -pkin(4) * t356 + t305 + t352 + t48;
t30 = t232 * t51 + t233 * t50;
t29 = t232 * t50 - t233 * t51;
t26 = t184 + (t361 + t263) * qJ(5) + t375 - pkin(4) * t147 - t241 + t352;
t25 = t245 - 0.2e1 * t273 + t333 - t349 - t384;
t24 = -pkin(4) * t368 + t341 * t355 + t243;
t23 = t232 * t47 + t233 * t46;
t22 = t232 * t46 - t233 * t47;
t21 = pkin(4) * t35 - qJ(4) * t44;
t20 = t232 * t40 + t233 * t38;
t19 = t232 * t38 - t233 * t40;
t18 = t232 * t35 + t233 * t36;
t17 = t232 * t36 - t233 * t35;
t13 = pkin(4) * t36 - t341 * t44;
t12 = -t341 * t367 + t345 * t53 - t336;
t11 = t340 * t367 + t345 * t54 + t339;
t10 = -t341 * t65 + t345 * t50 - t339;
t9 = t340 * t65 + t345 * t51 - t336;
t6 = -t341 * t74 + t345 * t38 + t7;
t5 = t340 * t74 + t345 * t40 + t8;
t4 = t232 * t8 + t233 * t7;
t3 = t232 * t7 - t233 * t8;
t2 = -t31 * t341 + t345 * t7;
t1 = t31 * t340 + t345 * t8;
t14 = [0, 0, 0, 0, 0, qJDD(1), t271, t261, 0, 0, t257 * t235, t208 * t238 + t211 * t235, t292 + t238 * (-t226 + t346), -t221 * t238 + t192, t235 * (t294 - t346) + t291, 0, t238 * t190 + pkin(1) * t211 + pkin(7) * (t238 * (-t294 - t346) - t292), -t235 * t190 - pkin(1) * t208 + pkin(7) * (-t291 - t235 * (-t226 - t346)), pkin(1) * (t226 + t294) + (t229 + t230) * t331 + t269, pkin(1) * t190 + pkin(7) * t269, t259, -t386, -t406, t260, t385, t192 + t235 * (t201 * t233 - t203 * t232) * t287, t235 * (t330 + t409) + t238 * (t71 + t413) + t395, t235 * (t329 + t407) + t238 * (t72 + t411) - t394, -t235 * t256 - t253 * t77 + t396, pkin(7) * (t110 * t235 + t238 * t42) + t253 * t256, t247, t406, -t385, t259, -t386, t260, t235 * (-t232 * t43 + t233 * t45 + t410) + t238 * (pkin(3) * t132 + t404) - t416, t235 * (-qJ(4) * t369 - t232 * t49 - t409) + t238 * (-pkin(3) * t358 + t245 + t405) - t395, t235 * (-pkin(3) * t323 + t233 * t48 - t407) + t238 * (-pkin(3) * t356 + t255 + t359 + t403) + t394, t235 * (-qJ(3) * t22 + (pkin(3) * t232 - qJ(4) * t233) * t52) + t238 * (-pkin(2) * t22 + pkin(3) * t47 - qJ(4) * t46) - pkin(1) * t22 + pkin(7) * (t23 * t238 + t235 * t52), t260, t385, t386, t247, t406, t259, t235 * (-t232 * t73 + t233 * t32 - t407) + t238 * (-t341 * t356 + t243 + t403) + t394, t235 * (-t232 * t24 + t233 * t25 - t410) + t238 * (-t132 * t341 - t404) + t416, t235 * (-t232 * t26 + t233 * t86 + t409) + t238 * (t341 * t358 + t35 - t405) + t395, t235 * (-qJ(3) * t17 - t13 * t232 + t21 * t233) + t238 * (-pkin(2) * t17 - qJ(4) * t36 + t341 * t35) - pkin(1) * t17 + pkin(7) * (t18 * t238 + t235 * t44), t235 * (-t232 * t58 + t233 * t57) - t277, t235 * (-t232 * t39 + t233 * t37) - t238 * t96, t235 * (-t232 * t61 + t233 * t59) - t238 * t70, t235 * (-t232 * t56 + t233 * t55) + t277, t235 * (-t232 * t62 + t233 * t60) + t238 * t66, t235 * (-t232 * t76 + t233 * t75) + t238 * t205, t235 * (-qJ(3) * t29 - t10 * t232 + t233 * t9) + t238 * (-pkin(2) * t29 + t340 * t50 + t341 * t51 + t15) - pkin(1) * t29 + pkin(7) * (t235 * t65 + t238 * t30), t235 * (-qJ(3) * t33 + t11 * t233 - t12 * t232) + t238 * (-pkin(2) * t33 + t340 * t53 + t341 * t54 + t16) - pkin(1) * t33 + pkin(7) * (t235 * t367 + t238 * t34), t235 * (-qJ(3) * t19 - t232 * t6 + t233 * t5) + t238 * (-pkin(2) * t19 + t340 * t38 + t341 * t40) - pkin(1) * t19 + pkin(7) * (t20 * t238 + t235 * t74), t235 * (-qJ(3) * t3 + t1 * t233 - t2 * t232) + t238 * (-pkin(2) * t3 + t340 * t7 + t341 * t8) - pkin(1) * t3 + pkin(7) * (t235 * t31 + t238 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t226 - t294, t225, t215, t279, qJDD(2), -t159, -t160, 0, 0, t113, t378, -t380, t112, t389, t289, -t329 + t402, t330 - t400, t399 + t42, -pkin(2) * t110 + qJ(3) * t42, t289, t380, -t389, t113, t378, t112, t232 * t45 + t233 * t43 + t399, t233 * t49 + t272 * t361 + t408, pkin(3) * t322 + t232 * t48 + t400, qJ(3) * t23 + (-pkin(3) * t233 + t272) * t52, t112, t389, -t378, t289, t380, t113, t232 * t32 + t233 * t73 + t400, t232 * t25 + t233 * t24 - t399, t232 * t86 + t233 * t26 + t402, -pkin(2) * t44 + qJ(3) * t18 + t13 * t233 + t21 * t232, t232 * t57 + t233 * t58, t232 * t37 + t233 * t39, t232 * t59 + t233 * t61, t232 * t55 + t233 * t56, t232 * t60 + t233 * t62, t232 * t75 + t233 * t76, -pkin(2) * t65 + qJ(3) * t30 + t10 * t233 + t232 * t9, -pkin(2) * t367 + qJ(3) * t34 + t11 * t232 + t12 * t233, -pkin(2) * t74 + qJ(3) * t20 + t232 * t5 + t233 * t6, -pkin(2) * t31 + qJ(3) * t4 + t1 * t232 + t2 * t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t361, t363, t355, t110, 0, 0, 0, 0, 0, 0, t355, t361, -t363, t52, 0, 0, 0, 0, 0, 0, -t363, -t355, -t361, t44, 0, 0, 0, 0, 0, 0, t65, t367, t74, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t358, -t356, t47, 0, 0, 0, 0, 0, 0, -t356, -t132, t358, t35, 0, 0, 0, 0, 0, 0, t51, t54, t40, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t357, -t368, -t147, t36, 0, 0, 0, 0, 0, 0, t50, t53, t38, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, t96, t70, -t307, -t66, -t205, -t15, -t16, 0, 0;];
tauJ_reg  = t14;
