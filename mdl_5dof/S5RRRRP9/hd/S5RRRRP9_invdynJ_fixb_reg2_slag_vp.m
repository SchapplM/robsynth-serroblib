% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:10
% EndTime: 2019-12-31 22:06:23
% DurationCPUTime: 6.73s
% Computational Cost: add. (7381->567), mult. (16569->711), div. (0->0), fcn. (11206->10), ass. (0->273)
t220 = cos(qJ(2));
t316 = qJD(1) * t220;
t187 = -qJD(3) + t316;
t217 = sin(qJ(2));
t305 = t217 * qJDD(1);
t200 = pkin(6) * t305;
t307 = qJD(1) * qJD(2);
t285 = t220 * t307;
t346 = qJDD(2) * pkin(2);
t120 = pkin(6) * t285 + t200 - t346;
t218 = sin(qJ(1));
t221 = cos(qJ(1));
t268 = g(1) * t221 + g(2) * t218;
t368 = g(3) * t220;
t238 = t217 * t268 - t368;
t232 = -t120 + t238;
t402 = -qJD(3) * pkin(7) * t187 - t232;
t171 = -qJD(4) + t187;
t216 = sin(qJ(3));
t317 = qJD(1) * t217;
t289 = t216 * t317;
t219 = cos(qJ(3));
t313 = qJD(2) * t219;
t140 = -t289 + t313;
t294 = t219 * t317;
t315 = qJD(2) * t216;
t141 = t294 + t315;
t215 = sin(qJ(4));
t376 = cos(qJ(4));
t250 = t215 * t140 + t376 * t141;
t310 = qJD(3) * t217;
t284 = qJD(1) * t310;
t393 = qJD(2) * qJD(3) + t285 + t305;
t280 = t393 * t216 + t219 * t284;
t245 = t219 * qJDD(2) - t280;
t73 = (-qJDD(2) + t284) * t216 - t393 * t219;
t27 = qJD(4) * t250 - t215 * t73 - t376 * t245;
t397 = t171 * t250 + t27;
t205 = t220 * qJDD(1);
t386 = -t217 * t307 + t205;
t137 = qJDD(3) - t386;
t133 = qJDD(4) + t137;
t214 = qJ(3) + qJ(4);
t206 = sin(t214);
t377 = pkin(8) + pkin(7);
t160 = t377 * t216;
t161 = t377 * t219;
t97 = -t215 * t160 + t376 * t161;
t401 = -t133 * t97 - t206 * t238;
t158 = -qJD(2) * pkin(2) + pkin(6) * t317;
t100 = -pkin(3) * t140 + t158;
t79 = -t376 * t140 + t141 * t215;
t32 = pkin(4) * t79 - qJ(5) * t250 + t100;
t400 = t32 * t79;
t399 = t100 * t79;
t398 = t79 * t250;
t337 = t215 * t219;
t143 = t376 * t216 + t337;
t385 = qJD(3) + qJD(4);
t93 = t385 * t143;
t350 = t143 * t316 - t93;
t288 = t376 * qJD(3);
t278 = t219 * t288;
t287 = t376 * qJD(4);
t295 = t216 * t316;
t296 = t376 * t219;
t338 = t215 * t216;
t349 = -t215 * t295 - t219 * t287 + t296 * t316 + t385 * t338 - t278;
t199 = pkin(3) * t219 + pkin(2);
t163 = t220 * t199;
t330 = t217 * t377;
t396 = t163 + t330;
t312 = qJD(2) * t220;
t293 = t216 * t312;
t309 = qJD(3) * t219;
t394 = t217 * t309 + t293;
t378 = t250 ^ 2;
t302 = t79 ^ 2 - t378;
t308 = qJD(4) * t215;
t26 = -t140 * t287 + t141 * t308 - t215 * t245 + t376 * t73;
t392 = -t171 * t79 - t26;
t45 = pkin(4) * t250 + qJ(5) * t79;
t297 = qJD(3) * t377;
t146 = t216 * t297;
t326 = t219 * t220;
t257 = pkin(3) * t217 - pkin(8) * t326;
t271 = pkin(2) * t217 - pkin(7) * t220;
t145 = t271 * qJD(1);
t98 = pkin(6) * t289 + t219 * t145;
t72 = qJD(1) * t257 + t98;
t128 = t216 * t145;
t332 = t217 * t219;
t334 = t216 * t220;
t86 = t128 + (-pkin(6) * t332 - pkin(8) * t334) * qJD(1);
t361 = -t97 * qJD(4) - t377 * t278 - t376 * t72 + (t146 + t86) * t215;
t119 = t386 * pkin(6) + qJDD(2) * pkin(7);
t272 = pkin(2) * t220 + pkin(7) * t217;
t153 = -pkin(1) - t272;
t134 = t153 * qJD(1);
t202 = pkin(6) * t316;
t159 = qJD(2) * pkin(7) + t202;
t311 = qJD(3) * t216;
t147 = t271 * qJD(2);
t94 = qJD(1) * t147 + qJDD(1) * t153;
t34 = t219 * t119 + t134 * t309 - t159 * t311 + t216 * t94;
t89 = t219 * t134 - t159 * t216;
t389 = t187 * t89 + t34;
t118 = t133 * qJ(5);
t152 = t171 * qJD(5);
t387 = t118 - t152;
t273 = -t202 + (-t295 + t311) * pkin(3);
t191 = pkin(6) * t326;
t103 = t216 * t153 + t191;
t121 = t133 * pkin(4);
t384 = t121 - qJDD(5);
t207 = cos(t214);
t327 = t218 * t220;
t110 = t206 * t327 + t207 * t221;
t323 = t221 * t206;
t329 = t218 * t207;
t112 = t220 * t323 - t329;
t85 = t219 * t94;
t90 = t134 * t216 + t159 * t219;
t35 = -qJD(3) * t90 - t216 * t119 + t85;
t19 = t137 * pkin(3) + t73 * pkin(8) + t35;
t23 = pkin(8) * t245 + t34;
t66 = -pkin(8) * t141 + t89;
t59 = -pkin(3) * t187 + t66;
t67 = pkin(8) * t140 + t90;
t282 = -t376 * t19 + t215 * t23 + t67 * t287 + t59 * t308;
t340 = t206 * t217;
t235 = g(1) * t112 + g(2) * t110 + g(3) * t340 - t282;
t226 = t250 * t32 - t235 - t384;
t383 = -t100 * t250 + t235;
t142 = -t296 + t338;
t382 = t133 * t142 + t171 * t350 - t317 * t79;
t369 = g(3) * t217;
t381 = t220 * t268 + t369;
t139 = t219 * t153;
t375 = pkin(6) * t216;
t87 = -pkin(8) * t332 + t139 + (-pkin(3) - t375) * t220;
t336 = t216 * t217;
t95 = -pkin(8) * t336 + t103;
t359 = t215 * t87 + t376 * t95;
t314 = qJD(2) * t217;
t321 = t219 * t147 + t314 * t375;
t46 = t257 * qJD(2) + (-t191 + (pkin(8) * t217 - t153) * t216) * qJD(3) + t321;
t64 = t216 * t147 + t153 * t309 + (-t217 * t313 - t220 * t311) * pkin(6);
t50 = -pkin(8) * t394 + t64;
t12 = -qJD(4) * t359 - t215 * t50 + t376 * t46;
t380 = t26 * t142 - t143 * t27 + t250 * t350 + t349 * t79;
t379 = -0.2e1 * pkin(1);
t374 = g(1) * t218;
t325 = t219 * t221;
t124 = t216 * t327 + t325;
t372 = g(2) * t124;
t370 = g(2) * t221;
t364 = -t350 * pkin(4) + t349 * qJ(5) - qJD(5) * t143 + t273;
t41 = t215 * t72 + t376 * t86;
t38 = qJ(5) * t317 + t41;
t249 = -t376 * t160 - t215 * t161;
t56 = t249 * qJD(4) - t376 * t146 - t297 * t337;
t363 = t56 - t38;
t362 = t56 - t41;
t360 = pkin(4) * t317 - t361;
t301 = t376 * t67;
t29 = t215 * t59 + t301;
t357 = t171 * t29;
t354 = t215 * t67;
t353 = t73 * t216;
t351 = t90 * t187;
t31 = t376 * t66 - t354;
t348 = pkin(3) * t287 + qJD(5) - t31;
t347 = pkin(6) * qJDD(1);
t345 = t140 * t187;
t344 = t140 * t216;
t343 = t141 * t140;
t342 = t141 * t187;
t341 = t141 * t219;
t339 = t207 * t217;
t335 = t216 * t218;
t333 = t216 * t221;
t331 = t217 * t221;
t328 = t218 * t219;
t324 = t220 * t221;
t28 = t376 * t59 - t354;
t322 = qJD(5) - t28;
t188 = pkin(3) * t336;
t148 = t217 * pkin(6) + t188;
t320 = t221 * pkin(1) + t218 * pkin(6);
t212 = t217 ^ 2;
t213 = t220 ^ 2;
t319 = t212 - t213;
t318 = t212 + t213;
t300 = t215 * t336;
t299 = t216 * t324;
t224 = qJD(1) ^ 2;
t298 = t217 * t224 * t220;
t101 = t394 * pkin(3) + pkin(6) * t312;
t292 = t216 * t310;
t290 = t171 * t317;
t3 = t215 * t19 + t376 * t23 + t59 * t287 - t67 * t308;
t113 = t206 * t218 + t207 * t324;
t281 = -t112 * pkin(4) + qJ(5) * t113;
t279 = t376 * t312;
t277 = t217 * t285;
t30 = t215 * t66 + t301;
t276 = pkin(3) * t308 - t30;
t275 = -pkin(4) * t340 + qJ(5) * t339;
t193 = t217 * t374;
t274 = -g(2) * t331 + t193;
t270 = g(1) * t110 - g(2) * t112;
t111 = t207 * t327 - t323;
t269 = g(1) * t111 - g(2) * t113;
t115 = t143 * t217;
t52 = t216 * t279 - t215 * t292 - qJD(4) * t300 + (t215 * t312 + (t288 + t287) * t217) * t219;
t266 = t115 * t27 + t52 * t79;
t265 = g(2) * (-t110 * pkin(4) + t111 * qJ(5));
t264 = pkin(4) * t207 + qJ(5) * t206;
t263 = pkin(6) * t140 - t158 * t216;
t262 = pkin(6) * t141 + t158 * t219;
t260 = -t216 * t90 - t219 * t89;
t259 = -pkin(7) * t137 + qJD(3) * t158;
t47 = -t215 * t95 + t376 * t87;
t253 = -pkin(6) * qJDD(2) + t307 * t379;
t11 = t215 * t46 + t87 * t287 - t308 * t95 + t376 * t50;
t251 = t27 * t142 - t350 * t79;
t248 = t216 * t137 - t187 * t309;
t247 = t219 * t137 + t187 * t311;
t246 = pkin(1) * t224 + t268;
t223 = qJD(2) ^ 2;
t244 = pkin(6) * t223 + qJDD(1) * t379 + t370;
t242 = pkin(3) * t335 + t199 * t324 + t221 * t330 + t320;
t241 = g(2) * t217 * t329 + t133 * t249 + (g(1) * t331 - t368) * t207;
t240 = t245 * t219;
t210 = t221 * pkin(6);
t239 = pkin(3) * t333 + t210 + (-pkin(1) - t396) * t218;
t234 = g(1) * t113 + g(2) * t111 + g(3) * t339 - t3;
t116 = t217 * t296 - t300;
t51 = t215 * t293 + t217 * t93 - t219 * t279;
t233 = t115 * t26 - t116 * t27 - t250 * t52 + t51 * t79;
t231 = t115 * t133 - t171 * t52 - t220 * t27 + t314 * t79;
t229 = -t171 * t28 + t234;
t227 = t249 * t26 - t97 * t27 - t56 * t79 - t381;
t58 = -pkin(3) * t245 + t120;
t198 = -t376 * pkin(3) - pkin(4);
t194 = pkin(3) * t215 + qJ(5);
t190 = pkin(3) * t328;
t127 = t219 * t324 + t335;
t126 = -t299 + t328;
t125 = -t218 * t326 + t333;
t102 = -pkin(6) * t334 + t139;
t99 = -pkin(6) * t294 + t128;
t88 = -t133 * t220 - t171 * t314;
t77 = pkin(4) * t142 - qJ(5) * t143 - t199;
t65 = -t103 * qJD(3) + t321;
t60 = pkin(4) * t115 - qJ(5) * t116 + t148;
t44 = t220 * pkin(4) - t47;
t43 = -qJ(5) * t220 + t359;
t37 = pkin(3) * t141 + t45;
t25 = -t171 * qJ(5) + t29;
t24 = t171 * pkin(4) + t322;
t20 = t143 * t133 + t349 * t171 - t250 * t317;
t13 = pkin(4) * t52 + qJ(5) * t51 - qJD(5) * t116 + t101;
t10 = -pkin(4) * t314 - t12;
t9 = qJ(5) * t314 - qJD(5) * t220 + t11;
t8 = -t116 * t26 - t250 * t51;
t7 = -t26 * t143 - t250 * t349;
t6 = t116 * t133 + t171 * t51 + t220 * t26 + t250 * t314;
t5 = t27 * pkin(4) + t26 * qJ(5) - qJD(5) * t250 + t58;
t2 = t282 - t384;
t1 = t3 + t387;
t4 = [0, 0, 0, 0, 0, qJDD(1), -t370 + t374, t268, 0, 0, qJDD(1) * t212 + 0.2e1 * t277, 0.2e1 * t205 * t217 - 0.2e1 * t307 * t319, qJDD(2) * t217 + t220 * t223, qJDD(1) * t213 - 0.2e1 * t277, qJDD(2) * t220 - t217 * t223, 0, t253 * t217 + (-t244 + t374) * t220, t217 * t244 + t220 * t253 - t193, 0.2e1 * t318 * t347 - t268, -g(1) * (-pkin(1) * t218 + t210) - g(2) * t320 + (pkin(6) ^ 2 * t318 + pkin(1) ^ 2) * qJDD(1), -t73 * t332 + (t219 * t312 - t292) * t141, (t140 * t219 - t141 * t216) * t312 + (t240 + t353 + (-t341 - t344) * qJD(3)) * t217, (-t187 * t313 + t73) * t220 + (qJD(2) * t141 + t247) * t217, -t140 * t394 - t245 * t336, (t187 * t315 - t245) * t220 + (qJD(2) * t140 - t248) * t217, -t137 * t220 - t187 * t314, -g(1) * t125 - g(2) * t127 + t102 * t137 - t65 * t187 + (-qJD(2) * t263 - t35) * t220 + (-pkin(6) * t245 + t89 * qJD(2) + t120 * t216 + t158 * t309) * t217, -g(1) * t124 - g(2) * t126 - t103 * t137 + t64 * t187 + (qJD(2) * t262 + t34) * t220 + (-pkin(6) * t73 - qJD(2) * t90 + t120 * t219 - t158 * t311) * t217, -t65 * t141 + t102 * t73 + t64 * t140 + t103 * t245 + t193 + t260 * t312 + (-t370 - t34 * t216 - t35 * t219 + (t216 * t89 - t219 * t90) * qJD(3)) * t217, t34 * t103 + t90 * t64 + t35 * t102 + t89 * t65 - g(1) * t210 - g(2) * (t221 * t272 + t320) - t153 * t374 + (t120 * t217 + t158 * t312) * pkin(6), t8, t233, t6, t266, -t231, t88, t100 * t52 + t101 * t79 + t115 * t58 - t12 * t171 + t133 * t47 + t148 * t27 + t220 * t282 + t28 * t314 + t269, -t100 * t51 + t101 * t250 + t11 * t171 + t116 * t58 - t133 * t359 - t148 * t26 + t220 * t3 - t29 * t314 - t270, -t11 * t79 - t115 * t3 + t116 * t282 - t12 * t250 + t26 * t47 - t27 * t359 + t28 * t51 - t29 * t52 + t274, -g(1) * t239 - g(2) * t242 + t100 * t101 + t29 * t11 + t28 * t12 + t58 * t148 - t282 * t47 + t3 * t359, t8, t6, -t233, t88, t231, t266, t10 * t171 + t115 * t5 + t13 * t79 - t133 * t44 + t2 * t220 - t24 * t314 + t27 * t60 + t32 * t52 + t269, -t1 * t115 + t10 * t250 + t116 * t2 - t24 * t51 - t25 * t52 - t26 * t44 - t27 * t43 - t79 * t9 + t274, -t1 * t220 - t116 * t5 - t13 * t250 + t133 * t43 - t171 * t9 + t25 * t314 + t26 * t60 + t32 * t51 + t270, t1 * t43 + t25 * t9 + t5 * t60 + t32 * t13 + t2 * t44 + t24 * t10 - g(1) * (-t111 * pkin(4) - t110 * qJ(5) + t239) - g(2) * (pkin(4) * t113 + qJ(5) * t112 + t242); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t298, t319 * t224, t305, t298, t205, qJDD(2), t217 * t246 - t200 - t368, t369 + (t246 - t347) * t220, 0, 0, -t187 * t341 - t353, (-t73 - t345) * t219 + (t245 + t342) * t216, (-t141 * t217 + t187 * t326) * qJD(1) + t248, t187 * t344 + t240, (-t140 * t217 - t187 * t334) * qJD(1) + t247, t187 * t317, -pkin(2) * t280 + t98 * t187 + t259 * t216 + (-t89 * t217 + t220 * t263) * qJD(1) + (t346 - t402) * t219, pkin(2) * t73 - t99 * t187 + t259 * t219 + (t217 * t90 - t220 * t262) * qJD(1) + t402 * t216, -t99 * t140 + t98 * t141 + ((qJD(3) * t141 + t245) * pkin(7) + t389) * t219 + (-t35 + t351 + (-qJD(3) * t140 - t73) * pkin(7)) * t216 - t381, -t158 * t202 - t89 * t98 - t90 * t99 + t232 * pkin(2) + (qJD(3) * t260 - t35 * t216 + t34 * t219 - t381) * pkin(7), t7, t380, t20, t251, -t382, t290, -t100 * t350 + t142 * t58 - t171 * t361 - t199 * t27 + t273 * t79 - t28 * t317 + t241, -t100 * t349 + t143 * t58 + t171 * t362 + t199 * t26 + t250 * t273 + t29 * t317 + t401, -t142 * t3 + t143 * t282 - t250 * t361 + t28 * t349 + t29 * t350 + t41 * t79 + t227, t3 * t97 - t282 * t249 - t58 * t199 - g(3) * t396 + t362 * t29 + t361 * t28 + t273 * t100 + t268 * (t199 * t217 - t220 * t377), t7, t20, -t380, t290, t382, t251, t142 * t5 + t171 * t360 + t24 * t317 + t27 * t77 - t32 * t350 + t364 * t79 + t241, -t1 * t142 + t143 * t2 - t24 * t349 + t25 * t350 + t250 * t360 + t38 * t79 + t227, -t143 * t5 - t171 * t363 - t25 * t317 - t250 * t364 + t26 * t77 + t32 * t349 - t401, -g(3) * t163 + t1 * t97 - t2 * t249 + t5 * t77 + t364 * t32 + t363 * t25 + t360 * t24 + (-g(3) * t264 - t268 * t377) * t220 + (-g(3) * t377 + t268 * (t199 + t264)) * t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t343, -t140 ^ 2 + t141 ^ 2, -t73 + t345, t343, t245 - t342, t137, -t159 * t309 - g(1) * t126 + t372 - t158 * t141 - t351 + t85 + (-qJD(3) * t134 - t119 + t369) * t216, g(1) * t127 - g(2) * t125 + g(3) * t332 - t140 * t158 - t389, 0, 0, t398, -t302, t392, -t398, -t397, t133, -t30 * t171 + (t376 * t133 - t141 * t79 + t171 * t308) * pkin(3) + t383, t399 - t31 * t171 + (-t133 * t215 - t141 * t250 + t171 * t287) * pkin(3) + t234, -t28 * t79 + t29 * t250 - t30 * t250 + t31 * t79 + (t376 * t26 - t215 * t27 + (t215 * t250 - t376 * t79) * qJD(4)) * pkin(3), -g(1) * t190 + t28 * t30 - t29 * t31 + (g(2) * t325 - t282 * t376 - t100 * t141 + t3 * t215 + t381 * t216 + (-t28 * t215 + t29 * t376) * qJD(4)) * pkin(3), t398, t392, t302, t133, t397, -t398, -t133 * t198 + t171 * t276 - t37 * t79 - t226, -t194 * t27 - t198 * t26 + (t25 + t276) * t250 + (t24 - t348) * t79, t133 * t194 - t171 * t348 + t250 * t37 - t234 + t387 - t400, t1 * t194 + t2 * t198 - t32 * t37 - t24 * t30 - g(1) * (t190 + t281) - t265 - g(3) * (-t188 + t275) + t348 * t25 + (g(1) * t299 + t24 * t308 + t372) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t398, -t302, t392, -t398, -t397, t133, -t357 + t383, t229 + t399, 0, 0, t398, t392, t302, t133, t397, -t398, -t45 * t79 + t121 - t226 - t357, pkin(4) * t26 - t27 * qJ(5) + (t25 - t29) * t250 + (t24 - t322) * t79, t250 * t45 + 0.2e1 * t118 - 0.2e1 * t152 - t229 - t400, -t2 * pkin(4) - g(1) * t281 - g(3) * t275 + t1 * qJ(5) - t24 * t29 + t25 * t322 - t32 * t45 - t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 + t398, t392, -t171 ^ 2 - t378, t171 * t25 + t226;];
tau_reg = t4;
