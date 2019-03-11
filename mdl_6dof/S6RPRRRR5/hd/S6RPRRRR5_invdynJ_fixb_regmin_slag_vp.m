% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:54
% EndTime: 2019-03-09 07:10:12
% DurationCPUTime: 8.02s
% Computational Cost: add. (10964->492), mult. (26662->635), div. (0->0), fcn. (22233->18), ass. (0->285)
t249 = cos(pkin(11));
t257 = cos(qJ(3));
t346 = t257 * t249;
t248 = sin(pkin(11));
t253 = sin(qJ(3));
t352 = t248 * t253;
t189 = -t346 + t352;
t177 = t189 * qJD(1);
t190 = t248 * t257 + t249 * t253;
t178 = t190 * qJD(1);
t252 = sin(qJ(4));
t393 = cos(qJ(4));
t144 = t393 * t177 + t178 * t252;
t401 = qJD(5) + qJD(6);
t451 = t144 + t401;
t417 = -qJD(5) - t144;
t139 = qJD(6) - t417;
t251 = sin(qJ(5));
t255 = cos(qJ(6));
t250 = sin(qJ(6));
t256 = cos(qJ(5));
t351 = t250 * t256;
t196 = t251 * t255 + t351;
t195 = t250 * t251 - t255 * t256;
t446 = t451 * t195;
t334 = qJD(1) * qJD(3);
t316 = t253 * t334;
t315 = t257 * t334;
t332 = t249 * qJDD(1);
t333 = t248 * qJDD(1);
t324 = t249 * t315 + t253 * t332 + t257 * t333;
t151 = -t248 * t316 + t324;
t325 = -t248 * t315 - t249 * t316 - t253 * t333;
t274 = t257 * t332 + t325;
t280 = -t252 * t177 + t178 * t393;
t75 = qJD(4) * t280 + t252 * t151 - t393 * t274;
t72 = qJDD(5) + t75;
t71 = qJDD(6) + t72;
t450 = -t139 * t446 + t196 * t71;
t449 = t451 * t196;
t340 = qJD(5) * t256;
t425 = t144 * t256;
t448 = t340 + t425;
t341 = qJD(5) * t251;
t426 = t144 * t251;
t447 = t341 + t426;
t245 = pkin(11) + qJ(3);
t238 = qJ(4) + t245;
t228 = sin(t238);
t254 = sin(qJ(1));
t258 = cos(qJ(1));
t295 = g(1) * t258 + g(2) * t254;
t445 = t295 * t228;
t444 = pkin(10) * t426;
t317 = qJD(4) * t393;
t385 = pkin(7) + qJ(2);
t210 = t385 * t248;
t191 = qJD(1) * t210;
t211 = t385 * t249;
t192 = qJD(1) * t211;
t285 = t191 * t253 - t192 * t257;
t131 = -pkin(8) * t177 - t285;
t125 = t252 * t131;
t404 = -t257 * t191 - t192 * t253;
t130 = -pkin(8) * t178 + t404;
t86 = t130 * t393 - t125;
t428 = pkin(3) * t317 - t86;
t297 = -t139 * t449 - t195 * t71;
t246 = qJD(3) + qJD(4);
t132 = -t256 * t246 + t251 * t280;
t134 = t246 * t251 + t256 * t280;
t88 = t255 * t132 + t134 * t250;
t380 = t280 * t88;
t443 = t297 + t380;
t442 = t447 * pkin(5);
t441 = pkin(5) * t280 + pkin(10) * t425;
t109 = pkin(4) * t280 + pkin(9) * t144;
t94 = pkin(3) * t178 + t109;
t440 = -t428 * t251 - t256 * t94;
t242 = qJDD(3) + qJDD(4);
t127 = qJD(3) * pkin(3) + t130;
t342 = qJD(4) * t252;
t335 = qJD(1) * qJD(2);
t395 = qJDD(1) * t385 + t335;
t166 = t395 * t248;
t167 = t395 * t249;
t306 = -t257 * t166 - t253 * t167;
t69 = qJDD(3) * pkin(3) - t151 * pkin(8) + qJD(3) * t285 + t306;
t286 = -t253 * t166 + t257 * t167;
t73 = t274 * pkin(8) + t404 * qJD(3) + t286;
t303 = -t127 * t342 - t131 * t317 - t252 * t73 + t393 * t69;
t19 = -pkin(4) * t242 - t303;
t74 = t393 * t151 - t177 * t317 - t178 * t342 + t252 * t274;
t309 = -t256 * t242 + t251 * t74;
t50 = qJD(5) * t134 + t309;
t10 = pkin(5) * t50 + t19;
t229 = cos(t238);
t247 = qJ(5) + qJ(6);
t240 = cos(t247);
t387 = g(3) * t240;
t82 = t127 * t393 - t125;
t76 = -t246 * pkin(4) - t82;
t54 = t132 * pkin(5) + t76;
t126 = t393 * t131;
t83 = t252 * t127 + t126;
t77 = pkin(9) * t246 + t83;
t230 = -pkin(2) * t249 - pkin(1);
t204 = qJD(1) * t230 + qJD(2);
t156 = t177 * pkin(3) + t204;
t84 = t144 * pkin(4) - pkin(9) * t280 + t156;
t35 = -t251 * t77 + t256 * t84;
t29 = -pkin(10) * t134 + t35;
t21 = -pkin(5) * t417 + t29;
t36 = t251 * t84 + t256 * t77;
t30 = -pkin(10) * t132 + t36;
t6 = t21 * t255 - t250 * t30;
t439 = t10 * t195 - t229 * t387 + t240 * t445 - t6 * t280 + t449 * t54;
t239 = sin(t247);
t388 = g(3) * t239;
t379 = t255 * t30;
t7 = t21 * t250 + t379;
t438 = t10 * t196 + t229 * t388 - t239 * t445 + t7 * t280 - t446 * t54;
t49 = t251 * t242 + t246 * t340 + t256 * t74 - t280 * t341;
t437 = -t132 * t448 - t251 * t50 + t49 * t256;
t47 = t49 * t251;
t436 = t134 * t448 + t47;
t287 = t132 * t250 - t255 * t134;
t378 = t287 * t280;
t435 = t378 + t450;
t338 = qJD(6) * t255;
t339 = qJD(6) * t250;
t14 = -t132 * t338 - t134 * t339 - t250 * t50 + t255 * t49;
t434 = t14 * t196 + t287 * t446;
t367 = t134 * t280;
t64 = t251 * t72;
t433 = -t417 * t448 - t367 + t64;
t267 = qJD(6) * t287 - t250 * t49 - t255 * t50;
t432 = -t14 * t195 + t196 * t267 + t287 * t449 + t446 * t88;
t431 = t144 * t76;
t430 = t287 * t88;
t389 = g(3) * t229;
t429 = t19 + t389;
t362 = t144 * t246;
t427 = t74 + t362;
t422 = t280 * t144;
t365 = t280 * t246;
t420 = -t75 + t365;
t416 = -t144 ^ 2 + t280 ^ 2;
t415 = t287 ^ 2 - t88 ^ 2;
t414 = t139 * t88 + t14;
t354 = t240 * t254;
t355 = t239 * t258;
t161 = -t229 * t354 + t355;
t353 = t240 * t258;
t356 = t239 * t254;
t163 = t229 * t353 + t356;
t25 = t30 * t339;
t413 = g(1) * t163 - g(2) * t161 + t228 * t387 + t54 * t88 + t25;
t160 = t229 * t356 + t353;
t162 = -t229 * t355 + t354;
t265 = t127 * t317 - t131 * t342 + t252 * t69 + t393 * t73;
t18 = t242 * pkin(9) + t265;
t135 = qJDD(2) - t325 * pkin(3) + (-pkin(1) + (-pkin(3) * t257 - pkin(2)) * t249) * qJDD(1);
t28 = t75 * pkin(4) - t74 * pkin(9) + t135;
t27 = t256 * t28;
t2 = t72 * pkin(5) - t49 * pkin(10) - qJD(5) * t36 - t251 * t18 + t27;
t282 = t256 * t18 + t251 * t28 + t84 * t340 - t341 * t77;
t3 = -pkin(10) * t50 + t282;
t323 = t255 * t2 - t250 * t3;
t412 = -g(1) * t162 + g(2) * t160 - qJD(6) * t7 + t228 * t388 + t54 * t287 + t323;
t411 = -t139 * t287 + t267;
t223 = g(3) * t228;
t410 = t156 * t144 + t229 * t295 + t223 - t265;
t85 = t252 * t130 + t126;
t300 = pkin(3) * t342 - t85;
t368 = t132 * t280;
t406 = t139 * t280;
t405 = t417 * t280;
t153 = -t252 * t189 + t190 * t393;
t110 = t196 * t153;
t345 = -t253 * t210 + t257 * t211;
t403 = t251 * t94 - t428 * t256;
t402 = qJ(2) * qJDD(1);
t66 = t76 * t341;
t398 = t256 * t445 - t35 * t280 + t66;
t397 = t429 * t251 + t36 * t280 + t76 * t340;
t396 = -t156 * t280 + t303 - t389 + t445;
t394 = -pkin(9) - pkin(10);
t180 = t190 * qJD(3);
t392 = pkin(3) * t180;
t386 = t256 * pkin(5);
t232 = pkin(3) * t252 + pkin(9);
t384 = -pkin(10) - t232;
t164 = pkin(3) * t189 + t230;
t279 = -t189 * t393 - t252 * t190;
t101 = -pkin(4) * t279 - pkin(9) * t153 + t164;
t187 = t257 * t210;
t305 = -t211 * t253 - t187;
t137 = -pkin(8) * t190 + t305;
t138 = -pkin(8) * t189 + t345;
t100 = t252 * t137 + t138 * t393;
t96 = t256 * t100;
t381 = t251 * t101 + t96;
t373 = t300 + t442;
t372 = t251 * t109 + t256 * t82;
t371 = qJDD(1) * pkin(1);
t179 = t189 * qJD(3);
t112 = qJD(4) * t279 - t179 * t393 - t252 * t180;
t370 = t112 * t251;
t369 = t112 * t256;
t366 = t134 * t251;
t359 = t153 * t251;
t358 = t153 * t256;
t350 = t251 * t254;
t349 = t251 * t258;
t348 = t254 * t256;
t347 = t256 * t258;
t344 = t248 ^ 2 + t249 ^ 2;
t343 = qJD(3) * t178;
t327 = qJD(5) * pkin(9) * t417;
t322 = qJD(5) * t394;
t321 = qJD(1) * t352;
t320 = t153 * t341;
t319 = t153 * t340;
t314 = qJD(6) * t21 + t3;
t312 = qJD(5) * t384;
t310 = -qJD(5) * t84 - t18;
t308 = t344 * qJD(1) ^ 2;
t304 = t417 * t251;
t301 = 0.2e1 * t344;
t233 = -pkin(3) * t393 - pkin(4);
t299 = t442 - t83;
t298 = -t340 * t77 + t27;
t296 = -pkin(9) * t72 + t431;
t294 = g(1) * t254 - g(2) * t258;
t184 = t384 * t251;
t293 = -qJD(6) * t184 - t251 * t312 + t403 + t444;
t241 = t256 * pkin(10);
t185 = t232 * t256 + t241;
t292 = qJD(6) * t185 - t256 * t312 - t440 + t441;
t213 = t394 * t251;
t291 = -qJD(6) * t213 - t251 * t322 + t372 + t444;
t107 = t256 * t109;
t214 = pkin(9) * t256 + t241;
t290 = qJD(6) * t214 - t251 * t82 - t256 * t322 + t107 + t441;
t289 = -t232 * t72 + t431;
t65 = t256 * t72;
t284 = t417 * t447 + t65;
t281 = t137 * t393 - t252 * t138;
t278 = t319 + t370;
t277 = -t320 + t369;
t269 = -qJD(3) * t187 + qJD(2) * t346 + (-qJD(2) * t248 - qJD(3) * t211) * t253;
t116 = -t180 * pkin(8) + t269;
t263 = -t190 * qJD(2) - t345 * qJD(3);
t117 = t179 * pkin(8) + t263;
t39 = qJD(4) * t281 + t116 * t393 + t252 * t117;
t113 = qJD(4) * t153 - t252 * t179 + t180 * t393;
t53 = pkin(4) * t113 - pkin(9) * t112 + t392;
t276 = -t100 * t341 + t101 * t340 + t251 * t53 + t256 * t39;
t272 = -t294 - t371;
t235 = qJDD(2) - t371;
t270 = -t235 - t272;
t266 = t301 * t335 - t295;
t40 = qJD(4) * t100 + t252 * t116 - t117 * t393;
t237 = cos(t245);
t236 = sin(t245);
t234 = -pkin(4) - t386;
t212 = t233 - t386;
t203 = qJDD(1) * t230 + qJDD(2);
t174 = t229 * t347 + t350;
t173 = -t229 * t349 + t348;
t172 = -t229 * t348 + t349;
t171 = t229 * t350 + t347;
t111 = t195 * t153;
t98 = t256 * t101;
t62 = pkin(5) * t359 - t281;
t52 = t256 * t53;
t38 = -pkin(10) * t359 + t381;
t32 = -pkin(5) * t279 - pkin(10) * t358 - t100 * t251 + t98;
t23 = t112 * t351 - t250 * t320 - t339 * t359 + (t401 * t358 + t370) * t255;
t22 = -t401 * t110 - t195 * t112;
t20 = pkin(5) * t278 + t40;
t5 = -pkin(10) * t278 + t276;
t4 = -pkin(10) * t369 + t113 * pkin(5) - t251 * t39 + t52 + (-t96 + (pkin(10) * t153 - t101) * t251) * qJD(5);
t1 = [qJDD(1), t294, t295, t270 * t249, -t270 * t248, t301 * t402 + t266 (-t235 + t294) * pkin(1) + (t344 * t402 + t266) * qJ(2), t151 * t190 - t178 * t179, -t151 * t189 + t179 * t177 - t178 * t180 + t190 * t274, -qJD(3) * t179 + qJDD(3) * t190, -qJD(3) * t180 - qJDD(3) * t189, 0, qJD(3) * t263 + qJDD(3) * t305 + t204 * t180 + t203 * t189 - t230 * t274 + t237 * t294, -qJD(3) * t269 - qJDD(3) * t345 + t230 * t151 - t204 * t179 + t203 * t190 - t236 * t294, t112 * t280 + t153 * t74, -t112 * t144 - t113 * t280 - t153 * t75 + t279 * t74, t112 * t246 + t153 * t242, -t113 * t246 + t242 * t279, 0, t156 * t113 - t135 * t279 + t144 * t392 + t164 * t75 + t229 * t294 + t242 * t281 - t40 * t246, -t100 * t242 + t156 * t112 + t135 * t153 + t164 * t74 - t228 * t294 - t39 * t246 + t280 * t392, t134 * t277 + t358 * t49 (-t132 * t256 - t366) * t112 + (-t47 - t256 * t50 + (t132 * t251 - t134 * t256) * qJD(5)) * t153, t134 * t113 - t277 * t417 - t279 * t49 + t358 * t72, -t132 * t113 + t278 * t417 + t279 * t50 - t359 * t72, -t113 * t417 - t279 * t72 -(-t100 * t340 + t52) * t417 + t98 * t72 - t298 * t279 + t35 * t113 + t40 * t132 - t281 * t50 + t76 * t319 - g(1) * t172 - g(2) * t174 + (-(-qJD(5) * t101 - t39) * t417 - t100 * t72 - t310 * t279 + t19 * t153 + t76 * t112) * t251, t276 * t417 - t381 * t72 + t282 * t279 - t36 * t113 + t40 * t134 - t281 * t49 + t76 * t369 - g(1) * t171 - g(2) * t173 + (t19 * t256 - t66) * t153, -t111 * t14 - t22 * t287, -t110 * t14 - t111 * t267 - t22 * t88 + t23 * t287, -t111 * t71 - t113 * t287 + t139 * t22 - t14 * t279, -t110 * t71 - t113 * t88 - t139 * t23 - t267 * t279, t113 * t139 - t279 * t71 (-t250 * t5 + t255 * t4) * t139 + (-t250 * t38 + t255 * t32) * t71 - t323 * t279 + t6 * t113 + t20 * t88 - t62 * t267 + t10 * t110 + t54 * t23 - g(1) * t161 - g(2) * t163 + ((-t250 * t32 - t255 * t38) * t139 + t7 * t279) * qJD(6), -g(1) * t160 - g(2) * t162 - t10 * t111 - t7 * t113 + t62 * t14 - t25 * t279 - t20 * t287 + t54 * t22 + (-(-qJD(6) * t38 + t4) * t139 - t32 * t71 + t2 * t279) * t250 + (-(qJD(6) * t32 + t5) * t139 - t38 * t71 + t314 * t279) * t255; 0, 0, 0, -t332, t333, -t308, -qJ(2) * t308 + qJDD(2) + t272, 0, 0, 0, 0, 0, -t274 + t343 (-t177 - t321) * qJD(3) + t324, 0, 0, 0, 0, 0, t75 + t365, t74 - t362, 0, 0, 0, 0, 0, t284 - t368, -t256 * t417 ^ 2 - t367 - t64, 0, 0, 0, 0, 0, t297 - t380, t378 - t450; 0, 0, 0, 0, 0, 0, 0, t178 * t177, -t177 ^ 2 + t178 ^ 2 (t177 - t321) * qJD(3) + t324, t274 + t343, qJDD(3), -g(3) * t237 - t204 * t178 + t236 * t295 + t306, g(3) * t236 + t204 * t177 + t237 * t295 - t286, t422, t416, t427, t420, t242, t85 * t246 + (-t144 * t178 + t242 * t393 - t246 * t342) * pkin(3) + t396, t86 * t246 + (-t178 * t280 - t242 * t252 - t246 * t317) * pkin(3) + t410, t436, t366 * t417 + t437, t433, t284 + t368, t405, t233 * t50 - t429 * t256 + t289 * t251 + t300 * t132 - (-t232 * t340 + t440) * t417 + t398, t233 * t49 + t289 * t256 - t251 * t445 + t300 * t134 - (t232 * t341 + t403) * t417 + t397, t434, t432, t435, t443, -t406 (t184 * t255 - t185 * t250) * t71 - t212 * t267 + t373 * t88 + (t250 * t293 - t255 * t292) * t139 + t439 -(t184 * t250 + t185 * t255) * t71 + t212 * t14 - t373 * t287 + (t250 * t292 + t255 * t293) * t139 + t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t422, t416, t427, t420, t242, t246 * t83 + t396, t82 * t246 + t410, t436, t134 * t304 + t437, t433, -t304 * t417 + t368 + t65, t405, -pkin(4) * t50 + t107 * t417 - t83 * t132 + (-t417 * t82 + t296) * t251 + (-t429 + t327) * t256 + t398, -pkin(4) * t49 - t372 * t417 - t83 * t134 + t296 * t256 + (-t445 - t327) * t251 + t397, t434, t432, t435, t443, -t406 (t213 * t255 - t214 * t250) * t71 - t234 * t267 + t299 * t88 + (t250 * t291 - t255 * t290) * t139 + t439 -(t213 * t250 + t214 * t255) * t71 + t234 * t14 - t299 * t287 + (t250 * t290 + t255 * t291) * t139 + t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134 * t132, -t132 ^ 2 + t134 ^ 2, -t132 * t417 + t49, -t309 + (-qJD(5) - t417) * t134, t72, -g(1) * t173 + g(2) * t171 - t76 * t134 - t36 * t417 + (t310 + t223) * t251 + t298, g(1) * t174 - g(2) * t172 + t132 * t76 + t223 * t256 - t35 * t417 - t282, -t430, t415, t414, t411, t71 -(-t250 * t29 - t379) * t139 + (-t134 * t88 - t139 * t339 + t255 * t71) * pkin(5) + t412 (-t139 * t30 - t2) * t250 + (t139 * t29 - t314) * t255 + (t134 * t287 - t139 * t338 - t250 * t71) * pkin(5) + t413; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t430, t415, t414, t411, t71, t7 * t139 + t412, t6 * t139 - t250 * t2 - t255 * t314 + t413;];
tau_reg  = t1;
