% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:30
% EndTime: 2019-03-09 13:53:43
% DurationCPUTime: 5.89s
% Computational Cost: add. (6880->472), mult. (14918->587), div. (0->0), fcn. (11056->12), ass. (0->259)
t231 = cos(qJ(6));
t315 = qJD(6) * t231;
t228 = sin(qJ(4));
t233 = cos(qJ(4));
t234 = cos(qJ(2));
t323 = qJD(1) * t234;
t229 = sin(qJ(2));
t324 = qJD(1) * t229;
t133 = t228 * t324 + t233 * t323;
t136 = -t228 * t323 + t233 * t324;
t227 = sin(qJ(5));
t232 = cos(qJ(5));
t377 = t232 * t133 + t227 * t136;
t395 = t231 * t377;
t417 = t315 + t395;
t205 = pkin(7) * t324;
t416 = -pkin(8) * t324 + qJD(3) + t205;
t217 = qJDD(2) - qJDD(4);
t209 = -qJDD(5) + t217;
t218 = qJD(2) - qJD(4);
t210 = -qJD(5) + t218;
t226 = sin(qJ(6));
t317 = qJD(5) * t232;
t318 = qJD(5) * t227;
t150 = t229 * t228 + t234 * t233;
t257 = t150 * qJD(4);
t314 = qJD(1) * qJD(2);
t305 = t229 * t314;
t313 = t234 * qJDD(1);
t392 = t305 - t313;
t211 = t229 * qJDD(1);
t304 = t234 * t314;
t393 = -t304 - t211;
t65 = -qJD(1) * t257 + t392 * t228 - t393 * t233;
t319 = qJD(4) * t233;
t320 = qJD(4) * t228;
t321 = qJD(2) * t234;
t409 = t228 * t321 + t229 * t319 - t234 * t320;
t66 = qJD(1) * t409 + qJDD(1) * t150 - t233 * t305;
t261 = t133 * t317 + t136 * t318 + t227 * t66 - t232 * t65;
t316 = qJD(6) * t226;
t77 = t227 * t133 - t232 * t136;
t14 = -t226 * t209 - t210 * t315 - t231 * t261 + t316 * t77;
t271 = t226 * t210 + t231 * t77;
t15 = -qJD(6) * t271 + t231 * t209 - t226 * t261;
t390 = qJD(6) + t377;
t415 = t226 * t390;
t67 = t231 * t210 - t226 * t77;
t414 = t14 * t231 - t226 * t15 + t271 * t415 - t417 * t67;
t245 = qJD(5) * t77 - t227 * t65 - t232 * t66;
t25 = qJDD(6) - t245;
t23 = t231 * t25;
t413 = t390 * t415 + t67 * t77 - t23;
t206 = pkin(7) * t323;
t160 = -pkin(8) * t323 + t206;
t236 = -pkin(2) - pkin(3);
t290 = -t228 * qJ(3) + t233 * t236;
t411 = t290 * qJD(4) - t228 * t160 + t233 * t416;
t162 = t233 * qJ(3) + t228 * t236;
t410 = qJD(4) * t162 + t233 * t160 + t228 * t416;
t384 = t77 * t210 + t245;
t309 = t236 * qJD(2);
t117 = t309 + t416;
t221 = qJD(2) * qJ(3);
t137 = t160 + t221;
t270 = t228 * t117 + t233 * t137;
t189 = pkin(7) * t304;
t200 = pkin(7) * t211;
t303 = qJDD(3) + t189 + t200;
t90 = t393 * pkin(8) + t236 * qJDD(2) + t303;
t201 = pkin(7) * t313;
t219 = qJDD(2) * qJ(3);
t220 = qJD(2) * qJD(3);
t115 = -pkin(7) * t305 + t201 + t219 + t220;
t92 = t392 * pkin(8) + t115;
t408 = -t270 * qJD(4) - t228 * t92 + t233 * t90;
t360 = t133 * pkin(9);
t58 = t270 - t360;
t348 = t232 * t58;
t288 = t233 * t117 - t228 * t137;
t358 = t136 * pkin(9);
t57 = t288 - t358;
t54 = -t218 * pkin(4) + t57;
t33 = t227 * t54 + t348;
t31 = -t210 * pkin(10) + t33;
t138 = -qJD(1) * pkin(1) - pkin(2) * t323 - qJ(3) * t324;
t113 = pkin(3) * t323 - t138;
t84 = t133 * pkin(4) + t113;
t36 = pkin(5) * t377 + pkin(10) * t77 + t84;
t10 = t226 * t36 + t231 * t31;
t18 = -t217 * pkin(4) - t65 * pkin(9) + t408;
t262 = -t117 * t319 + t137 * t320 - t228 * t90 - t233 * t92;
t21 = -t66 * pkin(9) - t262;
t298 = -t232 * t18 + t227 * t21 + t58 * t317 + t54 * t318;
t3 = t209 * pkin(5) + t298;
t352 = t227 * t58;
t32 = t232 * t54 - t352;
t30 = t210 * pkin(5) - t32;
t407 = t10 * t77 - t3 * t226 - t30 * t315;
t12 = t14 * t226;
t405 = t271 * t417 - t12;
t22 = t226 * t25;
t73 = t390 * t315;
t404 = -t271 * t77 + t390 * t395 + t22 + t73;
t402 = t30 * t377;
t401 = t390 * t77;
t399 = t77 * t377;
t398 = -t358 + t411;
t397 = t360 - t410;
t396 = -t133 * t218 + t65;
t230 = sin(qJ(1));
t267 = t234 * t228 - t229 * t233;
t124 = t267 * t230;
t235 = cos(qJ(1));
t335 = t234 * t235;
t338 = t229 * t235;
t126 = t228 * t335 - t233 * t338;
t394 = -g(1) * t126 - g(2) * t124 - g(3) * t150 + t113 * t136 - t408;
t334 = qJ(4) + qJ(5);
t213 = sin(t334);
t306 = cos(t334);
t275 = t229 * t306;
t376 = -t234 * t213 + t275;
t109 = t376 * t230;
t111 = t213 * t335 - t235 * t275;
t251 = t229 * t213 + t234 * t306;
t256 = g(1) * t111 - g(2) * t109 + g(3) * t251;
t389 = t377 ^ 2 - t77 ^ 2;
t9 = -t226 * t31 + t231 * t36;
t388 = t30 * t316 + t9 * t77;
t386 = -t256 * t226 - t407;
t385 = t84 * t77 + t256 - t298;
t112 = t251 * t235;
t110 = t251 * t230;
t277 = g(2) * t110 + g(3) * t376;
t297 = -t227 * t18 - t232 * t21 - t54 * t317 + t58 * t318;
t383 = g(1) * t112 + t84 * t377 + t277 + t297;
t381 = t210 * t377 + t261;
t49 = -t77 * pkin(5) + pkin(10) * t377;
t326 = t234 * pkin(2) + t229 * qJ(3);
t166 = -pkin(1) - t326;
t198 = t227 * pkin(4) + pkin(10);
t359 = t136 * pkin(4);
t380 = (qJD(6) * t198 + t359 + t49) * t390;
t196 = qJ(3) * t323;
t123 = t236 * t324 + t196;
t89 = t123 - t359;
t157 = -pkin(4) + t290;
t330 = t227 * t157 + t232 * t162;
t98 = -pkin(10) + t330;
t379 = (qJD(6) * t98 - t49 + t89) * t390;
t378 = t136 * t218 + t66;
t370 = pkin(7) - pkin(8);
t169 = t370 * t229;
t170 = t370 * t234;
t329 = t228 * t169 + t233 * t170;
t362 = g(2) * t235;
t366 = g(1) * t230;
t375 = -t362 + t366;
t363 = g(2) * t230;
t365 = g(1) * t235;
t279 = t363 + t365;
t346 = pkin(7) * qJDD(2);
t372 = (qJD(1) * t166 + t138) * qJD(2) - t346;
t302 = -t209 * pkin(10) + qJD(6) * t36 - t297;
t322 = qJD(2) * t229;
t101 = -t233 * t322 + t409;
t102 = qJD(2) * t150 - t257;
t93 = t232 * t150 - t227 * t267;
t40 = -qJD(5) * t93 - t227 * t101 + t232 * t102;
t146 = t234 * pkin(3) - t166;
t103 = t150 * pkin(4) + t146;
t94 = -t227 * t150 - t232 * t267;
t44 = t93 * pkin(5) - t94 * pkin(10) + t103;
t286 = t233 * t169 - t228 * t170;
t74 = pkin(9) * t267 + t286;
t75 = -t150 * pkin(9) + t329;
t48 = t227 * t74 + t232 * t75;
t159 = t370 * t322;
t161 = qJD(2) * t170;
t259 = -t233 * t159 + t228 * t161 + t169 * t319 - t170 * t320;
t45 = -t101 * pkin(9) + t259;
t243 = -t329 * qJD(4) + t228 * t159 + t233 * t161;
t46 = -t102 * pkin(9) + t243;
t47 = t227 * t75 - t232 * t74;
t6 = -qJD(5) * t47 + t227 * t46 + t232 * t45;
t371 = -t48 * t25 + t3 * t94 + t30 * t40 - (qJD(6) * t44 + t6) * t390 - t302 * t93 + t365;
t368 = g(1) * t110;
t357 = t30 * t94;
t269 = t232 * t157 - t227 * t162;
t354 = qJD(5) * t269 + t397 * t227 + t398 * t232;
t353 = t330 * qJD(5) + t398 * t227 - t397 * t232;
t149 = t227 * t228 - t232 * t233;
t347 = t210 * t149;
t345 = qJD(6) * t31;
t224 = qJDD(1) * pkin(1);
t344 = qJDD(2) * pkin(2);
t342 = t136 * t133;
t238 = qJD(1) ^ 2;
t337 = t229 * t238;
t152 = t227 * t233 + t232 * t228;
t332 = t210 * t152;
t212 = t229 * qJD(3);
t327 = qJ(3) * t321 + t212;
t222 = t229 ^ 2;
t223 = t234 ^ 2;
t325 = t222 - t223;
t311 = t94 * t316;
t310 = t234 * t337;
t307 = t133 ^ 2 - t136 ^ 2;
t291 = -qJD(2) * pkin(2) + qJD(3);
t284 = t218 ^ 2;
t283 = t229 * t309;
t34 = t227 * t57 + t348;
t282 = pkin(4) * t318 - t34;
t237 = qJD(2) ^ 2;
t280 = pkin(7) * t237 + t362;
t278 = t44 * t25 + t368;
t276 = t25 * t94 + t390 * t40;
t274 = pkin(2) * t229 - qJ(3) * t234;
t273 = -t345 - t362;
t272 = pkin(2) * t313 - t393 * qJ(3) + qJD(1) * t212 + t224;
t163 = t205 + t291;
t167 = t206 + t221;
t268 = t163 * t234 - t167 * t229;
t266 = g(1) * t338 - g(3) * t234 + t229 * t363 - t200;
t265 = qJD(6) * t152 + t324;
t263 = -0.2e1 * pkin(1) * t314 - t346;
t260 = -qJDD(3) + t266;
t107 = t283 + t327;
t254 = -t280 + 0.2e1 * t224;
t252 = t256 - t3;
t250 = t277 - t302;
t249 = -t98 * t25 - t354 * t390 - t402;
t61 = t101 * pkin(4) + t107;
t35 = t232 * t57 - t352;
t244 = -t198 * t25 + t402 + (-pkin(4) * t317 + t35) * t390;
t128 = pkin(2) * t322 - t327;
t88 = pkin(2) * t305 - t272;
t242 = -qJD(1) * t128 - qJDD(1) * t166 - t280 - t88;
t70 = pkin(3) * t313 + qJD(1) * t283 + t272;
t42 = t66 * pkin(4) + t70;
t125 = t150 * t230;
t127 = t150 * t235;
t240 = g(1) * t127 + g(2) * t125 - g(3) * t267 + t113 * t133 + t262;
t122 = t303 - t344;
t239 = qJD(2) * t268 + t115 * t234 + t122 * t229 - t279;
t199 = -t232 * pkin(4) - pkin(5);
t191 = t234 * t366;
t156 = pkin(2) * t324 - t196;
t97 = pkin(5) - t269;
t96 = t112 * t231 - t230 * t226;
t95 = -t112 * t226 - t230 * t231;
t41 = qJD(5) * t94 + t232 * t101 + t227 * t102;
t8 = t41 * pkin(5) - t40 * pkin(10) + t61;
t7 = qJD(5) * t48 + t227 * t45 - t232 * t46;
t5 = -pkin(5) * t245 + pkin(10) * t261 + t42;
t4 = t231 * t5;
t1 = [qJDD(1), t375, t279, t222 * qJDD(1) + 0.2e1 * t229 * t304, 0.2e1 * t229 * t313 - 0.2e1 * t314 * t325, qJDD(2) * t229 + t237 * t234, qJDD(2) * t234 - t237 * t229, 0, t229 * t263 + t234 * t254 + t191, t263 * t234 + (-t254 - t366) * t229, t372 * t229 + t242 * t234 + t191 (t222 + t223) * qJDD(1) * pkin(7) + t239, -t372 * t234 + (t242 + t366) * t229, pkin(7) * t239 + t138 * t128 + (-t375 + t88) * t166, t136 * t102 - t267 * t65, -t136 * t101 - t102 * t133 - t65 * t150 + t267 * t66, -t102 * t218 + t217 * t267, t101 * t218 + t150 * t217, 0, g(1) * t125 - g(2) * t127 + t113 * t101 + t107 * t133 + t146 * t66 + t70 * t150 - t217 * t286 - t218 * t243, -g(1) * t124 + g(2) * t126 + t113 * t102 + t107 * t136 + t146 * t65 + t217 * t329 + t218 * t259 - t267 * t70, -t261 * t94 - t40 * t77, t245 * t94 + t261 * t93 - t377 * t40 + t41 * t77, -t94 * t209 - t40 * t210, t93 * t209 + t41 * t210, 0, -g(2) * t112 - t103 * t245 + t47 * t209 + t7 * t210 + t377 * t61 + t84 * t41 + t42 * t93 + t368, g(1) * t109 + g(2) * t111 - t103 * t261 + t48 * t209 + t6 * t210 + t84 * t40 + t42 * t94 - t61 * t77, t271 * t311 + (t14 * t94 - t271 * t40) * t231 (t226 * t271 - t231 * t67) * t40 + (-t12 - t15 * t231 + (t226 * t67 + t231 * t271) * qJD(6)) * t94, t14 * t93 + t231 * t276 - t271 * t41 - t311 * t390, -t15 * t93 - t226 * t276 - t67 * t41 - t73 * t94, t25 * t93 + t390 * t41, -g(2) * t96 + t47 * t15 + t4 * t93 + t9 * t41 + t7 * t67 + (t8 * t390 + (-t31 * t93 - t390 * t48 + t357) * qJD(6) + t278) * t231 + t371 * t226, -g(2) * t95 - t10 * t41 + t47 * t14 - t7 * t271 + (-(-qJD(6) * t48 + t8) * t390 - (t5 - t345) * t93 - qJD(6) * t357 - t278) * t226 + t371 * t231; 0, 0, 0, -t310, t325 * t238, t211, t313, qJDD(2), pkin(1) * t337 + t266, g(3) * t229 - t201 + (pkin(1) * t238 + t279) * t234, 0.2e1 * t344 + (-t138 * t229 + t156 * t234) * qJD(1) + t260, -t274 * qJDD(1) + ((t167 - t221) * t229 + (-t163 + t291) * t234) * qJD(1), t201 + 0.2e1 * t219 + 0.2e1 * t220 + (qJD(1) * t156 - g(3)) * t229 + (qJD(1) * t138 - t279) * t234, -t268 * qJD(1) * pkin(7) - t122 * pkin(2) - g(3) * t326 + t115 * qJ(3) + t167 * qJD(3) - t138 * t156 + t274 * t279, -t342, t307, -t396, t378, t217, -t123 * t133 - t290 * t217 + t218 * t410 + t394, -t123 * t136 + t162 * t217 + t218 * t411 - t240, t399, t389, t381, -t384, t209, -t209 * t269 + t210 * t353 - t377 * t89 - t385, t209 * t330 + t210 * t354 + t77 * t89 - t383, t405, -t414, -t404, t413, -t401, t97 * t15 + t353 * t67 + t249 * t226 + (-t252 - t379) * t231 - t388, t97 * t14 - t353 * t271 + t249 * t231 + (t256 + t379) * t226 + t407; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t310, t211, -t222 * t238 - t237, -t167 * qJD(2) + t138 * t324 + t189 - t260 - t344, 0, 0, 0, 0, 0, -t133 * t324 - t233 * t217 - t228 * t284, -t136 * t324 + t228 * t217 - t233 * t284, 0, 0, 0, 0, 0, t149 * t209 - t210 * t332 - t324 * t377, t152 * t209 + t210 * t347 + t324 * t77, 0, 0, 0, 0, 0, -t152 * t22 + t149 * t15 - t332 * t67 + (-t226 * t347 - t231 * t265) * t390, -t152 * t23 + t149 * t14 + t332 * t271 + (t226 * t265 - t231 * t347) * t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, -t307, t396, -t378, -t217, -t218 * t270 - t394, -t218 * t288 + t240, -t399, -t389, -t381, t384, -t209, -t34 * t210 + (-t136 * t377 - t209 * t232 + t210 * t318) * pkin(4) + t385, -t35 * t210 + (t136 * t77 + t209 * t227 + t210 * t317) * pkin(4) + t383, -t405, t414, t404, -t413, t401, t199 * t15 + t282 * t67 + t244 * t226 + (t252 - t380) * t231 + t388, t199 * t14 + t226 * t380 + t231 * t244 - t271 * t282 + t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t399, -t389, -t381, t384, -t209, -t210 * t33 + t385, -t210 * t32 + t383, -t405, t414, t404, -t413, t401, -pkin(5) * t15 - t33 * t67 + (-pkin(10) * t25 + t32 * t390 + t402) * t226 + ((-pkin(10) * qJD(6) - t49) * t390 + t252) * t231 + t388, -pkin(5) * t14 + (t226 * t49 + t231 * t32) * t390 + t33 * t271 + t30 * t395 + (t316 * t390 - t23) * pkin(10) + t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271 * t67, t271 ^ 2 - t67 ^ 2, t390 * t67 + t14, -t271 * t390 - t15, t25, -g(1) * t95 + t10 * t390 + t226 * t250 + t231 * t273 + t271 * t30 + t4, g(1) * t96 + t30 * t67 + t9 * t390 + (-t273 - t5) * t226 + t250 * t231;];
tau_reg  = t1;
