% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR2
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
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:34
% EndTime: 2019-03-09 06:58:47
% DurationCPUTime: 5.22s
% Computational Cost: add. (6720->433), mult. (14084->585), div. (0->0), fcn. (10544->18), ass. (0->258)
t220 = qJD(3) + qJD(4);
t230 = sin(qJ(4));
t235 = cos(qJ(3));
t357 = cos(qJ(4));
t297 = qJD(1) * t357;
t231 = sin(qJ(3));
t319 = qJD(1) * t231;
t374 = -t230 * t319 + t235 * t297;
t377 = t220 * t374;
t225 = qJ(3) + qJ(4);
t215 = sin(t225);
t221 = qJ(1) + pkin(11);
t210 = sin(t221);
t211 = cos(t221);
t368 = g(1) * t211 + g(2) * t210;
t376 = t368 * t215;
t360 = qJD(5) + qJD(6);
t375 = -t374 + t360;
t296 = qJD(4) * t357;
t226 = sin(pkin(11));
t202 = pkin(1) * t226 + pkin(7);
t352 = pkin(8) + t202;
t289 = t352 * qJD(1);
t138 = t231 * qJD(2) + t289 * t235;
t125 = t230 * t138;
t137 = t235 * qJD(2) - t289 * t231;
t84 = t357 * t137 - t125;
t370 = pkin(3) * t296 - t84;
t324 = t230 * t235;
t156 = -qJD(1) * t324 - t231 * t297;
t113 = -pkin(4) * t156 - pkin(9) * t374;
t104 = pkin(3) * t319 + t113;
t229 = sin(qJ(5));
t234 = cos(qJ(5));
t373 = -t234 * t104 - t229 * t370;
t128 = -t156 * t229 - t234 * t220;
t228 = sin(qJ(6));
t233 = cos(qJ(6));
t261 = t156 * t234 - t220 * t229;
t262 = t128 * t228 + t233 * t261;
t72 = t233 * t128 - t228 * t261;
t372 = t262 * t72;
t325 = t228 * t234;
t166 = t229 * t233 + t325;
t371 = t375 * t166;
t164 = t228 * t229 - t233 * t234;
t322 = t375 * t164;
t367 = t262 ^ 2 - t72 ^ 2;
t151 = qJD(5) - t374;
t147 = qJD(6) + t151;
t314 = qJD(6) * t233;
t315 = qJD(6) * t228;
t219 = qJDD(3) + qJDD(4);
t316 = qJD(5) * t234;
t317 = qJD(5) * t229;
t290 = qJDD(1) * t357;
t310 = t235 * qJDD(1);
t93 = t230 * t310 + t231 * t290 + t377;
t54 = t156 * t317 + t229 * t219 + t220 * t316 + t234 * t93;
t55 = -qJD(5) * t261 - t234 * t219 + t229 * t93;
t17 = -t128 * t314 - t228 * t55 + t233 * t54 + t261 * t315;
t366 = t147 * t72 + t17;
t224 = qJ(5) + qJ(6);
t214 = sin(t224);
t216 = cos(t224);
t217 = cos(t225);
t328 = t216 * t217;
t132 = -t210 * t328 + t211 * t214;
t134 = t210 * t214 + t211 * t328;
t204 = g(3) * t215;
t126 = t357 * t138;
t346 = qJD(3) * pkin(3);
t127 = t137 + t346;
t79 = t230 * t127 + t126;
t69 = pkin(9) * t220 + t79;
t227 = cos(pkin(11));
t203 = -pkin(1) * t227 - pkin(2);
t178 = -pkin(3) * t235 + t203;
t157 = t178 * qJD(1);
t95 = -pkin(4) * t374 + t156 * pkin(9) + t157;
t44 = t229 * t95 + t234 * t69;
t32 = -pkin(10) * t128 + t44;
t30 = t32 * t315;
t78 = t357 * t127 - t125;
t68 = -t220 * pkin(4) - t78;
t52 = t128 * pkin(5) + t68;
t365 = g(1) * t134 - g(2) * t132 + t216 * t204 + t52 * t72 + t30;
t329 = t214 * t217;
t131 = t210 * t329 + t211 * t216;
t133 = t210 * t216 - t211 * t329;
t318 = qJD(4) * t230;
t212 = t235 * qJDD(2);
t182 = t202 * qJDD(1);
t286 = pkin(8) * qJDD(1) + t182;
t82 = qJDD(3) * pkin(3) - qJD(3) * t138 - t286 * t231 + t212;
t91 = qJD(3) * t137 + t231 * qJDD(2) + t286 * t235;
t243 = t127 * t296 - t138 * t318 + t230 * t82 + t357 * t91;
t21 = t219 * pkin(9) + t243;
t312 = qJD(1) * qJD(3);
t295 = t231 * t312;
t139 = pkin(3) * t295 + qJDD(1) * t178;
t167 = t357 * t231 + t324;
t120 = t220 * t167;
t311 = t231 * qJDD(1);
t270 = t230 * t311 - t235 * t290;
t94 = qJD(1) * t120 + t270;
t38 = pkin(4) * t94 - pkin(9) * t93 + t139;
t37 = t234 * t38;
t90 = qJDD(5) + t94;
t3 = pkin(5) * t90 - pkin(10) * t54 - t44 * qJD(5) - t229 * t21 + t37;
t257 = t234 * t21 + t229 * t38 + t95 * t316 - t69 * t317;
t4 = -pkin(10) * t55 + t257;
t303 = -t228 * t4 + t233 * t3;
t43 = -t229 * t69 + t234 * t95;
t31 = pkin(10) * t261 + t43;
t25 = pkin(5) * t151 + t31;
t344 = t233 * t32;
t9 = t228 * t25 + t344;
t364 = -g(1) * t133 + g(2) * t131 - t9 * qJD(6) + t214 * t204 + t52 * t262 + t303;
t247 = qJD(6) * t262 - t228 * t54 - t233 * t55;
t363 = -t147 * t262 + t247;
t83 = t230 * t137 + t126;
t282 = pkin(3) * t318 - t83;
t109 = t166 * t167;
t333 = t374 * t229;
t362 = (t317 - t333) * pkin(5);
t361 = t229 * t104 - t234 * t370;
t300 = t167 * t317;
t255 = -t230 * t231 + t357 * t235;
t119 = t220 * t255;
t338 = t119 * t234;
t253 = t300 - t338;
t330 = t167 * t234;
t359 = -t151 * t253 + t90 * t330;
t358 = -pkin(9) - pkin(10);
t354 = g(3) * t217;
t353 = t234 * pkin(5);
t218 = t234 * pkin(10);
t207 = pkin(3) * t230 + pkin(9);
t351 = -pkin(10) - t207;
t350 = -t120 * t262 - t17 * t255;
t331 = t167 * t229;
t339 = t119 * t229;
t29 = t119 * t325 - t228 * t300 - t315 * t331 + (t330 * t360 + t339) * t233;
t86 = qJDD(6) + t90;
t349 = -t109 * t86 - t29 * t147;
t348 = -t120 * t261 - t255 * t54;
t345 = t374 * t68;
t343 = t54 * t229;
t341 = t362 + t282;
t340 = t229 * t113 + t234 * t78;
t337 = t128 * t151;
t336 = t261 * t151;
t335 = t147 * t156;
t334 = t151 * t156;
t332 = t156 * t374;
t327 = t217 * t229;
t326 = t217 * t234;
t158 = t352 * t231;
t159 = t352 * t235;
t112 = -t230 * t158 + t357 * t159;
t105 = t234 * t112;
t323 = qJDD(2) - g(3);
t106 = -pkin(4) * t255 - pkin(9) * t167 + t178;
t321 = t229 * t106 + t105;
t222 = t231 ^ 2;
t320 = -t235 ^ 2 + t222;
t185 = qJD(1) * t203;
t309 = pkin(10) * t333;
t308 = t231 * t346;
t305 = qJD(5) * pkin(9) * t151;
t66 = t68 * t317;
t302 = qJD(5) * t358;
t299 = t167 * t316;
t284 = -t127 * t318 - t138 * t296 - t230 * t91 + t357 * t82;
t22 = -pkin(4) * t219 - t284;
t298 = -t22 - t354;
t294 = qJD(6) * t25 + t4;
t292 = qJD(3) * t352;
t291 = qJD(5) * t351;
t288 = -qJD(5) * t95 - t21;
t285 = t151 * t234;
t208 = -t357 * pkin(3) - pkin(4);
t281 = -t79 + t362;
t280 = -t156 * pkin(5) - t218 * t374;
t279 = -t69 * t316 + t37;
t278 = -pkin(9) * t90 - t345;
t276 = g(1) * t210 - g(2) * t211;
t232 = sin(qJ(1));
t236 = cos(qJ(1));
t275 = g(1) * t232 - g(2) * t236;
t160 = t351 * t229;
t274 = -qJD(6) * t160 - t229 * t291 - t309 + t361;
t161 = t207 * t234 + t218;
t273 = qJD(6) * t161 - t234 * t291 + t280 - t373;
t186 = t358 * t229;
t272 = -qJD(6) * t186 - t229 * t302 - t309 + t340;
t108 = t234 * t113;
t187 = pkin(9) * t234 + t218;
t271 = qJD(6) * t187 - t229 * t78 - t234 * t302 + t108 + t280;
t110 = t164 * t167;
t28 = -t109 * t360 - t164 * t119;
t269 = t110 * t86 - t147 * t28;
t268 = -t120 * t72 - t247 * t255;
t267 = -t207 * t90 - t345;
t264 = -t120 * t128 + t255 * t55;
t263 = t119 * t220 + t167 * t219;
t260 = g(3) * t327 - t44 * t156 + t22 * t229 + t68 * t316;
t259 = t43 * t156 + t234 * t376 + t66;
t256 = -t357 * t158 - t230 * t159;
t254 = t299 + t339;
t148 = t231 * t292;
t149 = t235 * t292;
t57 = t256 * qJD(4) - t357 * t148 - t230 * t149;
t65 = pkin(4) * t120 - pkin(9) * t119 + t308;
t251 = t106 * t316 - t112 * t317 + t229 * t65 + t234 * t57;
t250 = -qJD(1) * t185 - t182 + t368;
t249 = 0.2e1 * t185 * qJD(3) - qJDD(3) * t202;
t248 = t157 * t156 + t284 - t354 + t376;
t237 = qJD(3) ^ 2;
t246 = -0.2e1 * qJDD(1) * t203 - t202 * t237 + t276;
t13 = pkin(5) * t55 + t22;
t8 = -t228 * t32 + t233 * t25;
t245 = -g(3) * t328 + t13 * t164 + t8 * t156 + t216 * t376 + t371 * t52;
t244 = -t151 * t254 - t90 * t331;
t58 = t112 * qJD(4) - t230 * t148 + t357 * t149;
t242 = g(3) * t329 + t13 * t166 - t9 * t156 - t214 * t376 - t322 * t52;
t240 = -t157 * t374 + t217 * t368 + t204 - t243;
t238 = qJD(1) ^ 2;
t209 = -pkin(4) - t353;
t181 = t208 - t353;
t180 = qJDD(3) * t235 - t231 * t237;
t179 = qJDD(3) * t231 + t235 * t237;
t143 = t210 * t229 + t211 * t326;
t142 = t210 * t234 - t211 * t327;
t141 = -t210 * t326 + t211 * t229;
t140 = t210 * t327 + t211 * t234;
t103 = t234 * t106;
t96 = t156 ^ 2 - t374 ^ 2;
t92 = -t120 * t220 + t219 * t255;
t88 = pkin(5) * t331 - t256;
t64 = -t270 + (-qJD(1) * t167 - t156) * t220;
t63 = t93 - t377;
t62 = t234 * t65;
t47 = -pkin(10) * t331 + t321;
t45 = -pkin(5) * t255 - pkin(10) * t330 - t112 * t229 + t103;
t35 = pkin(5) * t254 + t58;
t27 = t151 * t285 - t156 * t261 + t229 * t90;
t26 = -t151 ^ 2 * t229 - t128 * t156 + t234 * t90;
t24 = -t261 * t285 + t343;
t16 = -t147 * t371 - t72 * t156 - t164 * t86;
t15 = -t322 * t147 - t156 * t262 + t166 * t86;
t10 = -pkin(10) * t254 + t251;
t7 = -pkin(10) * t338 + pkin(5) * t120 - t229 * t57 + t62 + (-t105 + (pkin(10) * t167 - t106) * t229) * qJD(5);
t6 = (t54 - t337) * t234 + (-t55 + t336) * t229;
t5 = t17 * t166 + t262 * t322;
t1 = -t164 * t17 + t166 * t247 + t262 * t371 + t322 * t72;
t2 = [qJDD(1), t275, g(1) * t236 + g(2) * t232 (t275 + (t226 ^ 2 + t227 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t222 + 0.2e1 * t235 * t295, 0.2e1 * t231 * t310 - 0.2e1 * t320 * t312, t179, t180, 0, t231 * t249 + t235 * t246, -t231 * t246 + t235 * t249, -t119 * t156 + t167 * t93, t119 * t374 + t120 * t156 - t167 * t94 + t255 * t93, t263, t92, 0, t120 * t157 - t139 * t255 + t178 * t94 + t217 * t276 + t219 * t256 - t220 * t58 - t308 * t374, -t112 * t219 + t119 * t157 + t139 * t167 - t156 * t308 + t178 * t93 - t215 * t276 - t220 * t57, t253 * t261 + t54 * t330 (-t128 * t234 + t229 * t261) * t119 + (-t343 - t234 * t55 + (t128 * t229 + t234 * t261) * qJD(5)) * t167, t348 + t359, t244 + t264, t120 * t151 - t255 * t90 (-t112 * t316 + t62) * t151 + t103 * t90 - t279 * t255 + t43 * t120 + t58 * t128 - t256 * t55 + t68 * t299 - g(1) * t141 - g(2) * t143 + ((-qJD(5) * t106 - t57) * t151 - t112 * t90 - t288 * t255 + t22 * t167 + t68 * t119) * t229, -t251 * t151 - t321 * t90 + t257 * t255 - t44 * t120 - t58 * t261 - t256 * t54 + t68 * t338 - g(1) * t140 - g(2) * t142 + (t22 * t234 - t66) * t167, -t110 * t17 - t262 * t28, -t109 * t17 - t110 * t247 + t262 * t29 - t28 * t72, -t269 + t350, t268 + t349, t120 * t147 - t255 * t86 (-t228 * t10 + t233 * t7) * t147 + (-t228 * t47 + t233 * t45) * t86 - t303 * t255 + t8 * t120 + t35 * t72 - t88 * t247 + t13 * t109 + t52 * t29 - g(1) * t132 - g(2) * t134 + ((-t228 * t45 - t233 * t47) * t147 + t9 * t255) * qJD(6), -g(1) * t131 - g(2) * t133 - t13 * t110 - t9 * t120 - t30 * t255 + t88 * t17 + t52 * t28 - t35 * t262 + (-(-qJD(6) * t47 + t7) * t147 - t45 * t86 + t3 * t255) * t228 + (-(qJD(6) * t45 + t10) * t147 - t47 * t86 + t294 * t255) * t233; 0, 0, 0, t323, 0, 0, 0, 0, 0, t180, -t179, 0, 0, 0, 0, 0, t92, -t263, 0, 0, 0, 0, 0, t244 - t264, t348 - t359, 0, 0, 0, 0, 0, -t268 + t349, t269 + t350; 0, 0, 0, 0, -t231 * t238 * t235, t320 * t238, t311, t310, qJDD(3), -g(3) * t235 + t231 * t250 + t212, -t323 * t231 + t250 * t235, t332, t96, t63, t64, t219, t83 * t220 + (t357 * t219 - t220 * t318 + t319 * t374) * pkin(3) + t248, t84 * t220 + (t156 * t319 - t219 * t230 - t220 * t296) * pkin(3) + t240, t24, t6, t27, t26, t334, t208 * t55 + t298 * t234 + t267 * t229 + t282 * t128 + (-t207 * t316 + t373) * t151 + t259, t208 * t54 + t267 * t234 - t229 * t376 - t282 * t261 + (t207 * t317 + t361) * t151 + t260, t5, t1, t15, t16, t335 (t160 * t233 - t161 * t228) * t86 - t181 * t247 + t341 * t72 + (t228 * t274 - t233 * t273) * t147 + t245 -(t160 * t228 + t161 * t233) * t86 + t181 * t17 - t341 * t262 + (t228 * t273 + t233 * t274) * t147 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, t96, t63, t64, t219, t220 * t79 + t248, t78 * t220 + t240, t24, t6, t27, t26, t334, -pkin(4) * t55 - t108 * t151 - t79 * t128 + (t151 * t78 + t278) * t229 + (t298 - t305) * t234 + t259, -pkin(4) * t54 + t340 * t151 + t79 * t261 + t278 * t234 + (-t376 + t305) * t229 + t260, t5, t1, t15, t16, t335 (t186 * t233 - t187 * t228) * t86 - t209 * t247 + t281 * t72 + (t228 * t272 - t233 * t271) * t147 + t245 -(t186 * t228 + t187 * t233) * t86 + t209 * t17 - t281 * t262 + (t228 * t271 + t233 * t272) * t147 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261 * t128, -t128 ^ 2 + t261 ^ 2, t54 + t337, -t55 - t336, t90, -g(1) * t142 + g(2) * t140 + t261 * t68 + t151 * t44 + (t288 + t204) * t229 + t279, g(1) * t143 - g(2) * t141 + t128 * t68 + t151 * t43 + t234 * t204 - t257, -t372, t367, t366, t363, t86 -(-t228 * t31 - t344) * t147 + (-t147 * t315 + t233 * t86 + t261 * t72) * pkin(5) + t364 (-t147 * t32 - t3) * t228 + (t147 * t31 - t294) * t233 + (-t147 * t314 - t228 * t86 - t261 * t262) * pkin(5) + t365; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t372, t367, t366, t363, t86, t147 * t9 + t364, t147 * t8 - t228 * t3 - t233 * t294 + t365;];
tau_reg  = t2;
