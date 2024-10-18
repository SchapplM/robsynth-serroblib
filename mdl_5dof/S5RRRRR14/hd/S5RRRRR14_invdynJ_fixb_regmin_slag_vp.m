% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR14
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR14_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:38
% EndTime: 2024-09-27 18:43:42
% DurationCPUTime: 3.39s
% Computational Cost: add. (6254->381), mult. (11210->523), div. (0->0), fcn. (8397->26), ass. (0->273)
t230 = sin(pkin(5));
t238 = cos(qJ(4));
t239 = cos(qJ(3));
t339 = t238 * t239;
t310 = t230 * t339;
t233 = sin(qJ(4));
t234 = sin(qJ(3));
t344 = t233 * t234;
t398 = -t230 * t344 + t310;
t272 = t233 * t239 + t234 * t238;
t397 = t272 * t230;
t235 = sin(qJ(2));
t240 = cos(qJ(2));
t231 = cos(pkin(5));
t330 = qJD(3) * t239;
t305 = t231 * t330;
t346 = t231 * t234;
t362 = pkin(1) * qJD(1);
t396 = (-t235 * t346 + t239 * t240) * t362 - pkin(2) * t305;
t228 = qJ(3) + qJ(4);
t304 = pkin(5) - t228;
t285 = -qJ(5) + t304;
t267 = sin(t285);
t214 = pkin(5) + t228;
t208 = qJ(5) + t214;
t321 = sin(t208) / 0.2e1;
t159 = t321 - t267 / 0.2e1;
t221 = qJ(5) + t228;
t210 = cos(t221);
t229 = qJ(1) + qJ(2);
t216 = sin(t229);
t218 = cos(t229);
t107 = t159 * t216 - t210 * t218;
t268 = cos(t285);
t274 = -t159 * t218 - t210 * t216;
t225 = qJD(1) + qJD(2);
t350 = t225 * t231;
t189 = qJD(3) + t350;
t375 = pkin(2) * t225;
t185 = t240 * t362 + t375;
t345 = t231 * t239;
t156 = t185 * t345;
t319 = t235 * t362;
t351 = t225 * t230;
t166 = pkin(8) * t351 + t319;
t284 = pkin(9) * t351 + t166;
t85 = -t234 * t284 + t156;
t74 = t189 * pkin(3) + t85;
t311 = t185 * t346;
t86 = t239 * t284 + t311;
t82 = t238 * t86;
t276 = -t233 * t74 - t82;
t389 = t398 * t225;
t372 = pkin(10) * t389;
t32 = -t276 + t372;
t232 = sin(qJ(5));
t327 = qJD(5) * t232;
t30 = t32 * t327;
t322 = cos(t208) / 0.2e1;
t132 = t397 * t225;
t237 = cos(qJ(5));
t72 = -t132 * t232 + t237 * t389;
t349 = t225 * t239;
t129 = (-pkin(3) * t349 - t185) * t230;
t77 = -pkin(4) * t389 + t129;
t395 = -g(1) * t107 - g(2) * t274 - g(3) * (t322 - t268 / 0.2e1) + t30 - t77 * t72;
t209 = sin(t221);
t252 = t268 / 0.2e1 + t322;
t102 = t209 * t216 - t218 * t252;
t105 = -t218 * t209 - t216 * t252;
t275 = -t237 * t132 - t232 * t389;
t394 = t77 * t275 + g(2) * t102 - g(3) * (t321 + t267 / 0.2e1) - g(1) * t105;
t223 = qJDD(1) + qJDD(2);
t382 = qJD(3) + qJD(4);
t244 = t382 * t272;
t53 = -t223 * t310 + (t223 * t344 + t225 * t244) * t230;
t317 = t230 * (-pkin(8) - pkin(9));
t288 = t234 * t317;
t393 = -qJD(3) * t288 + t396;
t257 = pkin(1) * (-t234 * t240 - t235 * t345);
t143 = qJD(1) * t257;
t202 = pkin(2) * t346;
t392 = (t239 * t317 - t202) * qJD(3) - t143;
t186 = qJD(4) + t189;
t125 = t132 * pkin(10);
t80 = t233 * t86;
t298 = t238 * t74 - t80;
t31 = -t125 + t298;
t29 = pkin(4) * t186 + t31;
t361 = t237 * t32;
t279 = -t232 * t29 - t361;
t354 = t223 * t231;
t187 = qJDD(3) + t354;
t184 = qJDD(4) + t187;
t377 = pkin(1) * t240;
t213 = qJDD(1) * t377;
t378 = pkin(1) * t235;
t318 = qJD(2) * t378;
t376 = pkin(2) * t223;
t152 = -qJD(1) * t318 + t213 + t376;
t136 = t152 * t345;
t324 = qJDD(1) * t235;
t332 = qJD(2) * t240;
t355 = t223 * t230;
t141 = pkin(8) * t355 + (qJD(1) * t332 + t324) * pkin(1);
t40 = t187 * pkin(3) + t136 + (-pkin(9) * t355 - t141) * t234 - t86 * qJD(3);
t331 = qJD(3) * t234;
t307 = t225 * t331;
t353 = t223 * t239;
t373 = pkin(9) * t230;
t309 = t239 * t141 + t152 * t346 + t185 * t305;
t48 = -t166 * t331 + t309;
t44 = (-t307 + t353) * t373 + t48;
t250 = qJD(4) * t276 - t233 * t44 + t238 * t40;
t52 = t223 * t397 + t389 * t382;
t6 = t184 * pkin(4) - t52 * pkin(10) + t250;
t328 = qJD(4) * t238;
t329 = qJD(4) * t233;
t296 = -t233 * t40 - t238 * t44 - t74 * t328 + t86 * t329;
t9 = -pkin(10) * t53 - t296;
t308 = -t232 * t9 + t237 * t6;
t251 = t279 * qJD(5) + t308;
t391 = t251 + t394;
t180 = qJD(5) + t186;
t390 = (-t180 * t32 - t6) * t232 + t395;
t366 = t275 * t72;
t17 = t275 ^ 2 - t72 ^ 2;
t326 = qJD(5) * t237;
t18 = -t132 * t327 - t232 * t53 + t237 * t52 + t326 * t389;
t13 = -t180 * t72 + t18;
t248 = qJD(5) * t275 - t232 * t52 - t237 * t53;
t14 = -t180 * t275 + t248;
t203 = pkin(2) * t345;
t220 = t231 * pkin(3);
t126 = t203 + t220 + t288;
t347 = t230 * t239;
t200 = pkin(9) * t347;
t335 = pkin(8) * t347 + t202;
t139 = t200 + t335;
t385 = -t126 * t328 + t139 * t329 - t392 * t233 + t393 * t238;
t338 = t233 * t126 + t238 * t139;
t384 = -t338 * qJD(4) + t393 * t233 + t392 * t238;
t383 = t230 * (-pkin(3) * t331 + t319);
t212 = pkin(2) + t377;
t182 = t212 * t346;
t188 = pkin(8) * t230 + t378;
t336 = t239 * t188 + t182;
t325 = qJD(3) - t189;
t380 = g(3) * t230 + t166 * t325;
t92 = t382 * (t339 - t344) * t230;
t379 = t92 * pkin(10);
t374 = pkin(3) * t239;
t367 = t398 * pkin(4);
t363 = t238 * t85 - t80;
t183 = t212 * t345;
t302 = -t188 - t373;
t100 = t234 * t302 + t183 + t220;
t108 = t200 + t336;
t360 = t233 * t100 + t238 * t108;
t359 = t132 * t389;
t171 = qJDD(5) + t184;
t358 = t171 * t233;
t211 = pkin(3) * t238 + pkin(4);
t357 = t211 * t171;
t224 = t230 ^ 2;
t356 = t225 ^ 2 * t224;
t352 = t224 * t239;
t348 = t230 * t234;
t343 = t233 * t237;
t340 = t234 * t239;
t337 = t239 * pkin(1) * t332 + t212 * t305;
t334 = g(1) * t218 + g(2) * t216;
t226 = t234 ^ 2;
t333 = -t239 ^ 2 + t226;
t199 = cos(t214) / 0.2e1;
t204 = cos(t304);
t323 = t204 / 0.2e1 + t199;
t320 = sin(t214) / 0.2e1;
t315 = t225 * t348;
t306 = t230 * t331;
t193 = pkin(3) * t306;
t93 = t244 * t230;
t78 = pkin(4) * t93 + t193;
t303 = -t185 - t375;
t301 = t231 * pkin(4) - pkin(10) * t397;
t300 = qJD(5) * t29 + t9;
t297 = -t233 * t85 - t82;
t295 = t238 * t100 - t108 * t233;
t294 = t238 * t126 - t139 * t233;
t293 = t187 + t354;
t292 = t189 + t350;
t291 = qJD(1) * (-qJD(2) + t225);
t290 = qJD(2) * (-qJD(1) - t225);
t289 = t231 * t318;
t283 = g(1) * t216 - g(2) * t218 + t213;
t54 = t294 + t301;
t91 = t93 * pkin(10);
t282 = -qJD(5) * t54 + t385 + t91;
t145 = t398 * pkin(10);
t56 = t145 + t338;
t281 = qJD(5) * t56 + t379 - t384;
t280 = sin(t304);
t176 = (-pkin(2) - t374) * t230;
t157 = (-t212 - t374) * t230;
t45 = t295 + t301;
t46 = t145 + t360;
t278 = -t232 * t46 + t237 * t45;
t277 = t232 * t45 + t237 * t46;
t89 = t232 * t397 - t237 * t398;
t90 = t232 * t398 + t237 * t397;
t164 = t320 - t280 / 0.2e1;
t217 = cos(t228);
t273 = -t164 * t218 - t216 * t217;
t117 = t164 * t216 - t217 * t218;
t270 = qJD(3) * (-t212 * t225 - t185);
t269 = -t230 * t319 + t78;
t75 = (qJD(3) * t302 - t289) * t234 + t337;
t255 = qJD(2) * t257;
t76 = t255 + (t239 * t302 - t182) * qJD(3);
t266 = t100 * t328 - t108 * t329 + t233 * t76 + t238 * t75;
t148 = -t216 * t239 - t218 * t346;
t150 = -t216 * t346 + t218 * t239;
t265 = -g(1) * t148 - g(2) * t150 + t152 * t352 + (-t234 * t141 + t136 + (-t166 * t239 - t311) * qJD(3)) * t231;
t264 = t225 * t319 + t376;
t261 = t212 * t223 - t225 * t318;
t147 = t216 * t234 - t218 * t345;
t149 = -t216 * t345 - t218 * t234;
t260 = -g(1) * t147 - g(2) * t149 - t48 * t231;
t34 = qJD(5) * t90 + t232 * t92 + t237 * t93;
t87 = t225 * t193 + (-pkin(3) * t353 - t152) * t230;
t35 = t53 * pkin(4) + t87;
t259 = -g(1) * t274 + g(2) * t107 + t251 * t231 + t77 * t34 + t35 * t89;
t258 = -g(1) * t273 + g(2) * t117 + t129 * t93 + t250 * t231 - t398 * t87;
t33 = -qJD(5) * t89 - t232 * t93 + t237 * t92;
t254 = -g(1) * t102 - g(2) * t105 - (t232 * t6 + t300 * t237 - t30) * t231 + t77 * t33 + t35 * t90;
t215 = sin(t228);
t112 = t215 * t216 - t218 * t323;
t115 = -t218 * t215 - t216 * t323;
t253 = -g(1) * t112 - g(2) * t115 + t129 * t92 + t231 * t296 + t397 * t87;
t249 = -qJD(4) * t360 - t233 * t75 + t238 * t76;
t246 = -g(1) * t117 - g(2) * t273 - g(3) * (t199 - t204 / 0.2e1) - t129 * t389 + t296;
t242 = -g(1) * t115 + g(2) * t112 - g(3) * (t320 + t280 / 0.2e1) - t129 * t132 + t250;
t241 = cos(qJ(1));
t236 = sin(qJ(1));
t194 = t230 * t318;
t170 = t187 * t231;
t162 = t184 * t231;
t161 = t171 * t231;
t155 = t194 + t193;
t138 = (t223 * t226 + 0.2e1 * t239 * t307) * t224;
t109 = t176 - t367;
t99 = t157 - t367;
t98 = 0.2e1 * (-qJD(3) * t225 * t333 + t223 * t340) * t224;
t97 = pkin(3) * t315 + pkin(4) * t132;
t64 = (t234 * t293 + t292 * t330) * t230;
t63 = (t239 * t293 - t292 * t331) * t230;
t62 = t194 + t78;
t55 = t132 ^ 2 - t389 ^ 2;
t42 = t132 * t186 - t53;
t41 = -t186 * t389 + t52;
t37 = -t125 + t363;
t36 = t297 - t372;
t26 = t132 * t92 + t397 * t52;
t25 = t184 * t398 - t186 * t93 - t231 * t53;
t24 = t184 * t397 + t186 * t92 + t231 * t52;
t16 = t249 - t379;
t15 = t266 - t91;
t12 = -t132 * t93 + t389 * t92 - t397 * t53 + t398 * t52;
t8 = -t171 * t89 - t180 * t34 + t231 * t248;
t7 = t171 * t90 + t18 * t231 + t180 * t33;
t4 = t18 * t90 - t275 * t33;
t3 = -t18 * t89 + t248 * t90 + t275 * t34 + t33 * t72;
t1 = [qJDD(1), g(1) * t236 - g(2) * t241, g(1) * t241 + g(2) * t236, t223, (t223 * t240 + t235 * t290) * pkin(1) + t283, ((-qJDD(1) - t223) * t235 + t240 * t290) * pkin(1) + t334, t138, t98, t64, t63, t170, (-t336 * qJD(3) + t255) * t189 + (-t188 * t234 + t183) * t187 + (t234 * t270 + t239 * t261) * t224 + t265, -((-qJD(3) * t188 - t289) * t234 + t337) * t189 - t336 * t187 + (t239 * t270 + (-t152 - t261) * t234) * t224 + t260, t26, t12, t24, t25, t162, -t155 * t389 + t157 * t53 + t184 * t295 + t186 * t249 + t258, t155 * t132 + t157 * t52 - t184 * t360 - t186 * t266 + t253, t4, t3, t7, t8, t161, (-qJD(5) * t277 - t232 * t15 + t237 * t16) * t180 + t278 * t171 - t62 * t72 - t99 * t248 + t259, -(qJD(5) * t278 + t237 * t15 + t232 * t16) * t180 - t277 * t171 - t62 * t275 + t99 * t18 + t254; 0, 0, 0, t223, t291 * t378 + t283, (t240 * t291 - t324) * pkin(1) + t334, t138, t98, t64, t63, t170, (-pkin(8) * t348 + t203) * t187 - t143 * t189 + t264 * t352 + (t224 * t234 * t303 - t189 * t335) * qJD(3) + t265, -t335 * t187 + (pkin(8) * t306 + t396) * t189 + (t303 * t330 + (-t152 - t264) * t234) * t224 + t260, t26, t12, t24, t25, t162, t176 * t53 + t294 * t184 + t384 * t186 + t383 * t389 + t258, -t132 * t383 + t176 * t52 - t338 * t184 + t385 * t186 + t253, t4, t3, t7, t8, t161, (-t232 * t56 + t237 * t54) * t171 - t109 * t248 - t269 * t72 + (t232 * t282 - t237 * t281) * t180 + t259, -(t232 * t54 + t237 * t56) * t171 + t109 * t18 - t269 * t275 + (t232 * t281 + t237 * t282) * t180 + t254; 0, 0, 0, 0, 0, 0, -t340 * t356, t333 * t356, (t223 * t234 + t325 * t349) * t230, (-t225 * t234 * t325 + t353) * t230, t187, -g(1) * t149 + g(2) * t147 + t136 - t380 * t239 + (-t141 + (t224 * t225 - t231 * t325) * t185) * t234, t224 * t185 * t349 + g(1) * t150 - g(2) * t148 + t156 * t189 + t380 * t234 - t309, -t359, t55, t41, t42, t184, -t297 * t186 + (t184 * t238 - t186 * t329 + t315 * t389) * pkin(3) + t242, t363 * t186 + (-t132 * t315 - t184 * t233 - t186 * t328) * pkin(3) + t246, t366, t17, t13, t14, t171, t237 * t357 - (-t232 * t37 + t237 * t36) * t180 + t97 * t72 + (-t232 * t358 + (-t232 * t238 - t343) * t180 * qJD(4)) * pkin(3) + ((-pkin(3) * t343 - t211 * t232) * t180 + t279) * qJD(5) + t308 + t394, t97 * t275 + (-t357 - t6 + (t36 - (-qJD(4) - qJD(5)) * t233 * pkin(3)) * t180) * t232 + (-pkin(3) * t358 + (-pkin(3) * t328 - qJD(5) * t211 + t37) * t180 - t300) * t237 + t395; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t359, t55, t41, t42, t184, -t186 * t276 + t242, t186 * t298 + t246, t366, t17, t13, t14, t171, -(-t232 * t31 - t361) * t180 + (t132 * t72 + t171 * t237 - t180 * t327) * pkin(4) + t391, (t180 * t31 - t300) * t237 + (t132 * t275 - t171 * t232 - t180 * t326) * pkin(4) + t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t366, t17, t13, t14, t171, -t180 * t279 + t391, (-t9 + (-qJD(5) + t180) * t29) * t237 + t390;];
tau_reg = t1;
