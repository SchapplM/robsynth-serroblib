% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:43:44
% EndTime: 2019-05-05 21:44:02
% DurationCPUTime: 7.39s
% Computational Cost: add. (16880->378), mult. (34219->464), div. (0->0), fcn. (22687->8), ass. (0->269)
t328 = pkin(7) + pkin(1);
t233 = sin(qJ(3));
t236 = cos(qJ(3));
t232 = sin(qJ(4));
t235 = cos(qJ(4));
t288 = qJD(1) * t236;
t210 = t232 * qJD(3) + t235 * t288;
t223 = t236 * qJDD(1);
t282 = qJD(1) * qJD(3);
t275 = t233 * t282;
t213 = t223 - t275;
t269 = -t235 * qJDD(3) + t232 * t213;
t177 = -t210 * qJD(4) - t269;
t208 = -t235 * qJD(3) + t232 * t288;
t255 = -t232 * qJDD(3) - t235 * t213;
t178 = -t208 * qJD(4) - t255;
t230 = sin(pkin(9));
t231 = cos(pkin(9));
t130 = t230 * t177 + t231 * t178;
t183 = t231 * t208 + t230 * t210;
t221 = t233 * qJD(1) + qJD(4);
t315 = t183 * t221;
t344 = t130 - t315;
t274 = t236 * t282;
t280 = t233 * qJDD(1);
t212 = -t274 - t280;
t207 = qJDD(4) - t212;
t185 = -t230 * t208 + t231 * t210;
t314 = t185 * t183;
t128 = -t314 - t207;
t304 = t230 * t128;
t182 = t185 ^ 2;
t329 = t221 ^ 2;
t340 = -t182 - t329;
t79 = t231 * t340 + t304;
t300 = t231 * t128;
t82 = -t230 * t340 + t300;
t46 = t232 * t79 - t235 * t82;
t37 = t233 * t46 + t236 * t344;
t51 = t232 * t82 + t235 * t79;
t405 = -qJ(2) * t51 - t328 * t37;
t404 = pkin(3) * t51;
t403 = pkin(8) * t51;
t401 = pkin(3) * t344 + pkin(8) * t46;
t165 = t221 * t185;
t270 = -t231 * t177 + t230 * t178;
t105 = -t270 + t165;
t330 = t183 ^ 2;
t158 = t330 - t329;
t87 = -t230 * t158 + t300;
t91 = -t231 * t158 - t304;
t400 = t233 * t105 + t236 * (t232 * t87 - t235 * t91);
t399 = pkin(4) * t79;
t398 = qJ(5) * t79;
t397 = qJ(5) * t82;
t339 = t182 - t330;
t342 = t165 + t270;
t61 = -t230 * t342 + t231 * t344;
t306 = t230 * t344;
t63 = t231 * t342 + t306;
t396 = -t233 * t339 + t236 * (t232 * t61 + t235 * t63);
t343 = t130 + t315;
t361 = t231 * t105 + t230 * t343;
t362 = t230 * t105 - t231 * t343;
t374 = t232 * t361 + t235 * t362;
t114 = -t330 - t182;
t375 = -t232 * t362 + t235 * t361;
t386 = -t236 * t114 + t233 * t375;
t395 = qJ(2) * t374 - t328 * t386;
t394 = t232 * t91 + t235 * t87;
t392 = pkin(3) * t374;
t391 = pkin(8) * t374;
t389 = -pkin(3) * t114 + pkin(8) * t375;
t388 = t232 * t63 - t235 * t61;
t337 = -t314 + t207;
t303 = t230 * t337;
t336 = -t329 - t330;
t346 = t231 * t336 - t303;
t120 = t231 * t337;
t347 = t230 * t336 + t120;
t359 = t232 * t346 + t235 * t347;
t360 = -t232 * t347 + t235 * t346;
t376 = t233 * t360 - t236 * t342;
t385 = qJ(2) * t359 - t328 * t376;
t384 = pkin(3) * t359;
t326 = pkin(4) * t362;
t383 = pkin(8) * t359;
t379 = qJ(5) * t362;
t378 = -pkin(3) * t342 + pkin(8) * t360;
t377 = -pkin(4) * t114 + qJ(5) * t361;
t160 = -t182 + t329;
t363 = t231 * t160 + t303;
t364 = -t230 * t160 + t120;
t373 = t232 * t364 + t235 * t363;
t372 = t236 * (-t232 * t363 + t235 * t364) + t233 * t343;
t325 = pkin(4) * t347;
t369 = qJ(5) * t346;
t368 = qJ(5) * t347;
t356 = qJ(6) * t344;
t190 = t210 * t208;
t338 = -t190 + t207;
t353 = t232 * t338;
t350 = t235 * t338;
t239 = qJD(1) ^ 2;
t345 = t328 * t239;
t287 = qJD(5) * t183;
t174 = -0.2e1 * t287;
t285 = qJD(6) * t221;
t341 = t174 + 0.2e1 * t285;
t198 = t221 * t208;
t149 = t178 + t198;
t135 = t183 * pkin(5) - t185 * qJ(6);
t281 = qJD(2) * qJD(1);
t225 = 0.2e1 * t281;
t227 = qJDD(1) * qJ(2);
t234 = sin(qJ(1));
t237 = cos(qJ(1));
t263 = t237 * g(1) + t234 * g(2);
t253 = -t227 + t263;
t249 = t225 - t253;
t257 = -t213 + t275;
t258 = -t212 + t274;
t144 = t258 * pkin(3) + t257 * pkin(8) + t249 - t345;
t272 = t234 * g(1) - t237 * g(2);
t260 = qJDD(2) - t272;
t291 = t239 * qJ(2);
t244 = t260 - t291;
t197 = -t328 * qJDD(1) + t244;
t187 = t236 * g(3) - t233 * t197;
t238 = qJD(3) ^ 2;
t264 = pkin(3) * t233 - pkin(8) * t236;
t250 = t239 * t264;
t153 = -t238 * pkin(3) + qJDD(3) * pkin(8) - t233 * t250 - t187;
t100 = -t235 * t144 + t232 * t153;
t69 = pkin(4) * t338 - qJ(5) * t149 - t100;
t101 = t232 * t144 + t235 * t153;
t205 = t208 ^ 2;
t259 = t221 * pkin(4) - t210 * qJ(5);
t71 = -t205 * pkin(4) + t177 * qJ(5) - t221 * t259 + t101;
t322 = t230 * t69 + t231 * t71;
t267 = t207 * qJ(6) - t183 * t135 + t322;
t335 = -t399 - pkin(5) * (t340 + t329) - qJ(6) * t128 + t267;
t145 = (qJD(4) - t221) * t210 + t269;
t247 = (-t183 * t230 - t185 * t231) * t221;
t313 = t221 * t230;
t157 = t185 * t313;
t312 = t221 * t231;
t279 = t183 * t312;
t261 = t157 - t279;
t334 = t232 * t261 + t235 * t247;
t251 = t230 * t270 + t279;
t262 = t183 * t313 - t231 * t270;
t333 = t232 * t251 + t235 * t262;
t296 = t233 * t207;
t332 = t236 * (-t232 * t247 + t235 * t261) + t296;
t278 = t233 * t314;
t331 = t236 * (-t232 * t262 + t235 * t251) - t278;
t206 = t210 ^ 2;
t271 = t230 * t71 - t231 * t69;
t286 = qJD(5) * t185;
t35 = t271 + 0.2e1 * t286;
t36 = t174 + t322;
t17 = t230 * t36 - t231 * t35;
t327 = pkin(4) * t17;
t324 = pkin(5) * t231;
t323 = t270 * pkin(5);
t186 = t233 * g(3) + t236 * t197;
t152 = qJDD(3) * pkin(3) + t238 * pkin(8) - t236 * t250 + t186;
t78 = t177 * pkin(4) + t205 * qJ(5) - t210 * t259 - qJDD(5) + t152;
t321 = t230 * t78;
t320 = t231 * t78;
t319 = t232 * t17;
t318 = t235 * t17;
t317 = qJ(6) * t231;
t316 = qJDD(1) * pkin(1);
t311 = t221 * t232;
t310 = t221 * t235;
t228 = t233 ^ 2;
t309 = t228 * t239;
t229 = t236 ^ 2;
t308 = t229 * t239;
t299 = t232 * t152;
t169 = t190 + t207;
t298 = t232 * t169;
t276 = t233 * t239 * t236;
t295 = t233 * (qJDD(3) + t276);
t294 = t235 * t152;
t293 = t235 * t169;
t292 = t236 * (qJDD(3) - t276);
t289 = t228 + t229;
t283 = qJD(4) + t221;
t277 = t233 * t190;
t273 = -qJ(6) * t230 - pkin(4);
t18 = t230 * t35 + t231 * t36;
t58 = t232 * t100 + t235 * t101;
t268 = (0.2e1 * qJD(5) + t135) * t185;
t266 = -t322 + t399;
t96 = t230 * t130 + t185 * t312;
t97 = t231 * t130 - t157;
t265 = t236 * (-t232 * t96 + t235 * t97) + t278;
t256 = t235 * t100 - t232 * t101;
t143 = t236 * t186 - t233 * t187;
t254 = qJ(2) + t264;
t252 = t267 + t341;
t248 = -t207 * pkin(5) - qJ(6) * t329 + qJDD(6) + t271;
t26 = -pkin(5) * t329 + t252;
t27 = t268 + t248;
t12 = t230 * t26 - t231 * t27;
t246 = pkin(4) * t12 - pkin(5) * t27 + qJ(6) * t26;
t245 = -pkin(5) * t343 + qJ(6) * t105 + t326;
t242 = -pkin(5) * t337 - qJ(6) * t336 + t248 - t325;
t241 = -pkin(5) * t165 + 0.2e1 * qJD(6) * t185 + t78;
t240 = t241 + t356;
t215 = t289 * qJDD(1);
t214 = t223 - 0.2e1 * t275;
t211 = 0.2e1 * t274 + t280;
t199 = -t244 + t316;
t196 = -t206 + t329;
t195 = t205 - t329;
t194 = t253 - 0.2e1 * t281 + t345;
t192 = -t295 + t236 * (-t238 - t308);
t191 = t233 * (-t238 - t309) + t292;
t189 = t206 - t205;
t188 = -t206 - t329;
t179 = -t329 - t205;
t176 = -0.2e1 * t286;
t175 = 0.2e1 * t287;
t167 = t205 + t206;
t150 = t283 * t208 + t255;
t148 = t178 - t198;
t146 = -t283 * t210 - t269;
t134 = -t232 * t188 - t293;
t133 = t235 * t188 - t298;
t123 = t235 * t179 - t353;
t122 = t232 * t179 + t350;
t112 = -t145 * t235 + t232 * t149;
t81 = t233 * t134 + t236 * t150;
t77 = t233 * t123 + t236 * t146;
t72 = t233 * t112 + t236 * t167;
t56 = t232 * t97 + t235 * t96;
t50 = -t320 - t398;
t49 = t236 * t152 + t233 * t58;
t48 = -t321 - t368;
t43 = -pkin(4) * t344 - t321 + t397;
t40 = t240 - t323;
t39 = -pkin(4) * t342 + t320 + t369;
t29 = (-t342 - t270) * pkin(5) + t240;
t28 = t241 - t323 + 0.2e1 * t356;
t23 = -qJ(6) * t114 + t27;
t22 = (-t114 - t329) * pkin(5) + t252;
t21 = -t230 * t29 - t317 * t342 - t368;
t20 = -pkin(5) * t306 + t231 * t28 + t398;
t19 = t231 * t29 + t273 * t342 + t369;
t16 = -t397 + t230 * t28 + (pkin(4) + t324) * t344;
t15 = pkin(4) * t78 + qJ(5) * t18;
t14 = -t17 - t379;
t13 = t230 * t27 + t231 * t26;
t11 = t18 + t377;
t10 = -t230 * t22 + t231 * t23 - t379;
t9 = t231 * t22 + t230 * t23 + t377;
t8 = t235 * t18 - t319;
t7 = t232 * t18 + t318;
t6 = t233 * t8 + t236 * t78;
t5 = -qJ(5) * t12 + (-pkin(5) * t230 + t317) * t40;
t4 = -t232 * t12 + t235 * t13;
t3 = t235 * t12 + t232 * t13;
t2 = qJ(5) * t13 + (-t273 + t324) * t40;
t1 = t233 * t4 + t236 * t40;
t24 = [0, 0, 0, 0, 0, qJDD(1), t272, t263, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t260 - 0.2e1 * t316, t225 + 0.2e1 * t227 - t263, pkin(1) * t199 + qJ(2) * (-t239 * pkin(1) + t249), -t257 * t236, -t236 * t211 - t233 * t214, t292 - t233 * (t238 - t308), t258 * t233, t236 * (-t238 + t309) - t295, 0, qJ(2) * t211 - t328 * t191 - t233 * t194, qJ(2) * t214 - t328 * t192 - t236 * t194, t328 * t215 - t289 * t291 - t143, -qJ(2) * t194 - t328 * t143, t236 * (t235 * t178 - t210 * t311) + t277, t236 * (t235 * t146 - t232 * t148) + t233 * t189, t236 * (-t232 * t196 + t350) + t233 * t149, t236 * (-t232 * t177 + t208 * t310) - t277, t236 * (t235 * t195 - t298) - t233 * t145, t296 + t236 * (-t208 * t235 + t210 * t232) * t221, t236 * (-pkin(8) * t122 - t299) - t233 * (-pkin(3) * t122 + t100) + qJ(2) * t122 - t328 * t77, t236 * (-pkin(8) * t133 - t294) - t233 * (-pkin(3) * t133 + t101) + qJ(2) * t133 - t328 * t81, t236 * t256 - t328 * t72 + t254 * (-t145 * t232 - t235 * t149), -t254 * t256 - t328 * t49, t265, -t396, t372, t331, t400, t332, t236 * (-t232 * t39 + t235 * t48 - t383) - t233 * (-t325 + t35 - t384) + t385, t236 * (-t232 * t43 + t235 * t50 - t403) - t233 * (t174 - t266 - t404) - t405, t236 * (-t232 * t11 + t235 * t14 - t391) - t233 * (-t326 - t392) + t395, t236 * (-pkin(8) * t7 - qJ(5) * t318 - t232 * t15) - t233 * (-pkin(3) * t7 - t327) + qJ(2) * t7 - t328 * t6, t265, t372, t396, t332, -t400, t331, t236 * (-t232 * t19 + t235 * t21 - t383) - t233 * (t242 + t268 - t384) + t385, t236 * (t235 * t10 - t232 * t9 - t391) - t233 * (-t245 - t392) + t395, t236 * (-t232 * t16 + t235 * t20 + t403) - t233 * (t175 - 0.2e1 * t285 - t335 + t404) + t405, t236 * (-pkin(8) * t3 - t232 * t2 + t235 * t5) - t233 * (-pkin(3) * t3 - t246) + qJ(2) * t3 - t328 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t239, -t199, 0, 0, 0, 0, 0, 0, t191, t192, -t215, t143, 0, 0, 0, 0, 0, 0, t77, t81, t72, t49, 0, 0, 0, 0, 0, 0, t376, -t37, t386, t6, 0, 0, 0, 0, 0, 0, t376, t386, t37, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, (-t228 + t229) * t239, t223, -t276, -t280, qJDD(3), t186, t187, 0, 0, t232 * t178 + t210 * t310, t232 * t146 + t235 * t148, t235 * t196 + t353, t235 * t177 + t208 * t311, t232 * t195 + t293, (-t208 * t232 - t210 * t235) * t221, pkin(3) * t146 + pkin(8) * t123 + t294, pkin(3) * t150 + pkin(8) * t134 - t299, pkin(3) * t167 + pkin(8) * t112 + t58, pkin(3) * t152 + pkin(8) * t58, t56, -t388, t373, t333, -t394, t334, t232 * t48 + t235 * t39 + t378, t232 * t50 + t235 * t43 - t401, t235 * t11 + t232 * t14 + t389, pkin(3) * t78 + pkin(8) * t8 - qJ(5) * t319 + t235 * t15, t56, t373, t388, t334, t394, t333, t235 * t19 + t232 * t21 + t378, t232 * t10 + t235 * t9 + t389, t235 * t16 + t232 * t20 + t401, pkin(3) * t40 + pkin(8) * t4 + t235 * t2 + t232 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, t189, t149, -t190, -t145, t207, -t100, -t101, 0, 0, t314, t339, t343, -t314, t105, t207, t176 - t271 + t325, t175 + t266, t326, t327, t314, t343, -t339, t207, -t105, -t314, -t185 * t135 + t176 - t242, t245, t335 + t341, t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, t344, t114, -t78, 0, 0, 0, 0, 0, 0, t342, t114, -t344, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, t343, t340, t27;];
tauJ_reg  = t24;
