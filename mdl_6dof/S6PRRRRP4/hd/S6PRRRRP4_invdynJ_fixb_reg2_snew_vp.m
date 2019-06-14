% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:55:48
% EndTime: 2019-05-05 09:56:11
% DurationCPUTime: 9.31s
% Computational Cost: add. (21659->429), mult. (42270->567), div. (0->0), fcn. (31056->12), ass. (0->280)
t249 = sin(pkin(6));
t250 = cos(pkin(6));
t255 = sin(qJ(3));
t256 = sin(qJ(2));
t259 = cos(qJ(3));
t260 = cos(qJ(2));
t254 = sin(qJ(4));
t258 = cos(qJ(4));
t308 = qJD(2) * t255;
t216 = -t258 * qJD(3) + t254 * t308;
t304 = qJD(2) * qJD(3);
t297 = t259 * t304;
t303 = t255 * qJDD(2);
t220 = t297 + t303;
t279 = -t254 * qJDD(3) - t258 * t220;
t189 = -t216 * qJD(4) - t279;
t217 = t254 * qJD(3) + t258 * t308;
t253 = sin(qJ(5));
t257 = cos(qJ(5));
t194 = t257 * t216 + t253 * t217;
t280 = t258 * qJDD(3) - t254 * t220;
t270 = -t217 * qJD(4) + t280;
t263 = -t194 * qJD(5) + t257 * t189 + t253 * t270;
t238 = t259 * qJD(2) - qJD(4);
t232 = -qJD(5) + t238;
t326 = t194 * t232;
t356 = t263 + t326;
t196 = -t253 * t216 + t257 * t217;
t157 = t196 * t194;
t241 = t255 * t304;
t302 = t259 * qJDD(2);
t221 = -t241 + t302;
t215 = -qJDD(4) + t221;
t212 = -qJDD(5) + t215;
t360 = -t157 + t212;
t313 = t257 * t360;
t193 = t196 ^ 2;
t345 = t232 ^ 2;
t358 = -t193 - t345;
t105 = -t253 * t358 + t313;
t319 = t253 * t360;
t79 = -t257 * t358 - t319;
t70 = t258 * t105 + t254 * t79;
t45 = t255 * t356 + t259 * t70;
t68 = t254 * t105 - t258 * t79;
t420 = (t256 * t45 - t260 * t68) * t249 + t250 * (t255 * t70 - t259 * t356);
t419 = pkin(2) * t68 - pkin(8) * t45;
t415 = pkin(3) * t68;
t414 = pkin(9) * t68;
t412 = pkin(3) * t356 - pkin(9) * t70;
t347 = t194 ^ 2;
t172 = t347 - t345;
t111 = t253 * t172 - t313;
t115 = t257 * t172 + t319;
t179 = t196 * t232;
t294 = -t253 * t189 + t257 * t270;
t275 = t196 * qJD(5) - t294;
t94 = t179 + t275;
t410 = t255 * (t254 * t111 - t258 * t115) - t259 * t94;
t357 = t193 - t347;
t359 = -t179 + t275;
t60 = -t253 * t359 + t257 * t356;
t335 = t253 * t356;
t61 = t257 * t359 + t335;
t409 = t255 * (t254 * t60 + t258 * t61) + t259 * t357;
t408 = pkin(4) * t79;
t407 = pkin(10) * t79;
t406 = pkin(10) * t105;
t139 = t212 + t157;
t318 = t253 * t139;
t354 = -t345 - t347;
t363 = t257 * t354 + t318;
t135 = t257 * t139;
t364 = t253 * t354 - t135;
t376 = t254 * t363 + t258 * t364;
t377 = -t254 * t364 + t258 * t363;
t388 = t255 * t359 + t259 * t377;
t405 = -pkin(2) * t376 + pkin(8) * t388;
t404 = t254 * t61 - t258 * t60;
t402 = t258 * t111 + t254 * t115;
t173 = -t193 + t345;
t390 = t257 * t173 - t318;
t391 = -t253 * t173 - t135;
t400 = t254 * t391 + t258 * t390;
t399 = t250 * (t255 * t377 - t259 * t359) + (t256 * t388 - t260 * t376) * t249;
t355 = -t326 + t263;
t398 = t255 * (-t254 * t390 + t258 * t391) - t259 * t355;
t396 = pkin(3) * t376;
t395 = pkin(9) * t376;
t389 = -pkin(3) * t359 + pkin(9) * t377;
t123 = -t347 - t193;
t386 = pkin(3) * t123;
t385 = pkin(4) * t123;
t384 = pkin(4) * t364;
t383 = pkin(10) * t363;
t382 = pkin(10) * t364;
t380 = qJ(6) * t356;
t328 = sin(pkin(11));
t329 = cos(pkin(11));
t273 = g(1) * t328 - g(2) * t329;
t309 = -g(3) + qJDD(1);
t265 = -t249 * t273 + t250 * t309;
t201 = t259 * t265;
t227 = -g(1) * t329 - g(2) * t328;
t269 = t250 * t273;
t365 = t249 * t309 + t269;
t185 = t260 * t227 + t256 * t365;
t261 = qJD(2) ^ 2;
t169 = -t261 * pkin(2) + qJDD(2) * pkin(8) + t185;
t287 = -t259 * pkin(3) - t255 * pkin(9);
t292 = t261 * t287 + t169;
t344 = qJD(3) ^ 2;
t126 = -qJDD(3) * pkin(3) - t344 * pkin(9) + t255 * t292 - t201;
t202 = -t238 * pkin(4) - t217 * pkin(10);
t346 = t216 ^ 2;
t77 = -t270 * pkin(4) - t346 * pkin(10) + t217 * t202 + t126;
t381 = pkin(5) * t275 - t380 + t77;
t379 = t255 * t123;
t378 = t259 * t123;
t324 = t217 * t216;
t272 = -t215 - t324;
t362 = t254 * t272;
t361 = t258 * t272;
t205 = t216 * t238;
t162 = t189 - t205;
t158 = (qJD(4) + t238) * t217 - t280;
t274 = (t194 * t253 + t196 * t257) * t232;
t323 = t232 * t253;
t170 = t196 * t323;
t322 = t232 * t257;
t301 = t194 * t322;
t285 = -t170 + t301;
t351 = t254 * t285 + t258 * t274;
t276 = t253 * t275 - t301;
t286 = -t194 * t323 - t257 * t275;
t350 = t254 * t276 + t258 * t286;
t349 = t255 * (-t254 * t274 + t258 * t285) + t259 * t212;
t300 = t259 * t157;
t348 = t255 * (-t254 * t286 + t258 * t276) + t300;
t214 = t217 ^ 2;
t236 = t238 ^ 2;
t264 = t255 * t265;
t127 = -pkin(3) * t344 + qJDD(3) * pkin(9) + t259 * t292 + t264;
t284 = t256 * t227 - t260 * t365;
t168 = -qJDD(2) * pkin(2) - t261 * pkin(8) + t284;
t282 = -t221 + t241;
t283 = t220 + t297;
t134 = pkin(3) * t282 - pkin(9) * t283 + t168;
t75 = t254 * t127 - t258 * t134;
t67 = t272 * pkin(4) - t162 * pkin(10) - t75;
t76 = t258 * t127 + t254 * t134;
t73 = -pkin(4) * t346 + pkin(10) * t270 + t238 * t202 + t76;
t37 = t253 * t73 - t257 * t67;
t38 = t253 * t67 + t257 * t73;
t18 = t253 * t38 - t257 * t37;
t343 = pkin(4) * t18;
t89 = t257 * t355;
t95 = (-qJD(5) - t232) * t196 + t294;
t58 = t253 * t95 - t89;
t342 = pkin(4) * t58;
t341 = pkin(5) * t257;
t307 = qJD(6) * t232;
t223 = -0.2e1 * t307;
t153 = t194 * pkin(5) - t196 * qJ(6);
t293 = -t212 * qJ(6) - t194 * t153 + t38;
t277 = -pkin(5) * t345 + t293;
t27 = t223 + t277;
t29 = t212 * pkin(5) - qJ(6) * t345 + t196 * t153 + qJDD(6) + t37;
t340 = -pkin(5) * t29 + qJ(6) * t27;
t339 = -pkin(5) * t355 - qJ(6) * t94;
t338 = t253 * t77;
t336 = t253 * t355;
t334 = t254 * t18;
t333 = t257 * t77;
t331 = t258 * t18;
t327 = qJ(6) * t257;
t321 = t238 * t254;
t320 = t238 * t258;
t316 = t254 * t126;
t181 = t215 - t324;
t315 = t254 * t181;
t237 = t255 * t261 * t259;
t229 = qJDD(3) + t237;
t314 = t255 * t229;
t312 = t258 * t126;
t311 = t258 * t181;
t228 = -t237 + qJDD(3);
t310 = t259 * t228;
t306 = qJD(4) - t238;
t299 = t259 * t324;
t296 = -qJ(6) * t253 - pkin(4);
t19 = t253 * t37 + t257 * t38;
t48 = t254 * t75 + t258 * t76;
t151 = t255 * t169 - t201;
t152 = t259 * t169 + t264;
t101 = t255 * t151 + t259 * t152;
t13 = t253 * t27 - t257 * t29;
t291 = pkin(4) * t13 + t340;
t59 = -t253 * t94 - t89;
t290 = pkin(4) * t59 + t339;
t289 = -t38 - t408;
t87 = -t196 * t322 + t253 * t263;
t88 = t257 * t263 + t170;
t288 = t255 * (-t254 * t87 + t258 * t88) - t300;
t47 = t254 * t76 - t258 * t75;
t281 = -pkin(2) + t287;
t278 = -t37 + t384;
t271 = -pkin(5) * t358 - qJ(6) * t360 + t277;
t268 = t271 + t408;
t267 = -pkin(5) * t139 + qJ(6) * t354 - t29;
t266 = t267 + t384;
t262 = 0.2e1 * qJD(6) * t196 - t381;
t246 = t259 ^ 2;
t245 = t255 ^ 2;
t244 = t246 * t261;
t242 = t245 * t261;
t235 = -t244 - t344;
t234 = -t242 - t344;
t225 = t242 + t244;
t224 = (t245 + t246) * qJDD(2);
t222 = -0.2e1 * t241 + t302;
t219 = 0.2e1 * t297 + t303;
t204 = -t214 + t236;
t203 = -t236 + t346;
t200 = -t255 * t234 - t310;
t199 = t259 * t235 - t314;
t198 = t214 - t346;
t197 = -t214 - t236;
t190 = -t236 - t346;
t180 = t214 + t346;
t163 = t216 * t306 + t279;
t161 = t189 + t205;
t159 = -t217 * t306 + t280;
t150 = -t254 * t197 + t311;
t149 = t258 * t197 + t315;
t144 = t258 * t190 - t362;
t143 = t254 * t190 + t361;
t117 = -t158 * t258 + t254 * t162;
t116 = -t158 * t254 - t258 * t162;
t107 = t259 * t150 - t255 * t163;
t102 = t259 * t144 - t255 * t159;
t78 = t259 * t117 - t255 * t180;
t63 = -t257 * t94 + t336;
t62 = t257 * t95 + t336;
t56 = t333 + t407;
t55 = t254 * t88 + t258 * t87;
t49 = t338 - t382;
t44 = -pkin(4) * t356 + t338 + t406;
t43 = t255 * t126 + t259 * t48;
t40 = -pkin(4) * t359 - t333 + t383;
t39 = (-pkin(5) * t232 - 0.2e1 * qJD(6)) * t196 + t381;
t35 = -t254 * t59 + t258 * t63;
t34 = -t254 * t58 + t258 * t62;
t33 = t254 * t63 + t258 * t59;
t32 = t254 * t62 + t258 * t58;
t31 = (-t359 + t179) * pkin(5) + t262;
t30 = pkin(5) * t179 + t262 + t380;
t25 = t259 * t35 + t379;
t24 = t259 * t34 + t379;
t23 = -qJ(6) * t123 + t29;
t22 = t223 + (-t123 - t345) * pkin(5) + t293;
t21 = -t253 * t31 - t327 * t359 - t382;
t20 = -pkin(5) * t335 + t257 * t30 - t407;
t17 = t257 * t31 + t296 * t359 + t383;
t16 = -t406 + t253 * t30 + (pkin(4) + t341) * t356;
t15 = -pkin(4) * t77 + pkin(10) * t19;
t14 = t253 * t29 + t257 * t27;
t12 = -pkin(10) * t58 - t18;
t11 = pkin(10) * t62 + t19 - t385;
t10 = -pkin(10) * t59 - t253 * t22 + t257 * t23;
t9 = pkin(10) * t63 + t257 * t22 + t253 * t23 - t385;
t8 = t258 * t19 - t334;
t7 = t254 * t19 + t331;
t6 = t255 * t77 + t259 * t8;
t5 = -pkin(10) * t13 + (pkin(5) * t253 - t327) * t39;
t4 = -t254 * t13 + t258 * t14;
t3 = t258 * t13 + t254 * t14;
t2 = pkin(10) * t14 + (t296 - t341) * t39;
t1 = t255 * t39 + t259 * t4;
t26 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t309, 0, 0, 0, 0, 0, 0, (qJDD(2) * t260 - t256 * t261) * t249, (-qJDD(2) * t256 - t260 * t261) * t249, 0, t250 ^ 2 * t309 + (t256 * t185 - t260 * t284 - t269) * t249, 0, 0, 0, 0, 0, 0, t250 * (t259 * t229 + t255 * t235) + (t256 * t199 + t260 * t222) * t249, t250 * (-t255 * t228 + t259 * t234) + (t256 * t200 - t260 * t219) * t249, (t224 * t256 + t225 * t260) * t249, t250 * (-t259 * t151 + t255 * t152) + (t256 * t101 - t260 * t168) * t249, 0, 0, 0, 0, 0, 0, t250 * (t255 * t144 + t259 * t159) + (t256 * t102 - t260 * t143) * t249, t250 * (t255 * t150 + t259 * t163) + (t256 * t107 - t260 * t149) * t249, t250 * (t255 * t117 + t259 * t180) + (-t260 * t116 + t256 * t78) * t249, t250 * (-t259 * t126 + t255 * t48) + (t256 * t43 - t260 * t47) * t249, 0, 0, 0, 0, 0, 0, t399, t420, t250 * (t255 * t34 - t378) + (t256 * t24 - t260 * t32) * t249, t250 * (t255 * t8 - t259 * t77) + (t256 * t6 - t260 * t7) * t249, 0, 0, 0, 0, 0, 0, t399, t250 * (t255 * t35 - t378) + (t256 * t25 - t260 * t33) * t249, -t420, t250 * (t255 * t4 - t259 * t39) + (t256 * t1 - t260 * t3) * t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t284, -t185, 0, 0, t283 * t255, t259 * t219 + t255 * t222, t314 + t259 * (-t242 + t344), -t282 * t259, t255 * (t244 - t344) + t310, 0, pkin(2) * t222 + pkin(8) * t199 - t259 * t168, -pkin(2) * t219 + pkin(8) * t200 + t255 * t168, pkin(2) * t225 + pkin(8) * t224 + t101, -pkin(2) * t168 + pkin(8) * t101, t255 * (t258 * t189 + t217 * t321) - t299, t255 * (t258 * t159 - t254 * t161) - t259 * t198, t255 * (-t254 * t204 + t361) - t259 * t162, t255 * (-t216 * t320 - t254 * t270) + t299, t255 * (t258 * t203 + t315) + t259 * t158, t259 * t215 + t255 * (t216 * t258 - t217 * t254) * t238, t255 * (-pkin(9) * t143 + t316) + t259 * (-pkin(3) * t143 + t75) - pkin(2) * t143 + pkin(8) * t102, t255 * (-pkin(9) * t149 + t312) + t259 * (-pkin(3) * t149 + t76) - pkin(2) * t149 + pkin(8) * t107, pkin(8) * t78 + t116 * t281 - t255 * t47, pkin(8) * t43 + t281 * t47, t288, -t409, t398, t348, -t410, t349, t255 * (-t254 * t40 + t258 * t49 - t395) + t259 * (-t278 - t396) + t405, t255 * (-t254 * t44 + t258 * t56 - t414) + t259 * (-t289 - t415) - t419, t255 * (-pkin(9) * t32 - t254 * t11 + t258 * t12) + t259 * (-pkin(3) * t32 - t342) - pkin(2) * t32 + pkin(8) * t24, t255 * (-pkin(9) * t7 - pkin(10) * t331 - t254 * t15) + t259 * (-pkin(3) * t7 - t343) - pkin(2) * t7 + pkin(8) * t6, t288, t398, t409, t349, t410, t348, t255 * (-t254 * t17 + t258 * t21 - t395) + t259 * (-t266 - t396) + t405, t255 * (-pkin(9) * t33 + t258 * t10 - t254 * t9) + t259 * (-pkin(3) * t33 - t290) - pkin(2) * t33 + pkin(8) * t25, t255 * (-t254 * t16 + t258 * t20 + t414) + t259 * (-t268 + 0.2e1 * t307 + t415) + t419, t255 * (-pkin(9) * t3 - t254 * t2 + t258 * t5) + t259 * (-pkin(3) * t3 - t291) - pkin(2) * t3 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, t242 - t244, t303, t237, t302, qJDD(3), -t151, -t152, 0, 0, t254 * t189 - t217 * t320, t254 * t159 + t258 * t161, t258 * t204 + t362, -t216 * t321 + t258 * t270, t254 * t203 - t311, (t216 * t254 + t217 * t258) * t238, pkin(3) * t159 + pkin(9) * t144 - t312, pkin(3) * t163 + pkin(9) * t150 + t316, pkin(3) * t180 + pkin(9) * t117 + t48, -pkin(3) * t126 + pkin(9) * t48, t55, -t404, t400, t350, t402, t351, t254 * t49 + t258 * t40 + t389, t254 * t56 + t258 * t44 - t412, pkin(9) * t34 + t258 * t11 + t254 * t12 - t386, -pkin(3) * t77 + pkin(9) * t8 - pkin(10) * t334 + t258 * t15, t55, t400, t404, t351, -t402, t350, t258 * t17 + t254 * t21 + t389, pkin(9) * t35 + t254 * t10 + t258 * t9 - t386, t258 * t16 + t254 * t20 + t412, -pkin(3) * t39 + pkin(9) * t4 + t258 * t2 + t254 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t198, t162, -t324, -t158, -t215, -t75, -t76, 0, 0, t157, t357, t355, -t157, -t94, -t212, t278, t289, t342, t343, t157, t355, -t357, -t212, t94, -t157, t266, t290, t223 + t268, t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t357, t355, -t157, -t94, -t212, -t37, -t38, 0, 0, t157, t355, -t357, -t212, t94, -t157, t267, t339, t223 + t271, t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t355, t358, t29;];
tauJ_reg  = t26;
