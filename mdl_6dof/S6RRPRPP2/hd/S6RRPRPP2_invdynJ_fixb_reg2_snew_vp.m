% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRPP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:25:56
% EndTime: 2019-05-06 12:26:19
% DurationCPUTime: 9.89s
% Computational Cost: add. (19384->494), mult. (44432->589), div. (0->0), fcn. (31185->8), ass. (0->289)
t274 = sin(qJ(2));
t276 = cos(qJ(2));
t270 = sin(pkin(9));
t271 = cos(pkin(9));
t327 = qJD(1) * t276;
t328 = qJD(1) * t274;
t248 = t270 * t327 + t271 * t328;
t273 = sin(qJ(4));
t275 = cos(qJ(4));
t227 = qJD(2) * t273 + t248 * t275;
t224 = t227 ^ 2;
t225 = -t275 * qJD(2) + t248 * t273;
t384 = t225 ^ 2;
t150 = -t384 - t224;
t438 = t150 * t271;
t246 = t270 * t328 - t271 * t327;
t242 = qJD(4) + t246;
t261 = t274 * qJDD(1);
t318 = qJD(1) * qJD(2);
t313 = t276 * t318;
t253 = t261 + t313;
t314 = t274 * t318;
t317 = t276 * qJDD(1);
t293 = -t314 + t317;
t215 = t271 * t253 + t270 * t293;
t309 = -t275 * qJDD(2) + t273 * t215;
t125 = (qJD(4) - t242) * t227 + t309;
t298 = -t273 * qJDD(2) - t275 * t215;
t288 = qJD(4) * t225 + t298;
t340 = t225 * t242;
t397 = t340 - t288;
t361 = t397 * t273;
t78 = t125 * t275 - t361;
t53 = t270 * t78 + t438;
t439 = t150 * t270;
t56 = t271 * t78 - t439;
t490 = pkin(7) * (t274 * t53 - t276 * t56);
t489 = qJ(3) * t53;
t488 = qJ(3) * t56;
t134 = t340 + t288;
t383 = t242 ^ 2;
t159 = -t383 - t224;
t183 = t227 * t225;
t307 = t253 * t270 - t271 * t293;
t212 = qJDD(4) + t307;
t402 = t183 + t212;
t415 = t275 * t402;
t454 = t159 * t273 + t415;
t463 = t271 * t454;
t457 = t134 * t270 + t463;
t414 = t402 * t273;
t97 = t159 * t275 - t414;
t469 = pkin(2) * t97;
t487 = qJ(3) * t457 + t469;
t464 = t270 * t454;
t458 = -t271 * t134 + t464;
t470 = pkin(1) * t97;
t486 = t470 + pkin(7) * (-t274 * t458 + t276 * t457);
t399 = t224 - t384;
t358 = t134 * t273;
t164 = qJD(4) * t227 + t309;
t338 = t242 * t227;
t404 = t164 + t338;
t410 = t404 * t275;
t81 = -t410 + t358;
t472 = t276 * (t270 * t81 - t271 * t399) + t274 * (t399 * t270 + t271 * t81);
t444 = pkin(3) * t150;
t485 = -pkin(2) * t53 - pkin(8) * t78 - t444;
t395 = t384 - t383;
t108 = -t275 * t395 + t414;
t405 = t164 - t338;
t483 = t274 * (t108 * t271 + t405 * t270) + t276 * (t108 * t270 - t405 * t271);
t193 = t224 - t383;
t403 = -t183 + t212;
t412 = t403 * t275;
t107 = -t193 * t273 - t412;
t482 = t274 * (t107 * t271 - t270 * t397) + t276 * (t107 * t270 + t271 * t397);
t468 = pkin(3) * t97;
t466 = pkin(8) * t97;
t480 = qJ(3) * t458;
t465 = pkin(8) * t454;
t477 = pkin(2) * t458 - pkin(3) * t134 + t465;
t360 = t397 * t275;
t74 = t125 * t273 + t360;
t467 = pkin(8) * t74;
t396 = -t383 - t384;
t413 = t403 * t273;
t433 = t275 * t396 - t413;
t451 = t270 * t433 - t271 * t404;
t456 = qJ(3) * t451;
t453 = pkin(2) * t451 - pkin(3) * t404 + pkin(8) * t433;
t434 = t273 * t396 + t412;
t450 = t404 * t270 + t271 * t433;
t452 = -pkin(2) * t434 + qJ(3) * t450;
t449 = pkin(7) * (-t274 * t451 + t276 * t450) - pkin(1) * t434;
t443 = pkin(3) * t434;
t442 = pkin(8) * t434;
t455 = t193 * t275 - t413;
t357 = t134 * t275;
t411 = t404 * t273;
t428 = t357 + t411;
t431 = t273 * t395 + t415;
t441 = qJ(5) * t134;
t440 = qJ(5) * t150;
t370 = qJ(5) * t396;
t435 = pkin(4) * t403 + t370;
t425 = qJ(5) * t402;
t424 = qJ(6) * t397;
t213 = t248 * t246;
t393 = -t213 + qJDD(2);
t423 = t270 * t393;
t421 = t271 * t393;
t326 = qJD(2) * t248;
t184 = t307 + t326;
t176 = t270 * t183;
t177 = t271 * t183;
t337 = t242 * t273;
t189 = t227 * t337;
t330 = -t275 * t288 - t189;
t387 = t276 * (t270 * t330 - t177) + t274 * (t271 * t330 + t176);
t321 = qJD(5) * t242;
t230 = 0.2e1 * t321;
t320 = qJD(6) * t225;
t401 = 0.2e1 * t320 + t230;
t231 = -0.2e1 * t321;
t400 = -0.2e1 * t320 + t231;
t239 = qJD(2) * t246;
t186 = -t239 + t215;
t268 = t276 ^ 2;
t277 = qJD(1) ^ 2;
t265 = t268 * t277;
t381 = qJD(2) ^ 2;
t398 = -t265 - t381;
t378 = sin(qJ(1));
t379 = cos(qJ(1));
t292 = g(1) * t379 + g(2) * t378;
t369 = qJDD(1) * pkin(7);
t283 = -t277 * pkin(1) - t292 + t369;
t229 = -t274 * g(3) + t276 * t283;
t174 = pkin(2) * t398 + qJ(3) * t317 + t229;
t281 = t274 * t283;
t334 = t274 * t277;
t278 = -t281 - t253 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t334 + qJ(3) * t318 - g(3)) * t276;
t104 = -0.2e1 * qJD(3) * t246 + t271 * t174 + t270 * t278;
t300 = -pkin(5) * t242 - qJ(6) * t227;
t394 = t227 * t300 + qJDD(6);
t392 = t164 * pkin(4) + t441;
t391 = -t164 * pkin(5) + t394;
t203 = pkin(3) * t246 - pkin(8) * t248;
t86 = -pkin(3) * t381 + qJDD(2) * pkin(8) - t203 * t246 + t104;
t306 = g(1) * t378 - t379 * g(2);
t291 = qJDD(1) * pkin(1) + t306;
t180 = (qJ(3) * t268 + pkin(7)) * t277 + pkin(2) * t293 - qJDD(3) - (qJD(2) * pkin(2) - qJ(3) * t328) * t328 + t291;
t95 = pkin(3) * t184 - t186 * pkin(8) - t180;
t50 = t273 * t86 - t275 * t95;
t295 = -t212 * pkin(4) - qJ(5) * t383 + qJDD(5) + t50;
t286 = -t212 * pkin(5) + t295 - t424;
t175 = pkin(4) * t225 - qJ(5) * t227;
t312 = -pkin(5) * t225 - t175;
t34 = (-0.2e1 * qJD(6) - t312) * t227 + t286;
t51 = t273 * t95 + t275 * t86;
t302 = -pkin(4) * t383 + t212 * qJ(5) - t225 * t175 + t51;
t289 = -pkin(5) * t384 + t164 * qJ(6) + t242 * t300 + t302;
t35 = t289 + t401;
t380 = pkin(4) + pkin(5);
t390 = qJ(5) * t35 - t34 * t380;
t389 = qJ(5) * t125 + t380 * t397;
t388 = -t159 * t380 + t289 + t425;
t336 = t242 * t275;
t316 = t225 * t336;
t296 = t164 * t273 + t316;
t386 = t274 * (t271 * t296 - t176) + t276 * (t270 * t296 + t177);
t202 = t271 * t212;
t303 = t189 - t316;
t341 = t212 * t270;
t385 = t274 * (t271 * t303 + t341) + t276 * (t270 * t303 - t202);
t244 = t246 ^ 2;
t245 = t248 ^ 2;
t382 = 0.2e1 * t227;
t377 = pkin(3) * t270;
t310 = t270 * t174 - t271 * t278;
t297 = -qJDD(2) * pkin(3) - t381 * pkin(8) + t310;
t85 = (0.2e1 * qJD(3) + t203) * t248 + t297;
t374 = t273 * t85;
t322 = qJD(3) * t248;
t103 = t310 + 0.2e1 * t322;
t60 = -t103 * t271 + t104 * t270;
t373 = t274 * t60;
t372 = t275 * t85;
t346 = t180 * t270;
t345 = t180 * t271;
t206 = qJDD(2) + t213;
t343 = t206 * t270;
t342 = t206 * t271;
t339 = t227 * t175;
t260 = t276 * t334;
t335 = t274 * (qJDD(2) + t260);
t333 = t276 * (qJDD(2) - t260);
t325 = qJD(2) * t270;
t324 = qJD(2) * t271;
t315 = -pkin(3) * t271 - pkin(2);
t33 = t273 * t50 + t275 * t51;
t61 = t103 * t270 + t271 * t104;
t228 = t276 * g(3) + t281;
t308 = t274 * t228 + t276 * t229;
t44 = t230 + t302;
t45 = t295 + t339;
t305 = -pkin(4) * t45 + qJ(5) * t44;
t116 = t227 * t336 - t273 * t288;
t304 = -t275 * t164 + t225 * t337;
t301 = -pkin(4) * t397 - qJ(5) * t405;
t32 = t273 * t51 - t275 * t50;
t185 = -t307 + t326;
t294 = (-t225 * t273 - t227 * t275) * t242;
t254 = -0.2e1 * t314 + t317;
t287 = -pkin(4) * t159 + t302 + t425;
t284 = qJD(6) * t382 - t286;
t238 = -0.2e1 * t322;
t282 = qJD(5) * t382 - t248 * t203 + t238 - t297 - t392;
t280 = -t45 + t435;
t47 = (pkin(4) * t242 - 0.2e1 * qJD(5)) * t227 + t85 + t392;
t41 = (-t404 - t338) * pkin(4) + t282;
t40 = -pkin(4) * t338 + t282 - t441;
t267 = t274 ^ 2;
t263 = t267 * t277;
t252 = t261 + 0.2e1 * t313;
t250 = t277 * pkin(7) + t291;
t234 = -t245 - t381;
t233 = -t245 + t381;
t232 = t244 - t381;
t204 = -t381 - t244;
t187 = t239 + t215;
t181 = -t244 - t245;
t167 = -t234 * t270 - t342;
t166 = t234 * t271 - t343;
t149 = t204 * t271 - t423;
t148 = t204 * t270 + t421;
t139 = (-t225 * t275 + t227 * t273) * t242;
t137 = t185 * t271 + t187 * t270;
t136 = t185 * t270 - t187 * t271;
t135 = (qJD(4) + t242) * t225 + t298;
t84 = -qJ(5) * t404 + qJ(6) * t403;
t80 = -t275 * t405 + t361;
t75 = -t273 * t405 - t360;
t72 = -t135 * t270 - t463;
t69 = t135 * t271 - t464;
t59 = -qJ(6) * t402 - t134 * t380;
t57 = t271 * t80 + t439;
t54 = t270 * t80 - t438;
t52 = t372 - t466;
t48 = t374 - t442;
t46 = -pkin(3) * t75 - t301;
t43 = t51 - t468;
t42 = t50 - t443;
t39 = -pkin(3) * t74 - t389;
t38 = t45 - t440;
t37 = -pkin(4) * t150 + t44;
t36 = qJ(6) * t384 - t391 + t47;
t31 = (-t159 - t384) * qJ(6) + t40 + t391;
t30 = -t280 - t443;
t29 = t231 - t287 + t468;
t28 = -qJ(5) * t410 - t273 * t41 - t442;
t27 = pkin(4) * t358 + t275 * t40 + t466;
t26 = t227 * t312 + t284 + t424 + t440;
t25 = (-t396 - t384) * qJ(6) + (-t404 - t164) * pkin(5) + t41 + t394;
t23 = t270 * t33 - t271 * t85;
t22 = -t32 + t467;
t21 = -qJ(6) * t125 + t150 * t380 - t289 + t400;
t20 = t273 * t45 + t275 * t44;
t19 = t273 * t44 - t275 * t45;
t18 = -t388 + t400 + t468;
t17 = -t380 * t403 + t34 - t370 - t443;
t16 = -t273 * t59 + t275 * t31 + t466;
t15 = -t25 * t273 + t275 * t84 - t442;
t14 = -qJ(5) * t36 - qJ(6) * t34;
t13 = -pkin(8) * t75 - t273 * t37 + t275 * t38;
t12 = t273 * t34 + t275 * t35;
t11 = t273 * t35 - t275 * t34;
t10 = t20 * t271 + t270 * t47;
t9 = t20 * t270 - t271 * t47;
t8 = -pkin(8) * t19 + (pkin(4) * t273 - qJ(5) * t275) * t47;
t7 = -qJ(6) * t35 - t36 * t380;
t6 = -pkin(3) * t19 - t305;
t5 = t12 * t271 + t270 * t36;
t4 = t12 * t270 - t271 * t36;
t3 = -t21 * t273 + t26 * t275 - t467;
t2 = -pkin(3) * t11 - t390;
t1 = -pkin(8) * t11 + t14 * t275 - t273 * t7;
t24 = [0, 0, 0, 0, 0, qJDD(1), t306, t292, 0, 0, (t253 + t313) * t274, t252 * t276 + t254 * t274, t335 + t276 * (-t263 + t381), t254 * t276, t274 * (t265 - t381) + t333, 0, t276 * t250 + pkin(1) * t254 + pkin(7) * (t276 * t398 - t335), -t274 * t250 - pkin(1) * t252 + pkin(7) * (-t333 - t274 * (-t263 - t381)), pkin(1) * (t263 + t265) + (t267 + t268) * t369 + t308, pkin(1) * t250 + pkin(7) * t308, t274 * (t215 * t271 - t248 * t325) + t276 * (t215 * t270 + t248 * t324), t274 * (-t184 * t271 - t186 * t270) + t276 * (-t184 * t270 + t186 * t271), t274 * (-t233 * t270 + t421) + t276 * (t233 * t271 + t423), t274 * (t246 * t324 + t270 * t307) + t276 * (t246 * t325 - t271 * t307), t274 * (t232 * t271 - t343) + t276 * (t232 * t270 + t342), (t274 * (-t246 * t271 + t248 * t270) + t276 * (-t246 * t270 - t248 * t271)) * qJD(2), t274 * (-qJ(3) * t148 - t346) + t276 * (-pkin(2) * t184 + qJ(3) * t149 + t345) - pkin(1) * t184 + pkin(7) * (-t148 * t274 + t149 * t276), t274 * (-qJ(3) * t166 - t345) + t276 * (-pkin(2) * t186 + qJ(3) * t167 - t346) - pkin(1) * t186 + pkin(7) * (-t166 * t274 + t167 * t276), t274 * (-qJ(3) * t136 - t60) + t276 * (-pkin(2) * t181 + qJ(3) * t137 + t61) - pkin(1) * t181 + pkin(7) * (-t136 * t274 + t137 * t276), -qJ(3) * t373 + t276 * (pkin(2) * t180 + qJ(3) * t61) + pkin(1) * t180 + pkin(7) * (t276 * t61 - t373), t387, t472, -t482, t386, -t483, t385, t274 * (-t270 * t42 + t271 * t48 - t456) + t276 * (t270 * t48 + t271 * t42 + t452) + t449, t274 * (-qJ(3) * t69 - t270 * t43 + t271 * t52) + t276 * (qJ(3) * t72 + t270 * t52 + t271 * t43 - t469) - t470 + pkin(7) * (-t274 * t69 + t276 * t72), t274 * (t22 * t271 + t489) + t276 * (t270 * t22 - t488) + t490 - (t274 * t377 + t276 * t315 - pkin(1)) * t74, (t274 * (-pkin(8) * t271 + t377) + t276 * (-pkin(8) * t270 + t315) - pkin(1)) * t32 + (pkin(7) + qJ(3)) * (-t23 * t274 + (t270 * t85 + t271 * t33) * t276), t387, -t482, -t472, t385, t483, t386, t274 * (-t270 * t30 + t271 * t28 - t456) + t276 * (t270 * t28 + t271 * t30 + t452) + t449, t274 * (-qJ(3) * t54 + t13 * t271 - t270 * t46) + t276 * (-pkin(2) * t75 + qJ(3) * t57 + t13 * t270 + t271 * t46) - pkin(1) * t75 + pkin(7) * (-t274 * t54 + t276 * t57), t274 * (t27 * t271 - t270 * t29 - t480) + t276 * (t27 * t270 + t271 * t29 + t487) + t486, t274 * (-qJ(3) * t9 - t270 * t6 + t271 * t8) + t276 * (-pkin(2) * t19 + qJ(3) * t10 + t270 * t8 + t271 * t6) - pkin(1) * t19 + pkin(7) * (t10 * t276 - t274 * t9), t387, -t472, t482, t386, -t483, t274 * (t139 * t271 + t341) + t276 * (t139 * t270 - t202), t274 * (t15 * t271 - t17 * t270 - t456) + t276 * (t15 * t270 + t17 * t271 + t452) + t449, t274 * (t16 * t271 - t18 * t270 - t480) + t276 * (t16 * t270 + t18 * t271 + t487) + t486, t274 * (-t270 * t39 + t271 * t3 - t489) + t276 * (-pkin(2) * t74 + t270 * t3 + t271 * t39 + t488) - pkin(1) * t74 - t490, t274 * (-qJ(3) * t4 + t1 * t271 - t2 * t270) + t276 * (-pkin(2) * t11 + qJ(3) * t5 + t1 * t270 + t2 * t271) - pkin(1) * t11 + pkin(7) * (-t274 * t4 + t276 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t260, t263 - t265, t261, t260, t317, qJDD(2), -t228, -t229, 0, 0, t213, t245 - t244, t187, -t213, t185, qJDD(2), pkin(2) * t148 + t238 - t310, pkin(2) * t166 - t104, pkin(2) * t136, pkin(2) * t60, t116, -t428, -t455, t304, t431, t294, -t372 + t453, pkin(2) * t69 + pkin(3) * t135 + t374 - t465, t33 + t485, pkin(2) * t23 - pkin(3) * t85 + pkin(8) * t33, t116, -t455, t428, t294, -t431, t304, -qJ(5) * t411 + t275 * t41 + t453, pkin(2) * t54 + pkin(8) * t80 + t273 * t38 + t275 * t37 - t444, -pkin(4) * t357 + t273 * t40 + t477, pkin(2) * t9 + pkin(8) * t20 + (-pkin(4) * t275 - qJ(5) * t273 - pkin(3)) * t47, t116, t428, t455, t304, t431, t294, t25 * t275 + t273 * t84 + t453, t273 * t31 + t275 * t59 + t477, t21 * t275 + t26 * t273 - t485, pkin(2) * t4 - pkin(3) * t36 + pkin(8) * t12 + t14 * t273 + t275 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t186, t181, -t180, 0, 0, 0, 0, 0, 0, t434, t97, -t74, t32, 0, 0, 0, 0, 0, 0, t434, t75, -t97, t19, 0, 0, 0, 0, 0, 0, t434, -t97, t74, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t399, t397, -t183, -t405, t212, -t50, -t51, 0, 0, t183, t397, -t399, t212, t405, -t183, t280, t301, t230 + t287, t305, t183, -t399, -t397, -t183, -t405, t212, -t339 + (t403 - t183) * pkin(5) + t284 + t435, t388 + t401, t389, t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t403, t397, t159, t45, 0, 0, 0, 0, 0, 0, -t403, t159, -t397, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t404, -t134, t150, -t36;];
tauJ_reg  = t24;
