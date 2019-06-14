% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:38:27
% EndTime: 2019-05-05 21:38:49
% DurationCPUTime: 9.39s
% Computational Cost: add. (16936->464), mult. (40262->558), div. (0->0), fcn. (29874->8), ass. (0->279)
t258 = sin(pkin(9));
t259 = cos(pkin(9));
t262 = sin(qJ(3));
t264 = cos(qJ(3));
t324 = t259 * t262;
t286 = t258 * t264 + t324;
t245 = t286 * qJD(1);
t261 = sin(qJ(4));
t263 = cos(qJ(4));
t225 = qJD(3) * t261 + t245 * t263;
t222 = t225 ^ 2;
t223 = -t263 * qJD(3) + t245 * t261;
t370 = t223 ^ 2;
t150 = -t370 - t222;
t427 = t150 * t264;
t312 = qJD(1) * t259;
t325 = t258 * t262;
t243 = qJD(1) * t325 - t264 * t312;
t237 = qJD(4) + t243;
t242 = t286 * qJDD(1);
t308 = t243 * qJD(3);
t211 = t242 - t308;
t298 = -t263 * qJDD(3) + t211 * t261;
t122 = (qJD(4) - t237) * t225 + t298;
t285 = -qJDD(3) * t261 - t211 * t263;
t278 = qJD(4) * t223 + t285;
t329 = t223 * t237;
t384 = t329 - t278;
t347 = t384 * t261;
t76 = t122 * t263 - t347;
t50 = t262 * t76 + t427;
t426 = t262 * t150;
t53 = t264 * t76 - t426;
t474 = qJ(2) * (t258 * t50 - t259 * t53);
t473 = pkin(7) * t50;
t472 = pkin(7) * t53;
t131 = t329 + t278;
t369 = t237 ^ 2;
t164 = -t369 - t222;
t179 = t225 * t223;
t304 = t259 * qJDD(1);
t305 = t258 * qJDD(1);
t288 = t262 * t305 - t264 * t304;
t307 = t245 * qJD(3);
t209 = -t288 - t307;
t200 = qJDD(4) - t209;
t388 = t179 + t200;
t405 = t263 * t388;
t441 = t164 * t261 + t405;
t450 = t264 * t441;
t444 = t262 * t131 + t450;
t402 = t388 * t261;
t92 = t164 * t263 - t402;
t456 = pkin(2) * t92;
t471 = pkin(7) * t444 + t456;
t451 = t262 * t441;
t445 = -t264 * t131 + t451;
t457 = pkin(1) * t92;
t470 = t457 + qJ(2) * (-t258 * t445 + t259 * t444);
t385 = t222 - t370;
t345 = t131 * t261;
t161 = qJD(4) * t225 + t298;
t328 = t237 * t225;
t390 = t161 + t328;
t398 = t390 * t263;
t79 = t345 - t398;
t459 = t259 * (t262 * t79 - t264 * t385) + t258 * (t262 * t385 + t264 * t79);
t186 = t222 - t369;
t389 = -t179 + t200;
t400 = t389 * t263;
t424 = t186 * t261 + t400;
t435 = t258 * (t262 * t384 + t264 * t424) + t259 * (t262 * t424 - t264 * t384);
t382 = t370 - t369;
t105 = -t263 * t382 + t402;
t391 = t161 - t328;
t468 = t258 * (t105 * t264 + t262 * t391) + t259 * (t262 * t105 - t391 * t264);
t455 = pkin(3) * t92;
t453 = pkin(8) * t92;
t466 = pkin(7) * t445;
t346 = t384 * t263;
t72 = t122 * t261 + t346;
t454 = pkin(8) * t72;
t383 = -t369 - t370;
t401 = t389 * t261;
t421 = t263 * t383 - t401;
t439 = t262 * t421 - t264 * t390;
t443 = pkin(7) * t439;
t452 = pkin(8) * t441;
t432 = pkin(3) * t150;
t449 = -pkin(8) * t76 - t432;
t422 = t261 * t383 + t400;
t438 = t262 * t390 + t264 * t421;
t440 = -pkin(2) * t422 + pkin(7) * t438;
t448 = -pkin(3) * t131 + t452;
t437 = qJ(2) * (-t258 * t439 + t259 * t438) - pkin(1) * t422;
t431 = pkin(3) * t422;
t430 = pkin(8) * t422;
t442 = t186 * t263 - t401;
t344 = t131 * t263;
t399 = t390 * t261;
t416 = t344 + t399;
t419 = t261 * t382 + t405;
t418 = -pkin(3) * t390 + pkin(8) * t421;
t429 = qJ(5) * t150;
t428 = t131 * qJ(5);
t354 = qJ(5) * t383;
t423 = pkin(4) * t389 + t354;
t413 = qJ(5) * t388;
t412 = qJ(6) * t384;
t213 = t245 * t243;
t380 = -t213 + qJDD(3);
t409 = t262 * t380;
t404 = t264 * t380;
t172 = t262 * t179;
t173 = t264 * t179;
t327 = t237 * t261;
t182 = t225 * t327;
t314 = -t263 * t278 - t182;
t373 = t259 * (t262 * t314 - t173) + t258 * (t264 * t314 + t172);
t362 = sin(qJ(1));
t363 = cos(qJ(1));
t295 = g(1) * t362 - t363 * g(2);
t284 = -qJDD(2) + t295;
t367 = qJD(1) ^ 2;
t311 = t367 * qJ(2);
t353 = qJDD(1) * pkin(1);
t239 = t284 + t311 + t353;
t255 = t258 ^ 2;
t256 = t259 ^ 2;
t313 = t255 + t256;
t395 = t311 * t313 - t239 - t353;
t302 = pkin(2) * t259 + pkin(1);
t204 = t302 * qJDD(1) + (pkin(7) * t313 + qJ(2)) * t367 + t284;
t394 = pkin(7) + qJ(2);
t310 = qJD(5) * t237;
t226 = 0.2e1 * t310;
t309 = qJD(6) * t223;
t387 = 0.2e1 * t309 + t226;
t227 = -0.2e1 * t310;
t386 = -0.2e1 * t309 + t227;
t289 = -pkin(5) * t237 - qJ(6) * t225;
t381 = t225 * t289 + qJDD(6);
t379 = -pkin(5) * t161 + t381;
t101 = (-t211 + t308) * pkin(8) + (-t209 + t307) * pkin(3) - t204;
t201 = pkin(3) * t243 - pkin(8) * t245;
t280 = g(1) * t363 + g(2) * t362;
t268 = (-t394 * qJDD(1) + (qJD(1) * t302 - (2 * qJD(2))) * qJD(1) + t280) * t258;
t360 = t259 * g(3);
t267 = t268 - t360;
t276 = qJDD(1) * qJ(2) - t280;
t365 = 2 * qJD(2);
t296 = -g(3) * t258 + t259 * (-pkin(1) * t367 + t276) + t312 * t365;
t190 = -pkin(2) * t256 * t367 + pkin(7) * t304 + t296;
t317 = t264 * t190;
t366 = qJD(3) ^ 2;
t98 = -pkin(3) * t366 + qJDD(3) * pkin(8) - t243 * t201 + t262 * t267 + t317;
t58 = -t263 * t101 + t261 * t98;
t281 = -t200 * pkin(4) - qJ(5) * t369 + qJDD(5) + t58;
t275 = -t200 * pkin(5) + t281 - t412;
t171 = pkin(4) * t223 - qJ(5) * t225;
t300 = -pkin(5) * t223 - t171;
t32 = (-0.2e1 * qJD(6) - t300) * t225 + t275;
t59 = t261 * t101 + t263 * t98;
t291 = -pkin(4) * t369 + t200 * qJ(5) - t223 * t171 + t59;
t279 = -pkin(5) * t370 + t161 * qJ(6) + t237 * t289 + t291;
t33 = t279 + t387;
t364 = pkin(4) + pkin(5);
t377 = qJ(5) * t33 - t32 * t364;
t376 = qJ(5) * t122 + t364 * t384;
t375 = -t164 * t364 + t279 + t413;
t326 = t237 * t263;
t303 = t223 * t326;
t283 = t161 * t261 + t303;
t372 = t258 * (t264 * t283 - t172) + t259 * (t262 * t283 + t173);
t196 = t264 * t200;
t292 = t182 - t303;
t319 = t262 * t200;
t371 = t258 * (t264 * t292 + t319) + t259 * (t262 * t292 - t196);
t240 = t243 ^ 2;
t241 = t245 ^ 2;
t368 = 0.2e1 * t225;
t134 = t190 * t262 - t264 * t267;
t135 = -g(3) * t324 + t262 * t268 + t317;
t82 = -t134 * t264 + t262 * t135;
t358 = t258 * t82;
t97 = -qJDD(3) * pkin(3) - t366 * pkin(8) + t201 * t245 + t134;
t357 = t261 * t97;
t356 = t263 * t97;
t334 = t171 * t225;
t332 = t204 * t262;
t331 = t204 * t264;
t205 = qJDD(3) + t213;
t330 = t205 * t264;
t318 = t262 * t205;
t301 = -pkin(3) * t264 - pkin(2);
t35 = t261 * t58 + t263 * t59;
t83 = t134 * t262 + t264 * t135;
t297 = t258 * (t360 + ((-pkin(1) * qJD(1) + t365) * qJD(1) + t276) * t258) + t259 * t296;
t45 = t226 + t291;
t46 = t281 + t334;
t294 = -pkin(4) * t46 + qJ(5) * t45;
t113 = t225 * t326 - t261 * t278;
t293 = -t263 * t161 + t223 * t327;
t290 = -pkin(4) * t384 - qJ(5) * t391;
t34 = t261 * t59 - t263 * t58;
t282 = (-t223 * t261 - t225 * t263) * t237;
t277 = -pkin(4) * t164 + t291 + t413;
t273 = qJD(6) * t368 - t275;
t272 = t161 * pkin(4) + t428 + t97;
t270 = qJD(5) * t368 - t272;
t269 = -t46 + t423;
t47 = (pkin(4) * t237 - 0.2e1 * qJD(5)) * t225 + t272;
t40 = (-t390 - t328) * pkin(4) + t270;
t39 = -pkin(4) * t328 + t270 - t428;
t251 = t256 * qJDD(1);
t250 = t255 * qJDD(1);
t246 = t313 * t367;
t230 = -t241 - t366;
t229 = -t241 + t366;
t228 = t240 - t366;
t210 = t242 - 0.2e1 * t308;
t208 = t288 + 0.2e1 * t307;
t202 = -t366 - t240;
t177 = -t240 - t241;
t167 = -t262 * t230 - t330;
t166 = t230 * t264 - t318;
t149 = t262 * t242 - t264 * t288;
t148 = -t242 * t264 - t262 * t288;
t147 = t202 * t264 - t409;
t146 = t262 * t202 + t404;
t137 = (-t223 * t263 + t225 * t261) * t237;
t132 = (qJD(4) + t237) * t223 + t285;
t81 = -qJ(5) * t390 + qJ(6) * t389;
t78 = -t263 * t391 + t347;
t73 = -t261 * t391 - t346;
t70 = -t262 * t132 - t450;
t67 = t132 * t264 - t451;
t56 = t356 - t453;
t54 = t264 * t78 + t426;
t51 = t262 * t78 - t427;
t49 = -qJ(6) * t388 - t131 * t364;
t48 = t357 - t430;
t44 = t59 - t455;
t43 = -pkin(3) * t73 - t290;
t42 = t58 - t431;
t41 = t46 - t429;
t38 = -pkin(4) * t150 + t45;
t37 = -pkin(3) * t72 - t376;
t36 = qJ(6) * t370 - t379 + t47;
t31 = -t269 - t431;
t30 = t39 + (-t164 - t370) * qJ(6) + t379;
t27 = t227 - t277 + t455;
t26 = -qJ(5) * t398 - t261 * t40 - t430;
t25 = pkin(4) * t345 + t263 * t39 + t453;
t24 = -t34 + t454;
t23 = t225 * t300 + t273 + t412 + t429;
t22 = t40 + (-t383 - t370) * qJ(6) + (-t390 - t161) * pkin(5) + t381;
t21 = t261 * t46 + t263 * t45;
t20 = t261 * t45 - t263 * t46;
t19 = -qJ(6) * t122 + t150 * t364 - t279 + t386;
t18 = -t375 + t386 + t455;
t17 = -t364 * t389 + t32 - t354 - t431;
t16 = -t261 * t49 + t263 * t30 + t453;
t15 = -t22 * t261 + t263 * t81 - t430;
t14 = -qJ(5) * t36 - qJ(6) * t32;
t13 = -pkin(8) * t73 - t261 * t38 + t263 * t41;
t12 = t21 * t264 + t262 * t47;
t11 = t262 * t21 - t264 * t47;
t10 = t261 * t32 + t263 * t33;
t9 = t261 * t33 - t263 * t32;
t8 = -pkin(8) * t20 + (pkin(4) * t261 - qJ(5) * t263) * t47;
t7 = -pkin(3) * t20 - t294;
t6 = -qJ(6) * t33 - t36 * t364;
t5 = t10 * t264 + t262 * t36;
t4 = t262 * t10 - t264 * t36;
t3 = -t19 * t261 + t23 * t263 - t454;
t2 = -pkin(3) * t9 - t377;
t1 = -pkin(8) * t9 + t14 * t263 - t261 * t6;
t28 = [0, 0, 0, 0, 0, qJDD(1), t295, t280, 0, 0, t250, 0.2e1 * t258 * t304, 0, t251, 0, 0, -t395 * t259, t395 * t258, pkin(1) * t246 + qJ(2) * (t251 + t250) + t297, pkin(1) * t239 + qJ(2) * t297, t258 * (t211 * t264 - t262 * t307) + t259 * (t262 * t211 + t264 * t307), t258 * (-t208 * t264 - t262 * t210) + t259 * (-t262 * t208 + t210 * t264), t258 * (-t262 * t229 + t404) + t259 * (t229 * t264 + t409), t258 * (-t262 * t209 + t264 * t308) + t259 * (t209 * t264 + t262 * t308), t258 * (t228 * t264 - t318) + t259 * (t262 * t228 + t330), (t258 * (-t243 * t264 + t245 * t262) + t259 * (-t243 * t262 - t245 * t264)) * qJD(3), t258 * (-pkin(7) * t146 - t332) + t259 * (-pkin(2) * t208 + pkin(7) * t147 + t331) - pkin(1) * t208 + qJ(2) * (-t146 * t258 + t147 * t259), t258 * (-pkin(7) * t166 - t331) + t259 * (-pkin(2) * t210 + pkin(7) * t167 - t332) - pkin(1) * t210 + qJ(2) * (-t166 * t258 + t167 * t259), t258 * (-pkin(7) * t148 - t82) + t259 * (-pkin(2) * t177 + pkin(7) * t149 + t83) - pkin(1) * t177 + qJ(2) * (-t148 * t258 + t149 * t259), -pkin(7) * t358 + t259 * (pkin(2) * t204 + pkin(7) * t83) + pkin(1) * t204 + qJ(2) * (t259 * t83 - t358), t373, t459, t435, t372, -t468, t371, t258 * (-t262 * t42 + t264 * t48 - t443) + t259 * (t262 * t48 + t264 * t42 + t440) + t437, t258 * (-pkin(7) * t67 - t262 * t44 + t264 * t56) + t259 * (pkin(7) * t70 + t262 * t56 + t264 * t44 - t456) - t457 + qJ(2) * (-t258 * t67 + t259 * t70), t258 * (t24 * t264 + t473) + t259 * (t262 * t24 - t472) + t474 - (pkin(3) * t325 + t259 * t301 - pkin(1)) * t72, (t258 * (pkin(3) * t262 - pkin(8) * t264) + t259 * (-pkin(8) * t262 + t301) - pkin(1)) * t34 + t394 * (-t258 * (t262 * t35 - t264 * t97) + t259 * (t262 * t97 + t264 * t35)), t373, t435, -t459, t371, t468, t372, t258 * (t26 * t264 - t262 * t31 - t443) + t259 * (t262 * t26 + t264 * t31 + t440) + t437, t258 * (-pkin(7) * t51 + t13 * t264 - t262 * t43) + t259 * (-pkin(2) * t73 + pkin(7) * t54 + t262 * t13 + t264 * t43) - pkin(1) * t73 + qJ(2) * (-t258 * t51 + t259 * t54), t258 * (t25 * t264 - t262 * t27 - t466) + t259 * (t262 * t25 + t264 * t27 + t471) + t470, t258 * (-pkin(7) * t11 - t262 * t7 + t264 * t8) + t259 * (-pkin(2) * t20 + pkin(7) * t12 + t262 * t8 + t264 * t7) - pkin(1) * t20 + qJ(2) * (-t11 * t258 + t12 * t259), t373, -t459, -t435, t372, -t468, t258 * (t137 * t264 + t319) + t259 * (t137 * t262 - t196), t258 * (t15 * t264 - t262 * t17 - t443) + t259 * (t262 * t15 + t17 * t264 + t440) + t437, t258 * (t16 * t264 - t262 * t18 - t466) + t259 * (t262 * t16 + t18 * t264 + t471) + t470, t258 * (-t262 * t37 + t264 * t3 - t473) + t259 * (-pkin(2) * t72 + t262 * t3 + t264 * t37 + t472) - pkin(1) * t72 - t474, t258 * (-pkin(7) * t4 + t1 * t264 - t262 * t2) + t259 * (-pkin(2) * t9 + pkin(7) * t5 + t262 * t1 + t2 * t264) - pkin(1) * t9 + qJ(2) * (-t258 * t4 + t259 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t304, t305, -t246, -t239, 0, 0, 0, 0, 0, 0, t208, t210, t177, -t204, 0, 0, 0, 0, 0, 0, t422, t92, -t72, t34, 0, 0, 0, 0, 0, 0, t422, t73, -t92, t20, 0, 0, 0, 0, 0, 0, t422, -t92, t72, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t241 - t240, t242, -t213, -t288, qJDD(3), -t134, -t135, 0, 0, t113, -t416, -t442, t293, t419, t282, -t356 + t418, pkin(3) * t132 + t357 - t452, t35 + t449, -pkin(3) * t97 + pkin(8) * t35, t113, -t442, t416, t282, -t419, t293, -qJ(5) * t399 + t263 * t40 + t418, pkin(8) * t78 + t261 * t41 + t263 * t38 - t432, -pkin(4) * t344 + t261 * t39 + t448, pkin(8) * t21 + (-pkin(4) * t263 - qJ(5) * t261 - pkin(3)) * t47, t113, t416, t442, t293, t419, t282, t22 * t263 + t261 * t81 + t418, t261 * t30 + t263 * t49 + t448, t19 * t263 + t23 * t261 - t449, -pkin(3) * t36 + pkin(8) * t10 + t14 * t261 + t263 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t385, t384, -t179, -t391, t200, -t58, -t59, 0, 0, t179, t384, -t385, t200, t391, -t179, t269, t290, t226 + t277, t294, t179, -t385, -t384, -t179, -t391, t200, -t334 + (t389 - t179) * pkin(5) + t273 + t423, t375 + t387, t376, t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t389, t384, t164, t46, 0, 0, 0, 0, 0, 0, -t389, t164, -t384, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t390, -t131, t150, -t36;];
tauJ_reg  = t28;
