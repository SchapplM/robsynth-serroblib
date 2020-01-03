% Calculate vector of inverse dynamics joint torques for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:47
% DurationCPUTime: 13.43s
% Computational Cost: add. (5987->425), mult. (7074->525), div. (0->0), fcn. (5487->6), ass. (0->236)
t448 = Icges(4,3) + Icges(5,3);
t205 = qJ(1) + pkin(6);
t198 = sin(t205);
t209 = cos(qJ(3));
t310 = t198 * t209;
t207 = sin(qJ(3));
t311 = t198 * t207;
t199 = cos(t205);
t320 = Icges(5,6) * t199;
t91 = Icges(5,4) * t310 - Icges(5,2) * t311 - t320;
t321 = Icges(4,6) * t199;
t93 = Icges(4,4) * t310 - Icges(4,2) * t311 - t321;
t442 = t91 + t93;
t150 = Icges(5,5) * t209 - Icges(5,6) * t207;
t152 = Icges(4,5) * t209 - Icges(4,6) * t207;
t430 = t150 + t152;
t171 = Icges(5,4) * t311;
t324 = Icges(5,5) * t199;
t95 = Icges(5,1) * t310 - t171 - t324;
t172 = Icges(4,4) * t311;
t325 = Icges(4,5) * t199;
t97 = Icges(4,1) * t310 - t172 - t325;
t446 = t95 + t97;
t326 = Icges(5,4) * t207;
t158 = Icges(5,1) * t209 - t326;
t96 = Icges(5,5) * t198 + t158 * t199;
t327 = Icges(4,4) * t207;
t160 = Icges(4,1) * t209 - t327;
t98 = Icges(4,5) * t198 + t160 * t199;
t434 = t96 + t98;
t447 = t448 * t199;
t412 = t447 + (Icges(4,6) + Icges(5,6)) * t311 + (-Icges(4,5) - Icges(5,5)) * t310;
t420 = t448 * t198 + t430 * t199;
t153 = Icges(5,2) * t209 + t326;
t155 = Icges(4,2) * t209 + t327;
t445 = t153 + t155;
t200 = Icges(5,4) * t209;
t157 = Icges(5,1) * t207 + t200;
t201 = Icges(4,4) * t209;
t159 = Icges(4,1) * t207 + t201;
t440 = t157 + t159;
t444 = t442 * t207;
t443 = t434 * t310;
t239 = -Icges(5,2) * t207 + t200;
t92 = Icges(5,6) * t198 + t199 * t239;
t240 = -Icges(4,2) * t207 + t201;
t94 = Icges(4,6) * t198 + t199 * t240;
t435 = t92 + t94;
t439 = t158 + t160;
t438 = t239 + t240;
t417 = -t446 * t209 + t444;
t427 = t445 * t207 - t440 * t209;
t437 = t420 * t199 - t443;
t308 = t199 * t209;
t436 = t412 * t198 - t446 * t308;
t397 = t420 * t198 + t434 * t308;
t391 = -t417 * t198 + t412 * t199;
t390 = -t435 * t311 - t437;
t309 = t199 * t207;
t389 = -t442 * t309 - t436;
t388 = -t435 * t309 + t397;
t433 = t438 * qJD(3);
t432 = t439 * qJD(3);
t149 = Icges(5,5) * t207 + Icges(5,6) * t209;
t151 = Icges(4,5) * t207 + Icges(4,6) * t209;
t431 = t149 + t151;
t429 = t445 * qJD(3);
t428 = t440 * qJD(3);
t426 = -t440 * t207 - t209 * t445;
t312 = t151 * t199;
t314 = t149 * t199;
t401 = -t427 * t198 - t312 - t314;
t313 = t151 * t198;
t315 = t149 * t198;
t400 = -t427 * t199 + t313 + t315;
t425 = t435 * t207;
t424 = -t429 * t199 + (-t438 * t198 + t320 + t321) * qJD(1);
t423 = qJD(1) * t435 - t198 * t429;
t422 = -t428 * t199 + (-t439 * t198 + t324 + t325) * qJD(1);
t421 = -qJD(1) * t434 + t198 * t428;
t419 = qJD(1) * t431 + qJD(3) * t426 - t207 * t433 + t209 * t432;
t418 = t209 * t434 - t425;
t416 = t388 * t198 - t199 * t389;
t415 = t390 * t198 - t391 * t199;
t414 = t427 * qJD(1) + qJD(3) * t430;
t210 = cos(qJ(1));
t204 = t210 * pkin(1);
t387 = t446 * t207 + t442 * t209;
t386 = t207 * t434 + t209 * t435;
t413 = t400 * qJD(1);
t129 = rSges(3,1) * t198 + rSges(3,2) * t199;
t208 = sin(qJ(1));
t365 = t208 * pkin(1);
t122 = -t129 - t365;
t411 = t431 * qJD(3);
t410 = t401 * qJD(1);
t132 = t199 * pkin(2) + t198 * pkin(5);
t124 = t132 * qJD(1);
t287 = qJD(1) * qJD(3);
t126 = -qJDD(3) * t199 + t198 * t287;
t202 = t209 * rSges(5,1);
t381 = -rSges(5,2) * t207 + t202;
t139 = t381 * qJD(3);
t342 = t209 * rSges(5,2);
t161 = rSges(5,1) * t207 + t342;
t182 = qJD(4) * t199;
t195 = t199 * pkin(5);
t131 = pkin(2) * t198 - t195;
t206 = -qJ(4) - pkin(5);
t177 = t199 * t206;
t203 = t209 * pkin(3);
t362 = t203 + pkin(2);
t385 = rSges(5,1) * t310 - rSges(5,2) * t311 - t199 * rSges(5,3) + t198 * t362 + t177;
t355 = -t131 + t385;
t264 = -t355 - t365;
t243 = -t131 + t264;
t212 = qJD(1) ^ 2;
t285 = t212 * t204;
t289 = qJD(3) * t199;
t307 = t209 * qJD(3) ^ 2;
t117 = t161 * t198;
t189 = t198 * rSges(5,3);
t286 = pkin(3) * t311;
t300 = qJD(3) * t286 + t182;
t383 = t381 + t203;
t360 = -qJD(3) * t117 - t300 + ((-pkin(5) - t206) * t198 + t189 + t383 * t199) * qJD(1);
t10 = -t285 - t139 * t289 + qJDD(4) * t198 + t126 * t161 + (t126 * t207 - t199 * t307) * pkin(3) + t243 * qJDD(1) + (-t124 + t182 - t360) * qJD(1);
t409 = -g(1) + t10;
t125 = qJDD(3) * t198 + t199 * t287;
t181 = qJD(4) * t198;
t291 = qJD(1) * t199;
t180 = pkin(5) * t291;
t263 = qJDD(1) * t204 - t212 * t365;
t292 = qJD(1) * t198;
t232 = qJD(1) * (-pkin(2) * t292 + t180) + qJDD(1) * t132 + t263;
t290 = qJD(3) * t198;
t384 = rSges(5,1) * t308 - rSges(5,2) * t309 - t198 * t206 + t199 * t362 + t189;
t338 = -t132 + t384;
t282 = t207 * t289;
t225 = -t209 * t292 - t282;
t288 = qJD(3) * t209;
t281 = t199 * t288;
t284 = t207 * t292;
t382 = rSges(5,2) * t284 + rSges(5,3) * t291 + t181;
t361 = -pkin(3) * t282 - t180 + (-t198 * t203 - t177) * qJD(1) + rSges(5,1) * t225 - rSges(5,2) * t281 + t382;
t11 = -t139 * t290 - qJDD(4) * t199 - t125 * t161 + t338 * qJDD(1) + (-t125 * t207 - t198 * t307) * pkin(3) + (t181 + t361) * qJD(1) + t232;
t408 = g(2) - t11;
t407 = qJD(3) * t415 + t410;
t406 = qJD(3) * t416 + t413;
t405 = qJD(3) * t417 + t207 * t421 - t209 * t423;
t404 = t418 * qJD(3) + t422 * t207 + t209 * t424;
t403 = t198 * t414 + t199 * t419;
t402 = t198 * t419 - t199 * t414;
t399 = t412 + t425;
t398 = t420 * qJD(1);
t396 = t199 ^ 2;
t395 = -t386 * qJD(3) - t207 * t424 + t422 * t209 + t398;
t394 = t412 * qJD(1) + t387 * qJD(3) + t207 * t423 + t209 * t421;
t393 = qJD(1) * t417 - t411 * t198 + t398;
t392 = -t411 * t199 + (-t198 * t430 - t418 + t447) * qJD(1);
t190 = t198 * rSges(4,3);
t102 = rSges(4,1) * t308 - rSges(4,2) * t309 + t190;
t272 = t132 + t204;
t74 = t102 + t272;
t130 = t199 * rSges(3,1) - rSges(3,2) * t198;
t123 = t130 + t204;
t380 = t393 * t396 + (t395 * t198 + (-t392 + t394) * t199) * t198;
t379 = t394 * t396 + (t392 * t198 + (-t393 + t395) * t199) * t198;
t331 = -t159 * t198 - t93;
t335 = -Icges(4,2) * t310 - t172 + t97;
t371 = -t207 * t335 + t209 * t331;
t333 = -t157 * t198 - t91;
t337 = -Icges(5,2) * t310 - t171 + t95;
t370 = -t207 * t337 + t209 * t333;
t369 = t125 / 0.2e1;
t368 = t126 / 0.2e1;
t354 = rSges(4,1) * t209;
t162 = rSges(4,1) * t207 + rSges(4,2) * t209;
t120 = t162 * t199;
t42 = qJD(1) * t74 - t162 * t290;
t351 = t120 * t42;
t295 = rSges(4,2) * t311 + t199 * rSges(4,3);
t100 = rSges(4,1) * t310 - t295;
t273 = -t100 - t365;
t261 = -t131 + t273;
t283 = t162 * t289;
t41 = qJD(1) * t261 - t283;
t350 = t198 * t41;
t349 = t199 * t41;
t271 = -pkin(3) * t207 - t161;
t224 = t271 * t289 + t181;
t30 = qJD(1) * t243 + t224;
t340 = t30 * t161;
t336 = -t153 * t199 + t96;
t334 = -t155 * t199 + t98;
t332 = -t157 * t199 - t92;
t330 = -t159 * t199 - t94;
t301 = rSges(4,2) * t284 + rSges(4,3) * t291;
t299 = -t153 + t158;
t298 = t157 + t239;
t297 = -t155 + t160;
t296 = t159 + t240;
t294 = qJD(1) * t150;
t293 = qJD(1) * t152;
t280 = -pkin(2) - t354;
t277 = -t290 / 0.2e1;
t276 = t290 / 0.2e1;
t275 = -t289 / 0.2e1;
t274 = t289 / 0.2e1;
t262 = -pkin(3) * t309 - t161 * t199;
t260 = -pkin(3) * t288 - t139;
t166 = rSges(2,1) * t210 - rSges(2,2) * t208;
t163 = rSges(2,1) * t208 + rSges(2,2) * t210;
t165 = -rSges(4,2) * t207 + t354;
t249 = -t198 * t42 - t349;
t68 = rSges(4,1) * t225 - rSges(4,2) * t281 + t301;
t118 = t162 * t198;
t70 = -qJD(3) * t118 + (t165 * t199 + t190) * qJD(1);
t248 = t198 * t70 + t199 * t68;
t238 = t100 * t198 + t102 * t199;
t233 = -t342 + (-rSges(5,1) - pkin(3)) * t207;
t223 = t199 * t233;
t222 = -t207 * t336 + t209 * t332;
t221 = -t207 * t334 + t209 * t330;
t220 = (-t207 * t298 + t209 * t299) * qJD(1);
t219 = (-t207 * t296 + t209 * t297) * qJD(1);
t140 = t165 * qJD(3);
t128 = qJD(1) * t131;
t40 = qJD(3) * t238 + qJD(2);
t31 = -t161 * t290 + (t272 + t338) * qJD(1) - t300;
t27 = qJD(2) + (t198 * t355 + t199 * t338) * qJD(3);
t22 = qJD(1) * t68 + qJDD(1) * t102 - t125 * t162 - t140 * t290 + t232;
t21 = -t285 - t140 * t289 + t126 * t162 + (-t124 - t70) * qJD(1) + t261 * qJDD(1);
t16 = qJD(3) * t248 + t100 * t125 - t102 * t126 + qJDD(2);
t1 = qJDD(2) - t338 * t126 + t355 * t125 + (t198 * t360 + t199 * t361) * qJD(3);
t2 = [-m(2) * (-g(1) * t163 + g(2) * t166) + ((t397 * t198 + ((t420 + t444) * t199 + t390 + t436 - t443) * t199) * qJD(3) + t413) * t274 + (-t427 * qJD(3) + t432 * t207 + t433 * t209) * qJD(1) + ((-t129 * t212 - g(2) + t263) * t123 + (-g(1) - t285 + (-0.2e1 * t130 - t204 + t123) * t212) * t122) * m(3) + (-(qJD(1) * t264 - t128 + t224 - t30) * t31 + t30 * t300 + t31 * t382 + (t198 * t340 + t223 * t31) * qJD(3) + ((-t208 * t31 - t210 * t30) * pkin(1) + (t30 * (-t381 - t362) - t31 * t206) * t199 + (t30 * (t206 - rSges(5,3)) + t31 * (-t202 - t362)) * t198) * qJD(1) - t408 * (t204 + t384) + t409 * (-t365 - t385)) * m(5) + (t42 * (t180 + t301) + (t162 * t350 - t351) * qJD(3) + ((-t208 * t42 - t210 * t41) * pkin(1) + (-pkin(2) - t165) * t349 + (t41 * (-rSges(4,3) - pkin(5)) + t42 * t280) * t198) * qJD(1) - (qJD(1) * t273 - t128 - t283 - t41) * t42 + (t22 - g(2)) * t74 + (t21 - g(1)) * (t280 * t198 + t195 + t295 - t365)) * m(4) + (t400 + t386) * t369 + (t401 + t387) * t368 + (((t199 * t399 + t388 - t397) * t199 + (t198 * t399 + t389 + t437) * t198) * qJD(3) + t407 - t410) * t277 + (t403 + t404) * t276 + (m(3) * (t122 ^ 2 + t130 * t123) + m(2) * (t163 ^ 2 + t166 ^ 2) + Icges(2,3) + Icges(3,3) - t426) * qJDD(1) + (t402 - t405 + t406) * t275; m(3) * qJDD(2) + m(4) * t16 + m(5) * t1 + (-m(3) - m(4) - m(5)) * g(3); t416 * t369 + t415 * t368 + (qJD(1) * t403 + t379 * qJD(3) + qJDD(1) * t400 + t388 * t125 + t389 * t126) * t198 / 0.2e1 - (qJD(1) * t402 + t380 * qJD(3) + qJDD(1) * t401 + t390 * t125 + t391 * t126) * t199 / 0.2e1 - ((((-t335 - t337) * t199 + (t334 + t336) * t198) * t209 + ((-t331 - t333) * t199 + (t330 + t332) * t198) * t207) * qJD(3) + ((t296 + t298) * t209 + (t297 + t299) * t207) * qJD(1)) * qJD(1) / 0.2e1 + (t405 * t199 + t404 * t198 + (t387 * t198 + t386 * t199) * qJD(1)) * qJD(1) / 0.2e1 + (t198 * t386 - t199 * t387) * qJDD(1) / 0.2e1 + t407 * t292 / 0.2e1 + t406 * t291 / 0.2e1 + ((-t290 * t312 + t293) * t198 + (t219 + (-t371 * t199 + (t313 + t221) * t198) * qJD(3)) * t199 + (-t290 * t314 + t294) * t198 + (t220 + (-t370 * t199 + (t315 + t222) * t198) * qJD(3)) * t199) * t277 + ((t389 * t198 + t199 * t388) * qJD(1) + t379) * t276 + ((t198 * t391 + t390 * t199) * qJD(1) + t380) * t275 + ((-t289 * t313 - t293) * t199 + (t219 + (t221 * t198 + (t312 - t371) * t199) * qJD(3)) * t198 + (-t289 * t315 - t294) * t199 + (t220 + (t222 * t198 + (t314 - t370) * t199) * qJD(3)) * t198) * t274 + (g(1) * t120 + g(2) * t118 - g(3) * t165 + t16 * t238 + t40 * ((t100 * t199 - t102 * t198) * qJD(1) + t248) + t249 * t140 + (-t22 * t198 - t21 * t199 + (-t199 * t42 + t350) * qJD(1)) * t162 - (t118 * t41 - t351) * qJD(1) - (t40 * (-t118 * t198 - t120 * t199) + t249 * t165) * qJD(3)) * m(4) + (-(t30 * t117 + t262 * t31) * qJD(1) - g(3) * t383 - g(1) * t223 + (-(t262 * t27 - t30 * t383) * qJD(3) + t10 * t271 + t30 * t260 + t1 * t338 + t27 * t361 + (t27 * t355 + t271 * t31) * qJD(1)) * t199 + (-(-t31 * t383 + (-t117 - t286) * t27) * qJD(3) + t11 * t271 + t31 * t260 + t1 * t355 + t27 * t360 + (-t27 * t338 + t340) * qJD(1) - g(2) * t233) * t198) * m(5); (t198 * t409 + t408 * t199) * m(5);];
tau = t2;
