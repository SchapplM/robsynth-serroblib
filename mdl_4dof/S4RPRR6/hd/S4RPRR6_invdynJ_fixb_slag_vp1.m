% Calculate vector of inverse dynamics joint torques for
% S4RPRR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:39
% DurationCPUTime: 7.84s
% Computational Cost: add. (9882->581), mult. (10181->761), div. (0->0), fcn. (8007->8), ass. (0->323)
t259 = pkin(7) + qJ(3);
t235 = sin(t259);
t457 = rSges(4,2) * t235;
t239 = qJ(4) + t259;
t225 = cos(t239);
t264 = sin(qJ(1));
t393 = t225 * t264;
t224 = sin(t239);
t395 = t224 * t264;
t265 = cos(qJ(1));
t405 = Icges(5,6) * t265;
t122 = Icges(5,4) * t393 - Icges(5,2) * t395 - t405;
t216 = Icges(5,4) * t225;
t173 = Icges(5,1) * t224 + t216;
t456 = -t173 * t264 - t122;
t301 = -Icges(5,2) * t224 + t216;
t123 = Icges(5,6) * t264 + t265 * t301;
t455 = -t173 * t265 - t123;
t411 = Icges(5,4) * t224;
t174 = Icges(5,1) * t225 - t411;
t125 = Icges(5,5) * t264 + t174 * t265;
t171 = Icges(5,2) * t225 + t411;
t454 = -t171 * t265 + t125;
t453 = -t171 + t174;
t452 = t173 + t301;
t200 = rSges(5,2) * t395;
t126 = rSges(5,1) * t393 - t265 * rSges(5,3) - t200;
t251 = t264 * rSges(5,3);
t392 = t225 * t265;
t394 = t224 * t265;
t127 = rSges(5,1) * t392 - rSges(5,2) * t394 + t251;
t421 = rSges(5,2) * t225;
t175 = rSges(5,1) * t224 + t421;
t151 = t175 * t264;
t152 = t175 * t265;
t260 = qJD(3) + qJD(4);
t202 = t260 * t264;
t203 = t260 * t265;
t262 = cos(pkin(7));
t226 = t262 * pkin(2) + pkin(1);
t263 = -pkin(5) - qJ(2);
t258 = -pkin(6) + t263;
t236 = cos(t259);
t223 = pkin(3) * t236;
t344 = t223 + t226;
t368 = -t265 * t258 - t264 * t344;
t383 = t265 * t263;
t95 = t264 * t226 + t368 + t383;
t178 = t265 * t344;
t211 = t265 * t226;
t356 = t263 - t258;
t96 = t264 * t356 + t178 - t211;
t38 = t126 * t202 + t127 * t203 + (-t264 * t95 + t265 * t96) * qJD(3);
t240 = qJD(2) * t264;
t349 = qJD(3) * t265;
t339 = t235 * t349;
t316 = pkin(3) * t339;
t291 = t240 - t316;
t281 = -t203 * t175 + t291;
t244 = t265 * qJ(2);
t427 = pkin(1) - t226;
t289 = t264 * t427 - t383;
t128 = -t244 + t289;
t204 = pkin(1) * t264 - t244;
t373 = t128 - t204;
t416 = -t126 + t95;
t315 = t373 + t416;
t41 = qJD(1) * t315 + t281;
t350 = qJD(3) * t264;
t430 = pkin(3) * t235;
t198 = t350 * t430;
t241 = qJD(2) * t265;
t363 = t198 + t241;
t242 = t264 * qJ(2);
t206 = t265 * pkin(1) + t242;
t319 = -t263 * t264 + t211;
t129 = t319 - t206;
t372 = t129 + t206;
t415 = t127 + t96;
t42 = -t175 * t202 + (t372 + t415) * qJD(1) - t363;
t217 = t225 * rSges(5,1);
t422 = rSges(5,2) * t224;
t449 = t217 - t422;
t347 = qJD(1) * qJD(3);
t190 = qJDD(3) * t264 + t265 * t347;
t346 = qJD(1) * qJD(4);
t143 = qJDD(4) * t264 + t265 * t346 + t190;
t227 = t264 * t347;
t144 = t264 * t346 + t227 + (-qJDD(3) - qJDD(4)) * t265;
t191 = -qJDD(3) * t265 + t227;
t352 = qJD(1) * t264;
t351 = qJD(1) * t265;
t365 = rSges(5,3) * t351 + qJD(1) * t200;
t76 = -t203 * t421 + (-t203 * t224 - t225 * t352) * rSges(5,1) + t365;
t77 = -t260 * t151 + (t265 * t449 + t251) * qJD(1);
t318 = t226 - t344;
t79 = -t316 + (t264 * t318 + t265 * t356) * qJD(1);
t218 = t263 * t352;
t386 = t264 * t258;
t80 = t218 - t198 + (-t265 * t318 - t386) * qJD(1);
t8 = t126 * t143 - t127 * t144 - t190 * t95 - t191 * t96 + t202 * t77 + t203 * t76 + (t264 * t80 + t265 * t79) * qJD(3);
t451 = -t41 * (qJD(1) * t151 - t203 * t449) - t38 * (-t202 * t151 - t152 * t203) - t42 * (-qJD(1) * t152 - t202 * t449) + t8 * (t264 * t126 + t265 * t127);
t193 = qJD(1) * t204;
t450 = qJD(1) * t128 - t193;
t222 = Icges(4,4) * t236;
t302 = -Icges(4,2) * t235 + t222;
t183 = Icges(4,1) * t235 + t222;
t261 = sin(pkin(7));
t424 = rSges(3,2) * t261;
t426 = rSges(3,1) * t262;
t154 = t264 * rSges(3,3) + (-t424 + t426) * t265;
t180 = Icges(4,5) * t236 - Icges(4,6) * t235;
t179 = Icges(4,5) * t235 + Icges(4,6) * t236;
t285 = qJD(3) * t179;
t412 = Icges(4,4) * t235;
t184 = Icges(4,1) * t236 - t412;
t135 = Icges(4,5) * t264 + t184 * t265;
t133 = Icges(4,6) * t264 + t265 * t302;
t400 = t133 * t235;
t297 = -t135 * t236 + t400;
t404 = Icges(4,3) * t265;
t448 = -t265 * t285 + (-t180 * t264 + t297 + t404) * qJD(1);
t391 = t235 * t264;
t210 = Icges(4,4) * t391;
t389 = t236 * t264;
t410 = Icges(4,5) * t265;
t134 = Icges(4,1) * t389 - t210 - t410;
t406 = Icges(4,6) * t265;
t132 = Icges(4,4) * t389 - Icges(4,2) * t391 - t406;
t401 = t132 * t235;
t298 = -t134 * t236 + t401;
t131 = Icges(4,3) * t264 + t180 * t265;
t354 = qJD(1) * t131;
t447 = qJD(1) * t298 - t264 * t285 + t354;
t170 = Icges(5,5) * t225 - Icges(5,6) * t224;
t169 = Icges(5,5) * t224 + Icges(5,6) * t225;
t398 = t169 * t265;
t402 = t123 * t224;
t403 = Icges(5,3) * t265;
t446 = -t260 * t398 + (-t125 * t225 - t170 * t264 + t402 + t403) * qJD(1);
t197 = Icges(5,4) * t395;
t409 = Icges(5,5) * t265;
t124 = Icges(5,1) * t393 - t197 - t409;
t300 = t122 * t224 - t124 * t225;
t121 = Icges(5,3) * t264 + t170 * t265;
t355 = qJD(1) * t121;
t399 = t169 * t264;
t445 = qJD(1) * t300 - t260 * t399 + t355;
t130 = Icges(4,5) * t389 - Icges(4,6) * t391 - t404;
t50 = -t265 * t130 - t264 * t298;
t295 = t171 * t224 - t173 * t225;
t444 = qJD(1) * t295 + t170 * t260;
t181 = Icges(4,2) * t236 + t412;
t293 = t181 * t235 - t183 * t236;
t443 = t293 * qJD(1) + t180 * qJD(3);
t442 = t264 * (-t181 * t265 + t135) - t265 * (-Icges(4,2) * t389 + t134 - t210);
t441 = qJD(1) * t452 + t202 * t454 - t203 * (-Icges(5,2) * t393 + t124 - t197);
t440 = t143 / 0.2e1;
t439 = t144 / 0.2e1;
t438 = t190 / 0.2e1;
t437 = t191 / 0.2e1;
t436 = -t202 / 0.2e1;
t435 = t202 / 0.2e1;
t434 = -t203 / 0.2e1;
t433 = t203 / 0.2e1;
t432 = t264 / 0.2e1;
t431 = -t265 / 0.2e1;
t429 = -qJD(1) / 0.2e1;
t428 = qJD(1) / 0.2e1;
t425 = rSges(4,1) * t236;
t423 = rSges(4,2) * t236;
t185 = rSges(4,1) * t235 + t423;
t162 = t185 * t265;
t252 = t264 * rSges(4,3);
t388 = t236 * t265;
t390 = t235 * t265;
t141 = rSges(4,1) * t388 - rSges(4,2) * t390 + t252;
t57 = -t185 * t350 - t241 + (t141 + t372) * qJD(1);
t420 = t162 * t57;
t313 = -t185 * t349 + t240;
t361 = rSges(4,2) * t391 + t265 * rSges(4,3);
t140 = rSges(4,1) * t389 - t361;
t341 = -t140 + t373;
t56 = qJD(1) * t341 + t313;
t419 = t264 * t56;
t418 = t41 * t175;
t417 = qJDD(1) / 0.2e1;
t397 = t179 * t264;
t396 = t179 * t265;
t387 = t236 * qJD(3) ^ 2;
t120 = Icges(5,5) * t393 - Icges(5,6) * t395 - t403;
t385 = t265 * t120;
t59 = -t264 * t295 - t398;
t382 = t59 * qJD(1);
t64 = -t264 * t293 - t396;
t381 = t64 * qJD(1);
t380 = -t264 * t120 - t124 * t392;
t379 = t264 * t121 + t125 * t392;
t378 = -t264 * t130 - t134 * t388;
t377 = t264 * t131 + t135 * t388;
t167 = qJD(1) * t206 - t241;
t376 = t218 - (-t265 * t427 - t242) * qJD(1) - t167;
t114 = t154 + t206;
t367 = -t181 + t184;
t366 = t183 + t302;
t364 = rSges(4,3) * t351 + t352 * t457;
t343 = t264 * t426;
t219 = t264 * t424;
t359 = t265 * rSges(3,3) + t219;
t153 = t343 - t359;
t362 = -t204 - t153;
t360 = rSges(3,3) * t351 + qJD(1) * t219;
t348 = qJD(1) * qJD(2);
t358 = qJDD(2) * t264 + t265 * t348;
t231 = qJ(2) * t351;
t357 = t231 + t240;
t353 = qJD(1) * t180;
t345 = t126 * t351 + t264 * t77 + t265 * t76;
t338 = -pkin(1) - t426;
t336 = t352 / 0.2e1;
t335 = t351 / 0.2e1;
t334 = -t350 / 0.2e1;
t333 = t350 / 0.2e1;
t332 = -t349 / 0.2e1;
t331 = t349 / 0.2e1;
t284 = -t175 - t430;
t330 = -t226 - t425;
t99 = t125 * t393;
t329 = t265 * t121 - t99;
t328 = qJD(1) * t125 + t456 * t260;
t327 = (-t174 * t264 + t409) * qJD(1) + t455 * t260;
t326 = qJD(1) * t123 + t124 * t260 - t171 * t202;
t325 = (-t264 * t301 + t405) * qJD(1) + t454 * t260;
t104 = t135 * t389;
t324 = t265 * t131 - t104;
t323 = -t120 + t402;
t322 = -t130 + t400;
t321 = t452 * t260;
t320 = t453 * t260;
t142 = t449 * t260;
t314 = -qJD(3) * t223 - t142;
t207 = rSges(2,1) * t265 - rSges(2,2) * t264;
t205 = rSges(2,1) * t264 + rSges(2,2) * t265;
t186 = t425 - t457;
t69 = t133 * t236 + t135 * t235;
t286 = qJD(3) * t181;
t83 = -t265 * t286 + (-t264 * t302 + t406) * qJD(1);
t287 = qJD(3) * t183;
t85 = -t265 * t287 + (-t184 * t264 + t410) * qJD(1);
t271 = -qJD(3) * t69 - t235 * t83 + t236 * t85 + t354;
t68 = t132 * t236 + t134 * t235;
t84 = qJD(1) * t133 - t264 * t286;
t86 = qJD(1) * t135 - t264 * t287;
t272 = qJD(1) * t130 - qJD(3) * t68 - t235 * t84 + t236 * t86;
t311 = -(t264 * t447 + t272 * t265) * t265 + (t264 * t448 + t271 * t265) * t264;
t310 = -(t272 * t264 - t265 * t447) * t265 + (t271 * t264 - t265 * t448) * t264;
t309 = -t264 * t42 - t265 * t41;
t51 = -t133 * t391 - t324;
t308 = t264 * t51 - t265 * t50;
t52 = -t132 * t390 - t378;
t53 = -t133 * t390 + t377;
t307 = t264 * t53 - t265 * t52;
t306 = -t264 * t57 - t265 * t56;
t87 = -t349 * t423 + (-t236 * t352 - t339) * rSges(4,1) + t364;
t161 = t185 * t264;
t88 = -qJD(3) * t161 + (t186 * t265 + t252) * qJD(1);
t305 = t264 * t88 + t265 * t87;
t61 = t122 * t225 + t124 * t224;
t296 = t140 * t264 + t141 * t265;
t294 = t181 * t236 + t183 * t235;
t292 = -t344 - t217;
t290 = -qJDD(2) * t265 + qJD(1) * (-pkin(1) * t352 + t357) + qJDD(1) * t206 + t264 * t348;
t288 = t300 * t264;
t283 = qJD(1) * t170 - t202 * t398 + t203 * t399;
t282 = qJD(1) * (qJD(1) * t289 - t231) + qJDD(1) * t129 + t290;
t280 = t132 * t265 - t133 * t264;
t279 = (-t235 * t366 + t236 * t367) * qJD(1);
t275 = qJD(1) * t120 - t224 * t326 + t225 * t328;
t11 = t264 * t445 + t275 * t265;
t274 = -t224 * t325 + t225 * t327 + t355;
t12 = t264 * t446 + t274 * t265;
t13 = t275 * t264 - t265 * t445;
t14 = t274 * t264 - t265 * t446;
t44 = -t288 - t385;
t45 = -t123 * t395 - t329;
t22 = t202 * t45 - t203 * t44 + t382;
t46 = -t122 * t394 - t380;
t47 = -t123 * t394 + t379;
t60 = -t265 * t295 + t399;
t58 = t60 * qJD(1);
t23 = t202 * t47 - t203 * t46 + t58;
t276 = t453 * qJD(1) + t455 * t202 - t456 * t203;
t268 = -t224 * t441 + t276 * t225;
t273 = qJD(1) * t169 - t224 * t321 + t225 * t320;
t28 = t264 * t444 + t273 * t265;
t29 = t273 * t264 - t265 * t444;
t32 = t224 * t328 + t225 * t326;
t33 = t224 * t327 + t225 * t325;
t62 = t123 * t225 + t125 * t224;
t278 = (qJD(1) * t28 + qJDD(1) * t60 - t11 * t203 + t12 * t202 + t143 * t47 + t144 * t46) * t432 + (t276 * t224 + t225 * t441) * t429 + (qJD(1) * t29 + qJDD(1) * t59 - t13 * t203 + t14 * t202 + t143 * t45 + t144 * t44) * t431 + t22 * t336 + t23 * t335 + (t264 * t47 - t265 * t46) * t440 + (t264 * t45 - t265 * t44) * t439 + (-t11 * t265 + t12 * t264 + (t264 * t46 + t265 * t47) * qJD(1)) * t435 + (-t13 * t265 + t14 * t264 + (t264 * t44 + t265 * t45) * qJD(1)) * t434 + (t264 * t62 - t265 * t61) * t417 + (t264 * t33 - t265 * t32 + (t264 * t61 + t265 * t62) * qJD(1)) * t428 + (t264 * t283 + t265 * t268) * t436 + (t264 * t268 - t265 * t283) * t433;
t165 = t302 * qJD(3);
t166 = t184 * qJD(3);
t270 = qJD(1) * t179 - qJD(3) * t294 - t165 * t235 + t166 * t236;
t269 = -t235 * t442 + t280 * t236;
t168 = t186 * qJD(3);
t98 = qJD(1) * t114 - t241;
t97 = qJD(1) * t362 + t240;
t78 = t296 * qJD(3);
t65 = -t265 * t293 + t397;
t63 = t65 * qJD(1);
t55 = qJDD(1) * t154 + qJD(1) * (-qJD(1) * t343 + t360) + t290;
t54 = t362 * qJDD(1) + (-qJD(1) * t154 - t167) * qJD(1) + t358;
t37 = -qJD(3) * t297 + t235 * t85 + t236 * t83;
t36 = -t298 * qJD(3) + t235 * t86 + t236 * t84;
t35 = t270 * t264 - t265 * t443;
t34 = t264 * t443 + t270 * t265;
t31 = qJD(1) * t87 + qJDD(1) * t141 - t168 * t350 - t185 * t190 + t282;
t30 = -t168 * t349 + t185 * t191 + t341 * qJDD(1) + (-t88 + t376) * qJD(1) + t358;
t25 = qJD(3) * t307 + t63;
t24 = qJD(3) * t308 + t381;
t10 = -t142 * t202 - t143 * t175 + t415 * qJDD(1) + (t76 + t79) * qJD(1) + (-t190 * t235 - t264 * t387) * pkin(3) + t282;
t9 = -t142 * t203 + t144 * t175 + (t191 * t235 - t265 * t387) * pkin(3) + t315 * qJDD(1) + (-t77 - t80 + t376) * qJD(1) + t358;
t1 = [(t63 + ((t51 - t104 + (t131 + t401) * t265 + t378) * t265 + t377 * t264) * qJD(3)) * t331 - m(2) * (-g(1) * t205 + g(2) * t207) + (t58 + (t45 + (t122 * t265 + t123 * t264) * t224 + t329 + t380) * t203 + (-t124 * t393 + t385 + t44 + (t122 * t264 - t123 * t265) * t224 + t379) * t202) * t433 + (t62 + t60) * t440 + (t61 + t59) * t439 + (t69 + t65) * t438 + (t68 + t64) * t437 + (-t382 + (t47 - t288 - t379) * t203 + (t323 * t264 + t46 - t99) * t202 + ((t121 + t300) * t202 + t323 * t203) * t265 + t22) * t436 + (t33 + t28) * t435 + (t24 - t381 + ((t265 * t322 - t377 + t53) * t265 + (t264 * t322 + t324 + t52) * t264) * qJD(3)) * t334 + (t34 + t37) * t333 + (-qJD(3) * t293 + t165 * t236 + t166 * t235 + t224 * t320 + t225 * t321) * qJD(1) + (t32 + t29 + t23) * t434 + (t36 + t35 + t25) * t332 + (-(qJD(1) * t416 + t281 - t41 + t450) * t42 + t41 * t363 + t42 * (t291 + t365) + (-t152 * t42 + t264 * t418) * t260 + ((t41 * (t292 + t422) - t42 * t258) * t265 + (t41 * (t258 - rSges(5,3)) + t42 * t292) * t264) * qJD(1) + (t10 - g(2)) * (t127 + t178 - t386) + (t9 - g(1)) * (-t126 + t368)) * m(5) + (t56 * (t218 + t241) + t57 * (t240 + t364) + (t185 * t419 - t420) * qJD(3) + ((-t56 * rSges(4,3) + t330 * t57) * t264 + (t56 * (-t186 - t226) - t57 * t263) * t265) * qJD(1) - (-qJD(1) * t140 + t313 + t450 - t56) * t57 + (t31 - g(2)) * (t141 + t319) + (t30 - g(1)) * (t330 * t264 + t361 - t383)) * m(4) + (t97 * t241 + t98 * (t357 + t360) + (t97 * (t338 + t424) * t265 + (t97 * (-rSges(3,3) - qJ(2)) + t98 * t338) * t264) * qJD(1) - (-qJD(1) * t153 - t193 + t240 - t97) * t98 + (t55 - g(2)) * t114 + (t54 - g(1)) * (t338 * t264 + t244 + t359)) * m(3) + (m(2) * (t205 ^ 2 + t207 ^ 2) + t294 + t171 * t225 + t173 * t224 + Icges(2,3) + Icges(3,2) * t262 ^ 2 + (Icges(3,1) * t261 + 0.2e1 * Icges(3,4) * t262) * t261) * qJDD(1); (-m(3) - m(4) - m(5)) * (g(1) * t264 - g(2) * t265) + 0.2e1 * (t10 * t431 + t432 * t9) * m(5) + 0.2e1 * (t30 * t432 + t31 * t431) * m(4) + 0.2e1 * (t431 * t55 + t432 * t54) * m(3); ((t235 * t367 + t236 * t366) * qJD(1) + (t280 * t235 + t236 * t442) * qJD(3)) * t429 + ((t52 * t264 + t53 * t265) * qJD(1) + t311) * t333 + ((t50 * t264 + t51 * t265) * qJD(1) + t310) * t332 + t278 + (t264 * t37 - t265 * t36 + (t68 * t264 + t265 * t69) * qJD(1)) * t428 + ((-t349 * t397 - t353) * t265 + (t279 + (t265 * t396 + t269) * qJD(3)) * t264) * t331 + ((-t350 * t396 + t353) * t264 + (t279 + (t264 * t397 + t269) * qJD(3)) * t265) * t334 + t24 * t336 + t25 * t335 + (qJD(1) * t34 + qJD(3) * t311 + qJDD(1) * t65 + t190 * t53 + t191 * t52) * t432 + (qJD(1) * t35 + qJD(3) * t310 + qJDD(1) * t64 + t190 * t51 + t191 * t50) * t431 + (t264 * t69 - t265 * t68) * t417 + t307 * t438 + t308 * t437 + (t38 * t345 + (t9 * t284 + t41 * t314 + t8 * t96 + t38 * t79 + (t284 * t42 - t38 * t95) * qJD(1)) * t265 + (t10 * t284 + t42 * t314 - t8 * t95 + t38 * t80 + (-t38 * t415 + t418) * qJD(1)) * t264 - g(3) * (t449 + t223) - (g(1) * t265 + g(2) * t264) * t284 - (-t42 * t235 * t351 + (t309 * t236 + t38 * (-t264 ^ 2 - t265 ^ 2) * t235) * qJD(3)) * pkin(3) + t451) * m(5) + (-(t161 * t56 - t420) * qJD(1) - (t78 * (-t161 * t264 - t162 * t265) + t306 * t186) * qJD(3) + (qJD(3) * t305 + t140 * t190 - t141 * t191) * t296 + t78 * ((t140 * t265 - t141 * t264) * qJD(1) + t305) + t306 * t168 + (-t264 * t31 - t265 * t30 + (-t265 * t57 + t419) * qJD(1)) * t185 + g(1) * t162 + g(2) * t161 - g(3) * t186) * m(4); t278 + (t38 * (-t127 * t352 + t345) + t309 * t142 + (-t10 * t264 - t9 * t265 + (t264 * t41 - t265 * t42) * qJD(1)) * t175 + g(1) * t152 + g(2) * t151 - g(3) * t449 + t451) * m(5);];
tau = t1;
