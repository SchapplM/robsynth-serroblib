% Calculate vector of inverse dynamics joint torques for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:52
% DurationCPUTime: 12.61s
% Computational Cost: add. (6007->450), mult. (7444->551), div. (0->0), fcn. (5785->6), ass. (0->246)
t213 = sin(qJ(3));
t215 = cos(qJ(3));
t153 = Icges(4,5) * t215 - Icges(4,6) * t213;
t212 = qJ(1) + pkin(6);
t203 = sin(t212);
t204 = cos(t212);
t89 = Icges(4,3) * t203 + t153 * t204;
t155 = Icges(5,4) * t215 + Icges(5,6) * t213;
t91 = Icges(5,2) * t203 + t155 * t204;
t446 = t89 + t91;
t206 = Icges(5,5) * t213;
t248 = Icges(5,1) * t215 + t206;
t95 = Icges(5,4) * t203 + t204 * t248;
t339 = Icges(4,4) * t213;
t161 = Icges(4,1) * t215 - t339;
t97 = Icges(4,5) * t203 + t161 * t204;
t442 = t95 + t97;
t156 = Icges(4,2) * t215 + t339;
t246 = Icges(5,3) * t215 - t206;
t445 = t156 + t246;
t336 = Icges(5,5) * t215;
t158 = Icges(5,1) * t213 - t336;
t207 = Icges(4,4) * t215;
t160 = Icges(4,1) * t213 + t207;
t444 = t158 + t160;
t319 = t204 * t215;
t320 = t204 * t213;
t176 = Icges(5,5) * t319;
t332 = Icges(5,6) * t203;
t87 = Icges(5,3) * t320 + t176 + t332;
t443 = t446 * t203 + t442 * t319 + t87 * t320;
t436 = -t445 * t213 + t444 * t215;
t321 = t203 * t215;
t322 = t203 * t213;
t333 = Icges(4,6) * t204;
t92 = Icges(4,4) * t321 - Icges(4,2) * t322 - t333;
t355 = t213 * t92;
t177 = Icges(4,4) * t322;
t337 = Icges(4,5) * t204;
t96 = Icges(4,1) * t321 - t177 - t337;
t253 = -t215 * t96 + t355;
t90 = -Icges(5,2) * t204 + t155 * t203;
t356 = t204 * t90;
t151 = Icges(5,3) * t213 + t336;
t86 = -Icges(5,6) * t204 + t151 * t203;
t94 = -Icges(5,4) * t204 + t203 * t248;
t255 = t213 * t86 + t215 * t94;
t394 = t203 * t255;
t30 = -t356 + t394;
t330 = Icges(4,3) * t204;
t88 = Icges(4,5) * t321 - Icges(4,6) * t322 - t330;
t405 = -t203 * t253 - t204 * t88 + t30;
t247 = -Icges(4,2) * t213 + t207;
t93 = Icges(4,6) * t203 + t204 * t247;
t403 = -t320 * t93 + t443;
t271 = t204 * t91 - t95 * t321 - t87 * t322;
t75 = t97 * t321;
t276 = t204 * t89 - t75;
t33 = -t322 * t93 - t276;
t399 = -t271 + t33;
t441 = (t151 - t247) * qJD(3);
t440 = (t161 + t248) * qJD(3);
t152 = Icges(4,5) * t213 + Icges(4,6) * t215;
t154 = Icges(5,4) * t213 - Icges(5,6) * t215;
t439 = t152 + t154;
t438 = t445 * qJD(3);
t437 = t444 * qJD(3);
t435 = -t444 * t213 - t445 * t215;
t323 = t154 * t204;
t325 = t152 * t204;
t413 = t436 * t203 - t323 - t325;
t324 = t154 * t203;
t326 = t152 * t203;
t412 = t436 * t204 + t324 + t326;
t402 = (-t86 + t92) * t215 + (t94 + t96) * t213;
t401 = (-t87 + t93) * t215 + t442 * t213;
t434 = t438 * t204 + (t203 * t247 - t333 - t86) * qJD(1);
t433 = t438 * t203 + (t151 * t204 + t332 - t93) * qJD(1);
t432 = -t437 * t204 + (-t161 * t203 + t337 - t94) * qJD(1);
t431 = -qJD(1) * t442 + t203 * t437;
t430 = qJD(1) * t439 + qJD(3) * t435 + t213 * t441 + t215 * t440;
t429 = t439 * qJD(3);
t354 = t213 * t93;
t428 = t213 * t87 + t215 * t442 - t354;
t427 = t253 - t255;
t366 = -t203 * t88 - t96 * t319;
t36 = -t320 * t92 - t366;
t84 = t203 * t90;
t34 = t94 * t319 + t86 * t320 + t84;
t400 = t204 * t34;
t426 = t403 * t203 - t204 * t36 - t400;
t425 = t399 * t203 - t405 * t204;
t424 = (-t153 - t155) * qJD(3) + t436 * qJD(1);
t216 = cos(qJ(1));
t211 = t216 * pkin(1);
t423 = t412 * qJD(1);
t133 = rSges(3,1) * t203 + rSges(3,2) * t204;
t214 = sin(qJ(1));
t373 = pkin(1) * t214;
t125 = -t133 - t373;
t422 = t413 * qJD(1);
t421 = t446 * qJD(1);
t420 = t204 ^ 2;
t217 = qJD(1) ^ 2;
t292 = t217 * t211;
t419 = qJD(3) * t425 + t422;
t418 = qJD(3) * t426 + t423;
t417 = qJD(3) * t427 + t213 * t431 + t215 * t433;
t416 = qJD(3) * t428 + t213 * t432 - t215 * t434;
t415 = -t203 * t424 + t204 * t430;
t414 = t203 * t430 + t204 * t424;
t404 = t34 + t36;
t411 = -qJD(3) * t401 + t213 * t434 + t215 * t432 + t421;
t395 = qJD(1) * t90;
t410 = -qJD(1) * t88 + qJD(3) * t402 - t433 * t213 + t431 * t215 - t395;
t409 = t356 + t443;
t408 = qJD(1) * t427 - t203 * t429 + t421;
t407 = -t395 - t429 * t204 + (-t153 * t203 + t330 - t428) * qJD(1);
t406 = rSges(5,3) + qJ(4);
t374 = rSges(5,1) + pkin(3);
t194 = t203 * rSges(4,3);
t101 = rSges(4,1) * t319 - rSges(4,2) * t320 + t194;
t136 = t204 * pkin(2) + t203 * pkin(5);
t390 = t211 + t136;
t72 = t101 + t390;
t205 = t213 * qJ(4);
t263 = t215 * pkin(3) + t205;
t264 = t215 * rSges(5,1) + t213 * rSges(5,3);
t398 = t263 + t264;
t397 = t408 * t420 + (t411 * t203 + (-t407 + t410) * t204) * t203;
t396 = t410 * t420 + (t407 * t203 + (-t408 + t411) * t204) * t203;
t195 = t203 * rSges(5,2);
t315 = rSges(5,3) * t320 + t204 * t205 + t319 * t374 + t195;
t42 = t390 + t315;
t393 = t213 * t374;
t134 = t204 * rSges(3,1) - rSges(3,2) * t203;
t126 = t134 + t211;
t392 = t406 * t321;
t391 = t406 * t319;
t295 = qJD(4) * t213;
t162 = t204 * t295;
t296 = qJD(3) * t204;
t287 = t215 * t296;
t298 = qJD(1) * t204;
t389 = rSges(5,2) * t298 + t287 * t406 + t162;
t388 = -g(1) * t204 - g(2) * t203;
t344 = -t160 * t203 - t92;
t348 = -Icges(4,2) * t321 - t177 + t96;
t380 = -t213 * t348 + t215 * t344;
t346 = -t158 * t203 + t86;
t350 = t246 * t203 - t94;
t379 = t213 * t350 + t215 * t346;
t293 = qJD(1) * qJD(3);
t129 = qJDD(3) * t203 + t204 * t293;
t378 = t129 / 0.2e1;
t130 = -qJDD(3) * t204 + t203 * t293;
t377 = t130 / 0.2e1;
t288 = t213 * t296;
t299 = qJD(1) * t203;
t230 = -t215 * t299 - t288;
t290 = t213 * t299;
t368 = t230 * t374 - t290 * t406 + t389;
t163 = pkin(3) * t213 - qJ(4) * t215;
t164 = rSges(5,1) * t213 - rSges(5,3) * t215;
t297 = qJD(3) * t203;
t367 = t263 * t298 + (-qJD(3) * t163 + t295) * t203 - t164 * t297 + (t204 * t264 + t195) * qJD(1);
t364 = rSges(4,1) * t215;
t165 = rSges(4,1) * t213 + rSges(4,2) * t215;
t122 = t165 * t204;
t40 = qJD(1) * t72 - t165 * t297;
t362 = t122 * t40;
t200 = t204 * pkin(5);
t135 = pkin(2) * t203 - t200;
t302 = rSges(4,2) * t322 + t204 * rSges(4,3);
t99 = rSges(4,1) * t321 - t302;
t285 = -t99 - t373;
t272 = -t135 + t285;
t289 = t165 * t296;
t39 = qJD(1) * t272 - t289;
t360 = t203 * t39;
t304 = t163 + t164;
t274 = t304 * qJD(3);
t240 = -t204 * t274 + t162;
t197 = t204 * rSges(5,2);
t342 = t203 * t398 - t197;
t273 = -t342 - t373;
t241 = -t135 + t273;
t28 = qJD(1) * t241 + t240;
t359 = t204 * t28;
t358 = t204 * t39;
t294 = qJD(4) * t215;
t27 = -t294 + qJD(2) + (t203 * t342 + t204 * t315) * qJD(3);
t353 = t27 * t213;
t349 = t246 * t204 - t95;
t347 = -t156 * t204 + t97;
t345 = -Icges(5,1) * t320 + t176 + t87;
t343 = -t160 * t204 - t93;
t314 = -t322 * t374 + t392;
t313 = -t320 * t374 + t391;
t312 = -qJD(3) * t398 + t294;
t309 = rSges(4,2) * t290 + rSges(4,3) * t298;
t308 = -t246 + t248;
t307 = t151 - t158;
t306 = -t156 + t161;
t305 = t160 + t247;
t301 = qJD(1) * t153;
t300 = qJD(1) * t155;
t286 = -pkin(2) - t364;
t282 = -t297 / 0.2e1;
t281 = t297 / 0.2e1;
t280 = -t296 / 0.2e1;
t279 = t296 / 0.2e1;
t277 = t200 - t373;
t275 = -t88 + t354;
t270 = qJDD(1) * t211 - t217 * t373;
t170 = rSges(2,1) * t216 - rSges(2,2) * t214;
t166 = rSges(2,1) * t214 + rSges(2,2) * t216;
t169 = -rSges(4,2) * t213 + t364;
t29 = (-t274 + t295) * t203 + t42 * qJD(1);
t262 = t203 * t29 + t359;
t257 = -t203 * t40 - t358;
t68 = rSges(4,1) * t230 - rSges(4,2) * t287 + t309;
t118 = t165 * t203;
t70 = -qJD(3) * t118 + (t169 * t204 + t194) * qJD(1);
t256 = t203 * t70 + t204 * t68;
t251 = t101 * t204 + t203 * t99;
t187 = pkin(5) * t298;
t239 = qJD(1) * (-pkin(2) * t299 + t187) + qJDD(1) * t136 + t270;
t238 = t262 * t215;
t229 = t213 * t349 + t215 * t345;
t228 = -t213 * t347 + t215 * t343;
t227 = -t213 * t406 - t215 * t374 - pkin(2);
t226 = qJDD(4) * t213 + (t294 + t312) * qJD(3);
t225 = (t213 * t307 + t215 * t308) * qJD(1);
t224 = (-t213 * t305 + t215 * t306) * qJD(1);
t144 = t169 * qJD(3);
t132 = qJD(1) * t135;
t128 = t136 * qJD(1);
t38 = qJD(3) * t251 + qJD(2);
t22 = qJD(1) * t68 + qJDD(1) * t101 - t129 * t165 - t144 * t297 + t239;
t21 = -t292 - t144 * t296 + t130 * t165 + (-t128 - t70) * qJD(1) + t272 * qJDD(1);
t16 = qJD(3) * t256 - t101 * t130 + t129 * t99 + qJDD(2);
t11 = -t304 * t129 + t315 * qJDD(1) + (t162 + t368) * qJD(1) + t226 * t203 + t239;
t10 = -t292 + t304 * t130 + t226 * t204 + t241 * qJDD(1) + (-t203 * t295 - t128 - t367) * qJD(1);
t1 = -qJDD(4) * t215 + qJDD(2) - t315 * t130 + t342 * t129 + (t203 * t367 + t204 * t368 + t295) * qJD(3);
t2 = [-m(2) * (-g(1) * t166 + g(2) * t170) + ((-t400 + (t33 - t75 + (t89 + t355) * t204 + t366) * t204 + (t30 - t394 + t409) * t203) * qJD(3) + t423) * t279 + (t436 * qJD(3) + t440 * t213 - t441 * t215) * qJD(1) + ((-t133 * t217 - g(2) + t270) * t126 + (-t292 + (-0.2e1 * t134 - t211 + t126) * t217 - g(1)) * t125) * m(3) + (t29 * (-t374 * t288 + t187 + t389) + (-t295 + (-t215 * t406 + t393) * qJD(3)) * t203 * t28 + ((-t214 * t29 - t216 * t28) * pkin(1) + t227 * t359 + (t28 * (-rSges(5,2) - pkin(5)) + t29 * (-pkin(2) - t398)) * t203) * qJD(1) - (qJD(1) * t273 - t132 + t240 - t28) * t29 + (t11 - g(2)) * t42 + (t10 - g(1)) * (t203 * t227 + t197 + t277)) * m(5) + (t40 * (t187 + t309) + (t165 * t360 - t362) * qJD(3) + ((-t214 * t40 - t216 * t39) * pkin(1) + (-pkin(2) - t169) * t358 + (t39 * (-rSges(4,3) - pkin(5)) + t40 * t286) * t203) * qJD(1) - (qJD(1) * t285 - t132 - t289 - t39) * t40 + (t22 - g(2)) * t72 + (t21 - g(1)) * (t286 * t203 + t277 + t302)) * m(4) + (t412 + t401) * t378 + (t413 + t402) * t377 + (((t204 * t275 + t403 - t409) * t204 + (t203 * t275 + t271 + t276 + t404 - t84) * t203) * qJD(3) + t419 - t422) * t282 + (t415 + t416) * t281 + (m(3) * (t125 ^ 2 + t134 * t126) + m(2) * (t166 ^ 2 + t170 ^ 2) + Icges(2,3) + Icges(3,3) - t435) * qJDD(1) + (t414 - t417 + t418) * t280; m(3) * qJDD(2) + m(4) * t16 + m(5) * t1 + (-m(3) - m(4) - m(5)) * g(3); t426 * t378 + t425 * t377 + (qJD(1) * t415 + t396 * qJD(3) + qJDD(1) * t412 + t403 * t129 + t404 * t130) * t203 / 0.2e1 - (qJD(1) * t414 + t397 * qJD(3) + qJDD(1) * t413 + t399 * t129 + t405 * t130) * t204 / 0.2e1 - ((((-t348 + t350) * t204 + (t347 - t349) * t203) * t215 + ((-t344 - t346) * t204 + (t343 + t345) * t203) * t213) * qJD(3) + ((t305 - t307) * t215 + (t306 + t308) * t213) * qJD(1)) * qJD(1) / 0.2e1 + (t417 * t204 + t416 * t203 + (t402 * t203 + t401 * t204) * qJD(1)) * qJD(1) / 0.2e1 + (t401 * t203 - t402 * t204) * qJDD(1) / 0.2e1 + t419 * t299 / 0.2e1 + t418 * t298 / 0.2e1 + ((-t297 * t323 + t300) * t203 + (t225 + (-t379 * t204 + (t324 + t229) * t203) * qJD(3)) * t204 + (-t297 * t325 + t301) * t203 + (t224 + (-t380 * t204 + (t326 + t228) * t203) * qJD(3)) * t204) * t282 + ((t404 * t203 + t204 * t403) * qJD(1) + t396) * t281 + ((t203 * t405 + t399 * t204) * qJD(1) + t397) * t280 + ((-t296 * t324 - t300) * t204 + (t225 + (t229 * t203 + (t323 - t379) * t204) * qJD(3)) * t203 + (-t296 * t326 - t301) * t204 + (t224 + (t228 * t203 + (t325 - t380) * t204) * qJD(3)) * t203) * t279 + (-g(1) * t391 - g(2) * t392 - g(3) * t398 - t388 * t393 + (-t10 * t304 + t28 * t312 + t1 * t315 + t27 * t368 + (t27 * t342 - t29 * t304) * qJD(1)) * t204 + (-t11 * t304 + t29 * t312 + t1 * t342 + t27 * t367 + (-t27 * t315 + t28 * t304) * qJD(1)) * t203 - (t238 + t353) * qJD(4) - (-t28 * t314 + t29 * t313) * qJD(1) - ((t27 * t313 - t28 * t398) * t204 + (t27 * t314 - t29 * t398) * t203) * qJD(3)) * m(5) + (g(1) * t122 + g(2) * t118 - g(3) * t169 - (t118 * t39 - t362) * qJD(1) - (t38 * (-t118 * t203 - t122 * t204) + t257 * t169) * qJD(3) + t16 * t251 + t38 * ((-t101 * t203 + t204 * t99) * qJD(1) + t256) + t257 * t144 + (-t22 * t203 - t21 * t204 + (-t204 * t40 + t360) * qJD(1)) * t165) * m(4); (-(t238 + (t203 ^ 2 + t420) * t353) * qJD(3) + (qJD(3) * t262 + g(3) - t1) * t215 + (qJD(3) * t27 + t10 * t204 + t11 * t203 + t388) * t213) * m(5);];
tau = t2;
