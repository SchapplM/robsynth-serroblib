% Calculate vector of inverse dynamics joint torques for
% S4PRRP4
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:45
% EndTime: 2019-12-31 16:27:57
% DurationCPUTime: 10.50s
% Computational Cost: add. (5906->419), mult. (7224->522), div. (0->0), fcn. (5677->4), ass. (0->227)
t203 = pkin(6) + qJ(2);
t195 = sin(t203);
t196 = cos(t203);
t204 = sin(qJ(3));
t198 = Icges(5,5) * t204;
t205 = cos(qJ(3));
t235 = Icges(5,1) * t205 + t198;
t94 = -Icges(5,4) * t196 + t195 * t235;
t301 = t195 * t204;
t171 = Icges(4,4) * t301;
t300 = t195 * t205;
t316 = Icges(4,5) * t196;
t96 = Icges(4,1) * t300 - t171 - t316;
t425 = t94 + t96;
t318 = Icges(4,4) * t204;
t152 = Icges(4,2) * t205 + t318;
t233 = Icges(5,3) * t205 - t198;
t427 = t152 + t233;
t315 = Icges(5,5) * t205;
t154 = Icges(5,1) * t204 - t315;
t199 = Icges(4,4) * t205;
t156 = Icges(4,1) * t204 + t199;
t423 = t154 + t156;
t149 = Icges(4,5) * t205 - Icges(4,6) * t204;
t89 = Icges(4,3) * t195 + t149 * t196;
t151 = Icges(5,4) * t205 + Icges(5,6) * t204;
t91 = Icges(5,2) * t195 + t151 * t196;
t426 = t89 + t91;
t95 = Icges(5,4) * t195 + t196 * t235;
t157 = Icges(4,1) * t205 - t318;
t97 = Icges(4,5) * t195 + t157 * t196;
t421 = t95 + t97;
t312 = Icges(4,6) * t196;
t92 = Icges(4,4) * t300 - Icges(4,2) * t301 - t312;
t335 = t204 * t92;
t147 = Icges(5,3) * t204 + t315;
t86 = -Icges(5,6) * t196 + t147 * t195;
t399 = -t204 * t86 - t425 * t205 + t335;
t309 = Icges(4,3) * t196;
t88 = Icges(4,5) * t300 - Icges(4,6) * t301 - t309;
t422 = -t399 * t195 - t196 * t88;
t148 = Icges(4,5) * t204 + Icges(4,6) * t205;
t150 = Icges(5,4) * t204 - Icges(5,6) * t205;
t419 = t148 + t150;
t418 = t427 * qJD(3);
t417 = t423 * qJD(3);
t298 = t196 * t205;
t299 = t196 * t204;
t170 = Icges(5,5) * t298;
t311 = Icges(5,6) * t195;
t87 = Icges(5,3) * t299 + t170 + t311;
t416 = t426 * t195 + t421 * t298 + t87 * t299;
t90 = -Icges(5,2) * t196 + t151 * t195;
t84 = t195 * t90;
t415 = -t195 * t88 - t425 * t298 - t86 * t299 - t84;
t407 = -t204 * t427 + t423 * t205;
t336 = t196 * t90;
t384 = -t336 + t422;
t383 = -t299 * t92 - t415;
t234 = -Icges(4,2) * t204 + t199;
t93 = Icges(4,6) * t195 + t196 * t234;
t382 = -t299 * t93 + t416;
t381 = (-t86 + t92) * t205 + t425 * t204;
t380 = (-t87 + t93) * t205 + t421 * t204;
t414 = t418 * t196 + (t195 * t234 - t312 - t86) * qJD(2);
t413 = t418 * t195 + (t147 * t196 + t311 - t93) * qJD(2);
t412 = -t417 * t196 + (-t157 * t195 + t316 - t94) * qJD(2);
t411 = -t421 * qJD(2) + t417 * t195;
t256 = t196 * t91 - t95 * t300 - t87 * t301;
t75 = t97 * t300;
t259 = t196 * t89 - t75;
t33 = -t301 * t93 - t259;
t377 = -t256 + t33;
t410 = (t147 - t234) * qJD(3);
t409 = (t157 + t235) * qJD(3);
t408 = t419 * qJD(3);
t406 = -t423 * t204 - t205 * t427;
t334 = t204 * t93;
t405 = t204 * t87 + t421 * t205 - t334;
t302 = t150 * t196;
t304 = t148 * t196;
t379 = t407 * t195 - t302 - t304;
t303 = t150 * t195;
t305 = t148 * t195;
t378 = t407 * t196 + t303 + t305;
t404 = t426 * qJD(2);
t403 = t196 ^ 2;
t402 = t419 * qJD(2) + t406 * qJD(3) + t410 * t204 + t409 * t205;
t401 = -t380 * qJD(3) + t414 * t204 + t412 * t205 + t404;
t373 = qJD(2) * t90;
t400 = -qJD(2) * t88 + t381 * qJD(3) - t413 * t204 + t411 * t205 - t373;
t398 = t382 * t195 - t383 * t196;
t397 = t377 * t195 - t384 * t196;
t396 = (-t149 - t151) * qJD(3) + t407 * qJD(2);
t395 = t399 * qJD(2) - t408 * t195 + t404;
t394 = -t373 - t408 * t196 + (-t149 * t195 + t309 - t405) * qJD(2);
t393 = t378 * qJD(2);
t392 = rSges(5,3) + qJ(4);
t391 = t379 * qJD(2);
t352 = rSges(5,1) + pkin(3);
t390 = t397 * qJD(3) + t391;
t389 = t398 * qJD(3) + t393;
t388 = t399 * qJD(3) + t411 * t204 + t413 * t205;
t387 = t405 * qJD(3) + t412 * t204 - t414 * t205;
t386 = -t396 * t195 + t402 * t196;
t385 = t402 * t195 + t396 * t196;
t197 = t204 * qJ(4);
t250 = t205 * pkin(3) + t197;
t251 = t205 * rSges(5,1) + t204 * rSges(5,3);
t376 = t250 + t251;
t375 = t395 * t403 + (t401 * t195 + (-t394 + t400) * t196) * t195;
t374 = t400 * t403 + (t394 * t195 + (-t395 + t401) * t196) * t195;
t372 = t204 * t352;
t134 = t196 * pkin(2) + t195 * pkin(5);
t189 = t195 * rSges(5,2);
t295 = rSges(5,3) * t299 + t196 * t197 + t352 * t298 + t189;
t371 = t134 + t295;
t370 = t392 * t300;
t369 = t392 * t298;
t368 = t336 + t416;
t274 = qJD(4) * t204;
t158 = t196 * t274;
t275 = qJD(3) * t196;
t267 = t205 * t275;
t277 = qJD(2) * t196;
t367 = rSges(5,2) * t277 + t392 * t267 + t158;
t366 = -g(1) * t196 - g(2) * t195;
t324 = -t156 * t195 - t92;
t328 = -Icges(4,2) * t300 - t171 + t96;
t359 = -t204 * t328 + t205 * t324;
t326 = -t154 * t195 + t86;
t330 = t233 * t195 - t94;
t358 = t204 * t330 + t205 * t326;
t357 = m(2) + m(3);
t272 = qJD(2) * qJD(3);
t127 = qJDD(3) * t195 + t196 * t272;
t356 = t127 / 0.2e1;
t128 = -qJDD(3) * t196 + t195 * t272;
t355 = t128 / 0.2e1;
t268 = t204 * t275;
t278 = qJD(2) * t195;
t219 = -t205 * t278 - t268;
t270 = t204 * t278;
t347 = t352 * t219 - t270 * t392 + t367;
t159 = pkin(3) * t204 - qJ(4) * t205;
t160 = rSges(5,1) * t204 - rSges(5,3) * t205;
t276 = qJD(3) * t195;
t346 = t250 * t277 + (-qJD(3) * t159 + t274) * t195 - t160 * t276 + (t196 * t251 + t189) * qJD(2);
t343 = rSges(4,1) * t205;
t161 = rSges(4,1) * t204 + rSges(4,2) * t205;
t120 = t161 * t196;
t188 = t195 * rSges(4,3);
t101 = rSges(4,1) * t298 - rSges(4,2) * t299 + t188;
t72 = t101 + t134;
t40 = qJD(2) * t72 - t161 * t276;
t342 = t120 * t40;
t269 = t161 * t275;
t193 = t196 * pkin(5);
t133 = pkin(2) * t195 - t193;
t281 = rSges(4,2) * t301 + t196 * rSges(4,3);
t99 = rSges(4,1) * t300 - t281;
t321 = -t133 - t99;
t39 = qJD(2) * t321 - t269;
t340 = t195 * t39;
t283 = t159 + t160;
t257 = t283 * qJD(3);
t228 = -t196 * t257 + t158;
t191 = t196 * rSges(5,2);
t322 = t376 * t195 - t191;
t271 = -t133 - t322;
t28 = qJD(2) * t271 + t228;
t339 = t196 * t28;
t338 = t196 * t39;
t273 = qJD(4) * t205;
t27 = -t273 + qJD(1) + (t195 * t322 + t196 * t295) * qJD(3);
t333 = t27 * t204;
t329 = t233 * t196 - t95;
t327 = -t152 * t196 + t97;
t325 = -Icges(5,1) * t299 + t170 + t87;
t323 = -t156 * t196 - t93;
t294 = -t352 * t301 + t370;
t293 = -t352 * t299 + t369;
t181 = pkin(5) * t277;
t292 = qJD(2) * (-pkin(2) * t278 + t181) + qJDD(2) * t134;
t291 = -t376 * qJD(3) + t273;
t288 = rSges(4,2) * t270 + rSges(4,3) * t277;
t287 = -t233 + t235;
t286 = t147 - t154;
t285 = -t152 + t157;
t284 = t156 + t234;
t280 = qJD(2) * t149;
t279 = qJD(2) * t151;
t266 = -pkin(2) - t343;
t263 = -t276 / 0.2e1;
t262 = t276 / 0.2e1;
t261 = -t275 / 0.2e1;
t260 = t275 / 0.2e1;
t258 = -t88 + t334;
t132 = rSges(3,1) * t196 - rSges(3,2) * t195;
t131 = rSges(3,1) * t195 + rSges(3,2) * t196;
t164 = -rSges(4,2) * t204 + t343;
t29 = (-t257 + t274) * t195 + t371 * qJD(2);
t249 = t195 * t29 + t339;
t244 = -t195 * t40 - t338;
t68 = rSges(4,1) * t219 - rSges(4,2) * t267 + t288;
t116 = t161 * t195;
t70 = -qJD(3) * t116 + (t164 * t196 + t188) * qJD(2);
t243 = t195 * t70 + t196 * t68;
t238 = t101 * t196 + t195 * t99;
t227 = t249 * t205;
t218 = t204 * t329 + t205 * t325;
t217 = -t204 * t327 + t205 * t323;
t216 = -t204 * t392 - t205 * t352 - pkin(2);
t215 = qJDD(4) * t204 + (t273 + t291) * qJD(3);
t214 = (t204 * t286 + t205 * t287) * qJD(2);
t213 = (-t204 * t284 + t205 * t285) * qJD(2);
t142 = t164 * qJD(3);
t130 = qJD(2) * t133;
t126 = t134 * qJD(2);
t38 = qJD(3) * t238 + qJD(1);
t22 = qJD(2) * t68 + qJDD(2) * t101 - t127 * t161 - t142 * t276 + t292;
t21 = -t142 * t275 + t128 * t161 + t321 * qJDD(2) + (-t126 - t70) * qJD(2);
t16 = qJD(3) * t243 - t101 * t128 + t127 * t99 + qJDD(1);
t11 = -t283 * t127 + t295 * qJDD(2) + (t158 + t347) * qJD(2) + t215 * t195 + t292;
t10 = t283 * t128 + t271 * qJDD(2) + t215 * t196 + (-t195 * t274 - t126 - t346) * qJD(2);
t1 = -qJDD(4) * t205 + qJDD(1) - t295 * t128 + t322 * t127 + (t195 * t346 + t196 * t347 + t274) * qJD(3);
t2 = [t357 * qJDD(1) + m(4) * t16 + m(5) * t1 + (-m(4) - m(5) - t357) * g(3); -m(3) * (-g(1) * t131 + g(2) * t132) + (((t33 - t75 + (t89 + t335) * t196 + t415) * t196 + (t368 + t384 - t422) * t195) * qJD(3) + t393) * t260 + (t407 * qJD(3) + t409 * t204 - t410 * t205) * qJD(2) + (t29 * (-t352 * t268 + t181 + t367) + (-t274 + (-t205 * t392 + t372) * qJD(3)) * t195 * t28 + (t216 * t339 + (t28 * (-rSges(5,2) - pkin(5)) + t29 * (-pkin(2) - t376)) * t195) * qJD(2) - (-qJD(2) * t322 - t130 + t228 - t28) * t29 + (-g(2) + t11) * t371 + (-g(1) + t10) * (t195 * t216 + t191 + t193)) * m(5) + (t40 * (t181 + t288) + (t161 * t340 - t342) * qJD(3) + ((-pkin(2) - t164) * t338 + (t39 * (-rSges(4,3) - pkin(5)) + t40 * t266) * t195) * qJD(2) - (-qJD(2) * t99 - t130 - t269 - t39) * t40 + (-g(2) + t22) * t72 + (-g(1) + t21) * (t266 * t195 + t193 + t281)) * m(4) + (t378 + t380) * t356 + (t379 + t381) * t355 + (((t196 * t258 - t368 + t382) * t196 + (t195 * t258 + t256 + t259 + t383 - t84) * t195) * qJD(3) + t390 - t391) * t263 + (t386 + t387) * t262 + (m(3) * (t131 ^ 2 + t132 ^ 2) + Icges(3,3) - t406) * qJDD(2) + (t385 - t388 + t389) * t261; t398 * t356 + t397 * t355 + (t386 * qJD(2) + t374 * qJD(3) + t378 * qJDD(2) + t382 * t127 + t383 * t128) * t195 / 0.2e1 - (t385 * qJD(2) + t375 * qJD(3) + t379 * qJDD(2) + t377 * t127 + t384 * t128) * t196 / 0.2e1 - ((((-t328 + t330) * t196 + (t327 - t329) * t195) * t205 + ((-t324 - t326) * t196 + (t323 + t325) * t195) * t204) * qJD(3) + ((t284 - t286) * t205 + (t285 + t287) * t204) * qJD(2)) * qJD(2) / 0.2e1 + (t388 * t196 + t387 * t195 + (t381 * t195 + t380 * t196) * qJD(2)) * qJD(2) / 0.2e1 + (t380 * t195 - t381 * t196) * qJDD(2) / 0.2e1 + t390 * t278 / 0.2e1 + t389 * t277 / 0.2e1 + ((-t276 * t302 + t279) * t195 + (t214 + (-t358 * t196 + (t303 + t218) * t195) * qJD(3)) * t196 + (-t276 * t304 + t280) * t195 + (t213 + (-t359 * t196 + (t305 + t217) * t195) * qJD(3)) * t196) * t263 + ((t383 * t195 + t382 * t196) * qJD(2) + t374) * t262 + ((t384 * t195 + t377 * t196) * qJD(2) + t375) * t261 + ((-t275 * t303 - t279) * t196 + (t214 + (t218 * t195 + (t302 - t358) * t196) * qJD(3)) * t195 + (-t275 * t305 - t280) * t196 + (t213 + (t217 * t195 + (t304 - t359) * t196) * qJD(3)) * t195) * t260 + (-g(1) * t369 - g(2) * t370 - g(3) * t376 - t366 * t372 + (-t10 * t283 + t28 * t291 + t1 * t295 + t27 * t347 + (t27 * t322 - t283 * t29) * qJD(2)) * t196 + (-t11 * t283 + t29 * t291 + t1 * t322 + t27 * t346 + (-t27 * t295 + t28 * t283) * qJD(2)) * t195 - (t227 + t333) * qJD(4) - (-t28 * t294 + t29 * t293) * qJD(2) - ((t27 * t293 - t28 * t376) * t196 + (t27 * t294 - t29 * t376) * t195) * qJD(3)) * m(5) + (t16 * t238 + t38 * ((-t101 * t195 + t196 * t99) * qJD(2) + t243) + t244 * t142 + (-t22 * t195 - t21 * t196 + (-t196 * t40 + t340) * qJD(2)) * t161 - (t116 * t39 - t342) * qJD(2) - (t38 * (-t116 * t195 - t120 * t196) + t244 * t164) * qJD(3) + g(1) * t120 + g(2) * t116 - g(3) * t164) * m(4); (-(t227 + (t195 ^ 2 + t403) * t333) * qJD(3) + (qJD(3) * t249 + g(3) - t1) * t205 + (qJD(3) * t27 + t10 * t196 + t11 * t195 + t366) * t204) * m(5);];
tau = t2;
