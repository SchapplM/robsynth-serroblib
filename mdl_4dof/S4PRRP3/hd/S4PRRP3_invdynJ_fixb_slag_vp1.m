% Calculate vector of inverse dynamics joint torques for
% S4PRRP3
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:55
% DurationCPUTime: 9.72s
% Computational Cost: add. (5888->397), mult. (6858->500), div. (0->0), fcn. (5381->4), ass. (0->223)
t430 = Icges(4,3) + Icges(5,3);
t196 = pkin(6) + qJ(2);
t190 = sin(t196);
t199 = cos(qJ(3));
t292 = t190 * t199;
t198 = sin(qJ(3));
t293 = t190 * t198;
t191 = cos(t196);
t302 = Icges(5,6) * t191;
t91 = Icges(5,4) * t292 - Icges(5,2) * t293 - t302;
t303 = Icges(4,6) * t191;
t93 = Icges(4,4) * t292 - Icges(4,2) * t293 - t303;
t422 = t91 + t93;
t146 = Icges(5,5) * t199 - Icges(5,6) * t198;
t148 = Icges(4,5) * t199 - Icges(4,6) * t198;
t401 = t146 + t148;
t165 = Icges(5,4) * t293;
t306 = Icges(5,5) * t191;
t95 = Icges(5,1) * t292 - t165 - t306;
t166 = Icges(4,4) * t293;
t307 = Icges(4,5) * t191;
t97 = Icges(4,1) * t292 - t166 - t307;
t427 = t95 + t97;
t308 = Icges(5,4) * t198;
t154 = Icges(5,1) * t199 - t308;
t96 = Icges(5,5) * t190 + t154 * t191;
t309 = Icges(4,4) * t198;
t156 = Icges(4,1) * t199 - t309;
t98 = Icges(4,5) * t190 + t156 * t191;
t420 = t96 + t98;
t149 = Icges(5,2) * t199 + t308;
t151 = Icges(4,2) * t199 + t309;
t426 = t149 + t151;
t192 = Icges(5,4) * t199;
t153 = Icges(5,1) * t198 + t192;
t193 = Icges(4,4) * t199;
t155 = Icges(4,1) * t198 + t193;
t429 = t153 + t155;
t428 = t430 * t191;
t404 = t428 + (Icges(4,6) + Icges(5,6)) * t293 + (-Icges(4,5) - Icges(5,5)) * t292;
t413 = t430 * t190 + t401 * t191;
t227 = -Icges(5,2) * t198 + t192;
t92 = Icges(5,6) * t190 + t191 * t227;
t228 = -Icges(4,2) * t198 + t193;
t94 = Icges(4,6) * t190 + t191 * t228;
t421 = t92 + t94;
t424 = t422 * t198;
t391 = -t427 * t199 + t424;
t375 = -t391 * t190 + t404 * t191;
t423 = t420 * t292;
t145 = Icges(5,5) * t198 + Icges(5,6) * t199;
t147 = Icges(4,5) * t198 + Icges(4,6) * t199;
t418 = t145 + t147;
t417 = t154 + t156;
t416 = t426 * qJD(3);
t415 = t429 * qJD(3);
t414 = t227 + t228;
t399 = t426 * t198 - t199 * t429;
t412 = t413 * t191 - t423;
t411 = t421 * t198;
t290 = t191 * t199;
t359 = -t413 * t190 - t420 * t290;
t409 = t404 * t190 - t427 * t290;
t374 = -t421 * t293 - t412;
t291 = t191 * t198;
t373 = -t422 * t291 - t409;
t372 = -t421 * t291 - t359;
t371 = t427 * t198 + t422 * t199;
t370 = t420 * t198 + t421 * t199;
t408 = -t416 * t191 + (-t414 * t190 + t302 + t303) * qJD(2);
t407 = t421 * qJD(2) - t416 * t190;
t406 = -t415 * t191 + (-t417 * t190 + t306 + t307) * qJD(2);
t405 = -t420 * qJD(2) + t415 * t190;
t403 = t414 * qJD(3);
t402 = t417 * qJD(3);
t400 = t418 * qJD(3);
t398 = -t198 * t429 - t426 * t199;
t397 = t420 * t199 - t411;
t294 = t147 * t191;
t296 = t145 * t191;
t369 = -t399 * t190 - t294 - t296;
t295 = t147 * t190;
t297 = t145 * t190;
t368 = -t399 * t191 + t295 + t297;
t396 = t413 * qJD(2);
t395 = t191 ^ 2;
t394 = t418 * qJD(2) + t398 * qJD(3) - t403 * t198 + t402 * t199;
t393 = -qJD(3) * t370 - t198 * t408 + t199 * t406 + t396;
t392 = qJD(2) * t404 + qJD(3) * t371 + t198 * t407 + t199 * t405;
t390 = t372 * t190 - t191 * t373;
t389 = t374 * t190 - t375 * t191;
t388 = t399 * qJD(2) + qJD(3) * t401;
t387 = qJD(2) * t391 - t190 * t400 + t396;
t386 = -t400 * t191 + (-t190 * t401 - t397 + t428) * qJD(2);
t385 = t368 * qJD(2);
t384 = t369 * qJD(2);
t195 = t199 * pkin(3);
t194 = t199 * rSges(5,1);
t363 = -rSges(5,2) * t198 + t194;
t365 = t363 + t195;
t130 = t191 * pkin(2) + t190 * pkin(5);
t122 = t130 * qJD(2);
t268 = qJD(2) * qJD(3);
t124 = -qJDD(3) * t191 + t190 * t268;
t137 = t363 * qJD(3);
t324 = t199 * rSges(5,2);
t157 = rSges(5,1) * t198 + t324;
t176 = qJD(4) * t191;
t188 = t191 * pkin(5);
t129 = pkin(2) * t190 - t188;
t197 = -qJ(4) - pkin(5);
t171 = t191 * t197;
t343 = t195 + pkin(2);
t367 = -rSges(5,1) * t292 + rSges(5,2) * t293 + t191 * rSges(5,3) - t190 * t343 - t171;
t336 = t129 + t367;
t266 = -t129 + t336;
t270 = qJD(3) * t191;
t289 = t199 * qJD(3) ^ 2;
t115 = t157 * t190;
t183 = t190 * rSges(5,3);
t267 = pkin(3) * t293;
t281 = qJD(3) * t267 + t176;
t341 = -qJD(3) * t115 - t281 + ((-pkin(5) - t197) * t190 + t183 + t365 * t191) * qJD(2);
t10 = -t137 * t270 + qJDD(4) * t190 + t124 * t157 + (t124 * t198 - t191 * t289) * pkin(3) + t266 * qJDD(2) + (-t122 + t176 - t341) * qJD(2);
t383 = -g(1) + t10;
t123 = qJDD(3) * t190 + t191 * t268;
t175 = qJD(4) * t190;
t271 = qJD(3) * t190;
t272 = qJD(2) * t191;
t174 = pkin(5) * t272;
t273 = qJD(2) * t190;
t285 = qJD(2) * (-pkin(2) * t273 + t174) + qJDD(2) * t130;
t366 = rSges(5,1) * t290 - rSges(5,2) * t291 - t190 * t197 + t191 * t343 + t183;
t320 = -t130 + t366;
t263 = t198 * t270;
t214 = -t199 * t273 - t263;
t269 = qJD(3) * t199;
t262 = t191 * t269;
t265 = t198 * t273;
t364 = rSges(5,2) * t265 + rSges(5,3) * t272 + t175;
t342 = -pkin(3) * t263 - t174 + (-t190 * t195 - t171) * qJD(2) + rSges(5,1) * t214 - rSges(5,2) * t262 + t364;
t11 = -t137 * t271 - qJDD(4) * t191 - t123 * t157 + t320 * qJDD(2) + (-t123 * t198 - t190 * t289) * pkin(3) + (t175 + t342) * qJD(2) + t285;
t382 = g(2) - t11;
t381 = t389 * qJD(3) + t384;
t380 = t390 * qJD(3) + t385;
t379 = t391 * qJD(3) + t198 * t405 - t199 * t407;
t378 = qJD(3) * t397 + t198 * t406 + t199 * t408;
t377 = t388 * t190 + t394 * t191;
t376 = t394 * t190 - t388 * t191;
t362 = t387 * t395 + (t393 * t190 + (-t386 + t392) * t191) * t190;
t361 = t392 * t395 + (t386 * t190 + (-t387 + t393) * t191) * t190;
t360 = t404 + t411;
t313 = -t155 * t190 - t93;
t317 = -Icges(4,2) * t292 - t166 + t97;
t352 = -t198 * t317 + t199 * t313;
t315 = -t153 * t190 - t91;
t319 = -Icges(5,2) * t292 - t165 + t95;
t351 = -t198 * t319 + t199 * t315;
t350 = m(2) + m(3);
t349 = t123 / 0.2e1;
t348 = t124 / 0.2e1;
t335 = rSges(4,1) * t199;
t158 = rSges(4,1) * t198 + rSges(4,2) * t199;
t118 = t158 * t191;
t184 = t190 * rSges(4,3);
t102 = rSges(4,1) * t290 - rSges(4,2) * t291 + t184;
t74 = t102 + t130;
t42 = qJD(2) * t74 - t158 * t271;
t333 = t118 * t42;
t264 = t158 * t270;
t276 = rSges(4,2) * t293 + t191 * rSges(4,3);
t100 = rSges(4,1) * t292 - t276;
t286 = -t100 - t129;
t41 = qJD(2) * t286 - t264;
t332 = t190 * t41;
t331 = t191 * t41;
t254 = -pkin(3) * t198 - t157;
t213 = t254 * t270 + t175;
t30 = qJD(2) * t266 + t213;
t322 = t30 * t157;
t318 = -t149 * t191 + t96;
t316 = -t151 * t191 + t98;
t314 = -t153 * t191 - t92;
t312 = -t155 * t191 - t94;
t282 = rSges(4,2) * t265 + rSges(4,3) * t272;
t280 = -t149 + t154;
t279 = t153 + t227;
t278 = -t151 + t156;
t277 = t155 + t228;
t275 = qJD(2) * t146;
t274 = qJD(2) * t148;
t261 = -pkin(2) - t335;
t258 = -t271 / 0.2e1;
t257 = t271 / 0.2e1;
t256 = -t270 / 0.2e1;
t255 = t270 / 0.2e1;
t247 = -pkin(3) * t291 - t157 * t191;
t246 = -pkin(3) * t269 - t137;
t128 = rSges(3,1) * t191 - rSges(3,2) * t190;
t127 = rSges(3,1) * t190 + rSges(3,2) * t191;
t160 = -rSges(4,2) * t198 + t335;
t236 = -t190 * t42 - t331;
t68 = rSges(4,1) * t214 - rSges(4,2) * t262 + t282;
t116 = t158 * t190;
t70 = -qJD(3) * t116 + (t160 * t191 + t184) * qJD(2);
t235 = t190 * t70 + t191 * t68;
t226 = t100 * t190 + t102 * t191;
t221 = -t324 + (-rSges(5,1) - pkin(3)) * t198;
t212 = t191 * t221;
t211 = -t198 * t318 + t199 * t314;
t210 = -t198 * t316 + t199 * t312;
t209 = (-t198 * t279 + t199 * t280) * qJD(2);
t208 = (-t198 * t277 + t199 * t278) * qJD(2);
t138 = t160 * qJD(3);
t126 = qJD(2) * t129;
t40 = qJD(3) * t226 + qJD(1);
t31 = -t157 * t271 + (t130 + t320) * qJD(2) - t281;
t27 = qJD(1) + (-t190 * t336 + t191 * t320) * qJD(3);
t22 = qJD(2) * t68 + qJDD(2) * t102 - t123 * t158 - t138 * t271 + t285;
t21 = -t138 * t270 + t124 * t158 + t286 * qJDD(2) + (-t122 - t70) * qJD(2);
t16 = qJD(3) * t235 + t100 * t123 - t102 * t124 + qJDD(1);
t1 = qJDD(1) - t320 * t124 - t336 * t123 + (t190 * t341 + t191 * t342) * qJD(3);
t2 = [t350 * qJDD(1) + m(4) * t16 + m(5) * t1 + (-m(4) - m(5) - t350) * g(3); -m(3) * (-g(1) * t127 + g(2) * t128) + ((((t413 + t424) * t191 + t374 + t409 - t423) * t191 - t359 * t190) * qJD(3) + t385) * t255 + (-t399 * qJD(3) + t402 * t198 + t403 * t199) * qJD(2) + (t30 * t281 + t31 * t364 + (t190 * t322 + t212 * t31) * qJD(3) + ((t30 * (-t363 - t343) - t31 * t197) * t191 + (t30 * (t197 - rSges(5,3)) + t31 * (-t194 - t343)) * t190) * qJD(2) - (qJD(2) * t336 - t126 + t213 - t30) * t31 - t382 * t366 + t383 * t367) * m(5) + (t42 * (t174 + t282) + (t158 * t332 - t333) * qJD(3) + ((-pkin(2) - t160) * t331 + (t41 * (-rSges(4,3) - pkin(5)) + t42 * t261) * t190) * qJD(2) - (-qJD(2) * t100 - t126 - t264 - t41) * t42 + (-g(2) + t22) * t74 + (-g(1) + t21) * (t261 * t190 + t188 + t276)) * m(4) + (t368 + t370) * t349 + (t369 + t371) * t348 + (((t191 * t360 + t359 + t372) * t191 + (t190 * t360 + t373 + t412) * t190) * qJD(3) + t381 - t384) * t258 + (t377 + t378) * t257 + (m(3) * (t127 ^ 2 + t128 ^ 2) + Icges(3,3) - t398) * qJDD(2) + (t376 - t379 + t380) * t256; t390 * t349 + t389 * t348 + (t377 * qJD(2) + t361 * qJD(3) + t368 * qJDD(2) + t372 * t123 + t373 * t124) * t190 / 0.2e1 - (t376 * qJD(2) + t362 * qJD(3) + t369 * qJDD(2) + t374 * t123 + t375 * t124) * t191 / 0.2e1 - ((((-t317 - t319) * t191 + (t316 + t318) * t190) * t199 + ((-t313 - t315) * t191 + (t312 + t314) * t190) * t198) * qJD(3) + ((t277 + t279) * t199 + (t278 + t280) * t198) * qJD(2)) * qJD(2) / 0.2e1 + (t379 * t191 + t378 * t190 + (t371 * t190 + t370 * t191) * qJD(2)) * qJD(2) / 0.2e1 + (t370 * t190 - t371 * t191) * qJDD(2) / 0.2e1 + t381 * t273 / 0.2e1 + t380 * t272 / 0.2e1 + ((-t271 * t296 + t275) * t190 + (t209 + (-t351 * t191 + (t297 + t211) * t190) * qJD(3)) * t191 + (-t271 * t294 + t274) * t190 + (t208 + (-t352 * t191 + (t295 + t210) * t190) * qJD(3)) * t191) * t258 + ((t373 * t190 + t191 * t372) * qJD(2) + t361) * t257 + ((t190 * t375 + t374 * t191) * qJD(2) + t362) * t256 + ((-t270 * t297 - t275) * t191 + (t209 + (t211 * t190 + (t296 - t351) * t191) * qJD(3)) * t190 + (-t270 * t295 - t274) * t191 + (t208 + (t210 * t190 + (t294 - t352) * t191) * qJD(3)) * t190) * t255 + (g(1) * t118 + g(2) * t116 - g(3) * t160 - (t116 * t41 - t333) * qJD(2) - (t40 * (-t116 * t190 - t118 * t191) + t236 * t160) * qJD(3) + t16 * t226 + t40 * ((t100 * t191 - t102 * t190) * qJD(2) + t235) + t236 * t138 + (-t22 * t190 - t21 * t191 + (-t191 * t42 + t332) * qJD(2)) * t158) * m(4) + (-g(3) * t365 - g(1) * t212 - (t30 * t115 + t247 * t31) * qJD(2) + (t10 * t254 + t30 * t246 + t1 * t320 + t27 * t342 + (t254 * t31 - t27 * t336) * qJD(2) - (t247 * t27 - t30 * t365) * qJD(3)) * t191 + (-g(2) * t221 + t11 * t254 + t31 * t246 - t1 * t336 + t27 * t341 + (-t27 * t320 + t322) * qJD(2) - (-t31 * t365 + (-t115 - t267) * t27) * qJD(3)) * t190) * m(5); (t190 * t383 + t382 * t191) * m(5);];
tau = t2;
