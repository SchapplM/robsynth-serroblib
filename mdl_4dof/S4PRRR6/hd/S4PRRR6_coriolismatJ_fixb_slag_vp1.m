% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:41
% EndTime: 2019-12-31 16:34:53
% DurationCPUTime: 10.32s
% Computational Cost: add. (32409->595), mult. (49093->907), div. (0->0), fcn. (53863->8), ass. (0->349)
t322 = sin(qJ(2));
t317 = t322 ^ 2;
t324 = cos(qJ(2));
t481 = t324 ^ 2;
t492 = t322 * t324;
t323 = cos(qJ(3));
t458 = pkin(3) * t323;
t254 = -pkin(6) * t324 + t322 * t458;
t319 = sin(pkin(7));
t229 = t254 * t319;
t490 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t492 + (0.2e1 * t481 - 0.2e1 * t317) * Icges(3,4);
t465 = t319 / 0.2e1;
t320 = cos(pkin(7));
t464 = -t320 / 0.2e1;
t462 = t322 / 0.2e1;
t461 = -t324 / 0.2e1;
t362 = -Icges(3,5) * t322 - Icges(3,6) * t324;
t298 = t362 * t319;
t299 = t362 * t320;
t479 = t320 ^ 2;
t480 = t319 ^ 2;
t486 = t479 + t480;
t318 = qJ(3) + qJ(4);
t315 = sin(t318);
t316 = cos(t318);
t417 = t319 * t324;
t270 = -t315 * t417 - t316 * t320;
t271 = -t315 * t320 + t316 * t417;
t418 = t319 * t322;
t191 = rSges(5,1) * t271 + rSges(5,2) * t270 + rSges(5,3) * t418;
t415 = t320 * t322;
t181 = t191 * t415;
t321 = sin(qJ(3));
t416 = t320 * t321;
t255 = pkin(6) * t322 + t324 * t458;
t425 = t255 * t319;
t213 = -pkin(3) * t416 + t425;
t414 = t320 * t324;
t272 = -t315 * t414 + t316 * t319;
t273 = t315 * t319 + t316 * t414;
t192 = rSges(5,1) * t273 + rSges(5,2) * t272 + rSges(5,3) * t415;
t419 = t319 * t321;
t384 = pkin(3) * t419;
t214 = t255 * t320 + t384;
t403 = t192 + t214;
t119 = t181 + (t213 * t320 - t319 * t403) * t322;
t184 = t324 * t191;
t372 = rSges(5,1) * t316 - rSges(5,2) * t315;
t262 = -rSges(5,3) * t324 + t322 * t372;
t167 = t262 * t418 + t184;
t426 = t213 * t324;
t135 = t254 * t418 + t167 + t426;
t375 = t403 * t324;
t396 = -t254 - t262;
t136 = t396 * t415 - t375;
t209 = rSges(5,1) * t270 - rSges(5,2) * t271;
t193 = t209 * t415;
t210 = rSges(5,1) * t272 - rSges(5,2) * t273;
t158 = -t210 * t418 + t193;
t289 = (-rSges(5,1) * t315 - rSges(5,2) * t316) * t322;
t173 = t209 * t324 + t289 * t418;
t174 = -t210 * t324 - t289 * t415;
t197 = Icges(5,5) * t270 - Icges(5,6) * t271;
t264 = Icges(5,4) * t270;
t189 = Icges(5,1) * t271 + Icges(5,5) * t418 + t264;
t406 = -Icges(5,2) * t271 + t189 + t264;
t443 = Icges(5,4) * t271;
t187 = Icges(5,2) * t270 + Icges(5,6) * t418 + t443;
t408 = Icges(5,1) * t270 - t187 - t443;
t105 = -t197 * t324 + (-t315 * t406 + t316 * t408) * t322;
t198 = Icges(5,5) * t272 - Icges(5,6) * t273;
t265 = Icges(5,4) * t272;
t190 = Icges(5,1) * t273 + Icges(5,5) * t415 + t265;
t405 = -Icges(5,2) * t273 + t190 + t265;
t442 = Icges(5,4) * t273;
t188 = Icges(5,2) * t272 + Icges(5,6) * t415 + t442;
t407 = Icges(5,1) * t272 - t188 - t442;
t106 = -t198 * t324 + (-t315 * t405 + t316 * t407) * t322;
t282 = (-Icges(5,5) * t315 - Icges(5,6) * t316) * t322;
t441 = Icges(5,4) * t315;
t367 = Icges(5,1) * t316 - t441;
t260 = -Icges(5,5) * t324 + t322 * t367;
t393 = t260 + (-Icges(5,2) * t316 - t441) * t322;
t440 = Icges(5,4) * t316;
t363 = -Icges(5,2) * t315 + t440;
t258 = -Icges(5,6) * t324 + t322 * t363;
t394 = -t258 + (-Icges(5,1) * t315 - t440) * t322;
t422 = t282 * t324;
t89 = t197 * t418 + t270 * t406 + t271 * t408;
t90 = t198 * t418 + t270 * t405 + t271 * t407;
t36 = -(t270 * t393 + t271 * t394) * t324 + (t90 * t320 + (t89 - t422) * t319) * t322;
t91 = t197 * t415 + t272 * t406 + t273 * t408;
t92 = t198 * t415 + t272 * t405 + t273 * t407;
t37 = -(t272 * t393 + t273 * t394) * t324 + (t91 * t319 + (t92 - t422) * t320) * t322;
t377 = t415 / 0.2e1;
t379 = t418 / 0.2e1;
t385 = t36 * t379 + t37 * t377 + (t282 * t481 + (t106 * t320 + t105 * t319 - (-t315 * t393 + t316 * t394) * t324) * t322) * t461;
t7 = t385 + m(5) * (t119 * t158 + t135 * t173 + t136 * t174);
t487 = t7 * qJD(4);
t421 = (-Icges(4,5) * t321 - Icges(4,6) * t323) * t492;
t444 = Icges(4,4) * t323;
t364 = -Icges(4,2) * t321 + t444;
t276 = -Icges(4,6) * t324 + t322 * t364;
t392 = -t276 + (-Icges(4,1) * t321 - t444) * t322;
t445 = Icges(4,4) * t321;
t368 = Icges(4,1) * t323 - t445;
t278 = -Icges(4,5) * t324 + t322 * t368;
t391 = t278 + (-Icges(4,2) * t323 - t445) * t322;
t311 = pkin(2) * t324 + pkin(5) * t322;
t387 = t486 * t311;
t404 = t191 + t213;
t113 = t319 * t404 + t320 * t403 + t387;
t163 = t209 * t319 + t210 * t320;
t413 = t321 * t324;
t294 = -t319 * t413 - t320 * t323;
t292 = t294 * pkin(3);
t296 = t319 * t323 - t320 * t413;
t293 = t296 * pkin(3);
t147 = t292 * t319 + t293 * t320 + t163;
t412 = t323 * t324;
t295 = t319 * t412 - t416;
t211 = rSges(4,1) * t295 + rSges(4,2) * t294 + rSges(4,3) * t418;
t297 = t320 * t412 + t419;
t212 = rSges(4,1) * t297 + rSges(4,2) * t296 + rSges(4,3) * t415;
t148 = t211 * t319 + t212 * t320 + t387;
t227 = rSges(4,1) * t294 - rSges(4,2) * t295;
t228 = rSges(4,1) * t296 - rSges(4,2) * t297;
t169 = t227 * t319 + t228 * t320;
t310 = t322 * pkin(2) - t324 * pkin(5);
t381 = -t310 + t396;
t175 = t381 * t319;
t177 = t381 * t320;
t373 = rSges(4,1) * t323 - rSges(4,2) * t321;
t280 = -rSges(4,3) * t324 + t322 * t373;
t390 = -t280 - t310;
t223 = t390 * t319;
t225 = t390 * t320;
t459 = pkin(3) * t321;
t374 = t322 * t459 - t289;
t232 = t374 * t319;
t233 = t374 * t320;
t307 = (-rSges(4,1) * t321 - rSges(4,2) * t323) * t322;
t217 = Icges(4,5) * t294 - Icges(4,6) * t295;
t290 = Icges(4,4) * t294;
t207 = Icges(4,1) * t295 + Icges(4,5) * t418 + t290;
t400 = -Icges(4,2) * t295 + t207 + t290;
t447 = Icges(4,4) * t295;
t205 = Icges(4,2) * t294 + Icges(4,6) * t418 + t447;
t402 = Icges(4,1) * t294 - t205 - t447;
t107 = t217 * t418 + t294 * t400 + t295 * t402;
t218 = Icges(4,5) * t296 - Icges(4,6) * t297;
t291 = Icges(4,4) * t296;
t208 = Icges(4,1) * t297 + Icges(4,5) * t415 + t291;
t399 = -Icges(4,2) * t297 + t208 + t291;
t446 = Icges(4,4) * t297;
t206 = Icges(4,2) * t296 + Icges(4,6) * t415 + t446;
t401 = Icges(4,1) * t296 - t206 - t446;
t108 = t218 * t418 + t294 * t399 + t295 * t401;
t63 = -t107 * t320 + t108 * t319;
t109 = t217 * t415 + t296 * t400 + t297 * t402;
t110 = t218 * t415 + t296 * t399 + t297 * t401;
t64 = -t109 * t320 + t110 * t319;
t485 = t63 * t464 + t64 * t465 + m(5) * (t113 * t147 + t175 * t232 + t177 * t233) + m(4) * (t148 * t169 + (-t223 * t319 - t225 * t320) * t307);
t186 = Icges(5,5) * t273 + Icges(5,6) * t272 + Icges(5,3) * t415;
t237 = t258 * t320;
t239 = t260 * t320;
t360 = Icges(5,5) * t316 - Icges(5,6) * t315;
t256 = -Icges(5,3) * t324 + t322 * t360;
t356 = -t188 * t315 + t190 * t316;
t347 = -t256 * t320 - t356;
t100 = -t347 * t324 + (t237 * t315 - t239 * t316 + t186) * t322;
t353 = -t258 * t315 + t260 * t316;
t424 = t256 * t324;
t164 = t322 * t353 - t424;
t259 = Icges(5,6) * t322 + t324 * t363;
t261 = Icges(5,5) * t322 + t324 * t367;
t344 = (Icges(5,3) * t322 + t324 * t360 - t353) * t324;
t357 = -t187 * t315 + t189 * t316;
t185 = Icges(5,5) * t271 + Icges(5,6) * t270 + Icges(5,3) * t418;
t433 = t185 * t324;
t127 = t322 * t357 - t433;
t432 = t186 * t324;
t128 = t322 * t356 - t432;
t359 = t127 * t319 + t128 * t320;
t236 = t258 * t319;
t238 = t260 * t319;
t348 = -t256 * t319 - t357;
t99 = -t348 * t324 + (t236 * t315 - t238 * t316 + t185) * t322;
t484 = ((t344 + t359) * t324 + (t100 * t320 + t99 * t319 - (-t259 * t315 + t261 * t316 + t256) * t324 + t164) * t322) * t461 + (-t164 * t324 + t322 * t359) * t462;
t56 = t319 * t90 - t320 * t89;
t57 = t319 * t92 - t320 * t91;
t456 = t464 * t56 + t465 * t57;
t361 = Icges(4,5) * t323 - Icges(4,6) * t321;
t274 = -Icges(4,3) * t324 + t322 * t361;
t478 = 2 * qJD(2);
t477 = 2 * qJD(3);
t476 = 4 * qJD(3);
t475 = m(4) / 0.2e1;
t474 = m(5) / 0.2e1;
t251 = t280 * t319;
t252 = t280 * t320;
t427 = t212 * t324;
t428 = t211 * t324;
t134 = (-t251 * t322 + t428) * t320 + (t252 * t322 - t427) * t319;
t281 = rSges(4,3) * t322 + t324 * t373;
t351 = t280 * t324 + t281 * t322;
t149 = -t211 * t322 - t251 * t324 + t319 * t351;
t150 = t212 * t322 + t252 * t324 - t320 * t351;
t160 = (t211 * t320 - t212 * t319) * t322;
t171 = t280 * t418 + t428;
t172 = -t280 * t415 - t427;
t473 = m(4) * (t134 * t160 + t149 * t171 + t150 * t172);
t241 = t262 * t320;
t240 = t262 * t319;
t409 = t184 * t320 - t240 * t415;
t431 = t192 * t324;
t125 = (t241 * t322 - t431) * t319 + t409;
t263 = rSges(5,3) * t322 + t324 * t372;
t382 = -t240 * t324 + t262 * t417 + t263 * t418;
t139 = -t191 * t322 + t382;
t183 = t322 * t192;
t141 = t241 * t324 + t183 + (-t262 * t324 - t263 * t322) * t320;
t153 = -t183 * t319 + t181;
t168 = -t262 * t415 - t431;
t397 = -t254 * t320 - t241;
t82 = (-t229 * t322 + t426) * t320 + (-t322 * t397 - t375) * t319 + t409;
t97 = (-t404 + t425) * t322 + t382;
t395 = -t255 - t263;
t98 = t214 * t322 + t183 - t397 * t324 + (t322 * t395 + t324 * t396) * t320;
t471 = m(5) * (t119 * t125 + t135 * t139 + t136 * t141 + t153 * t82 + t167 * t97 + t168 * t98);
t470 = m(5) * (t119 * t82 + t135 * t97 + t136 * t98);
t383 = t113 * t158 + t173 * t177 + t174 * t175;
t469 = m(5) * (t119 * t163 + (-t135 * t320 - t136 * t319) * t289 + t383);
t468 = m(5) * (t147 * t153 + t167 * t233 + t168 * t232 + t383);
t463 = t320 / 0.2e1;
t460 = m(5) * t158;
t203 = Icges(4,5) * t295 + Icges(4,6) * t294 + Icges(4,3) * t418;
t430 = t203 * t324;
t204 = Icges(4,5) * t297 + Icges(4,6) * t296 + Icges(4,3) * t415;
t429 = t204 * t324;
t423 = t274 * t324;
t39 = 0.2e1 * (t82 / 0.4e1 - t147 / 0.4e1) * m(5) + 0.2e1 * (t134 / 0.4e1 - t169 / 0.4e1) * m(4);
t411 = t39 * qJD(1);
t95 = 0.2e1 * (t125 / 0.4e1 - t163 / 0.4e1) * m(5);
t410 = t95 * qJD(1);
t398 = -t210 - t293;
t389 = -t281 - t311;
t388 = t486 * t310;
t386 = qJD(3) * t322;
t380 = -t311 + t395;
t378 = t417 / 0.2e1;
t376 = t414 / 0.2e1;
t308 = rSges(3,1) * t322 + rSges(3,2) * t324;
t120 = t185 * t418 + t187 * t270 + t189 * t271;
t121 = t186 * t418 + t188 * t270 + t190 * t271;
t151 = t256 * t418 + t258 * t270 + t260 * t271;
t329 = t322 * t348 + t433;
t85 = -t236 * t270 - t238 * t271 + t319 * t329;
t328 = t322 * t347 + t432;
t86 = -t237 * t270 - t239 * t271 + t319 * t328;
t20 = -(t259 * t270 + t261 * t271) * t324 + t151 * t322 + (t121 * t324 + t322 * t86) * t320 + ((t120 - t424) * t324 + (t85 - t344) * t322) * t319;
t122 = t185 * t415 + t187 * t272 + t189 * t273;
t123 = t186 * t415 + t188 * t272 + t190 * t273;
t152 = t256 * t415 + t258 * t272 + t260 * t273;
t87 = -t236 * t272 - t238 * t273 + t320 * t329;
t88 = -t237 * t272 - t239 * t273 + t320 * t328;
t21 = -(t259 * t272 + t261 * t273) * t324 + t152 * t322 + (t122 * t324 + t322 * t87) * t319 + ((t123 - t424) * t324 + (t88 - t344) * t322) * t320;
t69 = -t151 * t324 + (t120 * t319 + t121 * t320) * t322;
t70 = -t152 * t324 + (t122 * t319 + t123 * t320) * t322;
t371 = t20 * t379 + t21 * t377 + t376 * t70 + t378 * t69 + t484;
t355 = -t205 * t321 + t207 * t323;
t142 = t322 * t355 - t430;
t354 = -t206 * t321 + t208 * t323;
t143 = t322 * t354 - t429;
t358 = t142 * t319 + t143 * t320;
t352 = -t276 * t321 + t278 * t323;
t349 = t471 / 0.2e1 + t371;
t346 = -t274 * t319 - t355;
t345 = -t274 * t320 - t354;
t343 = (Icges(4,3) * t322 + t324 * t361 - t352) * t324;
t341 = t319 * t490 + t299;
t340 = -t320 * t490 + t298;
t50 = t319 * t86 - t320 * t85;
t51 = t319 * t88 - t320 * t87;
t339 = t20 * t464 + t21 * t465 + t51 * t377 + t50 * t379 + (t100 * t319 - t320 * t99) * t461 + (-t120 * t320 + t121 * t319) * t378 + (-t122 * t320 + t123 * t319) * t376 + (-t127 * t320 + t128 * t319) * t462 - t456;
t336 = t37 * t465 + t36 * t464 + (-t105 * t320 + t106 * t319) * t461 + t57 * t377 + t56 * t379 - t70 * t414 / 0.2e1 - t21 * t415 / 0.2e1 - t69 * t417 / 0.2e1 - t20 * t418 / 0.2e1 - t484;
t327 = t322 * t346 + t430;
t326 = t322 * t345 + t429;
t279 = Icges(4,5) * t322 + t324 * t368;
t277 = Icges(4,6) * t322 + t324 * t364;
t250 = t278 * t320;
t249 = t278 * t319;
t248 = t276 * t320;
t247 = t276 * t319;
t226 = t389 * t320;
t224 = t389 * t319;
t215 = t486 * t308;
t180 = -t228 * t324 - t307 * t415;
t179 = t227 * t324 + t307 * t418;
t178 = t380 * t320;
t176 = t380 * t319;
t170 = t322 * t352 - t423;
t165 = (t227 * t320 - t228 * t319) * t322;
t162 = t274 * t415 + t276 * t296 + t278 * t297;
t161 = t274 * t418 + t276 * t294 + t278 * t295;
t157 = -t251 * t319 - t252 * t320 - t388;
t156 = t398 * t324 + (-t289 * t322 + t317 * t459) * t320;
t155 = t292 * t324 - t317 * t384 + t173;
t154 = qJD(4) * t460;
t140 = t193 + (t292 * t320 + t319 * t398) * t322;
t133 = t204 * t415 + t206 * t296 + t208 * t297;
t132 = t203 * t415 + t205 * t296 + t207 * t297;
t131 = t204 * t418 + t206 * t294 + t208 * t295;
t130 = t203 * t418 + t205 * t294 + t207 * t295;
t126 = t397 * t320 + (-t229 - t240) * t319 - t388;
t117 = -t218 * t324 + (-t321 * t399 + t323 * t401) * t322;
t116 = -t217 * t324 + (-t321 * t400 + t323 * t402) * t322;
t112 = -t345 * t324 + (t248 * t321 - t250 * t323 + t204) * t322;
t111 = -t346 * t324 + (t247 * t321 - t249 * t323 + t203) * t322;
t104 = -t248 * t296 - t250 * t297 + t320 * t326;
t103 = -t247 * t296 - t249 * t297 + t320 * t327;
t102 = -t248 * t294 - t250 * t295 + t319 * t326;
t101 = -t247 * t294 - t249 * t295 + t319 * t327;
t96 = (t125 + t163) * t474;
t80 = -t170 * t324 + t322 * t358;
t76 = -t162 * t324 + (t132 * t319 + t133 * t320) * t322;
t75 = -t161 * t324 + (t130 * t319 + t131 * t320) * t322;
t72 = t113 * t163 + (-t175 * t319 - t177 * t320) * t289;
t61 = -t103 * t320 + t104 * t319;
t60 = -t101 * t320 + t102 * t319;
t43 = -(t296 * t391 + t297 * t392) * t324 + (t109 * t319 + (t110 - t421) * t320) * t322;
t42 = -(t294 * t391 + t295 * t392) * t324 + (t108 * t320 + (t107 - t421) * t319) * t322;
t41 = t125 * t153 + t139 * t167 + t141 * t168;
t40 = (t134 + t169) * t475 + (t147 + t82) * t474;
t31 = (t343 + t358) * t324 + (t112 * t320 + t111 * t319 - (-t277 * t321 + t279 * t323 + t274) * t324 + t170) * t322;
t29 = t468 / 0.2e1;
t28 = -(t277 * t296 + t279 * t297) * t324 + t162 * t322 + (t103 * t322 + t132 * t324) * t319 + ((t133 - t423) * t324 + (t104 - t343) * t322) * t320;
t27 = -(t277 * t294 + t279 * t295) * t324 + t161 * t322 + (t102 * t322 + t131 * t324) * t320 + ((t130 - t423) * t324 + (t101 - t343) * t322) * t319;
t22 = t469 / 0.2e1;
t12 = m(5) * t72 + t456;
t9 = m(5) * (t153 * t158 + t167 * t173 + t168 * t174) + t385;
t8 = t9 * qJD(4);
t6 = t456 + t485;
t5 = m(5) * t41 + t371;
t4 = t29 - t469 / 0.2e1 + t349;
t3 = t22 - t468 / 0.2e1 + t349;
t2 = t473 + t470 + (-t31 / 0.2e1 + t76 * t463 + t75 * t465) * t324 + (t80 / 0.2e1 + t28 * t463 + t27 * t465) * t322 + t371;
t1 = t336 + t22 - t471 / 0.2e1 + t29;
t10 = [0, t40 * qJD(3) + t96 * qJD(4) + (-m(3) * t215 / 0.2e1 + t157 * t475 + t126 * t474) * t478, t40 * qJD(2) + (t140 * t474 + t165 * t475) * t477 + t154, qJD(2) * t96 + qJD(3) * t460 + t154; -qJD(3) * t39 - qJD(4) * t95, t6 * qJD(3) + t12 * qJD(4) + (m(4) * (t148 * t157 + t223 * t224 + t225 * t226) + m(5) * (t113 * t126 + t175 * t176 + t177 * t178) + m(3) * (-t215 + t308) * t486 * (rSges(3,1) * t324 - rSges(3,2) * t322) + (t61 + t480 * t299 + (t341 * t320 + (-t298 + t340) * t319) * t320 + t51) * t465 + (t60 + t479 * t298 + (t340 * t319 + (-t299 + t341) * t320) * t319 + t50) * t464) * qJD(2), -t411 + t6 * qJD(2) + t1 * qJD(4) + (-t470 / 0.4e1 - t473 / 0.4e1) * t476 + ((t148 * t165 + t160 * t169 + t179 * t225 + t180 * t223 + (-t171 * t320 - t172 * t319) * t307) * t475 + (t113 * t140 + t119 * t147 + t135 * t233 + t136 * t232 + t155 * t177 + t156 * t175) * t474) * t477 + (-t80 / 0.2e1 + (-t28 / 0.2e1 + t64 / 0.2e1) * t320 + (-t27 / 0.2e1 + t63 / 0.2e1) * t319) * t386 + (t42 * t464 + t43 * t465 + t336 + (t31 / 0.2e1 + (-t76 / 0.2e1 + t116 / 0.2e1) * t320 + (-t75 / 0.2e1 - t117 / 0.2e1) * t319) * t324) * qJD(3), -t410 + t12 * qJD(2) + t1 * qJD(3) + ((t153 * t163 + (-t167 * t320 - t168 * t319) * t289 - t41 + t383) * m(5) + t336) * qJD(4); t39 * qJD(2), t411 + t2 * qJD(3) + t3 * qJD(4) + ((t113 * t82 + t119 * t126 + t135 * t178 + t136 * t176 + t175 * t98 + t177 * t97) * t474 + (t134 * t148 + t149 * t225 + t150 * t223 + t157 * t160 + t171 * t226 + t172 * t224) * t475) * t478 + (t339 + (-t132 * t320 + t133 * t319) * t376 + (-t130 * t320 + t131 * t319) * t378 + (-t111 * t320 + t112 * t319) * t461 + (-t142 * t320 + t143 * t319) * t462 + t27 * t464 + t28 * t465 + t61 * t377 + t60 * t379 - t485) * qJD(2), t2 * qJD(2) + (-t481 * t421 / 0.2e1 + t385) * qJD(3) + t487 + (m(5) * (t119 * t140 + t135 * t155 + t136 * t156) / 0.4e1 + m(4) * (t160 * t165 + t171 * t179 + t172 * t180) / 0.4e1) * t476 + (t43 * t463 + t42 * t465 + (t117 * t320 + t116 * t319 - (-t321 * t391 + t323 * t392) * t324) * t461) * t386, qJD(2) * t3 + qJD(3) * t7 + t487; t95 * qJD(2), t410 + ((t113 * t125 + t126 * t153 + t139 * t177 + t141 * t175 + t167 * t178 + t168 * t176 - t72) * m(5) + t339) * qJD(2) + t4 * qJD(3) + t5 * qJD(4), t4 * qJD(2) + ((t140 * t153 + t155 * t167 + t156 * t168) * m(5) + t385) * qJD(3) + t8, qJD(2) * t5 + qJD(3) * t9 + t8;];
Cq = t10;
