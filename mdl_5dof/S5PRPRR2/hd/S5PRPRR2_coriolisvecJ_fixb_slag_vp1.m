% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:46
% EndTime: 2019-12-05 15:45:05
% DurationCPUTime: 12.43s
% Computational Cost: add. (16564->550), mult. (16791->857), div. (0->0), fcn. (15759->10), ass. (0->315)
t516 = -Icges(5,1) + Icges(5,2);
t281 = qJ(2) + pkin(9);
t273 = sin(t281);
t274 = cos(t281);
t286 = sin(qJ(2));
t288 = cos(qJ(2));
t510 = -Icges(3,5) * t286 - Icges(4,5) * t273 - Icges(3,6) * t288 - Icges(4,6) * t274;
t275 = qJ(4) + t281;
t268 = cos(t275);
t283 = cos(pkin(8));
t287 = cos(qJ(5));
t439 = t283 * t287;
t282 = sin(pkin(8));
t285 = sin(qJ(5));
t442 = t282 * t285;
t230 = t268 * t439 + t442;
t267 = sin(t275);
t280 = qJD(2) + qJD(4);
t448 = t267 * t280;
t413 = t285 * t448;
t126 = -qJD(5) * t230 + t283 * t413;
t440 = t283 * t285;
t441 = t282 * t287;
t229 = -t268 * t440 + t441;
t412 = t287 * t448;
t127 = qJD(5) * t229 - t283 * t412;
t445 = t268 * t280;
t443 = t280 * t283;
t410 = t268 * t443;
t68 = Icges(6,5) * t127 + Icges(6,6) * t126 + Icges(6,3) * t410;
t446 = t267 * t283;
t88 = Icges(6,5) * t230 + Icges(6,6) * t229 + Icges(6,3) * t446;
t344 = t267 * t68 + t445 * t88;
t70 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t410;
t72 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t410;
t455 = Icges(6,4) * t230;
t90 = Icges(6,2) * t229 + Icges(6,6) * t446 + t455;
t214 = Icges(6,4) * t229;
t92 = Icges(6,1) * t230 + Icges(6,5) * t446 + t214;
t20 = t126 * t90 + t127 * t92 + t229 * t70 + t230 * t72 + t283 * t344;
t272 = qJD(2) * t282;
t262 = qJD(4) * t282 + t272;
t421 = qJD(5) * t267;
t204 = t283 * t421 + t262;
t478 = t204 / 0.2e1;
t515 = t20 * t478;
t388 = rSges(6,1) * t287 - rSges(6,2) * t285;
t158 = -rSges(6,3) * t268 + t267 * t388;
t418 = qJD(5) * t282;
t205 = t267 * t418 - t443;
t420 = qJD(5) * t268;
t227 = -t268 * t442 - t439;
t228 = t268 * t441 - t440;
t447 = t267 * t282;
t93 = rSges(6,1) * t228 + rSges(6,2) * t227 + rSges(6,3) * t447;
t514 = t158 * t205 + t93 * t420;
t476 = t205 / 0.2e1;
t511 = t510 * qJD(2);
t509 = 0.2e1 * Icges(5,4) * t267 + t516 * t268;
t508 = -0.2e1 * Icges(5,4) * t268 + t516 * t267;
t462 = Icges(3,4) * t286;
t364 = -Icges(3,2) * t288 - t462;
t461 = Icges(3,4) * t288;
t372 = -Icges(3,1) * t286 - t461;
t459 = Icges(4,4) * t274;
t460 = Icges(4,4) * t273;
t507 = (t286 * t364 - t288 * t372 + t273 * (-Icges(4,2) * t274 - t460) - t274 * (-Icges(4,1) * t273 - t459)) * qJD(2);
t124 = -qJD(5) * t228 + t282 * t413;
t125 = qJD(5) * t227 - t282 * t412;
t444 = t280 * t282;
t411 = t268 * t444;
t67 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t411;
t87 = Icges(6,5) * t228 + Icges(6,6) * t227 + Icges(6,3) * t447;
t345 = t267 * t67 + t445 * t87;
t69 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t411;
t71 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t411;
t456 = Icges(6,4) * t228;
t89 = Icges(6,2) * t227 + Icges(6,6) * t447 + t456;
t213 = Icges(6,4) * t227;
t91 = Icges(6,1) * t228 + Icges(6,5) * t447 + t213;
t19 = t126 * t89 + t127 * t91 + t229 * t69 + t230 * t71 + t283 * t345;
t453 = Icges(6,4) * t287;
t359 = -Icges(6,2) * t285 + t453;
t148 = -Icges(6,6) * t268 + t267 * t359;
t454 = Icges(6,4) * t285;
t367 = Icges(6,1) * t287 - t454;
t150 = -Icges(6,5) * t268 + t267 * t367;
t354 = Icges(6,5) * t287 - Icges(6,6) * t285;
t146 = -Icges(6,3) * t268 + t267 * t354;
t327 = t354 * t268;
t353 = -Icges(6,5) * t285 - Icges(6,6) * t287;
t77 = t280 * t327 + (Icges(6,3) * t280 + qJD(5) * t353) * t267;
t340 = t146 * t445 + t267 * t77;
t43 = t229 * t90 + t230 * t92 + t446 * t88;
t42 = t229 * t89 + t230 * t91 + t446 * t87;
t467 = t282 * t42;
t383 = t283 * t43 + t467;
t328 = t359 * t268;
t358 = -Icges(6,2) * t287 - t454;
t78 = t280 * t328 + (Icges(6,6) * t280 + qJD(5) * t358) * t267;
t329 = t367 * t268;
t366 = -Icges(6,1) * t285 - t453;
t79 = t280 * t329 + (Icges(6,5) * t280 + qJD(5) * t366) * t267;
t310 = (-t126 * t148 - t127 * t150 - t229 * t78 - t230 * t79 + t280 * t383 - t283 * t340) * t268;
t58 = t146 * t446 + t148 * t229 + t150 * t230;
t506 = t19 * t476 + t515 + (t448 * t58 + t310) * qJD(5) / 0.2e1;
t355 = -Icges(5,5) * t267 - Icges(5,6) * t268;
t197 = t355 * t282;
t198 = t355 * t283;
t505 = t197 * t443 / 0.2e1 - t198 * t262 / 0.2e1;
t278 = t282 ^ 2;
t279 = t283 ^ 2;
t426 = t278 + t279;
t152 = t280 * t197;
t153 = t280 * t198;
t500 = -Icges(5,5) * t282 + t509 * t283;
t502 = -Icges(5,6) * t282 + t508 * t283;
t302 = (t500 * t267 + t502 * t268) * t280;
t501 = Icges(5,5) * t283 + t509 * t282;
t503 = Icges(5,6) * t283 + t508 * t282;
t303 = (t501 * t267 + t503 * t268) * t280;
t504 = -(t302 + t152) * t282 / 0.2e1 + (t153 / 0.2e1 - t303 / 0.2e1) * t283;
t499 = t510 * t282;
t498 = t511 * t282;
t497 = t511 * t283;
t496 = t510 * t283;
t365 = -Icges(3,2) * t286 + t461;
t216 = Icges(3,6) * t282 + t283 * t365;
t373 = Icges(3,1) * t288 - t462;
t218 = Icges(3,5) * t282 + t283 * t373;
t363 = -Icges(4,2) * t273 + t459;
t371 = Icges(4,1) * t274 - t460;
t489 = -(Icges(4,6) * t282 + t283 * t363) * t274 - (Icges(4,5) * t282 + t283 * t371) * t273;
t495 = t507 * t283 + (t216 * t288 + t218 * t286 - t489) * qJD(2);
t215 = -Icges(3,6) * t283 + t282 * t365;
t217 = -Icges(3,5) * t283 + t282 * t373;
t490 = (-Icges(4,6) * t283 + t282 * t363) * t274 + (-Icges(4,5) * t283 + t282 * t371) * t273;
t494 = t507 * t282 + (t215 * t288 + t217 * t286 + t490) * qJD(2);
t40 = t227 * t89 + t228 * t91 + t447 * t87;
t41 = t227 * t90 + t228 * t92 + t447 * t88;
t57 = t146 * t447 + t148 * t227 + t150 * t228;
t491 = (t204 * t43 + t205 * t42 - t420 * t58) * t283 + (t204 * t41 + t205 * t40 - t420 * t57) * t282;
t264 = rSges(3,1) * t286 + rSges(3,2) * t288;
t332 = qJD(2) * t264;
t377 = -t285 * t90 + t287 * t92;
t378 = -t285 * t89 + t287 * t91;
t481 = -(-t146 * t283 - t377) * t204 - (-t146 * t282 - t378) * t205;
t480 = -t286 * (t364 * t282 + t217) - t288 * (-t372 * t282 + t215);
t296 = t204 * (-Icges(6,2) * t230 + t214 + t92) + t205 * (-Icges(6,2) * t228 + t213 + t91) - t420 * (t358 * t267 + t150);
t289 = qJD(2) ^ 2;
t479 = -t204 / 0.2e1;
t477 = -t205 / 0.2e1;
t473 = pkin(2) * t286;
t472 = pkin(3) * t274;
t276 = t288 * pkin(2);
t468 = pkin(2) * qJD(2);
t466 = t283 * t41;
t247 = pkin(4) * t268 + pkin(7) * t267;
t194 = t247 * t280;
t334 = t388 * t268;
t387 = -rSges(6,1) * t285 - rSges(6,2) * t287;
t80 = t280 * t334 + (rSges(6,3) * t280 + qJD(5) * t387) * t267;
t463 = -t194 - t80;
t449 = t146 * t268;
t244 = rSges(5,1) * t267 + rSges(5,2) * t268;
t335 = t282 * t244;
t160 = t280 * t335;
t203 = t244 * t283;
t161 = t280 * t203;
t438 = -t282 * t160 - t283 * t161;
t245 = rSges(5,1) * t268 - rSges(5,2) * t267;
t173 = -rSges(5,3) * t283 + t245 * t282;
t174 = rSges(5,3) * t282 + t245 * t283;
t437 = t282 * t173 + t283 * t174;
t246 = pkin(4) * t267 - pkin(7) * t268;
t435 = -t158 - t246;
t177 = -qJ(3) * t283 + t276 * t282;
t178 = qJ(3) * t282 + t276 * t283;
t434 = t282 * t177 + t283 * t178;
t261 = -pkin(3) * t273 - t473;
t433 = t426 * qJD(2) * (t261 + t473);
t417 = t286 * t468;
t422 = qJD(3) * t283;
t254 = -t282 * t417 - t422;
t271 = qJD(3) * t282;
t255 = -t283 * t417 + t271;
t424 = qJD(2) * t283;
t430 = t254 * t272 + t255 * t424;
t429 = t282 * t254 + t283 * t255;
t419 = qJD(5) * t280;
t416 = t288 * t468;
t409 = t268 * t418;
t407 = t177 * t272 + t178 * t424 + qJD(1);
t405 = t424 / 0.2e1;
t404 = -t420 / 0.2e1;
t403 = t420 / 0.2e1;
t402 = t419 / 0.2e1;
t257 = rSges(4,1) * t273 + rSges(4,2) * t274;
t401 = -t257 - t473;
t258 = rSges(4,1) * t274 - rSges(4,2) * t273;
t400 = -t258 - t276;
t399 = t426 * t286;
t398 = -t203 * t443 - t262 * t335;
t339 = t282 * t246;
t164 = t280 * t339;
t207 = t246 * t283;
t165 = t280 * t207;
t74 = rSges(6,1) * t125 + rSges(6,2) * t124 + rSges(6,3) * t411;
t75 = rSges(6,1) * t127 + rSges(6,2) * t126 + rSges(6,3) * t410;
t397 = (-t165 + t75) * t283 + (-t164 + t74) * t282;
t206 = t247 * t282;
t208 = t247 * t283;
t94 = rSges(6,1) * t230 + rSges(6,2) * t229 + rSges(6,3) * t446;
t396 = (t208 + t94) * t283 + (t206 + t93) * t282;
t114 = -pkin(6) * t283 + t282 * t472;
t115 = pkin(6) * t282 + t283 * t472;
t395 = t282 * t114 + t283 * t115 + t434;
t256 = t261 * qJD(2);
t389 = t256 + t417;
t175 = t389 * t282;
t176 = t389 * t283;
t394 = t175 * t272 + t176 * t424 + t430;
t393 = t282 * t175 + t283 * t176 + t429;
t392 = t267 * t402;
t391 = t268 * t402;
t243 = t258 * qJD(2);
t390 = -t243 - t416;
t265 = rSges(3,1) * t288 - rSges(3,2) * t286;
t386 = t88 * t204 + t87 * t205;
t338 = (-t276 - t472) * t289;
t319 = t282 * t338;
t37 = -t194 * t262 - t204 * t80 + t319 + (t94 * t448 + (-t158 * t443 - t75) * t268) * qJD(5);
t318 = t283 * t338;
t38 = -t194 * t443 + t205 * t80 + t318 + (-t93 * t448 + (t158 * t444 + t74) * t268) * qJD(5);
t385 = t282 * t38 - t283 * t37;
t384 = t282 * t40 + t466;
t48 = t267 * t378 - t268 * t87;
t49 = t267 * t377 - t268 * t88;
t382 = t48 * t282 + t49 * t283;
t304 = t256 * t282 - t422;
t51 = -t158 * t204 - t246 * t262 - t420 * t94 + t304;
t314 = t256 * t283 + t271;
t52 = -t246 * t443 + t314 + t514;
t381 = -t282 * t51 - t283 * t52;
t380 = -t282 * t94 + t283 * t93;
t95 = -t244 * t262 + t304;
t96 = -t244 * t443 + t314;
t379 = -t282 * t95 - t283 * t96;
t376 = qJD(2) * t401;
t374 = t114 * t272 + t115 * t424 + t407;
t352 = -t148 * t285 + t150 * t287;
t350 = t426 * t265;
t349 = t426 * t332;
t348 = t282 * t391;
t347 = t283 * t391;
t346 = -t244 + t261;
t341 = -qJD(2) * t243 - t276 * t289;
t337 = Icges(6,3) * t267 + t327 - t352;
t185 = -rSges(4,3) * t283 + t258 * t282;
t186 = rSges(4,3) * t282 + t258 * t283;
t63 = (t185 * t282 + t186 * t283) * qJD(2) + t407;
t336 = t63 * t257;
t333 = (-t262 * t95 - t443 * t96) * t245;
t331 = qJD(2) * t257;
t330 = t261 + t435;
t320 = -qJD(2) * t472 - t416;
t187 = t245 * t280;
t317 = -t187 + t320;
t316 = (t372 * t283 - t216) * t288 + (-t364 * t283 - t218) * t286;
t315 = t320 + t463;
t313 = -(Icges(6,5) * t227 - Icges(6,6) * t228) * t205 - (Icges(6,5) * t229 - Icges(6,6) * t230) * t204 + t353 * t267 * t420;
t312 = (t280 * t382 - (t280 * t352 - t77) * t268 - (t146 * t280 - t285 * t78 + t287 * t79 + (-t148 * t287 - t150 * t285) * qJD(5)) * t267) * t268;
t311 = (-t124 * t148 - t125 * t150 - t227 * t78 - t228 * t79 + t280 * t384 - t282 * t340) * t268;
t305 = t267 * t313;
t122 = t282 * t158;
t301 = -t204 * t122 - t207 * t443 - t262 * t339 + t283 * t514 - t409 * t94;
t295 = (Icges(6,1) * t229 - t455 - t90) * t204 + (Icges(6,1) * t227 - t456 - t89) * t205 - (t366 * t267 - t148) * t420;
t11 = (t280 * t378 - t67) * t268 + (t280 * t87 - t285 * t69 + t287 * t71 + (-t285 * t91 - t287 * t89) * qJD(5)) * t267;
t118 = t148 * t282;
t119 = t148 * t283;
t12 = (t280 * t377 - t68) * t268 + (t280 * t88 - t285 * t70 + t287 * t72 + (-t285 * t92 - t287 * t90) * qJD(5)) * t267;
t120 = t150 * t282;
t121 = t150 * t283;
t149 = Icges(6,6) * t267 + t328;
t151 = Icges(6,5) * t267 + t329;
t17 = t124 * t89 + t125 * t91 + t227 * t69 + t228 * t71 + t282 * t345;
t18 = t124 * t90 + t125 * t92 + t227 * t70 + t228 * t72 + t282 * t344;
t64 = t267 * t352 - t449;
t23 = t204 * t49 + t205 * t48 - t420 * t64;
t290 = (-t337 * t420 - t481) * t267;
t293 = (t502 * t262 - t503 * t443) * t268 + (t500 * t262 - t501 * t443) * t267;
t3 = t17 * t205 + t18 * t204 + (t448 * t57 + t311) * qJD(5);
t294 = -t23 * t421 / 0.2e1 + ((-t119 * t229 - t121 * t230) * t204 + (-t118 * t229 - t120 * t230) * t205 + (t58 * t267 + (-t149 * t229 - t151 * t230 + t467) * t268) * qJD(5)) * t479 + ((-t119 * t227 - t121 * t228) * t204 + (-t118 * t227 - t120 * t228) * t205 + (t57 * t267 + (-t149 * t227 - t151 * t228 + t466) * t268) * qJD(5)) * t477 + (((t119 * t285 - t121 * t287 + t88) * t204 + (t118 * t285 - t120 * t287 + t87) * t205 + t64 * qJD(5)) * t267 + ((t337 * t268 + (t149 * t285 - t151 * t287 - t146) * t267 + t382) * qJD(5) + t481) * t268) * t403 + t491 * t404 + (t18 * t476 + t41 * t348 + t43 * t347 + t49 * t392 + (((t40 - t449) * qJD(5) + t386) * t268 + t290) * t477 + t515 + t506 + t12 * t404 + (t293 / 0.2e1 + t504) * t443) * t282 + (-t17 * t476 - t40 * t348 - t42 * t347 - t48 * t392 + (((t43 - t449) * qJD(5) + t386) * t268 + t290) * t479 - t19 * t478 - t3 / 0.2e1 - t11 * t404 + (-t152 * t283 + t282 * t303 + t505) * t443) * t283 + ((t153 * t282 + t283 * t302 + t505) * t282 + (-t293 / 0.2e1 + t504) * t283) * t262;
t292 = t489 * t282 + t490 * t283;
t159 = rSges(6,3) * t267 + t334;
t291 = t51 * (-t159 * t204 - t247 * t262 + t94 * t421) + t52 * (-t122 * t420 + t158 * t409 + t205 * t159 - t247 * t443 - t421 * t93);
t212 = t387 * t267;
t196 = t283 * t331;
t195 = t282 * t331;
t139 = t341 * t283;
t138 = t341 * t282;
t131 = t283 * t376 + t271;
t130 = t282 * t376 - t422;
t112 = rSges(6,1) * t229 - rSges(6,2) * t230;
t111 = rSges(6,1) * t227 - rSges(6,2) * t228;
t99 = t349 * qJD(2);
t97 = qJD(2) * t350 + qJD(1);
t86 = -t187 * t443 + t318;
t85 = -t187 * t262 + t319;
t73 = (-t195 * t282 - t196 * t283) * qJD(2) + t430;
t50 = -t160 * t262 - t161 * t443 + t394;
t39 = t173 * t262 + t174 * t443 + t374;
t34 = t204 * t93 - t205 * t94 + t206 * t262 + t208 * t443 + t374;
t22 = t268 * t380 * t419 - t164 * t262 - t165 * t443 + t204 * t74 - t205 * t75 + t394;
t1 = [-m(3) * t99 + m(4) * t73 + m(5) * t50 + m(6) * t22; t294 + (-t498 * t279 + (t495 * t282 + (-t494 + t497) * t283) * t282) * t424 - (t496 * qJD(2) * t278 + (-t480 * t283 + t292 + (t316 - t499) * t282) * t424) * t272 / 0.2e1 + t499 * qJD(2) * t279 * t405 + (t292 * t405 + (t494 * t283 + (-t480 - t496) * t405) * t283 + (t316 * t405 + t497 * t282 + (-t495 - t498) * t283) * t282) * t272 + (-t34 * (t301 + t433) - (t381 * t472 + (t288 * t381 - t34 * t399) * pkin(2)) * qJD(2) - t291 + t22 * (t395 + t396) + t34 * (t393 + t397) + (t315 * t52 + t330 * t38) * t283 + (t315 * t51 + t330 * t37) * t282) * m(6) + (-t39 * (t398 + t433) - t333 - (t379 * t472 + (t288 * t379 - t39 * t399) * pkin(2)) * qJD(2) + t50 * (t395 + t437) + t39 * (t393 + t438) + (t317 * t96 + t346 * t86) * t283 + (t317 * t95 + t346 * t85) * t282) * m(5) + (-(-t63 * pkin(2) * t399 + (t131 * t400 - t283 * t336) * t283 + (t130 * t400 - t282 * t336) * t282) * qJD(2) + t73 * t434 + t63 * t429 + (t131 * t390 + t139 * t401 + t73 * t186 - t63 * t196) * t283 + (t130 * t390 + t138 * t401 + t73 * t185 - t63 * t195) * t282) * m(4) + (-t349 * t97 - t350 * t99 + (t264 * t265 * t289 + t332 * t97) * t426) * m(3); m(4) * (-t138 * t283 + t139 * t282) + m(5) * (t282 * t86 - t283 * t85) + m(6) * t385; t294 + (-t291 + t22 * t396 + (t38 * t435 + t463 * t52) * t283 + (t37 * t435 + t463 * t51) * t282 + (-t301 + t397) * t34) * m(6) + (t50 * t437 + (-t282 * t85 - t283 * t86) * t244 + t379 * t187 - t333 + (t438 - t398) * t39) * m(5); t446 * t506 + (t267 * t383 - t268 * t58) * t347 + (t310 + (t19 * t282 + t20 * t283 + t280 * t58) * t267) * t478 + t3 * t447 / 0.2e1 + (t267 * t384 - t268 * t57) * t348 + (t311 + (t17 * t282 + t18 * t283 + t280 * t57) * t267) * t476 + t23 * t448 / 0.2e1 - t268 * (t11 * t205 + t12 * t204 + (t448 * t64 + t312) * qJD(5)) / 0.2e1 + (t267 * t382 - t268 * t64) * t392 + (t312 + (t11 * t282 + t12 * t283 + t280 * t64) * t267) * t404 + (t229 * t296 + t230 * t295 - t283 * t305) * t479 + (t227 * t296 + t228 * t295 - t282 * t305) * t477 + (t313 * t268 + (-t285 * t296 + t295 * t287) * t267) * t403 + t491 * t445 / 0.2e1 + ((-t37 * t94 + t38 * t93 - t51 * t75 + t52 * t74 + (t34 * t380 + (t282 * t52 - t283 * t51) * t158) * t280) * t268 + (t52 * (-t280 * t93 + t282 * t80) + t51 * (t280 * t94 - t283 * t80) + t22 * t380 + t34 * (-t282 * t75 + t283 * t74) + t385 * t158) * t267 - t52 * (t111 * t420 + t205 * t212) - t51 * (-t112 * t420 - t204 * t212) - t34 * (t111 * t204 - t112 * t205)) * m(6);];
tauc = t1(:);
