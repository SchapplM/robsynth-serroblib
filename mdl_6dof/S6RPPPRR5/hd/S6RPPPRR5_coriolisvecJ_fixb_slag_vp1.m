% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:14
% EndTime: 2019-03-09 01:37:36
% DurationCPUTime: 19.20s
% Computational Cost: add. (14440->755), mult. (34441->1005), div. (0->0), fcn. (38337->8), ass. (0->359)
t287 = cos(pkin(9));
t292 = cos(qJ(1));
t467 = sin(pkin(9));
t486 = sin(qJ(1));
t373 = t486 * t467;
t236 = t287 * t292 - t373;
t289 = sin(qJ(5));
t408 = qJD(6) * t289;
t237 = t287 * t486 + t292 * t467;
t413 = qJD(5) * t237;
t190 = t236 * t408 + t413;
t414 = qJD(5) * t236;
t192 = t237 * t408 - t414;
t290 = cos(qJ(6));
t288 = sin(qJ(6));
t291 = cos(qJ(5));
t442 = t288 * t291;
t172 = -t236 * t290 - t237 * t442;
t441 = t290 * t291;
t173 = -t236 * t288 + t237 * t441;
t444 = t237 * t289;
t91 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t444;
t464 = Icges(7,4) * t173;
t94 = Icges(7,2) * t172 + Icges(7,6) * t444 + t464;
t165 = Icges(7,4) * t172;
t97 = Icges(7,1) * t173 + Icges(7,5) * t444 + t165;
t23 = t172 * t94 + t173 * t97 + t91 * t444;
t174 = t236 * t442 - t237 * t290;
t448 = t236 * t289;
t171 = t236 * t441 + t237 * t288;
t463 = Icges(7,4) * t171;
t93 = -Icges(7,2) * t174 + Icges(7,6) * t448 + t463;
t166 = Icges(7,4) * t174;
t96 = Icges(7,1) * t171 + Icges(7,5) * t448 - t166;
t540 = t172 * t93 + t173 * t96;
t90 = Icges(7,5) * t171 - Icges(7,6) * t174 + Icges(7,3) * t448;
t24 = -t444 * t90 - t540;
t407 = qJD(6) * t291;
t262 = qJD(1) - t407;
t349 = Icges(7,5) * t290 - Icges(7,6) * t288;
t204 = -Icges(7,3) * t291 + t289 * t349;
t461 = Icges(7,4) * t290;
t350 = -Icges(7,2) * t288 + t461;
t206 = -Icges(7,6) * t291 + t289 * t350;
t462 = Icges(7,4) * t288;
t352 = Icges(7,1) * t290 - t462;
t208 = -Icges(7,5) * t291 + t289 * t352;
t57 = t172 * t206 + t173 * t208 + t204 * t444;
t10 = -t190 * t24 + t192 * t23 + t57 * t262;
t358 = t288 * t93 - t290 * t96;
t32 = t289 * t358 + t291 * t90;
t26 = t171 * t96 - t174 * t93 + t448 * t90;
t365 = rSges(7,1) * t171 - t174 * rSges(7,2);
t101 = -rSges(7,3) * t448 - t365;
t364 = rSges(7,1) * t290 - rSges(7,2) * t288;
t210 = -rSges(7,3) * t291 + t289 * t364;
t415 = qJD(1) * t292;
t273 = pkin(3) * t415;
t274 = qJD(3) * t292;
t275 = qJD(2) * t486;
t419 = t274 + t275;
t397 = t273 + t419;
t485 = pkin(5) * t289;
t259 = -pkin(8) * t291 + t485;
t411 = qJD(5) * t259;
t537 = -t262 * t101 - t190 * t210 - t237 * t411 + t397;
t518 = t171 * t208 - t174 * t206 + t204 * t448;
t536 = t190 * t26 + t518 * t262;
t246 = Icges(6,5) * t291 - Icges(6,6) * t289;
t143 = -Icges(6,3) * t236 + t237 * t246;
t534 = t237 * t143;
t225 = t237 * qJD(1);
t226 = -qJD(1) * t373 + t287 * t415;
t533 = t226 * pkin(4) + t225 * pkin(7);
t278 = t292 * qJ(3);
t284 = t486 * pkin(3);
t333 = -t284 - t278;
t277 = t486 * qJ(2);
t509 = t292 * pkin(1) + t277;
t320 = t509 - t333;
t276 = qJD(2) * t292;
t388 = qJD(3) * t486;
t368 = t388 - t276;
t532 = -t320 * qJD(1) - t368;
t142 = Icges(6,3) * t237 + t236 * t246;
t280 = Icges(6,4) * t291;
t351 = -Icges(6,2) * t289 + t280;
t145 = Icges(6,6) * t237 + t236 * t351;
t465 = Icges(6,4) * t289;
t250 = Icges(6,1) * t291 - t465;
t148 = Icges(6,5) * t237 + t236 * t250;
t346 = -t145 * t289 + t148 * t291;
t55 = t142 * t237 + t236 * t346;
t75 = -t145 * t291 - t148 * t289;
t251 = rSges(6,1) * t289 + rSges(6,2) * t291;
t412 = qJD(5) * t251;
t529 = t236 * t412;
t394 = t486 * qJ(3);
t405 = t486 * pkin(1);
t319 = -t405 - t394;
t528 = t319 + t394;
t410 = qJD(5) * t289;
t325 = t226 * t291 - t237 * t410;
t409 = qJD(5) * t291;
t451 = t226 * t289;
t326 = t237 * t409 + t451;
t89 = rSges(6,1) * t325 - rSges(6,2) * t326 + t225 * rSges(6,3);
t395 = t278 + t509;
t282 = t486 * rSges(4,1);
t417 = t292 * rSges(4,3) + t282;
t527 = t395 + t417;
t510 = -t292 * rSges(3,2) + t486 * rSges(3,3);
t512 = t509 + t510;
t511 = Icges(6,1) * t289 + t280;
t423 = t511 + t351;
t247 = Icges(6,2) * t291 + t465;
t424 = t247 - t250;
t525 = (t289 * t423 + t291 * t424) * qJD(1);
t220 = t226 * pkin(7);
t521 = t220 + t532;
t268 = qJ(2) * t415;
t377 = t268 + t397;
t233 = t237 * pkin(7);
t513 = t236 * pkin(4) + t233;
t520 = -qJD(1) * t513 + t377 + t533;
t519 = t171 * t97 - t174 * t94;
t517 = 0.2e1 * qJD(5);
t260 = pkin(5) * t291 + pkin(8) * t289;
t515 = t260 * t236;
t205 = Icges(7,3) * t289 + t291 * t349;
t344 = -t206 * t288 + t208 * t290;
t359 = -t288 * t94 + t290 * t97;
t294 = t190 * (t204 * t236 - t358) - t192 * (-t204 * t237 - t359) - t262 * (t205 - t344);
t514 = t294 * t289;
t508 = t236 * rSges(5,1) - t237 * rSges(5,2);
t341 = t247 * t289 - t291 * t511;
t245 = Icges(6,5) * t289 + Icges(6,6) * t291;
t445 = t237 * t245;
t115 = t236 * t341 - t445;
t507 = qJD(1) * t115;
t366 = rSges(6,1) * t291 - rSges(6,2) * t289;
t476 = t237 * rSges(6,3);
t151 = t236 * t366 + t476;
t239 = t351 * qJD(5);
t240 = t250 * qJD(5);
t504 = qJD(5) * (t247 * t291 + t289 * t511) + t239 * t289 - t240 * t291;
t228 = (-Icges(7,2) * t290 - t462) * t289;
t502 = t190 * (Icges(7,2) * t171 + t166 - t96) - t192 * (-Icges(7,2) * t173 + t165 + t97) - t262 * (t208 + t228);
t293 = qJD(1) ^ 2;
t501 = -t10 / 0.2e1;
t25 = -t91 * t448 - t519;
t11 = t192 * t25 - t536;
t500 = t11 / 0.2e1;
t390 = t237 * t407;
t136 = t226 * t408 + (t225 + t390) * qJD(5);
t499 = t136 / 0.2e1;
t391 = t236 * t407;
t137 = t225 * t408 + (-t226 - t391) * qJD(5);
t498 = t137 / 0.2e1;
t497 = t190 / 0.2e1;
t496 = -t190 / 0.2e1;
t495 = -t192 / 0.2e1;
t494 = t192 / 0.2e1;
t493 = t225 / 0.2e1;
t492 = -t226 / 0.2e1;
t489 = -t262 / 0.2e1;
t488 = t262 / 0.2e1;
t487 = -t292 / 0.2e1;
t307 = qJD(6) * t236 - t325;
t370 = t225 - t390;
t80 = t288 * t307 + t290 * t370;
t81 = t288 * t370 - t290 * t307;
t41 = Icges(7,5) * t81 + Icges(7,6) * t80 + Icges(7,3) * t326;
t43 = Icges(7,4) * t81 + Icges(7,2) * t80 + Icges(7,6) * t326;
t45 = Icges(7,1) * t81 + Icges(7,4) * t80 + Icges(7,5) * t326;
t7 = (qJD(5) * t359 - t41) * t291 + (qJD(5) * t91 - t288 * t43 + t290 * t45 + (-t288 * t97 - t290 * t94) * qJD(6)) * t289;
t483 = t7 * t192;
t452 = t225 * t289;
t328 = t236 * t409 - t452;
t327 = t225 * t291 + t236 * t410;
t308 = -qJD(6) * t237 + t327;
t369 = -t226 + t391;
t78 = -t288 * t308 + t290 * t369;
t79 = t288 * t369 + t290 * t308;
t40 = Icges(7,5) * t79 + Icges(7,6) * t78 - Icges(7,3) * t328;
t42 = Icges(7,4) * t79 + Icges(7,2) * t78 - Icges(7,6) * t328;
t44 = Icges(7,1) * t79 + Icges(7,4) * t78 - Icges(7,5) * t328;
t8 = (qJD(5) * t358 - t40) * t291 + (-qJD(5) * t90 - t288 * t42 + t290 * t44 + (t288 * t96 + t290 * t93) * qJD(6)) * t289;
t482 = t8 * t190;
t481 = rSges(7,3) * t289;
t479 = t226 * rSges(6,3);
t475 = t292 * rSges(4,1);
t31 = t289 * t359 - t291 * t91;
t473 = t31 * t136;
t472 = t32 * t137;
t118 = -t204 * t291 + t289 * t344;
t387 = qJD(5) * t408;
t227 = (-Icges(7,5) * t288 - Icges(7,6) * t290) * t289;
t157 = qJD(5) * t205 + qJD(6) * t227;
t207 = Icges(7,6) * t289 + t291 * t350;
t158 = qJD(5) * t207 + qJD(6) * t228;
t209 = Icges(7,5) * t289 + t291 * t352;
t229 = (-Icges(7,1) * t288 - t461) * t289;
t159 = qJD(5) * t209 + qJD(6) * t229;
t39 = (qJD(5) * t344 - t157) * t291 + (qJD(5) * t204 - t158 * t288 + t159 * t290 + (-t206 * t290 - t208 * t288) * qJD(6)) * t289;
t471 = t118 * t387 + t39 * t262;
t107 = pkin(5) * t325 + pkin(8) * t326;
t47 = t81 * rSges(7,1) + t80 * rSges(7,2) + rSges(7,3) * t326;
t470 = t107 + t47;
t146 = -Icges(6,6) * t236 + t237 * t351;
t453 = t146 * t289;
t243 = qJD(5) * t260;
t449 = t236 * t243;
t177 = t236 * t245;
t447 = t236 * t291;
t446 = t237 * t243;
t443 = t237 * t291;
t440 = t292 * t293;
t100 = t173 * rSges(7,1) + t172 * rSges(7,2) + rSges(7,3) * t444;
t187 = pkin(5) * t443 + pkin(8) * t444;
t439 = t100 + t187;
t438 = t101 - t515;
t149 = -Icges(6,5) * t236 + t237 * t250;
t437 = -t149 * t447 - t534;
t436 = t236 * t143 - t149 * t443;
t435 = t237 * t511 + t146;
t434 = t236 * t511 + t145;
t433 = -t247 * t237 + t149;
t432 = t247 * t236 - t148;
t211 = t291 * t364 + t481;
t235 = (-rSges(7,1) * t288 - rSges(7,2) * t290) * t289;
t162 = qJD(5) * t211 + qJD(6) * t235;
t431 = t162 + t243;
t224 = qJD(1) * t509 - t276;
t430 = -t225 * pkin(4) + t220 - t224;
t428 = t210 + t259;
t389 = qJD(1) * t486;
t421 = t268 + t275;
t427 = (-pkin(1) * t389 + t275 + t421) * qJD(1);
t426 = t226 * rSges(5,1) - t225 * rSges(5,2);
t425 = t237 * rSges(5,1) + t236 * rSges(5,2);
t420 = rSges(3,2) * t389 + rSges(3,3) * t415;
t279 = t292 * qJ(2);
t252 = t405 - t279;
t244 = qJD(1) * t252;
t418 = t275 - t244;
t416 = qJD(1) * t246;
t406 = t486 / 0.2e1;
t402 = t486 * rSges(4,3);
t396 = t274 + t418;
t384 = -t414 / 0.2e1;
t383 = t414 / 0.2e1;
t381 = t413 / 0.2e1;
t380 = t410 / 0.2e1;
t267 = qJD(1) * t276;
t378 = -qJD(1) * t224 + t267;
t376 = t273 + t396;
t375 = t284 + t395;
t374 = qJD(6) * t380;
t372 = t79 * rSges(7,1) + t78 * rSges(7,2);
t371 = -t394 - t402 + t475;
t367 = -t225 * rSges(5,1) - t226 * rSges(5,2);
t363 = t23 * t237 - t236 * t24;
t362 = -t236 * t26 + t237 * t25;
t361 = -t236 * t32 + t237 * t31;
t339 = -t252 - t394 + t513;
t62 = -t237 * t412 + (t151 + t339) * qJD(1) + t397;
t152 = rSges(6,1) * t443 - rSges(6,2) * t444 - t236 * rSges(6,3);
t234 = t237 * pkin(4);
t313 = -pkin(7) * t236 + t234 + t320;
t63 = t529 + (t152 + t313) * qJD(1) + t368;
t360 = t236 * t63 - t237 * t62;
t354 = t234 + t375;
t348 = t100 * t236 + t101 * t237;
t74 = t146 * t291 + t149 * t289;
t347 = t149 * t291 - t453;
t345 = t151 * t236 + t152 * t237;
t340 = pkin(4) + t366;
t53 = t142 * t236 + t145 * t444 - t148 * t443;
t253 = rSges(3,2) * t486 + t292 * rSges(3,3);
t331 = t289 * t40 - t409 * t90;
t330 = t289 * t41 + t409 * t91;
t52 = -t146 * t444 - t436;
t324 = (-t236 * t52 - t237 * t53) * qJD(5);
t54 = t146 * t448 + t437;
t323 = (-t236 * t54 - t237 * t55) * qJD(5);
t318 = t190 * t90 + t192 * t91 + t204 * t262;
t316 = 0.2e1 * qJD(1) * t274 - t293 * t394 + t427;
t315 = (Icges(7,5) * t172 - Icges(7,6) * t173) * t192 - (Icges(7,5) * t174 + Icges(7,6) * t171) * t190 + t227 * t262;
t314 = pkin(4) + t260 + t481;
t312 = t289 * t433 + t291 * t435;
t311 = -t289 * t432 + t291 * t434;
t310 = -t394 + t508;
t309 = pkin(3) * t440 + t316;
t306 = t292 * pkin(3) + t279 + t319;
t305 = qJD(1) * t533 + t309;
t88 = rSges(6,1) * t327 + rSges(6,2) * t328 - t479;
t304 = -t151 * t225 + t152 * t226 - t236 * t88 + t237 * t89;
t263 = -0.2e1 * qJD(1) * t388;
t303 = t293 * t333 + t263 + t267;
t302 = t233 + t306;
t301 = -t402 + t319;
t296 = (Icges(7,1) * t172 - t464 - t94) * t192 - (Icges(7,1) * t174 + t463 + t93) * t190 + (-t206 + t229) * t262;
t28 = t100 * t190 + t101 * t192 + qJD(4) + (t187 * t237 + t236 * t515) * qJD(5);
t36 = (t515 + t339) * qJD(1) + t537;
t37 = t236 * t411 + t262 * t100 - t192 * t210 + (t187 + t313) * qJD(1) + t368;
t295 = t28 * t348 + (-t236 * t36 - t237 * t37) * t210;
t271 = rSges(4,1) * t415;
t241 = t366 * qJD(5);
t238 = t246 * qJD(5);
t188 = t259 * t236;
t186 = t259 * t237;
t183 = t251 * t236;
t182 = t251 * t237;
t163 = (-t252 + t371) * qJD(1) + t419;
t161 = -t293 * t510 + t378;
t160 = qJD(1) * t420 + t427;
t141 = -qJ(3) * t440 - t293 * t417 + t263 + t378;
t140 = qJD(1) * (-rSges(4,3) * t389 + t271) + t316;
t139 = t210 * t236;
t138 = t210 * t237;
t135 = t208 * t236;
t134 = t208 * t237;
t133 = t206 * t236;
t132 = t206 * t237;
t125 = (-t252 + t310) * qJD(1) + t397;
t117 = rSges(7,1) * t174 + rSges(7,2) * t171;
t116 = rSges(7,1) * t172 - rSges(7,2) * t173;
t114 = -t237 * t341 - t177;
t106 = pkin(5) * t327 - pkin(8) * t328;
t105 = t114 * qJD(1);
t103 = (-t224 + t367) * qJD(1) + t303;
t102 = qJD(1) * t426 + t309;
t83 = Icges(6,5) * t325 - Icges(6,6) * t326 + Icges(6,3) * t225;
t82 = Icges(6,5) * t327 + Icges(6,6) * t328 - Icges(6,3) * t226;
t64 = qJD(5) * t345 + qJD(4);
t51 = t225 * t245 - t341 * t226 - t236 * t238 - t237 * t504;
t50 = -t341 * t225 - t226 * t245 + t236 * t504 - t237 * t238;
t49 = (-t226 * t251 - t237 * t241) * qJD(5) + (-t88 + t430) * qJD(1) + t303;
t48 = qJD(1) * t89 + (-t225 * t251 + t236 * t241) * qJD(5) + t305;
t46 = -rSges(7,3) * t328 + t372;
t30 = -qJD(5) * t346 + t289 * (Icges(6,1) * t327 + Icges(6,4) * t328 - Icges(6,5) * t226) + t291 * (Icges(6,4) * t327 + Icges(6,2) * t328 - Icges(6,6) * t226);
t29 = qJD(5) * t347 + t289 * (Icges(6,1) * t325 - Icges(6,4) * t326 + Icges(6,5) * t225) + t291 * (Icges(6,4) * t325 - Icges(6,2) * t326 + Icges(6,6) * t225);
t27 = t304 * qJD(5);
t22 = t323 + t507;
t21 = t105 + t324;
t20 = t157 * t444 + t158 * t172 + t159 * t173 + t204 * t326 + t206 * t80 + t208 * t81;
t19 = -t157 * t448 + t158 * t174 - t159 * t171 - t204 * t328 + t206 * t78 + t208 * t79;
t14 = t137 * t210 - t190 * t162 - t262 * t46 + (-t101 * t408 - t226 * t259 - t446) * qJD(5) + (-t106 + t430) * qJD(1) + t303;
t13 = qJD(1) * t107 - t136 * t210 - t192 * t162 + t262 * t47 + (t100 * t408 - t225 * t259 + t449) * qJD(5) + t305;
t12 = t118 * t262 - t190 * t32 + t192 * t31;
t9 = -t100 * t137 + t101 * t136 + t190 * t47 + t192 * t46 + (-t106 * t236 + t107 * t237 + t187 * t226 - t225 * t515) * qJD(5);
t6 = t172 * t42 + t173 * t44 + t237 * t331 - t451 * t90 - t80 * t93 - t81 * t96;
t5 = t172 * t43 + t173 * t45 + t237 * t330 + t451 * t91 + t80 * t94 + t81 * t97;
t4 = -t171 * t44 + t174 * t42 - t236 * t331 - t452 * t90 - t78 * t93 - t79 * t96;
t3 = -t171 * t45 + t174 * t43 - t236 * t330 + t452 * t91 + t78 * t94 + t79 * t97;
t2 = t136 * t23 + t137 * t24 - t190 * t6 + t192 * t5 + t20 * t262 + t387 * t57;
t1 = t136 * t25 + t137 * t26 + t19 * t262 - t190 * t4 + t192 * t3 - t387 * t518;
t15 = [t20 * t494 + (t14 * (t302 + t365) + t13 * (t354 + t439) + (-t13 * pkin(7) + t14 * t314) * t236 + (-t372 - t314 * t225 + (-t485 + (rSges(7,3) + pkin(8)) * t291) * t414 + t521) * t36 + (t244 + t36 + t470 + (-t515 + t528) * qJD(1) + t520 - t537) * t37) * m(7) + t19 * t496 - t518 * t498 + t57 * t499 + t190 * t501 + ((t24 + (t236 * t91 + t237 * t90) * t289 + t519 + t540) * t192 + t11 + t536) * t495 + (t239 * t291 + t240 * t289) * qJD(1) + t483 / 0.2e1 - t482 / 0.2e1 + (-t507 + t22) * t383 + (t29 + t51) * t384 + (t103 * (t306 + t508) + t102 * (t375 + t425) + (t367 + t532) * t125 + (t125 - t376 + t377 + t426 + (-t310 + t319) * qJD(1)) * ((t320 + t425) * qJD(1) + t368)) * m(5) + t10 * t497 + (t141 * (t279 + t301 + t475) + t140 * t527 + (-t368 + (-t282 - t277 + (-rSges(4,3) - pkin(1) - qJ(3)) * t292) * qJD(1)) * t163 + (t163 + t268 + t271 - t396 + t419 + (t301 - t371) * qJD(1)) * (t527 * qJD(1) + t368)) * m(4) + t105 * t381 - (t30 + t50 + t21) * t413 / 0.2e1 + (t160 * t512 + t161 * (-t252 + t253) + (-t418 + t420 + t421 + (-t253 - t405) * qJD(1)) * (qJD(1) * t512 - t276)) * m(3) + (t49 * (t302 + t476) + t48 * (t152 + t354) + (-t48 * pkin(7) + t340 * t49) * t236 + (-t225 * t340 + t479 + t521 - t529) * t62 + (t251 * t413 - t376 + t62 + (-t151 + t528) * qJD(1) + t520 + t89) * t63) * m(6) + t471 + t472 / 0.2e1 + t473 / 0.2e1 + ((t114 + t74) * t493 + (t115 + t75) * t492 - t341 * qJD(1) + ((t52 + (t142 + t453) * t237 + t436) * t237 + (-t53 + (t142 - t347) * t236 - t534) * t236) * t383 + ((-t53 + t54 - t437) * t237 + t436 * t236) * t381) * qJD(5); 0.2e1 * (t13 * t487 + t14 * t406) * m(7) + 0.2e1 * (t406 * t49 + t48 * t487) * m(6) + 0.2e1 * (t102 * t487 + t103 * t406) * m(5) + 0.2e1 * (t140 * t487 + t141 * t406) * m(4) + 0.2e1 * (t160 * t487 + t161 * t406) * m(3); m(4) * (t140 * t486 + t141 * t292) + m(5) * (t102 * t486 + t103 * t292) + m(6) * (t49 * t292 + t48 * t486) + m(7) * (t13 * t486 + t14 * t292); m(6) * t27 + m(7) * t9; (t225 * t31 - t226 * t32 - t236 * t7 - t237 * t8) * t488 + (((t132 * t288 - t134 * t290 + t91) * t192 - (-t133 * t288 + t135 * t290 - t90) * t190 + (-t207 * t288 + t209 * t290 + t204) * t262 + t118 * qJD(6)) * t289 + (qJD(6) * t361 + t294) * t291) * t489 + (t225 * t23 - t226 * t24 - t236 * t5 - t237 * t6) * t494 - qJD(1) * ((-t289 * t424 + t291 * t423) * qJD(1) + ((-t236 * t433 - t237 * t432) * t291 + (t236 * t435 - t237 * t434) * t289) * qJD(5)) / 0.2e1 + ((-t132 * t172 - t134 * t173) * t192 - (t133 * t172 + t135 * t173) * t190 + (t172 * t207 + t173 * t209) * t262 + (-t24 * t447 + t289 * t57) * qJD(6) + ((qJD(6) * t23 + t318) * t291 - t514) * t237) * t495 + (t225 * t25 - t226 * t26 - t236 * t3 - t237 * t4) * t496 + ((-t132 * t174 + t134 * t171) * t192 - (t133 * t174 - t135 * t171) * t190 + (-t171 * t209 + t174 * t207) * t262 + (t25 * t443 - t289 * t518) * qJD(6) + ((-qJD(6) * t26 - t318) * t291 + t514) * t236) * t497 + (-t236 * t25 - t237 * t26) * t498 + (-t23 * t236 - t237 * t24) * t499 + t391 * t500 + t390 * t501 + ((t177 * t413 - t416) * t237 + (t525 + (-t312 * t236 + (-t445 + t311) * t237) * qJD(5)) * t236) * t381 + ((-t414 * t445 - t416) * t236 + (-t525 + (-t311 * t237 + (t177 + t312) * t236) * qJD(5)) * t237) * t383 + (-t236 * t31 - t237 * t32) * t374 - t12 * t408 / 0.2e1 + qJD(1) * (t225 * t74 - t226 * t75 - t236 * t29 - t237 * t30) / 0.2e1 + (-t36 * (-qJD(1) * t188 - t139 * t262 - t190 * t211 - t446) - t37 * (-qJD(1) * t186 - t138 * t262 - t192 * t211 + t449) - t28 * (-t138 * t190 + t139 * t192 - t186 * t413 - t188 * t414) - ((t100 * t37 - t101 * t36) * t289 + t295 * t291) * qJD(6) + (t28 * t439 - t36 * t428) * t226 + (t28 * t438 - t37 * t428) * t225 + (-t14 * t428 + t28 * t470 - t36 * t431 + t439 * t9) * t237 + (t13 * t428 + t37 * t431 - t9 * t438 + t28 * (-t106 - t46)) * t236) * m(7) + (-(-t182 * t63 - t183 * t62) * qJD(1) - (t64 * (-t182 * t237 - t183 * t236) + t360 * t366) * qJD(5) + t27 * t345 + t64 * t304 + t360 * t241 + (-t225 * t63 - t226 * t62 + t236 * t48 - t237 * t49) * t251) * m(6) - (qJD(1) * t51 + t2 + (-(t143 * t225 + t347 * t226 - t236 * t83) * t236 - (-t142 * t225 - t346 * t226 - t236 * t82) * t237 + t225 * t52 - t226 * t53) * t517) * t236 / 0.2e1 - (qJD(1) * t50 + t1 + (-(-t143 * t226 + t347 * t225 - t237 * t83) * t236 - (t142 * t226 - t346 * t225 - t237 * t82) * t237 + t225 * t54 - t226 * t55) * t517) * t237 / 0.2e1 + (t324 + t21 + t10) * t493 + (t323 + t22 + t11) * t492; t2 * t444 / 0.2e1 + (t289 * t363 - t291 * t57) * t499 + ((qJD(5) * t363 - t20) * t291 + (qJD(5) * t57 + t225 * t24 + t226 * t23 - t236 * t6 + t237 * t5) * t289) * t494 + t452 * t500 + t291 * t11 * t384 - t1 * t448 / 0.2e1 + (t289 * t362 + t291 * t518) * t498 + ((qJD(5) * t362 - t19) * t291 + (-qJD(5) * t518 + t225 * t26 + t226 * t25 - t236 * t4 + t237 * t3) * t289) * t496 + t12 * t380 - t291 * (t471 + t472 + t473 - t482 + t483) / 0.2e1 + (-t118 * t291 + t289 * t361) * t374 + ((qJD(5) * t361 - t39) * t291 + (qJD(5) * t118 + t225 * t32 + t226 * t31 - t236 * t8 + t237 * t7) * t289) * t488 + (-t172 * t502 + t173 * t296 + t315 * t444) * t495 + (-t171 * t296 - t174 * t502 - t315 * t448) * t497 + (-t315 * t291 + (t288 * t502 + t296 * t290) * t289) * t489 + (t451 / 0.2e1 + t291 * t381) * t10 + ((qJD(5) * t295 - t13 * t100 + t14 * t101 + t36 * t46 - t37 * t47) * t291 + (t36 * (-qJD(5) * t101 - t162 * t236) + t37 * (qJD(5) * t100 - t162 * t237) + t9 * t348 + t28 * (-t100 * t225 + t101 * t226 + t236 * t47 + t237 * t46) + (-t13 * t237 - t14 * t236 + t225 * t36 - t226 * t37) * t210) * t289 - t36 * (-t117 * t262 - t190 * t235) - t37 * (t116 * t262 - t192 * t235) - t28 * (t116 * t190 + t117 * t192)) * m(7);];
tauc  = t15(:);
