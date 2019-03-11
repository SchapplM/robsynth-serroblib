% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:51
% EndTime: 2019-03-09 01:30:11
% DurationCPUTime: 17.38s
% Computational Cost: add. (14792->738), mult. (19479->1013), div. (0->0), fcn. (17727->8), ass. (0->358)
t261 = qJ(1) + pkin(9);
t258 = sin(t261);
t259 = cos(t261);
t266 = cos(qJ(5));
t262 = sin(qJ(6));
t263 = sin(qJ(5));
t418 = t262 * t263;
t265 = cos(qJ(6));
t420 = t259 * t265;
t165 = t258 * t418 - t420;
t417 = t263 * t265;
t422 = t259 * t262;
t166 = t258 * t417 + t422;
t161 = Icges(7,4) * t166;
t423 = t258 * t266;
t94 = -Icges(7,2) * t165 - Icges(7,6) * t423 + t161;
t160 = Icges(7,4) * t165;
t98 = -Icges(7,1) * t166 + Icges(7,5) * t423 + t160;
t335 = t165 * t94 + t166 * t98;
t424 = t258 * t265;
t167 = -t259 * t418 - t424;
t426 = t258 * t262;
t168 = t259 * t417 - t426;
t419 = t259 * t266;
t435 = Icges(7,4) * t168;
t96 = Icges(7,2) * t167 - Icges(7,6) * t419 + t435;
t162 = Icges(7,4) * t167;
t99 = Icges(7,1) * t168 - Icges(7,5) * t419 + t162;
t455 = t167 * t96 + t168 * t99;
t92 = -Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t423;
t93 = Icges(7,5) * t168 + Icges(7,6) * t167 - Icges(7,3) * t419;
t505 = t335 + t455 + (-t258 * t92 - t259 * t93) * t266;
t456 = t167 * t94 - t168 * t98;
t457 = -t165 * t96 + t166 * t99;
t504 = t456 + t266 * (-t258 * t93 + t259 * t92) + t457;
t320 = Icges(6,5) * t263 + Icges(6,6) * t266;
t149 = -Icges(6,3) * t258 + t259 * t320;
t437 = Icges(6,4) * t263;
t322 = Icges(6,2) * t266 + t437;
t150 = Icges(6,6) * t259 + t258 * t322;
t436 = Icges(6,4) * t266;
t324 = Icges(6,1) * t263 + t436;
t152 = Icges(6,5) * t259 + t258 * t324;
t317 = t150 * t266 + t152 * t263;
t502 = t149 + t317;
t351 = -rSges(4,2) * t259 + t258 * rSges(4,3);
t267 = cos(qJ(1));
t260 = t267 * pkin(1);
t492 = t259 * pkin(2) + t258 * qJ(3);
t371 = t260 + t492;
t497 = t351 + t371;
t326 = t262 * t94 + t265 * t98;
t36 = -t263 * t92 - t266 * t326;
t268 = qJD(1) ^ 2;
t499 = t258 * t259;
t385 = qJD(6) * t263;
t249 = qJD(1) + t385;
t319 = Icges(7,5) * t265 - Icges(7,6) * t262;
t170 = Icges(7,3) * t263 + t266 * t319;
t433 = Icges(7,4) * t265;
t321 = -Icges(7,2) * t262 + t433;
t172 = Icges(7,6) * t263 + t266 * t321;
t434 = Icges(7,4) * t262;
t323 = Icges(7,1) * t265 - t434;
t174 = Icges(7,5) * t263 + t266 * t323;
t285 = t165 * t172 - t166 * t174 + t170 * t423;
t498 = t285 * t249;
t387 = qJD(5) * t263;
t368 = t259 * t387;
t390 = qJD(1) * t266;
t370 = t258 * t390;
t290 = t368 + t370;
t100 = rSges(7,1) * t166 - rSges(7,2) * t165 - rSges(7,3) * t423;
t340 = rSges(7,1) * t265 - rSges(7,2) * t262;
t184 = rSges(7,3) * t263 + t266 * t340;
t464 = pkin(5) * t266;
t231 = pkin(8) * t263 + t464;
t384 = qJD(6) * t266;
t388 = qJD(5) * t259;
t301 = t258 * t384 - t388;
t495 = -t100 * t249 - t184 * t301 + t231 * t388;
t493 = 0.2e1 * qJD(5);
t396 = t258 * rSges(5,2) + t259 * rSges(5,3);
t352 = t259 * rSges(3,1) - rSges(3,2) * t258;
t491 = t260 + t352;
t389 = qJD(5) * t258;
t189 = -t259 * t384 - t389;
t25 = t423 * t92 - t335;
t26 = -t423 * t93 + t457;
t10 = t189 * t26 - t25 * t301 - t498;
t27 = t419 * t92 + t456;
t28 = -t419 * t93 + t455;
t54 = t167 * t172 + t168 * t174 - t170 * t419;
t52 = t54 * t249;
t11 = t189 * t28 - t27 * t301 + t52;
t490 = t10 * t258 + t11 * t259;
t148 = Icges(6,3) * t259 + t258 * t320;
t341 = rSges(6,1) * t263 + rSges(6,2) * t266;
t449 = rSges(6,3) * t259;
t154 = t258 * t341 + t449;
t216 = Icges(6,5) * t266 - Icges(6,6) * t263;
t294 = qJD(5) * t216;
t151 = -Icges(6,6) * t258 + t259 * t322;
t226 = Icges(6,4) * t419;
t421 = t259 * t263;
t432 = Icges(6,5) * t258;
t153 = Icges(6,1) * t421 + t226 - t432;
t316 = t151 * t266 + t153 * t263;
t486 = t259 * t294 + (-t148 + t316) * qJD(1);
t485 = qJD(1) * t502 + t258 * t294;
t218 = -Icges(6,2) * t263 + t436;
t220 = Icges(6,1) * t266 - t437;
t313 = t218 * t266 + t220 * t263;
t484 = t313 * qJD(1) - t320 * qJD(5);
t408 = -Icges(6,2) * t421 + t153 + t226;
t410 = -t220 * t259 + t151;
t483 = -t263 * t410 + t266 * t408;
t169 = Icges(7,3) * t266 - t263 * t319;
t314 = t172 * t262 - t174 * t265;
t325 = t262 * t96 - t265 * t99;
t482 = -t189 * (t170 * t259 + t325) + t301 * (t170 * t258 + t326) - t249 * (t169 + t314);
t195 = (-Icges(7,2) * t265 - t434) * t266;
t481 = t189 * (-Icges(7,2) * t168 + t162 + t99) - t301 * (-Icges(7,2) * t166 - t160 - t98) + t249 * (t174 + t195);
t380 = qJD(5) * qJD(6);
t366 = t263 * t380;
t146 = qJD(1) * t189 + t258 * t366;
t480 = t146 / 0.2e1;
t147 = qJD(1) * t301 + t259 * t366;
t479 = t147 / 0.2e1;
t478 = -t189 / 0.2e1;
t477 = t189 / 0.2e1;
t476 = t301 / 0.2e1;
t475 = -t301 / 0.2e1;
t474 = -t249 / 0.2e1;
t473 = t249 / 0.2e1;
t471 = t258 / 0.2e1;
t470 = -t259 / 0.2e1;
t468 = -rSges(6,3) - pkin(7);
t264 = sin(qJ(1));
t467 = pkin(1) * t264;
t466 = pkin(2) * t258;
t465 = pkin(5) * t263;
t463 = pkin(7) * t259;
t369 = t259 * t390;
t291 = t258 * t387 - t369;
t386 = qJD(5) * t266;
t282 = -t249 * t265 - t262 * t386;
t391 = qJD(1) * t263;
t350 = qJD(6) + t391;
t88 = t258 * t282 - t350 * t422;
t281 = -t249 * t262 + t265 * t386;
t89 = t258 * t281 + t350 * t420;
t45 = Icges(7,5) * t89 + Icges(7,6) * t88 + Icges(7,3) * t291;
t47 = Icges(7,4) * t89 + Icges(7,2) * t88 + Icges(7,6) * t291;
t49 = Icges(7,1) * t89 + Icges(7,4) * t88 + Icges(7,5) * t291;
t7 = (qJD(5) * t326 + t45) * t263 + (-qJD(5) * t92 - t262 * t47 + t265 * t49 + (t262 * t98 - t265 * t94) * qJD(6)) * t266;
t462 = t7 * t301;
t86 = t259 * t282 + t350 * t426;
t87 = t259 * t281 - t350 * t424;
t44 = Icges(7,5) * t87 + Icges(7,6) * t86 + Icges(7,3) * t290;
t46 = Icges(7,4) * t87 + Icges(7,2) * t86 + Icges(7,6) * t290;
t48 = Icges(7,1) * t87 + Icges(7,4) * t86 + Icges(7,5) * t290;
t8 = (qJD(5) * t325 + t44) * t263 + (qJD(5) * t93 - t262 * t46 + t265 * t48 + (-t262 * t99 - t265 * t96) * qJD(6)) * t266;
t461 = t8 * t189;
t460 = -qJD(1) / 0.2e1;
t459 = -pkin(2) - qJ(4);
t194 = (-Icges(7,5) * t262 - Icges(7,6) * t265) * t266;
t127 = qJD(5) * t169 + qJD(6) * t194;
t171 = Icges(7,6) * t266 - t263 * t321;
t128 = qJD(5) * t171 + qJD(6) * t195;
t173 = Icges(7,5) * t266 - t263 * t323;
t196 = (-Icges(7,1) * t262 - t433) * t266;
t129 = qJD(5) * t173 + qJD(6) * t196;
t24 = (qJD(5) * t314 + t127) * t263 + (qJD(5) * t170 - t128 * t262 + t129 * t265 + (-t172 * t265 - t174 * t262) * qJD(6)) * t266;
t365 = t266 * t380;
t67 = t170 * t263 - t266 * t314;
t458 = t24 * t249 + t67 * t365;
t452 = rSges(5,2) * t259;
t451 = rSges(4,3) * t259;
t450 = rSges(6,3) * t258;
t448 = rSges(7,3) * t266;
t102 = t168 * rSges(7,1) + t167 * rSges(7,2) - rSges(7,3) * t419;
t367 = t259 * t386;
t374 = pkin(5) * t367 + pkin(8) * t290;
t425 = t258 * t263;
t378 = pkin(5) * t425;
t125 = -qJD(1) * t378 + t374;
t183 = -t263 * t340 + t448;
t197 = (-rSges(7,1) * t262 - rSges(7,2) * t265) * t266;
t130 = qJD(5) * t183 + qJD(6) * t197;
t230 = pkin(8) * t266 - t465;
t208 = qJD(5) * t230;
t338 = -qJ(4) * t258 - t467;
t292 = t338 - t463;
t381 = qJD(1) * qJD(4);
t382 = qJD(1) * qJD(3);
t393 = qJD(1) * t258;
t392 = qJD(1) * t259;
t238 = qJ(3) * t392;
t247 = qJD(3) * t258;
t401 = t238 + t247;
t407 = qJD(1) * (-pkin(2) * t393 + t401) + t258 * t382;
t375 = 0.2e1 * t259 * t381 + t407;
t277 = t268 * t292 + t375;
t50 = t87 * rSges(7,1) + t86 * rSges(7,2) + rSges(7,3) * t290;
t13 = qJD(1) * t125 - t130 * t189 - t147 * t184 + t249 * t50 + (t102 * t384 + t208 * t258 + t231 * t392) * qJD(5) + t277;
t445 = t13 * t259;
t126 = t291 * pkin(8) + (t258 * t386 + t259 * t391) * pkin(5);
t248 = qJD(3) * t259;
t157 = qJD(1) * t492 - t248;
t244 = t268 * t258 * pkin(7);
t232 = -0.2e1 * t258 * t381;
t237 = t259 * t382;
t251 = t259 * qJ(4);
t337 = -t251 - t260;
t284 = t268 * t337 + t232 + t237;
t344 = rSges(7,1) * t89 + rSges(7,2) * t88;
t51 = rSges(7,3) * t291 + t344;
t14 = -t130 * t301 + t146 * t184 - t249 * t51 + t244 + (-t100 * t384 + t208 * t259) * qJD(5) + (-t231 * t389 - t126 - t157) * qJD(1) + t284;
t444 = t14 * t258;
t441 = t36 * t146;
t37 = t263 * t93 - t266 * t325;
t440 = t37 * t147;
t428 = t148 * t258;
t427 = t148 * t259;
t175 = t216 * t258;
t176 = t216 * t259;
t145 = t259 * t149;
t234 = pkin(8) * t423;
t185 = -t234 + t378;
t416 = -t100 - t185;
t187 = pkin(5) * t421 - pkin(8) * t419;
t415 = t102 + t187;
t414 = t130 + t208;
t413 = t150 * t419 + t152 * t421;
t412 = t151 * t419 + t153 * t421;
t411 = t220 * t258 - t150;
t409 = t218 * t258 + t152;
t405 = t184 + t231;
t404 = -t322 + t220;
t403 = -t218 - t324;
t402 = rSges(6,1) * t421 + rSges(6,2) * t419;
t400 = rSges(4,2) * t393 + rSges(4,3) * t392;
t399 = pkin(7) * t393 + t248;
t246 = qJD(4) * t259;
t398 = t246 + t247;
t252 = t259 * qJ(3);
t198 = -t252 + t466;
t397 = -qJD(1) * t198 + t247;
t394 = qJD(1) * t320;
t111 = t258 * t313 + t176;
t383 = t111 * qJD(1);
t379 = -rSges(5,3) + t459;
t377 = t268 * t467;
t376 = t268 * t260;
t59 = t151 * t423 + t153 * t425 + t145;
t373 = t238 + t398;
t372 = t246 + t397;
t224 = rSges(6,1) * t266 - rSges(6,2) * t263;
t191 = t224 * t388;
t363 = -t392 / 0.2e1;
t361 = t389 / 0.2e1;
t360 = -t388 / 0.2e1;
t354 = t252 - t467;
t349 = t251 + t371;
t348 = t459 - t465;
t347 = t237 - t376;
t200 = rSges(3,1) * t258 + rSges(3,2) * t259;
t334 = t25 * t259 - t258 * t26;
t333 = t25 * t258 + t259 * t26;
t332 = t258 * t28 - t259 * t27;
t331 = t258 * t27 + t259 * t28;
t330 = t258 * t37 - t259 * t36;
t329 = t258 * t36 + t259 * t37;
t312 = rSges(6,1) * t367 - rSges(6,2) * t368;
t113 = -qJD(1) * t154 + t312;
t206 = t341 * qJD(5);
t40 = -t206 * t389 + (t113 + t191) * qJD(1) + t277;
t181 = t224 * t258;
t114 = qJD(5) * t181 + (t259 * t341 - t450) * qJD(1);
t41 = t232 + t244 + (-qJ(4) * t268 - qJD(5) * t206) * t259 + (-t224 * t389 - t114 - t157) * qJD(1) + t347;
t328 = t40 * t258 + t41 * t259;
t287 = -t198 + t292;
t62 = t191 + (-t154 + t287) * qJD(1) + t398;
t155 = t402 - t450;
t310 = t492 - t337;
t63 = (qJD(5) * t224 + qJD(4)) * t258 + (t155 + t310) * qJD(1) - t399;
t327 = t258 * t63 + t259 * t62;
t318 = t100 * t259 - t102 * t258;
t80 = -t150 * t263 + t152 * t266;
t81 = -t151 * t263 + t153 * t266;
t315 = -t154 * t258 - t155 * t259;
t311 = t354 - t463;
t309 = -rSges(5,3) * t258 + t338 + t452;
t306 = -t266 * t44 + t387 * t93;
t305 = -t266 * t45 - t387 * t92;
t58 = t258 * t317 + t427;
t299 = (-t258 * t59 + t259 * t58) * qJD(5);
t60 = t413 - t428;
t61 = -t149 * t258 + t412;
t298 = (-t258 * t61 + t259 * t60) * qJD(5);
t297 = -t341 + t459;
t296 = qJD(5) * t220;
t295 = qJD(5) * t218;
t289 = t170 * t249 + t189 * t93 + t301 * t92;
t288 = -(-Icges(7,5) * t165 - Icges(7,6) * t166) * t301 + (Icges(7,5) * t167 - Icges(7,6) * t168) * t189 + t194 * t249;
t286 = t263 * t411 + t266 * t409;
t283 = t266 * t288;
t280 = (t263 * t403 + t266 * t404) * qJD(1);
t278 = (Icges(7,1) * t167 - t435 - t96) * t189 - (-Icges(7,1) * t165 - t161 - t94) * t301 + (-t172 + t196) * t249;
t276 = -t113 * t259 - t114 * t258 + (-t154 * t259 + t155 * t258) * qJD(1);
t29 = t100 * t189 + t102 * t301 + qJD(2) + (-t185 * t258 - t187 * t259) * qJD(5);
t33 = (-t185 + t287) * qJD(1) + t398 + t495;
t34 = t102 * t249 - t184 * t189 + (qJD(5) * t231 + qJD(4)) * t258 + (t187 + t310) * qJD(1) - t399;
t273 = t29 * t318 + (t258 * t33 - t259 * t34) * t184;
t204 = t322 * qJD(5);
t205 = t324 * qJD(5);
t270 = -qJD(1) * t216 - t204 * t266 - t205 * t263 + (-t218 * t263 + t220 * t266) * qJD(5);
t269 = t482 * t266;
t243 = rSges(5,2) * t392;
t188 = t231 * t259;
t186 = t231 * t258;
t182 = t224 * t259;
t144 = t184 * t259;
t143 = t184 * t258;
t142 = t174 * t259;
t141 = t174 * t258;
t140 = t172 * t259;
t139 = t172 * t258;
t124 = rSges(7,1) * t167 - rSges(7,2) * t168;
t123 = -rSges(7,1) * t165 - rSges(7,2) * t166;
t112 = t259 * t313 - t175;
t90 = t112 * qJD(1);
t83 = (-qJD(1) * t351 - t157) * qJD(1) + t347;
t82 = qJD(1) * t400 - t377 + t407;
t79 = qJD(4) * t258 - t248 + (t310 + t396) * qJD(1);
t78 = (-t198 + t309) * qJD(1) + t398;
t74 = qJD(5) * t315 + qJD(2);
t66 = (-qJD(1) * t396 - t157) * qJD(1) + t284;
t65 = qJD(1) * (-rSges(5,3) * t393 + t243) + t338 * t268 + t375;
t43 = t270 * t258 + t259 * t484;
t42 = -t258 * t484 + t270 * t259;
t39 = -qJD(5) * t316 - (-qJD(1) * t150 + t259 * t295) * t263 + (-qJD(1) * t152 + t259 * t296) * t266;
t38 = -qJD(5) * t317 - (qJD(1) * t151 + t258 * t295) * t263 + (t258 * t296 + (t259 * t324 - t432) * qJD(1)) * t266;
t35 = t276 * qJD(5);
t23 = t90 + t298;
t22 = t299 + t383;
t16 = -t127 * t423 - t128 * t165 + t129 * t166 + t170 * t291 + t172 * t88 + t174 * t89;
t15 = -t127 * t419 + t128 * t167 + t129 * t168 + t170 * t290 + t172 * t86 + t174 * t87;
t12 = t189 * t37 + t249 * t67 - t301 * t36;
t9 = t100 * t147 - t102 * t146 + t189 * t51 + t301 * t50 + (-t125 * t259 - t126 * t258 + (-t185 * t259 + t187 * t258) * qJD(1)) * qJD(5);
t6 = -t165 * t46 + t166 * t48 + t258 * t306 - t369 * t93 + t88 * t96 + t89 * t99;
t5 = -t165 * t47 + t166 * t49 + t258 * t305 + t369 * t92 + t88 * t94 - t89 * t98;
t4 = t167 * t46 + t168 * t48 + t259 * t306 + t370 * t93 + t86 * t96 + t87 * t99;
t3 = t167 * t47 + t168 * t49 + t259 * t305 - t370 * t92 + t86 * t94 - t87 * t98;
t2 = t146 * t25 + t147 * t26 + t16 * t249 + t189 * t6 - t285 * t365 - t301 * t5;
t1 = t146 * t27 + t147 * t28 + t15 * t249 + t189 * t4 - t3 * t301 + t365 * t54;
t17 = [(t52 - (-t26 + t504) * t301 + (t25 + t505) * t189) * t476 + (t90 + (t413 * t259 + (t258 * t502 - t412 - t58) * t258) * qJD(5)) * t360 + t461 / 0.2e1 - t462 / 0.2e1 + t54 * t479 + t15 * t477 + t458 + (-qJD(5) * t313 + t204 * t263 - t205 * t266) * qJD(1) + t440 / 0.2e1 + t441 / 0.2e1 - t285 * t480 + m(3) * ((-t200 * t268 - t377) * t491 + (-t376 + (-0.2e1 * t352 - t260 + t491) * t268) * (-t200 - t467)) + (t498 - (-t28 + t505) * t301 + (t27 - t504) * t189 + t10) * t478 + (t16 + t11) * t475 - (t42 + t39) * t389 / 0.2e1 + (t22 - t383 + ((t412 - t61 - t427) * t259 + (-t145 + t59 - t60 - t428) * t258) * qJD(5)) * t361 + (-(-t33 + t372 + (-t185 + t292) * qJD(1) + t495) * t34 + t14 * (-t100 + t234 + t311) + t33 * (-t344 + t399) + t13 * (t349 + t415) + t34 * (t50 + t373 + t374) + (-t13 * pkin(7) + t14 * t348 + (-qJD(4) + (-t464 + (-rSges(7,3) - pkin(8)) * t263) * qJD(5)) * t33) * t258 + ((-t264 * t34 - t267 * t33) * pkin(1) + (-t33 * qJ(3) + t34 * t348) * t258 + (t33 * (t230 + t448 + t459) - t34 * pkin(7)) * t259) * qJD(1)) * m(7) + (-(-t62 + t191 + t372 + (-t154 + t292) * qJD(1)) * t63 + t41 * (t311 - t449) + t62 * t399 + t40 * (t349 + t402) + t63 * (t312 + t373) + (t41 * t297 + t62 * (-rSges(6,1) * t386 + rSges(6,2) * t387 - qJD(4)) + t40 * t468) * t258 + ((-t264 * t63 - t267 * t62) * pkin(1) + (t297 * t62 + t468 * t63) * t259 + (t62 * (rSges(6,3) - qJ(3)) + t63 * t297) * t258) * qJD(1)) * m(6) + (-(qJD(1) * t309 + t372 - t78) * t79 + t66 * (t354 + t452) + t78 * t248 + t65 * (t349 + t396) + t79 * (t243 + t373) + (-t78 * qJD(4) + t379 * t66) * t258 + ((-t264 * t79 - t267 * t78) * pkin(1) + t78 * t379 * t259 + (t78 * (-rSges(5,2) - qJ(3)) + t79 * t379) * t258) * qJD(1)) * m(5) + (t83 * (t451 + (rSges(4,2) - pkin(2)) * t258 + t354) + t82 * t497 + (t400 + t401 - t397 + (-rSges(4,2) * t258 - t451 - t466) * qJD(1)) * (qJD(1) * t497 - t248)) * m(4) + (t38 + t43 + t23) * t388 / 0.2e1 + ((t111 + t80) * t258 + (t112 + t81) * t259) * qJD(5) * t460; m(6) * t35 + m(7) * t9; 0.2e1 * (-t445 / 0.2e1 + t444 / 0.2e1) * m(7) + 0.2e1 * (t40 * t470 + t41 * t471) * m(6) + 0.2e1 * (t470 * t65 + t471 * t66) * m(5) + 0.2e1 * (t470 * t82 + t471 * t83) * m(4); m(5) * (t258 * t65 + t259 * t66) + m(6) * t328 + m(7) * (t13 * t258 + t14 * t259); t334 * t480 + (-qJD(1) * t329 - t258 * t8 + t259 * t7) * t473 + (((-t140 * t262 + t142 * t265 + t93) * t189 - (-t139 * t262 + t141 * t265 - t92) * t301 + (-t171 * t262 + t173 * t265 + t170) * t249 + t67 * qJD(6)) * t266 + (qJD(6) * t329 - t482) * t263) * t474 + (-qJD(1) * t333 - t258 * t6 + t259 * t5) * t475 + ((-t140 * t165 + t142 * t166) * t189 - (-t139 * t165 + t141 * t166) * t301 + (-t165 * t171 + t166 * t173) * t249 + (t26 * t421 - t266 * t285) * qJD(6) + ((qJD(6) * t25 + t289) * t263 + t269) * t258) * t476 + (-qJD(1) * t331 - t258 * t4 + t259 * t3) * t477 + ((t140 * t167 + t142 * t168) * t189 - (t139 * t167 + t141 * t168) * t301 + (t167 * t171 + t168 * t173) * t249 + (t266 * t54 + t27 * t425) * qJD(6) + ((qJD(6) * t28 + t289) * t263 + t269) * t259) * t478 + ((-t263 * t404 + t266 * t403) * qJD(1) + ((t258 * t410 + t259 * t411) * t266 + (t258 * t408 - t259 * t409) * t263) * qJD(5)) * t460 - t147 * t332 / 0.2e1 - t12 * t384 / 0.2e1 + ((t176 * t389 + t394) * t258 + (t280 + (t286 * t259 + (-t175 - t483) * t258) * qJD(5)) * t259) * t361 + ((t175 * t388 - t394) * t259 + (t280 + (-t483 * t258 + (-t176 + t286) * t259) * qJD(5)) * t258) * t360 + qJD(1) * (-t258 * t39 + t259 * t38 + (-t258 * t80 - t81 * t259) * qJD(1)) / 0.2e1 - t330 * t365 / 0.2e1 - t490 * t385 / 0.2e1 + ((t14 * t405 + t33 * t414 - t9 * t415 + t29 * (-t125 - t50) + (t29 * t416 + t34 * t405) * qJD(1)) * t259 + (t13 * t405 + t34 * t414 + t9 * t416 + t29 * (-t126 - t51) + (t29 * t415 - t33 * t405) * qJD(1)) * t258 - t33 * (-qJD(1) * t186 - t143 * t249 - t183 * t301 + t230 * t388) - t34 * (qJD(1) * t188 + t144 * t249 - t183 * t189 + t230 * t389) - t29 * (t143 * t189 + t144 * t301 - t186 * t389 - t188 * t388) - ((-t100 * t33 + t102 * t34) * t266 + t273 * t263) * qJD(6)) * m(7) + (-(-t181 * t62 + t182 * t63) * qJD(1) - (t74 * (-t181 * t258 - t182 * t259) - t327 * t341) * qJD(5) + t35 * t315 + t74 * t276 - t327 * t206 + ((-t258 * t62 + t259 * t63) * qJD(1) + t328) * t224) * m(6) - (qJD(1) * t42 + t1 + (-t485 * t499 + t486 * t258 ^ 2 + (-t60 * t258 - t61 * t259) * qJD(1)) * t493) * t258 / 0.2e1 + (qJD(1) * t43 + t2 + (t485 * t259 ^ 2 - t486 * t499 + (-t58 * t258 - t59 * t259) * qJD(1)) * t493) * t259 / 0.2e1 - (t299 + t22 + t10) * t393 / 0.2e1 + (t23 + t11 + t298) * t363; t11 * t370 / 0.2e1 - t1 * t419 / 0.2e1 + (t263 * t54 - t266 * t331) * t479 + ((qJD(5) * t331 + t15) * t263 + (qJD(1) * t332 + qJD(5) * t54 - t258 * t3 - t259 * t4) * t266) * t477 + t266 * t10 * t363 - t2 * t423 / 0.2e1 + (-t263 * t285 - t266 * t333) * t480 + ((qJD(5) * t333 + t16) * t263 + (-qJD(1) * t334 - qJD(5) * t285 - t258 * t5 - t259 * t6) * t266) * t475 + t263 * (t440 + t441 + t458 + t461 - t462) / 0.2e1 + ((qJD(5) * t329 + t24) * t263 + (qJD(1) * t330 + qJD(5) * t67 - t258 * t7 - t259 * t8) * t266) * t473 + (t167 * t481 + t278 * t168 - t259 * t283) * t478 + (-t165 * t481 + t166 * t278 - t258 * t283) * t476 + (t288 * t263 + (-t262 * t481 + t265 * t278) * t266) * t474 + t490 * t387 / 0.2e1 + (t12 + qJD(6) * (t263 * t67 - t266 * t329)) * t386 / 0.2e1 + ((qJD(5) * t273 - t14 * t100 + t13 * t102 - t33 * t51 + t34 * t50) * t263 + (t33 * (-qJD(5) * t100 - t130 * t258) + t34 * (qJD(5) * t102 + t130 * t259) - t9 * t318 + t29 * (t100 * t393 + t102 * t392 + t258 * t50 - t259 * t51) + (t445 - t444 + (-t258 * t34 - t259 * t33) * qJD(1)) * t184) * t266 - t33 * (-t123 * t249 - t197 * t301) - t34 * (t124 * t249 - t189 * t197) - t29 * (t123 * t189 + t124 * t301)) * m(7);];
tauc  = t17(:);
