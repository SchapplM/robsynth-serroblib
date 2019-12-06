% Calculate vector of inverse dynamics joint torques for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:14
% DurationCPUTime: 9.12s
% Computational Cost: add. (6732->607), mult. (10817->786), div. (0->0), fcn. (8203->6), ass. (0->336)
t254 = qJ(4) + qJ(5);
t243 = sin(t254);
t244 = cos(t254);
t176 = rSges(6,1) * t244 - rSges(6,2) * t243;
t253 = qJD(4) + qJD(5);
t256 = sin(qJ(1));
t397 = t253 * t256;
t471 = t176 * t397;
t258 = cos(qJ(1));
t414 = Icges(6,4) * t243;
t305 = Icges(6,2) * t244 + t414;
t122 = Icges(6,6) * t258 + t256 * t305;
t174 = Icges(6,1) * t244 - t414;
t470 = t174 * t256 - t122;
t123 = -Icges(6,6) * t256 + t258 * t305;
t469 = t174 * t258 - t123;
t413 = Icges(6,4) * t244;
t307 = Icges(6,1) * t243 + t413;
t124 = Icges(6,5) * t258 + t256 * t307;
t172 = -Icges(6,2) * t243 + t413;
t468 = t172 * t256 + t124;
t467 = -t172 - t307;
t462 = -t256 * rSges(4,2) - t258 * rSges(4,3);
t202 = t258 * pkin(1) + t256 * qJ(2);
t246 = t258 * qJ(3);
t463 = t246 + t202;
t466 = -t462 + t463;
t398 = t244 * t258;
t400 = t243 * t258;
t154 = rSges(6,1) * t398 - rSges(6,2) * t400;
t465 = -t305 + t174;
t464 = t122 * t398 + t124 * t400;
t204 = -rSges(3,2) * t258 + t256 * rSges(3,3);
t303 = Icges(6,5) * t243 + Icges(6,6) * t244;
t120 = Icges(6,3) * t258 + t256 * t303;
t461 = qJD(1) * t120;
t319 = rSges(6,1) * t243 + rSges(6,2) * t244;
t126 = rSges(6,3) * t258 + t256 * t319;
t460 = qJD(1) * t126;
t255 = sin(qJ(4));
t257 = cos(qJ(4));
t304 = Icges(5,5) * t255 + Icges(5,6) * t257;
t138 = Icges(5,3) * t258 + t256 * t304;
t459 = qJD(1) * t138;
t399 = t244 * t256;
t401 = t243 * t256;
t153 = rSges(6,1) * t399 - rSges(6,2) * t401;
t189 = t253 * t258;
t424 = rSges(6,3) * t256;
t127 = rSges(6,1) * t400 + rSges(6,2) * t398 - t424;
t259 = -pkin(7) - pkin(6);
t395 = t255 * t258;
t371 = pkin(4) * t395 + t256 * t259;
t434 = pkin(6) * t256;
t157 = t371 + t434;
t396 = t255 * t256;
t222 = pkin(4) * t396;
t158 = -t222 + (pkin(6) + t259) * t258;
t42 = -t126 * t397 - t127 * t189 + (-t157 * t258 + t158 * t256) * qJD(4);
t247 = t258 * qJ(2);
t197 = pkin(1) * t256 - t247;
t404 = qJ(3) * t256;
t433 = pkin(6) * t258;
t318 = -t404 - t433;
t381 = -t126 + t158;
t283 = t318 + t381;
t278 = -t197 + t283;
t356 = qJD(4) * t258;
t343 = t257 * t356;
t214 = pkin(4) * t343;
t240 = qJD(3) * t258;
t241 = qJD(2) * t256;
t366 = t240 + t241;
t290 = t176 * t189 + t214 + t366;
t43 = qJD(1) * t278 + t290;
t357 = qJD(4) * t257;
t322 = pkin(4) * t357 + qJD(3);
t361 = qJD(1) * t256;
t235 = pkin(6) * t361;
t242 = qJD(2) * t258;
t367 = t235 + t242;
t380 = t127 + t157;
t44 = t471 + t322 * t256 + (t463 + t380) * qJD(1) - t367;
t458 = -t44 * (qJD(1) * t154 - t319 * t397) - t42 * (-t153 * t397 - t154 * t189) - t43 * (-qJD(1) * t153 - t189 * t319);
t191 = Icges(5,5) * t257 - Icges(5,6) * t255;
t284 = qJD(4) * t191;
t416 = Icges(5,4) * t255;
t306 = Icges(5,2) * t257 + t416;
t141 = -Icges(5,6) * t256 + t258 * t306;
t390 = t257 * t258;
t217 = Icges(5,4) * t390;
t412 = Icges(5,5) * t256;
t143 = Icges(5,1) * t395 + t217 - t412;
t299 = t141 * t257 + t143 * t255;
t457 = qJD(1) * t299 + t258 * t284 - t459;
t140 = Icges(5,6) * t258 + t256 * t306;
t415 = Icges(5,4) * t257;
t308 = Icges(5,1) * t255 + t415;
t142 = Icges(5,5) * t258 + t256 * t308;
t300 = t140 * t257 + t142 * t255;
t139 = -Icges(5,3) * t256 + t258 * t304;
t363 = qJD(1) * t139;
t456 = qJD(1) * t300 + t256 * t284 + t363;
t208 = Icges(6,4) * t398;
t410 = Icges(6,5) * t256;
t125 = Icges(6,1) * t400 + t208 - t410;
t170 = Icges(6,5) * t244 - Icges(6,6) * t243;
t148 = t170 * t258;
t455 = qJD(1) * (t123 * t244 + t125 * t243) + t253 * t148 - t461;
t147 = t170 * t256;
t302 = t122 * t244 + t124 * t243;
t121 = -Icges(6,3) * t256 + t258 * t303;
t364 = qJD(1) * t121;
t454 = qJD(1) * t302 + t147 * t253 + t364;
t320 = rSges(5,1) * t255 + rSges(5,2) * t257;
t144 = rSges(5,3) * t258 + t256 * t320;
t297 = t172 * t244 + t174 * t243;
t453 = qJD(1) * t297 - t303 * t253;
t193 = -Icges(5,2) * t255 + t415;
t195 = Icges(5,1) * t257 - t416;
t296 = t257 * t193 + t255 * t195;
t452 = t296 * qJD(1) - t304 * qJD(4);
t451 = t256 * (-Icges(5,2) * t395 + t143 + t217) - t258 * (t193 * t256 + t142);
t450 = qJD(1) * t465 - t397 * (-Icges(6,2) * t400 + t125 + t208) + t189 * t468;
t325 = qJD(1) * t253;
t131 = (-qJDD(4) - qJDD(5)) * t256 - t258 * t325;
t449 = t131 / 0.2e1;
t236 = qJDD(4) * t258;
t132 = qJDD(5) * t258 - t256 * t325 + t236;
t448 = t132 / 0.2e1;
t352 = qJD(1) * qJD(4);
t184 = -qJDD(4) * t256 - t258 * t352;
t447 = t184 / 0.2e1;
t185 = -t256 * t352 + t236;
t446 = t185 / 0.2e1;
t445 = t397 / 0.2e1;
t444 = -t397 / 0.2e1;
t443 = -t189 / 0.2e1;
t442 = t189 / 0.2e1;
t441 = -t256 / 0.2e1;
t440 = t256 / 0.2e1;
t439 = -t258 / 0.2e1;
t438 = t258 / 0.2e1;
t437 = rSges(3,2) - pkin(1);
t436 = -rSges(5,3) - pkin(6);
t435 = pkin(4) * t255;
t432 = -qJD(1) / 0.2e1;
t431 = qJD(1) / 0.2e1;
t430 = -pkin(1) - qJ(3);
t428 = rSges(5,2) * t255;
t427 = rSges(3,3) * t258;
t426 = rSges(5,3) * t256;
t422 = t256 * t43;
t201 = rSges(5,1) * t257 - t428;
t168 = t201 * t356;
t289 = -t144 + t318;
t282 = -t197 + t289;
t58 = qJD(1) * t282 + t168 + t366;
t421 = t256 * t58;
t420 = qJDD(1) / 0.2e1;
t419 = -rSges(6,3) + t259;
t359 = qJD(1) * t259;
t360 = qJD(1) * t258;
t106 = t256 * t359 + t235 + (t255 * t360 + t256 * t357) * pkin(4);
t74 = t471 + (t258 * t319 - t424) * qJD(1);
t418 = -t106 - t74;
t373 = t258 * t359 + t214;
t107 = (-t222 + t433) * qJD(1) + t373;
t293 = t154 * t253;
t73 = t293 - t460;
t417 = t107 + t73;
t403 = t138 * t256;
t402 = t138 * t258;
t160 = t191 * t256;
t161 = t191 * t258;
t394 = t255 * qJD(4) ^ 2;
t393 = t256 * t120;
t392 = t256 * t257;
t261 = qJD(1) ^ 2;
t391 = t256 * t261;
t389 = t258 * t120;
t128 = t258 * t139;
t63 = t256 * t297 + t148;
t388 = t63 * qJD(1);
t81 = t256 * t296 + t161;
t387 = t81 * qJD(1);
t386 = t123 * t398 + t125 * t400;
t385 = t140 * t390 + t142 * t395;
t384 = t141 * t390 + t143 * t395;
t376 = -t306 + t195;
t375 = -t193 - t308;
t199 = rSges(3,2) * t256 + t427;
t374 = -t197 + t199;
t156 = t202 + t204;
t372 = rSges(5,1) * t395 + rSges(5,2) * t390;
t354 = qJD(1) * qJD(2);
t370 = qJDD(2) * t256 + t258 * t354;
t229 = qJ(2) * t360;
t369 = t229 + t241;
t368 = rSges(3,2) * t361 + rSges(3,3) * t360;
t187 = qJD(1) * t197;
t365 = t241 - t187;
t362 = qJD(1) * t304;
t358 = qJD(4) * t256;
t355 = -m(4) - m(5) - m(6);
t353 = qJD(1) * qJD(3);
t351 = qJDD(2) * t258;
t350 = -rSges(4,3) + t430;
t225 = pkin(4) * t390;
t348 = qJ(3) * t391;
t48 = t258 * t121 + t123 * t399 + t125 * t401;
t52 = t141 * t392 + t143 * t396 + t128;
t347 = qJD(1) * (-pkin(1) * t361 + t369) + qJDD(1) * t202 + t256 * t354;
t346 = t229 + t366;
t345 = t240 + t365;
t342 = -t361 / 0.2e1;
t341 = -t360 / 0.2e1;
t340 = -t358 / 0.2e1;
t339 = t358 / 0.2e1;
t338 = -t356 / 0.2e1;
t337 = t356 / 0.2e1;
t203 = rSges(4,2) * t258 - t256 * rSges(4,3);
t336 = (t258 * t307 - t410) * qJD(1) + t470 * t253;
t335 = -qJD(1) * t124 + t469 * t253;
t334 = qJD(1) * t123 + t468 * t253;
t333 = -qJD(1) * t122 + t125 * t253 + t172 * t189;
t331 = t203 - t404;
t330 = t465 * t253;
t329 = t467 * t253;
t327 = -qJD(3) * t256 + t242;
t326 = qJD(4) * t201 + qJD(3);
t324 = qJDD(3) * t258 - 0.2e1 * t256 * t353 + t370;
t137 = t319 * t253;
t323 = -qJD(4) * t435 - t137;
t321 = -t197 + t331;
t205 = rSges(2,1) * t258 - rSges(2,2) * t256;
t200 = rSges(2,1) * t256 + rSges(2,2) * t258;
t159 = qJD(1) * t202 - t242;
t234 = pkin(6) * t391;
t280 = -t246 * t261 + t324;
t14 = t132 * t176 - t137 * t189 + t234 + (t185 * t257 - t258 * t394) * pkin(4) + (-t159 + t418) * qJD(1) + t278 * qJDD(1) + t280;
t288 = qJDD(1) * t246 + qJDD(3) * t256 + 0.2e1 * t258 * t353 + t347;
t272 = (-pkin(6) * t261 - qJDD(2)) * t258 + t288;
t15 = -t348 - t131 * t176 - t137 * t397 + t417 * qJD(1) + (-t184 * t257 - t256 * t394) * pkin(4) + (t380 - t434) * qJDD(1) + t272;
t317 = t14 * t258 + t15 * t256;
t80 = -t141 * t255 + t143 * t257;
t285 = qJD(4) * t193;
t88 = -qJD(1) * t140 + t258 * t285;
t286 = qJD(4) * t195;
t90 = -qJD(1) * t142 + t258 * t286;
t265 = qJD(4) * t80 + t255 * t90 + t257 * t88 - t363;
t79 = -t140 * t255 + t142 * t257;
t89 = qJD(1) * t141 + t256 * t285;
t91 = t256 * t286 + (t258 * t308 - t412) * qJD(1);
t266 = qJD(4) * t79 + t255 * t91 + t257 * t89 - t459;
t316 = (-t256 * t456 + t266 * t258) * t258 - (-t256 * t457 + t265 * t258) * t256;
t315 = (t266 * t256 + t258 * t456) * t258 - (t265 * t256 + t258 * t457) * t256;
t181 = t320 * qJD(4);
t294 = -qJ(3) * t261 - qJD(4) * t181;
t166 = t201 * t256;
t93 = qJD(4) * t166 + (t258 * t320 - t426) * qJD(1);
t32 = t185 * t201 + t234 + t294 * t258 + (-t159 - t93) * qJD(1) + t282 * qJDD(1) + t324;
t145 = t372 - t426;
t292 = rSges(5,1) * t343 - t356 * t428;
t92 = -qJD(1) * t144 + t292;
t33 = qJD(1) * t92 + qJDD(1) * t145 - t184 * t201 + (-pkin(6) * qJDD(1) + t294) * t256 + t272;
t314 = t33 * t256 + t32 * t258;
t313 = t256 * t44 + t258 * t43;
t51 = t256 * t300 + t402;
t312 = -t256 * t52 + t258 * t51;
t53 = t385 - t403;
t54 = -t139 * t256 + t384;
t311 = -t256 * t54 + t258 * t53;
t59 = t326 * t256 + (t145 + t463) * qJD(1) - t367;
t310 = t256 * t59 + t258 * t58;
t309 = -t256 * t93 - t258 * t92;
t65 = -t122 * t243 + t124 * t244;
t298 = -t144 * t256 - t145 * t258;
t295 = -t193 * t255 + t195 * t257;
t291 = -t48 + t393;
t287 = -t320 + t430;
t281 = -t319 - t435;
t279 = -qJD(1) * t303 + t147 * t189 - t148 * t397;
t277 = -t140 * t258 + t141 * t256;
t276 = t281 + t430;
t275 = (t255 * t375 + t257 * t376) * qJD(1);
t268 = t243 * t335 + t244 * t333 - t364;
t10 = -t256 * t455 + t268 * t258;
t269 = t243 * t336 + t244 * t334 - t461;
t11 = t269 * t256 + t258 * t454;
t12 = t268 * t256 + t258 * t455;
t47 = t256 * t302 + t389;
t22 = t189 * t47 - t397 * t48 + t388;
t49 = -t393 + t464;
t50 = -t121 * t256 + t386;
t64 = t258 * t297 - t147;
t62 = t64 * qJD(1);
t23 = t189 * t49 - t397 * t50 + t62;
t271 = t467 * qJD(1) + t470 * t189 - t469 * t397;
t262 = t271 * t243 + t244 * t450;
t267 = -qJD(1) * t170 + t243 * t329 + t244 * t330;
t28 = -t256 * t453 + t267 * t258;
t29 = t267 * t256 + t258 * t453;
t30 = -t243 * t334 + t244 * t336;
t31 = -t243 * t333 + t244 * t335;
t66 = -t123 * t243 + t125 * t244;
t9 = -t256 * t454 + t269 * t258;
t274 = (qJD(1) * t28 + qJDD(1) * t64 - t10 * t397 + t131 * t50 + t132 * t49 + t189 * t9) * t441 + (-t243 * t450 + t271 * t244) * t432 + t22 * t342 + t23 * t341 + (qJD(1) * t29 + qJDD(1) * t63 + t11 * t189 - t12 * t397 + t131 * t48 + t132 * t47) * t438 + (-t256 * t50 + t258 * t49) * t449 + (-t256 * t48 + t258 * t47) * t448 + (-t10 * t256 + t258 * t9 + (-t256 * t49 - t258 * t50) * qJD(1)) * t444 + (-t256 * t66 + t258 * t65) * t420 + (t11 * t258 - t12 * t256 + (-t256 * t47 - t258 * t48) * qJD(1)) * t442 + (-t256 * t31 + t258 * t30 + (-t256 * t65 - t258 * t66) * qJD(1)) * t431 + (-t256 * t279 + t258 * t262) * t445 + (t256 * t262 + t258 * t279) * t443;
t270 = -t389 + (-t121 - t302) * t256 + t386;
t178 = t306 * qJD(4);
t179 = t308 * qJD(4);
t264 = -qJD(1) * t191 + qJD(4) * t295 - t178 * t257 - t179 * t255;
t263 = t277 * t255 - t257 * t451;
t232 = rSges(4,2) * t360;
t223 = pkin(4) * t392;
t167 = t201 * t258;
t111 = t127 * t361;
t110 = qJD(1) * t156 - t242;
t109 = qJD(1) * t374 + t241;
t95 = t466 * qJD(1) - t327;
t94 = qJD(1) * t321 + t366;
t82 = t258 * t296 - t160;
t78 = t82 * qJD(1);
t77 = t298 * qJD(4);
t61 = qJD(1) * t368 + qJDD(1) * t204 + t347 - t351;
t60 = t374 * qJDD(1) + (-qJD(1) * t204 - t159) * qJD(1) + t370;
t46 = -t351 - t348 - qJDD(1) * t462 + qJD(1) * (-rSges(4,3) * t361 + t232) + t288;
t45 = t321 * qJDD(1) + (qJD(1) * t462 - t159) * qJD(1) + t280;
t37 = t264 * t256 + t258 * t452;
t36 = -t256 * t452 + t264 * t258;
t35 = -qJD(4) * t299 - t255 * t88 + t257 * t90;
t34 = -qJD(4) * t300 - t255 * t89 + t257 * t91;
t27 = qJD(4) * t311 + t78;
t26 = qJD(4) * t312 + t387;
t8 = t126 * t131 - t127 * t132 - t157 * t185 - t158 * t184 - t397 * t74 - t189 * t73 + (-t106 * t256 - t107 * t258) * qJD(4);
t1 = [(t62 + (-t291 - t48 + t464) * t189 - (t47 + t270) * t397) * t443 + (t78 + (t385 * t258 + (-t51 + (t139 + t300) * t256 - t384) * t256) * qJD(4)) * t338 - m(2) * (-g(1) * t200 + g(2) * t205) + (t66 + t64) * t449 + (t63 + t65) * t448 + (t80 + t82) * t447 + (t79 + t81) * t446 + (-t388 + (-t50 + t270) * t189 - (-t258 * t302 + t291 + t49) * t397 + t22) * t445 + (t31 + t28) * t444 + (t35 + t36) * t340 + (-t387 + ((t384 - t54 - t402) * t258 + (-t128 + t52 - t53 - t403) * t256) * qJD(4) + t26) * t339 + (-t296 * qJD(4) + t178 * t255 - t179 * t257 - t330 * t243 + t329 * t244) * qJD(1) + (t30 + t29 + t23) * t442 + (t34 + t37 + t27) * t337 + (t43 * t242 + t44 * (t293 + t346 + t373) + (-t176 * t253 - t322) * t422 + ((-t44 * rSges(6,3) + t276 * t43) * t258 + (t43 * (-qJ(2) - t419) + t44 * t276) * t256) * qJD(1) - (t283 * qJD(1) - t187 + t290 - t43) * t44 + (t15 - g(2)) * (t127 + t463 + t371) + (t14 - g(1)) * (-t222 + t247 + t419 * t258 + (-t319 + t430) * t256)) * m(6) + (t58 * t367 + t59 * (t292 + t346) - t326 * t421 + ((t287 * t58 + t436 * t59) * t258 + (t58 * (rSges(5,3) - qJ(2)) + t59 * t287) * t256) * qJD(1) - (qJD(1) * t289 + t168 + t345 - t58) * t59 + (t33 - g(2)) * (t256 * t436 + t372 + t463) + (t32 - g(1)) * (t287 * t256 + t436 * t258 + t247)) * m(5) + (t94 * t327 + t95 * (t232 + t346) + (t94 * t350 * t258 + (t94 * (-rSges(4,2) - qJ(2)) + t95 * t350) * t256) * qJD(1) - (qJD(1) * t331 + t345 - t94) * t95 + (t46 - g(2)) * t466 + (t45 - g(1)) * (t256 * t430 + t203 + t247)) * m(4) + (t109 * t242 + t110 * (t368 + t369) + (t109 * t437 * t258 + (t109 * (-rSges(3,3) - qJ(2)) - t110 * pkin(1)) * t256) * qJD(1) - (qJD(1) * t199 - t109 + t365) * t110 + (t61 - g(2)) * t156 + (t60 - g(1)) * (t256 * t437 + t247 + t427)) * m(3) + (t295 + m(2) * (t200 ^ 2 + t205 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1) - t172 * t243 + t174 * t244) * qJDD(1); (-m(3) + t355) * (g(1) * t256 - g(2) * t258) + 0.2e1 * (t14 * t440 + t15 * t439) * m(6) + 0.2e1 * (t32 * t440 + t33 * t439) * m(5) + 0.2e1 * (t439 * t46 + t440 * t45) * m(4) + 0.2e1 * (t439 * t61 + t440 * t60) * m(3); t355 * (g(1) * t258 + g(2) * t256) + m(4) * (t256 * t46 + t258 * t45) + m(5) * t314 + m(6) * t317; ((t160 * t356 - t362) * t258 + (t275 + (-t258 * t161 + t263) * qJD(4)) * t256) * t338 + ((t161 * t358 + t362) * t256 + (t275 + (-t256 * t160 + t263) * qJD(4)) * t258) * t339 + (-t256 * t35 + t258 * t34 + (-t79 * t256 - t258 * t80) * qJD(1)) * t431 + (qJD(1) * t37 + qJD(4) * t315 + qJDD(1) * t81 + t184 * t52 + t185 * t51) * t438 + (qJD(1) * t36 + qJD(4) * t316 + qJDD(1) * t82 + t184 * t54 + t185 * t53) * t441 + ((-t51 * t256 - t52 * t258) * qJD(1) + t315) * t337 + ((-t53 * t256 - t54 * t258) * qJD(1) + t316) * t340 + ((-t255 * t376 + t257 * t375) * qJD(1) + (t255 * t451 + t277 * t257) * qJD(4)) * t432 + t27 * t341 + t274 + t26 * t342 + (-t256 * t80 + t258 * t79) * t420 + t311 * t447 + t312 * t446 + (-(-t43 * t257 * t361 + (t42 * (-t256 ^ 2 - t258 ^ 2) * t257 - t313 * t255) * qJD(4)) * pkin(4) + t42 * t111 + t14 * t225 + t15 * t223 + (t14 * t176 + t43 * t323 - t8 * t380 - t42 * t417 + (t44 * t176 + t381 * t42) * qJD(1)) * t258 + (t15 * t176 + t44 * t323 + t8 * t381 + t42 * t418 + (t43 * (-pkin(4) * t257 - t176) + t42 * t157) * qJD(1)) * t256 - g(1) * (t154 + t225) - g(2) * (t153 + t223) - g(3) * t281 + t458) * m(6) + ((qJD(4) * t309 + t144 * t184 - t145 * t185) * t298 + t77 * ((-t144 * t258 + t145 * t256) * qJD(1) + t309) - t310 * t181 + ((t258 * t59 - t421) * qJD(1) + t314) * t201 - (-t166 * t58 + t167 * t59) * qJD(1) - (t77 * (-t166 * t256 - t167 * t258) - t310 * t320) * qJD(4) - g(1) * t167 - g(2) * t166 + g(3) * t320) * m(5); t274 + (t8 * (-t126 * t256 - t127 * t258) + t42 * (-t256 * t74 + t111 + (-t73 - t460) * t258) - t313 * t137 + ((t258 * t44 - t422) * qJD(1) + t317) * t176 - g(1) * t154 - g(2) * t153 + g(3) * t319 + t458) * m(6);];
tau = t1;
