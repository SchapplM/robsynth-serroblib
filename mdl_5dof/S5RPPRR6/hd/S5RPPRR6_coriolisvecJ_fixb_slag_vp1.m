% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:58:04
% DurationCPUTime: 18.36s
% Computational Cost: add. (19736->707), mult. (18541->970), div. (0->0), fcn. (17225->10), ass. (0->353)
t262 = qJ(1) + pkin(8);
t259 = cos(t262);
t265 = -pkin(6) - qJ(3);
t236 = t259 * t265;
t264 = cos(pkin(9));
t255 = pkin(3) * t264 + pkin(2);
t257 = sin(t262);
t395 = -t257 * t255 - t236;
t267 = sin(qJ(1));
t459 = pkin(1) * t267;
t493 = t395 - t459;
t261 = pkin(9) + qJ(4);
t258 = cos(t261);
t417 = t257 * t258;
t374 = rSges(5,1) * t417;
t492 = -t374 + t493;
t256 = sin(t261);
t266 = sin(qJ(5));
t413 = t259 * t266;
t268 = cos(qJ(5));
t415 = t257 * t268;
t184 = t258 * t415 - t413;
t170 = Icges(6,4) * t184;
t412 = t259 * t268;
t416 = t257 * t266;
t183 = t258 * t416 + t412;
t419 = t256 * t257;
t100 = -Icges(6,2) * t183 + Icges(6,6) * t419 + t170;
t169 = Icges(6,4) * t183;
t104 = -Icges(6,1) * t184 - Icges(6,5) * t419 + t169;
t490 = t100 * t266 + t104 * t268;
t97 = Icges(6,5) * t184 - Icges(6,6) * t183 + Icges(6,3) * t419;
t39 = -t256 * t490 - t258 * t97;
t334 = rSges(6,1) * t184 - rSges(6,2) * t183;
t106 = rSges(6,3) * t419 + t334;
t333 = rSges(6,1) * t268 - rSges(6,2) * t266;
t160 = -rSges(6,3) * t258 + t256 * t333;
t382 = qJD(5) * t256;
t383 = qJD(4) * t259;
t191 = -t257 * t382 + t383;
t458 = pkin(4) * t256;
t210 = -pkin(7) * t258 + t458;
t381 = qJD(5) * t258;
t233 = qJD(1) - t381;
t241 = qJD(3) * t257;
t491 = -t106 * t233 - t160 * t191 - t210 * t383 + t241;
t384 = qJD(4) * t257;
t190 = t259 * t382 + t384;
t25 = -t100 * t183 - t104 * t184 + t419 * t97;
t185 = -t258 * t413 + t415;
t418 = t256 * t259;
t186 = t258 * t412 + t416;
t434 = Icges(6,4) * t186;
t102 = Icges(6,2) * t185 + Icges(6,6) * t418 + t434;
t171 = Icges(6,4) * t185;
t105 = Icges(6,1) * t186 + Icges(6,5) * t418 + t171;
t99 = Icges(6,5) * t186 + Icges(6,6) * t185 + Icges(6,3) * t418;
t26 = -t183 * t102 + t184 * t105 + t99 * t419;
t321 = Icges(6,5) * t268 - Icges(6,6) * t266;
t154 = -Icges(6,3) * t258 + t256 * t321;
t432 = Icges(6,4) * t268;
t322 = -Icges(6,2) * t266 + t432;
t156 = -Icges(6,6) * t258 + t256 * t322;
t433 = Icges(6,4) * t266;
t324 = Icges(6,1) * t268 - t433;
t158 = -Icges(6,5) * t258 + t256 * t324;
t52 = t154 * t419 - t156 * t183 + t158 * t184;
t10 = t190 * t26 - t191 * t25 + t233 * t52;
t27 = t185 * t100 - t104 * t186 + t97 * t418;
t28 = t185 * t102 + t186 * t105 + t99 * t418;
t53 = t154 * t418 + t156 * t185 + t158 * t186;
t11 = t190 * t28 - t191 * t27 + t53 * t233;
t270 = qJD(1) ^ 2;
t449 = rSges(4,2) * sin(pkin(9));
t451 = rSges(4,1) * t264;
t311 = t257 * rSges(4,3) + (-t449 + t451) * t259;
t243 = t257 * qJ(3);
t209 = t259 * pkin(2) + t243;
t269 = cos(qJ(1));
t260 = t269 * pkin(1);
t353 = t209 + t260;
t485 = t311 + t353;
t483 = 0.2e1 * qJD(4);
t244 = t259 * qJ(3);
t206 = pkin(2) * t257 - t244;
t143 = t206 + t395;
t198 = qJD(1) * t206;
t481 = qJD(1) * t143 - t198;
t480 = -rSges(5,2) * t419 - t259 * rSges(5,3);
t248 = Icges(5,4) * t258;
t323 = -Icges(5,2) * t256 + t248;
t203 = Icges(5,1) * t256 + t248;
t351 = t259 * rSges(3,1) - rSges(3,2) * t257;
t479 = t260 + t351;
t200 = Icges(5,5) * t258 - Icges(5,6) * t256;
t199 = Icges(5,5) * t256 + Icges(5,6) * t258;
t298 = qJD(4) * t199;
t435 = Icges(5,4) * t256;
t204 = Icges(5,1) * t258 - t435;
t150 = Icges(5,5) * t257 + t204 * t259;
t148 = Icges(5,6) * t257 + t259 * t323;
t422 = t148 * t256;
t316 = -t150 * t258 + t422;
t426 = Icges(5,3) * t259;
t476 = -t259 * t298 + (-t200 * t257 + t316 + t426) * qJD(1);
t219 = Icges(5,4) * t419;
t431 = Icges(5,5) * t259;
t149 = Icges(5,1) * t417 - t219 - t431;
t428 = Icges(5,6) * t259;
t147 = Icges(5,4) * t417 - Icges(5,2) * t419 - t428;
t423 = t147 * t256;
t317 = -t149 * t258 + t423;
t146 = Icges(5,3) * t257 + t200 * t259;
t390 = qJD(1) * t146;
t475 = qJD(1) * t317 - t257 * t298 + t390;
t145 = Icges(5,5) * t417 - Icges(5,6) * t419 - t426;
t56 = -t145 * t259 - t257 * t317;
t201 = Icges(5,2) * t258 + t435;
t313 = t201 * t256 - t203 * t258;
t474 = t313 * qJD(1) + t200 * qJD(4);
t473 = t257 * (-t201 * t259 + t150) - t259 * (-Icges(5,2) * t417 + t149 - t219);
t155 = Icges(6,3) * t256 + t258 * t321;
t314 = -t156 * t266 + t158 * t268;
t319 = -t102 * t266 + t105 * t268;
t472 = t190 * (-t154 * t259 - t319) - t191 * (-t154 * t257 + t490) + t233 * (t155 - t314);
t188 = (-Icges(6,2) * t268 - t433) * t256;
t471 = t190 * (-Icges(6,2) * t186 + t105 + t171) - t191 * (-Icges(6,2) * t184 - t104 - t169) + t233 * (t158 + t188);
t379 = qJD(4) * qJD(5);
t363 = t258 * t379;
t141 = qJD(1) * t190 + t257 * t363;
t470 = t141 / 0.2e1;
t142 = qJD(1) * t191 + t259 * t363;
t469 = t142 / 0.2e1;
t468 = -t190 / 0.2e1;
t467 = t190 / 0.2e1;
t466 = -t191 / 0.2e1;
t465 = t191 / 0.2e1;
t464 = -t233 / 0.2e1;
t463 = t233 / 0.2e1;
t462 = t257 / 0.2e1;
t461 = -t259 / 0.2e1;
t460 = -rSges(6,3) - pkin(7);
t457 = pkin(4) * t258;
t368 = t258 * t384;
t386 = qJD(1) * t259;
t294 = t256 * t386 + t368;
t385 = qJD(4) * t256;
t287 = t233 * t268 + t266 * t385;
t387 = qJD(1) * t258;
t347 = -qJD(5) + t387;
t87 = t257 * t287 - t347 * t413;
t286 = t233 * t266 - t268 * t385;
t88 = t257 * t286 + t347 * t412;
t44 = Icges(6,5) * t88 + Icges(6,6) * t87 + Icges(6,3) * t294;
t46 = Icges(6,4) * t88 + Icges(6,2) * t87 + Icges(6,6) * t294;
t48 = Icges(6,1) * t88 + Icges(6,4) * t87 + Icges(6,5) * t294;
t7 = (-qJD(4) * t490 - t44) * t258 + (qJD(4) * t97 - t266 * t46 + t268 * t48 + (-t100 * t268 + t104 * t266) * qJD(5)) * t256;
t456 = t7 * t191;
t367 = t258 * t383;
t388 = qJD(1) * t257;
t372 = t256 * t388;
t293 = t367 - t372;
t85 = t259 * t287 + t347 * t416;
t86 = t259 * t286 - t347 * t415;
t43 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t293;
t45 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t293;
t47 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t293;
t8 = (qJD(4) * t319 - t43) * t258 + (qJD(4) * t99 - t266 * t45 + t268 * t47 + (-t102 * t268 - t105 * t266) * qJD(5)) * t256;
t455 = t8 * t190;
t454 = qJD(1) / 0.2e1;
t453 = pkin(2) - t255;
t187 = (-Icges(6,5) * t266 - Icges(6,6) * t268) * t256;
t111 = qJD(4) * t155 + qJD(5) * t187;
t157 = Icges(6,6) * t256 + t258 * t322;
t112 = qJD(4) * t157 + qJD(5) * t188;
t159 = Icges(6,5) * t256 + t258 * t324;
t189 = (-Icges(6,1) * t266 - t432) * t256;
t113 = qJD(4) * t159 + qJD(5) * t189;
t22 = (qJD(4) * t314 - t111) * t258 + (qJD(4) * t154 - t112 * t266 + t113 * t268 + (-t156 * t268 - t158 * t266) * qJD(5)) * t256;
t364 = t256 * t379;
t63 = -t154 * t258 + t256 * t314;
t452 = t22 * t233 + t63 * t364;
t448 = rSges(6,3) * t256;
t108 = t186 * rSges(6,1) + t185 * rSges(6,2) + rSges(6,3) * t418;
t214 = pkin(7) * t367;
t369 = t256 * t383;
t295 = -t257 * t387 - t369;
t114 = pkin(4) * t295 - pkin(7) * t372 + t214;
t161 = t258 * t333 + t448;
t192 = (-rSges(6,1) * t266 - rSges(6,2) * t268) * t256;
t116 = qJD(4) * t161 + qJD(5) * t192;
t211 = pkin(7) * t256 + t457;
t197 = qJD(4) * t211;
t237 = qJ(3) * t386;
t378 = t270 * t459;
t380 = qJD(1) * qJD(3);
t391 = t237 + t241;
t312 = qJD(1) * (-pkin(2) * t388 + t391) + t257 * t380 - t378;
t304 = qJD(1) * (-t237 + (t257 * t453 - t236) * qJD(1)) + t312;
t376 = t86 * rSges(6,1) + t85 * rSges(6,2) + rSges(6,3) * t367;
t49 = -rSges(6,3) * t372 + t376;
t13 = qJD(1) * t114 - t116 * t190 - t142 * t160 + t233 * t49 + (t108 * t382 - t197 * t257 - t210 * t386) * qJD(4) + t304;
t446 = t13 * t259;
t115 = t294 * pkin(7) + (-t256 * t384 + t258 * t386) * pkin(4);
t377 = t270 * t260;
t346 = t259 * t380 - t377;
t370 = t210 * t384;
t242 = qJD(3) * t259;
t172 = qJD(1) * t209 - t242;
t228 = t265 * t388;
t403 = t228 - (-t259 * t453 - t243) * qJD(1) - t172;
t340 = rSges(6,1) * t88 + rSges(6,2) * t87;
t50 = rSges(6,3) * t294 + t340;
t14 = -t116 * t191 + t141 * t160 - t233 * t50 + (-t106 * t382 - t197 * t259) * qJD(4) + (-t115 + t370 + t403) * qJD(1) + t346;
t445 = t14 * t257;
t205 = rSges(5,1) * t256 + rSges(5,2) * t258;
t174 = t205 * t259;
t249 = t257 * rSges(5,3);
t414 = t258 * t259;
t152 = rSges(5,1) * t414 - rSges(5,2) * t418 + t249;
t221 = t259 * t255;
t348 = -t257 * t265 + t221;
t344 = t348 - t209 + t353;
t371 = t205 * t384;
t61 = -t371 - t242 + (t152 + t344) * qJD(1);
t444 = t174 * t61;
t180 = pkin(4) * t414 + pkin(7) * t418;
t32 = -t370 + t108 * t233 - t160 * t190 - t242 + (t180 + t344) * qJD(1);
t440 = t257 * t32;
t151 = t374 + t480;
t366 = t205 * t383;
t339 = t241 - t366;
t354 = -t206 - t459;
t345 = t143 + t354;
t60 = (-t151 + t345) * qJD(1) + t339;
t439 = t257 * t60;
t438 = t39 * t141;
t40 = t256 * t319 - t258 * t99;
t437 = t40 * t142;
t421 = t199 * t257;
t420 = t199 * t259;
t79 = -t257 * t313 - t420;
t411 = t79 * qJD(1);
t178 = t211 * t257;
t408 = t106 + t178;
t407 = t108 + t180;
t406 = -t116 - t197;
t405 = -t257 * t145 - t149 * t414;
t404 = t257 * t146 + t150 * t414;
t399 = t160 + t210;
t398 = -t201 + t204;
t397 = t203 + t323;
t396 = rSges(5,2) * t372 + rSges(5,3) * t386;
t229 = t257 * t449;
t394 = rSges(4,3) * t386 + qJD(1) * t229;
t393 = t228 + t242;
t392 = t259 * rSges(4,3) + t229;
t389 = qJD(1) * t200;
t375 = t257 * t451;
t365 = -pkin(2) - t451;
t361 = t386 / 0.2e1;
t360 = t385 / 0.2e1;
t359 = -t384 / 0.2e1;
t358 = t384 / 0.2e1;
t356 = t383 / 0.2e1;
t126 = t150 * t417;
t350 = t146 * t259 - t126;
t349 = -t145 + t422;
t341 = qJD(5) * t360;
t207 = rSges(3,1) * t257 + rSges(3,2) * t259;
t335 = rSges(5,1) * t258 - rSges(5,2) * t256;
t332 = t25 * t259 - t257 * t26;
t331 = t25 * t257 + t259 * t26;
t330 = t257 * t28 - t259 * t27;
t329 = t257 * t27 + t259 * t28;
t328 = t257 * t40 - t259 * t39;
t327 = t257 * t39 + t259 * t40;
t326 = -t257 * t61 - t259 * t60;
t318 = t106 * t259 - t108 * t257;
t71 = t147 * t258 + t149 * t256;
t72 = t148 * t258 + t150 * t256;
t315 = t151 * t257 + t152 * t259;
t173 = t205 * t257;
t57 = -t148 * t419 - t350;
t302 = (t257 * t57 - t259 * t56) * qJD(4);
t58 = -t147 * t418 - t405;
t59 = -t148 * t418 + t404;
t301 = (t257 * t59 - t259 * t58) * qJD(4);
t300 = qJD(4) * t203;
t299 = qJD(4) * t201;
t297 = -t211 - t448;
t292 = t154 * t233 + t190 * t99 - t191 * t97;
t291 = (-Icges(6,5) * t183 - Icges(6,6) * t184) * t191 - (Icges(6,5) * t185 - Icges(6,6) * t186) * t190 - t187 * t233;
t290 = t147 * t259 - t148 * t257;
t289 = t256 * t291;
t282 = (-t256 * t397 + t258 * t398) * qJD(1);
t95 = rSges(5,1) * t295 - rSges(5,2) * t367 + t396;
t96 = -qJD(4) * t173 + (t259 * t335 + t249) * qJD(1);
t281 = t257 * t96 + t259 * t95 + (t151 * t259 - t152 * t257) * qJD(1);
t280 = (Icges(6,1) * t185 - t102 - t434) * t190 - (-Icges(6,1) * t183 - t100 - t170) * t191 + (-t156 + t189) * t233;
t31 = (-t178 + t345) * qJD(1) + t491;
t33 = t106 * t190 + t108 * t191 + qJD(2) + (t178 * t257 + t180 * t259) * qJD(4);
t276 = t33 * t318 + (t257 * t31 - t259 * t32) * t160;
t92 = qJD(1) * t148 - t257 * t299;
t94 = qJD(1) * t150 - t257 * t300;
t275 = qJD(1) * t145 - qJD(4) * t71 - t256 * t92 + t258 * t94;
t91 = -t259 * t299 + (-t257 * t323 + t428) * qJD(1);
t93 = -t259 * t300 + (-t204 * t257 + t431) * qJD(1);
t274 = -qJD(4) * t72 - t256 * t91 + t258 * t93 + t390;
t194 = t323 * qJD(4);
t195 = t204 * qJD(4);
t273 = qJD(1) * t199 - t194 * t256 + t195 * t258 + (-t201 * t258 - t203 * t256) * qJD(4);
t272 = -t256 * t473 + t290 * t258;
t271 = t472 * t256;
t196 = t335 * qJD(4);
t179 = t210 * t259;
t177 = t210 * t257;
t153 = t375 - t392;
t136 = t160 * t259;
t135 = t160 * t257;
t134 = t158 * t259;
t133 = t158 * t257;
t132 = t156 * t259;
t131 = t156 * t257;
t124 = rSges(6,1) * t185 - rSges(6,2) * t186;
t123 = -rSges(6,1) * t183 - rSges(6,2) * t184;
t110 = qJD(1) * t485 - t242;
t109 = t241 + (-t153 + t354) * qJD(1);
t80 = -t259 * t313 + t421;
t75 = t80 * qJD(1);
t74 = (-qJD(1) * t311 - t172) * qJD(1) + t346;
t73 = qJD(1) * (-qJD(1) * t375 + t394) + t312;
t70 = qJD(4) * t315 + qJD(2);
t42 = -t196 * t383 + (-t96 + t371 + t403) * qJD(1) + t346;
t41 = -t196 * t384 + (t95 - t366) * qJD(1) + t304;
t38 = t273 * t257 - t259 * t474;
t37 = t257 * t474 + t273 * t259;
t36 = -qJD(4) * t316 + t256 * t93 + t258 * t91;
t35 = -t317 * qJD(4) + t256 * t94 + t258 * t92;
t34 = t281 * qJD(4);
t24 = t75 + t301;
t23 = t302 + t411;
t16 = t111 * t419 - t112 * t183 + t113 * t184 + t154 * t294 + t156 * t87 + t158 * t88;
t15 = t111 * t418 + t112 * t185 + t113 * t186 + t154 * t293 + t156 * t85 + t158 * t86;
t12 = t190 * t40 - t191 * t39 + t233 * t63;
t9 = t106 * t142 - t108 * t141 + t190 * t50 + t191 * t49 + (t114 * t259 + t115 * t257 + (t178 * t259 - t180 * t257) * qJD(1)) * qJD(4);
t6 = t99 * t368 + t102 * t87 + t105 * t88 - t183 * t45 + t184 * t47 + (t257 * t43 + t386 * t99) * t256;
t5 = t97 * t368 + t100 * t87 - t104 * t88 - t183 * t46 + t184 * t48 + (t257 * t44 + t386 * t97) * t256;
t4 = t99 * t367 + t102 * t85 + t105 * t86 + t185 * t45 + t186 * t47 + (t259 * t43 - t388 * t99) * t256;
t3 = t97 * t367 + t100 * t85 - t104 * t86 + t185 * t46 + t186 * t48 + (t259 * t44 - t388 * t97) * t256;
t2 = t141 * t25 + t142 * t26 + t16 * t233 + t190 * t6 - t191 * t5 + t364 * t52;
t1 = t141 * t27 + t142 * t28 + t15 * t233 + t190 * t4 - t191 * t3 + t364 * t53;
t17 = [(-qJD(4) * t313 + t194 * t258 + t195 * t256) * qJD(1) + t452 + t16 * t466 + t15 * t467 + t53 * t469 + t52 * t470 + t437 / 0.2e1 + t438 / 0.2e1 + t455 / 0.2e1 - t456 / 0.2e1 + m(3) * ((-t207 * t270 - t378) * t479 + (-t377 + (-0.2e1 * t351 - t260 + t479) * t270) * (-t207 - t459)) + (t75 + ((t57 - t126 + (t146 + t423) * t259 + t405) * t259 + t404 * t257) * qJD(4)) * t356 + (t23 - t411 + ((t259 * t349 - t404 + t59) * t259 + (t257 * t349 + t350 + t58) * t257) * qJD(4)) * t359 + (t36 + t37) * t358 + (t466 + t465) * t11 + (-(-t31 + (-t178 - t459) * qJD(1) + t481 + t491) * t32 + t14 * (-t334 + t493) + t31 * (-t340 + t393) + t13 * (t221 + t260 + t407) + t32 * (-pkin(4) * t369 + t214 + t241 + t376) + (t14 * t297 - t13 * t265 + t31 * (t258 * t460 + t458) * qJD(4)) * t257 + ((-t267 * t32 - t269 * t31) * pkin(1) + (t31 * (-t255 + t297) - t32 * t265) * t259 + (t256 * t460 - t255 - t457) * t440) * qJD(1)) * m(6) + (t42 * (-t480 + t492) + t41 * (t152 + t260 + t348) + (t205 * t439 - t444) * qJD(4) + (t393 + (-t249 - t260 + (-t255 - t335) * t259) * qJD(1)) * t60 + (t60 - t339 - t481 + t241 + t396 + (t151 + t459 + t492) * qJD(1)) * t61) * m(5) + (-(-t109 - t198 + t241 + (-t153 - t459) * qJD(1)) * t110 + t74 * (t257 * t365 + t244 + t392 - t459) + t109 * t242 + t73 * t485 + t110 * (t391 + t394) + ((-t109 * t269 - t110 * t267) * pkin(1) + t109 * (t365 + t449) * t259 + (t109 * (-rSges(4,3) - qJ(3)) + t110 * t365) * t257) * qJD(1)) * m(4) - (t35 + t38 + t24) * t383 / 0.2e1 + ((t71 + t79) * t257 + (t72 + t80) * t259) * qJD(4) * t454; m(5) * t34 + m(6) * t9; 0.2e1 * (-t446 / 0.2e1 + t445 / 0.2e1) * m(6) + 0.2e1 * (t41 * t461 + t42 * t462) * m(5) + 0.2e1 * (t461 * t73 + t462 * t74) * m(4); -qJD(1) * ((t256 * t398 + t258 * t397) * qJD(1) + (t290 * t256 + t258 * t473) * qJD(4)) / 0.2e1 + ((-t383 * t421 - t389) * t259 + (t282 + (t259 * t420 + t272) * qJD(4)) * t257) * t356 + ((-t384 * t420 + t389) * t257 + (t282 + (t257 * t421 + t272) * qJD(4)) * t259) * t359 - t12 * t382 / 0.2e1 + ((t132 * t183 - t134 * t184) * t190 - (t131 * t183 - t133 * t184) * t191 + (-t157 * t183 + t159 * t184) * t233 + (t256 * t52 + t26 * t414) * qJD(5) + ((qJD(5) * t25 + t292) * t258 + t271) * t257) * t465 + (qJD(1) * t331 + t257 * t6 - t259 * t5) * t466 + (qJD(1) * t329 + t257 * t4 - t259 * t3) * t467 + ((-t132 * t185 - t134 * t186) * t190 - (-t131 * t185 - t133 * t186) * t191 + (t157 * t185 + t159 * t186) * t233 + (t256 * t53 + t27 * t417) * qJD(5) + ((qJD(5) * t28 + t292) * t258 + t271) * t259) * t468 + t330 * t469 + (t257 * t36 - t259 * t35 + (t71 * t257 + t72 * t259) * qJD(1)) * t454 + (qJD(1) * t327 + t257 * t8 - t259 * t7) * t463 + (((t132 * t266 - t134 * t268 + t99) * t190 - (t131 * t266 - t133 * t268 + t97) * t191 + (-t157 * t266 + t159 * t268 + t154) * t233 + t63 * qJD(5)) * t256 + (qJD(5) * t327 - t472) * t258) * t464 + t328 * t341 - t141 * t332 / 0.2e1 - (t257 * t10 + t259 * t11) * t381 / 0.2e1 + ((-t14 * t399 + t31 * t406 + t9 * t407 + t33 * (t114 + t49) + (-t32 * t399 + t33 * t408) * qJD(1)) * t259 + (-t13 * t399 + t32 * t406 + t9 * t408 + t33 * (t115 + t50) + (t31 * t399 - t33 * t407) * qJD(1)) * t257 - t31 * (qJD(1) * t177 + t135 * t233 - t161 * t191 - t211 * t383) - t32 * (-qJD(1) * t179 - t136 * t233 - t161 * t190 - t211 * t384) - t33 * (-t135 * t190 - t136 * t191 - t177 * t384 - t179 * t383) - ((-t106 * t31 + t108 * t32) * t256 + t276 * t258) * qJD(5)) * m(6) + (t34 * t315 + t70 * t281 + t326 * t196 + (-t41 * t257 - t42 * t259 + (-t259 * t61 + t439) * qJD(1)) * t205 - (t173 * t60 - t444) * qJD(1) - (t70 * (-t173 * t257 - t174 * t259) + t326 * t335) * qJD(4)) * m(5) + (qJD(1) * t37 + t1 + (-(t257 * t475 + t275 * t259) * t259 + (t257 * t476 + t274 * t259) * t257 + (t58 * t257 + t59 * t259) * qJD(1)) * t483) * t462 + (qJD(1) * t38 + t2 + (-(t275 * t257 - t259 * t475) * t259 + (t274 * t257 - t259 * t476) * t257 + (t56 * t257 + t57 * t259) * qJD(1)) * t483) * t461 + (t23 + t10 + t302) * t388 / 0.2e1 + (t301 + t24 + t11) * t361; t1 * t418 / 0.2e1 + (t256 * t329 - t258 * t53) * t469 + ((qJD(4) * t329 - t15) * t258 + (-qJD(1) * t330 + qJD(4) * t53 + t257 * t3 + t259 * t4) * t256) * t467 + t2 * t419 / 0.2e1 + (t256 * t331 - t258 * t52) * t470 + ((qJD(4) * t331 - t16) * t258 + (qJD(1) * t332 + qJD(4) * t52 + t257 * t5 + t259 * t6) * t256) * t466 + t12 * t360 - t258 * (t437 + t438 + t452 + t455 - t456) / 0.2e1 + (t256 * t327 - t258 * t63) * t341 + ((qJD(4) * t327 - t22) * t258 + (-qJD(1) * t328 + qJD(4) * t63 + t257 * t7 + t259 * t8) * t256) * t463 + (t185 * t471 + t280 * t186 - t259 * t289) * t468 + (-t183 * t471 + t184 * t280 - t257 * t289) * t465 + (t291 * t258 + (-t266 * t471 + t268 * t280) * t256) * t464 + (-t372 / 0.2e1 + t258 * t356) * t11 + (t256 * t361 + t258 * t358) * t10 + ((qJD(4) * t276 + t14 * t106 - t13 * t108 + t31 * t50 - t32 * t49) * t258 + (t31 * (-qJD(4) * t106 + t116 * t257) + t32 * (qJD(4) * t108 - t116 * t259) + t9 * t318 + t33 * (-t106 * t388 - t108 * t386 - t257 * t49 + t259 * t50) + (-t446 + t445 + (t259 * t31 + t440) * qJD(1)) * t160) * t256 - t31 * (-t123 * t233 - t191 * t192) - t32 * (t124 * t233 - t190 * t192) - t33 * (t123 * t190 + t124 * t191)) * m(6);];
tauc = t17(:);
