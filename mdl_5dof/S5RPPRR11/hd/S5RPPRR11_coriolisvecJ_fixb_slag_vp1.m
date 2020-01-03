% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR11_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:32
% EndTime: 2019-12-31 18:05:50
% DurationCPUTime: 15.06s
% Computational Cost: add. (7233->716), mult. (18939->998), div. (0->0), fcn. (17421->6), ass. (0->336)
t255 = sin(qJ(1));
t257 = cos(qJ(4));
t258 = cos(qJ(1));
t254 = sin(qJ(4));
t256 = cos(qJ(5));
t388 = t258 * t256;
t253 = sin(qJ(5));
t393 = t255 * t253;
t170 = t254 * t393 - t388;
t392 = t255 * t256;
t171 = t253 * t258 + t254 * t392;
t165 = Icges(6,4) * t171;
t391 = t255 * t257;
t91 = -Icges(6,2) * t170 - Icges(6,6) * t391 + t165;
t164 = Icges(6,4) * t170;
t95 = -Icges(6,1) * t171 + Icges(6,5) * t391 + t164;
t318 = t170 * t91 + t171 * t95;
t394 = t254 * t258;
t172 = -t253 * t394 - t392;
t173 = t254 * t388 - t393;
t389 = t257 * t258;
t405 = Icges(6,4) * t173;
t93 = Icges(6,2) * t172 - Icges(6,6) * t389 + t405;
t166 = Icges(6,4) * t172;
t96 = Icges(6,1) * t173 - Icges(6,5) * t389 + t166;
t424 = t172 * t93 + t173 * t96;
t89 = -Icges(6,5) * t171 + Icges(6,6) * t170 + Icges(6,3) * t391;
t90 = Icges(6,5) * t173 + Icges(6,6) * t172 - Icges(6,3) * t389;
t471 = t318 + t424 + (-t255 * t89 - t258 * t90) * t257;
t425 = t172 * t91 - t173 * t95;
t426 = -t170 * t93 + t171 * t96;
t470 = t425 + t257 * (-t255 * t90 + t258 * t89) + t426;
t302 = Icges(5,5) * t254 + Icges(5,6) * t257;
t150 = -Icges(5,3) * t255 + t258 * t302;
t407 = Icges(5,4) * t254;
t304 = Icges(5,2) * t257 + t407;
t153 = Icges(5,6) * t258 + t255 * t304;
t406 = Icges(5,4) * t257;
t306 = Icges(5,1) * t254 + t406;
t157 = Icges(5,5) * t258 + t255 * t306;
t299 = t153 * t257 + t157 * t254;
t469 = t150 + t299;
t315 = t253 * t91 + t256 * t95;
t35 = -t254 * t89 - t257 * t315;
t466 = t255 * t258;
t360 = qJD(5) * t254;
t228 = qJD(1) + t360;
t301 = Icges(6,5) * t256 - Icges(6,6) * t253;
t148 = Icges(6,3) * t254 + t257 * t301;
t403 = Icges(6,4) * t256;
t303 = -Icges(6,2) * t253 + t403;
t152 = Icges(6,6) * t254 + t257 * t303;
t404 = Icges(6,4) * t253;
t305 = Icges(6,1) * t256 - t404;
t156 = Icges(6,5) * t254 + t257 * t305;
t273 = t148 * t391 + t152 * t170 - t156 * t171;
t465 = t273 * t228;
t361 = qJD(4) * t258;
t344 = t254 * t361;
t366 = qJD(1) * t255;
t346 = t257 * t366;
t277 = t344 + t346;
t368 = t255 * rSges(4,2) + t258 * rSges(4,3);
t212 = t258 * pkin(1) + t255 * qJ(2);
t458 = t258 * qJ(3) + t212;
t464 = t458 + t368;
t320 = rSges(6,1) * t256 - rSges(6,2) * t253;
t161 = rSges(6,3) * t254 + t257 * t320;
t432 = t257 * pkin(4);
t216 = pkin(7) * t254 + t432;
t359 = qJD(5) * t257;
t286 = t255 * t359 - t361;
t97 = rSges(6,1) * t171 - rSges(6,2) * t170 - rSges(6,3) * t391;
t463 = -t161 * t286 + t216 * t361 - t228 * t97;
t460 = 0.2e1 * qJD(4);
t330 = -rSges(3,2) * t258 + t255 * rSges(3,3);
t459 = t212 + t330;
t149 = Icges(5,3) * t258 + t255 * t302;
t321 = rSges(5,1) * t254 + rSges(5,2) * t257;
t160 = rSges(5,3) * t258 + t255 * t321;
t327 = qJD(1) * t254 + qJD(5);
t342 = t257 * t361;
t455 = t255 * t327 - t342;
t202 = Icges(5,5) * t257 - Icges(5,6) * t254;
t279 = qJD(4) * t202;
t154 = -Icges(5,6) * t255 + t258 * t304;
t224 = Icges(5,4) * t389;
t402 = Icges(5,5) * t255;
t158 = Icges(5,1) * t394 + t224 - t402;
t298 = t154 * t257 + t158 * t254;
t454 = t258 * t279 + (-t149 + t298) * qJD(1);
t453 = t469 * qJD(1) + t255 * t279;
t204 = -Icges(5,2) * t254 + t406;
t206 = Icges(5,1) * t257 - t407;
t297 = t204 * t257 + t206 * t254;
t452 = t297 * qJD(1) - t302 * qJD(4);
t451 = t255 * (-Icges(5,2) * t394 + t158 + t224) - t258 * (t204 * t255 + t157);
t147 = Icges(6,3) * t257 - t254 * t301;
t363 = qJD(4) * t255;
t192 = -t258 * t359 - t363;
t300 = t152 * t253 - t156 * t256;
t314 = t253 * t93 - t256 * t96;
t450 = -t192 * (t148 * t258 + t314) + t286 * (t148 * t255 + t315) - t228 * (t147 + t300);
t177 = (-Icges(6,2) * t256 - t404) * t257;
t449 = t192 * (-Icges(6,2) * t173 + t166 + t96) - t286 * (-Icges(6,2) * t171 - t164 - t95) + t228 * (t156 + t177);
t355 = qJD(4) * qJD(5);
t341 = t254 * t355;
t144 = qJD(1) * t192 + t255 * t341;
t448 = t144 / 0.2e1;
t145 = qJD(1) * t286 + t258 * t341;
t447 = t145 / 0.2e1;
t446 = -t192 / 0.2e1;
t445 = t192 / 0.2e1;
t444 = t286 / 0.2e1;
t443 = -t286 / 0.2e1;
t442 = -t228 / 0.2e1;
t441 = t228 / 0.2e1;
t439 = t255 / 0.2e1;
t438 = -t258 / 0.2e1;
t436 = rSges(3,2) - pkin(1);
t435 = -rSges(5,3) - pkin(6);
t434 = pkin(4) * t254;
t433 = pkin(6) * t258;
t365 = qJD(1) * t258;
t345 = t257 * t365;
t278 = t254 * t363 - t345;
t294 = t228 * t256;
t362 = qJD(4) * t257;
t343 = t255 * t362;
t81 = -t255 * t294 + (-t258 * t327 - t343) * t253;
t295 = t228 * t253;
t82 = t327 * t388 + (t256 * t362 - t295) * t255;
t44 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t278;
t46 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t278;
t48 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t278;
t7 = (qJD(4) * t315 + t44) * t254 + (-qJD(4) * t89 - t253 * t46 + t256 * t48 + (t253 * t95 - t256 * t91) * qJD(5)) * t257;
t431 = t7 * t286;
t79 = t253 * t455 - t258 * t294;
t80 = -t256 * t455 - t258 * t295;
t43 = Icges(6,5) * t80 + Icges(6,6) * t79 + Icges(6,3) * t277;
t45 = Icges(6,4) * t80 + Icges(6,2) * t79 + Icges(6,6) * t277;
t47 = Icges(6,1) * t80 + Icges(6,4) * t79 + Icges(6,5) * t277;
t8 = (qJD(4) * t314 + t43) * t254 + (qJD(4) * t90 - t253 * t45 + t256 * t47 + (-t253 * t96 - t256 * t93) * qJD(5)) * t257;
t430 = t8 * t192;
t429 = -qJD(1) / 0.2e1;
t428 = -pkin(1) - qJ(3);
t174 = (-Icges(6,5) * t253 - Icges(6,6) * t256) * t257;
t100 = qJD(4) * t147 + qJD(5) * t174;
t151 = Icges(6,6) * t257 - t254 * t303;
t103 = qJD(4) * t151 + qJD(5) * t177;
t155 = Icges(6,5) * t257 - t254 * t305;
t180 = (-Icges(6,1) * t253 - t403) * t257;
t106 = qJD(4) * t155 + qJD(5) * t180;
t18 = (qJD(4) * t300 + t100) * t254 + (qJD(4) * t148 - t103 * t253 + t106 * t256 + (-t152 * t256 - t156 * t253) * qJD(5)) * t257;
t340 = t257 * t355;
t57 = t148 * t254 - t257 * t300;
t427 = t18 * t228 + t57 * t340;
t422 = rSges(4,2) * t258;
t421 = rSges(3,3) * t258;
t420 = rSges(5,3) * t255;
t418 = rSges(6,3) * t257;
t159 = -t254 * t320 + t418;
t184 = (-rSges(6,1) * t253 - rSges(6,2) * t256) * t257;
t111 = qJD(4) * t159 + qJD(5) * t184;
t350 = pkin(4) * t342 + pkin(7) * t277;
t395 = t254 * t255;
t353 = pkin(4) * t395;
t124 = -qJD(1) * t353 + t350;
t215 = pkin(7) * t257 - t434;
t199 = qJD(4) * t215;
t259 = qJD(1) ^ 2;
t398 = qJ(3) * t255;
t319 = -t398 - t433;
t356 = qJD(1) * qJD(3);
t357 = qJD(1) * qJD(2);
t235 = qJ(2) * t365;
t244 = qJD(2) * t255;
t373 = t235 + t244;
t378 = qJD(1) * (-pkin(1) * t366 + t373) + t255 * t357;
t351 = 0.2e1 * t258 * t356 + t378;
t49 = t80 * rSges(6,1) + t79 * rSges(6,2) + rSges(6,3) * t277;
t99 = t173 * rSges(6,1) + t172 * rSges(6,2) - rSges(6,3) * t389;
t13 = qJD(1) * t124 - t111 * t192 - t145 * t161 + t228 * t49 + t319 * t259 + (t199 * t255 + t216 * t365 + t359 * t99) * qJD(4) + t351;
t417 = t13 * t258;
t125 = t278 * pkin(7) + (t254 * t365 + t343) * pkin(4);
t245 = qJD(2) * t258;
t167 = qJD(1) * t212 - t245;
t390 = t255 * t259;
t241 = pkin(6) * t390;
t234 = t258 * t357;
t374 = -0.2e1 * t255 * t356 + t234;
t387 = t258 * t259;
t293 = -qJ(3) * t387 + t374;
t323 = rSges(6,1) * t82 + rSges(6,2) * t81;
t50 = rSges(6,3) * t278 + t323;
t14 = -t111 * t286 + t144 * t161 - t228 * t50 + t241 + (t199 * t258 - t359 * t97) * qJD(4) + (-t216 * t363 - t125 - t167) * qJD(1) + t293;
t416 = t14 * t255;
t413 = t35 * t144;
t36 = t254 * t90 - t257 * t314;
t412 = t36 * t145;
t231 = pkin(7) * t391;
t186 = -t231 + t353;
t409 = -t186 - t97;
t188 = pkin(4) * t394 - pkin(7) * t389;
t408 = t188 + t99;
t397 = t149 * t255;
t396 = t149 * t258;
t175 = t202 * t255;
t176 = t202 * t258;
t146 = t258 * t150;
t83 = t255 * t297 + t176;
t386 = t83 * qJD(1);
t385 = t111 + t199;
t384 = t153 * t389 + t157 * t394;
t383 = t154 * t389 + t158 * t394;
t379 = t161 + t216;
t377 = -t304 + t206;
t376 = -t204 - t306;
t375 = rSges(5,1) * t394 + rSges(5,2) * t389;
t372 = rSges(3,2) * t366 + rSges(3,3) * t365;
t371 = pkin(6) * t366 + t245;
t243 = qJD(3) * t258;
t370 = t243 + t244;
t248 = t258 * qJ(2);
t208 = pkin(1) * t255 - t248;
t369 = -qJD(1) * t208 + t244;
t364 = qJD(4) * t254;
t358 = t302 * qJD(1);
t354 = -rSges(4,3) + t428;
t352 = rSges(5,2) * t364;
t56 = t154 * t391 + t158 * t395 + t146;
t349 = t235 + t370;
t348 = t243 + t369;
t211 = rSges(5,1) * t257 - rSges(5,2) * t254;
t190 = t211 * t361;
t338 = -t365 / 0.2e1;
t336 = t363 / 0.2e1;
t334 = -t361 / 0.2e1;
t333 = t361 / 0.2e1;
t328 = -rSges(4,3) * t255 - t398 + t422;
t326 = t428 - t434;
t25 = t391 * t89 - t318;
t26 = -t391 * t90 + t426;
t317 = t25 * t258 - t255 * t26;
t316 = t25 * t255 + t258 * t26;
t27 = t389 * t89 + t425;
t28 = -t389 * t90 + t424;
t313 = t255 * t28 - t258 * t27;
t312 = t255 * t27 + t258 * t28;
t311 = t255 * t36 - t258 * t35;
t310 = t255 * t35 + t258 * t36;
t218 = rSges(5,1) * t342;
t112 = -rSges(5,2) * t344 - t160 * qJD(1) + t218;
t197 = t321 * qJD(4);
t296 = -qJ(3) * t259 - qJD(4) * t197;
t41 = -pkin(6) * t387 + t296 * t255 + (t112 + t190) * qJD(1) + t351;
t183 = t211 * t255;
t113 = qJD(4) * t183 + (t258 * t321 - t420) * qJD(1);
t42 = t241 + t296 * t258 + (-t211 * t363 - t113 - t167) * qJD(1) + t374;
t309 = t255 * t41 + t258 * t42;
t292 = -t208 + t319;
t63 = t190 + (-t160 + t292) * qJD(1) + t370;
t162 = t375 - t420;
t64 = (qJD(4) * t211 + qJD(3)) * t255 + (t162 + t458) * qJD(1) - t371;
t308 = t255 * t64 + t258 * t63;
t307 = t255 * t99 - t258 * t97;
t77 = -t153 * t254 + t157 * t257;
t78 = -t154 * t254 + t158 * t257;
t55 = t255 * t299 + t396;
t284 = (-t255 * t56 + t258 * t55) * qJD(4);
t58 = t384 - t397;
t59 = -t150 * t255 + t383;
t283 = (-t255 * t59 + t258 * t58) * qJD(4);
t282 = -t321 + t428;
t281 = qJD(4) * t206;
t280 = qJD(4) * t204;
t73 = (-t160 * t255 - t162 * t258) * qJD(4);
t276 = t148 * t228 + t192 * t90 + t286 * t89;
t275 = -(-Icges(6,5) * t170 - Icges(6,6) * t171) * t286 + (Icges(6,5) * t172 - Icges(6,6) * t173) * t192 + t174 * t228;
t274 = -t153 * t258 + t154 * t255;
t272 = t257 * t275;
t271 = (t254 * t376 + t257 * t377) * qJD(1);
t269 = (Icges(6,1) * t172 - t405 - t93) * t192 - (-Icges(6,1) * t170 - t165 - t91) * t286 + (-t152 + t180) * t228;
t32 = t192 * t97 + t286 * t99 + (-t186 * t255 - t188 * t258) * qJD(4);
t33 = (-t186 + t292) * qJD(1) + t370 + t463;
t34 = -t161 * t192 + t228 * t99 + (qJD(4) * t216 + qJD(3)) * t255 + (t188 + t458) * qJD(1) - t371;
t265 = -t32 * t307 + (t255 * t33 - t258 * t34) * t161;
t195 = t304 * qJD(4);
t196 = t306 * qJD(4);
t262 = -qJD(1) * t202 - t195 * t257 - t196 * t254 + (-t204 * t254 + t206 * t257) * qJD(4);
t261 = t274 * t254 - t257 * t451;
t260 = t450 * t257;
t240 = rSges(4,2) * t365;
t209 = rSges(3,2) * t255 + t421;
t189 = t216 * t258;
t187 = t216 * t255;
t185 = t211 * t258;
t137 = t161 * t258;
t136 = t161 * t255;
t135 = t156 * t258;
t134 = t156 * t255;
t133 = t152 * t258;
t132 = t152 * t255;
t129 = qJD(1) * t459 - t245;
t128 = t244 + (-t208 + t209) * qJD(1);
t123 = rSges(6,1) * t172 - rSges(6,2) * t173;
t122 = -rSges(6,1) * t170 - rSges(6,2) * t171;
t115 = qJD(1) * t464 + qJD(3) * t255 - t245;
t114 = (-t208 + t328) * qJD(1) + t370;
t110 = t234 + (-qJD(1) * t330 - t167) * qJD(1);
t109 = qJD(1) * t372 + t378;
t84 = t258 * t297 - t175;
t74 = t84 * qJD(1);
t72 = (-qJD(1) * t368 - t167) * qJD(1) + t293;
t71 = -qJ(3) * t390 + qJD(1) * (-rSges(4,3) * t366 + t240) + t351;
t53 = -t148 * t389 + t152 * t172 + t156 * t173;
t51 = t53 * t228;
t40 = t262 * t255 + t258 * t452;
t39 = -t255 * t452 + t262 * t258;
t38 = -qJD(4) * t298 - (-qJD(1) * t153 + t258 * t280) * t254 + (-qJD(1) * t157 + t258 * t281) * t257;
t37 = -qJD(4) * t299 - (qJD(1) * t154 + t255 * t280) * t254 + (t255 * t281 + (t258 * t306 - t402) * qJD(1)) * t257;
t24 = t74 + t283;
t23 = t284 + t386;
t16 = -t100 * t391 - t103 * t170 + t106 * t171 + t148 * t278 + t152 * t81 + t156 * t82;
t15 = -t100 * t389 + t103 * t172 + t106 * t173 + t148 * t277 + t152 * t79 + t156 * t80;
t12 = t192 * t36 + t228 * t57 - t286 * t35;
t11 = t192 * t28 - t27 * t286 + t51;
t10 = t192 * t26 - t25 * t286 - t465;
t9 = -t144 * t99 + t145 * t97 + t192 * t50 + t286 * t49 + (-t124 * t258 - t125 * t255 + (-t186 * t258 + t188 * t255) * qJD(1)) * qJD(4);
t6 = -t90 * t345 - t170 * t45 + t171 * t47 + t81 * t93 + t82 * t96 + (-t257 * t43 + t364 * t90) * t255;
t5 = t89 * t345 - t170 * t46 + t171 * t48 + t81 * t91 - t82 * t95 + (-t257 * t44 - t364 * t89) * t255;
t4 = t90 * t344 + t172 * t45 + t173 * t47 + t79 * t93 + t80 * t96 + (-t258 * t43 + t366 * t90) * t257;
t3 = -t89 * t344 + t172 * t46 + t173 * t48 + t79 * t91 - t80 * t95 + (-t258 * t44 - t366 * t89) * t257;
t2 = t144 * t25 + t145 * t26 + t16 * t228 + t192 * t6 - t273 * t340 - t286 * t5;
t1 = t144 * t27 + t145 * t28 + t15 * t228 + t192 * t4 - t286 * t3 + t340 * t53;
t17 = [(t51 - (-t26 + t470) * t286 + (t25 + t471) * t192) * t444 + (-qJD(4) * t297 + t195 * t254 - t196 * t257) * qJD(1) + (t74 + (t384 * t258 + (t469 * t255 - t383 - t55) * t255) * qJD(4)) * t334 + t427 + t430 / 0.2e1 + t412 / 0.2e1 + t413 / 0.2e1 - t431 / 0.2e1 + t15 * t445 + t53 * t447 - t273 * t448 + (t465 - (-t28 + t471) * t286 + (t27 - t470) * t192 + t10) * t446 + (t16 + t11) * t443 - (t38 + t39) * t363 / 0.2e1 + (-t386 + ((t383 - t59 - t396) * t258 + (-t146 + t56 - t58 - t397) * t255) * qJD(4) + t23) * t336 + (-(-t33 + t348 + (-t186 + t319) * qJD(1) + t463) * t34 + t14 * (t231 + t248 - t97 - t433) + t33 * (-t323 + t371) + t13 * (t458 + t408) + t34 * (t49 + t349 + t350) + (-t13 * pkin(6) + t14 * t326 + (-qJD(3) + (-t432 + (-rSges(6,3) - pkin(7)) * t254) * qJD(4)) * t33) * t255 + ((-t33 * qJ(2) + t326 * t34) * t255 + (t33 * (t215 + t418 + t428) - t34 * pkin(6)) * t258) * qJD(1)) * m(6) + (t42 * t248 + t63 * t371 + t41 * (t458 + t375) + t64 * (t218 + t349) + (t42 * t435 - t64 * t352 + (t282 * t63 + t435 * t64) * qJD(1)) * t258 + (t42 * t282 + t63 * (-rSges(5,1) * t362 - qJD(3) + t352) + t41 * t435 + (t63 * (rSges(5,3) - qJ(2)) + t64 * t282) * qJD(1)) * t255 - (t190 - t63 + (-t160 + t319) * qJD(1) + t348) * t64) * m(5) + (t72 * (t248 + t422) + t114 * t245 + t71 * t464 + t115 * (t240 + t349) + (-t114 * qJD(3) + t72 * t354) * t255 + (t114 * t354 * t258 + (t114 * (-rSges(4,2) - qJ(2)) + t115 * t354) * t255) * qJD(1) - (qJD(1) * t328 - t114 + t348) * t115) * m(4) + (t110 * (t255 * t436 + t248 + t421) + t128 * t245 + t109 * t459 + t129 * (t372 + t373) + (t128 * t436 * t258 + (t128 * (-rSges(3,3) - qJ(2)) - t129 * pkin(1)) * t255) * qJD(1) - (qJD(1) * t209 - t128 + t369) * t129) * m(3) + (t37 + t40 + t24) * t333 + ((t77 + t83) * t255 + (t78 + t84) * t258) * qJD(4) * t429; 0.2e1 * (-t417 / 0.2e1 + t416 / 0.2e1) * m(6) + 0.2e1 * (t41 * t438 + t42 * t439) * m(5) + 0.2e1 * (t438 * t71 + t439 * t72) * m(4) + 0.2e1 * (t109 * t438 + t110 * t439) * m(3); m(4) * (t255 * t71 + t258 * t72) + m(5) * t309 + m(6) * (t13 * t255 + t14 * t258); qJD(1) * (-t255 * t38 + t258 * t37 + (-t77 * t255 - t78 * t258) * qJD(1)) / 0.2e1 - t145 * t313 / 0.2e1 - t12 * t359 / 0.2e1 + ((t175 * t361 - t358) * t258 + (t271 + (-t258 * t176 + t261) * qJD(4)) * t255) * t334 + ((t176 * t363 + t358) * t255 + (t271 + (-t255 * t175 + t261) * qJD(4)) * t258) * t336 + (-qJD(1) * t312 - t255 * t4 + t258 * t3) * t445 + ((t133 * t172 + t135 * t173) * t192 - (t132 * t172 + t134 * t173) * t286 + (t151 * t172 + t155 * t173) * t228 + (t257 * t53 + t27 * t395) * qJD(5) + ((qJD(5) * t28 + t276) * t254 + t260) * t258) * t446 + t317 * t448 + ((-t254 * t377 + t257 * t376) * qJD(1) + (t254 * t451 + t274 * t257) * qJD(4)) * t429 + (-qJD(1) * t310 - t255 * t8 + t258 * t7) * t441 + (((-t133 * t253 + t135 * t256 + t90) * t192 - (-t132 * t253 + t134 * t256 - t89) * t286 + (-t151 * t253 + t155 * t256 + t148) * t228 + t57 * qJD(5)) * t257 + (qJD(5) * t310 - t450) * t254) * t442 + (-qJD(1) * t316 - t255 * t6 + t258 * t5) * t443 + ((-t133 * t170 + t135 * t171) * t192 - (-t132 * t170 + t134 * t171) * t286 + (-t151 * t170 + t155 * t171) * t228 + (-t257 * t273 + t26 * t394) * qJD(5) + ((qJD(5) * t25 + t276) * t254 + t260) * t255) * t444 - t311 * t340 / 0.2e1 - (t255 * t10 + t258 * t11) * t360 / 0.2e1 + ((t14 * t379 + t33 * t385 - t9 * t408 + t32 * (-t124 - t49) + (t32 * t409 + t34 * t379) * qJD(1)) * t258 + (t13 * t379 + t34 * t385 + t9 * t409 + t32 * (-t125 - t50) + (t32 * t408 - t33 * t379) * qJD(1)) * t255 - t33 * (-qJD(1) * t187 - t136 * t228 - t159 * t286 + t215 * t361) - t34 * (qJD(1) * t189 + t137 * t228 - t159 * t192 + t215 * t363) - t32 * (t136 * t192 + t137 * t286 - t187 * t363 - t189 * t361) - ((-t33 * t97 + t34 * t99) * t257 + t265 * t254) * qJD(5)) * m(6) + (0.2e1 * t73 * (-t112 * t258 - t113 * t255 + (-t160 * t258 + t162 * t255) * qJD(1)) - t308 * t197 + ((-t255 * t63 + t258 * t64) * qJD(1) + t309) * t211 - (-t183 * t63 + t185 * t64) * qJD(1) - (t73 * (-t183 * t255 - t185 * t258) - t308 * t321) * qJD(4)) * m(5) - (qJD(1) * t39 + t1 + (-t453 * t466 + t454 * t255 ^ 2 + (-t58 * t255 - t59 * t258) * qJD(1)) * t460) * t255 / 0.2e1 + (qJD(1) * t40 + t2 + (t453 * t258 ^ 2 - t454 * t466 + (-t55 * t255 - t56 * t258) * qJD(1)) * t460) * t258 / 0.2e1 - (t10 + t284 + t23) * t366 / 0.2e1 + (t283 + t24 + t11) * t338; -t1 * t389 / 0.2e1 + (t254 * t53 - t257 * t312) * t447 + ((qJD(4) * t312 + t15) * t254 + (qJD(1) * t313 + qJD(4) * t53 - t255 * t3 - t258 * t4) * t257) * t445 - t2 * t391 / 0.2e1 + (-t254 * t273 - t257 * t316) * t448 + ((qJD(4) * t316 + t16) * t254 + (-qJD(1) * t317 - qJD(4) * t273 - t255 * t5 - t258 * t6) * t257) * t443 + t254 * (t412 + t413 + t427 + t430 - t431) / 0.2e1 + ((qJD(4) * t310 + t18) * t254 + (qJD(1) * t311 + qJD(4) * t57 - t255 * t7 - t258 * t8) * t257) * t441 + (t172 * t449 + t269 * t173 - t258 * t272) * t446 + (-t170 * t449 + t171 * t269 - t255 * t272) * t444 + (t275 * t254 + (-t253 * t449 + t256 * t269) * t257) * t442 + (t12 + qJD(5) * (t254 * t57 - t257 * t310)) * t362 / 0.2e1 + (t346 / 0.2e1 + t254 * t333) * t11 + (t254 * t336 + t257 * t338) * t10 + ((qJD(4) * t265 + t13 * t99 - t14 * t97 - t33 * t50 + t34 * t49) * t254 + (t33 * (-qJD(4) * t97 - t111 * t255) + t34 * (qJD(4) * t99 + t111 * t258) + t9 * t307 + t32 * (t255 * t49 - t258 * t50 + t365 * t99 + t366 * t97) + (t417 - t416 + (-t255 * t34 - t258 * t33) * qJD(1)) * t161) * t257 - t33 * (-t122 * t228 - t184 * t286) - t34 * (t123 * t228 - t184 * t192) - t32 * (t122 * t192 + t123 * t286)) * m(6);];
tauc = t17(:);
