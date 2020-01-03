% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:16
% DurationCPUTime: 5.93s
% Computational Cost: add. (12381->512), mult. (10652->675), div. (0->0), fcn. (8372->8), ass. (0->313)
t263 = qJ(3) + qJ(4);
t257 = cos(t263);
t245 = Icges(5,4) * t257;
t255 = sin(t263);
t192 = Icges(5,1) * t255 + t245;
t324 = -Icges(5,2) * t255 + t245;
t470 = t192 + t324;
t266 = sin(qJ(1));
t436 = pkin(1) * qJD(1);
t370 = t266 * t436;
t264 = qJ(1) + qJ(2);
t256 = sin(t264);
t258 = cos(t264);
t195 = rSges(3,1) * t256 + rSges(3,2) * t258;
t262 = qJD(1) + qJD(2);
t417 = t195 * t262;
t153 = -t370 - t417;
t261 = qJD(3) + qJD(4);
t437 = rSges(5,2) * t257;
t371 = t261 * t437;
t440 = rSges(5,1) * t255;
t469 = -t261 * t440 - t371;
t468 = 2 * qJD(3);
t406 = t257 * t258;
t385 = rSges(5,1) * t406 + rSges(5,3) * t256;
t267 = cos(qJ(3));
t403 = t258 * t267;
t467 = rSges(4,1) * t403 + t256 * rSges(4,3);
t251 = t256 * pkin(6);
t199 = pkin(2) * t258 + t251;
t259 = Icges(4,4) * t267;
t265 = sin(qJ(3));
t325 = -Icges(4,2) * t265 + t259;
t223 = Icges(4,1) * t265 + t259;
t409 = t256 * t265;
t382 = rSges(4,2) * t409 + rSges(4,3) * t258;
t408 = t256 * t267;
t150 = rSges(4,1) * t408 - t382;
t138 = t262 * t150;
t252 = t258 * pkin(6);
t198 = pkin(2) * t256 - t252;
t180 = t262 * t198;
t405 = t258 * t262;
t228 = pkin(6) * t405;
t378 = qJD(3) * t267;
t368 = rSges(4,2) * t378;
t402 = t262 * t265;
t313 = rSges(4,2) * t256 * t402 + rSges(4,3) * t405 - t258 * t368;
t379 = qJD(3) * t265;
t359 = t258 * t379;
t466 = -rSges(4,1) * t359 + t138 + t180 + t228 + t313;
t220 = Icges(4,5) * t267 - Icges(4,6) * t265;
t219 = Icges(4,5) * t265 + Icges(4,6) * t267;
t290 = Icges(4,3) * t262 - qJD(3) * t219;
t430 = Icges(4,4) * t265;
t224 = Icges(4,1) * t267 - t430;
t312 = t224 * t258;
t149 = Icges(4,5) * t256 + t312;
t310 = t325 * t258;
t147 = Icges(4,6) * t256 + t310;
t421 = t147 * t265;
t320 = -t149 * t267 + t421;
t410 = t256 * t262;
t465 = -t220 * t410 + t258 * t290 + t262 * t320;
t429 = Icges(5,4) * t255;
t193 = Icges(5,1) * t257 - t429;
t311 = t193 * t258;
t137 = Icges(5,5) * t256 + t311;
t189 = Icges(5,5) * t257 - Icges(5,6) * t255;
t188 = Icges(5,5) * t255 + Icges(5,6) * t257;
t293 = Icges(5,3) * t262 - t188 * t261;
t309 = t324 * t258;
t135 = Icges(5,6) * t256 + t309;
t424 = t135 * t255;
t464 = -t189 * t410 + t258 * t293 + t262 * (-t137 * t257 + t424);
t307 = t189 * t258;
t411 = t256 * t257;
t413 = t255 * t256;
t134 = Icges(5,4) * t411 - Icges(5,2) * t413 - Icges(5,6) * t258;
t211 = Icges(5,4) * t413;
t136 = Icges(5,1) * t411 - Icges(5,5) * t258 - t211;
t323 = t134 * t255 - t136 * t257;
t463 = t256 * t293 + (t307 + t323) * t262;
t308 = t220 * t258;
t235 = Icges(4,4) * t409;
t148 = Icges(4,1) * t408 - Icges(4,5) * t258 - t235;
t146 = Icges(4,4) * t408 - Icges(4,2) * t409 - Icges(4,6) * t258;
t422 = t146 * t265;
t321 = -t148 * t267 + t422;
t462 = t256 * t290 + (t308 + t321) * t262;
t190 = Icges(5,2) * t257 + t429;
t319 = t190 * t255 - t192 * t257;
t461 = t189 * t261 + t262 * t319;
t229 = rSges(4,1) * t265 + rSges(4,2) * t267;
t381 = qJD(3) * t256;
t181 = t229 * t381;
t404 = t258 * t265;
t151 = -rSges(4,2) * t404 + t467;
t302 = t151 + t199;
t460 = t262 * t302 - t181;
t221 = Icges(4,2) * t267 + t430;
t318 = t221 * t265 - t223 * t267;
t459 = qJD(3) * t220 + t262 * t318;
t144 = Icges(4,5) * t408 - Icges(4,6) * t409 - Icges(4,3) * t258;
t66 = -t144 * t258 - t256 * t321;
t269 = -pkin(7) - pkin(6);
t238 = t258 * t269;
t443 = pkin(3) * t267;
t254 = pkin(2) + t443;
t386 = -t254 * t256 - t238;
t126 = t198 + t386;
t117 = t262 * t126;
t214 = rSges(5,2) * t413;
t139 = rSges(5,1) * t411 - rSges(5,3) * t258 - t214;
t125 = t262 * t139;
t367 = t257 * t410;
t389 = rSges(5,3) * t405 + t214 * t262;
t458 = -rSges(5,1) * t367 - t254 * t410 - t117 + t125 + t180 + t389;
t336 = pkin(3) * t359;
t442 = pkin(2) - t254;
t107 = -t336 - t228 + (t256 * t442 - t238) * t262;
t360 = t256 * t379;
t218 = pkin(3) * t360;
t407 = t256 * t269;
t387 = -t262 * t407 - t218;
t108 = (-t258 * t442 - t251) * t262 + t387;
t337 = t254 * t258 - t407;
t127 = t337 - t199;
t412 = t255 * t258;
t372 = rSges(5,2) * t412;
t140 = -t372 + t385;
t186 = t256 * t261;
t163 = t262 * t186;
t187 = t258 * t261;
t164 = t262 * t187;
t89 = -t258 * t371 + (-t187 * t255 - t367) * rSges(5,1) + t389;
t363 = t256 * t469 - t262 * t372;
t90 = t262 * t385 + t363;
t12 = t139 * t164 - t140 * t163 + t186 * t90 + t187 * t89 + ((t107 - t117) * t258 + (-t127 * t262 + t108) * t256) * qJD(3);
t194 = t437 + t440;
t161 = t194 * t256;
t162 = t194 * t258;
t438 = rSges(5,2) * t255;
t439 = rSges(5,1) * t257;
t196 = -t438 + t439;
t52 = t139 * t186 + t140 * t187 + (-t126 * t256 + t127 * t258) * qJD(3);
t305 = -t187 * t194 - t336;
t288 = t305 - t370;
t56 = (t126 - t139 - t198) * t262 + t288;
t397 = -t127 - t140;
t364 = t199 - t397;
t268 = cos(qJ(1));
t369 = t268 * t436;
t394 = t186 * t194 + t218;
t57 = t262 * t364 + t369 - t394;
t457 = -t56 * (t161 * t262 - t187 * t196) - t57 * (-t162 * t262 - t186 * t196) - t52 * (-t161 * t186 - t162 * t187) + t12 * (t139 * t256 + t140 * t258);
t391 = -Icges(4,2) * t408 + t148 - t235;
t393 = t223 * t256 + t146;
t456 = -t265 * t391 - t267 * t393;
t455 = t186 * (-t190 * t258 + t137) - t187 * (-Icges(5,2) * t411 + t136 - t211) + t262 * t470;
t454 = t163 / 0.2e1;
t453 = t164 / 0.2e1;
t452 = -t186 / 0.2e1;
t451 = t186 / 0.2e1;
t450 = -t187 / 0.2e1;
t449 = t187 / 0.2e1;
t448 = t256 / 0.2e1;
t447 = -t258 / 0.2e1;
t446 = -t262 / 0.2e1;
t445 = t262 / 0.2e1;
t444 = pkin(1) * t266;
t260 = t268 * pkin(1);
t441 = rSges(4,1) * t267;
t380 = qJD(3) * t258;
t361 = t229 * t380;
t304 = -t361 - t370;
t87 = (-t150 - t198) * t262 + t304;
t435 = t258 * t87;
t434 = t262 * t56;
t419 = t188 * t258;
t91 = -t256 * t319 - t419;
t433 = t91 * t262;
t415 = t219 * t258;
t110 = -t256 * t318 - t415;
t426 = t110 * t262;
t132 = Icges(5,5) * t411 - Icges(5,6) * t413 - Icges(5,3) * t258;
t425 = t132 * t258;
t420 = t188 * t256;
t418 = t190 * t261;
t416 = t219 * t256;
t414 = t220 * t262;
t401 = -t132 * t256 - t136 * t406;
t133 = Icges(5,3) * t256 + t307;
t400 = t133 * t256 + t137 * t406;
t399 = -t144 * t256 - t148 * t403;
t145 = Icges(4,3) * t256 + t308;
t398 = t145 * t256 + t149 * t403;
t392 = -t223 * t258 - t147;
t390 = -t221 * t258 + t149;
t384 = -t221 + t224;
t383 = t223 + t325;
t271 = qJD(1) ^ 2;
t377 = t271 * t444;
t376 = t271 * t260;
t375 = (qJD(3) ^ 2) * t443;
t374 = t256 * t90 + (t125 + t89) * t258;
t365 = t258 * t402;
t362 = rSges(4,1) * t360 + rSges(4,2) * t365 + t256 * t368;
t358 = t410 / 0.2e1;
t357 = t405 / 0.2e1;
t356 = -pkin(2) - t441;
t355 = -t381 / 0.2e1;
t352 = t380 / 0.2e1;
t350 = -pkin(3) * t265 - t194;
t197 = rSges(3,1) * t258 - rSges(3,2) * t256;
t295 = Icges(5,5) * t262 - t192 * t261;
t348 = -t134 * t261 + t256 * t295 + t262 * t311;
t347 = -t135 * t261 - t193 * t410 + t258 * t295;
t294 = Icges(5,6) * t262 - t418;
t346 = t136 * t261 + t256 * t294 + t262 * t309;
t345 = t137 * t261 + t258 * t294 - t324 * t410;
t112 = t137 * t411;
t344 = t133 * t258 - t112;
t122 = t149 * t408;
t343 = t145 * t258 - t122;
t342 = -t132 + t424;
t340 = -t144 + t421;
t339 = t470 * t261;
t338 = t193 * t261 - t418;
t333 = t262 * (-pkin(2) * t410 + t228) - t377;
t170 = rSges(3,1) * t405 - rSges(3,2) * t410;
t168 = t196 * t261;
t332 = -pkin(3) * t378 - t168;
t330 = -rSges(4,2) * t265 + t441;
t329 = -t256 * t57 - t258 * t56;
t88 = t369 + t460;
t328 = -t256 * t88 - t435;
t77 = t134 * t257 + t136 * t255;
t97 = t146 * t267 + t148 * t265;
t98 = t147 * t267 + t149 * t265;
t317 = t337 + t385;
t67 = -t147 * t409 - t343;
t315 = (t256 * t67 - t258 * t66) * qJD(3);
t68 = -t146 * t404 - t399;
t69 = -t147 * t404 + t398;
t314 = (t256 * t69 - t258 * t68) * qJD(3);
t306 = t323 * t256;
t94 = (t150 * t256 + t151 * t258) * qJD(3);
t303 = -t139 + t386;
t301 = t186 * t419 - t187 * t420 - t189 * t262;
t300 = -t265 * t390 + t267 * t392;
t299 = t256 * t356 + t252 + t382;
t298 = -rSges(5,3) * t410 - t363 - t387;
t281 = t132 * t262 - t255 * t346 + t257 * t348;
t13 = t256 * t463 + t258 * t281;
t280 = t133 * t262 - t255 * t345 + t257 * t347;
t14 = t256 * t464 + t258 * t280;
t15 = t256 * t281 - t258 * t463;
t16 = t256 * t280 - t258 * t464;
t284 = (-t192 * t258 - t135) * t186 - (-t192 * t256 - t134) * t187 + (-t190 + t193) * t262;
t272 = -t255 * t455 + t257 * t284;
t61 = -t306 - t425;
t62 = -t135 * t413 - t344;
t28 = t186 * t62 - t187 * t61 + t433;
t63 = -t134 * t412 - t401;
t64 = -t135 * t412 + t400;
t92 = -t258 * t319 + t420;
t76 = t92 * t262;
t29 = t186 * t64 - t187 * t63 + t76;
t42 = t255 * t348 + t257 * t346;
t43 = t255 * t347 + t257 * t345;
t279 = t188 * t262 - t255 * t339 + t257 * t338;
t44 = t256 * t461 + t258 * t279;
t45 = t256 * t279 - t258 * t461;
t78 = t135 * t257 + t137 * t255;
t297 = (-t13 * t187 + t14 * t186 + t163 * t63 + t164 * t64 + t262 * t44) * t448 + (-t256 * t301 + t258 * t272) * t452 + (t256 * t272 + t258 * t301) * t449 + (-t15 * t187 + t16 * t186 + t163 * t61 + t164 * t62 + t262 * t45) * t447 + (t255 * t284 + t257 * t455) * t446 + t28 * t358 + t29 * t357 + ((t262 * t64 - t13) * t258 + (t262 * t63 + t14) * t256) * t451 + (t256 * t62 - t258 * t61) * t454 + (t256 * t64 - t258 * t63) * t453 + ((t262 * t62 - t15) * t258 + (t262 * t61 + t16) * t256) * t450 + ((t262 * t78 - t42) * t258 + (t262 * t77 + t43) * t256) * t445;
t296 = (-t265 * t383 + t267 * t384) * t262;
t292 = Icges(4,5) * t262 - qJD(3) * t223;
t291 = Icges(4,6) * t262 - qJD(3) * t221;
t101 = t258 * t291 - t325 * t410;
t103 = -t224 * t410 + t258 * t292;
t278 = -qJD(3) * t98 - t101 * t265 + t103 * t267 + t145 * t262;
t102 = t256 * t291 + t262 * t310;
t104 = t256 * t292 + t262 * t312;
t277 = -qJD(3) * t97 - t102 * t265 + t104 * t267 + t144 * t262;
t201 = t325 * qJD(3);
t202 = t224 * qJD(3);
t276 = -t201 * t265 + t202 * t267 + t219 * t262 + (-t221 * t267 - t223 * t265) * qJD(3);
t275 = (t356 * t435 + (t87 * (-rSges(4,3) - pkin(6)) + t88 * t356) * t256) * t262;
t111 = -t258 * t318 + t416;
t109 = t111 * t262;
t32 = t315 + t426;
t33 = t109 + t314;
t48 = -qJD(3) * t321 + t102 * t267 + t104 * t265;
t49 = -qJD(3) * t320 + t101 * t267 + t103 * t265;
t53 = t256 * t459 + t258 * t276;
t54 = t256 * t276 - t258 * t459;
t274 = (t76 + (t62 + (t134 * t258 + t135 * t256) * t255 + t344 + t401) * t187 + (-t136 * t411 + t425 + t61 + (t134 * t256 - t135 * t258) * t255 + t400) * t186) * t449 + (t109 + ((t67 - t122 + (t145 + t422) * t258 + t399) * t258 + t398 * t256) * qJD(3)) * t352 + (t77 + t91) * t454 + (t78 + t92) * t453 + (t28 - t433 + (t64 - t306 - t400) * t187 + (t342 * t256 - t112 + t63) * t186 + ((t133 + t323) * t186 + t342 * t187) * t258) * t452 + (t43 + t44) * t451 + (-t426 + ((t258 * t340 - t398 + t69) * t258 + (t256 * t340 + t343 + t68) * t256) * qJD(3) + t32) * t355 + (t49 + t53) * t381 / 0.2e1 + (-qJD(3) * t318 + t201 * t267 + t202 * t265 + t255 * t338 + t257 * t339) * t262 + (t42 + t45 + t29) * t450 - (t48 + t54 + t33) * t380 / 0.2e1 + ((t110 + t97) * t256 + (t111 + t98) * t258) * qJD(3) * t445;
t36 = -t256 * t375 - t164 * t194 - t168 * t186 + (t107 + t89 - t336) * t262 + t333;
t273 = (-t36 * t438 + t57 * (-pkin(3) * t379 + t469) + (t56 * (-t254 - t439) - t57 * t269) * t262) * t258;
t205 = t330 * qJD(3);
t179 = t229 * t258;
t178 = t229 * t256;
t171 = t199 * t262;
t154 = t197 * t262 + t369;
t129 = -t170 * t262 - t376;
t128 = -t262 * t417 - t377;
t106 = t262 * t467 - t362;
t105 = (-t262 * t408 - t359) * rSges(4,1) + t313;
t59 = -t376 - t205 * t380 + (-t106 - t171 + t181) * t262;
t58 = t105 * t262 + (-t205 * t256 - t229 * t405) * qJD(3) + t333;
t37 = -t258 * t375 - t376 + t163 * t194 - t168 * t187 + (-t108 - t171 - t90 + t218) * t262;
t1 = [m(3) * (t129 * (-t195 - t444) + t128 * (t197 + t260) + (-t170 - t369 + t154) * t153) + t274 + (t37 * (t303 - t444) + t56 * (t298 - t369) + t36 * (t260 + t317) + t273 + (-t288 + t56 - t370 + t458) * t57) * m(5) + (t59 * (t299 - t444) + t87 * (t362 - t369) + t58 * (t260 + t302) + t275 + (-t370 - t304 + t87 + t466) * t88) * m(4); t274 + (t37 * t303 + t36 * t317 + t364 * t434 + t273 + (-t305 + t458) * t57 + (-t394 + t298) * t56) * m(5) + (t59 * t299 + t58 * t302 + t275 + (t361 + t466) * t88 + (t362 + t460) * t87) * m(4) + (-(-t153 * t197 - t154 * t195) * t262 + t128 * t197 - t129 * t195 - t153 * t170 - t154 * t417) * m(3); ((t265 * t384 + t267 * t383) * t262 + ((t256 * t390 - t258 * t391) * t267 + (t256 * t392 + t258 * t393) * t265) * qJD(3)) * t446 + ((t262 * t98 - t48) * t258 + (t262 * t97 + t49) * t256) * t445 + ((-t380 * t416 - t414) * t258 + (t296 + (t300 * t256 + (t415 - t456) * t258) * qJD(3)) * t256) * t352 + ((-t381 * t415 + t414) * t256 + (t296 + (-t456 * t258 + (t416 + t300) * t256) * qJD(3)) * t258) * t355 + t297 + (t262 * t53 + ((-t256 * t462 - t258 * t277 + t262 * t69) * t258 + (t256 * t465 + t258 * t278 + t262 * t68) * t256) * t468) * t448 + (t262 * t54 + ((-t256 * t277 + t258 * t462 + t262 * t67) * t258 + (t256 * t278 - t258 * t465 + t262 * t66) * t256) * t468) * t447 + (t315 + t32) * t358 + (t314 + t33) * t357 + (-(-t57 * t365 + (t329 * t267 + t52 * (-t256 ^ 2 - t258 ^ 2) * t265) * qJD(3)) * pkin(3) + t52 * t374 + (t37 * t350 + t56 * t332 + t12 * t127 + t52 * t107 + (-t126 * t52 + t350 * t57) * t262) * t258 + (t36 * t350 + t57 * t332 - t12 * t126 + t52 * t108 + (t194 * t56 + t397 * t52) * t262) * t256 + t457) * m(5) + (0.2e1 * t94 * ((t105 + t138) * t258 + (-t151 * t262 + t106) * t256) + t328 * t205 + ((-t262 * t88 - t59) * t258 + (t262 * t87 - t58) * t256) * t229 - (t178 * t87 - t179 * t88) * t262 - (t94 * (-t178 * t256 - t179 * t258) + t328 * t330) * qJD(3)) * m(4); t297 + (t52 * (-t140 * t410 + t374) + t329 * t168 + ((-t262 * t57 - t37) * t258 + (-t36 + t434) * t256) * t194 + t457) * m(5);];
tauc = t1(:);
