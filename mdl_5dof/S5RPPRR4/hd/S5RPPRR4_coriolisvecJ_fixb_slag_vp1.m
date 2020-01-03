% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR4
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:27
% EndTime: 2020-01-03 11:31:09
% DurationCPUTime: 18.15s
% Computational Cost: add. (18911->726), mult. (21699->989), div. (0->0), fcn. (20721->10), ass. (0->344)
t326 = cos(qJ(1));
t325 = sin(qJ(1));
t469 = qJ(2) * t325;
t272 = pkin(1) * t326 + t469;
t420 = qJD(1) * t326;
t421 = qJD(1) * t325;
t426 = pkin(1) * t420 + qJ(2) * t421;
t535 = -qJD(1) * t272 + t426;
t323 = cos(pkin(8));
t299 = -qJD(4) * t323 + qJD(1);
t318 = pkin(9) + qJ(4);
t309 = sin(t318);
t310 = cos(t318);
t321 = sin(pkin(8));
t201 = -Icges(5,3) * t323 + (Icges(5,5) * t310 - Icges(5,6) * t309) * t321;
t473 = Icges(5,4) * t310;
t202 = -Icges(5,6) * t323 + (-Icges(5,2) * t309 + t473) * t321;
t474 = Icges(5,4) * t309;
t203 = -Icges(5,5) * t323 + (Icges(5,1) * t310 - t474) * t321;
t454 = t326 * t309;
t456 = t325 * t310;
t237 = t323 * t454 - t456;
t458 = t323 * t326;
t238 = t309 * t325 + t310 * t458;
t460 = t321 * t326;
t336 = t201 * t460 - t202 * t237 + t203 * t238;
t534 = t336 * t299;
t468 = qJ(3) * t321;
t491 = pkin(2) * t323;
t367 = t468 + t491;
t258 = t367 * t326;
t418 = qJD(3) * t321;
t289 = t325 * t418;
t533 = -qJD(1) * t258 + t289 + t535;
t322 = cos(pkin(9));
t303 = pkin(3) * t322 + pkin(2);
t263 = pkin(4) * t310 + t303;
t489 = pkin(4) * t309;
t320 = sin(pkin(9));
t490 = pkin(3) * t320;
t267 = t489 + t490;
t403 = t323 * t420;
t484 = pkin(4) * qJD(4);
t507 = t484 * t309 * t323;
t532 = t263 * t403 + t267 * t421 - t325 * t507;
t146 = Icges(5,5) * t238 - Icges(5,6) * t237 + Icges(5,3) * t460;
t224 = Icges(5,4) * t238;
t150 = Icges(5,2) * t237 - Icges(5,6) * t460 - t224;
t223 = Icges(5,4) * t237;
t152 = Icges(5,1) * t238 + Icges(5,5) * t460 - t223;
t64 = t323 * t146 - (t150 * t309 + t152 * t310) * t321;
t311 = qJ(5) + t318;
t301 = sin(t311);
t455 = t326 * t301;
t302 = cos(t311);
t457 = t325 * t302;
t219 = t323 * t455 - t457;
t220 = t301 * t325 + t302 * t458;
t127 = Icges(6,5) * t220 - Icges(6,6) * t219 + Icges(6,3) * t460;
t207 = Icges(6,4) * t220;
t131 = Icges(6,2) * t219 - Icges(6,6) * t460 - t207;
t206 = Icges(6,4) * t219;
t133 = Icges(6,1) * t220 + Icges(6,5) * t460 - t206;
t56 = t323 * t127 - (t131 * t301 + t133 * t302) * t321;
t361 = t150 * t237 + t152 * t238;
t459 = t323 * t325;
t235 = -t309 * t459 - t310 * t326;
t236 = t323 * t456 - t454;
t461 = t321 * t325;
t145 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t461;
t475 = Icges(5,4) * t236;
t148 = Icges(5,2) * t235 + Icges(5,6) * t461 + t475;
t222 = Icges(5,4) * t235;
t151 = Icges(5,1) * t236 + Icges(5,5) * t461 + t222;
t51 = t145 * t461 + t148 * t235 + t151 * t236;
t531 = -t361 + t51;
t412 = t310 * t484;
t324 = -pkin(6) - qJ(3);
t278 = t324 * t460;
t463 = t320 * t325;
t375 = -pkin(3) * t463 + t278;
t502 = t323 * (pkin(2) - t303) + t468;
t184 = t326 * t502 + t375;
t529 = qJD(1) * t184 + t533;
t436 = -t263 * t458 - t267 * t325;
t316 = -pkin(7) + t324;
t464 = t316 * t321;
t465 = t303 * t323;
t109 = (t464 + t465) * t326 - t375 + t436;
t136 = rSges(6,1) * t220 - rSges(6,2) * t219 + rSges(6,3) * t460;
t199 = -rSges(6,3) * t323 + (rSges(6,1) * t302 - rSges(6,2) * t301) * t321;
t319 = qJD(4) + qJD(5);
t386 = t319 * t321;
t256 = t326 * t386;
t265 = -t319 * t323 + qJD(1);
t424 = t316 - t324;
t432 = t263 - t303;
t186 = t321 * t432 + t323 * t424;
t416 = qJD(4) * t321;
t346 = (-t186 * t416 - qJD(2)) * t326;
t528 = -t109 * t299 + t265 * t136 - t199 * t256 + t289 + t346;
t255 = t325 * t386;
t217 = -t301 * t459 - t302 * t326;
t218 = t323 * t457 - t455;
t126 = Icges(6,5) * t218 + Icges(6,6) * t217 + Icges(6,3) * t461;
t472 = Icges(6,4) * t218;
t129 = Icges(6,2) * t217 + Icges(6,6) * t461 + t472;
t205 = Icges(6,4) * t217;
t132 = Icges(6,1) * t218 + Icges(6,5) * t461 + t205;
t46 = t126 * t461 + t129 * t217 + t132 * t218;
t527 = -t131 * t217 + t218 * t133;
t47 = -t127 * t461 - t527;
t196 = -Icges(6,3) * t323 + (Icges(6,5) * t302 - Icges(6,6) * t301) * t321;
t470 = Icges(6,4) * t302;
t197 = -Icges(6,6) * t323 + (-Icges(6,2) * t301 + t470) * t321;
t471 = Icges(6,4) * t301;
t198 = -Icges(6,5) * t323 + (Icges(6,1) * t302 - t471) * t321;
t68 = t196 * t461 + t197 * t217 + t198 * t218;
t17 = t255 * t46 - t256 * t47 + t265 * t68;
t525 = -t150 * t235 + t152 * t236;
t49 = t127 * t460 + t131 * t219 + t133 * t220;
t204 = -rSges(5,3) * t323 + (rSges(5,1) * t310 - rSges(5,2) * t309) * t321;
t521 = t204 * t416;
t404 = t321 * t420;
t419 = qJD(2) * t326;
t516 = -rSges(4,2) * (t320 * t458 - t322 * t325) + rSges(4,1) * (t322 * t458 + t463);
t517 = rSges(4,3) * t404 + qJD(1) * t516 - t419;
t462 = t320 * t326;
t515 = -(t322 * t459 - t462) * rSges(4,1) - (-t320 * t459 - t322 * t326) * rSges(4,2);
t157 = rSges(5,1) * t238 - rSges(5,2) * t237 + rSges(5,3) * t460;
t513 = t157 * t299 + t289;
t337 = t196 * t460 - t197 * t219 + t198 * t220;
t512 = t256 * t49 + t337 * t265;
t422 = qJD(1) * t321;
t511 = t326 * (-rSges(3,2) * t422 - qJD(2)) + rSges(3,1) * t403 + rSges(3,3) * t421;
t506 = -t129 * t219 + t132 * t220;
t452 = -t148 * t237 + t151 * t238;
t295 = rSges(3,1) * t459;
t504 = rSges(3,2) * t461 + rSges(3,3) * t326 - t295;
t338 = t325 * (-Icges(5,2) * t236 + t151 + t222) - t326 * (Icges(5,2) * t238 - t152 + t223);
t501 = t325 * (-Icges(5,1) * t235 + t148 + t475) - t326 * (-Icges(5,1) * t237 + t150 - t224);
t215 = (-Icges(6,2) * t302 - t471) * t321;
t332 = t255 * (-Icges(6,2) * t218 + t132 + t205) - t256 * (Icges(6,2) * t220 - t133 + t206) + t265 * (t198 + t215);
t216 = (-Icges(6,1) * t301 - t470) * t321;
t500 = t255 * (-Icges(6,1) * t217 + t129 + t472) - t256 * (-Icges(6,1) * t219 + t131 - t207) + t265 * (t197 - t216);
t362 = qJD(1) * t386;
t244 = t325 * t362;
t499 = t244 / 0.2e1;
t245 = t326 * t362;
t498 = t245 / 0.2e1;
t497 = -t255 / 0.2e1;
t496 = t255 / 0.2e1;
t495 = -t256 / 0.2e1;
t494 = t256 / 0.2e1;
t492 = -t323 / 0.2e1;
t486 = rSges(3,2) * t321;
t142 = -qJD(1) * t219 - t218 * t319;
t143 = qJD(1) * t220 + t217 * t319;
t78 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t404;
t80 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t404;
t82 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t404;
t27 = -t323 * t78 + ((-t129 * t319 + t82) * t302 + (-t132 * t319 - t80) * t301) * t321;
t481 = t27 * t255;
t140 = qJD(1) * t217 + t220 * t319;
t141 = qJD(1) * t218 + t219 * t319;
t405 = t321 * t421;
t77 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t405;
t79 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t405;
t81 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t405;
t28 = -t323 * t77 + ((-t131 * t319 + t81) * t302 + (t133 * t319 - t79) * t301) * t321;
t480 = t28 * t256;
t55 = -t126 * t323 + (-t129 * t301 + t132 * t302) * t321;
t479 = t55 * t245;
t478 = t56 * t244;
t477 = rSges(3,3) + qJ(2);
t368 = rSges(6,1) * t141 + rSges(6,2) * t140;
t83 = rSges(6,3) * t405 + t368;
t476 = t199 * t405 + t323 * t83;
t467 = t145 * t326;
t466 = t146 * t325;
t246 = t263 * t459;
t264 = t303 * t459;
t392 = t424 * t321;
t108 = t246 - t264 + (-t267 + t490) * t326 - t325 * t392;
t135 = rSges(6,1) * t218 + rSges(6,2) * t217 + rSges(6,3) * t461;
t453 = -t108 - t135;
t443 = -t186 - t199;
t312 = qJD(2) * t325;
t425 = qJ(2) * t420 + t312;
t251 = pkin(1) * t421 - t425;
t290 = t326 * t418;
t440 = -t367 * t421 - t251 + t290;
t233 = (-Icges(5,1) * t309 - t473) * t321;
t439 = t202 - t233;
t232 = (-Icges(5,2) * t310 - t474) * t321;
t438 = t203 + t232;
t257 = pkin(2) * t459 + qJ(3) * t461;
t315 = t325 * pkin(1);
t270 = -qJ(2) * t326 + t315;
t435 = t257 + t270;
t434 = t258 + t272;
t414 = qJD(1) * t490;
t433 = t303 * t403 + t325 * t414;
t431 = t324 * t405 + t326 * t414;
t304 = qJD(1) * t312;
t399 = qJD(1) * t418;
t430 = t326 * t399 + t304;
t429 = pkin(2) * t403 + qJ(3) * t404;
t427 = -t290 - t312;
t371 = rSges(3,1) * t323 - t486;
t240 = rSges(3,3) * t325 + t326 * t371;
t189 = -t419 + (t240 + t272) * qJD(1);
t423 = qJD(1) * t189;
t417 = qJD(3) * t323;
t84 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t404;
t415 = -t136 * t404 + t460 * t84 + t461 * t83;
t52 = -t146 * t461 - t525;
t411 = t421 * t502 + t431 + t440;
t173 = -qJD(1) * t237 - qJD(4) * t236;
t174 = qJD(1) * t238 + qJD(4) * t235;
t96 = rSges(5,1) * t174 + rSges(5,2) * t173 + rSges(5,3) * t404;
t243 = qJD(1) * (-t419 + t426);
t410 = qJD(1) * (t289 + t429) + t243 + t325 * t399;
t156 = rSges(5,1) * t236 + rSges(5,2) * t235 + rSges(5,3) * t461;
t408 = rSges(4,3) * t461 - t515;
t406 = t290 + t425;
t402 = t325 * t416;
t401 = t461 / 0.2e1;
t400 = -t460 / 0.2e1;
t398 = t422 / 0.2e1;
t397 = -t416 / 0.2e1;
t396 = t416 / 0.2e1;
t395 = pkin(1) + t465;
t394 = qJ(2) + t490;
t340 = t237 * qJD(4);
t87 = pkin(4) * t340 + (-t267 * t326 + (t323 * t432 - t464) * t325) * qJD(1) + t431;
t88 = (-qJD(1) * t392 - t412) * t326 - t433 + t532;
t10 = -t135 * t244 - t136 * t245 + t255 * t83 + t256 * t84 + (t325 * t87 + t326 * t88 + (-t108 * t325 + t109 * t326) * qJD(1)) * t416;
t393 = t10 * (t135 * t460 - t136 * t461);
t166 = rSges(6,1) * t217 - rSges(6,2) * t218;
t167 = rSges(6,1) * t219 + rSges(6,2) * t220;
t391 = t166 * t256 + t167 * t255;
t221 = (-rSges(6,1) * t301 - rSges(6,2) * t302) * t321;
t390 = t166 * t265 - t221 * t255;
t389 = -t167 * t265 - t221 * t256;
t388 = t321 ^ 2 * qJD(4) ^ 2 * t489;
t387 = t289 - t419;
t384 = qJD(1) * (-t324 * t404 - t429 + t433) + t410;
t383 = t186 * t402;
t382 = t204 * t402;
t381 = t325 * t398;
t380 = t326 * t398;
t379 = t325 * t397;
t378 = t325 * t396;
t377 = t326 * t397;
t376 = t326 * t396;
t374 = qJD(1) * t396;
t373 = (-t184 + t434) * qJD(1);
t370 = t515 * qJD(1);
t171 = qJD(1) * t235 + qJD(4) * t238;
t172 = qJD(1) * t236 + t340;
t369 = rSges(5,1) * t172 + rSges(5,2) * t171;
t347 = (-rSges(5,1) * t309 - rSges(5,2) * t310) * t321;
t213 = qJD(4) * t347;
t95 = rSges(5,3) * t405 + t369;
t44 = -t213 * t326 * t416 - t299 * t95 + (t382 + t411) * qJD(1) + t430;
t345 = (-qJD(2) - t521) * t326;
t45 = qJD(1) * t345 - t213 * t402 + t299 * t96 + t384;
t366 = -t325 * t45 - t326 * t44;
t365 = t52 * t325 + t51 * t326;
t342 = (-pkin(3) * t462 - t324 * t461 - t257 + t264 + t435) * qJD(1) + t427;
t61 = t156 * t299 + t342 - t382;
t62 = t345 + t373 + t513;
t364 = -t325 * t61 - t326 * t62;
t363 = t325 * t95 + t326 * t96;
t360 = t156 * t326 - t157 * t325;
t359 = (Icges(5,5) * t235 - Icges(5,6) * t236) * t325 - (Icges(5,5) * t237 + Icges(5,6) * t238) * t326;
t358 = t326 * t374;
t357 = t325 * t374;
t356 = pkin(1) + t371;
t48 = -t126 * t460 - t506;
t354 = rSges(4,3) * t460 + t516;
t344 = (t325 * t51 - t326 * t52) * t321;
t53 = -t145 * t460 - t452;
t54 = t146 * t460 + t361;
t343 = (t325 * t53 - t326 * t54) * t321;
t231 = (-Icges(5,5) * t309 - Icges(5,6) * t310) * t321;
t214 = (-Icges(6,5) * t301 - Icges(6,6) * t302) * t321;
t341 = (Icges(6,5) * t217 - Icges(6,6) * t218) * t255 - (Icges(6,5) * t219 + Icges(6,6) * t220) * t256 + t214 * t265;
t11 = t129 * t140 + t132 * t141 + t219 * t80 - t220 * t82 + (t126 * t421 - t326 * t78) * t321;
t12 = t131 * t140 - t133 * t141 + t219 * t79 - t220 * t81 + (-t127 * t421 - t326 * t77) * t321;
t13 = t129 * t142 + t132 * t143 + t217 * t80 + t218 * t82 + (t126 * t420 + t325 * t78) * t321;
t14 = t131 * t142 - t133 * t143 + t217 * t79 + t218 * t81 + (-t127 * t420 + t325 * t77) * t321;
t18 = t255 * t48 - t512;
t191 = t319 * t214;
t192 = t319 * t215;
t193 = t319 * t216;
t34 = t140 * t197 + t141 * t198 + t192 * t219 - t193 * t220 + (-t191 * t326 + t196 * t421) * t321;
t35 = t142 * t197 + t143 * t198 + t192 * t217 + t193 * t218 + (t191 * t325 + t196 * t420) * t321;
t57 = -t191 * t323 + ((-t197 * t319 + t193) * t302 + (-t198 * t319 - t192) * t301) * t321;
t50 = t57 * t265;
t335 = (t13 * t255 - t14 * t256 + t244 * t47 + t245 * t46 + t265 * t35) * t401 + t18 * t381 + t17 * t380 + (t323 * t337 + (t325 * t48 - t326 * t49) * t321) * t499 + (t11 * t255 - t12 * t256 + t244 * t49 + t245 * t48 + t265 * t34) * t400 + (-t323 * t68 + (t325 * t46 - t326 * t47) * t321) * t498 + (-t323 * t35 + (t13 * t325 - t14 * t326 + (t325 * t47 + t326 * t46) * qJD(1)) * t321) * t496 + (-t323 * t34 + (t11 * t325 - t12 * t326 + (t325 * t49 + t326 * t48) * qJD(1)) * t321) * t495 + (t478 + t479 - t480 + t50 + t481) * t492 + t265 * (-t323 * t57 + (t27 * t325 - t28 * t326 + (t325 * t56 + t326 * t55) * qJD(1)) * t321) / 0.2e1 + (t217 * t332 - t218 * t500 + t341 * t461) * t497 + (t219 * t332 + t220 * t500 - t341 * t460) * t494 - (-t341 * t323 + (-t301 * t332 - t302 * t500) * t321) * t265 / 0.2e1;
t89 = Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t405;
t90 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t404;
t91 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t405;
t92 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t404;
t93 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t405;
t94 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t404;
t331 = ((t148 * t171 + t151 * t172 + t237 * t92 - t238 * t94 + (t145 * t421 - t326 * t90) * t321) * t325 - (t150 * t171 - t152 * t172 + t237 * t91 - t238 * t93 + (-t146 * t421 - t326 * t89) * t321) * t326 + (t54 * t325 + t53 * t326) * qJD(1)) * t321;
t330 = (qJD(1) * t365 + (t148 * t173 + t151 * t174 + t235 * t92 + t236 * t94 + (t145 * t420 + t325 * t90) * t321) * t325 - (t150 * t173 - t152 * t174 + t235 * t91 + t236 * t93 + (-t146 * t420 + t325 * t89) * t321) * t326) * t321;
t31 = -t323 * t90 + (-t309 * t92 + t310 * t94 + (-t148 * t310 - t151 * t309) * qJD(4)) * t321;
t32 = -t323 * t89 + (-t309 * t91 + t310 * t93 + (-t150 * t310 + t152 * t309) * qJD(4)) * t321;
t63 = -t145 * t323 + (-t148 * t309 + t151 * t310) * t321;
t329 = (t31 * t325 - t32 * t326 + (t325 * t64 + t326 * t63) * qJD(1)) * t321;
t226 = t237 * pkin(4);
t225 = t235 * pkin(4);
t212 = qJD(4) * t233;
t211 = qJD(4) * t232;
t210 = qJD(4) * t231;
t194 = t319 * t221;
t182 = rSges(5,1) * t237 + rSges(5,2) * t238;
t181 = rSges(5,1) * t235 - rSges(5,2) * t236;
t155 = qJD(1) * t511 + t243;
t154 = t304 + (qJD(1) * t504 - t251) * qJD(1);
t123 = t323 * t136;
t101 = (t354 + t434) * qJD(1) + t387;
t86 = qJD(1) * t517 + t410;
t85 = (-rSges(4,3) * t405 + t370 + t440) * qJD(1) + t430;
t75 = t201 * t461 + t202 * t235 + t203 * t236;
t74 = t360 * t416 - t417;
t70 = t75 * t299;
t66 = -t210 * t323 + (-t211 * t309 + t212 * t310 + (-t202 * t310 - t203 * t309) * qJD(4)) * t321;
t65 = t66 * t299;
t43 = t173 * t202 + t174 * t203 + t211 * t235 + t212 * t236 + (t201 * t420 + t210 * t325) * t321;
t42 = t171 * t202 + t172 * t203 + t211 * t237 - t212 * t238 + (t201 * t421 - t210 * t326) * t321;
t41 = -t417 + t135 * t256 - t136 * t255 + (t108 * t326 + t109 * t325) * t416;
t40 = t373 + t528;
t39 = t108 * t299 + t135 * t265 - t199 * t255 + t342 - t383;
t33 = ((-t156 * t325 - t157 * t326) * qJD(1) + t363) * t416;
t30 = qJD(1) * t346 - t194 * t255 - t199 * t245 + t265 * t84 + t299 * t88 + t325 * t388 + t384;
t29 = t326 * t388 - t194 * t256 + t199 * t244 - t265 * t83 - t299 * t87 + (t383 + t411) * qJD(1) + t430;
t26 = qJD(4) * t343 - t534;
t25 = qJD(4) * t344 + t70;
t1 = [t478 / 0.2e1 + t479 / 0.2e1 + t481 / 0.2e1 + t65 + t35 * t496 + t68 * t498 - t337 * t499 + (t70 + ((t54 + t531) * t325 + (t53 + (-t466 + t467) * t321 - t52 + t452) * t326) * t416) * t376 + t50 - t480 / 0.2e1 + t17 * t494 + ((t47 + (t126 * t326 + t127 * t325) * t321 + t506 + t527) * t255 + t18 + t512) * t497 + (t34 + t17) * t495 + (t31 + t43) * t378 + (t63 + t75) * t358 + (t64 - t336) * t357 + (t29 * (t136 - t436) + t30 * (t246 + t315 + t135) + (t29 * qJ(2) - t30 * t464) * t325 + (t29 * (pkin(1) - t464) + t30 * (-qJ(2) - t267)) * t326 + (-t368 + t406 + (t412 + (-t263 * t323 - pkin(1) + (-rSges(6,3) + t316) * t321) * qJD(1)) * t325 + (qJD(1) * t267 - t507) * t326) * t40 + (t40 + t84 + (-t316 * t422 - qJD(2) - t412) * t326 - t528 + t529 + t532) * t39) * m(6) + (t44 * (t157 - t278) + t62 * (-t369 + t406 + t431) + t45 * (t264 + t315 + t156) + (-t45 * t394 + t44 * t395) * t326 + (t44 * t394 - t45 * t321 * t324 + t62 * (-rSges(5,3) * t321 - t395) * qJD(1)) * t325 + (t433 - t513 + t62 + t96 + (-t324 * t422 + t521) * t326 + t529) * t61) * m(5) + (t85 * (t354 + t469) + t86 * (t315 + t408 + t257) + (t85 * (pkin(1) + t367) - t86 * qJ(2)) * t326 + (t370 + t406 + (-t491 - pkin(1) + (-rSges(4,3) - qJ(3)) * t321) * t421) * t101 + (-qJD(1) * t354 + t101 - t387 + t429 + t517 + t533) * ((t408 + t435) * qJD(1) + t427)) * m(4) + (t189 * t425 + t155 * (t295 + t315) + (t154 * t477 - t155 * t486 - t356 * t423) * t325 + (rSges(3,3) * t423 + t154 * t356 - t155 * t477) * t326 + (-qJD(1) * t240 + t189 + t419 + t511 + t535) * (-t312 + (t270 - t504) * qJD(1))) * m(3) + (((t452 + t525) * t325 - t531 * t326 + ((t466 + t467) * t325 + t146 * t326 ^ 2) * t321 + t365) * t416 + t26 + t534) * t379 + (t32 + t42 + t25) * t377; m(3) * (-t154 * t326 - t155 * t325) + m(4) * (-t325 * t86 - t326 * t85) + m(5) * t366 + m(6) * (-t29 * t326 - t30 * t325); 0.2e1 * (-m(5) * t33 / 0.2e1 - m(6) * t10 / 0.2e1) * t323 + 0.2e1 * (m(4) * (t325 * t85 - t326 * t86) / 0.2e1 + m(5) * (t325 * t44 - t326 * t45) / 0.2e1 + m(6) * (t29 * t325 - t30 * t326) / 0.2e1) * t321; (qJD(4) * t329 + t65) * t492 + t299 * (-t323 * t66 + t329) / 0.2e1 + t25 * t380 + t335 + (-t323 * t43 + t330) * t378 + (-t323 * t42 + t331) * t377 + t26 * t381 - t299 * (-t323 * t231 * t299 + ((-t309 * t438 - t310 * t439) * t299 + ((-t309 * t338 - t310 * t501) * t321 - t359 * t323) * qJD(4)) * t321) / 0.2e1 + (qJD(4) * t330 + t299 * t43) * t401 + (qJD(4) * t331 + t299 * t42) * t400 + ((t231 * t461 + t235 * t438 - t236 * t439) * t299 + (t235 * t338 - t236 * t501 + t359 * t461) * t416) * t379 + ((-t231 * t460 + t237 * t438 + t238 * t439) * t299 + (t338 * t237 + t238 * t501 - t359 * t460) * t416) * t376 + (-t323 * t75 + t344) * t358 + (t323 * t336 + t343) * t357 + (-t41 * ((t225 * t326 + t226 * t325) * t416 + t391) - t40 * (-t226 * t299 + t389) - t39 * (t225 * t299 + t390) + t393 + t41 * t415 - t29 * t123 + t40 * t476 + (t29 * t109 + t40 * t87 + t30 * t453 + t39 * (-t84 - t88)) * t323 + ((t10 * t108 + t41 * t88 + t29 * t443 - t40 * t194 + (t109 * t41 + t39 * t443) * qJD(1)) * t326 + (t10 * t109 + t41 * t87 + t30 * t443 - t39 * t194 + (t186 * t40 + t41 * t453) * qJD(1)) * t325) * t321) * m(6) + ((-t156 * t45 - t157 * t44 - t61 * t96 + t62 * t95) * t323 + (t33 * t360 + t74 * (-t156 * t421 - t157 * t420 + t363) + t364 * t213 + ((t325 * t62 - t326 * t61) * qJD(1) + t366) * t204) * t321 - (t181 * t61 - t182 * t62) * t299 - (t74 * (t181 * t326 + t182 * t325) + t364 * t347) * t416) * m(5); t335 + (t393 + t29 * (-t199 * t460 - t123) + t30 * (-t135 * t323 - t199 * t461) + (-t135 * t405 - t391 + t415) * t41 + (-t194 * t460 - t389 + t476) * t40 + (-t323 * t84 + (-t194 * t325 - t199 * t420) * t321 - t390) * t39) * m(6);];
tauc = t1(:);
