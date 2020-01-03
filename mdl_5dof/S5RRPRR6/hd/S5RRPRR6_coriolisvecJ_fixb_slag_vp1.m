% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:30
% EndTime: 2020-01-03 12:05:49
% DurationCPUTime: 12.83s
% Computational Cost: add. (24245->686), mult. (23385->935), div. (0->0), fcn. (22651->10), ass. (0->345)
t338 = qJ(1) + qJ(2);
t328 = sin(t338);
t330 = cos(t338);
t343 = cos(qJ(4));
t340 = cos(pkin(9));
t341 = sin(qJ(4));
t480 = t340 * t341;
t259 = -t328 * t343 + t330 * t480;
t479 = t340 * t343;
t489 = t328 * t341;
t260 = t330 * t479 + t489;
t339 = sin(pkin(9));
t487 = t330 * t339;
t167 = Icges(5,5) * t260 - Icges(5,6) * t259 + Icges(5,3) * t487;
t241 = Icges(5,4) * t260;
t171 = Icges(5,2) * t259 - Icges(5,6) * t487 - t241;
t240 = Icges(5,4) * t259;
t173 = Icges(5,1) * t260 + Icges(5,5) * t487 - t240;
t79 = t340 * t167 - (t171 * t341 + t173 * t343) * t339;
t337 = qJ(4) + qJ(5);
t327 = sin(t337);
t329 = cos(t337);
t486 = t330 * t340;
t236 = t327 * t486 - t328 * t329;
t237 = t327 * t328 + t329 * t486;
t149 = Icges(6,5) * t237 - Icges(6,6) * t236 + Icges(6,3) * t487;
t225 = Icges(6,4) * t237;
t153 = Icges(6,2) * t236 - Icges(6,6) * t487 - t225;
t224 = Icges(6,4) * t236;
t155 = Icges(6,1) * t237 + Icges(6,5) * t487 - t224;
t74 = t340 * t149 - (t153 * t327 + t155 * t329) * t339;
t390 = t171 * t259 + t173 * t260;
t72 = t167 * t487 + t390;
t555 = t330 * t72;
t336 = qJD(1) + qJD(2);
t534 = t328 * rSges(3,1) + t330 * rSges(3,2);
t249 = t534 * t336;
t342 = sin(qJ(1));
t512 = pkin(1) * qJD(1);
t443 = t342 * t512;
t232 = t443 + t249;
t335 = qJD(4) + qJD(5);
t420 = t335 * t339;
t252 = t328 * t420;
t253 = t330 * t420;
t277 = -t335 * t340 + t336;
t490 = t328 * t340;
t234 = -t327 * t490 - t329 * t330;
t235 = -t330 * t327 + t329 * t490;
t491 = t328 * t339;
t148 = Icges(6,5) * t235 + Icges(6,6) * t234 + Icges(6,3) * t491;
t500 = Icges(6,4) * t235;
t151 = Icges(6,2) * t234 + Icges(6,6) * t491 + t500;
t223 = Icges(6,4) * t234;
t154 = Icges(6,1) * t235 + Icges(6,5) * t491 + t223;
t63 = t148 * t491 + t234 * t151 + t235 * t154;
t554 = -t234 * t153 + t155 * t235;
t64 = -t149 * t491 - t554;
t228 = -Icges(6,3) * t340 + (Icges(6,5) * t329 - Icges(6,6) * t327) * t339;
t498 = Icges(6,4) * t329;
t229 = -Icges(6,6) * t340 + (-Icges(6,2) * t327 + t498) * t339;
t499 = Icges(6,4) * t327;
t230 = -Icges(6,5) * t340 + (Icges(6,1) * t329 - t499) * t339;
t99 = t228 * t491 + t229 * t234 + t230 * t235;
t24 = t252 * t63 - t253 * t64 + t99 * t277;
t257 = -t328 * t480 - t330 * t343;
t485 = t330 * t341;
t258 = t328 * t479 - t485;
t552 = -t257 * t171 + t173 * t258;
t66 = t149 * t487 + t153 * t236 + t155 * t237;
t314 = -qJD(4) * t340 + t336;
t246 = -Icges(5,3) * t340 + (Icges(5,5) * t343 - Icges(5,6) * t341) * t339;
t501 = Icges(5,4) * t343;
t247 = -Icges(5,6) * t340 + (-Icges(5,2) * t341 + t501) * t339;
t502 = Icges(5,4) * t341;
t248 = -Icges(5,5) * t340 + (Icges(5,1) * t343 - t502) * t339;
t367 = t246 * t487 - t247 * t259 + t248 * t260;
t546 = t367 * t314;
t324 = pkin(4) * t343 + pkin(3);
t545 = pkin(4) * t489 + t324 * t486;
t345 = -pkin(8) - pkin(7);
t481 = t339 * t345;
t437 = t336 * t481;
t511 = pkin(4) * qJD(4);
t544 = -t343 * t511 - t437;
t160 = rSges(6,1) * t237 - rSges(6,2) * t236 + rSges(6,3) * t487;
t513 = pkin(7) + t345;
t517 = pkin(3) * t340;
t176 = (t339 * t513 + t517) * t330 - t545;
t231 = -rSges(6,3) * t340 + (rSges(6,1) * t329 - rSges(6,2) * t327) * t339;
t514 = -pkin(3) + t324;
t222 = t339 * t514 + t340 * t513;
t447 = qJD(4) * t339;
t380 = (-t222 * t447 - qJD(3)) * t330;
t543 = t277 * t160 - t314 * t176 - t253 * t231 + t380;
t178 = rSges(5,1) * t260 - rSges(5,2) * t259 + rSges(5,3) * t487;
t251 = -rSges(5,3) * t340 + (rSges(5,1) * t343 - rSges(5,2) * t341) * t339;
t542 = t314 * t178 + (-t251 * t447 - qJD(3)) * t330;
t368 = t228 * t487 - t229 * t236 + t230 * t237;
t541 = t253 * t66 + t368 * t277;
t418 = t511 * t480;
t446 = pkin(4) * t485;
t435 = t328 * t544 - t336 * t446;
t492 = t328 * t336;
t515 = pkin(7) * t339;
t116 = t330 * t418 + (t340 * t514 - t515) * t492 + t435;
t434 = t328 * t447;
t205 = t222 * t434;
t213 = t336 * t252;
t261 = (-rSges(6,1) * t327 - rSges(6,2) * t329) * t339;
t218 = t335 * t261;
t332 = t342 * pkin(1);
t347 = qJD(1) ^ 2;
t445 = t347 * t332;
t448 = qJD(3) * t336;
t408 = t328 * t448 - t445;
t516 = pkin(4) * t341;
t421 = t339 ^ 2 * qJD(4) ^ 2 * t516;
t317 = qJD(3) * t328;
t488 = t330 * t336;
t451 = qJ(3) * t488 + t317;
t212 = pkin(2) * t492 - t451;
t404 = t515 + t517;
t467 = -t404 * t492 - t212;
t393 = t277 * t330;
t482 = t336 * t340;
t423 = -t335 + t482;
t396 = t327 * t423;
t162 = -t328 * t396 - t329 * t393;
t394 = t423 * t329;
t163 = -t327 * t393 + t328 * t394;
t399 = rSges(6,1) * t163 + rSges(6,2) * t162;
t483 = t336 * t339;
t440 = t328 * t483;
t97 = rSges(6,3) * t440 + t399;
t43 = t330 * t421 - t116 * t314 + t213 * t231 - t218 * t253 - t277 * t97 + (t205 + t467) * t336 + t408;
t373 = t257 * qJD(4);
t438 = t330 * t482;
t439 = t330 * t483;
t453 = pkin(3) * t438 + pkin(7) * t439;
t457 = t545 * t336;
t117 = pkin(4) * t373 - t330 * t437 - t453 + t457;
t214 = t336 * t253;
t344 = cos(qJ(1));
t333 = t344 * pkin(1);
t325 = t347 * t333;
t449 = qJD(3) * t330;
t452 = pkin(2) * t488 + qJ(3) * t492;
t402 = -t449 + t452;
t468 = t336 * t402 + t325;
t436 = t336 * t453 + t468;
t395 = t328 * t277;
t164 = t329 * t395 - t330 * t396;
t165 = t327 * t395 + t330 * t394;
t98 = t165 * rSges(6,1) + t164 * rSges(6,2) + rSges(6,3) * t439;
t44 = t117 * t314 - t214 * t231 - t218 * t252 + t277 * t98 + t328 * t421 + t336 * t380 + t436;
t406 = -t317 + t443;
t262 = pkin(3) * t490 + pkin(7) * t491;
t322 = t328 * pkin(2);
t278 = -qJ(3) * t330 + t322;
t459 = t262 + t278;
t537 = t336 * t459;
t370 = t406 + t537;
t159 = t235 * rSges(6,1) + t234 * rSges(6,2) + rSges(6,3) * t491;
t275 = t324 * t490;
t175 = -t328 * t481 - t262 + t275 - t446;
t528 = t159 * t277 + t175 * t314 - t252 * t231 - t205;
t58 = t370 + t528;
t326 = t344 * t512;
t263 = t404 * t330;
t280 = t330 * pkin(2) + t328 * qJ(3);
t458 = t263 + t280;
t391 = t336 * t458 + t326;
t59 = t391 + t543;
t540 = (-t58 * t418 - t44 * t481 + t59 * (-rSges(6,3) * t339 - t324 * t340 - pkin(2)) * t336) * t328 + (-t43 * t481 - t59 * t418 + t44 * (-qJ(3) - t516) + t58 * (-qJD(3) + t544)) * t330;
t85 = t391 + t542;
t539 = t85 * (-t517 - pkin(2) + (-rSges(5,3) - pkin(7)) * t339) * t492;
t307 = rSges(4,1) * t490;
t407 = -rSges(4,2) * t491 + t307;
t536 = t336 * (-rSges(4,3) * t330 + t278 + t407);
t535 = -t236 * t151 + t237 * t154;
t503 = Icges(5,4) * t258;
t169 = Icges(5,2) * t257 + Icges(5,6) * t491 + t503;
t239 = Icges(5,4) * t257;
t172 = Icges(5,1) * t258 + Icges(5,5) * t491 + t239;
t478 = -t259 * t169 + t260 * t172;
t269 = t336 * t280;
t462 = t336 * t263 + t269;
t533 = t452 + t457 - t462 - t543 + t98;
t192 = -qJD(4) * t258 - t259 * t336;
t193 = t260 * t336 + t373;
t113 = t193 * rSges(5,1) + t192 * rSges(5,2) + rSges(5,3) * t439;
t532 = t113 + t402 + t453 - t462 - t542;
t450 = rSges(4,1) * t486 + t328 * rSges(4,3);
t221 = -rSges(4,2) * t487 + t450;
t454 = rSges(4,1) * t438 + rSges(4,3) * t492;
t530 = -t336 * t221 - t269 + t449 + t452 + t454;
t177 = t258 * rSges(5,1) + t257 * rSges(5,2) + rSges(5,3) * t491;
t210 = t251 * t434;
t529 = t177 * t314 - t210;
t371 = t328 * (-Icges(5,2) * t258 + t172 + t239) - t330 * (Icges(5,2) * t260 - t173 + t240);
t527 = t328 * (-Icges(5,1) * t257 + t169 + t503) - t330 * (-Icges(5,1) * t259 + t171 - t241);
t255 = (-Icges(6,2) * t329 - t499) * t339;
t353 = t252 * (-Icges(6,2) * t235 + t154 + t223) - t253 * (Icges(6,2) * t237 - t155 + t224) + t277 * (t230 + t255);
t256 = (-Icges(6,1) * t327 - t498) * t339;
t526 = t252 * (-Icges(6,1) * t234 + t151 + t500) - t253 * (-Icges(6,1) * t236 + t153 - t225) + t277 * (t229 - t256);
t525 = t213 / 0.2e1;
t524 = t214 / 0.2e1;
t523 = -t252 / 0.2e1;
t522 = t252 / 0.2e1;
t521 = -t253 / 0.2e1;
t520 = t253 / 0.2e1;
t518 = -t340 / 0.2e1;
t92 = Icges(6,5) * t165 + Icges(6,6) * t164 + Icges(6,3) * t439;
t94 = Icges(6,4) * t165 + Icges(6,2) * t164 + Icges(6,6) * t439;
t96 = Icges(6,1) * t165 + Icges(6,4) * t164 + Icges(6,5) * t439;
t34 = -t340 * t92 + ((-t151 * t335 + t96) * t329 + (-t154 * t335 - t94) * t327) * t339;
t508 = t34 * t252;
t91 = Icges(6,5) * t163 + Icges(6,6) * t162 + Icges(6,3) * t440;
t93 = Icges(6,4) * t163 + Icges(6,2) * t162 + Icges(6,6) * t440;
t95 = Icges(6,1) * t163 + Icges(6,4) * t162 + Icges(6,5) * t440;
t35 = -t340 * t91 + ((-t153 * t335 + t95) * t329 + (t155 * t335 - t93) * t327) * t339;
t507 = t35 * t253;
t73 = -t148 * t340 + (-t151 * t327 + t154 * t329) * t339;
t506 = t73 * t214;
t505 = t74 * t213;
t504 = t231 * t440 + t340 * t97;
t166 = Icges(5,5) * t258 + Icges(5,6) * t257 + Icges(5,3) * t491;
t496 = t166 * t330;
t495 = t167 * t328;
t473 = -t159 - t175;
t465 = -t222 - t231;
t272 = (-Icges(5,1) * t341 - t501) * t339;
t461 = t247 - t272;
t271 = (-Icges(5,2) * t343 - t502) * t339;
t460 = t248 + t271;
t455 = rSges(4,2) * t440 + rSges(4,3) * t488;
t444 = -t160 * t439 + t98 * t487 + t97 * t491;
t69 = t166 * t491 + t257 * t169 + t258 * t172;
t70 = -t167 * t491 - t552;
t432 = t491 / 0.2e1;
t431 = -t487 / 0.2e1;
t430 = t483 / 0.2e1;
t429 = -t447 / 0.2e1;
t428 = t447 / 0.2e1;
t15 = -t159 * t213 - t160 * t214 + t252 * t97 + t253 * t98 + ((t176 * t336 + t117) * t330 + (-t175 * t336 + t116) * t328) * t447;
t427 = t15 * (t159 * t487 - t160 * t491);
t186 = rSges(6,1) * t234 - rSges(6,2) * t235;
t187 = rSges(6,1) * t236 + rSges(6,2) * t237;
t426 = t186 * t253 + t252 * t187;
t425 = t277 * t186 - t252 * t261;
t424 = -t187 * t277 - t253 * t261;
t281 = rSges(3,1) * t330 - rSges(3,2) * t328;
t233 = t281 * t336 + t326;
t416 = t280 + t450;
t415 = t328 * t430;
t414 = t330 * t430;
t413 = t328 * t429;
t412 = t328 * t428;
t411 = t330 * t429;
t410 = t330 * t428;
t409 = t336 * t428;
t250 = rSges(3,1) * t488 - rSges(3,2) * t492;
t405 = -rSges(4,2) * t483 - qJD(3);
t190 = qJD(4) * t260 + t257 * t336;
t191 = qJD(4) * t259 + t258 * t336;
t400 = rSges(5,1) * t191 + rSges(5,2) * t190;
t84 = t370 + t529;
t398 = -t328 * t84 - t330 * t85;
t397 = t275 + t322 + t159;
t392 = t317 - t537;
t389 = (Icges(5,5) * t257 - Icges(5,6) * t258) * t328 - (Icges(5,5) * t259 + Icges(5,6) * t260) * t330;
t388 = t328 * t409;
t387 = t330 * t409;
t386 = t322 + t407;
t65 = -t148 * t487 - t535;
t381 = (-rSges(5,1) * t341 - rSges(5,2) * t343) * t339;
t378 = (t328 * t69 - t330 * t70) * t339;
t71 = -t166 * t487 - t478;
t377 = (t328 * t71 - t555) * t339;
t270 = (-Icges(5,5) * t341 - Icges(5,6) * t343) * t339;
t254 = (-Icges(6,5) * t327 - Icges(6,6) * t329) * t339;
t376 = -t400 + t451;
t374 = (Icges(6,5) * t234 - Icges(6,6) * t235) * t252 - (Icges(6,5) * t236 + Icges(6,6) * t237) * t253 + t254 * t277;
t104 = (t177 * t330 - t178 * t328) * t447;
t369 = t177 + t459;
t364 = -t399 - t435 + t451;
t363 = t160 + t280 + t545;
t360 = (-rSges(4,1) * t340 - pkin(2)) * t492 + t451 + t455;
t16 = t151 * t162 + t154 * t163 + t236 * t94 - t237 * t96 + (t148 * t492 - t330 * t92) * t339;
t17 = t153 * t162 - t155 * t163 + t236 * t93 - t237 * t95 + (-t149 * t492 - t330 * t91) * t339;
t18 = t151 * t164 + t154 * t165 + t234 * t94 + t235 * t96 + (t148 * t488 + t328 * t92) * t339;
t19 = t153 * t164 - t155 * t165 + t234 * t93 + t235 * t95 + (-t149 * t488 + t328 * t91) * t339;
t25 = t252 * t65 - t541;
t215 = t335 * t254;
t216 = t335 * t255;
t217 = t335 * t256;
t50 = t162 * t229 + t163 * t230 + t216 * t236 - t217 * t237 + (-t215 * t330 + t228 * t492) * t339;
t51 = t164 * t229 + t165 * t230 + t216 * t234 + t217 * t235 + (t215 * t328 + t228 * t488) * t339;
t80 = -t215 * t340 + ((-t229 * t335 + t217) * t329 + (-t230 * t335 - t216) * t327) * t339;
t77 = t80 * t277;
t359 = (t18 * t252 - t19 * t253 + t213 * t64 + t214 * t63 + t277 * t51) * t432 + (t234 * t353 - t235 * t526 + t374 * t491) * t523 + (t353 * t236 + t237 * t526 - t374 * t487) * t520 - (-t374 * t340 + (-t327 * t353 - t329 * t526) * t339) * t277 / 0.2e1 + (t16 * t252 - t17 * t253 + t213 * t66 + t214 * t65 + t277 * t50) * t431 + t25 * t415 + t24 * t414 + (t368 * t340 + (t328 * t65 - t330 * t66) * t339) * t525 + (-t340 * t99 + (t328 * t63 - t330 * t64) * t339) * t524 + (-t340 * t51 + ((t336 * t63 - t19) * t330 + (t336 * t64 + t18) * t328) * t339) * t522 + (-t340 * t50 + ((t336 * t65 - t17) * t330 + (t336 * t66 + t16) * t328) * t339) * t521 + (t505 + t506 - t507 + t77 + t508) * t518 + t277 * (-t340 * t80 + ((t336 * t73 - t35) * t330 + (t336 * t74 + t34) * t328) * t339) / 0.2e1;
t106 = Icges(5,5) * t191 + Icges(5,6) * t190 + Icges(5,3) * t440;
t107 = Icges(5,5) * t193 + Icges(5,6) * t192 + Icges(5,3) * t439;
t108 = Icges(5,4) * t191 + Icges(5,2) * t190 + Icges(5,6) * t440;
t109 = Icges(5,4) * t193 + Icges(5,2) * t192 + Icges(5,6) * t439;
t110 = Icges(5,1) * t191 + Icges(5,4) * t190 + Icges(5,5) * t440;
t111 = Icges(5,1) * t193 + Icges(5,4) * t192 + Icges(5,5) * t439;
t357 = ((t336 * t71 - t108 * t259 + t110 * t260 - t171 * t190 + t173 * t191 - (-t106 * t330 - t167 * t492) * t339) * t330 + (t336 * t72 + t109 * t259 - t111 * t260 + t169 * t190 + t172 * t191 + (-t107 * t330 + t166 * t492) * t339) * t328) * t339;
t356 = ((t336 * t69 - t108 * t257 - t110 * t258 - t171 * t192 + t173 * t193 - (t106 * t328 - t167 * t488) * t339) * t330 + (t336 * t70 + t109 * t257 + t111 * t258 + t169 * t192 + t172 * t193 + (t107 * t328 + t166 * t488) * t339) * t328) * t339;
t45 = -t107 * t340 + (-t109 * t341 + t111 * t343 + (-t169 * t343 - t172 * t341) * qJD(4)) * t339;
t46 = -t106 * t340 + (-t108 * t341 + t110 * t343 + (-t171 * t343 + t173 * t341) * qJD(4)) * t339;
t78 = -t166 * t340 + (-t169 * t341 + t172 * t343) * t339;
t355 = ((t336 * t78 - t46) * t330 + (t336 * t79 + t45) * t328) * t339;
t351 = t178 + t458;
t118 = (-t307 * t336 - t212 + t455) * t336 + t408;
t119 = (t330 * t405 + t454) * t336 + t468;
t157 = t406 + t536;
t350 = (-t118 * rSges(4,2) * t339 + t119 * (-rSges(4,3) - qJ(3)) + t157 * t405) * t330;
t114 = t246 * t491 + t247 * t257 + t248 * t258;
t105 = t114 * t314;
t39 = qJD(4) * t378 + t105;
t40 = qJD(4) * t377 - t546;
t264 = qJD(4) * t270;
t265 = qJD(4) * t271;
t266 = qJD(4) * t272;
t56 = t190 * t247 + t191 * t248 + t259 * t265 - t260 * t266 + (t246 * t492 - t264 * t330) * t339;
t57 = t192 * t247 + t193 * t248 + t257 * t265 + t258 * t266 + (t246 * t488 + t264 * t328) * t339;
t101 = -t264 * t340 + (-t265 * t341 + t266 * t343 + (-t247 * t343 - t248 * t341) * qJD(4)) * t339;
t90 = t101 * t314;
t349 = t508 / 0.2e1 + t77 + t51 * t522 + t99 * t524 + (t105 + ((-t390 + t69 + t72) * t328 + (t71 + (-t495 + t496) * t339 - t70 + t478) * t330) * t447) * t410 - t368 * t525 - t507 / 0.2e1 + t505 / 0.2e1 + t506 / 0.2e1 + t90 + t24 * t520 + (t25 + (t64 + (t148 * t330 + t149 * t328) * t339 + t535 + t554) * t252 + t541) * t523 + (t50 + t24) * t521 + (t45 + t57) * t412 + (-t367 + t79) * t388 + (t114 + t78) * t387 + (t40 + t546 + (t555 + (t478 + t70 + (t495 + t496) * t339 + t552) * t328) * t447) * t413 + (t46 + t56 + t39) * t411;
t267 = qJD(4) * t381;
t243 = t259 * pkin(4);
t242 = t257 * pkin(4);
t208 = t250 * t336 + t325;
t207 = -t249 * t336 - t445;
t202 = rSges(5,1) * t259 + rSges(5,2) * t260;
t201 = rSges(5,1) * t257 - rSges(5,2) * t258;
t158 = -t449 + t326 + (t221 + t280) * t336;
t141 = t340 * t160;
t112 = rSges(5,3) * t440 + t400;
t61 = -t330 * t448 + t113 * t314 + (-t251 * t488 - t267 * t328) * t447 + t436;
t60 = -t267 * t330 * t447 - t112 * t314 + (t210 + t467) * t336 + t408;
t54 = t159 * t253 - t160 * t252 + (t175 * t330 + t176 * t328) * t447;
t1 = [t349 + m(3) * (t207 * (t281 + t333) + t208 * (t332 + t534) + (-t233 + t250 + t326) * t232) + (t43 * (t333 + t363) + t59 * (t364 - t443) + t44 * (t332 + t397) + (t59 + t533) * t58 + t540) * m(6) + (t60 * (t333 + t351) + t85 * (t376 - t443) + t61 * (t332 + t369) + (t85 + t532) * t84 + t539) * m(5) + (t118 * (t333 + t416) + t158 * (t360 - t443) + t119 * (t332 + t386) + t350 + (t158 + t530) * t157) * m(4); t349 + (t43 * t363 + t44 * t397 + (-t392 + t364 + t528) * t59 + t533 * t58 + t540) * m(6) + (t60 * t351 + t61 * t369 + (-t392 + t376 + t529) * t85 + t532 * t84 + t539) * m(5) + (t118 * t416 + t119 * t386 + t350 + (-t317 + t360 + t536) * t158 + t530 * t157) * m(4) + (t207 * t281 + t208 * t534 + t232 * t250 - t233 * t249 - (t232 * t281 - t233 * t534) * t336) * m(3); m(4) * (-t118 * t330 - t119 * t328) + m(5) * (-t61 * t328 - t60 * t330) + m(6) * (-t328 * t44 - t330 * t43); ((t259 * t460 + t260 * t461 - t270 * t487) * t314 + (t371 * t259 + t260 * t527 - t389 * t487) * t447) * t410 + ((t257 * t460 - t258 * t461 + t270 * t491) * t314 + (t257 * t371 - t258 * t527 + t389 * t491) * t447) * t413 + (t340 * t367 + t377) * t388 + (-t114 * t340 + t378) * t387 + t40 * t415 + t39 * t414 + (qJD(4) * t356 + t314 * t57) * t432 + (qJD(4) * t355 + t90) * t518 + t314 * (-t101 * t340 + t355) / 0.2e1 + (-t340 * t57 + t356) * t412 + (-t340 * t56 + t357) * t411 + (qJD(4) * t357 + t314 * t56) * t431 - t314 * (-t270 * t340 * t314 + ((-t341 * t460 - t343 * t461) * t314 + ((-t341 * t371 - t343 * t527) * t339 - t389 * t340) * qJD(4)) * t339) / 0.2e1 + t359 + (-t54 * ((t242 * t330 + t243 * t328) * t447 + t426) - t59 * (-t243 * t314 + t424) - t58 * (t242 * t314 + t425) + t427 + t54 * t444 - t43 * t141 + t59 * t504 + (t43 * t176 + t59 * t116 + t44 * t473 + t58 * (-t117 - t98)) * t340 + ((t15 * t175 + t54 * t117 + t43 * t465 - t59 * t218 + (t54 * t176 + t465 * t58) * t336) * t330 + (t15 * t176 + t54 * t116 + t44 * t465 - t58 * t218 + (t59 * t222 + t473 * t54) * t336) * t328) * t339) * m(6) + ((t112 * t85 - t113 * t84 - t177 * t61 - t178 * t60) * t340 + (0.2e1 * t104 * ((-t178 * t336 + t113) * t330 + (-t177 * t336 + t112) * t328) + t398 * t267 + ((-t336 * t84 - t60) * t330 + (t336 * t85 - t61) * t328) * t251) * t339 - (t201 * t84 - t202 * t85) * t314 - (t104 * (t201 * t330 + t202 * t328) + t398 * t381) * t447) * m(5); t359 + (t427 + t43 * (-t231 * t487 - t141) + t44 * (-t159 * t340 - t231 * t491) + (-t218 * t487 - t424 + t504) * t59 + (-t425 - t340 * t98 + (-t218 * t328 - t231 * t488) * t339) * t58 + (-t159 * t440 - t426 + t444) * t54) * m(6);];
tauc = t1(:);
