% Calculate time derivative of joint inertia matrix for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:31
% EndTime: 2019-03-09 03:04:57
% DurationCPUTime: 15.85s
% Computational Cost: add. (29812->791), mult. (29290->1086), div. (0->0), fcn. (27832->10), ass. (0->388)
t543 = Icges(4,3) + Icges(5,3);
t286 = qJ(3) + pkin(10);
t281 = sin(t286);
t283 = cos(t286);
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t542 = Icges(4,5) * t293 + Icges(5,5) * t283 - Icges(4,6) * t290 - Icges(5,6) * t281;
t287 = qJ(1) + pkin(9);
t282 = sin(t287);
t284 = cos(t287);
t535 = t282 * t543 + t542 * t284;
t289 = sin(qJ(5));
t292 = cos(qJ(5));
t341 = Icges(6,5) * t292 - Icges(6,6) * t289;
t205 = -Icges(6,3) * t283 + t281 * t341;
t344 = Icges(7,4) * t292 + Icges(7,6) * t289;
t208 = -Icges(7,2) * t283 + t281 * t344;
t541 = t205 + t208;
t433 = t284 * t292;
t438 = t282 * t289;
t227 = t283 * t438 + t433;
t434 = t284 * t289;
t437 = t282 * t292;
t228 = t283 * t437 - t434;
t529 = rSges(7,3) + qJ(6);
t530 = rSges(7,1) + pkin(5);
t540 = -t529 * t227 - t228 * t530;
t413 = qJD(3) * t289;
t390 = t281 * t413;
t407 = qJD(5) * t292;
t408 = qJD(5) * t289;
t419 = qJD(1) * t284;
t421 = qJD(1) * t282;
t126 = -t282 * t390 - t284 * t408 - t292 * t421 + (t282 * t407 + t289 * t419) * t283;
t411 = qJD(3) * t292;
t299 = -t281 * t411 + (-qJD(5) * t283 + qJD(1)) * t289;
t420 = qJD(1) * t283;
t376 = -qJD(5) + t420;
t127 = t282 * t299 + t376 * t433;
t539 = -t227 * qJD(6) - t529 * t126 - t127 * t530;
t462 = Icges(5,4) * t283;
t347 = -Icges(5,2) * t281 + t462;
t194 = Icges(5,6) * t282 + t284 * t347;
t463 = Icges(5,4) * t281;
t353 = Icges(5,1) * t283 - t463;
t196 = Icges(5,5) * t282 + t284 * t353;
t464 = Icges(4,4) * t293;
t349 = -Icges(4,2) * t290 + t464;
t211 = Icges(4,6) * t282 + t284 * t349;
t465 = Icges(4,4) * t290;
t355 = Icges(4,1) * t293 - t465;
t215 = Icges(4,5) * t282 + t284 * t355;
t504 = t194 * t281 - t196 * t283 + t211 * t290 - t215 * t293;
t538 = t282 * t504;
t537 = t504 * t284;
t536 = t542 * t282 - t284 * t543;
t534 = (-Icges(4,5) * t290 - Icges(5,5) * t281 - Icges(4,6) * t293 - Icges(5,6) * t283) * qJD(3);
t193 = -Icges(5,6) * t284 + t282 * t347;
t195 = -Icges(5,5) * t284 + t282 * t353;
t210 = -Icges(4,6) * t284 + t282 * t349;
t214 = -Icges(4,5) * t284 + t282 * t355;
t533 = t193 * t281 - t195 * t283 + t210 * t290 - t214 * t293;
t442 = t281 * t282;
t143 = Icges(6,5) * t228 - Icges(6,6) * t227 + Icges(6,3) * t442;
t147 = Icges(6,4) * t228 - Icges(6,2) * t227 + Icges(6,6) * t442;
t151 = Icges(6,1) * t228 - Icges(6,4) * t227 + Icges(6,5) * t442;
t229 = t283 * t434 - t437;
t399 = t283 * t433;
t230 = t399 + t438;
t441 = t281 * t284;
t59 = t143 * t441 - t147 * t229 + t151 * t230;
t144 = Icges(6,5) * t230 - Icges(6,6) * t229 + Icges(6,3) * t441;
t148 = Icges(6,4) * t230 - Icges(6,2) * t229 + Icges(6,6) * t441;
t152 = Icges(6,1) * t230 - Icges(6,4) * t229 + Icges(6,5) * t441;
t60 = t144 * t441 - t148 * t229 + t152 * t230;
t356 = t282 * t59 + t284 * t60;
t141 = Icges(7,5) * t228 + Icges(7,6) * t442 + Icges(7,3) * t227;
t145 = Icges(7,4) * t228 + Icges(7,2) * t442 + Icges(7,6) * t227;
t149 = Icges(7,1) * t228 + Icges(7,4) * t442 + Icges(7,5) * t227;
t57 = t141 * t229 + t145 * t441 + t149 * t230;
t142 = Icges(7,5) * t230 + Icges(7,6) * t441 + Icges(7,3) * t229;
t146 = Icges(7,4) * t230 + Icges(7,2) * t441 + Icges(7,6) * t229;
t150 = Icges(7,1) * t230 + Icges(7,4) * t441 + Icges(7,5) * t229;
t58 = t142 * t229 + t146 * t441 + t150 * t230;
t357 = t282 * t57 + t284 * t58;
t532 = t356 + t357;
t55 = t143 * t442 - t147 * t227 + t151 * t228;
t56 = t144 * t442 - t148 * t227 + t152 * t228;
t358 = t282 * t55 + t284 * t56;
t53 = t141 * t227 + t145 * t442 + t149 * t228;
t54 = t142 * t227 + t146 * t442 + t150 * t228;
t359 = t282 * t53 + t284 * t54;
t531 = t358 + t359;
t456 = Icges(7,5) * t292;
t340 = Icges(7,3) * t289 + t456;
t204 = -Icges(7,6) * t283 + t281 * t340;
t460 = Icges(6,4) * t292;
t345 = -Icges(6,2) * t289 + t460;
t209 = -Icges(6,6) * t283 + t281 * t345;
t457 = Icges(7,5) * t289;
t350 = Icges(7,1) * t292 + t457;
t212 = -Icges(7,4) * t283 + t281 * t350;
t461 = Icges(6,4) * t289;
t351 = Icges(6,1) * t292 - t461;
t213 = -Icges(6,5) * t283 + t281 * t351;
t528 = t541 * t283 + ((-t212 - t213) * t292 + (-t204 + t209) * t289) * t281;
t387 = t283 * t413;
t527 = -t281 * t407 - t387;
t526 = t535 * qJD(1);
t525 = t536 * t282 - t538 + (-t533 - t535) * t284;
t336 = t142 * t289 + t150 * t292;
t62 = -t146 * t283 + t281 * t336;
t334 = -t148 * t289 + t152 * t292;
t64 = -t144 * t283 + t281 * t334;
t474 = t62 + t64;
t337 = t141 * t289 + t149 * t292;
t61 = -t145 * t283 + t281 * t337;
t335 = -t147 * t289 + t151 * t292;
t63 = -t143 * t283 + t281 * t335;
t475 = t61 + t63;
t524 = t282 * t475 + t284 * t474;
t496 = t282 ^ 2;
t495 = t284 ^ 2;
t124 = qJD(1) * t227 - qJD(5) * t399 - t282 * t408 + t284 * t390;
t125 = t284 * t299 - t376 * t437;
t414 = qJD(3) * t284;
t388 = t283 * t414;
t392 = t281 * t421;
t302 = t388 - t392;
t415 = qJD(3) * t283;
t389 = t282 * t415;
t303 = t281 * t419 + t389;
t69 = Icges(7,5) * t127 + Icges(7,6) * t303 + Icges(7,3) * t126;
t73 = Icges(7,4) * t127 + Icges(7,2) * t303 + Icges(7,6) * t126;
t77 = Icges(7,1) * t127 + Icges(7,4) * t303 + Icges(7,5) * t126;
t11 = -t124 * t141 + t125 * t149 + t145 * t302 + t229 * t69 + t230 * t77 + t441 * t73;
t68 = Icges(7,5) * t125 + Icges(7,6) * t302 - Icges(7,3) * t124;
t72 = Icges(7,4) * t125 + Icges(7,2) * t302 - Icges(7,6) * t124;
t76 = Icges(7,1) * t125 + Icges(7,4) * t302 - Icges(7,5) * t124;
t12 = -t124 * t142 + t125 * t150 + t146 * t302 + t229 * t68 + t230 * t76 + t441 * t72;
t71 = Icges(6,5) * t127 - Icges(6,6) * t126 + Icges(6,3) * t303;
t75 = Icges(6,4) * t127 - Icges(6,2) * t126 + Icges(6,6) * t303;
t79 = Icges(6,1) * t127 - Icges(6,4) * t126 + Icges(6,5) * t303;
t13 = t124 * t147 + t125 * t151 + t143 * t302 - t229 * t75 + t230 * t79 + t441 * t71;
t70 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t302;
t74 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t302;
t78 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t302;
t14 = t124 * t148 + t125 * t152 + t144 * t302 - t229 * t74 + t230 * t78 + t441 * t70;
t523 = (-t11 - t13) * t284 + (t12 + t14) * t282 + t532 * qJD(1);
t15 = t126 * t141 + t127 * t149 + t145 * t303 + t227 * t69 + t228 * t77 + t442 * t73;
t16 = t126 * t142 + t127 * t150 + t146 * t303 + t227 * t68 + t228 * t76 + t442 * t72;
t17 = -t126 * t147 + t127 * t151 + t143 * t303 - t227 * t75 + t228 * t79 + t442 * t71;
t18 = -t126 * t148 + t127 * t152 + t144 * t303 - t227 * t74 + t228 * t78 + t442 * t70;
t522 = (-t15 - t17) * t284 + (t16 + t18) * t282 + t531 * qJD(1);
t19 = (qJD(3) * t337 - t73) * t283 + (qJD(3) * t145 + t289 * t69 + t292 * t77 + (t141 * t292 - t149 * t289) * qJD(5)) * t281;
t21 = (qJD(3) * t335 - t71) * t283 + (qJD(3) * t143 - t289 * t75 + t292 * t79 + (-t147 * t292 - t151 * t289) * qJD(5)) * t281;
t521 = -t19 - t21;
t20 = (qJD(3) * t336 - t72) * t283 + (qJD(3) * t146 + t289 * t68 + t292 * t76 + (t142 * t292 - t150 * t289) * qJD(5)) * t281;
t22 = (qJD(3) * t334 - t70) * t283 + (qJD(3) * t144 - t289 * t74 + t292 * t78 + (-t148 * t292 - t152 * t289) * qJD(5)) * t281;
t520 = t20 + t22;
t91 = t204 * t227 + t208 * t442 + t212 * t228;
t92 = t205 * t442 - t209 * t227 + t213 * t228;
t519 = (-t91 - t92) * t283 + t531 * t281;
t93 = t204 * t229 + t208 * t441 + t212 * t230;
t94 = t205 * t441 - t209 * t229 + t213 * t230;
t518 = (-t93 - t94) * t283 + t532 * t281;
t280 = pkin(3) * t293 + pkin(2);
t259 = t284 * t280;
t288 = -qJ(4) - pkin(7);
t478 = -pkin(7) - t288;
t188 = -pkin(2) * t284 + t282 * t478 + t259;
t470 = rSges(5,1) * t283;
t366 = -rSges(5,2) * t281 + t470;
t197 = -rSges(5,3) * t284 + t282 * t366;
t436 = t283 * t284;
t277 = t282 * rSges(5,3);
t506 = -rSges(5,2) * t441 + t277;
t198 = rSges(5,1) * t436 + t506;
t243 = rSges(5,1) * t281 + rSges(5,2) * t283;
t314 = qJD(3) * t243;
t279 = t284 * pkin(7);
t435 = t284 * t288;
t479 = pkin(2) - t280;
t509 = t282 * t479;
t187 = t279 + t435 - t509;
t275 = qJD(4) * t282;
t412 = qJD(3) * t290;
t403 = pkin(3) * t412;
t394 = qJD(4) * t284 + t282 * t403 + t288 * t421;
t482 = pkin(7) * t282;
t396 = t282 * ((-t284 * t479 - t482) * qJD(1) - t394) + t284 * (-t284 * t403 + t275 + (t284 * t478 + t509) * qJD(1)) + t187 * t419;
t425 = rSges(5,2) * t392 + rSges(5,3) * t419;
t39 = (qJD(1) * t197 - t284 * t314 + t425) * t284 + (-t282 * t314 + (-t188 - t198 + t506) * qJD(1)) * t282 + t396;
t499 = 2 * m(5);
t517 = t39 * t499;
t500 = 2 * m(4);
t469 = rSges(4,2) * t290;
t471 = rSges(4,1) * t293;
t367 = -t469 + t471;
t468 = rSges(4,3) * t284;
t218 = t282 * t367 - t468;
t405 = t284 * t469;
t278 = t282 * rSges(4,3);
t424 = t284 * t471 + t278;
t219 = -t405 + t424;
t270 = rSges(4,1) * t290 + rSges(4,2) * t293;
t315 = qJD(3) * t270;
t418 = qJD(1) * t290;
t391 = t282 * t418;
t297 = rSges(4,2) * t391 + rSges(4,3) * t419 - t284 * t315;
t67 = (qJD(1) * t218 + t297) * t284 + (-t282 * t315 + (-t219 - t405 + t278) * qJD(1)) * t282;
t516 = t500 * t67;
t515 = t282 * t533 + t284 * t536;
t409 = qJD(5) * t281;
t162 = (Icges(7,3) * t292 - t457) * t409 + (Icges(7,6) * t281 + t283 * t340) * qJD(3);
t166 = (-Icges(7,4) * t289 + Icges(7,6) * t292) * t409 + (Icges(7,2) * t281 + t283 * t344) * qJD(3);
t170 = (-Icges(7,1) * t289 + t456) * t409 + (Icges(7,4) * t281 + t283 * t350) * qJD(3);
t171 = (-Icges(6,1) * t289 - t460) * t409 + (Icges(6,5) * t281 + t283 * t351) * qJD(3);
t385 = t281 * t408;
t386 = t283 * t411;
t417 = qJD(3) * t281;
t440 = t281 * t289;
t512 = t213 * t386 + t162 * t440 - t283 * t166 + (-t385 + t386) * t212 - t527 * t204 + (t171 + t170) * t281 * t292 + t541 * t417;
t511 = rSges(7,2) * t388 + qJD(6) * t229 - t529 * t124 + t125 * t530;
t510 = t282 * t535 - t537;
t508 = -qJD(1) * t536 + t284 * t534;
t507 = -t282 * t534 - t526;
t486 = sin(qJ(1)) * pkin(1);
t505 = t279 - t486;
t502 = t283 * t474 - t518;
t431 = rSges(7,2) * t441 + t529 * t229 + t230 * t530;
t432 = rSges(7,2) * t442 - t540;
t501 = -t282 * t432 - t284 * t431;
t498 = 2 * m(6);
t497 = 2 * m(7);
t494 = t282 / 0.2e1;
t489 = -rSges(7,2) - pkin(8);
t488 = -rSges(6,3) - pkin(8);
t487 = m(4) * t270;
t485 = pkin(3) * t290;
t484 = pkin(4) * t281;
t483 = pkin(4) * t283;
t285 = cos(qJ(1)) * pkin(1);
t163 = (-Icges(6,5) * t289 - Icges(6,6) * t292) * t409 + (Icges(6,3) * t281 + t283 * t341) * qJD(3);
t167 = (-Icges(6,2) * t292 - t461) * t409 + (Icges(6,6) * t281 + t283 * t345) * qJD(3);
t476 = (-t209 * t413 - t163) * t283 + (-t167 * t289 + (-t209 * t292 - t213 * t289) * qJD(5)) * t281 + t512;
t473 = -rSges(7,2) * t392 + t511;
t472 = rSges(7,2) * t303 - t539;
t466 = t528 * t417;
t364 = -rSges(6,1) * t228 + rSges(6,2) * t227;
t154 = rSges(6,3) * t442 - t364;
t451 = t154 * t284;
t450 = t193 * t283;
t449 = t194 * t283;
t448 = t195 * t281;
t447 = t196 * t281;
t446 = t210 * t293;
t445 = t211 * t293;
t444 = t214 * t290;
t443 = t215 * t290;
t360 = pkin(5) * t292 + qJ(6) * t289;
t362 = rSges(7,1) * t292 + rSges(7,3) * t289;
t430 = t360 * t415 + (qJD(6) * t289 + (-pkin(5) * t289 + qJ(6) * t292) * qJD(5)) * t281 + (-rSges(7,1) * t289 + rSges(7,3) * t292) * t409 + (rSges(7,2) * t281 + t283 * t362) * qJD(3);
t429 = t282 * t187 + t284 * t188;
t258 = pkin(4) * t436;
t226 = pkin(8) * t441 + t258;
t428 = -t188 - t226;
t427 = -rSges(7,2) * t283 + (t360 + t362) * t281;
t244 = -pkin(8) * t283 + t484;
t262 = pkin(3) * t391;
t426 = t244 * t421 + t262;
t416 = qJD(3) * t282;
t410 = qJD(3) * t293;
t402 = pkin(3) * t410;
t31 = t282 * t54 - t284 * t53;
t32 = t282 * t56 - t284 * t55;
t401 = t32 / 0.2e1 + t31 / 0.2e1;
t33 = t282 * t58 - t284 * t57;
t34 = t282 * t60 - t284 * t59;
t400 = t33 / 0.2e1 + t34 / 0.2e1;
t398 = t125 * rSges(6,1) + t124 * rSges(6,2) + rSges(6,3) * t388;
t156 = t230 * rSges(6,1) - t229 * rSges(6,2) + rSges(6,3) * t441;
t363 = rSges(6,1) * t292 - rSges(6,2) * t289;
t217 = -rSges(6,3) * t283 + t281 * t363;
t393 = t217 * t421;
t383 = -t243 - t485;
t382 = -t244 - t485;
t381 = -t280 - t483;
t380 = t432 * t284;
t379 = t427 * t284;
t378 = t427 * t282;
t377 = qJD(1) * t427;
t371 = pkin(8) * t281 + t483;
t225 = t371 * t282;
t375 = t282 * t225 + t284 * t226 + t429;
t248 = t416 * t484;
t374 = t248 + t394;
t373 = -t217 + t382;
t372 = -t371 * qJD(3) - t402;
t370 = -t282 * t288 + t259 + t285;
t369 = t282 * t377;
t201 = t383 * t284;
t365 = t127 * rSges(6,1) - t126 * rSges(6,2);
t361 = -t435 - t486;
t354 = Icges(4,1) * t290 + t464;
t352 = Icges(5,1) * t281 + t462;
t348 = Icges(4,2) * t293 + t465;
t346 = Icges(5,2) * t283 + t463;
t333 = -t156 * t282 + t451;
t332 = -t154 * t282 - t156 * t284;
t325 = t382 - t427;
t175 = (-rSges(6,1) * t289 - rSges(6,2) * t292) * t409 + (rSges(6,3) * t281 + t283 * t363) * qJD(3);
t324 = -t175 + t372;
t323 = -t283 * t475 + t519;
t322 = -pkin(2) - t367;
t321 = t63 / 0.2e1 + t61 / 0.2e1 + t91 / 0.2e1 + t92 / 0.2e1;
t320 = -t93 / 0.2e1 - t94 / 0.2e1 - t64 / 0.2e1 - t62 / 0.2e1;
t137 = t373 * t284;
t319 = -t280 - t366;
t249 = pkin(8) * t388;
t317 = t282 * (pkin(8) * t303 + qJD(1) * t258 - t248) + t284 * (-pkin(8) * t392 + t249 + (-t281 * t414 - t282 * t420) * pkin(4)) + t225 * t419 + t396;
t316 = t372 - t430;
t313 = t370 + t226;
t312 = qJD(3) * t354;
t311 = qJD(3) * t352;
t310 = qJD(3) * t348;
t309 = qJD(3) * t346;
t114 = t325 * t284;
t306 = t281 * t489 + t381;
t305 = t281 * t488 + t381;
t300 = -t282 * t431 + t380;
t298 = t249 + t275 + (-t484 - t485) * t414;
t296 = t282 * t306 + t361;
t295 = t282 * t305 + t361;
t253 = t367 * qJD(3);
t236 = t366 * qJD(3);
t200 = t383 * t282;
t179 = t482 + t285 + (pkin(2) - t469) * t284 + t424;
t178 = t282 * t322 + t468 + t505;
t160 = t198 + t370;
t159 = -t486 + (rSges(5,3) - t288) * t284 + t319 * t282;
t136 = t373 * t282;
t122 = -t243 * t419 - t236 * t282 + (-t282 * t410 - t284 * t418) * pkin(3);
t121 = t243 * t421 + t262 + (-t236 - t402) * t284;
t116 = t270 * t416 + (-t285 + (-rSges(4,3) - pkin(7)) * t282 + t322 * t284) * qJD(1);
t115 = ((-pkin(2) - t471) * t282 + t505) * qJD(1) + t297;
t113 = t325 * t282;
t106 = -t156 * t283 - t217 * t441;
t105 = t154 * t283 + t217 * t442;
t104 = t243 * t416 + (t284 * t319 - t277 - t285) * qJD(1) + t394;
t103 = t275 + qJD(3) * t201 + ((-t280 - t470) * t282 + t361) * qJD(1) + t425;
t102 = t313 + t156;
t101 = t295 + t364;
t90 = t333 * t281;
t89 = t313 + t431;
t88 = t296 + t540;
t87 = qJD(1) * t137 + t282 * t324;
t86 = t284 * t324 + t393 + t426;
t83 = rSges(6,3) * t303 + t365;
t81 = -rSges(6,3) * t392 + t398;
t66 = -t281 * t379 - t283 * t431;
t65 = t281 * t378 + t283 * t432;
t52 = qJD(1) * t114 + t282 * t316;
t51 = t284 * t316 + t369 + t426;
t50 = t488 * t389 + (t284 * t305 - t285) * qJD(1) - t365 + t374;
t49 = qJD(1) * t295 + t298 + t398;
t48 = t300 * t281;
t47 = -t332 + t375;
t44 = (t217 * t416 + t83) * t283 + (-qJD(3) * t154 + t175 * t282 + t217 * t419) * t281;
t43 = (-t217 * t414 - t81) * t283 + (qJD(3) * t156 - t175 * t284 + t393) * t281;
t42 = t375 - t501;
t41 = t489 * t389 + (t284 * t306 - t285) * qJD(1) + t374 + t539;
t40 = qJD(1) * t296 + t298 + t511;
t38 = -t126 * t209 + t127 * t213 + t163 * t442 - t167 * t227 + t171 * t228 + t205 * t303;
t37 = t126 * t204 + t127 * t212 + t162 * t227 + t166 * t442 + t170 * t228 + t208 * t303;
t36 = t124 * t209 + t125 * t213 + t163 * t441 - t167 * t229 + t171 * t230 + t205 * t302;
t35 = -t124 * t204 + t125 * t212 + t162 * t229 + t166 * t441 + t170 * t230 + t208 * t302;
t30 = t333 * t415 + (qJD(1) * t332 - t282 * t81 + t284 * t83) * t281;
t25 = (qJD(3) * t378 + t472) * t283 + (-qJD(3) * t432 + t282 * t430 + t284 * t377) * t281;
t24 = (-qJD(3) * t379 - t473) * t283 + (qJD(3) * t431 - t284 * t430 + t369) * t281;
t23 = t282 * t83 + t284 * t81 + (t451 + (-t156 + t428) * t282) * qJD(1) + t317;
t10 = t300 * t415 + (qJD(1) * t501 - t473 * t282 + t472 * t284) * t281;
t9 = t473 * t284 + t472 * t282 + (t380 + (t428 - t431) * t282) * qJD(1) + t317;
t4 = (qJD(3) * t358 - t38) * t283 + (-qJD(1) * t32 + qJD(3) * t92 + t17 * t282 + t18 * t284) * t281;
t3 = (qJD(3) * t359 - t37) * t283 + (-qJD(1) * t31 + qJD(3) * t91 + t15 * t282 + t16 * t284) * t281;
t2 = (qJD(3) * t356 - t36) * t283 + (-qJD(1) * t34 + qJD(3) * t94 + t13 * t282 + t14 * t284) * t281;
t1 = (qJD(3) * t357 - t35) * t283 + (-qJD(1) * t33 + qJD(3) * t93 + t11 * t282 + t12 * t284) * t281;
t5 = [(t115 * t179 + t116 * t178) * t500 + (t103 * t160 + t104 * t159) * t499 + (t101 * t50 + t102 * t49) * t498 + (t40 * t89 + t41 * t88) * t497 - t213 * t385 - t283 * t163 - t167 * t440 + t512 + (-t346 + t353) * t417 + (t347 + t352) * t415 + (-t348 + t355) * t412 + (t349 + t354) * t410 + t527 * t209; 0; 0; m(6) * (t101 * t86 + t102 * t87 + t136 * t49 + t137 * t50) + m(7) * (t113 * t40 + t114 * t41 + t51 * t88 + t52 * t89) + m(5) * (t103 * t200 + t104 * t201 + t121 * t159 + t122 * t160) + m(4) * ((-t115 * t282 - t116 * t284) * t270 + (-t178 * t284 - t179 * t282) * t253) + ((-t179 * t487 + t449 / 0.2e1 + t447 / 0.2e1 + t445 / 0.2e1 + t443 / 0.2e1 - t320) * t284 + (t450 / 0.2e1 + t448 / 0.2e1 + t178 * t487 + t446 / 0.2e1 + t444 / 0.2e1 + t321) * t282) * qJD(1) + t542 * qJD(3) * (t495 / 0.2e1 + t496 / 0.2e1) + (-qJD(3) * t504 + (-qJD(1) * t193 - t284 * t309) * t283 + (-qJD(1) * t195 - t284 * t311) * t281 + (-qJD(1) * t210 - t284 * t310) * t293 + (-qJD(1) * t214 - t284 * t312) * t290 + t35 + t36 + t520) * t494 - (-t533 * qJD(3) + (qJD(1) * t194 - t282 * t309) * t283 + (qJD(1) * t196 - t282 * t311) * t281 + (qJD(1) * t211 - t282 * t310) * t293 + (qJD(1) * t215 - t282 * t312) * t290 + t37 + t38 - t521) * t284 / 0.2e1; m(4) * t67 + m(5) * t39 + m(6) * t23 + m(7) * t9; (t113 * t52 + t114 * t51 + t42 * t9) * t497 + (t136 * t87 + t137 * t86 + t23 * t47) * t498 + (t201 * t121 + t200 * t122 + t39 * t429) * t499 + (t495 + t496) * t270 * t253 * t500 + (t32 + t31) * t421 + (t34 + t33) * t419 + (t198 * t517 + t219 * t516 + t515 * t421 + t507 * t495 - t522 + (-t284 * t533 - t525) * t419) * t284 + (t197 * t517 + t218 * t516 + t508 * t496 + t510 * t419 + ((t194 * t415 + t196 * t417 + t211 * t410 + t215 * t412 + t507 - t526) * t282 + (t193 * t415 + t195 * t417 + t210 * t410 + t214 * t412 + t508) * t284 + ((-t447 - t449 - t443 - t445) * t282 + (-t444 - t446 - t448 - t450) * t284) * qJD(3) + ((-t533 + t535) * t282 + t537 + t510 + t515) * qJD(1)) * t284 + t523 + (t525 + t538) * t421) * t282; m(7) * (t282 * t41 - t284 * t40 + (t282 * t89 + t284 * t88) * qJD(1)) + m(6) * (t282 * t50 - t284 * t49 + (t101 * t284 + t102 * t282) * qJD(1)) + m(5) * (-t103 * t284 + t104 * t282 + (t159 * t284 + t160 * t282) * qJD(1)); 0; m(7) * (t282 * t51 - t284 * t52 + (t113 * t282 + t114 * t284) * qJD(1)) + m(6) * (t282 * t86 - t284 * t87 + (t136 * t282 + t137 * t284) * qJD(1)) + m(5) * (t121 * t282 - t122 * t284 + (t200 * t282 + t201 * t284) * qJD(1)); 0; m(6) * (t101 * t44 + t102 * t43 + t105 * t50 + t106 * t49) + m(7) * (t24 * t89 + t25 * t88 + t40 * t66 + t41 * t65) + ((t282 * t321 - t284 * t320) * qJD(3) - t476) * t283 + ((t22 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1) * t284 + (t21 / 0.2e1 + t19 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1) * t282 + (t282 * t320 + t284 * t321) * qJD(1)) * t281 - t466; m(6) * t30 + m(7) * t10; m(6) * (t105 * t86 + t106 * t87 + t136 * t43 + t137 * t44 + t23 * t90 + t30 * t47) + m(7) * (t10 * t42 + t113 * t24 + t114 * t25 + t48 * t9 + t51 * t65 + t52 * t66) + (-t3 / 0.2e1 - t4 / 0.2e1 + t400 * t415) * t284 + (t2 / 0.2e1 + t1 / 0.2e1 + t401 * t415) * t282 + ((-t282 * t400 + t284 * t401) * qJD(1) + t522 * t494 + t523 * t284 / 0.2e1 + (t282 * t474 - t284 * t475) * qJD(3) / 0.2e1) * t281 - (qJD(1) * t524 + t520 * t282 + t521 * t284) * t283 / 0.2e1 + (t519 * t282 + t518 * t284) * qJD(1) / 0.2e1; m(6) * (t282 * t44 - t284 * t43 + (t105 * t284 + t106 * t282) * qJD(1)) + m(7) * (-t24 * t284 + t25 * t282 + (t282 * t66 + t284 * t65) * qJD(1)); (t10 * t48 + t24 * t66 + t25 * t65) * t497 + (t105 * t44 + t106 * t43 + t30 * t90) * t498 + (t476 * t283 + (t323 * t282 - t284 * t502) * qJD(3) + t466) * t283 + ((-t283 * t520 + t1 + t2) * t284 + (t283 * t521 + t3 + t4) * t282 + (t281 * t524 + t283 * t528) * qJD(3) + (t282 * t502 + t323 * t284) * qJD(1)) * t281; m(7) * (-t124 * t88 + t126 * t89 + t227 * t40 + t229 * t41); -t527 * m(7); m(7) * (t42 * t387 + t113 * t126 - t114 * t124 + t227 * t52 + t229 * t51 + (t289 * t9 + t407 * t42) * t281); m(7) * (-t124 * t282 - t126 * t284 + (t227 * t282 + t229 * t284) * qJD(1)); m(7) * (t48 * t387 - t124 * t65 + t126 * t66 + t227 * t24 + t229 * t25 + (t10 * t289 + t407 * t48) * t281); (-t124 * t229 + t126 * t227 - t440 * t527) * t497;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
