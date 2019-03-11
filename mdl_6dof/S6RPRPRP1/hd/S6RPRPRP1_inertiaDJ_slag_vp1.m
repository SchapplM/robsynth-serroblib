% Calculate time derivative of joint inertia matrix for
% S6RPRPRP1
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:22
% EndTime: 2019-03-09 03:01:47
% DurationCPUTime: 16.26s
% Computational Cost: add. (30421->815), mult. (29311->1105), div. (0->0), fcn. (27403->10), ass. (0->402)
t556 = Icges(4,3) + Icges(5,3);
t284 = qJ(3) + pkin(10);
t279 = sin(t284);
t281 = cos(t284);
t289 = sin(qJ(3));
t292 = cos(qJ(3));
t555 = Icges(4,5) * t292 + Icges(5,5) * t281 - Icges(4,6) * t289 - Icges(5,6) * t279;
t288 = sin(qJ(5));
t291 = cos(qJ(5));
t340 = Icges(7,5) * t291 - Icges(7,6) * t288;
t201 = -Icges(7,3) * t281 + t279 * t340;
t341 = Icges(6,5) * t291 - Icges(6,6) * t288;
t202 = -Icges(6,3) * t281 + t279 * t341;
t554 = t201 + t202;
t285 = qJ(1) + pkin(9);
t280 = sin(t285);
t282 = cos(t285);
t549 = t280 * t556 + t555 * t282;
t468 = Icges(7,4) * t291;
t344 = -Icges(7,2) * t288 + t468;
t205 = -Icges(7,6) * t281 + t279 * t344;
t470 = Icges(6,4) * t291;
t345 = -Icges(6,2) * t288 + t470;
t206 = -Icges(6,6) * t281 + t279 * t345;
t548 = t206 + t205;
t469 = Icges(7,4) * t288;
t350 = Icges(7,1) * t291 - t469;
t209 = -Icges(7,5) * t281 + t279 * t350;
t471 = Icges(6,4) * t288;
t351 = Icges(6,1) * t291 - t471;
t210 = -Icges(6,5) * t281 + t279 * t351;
t540 = -t210 - t209;
t286 = -qJ(6) - pkin(8);
t553 = rSges(7,3) - t286;
t472 = Icges(5,4) * t281;
t347 = -Icges(5,2) * t279 + t472;
t191 = Icges(5,6) * t280 + t282 * t347;
t473 = Icges(5,4) * t279;
t353 = Icges(5,1) * t281 - t473;
t193 = Icges(5,5) * t280 + t282 * t353;
t474 = Icges(4,4) * t292;
t349 = -Icges(4,2) * t289 + t474;
t208 = Icges(4,6) * t280 + t282 * t349;
t475 = Icges(4,4) * t289;
t355 = Icges(4,1) * t292 - t475;
t212 = Icges(4,5) * t280 + t282 * t355;
t515 = t191 * t279 - t193 * t281 + t208 * t289 - t212 * t292;
t552 = t280 * t515;
t551 = t515 * t282;
t550 = t555 * t280 - t282 * t556;
t547 = (-Icges(4,5) * t289 - Icges(5,5) * t279 - Icges(4,6) * t292 - Icges(5,6) * t281) * qJD(3);
t190 = -Icges(5,6) * t282 + t280 * t347;
t192 = -Icges(5,5) * t282 + t280 * t353;
t207 = -Icges(4,6) * t282 + t280 * t349;
t211 = -Icges(4,5) * t282 + t280 * t355;
t546 = t190 * t279 - t192 * t281 + t207 * t289 - t211 * t292;
t442 = t282 * t291;
t447 = t280 * t288;
t223 = -t281 * t447 - t442;
t443 = t282 * t288;
t446 = t280 * t291;
t224 = t281 * t446 - t443;
t451 = t279 * t280;
t143 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t451;
t147 = Icges(6,4) * t224 + Icges(6,2) * t223 + Icges(6,6) * t451;
t151 = Icges(6,1) * t224 + Icges(6,4) * t223 + Icges(6,5) * t451;
t225 = -t281 * t443 + t446;
t226 = t281 * t442 + t447;
t450 = t279 * t282;
t59 = t143 * t450 + t147 * t225 + t151 * t226;
t144 = Icges(6,5) * t226 + Icges(6,6) * t225 + Icges(6,3) * t450;
t148 = Icges(6,4) * t226 + Icges(6,2) * t225 + Icges(6,6) * t450;
t152 = Icges(6,1) * t226 + Icges(6,4) * t225 + Icges(6,5) * t450;
t60 = t144 * t450 + t148 * t225 + t152 * t226;
t358 = t280 * t59 + t282 * t60;
t141 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t451;
t145 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t451;
t149 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t451;
t57 = t141 * t450 + t145 * t225 + t149 * t226;
t142 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t450;
t146 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t450;
t150 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t450;
t58 = t142 * t450 + t146 * t225 + t150 * t226;
t359 = t280 * t57 + t282 * t58;
t545 = t358 + t359;
t55 = t143 * t451 + t147 * t223 + t151 * t224;
t56 = t144 * t451 + t148 * t223 + t152 * t224;
t360 = t280 * t55 + t282 * t56;
t53 = t141 * t451 + t145 * t223 + t149 * t224;
t54 = t142 * t451 + t146 * t223 + t150 * t224;
t361 = t280 * t53 + t282 * t54;
t544 = t360 + t361;
t543 = t554 * t281 + (t288 * t548 + t291 * t540) * t279;
t417 = qJD(5) * t279;
t162 = (-Icges(7,5) * t288 - Icges(7,6) * t291) * t417 + (Icges(7,3) * t279 + t281 * t340) * qJD(3);
t163 = (-Icges(6,5) * t288 - Icges(6,6) * t291) * t417 + (Icges(6,3) * t279 + t281 * t341) * qJD(3);
t542 = -t163 - t162;
t538 = t549 * qJD(1);
t537 = t550 * t280 - t552 + (-t546 - t549) * t282;
t166 = (-Icges(7,2) * t291 - t469) * t417 + (Icges(7,6) * t279 + t281 * t344) * qJD(3);
t167 = (-Icges(6,2) * t291 - t471) * t417 + (Icges(6,6) * t279 + t281 * t345) * qJD(3);
t536 = (-t167 - t166) * t288;
t335 = -t146 * t288 + t150 * t291;
t64 = -t142 * t281 + t279 * t335;
t333 = -t148 * t288 + t152 * t291;
t66 = -t144 * t281 + t279 * t333;
t485 = t64 + t66;
t336 = -t145 * t288 + t149 * t291;
t63 = -t141 * t281 + t279 * t336;
t334 = -t147 * t288 + t151 * t291;
t65 = -t143 * t281 + t279 * t334;
t486 = t63 + t65;
t535 = t280 * t486 + t282 * t485;
t277 = t280 ^ 2;
t278 = t282 ^ 2;
t380 = -qJD(5) * t281 + qJD(1);
t421 = qJD(3) * t288;
t297 = t279 * t421 + t291 * t380;
t428 = qJD(1) * t281;
t379 = -qJD(5) + t428;
t123 = t282 * t297 + t379 * t447;
t323 = t380 * t288;
t419 = qJD(3) * t291;
t296 = -t279 * t419 + t323;
t124 = t282 * t296 - t379 * t446;
t422 = qJD(3) * t282;
t393 = t281 * t422;
t429 = qJD(1) * t280;
t396 = t279 * t429;
t301 = t393 - t396;
t125 = t280 * t297 - t379 * t443;
t126 = t280 * t296 + t379 * t442;
t423 = qJD(3) * t281;
t394 = t280 * t423;
t427 = qJD(1) * t282;
t302 = t279 * t427 + t394;
t71 = Icges(7,5) * t126 + Icges(7,6) * t125 + Icges(7,3) * t302;
t75 = Icges(7,4) * t126 + Icges(7,2) * t125 + Icges(7,6) * t302;
t79 = Icges(7,1) * t126 + Icges(7,4) * t125 + Icges(7,5) * t302;
t11 = t123 * t145 + t124 * t149 + t141 * t301 + t225 * t75 + t226 * t79 + t450 * t71;
t70 = Icges(7,5) * t124 + Icges(7,6) * t123 + Icges(7,3) * t301;
t74 = Icges(7,4) * t124 + Icges(7,2) * t123 + Icges(7,6) * t301;
t78 = Icges(7,1) * t124 + Icges(7,4) * t123 + Icges(7,5) * t301;
t12 = t123 * t146 + t124 * t150 + t142 * t301 + t225 * t74 + t226 * t78 + t450 * t70;
t73 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t302;
t77 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t302;
t81 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t302;
t13 = t123 * t147 + t124 * t151 + t143 * t301 + t225 * t77 + t226 * t81 + t450 * t73;
t72 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t301;
t76 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t301;
t80 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t301;
t14 = t123 * t148 + t124 * t152 + t144 * t301 + t225 * t76 + t226 * t80 + t450 * t72;
t534 = (-t11 - t13) * t282 + (t12 + t14) * t280 + t545 * qJD(1);
t15 = t125 * t145 + t126 * t149 + t141 * t302 + t223 * t75 + t224 * t79 + t451 * t71;
t16 = t125 * t146 + t126 * t150 + t142 * t302 + t223 * t74 + t224 * t78 + t451 * t70;
t17 = t125 * t147 + t126 * t151 + t143 * t302 + t223 * t77 + t224 * t81 + t451 * t73;
t18 = t125 * t148 + t126 * t152 + t144 * t302 + t223 * t76 + t224 * t80 + t451 * t72;
t533 = (-t15 - t17) * t282 + (t16 + t18) * t280 + t544 * qJD(1);
t19 = (qJD(3) * t336 - t71) * t281 + (qJD(3) * t141 - t288 * t75 + t291 * t79 + (-t145 * t291 - t149 * t288) * qJD(5)) * t279;
t21 = (qJD(3) * t334 - t73) * t281 + (qJD(3) * t143 - t288 * t77 + t291 * t81 + (-t147 * t291 - t151 * t288) * qJD(5)) * t279;
t532 = -t19 - t21;
t20 = (qJD(3) * t335 - t70) * t281 + (qJD(3) * t142 - t288 * t74 + t291 * t78 + (-t146 * t291 - t150 * t288) * qJD(5)) * t279;
t22 = (qJD(3) * t333 - t72) * t281 + (qJD(3) * t144 - t288 * t76 + t291 * t80 + (-t148 * t291 - t152 * t288) * qJD(5)) * t279;
t531 = t20 + t22;
t89 = t201 * t451 + t205 * t223 + t209 * t224;
t90 = t202 * t451 + t206 * t223 + t210 * t224;
t530 = (-t89 - t90) * t281 + t544 * t279;
t91 = t201 * t450 + t205 * t225 + t209 * t226;
t92 = t202 * t450 + t206 * t225 + t210 * t226;
t529 = (-t91 - t92) * t281 + t545 * t279;
t276 = pkin(3) * t292 + pkin(2);
t255 = t282 * t276;
t287 = -qJ(4) - pkin(7);
t490 = -pkin(7) - t287;
t185 = -pkin(2) * t282 + t280 * t490 + t255;
t481 = rSges(5,1) * t281;
t369 = -rSges(5,2) * t279 + t481;
t194 = -rSges(5,3) * t282 + t280 * t369;
t445 = t281 * t282;
t272 = t280 * rSges(5,3);
t517 = -rSges(5,2) * t450 + t272;
t195 = rSges(5,1) * t445 + t517;
t238 = rSges(5,1) * t279 + rSges(5,2) * t281;
t311 = qJD(3) * t238;
t274 = t282 * pkin(7);
t444 = t282 * t287;
t492 = pkin(2) - t276;
t522 = t280 * t492;
t184 = t274 + t444 - t522;
t270 = qJD(4) * t280;
t420 = qJD(3) * t289;
t409 = pkin(3) * t420;
t398 = qJD(4) * t282 + t280 * t409 + t287 * t429;
t495 = pkin(7) * t280;
t401 = t280 * ((-t282 * t492 - t495) * qJD(1) - t398) + t282 * (-t282 * t409 + t270 + (t282 * t490 + t522) * qJD(1)) + t184 * t427;
t434 = rSges(5,2) * t396 + rSges(5,3) * t427;
t39 = (qJD(1) * t194 - t282 * t311 + t434) * t282 + (-t280 * t311 + (-t185 - t195 + t517) * qJD(1)) * t280 + t401;
t508 = 2 * m(5);
t528 = t39 * t508;
t509 = 2 * m(4);
t480 = rSges(4,2) * t289;
t482 = rSges(4,1) * t292;
t370 = -t480 + t482;
t479 = rSges(4,3) * t282;
t215 = t280 * t370 - t479;
t411 = t282 * t480;
t273 = t280 * rSges(4,3);
t433 = t282 * t482 + t273;
t216 = -t411 + t433;
t265 = rSges(4,1) * t289 + rSges(4,2) * t292;
t312 = qJD(3) * t265;
t426 = qJD(1) * t289;
t395 = t280 * t426;
t295 = rSges(4,2) * t395 + rSges(4,3) * t427 - t282 * t312;
t67 = (qJD(1) * t215 + t295) * t282 + (-t280 * t312 + (-t216 - t411 + t273) * qJD(1)) * t280;
t527 = t509 * t67;
t526 = t280 * t546 + t282 * t550;
t523 = t280 * t549 - t551;
t489 = pkin(8) + t286;
t521 = t281 * t489;
t520 = -qJD(1) * t550 + t282 * t547;
t519 = -t280 * t547 - t538;
t267 = pkin(5) * t447;
t275 = pkin(5) * t291 + pkin(4);
t518 = t226 * rSges(7,1) + t225 * rSges(7,2) + t275 * t445 + t450 * t553 + t267;
t499 = sin(qJ(1)) * pkin(1);
t516 = t274 - t499;
t415 = qJD(5) * t291;
t406 = pkin(5) * t415;
t413 = pkin(5) * t443;
t414 = qJD(6) * t279;
t513 = t124 * rSges(7,1) + t123 * rSges(7,2) + rSges(7,3) * t393 + qJD(1) * t413 + t280 * t406 + t282 * t414 + t286 * t396;
t170 = (-Icges(7,1) * t288 - t468) * t417 + (Icges(7,5) * t279 + t281 * t350) * qJD(3);
t171 = (-Icges(6,1) * t288 - t470) * t417 + (Icges(6,5) * t279 + t281 * t351) * qJD(3);
t425 = qJD(3) * t279;
t512 = (t170 + t171) * t279 * t291 + t554 * t425 - t540 * t281 * t419;
t511 = t281 * t485 - t529;
t254 = pkin(4) * t445;
t222 = pkin(8) * t450 + t254;
t440 = -t222 + t518;
t491 = pkin(4) - t275;
t299 = -t279 * t489 - t281 * t491;
t364 = -rSges(7,1) * t224 - rSges(7,2) * t223;
t441 = rSges(7,3) * t451 + t280 * t299 - t364 - t413;
t510 = -t280 * t441 - t282 * t440;
t507 = 2 * m(6);
t506 = 2 * m(7);
t505 = t280 / 0.2e1;
t501 = -rSges(6,3) - pkin(8);
t500 = m(4) * t265;
t498 = pkin(3) * t289;
t497 = pkin(4) * t279;
t496 = pkin(4) * t281;
t283 = cos(qJ(1)) * pkin(1);
t487 = t512 + (-t421 * t548 + t542) * t281 + (t536 + (t288 * t540 - t291 * t548) * qJD(5)) * t279;
t244 = pkin(8) * t393;
t416 = qJD(5) * t288;
t407 = pkin(5) * t416;
t484 = -t244 + (pkin(8) * t429 + t422 * t491) * t279 + ((-qJD(3) * t286 - t407) * t282 + t491 * t429) * t281 - rSges(7,3) * t396 + t513;
t424 = qJD(3) * t280;
t243 = t424 * t497;
t365 = t126 * rSges(7,1) + t125 * rSges(7,2);
t452 = t275 * t279;
t483 = t243 + (qJD(1) * t299 - t406) * t282 + (t414 + pkin(5) * t323 + (-t452 - t521) * qJD(3)) * t280 + rSges(7,3) * t302 + t365;
t478 = rSges(7,3) * t279;
t476 = t543 * t425;
t367 = -rSges(6,1) * t224 - rSges(6,2) * t223;
t154 = rSges(6,3) * t451 - t367;
t461 = t154 * t282;
t460 = t190 * t281;
t459 = t191 * t281;
t458 = t192 * t279;
t457 = t193 * t279;
t456 = t207 * t292;
t455 = t208 * t292;
t454 = t211 * t289;
t453 = t212 * t289;
t363 = rSges(7,1) * t291 - rSges(7,2) * t288;
t390 = t279 * t416;
t439 = -pkin(5) * t390 - qJD(6) * t281 + (-rSges(7,1) * t288 - rSges(7,2) * t291) * t417 + (t281 * t363 + t299 + t478) * qJD(3);
t438 = t280 * t184 + t282 * t185;
t437 = -rSges(7,3) * t281 + t521 + (t363 - t491) * t279;
t436 = -t185 - t222;
t239 = -pkin(8) * t281 + t497;
t259 = pkin(3) * t395;
t435 = t239 * t429 + t259;
t432 = t277 + t278;
t418 = qJD(3) * t292;
t412 = m(7) * t425;
t408 = pkin(3) * t418;
t31 = t280 * t54 - t282 * t53;
t32 = t280 * t56 - t282 * t55;
t405 = t31 / 0.2e1 + t32 / 0.2e1;
t33 = t280 * t58 - t282 * t57;
t34 = t280 * t60 - t282 * t59;
t404 = -t34 / 0.2e1 - t33 / 0.2e1;
t402 = t124 * rSges(6,1) + t123 * rSges(6,2) + rSges(6,3) * t393;
t156 = t226 * rSges(6,1) + t225 * rSges(6,2) + rSges(6,3) * t450;
t366 = rSges(6,1) * t291 - rSges(6,2) * t288;
t214 = -rSges(6,3) * t281 + t279 * t366;
t397 = t214 * t429;
t388 = -t238 - t498;
t387 = -t239 - t498;
t386 = t280 * t437;
t385 = t282 * t437;
t384 = t441 * t282;
t383 = -t275 * t281 - t276;
t382 = qJD(1) * t437;
t381 = t281 * t407;
t374 = pkin(8) * t279 + t496;
t221 = t374 * t280;
t378 = t280 * t221 + t282 * t222 + t438;
t376 = -t214 + t387;
t375 = -t374 * qJD(3) - t408;
t373 = -t280 * t287 + t255 + t283;
t372 = t280 * t382;
t198 = t388 * t282;
t368 = t126 * rSges(6,1) + t125 * rSges(6,2);
t362 = -t444 - t499;
t61 = t279 * t386 + t281 * t441;
t62 = -t279 * t385 - t281 * t440;
t357 = t280 * t62 + t282 * t61;
t300 = -t279 * t553 + t383;
t93 = -t499 + (pkin(5) * t288 - t287) * t282 + t300 * t280 + t364;
t94 = t373 + t518;
t356 = t280 * t94 + t282 * t93;
t354 = Icges(4,1) * t289 + t474;
t352 = Icges(5,1) * t279 + t472;
t348 = Icges(4,2) * t292 + t475;
t346 = Icges(5,2) * t281 + t473;
t324 = t387 - t437;
t113 = t324 * t280;
t114 = t324 * t282;
t339 = t113 * t280 + t114 * t282;
t332 = -t156 * t280 + t461;
t331 = -t154 * t280 - t156 * t282;
t175 = (-rSges(6,1) * t288 - rSges(6,2) * t291) * t417 + (rSges(6,3) * t279 + t281 * t366) * qJD(3);
t322 = -t175 + t375;
t321 = -t281 * t486 + t530;
t320 = -pkin(2) - t370;
t318 = t65 / 0.2e1 + t63 / 0.2e1 + t89 / 0.2e1 + t90 / 0.2e1;
t317 = t66 / 0.2e1 + t64 / 0.2e1 + t91 / 0.2e1 + t92 / 0.2e1;
t138 = t376 * t282;
t316 = -t276 - t369;
t314 = t280 * (pkin(8) * t302 + qJD(1) * t254 - t243) + t282 * (-pkin(8) * t396 + t244 + (-t279 * t422 - t280 * t428) * pkin(4)) + t221 * t427 + t401;
t313 = t375 - t439;
t310 = qJD(3) * t354;
t309 = qJD(3) * t352;
t308 = qJD(3) * t348;
t307 = qJD(3) * t346;
t304 = t279 * t501 - t276 - t496;
t298 = -t280 * t440 + t384;
t294 = t280 * t304 + t362;
t248 = t370 * qJD(3);
t231 = t369 * qJD(3);
t197 = t388 * t280;
t177 = t495 + t283 + (pkin(2) - t480) * t282 + t433;
t176 = t280 * t320 + t479 + t516;
t161 = t195 + t373;
t160 = -t499 + (rSges(5,3) - t287) * t282 + t316 * t280;
t137 = t376 * t280;
t122 = -t238 * t427 - t231 * t280 + (-t280 * t418 - t282 * t426) * pkin(3);
t121 = t238 * t429 + t259 + (-t231 - t408) * t282;
t116 = t265 * t424 + (-t283 + (-rSges(4,3) - pkin(7)) * t280 + t320 * t282) * qJD(1);
t115 = ((-pkin(2) - t482) * t280 + t516) * qJD(1) + t295;
t106 = -t156 * t281 - t214 * t450;
t105 = t154 * t281 + t214 * t451;
t104 = t238 * t424 + (t282 * t316 - t272 - t283) * qJD(1) + t398;
t103 = t270 + qJD(3) * t198 + ((-t276 - t481) * t280 + t362) * qJD(1) + t434;
t102 = t373 + t156 + t222;
t101 = t294 + t367;
t88 = t332 * t279;
t87 = qJD(1) * t138 + t280 * t322;
t86 = t282 * t322 + t397 + t435;
t85 = rSges(6,3) * t302 + t368;
t83 = -rSges(6,3) * t396 + t402;
t52 = t243 + t501 * t394 + (t282 * t304 - t283) * qJD(1) - t368 + t398;
t51 = t244 + t270 + (-t497 - t498) * t422 + t294 * qJD(1) + t402;
t50 = qJD(1) * t114 + t280 * t313;
t49 = t282 * t313 + t372 + t435;
t48 = -t331 + t378;
t47 = t298 * t279;
t46 = t282 * t406 + (t381 - t414 + (-t281 * t553 + t452) * qJD(3)) * t280 + (t282 * t300 - t267 - t283) * qJD(1) - t365 + t398;
t45 = t270 + (-t381 + (-t281 * t286 - t452 - t498) * qJD(3)) * t282 + ((t383 - t478) * t280 + t362) * qJD(1) + t513;
t42 = (t214 * t424 + t85) * t281 + (-qJD(3) * t154 + t175 * t280 + t214 * t427) * t279;
t41 = (-t214 * t422 - t83) * t281 + (qJD(3) * t156 - t175 * t282 + t397) * t279;
t40 = t378 - t510;
t38 = t125 * t206 + t126 * t210 + t163 * t451 + t167 * t223 + t171 * t224 + t202 * t302;
t37 = t125 * t205 + t126 * t209 + t162 * t451 + t166 * t223 + t170 * t224 + t201 * t302;
t36 = t123 * t206 + t124 * t210 + t163 * t450 + t167 * t225 + t171 * t226 + t202 * t301;
t35 = t123 * t205 + t124 * t209 + t162 * t450 + t166 * t225 + t170 * t226 + t201 * t301;
t30 = t332 * t423 + (qJD(1) * t331 - t280 * t83 + t282 * t85) * t279;
t25 = (qJD(3) * t386 + t483) * t281 + (-qJD(3) * t441 + t280 * t439 + t282 * t382) * t279;
t24 = (-qJD(3) * t385 - t484) * t281 + (qJD(3) * t440 - t282 * t439 + t372) * t279;
t23 = t280 * t85 + t282 * t83 + (t461 + (-t156 + t436) * t280) * qJD(1) + t314;
t10 = t298 * t423 + (qJD(1) * t510 - t484 * t280 + t483 * t282) * t279;
t9 = t484 * t282 + t483 * t280 + (t384 + (t436 - t440) * t280) * qJD(1) + t314;
t4 = (qJD(3) * t360 - t38) * t281 + (-qJD(1) * t32 + qJD(3) * t90 + t17 * t280 + t18 * t282) * t279;
t3 = (qJD(3) * t361 - t37) * t281 + (-qJD(1) * t31 + qJD(3) * t89 + t15 * t280 + t16 * t282) * t279;
t2 = (qJD(3) * t358 - t36) * t281 + (-qJD(1) * t34 + qJD(3) * t92 + t13 * t280 + t14 * t282) * t279;
t1 = (qJD(3) * t359 - t35) * t281 + (-qJD(1) * t33 + qJD(3) * t91 + t11 * t280 + t12 * t282) * t279;
t5 = [t512 + (t45 * t94 + t46 * t93) * t506 + (t101 * t52 + t102 * t51) * t507 + (t103 * t161 + t104 * t160) * t508 + (t115 * t177 + t116 * t176) * t509 + t279 * t536 + (-t346 + t353) * t425 + (t347 + t352) * t423 + (-t348 + t355) * t420 + (t349 + t354) * t418 + t540 * t390 + t542 * t281 + t548 * (-t279 * t415 - t281 * t421); 0; 0; m(5) * (t103 * t197 + t104 * t198 + t121 * t160 + t122 * t161) + m(6) * (t101 * t86 + t102 * t87 + t137 * t51 + t138 * t52) + m(7) * (t113 * t45 + t114 * t46 + t49 * t93 + t50 * t94) + m(4) * ((-t115 * t280 - t116 * t282) * t265 + (-t176 * t282 - t177 * t280) * t248) + ((t459 / 0.2e1 + t457 / 0.2e1 - t177 * t500 + t455 / 0.2e1 + t453 / 0.2e1 + t317) * t282 + (t460 / 0.2e1 + t458 / 0.2e1 + t176 * t500 + t456 / 0.2e1 + t454 / 0.2e1 + t318) * t280) * qJD(1) + t555 * qJD(3) * (t277 / 0.2e1 + t278 / 0.2e1) + (-qJD(3) * t515 + (-qJD(1) * t190 - t282 * t307) * t281 + (-qJD(1) * t192 - t282 * t309) * t279 + (-qJD(1) * t207 - t282 * t308) * t292 + (-qJD(1) * t211 - t282 * t310) * t289 + t35 + t36 + t531) * t505 - (-qJD(3) * t546 + (qJD(1) * t191 - t280 * t307) * t281 + (qJD(1) * t193 - t280 * t309) * t279 + (qJD(1) * t208 - t280 * t308) * t292 + (qJD(1) * t212 - t280 * t310) * t289 + t37 + t38 - t532) * t282 / 0.2e1; m(4) * t67 + m(5) * t39 + m(6) * t23 + m(7) * t9; (t113 * t50 + t114 * t49 + t40 * t9) * t506 + (t137 * t87 + t138 * t86 + t23 * t48) * t507 + (t198 * t121 + t197 * t122 + t39 * t438) * t508 + t432 * t265 * t248 * t509 + (t32 + t31) * t429 + (t34 + t33) * t427 + (t195 * t528 + t216 * t527 + t519 * t278 + t526 * t429 - t533 + (-t282 * t546 - t537) * t427) * t282 + (t194 * t528 + t215 * t527 + t520 * t277 + t523 * t427 + ((t191 * t423 + t193 * t425 + t208 * t418 + t212 * t420 + t519 - t538) * t280 + (t190 * t423 + t192 * t425 + t207 * t418 + t211 * t420 + t520) * t282 + ((-t453 - t455 - t457 - t459) * t280 + (-t454 - t456 - t458 - t460) * t282) * qJD(3) + ((-t546 + t549) * t280 + t551 + t523 + t526) * qJD(1)) * t282 + t534 + (t537 + t552) * t429) * t280; m(7) * (qJD(1) * t356 + t280 * t46 - t282 * t45) + m(6) * (t280 * t52 - t282 * t51 + (t101 * t282 + t102 * t280) * qJD(1)) + m(5) * (-t103 * t282 + t104 * t280 + (t160 * t282 + t161 * t280) * qJD(1)); 0; m(7) * (qJD(1) * t339 + t280 * t49 - t282 * t50) + m(6) * (t280 * t86 - t282 * t87 + (t137 * t280 + t138 * t282) * qJD(1)) + m(5) * (t121 * t280 - t122 * t282 + (t197 * t280 + t198 * t282) * qJD(1)); 0; m(6) * (t101 * t42 + t102 * t41 + t105 * t52 + t106 * t51) + m(7) * (t24 * t94 + t25 * t93 + t45 * t62 + t46 * t61) + ((t280 * t318 + t282 * t317) * qJD(3) - t487) * t281 + ((t22 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1) * t282 + (t21 / 0.2e1 + t19 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1) * t280 + (-t280 * t317 + t282 * t318) * qJD(1)) * t279 - t476; m(6) * t30 + m(7) * t10; m(6) * (t105 * t86 + t106 * t87 + t137 * t41 + t138 * t42 + t23 * t88 + t30 * t48) + m(7) * (t10 * t40 + t113 * t24 + t114 * t25 + t47 * t9 + t49 * t61 + t50 * t62) + (-t4 / 0.2e1 - t3 / 0.2e1 - t404 * t423) * t282 + (t2 / 0.2e1 + t1 / 0.2e1 + t405 * t423) * t280 + ((t280 * t404 + t282 * t405) * qJD(1) + t533 * t505 + t534 * t282 / 0.2e1 + (t280 * t485 - t282 * t486) * qJD(3) / 0.2e1) * t279 - (qJD(1) * t535 + t531 * t280 + t532 * t282) * t281 / 0.2e1 + (t280 * t530 + t282 * t529) * qJD(1) / 0.2e1; m(6) * (t280 * t42 - t282 * t41 + (t105 * t282 + t106 * t280) * qJD(1)) + m(7) * (qJD(1) * t357 - t24 * t282 + t25 * t280); (t10 * t47 + t24 * t62 + t25 * t61) * t506 + (t105 * t42 + t106 * t41 + t30 * t88) * t507 + (t487 * t281 + (t321 * t280 - t282 * t511) * qJD(3) + t476) * t281 + ((-t281 * t531 + t1 + t2) * t282 + (t281 * t532 + t3 + t4) * t280 + (t535 * t279 + t281 * t543) * qJD(3) + (t280 * t511 + t321 * t282) * qJD(1)) * t279; m(7) * (t356 * t423 + (t280 * t45 + t282 * t46 + (-t280 * t93 + t282 * t94) * qJD(1)) * t279); t412; m(7) * ((qJD(3) * t339 - t9) * t281 + (qJD(3) * t40 + t280 * t50 + t282 * t49 + (t113 * t282 - t114 * t280) * qJD(1)) * t279); 0; m(7) * ((qJD(3) * t357 - t10) * t281 + (qJD(3) * t47 + t24 * t280 + t25 * t282 + (-t280 * t61 + t282 * t62) * qJD(1)) * t279); 0.2e1 * (-0.1e1 + t432) * t281 * t412;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
