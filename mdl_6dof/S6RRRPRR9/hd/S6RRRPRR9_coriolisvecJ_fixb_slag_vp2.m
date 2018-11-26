% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:57:01
% EndTime: 2018-11-23 17:57:49
% DurationCPUTime: 48.44s
% Computational Cost: add. (46224->1046), mult. (143760->1486), div. (0->0), fcn. (120646->14), ass. (0->460)
t401 = cos(qJ(2));
t393 = cos(pkin(6));
t531 = pkin(1) * t393;
t384 = t401 * t531;
t376 = qJD(1) * t384;
t397 = sin(qJ(2));
t390 = sin(pkin(6));
t392 = cos(pkin(7));
t464 = pkin(10) * t392 + pkin(9);
t446 = t390 * t464;
t419 = t397 * t446;
t314 = -qJD(1) * t419 + t376;
t383 = t397 * t531;
t408 = -t401 * t446 - t383;
t315 = t408 * qJD(1);
t389 = sin(pkin(7));
t529 = pkin(10) * t389;
t411 = (pkin(2) * t397 - t401 * t529) * t390;
t343 = qJD(1) * t411;
t400 = cos(qJ(3));
t467 = qJD(3) * t400;
t458 = t392 * t467;
t396 = sin(qJ(3));
t489 = t392 * t396;
t494 = t389 * t396;
t630 = -pkin(2) * t458 + t400 * t314 + t315 * t489 + t343 * t494;
t488 = t392 * t400;
t493 = t389 * t400;
t217 = -t314 * t396 + t315 * t488 + t343 * t493;
t480 = t400 * t401;
t484 = t396 * t397;
t415 = -t392 * t484 + t480;
t471 = qJD(1) * t390;
t336 = t415 * t471;
t381 = pkin(2) * t489;
t463 = t397 * t471;
t449 = t389 * t463;
t524 = -pkin(10) - qJ(4);
t629 = -pkin(3) * t449 + qJ(4) * t336 - t217 - qJD(4) * t494 + (t493 * t524 - t381) * qJD(3);
t482 = t397 * t400;
t483 = t396 * t401;
t417 = -t392 * t482 - t483;
t335 = t417 * t471;
t454 = t524 * t396;
t628 = qJ(4) * t335 - (qJD(3) * t454 + qJD(4) * t400) * t389 + t630;
t388 = sin(pkin(13));
t391 = cos(pkin(13));
t265 = -t391 * t335 + t336 * t388;
t349 = (t388 * t400 + t391 * t396) * t389;
t341 = qJD(3) * t349;
t474 = t265 - t341;
t266 = t335 * t388 + t336 * t391;
t469 = qJD(3) * t389;
t490 = t391 * t400;
t342 = (-t388 * t396 + t490) * t469;
t473 = t266 - t342;
t605 = t388 * t629 - t628 * t391;
t258 = -t315 * t389 + t392 * t343;
t223 = -pkin(3) * t335 + t258;
t468 = qJD(3) * t396;
t460 = t389 * t468;
t627 = pkin(3) * t460 - t223;
t626 = pkin(11) * t449 - t605;
t625 = -pkin(4) * t474 + t473 * pkin(11) + t627;
t606 = t628 * t388 + t391 * t629;
t395 = sin(qJ(5));
t399 = cos(qJ(5));
t234 = t266 * t395 - t399 * t449;
t313 = t349 * t399 + t392 * t395;
t261 = qJD(5) * t313 + t342 * t395;
t624 = t234 - t261;
t235 = t266 * t399 + t395 * t449;
t312 = t349 * t395 - t399 * t392;
t260 = -qJD(5) * t312 + t342 * t399;
t475 = t235 - t260;
t378 = qJD(1) * t393 + qJD(2);
t418 = t392 * t480 - t484;
t410 = t418 * t390;
t281 = qJD(1) * t410 + t378 * t493;
t416 = t392 * t483 + t482;
t409 = t416 * t390;
t282 = qJD(1) * t409 + t378 * t494;
t451 = t391 * t281 - t282 * t388;
t222 = qJD(5) - t451;
t406 = (qJD(2) * t415 + qJD(3) * t418) * t390;
t459 = t389 * t467;
t242 = qJD(1) * t406 + t378 * t459;
t405 = (qJD(2) * t417 - qJD(3) * t416) * t390;
t243 = qJD(1) * t405 - t378 * t460;
t182 = t242 * t391 + t243 * t388;
t462 = t401 * t471;
t337 = t392 * t378 - t389 * t462 + qJD(3);
t422 = t281 * t388 + t391 * t282;
t207 = t337 * t399 - t395 * t422;
t470 = qJD(2) * t390;
t453 = qJD(1) * t470;
t447 = t397 * t453;
t421 = t389 * t447;
t121 = qJD(5) * t207 + t182 * t399 + t395 * t421;
t569 = t121 / 0.2e1;
t623 = Ifges(6,4) * t569;
t382 = pkin(2) * t488;
t322 = pkin(3) * t392 + t389 * t454 + t382;
t361 = pkin(10) * t493 + t381;
t339 = qJ(4) * t493 + t361;
t263 = t388 * t322 + t391 * t339;
t249 = pkin(11) * t392 + t263;
t348 = t388 * t494 - t389 * t490;
t369 = (-pkin(3) * t400 - pkin(2)) * t389;
t270 = pkin(4) * t348 - pkin(11) * t349 + t369;
t465 = qJD(5) * t399;
t466 = qJD(5) * t395;
t608 = -t249 * t466 + t270 * t465 + t625 * t395 - t399 * t626;
t607 = pkin(4) * t449 - t606;
t279 = pkin(2) * t378 + t314;
t340 = (-pkin(2) * t401 - t397 * t529 - pkin(1)) * t390;
t329 = qJD(1) * t340;
t236 = -t279 * t389 + t392 * t329;
t209 = -pkin(3) * t281 + qJD(4) + t236;
t117 = -pkin(4) * t451 - pkin(11) * t422 + t209;
t492 = t390 * t401;
t275 = t378 * t529 + (t464 * t492 + t383) * qJD(1);
t201 = -t275 * t396 + t279 * t488 + t329 * t493;
t166 = -qJ(4) * t282 + t201;
t157 = pkin(3) * t337 + t166;
t202 = t396 * (t279 * t392 + t329 * t389) + t275 * t400;
t167 = qJ(4) * t281 + t202;
t491 = t391 * t167;
t87 = t388 * t157 + t491;
t82 = pkin(11) * t337 + t87;
t44 = t117 * t395 + t399 * t82;
t34 = pkin(12) * t222 + t44;
t394 = sin(qJ(6));
t398 = cos(qJ(6));
t208 = t337 * t395 + t399 * t422;
t163 = t388 * t167;
t86 = t157 * t391 - t163;
t81 = -pkin(4) * t337 - t86;
t56 = -pkin(5) * t207 - pkin(12) * t208 + t81;
t16 = -t34 * t394 + t398 * t56;
t17 = t34 * t398 + t394 * t56;
t599 = t16 * mrSges(7,1) - t17 * mrSges(7,2);
t596 = t81 * mrSges(6,1) - t44 * mrSges(6,3) + t599;
t206 = qJD(6) - t207;
t505 = Ifges(7,3) * t206;
t144 = -t208 * t394 + t222 * t398;
t506 = Ifges(7,6) * t144;
t145 = t208 * t398 + t222 * t394;
t509 = Ifges(7,5) * t145;
t68 = t505 + t506 + t509;
t507 = Ifges(6,6) * t222;
t508 = Ifges(6,2) * t207;
t517 = Ifges(6,4) * t208;
t98 = t507 + t508 + t517;
t616 = t98 / 0.2e1 - t68 / 0.2e1;
t622 = t616 - t596;
t122 = qJD(5) * t208 + t182 * t395 - t399 * t421;
t568 = -t122 / 0.2e1;
t181 = t242 * t388 - t391 * t243;
t560 = t181 / 0.2e1;
t621 = -pkin(12) * t474 + t608;
t620 = -pkin(5) * t624 + pkin(12) * t475 + t607;
t93 = t166 * t388 + t491;
t619 = -t93 + t222 * (pkin(5) * t395 - pkin(12) * t399);
t386 = pkin(3) * t388 + pkin(11);
t530 = pkin(3) * t282;
t142 = pkin(4) * t422 - pkin(11) * t451 + t530;
t94 = t166 * t391 - t163;
t63 = t395 * t142 + t399 * t94;
t618 = pkin(12) * t422 + t386 * t466 + t63;
t587 = t399 * t249 + t395 * t270;
t609 = -qJD(5) * t587 + t395 * t626 + t625 * t399;
t43 = t117 * t399 - t395 * t82;
t33 = -pkin(5) * t222 - t43;
t434 = t16 * t398 + t17 * t394;
t436 = Ifges(7,5) * t398 - Ifges(7,6) * t394;
t512 = Ifges(7,4) * t398;
t438 = -Ifges(7,2) * t394 + t512;
t513 = Ifges(7,4) * t394;
t440 = Ifges(7,1) * t398 - t513;
t441 = mrSges(7,1) * t394 + mrSges(7,2) * t398;
t534 = t398 / 0.2e1;
t535 = -t394 / 0.2e1;
t557 = t206 / 0.2e1;
t563 = t145 / 0.2e1;
t565 = t144 / 0.2e1;
t514 = Ifges(7,4) * t145;
t69 = Ifges(7,2) * t144 + Ifges(7,6) * t206 + t514;
t143 = Ifges(7,4) * t144;
t70 = Ifges(7,1) * t145 + Ifges(7,5) * t206 + t143;
t617 = -t434 * mrSges(7,3) + t33 * t441 + t436 * t557 + t438 * t565 + t440 * t563 + t534 * t70 + t535 * t69;
t372 = qJD(2) * t376;
t412 = qJD(2) * t419;
t296 = -qJD(1) * t412 + t372;
t317 = t408 * qJD(2);
t297 = qJD(1) * t317;
t344 = qJD(2) * t411;
t338 = qJD(1) * t344;
t138 = -qJD(3) * t202 - t296 * t396 + t297 * t488 + t338 * t493;
t88 = pkin(3) * t421 - qJ(4) * t242 - qJD(4) * t282 + t138;
t137 = -t275 * t468 + t279 * t458 + t400 * t296 + t297 * t489 + t329 * t459 + t338 * t494;
t91 = qJ(4) * t243 + qJD(4) * t281 + t137;
t41 = -t388 * t91 + t391 * t88;
t36 = -pkin(4) * t421 - t41;
t567 = t122 / 0.2e1;
t48 = -qJD(6) * t145 - t121 * t394 + t181 * t398;
t578 = t48 / 0.2e1;
t47 = qJD(6) * t144 + t121 * t398 + t181 * t394;
t579 = t47 / 0.2e1;
t19 = pkin(5) * t122 - pkin(12) * t121 + t36;
t42 = t388 * t88 + t391 * t91;
t37 = pkin(11) * t421 + t42;
t247 = -t297 * t389 + t392 * t338;
t191 = -pkin(3) * t243 + t247;
t78 = pkin(4) * t181 - pkin(11) * t182 + t191;
t7 = t117 * t465 + t399 * t37 + t395 * t78 - t466 * t82;
t5 = pkin(12) * t181 + t7;
t1 = qJD(6) * t16 + t19 * t394 + t398 * t5;
t2 = -qJD(6) * t17 + t19 * t398 - t394 * t5;
t600 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t615 = -t36 * mrSges(6,1) - Ifges(7,5) * t579 + 0.2e1 * Ifges(6,2) * t568 + 0.2e1 * Ifges(6,6) * t560 - Ifges(7,6) * t578 - Ifges(7,3) * t567 - t600 + t623;
t591 = Ifges(4,3) + Ifges(5,3);
t184 = pkin(12) * t348 + t587;
t262 = t322 * t391 - t388 * t339;
t248 = -pkin(4) * t392 - t262;
t198 = pkin(5) * t312 - pkin(12) * t313 + t248;
t125 = t184 * t398 + t198 * t394;
t614 = -qJD(6) * t125 - t394 * t621 + t398 * t620;
t124 = -t184 * t394 + t198 * t398;
t613 = qJD(6) * t124 + t394 * t620 + t398 * t621;
t387 = -pkin(3) * t391 - pkin(4);
t364 = -pkin(5) * t399 - pkin(12) * t395 + t387;
t486 = t394 * t399;
t305 = t364 * t398 - t386 * t486;
t612 = qJD(6) * t305 + t394 * t619 - t398 * t618;
t481 = t398 * t399;
t306 = t364 * t394 + t386 * t481;
t611 = -qJD(6) * t306 + t394 * t618 + t398 * t619;
t610 = pkin(5) * t474 - t609;
t177 = Ifges(5,6) * t181;
t178 = Ifges(5,5) * t182;
t104 = Ifges(5,3) * t421 - t177 + t178;
t238 = Ifges(4,6) * t243;
t239 = Ifges(4,5) * t242;
t171 = Ifges(4,3) * t421 + t238 + t239;
t604 = t171 + t104;
t603 = -pkin(10) * t460 - t630;
t602 = -t361 * qJD(3) - t217;
t205 = Ifges(6,4) * t207;
t510 = Ifges(6,5) * t222;
t522 = Ifges(6,1) * t208;
t99 = t205 + t510 + t522;
t597 = -t81 * mrSges(6,2) + t43 * mrSges(6,3) - t99 / 0.2e1;
t601 = -t205 / 0.2e1 - t510 / 0.2e1 + t597 - t617;
t598 = Ifges(6,1) * t569 + Ifges(6,5) * t560;
t595 = mrSges(5,1) * t209 + mrSges(6,1) * t43 - t44 * mrSges(6,2) - t87 * mrSges(5,3);
t11 = Ifges(7,5) * t47 + Ifges(7,6) * t48 + Ifges(7,3) * t122;
t593 = -mrSges(6,3) * t7 - t623 + t11 / 0.2e1 - t615;
t561 = -t181 / 0.2e1;
t559 = t182 / 0.2e1;
t545 = t242 / 0.2e1;
t544 = t243 / 0.2e1;
t592 = t421 / 0.2e1;
t547 = -t422 / 0.2e1;
t546 = t422 / 0.2e1;
t549 = -t451 / 0.2e1;
t590 = Ifges(5,4) * t422;
t589 = Ifges(5,4) * t451;
t472 = pkin(9) * t492 + t383;
t299 = (t389 * t393 + t392 * t492) * pkin(10) + t472;
t311 = pkin(2) * t393 + t384 - t419;
t215 = -t299 * t396 + t311 * t488 + t340 * t493;
t303 = t393 * t494 + t409;
t358 = -t389 * t492 + t393 * t392;
t174 = pkin(3) * t358 - qJ(4) * t303 + t215;
t287 = t400 * t299;
t216 = t311 * t489 + t340 * t494 + t287;
t302 = t393 * t493 + t410;
t192 = qJ(4) * t302 + t216;
t111 = t388 * t174 + t391 * t192;
t102 = pkin(11) * t358 + t111;
t253 = -t311 * t389 + t392 * t340;
t221 = -pkin(3) * t302 + t253;
t240 = -t391 * t302 + t303 * t388;
t241 = t302 * t388 + t303 * t391;
t136 = pkin(4) * t240 - pkin(11) * t241 + t221;
t588 = t399 * t102 + t395 * t136;
t444 = t1 * t398 - t2 * t394;
t8 = -qJD(5) * t44 - t37 * t395 + t399 * t78;
t461 = t397 * t470;
t448 = t389 * t461;
t377 = qJD(2) * t384;
t316 = t377 - t412;
t150 = -t316 * t396 + t317 * t488 + t344 * t493 + (-t287 + (-t311 * t392 - t340 * t389) * t396) * qJD(3);
t256 = t393 * t459 + t406;
t107 = pkin(3) * t448 - qJ(4) * t256 - qJD(4) * t303 + t150;
t149 = -t299 * t468 + t311 * t458 + t400 * t316 + t317 * t489 + t340 * t459 + t344 * t494;
t257 = -t393 * t460 + t405;
t112 = qJ(4) * t257 + qJD(4) * t302 + t149;
t58 = t388 * t107 + t391 * t112;
t55 = pkin(11) * t448 + t58;
t196 = t256 * t388 - t391 * t257;
t197 = t256 * t391 + t257 * t388;
t259 = -t317 * t389 + t392 * t344;
t204 = -pkin(3) * t257 + t259;
t84 = pkin(4) * t196 - pkin(11) * t197 + t204;
t15 = -qJD(5) * t588 - t395 * t55 + t399 * t84;
t12 = Ifges(7,4) * t47 + Ifges(7,2) * t48 + Ifges(7,6) * t122;
t584 = t12 / 0.2e1;
t583 = Ifges(7,1) * t579 + Ifges(7,4) * t578 + Ifges(7,5) * t567;
t38 = Ifges(6,5) * t121 - Ifges(6,6) * t122 + Ifges(6,3) * t181;
t582 = t38 / 0.2e1;
t580 = Ifges(6,4) * t568 + t598;
t576 = -t69 / 0.2e1;
t573 = pkin(1) * mrSges(3,1);
t572 = pkin(1) * mrSges(3,2);
t571 = -t182 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t560 - Ifges(5,6) * t421 / 0.2e1;
t570 = Ifges(5,1) * t559 + Ifges(5,4) * t561 + Ifges(5,5) * t592;
t566 = -t144 / 0.2e1;
t564 = -t145 / 0.2e1;
t562 = Ifges(4,1) * t545 + Ifges(4,4) * t544 + Ifges(4,5) * t592;
t558 = -t206 / 0.2e1;
t556 = -t207 / 0.2e1;
t555 = t207 / 0.2e1;
t554 = -t208 / 0.2e1;
t553 = t208 / 0.2e1;
t520 = Ifges(4,4) * t282;
t213 = t281 * Ifges(4,2) + t337 * Ifges(4,6) + t520;
t552 = -t213 / 0.2e1;
t551 = -t222 / 0.2e1;
t550 = t222 / 0.2e1;
t548 = t451 / 0.2e1;
t543 = -t281 / 0.2e1;
t542 = t281 / 0.2e1;
t541 = -t282 / 0.2e1;
t540 = t282 / 0.2e1;
t539 = -t337 / 0.2e1;
t538 = t337 / 0.2e1;
t533 = -t399 / 0.2e1;
t532 = t400 / 0.2e1;
t526 = t399 * t7;
t18 = -mrSges(7,1) * t48 + mrSges(7,2) * t47;
t76 = mrSges(6,1) * t181 - mrSges(6,3) * t121;
t523 = t18 - t76;
t521 = Ifges(3,4) * t397;
t519 = Ifges(4,4) * t396;
t518 = Ifges(4,4) * t400;
t516 = Ifges(6,4) * t395;
t515 = Ifges(6,4) * t399;
t511 = Ifges(3,5) * t401;
t504 = t451 * Ifges(5,6);
t503 = t422 * Ifges(5,5);
t502 = t281 * Ifges(4,6);
t501 = t282 * Ifges(4,5);
t500 = t337 * Ifges(5,3);
t499 = t378 * Ifges(3,5);
t148 = mrSges(6,1) * t222 - mrSges(6,3) * t208;
t80 = -mrSges(7,1) * t144 + mrSges(7,2) * t145;
t498 = t148 - t80;
t497 = t451 * t395;
t495 = t386 * t395;
t487 = t394 * t395;
t485 = t395 * t398;
t139 = -mrSges(6,1) * t207 + mrSges(6,2) * t208;
t211 = mrSges(5,1) * t337 - mrSges(5,3) * t422;
t479 = -t139 + t211;
t268 = -t313 * t394 + t348 * t398;
t179 = qJD(6) * t268 + t260 * t398 + t341 * t394;
t187 = t235 * t398 + t265 * t394;
t478 = t179 - t187;
t269 = t313 * t398 + t348 * t394;
t180 = -qJD(6) * t269 - t260 * t394 + t341 * t398;
t186 = -t235 * t394 + t265 * t398;
t477 = t180 - t186;
t113 = t181 * mrSges(5,1) + t182 * mrSges(5,2);
t57 = t107 * t391 - t388 * t112;
t110 = t174 * t391 - t388 * t192;
t345 = -pkin(9) * t447 + t372;
t357 = t472 * qJD(2);
t346 = qJD(1) * t357;
t443 = -t346 * mrSges(3,1) - t345 * mrSges(3,2);
t442 = mrSges(7,1) * t398 - mrSges(7,2) * t394;
t439 = Ifges(7,1) * t394 + t512;
t437 = Ifges(7,2) * t398 + t513;
t435 = Ifges(7,5) * t394 + Ifges(7,6) * t398;
t433 = t16 * t394 - t17 * t398;
t27 = mrSges(7,1) * t122 - mrSges(7,3) * t47;
t28 = -mrSges(7,2) * t122 + mrSges(7,3) * t48;
t432 = -t394 * t27 + t398 * t28;
t52 = pkin(12) * t240 + t588;
t101 = -pkin(4) * t358 - t110;
t219 = t241 * t395 - t399 * t358;
t220 = t241 * t399 + t358 * t395;
t72 = pkin(5) * t219 - pkin(12) * t220 + t101;
t23 = t394 * t72 + t398 * t52;
t22 = -t394 * t52 + t398 * t72;
t95 = -mrSges(7,2) * t206 + mrSges(7,3) * t144;
t96 = mrSges(7,1) * t206 - mrSges(7,3) * t145;
t430 = -t394 * t95 - t398 * t96;
t429 = t395 * t44 + t399 * t43;
t62 = t142 * t399 - t395 * t94;
t60 = -t102 * t395 + t136 * t399;
t162 = t220 * t398 + t240 * t394;
t161 = -t220 * t394 + t240 * t398;
t199 = -t249 * t395 + t270 * t399;
t14 = -t102 * t466 + t136 * t465 + t395 * t84 + t399 * t55;
t54 = -pkin(4) * t448 - t57;
t407 = t138 * mrSges(4,1) + t41 * mrSges(5,1) - t137 * mrSges(4,2) - t42 * mrSges(5,2);
t404 = -t509 / 0.2e1 - t505 / 0.2e1 - t506 / 0.2e1 + t507 / 0.2e1 + t517 / 0.2e1 + t622;
t374 = Ifges(3,4) * t462;
t371 = t453 * t511;
t360 = -pkin(9) * t390 * t397 + t384;
t359 = -pkin(10) * t494 + t382;
t356 = -pkin(9) * t461 + t377;
t353 = t472 * qJD(1);
t352 = -pkin(9) * t463 + t376;
t351 = -t378 * mrSges(3,2) + mrSges(3,3) * t462;
t350 = mrSges(3,1) * t378 - mrSges(3,3) * t463;
t324 = Ifges(3,1) * t463 + t374 + t499;
t323 = Ifges(3,6) * t378 + (Ifges(3,2) * t401 + t521) * t471;
t278 = Ifges(4,4) * t281;
t245 = mrSges(4,1) * t337 - mrSges(4,3) * t282;
t244 = -mrSges(4,2) * t337 + mrSges(4,3) * t281;
t231 = -mrSges(4,1) * t281 + mrSges(4,2) * t282;
t225 = -mrSges(4,2) * t421 + mrSges(4,3) * t243;
t224 = mrSges(4,1) * t421 - mrSges(4,3) * t242;
t214 = t282 * Ifges(4,1) + t337 * Ifges(4,5) + t278;
t212 = t337 * Ifges(4,3) + t501 + t502;
t210 = -mrSges(5,2) * t337 + mrSges(5,3) * t451;
t185 = -mrSges(4,1) * t243 + mrSges(4,2) * t242;
t183 = -pkin(5) * t348 - t199;
t172 = t242 * Ifges(4,4) + t243 * Ifges(4,2) + Ifges(4,6) * t421;
t169 = mrSges(5,1) * t421 - mrSges(5,3) * t182;
t168 = -mrSges(5,2) * t421 - mrSges(5,3) * t181;
t160 = -mrSges(5,1) * t451 + mrSges(5,2) * t422;
t159 = t394 * t422 + t451 * t481;
t158 = t398 * t422 - t451 * t486;
t154 = Ifges(5,1) * t422 + t337 * Ifges(5,5) + t589;
t153 = Ifges(5,2) * t451 + t337 * Ifges(5,6) + t590;
t152 = t500 + t503 + t504;
t147 = -mrSges(6,2) * t222 + mrSges(6,3) * t207;
t140 = pkin(5) * t208 - pkin(12) * t207;
t135 = qJD(5) * t220 + t197 * t395 - t399 * t448;
t134 = -qJD(5) * t219 + t197 * t399 + t395 * t448;
t97 = t208 * Ifges(6,5) + t207 * Ifges(6,6) + t222 * Ifges(6,3);
t77 = -mrSges(6,2) * t181 - mrSges(6,3) * t122;
t65 = -qJD(6) * t162 - t134 * t394 + t196 * t398;
t64 = qJD(6) * t161 + t134 * t398 + t196 * t394;
t59 = mrSges(6,1) * t122 + mrSges(6,2) * t121;
t51 = -pkin(5) * t240 - t60;
t49 = -pkin(5) * t422 - t62;
t30 = t140 * t394 + t398 * t43;
t29 = t140 * t398 - t394 * t43;
t24 = pkin(5) * t135 - pkin(12) * t134 + t54;
t10 = -pkin(5) * t196 - t15;
t9 = pkin(12) * t196 + t14;
t6 = -pkin(5) * t181 - t8;
t4 = -qJD(6) * t23 + t24 * t398 - t394 * t9;
t3 = qJD(6) * t22 + t24 * t394 + t398 * t9;
t13 = [(t220 * t36 - t240 * t7) * mrSges(6,2) + t593 * t219 + (Ifges(6,5) * t220 + Ifges(6,3) * t240) * t560 + (t212 + t152) * t448 / 0.2e1 + (t401 * t324 + t378 * (-Ifges(3,6) * t397 + t511)) * t470 / 0.2e1 + m(6) * (t101 * t36 + t14 * t44 + t15 * t43 + t54 * t81 + t588 * t7 + t60 * t8) + t588 * t77 + m(3) * (t345 * t472 - t346 * t360 - t352 * t357 + t353 * t356) + (t97 / 0.2e1 - t153 / 0.2e1 - Ifges(5,6) * t538 - Ifges(5,4) * t546 - Ifges(5,2) * t548 + Ifges(6,3) * t550 + Ifges(6,5) * t553 + Ifges(6,6) * t555 + t595) * t196 + (Ifges(7,5) * t64 + Ifges(7,6) * t65) * t557 + (Ifges(7,5) * t162 + Ifges(7,6) * t161) * t567 + (Ifges(6,4) * t220 + Ifges(6,6) * t240) * t568 + (-Ifges(6,4) * t553 + Ifges(7,5) * t563 - Ifges(6,2) * t555 - Ifges(6,6) * t550 + Ifges(7,6) * t565 + Ifges(7,3) * t557 - t622) * t135 + (t345 * t401 + t346 * t397 + (-t352 * t401 - t353 * t397) * qJD(2)) * t390 * mrSges(3,3) + (Ifges(7,1) * t64 + Ifges(7,4) * t65) * t563 + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t579 + (t191 * t241 + t197 * t209 - t358 * t42 - t448 * t87) * mrSges(5,2) + (t1 * t161 - t16 * t64 - t162 * t2 + t17 * t65) * mrSges(7,3) + (Ifges(5,4) * t197 + Ifges(5,6) * t448) * t548 + (Ifges(5,1) * t197 + Ifges(5,5) * t448) * t546 + (Ifges(6,1) * t553 + Ifges(6,4) * t555 + Ifges(6,5) * t550 - t597) * t134 + (Ifges(7,4) * t64 + Ifges(7,2) * t65) * t565 + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t578 + (Ifges(6,1) * t220 + Ifges(6,5) * t240) * t569 + m(7) * (t1 * t23 + t10 * t33 + t16 * t4 + t17 * t3 + t2 * t22 + t51 * t6) + m(4) * (t137 * t216 + t138 * t215 + t149 * t202 + t150 * t201 + t236 * t259 + t247 * t253) + m(5) * (t110 * t41 + t111 * t42 + t191 * t221 + t204 * t209 + t57 * t86 + t58 * t87) + t604 * t358 / 0.2e1 - t323 * t461 / 0.2e1 + t138 * (mrSges(4,1) * t358 - mrSges(4,3) * t303) + t137 * (-mrSges(4,2) * t358 + mrSges(4,3) * t302) + t41 * (mrSges(5,1) * t358 - mrSges(5,3) * t241) + t86 * (mrSges(5,1) * t448 - mrSges(5,3) * t197) + t201 * (mrSges(4,1) * t448 - mrSges(4,3) * t256) + t356 * t351 - t357 * t350 + t202 * (-mrSges(4,2) * t448 + mrSges(4,3) * t257) + t247 * (-mrSges(4,1) * t302 + mrSges(4,2) * t303) + t302 * t172 / 0.2e1 + t256 * t214 / 0.2e1 + t236 * (-mrSges(4,1) * t257 + mrSges(4,2) * t256) + t257 * t213 / 0.2e1 + t259 * t231 + t253 * t185 + t149 * t244 + t150 * t245 + ((Ifges(3,5) * t393 / 0.2e1 - t360 * mrSges(3,3) + (-0.2e1 * t572 + 0.3e1 / 0.2e1 * Ifges(3,4) * t401) * t390) * t401 + (-Ifges(3,6) * t393 - t472 * mrSges(3,3) + (-0.2e1 * t573 - 0.3e1 / 0.2e1 * t521 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t401) * t390 + (Ifges(4,5) * t303 + Ifges(5,5) * t241 + Ifges(4,6) * t302 - Ifges(5,6) * t240 + t358 * t591) * t389 / 0.2e1) * t397) * t453 + (Ifges(4,5) * t256 + Ifges(5,5) * t197 + Ifges(4,6) * t257 + t448 * t591) * t538 + t22 * t27 + t23 * t28 + t51 * t18 + t33 * (-mrSges(7,1) * t65 + mrSges(7,2) * t64) + t65 * t69 / 0.2e1 + t64 * t70 / 0.2e1 + t60 * t76 + t10 * t80 + t3 * t95 + t4 * t96 + t101 * t59 + t54 * t139 - t42 * mrSges(5,3) * t240 + t191 * mrSges(5,1) * t240 + t14 * t147 + t15 * t148 + t6 * (-mrSges(7,1) * t161 + mrSges(7,2) * t162) + t111 * t168 + t110 * t169 + t197 * t154 / 0.2e1 + t204 * t160 + t58 * t210 + t57 * t211 + (Ifges(4,1) * t256 + Ifges(4,4) * t257 + Ifges(4,5) * t448) * t540 + (Ifges(4,4) * t256 + Ifges(4,2) * t257 + Ifges(4,6) * t448) * t542 + (Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * t358) * t544 + (Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * t358) * t545 + (Ifges(5,1) * t241 - Ifges(5,4) * t240 + Ifges(5,5) * t358) * t559 + (Ifges(5,4) * t241 - Ifges(5,2) * t240 + Ifges(5,6) * t358) * t561 + t303 * t562 + t241 * t570 + t240 * t571 + t220 * t580 + t240 * t582 + t162 * t583 + t161 * t584 + t221 * t113 + t215 * t224 + t216 * t225 + (t371 / 0.2e1 + t443) * t393 + t8 * (mrSges(6,1) * t240 - mrSges(6,3) * t220); t593 * t312 + (t313 * t36 - t348 * t7 + t44 * t474 - t475 * t81) * mrSges(6,2) + (t97 - t153) * (t341 / 0.2e1 - t265 / 0.2e1) + (t68 - t98) * (t261 / 0.2e1 - t234 / 0.2e1) + t587 * t77 + (Ifges(6,1) * t313 + Ifges(6,5) * t348) * t569 + t371 - t596 * t624 + (t342 / 0.2e1 - t266 / 0.2e1) * t154 + (t214 * t532 + t518 * t542 - t519 * t540 + (t236 * mrSges(4,2) - t201 * mrSges(4,3) + Ifges(4,1) * t540 + Ifges(4,5) * t538) * t400 + (t236 * mrSges(4,1) - t202 * mrSges(4,3) - Ifges(4,2) * t542 - Ifges(4,6) * t538 + pkin(3) * t160 + t552) * t396) * t469 + (Ifges(6,4) * t313 + Ifges(6,6) * t348) * t568 + (t260 / 0.2e1 - t235 / 0.2e1) * t99 + (t191 * t369 + t209 * t627 + t262 * t41 + t263 * t42 + t605 * t87 + t606 * t86) * m(5) + (t179 / 0.2e1 - t187 / 0.2e1) * t70 + t443 + (t180 / 0.2e1 - t186 / 0.2e1) * t69 + (t201 * t336 - t202 * t335) * mrSges(4,3) + ((Ifges(4,1) * t396 + t518) * t545 + (Ifges(4,2) * t400 + t519) * t544 + t247 * (-mrSges(4,1) * t400 + mrSges(4,2) * t396) + t396 * t562 - pkin(2) * t185 + t172 * t532 + (t137 * t400 - t138 * t396) * mrSges(4,3)) * t389 + (t239 / 0.2e1 + t104 / 0.2e1 + t171 / 0.2e1 - t177 / 0.2e1 + t178 / 0.2e1 + t238 / 0.2e1 + t407) * t392 + t602 * t245 + (-pkin(2) * t247 * t389 + t137 * t361 + t138 * t359 + t201 * t602 + t202 * t603 - t236 * t258) * m(4) + t603 * t244 + t605 * t210 + t606 * t211 + t607 * t139 + (-t348 * t42 - t349 * t41 + t473 * t86 + t474 * t87) * mrSges(5,3) + (-mrSges(5,1) * t474 - mrSges(5,2) * t473) * t209 + (-mrSges(6,1) * t474 + mrSges(6,3) * t475) * t43 + (-mrSges(7,1) * t477 + mrSges(7,2) * t478) * t33 + t369 * t113 + t359 * t224 + t361 * t225 - t352 * t351 + t353 * t350 + t191 * (mrSges(5,1) * t348 + mrSges(5,2) * t349) + t8 * (mrSges(6,1) * t348 - mrSges(6,3) * t313) + (Ifges(7,5) * t269 + Ifges(7,6) * t268) * t567 - t236 * (-mrSges(4,1) * t335 + mrSges(4,2) * t336) - t336 * t214 / 0.2e1 + (Ifges(7,1) * t269 + Ifges(7,4) * t268) * t579 + t6 * (-mrSges(7,1) * t268 + mrSges(7,2) * t269) + t262 * t169 + t263 * t168 - t258 * t231 + t248 * t59 + t608 * t147 + t609 * t148 + (t199 * t8 + t248 * t36 + t43 * t609 + t44 * t608 + t587 * t7 + t607 * t81) * m(6) + t610 * t80 + t613 * t95 + t614 * t96 + (t1 * t125 + t124 * t2 + t16 * t614 + t17 * t613 + t183 * t6 + t33 * t610) * m(7) + ((-t374 / 0.2e1 - t324 / 0.2e1 + t352 * mrSges(3,3) - t499 / 0.2e1 + t471 * t572) * t401 + (t323 / 0.2e1 + t353 * mrSges(3,3) + (t573 + t521 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t401) * t471 + (t378 / 0.2e1 - qJD(2)) * Ifges(3,6) + (-t500 / 0.2e1 + t202 * mrSges(4,2) - t502 / 0.2e1 - t501 / 0.2e1 - t201 * mrSges(4,1) - t503 / 0.2e1 - t86 * mrSges(5,1) + t87 * mrSges(5,2) - t504 / 0.2e1 - t212 / 0.2e1 - t152 / 0.2e1 + t539 * Ifges(4,3) + (Ifges(5,5) * t349 - Ifges(5,6) * t348 + (Ifges(4,5) * t396 + Ifges(4,6) * t400) * t389 + t591 * t392) * qJD(2) / 0.2e1) * t389) * t397) * t471 + (Ifges(6,5) * t313 + Ifges(6,3) * t348) * t560 + (Ifges(7,4) * t269 + Ifges(7,2) * t268) * t578 + t124 * t27 + t125 * t28 + (t1 * t268 - t16 * t478 + t17 * t477 - t2 * t269) * mrSges(7,3) + t183 * t18 + t199 * t76 + (Ifges(5,5) * t342 - Ifges(5,6) * t341) * t538 + (Ifges(5,5) * t266 - Ifges(5,6) * t265) * t539 + (Ifges(4,5) * t336 + Ifges(4,6) * t335) * t539 + (Ifges(4,1) * t336 + Ifges(4,4) * t335) * t541 + (Ifges(4,4) * t336 + Ifges(4,2) * t335) * t543 + (Ifges(5,1) * t342 - Ifges(5,4) * t341) * t546 + (Ifges(5,1) * t266 - Ifges(5,4) * t265) * t547 + (Ifges(5,4) * t342 - Ifges(5,2) * t341) * t548 + (Ifges(5,4) * t266 - Ifges(5,2) * t265) * t549 + (Ifges(6,5) * t260 - Ifges(6,6) * t261 + Ifges(6,3) * t341) * t550 + (Ifges(6,5) * t235 - Ifges(6,6) * t234 + Ifges(6,3) * t265) * t551 + t335 * t552 + (Ifges(6,1) * t260 - Ifges(6,4) * t261 + Ifges(6,5) * t341) * t553 + (Ifges(6,1) * t235 - Ifges(6,4) * t234 + Ifges(6,5) * t265) * t554 + (Ifges(6,4) * t260 - Ifges(6,2) * t261 + Ifges(6,6) * t341) * t555 + (Ifges(6,4) * t235 - Ifges(6,2) * t234 + Ifges(6,6) * t265) * t556 + (Ifges(7,5) * t179 + Ifges(7,6) * t180 + Ifges(7,3) * t261) * t557 + (Ifges(7,5) * t187 + Ifges(7,6) * t186 + Ifges(7,3) * t234) * t558 + (Ifges(5,1) * t349 - Ifges(5,4) * t348) * t559 + (Ifges(5,4) * t349 - Ifges(5,2) * t348) * t561 + (Ifges(7,1) * t179 + Ifges(7,4) * t180 + Ifges(7,5) * t261) * t563 + (Ifges(7,1) * t187 + Ifges(7,4) * t186 + Ifges(7,5) * t234) * t564 + (Ifges(7,4) * t179 + Ifges(7,2) * t180 + Ifges(7,6) * t261) * t565 + (Ifges(7,4) * t187 + Ifges(7,2) * t186 + Ifges(7,6) * t234) * t566 + t349 * t570 + t348 * t571 + t313 * t580 + t348 * t582 + t269 * t583 + t268 * t584 - t223 * t160; (t99 * t533 - t81 * (mrSges(6,1) * t395 + mrSges(6,2) * t399) + t86 * mrSges(5,3) + Ifges(5,5) * t539 + (Ifges(6,5) * t399 - Ifges(6,6) * t395) * t551 + (Ifges(6,1) * t399 - t516) * t554 + (-Ifges(6,2) * t395 + t515) * t556 - t209 * mrSges(5,2) + t429 * mrSges(6,3) + Ifges(5,1) * t547) * t451 + t516 * t568 + (Ifges(6,5) * t554 - Ifges(5,2) * t549 - Ifges(5,6) * t539 + Ifges(6,6) * t556 + Ifges(6,3) * t551 - t595) * t422 + (t36 * t387 - t43 * t62 - t44 * t63 - t495 * t8 - t81 * t93 + (-qJD(5) * t429 + t526) * t386) * m(6) + t604 + (-t1 * t487 - t158 * t17 + t159 * t16 - t2 * t485) * mrSges(7,3) + (-t590 + t97) * t547 + t526 * mrSges(6,3) + t515 * t569 + (-t160 * t282 + t168 * t388 + t169 * t391) * pkin(3) + (t589 + t154) * t549 + (Ifges(7,1) * t159 + Ifges(7,4) * t158) * t564 + (-Ifges(4,2) * t282 + t214 + t278) * t543 + (-t8 * mrSges(6,3) + (mrSges(7,3) * t433 + t33 * t442 + t398 * t576 + t435 * t558 + t437 * t566 + t439 * t564 + t535 * t70) * qJD(6) + t6 * t441 + t523 * t386 + (-t508 / 0.2e1 - t404 - t147 * t386) * qJD(5) + t36 * mrSges(6,2) + t436 * t567 + t438 * t578 + t440 * t579 + t580 + t598) * t395 + t407 + (-t209 * t530 + t86 * t93 - t87 * t94 + (t388 * t42 + t391 * t41) * pkin(3)) * m(5) - t12 * t487 / 0.2e1 + t479 * t93 + t387 * t59 + t305 * t27 + t306 * t28 - t236 * (mrSges(4,1) * t282 + mrSges(4,2) * t281) - t201 * t244 + t202 * t245 + (Ifges(7,4) * t159 + Ifges(7,2) * t158) * t566 + t611 * t96 + t612 * t95 + (t1 * t306 + t16 * t611 + t17 * t612 + t2 * t305 - t33 * t49 + t495 * t6) * m(7) + t153 * t546 + (t201 * t281 + t202 * t282) * mrSges(4,3) + (Ifges(7,5) * t159 + Ifges(7,6) * t158) * t558 + (Ifges(7,5) * t564 + Ifges(7,6) * t566 + Ifges(7,3) * t558 - t599 + t616) * t497 + ((t522 / 0.2e1 + (m(7) * t33 - t498) * t386 - t601) * qJD(5) + t77 * t386 + t615) * t399 - t49 * t80 - t63 * t147 - t62 * t148 - t159 * t70 / 0.2e1 - t33 * (-mrSges(7,1) * t158 + mrSges(7,2) * t159) - t94 * t210 + t11 * t533 + (Ifges(4,5) * t281 - Ifges(4,6) * t282) * t539 + t213 * t540 + (Ifges(4,1) * t281 - t520) * t541 + t158 * t576 + t485 * t583; -t158 * t96 - t159 * t95 - t451 * t210 + t479 * t422 + (-t451 * t147 + (-t394 * t96 + t398 * t95 + t147) * qJD(5) - t523) * t399 + (qJD(6) * t430 - t222 * t498 + t432 + t77) * t395 + t113 + (-t158 * t16 - t159 * t17 - t33 * t497 + (-qJD(5) * t433 - t6) * t399 + (qJD(5) * t33 - qJD(6) * t434 + t444) * t395) * m(7) + (-t422 * t81 + t395 * t7 + t399 * t8 + t222 * (-t395 * t43 + t399 * t44)) * m(6) + (t422 * t86 - t451 * t87 + t191) * m(5); t38 + t12 * t534 + t394 * t583 - t7 * mrSges(6,2) + t8 * mrSges(6,1) + t404 * t208 - pkin(5) * t18 - t30 * t95 - t29 * t96 - t43 * t147 + t498 * t44 + ((Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t208 + t601) * t207 + t617 * qJD(6) + t435 * t567 + t437 * t578 + t439 * t579 - t6 * t442 + t444 * mrSges(7,3) + (-pkin(5) * t6 - t16 * t29 - t17 * t30 - t33 * t44) * m(7) + (m(7) * t444 + t432 + (-m(7) * t434 + t430) * qJD(6)) * pkin(12); t17 * t96 - t16 * t95 + t69 * t563 + (Ifges(7,5) * t144 - Ifges(7,6) * t145) * t558 + (Ifges(7,1) * t144 - t514) * t564 - t33 * (mrSges(7,1) * t145 + mrSges(7,2) * t144) + (t144 * t16 + t145 * t17) * mrSges(7,3) + t11 + (-Ifges(7,2) * t145 + t143 + t70) * t566 + t600;];
tauc  = t13(:);
