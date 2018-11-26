% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:32:50
% EndTime: 2018-11-23 18:33:20
% DurationCPUTime: 30.88s
% Computational Cost: add. (24821->900), mult. (63618->1235), div. (0->0), fcn. (50037->10), ass. (0->365)
t363 = sin(qJ(4));
t365 = sin(qJ(2));
t367 = cos(qJ(4));
t360 = sin(pkin(6));
t423 = qJD(1) * t360;
t368 = cos(qJ(3));
t369 = cos(qJ(2));
t432 = t368 * t369;
t275 = (-t363 * t432 + t365 * t367) * t423;
t364 = sin(qJ(3));
t417 = qJD(4) * t367;
t419 = qJD(3) * t368;
t570 = t363 * t419 + t364 * t417 + t275;
t361 = cos(pkin(6));
t352 = t361 * t365 * pkin(1);
t435 = t360 * t369;
t424 = pkin(8) * t435 + t352;
t309 = t424 * qJD(1);
t348 = qJD(1) * t361 + qJD(2);
t267 = t348 * pkin(9) + t309;
t303 = (-pkin(2) * t369 - pkin(9) * t365 - pkin(1)) * t360;
t280 = qJD(1) * t303;
t201 = -t364 * t267 + t280 * t368;
t407 = t369 * t423;
t340 = qJD(3) - t407;
t180 = -pkin(3) * t340 - t201;
t408 = t365 * t423;
t288 = t348 * t364 + t368 * t408;
t233 = -t288 * t363 + t340 * t367;
t131 = -pkin(4) * t233 + t180;
t234 = t288 * t367 + t340 * t363;
t362 = sin(qJ(5));
t366 = cos(qJ(5));
t163 = t233 * t362 + t234 * t366;
t287 = t348 * t368 - t364 * t408;
t282 = qJD(4) - t287;
t470 = pkin(1) * t369;
t414 = t361 * t470;
t306 = -pkin(8) * t408 + qJD(1) * t414;
t266 = -t348 * pkin(2) - t306;
t176 = -t287 * pkin(3) - t288 * pkin(10) + t266;
t202 = t368 * t267 + t364 * t280;
t181 = pkin(10) * t340 + t202;
t105 = t367 * t176 - t181 * t363;
t90 = -pkin(11) * t234 + t105;
t76 = pkin(4) * t282 + t90;
t106 = t176 * t363 + t181 * t367;
t91 = pkin(11) * t233 + t106;
t87 = t362 * t91;
t33 = t366 * t76 - t87;
t552 = qJ(6) * t163;
t23 = t33 - t552;
t271 = qJD(5) + t282;
t19 = pkin(5) * t271 + t23;
t89 = t366 * t91;
t34 = t362 * t76 + t89;
t396 = t366 * t233 - t234 * t362;
t518 = qJ(6) * t396;
t24 = t34 + t518;
t486 = -t271 / 0.2e1;
t495 = t163 / 0.2e1;
t496 = -t163 / 0.2e1;
t421 = qJD(2) * t369;
t406 = t364 * t421;
t420 = qJD(3) * t364;
t243 = t348 * t420 + (t365 * t419 + t406) * t423;
t405 = t368 * t421;
t242 = t348 * t419 + (-t365 * t420 + t405) * t423;
t422 = qJD(2) * t360;
t397 = qJD(1) * t422;
t394 = t365 * t397;
t149 = qJD(4) * t233 + t242 * t367 + t363 * t394;
t150 = -qJD(4) * t234 - t242 * t363 + t367 * t394;
t65 = qJD(5) * t396 + t149 * t366 + t150 * t362;
t381 = (pkin(2) * t365 - pkin(9) * t369) * t360;
t308 = qJD(2) * t381;
t297 = qJD(1) * t308;
t436 = t360 * t365;
t349 = pkin(8) * t436;
t321 = -t349 + t414;
t310 = t321 * qJD(2);
t298 = qJD(1) * t310;
t132 = -t267 * t420 + t280 * t419 + t364 * t297 + t368 * t298;
t121 = pkin(10) * t394 + t132;
t311 = t424 * qJD(2);
t299 = qJD(1) * t311;
t146 = t243 * pkin(3) - t242 * pkin(10) + t299;
t44 = -qJD(4) * t106 - t121 * t363 + t367 * t146;
t22 = pkin(4) * t243 - pkin(11) * t149 + t44;
t418 = qJD(4) * t363;
t43 = t367 * t121 + t363 * t146 + t176 * t417 - t181 * t418;
t29 = pkin(11) * t150 + t43;
t8 = -qJD(5) * t34 + t366 * t22 - t29 * t362;
t2 = pkin(5) * t243 - qJ(6) * t65 - qJD(6) * t163 + t8;
t66 = -qJD(5) * t163 - t149 * t362 + t150 * t366;
t415 = qJD(5) * t366;
t416 = qJD(5) * t362;
t7 = t362 * t22 + t366 * t29 + t76 * t415 - t416 * t91;
t3 = qJ(6) * t66 + qJD(6) * t396 + t7;
t513 = t8 * mrSges(6,1) + t2 * mrSges(7,1) - t7 * mrSges(6,2) - t3 * mrSges(7,2);
t540 = Ifges(6,5) + Ifges(7,5);
t542 = Ifges(6,1) + Ifges(7,1);
t541 = Ifges(6,4) + Ifges(7,4);
t562 = t541 * t396;
t532 = t163 * t542 + t540 * t271 + t562;
t538 = Ifges(6,6) + Ifges(7,6);
t539 = Ifges(6,2) + Ifges(7,2);
t559 = t163 * t541;
t533 = t271 * t538 + t396 * t539 + t559;
t11 = Ifges(7,5) * t65 + Ifges(7,6) * t66 + Ifges(7,3) * t243;
t12 = Ifges(6,5) * t65 + Ifges(6,6) * t66 + Ifges(6,3) * t243;
t536 = t11 + t12;
t544 = -t396 / 0.2e1;
t86 = -pkin(5) * t396 + qJD(6) + t131;
t569 = t513 + t536 + (-t163 * t538 + t396 * t540) * t486 + (t163 * t24 + t19 * t396) * mrSges(7,3) + (t163 * t34 + t33 * t396) * mrSges(6,3) - t131 * (mrSges(6,1) * t163 + mrSges(6,2) * t396) - t86 * (mrSges(7,1) * t163 + mrSges(7,2) * t396) + t533 * t495 + (-t163 * t539 + t532 + t562) * t544 + (t542 * t396 - t559) * t496;
t307 = qJD(1) * t381;
t227 = t368 * t306 + t364 * t307;
t213 = pkin(10) * t408 + t227;
t392 = pkin(3) * t364 - pkin(10) * t368;
t230 = (t352 + (pkin(8) + t392) * t435) * qJD(1);
t143 = t367 * t213 + t363 * t230;
t334 = t392 * qJD(3);
t339 = -pkin(3) * t368 - pkin(10) * t364 - pkin(2);
t219 = t363 * t334 + t339 * t417 + (-t367 * t420 - t368 * t418) * pkin(9);
t568 = t219 - t143;
t142 = -t363 * t213 + t367 * t230;
t276 = (t363 * t365 + t367 * t432) * t423;
t433 = t367 * t368;
t354 = pkin(9) * t433;
t395 = t364 * t407;
t469 = pkin(9) * t363;
t425 = t367 * t334 + t420 * t469;
t567 = -pkin(4) * t395 + t276 * pkin(11) - t142 + (pkin(4) * t364 - pkin(11) * t433) * qJD(3) + (-t354 + (pkin(11) * t364 - t339) * t363) * qJD(4) + t425;
t566 = pkin(11) * t570 - t568;
t221 = pkin(3) * t288 - pkin(10) * t287;
t129 = -t201 * t363 + t367 * t221;
t506 = -pkin(11) - pkin(10);
t409 = qJD(4) * t506;
t467 = pkin(11) * t367;
t565 = -pkin(4) * t288 + t287 * t467 + t367 * t409 - t129;
t130 = t367 * t201 + t363 * t221;
t437 = t287 * t363;
t564 = -pkin(11) * t437 - t363 * t409 + t130;
t328 = t362 * t367 + t363 * t366;
t514 = qJD(4) + qJD(5);
t253 = t514 * t328;
t382 = t362 * t363 - t366 * t367;
t191 = -t253 * t364 - t382 * t419;
t207 = t275 * t362 + t276 * t366;
t431 = t191 - t207;
t313 = t382 * t364;
t192 = t313 * t514 - t328 * t419;
t206 = t275 * t366 - t276 * t362;
t430 = t192 - t206;
t535 = t243 * t538 + t539 * t66 + t541 * t65;
t563 = -t535 / 0.2e1;
t326 = t367 * t339;
t248 = -t364 * t467 + t326 + (-pkin(4) - t469) * t368;
t296 = t363 * t339 + t354;
t434 = t363 * t364;
t262 = -pkin(11) * t434 + t296;
t183 = t362 * t248 + t366 * t262;
t554 = -qJD(5) * t183 + t362 * t566 + t366 * t567;
t553 = t248 * t415 - t262 * t416 + t362 * t567 - t366 * t566;
t226 = -t364 * t306 + t307 * t368;
t212 = -pkin(3) * t408 - t226;
t550 = pkin(4) * t570 + pkin(9) * t419 - t212;
t204 = t328 * t287;
t561 = t204 - t253;
t205 = t382 * t287;
t252 = t514 * t382;
t428 = -t205 + t252;
t463 = Ifges(5,4) * t234;
t135 = Ifges(5,2) * t233 + Ifges(5,6) * t282 + t463;
t229 = Ifges(5,4) * t233;
t136 = t234 * Ifges(5,1) + t282 * Ifges(5,5) + t229;
t383 = t105 * t367 + t106 * t363;
t385 = Ifges(5,5) * t367 - Ifges(5,6) * t363;
t461 = Ifges(5,4) * t367;
t387 = -Ifges(5,2) * t363 + t461;
t462 = Ifges(5,4) * t363;
t389 = Ifges(5,1) * t367 - t462;
t390 = mrSges(5,1) * t363 + mrSges(5,2) * t367;
t472 = t367 / 0.2e1;
t473 = -t363 / 0.2e1;
t482 = t282 / 0.2e1;
t489 = t234 / 0.2e1;
t491 = t233 / 0.2e1;
t560 = t383 * mrSges(5,3) - t135 * t473 - t136 * t472 - t180 * t390 - t385 * t482 - t387 * t491 - t389 * t489;
t534 = t540 * t243 + t541 * t66 + t542 * t65;
t557 = t534 / 0.2e1;
t556 = qJD(6) * t313 + t554 - t431 * qJ(6) + (-t395 + t420) * pkin(5);
t312 = t328 * t364;
t555 = qJ(6) * t430 - qJD(6) * t312 + t553;
t341 = t506 * t363;
t342 = t506 * t367;
t521 = t341 * t415 + t342 * t416 + t362 * t565 - t366 * t564;
t270 = t362 * t341 - t366 * t342;
t520 = -qJD(5) * t270 + t362 * t564 + t366 * t565;
t551 = -pkin(5) * t430 + t550;
t281 = Ifges(4,4) * t287;
t444 = t340 * Ifges(4,5);
t447 = t288 * Ifges(4,1);
t198 = t281 + t444 + t447;
t376 = t201 * mrSges(4,3) - t198 / 0.2e1 - t266 * mrSges(4,2) - t444 / 0.2e1;
t549 = t376 + t560;
t543 = -Ifges(3,6) * t348 / 0.2e1;
t537 = Ifges(6,3) + Ifges(7,3);
t523 = qJ(6) * t561 - qJD(6) * t382 + t521;
t522 = -pkin(5) * t288 + qJ(6) * t428 - qJD(6) * t328 + t520;
t95 = -mrSges(6,1) * t396 + mrSges(6,2) * t163;
t519 = m(6) * t131 + t95;
t301 = t349 + (-pkin(2) - t470) * t361;
t314 = -t361 * t368 + t364 * t436;
t315 = t361 * t364 + t368 * t436;
t209 = t314 * pkin(3) - t315 * pkin(10) + t301;
t302 = pkin(9) * t361 + t424;
t223 = t368 * t302 + t364 * t303;
t211 = -pkin(10) * t435 + t223;
t124 = t363 * t209 + t367 * t211;
t165 = pkin(4) * t437 + t202;
t517 = pkin(4) * t418 - pkin(5) * t561 - t165;
t516 = -t363 * t44 + t367 * t43;
t515 = -t44 * mrSges(5,1) + t43 * mrSges(5,2);
t450 = t282 * Ifges(5,3);
t454 = t234 * Ifges(5,5);
t455 = t233 * Ifges(5,6);
t134 = t450 + t454 + t455;
t443 = t340 * Ifges(4,6);
t446 = t288 * Ifges(4,4);
t449 = t287 * Ifges(4,2);
t197 = t443 + t446 + t449;
t411 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t412 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t413 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t80 = t163 * Ifges(7,5) + Ifges(7,6) * t396 + t271 * Ifges(7,3);
t81 = t163 * Ifges(6,5) + Ifges(6,6) * t396 + t271 * Ifges(6,3);
t370 = t413 * t163 - t411 * t271 + t412 * t396 + t106 * mrSges(5,2) + t202 * mrSges(4,3) + t24 * mrSges(7,2) + t34 * mrSges(6,2) - t134 / 0.2e1 + t197 / 0.2e1 - t80 / 0.2e1 - t81 / 0.2e1 - t105 * mrSges(5,1) - t19 * mrSges(7,1) - t455 / 0.2e1 - t454 / 0.2e1 - t266 * mrSges(4,1) - t450 / 0.2e1 + t446 / 0.2e1 - t33 * mrSges(6,1) + t443 / 0.2e1;
t512 = t449 / 0.2e1 + t370;
t511 = t65 / 0.2e1;
t510 = t66 / 0.2e1;
t71 = t149 * Ifges(5,1) + t150 * Ifges(5,4) + t243 * Ifges(5,5);
t509 = t71 / 0.2e1;
t505 = pkin(1) * mrSges(3,1);
t504 = pkin(1) * mrSges(3,2);
t503 = -t135 / 0.2e1;
t502 = t149 / 0.2e1;
t501 = t150 / 0.2e1;
t498 = t396 / 0.2e1;
t492 = -t233 / 0.2e1;
t490 = -t234 / 0.2e1;
t488 = t243 / 0.2e1;
t485 = t271 / 0.2e1;
t484 = -t281 / 0.2e1;
t483 = -t282 / 0.2e1;
t479 = -t314 / 0.2e1;
t477 = t315 / 0.2e1;
t474 = t361 / 0.2e1;
t468 = pkin(9) * t368;
t39 = t366 * t90 - t87;
t464 = Ifges(3,4) * t365;
t148 = Ifges(5,5) * t149;
t147 = Ifges(5,6) * t150;
t457 = t132 * mrSges(4,2);
t133 = -t267 * t419 - t280 * t420 + t297 * t368 - t364 * t298;
t456 = t133 * mrSges(4,1);
t453 = t242 * Ifges(4,1);
t452 = t242 * Ifges(4,4);
t451 = t243 * Ifges(4,4);
t123 = t367 * t209 - t211 * t363;
t380 = -t367 * t315 + t363 * t435;
t101 = pkin(4) * t314 + pkin(11) * t380 + t123;
t258 = -t363 * t315 - t367 * t435;
t107 = pkin(11) * t258 + t124;
t53 = t362 * t101 + t366 * t107;
t427 = -mrSges(3,1) * t348 - mrSges(4,1) * t287 + mrSges(4,2) * t288 + mrSges(3,3) * t408;
t168 = -mrSges(5,1) * t233 + mrSges(5,2) * t234;
t245 = mrSges(4,1) * t340 - mrSges(4,3) * t288;
t426 = t245 - t168;
t335 = pkin(4) * t434 + t364 * pkin(9);
t69 = Ifges(5,3) * t243 + t147 + t148;
t410 = Ifges(4,5) * t242 - Ifges(4,6) * t243 + Ifges(4,3) * t394;
t357 = -pkin(4) * t367 - pkin(3);
t402 = t365 * t422;
t17 = -t66 * mrSges(7,1) + t65 * mrSges(7,2);
t38 = -t362 * t90 - t89;
t52 = t366 * t101 - t107 * t362;
t182 = t366 * t248 - t262 * t362;
t222 = -t364 * t302 + t303 * t368;
t269 = t366 * t341 + t342 * t362;
t210 = pkin(3) * t435 - t222;
t391 = mrSges(5,1) * t367 - mrSges(5,2) * t363;
t388 = Ifges(5,1) * t363 + t461;
t386 = Ifges(5,2) * t367 + t462;
t384 = Ifges(5,5) * t363 + Ifges(5,6) * t367;
t187 = t258 * t366 + t362 * t380;
t188 = t258 * t362 - t366 * t380;
t152 = -t302 * t419 - t303 * t420 + t308 * t368 - t364 * t310;
t256 = -qJD(3) * t314 + t360 * t405;
t174 = qJD(4) * t258 + t367 * t256 + t363 * t402;
t257 = qJD(3) * t315 + t360 * t406;
t151 = -t302 * t420 + t303 * t419 + t364 * t308 + t368 * t310;
t140 = pkin(10) * t402 + t151;
t169 = t257 * pkin(3) - t256 * pkin(10) + t311;
t57 = -qJD(4) * t124 - t140 * t363 + t367 * t169;
t32 = pkin(4) * t257 - pkin(11) * t174 + t57;
t175 = qJD(4) * t380 - t363 * t256 + t367 * t402;
t56 = t367 * t140 + t363 * t169 + t209 * t417 - t211 * t418;
t42 = pkin(11) * t175 + t56;
t9 = t101 * t415 - t107 * t416 + t362 * t32 + t366 * t42;
t164 = -pkin(4) * t258 + t210;
t343 = Ifges(3,4) * t407;
t377 = -t306 * mrSges(3,3) + Ifges(3,1) * t408 / 0.2e1 + t343 / 0.2e1 + t348 * Ifges(3,5);
t141 = -pkin(3) * t402 - t152;
t10 = -qJD(5) * t53 + t366 * t32 - t362 * t42;
t122 = -pkin(3) * t394 - t133;
t93 = -pkin(4) * t175 + t141;
t77 = -pkin(4) * t150 + t122;
t374 = t201 * mrSges(4,1) + t340 * Ifges(4,3) + t288 * Ifges(4,5) + t287 * Ifges(4,6) + t543 - (Ifges(3,2) * t369 + t464) * t423 / 0.2e1 - t202 * mrSges(4,2) - t309 * mrSges(3,3);
t356 = pkin(4) * t366 + pkin(5);
t338 = Ifges(3,5) * t369 * t397;
t305 = -t348 * mrSges(3,2) + mrSges(3,3) * t407;
t300 = pkin(5) * t382 + t357;
t295 = -t363 * t468 + t326;
t250 = pkin(5) * t312 + t335;
t244 = -mrSges(4,2) * t340 + mrSges(4,3) * t287;
t232 = -qJ(6) * t382 + t270;
t231 = -qJ(6) * t328 + t269;
t220 = -qJD(4) * t296 + t425;
t217 = -mrSges(4,2) * t394 - mrSges(4,3) * t243;
t216 = mrSges(4,1) * t394 - mrSges(4,3) * t242;
t186 = mrSges(5,1) * t282 - mrSges(5,3) * t234;
t185 = -mrSges(5,2) * t282 + mrSges(5,3) * t233;
t172 = mrSges(4,1) * t243 + mrSges(4,2) * t242;
t156 = Ifges(4,5) * t394 - t451 + t453;
t155 = -t243 * Ifges(4,2) + Ifges(4,6) * t394 + t452;
t154 = -qJ(6) * t312 + t183;
t153 = -pkin(5) * t368 + qJ(6) * t313 + t182;
t128 = mrSges(6,1) * t271 - mrSges(6,3) * t163;
t127 = mrSges(7,1) * t271 - mrSges(7,3) * t163;
t126 = -mrSges(6,2) * t271 + mrSges(6,3) * t396;
t125 = -mrSges(7,2) * t271 + mrSges(7,3) * t396;
t118 = pkin(4) * t234 + pkin(5) * t163;
t116 = -mrSges(5,2) * t243 + mrSges(5,3) * t150;
t115 = mrSges(5,1) * t243 - mrSges(5,3) * t149;
t103 = -pkin(5) * t187 + t164;
t94 = -mrSges(7,1) * t396 + mrSges(7,2) * t163;
t92 = -mrSges(5,1) * t150 + mrSges(5,2) * t149;
t73 = -qJD(5) * t188 - t174 * t362 + t175 * t366;
t72 = qJD(5) * t187 + t174 * t366 + t175 * t362;
t70 = t149 * Ifges(5,4) + t150 * Ifges(5,2) + t243 * Ifges(5,6);
t51 = -mrSges(6,2) * t243 + mrSges(6,3) * t66;
t50 = -mrSges(7,2) * t243 + mrSges(7,3) * t66;
t49 = mrSges(6,1) * t243 - mrSges(6,3) * t65;
t48 = mrSges(7,1) * t243 - mrSges(7,3) * t65;
t40 = qJ(6) * t187 + t53;
t36 = -pkin(5) * t73 + t93;
t35 = pkin(5) * t314 - qJ(6) * t188 + t52;
t27 = t39 - t552;
t26 = t38 - t518;
t25 = -pkin(5) * t66 + t77;
t18 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t5 = qJ(6) * t73 + qJD(6) * t187 + t9;
t4 = pkin(5) * t257 - qJ(6) * t72 - qJD(6) * t188 + t10;
t1 = [(-Ifges(5,5) * t380 + Ifges(5,6) * t258 + t540 * t188 + t538 * t187 + (Ifges(5,3) + t537) * t314) * t488 + t122 * (-mrSges(5,1) * t258 - mrSges(5,2) * t380) + t44 * (mrSges(5,1) * t314 + mrSges(5,3) * t380) + (-Ifges(5,4) * t380 + Ifges(5,2) * t258 + Ifges(5,6) * t314) * t501 + (-Ifges(5,1) * t380 + Ifges(5,4) * t258 + Ifges(5,5) * t314) * t502 - t380 * t509 + m(4) * (t132 * t223 + t133 * t222 + t151 * t202 + t152 * t201 + t266 * t311 + t299 * t301) + m(5) * (t105 * t57 + t106 * t56 + t122 * t210 + t123 * t44 + t124 * t43 + t141 * t180) + m(6) * (t10 * t33 + t131 * t93 + t164 * t77 + t34 * t9 + t52 * t8 + t53 * t7) + m(7) * (t103 * t25 + t19 * t4 + t2 * t35 + t24 * t5 + t3 * t40 + t36 * t86) + (t187 * t539 + t188 * t541 + t314 * t538) * t510 + (t257 * t538 + t539 * t73 + t541 * t72) * t498 + (t187 * t541 + t188 * t542 + t314 * t540) * t511 + (t257 * t540 + t541 * t73 + t542 * t72) * t495 + (t69 + t536) * t314 / 0.2e1 + (t257 * t537 + t538 * t73 + t540 * t72) * t485 + t532 * t72 / 0.2e1 + t533 * t73 / 0.2e1 + t535 * t187 / 0.2e1 + (-t132 * t314 - t133 * t315 - t201 * t256 - t202 * t257) * mrSges(4,3) + t188 * t557 + (t377 * t369 + (t543 + t374) * t365) * t422 + t288 * (Ifges(4,1) * t256 - Ifges(4,4) * t257) / 0.2e1 + t340 * (Ifges(4,5) * t256 - Ifges(4,6) * t257) / 0.2e1 + m(3) * (t298 * t424 - t299 * t321 - t306 * t311 + t309 * t310) + ((Ifges(3,5) * t474 - t321 * mrSges(3,3) + (-0.2e1 * t504 + 0.3e1 / 0.2e1 * Ifges(3,4) * t369) * t360) * t369 + (-Ifges(3,6) * t361 + Ifges(4,5) * t477 + Ifges(4,6) * t479 - t424 * mrSges(3,3) + (-0.2e1 * t505 - 0.3e1 / 0.2e1 * t464 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1) * t369) * t360) * t365) * t397 + t287 * (Ifges(4,4) * t256 - Ifges(4,2) * t257) / 0.2e1 + t266 * (mrSges(4,1) * t257 + mrSges(4,2) * t256) + t258 * t70 / 0.2e1 - t257 * t197 / 0.2e1 + t106 * (-mrSges(5,2) * t257 + mrSges(5,3) * t175) + t105 * (mrSges(5,1) * t257 - mrSges(5,3) * t174) + t19 * (mrSges(7,1) * t257 - mrSges(7,3) * t72) + t33 * (mrSges(6,1) * t257 - mrSges(6,3) * t72) + t24 * (-mrSges(7,2) * t257 + mrSges(7,3) * t73) + t34 * (-mrSges(6,2) * t257 + mrSges(6,3) * t73) + t256 * t198 / 0.2e1 + t151 * t244 + t152 * t245 + t222 * t216 + t223 * t217 + t210 * t92 + t77 * (-mrSges(6,1) * t187 + mrSges(6,2) * t188) - t435 * t456 + t25 * (-mrSges(7,1) * t187 + mrSges(7,2) * t188) + t56 * t185 + t57 * t186 + t175 * t135 / 0.2e1 + t180 * (-mrSges(5,1) * t175 + mrSges(5,2) * t174) + t174 * t136 / 0.2e1 + t164 * t18 + t141 * t168 - t410 * t435 / 0.2e1 - t243 * (Ifges(4,4) * t315 - Ifges(4,2) * t314 - Ifges(4,6) * t435) / 0.2e1 + t242 * (Ifges(4,1) * t315 - Ifges(4,4) * t314 - Ifges(4,5) * t435) / 0.2e1 + t298 * (-t361 * mrSges(3,2) + mrSges(3,3) * t435) + (-mrSges(3,1) * t361 + mrSges(4,1) * t314 + mrSges(4,2) * t315 + mrSges(3,3) * t436) * t299 + t131 * (-mrSges(6,1) * t73 + mrSges(6,2) * t72) + t9 * t126 + t4 * t127 + t427 * t311 + t10 * t128 + t5 * t125 + t123 * t115 + t124 * t116 + t103 * t17 + t36 * t94 + t93 * t95 + t86 * (-mrSges(7,1) * t73 + mrSges(7,2) * t72) + t301 * t172 + t310 * t305 + t43 * (-mrSges(5,2) * t314 + mrSges(5,3) * t258) + t2 * (mrSges(7,1) * t314 - mrSges(7,3) * t188) + t8 * (mrSges(6,1) * t314 - mrSges(6,3) * t188) + t3 * (-mrSges(7,2) * t314 + mrSges(7,3) * t187) + t7 * (-mrSges(6,2) * t314 + mrSges(6,3) * t187) + t435 * t457 + t338 * t474 + t156 * t477 + t155 * t479 + (Ifges(5,5) * t174 + Ifges(5,6) * t175 + Ifges(5,3) * t257) * t482 + (Ifges(5,1) * t174 + Ifges(5,4) * t175 + Ifges(5,5) * t257) * t489 + (Ifges(5,4) * t174 + Ifges(5,2) * t175 + Ifges(5,6) * t257) * t491 + t35 * t48 + t40 * t50 + t52 * t49 + t53 * t51 + (t134 + t81 + t80) * t257 / 0.2e1; (t220 - t142) * t186 - t180 * (-mrSges(5,1) * t275 + mrSges(5,2) * t276) - t276 * t136 / 0.2e1 + ((qJD(2) * (Ifges(4,5) * t364 + Ifges(4,6) * t368) / 0.2e1 + (t505 + t464 / 0.2e1) * t423 + (t348 / 0.2e1 - qJD(2)) * Ifges(3,6) - t374) * t365 + (-t343 / 0.2e1 + (t504 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t365) * t423 + (t484 - t447 / 0.2e1 + t376) * t368 + t512 * t364 - t377) * t369) * t423 + (-t299 * mrSges(4,1) - t147 / 0.2e1 - t148 / 0.2e1 + t452 / 0.2e1 + pkin(9) * t217 + t132 * mrSges(4,3) + t155 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1 - t69 / 0.2e1 + t412 * t66 + t413 * t65 + (-Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1 - t411) * t243 - t513 + t515) * t368 + t312 * t563 - m(4) * (t201 * t226 + t202 * t227 + t266 * t309) - m(5) * (t105 * t142 + t106 * t143 + t180 * t212) + t338 + (-t312 * t7 + t313 * t8 - t33 * t431 + t34 * t430) * mrSges(6,3) + (-t19 * t431 + t2 * t313 + t24 * t430 - t3 * t312) * mrSges(7,3) + t25 * (mrSges(7,1) * t312 - mrSges(7,2) * t313) + t77 * (mrSges(6,1) * t312 - mrSges(6,2) * t313) + m(5) * (t105 * t220 + t106 * t219 + t295 * t44 + t296 * t43) - t227 * t244 - t226 * t245 + t250 * t17 + m(4) * (-pkin(2) * t299 + t132 * t468) - t212 * t168 + t182 * t49 + t183 * t51 - pkin(2) * t172 + t154 * t50 + t153 * t48 + (-mrSges(7,1) * t430 + mrSges(7,2) * t431) * t86 + (-mrSges(6,1) * t430 + mrSges(6,2) * t431) * t131 + (t156 / 0.2e1 + t122 * t390 + t453 / 0.2e1 - t451 / 0.2e1 + t299 * mrSges(4,2) + t387 * t501 + t389 * t502 + t70 * t473 - t133 * mrSges(4,3) + t71 * t472 + (-t363 * t43 - t367 * t44) * mrSges(5,3) + (-m(4) * t133 + m(5) * t122 - t216 + t92) * pkin(9) + (t136 * t473 + t367 * t503 + t180 * t391 + t386 * t492 + t388 * t490 + t384 * t483 + (t105 * t363 - t106 * t367) * mrSges(5,3)) * qJD(4)) * t364 - t427 * t309 + ((t447 / 0.2e1 + t281 / 0.2e1 + (-m(4) * t201 + m(5) * t180 - t426) * pkin(9) - t549) * t368 + ((-m(4) * t202 - t244) * pkin(9) - t512) * t364) * qJD(3) + t295 * t115 + t296 * t116 + t550 * t95 - t298 * mrSges(3,2) + t551 * t94 - t299 * mrSges(3,1) - t306 * t305 + t568 * t185 + t335 * t18 + t532 * (-t207 / 0.2e1 + t191 / 0.2e1) + t533 * (-t206 / 0.2e1 + t192 / 0.2e1) + t553 * t126 + t554 * t128 + (t131 * t550 + t182 * t8 + t183 * t7 + t33 * t554 + t335 * t77 + t34 * t553) * m(6) + t555 * t125 + t556 * t127 + (t153 * t2 + t154 * t3 + t19 * t556 + t24 * t555 + t25 * t250 + t551 * t86) * m(7) - t534 * t313 / 0.2e1 + (t191 * t540 + t192 * t538) * t485 + (t206 * t538 + t207 * t540) * t486 + (-t312 * t538 - t313 * t540 + t364 * t385) * t488 + (t191 * t541 + t192 * t539) * t498 + (t206 * t539 + t207 * t541) * t544 + (-t312 * t539 - t313 * t541) * t510 + (t206 * t541 + t207 * t542) * t496 + (t191 * t542 + t192 * t541) * t495 + (-t312 * t541 - t313 * t542) * t511 + (t105 * t276 - t106 * t275) * mrSges(5,3) + (Ifges(5,5) * t276 + Ifges(5,6) * t275) * t483 + (Ifges(5,1) * t276 + Ifges(5,4) * t275) * t490 + (Ifges(5,4) * t276 + Ifges(5,2) * t275) * t492 + t275 * t503; (t328 * t541 - t382 * t539) * t510 + (t328 * t542 - t382 * t541) * t511 + (t328 * t540 - t382 * t538 + t384) * t488 + t25 * (mrSges(7,1) * t382 + mrSges(7,2) * t328) + t77 * (mrSges(6,1) * t382 + mrSges(6,2) * t328) + (-t204 * t539 - t205 * t541) * t544 + t532 * (-t252 / 0.2e1 + t205 / 0.2e1) + (-t252 * t541 - t253 * t539) * t498 + (pkin(4) * t363 * t519 - t560) * qJD(4) + (-t252 * t542 - t253 * t541) * t495 + (-t252 * t540 - t253 * t538) * t485 + t533 * (-t253 / 0.2e1 + t204 / 0.2e1) + t516 * mrSges(5,3) + ((-m(5) * t383 - t363 * t185 - t367 * t186) * qJD(4) + m(5) * t516 - t115 * t363 + t116 * t367) * pkin(10) + t517 * t94 + (-t204 * t541 - t205 * t542) * t496 + (-t204 * t538 - t205 * t540) * t486 + t328 * t557 + t382 * t563 + t456 - t457 + t370 * t288 + t520 * t128 + (-pkin(3) * t122 - t105 * t129 - t106 * t130 - t180 * t202) * m(5) + (-t131 * t165 + t269 * t8 + t270 * t7 + t33 * t520 + t34 * t521 + t357 * t77) * m(6) + t521 * t126 + t522 * t127 + t523 * t125 + (t19 * t522 + t2 * t231 + t232 * t3 + t24 * t523 + t25 * t300 + t517 * t86) * m(7) + (-t328 * t8 + t33 * t428 + t34 * t561 - t382 * t7) * mrSges(6,3) + (t19 * t428 - t2 * t328 + t24 * t561 - t3 * t382) * mrSges(7,3) + (-mrSges(7,1) * t561 - mrSges(7,2) * t428) * t86 + (-mrSges(6,1) * t561 - mrSges(6,2) * t428) * t131 + t269 * t49 + t270 * t51 - t201 * t244 + t231 * t48 + t232 * t50 - t130 * t185 - t129 * t186 - t165 * t95 + t426 * t202 - pkin(3) * t92 + (t484 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t288 + t549) * t287 + t300 * t17 + t410 + t357 * t18 - t122 * t391 + t70 * t472 + t386 * t501 + t388 * t502 + t363 * t509; -t515 + t69 + (t366 * t49 + (t50 + t51) * t362 + ((t125 + t126) * t366 + (-t127 - t128) * t362) * qJD(5) + m(6) * (-t33 * t416 + t34 * t415 + t362 * t7 + t366 * t8) - t519 * t234) * pkin(4) + ((-t19 * t416 + t24 * t415 + t3 * t362) * pkin(4) + t2 * t356 - t118 * t86 - t19 * t26 - t24 * t27) * m(7) - m(6) * (t33 * t38 + t34 * t39) + (t105 * t233 + t106 * t234) * mrSges(5,3) + (-Ifges(5,2) * t234 + t136 + t229) * t492 - t180 * (mrSges(5,1) * t234 + mrSges(5,2) * t233) - t105 * t185 + t106 * t186 - t26 * t127 - t38 * t128 - t39 * t126 - t27 * t125 - t118 * t94 + t356 * t48 + (Ifges(5,5) * t233 - Ifges(5,6) * t234) * t483 + t135 * t489 + (Ifges(5,1) * t233 - t463) * t490 + t569; (-t163 * t94 + t48) * pkin(5) + t24 * t127 + t34 * t128 - t33 * t126 - t23 * t125 + (-(-t19 + t23) * t24 + (-t163 * t86 + t2) * pkin(5)) * m(7) + t569; -t396 * t125 + t163 * t127 + 0.2e1 * (t25 / 0.2e1 + t24 * t544 + t19 * t495) * m(7) + t17;];
tauc  = t1(:);
