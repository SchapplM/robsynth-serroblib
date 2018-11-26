% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:16:03
% EndTime: 2018-11-23 18:16:37
% DurationCPUTime: 34.82s
% Computational Cost: add. (34296->956), mult. (89253->1339), div. (0->0), fcn. (71862->12), ass. (0->426)
t369 = sin(qJ(4));
t370 = sin(qJ(3));
t373 = cos(qJ(4));
t374 = cos(qJ(3));
t333 = t369 * t374 + t370 * t373;
t556 = qJD(3) + qJD(4);
t270 = t556 * t333;
t375 = cos(qJ(2));
t366 = sin(pkin(6));
t442 = qJD(1) * t366;
t420 = t375 * t442;
t280 = t333 * t420;
t588 = t270 - t280;
t332 = -t369 * t370 + t373 * t374;
t269 = t556 * t332;
t365 = sin(pkin(12));
t455 = cos(pkin(12));
t213 = t269 * t365 + t270 * t455;
t281 = t332 * t420;
t222 = t280 * t455 + t281 * t365;
t447 = t213 - t222;
t214 = t269 * t455 - t270 * t365;
t223 = -t280 * t365 + t281 * t455;
t446 = t214 - t223;
t371 = sin(qJ(2));
t367 = cos(pkin(6));
t441 = qJD(1) * t367;
t429 = pkin(1) * t441;
t318 = pkin(8) * t420 + t371 * t429;
t407 = t370 * t420;
t271 = pkin(3) * t407 + t318;
t438 = qJD(3) * t370;
t587 = pkin(3) * t438 - t271;
t421 = t371 * t442;
t314 = -pkin(8) * t421 + t375 * t429;
t389 = (pkin(2) * t371 - pkin(9) * t375) * t366;
t315 = qJD(1) * t389;
t251 = -t314 * t370 + t315 * t374;
t221 = (-pkin(10) * t374 * t375 + pkin(3) * t371) * t442 + t251;
t252 = t314 * t374 + t315 * t370;
t234 = -pkin(10) * t407 + t252;
t161 = t221 * t373 - t234 * t369;
t133 = pkin(4) * t421 - qJ(5) * t281 + t161;
t162 = t221 * t369 + t234 * t373;
t139 = -qJ(5) * t280 + t162;
t544 = -pkin(10) - pkin(9);
t423 = qJD(3) * t544;
t337 = t370 * t423;
t338 = t374 * t423;
t347 = t544 * t370;
t348 = t544 * t374;
t435 = qJD(4) * t373;
t436 = qJD(4) * t369;
t224 = t337 * t373 + t338 * t369 + t347 * t435 + t348 * t436;
t182 = -qJ(5) * t270 + qJD(5) * t332 + t224;
t283 = t347 * t369 - t348 * t373;
t225 = -qJD(4) * t283 - t337 * t369 + t338 * t373;
t378 = -qJ(5) * t269 - qJD(5) * t333 + t225;
t581 = (-t139 + t182) * t455 + (-t133 + t378) * t365;
t578 = pkin(4) * t588 + t587;
t586 = -pkin(11) * t421 + t581;
t585 = pkin(5) * t447 - pkin(11) * t446 + t578;
t584 = Ifges(5,3) + Ifges(6,3);
t266 = -t332 * t455 + t333 * t365;
t267 = t332 * t365 + t333 * t455;
t363 = -pkin(3) * t374 - pkin(2);
t305 = -pkin(4) * t332 + t363;
t192 = pkin(5) * t266 - pkin(11) * t267 + t305;
t254 = qJ(5) * t332 + t283;
t282 = t347 * t373 + t348 * t369;
t386 = -qJ(5) * t333 + t282;
t194 = t254 * t455 + t365 * t386;
t368 = sin(qJ(6));
t372 = cos(qJ(6));
t124 = t192 * t368 + t194 * t372;
t583 = -qJD(6) * t124 - t368 * t586 + t372 * t585;
t123 = t192 * t372 - t194 * t368;
t582 = qJD(6) * t123 + t368 * t585 + t372 * t586;
t580 = -t161 + t225;
t579 = -t162 + t224;
t355 = qJD(2) + t441;
t294 = t355 * t374 - t370 * t421;
t295 = t355 * t370 + t374 * t421;
t242 = t294 * t369 + t295 * t373;
t408 = t294 * t373 - t295 * t369;
t572 = t242 * t455 + t365 * t408;
t573 = -t242 * t365 + t408 * t455;
t577 = pkin(5) * t572 - pkin(11) * t573;
t439 = qJD(2) * t375;
t418 = t374 * t439;
t437 = qJD(3) * t374;
t261 = t355 * t437 + (-t371 * t438 + t418) * t442;
t419 = t370 * t439;
t262 = -t355 * t438 + (-t371 * t437 - t419) * t442;
t153 = qJD(4) * t408 + t261 * t373 + t262 * t369;
t432 = qJD(1) * qJD(2);
t414 = t366 * t432;
t406 = t371 * t414;
t279 = pkin(9) * t355 + t318;
t308 = (-pkin(2) * t375 - pkin(9) * t371 - pkin(1)) * t366;
t289 = qJD(1) * t308;
t230 = -t279 * t370 + t289 * t374;
t202 = -pkin(10) * t295 + t230;
t345 = qJD(3) - t420;
t189 = pkin(3) * t345 + t202;
t231 = t279 * t374 + t289 * t370;
t203 = pkin(10) * t294 + t231;
t201 = t373 * t203;
t128 = t189 * t369 + t201;
t316 = qJD(2) * t389;
t302 = qJD(1) * t316;
t451 = t366 * t371;
t356 = pkin(8) * t451;
t503 = pkin(1) * t375;
t329 = t367 * t503 - t356;
t319 = t329 * qJD(2);
t303 = qJD(1) * t319;
t179 = -qJD(3) * t231 + t302 * t374 - t303 * t370;
t138 = pkin(3) * t406 - pkin(10) * t261 + t179;
t178 = -t279 * t438 + t289 * t437 + t302 * t370 + t303 * t374;
t141 = pkin(10) * t262 + t178;
t54 = -qJD(4) * t128 + t138 * t373 - t141 * t369;
t32 = pkin(4) * t406 - qJ(5) * t153 - qJD(5) * t242 + t54;
t154 = -qJD(4) * t242 - t261 * t369 + t262 * t373;
t53 = t138 * t369 + t141 * t373 + t189 * t435 - t203 * t436;
t34 = qJ(5) * t154 + qJD(5) * t408 + t53;
t10 = t32 * t455 - t34 * t365;
t11 = t32 * t365 + t34 * t455;
t238 = Ifges(5,4) * t408;
t339 = qJD(4) + t345;
t168 = Ifges(5,1) * t242 + Ifges(5,5) * t339 + t238;
t157 = t339 * t372 - t368 * t572;
t97 = t153 * t455 + t154 * t365;
t65 = qJD(6) * t157 + t368 * t406 + t372 * t97;
t158 = t339 * t368 + t372 * t572;
t66 = -qJD(6) * t158 - t368 * t97 + t372 * t406;
t96 = t153 * t365 - t154 * t455;
t19 = Ifges(7,4) * t65 + Ifges(7,2) * t66 + Ifges(7,6) * t96;
t199 = t369 * t203;
t127 = t189 * t373 - t199;
t575 = qJ(5) * t242;
t110 = t127 - t575;
t102 = pkin(4) * t339 + t110;
t562 = qJ(5) * t408;
t111 = t128 + t562;
t105 = t455 * t111;
t56 = t102 * t365 + t105;
t51 = pkin(11) * t339 + t56;
t278 = -pkin(2) * t355 - t314;
t243 = -pkin(3) * t294 + t278;
t184 = -pkin(4) * t408 + qJD(5) + t243;
t82 = -pkin(5) * t573 - pkin(11) * t572 + t184;
t22 = -t368 * t51 + t372 * t82;
t450 = t366 * t375;
t330 = pkin(1) * t367 * t371 + pkin(8) * t450;
t304 = t330 * t432;
t229 = -pkin(3) * t262 + t304;
t120 = -pkin(4) * t154 + t229;
t28 = pkin(5) * t96 - pkin(11) * t97 + t120;
t8 = pkin(11) * t406 + t11;
t2 = qJD(6) * t22 + t28 * t368 + t372 * t8;
t20 = Ifges(7,1) * t65 + Ifges(7,4) * t66 + Ifges(7,5) * t96;
t23 = t368 * t82 + t372 * t51;
t3 = -qJD(6) * t23 + t28 * t372 - t368 * t8;
t170 = qJD(6) - t573;
t396 = Ifges(7,5) * t372 - Ifges(7,6) * t368;
t381 = t170 * t396;
t481 = Ifges(7,4) * t368;
t400 = Ifges(7,1) * t372 - t481;
t382 = t158 * t400;
t480 = Ifges(7,4) * t372;
t398 = -Ifges(7,2) * t368 + t480;
t383 = t157 * t398;
t401 = mrSges(7,1) * t368 + mrSges(7,2) * t372;
t453 = t365 * t111;
t55 = t102 * t455 - t453;
t50 = -pkin(5) * t339 - t55;
t384 = t50 * t401;
t395 = Ifges(7,5) * t368 + Ifges(7,6) * t372;
t397 = Ifges(7,2) * t372 + t481;
t399 = Ifges(7,1) * t368 + t480;
t402 = mrSges(7,1) * t372 - mrSges(7,2) * t368;
t433 = qJD(6) * t372;
t434 = qJD(6) * t368;
t485 = Ifges(5,4) * t242;
t490 = mrSges(7,3) * t372;
t491 = mrSges(7,3) * t368;
t505 = t372 / 0.2e1;
t507 = t368 / 0.2e1;
t511 = -t339 / 0.2e1;
t523 = -t242 / 0.2e1;
t525 = -t408 / 0.2e1;
t545 = t96 / 0.2e1;
t156 = Ifges(7,4) * t157;
t77 = Ifges(7,1) * t158 + Ifges(7,5) * t170 + t156;
t546 = t77 / 0.2e1;
t548 = t66 / 0.2e1;
t549 = t65 / 0.2e1;
t559 = Ifges(5,5) * t153 + Ifges(6,5) * t97 + Ifges(5,6) * t154 - Ifges(6,6) * t96 + t406 * t584;
t7 = -pkin(5) * t406 - t10;
t482 = Ifges(7,4) * t158;
t76 = Ifges(7,2) * t157 + Ifges(7,6) * t170 + t482;
t576 = t54 * mrSges(5,1) + t10 * mrSges(6,1) - t53 * mrSges(5,2) - t11 * mrSges(6,2) + qJD(6) * t384 + t19 * t505 + t2 * t490 + t20 * t507 - t3 * t491 + t399 * t549 - t7 * t402 + t397 * t548 - t76 * t434 / 0.2e1 + t433 * t546 + t395 * t545 + t559 + (t383 + t382 + t381) * qJD(6) / 0.2e1 + (Ifges(5,5) * t408 - Ifges(5,6) * t242) * t511 + (-Ifges(5,2) * t242 + t168 + t238) * t525 - t243 * (mrSges(5,1) * t242 + mrSges(5,2) * t408) + (Ifges(5,1) * t408 - t485) * t523;
t502 = pkin(4) * t242;
t63 = Ifges(7,6) * t66;
t64 = Ifges(7,5) * t65;
t18 = Ifges(7,3) * t96 + t63 + t64;
t404 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t496 = t97 * Ifges(6,4);
t574 = t404 + t120 * mrSges(6,1) - t11 * mrSges(6,3) + t18 / 0.2e1 - t496 / 0.2e1;
t475 = Ifges(6,6) * t339;
t476 = Ifges(6,2) * t573;
t483 = Ifges(6,4) * t572;
t108 = t475 + t476 + t483;
t473 = Ifges(7,3) * t170;
t474 = Ifges(7,6) * t157;
t477 = Ifges(7,5) * t158;
t75 = t473 + t474 + t477;
t551 = t56 * mrSges(6,3) + t108 / 0.2e1 - t75 / 0.2e1 - t184 * mrSges(6,1) - t22 * mrSges(7,1) + t23 * mrSges(7,2);
t167 = Ifges(5,2) * t408 + Ifges(5,6) * t339 + t485;
t570 = t167 / 0.2e1;
t568 = t96 * Ifges(6,2);
t394 = t22 * t372 + t23 * t368;
t564 = t394 * mrSges(7,3);
t131 = t202 * t373 - t199;
t116 = t131 - t575;
t130 = -t202 * t369 - t201;
t387 = t130 - t562;
t412 = t455 * t369;
t472 = pkin(3) * qJD(4);
t563 = t116 * t365 - t387 * t455 - (t365 * t373 + t412) * t472;
t164 = mrSges(6,1) * t339 - mrSges(6,3) * t572;
t99 = -mrSges(7,1) * t157 + mrSges(7,2) * t158;
t456 = t99 - t164;
t307 = pkin(9) * t367 + t330;
t248 = -t307 * t370 + t308 * t374;
t325 = t367 * t370 + t374 * t451;
t212 = -pkin(3) * t450 - pkin(10) * t325 + t248;
t249 = t307 * t374 + t308 * t370;
t324 = t367 * t374 - t370 * t451;
t220 = pkin(10) * t324 + t249;
t145 = t212 * t369 + t220 * t373;
t29 = mrSges(7,1) * t96 - mrSges(7,3) * t65;
t30 = -mrSges(7,2) * t96 + mrSges(7,3) * t66;
t558 = -t29 * t368 + t30 * t372;
t557 = -t22 * t368 + t23 * t372;
t463 = t295 * Ifges(4,4);
t227 = Ifges(4,2) * t294 + Ifges(4,6) * t345 + t463;
t290 = Ifges(4,4) * t294;
t228 = Ifges(4,1) * t295 + t345 * Ifges(4,5) + t290;
t390 = t230 * t374 + t231 * t370;
t486 = Ifges(4,4) * t374;
t487 = Ifges(4,4) * t370;
t504 = t374 / 0.2e1;
t509 = t345 / 0.2e1;
t514 = t295 / 0.2e1;
t516 = t294 / 0.2e1;
t555 = -t390 * mrSges(4,3) + t278 * (mrSges(4,1) * t370 + mrSges(4,2) * t374) + (-Ifges(4,2) * t370 + t486) * t516 + (Ifges(4,1) * t374 - t487) * t514 + (Ifges(4,5) * t374 - Ifges(4,6) * t370) * t509 - t370 * t227 / 0.2e1 + t228 * t504;
t478 = Ifges(6,5) * t339;
t484 = Ifges(6,4) * t573;
t489 = Ifges(6,1) * t572;
t109 = t478 + t484 + t489;
t506 = -t372 / 0.2e1;
t552 = -t184 * mrSges(6,2) + t55 * mrSges(6,3) + t77 * t506 + t76 * t507 - t384 - t109 / 0.2e1;
t543 = pkin(1) * mrSges(3,1);
t542 = pkin(1) * mrSges(3,2);
t539 = t153 / 0.2e1;
t538 = t154 / 0.2e1;
t537 = -t157 / 0.2e1;
t536 = t157 / 0.2e1;
t535 = -t158 / 0.2e1;
t534 = t158 / 0.2e1;
t533 = -t170 / 0.2e1;
t532 = t170 / 0.2e1;
t531 = -t573 / 0.2e1;
t530 = t573 / 0.2e1;
t529 = -t572 / 0.2e1;
t528 = t572 / 0.2e1;
t255 = t324 * t373 - t325 * t369;
t256 = t324 * t369 + t325 * t373;
t197 = -t255 * t455 + t256 * t365;
t527 = -t197 / 0.2e1;
t198 = t255 * t365 + t256 * t455;
t526 = t198 / 0.2e1;
t524 = t408 / 0.2e1;
t522 = t242 / 0.2e1;
t521 = t255 / 0.2e1;
t520 = t256 / 0.2e1;
t519 = t261 / 0.2e1;
t518 = t262 / 0.2e1;
t515 = -t295 / 0.2e1;
t513 = t324 / 0.2e1;
t512 = t325 / 0.2e1;
t510 = t339 / 0.2e1;
t508 = -t368 / 0.2e1;
t501 = pkin(4) * t365;
t498 = t96 * Ifges(6,4);
t497 = t97 * Ifges(6,1);
t272 = qJD(3) * t324 + t366 * t418;
t273 = -qJD(3) * t325 - t366 * t419;
t180 = qJD(4) * t255 + t272 * t373 + t273 * t369;
t440 = qJD(2) * t366;
t417 = t371 * t440;
t191 = -qJD(3) * t249 + t316 * t374 - t319 * t370;
t155 = pkin(3) * t417 - pkin(10) * t272 + t191;
t190 = -t307 * t438 + t308 * t437 + t316 * t370 + t319 * t374;
t165 = pkin(10) * t273 + t190;
t68 = -qJD(4) * t145 + t155 * t373 - t165 * t369;
t41 = pkin(4) * t417 - qJ(5) * t180 - qJD(5) * t256 + t68;
t181 = -qJD(4) * t256 - t272 * t369 + t273 * t373;
t67 = t155 * t369 + t165 * t373 + t212 * t435 - t220 * t436;
t45 = qJ(5) * t181 + qJD(5) * t255 + t67;
t15 = t365 * t41 + t45 * t455;
t493 = mrSges(5,3) * t408;
t492 = mrSges(5,3) * t242;
t488 = Ifges(3,4) * t371;
t479 = Ifges(3,5) * t375;
t471 = t573 * Ifges(6,6);
t470 = t572 * Ifges(6,5);
t467 = t408 * Ifges(5,6);
t466 = t242 * Ifges(5,5);
t464 = t294 * Ifges(4,6);
t462 = t295 * Ifges(4,5);
t460 = t303 * mrSges(3,2);
t459 = t345 * Ifges(4,3);
t458 = t355 * Ifges(3,5);
t452 = t365 * t369;
t449 = t368 * t214;
t448 = t372 * t214;
t144 = t212 * t373 - t220 * t369;
t121 = -pkin(4) * t450 - qJ(5) * t256 + t144;
t125 = qJ(5) * t255 + t145;
t72 = t121 * t365 + t125 * t455;
t445 = -mrSges(3,1) * t355 - mrSges(4,1) * t294 + mrSges(4,2) * t295 + mrSges(3,3) * t421;
t444 = t269 - t281;
t320 = t330 * qJD(2);
t362 = pkin(3) * t373 + pkin(4);
t322 = pkin(3) * t412 + t362 * t365;
t428 = -Ifges(6,3) / 0.2e1 - Ifges(5,3) / 0.2e1;
t424 = Ifges(4,5) * t261 + Ifges(4,6) * t262 + Ifges(4,3) * t406;
t422 = t455 * pkin(4);
t46 = mrSges(6,1) * t96 + mrSges(6,2) * t97;
t208 = -t223 * t368 + t372 * t421;
t411 = t208 + t449;
t209 = t223 * t372 + t368 * t421;
t410 = -t209 + t448;
t246 = -pkin(3) * t273 + t320;
t207 = pkin(3) * t295 + t502;
t403 = -t2 * t368 - t3 * t372;
t306 = t356 + (-pkin(2) - t503) * t367;
t260 = -pkin(3) * t324 + t306;
t204 = -pkin(4) * t255 + t260;
t100 = pkin(5) * t197 - pkin(11) * t198 + t204;
t70 = -pkin(11) * t450 + t72;
t35 = t100 * t372 - t368 * t70;
t36 = t100 * t368 + t372 * t70;
t103 = -mrSges(7,2) * t170 + mrSges(7,3) * t157;
t104 = mrSges(7,1) * t170 - mrSges(7,3) * t158;
t392 = t103 * t372 - t104 * t368;
t391 = t178 * t374 - t179 * t370;
t185 = -t198 * t368 - t372 * t450;
t388 = -t198 * t372 + t368 * t450;
t163 = -mrSges(6,2) * t339 + mrSges(6,3) * t573;
t385 = t163 + t392;
t14 = -t365 * t45 + t41 * t455;
t71 = t121 * t455 - t125 * t365;
t83 = t133 * t455 - t139 * t365;
t136 = -pkin(4) * t181 + t246;
t321 = -pkin(3) * t452 + t362 * t455;
t379 = -qJD(6) * t394 + t2 * t372 - t3 * t368;
t361 = -t422 - pkin(5);
t349 = Ifges(3,4) * t420;
t344 = t414 * t479;
t317 = (t373 * t455 - t452) * t472;
t312 = pkin(11) + t322;
t311 = -pkin(5) - t321;
t310 = -mrSges(3,2) * t355 + mrSges(3,3) * t420;
t276 = Ifges(3,1) * t421 + t349 + t458;
t275 = Ifges(3,6) * t355 + (Ifges(3,2) * t375 + t488) * t442;
t264 = mrSges(4,1) * t345 - mrSges(4,3) * t295;
t263 = -mrSges(4,2) * t345 + mrSges(4,3) * t294;
t245 = -mrSges(4,2) * t406 + mrSges(4,3) * t262;
t244 = mrSges(4,1) * t406 - mrSges(4,3) * t261;
t226 = t459 + t462 + t464;
t219 = mrSges(5,1) * t339 - t492;
t218 = -mrSges(5,2) * t339 + t493;
t205 = -mrSges(4,1) * t262 + mrSges(4,2) * t261;
t196 = Ifges(4,1) * t261 + Ifges(4,4) * t262 + Ifges(4,5) * t406;
t195 = Ifges(4,4) * t261 + Ifges(4,2) * t262 + Ifges(4,6) * t406;
t193 = t254 * t365 - t386 * t455;
t183 = -mrSges(5,1) * t408 + mrSges(5,2) * t242;
t166 = t339 * Ifges(5,3) + t466 + t467;
t143 = -mrSges(5,2) * t406 + mrSges(5,3) * t154;
t142 = mrSges(5,1) * t406 - mrSges(5,3) * t153;
t117 = t182 * t365 - t378 * t455;
t115 = t180 * t455 + t181 * t365;
t114 = t180 * t365 - t181 * t455;
t113 = -mrSges(6,1) * t573 + mrSges(6,2) * t572;
t107 = t339 * Ifges(6,3) + t470 + t471;
t98 = -mrSges(5,1) * t154 + mrSges(5,2) * t153;
t95 = t502 + t577;
t90 = Ifges(5,1) * t153 + Ifges(5,4) * t154 + Ifges(5,5) * t406;
t89 = Ifges(5,4) * t153 + Ifges(5,2) * t154 + Ifges(5,6) * t406;
t87 = mrSges(6,1) * t406 - mrSges(6,3) * t97;
t86 = -mrSges(6,2) * t406 - mrSges(6,3) * t96;
t85 = t207 + t577;
t80 = -pkin(5) * t421 - t83;
t79 = qJD(6) * t388 - t115 * t368 + t372 * t417;
t78 = qJD(6) * t185 + t115 * t372 + t368 * t417;
t69 = pkin(5) * t450 - t71;
t60 = t116 * t455 + t365 * t387;
t58 = t110 * t455 - t453;
t57 = t110 * t365 + t105;
t44 = Ifges(6,5) * t406 + t497 - t498;
t43 = Ifges(6,6) * t406 + t496 - t568;
t37 = pkin(5) * t114 - pkin(11) * t115 + t136;
t27 = t368 * t95 + t372 * t58;
t26 = -t368 * t58 + t372 * t95;
t25 = t368 * t85 + t372 * t60;
t24 = -t368 * t60 + t372 * t85;
t21 = -mrSges(7,1) * t66 + mrSges(7,2) * t65;
t13 = pkin(11) * t417 + t15;
t12 = -pkin(5) * t417 - t14;
t5 = -qJD(6) * t36 - t13 * t368 + t37 * t372;
t4 = qJD(6) * t35 + t13 * t372 + t368 * t37;
t1 = [-(t424 + t559) * t450 / 0.2e1 + (Ifges(7,5) * t78 + Ifges(7,6) * t79) * t532 + (Ifges(7,4) * t78 + Ifges(7,2) * t79) * t536 - t96 * (Ifges(6,4) * t198 - Ifges(6,6) * t450) / 0.2e1 + (Ifges(6,4) * t115 + Ifges(6,6) * t417) * t530 + (Ifges(7,1) * t78 + Ifges(7,4) * t79) * t534 + t97 * (Ifges(6,1) * t198 - Ifges(6,5) * t450) / 0.2e1 + (Ifges(6,1) * t115 + Ifges(6,5) * t417) * t528 + (t11 * t450 + t115 * t184 + t120 * t198 - t417 * t56) * mrSges(6,2) + (t568 / 0.2e1 + Ifges(7,3) * t545 + Ifges(7,6) * t548 + Ifges(7,5) * t549 + t574) * t197 + (-Ifges(6,4) * t528 + Ifges(7,5) * t534 - Ifges(6,2) * t530 - Ifges(6,6) * t510 + Ifges(7,6) * t536 + Ifges(7,3) * t532 - t551) * t114 + (t344 / 0.2e1 - t304 * mrSges(3,1) - t460) * t367 + (t303 * t375 + t304 * t371 + (-t314 * t375 - t318 * t371) * qJD(2)) * t366 * mrSges(3,3) + t15 * t163 + t14 * t164 + t10 * (-mrSges(6,1) * t450 - mrSges(6,3) * t198) + t54 * (-mrSges(5,1) * t450 - mrSges(5,3) * t256) + t179 * (-mrSges(4,1) * t450 - mrSges(4,3) * t325) + t53 * (mrSges(5,2) * t450 + mrSges(5,3) * t255) + t178 * (mrSges(4,2) * t450 + mrSges(4,3) * t324) + t445 * t320 + ((-t329 * mrSges(3,3) + Ifges(3,5) * t367 / 0.2e1 + (-0.2e1 * t542 + 0.3e1 / 0.2e1 * Ifges(3,4) * t375) * t366) * t375 + (-t330 * mrSges(3,3) - Ifges(3,6) * t367 + Ifges(4,5) * t512 + Ifges(4,6) * t513 + Ifges(5,5) * t520 + Ifges(5,6) * t521 + Ifges(6,5) * t526 + Ifges(6,6) * t527 + (-0.2e1 * t543 - 0.3e1 / 0.2e1 * t488) * t366 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - Ifges(4,3) / 0.2e1 + t428) * t450) * t371) * t414 + t145 * t143 + t144 * t142 + t136 * t113 + t115 * t109 / 0.2e1 + t55 * (mrSges(6,1) * t417 - mrSges(6,3) * t115) + t127 * (mrSges(5,1) * t417 - mrSges(5,3) * t180) + t230 * (mrSges(4,1) * t417 - mrSges(4,3) * t272) + t4 * t103 + t5 * t104 - t275 * t417 / 0.2e1 + t128 * (-mrSges(5,2) * t417 + mrSges(5,3) * t181) + t231 * (-mrSges(4,2) * t417 + mrSges(4,3) * t273) + t12 * t99 + t72 * t86 + t71 * t87 + t50 * (-mrSges(7,1) * t79 + mrSges(7,2) * t78) + t79 * t76 / 0.2e1 + t69 * t21 + m(4) * (t178 * t249 + t179 * t248 + t190 * t231 + t191 * t230 + t278 * t320 + t304 * t306) + m(3) * (t303 * t330 - t304 * t329 - t314 * t320 + t318 * t319) + m(5) * (t127 * t68 + t128 * t67 + t144 * t54 + t145 * t53 + t229 * t260 + t243 * t246) + m(7) * (t12 * t50 + t2 * t36 + t22 * t5 + t23 * t4 + t3 * t35 + t69 * t7) + m(6) * (t10 * t71 + t11 * t72 + t120 * t204 + t136 * t184 + t14 * t55 + t15 * t56) + t35 * t29 + t36 * t30 + (Ifges(5,5) * t180 + Ifges(6,5) * t115 + Ifges(5,6) * t181 + t417 * t584) * t510 + t180 * t168 / 0.2e1 + t185 * t19 / 0.2e1 + t204 * t46 + t67 * t218 + t68 * t219 + t181 * t570 + (Ifges(4,5) * t272 + Ifges(4,6) * t273 + Ifges(4,3) * t417) * t509 + t196 * t512 + t195 * t513 + (Ifges(4,1) * t272 + Ifges(4,4) * t273 + Ifges(4,5) * t417) * t514 + (Ifges(4,4) * t272 + Ifges(4,2) * t273 + Ifges(4,6) * t417) * t516 + (Ifges(4,4) * t325 + Ifges(4,2) * t324 - Ifges(4,6) * t450) * t518 + (Ifges(4,1) * t325 + Ifges(4,4) * t324 - Ifges(4,5) * t450) * t519 + t90 * t520 + t89 * t521 + (Ifges(5,1) * t180 + Ifges(5,4) * t181 + Ifges(5,5) * t417) * t522 + (Ifges(5,4) * t180 + Ifges(5,2) * t181 + Ifges(5,6) * t417) * t524 + t243 * (-mrSges(5,1) * t181 + mrSges(5,2) * t180) + t246 * t183 + t248 * t244 + t249 * t245 + t229 * (-mrSges(5,1) * t255 + mrSges(5,2) * t256) + t260 * t98 + t190 * t263 + t191 * t264 + ((t226 + t166 + t107) * t371 + t375 * t276 + t355 * (-Ifges(3,6) * t371 + t479)) * t440 / 0.2e1 + (-Ifges(7,5) * t388 + Ifges(7,6) * t185) * t545 + (-Ifges(7,4) * t388 + Ifges(7,2) * t185) * t548 + (-Ifges(7,1) * t388 + Ifges(7,4) * t185) * t549 + (t185 * t2 - t22 * t78 + t23 * t79 + t3 * t388) * mrSges(7,3) + t7 * (-mrSges(7,1) * t185 - mrSges(7,2) * t388) - t388 * t20 / 0.2e1 + t272 * t228 / 0.2e1 + t273 * t227 / 0.2e1 + t278 * (-mrSges(4,1) * t273 + mrSges(4,2) * t272) + t306 * t205 + t319 * t310 + t304 * (-mrSges(4,1) * t324 + mrSges(4,2) * t325) + t44 * t526 + t43 * t527 + (Ifges(5,4) * t256 + Ifges(5,2) * t255 - Ifges(5,6) * t450) * t538 + (Ifges(5,1) * t256 + Ifges(5,4) * t255 - Ifges(5,5) * t450) * t539 + t78 * t546; (t20 * t505 + t19 * t508 - t10 * mrSges(6,3) + t44 / 0.2e1 + t120 * mrSges(6,2) + t497 / 0.2e1 - t498 / 0.2e1 + t7 * t401 + t396 * t545 + t400 * t549 + t398 * t548 + t403 * mrSges(7,3) + (-mrSges(7,3) * t557 + t395 * t533 + t397 * t537 + t399 * t535 + t402 * t50 + t506 * t76 + t508 * t77) * qJD(6)) * t267 + (t75 - t108) * (t213 / 0.2e1 - t222 / 0.2e1) + (t580 * t127 + t579 * t128 + t229 * t363 + t243 * t587 + t282 * t54 + t283 * t53) * m(5) + (t269 / 0.2e1 - t281 / 0.2e1) * t168 + (t64 / 0.2e1 + t63 / 0.2e1 - t43 / 0.2e1 + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t96 + t574) * t266 + (-t270 / 0.2e1 + t280 / 0.2e1) * t167 + (Ifges(5,1) * t269 - Ifges(5,4) * t270) * t522 + (Ifges(5,4) * t269 - Ifges(5,2) * t270) * t524 + (Ifges(5,5) * t269 + Ifges(6,5) * t214 - Ifges(5,6) * t270 - Ifges(6,6) * t213) * t510 + (-t244 * t370 + t245 * t374) * pkin(9) + (t21 - t87) * t193 + t344 - t460 - t83 * t164 + (t370 * pkin(3) * t183 + (-t263 * t370 - t264 * t374) * pkin(9) + t555) * qJD(3) + (mrSges(5,1) * t588 + mrSges(5,2) * t444) * t243 + (-t127 * t444 - t128 * t588 + t332 * t53 - t333 * t54) * mrSges(5,3) + (t214 / 0.2e1 - t223 / 0.2e1) * t109 + t456 * t117 + (-t449 / 0.2e1 - t208 / 0.2e1) * t76 + (t448 / 0.2e1 - t209 / 0.2e1) * t77 + (mrSges(7,1) * t447 - mrSges(7,3) * t410) * t22 + (-t446 * t55 - t447 * t56) * mrSges(6,3) + (-mrSges(7,2) * t447 - mrSges(7,3) * t411) * t23 - t445 * t318 + (mrSges(6,1) * t447 + mrSges(6,2) * t446) * t184 + t123 * t29 + t124 * t30 - t80 * t99 + (mrSges(7,1) * t411 + mrSges(7,2) * t410) * t50 + t582 * t103 + (t123 * t3 + t124 * t2 + t193 * t7 + (t117 - t80) * t50 + t582 * t23 + t583 * t22) * m(7) + t583 * t104 + t581 * t163 + (-t10 * t193 + t11 * t194 + t120 * t305 + t581 * t56 + (-t117 - t83) * t55 + t578 * t184) * m(6) + (-t230 * t251 - t231 * t252 - t278 * t318 - pkin(2) * t304 + (-qJD(3) * t390 + t391) * pkin(9)) * m(4) + t578 * t113 + t579 * t218 + t580 * t219 + (-mrSges(4,1) * t374 + mrSges(4,2) * t370 - mrSges(3,1)) * t304 + (Ifges(5,1) * t281 - Ifges(5,4) * t280) * t523 + (Ifges(5,4) * t281 - Ifges(5,2) * t280) * t525 + (Ifges(5,5) * t281 + Ifges(6,5) * t223 - Ifges(5,6) * t280 - Ifges(6,6) * t222) * t511 + t194 * t86 - pkin(2) * t205 + t195 * t504 + (Ifges(4,2) * t374 + t487) * t518 + (Ifges(4,1) * t370 + t486) * t519 - t252 * t263 - t251 * t264 - t271 * t183 + t282 * t142 + t283 * t143 + t305 * t46 - t314 * t310 + t332 * t89 / 0.2e1 + t333 * t90 / 0.2e1 + t229 * (-mrSges(5,1) * t332 + mrSges(5,2) * t333) + ((-t458 / 0.2e1 - t349 / 0.2e1 - t276 / 0.2e1 + t314 * mrSges(3,3) + t442 * t542 - t555) * t375 + (-t166 / 0.2e1 - t226 / 0.2e1 + t318 * mrSges(3,3) + t56 * mrSges(6,2) - t55 * mrSges(6,1) - t471 / 0.2e1 - t470 / 0.2e1 - t127 * mrSges(5,1) + t128 * mrSges(5,2) - t467 / 0.2e1 - t466 / 0.2e1 + t231 * mrSges(4,2) - t230 * mrSges(4,1) - t464 / 0.2e1 - t462 / 0.2e1 - t459 / 0.2e1 + t275 / 0.2e1 - t107 / 0.2e1 + t428 * t339 + (t543 + t488 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t375) * t442 + (-qJD(2) + t355 / 0.2e1) * Ifges(3,6) + (Ifges(4,5) * t370 + Ifges(5,5) * t333 + Ifges(6,5) * t267 + Ifges(4,6) * t374 + Ifges(5,6) * t332 - Ifges(6,6) * t266) * qJD(2) / 0.2e1) * t371) * t442 + t363 * t98 + t370 * t196 / 0.2e1 + (Ifges(6,1) * t214 - Ifges(6,4) * t213) * t528 + (Ifges(6,1) * t223 - Ifges(6,4) * t222) * t529 + (Ifges(6,4) * t214 - Ifges(6,2) * t213) * t530 + (Ifges(6,4) * t223 - Ifges(6,2) * t222) * t531 + (Ifges(7,5) * t448 - Ifges(7,6) * t449 + Ifges(7,3) * t213) * t532 + (Ifges(7,5) * t209 + Ifges(7,6) * t208 + Ifges(7,3) * t222) * t533 + (Ifges(7,1) * t448 - Ifges(7,4) * t449 + Ifges(7,5) * t213) * t534 + (Ifges(7,1) * t209 + Ifges(7,4) * t208 + Ifges(7,5) * t222) * t535 + (Ifges(7,4) * t448 - Ifges(7,2) * t449 + Ifges(7,6) * t213) * t536 + (Ifges(7,4) * t209 + Ifges(7,2) * t208 + Ifges(7,6) * t222) * t537 + (Ifges(5,4) * t333 + Ifges(5,2) * t332) * t538 + (Ifges(5,1) * t333 + Ifges(5,4) * t332) * t539 + t391 * mrSges(4,3); (t476 / 0.2e1 + t483 / 0.2e1 + t475 / 0.2e1 - t473 / 0.2e1 - t474 / 0.2e1 - t477 / 0.2e1 + t551) * t572 + (-t484 / 0.2e1 - t489 / 0.2e1 - t478 / 0.2e1 - t381 / 0.2e1 - t383 / 0.2e1 - t382 / 0.2e1 + t564 + t552) * t573 - t456 * t563 + (t10 * t321 + t11 * t322 - t184 * t207 + (t317 - t60) * t56 + t563 * t55) * m(6) + (-t22 * t24 - t23 * t25 + t311 * t7 + t312 * t379 + t317 * t557 - t50 * t563) * m(7) + ((-t103 * t368 - t104 * t372) * t312 - t564) * qJD(6) + t558 * t312 + (t373 * t142 + t369 * t143 - t295 * t183 + (t218 * t373 - t219 * t369) * qJD(4) + (-t127 * t436 + t128 * t435 + 0.2e1 * t243 * t515 + t369 * t53 + t373 * t54) * m(5)) * pkin(3) + t576 + (t230 * t294 + t231 * t295) * mrSges(4,3) - m(5) * (t127 * t130 + t128 * t131) - (-Ifges(4,2) * t295 + t228 + t290) * t294 / 0.2e1 - t60 * t163 - t25 * t103 - t24 * t104 - t178 * mrSges(4,2) + t179 * mrSges(4,1) + t424 - t207 * t113 - t131 * t218 - t130 * t219 + t227 * t514 + (Ifges(4,1) * t294 - t463) * t515 - t230 * t263 + t231 * t264 + (t127 * t408 + t128 * t242) * mrSges(5,3) + t242 * t570 - t278 * (mrSges(4,1) * t295 + mrSges(4,2) * t294) + t311 * t21 + t321 * t87 + t322 * t86 - t345 * (Ifges(4,5) * t294 - Ifges(4,6) * t295) / 0.2e1 + t385 * t317; (Ifges(6,1) * t529 + Ifges(6,4) * t531 + Ifges(6,5) * t511 + t22 * t490 + t23 * t491 + t396 * t533 + t398 * t537 + t400 * t535 + t552) * t573 - t456 * t57 + (m(7) * t379 - t103 * t434 - t104 * t433 + t558) * (pkin(11) + t501) + (-t22 * t26 - t23 * t27 + t361 * t7 - t50 * t57) * m(7) + (-t22 * t433 - t23 * t434) * mrSges(7,3) + (-t218 + t493) * t127 + (-t184 * t502 + t55 * t57 - t56 * t58 + (t10 * t455 + t11 * t365) * pkin(4)) * m(6) - t113 * t502 - t58 * t163 + (t219 + t492) * t128 + (-Ifges(6,4) * t529 + Ifges(7,5) * t535 - Ifges(6,2) * t531 - Ifges(6,6) * t511 + Ifges(7,6) * t537 + Ifges(7,3) * t533 + t551) * t572 - t27 * t103 - t26 * t104 + t87 * t422 + t86 * t501 + t167 * t522 + t361 * t21 + t576; t392 * qJD(6) - t456 * t572 - t385 * t573 + t372 * t29 + t368 * t30 + t46 + (t170 * t557 - t50 * t572 - t403) * m(7) + (t55 * t572 - t56 * t573 + t120) * m(6); -t50 * (mrSges(7,1) * t158 + mrSges(7,2) * t157) + (Ifges(7,1) * t157 - t482) * t535 + t76 * t534 + (Ifges(7,5) * t157 - Ifges(7,6) * t158) * t533 - t22 * t103 + t23 * t104 + (t157 * t22 + t158 * t23) * mrSges(7,3) + t404 + t18 + (-Ifges(7,2) * t158 + t156 + t77) * t537;];
tauc  = t1(:);
