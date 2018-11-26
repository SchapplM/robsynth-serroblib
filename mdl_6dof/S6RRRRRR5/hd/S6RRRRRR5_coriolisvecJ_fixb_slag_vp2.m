% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:40:22
% EndTime: 2018-11-23 18:40:57
% DurationCPUTime: 35.08s
% Computational Cost: add. (43815->969), mult. (111685->1366), div. (0->0), fcn. (90406->12), ass. (0->432)
t374 = sin(qJ(4));
t375 = sin(qJ(3));
t379 = cos(qJ(4));
t380 = cos(qJ(3));
t338 = t374 * t380 + t375 * t379;
t562 = qJD(3) + qJD(4);
t275 = t562 * t338;
t381 = cos(qJ(2));
t370 = sin(pkin(6));
t457 = qJD(1) * t370;
t434 = t381 * t457;
t286 = t338 * t434;
t595 = t275 - t286;
t337 = -t374 * t375 + t379 * t380;
t274 = t562 * t337;
t373 = sin(qJ(5));
t378 = cos(qJ(5));
t400 = t378 * t337 - t338 * t373;
t185 = qJD(5) * t400 + t274 * t378 - t275 * t373;
t287 = t337 * t434;
t226 = -t286 * t373 + t287 * t378;
t462 = t185 - t226;
t273 = t337 * t373 + t338 * t378;
t186 = qJD(5) * t273 + t274 * t373 + t378 * t275;
t225 = t378 * t286 + t287 * t373;
t461 = t186 - t225;
t376 = sin(qJ(2));
t371 = cos(pkin(6));
t456 = qJD(1) * t371;
t441 = pkin(1) * t456;
t320 = pkin(8) * t434 + t376 * t441;
t422 = t375 * t434;
t277 = pkin(3) * t422 + t320;
t453 = qJD(3) * t375;
t594 = pkin(3) * t453 - t277;
t435 = t376 * t457;
t317 = -pkin(8) * t435 + t381 * t441;
t398 = (pkin(2) * t376 - pkin(9) * t381) * t370;
t318 = qJD(1) * t398;
t253 = -t375 * t317 + t380 * t318;
t224 = (-pkin(10) * t380 * t381 + pkin(3) * t376) * t457 + t253;
t254 = t380 * t317 + t375 * t318;
t237 = -pkin(10) * t422 + t254;
t164 = t379 * t224 - t237 * t374;
t135 = pkin(4) * t435 - pkin(11) * t287 + t164;
t165 = t374 * t224 + t379 * t237;
t141 = -pkin(11) * t286 + t165;
t551 = -pkin(10) - pkin(9);
t436 = qJD(3) * t551;
t342 = t375 * t436;
t343 = t380 * t436;
t351 = t551 * t375;
t352 = t551 * t380;
t450 = qJD(4) * t379;
t451 = qJD(4) * t374;
t227 = t379 * t342 + t374 * t343 + t351 * t450 + t352 * t451;
t196 = -pkin(11) * t275 + t227;
t289 = t374 * t351 - t379 * t352;
t228 = -qJD(4) * t289 - t342 * t374 + t379 * t343;
t390 = -pkin(11) * t274 + t228;
t288 = t379 * t351 + t352 * t374;
t257 = -pkin(11) * t338 + t288;
t258 = pkin(11) * t337 + t289;
t402 = t378 * t257 - t258 * t373;
t588 = qJD(5) * t402 + (-t141 + t196) * t378 + (-t135 + t390) * t373;
t585 = pkin(4) * t595 + t594;
t593 = -pkin(12) * t435 + t588;
t592 = t461 * pkin(5) - pkin(12) * t462 + t585;
t368 = -pkin(3) * t380 - pkin(2);
t311 = -pkin(4) * t337 + t368;
t197 = -pkin(5) * t400 - pkin(12) * t273 + t311;
t199 = t257 * t373 + t258 * t378;
t372 = sin(qJ(6));
t377 = cos(qJ(6));
t128 = t197 * t372 + t199 * t377;
t591 = -qJD(6) * t128 - t372 * t593 + t592 * t377;
t127 = t197 * t377 - t199 * t372;
t590 = qJD(6) * t127 + t592 * t372 + t377 * t593;
t350 = qJD(3) - t434;
t344 = qJD(4) + t350;
t336 = qJD(5) + t344;
t360 = qJD(2) + t456;
t300 = t360 * t380 - t375 * t435;
t301 = t360 * t375 + t380 * t435;
t245 = t300 * t374 + t301 * t379;
t423 = t379 * t300 - t301 * t374;
t578 = t378 * t245 + t373 * t423;
t159 = t336 * t377 - t372 * t578;
t444 = qJD(1) * qJD(2);
t428 = t370 * t444;
t421 = t376 * t428;
t454 = qJD(2) * t381;
t431 = t380 * t454;
t452 = qJD(3) * t380;
t265 = t360 * t452 + (-t376 * t453 + t431) * t457;
t432 = t375 * t454;
t266 = -t360 * t453 + (-t376 * t452 - t432) * t457;
t155 = qJD(4) * t423 + t265 * t379 + t266 * t374;
t156 = -qJD(4) * t245 - t265 * t374 + t266 * t379;
t579 = -t245 * t373 + t378 * t423;
t81 = qJD(5) * t579 + t155 * t378 + t156 * t373;
t49 = qJD(6) * t159 + t372 * t421 + t377 * t81;
t160 = t336 * t372 + t377 * t578;
t50 = -qJD(6) * t160 - t372 * t81 + t377 * t421;
t21 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t285 = pkin(9) * t360 + t320;
t314 = (-pkin(2) * t381 - pkin(9) * t376 - pkin(1)) * t370;
t295 = qJD(1) * t314;
t233 = -t285 * t375 + t380 * t295;
t207 = -pkin(10) * t301 + t233;
t193 = pkin(3) * t350 + t207;
t234 = t285 * t380 + t295 * t375;
t208 = pkin(10) * t300 + t234;
t206 = t379 * t208;
t131 = t193 * t374 + t206;
t319 = qJD(2) * t398;
t308 = qJD(1) * t319;
t470 = t370 * t376;
t361 = pkin(8) * t470;
t513 = pkin(1) * t381;
t333 = t371 * t513 - t361;
t321 = t333 * qJD(2);
t309 = qJD(1) * t321;
t175 = -qJD(3) * t234 + t380 * t308 - t309 * t375;
t139 = pkin(3) * t421 - pkin(10) * t265 + t175;
t174 = -t285 * t453 + t295 * t452 + t375 * t308 + t380 * t309;
t143 = pkin(10) * t266 + t174;
t62 = -qJD(4) * t131 + t379 * t139 - t143 * t374;
t42 = pkin(4) * t421 - pkin(11) * t155 + t62;
t61 = t374 * t139 + t379 * t143 + t193 * t450 - t208 * t451;
t44 = pkin(11) * t156 + t61;
t204 = t374 * t208;
t130 = t379 * t193 - t204;
t582 = pkin(11) * t245;
t115 = t130 - t582;
t108 = pkin(4) * t344 + t115;
t574 = pkin(11) * t423;
t116 = t131 + t574;
t464 = t378 * t116;
t64 = t108 * t373 + t464;
t11 = -qJD(5) * t64 - t373 * t44 + t378 * t42;
t8 = -pkin(5) * t421 - t11;
t589 = m(7) * t8 + t21;
t587 = -t164 + t228;
t586 = -t165 + t227;
t120 = pkin(5) * t578 - pkin(12) * t579;
t241 = Ifges(5,4) * t423;
t170 = t245 * Ifges(5,1) + t344 * Ifges(5,5) + t241;
t284 = -t360 * pkin(2) - t317;
t246 = -t300 * pkin(3) + t284;
t173 = qJD(6) - t579;
t490 = Ifges(7,3) * t173;
t491 = Ifges(7,6) * t159;
t492 = Ifges(6,6) * t336;
t494 = Ifges(7,5) * t160;
t500 = Ifges(6,4) * t578;
t493 = Ifges(6,2) * t579;
t113 = t492 + t493 + t500;
t187 = -pkin(4) * t423 + t246;
t59 = pkin(12) * t336 + t64;
t92 = -pkin(5) * t579 - pkin(12) * t578 + t187;
t25 = -t372 * t59 + t377 * t92;
t26 = t372 * t92 + t377 * t59;
t87 = t490 + t491 + t494;
t559 = t26 * mrSges(7,2) + t64 * mrSges(6,3) + t113 / 0.2e1 - t87 / 0.2e1 - t187 * mrSges(6,1) - t25 * mrSges(7,1);
t386 = -t494 / 0.2e1 - t490 / 0.2e1 + t492 / 0.2e1 - t491 / 0.2e1 + t500 / 0.2e1 + t559;
t448 = qJD(5) * t378;
t449 = qJD(5) * t373;
t10 = t108 * t448 - t116 * t449 + t373 * t42 + t378 * t44;
t82 = qJD(5) * t578 + t155 * t373 - t378 * t156;
t15 = t49 * Ifges(7,4) + t50 * Ifges(7,2) + t82 * Ifges(7,6);
t16 = t49 * Ifges(7,1) + t50 * Ifges(7,4) + t82 * Ifges(7,5);
t412 = Ifges(7,5) * t377 - Ifges(7,6) * t372;
t392 = t173 * t412;
t498 = Ifges(7,4) * t372;
t416 = Ifges(7,1) * t377 - t498;
t393 = t160 * t416;
t497 = Ifges(7,4) * t377;
t414 = -Ifges(7,2) * t372 + t497;
t394 = t159 * t414;
t417 = mrSges(7,1) * t372 + mrSges(7,2) * t377;
t473 = t116 * t373;
t63 = t108 * t378 - t473;
t58 = -pkin(5) * t336 - t63;
t395 = t58 * t417;
t411 = Ifges(7,5) * t372 + Ifges(7,6) * t377;
t413 = Ifges(7,2) * t377 + t498;
t415 = Ifges(7,1) * t372 + t497;
t418 = mrSges(7,1) * t377 - mrSges(7,2) * t372;
t442 = Ifges(6,5) * t81 - Ifges(6,6) * t82 + Ifges(6,3) * t421;
t446 = qJD(6) * t377;
t447 = qJD(6) * t372;
t469 = t370 * t381;
t334 = t371 * t376 * pkin(1) + pkin(8) * t469;
t310 = t334 * t444;
t232 = -t266 * pkin(3) + t310;
t123 = -t156 * pkin(4) + t232;
t22 = t82 * pkin(5) - t81 * pkin(12) + t123;
t7 = pkin(12) * t421 + t10;
t2 = qJD(6) * t25 + t22 * t372 + t377 * t7;
t512 = t2 * t377;
t515 = t377 / 0.2e1;
t517 = t372 / 0.2e1;
t158 = Ifges(7,4) * t159;
t89 = Ifges(7,1) * t160 + Ifges(7,5) * t173 + t158;
t552 = t89 / 0.2e1;
t553 = t82 / 0.2e1;
t554 = t50 / 0.2e1;
t555 = t49 / 0.2e1;
t499 = Ifges(7,4) * t160;
t88 = Ifges(7,2) * t159 + Ifges(7,6) * t173 + t499;
t389 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + mrSges(7,3) * t512 + qJD(6) * t395 + t15 * t515 + t16 * t517 + t411 * t553 + t413 * t554 + t415 * t555 - t8 * t418 + t442 - t88 * t447 / 0.2e1 + t446 * t552 + (t394 + t393 + t392) * qJD(6) / 0.2e1;
t438 = Ifges(5,5) * t155 + Ifges(5,6) * t156 + Ifges(5,3) * t421;
t501 = Ifges(5,4) * t245;
t521 = -t344 / 0.2e1;
t534 = -t245 / 0.2e1;
t536 = -t423 / 0.2e1;
t584 = t62 * mrSges(5,1) - t61 * mrSges(5,2) + t389 + t438 + (Ifges(5,5) * t423 - Ifges(5,6) * t245) * t521 + (t130 * t423 + t131 * t245) * mrSges(5,3) - t246 * (mrSges(5,1) * t245 + mrSges(5,2) * t423) + (t493 / 0.2e1 + t386) * t578 + (-Ifges(5,2) * t245 + t170 + t241) * t536 + (Ifges(5,1) * t423 - t501) * t534;
t583 = pkin(4) * t245;
t47 = Ifges(7,6) * t50;
t48 = Ifges(7,5) * t49;
t14 = Ifges(7,3) * t82 + t47 + t48;
t3 = -qJD(6) * t26 + t22 * t377 - t372 * t7;
t419 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t509 = t81 * Ifges(6,4);
t581 = t419 + t123 * mrSges(6,1) - t10 * mrSges(6,3) + t14 / 0.2e1 - t509 / 0.2e1;
t169 = Ifges(5,2) * t423 + t344 * Ifges(5,6) + t501;
t576 = t169 / 0.2e1;
t573 = t82 * Ifges(6,2);
t484 = t25 * t377;
t409 = t26 * t372 + t484;
t569 = t409 * mrSges(7,3);
t133 = t379 * t207 - t204;
t118 = t133 - t582;
t367 = pkin(3) * t379 + pkin(4);
t132 = -t207 * t374 - t206;
t399 = t132 - t574;
t466 = t374 * t378;
t568 = t118 * t373 - t378 * t399 - t367 * t449 - (t374 * t448 + (t373 * t379 + t466) * qJD(4)) * pkin(3);
t104 = -mrSges(7,1) * t159 + mrSges(7,2) * t160;
t163 = mrSges(6,1) * t336 - mrSges(6,3) * t578;
t463 = t104 - t163;
t313 = pkin(9) * t371 + t334;
t251 = -t375 * t313 + t380 * t314;
t328 = t371 * t375 + t380 * t470;
t217 = -pkin(3) * t469 - t328 * pkin(10) + t251;
t252 = t380 * t313 + t375 * t314;
t327 = t371 * t380 - t375 * t470;
t223 = pkin(10) * t327 + t252;
t146 = t379 * t217 - t374 * t223;
t260 = t327 * t374 + t328 * t379;
t124 = -pkin(4) * t469 - t260 * pkin(11) + t146;
t147 = t374 * t217 + t379 * t223;
t259 = t327 * t379 - t328 * t374;
t126 = pkin(11) * t259 + t147;
t565 = t373 * t124 + t378 * t126;
t172 = Ifges(6,4) * t579;
t439 = -t172 / 0.2e1;
t505 = Ifges(6,1) * t578;
t564 = t439 - t505 / 0.2e1;
t563 = -t25 * t372 + t26 * t377;
t481 = t301 * Ifges(4,4);
t230 = t300 * Ifges(4,2) + t350 * Ifges(4,6) + t481;
t296 = Ifges(4,4) * t300;
t231 = t301 * Ifges(4,1) + t350 * Ifges(4,5) + t296;
t404 = t233 * t380 + t234 * t375;
t502 = Ifges(4,4) * t380;
t503 = Ifges(4,4) * t375;
t514 = t380 / 0.2e1;
t519 = t350 / 0.2e1;
t525 = t301 / 0.2e1;
t527 = t300 / 0.2e1;
t561 = -t404 * mrSges(4,3) + t284 * (mrSges(4,1) * t375 + mrSges(4,2) * t380) + (-Ifges(4,2) * t375 + t502) * t527 + (Ifges(4,1) * t380 - t503) * t525 + (Ifges(4,5) * t380 - Ifges(4,6) * t375) * t519 - t375 * t230 / 0.2e1 + t231 * t514;
t278 = qJD(3) * t327 + t370 * t431;
t279 = -qJD(3) * t328 - t370 * t432;
t182 = qJD(4) * t259 + t278 * t379 + t279 * t374;
t455 = qJD(2) * t370;
t433 = t376 * t455;
t195 = -qJD(3) * t252 + t380 * t319 - t321 * t375;
t157 = pkin(3) * t433 - pkin(10) * t278 + t195;
t194 = -t313 * t453 + t314 * t452 + t375 * t319 + t380 * t321;
t167 = pkin(10) * t279 + t194;
t72 = -qJD(4) * t147 + t379 * t157 - t167 * t374;
t52 = pkin(4) * t433 - pkin(11) * t182 + t72;
t183 = -qJD(4) * t260 - t278 * t374 + t279 * t379;
t71 = t374 * t157 + t379 * t167 + t217 * t450 - t223 * t451;
t56 = pkin(11) * t183 + t71;
t20 = -qJD(5) * t565 - t373 * t56 + t378 * t52;
t557 = Ifges(6,2) / 0.2e1;
t550 = pkin(1) * mrSges(3,1);
t549 = pkin(1) * mrSges(3,2);
t548 = t155 / 0.2e1;
t547 = t156 / 0.2e1;
t546 = -t159 / 0.2e1;
t545 = t159 / 0.2e1;
t544 = -t160 / 0.2e1;
t543 = t160 / 0.2e1;
t542 = -t173 / 0.2e1;
t541 = t173 / 0.2e1;
t540 = t579 / 0.2e1;
t539 = t578 / 0.2e1;
t401 = t378 * t259 - t260 * t373;
t538 = t401 / 0.2e1;
t203 = t259 * t373 + t260 * t378;
t537 = t203 / 0.2e1;
t535 = t423 / 0.2e1;
t533 = t245 / 0.2e1;
t532 = t259 / 0.2e1;
t531 = t260 / 0.2e1;
t530 = t265 / 0.2e1;
t529 = t266 / 0.2e1;
t526 = -t301 / 0.2e1;
t524 = t327 / 0.2e1;
t523 = t328 / 0.2e1;
t522 = t336 / 0.2e1;
t520 = t344 / 0.2e1;
t518 = -t372 / 0.2e1;
t516 = -t377 / 0.2e1;
t511 = t3 * t372;
t510 = t81 * Ifges(6,1);
t508 = t82 * Ifges(6,4);
t504 = Ifges(3,4) * t376;
t496 = Ifges(3,5) * t381;
t495 = Ifges(6,5) * t336;
t489 = t579 * Ifges(6,6);
t488 = t578 * Ifges(6,5);
t487 = t423 * Ifges(5,6);
t486 = t245 * Ifges(5,5);
t482 = t300 * Ifges(4,6);
t480 = t301 * Ifges(4,5);
t479 = t309 * mrSges(3,2);
t478 = t336 * Ifges(6,3);
t477 = t344 * Ifges(5,3);
t476 = t350 * Ifges(4,3);
t475 = t360 * Ifges(3,5);
t468 = t372 * t185;
t467 = t373 * t374;
t465 = t377 * t185;
t460 = -mrSges(3,1) * t360 - mrSges(4,1) * t300 + mrSges(4,2) * t301 + mrSges(3,3) * t435;
t459 = t274 - t287;
t322 = t334 * qJD(2);
t326 = pkin(3) * t466 + t373 * t367;
t437 = Ifges(4,5) * t265 + Ifges(4,6) * t266 + Ifges(4,3) * t421;
t213 = -t226 * t372 + t377 * t435;
t426 = t213 + t468;
t214 = t226 * t377 + t372 * t435;
t425 = -t214 + t465;
t249 = -pkin(3) * t279 + t322;
t212 = pkin(3) * t301 + t583;
t23 = mrSges(7,1) * t82 - mrSges(7,3) * t49;
t24 = -mrSges(7,2) * t82 + mrSges(7,3) * t50;
t410 = -t372 * t23 + t377 * t24;
t312 = t361 + (-pkin(2) - t513) * t371;
t264 = -t327 * pkin(3) + t312;
t209 = -t259 * pkin(4) + t264;
t105 = -pkin(5) * t401 - t203 * pkin(12) + t209;
t77 = -pkin(12) * t469 + t565;
t39 = t105 * t377 - t372 * t77;
t40 = t105 * t372 + t377 * t77;
t83 = t378 * t124 - t373 * t126;
t93 = t135 * t378 - t141 * t373;
t405 = t174 * t380 - t175 * t375;
t325 = -pkin(3) * t467 + t367 * t378;
t188 = -t372 * t203 - t377 * t469;
t397 = -t377 * t203 + t372 * t469;
t109 = -mrSges(7,2) * t173 + mrSges(7,3) * t159;
t110 = mrSges(7,1) * t173 - mrSges(7,3) * t160;
t162 = -mrSges(6,2) * t336 + mrSges(6,3) * t579;
t396 = t109 * t377 - t110 * t372 + t162;
t19 = t124 * t448 - t126 * t449 + t373 * t52 + t378 * t56;
t138 = -pkin(4) * t183 + t249;
t391 = -qJD(6) * t409 - t511;
t387 = -t110 * t446 - t109 * t447 + m(7) * (-t25 * t446 - t26 * t447 - t511 + t512) + t410;
t114 = t172 + t495 + t505;
t383 = t63 * mrSges(6,3) - t114 / 0.2e1 - t495 / 0.2e1 - t394 / 0.2e1 - t393 / 0.2e1 - t392 / 0.2e1 - t187 * mrSges(6,2) + t88 * t517 + t89 * t516 - t395;
t382 = t383 + t569;
t354 = Ifges(3,4) * t434;
t349 = t428 * t496;
t324 = pkin(12) + t326;
t323 = -pkin(5) - t325;
t316 = -t360 * mrSges(3,2) + mrSges(3,3) * t434;
t282 = Ifges(3,1) * t435 + t354 + t475;
t281 = Ifges(3,6) * t360 + (Ifges(3,2) * t381 + t504) * t457;
t270 = mrSges(4,1) * t350 - mrSges(4,3) * t301;
t269 = -mrSges(4,2) * t350 + mrSges(4,3) * t300;
t267 = t367 * t448 + (-t374 * t449 + (t378 * t379 - t467) * qJD(4)) * pkin(3);
t248 = -mrSges(4,2) * t421 + mrSges(4,3) * t266;
t247 = mrSges(4,1) * t421 - mrSges(4,3) * t265;
t229 = t476 + t480 + t482;
t222 = mrSges(5,1) * t344 - mrSges(5,3) * t245;
t221 = -mrSges(5,2) * t344 + mrSges(5,3) * t423;
t210 = -mrSges(4,1) * t266 + mrSges(4,2) * t265;
t201 = Ifges(4,1) * t265 + Ifges(4,4) * t266 + Ifges(4,5) * t421;
t200 = t265 * Ifges(4,4) + t266 * Ifges(4,2) + Ifges(4,6) * t421;
t184 = -mrSges(5,1) * t423 + mrSges(5,2) * t245;
t168 = t477 + t486 + t487;
t145 = -mrSges(5,2) * t421 + mrSges(5,3) * t156;
t144 = mrSges(5,1) * t421 - mrSges(5,3) * t155;
t119 = -mrSges(6,1) * t579 + mrSges(6,2) * t578;
t112 = t478 + t488 + t489;
t103 = -mrSges(5,1) * t156 + mrSges(5,2) * t155;
t102 = t120 + t583;
t101 = qJD(5) * t199 + t196 * t373 - t378 * t390;
t99 = t155 * Ifges(5,1) + t156 * Ifges(5,4) + Ifges(5,5) * t421;
t98 = t155 * Ifges(5,4) + t156 * Ifges(5,2) + Ifges(5,6) * t421;
t97 = t120 + t212;
t96 = qJD(5) * t203 + t182 * t373 - t378 * t183;
t95 = qJD(5) * t401 + t182 * t378 + t183 * t373;
t90 = -pkin(5) * t435 - t93;
t76 = pkin(5) * t469 - t83;
t74 = -mrSges(6,2) * t421 - mrSges(6,3) * t82;
t73 = mrSges(6,1) * t421 - mrSges(6,3) * t81;
t70 = qJD(6) * t397 - t372 * t95 + t377 * t433;
t69 = qJD(6) * t188 + t372 * t433 + t377 * t95;
t68 = t378 * t118 + t373 * t399;
t66 = t115 * t378 - t473;
t65 = t115 * t373 + t464;
t38 = t120 * t372 + t377 * t63;
t37 = t120 * t377 - t372 * t63;
t34 = t102 * t372 + t377 * t66;
t33 = t102 * t377 - t372 * t66;
t32 = t372 * t97 + t377 * t68;
t31 = -t372 * t68 + t377 * t97;
t30 = mrSges(6,1) * t82 + mrSges(6,2) * t81;
t29 = Ifges(6,5) * t421 - t508 + t510;
t28 = Ifges(6,6) * t421 + t509 - t573;
t27 = pkin(5) * t96 - pkin(12) * t95 + t138;
t18 = -pkin(5) * t433 - t20;
t17 = pkin(12) * t433 + t19;
t5 = -qJD(6) * t40 - t17 * t372 + t27 * t377;
t4 = qJD(6) * t39 + t17 * t377 + t27 * t372;
t1 = [(-Ifges(6,4) * t539 + Ifges(7,5) * t543 - Ifges(6,2) * t540 - Ifges(6,6) * t522 + Ifges(7,6) * t545 + Ifges(7,3) * t541 - t559) * t96 + (Ifges(7,1) * t69 + Ifges(7,4) * t70) * t543 + t81 * (Ifges(6,1) * t203 - Ifges(6,5) * t469) / 0.2e1 + (Ifges(6,1) * t95 + Ifges(6,5) * t433) * t539 - (t437 + t438 + t442) * t469 / 0.2e1 + (t381 * t282 + t360 * (-Ifges(3,6) * t376 + t496) + (t229 + t168 + t112) * t376) * t455 / 0.2e1 + m(6) * (t10 * t565 + t11 * t83 + t123 * t209 + t138 * t187 + t19 * t64 + t20 * t63) + t565 * t74 + (-Ifges(7,1) * t397 + Ifges(7,4) * t188) * t555 + (-Ifges(7,5) * t397 + Ifges(7,6) * t188) * t553 + (-Ifges(7,4) * t397 + Ifges(7,2) * t188) * t554 + (t188 * t2 - t25 * t69 + t26 * t70 + t3 * t397) * mrSges(7,3) + t8 * (-mrSges(7,1) * t188 - mrSges(7,2) * t397) - t397 * t16 / 0.2e1 + (Ifges(7,5) * t69 + Ifges(7,6) * t70) * t541 + t11 * (-mrSges(6,1) * t469 - t203 * mrSges(6,3)) + t62 * (-mrSges(5,1) * t469 - t260 * mrSges(5,3)) + t175 * (-mrSges(4,1) * t469 - t328 * mrSges(4,3)) + t61 * (mrSges(5,2) * t469 + t259 * mrSges(5,3)) + t174 * (mrSges(4,2) * t469 + t327 * mrSges(4,3)) + t183 * t576 + (t309 * t381 + t310 * t376 + (-t317 * t381 - t320 * t376) * qJD(2)) * t370 * mrSges(3,3) + t460 * t322 + ((-t333 * mrSges(3,3) + Ifges(3,5) * t371 / 0.2e1 + (-0.2e1 * t549 + 0.3e1 / 0.2e1 * Ifges(3,4) * t381) * t370) * t381 + (-t334 * mrSges(3,3) - Ifges(3,6) * t371 + Ifges(6,5) * t537 + Ifges(6,6) * t538 + Ifges(5,5) * t531 + Ifges(5,6) * t532 + Ifges(4,5) * t523 + Ifges(4,6) * t524 + (-0.2e1 * t550 - 0.3e1 / 0.2e1 * t504) * t370 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(6,3) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(4,3) / 0.2e1) * t469) * t376) * t428 - (t573 / 0.2e1 + Ifges(7,3) * t553 + Ifges(7,6) * t554 + Ifges(7,5) * t555 + t581) * t401 + (Ifges(6,5) * t95 + Ifges(6,3) * t433) * t522 + m(7) * (t18 * t58 + t2 * t40 + t25 * t5 + t26 * t4 + t3 * t39 + t76 * t8) + m(5) * (t130 * t72 + t131 * t71 + t146 * t62 + t147 * t61 + t232 * t264 + t246 * t249) + m(4) * (t174 * t252 + t175 * t251 + t194 * t234 + t195 * t233 + t284 * t322 + t310 * t312) + m(3) * (t309 * t334 - t310 * t333 - t317 * t322 + t320 * t321) + (t349 / 0.2e1 - t310 * mrSges(3,1) - t479) * t371 + (t10 * t469 + t123 * t203 + t187 * t95 - t433 * t64) * mrSges(6,2) + (Ifges(7,4) * t69 + Ifges(7,2) * t70) * t545 - t281 * t433 / 0.2e1 + t131 * (-mrSges(5,2) * t433 + mrSges(5,3) * t183) + t234 * (-mrSges(4,2) * t433 + mrSges(4,3) * t279) + t63 * (mrSges(6,1) * t433 - mrSges(6,3) * t95) + t130 * (mrSges(5,1) * t433 - mrSges(5,3) * t182) + t233 * (mrSges(4,1) * t433 - mrSges(4,3) * t278) - t82 * (Ifges(6,4) * t203 - Ifges(6,6) * t469) / 0.2e1 + (Ifges(6,4) * t95 + Ifges(6,6) * t433) * t540 + t39 * t23 + t40 * t24 + (Ifges(4,5) * t278 + Ifges(4,6) * t279 + Ifges(4,3) * t433) * t519 + (Ifges(5,5) * t182 + Ifges(5,6) * t183 + Ifges(5,3) * t433) * t520 + t201 * t523 + t200 * t524 + (Ifges(4,1) * t278 + Ifges(4,4) * t279 + Ifges(4,5) * t433) * t525 + (Ifges(4,4) * t278 + Ifges(4,2) * t279 + Ifges(4,6) * t433) * t527 + t58 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t76 * t21 + t83 * t73 + t70 * t88 / 0.2e1 + t18 * t104 + (Ifges(4,4) * t328 + Ifges(4,2) * t327 - Ifges(4,6) * t469) * t529 + (Ifges(4,1) * t328 + Ifges(4,4) * t327 - Ifges(4,5) * t469) * t530 + t99 * t531 + t98 * t532 + (Ifges(5,1) * t182 + Ifges(5,4) * t183 + Ifges(5,5) * t433) * t533 + (Ifges(5,4) * t182 + Ifges(5,2) * t183 + Ifges(5,6) * t433) * t535 + t29 * t537 + t28 * t538 + (Ifges(5,4) * t260 + Ifges(5,2) * t259 - Ifges(5,6) * t469) * t547 + (Ifges(5,1) * t260 + Ifges(5,4) * t259 - Ifges(5,5) * t469) * t548 + t69 * t552 + t4 * t109 + t5 * t110 + t95 * t114 / 0.2e1 + t138 * t119 + t146 * t144 + t147 * t145 + t19 * t162 + t20 * t163 + t182 * t170 / 0.2e1 + t188 * t15 / 0.2e1 + t209 * t30 + t71 * t221 + t72 * t222 + t246 * (-mrSges(5,1) * t183 + mrSges(5,2) * t182) + t249 * t184 + t251 * t247 + t252 * t248 + t232 * (-mrSges(5,1) * t259 + mrSges(5,2) * t260) + t264 * t103 + t194 * t269 + t195 * t270 + t278 * t231 / 0.2e1 + t279 * t230 / 0.2e1 + t284 * (-mrSges(4,1) * t279 + mrSges(4,2) * t278) + t312 * t210 + t321 * t316 + t310 * (-mrSges(4,1) * t327 + mrSges(4,2) * t328); (t375 * pkin(3) * t184 + (-t269 * t375 - t270 * t380) * pkin(9) + t561) * qJD(3) + (Ifges(5,5) * t287 - Ifges(5,6) * t286) * t521 + (Ifges(5,1) * t287 - Ifges(5,4) * t286) * t534 + (Ifges(5,4) * t287 - Ifges(5,2) * t286) * t536 - (t21 - t73) * t402 + ((t317 * mrSges(3,3) + t457 * t549 - t354 / 0.2e1 - t475 / 0.2e1 - t282 / 0.2e1 - t561) * t381 + (t281 / 0.2e1 - t112 / 0.2e1 - t488 / 0.2e1 - t487 / 0.2e1 - t486 / 0.2e1 - t482 / 0.2e1 - t489 / 0.2e1 + (t504 / 0.2e1 + t550 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t381) * t457 - t168 / 0.2e1 - t476 / 0.2e1 - t477 / 0.2e1 - t229 / 0.2e1 + (-qJD(2) + t360 / 0.2e1) * Ifges(3,6) + t234 * mrSges(4,2) - t233 * mrSges(4,1) + t320 * mrSges(3,3) - t63 * mrSges(6,1) + t64 * mrSges(6,2) - t130 * mrSges(5,1) + t131 * mrSges(5,2) - t480 / 0.2e1 - t478 / 0.2e1 + (Ifges(4,5) * t375 + Ifges(5,5) * t338 + Ifges(6,5) * t273 + Ifges(4,6) * t380 + Ifges(5,6) * t337 + Ifges(6,6) * t400) * qJD(2) / 0.2e1) * t376) * t457 + (-t275 / 0.2e1 + t286 / 0.2e1) * t169 + (Ifges(5,5) * t274 - Ifges(5,6) * t275) * t520 + (-t213 / 0.2e1 - t468 / 0.2e1) * t88 + (t10 * t199 + t11 * t402 + t123 * t311 + t588 * t64 + (-t101 - t93) * t63 + t585 * t187) * m(6) + t588 * t162 + (-t247 * t375 + t248 * t380) * pkin(9) + (t87 - t113) * (t186 / 0.2e1 - t225 / 0.2e1) - t578 * (Ifges(6,1) * t226 - Ifges(6,4) * t225) / 0.2e1 + (-t11 * mrSges(6,3) + t29 / 0.2e1 + t510 / 0.2e1 - t508 / 0.2e1 + t123 * mrSges(6,2) + t8 * t417 + t416 * t555 + t414 * t554 + t412 * t553 + t15 * t518 + t16 * t515 + (-t2 * t372 - t3 * t377) * mrSges(7,3) + (-mrSges(7,3) * t563 + t411 * t542 + t413 * t546 + t415 * t544 + t58 * t418 + t88 * t516 + t89 * t518) * qJD(6)) * t273 + (-t214 / 0.2e1 + t465 / 0.2e1) * t89 + t586 * t221 + t587 * t222 - t336 * (Ifges(6,5) * t226 - Ifges(6,6) * t225) / 0.2e1 - t460 * t320 + t405 * mrSges(4,3) + (Ifges(5,1) * t274 - Ifges(5,4) * t275) * t533 + (Ifges(5,4) * t274 - Ifges(5,2) * t275) * t535 - (t48 / 0.2e1 + t47 / 0.2e1 - t28 / 0.2e1 + (Ifges(7,3) / 0.2e1 + t557) * t82 + t581) * t400 + t585 * t119 + t200 * t514 + t349 + (t587 * t130 + t586 * t131 + t232 * t368 + t246 * t594 + t288 * t62 + t289 * t61) * m(5) - t579 * (Ifges(6,4) * t226 - Ifges(6,2) * t225) / 0.2e1 + t590 * t109 + t591 * t110 + (t127 * t3 + t128 * t2 - t402 * t8 + (t101 - t90) * t58 + t590 * t26 + t591 * t25) * m(7) + (mrSges(7,1) * t461 - mrSges(7,3) * t425) * t25 + (-mrSges(7,2) * t461 - mrSges(7,3) * t426) * t26 + (mrSges(6,1) * t461 + mrSges(6,2) * t462) * t187 + (-t461 * t64 - t462 * t63) * mrSges(6,3) + t463 * t101 + (mrSges(5,1) * t595 + mrSges(5,2) * t459) * t246 + (-t130 * t459 - t131 * t595 + t337 * t61 - t338 * t62) * mrSges(5,3) + (-mrSges(4,1) * t380 + mrSges(4,2) * t375 - mrSges(3,1)) * t310 - t479 + (t185 / 0.2e1 - t226 / 0.2e1) * t114 + (Ifges(6,5) * t185 - Ifges(6,6) * t186) * t522 + (mrSges(7,1) * t426 + mrSges(7,2) * t425) * t58 - t90 * t104 + (Ifges(4,2) * t380 + t503) * t529 + (Ifges(4,1) * t375 + t502) * t530 + (Ifges(6,1) * t185 - Ifges(6,4) * t186) * t539 + (Ifges(6,4) * t185 - Ifges(6,2) * t186) * t540 + (Ifges(7,5) * t465 - Ifges(7,6) * t468 + Ifges(7,3) * t186) * t541 + (Ifges(7,5) * t214 + Ifges(7,6) * t213 + Ifges(7,3) * t225) * t542 + (Ifges(7,1) * t465 - Ifges(7,4) * t468 + Ifges(7,5) * t186) * t543 + (Ifges(7,1) * t214 + Ifges(7,4) * t213 + Ifges(7,5) * t225) * t544 + (Ifges(7,4) * t465 - Ifges(7,2) * t468 + Ifges(7,6) * t186) * t545 + (Ifges(7,4) * t214 + Ifges(7,2) * t213 + Ifges(7,6) * t225) * t546 + (Ifges(5,4) * t338 + Ifges(5,2) * t337) * t547 + (Ifges(5,1) * t338 + Ifges(5,4) * t337) * t548 + (-pkin(2) * t310 - t233 * t253 - t234 * t254 - t284 * t320 + (-qJD(3) * t404 + t405) * pkin(9)) * m(4) + t127 * t23 + t128 * t24 - t93 * t163 + t199 * t74 - pkin(2) * t210 + (t274 / 0.2e1 - t287 / 0.2e1) * t170 - t254 * t269 - t253 * t270 - t277 * t184 + t288 * t144 + t289 * t145 + t311 * t30 - t317 * t316 + t337 * t98 / 0.2e1 + t338 * t99 / 0.2e1 + t232 * (-mrSges(5,1) * t337 + mrSges(5,2) * t338) + t368 * t103 + t375 * t201 / 0.2e1; -(-Ifges(4,2) * t301 + t231 + t296) * t300 / 0.2e1 + ((-t109 * t372 - t110 * t377) * t324 - t569) * qJD(6) + t245 * t576 + (t382 + t564) * t579 - mrSges(7,3) * t511 + t437 + t410 * t324 + (t233 * t300 + t234 * t301) * mrSges(4,3) + t584 + t396 * t267 + (t144 * t379 + t145 * t374 - t184 * t301 + (t221 * t379 - t222 * t374) * qJD(4) + (-t130 * t451 + t131 * t450 + 0.2e1 * t246 * t526 + t374 * t61 + t379 * t62) * m(5)) * pkin(3) - m(5) * (t130 * t132 + t131 * t133) + t230 * t525 + (Ifges(4,1) * t300 - t481) * t526 - t463 * t568 + (t10 * t326 + t11 * t325 - t187 * t212 + (t267 - t68) * t64 + t568 * t63) * m(6) + (t323 * t8 + (t391 + t512) * t324 - t25 * t31 - t26 * t32 - t568 * t58 + t563 * t267) * m(7) - t32 * t109 - t31 * t110 - t68 * t162 - t174 * mrSges(4,2) + t175 * mrSges(4,1) - t212 * t119 - t133 * t221 - t132 * t222 - t233 * t269 + t234 * t270 - t284 * (mrSges(4,1) * t301 + mrSges(4,2) * t300) + t323 * t21 + t325 * t73 + t326 * t74 - t350 * (Ifges(4,5) * t300 - Ifges(4,6) * t301) / 0.2e1; (-t245 * t119 + t373 * t74 + t378 * t73 + ((m(7) * t58 + t463) * t373 + (m(7) * t563 + t396) * t378) * qJD(5) + (0.2e1 * t187 * t534 + t10 * t373 + t11 * t378 + (-t373 * t63 + t378 * t64) * qJD(5)) * m(6)) * pkin(4) - m(6) * (-t63 * t65 + t64 * t66) + (t383 + t564) * t579 + (-t173 * t484 + (-t173 * t26 - t3) * t372) * mrSges(7,3) - m(7) * (t25 * t33 + t26 * t34 + t58 * t65) - t463 * t65 + t387 * (pkin(4) * t373 + pkin(12)) + t169 * t533 - t34 * t109 - t33 * t110 - t66 * t162 - t130 * t221 + t131 * t222 + t589 * (-pkin(4) * t378 - pkin(5)) + t584; t391 * mrSges(7,3) - m(7) * (t25 * t37 + t26 * t38 + t58 * t64) - t463 * t64 + t386 * t578 + t387 * pkin(12) + (t439 + t382 + (t557 - Ifges(6,1) / 0.2e1) * t578) * t579 + t389 - t38 * t109 - t37 * t110 - t63 * t162 - t589 * pkin(5); -t58 * (mrSges(7,1) * t160 + mrSges(7,2) * t159) - t25 * t109 + t26 * t110 + (Ifges(7,1) * t159 - t499) * t544 + t88 * t543 + (Ifges(7,5) * t159 - Ifges(7,6) * t160) * t542 + (t159 * t25 + t160 * t26) * mrSges(7,3) + t419 + t14 + (-Ifges(7,2) * t160 + t158 + t89) * t546;];
tauc  = t1(:);
