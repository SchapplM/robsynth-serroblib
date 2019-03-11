% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRR7
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:29:03
% EndTime: 2019-03-10 04:30:15
% DurationCPUTime: 38.07s
% Computational Cost: add. (42082->1045), mult. (105442->1464), div. (0->0), fcn. (84633->12), ass. (0->435)
t417 = sin(qJ(4));
t419 = sin(qJ(2));
t422 = cos(qJ(4));
t413 = sin(pkin(6));
t474 = qJD(1) * t413;
t423 = cos(qJ(3));
t424 = cos(qJ(2));
t487 = t423 * t424;
t326 = (-t417 * t487 + t419 * t422) * t474;
t418 = sin(qJ(3));
t468 = qJD(4) * t422;
t470 = qJD(3) * t423;
t645 = t417 * t470 + t418 * t468 + t326;
t460 = t419 * t474;
t414 = cos(pkin(6));
t537 = pkin(1) * t424;
t463 = t414 * t537;
t357 = -pkin(8) * t460 + qJD(1) * t463;
t436 = t413 * (pkin(2) * t419 - pkin(9) * t424);
t358 = qJD(1) * t436;
t273 = t423 * t357 + t418 * t358;
t257 = pkin(10) * t460 + t273;
t405 = t414 * t419 * pkin(1);
t448 = pkin(3) * t418 - pkin(10) * t423;
t492 = t413 * t424;
t276 = (t405 + (pkin(8) + t448) * t492) * qJD(1);
t183 = t422 * t257 + t417 * t276;
t387 = t448 * qJD(3);
t392 = -pkin(3) * t423 - pkin(10) * t418 - pkin(2);
t469 = qJD(4) * t417;
t471 = qJD(3) * t418;
t263 = t417 * t387 + t392 * t468 + (-t422 * t471 - t423 * t469) * pkin(9);
t644 = -t183 + t263;
t182 = -t417 * t257 + t422 * t276;
t327 = (t417 * t419 + t422 * t487) * t474;
t488 = t422 * t423;
t407 = pkin(9) * t488;
t459 = t424 * t474;
t450 = t418 * t459;
t536 = pkin(9) * t417;
t476 = t422 * t387 + t471 * t536;
t643 = -pkin(4) * t450 + t327 * pkin(11) - t182 + (pkin(4) * t418 - pkin(11) * t488) * qJD(3) + (-t407 + (pkin(11) * t418 - t392) * t417) * qJD(4) + t476;
t642 = pkin(11) * t645 - t644;
t475 = pkin(8) * t492 + t405;
t360 = t475 * qJD(1);
t401 = qJD(1) * t414 + qJD(2);
t317 = t401 * pkin(9) + t360;
t354 = (-pkin(2) * t424 - pkin(9) * t419 - pkin(1)) * t413;
t331 = qJD(1) * t354;
t245 = -t418 * t317 + t331 * t423;
t393 = qJD(3) - t459;
t224 = -pkin(3) * t393 - t245;
t339 = t401 * t418 + t423 * t460;
t279 = -t339 * t417 + t393 * t422;
t171 = -pkin(4) * t279 + t224;
t280 = t339 * t422 + t393 * t417;
t416 = sin(qJ(5));
t421 = cos(qJ(5));
t451 = t421 * t279 - t280 * t416;
t112 = -pkin(5) * t451 + t171;
t207 = t279 * t416 + t280 * t421;
t415 = sin(qJ(6));
t420 = cos(qJ(6));
t609 = -t207 * t415 + t420 * t451;
t121 = Ifges(7,4) * t609;
t127 = t207 * t420 + t415 * t451;
t338 = t401 * t423 - t418 * t460;
t333 = qJD(4) - t338;
t322 = qJD(5) + t333;
t316 = -t401 * pkin(2) - t357;
t220 = -t338 * pkin(3) - t339 * pkin(10) + t316;
t246 = t423 * t317 + t418 * t331;
t225 = pkin(10) * t393 + t246;
t141 = t422 * t220 - t225 * t417;
t116 = -pkin(11) * t280 + t141;
t103 = pkin(4) * t333 + t116;
t142 = t220 * t417 + t225 * t422;
t117 = pkin(11) * t279 + t142;
t113 = t416 * t117;
t59 = t421 * t103 - t113;
t621 = pkin(12) * t207;
t49 = t59 - t621;
t45 = pkin(5) * t322 + t49;
t115 = t421 * t117;
t60 = t103 * t416 + t115;
t604 = pkin(12) * t451;
t50 = t60 + t604;
t500 = t415 * t50;
t19 = t420 * t45 - t500;
t498 = t420 * t50;
t20 = t415 * t45 + t498;
t526 = Ifges(7,4) * t127;
t318 = qJD(6) + t322;
t551 = -t318 / 0.2e1;
t565 = t127 / 0.2e1;
t566 = -t127 / 0.2e1;
t568 = -t609 / 0.2e1;
t472 = qJD(2) * t424;
t457 = t423 * t472;
t288 = t401 * t470 + (-t419 * t471 + t457) * t474;
t473 = qJD(2) * t413;
t453 = qJD(1) * t473;
t449 = t419 * t453;
t189 = qJD(4) * t279 + t288 * t422 + t417 * t449;
t458 = t418 * t472;
t289 = t401 * t471 + (t419 * t470 + t458) * t474;
t359 = qJD(2) * t436;
t348 = qJD(1) * t359;
t493 = t413 * t419;
t402 = pkin(8) * t493;
t374 = -t402 + t463;
t361 = t374 * qJD(2);
t349 = qJD(1) * t361;
t172 = -t317 * t471 + t331 * t470 + t418 * t348 + t423 * t349;
t161 = pkin(10) * t449 + t172;
t362 = t475 * qJD(2);
t350 = qJD(1) * t362;
t186 = t289 * pkin(3) - t288 * pkin(10) + t350;
t73 = -qJD(4) * t142 - t161 * t417 + t422 * t186;
t48 = pkin(4) * t289 - pkin(11) * t189 + t73;
t190 = -qJD(4) * t280 - t288 * t417 + t422 * t449;
t72 = t422 * t161 + t417 * t186 + t220 * t468 - t225 * t469;
t55 = pkin(11) * t190 + t72;
t14 = -qJD(5) * t60 - t416 * t55 + t421 * t48;
t88 = qJD(5) * t451 + t189 * t421 + t190 * t416;
t6 = pkin(5) * t289 - pkin(12) * t88 + t14;
t466 = qJD(5) * t421;
t467 = qJD(5) * t416;
t13 = t103 * t466 - t117 * t467 + t416 * t48 + t421 * t55;
t89 = -qJD(5) * t207 - t189 * t416 + t190 * t421;
t7 = pkin(12) * t89 + t13;
t2 = qJD(6) * t19 + t415 * t6 + t420 * t7;
t3 = -qJD(6) * t20 - t415 * t7 + t420 * t6;
t603 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t67 = Ifges(7,2) * t609 + Ifges(7,6) * t318 + t526;
t68 = Ifges(7,1) * t127 + Ifges(7,5) * t318 + t121;
t34 = -qJD(6) * t127 - t415 * t88 + t420 * t89;
t31 = Ifges(7,6) * t34;
t33 = qJD(6) * t609 + t415 * t89 + t420 * t88;
t32 = Ifges(7,5) * t33;
t8 = Ifges(7,3) * t289 + t31 + t32;
t641 = (Ifges(7,1) * t609 - t526) * t566 + (Ifges(7,5) * t609 - Ifges(7,6) * t127) * t551 + (t127 * t20 + t19 * t609) * mrSges(7,3) - t112 * (mrSges(7,1) * t127 + mrSges(7,2) * t609) + t67 * t565 + t603 + t8 + (-Ifges(7,2) * t127 + t121 + t68) * t568;
t265 = pkin(3) * t339 - pkin(10) * t338;
t169 = -t245 * t417 + t422 * t265;
t572 = -pkin(11) - pkin(10);
t461 = qJD(4) * t572;
t534 = pkin(11) * t422;
t640 = -pkin(4) * t339 + t338 * t534 + t422 * t461 - t169;
t170 = t422 * t245 + t417 * t265;
t494 = t338 * t417;
t639 = -pkin(11) * t494 - t417 * t461 + t170;
t381 = t416 * t422 + t417 * t421;
t584 = qJD(4) + qJD(5);
t303 = t584 * t381;
t438 = t416 * t417 - t421 * t422;
t235 = -t303 * t418 - t438 * t470;
t251 = t326 * t416 + t327 * t421;
t481 = t235 - t251;
t364 = t438 * t418;
t236 = t364 * t584 - t381 * t470;
t250 = t326 * t421 - t327 * t416;
t480 = t236 - t250;
t201 = Ifges(6,4) * t451;
t111 = Ifges(6,1) * t207 + Ifges(6,5) * t322 + t201;
t86 = Ifges(6,6) * t89;
t87 = Ifges(6,5) * t88;
t39 = Ifges(6,3) * t289 + t86 + t87;
t527 = Ifges(6,4) * t207;
t549 = -t322 / 0.2e1;
t559 = -t207 / 0.2e1;
t561 = -t451 / 0.2e1;
t592 = t14 * mrSges(6,1) - t13 * mrSges(6,2);
t638 = t39 + t592 + (Ifges(6,5) * t451 - Ifges(6,6) * t207) * t549 + (t207 * t60 + t451 * t59) * mrSges(6,3) + (-Ifges(6,2) * t207 + t111 + t201) * t561 - t171 * (mrSges(6,1) * t207 + mrSges(6,2) * t451) + (Ifges(6,1) * t451 - t527) * t559 + t641;
t379 = t422 * t392;
t298 = -t418 * t534 + t379 + (-pkin(4) - t536) * t423;
t347 = t417 * t392 + t407;
t489 = t417 * t418;
t312 = -pkin(11) * t489 + t347;
t618 = t298 * t466 - t312 * t467 + t416 * t643 - t642 * t421;
t227 = t416 * t298 + t421 * t312;
t617 = -qJD(5) * t227 + t642 * t416 + t421 * t643;
t248 = t381 * t338;
t636 = t248 - t303;
t249 = t438 * t338;
t302 = t584 * t438;
t478 = -t249 + t302;
t634 = t617 - t481 * pkin(12) + (-t450 + t471) * pkin(5);
t633 = -pkin(12) * t480 - t618;
t394 = t572 * t417;
t395 = t572 * t422;
t321 = t416 * t394 - t421 * t395;
t591 = -qJD(5) * t321 + t639 * t416 + t421 * t640;
t590 = t394 * t466 + t395 * t467 + t416 * t640 - t639 * t421;
t626 = pkin(12) * t636 + t590;
t625 = -pkin(5) * t339 + pkin(12) * t478 + t591;
t271 = -t418 * t357 + t358 * t423;
t256 = -pkin(3) * t460 - t271;
t613 = pkin(4) * t645 + pkin(9) * t470 - t256;
t517 = t280 * Ifges(5,4);
t175 = t279 * Ifges(5,2) + t333 * Ifges(5,6) + t517;
t275 = Ifges(5,4) * t279;
t176 = t280 * Ifges(5,1) + t333 * Ifges(5,5) + t275;
t439 = t141 * t422 + t142 * t417;
t441 = Ifges(5,5) * t422 - Ifges(5,6) * t417;
t528 = Ifges(5,4) * t422;
t443 = -Ifges(5,2) * t417 + t528;
t529 = Ifges(5,4) * t417;
t445 = Ifges(5,1) * t422 - t529;
t446 = mrSges(5,1) * t417 + mrSges(5,2) * t422;
t539 = t422 / 0.2e1;
t540 = -t417 / 0.2e1;
t545 = t333 / 0.2e1;
t554 = t280 / 0.2e1;
t556 = t279 / 0.2e1;
t624 = t439 * mrSges(5,3) - t175 * t540 - t176 * t539 - t224 * t446 - t441 * t545 - t443 * t556 - t445 * t554;
t226 = t421 * t298 - t312 * t416;
t195 = -pkin(5) * t423 + pkin(12) * t364 + t226;
t363 = t381 * t418;
t202 = -pkin(12) * t363 + t227;
t119 = t195 * t420 - t202 * t415;
t620 = qJD(6) * t119 + t415 * t634 - t420 * t633;
t120 = t195 * t415 + t202 * t420;
t619 = -qJD(6) * t120 + t415 * t633 + t420 * t634;
t409 = pkin(4) * t421 + pkin(5);
t464 = qJD(6) * t420;
t465 = qJD(6) * t415;
t490 = t416 * t420;
t63 = -t116 * t416 - t115;
t52 = t63 - t604;
t64 = t421 * t116 - t113;
t53 = t64 - t621;
t616 = t415 * t53 - t420 * t52 - t409 * t465 + (-t416 * t464 + (-t415 * t421 - t490) * qJD(5)) * pkin(4);
t491 = t415 * t416;
t615 = -t415 * t52 - t420 * t53 + t409 * t464 + (-t416 * t465 + (t420 * t421 - t491) * qJD(5)) * pkin(4);
t614 = -pkin(5) * t480 + t613;
t332 = Ifges(4,4) * t338;
t504 = t393 * Ifges(4,5);
t507 = t339 * Ifges(4,1);
t242 = t332 + t504 + t507;
t431 = t245 * mrSges(4,3) - t242 / 0.2e1 - t316 * mrSges(4,2) - t504 / 0.2e1;
t612 = t431 + t624;
t581 = t33 / 0.2e1;
t580 = t34 / 0.2e1;
t575 = t88 / 0.2e1;
t574 = t89 / 0.2e1;
t110 = Ifges(6,2) * t451 + Ifges(6,6) * t322 + t527;
t606 = t110 / 0.2e1;
t553 = t289 / 0.2e1;
t605 = -Ifges(3,6) * t401 / 0.2e1;
t320 = t421 * t394 + t395 * t416;
t277 = -pkin(12) * t381 + t320;
t278 = -pkin(12) * t438 + t321;
t197 = t277 * t420 - t278 * t415;
t602 = qJD(6) * t197 + t415 * t625 + t420 * t626;
t198 = t277 * t415 + t278 * t420;
t601 = -qJD(6) * t198 - t415 * t626 + t420 * t625;
t131 = -mrSges(6,1) * t451 + mrSges(6,2) * t207;
t588 = m(6) * t171 + t131;
t352 = t402 + (-pkin(2) - t537) * t414;
t367 = -t414 * t423 + t418 * t493;
t368 = t414 * t418 + t423 * t493;
t253 = t367 * pkin(3) - t368 * pkin(10) + t352;
t353 = pkin(9) * t414 + t475;
t267 = t423 * t353 + t418 * t354;
t255 = -pkin(10) * t492 + t267;
t165 = t422 * t253 - t255 * t417;
t435 = -t368 * t422 + t417 * t492;
t137 = pkin(4) * t367 + pkin(11) * t435 + t165;
t166 = t417 * t253 + t422 * t255;
t308 = -t368 * t417 - t422 * t492;
t143 = pkin(11) * t308 + t166;
t81 = t416 * t137 + t421 * t143;
t209 = pkin(4) * t494 + t246;
t587 = pkin(4) * t469 - pkin(5) * t636 - t209;
t586 = -t417 * t73 + t422 * t72;
t585 = -t73 * mrSges(5,1) + t72 * mrSges(5,2);
t583 = Ifges(7,4) * t581 + Ifges(7,2) * t580 + Ifges(7,6) * t553;
t582 = Ifges(7,1) * t581 + Ifges(7,4) * t580 + Ifges(7,5) * t553;
t579 = Ifges(6,4) * t575 + Ifges(6,2) * t574 + Ifges(6,6) * t553;
t578 = Ifges(6,1) * t575 + Ifges(6,4) * t574 + Ifges(6,5) * t553;
t98 = t189 * Ifges(5,1) + t190 * Ifges(5,4) + t289 * Ifges(5,5);
t573 = t98 / 0.2e1;
t571 = pkin(1) * mrSges(3,1);
t570 = pkin(1) * mrSges(3,2);
t567 = t609 / 0.2e1;
t564 = -t175 / 0.2e1;
t563 = t189 / 0.2e1;
t562 = t190 / 0.2e1;
t560 = t451 / 0.2e1;
t558 = t207 / 0.2e1;
t557 = -t279 / 0.2e1;
t555 = -t280 / 0.2e1;
t550 = t318 / 0.2e1;
t548 = t322 / 0.2e1;
t547 = -t332 / 0.2e1;
t546 = -t333 / 0.2e1;
t544 = -t367 / 0.2e1;
t542 = t368 / 0.2e1;
t541 = t414 / 0.2e1;
t535 = pkin(9) * t423;
t530 = Ifges(3,4) * t419;
t188 = Ifges(5,5) * t189;
t187 = Ifges(5,6) * t190;
t524 = t609 * Ifges(7,6);
t523 = t127 * Ifges(7,5);
t522 = t172 * mrSges(4,2);
t173 = -t317 * t470 - t331 * t471 + t348 * t423 - t418 * t349;
t521 = t173 * mrSges(4,1);
t520 = t451 * Ifges(6,6);
t519 = t207 * Ifges(6,5);
t518 = t279 * Ifges(5,6);
t516 = t280 * Ifges(5,5);
t515 = t288 * Ifges(4,1);
t514 = t288 * Ifges(4,4);
t513 = t289 * Ifges(4,4);
t512 = t318 * Ifges(7,3);
t511 = t322 * Ifges(6,3);
t510 = t333 * Ifges(5,3);
t509 = t338 * Ifges(4,2);
t506 = t339 * Ifges(4,4);
t503 = t393 * Ifges(4,6);
t270 = -t363 * t420 + t364 * t415;
t129 = qJD(6) * t270 + t235 * t420 + t236 * t415;
t164 = t250 * t415 + t251 * t420;
t486 = t129 - t164;
t272 = -t363 * t415 - t364 * t420;
t130 = -qJD(6) * t272 - t235 * t415 + t236 * t420;
t163 = t250 * t420 - t251 * t415;
t485 = t130 - t163;
t159 = -t248 * t420 + t249 * t415;
t297 = t381 * t420 - t415 * t438;
t179 = -qJD(6) * t297 + t302 * t415 - t303 * t420;
t484 = t159 - t179;
t160 = -t248 * t415 - t249 * t420;
t296 = -t381 * t415 - t420 * t438;
t178 = qJD(6) * t296 - t302 * t420 - t303 * t415;
t483 = t160 - t178;
t212 = -mrSges(5,1) * t279 + mrSges(5,2) * t280;
t293 = mrSges(4,1) * t393 - mrSges(4,3) * t339;
t482 = -t212 + t293;
t477 = -mrSges(3,1) * t401 - mrSges(4,1) * t338 + mrSges(4,2) * t339 + mrSges(3,3) * t460;
t388 = pkin(4) * t489 + t418 * pkin(9);
t96 = Ifges(5,3) * t289 + t187 + t188;
t462 = Ifges(4,5) * t288 - Ifges(4,6) * t289 + Ifges(4,3) * t449;
t410 = -pkin(4) * t422 - pkin(3);
t454 = t419 * t473;
t80 = t421 * t137 - t143 * t416;
t266 = -t418 * t353 + t354 * t423;
t254 = pkin(3) * t492 - t266;
t447 = mrSges(5,1) * t422 - mrSges(5,2) * t417;
t444 = Ifges(5,1) * t417 + t528;
t442 = Ifges(5,2) * t422 + t529;
t440 = Ifges(5,5) * t417 + Ifges(5,6) * t422;
t232 = t308 * t416 - t421 * t435;
t61 = pkin(5) * t367 - pkin(12) * t232 + t80;
t231 = t308 * t421 + t416 * t435;
t69 = pkin(12) * t231 + t81;
t25 = -t415 * t69 + t420 * t61;
t26 = t415 * t61 + t420 * t69;
t149 = t231 * t420 - t232 * t415;
t150 = t231 * t415 + t232 * t420;
t192 = -t353 * t470 - t354 * t471 + t359 * t423 - t418 * t361;
t307 = -qJD(3) * t367 + t413 * t457;
t219 = qJD(4) * t308 + t307 * t422 + t417 * t454;
t306 = qJD(3) * t368 + t413 * t458;
t191 = -t353 * t471 + t354 * t470 + t418 * t359 + t423 * t361;
t180 = pkin(10) * t454 + t191;
t213 = t306 * pkin(3) - t307 * pkin(10) + t362;
t85 = -qJD(4) * t166 - t180 * t417 + t422 * t213;
t58 = pkin(4) * t306 - pkin(11) * t219 + t85;
t218 = qJD(4) * t435 - t307 * t417 + t422 * t454;
t84 = t422 * t180 + t417 * t213 + t253 * t468 - t255 * t469;
t71 = pkin(11) * t218 + t84;
t17 = t137 * t466 - t143 * t467 + t416 * t58 + t421 * t71;
t208 = -pkin(4) * t308 + t254;
t396 = Ifges(3,4) * t459;
t432 = t401 * Ifges(3,5) - t357 * mrSges(3,3) + Ifges(3,1) * t460 / 0.2e1 + t396 / 0.2e1;
t181 = -pkin(3) * t454 - t192;
t18 = -qJD(5) * t81 - t416 * t71 + t421 * t58;
t162 = -pkin(3) * t449 - t173;
t128 = -pkin(4) * t218 + t181;
t104 = -pkin(4) * t190 + t162;
t429 = t245 * mrSges(4,1) + t393 * Ifges(4,3) + t339 * Ifges(4,5) + t338 * Ifges(4,6) + t605 - (Ifges(3,2) * t424 + t530) * t474 / 0.2e1 - t246 * mrSges(4,2) - t360 * mrSges(3,3);
t109 = t511 + t519 + t520;
t174 = t510 + t516 + t518;
t241 = t503 + t506 + t509;
t66 = t512 + t523 + t524;
t426 = -t316 * mrSges(4,1) - t518 / 0.2e1 - t516 / 0.2e1 - t510 / 0.2e1 + t503 / 0.2e1 - t141 * mrSges(5,1) + t142 * mrSges(5,2) - t59 * mrSges(6,1) + t60 * mrSges(6,2) + t20 * mrSges(7,2) - t19 * mrSges(7,1) - t512 / 0.2e1 - t511 / 0.2e1 + t241 / 0.2e1 - t174 / 0.2e1 - t523 / 0.2e1 - t524 / 0.2e1 - t109 / 0.2e1 - t66 / 0.2e1 - t519 / 0.2e1 - t520 / 0.2e1 + t246 * mrSges(4,3) + t506 / 0.2e1;
t425 = t426 + t509 / 0.2e1;
t391 = Ifges(3,5) * t424 * t453;
t366 = pkin(4) * t490 + t409 * t415;
t365 = -pkin(4) * t491 + t409 * t420;
t356 = -t401 * mrSges(3,2) + mrSges(3,3) * t459;
t351 = pkin(5) * t438 + t410;
t346 = -t417 * t535 + t379;
t300 = pkin(5) * t363 + t388;
t292 = -mrSges(4,2) * t393 + mrSges(4,3) * t338;
t264 = -qJD(4) * t347 + t476;
t261 = -mrSges(4,2) * t449 - mrSges(4,3) * t289;
t260 = mrSges(4,1) * t449 - mrSges(4,3) * t288;
t230 = mrSges(5,1) * t333 - mrSges(5,3) * t280;
t229 = -mrSges(5,2) * t333 + mrSges(5,3) * t279;
t216 = mrSges(4,1) * t289 + mrSges(4,2) * t288;
t200 = Ifges(4,5) * t449 - t513 + t515;
t199 = -t289 * Ifges(4,2) + Ifges(4,6) * t449 + t514;
t168 = mrSges(6,1) * t322 - mrSges(6,3) * t207;
t167 = -mrSges(6,2) * t322 + mrSges(6,3) * t451;
t156 = pkin(4) * t280 + pkin(5) * t207;
t154 = -mrSges(5,2) * t289 + mrSges(5,3) * t190;
t153 = mrSges(5,1) * t289 - mrSges(5,3) * t189;
t139 = -pkin(5) * t231 + t208;
t118 = -mrSges(5,1) * t190 + mrSges(5,2) * t189;
t108 = mrSges(7,1) * t318 - mrSges(7,3) * t127;
t107 = -mrSges(7,2) * t318 + mrSges(7,3) * t609;
t100 = -qJD(5) * t232 + t218 * t421 - t219 * t416;
t99 = qJD(5) * t231 + t218 * t416 + t219 * t421;
t97 = t189 * Ifges(5,4) + t190 * Ifges(5,2) + t289 * Ifges(5,6);
t79 = -mrSges(6,2) * t289 + mrSges(6,3) * t89;
t78 = mrSges(6,1) * t289 - mrSges(6,3) * t88;
t75 = -mrSges(7,1) * t609 + mrSges(7,2) * t127;
t62 = -pkin(5) * t100 + t128;
t51 = -pkin(5) * t89 + t104;
t44 = -mrSges(6,1) * t89 + mrSges(6,2) * t88;
t43 = -qJD(6) * t150 + t100 * t420 - t415 * t99;
t42 = qJD(6) * t149 + t100 * t415 + t420 * t99;
t28 = -mrSges(7,2) * t289 + mrSges(7,3) * t34;
t27 = mrSges(7,1) * t289 - mrSges(7,3) * t33;
t22 = t420 * t49 - t500;
t21 = -t415 * t49 - t498;
t16 = pkin(12) * t100 + t17;
t15 = pkin(5) * t306 - pkin(12) * t99 + t18;
t11 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t5 = -qJD(6) * t26 + t15 * t420 - t16 * t415;
t4 = qJD(6) * t25 + t15 * t415 + t16 * t420;
t1 = [(t432 * t424 + (t605 + t429) * t419) * t473 + (-Ifges(5,5) * t435 + Ifges(6,5) * t232 + Ifges(7,5) * t150 + Ifges(5,6) * t308 + Ifges(6,6) * t231 + Ifges(7,6) * t149 + (Ifges(5,3) + Ifges(6,3) + Ifges(7,3)) * t367) * t553 + t73 * (mrSges(5,1) * t367 + mrSges(5,3) * t435) + (-Ifges(5,4) * t435 + Ifges(5,2) * t308 + Ifges(5,6) * t367) * t562 + (-Ifges(5,1) * t435 + Ifges(5,4) * t308 + Ifges(5,5) * t367) * t563 + t162 * (-mrSges(5,1) * t308 - mrSges(5,2) * t435) - t435 * t573 + (t96 + t39 + t8) * t367 / 0.2e1 + (t174 + t109 + t66) * t306 / 0.2e1 + t339 * (Ifges(4,1) * t307 - Ifges(4,4) * t306) / 0.2e1 + t100 * t606 + t224 * (-mrSges(5,1) * t218 + mrSges(5,2) * t219) + m(6) * (t104 * t208 + t128 * t171 + t13 * t81 + t14 * t80 + t17 * t60 + t18 * t59) + m(7) * (t112 * t62 + t139 * t51 + t19 * t5 + t2 * t26 + t20 * t4 + t25 * t3) + m(5) * (t141 * t85 + t142 * t84 + t162 * t254 + t165 * t73 + t166 * t72 + t181 * t224) + m(4) * (t172 * t267 + t173 * t266 + t191 * t246 + t192 * t245 + t316 * t362 + t350 * t352) + (-t172 * t367 - t173 * t368 - t245 * t307 - t246 * t306) * mrSges(4,3) + t393 * (Ifges(4,5) * t307 - Ifges(4,6) * t306) / 0.2e1 + t338 * (Ifges(4,4) * t307 - Ifges(4,2) * t306) / 0.2e1 + t219 * t176 / 0.2e1 + t72 * (-mrSges(5,2) * t367 + mrSges(5,3) * t308) + t2 * (-mrSges(7,2) * t367 + mrSges(7,3) * t149) + t3 * (mrSges(7,1) * t367 - mrSges(7,3) * t150) + t13 * (-mrSges(6,2) * t367 + mrSges(6,3) * t231) + t14 * (mrSges(6,1) * t367 - mrSges(6,3) * t232) + t361 * t356 + t352 * t216 + t218 * t175 / 0.2e1 + t208 * t44 + t181 * t212 + t171 * (-mrSges(6,1) * t100 + mrSges(6,2) * t99) + t165 * t153 + t166 * t154 + t17 * t167 + t18 * t168 + (Ifges(6,4) * t99 + Ifges(6,2) * t100 + Ifges(6,6) * t306) * t560 + (Ifges(7,1) * t42 + Ifges(7,4) * t43 + Ifges(7,5) * t306) * t565 + t51 * (-mrSges(7,1) * t149 + mrSges(7,2) * t150) + t139 * t11 + t128 * t131 - t492 * t521 + t112 * (-mrSges(7,1) * t43 + mrSges(7,2) * t42) + t99 * t111 / 0.2e1 + t5 * t108 + t4 * t107 - t462 * t492 / 0.2e1 - t289 * (Ifges(4,4) * t368 - Ifges(4,2) * t367 - Ifges(4,6) * t492) / 0.2e1 + t288 * (Ifges(4,1) * t368 - Ifges(4,4) * t367 - Ifges(4,5) * t492) / 0.2e1 + t349 * (-t414 * mrSges(3,2) + mrSges(3,3) * t492) + (-mrSges(3,1) * t414 + mrSges(4,1) * t367 + mrSges(4,2) * t368 + mrSges(3,3) * t493) * t350 + t81 * t79 + t80 * t78 + t62 * t75 + t477 * t362 + t43 * t67 / 0.2e1 + t42 * t68 / 0.2e1 + m(3) * (t349 * t475 - t350 * t374 - t357 * t362 + t360 * t361) + ((Ifges(3,5) * t541 - t374 * mrSges(3,3) + (-0.2e1 * t570 + 0.3e1 / 0.2e1 * Ifges(3,4) * t424) * t413) * t424 + (-Ifges(3,6) * t414 + Ifges(4,5) * t542 + Ifges(4,6) * t544 - t475 * mrSges(3,3) + (-0.2e1 * t571 - 0.3e1 / 0.2e1 * t530 + (-Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t424) * t413) * t419) * t453 + (Ifges(6,4) * t232 + Ifges(6,2) * t231 + Ifges(6,6) * t367) * t574 + (Ifges(6,1) * t232 + Ifges(6,4) * t231 + Ifges(6,5) * t367) * t575 + t232 * t578 + t231 * t579 + (Ifges(7,4) * t150 + Ifges(7,2) * t149 + Ifges(7,6) * t367) * t580 + t84 * t229 + t85 * t230 + t104 * (-mrSges(6,1) * t231 + mrSges(6,2) * t232) + t254 * t118 + (Ifges(7,4) * t42 + Ifges(7,2) * t43 + Ifges(7,6) * t306) * t567 + t266 * t260 + t267 * t261 + t191 * t292 + t192 * t293 + (Ifges(7,1) * t150 + Ifges(7,4) * t149 + Ifges(7,5) * t367) * t581 + t141 * (mrSges(5,1) * t306 - mrSges(5,3) * t219) + t142 * (-mrSges(5,2) * t306 + mrSges(5,3) * t218) + t59 * (mrSges(6,1) * t306 - mrSges(6,3) * t99) + t60 * (-mrSges(6,2) * t306 + mrSges(6,3) * t100) + t19 * (mrSges(7,1) * t306 - mrSges(7,3) * t42) + t20 * (-mrSges(7,2) * t306 + mrSges(7,3) * t43) - t306 * t241 / 0.2e1 + t307 * t242 / 0.2e1 + t308 * t97 / 0.2e1 + t316 * (mrSges(4,1) * t306 + mrSges(4,2) * t307) + t25 * t27 + t26 * t28 + t150 * t582 + t149 * t583 + t492 * t522 + t391 * t541 + t200 * t542 + t199 * t544 + (Ifges(5,5) * t219 + Ifges(5,6) * t218 + Ifges(5,3) * t306) * t545 + (Ifges(6,5) * t99 + Ifges(6,6) * t100 + Ifges(6,3) * t306) * t548 + (Ifges(7,5) * t42 + Ifges(7,6) * t43 + Ifges(7,3) * t306) * t550 + (Ifges(5,1) * t219 + Ifges(5,4) * t218 + Ifges(5,5) * t306) * t554 + (Ifges(5,4) * t219 + Ifges(5,2) * t218 + Ifges(5,6) * t306) * t556 + (Ifges(6,1) * t99 + Ifges(6,4) * t100 + Ifges(6,5) * t306) * t558; (-t350 * mrSges(4,1) - t86 / 0.2e1 + t199 / 0.2e1 - t96 / 0.2e1 - t32 / 0.2e1 - t31 / 0.2e1 - t8 / 0.2e1 - t39 / 0.2e1 - t87 / 0.2e1 - t187 / 0.2e1 - t188 / 0.2e1 + t514 / 0.2e1 + t172 * mrSges(4,3) + pkin(9) * t261 + (-Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1) * t289 + t585 - t592 - t603) * t423 + (-Ifges(6,5) * t364 + Ifges(7,5) * t272 - Ifges(6,6) * t363 + Ifges(7,6) * t270 + t418 * t441) * t553 + t617 * t168 + (t104 * t388 + t13 * t227 + t14 * t226 + t171 * t613 + t59 * t617 + t60 * t618) * m(6) + t618 * t167 - m(5) * (t141 * t182 + t142 * t183 + t224 * t256) - m(4) * (t245 * t271 + t246 * t273 + t316 * t360) + t104 * (mrSges(6,1) * t363 - mrSges(6,2) * t364) + (-t13 * t363 + t14 * t364 + t480 * t60 - t481 * t59) * mrSges(6,3) + (-Ifges(6,4) * t364 - Ifges(6,2) * t363) * t574 + (-Ifges(6,1) * t364 - Ifges(6,4) * t363) * t575 + (t236 / 0.2e1 - t250 / 0.2e1) * t110 + (-t163 / 0.2e1 + t130 / 0.2e1) * t67 + (t141 * t327 - t142 * t326) * mrSges(5,3) + (-t182 + t264) * t230 + ((t507 / 0.2e1 + (-m(4) * t245 + m(5) * t224 - t482) * pkin(9) + t332 / 0.2e1 - t612) * t423 + (-t425 + (-m(4) * t246 - t292) * pkin(9)) * t418) * qJD(3) + t613 * t131 + t614 * t75 + (-t164 / 0.2e1 + t129 / 0.2e1) * t68 + (t350 * mrSges(4,2) + t162 * t446 - t173 * mrSges(4,3) + t97 * t540 + t98 * t539 + t445 * t563 + t443 * t562 + t515 / 0.2e1 - t513 / 0.2e1 + t200 / 0.2e1 + (-t417 * t72 - t422 * t73) * mrSges(5,3) + (-m(4) * t173 + m(5) * t162 + t118 - t260) * pkin(9) + (t440 * t546 + t442 * t557 + t444 * t555 + t224 * t447 + t176 * t540 + t422 * t564 + (t141 * t417 - t142 * t422) * mrSges(5,3)) * qJD(4)) * t418 + t619 * t108 + t620 * t107 + (t112 * t614 + t119 * t3 + t120 * t2 + t19 * t619 + t20 * t620 + t300 * t51) * m(7) + (t235 / 0.2e1 - t251 / 0.2e1) * t111 + ((qJD(2) * (Ifges(4,5) * t418 + Ifges(4,6) * t423) / 0.2e1 + (t571 + t530 / 0.2e1) * t474 + (t401 / 0.2e1 - qJD(2)) * Ifges(3,6) - t429) * t419 + (-t396 / 0.2e1 + (t570 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t419) * t474 + (t547 - t507 / 0.2e1 + t431) * t423 + t425 * t418 - t432) * t424) * t474 + m(5) * (t141 * t264 + t142 * t263 + t346 * t73 + t347 * t72) + t388 * t44 - t357 * t356 + t346 * t153 + t347 * t154 - t349 * mrSges(3,2) - t350 * mrSges(3,1) - pkin(2) * t216 + m(4) * (-pkin(2) * t350 + t172 * t535) + (Ifges(6,4) * t235 + Ifges(6,2) * t236) * t560 + (Ifges(6,4) * t251 + Ifges(6,2) * t250) * t561 + t326 * t564 + (Ifges(7,1) * t129 + Ifges(7,4) * t130) * t565 + t120 * t28 + t119 * t27 + (-mrSges(7,1) * t485 + mrSges(7,2) * t486) * t112 + (-t19 * t486 + t2 * t270 + t20 * t485 - t272 * t3) * mrSges(7,3) + (-mrSges(6,1) * t480 + mrSges(6,2) * t481) * t171 - t477 * t360 - t364 * t578 - t363 * t579 + t226 * t78 + t227 * t79 + (Ifges(7,1) * t164 + Ifges(7,4) * t163) * t566 + (Ifges(7,4) * t129 + Ifges(7,2) * t130) * t567 + (Ifges(7,4) * t164 + Ifges(7,2) * t163) * t568 - t256 * t212 + t51 * (-mrSges(7,1) * t270 + mrSges(7,2) * t272) + t644 * t229 - t273 * t292 - t271 * t293 + t300 * t11 + (Ifges(7,4) * t272 + Ifges(7,2) * t270) * t580 + (Ifges(7,1) * t272 + Ifges(7,4) * t270) * t581 - t327 * t176 / 0.2e1 - t224 * (-mrSges(5,1) * t326 + mrSges(5,2) * t327) + t391 + t272 * t582 + t270 * t583 + (Ifges(5,5) * t327 + Ifges(5,6) * t326) * t546 + (Ifges(6,5) * t235 + Ifges(6,6) * t236) * t548 + (Ifges(6,5) * t251 + Ifges(6,6) * t250) * t549 + (Ifges(7,5) * t129 + Ifges(7,6) * t130) * t550 + (Ifges(7,5) * t164 + Ifges(7,6) * t163) * t551 + (Ifges(5,1) * t327 + Ifges(5,4) * t326) * t555 + (Ifges(5,4) * t327 + Ifges(5,2) * t326) * t557 + (Ifges(6,1) * t235 + Ifges(6,4) * t236) * t558 + (Ifges(6,1) * t251 + Ifges(6,4) * t250) * t559; (-Ifges(6,4) * t302 - Ifges(6,2) * t303) * t560 + (-Ifges(6,5) * t302 - Ifges(6,6) * t303) * t548 + (-Ifges(6,1) * t302 - Ifges(6,4) * t303) * t558 + (t248 / 0.2e1 - t303 / 0.2e1) * t110 + t104 * (mrSges(6,1) * t438 + mrSges(6,2) * t381) + (Ifges(6,5) * t381 + Ifges(7,5) * t297 - Ifges(6,6) * t438 + Ifges(7,6) * t296 + t440) * t553 + (Ifges(6,4) * t381 - Ifges(6,2) * t438) * t574 + (Ifges(6,1) * t381 - Ifges(6,4) * t438) * t575 - t438 * t579 + (t249 / 0.2e1 - t302 / 0.2e1) * t111 + t601 * t108 + (t112 * t587 + t19 * t601 + t197 * t3 + t198 * t2 + t20 * t602 + t351 * t51) * m(7) + t602 * t107 + (-t13 * t438 - t14 * t381 + t478 * t59 + t60 * t636) * mrSges(6,3) + (-mrSges(6,1) * t636 - mrSges(6,2) * t478) * t171 + t590 * t167 + (t104 * t410 + t13 * t321 + t14 * t320 - t171 * t209 + t59 * t591 + t590 * t60) * m(6) + t591 * t168 + t587 * t75 + (-t153 * t417 + t154 * t422 + (-m(5) * t439 - t417 * t229 - t422 * t230) * qJD(4) + m(5) * t586) * pkin(10) + t586 * mrSges(5,3) + (-Ifges(6,4) * t249 - Ifges(6,2) * t248) * t561 + (-Ifges(6,5) * t249 - Ifges(6,6) * t248) * t549 + (-Ifges(6,1) * t249 - Ifges(6,4) * t248) * t559 + t426 * t339 + t521 - t522 + (t547 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t339 + t612) * t338 + t410 * t44 + t351 * t11 - t209 * t131 + t197 * t27 + t198 * t28 + t442 * t562 + t444 * t563 + (Ifges(7,1) * t178 + Ifges(7,4) * t179) * t565 - pkin(3) * t118 + (mrSges(7,1) * t484 - mrSges(7,2) * t483) * t112 + (t19 * t483 + t2 * t296 - t20 * t484 - t297 * t3) * mrSges(7,3) + t482 * t246 + t417 * t573 + t381 * t578 - t170 * t229 - t169 * t230 + (t588 * t417 * pkin(4) - t624) * qJD(4) + t462 + (Ifges(7,1) * t160 + Ifges(7,4) * t159) * t566 + (Ifges(7,4) * t178 + Ifges(7,2) * t179) * t567 + (Ifges(7,4) * t160 + Ifges(7,2) * t159) * t568 - t245 * t292 + t51 * (-mrSges(7,1) * t296 + mrSges(7,2) * t297) + (Ifges(7,4) * t297 + Ifges(7,2) * t296) * t580 + (Ifges(7,1) * t297 + Ifges(7,4) * t296) * t581 + t297 * t582 + t320 * t78 + t321 * t79 + (t179 / 0.2e1 - t159 / 0.2e1) * t67 + (-pkin(3) * t162 - t141 * t169 - t142 * t170 - t224 * t246) * m(5) - t162 * t447 + (t178 / 0.2e1 - t160 / 0.2e1) * t68 + t296 * t583 + t97 * t539 + (Ifges(7,5) * t178 + Ifges(7,6) * t179) * t550 + (Ifges(7,5) * t160 + Ifges(7,6) * t159) * t551; -t585 + t207 * t606 + (t416 * t79 + t421 * t78 + (t167 * t421 - t168 * t416) * qJD(5) + m(6) * (t13 * t416 + t14 * t421 + t466 * t60 - t467 * t59) - t588 * t280) * pkin(4) + t615 * t107 + (-t112 * t156 + t19 * t616 + t2 * t366 + t20 * t615 + t3 * t365) * m(7) + t616 * t108 + t638 - m(6) * (t59 * t63 + t60 * t64) + (-Ifges(5,2) * t280 + t176 + t275) * t557 + t96 + t365 * t27 + t366 * t28 - t64 * t167 - t63 * t168 - t156 * t75 - t141 * t229 + t142 * t230 - t224 * (mrSges(5,1) * t280 + mrSges(5,2) * t279) + (t141 * t279 + t142 * t280) * mrSges(5,3) + (Ifges(5,5) * t279 - Ifges(5,6) * t280) * t546 + t175 * t554 + (Ifges(5,1) * t279 - t517) * t555; (-t207 * t75 + t420 * t27 + t415 * t28 + (t107 * t420 - t108 * t415) * qJD(6) + (-t112 * t207 - t19 * t465 + t2 * t415 + t20 * t464 + t3 * t420) * m(7)) * pkin(5) - m(7) * (t19 * t21 + t20 * t22) - t59 * t167 + t60 * t168 - t22 * t107 - t21 * t108 + t110 * t558 + t638; -t19 * t107 + t20 * t108 + t641;];
tauc  = t1(:);
