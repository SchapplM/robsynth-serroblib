% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR14_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:54
% EndTime: 2019-12-31 20:36:15
% DurationCPUTime: 9.54s
% Computational Cost: add. (23069->687), mult. (55859->983), div. (0->0), fcn. (61923->10), ass. (0->355)
t376 = sin(pkin(10));
t378 = cos(pkin(10));
t532 = sin(qJ(4));
t533 = cos(qJ(4));
t348 = -t533 * t376 - t532 * t378;
t382 = cos(qJ(5));
t483 = t382 * mrSges(6,2);
t380 = sin(qJ(5));
t485 = t380 * mrSges(6,1);
t354 = t483 + t485;
t276 = t354 * t348;
t610 = t348 * mrSges(5,3) + t276;
t379 = cos(pkin(5));
t377 = sin(pkin(5));
t381 = sin(qJ(2));
t466 = t377 * t381;
t338 = t376 * t379 + t378 * t466;
t451 = t376 * t466;
t464 = t378 * t379;
t405 = -t451 + t464;
t255 = t532 * t338 - t533 * t405;
t565 = t255 / 0.2e1;
t346 = t532 * t376 - t533 * t378;
t383 = cos(qJ(2));
t465 = t377 * t383;
t312 = t346 * t465;
t272 = t312 * t380 + t382 * t466;
t561 = t272 / 0.2e1;
t273 = -t312 * t382 + t380 * t466;
t560 = t273 / 0.2e1;
t311 = t348 * t465;
t556 = -t311 / 0.2e1;
t609 = pkin(9) * t255;
t395 = t533 * t338 + t532 * t405;
t602 = -t395 / 0.2e1;
t608 = mrSges(5,1) * t602;
t607 = m(6) * pkin(9) + mrSges(6,3);
t605 = -Ifges(5,6) / 0.2e1;
t453 = pkin(7) * t465;
t531 = pkin(1) * t381;
t319 = t453 + (qJ(3) + t531) * t379;
t320 = (-pkin(2) * t383 - qJ(3) * t381 - pkin(1)) * t377;
t236 = -t376 * t319 + t378 * t320;
t174 = -pkin(3) * t465 - t338 * pkin(8) + t236;
t237 = t378 * t319 + t376 * t320;
t188 = t405 * pkin(8) + t237;
t89 = t533 * t174 - t532 * t188;
t79 = pkin(4) * t465 - t89;
t604 = t79 / 0.2e1;
t603 = t311 / 0.2e1;
t555 = -t312 / 0.2e1;
t469 = t348 * t380;
t452 = mrSges(6,3) * t469;
t284 = -mrSges(6,2) * t346 + t452;
t457 = t382 * t284;
t601 = -t457 / 0.2e1;
t600 = pkin(4) * t395;
t599 = mrSges(5,3) * t395;
t374 = t380 ^ 2;
t375 = t382 ^ 2;
t424 = (t375 / 0.2e1 + t374 / 0.2e1) * mrSges(6,3);
t370 = Ifges(6,5) * t382;
t508 = Ifges(6,6) * t380;
t410 = -t370 / 0.2e1 + t508 / 0.2e1;
t598 = t410 * t255;
t219 = -t380 * t395 - t382 * t465;
t220 = -t380 * t465 + t382 * t395;
t118 = -mrSges(6,1) * t219 + mrSges(6,2) * t220;
t232 = -mrSges(5,1) * t465 - t599;
t597 = t118 - t232;
t371 = Ifges(6,4) * t382;
t359 = Ifges(6,1) * t380 + t371;
t596 = t376 * t338 + t378 * t405;
t468 = t348 * t382;
t286 = mrSges(6,1) * t346 + mrSges(6,3) * t468;
t460 = t380 * t286;
t594 = t346 * mrSges(5,3) + t460;
t420 = t508 - t370;
t421 = Ifges(6,2) * t380 - t371;
t486 = t378 * mrSges(4,2);
t488 = t376 * mrSges(4,1);
t318 = (t486 + t488) * t465;
t353 = -mrSges(6,1) * t382 + mrSges(6,2) * t380;
t593 = -m(6) * pkin(4) - mrSges(5,1) + t353;
t592 = Ifges(6,5) * t560 + Ifges(6,6) * t561 + Ifges(6,3) * t556;
t591 = m(4) / 0.2e1;
t590 = m(5) / 0.2e1;
t589 = -m(6) / 0.2e1;
t588 = m(6) / 0.2e1;
t587 = -pkin(4) / 0.2e1;
t586 = -mrSges(5,2) / 0.2e1;
t585 = mrSges(6,2) / 0.2e1;
t583 = Ifges(4,5) / 0.2e1;
t340 = (pkin(2) * t381 - qJ(3) * t383) * t377;
t365 = pkin(7) * t466;
t530 = pkin(1) * t383;
t341 = t379 * t530 - t365;
t264 = t378 * t340 - t376 * t341;
t463 = t378 * t383;
t210 = (pkin(3) * t381 - pkin(8) * t463) * t377 + t264;
t265 = t376 * t340 + t378 * t341;
t450 = t376 * t465;
t235 = -pkin(8) * t450 + t265;
t116 = t532 * t210 + t533 * t235;
t107 = pkin(9) * t466 + t116;
t342 = t379 * t531 + t453;
t308 = pkin(3) * t450 + t342;
t167 = -t311 * pkin(4) + t312 * pkin(9) + t308;
t57 = -t107 * t380 + t167 * t382;
t582 = -t57 / 0.2e1;
t58 = t107 * t382 + t167 * t380;
t581 = t58 / 0.2e1;
t503 = t219 * mrSges(6,3);
t141 = -mrSges(6,2) * t255 + t503;
t580 = -t141 / 0.2e1;
t369 = -pkin(3) * t378 - pkin(2);
t280 = pkin(4) * t346 + pkin(9) * t348 + t369;
t524 = pkin(8) + qJ(3);
t351 = t524 * t378;
t431 = t524 * t376;
t296 = t533 * t351 - t532 * t431;
t179 = t280 * t382 - t296 * t380;
t579 = -t179 / 0.2e1;
t180 = t280 * t380 + t296 * t382;
t578 = -t180 / 0.2e1;
t189 = mrSges(6,2) * t311 + mrSges(6,3) * t272;
t577 = t189 / 0.2e1;
t190 = -mrSges(6,1) * t311 - mrSges(6,3) * t273;
t576 = -t190 / 0.2e1;
t527 = pkin(9) * t346;
t528 = pkin(4) * t348;
t291 = t527 - t528;
t295 = t532 * t351 + t533 * t431;
t194 = t291 * t382 + t295 * t380;
t575 = -t194 / 0.2e1;
t574 = t219 / 0.2e1;
t573 = t220 / 0.2e1;
t491 = t346 * Ifges(6,6);
t228 = t421 * t348 + t491;
t572 = -t228 / 0.4e1;
t517 = Ifges(6,4) * t380;
t360 = Ifges(6,1) * t382 - t517;
t492 = t346 * Ifges(6,5);
t230 = -t348 * t360 + t492;
t571 = t230 / 0.4e1;
t570 = t232 / 0.2e1;
t568 = -t255 / 0.2e1;
t563 = t395 / 0.2e1;
t559 = -t284 / 0.2e1;
t558 = -t286 / 0.2e1;
t557 = t295 / 0.2e1;
t554 = -t346 / 0.2e1;
t553 = -t346 / 0.4e1;
t552 = t346 / 0.4e1;
t551 = -t348 / 0.2e1;
t550 = -t348 / 0.4e1;
t549 = t348 / 0.4e1;
t548 = -t353 / 0.2e1;
t547 = t353 / 0.2e1;
t546 = t354 / 0.2e1;
t355 = Ifges(6,5) * t380 + Ifges(6,6) * t382;
t545 = t355 / 0.4e1;
t544 = -t420 / 0.4e1;
t543 = t359 / 0.4e1;
t542 = t369 / 0.2e1;
t541 = -t376 / 0.2e1;
t540 = t378 / 0.2e1;
t539 = -t380 / 0.2e1;
t538 = t380 / 0.2e1;
t537 = t380 / 0.4e1;
t536 = -t382 / 0.2e1;
t535 = t382 / 0.2e1;
t534 = t382 / 0.4e1;
t115 = t533 * t210 - t532 * t235;
t106 = -pkin(4) * t466 - t115;
t529 = pkin(4) * t106;
t526 = pkin(9) * t380;
t525 = pkin(9) * t382;
t522 = Ifges(4,4) * t376;
t521 = Ifges(4,4) * t378;
t520 = Ifges(5,4) * t395;
t519 = Ifges(5,4) * t348;
t518 = Ifges(6,4) * t220;
t516 = Ifges(5,5) * t312;
t515 = Ifges(6,5) * t255;
t512 = Ifges(4,6) * t376;
t511 = Ifges(5,6) * t348;
t510 = Ifges(6,6) * t255;
t507 = Ifges(6,3) * t395;
t505 = Ifges(6,3) * t348;
t504 = t180 * mrSges(6,3);
t502 = t220 * mrSges(6,3);
t501 = t255 * mrSges(5,2);
t500 = t395 * mrSges(5,1);
t499 = t295 * mrSges(5,3);
t498 = t296 * mrSges(5,3);
t131 = Ifges(6,4) * t273 + Ifges(6,2) * t272 - Ifges(6,6) * t311;
t132 = Ifges(6,1) * t273 + Ifges(6,4) * t272 - Ifges(6,5) * t311;
t142 = mrSges(6,1) * t255 - t502;
t251 = Ifges(5,4) * t255;
t149 = Ifges(5,1) * t395 - Ifges(5,5) * t465 - t251;
t162 = -mrSges(6,1) * t272 + mrSges(6,2) * t273;
t495 = t312 * mrSges(5,2);
t496 = t311 * mrSges(5,1);
t221 = -t495 - t496;
t231 = mrSges(5,2) * t465 - t255 * mrSges(5,3);
t271 = pkin(3) * t451 + t365 + (t369 - t530) * t379;
t281 = -mrSges(5,2) * t466 + mrSges(5,3) * t311;
t282 = mrSges(5,1) * t466 + mrSges(5,3) * t312;
t487 = t376 * Ifges(4,2);
t293 = (Ifges(4,6) * t381 + (-t487 + t521) * t383) * t377;
t306 = mrSges(4,2) * t465 + t405 * mrSges(4,3);
t307 = -mrSges(4,1) * t465 - t338 * mrSges(4,3);
t326 = t365 + (-pkin(2) - t530) * t379;
t336 = (-mrSges(4,3) * t376 * t383 - mrSges(4,2) * t381) * t377;
t337 = (mrSges(4,1) * t381 - mrSges(4,3) * t463) * t377;
t363 = Ifges(3,5) * t465;
t411 = t308 * mrSges(5,1) + Ifges(5,4) * t312 / 0.2e1 + Ifges(5,2) * t556 + t466 * t605 + t592;
t427 = t342 * mrSges(4,2) + (t381 * Ifges(4,5) + (t378 * Ifges(4,1) - t522) * t383) * t377 / 0.2e1;
t428 = t342 * mrSges(4,1) - t293 / 0.2e1;
t429 = t308 * mrSges(5,2) + Ifges(5,1) * t555 + Ifges(5,4) * t603 + Ifges(5,5) * t466 / 0.2e1;
t104 = t255 * pkin(4) - pkin(9) * t395 + t271;
t90 = t532 * t174 + t533 * t188;
t80 = -pkin(9) * t465 + t90;
t43 = t104 * t382 - t380 * t80;
t430 = Ifges(4,6) * t540 - Ifges(3,6);
t44 = t104 * t380 + t382 * t80;
t148 = -Ifges(5,2) * t255 - Ifges(5,6) * t465 + t520;
t82 = Ifges(6,5) * t220 + Ifges(6,6) * t219 + t255 * Ifges(6,3);
t441 = t82 / 0.2e1 - t148 / 0.2e1;
t494 = t341 * mrSges(3,2);
t83 = Ifges(6,2) * t219 + t510 + t518;
t216 = Ifges(6,4) * t219;
t84 = Ifges(6,1) * t220 + t216 + t515;
t3 = ((Ifges(5,5) * t563 + Ifges(5,6) * t568 + t338 * t583 + t430 * t379 + t428 * t376 + (-pkin(1) * mrSges(3,1) + (-t512 / 0.2e1 - Ifges(3,4)) * t381) * t377) * t381 + (t516 / 0.2e1 + Ifges(5,6) * t556 + Ifges(3,5) * t379 / 0.2e1 + Ifges(3,4) * t465 + (Ifges(4,4) * t338 + Ifges(4,2) * t464) * t541 + (Ifges(4,1) * t338 + Ifges(4,4) * t464) * t540 + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t378 + t512) * t383) * t377 + (-Ifges(5,3) - Ifges(4,3) - Ifges(3,2) + Ifges(3,1) + (-t521 / 0.2e1 + t487 / 0.2e1) * t376) * t466) * t383) * t377 + t429 * t395 + t149 * t555 + t84 * t560 + t83 * t561 + m(4) * (t236 * t264 + t237 * t265 + t326 * t342) + m(6) * (t106 * t79 + t43 * t57 + t44 * t58) + m(5) * (t115 * t89 + t116 * t90 + t271 * t308) + t237 * t336 + t236 * t337 + t326 * t318 + t265 * t306 + t264 * t307 + t89 * t282 + t90 * t281 + t271 * t221 + t116 * t231 + t115 * t232 + t43 * t190 + t44 * t189 + t79 * t162 + t58 * t141 + t57 * t142 + (-t494 + t363 / 0.2e1 + t293 * t540 + (-t378 * mrSges(4,1) - mrSges(3,1)) * t342) * t379 + t106 * t118 - t441 * t311 + t132 * t573 + t131 * t574 + t411 * t255 + t427 * t338;
t497 = t3 * qJD(1);
t489 = t369 * mrSges(5,1);
t484 = t380 * t43;
t482 = t382 * mrSges(6,3);
t481 = t382 * t44;
t154 = t354 * t255;
t156 = mrSges(6,3) * t255 * t380 - mrSges(6,2) * t395;
t157 = mrSges(6,1) * t395 + t255 * t482;
t426 = -Ifges(6,3) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1;
t456 = -Ifges(5,5) * t255 - Ifges(5,6) * t395;
t158 = t600 + t609;
t51 = t158 * t382 - t380 * t89;
t52 = t158 * t380 + t382 * t89;
t95 = Ifges(6,6) * t395 + t421 * t255;
t96 = Ifges(6,5) * t395 - t360 * t255;
t4 = m(6) * (t43 * t51 + t44 * t52) + t52 * t141 + t44 * t156 + t51 * t142 + t43 * t157 - t79 * t154 + t95 * t574 + t96 * t573 + t89 * t231 - t456 * t465 / 0.2e1 + (t271 * mrSges(5,1) - t520 / 0.2e1 + t441) * t395 + (-t271 * mrSges(5,2) - t149 / 0.2e1 + t251 / 0.2e1 + t83 * t538 + t84 * t536 + t89 * mrSges(5,3) + t598 - t426 * t395) * t255 + (m(6) * t79 + t597 - t599) * t90;
t480 = t4 * qJD(1);
t479 = t57 * t380;
t478 = t58 * t382;
t117 = mrSges(6,1) * t220 + mrSges(6,2) * t219;
t119 = Ifges(6,5) * t219 - Ifges(6,6) * t220;
t120 = -Ifges(6,2) * t220 + t216;
t121 = Ifges(6,1) * t219 - t518;
t7 = t43 * t141 - t44 * t142 + t79 * t117 + t119 * t565 + (-t44 * mrSges(6,3) - t83 / 0.2e1 + t121 / 0.2e1) * t220 + (-t43 * mrSges(6,3) + t120 / 0.2e1 + t84 / 0.2e1) * t219;
t477 = t7 * qJD(1);
t419 = -t481 + t484;
t459 = t382 * t141;
t462 = t380 * t142;
t12 = t405 * t306 - t338 * t307 + t597 * t395 - (t231 + t459 - t462) * t255 + m(6) * (t255 * t419 + t395 * t79) + m(5) * (-t255 * t90 - t395 * t89) + m(4) * (-t236 * t338 + t237 * t405);
t476 = qJD(1) * t12;
t473 = t395 * t295;
t472 = t295 * t348;
t471 = t346 * t380;
t470 = t346 * t382;
t461 = t380 * t228;
t458 = t382 * t230;
t454 = t374 + t375;
t449 = -t519 / 0.2e1;
t446 = Ifges(5,3) * t381 / 0.2e1;
t445 = -t502 / 0.2e1;
t444 = t492 / 0.2e1;
t443 = -t491 / 0.2e1;
t442 = t120 / 0.4e1 + t84 / 0.4e1;
t440 = -t83 / 0.4e1 + t121 / 0.4e1;
t439 = t471 / 0.2e1;
t438 = -t470 / 0.2e1;
t278 = t359 * t348;
t436 = t572 + t278 / 0.4e1;
t357 = Ifges(6,2) * t382 + t517;
t277 = t348 * t357;
t435 = t571 + t277 / 0.4e1;
t434 = t543 - t421 / 0.4e1;
t433 = t360 / 0.4e1 - t357 / 0.4e1;
t425 = Ifges(6,3) / 0.4e1 - Ifges(5,1) / 0.4e1 + Ifges(5,2) / 0.4e1;
t423 = t580 + t503 / 0.2e1;
t422 = t445 - t142 / 0.2e1;
t195 = t291 * t380 - t295 * t382;
t226 = t346 * Ifges(6,3) + t420 * t348;
t227 = -Ifges(6,6) * t348 + t421 * t346;
t229 = -Ifges(6,5) * t348 - t360 * t346;
t275 = t354 * t346;
t283 = mrSges(6,2) * t348 + mrSges(6,3) * t471;
t285 = -mrSges(6,1) * t348 + mrSges(6,3) * t470;
t289 = -t346 * Ifges(5,2) - t519;
t344 = Ifges(5,4) * t346;
t290 = -t348 * Ifges(5,1) - t344;
t343 = t346 * mrSges(5,2);
t13 = t194 * t286 + t179 * t285 + t195 * t284 + t180 * t283 - t296 * t276 - t295 * t275 + m(6) * (t179 * t194 + t180 * t195 + t295 * t296) - t369 * t343 + (t229 * t536 + t227 * t538 - t226 / 0.2e1 - t489 + t289 / 0.2e1 + t449) * t348 + (-t458 / 0.2e1 + t461 / 0.2e1 - t290 / 0.2e1 + t344 / 0.2e1 + t410 * t346 + t426 * t348) * t346;
t386 = (t179 * t51 + t180 * t52 + t194 * t43 + t195 * t44 + t295 * t90 + t296 * t79) * t589 + t157 * t579 + t156 * t578 + t142 * t575 + t195 * t580 - t219 * t227 / 0.4e1 - t220 * t229 / 0.4e1 - t271 * (-t348 * mrSges(5,1) - t343) / 0.2e1 + t148 * t550 + t82 * t549 - t43 * t285 / 0.2e1 - t44 * t283 / 0.2e1 + t51 * t558 + t52 * t559 + t275 * t604 + t90 * t276 / 0.2e1;
t387 = -(t545 + t605) * t311 - t516 / 0.2e1 + t162 * t587 + t106 * t547 + t115 * mrSges(5,1) / 0.2e1 + t116 * t586 + t272 * t357 / 0.4e1 + t273 * t543;
t402 = t383 * (-Ifges(5,5) * t346 + t511);
t2 = t387 + t386 + (-t489 / 0.2e1 + t449 + t289 / 0.4e1 - t226 / 0.4e1 + t498 / 0.2e1 - t425 * t346) * t395 + (-t118 / 0.2e1 + t570) * t296 + (t154 / 0.2e1 + t231 / 0.2e1) * t295 + (t402 / 0.4e1 + t446) * t377 + (t132 / 0.4e1 + t95 * t550 + t83 * t553 + mrSges(6,3) * t582 + (t443 + t572) * t255 + (m(6) * t582 + t576) * pkin(9)) * t380 + (t131 / 0.4e1 + t96 * t549 + t84 * t552 + mrSges(6,3) * t581 + (t444 + t571) * t255 + (m(6) * t581 + t577) * pkin(9)) * t382 + t529 * t589 + (-t251 / 0.4e1 + t149 / 0.4e1) * t346 + (mrSges(5,2) * t542 - t344 / 0.4e1 + t290 / 0.4e1 + t499 / 0.2e1 + t425 * t348) * t255;
t417 = -t2 * qJD(1) + t13 * qJD(2);
t274 = t348 * t353;
t19 = t180 * t286 - t295 * t274 + ((-t230 / 0.2e1 - t277 / 0.2e1 - t492 / 0.2e1) * t380 + (t278 / 0.2e1 - t228 / 0.2e1 + t443 - t504) * t382) * t348 + (t452 - t284) * t179;
t385 = (mrSges(6,3) * t579 + t435) * t219 + (-t504 / 0.2e1 + t436) * t220 + (t255 * t545 + t83 * t534 - t382 * t121 / 0.4e1 + (-t484 / 0.2e1 + t481 / 0.2e1) * mrSges(6,3) + (t120 + t84) * t537) * t348 + t179 * t141 / 0.2e1 + t142 * t578 + t117 * t557 + t119 * t552 + t43 * t284 / 0.2e1 + t44 * t558 + t274 * t604;
t393 = t57 * mrSges(6,1) / 0.2e1 - t58 * mrSges(6,2) / 0.2e1 + t592;
t5 = t385 - t393;
t416 = t5 * qJD(1) - t19 * qJD(2);
t412 = t179 * t380 - t180 * t382;
t28 = t610 * t348 + (-t457 + t594) * t346 + m(6) * (t412 * t346 - t472) + m(5) * (-t296 * t346 - t472) + (m(4) * qJ(3) + mrSges(4,3)) * (t376 ^ 2 + t378 ^ 2);
t384 = (t596 * qJ(3) - t376 * t236 + t378 * t237) * t591 + (-t255 * t296 - t346 * t90 + t348 * t89 + t473) * t590 + (t255 * t412 + t419 * t346 - t348 * t79 + t473) * t588 + t231 * t554 + t118 * t551 + t348 * t570 + t307 * t541 + t306 * t540 + t142 * t439 + t255 * t601 + t141 * t438 + t596 * mrSges(4,3) / 0.2e1 + t594 * t565 + t610 * t602;
t388 = t342 * t591 + t308 * t590 + (t380 * t58 + t382 * t57) * t588 - t496 / 0.2e1 - t495 / 0.2e1 + t189 * t538 + t190 * t535;
t9 = -t384 + t388 + t318 / 0.2e1;
t415 = qJD(1) * t9 - qJD(2) * t28;
t390 = (t380 * t52 + t382 * t51) * t589 + t501 / 0.2e1 + t156 * t539 + t157 * t536;
t396 = (-t454 * t609 - t600) * t588 + t395 * t547;
t16 = 0.2e1 * t608 - (t586 + t424) * t255 + t390 + t396;
t392 = -t346 * t424 + (-t454 * t527 + t528) * t588;
t394 = (t194 * t382 + t195 * t380) * t588 + t283 * t538 + t285 * t535;
t32 = t343 + (t548 + mrSges(5,1)) * t348 + t392 - t394;
t414 = qJD(1) * t16 + qJD(2) * t32;
t22 = (mrSges(6,2) * t565 + t423) * t382 + (mrSges(6,1) * t565 - t422) * t380;
t408 = -t483 / 0.2e1 - t485 / 0.2e1;
t401 = t408 * t346;
t404 = t601 + t460 / 0.2e1;
t47 = -t401 + t404;
t413 = qJD(1) * t22 + qJD(2) * t47;
t409 = mrSges(6,1) * t575 + t195 * t585;
t407 = pkin(9) * t558 + t435;
t406 = pkin(9) * t559 + t436;
t403 = t357 * t538 + t359 * t536;
t389 = t117 * t587 + t434 * t219 + t433 * t220 + t255 * t544 + t79 * t546;
t398 = -t507 / 0.2e1 - t51 * mrSges(6,1) / 0.2e1 + t52 * t585;
t10 = (t515 / 0.2e1 + t422 * pkin(9) + t442) * t382 + (-t510 / 0.2e1 + t423 * pkin(9) + t440) * t380 + t389 + t398;
t193 = -pkin(4) * t354 + (t359 / 0.2e1 - t421 / 0.2e1) * t382 + (t360 / 0.2e1 - t357 / 0.2e1) * t380;
t391 = pkin(9) * t424 + t434 * t380 - t433 * t382;
t397 = t274 * t587 + t295 * t546 + t346 * t544;
t21 = (t444 + t407) * t382 + (t443 + t406) * t380 + (Ifges(6,3) / 0.2e1 + t391) * t348 + t397 + t409;
t399 = t10 * qJD(1) + t21 * qJD(2) + t193 * qJD(4);
t48 = -t401 - t404;
t33 = t348 * t548 + t392 + t394;
t23 = t459 / 0.2e1 + t380 * t445 - t462 / 0.2e1 - t219 * t482 / 0.2e1 - t408 * t255;
t20 = Ifges(6,5) * t438 + Ifges(6,6) * t439 - t505 / 0.2e1 + t407 * t382 + t406 * t380 + t391 * t348 + t397 - t409;
t17 = t500 / 0.2e1 + mrSges(5,2) * t565 + t608 - t255 * t424 - t390 + t396;
t11 = t442 * t382 + t440 * t380 + (t141 * t539 + t142 * t536 + (t219 * t538 + t220 * t536) * mrSges(6,3)) * pkin(9) + t389 - t398 + t598;
t8 = t384 + (t486 / 0.2e1 + t488 / 0.2e1) * t465 + t388;
t6 = t385 + t393;
t1 = (t346 * t420 + t461 - t505) * t255 / 0.4e1 - (Ifges(5,2) * t348 + t290 - t344 + t458) * t255 / 0.4e1 + (-Ifges(5,1) * t346 + t226 + t519) * t395 / 0.4e1 - t395 * t289 / 0.4e1 + (-Ifges(5,2) * t395 + t149 - t251) * t553 - t154 * t557 + t499 * t568 + (-Ifges(5,1) * t255 - t520) * t550 + (t255 * t420 + t507) * t552 + t387 - t386 + (t500 - t501) * t542 + t498 * t602 + t132 * t537 + t131 * t534 - t295 * t231 / 0.2e1 + (t478 / 0.2e1 - t479 / 0.2e1) * mrSges(6,3) + t83 * t471 / 0.4e1 - t96 * t468 / 0.4e1 + t95 * t469 / 0.4e1 - t84 * t470 / 0.4e1 + (-t402 / 0.4e1 + t446) * t377 + t526 * t576 + t525 * t577 + (-t529 + (t478 - t479) * pkin(9)) * t588 + (t118 / 0.2e1 - t232 / 0.2e1) * t296;
t14 = [qJD(2) * t3 + qJD(3) * t12 + qJD(4) * t4 + qJD(5) * t7, t8 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + t497 + (t289 * t603 + t226 * t556 + t230 * t560 + t228 * t561 + t290 * t555 + t363 + (t265 * mrSges(4,3) + qJ(3) * t336 - t428) * t378 + (-t264 * mrSges(4,3) - qJ(3) * t337 + t427) * t376 + (-t116 * mrSges(5,3) + t411) * t346 + (t115 * mrSges(5,3) + t131 * t538 + t132 * t536 - t429) * t348 + ((Ifges(5,5) * t551 + Ifges(5,6) * t554 + t376 * t583 + t430) * t381 + ((Ifges(4,2) * t378 + t522) * t541 + (Ifges(4,1) * t376 + t521) * t540) * t383) * t377 + t369 * t221 - t342 * mrSges(3,1) - pkin(2) * t318 + t295 * t162 - t295 * t282 + t296 * t281 + t58 * t284 + t57 * t286 - t106 * t276 - t494 + t179 * t190 + t180 * t189 + 0.2e1 * (t106 * t295 + t179 * t57 + t180 * t58) * t588 + 0.2e1 * (-t115 * t295 + t116 * t296 + t308 * t369) * t590 + 0.2e1 * (-pkin(2) * t342 + (-t264 * t376 + t265 * t378) * qJ(3)) * t591) * qJD(2), qJD(2) * t8 + qJD(4) * t17 + qJD(5) * t23 + t476, t480 + t1 * qJD(2) + t17 * qJD(3) + (-t89 * mrSges(5,2) + pkin(4) * t154 + t156 * t525 - t157 * t526 + t403 * t255 + t355 * t563 + t95 * t535 + t96 * t538 + t456 + t607 * (-t380 * t51 + t382 * t52) + t593 * t90) * qJD(4) + t11 * qJD(5), t477 + t6 * qJD(2) + t23 * qJD(3) + t11 * qJD(4) + (-mrSges(6,1) * t44 - mrSges(6,2) * t43 + t119) * qJD(5); -qJD(3) * t9 - qJD(4) * t2 + qJD(5) * t5 - t497, qJD(3) * t28 + qJD(4) * t13 - qJD(5) * t19, qJD(4) * t33 + qJD(5) * t48 - t415, t33 * qJD(3) + t20 * qJD(5) + t417 + (t355 * t551 + t229 * t538 + t227 * t535 + pkin(4) * t275 + t283 * t525 - t285 * t526 + t295 * mrSges(5,2) + t511 + (-Ifges(5,5) + t403) * t346 + t593 * t296 + t607 * (-t194 * t380 + t195 * t382)) * qJD(4), t48 * qJD(3) + t20 * qJD(4) + (-mrSges(6,1) * t180 - mrSges(6,2) * t179 + t348 * t355) * qJD(5) + t416; qJD(2) * t9 - qJD(4) * t16 - qJD(5) * t22 - t476, -qJD(4) * t32 - qJD(5) * t47 + t415, 0, -t414, -qJD(5) * t354 - t413; qJD(2) * t2 + qJD(3) * t16 + qJD(5) * t10 - t480, qJD(3) * t32 + qJD(5) * t21 - t417, t414, t193 * qJD(5), (pkin(9) * t353 - t420) * qJD(5) + t399; -qJD(2) * t5 + qJD(3) * t22 - qJD(4) * t10 - t477, qJD(3) * t47 - qJD(4) * t21 - t416, t413, -t399, 0;];
Cq = t14;
