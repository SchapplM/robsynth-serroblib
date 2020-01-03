% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:24
% EndTime: 2019-12-31 21:43:48
% DurationCPUTime: 11.28s
% Computational Cost: add. (13692->707), mult. (33209->985), div. (0->0), fcn. (33161->8), ass. (0->362)
t657 = -mrSges(4,1) + mrSges(5,2);
t656 = mrSges(4,2) - mrSges(5,3);
t655 = Ifges(4,1) + Ifges(6,3);
t654 = Ifges(5,4) - Ifges(4,5);
t648 = Ifges(5,5) - Ifges(4,6);
t631 = m(6) / 0.2e1;
t653 = 0.2e1 * t631;
t633 = m(5) / 0.2e1;
t652 = 0.2e1 * t633;
t403 = sin(pkin(5));
t651 = -t403 / 0.2e1;
t650 = Ifges(5,1) + Ifges(4,3);
t405 = sin(qJ(3));
t408 = cos(qJ(3));
t406 = sin(qJ(2));
t526 = t403 * t406;
t531 = cos(pkin(5));
t334 = t405 * t526 - t408 * t531;
t321 = Ifges(4,4) * t334;
t335 = t405 * t531 + t408 * t526;
t409 = cos(qJ(2));
t525 = t403 * t409;
t172 = t335 * Ifges(4,1) - Ifges(4,5) * t525 - t321;
t404 = sin(qJ(5));
t407 = cos(qJ(5));
t236 = t334 * t407 + t404 * t525;
t511 = t407 * t409;
t237 = -t334 * t404 + t403 * t511;
t77 = -t237 * Ifges(6,5) + t236 * Ifges(6,6) + t335 * Ifges(6,3);
t647 = t172 + t77;
t488 = pkin(1) * t531;
t472 = t406 * t488;
t338 = pkin(7) * t525 + t472;
t296 = pkin(8) * t531 + t338;
t297 = (-pkin(2) * t409 - pkin(8) * t406 - pkin(1)) * t403;
t158 = t408 * t296 + t405 * t297;
t498 = qJ(4) * t525;
t129 = t498 - t158;
t646 = t129 + t158;
t476 = -qJ(4) * t405 - pkin(2);
t369 = -pkin(3) * t408 + t476;
t373 = mrSges(5,2) * t408 - mrSges(5,3) * t405;
t645 = -m(5) * t369 - t373;
t554 = Ifges(6,6) * t407;
t562 = Ifges(6,5) * t404;
t457 = t554 + t562;
t567 = Ifges(4,4) * t408;
t644 = t655 * t405 - t408 * t457 + t567;
t577 = pkin(8) * t408;
t502 = -t577 / 0.2e1;
t578 = pkin(8) * t405;
t503 = -t578 / 0.2e1;
t643 = t334 * t503 + t335 * t502;
t287 = (-t404 * t406 + t405 * t511) * t403;
t517 = t405 * t409;
t288 = (t404 * t517 + t406 * t407) * t403;
t510 = t408 * t409;
t499 = t403 * t510;
t122 = Ifges(6,1) * t288 + Ifges(6,4) * t287 + Ifges(6,5) * t499;
t214 = mrSges(6,1) * t499 - mrSges(6,3) * t288;
t622 = pkin(3) + pkin(9);
t581 = -t622 / 0.2e1;
t642 = t214 * t581 + t122 / 0.4e1;
t121 = Ifges(6,4) * t288 + Ifges(6,2) * t287 + Ifges(6,6) * t499;
t213 = -mrSges(6,2) * t499 + mrSges(6,3) * t287;
t641 = t213 * t581 - t121 / 0.4e1;
t640 = t654 * t334 + t648 * t335;
t529 = qJ(4) * t408;
t358 = t405 * t622 - t529;
t621 = pkin(4) + pkin(8);
t382 = t621 * t408;
t234 = -t358 * t404 + t382 * t407;
t515 = t407 * t234;
t235 = t358 * t407 + t382 * t404;
t523 = t404 * t235;
t639 = -t523 - t515;
t336 = -pkin(7) * t526 + t409 * t488;
t638 = -t338 * mrSges(3,1) - t336 * mrSges(3,2);
t337 = (pkin(2) * t406 - pkin(8) * t409) * t403;
t201 = t408 * t336 + t405 * t337;
t128 = (-pkin(4) * t517 + qJ(4) * t406) * t403 + t201;
t159 = -mrSges(6,1) * t287 + mrSges(6,2) * t288;
t162 = -qJ(4) * t526 - t201;
t200 = -t405 * t336 + t337 * t408;
t164 = -pkin(3) * t526 - t200;
t313 = (mrSges(5,1) * t517 - mrSges(5,3) * t406) * t403;
t370 = mrSges(5,1) * t499;
t314 = mrSges(5,2) * t526 + t370;
t535 = t407 * mrSges(6,2);
t541 = t404 * mrSges(6,1);
t374 = t535 + t541;
t112 = (pkin(4) * t510 - t406 * t622) * t403 - t200;
t399 = pkin(3) * t525;
t442 = t405 * t399 + t472;
t475 = pkin(7) - t529;
t163 = (pkin(9) * t405 + t475) * t525 + t442;
t57 = t112 * t407 - t163 * t404;
t58 = t112 * t404 + t163 * t407;
t449 = t58 * t404 + t57 * t407;
t566 = Ifges(6,4) * t404;
t380 = Ifges(6,1) * t407 - t566;
t595 = t380 / 0.4e1;
t565 = Ifges(6,4) * t407;
t378 = -Ifges(6,2) * t404 + t565;
t596 = t378 / 0.4e1;
t637 = (-pkin(3) * t164 - qJ(4) * t162) * t633 + (qJ(4) * t128 - t449 * t622) * t631 - pkin(3) * t314 / 0.2e1 + t128 * t374 / 0.2e1 - t162 * mrSges(5,3) / 0.2e1 + t164 * mrSges(5,2) / 0.2e1 + t200 * mrSges(4,1) / 0.2e1 - t201 * mrSges(4,2) / 0.2e1 + t287 * t596 + t288 * t595 + (-t313 / 0.2e1 + t159 / 0.2e1) * qJ(4);
t635 = t407 ^ 2;
t634 = -m(5) / 0.2e1;
t632 = -m(6) / 0.2e1;
t630 = -mrSges(6,1) / 0.2e1;
t629 = mrSges(6,2) / 0.2e1;
t628 = -mrSges(6,3) / 0.2e1;
t627 = Ifges(6,5) / 0.2e1;
t626 = -t57 / 0.2e1;
t547 = t237 * Ifges(6,4);
t556 = Ifges(6,6) * t335;
t78 = t236 * Ifges(6,2) - t547 + t556;
t625 = t78 / 0.4e1;
t233 = Ifges(6,4) * t236;
t563 = Ifges(6,5) * t335;
t79 = -t237 * Ifges(6,1) + t233 + t563;
t624 = t79 / 0.2e1;
t579 = pkin(4) * t334;
t90 = -t129 - t579;
t623 = t90 / 0.2e1;
t544 = t335 * mrSges(6,1);
t548 = t237 * mrSges(6,3);
t147 = t544 + t548;
t618 = -t147 / 0.2e1;
t343 = -t408 * t622 + t476;
t381 = t621 * t405;
t229 = -t343 * t404 + t381 * t407;
t616 = t229 / 0.2e1;
t230 = t343 * t407 + t381 * t404;
t615 = t230 / 0.2e1;
t614 = t235 / 0.2e1;
t613 = t236 / 0.2e1;
t612 = -t237 / 0.2e1;
t559 = Ifges(5,6) * t405;
t456 = -Ifges(5,2) * t408 + t559;
t611 = (t406 * Ifges(5,4) + t409 * t456) * t651;
t610 = t287 / 0.2e1;
t571 = Ifges(6,1) * t404;
t461 = t565 + t571;
t561 = Ifges(6,5) * t405;
t306 = -t408 * t461 + t561;
t609 = t306 / 0.2e1;
t607 = -t334 / 0.2e1;
t603 = t335 / 0.2e1;
t602 = t335 / 0.4e1;
t340 = t374 * t408;
t601 = -t340 / 0.2e1;
t519 = t404 * t408;
t505 = mrSges(6,3) * t519;
t538 = t405 * mrSges(6,1);
t359 = t505 + t538;
t600 = -t359 / 0.2e1;
t599 = t359 / 0.2e1;
t512 = t407 * t408;
t504 = mrSges(6,3) * t512;
t537 = t405 * mrSges(6,2);
t360 = -t504 - t537;
t598 = t360 / 0.2e1;
t376 = Ifges(6,5) * t407 - Ifges(6,6) * t404;
t597 = t376 / 0.4e1;
t594 = t382 / 0.2e1;
t593 = -t404 / 0.2e1;
t592 = t404 / 0.2e1;
t591 = -t405 / 0.2e1;
t589 = t405 / 0.2e1;
t588 = t405 / 0.4e1;
t587 = -t407 / 0.2e1;
t586 = t407 / 0.2e1;
t584 = t408 / 0.2e1;
t582 = -t409 / 0.2e1;
t157 = t296 * t405 - t408 * t297;
t102 = -pkin(4) * t335 - t157;
t81 = pkin(9) * t525 - t102 + t399;
t295 = -pkin(2) * t531 - t336;
t422 = -t335 * qJ(4) + t295;
t86 = t334 * t622 + t422;
t38 = -t404 * t86 + t407 * t81;
t576 = t38 * mrSges(6,3);
t39 = t404 * t81 + t407 * t86;
t575 = t39 * mrSges(6,3);
t574 = Ifges(4,4) + Ifges(5,6);
t573 = mrSges(6,1) * t622;
t572 = mrSges(6,2) * t622;
t570 = Ifges(3,4) * t406;
t569 = Ifges(4,4) * t335;
t568 = Ifges(4,4) * t405;
t564 = Ifges(6,5) * t288;
t560 = Ifges(5,6) * t335;
t558 = Ifges(5,6) * t408;
t557 = Ifges(6,6) * t287;
t555 = Ifges(6,6) * t405;
t553 = Ifges(6,3) * t334;
t552 = Ifges(6,3) * t408;
t551 = t236 * mrSges(6,1);
t550 = t236 * mrSges(6,3);
t549 = t237 * mrSges(6,2);
t108 = -t549 - t551;
t120 = Ifges(6,3) * t499 + t557 + t564;
t127 = t334 * pkin(3) + t422;
t130 = t157 + t399;
t543 = t335 * mrSges(6,2);
t146 = -t543 + t550;
t169 = -Ifges(5,5) * t525 + t334 * Ifges(5,3) - t560;
t316 = Ifges(5,6) * t334;
t170 = -Ifges(5,4) * t525 - t335 * Ifges(5,2) + t316;
t171 = -t334 * Ifges(4,2) - Ifges(4,6) * t525 + t569;
t209 = -mrSges(5,2) * t334 - mrSges(5,3) * t335;
t210 = t475 * t525 + t442;
t460 = -Ifges(4,2) * t405 + t567;
t238 = (t406 * Ifges(4,6) + t409 * t460) * t403;
t463 = Ifges(4,1) * t408 - t568;
t239 = (t406 * Ifges(4,5) + t409 * t463) * t403;
t454 = Ifges(5,3) * t405 - t558;
t240 = (t406 * Ifges(5,5) + t409 * t454) * t403;
t395 = mrSges(5,3) * t525;
t545 = t334 * mrSges(5,1);
t258 = t395 + t545;
t259 = t335 * mrSges(5,1) - mrSges(5,2) * t525;
t260 = mrSges(4,2) * t525 - t334 * mrSges(4,3);
t261 = -mrSges(4,1) * t525 - t335 * mrSges(4,3);
t464 = -mrSges(5,2) * t405 - mrSges(5,3) * t408;
t290 = t464 * t525;
t465 = mrSges(4,1) * t405 + mrSges(4,2) * t408;
t291 = t465 * t525;
t311 = (-mrSges(4,2) * t406 - mrSges(4,3) * t517) * t403;
t312 = (mrSges(4,1) * t406 - mrSges(4,3) * t510) * t403;
t396 = Ifges(3,5) * t525;
t433 = t409 * (Ifges(4,5) * t408 - Ifges(4,6) * t405);
t434 = t409 * (-Ifges(5,4) * t408 + Ifges(5,5) * t405);
t485 = t525 / 0.2e1;
t468 = t408 * t485;
t486 = -t525 / 0.2e1;
t469 = t408 * t486;
t470 = t405 * t485;
t471 = t405 * t486;
t487 = t526 / 0.2e1;
t3 = (0.2e1 * Ifges(3,4) * t525 + (Ifges(3,1) - Ifges(3,2)) * t526) * t485 + (-Ifges(3,6) * t526 + t396 / 0.2e1 + Ifges(3,5) * t485 + t638) * t531 + t238 * t607 + t78 * t610 + t335 * t611 + t122 * t612 + t121 * t613 + t288 * t624 + t171 * t471 + t170 * t469 + t169 * t470 + (Ifges(3,2) * t409 + t570) * t526 * t651 + (t239 + t120) * t603 + (t334 * t648 - t335 * t654 - t525 * t650) * t487 + t338 * (mrSges(4,1) * t334 + mrSges(4,2) * t335) + t334 * t240 / 0.2e1 + t158 * t311 - t157 * t312 + t129 * t313 + t130 * t314 + t127 * t290 + t295 * t291 + t162 * t258 + t164 * t259 + t201 * t260 + t200 * t261 + t210 * t209 + t39 * t213 + t38 * t214 + t90 * t159 + t58 * t146 + t57 * t147 + t128 * t108 + t647 * t468 + ((t406 * t650 + t433 + t434) * t582 - pkin(1) * (mrSges(3,1) * t406 + mrSges(3,2) * t409) + t406 * (Ifges(3,1) * t409 - t570) / 0.2e1) * t403 ^ 2 + m(4) * (-t157 * t200 + t158 * t201 + t295 * t338) + m(6) * (t128 * t90 + t38 * t57 + t39 * t58) + m(5) * (t127 * t210 + t129 * t162 + t130 * t164);
t546 = t3 * qJD(1);
t458 = Ifges(6,2) * t407 + t566;
t117 = -t334 * Ifges(6,6) + t335 * t458;
t118 = -t334 * Ifges(6,5) + t335 * t461;
t536 = t407 * mrSges(6,1);
t540 = t404 * mrSges(6,2);
t372 = t536 - t540;
t177 = t372 * t335;
t527 = t335 * t404;
t194 = -mrSges(6,1) * t334 - mrSges(6,3) * t527;
t534 = t407 * mrSges(6,3);
t195 = mrSges(6,2) * t334 + t335 * t534;
t530 = qJ(4) * t334;
t208 = pkin(3) * t335 + t530;
t439 = -t562 / 0.2e1 - t554 / 0.2e1;
t103 = t158 - t579;
t133 = t335 * t622 + t530;
t51 = t103 * t407 - t133 * t404;
t52 = t103 * t404 + t133 * t407;
t4 = t117 * t613 + t118 * t612 + t208 * t209 + t38 * t194 + t39 * t195 - t90 * t177 + t52 * t146 + t51 * t147 + t102 * t108 + (t259 - t261) * t158 + (t258 - t260) * t157 + m(6) * (t102 * t90 + t38 * t51 + t39 * t52) + m(5) * (t127 * t208 + t129 * t157 + t130 * t158) + (-t295 * mrSges(4,2) + t127 * mrSges(5,3) + t321 / 0.2e1 + t316 / 0.2e1 + t170 / 0.2e1 - t172 / 0.2e1 - t77 / 0.2e1 - t130 * mrSges(5,1) - t157 * mrSges(4,3)) * t334 + (t295 * mrSges(4,1) - t171 / 0.2e1 - t127 * mrSges(5,2) + t169 / 0.2e1 + t78 * t586 + t79 * t592 - t158 * mrSges(4,3) + t129 * mrSges(5,1) + (-Ifges(5,6) / 0.2e1 - Ifges(4,4) / 0.2e1 - t439) * t335 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t334) * t335 + t640 * t403 * t582;
t542 = t4 * qJD(1);
t539 = t404 * mrSges(6,3);
t107 = -mrSges(6,1) * t237 + mrSges(6,2) * t236;
t109 = Ifges(6,5) * t236 + Ifges(6,6) * t237;
t110 = Ifges(6,2) * t237 + t233;
t111 = Ifges(6,1) * t236 + t547;
t7 = t109 * t603 + t38 * t146 - t39 * t147 + t90 * t107 + (t575 - t111 / 0.2e1 + t78 / 0.2e1) * t237 + (-t576 + t624 + t110 / 0.2e1) * t236;
t532 = t7 * qJD(1);
t451 = t38 * t404 - t39 * t407;
t516 = t407 * t146;
t524 = t404 * t147;
t12 = (-t108 + t258) * t525 + (-t209 - t516 + t524) * t335 + m(6) * (t335 * t451 - t525 * t90) + m(5) * (-t127 * t335 + t129 * t525);
t528 = qJD(1) * t12;
t522 = t404 * t306;
t521 = t404 * t359;
t520 = t404 * t405;
t518 = t405 * t407;
t305 = -t408 * t458 + t555;
t514 = t407 * t305;
t513 = t407 * t360;
t501 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t500 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t496 = t561 / 0.2e1;
t495 = t554 / 0.2e1;
t494 = -t539 / 0.2e1;
t493 = t539 / 0.2e1;
t491 = -t534 / 0.2e1;
t490 = -t110 / 0.4e1 - t79 / 0.4e1;
t489 = -t111 / 0.4e1 + t625;
t339 = t372 * t408;
t484 = t339 * t582;
t342 = t408 * t380;
t480 = t305 / 0.4e1 + t342 / 0.4e1;
t341 = t408 * t378;
t479 = t306 / 0.4e1 - t341 / 0.4e1;
t478 = t461 / 0.4e1 + t596;
t477 = t595 - t458 / 0.4e1;
t467 = t146 / 0.2e1 - t550 / 0.2e1;
t466 = t548 / 0.2e1 + t618;
t459 = Ifges(4,2) * t408 + t568;
t455 = -Ifges(5,2) * t405 - t558;
t453 = -Ifges(5,3) * t408 - t559;
t452 = pkin(3) * t405 - t529;
t450 = t404 * t52 + t407 * t51;
t432 = t372 * t405;
t440 = -mrSges(6,2) * t408 + mrSges(6,3) * t518;
t441 = mrSges(6,1) * t408 - mrSges(6,3) * t520;
t411 = -(Ifges(5,3) * t335 + t170 + t316) * t408 / 0.4e1 - (Ifges(5,2) * t334 + t171 + t560) * t405 / 0.4e1 - (t456 + t459) * t335 / 0.4e1 + (t455 + t454) * t334 / 0.4e1 + t261 * t502 + t260 * t503 + t146 * t614 + t195 * t615 + t194 * t616 - t432 * t623 + t518 * t625 - t177 * t594 + t52 * t598 + t51 * t599 - t118 * t519 / 0.4e1 + t79 * t520 / 0.4e1 - t117 * t512 / 0.4e1 + t643 * mrSges(4,3) - (t460 + t644) * t334 / 0.4e1 + t259 * t577 / 0.2e1 + t258 * t578 / 0.2e1 + (t405 * t457 + t453 + t463 + t514 + t522 + t552) * t602 + (t643 + t646 * t589 + (-t157 / 0.2e1 + t130 / 0.2e1) * t408) * mrSges(5,1) + (-Ifges(4,1) * t334 + t335 * t457 + t169 - t553 - t569) * t588 - t381 * t108 / 0.2e1 + t208 * t373 / 0.2e1 + t369 * (-mrSges(5,2) * t335 + mrSges(5,3) * t334) / 0.2e1 + t102 * t339 / 0.2e1 - pkin(2) * (mrSges(4,1) * t335 - mrSges(4,2) * t334) / 0.2e1 + t295 * t465 / 0.2e1 - t237 * (Ifges(6,5) * t408 + t405 * t461) / 0.4e1 + t127 * t464 / 0.2e1 + t452 * t209 / 0.2e1 + t236 * (Ifges(6,6) * t408 + t405 * t458) / 0.4e1 + t234 * t147 / 0.2e1 + t39 * t440 / 0.2e1 + t38 * t441 / 0.2e1 + (t452 * t127 + t369 * t208 + ((t130 - t157) * t408 + t646 * t405) * pkin(8)) * t633 + (-Ifges(4,2) * t335 - t321 + t647) * t408 / 0.4e1 + (t102 * t382 + t229 * t51 + t230 * t52 + t234 * t38 + t235 * t39 - t381 * t90) * t631 + (-t434 / 0.4e1 - t433 / 0.4e1) * t403;
t1 = Ifges(5,4) * t469 + Ifges(4,5) * t468 + Ifges(5,5) * t470 + Ifges(4,6) * t471 + t404 * t641 + t407 * t642 + t487 * t650 + t491 * t57 + t494 * t58 + t499 * t597 - t411 + t637;
t15 = -t235 * t360 - t234 * t359 - t369 * t464 - t230 * t440 - t229 * t441 + t381 * t339 + t382 * t432 - m(6) * (t229 * t234 + t230 * t235 - t381 * t382) + (pkin(2) * mrSges(4,1) - t514 / 0.2e1 - t522 / 0.2e1 + (t439 + t574) * t405) * t405 + (pkin(2) * mrSges(4,2) + (t457 - t574) * t408 + (Ifges(4,2) - Ifges(5,2) + Ifges(5,3) + Ifges(6,2) * t635 / 0.2e1 + (t565 + t571 / 0.2e1) * t404 - t655) * t405) * t408 + t645 * t452;
t448 = -t1 * qJD(1) - t15 * qJD(2);
t400 = Ifges(6,6) * t519;
t23 = t400 * t591 + t382 * t340 + ((-t342 / 0.2e1 - t305 / 0.2e1) * t404 + (t609 - t341 / 0.2e1 + t496) * t407) * t408 + (-t505 + t359) * t230 + (-t360 - t504) * t229;
t414 = (t229 * t628 + t479) * t236 + (mrSges(6,3) * t615 + t480) * t237 + t146 * t616 + t230 * t618 + t400 * t602 + t38 * t598 + t107 * t594 + t39 * t600 + t109 * t588 + t90 * t601;
t417 = (t575 / 0.2e1 + t489) * t404 + (-t563 / 0.4e1 + t576 / 0.2e1 + t490) * t407;
t423 = -t564 / 0.2e1 - t557 / 0.2e1 + mrSges(6,1) * t626 + t58 * t629;
t5 = (Ifges(6,3) * t486 + t417) * t408 + t414 + t423;
t447 = t5 * qJD(1) - t23 * qJD(2);
t436 = t521 / 0.2e1 - t513 / 0.2e1;
t437 = t524 / 0.2e1 - t516 / 0.2e1;
t444 = t229 * t404 - t230 * t407;
t412 = (-t209 / 0.2e1 + t437) * t405 + (-t373 / 0.2e1 + t436) * t335 + (-pkin(8) * t499 - t127 * t405 - t369 * t335) * t633 + (t335 * t444 - t382 * t525 + t405 * t451) * t631;
t418 = t164 * t634 + t213 * t593 + t214 * t587 + t449 * t632;
t10 = -t370 + (t484 - t406 * mrSges(5,2) / 0.2e1) * t403 + t412 + t418;
t42 = (m(6) * t444 - t513 + t521 + t645) * t405;
t446 = qJD(1) * t10 + qJD(2) * t42;
t18 = (t543 / 0.2e1 - t467) * t407 + (t544 / 0.2e1 - t466) * t404;
t45 = (t537 / 0.2e1 - t360 / 0.2e1 + t408 * t491) * t407 + (t538 / 0.2e1 + t599 + t408 * t494) * t404;
t445 = qJD(1) * t18 + qJD(2) * t45;
t438 = t535 / 0.2e1 + t541 / 0.2e1;
t435 = t378 * t586 + t380 * t592;
t100 = -qJ(4) * t372 + (t461 / 0.2e1 + t378 / 0.2e1) * t407 + (t380 / 0.2e1 - t458 / 0.2e1) * t404;
t420 = t360 * t581 + (t534 * t581 - t477) * t408 - t480;
t421 = -t622 * t600 + (-t493 * t622 + t478) * t408 - t479;
t426 = qJ(4) * t601 + t372 * t594 - t457 * t588;
t428 = -t552 / 0.2e1 + t234 * t630 + mrSges(6,2) * t614;
t21 = (-t555 / 0.2e1 + t420) * t407 + (-t561 / 0.2e1 + t421) * t404 + t426 + t428;
t416 = t478 * t237 + t477 * t236 + qJ(4) * t107 / 0.2e1 - t457 * t602 + t372 * t623;
t424 = -t466 * t622 + t490;
t425 = -t467 * t622 - t489;
t429 = t553 / 0.2e1 + t51 * t630 + t52 * t629;
t9 = (-t556 / 0.2e1 + t425) * t407 + (-t563 / 0.2e1 + t424) * t404 + t416 + t429;
t431 = t9 * qJD(1) + t21 * qJD(2) - t100 * qJD(3);
t394 = -0.2e1 * t498;
t415 = -t395 + (t394 + t158) * t633 + (t103 + t394) * t631 - t551 / 0.2e1 - t549 / 0.2e1 + t374 * t486;
t419 = t158 * t634 + t194 * t587 + t195 * t593 + t450 * t632;
t17 = t415 + t419;
t333 = mrSges(5,3) + (m(5) + m(6)) * qJ(4) + t374;
t47 = 0.2e1 * (t382 / 0.4e1 - t515 / 0.4e1 - t523 / 0.4e1) * m(6);
t430 = qJD(1) * t17 + qJD(2) * t47 + qJD(3) * t333;
t46 = t405 * t438 - t436 + (t404 ^ 2 + t635) * mrSges(6,3) * t584;
t41 = m(6) * t594 - t639 * t631 + t440 * t592 + t441 * t586 + (mrSges(5,1) - t540 / 0.2e1 + t536 / 0.2e1 + m(5) * pkin(8)) * t408;
t20 = t405 * t495 + t407 * t420 + t426 - t428 + (t421 + t496) * t404;
t19 = t236 * t491 + t237 * t493 + t335 * t438 - t437;
t16 = t415 - t419 - t545;
t11 = mrSges(5,2) * t487 + t403 * t484 + t412 - t418;
t8 = t335 * t495 + t404 * t424 + t407 * t425 + t527 * t627 + t416 - t429;
t6 = Ifges(6,3) * t468 + t408 * t417 + t414 - t423;
t2 = t411 + (mrSges(6,3) * t626 + t642) * t407 + (t58 * t628 + t641) * t404 + ((Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t406 + (t500 * t405 + (t597 + t501) * t408) * t409) * t403 + t637;
t13 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t12 + qJD(5) * t7, t2 * qJD(3) + t11 * qJD(4) + t6 * qJD(5) + t546 + (-pkin(2) * t291 + t128 * t339 + t382 * t159 + t210 * t373 + t230 * t213 + t229 * t214 + t288 * t609 + t369 * t290 + t305 * t610 + t57 * t359 + t58 * t360 + t396 + t638 + (t338 * mrSges(4,2) + t239 / 0.2e1 + t611 + t120 / 0.2e1 + t164 * mrSges(5,1) - t200 * mrSges(4,3) + (-t312 + t314) * pkin(8)) * t405 + (-t338 * mrSges(4,1) + t238 / 0.2e1 - t240 / 0.2e1 + t121 * t587 + t122 * t593 - t162 * mrSges(5,1) + t201 * mrSges(4,3) + (t311 - t313) * pkin(8)) * t408 + (t210 * t369 + (-t162 * t408 + t164 * t405) * pkin(8)) * t652 + m(4) * (-pkin(2) * t338 + (-t200 * t405 + t201 * t408) * pkin(8)) + (t128 * t382 + t229 * t57 + t230 * t58) * t653 + ((t405 * t501 - t408 * t500 - Ifges(3,6)) * t406 + (t459 * t591 - t408 * t455 / 0.2e1 + t453 * t589 + t644 * t584) * t409) * t403) * qJD(2), t2 * qJD(2) + t16 * qJD(4) + t8 * qJD(5) + t542 + (pkin(3) * t545 - qJ(4) * t177 + t102 * t374 + t376 * t607 + (-t622 * t194 - t51 * mrSges(6,3) + t118 / 0.2e1) * t407 + (-t622 * t195 - t52 * mrSges(6,3) - t117 / 0.2e1) * t404 + (qJ(4) * t102 - t450 * t622) * t653 + (-qJ(4) * mrSges(5,1) + t435) * t335 + t640 + (-pkin(3) * t652 + t657) * t158 + (-qJ(4) * t652 + t656) * t157) * qJD(3), qJD(2) * t11 + qJD(3) * t16 + qJD(5) * t19 + t528, t532 + t6 * qJD(2) + t8 * qJD(3) + t19 * qJD(4) + (-mrSges(6,1) * t39 - mrSges(6,2) * t38 + t109) * qJD(5); -qJD(3) * t1 + qJD(4) * t10 + qJD(5) * t5 - t546, -qJD(3) * t15 + qJD(4) * t42 - qJD(5) * t23, t41 * qJD(4) + t20 * qJD(5) + t448 + (-t381 * t374 + m(6) * (-qJ(4) * t381 + t622 * t639) + t639 * mrSges(6,3) + (-pkin(3) * mrSges(5,1) + t376 / 0.2e1 + (-t573 + t627) * t407 + (t572 - Ifges(6,6) / 0.2e1) * t404 + (-m(5) * pkin(3) + t657) * pkin(8) - t654) * t408 + (t461 * t586 + t458 * t593 + (-mrSges(5,1) - t372) * qJ(4) + (-m(5) * qJ(4) + t656) * pkin(8) + t435 + t648) * t405) * qJD(3), qJD(3) * t41 + qJD(5) * t46 + t446, t20 * qJD(3) + t46 * qJD(4) + (-mrSges(6,1) * t230 - mrSges(6,2) * t229 - Ifges(6,5) * t512 + t400) * qJD(5) + t447; qJD(2) * t1 + qJD(4) * t17 + qJD(5) * t9 - t542, qJD(4) * t47 + qJD(5) * t21 - t448, qJD(4) * t333 - qJD(5) * t100, t430, ((-Ifges(6,6) + t572) * t407 + (-Ifges(6,5) + t573) * t404) * qJD(5) + t431; -qJD(2) * t10 - qJD(3) * t17 - qJD(5) * t18 - t528, -qJD(3) * t47 - qJD(5) * t45 - t446, -t430, 0, -qJD(5) * t374 - t445; -qJD(2) * t5 - qJD(3) * t9 + qJD(4) * t18 - t532, -qJD(3) * t21 + qJD(4) * t45 - t447, -t431, t445, 0;];
Cq = t13;
