% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR14_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:17:05
% EndTime: 2019-12-31 19:17:27
% DurationCPUTime: 11.01s
% Computational Cost: add. (35162->700), mult. (97430->1041), div. (0->0), fcn. (109677->12), ass. (0->366)
t430 = sin(qJ(3));
t427 = sin(pkin(6));
t539 = sin(pkin(5));
t540 = cos(pkin(11));
t476 = t540 * t539;
t541 = cos(pkin(6));
t542 = cos(pkin(5));
t449 = t542 * t427 + t541 * t476;
t426 = sin(pkin(11));
t488 = t426 * t539;
t586 = cos(qJ(3));
t308 = t449 * t430 + t586 * t488;
t429 = sin(qJ(4));
t432 = cos(qJ(4));
t448 = t427 * t476 - t542 * t541;
t260 = t308 * t429 + t448 * t432;
t428 = sin(qJ(5));
t431 = cos(qJ(5));
t652 = Ifges(6,5) * t431 - Ifges(6,6) * t428;
t661 = t260 * t652;
t660 = -Ifges(6,5) / 0.2e1;
t659 = -Ifges(6,3) / 0.2e1;
t475 = mrSges(6,1) * t431 - mrSges(6,2) * t428;
t658 = mrSges(5,1) + t475;
t657 = t661 / 0.4e1;
t477 = t541 * t539;
t466 = t426 * t477;
t337 = -t430 * t466 + t586 * t476;
t478 = t427 * t488;
t280 = t337 * t429 - t432 * t478;
t624 = -t280 / 0.2e1;
t398 = Ifges(6,5) * t428 + Ifges(6,6) * t431;
t656 = -t398 / 0.4e1;
t655 = t448 / 0.2e1;
t307 = t430 * t488 - t449 * t586;
t530 = t307 * t429;
t654 = t530 * t659;
t261 = t308 * t432 - t448 * t429;
t558 = t261 * mrSges(5,3);
t205 = -t261 * t428 + t307 * t431;
t206 = t261 * t431 + t307 * t428;
t101 = -mrSges(6,1) * t205 + mrSges(6,2) * t206;
t208 = mrSges(5,1) * t307 - t558;
t653 = t101 - t208;
t420 = Ifges(6,4) * t431;
t651 = -Ifges(6,2) * t428 + t420;
t403 = Ifges(6,1) * t428 + t420;
t650 = -t558 / 0.2e1 - t208 / 0.2e1;
t409 = t429 * pkin(4) - pkin(10) * t432;
t523 = t428 * t429;
t349 = pkin(9) * t523 + t409 * t431;
t521 = t429 * t431;
t350 = -pkin(9) * t521 + t409 * t428;
t649 = -t349 * t428 + t350 * t431;
t421 = Ifges(5,4) * t432;
t647 = Ifges(5,2) * t429 - t421;
t646 = m(6) * pkin(10) + mrSges(6,3);
t522 = t428 * t432;
t227 = t307 * t522 + t308 * t431;
t520 = t431 * t432;
t228 = -t307 * t520 + t308 * t428;
t645 = -Ifges(6,6) * t227 / 0.2e1 + t228 * t660;
t644 = -m(6) * pkin(4) - t658;
t643 = -m(5) / 0.2e1;
t642 = m(5) / 0.2e1;
t641 = -m(6) / 0.2e1;
t640 = m(6) / 0.2e1;
t639 = -mrSges(6,1) / 0.2e1;
t638 = mrSges(5,2) / 0.2e1;
t637 = -mrSges(6,2) / 0.2e1;
t636 = mrSges(6,2) / 0.2e1;
t634 = -mrSges(6,3) / 0.2e1;
t504 = pkin(1) * t542;
t415 = t540 * t504;
t485 = t539 * qJ(2);
t309 = t542 * pkin(2) + t415 + (-pkin(8) * t477 - t485) * t426;
t333 = -t539 * pkin(1) - pkin(2) * t476 - pkin(8) * t478;
t249 = -t309 * t427 + t541 * t333;
t151 = pkin(3) * t307 - pkin(9) * t308 + t249;
t516 = qJ(2) * t476 + t426 * t504;
t291 = t449 * pkin(8) + t516;
t486 = t541 * t309;
t194 = t586 * t291 + (t333 * t427 + t486) * t430;
t156 = -t448 * pkin(9) + t194;
t65 = t151 * t432 - t429 * t156;
t60 = -t307 * pkin(4) - t65;
t633 = t60 / 0.2e1;
t203 = Ifges(6,4) * t205;
t573 = Ifges(6,5) * t260;
t76 = Ifges(6,1) * t206 + t203 + t573;
t632 = t76 / 0.2e1;
t560 = t206 * mrSges(6,3);
t136 = mrSges(6,1) * t260 - t560;
t631 = -t136 / 0.2e1;
t259 = Ifges(5,4) * t260;
t557 = t261 * Ifges(5,1);
t576 = Ifges(5,5) * t307;
t144 = -t259 + t557 + t576;
t630 = -t144 / 0.2e1;
t629 = t205 / 0.2e1;
t628 = t206 / 0.2e1;
t625 = t261 / 0.2e1;
t623 = t307 / 0.2e1;
t524 = t427 * t430;
t372 = t429 * t541 + t432 * t524;
t503 = t427 * t586;
t310 = -t428 * t372 - t431 * t503;
t622 = t310 / 0.2e1;
t311 = t372 * t431 - t428 * t503;
t621 = -t311 / 0.2e1;
t620 = t311 / 0.2e1;
t501 = t432 * t586;
t344 = (-t428 * t501 + t430 * t431) * t427;
t619 = t344 / 0.2e1;
t392 = -pkin(4) * t432 - t429 * pkin(10) - pkin(3);
t346 = -pkin(9) * t522 + t431 * t392;
t618 = t346 / 0.2e1;
t617 = t350 / 0.2e1;
t564 = Ifges(6,6) * t432;
t356 = t429 * t651 - t564;
t616 = t356 / 0.2e1;
t577 = Ifges(6,4) * t428;
t404 = Ifges(6,1) * t431 - t577;
t572 = Ifges(6,5) * t432;
t358 = t429 * t404 - t572;
t615 = t358 / 0.2e1;
t371 = t429 * t524 - t432 * t541;
t614 = t371 / 0.2e1;
t613 = t372 / 0.2e1;
t374 = t475 * t429;
t612 = -t374 / 0.2e1;
t546 = t431 * mrSges(6,2);
t548 = t428 * mrSges(6,1);
t396 = t546 + t548;
t375 = t396 * t429;
t611 = t375 / 0.2e1;
t376 = t396 * t432;
t610 = t376 / 0.2e1;
t514 = mrSges(6,3) * t523;
t386 = mrSges(6,2) * t432 - t514;
t609 = -t386 / 0.2e1;
t608 = t386 / 0.2e1;
t387 = -t429 * mrSges(6,2) - mrSges(6,3) * t522;
t607 = t387 / 0.2e1;
t388 = -mrSges(6,1) * t432 - mrSges(6,3) * t521;
t606 = -t388 / 0.2e1;
t605 = t388 / 0.2e1;
t389 = t429 * mrSges(6,1) - mrSges(6,3) * t520;
t604 = t389 / 0.2e1;
t603 = t475 / 0.2e1;
t602 = t396 / 0.2e1;
t601 = t398 / 0.2e1;
t400 = Ifges(6,2) * t431 + t577;
t600 = -t400 / 0.4e1;
t599 = -t428 / 0.2e1;
t598 = -t428 / 0.4e1;
t597 = t428 / 0.2e1;
t596 = t428 / 0.4e1;
t595 = -t429 / 0.2e1;
t594 = t429 / 0.2e1;
t593 = -t431 / 0.2e1;
t592 = -t431 / 0.4e1;
t591 = t431 / 0.2e1;
t590 = t431 / 0.4e1;
t589 = -t432 / 0.2e1;
t588 = -t432 / 0.4e1;
t587 = t432 / 0.2e1;
t585 = pkin(9) * t429;
t584 = pkin(9) * t432;
t583 = pkin(10) * t428;
t582 = pkin(10) * t431;
t581 = mrSges(5,2) * t432;
t579 = Ifges(5,4) * t429;
t578 = Ifges(6,4) * t206;
t575 = Ifges(5,5) * t308;
t569 = Ifges(5,6) * t307;
t568 = Ifges(5,6) * t308;
t566 = Ifges(6,6) * t260;
t563 = Ifges(6,3) * t261;
t562 = Ifges(6,3) * t429;
t561 = t205 * mrSges(6,3);
t559 = t260 * mrSges(5,3);
t556 = t261 * Ifges(5,4);
t109 = Ifges(6,4) * t228 + Ifges(6,2) * t227 - Ifges(6,6) * t530;
t110 = Ifges(6,1) * t228 + Ifges(6,4) * t227 - Ifges(6,5) * t530;
t193 = -t430 * t291 + t333 * t503 + t586 * t486;
t246 = pkin(3) * t308 + pkin(9) * t307;
t118 = -t429 * t193 + t246 * t432;
t119 = t432 * t193 + t429 * t246;
t135 = -mrSges(6,2) * t260 + t561;
t140 = -mrSges(6,1) * t227 + mrSges(6,2) * t228;
t155 = t448 * pkin(3) - t193;
t182 = mrSges(6,2) * t530 + mrSges(6,3) * t227;
t183 = -mrSges(6,1) * t530 - mrSges(6,3) * t228;
t474 = Ifges(5,1) * t432 - t579;
t192 = -t474 * t307 + t575;
t207 = -mrSges(5,2) * t307 - t559;
t397 = t429 * mrSges(5,1) + t581;
t230 = t397 * t307;
t236 = -mrSges(5,2) * t308 + mrSges(5,3) * t530;
t237 = mrSges(5,3) * t307 * t432 + t308 * mrSges(5,1);
t245 = mrSges(4,1) * t308 - mrSges(4,2) * t307;
t554 = t307 * mrSges(4,3);
t264 = t448 * mrSges(4,2) - t554;
t66 = t429 * t151 + t156 * t432;
t61 = t307 * pkin(10) + t66;
t77 = t260 * pkin(4) - t261 * pkin(10) + t155;
t36 = -t428 * t61 + t431 * t77;
t37 = t428 * t77 + t431 * t61;
t493 = t654 - t568 / 0.2e1 - t647 * t307 / 0.2e1 - t645;
t129 = -t409 * t307 + t194;
t83 = pkin(10) * t308 + t119;
t50 = t129 * t431 - t428 * t83;
t143 = -t260 * Ifges(5,2) + t556 + t569;
t74 = Ifges(6,5) * t206 + Ifges(6,6) * t205 + t260 * Ifges(6,3);
t506 = t143 / 0.2e1 - t74 / 0.2e1;
t51 = t129 * t428 + t431 * t83;
t511 = t569 / 0.2e1;
t513 = -t576 / 0.2e1;
t517 = -Ifges(4,5) * t307 - Ifges(4,6) * t308;
t171 = mrSges(5,1) * t260 + mrSges(5,2) * t261;
t553 = t308 * mrSges(4,3);
t265 = -t448 * mrSges(4,1) - t553;
t519 = t171 - t265;
t75 = Ifges(6,2) * t205 + t566 + t578;
t82 = -t308 * pkin(4) - t118;
t3 = t36 * t183 + t37 * t182 + m(5) * (t118 * t65 + t119 * t66 + t155 * t194) + m(6) * (t36 * t50 + t37 * t51 + t60 * t82) + (Ifges(5,5) * t625 - Ifges(5,6) * t260 / 0.2e1 + Ifges(4,6) * t655 - t194 * mrSges(4,3) - Ifges(4,4) * t308) * t308 + (Ifges(4,4) * t307 + Ifges(4,5) * t655 + t193 * mrSges(4,3) + (t630 + t513) * t432 + (t511 + t506) * t429 + (-Ifges(4,1) + Ifges(4,2) + Ifges(5,3)) * t308) * t307 + t50 * t136 + t60 * t140 + t51 * t135 + t82 * t101 + t519 * t194 - t448 * t517 / 0.2e1 + t493 * t260 + t192 * t625 + t110 * t628 + t109 * t629 + t228 * t632 + t119 * t207 + t118 * t208 + t227 * t75 / 0.2e1 - t155 * t230 + t66 * t236 + t65 * t237 + t249 * t245 + t193 * t264;
t555 = t3 * qJD(1);
t347 = pkin(9) * t520 + t428 * t392;
t552 = t347 * mrSges(6,3);
t551 = t36 * t428;
t550 = t37 * t431;
t120 = t563 - t661;
t121 = Ifges(6,6) * t261 - t260 * t651;
t122 = Ifges(6,5) * t261 - t404 * t260;
t159 = t396 * t260;
t534 = t260 * t428;
t165 = -mrSges(6,2) * t261 + mrSges(6,3) * t534;
t533 = t260 * t431;
t166 = mrSges(6,1) * t261 + mrSges(6,3) * t533;
t170 = mrSges(5,1) * t261 - mrSges(5,2) * t260;
t172 = -Ifges(5,2) * t261 - t259;
t173 = -Ifges(5,1) * t260 - t556;
t518 = -Ifges(5,5) * t260 - Ifges(5,6) * t261;
t178 = pkin(4) * t261 + pkin(10) * t260;
t54 = t178 * t431 - t428 * t65;
t545 = t431 * t76;
t547 = t428 * t75;
t55 = t178 * t428 + t431 * t65;
t4 = t122 * t628 + t121 * t629 + t55 * t135 + t37 * t165 + t518 * t623 + t65 * t207 + t155 * t170 + t54 * t136 + t36 * t166 - t60 * t159 + m(6) * (t36 * t54 + t37 * t55) + (t173 / 0.2e1 - t506) * t261 + (t65 * mrSges(5,3) - t545 / 0.2e1 + t547 / 0.2e1 + t630 - t172 / 0.2e1 + t120 / 0.2e1) * t260 + (m(6) * t60 - t558 + t653) * t66;
t549 = t4 * qJD(1);
t544 = t432 * mrSges(5,1);
t100 = mrSges(6,1) * t206 + mrSges(6,2) * t205;
t102 = Ifges(6,5) * t205 - Ifges(6,6) * t206;
t103 = -Ifges(6,2) * t206 + t203;
t104 = Ifges(6,1) * t205 - t578;
t7 = -t37 * t136 + t60 * t100 + t260 * t102 / 0.2e1 + t36 * t135 + (t104 / 0.2e1 - t75 / 0.2e1 - t37 * mrSges(6,3)) * t206 + (-t36 * mrSges(6,3) + t632 + t103 / 0.2e1) * t205;
t543 = t7 * qJD(1);
t281 = t432 * t337 + t429 * t478;
t336 = t430 * t476 + t586 * t466;
t238 = -t281 * t428 + t336 * t431;
t239 = t281 * t431 + t336 * t428;
t14 = t238 * t136 + t239 * t135 + t281 * t207 + t337 * t264 + (mrSges(4,1) * t307 + mrSges(4,2) * t308) * t478 + (-t542 * mrSges(3,2) + mrSges(3,3) * t476) * t476 - (t542 * mrSges(3,1) - mrSges(3,3) * t488) * t488 + t519 * t336 + t653 * t280 + m(6) * (t238 * t36 + t239 * t37 + t280 * t60) + m(5) * (t155 * t336 - t280 * t65 + t281 * t66) + m(4) * (-t193 * t336 + t194 * t337 + t249 * t478) + m(3) * (t516 * t476 - (-t426 * t485 + t415) * t488);
t538 = t14 * qJD(1);
t446 = (-t310 * t205 / 0.2e1 + t206 * t621) * mrSges(6,3) + t135 * t622 + t136 * t621 + t100 * t614;
t464 = t238 * mrSges(6,1) / 0.2e1 + t239 * t637;
t15 = t446 - t464;
t537 = t15 * qJD(1);
t536 = t238 * t428;
t535 = t239 * t431;
t532 = t280 * t429;
t529 = t310 * t428;
t528 = t311 * t431;
t525 = t371 * t280;
t512 = -t572 / 0.2e1;
t509 = mrSges(6,3) * t599;
t508 = mrSges(6,3) * t591;
t507 = t104 / 0.4e1 - t75 / 0.4e1;
t505 = t76 / 0.4e1 + t103 / 0.4e1;
t502 = t429 * t586;
t500 = t586 * t336;
t498 = t523 / 0.2e1;
t497 = -t522 / 0.2e1;
t496 = -t521 / 0.2e1;
t495 = t520 / 0.2e1;
t405 = t429 * Ifges(5,1) + t421;
t494 = t405 * t589;
t378 = t403 * t429;
t492 = -t356 / 0.4e1 - t378 / 0.4e1;
t377 = t429 * t400;
t491 = t358 / 0.4e1 - t377 / 0.4e1;
t490 = t600 + t404 / 0.4e1;
t489 = t651 / 0.4e1 + t403 / 0.4e1;
t483 = t429 ^ 2 * t503;
t482 = t427 * t502;
t481 = t371 * t502;
t480 = -t503 / 0.2e1;
t479 = t503 / 0.2e1;
t345 = (t428 * t430 + t431 * t501) * t427;
t467 = t429 * t479;
t468 = t429 * t480;
t433 = (-t371 * t118 + t372 * t119 + (t155 * t430 - t586 * t194 + t66 * t501 - t65 * t502) * t427) * t642 + (t310 * t50 + t311 * t51 + t344 * t36 + t345 * t37 + t371 * t82 + t60 * t482) * t640 + t183 * t622 + t182 * t620 + t136 * t619 + t345 * t135 / 0.2e1 + t140 * t614 - t371 * t237 / 0.2e1 + t236 * t613 + t541 * t245 / 0.2e1 + t171 * t524 / 0.2e1 - t230 * t480 + t101 * t467 + t208 * t468 - (t265 + t553) * t524 / 0.2e1 + (t432 * t207 + t264 + t554) * t479;
t395 = t429 * mrSges(5,2) - t544;
t435 = (pkin(9) * t532 + t238 * t346 + t239 * t347) * t641 + t238 * t606 + t239 * t609 + t375 * t624 + t337 * mrSges(4,2) / 0.2e1 + (-pkin(3) * t643 + mrSges(4,1) / 0.2e1 - t395 / 0.2e1) * t336 + (pkin(9) * t643 - mrSges(5,3) / 0.2e1) * (t281 * t432 + t532);
t9 = t433 + t435;
t99 = m(5) * (t372 * t501 - t430 * t503 + t481) * t427 + m(6) * (t310 * t344 + t311 * t345 + t427 * t481);
t470 = t9 * qJD(1) + t99 * qJD(2);
t457 = -t159 / 0.2e1 - t207 / 0.2e1 - t559 / 0.2e1;
t437 = (t136 * t597 + (-t550 + t66 + t551) * t640 + t135 * t593 + t457) * t371 + (t310 * t54 + t311 * t55 + t372 * t60) * t640 + t166 * t622 + t165 * t620 + t170 * t480;
t447 = (-pkin(4) * t280 + (t535 - t536) * pkin(10)) * t641 + t281 * t638;
t458 = t101 / 0.2e1 + t650;
t11 = (t603 + mrSges(5,1) / 0.2e1) * t280 + (-t535 / 0.2e1 + t536 / 0.2e1) * mrSges(6,3) + t458 * t372 + t437 + t447;
t89 = m(6) * (t372 - t528 + t529) * t371;
t469 = t11 * qJD(1) + t89 * qJD(2);
t465 = pkin(4) * t612 + t588 * t652;
t463 = mrSges(6,1) * t619 + t345 * t637;
t462 = mrSges(6,2) * t617 + t349 * t639;
t461 = t546 / 0.2e1 + t548 / 0.2e1;
t460 = pkin(10) * t606 + t491;
t459 = pkin(10) * t609 + t492;
t355 = -t432 * Ifges(6,3) + t429 * t652;
t357 = Ifges(6,6) * t429 + t432 * t651;
t359 = Ifges(6,5) * t429 + t404 * t432;
t402 = t432 * Ifges(5,2) + t579;
t419 = Ifges(5,5) * t432;
t434 = (t355 / 0.4e1 - t402 / 0.4e1) * t261 + (-t421 / 0.4e1 + t358 * t592 + t356 * t596 - t405 / 0.4e1) * t260 + (t346 * t54 + t347 * t55 + t349 * t36 + t350 * t37) * t640 - pkin(3) * t170 / 0.2e1 + t155 * t397 / 0.2e1 + t205 * t357 / 0.4e1 + t206 * t359 / 0.4e1 + t307 * t419 / 0.4e1 + t166 * t618 + t347 * t165 / 0.2e1 + t349 * t136 / 0.2e1 + t135 * t617 + t36 * t604 + t37 * t607 + t54 * t605 + t55 * t608 + t60 * t610 + t66 * t611;
t440 = (-pkin(4) * t82 + (-t50 * t428 + t51 * t431) * pkin(10)) * t641 - Ifges(5,3) * t308 / 0.2e1 + pkin(4) * t140 / 0.2e1 - t118 * mrSges(5,1) / 0.2e1 + t119 * t638 + t227 * t600 - t228 * t403 / 0.4e1 + t82 * t603;
t442 = t545 / 0.4e1 - t547 / 0.4e1 + t557 / 0.4e1 + t657 + t172 / 0.4e1 + t144 / 0.4e1 - t120 / 0.4e1;
t445 = (Ifges(5,2) / 0.4e1 + Ifges(6,3) / 0.4e1) * t260 - t143 / 0.4e1 + t173 / 0.4e1 + t74 / 0.4e1 - t556 / 0.4e1 + t121 * t598 + t122 * t590;
t1 = t434 + ((-0.3e1 / 0.4e1 * Ifges(5,6) + t398 / 0.4e1) * t307 + (t66 * t640 + t457) * pkin(9) + t445) * t429 + (t576 / 0.2e1 + (m(6) * t633 + t458) * pkin(9) + t442) * t432 + (-pkin(10) * t182 / 0.2e1 + t51 * t634 - t109 / 0.4e1) * t431 + (t50 * mrSges(6,3) / 0.2e1 + pkin(10) * t183 / 0.2e1 - t110 / 0.4e1) * t428 + t440;
t436 = (t349 * t310 + t350 * t311 + t372 * t585) * t640 + t310 * t604 + t311 * t607 + t372 * t611 + t397 * t480 + ((t346 * t428 - t347 * t431 + t584) * t640 + t610 + t388 * t597 + t386 * t593) * t371;
t438 = (-pkin(4) * t482 + (-t344 * t428 + t345 * t431) * pkin(10)) * t640 + t344 * t509 + t345 * t508 + mrSges(5,1) * t468 - t475 * t467 + t480 * t581;
t26 = -t436 + t438;
t39 = t522 * t616 + t357 * t498 - t375 * t584 - t376 * t585 - t358 * t520 / 0.2e1 + t359 * t496 - t350 * t386 - t347 * t387 - t349 * t388 - t346 * t389 + pkin(3) * t397 + t494 + t402 * t594 - t647 * t589 - m(6) * (pkin(9) ^ 2 * t429 * t432 + t346 * t349 + t347 * t350) + (t432 * t652 + t562) * t587 + (t474 + t355) * t595;
t455 = t1 * qJD(1) - t26 * qJD(2) - t39 * qJD(3);
t452 = t310 * t609 + t311 * t605 + t371 * t612;
t40 = (t528 / 0.2e1 - t529 / 0.2e1) * t429 * mrSges(6,3) + t452 + t463;
t439 = (t346 * t634 + t491) * t205 + (-t552 / 0.2e1 + t492) * t206 + t135 * t618 + t347 * t631 + t36 * t608 + t37 * t606 + t102 * t588 + t374 * t633;
t441 = t104 * t590 + t75 * t592 + pkin(9) * t100 / 0.2e1 + t260 * t656 + (t551 / 0.2e1 - t550 / 0.2e1) * mrSges(6,3) + (t76 + t103) * t598;
t450 = t50 * t639 + t51 * t636 + t645;
t5 = (Ifges(6,3) * t623 + t441) * t429 + t439 + t450;
t56 = -t374 * t585 + t347 * t388 + ((-t377 / 0.2e1 + t615 + t512) * t428 + (t616 + t378 / 0.2e1 + t552 - t564 / 0.2e1) * t431) * t429 + (-t386 - t514) * t346;
t454 = t5 * qJD(1) - t40 * qJD(2) - t56 * qJD(3);
t453 = -t563 / 0.2e1 + t54 * t639 + t55 * t636;
t443 = t490 * t206 + t489 * t205 - pkin(4) * t100 / 0.2e1 + t657 + t60 * t602;
t13 = (t573 / 0.2e1 + (t631 - t560 / 0.2e1) * pkin(10) + t505) * t431 + (-t566 / 0.2e1 + (t561 / 0.2e1 - t135 / 0.2e1) * pkin(10) + t507) * t428 + t443 + t453;
t248 = -pkin(4) * t396 + (t403 / 0.2e1 + t651 / 0.2e1) * t431 + (t404 / 0.2e1 - t400 / 0.2e1) * t428;
t422 = t428 ^ 2;
t424 = t431 ^ 2;
t444 = pkin(9) * t602 + t490 * t431 - t489 * t428 + (-t424 / 0.2e1 - t422 / 0.2e1) * pkin(10) * mrSges(6,3);
t42 = (t512 + t460) * t431 + (t564 / 0.2e1 + t459) * t428 + (t659 + t444) * t429 + t462 + t465;
t86 = (-t396 / 0.2e1 + t461) * t371;
t451 = t13 * qJD(1) - t86 * qJD(2) + t42 * qJD(3) + t248 * qJD(4);
t425 = t432 ^ 2;
t393 = pkin(9) * t483;
t87 = (t461 + t602) * t371;
t43 = Ifges(6,6) * t497 + Ifges(6,5) * t495 + t562 / 0.2e1 + t460 * t431 + t459 * t428 + t444 * t429 - t462 + t465;
t41 = -t452 + t463 + (t310 * t498 + t311 * t496) * mrSges(6,3);
t27 = t436 + t438;
t16 = t446 + t464;
t12 = t533 * t660 + Ifges(6,6) * t534 / 0.2e1 + t505 * t431 + t507 * t428 + (t136 * t593 + t135 * t599 + (t205 * t597 + t206 * t593) * mrSges(6,3)) * pkin(10) + t443 - t453;
t10 = t101 * t613 + t238 * t509 + t239 * t508 + t650 * t372 + t658 * t624 + t437 - t447;
t8 = t433 - t435;
t6 = t441 * t429 + t439 - t450 + t654;
t2 = t434 + (t101 * t587 + (t429 * t66 + t432 * t60) * t640 + t208 * t589 + t207 * t595 - t159 * t594 + (t260 * t595 + t261 * t589) * mrSges(5,3)) * pkin(9) + t432 * t513 + t182 * t582 / 0.2e1 + t51 * t508 + t429 * t511 + t530 * t656 + t50 * t509 - t183 * t583 / 0.2e1 + (-t569 / 0.4e1 + t445) * t429 - t440 + t442 * t432 + t110 * t596 + t109 * t590;
t17 = [qJD(2) * t14 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t7, t538 + 0.2e1 * ((t238 * t310 + t239 * t311 + t525) * t640 + (t372 * t281 + t525) * t642 + (t500 * t643 + m(4) * (t337 * t430 + t466 - t500) / 0.2e1) * t427) * qJD(2) + t8 * qJD(3) + t10 * qJD(4) + t16 * qJD(5), t8 * qJD(2) + t2 * qJD(4) + t6 * qJD(5) + t555 + (-t194 * mrSges(4,1) - t193 * mrSges(4,2) + pkin(3) * t230 + t347 * t182 + t346 * t183 + t194 * t395 + t227 * t616 + t228 * t615 + t82 * t375 + t51 * t386 + t50 * t388 + t517 + (pkin(9) * t236 + t119 * mrSges(5,3) + t568 / 0.2e1 - t493) * t432 + (t110 * t591 + t109 * t599 - t118 * mrSges(5,3) + t192 / 0.2e1 + t575 / 0.2e1 + (-t237 + t140) * pkin(9)) * t429 + 0.2e1 * (-pkin(3) * t194 + (-t118 * t429 + t119 * t432) * pkin(9)) * t642 + 0.2e1 * (t346 * t50 + t347 * t51 + t82 * t585) * t640 + (t494 + (t402 / 0.2e1 - t355 / 0.2e1) * t429) * t307) * qJD(3), t549 + t10 * qJD(2) + t2 * qJD(3) + (-t166 * t583 + pkin(4) * t159 + t122 * t597 + t121 * t591 + t165 * t582 + t261 * t601 - t65 * mrSges(5,2) + (t400 * t597 + t403 * t593) * t260 + t518 + t646 * (-t428 * t54 + t431 * t55) + t644 * t66) * qJD(4) + t12 * qJD(5), t543 + t16 * qJD(2) + t6 * qJD(3) + t12 * qJD(4) + (-mrSges(6,1) * t37 - mrSges(6,2) * t36 + t102) * qJD(5); qJD(3) * t9 + qJD(4) * t11 + qJD(5) * t15 - t538, t99 * qJD(3) + t89 * qJD(4), (-mrSges(4,2) * t503 + m(5) * (t393 + (t586 * pkin(9) * t425 - pkin(3) * t430) * t427) + t375 * t482 + m(6) * (t344 * t346 + t345 * t347 + t393) + t345 * t386 + t344 * t388 + (-mrSges(4,1) + t395) * t524 + (t425 * t503 + t483) * mrSges(5,3)) * qJD(3) + t27 * qJD(4) + t41 * qJD(5) + t470, t27 * qJD(3) + (t644 * t372 + (mrSges(5,2) + t646 * (-t422 - t424)) * t371) * qJD(4) + t87 * qJD(5) + t469, t537 + t41 * qJD(3) + t87 * qJD(4) + (-mrSges(6,1) * t311 - mrSges(6,2) * t310) * qJD(5); -qJD(2) * t9 + qJD(4) * t1 + qJD(5) * t5 - t555, -qJD(4) * t26 - qJD(5) * t40 - t470, -qJD(4) * t39 - qJD(5) * t56, t43 * qJD(5) + t455 + (t400 * t497 - t475 * t584 + t403 * t495 - pkin(9) * t544 + m(6) * (-pkin(4) * t584 + t649 * pkin(10)) + t359 * t597 + t357 * t591 + t387 * t582 - t389 * t583 - pkin(4) * t376 + t419 + (pkin(9) * mrSges(5,2) - Ifges(5,6) + t601) * t429 + t649 * mrSges(6,3)) * qJD(4), t43 * qJD(4) + (-mrSges(6,1) * t347 - mrSges(6,2) * t346 - t429 * t398) * qJD(5) + t454; -qJD(2) * t11 - qJD(3) * t1 + qJD(5) * t13 - t549, qJD(3) * t26 - qJD(5) * t86 - t469, qJD(5) * t42 - t455, t248 * qJD(5), (-pkin(10) * t475 + t652) * qJD(5) + t451; -qJD(2) * t15 - qJD(3) * t5 - qJD(4) * t13 - t543, t40 * qJD(3) + t86 * qJD(4) - t537, -qJD(4) * t42 - t454, -t451, 0;];
Cq = t17;
