% Calculate time derivative of joint inertia matrix for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:53
% EndTime: 2019-03-09 09:59:39
% DurationCPUTime: 33.73s
% Computational Cost: add. (28227->1176), mult. (50012->1597), div. (0->0), fcn. (47029->8), ass. (0->541)
t418 = cos(qJ(2));
t415 = sin(qJ(2));
t647 = Icges(4,6) * t415;
t659 = Icges(3,4) * t415;
t734 = -t647 - t659 + (-Icges(3,2) - Icges(4,3)) * t418;
t646 = Icges(4,6) * t418;
t658 = Icges(3,4) * t418;
t733 = -t646 - t658 + (-Icges(3,1) - Icges(4,2)) * t415;
t417 = cos(qJ(4));
t414 = sin(qJ(4));
t656 = Icges(5,4) * t414;
t501 = Icges(5,2) * t417 + t656;
t306 = Icges(5,6) * t415 - t418 * t501;
t655 = Icges(5,4) * t417;
t507 = Icges(5,1) * t414 + t655;
t309 = Icges(5,5) * t415 - t418 * t507;
t408 = qJ(4) + pkin(9);
t399 = sin(t408);
t400 = cos(t408);
t649 = Icges(7,5) * t400;
t505 = Icges(7,1) * t399 - t649;
t288 = Icges(7,4) * t415 - t418 * t505;
t653 = Icges(6,4) * t400;
t506 = Icges(6,1) * t399 + t653;
t289 = Icges(6,5) * t415 - t418 * t506;
t718 = -t289 - t288;
t732 = t306 * t417 + t309 * t414 - t718 * t399;
t731 = t734 * qJD(2);
t730 = t733 * qJD(2);
t496 = Icges(6,5) * t399 + Icges(6,6) * t400;
t285 = Icges(6,3) * t415 - t418 * t496;
t499 = Icges(7,4) * t399 - Icges(7,6) * t400;
t286 = Icges(7,2) * t415 - t418 * t499;
t497 = Icges(5,5) * t414 + Icges(5,6) * t417;
t303 = Icges(5,3) * t415 - t418 * t497;
t729 = t286 + t285 + t303;
t416 = sin(qJ(1));
t591 = qJD(4) * t415;
t542 = qJD(1) + t591;
t470 = t416 * t542;
t419 = cos(qJ(1));
t636 = t400 * t419;
t594 = qJD(2) * t418;
t555 = t416 * t594;
t597 = qJD(1) * t419;
t699 = t415 * t597 + t555;
t180 = -qJD(4) * t636 + t399 * t470 - t699 * t400;
t541 = qJD(1) * t415 + qJD(4);
t467 = t541 * t419;
t430 = t467 + t555;
t181 = t399 * t430 + t400 * t470;
t632 = t416 * t400;
t301 = t399 * t419 + t415 * t632;
t661 = rSges(7,3) + qJ(6);
t720 = rSges(7,1) + pkin(5);
t728 = t301 * qJD(6) - t661 * t180 - t720 * t181;
t634 = t415 * t416;
t302 = t399 * t634 - t636;
t727 = -t301 * t661 + t720 * t302;
t504 = -Icges(3,2) * t415 + t658;
t307 = -Icges(3,6) * t419 + t416 * t504;
t492 = -Icges(4,3) * t415 + t646;
t694 = Icges(4,5) * t419 + t416 * t492;
t725 = -t307 - t694;
t308 = Icges(3,6) * t416 + t419 * t504;
t312 = Icges(4,5) * t416 - t419 * t492;
t724 = t308 - t312;
t509 = Icges(3,1) * t418 - t659;
t310 = -Icges(3,5) * t419 + t416 * t509;
t494 = Icges(4,2) * t418 - t647;
t693 = Icges(4,4) * t419 + t416 * t494;
t723 = t310 + t693;
t311 = Icges(3,5) * t416 + t419 * t509;
t314 = Icges(4,4) * t416 - t419 * t494;
t722 = t311 - t314;
t589 = qJD(4) * t418;
t650 = Icges(7,5) * t399;
t209 = (-Icges(7,1) * t400 - t650) * t589 + (Icges(7,4) * t418 + t415 * t505) * qJD(2);
t654 = Icges(6,4) * t399;
t210 = (-Icges(6,1) * t400 + t654) * t589 + (Icges(6,5) * t418 + t415 * t506) * qJD(2);
t719 = -t210 - t209;
t593 = qJD(2) * t419;
t556 = t415 * t593;
t598 = qJD(1) * t418;
t717 = t416 * t598 + t556;
t630 = t416 * t417;
t334 = t414 * t419 + t415 * t630;
t628 = t417 * t419;
t631 = t416 * t414;
t335 = t415 * t631 - t628;
t629 = t416 * t418;
t234 = Icges(5,5) * t335 + Icges(5,6) * t334 + Icges(5,3) * t629;
t236 = Icges(5,4) * t335 + Icges(5,2) * t334 + Icges(5,6) * t629;
t238 = Icges(5,1) * t335 + Icges(5,4) * t334 + Icges(5,5) * t629;
t483 = t236 * t417 + t238 * t414;
t103 = t234 * t415 - t418 * t483;
t194 = Icges(7,4) * t302 + Icges(7,2) * t629 - Icges(7,6) * t301;
t190 = Icges(7,5) * t302 + Icges(7,6) * t629 - Icges(7,3) * t301;
t198 = Icges(7,1) * t302 + Icges(7,4) * t629 - Icges(7,5) * t301;
t487 = t190 * t400 - t198 * t399;
t90 = t194 * t415 + t418 * t487;
t192 = Icges(6,5) * t302 + Icges(6,6) * t301 + Icges(6,3) * t629;
t196 = Icges(6,4) * t302 + Icges(6,2) * t301 + Icges(6,6) * t629;
t200 = Icges(6,1) * t302 + Icges(6,4) * t301 + Icges(6,5) * t629;
t485 = t196 * t400 + t200 * t399;
t92 = t192 * t415 - t418 * t485;
t579 = t103 + t90 + t92;
t332 = t415 * t628 - t631;
t633 = t415 * t419;
t574 = t414 * t633;
t333 = t574 + t630;
t627 = t418 * t419;
t233 = Icges(5,5) * t333 + Icges(5,6) * t332 + Icges(5,3) * t627;
t235 = Icges(5,4) * t333 + Icges(5,2) * t332 + Icges(5,6) * t627;
t237 = Icges(5,1) * t333 + Icges(5,4) * t332 + Icges(5,5) * t627;
t484 = t235 * t417 + t237 * t414;
t102 = t233 * t415 - t418 * t484;
t299 = t399 * t416 - t400 * t633;
t300 = t399 * t633 + t632;
t193 = Icges(7,4) * t300 + Icges(7,2) * t627 + Icges(7,6) * t299;
t189 = Icges(7,5) * t300 + Icges(7,6) * t627 + Icges(7,3) * t299;
t197 = Icges(7,1) * t300 + Icges(7,4) * t627 + Icges(7,5) * t299;
t488 = t189 * t400 - t197 * t399;
t89 = t193 * t415 + t418 * t488;
t191 = Icges(6,5) * t300 - Icges(6,6) * t299 + Icges(6,3) * t627;
t195 = Icges(6,4) * t300 - Icges(6,2) * t299 + Icges(6,6) * t627;
t199 = Icges(6,1) * t300 - Icges(6,4) * t299 + Icges(6,5) * t627;
t486 = t195 * t400 + t199 * t399;
t91 = t191 * t415 - t418 * t486;
t580 = t102 + t89 + t91;
t715 = t416 * t579 + t419 * t580;
t714 = -t416 * t580 + t419 * t579;
t495 = -Icges(7,3) * t400 + t650;
t284 = Icges(7,6) * t415 - t418 * t495;
t500 = Icges(6,2) * t400 + t654;
t287 = Icges(6,6) * t415 - t418 * t500;
t713 = ((t284 - t287) * t400 - t732) * t418 + t729 * t415;
t205 = (-Icges(7,3) * t399 - t649) * t589 + (Icges(7,6) * t418 + t415 * t495) * qJD(2);
t206 = (-Icges(6,5) * t400 + Icges(6,6) * t399) * t589 + (Icges(6,3) * t418 + t415 * t496) * qJD(2);
t207 = (-Icges(7,4) * t400 - Icges(7,6) * t399) * t589 + (Icges(7,2) * t418 + t415 * t499) * qJD(2);
t243 = (-Icges(5,5) * t417 + Icges(5,6) * t414) * t589 + (Icges(5,3) * t418 + t415 * t497) * qJD(2);
t553 = t399 * t589;
t596 = qJD(2) * t415;
t558 = t400 * t596;
t637 = t400 * t418;
t712 = t414 * t306 * t589 + t205 * t637 - t284 * t558 + (t553 + t558) * t287 + t729 * t594 + (t207 + t243 + t206) * t415 + t732 * t596;
t602 = t416 ^ 2 + t419 ^ 2;
t711 = -0.1e1 + t602;
t473 = t312 * t415 - t314 * t418;
t710 = t416 * t473;
t474 = t308 * t415 - t311 * t418;
t709 = t416 * t474;
t472 = -t415 * t694 + t418 * t693;
t706 = t419 * t472;
t475 = t307 * t415 - t310 * t418;
t705 = t419 * t475;
t554 = t418 * t593;
t182 = qJD(1) * t301 + qJD(4) * t300 - t400 * t554;
t429 = -t416 * t541 + t554;
t469 = t419 * t542;
t183 = t399 * t429 + t400 * t469;
t702 = t299 * qJD(6) + t661 * t182 + t720 * t183;
t701 = t416 * rSges(4,1) - rSges(4,2) * t627;
t352 = t416 * pkin(3) + pkin(8) * t627;
t411 = t418 ^ 2;
t700 = -t415 ^ 2 + t411;
t698 = -rSges(3,2) * t633 + t416 * rSges(3,3);
t107 = Icges(7,5) * t183 - Icges(7,6) * t717 + Icges(7,3) * t182;
t111 = Icges(7,4) * t183 - Icges(7,2) * t717 + Icges(7,6) * t182;
t115 = Icges(7,1) * t183 - Icges(7,4) * t717 + Icges(7,5) * t182;
t27 = (-qJD(2) * t488 + t111) * t415 + (qJD(2) * t193 + t107 * t400 - t115 * t399 + (-t189 * t399 - t197 * t400) * qJD(4)) * t418;
t109 = Icges(6,5) * t183 - Icges(6,6) * t182 - Icges(6,3) * t717;
t113 = Icges(6,4) * t183 - Icges(6,2) * t182 - Icges(6,6) * t717;
t117 = Icges(6,1) * t183 - Icges(6,4) * t182 - Icges(6,5) * t717;
t29 = (qJD(2) * t486 + t109) * t415 + (qJD(2) * t191 - t113 * t400 - t117 * t399 + (t195 * t399 - t199 * t400) * qJD(4)) * t418;
t468 = t542 * t414;
t218 = t417 * t429 - t419 * t468;
t219 = t414 * t429 + t417 * t469;
t129 = Icges(5,5) * t219 + Icges(5,6) * t218 - Icges(5,3) * t717;
t131 = Icges(5,4) * t219 + Icges(5,2) * t218 - Icges(5,6) * t717;
t133 = Icges(5,1) * t219 + Icges(5,4) * t218 - Icges(5,5) * t717;
t35 = (qJD(2) * t484 + t129) * t415 + (qJD(2) * t233 - t131 * t417 - t133 * t414 + (t235 * t414 - t237 * t417) * qJD(4)) * t418;
t697 = t27 + t29 + t35;
t595 = qJD(2) * t416;
t557 = t415 * t595;
t560 = t418 * t597;
t448 = -t557 + t560;
t106 = Icges(7,5) * t181 + Icges(7,6) * t448 + Icges(7,3) * t180;
t110 = Icges(7,4) * t181 + Icges(7,2) * t448 + Icges(7,6) * t180;
t114 = Icges(7,1) * t181 + Icges(7,4) * t448 + Icges(7,5) * t180;
t28 = (-qJD(2) * t487 + t110) * t415 + (qJD(2) * t194 + t106 * t400 - t114 * t399 + (-t190 * t399 - t198 * t400) * qJD(4)) * t418;
t108 = Icges(6,5) * t181 - Icges(6,6) * t180 + Icges(6,3) * t448;
t112 = Icges(6,4) * t181 - Icges(6,2) * t180 + Icges(6,6) * t448;
t116 = Icges(6,1) * t181 - Icges(6,4) * t180 + Icges(6,5) * t448;
t30 = (qJD(2) * t485 + t108) * t415 + (qJD(2) * t192 - t112 * t400 - t116 * t399 + (t196 * t399 - t200 * t400) * qJD(4)) * t418;
t216 = t417 * t467 + (t417 * t594 - t468) * t416;
t217 = t414 * t430 + t417 * t470;
t128 = Icges(5,5) * t217 + Icges(5,6) * t216 + Icges(5,3) * t448;
t130 = Icges(5,4) * t217 + Icges(5,2) * t216 + Icges(5,6) * t448;
t132 = Icges(5,1) * t217 + Icges(5,4) * t216 + Icges(5,5) * t448;
t36 = (qJD(2) * t483 + t128) * t415 + (qJD(2) * t234 - t130 * t417 - t132 * t414 + (t236 * t414 - t238 * t417) * qJD(4)) * t418;
t696 = t28 + t30 + t36;
t122 = t299 * t284 + t286 * t627 + t300 * t288;
t123 = t285 * t627 - t299 * t287 + t300 * t289;
t146 = t303 * t627 + t332 * t306 + t333 * t309;
t93 = t233 * t627 + t332 * t235 + t333 * t237;
t94 = t234 * t627 + t332 * t236 + t333 * t238;
t513 = t416 * t94 + t419 * t93;
t79 = t191 * t627 - t299 * t195 + t300 * t199;
t80 = t192 * t627 - t299 * t196 + t300 * t200;
t516 = t416 * t80 + t419 * t79;
t77 = t299 * t189 + t193 * t627 + t300 * t197;
t78 = t299 * t190 + t194 * t627 + t300 * t198;
t517 = t416 * t78 + t419 * t77;
t695 = (t513 + t516 + t517) * t418 + (t122 + t123 + t146) * t415;
t124 = -t284 * t301 + t286 * t629 + t288 * t302;
t125 = t285 * t629 + t287 * t301 + t289 * t302;
t147 = t303 * t629 + t306 * t334 + t309 * t335;
t95 = t233 * t629 + t235 * t334 + t237 * t335;
t96 = t234 * t629 + t236 * t334 + t238 * t335;
t512 = t416 * t96 + t419 * t95;
t83 = t191 * t629 + t195 * t301 + t199 * t302;
t84 = t192 * t629 + t196 * t301 + t200 * t302;
t514 = t416 * t84 + t419 * t83;
t81 = -t189 * t301 + t193 * t629 + t197 * t302;
t82 = -t190 * t301 + t194 * t629 + t198 * t302;
t515 = t416 * t82 + t419 * t81;
t583 = (t512 + t514 + t515) * t418 + (t124 + t125 + t147) * t415;
t498 = Icges(3,5) * t418 - Icges(3,6) * t415;
t304 = -Icges(3,3) * t419 + t416 * t498;
t502 = Icges(4,4) * t418 - Icges(4,5) * t415;
t692 = Icges(4,1) * t419 + t416 * t502;
t670 = pkin(4) * t414;
t691 = -pkin(8) * t598 + t541 * t670;
t690 = 2 * m(3);
t689 = 2 * m(4);
t688 = 2 * m(5);
t687 = 2 * m(6);
t686 = 2 * m(7);
t685 = m(4) / 0.2e1;
t684 = m(5) / 0.2e1;
t683 = -m(6) / 0.2e1;
t682 = m(6) / 0.2e1;
t681 = -m(7) / 0.2e1;
t680 = m(7) / 0.2e1;
t678 = t416 / 0.2e1;
t677 = -t419 / 0.2e1;
t674 = -rSges(7,2) - pkin(2);
t673 = -rSges(6,3) - pkin(2);
t370 = rSges(3,1) * t415 + rSges(3,2) * t418;
t672 = m(3) * t370;
t671 = pkin(2) * t418;
t413 = -qJ(5) - pkin(8);
t667 = -pkin(8) - t413;
t666 = rSges(4,1) * t419;
t665 = rSges(4,2) * t415;
t664 = rSges(3,3) * t419;
t663 = rSges(5,3) * t415;
t662 = -rSges(4,3) - qJ(3);
t644 = qJ(3) * t415;
t643 = qJ(3) * t418;
t523 = -rSges(6,1) * t302 - rSges(6,2) * t301;
t204 = rSges(6,3) * t629 - t523;
t642 = t204 * t419;
t526 = -rSges(5,1) * t335 - rSges(5,2) * t334;
t240 = rSges(5,3) * t629 - t526;
t641 = t240 * t419;
t635 = t414 * t418;
t626 = rSges(7,2) * t448 - t728;
t625 = -rSges(7,2) * t717 + t702;
t360 = t413 * t560;
t380 = pkin(8) * t557;
t394 = pkin(4) * t417 + pkin(3);
t590 = qJD(4) * t417;
t551 = t415 * t590;
t588 = qJD(5) * t418;
t140 = -t360 + t380 + t691 * t419 + (t413 * t596 + t588 + (-pkin(3) + t394) * qJD(1) + (t414 * t594 + t551) * pkin(4)) * t416;
t462 = pkin(4) * t574 + t416 * t394 - t413 * t627;
t261 = t462 - t352;
t624 = t140 * t627 + t261 * t557;
t622 = t183 * rSges(6,1) - t182 * rSges(6,2);
t521 = rSges(7,1) * t399 - rSges(7,3) * t400;
t621 = (pkin(5) * t596 - qJ(6) * t589) * t399 + (-qJ(6) * t596 + (-pkin(5) * qJD(4) + qJD(6)) * t418) * t400 + (-rSges(7,1) * t400 - rSges(7,3) * t399) * t589 + (rSges(7,2) * t418 + t415 * t521) * qJD(2);
t620 = rSges(7,2) * t627 + t661 * t299 + t720 * t300;
t202 = t300 * rSges(6,1) - t299 * rSges(6,2) + rSges(6,3) * t627;
t619 = -t202 - t261;
t618 = rSges(7,2) * t629 + t727;
t406 = t419 * pkin(3);
t581 = t415 * t670;
t606 = -t419 * t394 - t413 * t629;
t262 = t406 + (-pkin(8) * t418 + t581) * t416 + t606;
t617 = -t204 - t262;
t616 = t219 * rSges(5,1) + t218 * rSges(5,2);
t550 = t417 * t589;
t263 = -pkin(4) * t550 + qJD(5) * t415 + (t418 * t667 + t581) * qJD(2);
t327 = -pkin(4) * t635 + t415 * t667;
t615 = t263 * t629 + t327 * t560;
t614 = rSges(7,2) * t415 + (-pkin(5) * t399 + qJ(6) * t400 - t521) * t418;
t522 = rSges(6,1) * t399 + rSges(6,2) * t400;
t292 = rSges(6,3) * t415 - t418 * t522;
t613 = -t292 - t327;
t518 = t644 + t671;
t338 = t518 * t416;
t339 = pkin(2) * t627 + qJ(3) * t633;
t612 = t416 * t338 + t419 * t339;
t328 = qJD(2) * t518 - qJD(3) * t418;
t520 = -rSges(4,2) * t418 + rSges(4,3) * t415;
t611 = -t520 * qJD(2) - t328;
t610 = -t339 - t352;
t368 = pkin(2) * t415 - t643;
t599 = qJD(1) * t416;
t340 = t368 * t599;
t565 = t415 * t599;
t609 = pkin(8) * t565 + t340;
t381 = pkin(2) * t557;
t608 = t360 + t381;
t519 = rSges(4,3) * t418 + t665;
t607 = -t368 + t519;
t592 = qJD(3) * t415;
t605 = qJ(3) * t554 + t419 * t592;
t604 = rSges(3,2) * t565 + rSges(3,3) * t597;
t603 = t419 * pkin(1) + t416 * pkin(7);
t305 = Icges(3,3) * t416 + t419 * t498;
t601 = qJD(1) * t305;
t316 = Icges(4,1) * t416 - t419 * t502;
t600 = qJD(1) * t316;
t586 = -rSges(5,3) - pkin(2) - pkin(8);
t585 = t682 + t680;
t576 = qJD(4) * t670;
t573 = -t261 - t620;
t572 = -t262 - t618;
t571 = t416 * (pkin(2) * t560 + t699 * qJ(3) + t416 * t592 - t381) + t419 * (-pkin(2) * t717 - qJ(3) * t565 + t605) + t338 * t597;
t570 = -t261 + t610;
t569 = -t327 - t614;
t568 = t327 * t599 + t609;
t239 = t333 * rSges(5,1) + t332 * rSges(5,2) + rSges(5,3) * t627;
t405 = t419 * pkin(7);
t567 = t405 - t606;
t397 = pkin(7) * t597;
t566 = t397 + t605;
t563 = t292 * t599;
t525 = rSges(5,1) * t414 + rSges(5,2) * t417;
t321 = -t418 * t525 + t663;
t562 = t321 * t599;
t549 = t411 * qJD(4) * t399;
t548 = t415 * t594;
t547 = (-t502 / 0.2e1 + t498 / 0.2e1) * qJD(2);
t546 = -qJ(3) - t670;
t545 = -pkin(8) * t415 - t368;
t544 = t618 * t419;
t277 = t607 * t419;
t543 = qJD(1) * t614;
t398 = pkin(3) * t597;
t459 = t394 * t597 + t554 * t670 + (pkin(4) * t551 + t588) * t419 + t717 * t413;
t141 = pkin(8) * t556 - t691 * t416 - t398 + t459;
t540 = t415 * t141 + t261 * t594 + t717 * t327;
t353 = pkin(8) * t629 - t406;
t539 = t419 * t352 + t416 * t353 + t612;
t538 = rSges(4,1) * t597 + t717 * rSges(4,2) + rSges(4,3) * t554;
t537 = t603 + t339;
t51 = t77 * t416 - t419 * t78;
t52 = t79 * t416 - t419 * t80;
t62 = t93 * t416 - t419 * t94;
t536 = -t51 / 0.2e1 - t52 / 0.2e1 - t62 / 0.2e1;
t53 = t81 * t416 - t419 * t82;
t54 = t83 * t416 - t419 * t84;
t63 = t95 * t416 - t419 * t96;
t535 = t53 / 0.2e1 + t54 / 0.2e1 + t63 / 0.2e1;
t534 = -t321 + t545;
t533 = -t327 + t545;
t532 = -pkin(8) * t594 - t328;
t531 = t416 * t543;
t530 = t415 * t662 - pkin(1);
t529 = t546 * t418;
t528 = rSges(3,1) * t418 - rSges(3,2) * t415;
t527 = t217 * rSges(5,1) + t216 * rSges(5,2);
t524 = t181 * rSges(6,1) - t180 * rSges(6,2);
t458 = t415 * t546 - pkin(1);
t428 = t418 * t674 + t458;
t100 = t416 * t428 + t567 - t727;
t431 = t462 + t537;
t101 = t431 + t620;
t490 = t100 * t419 + t101 * t416;
t427 = t418 * t673 + t458;
t142 = t416 * t427 + t523 + t567;
t143 = t431 + t202;
t489 = t142 * t419 + t143 * t416;
t482 = -t239 * t419 - t416 * t240;
t481 = t239 * t416 - t641;
t476 = t299 * t419 - t301 * t416;
t471 = -t292 + t533;
t322 = rSges(3,1) * t627 + t698;
t323 = rSges(4,3) * t633 + t701;
t259 = (-rSges(5,1) * t417 + rSges(5,2) * t414) * t589 + (rSges(5,3) * t418 + t415 * t525) * qJD(2);
t466 = -t259 + t532;
t465 = -t263 + t532;
t208 = (Icges(6,2) * t399 - t653) * t589 + (Icges(6,6) * t418 + t415 * t500) * qJD(2);
t246 = (Icges(5,2) * t414 - t655) * t589 + (Icges(5,6) * t418 + t415 * t501) * qJD(2);
t249 = (-Icges(5,1) * t417 + t656) * t589 + (Icges(5,5) * t418 + t415 * t507) * qJD(2);
t464 = t713 * t594 + ((-t414 * t249 + (-qJD(4) * t309 - t246) * t417 + (t718 * qJD(4) - t208) * t400 + (-t284 * qJD(4) + t719) * t399) * t418 + t712) * t415;
t463 = -pkin(1) - t528;
t229 = t534 * t419;
t461 = t416 * (t352 * qJD(1) - t380) + t419 * (-pkin(8) * t717 + t398) + t353 * t597 + t571;
t460 = t419 * t261 + t416 * t262 + t539;
t457 = t533 - t614;
t212 = (-rSges(6,1) * t400 + rSges(6,2) * t399) * t589 + (rSges(6,3) * t418 + t415 * t522) * qJD(2);
t456 = -t212 + t465;
t455 = qJD(2) * t370;
t452 = qJD(2) * (Icges(4,4) * t415 + Icges(4,5) * t418);
t451 = qJD(2) * (-Icges(3,5) * t415 - Icges(3,6) * t418);
t171 = t471 * t419;
t446 = t465 - t621;
t445 = t418 * t586 - pkin(1) - t644;
t154 = t457 * t419;
t13 = (t416 * t620 + t419 * t572) * t596 + (t626 * t419 + (-t141 - t625) * t416 + (t416 * t572 + t419 * t573) * qJD(1)) * t418 + t624;
t290 = t327 * t629;
t74 = t415 * t572 + t614 * t629 + t290;
t230 = t415 * t261;
t75 = t415 * t620 + t569 * t627 + t230;
t444 = t593 * t74 + t595 * t75 - t13;
t119 = rSges(6,3) * t448 + t524;
t121 = -rSges(6,3) * t717 + t622;
t15 = (t202 * t416 + t419 * t617) * t596 + (t119 * t419 + (-t121 - t141) * t416 + (t416 * t617 + t419 * t619) * qJD(1)) * t418 + t624;
t98 = t292 * t629 + t415 * t617 + t290;
t99 = t415 * t202 + t613 * t627 + t230;
t443 = t593 * t98 + t595 * t99 - t15;
t436 = t416 * t140 + t419 * t141 + t262 * t597 + t461;
t14 = t625 * t419 + t626 * t416 + (t544 + (t570 - t620) * t416) * qJD(1) + t436;
t153 = t457 * t416;
t442 = t153 * t595 + t154 * t593 - t14;
t16 = t416 * t119 + t121 * t419 + (t642 + (-t202 + t570) * t416) * qJD(1) + t436;
t170 = t471 * t416;
t441 = t170 * t595 + t171 * t593 - t16;
t440 = (rSges(4,2) - pkin(2)) * t418 + t530;
t49 = t182 * t284 + t183 * t288 + t299 * t205 + t207 * t627 + t300 * t209 - t286 * t717;
t50 = -t182 * t287 + t183 * t289 + t206 * t627 - t299 * t208 + t300 * t210 - t285 * t717;
t59 = t218 * t306 + t219 * t309 + t243 * t627 + t332 * t246 + t333 * t249 - t303 * t717;
t439 = t35 / 0.2e1 + t27 / 0.2e1 + t29 / 0.2e1 + t59 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1;
t47 = t180 * t284 + t181 * t288 - t301 * t205 + t207 * t629 + t302 * t209 + t286 * t448;
t48 = -t180 * t287 + t181 * t289 + t206 * t629 + t301 * t208 + t302 * t210 + t285 * t448;
t58 = t216 * t306 + t217 * t309 + t243 * t629 + t334 * t246 + t335 * t249 + t303 * t448;
t438 = t36 / 0.2e1 + t28 / 0.2e1 + t30 / 0.2e1 + t58 / 0.2e1 + t47 / 0.2e1 + t48 / 0.2e1;
t435 = t459 + t566;
t434 = -t122 / 0.2e1 - t123 / 0.2e1 - t146 / 0.2e1 - t102 / 0.2e1 - t89 / 0.2e1 - t91 / 0.2e1;
t433 = t124 / 0.2e1 + t125 / 0.2e1 + t147 / 0.2e1 + t103 / 0.2e1 + t90 / 0.2e1 + t92 / 0.2e1;
t432 = t445 * t416;
t426 = t180 * t416 + t182 * t419 + (-t299 * t416 - t301 * t419) * qJD(1);
t425 = (-pkin(7) - t394) * qJD(1) + (-pkin(4) * t590 - qJD(3)) * t415 - t588;
t424 = qJD(1) * t428 - t576;
t423 = qJD(1) * t427 - t576;
t55 = t416 * t424 + t556 * t674 + t435 + t702;
t56 = t424 * t419 + ((t529 + (rSges(7,2) - t413) * t415) * qJD(2) + t425) * t416 + t608 + t728;
t68 = t416 * t423 + t556 * t673 + t435 + t622;
t69 = t423 * t419 + ((t529 + (rSges(6,3) - t413) * t415) * qJD(2) + t425) * t416 - t524 + t608;
t422 = (-t100 * t599 + t101 * t597 + t416 * t55 + t419 * t56) * t680 + (-t142 * t599 + t143 * t597 + t416 * t68 + t419 * t69) * t682;
t17 = (t593 * t614 + t625) * t415 + (t620 * qJD(2) + t531 + (-t263 - t621) * t419) * t418 + t540;
t18 = (t569 * t595 - t140 - t626) * t415 + (qJD(2) * t572 + t416 * t621 + t419 * t543) * t418 + t615;
t37 = (t292 * t593 + t121) * t415 + (t563 + qJD(2) * t202 + (-t212 - t263) * t419) * t418 + t540;
t38 = (t416 * t212 + t292 * t597) * t418 + (-t119 - t140) * t415 + (t418 * t617 + t613 * t634) * qJD(2) + t615;
t220 = t262 * t627;
t64 = t220 + (t416 * t573 + t544) * t418;
t76 = t220 + (t416 * t619 + t642) * t418;
t421 = (qJD(2) * t64 + t17 * t416 + t18 * t419 + t597 * t75 - t599 * t74) * t680 + (qJD(2) * t76 + t37 * t416 + t38 * t419 + t597 * t99 - t599 * t98) * t682;
t57 = t416 * t618 + t419 * t620 + t460;
t66 = t419 * t446 + t531 + t568;
t67 = qJD(1) * t154 + t416 * t446;
t73 = t202 * t419 + t416 * t204 + t460;
t85 = t419 * t456 + t563 + t568;
t86 = qJD(1) * t171 + t416 * t456;
t420 = (qJD(2) * t57 + t153 * t597 - t154 * t599 + t416 * t67 + t419 * t66) * t680 + (qJD(2) * t73 + t170 * t597 - t171 * t599 + t416 * t86 + t419 * t85) * t682;
t348 = t528 * qJD(2);
t324 = t416 * t520 - t666;
t320 = t416 * t528 - t664;
t276 = t607 * t416;
t272 = t322 + t603;
t271 = t416 * t463 + t405 + t664;
t270 = t711 * t548;
t257 = t692 * qJD(1) + t419 * t452;
t256 = t416 * t452 + t600;
t245 = t416 * t451 + t601;
t244 = -qJD(1) * t304 + t419 * t451;
t242 = t323 + t537;
t241 = t416 * t440 + t405 + t666;
t228 = t534 * t416;
t188 = t370 * t595 + ((-rSges(3,3) - pkin(7)) * t416 + t463 * t419) * qJD(1);
t187 = -rSges(3,1) * t717 - rSges(3,2) * t554 - pkin(1) * t599 + t397 + t604;
t173 = qJD(1) * t277 + t416 * t611;
t172 = t419 * t611 - t519 * t599 + t340;
t169 = t415 * t239 - t321 * t627;
t168 = -t240 * t415 + t321 * t629;
t167 = t416 * t472 + t419 * t692;
t166 = -t316 * t419 + t710;
t165 = -t416 * t692 + t706;
t164 = t416 * t316 + t419 * t473;
t163 = t416 * t305 - t419 * t474;
t162 = t416 * t304 - t705;
t160 = -t305 * t419 - t709;
t159 = -t304 * t419 - t416 * t475;
t158 = t537 + t239 + t352;
t157 = t405 + t406 + t432 + t526;
t155 = t323 * t419 + t416 * t324 + t612;
t152 = t381 + (-t592 + (t418 * t662 - t665) * qJD(2)) * t416 + ((-rSges(4,1) - pkin(7)) * t416 + t440 * t419) * qJD(1);
t151 = -pkin(2) * t556 + (t530 - t671) * t599 + t538 + t566;
t148 = t481 * t418;
t135 = -rSges(5,3) * t717 + t616;
t134 = rSges(5,3) * t448 + t527;
t127 = qJD(1) * t229 + t416 * t466;
t126 = t419 * t466 + t562 + t609;
t97 = -t482 + t539;
t88 = t380 + t381 + (-t592 + (-t643 + t663) * qJD(2)) * t416 + ((-pkin(3) - pkin(7)) * t416 + t445 * t419) * qJD(1) - t527;
t87 = qJD(1) * t432 + t556 * t586 + t398 + t566 + t616;
t72 = (-t321 * t595 - t134) * t415 + (-qJD(2) * t240 + t416 * t259 + t321 * t597) * t418;
t71 = (t321 * t593 + t135) * t415 + (qJD(2) * t239 - t259 * t419 + t562) * t418;
t70 = (qJD(1) * t324 + t538) * t419 + (t519 * t595 + (-t323 - t339 + t701) * qJD(1)) * t416 + t571;
t46 = t481 * t596 + (qJD(1) * t482 + t134 * t419 - t135 * t416) * t418;
t43 = t416 * t134 + t135 * t419 + (t641 + (-t239 + t610) * t416) * qJD(1) + t461;
t34 = t128 * t627 + t332 * t130 + t333 * t132 + t218 * t236 + t219 * t238 - t234 * t717;
t33 = t129 * t627 + t332 * t131 + t333 * t133 + t218 * t235 + t219 * t237 - t233 * t717;
t32 = t128 * t629 + t334 * t130 + t335 * t132 + t216 * t236 + t217 * t238 + t234 * t448;
t31 = t129 * t629 + t334 * t131 + t335 * t133 + t216 * t235 + t217 * t237 + t233 * t448;
t26 = t108 * t627 - t299 * t112 + t300 * t116 - t182 * t196 + t183 * t200 - t192 * t717;
t25 = t109 * t627 - t299 * t113 + t300 * t117 - t182 * t195 + t183 * t199 - t191 * t717;
t24 = t299 * t106 + t110 * t627 + t300 * t114 + t182 * t190 + t183 * t198 - t194 * t717;
t23 = t299 * t107 + t111 * t627 + t300 * t115 + t182 * t189 + t183 * t197 - t193 * t717;
t22 = t108 * t629 + t301 * t112 + t302 * t116 - t180 * t196 + t181 * t200 + t192 * t448;
t21 = t109 * t629 + t301 * t113 + t302 * t117 - t180 * t195 + t181 * t199 + t191 * t448;
t20 = -t301 * t106 + t110 * t629 + t302 * t114 + t180 * t190 + t181 * t198 + t194 * t448;
t19 = -t301 * t107 + t111 * t629 + t302 * t115 + t180 * t189 + t181 * t197 + t193 * t448;
t12 = qJD(1) * t513 + t33 * t416 - t34 * t419;
t11 = qJD(1) * t512 + t31 * t416 - t32 * t419;
t10 = qJD(1) * t516 + t25 * t416 - t26 * t419;
t9 = qJD(1) * t517 + t23 * t416 - t24 * t419;
t8 = qJD(1) * t514 + t21 * t416 - t22 * t419;
t7 = qJD(1) * t515 + t19 * t416 - t20 * t419;
t6 = (-qJD(2) * t513 + t59) * t415 + (-qJD(1) * t62 + qJD(2) * t146 + t33 * t419 + t34 * t416) * t418;
t5 = (-qJD(2) * t512 + t58) * t415 + (-qJD(1) * t63 + qJD(2) * t147 + t31 * t419 + t32 * t416) * t418;
t4 = (-qJD(2) * t516 + t50) * t415 + (-qJD(1) * t52 + qJD(2) * t123 + t25 * t419 + t26 * t416) * t418;
t3 = (-qJD(2) * t517 + t49) * t415 + (-qJD(1) * t51 + qJD(2) * t122 + t23 * t419 + t24 * t416) * t418;
t2 = (-qJD(2) * t514 + t48) * t415 + (-qJD(1) * t54 + qJD(2) * t125 + t21 * t419 + t22 * t416) * t418;
t1 = (-qJD(2) * t515 + t47) * t415 + (-qJD(1) * t53 + qJD(2) * t124 + t19 * t419 + t20 * t416) * t418;
t39 = [-t284 * t553 - t309 * t550 + (t100 * t56 + t101 * t55) * t686 + (t142 * t69 + t143 * t68) * t687 + (t157 * t88 + t158 * t87) * t688 + (t151 * t242 + t152 * t241) * t689 + (t187 * t272 + t188 * t271) * t690 - t249 * t635 - t208 * t637 + t718 * t400 * t589 + (t509 + t494 + t734) * t596 + (t504 + t492 - t733) * t594 + (-t246 * t417 + t719 * t399) * t418 + t712; m(4) * (t151 * t276 + t152 * t277 + t172 * t241 + t173 * t242) + m(5) * (t126 * t157 + t127 * t158 + t228 * t87 + t229 * t88) + m(6) * (t142 * t85 + t143 * t86 + t170 * t68 + t171 * t69) + m(7) * (t100 * t66 + t101 * t67 + t153 * t55 + t154 * t56) + (m(3) * (-t188 * t370 - t271 * t348) + t547 * t419 - t438) * t419 + (m(3) * (-t187 * t370 - t272 * t348) + t547 * t416 + t439) * t416 + ((-t724 * qJD(2) + t730 * t419) * t678 + (t725 * qJD(2) + t730 * t416) * t677 + (t722 * t677 - t723 * t678) * qJD(1)) * t415 + ((t722 * qJD(2) + t731 * t419) * t678 + (t723 * qJD(2) + t731 * t416) * t677 + (t724 * t677 + t725 * t678) * qJD(1)) * t418 + ((-t272 * t672 + (t308 / 0.2e1 - t312 / 0.2e1) * t418 + (t311 / 0.2e1 - t314 / 0.2e1) * t415 - t434) * t419 + (t271 * t672 + (t307 / 0.2e1 + t694 / 0.2e1) * t418 + (t310 / 0.2e1 + t693 / 0.2e1) * t415 + t433) * t416) * qJD(1); -t419 * ((t419 * t256 + (t166 - t706) * qJD(1)) * t419 + (t167 * qJD(1) + (t312 * t594 + t314 * t596 + t600) * t416 + (-t257 + (t415 * t693 + t418 * t694) * qJD(2) + t473 * qJD(1)) * t419) * t416) - t419 * ((t419 * t245 + (t160 + t705) * qJD(1)) * t419 + (t159 * qJD(1) + (-t308 * t594 - t311 * t596 + t601) * t416 + (-t244 + (t307 * t418 + t310 * t415) * qJD(2) - t474 * qJD(1)) * t419) * t416) - t419 * t11 - t419 * t8 - t419 * t7 + t416 * t12 + t416 * t10 + t416 * t9 + t416 * ((t416 * t257 + (t165 - t710) * qJD(1)) * t416 + (t164 * qJD(1) + (t594 * t694 + t596 * t693) * t419 + (-t256 + (t312 * t418 + t314 * t415) * qJD(2) + (t316 + t472) * qJD(1)) * t416) * t419) + t416 * ((t416 * t244 + (t162 + t709) * qJD(1)) * t416 + (t163 * qJD(1) + (t307 * t594 + t310 * t596) * t419 + (-t245 + (-t308 * t418 - t311 * t415) * qJD(2) + (t305 - t475) * qJD(1)) * t416) * t419) + (t14 * t57 + t153 * t67 + t154 * t66) * t686 + (t16 * t73 + t170 * t86 + t171 * t85) * t687 + (t126 * t229 + t127 * t228 + t43 * t97) * t688 + (t155 * t70 + t172 * t277 + t173 * t276) * t689 + ((t416 * t320 + t322 * t419) * ((qJD(1) * t320 - t419 * t455 + t604) * t419 + (-t416 * t455 + (-t322 + t698) * qJD(1)) * t416) + t602 * t370 * t348) * t690 + (t53 + t54 + t63 + (-t159 - t167) * t419 + (t160 + t166) * t416) * t599 + (t51 + t52 + t62 + (-t162 - t165) * t419 + (t163 + t164) * t416) * t597; 0.2e1 * (t490 * t680 + t489 * t682 + (t157 * t419 + t158 * t416) * t684 + (t241 * t419 + t242 * t416) * t685) * t594 + 0.2e1 * ((-t157 * t599 + t158 * t597 + t416 * t87 + t419 * t88) * t684 + (t151 * t416 + t152 * t419 - t241 * t599 + t242 * t597) * t685 + t422) * t415; 0.2e1 * (t442 * t680 + t441 * t682 + (t228 * t595 + t229 * t593 - t43) * t684 + (t276 * t595 + t277 * t593 - t70) * t685) * t418 + 0.2e1 * ((qJD(2) * t97 + t126 * t419 + t127 * t416 + t228 * t597 - t229 * t599) * t684 + (qJD(2) * t155 + t172 * t419 + t173 * t416 + t276 * t597 - t277 * t599) * t685 + t420) * t415; 0.4e1 * (t685 + t684 + t585) * t270; m(5) * (t157 * t72 + t158 * t71 + t168 * t88 + t169 * t87) + m(6) * (t142 * t38 + t143 * t37 + t68 * t99 + t69 * t98) + m(7) * (t100 * t18 + t101 * t17 + t55 * t75 + t56 * t74) + (-t416 * t433 + t419 * t434) * t596 + (t439 * t419 + t438 * t416 + (t416 * t434 + t419 * t433) * qJD(1)) * t418 + t464; m(5) * (t126 * t168 + t127 * t169 - t148 * t43 + t228 * t71 + t229 * t72 + t46 * t97) + m(6) * (t15 * t73 + t16 * t76 + t170 * t37 + t171 * t38 + t85 * t98 + t86 * t99) + m(7) * (t13 * t57 + t14 * t64 + t153 * t17 + t154 * t18 + t66 * t74 + t67 * t75) + (-t5 / 0.2e1 - t2 / 0.2e1 - t1 / 0.2e1 + t536 * t596) * t419 + (t6 / 0.2e1 + t4 / 0.2e1 + t3 / 0.2e1 - t535 * t596) * t416 + (t416 * t536 + t419 * t535) * t598 + (t7 + t8 + t11) * t418 * t678 + (t9 + t10 + t12) * t627 / 0.2e1 - t714 * t594 / 0.2e1 + (t715 * qJD(1) + t697 * t416 - t696 * t419) * t415 / 0.2e1 + (t583 * t416 + t419 * t695) * qJD(1) / 0.2e1; 0.2e1 * ((t168 * t593 + t169 * t595 - t46) * t684 + t443 * t682 + t444 * t680) * t418 + 0.2e1 * ((-qJD(2) * t148 - t168 * t599 + t169 * t597 + t416 * t71 + t419 * t72) * t684 + t421) * t415; (-t148 * t46 + t168 * t72 + t169 * t71) * t688 + (t15 * t76 + t37 * t99 + t38 * t98) * t687 + (t13 * t64 + t17 * t75 + t18 * t74) * t686 + (((-t415 * t580 - t695) * t419 + (-t415 * t579 - t583) * t416) * qJD(2) + t464) * t415 + ((t3 + t4 + t6) * t419 + (t1 + t2 + t5) * t416 + (t416 * t696 + t419 * t697) * t415 + (t713 * t415 + t715 * t418) * qJD(2) + (t714 * t415 - t416 * t695 + t583 * t419) * qJD(1)) * t418; 0.2e1 * (t489 * t683 + t490 * t681) * t596 + 0.2e1 * t422 * t418; 0.2e1 * (t441 * t683 + t442 * t681) * t415 + 0.2e1 * t420 * t418; 0.2e1 * t585 * t700 * t711 * qJD(2); 0.2e1 * (t443 * t683 + t444 * t681) * t415 + 0.2e1 * t421 * t418; -0.4e1 * t585 * t270; m(7) * (t100 * t182 + t101 * t180 + t299 * t56 - t301 * t55); m(7) * (-t57 * t553 + t153 * t180 + t154 * t182 + t299 * t66 - t301 * t67 + (t418 * t14 - t57 * t596) * t400); m(7) * (t549 + t476 * t594 + (0.2e1 * t400 * t594 + t426) * t415); m(7) * (-t64 * t553 - t17 * t301 + t18 * t299 + t180 * t75 + t182 * t74 + (t13 * t418 - t596 * t64) * t400); m(7) * ((t400 * t700 - t476 * t415) * qJD(2) + (-t399 * t591 + t426) * t418); (-t180 * t301 + t182 * t299 + (-t400 * t548 - t549) * t400) * t686;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t39(1) t39(2) t39(4) t39(7) t39(11) t39(16); t39(2) t39(3) t39(5) t39(8) t39(12) t39(17); t39(4) t39(5) t39(6) t39(9) t39(13) t39(18); t39(7) t39(8) t39(9) t39(10) t39(14) t39(19); t39(11) t39(12) t39(13) t39(14) t39(15) t39(20); t39(16) t39(17) t39(18) t39(19) t39(20) t39(21);];
Mq  = res;
