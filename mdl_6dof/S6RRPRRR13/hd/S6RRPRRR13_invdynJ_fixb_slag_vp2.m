% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:59
% EndTime: 2019-03-09 14:46:19
% DurationCPUTime: 48.12s
% Computational Cost: add. (21140->1080), mult. (49389->1423), div. (0->0), fcn. (37906->14), ass. (0->475)
t384 = cos(qJ(2));
t375 = cos(pkin(6));
t503 = qJD(1) * t375;
t488 = pkin(1) * t503;
t350 = t384 * t488;
t683 = qJD(3) - t350;
t378 = sin(qJ(4));
t383 = cos(qJ(4));
t436 = pkin(4) * t383 + pkin(10) * t378;
t379 = sin(qJ(2));
t374 = sin(pkin(6));
t504 = qJD(1) * t374;
t473 = t379 * t504;
t591 = pkin(3) + pkin(8);
t682 = -(-t436 - t591) * t473 + qJD(4) * t436 + t683;
t377 = sin(qJ(5));
t382 = cos(qJ(5));
t518 = t378 * t379;
t256 = (-t377 * t518 + t382 * t384) * t504;
t498 = qJD(4) * t378;
t469 = t377 * t498;
t493 = qJD(5) * t383;
t633 = t382 * t493 + t256 - t469;
t381 = cos(qJ(6));
t356 = qJD(2) + t503;
t472 = t384 * t504;
t263 = -t378 * t356 - t383 * t472;
t260 = qJD(5) - t263;
t444 = t378 * t472;
t264 = t356 * t383 - t444;
t328 = qJD(4) + t473;
t194 = t264 * t382 + t328 * t377;
t348 = pkin(8) * t472;
t284 = t379 * t488 + t348;
t246 = pkin(3) * t472 + t284;
t337 = t356 * qJ(3);
t203 = t337 + t246;
t111 = -pkin(4) * t263 - pkin(10) * t264 + t203;
t592 = pkin(2) + pkin(9);
t190 = -t356 * t592 + t473 * t591 + t683;
t458 = -qJ(3) * t379 - pkin(1);
t211 = (-t384 * t592 + t458) * t504;
t105 = t190 * t378 + t211 * t383;
t96 = pkin(10) * t328 + t105;
t60 = t382 * t111 - t377 * t96;
t45 = -pkin(11) * t194 + t60;
t42 = pkin(5) * t260 + t45;
t376 = sin(qJ(6));
t193 = -t264 * t377 + t328 * t382;
t61 = t111 * t377 + t382 * t96;
t46 = pkin(11) * t193 + t61;
t534 = t376 * t46;
t15 = t381 * t42 - t534;
t531 = t381 * t46;
t16 = t376 * t42 + t531;
t645 = t328 * Ifges(5,6);
t661 = t60 * mrSges(6,1) - t61 * mrSges(6,2) - t105 * mrSges(5,3);
t681 = t661 - t645 / 0.2e1 + t203 * mrSges(5,1) + t15 * mrSges(7,1) - t16 * mrSges(7,2);
t345 = pkin(2) * t473;
t417 = pkin(9) * t379 - qJ(3) * t384;
t244 = t417 * t504 + t345;
t149 = t383 * t244 + t378 * t246;
t131 = pkin(10) * t472 + t149;
t455 = -pkin(10) * t383 + qJ(3);
t320 = pkin(4) * t378 + t455;
t496 = qJD(4) * t592;
t467 = t383 * t496;
t517 = t378 * t592;
t478 = t377 * t517;
t494 = qJD(5) * t382;
t642 = qJD(5) * t478 + t320 * t494 + (-t131 - t467) * t382 + t682 * t377;
t680 = t131 * t377 + t382 * t682;
t367 = pkin(5) * t382 + pkin(4);
t373 = qJ(5) + qJ(6);
t370 = sin(t373);
t371 = cos(t373);
t430 = -mrSges(6,1) * t382 + mrSges(6,2) * t377;
t629 = m(6) * pkin(4) + m(7) * t367 + mrSges(7,1) * t371 - mrSges(7,2) * t370 - t430;
t676 = -m(5) - m(7);
t492 = m(6) - t676;
t673 = -mrSges(3,1) + mrSges(4,2);
t618 = t492 * pkin(9) + mrSges(5,3) - t673;
t428 = t370 * mrSges(7,1) + t371 * mrSges(7,2);
t530 = t382 * mrSges(6,2);
t429 = t377 * mrSges(6,1) + t530;
t559 = pkin(5) * t377;
t490 = m(7) * t559;
t679 = -t490 - t428 - t429;
t104 = t190 * t383 - t378 * t211;
t542 = t104 * mrSges(5,3);
t646 = t328 * Ifges(5,5);
t677 = t203 * mrSges(5,2) - t542 + t646 / 0.2e1;
t523 = t374 * t379;
t471 = qJD(2) * t523;
t521 = t374 * t384;
t291 = qJD(1) * t471 - qJDD(1) * t521;
t491 = qJDD(1) * t375;
t354 = qJDD(2) + t491;
t497 = qJD(4) * t383;
t166 = qJD(4) * t444 + t291 * t383 - t378 * t354 - t356 * t497;
t159 = qJDD(5) - t166;
t165 = t263 * qJD(4) + t291 * t378 + t354 * t383;
t500 = qJD(2) * t384;
t292 = (qJD(1) * t500 + qJDD(1) * t379) * t374;
t273 = qJDD(4) + t292;
t85 = qJD(5) * t193 + t165 * t382 + t273 * t377;
t86 = -qJD(5) * t194 - t165 * t377 + t273 * t382;
t35 = Ifges(6,5) * t85 + Ifges(6,6) * t86 + Ifges(6,3) * t159;
t151 = qJDD(6) + t159;
t454 = t381 * t193 - t194 * t376;
t31 = qJD(6) * t454 + t376 * t86 + t381 * t85;
t110 = t193 * t376 + t194 * t381;
t32 = -qJD(6) * t110 - t376 * t85 + t381 * t86;
t8 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t151;
t675 = t35 / 0.2e1 + t8 / 0.2e1;
t258 = Ifges(5,4) * t263;
t671 = t263 * Ifges(5,2);
t515 = t379 * t382;
t257 = (t377 * t384 + t378 * t515) * t504;
t352 = t382 * t517;
t445 = t383 * t473;
t456 = t377 * t592 + pkin(5);
t670 = pkin(5) * t445 + pkin(11) * t257 + (t352 + (pkin(11) * t383 - t320) * t377) * qJD(5) + (pkin(11) * t378 * t382 + t383 * t456) * qJD(4) + t680;
t669 = -pkin(11) * t633 + t642;
t386 = -pkin(11) - pkin(10);
t474 = qJD(5) * t386;
t525 = t263 * t377;
t189 = pkin(4) * t264 - pkin(10) * t263;
t75 = t382 * t104 + t377 * t189;
t668 = pkin(11) * t525 + t377 * t474 - t75;
t524 = t263 * t382;
t74 = -t104 * t377 + t382 * t189;
t667 = -pkin(5) * t264 + pkin(11) * t524 + t382 * t474 - t74;
t495 = qJD(5) * t377;
t666 = t495 - t525;
t665 = -m(7) * t386 + mrSges(6,3) + mrSges(7,3);
t357 = pkin(8) * t523;
t563 = pkin(1) * t375;
t487 = qJD(2) * t563;
t448 = qJD(1) * t487;
t480 = pkin(1) * t491;
t196 = -qJD(2) * t348 - qJDD(1) * t357 - t379 * t448 + t384 * t480;
t396 = qJDD(3) - t196;
t118 = pkin(3) * t292 - t354 * t592 + t396;
t499 = qJD(3) * t379;
t393 = -qJ(3) * t292 + (-pkin(1) * qJDD(1) - qJD(1) * t499) * t374;
t124 = t291 * t592 + t393;
t43 = t378 * t118 + t383 * t124 + t190 * t497 - t211 * t498;
t39 = pkin(10) * t273 + t43;
t195 = -pkin(8) * t291 + t379 * t480 + t384 * t448;
t164 = -t354 * qJ(3) - t356 * qJD(3) - t195;
t120 = -pkin(3) * t291 - t164;
t57 = -pkin(4) * t166 - pkin(10) * t165 + t120;
t13 = -qJD(5) * t61 - t377 * t39 + t382 * t57;
t6 = pkin(5) * t159 - pkin(11) * t85 + t13;
t12 = t111 * t494 + t377 * t57 + t382 * t39 - t495 * t96;
t7 = pkin(11) * t86 + t12;
t2 = qJD(6) * t15 + t376 * t6 + t381 * t7;
t3 = -qJD(6) * t16 - t376 * t7 + t381 * t6;
t664 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t663 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t569 = t273 / 0.2e1;
t582 = t166 / 0.2e1;
t583 = t165 / 0.2e1;
t599 = Ifges(5,1) * t583 + Ifges(5,4) * t582 + Ifges(5,5) * t569;
t537 = t264 * Ifges(5,4);
t155 = t537 + t645 + t671;
t251 = qJD(6) + t260;
t575 = t251 / 0.2e1;
t587 = t110 / 0.2e1;
t589 = t454 / 0.2e1;
t647 = t194 * Ifges(6,5) + t110 * Ifges(7,5) + t193 * Ifges(6,6) + Ifges(7,6) * t454 + t260 * Ifges(6,3) + t251 * Ifges(7,3);
t660 = Ifges(7,5) * t587 + Ifges(7,6) * t589 + Ifges(7,3) * t575 - t155 / 0.2e1 + t647 / 0.2e1;
t414 = t376 * t377 - t381 * t382;
t623 = qJD(5) + qJD(6);
t390 = t623 * t414;
t148 = -t378 * t244 + t246 * t383;
t130 = -pkin(4) * t472 - t148;
t468 = t378 * t496;
t657 = -t130 - t468;
t254 = t356 * t382 + t377 * t445;
t255 = t356 * t377 - t382 * t445;
t314 = t376 * t382 + t377 * t381;
t288 = t314 * t383;
t656 = -qJD(4) * t288 - t254 * t381 + t255 * t376 + t378 * t390;
t218 = t623 * t314;
t290 = t414 * t383;
t655 = -qJD(4) * t290 - t218 * t378 - t254 * t376 - t255 * t381;
t482 = m(4) + t492;
t653 = pkin(2) * t482 + t618 - t679;
t432 = mrSges(5,1) * t378 + mrSges(5,2) * t383;
t652 = -t378 * t629 - t432;
t608 = t31 / 0.2e1;
t607 = t32 / 0.2e1;
t598 = t85 / 0.2e1;
t597 = t86 / 0.2e1;
t586 = t151 / 0.2e1;
t584 = t159 / 0.2e1;
t556 = mrSges(4,3) - mrSges(3,2);
t651 = -Ifges(4,4) + Ifges(3,5);
t650 = Ifges(4,5) - Ifges(3,6);
t310 = t382 * t320;
t512 = t382 * t383;
t210 = -pkin(11) * t512 + t378 * t456 + t310;
t262 = t377 * t320 - t352;
t519 = t377 * t383;
t223 = -pkin(11) * t519 + t262;
t128 = t210 * t381 - t223 * t376;
t649 = qJD(6) * t128 + t376 * t670 + t381 * t669;
t129 = t210 * t376 + t223 * t381;
t648 = -qJD(6) * t129 - t376 * t669 + t381 * t670;
t335 = t386 * t377;
t336 = t386 * t382;
t249 = t335 * t381 + t336 * t376;
t644 = qJD(6) * t249 + t376 * t667 + t381 * t668;
t250 = t335 * t376 - t336 * t381;
t643 = -qJD(6) * t250 - t376 * t668 + t381 * t667;
t641 = -qJD(5) * t262 + t377 * t467 + t680;
t640 = pkin(5) * t633 + t657;
t114 = mrSges(5,1) * t273 - mrSges(5,3) * t165;
t41 = -mrSges(6,1) * t86 + mrSges(6,2) * t85;
t639 = -t41 + t114;
t610 = m(7) * pkin(5);
t638 = -t610 - mrSges(6,1);
t637 = pkin(5) * t666 - t105;
t562 = pkin(1) * t384;
t475 = -pkin(2) - t562;
t214 = pkin(3) * t523 + t357 + (-pkin(9) + t475) * t375;
t506 = pkin(2) * t521 + qJ(3) * t523;
t564 = pkin(1) * t374;
t267 = -t506 - t564;
t360 = pkin(9) * t521;
t238 = t267 - t360;
t135 = t378 * t214 + t383 * t238;
t126 = pkin(10) * t523 + t135;
t366 = t379 * t563;
t307 = pkin(8) * t521 + t366;
t266 = -t375 * qJ(3) - t307;
t237 = pkin(3) * t521 - t266;
t299 = t375 * t378 + t383 * t521;
t479 = t378 * t521;
t300 = t375 * t383 - t479;
t142 = pkin(4) * t299 - pkin(10) * t300 + t237;
t73 = t382 * t126 + t377 * t142;
t145 = -t218 * t383 + t414 * t498;
t172 = t256 * t376 + t257 * t381;
t636 = t145 - t172;
t147 = t314 * t498 + t383 * t390;
t171 = t256 * t381 - t257 * t376;
t635 = t147 - t171;
t112 = -mrSges(6,1) * t193 + mrSges(6,2) * t194;
t205 = mrSges(5,1) * t328 - mrSges(5,3) * t264;
t634 = t205 - t112;
t632 = t377 * t493 + t382 * t498 + t257;
t188 = -mrSges(5,1) * t263 + mrSges(5,2) * t264;
t451 = mrSges(4,1) * t472;
t279 = -mrSges(4,3) * t356 - t451;
t631 = -t279 + t188;
t450 = mrSges(3,3) * t473;
t452 = mrSges(4,1) * t473;
t630 = t356 * t673 + t450 + t452;
t628 = t445 + t497;
t627 = -m(6) * pkin(10) - t665;
t44 = t118 * t383 - t378 * t124 - t190 * t498 - t211 * t497;
t626 = -t378 * t43 - t383 * t44;
t62 = mrSges(6,1) * t159 - mrSges(6,3) * t85;
t63 = -mrSges(6,2) * t159 + mrSges(6,3) * t86;
t625 = -t377 * t62 + t382 * t63;
t624 = t12 * t382 - t13 * t377;
t552 = Ifges(3,4) * t379;
t622 = -t379 * (Ifges(3,1) * t384 - t552) / 0.2e1 + pkin(1) * (mrSges(3,1) * t379 + mrSges(3,2) * t384);
t621 = mrSges(5,1) + t629;
t403 = mrSges(5,2) + t627;
t95 = -pkin(4) * t328 - t104;
t78 = -pkin(5) * t193 + t95;
t620 = -mrSges(7,1) * t78 + mrSges(7,3) * t16;
t619 = mrSges(7,2) * t78 - mrSges(7,3) * t15;
t617 = qJ(3) * t482 + t556;
t616 = t44 * mrSges(5,1) - t43 * mrSges(5,2) + Ifges(5,5) * t165 + Ifges(5,6) * t166 + Ifges(5,3) * t273;
t612 = -m(6) * t455 - t556 + t652 + t665 * t383 + (-m(4) + t676) * qJ(3);
t611 = Ifges(7,4) * t608 + Ifges(7,2) * t607 + Ifges(7,6) * t586;
t609 = Ifges(7,1) * t608 + Ifges(7,4) * t607 + Ifges(7,5) * t586;
t36 = Ifges(6,4) * t85 + Ifges(6,2) * t86 + Ifges(6,6) * t159;
t606 = t36 / 0.2e1;
t605 = Ifges(6,1) * t598 + Ifges(6,4) * t597 + Ifges(6,5) * t584;
t546 = Ifges(7,4) * t110;
t53 = Ifges(7,2) * t454 + Ifges(7,6) * t251 + t546;
t604 = -t53 / 0.2e1;
t603 = t53 / 0.2e1;
t103 = Ifges(7,4) * t454;
t54 = Ifges(7,1) * t110 + Ifges(7,5) * t251 + t103;
t602 = -t54 / 0.2e1;
t601 = t54 / 0.2e1;
t600 = -Ifges(5,4) * t165 / 0.2e1 - Ifges(5,2) * t166 / 0.2e1 - Ifges(5,6) * t273 / 0.2e1;
t549 = Ifges(6,4) * t194;
t92 = Ifges(6,2) * t193 + Ifges(6,6) * t260 + t549;
t596 = -t92 / 0.2e1;
t595 = t92 / 0.2e1;
t192 = Ifges(6,4) * t193;
t93 = Ifges(6,1) * t194 + Ifges(6,5) * t260 + t192;
t594 = -t93 / 0.2e1;
t593 = t93 / 0.2e1;
t590 = -t454 / 0.2e1;
t588 = -t110 / 0.2e1;
t581 = -t193 / 0.2e1;
t580 = t193 / 0.2e1;
t579 = -t194 / 0.2e1;
t578 = t194 / 0.2e1;
t576 = -t251 / 0.2e1;
t574 = -t260 / 0.2e1;
t573 = t260 / 0.2e1;
t572 = -t263 / 0.2e1;
t570 = t264 / 0.2e1;
t561 = pkin(5) * t194;
t228 = -t300 * t377 + t374 * t515;
t560 = pkin(5) * t228;
t555 = -Ifges(4,6) - Ifges(3,4);
t554 = mrSges(6,3) * t193;
t553 = mrSges(6,3) * t194;
t551 = Ifges(5,4) * t378;
t550 = Ifges(5,4) * t383;
t548 = Ifges(6,4) * t377;
t547 = Ifges(6,4) * t382;
t545 = Ifges(4,6) * t379;
t544 = Ifges(4,6) * t384;
t538 = t263 * Ifges(5,6);
t536 = t264 * Ifges(5,5);
t535 = t328 * Ifges(5,3);
t40 = -pkin(4) * t273 - t44;
t528 = t383 * t40;
t380 = sin(qJ(1));
t522 = t374 * t380;
t385 = cos(qJ(1));
t520 = t374 * t385;
t516 = t379 * t380;
t514 = t379 * t385;
t513 = t380 * t384;
t510 = t384 * t385;
t303 = t375 * t513 + t514;
t231 = t303 * t378 + t383 * t522;
t304 = -t375 * t516 + t510;
t173 = -t231 * t370 + t304 * t371;
t174 = t231 * t371 + t304 * t370;
t509 = t173 * mrSges(7,1) - t174 * mrSges(7,2);
t301 = -t375 * t510 + t516;
t234 = -t301 * t378 + t383 * t520;
t302 = t375 * t514 + t513;
t508 = (t234 * t370 + t302 * t371) * mrSges(7,1) + (t234 * t371 - t302 * t370) * mrSges(7,2);
t507 = (-t300 * t370 + t371 * t523) * mrSges(7,1) + (-t300 * t371 - t370 * t523) * mrSges(7,2);
t505 = t385 * pkin(1) + pkin(8) * t522;
t501 = qJD(1) ^ 2 * t374 ^ 2;
t485 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t484 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t477 = t304 * pkin(2) + t505;
t470 = t374 * t500;
t461 = -t498 / 0.2e1;
t459 = -t493 / 0.2e1;
t457 = -pkin(1) * t380 + pkin(8) * t520;
t225 = t292 * mrSges(4,1) + t354 * mrSges(4,2);
t72 = -t126 * t377 + t382 * t142;
t134 = t214 * t383 - t378 * t238;
t453 = t591 * t523;
t449 = mrSges(3,3) * t472;
t447 = pkin(3) * t522 + t477;
t446 = t378 * t473;
t438 = -t302 * pkin(2) + t457;
t434 = mrSges(5,1) * t299 + mrSges(5,2) * t300;
t229 = t300 * t382 + t377 * t523;
t431 = mrSges(6,1) * t228 - mrSges(6,2) * t229;
t427 = mrSges(4,2) * t384 - mrSges(4,3) * t379;
t426 = Ifges(5,1) * t378 + t550;
t425 = Ifges(6,1) * t382 - t548;
t424 = Ifges(6,1) * t377 + t547;
t423 = Ifges(5,2) * t383 + t551;
t422 = -Ifges(6,2) * t377 + t547;
t421 = Ifges(6,2) * t382 + t548;
t420 = Ifges(5,5) * t378 + Ifges(5,6) * t383;
t419 = Ifges(6,5) * t382 - Ifges(6,6) * t377;
t418 = Ifges(6,5) * t377 + Ifges(6,6) * t382;
t58 = pkin(5) * t299 - pkin(11) * t229 + t72;
t65 = pkin(11) * t228 + t73;
t25 = -t376 * t65 + t381 * t58;
t26 = t376 * t58 + t381 * t65;
t416 = t104 * t378 - t105 * t383;
t140 = t228 * t381 - t229 * t376;
t141 = t228 * t376 + t229 * t381;
t180 = -t231 * t377 + t304 * t382;
t283 = pkin(8) * t473 - t350;
t351 = t384 * t487;
t285 = -pkin(8) * t471 + t351;
t347 = pkin(2) * t471;
t207 = t347 + (qJD(2) * t417 - t499) * t374;
t247 = (t521 * t591 + t366) * qJD(2);
t80 = -t378 * t207 - t214 * t498 - t238 * t497 + t247 * t383;
t411 = t664 + t8;
t232 = t301 * t383 + t378 * t520;
t409 = t95 * t429;
t407 = t379 * (-Ifges(4,2) * t384 + t545);
t406 = t384 * (Ifges(4,3) * t379 - t544);
t368 = t375 * qJD(3);
t213 = -qJD(2) * t453 + t351 + t368;
t226 = -qJD(4) * t299 + t378 * t471;
t227 = -qJD(4) * t479 + t375 * t497 - t383 * t471;
t100 = pkin(4) * t227 - pkin(10) * t226 + t213;
t79 = t383 * t207 + t214 * t497 - t238 * t498 + t378 * t247;
t70 = pkin(10) * t470 + t79;
t23 = t377 * t100 - t126 * t495 + t142 * t494 + t382 * t70;
t125 = -pkin(4) * t523 - t134;
t186 = -pkin(2) * t354 + t396;
t394 = t196 * mrSges(3,1) - t195 * mrSges(3,2) + t186 * mrSges(4,2) - t164 * mrSges(4,3);
t71 = -pkin(4) * t470 - t80;
t24 = -qJD(5) * t73 + t382 * t100 - t377 * t70;
t392 = (-t377 * t61 - t382 * t60) * qJD(5) + t624;
t391 = -qJD(4) * t416 - t626;
t344 = Ifges(3,4) * t472;
t334 = Ifges(4,1) * t354;
t333 = Ifges(3,3) * t354;
t312 = (t592 + t559) * t383;
t306 = t375 * t562 - t357;
t305 = (-mrSges(3,1) * t384 + mrSges(3,2) * t379) * t374;
t289 = t414 * t378;
t287 = t314 * t378;
t286 = t307 * qJD(2);
t282 = -qJ(3) * t472 + t345;
t281 = t427 * t504;
t278 = -mrSges(3,2) * t356 + t449;
t272 = Ifges(4,4) * t292;
t271 = Ifges(3,5) * t292;
t270 = Ifges(4,5) * t291;
t269 = Ifges(3,6) * t291;
t268 = t375 * t475 + t357;
t261 = t310 + t478;
t259 = -t285 - t368;
t253 = (-pkin(2) * t384 + t458) * t504;
t248 = t347 + (-qJ(3) * t500 - t499) * t374;
t245 = -qJD(1) * t453 + t350;
t243 = -t337 - t284;
t242 = t356 * Ifges(4,4) + (-Ifges(4,2) * t379 - t544) * t504;
t241 = t356 * Ifges(4,5) + (-t384 * Ifges(4,3) - t545) * t504;
t240 = Ifges(3,1) * t473 + t356 * Ifges(3,5) + t344;
t239 = t356 * Ifges(3,6) + (t384 * Ifges(3,2) + t552) * t504;
t236 = -pkin(2) * t356 + qJD(3) + t283;
t230 = -t303 * t383 + t378 * t522;
t224 = mrSges(4,1) * t291 - mrSges(4,3) * t354;
t204 = -mrSges(5,2) * t328 + mrSges(5,3) * t263;
t181 = t231 * t382 + t304 * t377;
t175 = pkin(2) * t291 + t393;
t168 = t414 * t263;
t167 = t314 * t263;
t156 = t264 * Ifges(5,1) + t258 + t646;
t154 = t535 + t536 + t538;
t137 = mrSges(6,1) * t260 - t553;
t136 = -mrSges(6,2) * t260 + t554;
t122 = qJD(5) * t228 + t226 * t382 + t377 * t470;
t121 = -qJD(5) * t229 - t226 * t377 + t382 * t470;
t115 = -mrSges(5,2) * t273 + mrSges(5,3) * t166;
t94 = t125 - t560;
t89 = mrSges(7,1) * t251 - mrSges(7,3) * t110;
t88 = -mrSges(7,2) * t251 + mrSges(7,3) * t454;
t87 = -mrSges(5,1) * t166 + mrSges(5,2) * t165;
t64 = -mrSges(7,1) * t454 + mrSges(7,2) * t110;
t51 = -pkin(5) * t121 + t71;
t48 = -qJD(6) * t141 + t121 * t381 - t122 * t376;
t47 = qJD(6) * t140 + t121 * t376 + t122 * t381;
t22 = -pkin(5) * t86 + t40;
t21 = -mrSges(7,2) * t151 + mrSges(7,3) * t32;
t20 = mrSges(7,1) * t151 - mrSges(7,3) * t31;
t19 = t381 * t45 - t534;
t18 = -t376 * t45 - t531;
t17 = pkin(11) * t121 + t23;
t14 = pkin(5) * t227 - pkin(11) * t122 + t24;
t11 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t5 = -qJD(6) * t26 + t14 * t381 - t17 * t376;
t4 = qJD(6) * t25 + t14 * t376 + t17 * t381;
t1 = [(Ifges(7,1) * t47 + Ifges(7,4) * t48) * t587 + (Ifges(6,4) * t229 + Ifges(6,2) * t228) * t597 + (Ifges(6,4) * t122 + Ifges(6,2) * t121) * t580 + (t140 * t2 - t141 * t3 - t15 * t47 + t16 * t48) * mrSges(7,3) + (t12 * t228 + t121 * t61 - t122 * t60 - t13 * t229) * mrSges(6,3) + (-mrSges(5,3) * t44 + 0.2e1 * t599) * t300 + (Ifges(7,1) * t141 + Ifges(7,4) * t140) * t608 + (Ifges(6,1) * t229 + Ifges(6,4) * t228) * t598 + (Ifges(6,1) * t122 + Ifges(6,4) * t121) * t578 + (Ifges(7,4) * t141 + Ifges(7,2) * t140) * t607 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t589 + (-t671 / 0.2e1 + Ifges(6,3) * t573 + Ifges(6,5) * t578 + Ifges(6,6) * t580 - Ifges(5,4) * t570 + t660 + t681) * t227 + (Ifges(7,5) * t141 + Ifges(7,6) * t140) * t586 + t259 * t279 + t248 * t281 + (t271 / 0.2e1 - t269 / 0.2e1 + t333 / 0.2e1 + t334 / 0.2e1 - t272 / 0.2e1 + t270 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t354 + t485 * t292 + t484 * t291 + t394) * t375 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t575 + t267 * (-mrSges(4,2) * t291 - mrSges(4,3) * t292) + t285 * t278 + t48 * t603 + t229 * t605 + t228 * t606 + t141 * t609 + t140 * t611 + m(6) * (t12 * t73 + t125 * t40 + t13 * t72 + t23 * t61 + t24 * t60 + t71 * t95) + m(7) * (t15 * t5 + t16 * t4 + t2 * t26 + t22 * t94 + t25 * t3 + t51 * t78) + m(5) * (t104 * t80 + t105 * t79 + t120 * t237 + t134 * t44 + t135 * t43 + t203 * t213) + m(4) * (t164 * t266 + t175 * t267 + t186 * t268 + t236 * t286 + t243 * t259 + t248 * t253) + (t258 / 0.2e1 + t156 / 0.2e1 + Ifges(5,1) * t570 + t677) * t226 + m(3) * (t195 * t307 + t196 * t306 + t283 * t286 + t284 * t285) + (Ifges(6,5) * t229 + Ifges(6,6) * t228) * t584 + (Ifges(6,5) * t122 + Ifges(6,6) * t121) * t573 + Ifges(2,3) * qJDD(1) + t266 * t224 + t268 * t225 + t237 * t87 + t213 * t188 + t79 * t204 + t80 * t205 + t22 * (-mrSges(7,1) * t140 + mrSges(7,2) * t141) + t134 * t114 + t135 * t115 + t23 * t136 + t24 * t137 + t95 * (-mrSges(6,1) * t121 + mrSges(6,2) * t122) + t125 * t41 + t71 * t112 + t4 * t88 + t5 * t89 + t94 * t11 + t78 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + t72 * t62 + t73 * t63 + t51 * t64 + t25 * t20 + t26 * t21 + ((-mrSges(3,1) * t291 - mrSges(3,2) * t292 + (m(3) * t564 - t305) * qJDD(1)) * pkin(1) + (-t164 * mrSges(4,1) + t175 * mrSges(4,2) + t195 * mrSges(3,3) - t650 * t354 - t555 * t292 + (-Ifges(3,2) - Ifges(4,3)) * t291) * t384 + (-t196 * mrSges(3,3) + t186 * mrSges(4,1) - t175 * mrSges(4,3) + t651 * t354 + (Ifges(3,1) + Ifges(4,2)) * t292 + t555 * t291 + t616) * t379 + ((-t284 * mrSges(3,3) + t243 * mrSges(4,1) - t239 / 0.2e1 + t241 / 0.2e1 - t253 * mrSges(4,2) + t484 * t356) * t379 + (t283 * mrSges(3,3) + t236 * mrSges(4,1) + t536 / 0.2e1 + t538 / 0.2e1 - t105 * mrSges(5,2) + t240 / 0.2e1 - t242 / 0.2e1 + t154 / 0.2e1 - t253 * mrSges(4,3) + t104 * mrSges(5,1) + t535 / 0.2e1 + t485 * t356) * t384 + (t384 * (Ifges(3,4) * t384 - Ifges(3,2) * t379) / 0.2e1 - t407 / 0.2e1 - t406 / 0.2e1 - t622) * t504) * qJD(2) + (mrSges(3,3) + mrSges(4,1)) * (-g(1) * t385 - g(2) * t380)) * t374 + (-m(3) * t457 - m(4) * t438 + mrSges(2,1) * t380 + mrSges(2,2) * t385 + t617 * t301 + t403 * t232 - t621 * t234 + (-t377 * t638 + t428 + t530 + t618) * t302 + t492 * (-pkin(3) * t520 - t438)) * g(1) - t40 * t431 + t120 * t434 + t47 * t601 + t122 * t593 + t121 * t595 + (-m(3) * t505 - m(4) * t477 - mrSges(2,1) * t385 + mrSges(2,2) * t380 - m(7) * (t231 * t367 + t447) - t174 * mrSges(7,1) - t173 * mrSges(7,2) - m(6) * (pkin(4) * t231 + t447) - t181 * mrSges(6,1) - t180 * mrSges(6,2) - m(5) * t447 - t231 * mrSges(5,1) - t617 * t303 + t403 * t230 + (-t490 - t618) * t304) * g(2) + t306 * (mrSges(3,1) * t354 - mrSges(3,3) * t292) + t307 * (-mrSges(3,2) * t354 - mrSges(3,3) * t291) + (-t43 * mrSges(5,3) - Ifges(5,4) * t583 + Ifges(6,5) * t598 + Ifges(7,5) * t608 - Ifges(5,2) * t582 - Ifges(5,6) * t569 + Ifges(6,6) * t597 + Ifges(7,6) * t607 + Ifges(6,3) * t584 + Ifges(7,3) * t586 + t600 + t663 + t664 + t675) * t299 + t630 * t286; -t639 * t383 * t592 + t328 * (mrSges(5,1) * t383 - mrSges(5,2) * t378) * t203 + (-m(4) * t243 + t278 - t279 - t449) * t283 + ((t383 * t647 + t239) * t379 + t384 * t242) * t504 / 0.2e1 + (-t224 + t87) * qJ(3) - ((-Ifges(3,2) * t473 + t154 + t240 + t344) * t384 + (t383 * t155 + t378 * t156 + t241) * t379 + (t379 * t650 + t384 * t651) * t356 + t328 * (Ifges(5,3) * t384 + t379 * t420) + t264 * (Ifges(5,5) * t384 + t379 * t426) + t263 * (Ifges(5,6) * t384 + t379 * t423)) * t504 / 0.2e1 + (t120 * qJ(3) - t104 * t148 - t105 * t149 - t592 * t391 + (qJD(3) - t245) * t203) * m(5) + (-t130 * t95 + t12 * t262 + t13 * t261 - (t498 * t95 - t528) * t592 + t642 * t61 + t641 * t60) * m(6) + (-t149 - t467) * t204 + (t468 - t148) * t205 + (-pkin(2) * t186 - qJ(3) * t164 - qJD(3) * t243 - t253 * t282) * m(4) + (t661 + t660) * t497 - t282 * t281 - t243 * t452 - t236 * t451 + t2 * (-mrSges(7,2) * t378 - mrSges(7,3) * t288) + (-Ifges(7,4) * t290 - Ifges(7,2) * t288 + Ifges(7,6) * t378) * t607 + (-Ifges(7,1) * t290 - Ifges(7,4) * t288 + Ifges(7,5) * t378) * t608 + (-Ifges(7,5) * t290 - Ifges(7,6) * t288 + Ifges(7,3) * t378) * t586 + t22 * (mrSges(7,1) * t288 - mrSges(7,2) * t290) + t3 * (mrSges(7,1) * t378 + mrSges(7,3) * t290) - (t263 * t423 + t264 * t426 + t328 * t420) * qJD(4) / 0.2e1 + (t301 * t653 + t302 * t612) * g(2) + (t303 * t653 + t304 * t612) * g(1) + t334 - t36 * t519 / 0.2e1 + t271 - t272 + t333 - t269 + t270 + (-m(4) * t506 + t305 - t492 * (t360 + t506) + (t427 + (-mrSges(5,3) + t679) * t384 + (-t383 * t627 + t652) * t379) * t374) * g(3) + t147 * t603 + t171 * t604 + t512 * t605 - t290 * t609 - t288 * t611 + (t407 + t406) * t501 / 0.2e1 + (Ifges(7,5) * t145 + Ifges(7,6) * t147) * t575 + t657 * t112 + (Ifges(7,4) * t145 + Ifges(7,2) * t147) * t589 + (Ifges(7,1) * t145 + Ifges(7,4) * t147) * t587 + t394 + t156 * t461 + t261 * t62 + t262 * t63 - t245 * t188 - pkin(2) * t225 + t128 * t20 + t129 * t21 + (-t105 * (mrSges(5,3) * t379 * t383 - mrSges(5,2) * t384) - t253 * (-mrSges(4,2) * t379 - mrSges(4,3) * t384) - t104 * (mrSges(5,1) * t384 - mrSges(5,3) * t518)) * t504 + t648 * t89 + (t128 * t3 + t129 * t2 + t15 * t648 + t16 * t649 + t22 * t312 + t640 * t78) * m(7) + t649 * t88 + t640 * t64 + t641 * t137 + t642 * t136 + t120 * t432 + (Ifges(6,6) * t378 + t383 * t422) * t597 + (Ifges(6,5) * t378 + t383 * t425) * t598 + t383 * t599 + t378 * t600 + t145 * t601 + t172 * t602 + (Ifges(7,4) * t172 + Ifges(7,2) * t171 - Ifges(7,6) * t445) * t590 + t257 * t594 + t469 * t595 + t256 * t596 + (-Ifges(5,2) * t378 + t550) * t582 + (Ifges(5,1) * t383 - t551) * t583 + (Ifges(6,3) * t378 + t383 * t419) * t584 + (Ifges(7,1) * t172 + Ifges(7,4) * t171 - Ifges(7,5) * t445) * t588 + (-t418 * t493 + (Ifges(6,3) * t383 - t378 * t419) * qJD(4)) * t573 + (Ifges(6,5) * t257 + Ifges(6,6) * t256 - Ifges(6,3) * t445) * t574 + (Ifges(7,5) * t172 + Ifges(7,6) * t171 - Ifges(7,3) * t445) * t576 + (-t424 * t493 + (Ifges(6,5) * t383 - t378 * t425) * qJD(4)) * t578 + (Ifges(6,1) * t257 + Ifges(6,4) * t256 - Ifges(6,5) * t445) * t579 + (-t421 * t493 + (Ifges(6,6) * t383 - t378 * t422) * qJD(4)) * t580 + (Ifges(6,4) * t257 + Ifges(6,2) * t256 - Ifges(6,6) * t445) * t581 + (Ifges(5,5) * t383 - Ifges(5,6) * t378) * t569 + t498 * t542 + t429 * t528 + t377 * t93 * t459 + t312 * t11 + t378 * t675 - t115 * t517 + t622 * t501 + t626 * mrSges(5,3) + (-m(4) * t236 + t450 - t630) * t284 + t631 * qJD(3) + (-t12 * t378 - t445 * t61 - t632 * t95) * mrSges(6,2) + (t13 * t378 + t445 * t60 + t633 * t95) * mrSges(6,1) + (-t12 * t519 - t13 * t512 + t60 * t632 - t61 * t633) * mrSges(6,3) + (-mrSges(7,2) * t628 + mrSges(7,3) * t635) * t16 + (-mrSges(7,1) * t635 + mrSges(7,2) * t636) * t78 + (mrSges(7,1) * t628 - mrSges(7,3) * t636) * t15 + (t459 * t92 + t461 * t93) * t382; t655 * t88 + (t204 * t473 - t11 + (t136 * t382 - t137 * t377 + t204) * qJD(4) + t639) * t383 + t281 * t473 - t287 * t20 - t289 * t21 + (t115 + (-t377 * t136 - t382 * t137) * qJD(5) + t328 * (t64 - t634) + t625) * t378 + t656 * t89 - t631 * t356 + t225 - t254 * t137 - t255 * t136 + (-g(1) * t303 - g(2) * t301 + g(3) * t521) * t482 + (-t2 * t289 - t22 * t383 - t287 * t3 + (t446 + t498) * t78 + t655 * t16 + t656 * t15) * m(7) + ((-t40 + (-t377 * t60 + t382 * t61) * qJD(4)) * t383 + (qJD(4) * t95 + t392) * t378 - t254 * t60 - t255 * t61 + t446 * t95) * m(6) + (-t203 * t356 - t416 * t473 + t391) * m(5) + (t243 * t356 + t253 * t473 + t186) * m(4); (-pkin(4) * t40 - t105 * t95 - t60 * t74 - t61 * t75) * m(6) + (-t15 * t168 + t16 * t167 - t2 * t414 - t3 * t314) * mrSges(7,3) + (Ifges(7,4) * t314 - Ifges(7,2) * t414) * t607 + (Ifges(7,1) * t314 - Ifges(7,4) * t414) * t608 + (Ifges(7,5) * t314 - Ifges(7,6) * t414) * t586 + t22 * (mrSges(7,1) * t414 + mrSges(7,2) * t314) - t414 * t611 + (t258 + t156) * t572 + (Ifges(6,5) * t579 + Ifges(7,5) * t588 - Ifges(5,2) * t572 + Ifges(6,6) * t581 + Ifges(7,6) * t590 + Ifges(6,3) * t574 + Ifges(7,3) * t576 - t681) * t264 + (-t666 * t61 + (-t494 + t524) * t60 + t624) * mrSges(6,3) + qJD(5) * t409 + (t193 * t422 + t194 * t425 + t260 * t419) * qJD(5) / 0.2e1 + (-Ifges(7,4) * t168 - Ifges(7,2) * t167) * t590 + (-Ifges(7,1) * t168 - Ifges(7,4) * t167) * t588 + (-Ifges(7,5) * t168 - Ifges(7,6) * t167) * t576 - t78 * (mrSges(7,1) * t167 - mrSges(7,2) * t168) - t167 * t604 + t377 * t605 + t382 * t606 + t314 * t609 - (Ifges(7,1) * t587 + Ifges(7,4) * t589 + Ifges(7,5) * t575 + t601 + t619) * t390 + (t419 * t574 + t422 * t581 + t425 * t579 - t409 - t677) * t263 + t616 + t249 * t20 + t250 * t21 - t104 * t204 - t75 * t136 - t74 * t137 - pkin(4) * t41 - (Ifges(5,1) * t263 - t537 + t647) * t264 / 0.2e1 + t40 * t430 + t643 * t89 + t644 * t88 + (t15 * t643 + t16 * t644 + t2 * t250 - t22 * t367 + t249 * t3 + t637 * t78) * m(7) + t421 * t597 + t424 * t598 - t168 * t602 + t494 * t593 + t524 * t594 + t525 * t595 + t495 * t596 + t418 * t584 + t155 * t570 - (Ifges(7,4) * t587 + Ifges(7,2) * t589 + Ifges(7,6) * t575 + t603 + t620) * t218 - t367 * t11 + (t230 * t621 + t231 * t403) * g(1) + (-t232 * t621 - t234 * t403) * g(2) + (m(6) * t392 - t136 * t495 - t137 * t494 + t625) * pkin(10) + (t299 * t629 + t300 * t627 + t434) * g(3) + t634 * t105 + t637 * t64; (t137 + t553) * t61 + (-t136 + t554) * t60 + (-Ifges(6,2) * t194 + t192 + t93) * t581 + (-m(7) * t560 - t431 - t507) * g(3) + ((-t376 * t89 + t381 * t88) * qJD(6) + t20 * t381 + t21 * t376) * pkin(5) - t64 * t561 - m(7) * (t15 * t18 + t16 * t19 + t561 * t78) + t35 + (t2 * t376 + t3 * t381 + (-t15 * t376 + t16 * t381) * qJD(6)) * t610 + t411 + t663 - t95 * (mrSges(6,1) * t194 + mrSges(6,2) * t193) - t19 * t88 - t18 * t89 + (-t508 - (t234 * t382 - t302 * t377) * mrSges(6,2) + t638 * (t234 * t377 + t302 * t382)) * g(2) + (mrSges(6,2) * t181 + t180 * t638 - t509) * g(1) + (Ifges(6,5) * t193 - Ifges(6,6) * t194) * t574 + t92 * t578 + (Ifges(6,1) * t193 - t549) * t579 + (Ifges(7,1) * t588 + Ifges(7,4) * t590 + Ifges(7,5) * t576 + t602 - t619) * t454 - (Ifges(7,4) * t588 + Ifges(7,2) * t590 + Ifges(7,6) * t576 + t604 - t620) * t110; -t78 * (mrSges(7,1) * t110 + mrSges(7,2) * t454) + (Ifges(7,1) * t454 - t546) * t588 + t53 * t587 + (Ifges(7,5) * t454 - Ifges(7,6) * t110) * t576 - t15 * t88 + t16 * t89 - g(1) * t509 - g(2) * t508 - g(3) * t507 + (t110 * t16 + t15 * t454) * mrSges(7,3) + t411 + (-Ifges(7,2) * t110 + t103 + t54) * t590;];
tau  = t1;
