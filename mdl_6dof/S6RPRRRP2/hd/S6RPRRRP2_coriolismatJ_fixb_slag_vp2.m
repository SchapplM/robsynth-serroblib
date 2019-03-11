% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:59
% EndTime: 2019-03-09 05:59:21
% DurationCPUTime: 13.43s
% Computational Cost: add. (20276->727), mult. (42197->950), div. (0->0), fcn. (41502->8), ass. (0->376)
t752 = -mrSges(7,3) / 0.2e1;
t418 = sin(qJ(5));
t419 = sin(qJ(4));
t421 = cos(qJ(4));
t642 = cos(qJ(5));
t495 = t642 * t421;
t372 = -t418 * t419 + t495;
t373 = -t418 * t421 - t419 * t642;
t259 = -mrSges(7,1) * t372 - mrSges(7,2) * t373;
t634 = pkin(4) * t421;
t407 = -pkin(3) - t634;
t312 = -pkin(5) * t372 + t407;
t639 = m(7) * t312;
t751 = t259 + t639;
t736 = Ifges(6,6) + Ifges(7,6);
t737 = Ifges(6,5) + Ifges(7,5);
t480 = t372 * t737 + t373 * t736;
t682 = -pkin(9) - pkin(8);
t390 = t682 * t419;
t391 = t682 * t421;
t714 = t642 * t390 + t418 * t391;
t728 = t714 * mrSges(6,2);
t726 = t373 * qJ(6) + t714;
t733 = t726 * mrSges(7,2);
t281 = t418 * t390 - t642 * t391;
t734 = t281 * mrSges(6,1);
t212 = t372 * qJ(6) + t281;
t747 = t212 * mrSges(7,1);
t750 = t480 - t728 - t733 - t734 - t747;
t749 = -t728 / 0.2e1 - t733 / 0.2e1 - t734 / 0.2e1 - t747 / 0.2e1;
t420 = sin(qJ(3));
t345 = t373 * t420;
t565 = t345 * qJ(6);
t422 = cos(qJ(3));
t496 = -cos(pkin(10)) * pkin(1) - pkin(2);
t361 = -pkin(3) * t422 - t420 * pkin(8) + t496;
t348 = t421 * t361;
t403 = sin(pkin(10)) * pkin(1) + pkin(7);
t552 = t420 * t421;
t529 = pkin(9) * t552;
t202 = -t529 + t348 + (-t403 * t419 - pkin(4)) * t422;
t387 = t422 * t403;
t249 = t419 * t361 + t387 * t421;
t554 = t419 * t420;
t229 = -pkin(9) * t554 + t249;
t207 = t642 * t229;
t86 = t418 * t202 + t207;
t65 = t86 + t565;
t248 = -t387 * t419 + t348;
t228 = t248 - t529;
t90 = -t228 * t418 - t207;
t75 = t90 - t565;
t748 = t65 + t75;
t633 = pkin(5) * t422;
t344 = t418 * t554 - t420 * t495;
t313 = t344 * qJ(6);
t205 = t418 * t229;
t85 = t642 * t202 - t205;
t64 = t313 + t85;
t61 = -t633 + t64;
t592 = t212 * t61;
t746 = Ifges(7,3) + Ifges(6,3);
t91 = t642 * t228 - t205;
t745 = t85 - t91;
t528 = t642 * pkin(4);
t406 = t528 + pkin(5);
t347 = t372 * t422;
t469 = t420 * mrSges(7,1) - t347 * mrSges(7,3);
t471 = t420 * mrSges(6,1) - t347 * mrSges(6,3);
t484 = t528 / 0.2e1;
t636 = pkin(4) * t418;
t522 = t636 / 0.2e1;
t549 = t421 * t422;
t553 = t419 * t422;
t631 = t420 * pkin(3);
t632 = pkin(8) * t422;
t392 = t631 - t632;
t283 = t419 * t392 - t403 * t552;
t585 = t283 * mrSges(5,2);
t386 = t420 * t403;
t282 = t419 * t386 + t421 * t392;
t586 = t282 * mrSges(5,1);
t648 = t420 / 0.2e1;
t689 = m(6) * pkin(4);
t691 = m(7) / 0.2e1;
t649 = -t420 / 0.2e1;
t738 = -t347 / 0.2e1;
t234 = t420 * pkin(4) - pkin(9) * t549 + t282;
t246 = -pkin(9) * t553 + t283;
t96 = t642 * t234 - t246 * t418;
t77 = pkin(5) * t420 - qJ(6) * t347 + t96;
t740 = -t77 / 0.2e1;
t346 = t373 * t422;
t97 = t418 * t234 + t642 * t246;
t81 = qJ(6) * t346 + t97;
t699 = t97 * mrSges(6,2) / 0.2e1 - t96 * mrSges(6,1) / 0.2e1 + t81 * mrSges(7,2) / 0.2e1 + mrSges(7,1) * t740 + t737 * t738 + t746 * t649;
t739 = t346 / 0.2e1;
t698 = t736 * t739 - t699;
t467 = -mrSges(7,2) * t420 + mrSges(7,3) * t346;
t468 = -t420 * mrSges(6,2) + t346 * mrSges(6,3);
t712 = t468 + t467;
t744 = (t406 * t77 + t636 * t81) * t691 + Ifges(5,3) * t648 + t586 / 0.2e1 - t585 / 0.2e1 + t406 * t469 / 0.2e1 + (t418 * t97 + t642 * t96) * t689 / 0.2e1 + Ifges(5,5) * t549 / 0.2e1 - Ifges(5,6) * t553 / 0.2e1 + t471 * t484 + t712 * t522 + t698;
t287 = -mrSges(6,1) * t422 + t344 * mrSges(6,3);
t673 = -t287 / 0.2e1;
t286 = -mrSges(7,1) * t422 + t344 * mrSges(7,3);
t675 = -t286 / 0.2e1;
t730 = t675 + t673;
t363 = t372 * mrSges(7,2);
t257 = -t373 * mrSges(7,1) + t363;
t258 = -mrSges(6,1) * t373 + mrSges(6,2) * t372;
t368 = Ifges(7,4) * t372;
t261 = Ifges(7,2) * t373 + t368;
t369 = Ifges(6,4) * t372;
t263 = Ifges(6,2) * t373 + t369;
t266 = -Ifges(7,1) * t373 + t368;
t268 = -Ifges(6,1) * t373 + t369;
t741 = t312 * t257 + t407 * t258 + (t261 / 0.2e1 + t263 / 0.2e1 + t266 / 0.2e1 + t268 / 0.2e1) * t372;
t643 = t422 / 0.2e1;
t721 = t86 + t90;
t537 = t419 ^ 2 + t421 ^ 2;
t735 = mrSges(5,3) * t537;
t477 = t528 - t406;
t214 = -t344 * mrSges(6,1) + t345 * mrSges(6,2);
t352 = pkin(4) * t554 + t386;
t230 = -pkin(5) * t345 + t352;
t644 = -t422 / 0.4e1;
t651 = t407 / 0.2e1;
t583 = t345 * mrSges(6,3);
t285 = mrSges(6,2) * t422 + t583;
t676 = t285 / 0.2e1;
t582 = t345 * mrSges(7,3);
t284 = mrSges(7,2) * t422 + t582;
t677 = t284 / 0.2e1;
t679 = t258 / 0.2e1;
t680 = t257 / 0.2e1;
t315 = t345 * mrSges(7,2);
t584 = t344 * mrSges(7,1);
t489 = -t315 + t584;
t681 = -t489 / 0.2e1;
t731 = t214 * t651 + t230 * t680 + t281 * t673 + t312 * t681 + t352 * t679 + t480 * t644 + t714 * t676 + t726 * t677;
t729 = mrSges(6,3) + mrSges(7,3);
t725 = (t257 + t258) * t643;
t608 = Ifges(7,4) * t373;
t265 = Ifges(7,1) * t372 + t608;
t610 = Ifges(6,4) * t373;
t267 = Ifges(6,1) * t372 + t610;
t686 = mrSges(7,3) / 0.2e1;
t687 = mrSges(6,3) / 0.2e1;
t722 = t212 * t686 + t281 * t687 - t265 / 0.4e1 - t267 / 0.4e1;
t412 = Ifges(5,5) * t421;
t604 = Ifges(5,6) * t419;
t720 = Ifges(4,4) - t412 / 0.2e1 + t604 / 0.2e1;
t262 = Ifges(7,2) * t372 - t608;
t264 = Ifges(6,2) * t372 - t610;
t719 = t264 + t262;
t718 = t268 + t266;
t716 = t284 + t285;
t413 = Ifges(5,4) * t421;
t713 = -Ifges(5,2) * t419 + t413;
t384 = Ifges(5,1) * t419 + t413;
t710 = (mrSges(6,2) + mrSges(7,2)) * t642;
t216 = -mrSges(6,1) * t345 - mrSges(6,2) * t344;
t360 = t420 * t384;
t709 = pkin(4) * t216 / 0.2e1 - t360 / 0.4e1;
t316 = Ifges(7,6) * t344;
t317 = Ifges(6,6) * t344;
t318 = Ifges(7,5) * t345;
t319 = Ifges(6,5) * t345;
t481 = t319 + t317 + t318 + t316;
t707 = -t282 * t419 + t283 * t421;
t706 = -mrSges(5,1) * t421 + mrSges(5,2) * t419;
t704 = t91 / 0.2e1 - t85 / 0.2e1;
t703 = t344 * t729;
t615 = mrSges(7,3) * t372;
t509 = -t615 / 0.2e1;
t702 = (-t212 * t691 + t509) * pkin(5) + t749;
t566 = t344 * t373;
t493 = t566 / 0.2e1;
t240 = mrSges(7,3) * t493;
t494 = -t566 / 0.2e1;
t701 = mrSges(7,3) * t494 + t240 - t725 + (t493 + t494) * mrSges(6,3);
t215 = -mrSges(7,1) * t345 - mrSges(7,2) * t344;
t531 = pkin(4) * t552;
t276 = -pkin(5) * t344 + t531;
t635 = pkin(4) * t419;
t342 = -pkin(5) * t373 + t635;
t647 = -t421 / 0.2e1;
t650 = -t419 / 0.2e1;
t672 = t342 / 0.2e1;
t674 = t286 / 0.2e1;
t678 = t276 / 0.2e1;
t76 = t313 + t91;
t700 = (t212 * t76 + t230 * t342 + t276 * t312 + t726 * t748 - t592) * t691 - t212 * t674 + t259 * t678 + t215 * t672 + ((-mrSges(5,1) * t422 - mrSges(5,3) * t552) * t647 + (mrSges(5,2) * t422 - mrSges(5,3) * t554) * t650 + t649 * t735) * pkin(8) + t731;
t697 = 0.2e1 * m(7);
t695 = 2 * qJD(3);
t694 = -m(6) / 0.2e1;
t693 = m(6) / 0.2e1;
t692 = -m(7) / 0.2e1;
t690 = -pkin(3) / 0.2e1;
t688 = m(7) * pkin(5);
t685 = -t61 / 0.2e1;
t669 = t344 / 0.2e1;
t612 = Ifges(5,4) * t419;
t382 = Ifges(5,2) * t421 + t612;
t359 = t420 * t382;
t661 = -t359 / 0.4e1;
t659 = t372 / 0.2e1;
t575 = t421 * mrSges(5,2);
t576 = t419 * mrSges(5,1);
t381 = t575 + t576;
t654 = t381 / 0.2e1;
t653 = -t382 / 0.4e1;
t385 = Ifges(5,1) * t421 - t612;
t652 = t385 / 0.4e1;
t646 = t421 / 0.2e1;
t645 = -t422 / 0.2e1;
t640 = m(7) * t230;
t638 = m(7) * t344;
t637 = pkin(3) * t381;
t630 = t64 * mrSges(7,2);
t629 = t65 * mrSges(7,1);
t628 = t75 * mrSges(7,1);
t627 = t76 * mrSges(7,2);
t624 = t85 * mrSges(6,2);
t623 = t86 * mrSges(6,1);
t622 = t90 * mrSges(6,1);
t621 = t91 * mrSges(6,2);
t618 = -t61 + t64;
t614 = mrSges(7,3) * t373;
t611 = Ifges(6,4) * t344;
t609 = Ifges(7,4) * t344;
t591 = t212 * t64;
t581 = t346 * mrSges(7,1);
t580 = t347 * mrSges(7,2);
t579 = t372 * mrSges(6,3);
t578 = t373 * mrSges(6,3);
t577 = t406 * mrSges(7,3);
t574 = t422 * Ifges(5,6);
t573 = t706 - mrSges(4,1);
t343 = t345 ^ 2;
t464 = t61 * t344 + t65 * t345;
t485 = t531 / 0.2e1;
t516 = t687 + t686;
t525 = m(7) * t678;
t12 = t516 * t343 + t464 * t692 + (t681 + t214 / 0.2e1 + t525 + m(6) * t485) * t422 + (-t284 / 0.2e1 - t285 / 0.2e1 + t75 * t692 + t721 * t694) * t345 + (t516 * t344 - t76 * t692 + t694 * t745 + t730) * t344;
t572 = t12 * qJD(1);
t564 = t345 * t284;
t568 = t344 * t286;
t13 = t287 * t669 + t564 / 0.2e1 + t568 / 0.2e1 + t345 * t676 + (-t618 + t633) * t638 / 0.2e1 + (-t489 + t214) * t645 + t729 * (-t343 / 0.2e1 - t344 ^ 2 / 0.2e1);
t571 = t13 * qJD(1);
t567 = t344 * t372;
t563 = t345 * t373;
t560 = t352 * t419;
t36 = m(7) * t464 + t564 + t568;
t559 = t36 * qJD(1);
t557 = t406 * t373;
t338 = t420 * t713 - t574;
t555 = t419 * t338;
t551 = t420 * t422;
t340 = -t422 * Ifges(5,5) + t420 * t385;
t550 = t421 * t340;
t402 = pkin(4) * t553;
t353 = t387 + t402;
t536 = qJD(3) * t422;
t534 = t688 / 0.2e1;
t532 = t373 * t636;
t530 = t372 * t636;
t295 = t345 * t636;
t154 = t406 * t344 + t295;
t527 = t154 * t691;
t517 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t515 = -Ifges(6,5) / 0.2e1 - Ifges(7,5) / 0.2e1;
t514 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t513 = mrSges(7,1) + t688;
t508 = t615 / 0.2e1;
t501 = -t714 * mrSges(6,3) / 0.2e1;
t499 = -t582 / 0.2e1;
t498 = t582 / 0.2e1;
t492 = t403 * t654;
t320 = Ifges(7,4) * t345;
t184 = -Ifges(7,1) * t344 - Ifges(7,5) * t422 + t320;
t321 = Ifges(6,4) * t345;
t186 = -Ifges(6,1) * t344 - Ifges(6,5) * t422 + t321;
t490 = t184 / 0.2e1 + t186 / 0.2e1;
t488 = t412 - t604;
t487 = t372 * t528;
t479 = (-mrSges(7,2) / 0.2e1 - mrSges(6,2) / 0.2e1) * t347;
t478 = t514 * t346;
t476 = -t214 - t315;
t475 = t528 * t583;
t474 = t249 * mrSges(5,1) + t248 * mrSges(5,2);
t472 = -t346 * mrSges(6,1) + t347 * mrSges(6,2);
t470 = t580 - t581;
t465 = Ifges(5,5) * t419 + Ifges(5,6) * t421;
t231 = -pkin(5) * t346 + t353;
t10 = (t230 * t420 - t231 * t422 - t81 * t344 + t77 * t345 + t61 * t346 + t65 * t347) * t692 + (-t97 * t344 + t96 * t345 + t85 * t346 + t86 * t347 + t352 * t420 - t353 * t422) * t694 - m(5) * (t249 * t549 + t283 * t552 + (t420 ^ 2 - t422 ^ 2) * t403) / 0.2e1 + m(5) * (-t248 * t422 - t282 * t420) * t650 + t712 * t669 - (t471 + t469) * t345 / 0.2e1 + t346 * t730 + t716 * t738 + (t215 + t216) * t649 + (t470 + t472) * t643;
t181 = Ifges(7,4) * t347 + Ifges(7,2) * t346 + Ifges(7,6) * t420;
t183 = Ifges(6,4) * t347 + Ifges(6,2) * t346 + Ifges(6,6) * t420;
t185 = Ifges(7,1) * t347 + Ifges(7,4) * t346 + Ifges(7,5) * t420;
t187 = Ifges(6,1) * t347 + Ifges(6,4) * t346 + Ifges(6,5) * t420;
t339 = Ifges(5,6) * t420 + t422 * t713;
t341 = Ifges(5,5) * t420 + t385 * t422;
t180 = Ifges(7,2) * t345 - Ifges(7,6) * t422 - t609;
t182 = Ifges(6,2) * t345 - Ifges(6,6) * t422 - t611;
t446 = t65 * mrSges(7,3) + t86 * mrSges(6,3) + t180 / 0.2e1 + t182 / 0.2e1;
t6 = t231 * t215 + t353 * t216 + t81 * t284 + t97 * t285 + t77 * t286 + t96 * t287 + (t181 / 0.2e1 + t183 / 0.2e1) * t345 + (-t185 / 0.2e1 - t187 / 0.2e1) * t344 + m(6) * (t352 * t353 + t85 * t96 + t86 * t97) + m(7) * (t230 * t231 + t61 * t77 + t65 * t81) + m(5) * (t248 * t282 + t249 * t283) + (t352 * mrSges(6,2) + t230 * mrSges(7,2) - t85 * mrSges(6,3) - t61 * mrSges(7,3) + t490) * t347 + (-t352 * mrSges(6,1) - t230 * mrSges(7,1) + t446) * t346 + (-t586 + t585 + t550 / 0.2e1 - t555 / 0.2e1 + t496 * mrSges(4,2) + t515 * t347 + t478 + (-t248 * t421 - t249 * t419) * mrSges(5,3) + t720 * t422) * t422 + (t496 * mrSges(4,1) - t86 * mrSges(6,2) - t65 * mrSges(7,2) + t85 * mrSges(6,1) + t61 * mrSges(7,1) + t341 * t646 + t339 * t650 + t248 * mrSges(5,1) - t249 * mrSges(5,2) - t720 * t420 - t514 * t345 + t515 * t344 + (-t282 * t421 - t283 * t419) * mrSges(5,3) + (Ifges(4,1) - Ifges(4,2) - Ifges(5,3) + (m(5) * t403 + 0.2e1 * t381) * t403 - t746) * t422) * t420;
t463 = t6 * qJD(1) - t10 * qJD(2);
t217 = Ifges(7,2) * t344 + t320;
t218 = Ifges(6,2) * t344 + t321;
t433 = -t230 * t489 + t352 * t214 + (t217 / 0.2e1 + t218 / 0.2e1 + t490) * t345 - t61 * t582 - t85 * t583;
t219 = Ifges(7,1) * t345 + t609;
t220 = Ifges(6,1) * t345 + t611;
t438 = -t219 / 0.2e1 - t220 / 0.2e1 + t446;
t449 = t706 * t420;
t7 = t276 * t215 + t76 * t284 + t91 * t285 + t75 * t286 + t90 * t287 + (-t319 / 0.2e1 - t317 / 0.2e1 - t318 / 0.2e1 - t316 / 0.2e1 + t474) * t422 + m(7) * (t230 * t276 + t61 * t75 + t65 * t76) + m(6) * (t85 * t90 + t86 * t91) + t438 * t344 + (-t360 * t646 + t338 * t647 + t465 * t643 - t403 * t449 + (m(6) * t352 + t216) * t634 + (-t359 + t340) * t650) * t420 + t433;
t462 = t7 * qJD(1) - t12 * qJD(2);
t8 = -t86 * t287 + t85 * t285 + t64 * t284 + (m(7) * t618 - t286) * t65 + ((-t215 - t640) * pkin(5) + t438) * t344 + t433 + t481 * t645;
t461 = t8 * qJD(1) + t13 * qJD(2);
t41 = 0.4e1 * (m(7) / 0.4e1 + m(6) / 0.4e1) * (-t344 * t347 + t345 * t346 - t551) + m(5) * (-0.1e1 + t537) * t551;
t460 = -t10 * qJD(1) + t41 * qJD(2);
t54 = (t154 / 0.4e1 - t276 / 0.4e1) * t697 + t489;
t82 = (t530 / 0.4e1 + t557 / 0.4e1 - t342 / 0.4e1) * t697 - t257;
t459 = qJD(1) * t54 + qJD(3) * t82;
t118 = t344 * t513 - t315;
t179 = t373 * t513 - t363;
t458 = qJD(1) * t118 + qJD(3) * t179;
t456 = -t576 / 0.2e1 - t575 / 0.2e1;
t455 = pkin(5) * t469;
t453 = t180 / 0.4e1 + t182 / 0.4e1 - t219 / 0.4e1 - t220 / 0.4e1;
t452 = t184 / 0.4e1 + t186 / 0.4e1 + t217 / 0.4e1 + t218 / 0.4e1;
t450 = t262 / 0.2e1 + t264 / 0.2e1 - t265 / 0.2e1 - t267 / 0.2e1;
t260 = -mrSges(6,1) * t372 - mrSges(6,2) * t373;
t14 = -t637 + (t713 / 0.2e1 + t384 / 0.2e1) * t421 + (pkin(4) * t260 - t382 / 0.2e1 + t385 / 0.2e1) * t419 + m(6) * t407 * t635 + t450 * t373 + t741 + t751 * t342;
t428 = -t342 * t422 * t692 - t402 * t694;
t296 = t347 * t636;
t429 = t517 * t346 + t479 + (t346 * t528 + t296) * t693 + (t346 * t406 + t296) * t691;
t16 = (t680 + t679 + t654 + t456) * t422 + t428 + t429;
t440 = -t281 * t745 + t721 * t714;
t2 = (t182 + t180) * t373 / 0.4e1 - (t220 + t219) * t373 / 0.4e1 + (t218 + t217 + t186 + t184) * t372 / 0.4e1 + t550 / 0.4e1 - t744 + t748 * t614 / 0.2e1 + t76 * t508 + t61 * t509 + t420 * t492 + t345 * t501 + (t263 + t261 + t718) * t345 / 0.4e1 + t721 * t578 / 0.2e1 - (t384 + t713) * t554 / 0.4e1 + t726 * t499 + (t719 / 0.4e1 + t722) * t344 + t709 * t419 + t704 * t579 - t555 / 0.4e1 + t488 * t644 + t552 * t652 + t552 * t653 + t421 * t661 - t449 * t690 + ((t407 * t552 + t560) * pkin(4) + t440) * t693 + t700 + t260 * t485;
t448 = t2 * qJD(1) - t16 * qJD(2) + t14 * qJD(3);
t15 = (-pkin(5) * t751 + t450) * t373 + t741;
t435 = t479 + (t534 + t517) * t346;
t442 = m(7) * t373 * t633;
t19 = -t442 / 0.2e1 + t435 + t725;
t434 = (t501 + t726 * t752 + t268 / 0.4e1 + t266 / 0.4e1 + t263 / 0.4e1 + t261 / 0.4e1) * t345;
t436 = t264 / 0.4e1 + t262 / 0.4e1 + t722;
t423 = ((-t640 / 0.2e1 - t215 / 0.2e1) * pkin(5) + t453) * t373 + ((t64 / 0.2e1 + t685) * mrSges(7,3) + t452) * t372 + t434 + ((-t639 / 0.2e1 - t259 / 0.2e1) * pkin(5) + t436) * t344 + t212 * t675 + t731;
t5 = t423 - t455 / 0.2e1 + (pkin(5) * t740 - t592 / 0.2e1 + t591 / 0.2e1) * m(7) + t478 + t699;
t447 = t5 * qJD(1) - t19 * qJD(2) + t15 * qJD(3);
t101 = (-t567 / 0.2e1 + t563 / 0.2e1 + t649) * m(7);
t432 = t240 + (t212 * t345 + t344 * t726 + t372 * t65 + t373 * t61) * t691 + t373 * t674;
t441 = t231 * t692 + t581 / 0.2e1 - t580 / 0.2e1;
t25 = (t498 + t677) * t372 + t432 + t441;
t51 = m(7) * (t212 * t372 + t373 * t726) + (t372 ^ 2 + t373 ^ 2) * mrSges(7,3);
t445 = qJD(1) * t25 + qJD(2) * t101 + qJD(3) * t51;
t425 = (-t406 * t65 + (t418 * t618 + t642 * t65) * pkin(4)) * t691 - t630 / 0.2e1 - t629 / 0.2e1 - t624 / 0.2e1 - t623 / 0.2e1 + t406 * t499 - t475 / 0.2e1 + t636 * t730 + t522 * t703 + t716 * t484;
t431 = -t628 / 0.2e1 + t627 / 0.2e1 - t622 / 0.2e1 + t621 / 0.2e1 - t75 * t688 / 0.2e1 + pkin(5) * t498;
t11 = t425 + t431;
t197 = t710 * pkin(4) + (-m(7) * t477 + mrSges(6,1) + mrSges(7,1)) * t636;
t427 = t477 * t692 * t212 + t406 * t508 + t487 * t752 - t749;
t21 = t427 + t702;
t55 = (pkin(5) / 0.2e1 + t484 - t406 / 0.2e1) * t638;
t439 = t11 * qJD(1) - t55 * qJD(2) - t21 * qJD(3) - t197 * qJD(4);
t109 = (t530 + t557) * t691 + m(7) * t672;
t100 = (t563 - t567) * t691 + m(7) * t648;
t89 = t527 + t525;
t40 = -t477 * t638 / 0.2e1 + (mrSges(7,1) + t534) * t344 + t476;
t24 = t284 * t659 + t372 * t498 + t432 - t441;
t20 = t442 / 0.2e1 + t435 + t701;
t18 = -t427 + t480 + t702;
t17 = t381 * t645 + t422 * t456 - t428 + t429 + t701;
t9 = t425 - t431 + t481;
t4 = t423 + t698 + t455 / 0.2e1 + t77 * t534 + (t591 - t592) * t691;
t3 = -t10 * qJD(3) - t12 * qJD(4) + t13 * qJD(5);
t1 = t434 + t436 * t344 + (t492 + (-t384 / 0.4e1 - t713 / 0.4e1 + pkin(3) * mrSges(5,2) / 0.2e1) * t419 + (t652 + t653 + mrSges(5,1) * t690 + (m(6) * t651 + t260 / 0.2e1) * pkin(4)) * t421) * t420 + ((t76 / 0.2e1 + t685) * mrSges(7,3) + t704 * mrSges(6,3) + t452) * t372 + (t574 / 0.4e1 - t338 / 0.4e1 + t709) * t419 + (t661 + t340 / 0.4e1) * t421 + ((t65 / 0.2e1 + t75 / 0.2e1) * mrSges(7,3) + (t86 / 0.2e1 + t90 / 0.2e1) * mrSges(6,3) + t453) * t373 + t412 * t644 + (pkin(4) * t560 + t440) * t693 + t700 + t744;
t22 = [qJD(3) * t6 + qJD(4) * t7 + qJD(5) * t8 + qJD(6) * t36, t3, t1 * qJD(4) + t4 * qJD(5) + t24 * qJD(6) + ((t281 * t97 + t353 * t407 + t714 * t96) * t693 + (t212 * t81 + t231 * t312 + t726 * t77) * t691) * t695 + (t384 * t646 + t382 * t650 + Ifges(4,5) - t637 + (-m(5) * pkin(3) + t573) * t403) * t536 + t463 + (t77 * t614 + t81 * t615 + t96 * t578 + t97 * t579 + mrSges(4,2) * t386 - Ifges(4,6) * t420 + t419 * t341 / 0.2e1 + t353 * t260 + t231 * t259 + (m(5) * t707 - t381 * t420) * pkin(8) + t726 * t469 + t312 * t470 + t714 * t471 + t407 * t472 + t212 * t467 + t281 * t468 + t339 * t646 + t719 * t739 + t718 * t347 / 0.2e1 + (t183 + t181) * t659 - (t187 + t185) * t373 / 0.2e1 + t707 * mrSges(5,3) + (t372 * t736 - t373 * t737 + t465) * t648) * qJD(3), t1 * qJD(3) + (-t475 + (t418 * t91 + t642 * t90) * t689 + m(7) * (t406 * t75 + t636 * t76) - t345 * t577 - t627 + t628 - t621 + t622 - Ifges(5,5) * t554 - Ifges(5,6) * t552 - t474 + t481 + t636 * t703) * qJD(4) + t9 * qJD(5) + t89 * qJD(6) + t462, t4 * qJD(3) + t9 * qJD(4) + (-t623 - t629 - t624 - t630 + (-m(7) * t65 - t582) * pkin(5) + t481) * qJD(5) + t461, t24 * qJD(3) + t89 * qJD(4) + t559; t3, qJD(3) * t41, t17 * qJD(4) + t20 * qJD(5) + t100 * qJD(6) + (-mrSges(4,2) + t735) * t536 + ((t212 * t347 + t312 * t420 + t346 * t726) * t691 + (t281 * t347 + t346 * t714 + t407 * t420) * t693 + m(5) * (t537 * t632 - t631) / 0.2e1) * t695 + t460 + ((t259 + t260 + t573) * t420 + t729 * (t346 * t373 + t347 * t372)) * qJD(3), -t572 + t17 * qJD(3) + (-mrSges(5,1) * t552 + mrSges(5,2) * t554 + t476 + t584) * qJD(4) + t40 * qJD(5) + 0.2e1 * (t527 + (t344 * t528 + t295) * t693) * qJD(4), t571 + t20 * qJD(3) + t40 * qJD(4) + (t118 - t214) * qJD(5), t100 * qJD(3); qJD(4) * t2 + qJD(5) * t5 + qJD(6) * t25 - t463, -qJD(4) * t16 - qJD(5) * t19 + qJD(6) * t101 - t460, qJD(4) * t14 + qJD(5) * t15 + qJD(6) * t51 (m(7) * (-t406 * t212 + t636 * t726) + (-t281 * t642 + t418 * t714) * t689 - t372 * t577 + mrSges(7,3) * t532 + t488 + t706 * pkin(8) + (t532 - t487) * mrSges(6,3) + t750) * qJD(4) + t18 * qJD(5) + t109 * qJD(6) + t448, t18 * qJD(4) + ((-m(7) * t212 - t615) * pkin(5) + t750) * qJD(5) + t447, qJD(4) * t109 + t445; -qJD(3) * t2 + qJD(5) * t11 + qJD(6) * t54 - t462, qJD(3) * t16 - qJD(5) * t55 + t572, -qJD(5) * t21 + qJD(6) * t82 - t448, -t197 * qJD(5) ((-mrSges(6,1) - t513) * t418 - t710) * qJD(5) * pkin(4) + t439, t459; -qJD(3) * t5 - qJD(4) * t11 + qJD(6) * t118 - t461, qJD(3) * t19 + qJD(4) * t55 - t571, qJD(4) * t21 + qJD(6) * t179 - t447, -t439, 0, t458; -t25 * qJD(3) - t54 * qJD(4) - t118 * qJD(5) - t559, -t101 * qJD(3), -qJD(4) * t82 - qJD(5) * t179 - t445, -t459, -t458, 0;];
Cq  = t22;
