% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:31
% EndTime: 2019-03-09 05:34:55
% DurationCPUTime: 14.44s
% Computational Cost: add. (14085->772), mult. (27980->1014), div. (0->0), fcn. (25527->6), ass. (0->392)
t435 = sin(qJ(4));
t738 = pkin(8) - pkin(9);
t387 = t738 * t435;
t437 = cos(qJ(4));
t389 = t738 * t437;
t434 = sin(qJ(6));
t663 = cos(qJ(6));
t227 = -t387 * t663 + t434 * t389;
t478 = t434 * t387 + t389 * t663;
t543 = t663 * t435;
t358 = -t434 * t437 + t543;
t477 = t434 * t435 + t437 * t663;
t573 = -Ifges(7,5) * t477 - Ifges(7,6) * t358;
t748 = -t478 * mrSges(7,1) + t227 * mrSges(7,2) + t573;
t752 = qJD(6) * t748;
t438 = cos(qJ(3));
t657 = pkin(8) * t438;
t436 = sin(qJ(3));
t660 = pkin(3) * t436;
t374 = qJ(2) - t657 + t660;
t440 = -pkin(1) - pkin(7);
t584 = t436 * t440;
t568 = -t437 * t374 + t435 * t584;
t583 = t437 * t438;
t190 = pkin(9) * t583 - t568;
t701 = pkin(4) + pkin(5);
t163 = -t436 * t701 - t190;
t245 = t435 * t374 + t437 * t584;
t427 = t436 * qJ(5);
t216 = t245 + t427;
t587 = t435 * t438;
t412 = pkin(9) * t587;
t173 = t216 + t412;
t77 = t434 * t163 + t173 * t663;
t191 = t245 + t412;
t85 = -t434 * t190 + t191 * t663;
t751 = t77 - t85;
t432 = t435 ^ 2;
t433 = t437 ^ 2;
t566 = t432 + t433;
t567 = t566 * t657;
t361 = -t436 * mrSges(5,2) - mrSges(5,3) * t587;
t563 = mrSges(6,2) * t587;
t627 = t436 * mrSges(6,3);
t366 = -t563 + t627;
t572 = t361 + t366;
t618 = qJ(5) * t435;
t501 = pkin(4) * t437 + t618;
t632 = t358 * mrSges(7,2);
t634 = t477 * mrSges(7,1);
t574 = t634 / 0.2e1 + t632 / 0.2e1;
t710 = -m(7) / 0.2e1;
t712 = -m(6) / 0.2e1;
t368 = -t434 * qJ(5) - t663 * t701;
t369 = qJ(5) * t663 - t434 * t701;
t722 = t358 * t369 - t368 * t477;
t750 = -t501 * t712 - t710 * t722 + t574;
t308 = t434 * t583 - t438 * t543;
t309 = t477 * t438;
t155 = mrSges(7,1) * t309 - mrSges(7,2) * t308;
t428 = t438 * qJ(5);
t408 = t437 * t428;
t553 = t701 * t435;
t215 = t408 + (t440 - t553) * t438;
t267 = Ifges(7,6) * t309;
t268 = Ifges(7,5) * t308;
t503 = t268 + t267;
t672 = -t436 / 0.2e1;
t747 = t215 * t155 - t503 * t672;
t635 = t309 * mrSges(7,3);
t232 = -mrSges(7,1) * t436 - t635;
t746 = -t232 / 0.2e1;
t717 = -t437 * t701 - t618;
t343 = pkin(3) - t717;
t745 = m(7) * t343;
t76 = t163 * t663 - t434 * t173;
t86 = t190 * t663 + t434 * t191;
t744 = t86 + t76;
t735 = mrSges(6,2) + mrSges(5,3);
t741 = (t433 / 0.2e1 + t432 / 0.2e1) * t735;
t617 = qJ(5) * t437;
t354 = -t553 + t617;
t484 = -t434 * mrSges(7,1) - mrSges(7,2) * t663;
t740 = qJD(6) * t484;
t192 = mrSges(7,1) * t358 - mrSges(7,2) * t477;
t700 = t215 / 0.2e1;
t739 = t343 * t155 / 0.2e1 + t192 * t700 + t478 * t746;
t479 = t434 * t227 + t478 * t663;
t348 = Ifges(7,4) * t477;
t198 = t358 * Ifges(7,1) - t348;
t737 = t198 / 0.2e1;
t734 = Ifges(6,4) + Ifges(5,5);
t733 = Ifges(5,6) - Ifges(6,6);
t633 = t477 * mrSges(7,3);
t349 = Ifges(7,4) * t358;
t196 = -Ifges(7,2) * t477 + t349;
t522 = Ifges(7,1) * t477 + t349;
t730 = t196 + t522;
t218 = -pkin(4) * t436 + t568;
t579 = t218 - t568;
t586 = t436 * t437;
t306 = t434 * t586 - t436 * t543;
t307 = t477 * t436;
t729 = t434 * t306 + t663 * t307;
t270 = Ifges(7,4) * t309;
t523 = Ifges(7,1) * t308 + t270;
t588 = t435 * t436;
t360 = -mrSges(5,2) * t438 + mrSges(5,3) * t588;
t367 = mrSges(6,2) * t588 + mrSges(6,3) * t438;
t728 = t360 + t367;
t429 = Ifges(6,5) * t435;
t727 = Ifges(6,1) * t437 + t429;
t430 = Ifges(5,4) * t437;
t726 = -Ifges(5,2) * t435 + t430;
t382 = Ifges(5,1) * t435 + t430;
t725 = -t567 / 0.2e1;
t148 = -t308 * Ifges(7,2) - t436 * Ifges(7,6) + t270;
t724 = t523 + t148;
t721 = t437 * Ifges(6,3) - t429;
t723 = t727 - t721;
t658 = pkin(8) * t436;
t388 = pkin(3) * t438 + t658;
t582 = t438 * t440;
t405 = t435 * t582;
t247 = t388 * t437 - t405;
t248 = t435 * t388 + t437 * t582;
t493 = -t247 * t435 + t248 * t437;
t223 = t428 + t248;
t224 = -pkin(4) * t438 - t247;
t496 = t223 * t437 + t224 * t435;
t505 = Ifges(7,2) * t358 + t348;
t269 = Ifges(7,4) * t308;
t506 = Ifges(7,2) * t309 + t269;
t646 = Ifges(6,5) * t437;
t381 = Ifges(6,1) * t435 - t646;
t720 = t726 + t381 + t382;
t487 = t268 / 0.2e1 + t267 / 0.2e1;
t231 = mrSges(7,2) * t436 - t308 * mrSges(7,3);
t546 = t663 * t231;
t581 = t434 * t746 + t546 / 0.2e1;
t665 = t438 / 0.2e1;
t602 = t307 * t309;
t603 = t306 * t308;
t691 = -t307 / 0.2e1;
t693 = -t306 / 0.2e1;
t718 = (t602 / 0.2e1 + t603 / 0.2e1) * mrSges(7,3) - t231 * t693 - t232 * t691;
t448 = t155 * t665 - t718;
t16 = t448 + t574;
t719 = t16 * qJD(1);
t364 = mrSges(5,1) * t436 - mrSges(5,3) * t583;
t365 = -mrSges(6,1) * t436 + mrSges(6,2) * t583;
t668 = t437 / 0.2e1;
t669 = -t437 / 0.2e1;
t675 = -t435 / 0.2e1;
t715 = t364 * t669 + t365 * t668 + t572 * t675;
t714 = 2 * qJD(3);
t713 = m(5) / 0.2e1;
t711 = m(6) / 0.2e1;
t709 = m(7) / 0.2e1;
t708 = -mrSges(6,1) / 0.2e1;
t707 = -mrSges(7,1) / 0.2e1;
t706 = -mrSges(5,2) / 0.2e1;
t705 = mrSges(6,3) / 0.2e1;
t704 = -mrSges(7,3) / 0.2e1;
t703 = mrSges(7,3) / 0.2e1;
t702 = -Ifges(6,5) / 0.2e1;
t699 = -t227 / 0.2e1;
t230 = -mrSges(7,1) * t438 + mrSges(7,3) * t307;
t698 = t230 / 0.2e1;
t697 = -t231 / 0.2e1;
t695 = -t245 / 0.2e1;
t659 = pkin(4) * t435;
t265 = -t408 + (-t440 + t659) * t438;
t694 = t265 / 0.2e1;
t692 = t306 / 0.2e1;
t690 = -t308 / 0.2e1;
t689 = -t308 / 0.4e1;
t688 = t308 / 0.2e1;
t687 = -t309 / 0.2e1;
t685 = t309 / 0.2e1;
t337 = t501 * t438;
t684 = -t337 / 0.2e1;
t624 = t437 * mrSges(6,1);
t628 = t435 * mrSges(6,3);
t511 = t624 + t628;
t338 = t511 * t438;
t683 = -t338 / 0.2e1;
t510 = mrSges(6,1) * t435 - mrSges(6,3) * t437;
t342 = t438 * t510;
t682 = -t342 / 0.2e1;
t681 = -t477 / 0.2e1;
t680 = -t477 / 0.4e1;
t678 = t358 / 0.2e1;
t562 = mrSges(6,2) * t586;
t622 = t438 * mrSges(6,1);
t363 = -t562 - t622;
t677 = -t363 / 0.2e1;
t676 = t369 / 0.2e1;
t674 = t435 / 0.2e1;
t671 = t436 / 0.2e1;
t670 = t436 / 0.4e1;
t666 = -t438 / 0.2e1;
t664 = t440 / 0.2e1;
t375 = -pkin(3) - t501;
t662 = m(6) * t375;
t661 = m(6) * t438;
t656 = t76 * mrSges(7,2);
t655 = t77 * mrSges(7,1);
t654 = t85 * mrSges(7,1);
t653 = t86 * mrSges(7,2);
t647 = Ifges(5,4) * t435;
t147 = -Ifges(7,4) * t307 + Ifges(7,2) * t306 - Ifges(7,6) * t438;
t149 = -Ifges(7,1) * t307 + Ifges(7,4) * t306 - Ifges(7,5) * t438;
t150 = t309 * Ifges(7,1) - t436 * Ifges(7,5) - t269;
t153 = -mrSges(7,1) * t306 - mrSges(7,2) * t307;
t636 = t309 * mrSges(7,2);
t637 = t308 * mrSges(7,1);
t156 = t636 + t637;
t214 = (-t440 - t354) * t436;
t229 = mrSges(7,2) * t438 + mrSges(7,3) * t306;
t378 = -t617 + t659;
t264 = (-t378 + t440) * t436;
t504 = Ifges(6,3) * t435 + t646;
t294 = Ifges(6,6) * t438 - t436 * t504;
t296 = Ifges(5,6) * t438 - t436 * t726;
t298 = Ifges(6,4) * t438 - t436 * t727;
t509 = Ifges(5,1) * t437 - t647;
t300 = Ifges(5,5) * t438 - t436 * t509;
t340 = t510 * t436;
t512 = mrSges(5,1) * t435 + mrSges(5,2) * t437;
t341 = t436 * t512;
t362 = mrSges(5,1) * t438 + mrSges(5,3) * t586;
t299 = t436 * Ifges(6,4) + t438 * t727;
t474 = t509 * t438;
t301 = t436 * Ifges(5,5) + t474;
t557 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t464 = -t299 / 0.2e1 - t301 / 0.2e1 - t557 * t436;
t488 = Ifges(7,5) * t307 / 0.2e1 + Ifges(7,6) * t693;
t410 = Ifges(6,5) * t583;
t295 = t436 * Ifges(6,6) + Ifges(6,3) * t587 + t410;
t626 = t436 * Ifges(5,6);
t297 = t438 * t726 + t626;
t528 = -t295 / 0.2e1 + t297 / 0.2e1;
t556 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t169 = t405 - t701 * t438 + (pkin(9) * t436 - t388) * t437;
t183 = -pkin(9) * t588 + t223;
t78 = t169 * t663 - t434 * t183;
t79 = t434 * t169 + t183 * t663;
t3 = (t440 * t341 + qJ(2) * mrSges(4,1) + Ifges(7,5) * t687 + Ifges(7,6) * t688 - Ifges(4,4) * t438 + (t298 / 0.2e1 + t300 / 0.2e1 + t557 * t438) * t437 + (-t296 / 0.2e1 + t294 / 0.2e1 + t556 * t438) * t435 + (Ifges(6,2) + Ifges(5,3) + Ifges(4,2) - Ifges(4,1) + Ifges(7,3) + (-m(5) * t440 + t512) * t440) * t436) * t438 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t436 + t464 * t437 + (-t556 * t436 + t528) * t435 + t488) * t436 + m(6) * (t216 * t223 + t218 * t224 + t264 * t265) + m(7) * (t214 * t215 + t76 * t78 + t77 * t79) + t147 * t690 + t150 * t691 + t148 * t692 + t149 * t685 + m(5) * (t245 * t248 - t247 * t568) - t568 * t362 + t247 * t364 + t224 * t365 + t223 * t366 + t216 * t367 + t245 * t360 + t248 * t361 + t218 * t363 - t265 * t340 + t264 * t342 + t76 * t230 + t79 * t231 + t78 * t232 + t77 * t229 + t214 * t156 + t215 * t153;
t638 = t3 * qJD(1);
t631 = t358 * mrSges(7,3);
t246 = t717 * t438;
t625 = t437 * mrSges(5,1);
t629 = t435 * mrSges(5,2);
t513 = t625 - t629;
t339 = t513 * t438;
t409 = Ifges(6,6) * t583;
t500 = -t308 * t76 + t309 * t77;
t532 = -t583 / 0.2e1;
t571 = t364 - t365;
t4 = t409 * t671 + t337 * t342 + t265 * t338 + t506 * t690 + t150 * t688 + t246 * t156 + t86 * t231 + t85 * t232 - t571 * t245 - t572 * t568 + t500 * mrSges(7,3) + m(6) * (-t216 * t568 + t218 * t245 + t265 * t337) + m(7) * (t215 * t246 + t76 * t85 + t77 * t86) + (-t440 * t339 + (t410 / 0.2e1 - t626 / 0.2e1 - t245 * mrSges(5,3) - t216 * mrSges(6,2) + Ifges(5,4) * t532 - t528) * t437 + (-t568 * mrSges(5,3) - t218 * mrSges(6,2) + (t702 + Ifges(5,4) / 0.2e1) * t587 + (Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t583 + t464) * t435) * t438 + t724 * t685 - t747;
t630 = t4 * qJD(1);
t621 = t77 * t358;
t9 = -t77 * t232 + t76 * t231 + (-t77 * mrSges(7,3) - t148 / 0.2e1 - t523 / 0.2e1) * t309 + (t76 * mrSges(7,3) + t506 / 0.2e1 - t150 / 0.2e1) * t308 + t747;
t620 = t9 * qJD(1);
t619 = -t513 - mrSges(4,1);
t552 = t663 * t77;
t589 = t434 * t436;
t27 = -t232 * t589 + (-m(6) * t265 + m(7) * t215 + t156 - t342) * t583 + (t546 + m(7) * (-t434 * t76 + t552) + m(6) * t216 + t366) * t436;
t616 = qJD(1) * t27;
t417 = t438 ^ 2 * t437;
t453 = (t436 ^ 2 * t437 + t417) * t711 + (t436 * t729 + t417) * t709;
t544 = t663 * t477;
t590 = t434 * t358;
t458 = m(6) * t668 + (-t544 + t590) * t710;
t58 = t453 + t458;
t615 = qJD(1) * t58;
t80 = t369 * mrSges(7,1) + mrSges(7,2) * t368;
t614 = qJD(6) * t80;
t612 = t216 * t435;
t609 = t227 * t308;
t608 = t478 * t309;
t606 = t245 * t435;
t495 = -t437 * t568 + t606;
t497 = -t218 * t437 + t612;
t24 = t436 * mrSges(4,1) + t438 * mrSges(4,2) + t358 * t231 - t477 * t232 + mrSges(3,3) + t571 * t437 + t572 * t435 + (m(4) + m(3)) * qJ(2) + m(7) * (-t477 * t76 + t621) + m(6) * t497 + m(5) * t495;
t607 = t24 * qJD(1);
t601 = t309 * t477;
t600 = t343 * t192;
t598 = t358 * t308;
t596 = t368 * t308;
t595 = t369 * t309;
t39 = t307 * mrSges(7,1) - t306 * mrSges(7,2);
t594 = t39 * qJD(6);
t591 = t434 * t309;
t585 = t436 * t438;
t580 = -t216 + t245;
t193 = t632 + t634;
t570 = -t511 - t193;
t547 = t436 * t663;
t516 = mrSges(7,2) * t547;
t569 = t589 * t707 - t516 / 0.2e1;
t565 = qJD(3) * t436;
t564 = mrSges(7,3) * t596;
t559 = t708 - mrSges(5,1) / 0.2e1;
t558 = t705 + t706;
t555 = t85 / 0.2e1 - t77 / 0.2e1;
t554 = -t86 / 0.2e1 - t76 / 0.2e1;
t550 = -t635 / 0.2e1;
t545 = t663 * t308;
t541 = t589 / 0.2e1;
t537 = t587 / 0.2e1;
t531 = t583 / 0.2e1;
t529 = t192 * t665;
t525 = -mrSges(6,2) * qJ(5) - Ifges(5,6);
t517 = pkin(4) * mrSges(6,2) - t734;
t515 = -t544 / 0.2e1;
t514 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t380 = t437 * Ifges(5,2) + t647;
t53 = m(7) * (-t585 + t602 + t603) + 0.4e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * (-0.1e1 + t566) * t585;
t456 = t229 * t691 + t230 * t692 + t231 * t687 + t232 * t688;
t459 = -t306 * t78 + t307 * t79 + t500;
t471 = (t227 * t477 + t358 * t478) * t709;
t476 = m(5) * t493;
t482 = t265 + t496;
t483 = t216 * t437 + t218 * t435 - t264;
t494 = t245 * t437 + t435 * t568;
t7 = t471 + t459 * t710 + (t156 / 0.2e1 + t682 + (-t367 / 0.2e1 - t360 / 0.2e1) * t437 + (t677 + t362 / 0.2e1) * t435 + m(7) * t700 + t482 * t712 - t476 / 0.2e1) * t436 + (-t153 / 0.2e1 - t340 / 0.2e1 + t512 * t672 - t341 / 0.2e1 + t365 * t675 + t364 * t674 + t214 * t710 + t483 * t712 - m(5) * (t494 - 0.2e1 * t584) / 0.2e1 + t572 * t669) * t438 + t456;
t499 = -t7 * qJD(1) + t53 * qJD(2);
t498 = mrSges(7,3) * t515;
t492 = t265 * t378 + t337 * t375;
t491 = -Ifges(5,2) / 0.4e1 - Ifges(6,3) / 0.4e1 + Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1;
t442 = (-t155 / 0.2e1 - t339 / 0.2e1 + t683) * t438 + ((t495 - t497) * t711 - t438 * t741 + t715) * t436 + (t246 * t438 + t306 * t751 + t744 * t307) * t709 + t661 * t684 + t718;
t11 = -t435 * t558 + t437 * t559 + t442 - t750;
t489 = t11 * qJD(1);
t485 = t637 / 0.2e1 + t636 / 0.2e1;
t475 = t438 * t512;
t473 = t591 / 0.2e1 - t545 / 0.2e1;
t472 = m(7) * t479;
t470 = -t475 / 0.2e1;
t441 = (t215 * t354 + t246 * t343 + t744 * t478) * t710 + pkin(3) * t339 / 0.2e1 - t246 * t193 / 0.2e1 + t198 * t689 - t511 * t684 - t354 * t156 / 0.2e1 + t150 * t680 + t375 * t683 + t378 * t682 - t573 * t670 - t730 * t309 / 0.4e1 - t724 * t358 / 0.4e1 + t739 + (t710 * t751 + t697) * t227;
t457 = t78 * t707 + t79 * mrSges(7,2) / 0.2e1 + t488;
t444 = (-pkin(4) * t224 + qJ(5) * t223) * t711 + (t368 * t78 + t369 * t79) * t709 + pkin(4) * t677 + qJ(5) * t367 / 0.2e1 + t223 * t705 + t224 * t708 + t247 * mrSges(5,1) / 0.2e1 + t248 * t706 + t368 * t698 + t229 * t676 + t457;
t1 = (t601 / 0.4e1 + t598 / 0.4e1) * Ifges(7,2) + t441 + (-t301 / 0.4e1 - t299 / 0.4e1 + mrSges(6,3) * t694 + (-0.3e1 / 0.4e1 * Ifges(6,4) - 0.3e1 / 0.4e1 * Ifges(5,5)) * t436 + (t568 / 0.2e1 - t218 / 0.2e1) * mrSges(6,2) + (t579 * t712 - t365 / 0.2e1 + t364 / 0.2e1) * pkin(8) + (mrSges(5,2) * t664 - t429 / 0.4e1 + t380 / 0.4e1 + t721 / 0.4e1 - t491 * t437 + (0.3e1 / 0.4e1 * Ifges(5,4) + t702) * t435) * t438) * t437 + (-t410 / 0.4e1 - t295 / 0.4e1 + t297 / 0.4e1 + t265 * t708 + (0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(6,6)) * t436 + (t695 + t216 / 0.2e1) * mrSges(6,2) + (t366 / 0.2e1 + t361 / 0.2e1 + t580 * t712) * pkin(8) + (t382 / 0.4e1 + t381 / 0.4e1 + t430 / 0.4e1 + mrSges(5,1) * t664 + t491 * t435) * t438) * t435 + t492 * t712 + (pkin(8) * t741 + t514) * t438 + (-t608 / 0.2e1 - t609 / 0.2e1 + t555 * t358 - t554 * t477) * mrSges(7,3) - t348 * t690 + t444;
t14 = -t600 + t505 * t681 + t477 * t737 - t378 * t511 + t380 * t675 - pkin(3) * t512 + t504 * t669 + (m(6) * t378 + t510) * t375 + t730 * t678 + (t509 + t723) * t674 + t720 * t668 + (t745 + t193) * t354;
t447 = t354 * t438 * t710 + t529 + t378 * t661 / 0.2e1;
t451 = (-pkin(4) * t587 + t408) * t711 + (t595 - t596) * t709 + t485;
t21 = t447 + t451;
t469 = -t1 * qJD(1) - t21 * qJD(2) + t14 * qJD(3);
t28 = t600 + (-t196 / 0.2e1 - t522 / 0.2e1) * t358 - (-t505 / 0.2e1 + t737) * t477;
t34 = t529 + t485;
t445 = -(t150 / 0.4e1 - t506 / 0.4e1) * t477 + (-t523 / 0.4e1 - t148 / 0.4e1) * t358 + (t505 / 0.4e1 - t198 / 0.4e1 + mrSges(7,3) * t699) * t308 + (-t196 / 0.4e1 - t522 / 0.4e1 + t478 * t704) * t309 + t231 * t699 - t436 * t573 / 0.4e1 + t739;
t455 = Ifges(7,3) * t665 + t457;
t6 = t445 + t455;
t468 = t6 * qJD(1) + t34 * qJD(2) + t28 * qJD(3);
t460 = t232 * t676 + t368 * t697 - t487;
t12 = (t595 / 0.2e1 - t596 / 0.2e1) * mrSges(7,3) + t554 * mrSges(7,2) + t555 * mrSges(7,1) + t460 + t487;
t32 = (t699 + t227 / 0.2e1) * mrSges(7,2);
t467 = -t12 * qJD(1) + t32 * qJD(3) + t80 * qJD(4);
t105 = (t537 - t473) * m(7);
t443 = (-t265 * t435 + (-t375 * t438 + t658) * t437) * t711 + (t435 * t215 + t343 * t583 + t436 * t479) * t709 + t156 * t674 + t342 * t675 + t193 * t531 - t511 * t532 + t541 * t631 + t436 * t498;
t449 = t224 * t711 + (t434 * t79 + t663 * t78) * t709 + t434 * t229 / 0.2e1 - t622 / 0.2e1 + t663 * t698;
t17 = t443 - t449 + t562;
t75 = (-t570 - t662 + t745) * t435;
t466 = qJD(1) * t17 + qJD(2) * t105 + qJD(3) * t75;
t129 = mrSges(6,3) + m(6) * qJ(5) + m(7) * (-t368 * t434 + t369 * t663) - t484;
t446 = mrSges(7,1) * t541 + t516 / 0.2e1 + t627 + (t245 + 0.2e1 * t427) * t711 + (t369 * t547 + t552 + (-t368 * t436 - t76) * t434) * t709 + t581;
t462 = t434 * t550 + t545 * t703;
t450 = m(6) * t695 + (t434 * t86 + t663 * t85) * t710 + t462;
t20 = t446 + t450;
t463 = t479 * t709;
t49 = t463 - t472 / 0.2e1;
t465 = t20 * qJD(1) - t49 * qJD(3) + t129 * qJD(4);
t36 = t473 * mrSges(7,3) + t569 - t581;
t461 = t36 * qJD(1) + qJD(4) * t484;
t452 = qJD(4) * (-m(6) * t501 - t511 - t513);
t92 = (-t545 + t591) * t709 + (m(6) + t709) * t587;
t59 = m(6) * t586 + m(7) * t729;
t57 = t453 - t458;
t42 = t472 / 0.2e1 + t498 + t463 + (m(6) * pkin(8) + mrSges(6,2)) * t437 + (t515 + t590) * mrSges(7,3);
t37 = t462 + t569 + t581;
t35 = t529 - t485;
t22 = t470 + t682 + (t435 * t559 + t437 * t558) * t438 - t447 + t451;
t19 = t446 - t450 - t563;
t18 = t443 + t449;
t15 = t448 - t574;
t13 = t655 / 0.2e1 + t656 / 0.2e1 + t369 * t550 + t564 / 0.2e1 + t654 / 0.2e1 - t653 / 0.2e1 - t460 + t487;
t10 = t628 / 0.2e1 - t629 / 0.2e1 + t624 / 0.2e1 + t625 / 0.2e1 + t442 + t750;
t8 = (t214 * t438 - t215 * t436 + t459) * t709 + (t436 * t482 + t438 * t483) * t711 + (t494 * t438 + (t493 - 0.2e1 * t582) * t436) * t713 + t365 * t537 - t364 * t587 / 0.2e1 + t471 - t456 + t153 * t665 + t156 * t672 + (t475 + t342) * t671 + (-t341 - t340) * t666 + (t363 / 0.2e1 - t362 / 0.2e1) * t588 + t728 * t586 / 0.2e1 + t572 * t531;
t5 = t445 - t455;
t2 = -t744 * t633 / 0.2e1 + (t474 + t301 + t299) * t437 / 0.4e1 + (-Ifges(6,1) * t587 + t295 + t410) * t435 / 0.4e1 + t725 * mrSges(5,3) + (t504 - t382) * t587 / 0.4e1 + t380 * t532 + t440 * t470 - t441 + (-t435 * t556 - t437 * t557) * t436 + (-t435 * t733 + t437 * t734) * t670 + (t621 + t608) * t703 + t492 * t711 - t85 * t631 / 0.2e1 - t720 * t587 / 0.4e1 + ((t435 * t580 + t437 * t579) * t711 + t715) * pkin(8) - t609 * t704 + t510 * t694 + t506 * t680 + t505 * t689 + t514 * t438 - t435 * t297 / 0.4e1 + t444 + (t723 / 0.4e1 - t721 / 0.4e1) * t583 + (-t612 / 0.2e1 + t606 / 0.2e1 + t725 + t579 * t668) * mrSges(6,2);
t23 = [qJD(2) * t24 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t27 + qJD(6) * t9, t607 + (t306 * t477 + t307 * t358) * m(7) * qJD(2) + t8 * qJD(3) + t10 * qJD(4) + t57 * qJD(5) + t15 * qJD(6), t638 + t8 * qJD(2) + t2 * qJD(4) + t18 * qJD(5) + t5 * qJD(6) + ((t214 * t343 - t227 * t78 + t478 * t79) * t709 + t264 * t662 / 0.2e1) * t714 + (-Ifges(4,5) + (-t381 / 0.2e1 - t382 / 0.2e1) * t437 + (t721 / 0.2e1 + t380 / 0.2e1) * t435 + (-m(5) * pkin(3) + t619) * t440) * t565 + (-mrSges(4,2) * t582 - t78 * t631 + t198 * t691 + t196 * t692 + t149 * t678 + t147 * t681 + t296 * t668 + t294 * t669 + (Ifges(7,5) * t358 - Ifges(7,6) * t477) * t666 - t79 * t633 - Ifges(4,6) * t438 - t375 * t340 - t264 * t511 + pkin(3) * t341 + t343 * t153 - t227 * t230 + t478 * t229 + t214 * t193 + (t728 * t437 + (-t362 + t363) * t435 + t476 + m(6) * t496) * pkin(8) + (t300 + t298) * t674 + (t435 * t734 + t437 * t733) * t665 + t493 * mrSges(5,3) + t496 * mrSges(6,2)) * qJD(3), t10 * qJD(2) + t2 * qJD(3) + t19 * qJD(5) + t13 * qJD(6) + t630 + (-t245 * mrSges(5,1) - t245 * mrSges(6,1) + t568 * mrSges(5,2) - t568 * mrSges(6,3) + mrSges(7,3) * t595 + t409 - t503 - t564 + t653 - t654 + 0.2e1 * (t368 * t85 + t369 * t86) * t709 + 0.2e1 * (-pkin(4) * t245 - qJ(5) * t568) * t711 + (t435 * t517 + t437 * t525) * t438) * qJD(4), qJD(2) * t57 + qJD(3) * t18 + qJD(4) * t19 + qJD(6) * t37 + t616, t620 + t15 * qJD(2) + t5 * qJD(3) + t13 * qJD(4) + t37 * qJD(5) + (-t503 - t655 - t656) * qJD(6); -qJD(3) * t7 + qJD(4) * t11 + qJD(5) * t58 + qJD(6) * t16 - t607, t53 * qJD(3), t22 * qJD(4) + t92 * qJD(5) + t35 * qJD(6) + (t570 + t619) * t565 + ((-t343 * t436 + t608 + t609) * t709 + (t375 * t436 + t567) * t711 + (t567 - t660) * t713) * t714 + t499 + ((t598 - t601) * mrSges(7,3) + (t566 * t735 - mrSges(4,2)) * t438) * qJD(3), t22 * qJD(3) + (m(7) * (t306 * t369 + t307 * t368) - t39) * qJD(4) + t59 * qJD(5) + t594 + t436 * t452 + t489, qJD(3) * t92 + qJD(4) * t59 + t615, t35 * qJD(3) + t39 * qJD(4) - t594 + t719; qJD(2) * t7 - qJD(4) * t1 + qJD(5) * t17 + qJD(6) * t6 - t638, -qJD(4) * t21 + qJD(5) * t105 + qJD(6) * t34 - t499, qJD(4) * t14 + qJD(5) * t75 + qJD(6) * t28, pkin(8) * t452 + t42 * qJD(5) - t752 + t469 + (m(7) * (t227 * t369 + t368 * t478) + t722 * mrSges(7,3) - t517 * t437 + (Ifges(6,6) + t525) * t435 + t748) * qJD(4), qJD(4) * t42 + t466, -qJD(4) * t748 + t468 + t752; -qJD(2) * t11 + qJD(3) * t1 + qJD(5) * t20 - qJD(6) * t12 - t630, t21 * qJD(3) - t489, -qJD(5) * t49 + qJD(6) * t32 - t469, qJD(5) * t129 + t614, t465, t467 - t614; -qJD(2) * t58 - qJD(3) * t17 - qJD(4) * t20 - qJD(6) * t36 - t616, -qJD(3) * t105 - t615, qJD(4) * t49 - t466, -t465 - t740, 0, -t461 + t740; -qJD(2) * t16 - qJD(3) * t6 + qJD(4) * t12 + qJD(5) * t36 - t620, -t34 * qJD(3) - t719, -qJD(4) * t32 - t468, qJD(5) * t484 - t467, t461, 0;];
Cq  = t23;
