% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:28
% EndTime: 2019-03-09 11:47:58
% DurationCPUTime: 16.58s
% Computational Cost: add. (34391->751), mult. (67680->974), div. (0->0), fcn. (76733->8), ass. (0->385)
t405 = sin(qJ(4));
t577 = sin(pkin(10));
t494 = t577 * pkin(2);
t396 = t494 + pkin(8);
t624 = pkin(9) + t396;
t372 = t624 * t405;
t407 = cos(qJ(4));
t373 = t624 * t407;
t404 = sin(qJ(5));
t406 = cos(qJ(5));
t316 = -t372 * t404 + t406 * t373;
t711 = -t404 * t405 + t406 * t407;
t244 = qJ(6) * t711 + t316;
t578 = cos(pkin(10));
t646 = sin(qJ(2));
t647 = cos(qJ(2));
t382 = -t577 * t647 - t578 * t646;
t294 = t711 * t382;
t446 = t404 * t407 + t406 * t405;
t296 = t446 * t382;
t380 = t577 * t646 - t578 * t647;
t286 = t294 * qJ(6);
t401 = -pkin(2) * t647 - pkin(1);
t317 = t380 * pkin(3) + t382 * pkin(8) + t401;
t516 = t646 * pkin(7);
t388 = -qJ(3) * t646 - t516;
t518 = t647 * pkin(7);
t389 = qJ(3) * t647 + t518;
t714 = t577 * t388 + t578 * t389;
t194 = t407 * t317 - t405 * t714;
t554 = t382 * t407;
t159 = pkin(9) * t554 + t194;
t127 = t380 * pkin(4) + t159;
t195 = t405 * t317 + t407 * t714;
t555 = t382 * t405;
t160 = pkin(9) * t555 + t195;
t152 = t404 * t160;
t82 = t406 * t127 - t152;
t60 = t286 + t82;
t52 = pkin(5) * t380 + t60;
t655 = -t446 / 0.2e1;
t659 = t711 / 0.2e1;
t753 = t294 * t655 + t296 * t659;
t532 = mrSges(7,3) * t753;
t595 = t296 * mrSges(7,3);
t214 = -mrSges(7,2) * t380 + t595;
t218 = mrSges(7,1) * t380 + t294 * mrSges(7,3);
t533 = t214 * t659 + t218 * t655;
t576 = qJ(6) * t296;
t154 = t406 * t160;
t83 = t127 * t404 + t154;
t61 = t83 + t576;
t699 = -m(7) / 0.2e1;
t716 = -t406 * t372 - t404 * t373;
t726 = -t446 * qJ(6) + t716;
t768 = (t244 * t296 + t294 * t726 - t446 * t52 + t61 * t711) * t699 - t533 - t532;
t89 = t406 * t159 - t152;
t70 = t286 + t89;
t767 = t52 - t70;
t625 = Ifges(7,4) + Ifges(6,4);
t766 = t625 * t294;
t765 = t625 * t296;
t764 = t625 * t446;
t744 = t625 * t711;
t324 = -mrSges(7,1) * t711 + mrSges(7,2) * t446;
t495 = t578 * pkin(2);
t397 = -t495 - pkin(3);
t635 = t407 * pkin(4);
t387 = t397 - t635;
t338 = -pkin(5) * t711 + t387;
t763 = m(7) * t338 + t324;
t376 = Ifges(7,6) * t446;
t377 = Ifges(6,6) * t446;
t378 = Ifges(7,5) * t711;
t379 = Ifges(6,5) * t711;
t474 = t379 - t377 + t378 - t376;
t762 = -t316 * mrSges(6,1) - t244 * mrSges(7,1) - t716 * mrSges(6,2) - t726 * mrSges(7,2) + t474;
t219 = mrSges(6,1) * t380 + t294 * mrSges(6,3);
t665 = t316 / 0.2e1;
t760 = t244 / 0.2e1;
t761 = t218 * t760 + t219 * t665;
t737 = Ifges(6,1) + Ifges(7,1);
t721 = Ifges(6,5) + Ifges(7,5);
t736 = Ifges(6,2) + Ifges(7,2);
t720 = Ifges(7,6) + Ifges(6,6);
t88 = -t159 * t404 - t154;
t69 = t88 - t576;
t757 = t61 + t69;
t756 = t82 - t89;
t517 = t646 * pkin(2);
t752 = m(4) * t517;
t636 = t406 * pkin(4);
t400 = pkin(5) + t636;
t751 = (t400 - t636) * m(7);
t750 = t83 + t88;
t676 = -t726 / 0.2e1;
t741 = t726 / 0.2e1 + t676;
t740 = -t382 * mrSges(4,1) - t380 * mrSges(4,2);
t739 = t387 * mrSges(6,2) + t338 * mrSges(7,2) + mrSges(7,3) * t726 + t744;
t295 = t446 * t380;
t297 = t711 * t380;
t734 = t295 * t736 - t297 * t625 - t382 * t720;
t733 = t625 * t295 - t297 * t737 - t721 * t382;
t731 = t736 * t711 + t764;
t730 = t736 * t296 - t766;
t729 = t737 * t446 + t744;
t728 = -t737 * t294 + t765;
t691 = mrSges(7,3) / 0.2e1;
t692 = mrSges(6,3) / 0.2e1;
t515 = t691 + t692;
t429 = -t218 / 0.2e1 - t219 / 0.2e1 + t515 * t294;
t623 = t52 - t60;
t723 = m(7) * t623;
t724 = -t723 / 0.2e1 + t429;
t640 = pkin(4) * t404;
t722 = -mrSges(6,1) - mrSges(7,1);
t402 = Ifges(5,4) * t407;
t715 = Ifges(5,1) * t405 + t402;
t318 = -t382 * pkin(3) + t380 * pkin(8) + t517;
t331 = -t578 * t388 + t389 * t577;
t196 = t407 * t318 + t331 * t405;
t197 = t405 * t318 - t331 * t407;
t693 = -mrSges(5,2) / 0.2e1;
t712 = t197 * t693 + t196 * mrSges(5,1) / 0.2e1;
t475 = t294 * t720 + t296 * t721;
t459 = Ifges(5,2) * t405 - t402;
t710 = t405 ^ 2 + t407 ^ 2;
t508 = -t89 / 0.2e1 + t82 / 0.2e1;
t287 = t294 * mrSges(7,1);
t481 = t296 * mrSges(7,2) - t287;
t709 = -t387 * (-mrSges(6,1) * t294 + mrSges(6,2) * t296) / 0.2e1 - t338 * t481 / 0.2e1;
t708 = -t474 * t380 / 0.4e1;
t593 = t297 * mrSges(7,2);
t597 = t295 * mrSges(7,1);
t707 = t593 / 0.2e1 + t597 / 0.2e1;
t596 = t296 * mrSges(6,3);
t215 = -mrSges(6,2) * t380 + t596;
t706 = -mrSges(6,3) * t753 + t215 * t659 + t219 * t655 - t532 + t533;
t705 = 0.2e1 * m(7);
t702 = m(5) / 0.2e1;
t701 = -m(6) / 0.2e1;
t700 = m(6) / 0.2e1;
t698 = m(7) / 0.2e1;
t697 = pkin(5) / 0.2e1;
t696 = m(4) * pkin(2);
t695 = m(7) * pkin(5);
t694 = mrSges(6,1) / 0.2e1;
t690 = t52 / 0.2e1;
t556 = t380 * t407;
t128 = -t382 * pkin(4) + pkin(9) * t556 + t196;
t557 = t380 * t405;
t161 = pkin(9) * t557 + t197;
t84 = t406 * t128 - t161 * t404;
t53 = -pkin(5) * t382 + qJ(6) * t297 + t84;
t689 = t53 / 0.2e1;
t688 = -t60 / 0.2e1;
t687 = -t61 / 0.2e1;
t686 = -t69 / 0.2e1;
t685 = -t83 / 0.2e1;
t684 = m(7) * t61;
t255 = -pkin(4) * t555 + t331;
t164 = -t296 * pkin(5) + t255;
t683 = t164 / 0.2e1;
t167 = -mrSges(7,1) * t296 - mrSges(7,2) * t294;
t682 = -t167 / 0.2e1;
t681 = -t214 / 0.2e1;
t680 = -t215 / 0.2e1;
t216 = -mrSges(7,1) * t382 + t297 * mrSges(7,3);
t679 = t216 / 0.2e1;
t217 = -mrSges(6,1) * t382 + t297 * mrSges(6,3);
t678 = t217 / 0.2e1;
t675 = -t244 / 0.2e1;
t673 = t255 / 0.2e1;
t669 = t295 / 0.2e1;
t668 = t296 / 0.2e1;
t666 = -t297 / 0.2e1;
t664 = -t716 / 0.2e1;
t638 = pkin(5) * t446;
t368 = m(7) * t638;
t663 = -t368 / 0.2e1;
t662 = -t380 / 0.2e1;
t661 = t380 / 0.2e1;
t660 = -t382 / 0.2e1;
t657 = t446 / 0.2e1;
t652 = t397 / 0.2e1;
t651 = -t400 / 0.2e1;
t650 = t405 / 0.2e1;
t649 = -t407 / 0.2e1;
t648 = t407 / 0.2e1;
t645 = m(7) * t164;
t642 = m(7) * t400;
t641 = m(7) * t404;
t639 = pkin(4) * t405;
t637 = t294 * pkin(5);
t634 = t60 * mrSges(7,2);
t633 = t61 * mrSges(7,1);
t632 = t69 * mrSges(7,1);
t631 = t70 * mrSges(7,2);
t630 = t82 * mrSges(6,2);
t629 = t83 * mrSges(6,1);
t628 = t88 * mrSges(6,1);
t627 = t89 * mrSges(6,2);
t626 = mrSges(6,2) + mrSges(7,2);
t622 = mrSges(4,3) * t380;
t621 = mrSges(6,3) * t446;
t620 = mrSges(7,3) * t711;
t619 = mrSges(7,3) * t446;
t617 = Ifges(5,4) * t405;
t608 = Ifges(5,5) * t407;
t606 = Ifges(5,6) * t405;
t605 = pkin(4) * qJD(4);
t604 = t164 * mrSges(7,2);
t599 = t255 * mrSges(6,2);
t598 = t295 * mrSges(6,1);
t594 = t297 * mrSges(6,2);
t254 = -pkin(4) * t557 + t714;
t163 = -t295 * pkin(5) + t254;
t165 = -t593 - t597;
t166 = -t594 - t598;
t168 = -mrSges(6,1) * t296 - mrSges(6,2) * t294;
t212 = mrSges(7,2) * t382 + t295 * mrSges(7,3);
t213 = mrSges(6,2) * t382 + t295 * mrSges(6,3);
t236 = -Ifges(5,6) * t382 + t380 * t459;
t465 = Ifges(5,1) * t407 - t617;
t238 = -Ifges(5,5) * t382 - t380 * t465;
t581 = t407 * mrSges(5,2);
t582 = t405 * mrSges(5,1);
t466 = t581 + t582;
t311 = t466 * t380;
t312 = t466 * t382;
t319 = mrSges(5,2) * t382 + mrSges(5,3) * t557;
t320 = -mrSges(5,2) * t380 + mrSges(5,3) * t555;
t321 = -t382 * mrSges(5,1) + mrSges(5,3) * t556;
t322 = t380 * mrSges(5,1) + mrSges(5,3) * t554;
t438 = -t608 / 0.2e1 + t606 / 0.2e1;
t433 = t438 * t380;
t512 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t513 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t239 = Ifges(5,5) * t380 - t382 * t465;
t537 = t407 * t239;
t237 = Ifges(5,6) * t380 + t382 * t459;
t542 = t405 * t237;
t85 = t404 * t128 + t406 * t161;
t63 = qJ(6) * t295 + t85;
t3 = -pkin(1) * (mrSges(3,1) * t646 + mrSges(3,2) * t647) + (Ifges(4,4) * t380 - t537 / 0.2e1 + t542 / 0.2e1 + mrSges(4,1) * t517 + t433 - t721 * t297 + t720 * t295 + (-Ifges(4,2) + Ifges(4,1) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3)) * t382) * t380 + m(5) * (t194 * t196 + t195 * t197 + t331 * t714) + (t740 + t752) * t401 - t714 * t312 + (-Ifges(3,2) + Ifges(3,1)) * t647 * t646 - t733 * t294 / 0.2e1 + t734 * t668 + t728 * t666 + t730 * t669 + m(6) * (t254 * t255 + t82 * t84 + t83 * t85) + m(7) * (t163 * t164 + t52 * t53 + t61 * t63) + (-t646 ^ 2 + t647 ^ 2) * Ifges(3,4) + (t238 * t649 + t236 * t650 - mrSges(4,2) * t517 + (-Ifges(4,4) - t438) * t382 + t512 * t296 + t513 * t294) * t382 - t331 * t311 + t197 * t320 + t194 * t321 + t196 * t322 + t195 * t319 + t254 * t168 + t255 * t166 + t84 * t219 + t61 * t212 + t83 * t213 + t63 * t214 + t85 * t215 + t52 * t216 + t82 * t217 + t53 * t218 + t164 * t165 + t163 * t167;
t592 = t3 * qJD(1);
t587 = t382 * mrSges(4,3);
t586 = t711 * mrSges(6,3);
t585 = t396 * mrSges(5,3);
t584 = t400 * mrSges(7,3);
t583 = t404 * t63;
t210 = -pkin(4) * t554 - t637;
t413 = (t380 * t513 + t599 + t604 + t765) * t296 - t52 * t595 - t82 * t596 - t164 * t287 + t475 * t661;
t483 = t736 - t737;
t419 = -t255 * mrSges(6,1) + t83 * mrSges(6,3) + t61 * mrSges(7,3) + t296 * t483 - t380 * t512 - t766;
t447 = t194 * t405 - t195 * t407;
t453 = Ifges(5,5) * t405 + Ifges(5,6) * t407;
t467 = -t407 * mrSges(5,1) + t405 * mrSges(5,2);
t390 = Ifges(5,2) * t407 + t617;
t539 = t405 * t390;
t6 = t194 * t320 - t195 * t322 + t88 * t219 + t210 * t167 + t70 * t214 + t89 * t215 + t69 * t218 + m(6) * (t82 * t88 + t83 * t89) + m(7) * (t164 * t210 + t52 * t69 + t61 * t70) + (t237 * t648 + t239 * t650 + t453 * t661 + t331 * t467 + (t539 / 0.2e1 + t715 * t649) * t382 + (-m(6) * t255 - t168) * t635 - t447 * mrSges(5,3)) * t382 + t419 * t294 + t413;
t580 = t6 * qJD(1);
t7 = t60 * t214 + t82 * t215 - t83 * t219 + (-t218 - t723) * t61 + ((-t167 - t645) * pkin(5) + t419) * t294 + t413;
t579 = t7 * qJD(1);
t436 = t382 * t467;
t325 = -mrSges(6,1) * t711 + mrSges(6,2) * t446;
t489 = t325 * t660;
t408 = (-t380 * t396 * t710 - t397 * t382) * t702 + (t295 * t716 - t297 * t316 - t382 * t387) * t700 + (-t244 * t297 + t295 * t726 - t338 * t382) * t698 + t324 * t660 + t489 - t436 / 0.2e1 + (-t380 * t577 + t382 * t578) * t696 / 0.2e1 + t710 * mrSges(5,3) * t662 + (mrSges(7,3) + mrSges(6,3)) * (t295 * t655 - t297 * t659);
t412 = (t196 * t407 + t405 * t197) * t702 + (t446 * t85 + t711 * t84) * t700 + (t446 * t63 + t53 * t711) * t698 + t319 * t650 + t321 * t648 + t752 / 0.2e1;
t14 = -t408 + t412 + (t216 + t217) * t659 + (t212 + t213) * t657 + t740;
t574 = qJD(1) * t14;
t32 = t296 * t214 + t294 * t218 + m(7) * (t294 * t52 + t296 * t61);
t573 = qJD(1) * t32;
t536 = t407 * t320;
t540 = t405 * t322;
t559 = t331 * t382;
t16 = -(t214 + t215) * t297 + (t218 + t219) * t295 + (-t536 + t540 + t622) * t380 + (-t167 - t168 + t312 + t587) * t382 + m(7) * (-t164 * t382 + t295 * t52 - t297 * t61) + m(6) * (-t255 * t382 + t295 * t82 - t297 * t83) + m(5) * (t380 * t447 - t559) + m(4) * (-t380 * t714 - t559);
t572 = t16 * qJD(1);
t569 = t244 * t294;
t568 = t244 * t400;
t567 = t255 * t405;
t566 = t294 * t324;
t565 = t294 * t338;
t563 = t294 * t404;
t561 = t316 * t294;
t552 = t711 * t294;
t550 = t446 * t296;
t549 = t446 * t404;
t548 = t396 * t405;
t547 = t396 * t407;
t546 = t400 * t294;
t545 = t400 * t446;
t543 = t405 * t168;
t541 = t405 * t320;
t535 = t407 * t322;
t534 = t407 * t715;
t101 = (-t382 / 0.4e1 - t552 / 0.4e1 - t550 / 0.4e1) * t705;
t523 = t101 * qJD(1);
t521 = t695 / 0.2e1;
t520 = t296 * t640;
t514 = Ifges(5,1) / 0.4e1 - Ifges(5,2) / 0.4e1;
t511 = -Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1;
t510 = t687 + t686;
t509 = t685 - t88 / 0.2e1;
t504 = -t621 / 0.2e1;
t503 = -t620 / 0.2e1;
t502 = t620 / 0.2e1;
t501 = -t619 / 0.2e1;
t498 = -t596 / 0.2e1;
t497 = -t595 / 0.2e1;
t496 = -t585 / 0.2e1;
t487 = t213 / 0.2e1 + t212 / 0.2e1;
t486 = t681 + t680;
t375 = t446 * mrSges(6,1);
t482 = mrSges(6,2) * t711 + t375;
t374 = t446 * mrSges(7,1);
t480 = mrSges(7,2) * t711 + t374;
t472 = m(7) * t689 + t679;
t471 = -Ifges(5,3) / 0.2e1 + t511;
t470 = t511 * t382;
t469 = t651 + t636 / 0.2e1;
t468 = t163 * t698 - t707;
t454 = -t606 + t608;
t451 = t594 / 0.2e1 + t707;
t269 = t297 * t640;
t414 = (t582 / 0.2e1 + t581 / 0.2e1) * t380 + (t295 * t636 - t269) * t700 + (t295 * t400 - t269) * t698 + t598 / 0.2e1 + t451;
t416 = (-t446 * t756 + t750 * t711) * t701 + (-t446 * t767 + t757 * t711) * t699 + t540 / 0.2e1 - t536 / 0.2e1;
t426 = (t296 * t515 + t486) * t711;
t10 = -t429 * t446 + t414 + t416 + t426;
t450 = t10 * qJD(1);
t423 = (t521 + t694) * t295 + t451;
t12 = -t446 * t724 + t423 + t426;
t449 = t12 * qJD(1);
t21 = t468 + t768;
t91 = m(7) * (t244 * t711 - t446 * t726) + (t446 ^ 2 + t711 ^ 2) * mrSges(7,3);
t448 = -t21 * qJD(1) + t91 * qJD(2);
t348 = t638 + t639;
t362 = t711 * t640;
t148 = (t348 / 0.4e1 + t545 / 0.4e1 - t362 / 0.4e1) * t705 + t480;
t56 = (t210 / 0.4e1 - t520 / 0.4e1 - t546 / 0.4e1) * t705 + t481;
t444 = qJD(1) * t56 + qJD(2) * t148;
t103 = -m(7) * t637 + t481;
t256 = t368 + t480;
t443 = qJD(1) * t103 + qJD(2) * t256;
t442 = t446 * t751;
t441 = Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 - Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1;
t440 = -t626 * t711 - t374 - t375;
t439 = (t362 - t545) * t698;
t437 = t338 * t374 + t387 * t375 - t620 * t726;
t434 = (t404 * t85 + t406 * t84) * t700;
t418 = -t513 * t297 - t512 * t295 + mrSges(7,1) * t689 - t63 * mrSges(7,2) / 0.2e1 + t84 * t694 - t85 * mrSges(6,2) / 0.2e1;
t409 = (t726 * t668 - t569 / 0.2e1) * mrSges(7,3) + (-t561 / 0.2e1 + t716 * t668) * mrSges(6,3) - t374 * t683 - t375 * t673 + t418 + t709;
t411 = (t164 * t348 + t210 * t338 - t244 * t767) * t699 - t210 * t324 / 0.2e1 + t716 * t680 - t331 * t466 / 0.2e1 + t348 * t682 + t542 / 0.4e1 - t537 / 0.4e1 + (t699 * t757 + t681) * t726 + t761;
t417 = -t765 + (-Ifges(6,5) / 0.4e1 - Ifges(7,5) / 0.4e1) * t380 + t441 * t294 - t604 / 0.2e1 - t599 / 0.2e1;
t420 = t766 + (-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t380 + t441 * t296;
t425 = -t316 * t756 + t750 * t716;
t428 = t404 * t487 + t406 * t678;
t1 = t411 + t409 - (mrSges(6,3) * t509 + mrSges(7,3) * t510 + t420) * t446 + pkin(4) * t434 + ((-t70 / 0.2e1 + t690) * mrSges(7,3) + t508 * mrSges(6,3) + t417) * t711 + (-t543 / 0.2e1 + t428) * pkin(4) + ((-t715 / 0.4e1 - t402 / 0.4e1 + t397 * t693 + (t496 - t514) * t405) * t405 + (-0.3e1 / 0.4e1 * t617 - t390 / 0.4e1 + mrSges(5,1) * t652 + (t496 + t514) * t407 + (t387 * t700 + t325 / 0.2e1) * pkin(4)) * t407 + t471) * t382 + t400 * t679 + (pkin(4) * t567 + t425) * t701 + (pkin(4) * t583 + t400 * t53) * t698 + (t535 / 0.2e1 + t541 / 0.2e1) * t396 + (-0.3e1 / 0.4e1 * t608 + 0.3e1 / 0.4e1 * t606 - t379 / 0.4e1 + t377 / 0.4e1 - t378 / 0.4e1 + t376 / 0.4e1) * t380 + t712;
t17 = t465 * t650 - t459 * t648 + t397 * t466 + t534 / 0.2e1 - t539 / 0.2e1 + t739 * t711 - (t483 * t711 + t764) * t446 + t437 + (m(6) * t387 + t325) * t639 + t763 * t348;
t432 = -t1 * qJD(1) + t17 * qJD(2);
t20 = -(-pkin(5) * t763 + t764) * t446 + (-t446 * t483 + t739) * t711 + t437;
t421 = t214 * t676 + t215 * t664 + t761;
t4 = t409 + ((t688 + t690) * mrSges(7,3) + t417) * t711 + (t566 / 0.2e1 + t679) * pkin(5) + t470 + (pkin(5) * t689 + t52 * t760 + t565 * t697 + t60 * t675 + t61 * t741) * m(7) - ((t167 / 0.2e1 + t645 / 0.2e1) * pkin(5) + t420) * t446 + t421 + t708;
t431 = -t4 * qJD(1) + t20 * qJD(2);
t430 = t244 * t636;
t140 = t663 + t442 / 0.2e1;
t268 = t626 * t636 + (-t722 + t751) * t640;
t29 = t741 * mrSges(7,2) + (t664 + t716 / 0.2e1) * mrSges(6,2) + (t675 + t760) * mrSges(7,1) + (-t316 / 0.2e1 + t665) * mrSges(6,1) + (-t568 / 0.2e1 + t430 / 0.2e1 + pkin(5) * t760) * m(7) + (t697 + t469) * t620;
t415 = ((t684 / 0.2e1 + t498 - t486) * t406 + t724 * t404) * pkin(4);
t9 = (t697 + t651) * t595 + (t688 + t70 / 0.2e1) * mrSges(7,2) - t508 * mrSges(6,2) + t510 * mrSges(7,1) + t509 * mrSges(6,1) + (pkin(5) * t686 + t61 * t651) * m(7) + t415;
t424 = t9 * qJD(1) + t29 * qJD(2) - t140 * qJD(3) - t268 * qJD(4);
t410 = t726 * t497 + t716 * t498 + t480 * t683 + t482 * t673 + t61 * t501 + t52 * t503 + t83 * t504 + t561 * t692 + t569 * t691 + t418 - t708 - t709 - (t711 * t737 - t764) * t294 / 0.4e1 + t731 * t294 / 0.4e1 + (t737 * t296 + t766) * t446 / 0.4e1 - (t720 * t380 + t730) * t446 / 0.4e1 + (-t736 * t446 + t729 + t744) * t296 / 0.4e1 + (t736 * t294 + t721 * t380 + t728 + t765) * t711 / 0.4e1;
t199 = t348 * t698 + t439;
t102 = (t550 + t552) * t698 + m(7) * t660;
t94 = -t442 / 0.2e1 + t663 + t440;
t93 = (t520 + t546 + t210) * t698;
t28 = (t430 - t568) * t698 - t244 * t521 + pkin(5) * t503 + t469 * t620 + t762;
t22 = t468 - t768;
t15 = t487 * t446 + (t679 + t678) * t711 + t408 + t412;
t13 = t655 * t723 + t423 + t706;
t11 = t414 - t416 + t706;
t8 = pkin(5) * t497 + t69 * t521 - t634 / 0.2e1 - t633 / 0.2e1 - t630 / 0.2e1 - t629 / 0.2e1 - t631 / 0.2e1 + t632 / 0.2e1 - t627 / 0.2e1 + t628 / 0.2e1 + (t497 - t684 / 0.2e1) * t400 + t415 + t475;
t5 = t410 - t638 * t682 - t621 * t685 - t619 * t687 + t60 * t502 + t470 - t623 * t244 * t698 - t421 + (t472 - t566 / 0.2e1 + (t164 * t446 - t565) * t698) * pkin(5);
t2 = t410 - t411 + t425 * t700 + t88 * t504 + t436 * t652 + t489 * t635 + t69 * t501 + t70 * t502 + t472 * t400 + t380 * t454 / 0.4e1 + t433 - (t541 + t535) * t396 / 0.2e1 - t508 * t586 + (-t465 / 0.4e1 + t390 / 0.2e1) * t554 + (t583 * t698 + t428 + t434 + t543 / 0.2e1 + (-t387 * t554 + t567) * t700) * pkin(4) + (0.2e1 * t715 - t459) * t555 / 0.4e1 + t712 + (t471 + t710 * t585 / 0.2e1) * t382;
t18 = [qJD(2) * t3 + qJD(3) * t16 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t32, t15 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t22 * qJD(6) + t592 + (-Ifges(3,6) * t646 + Ifges(3,5) * t647 - t321 * t548 - mrSges(3,1) * t518 + m(7) * (t163 * t338 + t244 * t63 + t53 * t726) + t726 * t216 + (t446 * t721 + t711 * t720 + t453) * t660 + (m(5) * t396 + mrSges(5,3)) * (-t196 * t405 + t197 * t407) + m(6) * (t254 * t387 + t316 * t85 + t716 * t84) + t716 * t217 + (m(5) * t397 - t578 * t696 - mrSges(4,1) + t467) * t714 + (-t577 * t696 + mrSges(4,2)) * t331 + t534 * t662 + t539 * t661 + t236 * t648 + t238 * t650 - t53 * t619 + t63 * t620 - t84 * t621 + t495 * t622 + t85 * t586 + t494 * t587 + t733 * t657 + t734 * t659 + t729 * t666 + t731 * t669 - t397 * t311 + t387 * t166 - Ifges(4,5) * t380 + Ifges(4,6) * t382 + t338 * t165 + t163 * t324 + t254 * t325 + t316 * t213 + t244 * t212 + mrSges(3,2) * t516 + t319 * t547) * qJD(2), t572 + t15 * qJD(2) + 0.2e1 * (t698 + t700) * (t295 * t711 - t297 * t446) * qJD(3) + t11 * qJD(4) + t13 * qJD(5) + t102 * qJD(6), t580 + t2 * qJD(2) + t11 * qJD(3) + (-t195 * mrSges(5,1) - t194 * mrSges(5,2) + Ifges(5,5) * t555 + Ifges(5,6) * t554 - t296 * t584 + t642 * t69 + t475 - t627 + t628 - t631 + t632) * qJD(4) + t8 * qJD(5) + t93 * qJD(6) + (mrSges(7,3) * t563 + (-t296 * t406 + t563) * mrSges(6,3) + t70 * t641 + m(6) * (t404 * t89 + t406 * t88)) * t605, t579 + t5 * qJD(2) + t13 * qJD(3) + t8 * qJD(4) + (-t629 - t633 - t630 - t634 + (-t595 - t684) * pkin(5) + t475) * qJD(5), qJD(2) * t22 + qJD(3) * t102 + qJD(4) * t93 + t573; -qJD(3) * t14 - qJD(4) * t1 - qJD(5) * t4 - qJD(6) * t21 - t592, qJD(4) * t17 + qJD(5) * t20 + qJD(6) * t91, -t574 (-mrSges(5,1) * t547 + mrSges(5,2) * t548 - t244 * t642 - t584 * t711 + t454 + t762) * qJD(4) + t28 * qJD(5) + t199 * qJD(6) + (-mrSges(7,3) * t549 + (-t406 * t711 - t549) * mrSges(6,3) + m(6) * (-t316 * t406 + t404 * t716) + t726 * t641) * t605 + t432, t28 * qJD(4) + ((-m(7) * t244 - t620) * pkin(5) + t762) * qJD(5) + t431, t199 * qJD(4) + t448; qJD(2) * t14 - qJD(4) * t10 - qJD(5) * t12 - qJD(6) * t101 - t572, t574, 0 (-t466 - t480 - t482) * qJD(4) + t94 * qJD(5) + 0.2e1 * ((-t446 * t636 + t362) * t700 + t439) * qJD(4) - t450, t94 * qJD(4) + (-t368 + t440) * qJD(5) - t449, -t523; qJD(2) * t1 + qJD(3) * t10 + qJD(5) * t9 - qJD(6) * t56 - t580, qJD(5) * t29 - qJD(6) * t148 - t432, -qJD(5) * t140 + t450, -t268 * qJD(5) (-t626 * t406 + (-t695 + t722) * t404) * qJD(5) * pkin(4) + t424, -t444; qJD(2) * t4 + qJD(3) * t12 - qJD(4) * t9 - qJD(6) * t103 - t579, -qJD(4) * t29 - qJD(6) * t256 - t431, qJD(4) * t140 + t449, -t424, 0, -t443; qJD(2) * t21 + qJD(3) * t101 + qJD(4) * t56 + qJD(5) * t103 - t573, t148 * qJD(4) + t256 * qJD(5) - t448, t523, t444, t443, 0;];
Cq  = t18;
