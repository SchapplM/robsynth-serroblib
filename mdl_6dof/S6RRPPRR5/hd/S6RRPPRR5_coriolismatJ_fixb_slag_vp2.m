% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:23
% EndTime: 2019-03-09 09:07:43
% DurationCPUTime: 11.42s
% Computational Cost: add. (17152->831), mult. (39690->1132), div. (0->0), fcn. (38641->8), ass. (0->419)
t466 = sin(qJ(6));
t459 = t466 ^ 2;
t469 = cos(qJ(6));
t461 = t469 ^ 2;
t737 = (t461 / 0.2e1 + t459 / 0.2e1) * mrSges(7,3);
t736 = Ifges(7,3) / 0.2e1;
t463 = sin(pkin(6));
t471 = cos(qJ(2));
t611 = t463 * t471;
t468 = sin(qJ(2));
t612 = t463 * t468;
t319 = -t463 * pkin(1) - pkin(2) * t611 - qJ(3) * t612;
t267 = pkin(3) * t611 - t319;
t188 = (pkin(4) * t471 - pkin(9) * t468) * t463 + t267;
t429 = qJ(4) * t611;
t464 = cos(pkin(6));
t450 = t464 * qJ(3);
t578 = pkin(8) * t611;
t673 = pkin(1) * t468;
t241 = t578 - t429 + t450 + (-pkin(9) + t673) * t464;
t467 = sin(qJ(5));
t470 = cos(qJ(5));
t97 = t188 * t470 - t241 * t467;
t85 = -pkin(5) * t611 - t97;
t735 = -t85 / 0.2e1;
t734 = m(7) * t85;
t594 = t467 * t471;
t572 = t463 * t594;
t733 = t572 * t736;
t732 = Ifges(4,5) + Ifges(5,4);
t365 = -t464 * t467 + t470 * t612;
t638 = t365 * mrSges(6,3);
t253 = -t365 * t466 + t469 * t611;
t254 = t365 * t469 + t466 * t611;
t138 = -mrSges(7,1) * t253 + mrSges(7,2) * t254;
t276 = mrSges(6,1) * t611 - t638;
t731 = -t138 + t276;
t457 = Ifges(7,4) * t469;
t411 = Ifges(7,1) * t466 + t457;
t364 = t464 * t470 + t467 * t612;
t640 = t364 * mrSges(7,2);
t649 = t253 * mrSges(7,3);
t169 = -t640 + t649;
t641 = t364 * mrSges(7,1);
t648 = t254 * mrSges(7,3);
t170 = t641 - t648;
t681 = -t469 / 0.2e1;
t687 = -t466 / 0.2e1;
t510 = t169 * t687 + t170 * t681;
t456 = Ifges(7,5) * t469;
t570 = -t456 / 0.2e1;
t653 = Ifges(7,6) * t466;
t729 = t653 / 0.2e1 + t570;
t669 = pkin(10) * t470;
t670 = pkin(5) * t467;
t415 = t669 - t670;
t465 = qJ(3) - pkin(9);
t602 = t466 * t467;
t300 = t415 * t469 + t465 * t602;
t596 = t467 * t469;
t301 = t415 * t466 - t465 * t596;
t728 = -t300 * t466 + t301 * t469;
t530 = t653 - t456;
t531 = Ifges(7,2) * t466 - t457;
t584 = t470 * t471;
t313 = (t466 * t584 + t468 * t469) * t463;
t571 = t463 * t584;
t314 = -t466 * t612 + t469 * t571;
t727 = Ifges(7,6) * t313 / 0.2e1 - Ifges(7,5) * t314 / 0.2e1;
t679 = t469 / 0.2e1;
t726 = mrSges(7,3) * (t253 * t687 + t254 * t679) - t510;
t430 = qJ(3) * t611;
t472 = -pkin(2) - pkin(3);
t298 = t472 * t612 + t430;
t458 = pkin(4) - t472;
t614 = t458 * t468;
t215 = t430 + (-pkin(9) * t471 - t614) * t463;
t369 = t464 * t471 * pkin(1) - pkin(8) * t612;
t427 = qJ(4) * t612;
t297 = t427 + t369;
t119 = t215 * t470 - t297 * t467;
t105 = pkin(5) * t612 - t119;
t120 = t467 * t215 + t470 * t297;
t692 = -t411 / 0.4e1;
t659 = Ifges(7,4) * t466;
t409 = Ifges(7,2) * t469 + t659;
t693 = t409 / 0.4e1;
t405 = -mrSges(7,1) * t469 + mrSges(7,2) * t466;
t694 = -t405 / 0.2e1;
t642 = t314 * mrSges(7,2);
t643 = t313 * mrSges(7,1);
t175 = t642 + t643;
t712 = t175 / 0.2e1;
t719 = m(7) / 0.2e1;
t725 = t105 * t694 - t119 * mrSges(6,1) / 0.2e1 + t120 * mrSges(6,2) / 0.2e1 + t313 * t693 + t314 * t692 + (t105 * t719 + t712) * pkin(5);
t724 = -m(7) * pkin(5) - mrSges(6,1) + t405;
t723 = m(4) / 0.2e1;
t722 = m(5) / 0.2e1;
t721 = m(6) / 0.2e1;
t720 = -m(7) / 0.2e1;
t718 = -mrSges(7,1) / 0.2e1;
t106 = -pkin(10) * t612 + t120;
t404 = pkin(5) * t572;
t579 = t464 * t673;
t187 = -t579 + t404 + t429 + (-pkin(8) - t669) * t611;
t75 = -t106 * t466 + t187 * t469;
t717 = t75 / 0.2e1;
t76 = t106 * t469 + t187 * t466;
t716 = -t76 / 0.2e1;
t245 = Ifges(7,4) * t253;
t101 = Ifges(7,1) * t254 + Ifges(7,5) * t364 + t245;
t715 = t101 / 0.2e1;
t647 = t254 * Ifges(7,4);
t141 = Ifges(7,1) * t253 - t647;
t714 = -t141 / 0.4e1;
t713 = -t170 / 0.2e1;
t323 = -t464 * pkin(2) - t369;
t242 = -t464 * pkin(3) + t323 - t427;
t216 = t464 * pkin(4) - t242;
t711 = -t216 / 0.2e1;
t710 = t253 / 0.2e1;
t709 = t254 / 0.2e1;
t396 = pkin(5) * t470 + pkin(10) * t467 + t458;
t601 = t466 * t470;
t281 = t396 * t469 - t465 * t601;
t708 = t281 / 0.2e1;
t706 = t301 / 0.2e1;
t705 = -t313 / 0.2e1;
t445 = Ifges(7,4) * t602;
t656 = Ifges(7,5) * t470;
t338 = -Ifges(7,1) * t596 + t445 + t656;
t704 = t338 / 0.2e1;
t703 = t364 / 0.2e1;
t702 = t364 / 0.4e1;
t701 = t365 / 0.2e1;
t371 = t467 * t405;
t700 = t371 / 0.2e1;
t630 = t469 * mrSges(7,2);
t633 = t466 * mrSges(7,1);
t407 = t630 + t633;
t372 = t407 * t467;
t699 = -t372 / 0.2e1;
t626 = t470 * mrSges(7,2);
t392 = mrSges(7,3) * t602 - t626;
t698 = -t392 / 0.2e1;
t697 = t392 / 0.2e1;
t628 = t470 * mrSges(7,1);
t394 = mrSges(7,3) * t596 + t628;
t696 = -t394 / 0.2e1;
t695 = t394 / 0.2e1;
t691 = -t456 / 0.4e1;
t688 = t465 / 0.2e1;
t686 = t466 / 0.2e1;
t685 = t466 / 0.4e1;
t684 = -t467 / 0.2e1;
t683 = t467 / 0.2e1;
t682 = t468 / 0.2e1;
t680 = -t469 / 0.4e1;
t678 = t469 / 0.4e1;
t677 = -t470 / 0.2e1;
t676 = t470 / 0.2e1;
t675 = t470 / 0.4e1;
t582 = t459 + t461;
t595 = t467 * t470;
t289 = (-0.1e1 + t582) * t595;
t674 = m(7) * t289;
t671 = pkin(5) * t365;
t109 = pkin(5) * t364 - pkin(10) * t365 + t216;
t98 = t188 * t467 + t241 * t470;
t86 = pkin(10) * t611 + t98;
t56 = t109 * t469 - t466 * t86;
t668 = t56 * mrSges(7,3);
t57 = t109 * t466 + t469 * t86;
t667 = t57 * mrSges(7,3);
t666 = mrSges(5,2) + mrSges(4,3);
t665 = -Ifges(3,6) - Ifges(5,6);
t664 = m(7) * qJD(4);
t663 = Ifges(6,1) * t470;
t661 = Ifges(6,4) * t467;
t660 = Ifges(6,4) * t470;
t658 = Ifges(6,5) * t470;
t455 = Ifges(6,6) * t467;
t652 = Ifges(7,6) * t470;
t651 = Ifges(7,3) * t365;
t650 = Ifges(7,3) * t467;
t646 = t281 * mrSges(7,3);
t587 = t469 * t470;
t282 = t396 * t466 + t465 * t587;
t645 = t282 * mrSges(7,3);
t100 = t253 * Ifges(7,2) + t364 * Ifges(7,6) + t647;
t149 = Ifges(7,4) * t314 - Ifges(7,2) * t313 + Ifges(7,6) * t572;
t150 = Ifges(7,1) * t314 - Ifges(7,4) * t313 + Ifges(7,5) * t572;
t352 = Ifges(6,4) * t364;
t574 = Ifges(6,5) * t611;
t637 = t365 * Ifges(6,1);
t183 = -t352 + t574 + t637;
t230 = -mrSges(7,2) * t572 - mrSges(7,3) * t313;
t231 = mrSges(7,1) * t572 - mrSges(7,3) * t314;
t534 = -t661 + t663;
t256 = (-Ifges(6,5) * t468 + t471 * t534) * t463;
t370 = t578 + t579;
t318 = t370 + t450;
t257 = t318 - t429;
t639 = t364 * mrSges(6,3);
t275 = -mrSges(6,2) * t611 - t639;
t299 = t370 - t429;
t627 = t470 * mrSges(6,2);
t631 = t467 * mrSges(6,1);
t535 = t627 + t631;
t317 = t535 * t611;
t345 = (mrSges(6,2) * t468 - mrSges(6,3) * t594) * t463;
t346 = (-mrSges(6,1) * t468 - mrSges(6,3) * t584) * t463;
t366 = (-mrSges(4,1) * t471 - mrSges(4,3) * t468) * t463;
t367 = (mrSges(5,1) * t471 + mrSges(5,2) * t468) * t463;
t368 = pkin(2) * t612 - t430;
t389 = -mrSges(5,1) * t464 - mrSges(5,3) * t612;
t575 = mrSges(5,3) * t611;
t390 = t464 * mrSges(5,2) - t575;
t432 = mrSges(4,2) * t611;
t391 = t464 * mrSges(4,3) + t432;
t431 = mrSges(5,2) * t611;
t433 = Ifges(4,6) * t612;
t435 = Ifges(3,5) * t611;
t437 = Ifges(4,4) * t611;
t536 = mrSges(6,1) * t364 + mrSges(6,2) * t365;
t532 = -Ifges(6,2) * t467 + t660;
t548 = -(-Ifges(6,6) * t468 + t471 * t532) * t463 / 0.2e1 + t733 - t727;
t636 = t365 * Ifges(6,4);
t182 = -t364 * Ifges(6,2) + Ifges(6,6) * t611 + t636;
t99 = Ifges(7,5) * t254 + Ifges(7,6) * t253 + t364 * Ifges(7,3);
t561 = t99 / 0.2e1 - t182 / 0.2e1;
t635 = t369 * mrSges(3,2);
t3 = t297 * t390 + t369 * t391 + t298 * t367 + t368 * t366 + t98 * t345 + t97 * t346 + t216 * t317 + t120 * t275 + t119 * t276 + t57 * t230 + t56 * t231 + t76 * t169 + t75 * t170 + t85 * t175 + t105 * t138 + m(5) * (t242 * t299 + t257 * t297 + t267 * t298) + m(6) * (t119 * t97 + t120 * t98 - t216 * t299) + m(7) * (t105 * t85 + t56 * t75 + t57 * t76) + m(4) * (t318 * t369 + t319 * t368 + t323 * t370) + (t437 / 0.2e1 + t433 / 0.2e1 - t635 + t435 / 0.2e1 + (-mrSges(4,1) - mrSges(3,1)) * t370) * t464 + t267 * t431 + t548 * t364 + t256 * t701 + t100 * t705 + t150 * t709 + t149 * t710 + t314 * t715 + (t389 - t536) * t299 + ((Ifges(3,4) * t611 - t319 * mrSges(4,3) + t323 * mrSges(4,2) - t242 * mrSges(5,3) + t183 * t676 + t561 * t467 + (-pkin(1) * mrSges(3,2) + (t658 / 0.2e1 - t455 / 0.2e1 - t732) * t471) * t463 + (Ifges(3,5) / 0.2e1 - Ifges(5,5) + Ifges(4,4) / 0.2e1) * t464) * t471 + (-Ifges(6,5) * t365 / 0.2e1 + Ifges(6,6) * t703 + t319 * mrSges(4,1) + t257 * mrSges(5,3) - t267 * mrSges(5,1) + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t468) * t463 + (-t318 + t370) * mrSges(4,2) + (Ifges(4,6) / 0.2e1 + t665) * t464 + (Ifges(3,1) + Ifges(4,1) + Ifges(5,1) - Ifges(3,2) - Ifges(5,2) - Ifges(4,3) - Ifges(6,3)) * t611 + t732 * t612) * t468) * t463;
t644 = t3 * qJD(1);
t145 = t364 * t530 + t651;
t146 = Ifges(7,6) * t365 + t364 * t531;
t533 = Ifges(7,1) * t469 - t659;
t147 = Ifges(7,5) * t365 - t364 * t533;
t190 = t407 * t364;
t632 = t466 * mrSges(7,3);
t204 = -mrSges(7,2) * t365 + t364 * t632;
t629 = t469 * mrSges(7,3);
t205 = mrSges(7,1) * t365 + t364 * t629;
t217 = t365 * mrSges(6,1) - t364 * mrSges(6,2);
t218 = -Ifges(6,2) * t365 - t352;
t219 = -Ifges(6,1) * t364 - t636;
t558 = t611 / 0.2e1;
t583 = -Ifges(6,5) * t364 - Ifges(6,6) * t365;
t593 = t469 * t101;
t608 = t466 * t100;
t220 = pkin(10) * t364 + t671;
t73 = t220 * t469 - t466 * t97;
t74 = t220 * t466 + t469 * t97;
t4 = -t85 * t190 + t74 * t169 + t57 * t204 + t146 * t710 + t147 * t709 + t73 * t170 + t56 * t205 + m(7) * (t56 * t73 + t57 * t74) + t216 * t217 + t97 * t275 + t583 * t558 + (t219 / 0.2e1 + t561) * t365 + (t145 / 0.2e1 - t183 / 0.2e1 - t218 / 0.2e1 + t608 / 0.2e1 - t593 / 0.2e1 + t97 * mrSges(6,3)) * t364 + (-t638 - t731 + t734) * t98;
t634 = t4 * qJD(1);
t137 = mrSges(7,1) * t254 + mrSges(7,2) * t253;
t139 = Ifges(7,5) * t253 - Ifges(7,6) * t254;
t140 = -Ifges(7,2) * t254 + t245;
t7 = t85 * t137 + t56 * t169 + t139 * t703 - t57 * t170 + (-t667 - t100 / 0.2e1 + t141 / 0.2e1) * t254 + (t140 / 0.2e1 + t715 - t668) * t253;
t625 = t7 * qJD(1);
t624 = t75 * t466;
t529 = -t466 * t56 + t469 * t57;
t592 = t469 * t169;
t607 = t466 * t170;
t14 = (t390 + t391) * t464 + (-t366 + t367) * t612 + t731 * t365 + (t275 + t592 - t607) * t364 + m(7) * (t364 * t529 - t365 * t85) + m(6) * (t364 * t98 + t365 * t97) + m(5) * (t257 * t464 + t267 * t612) + m(4) * (t318 * t464 - t319 * t612);
t623 = qJD(1) * t14;
t406 = mrSges(6,1) * t470 - mrSges(6,2) * t467;
t460 = t467 ^ 2;
t462 = t470 ^ 2;
t581 = -t460 - t462;
t609 = t465 * t471;
t477 = -t298 * t722 + (-t460 * t463 * t609 + t281 * t313 - t282 * t314) * t719 + t313 * t695 + t314 * t698 + (t581 * t609 + t614) * t463 * t721;
t527 = t76 * t469 - t624;
t479 = -m(5) * t298 / 0.2e1 - m(6) * (t119 * t470 + t120 * t467) / 0.2e1 + (-t105 * t470 + t467 * t527) * t720;
t598 = t467 * t372;
t490 = (t598 / 0.2e1 + (t460 / 0.2e1 + t462 / 0.2e1) * mrSges(6,3)) * t471;
t590 = t469 * t230;
t605 = t466 * t231;
t16 = -t431 + (t712 - t346 / 0.2e1) * t470 + (-t590 / 0.2e1 + t605 / 0.2e1 - t345 / 0.2e1) * t467 + ((t406 / 0.2e1 + mrSges(5,1)) * t468 + t490) * t463 + t477 + t479;
t622 = qJD(1) * t16;
t586 = t470 * t275;
t600 = t467 * t138;
t15 = t314 * t169 - t313 * t170 - m(7) * (t313 * t56 - t314 * t57) - t276 * t572 - t536 * t612 + (t468 * t389 + (t390 + t586 + t600) * t471 + t594 * t734 - m(6) * (t216 * t468 - t584 * t98 + t594 * t97) - m(5) * (-t242 * t468 - t257 * t471)) * t463;
t621 = t15 * qJD(1);
t618 = t313 * t466;
t617 = t314 * t469;
t616 = t364 * t470;
t615 = t365 * t467;
t613 = t462 * t465;
t610 = t465 * t467;
t606 = t466 * t205;
t604 = t466 * t394;
t395 = -mrSges(7,1) * t467 + mrSges(7,3) * t587;
t603 = t466 * t395;
t599 = t467 * t364;
t408 = Ifges(7,5) * t466 + Ifges(7,6) * t469;
t597 = t467 * t408;
t591 = t469 * t204;
t589 = t469 * t392;
t393 = mrSges(7,2) * t467 + mrSges(7,3) * t601;
t588 = t469 * t393;
t585 = t470 * t365;
t374 = Ifges(7,5) * t602 + Ifges(7,6) * t596;
t580 = mrSges(6,1) + t694;
t373 = t470 * t407;
t447 = t460 * t465;
t505 = m(7) * t728;
t508 = t604 / 0.2e1 - t589 / 0.2e1;
t522 = -t281 * t466 + t282 * t469;
t41 = (t447 - t613) * t719 + (t373 / 0.2e1 + t522 * t719 - t508) * t470 + (t699 + t505 / 0.2e1 + t588 / 0.2e1 - t603 / 0.2e1) * t467;
t573 = t674 / 0.2e1;
t63 = t371 * t676 + (t392 * t686 + t394 * t679) * t467 - t460 * t737;
t577 = qJD(3) * t573 + t41 * qJD(5) - t63 * qJD(6);
t576 = t289 * t664;
t568 = -t648 / 0.2e1;
t565 = t632 / 0.2e1;
t564 = t629 / 0.2e1;
t563 = t627 / 0.2e1;
t562 = t694 + mrSges(6,1) / 0.2e1;
t560 = -t612 / 0.2e1;
t559 = -t611 / 0.2e1;
t555 = -t602 / 0.2e1;
t554 = -t601 / 0.2e1;
t553 = t596 / 0.2e1;
t551 = t590 / 0.2e1;
t550 = t587 / 0.2e1;
t412 = -Ifges(6,1) * t467 - t660;
t549 = t412 * t676;
t336 = t467 * t531 + t652;
t376 = t411 * t467;
t547 = -t336 / 0.4e1 + t376 / 0.4e1;
t375 = Ifges(7,2) * t596 + t445;
t546 = t338 / 0.4e1 + t375 / 0.4e1;
t335 = Ifges(7,3) * t470 + t467 * t530;
t410 = -Ifges(6,2) * t470 - t661;
t545 = -t410 / 0.4e1 + t335 / 0.4e1;
t541 = t582 * t364;
t540 = t582 * t470;
t538 = t467 * t559;
t528 = -t466 * t73 + t469 * t74;
t526 = t600 / 0.2e1 + t170 * t554 + t169 * t550 + t276 * t684 + t586 / 0.2e1;
t480 = t137 * t677 - t467 * t726;
t515 = -t643 / 0.2e1 - t642 / 0.2e1;
t22 = t480 + t515;
t525 = -qJD(1) * t22 + qJD(2) * t63;
t23 = (-t640 / 0.2e1 + t169 / 0.2e1 - t649 / 0.2e1) * t469 + (-t641 / 0.2e1 + t568 + t713) * t466;
t78 = (-t626 / 0.2e1 + t697) * t469 + (-t628 / 0.2e1 + t696) * t466;
t524 = qJD(1) * t23 + qJD(2) * t78;
t487 = (-t466 * t74 - t469 * t73) * t719 + t204 * t687 + t205 * t681;
t495 = m(7) * (pkin(10) * t541 + t671);
t518 = mrSges(6,2) - t737;
t27 = -t580 * t365 + t518 * t364 - t495 / 0.2e1 + t487;
t486 = (-t300 * t469 - t301 * t466) * t719 + t393 * t687 + t395 * t681;
t494 = m(7) * (pkin(10) * t540 - t670);
t59 = t580 * t467 + t518 * t470 - t494 / 0.2e1 + t486;
t523 = qJD(1) * t27 + qJD(2) * t59;
t521 = -t617 - t618;
t520 = t529 - t98;
t519 = t528 + t85;
t517 = -pkin(5) * t371 / 0.2e1 + t456 * t675;
t516 = mrSges(7,2) * t706 + t300 * t718;
t514 = -t630 / 0.2e1 - t633 / 0.2e1;
t513 = pkin(10) * t698 + t547;
t512 = pkin(10) * t696 + t546;
t511 = t608 / 0.4e1 - t593 / 0.4e1;
t509 = t607 / 0.2e1 - t592 / 0.2e1;
t507 = t409 * t686 + t411 * t681;
t504 = m(7) * (-t313 * t469 + t314 * t466);
t484 = (t585 + t599) * t721 + (t467 * t541 + t585) * t719;
t44 = (m(5) + t721) * t612 - t504 / 0.2e1 + t484;
t503 = t44 * qJD(1) + qJD(2) * t573;
t502 = (-Ifges(6,2) / 0.4e1 - Ifges(7,3) / 0.4e1) * t467 - t412 / 0.4e1;
t501 = t514 * t470;
t500 = mrSges(6,1) * t711 - t219 / 0.4e1 + t182 / 0.4e1 - t99 / 0.4e1;
t337 = -Ifges(7,6) * t467 + t470 * t531;
t339 = -Ifges(7,5) * t467 - t470 * t533;
t474 = (t276 * t677 + t138 * t676 + (t467 * t98 + t470 * t85) * t719 + t275 * t684 - t190 * t683 + (-t585 / 0.2e1 - t599 / 0.2e1) * mrSges(6,3)) * t465 + (t281 * t73 + t282 * t74 + t300 * t56 + t301 * t57) * t719 + t253 * t337 / 0.4e1 + t254 * t339 / 0.4e1 + t205 * t708 + t282 * t204 / 0.2e1 + t300 * t170 / 0.2e1 + t169 * t706 + t458 * t217 / 0.2e1 + t56 * t395 / 0.2e1 + t57 * t393 / 0.2e1 + t73 * t695 + t74 * t697 + t373 * t735 + t98 * t699;
t496 = -t218 / 0.4e1 + mrSges(6,2) * t711 - t183 / 0.4e1 + t145 / 0.4e1;
t2 = t474 + (-t663 / 0.4e1 + t661 / 0.4e1 + t545) * t365 + (Ifges(6,3) * t682 + (-0.3e1 / 0.4e1 * t658 + t455 / 0.4e1 + (Ifges(6,6) / 0.2e1 - t408 / 0.4e1) * t467) * t471) * t463 + (-t149 / 0.4e1 - t470 * t101 / 0.4e1 - t467 * t147 / 0.4e1 + mrSges(7,3) * t716 + (-t656 / 0.4e1 - t338 / 0.4e1) * t364 + (m(7) * t716 - t230 / 0.2e1) * pkin(10)) * t469 + (-t150 / 0.4e1 + t100 * t675 + t467 * t146 / 0.4e1 + mrSges(7,3) * t717 + (t652 / 0.4e1 + t336 / 0.4e1) * t364 + (m(7) * t717 + t231 / 0.2e1) * pkin(10)) * t466 + t500 * t467 + t496 * t470 + (t660 / 0.4e1 + t502) * t364 + t725;
t28 = t338 * t550 + t339 * t553 + t336 * t554 + t337 * t555 + t470 * t465 * t372 + t373 * t610 + t335 * t683 - m(7) * (t465 ^ 2 * t595 + t281 * t300 + t282 * t301) - t301 * t392 - t282 * t393 - t300 * t394 - t281 * t395 + t549 + t458 * t535 + (t534 + t410) * t684 + (t470 * t530 + t532 - t650) * t677;
t499 = t2 * qJD(1) - t28 * qJD(2) + t41 * qJD(4);
t32 = t374 * t676 + t281 * t392 - t282 * t394 + (t465 * t371 + (-t376 / 0.2e1 + t336 / 0.2e1 + t645) * t469 + (t704 + t375 / 0.2e1 - t646) * t466) * t467;
t476 = (-t646 / 0.2e1 + t546) * t253 + (-t645 / 0.2e1 + t547) * t254 + t169 * t708 + t282 * t713 + t374 * t702 + t139 * t675 + t56 * t697 + t57 * t696 + t85 * t700;
t481 = t137 * t688 + (t714 + t100 / 0.4e1 + t667 / 0.2e1) * t469 + (t140 / 0.4e1 + t101 / 0.4e1 - t668 / 0.2e1) * t466;
t488 = t75 * t718 + t76 * mrSges(7,2) / 0.2e1 + t727;
t6 = (Ifges(7,3) * t559 + t481) * t467 + t476 + t488;
t498 = t6 * qJD(1) + t32 * qJD(2) - t63 * qJD(4);
t497 = t370 + 0.2e1 * t450;
t482 = (-t618 / 0.2e1 - t617 / 0.2e1) * mrSges(7,3) + (pkin(10) * t521 + t404) * t719;
t12 = (mrSges(6,2) * t558 - t190 / 0.2e1 + t520 * t720 - t275 / 0.2e1 - t639 / 0.2e1 + t509) * t470 + (-t138 / 0.2e1 - t591 / 0.2e1 + t519 * t720 + t606 / 0.2e1 + t638 / 0.2e1 + t276 / 0.2e1 + t562 * t611) * t467 + t482;
t493 = -t12 * qJD(1) + t41 * qJD(2) + t576;
t492 = t651 / 0.2e1 + t73 * mrSges(7,1) / 0.2e1 - t74 * mrSges(7,2) / 0.2e1;
t473 = t666 * t464 + (t372 / 0.2e1 + mrSges(6,3) * t683) * t365 + (mrSges(6,3) * t677 - t508) * t364 + t497 * t723 + (-t429 + t497) * t722 + (-t467 * t97 + t470 * t98 + (-t615 + t616) * t465) * t721 + (t529 * t470 + (-t365 * t465 + t85) * t467 + t522 * t364) * t719 + t526;
t485 = -m(4) * t370 / 0.2e1 + (-t466 * t76 - t469 * t75) * t720 + t230 * t686 + t231 * t679;
t11 = t473 + 0.2e1 * (-m(6) / 0.4e1 - m(5) / 0.4e1) * t299 + (t563 + t631 / 0.2e1) * t611 + t485;
t58 = -t598 + (t589 - t604) * t470 + (m(5) + m(4)) * qJ(3) + t581 * mrSges(6,3) + m(7) * (t470 * t522 + t447) + m(6) * (t447 + t613) + t666;
t491 = -t11 * qJD(1) - t58 * qJD(2) - t576 / 0.2e1;
t130 = pkin(5) * t407 - t531 * t681 + t533 * t687 + t507;
t251 = (t407 / 0.2e1 + t514) * t470;
t478 = pkin(10) * t737 + t407 * t688 + t409 * t678 + t533 * t680 + (t411 - t531) * t685;
t30 = (t656 / 0.2e1 + t512) * t469 + (-0.3e1 / 0.4e1 * t652 + t513) * t466 + (t736 + t478) * t467 + t516 + t517;
t475 = pkin(5) * t137 / 0.2e1 + t466 * t714 + t140 * t680 + t407 * t735 + t511 + (t693 - t533 / 0.4e1) * t254 + (t692 + t531 / 0.4e1) * t253;
t8 = t475 + t726 * pkin(10) + (0.3e1 / 0.4e1 * t653 + t570 + t691) * t364 + t492;
t489 = t8 * qJD(1) - t30 * qJD(2) + t251 * qJD(4) + t130 * qJD(5);
t252 = t407 * t677 + t501;
t79 = t501 + t508;
t62 = t700 + t494 / 0.2e1 + t486 + t737 * t470;
t53 = t504 / 0.2e1 + m(6) * t560 + t484;
t31 = -t650 / 0.2e1 + (-t652 / 0.4e1 + t513) * t466 + t512 * t469 + t478 * t467 - t516 + t517 + t729 * t470;
t26 = t365 * t694 + t495 / 0.2e1 + t487 + t737 * t364;
t24 = t253 * t564 + t254 * t565 + t364 * t514 + t509;
t21 = t480 - t515;
t17 = t175 * t677 + t467 * t551 + t231 * t555 + t345 * t683 + t346 * t676 + (t406 * t682 + t490) * t463 + t477 - t479;
t13 = -t190 * t677 + t204 * t553 + (t467 * t519 + t470 * t520) * t719 + t205 * t555 + (t467 * t562 + t563) * t611 + t482 + t526 + (-t615 / 0.2e1 + t616 / 0.2e1) * mrSges(6,3);
t10 = mrSges(6,1) * t538 + t559 * t627 + t432 + t473 - t485 - t575 + (m(5) + m(6)) * t299 / 0.2e1;
t9 = -t475 - t530 * t702 + t729 * t364 + t492 + (t253 * t565 + t469 * t568 + t510) * pkin(10);
t5 = t481 * t467 + t476 - t488 + t733;
t1 = t474 + t149 * t678 + t150 * t685 + (t338 * t680 + t336 * t685 + (Ifges(6,4) / 0.4e1 + t691 + t653 / 0.4e1) * t470 + t502) * t364 + t76 * t564 + Ifges(6,3) * t560 - mrSges(7,3) * t624 / 0.2e1 + t545 * t365 + (-t637 / 0.4e1 - t574 / 0.4e1 + t496 + t511) * t470 + (t636 / 0.4e1 + t147 * t680 + t146 * t685 + t500) * t467 + t558 * t658 + Ifges(6,6) * t538 + (t455 + t597) * t611 / 0.4e1 + (t551 - t605 / 0.2e1 + t527 * t719) * pkin(10) - t725;
t18 = [qJD(2) * t3 + qJD(3) * t14 - qJD(4) * t15 + qJD(5) * t4 + qJD(6) * t7, t10 * qJD(3) + t17 * qJD(4) + t1 * qJD(5) + t5 * qJD(6) + t644 + (t433 + t435 + t437 + t458 * t317 - t299 * t406 + t76 * t392 + t75 * t394 - t370 * mrSges(4,1) - t370 * mrSges(3,1) - t105 * t372 + t369 * mrSges(4,3) + t297 * mrSges(5,2) - t299 * mrSges(5,1) + t281 * t231 + t282 * t230 + t314 * t704 + t336 * t705 + (-t120 * mrSges(6,3) + t465 * t345 + t548) * t470 + (-t256 / 0.2e1 + t150 * t681 + t149 * t686 + t119 * mrSges(6,3) + (-t346 + t175) * t465) * t467 + ((Ifges(6,5) * t683 + Ifges(6,6) * t676 + (-mrSges(4,2) + mrSges(5,3)) * qJ(3) + t665) * t468 + (-Ifges(5,5) - pkin(2) * mrSges(4,2) + t549 - t472 * mrSges(5,3) + (t335 / 0.2e1 - t410 / 0.2e1) * t467) * t471) * t463 + 0.2e1 * (t105 * t610 + t281 * t75 + t282 * t76) * t719 + 0.2e1 * (-t299 * t458 + (-t119 * t467 + t120 * t470) * t465) * t721 + 0.2e1 * (qJ(3) * t297 + t299 * t472) * t722 + 0.2e1 * (-pkin(2) * t370 + qJ(3) * t369) * t723 - t635) * qJD(2), qJD(2) * t10 + qJD(4) * t53 + qJD(5) * t26 + qJD(6) * t24 + t623, -t621 + t17 * qJD(2) + t53 * qJD(3) + t13 * qJD(5) + t21 * qJD(6) + (t521 + t571) * t467 * t664, t634 + t1 * qJD(2) + t26 * qJD(3) + t13 * qJD(4) + (-t97 * mrSges(6,2) + t528 * mrSges(7,3) + pkin(5) * t190 + t146 * t679 + t147 * t686 + t507 * t364 + t408 * t701 + t583 + t724 * t98 + (m(7) * t528 + t591 - t606) * pkin(10)) * qJD(5) + t9 * qJD(6), t625 + t5 * qJD(2) + t24 * qJD(3) + t21 * qJD(4) + t9 * qJD(5) + (-mrSges(7,1) * t57 - mrSges(7,2) * t56 + t139) * qJD(6); qJD(3) * t11 + qJD(4) * t16 + qJD(5) * t2 + qJD(6) * t6 - t644, qJD(3) * t58 - qJD(5) * t28 + qJD(6) * t32, t62 * qJD(5) + t79 * qJD(6) - t491, t577 + t622, t62 * qJD(3) + t31 * qJD(6) + t499 + (mrSges(6,2) * t610 - t597 / 0.2e1 + t339 * t686 + t337 * t679 + pkin(5) * t373 + t455 + (t505 + t588 - t603) * pkin(10) + (t465 * t724 - Ifges(6,5) + t507) * t470 + t728 * mrSges(7,3)) * qJD(5), t79 * qJD(3) + t31 * qJD(5) + (-mrSges(7,1) * t282 - mrSges(7,2) * t281 + t374) * qJD(6) + t498; -qJD(2) * t11 - qJD(4) * t44 + qJD(5) * t27 - qJD(6) * t23 - t623, t59 * qJD(5) - t78 * qJD(6) + t491, 0, -t503, t523, qJD(6) * t407 - t524; -qJD(2) * t16 + qJD(3) * t44 - qJD(5) * t12 + qJD(6) * t22 + t621, t577 - t622, t503, qJD(5) * t674 (mrSges(7,3) * t540 + t371 + t494 - t535) * qJD(5) + t252 * qJD(6) + t493, qJD(5) * t252 + qJD(6) * t371 - t525; -qJD(2) * t2 - qJD(3) * t27 + qJD(4) * t12 - qJD(6) * t8 - t634, -qJD(3) * t59 + qJD(6) * t30 - t499, -t523, -t251 * qJD(6) - t493, -t130 * qJD(6) (pkin(10) * t405 - t530) * qJD(6) - t489; -qJD(2) * t6 + qJD(3) * t23 - qJD(4) * t22 + qJD(5) * t8 - t625, qJD(3) * t78 - qJD(5) * t30 - t498, t524, qJD(5) * t251 + t525, t489, 0;];
Cq  = t18;
