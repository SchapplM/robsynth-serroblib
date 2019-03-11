% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP1
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:39:12
% EndTime: 2019-03-09 11:39:38
% DurationCPUTime: 16.81s
% Computational Cost: add. (33012->708), mult. (64161->914), div. (0->0), fcn. (75188->8), ass. (0->368)
t423 = sin(qJ(4));
t420 = sin(pkin(10));
t421 = cos(pkin(10));
t424 = sin(qJ(2));
t604 = -qJ(3) - pkin(7);
t493 = t604 * t424;
t426 = cos(qJ(2));
t689 = t604 * t426;
t340 = -t420 * t493 + t421 * t689;
t380 = -t420 * t424 + t421 * t426;
t449 = -t380 * pkin(8) + t340;
t625 = cos(qJ(4));
t381 = -t420 * t426 - t421 * t424;
t708 = t420 * t689 + t421 * t493;
t737 = t381 * pkin(8) + t708;
t183 = t423 * t449 + t625 * t737;
t425 = cos(qJ(5));
t422 = sin(qJ(5));
t598 = mrSges(7,2) * t422;
t388 = -mrSges(7,1) * t425 + t598;
t330 = t625 * t380 + t381 * t423;
t720 = t422 * t330;
t758 = t423 * t737 - t625 * t449;
t760 = pkin(5) * t720 + t758;
t771 = -t183 * mrSges(5,2) + t760 * t388;
t459 = t423 * t380 - t381 * t625;
t563 = t459 * t422;
t107 = pkin(5) * t563 - t183;
t769 = t107 * t760;
t620 = pkin(2) * t421;
t407 = pkin(3) + t620;
t621 = pkin(2) * t420;
t360 = t423 * t407 + t621 * t625;
t768 = t183 * t360;
t767 = t183 * t422;
t567 = t183 * t758;
t359 = t407 * t625 - t423 * t621;
t357 = -pkin(4) - t359;
t610 = t425 * pkin(5);
t353 = t357 - t610;
t766 = t353 * t760;
t408 = -pkin(4) - t610;
t765 = t408 * t760;
t764 = t425 * t183;
t597 = mrSges(7,2) * t425;
t600 = mrSges(7,1) * t422;
t391 = t597 + t600;
t222 = t391 * t459;
t599 = mrSges(6,2) * t425;
t392 = mrSges(6,1) * t422 + t599;
t223 = t392 * t459;
t562 = t459 * t425;
t409 = -pkin(2) * t426 - pkin(1);
t355 = -t380 * pkin(3) + t409;
t613 = t459 * pkin(9);
t179 = -pkin(4) * t330 + t355 - t613;
t82 = t425 * t179 - t422 * t758;
t67 = -qJ(6) * t562 + t82;
t55 = -pkin(5) * t330 + t67;
t83 = t179 * t422 + t425 * t758;
t68 = -qJ(6) * t563 + t83;
t721 = t392 * t330;
t722 = t391 * t330;
t719 = t425 * t330;
t727 = mrSges(7,1) * t459 - mrSges(7,3) * t719;
t728 = -mrSges(7,2) * t459 - mrSges(7,3) * t720;
t729 = mrSges(6,1) * t459 - mrSges(6,3) * t719;
t730 = -mrSges(6,2) * t459 - mrSges(6,3) * t720;
t762 = t107 * t722 - t183 * t721 + t760 * t222 + t758 * t223 + t55 * t727 + t68 * t728 + t82 * t729 + t83 * t730;
t761 = t357 * t758;
t757 = (t730 / 0.2e1 + t728 / 0.2e1) * t422;
t756 = -t381 * mrSges(4,1) + t380 * mrSges(4,2);
t670 = -m(7) / 0.2e1;
t755 = -t727 / 0.2e1;
t754 = t727 / 0.2e1;
t358 = pkin(9) + t360;
t548 = qJ(6) + t358;
t345 = t548 * t425;
t753 = t345 * t728;
t752 = t358 * t730;
t418 = t422 ^ 2;
t419 = t425 ^ 2;
t544 = t418 + t419;
t750 = t544 * (mrSges(6,3) + mrSges(7,3));
t415 = Ifges(6,4) * t425;
t684 = -Ifges(6,2) * t422 + t415;
t414 = Ifges(7,4) * t425;
t685 = -Ifges(7,2) * t422 + t414;
t401 = Ifges(7,1) * t422 + t414;
t403 = Ifges(6,1) * t422 + t415;
t716 = t403 + t401;
t748 = t684 + t685 + t716;
t703 = t459 / 0.2e1;
t704 = -t459 / 0.2e1;
t413 = Ifges(6,5) * t425;
t686 = -Ifges(6,6) * t422 + t413;
t412 = Ifges(7,5) * t425;
t687 = -Ifges(7,6) * t422 + t412;
t713 = t686 + t687;
t723 = Ifges(6,3) + Ifges(7,3);
t724 = -t330 / 0.2e1;
t592 = Ifges(6,4) * t422;
t404 = Ifges(6,1) * t425 - t592;
t711 = Ifges(6,5) * t459 + t330 * t404;
t591 = Ifges(7,4) * t422;
t402 = Ifges(7,1) * t425 - t591;
t712 = Ifges(7,5) * t459 + t330 * t402;
t731 = t712 / 0.2e1 + t711 / 0.2e1;
t709 = Ifges(7,6) * t459 + t330 * t685;
t710 = Ifges(6,6) * t459 + t330 * t684;
t733 = -t710 / 0.2e1 - t709 / 0.2e1;
t746 = t422 * t733 + t425 * t731 + 0.2e1 * Ifges(5,4) * t704 + t713 * t703 + (Ifges(5,2) + t723) * t724;
t397 = Ifges(7,2) * t425 + t591;
t399 = Ifges(6,2) * t425 + t592;
t717 = t399 + t397;
t744 = -t717 / 0.4e1;
t743 = pkin(4) * t721;
t740 = qJ(6) * t720;
t739 = t404 + t402;
t603 = -qJ(6) - pkin(9);
t387 = t603 * t422;
t390 = t603 * t425;
t536 = mrSges(7,3) * t563;
t234 = mrSges(7,2) * t330 - t536;
t240 = -mrSges(7,1) * t330 - mrSges(7,3) * t562;
t629 = t422 / 0.2e1;
t460 = t240 * t629 - t425 * t234 / 0.2e1;
t469 = -t422 * t55 + t425 * t68;
t736 = t460 + ((-t387 * t425 + t390 * t422) * t459 + t469) * t670;
t344 = t548 * t422;
t735 = t460 + ((t344 * t425 - t345 * t422) * t459 + t469) * t670;
t734 = t709 / 0.4e1 + t710 / 0.4e1;
t732 = t711 / 0.4e1 + t712 / 0.4e1;
t726 = pkin(5) * t459 - qJ(6) * t719;
t645 = t330 / 0.2e1;
t601 = mrSges(6,1) * t425;
t474 = -mrSges(6,2) * t422 + t601;
t577 = -t474 - mrSges(5,1);
t531 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t532 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t718 = (t422 * t531 - t425 * t532) * t330;
t144 = -Ifges(7,5) * t330 + t402 * t459;
t147 = -Ifges(6,5) * t330 + t404 * t459;
t482 = t532 * t330;
t715 = t482 - t144 / 0.2e1 - t147 / 0.2e1;
t618 = pkin(4) * t459;
t245 = -pkin(9) * t330 + t618;
t669 = m(7) / 0.2e1;
t707 = 0.2e1 * t669;
t671 = m(6) / 0.2e1;
t706 = 0.2e1 * t671;
t700 = Ifges(6,5) + Ifges(7,5);
t699 = Ifges(6,6) + Ifges(7,6);
t495 = -t419 / 0.2e1 - t418 / 0.2e1;
t478 = t495 * mrSges(6,3);
t666 = m(7) * pkin(5);
t528 = mrSges(7,1) + t666;
t393 = Ifges(7,5) * t422 + Ifges(7,6) * t425;
t395 = Ifges(6,5) * t422 + Ifges(6,6) * t425;
t680 = t395 / 0.4e1 + t393 / 0.4e1;
t690 = (-Ifges(5,6) / 0.2e1 + t680) * t459;
t337 = t345 * t425;
t561 = t344 * t422;
t466 = t337 + t561;
t385 = t390 * t425;
t465 = -t387 * t422 - t385;
t683 = (-mrSges(5,2) + t750) * t359 + (t388 + t577) * t360;
t640 = t357 / 0.2e1;
t682 = (-t561 / 0.2e1 - t337 / 0.2e1) * mrSges(7,3) + t474 * t640;
t537 = mrSges(6,3) * t563;
t235 = mrSges(6,2) * t330 - t537;
t630 = -t422 / 0.2e1;
t241 = -mrSges(6,1) * t330 - mrSges(6,3) * t562;
t652 = -t241 / 0.2e1;
t681 = t235 * t630 + t425 * t652;
t679 = t724 + t645;
t678 = -t355 * mrSges(5,2) + Ifges(5,1) * t704 - Ifges(5,4) * t330;
t677 = t183 * t392 / 0.2e1 - t107 * t391 / 0.2e1;
t313 = mrSges(7,1) * t562;
t479 = -mrSges(7,2) * t563 + t313;
t538 = pkin(5) * t562;
t126 = -m(7) * t538 - t479;
t370 = t422 * t528 + t597;
t676 = -qJD(1) * t126 + (qJD(2) + qJD(4)) * t370;
t496 = t403 / 0.2e1 + t401 / 0.2e1;
t499 = t399 / 0.2e1 + t397 / 0.2e1;
t675 = t422 * t499 - t425 * t496 - Ifges(5,5);
t627 = t425 / 0.2e1;
t674 = t627 * t748 + t739 * t629 + t717 * t630;
t673 = m(5) / 0.2e1;
t672 = -m(6) / 0.2e1;
t668 = m(4) * pkin(2);
t667 = m(6) * pkin(9);
t665 = -mrSges(6,1) / 0.2e1;
t664 = mrSges(6,1) / 0.2e1;
t663 = -mrSges(7,1) / 0.2e1;
t662 = -mrSges(6,2) / 0.2e1;
t661 = -mrSges(7,2) / 0.2e1;
t91 = t425 * t245 - t767;
t60 = t726 + t91;
t660 = t60 / 0.2e1;
t417 = t424 * pkin(2);
t356 = -pkin(3) * t381 + t417;
t185 = t245 + t356;
t85 = t425 * t185 - t767;
t659 = t85 / 0.2e1;
t86 = t422 * t185 + t764;
t658 = -t86 / 0.2e1;
t56 = t726 + t85;
t657 = pkin(5) * t56;
t656 = t760 / 0.2e1;
t655 = t235 / 0.2e1;
t654 = -t729 / 0.2e1;
t653 = -t240 / 0.2e1;
t643 = -t344 / 0.2e1;
t642 = -t353 / 0.2e1;
t641 = t353 / 0.2e1;
t639 = -t387 / 0.2e1;
t638 = -t388 / 0.2e1;
t637 = t388 / 0.2e1;
t636 = -t474 / 0.2e1;
t635 = t390 / 0.2e1;
t632 = -t408 / 0.2e1;
t631 = t408 / 0.2e1;
t628 = t422 / 0.4e1;
t626 = t425 / 0.4e1;
t624 = m(7) * t107;
t623 = m(7) * t459;
t622 = m(7) * t359;
t619 = pkin(4) * t758;
t617 = pkin(4) * t392;
t616 = pkin(5) * t422;
t615 = pkin(9) * t730;
t614 = pkin(9) * t729;
t609 = t60 * mrSges(7,3);
t92 = t422 * t245 + t764;
t73 = t92 - t740;
t608 = t73 * mrSges(7,3);
t607 = t91 * mrSges(6,3);
t606 = t92 * mrSges(6,3);
t602 = -t55 + t67;
t596 = mrSges(7,3) * t425;
t586 = t758 * mrSges(5,1);
t326 = t459 * mrSges(5,1);
t530 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t483 = Ifges(5,2) / 0.2e1 + t530;
t552 = t425 * t147;
t553 = t425 * t144;
t141 = -Ifges(6,6) * t330 + t459 * t684;
t557 = t422 * t141;
t138 = -Ifges(7,6) * t330 + t459 * t685;
t558 = t422 * t138;
t70 = t86 - t740;
t3 = m(6) * (t82 * t85 + t83 * t86 - t567) + t56 * t240 + t85 * t241 + t70 * t234 + t86 * t235 + m(7) * (t55 * t56 + t68 * t70 + t769) + (t356 * mrSges(5,2) + Ifges(5,1) * t645 + t746) * t459 + (-mrSges(4,2) * t417 - Ifges(4,4) * t381) * t381 + (-mrSges(4,1) * t417 + Ifges(4,4) * t380 + (-Ifges(4,1) + Ifges(4,2)) * t381) * t380 + (t553 / 0.2e1 + t552 / 0.2e1 - t558 / 0.2e1 - t557 / 0.2e1 + t713 * t724 - t678 - t356 * mrSges(5,1) - t483 * t459) * t330 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t426 + (Ifges(3,1) - Ifges(3,2)) * t424) * t426 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t424) * t424 + (m(4) * t417 + t756) * t409 + t762 + (m(5) * t356 + t326) * t355;
t583 = t3 * qJD(1);
t481 = t531 * t330;
t450 = t138 / 0.2e1 + t141 / 0.2e1 - t481;
t4 = t73 * t234 + t92 * t235 + t60 * t240 + t91 * t241 + m(6) * (t82 * t91 + t83 * t92 - t567) + m(7) * (t55 * t60 + t68 * t73 + t769) + (t355 * mrSges(5,1) + t746) * t459 - (t715 * t425 + t450 * t422 + (-Ifges(5,1) / 0.2e1 + t483) * t459 + t678) * t330 + t762;
t582 = t4 * qJD(1);
t581 = t422 * mrSges(6,3);
t580 = t422 * mrSges(7,3);
t579 = t425 * t86;
t224 = t397 * t459;
t225 = t399 * t459;
t226 = t401 * t459;
t227 = t403 * t459;
t456 = t474 * t459;
t527 = m(7) * t602;
t9 = t83 * t241 + t183 * t456 - t222 * t538 - t107 * t479 - t67 * t234 + t68 * t240 - t55 * t536 - t68 * t527 + ((-t224 / 0.2e1 - t225 / 0.2e1 - t715) * t422 + (t226 / 0.2e1 + t227 / 0.2e1 + t83 * mrSges(6,3) + t68 * mrSges(7,3) - pkin(5) * t624 + t450) * t425) * t459 + (-t537 - t235) * t82;
t578 = t9 * qJD(1);
t565 = t183 * t459;
t10 = (t380 ^ 2 + t381 ^ 2) * mrSges(4,3) + (mrSges(5,3) * t459 + t222 + t223) * t459 + (mrSges(5,3) * t330 + (t234 + t235) * t425 + (-t240 - t241) * t422) * t330 + m(6) * (-t565 + (-t422 * t82 + t425 * t83) * t330) + m(7) * (t107 * t459 + t330 * t469) + m(5) * (t330 * t758 - t565) + m(4) * (-t340 * t380 + t381 * t708);
t576 = qJD(1) * t10;
t490 = t544 * t330;
t429 = -m(5) * (t330 * t360 - t359 * t459) / 0.2e1 + (t357 * t459 + t358 * t490) * t672 + (t330 * t466 + t353 * t459) * t670 - (t380 * t420 + t381 * t421) * t668 / 0.2e1;
t539 = t668 / 0.2e1;
t430 = t757 + (t729 / 0.2e1 + t754) * t425 + t356 * t673 + (t422 * t86 + t425 * t85) * t671 + (t422 * t70 + t425 * t56) * t669 + t424 * t539;
t11 = t326 + (t638 + t474 / 0.2e1) * t459 + (mrSges(7,3) * t495 + mrSges(5,2) + t478) * t330 + t429 + t430 + t756;
t575 = qJD(1) * t11;
t461 = (t636 + t637) * t459 + t645 * t750;
t431 = -t326 / 0.2e1 + (pkin(9) * t490 - t618) * t671 + (t330 * t465 + t408 * t459) * t669 + t461;
t437 = (t422 * t92 + t425 * t91) * t672 + (t422 * t73 + t425 * t60) * t670 + mrSges(5,1) * t704;
t15 = t431 + 0.2e1 * t724 * mrSges(5,2) + (t654 + t755) * t425 - t757 + t437;
t574 = qJD(1) * t15;
t27 = (t422 * t234 - m(7) * (-t422 * t68 - t425 * t55) + t425 * t240) * t459;
t573 = qJD(1) * t27;
t570 = t107 * t422;
t500 = t652 + t653;
t501 = t655 + t234 / 0.2e1;
t529 = -t67 / 0.2e1 + t55 / 0.2e1;
t512 = -t720 / 0.2e1;
t524 = -t597 / 0.2e1;
t546 = mrSges(7,1) * t512 + t330 * t524;
t13 = (mrSges(6,2) * t724 - t501) * t425 + (mrSges(6,1) * t724 + (pkin(5) * t724 + t529) * m(7) - t500) * t422 + t546;
t568 = t13 * qJD(1);
t566 = t758 * t474;
t560 = t387 * t234;
t545 = t544 * mrSges(7,3);
t542 = qJD(5) * t422;
t462 = m(7) * t544 * t704;
t112 = -t623 / 0.2e1 + t462;
t541 = t112 * qJD(1);
t535 = t360 * t669;
t534 = t616 / 0.2e1;
t533 = t661 + t662;
t523 = -t596 / 0.2e1;
t522 = t596 / 0.2e1;
t521 = -t581 / 0.2e1;
t520 = -t580 / 0.2e1;
t519 = t580 / 0.2e1;
t509 = t719 / 0.2e1;
t498 = t399 / 0.4e1 + t397 / 0.4e1;
t497 = -t403 / 0.4e1 - t401 / 0.4e1;
t491 = t602 * t345;
t489 = t544 * t359;
t487 = t665 - t666 / 0.2e1;
t485 = t393 / 0.2e1 + t395 / 0.2e1 - Ifges(5,6);
t475 = t760 * t669 + mrSges(7,1) * t720 / 0.2e1 + mrSges(7,2) * t509;
t468 = -t422 * t85 + t579;
t428 = (t223 / 0.2e1 + t222 / 0.2e1) * t360 + t690 + (t761 - t768) * t671 + (t107 * t360 - t344 * t60 + t345 * t73 + t766) * t669 + t727 * t643 + t753 / 0.2e1 + t722 * t641 + t721 * t640;
t433 = t500 * t359 + (-t358 * t91 - t359 * t82) * t671 + t358 * t654 - t55 * t622 / 0.2e1 + t732;
t434 = t501 * t359 + (t358 * t92 + t359 * t83) * t671 + t752 / 0.2e1 + t68 * t622 / 0.2e1 + t734;
t435 = (t387 * t56 - t390 * t70 + t765) * t670 + t743 / 0.2e1 + t727 * t639 + t728 * t635 + t722 * t632;
t1 = t428 + (t656 - t760 / 0.2e1) * t388 + t619 * t671 - t690 + t679 * Ifges(5,5) + (-t615 / 0.2e1 + (t73 / 0.2e1 - t70 / 0.2e1) * mrSges(7,3) + (t92 / 0.2e1 + t658) * mrSges(6,3) + t658 * t667 + t434 - t734) * t425 + (t614 / 0.2e1 + (-t60 / 0.2e1 + t56 / 0.2e1) * mrSges(7,3) + (-t91 / 0.2e1 + t659) * mrSges(6,3) + t659 * t667 + t433 - t732) * t422 + t435;
t35 = m(7) * (t353 * t360 + t359 * t466) + m(6) * (t357 * t360 + t358 * t489) + t683;
t467 = t1 * qJD(1) + t35 * qJD(2);
t464 = -pkin(5) * t596 + t412 + t413;
t22 = t475 + t735;
t242 = m(7) * t466 + t545;
t463 = -qJD(1) * t22 + qJD(2) * t242;
t458 = -t402 / 0.2e1 - t404 / 0.2e1 + t499;
t457 = -t685 / 0.2e1 - t684 / 0.2e1 - t496;
t454 = (t353 * t562 + t570) * pkin(5);
t444 = (t337 - t385 + (t344 - t387) * t422) * t670 - t545;
t109 = t535 + t444;
t443 = -(t524 - t600 / 0.2e1) * t330 + m(7) * t656;
t25 = t443 + t736;
t289 = m(7) * t465 + t545;
t451 = -qJD(1) * t25 - qJD(2) * t109 + qJD(4) * t289;
t336 = t353 * t391;
t339 = t357 * t392;
t379 = t388 * t616;
t44 = -t379 - t336 - t339 + t457 * t425 + (-t353 * t666 + t458) * t422;
t432 = -t677 - t557 / 0.4e1 - t558 / 0.4e1 + t552 / 0.4e1 + t553 / 0.4e1 + t538 * t637 + t67 * t522 + t55 * t523 + t222 * t534 + (t581 / 0.2e1 + t521) * t83 + (t519 + t520) * t68 - t713 * t330 / 0.4e1 + (-t227 - t226) * t628 + (-t225 - t224) * t626 - t748 * t563 / 0.4e1 + (t744 + t739 / 0.4e1) * t562;
t427 = t432 + t234 * t643 + t345 * t653 + (t459 * t478 + t681) * t358;
t438 = pkin(5) * t755 + t56 * t663 + t70 * mrSges(7,2) / 0.2e1 + t85 * t665 + t86 * mrSges(6,2) / 0.2e1;
t6 = t427 + t313 * t641 + (t491 / 0.2e1 + t454 / 0.2e1 - t657 / 0.2e1) * m(7) + t718 + t438 + (t598 * t642 - t530 + t682) * t459;
t448 = -t6 * qJD(1) + t44 * qJD(2);
t446 = -t336 / 0.2e1 - t339 / 0.2e1 - t379 + t617 / 0.2e1 + t391 * t632;
t31 = (t359 * t533 + t457) * t425 + ((t663 + t665) * t359 + (-t359 / 0.2e1 + t642 + t632) * t666 + t458) * t422 + t446;
t65 = -t617 + t379 + (m(7) * t616 + t391) * t408 + t674;
t436 = (m(7) * t660 + t754) * pkin(5) + mrSges(7,1) * t660 + t73 * t661 + t91 * t664 + t92 * t662;
t7 = -t560 / 0.2e1 + t313 * t632 - (-t687 / 0.4e1 - t686 / 0.4e1) * t330 - (t240 / 0.2e1 - t527 / 0.2e1) * t390 + (t227 / 0.4e1 + t226 / 0.4e1 + t141 / 0.4e1 + t138 / 0.4e1 + pkin(9) * t655 - t481 + (-t222 / 0.2e1 - t624 / 0.2e1) * pkin(5)) * t422 + (t225 / 0.4e1 + t224 / 0.4e1 - t147 / 0.4e1 - t144 / 0.4e1 + pkin(9) * t241 / 0.2e1 + t482 + t529 * mrSges(7,3)) * t425 + (-pkin(9) * t478 + (pkin(4) * t662 + mrSges(7,2) * t631 + t684 / 0.4e1 + t685 / 0.4e1 + mrSges(7,3) * t639 - t497) * t422 + (pkin(4) * t664 - t404 / 0.4e1 - t402 / 0.4e1 - t390 * mrSges(7,3) / 0.2e1 + (m(7) * t632 + t638) * pkin(5) + t498) * t425 + t530) * t459 + t436 + t677;
t442 = t7 * qJD(1) + t31 * qJD(2) - t65 * qJD(4);
t365 = t370 * qJD(6);
t364 = t370 * qJD(5);
t111 = t623 / 0.2e1 + t462;
t110 = t535 - t444;
t32 = m(7) * (t353 + t408) * t534 + (t533 * t425 + (t663 + t487) * t422) * t359 - t446 + t674;
t26 = t443 - t736;
t23 = t475 - t735;
t16 = t431 - t437 + (t728 + t730) * t629 + (t727 + t729) * t627 + t679 * mrSges(5,2);
t14 = t241 * t630 + t235 * t627 + t527 * t629 + (-t599 / 0.2e1 + t487 * t422) * t330 + t546 - t460;
t12 = -t429 + t430 + t461;
t8 = t432 + t479 * t631 + t560 / 0.2e1 + t240 * t635 - pkin(4) * t456 / 0.2e1 + (-t602 * t390 + (t408 * t562 + t570) * pkin(5)) * t669 - t718 + t436 - t544 * mrSges(6,3) * t613 / 0.2e1 + t681 * pkin(9) + (t387 * t519 - t390 * t523 + t530) * t459;
t5 = t479 * t641 + t427 - t438 + (t657 + t491 + t454) * t669 + t723 * t703 + t699 * t512 + t700 * t509 + t682 * t459;
t2 = -t435 + t771 + t680 * t459 + t428 + (t606 / 0.2e1 + t608 / 0.2e1 + t434) * t425 + (-t607 / 0.2e1 - t609 / 0.2e1 + t433) * t422 + mrSges(6,3) * t579 / 0.2e1 + (t636 - mrSges(5,1) / 0.2e1) * t758 + Ifges(5,5) * t645 + t614 * t630 + t615 * t627 + t720 * t744 + t716 * t719 / 0.4e1 + (t710 + t709) * t626 + (t711 + t712) * t628 - (-Ifges(5,5) / 0.2e1 + t497 * t425 + t498 * t422) * t330 + Ifges(5,6) * t704 - t586 / 0.2e1 - t566 / 0.2e1 + t56 * t520 + t85 * t521 + t70 * t522 + (pkin(9) * t468 - t619) * t671;
t17 = [qJD(2) * t3 + qJD(3) * t10 + qJD(4) * t4 - qJD(5) * t9 - qJD(6) * t27, t12 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t23 * qJD(6) + t583 + (t753 + (-t360 * mrSges(5,3) + t485) * t459 + t771 + Ifges(3,5) * t426 - Ifges(3,6) * t424 + Ifges(4,5) * t380 + Ifges(4,6) * t381 + t340 * mrSges(4,1) + 0.2e1 * (-t359 * t758 + t768) * t673 + (-t344 * t56 + t345 * t70 + t766) * t707 - t344 * t727 + (-t85 * mrSges(6,3) - t56 * mrSges(7,3) - t358 * t729 + t731) * t422 + (t358 * t468 + t761) * t706 + (-mrSges(3,1) * t426 + mrSges(3,2) * t424) * pkin(7) + t353 * t722 + t357 * t721 - t708 * mrSges(4,2) + 0.2e1 * (t340 * t421 + t420 * t708) * t539 + (t86 * mrSges(6,3) + t70 * mrSges(7,3) - t733 + t752) * t425 + (-t380 * t620 + t381 * t621) * mrSges(4,3) - t586 - t566 + (-t359 * mrSges(5,3) - t675) * t330) * qJD(2), qJD(2) * t12 + qJD(4) * t16 + qJD(5) * t14 + qJD(6) * t111 + t576, t2 * qJD(2) + t16 * qJD(3) + t8 * qJD(5) + t26 * qJD(6) + t582 + (t408 * t722 + t387 * t727 - t390 * t728 - t743 + m(7) * (t387 * t60 - t390 * t73 + t765) + (-m(6) * pkin(4) + t577) * t758 + (t606 + t608 + (m(6) * t92 + t730) * pkin(9) - t733) * t425 + (-t607 - t609 + (-m(6) * t91 - t729) * pkin(9) + t731) * t422 + t485 * t459 - t675 * t330 + t771) * qJD(4), t5 * qJD(2) + t14 * qJD(3) + t8 * qJD(4) - t578 + (-mrSges(6,1) * t83 - mrSges(6,2) * t82 - mrSges(7,2) * t67 + (-t699 * t425 + (mrSges(7,3) * pkin(5) - t700) * t422) * t459 - t528 * t68) * qJD(5), qJD(2) * t23 + qJD(3) * t111 + qJD(4) * t26 - t573; -qJD(3) * t11 + qJD(4) * t1 + qJD(5) * t6 - qJD(6) * t22 - t583, qJD(4) * t35 - qJD(5) * t44 + qJD(6) * t242, -t575, t32 * qJD(5) + t110 * qJD(6) + t467 + ((t359 * t465 + t360 * t408) * t707 + (-pkin(4) * t360 + pkin(9) * t489) * t706 + t683) * qJD(4), t32 * qJD(4) + (mrSges(7,2) * t344 - t345 * t528 - t358 * t601 + t464) * qJD(5) + (mrSges(6,2) * t358 - t699) * t542 - t448, qJD(4) * t110 + t463; qJD(2) * t11 - qJD(4) * t15 - qJD(5) * t13 + qJD(6) * t112 - t576, t575, 0, -t574, -t568 + (-mrSges(6,2) - mrSges(7,2)) * qJD(5) * t425 + (-mrSges(6,1) - t528) * t542, t541; -qJD(2) * t1 + qJD(3) * t15 - qJD(5) * t7 - qJD(6) * t25 - t582, -qJD(5) * t31 - qJD(6) * t109 - t467, t574, qJD(5) * t65 + qJD(6) * t289 (-mrSges(7,2) * t387 - pkin(9) * t601 + t390 * t528 + t464) * qJD(5) + (mrSges(6,2) * pkin(9) - t699) * t542 - t442, t451; -qJD(2) * t6 + qJD(3) * t13 + qJD(4) * t7 + qJD(6) * t126 + t578, t31 * qJD(4) - t365 + t448, t568, -t365 + t442, 0, -t676; qJD(2) * t22 - qJD(3) * t112 + qJD(4) * t25 - qJD(5) * t126 + t573, qJD(4) * t109 + t364 - t463, -t541, t364 - t451, t676, 0;];
Cq  = t17;
