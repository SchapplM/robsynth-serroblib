% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:26
% EndTime: 2019-03-09 18:02:57
% DurationCPUTime: 21.99s
% Computational Cost: add. (62040->633), mult. (118914->832), div. (0->0), fcn. (145966->10), ass. (0->370)
t646 = sin(qJ(3));
t647 = sin(qJ(2));
t649 = cos(qJ(3));
t650 = cos(qJ(2));
t443 = t646 * t647 - t649 * t650;
t444 = -t646 * t650 - t649 * t647;
t591 = sin(pkin(11));
t592 = cos(pkin(11));
t333 = t591 * t443 + t592 * t444;
t334 = -t592 * t443 + t591 * t444;
t392 = sin(qJ(5));
t648 = cos(qJ(5));
t274 = t392 * t333 + t648 * t334;
t385 = -t650 * pkin(2) - pkin(1);
t353 = t443 * pkin(3) + t385;
t283 = -pkin(4) * t334 + t353;
t713 = -t648 * t333 + t334 * t392;
t795 = m(6) * t283 - mrSges(6,1) * t274 + mrSges(6,2) * t713;
t391 = sin(qJ(6));
t393 = cos(qJ(6));
t486 = mrSges(7,1) * t393 - t391 * mrSges(7,2);
t545 = t647 * pkin(7);
t462 = -t647 * pkin(8) - t545;
t547 = t650 * pkin(7);
t463 = t650 * pkin(8) + t547;
t345 = t646 * t462 + t649 * t463;
t406 = -t443 * qJ(4) + t345;
t706 = t649 * t462 - t646 * t463;
t719 = t444 * qJ(4) + t706;
t227 = t592 * t406 + t591 * t719;
t199 = t334 * pkin(9) + t227;
t748 = -t591 * t406 + t592 * t719;
t762 = t333 * pkin(9) + t748;
t782 = t648 * t199 + t392 * t762;
t783 = t782 * t486;
t786 = t782 * mrSges(6,1);
t94 = t199 * t392 - t648 * t762;
t791 = t94 * mrSges(6,2);
t793 = -t783 - t786 + t791;
t676 = mrSges(6,1) / 0.2e1;
t792 = (t486 / 0.2e1 + t676) * t782 - t791 / 0.2e1;
t685 = t791 / 0.2e1 - t783 / 0.2e1 - t786 / 0.2e1;
t546 = t649 * pkin(2);
t384 = t546 + pkin(3);
t487 = t591 * t646;
t361 = -pkin(2) * t487 + t592 * t384;
t358 = pkin(4) + t361;
t488 = t592 * t646;
t363 = pkin(2) * t488 + t591 * t384;
t317 = t392 * t358 + t648 * t363;
t790 = t317 * t94;
t524 = t592 * pkin(3);
t383 = t524 + pkin(4);
t525 = pkin(3) * t591;
t364 = t392 * t383 + t648 * t525;
t789 = t364 * t94;
t771 = t391 * t94;
t769 = t393 * t94;
t774 = t782 * t94;
t316 = t648 * t358 - t392 * t363;
t304 = -pkin(5) - t316;
t773 = t304 * t782;
t362 = t648 * t383 - t392 * t525;
t359 = -pkin(5) - t362;
t772 = t359 * t782;
t605 = t393 * mrSges(7,2);
t608 = t391 * mrSges(7,1);
t378 = t605 + t608;
t174 = t378 * t713;
t784 = t782 * t174;
t270 = Ifges(6,5) * t274;
t630 = Ifges(7,4) * t393;
t381 = Ifges(7,1) * t391 + t630;
t551 = t393 * t381;
t631 = Ifges(7,4) * t391;
t380 = Ifges(7,2) * t393 + t631;
t559 = t391 * t380;
t469 = t559 / 0.2e1 - t551 / 0.2e1;
t379 = Ifges(7,5) * t391 + Ifges(7,6) * t393;
t722 = t713 * t379;
t725 = Ifges(6,6) * t713;
t629 = Ifges(7,2) * t391;
t484 = -t629 + t630;
t732 = Ifges(7,6) * t713 + t484 * t274;
t750 = t393 * t732;
t633 = Ifges(7,1) * t393;
t485 = -t631 + t633;
t731 = Ifges(7,5) * t713 + t485 * t274;
t751 = t391 * t731;
t758 = t722 / 0.2e1 - t725 + t750 / 0.2e1 + t751 / 0.2e1;
t437 = -Ifges(4,5) * t443 + Ifges(5,5) * t334 + Ifges(4,6) * t444 + Ifges(5,6) * t333 - t469 * t274 + t270 + t758;
t718 = t706 * mrSges(4,2);
t723 = t345 * mrSges(4,1);
t752 = t227 * mrSges(5,1);
t753 = t748 * mrSges(5,2);
t780 = t437 - t718 - t723 - t752 - t753 + t793;
t438 = t444 * pkin(3);
t779 = m(5) * t438;
t778 = pkin(5) * t782;
t766 = -t316 * t782 - t790;
t765 = m(6) * (-t362 * t782 - t789);
t607 = t391 * mrSges(7,3);
t728 = mrSges(7,2) * t713;
t178 = -t274 * t607 - t728;
t604 = t393 * mrSges(7,3);
t729 = mrSges(7,1) * t713;
t181 = -t274 * t604 + t729;
t268 = t274 * mrSges(6,2);
t621 = t713 * mrSges(6,1);
t183 = t268 + t621;
t277 = -t333 * mrSges(5,1) + t334 * mrSges(5,2);
t412 = t334 ^ 2;
t413 = t333 ^ 2;
t675 = Ifges(7,3) / 0.2e1;
t538 = t675 + Ifges(6,2) / 0.2e1;
t124 = -pkin(5) * t274 - pkin(10) * t713 + t283;
t55 = t124 * t393 - t391 * t782;
t56 = t391 * t124 + t393 * t782;
t619 = t274 * mrSges(6,3);
t668 = -t274 / 0.2e1;
t465 = t485 * t713;
t115 = -Ifges(7,5) * t274 + t465;
t556 = t393 * t115;
t464 = t484 * t713;
t112 = -Ifges(7,6) * t274 + t464;
t566 = t391 * t112;
t703 = t713 / 0.2e1;
t687 = Ifges(6,1) * t703 + Ifges(6,4) * t274 - t566 / 0.2e1 + t556 / 0.2e1;
t386 = Ifges(7,5) * t393;
t627 = Ifges(7,6) * t391;
t692 = t386 - t627;
t739 = t378 * t274;
t741 = t274 / 0.2e1;
t652 = t393 / 0.2e1;
t654 = -t391 / 0.2e1;
t730 = -t713 / 0.2e1;
t763 = 0.2e1 * t730;
t743 = t652 * t731 + t654 * t732 + t692 * t703 + Ifges(6,4) * t763 + (Ifges(6,2) + Ifges(7,3)) * t668;
t764 = t283 * t183 + t353 * t277 + t385 * (-t444 * mrSges(4,1) - t443 * mrSges(4,2)) + t55 * t181 + t56 * t178 + (t412 - t413) * Ifges(5,4) - t444 ^ 2 * Ifges(4,4) + (Ifges(4,4) * t443 + (Ifges(4,1) - Ifges(4,2)) * t444) * t443 + (t619 + t739) * t94 + (-Ifges(5,1) + Ifges(5,2)) * t333 * t334 + t784 + (-mrSges(6,3) * t94 + t668 * t692 + t687) * t274 + (Ifges(6,1) * t741 - t274 * t538 + t743) * t713;
t760 = t750 / 0.4e1 + t751 / 0.4e1;
t759 = -t718 / 0.2e1 - t723 / 0.2e1 - t752 / 0.2e1 - t753 / 0.2e1;
t755 = -t731 / 0.4e1;
t754 = -t732 / 0.4e1;
t745 = -t227 * t361 + t363 * t748;
t677 = m(5) * pkin(3);
t744 = (-t227 * t592 + t591 * t748) * t677;
t721 = t713 * t486;
t740 = -t721 / 0.2e1;
t697 = (t386 / 0.2e1 - t627 / 0.2e1) * t274;
t653 = -t393 / 0.2e1;
t493 = -t484 * t653 - t485 * t654 - t469;
t733 = -t722 / 0.4e1 + t725 / 0.2e1;
t644 = pkin(5) * t713;
t187 = -t274 * pkin(10) + t644;
t727 = mrSges(6,3) * t713;
t390 = t393 ^ 2;
t609 = t390 * mrSges(7,3);
t389 = t391 ^ 2;
t610 = t389 * mrSges(7,3);
t714 = -t609 / 0.2e1 - t610 / 0.2e1;
t678 = pkin(5) / 0.2e1;
t549 = t389 + t390;
t291 = -t333 * pkin(4) - t438;
t388 = t647 * pkin(2);
t286 = t388 + t291;
t356 = t388 - t438;
t550 = t391 * t178 / 0.2e1 + t181 * t652;
t138 = t187 + t291;
t125 = t388 + t138;
t57 = t125 * t393 + t771;
t58 = t391 * t125 - t769;
t679 = m(7) / 0.2e1;
t681 = m(6) / 0.2e1;
t682 = m(5) / 0.2e1;
t711 = t356 * t682 + t286 * t681 + (t391 * t58 + t393 * t57) * t679 + t550;
t710 = t178 * t652 + t181 * t654;
t709 = -t559 / 0.4e1 + t551 / 0.4e1;
t305 = pkin(10) + t317;
t367 = (-t591 * t649 - t488) * pkin(2);
t368 = (t592 * t649 - t487) * pkin(2);
t326 = t392 * t367 + t648 * t368;
t578 = t713 * t391;
t179 = mrSges(7,2) * t274 - mrSges(7,3) * t578;
t555 = t393 * t179;
t577 = t713 * t393;
t182 = -mrSges(7,1) * t274 - mrSges(7,3) * t577;
t562 = t391 * t182;
t468 = t555 / 0.2e1 - t562 / 0.2e1;
t61 = t138 * t393 + t771;
t62 = t391 * t138 - t769;
t480 = -t61 * t391 + t62 * t393;
t482 = t391 * t55 - t393 * t56;
t528 = t604 / 0.2e1;
t530 = -t607 / 0.2e1;
t490 = t62 * t528 + t61 * t530 + t685;
t567 = t363 * t333;
t568 = t361 * t334;
t572 = t317 * t713;
t573 = t316 * t274;
t575 = t304 * t739;
t325 = -t648 * t367 + t368 * t392;
t617 = t325 * t94;
t661 = t325 / 0.2e1;
t708 = (t227 * t368 + t367 * t748 + t745) * t682 + (t326 * t782 + t617 + t766) * t681 + (t480 * t305 - t482 * t326 + t617 + t773) * t679 + t575 / 0.2e1 + t174 * t661 + t490 + (-t568 / 0.2e1 + t567 / 0.2e1 + t367 * t333 / 0.2e1 + t368 * t334 / 0.2e1) * mrSges(5,3) + (t713 * t661 + t326 * t741 - t572 / 0.2e1 - t573 / 0.2e1) * mrSges(6,3) + t468 * t326 + t759;
t297 = t325 * t486;
t318 = t326 * t610;
t319 = t326 * t609;
t324 = t325 * mrSges(6,1);
t616 = t326 * mrSges(6,2);
t705 = t367 * mrSges(5,1) - t368 * mrSges(5,2) + (-t646 * mrSges(4,1) - t649 * mrSges(4,2)) * pkin(2) - t297 + t318 + t319 - t324 - t616;
t660 = -t359 / 0.2e1;
t489 = mrSges(7,3) * (-t390 / 0.2e1 - t389 / 0.2e1);
t700 = -t486 - mrSges(6,1);
t434 = t646 * t444;
t691 = t714 * t713;
t673 = -t62 / 0.2e1;
t674 = t61 / 0.2e1;
t689 = mrSges(7,1) * t674 + mrSges(7,2) * t673;
t688 = -t58 * mrSges(7,2) / 0.2e1 + t57 * mrSges(7,1) / 0.2e1;
t68 = t187 * t393 + t771;
t69 = t391 * t187 - t769;
t479 = -t391 * t68 + t393 * t69;
t684 = -t270 / 0.2e1 + t739 * t678 + t733;
t680 = -m(7) / 0.2e1;
t672 = mrSges(4,3) * pkin(2);
t669 = t274 / 0.4e1;
t667 = -t274 / 0.4e1;
t664 = -t304 / 0.2e1;
t663 = t305 / 0.2e1;
t662 = -t316 / 0.2e1;
t360 = pkin(10) + t364;
t659 = t360 / 0.2e1;
t658 = -t362 / 0.2e1;
t657 = t364 / 0.2e1;
t656 = t380 / 0.4e1;
t655 = -t381 / 0.4e1;
t643 = pkin(5) * t378;
t278 = -mrSges(5,1) * t334 - mrSges(5,2) * t333;
t2 = m(7) * (t55 * t57 + t56 * t58 + t774) - pkin(1) * (t647 * mrSges(3,1) + t650 * mrSges(3,2)) + t58 * t179 + t57 * t182 + (-Ifges(3,2) + Ifges(3,1)) * t650 * t647 + (-t647 ^ 2 + t650 ^ 2) * Ifges(3,4) + (m(4) * t385 + t443 * mrSges(4,1) - t444 * mrSges(4,2)) * t388 + t764 + (m(5) * t353 + t278) * t356 + t795 * t286;
t626 = t2 * qJD(1);
t3 = m(7) * (t55 * t61 + t56 * t62 + t774) - t353 * t779 - t278 * t438 + t62 * t179 + t61 * t182 + t764 + t795 * t291;
t618 = t3 * qJD(1);
t602 = t57 * t391;
t601 = t58 * t393;
t581 = t274 * t391;
t177 = -mrSges(7,3) * t581 - t728;
t580 = t274 * t393;
t180 = -mrSges(7,3) * t580 + t729;
t9 = t69 * t179 + t56 * t177 + t68 * t182 + t55 * t180 + m(7) * (t55 * t68 + t56 * t69 + t774) + t784 + t94 * t739 + (t283 * mrSges(6,1) + t743) * t713 - (-t283 * mrSges(6,2) + t697 + (-Ifges(6,1) / 0.2e1 + t538) * t713 - t687) * t274;
t598 = t9 * qJD(1);
t595 = t94 * t713;
t594 = t94 * t378;
t17 = (t174 + t727) * t713 + (t413 + t412) * mrSges(5,3) + (t555 - t562 + t619) * t274 + m(7) * (-t482 * t274 + t595) + m(6) * (t274 * t782 + t595) + m(5) * (t227 * t334 + t333 * t748);
t590 = qJD(1) * t17;
t175 = t713 * t380;
t176 = t381 * t713;
t10 = t94 * t721 + t55 * t179 - t56 * t182 + (t482 * mrSges(7,3) + t112 * t653 - t176 * t652 + t379 * t741 + (t115 - t175) * t654) * t713;
t589 = t10 * qJD(1);
t500 = t549 * t274;
t526 = mrSges(7,3) * t549 * t741 + t740;
t411 = (t333 * t361 + t334 * t363) * t682 + (t274 * t317 - t316 * t713) * t681 + (t304 * t713 + t305 * t500) * t679 + t526;
t461 = -t183 - t277;
t20 = t411 + t461 - t711;
t588 = t20 * qJD(1);
t497 = t549 * t360;
t409 = (t274 * t364 - t362 * t713) * t681 + (t274 * t497 + t359 * t713) * t679 + (t592 * t333 + t591 * t334) * t677 / 0.2e1;
t418 = t291 * t681 + (t391 * t62 + t393 * t61) * t679 - t779 / 0.2e1 + t550;
t22 = t721 / 0.2e1 - t409 + t418 - t461 + t714 * t274;
t587 = t22 * qJD(1);
t423 = (t391 * t69 + t393 * t68) * t680 + mrSges(6,2) * t668 + t177 * t654 + t180 * t653;
t426 = -t268 / 0.2e1 + (pkin(10) * t500 - t644) * t679 + t526;
t25 = mrSges(6,1) * t763 + t423 + t426;
t585 = t25 * qJD(1);
t471 = -t605 / 0.2e1 - t608 / 0.2e1;
t460 = t471 * t274;
t27 = t460 - t468;
t584 = t27 * qJD(1);
t576 = t304 * t721;
t574 = t304 * t378;
t571 = t359 * t721;
t570 = t359 * t739;
t569 = t359 * t378;
t561 = t391 * t305;
t560 = t391 * t360;
t553 = t393 * t305;
t552 = t393 * t360;
t548 = mrSges(7,3) * t601;
t544 = t362 * t619;
t543 = t364 * t727;
t535 = t181 * t560;
t533 = t713 * t675;
t529 = t607 / 0.2e1;
t519 = t177 * t659;
t514 = -t561 / 0.2e1;
t513 = -t560 / 0.2e1;
t509 = t182 * t653;
t508 = t553 / 0.2e1;
t507 = -t552 / 0.2e1;
t504 = t662 + t658;
t499 = t549 * t316;
t498 = t549 * t326;
t496 = t549 * t362;
t491 = t58 * t528 + t57 * t530 + t685;
t481 = t601 - t602;
t478 = t524 * t334 * mrSges(5,3);
t477 = t333 * mrSges(5,3) * t525;
t476 = -t643 / 0.2e1 + t493;
t396 = -t765 / 0.2e1 + (t481 * t360 + t772) * t680 - t570 / 0.2e1 - t744 / 0.2e1 + t544 / 0.2e1 + t543 / 0.2e1 + t535 / 0.2e1 + t178 * t507 + t478 / 0.2e1 - t477 / 0.2e1 - t759;
t427 = t434 * t672;
t428 = -t434 / 0.2e1;
t4 = t427 / 0.2e1 + t428 * t672 - t548 / 0.2e1 + t178 * t508 + t181 * t514 + t57 * t529 + t396 - t685 + t708;
t44 = -m(7) * (t304 * t325 + t305 * t498) - m(6) * (-t316 * t325 + t317 * t326) - m(5) * (t361 * t367 + t363 * t368) - t705;
t475 = t4 * qJD(1) - t44 * qJD(2);
t433 = -t316 * mrSges(6,2) + mrSges(7,3) * t499 + t700 * t317;
t45 = m(7) * (t304 * t317 + t305 * t499) + t433;
t414 = (t479 * t305 - t482 * t316 + t773 + t790) * t680 + t739 * t664 - t317 * t174 / 0.2e1;
t441 = t710 * pkin(10) + t709 * t274 - t684 + t760;
t417 = t441 + (t481 * pkin(10) - t778) * t679 + t491;
t6 = (-t69 * mrSges(7,3) / 0.2e1 + t179 * t662 - t305 * t177 / 0.2e1 + t754) * t393 + (t180 * t663 + t68 * mrSges(7,3) / 0.2e1 + t316 * t182 / 0.2e1 + t755) * t391 + (-t379 / 0.4e1 + Ifges(6,6) / 0.2e1) * t713 - (Ifges(6,5) / 0.2e1 + t709) * t274 + t417 + t414 + t792;
t474 = -t6 * qJD(1) + t45 * qJD(2);
t424 = (t656 - t633 / 0.4e1) * t393 + (t381 / 0.4e1 + t630 / 0.2e1 - t629 / 0.4e1) * t391;
t451 = t533 + t386 * t669 - t594 / 0.2e1;
t453 = (t668 + t667) * Ifges(7,6) + t112 / 0.4e1 + t176 / 0.4e1;
t470 = -t115 / 0.4e1 + t175 / 0.4e1 + Ifges(7,5) * t741;
t11 = -t576 / 0.2e1 + (t179 * t663 + t453) * t391 + (t182 * t663 + t470) * t393 + (-t305 * t489 + t424) * t713 + t451 + t688;
t206 = t493 + t574;
t473 = -t11 * qJD(1) + t206 * qJD(2);
t459 = t471 * t316;
t458 = t471 * t326;
t457 = t471 * t362;
t400 = (t304 * t364 + t305 * t496 + t316 * t497 + t317 * t359) * t679 + t700 * (t317 / 0.2e1 + t657) + (-t316 - t362) * t489;
t425 = t297 / 0.2e1 - t318 / 0.2e1 - t319 / 0.2e1 + t324 / 0.2e1 + (-pkin(5) * t325 + pkin(10) * t498) * t680;
t30 = (t326 / 0.2e1 + t504) * mrSges(6,2) + t400 + t425;
t432 = -t362 * mrSges(6,2) + mrSges(7,3) * t496 + t700 * t364;
t79 = m(7) * (t359 * t364 + t360 * t496) + t432;
t436 = Ifges(6,5) * t741 + t69 * t528 + t68 * t530 + t551 * t669 + t559 * t667 + t685 - t733 + t760;
t401 = t436 + (t479 * t360 - t482 * t362 + t772 + t789) * t679 + t174 * t657 - t739 * t660;
t440 = m(7) * (t480 * pkin(10) - t778);
t8 = (mrSges(7,3) * t673 + t274 * t655 - pkin(10) * t178 / 0.2e1 + t362 * t179 / 0.2e1 + t519 + t754) * t393 + (mrSges(7,3) * t674 + t274 * t656 + pkin(10) * t181 / 0.2e1 + t182 * t658 - t360 * t180 / 0.2e1 + t755) * t391 - t440 / 0.2e1 + t401 + t684 + t792;
t456 = t8 * qJD(1) + t30 * qJD(2) + t79 * qJD(3);
t120 = (t664 + t660) * t378 + t458 - t493;
t13 = -t571 / 0.2e1 + (t179 * t659 + t453) * t391 + (t182 * t659 + t470) * t393 + (-t360 * t489 + t424) * t713 + t451 + t689;
t244 = t493 + t569;
t454 = -t13 * qJD(1) - t120 * qJD(2) + t244 * qJD(3);
t452 = Ifges(7,3) * t730 - t68 * mrSges(7,1) / 0.2e1 + t69 * mrSges(7,2) / 0.2e1;
t122 = (t678 + t664) * t378 + t459 - t493;
t435 = -t566 / 0.4e1 + t556 / 0.4e1 + t578 * t655 - t380 * t577 / 0.4e1 + t692 * t667 + t594 / 0.2e1 - (t176 + t464) * t391 / 0.4e1 + (-t175 + t465) * t393 / 0.4e1 + (t529 + t530) * t56;
t410 = t435 + (t179 * t654 + t489 * t713 + t509) * pkin(10) + pkin(5) * t740;
t16 = t410 + t452 - t697;
t188 = (t678 + t660) * t378 + t457 - t493;
t279 = -t493 + t643;
t439 = t16 * qJD(1) - t122 * qJD(2) - t188 * qJD(3) - t279 * qJD(5);
t419 = t435 + t533 + t697;
t335 = t569 / 0.2e1;
t282 = t574 / 0.2e1;
t189 = t335 + t457 + t476;
t123 = t282 + t459 + t476;
t121 = t282 + t335 + t458 + t493;
t29 = -t616 / 0.2e1 + t504 * mrSges(6,2) + t400 - t425;
t28 = t460 + t468;
t26 = t713 * t676 - t621 / 0.2e1 - t423 + t426;
t24 = t409 + t418 + t526;
t23 = t411 + t711;
t15 = t410 + Ifges(7,5) * t580 / 0.2e1 - Ifges(7,6) * t581 / 0.2e1 - t452;
t14 = t182 * t507 + t179 * t513 + t419 + t571 / 0.2e1 + t691 * t360 + t689;
t12 = t179 * t514 + t419 + t576 / 0.2e1 + (t509 + t691) * t305 + t688;
t7 = t393 * t519 + t180 * t513 + t440 / 0.2e1 + t441 + t401 + t490 + t468 * t362;
t5 = t177 * t508 + t180 * t514 + t468 * t316 - t414 + t417 + t436;
t1 = t710 * t305 + (t428 + t434 / 0.2e1) * t672 - t396 + t437 + t491 + t708;
t18 = [qJD(2) * t2 + qJD(3) * t3 + qJD(4) * t17 + qJD(5) * t9 + qJD(6) * t10, t626 + (m(4) * (-t345 * t649 + t646 * t706) * pkin(2) + t780 + mrSges(3,2) * t545 + t178 * t553 + m(5) * t745 + t548 + (-t572 - t573) * mrSges(6,3) + (t567 - t568) * mrSges(5,3) + t575 + m(7) * (t305 * t481 + t773) + m(6) * t766 + t443 * mrSges(4,3) * t546 + t427 - Ifges(3,6) * t647 + Ifges(3,5) * t650 - mrSges(7,3) * t602 - t181 * t561 - mrSges(3,1) * t547) * qJD(2) + t1 * qJD(3) + t23 * qJD(4) + t5 * qJD(5) + t12 * qJD(6), t618 + t1 * qJD(2) + (t744 + t765 + t178 * t552 - t535 - t544 - t543 + m(7) * (t360 * t480 + t772) + t570 + t477 - t478 + t480 * mrSges(7,3) + t780) * qJD(3) + t24 * qJD(4) + t7 * qJD(5) + t14 * qJD(6), qJD(2) * t23 + qJD(3) * t24 + qJD(5) * t26 + qJD(6) * t28 + t590, t5 * qJD(2) + t7 * qJD(3) + t26 * qJD(4) + t15 * qJD(6) + t598 + (-(-Ifges(6,5) + t469) * t274 + (-m(7) * t782 - t739) * pkin(5) + (m(7) * t479 + t393 * t177 - t391 * t180) * pkin(10) + t479 * mrSges(7,3) + t758 + t793) * qJD(5), t589 + t12 * qJD(2) + t14 * qJD(3) + t28 * qJD(4) + t15 * qJD(5) + (-t56 * mrSges(7,1) - t55 * mrSges(7,2) - t722) * qJD(6); qJD(3) * t4 + qJD(4) * t20 - qJD(5) * t6 - qJD(6) * t11 - t626, -qJD(3) * t44 + qJD(5) * t45 + qJD(6) * t206 (m(7) * (t325 * t359 + t326 * t497) + m(6) * (-t325 * t362 + t326 * t364) + (t367 * t592 + t368 * t591) * t677 + t705) * qJD(3) + t29 * qJD(5) + t121 * qJD(6) + t475, t588, t29 * qJD(3) + (m(7) * (-pkin(5) * t317 + pkin(10) * t499) + t433) * qJD(5) + t123 * qJD(6) + t474, t121 * qJD(3) + t123 * qJD(5) + (-t305 * t486 + t692) * qJD(6) + t473; -qJD(2) * t4 - qJD(4) * t22 + qJD(5) * t8 - qJD(6) * t13 - t618, qJD(5) * t30 - qJD(6) * t120 - t475, qJD(5) * t79 + qJD(6) * t244, -t587 (m(7) * (-pkin(5) * t364 + pkin(10) * t496) + t432) * qJD(5) + t189 * qJD(6) + t456, t189 * qJD(5) + (-t360 * t486 + t692) * qJD(6) + t454; -qJD(2) * t20 + qJD(3) * t22 - qJD(5) * t25 - qJD(6) * t27 - t590, -t588, t587, 0, -t585, -qJD(6) * t378 - t584; qJD(2) * t6 - qJD(3) * t8 + qJD(4) * t25 + qJD(6) * t16 - t598, -qJD(3) * t30 - qJD(6) * t122 - t474, -qJD(6) * t188 - t456, t585, -t279 * qJD(6) (-pkin(10) * t486 + t692) * qJD(6) + t439; qJD(2) * t11 + qJD(3) * t13 + qJD(4) * t27 - qJD(5) * t16 - t589, qJD(3) * t120 + qJD(5) * t122 - t473, qJD(5) * t188 - t454, t584, -t439, 0;];
Cq  = t18;
