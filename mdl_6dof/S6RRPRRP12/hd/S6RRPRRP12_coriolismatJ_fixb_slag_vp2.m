% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP12_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:29
% EndTime: 2019-03-09 12:50:55
% DurationCPUTime: 14.41s
% Computational Cost: add. (20953->743), mult. (39436->940), div. (0->0), fcn. (37874->6), ass. (0->382)
t438 = cos(qJ(2));
t437 = sin(qJ(2));
t579 = t437 * qJ(3);
t710 = pkin(2) + pkin(8);
t368 = -t438 * t710 - pkin(1) - t579;
t709 = pkin(3) + pkin(7);
t399 = t709 * t437;
t681 = cos(qJ(4));
t385 = t681 * t399;
t436 = sin(qJ(4));
t250 = -t368 * t436 + t385;
t580 = t436 * t438;
t215 = pkin(9) * t580 + t250;
t435 = sin(qJ(5));
t251 = t368 * t681 + t436 * t399;
t542 = t438 * t681;
t216 = -pkin(9) * t542 + t251;
t680 = cos(qJ(5));
t539 = t680 * t216;
t105 = t215 * t435 + t539;
t582 = t435 * t216;
t106 = t215 * t680 - t582;
t500 = t680 * t681;
t382 = t435 * t436 - t500;
t537 = t680 * t436;
t383 = t435 * t681 + t537;
t714 = m(7) / 0.2e1;
t716 = m(6) / 0.2e1;
t349 = t383 * t438;
t633 = t349 * mrSges(6,3);
t297 = mrSges(6,1) * t437 + t633;
t634 = t349 * mrSges(7,2);
t298 = -mrSges(7,1) * t437 - t634;
t347 = t435 * t580 - t438 * t500;
t691 = t382 / 0.2e1;
t627 = t437 * mrSges(7,3);
t637 = t347 * mrSges(7,2);
t291 = t627 + t637;
t636 = t347 * mrSges(6,3);
t294 = -mrSges(6,2) * t437 + t636;
t736 = t291 + t294;
t740 = mrSges(7,2) + mrSges(6,3);
t767 = t383 / 0.2e1;
t768 = -t383 / 0.2e1;
t722 = -t297 * t767 - t298 * t768 - t691 * t736 + t740 * (t347 * t691 + t349 * t767);
t560 = mrSges(5,3) * t580;
t388 = t437 * mrSges(5,1) + t560;
t507 = mrSges(5,3) * t542;
t390 = -t437 * mrSges(5,2) - t507;
t540 = t681 * t390;
t682 = t438 / 0.2e1;
t688 = -t436 / 0.2e1;
t727 = -t540 / 0.2e1 - t388 * t688 + (-t436 ^ 2 - t681 ^ 2) * mrSges(5,3) * t682;
t429 = t437 * qJ(6);
t189 = pkin(4) * t437 + t215;
t583 = t435 * t189;
t97 = t539 + t583;
t85 = t429 + t97;
t96 = t189 * t680 - t582;
t86 = -t437 * pkin(5) - t96;
t773 = -((t106 - t96) * t383 + (t105 - t97) * t382) * t716 - ((t106 + t86) * t383 + (t105 - t85) * t382) * t714 + t727 - t722;
t765 = Ifges(7,4) + Ifges(6,5);
t764 = Ifges(6,6) - Ifges(7,6);
t697 = -t347 / 0.2e1;
t502 = t681 * t710;
t463 = -pkin(9) * t681 - t502;
t426 = t436 * t710;
t515 = pkin(9) * t436 + t426;
t728 = t435 * t463 - t515 * t680;
t706 = t728 / 0.2e1;
t329 = Ifges(6,4) * t347;
t203 = Ifges(6,2) * t349 + t329;
t656 = Ifges(7,5) * t347;
t178 = -Ifges(7,1) * t349 + t437 * Ifges(7,4) - t656;
t180 = -Ifges(6,1) * t349 + t437 * Ifges(6,5) + t329;
t737 = t180 + t178;
t772 = t203 + t737;
t265 = t383 * mrSges(7,1) + t382 * mrSges(7,3);
t419 = t436 * pkin(4) + qJ(3);
t491 = pkin(5) * t383 + qJ(6) * t382;
t221 = t491 + t419;
t677 = m(7) * t221;
t749 = t265 + t677;
t400 = t709 * t438;
t363 = pkin(4) * t542 + t400;
t620 = qJ(6) * t349;
t672 = pkin(5) * t347;
t492 = t620 - t672;
t142 = t492 + t363;
t660 = Ifges(6,4) * t349;
t176 = Ifges(6,2) * t347 + t437 * Ifges(6,6) - t660;
t196 = -mrSges(7,1) * t349 - mrSges(7,3) * t347;
t197 = -mrSges(6,1) * t349 + mrSges(6,2) * t347;
t202 = -Ifges(7,3) * t349 + t656;
t262 = -mrSges(7,1) * t382 + mrSges(7,3) * t383;
t263 = -mrSges(6,1) * t382 - mrSges(6,2) * t383;
t655 = Ifges(7,5) * t383;
t267 = -Ifges(7,3) * t382 - t655;
t373 = Ifges(7,5) * t382;
t268 = Ifges(7,3) * t383 - t373;
t376 = Ifges(6,4) * t383;
t271 = Ifges(6,2) * t382 - t376;
t659 = Ifges(6,4) * t382;
t272 = -Ifges(6,2) * t383 - t659;
t273 = -Ifges(7,1) * t383 - t373;
t275 = -Ifges(6,1) * t383 + t659;
t280 = t435 * t515 + t680 * t463;
t503 = t764 * t382 - t765 * t383;
t326 = Ifges(7,5) * t349;
t174 = t437 * Ifges(7,6) - Ifges(7,3) * t347 - t326;
t204 = Ifges(7,1) * t347 - t326;
t205 = Ifges(6,1) * t347 + t660;
t759 = t205 + t204 + t174;
t274 = -Ifges(7,1) * t382 + t655;
t276 = -Ifges(6,1) * t382 - t376;
t763 = t274 + t276;
t770 = -t740 * (t280 * t697 + t349 * t706) - t419 * t197 / 0.2e1 - t363 * t263 / 0.2e1 - t349 * t272 / 0.4e1 + t347 * t267 / 0.4e1 - t221 * t196 / 0.2e1 - t142 * t262 / 0.2e1 - t503 * t437 / 0.4e1 - (t271 + t763) * t347 / 0.4e1 + (t275 + t273 + t268) * t349 / 0.4e1 + (-t176 / 0.4e1 + t759 / 0.4e1) * t382 + (-t202 / 0.4e1 + t772 / 0.4e1) * t383;
t769 = -mrSges(6,1) / 0.2e1;
t766 = m(6) * t419;
t630 = t383 * mrSges(7,2);
t544 = t630 / 0.2e1;
t545 = -t630 / 0.2e1;
t252 = t545 + t544;
t761 = qJD(2) * t252;
t760 = qJD(6) * t252;
t757 = (-t276 / 0.2e1 - t274 / 0.2e1 + t267 / 0.2e1 - t271 / 0.2e1) * t383 + (-t275 / 0.2e1 - t273 / 0.2e1 - t268 / 0.2e1 + t272 / 0.2e1) * t382 + t221 * t262 + t419 * t263;
t756 = 0.2e1 * t714;
t750 = mrSges(7,3) - mrSges(6,2);
t663 = t86 + t96;
t428 = m(7) * qJ(6) + mrSges(7,3);
t747 = qJD(5) * t428;
t746 = t428 * qJD(6);
t396 = t437 * pkin(2) - qJ(3) * t438;
t381 = t437 * pkin(8) + t396;
t386 = t681 * t400;
t190 = pkin(4) * t438 + t386 + (-pkin(9) * t437 - t381) * t436;
t260 = t681 * t381 + t436 * t400;
t543 = t437 * t681;
t218 = pkin(9) * t543 + t260;
t103 = t190 * t680 - t435 * t218;
t104 = t435 * t190 + t680 * t218;
t91 = qJ(6) * t438 + t104;
t92 = -t438 * pkin(5) - t103;
t743 = t104 * mrSges(6,2) / 0.2e1 + t103 * t769 + t92 * mrSges(7,1) / 0.2e1 - t91 * mrSges(7,3) / 0.2e1;
t673 = pkin(4) * t435;
t557 = t673 / 0.2e1;
t741 = mrSges(7,1) + mrSges(6,1);
t478 = Ifges(5,5) * t436 + Ifges(5,6) * t681;
t739 = t478 * t437;
t581 = t436 * t437;
t346 = t435 * t581 - t437 * t500;
t348 = -t435 * t543 - t437 * t537;
t738 = -t280 * t348 + t346 * t728;
t493 = -pkin(2) * t438 - t579;
t394 = -pkin(1) + t493;
t395 = t438 * mrSges(4,2) - t437 * mrSges(4,3);
t735 = m(4) * t394 + t395;
t551 = Ifges(5,4) * t681;
t397 = -Ifges(5,2) * t436 + t551;
t661 = Ifges(5,4) * t436;
t398 = Ifges(5,1) * t681 - t661;
t564 = t681 / 0.2e1;
t686 = t436 / 0.2e1;
t469 = t397 * t564 + t398 * t686;
t504 = t347 * t765 + t349 * t764;
t562 = t680 * pkin(4);
t421 = -t562 - pkin(5);
t584 = t421 * t383;
t418 = qJ(6) + t673;
t586 = t418 * t382;
t483 = t584 - t586;
t731 = (Ifges(3,4) + Ifges(4,6)) * t437;
t516 = -t297 / 0.2e1 + t298 / 0.2e1;
t200 = -t347 * mrSges(7,1) + t349 * mrSges(7,3);
t729 = m(7) * t142 + t200;
t726 = m(7) * t85 + t736;
t725 = m(7) * t86 - t297 + t298;
t724 = m(6) * t97 + t726;
t723 = -m(6) * t96 + t725;
t648 = Ifges(6,3) * t438;
t650 = Ifges(7,6) * t346;
t651 = Ifges(6,6) * t346;
t654 = Ifges(7,2) * t438;
t657 = Ifges(6,5) * t348;
t658 = Ifges(7,4) * t348;
t721 = t650 / 0.2e1 - t651 / 0.2e1 - t657 / 0.2e1 - t658 / 0.2e1 + t648 / 0.2e1 + t654 / 0.2e1 - t743;
t717 = m(5) / 0.2e1;
t715 = -m(7) / 0.2e1;
t713 = m(6) * pkin(4);
t712 = m(7) * pkin(4);
t711 = mrSges(6,3) / 0.2e1;
t708 = t200 / 0.2e1;
t707 = t265 / 0.2e1;
t705 = t280 / 0.2e1;
t292 = -t346 * mrSges(7,2) + mrSges(7,3) * t438;
t704 = -t292 / 0.2e1;
t295 = mrSges(6,1) * t438 + t348 * mrSges(6,3);
t703 = t295 / 0.2e1;
t626 = t438 * mrSges(7,1);
t635 = t348 * mrSges(7,2);
t296 = -t626 - t635;
t702 = t296 / 0.2e1;
t699 = -t346 / 0.2e1;
t698 = t346 / 0.2e1;
t696 = t347 / 0.2e1;
t695 = -t348 / 0.2e1;
t694 = -t349 / 0.2e1;
t693 = t349 / 0.2e1;
t692 = -t382 / 0.2e1;
t689 = t418 / 0.2e1;
t684 = t437 / 0.2e1;
t678 = m(7) * t105;
t676 = m(7) * t491;
t675 = m(7) * t728;
t674 = m(7) * t383;
t669 = t96 * mrSges(6,2);
t668 = t96 * mrSges(7,3);
t667 = t97 * mrSges(6,1);
t666 = t97 * mrSges(7,1);
t664 = -t85 + t97;
t649 = Ifges(5,3) * t438;
t645 = t105 * mrSges(6,1);
t644 = t105 * mrSges(7,1);
t643 = t106 * mrSges(6,2);
t642 = t106 * mrSges(7,3);
t641 = t728 * mrSges(6,1);
t640 = t728 * mrSges(7,1);
t639 = t280 * mrSges(6,2);
t638 = t280 * mrSges(7,3);
t632 = t382 * mrSges(7,2);
t631 = t382 * t97;
t629 = t383 * t96;
t628 = t436 * mrSges(5,2);
t625 = t438 * mrSges(5,2);
t563 = t681 * pkin(4);
t362 = (-t563 - t709) * t437;
t141 = t346 * pkin(5) + t348 * qJ(6) + t362;
t173 = -Ifges(7,5) * t348 + Ifges(7,6) * t438 + Ifges(7,3) * t346;
t175 = -Ifges(6,4) * t348 - Ifges(6,2) * t346 + Ifges(6,6) * t438;
t177 = -Ifges(7,1) * t348 + Ifges(7,4) * t438 + Ifges(7,5) * t346;
t179 = -Ifges(6,1) * t348 - Ifges(6,4) * t346 + Ifges(6,5) * t438;
t198 = mrSges(7,1) * t346 + mrSges(7,3) * t348;
t199 = mrSges(6,1) * t346 - mrSges(6,2) * t348;
t201 = -t347 * mrSges(6,1) - t349 * mrSges(6,2);
t259 = -t436 * t381 + t386;
t293 = -mrSges(6,2) * t438 - t346 * mrSges(6,3);
t479 = Ifges(5,2) * t681 + t661;
t338 = Ifges(5,6) * t438 + t437 * t479;
t466 = t479 * t438;
t339 = Ifges(5,6) * t437 - t466;
t480 = Ifges(5,1) * t436 + t551;
t340 = Ifges(5,5) * t438 + t437 * t480;
t467 = t480 * t438;
t341 = Ifges(5,5) * t437 - t467;
t550 = t681 * mrSges(5,1);
t476 = -t550 + t628;
t367 = t476 * t437;
t387 = mrSges(5,1) * t438 - mrSges(5,3) * t581;
t389 = mrSges(5,3) * t543 - t625;
t498 = -t542 / 0.2e1;
t499 = t543 / 0.2e1;
t521 = -t580 / 0.2e1;
t522 = t581 / 0.2e1;
t5 = (t179 + t177) * t694 + t340 * t521 + t341 * t522 + t338 * t498 + t339 * t499 + m(5) * (t250 * t259 + t251 * t260 - t399 * t400) + m(6) * (t103 * t96 + t104 * t97 + t362 * t363) + m(7) * (t141 * t142 + t85 * t91 + t86 * t92) + ((-t478 + 0.2e1 * Ifges(3,4)) * t682 + (Ifges(3,1) - Ifges(4,3)) * t684 + t399 * t476 - pkin(1) * mrSges(3,2) - t394 * mrSges(4,3) + Ifges(4,6) * t438 + (Ifges(4,2) - Ifges(4,3) / 0.2e1 - Ifges(3,2) / 0.2e1) * t437) * t438 + ((Ifges(3,1) - Ifges(3,2) + Ifges(7,2) + Ifges(5,3) + Ifges(6,3)) * t682 - t731 / 0.2e1 - pkin(1) * mrSges(3,1) - t394 * mrSges(4,2)) * t437 + t735 * t396 + t737 * t695 + (t739 + t648 + t649 + t650 - t651 + t654 - t657 - t658 - t731) * t684 + t400 * t367 + t260 * t390 + t250 * t387 + t259 * t388 + t251 * t389 + t362 * t201 + t363 * t199 + t104 * t294 + t96 * t295 + t86 * t296 + t103 * t297 + t92 * t298 + t91 * t291 + t85 * t292 + t97 * t293 + t142 * t198 + t141 * t200 + t175 * t696 + t173 * t697 + t174 * t698 + t176 * t699 + (t764 * t347 - t765 * t349) * t682;
t624 = t5 * qJD(1);
t195 = -t349 * pkin(5) - t347 * qJ(6);
t567 = pkin(4) * t580;
t167 = t195 - t567;
t452 = t142 * t196 + t363 * t197 + t504 * t684 + t85 * t634 - t96 * t636;
t477 = -t436 * mrSges(5,1) - mrSges(5,2) * t681;
t468 = t438 * t477;
t482 = -Ifges(5,5) * t542 + Ifges(5,6) * t580;
t520 = t580 / 0.2e1;
t6 = t176 * t693 + t202 * t697 + t339 * t520 + t341 * t498 + t400 * t468 + t482 * t684 + t97 * t633 + t86 * t637 + t452 + t469 * t438 ^ 2 + (-m(6) * t363 - t201) * t567 + (-t388 + t560) * t251 + (t507 + t390) * t250 + t729 * t167 + t724 * t106 + t723 * t105 + t772 * t696 + t759 * t694;
t623 = t6 * qJD(1);
t7 = (-t204 / 0.2e1 + t176 / 0.2e1 - t205 / 0.2e1 - t174 / 0.2e1) * t349 + (t178 / 0.2e1 + t203 / 0.2e1 + t86 * mrSges(7,2) + t180 / 0.2e1 - t202 / 0.2e1) * t347 + t452 + t729 * t195 + (t725 + t633) * t97 + t726 * t96;
t622 = t7 * qJD(1);
t621 = qJD(6) * t674;
t36 = t437 * t291 + t349 * t200 + m(7) * (t142 * t349 + t437 * t85);
t619 = qJD(1) * t36;
t457 = mrSges(6,2) * t699 + mrSges(7,3) * t698 + t695 * t741;
t448 = (-pkin(5) * t348 + qJ(6) * t346) * t714 + t457;
t517 = -t294 / 0.2e1 - t291 / 0.2e1;
t555 = -mrSges(7,2) / 0.2e1 - mrSges(6,3) / 0.2e1;
t13 = (t349 * t555 + t663 * t715 - t516) * t383 + (t347 * t555 + t664 * t715 - t517) * t382 + t448;
t618 = t13 * qJD(1);
t18 = t388 * t581 + t723 * t348 + t724 * t346 + (-t540 + m(5) * (t250 * t436 - t251 * t681) - t735) * t437;
t616 = t18 * qJD(1);
t606 = t348 * t382;
t593 = t382 * t435;
t249 = t383 * t346;
t571 = 0.2e1 * t545;
t525 = t383 * t684;
t494 = t525 + t695;
t193 = t494 * m(7);
t570 = t193 * qJD(1);
t569 = mrSges(6,3) * t673;
t568 = t713 / 0.2e1;
t566 = -t681 / 0.2e1;
t561 = t418 * t634;
t559 = -t676 / 0.2e1;
t558 = -t673 / 0.2e1;
t556 = -mrSges(7,1) / 0.2e1 + t769;
t554 = mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t553 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t552 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t548 = -t637 / 0.2e1;
t547 = t632 / 0.2e1;
t546 = -t631 / 0.2e1;
t541 = t681 * t259;
t538 = t680 * t383;
t536 = t620 / 0.2e1;
t266 = t383 * mrSges(6,1) - t382 * mrSges(6,2);
t514 = t349 * t569;
t509 = mrSges(6,3) * t562;
t508 = mrSges(6,3) * t546;
t506 = -t562 / 0.2e1;
t501 = t550 / 0.2e1;
t495 = t347 * t509;
t261 = -pkin(5) * t382 + t383 * qJ(6);
t17 = t749 * t261 + t757;
t446 = t86 * t545 + t85 * t547 + t631 * t711 - t770;
t440 = t446 - t517 * t280 + (t142 * t261 + t195 * t221 - t280 * t664) * t714 + (-t629 / 0.2e1 + t546) * mrSges(7,2) + t261 * t708 + t195 * t707 + t508 + (t663 * t714 + t516) * t728;
t455 = (-pkin(5) * t92 + qJ(6) * t91) * t715 + pkin(5) * t702 + qJ(6) * t704;
t3 = (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t438 + t552 * t346 + t440 - t553 * t348 + t455 + t743;
t490 = t3 * qJD(1) + t17 * qJD(2);
t443 = (t346 * t418 + t348 * t421) * t714 + (t346 * t435 - t348 * t680) * t568 + mrSges(5,1) * t522 + mrSges(5,2) * t499 + t457;
t10 = t443 + t773;
t488 = t10 * qJD(1);
t445 = t556 * t347 + t554 * t349 + t400 * t717 + (t363 + t738) * t716 + (t142 + t738) * t714;
t447 = -m(5) * (t436 * t260 + t541) / 0.2e1 - m(6) * (-t103 * t382 + t104 * t383) / 0.2e1 + (t382 * t92 + t383 * t91) * t715 + t387 * t566;
t16 = t438 * t501 + (-t389 / 0.2e1 - t625 / 0.2e1) * t436 + (t704 - t293 / 0.2e1 + t555 * t346) * t383 + (-t296 / 0.2e1 + t703 + t555 * t348) * t382 + t445 + t447;
t470 = t265 + t266;
t456 = -t470 + t477;
t57 = mrSges(4,3) + (m(5) + m(4)) * qJ(3) + t766 + t677 - t456;
t487 = -qJD(1) * t16 - qJD(2) * t57;
t451 = (t142 * t382 + t221 * t349 + t437 * t728) * t714 + t265 * t693 + t200 * t691;
t481 = t92 * t715 + t626 / 0.2e1;
t34 = -mrSges(7,2) * t494 + t451 + t481;
t83 = t749 * t382;
t486 = qJD(1) * t34 + qJD(2) * t83;
t485 = -t105 * t280 + t106 * t728;
t475 = m(7) * (-pkin(5) * t728 + qJ(6) * t280);
t473 = t538 + t593;
t225 = t563 + t261;
t12 = -qJ(3) * t476 + t479 * t686 + t480 * t566 - t469 + t757 + (t266 + t766) * t563 + t749 * t225;
t439 = pkin(4) * t266 * t521 + t397 * t520 + t398 * t498 + (t280 * t97 + (t363 * t681 - t419 * t580) * pkin(4) + t485) * t716 + t167 * t707 + t225 * t708 + t629 * t711 + (t142 * t225 + t167 * t221 + t280 * t85 + t485) * t714 + t201 * t563 / 0.2e1 - t739 / 0.4e1 - t400 * t476 / 0.2e1 + qJ(3) * t468 / 0.2e1 + t736 * t705 - (-t467 + t341) * t436 / 0.4e1 - (-t466 + t339) * t681 / 0.4e1 + (t714 * t86 - t716 * t96 + t516) * t728 + t727 * t710 + t740 * (t105 * t692 + t106 * t768);
t442 = t721 + t293 * t557 + Ifges(5,5) * t522 + Ifges(5,6) * t499 + (t103 * t680 + t104 * t435) * t568 + (t418 * t91 + t421 * t92) * t714 + t649 / 0.2e1 + t259 * mrSges(5,1) / 0.2e1 - t260 * mrSges(5,2) / 0.2e1 + t292 * t689 + t421 * t702 + t562 * t703;
t2 = t86 * t544 - t85 * t632 / 0.2e1 - t439 + t508 + t442 + t770;
t465 = -t2 * qJD(1) + t12 * qJD(2);
t453 = t627 + (t418 * t437 + t85) * t714;
t38 = -t678 / 0.2e1 + t453;
t391 = m(7) * t418 + mrSges(7,3);
t464 = -qJD(1) * t38 - qJD(4) * t391 + t761;
t458 = t750 * t562 - t741 * t673;
t166 = -(t418 * t680 + t421 * t435) * t712 - t458;
t450 = (-(-t418 + t673) * t280 + (t421 + t562) * t728) * t714;
t20 = -t475 / 0.2e1 + ((-t421 / 0.2e1 - pkin(5) / 0.2e1 + t506) * t383 + (t689 - qJ(6) / 0.2e1 + t558) * t382) * mrSges(7,2) + t450 + t741 * (t706 - t728 / 0.2e1) - t750 * (t705 - t280 / 0.2e1);
t454 = m(7) * (pkin(4) * t473 + t483);
t55 = t559 - t454 / 0.2e1;
t441 = (t418 * t96 + t421 * t97 + (t435 * t86 + t680 * t85) * pkin(4)) * t715 + t669 / 0.2e1 - t668 / 0.2e1 + t667 / 0.2e1 + t666 / 0.2e1 - t561 / 0.2e1 + t421 * t548 + t297 * t557 + t298 * t558 - t514 / 0.2e1 + t495 / 0.2e1 + t736 * t506;
t449 = (-pkin(5) * t105 + qJ(6) * t106) * t714 - t645 / 0.2e1 - t644 / 0.2e1 - t643 / 0.2e1 + t642 / 0.2e1;
t9 = mrSges(7,2) * t536 + pkin(5) * t548 + t441 + t449;
t460 = -t9 * qJD(1) + t20 * qJD(2) - t55 * qJD(3) - t166 * qJD(4);
t40 = t627 + 0.2e1 * (t539 / 0.4e1 + t583 / 0.4e1 + t429 / 0.2e1 - t97 / 0.4e1) * m(7);
t459 = qJD(1) * t40 + qJD(4) * t428 + t747 + t761;
t364 = mrSges(7,3) + (qJ(6) + 0.2e1 * t557) * m(7);
t192 = m(7) * t525 + t348 * t714;
t122 = t571 + t675;
t78 = t675 / 0.2e1 + m(7) * t706 + t571;
t42 = t454 / 0.2e1 + t559 - t470;
t39 = t85 * t756 + t291;
t37 = t637 + t678 / 0.2e1 + t453;
t33 = t437 * t545 - t635 / 0.2e1 + t451 - t481;
t19 = pkin(5) * t544 + qJ(6) * t547 - t640 / 0.2e1 + t638 / 0.2e1 - t639 / 0.2e1 - t641 / 0.2e1 + t475 / 0.2e1 + t556 * t728 + t554 * t280 + (-t584 / 0.2e1 + t586 / 0.2e1 + (-t593 / 0.2e1 - t538 / 0.2e1) * pkin(4)) * mrSges(7,2) + t450 + t503;
t15 = t296 * t691 + t295 * t692 + t389 * t686 + (m(4) * pkin(7) + t501 - t628 / 0.2e1 + mrSges(4,1)) * t438 + t445 - t447 + (t293 + t292) * t767 + t740 * (-t249 / 0.2e1 - t606 / 0.2e1);
t14 = (t382 * t664 + t383 * t663) * t714 + t448 + t722;
t11 = t443 - t773;
t8 = (t536 - t672 / 0.2e1) * mrSges(7,2) - t441 + t449 + t504;
t4 = t440 - t455 + t721;
t1 = t446 + t439 + t442;
t21 = [qJD(2) * t5 + qJD(3) * t18 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t36, t15 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + t33 * qJD(6) + t624 + (0.2e1 * (-qJ(3) * t399 - t259 * t502 - t260 * t426) * t717 - mrSges(5,2) * t385 - mrSges(5,3) * t541 + t340 * t564 + t728 * t293 + t728 * t292 - t387 * t502 + (t141 * t221 - t280 * t92 + t728 * t91) * t756 + 0.2e1 * (t103 * t280 + t104 * t728 + t362 * t419) * t716 + t280 * t295 - t280 * t296 + (-t177 / 0.2e1 - t179 / 0.2e1 + t103 * mrSges(6,3) - t92 * mrSges(7,2)) * t382 + t419 * t199 + t362 * t266 + qJ(3) * t367 + t141 * t265 + t221 * t198 + (-t399 * mrSges(5,1) - t338 / 0.2e1 - t710 * t389 - t260 * mrSges(5,3)) * t436 + (-t175 / 0.2e1 + t173 / 0.2e1 - t104 * mrSges(6,3) - t91 * mrSges(7,2)) * t383 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t564 + Ifges(5,6) * t688 + t382 * t553 - t383 * t552 - Ifges(4,4) + Ifges(3,5)) * t438 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t469) * t437 + (m(4) * t493 - t438 * mrSges(3,1) + t437 * mrSges(3,2) + t395) * pkin(7) + t268 * t698 + t272 * t699 + t763 * t695) * qJD(2), t15 * qJD(2) + t11 * qJD(4) + t14 * qJD(5) + t192 * qJD(6) + t616 + 0.2e1 * (t716 + t714) * qJD(3) * (t249 + t606) t623 + t1 * qJD(2) + t11 * qJD(3) + (t642 - t643 + t421 * t637 + t561 + m(7) * (t105 * t421 + t106 * t418) + (-t105 * t680 + t106 * t435) * t713 - t645 - t644 - t250 * mrSges(5,2) - t251 * mrSges(5,1) + t514 - t495 + t482 + t504) * qJD(4) + t8 * qJD(5) + t37 * qJD(6), t622 + t4 * qJD(2) + t14 * qJD(3) + t8 * qJD(4) + (t668 - t669 + m(7) * (-pkin(5) * t97 + qJ(6) * t96) - t667 - t666 + t492 * mrSges(7,2) + t504) * qJD(5) + t39 * qJD(6), qJD(2) * t33 + qJD(3) * t192 + qJD(4) * t37 + qJD(5) * t39 + t619; qJD(3) * t16 - qJD(4) * t2 + qJD(5) * t3 + qJD(6) * t34 - t624, qJD(3) * t57 + qJD(4) * t12 + qJD(5) * t17 + qJD(6) * t83, -t487 (m(7) * (t280 * t418 + t421 * t728) + (t280 * t435 - t680 * t728) * t713 - t640 + t638 - t639 - t641 + t382 * t569 + t383 * t509 + mrSges(5,2) * t502 + mrSges(5,1) * t426 - t478 + t503 - t483 * mrSges(7,2)) * qJD(4) + t19 * qJD(5) + t78 * qJD(6) + t465, t19 * qJD(4) + (t491 * mrSges(7,2) + t503 + (-m(7) * pkin(5) - t741) * t728 - (mrSges(6,2) - t428) * t280) * qJD(5) + t122 * qJD(6) + t490, qJD(4) * t78 + qJD(5) * t122 + t486; -qJD(2) * t16 - qJD(4) * t10 - qJD(5) * t13 + qJD(6) * t193 - t616, t487, 0 (m(7) * t483 - t473 * t713 + t456) * qJD(4) + t42 * qJD(5) - t488 + t621, -t618 + t42 * qJD(4) + (-t470 - t676) * qJD(5) + t621, t570 + 0.2e1 * (qJD(4) / 0.2e1 + qJD(5) / 0.2e1) * t674; qJD(2) * t2 + qJD(3) * t10 - qJD(5) * t9 + qJD(6) * t38 - t623, qJD(5) * t20 - t465 - t760, -qJD(5) * t55 + t488, -qJD(5) * t166 + qJD(6) * t391 ((-pkin(5) * t435 + qJ(6) * t680) * t712 + t458) * qJD(5) + t364 * qJD(6) + t460, qJD(5) * t364 - t464; -qJD(2) * t3 + qJD(3) * t13 + qJD(4) * t9 + qJD(6) * t40 - t622, -qJD(4) * t20 - t490 + t760, qJD(4) * t55 + t618, -t460 + t746, t746, t459; -qJD(2) * t34 - qJD(3) * t193 - qJD(4) * t38 - qJD(5) * t40 - t619, -t486 + (qJD(4) - qJD(5)) * t252, -t570, t464 - t747, -t459, 0;];
Cq  = t21;
