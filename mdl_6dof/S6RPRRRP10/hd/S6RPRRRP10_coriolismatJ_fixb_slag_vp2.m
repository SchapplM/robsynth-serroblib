% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:23
% EndTime: 2019-03-09 06:30:47
% DurationCPUTime: 14.16s
% Computational Cost: add. (19456->815), mult. (39794->1052), div. (0->0), fcn. (38311->6), ass. (0->424)
t760 = m(6) / 0.2e1;
t808 = Ifges(7,4) + Ifges(6,5);
t814 = Ifges(6,3) + Ifges(7,2);
t462 = cos(qJ(4));
t698 = pkin(4) * t462;
t443 = -pkin(3) - t698;
t460 = sin(qJ(4));
t705 = sin(qJ(5));
t558 = t705 * t460;
t706 = cos(qJ(5));
t405 = -t706 * t462 + t558;
t562 = t706 * t460;
t406 = -t462 * t705 - t562;
t516 = pkin(5) * t405 + qJ(6) * t406;
t239 = t443 + t516;
t649 = t406 * mrSges(7,3);
t654 = t405 * mrSges(7,1);
t272 = t649 + t654;
t790 = m(7) * t239 + t272;
t461 = sin(qJ(3));
t463 = cos(qJ(3));
t427 = pkin(3) * t463 + pkin(8) * t461;
t410 = t462 * t427;
t464 = -pkin(1) - pkin(7);
t541 = -t460 * t464 + pkin(4);
t612 = t461 * t462;
t253 = pkin(9) * t612 + t463 * t541 + t410;
t603 = t463 * t464;
t324 = t460 * t427 + t462 * t603;
t614 = t460 * t461;
t292 = pkin(9) * t614 + t324;
t112 = t253 * t706 - t292 * t705;
t100 = -t463 * pkin(5) - t112;
t113 = t705 * t253 + t706 * t292;
t367 = t406 * t461;
t304 = -mrSges(6,2) * t463 - mrSges(6,3) * t367;
t369 = t405 * t461;
t306 = mrSges(6,1) * t463 - mrSges(6,3) * t369;
t323 = -t460 * t603 + t410;
t584 = t705 * pkin(4);
t438 = t584 + qJ(6);
t585 = t706 * pkin(4);
t442 = -t585 - pkin(5);
t529 = t584 / 0.2e1;
t531 = t585 / 0.2e1;
t646 = t461 * Ifges(5,5);
t566 = -t646 / 0.2e1;
t756 = m(6) * pkin(4);
t588 = t756 / 0.2e1;
t707 = t463 / 0.2e1;
t643 = t463 * mrSges(7,1);
t656 = t369 * mrSges(7,2);
t307 = -t643 + t656;
t733 = t307 / 0.2e1;
t303 = -mrSges(7,2) * t367 + mrSges(7,3) * t463;
t734 = t303 / 0.2e1;
t754 = -mrSges(5,2) / 0.2e1;
t758 = m(7) / 0.2e1;
t727 = t369 / 0.2e1;
t729 = t367 / 0.2e1;
t810 = mrSges(7,1) / 0.2e1;
t99 = qJ(6) * t463 + t113;
t766 = t113 * mrSges(6,2) / 0.2e1 - t112 * mrSges(6,1) / 0.2e1 + t100 * t810 - t99 * mrSges(7,3) / 0.2e1;
t764 = Ifges(7,6) * t729 - Ifges(6,6) * t367 / 0.2e1 + t808 * t727 + t814 * t707 - t766;
t813 = (t100 * t442 + t438 * t99) * t758 + Ifges(5,3) * t707 + t323 * mrSges(5,1) / 0.2e1 + t324 * t754 + t438 * t734 + t442 * t733 + (t112 * t706 + t113 * t705) * t588 + t462 * t566 + Ifges(5,6) * t614 / 0.2e1 + t306 * t531 + t304 * t529 + t764;
t807 = Ifges(6,6) - Ifges(7,6);
t523 = -t808 * t405 + t807 * t406;
t796 = pkin(9) + pkin(8);
t592 = t796 * t462;
t769 = -t558 * t796 + t592 * t706;
t791 = t769 * mrSges(7,1);
t792 = t769 * mrSges(6,1);
t299 = t562 * t796 + t592 * t705;
t804 = t299 * mrSges(7,3);
t805 = t299 * mrSges(6,2);
t812 = t523 - t791 - t792 + t805 - t804;
t811 = t791 / 0.2e1 + t792 / 0.2e1 + t804 / 0.2e1 - t805 / 0.2e1;
t779 = mrSges(7,2) + mrSges(6,3);
t458 = t460 ^ 2;
t459 = t462 ^ 2;
t591 = t458 + t459;
t806 = mrSges(5,3) * t591;
t647 = t461 * mrSges(7,3);
t368 = t406 * t463;
t657 = t368 * mrSges(7,2);
t302 = t647 + t657;
t305 = -mrSges(6,2) * t461 + t368 * mrSges(6,3);
t543 = -t305 / 0.2e1 - t302 / 0.2e1;
t803 = t543 * t299;
t653 = t405 * mrSges(7,2);
t570 = -t653 / 0.2e1;
t257 = t570 + t653 / 0.2e1;
t802 = qJD(3) * t257;
t801 = qJD(6) * t257;
t798 = -pkin(5) * t769 - qJ(6) * t299;
t270 = -mrSges(7,1) * t406 + mrSges(7,3) * t405;
t271 = -mrSges(6,1) * t406 - mrSges(6,2) * t405;
t688 = Ifges(7,5) * t405;
t274 = -Ifges(7,3) * t406 - t688;
t396 = Ifges(7,5) * t406;
t275 = Ifges(7,3) * t405 - t396;
t399 = Ifges(6,4) * t405;
t276 = Ifges(6,2) * t406 - t399;
t278 = -Ifges(7,1) * t405 - t396;
t279 = -Ifges(7,1) * t406 + t688;
t692 = Ifges(6,4) * t406;
t280 = -Ifges(6,1) * t405 + t692;
t281 = -Ifges(6,1) * t406 - t399;
t277 = -Ifges(6,2) * t405 - t692;
t737 = t277 / 0.2e1;
t797 = (t274 / 0.2e1 - t281 / 0.2e1 - t276 / 0.2e1 - t279 / 0.2e1) * t405 + (-t275 / 0.2e1 - t280 / 0.2e1 + t737 - t278 / 0.2e1) * t406 + t239 * t270 + t443 * t271;
t794 = -mrSges(6,1) - mrSges(7,1);
t793 = -mrSges(6,2) + mrSges(7,3);
t697 = pkin(8) * t463;
t700 = pkin(3) * t461;
t418 = qJ(2) - t697 + t700;
t402 = t462 * t418;
t607 = t462 * t463;
t586 = pkin(9) * t607;
t245 = t461 * t541 + t402 - t586;
t450 = t461 * t464;
t576 = t462 * t450;
t268 = t576 + (-pkin(9) * t463 + t418) * t460;
t559 = t705 * t268;
t106 = t245 * t706 - t559;
t94 = -t461 * pkin(5) - t106;
t638 = t106 + t94;
t789 = t281 + t279;
t599 = t305 + t302;
t613 = t460 * t463;
t436 = pkin(4) * t613;
t404 = t436 - t603;
t366 = t405 * t463;
t517 = -pkin(5) * t368 + t366 * qJ(6);
t142 = t404 + t517;
t218 = -pkin(5) * t366 - qJ(6) * t368;
t587 = pkin(4) * t607;
t167 = t218 + t587;
t269 = -pkin(5) * t406 + qJ(6) * t405;
t699 = pkin(4) * t460;
t246 = t269 + t699;
t644 = t462 * mrSges(5,1);
t648 = t460 * mrSges(5,2);
t519 = t644 - t648;
t389 = t519 * t463;
t220 = -mrSges(6,1) * t366 + mrSges(6,2) * t368;
t710 = t461 / 0.4e1;
t714 = t443 / 0.2e1;
t219 = -mrSges(7,1) * t366 - mrSges(7,3) * t368;
t743 = t219 / 0.2e1;
t479 = t142 * t270 / 0.2e1 + t239 * t743 + t404 * t271 / 0.2e1 + t220 * t714 + t523 * t710;
t317 = -t450 * t460 + t402;
t267 = t317 - t586;
t565 = t706 * t268;
t114 = t267 * t705 + t565;
t115 = t267 * t706 - t559;
t539 = t114 * t299 + t115 * t769;
t308 = mrSges(6,1) * t461 + t366 * mrSges(6,3);
t309 = -mrSges(7,1) * t461 - t366 * mrSges(7,2);
t542 = -t308 / 0.2e1 + t309 / 0.2e1;
t740 = t272 / 0.2e1;
t658 = t368 * mrSges(7,1);
t661 = t366 * mrSges(7,3);
t223 = -t658 + t661;
t742 = t223 / 0.2e1;
t560 = t705 * t245;
t107 = t565 + t560;
t453 = t461 * qJ(6);
t93 = t453 + t107;
t788 = (t142 * t246 + t167 * t239 - t299 * t93 + t769 * t94 + t539) * t758 - pkin(3) * t389 / 0.2e1 + t167 * t740 + t246 * t742 + t479 + t542 * t769;
t652 = t406 * mrSges(6,2);
t655 = t405 * mrSges(6,1);
t524 = -t654 / 0.2e1 - t655 / 0.2e1 - t649 / 0.2e1 + t652 / 0.2e1;
t563 = t706 * t405;
t774 = t442 * t405 - t438 * t406;
t787 = t774 * t758 - t648 / 0.2e1 + t644 / 0.2e1 + (-t406 * t705 - t563) * t588 + t524;
t511 = -t323 * t460 + t324 * t462;
t452 = m(7) * qJ(6) + mrSges(7,3);
t785 = qJD(5) * t452;
t784 = t452 * qJD(6);
t328 = t657 / 0.2e1;
t781 = 0.2e1 * t328;
t598 = -t308 + t309;
t455 = Ifges(5,4) * t462;
t777 = -Ifges(5,2) * t460 + t455;
t423 = Ifges(5,1) * t460 + t455;
t414 = mrSges(5,1) * t461 - mrSges(5,3) * t607;
t608 = t462 * t414;
t412 = -mrSges(5,2) * t461 - mrSges(5,3) * t613;
t616 = t460 * t412;
t776 = -t608 / 0.2e1 - t616 / 0.2e1;
t215 = t769 * t366;
t629 = t299 * t368;
t775 = t215 / 0.2e1 + t629 / 0.2e1;
t525 = t807 * t366 + t808 * t368;
t393 = t463 * t423;
t659 = t368 * mrSges(6,1);
t662 = t366 * mrSges(6,2);
t224 = -t659 - t662;
t741 = t224 / 0.2e1;
t772 = pkin(4) * t741 - t393 / 0.4e1;
t771 = t760 + t758;
t770 = (t270 + t271) * t707;
t708 = -t463 / 0.2e1;
t768 = t708 * t806 + t776;
t577 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t578 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t765 = -t577 * t367 - t578 * t369;
t763 = 0.2e1 * m(7);
t762 = 0.2e1 * qJD(3);
t761 = m(5) / 0.2e1;
t759 = -m(7) / 0.2e1;
t757 = -pkin(8) / 0.2e1;
t755 = m(7) * pkin(4);
t753 = -mrSges(6,3) / 0.2e1;
t752 = mrSges(6,3) / 0.2e1;
t751 = t93 / 0.2e1;
t750 = -t94 / 0.2e1;
t748 = -qJ(6) / 0.2e1;
t747 = t107 / 0.2e1;
t746 = -t114 / 0.2e1;
t745 = -t115 / 0.2e1;
t693 = Ifges(6,4) * t366;
t199 = Ifges(6,2) * t368 + t461 * Ifges(6,6) - t693;
t744 = t199 / 0.4e1;
t273 = -t652 + t655;
t739 = t273 / 0.2e1;
t738 = t274 / 0.4e1;
t736 = t277 / 0.4e1;
t726 = -t369 / 0.2e1;
t420 = mrSges(5,1) * t460 + mrSges(5,2) * t462;
t391 = t463 * t420;
t725 = t391 / 0.2e1;
t694 = Ifges(5,4) * t460;
t421 = Ifges(5,2) * t462 + t694;
t392 = t463 * t421;
t724 = -t392 / 0.4e1;
t722 = -t405 / 0.2e1;
t720 = t405 / 0.2e1;
t719 = t405 / 0.4e1;
t718 = -t406 / 0.2e1;
t716 = -t421 / 0.4e1;
t424 = Ifges(5,1) * t462 - t694;
t715 = t424 / 0.4e1;
t713 = -t460 / 0.2e1;
t712 = t460 / 0.2e1;
t709 = t462 / 0.2e1;
t704 = m(7) * t114;
t703 = m(7) * (pkin(5) * t369 + qJ(6) * t367);
t702 = m(7) * t769;
t701 = m(7) * t369;
t454 = Ifges(5,5) * t462;
t689 = Ifges(7,5) * t368;
t685 = Ifges(5,6) * t460;
t680 = t106 * mrSges(6,2);
t679 = t106 * mrSges(7,3);
t678 = t107 * mrSges(6,1);
t677 = t107 * mrSges(7,1);
t674 = t114 * mrSges(6,1);
t673 = t114 * mrSges(7,1);
t672 = t115 * mrSges(6,2);
t671 = t115 * mrSges(7,3);
t660 = t366 * t93;
t651 = t406 * mrSges(7,2);
t650 = t406 * mrSges(6,3);
t645 = t461 * Ifges(5,6);
t403 = -pkin(4) * t614 + t450;
t141 = pkin(5) * t367 - qJ(6) * t369 + t403;
t196 = Ifges(7,5) * t369 + Ifges(7,6) * t463 + Ifges(7,3) * t367;
t198 = Ifges(6,4) * t369 - Ifges(6,2) * t367 + Ifges(6,6) * t463;
t200 = Ifges(7,1) * t369 + Ifges(7,4) * t463 + Ifges(7,5) * t367;
t202 = Ifges(6,1) * t369 - Ifges(6,4) * t367 + Ifges(6,5) * t463;
t221 = mrSges(7,1) * t367 - mrSges(7,3) * t369;
t222 = mrSges(6,1) * t367 + mrSges(6,2) * t369;
t318 = t460 * t418 + t576;
t360 = Ifges(5,6) * t463 - t461 * t777;
t362 = t463 * Ifges(5,5) - t424 * t461;
t390 = t420 * t461;
t411 = -mrSges(5,2) * t463 + mrSges(5,3) * t614;
t413 = mrSges(5,1) * t463 + mrSges(5,3) * t612;
t499 = t454 / 0.2e1 - t685 / 0.2e1 - Ifges(4,4);
t201 = -Ifges(7,1) * t366 + Ifges(7,4) * t461 - t689;
t344 = Ifges(6,4) * t368;
t203 = -Ifges(6,1) * t366 + Ifges(6,5) * t461 + t344;
t544 = -t201 / 0.2e1 - t203 / 0.2e1;
t341 = Ifges(7,5) * t366;
t197 = t461 * Ifges(7,6) - Ifges(7,3) * t368 - t341;
t545 = -t197 / 0.2e1 + t199 / 0.2e1;
t363 = t463 * t424 + t646;
t610 = t462 * t363;
t361 = t463 * t777 + t645;
t617 = t460 * t361;
t5 = (t464 * t391 - t610 / 0.2e1 + t617 / 0.2e1 - qJ(2) * mrSges(4,2) - t499 * t461 + (-m(5) * t464 ^ 2 - Ifges(4,1) + Ifges(4,2) + Ifges(5,3) + t814) * t463 + t765) * t461 + (qJ(2) * mrSges(4,1) + t360 * t713 + t362 * t709 + t578 * t366 + t577 * t368 + t464 * t390 + t499 * t463) * t463 + (-t202 / 0.2e1 - t200 / 0.2e1) * t366 + m(6) * (t106 * t112 + t107 * t113 + t403 * t404) + m(7) * (t100 * t94 + t141 * t142 + t93 * t99) + m(5) * (t317 * t323 + t318 * t324) - t544 * t369 - t545 * t367 - (-t198 / 0.2e1 + t196 / 0.2e1) * t368 + t318 * t411 + t324 * t412 + t317 * t413 + t323 * t414 + t403 * t224 + t404 * t222 + t100 * t309 + t113 * t305 + t106 * t306 + t94 * t307 + t112 * t308 + t99 * t302 + t93 * t303 + t107 * t304 + t142 * t221 + t141 * t223;
t642 = t5 * qJD(1);
t226 = -Ifges(7,3) * t366 + t689;
t227 = Ifges(6,2) * t366 + t344;
t228 = Ifges(7,1) * t368 - t341;
t229 = Ifges(6,1) * t368 + t693;
t634 = t106 * t368;
t472 = mrSges(7,2) * t660 - mrSges(6,3) * t634 + t142 * t219 + t404 * t220 - (-t227 / 0.2e1 - t94 * mrSges(7,2) + t226 / 0.2e1 + t544) * t368 + (t107 * mrSges(6,3) - t228 / 0.2e1 - t229 / 0.2e1 + t545) * t366 + t525 * t461 / 0.2e1;
t6 = m(7) * (t114 * t94 + t115 * t93 + t142 * t167) + (-t464 * t389 + (t392 / 0.2e1 - t363 / 0.2e1 + t566 + t317 * mrSges(5,3)) * t460 + (-t393 / 0.2e1 - t361 / 0.2e1 - t645 / 0.2e1 - t318 * mrSges(5,3) + (m(6) * t404 + t224) * pkin(4)) * t462) * t463 + t317 * t412 - t318 * t414 + t167 * t223 + t472 + m(6) * (-t106 * t114 + t107 * t115) + t598 * t114 + t599 * t115;
t641 = t6 * qJD(1);
t7 = t472 + (m(7) * t142 + t223) * t218 + (m(7) * t94 + t598) * t107 + (m(7) * t93 + t599) * t106;
t640 = t7 * qJD(1);
t639 = -mrSges(4,1) - t519;
t637 = t107 - t93;
t636 = qJD(6) * t701;
t39 = m(7) * (t142 * t366 + t461 * t93) + t366 * t223 + t461 * t302;
t635 = qJD(1) * t39;
t633 = t106 * t405;
t237 = t369 * t366;
t625 = t367 * t368;
t489 = t308 * t727 + t309 * t726 + t599 * t729 + (t219 + t220) * t708 + t779 * (-t237 / 0.2e1 - t625 / 0.2e1);
t77 = t369 * t115;
t471 = (-t463 ^ 2 * t698 + t106 * t369 - t77 - (-t107 + t114) * t367) * t760 + (-t167 * t463 - t369 * t94 - t77 - (t114 - t93) * t367) * t758 + t389 * t708 + t489;
t12 = t461 * t768 + t471 - t787;
t632 = t12 * qJD(1);
t486 = t516 * t759 + t524;
t579 = t752 + mrSges(7,2) / 0.2e1;
t14 = (t220 / 0.2e1 + t218 * t758 + t743) * t463 - (-t366 * t579 + t638 * t759 - t542) * t369 - (-t368 * t579 + t637 * t759 - t543) * t367 + t486;
t631 = t14 * qJD(1);
t24 = t461 * mrSges(4,1) + t463 * mrSges(4,2) + t616 + t608 + mrSges(3,3) - t599 * t406 + t598 * t405 + (m(4) + m(3)) * qJ(2) + m(6) * (-t107 * t406 - t633) + m(7) * (t94 * t405 - t93 * t406) + m(5) * (t317 * t462 + t318 * t460);
t630 = t24 * qJD(1);
t256 = t369 * t406;
t623 = t369 * t461;
t622 = t404 * t460;
t621 = t438 * t366;
t619 = t442 * t368;
t615 = t460 * t413;
t611 = t461 * t463;
t609 = t462 * t411;
t604 = t463 * t366;
t595 = 0.2e1 * t570;
t590 = qJD(3) * t461;
t117 = (t719 + t604 / 0.4e1 + t623 / 0.4e1) * t763;
t589 = t117 * qJD(1);
t583 = mrSges(7,2) * t621;
t582 = t703 / 0.2e1;
t575 = mrSges(7,2) * t748;
t568 = t651 / 0.2e1;
t567 = t650 / 0.2e1;
t564 = t706 * t368;
t551 = t406 * t708;
t548 = t461 * t722;
t546 = -t464 * t420 / 0.2e1;
t540 = t454 - t685;
t536 = -t215 - t629;
t533 = mrSges(6,3) * t584;
t532 = pkin(4) * t563;
t530 = -t584 / 0.2e1;
t526 = t658 / 0.2e1 + t659 / 0.2e1 - t661 / 0.2e1 + t662 / 0.2e1;
t522 = mrSges(5,3) * (-t459 / 0.2e1 - t458 / 0.2e1);
t521 = pkin(4) * mrSges(6,3) * t564;
t520 = t366 * t533;
t46 = m(5) * (-0.1e1 + t591) * t611 + (m(6) + m(7)) * (t237 - t611 + t625);
t466 = -(t304 / 0.2e1 + t734) * t369 + t543 * t366 - (-t306 / 0.2e1 + t733) * t367 - t542 * t368 + (-t222 / 0.2e1 - t221 / 0.2e1 + t390 / 0.2e1 + t412 * t709 + t414 * t713) * t463 + (t741 + t742 + t725 + t609 / 0.2e1 - t615 / 0.2e1) * t461 + ((-t317 * t460 + t318 * t462) * t463 + (t511 - 0.2e1 * t603) * t461) * t761 + (-t107 * t366 + t112 * t367 - t113 * t369 - t403 * t463 + t404 * t461 + t634) * t760 + (-t100 * t367 - t141 * t463 + t142 * t461 - t368 * t94 - t369 * t99 - t660) * t758;
t483 = t771 * (t299 * t405 - t406 * t769);
t9 = t466 - t483;
t515 = t9 * qJD(1) + t46 * qJD(2);
t510 = t438 * t367 - t442 * t369;
t509 = t100 * t759 + t643 / 0.2e1;
t508 = -t256 * t753 - t257 * t367 - t369 * t567 - t770;
t507 = m(7) * t517;
t505 = -t227 / 0.4e1 - t203 / 0.4e1 - t201 / 0.4e1 + t226 / 0.4e1;
t504 = -t229 / 0.4e1 - t228 / 0.4e1 - t197 / 0.4e1 + t744;
t502 = t738 - t276 / 0.4e1 - t279 / 0.4e1 - t281 / 0.4e1;
t500 = -t275 / 0.4e1 + t736 - t278 / 0.4e1 - t280 / 0.4e1;
t498 = t793 * t367 - t794 * t369;
t496 = -t367 * t705 - t369 * t706;
t16 = -pkin(3) * t420 + (t777 / 0.2e1 + t423 / 0.2e1) * t462 + (pkin(4) * t273 - t421 / 0.2e1 + t424 / 0.2e1) * t460 + m(6) * t443 * t699 + t797 + t790 * t246;
t473 = (-t619 - t621) * t758 + (-t366 * t705 + t564) * t588 - mrSges(5,1) * t613 / 0.2e1 + t607 * t754 + t526;
t477 = -t246 * t463 * t759 + t436 * t760 + t725;
t18 = t473 + t477 + t770;
t490 = -t106 * t769 - t107 * t299 + t539;
t2 = -t813 - (t229 + t228 + t197) * t406 / 0.4e1 - (t227 + t203 + t201) * t405 / 0.4e1 - (t280 + t278 + t275) * t366 / 0.4e1 + t788 - (t423 + t777) * t613 / 0.4e1 - t599 * t299 / 0.2e1 + ((t443 * t607 + t622) * pkin(4) + t490) * t760 + t366 * t736 - t368 * t738 + t587 * t739 + t406 * t744 + t633 * t752 + t607 * t715 + t607 * t716 + t226 * t719 + t462 * t724 + t540 * t710 + t107 * t567 + t93 * t568 + t94 * t570 + t610 / 0.4e1 - t617 / 0.4e1 + t463 * t546 + (t276 + t789) * t368 / 0.4e1 + t779 * (t114 * t718 + t115 * t722 + t775) + t768 * pkin(8) + t772 * t460;
t495 = t2 * qJD(1) - t18 * qJD(2) + t16 * qJD(3);
t17 = t790 * t269 + t797;
t474 = -t269 * t463 * t758 + t508;
t20 = -(mrSges(6,1) / 0.2e1 + t810) * t368 + (-mrSges(6,2) / 0.2e1 + mrSges(7,3) / 0.2e1) * t366 + t507 / 0.2e1 + t474;
t465 = t803 + t505 * t405 + t504 * t406 + ((-t107 / 0.2e1 + t751) * t406 + (t750 - t106 / 0.2e1) * t405 + t775) * mrSges(7,2) - (t299 * t753 + t502) * t368 + t500 * t366 + (t142 * t269 + t218 * t239 + t299 * t637) * t758 + t218 * t740 + t269 * t742 + t479 + (t752 * t366 + t638 * t758 + t542) * t769;
t485 = (-pkin(5) * t100 + qJ(6) * t99) * t759 + pkin(5) * t733 + t303 * t748;
t3 = t465 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t463 + t485 - t765 + t766;
t494 = t3 * qJD(1) + t20 * qJD(2) + t17 * qJD(3);
t217 = (t551 + t368 / 0.2e1) * m(7);
t480 = (t142 * t406 + t239 * t366 + t461 * t769) * t758 + t366 * t740 + t406 * t742;
t34 = (t548 + t726) * mrSges(7,2) + t480 + t509;
t65 = t790 * t406;
t492 = qJD(1) * t34 + qJD(2) * t217 + qJD(3) * t65;
t415 = m(7) * t438 + mrSges(7,3);
t481 = t647 + (t438 * t461 + t93) * t758;
t43 = -t704 / 0.2e1 + t481;
t491 = -qJD(1) * t43 - qJD(4) * t415 + t802;
t467 = (t438 * t106 + t442 * t107 + (t705 * t94 + t706 * t93) * pkin(4)) * t758 - t680 / 0.2e1 + t679 / 0.2e1 - t678 / 0.2e1 - t677 / 0.2e1 + t583 / 0.2e1 + t442 * t328 + t308 * t530 + t309 * t529 - t521 / 0.2e1 + t520 / 0.2e1 + t599 * t531;
t476 = (-pkin(5) * t114 + qJ(6) * t115) * t759 + t674 / 0.2e1 + t673 / 0.2e1 + t672 / 0.2e1 - t671 / 0.2e1 + pkin(5) * t328 + t366 * t575;
t13 = t467 + t476;
t484 = t794 * t584 + t793 * t585;
t166 = -(t438 * t706 + t442 * t705) * t755 - t484;
t468 = ((t584 - t438) * t299 + (t442 + t585) * t769) * t758 + t438 * t568 + t442 * t570 - mrSges(7,2) * t532 / 0.2e1 + t530 * t651 - t811;
t475 = pkin(5) * t570 + t406 * t575 + t798 * t759 + t811;
t23 = t468 + t475;
t482 = m(7) * (pkin(4) * t496 + t510);
t47 = t582 - t482 / 0.2e1;
t488 = t13 * qJD(1) - t47 * qJD(2) + t23 * qJD(3) - t166 * qJD(4);
t45 = t647 + (t565 / 0.4e1 + t560 / 0.4e1 + t453 / 0.2e1 - t107 / 0.4e1) * t763;
t487 = qJD(1) * t45 + qJD(4) * t452 + t785 + t802;
t382 = mrSges(7,3) + (0.2e1 * t529 + qJ(6)) * m(7);
t216 = m(7) * t551 - t368 * t758;
t124 = t595 + t702;
t118 = (-t604 - t623) * t758 + m(7) * t720;
t69 = t702 / 0.2e1 + t769 * t758 + t595;
t44 = (0.2e1 * t453 + t107) * t758 + t647 + t781 + m(7) * t747;
t42 = t781 + t704 / 0.2e1 + t481;
t38 = t482 / 0.2e1 + t582 + t498;
t33 = mrSges(7,2) * t548 + t656 / 0.2e1 + t480 - t509;
t22 = t468 - t475 + t523;
t21 = -t507 / 0.2e1 + t474 + t526;
t19 = t473 - t477 + t508;
t15 = (-t218 * t463 - t367 * t637 - t369 * t638) * t758 + t486 + t489;
t11 = (t463 * t522 + t776) * t461 + t471 + t787;
t10 = t467 - t476 + t525;
t8 = t466 + t483;
t4 = t465 - t485 + t764;
t1 = (t579 * t769 + t500) * t366 - (-t299 * t579 + t502) * t368 + (pkin(4) * t622 + t490) * t760 + t454 * t710 + t803 + ((t106 / 0.2e1 + t745) * mrSges(6,3) + (t750 + t745) * mrSges(7,2) + t505) * t405 + ((t747 + t746) * mrSges(6,3) + (t746 + t751) * mrSges(7,2) + t504) * t406 + (t414 * t757 + t724 + t363 / 0.4e1) * t462 + (t412 * t757 - t361 / 0.4e1 - t645 / 0.4e1 + t772) * t460 + (t546 + (-t423 / 0.4e1 - t777 / 0.4e1) * t460 + pkin(8) * t522 + (t715 + t716 + (m(6) * t714 + t739) * pkin(4)) * t462) * t463 + t788 + t813;
t25 = [qJD(2) * t24 + qJD(3) * t5 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t39, t8 * qJD(3) + t11 * qJD(4) + t15 * qJD(5) + t118 * qJD(6) + t630 + 0.2e1 * t771 * qJD(2) * (-t367 * t405 + t256) t642 + t8 * qJD(2) + t1 * qJD(4) + t4 * qJD(5) + t33 * qJD(6) + (-t462 * t423 / 0.2e1 + t421 * t712 - Ifges(4,5) + t639 * t464) * t590 + ((-pkin(3) * t450 + pkin(8) * t511) * t761 + (-t112 * t299 + t113 * t769 + t403 * t443) * t760 + (t100 * t299 + t141 * t239 + t769 * t99) * t758) * t762 + (-mrSges(4,2) * t603 + t275 * t729 - t367 * t737 + t362 * t712 + t196 * t720 + t198 * t722 + t360 * t709 + t112 * t650 - t113 * t405 * mrSges(6,3) - t100 * t651 - t99 * t653 - Ifges(4,6) * t463 + t443 * t222 + t403 * t273 + pkin(3) * t390 + t141 * t272 + t239 * t221 + (Ifges(5,5) * t460 + Ifges(5,6) * t462 - t807 * t405 - t808 * t406) * t707 + t789 * t727 + (t202 + t200) * t718 + (t303 + t304) * t769 + (-t306 + t307) * t299 + (-t615 + t609) * pkin(8) + t511 * mrSges(5,3)) * qJD(3), t641 + t11 * qJD(2) + t1 * qJD(3) + (-t521 + t520 - Ifges(5,5) * t613 - Ifges(5,6) * t607 + (-t114 * t706 + t115 * t705) * t756 + mrSges(7,2) * t619 + m(7) * (t114 * t442 + t115 * t438) + t583 + t671 - t672 - t674 - t673 - t317 * mrSges(5,2) - t318 * mrSges(5,1) + t525) * qJD(4) + t10 * qJD(5) + t42 * qJD(6), t640 + t15 * qJD(2) + t4 * qJD(3) + t10 * qJD(4) + (m(7) * (-pkin(5) * t107 + qJ(6) * t106) + t679 - t680 - t678 - t677 + t517 * mrSges(7,2) + t525) * qJD(5) + t44 * qJD(6), qJD(2) * t118 + qJD(3) * t33 + qJD(4) * t42 + qJD(5) * t44 + t635; qJD(3) * t9 + qJD(4) * t12 - qJD(5) * t14 - qJD(6) * t117 - t630, qJD(3) * t46, t19 * qJD(4) + t21 * qJD(5) + t216 * qJD(6) + (t272 + t273 + t639) * t590 + ((t239 * t461 + t536) * t758 + (t443 * t461 + t536) * t760 + (t591 * t697 - t700) * t761) * t762 + t515 + ((-mrSges(4,2) + t806) * t463 + t779 * (t366 * t405 + t368 * t406)) * qJD(3), t632 + t19 * qJD(3) + (m(7) * t510 - mrSges(5,1) * t612 + mrSges(5,2) * t614 - t496 * t756 + t498) * qJD(4) + t38 * qJD(5) - t636, -t631 + t21 * qJD(3) + t38 * qJD(4) + (t498 + t703) * qJD(5) - t636, -t589 + t216 * qJD(3) - 0.2e1 * (qJD(4) / 0.2e1 + qJD(5) / 0.2e1) * t701; -qJD(2) * t9 + qJD(4) * t2 + qJD(5) * t3 + qJD(6) * t34 - t642, -qJD(4) * t18 + qJD(5) * t20 + qJD(6) * t217 - t515, qJD(4) * t16 + qJD(5) * t17 + qJD(6) * t65 (m(7) * (-t299 * t438 + t442 * t769) + (-t299 * t705 - t706 * t769) * t756 + t406 * t533 + mrSges(6,3) * t532 + t540 - t519 * pkin(8) - t774 * mrSges(7,2) + t812) * qJD(4) + t22 * qJD(5) + t69 * qJD(6) + t495, t22 * qJD(4) + (m(7) * t798 + t516 * mrSges(7,2) + t812) * qJD(5) + t124 * qJD(6) + t494, qJD(4) * t69 + qJD(5) * t124 + t492; -qJD(2) * t12 - qJD(3) * t2 + qJD(5) * t13 + qJD(6) * t43 - t641, qJD(3) * t18 - qJD(5) * t47 - t632, qJD(5) * t23 - t495 - t801, -qJD(5) * t166 + qJD(6) * t415 ((-pkin(5) * t705 + qJ(6) * t706) * t755 + t484) * qJD(5) + t382 * qJD(6) + t488, qJD(5) * t382 - t491; qJD(2) * t14 - qJD(3) * t3 - qJD(4) * t13 + qJD(6) * t45 - t640, -qJD(3) * t20 + qJD(4) * t47 + t631, -qJD(4) * t23 - t494 + t801, -t488 + t784, t784, t487; qJD(2) * t117 - qJD(3) * t34 - qJD(4) * t43 - qJD(5) * t45 - t635, -qJD(3) * t217 + t589, -t492 + (qJD(4) - qJD(5)) * t257, t491 - t785, -t487, 0;];
Cq  = t25;
