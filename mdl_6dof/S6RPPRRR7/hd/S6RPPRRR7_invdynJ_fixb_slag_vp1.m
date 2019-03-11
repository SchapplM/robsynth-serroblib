% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR7_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:08
% EndTime: 2019-03-09 02:33:46
% DurationCPUTime: 31.19s
% Computational Cost: add. (33534->1149), mult. (34027->1460), div. (0->0), fcn. (30395->10), ass. (0->577)
t492 = sin(qJ(1));
t876 = pkin(1) * t492;
t494 = cos(qJ(1));
t875 = t494 * pkin(1);
t485 = qJD(4) + qJD(5);
t736 = t485 * t494;
t484 = pkin(10) + qJ(4);
t461 = qJ(5) + t484;
t441 = sin(t461);
t493 = cos(qJ(6));
t731 = t492 * t493;
t491 = sin(qJ(6));
t734 = t491 * t494;
t340 = t441 * t734 + t731;
t730 = t493 * t494;
t732 = t492 * t491;
t341 = t441 * t730 - t732;
t704 = rSges(7,1) * t341 - rSges(7,2) * t340;
t442 = cos(t461);
t739 = t442 * t494;
t190 = rSges(7,3) * t739 - t704;
t784 = rSges(7,2) * t491;
t787 = rSges(7,1) * t493;
t595 = -t784 + t787;
t782 = rSges(7,3) * t441;
t257 = t442 * t595 + t782;
t412 = t485 * t492;
t673 = qJD(6) * t494;
t631 = t442 * t673;
t312 = t412 + t631;
t793 = t442 * pkin(5);
t364 = pkin(9) * t441 + t793;
t407 = qJD(6) * t441 + qJD(1);
t874 = t190 * t407 - t257 * t312 - t364 * t412;
t338 = -t441 * t732 + t730;
t339 = t441 * t731 + t734;
t741 = t442 * t492;
t177 = Icges(7,5) * t339 + Icges(7,6) * t338 - Icges(7,3) * t741;
t179 = -Icges(7,5) * t341 + Icges(7,6) * t340 + Icges(7,3) * t739;
t320 = Icges(7,4) * t341;
t182 = Icges(7,2) * t340 + Icges(7,6) * t739 - t320;
t319 = Icges(7,4) * t340;
t184 = Icges(7,1) * t341 - Icges(7,5) * t739 - t319;
t568 = -t182 * t340 - t341 * t184;
t764 = Icges(7,4) * t339;
t180 = Icges(7,2) * t338 - Icges(7,6) * t741 + t764;
t318 = Icges(7,4) * t338;
t183 = Icges(7,1) * t339 - Icges(7,5) * t741 + t318;
t726 = t180 * t338 + t183 * t339;
t870 = t568 + t726 + (-t177 * t492 - t179 * t494) * t442;
t574 = Icges(7,5) * t493 - Icges(7,6) * t491;
t248 = Icges(7,3) * t441 + t442 * t574;
t762 = Icges(7,4) * t493;
t578 = -Icges(7,2) * t491 + t762;
t250 = Icges(7,6) * t441 + t442 * t578;
t763 = Icges(7,4) * t491;
t582 = Icges(7,1) * t493 - t763;
t252 = Icges(7,5) * t441 + t442 * t582;
t103 = t248 * t739 + t250 * t340 - t252 * t341;
t462 = qJD(5) * t494;
t463 = qJD(4) * t494;
t413 = t463 + t462;
t674 = qJD(6) * t492;
t313 = -t442 * t674 + t413;
t74 = t177 * t739 + t180 * t340 - t183 * t341;
t869 = -t103 * t407 - t313 * t74;
t567 = t182 * t491 + t184 * t493;
t81 = t179 * t441 - t442 * t567;
t394 = Icges(6,4) * t741;
t743 = t441 * t492;
t760 = Icges(6,5) * t494;
t270 = Icges(6,1) * t743 + t394 + t760;
t766 = Icges(6,4) * t441;
t359 = Icges(6,1) * t442 - t766;
t579 = Icges(6,2) * t442 + t766;
t856 = -t579 + t359;
t765 = Icges(6,4) * t442;
t583 = Icges(6,1) * t441 + t765;
t271 = -Icges(6,5) * t492 + t494 * t583;
t357 = -Icges(6,2) * t441 + t765;
t858 = t357 * t494 + t271;
t825 = qJD(1) * t856 - t412 * t858 + t413 * (-Icges(6,2) * t743 + t270 + t394);
t857 = t357 + t583;
t269 = -Icges(6,6) * t492 + t494 * t579;
t859 = -t359 * t494 + t269;
t268 = Icges(6,6) * t494 + t492 * t579;
t860 = -t359 * t492 + t268;
t826 = qJD(1) * t857 - t412 * t859 + t413 * t860;
t866 = t441 * t826 - t442 * t825;
t73 = -t179 * t741 + t182 * t338 - t184 * t339;
t452 = sin(t484);
t865 = rSges(5,2) * t452;
t411 = pkin(9) * t739;
t742 = t441 * t494;
t310 = pkin(5) * t742 - t411;
t862 = t190 - t310;
t861 = t257 + t364;
t655 = t442 * t412;
t679 = qJD(1) * t494;
t855 = t441 * t679 + t655;
t453 = cos(t484);
t737 = t453 * t492;
t420 = Icges(5,4) * t737;
t738 = t452 * t492;
t761 = Icges(5,5) * t494;
t281 = Icges(5,1) * t738 + t420 + t761;
t767 = Icges(5,4) * t453;
t584 = Icges(5,1) * t452 + t767;
t282 = -Icges(5,5) * t492 + t494 * t584;
t370 = -Icges(5,2) * t452 + t767;
t331 = t370 * t494;
t520 = t492 * (t282 + t331) - t494 * (-Icges(5,2) * t738 + t281 + t420);
t768 = Icges(5,4) * t452;
t580 = Icges(5,2) * t453 + t768;
t279 = Icges(5,6) * t494 + t492 * t580;
t280 = -Icges(5,6) * t492 + t494 * t580;
t372 = Icges(5,1) * t453 - t768;
t332 = t372 * t492;
t333 = t372 * t494;
t521 = t492 * (t280 - t333) - t494 * (t279 - t332);
t854 = -t521 * t452 + t520 * t453;
t698 = t370 + t584;
t699 = -t580 + t372;
t853 = (t452 * t698 - t453 * t699) * qJD(1);
t852 = -qJ(3) * qJD(1) ^ 2 + qJDD(3);
t305 = rSges(6,1) * t741 - rSges(6,2) * t743;
t598 = rSges(6,1) * t441 + rSges(6,2) * t442;
t479 = t494 * rSges(6,3);
t273 = rSges(6,1) * t743 + rSges(6,2) * t741 + t479;
t788 = rSges(6,1) * t442;
t361 = -rSges(6,2) * t441 + t788;
t798 = pkin(4) * t452;
t488 = sin(pkin(10));
t799 = pkin(3) * t488;
t387 = t798 + t799;
t365 = t492 * t387;
t735 = t488 * t492;
t437 = pkin(3) * t735;
t490 = -pkin(7) - qJ(3);
t483 = -pkin(8) + t490;
t229 = t365 - t437 + (-t483 + t490) * t494;
t728 = qJ(3) + t490;
t344 = -t494 * t728 + t437;
t417 = qJ(2) * t492 + t875;
t755 = qJ(3) * t494;
t621 = t417 + t755;
t604 = t344 + t621;
t549 = t229 + t604;
t634 = t453 * t463;
t397 = pkin(4) * t634;
t466 = qJD(2) * t494;
t611 = qJD(3) * t492 - t466;
t601 = -t397 + t611;
t89 = -t361 * t413 + (t273 + t549) * qJD(1) + t601;
t850 = (qJD(1) * t305 + t413 * t598) * t89;
t102 = -t248 * t741 + t250 * t338 + t252 * t339;
t849 = t102 * t407 + t312 * t73;
t800 = rSges(7,3) + pkin(9);
t848 = t441 * t800;
t562 = t269 * t442 + t271 * t441;
t847 = t494 * t562;
t532 = t574 * t441;
t564 = t250 * t491 - t252 * t493;
t569 = t180 * t491 - t183 * t493;
t497 = t312 * (-t248 * t494 + t567) + t313 * (t248 * t492 + t569) + t407 * (Icges(7,3) * t442 - t532 + t564);
t846 = t497 * t442;
t559 = t280 * t453 + t282 * t452;
t843 = t559 * t494;
t575 = Icges(6,5) * t441 + Icges(6,6) * t442;
t267 = -Icges(6,3) * t492 + t494 * t575;
t109 = -t267 * t494 - t269 * t741 - t271 * t743;
t557 = t357 * t442 + t359 * t441;
t355 = Icges(6,5) * t442 - Icges(6,6) * t441;
t748 = t355 * t494;
t122 = t492 * t557 + t748;
t842 = qJD(1) * t122 + t109 * t412;
t680 = qJD(1) * t492;
t638 = t488 * t679;
t693 = pkin(3) * t638 + t490 * t680;
t840 = qJD(1) * (qJ(3) * t680 + t693) + qJDD(1) * t344;
t489 = cos(pkin(10));
t785 = rSges(4,2) * t489;
t326 = rSges(4,1) * t735 + rSges(4,3) * t494 + t492 * t785;
t839 = t326 + t621;
t432 = t442 * pkin(9);
t796 = pkin(5) * t441;
t363 = t432 - t796;
t418 = -rSges(3,2) * t494 + rSges(3,3) * t492;
t594 = -rSges(7,1) * t491 - rSges(7,2) * t493;
t744 = t441 * t485;
t139 = -t595 * t744 + (rSges(7,3) * t485 + qJD(6) * t594) * t442;
t134 = t492 * t139;
t661 = t442 * t784;
t380 = t494 * t661;
t740 = t442 * t493;
t662 = rSges(7,1) * t740;
t227 = t380 + (-t662 - t782) * t494;
t431 = t442 * rSges(7,3);
t696 = t441 * t784 + t431;
t256 = -t441 * t787 + t696;
t291 = t363 * t485;
t311 = t364 * t494;
t838 = -qJD(1) * t311 + t227 * t407 - t312 * t256 - t412 * t363 - (-t190 * t442 - t257 * t742) * qJD(6) + t492 * t291 + t134 + t861 * t679;
t705 = rSges(7,1) * t339 + rSges(7,2) * t338;
t188 = -rSges(7,3) * t741 + t705;
t656 = t441 * t736;
t529 = -t442 * t680 - t656;
t654 = t442 * t736;
t191 = t529 * pkin(9) + (t441 * t680 - t654) * pkin(5);
t226 = rSges(7,3) * t743 + (-t661 + t662) * t492;
t309 = pkin(5) * t741 + pkin(9) * t743;
t632 = t441 * t673;
t633 = t441 * t674;
t553 = t494 * t407;
t610 = qJD(1) * t441 + qJD(6);
t836 = t492 * t610 - t654;
t162 = -t491 * t836 + t493 * t553;
t163 = t491 * t553 + t493 * t836;
t597 = rSges(7,1) * t163 + rSges(7,2) * t162;
t99 = rSges(7,3) * t529 + t597;
t837 = -t188 * t632 - t190 * t633 + t226 * t312 - t313 * t227 + t309 * t412 + t413 * t311 + (t191 + t99) * t494;
t150 = t280 * t452 - t282 * t453;
t196 = qJD(1) * t279 - qJD(4) * t331;
t198 = -qJD(4) * t333 + (t492 * t584 + t761) * qJD(1);
t576 = Icges(5,5) * t452 + Icges(5,6) * t453;
t278 = -Icges(5,3) * t492 + t494 * t576;
t683 = qJD(1) * t278;
t835 = qJD(4) * t150 + t196 * t453 + t198 * t452 + t683;
t350 = t580 * qJD(4);
t351 = t584 * qJD(4);
t368 = Icges(5,5) * t453 - Icges(5,6) * t452;
t555 = t370 * t452 - t372 * t453;
t834 = qJD(1) * t368 + qJD(4) * t555 + t350 * t453 + t351 * t452;
t677 = qJD(4) * t492;
t197 = qJD(1) * t280 + t370 * t677;
t199 = qJD(1) * t282 + qJD(4) * t332;
t560 = t279 * t452 - t281 * t453;
t277 = Icges(5,3) * t494 + t492 * t576;
t684 = qJD(1) * t277;
t833 = qJD(4) * t560 - t197 * t453 - t199 * t452 + t684;
t409 = pkin(5) * t743;
t308 = -pkin(9) * t741 + t409;
t63 = t188 * t407 - t257 * t313 - t364 * t413 + (t308 + t549) * qJD(1) + t601;
t676 = qJD(6) * t442;
t832 = (-qJD(1) * t309 - t188 * t676 - t226 * t407 + t256 * t313 + t257 * t633 + t363 * t413 + t680 * t861) * t63 - g(2) * t380;
t614 = t857 * t485;
t615 = t856 * t485;
t830 = qJD(1) * t355 + t441 * t614 - t442 * t615;
t617 = -qJD(1) * t271 + t485 * t860;
t619 = qJD(1) * t269 + t270 * t485 + t357 * t412;
t266 = Icges(6,3) * t494 + t492 * t575;
t686 = qJD(1) * t266;
t829 = t441 * t617 - t442 * t619 + t686;
t618 = (t492 * t583 + t760) * qJD(1) + t859 * t485;
t620 = -qJD(1) * t268 + t485 * t858;
t685 = qJD(1) * t267;
t828 = t441 * t618 - t442 * t620 + t685;
t577 = -Icges(7,2) * t493 - t763;
t506 = t312 * (Icges(7,2) * t341 - t184 + t319) + t313 * (-Icges(7,2) * t339 + t183 + t318) + t407 * (t442 * t577 + t252);
t581 = -Icges(7,1) * t491 - t762;
t827 = t312 * (-Icges(7,1) * t340 + t182 - t320) + t313 * (-Icges(7,1) * t338 + t180 + t764) + t407 * (-t442 * t581 + t250);
t486 = t492 ^ 2;
t72 = -t177 * t741 + t726;
t28 = t313 * t72 + t849;
t824 = -t28 / 0.2e1;
t668 = qJD(1) * qJD(4);
t385 = qJDD(4) * t492 + t494 * t668;
t292 = qJD(1) * t462 + qJDD(5) * t492 + t385;
t666 = qJDD(6) * t442;
t146 = qJD(6) * t529 + t494 * t666 + t292;
t823 = t146 / 0.2e1;
t609 = t485 * qJD(1);
t675 = qJD(6) * t485;
t457 = qJDD(4) * t494;
t689 = qJDD(5) * t494 + t457;
t147 = -qJD(1) * t631 + (t441 * t675 - t609 - t666) * t492 + t689;
t822 = t147 / 0.2e1;
t821 = t292 / 0.2e1;
t293 = -t492 * t609 + t689;
t820 = t293 / 0.2e1;
t295 = qJDD(6) * t441 + t442 * t675 + qJDD(1);
t819 = t295 / 0.2e1;
t818 = -t312 / 0.2e1;
t817 = t312 / 0.2e1;
t816 = -t313 / 0.2e1;
t815 = t313 / 0.2e1;
t814 = t385 / 0.2e1;
t386 = -t492 * t668 + t457;
t813 = t386 / 0.2e1;
t812 = -t407 / 0.2e1;
t811 = t407 / 0.2e1;
t810 = -t412 / 0.2e1;
t809 = t412 / 0.2e1;
t808 = -t413 / 0.2e1;
t807 = t413 / 0.2e1;
t806 = t492 / 0.2e1;
t805 = -t494 / 0.2e1;
t804 = t494 / 0.2e1;
t803 = rSges(3,2) - pkin(1);
t802 = -rSges(5,3) - pkin(1);
t801 = -rSges(6,3) - pkin(1);
t797 = pkin(4) * t453;
t794 = g(2) * t494;
t792 = -qJD(1) / 0.2e1;
t791 = qJD(1) / 0.2e1;
t790 = -pkin(1) - qJ(3);
t789 = rSges(5,1) * t452;
t783 = rSges(3,3) * t494;
t164 = -t407 * t731 + (-t494 * t610 - t655) * t491;
t165 = t610 * t730 + (-t407 * t491 + t485 * t740) * t492;
t640 = t442 * t679;
t657 = t441 * t412;
t530 = -t640 + t657;
t93 = Icges(7,5) * t165 + Icges(7,6) * t164 + Icges(7,3) * t530;
t95 = Icges(7,4) * t165 + Icges(7,2) * t164 + Icges(7,6) * t530;
t97 = Icges(7,1) * t165 + Icges(7,4) * t164 + Icges(7,5) * t530;
t22 = (t485 * t569 + t93) * t441 + (t177 * t485 - t491 * t95 + t493 * t97 + (-t180 * t493 - t183 * t491) * qJD(6)) * t442;
t781 = t22 * t313;
t92 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t529;
t94 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t529;
t96 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t529;
t23 = (t485 * t567 + t92) * t441 + (t179 * t485 - t491 * t94 + t493 * t96 + (-t182 * t493 + t184 * t491) * qJD(6)) * t442;
t780 = t23 * t312;
t160 = -t361 * t736 + (t492 * t598 + t479) * qJD(1);
t289 = t598 * t485;
t669 = qJD(1) * qJD(3);
t670 = qJD(1) * qJD(2);
t692 = qJDD(2) * t492 + t494 * t670;
t528 = -0.2e1 * t492 * t669 + t494 * t852 + t692;
t663 = qJD(4) ^ 2 * t798;
t501 = t385 * t797 - t492 * t663 + t528;
t274 = -rSges(6,3) * t492 + t494 * t598;
t664 = t494 * t799;
t700 = t387 * t494 + t483 * t492;
t228 = t490 * t492 + t664 - t700;
t343 = t492 * t728 + t664;
t468 = t494 * qJ(2);
t414 = -t468 + t876;
t756 = qJ(3) * t492;
t622 = -t414 - t756;
t605 = t343 + t622;
t550 = -t228 + t605;
t536 = t274 + t550;
t430 = t483 * t679;
t434 = t490 * t679;
t187 = -t397 - t430 + t434 + (t387 - t799) * t680;
t352 = qJD(1) * t417 - t466;
t706 = t434 - (t437 - t755) * qJD(1) - t352;
t652 = -t187 + t706;
t43 = -t289 * t412 + t292 * t361 + (-t160 + t652) * qJD(1) + t536 * qJDD(1) + t501;
t778 = t43 * t492;
t648 = rSges(6,1) * t855 + rSges(6,2) * t640;
t161 = (-rSges(6,2) * t744 - rSges(6,3) * qJD(1)) * t492 + t648;
t635 = t453 * t677;
t395 = pkin(4) * t635;
t649 = t387 * t679 + t483 * t680 + t395;
t186 = t649 - t693;
t448 = qJ(2) * t679;
t465 = qJD(2) * t492;
t691 = t448 + t465;
t650 = qJD(1) * (-pkin(1) * t680 + t691) + qJDD(1) * t417 + t492 * t670;
t508 = qJDD(1) * t755 + t492 * t852 + 0.2e1 * t494 * t669 + t650;
t667 = qJDD(2) * t494;
t505 = t508 - t667;
t499 = qJD(1) * t186 + qJDD(1) * t229 - t386 * t797 + t494 * t663 + t505 + t840;
t44 = qJD(1) * t161 + qJDD(1) * t273 + t289 * t413 - t293 * t361 + t499;
t777 = t44 * t494;
t480 = t494 * rSges(5,3);
t464 = qJD(3) * t494;
t688 = t464 + t465;
t645 = t395 + t688;
t746 = t412 * t361;
t88 = qJD(1) * t536 + t645 + t746;
t776 = t494 * t88;
t374 = rSges(5,1) * t453 - t865;
t336 = t374 * t494;
t599 = rSges(5,2) * t453 + t789;
t202 = -qJD(4) * t336 + (t492 * t599 + t480) * qJD(1);
t353 = t599 * qJD(4);
t476 = t492 * rSges(5,3);
t288 = t494 * t599 - t476;
t548 = t288 + t605;
t68 = -t353 * t677 + t374 * t385 + (-t202 + t706) * qJD(1) + t548 * qJDD(1) + t528;
t775 = t68 * t492;
t639 = t453 * t679;
t646 = rSges(5,1) * t635 + rSges(5,2) * t639 + t679 * t789;
t678 = qJD(4) * t452;
t203 = (-rSges(5,2) * t678 - rSges(5,3) * qJD(1)) * t492 + t646;
t287 = rSges(5,1) * t738 + rSges(5,2) * t737 + t480;
t69 = qJD(1) * t203 + qJDD(1) * t287 - t374 * t386 + (qJD(4) * t353 - qJDD(2)) * t494 + t508 + t840;
t774 = t69 * t494;
t80 = t177 * t441 - t442 * t569;
t773 = t80 * t147;
t772 = t81 * t146;
t771 = qJDD(1) / 0.2e1;
t104 = t248 * t441 - t442 * t564;
t573 = -Icges(7,5) * t491 - Icges(7,6) * t493;
t136 = -t485 * t532 + (Icges(7,3) * t485 + qJD(6) * t573) * t442;
t533 = t578 * t441;
t137 = -t485 * t533 + (Icges(7,6) * t485 + qJD(6) * t577) * t442;
t534 = t582 * t441;
t138 = -t485 * t534 + (Icges(7,5) * t485 + qJD(6) * t581) * t442;
t41 = (t485 * t564 + t136) * t441 + (-t137 * t491 + t138 * t493 + t248 * t485 + (-t250 * t493 - t252 * t491) * qJD(6)) * t442;
t770 = t104 * t295 + t407 * t41;
t602 = t228 * t463 - t229 * t677;
t58 = -t188 * t312 + t190 * t313 - t308 * t412 - t310 * t413 + t602;
t753 = qJD(1) * t58;
t342 = t374 * t677;
t119 = qJD(1) * t548 + t342 + t688;
t750 = t119 * t494;
t747 = t368 * t494;
t733 = t492 * t267;
t297 = t492 * t355;
t328 = t492 * t368;
t729 = -qJ(2) - t387;
t653 = rSges(7,1) * t165 + rSges(7,2) * t164 + rSges(7,3) * t657;
t100 = -rSges(7,3) * t640 + t653;
t647 = pkin(5) * t855 + pkin(9) * t657;
t192 = -pkin(9) * t640 + t647;
t727 = -t100 - t192;
t725 = t862 * t494;
t720 = -t188 - t308;
t719 = -t229 - t273;
t703 = t861 * t492;
t415 = rSges(3,2) * t492 + t783;
t695 = -t414 + t415;
t346 = t417 + t418;
t694 = rSges(4,1) * t638 + t679 * t785;
t690 = rSges(3,2) * t680 + rSges(3,3) * t679;
t687 = -qJD(1) * t414 + t465;
t681 = qJD(1) * t576;
t123 = t494 * t557 - t297;
t672 = t123 * qJD(1);
t556 = t370 * t453 + t372 * t452;
t141 = t494 * t556 - t328;
t671 = t141 * qJD(1);
t665 = -rSges(4,3) + t790;
t425 = pkin(4) * t737;
t659 = pkin(4) * t678;
t658 = -m(4) - m(5) - m(6) - m(7);
t651 = -t229 + t720;
t108 = t266 * t494 + t268 * t741 + t270 * t743;
t115 = t277 * t494 + t279 * t737 + t281 * t738;
t116 = -t278 * t494 - t280 * t737 - t282 * t738;
t644 = t448 + t688;
t643 = t464 + t687;
t642 = t442 * t800;
t637 = t452 * t677;
t630 = -pkin(5) - t787;
t629 = -t680 / 0.2e1;
t628 = t679 / 0.2e1;
t627 = -t677 / 0.2e1;
t626 = t677 / 0.2e1;
t625 = -t463 / 0.2e1;
t624 = t463 / 0.2e1;
t623 = t63 * (-t139 - t291);
t406 = rSges(6,2) * t742;
t307 = -rSges(6,1) * t739 + t406;
t616 = -t305 * t412 + t307 * t413;
t612 = -qJD(1) * t307 - t412 * t598;
t607 = qJD(1) * t343 + t643;
t600 = rSges(4,1) * t488 + t785;
t327 = -rSges(4,3) * t492 + t494 * t600;
t606 = t327 + t622;
t603 = -t431 + t796;
t419 = rSges(2,1) * t494 - rSges(2,2) * t492;
t416 = rSges(2,1) * t492 + rSges(2,2) * t494;
t535 = t310 + t550;
t20 = t312 * t139 + t146 * t257 - t295 * t190 + t412 * t291 + t292 * t364 - t407 * t99 + (-t191 + t652) * qJD(1) + t535 * qJDD(1) + t501;
t21 = qJD(1) * t192 + qJDD(1) * t308 + t100 * t407 - t139 * t313 - t147 * t257 + t188 * t295 - t291 * t413 - t293 * t364 + t499;
t593 = t20 * t494 + t21 * t492;
t561 = t279 * t453 + t281 * t452;
t514 = qJD(1) * t561 + qJD(4) * t328 + t683;
t515 = -qJD(1) * t559 - qJD(4) * t747 + t684;
t592 = t492 * (t492 * t515 - t494 * t835) + t494 * (t492 * t514 + t494 * t833);
t591 = t492 * (t492 * t835 + t494 * t515) + t494 * (-t492 * t833 + t494 * t514);
t590 = t492 * t73 + t494 * t72;
t589 = t492 * t72 - t494 * t73;
t75 = t179 * t739 - t568;
t588 = t492 * t75 + t494 * t74;
t587 = t492 * t74 - t494 * t75;
t586 = t492 * t81 + t494 * t80;
t585 = t492 * t80 - t494 * t81;
t572 = t115 * t494 + t116 * t492;
t261 = t492 * t277;
t117 = -t494 * t561 + t261;
t118 = -t278 * t492 + t843;
t571 = t117 * t494 + t118 * t492;
t120 = -t374 * t463 + (t287 + t604) * qJD(1) + t611;
t570 = t119 * t492 - t120 * t494;
t566 = t188 * t494 + t190 * t492;
t565 = t202 * t494 - t203 * t492;
t563 = t268 * t442 + t270 * t441;
t133 = t269 * t441 - t271 * t442;
t558 = -t287 * t492 - t288 * t494;
t554 = -t414 + t700;
t552 = (-t494 ^ 2 - t486) * qJD(4) * t797;
t400 = pkin(4) * t639;
t551 = -pkin(4) * t637 + t400;
t547 = t108 + t733;
t546 = -t483 * t494 + t365 + t417;
t545 = t430 - t601;
t544 = -qJD(1) * t228 + t395 + t607;
t543 = t644 + t649;
t531 = t599 + t799;
t527 = t177 * t313 + t179 * t312 + t248 * t407;
t526 = (Icges(7,5) * t338 - Icges(7,6) * t339) * t313 + (Icges(7,5) * t340 + Icges(7,6) * t341) * t312 + t573 * t442 * t407;
t523 = -qJD(1) * t575 + t297 * t413 - t412 * t748;
t253 = t492 * t266;
t110 = -t494 * t563 + t253;
t519 = t441 * t630 + t432 + t696;
t518 = -t186 * t677 + t187 * t463 + t228 * t386 - t229 * t385;
t517 = -qJD(1) * t562 - t485 * t748 + t686;
t516 = qJD(1) * t563 + t297 * t485 + t685;
t513 = qJD(1) * t557 - t485 * t575;
t512 = t556 * qJD(1) - qJD(4) * t576;
t511 = t226 + t309;
t510 = t442 * t630 - t848;
t111 = -t733 + t847;
t132 = -t268 * t441 + t270 * t442;
t16 = t162 * t180 + t163 * t183 + t177 * t529 + t340 * t95 - t341 * t97 + t739 * t93;
t17 = t162 * t182 - t163 * t184 + t179 * t529 + t340 * t94 - t341 * t96 + t739 * t92;
t18 = t164 * t180 + t165 * t183 + t177 * t530 + t338 * t95 + t339 * t97 - t741 * t93;
t19 = t164 * t182 - t165 * t184 + t179 * t530 + t338 * t94 + t339 * t96 - t741 * t92;
t222 = t250 * t492;
t223 = t250 * t494;
t224 = t252 * t492;
t225 = t252 * t494;
t249 = Icges(7,6) * t442 - t533;
t251 = Icges(7,5) * t442 - t534;
t29 = t312 * t75 - t869;
t33 = t136 * t739 + t137 * t340 - t138 * t341 + t162 * t250 + t163 * t252 + t248 * t529;
t3 = t103 * t295 + t146 * t75 + t147 * t74 + t16 * t313 + t17 * t312 + t33 * t407;
t32 = t104 * t407 + t312 * t81 + t313 * t80;
t34 = -t136 * t741 + t137 * t338 + t138 * t339 + t164 * t250 + t165 * t252 + t248 * t530;
t4 = t102 * t295 + t146 * t73 + t147 * t72 + t18 * t313 + t19 * t312 + t34 * t407;
t45 = t492 * t516 + t494 * t829;
t46 = t492 * t517 - t494 * t828;
t47 = -t492 * t829 + t494 * t516;
t48 = t492 * t828 + t494 * t517;
t56 = t108 * t413 + t842;
t57 = t110 * t413 + t111 * t412 - t672;
t70 = t492 * t513 + t494 * t830;
t71 = -t492 * t830 + t494 * t513;
t76 = -t441 * t619 - t442 * t617;
t77 = t441 * t620 + t442 * t618;
t498 = (t523 * t492 + t494 * t866) * t810 + (-t492 * t866 + t523 * t494) * t808 + t57 * t628 + (qJD(1) * t71 + qJDD(1) * t122 + t108 * t293 + t109 * t292 + t412 * t48 + t413 * t47 + t4) * t804 + (-t441 * t825 - t442 * t826) * t792 + (t628 + t632 / 0.2e1) * t29 + (qJD(1) * t70 - qJDD(1) * t123 + t110 * t293 + t111 * t292 + t412 * t46 + t413 * t45 + t3) * t806 + (t28 + t56) * t629 + t586 * t819 + (t108 * t494 + t109 * t492) * t820 + (t110 * t494 + t111 * t492) * t821 + t590 * t822 + t588 * t823 + t633 * t824 + (t47 * t494 + t48 * t492 + (-t108 * t492 + t109 * t494) * qJD(1)) * t807 + (t45 * t494 + t46 * t492 + (-t110 * t492 + t111 * t494) * qJD(1)) * t809 + (-qJD(1) * t585 + t22 * t494 + t23 * t492) * t811 + (-qJD(1) * t589 + t18 * t494 + t19 * t492) * t815 + (-qJD(1) * t587 + t16 * t494 + t17 * t492) * t817 + ((t222 * t338 + t224 * t339) * t313 + (-t223 * t338 - t225 * t339) * t312 + (t249 * t338 + t251 * t339) * t407 + (t102 * t442 - t73 * t742) * qJD(6) + ((qJD(6) * t72 + t527) * t441 - t846) * t492) * t816 + ((t222 * t340 - t224 * t341) * t313 + (-t223 * t340 + t225 * t341) * t312 + (t249 * t340 - t251 * t341) * t407 + (t103 * t442 + t74 * t743) * qJD(6) + ((-qJD(6) * t75 - t527) * t441 + t846) * t494) * t818 + (((-t222 * t491 + t224 * t493 + t177) * t313 + (t223 * t491 - t225 * t493 + t179) * t312 + (-t249 * t491 + t251 * t493 + t248) * t407 + t104 * qJD(6)) * t442 + (qJD(6) * t585 + t497) * t441) * t812 + (t132 * t494 + t133 * t492) * t771 + (t492 * t77 + t494 * t76 + (-t132 * t492 + t133 * t494) * qJD(1)) * t791 - t32 * t676 / 0.2e1;
t335 = t374 * t492;
t306 = t594 * t442;
t260 = qJD(1) * t346 - t466;
t259 = qJD(1) * t695 + t465;
t258 = t494 * t274;
t217 = t494 * t228;
t215 = rSges(7,1) * t340 + rSges(7,2) * t341;
t214 = rSges(7,1) * t338 - rSges(7,2) * t339;
t201 = qJD(1) * t839 + t611;
t200 = qJD(1) * t606 + t688;
t170 = t494 * t187;
t166 = t558 * qJD(4);
t145 = t494 * t160;
t140 = t492 * t556 + t747;
t135 = t140 * qJD(1);
t131 = qJD(1) * t690 + qJDD(1) * t418 + t650 - t667;
t130 = t695 * qJDD(1) + (-qJD(1) * t418 - t352) * qJD(1) + t692;
t106 = qJDD(1) * t326 + qJD(1) * (-rSges(4,3) * t680 + t694) + t505;
t105 = t606 * qJDD(1) + (-qJD(1) * t326 - t352) * qJD(1) + t528;
t87 = -t273 * t412 - t274 * t413 + t602;
t83 = qJD(4) * t559 - t196 * t452 + t198 * t453;
t82 = -qJD(4) * t561 - t197 * t452 + t199 * t453;
t79 = -t492 * t834 + t494 * t512;
t78 = t492 * t512 + t494 * t834;
t65 = qJD(4) * t571 - t671;
t64 = qJD(4) * t572 + t135;
t62 = qJD(1) * t535 + t645 - t874;
t39 = t160 * t413 - t161 * t412 - t273 * t292 - t274 * t293 + t518;
t15 = -t100 * t312 - t146 * t188 + t147 * t190 + t191 * t413 - t192 * t412 - t292 * t308 - t293 * t310 + t313 * t99 + t518;
t1 = [(t140 - t560) * t813 + ((t111 + t547 - t847) * t413 + t842) * t810 + (t83 + t78 + t64) * t626 + (t82 + t79) * t624 + (t135 + ((-t117 + t261 + t116) * t492 + (t118 - t843 + (t278 - t561) * t492 + t115) * t494) * qJD(4)) * t627 + ((t21 - g(2)) * (-t492 * t642 + t409 + t546 + t705) + (t545 - t597 + (t793 + t848) * t736 + (-t875 + (-t603 + t729 + t432) * t492) * qJD(1)) * t62 + (t543 + t647 + t653 + t62 - t544 + (-t494 * t642 - t310 + t756 - t876) * qJD(1) + t874) * t63 + (-g(1) + t20) * (t494 * t603 - t411 + t554 + t704)) * m(7) - t292 * t123 / 0.2e1 + (t33 + t28) * t817 + (-qJD(4) * t556 + t350 * t452 - t351 * t453 - t441 * t615 - t442 * t614) * qJD(1) + (-(qJD(1) * t415 - t259 + t687) * t260 + t259 * t466 + t260 * (t690 + t691) + (t259 * t803 * t494 + (t259 * (-rSges(3,3) - qJ(2)) - t260 * pkin(1)) * t492) * qJD(1) + (t131 - g(2)) * t346 + (t130 - g(1)) * (t492 * t803 + t468 + t783)) * m(3) + t780 / 0.2e1 + t781 / 0.2e1 + (t88 * (rSges(6,1) * t654 - rSges(6,2) * t656 + t545) + t89 * (-rSges(6,2) * t657 + t543 + t648) + (t801 * t776 + (t88 * (-t598 + t729) + t89 * t801) * t492) * qJD(1) - (t746 - t88 + (t274 - t756) * qJD(1) + t544) * t89 + (t44 - g(2)) * (t546 + t273) + (t43 - g(1)) * (t274 + t554)) * m(6) + (t122 + t132) * t820 + (t65 + t671 + (t486 * t278 + (-t261 + t116 + (t278 + t561) * t494) * t494) * qJD(4)) * t625 + (-t555 + Icges(4,1) * t489 ^ 2 + (-0.2e1 * Icges(4,4) * t489 + Icges(4,2) * t488) * t488 + m(2) * (t416 ^ 2 + t419 ^ 2) + Icges(3,1) + Icges(2,3) - t357 * t441 + t359 * t442) * qJDD(1) + (-(-t200 + (t327 - t756) * qJD(1) + t643) * t201 - t200 * t611 + t201 * (t644 + t694) + (t200 * t665 * t494 + (t200 * (-qJ(2) - t600) + t201 * t665) * t492) * qJD(1) + (-g(2) + t106) * t839 + (-g(1) + t105) * (t492 * t790 + t327 + t468)) * m(4) + t773 / 0.2e1 + t772 / 0.2e1 + t133 * t821 + t102 * t822 + t103 * t823 + t150 * t814 + t34 * t815 + t770 + (t76 + t71) * t807 + ((-t72 + t870) * t312 + t29 + t869) * t816 + ((t75 + t870) * t313 + t849) * t818 - t385 * t141 / 0.2e1 - m(2) * (-g(1) * t416 + g(2) * t419) + (t672 + (t492 * t562 + t109 - t253) * t413 + (-t108 + t547) * t412 + ((t267 + t563) * t413 - t562 * t412) * t494 + t57) * t808 + (t119 * (rSges(5,1) * t634 - t463 * t865 + t434 - t611) + t120 * (-rSges(5,2) * t637 + t644 + t646 + t693) + (t802 * t750 + (t119 * (-qJ(2) - t531) + t120 * t802) * t492) * qJD(1) - (-t119 + t342 + (t288 - t756) * qJD(1) + t607) * t120 + (t69 - g(2)) * (-t490 * t494 + t287 + t417 + t437) + (t68 - g(1)) * (t468 - t476 + (-pkin(1) + t490) * t492 + t531 * t494)) * m(5) + (t77 + t70 + t56) * t809; (-m(3) + t658) * (g(1) * t492 - t794) + 0.2e1 * (t20 * t806 + t21 * t805) * m(7) + 0.2e1 * (t778 / 0.2e1 - t777 / 0.2e1) * m(6) + 0.2e1 * (t775 / 0.2e1 - t774 / 0.2e1) * m(5) + 0.2e1 * (t105 * t806 + t106 * t805) * m(4) + 0.2e1 * (t130 * t806 + t131 * t805) * m(3); t658 * (g(1) * t494 + g(2) * t492) + m(4) * (t105 * t494 + t106 * t492) + m(5) * (t492 * t69 + t494 * t68) + m(6) * (t43 * t494 + t44 * t492) + m(7) * t593; ((-t117 * t492 + t118 * t494) * qJD(1) + t592) * t626 + ((-t115 * t492 + t116 * t494) * qJD(1) + t591) * t624 + ((-t677 * t747 - t681) * t492 + (t853 + (t492 * t328 + t854) * qJD(4)) * t494) * t627 + t65 * t628 + (qJD(1) * t78 + qJD(4) * t592 - qJDD(1) * t141 + t117 * t386 + t118 * t385) * t806 + t64 * t629 + ((t328 * t463 - t681) * t494 + (-t853 + (-t494 * t747 - t854) * qJD(4)) * t492) * t625 + t498 + (t492 * t83 + t494 * t82 + (t150 * t494 + t492 * t560) * qJD(1)) * t791 + ((-t452 * t699 - t453 * t698) * qJD(1) + (t452 * t520 + t453 * t521) * qJD(4)) * t792 + (qJD(1) * t79 + qJD(4) * t591 + qJDD(1) * t140 + t115 * t386 + t116 * t385) * t804 + (t150 * t492 - t494 * t560) * t771 + t571 * t814 + t572 * t813 + (-g(1) * (t425 + t511) - g(3) * (t519 - t798) - (t510 - t797) * t794 + t20 * (t425 + t703) + t15 * (t217 + t725) + (t21 * (-t861 - t797) + t623 + t651 * t753) * t494 + (t15 * t651 + (-t228 - t862) * t753) * t492 + (-t492 * t659 + t400 - t551 + t838) * t62 + (-t552 + t170 + (-t186 + t727) * t492 + t837) * t58 + t832) * m(7) + (t43 * t425 + t88 * t400 + t39 * (t217 - t258) + t87 * (t145 + t170) + (t44 * (-t361 - t797) + t89 * t289 + (t361 * t88 + t719 * t87) * qJD(1)) * t494 + (t43 * t361 + t88 * (-t289 - t659) + t39 * t719 + t87 * (-t161 - t186) + (t89 * t361 + t87 * (-t228 + t274)) * qJD(1)) * t492 - t88 * (t551 + t612) - t850 - t87 * (t552 + t616) - g(1) * (t305 + t425) - g(2) * (t406 + (-t788 - t797) * t494) - g(3) * (-t598 - t798)) * m(6) + (-(t119 * t336 + t120 * t335) * qJD(1) - (t166 * (-t335 * t492 - t336 * t494) - t570 * t599) * qJD(4) + (qJD(4) * t565 - t287 * t385 - t288 * t386) * t558 + t166 * ((-t287 * t494 + t288 * t492) * qJD(1) + t565) - t570 * t353 + (t775 - t774 + (t120 * t492 + t750) * qJD(1)) * t374 - g(1) * t335 + g(2) * t336 + g(3) * t599) * m(5); t498 + (t20 * t703 + (-t21 * t861 + t720 * t753 + t623) * t494 - g(1) * t511 - g(3) * t519 - t510 * t794 + t838 * t62 + (t492 * t720 + t725) * t15 + ((-qJD(1) * t862 + t727) * t492 + t837) * t58 + t832) * m(7) + (t39 * (-t273 * t492 - t258) - (t492 * t88 - t494 * t89) * t289 + (t778 - t777 + (t492 * t89 + t776) * qJD(1)) * t361 - t612 * t88 - t850 - g(1) * t305 - g(2) * t307 + g(3) * t598 + (-t161 * t492 + t145 + (-t273 * t494 + t274 * t492) * qJD(1) - t616) * t87) * m(6); t640 * t824 + t28 * t657 / 0.2e1 - t4 * t741 / 0.2e1 + (t102 * t441 - t442 * t589) * t822 + ((t485 * t589 + t34) * t441 + (-qJD(1) * t590 + t102 * t485 - t18 * t492 + t19 * t494) * t442) * t815 + t3 * t739 / 0.2e1 + (t103 * t441 - t442 * t587) * t823 + ((t485 * t587 + t33) * t441 + (-qJD(1) * t588 + t103 * t485 - t16 * t492 + t17 * t494) * t442) * t817 + t485 * t442 * t32 / 0.2e1 + t441 * (t770 + t772 + t773 + t780 + t781) / 0.2e1 + (t104 * t441 - t442 * t585) * t819 + ((t485 * t585 + t41) * t441 + (-qJD(1) * t586 + t104 * t485 - t22 * t492 + t23 * t494) * t442) * t811 + (t338 * t506 - t339 * t827 - t526 * t741) * t816 + (t506 * t340 + t341 * t827 + t526 * t739) * t818 + (t526 * t441 + (-t491 * t506 - t493 * t827) * t442) * t812 + (t442 * t629 - t656 / 0.2e1) * t29 + ((t63 * t100 + t21 * t188 - t20 * t190 - t62 * t99 + (t58 * t566 + (-t492 * t63 - t494 * t62) * t257) * t485) * t441 + (t62 * (t139 * t494 - t190 * t485) + t63 * (t188 * t485 + t134) - t15 * t566 + t58 * (-t100 * t494 + t188 * t680 - t190 * t679 - t492 * t99) + ((-t492 * t62 + t494 * t63) * qJD(1) + t593) * t257) * t442 - t62 * (-t215 * t407 + t306 * t312) - t63 * (t214 * t407 - t306 * t313) - t58 * (-t214 * t312 + t215 * t313) - g(1) * t214 - g(2) * t215 - g(3) * t306) * m(7);];
tau  = t1;
