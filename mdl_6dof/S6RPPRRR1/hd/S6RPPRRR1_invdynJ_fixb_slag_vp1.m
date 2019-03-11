% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR1_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:57
% EndTime: 2019-03-09 02:18:37
% DurationCPUTime: 34.71s
% Computational Cost: add. (46333->1109), mult. (32876->1426), div. (0->0), fcn. (29609->12), ass. (0->561)
t478 = pkin(11) + qJ(4);
t472 = cos(t478);
t463 = pkin(4) * t472;
t482 = cos(pkin(11));
t466 = t482 * pkin(3) + pkin(2);
t400 = t463 + t466;
t480 = qJ(1) + pkin(10);
t471 = sin(t480);
t473 = cos(t480);
t483 = -pkin(7) - qJ(3);
t477 = -pkin(8) + t483;
t698 = -t471 * t400 - t473 * t477;
t485 = sin(qJ(1));
t802 = pkin(1) * t485;
t867 = t698 - t802;
t487 = cos(qJ(1));
t476 = t487 * pkin(1);
t383 = rSges(3,1) * t471 + rSges(3,2) * t473;
t346 = -t383 - t802;
t474 = qJ(5) + t478;
t465 = cos(t474);
t734 = t465 * t471;
t662 = rSges(6,1) * t734;
t866 = -t662 + t867;
t722 = t473 * t483;
t865 = t466 * t471 + t722;
t486 = cos(qJ(6));
t720 = t473 * t486;
t484 = sin(qJ(6));
t726 = t471 * t484;
t332 = t465 * t726 + t720;
t721 = t473 * t484;
t725 = t471 * t486;
t333 = t465 * t725 - t721;
t594 = t333 * rSges(7,1) - t332 * rSges(7,2);
t464 = sin(t474);
t738 = t464 * t471;
t192 = -rSges(7,3) * t738 - t594;
t791 = rSges(7,1) * t486;
t593 = -rSges(7,2) * t484 + t791;
t291 = -rSges(7,3) * t465 + t464 * t593;
t479 = qJD(4) + qJD(5);
t373 = t473 * t479;
t674 = qJD(6) * t471;
t296 = -t464 * t674 + t373;
t800 = pkin(5) * t464;
t366 = -pkin(9) * t465 + t800;
t675 = qJD(6) * t465;
t416 = qJD(1) - t675;
t444 = qJD(3) * t471;
t470 = sin(t478);
t677 = qJD(4) * t473;
t640 = t470 * t677;
t610 = pkin(4) * t640;
t551 = t444 - t610;
t864 = t416 * t192 - t296 * t291 - t373 * t366 + t551;
t181 = Icges(7,5) * t333 - Icges(7,6) * t332 + Icges(7,3) * t738;
t313 = Icges(7,4) * t333;
t184 = -Icges(7,2) * t332 + Icges(7,6) * t738 + t313;
t312 = Icges(7,4) * t332;
t188 = -Icges(7,1) * t333 - Icges(7,5) * t738 + t312;
t861 = t184 * t484 + t188 * t486;
t78 = -t181 * t465 - t464 * t861;
t489 = qJD(1) ^ 2;
t863 = t489 * t476;
t573 = Icges(7,5) * t486 - Icges(7,6) * t484;
t285 = -Icges(7,3) * t465 + t464 * t573;
t764 = Icges(7,4) * t486;
t575 = -Icges(7,2) * t484 + t764;
t287 = -Icges(7,6) * t465 + t464 * t575;
t765 = Icges(7,4) * t484;
t579 = Icges(7,1) * t486 - t765;
t289 = -Icges(7,5) * t465 + t464 * t579;
t101 = t285 * t738 - t287 * t332 + t289 * t333;
t443 = qJD(4) * t471;
t372 = qJD(5) * t471 + t443;
t673 = qJD(6) * t473;
t295 = t464 * t673 + t372;
t69 = t181 * t738 - t184 * t332 - t188 * t333;
t334 = -t465 * t721 + t725;
t335 = t465 * t720 + t726;
t737 = t464 * t473;
t183 = Icges(7,5) * t335 + Icges(7,6) * t334 + Icges(7,3) * t737;
t766 = Icges(7,4) * t335;
t186 = Icges(7,2) * t334 + Icges(7,6) * t737 + t766;
t314 = Icges(7,4) * t334;
t189 = Icges(7,1) * t335 + Icges(7,5) * t737 + t314;
t70 = t183 * t738 - t332 * t186 + t333 * t189;
t28 = t101 * t416 + t295 * t70 - t296 * t69;
t102 = t285 * t737 + t287 * t334 + t289 * t335;
t71 = t181 * t737 + t334 * t184 - t188 * t335;
t72 = t183 * t737 + t334 * t186 + t335 * t189;
t29 = t102 * t416 + t295 * t72 - t296 * t71;
t818 = -t295 / 0.2e1;
t815 = t296 / 0.2e1;
t860 = rSges(5,2) * t470;
t665 = -qJDD(4) - qJDD(5);
t666 = qJDD(6) * t464;
t672 = qJD(6) * t479;
t676 = qJD(6) * t464;
t668 = qJD(1) * qJD(4);
t425 = t471 * t668;
t667 = qJD(1) * qJD(5);
t685 = t471 * t667 + t425;
t134 = (t465 * t672 + t666) * t471 + (qJD(1) * t676 + t665) * t473 + t685;
t592 = -rSges(7,1) * t484 - rSges(7,2) * t486;
t732 = t465 * t479;
t190 = t593 * t732 + (rSges(7,3) * t479 + qJD(6) * t592) * t464;
t727 = t471 * t479;
t652 = t465 * t727;
t678 = qJD(1) * t473;
t526 = t464 * t678 + t652;
t654 = t464 * t727;
t201 = t526 * pkin(9) + (t465 * t678 - t654) * pkin(5);
t272 = t473 * t665 + t685;
t438 = t465 * pkin(5);
t838 = t464 * pkin(9) + t438;
t336 = t838 * t479;
t339 = -qJDD(6) * t465 + t464 * t672 + qJDD(1);
t363 = -qJDD(4) * t473 + t425;
t669 = qJD(1) * qJD(3);
t555 = qJDD(3) * t471 + t473 * t669 - t863;
t723 = t472 * qJD(4) ^ 2;
t801 = pkin(4) * t470;
t509 = -pkin(4) * t473 * t723 + t363 * t801 + t555;
t308 = t838 * t471;
t227 = t698 + t865;
t447 = t473 * qJ(3);
t794 = pkin(2) - t466;
t538 = t471 * t794 - t722;
t260 = -t447 + t538;
t382 = pkin(2) * t471 - t447;
t625 = -t382 - t802;
t606 = t260 + t625;
t554 = t227 + t606;
t540 = -t308 + t554;
t785 = pkin(4) * qJD(4);
t657 = t470 * t785;
t393 = t471 * t657;
t679 = qJD(1) * t471;
t414 = t477 * t679;
t417 = t483 * t679;
t689 = t400 - t466;
t170 = t678 * t689 - t393 - t414 + t417;
t446 = t471 * qJ(3);
t385 = t473 * pkin(2) + t446;
t445 = qJD(3) * t473;
t325 = qJD(1) * t385 - t445;
t708 = t417 - (-t473 * t794 - t446) * qJD(1) - t325;
t647 = -t170 + t708;
t736 = t464 * t484;
t517 = t416 * t486 + t479 * t736;
t609 = qJD(1) * t465 - qJD(6);
t164 = t471 * t517 - t609 * t721;
t735 = t464 * t486;
t516 = t416 * t484 - t479 * t735;
t165 = t471 * t516 + t609 * t720;
t595 = rSges(7,1) * t165 + rSges(7,2) * t164;
t99 = rSges(7,3) * t526 + t595;
t21 = t134 * t291 - t296 * t190 + t339 * t192 + t272 * t366 - t373 * t336 - t416 * t99 + (-t201 + t647) * qJD(1) + t540 * qJDD(1) + t509;
t193 = t335 * rSges(7,1) + t334 * rSges(7,2) + rSges(7,3) * t737;
t733 = t465 * t473;
t310 = pkin(5) * t733 + pkin(9) * t737;
t349 = t473 * t400;
t410 = t473 * t466;
t683 = -t477 + t483;
t228 = t471 * t683 + t349 - t410;
t612 = -t471 * t483 + t410;
t261 = t612 - t385;
t624 = t385 + t476;
t605 = t261 + t624;
t553 = t228 + t605;
t690 = -t393 - t445;
t64 = t193 * t416 - t291 * t295 - t366 * t372 + (t310 + t553) * qJD(1) + t690;
t856 = t473 * (qJD(1) * t64 + t21);
t855 = -t192 + t308;
t854 = t193 + t310;
t757 = Icges(6,6) * t473;
t264 = Icges(6,4) * t734 - Icges(6,2) * t738 - t757;
t434 = Icges(6,4) * t465;
t358 = Icges(6,1) * t464 + t434;
t853 = -t358 * t471 - t264;
t576 = -Icges(6,2) * t464 + t434;
t265 = Icges(6,6) * t471 + t473 * t576;
t852 = -t358 * t473 - t265;
t767 = Icges(6,4) * t464;
t359 = Icges(6,1) * t465 - t767;
t267 = Icges(6,5) * t471 + t359 * t473;
t356 = Icges(6,2) * t465 + t767;
t851 = -t356 * t473 + t267;
t850 = t291 + t366;
t849 = -t356 + t359;
t848 = t358 + t576;
t661 = rSges(7,1) * t735;
t659 = rSges(7,2) * t736;
t696 = rSges(7,3) * t734 + t471 * t659;
t238 = -t471 * t661 + t696;
t695 = rSges(7,3) * t733 + t473 * t659;
t239 = -t473 * t661 + t695;
t435 = t464 * rSges(7,3);
t292 = t465 * t593 + t435;
t401 = pkin(9) * t734;
t307 = -pkin(5) * t738 + t401;
t403 = pkin(9) * t733;
t309 = -pkin(5) * t737 + t403;
t635 = -t227 * t443 + t228 * t677 + qJD(2);
t56 = -t192 * t295 + t193 * t296 + t308 * t372 + t310 * t373 + t635;
t637 = t465 * t673;
t638 = t465 * t674;
t836 = g(1) * t473 + g(2) * t471;
t847 = -t56 * (-t192 * t637 - t193 * t638 + t295 * t238 + t239 * t296 + t372 * t307 + t309 * t373) - t64 * (qJD(1) * t309 + t193 * t676 + t416 * t239 - t291 * t637 - t292 * t295 - t372 * t838) - t836 * (-pkin(5) - t791) * t464;
t846 = qJD(1) * t307 - t192 * t676 + t238 * t416 - t291 * t638 + t296 * t292 + t373 * t838 + t850 * t679;
t845 = t28 * t471 + t29 * t473;
t360 = rSges(6,1) * t464 + rSges(6,2) * t465;
t305 = t360 * t471;
t306 = t360 * t473;
t840 = -rSges(6,2) * t738 - t473 * rSges(6,3);
t268 = t662 + t840;
t455 = t471 * rSges(6,3);
t269 = rSges(6,1) * t733 - rSges(6,2) * t737 + t455;
t82 = t268 * t372 + t269 * t373 + t635;
t839 = t465 * rSges(6,1) - rSges(6,2) * t464;
t520 = -t373 * t360 + t551;
t541 = -t268 + t554;
t87 = qJD(1) * t541 + t520;
t88 = -t360 * t372 + (t269 + t553) * qJD(1) + t690;
t844 = -(qJD(1) * t305 - t373 * t839) * t87 - t82 * (-t372 * t305 - t306 * t373) - t88 * (-qJD(1) * t306 - t372 * t839);
t369 = qJD(1) * t382;
t841 = qJD(1) * t260 - t369;
t481 = sin(pkin(11));
t789 = rSges(4,2) * t481;
t793 = rSges(4,1) * t482;
t298 = t471 * rSges(4,3) + (-t789 + t793) * t473;
t241 = t298 + t624;
t386 = t473 * rSges(3,1) - rSges(3,2) * t471;
t347 = t386 + t476;
t454 = Icges(5,4) * t472;
t577 = -Icges(5,2) * t470 + t454;
t378 = Icges(5,1) * t470 + t454;
t834 = qJD(1) * t227 + t841;
t728 = t471 * t472;
t730 = t470 * t471;
t755 = Icges(5,3) * t473;
t277 = Icges(5,5) * t728 - Icges(5,6) * t730 - t755;
t409 = Icges(5,4) * t730;
t763 = Icges(5,5) * t473;
t281 = Icges(5,1) * t728 - t409 - t763;
t758 = Icges(5,6) * t473;
t279 = Icges(5,4) * t728 - Icges(5,2) * t730 - t758;
t744 = t279 * t470;
t562 = -t281 * t472 + t744;
t113 = -t277 * t473 - t471 * t562;
t355 = Icges(6,5) * t465 - Icges(6,6) * t464;
t354 = Icges(6,5) * t464 + Icges(6,6) * t465;
t535 = t354 * t479;
t746 = t265 * t464;
t754 = Icges(6,3) * t473;
t833 = -t473 * t535 + (-t267 * t465 - t355 * t471 + t746 + t754) * qJD(1);
t392 = Icges(6,4) * t738;
t762 = Icges(6,5) * t473;
t266 = Icges(6,1) * t734 - t392 - t762;
t564 = t264 * t464 - t266 * t465;
t263 = Icges(6,3) * t471 + t355 * t473;
t682 = qJD(1) * t263;
t832 = qJD(1) * t564 - t471 * t535 + t682;
t375 = Icges(5,5) * t472 - Icges(5,6) * t470;
t374 = Icges(5,5) * t470 + Icges(5,6) * t472;
t530 = qJD(4) * t374;
t768 = Icges(5,4) * t470;
t379 = Icges(5,1) * t472 - t768;
t282 = Icges(5,5) * t471 + t379 * t473;
t280 = Icges(5,6) * t471 + t473 * t577;
t743 = t280 * t470;
t561 = -t282 * t472 + t743;
t831 = -t473 * t530 + (-t375 * t471 + t561 + t755) * qJD(1);
t278 = Icges(5,3) * t471 + t375 * t473;
t681 = qJD(1) * t278;
t830 = qJD(1) * t562 - t471 * t530 + t681;
t558 = t356 * t464 - t358 * t465;
t829 = qJD(1) * t558 + t355 * t479;
t376 = Icges(5,2) * t472 + t768;
t556 = t376 * t470 - t378 * t472;
t828 = t556 * qJD(1) + t375 * qJD(4);
t827 = t471 * (-t376 * t473 + t282) - t473 * (-Icges(5,2) * t728 + t281 - t409);
t534 = t573 * t465;
t559 = -t287 * t484 + t289 * t486;
t567 = -t186 * t484 + t189 * t486;
t826 = t295 * (-t285 * t473 - t567) - t296 * (-t285 * t471 + t861) + t416 * (Icges(7,3) * t464 + t534 - t559);
t574 = -Icges(7,2) * t486 - t765;
t825 = t295 * (-Icges(7,2) * t335 + t189 + t314) - t296 * (-Icges(7,2) * t333 - t188 - t312) + t416 * (t574 * t464 + t289);
t824 = qJD(1) * t848 + t372 * t851 - t373 * (-Icges(6,2) * t734 + t266 - t392);
t362 = qJDD(4) * t471 + t473 * t668;
t271 = qJDD(5) * t471 + t473 * t667 + t362;
t643 = t464 * t679;
t651 = t465 * t373;
t525 = -t643 + t651;
t133 = qJD(6) * t525 + t473 * t666 + t271;
t93 = Icges(7,5) * t165 + Icges(7,6) * t164 + Icges(7,3) * t526;
t95 = Icges(7,4) * t165 + Icges(7,2) * t164 + Icges(7,6) * t526;
t97 = Icges(7,1) * t165 + Icges(7,4) * t164 + Icges(7,5) * t526;
t22 = (-t479 * t861 - t93) * t465 + (t181 * t479 - t484 * t95 + t486 * t97 + (-t184 * t486 + t188 * t484) * qJD(6)) * t464;
t162 = t473 * t517 + t609 * t726;
t163 = t473 * t516 - t609 * t725;
t92 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t525;
t94 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t525;
t96 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t525;
t23 = (t479 * t567 - t92) * t465 + (t183 * t479 - t484 * t94 + t486 * t96 + (-t186 * t486 - t189 * t484) * qJD(6)) * t464;
t79 = -t183 * t465 + t464 * t567;
t821 = t134 / 0.2e1;
t822 = t133 / 0.2e1;
t823 = -t79 * t133 / 0.2e1 - t78 * t134 / 0.2e1 + t23 * t818 + t22 * t815 - t101 * t821 - t102 * t822;
t820 = t271 / 0.2e1;
t819 = t272 / 0.2e1;
t817 = t295 / 0.2e1;
t816 = -t296 / 0.2e1;
t814 = t339 / 0.2e1;
t813 = t362 / 0.2e1;
t812 = t363 / 0.2e1;
t811 = -t372 / 0.2e1;
t810 = t372 / 0.2e1;
t809 = -t373 / 0.2e1;
t808 = t373 / 0.2e1;
t807 = -t416 / 0.2e1;
t806 = t416 / 0.2e1;
t805 = t471 / 0.2e1;
t804 = -t473 / 0.2e1;
t803 = -rSges(7,3) - pkin(9);
t799 = g(1) * t471;
t796 = -qJD(1) / 0.2e1;
t795 = qJD(1) / 0.2e1;
t792 = rSges(5,1) * t472;
t788 = rSges(5,2) * t472;
t368 = pkin(9) * t651;
t653 = t464 * t373;
t527 = -t465 * t679 - t653;
t200 = pkin(5) * t527 - pkin(9) * t643 + t368;
t169 = -t610 + (-t471 * t689 + t473 * t683) * qJD(1);
t429 = qJ(3) * t678;
t607 = qJDD(1) * t476 - t489 * t802;
t684 = t429 + t444;
t511 = -qJDD(3) * t473 + qJD(1) * (-pkin(2) * t679 + t684) + qJDD(1) * t385 + t471 * t669 + t607;
t506 = qJD(1) * (qJD(1) * t538 - t429) + qJDD(1) * t261 + t511;
t494 = qJD(1) * t169 + qJDD(1) * t228 + (-t362 * t470 - t471 * t723) * pkin(4) + t506;
t648 = t163 * rSges(7,1) + t162 * rSges(7,2) + rSges(7,3) * t651;
t98 = -rSges(7,3) * t643 + t648;
t20 = qJD(1) * t200 + qJDD(1) * t310 - t133 * t291 - t295 * t190 + t339 * t193 - t271 * t366 - t372 * t336 + t416 * t98 + t494;
t784 = t20 * t473;
t783 = t21 * t471;
t456 = t471 * rSges(5,3);
t776 = t471 * t64;
t775 = t473 * t88;
t772 = qJDD(1) / 0.2e1;
t119 = -t285 * t465 + t464 * t559;
t572 = -Icges(7,5) * t484 - Icges(7,6) * t486;
t176 = t479 * t534 + (Icges(7,3) * t479 + qJD(6) * t572) * t464;
t536 = t575 * t465;
t177 = t479 * t536 + (Icges(7,6) * t479 + qJD(6) * t574) * t464;
t537 = t579 * t465;
t578 = -Icges(7,1) * t484 - t764;
t178 = t479 * t537 + (Icges(7,5) * t479 + qJD(6) * t578) * t464;
t48 = (t479 * t559 - t176) * t465 + (-t177 * t484 + t178 * t486 + t285 * t479 + (-t287 * t486 - t289 * t484) * qJD(6)) * t464;
t771 = t119 * t339 + t48 * t416;
t752 = qJD(1) * t56;
t663 = rSges(5,1) * t728;
t688 = rSges(5,2) * t730 + t473 * rSges(5,3);
t283 = t663 - t688;
t552 = -t283 + t606;
t381 = rSges(5,1) * t470 + t788;
t598 = -t381 * t677 + t444;
t117 = qJD(1) * t552 + t598;
t749 = t117 * t471;
t724 = t472 * t473;
t729 = t470 * t473;
t284 = rSges(5,1) * t724 - rSges(5,2) * t729 + t456;
t118 = -t381 * t443 - t445 + (t284 + t605) * qJD(1);
t327 = t381 * t473;
t748 = t118 * t327;
t262 = Icges(6,5) * t734 - Icges(6,6) * t738 - t754;
t747 = t262 * t473;
t742 = t354 * t471;
t741 = t354 * t473;
t740 = t374 * t471;
t739 = t374 * t473;
t716 = -t190 - t336;
t713 = -t471 * t227 + t473 * t228;
t712 = -t471 * t262 - t266 * t733;
t711 = t471 * t263 + t267 * t733;
t710 = -t471 * t277 - t281 * t724;
t709 = t471 * t278 + t282 * t724;
t707 = t471 * t268 + t473 * t269;
t694 = -t376 + t379;
t693 = t378 + t577;
t692 = rSges(6,2) * t643 + rSges(6,3) * t678;
t691 = rSges(5,3) * t678 + t679 * t860;
t418 = t471 * t789;
t687 = rSges(4,3) * t678 + qJD(1) * t418;
t686 = t473 * rSges(4,3) + t418;
t680 = qJD(1) * t375;
t129 = -t471 * t558 - t741;
t671 = t129 * qJD(1);
t158 = -t471 * t556 - t739;
t670 = t158 * qJD(1);
t664 = t471 * t793;
t656 = t472 * t785;
t655 = -m(4) - m(5) - m(6) - m(7);
t160 = rSges(6,1) * t527 - rSges(6,2) * t651 + t692;
t539 = t360 * t479;
t161 = -t471 * t539 + (t473 * t839 + t455) * qJD(1);
t650 = t473 * t160 + t471 * t161 + t268 * t678;
t649 = t473 * t169 + t471 * t170 - t227 * t678;
t646 = t401 + t696;
t645 = t403 + t695;
t644 = t414 - t690;
t641 = t470 * t678;
t634 = -pkin(2) - t793;
t632 = t679 / 0.2e1;
t631 = t678 / 0.2e1;
t630 = -t443 / 0.2e1;
t629 = t443 / 0.2e1;
t628 = -t677 / 0.2e1;
t627 = t677 / 0.2e1;
t528 = -t360 - t801;
t623 = (-t471 ^ 2 - t473 ^ 2) * t470;
t622 = (-t471 * t576 + t757) * qJD(1) + t851 * t479;
t621 = qJD(1) * t265 + t266 * t479 - t356 * t727;
t620 = (-t359 * t471 + t762) * qJD(1) + t852 * t479;
t619 = qJD(1) * t267 + t479 * t853;
t229 = t267 * t734;
t618 = t263 * t473 - t229;
t245 = t282 * t728;
t617 = t278 * t473 - t245;
t616 = -t262 + t746;
t615 = -t277 + t743;
t614 = t848 * t479;
t613 = t849 * t479;
t608 = t471 * t855 + t473 * t854;
t297 = t664 - t686;
t604 = -t297 + t625;
t602 = -t850 - t801;
t331 = t839 * t479;
t601 = -t331 - t656;
t599 = -t471 * t477 + t349 + t476;
t421 = rSges(2,1) * t487 - rSges(2,2) * t485;
t420 = rSges(2,1) * t485 + rSges(2,2) * t487;
t384 = t792 - t860;
t141 = t280 * t472 + t282 * t470;
t531 = qJD(4) * t376;
t196 = -t473 * t531 + (-t471 * t577 + t758) * qJD(1);
t532 = qJD(4) * t378;
t198 = -t473 * t532 + (-t379 * t471 + t763) * qJD(1);
t497 = -qJD(4) * t141 - t196 * t470 + t198 * t472 + t681;
t140 = t279 * t472 + t281 * t470;
t197 = qJD(1) * t280 - t471 * t531;
t199 = qJD(1) * t282 - t471 * t532;
t498 = qJD(1) * t277 - qJD(4) * t140 - t197 * t470 + t199 * t472;
t591 = t471 * (t471 * t831 + t497 * t473) - t473 * (t471 * t830 + t498 * t473);
t590 = t471 * (t497 * t471 - t473 * t831) - t473 * (t498 * t471 - t473 * t830);
t63 = qJD(1) * t540 + t864;
t589 = t473 * t63 + t776;
t588 = t471 * t70 - t473 * t69;
t587 = t471 * t69 + t473 * t70;
t586 = t471 * t72 - t473 * t71;
t585 = t471 * t71 + t473 * t72;
t584 = t471 * t79 - t473 * t78;
t583 = t471 * t78 + t473 * t79;
t582 = -t471 * t88 - t473 * t87;
t114 = -t280 * t730 - t617;
t571 = -t113 * t473 + t114 * t471;
t115 = -t279 * t729 - t710;
t116 = -t280 * t729 + t709;
t570 = -t115 * t473 + t116 * t471;
t569 = -t117 * t473 - t118 * t471;
t566 = -t192 * t473 - t193 * t471;
t202 = -t677 * t788 + (-t472 * t679 - t640) * rSges(5,1) + t691;
t326 = t381 * t471;
t203 = -qJD(4) * t326 + (t384 * t473 + t456) * qJD(1);
t565 = t202 * t473 + t203 * t471;
t127 = t264 * t465 + t266 * t464;
t560 = t283 * t471 + t284 * t473;
t557 = t376 * t472 + t378 * t470;
t550 = -t656 + t716;
t549 = t855 * t678 + (t200 + t98) * t473 + (t201 + t99) * t471;
t547 = t464 * t803 - t438;
t533 = t564 * t471;
t529 = t169 * t677 + t170 * t443 - t362 * t227 - t228 * t363 + qJDD(2);
t524 = t292 + t838;
t523 = -t181 * t296 + t183 * t295 + t285 * t416;
t522 = (-Icges(7,5) * t332 - Icges(7,6) * t333) * t296 - (Icges(7,5) * t334 - Icges(7,6) * t335) * t295 - t572 * t464 * t416;
t521 = qJD(1) * t355 - t372 * t741 + t373 * t742;
t519 = t279 * t473 - t280 * t471;
t518 = t464 * t522;
t512 = (-t470 * t693 + t472 * t694) * qJD(1);
t505 = (Icges(7,1) * t334 - t186 - t766) * t295 - (-Icges(7,1) * t332 - t184 - t313) * t296 + (t578 * t464 - t287) * t416;
t503 = qJD(1) * t849 + t372 * t852 - t373 * t853;
t501 = qJD(1) * t262 - t464 * t621 + t465 * t619;
t500 = -t464 * t622 + t465 * t620 + t682;
t499 = qJD(1) * t354 - t464 * t614 + t465 * t613;
t351 = t577 * qJD(4);
t352 = t379 * qJD(4);
t496 = qJD(1) * t374 - qJD(4) * t557 - t351 * t470 + t352 * t472;
t495 = -t470 * t827 + t519 * t472;
t106 = -t533 - t747;
t107 = -t265 * t738 - t618;
t108 = -t264 * t737 - t712;
t109 = -t265 * t737 + t711;
t128 = t265 * t465 + t267 * t464;
t130 = -t473 * t558 + t742;
t16 = t162 * t184 - t163 * t188 + t181 * t525 + t334 * t95 + t335 * t97 + t737 * t93;
t17 = t162 * t186 + t163 * t189 + t183 * t525 + t334 * t94 + t335 * t96 + t737 * t92;
t18 = t164 * t184 - t165 * t188 + t181 * t526 - t332 * t95 + t333 * t97 + t738 * t93;
t19 = t164 * t186 + t165 * t189 + t183 * t526 - t332 * t94 + t333 * t96 + t738 * t92;
t234 = t287 * t471;
t235 = t287 * t473;
t236 = t289 * t471;
t237 = t289 * t473;
t288 = Icges(7,6) * t464 + t536;
t290 = Icges(7,5) * t464 + t537;
t38 = t162 * t287 + t163 * t289 + t176 * t737 + t177 * t334 + t178 * t335 + t285 * t525;
t3 = t102 * t339 + t133 * t72 + t134 * t71 - t16 * t296 + t17 * t295 + t38 * t416;
t32 = t119 * t416 + t295 * t79 - t296 * t78;
t39 = t164 * t287 + t165 * t289 + t176 * t738 - t177 * t332 + t178 * t333 + t285 * t526;
t4 = t101 * t339 + t133 * t70 + t134 * t69 - t18 * t296 + t19 * t295 + t39 * t416;
t43 = t471 * t832 + t501 * t473;
t44 = t471 * t833 + t500 * t473;
t45 = t501 * t471 - t473 * t832;
t46 = t500 * t471 - t473 * t833;
t491 = -t464 * t824 + t503 * t465;
t492 = t826 * t464;
t57 = -t106 * t373 + t107 * t372 + t671;
t126 = t130 * qJD(1);
t58 = -t108 * t373 + t109 * t372 + t126;
t73 = t471 * t829 + t499 * t473;
t74 = t499 * t471 - t473 * t829;
t75 = t464 * t619 + t465 * t621;
t76 = t464 * t620 + t465 * t622;
t493 = (-t127 * t473 + t128 * t471) * t772 + (t503 * t464 + t465 * t824) * t796 + (((t235 * t484 - t237 * t486 + t183) * t295 - (t234 * t484 - t236 * t486 + t181) * t296 + (-t288 * t484 + t290 * t486 + t285) * t416 + t119 * qJD(6)) * t464 + (qJD(6) * t583 - t826) * t465) * t807 + (qJD(1) * t583 - t22 * t473 + t23 * t471) * t806 + (t471 * t491 - t473 * t521) * t808 + (-t45 * t473 + t46 * t471 + (t106 * t471 + t107 * t473) * qJD(1)) * t809 + (t471 * t76 - t473 * t75 + (t127 * t471 + t128 * t473) * qJD(1)) * t795 + t586 * t822 + (-t43 * t473 + t44 * t471 + (t108 * t471 + t109 * t473) * qJD(1)) * t810 + (t471 * t521 + t473 * t491) * t811 + t584 * t814 + ((t235 * t332 - t237 * t333) * t295 - (t234 * t332 - t236 * t333) * t296 + (-t288 * t332 + t290 * t333) * t416 + (t101 * t464 + t70 * t733) * qJD(6) + ((qJD(6) * t69 + t523) * t465 + t492) * t471) * t815 + (qJD(1) * t587 - t18 * t473 + t19 * t471) * t816 + (qJD(1) * t585 - t16 * t473 + t17 * t471) * t817 + ((-t235 * t334 - t237 * t335) * t295 - (-t234 * t334 - t236 * t335) * t296 + (t288 * t334 + t290 * t335) * t416 + (t102 * t464 + t71 * t734) * qJD(6) + ((qJD(6) * t72 + t523) * t465 + t492) * t473) * t818 + (-t106 * t473 + t107 * t471) * t819 + (-t108 * t473 + t109 * t471) * t820 + t588 * t821 - t32 * t676 / 0.2e1 + (qJD(1) * t73 + qJDD(1) * t130 + t108 * t272 + t109 * t271 + t372 * t44 - t373 * t43 + t3) * t805 + (qJD(1) * t74 + qJDD(1) * t129 + t106 * t272 + t107 * t271 + t372 * t46 - t373 * t45 + t4) * t804 + (t57 + t28) * t632 + (t58 + t29) * t631 - t845 * t675 / 0.2e1;
t353 = t384 * qJD(4);
t345 = t592 * t464;
t224 = rSges(7,1) * t334 - rSges(7,2) * t335;
t223 = -rSges(7,1) * t332 - rSges(7,2) * t333;
t205 = qJD(1) * t241 - t445;
t204 = qJD(1) * t604 + t444;
t159 = -t473 * t556 + t740;
t145 = t159 * qJD(1);
t131 = qJD(4) * t560 + qJD(2);
t105 = qJDD(1) * t298 + qJD(1) * (-qJD(1) * t664 + t687) + t511;
t104 = t604 * qJDD(1) + (-qJD(1) * t298 - t325) * qJD(1) + t555;
t84 = t496 * t471 - t473 * t828;
t83 = t471 * t828 + t496 * t473;
t81 = -qJD(4) * t561 + t196 * t472 + t198 * t470;
t80 = -qJD(4) * t562 + t197 * t472 + t199 * t470;
t77 = qJD(4) * t565 + t283 * t362 - t284 * t363 + qJDD(2);
t68 = qJD(1) * t202 + qJDD(1) * t284 - t353 * t443 - t362 * t381 + t506;
t67 = -t353 * t677 + t363 * t381 + (-t203 + t708) * qJD(1) + t552 * qJDD(1) + t555;
t62 = qJD(4) * t570 + t145;
t61 = qJD(4) * t571 + t670;
t41 = qJD(1) * t160 + qJDD(1) * t269 - t271 * t360 - t331 * t372 + t494;
t40 = t272 * t360 - t331 * t373 + (-t161 + t647) * qJD(1) + t541 * qJDD(1) + t509;
t37 = t160 * t373 + t161 * t372 + t268 * t271 - t269 * t272 + t529;
t15 = -t133 * t192 - t134 * t193 + t200 * t373 + t201 * t372 + t271 * t308 - t272 * t310 + t295 * t99 + t296 * t98 + t529;
t1 = [(t204 * t445 + t205 * (t684 + t687) + ((-t204 * t487 - t205 * t485) * pkin(1) + t204 * (t634 + t789) * t473 + (t204 * (-rSges(4,3) - qJ(3)) + t205 * t634) * t471) * qJD(1) - (-t204 - t369 + t444 + (-t297 - t802) * qJD(1)) * t205 + (t105 - g(2)) * t241 + (t104 - g(1)) * (t471 * t634 + t447 + t686 - t802)) * m(4) + (t816 + t815) * t29 + (t130 + t128) * t820 + (-t671 + (t109 - t533 - t711) * t373 + (t471 * t616 + t108 - t229) * t372 + ((t263 + t564) * t372 + t616 * t373) * t473 + t57) * t811 + (t80 + t84 + t62) * t628 + (t75 + t74 + t58) * t809 + (t81 + t83) * t629 + (-qJD(4) * t556 + t351 * t472 + t352 * t470 + t464 * t613 + t465 * t614) * qJD(1) + (t159 + t141) * t813 + ((-t383 * t489 - g(2) + t607) * t347 + (-t863 + (-0.2e1 * t386 - t476 + t347) * t489 - g(1)) * t346) * m(3) + (t61 - t670 + ((t473 * t615 + t116 - t709) * t473 + (t471 * t615 + t115 + t617) * t471) * qJD(4)) * t630 + (t158 + t140) * t812 + (t145 + ((t114 - t245 + (t278 + t744) * t473 + t710) * t473 + t709 * t471) * qJD(4)) * t627 - m(2) * (-g(1) * t420 + g(2) * t421) + ((-t539 - t657) * t775 + (rSges(6,1) * t654 + rSges(6,2) * t652 + t644 + (-t455 - t476 + (-t400 - t839) * t473) * qJD(1)) * t87 + (t444 + t692 + t87 - t520 - t834 + (t268 + t802 + t866) * qJD(1)) * t88 + (t41 - g(2)) * (t269 + t599) + (t40 - g(1)) * (-t840 + t866)) * m(6) + (t557 + t356 * t465 + t358 * t464 + m(3) * (t346 ^ 2 + t386 * t347) + m(2) * (t420 ^ 2 + t421 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t482 ^ 2 + (Icges(4,1) * t481 + 0.2e1 * Icges(4,4) * t482) * t481) * qJDD(1) + ((t381 * t749 - t748) * qJD(4) + (t68 - g(2)) * (t284 + t476 + t612) + (t67 - g(1)) * (-t802 - t722 + (-t466 - t792) * t471 + t688) + (t417 + t445 + (-t456 - t476 + (-t384 - t466) * t473) * qJD(1)) * t117 + (t444 + t691 + t117 - t598 - t841 + (t283 - t663 - t865) * qJD(1)) * t118) * m(5) + (t126 + (t107 + (t264 * t473 + t265 * t471) * t464 + t618 + t712) * t373 + (-t266 * t734 + t747 + t106 + (t264 * t471 - t265 * t473) * t464 + t711) * t372) * t808 + (-(-t63 + (-t308 - t802) * qJD(1) + t834 + t864) * t64 + t63 * (-t595 + t644) + t64 * (-pkin(5) * t653 + t368 + t551 + t648) + (t21 * t547 + t63 * (t465 * t803 + t800) * t479) * t471 + ((-t485 * t64 - t487 * t63) * pkin(1) + (t63 * (-t400 - t838 - t435) - t64 * t477) * t473 + (-t400 + t547) * t776) * qJD(1) - t547 * t799 + (t21 - g(1)) * (-t594 + t867) + (t20 - g(2)) * (t599 + t854)) * m(7) + (t129 + t127) * t819 - t823 + t771 + t39 * t816 + t38 * t817 + (t76 + t73) * t810; (m(3) + m(4)) * qJDD(2) + m(5) * t77 + m(6) * t37 + m(7) * t15 + (-m(3) + t655) * g(3); t655 * (-g(2) * t473 + t799) + 0.2e1 * (-t784 / 0.2e1 + t783 / 0.2e1) * m(7) + 0.2e1 * (t40 * t805 + t41 * t804) * m(6) + 0.2e1 * (t67 * t805 + t68 * t804) * m(5) + 0.2e1 * (t104 * t805 + t105 * t804) * m(4); t62 * t631 + t61 * t632 + ((t115 * t471 + t116 * t473) * qJD(1) + t591) * t629 + ((t113 * t471 + t114 * t473) * qJD(1) + t590) * t628 + (qJD(1) * t84 + qJD(4) * t590 + qJDD(1) * t158 + t113 * t363 + t114 * t362) * t804 + t493 + t571 * t812 + t570 * t813 + (-t140 * t473 + t141 * t471) * t772 + ((t470 * t694 + t472 * t693) * qJD(1) + (t519 * t470 + t472 * t827) * qJD(4)) * t796 + ((-t677 * t740 - t680) * t473 + (t512 + (t473 * t739 + t495) * qJD(4)) * t471) * t627 + ((-t443 * t739 + t680) * t471 + (t512 + (t471 * t740 + t495) * qJD(4)) * t473) * t630 + (t471 * t81 - t473 * t80 + (t140 * t471 + t141 * t473) * qJD(1)) * t795 + (qJD(1) * t83 + qJD(4) * t591 + qJDD(1) * t159 + t115 * t363 + t116 * t362) * t805 + (-g(1) * (-pkin(4) * t729 + t645) - g(2) * (-pkin(4) * t730 + t646) - g(3) * (t463 + t524) - (-t64 * t641 + (-t472 * t589 + t56 * t623) * qJD(4)) * pkin(4) + t15 * (t608 + t713) + t56 * (t549 + t649) + t602 * t856 + (t20 * t602 + t64 * t550 + (-t228 - t854) * t752) * t471 + (t473 * t550 + t846) * t63 + t847) * m(7) + (-g(3) * (t839 + t463) - t836 * t528 + t37 * (t707 + t713) + t82 * (t649 + t650) + (t601 * t87 + (qJD(1) * t88 + t40) * t528) * t473 + (t41 * t528 + t88 * t601 + (t87 * t360 + t82 * (-t228 - t269)) * qJD(1)) * t471 - (-t88 * t641 + (t472 * t582 + t623 * t82) * qJD(4)) * pkin(4) + t844) * m(6) + (g(1) * t327 + g(2) * t326 - g(3) * t384 + t77 * t560 + t131 * ((t283 * t473 - t284 * t471) * qJD(1) + t565) + t569 * t353 + (-t68 * t471 - t67 * t473 + (-t118 * t473 + t749) * qJD(1)) * t381 - (t117 * t326 - t748) * qJD(1) - (t131 * (-t326 * t471 - t327 * t473) + t569 * t384) * qJD(4)) * m(5); t493 + (-g(1) * t645 - g(2) * t646 - g(3) * t524 + t15 * t608 + t56 * t549 + (t64 * t716 - t752 * t854) * t471 - (t20 * t471 + t856) * t850 + (t473 * t716 + t846) * t63 + t847) * m(7) + (g(1) * t306 + g(2) * t305 - g(3) * t839 + t37 * t707 + t82 * (-t269 * t679 + t650) + t582 * t331 + (-t40 * t473 - t41 * t471 + (t471 * t87 - t775) * qJD(1)) * t360 + t844) * m(6); -t29 * t643 / 0.2e1 + t3 * t737 / 0.2e1 + t4 * t738 / 0.2e1 + (t334 * t825 + t505 * t335 - t473 * t518) * t818 + (-t332 * t825 + t333 * t505 - t471 * t518) * t815 + t845 * t732 / 0.2e1 + ((t479 * t585 - t38) * t817 + (t479 * t587 - t39) * t816 - t771 / 0.2e1 - t119 * t814 + (t479 * t583 - t48) * t806 + t522 * t807 + t823) * t465 + (t585 * t822 + (-qJD(1) * t586 + t102 * t479 + t16 * t471 + t17 * t473) * t817 + t28 * t631 + t587 * t821 + (-qJD(1) * t588 + t101 * t479 + t18 * t471 + t19 * t473) * t816 + t479 * t32 / 0.2e1 + t583 * t814 + (-qJD(1) * t584 + t119 * t479 + t22 * t471 + t23 * t473) * t806 + (-t484 * t825 + t486 * t505) * t807) * t464 + ((-t21 * t192 - t20 * t193 + t63 * t99 - t64 * t98 + (t56 * t566 + (t471 * t63 - t473 * t64) * t291) * t479) * t465 + (t63 * (t190 * t471 + t192 * t479) + t64 * (-t190 * t473 + t193 * t479) + t15 * t566 + t56 * (t192 * t679 - t193 * t678 - t471 * t98 + t473 * t99) + (qJD(1) * t589 + t783 - t784) * t291) * t464 - t63 * (-t223 * t416 - t296 * t345) - t64 * (t224 * t416 - t295 * t345) - t56 * (t223 * t295 + t224 * t296) - g(1) * t224 - g(2) * t223 - g(3) * t345) * m(7);];
tau  = t1;
