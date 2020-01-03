% Calculate vector of inverse dynamics joint torques for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR9_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR9_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:24:36
% DurationCPUTime: 48.23s
% Computational Cost: add. (18069->958), mult. (23546->1218), div. (0->0), fcn. (21237->8), ass. (0->450)
t774 = Icges(5,1) + Icges(4,3);
t390 = sin(qJ(3));
t393 = cos(qJ(3));
t311 = Icges(4,5) * t393 - Icges(4,6) * t390;
t467 = Icges(5,4) * t393 - Icges(5,5) * t390;
t773 = t311 - t467;
t388 = qJ(1) + pkin(8);
t380 = cos(t388);
t772 = t774 * t380;
t379 = sin(t388);
t613 = t379 * t393;
t614 = t379 * t390;
t747 = t772 + (-Icges(5,5) + Icges(4,6)) * t614 + (Icges(5,4) - Icges(4,5)) * t613;
t753 = t774 * t379 + t773 * t380;
t629 = Icges(4,6) * t380;
t183 = Icges(4,4) * t613 - Icges(4,2) * t614 - t629;
t336 = Icges(5,6) * t614;
t637 = Icges(5,4) * t380;
t190 = Icges(5,2) * t613 - t336 + t637;
t771 = t183 * t390 - t190 * t393;
t639 = Icges(4,4) * t390;
t317 = Icges(4,1) * t393 - t639;
t186 = Icges(4,5) * t379 + t317 * t380;
t627 = Icges(5,6) * t393;
t461 = -Icges(5,3) * t390 + t627;
t187 = Icges(5,5) * t379 - t380 * t461;
t770 = -t186 * t613 - t187 * t614;
t314 = Icges(4,2) * t393 + t639;
t628 = Icges(5,6) * t390;
t460 = Icges(5,3) * t393 + t628;
t769 = -t314 - t460;
t383 = Icges(4,4) * t393;
t316 = Icges(4,1) * t390 + t383;
t462 = Icges(5,2) * t390 + t627;
t768 = t316 + t462;
t463 = Icges(5,2) * t393 - t628;
t767 = t317 + t463;
t341 = Icges(4,4) * t614;
t633 = Icges(4,5) * t380;
t185 = Icges(4,1) * t613 - t341 - t633;
t632 = Icges(5,5) * t380;
t188 = Icges(5,6) * t613 - Icges(5,3) * t614 + t632;
t744 = -t185 * t393 + t188 * t390 + t771;
t468 = -Icges(4,2) * t390 + t383;
t766 = t461 + t468;
t445 = t314 * t390 - t316 * t393;
t756 = -t390 * t460 + t393 * t462 - t445;
t765 = t380 * t753 + t770;
t611 = t380 * t393;
t612 = t380 * t390;
t723 = t186 * t611 + t187 * t612 + t379 * t753;
t764 = -t185 * t611 + t188 * t612 + t379 * t747;
t719 = -t744 * t379 + t747 * t380;
t184 = Icges(4,6) * t379 + t380 * t468;
t337 = Icges(5,6) * t612;
t638 = Icges(5,4) * t379;
t189 = -Icges(5,2) * t611 + t337 + t638;
t718 = -t184 * t614 - t189 * t613 - t765;
t717 = -t183 * t612 + t190 * t611 - t764;
t716 = -t184 * t612 - t189 * t611 + t723;
t763 = t184 - t187;
t762 = t766 * qJD(3);
t761 = t767 * qJD(3);
t310 = Icges(4,5) * t390 + Icges(4,6) * t393;
t466 = Icges(5,4) * t390 + Icges(5,5) * t393;
t760 = t310 - t466;
t758 = t769 * qJD(3);
t757 = t768 * qJD(3);
t755 = -t768 * t390 + t769 * t393;
t615 = t310 * t380;
t136 = -t379 * t445 - t615;
t240 = t466 * t380;
t139 = t460 * t614 - t462 * t613 - t240;
t754 = t136 - t139;
t239 = t466 * t379;
t616 = t310 * t379;
t728 = t756 * t380 - t239 + t616;
t752 = t184 * t390 + t189 * t393;
t715 = (t183 + t188) * t393 + (t185 + t190) * t390;
t714 = t763 * t393 + (t186 - t189) * t390;
t751 = t758 * t380 + (-t766 * t379 + t629 - t632) * qJD(1);
t750 = t763 * qJD(1) + t758 * t379;
t749 = -t757 * t380 + (-t767 * t379 + t633 - t637) * qJD(1);
t748 = t757 * t379 + (-t380 * t463 - t186 + t638) * qJD(1);
t746 = t760 * qJD(1) + t755 * qJD(3) - t762 * t390 + t761 * t393;
t745 = t760 * qJD(3);
t743 = t186 * t393 + t187 * t390 - t752;
t742 = t716 * t379 - t717 * t380;
t741 = t718 * t379 - t719 * t380;
t740 = t756 * qJD(1) - qJD(3) * t773;
t394 = cos(qJ(1));
t387 = t394 * pkin(1);
t739 = t728 * qJD(1);
t275 = rSges(3,1) * t379 + rSges(3,2) * t380;
t391 = sin(qJ(1));
t669 = pkin(1) * t391;
t256 = -t275 - t669;
t738 = t754 * qJD(1);
t737 = t753 * qJD(1);
t389 = sin(qJ(5));
t392 = cos(qJ(5));
t609 = t390 * t392;
t225 = t379 * t609 + t380 * t389;
t610 = t389 * t390;
t226 = -t379 * t610 + t380 * t392;
t584 = t226 * rSges(6,1) - t225 * rSges(6,2);
t123 = rSges(6,3) * t613 - t584;
t490 = rSges(6,1) * t389 + rSges(6,2) * t392;
t705 = t393 * t490;
t254 = rSges(6,3) * t390 - t705;
t552 = qJD(5) * t393;
t557 = qJD(3) * t380;
t260 = -t379 * t552 + t557;
t381 = qJD(4) * t390;
t322 = t380 * t381;
t553 = qJD(5) * t390;
t363 = qJD(1) + t553;
t324 = pkin(3) * t390 - qJ(4) * t393;
t512 = -pkin(7) * t390 - t324;
t736 = t123 * t363 + t260 * t254 - t512 * t557 - t322;
t735 = t380 ^ 2;
t396 = qJD(1) ^ 2;
t544 = t396 * t387;
t734 = qJD(3) * t741 + t738;
t733 = qJD(3) * t742 + t739;
t732 = qJD(3) * t744 + t390 * t748 - t393 * t750;
t731 = qJD(3) * t743 + t390 * t749 + t393 * t751;
t730 = -t379 * t740 + t380 * t746;
t729 = t379 * t746 + t380 * t740;
t727 = t747 * qJD(1) + qJD(3) * t715 + t750 * t390 + t748 * t393;
t726 = -qJD(3) * t714 - t390 * t751 + t393 * t749 + t737;
t725 = t747 + t752;
t382 = t390 * qJ(4);
t700 = t393 * pkin(3) + t382;
t248 = t700 * t379;
t373 = t380 * pkin(6);
t277 = pkin(2) * t379 - t373;
t270 = qJD(1) * t277;
t561 = qJD(1) * t380;
t358 = pkin(6) * t561;
t555 = qJD(3) * t393;
t528 = t380 * t555;
t576 = qJ(4) * t528 + t322;
t724 = qJD(1) * t248 + t270 + t358 + t576;
t722 = -t745 * t380 + (-t379 * t773 - t743 + t772) * qJD(1);
t721 = qJD(1) * t744 - t379 * t745 + t737;
t368 = t379 * rSges(5,1);
t203 = -rSges(5,2) * t611 + rSges(5,3) * t612 + t368;
t252 = pkin(3) * t611 + t380 * t382;
t278 = t380 * pkin(2) + t379 * pkin(6);
t699 = t387 + t278;
t710 = t252 + t699;
t92 = t203 + t710;
t266 = pkin(4) * t380 - pkin(7) * t613;
t514 = t373 - t669;
t720 = t266 + t514 + t584;
t558 = qJD(3) * t379;
t259 = t380 * t552 + t558;
t223 = -t379 * t389 + t380 * t609;
t224 = t379 * t392 + t380 * t610;
t110 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t611;
t636 = Icges(6,4) * t224;
t113 = Icges(6,2) * t223 + Icges(6,6) * t611 + t636;
t210 = Icges(6,4) * t223;
t116 = Icges(6,1) * t224 + Icges(6,5) * t611 + t210;
t35 = t110 * t611 + t223 * t113 + t224 * t116;
t112 = -Icges(6,5) * t226 + Icges(6,6) * t225 + Icges(6,3) * t613;
t212 = Icges(6,4) * t226;
t115 = Icges(6,2) * t225 + Icges(6,6) * t613 - t212;
t211 = Icges(6,4) * t225;
t117 = Icges(6,1) * t226 - Icges(6,5) * t613 - t211;
t36 = t112 * t611 + t223 * t115 - t117 * t224;
t464 = Icges(6,5) * t389 + Icges(6,6) * t392;
t412 = -Icges(6,3) * t390 + t393 * t464;
t635 = Icges(6,4) * t389;
t465 = Icges(6,2) * t392 + t635;
t413 = -Icges(6,6) * t390 + t393 * t465;
t634 = Icges(6,4) * t392;
t469 = Icges(6,1) * t389 + t634;
t414 = -Icges(6,5) * t390 + t393 * t469;
t66 = -t223 * t413 - t224 * t414 - t412 * t611;
t12 = t259 * t35 - t260 * t36 + t66 * t363;
t37 = t110 * t613 + t225 * t113 - t226 * t116;
t38 = t112 * t613 + t115 * t225 + t117 * t226;
t67 = -t225 * t413 + t226 * t414 - t412 * t613;
t13 = t259 * t37 - t260 * t38 + t363 * t67;
t458 = t115 * t392 - t117 * t389;
t43 = t112 * t390 - t393 * t458;
t367 = t379 * rSges(4,3);
t202 = rSges(4,1) * t611 - rSges(4,2) * t612 + t367;
t155 = t202 + t699;
t709 = t727 * t735 + (t722 * t379 + (-t721 + t726) * t380) * t379;
t708 = t721 * t735 + (t726 * t379 + (-t722 + t727) * t380) * t379;
t548 = qJD(3) * qJD(4);
t707 = qJDD(4) * t390 + t393 * t548;
t556 = qJD(3) * t390;
t529 = t380 * t556;
t559 = qJD(1) * t393;
t534 = t379 * t559;
t706 = t529 + t534;
t704 = t12 * t380 + t13 * t379;
t562 = qJD(1) * t379;
t276 = t380 * rSges(3,1) - rSges(3,2) * t379;
t257 = t276 + t387;
t701 = -pkin(7) * t393 - t700;
t265 = t379 * pkin(4) + pkin(7) * t611;
t384 = t390 * rSges(5,3);
t489 = -rSges(5,2) * t393 + t384;
t513 = -pkin(2) - t382;
t670 = -rSges(6,3) - pkin(3);
t545 = -pkin(7) + t670;
t698 = t545 * t393 + t513;
t121 = t224 * rSges(6,1) + t223 * rSges(6,2) + rSges(6,3) * t611;
t263 = t324 * t558;
t532 = t379 * t556;
t303 = pkin(7) * t532;
t527 = t379 * t381;
t41 = t527 + t121 * t363 - t254 * t259 - t263 - t303 + (t265 + t710) * qJD(1);
t516 = -t277 - t669;
t497 = -t248 + t516;
t443 = t266 + t497;
t40 = qJD(1) * t443 - t736;
t647 = t380 * t40;
t697 = t379 * t41 + t647;
t588 = -Icges(5,3) * t613 + t190 - t336;
t590 = t462 * t379 + t188;
t688 = -t390 * t588 - t393 * t590;
t593 = -Icges(4,2) * t613 + t185 - t341;
t595 = t316 * t379 + t183;
t687 = -t390 * t593 - t393 * t595;
t272 = (Icges(6,2) * t389 - t634) * t393;
t406 = t259 * (-Icges(6,2) * t224 + t116 + t210) - t260 * (Icges(6,2) * t226 - t117 + t211) + t363 * (-t414 + t272);
t273 = (-Icges(6,1) * t392 + t635) * t393;
t407 = t259 * (-Icges(6,1) * t223 + t113 + t636) - t260 * (-Icges(6,1) * t225 + t115 - t212) + t363 * (-t413 - t273);
t686 = m(5) / 0.2e1;
t685 = m(6) / 0.2e1;
t684 = -m(5) - m(6);
t549 = qJD(1) * qJD(3);
t267 = qJDD(3) * t379 + t380 * t549;
t546 = qJDD(5) * t393;
t144 = -qJD(5) * t706 + t380 * t546 + t267;
t683 = t144 / 0.2e1;
t268 = -qJDD(3) * t380 + t379 * t549;
t533 = t380 * t559;
t424 = -t532 + t533;
t145 = qJD(5) * t424 + t379 * t546 + t268;
t682 = t145 / 0.2e1;
t681 = -t259 / 0.2e1;
t680 = t259 / 0.2e1;
t679 = -t260 / 0.2e1;
t678 = t260 / 0.2e1;
t677 = t267 / 0.2e1;
t676 = t268 / 0.2e1;
t279 = qJD(3) * t552 + qJDD(5) * t390 + qJDD(1);
t675 = t279 / 0.2e1;
t674 = -t363 / 0.2e1;
t673 = t363 / 0.2e1;
t667 = g(2) * t379;
t459 = t113 * t392 + t116 * t389;
t415 = -t363 * t389 + t392 * t555;
t560 = qJD(1) * t390;
t501 = qJD(5) + t560;
t442 = t379 * t501;
t105 = t380 * t415 - t392 * t442;
t416 = t363 * t392 + t389 * t555;
t106 = t380 * t416 - t389 * t442;
t55 = Icges(6,5) * t106 + Icges(6,6) * t105 - Icges(6,3) * t706;
t57 = Icges(6,4) * t106 + Icges(6,2) * t105 - Icges(6,6) * t706;
t59 = Icges(6,1) * t106 + Icges(6,4) * t105 - Icges(6,5) * t706;
t8 = (qJD(3) * t459 + t55) * t390 + (qJD(3) * t110 - t389 * t59 - t392 * t57 + (t113 * t389 - t116 * t392) * qJD(5)) * t393;
t666 = t8 * t259;
t441 = t380 * t501;
t103 = t379 * t415 + t392 * t441;
t104 = t379 * t416 + t389 * t441;
t54 = Icges(6,5) * t104 + Icges(6,6) * t103 + Icges(6,3) * t424;
t56 = Icges(6,4) * t104 + Icges(6,2) * t103 + Icges(6,6) * t424;
t58 = Icges(6,1) * t104 + Icges(6,4) * t103 + Icges(6,5) * t424;
t9 = (qJD(3) * t458 + t54) * t390 + (qJD(3) * t112 - t389 * t58 - t392 * t56 + (t115 * t389 + t117 * t392) * qJD(5)) * t393;
t665 = t9 * t260;
t227 = Icges(6,3) * t393 + t390 * t464;
t271 = (-Icges(6,5) * t392 + Icges(6,6) * t389) * t393;
t156 = qJD(3) * t227 + qJD(5) * t271;
t229 = Icges(6,6) * t393 + t390 * t465;
t157 = qJD(3) * t229 + qJD(5) * t272;
t231 = Icges(6,5) * t393 + t390 * t469;
t158 = qJD(3) * t231 + qJD(5) * t273;
t449 = -t389 * t414 - t392 * t413;
t33 = (qJD(3) * t449 + t156) * t390 + (-qJD(3) * t412 - t157 * t392 - t158 * t389 + (-t389 * t413 + t392 * t414) * qJD(5)) * t393;
t78 = -t390 * t412 - t393 * t449;
t662 = t78 * t279 + t33 * t363;
t661 = rSges(4,1) * t393;
t658 = rSges(5,2) * t390;
t326 = rSges(4,1) * t390 + rSges(4,2) * t393;
t251 = t326 * t380;
t87 = qJD(1) * t155 - t326 * t558;
t653 = t251 * t87;
t568 = rSges(4,2) * t614 + t380 * rSges(4,3);
t201 = rSges(4,1) * t613 - t568;
t498 = -t201 + t516;
t530 = t326 * t557;
t86 = qJD(1) * t498 - t530;
t648 = t379 * t86;
t646 = t380 * t86;
t385 = t393 * rSges(6,3);
t42 = t110 * t390 - t393 * t459;
t645 = t42 * t144;
t644 = t43 * t145;
t642 = -rSges(5,3) - qJ(4);
t554 = qJD(4) * t393;
t437 = t248 * t558 + t252 * t557 + qJD(2) - t554;
t34 = t121 * t260 + t123 * t259 + (t265 * t380 - t266 * t379) * qJD(3) + t437;
t623 = qJD(3) * t34;
t608 = t393 * qJD(3) ^ 2;
t606 = t106 * rSges(6,1) + t105 * rSges(6,2);
t601 = t121 + t265;
t600 = t123 - t266;
t594 = -t316 * t380 - t184;
t592 = -t314 * t380 + t186;
t591 = -t462 * t380 + t187;
t589 = Icges(5,3) * t611 + t189 + t337;
t586 = -t203 - t252;
t585 = t379 * t248 + t380 * t252;
t334 = qJ(4) * t611;
t249 = -pkin(3) * t612 + t334;
t583 = qJD(1) * t249 + t379 * t554;
t580 = -t252 - t265;
t261 = qJD(3) * t700 - t554;
t579 = -t489 * qJD(3) - t261;
t578 = t705 * t379;
t577 = t705 * t380;
t535 = t379 * t560;
t575 = rSges(4,2) * t535 + rSges(4,3) * t561;
t574 = -t460 + t463;
t573 = -t461 - t462;
t572 = -t314 + t317;
t571 = t316 + t468;
t488 = rSges(5,3) * t393 + t658;
t570 = -t324 + t488;
t569 = -t700 - t489;
t246 = rSges(5,2) * t614 + rSges(5,3) * t613;
t250 = rSges(5,2) * t612 + rSges(5,3) * t611;
t567 = t379 ^ 2 + t735;
t564 = qJD(1) * t311;
t563 = qJD(1) * t467;
t119 = -pkin(3) * t706 - qJ(4) * t535 + t576;
t221 = t379 * t555 + t380 * t560;
t304 = pkin(3) * t532;
t120 = pkin(3) * t533 + qJ(4) * t221 - t304 + t527;
t540 = t380 * t119 + t379 * t120 + t248 * t561;
t539 = t41 * t561;
t332 = qJ(4) * t613;
t245 = -pkin(3) * t614 + t332;
t538 = t245 * t558 + t249 * t557 + t381;
t253 = rSges(6,1) * t610 + rSges(6,2) * t609 + t385;
t526 = -pkin(2) - t661;
t523 = t561 / 0.2e1;
t522 = -t558 / 0.2e1;
t521 = t558 / 0.2e1;
t520 = -t557 / 0.2e1;
t519 = t557 / 0.2e1;
t510 = rSges(5,1) * t380 - rSges(5,3) * t614;
t505 = qJD(3) * t579;
t504 = -qJD(1) * t245 + t380 * t554;
t503 = t380 * t545;
t500 = rSges(5,1) * t561 + rSges(5,2) * t706 + rSges(5,3) * t528;
t499 = qJDD(1) * t387 - t396 * t669;
t495 = -t254 + t512;
t494 = (rSges(5,2) - pkin(3)) * t393 - pkin(2);
t331 = rSges(2,1) * t394 - rSges(2,2) * t391;
t327 = rSges(2,1) * t391 + rSges(2,2) * t394;
t330 = -rSges(4,2) * t390 + t661;
t492 = rSges(6,1) * t104 + rSges(6,2) * t103;
t482 = t35 * t380 + t36 * t379;
t481 = t35 * t379 - t36 * t380;
t480 = t37 * t380 + t379 * t38;
t479 = t37 * t379 - t38 * t380;
t478 = t379 * t43 + t380 * t42;
t477 = t379 * t42 - t380 * t43;
t472 = -t379 * t87 - t646;
t457 = t121 * t379 - t123 * t380;
t140 = -rSges(4,1) * t706 - rSges(4,2) * t528 + t575;
t247 = t326 * t379;
t141 = -qJD(3) * t247 + (t330 * t380 + t367) * qJD(1);
t456 = t140 * t380 + t141 * t379;
t450 = t201 * t379 + t202 * t380;
t204 = rSges(5,2) * t613 + t510;
t444 = t204 + t497;
t274 = (-rSges(6,1) * t392 + rSges(6,2) * t389) * t393;
t159 = qJD(5) * t274 + (t390 * t490 + t385) * qJD(3);
t440 = -pkin(7) * t555 - t159 - t261;
t439 = -pkin(2) - t700;
t262 = t278 * qJD(1);
t438 = -t120 - t262 - t527;
t434 = t557 * t570 + t322;
t433 = t268 * t324 + t380 * t707 - t544;
t432 = qJD(1) * (-pkin(2) * t562 + t358) + qJDD(1) * t278 + t499;
t422 = -t110 * t259 + t112 * t260 + t363 * t412;
t421 = (Icges(6,5) * t223 - Icges(6,6) * t224) * t259 - (Icges(6,5) * t225 + Icges(6,6) * t226) * t260 + t271 * t363;
t420 = -qJDD(4) * t393 + t119 * t557 + t120 * t558 + t267 * t248 + t390 * t548 + qJDD(2);
t419 = -t390 * t592 + t393 * t594;
t418 = t390 * t589 + t393 * t591;
t417 = t393 * t421;
t411 = (t390 * t573 + t393 * t574) * qJD(1);
t410 = (-t390 * t571 + t393 * t572) * qJD(1);
t408 = qJDD(1) * t252 + t432 + t707 * t379 + (t119 + t322) * qJD(1);
t405 = t34 * t457 + (-t379 * t40 + t380 * t41) * t254;
t398 = (t412 * t380 + t459) * t259 - (t412 * t379 + t458) * t260 + (t227 + t449) * t363;
t397 = t398 * t393;
t359 = pkin(4) * t561;
t287 = t330 * qJD(3);
t264 = t324 * t562;
t222 = t528 - t535;
t220 = t567 * t556;
t180 = -pkin(7) * t706 + t359;
t179 = qJD(1) * t265 - t303;
t174 = -rSges(6,3) * t612 + t577;
t173 = -rSges(6,3) * t614 + t578;
t172 = t414 * t380;
t171 = t414 * t379;
t170 = t413 * t380;
t169 = t413 * t379;
t153 = rSges(6,1) * t225 + rSges(6,2) * t226;
t152 = rSges(6,1) * t223 - rSges(6,2) * t224;
t143 = -rSges(5,3) * t535 + t500;
t142 = t488 * t558 + (t380 * t489 + t368) * qJD(1);
t85 = qJD(3) * t450 + qJD(2);
t65 = -t263 + (qJD(3) * t488 + t381) * t379 + t92 * qJD(1);
t64 = qJD(1) * t444 + t434;
t62 = (t203 * t380 - t204 * t379) * qJD(3) + t437;
t61 = -rSges(6,3) * t706 + t606;
t60 = rSges(6,3) * t424 + t492;
t49 = qJD(1) * t140 + qJDD(1) * t202 - t267 * t326 - t287 * t558 + t432;
t48 = -t544 - t287 * t557 + t268 * t326 + (-t141 - t262) * qJD(1) + t498 * qJDD(1);
t39 = qJD(3) * t456 + t201 * t267 - t202 * t268 + qJDD(2);
t27 = qJD(1) * t143 + qJDD(1) * t203 + t267 * t570 + t379 * t505 + t408;
t26 = -t268 * t488 + t380 * t505 + t444 * qJDD(1) + (-t142 + t438) * qJD(1) + t433;
t17 = -t105 * t413 - t106 * t414 + t156 * t611 + t157 * t223 + t158 * t224 + t412 * t706;
t16 = -t103 * t413 - t104 * t414 + t156 * t613 + t157 * t225 - t158 * t226 - t412 * t424;
t15 = -t204 * t267 + t586 * t268 + (t142 * t379 + t143 * t380) * qJD(3) + t420;
t14 = t259 * t42 - t260 * t43 + t363 * t78;
t11 = t408 + (-t267 * t390 - t379 * t608) * pkin(7) - t261 * t558 + qJD(1) * t180 - t144 * t254 - t259 * t159 + qJDD(1) * t265 + t279 * t121 - t267 * t324 + t363 * t61;
t10 = -t261 * t557 - t279 * t123 + t145 * t254 - t260 * t159 - t363 * t60 + (t268 * t390 - t380 * t608) * pkin(7) + t443 * qJDD(1) + (-t179 + t438) * qJD(1) + t433;
t7 = t105 * t115 - t106 * t117 - t112 * t706 + t223 * t56 + t224 * t58 + t54 * t611;
t6 = t105 * t113 + t106 * t116 - t110 * t706 + t223 * t57 + t224 * t59 + t55 * t611;
t5 = t103 * t115 - t104 * t117 + t112 * t424 + t225 * t56 - t226 * t58 + t54 * t613;
t4 = t103 * t113 + t104 * t116 + t110 * t424 + t225 * t57 - t226 * t59 + t55 * t613;
t3 = -t121 * t145 + t123 * t144 + t259 * t60 + t260 * t61 - t266 * t267 + t580 * t268 + (t179 * t379 + t180 * t380) * qJD(3) + t420;
t2 = t144 * t35 + t145 * t36 + t17 * t363 + t259 * t6 - t260 * t7 + t279 * t66;
t1 = t144 * t37 + t145 * t38 + t16 * t363 + t259 * t4 - t260 * t5 + t279 * t67;
t18 = [t666 / 0.2e1 - t665 / 0.2e1 + t662 + t644 / 0.2e1 + t645 / 0.2e1 + t16 * t679 + t17 * t680 + t67 * t682 + t66 * t683 - t268 * t139 / 0.2e1 - m(2) * (-g(1) * t327 + g(2) * t331) + ((t723 * t379 + ((t753 + t771) * t380 + t718 + t764 + t770) * t380) * qJD(3) + t739) * t519 + (t679 + t678) * t12 + (t756 * qJD(3) + t761 * t390 + t762 * t393) * qJD(1) + ((-t275 * t396 - g(2) + t499) * t257 + (-t544 + (-0.2e1 * t276 - t387 + t257) * t396 - g(1)) * t256) * m(3) + (t136 + t715) * t676 + (-((t393 * t670 + t513) * t379 + t720) * g(1) + t698 * t647 * qJD(1) + (-g(2) + t11) * (t710 + t601) + (t720 + (t439 - t385) * t379) * t10 + (t303 + t304 - t492 + (rSges(6,3) * t556 - qJ(4) * t555 - t381) * t379 + (-t387 + (-pkin(4) - pkin(6)) * t379) * qJD(1)) * t40 + (t503 * t556 + t359 + t40 + t606 + (t698 * t379 - t266) * qJD(1) + t724 + t736) * t41) * m(6) + ((t304 + (-t381 + (t393 * t642 - t658) * qJD(3)) * t379 + (-t387 + (t390 * t642 + t494) * t380 + (-rSges(5,1) - pkin(6)) * t379) * qJD(1)) * t64 + (-g(2) + t27) * t92 + (-g(1) + t26) * ((t494 - t382) * t379 + t510 + t514) + (-pkin(3) * t529 - t434 + t500 + t64 + ((t439 - t384) * t379 - t204) * qJD(1) + t724) * t65) * m(5) + (t87 * (t358 + t575) + (t326 * t648 - t653) * qJD(3) + ((-t391 * t87 - t394 * t86) * pkin(1) + (-pkin(2) - t330) * t646 + (t86 * (-rSges(4,3) - pkin(6)) + t87 * t526) * t379) * qJD(1) - (-t530 - t270 - t86 + (-t201 - t669) * qJD(1)) * t87 + (t49 - g(2)) * t155 + (t48 - g(1)) * (t526 * t379 + t514 + t568)) * m(4) + (t728 + t714) * t677 + (((t380 * t725 + t716 - t723) * t380 + (t379 * t725 + t717 + t765) * t379) * qJD(3) + t734 - t738) * t522 + (t730 + t731) * t521 + (m(3) * (t256 ^ 2 + t276 * t257) + m(2) * (t327 ^ 2 + t331 ^ 2) + Icges(2,3) + Icges(3,3) - t755) * qJDD(1) + (t729 - t732 + t733) * t520; m(3) * qJDD(2) + m(4) * t39 + m(5) * t15 + m(6) * t3 + (-m(3) - m(4) + t684) * g(3); -((((-t588 - t593) * t380 + (-t589 + t592) * t379) * t393 + ((t590 + t595) * t380 + (t591 + t594) * t379) * t390) * qJD(3) + ((t571 - t573) * t393 + (t572 + t574) * t390) * qJD(1)) * qJD(1) / 0.2e1 + ((-t558 * t615 + t564) * t379 + (t410 + (-t687 * t380 + (t616 + t419) * t379) * qJD(3)) * t380 + (t240 * t558 - t563) * t379 + (t411 + (-t688 * t380 + (-t239 + t418) * t379) * qJD(3)) * t380) * t522 + ((-t557 * t616 - t564) * t380 + (t410 + (t419 * t379 + (t615 - t687) * t380) * qJD(3)) * t379 + (t239 * t557 + t563) * t380 + (t411 + (t418 * t379 + (-t240 - t688) * t380) * qJD(3)) * t379) * t519 + (((-t170 * t392 - t172 * t389 + t110) * t259 - (-t169 * t392 - t171 * t389 + t112) * t260 + (-t229 * t392 - t231 * t389 - t412) * t363 + t78 * qJD(5)) * t393 + (-qJD(5) * t478 + t398) * t390) * t674 + (g(1) * t251 + g(2) * t247 - g(3) * t330 - (t247 * t86 - t653) * qJD(1) - (t85 * (-t247 * t379 - t251 * t380) + t472 * t330) * qJD(3) + t39 * t450 + t85 * ((t201 * t380 - t202 * t379) * qJD(1) + t456) + t472 * t287 + (-t49 * t379 - t48 * t380 + (-t380 * t87 + t648) * qJD(1)) * t326) * m(4) + (t3 * t585 + (t11 * t495 + t3 * t600 + t41 * t440) * t379 + (t3 * t601 + (qJD(1) * t41 + t10) * t495) * t380 - g(1) * (t334 + t577) - g(2) * (t332 + t578) - g(3) * (t253 - t701) - t41 * (t121 * t552 + t174 * t363 - t253 * t259 + t583) - t697 * qJD(3) * t701 + (t123 * t552 + t173 * t363 + t253 * t260 + t254 * t562 + t380 * t440 + t264 - t504) * t40 + (-g(1) * t503 - t545 * t667 - (-t567 * t623 - t539) * pkin(7) - t405 * qJD(5)) * t390 + (t540 + (t179 + t60 + (-t121 + t580) * qJD(1)) * t379 + (qJD(1) * t600 + t180 + t61) * t380 - t173 * t259 - t174 * t260 - t538) * t34) * m(6) + (qJD(1) * t478 + t379 * t8 - t380 * t9) * t673 + t477 * t675 + t704 * t553 / 0.2e1 + ((t379 * t717 + t380 * t716) * qJD(1) + t709) * t521 + ((t379 * t719 + t380 * t718) * qJD(1) + t708) * t520 + (qJD(1) * t730 + t709 * qJD(3) + qJDD(1) * t728 + t716 * t267 + t717 * t268 + t2) * t379 / 0.2e1 + (t732 * t380 + t731 * t379 + (t715 * t379 + t714 * t380) * qJD(1)) * qJD(1) / 0.2e1 + (t12 + t733) * t523 + (t13 + t734) * t562 / 0.2e1 + (t64 * t264 + t15 * t585 + t62 * t540 + (t26 * t570 + t64 * t579 + t15 * t203 + t62 * t143 + (-t62 * t204 + t570 * t65) * qJD(1)) * t380 + (t27 * t570 + t65 * t579 - t15 * t204 + t62 * t142 + (-t488 * t64 + t586 * t62) * qJD(1)) * t379 - g(1) * (t249 + t250) - g(2) * (t245 + t246) + g(3) * t569 - t64 * (-qJD(1) * t246 + t504) - t65 * (qJD(1) * t250 + t583) - t62 * t538 - ((t62 * t250 + t569 * t64) * t380 + (t62 * t246 + t569 * t65) * t379) * qJD(3)) * m(5) - t14 * t552 / 0.2e1 + ((t170 * t225 - t172 * t226) * t259 - (t169 * t225 - t171 * t226) * t260 + (t225 * t229 - t226 * t231) * t363 + (-t37 * t612 + t393 * t67) * qJD(5) + ((-qJD(5) * t38 + t422) * t390 + t397) * t379) * t678 + (qJD(1) * t480 + t379 * t4 - t380 * t5) * t679 + (qJD(1) * t482 + t379 * t6 - t380 * t7) * t680 + ((t170 * t223 + t172 * t224) * t259 - (t169 * t223 + t171 * t224) * t260 + (t223 * t229 + t224 * t231) * t363 + (-t36 * t614 + t393 * t66) * qJD(5) + ((-qJD(5) * t35 + t422) * t390 + t397) * t380) * t681 + t479 * t682 + t481 * t683 - (t729 * qJD(1) + t708 * qJD(3) + t754 * qJDD(1) + t718 * t267 + t719 * t268 + t1) * t380 / 0.2e1 + t741 * t676 + t742 * t677 + (t714 * t379 - t715 * t380) * qJDD(1) / 0.2e1; t684 * (-g(3) * t393 + (g(1) * t380 + t667) * t390) - m(5) * (t220 * t62 + t221 * t65 + t222 * t64) - m(6) * (t220 * t34 + t221 * t41 + t222 * t40) + 0.2e1 * ((t557 * t64 + t558 * t65 - t15) * t686 + (t40 * t557 + t41 * t558 - t3) * t685) * t393 + 0.2e1 * ((qJD(3) * t62 + t26 * t380 + t27 * t379 + t561 * t65 - t562 * t64) * t686 + (t10 * t380 + t11 * t379 - t40 * t562 + t539 + t623) * t685) * t390; -t12 * t534 / 0.2e1 + t2 * t611 / 0.2e1 + (t390 * t66 + t393 * t482) * t683 + ((-qJD(3) * t482 + t17) * t390 + (-qJD(1) * t481 + qJD(3) * t66 + t379 * t7 + t380 * t6) * t393) * t680 + t393 * t13 * t523 + t1 * t613 / 0.2e1 + (t390 * t67 + t393 * t480) * t682 + ((-qJD(3) * t480 + t16) * t390 + (-qJD(1) * t479 + qJD(3) * t67 + t379 * t5 + t380 * t4) * t393) * t679 + t14 * t555 / 0.2e1 + t390 * (t644 + t645 + t662 - t665 + t666) / 0.2e1 + (t390 * t78 + t393 * t478) * t675 + ((-qJD(3) * t478 + t33) * t390 + (-qJD(1) * t477 + qJD(3) * t78 + t379 * t9 + t380 * t8) * t393) * t673 + (t406 * t223 - t224 * t407 + t380 * t417) * t681 + (t225 * t406 + t226 * t407 + t379 * t417) * t678 + (t421 * t390 + (t407 * t389 - t392 * t406) * t393) * t674 - t704 * t556 / 0.2e1 + ((qJD(3) * t405 - t10 * t123 + t11 * t121 - t40 * t60 + t41 * t61) * t390 + (t40 * (-qJD(3) * t123 + t159 * t379) + t41 * (qJD(3) * t121 - t159 * t380) - t3 * t457 + t34 * (-t121 * t561 - t123 * t562 - t379 * t61 + t380 * t60) + (qJD(1) * t697 + t10 * t379 - t11 * t380) * t254) * t393 - t40 * (-t153 * t363 - t260 * t274) - t41 * (t152 * t363 - t259 * t274) - t34 * (t152 * t260 + t153 * t259) - g(1) * t152 - g(2) * t153 - g(3) * t274) * m(6);];
tau = t18;
