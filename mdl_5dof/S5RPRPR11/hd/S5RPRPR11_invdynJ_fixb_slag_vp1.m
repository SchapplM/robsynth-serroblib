% Calculate vector of inverse dynamics joint torques for
% S5RPRPR11
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR11_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:53
% DurationCPUTime: 35.76s
% Computational Cost: add. (16633->838), mult. (22845->1037), div. (0->0), fcn. (20751->8), ass. (0->407)
t784 = Icges(5,4) + Icges(4,5);
t428 = sin(qJ(1));
t429 = cos(qJ(1));
t421 = pkin(8) + qJ(3);
t399 = sin(t421);
t384 = Icges(5,5) * t399;
t400 = cos(t421);
t488 = Icges(5,1) * t400 + t384;
t783 = -t428 * t488 + t784 * t429;
t608 = t399 * t428;
t362 = Icges(4,4) * t608;
t606 = t400 * t428;
t775 = -Icges(4,1) * t606 + t362 + t783;
t633 = Icges(4,4) * t399;
t314 = Icges(4,2) * t400 + t633;
t782 = Icges(5,3) * t400 + t314 - t384;
t628 = Icges(5,5) * t400;
t316 = Icges(5,1) * t399 - t628;
t385 = Icges(4,4) * t400;
t780 = Icges(4,1) * t399 + t316 + t385;
t311 = Icges(4,5) * t400 - Icges(4,6) * t399;
t225 = Icges(4,3) * t428 + t311 * t429;
t313 = Icges(5,4) * t400 + Icges(5,6) * t399;
t227 = Icges(5,2) * t428 + t313 * t429;
t781 = t225 + t227;
t231 = Icges(5,4) * t428 + t429 * t488;
t319 = Icges(4,1) * t400 - t633;
t233 = Icges(4,5) * t428 + t319 * t429;
t774 = t231 + t233;
t619 = Icges(4,3) * t429;
t224 = Icges(4,5) * t606 - Icges(4,6) * t608 - t619;
t309 = Icges(5,3) * t399 + t628;
t222 = -Icges(5,6) * t429 + t309 * t428;
t622 = Icges(4,6) * t429;
t228 = Icges(4,4) * t606 - Icges(4,2) * t608 - t622;
t614 = t228 * t399;
t746 = -t222 * t399 + t775 * t400 + t614;
t779 = -t224 * t429 - t746 * t428;
t777 = t222 - t228;
t605 = t400 * t429;
t361 = Icges(5,5) * t605;
t607 = t399 * t429;
t621 = Icges(5,6) * t428;
t223 = Icges(5,3) * t607 + t361 + t621;
t486 = -Icges(4,2) * t399 + t385;
t229 = Icges(4,6) * t428 + t429 * t486;
t776 = -t223 + t229;
t773 = t309 - t486;
t772 = (Icges(4,6) - Icges(5,6)) * t400 + t784 * t399;
t769 = t319 + t488;
t768 = t782 * qJD(3);
t767 = t780 * qJD(3);
t766 = t223 * t607 + t781 * t428 + t774 * t605;
t226 = -Icges(5,2) * t429 + t313 * t428;
t206 = t428 * t226;
t765 = -t222 * t607 - t428 * t224 + t775 * t605 - t206;
t756 = -t399 * t782 + t780 * t400;
t615 = t226 * t429;
t721 = -t615 + t779;
t720 = -t228 * t607 - t765;
t719 = -t229 * t607 + t766;
t510 = -t223 * t608 + t227 * t429 - t231 * t606;
t193 = t233 * t606;
t517 = t225 * t429 - t193;
t80 = -t229 * t608 - t517;
t718 = -t510 + t80;
t764 = t768 * t429 + (t428 * t486 - t222 - t622) * qJD(1);
t763 = t768 * t428 + (t309 * t429 - t229 + t621) * qJD(1);
t762 = -t767 * t429 + (-t319 * t428 + t783) * qJD(1);
t761 = -t774 * qJD(1) + t767 * t428;
t760 = t773 * qJD(3);
t759 = t769 * qJD(3);
t758 = -t311 - t313;
t757 = t772 * qJD(3);
t755 = -t399 * t780 - t400 * t782;
t613 = t229 * t399;
t754 = t223 * t399 + t774 * t400 - t613;
t710 = t774 * t399 + t776 * t400;
t709 = t775 * t399 + t777 * t400;
t704 = t772 * t429;
t703 = t772 * t428;
t714 = t428 * t756 - t704;
t713 = t429 * t756 + t703;
t753 = t781 * qJD(1);
t695 = qJD(3) - qJD(5);
t345 = t695 * t428;
t553 = qJD(3) * t429;
t346 = -qJD(5) * t429 + t553;
t427 = sin(qJ(5));
t661 = cos(qJ(5));
t533 = t399 * t661;
t258 = t427 * t606 - t428 * t533;
t468 = t399 * t427 + t400 * t661;
t259 = t468 * t428;
t124 = Icges(6,5) * t259 - Icges(6,6) * t258 + Icges(6,3) * t429;
t630 = Icges(6,4) * t259;
t127 = -Icges(6,2) * t258 + Icges(6,6) * t429 + t630;
t631 = Icges(6,4) * t258;
t131 = -Icges(6,1) * t259 - Icges(6,5) * t429 + t631;
t260 = t427 * t605 - t429 * t533;
t261 = t468 * t429;
t46 = -t124 * t428 - t260 * t127 - t131 * t261;
t126 = Icges(6,5) * t261 - Icges(6,6) * t260 - Icges(6,3) * t428;
t221 = Icges(6,4) * t261;
t129 = -Icges(6,2) * t260 - Icges(6,6) * t428 + t221;
t220 = Icges(6,4) * t260;
t132 = Icges(6,1) * t261 - Icges(6,5) * t428 - t220;
t511 = t126 * t428 + t260 * t129 - t261 * t132;
t304 = -t400 * t427 + t533;
t179 = Icges(6,5) * t304 - Icges(6,6) * t468;
t292 = Icges(6,4) * t304;
t182 = -Icges(6,2) * t468 + t292;
t291 = Icges(6,4) * t468;
t185 = Icges(6,1) * t304 - t291;
t58 = -t179 * t428 - t182 * t260 + t185 * t261;
t14 = t58 * qJD(1) - t345 * t511 - t346 * t46;
t175 = t695 * t468;
t558 = qJD(1) * t428;
t101 = t175 * t429 - t304 * t558;
t715 = t695 * t304;
t102 = -t429 * t715 - t468 * t558;
t557 = qJD(1) * t429;
t89 = Icges(6,5) * t175 + Icges(6,6) * t715;
t90 = Icges(6,4) * t175 + Icges(6,2) * t715;
t91 = Icges(6,1) * t175 + Icges(6,4) * t715;
t15 = t101 * t182 + t102 * t185 - t179 * t557 - t260 * t90 + t261 * t91 - t428 * t89;
t103 = t175 * t428 + t304 * t557;
t104 = qJD(1) * t261 - t428 * t715;
t16 = t103 * t182 + t104 * t185 - t179 * t558 - t258 * t90 + t259 * t91 + t429 * t89;
t62 = Icges(6,4) * t104 + Icges(6,2) * t103 - Icges(6,6) * t558;
t64 = Icges(6,1) * t104 + Icges(6,4) * t103 - Icges(6,5) * t558;
t17 = t127 * t715 - t131 * t175 + t304 * t64 - t468 * t62;
t61 = Icges(6,4) * t102 + Icges(6,2) * t101 - Icges(6,6) * t557;
t63 = Icges(6,1) * t102 + Icges(6,4) * t101 - Icges(6,5) * t557;
t18 = t129 * t715 + t132 * t175 + t304 * t63 - t468 * t61;
t546 = qJD(1) * qJD(3);
t330 = qJDD(3) * t428 + t429 * t546;
t545 = qJD(1) * qJD(5);
t254 = -qJDD(5) * t428 - t429 * t545 + t330;
t390 = t428 * t546;
t255 = -t428 * t545 + t390 + (-qJDD(3) + qJDD(5)) * t429;
t44 = t124 * t429 - t127 * t258 - t131 * t259;
t445 = qJD(1) * (Icges(6,1) * t468 + t182 + t292) + t345 * (Icges(6,1) * t260 + t129 + t221) - t346 * (Icges(6,1) * t258 + t127 + t630);
t45 = t429 * t126 - t258 * t129 + t259 * t132;
t455 = qJD(1) * (-Icges(6,5) * t468 - Icges(6,6) * t304) + (Icges(6,5) * t258 + Icges(6,6) * t259) * t346 - (Icges(6,5) * t260 + Icges(6,6) * t261) * t345;
t57 = t179 * t429 - t182 * t258 + t185 * t259;
t60 = Icges(6,5) * t104 + Icges(6,6) * t103 - Icges(6,3) * t558;
t6 = t101 * t127 - t102 * t131 - t124 * t557 - t260 * t62 + t261 * t64 - t428 * t60;
t654 = -qJD(1) / 0.2e1;
t666 = t346 / 0.2e1;
t67 = -t127 * t468 - t131 * t304;
t679 = qJD(1) * (Icges(6,2) * t304 - t185 + t291) - t345 * (-Icges(6,2) * t261 + t132 - t220) + t346 * (-Icges(6,2) * t259 - t131 - t631);
t68 = -t129 * t468 + t132 * t304;
t59 = Icges(6,5) * t102 + Icges(6,6) * t101 - Icges(6,3) * t557;
t7 = t101 * t129 + t102 * t132 - t126 * t557 - t260 * t61 + t261 * t63 - t428 * t59;
t702 = t46 / 0.2e1 + t45 / 0.2e1;
t8 = t103 * t127 - t104 * t131 - t124 * t558 - t258 * t62 + t259 * t64 + t429 * t60;
t9 = t103 * t129 + t104 * t132 - t126 * t558 - t258 * t61 + t259 * t63 + t429 * t59;
t752 = (-t17 * t429 + t18 * t428 + (t428 * t67 + t429 * t68) * qJD(1)) * t654 - qJDD(1) * (t428 * t68 - t429 * t67) / 0.2e1 - t428 * (qJD(1) * t15 + qJDD(1) * t58 + t345 * t7 - t346 * t6) / 0.2e1 + t429 * (qJD(1) * t16 + qJDD(1) * t57 + t345 * t9 - t346 * t8) / 0.2e1 - (qJD(1) * t57 + t345 * t45 - t346 * t44) * t558 / 0.2e1 - t14 * t557 / 0.2e1 + (-t428 * t702 + t429 * t44) * t255 + (t428 * t511 + t429 * t702) * t254 - (t260 * t679 - t261 * t445 + (-t511 * qJD(1) - t6) * t429 + (t46 * qJD(1) - t455 + t7) * t428) * t345 / 0.2e1 + (t258 * t679 - t259 * t445 + (t44 * qJD(1) + t9) * t428 + (t45 * qJD(1) + t455 - t8) * t429) * t666;
t751 = t429 ^ 2;
t504 = t259 * rSges(6,1) - t258 * rSges(6,2);
t133 = rSges(6,3) * t429 + t504;
t658 = pkin(7) * t429;
t320 = pkin(4) * t606 + t658;
t599 = t133 + t320;
t750 = t772 * qJD(1) + t755 * qJD(3) + t760 * t399 + t759 * t400;
t749 = -t710 * qJD(3) + t764 * t399 + t762 * t400 + t753;
t696 = qJD(1) * t226;
t748 = -qJD(1) * t224 - t709 * qJD(3) - t763 * t399 + t761 * t400 - t696;
t747 = -t777 * t429 + (-Icges(5,1) * t607 + t316 * t429 + t361 - t776) * t428;
t745 = t719 * t428 - t720 * t429;
t744 = t718 * t428 - t721 * t429;
t743 = t780 - t773;
t742 = t769 - t782;
t741 = (Icges(4,2) * t606 + t362 + t775) * t429 + (-t314 * t429 + t774) * t428;
t740 = t756 * qJD(1) + t758 * qJD(3);
t739 = t746 * qJD(1) - t757 * t428 + t753;
t738 = -t696 - t757 * t429 + (-t311 * t428 + t619 - t754) * qJD(1);
t737 = -t445 * t304 + t468 * t679;
t187 = -rSges(6,1) * t468 - t304 * rSges(6,2);
t734 = t187 * t345;
t733 = t187 * t346;
t730 = t713 * qJD(1);
t729 = t714 * qJD(1);
t167 = -t258 * rSges(6,1) - t259 * rSges(6,2);
t168 = t260 * rSges(6,1) + t261 * rSges(6,2);
t728 = t167 * t345 - t168 * t346;
t727 = t744 * qJD(3) + t729;
t726 = t745 * qJD(3) + t730;
t725 = -t740 * t428 + t750 * t429;
t724 = t750 * t428 + t740 * t429;
t723 = t746 * qJD(3) + t761 * t399 + t763 * t400;
t722 = t754 * qJD(3) + t762 * t399 - t764 * t400;
t712 = -t741 * t399 + t747 * t400;
t711 = (-t743 * t399 + t742 * t400) * qJD(1);
t708 = t739 * t751 + (t749 * t428 + (-t738 + t748) * t429) * t428;
t707 = t748 * t751 + (t738 * t428 + (-t739 + t749) * t429) * t428;
t544 = qJD(3) * qJD(4);
t706 = qJDD(4) * t399 + t400 * t544;
t705 = t758 * qJD(1);
t506 = t400 * rSges(5,1) + t399 * rSges(5,3);
t383 = t399 * qJ(4);
t698 = t400 * pkin(3) + t383;
t574 = -t698 - t506;
t459 = -pkin(4) * t400 - t698;
t407 = t429 * qJ(2);
t352 = pkin(1) * t428 - t407;
t425 = cos(pkin(8));
t389 = pkin(2) * t425 + pkin(1);
t426 = -pkin(6) - qJ(2);
t394 = t429 * t426;
t569 = -t428 * t389 - t394;
t214 = t352 + t569;
t333 = qJD(1) * t352;
t699 = qJD(1) * t214 - t333;
t372 = pkin(4) * t605;
t321 = -pkin(7) * t428 + t372;
t406 = t428 * qJ(2);
t354 = t429 * pkin(1) + t406;
t697 = t615 + t766;
t285 = t698 * t428;
t689 = -qJD(1) * t285 + t699;
t424 = sin(pkin(8));
t649 = rSges(3,2) * t424;
t651 = rSges(3,1) * t425;
t269 = t428 * rSges(3,3) + (-t649 + t651) * t429;
t655 = g(2) * t428;
t688 = (g(1) * t429 + t655) * t399;
t677 = m(5) / 0.2e1;
t676 = m(6) / 0.2e1;
t675 = -m(5) - m(6);
t674 = -pkin(3) - pkin(4);
t671 = t330 / 0.2e1;
t331 = -qJDD(3) * t429 + t390;
t670 = t331 / 0.2e1;
t665 = t428 / 0.2e1;
t664 = -t429 / 0.2e1;
t663 = -rSges(5,1) - pkin(3);
t662 = rSges(6,3) + pkin(7);
t657 = g(1) * t428;
t652 = pkin(1) - t389;
t650 = rSges(4,1) * t400;
t324 = rSges(4,1) * t399 + rSges(4,2) * t400;
t288 = t324 * t429;
t414 = t428 * rSges(4,3);
t253 = rSges(4,1) * t605 - rSges(4,2) * t607 + t414;
t405 = qJD(2) * t429;
t554 = qJD(3) * t428;
t364 = t429 * t389;
t515 = -t426 * t428 + t364;
t215 = t515 - t354;
t591 = t215 + t354;
t88 = -t324 * t554 - t405 + (t253 + t591) * qJD(1);
t647 = t288 * t88;
t645 = t399 * rSges(5,1);
t416 = t428 * rSges(5,2);
t404 = qJD(2) * t428;
t509 = -t324 * t553 + t404;
t251 = rSges(4,1) * t606 - rSges(4,2) * t608 - t429 * rSges(4,3);
t592 = t214 - t352;
t538 = -t251 + t592;
t87 = qJD(1) * t538 + t509;
t644 = t428 * t87;
t642 = -rSges(5,3) - qJ(4);
t641 = t102 * rSges(6,1) + t101 * rSges(6,2);
t457 = -t399 * t553 - t400 * t558;
t208 = pkin(4) * t457 - pkin(7) * t557;
t65 = -rSges(6,3) * t557 + t641;
t640 = t208 + t65;
t530 = t399 * t554;
t342 = pkin(4) * t530;
t209 = qJD(1) * t321 - t342;
t505 = rSges(6,1) * t104 + rSges(6,2) * t103;
t66 = -rSges(6,3) * t558 + t505;
t639 = t209 + t66;
t604 = t400 * qJD(3) ^ 2;
t586 = t261 * rSges(6,1) - t260 * rSges(6,2);
t135 = -rSges(6,3) * t428 + t586;
t598 = t135 + t321;
t305 = qJD(1) * t354 - t405;
t379 = t426 * t558;
t593 = t379 - (-t429 * t652 - t406) * qJD(1) - t305;
t552 = qJD(4) * t400;
t249 = qJD(3) * t698 - t552;
t584 = -t506 * qJD(3) - t249;
t252 = rSges(5,1) * t605 + rSges(5,3) * t607 + t416;
t373 = pkin(3) * t605;
t289 = qJ(4) * t607 + t373;
t583 = -t252 - t289;
t582 = t428 * t285 + t429 * t289;
t358 = qJ(4) * t605;
t286 = -pkin(3) * t607 + t358;
t551 = qJD(4) * t428;
t581 = qJD(1) * t286 + t400 * t551;
t203 = t269 + t354;
t580 = -t289 - t321;
t322 = pkin(3) * t399 - qJ(4) * t400;
t323 = -rSges(5,3) * t400 + t645;
t575 = -t322 - t323;
t529 = t400 * t553;
t573 = rSges(5,2) * t557 + rSges(5,3) * t529;
t531 = t399 * t558;
t572 = rSges(4,2) * t531 + rSges(4,3) * t557;
t550 = qJD(4) * t429;
t350 = t399 * t550;
t571 = t350 + t404;
t541 = t428 * t651;
t380 = t428 * t649;
t566 = t429 * rSges(3,3) + t380;
t268 = t541 - t566;
t570 = -t352 - t268;
t568 = rSges(3,3) * t557 + qJD(1) * t380;
t567 = t379 + t405;
t547 = qJD(1) * qJD(2);
t565 = qJDD(2) * t428 + t429 * t547;
t395 = qJ(2) * t557;
t564 = t395 + t404;
t563 = t428 ^ 2 + t751;
t556 = qJD(3) * t399;
t555 = qJD(3) * t400;
t382 = qJD(4) * t399;
t542 = -t426 - t662;
t337 = qJ(4) * t529;
t136 = pkin(3) * t457 - qJ(4) * t531 + t337 + t350;
t266 = t399 * t557 + t400 * t554;
t343 = pkin(3) * t530;
t527 = t399 * t551;
t137 = qJ(4) * t266 + qJD(1) * t373 - t343 + t527;
t539 = t429 * t136 + t428 * t137 + t285 * t557;
t537 = -t285 + t592;
t536 = t289 + t591;
t356 = qJ(4) * t606;
t282 = -pkin(3) * t608 + t356;
t535 = t282 * t554 + t286 * t553 + t382;
t534 = t337 + t571;
t532 = t663 * t429;
t526 = -pkin(1) - t651;
t522 = -t554 / 0.2e1;
t521 = t554 / 0.2e1;
t520 = -t553 / 0.2e1;
t519 = t553 / 0.2e1;
t518 = -pkin(4) * t399 - t322;
t516 = -t224 + t613;
t514 = qJD(3) * t584;
t419 = t429 * rSges(5,2);
t250 = t428 * t506 - t419;
t513 = -t250 + t537;
t188 = rSges(6,1) * t304 - rSges(6,2) * t468;
t512 = -t188 + t518;
t508 = -t405 + t527;
t507 = t285 * t554 + t289 * t553 - t552;
t355 = rSges(2,1) * t429 - rSges(2,2) * t428;
t353 = rSges(2,1) * t428 + rSges(2,2) * t429;
t327 = -rSges(4,2) * t399 + t650;
t467 = -t137 - t527 + t593;
t490 = t331 * t322 + t429 * t706 + t565;
t491 = t537 - t599;
t92 = rSges(6,1) * t175 + rSges(6,2) * t715;
t11 = -t249 * t553 + t255 * t188 - t346 * t92 + (t331 * t399 - t429 * t604) * pkin(4) + t491 * qJDD(1) + (t467 - t639) * qJD(1) + t490;
t469 = -qJDD(2) * t429 + qJD(1) * (-pkin(1) * t558 + t564) + qJDD(1) * t354 + t428 * t547;
t454 = qJD(1) * (-t395 + (t428 * t652 - t394) * qJD(1)) + qJDD(1) * t215 + t469;
t442 = qJDD(1) * t289 + t454 + t706 * t428 + (t136 + t350) * qJD(1);
t12 = t442 + t598 * qJDD(1) + (-t330 * t399 - t428 * t604) * pkin(4) + t640 * qJD(1) - t249 * t554 - t254 * t188 - t330 * t322 - t345 * t92;
t502 = t11 * t429 + t12 * t428;
t447 = -t346 * t188 + t518 * t553 + t571;
t41 = qJD(1) * t491 + t447;
t295 = t322 * t554;
t42 = -t188 * t345 - t295 - t342 + (t536 + t598) * qJD(1) + t508;
t497 = -t41 * t429 - t42 * t428;
t492 = -t428 * t88 - t429 * t87;
t163 = rSges(4,1) * t457 - rSges(4,2) * t529 + t572;
t284 = t324 * t428;
t165 = -qJD(3) * t284 + (t327 * t429 + t414) * qJD(1);
t483 = t163 * t429 + t165 * t428;
t476 = t251 * t428 + t253 * t429;
t471 = -pkin(4) * t555 - t249 - t92;
t460 = -qJDD(4) * t400 + t136 * t553 + t137 * t554 + t330 * t285 + t399 * t544;
t458 = t553 * t575 + t571;
t453 = -t389 + t459;
t448 = -t389 + t574;
t368 = rSges(5,3) * t605;
t366 = rSges(5,3) * t606;
t351 = t400 * t550;
t307 = t327 * qJD(3);
t296 = t322 * t558;
t287 = -rSges(5,1) * t607 + t368;
t283 = -rSges(5,1) * t608 + t366;
t267 = t529 - t531;
t265 = t563 * t556;
t177 = qJD(1) * t203 - t405;
t176 = qJD(1) * t570 + t404;
t164 = -t323 * t554 + (t429 * t506 + t416) * qJD(1);
t162 = rSges(5,1) * t457 - rSges(5,3) * t531 + t573;
t119 = t476 * qJD(3);
t86 = qJDD(1) * t269 + qJD(1) * (-qJD(1) * t541 + t568) + t469;
t85 = t570 * qJDD(1) + (-qJD(1) * t269 - t305) * qJD(1) + t565;
t73 = (t250 * t428 + t252 * t429) * qJD(3) + t507;
t72 = -t295 - t405 + (-qJD(3) * t323 + t382) * t428 + (t252 + t536) * qJD(1);
t71 = qJD(1) * t513 + t458;
t43 = t133 * t345 + t135 * t346 + (t320 * t428 + t321 * t429) * qJD(3) + t507;
t40 = qJD(1) * t163 + qJDD(1) * t253 - t307 * t554 - t324 * t330 + t454;
t39 = -t307 * t553 + t324 * t331 + t538 * qJDD(1) + (-t165 + t593) * qJD(1) + t565;
t23 = t250 * t330 + t583 * t331 + (t162 * t429 + t164 * t428) * qJD(3) + t460;
t22 = qJD(1) * t162 + qJDD(1) * t252 + t330 * t575 + t428 * t514 + t442;
t21 = t323 * t331 + t429 * t514 + t513 * qJDD(1) + (-t164 + t467) * qJD(1) + t490;
t10 = t133 * t254 - t135 * t255 + t320 * t330 + t345 * t66 + t346 * t65 + t580 * t331 + (t208 * t429 + t209 * t428) * qJD(3) + t460;
t1 = [t14 * t666 - m(2) * (-g(1) * t353 + g(2) * t355) + (t68 + t58) * t254 / 0.2e1 + (t67 + t57) * t255 / 0.2e1 + (t15 + t18) * t345 / 0.2e1 + (((t80 - t193 + (t225 + t614) * t429 + t765) * t429 + (t697 + t721 - t779) * t428) * qJD(3) + t730) * t519 - (t16 + t14 + t17) * t346 / 0.2e1 + (t756 * qJD(3) + t175 * t185 + t715 * t182 + t304 * t91 + t759 * t399 - t760 * t400 - t468 * t90) * qJD(1) + (-g(1) * (-t133 + t569 - t658) - (t400 * t674 - t383) * t657 + t11 * (-t504 + t569) + t41 * (t342 + t343 - t505 + t567) + t42 * (t534 + t641) + (t42 * t556 * t674 - t11 * t662) * t429 + (t11 * t459 + t41 * (-qJ(4) * t555 - t382)) * t428 + ((t41 * t662 + t42 * t453) * t428 + (t41 * t453 + t42 * t542) * t429) * qJD(1) - (-t599 * qJD(1) - t41 + t447 + t689) * t42 + (-g(2) + t12) * (t428 * t542 + t289 + t364 + t372 + t586)) * m(6) + (-(-qJD(1) * t250 + t458 + t689 - t71) * t72 + t71 * (t343 + t379 - t508) + t72 * (t534 + t573) + (t72 * t399 * t532 + t71 * (t400 * t642 + t645) * t428) * qJD(3) + ((-t72 * t426 + t448 * t71) * t429 + (-t71 * rSges(5,2) + t448 * t72) * t428) * qJD(1) + (t22 - g(2)) * (t515 - t583) + (t21 - g(1)) * (t419 + (t399 * t642 + t400 * t663) * t428 + t569)) * m(5) + (t87 * t567 + t88 * (t404 + t572) + (t324 * t644 - t647) * qJD(3) + ((-t87 * rSges(4,3) + t88 * (-t389 - t650)) * t428 + (t87 * (-t327 - t389) - t88 * t426) * t429) * qJD(1) - (-qJD(1) * t251 + t509 + t699 - t87) * t88 + (t40 - g(2)) * (t253 + t515) + (t39 - g(1)) * (-t251 + t569)) * m(4) + (t176 * t405 + t177 * (t564 + t568) + (t176 * (t526 + t649) * t429 + (t176 * (-rSges(3,3) - qJ(2)) + t177 * t526) * t428) * qJD(1) - (-qJD(1) * t268 - t176 - t333 + t404) * t177 + (t86 - g(2)) * t203 + (t85 - g(1)) * (t526 * t428 + t407 + t566)) * m(3) + (t710 + t713) * t671 + (-t709 + t714) * t670 + (t722 + t725) * t521 + (((t429 * t516 - t697 + t719) * t429 + (t428 * t516 - t206 + t510 + t517 + t720) * t428) * qJD(3) + t727 - t729) * t522 + (-t723 + t724 + t726) * t520 + (Icges(3,2) * t425 ^ 2 + (Icges(3,1) * t424 + 0.2e1 * Icges(3,4) * t425) * t424 + m(2) * (t353 ^ 2 + t355 ^ 2) + Icges(2,3) - t182 * t468 + t185 * t304 - t755) * qJDD(1); (-m(3) - m(4) + t675) * (-g(2) * t429 + t657) + 0.2e1 * (t11 * t665 + t12 * t664) * m(6) + 0.2e1 * (t21 * t665 + t22 * t664) * m(5) + 0.2e1 * (t39 * t665 + t40 * t664) * m(4) + 0.2e1 * (t664 * t86 + t665 * t85) * m(3); t745 * t671 + t744 * t670 + (qJD(1) * t725 + qJD(3) * t707 + qJDD(1) * t713 + t330 * t719 + t331 * t720) * t665 + (qJD(1) * t724 + qJD(3) * t708 + qJDD(1) * t714 + t330 * t718 + t331 * t721) * t664 + (t723 * t429 + t722 * t428 + (-t709 * t428 + t710 * t429) * qJD(1)) * qJD(1) / 0.2e1 + (t710 * t428 + t709 * t429) * qJDD(1) / 0.2e1 + t727 * t558 / 0.2e1 + t726 * t557 / 0.2e1 + ((-t554 * t704 - t705) * t428 + ((t428 * t703 + t712) * qJD(3) + t711) * t429) * t522 + ((t428 * t720 + t429 * t719) * qJD(1) + t707) * t521 + ((t428 * t721 + t429 * t718) * qJD(1) + t708) * t520 + ((-t553 * t703 + t705) * t429 + ((t429 * t704 + t712) * qJD(3) + t711) * t428) * t519 + ((t747 * t399 + t741 * t400) * qJD(3) + (t742 * t399 + t743 * t400) * qJD(1) - t737) * t654 + (-g(1) * (t358 + t168) - g(2) * (t356 - t167) - g(3) * (-t187 - t459) - t674 * t688 + t41 * t296 + t10 * t582 + t43 * t539 + (t12 * t512 + t42 * t471 + t10 * t599 + t43 * t639 + (t41 * t188 + t43 * (-t135 + t580)) * qJD(1)) * t428 + (t11 * t512 + t41 * t471 + t10 * t598 + t43 * t640 + (t42 * t512 + t43 * t599) * qJD(1)) * t429 - t41 * (t351 + t733) - t42 * (t581 + t734) - t43 * (t535 - t728) - (t41 * (t167 - t282) + t42 * (-pkin(4) * t607 + t168)) * qJD(1) - (t497 * t698 + (-t399 * t43 * t563 + t400 * t497) * pkin(4)) * qJD(3)) * m(6) + (-t71 * (t351 + (-t282 - t283) * qJD(1)) - t72 * (qJD(1) * t287 + t581) - t73 * t535 - ((t73 * t287 + t574 * t71) * t429 + (t73 * t283 + t574 * t72) * t428) * qJD(3) + t71 * t296 + t23 * t582 + t73 * t539 + (t21 * t575 + t71 * t584 + t23 * t252 + t73 * t162 + (t73 * t250 + t575 * t72) * qJD(1)) * t429 + (t22 * t575 + t72 * t584 + t23 * t250 + t73 * t164 + (t71 * t323 + t583 * t73) * qJD(1)) * t428 - g(1) * (t358 + t368) - g(2) * (t356 + t366) + g(3) * t574 - (g(1) * t532 + t655 * t663) * t399) * m(5) + (g(1) * t288 + g(2) * t284 - g(3) * t327 - (t284 * t87 - t647) * qJD(1) - (t119 * (-t284 * t428 - t288 * t429) + t492 * t327) * qJD(3) + (qJD(3) * t483 + t251 * t330 - t253 * t331) * t476 + t119 * ((t251 * t429 - t253 * t428) * qJD(1) + t483) + t492 * t307 + (-t39 * t429 - t40 * t428 + (-t429 * t88 + t644) * qJD(1)) * t324) * m(4) - t752; t675 * (-g(3) * t400 + t688) - m(5) * (t265 * t73 + t266 * t72 + t267 * t71) - m(6) * (t265 * t43 + t266 * t42 + t267 * t41) + 0.2e1 * ((t553 * t71 + t554 * t72 - t23) * t677 + (t41 * t553 + t42 * t554 - t10) * t676) * t400 + 0.2e1 * ((qJD(3) * t73 + t21 * t429 + t22 * t428 + t557 * t72 - t558 * t71) * t677 + (qJD(3) * t43 - t41 * t558 + t42 * t557 + t502) * t676) * t399; t737 * t654 + ((t41 * t92 - t10 * t135 + t43 * (-qJD(1) * t133 - t65)) * t429 + (t42 * t92 - t10 * t133 + t43 * (qJD(1) * t135 - t66)) * t428 + ((-t41 * t428 + t42 * t429) * qJD(1) + t502) * t188 - t41 * (-qJD(1) * t167 - t733) - t42 * (-qJD(1) * t168 - t734) - t43 * t728 + g(1) * t168 - g(2) * t167 - g(3) * t187) * m(6) + t752;];
tau = t1;
