% Calculate vector of inverse dynamics joint torques for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:10
% EndTime: 2019-12-31 17:21:03
% DurationCPUTime: 47.58s
% Computational Cost: add. (13282->928), mult. (34908->1212), div. (0->0), fcn. (33684->6), ass. (0->437)
t780 = Icges(4,1) + Icges(5,1);
t794 = Icges(5,4) + Icges(4,5);
t793 = Icges(4,6) - Icges(5,6);
t411 = sin(qJ(3));
t819 = (Icges(4,4) - Icges(5,5)) * t411;
t779 = Icges(4,2) + Icges(5,3);
t818 = Icges(5,2) + Icges(4,3);
t414 = cos(qJ(3));
t817 = -t793 * t411 + t794 * t414;
t816 = t780 * t414 - t819;
t412 = sin(qJ(2));
t415 = cos(qJ(2));
t654 = Icges(4,4) * t414;
t496 = -Icges(4,2) * t411 + t654;
t815 = t412 * t496 - t793 * t415;
t628 = t412 * t414;
t371 = Icges(5,5) * t628;
t631 = t411 * t412;
t792 = Icges(5,3) * t631 + t371 - t815;
t648 = Icges(5,5) * t414;
t493 = Icges(5,3) * t411 + t648;
t746 = (t493 - t496) * t415 - t793 * t412;
t814 = t818 * t412 + t817 * t415;
t784 = t816 * t412 - t794 * t415;
t745 = t794 * t412 + t816 * t415;
t813 = (t779 * t414 + t819) * t412;
t812 = (-t794 * t411 - t793 * t414) * t412;
t416 = cos(qJ(1));
t575 = qJD(3) * t416;
t413 = sin(qJ(1));
t582 = qJD(2) * t413;
t325 = t412 * t575 + t582;
t578 = qJD(3) * t413;
t580 = qJD(2) * t416;
t326 = -t412 * t578 + t580;
t576 = qJD(3) * t415;
t385 = qJD(1) - t576;
t621 = t416 * t411;
t296 = -t413 * t414 + t415 * t621;
t623 = t415 * t416;
t566 = t414 * t623;
t297 = t411 * t413 + t566;
t627 = t412 * t416;
t786 = t817 * t412 - t818 * t415;
t731 = t792 * t296 + t784 * t297 + t786 * t627;
t141 = Icges(4,5) * t297 - Icges(4,6) * t296 + Icges(4,3) * t627;
t144 = Icges(5,4) * t297 + Icges(5,2) * t627 + Icges(5,6) * t296;
t787 = t141 + t144;
t650 = Icges(5,5) * t296;
t150 = Icges(5,1) * t297 + Icges(5,4) * t627 + t650;
t274 = Icges(4,4) * t296;
t153 = Icges(4,1) * t297 + Icges(4,5) * t627 - t274;
t799 = t150 + t153;
t271 = Icges(5,5) * t297;
t138 = Icges(5,6) * t627 + Icges(5,3) * t296 + t271;
t656 = Icges(4,4) * t297;
t147 = -Icges(4,2) * t296 + Icges(4,6) * t627 + t656;
t801 = t138 - t147;
t756 = t801 * t296 + t799 * t297 + t787 * t627;
t624 = t414 * t416;
t626 = t413 * t415;
t294 = t411 * t626 + t624;
t625 = t414 * t415;
t295 = t413 * t625 - t621;
t629 = t412 * t413;
t139 = Icges(4,5) * t295 - Icges(4,6) * t294 + Icges(4,3) * t629;
t142 = Icges(5,4) * t295 + Icges(5,2) * t629 + Icges(5,6) * t294;
t753 = t139 + t142;
t270 = Icges(5,5) * t295;
t137 = -Icges(5,6) * t629 - Icges(5,3) * t294 - t270;
t273 = Icges(4,4) * t295;
t145 = -Icges(4,2) * t294 + Icges(4,6) * t629 + t273;
t754 = t137 + t145;
t269 = Icges(5,5) * t294;
t148 = Icges(5,1) * t295 + Icges(5,4) * t629 + t269;
t272 = Icges(4,4) * t294;
t152 = -Icges(4,1) * t295 - Icges(4,5) * t629 + t272;
t800 = t148 - t152;
t757 = -t296 * t754 + t800 * t297 + t627 * t753;
t763 = t325 * t756 - t757 * t326 + t731 * t385;
t732 = t792 * t294 + t784 * t295 + t786 * t629;
t758 = t801 * t294 + t799 * t295 + t787 * t629;
t759 = -t294 * t754 + t800 * t295 + t629 * t753;
t764 = t325 * t758 - t759 * t326 + t732 * t385;
t555 = t412 * t580;
t130 = qJD(1) * t294 - qJD(3) * t566 + (t555 - t578) * t411;
t476 = t385 * t411;
t585 = qJD(1) * t415;
t533 = -qJD(3) + t585;
t131 = t416 * t476 + (-t413 * t533 - t555) * t414;
t554 = t415 * t580;
t586 = qJD(1) * t413;
t559 = t412 * t586;
t452 = t554 - t559;
t67 = Icges(5,5) * t131 + Icges(5,6) * t452 - Icges(5,3) * t130;
t73 = Icges(4,4) * t131 + Icges(4,2) * t130 + Icges(4,6) * t452;
t807 = t67 - t73;
t556 = t412 * t582;
t577 = qJD(3) * t414;
t584 = qJD(1) * t416;
t132 = t577 * t626 - t414 * t586 + (t415 * t584 - t556 - t575) * t411;
t583 = qJD(2) * t412;
t133 = t533 * t624 + (-t414 * t583 + t476) * t413;
t581 = qJD(2) * t415;
t453 = t412 * t584 + t413 * t581;
t68 = Icges(5,5) * t133 + Icges(5,6) * t453 + Icges(5,3) * t132;
t74 = Icges(4,4) * t133 - Icges(4,2) * t132 + Icges(4,6) * t453;
t806 = t68 - t74;
t69 = Icges(4,5) * t131 + Icges(4,6) * t130 + Icges(4,3) * t452;
t71 = Icges(5,4) * t131 + Icges(5,2) * t452 - Icges(5,6) * t130;
t805 = t69 + t71;
t70 = Icges(4,5) * t133 - Icges(4,6) * t132 + Icges(4,3) * t453;
t72 = Icges(5,4) * t133 + Icges(5,2) * t453 + Icges(5,6) * t132;
t804 = t70 + t72;
t75 = Icges(5,1) * t131 + Icges(5,4) * t452 - Icges(5,5) * t130;
t77 = Icges(4,1) * t131 + Icges(4,4) * t130 + Icges(4,5) * t452;
t803 = t75 + t77;
t76 = Icges(5,1) * t133 + Icges(5,4) * t453 + Icges(5,5) * t132;
t78 = Icges(4,1) * t133 - Icges(4,4) * t132 + Icges(4,5) * t453;
t802 = t76 + t78;
t798 = t746 * qJD(2) + t813 * qJD(3);
t797 = t814 * qJD(2) + t812 * qJD(3);
t307 = (-Icges(4,1) * t411 - t654) * t412;
t579 = qJD(3) * t412;
t796 = (-Icges(5,1) * t411 + t648) * t579 + qJD(3) * t307 + t745 * qJD(2);
t795 = t792 * t411 + t784 * t414;
t268 = qJD(4) * t296;
t755 = rSges(5,3) + qJ(4);
t771 = rSges(5,1) + pkin(3);
t751 = -t294 * t755 - t295 * t771;
t619 = rSges(5,2) * t629 - t751;
t791 = t385 * t619 - t268;
t572 = qJD(1) * qJD(2);
t336 = -qJDD(2) * t416 + t413 * t572;
t571 = qJDD(3) * t412;
t182 = qJD(3) * t453 + t413 * t571 + t336;
t324 = qJD(2) * t579 - qJDD(3) * t415 + qJDD(1);
t360 = pkin(2) * t556;
t390 = pkin(2) * t623;
t206 = pkin(6) * t453 + qJD(1) * t390 - t360;
t408 = t415 * pkin(2);
t718 = t412 * pkin(6) + t408;
t333 = t718 * qJD(2);
t356 = t416 * pkin(1) + t413 * pkin(5);
t334 = qJD(1) * t356;
t353 = pkin(2) * t412 - pkin(6) * t415;
t318 = t718 * t413;
t409 = t416 * pkin(5);
t354 = pkin(1) * t413 - t409;
t596 = -t318 - t354;
t426 = t336 * t353 + (-t206 - t334) * qJD(1) + t596 * qJDD(1) - t333 * t580;
t518 = rSges(5,1) * t414 + rSges(5,3) * t411;
t734 = (-pkin(3) * t414 - qJ(4) * t411 - t518) * t412;
t607 = -rSges(5,2) * t415 - t734;
t574 = qJD(4) * t411;
t364 = t412 * t574;
t403 = t412 * rSges(5,2);
t451 = t411 * t581 + t412 * t577;
t620 = t364 + t451 * qJ(4) + (-t411 * t579 + t414 * t581) * pkin(3) + (-rSges(5,1) * t411 + rSges(5,3) * t414) * t579 + (t415 * t518 + t403) * qJD(2);
t267 = qJD(4) * t294;
t728 = -t755 * t132 - t133 * t771 - t267;
t675 = rSges(5,2) * t453 - t728;
t6 = -qJD(4) * t130 + qJDD(4) * t296 + t182 * t607 - t324 * t619 - t326 * t620 - t385 * t675 + t426;
t790 = -g(1) + t6;
t770 = t754 * t130 + t800 * t131 + t806 * t296 + t802 * t297 + t753 * t452 + t804 * t627;
t769 = -t801 * t130 + t799 * t131 + t807 * t296 + t803 * t297 + t787 * t452 + t805 * t627;
t768 = -t754 * t132 + t800 * t133 + t806 * t294 + t802 * t295 + t753 * t453 + t804 * t629;
t767 = t801 * t132 + t799 * t133 + t807 * t294 + t803 * t295 + t787 * t453 + t805 * t629;
t761 = -t792 * t130 + t784 * t131 + t798 * t296 + t796 * t297 + t786 * t452 + t797 * t627;
t760 = t792 * t132 + t784 * t133 + t798 * t294 + t796 * t295 + t786 * t453 + t797 * t629;
t491 = -t137 * t411 + t148 * t414;
t55 = -t142 * t415 + t412 * t491;
t489 = -t145 * t411 - t152 * t414;
t57 = -t139 * t415 + t412 * t489;
t789 = t55 + t57;
t490 = t138 * t411 + t150 * t414;
t56 = -t144 * t415 + t412 * t490;
t488 = -t147 * t411 + t153 * t414;
t58 = -t141 * t415 + t412 * t488;
t788 = t56 + t58;
t730 = t795 * t412 - t786 * t415;
t785 = -t412 * t493 + t815;
t783 = (-t795 + t814) * t385 + (t786 * t413 + t489 + t491) * t326 + (-t786 * t416 - t488 - t490) * t325;
t720 = -t295 * rSges(4,1) + t294 * rSges(4,2);
t156 = rSges(4,3) * t629 - t720;
t673 = rSges(4,1) * t414;
t521 = -rSges(4,2) * t411 + t673;
t253 = -rSges(4,3) * t415 + t412 * t521;
t782 = t156 * t385 + t253 * t326;
t597 = -t755 * t628 + t631 * t771;
t781 = (t795 * qJD(2) - t797) * t415 + (t796 * t414 + t798 * t411 + (-t784 * t411 + t792 * t414) * qJD(3) + t786 * qJD(2)) * t412;
t777 = -t812 * t385 + (-t794 * t294 - t793 * t295) * t326 + (t794 * t296 + t793 * t297) * t325;
t320 = pkin(6) * t627 + t390;
t532 = t320 + t356;
t478 = qJD(1) * t532 - t353 * t582;
t617 = rSges(5,2) * t627 + t296 * t755 + t297 * t771;
t43 = -t325 * t607 + t617 * t385 + t267 + t478;
t335 = qJDD(2) * t413 + t416 * t572;
t181 = qJD(3) * t452 + t416 * t571 + t335;
t361 = pkin(6) * t554;
t454 = -t413 * t585 - t555;
t205 = pkin(2) * t454 - pkin(6) * t559 + t361;
t396 = pkin(5) * t584;
t595 = qJD(1) * (-pkin(1) * t586 + t396) + qJDD(1) * t356;
t431 = qJD(1) * t205 + qJDD(1) * t320 - t333 * t582 - t335 * t353 + t595;
t727 = rSges(5,2) * t554 - t755 * t130 + t131 * t771 + t268;
t676 = -rSges(5,2) * t559 + t727;
t7 = qJD(4) * t132 + qJDD(4) * t294 - t181 * t607 + t324 * t617 - t325 * t620 + t385 * t676 + t431;
t774 = -g(2) + t7;
t773 = t181 * t756 + t182 * t757 + t324 * t731 + t325 * t769 - t326 * t770 + t385 * t761;
t772 = t181 * t758 + t182 * t759 + t324 * t732 + t325 * t767 - t326 * t768 + t385 * t760;
t16 = (qJD(2) * t491 - t72) * t415 + (qJD(2) * t142 + t411 * t68 + t414 * t76 + (-t137 * t414 - t148 * t411) * qJD(3)) * t412;
t18 = (qJD(2) * t489 - t70) * t415 + (qJD(2) * t139 - t411 * t74 + t414 * t78 + (-t145 * t414 + t152 * t411) * qJD(3)) * t412;
t766 = t16 + t18;
t17 = (qJD(2) * t490 - t71) * t415 + (qJD(2) * t144 + t411 * t67 + t414 * t75 + (t138 * t414 - t150 * t411) * qJD(3)) * t412;
t19 = (qJD(2) * t488 - t69) * t415 + (qJD(2) * t141 - t411 * t73 + t414 * t77 + (-t147 * t414 - t153 * t411) * qJD(3)) * t412;
t765 = t17 + t19;
t762 = t325 * t788 - t326 * t789 + t385 * t730;
t750 = t785 * t413;
t749 = t785 * t416;
t748 = t784 * t413;
t747 = t784 * t416;
t744 = t364 - t333 - t620;
t743 = t409 + t720;
t742 = t783 * t412;
t741 = t325 * t787 - t326 * t753 + t385 * t786;
t338 = qJD(1) * t354;
t557 = t353 * t580;
t740 = -pkin(2) * t555 + qJD(1) * t318 + t338 + t361 + t396 + t557;
t714 = t413 * t789 + t416 * t788;
t739 = t413 * t788 - t416 * t789;
t713 = t413 * t757 + t416 * t756;
t738 = t413 * t756 - t416 * t757;
t712 = t413 * t759 + t416 * t758;
t737 = t413 * t758 - t416 * t759;
t536 = t43 * t607;
t602 = t318 * t582 + t320 * t580;
t41 = t325 * t619 + t326 * t617 + t364 + t602;
t538 = t41 * t619;
t736 = t536 - t538;
t386 = pkin(6) * t626;
t317 = -pkin(2) * t629 + t386;
t389 = pkin(6) * t623;
t319 = -pkin(2) * t627 + t389;
t735 = t416 * t205 + t413 * t206 - t317 * t582 + t318 * t584 - t319 * t580;
t729 = pkin(5) * qJD(1);
t401 = Icges(3,4) * t415;
t497 = -Icges(3,2) * t412 + t401;
t345 = Icges(3,1) * t412 + t401;
t717 = (-t784 + t813) * t385 + (-t295 * t779 + t269 - t272 + t800) * t326 + (t297 * t779 + t274 - t650 - t799) * t325;
t716 = (-Icges(5,1) * t631 + t307 + t371 + t792) * t385 + (t294 * t780 - t270 + t273 + t754) * t326 + (-t296 * t780 + t271 - t656 + t801) * t325;
t715 = t777 * t412;
t707 = t324 * t730 + t385 * t781;
t706 = t412 * (g(1) * t416 + g(2) * t413);
t642 = Icges(3,3) * t416;
t238 = Icges(3,5) * t626 - Icges(3,6) * t629 - t642;
t372 = Icges(3,4) * t629;
t652 = Icges(3,5) * t416;
t250 = Icges(3,1) * t626 - t372 - t652;
t645 = Icges(3,6) * t416;
t244 = Icges(3,4) * t626 - Icges(3,2) * t629 - t645;
t636 = t244 * t412;
t483 = -t250 * t415 + t636;
t96 = -t416 * t238 - t413 * t483;
t342 = Icges(3,5) * t415 - Icges(3,6) * t412;
t341 = Icges(3,5) * t412 + Icges(3,6) * t415;
t455 = qJD(2) * t341;
t657 = Icges(3,4) * t412;
t346 = Icges(3,1) * t415 - t657;
t251 = Icges(3,5) * t413 + t346 * t416;
t245 = Icges(3,6) * t413 + t416 * t497;
t635 = t245 * t412;
t482 = -t251 * t415 + t635;
t705 = -t416 * t455 + (-t342 * t413 + t482 + t642) * qJD(1);
t239 = Icges(3,3) * t413 + t342 * t416;
t588 = qJD(1) * t239;
t704 = qJD(1) * t483 - t413 * t455 + t588;
t343 = Icges(3,2) * t415 + t657;
t479 = t343 * t412 - t345 * t415;
t703 = t479 * qJD(1) + t342 * qJD(2);
t475 = t205 * t580 + t206 * t582 + t335 * t318 - t320 * t336;
t5 = qJD(4) * t451 + qJDD(4) * t631 + t181 * t619 - t182 * t617 + t325 * t675 + t326 * t676 + t475;
t702 = t41 * t676 + t5 * t617;
t701 = t413 * (-t343 * t416 + t251) - t416 * (-Icges(3,2) * t626 + t250 - t372);
t698 = t181 / 0.2e1;
t697 = t182 / 0.2e1;
t696 = t324 / 0.2e1;
t695 = -t325 / 0.2e1;
t694 = t325 / 0.2e1;
t693 = -t326 / 0.2e1;
t692 = t326 / 0.2e1;
t691 = t335 / 0.2e1;
t690 = t336 / 0.2e1;
t689 = -t385 / 0.2e1;
t688 = t385 / 0.2e1;
t683 = -rSges(5,2) - pkin(6);
t682 = -rSges(4,3) - pkin(6);
t681 = g(1) * t413;
t674 = rSges(3,1) * t415;
t672 = t16 * t326;
t671 = t17 * t325;
t670 = t18 * t326;
t669 = t19 * t325;
t402 = t412 * rSges(4,3);
t404 = t413 * rSges(3,3);
t662 = t55 * t182;
t661 = t56 * t181;
t660 = t57 * t182;
t659 = t58 * t181;
t349 = rSges(3,1) * t412 + rSges(3,2) * t415;
t558 = t349 * t580;
t589 = rSges(3,2) * t629 + t416 * rSges(3,3);
t254 = rSges(3,1) * t626 - t589;
t605 = -t254 - t354;
t121 = qJD(1) * t605 - t558;
t639 = t121 * t413;
t638 = t121 * t416;
t257 = rSges(3,1) * t623 - rSges(3,2) * t627 + t404;
t208 = t257 + t356;
t122 = qJD(1) * t208 - t349 * t582;
t316 = t349 * t416;
t637 = t122 * t316;
t633 = t341 * t413;
t632 = t341 * t416;
t630 = t411 * t415;
t312 = (-rSges(4,1) * t411 - rSges(4,2) * t414) * t412;
t174 = qJD(3) * t312 + (t415 * t521 + t402) * qJD(2);
t616 = -t174 - t333;
t615 = -t294 * t771 + t295 * t755;
t614 = -t296 * t771 + t297 * t755;
t376 = rSges(5,2) * t626;
t613 = t413 * t734 + t376;
t383 = rSges(5,2) * t623;
t612 = t416 * t734 + t383;
t611 = -t413 * t238 - t250 * t623;
t610 = t413 * t239 + t251 * t623;
t606 = -t253 - t353;
t604 = -t625 * t771 - t630 * t755 - t403;
t600 = t413 * t318 + t416 * t320;
t594 = -t343 + t346;
t593 = t345 + t497;
t569 = rSges(4,2) * t631;
t592 = rSges(4,3) * t626 + t413 * t569;
t591 = rSges(4,3) * t623 + t416 * t569;
t590 = rSges(3,2) * t559 + rSges(3,3) * t584;
t587 = qJD(1) * t342;
t134 = -t413 * t479 - t632;
t573 = t134 * qJD(1);
t570 = rSges(4,1) * t628;
t565 = t131 * rSges(4,1) + t130 * rSges(4,2) + rSges(4,3) * t554;
t561 = -t353 - t607;
t160 = t297 * rSges(4,1) - t296 * rSges(4,2) + rSges(4,3) * t627;
t560 = -pkin(1) - t408;
t549 = -pkin(1) - t674;
t546 = t584 / 0.2e1;
t544 = -t582 / 0.2e1;
t543 = t582 / 0.2e1;
t542 = -t580 / 0.2e1;
t541 = t580 / 0.2e1;
t443 = qJD(1) * t596 - t557;
t42 = -t326 * t607 + t443 - t791;
t537 = t42 * t607;
t226 = t251 * t626;
t535 = t416 * t239 - t226;
t534 = -t238 + t635;
t524 = qJD(1) * t319 - t582 * t718;
t352 = rSges(2,1) * t416 - rSges(2,2) * t413;
t350 = rSges(2,1) * t413 + rSges(2,2) * t416;
t351 = -rSges(3,2) * t412 + t674;
t523 = rSges(4,1) * t133 - rSges(4,2) * t132;
t129 = t245 * t415 + t251 * t412;
t456 = qJD(2) * t343;
t167 = -t416 * t456 + (-t413 * t497 + t645) * qJD(1);
t457 = qJD(2) * t345;
t171 = -t416 * t457 + (-t346 * t413 + t652) * qJD(1);
t423 = -qJD(2) * t129 - t167 * t412 + t171 * t415 + t588;
t128 = t244 * t415 + t250 * t412;
t168 = qJD(1) * t245 - t413 * t456;
t172 = qJD(1) * t251 - t413 * t457;
t424 = qJD(1) * t238 - qJD(2) * t128 - t168 * t412 + t172 * t415;
t516 = -(t413 * t704 + t424 * t416) * t416 + (t413 * t705 + t423 * t416) * t413;
t515 = -(t424 * t413 - t416 * t704) * t416 + (t423 * t413 - t416 * t705) * t413;
t97 = -t245 * t629 - t535;
t502 = t413 * t97 - t416 * t96;
t98 = -t244 * t627 - t611;
t99 = -t245 * t627 + t610;
t501 = t413 * t99 - t416 * t98;
t492 = -t122 * t413 - t638;
t487 = t156 * t416 - t160 * t413;
t175 = rSges(3,1) * t454 - rSges(3,2) * t554 + t590;
t314 = t349 * t413;
t176 = -qJD(2) * t314 + (t351 * t416 + t404) * qJD(1);
t486 = t175 * t416 + t176 * t413;
t481 = t254 * t413 + t257 * t416;
t480 = t343 * t415 + t345 * t412;
t477 = -pkin(1) - t718;
t256 = rSges(4,1) * t625 - rSges(4,2) * t630 + t402;
t461 = -qJD(1) * t317 - t580 * t718;
t459 = t412 * t683 + t560;
t458 = t412 * t682 + t560;
t450 = t41 * t675 + t5 * t619;
t445 = -t42 * t619 + t43 * t617;
t444 = -t41 * t617 + t537;
t442 = t244 * t416 - t245 * t413;
t432 = (-t412 * t593 + t415 * t594) * qJD(1);
t54 = t156 * t325 + t160 * t326 + t602;
t65 = t443 - t782;
t66 = t160 * t385 - t253 * t325 + t478;
t425 = t54 * t487 + (t413 * t65 - t416 * t66) * t253;
t328 = t497 * qJD(2);
t329 = t346 * qJD(2);
t422 = qJD(1) * t341 - qJD(2) * t480 - t328 * t412 + t329 * t415;
t421 = t444 * t413 - t416 * t736;
t420 = -t412 * t701 + t442 * t415;
t330 = t351 * qJD(2);
t323 = t353 * t586;
t224 = -t416 * t570 + t591;
t222 = -t413 * t570 + t592;
t203 = -rSges(4,1) * t296 - rSges(4,2) * t297;
t198 = -rSges(4,1) * t294 - rSges(4,2) * t295;
t135 = -t416 * t479 + t633;
t120 = t135 * qJD(1);
t113 = t481 * qJD(2);
t84 = rSges(4,3) * t453 + t523;
t82 = -rSges(4,3) * t559 + t565;
t64 = qJD(1) * t175 + qJDD(1) * t257 - t330 * t582 - t335 * t349 + t595;
t63 = -t330 * t580 + t336 * t349 + t605 * qJDD(1) + (-t176 - t334) * qJD(1);
t62 = t422 * t413 - t416 * t703;
t61 = t413 * t703 + t422 * t416;
t60 = -qJD(2) * t482 + t167 * t415 + t171 * t412;
t59 = -qJD(2) * t483 + t168 * t415 + t172 * t412;
t45 = qJD(2) * t501 + t120;
t44 = qJD(2) * t502 + t573;
t26 = t160 * t324 - t174 * t325 - t181 * t253 + t385 * t82 + t431;
t25 = -t156 * t324 - t174 * t326 + t182 * t253 - t385 * t84 + t426;
t20 = t156 * t181 - t160 * t182 + t325 * t84 + t326 * t82 + t475;
t1 = [(t120 + ((t97 - t226 + (t239 + t636) * t416 + t611) * t416 + t610 * t413) * qJD(2)) * t541 + (t536 * t326 + t774 * (t532 + t617) + (t413 * t6 - t681) * t459 + (t727 + (t477 - t403) * t586 + t740 + t791) * t43 + (t360 + t459 * t584 + (t581 * t683 - t729) * t413 + t728 + t43) * t42 + t790 * (t409 + t751)) * m(5) + (t135 + t129) * t691 + (t44 - t573 + ((t416 * t534 - t610 + t99) * t416 + (t413 * t534 + t535 + t98) * t413) * qJD(2)) * t544 + (t134 + t128) * t690 + (t60 + t61) * t543 + (t122 * (t396 + t590) + (t349 * t639 - t637) * qJD(2) + ((-pkin(1) - t351) * t638 + (t121 * (-rSges(3,3) - pkin(5)) + t122 * t549) * t413) * qJD(1) - (-qJD(1) * t254 - t121 - t338 - t558) * t122 + (t64 - g(2)) * t208 + (t63 - g(1)) * (t413 * t549 + t409 + t589)) * m(3) + (t59 + t62 + t45) * t542 + t662 / 0.2e1 + t661 / 0.2e1 + t659 / 0.2e1 + t660 / 0.2e1 + (m(2) * (t350 ^ 2 + t352 ^ 2) + t480 + Icges(2,3)) * qJDD(1) + (t753 * t415 + (t411 * t754 - t414 * t800) * t412 + t789) * t325 * t689 + t707 + ((t360 - t523 + t458 * t584 + (t581 * t682 - t729) * t413) * t65 - g(1) * t743 - t458 * t681 + (t413 * t458 + t743) * t25 + (t65 + t565 + (t477 - t402) * t586 + t740 + t782) * t66 + (t26 - g(2)) * (t532 + t160)) * m(4) + t763 * t692 - m(2) * (-g(1) * t350 + g(2) * t352) + t761 * t694 + (t760 + t763) * t693 + t671 / 0.2e1 - t672 / 0.2e1 - t670 / 0.2e1 + t669 / 0.2e1 + t731 * t698 + t732 * t697 + (-qJD(2) * t479 + t328 * t415 + t329 * t412) * qJD(1); -qJD(1) * ((t412 * t594 + t415 * t593) * qJD(1) + (t442 * t412 + t415 * t701) * qJD(2)) / 0.2e1 + ((qJD(3) * t714 - t783) * t415 + ((t411 * t746 + t414 * t745 + t786) * t385 + (-t411 * t750 + t414 * t748 - t753) * t326 + (t411 * t749 - t414 * t747 + t787) * t325 + t730 * qJD(3)) * t412) * t689 + ((t412 * t732 + t623 * t758) * qJD(3) + ((qJD(3) * t759 + t741) * t415 + t742) * t413 + (t294 * t746 + t295 * t745) * t385 + (-t294 * t750 + t295 * t748) * t326 + (t294 * t749 - t295 * t747) * t325) * t692 + ((t412 * t731 + t626 * t757) * qJD(3) + ((qJD(3) * t756 + t741) * t415 + t742) * t416 + (t296 * t746 + t297 * t745) * t385 + (-t296 * t750 + t297 * t748) * t326 + (t296 * t749 - t297 * t747) * t325) * t695 + ((-t580 * t633 - t587) * t416 + (t432 + (t416 * t632 + t420) * qJD(2)) * t413) * t541 + ((-t582 * t632 + t587) * t413 + (t432 + (t413 * t633 + t420) * qJD(2)) * t416) * t544 + (-g(1) * (t383 + t389) - g(2) * (t376 + t386) - g(3) * (t718 - t604) - (-t411 * t755 - t414 * t771 - pkin(2)) * t706 - (t412 * t445 + t415 * t421) * qJD(3) + t5 * t600 + (t538 * qJD(1) + t6 * t561 + t702) * t416 + (t537 * qJD(1) + t7 * t561 + t450) * t413 + (-t604 * t325 - t612 * t385 + t413 * t744 + t561 * t584 - t524) * t43 + (-t415 * t574 - t612 * t326 - t613 * t325 + (-t320 - t617) * t586 + t735) * t41 + (-t604 * t326 + t613 * t385 + t416 * t744 + t323 - t461) * t42) * m(5) - (t413 * t764 + t416 * t763) * t576 / 0.2e1 + ((qJD(2) * t486 + t254 * t335 - t257 * t336) * t481 + t113 * ((t254 * t416 - t257 * t413) * qJD(1) + t486) + t492 * t330 + (-t64 * t413 - t63 * t416 + (-t122 * t416 + t639) * qJD(1)) * t349 + g(1) * t316 + g(2) * t314 - g(3) * t351 - (t121 * t314 - t637) * qJD(1) - (t113 * (-t314 * t413 - t316 * t416) + t492 * t351) * qJD(2)) * m(3) + ((t96 * t413 + t97 * t416) * qJD(1) + t515) * t542 + ((t98 * t413 + t99 * t416) * qJD(1) + t516) * t543 + t502 * t690 + t501 * t691 + (t65 * t323 + t20 * t600 + (t20 * t160 + t65 * t616 + (qJD(1) * t66 + t25) * t606) * t416 + (t65 * t253 * qJD(1) + t20 * t156 + t26 * t606 + t66 * t616) * t413 - g(1) * (t389 + t591) - g(2) * (t386 + t592) - g(3) * (t256 + t718) - (-pkin(2) - t673) * t706 - t65 * (-t222 * t385 - t256 * t326 + t461) - t66 * (t224 * t385 - t256 * t325 + t524) - ((-t156 * t65 + t160 * t66) * t412 + t425 * t415) * qJD(3) + ((t156 * qJD(1) + t82) * t416 + (t84 + (-t160 - t320) * qJD(1)) * t413 - t222 * t325 - t224 * t326 + t735) * t54) * m(4) + qJD(1) * (t413 * t60 - t416 * t59 + (t128 * t413 + t129 * t416) * qJD(1)) / 0.2e1 - (qJD(1) * t62 + qJD(2) * t515 + qJDD(1) * t134 + t335 * t97 + t336 * t96 + t772) * t416 / 0.2e1 + (qJD(1) * t61 + qJD(2) * t516 + qJDD(1) * t135 + t335 * t99 + t336 * t98 + t773) * t413 / 0.2e1 + (qJD(1) * t713 + t413 * t769 - t416 * t770) * t694 + (qJD(1) * t714 + t413 * t765 - t416 * t766) * t688 + (qJD(1) * t712 + t413 * t767 - t416 * t768) * t693 - t762 * t579 / 0.2e1 + (t45 + t763) * t546 + (t44 + t764) * t586 / 0.2e1 + qJDD(1) * (-t128 * t416 + t129 * t413) / 0.2e1 + t737 * t697 + t738 * t698 + t739 * t696; (t412 * t713 - t415 * t731) * t698 + (t712 * t412 - t415 * t732) * t697 + (t412 * t714 - t415 * t730) * t696 + (t296 * t717 + t297 * t716 - t416 * t715) * t695 + ((qJD(2) * t713 - t761) * t415 + (-qJD(1) * t738 + t731 * qJD(2) + t413 * t770 + t416 * t769) * t412) * t694 + ((qJD(2) * t712 - t760) * t415 + (-qJD(1) * t737 + t732 * qJD(2) + t413 * t768 + t416 * t767) * t412) * t693 + (t294 * t717 + t295 * t716 - t413 * t715) * t692 + (t777 * t415 + (t411 * t717 + t414 * t716) * t412) * t689 + ((qJD(2) * t714 - t781) * t415 + (-qJD(1) * t739 + t730 * qJD(2) + t413 * t766 + t416 * t765) * t412) * t688 - (t661 + t662 + t671 - t672 + t659 + t660 + t669 - t670 + t707) * t415 / 0.2e1 + t772 * t629 / 0.2e1 + t773 * t627 / 0.2e1 + t762 * t583 / 0.2e1 + ((qJD(2) * t421 + t42 * t675 - t43 * t676 + t6 * t619 - t617 * t7) * t415 + (t445 * qJD(2) + (qJD(1) * t444 - t43 * t620 - t607 * t7 + t450) * t416 + (qJD(1) * t736 + t42 * t620 + t6 * t607 - t702) * t413) * t412 - g(1) * t614 - g(2) * t615 + g(3) * t597 - (t295 * t43 + t297 * t42 + t41 * t628) * qJD(4) - (-t42 * t615 + t43 * t614) * t385 - (t41 * t614 + t42 * t597) * t326 - (t41 * t615 + t43 * t597) * t325) * m(5) + (-t65 * (-t198 * t385 - t312 * t326) - t66 * (t203 * t385 - t312 * t325) - t54 * (t198 * t325 + t203 * t326) - g(1) * t203 - g(2) * t198 - g(3) * t312 + (qJD(2) * t425 + t25 * t156 - t26 * t160 + t65 * t84 - t66 * t82) * t415 + (t65 * (-qJD(2) * t156 + t174 * t413) + t66 * (qJD(2) * t160 - t174 * t416) + t20 * t487 + t54 * (-t156 * t586 - t160 * t584 - t413 * t82 + t416 * t84) + (t25 * t413 - t26 * t416 + (t413 * t66 + t416 * t65) * qJD(1)) * t253) * t412) * m(4) + t763 * (t415 * t541 - t559 / 0.2e1) + t764 * (t412 * t546 + t415 * t543); (-t130 * t42 + t132 * t43 + t451 * t41 + (t43 * t325 + t42 * t326 - g(3) + t5) * t631 + (-t41 * t326 - t43 * t385 + t790) * t296 + (-t41 * t325 + t42 * t385 + t774) * t294) * m(5);];
tau = t1;
