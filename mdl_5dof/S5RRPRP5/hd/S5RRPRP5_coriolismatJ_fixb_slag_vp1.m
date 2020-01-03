% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:02
% EndTime: 2019-12-31 19:54:25
% DurationCPUTime: 13.28s
% Computational Cost: add. (34468->565), mult. (30223->725), div. (0->0), fcn. (27950->8), ass. (0->351)
t482 = qJ(2) + pkin(8);
t462 = qJ(4) + t482;
t457 = cos(t462);
t456 = sin(t462);
t648 = Icges(5,4) * t456;
t394 = Icges(5,1) * t457 - t648;
t487 = sin(qJ(1));
t489 = cos(qJ(1));
t304 = Icges(5,5) * t487 + t394 * t489;
t446 = Icges(6,5) * t456;
t740 = Icges(6,1) * t457 + t446;
t788 = Icges(6,4) * t487 + t489 * t740 + t304;
t447 = Icges(5,4) * t457;
t393 = Icges(5,1) * t456 + t447;
t647 = Icges(6,5) * t457;
t787 = Icges(6,1) * t456 + t393 - t647;
t513 = Icges(6,3) * t457 - t446;
t781 = Icges(5,2) * t457 + t513 + t648;
t398 = pkin(4) * t457 + qJ(5) * t456;
t399 = rSges(6,1) * t457 + rSges(6,3) * t456;
t769 = t398 + t399;
t244 = t769 * t487;
t625 = t456 * t487;
t695 = rSges(6,1) + pkin(4);
t562 = t695 * t625;
t619 = t457 * t487;
t655 = rSges(6,3) + qJ(5);
t618 = t457 * t489;
t624 = t456 * t489;
t237 = t655 * t618 - t624 * t695;
t774 = t237 * t489;
t152 = t774 + (t619 * t655 - t562) * t487;
t460 = sin(t482);
t486 = sin(qJ(2));
t671 = pkin(2) * t486;
t527 = pkin(3) * t460 + t671;
t549 = t655 * t457;
t566 = t456 * t695 - t549;
t501 = t527 + t566;
t212 = t501 * t487;
t214 = t501 * t489;
t247 = t769 * t489;
t479 = t489 * rSges(6,2);
t484 = t489 ^ 2;
t660 = t487 * rSges(6,2);
t141 = t489 * (t399 * t489 + t660) + t484 * t398 + (-t479 + t244) * t487;
t488 = cos(qJ(2));
t480 = t488 * pkin(2);
t459 = t480 + pkin(1);
t461 = cos(t482);
t670 = pkin(3) * t461;
t420 = t459 + t670;
t485 = -qJ(3) - pkin(6);
t557 = -pkin(7) + t485;
t453 = t487 * t557;
t458 = t487 * t485;
t563 = -t487 * t420 - t489 * t557;
t481 = t489 * pkin(6);
t613 = t485 * t489;
t669 = pkin(1) - t459;
t585 = -t487 * (t487 * t669 - t481 - t613) + t489 * (-t487 * pkin(6) - t489 * t669 - t458);
t529 = -t487 * (t487 * t459 + t563 + t613) + t489 * (-t453 + t458 + (t420 - t459) * t489) + t585;
t88 = t141 + t529;
t47 = t88 * t152 + t212 * t244 + t214 * t247;
t307 = rSges(5,1) * t619 - rSges(5,2) * t625 - t489 * rSges(5,3);
t545 = -rSges(5,2) * t624 + t487 * rSges(5,3);
t205 = t487 * t307 + t489 * (rSges(5,1) * t618 + t545);
t107 = t529 + t205;
t661 = rSges(5,1) * t457;
t400 = -rSges(5,2) * t456 + t661;
t397 = rSges(5,1) * t456 + rSges(5,2) * t457;
t368 = t397 * t487;
t370 = t397 * t489;
t743 = t487 * t368 + t489 * t370;
t505 = t397 + t527;
t254 = t505 * t489;
t760 = t505 * t487;
t772 = t487 * t760;
t744 = t254 * t489 + t772;
t64 = -t107 * t743 + t744 * t400;
t786 = -m(5) * t64 - m(6) * t47;
t207 = t487 * t212;
t724 = m(6) / 0.2e1;
t418 = rSges(4,1) * t460 + rSges(4,2) * t461;
t496 = t418 + t671;
t746 = t496 * t489;
t747 = t496 * t487;
t751 = t487 * t747 + t489 * t746;
t766 = -m(5) / 0.2e1;
t552 = (-t214 * t489 - t207) * t724 + t744 * t766 - m(4) * t751 / 0.2e1;
t412 = t489 * t527;
t208 = -t412 + t237;
t209 = (t527 - t549) * t487 + t562;
t257 = -t412 - t370;
t725 = m(5) / 0.2e1;
t726 = m(4) / 0.2e1;
t554 = (-t208 * t489 + t487 * t209) * t724 + (-t257 * t489 + t772) * t725 + t751 * t726;
t18 = t554 - t552;
t785 = t18 * qJD(1);
t243 = t566 * t487;
t246 = t566 * t489;
t483 = t487 ^ 2;
t559 = t483 + t484;
t773 = t397 * t559;
t594 = (-t243 * t487 - t246 * t489) * t724 + t766 * t773;
t236 = -t487 * t549 + t562;
t595 = (t487 * t236 - t774) * t724 + t743 * t725;
t60 = t595 - t594;
t784 = t60 * qJD(1);
t425 = Icges(6,5) * t618;
t294 = Icges(6,6) * t487 + Icges(6,3) * t624 + t425;
t387 = Icges(5,5) * t457 - Icges(5,6) * t456;
t633 = t387 * t489;
t296 = Icges(5,3) * t487 + t633;
t388 = Icges(6,4) * t457 + Icges(6,6) * t456;
t632 = t388 * t489;
t783 = t294 * t624 + t788 * t618 + (Icges(6,2) * t487 + t296 + t632) * t487;
t390 = -Icges(5,2) * t456 + t447;
t782 = t390 + t787;
t780 = (Icges(5,6) - Icges(6,6)) * t457 + (Icges(6,4) + Icges(5,5)) * t456;
t779 = -t781 * t489 + t788;
t301 = -Icges(6,4) * t489 + t487 * t740;
t426 = Icges(5,4) * t625;
t303 = Icges(5,1) * t619 - Icges(5,5) * t489 - t426;
t778 = -Icges(5,2) * t619 - t513 * t487 + t301 + t303 - t426;
t300 = Icges(5,6) * t487 + t390 * t489;
t777 = -Icges(6,1) * t624 - t393 * t489 + t294 - t300 + t425;
t386 = Icges(6,3) * t456 + t647;
t293 = -Icges(6,6) * t489 + t386 * t487;
t299 = Icges(5,4) * t619 - Icges(5,2) * t625 - Icges(5,6) * t489;
t776 = t787 * t487 - t293 + t299;
t775 = t740 + t394;
t771 = -t300 * t624 + t783;
t768 = -Icges(3,5) * t486 - Icges(4,5) * t460 - Icges(3,6) * t488 - Icges(4,6) * t461;
t748 = t487 * (t293 * t456 + t301 * t457);
t767 = t748 + t783;
t699 = t487 / 0.2e1;
t697 = -t489 / 0.2e1;
t762 = t489 / 0.2e1;
t606 = t487 * t388;
t274 = t487 * (-Icges(6,2) * t489 + t606);
t148 = t293 * t624 + t301 * t618 + t274;
t761 = t148 * t489;
t617 = t460 * t487;
t615 = t461 * t487;
t612 = t486 * t487;
t605 = t487 * t488;
t759 = t780 * t487;
t758 = t780 * t489;
t756 = (t775 - t781) * t457 + (t386 - t782) * t456;
t226 = -t307 + t563;
t227 = -t453 + (t420 + t661) * t489 + t545;
t745 = t226 * t489 + t227 * t487;
t499 = t745 * t400;
t735 = -t456 * t655 - t457 * t695;
t195 = t487 * t735 + t479 + t563;
t196 = t660 - t453 + (t420 - t735) * t489;
t596 = -t247 * t195 - t244 * t196;
t667 = (-t208 * t243 - t209 * t246 + t596) * t724 + (-t499 + (-t257 * t487 - t489 * t760) * t397) * t725;
t668 = (-t212 * t237 - t214 * t236 + t596) * t724 + (-t254 * t368 + t370 * t760 - t499) * t725;
t755 = t667 - t668;
t753 = (t777 * t487 + t776 * t489) * t457 + (-t779 * t487 + t778 * t489) * t456;
t649 = Icges(4,4) * t460;
t414 = Icges(4,2) * t461 + t649;
t417 = Icges(4,1) * t461 - t649;
t650 = Icges(3,4) * t486;
t433 = Icges(3,2) * t488 + t650;
t436 = Icges(3,1) * t488 - t650;
t752 = -(t417 / 0.2e1 - t414 / 0.2e1) * t460 - (t436 / 0.2e1 - t433 / 0.2e1) * t486;
t742 = t768 * t489;
t741 = t768 * t487;
t455 = Icges(4,4) * t461;
t415 = -Icges(4,2) * t460 + t455;
t416 = Icges(4,1) * t460 + t455;
t475 = Icges(3,4) * t488;
t434 = -Icges(3,2) * t486 + t475;
t435 = Icges(3,1) * t486 + t475;
t111 = -t205 * t743 + t400 * t773;
t63 = t141 * t152 + t243 * t244 + t246 * t247;
t739 = -m(5) * t111 - m(6) * t63;
t376 = Icges(3,5) * t487 + t436 * t489;
t567 = -t433 * t489 + t376;
t374 = Icges(3,6) * t487 + t434 * t489;
t569 = -t435 * t489 - t374;
t352 = Icges(4,5) * t487 + t417 * t489;
t571 = -t414 * t489 + t352;
t350 = Icges(4,6) * t487 + t415 * t489;
t573 = -t416 * t489 - t350;
t734 = -t567 * t612 + t569 * t605 - t571 * t617 + t573 * t615;
t450 = Icges(3,4) * t612;
t375 = Icges(3,1) * t605 - Icges(3,5) * t489 - t450;
t568 = -Icges(3,2) * t605 + t375 - t450;
t373 = Icges(3,4) * t605 - Icges(3,2) * t612 - Icges(3,6) * t489;
t570 = t435 * t487 + t373;
t443 = Icges(4,4) * t617;
t351 = Icges(4,1) * t615 - Icges(4,5) * t489 - t443;
t572 = -Icges(4,2) * t615 + t351 - t443;
t349 = Icges(4,4) * t615 - Icges(4,2) * t617 - Icges(4,6) * t489;
t574 = t416 * t487 + t349;
t733 = t460 * t572 + t461 * t574 + t486 * t568 + t488 * t570;
t495 = (-t386 / 0.2e1 + t782 / 0.2e1) * t457 + (-t781 / 0.2e1 + t775 / 0.2e1) * t456;
t262 = t304 * t619;
t536 = t296 * t489 - t262;
t147 = -t300 * t625 - t536;
t295 = Icges(5,5) * t619 - Icges(5,6) * t625 - Icges(5,3) * t489;
t590 = -t487 * t295 - t303 * t618;
t150 = -t299 * t624 - t590;
t534 = t300 * t456 - t295;
t638 = t299 * t456;
t497 = (-t150 * t489 + t771 * t487 - t761) * t762 + (-t761 + (t147 - t262 + (t296 + t638) * t489 + t590) * t489 + (-t748 + t767) * t487) * t697 + (((t534 + t295) * t489 - t767 + t771) * t489 + (t148 - t274 + t150 + t536 - (t303 * t457 - t638) * t489 + t147 + t487 * t534) * t487) * t699;
t730 = 0.4e1 * qJD(1);
t729 = 2 * qJD(2);
t727 = 2 * qJD(4);
t383 = t559 * t456;
t592 = -t457 * t207 - t214 * t618;
t56 = t88 * t383 + t592;
t591 = -t243 * t619 - t246 * t618;
t90 = t141 * t383 + t591;
t715 = m(6) * (t90 + t56);
t525 = -t152 * t457 - t244 * t625 - t247 * t624;
t555 = t456 * t88 + t592;
t714 = m(6) * (t525 + t555);
t52 = t141 * t456 + t525 + t591;
t712 = m(6) * t52;
t593 = t195 * t618 + t196 * t619;
t707 = m(6) * ((t208 * t487 + t209 * t489) * t456 + t593);
t706 = m(6) * (-t212 * t624 + t214 * t625 + t593);
t705 = m(6) * ((t236 * t489 + t237 * t487) * t456 + t593);
t704 = m(6) * (-t243 * t624 + t246 * t625 + t593);
t703 = m(6) * (t195 * t209 + t196 * t208);
t702 = m(6) * (t195 * t236 + t196 * t237);
t700 = -t487 / 0.2e1;
t663 = rSges(3,1) * t488;
t551 = pkin(1) + t663;
t560 = rSges(3,2) * t612 + t489 * rSges(3,3);
t289 = -t487 * t551 + t481 + t560;
t611 = t486 * t489;
t452 = rSges(3,2) * t611;
t290 = -t452 + t551 * t489 + (rSges(3,3) + pkin(6)) * t487;
t437 = rSges(3,1) * t486 + rSges(3,2) * t488;
t410 = t437 * t487;
t411 = t437 * t489;
t694 = m(3) * (t289 * t410 - t290 * t411);
t662 = rSges(4,1) * t461;
t547 = t459 + t662;
t561 = rSges(4,2) * t617 + t489 * rSges(4,3);
t249 = -t487 * t547 + t561 - t613;
t616 = t460 * t489;
t546 = -rSges(4,2) * t616 + t487 * rSges(4,3);
t250 = t489 * t547 - t458 + t546;
t693 = m(4) * (t249 * t747 - t250 * t746);
t692 = m(4) * (t249 * t489 + t250 * t487);
t687 = m(5) * (t226 * t760 + t227 * t257);
t686 = m(5) * (t226 * t368 - t227 * t370);
t685 = m(5) * t745;
t678 = m(6) * (t195 * t489 + t196 * t487);
t666 = m(6) * qJD(2);
t665 = m(6) * qJD(4);
t664 = m(6) * qJD(5);
t636 = t349 * t460;
t634 = t373 * t486;
t626 = t456 * t457;
t614 = t461 * t489;
t604 = t488 * t489;
t347 = Icges(4,5) * t615 - Icges(4,6) * t617 - Icges(4,3) * t489;
t588 = -t487 * t347 - t351 * t614;
t516 = Icges(4,5) * t461 - Icges(4,6) * t460;
t348 = Icges(4,3) * t487 + t489 * t516;
t587 = t487 * t348 + t352 * t614;
t371 = Icges(3,5) * t605 - Icges(3,6) * t612 - Icges(3,3) * t489;
t576 = -t487 * t371 - t375 * t604;
t518 = Icges(3,5) * t488 - Icges(3,6) * t486;
t372 = Icges(3,3) * t487 + t489 * t518;
t575 = t487 * t372 + t376 * t604;
t564 = t559 * t626;
t228 = m(6) * t383;
t558 = t228 * qJD(1);
t114 = -t195 * t625 + t196 * t624;
t556 = m(6) * t114 * qJD(1);
t550 = rSges(4,2) * t460 - t480 - t662;
t548 = ((t487 * t759 + t753) * t489 - t758 * t483) * t699 + ((t489 * t758 + t753) * t487 - t759 * t484) * t697;
t269 = t352 * t615;
t535 = t348 * t489 - t269;
t309 = t376 * t605;
t533 = t372 * t489 - t309;
t532 = t350 * t460 - t347;
t531 = t374 * t486 - t371;
t528 = t559 * t671;
t526 = -t480 - t670;
t504 = -t400 + t526;
t500 = t526 - t769;
t498 = t483 * (-t527 + t671) + t489 * (pkin(2) * t611 - t412) - t528;
t491 = -t497 + (t487 * t387 + t777 * t456 + t779 * t457 + t756 * t489 + t606) * t699 + (-t776 * t456 + t778 * t457 + t756 * t487 - t632 - t633) * t697;
t490 = -t495 + ((t293 + t299) * t457 + (-t301 + t303) * t456) * (t699 + t700);
t439 = -rSges(3,2) * t486 + t663;
t343 = t550 * t489;
t341 = t550 * t487;
t255 = t504 * t489;
t253 = t504 * t487;
t234 = t564 - t626;
t218 = t234 * t664;
t215 = t500 * t489;
t213 = t500 * t487;
t187 = -t374 * t611 + t575;
t186 = -t373 * t611 - t576;
t185 = -t374 * t612 - t533;
t164 = t244 * t489 - t487 * t247;
t160 = -t350 * t616 + t587;
t159 = -t349 * t616 - t588;
t158 = -t350 * t617 - t535;
t156 = t164 * t665;
t136 = (-t383 * t457 - t234 + t564) * t664;
t135 = t498 - t743;
t116 = -t186 * t489 + t187 * t487;
t115 = -(-t487 * (-t375 * t488 + t634) - t371 * t489) * t489 + t185 * t487;
t113 = t498 + t152;
t110 = -t159 * t489 + t160 * t487;
t109 = -(-t487 * (-t351 * t461 + t636) - t347 * t489) * t489 + t158 * t487;
t79 = t704 / 0.2e1;
t77 = t705 / 0.2e1;
t73 = t706 / 0.2e1;
t72 = t707 / 0.2e1;
t66 = t678 + t685 + t692;
t61 = t594 + t595;
t50 = t712 / 0.2e1;
t49 = (t185 - t309 + (t372 + t634) * t489 + t576) * t489 + t575 * t487;
t48 = (t489 * t531 + t187 - t575) * t489 + (t487 * t531 + t186 + t533) * t487;
t45 = (t158 - t269 + (t348 + t636) * t489 + t588) * t489 + t587 * t487;
t44 = (t489 * t532 + t160 - t587) * t489 + (t487 * t532 + t159 + t535) * t487;
t43 = t495 + t686 + t702;
t41 = t714 / 0.2e1;
t26 = t715 / 0.2e1;
t23 = t79 - t705 / 0.2e1;
t22 = t79 + t77;
t21 = t77 - t704 / 0.2e1;
t20 = t552 + t554;
t17 = t73 - t707 / 0.2e1;
t16 = t73 + t72;
t15 = t72 - t706 / 0.2e1;
t14 = t703 + t687 + t693 + t694 + (t416 / 0.2e1 + t415 / 0.2e1) * t461 + t495 + (t435 / 0.2e1 + t434 / 0.2e1) * t488 - t752;
t11 = t26 + t50 - t714 / 0.2e1;
t10 = t26 + t41 - t712 / 0.2e1;
t9 = t41 + t50 - t715 / 0.2e1;
t8 = t548 - t739;
t7 = t8 * qJD(4);
t6 = t548 - t786;
t4 = t497 + t755;
t3 = t497 - t755;
t2 = t491 + t667 + t668;
t1 = (t110 / 0.2e1 + t116 / 0.2e1 - t49 / 0.2e1 - t45 / 0.2e1) * t489 + (t44 / 0.2e1 + t48 / 0.2e1 + t115 / 0.2e1 + t109 / 0.2e1) * t487 + t497;
t5 = [t14 * qJD(2) + t66 * qJD(3) + t43 * qJD(4) + t114 * t664, t14 * qJD(1) + t20 * qJD(3) + t2 * qJD(4) + t16 * qJD(5) + (m(3) * ((-t289 * t489 - t290 * t487) * t439 + (-t410 * t489 + t411 * t487) * t437) / 0.2e1 + (t195 * t215 + t196 * t213 - t208 * t212 - t209 * t214) * t724 + (t226 * t255 + t227 * t253 + (-t254 - t257) * t760) * t725 + (t249 * t343 + t250 * t341) * t726) * t729 + (t491 + (t460 * t573 + t461 * t571 + t486 * t569 + t488 * t567) * t699 + (t49 + t45) * t762 + (t518 + t516) * (t484 / 0.2e1 + t483 / 0.2e1) + (t44 + t48 + t115 + t109) * t700 + (-t460 * t574 + t461 * t572 - t486 * t570 + t488 * t568 + t110 + t116) * t697) * qJD(2), qJD(1) * t66 + qJD(2) * t20 + qJD(4) * t61, t43 * qJD(1) + t2 * qJD(2) + t61 * qJD(3) + t491 * qJD(4) + t22 * qJD(5) + ((-t499 + (-t368 * t489 + t370 * t487) * t397) * t725 + (-t236 * t246 - t237 * t243 + t596) * t724) * t727, t16 * qJD(2) + t22 * qJD(4) + t556; (t490 - (t416 + t415) * t461 / 0.2e1 - (t435 + t434) * t488 / 0.2e1 + t752) * qJD(1) + t1 * qJD(2) - t18 * qJD(3) + t3 * qJD(4) + t17 * qJD(5) + (-t703 / 0.4e1 - t687 / 0.4e1 - t693 / 0.4e1 - t694 / 0.4e1) * t730, t1 * qJD(1) + (m(5) * (t107 * t135 - t253 * t760 - t254 * t255) + m(4) * (-t746 * t343 - t747 * t341 + (t487 * (rSges(4,1) * t615 - t561) + t489 * (rSges(4,1) * t614 + t546) + t585) * (-t559 * t418 - t528)) + m(3) * ((t487 * (rSges(3,1) * t605 - t560) + t489 * (rSges(3,1) * t604 + t487 * rSges(3,3) - t452)) * (-t487 * t410 - t411 * t489) + t559 * t439 * t437) + m(6) * (t113 * t88 - t212 * t213 - t214 * t215) + t548 + ((-t741 * t487 + t733 * t489 + t734) * t489 + t742 * t483) * t699 + (((t733 - t742) * t489 + t734) * t487 + t741 * t484) * t697) * qJD(2) + t6 * qJD(4) + t56 * t664, -t785, t3 * qJD(1) + t6 * qJD(2) + t10 * qJD(5) + ((t63 + t47) * t724 + (t64 + t111) * t725) * t727 + (t548 + t739) * qJD(4), t17 * qJD(1) + t10 * qJD(4) + t56 * t666 + t136; t18 * qJD(2) + t60 * qJD(4) - t228 * qJD(5) + (-t685 / 0.4e1 - t678 / 0.4e1 - t692 / 0.4e1) * t730, t785 + t156 + ((-t213 * t489 + t487 * t215) * t724 + (-t253 * t489 + t487 * t255) * t725 + (-t341 * t489 + t487 * t343) * t726) * t729, 0, t164 * t666 + t156 + t784, -t558; t490 * qJD(1) + t4 * qJD(2) - t60 * qJD(3) + t497 * qJD(4) + t23 * qJD(5) + (-t686 / 0.4e1 - t702 / 0.4e1) * t730, t4 * qJD(1) + t7 + t11 * qJD(5) + ((t113 * t141 - t213 * t243 - t215 * t246 + t47) * t724 + (t205 * t135 + (-t253 * t487 - t255 * t489) * t397 + t64) * t725) * t729 + (t548 + t786) * qJD(2), -t784, qJD(1) * t497 + t8 * qJD(2) + t664 * t90 + t7, t23 * qJD(1) + t11 * qJD(2) + t665 * t90 + t136; t15 * qJD(2) + t228 * qJD(3) + t21 * qJD(4) - t556, t15 * qJD(1) + (-t457 * t113 + (t213 * t487 + t215 * t489) * t456 - t56 + t555) * t666 + t9 * qJD(4) + t218, t558, t21 * qJD(1) + t9 * qJD(2) + (t52 - t90) * t665 + t218, 0.4e1 * (qJD(2) / 0.4e1 + qJD(4) / 0.4e1) * t234 * m(6);];
Cq = t5;
