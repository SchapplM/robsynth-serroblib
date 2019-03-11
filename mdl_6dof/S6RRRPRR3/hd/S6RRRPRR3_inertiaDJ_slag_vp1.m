% Calculate time derivative of joint inertia matrix for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:57
% EndTime: 2019-03-09 18:11:44
% DurationCPUTime: 29.78s
% Computational Cost: add. (66282->1090), mult. (90171->1447), div. (0->0), fcn. (92732->10), ass. (0->535)
t459 = qJ(2) + qJ(3);
t451 = cos(t459);
t450 = sin(t459);
t682 = Icges(5,5) * t450;
t401 = -Icges(5,3) * t451 + t682;
t687 = Icges(4,4) * t450;
t404 = Icges(4,2) * t451 + t687;
t779 = -t404 + t401;
t681 = Icges(5,5) * t451;
t405 = Icges(5,1) * t450 - t681;
t686 = Icges(4,4) * t451;
t406 = Icges(4,1) * t450 + t686;
t778 = t406 + t405;
t465 = cos(qJ(1));
t721 = sin(qJ(5));
t595 = t465 * t721;
t565 = t451 * t595;
t722 = cos(qJ(5));
t596 = t465 * t722;
t370 = -t450 * t596 + t565;
t501 = -t450 * t721 - t451 * t722;
t371 = t501 * t465;
t462 = sin(qJ(1));
t285 = -Icges(6,4) * t371 - Icges(6,2) * t370 - Icges(6,6) * t462;
t287 = -Icges(6,1) * t371 - Icges(6,4) * t370 - Icges(6,5) * t462;
t597 = t462 * t721;
t598 = t462 * t722;
t368 = -t450 * t598 + t451 * t597;
t369 = t501 * t462;
t283 = -Icges(6,5) * t371 - Icges(6,6) * t370 - Icges(6,3) * t462;
t670 = t283 * t465;
t128 = -t285 * t368 - t287 * t369 + t670;
t599 = t451 * t721;
t600 = t450 * t722;
t396 = t600 - t599;
t456 = qJD(2) + qJD(3);
t741 = t456 - qJD(5);
t470 = t741 * t396;
t625 = qJD(1) * t462;
t232 = -t465 * t470 + t501 * t625;
t460 = sin(qJ(6));
t463 = cos(qJ(6));
t328 = -t371 * t463 - t460 * t462;
t624 = qJD(1) * t465;
t184 = -qJD(6) * t328 - t232 * t460 - t463 * t624;
t327 = t371 * t460 - t462 * t463;
t185 = qJD(6) * t327 + t232 * t463 - t460 * t624;
t564 = qJD(1) * t600;
t587 = qJD(5) * t722;
t654 = t451 * t465;
t231 = -t587 * t654 - t462 * t564 + (qJD(1) * t597 + t456 * t596) * t451 + t741 * t450 * t595;
t107 = Icges(7,5) * t185 + Icges(7,6) * t184 - Icges(7,3) * t231;
t108 = Icges(7,4) * t185 + Icges(7,2) * t184 - Icges(7,6) * t231;
t109 = Icges(7,1) * t185 + Icges(7,4) * t184 - Icges(7,5) * t231;
t234 = -qJD(1) * t371 - t462 * t470;
t326 = -t369 * t463 + t460 * t465;
t186 = -qJD(6) * t326 - t234 * t460 - t463 * t625;
t325 = t369 * t460 + t463 * t465;
t187 = qJD(6) * t325 + t234 * t463 - t460 * t625;
t203 = Icges(7,5) * t328 + Icges(7,6) * t327 + Icges(7,3) * t370;
t205 = Icges(7,4) * t328 + Icges(7,2) * t327 + Icges(7,6) * t370;
t207 = Icges(7,1) * t328 + Icges(7,4) * t327 + Icges(7,5) * t370;
t655 = t451 * t456;
t657 = t450 * t456;
t233 = qJD(1) * t565 - qJD(5) * t369 - t465 * t564 - t597 * t657 - t598 * t655;
t31 = t107 * t368 + t108 * t325 + t109 * t326 + t186 * t205 + t187 * t207 + t203 * t233;
t202 = Icges(7,5) * t326 + Icges(7,6) * t325 + Icges(7,3) * t368;
t204 = Icges(7,4) * t326 + Icges(7,2) * t325 + Icges(7,6) * t368;
t206 = Icges(7,1) * t326 + Icges(7,4) * t325 + Icges(7,5) * t368;
t496 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t233;
t497 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t233;
t498 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t233;
t467 = t186 * t204 + t187 * t206 + t233 * t202 + t325 * t497 + t326 * t498 + t368 * t496;
t96 = t202 * t368 + t204 * t325 + t206 * t326;
t97 = t203 * t368 + t205 * t325 + t207 * t326;
t13 = t31 * t462 - t467 * t465 + (t462 * t96 + t465 * t97) * qJD(1);
t160 = Icges(6,5) * t232 + Icges(6,6) * t231 - Icges(6,3) * t624;
t161 = Icges(6,4) * t232 + Icges(6,2) * t231 - Icges(6,6) * t624;
t162 = Icges(6,1) * t232 + Icges(6,4) * t231 - Icges(6,5) * t624;
t282 = -Icges(6,5) * t369 - Icges(6,6) * t368 + Icges(6,3) * t465;
t284 = -Icges(6,4) * t369 - Icges(6,2) * t368 + Icges(6,6) * t465;
t286 = -Icges(6,1) * t369 - Icges(6,4) * t368 + Icges(6,5) * t465;
t485 = Icges(6,5) * t234 - Icges(6,6) * t233 - Icges(6,3) * t625;
t486 = -Icges(6,4) * t234 + Icges(6,2) * t233 + Icges(6,6) * t625;
t487 = Icges(6,1) * t234 - Icges(6,4) * t233 - Icges(6,5) * t625;
t495 = t282 * t465 - t284 * t368 - t286 * t369;
t671 = t283 * t462;
t771 = (t160 * t465 - t161 * t368 - t162 * t369 - t233 * t285 + t234 * t287) * t462 - (-t233 * t284 + t234 * t286 - t282 * t625 + t368 * t486 - t369 * t487 + t465 * t485) * t465 + (t128 * t465 + (t495 - t671) * t462) * qJD(1) + t13;
t744 = t771 * t465;
t30 = t107 * t370 + t108 * t327 + t109 * t328 + t184 * t205 + t185 * t207 - t203 * t231;
t468 = t184 * t204 + t185 * t206 - t231 * t202 + t327 * t497 + t328 * t498 + t370 * t496;
t98 = t202 * t370 + t204 * t327 + t206 * t328;
t99 = t203 * t370 + t205 * t327 + t207 * t328;
t12 = t30 * t462 - t468 * t465 + (t462 * t98 + t465 * t99) * qJD(1);
t129 = -t285 * t370 - t287 * t371 - t671;
t494 = -t282 * t462 - t284 * t370 - t286 * t371;
t753 = ((-t160 * t462 - t161 * t370 - t162 * t371 + t231 * t285 + t232 * t287) * t462 - (t231 * t284 + t232 * t286 - t282 * t624 + t370 * t486 - t371 * t487 - t462 * t485) * t465 + (t129 * t465 + (t494 - t670) * t462) * qJD(1) + t12) * t462;
t777 = t753 - t744;
t542 = Icges(5,3) * t450 + t681;
t346 = Icges(5,6) * t462 + t465 * t542;
t546 = -Icges(4,2) * t450 + t686;
t352 = Icges(4,6) * t462 + t465 * t546;
t776 = t346 - t352;
t543 = Icges(4,5) * t451 - Icges(4,6) * t450;
t348 = Icges(4,3) * t462 + t465 * t543;
t545 = Icges(5,4) * t451 + Icges(5,6) * t450;
t350 = Icges(5,2) * t462 + t465 * t545;
t775 = t348 + t350;
t549 = Icges(5,1) * t451 + t682;
t354 = Icges(5,4) * t462 + t465 * t549;
t550 = Icges(4,1) * t451 - t687;
t356 = Icges(4,5) * t462 + t465 * t550;
t774 = t354 + t356;
t773 = (-t542 + t546) * t456;
t772 = (t549 + t550) * t456;
t402 = Icges(4,5) * t450 + Icges(4,6) * t451;
t403 = Icges(5,4) * t450 - Icges(5,6) * t451;
t760 = t402 + t403;
t758 = t450 * t779 + t778 * t451;
t345 = -Icges(5,6) * t465 + t462 * t542;
t351 = -Icges(4,6) * t465 + t462 * t546;
t770 = -t345 + t351;
t353 = -Icges(5,4) * t465 + t462 * t549;
t355 = -Icges(4,5) * t465 + t462 * t550;
t769 = t353 + t355;
t768 = -t462 / 0.2e1;
t61 = t462 * t97 - t96 * t465;
t89 = t128 * t462 - t465 * t495;
t767 = t61 + t89;
t62 = t462 * t99 - t98 * t465;
t90 = t129 * t462 - t465 * t494;
t766 = t62 + t90;
t317 = rSges(6,1) * t396 + rSges(6,2) * t501;
t457 = t462 ^ 2;
t458 = t465 ^ 2;
t629 = t457 + t458;
t765 = t317 * t629;
t110 = t185 * rSges(7,1) + t184 * rSges(7,2) - t231 * rSges(7,3);
t764 = -t232 * pkin(5) + t231 * pkin(10) - t110;
t553 = -t326 * rSges(7,1) - t325 * rSges(7,2);
t208 = t368 * rSges(7,3) - t553;
t715 = t369 * pkin(5);
t763 = -t368 * pkin(10) - t208 + t715;
t209 = t328 * rSges(7,1) + t327 * rSges(7,2) + t370 * rSges(7,3);
t648 = -t371 * pkin(5) + t370 * pkin(10) + t209;
t347 = -Icges(4,3) * t465 + t462 * t543;
t349 = -Icges(5,2) * t465 + t462 * t545;
t535 = t345 * t450 + t353 * t451;
t745 = t465 * t535;
t533 = t351 * t450 - t355 * t451;
t746 = t465 * t533;
t762 = -t745 + t746 + (-t347 - t349) * t462;
t532 = t352 * t450 - t356 * t451;
t534 = t346 * t450 + t354 * t451;
t761 = (-t532 + t534) * t465 + t775 * t462;
t255 = -rSges(7,3) * t501 + (rSges(7,1) * t463 - rSges(7,2) * t460) * t396;
t645 = pkin(5) * t396 - pkin(10) * t501 + t255;
t658 = t406 * t456;
t659 = t405 * t456;
t660 = t404 * t456;
t661 = t401 * t456;
t759 = (t661 - t660 + t772) * t451 + (-t659 - t658 - t773) * t450 + t760 * qJD(1);
t278 = -qJD(1) * t353 - t465 * t659;
t280 = -qJD(1) * t355 - t465 * t658;
t757 = t456 * t776 + t278 + t280;
t270 = -qJD(1) * t345 - t465 * t661;
t276 = -qJD(1) * t351 - t465 * t660;
t756 = -t456 * t774 + t270 - t276;
t466 = -pkin(8) - pkin(7);
t461 = sin(qJ(2));
t622 = qJD(2) * t461;
t615 = pkin(2) * t622;
t755 = qJD(1) * t466 + t615;
t754 = (-t543 - t545) * t456 + t758 * qJD(1);
t464 = cos(qJ(2));
t430 = rSges(3,1) * t461 + rSges(3,2) * t464;
t514 = qJD(2) * t430;
t752 = t462 * t514;
t688 = Icges(3,4) * t464;
t548 = -Icges(3,2) * t461 + t688;
t383 = Icges(3,6) * t462 + t465 * t548;
t689 = Icges(3,4) * t461;
t552 = Icges(3,1) * t464 - t689;
t385 = Icges(3,5) * t462 + t465 * t552;
t530 = t383 * t461 - t385 * t464;
t751 = t462 * t530;
t750 = t462 * t532;
t749 = t462 * t534;
t445 = pkin(2) * t464 + pkin(1);
t711 = pkin(1) - t445;
t748 = t462 * t711;
t382 = -Icges(3,6) * t465 + t462 * t548;
t384 = -Icges(3,5) * t465 + t462 * t552;
t531 = t382 * t461 - t384 * t464;
t747 = t465 * t531;
t743 = qJD(1) * t347;
t742 = qJD(1) * t349;
t544 = Icges(3,5) * t464 - Icges(3,6) * t461;
t380 = -Icges(3,3) * t465 + t462 * t544;
t539 = -t204 * t460 + t206 * t463;
t102 = -t202 * t501 + t396 * t539;
t538 = -t205 * t460 + t207 * t463;
t103 = -t203 * t501 + t396 * t538;
t298 = -qJD(5) * t599 - t396 * t456 + t450 * t587;
t252 = -Icges(7,3) * t501 + (Icges(7,5) * t463 - Icges(7,6) * t460) * t396;
t253 = -Icges(7,6) * t501 + (Icges(7,4) * t463 - Icges(7,2) * t460) * t396;
t254 = -Icges(7,5) * t501 + (Icges(7,1) * t463 - Icges(7,4) * t460) * t396;
t120 = t252 * t370 + t253 * t327 + t254 * t328;
t619 = qJD(6) * t396;
t591 = t460 * t619;
t299 = t741 * t501;
t668 = t299 * t463;
t517 = -t591 - t668;
t590 = t463 * t619;
t669 = t299 * t460;
t518 = -t590 + t669;
t151 = Icges(7,5) * t517 + Icges(7,6) * t518 + Icges(7,3) * t298;
t152 = Icges(7,4) * t517 + Icges(7,2) * t518 + Icges(7,6) * t298;
t153 = Icges(7,1) * t517 + Icges(7,4) * t518 + Icges(7,5) * t298;
t36 = t151 * t370 + t152 * t327 + t153 * t328 + t184 * t253 + t185 * t254 - t231 * t252;
t3 = t120 * t298 - t99 * t231 + t98 * t233 + t30 * t370 - t36 * t501 + t368 * t468;
t119 = t252 * t368 + t253 * t325 + t254 * t326;
t37 = t151 * t368 + t152 * t325 + t153 * t326 + t186 * t253 + t187 * t254 + t233 * t252;
t4 = t119 * t298 - t97 * t231 + t96 * t233 + t31 * t370 + t368 * t467 - t37 * t501;
t27 = -t107 * t501 + t203 * t298 - t538 * t299 + (-t108 * t460 + t109 * t463 + (-t205 * t463 - t207 * t460) * qJD(6)) * t396;
t699 = t27 * t462;
t26 = t298 * t202 - t501 * t496 - t539 * t299 + (t463 * t498 - t460 * t497 + (-t204 * t463 - t206 * t460) * qJD(6)) * t396;
t700 = t26 * t465;
t738 = t501 * (-t700 + t699 + (t102 * t462 + t103 * t465) * qJD(1)) / 0.2e1 - t370 * t12 / 0.2e1 - t368 * t13 / 0.2e1 - t298 * (-t102 * t465 + t103 * t462) / 0.2e1 - t233 * t61 / 0.2e1 + t231 * t62 / 0.2e1 + t465 * t4 / 0.2e1 + t3 * t768;
t737 = 2 * m(3);
t736 = 2 * m(4);
t735 = 2 * m(5);
t734 = 2 * m(6);
t733 = 2 * m(7);
t732 = m(5) / 0.2e1;
t731 = m(6) / 0.2e1;
t730 = m(7) / 0.2e1;
t727 = -pkin(3) - pkin(4);
t726 = t462 / 0.2e1;
t725 = -t465 / 0.2e1;
t724 = -rSges(5,1) - pkin(3);
t723 = -rSges(7,3) - pkin(10);
t720 = m(3) * t430;
t409 = rSges(4,1) * t450 + rSges(4,2) * t451;
t719 = m(4) * t409;
t718 = m(6) * t317;
t717 = pkin(2) * t461;
t716 = pkin(5) * t234;
t714 = t462 * pkin(7);
t455 = t465 * pkin(7);
t710 = -pkin(7) - t466;
t709 = -pkin(9) - t466;
t707 = rSges(3,1) * t464;
t706 = rSges(4,1) * t451;
t705 = rSges(5,1) * t450;
t704 = rSges(3,2) * t461;
t703 = rSges(3,3) * t465;
t454 = t462 * rSges(5,2);
t453 = t462 * rSges(3,3);
t452 = t462 * rSges(4,3);
t694 = t462 * rSges(6,3);
t691 = -rSges(5,3) - qJ(4);
t122 = -t252 * t501 + (-t253 * t460 + t254 * t463) * t396;
t499 = t396 * t463 * t153 - t151 * t501 + t298 * t252 + t253 * t669 - t254 * t668;
t651 = t460 * t152;
t690 = t122 * t298 - ((-t651 + (-t253 * t463 - t254 * t460) * qJD(6)) * t396 + t499) * t501;
t674 = qJ(4) * t450;
t223 = -rSges(6,1) * t299 - rSges(6,2) * t298;
t673 = t223 * t462;
t672 = t223 * t465;
t558 = -rSges(4,2) * t450 + t706;
t379 = t558 * t456;
t666 = t379 * t462;
t665 = t382 * t464;
t664 = t383 * t464;
t663 = t384 * t461;
t662 = t385 * t461;
t656 = t450 * t465;
t653 = t456 * t462;
t652 = t456 * t465;
t650 = t465 * t466;
t154 = rSges(7,1) * t517 + rSges(7,2) * t518 + rSges(7,3) * t298;
t649 = -pkin(5) * t299 + pkin(10) * t298 + t154;
t647 = t232 * rSges(6,1) + t231 * rSges(6,2);
t646 = t645 * t625;
t555 = t369 * rSges(6,1) + t368 * rSges(6,2);
t288 = t465 * rSges(6,3) - t555;
t642 = -t371 * rSges(6,1) - t370 * rSges(6,2);
t289 = t642 - t694;
t188 = -t462 * t288 - t465 * t289;
t343 = t455 + t650 - t748;
t434 = t465 * t445;
t344 = -t465 * pkin(1) + t462 * t710 + t434;
t644 = t462 * t343 + t465 * t344;
t336 = qJ(4) * t657 + (pkin(3) * t456 - qJD(4)) * t451;
t557 = rSges(5,1) * t451 + rSges(5,3) * t450;
t378 = t557 * t456;
t643 = -t336 - t378;
t363 = -t465 * rSges(4,3) + t462 * t558;
t365 = rSges(4,1) * t654 - rSges(4,2) * t656 + t452;
t269 = t462 * t363 + t465 * t365;
t364 = rSges(5,1) * t654 + rSges(5,3) * t656 + t454;
t389 = pkin(3) * t654 + qJ(4) * t656;
t641 = -t364 - t389;
t388 = (pkin(3) * t451 + t674) * t462;
t640 = t462 * t388 + t465 * t389;
t438 = pkin(4) * t654;
t400 = -t462 * pkin(9) + t438;
t639 = -t389 - t400;
t407 = pkin(3) * t450 - qJ(4) * t451;
t391 = t407 * t625;
t408 = -rSges(5,3) * t451 + t705;
t638 = t408 * t625 + t391;
t594 = t450 * t625;
t637 = pkin(4) * t594 + t391;
t636 = -t407 - t408;
t610 = t451 * t652;
t620 = qJD(4) * t450;
t635 = qJ(4) * t610 + t465 * t620;
t634 = rSges(5,2) * t624 + rSges(5,3) * t610;
t611 = t450 * t653;
t633 = pkin(4) * t611 + pkin(9) * t625;
t632 = rSges(4,2) * t594 + rSges(4,3) * t624;
t631 = t755 * t462;
t630 = t465 * t707 + t453;
t628 = qJD(1) * t348;
t627 = qJD(1) * t350;
t381 = Icges(3,3) * t462 + t465 * t544;
t626 = qJD(1) * t381;
t621 = qJD(2) * t464;
t618 = -rSges(6,3) + t709;
t617 = pkin(4) * t655;
t616 = t465 * t704;
t614 = pkin(2) * t621;
t613 = t27 / 0.2e1 + t36 / 0.2e1;
t612 = t37 / 0.2e1 + t26 / 0.2e1;
t556 = -rSges(6,1) * t234 + rSges(6,2) * t233;
t607 = t462 * (-rSges(6,3) * t625 - t556) + t465 * (-rSges(6,3) * t624 + t647) + t288 * t624;
t417 = pkin(3) * t611;
t502 = -t450 * t652 - t451 * t625;
t593 = t451 * t624;
t606 = t462 * (pkin(3) * t593 + t462 * t620 - t417 + (t450 * t624 + t451 * t653) * qJ(4)) + t465 * (pkin(3) * t502 - qJ(4) * t594 + t635) + t388 * t624;
t515 = t409 * t456;
t605 = t462 * (-t462 * t515 + (t465 * t558 + t452) * qJD(1)) + t465 * (rSges(4,1) * t502 - rSges(4,2) * t610 + t632) + t363 * t624;
t604 = -t289 + t639;
t603 = t317 * t625 + t637;
t602 = t462 * ((-t465 * t711 - t714) * qJD(1) - t631) + t465 * (-t465 * t615 + (t465 * t710 + t748) * qJD(1)) + t343 * t624;
t601 = t417 + t631;
t592 = t461 * t625;
t589 = -t103 / 0.2e1 - t120 / 0.2e1;
t588 = t119 / 0.2e1 + t102 / 0.2e1;
t585 = t625 / 0.2e1;
t583 = t624 / 0.2e1;
t582 = -t409 - t717;
t581 = -pkin(4) * t450 - t407;
t194 = t645 * t465;
t320 = t636 * t465;
t271 = qJD(1) * t346 - t462 * t661;
t578 = t353 * t456 - t271;
t277 = qJD(1) * t352 - t462 * t660;
t576 = t355 * t456 + t277;
t279 = qJD(1) * t354 - t462 * t659;
t574 = t345 * t456 + t279;
t281 = qJD(1) * t356 - t462 * t658;
t572 = t351 * t456 - t281;
t571 = -t462 * t466 + t434;
t570 = t639 - t648;
t569 = t637 + t646;
t113 = t462 * t763 - t465 * t648;
t362 = -t465 * rSges(5,2) + t462 * t557;
t198 = t462 * t362 + t465 * t364 + t640;
t399 = t462 * t451 * pkin(4) + t465 * pkin(9);
t568 = t462 * t399 + t465 * t400 + t640;
t567 = t434 + t438 + t389;
t563 = -t336 - t617;
t562 = t636 - t717;
t561 = -t317 + t581;
t560 = -t336 - t614;
t559 = -t704 + t707;
t554 = -rSges(7,1) * t187 - rSges(7,2) * t186;
t551 = Icges(3,1) * t461 + t688;
t547 = Icges(3,2) * t464 + t689;
t442 = pkin(2) * t592;
t493 = t560 - t617;
t488 = -t223 + t493;
t124 = t465 * t488 + t442 + t603;
t523 = t581 - t717;
t513 = -t317 + t523;
t240 = t513 * t465;
t125 = qJD(1) * t240 + t462 * t488;
t541 = t124 * t465 + t125 * t462;
t527 = -t223 + t563;
t130 = t465 * t527 + t603;
t246 = t561 * t465;
t131 = qJD(1) * t246 + t462 * t527;
t540 = t130 * t465 + t131 * t462;
t526 = t581 - t645;
t525 = -t378 + t560;
t524 = -pkin(1) - t559;
t311 = t562 * t465;
t522 = -t445 - t558;
t111 = rSges(7,3) * t233 - t554;
t521 = t763 * t624 + t764 * t465 + (-pkin(10) * t233 - t111 - t716) * t462;
t520 = t462 * (-t408 * t653 + (t465 * t557 + t454) * qJD(1)) + t465 * (rSges(5,1) * t502 - rSges(5,3) * t594 + t634) + t362 * t624 + t606;
t519 = t462 * (pkin(4) * t593 - t633) + t465 * (pkin(4) * t502 - pkin(9) * t624) + t399 * t624 + t606;
t144 = t568 - t188;
t516 = t563 - t649;
t509 = t456 * t403;
t508 = t456 * t402;
t506 = qJD(2) * t551;
t505 = qJD(2) * t547;
t504 = qJD(2) * (-Icges(3,5) * t461 - Icges(3,6) * t464);
t181 = t526 * t465;
t503 = t451 * t727 - t445 - t674;
t500 = t523 - t645;
t93 = t568 - t113;
t41 = -t119 * t501 + t368 * t96 + t370 * t97;
t42 = -t120 * t501 + t368 * t98 + t370 * t99;
t492 = t41 * t585 + t42 * t583 - t738;
t491 = t503 * t462;
t490 = t503 * t465;
t489 = t450 * t691 + t451 * t724 - t445;
t484 = t519 + t607;
t177 = t500 * t465;
t483 = t493 - t649;
t482 = rSges(3,2) * t592 + rSges(3,3) * t624 - t465 * t514;
t481 = t489 * t462;
t212 = -t349 * t465 + t462 * t535;
t213 = -t350 * t465 + t749;
t214 = -t347 * t465 - t462 * t533;
t215 = -t348 * t465 - t750;
t272 = -t465 * t508 - t743;
t273 = -t462 * t508 + t628;
t274 = -t465 * t509 - t742;
t275 = -t462 * t509 + t627;
t59 = t61 * t625;
t60 = t62 * t624;
t86 = t89 * t625;
t87 = t90 * t624;
t480 = t59 + t60 + t753 + t86 + t87 + ((-t212 - t214) * t625 + t762 * t624) * t465 + (((t274 + t272) * t462 + (-t749 + t750 - t762) * qJD(1)) * t462 + (t213 + t215) * t625 + t761 * t624 + ((-t273 - t275) * t462 + (t770 * t655 + t769 * t657 - t742 - t743) * t465 + (t757 * t462 + (-t279 - t281) * t465) * t451 + (t756 * t462 + (-t271 + t277) * t465) * t450 + ((t535 - t533 + t775) * t462 + t761) * qJD(1)) * t465) * t462;
t479 = t519 - t521;
t478 = (t657 * t727 - t615) * t465 + t635;
t477 = t465 * t709 + t491;
t476 = -t59 / 0.2e1 - t60 / 0.2e1 - t86 / 0.2e1 - t87 / 0.2e1 - t767 * t625 / 0.2e1 - t766 * t624 / 0.2e1 - t777;
t475 = (-qJ(4) * t655 - t620) * t462 + t601 + t633;
t474 = t465 * t618 + t491;
t46 = (t275 * t465 + (t213 - t745) * qJD(1)) * t465 + (t212 * qJD(1) + (t270 * t450 + t278 * t451 + t346 * t655 - t354 * t657 + t627) * t462 + (-t274 - t574 * t451 + t578 * t450 + (-t349 + t534) * qJD(1)) * t465) * t462;
t47 = (t273 * t465 + (t215 + t746) * qJD(1)) * t465 + (t214 * qJD(1) + (-t276 * t450 + t280 * t451 - t352 * t655 - t356 * t657 + t628) * t462 + (-t272 + t572 * t451 + t576 * t450 + (-t347 - t532) * qJD(1)) * t465) * t462;
t471 = t480 + (-t47 - t46 - t771) * t465;
t314 = Icges(6,5) * t396 + Icges(6,6) * t501;
t315 = Icges(6,4) * t396 + Icges(6,2) * t501;
t316 = Icges(6,1) * t396 + Icges(6,4) * t501;
t165 = t314 * t465 - t315 * t368 - t316 * t369;
t166 = -t314 * t462 - t315 * t370 - t316 * t371;
t174 = t284 * t501 + t286 * t396;
t175 = t285 * t501 + t287 * t396;
t220 = -Icges(6,5) * t299 - Icges(6,6) * t298;
t221 = -Icges(6,4) * t299 - Icges(6,2) * t298;
t222 = -Icges(6,1) * t299 - Icges(6,4) * t298;
t70 = -t220 * t462 - t221 * t370 - t222 * t371 + t231 * t315 + t232 * t316 - t314 * t624;
t71 = t220 * t465 - t221 * t368 - t222 * t369 - t233 * t315 + t234 * t316 - t314 * t625;
t74 = -t298 * t284 - t299 * t286 + t396 * t487 - t486 * t501;
t75 = t161 * t501 + t162 * t396 - t285 * t298 - t287 * t299;
t469 = t699 / 0.2e1 - t700 / 0.2e1 + (t450 * t757 - t451 * t756 - t462 * t754 + t465 * t759 + t36 + t70 + t75) * t726 + (t37 + t71 + t74 + t754 * t465 + t759 * t462 + (t576 + t578) * t451 + (-t572 + t574) * t450) * t725 + (t450 * t769 + t451 * t770 + t758 * t462 - t760 * t465 + t102 + t119 + t165 + t174) * t585 + (t774 * t450 - t451 * t776 + t760 * t462 + t758 * t465 + t103 + t120 + t166 + t175) * t583;
t415 = t559 * qJD(2);
t387 = -t616 + t630;
t386 = t462 * t559 - t703;
t342 = t582 * t465;
t341 = t582 * t462;
t333 = t714 + (pkin(1) - t704) * t465 + t630;
t332 = t462 * t524 + t455 + t703;
t319 = t636 * t462;
t313 = t365 + t571;
t312 = (rSges(4,3) - t466) * t465 + t522 * t462;
t310 = t562 * t462;
t305 = t462 * t504 + t626;
t304 = -qJD(1) * t380 + t465 * t504;
t292 = t752 + ((-rSges(3,3) - pkin(7)) * t462 + t524 * t465) * qJD(1);
t291 = (t455 + (-pkin(1) - t707) * t462) * qJD(1) + t482;
t262 = t571 - t641;
t261 = (rSges(5,2) - t466) * t465 + t481;
t257 = -t409 * t624 - t666 + (-t461 * t624 - t462 * t621) * pkin(2);
t256 = t409 * t625 + t442 + (-t379 - t614) * t465;
t245 = t561 * t462;
t239 = t513 * t462;
t238 = t381 * t462 - t465 * t530;
t237 = t380 * t462 - t747;
t236 = -t381 * t465 - t751;
t235 = -t380 * t465 - t462 * t531;
t229 = t409 * t653 + (t465 * t522 - t452) * qJD(1) + t631;
t228 = (-t445 - t706) * t625 + (-t515 - t755) * t465 + t632;
t211 = qJD(1) * t320 + t462 * t643;
t210 = t465 * t643 + t638;
t197 = t462 * t618 + t567 + t642;
t196 = t474 + t555;
t193 = t645 * t462;
t191 = qJD(1) * t311 + t462 * t525;
t190 = t465 * t525 + t442 + t638;
t189 = t269 + t644;
t180 = t526 * t462;
t179 = (-t620 + (t451 * t691 + t705) * t456) * t462 + (t465 * t489 - t454) * qJD(1) + t601;
t178 = (t657 * t724 - t615) * t465 + (t481 - t650) * qJD(1) + t634 + t635;
t176 = t500 * t462;
t171 = t198 + t644;
t157 = -t365 * t625 + t605;
t150 = t462 * t709 + t567 + t648;
t149 = t368 * t723 + t477 + t553 + t715;
t143 = -t209 * t501 - t255 * t370;
t142 = t208 * t501 + t255 * t368;
t123 = t208 * t370 - t209 * t368;
t121 = t144 + t644;
t118 = (t490 + t694) * qJD(1) + t475 + t556;
t117 = qJD(1) * t474 + t478 + t647;
t112 = (-t344 - t365) * t625 + t602 + t605;
t104 = t625 * t641 + t520;
t92 = qJD(1) * t194 + t462 * t649;
t91 = t465 * t649 - t646;
t88 = t93 + t644;
t81 = t289 * t625 - t607;
t80 = qJD(1) * t181 + t462 * t516;
t79 = t465 * t516 + t569;
t78 = (-t344 + t641) * t625 + t520 + t602;
t77 = qJD(1) * t177 + t462 * t483;
t76 = t465 * t483 + t442 + t569;
t69 = qJD(1) * t490 + t233 * t723 + t475 + t554 - t716;
t68 = qJD(1) * t477 + t478 - t764;
t64 = t604 * t625 + t484;
t52 = t111 * t501 + t154 * t368 - t208 * t298 + t233 * t255;
t51 = -t110 * t501 - t154 * t370 + t209 * t298 + t231 * t255;
t50 = (-t344 + t604) * t625 + t484 + t602;
t43 = -t110 * t368 + t111 * t370 - t208 * t231 - t209 * t233;
t33 = t625 * t648 + t521;
t32 = t570 * t625 + t479;
t23 = t479 + (-t344 + t570) * t625 + t602;
t1 = [-t253 * t590 - t254 * t591 + (t178 * t262 + t179 * t261) * t735 + (t228 * t313 + t229 * t312) * t736 + (t291 * t333 + t292 * t332) * t737 + (t149 * t69 + t150 * t68) * t733 + (t117 * t197 + t118 * t196) * t734 + t501 * t221 - t298 * t315 - t299 * t316 + t499 + t779 * t657 + t778 * t655 + (-t547 + t552) * t622 + (t548 + t551) * t621 + t773 * t451 + t772 * t450 + (t222 - t651) * t396; (t457 / 0.2e1 + t458 / 0.2e1) * t544 * qJD(2) + t469 + m(3) * ((-t291 * t462 - t292 * t465) * t430 + (-t332 * t465 - t333 * t462) * t415) + m(7) * (t149 * t76 + t150 * t77 + t176 * t68 + t177 * t69) + ((t664 / 0.2e1 + t662 / 0.2e1 - t333 * t720) * t465 + (t332 * t720 + t665 / 0.2e1 + t663 / 0.2e1) * t462) * qJD(1) + m(4) * (t228 * t341 + t229 * t342 + t256 * t312 + t257 * t313) + m(5) * (t178 * t310 + t179 * t311 + t190 * t261 + t191 * t262) + m(6) * (t117 * t239 + t118 * t240 + t124 * t196 + t125 * t197) + (-qJD(2) * t530 + (-qJD(1) * t382 - t465 * t505) * t464 + (-qJD(1) * t384 - t465 * t506) * t461) * t726 + (-qJD(2) * t531 + (qJD(1) * t383 - t462 * t505) * t464 + (qJD(1) * t385 - t462 * t506) * t461) * t725; -t465 * ((t465 * t305 + (t236 + t747) * qJD(1)) * t465 + (t235 * qJD(1) + (-t383 * t621 - t385 * t622 + t626) * t462 + (-t304 + (t663 + t665) * qJD(2) - t530 * qJD(1)) * t465) * t462) + (t171 * t78 + t190 * t311 + t191 * t310) * t735 + (t121 * t50 + t124 * t240 + t125 * t239) * t734 + (t176 * t77 + t177 * t76 + t23 * t88) * t733 + (t112 * t189 + t256 * t342 + t257 * t341) * t736 + (-t235 * t465 + t236 * t462) * t625 + (-t237 * t465 + t238 * t462) * t624 + t480 - t465 * t47 - t465 * t46 + t462 * ((t462 * t304 + (t237 + t751) * qJD(1)) * t462 + (t238 * qJD(1) + (t382 * t621 + t384 * t622) * t465 + (-t305 + (-t662 - t664) * qJD(2) + (t381 - t531) * qJD(1)) * t462) * t465) + ((t386 * t462 + t387 * t465) * ((qJD(1) * t386 + t482) * t465 + (-t752 + (-t387 - t616 + t453) * qJD(1)) * t462) + t629 * t430 * t415) * t737 - t744; m(5) * (t178 * t319 + t179 * t320 + t210 * t261 + t211 * t262) + m(6) * (t117 * t245 + t118 * t246 + t130 * t196 + t131 * t197) + m(7) * (t149 * t79 + t150 * t80 + t180 * t68 + t181 * t69) + t469 + (-t228 * t462 - t229 * t465 + (t312 * t462 - t313 * t465) * qJD(1)) * t719 + m(4) * (-t312 * t465 - t313 * t462) * t379; t471 + m(7) * (t176 * t80 + t177 * t79 + t180 * t77 + t181 * t76 + t23 * t93 + t32 * t88) + m(6) * (t121 * t64 + t124 * t246 + t125 * t245 + t130 * t240 + t131 * t239 + t144 * t50) + m(5) * (t104 * t171 + t190 * t320 + t191 * t319 + t198 * t78 + t210 * t311 + t211 * t310) + (-t256 * t465 - t257 * t462 + (-t341 * t465 + t342 * t462) * qJD(1)) * t719 + m(4) * (-t342 * t379 * t465 + t112 * t269 + t157 * t189 - t341 * t666); t471 + (t180 * t80 + t181 * t79 + t32 * t93) * t733 + (t130 * t246 + t131 * t245 + t144 * t64) * t734 + (t104 * t198 + t210 * t320 + t211 * t319) * t735 + (t379 * t409 * t629 + t157 * t269) * t736; 0.2e1 * ((t149 * t465 + t150 * t462) * t730 + (t196 * t465 + t197 * t462) * t731 + (t261 * t465 + t262 * t462) * t732) * t655 + 0.2e1 * ((-t149 * t625 + t150 * t624 + t462 * t68 + t465 * t69) * t730 + (t117 * t462 + t118 * t465 - t196 * t625 + t197 * t624) * t731 + (t178 * t462 + t179 * t465 - t261 * t625 + t262 * t624) * t732) * t450; 0.2e1 * ((t176 * t653 + t177 * t652 - t23) * t730 + (t239 * t653 + t240 * t652 - t50) * t731 + (t310 * t653 + t311 * t652 - t78) * t732) * t451 + 0.2e1 * ((t176 * t624 - t177 * t625 + t456 * t88 + t462 * t77 + t465 * t76) * t730 + (t121 * t456 + t239 * t624 - t240 * t625 + t541) * t731 + (t171 * t456 + t190 * t465 + t191 * t462 + t310 * t624 - t311 * t625) * t732) * t450; 0.2e1 * ((t180 * t653 + t181 * t652 - t32) * t730 + (t245 * t653 + t246 * t652 - t64) * t731 + (t319 * t653 + t320 * t652 - t104) * t732) * t451 + 0.2e1 * ((t180 * t624 - t181 * t625 + t456 * t93 + t462 * t80 + t465 * t79) * t730 + (t144 * t456 + t245 * t624 - t246 * t625 + t540) * t731 + (t198 * t456 + t210 * t465 + t211 * t462 + t319 * t624 - t320 * t625) * t732) * t450; 0.4e1 * (t732 + t731 + t730) * (-0.1e1 + t629) * t450 * t655; m(7) * (t149 * t91 + t150 * t92 + t193 * t68 + t194 * t69) + (m(6) * (t118 * t317 + t196 * t223) + t74 / 0.2e1 + t71 / 0.2e1 + (t197 * t718 - t166 / 0.2e1 - t175 / 0.2e1 + t589) * qJD(1) + t612) * t465 + (m(6) * (t117 * t317 + t197 * t223) - t70 / 0.2e1 - t75 / 0.2e1 + (-t196 * t718 - t174 / 0.2e1 - t165 / 0.2e1 - t588) * qJD(1) - t613) * t462; m(7) * (t113 * t23 + t176 * t92 + t177 * t91 + t193 * t77 + t194 * t76 + t33 * t88) + m(6) * (t121 * t81 + t188 * t50 + t239 * t673 + t240 * t672) + t476 + ((t239 * t465 - t240 * t462) * qJD(1) + t541) * t718; m(6) * (t144 * t81 + t188 * t64 + t245 * t673 + t246 * t672) + t476 + m(7) * (t113 * t32 + t180 * t92 + t181 * t91 + t193 * t80 + t194 * t79 + t33 * t93) + ((t245 * t465 - t246 * t462) * qJD(1) + t540) * t718; 0.2e1 * ((t456 * t765 - t81) * t731 + (t193 * t653 + t194 * t652 - t33) * t730) * t451 + 0.2e1 * ((t188 * t456 + t223 * t629) * t731 + (t113 * t456 + t193 * t624 - t194 * t625 + t462 * t92 + t465 * t91) * t730) * t450; (t188 * t81 + t223 * t765) * t734 + (t113 * t33 + t193 * t92 + t194 * t91) * t733 + (t462 * t767 + t465 * t766) * qJD(1) + t777; m(7) * (t142 * t69 + t143 * t68 + t149 * t52 + t150 * t51) + t613 * t370 + t612 * t368 + t588 * t233 + t589 * t231 + t690; m(7) * (t123 * t23 + t142 * t76 + t143 * t77 + t176 * t51 + t177 * t52 + t43 * t88) + t492; m(7) * (t123 * t32 + t142 * t79 + t143 * t80 + t180 * t51 + t181 * t52 + t43 * t93) + t492; m(7) * ((-t43 + (t142 * t465 + t143 * t462) * t456) * t451 + (t123 * t456 + t462 * t51 + t465 * t52 + (-t142 * t462 + t143 * t465) * qJD(1)) * t450); m(7) * (t113 * t43 + t123 * t33 + t142 * t91 + t143 * t92 + t193 * t51 + t194 * t52) + (t41 * t768 + t42 * t725) * qJD(1) + t738; -t231 * t42 + t370 * t3 + t233 * t41 + t368 * t4 + t298 * (t102 * t368 + t103 * t370 - t122 * t501) - t501 * (t102 * t233 - t103 * t231 + t26 * t368 + t27 * t370 + t690) + (t123 * t43 + t142 * t52 + t143 * t51) * t733;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
