% Calculate time derivative of joint inertia matrix for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:35
% EndTime: 2019-03-09 10:50:26
% DurationCPUTime: 32.23s
% Computational Cost: add. (65377->1327), mult. (95485->1797), div. (0->0), fcn. (99800->10), ass. (0->583)
t449 = sin(qJ(2));
t452 = cos(qJ(2));
t441 = pkin(10) + qJ(4);
t435 = sin(t441);
t436 = cos(t441);
t517 = Icges(5,5) * t436 - Icges(5,6) * t435;
t343 = -Icges(5,3) * t452 + t449 * t517;
t520 = Icges(6,4) * t436 + Icges(6,6) * t435;
t344 = -Icges(6,2) * t452 + t449 * t520;
t739 = t343 + t344;
t738 = Icges(4,3) / 0.2e1;
t450 = sin(qJ(1));
t703 = t450 / 0.2e1;
t453 = cos(qJ(1));
t701 = -t453 / 0.2e1;
t694 = -qJD(1) / 0.2e1;
t675 = Icges(6,5) * t436;
t516 = Icges(6,3) * t435 + t675;
t342 = -Icges(6,6) * t452 + t449 * t516;
t678 = Icges(5,4) * t436;
t521 = -Icges(5,2) * t435 + t678;
t345 = -Icges(5,6) * t452 + t449 * t521;
t676 = Icges(6,5) * t435;
t525 = Icges(6,1) * t436 + t676;
t346 = -Icges(6,4) * t452 + t449 * t525;
t679 = Icges(5,4) * t435;
t526 = Icges(5,1) * t436 - t679;
t347 = -Icges(5,5) * t452 + t449 * t526;
t737 = t739 * t452 + ((-t346 - t347) * t436 + (-t342 + t345) * t435) * t449;
t446 = cos(pkin(10));
t432 = pkin(3) * t446 + pkin(2);
t691 = pkin(2) - t432;
t569 = t691 * t452;
t672 = qJ(3) * t449;
t736 = t569 + t672;
t616 = qJD(4) * t449;
t571 = t436 * t616;
t619 = qJD(2) * t452;
t579 = t435 * t619;
t468 = t571 + t579;
t572 = t435 * t616;
t578 = t436 * t619;
t735 = -t572 + t578;
t574 = t450 * t619;
t622 = qJD(1) * t453;
t470 = t449 * t622 + t574;
t448 = sin(qJ(6));
t451 = cos(qJ(6));
t496 = t435 * t451 - t436 * t448;
t352 = t496 * t449;
t495 = t435 * t448 + t436 * t451;
t353 = t495 * t449;
t245 = Icges(7,5) * t353 + Icges(7,6) * t352 + Icges(7,3) * t452;
t246 = Icges(7,4) * t353 + Icges(7,2) * t352 + Icges(7,6) * t452;
t247 = Icges(7,1) * t353 + Icges(7,4) * t352 + Icges(7,5) * t452;
t657 = t452 * t453;
t368 = t435 * t657 - t450 * t436;
t369 = t450 * t435 + t436 * t657;
t288 = t368 * t451 - t369 * t448;
t289 = t368 * t448 + t369 * t451;
t660 = t449 * t453;
t106 = -t245 * t660 + t288 * t246 + t289 * t247;
t658 = t450 * t452;
t366 = t435 * t658 + t436 * t453;
t614 = qJD(4) * t452;
t570 = t436 * t614;
t618 = qJD(2) * t453;
t575 = t449 * t618;
t615 = qJD(4) * t450;
t241 = qJD(1) * t366 - t453 * t570 + (t575 - t615) * t435;
t494 = t435 * (qJD(1) - t614);
t623 = qJD(1) * t452;
t562 = -qJD(4) + t623;
t242 = t453 * t494 + (-t450 * t562 - t575) * t436;
t141 = -qJD(6) * t289 - t241 * t451 - t242 * t448;
t142 = qJD(6) * t288 - t241 * t448 + t242 * t451;
t667 = t435 * t453;
t367 = t436 * t658 - t667;
t286 = t366 * t451 - t367 * t448;
t287 = t366 * t448 + t367 * t451;
t661 = t449 * t450;
t190 = Icges(7,5) * t287 + Icges(7,6) * t286 - Icges(7,3) * t661;
t192 = Icges(7,4) * t287 + Icges(7,2) * t286 - Icges(7,6) * t661;
t194 = Icges(7,1) * t287 + Icges(7,4) * t286 - Icges(7,5) * t661;
t573 = t452 * t618;
t624 = qJD(1) * t450;
t583 = t449 * t624;
t469 = -t573 + t583;
t620 = qJD(2) * t450;
t577 = t449 * t620;
t243 = -t435 * t577 - qJD(4) * t667 - t436 * t624 + (t435 * t622 + t436 * t615) * t452;
t244 = t450 * t494 + (t453 * t562 - t577) * t436;
t143 = -qJD(6) * t287 + t243 * t451 - t244 * t448;
t144 = qJD(6) * t286 + t243 * t448 + t244 * t451;
t86 = Icges(7,5) * t144 + Icges(7,6) * t143 - Icges(7,3) * t470;
t88 = Icges(7,4) * t144 + Icges(7,2) * t143 - Icges(7,6) * t470;
t90 = Icges(7,1) * t144 + Icges(7,4) * t143 - Icges(7,5) * t470;
t17 = t141 * t192 + t142 * t194 + t190 * t469 + t288 * t88 + t289 * t90 - t660 * t86;
t191 = Icges(7,5) * t289 + Icges(7,6) * t288 - Icges(7,3) * t660;
t193 = Icges(7,4) * t289 + Icges(7,2) * t288 - Icges(7,6) * t660;
t195 = Icges(7,1) * t289 + Icges(7,4) * t288 - Icges(7,5) * t660;
t85 = Icges(7,5) * t142 + Icges(7,6) * t141 + Icges(7,3) * t469;
t87 = Icges(7,4) * t142 + Icges(7,2) * t141 + Icges(7,6) * t469;
t89 = Icges(7,1) * t142 + Icges(7,4) * t141 + Icges(7,5) * t469;
t18 = t141 * t193 + t142 * t195 + t191 * t469 + t288 * t87 + t289 * t89 - t660 * t85;
t725 = qJD(4) - qJD(6);
t224 = t353 * t725 + t496 * t619;
t225 = -t352 * t725 + t495 * t619;
t621 = qJD(2) * t449;
t148 = Icges(7,5) * t225 + Icges(7,6) * t224 - Icges(7,3) * t621;
t149 = Icges(7,4) * t225 + Icges(7,2) * t224 - Icges(7,6) * t621;
t150 = Icges(7,1) * t225 + Icges(7,4) * t224 - Icges(7,5) * t621;
t29 = t141 * t246 + t142 * t247 - t148 * t660 + t288 * t149 + t289 * t150 + t245 * t469;
t74 = -t190 * t660 + t288 * t192 + t289 * t194;
t75 = -t191 * t660 + t288 * t193 + t289 * t195;
t534 = t450 * t74 + t453 * t75;
t1 = -t29 * t452 - t75 * t583 + t534 * t619 + (qJD(2) * t106 + t17 * t450 + (qJD(1) * t74 + t18) * t453) * t449;
t105 = -t245 * t661 + t246 * t286 + t247 * t287;
t19 = t143 * t192 + t144 * t194 - t190 * t470 + t286 * t88 + t287 * t90 - t661 * t86;
t20 = t143 * t193 + t144 * t195 - t191 * t470 + t286 * t87 + t287 * t89 - t661 * t85;
t30 = t143 * t246 + t144 * t247 - t148 * t661 + t286 * t149 + t287 * t150 - t245 * t470;
t72 = -t190 * t661 + t192 * t286 + t194 * t287;
t73 = -t191 * t661 + t193 * t286 + t195 * t287;
t535 = t450 * t72 + t453 * t73;
t2 = -t30 * t452 - t73 * t583 + t535 * t619 + (qJD(2) * t105 + t19 * t450 + (qJD(1) * t72 + t20) * t453) * t449;
t21 = -t190 * t621 + t192 * t224 + t194 * t225 + t352 * t88 + t353 * t90 + t452 * t86;
t22 = -t191 * t621 + t193 * t224 + t195 * t225 + t352 * t87 + t353 * t89 + t452 * t85;
t461 = t452 * t148 + t352 * t149 + t353 * t150 + t224 * t246 + t225 * t247 - t245 * t621;
t35 = t461 * t452;
t83 = t190 * t452 + t192 * t352 + t194 * t353;
t84 = t191 * t452 + t193 * t352 + t195 * t353;
t533 = t450 * t83 + t453 * t84;
t115 = t245 * t452 + t246 * t352 + t247 * t353;
t613 = t115 * qJD(2);
t82 = t84 * t583;
t3 = -t35 - t82 + t533 * t619 + (t613 + t21 * t450 + (qJD(1) * t83 + t22) * t453) * t449;
t31 = -t105 * t452 + t449 * t535;
t32 = -t106 * t452 + t449 * t534;
t34 = -t115 * t452 + t449 * t533;
t734 = -(qJD(1) * (t31 * t453 - t32 * t450) + qJD(2) * t34 + t453 * t1 + t450 * t2) * t449 - (qJD(2) * (t31 * t450 + t32 * t453) - t3) * t452;
t163 = Icges(6,5) * t244 + Icges(6,6) * t470 + Icges(6,3) * t243;
t167 = Icges(6,4) * t244 + Icges(6,2) * t470 + Icges(6,6) * t243;
t171 = Icges(6,1) * t244 + Icges(6,4) * t470 + Icges(6,5) * t243;
t252 = Icges(6,5) * t367 + Icges(6,6) * t661 + Icges(6,3) * t366;
t256 = Icges(6,4) * t367 + Icges(6,2) * t661 + Icges(6,6) * t366;
t260 = Icges(6,1) * t367 + Icges(6,4) * t661 + Icges(6,5) * t366;
t506 = t252 * t435 + t260 * t436;
t50 = (qJD(2) * t506 - t167) * t452 + (qJD(2) * t256 + t163 * t435 + t171 * t436 + (t252 * t436 - t260 * t435) * qJD(4)) * t449;
t165 = Icges(5,5) * t244 - Icges(5,6) * t243 + Icges(5,3) * t470;
t169 = Icges(5,4) * t244 - Icges(5,2) * t243 + Icges(5,6) * t470;
t173 = Icges(5,1) * t244 - Icges(5,4) * t243 + Icges(5,5) * t470;
t254 = Icges(5,5) * t367 - Icges(5,6) * t366 + Icges(5,3) * t661;
t258 = Icges(5,4) * t367 - Icges(5,2) * t366 + Icges(5,6) * t661;
t262 = Icges(5,1) * t367 - Icges(5,4) * t366 + Icges(5,5) * t661;
t504 = -t258 * t435 + t262 * t436;
t52 = (qJD(2) * t504 - t165) * t452 + (qJD(2) * t254 - t169 * t435 + t173 * t436 + (-t258 * t436 - t262 * t435) * qJD(4)) * t449;
t733 = -t50 - t52;
t162 = Icges(6,5) * t242 - Icges(6,6) * t469 - Icges(6,3) * t241;
t166 = Icges(6,4) * t242 - Icges(6,2) * t469 - Icges(6,6) * t241;
t170 = Icges(6,1) * t242 - Icges(6,4) * t469 - Icges(6,5) * t241;
t253 = Icges(6,5) * t369 + Icges(6,6) * t660 + Icges(6,3) * t368;
t257 = Icges(6,4) * t369 + Icges(6,2) * t660 + Icges(6,6) * t368;
t261 = Icges(6,1) * t369 + Icges(6,4) * t660 + Icges(6,5) * t368;
t505 = t253 * t435 + t261 * t436;
t51 = (qJD(2) * t505 - t166) * t452 + (qJD(2) * t257 + t162 * t435 + t170 * t436 + (t253 * t436 - t261 * t435) * qJD(4)) * t449;
t164 = Icges(5,5) * t242 + Icges(5,6) * t241 - Icges(5,3) * t469;
t168 = Icges(5,4) * t242 + Icges(5,2) * t241 - Icges(5,6) * t469;
t172 = Icges(5,1) * t242 + Icges(5,4) * t241 - Icges(5,5) * t469;
t255 = Icges(5,5) * t369 - Icges(5,6) * t368 + Icges(5,3) * t660;
t259 = Icges(5,4) * t369 - Icges(5,2) * t368 + Icges(5,6) * t660;
t263 = Icges(5,1) * t369 - Icges(5,4) * t368 + Icges(5,5) * t660;
t503 = -t259 * t435 + t263 * t436;
t53 = (qJD(2) * t503 - t164) * t452 + (qJD(2) * t255 - t168 * t435 + t172 * t436 + (-t259 * t436 - t263 * t435) * qJD(4)) * t449;
t732 = t51 + t53;
t185 = t342 * t366 + t344 * t661 + t346 * t367;
t186 = t343 * t661 - t345 * t366 + t347 * t367;
t120 = t254 * t661 - t258 * t366 + t262 * t367;
t121 = t255 * t661 - t259 * t366 + t263 * t367;
t512 = t120 * t450 + t121 * t453;
t118 = t252 * t366 + t256 * t661 + t260 * t367;
t119 = t253 * t366 + t257 * t661 + t261 * t367;
t514 = t118 * t450 + t119 * t453;
t731 = (-t185 - t186) * t452 + (t512 + t514) * t449;
t187 = t368 * t342 + t344 * t660 + t369 * t346;
t188 = t343 * t660 - t368 * t345 + t369 * t347;
t124 = t254 * t660 - t368 * t258 + t369 * t262;
t125 = t255 * t660 - t368 * t259 + t369 * t263;
t508 = t124 * t450 + t125 * t453;
t122 = t368 * t252 + t256 * t660 + t369 * t260;
t123 = t368 * t253 + t257 * t660 + t369 * t261;
t510 = t122 * t450 + t123 * t453;
t730 = (-t187 - t188) * t452 + (t508 + t510) * t449;
t269 = (Icges(6,3) * t436 - t676) * t616 + (Icges(6,6) * t449 + t452 * t516) * qJD(2);
t270 = (-Icges(5,5) * t435 - Icges(5,6) * t436) * t616 + (Icges(5,3) * t449 + t452 * t517) * qJD(2);
t271 = (-Icges(6,4) * t435 + Icges(6,6) * t436) * t616 + (Icges(6,2) * t449 + t452 * t520) * qJD(2);
t273 = (-Icges(6,1) * t435 + t675) * t616 + (Icges(6,4) * t449 + t452 * t525) * qJD(2);
t274 = (-Icges(5,1) * t435 - t678) * t616 + (Icges(5,5) * t449 + t452 * t526) * qJD(2);
t666 = t436 * t449;
t668 = t435 * t449;
t729 = t269 * t668 + t468 * t342 - t345 * t571 + t735 * t346 + t347 * t578 + (t273 + t274) * t666 + t739 * t621 + (-t270 - t271) * t452;
t728 = t449 * t691;
t680 = Icges(3,4) * t452;
t524 = -Icges(3,2) * t449 + t680;
t373 = Icges(3,6) * t450 + t453 * t524;
t681 = Icges(3,4) * t449;
t529 = Icges(3,1) * t452 - t681;
t375 = Icges(3,5) * t450 + t453 * t529;
t497 = t373 * t449 - t375 * t452;
t482 = t497 * t450;
t372 = -Icges(3,6) * t453 + t450 * t524;
t374 = -Icges(3,5) * t453 + t450 * t529;
t498 = t372 * t449 - t374 * t452;
t483 = t498 * t453;
t131 = -t256 * t452 + t449 * t506;
t133 = -t254 * t452 + t449 * t504;
t727 = t131 + t133;
t132 = -t257 * t452 + t449 * t505;
t134 = -t255 * t452 + t449 * t503;
t655 = t132 + t134;
t726 = -rSges(3,2) * t660 + t450 * rSges(3,3);
t519 = Icges(3,5) * t452 - Icges(3,6) * t449;
t370 = -Icges(3,3) * t453 + t450 * t519;
t723 = t452 * t655 - t32 - t730;
t721 = t450 * t727 + t453 * t655;
t720 = 2 * m(3);
t719 = 2 * m(4);
t718 = 2 * m(5);
t717 = 2 * m(6);
t716 = 2 * m(7);
t715 = m(4) / 0.2e1;
t714 = m(5) / 0.2e1;
t713 = m(6) / 0.2e1;
t712 = m(7) / 0.2e1;
t711 = -pkin(4) - pkin(5);
t710 = -t106 / 0.2e1;
t445 = sin(pkin(10));
t522 = Icges(4,4) * t446 - Icges(4,2) * t445;
t363 = -Icges(4,6) * t452 + t449 * t522;
t709 = t363 / 0.2e1;
t527 = Icges(4,1) * t446 - Icges(4,4) * t445;
t364 = -Icges(4,5) * t452 + t449 * t527;
t708 = t364 / 0.2e1;
t387 = -t445 * t657 + t450 * t446;
t707 = t387 / 0.2e1;
t659 = t450 * t445;
t388 = t446 * t657 + t659;
t706 = t388 / 0.2e1;
t705 = -t445 / 0.2e1;
t704 = t446 / 0.2e1;
t700 = t453 / 0.2e1;
t699 = -rSges(6,1) - pkin(4);
t698 = rSges(7,3) + pkin(9);
t411 = rSges(3,1) * t449 + rSges(3,2) * t452;
t697 = m(3) * t411;
t696 = pkin(2) * t452;
t438 = t450 * pkin(7);
t695 = t452 * (qJD(1) * t533 - t21 * t453 + t22 * t450);
t693 = qJD(1) / 0.2e1;
t272 = (-Icges(5,2) * t436 - t679) * t616 + (Icges(5,6) * t449 + t452 * t521) * qJD(2);
t690 = (-t345 * t619 + (-qJD(4) * t347 - t272) * t449) * t435 + t729;
t689 = rSges(6,2) * t449;
t688 = rSges(3,3) * t453;
t687 = rSges(5,3) * t449;
t686 = rSges(6,3) * t366;
t685 = t243 * rSges(6,3);
t684 = -rSges(4,3) - qJ(3);
t647 = t242 * pkin(5) + pkin(9) * t583;
t598 = t142 * rSges(7,1) + t141 * rSges(7,2) + rSges(7,3) * t583;
t91 = -rSges(7,3) * t573 + t598;
t683 = -pkin(9) * t573 + t647 + t91;
t540 = t144 * rSges(7,1) + t143 * rSges(7,2);
t92 = -rSges(7,3) * t470 + t540;
t682 = t244 * pkin(5) - pkin(9) * t470 + t92;
t265 = rSges(6,1) * t367 + rSges(6,2) * t661 + t686;
t671 = t265 * t453;
t543 = -rSges(5,1) * t367 + rSges(5,2) * t366;
t266 = rSges(5,3) * t661 - t543;
t670 = t266 * t453;
t669 = t432 * t449;
t665 = t445 * t452;
t664 = t445 * t453;
t663 = t446 * t452;
t447 = -pkin(8) - qJ(3);
t662 = t447 * t452;
t656 = qJ(3) + t447;
t151 = rSges(7,1) * t225 + rSges(7,2) * t224 - rSges(7,3) * t621;
t654 = t735 * pkin(5) - pkin(9) * t621 + t151;
t648 = -t243 * qJ(5) - t366 * qJD(5);
t159 = t244 * pkin(4) - t648;
t355 = t366 * qJ(5);
t294 = pkin(4) * t367 + t355;
t653 = t159 * t660 + t294 * t573;
t158 = t242 * pkin(4) - t241 * qJ(5) + t368 * qJD(5);
t594 = t242 * rSges(6,1) + rSges(6,2) * t573 - t241 * rSges(6,3);
t174 = -rSges(6,2) * t583 + t594;
t652 = -t158 - t174;
t539 = -rSges(7,1) * t287 - rSges(7,2) * t286;
t196 = -rSges(7,3) * t661 - t539;
t651 = pkin(5) * t367 - pkin(9) * t661 + t196;
t643 = t289 * rSges(7,1) + t288 * rSges(7,2);
t197 = -rSges(7,3) * t660 + t643;
t360 = t369 * pkin(5);
t650 = -pkin(9) * t660 + t197 + t360;
t649 = t737 * t621;
t248 = rSges(7,1) * t353 + rSges(7,2) * t352 + rSges(7,3) * t452;
t646 = pkin(5) * t666 + pkin(9) * t452 + t248;
t645 = -t265 - t294;
t267 = t369 * rSges(6,1) + rSges(6,2) * t660 + t368 * rSges(6,3);
t295 = t369 * pkin(4) + t368 * qJ(5);
t644 = -t267 - t295;
t377 = (pkin(4) * t436 + qJ(5) * t435) * t449;
t349 = t377 * t624;
t642 = t295 * t621 + t449 * t349;
t641 = t452 * t294 + t377 * t661;
t303 = t388 * rSges(4,1) + t387 * rSges(4,2) + rSges(4,3) * t660;
t431 = pkin(2) * t657;
t390 = qJ(3) * t660 + t431;
t640 = -t303 - t390;
t600 = t447 * t660;
t424 = pkin(3) * t659;
t632 = t432 * t657 + t424;
t492 = -t600 + t632;
t305 = t492 - t390;
t639 = -t305 - t390;
t538 = t672 + t696;
t384 = qJD(2) * t538 - qJD(3) * t452;
t638 = -(-t449 * t656 - t569) * qJD(2) - t384;
t338 = t452 * t656 - t728;
t410 = pkin(2) * t449 - qJ(3) * t452;
t391 = t410 * t624;
t637 = t338 * t624 + t391;
t636 = -t338 - t410;
t545 = rSges(4,1) * t446 - rSges(4,2) * t445;
t635 = -(rSges(4,3) * t449 + t452 * t545) * qJD(2) - t384;
t365 = -rSges(4,3) * t452 + t449 * t545;
t634 = -t365 - t410;
t389 = t538 * t450;
t633 = t450 * t389 + t453 * t390;
t425 = pkin(3) * t664;
t631 = qJD(1) * t425 + t447 * t583;
t630 = -t447 * t661 - t425;
t629 = rSges(3,2) * t583 + rSges(3,3) * t622;
t617 = qJD(3) * t449;
t423 = t453 * t617;
t434 = pkin(7) * t622;
t628 = t423 + t434;
t627 = t453 * pkin(1) + t438;
t626 = t450 ^ 2 + t453 ^ 2;
t371 = Icges(3,3) * t450 + t453 * t519;
t625 = qJD(1) * t371;
t612 = t35 + t82 / 0.2e1;
t611 = t713 + t712;
t610 = -t447 - t698;
t607 = -t21 / 0.2e1 - t30 / 0.2e1;
t606 = -t22 / 0.2e1 - t29 / 0.2e1;
t605 = -t158 - t683;
t386 = t446 * t658 - t664;
t487 = t445 * t658 + t446 * t453;
t296 = Icges(4,5) * t386 - Icges(4,6) * t487 + Icges(4,3) * t661;
t604 = t296 * t661;
t603 = t296 * t660;
t297 = Icges(4,5) * t388 + Icges(4,6) * t387 + Icges(4,3) * t660;
t602 = t297 * t661;
t601 = t297 * t660;
t599 = -t105 / 0.2e1 - t83 / 0.2e1;
t597 = -t294 - t651;
t596 = -t295 - t650;
t595 = t242 * rSges(5,1) + t241 * rSges(5,2) + rSges(5,3) * t573;
t249 = (pkin(4) * t619 + qJ(5) * t616) * t436 + (qJ(5) * t619 + (-pkin(4) * qJD(4) + qJD(5)) * t449) * t435;
t593 = -t249 + t638;
t542 = rSges(5,1) * t436 - rSges(5,2) * t435;
t277 = (-rSges(5,1) * t435 - rSges(5,2) * t436) * t616 + (t452 * t542 + t687) * qJD(2);
t592 = -t277 + t638;
t413 = qJ(3) * t573;
t420 = pkin(2) * t577;
t471 = -t450 * t623 - t575;
t591 = t450 * (qJ(3) * t470 + qJD(1) * t431 + t450 * t617 - t420) + t453 * (pkin(2) * t471 - qJ(3) * t583 + t413 + t423) + t389 * t622;
t590 = -t295 + t639;
t322 = qJD(1) * t487 + t445 * t575;
t323 = -qJD(1) * t386 - t446 * t575;
t589 = t323 * rSges(4,1) + t322 * rSges(4,2) + rSges(4,3) * t573;
t588 = t349 + t637;
t351 = -rSges(5,3) * t452 + t449 * t542;
t587 = -t351 + t636;
t586 = -t377 + t636;
t268 = t369 * rSges(5,1) - t368 * rSges(5,2) + rSges(5,3) * t660;
t585 = t432 * t577 + t470 * t447;
t439 = t453 * pkin(7);
t584 = t439 - t630;
t541 = rSges(6,1) * t436 + rSges(6,3) * t435;
t350 = -rSges(6,2) * t452 + t449 * t541;
t581 = t350 * t624;
t580 = t351 * t624;
t576 = t449 * t619;
t568 = -t619 / 0.2e1;
t567 = -t432 * t452 - pkin(1);
t566 = t450 * t646;
t565 = t453 * (-t350 - t377);
t564 = t651 * t453;
t313 = t634 * t453;
t563 = qJD(1) * t646;
t561 = t452 * t159 + t249 * t661 + t470 * t377;
t276 = (-rSges(6,1) * t435 + rSges(6,3) * t436) * t616 + (t452 * t541 + t689) * qJD(2);
t560 = -t276 + t593;
t304 = -t450 * t736 + t630;
t559 = t450 * t304 + t453 * t305 + t633;
t558 = -t350 + t586;
t557 = t628 + t631;
t556 = -t355 + t584;
t36 = t73 * t450 - t453 * t72;
t513 = t120 * t453 - t121 * t450;
t515 = t118 * t453 - t119 * t450;
t555 = -t515 / 0.2e1 - t513 / 0.2e1 + t36 / 0.2e1;
t37 = t75 * t450 - t453 * t74;
t509 = t124 * t453 - t125 * t450;
t511 = t122 * t453 - t123 * t450;
t554 = -t509 / 0.2e1 - t511 / 0.2e1 + t37 / 0.2e1;
t553 = -t1 / 0.2e1 + t31 * t694;
t552 = t2 / 0.2e1 + t32 * t694;
t551 = t453 * (-t377 - t646);
t227 = t587 * t453;
t550 = (-pkin(3) * t445 - pkin(7)) * t450;
t549 = t450 * t563;
t548 = rSges(3,1) * t452 - rSges(3,2) * t449;
t324 = qJD(1) * t387 + t445 * t577;
t325 = qJD(1) * t388 - t446 * t577;
t547 = -t325 * rSges(4,1) - t324 * rSges(4,2);
t546 = -rSges(4,1) * t386 + rSges(4,2) * t487;
t544 = t244 * rSges(5,1) - t243 * rSges(5,2);
t532 = t593 - t654;
t531 = t586 - t646;
t530 = t585 + t648;
t528 = Icges(3,1) * t449 + t680;
t523 = Icges(3,2) * t452 + t681;
t518 = Icges(4,5) * t446 - Icges(4,6) * t445;
t507 = -t196 * t453 + t197 * t450;
t502 = -t268 * t450 + t670;
t501 = -t450 * t266 - t268 * t453;
t379 = rSges(3,1) * t657 + t726;
t205 = t558 * t453;
t493 = -pkin(1) - t548;
t491 = t567 - t689;
t490 = t567 - t687;
t489 = t450 * (-qJ(3) * t574 + t420 + (-t453 * t736 + t424) * qJD(1) - t585) + t453 * (-t413 + (-t662 + t728) * t618 + t736 * t624 + t631) + t304 * t622 + t591;
t488 = t450 * t294 + t453 * t295 + t559;
t484 = qJD(2) * t411;
t161 = t531 * t453;
t481 = t450 * t644 + t671;
t480 = -t452 * t727 + t31 + t731;
t479 = qJD(2) * t528;
t478 = qJD(2) * t523;
t477 = qJD(2) * (-Icges(3,5) * t449 - Icges(3,6) * t452);
t476 = t449 * t684 - pkin(1) - t696;
t475 = t449 * t698 + t567;
t473 = t491 * t450;
t472 = t490 * t450;
t467 = t295 + t627 + t632;
t465 = -t187 / 0.2e1 - t188 / 0.2e1 + t710 - t132 / 0.2e1 - t134 / 0.2e1;
t464 = t476 * t450;
t70 = t243 * t342 + t244 * t346 + t366 * t269 + t271 * t661 + t367 * t273 + t344 * t470;
t71 = -t243 * t345 + t244 * t347 + t270 * t661 - t366 * t272 + t367 * t274 + t343 * t470;
t463 = t50 / 0.2e1 + t52 / 0.2e1 + t70 / 0.2e1 + t71 / 0.2e1 - t607;
t68 = -t241 * t342 + t242 * t346 + t368 * t269 + t271 * t660 + t369 * t273 - t344 * t469;
t69 = t241 * t345 + t242 * t347 + t270 * t660 - t368 * t272 + t369 * t274 - t343 * t469;
t462 = t51 / 0.2e1 + t53 / 0.2e1 + t68 / 0.2e1 + t69 / 0.2e1 - t606;
t4 = qJD(1) * t534 - t17 * t453 + t18 * t450;
t5 = qJD(1) * t535 - t19 * t453 + t20 * t450;
t460 = -qJD(2) * (t84 * t450 - t453 * t83) / 0.2e1 - t450 * t5 / 0.2e1 + t4 * t701;
t459 = t453 * t158 + t450 * t159 + t294 * t622 + t489;
t458 = t450 * t596 + t564;
t457 = t84 / 0.2e1 - t465;
t456 = t185 / 0.2e1 + t186 / 0.2e1 + t131 / 0.2e1 + t133 / 0.2e1 - t599;
t454 = (-t662 - t669) * t618 + t557;
t442 = t449 ^ 2;
t397 = t548 * qJD(2);
t394 = t519 * qJD(2);
t378 = t450 * t548 - t688;
t341 = (Icges(4,5) * t449 + t452 * t527) * qJD(2);
t340 = (Icges(4,6) * t449 + t452 * t522) * qJD(2);
t332 = t379 + t627;
t331 = t450 * t493 + t439 + t688;
t312 = t634 * t450;
t307 = t450 * t477 + t625;
t306 = -qJD(1) * t370 + t453 * t477;
t302 = rSges(4,3) * t661 - t546;
t301 = Icges(4,1) * t388 + Icges(4,4) * t387 + Icges(4,5) * t660;
t300 = Icges(4,1) * t386 - Icges(4,4) * t487 + Icges(4,5) * t661;
t299 = Icges(4,4) * t388 + Icges(4,2) * t387 + Icges(4,6) * t660;
t298 = Icges(4,4) * t386 - Icges(4,2) * t487 + Icges(4,6) * t661;
t275 = t294 * t660;
t251 = t411 * t620 + ((-rSges(3,3) - pkin(7)) * t450 + t493 * t453) * qJD(1);
t250 = rSges(3,1) * t471 - rSges(3,2) * t573 - pkin(1) * t624 + t434 + t629;
t229 = t627 - t640;
t228 = t439 + t464 + t546;
t226 = t587 * t450;
t223 = t450 * t371 - t453 * t497;
t222 = t450 * t370 - t483;
t221 = -t371 * t453 - t482;
t220 = -t370 * t453 - t450 * t498;
t217 = Icges(4,1) * t325 + Icges(4,4) * t324 + Icges(4,5) * t470;
t216 = Icges(4,1) * t323 + Icges(4,4) * t322 - Icges(4,5) * t469;
t215 = Icges(4,4) * t325 + Icges(4,2) * t324 + Icges(4,6) * t470;
t214 = Icges(4,4) * t323 + Icges(4,2) * t322 - Icges(4,6) * t469;
t211 = -t452 * t268 - t351 * t660;
t210 = t266 * t452 + t351 * t661;
t209 = qJD(1) * t313 + t450 * t635;
t208 = t365 * t624 + t453 * t635 + t391;
t207 = t492 + t268 + t627;
t206 = t472 + t543 + t584;
t204 = t558 * t450;
t189 = t502 * t449;
t182 = t450 * t302 + t303 * t453 + t633;
t181 = t267 + t467 - t600;
t180 = t367 * t699 + t473 + t556 - t686;
t179 = t420 + (t619 * t684 - t617) * t450 + (t453 * t476 - t438) * qJD(1) + t547;
t178 = -pkin(2) * t575 + qJD(1) * t464 + t413 + t589 + t628;
t177 = rSges(5,3) * t470 + t544;
t176 = t244 * rSges(6,1) + rSges(6,2) * t470 + t685;
t175 = -rSges(5,3) * t583 + t595;
t160 = t531 * t450;
t153 = t449 * t565 + t452 * t644;
t152 = t265 * t452 + t350 * t661 + t641;
t146 = t452 * t197 + t248 * t660;
t145 = -t196 * t452 - t248 * t661;
t140 = t387 * t299 + t388 * t301 + t601;
t139 = t387 * t298 + t388 * t300 + t603;
t138 = -t299 * t487 + t301 * t386 + t602;
t137 = -t298 * t487 + t300 * t386 + t604;
t130 = qJD(1) * t227 + t450 * t592;
t129 = t453 * t592 + t580 + t637;
t117 = (-rSges(5,3) * t619 - t617) * t450 + (t453 * t490 + t550) * qJD(1) - t544 + t585;
t116 = qJD(1) * t472 + t454 + t595;
t114 = t610 * t660 + t360 + t467 + t643;
t113 = t367 * t711 + t450 * t475 + t539 + t556;
t112 = t449 * t481 + t275;
t110 = t507 * t449;
t109 = -t501 + t559;
t108 = qJD(1) * t205 + t450 * t560;
t107 = t453 * t560 + t581 + t588;
t101 = t449 * t551 + t452 * t596;
t100 = t449 * t566 + t452 * t651 + t641;
t99 = (t351 * t620 + t177) * t452 + (-qJD(2) * t266 + t450 * t277 + t351 * t622) * t449;
t98 = (-t351 * t618 - t175) * t452 + (qJD(2) * t268 - t277 * t453 + t580) * t449;
t97 = t450 * t265 + t267 * t453 + t488;
t94 = -t685 + (-rSges(6,2) * t619 - t617) * t450 + t699 * t244 + (t453 * t491 + t550) * qJD(1) + t530;
t93 = qJD(1) * t473 + t158 + t454 + t594;
t81 = t450 * (rSges(4,3) * t574 - t547) + t453 * t589 + (t453 * t302 + t450 * t640) * qJD(1) + t591;
t80 = t449 * t458 + t275;
t65 = qJD(1) * t161 + t450 * t532;
t64 = t453 * t532 + t549 + t588;
t63 = t450 * t651 + t453 * t650 + t488;
t62 = t502 * t619 + (qJD(1) * t501 - t175 * t450 + t177 * t453) * t449;
t57 = (t350 * t620 + t176) * t452 + (qJD(2) * t645 + t450 * t276 + t350 * t622) * t449 + t561;
t56 = (qJD(2) * t565 + t652) * t452 + (t581 + qJD(2) * t267 + (-t249 - t276) * t453) * t449 + t642;
t55 = t711 * t244 + (t619 * t698 - t617) * t450 + (t453 * t475 + t550) * qJD(1) + t530 - t540;
t54 = t567 * t624 + (t452 * t610 - t669) * t618 + t158 + t557 + t598 + t647;
t48 = (-t248 * t620 - t92) * t452 + (qJD(2) * t196 - t450 * t151 - t248 * t622) * t449;
t47 = (t248 * t618 + t91) * t452 + (-qJD(2) * t197 + t151 * t453 - t248 * t624) * t449;
t46 = t164 * t661 - t366 * t168 + t367 * t172 - t243 * t259 + t244 * t263 + t255 * t470;
t45 = t165 * t661 - t366 * t169 + t367 * t173 - t243 * t258 + t244 * t262 + t254 * t470;
t44 = t366 * t162 + t166 * t661 + t367 * t170 + t243 * t253 + t244 * t261 + t257 * t470;
t43 = t366 * t163 + t167 * t661 + t367 * t171 + t243 * t252 + t244 * t260 + t256 * t470;
t42 = t164 * t660 - t368 * t168 + t369 * t172 + t241 * t259 + t242 * t263 - t255 * t469;
t41 = t165 * t660 - t368 * t169 + t369 * t173 + t241 * t258 + t242 * t262 - t254 * t469;
t40 = t368 * t162 + t166 * t660 + t369 * t170 - t241 * t253 + t242 * t261 - t257 * t469;
t39 = t368 * t163 + t167 * t660 + t369 * t171 - t241 * t252 + t242 * t260 - t256 * t469;
t38 = t175 * t453 + t450 * t177 + (t670 + (-t268 + t639) * t450) * qJD(1) + t489;
t33 = t481 * t619 + (t176 * t453 + t652 * t450 + (t450 * t645 + t453 * t644) * qJD(1)) * t449 + t653;
t26 = t174 * t453 + t450 * t176 + (t671 + (-t267 + t590) * t450) * qJD(1) + t459;
t25 = t507 * t619 + (t450 * t91 - t453 * t92 + (t450 * t196 + t197 * t453) * qJD(1)) * t449;
t24 = (qJD(2) * t566 + t682) * t452 + (qJD(2) * t597 + t450 * t654 + t453 * t563) * t449 + t561;
t23 = (qJD(2) * t551 + t605) * t452 + (t650 * qJD(2) + t549 + (-t249 - t654) * t453) * t449 + t642;
t16 = t683 * t453 + t682 * t450 + (t564 + (t590 - t650) * t450) * qJD(1) + t459;
t15 = t458 * t619 + (t682 * t453 + t605 * t450 + (t450 * t597 + t453 * t596) * qJD(1)) * t449 + t653;
t14 = qJD(1) * t512 - t45 * t453 + t46 * t450;
t13 = qJD(1) * t514 - t43 * t453 + t44 * t450;
t12 = qJD(1) * t508 - t41 * t453 + t42 * t450;
t11 = qJD(1) * t510 - t39 * t453 + t40 * t450;
t10 = (qJD(2) * t512 - t71) * t452 + (qJD(1) * t513 + qJD(2) * t186 + t45 * t450 + t453 * t46) * t449;
t9 = (qJD(2) * t514 - t70) * t452 + (qJD(1) * t515 + qJD(2) * t185 + t43 * t450 + t44 * t453) * t449;
t8 = (qJD(2) * t508 - t69) * t452 + (qJD(1) * t509 + qJD(2) * t188 + t41 * t450 + t42 * t453) * t449;
t7 = (qJD(2) * t510 - t68) * t452 + (qJD(1) * t511 + qJD(2) * t187 + t39 * t450 + t40 * t453) * t449;
t6 = [t461 + (t250 * t332 + t251 * t331) * t720 + (t178 * t229 + t179 * t228) * t719 + (t116 * t207 + t117 * t206) * t718 + (t180 * t94 + t181 * t93) * t717 + (t113 * t55 + t114 * t54) * t716 - t272 * t668 - t345 * t579 - t347 * t572 + (-t340 * t445 + t341 * t446) * t449 + (-Icges(4,3) * t452 + t449 * t518 - t523 + t529) * t621 + (-Icges(4,3) * t449 - t363 * t445 + t364 * t446 - t452 * t518 + t524 + t528) * t619 + t729; m(4) * (t178 * t312 + t179 * t313 + t208 * t228 + t209 * t229) + m(5) * (t116 * t226 + t117 * t227 + t129 * t206 + t130 * t207) + m(6) * (t107 * t180 + t108 * t181 + t204 * t93 + t205 * t94) + m(7) * (t113 * t64 + t114 * t65 + t160 * t54 + t161 * t55) + (-t324 * t363 / 0.2e1 - t325 * t364 / 0.2e1 + t487 * t340 / 0.2e1 - t386 * t341 / 0.2e1 + m(3) * (-t251 * t411 - t331 * t397) + t394 * t700 + (t373 * t694 + t478 * t703 + Icges(4,5) * t325 / 0.2e1 + Icges(4,6) * t324 / 0.2e1 + t470 * t738) * t452 - t463) * t453 + (t322 * t709 + t323 * t708 + t340 * t707 + t341 * t706 + m(3) * (-t250 * t411 - t332 * t397) + t394 * t703 + (t372 * t694 + t478 * t701 - Icges(4,5) * t323 / 0.2e1 - Icges(4,6) * t322 / 0.2e1 + t469 * t738) * t452 + t462) * t450 + ((-qJD(1) * t374 - t214 * t445 + t216 * t446 - t453 * t479) * t703 + (qJD(1) * t375 - t215 * t445 + t217 * t446 - t450 * t479) * t701) * t449 + (-t482 / 0.2e1 + (t297 * t449 - t299 * t665 + t301 * t663) * t703 + t483 / 0.2e1 + (t296 * t449 - t298 * t665 + t300 * t663) * t701) * qJD(2) + ((t363 * t707 + t364 * t706 - t332 * t697 + (t373 / 0.2e1 - t297 / 0.2e1) * t452 + (t375 / 0.2e1 + t299 * t705 + t301 * t704) * t449 + t457) * t453 + (-t487 * t709 + t386 * t708 + t331 * t697 + (-t296 / 0.2e1 + t372 / 0.2e1) * t452 + (t298 * t705 + t300 * t704 + t374 / 0.2e1) * t449 + t456) * t450) * qJD(1); (t109 * t38 + t129 * t227 + t130 * t226) * t718 + (t182 * t81 + t208 * t313 + t209 * t312) * t719 + ((t450 * t378 + t379 * t453) * ((qJD(1) * t378 - t453 * t484 + t629) * t453 + (-t450 * t484 + (-t379 + t726) * qJD(1)) * t450) + t626 * t411 * t397) * t720 + (t16 * t63 + t160 * t65 + t161 * t64) * t716 + (t107 * t205 + t108 * t204 + t26 * t97) * t717 - t453 * ((t307 * t453 + (t221 + t483) * qJD(1)) * t453 + (t220 * qJD(1) + (-t373 * t619 - t375 * t621 + t625) * t450 + (-t306 + (t372 * t452 + t374 * t449) * qJD(2) - t497 * qJD(1)) * t453) * t450) + t450 * ((t450 * t306 + (t222 + t482) * qJD(1)) * t450 + (t223 * qJD(1) + (t372 * t619 + t374 * t621) * t453 + (-t307 + (-t373 * t452 - t375 * t449) * qJD(2) + (t371 - t498) * qJD(1)) * t450) * t453) - t453 * ((t487 * t215 - t386 * t217 - t324 * t298 - t325 * t300 + (t138 - t603) * qJD(1)) * t453 + (-t487 * t214 + t386 * t216 + t324 * t299 + t325 * t301 + (t137 + t601) * qJD(1)) * t450) + t450 * ((t387 * t214 + t388 * t216 + t322 * t299 + t323 * t301 + (t139 - t602) * qJD(1)) * t450 + (-t387 * t215 - t388 * t217 - t322 * t298 - t323 * t300 + (t140 + t604) * qJD(1)) * t453) + t450 * t11 + t450 * t12 + t450 * t4 - t453 * t14 - t453 * t13 - t453 * t5 + (t36 - t515 - t513 + (-t137 - t220) * t453 + (t138 + t221) * t450) * t624 + (t37 - t511 - t509 + (-t139 - t222) * t453 + (t140 + t223) * t450) * t622; 0.2e1 * ((t113 * t453 + t114 * t450) * t712 + (t206 * t453 + t207 * t450) * t714 + (t180 * t453 + t181 * t450) * t713 + (t228 * t453 + t229 * t450) * t715) * t619 + 0.2e1 * ((-t113 * t624 + t114 * t622 + t450 * t54 + t453 * t55) * t712 + (t116 * t450 + t117 * t453 - t206 * t624 + t207 * t622) * t714 + (-t180 * t624 + t181 * t622 + t450 * t93 + t453 * t94) * t713 + (t178 * t450 + t179 * t453 - t228 * t624 + t229 * t622) * t715) * t449; 0.2e1 * ((t160 * t620 + t161 * t618 - t16) * t712 + (t204 * t620 + t205 * t618 - t26) * t713 + (t226 * t620 + t227 * t618 - t38) * t714 + (t312 * t620 + t313 * t618 - t81) * t715) * t452 + 0.2e1 * ((qJD(2) * t63 + t160 * t622 - t161 * t624 + t450 * t65 + t453 * t64) * t712 + (qJD(2) * t97 + t107 * t453 + t108 * t450 + t204 * t622 - t205 * t624) * t713 + (qJD(2) * t109 + t129 * t453 + t130 * t450 + t226 * t622 - t227 * t624) * t714 + (qJD(2) * t182 + t208 * t453 + t209 * t450 + t312 * t622 - t313 * t624) * t715) * t449; 0.4e1 * (t715 + t714 + t611) * (-0.1e1 + t626) * t576; m(5) * (t116 * t211 + t117 * t210 + t206 * t99 + t207 * t98) + m(6) * (t152 * t94 + t153 * t93 + t180 * t57 + t181 * t56) + m(7) * (t100 * t55 + t101 * t54 + t113 * t24 + t114 * t23) + ((t450 * t456 + t453 * t457) * qJD(2) - t690) * t452 + (t613 + (qJD(1) * t465 + t463) * t450 + (qJD(1) * t456 + t462) * t453) * t449 - t612 - t649; -t695 / 0.2e1 + m(5) * (t109 * t62 + t129 * t210 + t130 * t211 + t189 * t38 + t226 * t98 + t227 * t99) + m(6) * (t107 * t152 + t108 * t153 + t112 * t26 + t204 * t56 + t205 * t57 + t33 * t97) + m(7) * (t100 * t64 + t101 * t65 + t15 * t63 + t16 * t80 + t160 * t23 + t161 * t24) + (-t9 / 0.2e1 - t10 / 0.2e1 + t554 * t619 - t552) * t453 + (t8 / 0.2e1 + t7 / 0.2e1 + t555 * t619 - t553) * t450 + ((-t450 * t554 + t453 * t555) * qJD(1) - t460 + (t13 + t14) * t703 + (t11 + t12) * t700 + (t450 * t655 - t453 * t727) * qJD(2) / 0.2e1) * t449 - (t721 * qJD(1) + t732 * t450 + t453 * t733) * t452 / 0.2e1 + (t450 * t731 + t453 * t730) * t693; 0.2e1 * ((t210 * t618 + t211 * t620 - t62) * t714 + (t152 * t618 + t153 * t620 - t33) * t713 + (t100 * t618 + t101 * t620 - t15) * t712) * t452 + 0.2e1 * ((qJD(2) * t189 - t210 * t624 + t211 * t622 + t450 * t98 + t453 * t99) * t714 + (qJD(2) * t112 - t152 * t624 + t153 * t622 + t450 * t56 + t453 * t57) * t713 + (qJD(2) * t80 - t100 * t624 + t101 * t622 + t23 * t450 + t24 * t453) * t712) * t449; (t189 * t62 + t210 * t99 + t211 * t98) * t718 + (t112 * t33 + t152 * t57 + t153 * t56) * t717 + (t100 * t24 + t101 * t23 + t15 * t80) * t716 + (-t3 + t690 * t452 + (t480 * t450 - t453 * t723) * qJD(2) + t649) * t452 + ((-t732 * t452 + t1 + t7 + t8) * t453 + (t733 * t452 + t10 + t2 + t9) * t450 + (t721 * t449 + t452 * t737 + t34) * qJD(2) + (t450 * t723 + t480 * t453) * qJD(1)) * t449; m(7) * (-t113 * t241 + t114 * t243 + t366 * t54 + t368 * t55) + m(6) * (-t180 * t241 + t181 * t243 + t366 * t93 + t368 * t94); m(7) * (t63 * t571 + t160 * t243 - t161 * t241 + t366 * t65 + t368 * t64 + (t16 * t449 + t619 * t63) * t435) + m(6) * (t97 * t571 + t107 * t368 + t108 * t366 + t204 * t243 - t205 * t241 + (t26 * t449 + t619 * t97) * t435); 0.2e1 * t611 * ((t435 * t442 + (t366 * t450 + t368 * t453 - t435 * t452) * t452) * qJD(2) + (-t570 - t241 * t453 + t243 * t450 + (t366 * t453 - t368 * t450) * qJD(1)) * t449); m(7) * (t80 * t571 - t100 * t241 + t101 * t243 + t23 * t366 + t24 * t368 + (t15 * t449 + t619 * t80) * t435) + m(6) * (t112 * t468 - t152 * t241 + t153 * t243 + t33 * t668 + t366 * t56 + t368 * t57); 0.4e1 * t611 * (-t241 * t368 + t243 * t366 + (qJD(4) * t436 * t442 + t435 * t576) * t435); m(7) * (t113 * t48 + t114 * t47 + t145 * t55 + t146 * t54) + ((-t84 / 0.2e1 + t710) * t453 + t599 * t450) * t619 + (-t613 + (t106 * t693 + t607) * t450 + (qJD(1) * t599 + t606) * t453) * t449 + t612; m(7) * (t110 * t16 + t145 * t64 + t146 * t65 + t160 * t47 + t161 * t48 + t25 * t63) + t695 / 0.2e1 + (t37 * t568 + t552) * t453 + (t36 * t568 + t553) * t450 + ((t36 * t701 + t37 * t703) * qJD(1) + t460) * t449; m(7) * ((-t25 + (t145 * t453 + t146 * t450) * qJD(2)) * t452 + (qJD(2) * t110 + t450 * t47 + t453 * t48 + (-t145 * t450 + t146 * t453) * qJD(1)) * t449); m(7) * (t100 * t48 + t101 * t47 + t110 * t15 + t145 * t24 + t146 * t23 + t25 * t80) + t734; m(7) * (t110 * t468 - t145 * t241 + t146 * t243 + t25 * t668 + t366 * t47 + t368 * t48); (t110 * t25 + t145 * t48 + t146 * t47) * t716 - t734;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
