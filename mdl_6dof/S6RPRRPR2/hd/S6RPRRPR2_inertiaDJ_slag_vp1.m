% Calculate time derivative of joint inertia matrix for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:16
% EndTime: 2019-03-09 05:00:50
% DurationCPUTime: 21.33s
% Computational Cost: add. (62883->1101), mult. (56668->1482), div. (0->0), fcn. (53827->12), ass. (0->545)
t425 = sin(qJ(3));
t428 = cos(qJ(3));
t421 = qJ(4) + pkin(11);
t412 = sin(t421);
t414 = cos(t421);
t499 = Icges(6,5) * t414 - Icges(6,6) * t412;
t339 = -Icges(6,3) * t428 + t425 * t499;
t424 = sin(qJ(4));
t427 = cos(qJ(4));
t500 = Icges(5,5) * t427 - Icges(5,6) * t424;
t352 = -Icges(5,3) * t428 + t425 * t500;
t719 = t352 + t339;
t658 = Icges(6,4) * t412;
t508 = Icges(6,1) * t414 - t658;
t341 = -Icges(6,5) * t428 + t425 * t508;
t660 = Icges(5,4) * t424;
t509 = Icges(5,1) * t427 - t660;
t354 = -Icges(5,5) * t428 + t425 * t509;
t718 = t341 * t414 + t354 * t427;
t422 = qJ(1) + pkin(10);
t413 = sin(t422);
t415 = cos(t422);
t624 = t424 * t428;
t347 = -t413 * t624 - t415 * t427;
t623 = t427 * t428;
t631 = t415 * t424;
t348 = t413 * t623 - t631;
t634 = t413 * t425;
t261 = Icges(5,5) * t348 + Icges(5,6) * t347 + Icges(5,3) * t634;
t263 = Icges(5,4) * t348 + Icges(5,2) * t347 + Icges(5,6) * t634;
t265 = Icges(5,1) * t348 + Icges(5,4) * t347 + Icges(5,5) * t634;
t349 = t413 * t427 - t415 * t624;
t635 = t413 * t424;
t350 = t415 * t623 + t635;
t630 = t415 * t425;
t120 = t261 * t630 + t263 * t349 + t265 * t350;
t262 = Icges(5,5) * t350 + Icges(5,6) * t349 + Icges(5,3) * t630;
t264 = Icges(5,4) * t350 + Icges(5,2) * t349 + Icges(5,6) * t630;
t266 = Icges(5,1) * t350 + Icges(5,4) * t349 + Icges(5,5) * t630;
t121 = t262 * t630 + t264 * t349 + t266 * t350;
t484 = t120 * t413 + t121 * t415;
t633 = t413 * t428;
t334 = -t412 * t633 - t414 * t415;
t335 = -t412 * t415 + t414 * t633;
t235 = Icges(6,5) * t335 + Icges(6,6) * t334 + Icges(6,3) * t634;
t237 = Icges(6,4) * t335 + Icges(6,2) * t334 + Icges(6,6) * t634;
t239 = Icges(6,1) * t335 + Icges(6,4) * t334 + Icges(6,5) * t634;
t629 = t415 * t428;
t336 = -t412 * t629 + t413 * t414;
t337 = t412 * t413 + t414 * t629;
t108 = t235 * t630 + t237 * t336 + t239 * t337;
t236 = Icges(6,5) * t337 + Icges(6,6) * t336 + Icges(6,3) * t630;
t238 = Icges(6,4) * t337 + Icges(6,2) * t336 + Icges(6,6) * t630;
t240 = Icges(6,1) * t337 + Icges(6,4) * t336 + Icges(6,5) * t630;
t109 = t236 * t630 + t238 * t336 + t240 * t337;
t490 = t108 * t413 + t109 * t415;
t717 = t484 + t490;
t118 = t261 * t634 + t263 * t347 + t265 * t348;
t119 = t262 * t634 + t264 * t347 + t266 * t348;
t486 = t118 * t413 + t119 * t415;
t106 = t235 * t634 + t237 * t334 + t239 * t335;
t107 = t236 * t634 + t238 * t334 + t240 * t335;
t492 = t106 * t413 + t107 * t415;
t716 = t486 + t492;
t591 = qJD(1) * t425;
t715 = -t591 / 0.2e1;
t657 = Icges(6,4) * t414;
t503 = -Icges(6,2) * t412 + t657;
t340 = -Icges(6,6) * t428 + t425 * t503;
t659 = Icges(5,4) * t427;
t504 = -Icges(5,2) * t424 + t659;
t353 = -Icges(5,6) * t428 + t425 * t504;
t714 = t719 * t428 + (t340 * t412 + t353 * t424 - t718) * t425;
t586 = qJD(3) * t428;
t546 = t586 / 0.2e1;
t562 = t413 * t591;
t713 = t415 * t546 - t562 / 0.2e1;
t592 = qJD(1) * t415;
t547 = t592 / 0.2e1;
t712 = t413 * t546 + t425 * t547;
t417 = t427 * pkin(4);
t409 = t417 + pkin(3);
t671 = pkin(3) - t409;
t548 = t671 * t428;
t674 = pkin(8) * t425;
t711 = t548 + t674;
t554 = t413 * t586;
t560 = t415 * t591;
t443 = t554 + t560;
t423 = -qJ(5) - pkin(8);
t584 = qJD(4) * t427;
t710 = pkin(4) * t584 + t423 * t591;
t587 = qJD(3) * t425;
t553 = t415 * t587;
t590 = qJD(1) * t428;
t537 = -qJD(4) + t590;
t696 = t537 * t413;
t709 = t696 + t553;
t480 = -t238 * t412 + t240 * t414;
t123 = -t236 * t428 + t425 * t480;
t478 = -t264 * t424 + t266 * t427;
t142 = -t262 * t428 + t425 * t478;
t620 = t123 + t142;
t481 = -t237 * t412 + t239 * t414;
t122 = -t235 * t428 + t425 * t481;
t479 = -t263 * t424 + t265 * t427;
t141 = -t261 * t428 + t425 * t479;
t621 = t122 + t141;
t708 = t413 * t621 + t415 * t620;
t538 = -qJD(4) * t428 + qJD(1);
t469 = t414 * t538;
t555 = t413 * t587;
t699 = t415 * t537;
t224 = t413 * t469 + (-t699 + t555) * t412;
t467 = t538 * t412;
t225 = t414 * t699 + (-t414 * t587 + t467) * t413;
t145 = Icges(6,5) * t225 + Icges(6,6) * t224 + Icges(6,3) * t443;
t147 = Icges(6,4) * t225 + Icges(6,2) * t224 + Icges(6,6) * t443;
t149 = Icges(6,1) * t225 + Icges(6,4) * t224 + Icges(6,5) * t443;
t222 = t412 * t709 + t415 * t469;
t223 = -t414 * t709 + t415 * t467;
t552 = t415 * t586;
t442 = t552 - t562;
t35 = t145 * t630 + t147 * t336 + t149 * t337 + t222 * t237 + t223 * t239 + t235 * t442;
t144 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t442;
t146 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t442;
t148 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t442;
t36 = t144 * t630 + t146 * t336 + t148 * t337 + t222 * t238 + t223 * t240 + t236 * t442;
t440 = t424 * t587 + t427 * t538;
t258 = t413 * t440 - t424 * t699;
t439 = t424 * t538 - t427 * t587;
t259 = t413 * t439 + t427 * t699;
t166 = Icges(5,5) * t259 + Icges(5,6) * t258 + Icges(5,3) * t443;
t168 = Icges(5,4) * t259 + Icges(5,2) * t258 + Icges(5,6) * t443;
t170 = Icges(5,1) * t259 + Icges(5,4) * t258 + Icges(5,5) * t443;
t256 = t415 * t440 + t424 * t696;
t257 = t415 * t439 - t427 * t696;
t41 = t166 * t630 + t168 * t349 + t170 * t350 + t256 * t263 + t257 * t265 + t261 * t442;
t165 = Icges(5,5) * t257 + Icges(5,6) * t256 + Icges(5,3) * t442;
t167 = Icges(5,4) * t257 + Icges(5,2) * t256 + Icges(5,6) * t442;
t169 = Icges(5,1) * t257 + Icges(5,4) * t256 + Icges(5,5) * t442;
t42 = t165 * t630 + t167 * t349 + t169 * t350 + t256 * t264 + t257 * t266 + t262 * t442;
t707 = (-t35 - t41) * t415 + (t36 + t42) * t413 + t717 * qJD(1);
t37 = t145 * t634 + t147 * t334 + t149 * t335 + t224 * t237 + t225 * t239 + t235 * t443;
t38 = t144 * t634 + t146 * t334 + t148 * t335 + t224 * t238 + t225 * t240 + t236 * t443;
t43 = t166 * t634 + t168 * t347 + t170 * t348 + t258 * t263 + t259 * t265 + t261 * t443;
t44 = t165 * t634 + t167 * t347 + t169 * t348 + t258 * t264 + t259 * t266 + t262 * t443;
t706 = (-t37 - t43) * t415 + (t38 + t44) * t413 + t716 * qJD(1);
t39 = (qJD(3) * t481 - t145) * t428 + (qJD(3) * t235 - t147 * t412 + t149 * t414 + (-t237 * t414 - t239 * t412) * qJD(4)) * t425;
t45 = (qJD(3) * t479 - t166) * t428 + (qJD(3) * t261 - t168 * t424 + t170 * t427 + (-t263 * t427 - t265 * t424) * qJD(4)) * t425;
t705 = -t39 - t45;
t40 = (qJD(3) * t480 - t144) * t428 + (qJD(3) * t236 - t146 * t412 + t148 * t414 + (-t238 * t414 - t240 * t412) * qJD(4)) * t425;
t46 = (qJD(3) * t478 - t165) * t428 + (qJD(3) * t262 - t167 * t424 + t169 * t427 + (-t264 * t427 - t266 * t424) * qJD(4)) * t425;
t704 = t40 + t46;
t175 = t334 * t340 + t335 * t341 + t339 * t634;
t182 = t347 * t353 + t348 * t354 + t352 * t634;
t703 = (-t175 - t182) * t428 + t716 * t425;
t176 = t336 * t340 + t337 * t341 + t339 * t630;
t183 = t349 * t353 + t350 * t354 + t352 * t630;
t669 = (-t176 - t183) * t428 + t717 * t425;
t668 = rSges(4,1) * t428;
t522 = -rSges(4,2) * t425 + t668;
t667 = rSges(4,3) * t415;
t328 = t413 * t522 - t667;
t695 = -rSges(4,2) * t630 + rSges(4,3) * t413;
t329 = rSges(4,1) * t629 + t695;
t391 = rSges(4,1) * t425 + rSges(4,2) * t428;
t453 = qJD(3) * t391;
t436 = rSges(4,2) * t562 + rSges(4,3) * t592 - t415 * t453;
t140 = (qJD(1) * t328 + t436) * t415 + (-t413 * t453 + (-t329 + t695) * qJD(1)) * t413;
t689 = 2 * m(4);
t702 = t140 * t689;
t662 = Icges(4,4) * t425;
t511 = Icges(4,1) * t428 - t662;
t324 = Icges(4,5) * t413 + t415 * t511;
t641 = t324 * t428;
t661 = Icges(4,4) * t428;
t506 = -Icges(4,2) * t425 + t661;
t322 = Icges(4,6) * t413 + t415 * t506;
t646 = t322 * t425;
t472 = -t641 + t646;
t700 = t415 * t472;
t697 = t425 * t671;
t418 = cos(qJ(1)) * pkin(1);
t694 = pkin(7) * t413 + t418;
t585 = qJD(4) * t425;
t278 = (-Icges(6,5) * t412 - Icges(6,6) * t414) * t585 + (Icges(6,3) * t425 + t428 * t499) * qJD(3);
t280 = (-Icges(6,1) * t412 - t657) * t585 + (Icges(6,5) * t425 + t428 * t508) * qJD(3);
t288 = (-Icges(5,5) * t424 - Icges(5,6) * t427) * t585 + (Icges(5,3) * t425 + t428 * t500) * qJD(3);
t290 = (-Icges(5,1) * t424 - t659) * t585 + (Icges(5,5) * t425 + t428 * t509) * qJD(3);
t693 = (t290 * t427 - t353 * t584) * t425 + (t280 * t425 - t340 * t585) * t414 + t719 * t587 + t718 * t586 + (-t288 - t278) * t428;
t501 = Icges(4,5) * t428 - Icges(4,6) * t425;
t319 = -Icges(4,3) * t415 + t413 * t501;
t321 = -Icges(4,6) * t415 + t413 * t506;
t323 = -Icges(4,5) * t415 + t413 * t511;
t420 = qJD(4) + qJD(6);
t540 = -t420 + t590;
t692 = t413 * t540 + t553;
t691 = t415 * t540 - t555;
t690 = -t428 * t620 + t669;
t419 = -pkin(9) + t423;
t595 = t419 - t423;
t675 = pkin(5) * t414;
t373 = t409 + t675;
t598 = t373 - t409;
t294 = t425 * t598 + t428 * t595;
t688 = 2 * m(5);
t687 = 2 * m(6);
t686 = 2 * m(7);
t685 = m(6) / 0.2e1;
t684 = m(7) / 0.2e1;
t683 = t413 / 0.2e1;
t682 = t415 / 0.2e1;
t681 = -t428 / 0.2e1;
t680 = -rSges(5,3) - pkin(8);
t679 = m(4) * t391;
t678 = sin(qJ(1)) * pkin(1);
t677 = pkin(3) * t428;
t676 = pkin(4) * t424;
t670 = pkin(8) + t423;
t666 = rSges(6,3) * t425;
t665 = rSges(7,3) * t425;
t664 = -rSges(7,3) + t419;
t279 = (-Icges(6,2) * t414 - t658) * t585 + (Icges(6,6) * t425 + t428 * t503) * qJD(3);
t289 = (-Icges(5,2) * t427 - t660) * t585 + (Icges(5,6) * t425 + t428 * t504) * qJD(3);
t556 = t353 * t586;
t557 = t340 * t586;
t663 = (-t556 + (-qJD(4) * t354 - t289) * t425) * t424 + (-t557 + (-qJD(4) * t341 - t279) * t425) * t412 + t693;
t416 = qJ(6) + t421;
t407 = sin(t416);
t656 = Icges(7,4) * t407;
t408 = cos(t416);
t655 = Icges(7,4) * t408;
t517 = -rSges(6,1) * t335 - rSges(6,2) * t334;
t242 = rSges(6,3) * t634 - t517;
t651 = t242 * t415;
t520 = -rSges(5,1) * t348 - rSges(5,2) * t347;
t269 = rSges(5,3) * t634 - t520;
t650 = t269 * t415;
t507 = Icges(7,1) * t408 - t656;
t318 = -Icges(7,5) * t428 + t425 * t507;
t649 = t318 * t408;
t648 = t321 * t425;
t647 = t321 * t428;
t645 = t322 * t428;
t644 = t323 * t425;
t643 = t323 * t428;
t642 = t324 * t425;
t376 = pkin(5) * t412 + t676;
t362 = t376 * qJD(4);
t638 = t362 * t428;
t637 = t373 * t425;
t636 = t376 * t415;
t363 = (t417 + t675) * qJD(4);
t632 = t415 * t363;
t628 = t419 * t425;
t627 = t420 * t425;
t626 = t423 * t425;
t625 = t423 * t428;
t541 = -t420 * t428 + qJD(1);
t470 = t415 * t541;
t205 = t407 * t692 + t408 * t470;
t206 = t407 * t470 - t408 * t692;
t574 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t552;
t138 = -rSges(7,3) * t562 + t574;
t531 = -t362 * t629 + t363 * t413 + t376 * t592 + t419 * t562;
t581 = qJD(4) * t676;
t539 = t428 * t581;
t561 = t413 * t590;
t398 = pkin(4) * t631;
t566 = qJD(1) * t398 + t413 * t710;
t622 = -t598 * t561 + (-qJD(3) * t294 + t539) * t415 + t531 - t566 + t138;
t471 = t413 * t541;
t207 = -t407 * t691 + t408 * t471;
t208 = t407 * t471 + t408 * t691;
t515 = t208 * rSges(7,1) + t207 * rSges(7,2);
t139 = rSges(7,3) * t443 + t515;
t312 = -t407 * t633 - t408 * t415;
t313 = -t407 * t415 + t408 * t633;
t514 = -rSges(7,1) * t313 - rSges(7,2) * t312;
t232 = rSges(7,3) * t634 - t514;
t619 = t139 * t630 + t232 * t552;
t572 = rSges(6,1) * t223 + rSges(6,2) * t222 + rSges(6,3) * t552;
t150 = -rSges(6,3) * t562 + t572;
t385 = pkin(8) * t552;
t583 = qJD(5) * t425;
t390 = t415 * t583;
t530 = t390 + t566;
t593 = qJD(1) * t413;
t163 = -t385 + t711 * t593 + (-t539 + (-t625 + t697) * qJD(3)) * t415 + t530;
t618 = -t150 - t163;
t383 = pkin(3) * t555;
t397 = pkin(4) * t635;
t512 = t409 * t555 + t413 * t539 + t415 * t710 + t423 * t554;
t164 = t383 + (-pkin(8) * t586 + t583) * t413 + (-t415 * t711 + t397) * qJD(1) - t512;
t597 = t413 * t626 + t398;
t267 = -t413 * t711 - t597;
t617 = t164 * t630 + t267 * t552;
t616 = t714 * t587;
t543 = t598 * t428;
t451 = t543 - t628;
t200 = t413 * t451 + t597 - t636;
t615 = t200 + t232;
t542 = t595 * t425;
t599 = -t409 * t629 - t397;
t601 = t373 * t629 + t376 * t413;
t201 = -t415 * t542 + t599 + t601;
t314 = -t407 * t629 + t408 * t413;
t315 = t407 * t413 + t408 * t629;
t233 = rSges(7,1) * t315 + rSges(7,2) * t314 + rSges(7,3) * t630;
t614 = t201 + t233;
t513 = rSges(7,1) * t408 - rSges(7,2) * t407;
t325 = -rSges(7,3) * t428 + t425 * t513;
t613 = t233 * t587 + t325 * t562;
t612 = -t242 - t267;
t243 = rSges(6,1) * t337 + rSges(6,2) * t336 + rSges(6,3) * t630;
t400 = pkin(3) * t629;
t357 = pkin(8) * t630 + t400;
t462 = -t415 * t626 - t599;
t268 = t462 - t357;
t611 = -t243 - t268;
t333 = t428 * t670 - t697;
t610 = t268 * t587 + t333 * t562;
t609 = t267 * t428 + t333 * t634;
t270 = rSges(5,1) * t350 + rSges(5,2) * t349 + rSges(5,3) * t630;
t608 = -t270 - t357;
t516 = rSges(6,1) * t414 - rSges(6,2) * t412;
t281 = (-rSges(6,1) * t412 - rSges(6,2) * t414) * t585 + (t428 * t516 + t666) * qJD(3);
t551 = t424 * t585;
t287 = -pkin(4) * t551 - qJD(5) * t428 + (-t425 * t670 - t548) * qJD(3);
t607 = -t281 - t287;
t519 = rSges(5,1) * t427 - rSges(5,2) * t424;
t291 = (-rSges(5,1) * t424 - rSges(5,2) * t427) * t585 + (rSges(5,3) * t425 + t428 * t519) * qJD(3);
t525 = t674 + t677;
t375 = t525 * qJD(3);
t606 = -t291 - t375;
t605 = t294 + t325;
t187 = t232 * t428 + t325 * t634;
t396 = pkin(3) * t425 - pkin(8) * t428;
t360 = t396 * t593;
t604 = t333 * t593 + t360;
t356 = t525 * t413;
t603 = t356 * t413 + t357 * t415;
t342 = -rSges(6,3) * t428 + t425 * t516;
t602 = -t333 - t342;
t355 = -rSges(5,3) * t428 + t425 * t519;
t600 = -t355 - t396;
t596 = t413 ^ 2 + t415 ^ 2;
t320 = Icges(4,3) * t413 + t415 * t501;
t594 = qJD(1) * t320;
t589 = qJD(3) * t413;
t588 = qJD(3) * t415;
t487 = t118 * t415 - t119 * t413;
t493 = t106 * t415 - t107 * t413;
t579 = -t493 / 0.2e1 - t487 / 0.2e1;
t485 = t120 * t415 - t121 * t413;
t491 = t108 * t415 - t109 * t413;
t578 = -t491 / 0.2e1 - t485 / 0.2e1;
t577 = -t163 - t622;
t576 = -t267 - t615;
t575 = -t268 - t614;
t218 = (-t362 + t581) * t425 + (t543 - t542) * qJD(3);
t253 = (-rSges(7,1) * t407 - rSges(7,2) * t408) * t627 + (t428 * t513 + t665) * qJD(3);
t573 = -t218 - t253 - t287;
t571 = rSges(5,1) * t257 + rSges(5,2) * t256 + rSges(5,3) * t552;
t570 = -t375 + t607;
t569 = t413 * (pkin(8) * t443 + qJD(1) * t400 - t383) + t415 * (-pkin(8) * t562 + t385 + (-t553 - t561) * pkin(3)) + t356 * t592;
t568 = -t333 - t605;
t567 = -t396 + t602;
t565 = pkin(2) * t415 + t694;
t564 = t342 * t593;
t563 = t355 * t593;
t502 = -Icges(7,2) * t407 + t655;
t317 = -Icges(7,6) * t428 + t425 * t502;
t558 = t317 * t586;
t550 = t634 / 0.2e1;
t549 = t630 / 0.2e1;
t405 = t415 * pkin(7);
t545 = t405 - t678;
t544 = -t373 * t428 - pkin(2);
t293 = t600 * t415;
t536 = t428 * t139 + t253 * t634 + t325 * t443;
t535 = t428 * t164 + t287 * t634 + t333 * t443;
t534 = -t375 + t573;
t533 = t413 * t267 + t268 * t415 + t603;
t532 = -t396 + t568;
t213 = t567 * t415;
t498 = Icges(7,5) * t408 - Icges(7,6) * t407;
t316 = -Icges(7,3) * t428 + t425 * t498;
t189 = -t316 * t428 + (-t317 * t407 + t649) * t425;
t226 = Icges(7,5) * t313 + Icges(7,6) * t312 + Icges(7,3) * t634;
t228 = Icges(7,4) * t313 + Icges(7,2) * t312 + Icges(7,6) * t634;
t230 = Icges(7,1) * t313 + Icges(7,4) * t312 + Icges(7,5) * t634;
t483 = -t228 * t407 + t230 * t408;
t114 = -t226 * t428 + t425 * t483;
t227 = Icges(7,5) * t315 + Icges(7,6) * t314 + Icges(7,3) * t630;
t229 = Icges(7,4) * t315 + Icges(7,2) * t314 + Icges(7,6) * t630;
t231 = Icges(7,1) * t315 + Icges(7,4) * t314 + Icges(7,5) * t630;
t482 = -t229 * t407 + t231 * t408;
t115 = -t227 * t428 + t425 * t482;
t488 = t114 * t413 + t115 * t415;
t162 = t314 * t317 + t315 * t318 + t316 * t630;
t133 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t443;
t135 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t443;
t137 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t443;
t26 = t133 * t630 + t135 * t314 + t137 * t315 + t205 * t228 + t206 * t230 + t226 * t442;
t132 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t442;
t134 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t442;
t136 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t442;
t27 = t132 * t630 + t134 * t314 + t136 * t315 + t205 * t229 + t206 * t231 + t227 * t442;
t103 = t226 * t630 + t228 * t314 + t230 * t315;
t104 = t227 * t630 + t229 * t314 + t231 * t315;
t494 = t103 * t413 + t104 * t415;
t495 = t103 * t415 - t104 * t413;
t246 = (-Icges(7,5) * t407 - Icges(7,6) * t408) * t627 + (Icges(7,3) * t425 + t428 * t498) * qJD(3);
t247 = (-Icges(7,2) * t408 - t656) * t627 + (Icges(7,6) * t425 + t428 * t502) * qJD(3);
t248 = (-Icges(7,1) * t407 - t655) * t627 + (Icges(7,5) * t425 + t428 * t507) * qJD(3);
t72 = t205 * t317 + t206 * t318 + t246 * t630 + t247 * t314 + t248 * t315 + t316 * t442;
t5 = (qJD(3) * t494 - t72) * t428 + (qJD(1) * t495 + qJD(3) * t162 + t26 * t413 + t27 * t415) * t425;
t161 = t312 * t317 + t313 * t318 + t316 * t634;
t101 = t226 * t634 + t228 * t312 + t230 * t313;
t102 = t227 * t634 + t229 * t312 + t231 * t313;
t496 = t101 * t413 + t102 * t415;
t52 = -t161 * t428 + t425 * t496;
t53 = -t162 * t428 + t425 * t494;
t28 = t133 * t634 + t135 * t312 + t137 * t313 + t207 * t228 + t208 * t230 + t226 * t443;
t29 = t132 * t634 + t134 * t312 + t136 * t313 + t207 * t229 + t208 * t231 + t227 * t443;
t497 = t101 * t415 - t102 * t413;
t73 = t207 * t317 + t208 * t318 + t246 * t634 + t247 * t312 + t248 * t313 + t316 * t443;
t6 = (qJD(3) * t496 - t73) * t428 + (qJD(1) * t497 + qJD(3) * t161 + t28 * t413 + t29 * t415) * t425;
t524 = t5 * t630 + t53 * t552 + t6 * t634 + (-t189 * t428 + t425 * t488) * t587 + t443 * t52;
t521 = rSges(5,1) * t259 + rSges(5,2) * t258;
t518 = rSges(6,1) * t225 + rSges(6,2) * t224;
t505 = Icges(4,2) * t428 + t662;
t489 = t114 * t415 - t115 * t413;
t477 = -t270 * t413 + t650;
t476 = -t269 * t413 - t270 * t415;
t473 = -t643 + t648;
t179 = t532 * t415;
t463 = -pkin(2) - t522;
t79 = t224 * t340 + t225 * t341 + t278 * t634 + t279 * t334 + t280 * t335 + t339 * t443;
t87 = t258 * t353 + t259 * t354 + t288 * t634 + t289 * t347 + t290 * t348 + t352 * t443;
t461 = t45 / 0.2e1 + t39 / 0.2e1 + t79 / 0.2e1 + t87 / 0.2e1;
t78 = t222 * t340 + t223 * t341 + t278 * t630 + t279 * t336 + t280 * t337 + t339 * t442;
t86 = t256 * t353 + t257 * t354 + t288 * t630 + t289 * t349 + t290 * t350 + t352 * t442;
t460 = t46 / 0.2e1 + t40 / 0.2e1 + t78 / 0.2e1 + t86 / 0.2e1;
t459 = -t409 * t428 - pkin(2) - t666;
t458 = -t428 * t621 + t703;
t457 = t163 * t415 + t164 * t413 + t267 * t592 + t569;
t455 = t141 / 0.2e1 + t122 / 0.2e1 + t175 / 0.2e1 + t182 / 0.2e1;
t454 = -t183 / 0.2e1 - t142 / 0.2e1 - t123 / 0.2e1 - t176 / 0.2e1;
t452 = t413 * t611 + t651;
t450 = t425 * t680 - pkin(2) - t677;
t448 = qJD(3) * t505;
t447 = qJD(3) * (-Icges(4,5) * t425 - Icges(4,6) * t428);
t446 = (-t419 * t428 - t637) * qJD(3);
t444 = t425 * t664 + t544;
t441 = t200 * t415 + t413 * t575;
t14 = qJD(1) * t494 - t26 * t415 + t27 * t413;
t15 = qJD(1) * t496 - t28 * t415 + t29 * t413;
t33 = (qJD(3) * t483 - t133) * t428 + (qJD(3) * t226 + (-t228 * t420 + t137) * t408 + (-t230 * t420 - t135) * t407) * t425;
t34 = (qJD(3) * t482 - t132) * t428 + (qJD(3) * t227 + (-t229 * t420 + t136) * t408 + (-t231 * t420 - t134) * t407) * t425;
t438 = t5 * t683 + t15 * t550 + t14 * t549 + (qJD(1) * t488 - t33 * t415 + t34 * t413) * t681 - t415 * t6 / 0.2e1 + t52 * t593 / 0.2e1 + t53 * t547 - t489 * t587 / 0.2e1 - t713 * t495 - t712 * t497;
t437 = t413 * t459 - t678;
t435 = -t428 * t246 + t316 * t587 + t586 * t649 + (t248 * t425 - t317 * t627) * t408;
t432 = t413 * t450 - t678;
t184 = t189 * t587;
t85 = (-t558 + (-t318 * t420 - t247) * t425) * t407 + t435;
t9 = t184 + (qJD(3) * t488 - t85) * t428 + (qJD(1) * t489 + t33 * t413 + t34 * t415) * t425;
t431 = -t428 * t9 - t53 * t562 + t524;
t430 = t184 + (t33 + t73) * t550 + (t34 + t72) * t549 + (t115 + t162) * t713 + (t114 + t161) * t712;
t402 = pkin(7) * t592;
t374 = t522 * qJD(3);
t367 = t501 * qJD(3);
t292 = t600 * t413;
t286 = t329 + t565;
t285 = t413 * t463 + t545 + t667;
t272 = t413 * t447 + t594;
t271 = -qJD(1) * t319 + t415 * t447;
t244 = t267 * t630;
t215 = t391 * t589 + (-t418 + (-rSges(4,3) - pkin(7)) * t413 + t463 * t415) * qJD(1);
t214 = t402 + (-t678 + (-pkin(2) - t668) * t413) * qJD(1) + t436;
t212 = t567 * t413;
t210 = t232 * t630;
t199 = -t270 * t428 - t355 * t630;
t198 = t269 * t428 + t355 * t634;
t197 = t565 - t608;
t196 = t405 + t432 + t520;
t193 = t320 * t413 - t700;
t192 = t319 * t413 - t415 * t473;
t191 = -t320 * t415 - t413 * t472;
t190 = -t319 * t415 - t413 * t473;
t188 = -t233 * t428 - t325 * t630;
t186 = qJD(1) * t293 + t413 * t606;
t185 = t415 * t606 + t360 + t563;
t181 = t462 + t565 + t243;
t180 = t405 + t437 + t517 + t597;
t178 = t532 * t413;
t177 = t477 * t425;
t174 = -t415 * t628 + t233 + t565 + t601;
t173 = t413 * t444 + t514 + t545 + t636;
t172 = rSges(5,3) * t443 + t521;
t171 = -rSges(5,3) * t562 + t571;
t156 = -t233 * t634 + t210;
t151 = rSges(6,3) * t443 + t518;
t143 = -t476 + t603;
t130 = t428 * t611 + t602 * t630;
t129 = t242 * t428 + t342 * t634 + t609;
t127 = qJD(1) * t213 + t413 * t570;
t126 = t415 * t570 + t564 + t604;
t125 = t383 + t680 * t554 + (t415 * t450 - t694) * qJD(1) - t521;
t124 = -pkin(3) * t553 + qJD(1) * t432 + t385 + t402 + t571;
t117 = -t632 + (t446 - t638) * t413 + ((t376 - t676) * t413 + t451 * t415) * qJD(1) + t512;
t105 = t425 * t452 + t244;
t99 = (-rSges(6,3) * t586 - t583) * t413 + (-t418 + (-pkin(7) - t676) * t413 + t459 * t415) * qJD(1) + t512 - t518;
t98 = t402 + (-t539 + (-t409 * t425 - t625) * qJD(3)) * t415 + t437 * qJD(1) + t530 + t572;
t97 = t242 * t413 + t243 * t415 + t533;
t96 = (t355 * t589 + t172) * t428 + (-qJD(3) * t269 + t291 * t413 + t355 * t592) * t425;
t95 = (-t355 * t588 - t171) * t428 + (qJD(3) * t270 - t291 * t415 + t563) * t425;
t94 = t428 * t575 + t568 * t630;
t93 = t200 * t428 + t294 * t634 + t187 + t609;
t92 = qJD(1) * t179 + t413 * t534;
t91 = t415 * t534 + t593 * t605 + t604;
t89 = t632 + (-t583 + t638 + (t428 * t664 + t637) * qJD(3)) * t413 + (-t418 + (-pkin(7) - t376) * t413 + t444 * t415) * qJD(1) - t515;
t88 = t390 + t402 + t415 * t446 + (-t678 + (t544 - t665) * t413) * qJD(1) + t531 + t574;
t84 = -t232 * t587 + t536;
t83 = -t138 * t428 + (-t253 * t425 - t325 * t586) * t415 + t613;
t80 = t425 * t441 + t210 + t244;
t76 = t413 * t615 + t415 * t614 + t533;
t67 = t477 * t586 + (qJD(1) * t476 - t171 * t413 + t172 * t415) * t425;
t60 = t171 * t415 + t172 * t413 + (t413 * t608 + t650) * qJD(1) + t569;
t56 = (t342 * t589 + t151) * t428 + (qJD(3) * t612 + t281 * t413 + t342 * t592) * t425 + t535;
t55 = t618 * t428 + (qJD(3) * t243 + t564) * t425 + (t425 * t607 + t586 * t602) * t415 + t610;
t54 = -t233 * t560 + (-t233 * t586 + (-qJD(1) * t232 - t138) * t425) * t413 + t619;
t30 = t150 * t415 + t151 * t413 + (t651 + (-t357 + t611) * t413) * qJD(1) + t457;
t25 = t452 * t586 + (t151 * t415 + t618 * t413 + (t413 * t612 + t415 * t611) * qJD(1)) * t425 + t617;
t24 = (t294 * t589 + t117) * t428 + (qJD(3) * t576 + t218 * t413 + t294 * t592) * t425 + t535 + t536;
t23 = (qJD(3) * t201 + t294 * t593) * t425 + t577 * t428 + (t425 * t573 + t568 * t586) * t415 + t610 + t613;
t22 = t622 * t415 + (t117 + t139) * t413 + (t615 * t415 + (-t357 + t575) * t413) * qJD(1) + t457;
t21 = t441 * t586 + (t117 * t415 + t577 * t413 + (t413 * t576 + t415 * t575) * qJD(1)) * t425 + t617 + t619;
t11 = (qJD(3) * t486 - t87) * t428 + (qJD(1) * t487 + qJD(3) * t182 + t413 * t43 + t415 * t44) * t425;
t10 = (qJD(3) * t484 - t86) * t428 + (qJD(1) * t485 + qJD(3) * t183 + t41 * t413 + t415 * t42) * t425;
t8 = (qJD(3) * t492 - t79) * t428 + (qJD(1) * t493 + qJD(3) * t175 + t37 * t413 + t38 * t415) * t425;
t7 = (qJD(3) * t490 - t78) * t428 + (qJD(1) * t491 + qJD(3) * t176 + t35 * t413 + t36 * t415) * t425;
t1 = [t435 - t354 * t551 + (t173 * t89 + t174 * t88) * t686 + (t180 * t99 + t181 * t98) * t687 + (t124 * t197 + t125 * t196) * t688 + (t214 * t286 + t215 * t285) * t689 + t693 + (-t505 + t511) * t587 + (Icges(4,1) * t425 + t506 + t661) * t586 + (-t289 * t425 - t556) * t424 + (-t279 * t425 - t341 * t585 - t557) * t412 + (-t247 * t425 - t318 * t627 - t558) * t407; 0; 0; m(7) * (t173 * t91 + t174 * t92 + t178 * t88 + t179 * t89) + m(6) * (t126 * t180 + t127 * t181 + t212 * t98 + t213 * t99) + m(5) * (t124 * t292 + t125 * t293 + t185 * t196 + t186 * t197) + ((qJD(1) * t322 - t413 * t448) * t681 + t324 * t715 - t33 / 0.2e1 - t73 / 0.2e1 + m(4) * (-t215 * t391 - t285 * t374) + t367 * t682 + (t648 / 0.2e1 - t643 / 0.2e1) * qJD(3) - t461) * t415 + ((-qJD(1) * t321 - t415 * t448) * t428 / 0.2e1 + t323 * t715 + t72 / 0.2e1 + t34 / 0.2e1 + m(4) * (-t214 * t391 - t286 * t374) + t367 * t683 + (-t646 / 0.2e1 + t641 / 0.2e1) * qJD(3) + t460) * t413 + ((t645 / 0.2e1 + t642 / 0.2e1 - t286 * t679 + t162 / 0.2e1 + t115 / 0.2e1 - t454) * t415 + (t647 / 0.2e1 + t644 / 0.2e1 + t285 * t679 + t161 / 0.2e1 + t114 / 0.2e1 + t455) * t413) * qJD(1); m(4) * t140 + m(5) * t60 + m(6) * t30 + m(7) * t22; (t178 * t92 + t179 * t91 + t22 * t76) * t686 + (t126 * t213 + t127 * t212 + t30 * t97) * t687 + (t143 * t60 + t185 * t293 + t186 * t292) * t688 + t596 * t391 * t374 * t689 + (t329 * t702 - t15 + (-t191 * qJD(1) + (-qJD(1) * t473 - t272) * t415) * t415 - t706) * t415 + (t14 + t328 * t702 + (t192 * qJD(1) + (t472 * qJD(1) + t271) * t413) * t413 + ((-t272 + (-t642 - t645) * qJD(3) + t322 * t586 + t324 * t587 - t594) * t413 + (t321 * t586 + t323 * t587 + t271 - (t644 + t647) * qJD(3)) * t415 + (t193 - t190 + (t320 - t473) * t413 + t700) * qJD(1)) * t415 + t707) * t413 + (-t190 * t415 + t191 * t413 - t487 - t493 - t497) * t593 + (-t192 * t415 + t193 * t413 - t485 - t491 - t495) * t592; t430 + (t460 * t415 + t461 * t413 + (t413 * t454 + t415 * t455) * qJD(1)) * t425 + (-t85 + (t413 * t455 - t415 * t454) * qJD(3) - t663) * t428 + m(7) * (t173 * t24 + t174 * t23 + t88 * t94 + t89 * t93) + m(6) * (t129 * t99 + t130 * t98 + t180 * t56 + t181 * t55) + m(5) * (t124 * t199 + t125 * t198 + t196 * t96 + t197 * t95) - t616; m(5) * t67 + m(6) * t25 + m(7) * t21; t438 + ((-t413 * t578 + t415 * t579) * qJD(1) + t706 * t683 + t707 * t682 + (t413 * t620 - t415 * t621) * qJD(3) / 0.2e1) * t425 + (t10 / 0.2e1 + t7 / 0.2e1 + t579 * t586) * t413 + (-t11 / 0.2e1 - t8 / 0.2e1 + t578 * t586) * t415 + m(7) * (t178 * t23 + t179 * t24 + t21 * t76 + t22 * t80 + t91 * t93 + t92 * t94) + m(6) * (t105 * t30 + t126 * t129 + t127 * t130 + t212 * t55 + t213 * t56 + t25 * t97) + m(5) * (t143 * t67 + t177 * t60 + t185 * t198 + t186 * t199 + t292 * t95 + t293 * t96) + (qJD(1) * t708 + t704 * t413 + t705 * t415) * t681 + (t703 * t413 + t669 * t415) * qJD(1) / 0.2e1; (t21 * t80 + t23 * t94 + t24 * t93) * t686 + (t105 * t25 + t129 * t56 + t130 * t55) * t687 + (t177 * t67 + t198 * t96 + t199 * t95) * t688 + (-t9 + t663 * t428 + (t458 * t413 + t415 * t690) * qJD(3) + t616) * t428 + ((-t428 * t704 + t10 + t7) * t415 + (t705 * t428 + t11 + t8) * t413 + (t425 * t708 + t428 * t714) * qJD(3) + (t458 * t415 + (-t53 - t690) * t413) * qJD(1)) * t425 + t524; 0.2e1 * ((t173 * t415 + t174 * t413) * t684 + (t180 * t415 + t181 * t413) * t685) * t586 + 0.2e1 * ((-t173 * t593 + t174 * t592 + t413 * t88 + t415 * t89) * t684 + (-t180 * t593 + t181 * t592 + t413 * t98 + t415 * t99) * t685) * t425; (m(6) + m(7)) * t587; 0.2e1 * ((t178 * t589 + t179 * t588 - t22) * t684 + (t212 * t589 + t213 * t588 - t30) * t685) * t428 + 0.2e1 * ((qJD(3) * t76 + t178 * t592 - t179 * t593 + t413 * t92 + t415 * t91) * t684 + (qJD(3) * t97 + t126 * t415 + t127 * t413 + t212 * t592 - t213 * t593) * t685) * t425; 0.2e1 * ((t588 * t93 + t589 * t94 - t21) * t684 + (t129 * t588 + t130 * t589 - t25) * t685) * t428 + 0.2e1 * ((qJD(3) * t80 + t23 * t413 + t24 * t415 + t592 * t94 - t593 * t93) * t684 + (qJD(3) * t105 - t129 * t593 + t130 * t592 + t413 * t55 + t415 * t56) * t685) * t425; 0.4e1 * (t685 + t684) * (-0.1e1 + t596) * t425 * t586; t430 - t85 * t428 + m(7) * (t173 * t84 + t174 * t83 + t187 * t89 + t188 * t88); m(7) * t54; t438 + m(7) * (t156 * t22 + t178 * t83 + t179 * t84 + t187 * t91 + t188 * t92 + t54 * t76); m(7) * (t156 * t21 + t187 * t24 + t188 * t23 + t54 * t80 + t83 * t94 + t84 * t93) + t431; m(7) * ((-t54 + (t187 * t415 + t188 * t413) * qJD(3)) * t428 + (qJD(3) * t156 + t413 * t83 + t415 * t84 + (-t187 * t413 + t188 * t415) * qJD(1)) * t425); (t156 * t54 + t187 * t84 + t188 * t83) * t686 + t431;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
