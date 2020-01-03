% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:30
% EndTime: 2019-12-31 20:51:44
% DurationCPUTime: 9.74s
% Computational Cost: add. (26157->465), mult. (26886->577), div. (0->0), fcn. (24224->6), ass. (0->294)
t710 = Icges(5,1) + Icges(6,1);
t451 = sin(qJ(3));
t447 = Icges(5,5) * t451;
t448 = Icges(6,4) * t451;
t453 = cos(qJ(3));
t709 = t447 + t448 + (-Icges(6,2) - Icges(5,3)) * t453;
t449 = Icges(4,4) * t453;
t406 = Icges(4,1) * t451 + t449;
t581 = Icges(5,5) * t453;
t582 = Icges(6,4) * t453;
t708 = t710 * t451 + t406 - t581 - t582;
t583 = Icges(4,4) * t451;
t700 = Icges(4,2) * t453 + t583 - t709;
t450 = qJ(1) + qJ(2);
t445 = sin(t450);
t446 = cos(t450);
t625 = rSges(6,1) + pkin(4);
t501 = pkin(3) + t625;
t593 = rSges(6,2) + qJ(4);
t669 = t593 * t451 + t501 * t453 + pkin(2);
t681 = rSges(6,3) + qJ(5);
t215 = t669 * t446 + (pkin(7) - t681) * t445;
t604 = cos(qJ(1)) * pkin(1);
t213 = t215 + t604;
t558 = t446 * t451;
t203 = t213 * t558;
t592 = rSges(5,3) + qJ(4);
t626 = rSges(5,1) + pkin(3);
t668 = t592 * t451 + t626 * t453 + pkin(2);
t249 = (rSges(5,2) + pkin(7)) * t445 + t668 * t446;
t238 = t249 + t604;
t220 = t238 * t558;
t441 = t446 * rSges(5,2);
t442 = t446 * pkin(7);
t248 = -t668 * t445 + t441 + t442;
t605 = sin(qJ(1)) * pkin(1);
t237 = t248 - t605;
t494 = t249 * t558;
t495 = t215 * t558;
t214 = -t669 * t445 - t681 * t446 + t442;
t212 = t214 - t605;
t536 = -t212 - t214;
t561 = t445 * t451;
t659 = m(6) / 0.2e1;
t660 = m(5) / 0.2e1;
t600 = (t536 * t561 + t203 + t495) * t659 + (t220 + t494 + (-t237 - t248) * t561) * t660;
t696 = t237 - t248;
t698 = t212 - t214;
t601 = (-t561 * t698 + t203 - t495) * t659 + (-t696 * t561 + t220 - t494) * t660;
t6 = t601 - t600;
t707 = t6 * qJD(1);
t395 = Icges(5,3) * t451 + t581;
t313 = -Icges(5,6) * t446 + t395 * t445;
t398 = Icges(6,2) * t451 + t582;
t317 = Icges(6,6) * t446 + t398 * t445;
t706 = t313 + t317;
t675 = Icges(6,1) * t453 + t448;
t323 = Icges(6,5) * t446 + t445 * t675;
t676 = Icges(5,1) * t453 + t447;
t325 = -Icges(5,4) * t446 + t445 * t676;
t705 = -t323 - t325;
t704 = t395 + t398;
t401 = -Icges(4,2) * t451 + t449;
t702 = t401 + t708;
t407 = Icges(4,1) * t453 - t583;
t701 = t407 + t675 + t676;
t408 = pkin(3) * t451 - qJ(4) * t453;
t594 = rSges(5,1) * t451;
t410 = -rSges(5,3) * t453 + t594;
t507 = t408 + t410;
t304 = t507 * t445;
t306 = t507 * t446;
t487 = -t304 * t558 + t306 * t561;
t409 = rSges(6,1) * t451 - rSges(6,2) * t453;
t603 = pkin(4) * t451;
t490 = t408 + t409 + t603;
t271 = t490 * t445;
t273 = t490 * t446;
t488 = -t271 * t558 + t273 * t561;
t557 = t446 * t453;
t560 = t445 * t453;
t532 = t248 * t557 + t249 * t560;
t537 = t214 * t557 + t215 * t560;
t591 = (t487 + t532) * t660 + (t488 + t537) * t659;
t428 = pkin(3) * t561;
t280 = t428 + (-t592 * t453 + t594) * t445;
t418 = qJ(4) * t557;
t426 = rSges(5,3) * t557;
t281 = -t626 * t558 + t418 + t426;
t529 = t280 * t558 + t281 * t561;
t258 = t428 + (t625 * t451 - t593 * t453) * t445;
t427 = rSges(6,2) * t557;
t259 = -t501 * t558 + t418 + t427;
t531 = t258 * t558 + t259 * t561;
t598 = (t531 + t537) * t659 + (t529 + t532) * t660;
t12 = t598 - t591;
t533 = t237 * t557 + t238 * t560;
t538 = t212 * t557 + t213 * t560;
t597 = (t488 + t538) * t659 + (t487 + t533) * t660;
t599 = (t531 + t538) * t659 + (t529 + t533) * t660;
t9 = t599 - t597;
t699 = -t9 * qJD(1) - t12 * qJD(2);
t444 = t446 ^ 2;
t697 = -t213 + t215;
t443 = t445 ^ 2;
t503 = t443 + t444;
t694 = (-Icges(4,6) + Icges(5,6) - Icges(6,6)) * t453 + (-Icges(5,4) - Icges(4,5) + Icges(6,5)) * t451;
t324 = -Icges(6,5) * t445 + t446 * t675;
t326 = Icges(5,4) * t445 + t446 * t676;
t328 = Icges(4,5) * t445 + t407 * t446;
t693 = -t700 * t446 + t324 + t326 + t328;
t423 = Icges(4,4) * t561;
t327 = Icges(4,1) * t560 - Icges(4,5) * t446 - t423;
t692 = -Icges(4,2) * t560 + t709 * t445 + t327 - t423 - t705;
t421 = Icges(5,5) * t557;
t314 = Icges(5,6) * t445 + Icges(5,3) * t558 + t421;
t422 = Icges(6,4) * t557;
t318 = Icges(6,2) * t558 - Icges(6,6) * t445 + t422;
t322 = Icges(4,6) * t445 + t401 * t446;
t691 = -t406 * t446 - t710 * t558 + t314 + t318 - t322 + t421 + t422;
t321 = Icges(4,4) * t560 - Icges(4,2) * t561 - Icges(4,6) * t446;
t690 = t708 * t445 + t321 - t706;
t661 = m(4) / 0.2e1;
t623 = m(3) * (t604 * (-rSges(3,1) * t445 - rSges(3,2) * t446) + (t446 * rSges(3,1) - t445 * rSges(3,2)) * t605);
t399 = Icges(5,4) * t453 + Icges(5,6) * t451;
t566 = t399 * t445;
t319 = -Icges(5,2) * t446 + t566;
t301 = t445 * t319;
t176 = t313 * t558 + t325 * t557 + t301;
t685 = t176 * t446;
t684 = (-t700 + t701) * t453 + (-t702 + t704) * t451;
t680 = (t313 * t451 + t325 * t453) * t445;
t550 = t451 * t453;
t508 = t503 * t550;
t679 = (m(5) / 0.4e1 + m(6) / 0.4e1) * (t508 - t550);
t363 = t503 * t451;
t266 = m(6) * t363;
t502 = qJD(1) + qJD(2);
t678 = t502 * t266;
t69 = -t215 * t212 + t213 * t214;
t91 = -t249 * t237 + t238 * t248;
t595 = rSges(4,1) * t453;
t493 = pkin(2) + t595;
t505 = rSges(4,2) * t561 + t446 * rSges(4,3);
t278 = -t445 * t493 + t442 + t505;
t269 = t278 - t605;
t425 = rSges(4,2) * t558;
t279 = -t425 + t493 * t446 + (rSges(4,3) + pkin(7)) * t445;
t270 = t279 + t604;
t124 = -t279 * t269 + t270 * t278;
t414 = rSges(6,1) * t453 + rSges(6,2) * t451;
t677 = pkin(4) * t453 + t414;
t674 = t694 * t445;
t673 = t694 * t446;
t672 = t692 * t451 + t690 * t453;
t671 = -t693 * t451 + t691 * t453;
t411 = rSges(4,1) * t451 + rSges(4,2) * t453;
t499 = (t698 * t271 + t697 * t273) * t659 + ((-t238 + t249) * t306 + t696 * t304) * t660 + ((-t270 + t279) * t446 + (t269 - t278) * t445) * t411 * t661;
t118 = t237 * t280 + t238 * t281;
t120 = t248 * t280 + t249 * t281;
t358 = t411 * t445;
t360 = t411 * t446;
t141 = t269 * t358 - t270 * t360;
t151 = t278 * t358 - t279 * t360;
t87 = t212 * t258 + t213 * t259;
t90 = t214 * t258 + t215 * t259;
t670 = (t90 + t87) * t659 + (t120 + t118) * t660 + (t151 + t141) * t661;
t458 = (-t704 / 0.2e1 + t702 / 0.2e1) * t453 + (-t700 / 0.2e1 + t701 / 0.2e1) * t451;
t666 = 0.4e1 * qJD(1);
t664 = 0.4e1 * qJD(2);
t663 = 2 * qJD(3);
t651 = m(5) * t91;
t642 = m(6) * (t697 * t445 - t446 * t698);
t641 = m(6) * (t536 * t446 + (-t213 - t215) * t445);
t640 = m(6) * t69;
t413 = pkin(3) * t453 + qJ(4) * t451;
t510 = t503 * t413;
t130 = (pkin(4) * t560 + t414 * t445) * t445 + t510 + t677 * t444;
t530 = -t271 * t560 - t273 * t557;
t636 = m(6) * (t130 * t363 + t530);
t632 = m(6) * t87;
t631 = m(6) * t90;
t630 = -t445 / 0.2e1;
t629 = t445 / 0.2e1;
t628 = -t446 / 0.2e1;
t393 = Icges(6,5) * t453 + Icges(6,6) * t451;
t569 = t393 * t445;
t311 = Icges(6,3) * t446 + t569;
t474 = t317 * t451 + t323 * t453;
t168 = t311 * t446 + t445 * t474;
t568 = t393 * t446;
t312 = -Icges(6,3) * t445 + t568;
t169 = t446 * t312 + t318 * t561 + t324 * t560;
t108 = -t168 * t446 + t169 * t445;
t572 = t319 * t446;
t170 = -t572 + t680;
t565 = t399 * t446;
t320 = Icges(5,2) * t445 + t565;
t486 = -t314 * t561 + t320 * t446 - t326 * t560;
t109 = -t170 * t446 - t445 * t486;
t286 = t328 * t560;
t396 = Icges(4,5) * t453 - Icges(4,6) * t451;
t567 = t396 * t446;
t316 = Icges(4,3) * t445 + t567;
t492 = t316 * t446 - t286;
t173 = -t322 * t561 - t492;
t315 = Icges(4,5) * t560 - Icges(4,6) * t561 - Icges(4,3) * t446;
t571 = t321 * t451;
t110 = -(-t445 * (-t327 * t453 + t571) - t315 * t446) * t446 + t173 * t445;
t527 = -t317 * t558 - t323 * t557;
t574 = t311 * t445;
t174 = -t527 - t574;
t526 = t318 * t558 + t324 * t557;
t111 = -t174 * t446 + (-t312 * t445 + t526) * t445;
t177 = t314 * t558 + t445 * t320 + t326 * t557;
t112 = t177 * t445 - t685;
t525 = -t445 * t315 - t327 * t557;
t178 = -t321 * t558 - t525;
t524 = t445 * t316 + t328 * t557;
t179 = -t322 * t558 + t524;
t113 = -t178 * t446 + t179 * t445;
t29 = (-t169 + t174 + t574) * t445 + t311 * t444;
t470 = t177 + t572;
t30 = (t177 - t470) * t446 + (t176 + t486 - t301) * t445;
t491 = t322 * t451 - t315;
t31 = (t446 * t491 + t179 - t524) * t446 + (t445 * t491 + t178 + t492) * t445;
t32 = t527 * t446 + (t168 + (-t312 - t474) * t445 + t526) * t445;
t33 = -t685 + (t170 + t470 - t680) * t445;
t34 = (t173 - t286 + (t316 + t571) * t446 + t525) * t446 + t524 * t445;
t2 = (t113 / 0.2e1 - t32 / 0.2e1 - t34 / 0.2e1 - t33 / 0.2e1 + t111 / 0.2e1 + t112 / 0.2e1) * t446 + (t31 / 0.2e1 + t108 / 0.2e1 + t110 / 0.2e1 + t109 / 0.2e1 + t29 / 0.2e1 + t30 / 0.2e1) * t445;
t606 = m(6) * (t271 * t445 + t273 * t446);
t184 = t606 / 0.2e1;
t607 = m(6) * (-t258 * t445 + t259 * t446);
t67 = t184 - t607 / 0.2e1;
t624 = t2 * qJD(3) + t67 * qJD(5);
t621 = m(4) * t124;
t619 = m(4) * t141;
t618 = m(4) * t151;
t415 = rSges(5,1) * t453 + rSges(5,3) * t451;
t163 = t445 * (t415 * t445 - t441) + (t445 * rSges(5,2) + t415 * t446) * t446 + t510;
t528 = -t304 * t560 - t306 * t557;
t615 = m(5) * (t163 * t363 + t528);
t614 = m(5) * t118;
t612 = m(5) * t120;
t611 = m(5) * (-t237 * t561 + t220);
t610 = m(5) * (-t248 * t561 + t494);
t609 = m(6) * (-t212 * t561 + t203);
t608 = m(6) * (-t214 * t561 + t495);
t596 = m(6) * qJD(5);
t165 = t607 / 0.2e1;
t65 = t165 - t606 / 0.2e1;
t588 = t65 * qJD(3) + t266 * qJD(4);
t66 = t184 + t165;
t587 = t66 * qJD(3);
t539 = -t273 * t258 - t271 * t259;
t534 = -t306 * t280 - t304 * t281;
t511 = t445 * (qJ(4) * t560 - t428) + t446 * (-pkin(3) * t558 + t418);
t506 = -t413 - t415;
t123 = -t212 * t446 - t213 * t445;
t497 = m(6) * t123 * qJD(1);
t129 = -t214 * t446 - t215 * t445;
t496 = m(6) * t129 * qJD(2);
t489 = -t413 - t677;
t468 = (-t358 * t446 + t360 * t445) * t411;
t457 = t458 + t670;
t456 = -t458 + ((t321 + t706) * t453 + (t327 + t705) * t451) * (t629 + t630);
t455 = t66 * qJD(5) + ((t34 + t33 + t32) * t446 / 0.2e1 + (t110 + t109 + t108 + t31 + t30 + t29) * t630 + (t396 * t445 + t684 * t446 + t691 * t451 + t693 * t453 + t566 - t569) * t629 + (t684 * t445 - t690 * t451 + t692 * t453 + t111 + t112 + t113 - t565 - t567 + t568) * t628) * qJD(3);
t416 = -rSges(4,2) * t451 + t595;
t307 = t506 * t446;
t305 = t506 * t445;
t274 = t489 * t446;
t272 = t489 * t445;
t265 = t266 * qJD(5);
t227 = 0.4e1 * t679;
t182 = t446 * (-rSges(5,1) * t558 + t426) - t410 * t443 + t511;
t150 = t446 * (-rSges(6,1) * t558 + t427) - t409 * t443 - t503 * t603 + t511;
t68 = t608 + t610;
t60 = t609 + t611;
t57 = t641 / 0.2e1;
t56 = t642 / 0.2e1;
t36 = t615 + t636;
t22 = t458 + t612 + t618 + t631;
t19 = t458 + t614 + t619 + t632;
t18 = t621 + t623 + t640 + t651;
t17 = t57 - t642 / 0.2e1;
t16 = t57 + t56;
t15 = t56 - t641 / 0.2e1;
t13 = t591 + t598;
t10 = t597 + t599;
t7 = t600 + t601;
t5 = t457 + t499;
t4 = t457 - t499;
t3 = t456 + t499 - t670;
t1 = [t18 * qJD(2) + t19 * qJD(3) + t60 * qJD(4) + t123 * t596, t18 * qJD(1) + t5 * qJD(3) + t7 * qJD(4) + t16 * qJD(5) + 0.2e1 * (t623 / 0.2e1 + t124 * t661 + t69 * t659 + t91 * t660) * qJD(2), t19 * qJD(1) + t5 * qJD(2) + t10 * qJD(4) + (((-t269 * t446 - t270 * t445) * t416 + t468) * t661 + (t237 * t307 + t238 * t305 + t534) * t660 + (t212 * t274 + t213 * t272 + t539) * t659) * t663 + t455, qJD(1) * t60 + qJD(2) * t7 + qJD(3) * t10, t16 * qJD(2) + t497 + t587; t4 * qJD(3) - t6 * qJD(4) + t17 * qJD(5) + (-t640 / 0.4e1 - t651 / 0.4e1 - t621 / 0.4e1 - t623 / 0.4e1) * t666, t22 * qJD(3) + t68 * qJD(4) + t129 * t596, t4 * qJD(1) + t22 * qJD(2) + t13 * qJD(4) + (((-t278 * t446 - t279 * t445) * t416 + t468) * t661 + (t248 * t307 + t249 * t305 + t534) * t660 + (t214 * t274 + t215 * t272 + t539) * t659) * t663 + t455, qJD(2) * t68 + qJD(3) * t13 - t707, t17 * qJD(1) + t496 + t587; t456 * qJD(1) + t3 * qJD(2) - t9 * qJD(4) + (-t619 / 0.4e1 - t614 / 0.4e1 - t632 / 0.4e1) * t666 + t624, t3 * qJD(1) + t456 * qJD(2) - t12 * qJD(4) + (-t618 / 0.4e1 - t612 / 0.4e1 - t631 / 0.4e1) * t664 + t624, t36 * qJD(4) + t502 * t2 + (m(5) * (t163 * t182 - t304 * t305 - t306 * t307) + m(6) * (t130 * t150 - t271 * t272 - t273 * t274) + m(4) * ((t445 * (rSges(4,1) * t560 - t505) + t446 * (rSges(4,1) * t557 + t445 * rSges(4,3) - t425)) * (-t358 * t445 - t360 * t446) + t503 * t416 * t411) + (t673 * t443 + (t672 * t446 + (t671 - t674) * t445) * t446) * t629 + (t674 * t444 + (t671 * t445 + (t672 - t673) * t446) * t445) * t628) * qJD(3), t36 * qJD(3) + (-0.4e1 * t679 + 0.2e1 * (t660 + t659) * (-t363 * t453 + t508)) * qJD(4) + t699, t502 * t67; t6 * qJD(2) + t9 * qJD(3) - t265 + (-t609 / 0.4e1 - t611 / 0.4e1) * t666, t707 + t12 * qJD(3) - t265 + (-t610 / 0.4e1 - t608 / 0.4e1) * t664, t227 * qJD(4) + 0.4e1 * (-t615 / 0.4e1 - t636 / 0.4e1) * qJD(3) + ((-t182 * t453 + t528) * t660 + (-t150 * t453 + t530) * t659 + ((t305 * t445 + t307 * t446 + t163) * t660 + (t272 * t445 + t274 * t446 + t130) * t659) * t451) * t663 - t699, t227 * qJD(3), -t678; t15 * qJD(2) - t497 + t588, t15 * qJD(1) - t496 + t588, m(6) * (t272 * t446 - t274 * t445) * qJD(3) + t502 * t65, t678, 0;];
Cq = t1;
