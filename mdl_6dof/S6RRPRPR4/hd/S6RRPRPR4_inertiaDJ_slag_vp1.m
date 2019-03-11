% Calculate time derivative of joint inertia matrix for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:22
% EndTime: 2019-03-09 10:22:32
% DurationCPUTime: 42.93s
% Computational Cost: add. (120470->1362), mult. (275764->1802), div. (0->0), fcn. (324728->14), ass. (0->572)
t534 = sin(pkin(6));
t700 = sin(pkin(11));
t701 = cos(pkin(11));
t715 = sin(qJ(2));
t716 = cos(qJ(2));
t552 = t700 * t716 + t701 * t715;
t494 = t552 * t534;
t665 = qJ(4) + pkin(12);
t533 = sin(t665);
t612 = cos(t665);
t702 = cos(pkin(6));
t473 = t494 * t612 + t533 * t702;
t551 = t715 * t700 - t716 * t701;
t493 = t551 * t534;
t558 = -t494 * t533 + t612 * t702;
t350 = Icges(6,5) * t473 + Icges(6,6) * t558 + Icges(6,3) * t493;
t537 = sin(qJ(4));
t540 = cos(qJ(4));
t474 = -t494 * t537 + t540 * t702;
t622 = t702 * t537;
t571 = -t494 * t540 - t622;
t386 = -Icges(5,5) * t571 + Icges(5,6) * t474 + Icges(5,3) * t493;
t781 = t350 + t386;
t592 = t702 * t700;
t593 = t702 * t701;
t495 = t592 * t716 + t593 * t715;
t538 = sin(qJ(1));
t541 = cos(qJ(1));
t463 = t495 * t541 - t538 * t551;
t546 = -t592 * t715 + t593 * t716;
t485 = t546 * qJD(2);
t503 = t552 * qJD(2);
t374 = -qJD(1) * t463 - t538 * t485 - t503 * t541;
t465 = -t538 * t495 - t541 * t551;
t697 = t534 * t538;
t417 = t465 * t612 + t533 * t697;
t583 = t534 * t612;
t574 = qJD(1) * t583;
t279 = qJD(4) * t417 + t374 * t533 - t541 * t574;
t562 = -t465 * t533 + t538 * t583;
t669 = qJD(1) * t541;
t632 = t534 * t669;
t280 = qJD(4) * t562 + t374 * t612 + t533 * t632;
t462 = -t538 * t552 + t541 * t546;
t484 = t495 * qJD(2);
t549 = qJD(2) * t551;
t373 = qJD(1) * t462 - t538 * t484 - t541 * t549;
t167 = Icges(6,5) * t280 - Icges(6,6) * t279 + Icges(6,3) * t373;
t657 = t537 * t697;
t425 = t465 * t540 + t657;
t305 = -qJD(4) * t425 - t374 * t537 + t540 * t632;
t696 = t534 * t540;
t424 = -t465 * t537 + t538 * t696;
t604 = t537 * t632;
t306 = qJD(4) * t424 + t374 * t540 + t604;
t178 = Icges(5,5) * t306 + Icges(5,6) * t305 + Icges(5,3) * t373;
t780 = -t167 - t178;
t376 = qJD(1) * t465 + t485 * t541 - t538 * t503;
t695 = t534 * t541;
t415 = t463 * t612 - t533 * t695;
t281 = qJD(4) * t415 + t376 * t533 - t538 * t574;
t561 = -t463 * t533 - t541 * t583;
t670 = qJD(1) * t538;
t633 = t534 * t670;
t282 = qJD(4) * t561 + t376 * t612 + t533 * t633;
t544 = t538 * t546;
t375 = -qJD(1) * t544 - t541 * t484 + t538 * t549 - t552 * t669;
t168 = Icges(6,5) * t282 - Icges(6,6) * t281 - Icges(6,3) * t375;
t656 = t537 * t695;
t573 = -t463 * t540 + t656;
t307 = qJD(4) * t573 - t376 * t537 + t540 * t633;
t422 = -t463 * t537 - t540 * t695;
t603 = t537 * t633;
t308 = qJD(4) * t422 + t376 * t540 + t603;
t179 = Icges(5,5) * t308 + Icges(5,6) * t307 - Icges(5,3) * t375;
t779 = t168 + t179;
t288 = Icges(6,5) * t415 + Icges(6,6) * t561 - Icges(6,3) * t462;
t311 = -Icges(5,5) * t573 + Icges(5,6) * t422 - Icges(5,3) * t462;
t778 = t288 + t311;
t464 = -t541 * t552 - t544;
t289 = Icges(6,5) * t417 + Icges(6,6) * t562 - Icges(6,3) * t464;
t312 = Icges(5,5) * t425 + Icges(5,6) * t424 - Icges(5,3) * t464;
t777 = t289 + t312;
t482 = qJD(2) * t494;
t351 = Icges(6,4) * t473 + Icges(6,2) * t558 + Icges(6,6) * t493;
t352 = Icges(6,1) * t473 + Icges(6,4) * t558 + Icges(6,5) * t493;
t387 = -Icges(5,4) * t571 + Icges(5,2) * t474 + Icges(5,6) * t493;
t388 = -Icges(5,1) * t571 + Icges(5,4) * t474 + Icges(5,5) * t493;
t751 = t351 * t558 + t352 * t473 + t387 * t474 - t388 * t571 + t493 * t781;
t776 = t482 * t751;
t668 = qJD(2) * t534;
t483 = t551 * t668;
t404 = qJD(4) * t473 - t483 * t533;
t405 = qJD(4) * t558 - t483 * t612;
t300 = Icges(6,5) * t405 - Icges(6,6) * t404 + Icges(6,3) * t482;
t301 = Icges(6,4) * t405 - Icges(6,2) * t404 + Icges(6,6) * t482;
t302 = Icges(6,1) * t405 - Icges(6,4) * t404 + Icges(6,5) * t482;
t412 = qJD(4) * t571 + t483 * t537;
t567 = t474 * qJD(4);
t413 = -t483 * t540 + t567;
t321 = Icges(5,5) * t413 + Icges(5,6) * t412 + Icges(5,3) * t482;
t322 = Icges(5,4) * t413 + Icges(5,2) * t412 + Icges(5,6) * t482;
t323 = Icges(5,1) * t413 + Icges(5,4) * t412 + Icges(5,5) * t482;
t775 = t301 * t558 + t473 * t302 + t474 * t322 - t323 * t571 - t404 * t351 + t405 * t352 + t412 * t387 + t413 * t388 + (t300 + t321) * t493 + t781 * t482;
t170 = Icges(6,4) * t282 - Icges(6,2) * t281 - Icges(6,6) * t375;
t172 = Icges(6,1) * t282 - Icges(6,4) * t281 - Icges(6,5) * t375;
t181 = Icges(5,4) * t308 + Icges(5,2) * t307 - Icges(5,6) * t375;
t183 = Icges(5,1) * t308 + Icges(5,4) * t307 - Icges(5,5) * t375;
t290 = Icges(6,4) * t415 + Icges(6,2) * t561 - Icges(6,6) * t462;
t292 = Icges(6,1) * t415 + Icges(6,4) * t561 - Icges(6,5) * t462;
t313 = -Icges(5,4) * t573 + Icges(5,2) * t422 - Icges(5,6) * t462;
t315 = -Icges(5,1) * t573 + Icges(5,4) * t422 - Icges(5,5) * t462;
t774 = -t170 * t562 - t172 * t417 - t181 * t424 - t183 * t425 + t279 * t290 - t280 * t292 - t305 * t313 - t306 * t315 - t778 * t373 + t779 * t464;
t169 = Icges(6,4) * t280 - Icges(6,2) * t279 + Icges(6,6) * t373;
t171 = Icges(6,1) * t280 - Icges(6,4) * t279 + Icges(6,5) * t373;
t180 = Icges(5,4) * t306 + Icges(5,2) * t305 + Icges(5,6) * t373;
t182 = Icges(5,1) * t306 + Icges(5,4) * t305 + Icges(5,5) * t373;
t291 = Icges(6,4) * t417 + Icges(6,2) * t562 - Icges(6,6) * t464;
t293 = Icges(6,1) * t417 + Icges(6,4) * t562 - Icges(6,5) * t464;
t314 = Icges(5,4) * t425 + Icges(5,2) * t424 - Icges(5,6) * t464;
t316 = Icges(5,1) * t425 + Icges(5,4) * t424 - Icges(5,5) * t464;
t773 = t169 * t562 + t171 * t417 + t180 * t424 + t182 * t425 - t279 * t291 + t280 * t293 + t305 * t314 + t306 * t316 + t777 * t373 + t780 * t464;
t772 = -t170 * t561 - t172 * t415 - t181 * t422 + t183 * t573 + t281 * t290 - t282 * t292 - t307 * t313 - t308 * t315 + t778 * t375 + t779 * t462;
t771 = t169 * t561 + t171 * t415 + t180 * t422 - t182 * t573 - t281 * t291 + t282 * t293 + t307 * t314 + t308 * t316 - t777 * t375 + t780 * t462;
t64 = t168 * t493 + t170 * t558 + t172 * t473 + t288 * t482 - t290 * t404 + t292 * t405;
t67 = t179 * t493 + t181 * t474 - t183 * t571 + t311 * t482 + t313 * t412 + t315 * t413;
t770 = -t64 - t67;
t65 = t167 * t493 + t169 * t558 + t171 * t473 + t289 * t482 - t291 * t404 + t293 * t405;
t68 = t178 * t493 + t180 * t474 - t182 * t571 + t312 * t482 + t314 * t412 + t316 * t413;
t769 = t65 + t68;
t84 = -t279 * t351 + t280 * t352 - t300 * t464 + t301 * t562 + t302 * t417 + t350 * t373;
t86 = t305 * t387 + t306 * t388 - t321 * t464 + t322 * t424 + t323 * t425 + t373 * t386;
t768 = t84 + t86;
t85 = -t281 * t351 + t282 * t352 - t300 * t462 + t301 * t561 + t302 * t415 - t350 * t375;
t87 = t307 * t387 + t308 * t388 - t321 * t462 + t322 * t422 - t323 * t573 - t375 * t386;
t767 = t87 + t85;
t766 = t290 * t561 + t292 * t415 + t313 * t422 - t315 * t573 - t778 * t462;
t765 = t291 * t561 + t293 * t415 + t314 * t422 - t316 * t573 - t777 * t462;
t764 = t290 * t562 + t292 * t417 + t313 * t424 + t315 * t425 - t778 * t464;
t763 = t291 * t562 + t293 * t417 + t314 * t424 + t316 * t425 - t777 * t464;
t144 = t288 * t493 + t290 * t558 + t292 * t473;
t148 = t311 * t493 + t313 * t474 - t315 * t571;
t762 = t144 + t148;
t145 = t289 * t493 + t291 * t558 + t293 * t473;
t149 = t312 * t493 + t314 * t474 - t316 * t571;
t752 = t145 + t149;
t165 = -t350 * t462 + t351 * t561 + t352 * t415;
t186 = -t386 * t462 + t387 * t422 - t388 * t573;
t761 = t165 + t186;
t166 = -t350 * t464 + t351 * t562 + t352 * t417;
t187 = -t386 * t464 + t387 * t424 + t388 * t425;
t760 = t166 + t187;
t355 = Icges(4,5) * t463 + Icges(4,6) * t462 - Icges(4,3) * t695;
t597 = t702 * t715;
t555 = -t538 * t716 - t541 * t597;
t598 = t702 * t716;
t556 = t538 * t715 - t541 * t598;
t444 = -Icges(3,5) * t555 - Icges(3,6) * t556 - Icges(3,3) * t695;
t745 = t355 + t444;
t356 = Icges(4,5) * t465 + Icges(4,6) * t464 + Icges(4,3) * t697;
t506 = -t538 * t598 - t541 * t715;
t557 = t538 * t597 - t541 * t716;
t445 = -Icges(3,5) * t557 + Icges(3,6) * t506 + Icges(3,3) * t697;
t757 = t356 + t445;
t429 = -Icges(4,5) * t483 - Icges(4,6) * t482;
t430 = -Icges(4,4) * t483 - Icges(4,2) * t482;
t431 = -Icges(4,1) * t483 - Icges(4,4) * t482;
t436 = Icges(4,4) * t494 - Icges(4,2) * t493 + Icges(4,6) * t702;
t437 = Icges(4,1) * t494 - Icges(4,4) * t493 + Icges(4,5) * t702;
t634 = t715 * Icges(3,4);
t488 = Icges(3,6) * t702 + (Icges(3,2) * t716 + t634) * t534;
t635 = t716 * Icges(3,4);
t489 = Icges(3,5) * t702 + (Icges(3,1) * t715 + t635) * t534;
t496 = (Icges(3,5) * t716 - Icges(3,6) * t715) * t668;
t497 = (-Icges(3,2) * t715 + t635) * t668;
t498 = (Icges(3,1) * t716 - t634) * t668;
t637 = t534 * t716;
t602 = qJD(2) * t637;
t627 = qJD(2) * t715;
t636 = t534 * t715;
t759 = -t488 * t534 * t627 - t493 * t430 + t494 * t431 - t482 * t436 - t483 * t437 + t489 * t602 + t497 * t637 + t498 * t636 + (t429 + t496) * t702;
t758 = t775 * t702;
t251 = Icges(4,5) * t374 - Icges(4,6) * t373 + Icges(4,3) * t632;
t468 = qJD(1) * t556 + qJD(2) * t557;
t469 = qJD(1) * t555 + qJD(2) * t506;
t378 = Icges(3,5) * t469 + Icges(3,6) * t468 + Icges(3,3) * t632;
t756 = t378 + t251;
t755 = t775 * t493 + t776;
t532 = pkin(2) * t716 + pkin(1);
t710 = pkin(1) - t532;
t754 = t538 * t710;
t536 = sin(qJ(6));
t539 = cos(qJ(6));
t344 = t417 * t539 - t464 * t536;
t199 = -qJD(6) * t344 - t280 * t536 + t373 * t539;
t343 = -t417 * t536 - t464 * t539;
t200 = qJD(6) * t343 + t280 * t539 + t373 * t536;
t115 = t200 * rSges(7,1) + t199 * rSges(7,2) + t279 * rSges(7,3);
t691 = t280 * pkin(5) + t279 * pkin(10) + t115;
t342 = t415 * t539 - t462 * t536;
t201 = -qJD(6) * t342 - t282 * t536 - t375 * t539;
t341 = -t415 * t536 - t462 * t539;
t202 = qJD(6) * t341 + t282 * t539 - t375 * t536;
t586 = -t202 * rSges(7,1) - t201 * rSges(7,2);
t116 = t281 * rSges(7,3) - t586;
t713 = t282 * pkin(5);
t690 = t281 * pkin(10) + t116 + t713;
t357 = Icges(4,4) * t463 + Icges(4,2) * t462 - Icges(4,6) * t695;
t359 = Icges(4,1) * t463 + Icges(4,4) * t462 - Icges(4,5) * t695;
t446 = -Icges(3,4) * t555 - Icges(3,2) * t556 - Icges(3,6) * t695;
t448 = -Icges(3,1) * t555 - Icges(3,4) * t556 - Icges(3,5) * t695;
t750 = -t462 * t357 - t463 * t359 + t446 * t556 + t448 * t555 + t695 * t745;
t358 = Icges(4,4) * t465 + Icges(4,2) * t464 + Icges(4,6) * t697;
t360 = Icges(4,1) * t465 + Icges(4,4) * t464 + Icges(4,5) * t697;
t447 = -Icges(3,4) * t557 + Icges(3,2) * t506 + Icges(3,6) * t697;
t449 = -Icges(3,1) * t557 + Icges(3,4) * t506 + Icges(3,5) * t697;
t749 = t358 * t464 + t360 * t465 + t447 * t506 - t449 * t557 + t697 * t757;
t748 = t759 * t702;
t585 = -t342 * rSges(7,1) - t341 * rSges(7,2);
t228 = -rSges(7,3) * t561 - t585;
t712 = t415 * pkin(5);
t687 = -pkin(10) * t561 + t228 + t712;
t229 = t344 * rSges(7,1) + t343 * rSges(7,2) - rSges(7,3) * t562;
t686 = t417 * pkin(5) - pkin(10) * t562 + t229;
t747 = -t462 * t358 - t463 * t360 + t447 * t556 + t449 * t555 + t695 * t757;
t746 = t357 * t464 + t359 * t465 + t446 * t506 - t448 * t557 + t697 * t745;
t530 = pkin(8) * t695;
t744 = -t538 * pkin(1) + t530;
t406 = -t473 * t536 + t493 * t539;
t407 = t473 * t539 + t493 * t536;
t284 = Icges(7,5) * t407 + Icges(7,6) * t406 - Icges(7,3) * t558;
t285 = Icges(7,4) * t407 + Icges(7,2) * t406 - Icges(7,6) * t558;
t286 = Icges(7,1) * t407 + Icges(7,4) * t406 - Icges(7,5) * t558;
t126 = -t284 * t562 + t285 * t343 + t286 * t344;
t110 = Icges(7,5) * t202 + Icges(7,6) * t201 + Icges(7,3) * t281;
t112 = Icges(7,4) * t202 + Icges(7,2) * t201 + Icges(7,6) * t281;
t114 = Icges(7,1) * t202 + Icges(7,4) * t201 + Icges(7,5) * t281;
t222 = Icges(7,5) * t342 + Icges(7,6) * t341 - Icges(7,3) * t561;
t224 = Icges(7,4) * t342 + Icges(7,2) * t341 - Icges(7,6) * t561;
t226 = Icges(7,1) * t342 + Icges(7,4) * t341 - Icges(7,5) * t561;
t24 = -t110 * t562 + t112 * t343 + t114 * t344 + t199 * t224 + t200 * t226 + t222 * t279;
t109 = Icges(7,5) * t200 + Icges(7,6) * t199 + Icges(7,3) * t279;
t111 = Icges(7,4) * t200 + Icges(7,2) * t199 + Icges(7,6) * t279;
t113 = Icges(7,1) * t200 + Icges(7,4) * t199 + Icges(7,5) * t279;
t223 = Icges(7,5) * t344 + Icges(7,6) * t343 - Icges(7,3) * t562;
t225 = Icges(7,4) * t344 + Icges(7,2) * t343 - Icges(7,6) * t562;
t227 = Icges(7,1) * t344 + Icges(7,4) * t343 - Icges(7,5) * t562;
t25 = -t109 * t562 + t111 * t343 + t113 * t344 + t199 * t225 + t200 * t227 + t223 * t279;
t309 = -qJD(6) * t407 - t405 * t536 + t482 * t539;
t310 = qJD(6) * t406 + t405 * t539 + t482 * t536;
t190 = Icges(7,5) * t310 + Icges(7,6) * t309 + Icges(7,3) * t404;
t191 = Icges(7,4) * t310 + Icges(7,2) * t309 + Icges(7,6) * t404;
t192 = Icges(7,1) * t310 + Icges(7,4) * t309 + Icges(7,5) * t404;
t42 = -t190 * t562 + t191 * t343 + t192 * t344 + t199 * t285 + t200 * t286 + t279 * t284;
t92 = -t222 * t562 + t224 * t343 + t226 * t344;
t93 = -t223 * t562 + t225 * t343 + t227 * t344;
t3 = t126 * t482 - t24 * t462 - t25 * t464 + t373 * t93 - t375 * t92 + t42 * t493;
t743 = t763 * t373 - t764 * t375 + t462 * t774 - t773 * t464 + t760 * t482 + t768 * t493 + t3;
t125 = -t284 * t561 + t285 * t341 + t286 * t342;
t26 = -t110 * t561 + t112 * t341 + t114 * t342 + t201 * t224 + t202 * t226 + t222 * t281;
t27 = -t109 * t561 + t111 * t341 + t113 * t342 + t201 * t225 + t202 * t227 + t223 * t281;
t43 = -t190 * t561 + t191 * t341 + t192 * t342 + t201 * t285 + t202 * t286 + t281 * t284;
t90 = -t222 * t561 + t224 * t341 + t226 * t342;
t91 = -t223 * t561 + t225 * t341 + t227 * t342;
t4 = t125 * t482 - t26 * t462 - t27 * t464 + t373 * t91 - t375 * t90 + t43 * t493;
t742 = t373 * t765 - t375 * t766 + t462 * t772 - t464 * t771 + t482 * t761 + t493 * t767 + t4;
t5 = t42 * t702 + (-t24 * t541 + t25 * t538 + (t538 * t92 + t541 * t93) * qJD(1)) * t534;
t741 = t5 + t768 * t702 + (t774 * t541 + t773 * t538 + (t538 * t764 + t541 * t763) * qJD(1)) * t534;
t6 = t43 * t702 + (-t26 * t541 + t27 * t538 + (t538 * t90 + t541 * t91) * qJD(1)) * t534;
t740 = t6 + t767 * t702 + (t772 * t541 + t771 * t538 + (t538 * t766 + t541 * t765) * qJD(1)) * t534;
t102 = -t222 * t558 + t224 * t406 + t226 * t407;
t103 = -t223 * t558 + t225 * t406 + t227 * t407;
t28 = -t110 * t558 + t112 * t406 + t114 * t407 + t222 * t404 + t224 * t309 + t226 * t310;
t29 = -t109 * t558 + t111 * t406 + t113 * t407 + t223 * t404 + t225 * t309 + t227 * t310;
t132 = -t284 * t558 + t285 * t406 + t286 * t407;
t59 = -t190 * t558 + t406 * t191 + t407 * t192 + t404 * t284 + t309 * t285 + t310 * t286;
t704 = t132 * t482 + t59 * t493;
t8 = -t102 * t375 + t103 * t373 - t28 * t462 - t29 * t464 + t704;
t739 = t373 * t752 - t375 * t762 + t462 * t770 - t464 * t769 + t755 + t8;
t58 = t59 * t702;
t9 = t58 + (-t28 * t541 + t29 * t538 + (t102 * t538 + t103 * t541) * qJD(1)) * t534;
t738 = t9 + t758 + (t770 * t541 + t769 * t538 + (t538 * t762 + t541 * t752) * qJD(1)) * t534;
t38 = t126 * t493 - t462 * t92 - t464 * t93;
t737 = -t462 * t764 - t464 * t763 + t493 * t760 + t38;
t40 = t125 * t702 + (t538 * t91 - t541 * t90) * t534;
t736 = t40 + t761 * t702 + (t538 * t765 - t541 * t766) * t534;
t41 = t126 * t702 + (t538 * t93 - t541 * t92) * t534;
t735 = t41 + t760 * t702 + (t538 * t763 - t541 * t764) * t534;
t37 = t125 * t493 - t462 * t90 - t464 * t91;
t734 = -t462 * t766 - t464 * t765 + t493 * t761 + t37;
t252 = Icges(4,5) * t376 + Icges(4,6) * t375 + Icges(4,3) * t633;
t470 = qJD(1) * t506 + qJD(2) * t555;
t471 = -qJD(1) * t557 - qJD(2) * t556;
t379 = Icges(3,5) * t471 + Icges(3,6) * t470 + Icges(3,3) * t633;
t733 = (-t379 - t252) * t541;
t732 = 0.2e1 * t538;
t731 = m(7) / 0.2e1;
t730 = t279 / 0.2e1;
t729 = t281 / 0.2e1;
t728 = t373 / 0.2e1;
t727 = -t375 / 0.2e1;
t726 = t404 / 0.2e1;
t725 = -t561 / 0.2e1;
t724 = -t562 / 0.2e1;
t723 = -t462 / 0.2e1;
t722 = -t464 / 0.2e1;
t721 = -t558 / 0.2e1;
t720 = t482 / 0.2e1;
t719 = t493 / 0.2e1;
t718 = t538 / 0.2e1;
t717 = -rSges(7,3) - pkin(10);
t714 = pkin(1) * t541;
t531 = pkin(4) * t540 + pkin(3);
t709 = -pkin(3) + t531;
t535 = -qJ(5) - pkin(9);
t708 = -pkin(9) - t535;
t707 = pkin(4) * qJD(4);
t706 = rSges(6,3) - t535;
t705 = t132 * t404 - t558 * t59;
t699 = t375 * t535;
t698 = t462 * t535;
t694 = t538 * t532;
t500 = pkin(2) * t597 + (-pkin(8) - qJ(3)) * t534;
t693 = t541 * t500;
t523 = t541 * t532;
t266 = t374 * pkin(3) + t373 * pkin(9);
t609 = t696 * t707;
t660 = t537 * t707;
t550 = pkin(4) * t604 - t464 * qJD(5) - t373 * t535 + t374 * t531 - t465 * t660 + t538 * t609;
t159 = -t266 + t550;
t624 = -t465 * pkin(3) + pkin(9) * t464;
t520 = pkin(4) * t657;
t642 = t464 * t535 + t465 * t531 + t520;
t264 = t624 + t642;
t689 = t493 * t159 + t482 * t264;
t193 = rSges(7,1) * t310 + rSges(7,2) * t309 + rSges(7,3) * t404;
t688 = -pkin(5) * t405 - pkin(10) * t404 - t193;
t372 = t375 * pkin(9);
t641 = -t462 * qJD(5) - t463 * t660 - t541 * t609;
t160 = pkin(4) * t603 + t376 * t709 + t372 + t641 + t699;
t460 = t462 * pkin(9);
t521 = pkin(4) * t656;
t263 = t463 * t709 + t460 - t521 + t698;
t685 = -t464 * t160 + t373 * t263;
t304 = pkin(4) * t567 + t493 * qJD(5) + t482 * t708 - t483 * t709;
t349 = pkin(4) * t622 + t493 * t708 + t494 * t709;
t684 = -t462 * t304 - t375 * t349;
t587 = -t415 * rSges(6,1) - rSges(6,2) * t561;
t294 = -t462 * rSges(6,3) - t587;
t683 = -t263 - t294;
t295 = t417 * rSges(6,1) + rSges(6,2) * t562 - t464 * rSges(6,3);
t682 = t264 + t295;
t509 = pkin(2) * qJD(2) * t598 - t534 * qJD(3);
t608 = pkin(2) * t627;
t568 = -t538 * t509 - t541 * t608;
t625 = -t534 * pkin(8) - t500;
t397 = (t541 * t625 + t754) * qJD(1) + t568;
t394 = t702 * t397;
t681 = t702 * t266 + t394;
t287 = rSges(7,1) * t407 + rSges(7,2) * t406 - rSges(7,3) * t558;
t680 = pkin(5) * t473 - pkin(10) * t558 + t287;
t354 = rSges(6,1) * t473 + rSges(6,2) * t558 + rSges(6,3) * t493;
t679 = t349 + t354;
t455 = t538 * t625 + t523 - t714;
t442 = t702 * t455;
t678 = -t624 * t702 + t442;
t677 = t624 - t455;
t508 = pkin(2) * t602 + qJD(3) * t702;
t676 = rSges(4,1) * t483 + rSges(4,2) * t482 - t508;
t454 = t530 + t693 - t754;
t675 = t454 * t697 + t455 * t695;
t674 = t483 * pkin(3) - t482 * pkin(9) - t508;
t439 = t494 * rSges(4,1) - t493 * rSges(4,2) + rSges(4,3) * t702;
t512 = pkin(2) * t636 + qJ(3) * t702;
t673 = -t439 - t512;
t452 = t494 * pkin(3) + t493 * pkin(9);
t672 = -t452 - t512;
t671 = t500 * t670 + t538 * t608;
t664 = 0.2e1 * t702;
t663 = m(6) / 0.2e1 + t731;
t662 = 0.2e1 * qJD(1);
t661 = pkin(8) * t697;
t528 = rSges(4,3) * t697;
t659 = t28 / 0.2e1 + t43 / 0.2e1;
t658 = t42 / 0.2e1 + t29 / 0.2e1;
t654 = t702 / 0.2e1;
t653 = t702 * t159 + t681;
t652 = -t263 - t687;
t651 = t264 + t686;
t650 = t702 * t264 + t678;
t649 = -t264 + t677;
t173 = t280 * rSges(6,1) - t279 * rSges(6,2) + t373 * rSges(6,3);
t648 = t349 + t680;
t184 = t306 * rSges(5,1) + t305 * rSges(5,2) + t373 * rSges(5,3);
t647 = -t304 + t674;
t324 = rSges(5,1) * t413 + rSges(5,2) * t412 + rSges(5,3) * t482;
t646 = -t324 + t674;
t645 = -t349 + t672;
t257 = t374 * rSges(4,1) - t373 * rSges(4,2) + rSges(4,3) * t632;
t389 = -rSges(5,1) * t571 + rSges(5,2) * t474 + rSges(5,3) * t493;
t644 = -t389 + t672;
t594 = t509 * t541 - t671;
t398 = (-t541 * t710 - t661) * qJD(1) + t594;
t643 = t397 * t695 + t398 * t697 + t454 * t632;
t318 = t425 * rSges(5,1) + t424 * rSges(5,2) - t464 * rSges(5,3);
t362 = t465 * rSges(4,1) + t464 * rSges(4,2) + t528;
t384 = t469 * rSges(3,1) + t468 * rSges(3,2) + rSges(3,3) * t632;
t451 = -rSges(3,1) * t557 + t506 * rSges(3,2) + rSges(3,3) * t697;
t640 = m(5) * t702;
t639 = m(6) * t702;
t638 = m(7) * t702;
t629 = t102 / 0.2e1 + t125 / 0.2e1;
t628 = t126 / 0.2e1 + t103 / 0.2e1;
t267 = t376 * pkin(3) - t372;
t390 = t463 * pkin(3) - t460;
t623 = t702 * t454;
t621 = 2 * m(3);
t619 = 2 * m(4);
t617 = 0.2e1 * m(5);
t615 = 0.2e1 * m(6);
t613 = 0.2e1 * m(7);
t611 = t541 * t673;
t610 = -t538 * t500 + t523;
t303 = rSges(6,1) * t405 - rSges(6,2) * t404 + rSges(6,3) * t482;
t607 = -t303 + t647;
t606 = -t354 + t645;
t605 = t390 * t697 - t624 * t695 + t675;
t599 = t541 * t644;
t591 = -t471 * rSges(3,1) - t470 * rSges(3,2);
t590 = -t376 * rSges(4,1) - t375 * rSges(4,2);
t589 = -t463 * rSges(4,1) - t462 * rSges(4,2);
t588 = -t282 * rSges(6,1) + t281 * rSges(6,2);
t582 = -t398 * t702 + t512 * t633;
t581 = t647 + t688;
t580 = t645 - t680;
t578 = -t693 - t694;
t577 = t541 * t606;
t576 = t263 * t697 + t264 * t695 + t605;
t575 = t266 * t695 + t267 * t697 + t390 * t632 + t643;
t572 = t541 * t580;
t570 = t610 + t642;
t185 = t308 * rSges(5,1) + t307 * rSges(5,2) - t375 * rSges(5,3);
t317 = -rSges(5,1) * t573 + t422 * rSges(5,2) - t462 * rSges(5,3);
t569 = -t390 * t702 - t623;
t450 = -rSges(3,1) * t555 - rSges(3,2) * t556 - rSges(3,3) * t695;
t566 = -t87 / 0.2e1 - t85 / 0.2e1 - t67 / 0.2e1 - t64 / 0.2e1 - t659;
t565 = -t86 / 0.2e1 - t84 / 0.2e1 - t68 / 0.2e1 - t65 / 0.2e1 - t658;
t564 = -t463 * t531 + t521 + t578;
t563 = t159 * t695 + t160 * t697 + t263 * t632 + t575;
t560 = t187 / 0.2e1 + t166 / 0.2e1 + t149 / 0.2e1 + t145 / 0.2e1 + t628;
t559 = -t186 / 0.2e1 - t165 / 0.2e1 - t148 / 0.2e1 - t144 / 0.2e1 - t629;
t553 = -t267 * t702 + t452 * t633 + t582;
t548 = -t263 * t702 + t569;
t547 = -t160 * t702 + t349 * t633 + t553;
t545 = qJD(1) * t578 + t568;
t543 = (-t520 - t523) * qJD(1) - t376 * t531 - t594 - t641;
t542 = t545 + t550;
t499 = (rSges(3,1) * t716 - rSges(3,2) * t715) * t668;
t491 = t702 * rSges(3,3) + (rSges(3,1) * t715 + rSges(3,2) * t716) * t534;
t487 = Icges(3,3) * t702 + (Icges(3,5) * t715 + Icges(3,6) * t716) * t534;
t435 = Icges(4,5) * t494 - Icges(4,6) * t493 + Icges(4,3) * t702;
t421 = t451 + t661 + t714;
t420 = -t450 + t744;
t396 = -t450 * t702 - t491 * t695;
t395 = t451 * t702 - t491 * t697;
t385 = rSges(3,3) * t633 - t591;
t383 = Icges(3,1) * t471 + Icges(3,4) * t470 + Icges(3,5) * t633;
t382 = Icges(3,1) * t469 + Icges(3,4) * t468 + Icges(3,5) * t632;
t381 = Icges(3,4) * t471 + Icges(3,2) * t470 + Icges(3,6) * t633;
t380 = Icges(3,4) * t469 + Icges(3,2) * t468 + Icges(3,6) * t632;
t361 = -rSges(4,3) * t695 - t589;
t348 = (-t714 + (-rSges(3,3) - pkin(8)) * t697) * qJD(1) + t591;
t347 = qJD(1) * t744 + t384;
t346 = t487 * t697 + t488 * t506 - t489 * t557;
t345 = -t487 * t695 - t488 * t556 - t489 * t555;
t333 = t610 + t362;
t332 = -t694 + (rSges(4,3) * t534 - t500) * t541 + t589;
t329 = t462 * t349;
t328 = t702 * t384 + (-t491 * t669 - t499 * t538) * t534;
t327 = -t702 * t385 + (t491 * t670 - t499 * t541) * t534;
t320 = t702 * t445 + (t447 * t716 + t449 * t715) * t534;
t319 = t702 * t444 + (t446 * t716 + t448 * t715) * t534;
t258 = rSges(4,3) * t633 - t590;
t256 = Icges(4,1) * t376 + Icges(4,4) * t375 + Icges(4,5) * t633;
t255 = Icges(4,1) * t374 - Icges(4,4) * t373 + Icges(4,5) * t632;
t254 = Icges(4,4) * t376 + Icges(4,2) * t375 + Icges(4,6) * t633;
t253 = Icges(4,4) * t374 - Icges(4,2) * t373 + Icges(4,6) * t632;
t244 = t493 * t264;
t242 = -t361 * t702 + t534 * t611 - t623;
t241 = t362 * t702 + t673 * t697 + t442;
t240 = t464 * t263;
t238 = t435 * t697 + t436 * t464 + t437 * t465;
t237 = -t435 * t695 + t462 * t436 + t463 * t437;
t236 = t470 * t488 + t471 * t489 - t556 * t497 - t555 * t498 + (t487 * t670 - t496 * t541) * t534;
t235 = t468 * t488 + t469 * t489 + t506 * t497 - t557 * t498 + (t487 * t669 + t496 * t538) * t534;
t234 = (-t528 - t523) * qJD(1) + t590 - t594;
t233 = t545 + t257;
t232 = t610 - t624 + t318;
t231 = -t317 - t390 + t578;
t221 = t570 + t295;
t220 = t462 * t706 + t564 + t587;
t219 = t318 * t493 + t389 * t464;
t218 = -t317 * t493 - t389 * t462;
t213 = t356 * t702 - t493 * t358 + t494 * t360;
t212 = t355 * t702 - t493 * t357 + t494 * t359;
t195 = t702 * t378 + (t715 * t382 + t716 * t380 + (-t447 * t715 + t449 * t716) * qJD(2)) * t534;
t194 = t702 * t379 + (t715 * t383 + t716 * t381 + (-t446 * t715 + t448 * t716) * qJD(2)) * t534;
t189 = -t317 * t464 + t318 * t462;
t176 = t702 * t257 + t394 + (qJD(1) * t611 + t538 * t676) * t534;
t175 = -t702 * t258 + (t439 * t670 + t541 * t676) * t534 + t582;
t174 = -t375 * rSges(6,3) - t588;
t162 = -t317 * t702 + t534 * t599 + t569;
t161 = t318 * t702 + t644 * t697 + t678;
t153 = t570 + t686;
t152 = -t561 * t717 + t564 + t585 - t698 - t712;
t151 = -t229 * t558 + t287 * t562;
t150 = t228 * t558 - t287 * t561;
t147 = t375 * t436 + t376 * t437 + t462 * t430 + t463 * t431 + (-t429 * t541 + t435 * t670) * t534;
t146 = -t373 * t436 + t374 * t437 + t464 * t430 + t465 * t431 + (t429 * t538 + t435 * t669) * t534;
t143 = (t317 * t538 + t318 * t541) * t534 + t605;
t134 = (-qJD(1) * t532 - t509) * t541 - t185 - t267 + t671;
t133 = t266 + t545 + t184;
t131 = t295 * t493 + t464 * t679 + t244;
t130 = -t354 * t462 + t493 * t683 - t329;
t129 = -t228 * t562 + t229 * t561;
t124 = -t294 * t702 + t534 * t577 + t548;
t123 = t295 * t702 + t606 * t697 + t650;
t122 = t375 * t706 + t543 + t588;
t121 = t542 + t173;
t120 = t251 * t702 - t493 * t253 + t494 * t255 - t482 * t358 - t483 * t360;
t119 = t252 * t702 - t493 * t254 + t494 * t256 - t482 * t357 - t483 * t359;
t118 = (t257 * t541 + t258 * t538 + (t361 * t541 + (-t362 - t455) * t538) * qJD(1)) * t534 + t643;
t117 = -t294 * t464 + t462 * t682 - t240;
t105 = -t185 * t493 - t317 * t482 - t324 * t462 - t375 * t389;
t104 = t184 * t493 + t318 * t482 + t324 * t464 - t373 * t389;
t101 = (t294 * t538 + t295 * t541) * t534 + t576;
t97 = t702 * t184 + (qJD(1) * t599 + t538 * t646) * t534 + t681;
t96 = -t702 * t185 + (t389 * t670 + t541 * t646) * t534 + t553;
t95 = t464 * t648 + t493 * t686 + t244;
t94 = -t462 * t680 + t493 * t652 - t329;
t89 = t534 * t572 - t687 * t702 + t548;
t88 = t580 * t697 + t686 * t702 + t650;
t83 = t184 * t462 - t185 * t464 + t317 * t373 + t318 * t375;
t82 = t462 * t651 - t464 * t687 - t240;
t81 = (t538 * t687 + t541 * t686) * t534 + t576;
t80 = t281 * t717 + t543 + t586 - t699 - t713;
t79 = t542 + t691;
t70 = t702 * t173 + (qJD(1) * t577 + t538 * t607) * t534 + t653;
t69 = -t702 * t174 + (t354 * t670 + t541 * t607) * t534 + t547;
t66 = (t184 * t541 + t185 * t538 + (t317 * t541 + (-t318 + t677) * t538) * qJD(1)) * t534 + t575;
t63 = t116 * t558 - t193 * t561 - t228 * t404 + t281 * t287;
t62 = -t115 * t558 + t193 * t562 + t229 * t404 - t279 * t287;
t61 = -t303 * t462 - t354 * t375 + (-t160 - t174) * t493 + t683 * t482 + t684;
t60 = t173 * t493 + t295 * t482 + (t303 + t304) * t464 - t679 * t373 + t689;
t47 = t132 * t702 + (-t102 * t541 + t103 * t538) * t534;
t46 = -t102 * t462 - t103 * t464 + t132 * t493;
t45 = -t102 * t561 - t103 * t562 - t132 * t558;
t44 = t115 * t561 - t116 * t562 + t228 * t279 - t229 * t281;
t39 = -t174 * t464 + t294 * t373 + (t159 + t173) * t462 + t682 * t375 + t685;
t36 = -t126 * t558 - t561 * t92 - t562 * t93;
t35 = -t125 * t558 - t561 * t90 - t562 * t91;
t34 = (qJD(1) * t572 + t538 * t581) * t534 + t653 + t691 * t702;
t33 = (t541 * t581 + t670 * t680) * t534 + t547 - t690 * t702;
t32 = (t173 * t541 + t174 * t538 + (t294 * t541 + (-t295 + t649) * t538) * qJD(1)) * t534 + t563;
t31 = t688 * t462 - t680 * t375 + (-t160 - t690) * t493 + t652 * t482 + t684;
t30 = t691 * t493 + t686 * t482 + (t304 - t688) * t464 - t648 * t373 + t689;
t23 = -t690 * t464 + t687 * t373 + (t159 + t691) * t462 + t651 * t375 + t685;
t22 = (t691 * t541 + t690 * t538 + (t687 * t541 + (t649 - t686) * t538) * qJD(1)) * t534 + t563;
t7 = t102 * t281 + t103 * t279 - t28 * t561 - t29 * t562 + t705;
t2 = t125 * t404 - t26 * t561 - t27 * t562 + t279 * t91 + t281 * t90 - t43 * t558;
t1 = t126 * t404 - t24 * t561 - t25 * t562 + t279 * t93 + t281 * t92 - t42 * t558;
t10 = [t59 + (t152 * t80 + t153 * t79) * t613 + (t233 * t333 + t234 * t332) * t619 + (t347 * t421 + t348 * t420) * t621 + (t121 * t221 + t122 * t220) * t615 + (t133 * t232 + t134 * t231) * t617 + t759 + t775; m(3) * (t327 * t420 + t328 * t421 + t347 * t395 + t348 * t396) + (t175 * t332 + t176 * t333 + t233 * t241 + t234 * t242) * m(4) + (t133 * t161 + t134 * t162 + t231 * t96 + t232 * t97) * m(5) + (t121 * t123 + t122 * t124 + t220 * t69 + t221 * t70) * m(6) + (t152 * t33 + t153 * t34 + t79 * t88 + t80 * t89) * m(7) + t58 + ((-t236 / 0.2e1 - t147 / 0.2e1 - t119 / 0.2e1 - t194 / 0.2e1 + t566) * t541 + (t235 / 0.2e1 + t146 / 0.2e1 + t195 / 0.2e1 + t120 / 0.2e1 - t565) * t538 + ((t238 / 0.2e1 + t346 / 0.2e1 + t320 / 0.2e1 + t213 / 0.2e1 + t560) * t541 + (t237 / 0.2e1 + t345 / 0.2e1 + t319 / 0.2e1 + t212 / 0.2e1 - t559) * t538) * qJD(1)) * t534 + t748 + t758; (t242 * t175 + t241 * t176 + ((t361 * t538 + t362 * t541) * t534 + t675) * t118) * t619 + (t22 * t81 + t33 * t89 + t34 * t88) * t613 + (t101 * t32 + t123 * t70 + t124 * t69) * t615 + (t143 * t66 + t161 * t97 + t162 * t96) * t617 + (t396 * t327 + t395 * t328 + (t450 * t538 + t451 * t541) * (t384 * t541 + t385 * t538 + (t450 * t541 - t451 * t538) * qJD(1)) * t534 ^ 2) * t621 + (((-t119 - t194) * t541 + (t120 + t195) * t538 + ((t213 + t320) * t541 + (t212 + t319) * t538) * qJD(1)) * t534 + t738 + t748) * t702 + ((t235 + t146) * t702 + (t746 * t670 + t749 * t669 + (-t464 * t254 - t465 * t256 + t373 * t357 - t374 * t359 - t506 * t381 + t383 * t557 - t468 * t446 - t469 * t448 - t632 * t745) * t541 + (t464 * t253 + t465 * t255 - t373 * t358 + t374 * t360 + t506 * t380 - t557 * t382 + t468 * t447 + t469 * t449 + (t538 * t756 + t669 * t757 + t733) * t534) * t538) * t534 + t741) * t697 + ((-t236 - t147) * t702 + (t750 * t670 + t747 * t669 + (t462 * t254 + t463 * t256 + t375 * t357 + t376 * t359 - t556 * t381 - t555 * t383 + t470 * t446 + t471 * t448 + (t670 * t745 + t733) * t534) * t541 + (t556 * t380 + t555 * t382 - t470 * t447 - t471 * t449 - t462 * t253 - t463 * t255 - t375 * t358 - t376 * t360 + (t541 * t756 - t670 * t757) * t534) * t538) * t534 - t740) * t695 + ((t345 + t237) * t702 + (-t538 * t747 + t541 * t750) * t534 + t736) * t633 + ((t346 + t238) * t702 + (t538 * t749 - t541 * t746) * t534 + t735) * t632; ((t152 * t669 + t153 * t670 + t538 * t80 - t541 * t79) * m(7) + (-t121 * t541 + t122 * t538 + t220 * t669 + t221 * t670) * m(6) + (-t133 * t541 + t134 * t538 + t231 * t669 + t232 * t670) * m(5) + (-t233 * t541 + t234 * t538 + t332 * t669 + t333 * t670) * m(4)) * t534; m(4) * t702 * t118 + t66 * t640 + t32 * t639 + t22 * t638 + ((t33 * t538 - t34 * t541 + t669 * t89 + t670 * t88) * m(7) + (t123 * t670 + t124 * t669 + t538 * t69 - t541 * t70) * m(6) + (t161 * t670 + t162 * t669 + t538 * t96 - t541 * t97) * m(5) + (t175 * t538 - t176 * t541 + t241 * t670 + t242 * t669) * m(4)) * t534; 0; (t104 * t232 + t105 * t231 + t133 * t219 + t134 * t218) * m(5) + (t121 * t131 + t122 * t130 + t220 * t61 + t221 * t60) * m(6) + (t152 * t31 + t153 * t30 + t79 * t95 + t80 * t94) * m(7) + t565 * t464 + t566 * t462 + t559 * t375 + t560 * t373 + t704 + t755; (t22 * t82 + t23 * t81 + t30 * t88 + t31 * t89 + t33 * t94 + t34 * t95) * m(7) + (t101 * t39 + t117 * t32 + t123 * t60 + t124 * t61 + t130 * t69 + t131 * t70) * m(6) + (t104 * t161 + t105 * t162 + t143 * t83 + t189 * t66 + t218 * t96 + t219 * t97) * m(5) + t735 * t728 + t736 * t727 + t740 * t723 + t741 * t722 + (t47 + t751 * t702 + (t538 * t752 - t541 * t762) * t534) * t720 + t738 * t719 + t739 * t654 + t743 * t697 / 0.2e1 - t742 * t695 / 0.2e1 + (t538 * t734 + t541 * t737) * qJD(1) * t534 / 0.2e1; t83 * t640 + t39 * t639 + t23 * t638 + ((-t104 * t541 + t105 * t538 + t218 * t669 + t219 * t670) * m(5) + (t130 * t669 + t131 * t670 + t538 * t61 - t541 * t60) * m(6) + (-t30 * t541 + t31 * t538 + t669 * t94 + t670 * t95) * m(7)) * t534; t482 * t46 + (t23 * t82 + t30 * t95 + t31 * t94) * t613 + (t117 * t39 + t130 * t61 + t131 * t60) * t615 + (t104 * t219 + t105 * t218 + t189 * t83) * t617 + (t739 + t776) * t493 + (-t482 * t752 - t743) * t464 + (-t482 * t762 - t742) * t462 - t734 * t375 + t737 * t373; (-m(6) * t122 - m(7) * t80) * t464 + (-m(6) * t121 - m(7) * t79) * t462 + (-m(6) * t221 - m(7) * t153) * t375 + (m(6) * t220 + m(7) * t152) * t373; (m(6) * t32 + m(7) * t22) * t493 + (m(6) * t101 + m(7) * t81) * t482 + (-m(6) * t69 - m(7) * t33) * t464 + (-m(6) * t70 - m(7) * t34) * t462 + (-m(6) * t123 - m(7) * t88) * t375 + (m(6) * t124 + m(7) * t89) * t373; t663 * (t482 * t664 + (t373 * t732 + 0.2e1 * t375 * t541 + (-t462 * t538 - t464 * t541) * t662) * t534); (m(6) * t39 + m(7) * t23) * t493 + (m(6) * t117 + m(7) * t82) * t482 + (-m(6) * t61 - m(7) * t31) * t464 + (-m(6) * t60 - m(7) * t30) * t462 + (-m(6) * t131 - m(7) * t95) * t375 + (m(6) * t130 + m(7) * t94) * t373; 0.4e1 * t663 * (-t373 * t464 + t375 * t462 + t482 * t493); (t150 * t80 + t151 * t79 + t152 * t63 + t153 * t62) * m(7) - t658 * t562 - t659 * t561 + t629 * t281 + t628 * t279 + t705; (t129 * t22 + t150 * t33 + t151 * t34 + t44 * t81 + t62 * t88 + t63 * t89) * m(7) + t41 * t730 + t5 * t724 + t7 * t654 + t40 * t729 + t6 * t725 + t47 * t726 + t9 * t721 + (t1 * t718 - t541 * t2 / 0.2e1 + (t541 * t36 / 0.2e1 + t35 * t718) * qJD(1)) * t534; (t44 * t664 + (t63 * t732 - 0.2e1 * t541 * t62 + (t150 * t541 + t151 * t538) * t662) * t534) * t731; t45 * t720 + t7 * t719 + t36 * t728 + t1 * t722 + t38 * t730 + t3 * t724 + t35 * t727 + t2 * t723 + (t129 * t23 + t150 * t31 + t151 * t30 + t44 * t82 + t62 * t95 + t63 * t94) * m(7) + t46 * t726 + t8 * t721 + t37 * t729 + t4 * t725; (t129 * t482 + t150 * t373 - t151 * t375 + t44 * t493 - t462 * t62 - t464 * t63) * m(7); t279 * t36 - t562 * t1 + t281 * t35 - t561 * t2 + t404 * t45 - t558 * t7 + (t129 * t44 + t150 * t63 + t151 * t62) * t613;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
