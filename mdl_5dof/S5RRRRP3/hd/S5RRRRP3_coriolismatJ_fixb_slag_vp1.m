% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:08
% EndTime: 2019-12-31 21:49:18
% DurationCPUTime: 7.61s
% Computational Cost: add. (47835->423), mult. (31815->514), div. (0->0), fcn. (28123->8), ass. (0->276)
t422 = qJ(1) + qJ(2);
t421 = qJ(3) + t422;
t416 = cos(t421);
t411 = t416 * pkin(8);
t415 = sin(t421);
t423 = sin(qJ(4));
t425 = cos(qJ(4));
t541 = rSges(6,3) + qJ(5);
t587 = rSges(6,1) + pkin(4);
t633 = t541 * t423 + t587 * t425;
t625 = pkin(3) + t633;
t262 = t416 * rSges(6,2) - t625 * t415 + t411;
t417 = sin(t422);
t548 = pkin(2) * t417;
t255 = t262 - t548;
t550 = sin(qJ(1)) * pkin(1);
t239 = t255 - t550;
t458 = t587 * t423;
t660 = -t541 * t425 + t458;
t664 = t660 * t415;
t669 = t239 * t664;
t668 = t255 * t664;
t667 = t262 * t664;
t504 = t416 * t423;
t666 = t504 * t664;
t263 = (rSges(6,2) + pkin(8)) * t415 + t625 * t416;
t418 = cos(t422);
t547 = pkin(2) * t418;
t256 = t263 + t547;
t549 = cos(qJ(1)) * pkin(1);
t240 = t256 + t549;
t101 = -t263 * t239 + t240 * t262;
t542 = rSges(5,1) * t425;
t457 = pkin(3) + t542;
t508 = t415 * t423;
t468 = rSges(5,2) * t508 + t416 * rSges(5,3);
t281 = -t457 * t415 + t411 + t468;
t275 = t281 - t548;
t271 = t275 - t550;
t382 = rSges(5,2) * t504;
t282 = -t382 + t457 * t416 + (rSges(5,3) + pkin(8)) * t415;
t276 = t282 + t547;
t272 = t276 + t549;
t141 = -t282 * t271 + t272 * t281;
t354 = -rSges(4,1) * t415 - rSges(4,2) * t416;
t334 = t354 - t548;
t332 = t334 - t550;
t355 = t416 * rSges(4,1) - t415 * rSges(4,2);
t335 = t355 + t547;
t333 = t335 + t549;
t226 = -t355 * t332 + t333 * t354;
t480 = t355 * t334 - t335 * t354;
t488 = t282 * t275 - t276 * t281;
t493 = t263 * t255 - t256 * t262;
t615 = m(6) / 0.2e1;
t616 = m(5) / 0.2e1;
t617 = m(4) / 0.2e1;
t463 = (t480 + t226) * t617 + (t493 + t101) * t615 + (t488 + t141) * t616;
t464 = (-t480 + t226) * t617 + (-t493 + t101) * t615 + (-t488 + t141) * t616;
t4 = t464 - t463;
t665 = t4 * qJD(1);
t420 = Icges(5,4) * t425;
t392 = Icges(5,1) * t423 + t420;
t532 = Icges(6,5) * t425;
t663 = Icges(6,1) * t423 + t392 - t532;
t419 = Icges(6,5) * t423;
t447 = Icges(6,3) * t425 - t419;
t533 = Icges(5,4) * t423;
t661 = Icges(5,2) * t425 + t447 + t533;
t389 = -Icges(5,2) * t423 + t420;
t662 = t389 + t663;
t412 = t415 ^ 2;
t413 = t416 ^ 2;
t465 = t412 + t413;
t393 = Icges(5,1) * t425 - t533;
t632 = Icges(6,1) * t425 + t419;
t659 = t632 + t393;
t656 = (-Icges(5,6) + Icges(6,6)) * t425 + (-Icges(6,4) - Icges(5,5)) * t423;
t325 = Icges(6,4) * t415 + t416 * t632;
t327 = Icges(5,5) * t415 + t393 * t416;
t655 = -t661 * t416 + t325 + t327;
t324 = -Icges(6,4) * t416 + t415 * t632;
t380 = Icges(5,4) * t508;
t507 = t415 * t425;
t326 = Icges(5,1) * t507 - Icges(5,5) * t416 - t380;
t654 = -Icges(5,2) * t507 - t447 * t415 + t324 + t326 - t380;
t503 = t416 * t425;
t379 = Icges(6,5) * t503;
t317 = Icges(6,6) * t415 + Icges(6,3) * t504 + t379;
t323 = Icges(5,6) * t415 + t389 * t416;
t653 = -Icges(6,1) * t504 - t392 * t416 + t317 - t323 + t379;
t385 = Icges(6,3) * t423 + t532;
t316 = -Icges(6,6) * t416 + t385 * t415;
t322 = Icges(5,4) * t507 - Icges(5,2) * t508 - Icges(5,6) * t416;
t652 = t663 * t415 - t316 + t322;
t396 = rSges(5,1) * t423 + rSges(5,2) * t425;
t309 = t660 * t416;
t489 = t309 * t263 - t667;
t491 = -t309 * t256 + t668;
t538 = ((-t276 + t282) * t416 + (t275 - t281) * t415) * t396 * t616 + (t489 + t491) * t615;
t289 = -t416 * t458 + t541 * t503;
t134 = t256 * t289 + t668;
t138 = t263 * t289 + t667;
t348 = t396 * t415;
t350 = t396 * t416;
t162 = t275 * t348 - t276 * t350;
t173 = t281 * t348 - t282 * t350;
t638 = (t173 + t162) * t616 + (t138 + t134) * t615;
t651 = t538 - t638;
t257 = t271 * t348;
t492 = -t309 * t240 + t669;
t539 = (t257 + (-t281 * t415 + (-t272 + t282) * t416) * t396) * t616 + (t489 + t492) * t615;
t131 = t240 * t289 + t669;
t155 = -t272 * t350 + t257;
t639 = (t138 + t131) * t615 + (t173 + t155) * t616;
t650 = t539 - t639;
t586 = m(3) * (t549 * (-rSges(3,1) * t417 - rSges(3,2) * t418) + (t418 * rSges(3,1) - t417 * rSges(3,2)) * t550);
t460 = t256 * t504;
t153 = -t255 * t508 + t460;
t646 = t153 * m(6) * qJD(2);
t483 = t262 * t508 - t263 * t504;
t645 = t483 * m(6) * qJD(3);
t387 = Icges(6,4) * t425 + Icges(6,6) * t423;
t512 = t387 * t415;
t320 = -Icges(6,2) * t416 + t512;
t305 = t415 * t320;
t191 = t316 * t504 + t324 * t503 + t305;
t644 = t191 * t416;
t643 = (t659 - t661) * t425 + (t385 - t662) * t423;
t637 = (t316 * t423 + t324 * t425) * t415;
t96 = -t256 * t239 + t240 * t255;
t139 = -t276 * t271 + t272 * t275;
t214 = -t335 * t332 + t333 * t334;
t635 = t656 * t415;
t634 = t656 * t416;
t631 = -t655 * t423 + t653 * t425;
t630 = t654 * t423 + t652 * t425;
t629 = -m(5) * t173 - m(6) * t138;
t540 = (t257 + (-t275 * t415 + (-t272 + t276) * t416) * t396) * t616 + (-t491 + t492) * t615;
t626 = (t134 + t131) * t615 + (t162 + t155) * t616;
t440 = (-t385 / 0.2e1 + t662 / 0.2e1) * t425 + (-t661 / 0.2e1 + t659 / 0.2e1) * t423;
t623 = 0.4e1 * qJD(1);
t621 = 0.4e1 * qJD(2);
t620 = 0.2e1 * qJD(3);
t619 = 2 * qJD(4);
t229 = t240 * t504;
t599 = m(6) * (t229 - t460 + (-t239 + t255) * t508);
t598 = m(6) * (t229 + t460 + (-t239 - t255) * t508);
t147 = -t239 * t508 + t229;
t597 = m(6) * (t147 + t483);
t596 = m(6) * (t147 - t483);
t595 = m(6) * (t153 + t483);
t594 = m(6) * (t153 - t483);
t593 = m(6) * t96;
t591 = -t415 / 0.2e1;
t590 = t415 / 0.2e1;
t589 = -t416 / 0.2e1;
t582 = m(4) * t214;
t580 = m(4) * t226;
t579 = m(4) * t480;
t573 = m(5) * t139;
t571 = m(5) * t141;
t570 = m(5) * t488;
t568 = m(5) * t155;
t567 = m(5) * t162;
t565 = m(6) * t101;
t564 = m(6) * t493;
t482 = t289 * t508 + t666;
t487 = t239 * t503 + t240 * t507;
t563 = m(6) * (t482 + t487);
t486 = t255 * t503 + t256 * t507;
t562 = m(6) * (t482 + t486);
t454 = t309 * t508 - t666;
t560 = m(6) * (t454 + t487);
t484 = t262 * t503 + t263 * t507;
t559 = m(6) * (t482 + t484);
t558 = m(6) * (t454 + t486);
t557 = m(6) * t131;
t556 = m(6) * (t454 + t484);
t555 = m(6) * t134;
t544 = m(6) * qJD(4);
t543 = m(6) * qJD(5);
t519 = t322 * t423;
t386 = Icges(5,5) * t425 - Icges(5,6) * t423;
t513 = t386 * t416;
t511 = t387 * t416;
t505 = t416 * t320;
t498 = t423 * t425;
t485 = (-t289 - t309) * t664;
t481 = -t309 * t503 - t507 * t664;
t318 = Icges(5,5) * t507 - Icges(5,6) * t508 - Icges(5,3) * t416;
t479 = -t415 * t318 - t326 * t503;
t319 = Icges(5,3) * t415 + t513;
t478 = t415 * t319 + t327 * t503;
t469 = t465 * t498;
t461 = m(6) * t147 * qJD(1);
t321 = Icges(6,2) * t415 + t511;
t192 = t317 * t504 + t415 * t321 + t325 * t503;
t296 = t327 * t507;
t456 = t416 * t319 - t296;
t455 = t323 * t423 - t318;
t453 = -t317 * t508 + t321 * t416 - t325 * t507;
t442 = t192 + t505;
t441 = (-t348 * t416 + t350 * t415) * t396;
t432 = t440 + t626;
t429 = -t440 + ((t316 + t322) * t425 + (-t324 + t326) * t423) * (t590 + t591);
t187 = -t505 + t637;
t126 = -t187 * t416 - t415 * t453;
t190 = -t323 * t508 - t456;
t127 = -(-(-t326 * t425 + t519) * t415 - t416 * t318) * t416 + t190 * t415;
t128 = t192 * t415 - t644;
t193 = -t322 * t504 - t479;
t194 = -t323 * t504 + t478;
t129 = -t193 * t416 + t194 * t415;
t37 = (t192 - t442) * t416 + (t191 + t453 - t305) * t415;
t38 = (t455 * t416 + t194 - t478) * t416 + (t455 * t415 + t193 + t456) * t415;
t39 = -t644 + (t187 + t442 - t637) * t415;
t40 = (t190 - t296 + (t319 + t519) * t416 + t479) * t416 + t478 * t415;
t427 = ((t39 + t40) * t416 / 0.2e1 + (t126 + t127 + t37 + t38) * t591 + (t386 * t415 + t643 * t416 + t653 * t423 + t655 * t425 + t512) * t590 + (t643 * t415 - t652 * t423 + t654 * t425 + t128 + t129 - t511 - t513) * t589) * qJD(4);
t400 = -rSges(5,2) * t423 + t542;
t351 = t465 * t423;
t315 = t469 - t498;
t310 = t633 * t416;
t308 = t633 * t415;
t197 = t289 * t416 - t412 * t660;
t176 = t465 * t633;
t133 = t176 * t351 + t481;
t130 = t556 / 0.2e1;
t124 = t558 / 0.2e1;
t116 = t559 / 0.2e1;
t114 = t560 / 0.2e1;
t108 = t562 / 0.2e1;
t104 = t563 / 0.2e1;
t91 = t594 / 0.2e1;
t90 = t595 / 0.2e1;
t87 = t596 / 0.2e1;
t86 = t597 / 0.2e1;
t83 = t598 / 0.2e1;
t82 = t599 / 0.2e1;
t56 = t440 - t629;
t55 = t440 + t555 + t567;
t54 = t440 + t557 + t568;
t47 = -t564 - t570 - t579;
t44 = t565 + t571 + t580;
t41 = t573 + t582 + t586 + t593;
t36 = t130 - t559 / 0.2e1;
t35 = t130 + t116;
t34 = t116 - t556 / 0.2e1;
t29 = t124 - t562 / 0.2e1;
t28 = t124 + t108;
t27 = t108 - t558 / 0.2e1;
t26 = t114 - t563 / 0.2e1;
t25 = t114 + t104;
t24 = t104 - t560 / 0.2e1;
t23 = t91 - t595 / 0.2e1;
t22 = t91 + t90;
t21 = t90 - t594 / 0.2e1;
t20 = t87 - t597 / 0.2e1;
t19 = t87 + t86;
t18 = t86 - t596 / 0.2e1;
t17 = t83 - t599 / 0.2e1;
t16 = t83 + t82;
t15 = t82 - t598 / 0.2e1;
t14 = t440 - t651;
t13 = t440 + t538 + t638;
t12 = t440 + t539 + t639;
t11 = t440 - t650;
t10 = t432 + t540;
t9 = t432 - t540;
t8 = t429 + t651;
t7 = t429 + t650;
t6 = t429 + t540 - t626;
t3 = t463 + t464;
t2 = (t129 / 0.2e1 - t40 / 0.2e1 + t128 / 0.2e1 - t39 / 0.2e1) * t416 + (t38 / 0.2e1 + t127 / 0.2e1 + t37 / 0.2e1 + t126 / 0.2e1) * t415;
t1 = t2 * qJD(4);
t5 = [t41 * qJD(2) + t44 * qJD(3) + t54 * qJD(4) + t147 * t543, t41 * qJD(1) + t3 * qJD(3) + t10 * qJD(4) + t16 * qJD(5) + 0.2e1 * (t586 / 0.2e1 + t139 * t616 + t214 * t617 + t96 * t615) * qJD(2), t44 * qJD(1) + t3 * qJD(2) + t12 * qJD(4) + t19 * qJD(5) + (t101 * t615 + t141 * t616 + t226 * t617) * t620, t54 * qJD(1) + t10 * qJD(2) + t12 * qJD(3) + t427 + t25 * qJD(5) + ((-t239 * t310 - t240 * t308 + t485) * t615 + ((-t271 * t416 - t272 * t415) * t400 + t441) * t616) * t619, t16 * qJD(2) + t19 * qJD(3) + t25 * qJD(4) + t461; t4 * qJD(3) + t9 * qJD(4) + t17 * qJD(5) + (-t593 / 0.4e1 - t573 / 0.4e1 - t582 / 0.4e1 - t586 / 0.4e1) * t623, t47 * qJD(3) + t55 * qJD(4) + t153 * t543, t665 + t47 * qJD(2) + t13 * qJD(4) + t22 * qJD(5) + (-t480 * t617 - t488 * t616 - t493 * t615) * t620, t9 * qJD(1) + t55 * qJD(2) + t13 * qJD(3) + t427 + t28 * qJD(5) + ((-t255 * t310 - t256 * t308 + t485) * t615 + ((-t275 * t416 - t276 * t415) * t400 + t441) * t616) * t619, t17 * qJD(1) + t22 * qJD(3) + t28 * qJD(4) + t646; -t4 * qJD(2) + t11 * qJD(4) + t20 * qJD(5) + (-t565 / 0.4e1 - t571 / 0.4e1 - t580 / 0.4e1) * t623, -t665 + t14 * qJD(4) + t23 * qJD(5) + (t564 / 0.4e1 + t570 / 0.4e1 + t579 / 0.4e1) * t621, t56 * qJD(4) - t483 * t543, t11 * qJD(1) + t14 * qJD(2) + t56 * qJD(3) + t427 + t35 * qJD(5) + ((-t262 * t310 - t263 * t308 + t485) * t615 + ((-t281 * t416 - t282 * t415) * t400 + t441) * t616) * t619, t20 * qJD(1) + t23 * qJD(2) + t35 * qJD(4) - t645; t429 * qJD(1) + t6 * qJD(2) + t7 * qJD(3) + t1 + t26 * qJD(5) + (-t557 / 0.4e1 - t568 / 0.4e1) * t623, t6 * qJD(1) + t429 * qJD(2) + t8 * qJD(3) + t1 + t29 * qJD(5) + (-t555 / 0.4e1 - t567 / 0.4e1) * t621, t7 * qJD(1) + t8 * qJD(2) + t36 * qJD(5) + t1 + (t429 + t629) * qJD(3), (m(5) * ((t415 * (rSges(5,1) * t507 - t468) + t416 * (rSges(5,1) * t503 + t415 * rSges(5,3) - t382)) * (-t348 * t415 - t350 * t416) + t465 * t400 * t396) + m(6) * (t176 * t197 + t308 * t664 + t309 * t310) + (t634 * t412 + (t630 * t416 + (t631 - t635) * t415) * t416) * t590 + (t635 * t413 + (t631 * t415 + (t630 - t634) * t416) * t415) * t589) * qJD(4) + t133 * t543 + (qJD(1) + qJD(2) + qJD(3)) * t2, t26 * qJD(1) + t29 * qJD(2) + t36 * qJD(3) + t133 * t544 + (-t351 * t425 - t315 + t469) * t543; t15 * qJD(2) + t18 * qJD(3) + t24 * qJD(4) - t461, t15 * qJD(1) + t21 * qJD(3) + t27 * qJD(4) - t646, t18 * qJD(1) + t21 * qJD(2) + t34 * qJD(4) + t645, t24 * qJD(1) + t27 * qJD(2) + t34 * qJD(3) + (-t197 * t425 + (-t308 * t415 - t310 * t416 + t176) * t423 - t133 + t481) * t544 + t315 * t543, t315 * t544;];
Cq = t5;
