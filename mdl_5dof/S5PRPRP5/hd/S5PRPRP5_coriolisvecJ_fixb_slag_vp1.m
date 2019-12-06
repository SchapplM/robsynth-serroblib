% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:38:23
% DurationCPUTime: 35.58s
% Computational Cost: add. (15675->685), mult. (24801->998), div. (0->0), fcn. (23734->8), ass. (0->352)
t623 = Icges(6,4) + Icges(5,5);
t622 = Icges(5,6) - Icges(6,6);
t310 = pkin(8) + qJ(4);
t305 = sin(t310);
t306 = cos(t310);
t314 = cos(pkin(7));
t312 = sin(pkin(7));
t317 = cos(qJ(2));
t497 = t312 * t317;
t236 = t305 * t497 + t306 * t314;
t316 = sin(qJ(2));
t498 = t312 * t316;
t237 = -t314 * t305 + t306 * t497;
t519 = Icges(5,4) * t237;
t102 = -Icges(5,2) * t236 + Icges(5,6) * t498 + t519;
t228 = Icges(6,5) * t237;
t96 = Icges(6,6) * t498 + Icges(6,3) * t236 + t228;
t637 = t102 - t96;
t495 = t314 * t317;
t238 = t305 * t495 - t312 * t306;
t496 = t314 * t316;
t503 = t305 * t312;
t239 = t306 * t495 + t503;
t518 = Icges(5,4) * t239;
t103 = -Icges(5,2) * t238 + Icges(5,6) * t496 + t518;
t229 = Icges(6,5) * t239;
t97 = Icges(6,6) * t496 + Icges(6,3) * t238 + t229;
t636 = t103 - t97;
t100 = Icges(6,4) * t237 + Icges(6,2) * t498 + Icges(6,6) * t236;
t98 = Icges(5,5) * t237 - Icges(5,6) * t236 + Icges(5,3) * t498;
t612 = t98 + t100;
t101 = Icges(6,4) * t239 + Icges(6,2) * t496 + Icges(6,6) * t238;
t99 = Icges(5,5) * t239 - Icges(5,6) * t238 + Icges(5,3) * t496;
t611 = t99 + t101;
t513 = Icges(6,5) * t236;
t104 = Icges(6,1) * t237 + Icges(6,4) * t498 + t513;
t230 = Icges(5,4) * t236;
t106 = Icges(5,1) * t237 + Icges(5,5) * t498 - t230;
t635 = t104 + t106;
t512 = Icges(6,5) * t238;
t105 = Icges(6,1) * t239 + Icges(6,4) * t496 + t512;
t231 = Icges(5,4) * t238;
t107 = Icges(5,1) * t239 + Icges(5,5) * t496 - t231;
t634 = t105 + t107;
t599 = Icges(5,2) + Icges(6,3);
t656 = Icges(6,2) + Icges(5,3);
t501 = t306 * t316;
t294 = Icges(6,5) * t501;
t502 = t305 * t316;
t507 = Icges(6,6) * t317;
t211 = Icges(6,3) * t502 + t294 - t507;
t516 = Icges(5,4) * t306;
t401 = -Icges(5,2) * t305 + t516;
t217 = -Icges(5,6) * t317 + t316 * t401;
t631 = t211 - t217;
t398 = Icges(5,5) * t306 - Icges(5,6) * t305;
t213 = -Icges(5,3) * t317 + t316 * t398;
t400 = Icges(6,4) * t306 + Icges(6,6) * t305;
t215 = -Icges(6,2) * t317 + t316 * t400;
t655 = t213 + t215;
t511 = Icges(6,5) * t305;
t404 = Icges(6,1) * t306 + t511;
t219 = -Icges(6,4) * t317 + t316 * t404;
t517 = Icges(5,4) * t305;
t405 = Icges(5,1) * t306 - t517;
t221 = -Icges(5,5) * t317 + t316 * t405;
t597 = t219 + t221;
t654 = (-t623 * t305 - t622 * t306) * t316;
t620 = -t636 * t236 + t634 * t237 + t611 * t498;
t619 = -t637 * t238 + t635 * t239 + t612 * t496;
t621 = -t637 * t236 + t635 * t237 + t612 * t498;
t583 = -t636 * t238 + t634 * t239 + t611 * t496;
t653 = t631 * t236 + t597 * t237 + t655 * t498;
t652 = t631 * t238 + t597 * t239 + t655 * t496;
t467 = qJD(4) * t317;
t451 = t306 * t467;
t471 = qJD(2) * t316;
t456 = t312 * t471;
t166 = -t312 * t451 + (qJD(4) * t314 + t456) * t305;
t167 = -qJD(4) * t236 - t306 * t456;
t470 = qJD(2) * t317;
t455 = t312 * t470;
t651 = t622 * t166 + t623 * t167 + t656 * t455;
t454 = t314 * t471;
t168 = -qJD(4) * t503 + t305 * t454 - t314 * t451;
t169 = -qJD(4) * t238 - t306 * t454;
t453 = t314 * t470;
t650 = t622 * t168 + t623 * t169 + t656 * t453;
t649 = Icges(5,1) + Icges(6,1);
t648 = Icges(5,4) - Icges(6,5);
t214 = Icges(5,3) * t316 + t317 * t398;
t216 = Icges(6,2) * t316 + t317 * t400;
t647 = t654 * qJD(4) + (t214 + t216) * qJD(2);
t510 = Icges(6,5) * t306;
t397 = Icges(6,3) * t305 + t510;
t212 = Icges(6,6) * t316 + t317 * t397;
t218 = Icges(5,6) * t316 + t317 * t401;
t576 = -t212 + t218;
t220 = Icges(6,4) * t316 + t317 * t404;
t222 = Icges(5,5) * t316 + t317 * t405;
t575 = -t220 - t222;
t645 = (-t599 * t306 + t511 - t517) * t316;
t644 = t620 * t314;
t643 = t619 * t312;
t308 = t312 ^ 2;
t309 = t314 ^ 2;
t474 = t308 + t309;
t642 = m(3) * qJD(2) ^ 2 * (rSges(3,1) * t316 + rSges(3,2) * t317) * t474;
t641 = -t599 * t166 - t648 * t167 - t622 * t455;
t640 = -t599 * t168 - t648 * t169 - t622 * t453;
t639 = t648 * t166 + t649 * t167 + t623 * t455;
t638 = t648 * t168 + t649 * t169 + t623 * t453;
t633 = t576 * qJD(2) + t645 * qJD(4);
t253 = (-Icges(5,1) * t305 - t516) * t316;
t468 = qJD(4) * t316;
t632 = -(-Icges(6,1) * t305 + t510) * t468 - qJD(4) * t253 + t575 * qJD(2);
t630 = -t647 * t316 - t470 * t655;
t629 = t650 * t316 + t611 * t470;
t628 = t651 * t316 + t612 * t470;
t627 = t583 * t314 + t643;
t626 = t621 * t312 + t644;
t625 = t652 * t316;
t624 = t653 * t316;
t473 = qJD(2) * t312;
t283 = t314 * t468 + t473;
t544 = t283 / 0.2e1;
t472 = qJD(2) * t314;
t284 = t312 * t468 - t472;
t542 = t284 / 0.2e1;
t589 = t637 * t166 + t635 * t167 + t641 * t236 + t639 * t237 + t628 * t312;
t588 = t636 * t166 + t634 * t167 + t640 * t236 + t638 * t237 + t629 * t312;
t587 = t637 * t168 + t635 * t169 + t641 * t238 + t639 * t239 + t628 * t314;
t586 = t636 * t168 + t634 * t169 + t640 * t238 + t638 * t239 + t629 * t314;
t409 = t104 * t306 + t305 * t96;
t53 = -t100 * t317 + t316 * t409;
t395 = -t102 * t305 + t106 * t306;
t55 = t316 * t395 - t317 * t98;
t618 = t53 + t55;
t408 = t105 * t306 + t305 * t97;
t54 = -t101 * t317 + t316 * t408;
t394 = -t103 * t305 + t107 * t306;
t56 = t316 * t394 - t317 * t99;
t617 = t54 + t56;
t390 = t211 * t305 + t219 * t306;
t504 = t215 * t317;
t73 = t316 * t390 - t504;
t389 = -t217 * t305 + t221 * t306;
t505 = t213 * t317;
t74 = t316 * t389 - t505;
t602 = t73 + t74;
t421 = rSges(5,1) * t306 - rSges(5,2) * t305;
t225 = -rSges(5,3) * t317 + t316 * t421;
t610 = t225 * t283;
t609 = t225 * t284;
t606 = (t631 * t168 - t597 * t169 + t633 * t238 + t632 * t239 + t630 * t314) * t317 + (t627 * t317 + t625) * qJD(2);
t605 = (t631 * t166 - t597 * t167 + t633 * t236 + t632 * t237 + t630 * t312) * t317 + (t626 * t317 + t624) * qJD(2);
t311 = sin(pkin(8));
t313 = cos(pkin(8));
t520 = Icges(3,4) * t317;
t521 = Icges(3,4) * t316;
t604 = -(-Icges(3,1) * t316 - t520) * t470 + (-(-Icges(4,5) * t313 + Icges(4,6) * t311) * t316 - t521 + (-Icges(4,3) - Icges(3,2)) * t317) * t471;
t603 = t587 * t542 + t586 * t544 + t606 * qJD(4) / 0.2e1;
t596 = t654 * t467 + (t623 * t236 + t622 * t237) * t284 + (t623 * t238 + t622 * t239) * t283;
t412 = t55 * t312 + t56 * t314;
t413 = t53 * t312 + t54 * t314;
t595 = t412 + t413;
t307 = qJD(3) * t316;
t299 = t312 * t307;
t290 = pkin(2) * t316 - qJ(3) * t317;
t367 = qJD(2) * t290;
t205 = -t312 * t367 + t299;
t301 = t314 * t307;
t206 = -t314 * t367 + t301;
t594 = t312 * t205 + t314 * t206 + t367 * t474 - t307;
t593 = (t283 * t620 + t284 * t621 - t467 * t653) * t312 + (t283 * t583 + t284 * t619 - t467 * t652) * t314;
t554 = -rSges(4,1) * t313 + rSges(4,2) * t311;
t552 = rSges(4,3) * t317 + t316 * t554;
t592 = qJD(4) * t605 + t283 * t588 + t284 * t589;
t590 = rSges(6,1) + pkin(4);
t563 = ((t395 + t409) * qJD(2) - t651) * t317 + (t639 * t306 + t641 * t305 + (-t305 * t635 - t306 * t637) * qJD(4) + t612 * qJD(2)) * t316;
t562 = ((t394 + t408) * qJD(2) - t650) * t317 + (t638 * t306 + t640 * t305 + (-t305 * t634 - t306 * t636) * qJD(4) + t611 * qJD(2)) * t316;
t585 = t283 * t617 + t284 * t618 - t467 * t602;
t582 = rSges(6,3) + qJ(5);
t346 = -t316 * t397 + t507;
t176 = t346 * t312;
t182 = t217 * t312;
t580 = t176 + t182;
t177 = t346 * t314;
t183 = t217 * t314;
t579 = t177 + t183;
t184 = t219 * t312;
t186 = t221 * t312;
t578 = -t184 - t186;
t185 = t219 * t314;
t187 = t221 * t314;
t577 = -t185 - t187;
t374 = t214 - t389;
t375 = -t216 + t390;
t549 = (t215 * t314 + t408) * t283 + (t215 * t312 + t409) * t284;
t550 = -(-t213 * t314 - t394) * t283 - (-t213 * t312 - t395) * t284;
t574 = (-t549 - t550 + (-t374 + t375) * t467) * t316;
t272 = -t311 * t495 + t312 * t313;
t500 = t311 * t312;
t273 = t313 * t495 + t500;
t403 = -Icges(3,2) * t316 + t520;
t407 = Icges(3,1) * t317 - t521;
t573 = t604 * t314 + ((-Icges(4,5) * t273 + Icges(3,6) * t312 - Icges(4,6) * t272 - Icges(4,3) * t496 + t314 * t403) * t317 + (-(Icges(4,4) * t273 + Icges(4,2) * t272 + Icges(4,6) * t496) * t311 + (Icges(4,1) * t273 + Icges(4,4) * t272 + Icges(4,5) * t496) * t313 + Icges(3,5) * t312 + t314 * t407) * t316) * qJD(2);
t270 = -t311 * t497 - t313 * t314;
t499 = t311 * t314;
t271 = t313 * t497 - t499;
t572 = t604 * t312 + ((-Icges(4,5) * t271 - Icges(3,6) * t314 - Icges(4,6) * t270 - Icges(4,3) * t498 + t312 * t403) * t317 + (-(Icges(4,4) * t271 + Icges(4,2) * t270 + Icges(4,6) * t498) * t311 + (Icges(4,1) * t271 + Icges(4,4) * t270 + Icges(4,5) * t498) * t313 - Icges(3,5) * t314 + t312 * t407) * t316) * qJD(2);
t419 = pkin(4) * t306 + qJ(5) * t305;
t420 = rSges(6,1) * t306 + rSges(6,3) * t305;
t571 = rSges(6,2) * t317 + (-t419 - t420) * t316;
t570 = t283 * t611 + t284 * t612;
t538 = pkin(3) * t313;
t547 = pkin(6) * t317 - t316 * t538;
t353 = qJD(2) * t547;
t162 = t312 * t353;
t163 = t314 * t353;
t569 = t312 * t162 + t314 * t163 - t353 * t474 + t594;
t568 = -t504 - t505;
t564 = t474 * qJD(2) * t552;
t558 = t602 * t471 + (t647 * t317 + (t632 * t306 + t633 * t305 + (t305 * t597 - t306 * t631) * qJD(4)) * t316 + ((-t389 - t390) * t317 - t655 * t316 + t595) * qJD(2)) * t317;
t557 = (t597 + t645) * t467 + (t237 * t599 + t230 - t513 - t635) * t284 + (t239 * t599 + t231 - t512 - t634) * t283;
t556 = (Icges(6,1) * t502 - t253 - t294 - t631) * t467 + (-t236 * t649 + t228 - t519 - t637) * t284 + (-t238 * t649 + t229 - t518 - t636) * t283;
t555 = t596 * t316;
t464 = qJD(5) * t305;
t289 = t316 * t464;
t210 = pkin(6) * t316 + t317 * t538;
t160 = -pkin(3) * t499 + t210 * t312;
t161 = pkin(3) * t500 + t210 * t314;
t292 = pkin(2) * t317 + qJ(3) * t316;
t281 = t292 * t312;
t282 = t292 * t314;
t469 = qJD(3) * t317;
t384 = t281 * t473 + t282 * t472 + qJD(1) - t469;
t360 = t160 * t473 + t161 * t472 + t384;
t492 = rSges(6,2) * t496 + t238 * t582 + t239 * t590;
t493 = rSges(6,2) * t498 + t236 * t582 + t237 * t590;
t27 = t283 * t493 - t284 * t492 + t289 + t360;
t463 = qJD(2) * qJD(3);
t459 = t205 * t473 + t206 * t472 + t316 * t463;
t411 = t162 * t473 + t163 * t472 + t459;
t452 = t306 * t468;
t465 = qJD(5) * t238;
t536 = rSges(6,2) * t453 - t168 * t582 + t169 * t590 + t465;
t466 = qJD(5) * t236;
t537 = rSges(6,2) * t455 - t166 * t582 + t167 * t590 + t466;
t5 = qJD(5) * t452 - t536 * t284 + t537 * t283 + (t464 + (-t312 * t492 + t314 * t493) * qJD(4)) * t470 + t411;
t551 = t27 * t536 + t492 * t5;
t546 = -t474 * t471 / 0.2e1;
t545 = -t283 / 0.2e1;
t543 = -t284 / 0.2e1;
t539 = -t317 / 0.2e1;
t226 = rSges(6,2) * t316 + t317 * t420;
t255 = (-rSges(6,1) * t305 + rSges(6,3) * t306) * t316;
t362 = t305 * t470 + t452;
t522 = qJD(2) * t226 + qJD(4) * t255 + t289 + t362 * qJ(5) + (-t305 * t468 + t306 * t470) * pkin(4);
t491 = -t236 * t590 + t237 * t582;
t490 = t238 * t590 - t239 * t582;
t489 = t571 * t312;
t488 = t571 * t314;
t269 = qJD(2) * t292 - t469;
t487 = -t210 * qJD(2) - t269;
t485 = -t290 + t547;
t484 = -t210 - t292;
t247 = rSges(4,3) * t316 - t317 * t554;
t483 = -t247 * qJD(2) - t269;
t478 = -t290 + t552;
t477 = -t247 - t292;
t476 = (-pkin(4) * t305 + qJ(5) * t306) * t316 + t255;
t475 = t312 * t281 + t314 * t282;
t462 = qJD(2) * qJD(4);
t227 = rSges(5,3) * t316 + t317 * t421;
t256 = (-rSges(5,1) * t305 - rSges(5,2) * t306) * t316;
t122 = qJD(2) * t227 + qJD(4) * t256;
t460 = -t122 + t487;
t458 = -t225 + t485;
t448 = t317 * t463;
t442 = -t467 / 0.2e1;
t441 = t467 / 0.2e1;
t440 = t462 / 0.2e1;
t439 = t487 * t312;
t438 = t487 * t314;
t437 = t484 * t312;
t436 = t484 * t314;
t435 = qJD(2) * t485;
t434 = qJD(2) * t483;
t433 = qJD(2) * t478;
t432 = t487 - t522;
t431 = t312 * t160 + t314 * t161 + t475;
t429 = t571 + t485;
t424 = t316 * t440;
t423 = t317 * t440;
t399 = -Icges(3,5) * t316 - Icges(3,6) * t317;
t109 = rSges(5,1) * t237 - rSges(5,2) * t236 + rSges(5,3) * t498;
t111 = rSges(5,1) * t239 - rSges(5,2) * t238 + rSges(5,3) * t496;
t393 = t109 * t314 - t111 * t312;
t386 = t312 * t423;
t385 = t314 * t423;
t373 = t312 * t435 + t299;
t372 = t314 * t435 + t301;
t40 = t109 * t283 - t111 * t284 + t360;
t366 = t40 * t393;
t363 = qJD(2) * t399;
t361 = t27 * t537 + t493 * t5;
t43 = t283 * t571 - t467 * t492 + t373 + t466;
t44 = -t284 * t571 + t467 * t493 + t372 + t465;
t359 = t43 * t492 - t44 * t493;
t351 = Icges(4,5) * t317 + (-Icges(4,1) * t313 + Icges(4,4) * t311) * t316;
t348 = Icges(4,6) * t317 + (-Icges(4,4) * t313 + Icges(4,2) * t311) * t316;
t158 = rSges(4,1) * t271 + rSges(4,2) * t270 + rSges(4,3) * t498;
t159 = rSges(4,1) * t273 + rSges(4,2) * t272 + rSges(4,3) * t496;
t61 = (t158 * t312 + t159 * t314) * qJD(2) + t384;
t341 = t61 * t552;
t339 = qJD(2) * t351;
t338 = qJD(2) * t348;
t322 = (t27 * t493 + t43 * t571) * t314 + (-t27 * t492 - t44 * t571) * t312;
t302 = t314 * t469;
t300 = t312 * t469;
t296 = t314 * t448;
t295 = t312 * t448;
t276 = t399 * t314;
t275 = t399 * t312;
t262 = t314 * t363;
t261 = t312 * t363;
t204 = t351 * t314;
t203 = t351 * t312;
t202 = t348 * t314;
t201 = t348 * t312;
t175 = t314 * t339;
t174 = t312 * t339;
t173 = t314 * t338;
t172 = t312 * t338;
t156 = -rSges(5,1) * t238 - rSges(5,2) * t239;
t152 = -rSges(5,1) * t236 - rSges(5,2) * t237;
t126 = t314 * t433 + t301;
t125 = t312 * t433 + t299;
t92 = t314 * t434 + t296;
t91 = t312 * t434 + t295;
t90 = rSges(5,1) * t169 + rSges(5,2) * t168 + rSges(5,3) * t453;
t88 = rSges(5,1) * t167 + rSges(5,2) * t166 + rSges(5,3) * t455;
t62 = qJD(2) * t564 + t459;
t60 = t109 * t467 + t372 + t609;
t59 = -t111 * t467 + t373 - t610;
t35 = t88 * t467 + t122 * t284 + t296 + (t438 + (-t109 * t316 + t225 * t497) * qJD(4)) * qJD(2);
t34 = -t90 * t467 - t122 * t283 + t295 + (t439 + (t111 * t316 - t225 * t495) * qJD(4)) * qJD(2);
t26 = t317 * t393 * t462 + t283 * t88 - t284 * t90 + t411;
t7 = -qJD(5) * t168 + t296 + t522 * t284 + t537 * t467 + (t438 + (-t316 * t493 - t497 * t571) * qJD(4)) * qJD(2);
t6 = -qJD(5) * t166 + t295 - t522 * t283 - t536 * t467 + (t439 + (t316 * t492 + t495 * t571) * qJD(4)) * qJD(2);
t1 = [m(4) * t62 + m(5) * t26 + m(6) * t5 - t642; (((t238 * t576 + t239 * t575 + t643) * t317 + t625) * qJD(4) + (((t568 + t583) * qJD(4) + t570) * t317 + t574) * t314 + (t238 * t580 + t239 * t578) * t284 + (t238 * t579 + t239 * t577) * t283) * t545 + (t312 * t586 - t314 * t587) * t544 + (((t236 * t576 + t237 * t575 + t644) * t317 + t624) * qJD(4) + (((t568 + t621) * qJD(4) + t570) * t317 + t574) * t312 + (t236 * t580 + t237 * t578) * t284 + (t236 * t579 + t237 * t577) * t283) * t543 + (t312 * t588 - t314 * t589) * t542 + t312 * t603 - t592 * t314 / 0.2e1 + ((t173 * t272 + t175 * t273 + t262 * t312) * t312 + (-t172 * t272 - t174 * t273 + t572 * t314 + (-t261 - t573) * t312) * t314) * t473 + ((t172 * t270 + t174 * t271 - t261 * t314) * t314 + (-t173 * t270 - t175 * t271 + t573 * t312 + (t262 - t572) * t314) * t312) * t472 - (t276 * qJD(2) * t308 + (t202 * t272 + t204 * t273) * t473 + (-t272 * t201 - t273 * t203 - t275 * t312) * t472) * t473 / 0.2e1 + (t275 * qJD(2) * t309 - (t201 * t270 + t203 * t271) * t472 + (t270 * t202 + t271 * t204 - t276 * t314) * t473) * t472 / 0.2e1 - t585 * t468 / 0.2e1 + (((t177 * t305 - t185 * t306 + t101) * t283 + (t176 * t305 - t184 * t306 + t100) * t284 + t73 * qJD(4)) * t316 + ((-t375 * t317 + (-t212 * t305 - t220 * t306 - t215) * t316 + t413) * qJD(4) + t549) * t317 + ((t183 * t305 - t187 * t306 + t99) * t283 + (t182 * t305 - t186 * t306 + t98) * t284 + t74 * qJD(4)) * t316 + ((t374 * t317 + (t218 * t305 - t222 * t306 - t213) * t316 + t412) * qJD(4) + t550) * t317) * t441 + (t312 * t617 - t314 * t618) * t424 + (t312 * t620 - t314 * t621) * t386 + (t583 * t312 - t314 * t619) * t385 + (t5 * t431 + (t429 * t7 + t432 * t44 + t551) * t314 + (t429 * t6 + t43 * t432 + t361) * t312 - t44 * (-t289 * t314 + t302) - t43 * (-t289 * t312 + t300) - (t43 * t437 + t436 * t44) * qJD(2) - (t359 * t316 + (-t43 * t488 + t44 * t489 + t322) * t317) * qJD(4) + (t43 * t283 - t44 * t284) * (t419 * t317 + t226) + (-t489 * t283 + t488 * t284 - t317 * t464 + t569) * t27) * m(6) + (-t60 * (t227 * t284 + t302) - t59 * (-t227 * t283 + t300) - (t436 * t60 + t437 * t59) * qJD(2) - ((-t109 * t60 + t111 * t59) * t316 + t366 * t317) * qJD(4) + t26 * t431 + t569 * t40 + (t26 * t111 + t35 * t458 + t460 * t60 + (-t609 + t90) * t40) * t314 + (t26 * t109 + t34 * t458 + t460 * t59 + (t610 + t88) * t40) * t312) * m(5) + (t62 * t475 + (t126 * t483 + t62 * t159 + t478 * t92) * t314 + (t125 * t483 + t62 * t158 + t478 * t91) * t312 - t126 * t302 - t125 * t300 - ((t126 * t477 + t314 * t341) * t314 + (t125 * t477 + t312 * t341) * t312) * qJD(2) + (t564 + t594) * t61) * m(4) + (t312 * t562 - t314 * t563 + t593) * t442 + (-t474 + 0.1e1) * (rSges(3,1) * t317 - rSges(3,2) * t316) * t642; 0.2e1 * (t27 * t546 + t5 * t539) * m(6) + 0.2e1 * (t26 * t539 + t40 * t546) * m(5) + 0.2e1 * (t539 * t62 + t546 * t61) * m(4) + 0.2e1 * (m(4) * (qJD(2) * t61 + t312 * t91 + t314 * t92) / 0.2e1 + m(5) * (qJD(2) * t40 + t312 * t34 + t314 * t35) / 0.2e1 + m(6) * (qJD(2) * t27 + t312 * t6 + t314 * t7) / 0.2e1) * t316; (t238 * t557 + t239 * t556 - t314 * t555) * t545 + ((t312 * t587 + t314 * t586) * t316 + t606) * t544 + (t236 * t557 + t237 * t556 - t312 * t555) * t543 + ((t312 * t589 + t314 * t588) * t316 + t605) * t542 + (qJD(4) * t558 + t283 * t562 + t284 * t563) * t539 + t592 * t498 / 0.2e1 + t496 * t603 + t585 * t471 / 0.2e1 + ((t312 * t563 + t314 * t562) * t316 + t558) * t442 + (t596 * t317 + (t305 * t557 + t306 * t556) * t316) * t441 + (t316 * t595 - t317 * t602) * t424 + (t626 * t316 - t317 * t653) * t386 + (t627 * t316 - t317 * t652) * t385 + ((qJD(2) * t322 - t43 * t536 + t44 * t537 - t492 * t6 + t493 * t7) * t317 + (t359 * qJD(2) + (-t43 * t522 + t571 * t6 + t361) * t314 + (t44 * t522 - t571 * t7 - t551) * t312) * t316 - (t237 * t43 + t239 * t44 + t27 * t501) * qJD(5) - (t27 * t490 + t44 * t476) * t284 - (t27 * t491 - t43 * t476) * t283 - (t43 * t490 + t44 * t491) * t467) * m(6) + ((t35 * t109 - t34 * t111 - t59 * t90 + t60 * t88 + (t366 + (t312 * t60 - t314 * t59) * t225) * qJD(2)) * t317 + (t60 * (-qJD(2) * t109 + t122 * t312) + t59 * (qJD(2) * t111 - t122 * t314) + t26 * t393 + t40 * (-t312 * t90 + t314 * t88) + (t312 * t35 - t314 * t34) * t225) * t316 - t60 * (t152 * t467 + t256 * t284) - t59 * (-t156 * t467 - t256 * t283) - t40 * (t152 * t283 - t156 * t284)) * m(5) + t593 * t470 / 0.2e1; (t236 * t6 + t238 * t7 + t5 * t502 + (-t236 * t467 - t284 * t502 - t168) * t44 + (t238 * t467 + t283 * t502 - t166) * t43 + (-t236 * t283 + t284 * t238 + t362) * t27) * m(6);];
tauc = t1(:);
