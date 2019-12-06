% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:15
% EndTime: 2019-12-05 17:59:38
% DurationCPUTime: 10.03s
% Computational Cost: add. (17405->402), mult. (22437->531), div. (0->0), fcn. (20608->6), ass. (0->255)
t619 = Icges(5,4) + Icges(6,4);
t613 = Icges(5,2) + Icges(6,2);
t617 = Icges(5,1) + Icges(6,1);
t377 = qJ(3) + qJ(4);
t359 = sin(t377);
t618 = t619 * t359;
t614 = Icges(5,5) + Icges(6,5);
t616 = Icges(5,6) + Icges(6,6);
t360 = cos(t377);
t596 = t613 * t360 + t618;
t615 = t619 * t360;
t379 = sin(qJ(1));
t457 = t360 * t379;
t611 = t619 * t457;
t610 = -t617 * t359 - t615;
t381 = cos(qJ(1));
t600 = t614 * t359 + t616 * t360;
t609 = t381 * t600;
t608 = t600 * t379;
t607 = Icges(6,3) + Icges(5,3);
t605 = t596 * t379 + t616 * t381;
t463 = t359 * t379;
t603 = t614 * t381 + t617 * t463 + t611;
t604 = -t379 * t616 + t596 * t381;
t602 = t614 * t379 + t610 * t381;
t606 = -t613 * t359 + t615;
t598 = t617 * t360 - t618;
t585 = -t607 * t379 + t609;
t601 = -t603 * t359 - t605 * t360;
t302 = rSges(5,1) * t457 - rSges(5,2) * t463;
t380 = cos(qJ(3));
t455 = t379 * t380;
t357 = pkin(3) * t455;
t233 = t302 + t357;
t333 = rSges(5,1) * t360 - rSges(5,2) * t359;
t492 = pkin(3) * t380;
t415 = (t333 + t492) * t381;
t557 = t379 * t415;
t164 = -t233 * t381 + t557;
t378 = sin(qJ(3));
t345 = rSges(4,1) * t380 - rSges(4,2) * t378;
t316 = t345 * t379;
t317 = t345 * t381;
t188 = -t316 * t381 + t317 * t379;
t301 = rSges(6,1) * t457 - rSges(6,2) * t463;
t491 = pkin(4) * t360;
t337 = t491 + t492;
t198 = t379 * t337 + t301;
t486 = rSges(6,1) * t360;
t351 = t381 * t486;
t485 = rSges(6,2) * t359;
t199 = t351 + (t337 - t485) * t381;
t537 = m(6) / 0.2e1;
t538 = m(5) / 0.2e1;
t539 = m(4) / 0.2e1;
t422 = (-t198 * t381 + t199 * t379) * t537 + t164 * t538 + t188 * t539;
t332 = -t485 + t486;
t352 = pkin(4) * t457;
t226 = t379 * t332 + t352;
t195 = t357 + t226;
t197 = (t332 + t337) * t381;
t237 = t333 * t379 + t357;
t451 = (t195 * t381 - t197 * t379) * t537 + (t237 * t381 - t557) * t538;
t29 = t451 - t422;
t599 = t29 * qJD(1);
t597 = -t359 * t616 + t614 * t360;
t595 = t606 * t381 - t602;
t594 = t613 * t463 - t603 - t611;
t593 = -t598 * t381 + t604;
t592 = t598 * t379 - t605;
t591 = (t602 * t359 - t604 * t360) * t381;
t590 = t606 - t610;
t413 = rSges(5,1) * t359 + rSges(5,2) * t360;
t375 = t379 ^ 2;
t376 = t381 ^ 2;
t427 = t375 + t376;
t589 = t413 * t427;
t493 = pkin(3) * t378;
t566 = t413 + t493;
t238 = t566 * t381;
t588 = t603 * t463 + t605 * t457 + (t607 * t381 + t608) * t381;
t587 = -t585 * t381 - t604 * t457 + t602 * t463;
t586 = -t585 * t379 - t591;
t462 = t359 * t381;
t262 = t381 * (rSges(6,2) * t462 - t351);
t148 = -t301 * t379 - t427 * t491 + t262;
t412 = rSges(6,1) * t359 + rSges(6,2) * t360;
t308 = t412 * t379;
t225 = -pkin(4) * t463 - t308;
t353 = pkin(4) * t462;
t573 = t412 * t381;
t227 = t353 + t573;
t448 = t195 * t225 - t197 * t227;
t382 = -pkin(7) - pkin(6);
t424 = t381 * t493;
t265 = t381 * (t424 + (pkin(6) + t382) * t379);
t358 = t381 * t382;
t456 = t378 * t379;
t307 = pkin(3) * t456 - pkin(6) * t381 - t358;
t425 = -qJ(5) + t382;
t356 = t381 * t425;
t336 = pkin(4) * t359 + t493;
t569 = t336 + t412;
t445 = -rSges(6,3) * t381 + t356 - t358 + (t493 - t569) * t379;
t568 = t381 * t336 + t573 + (-rSges(6,3) + t425) * t379;
t446 = (t379 * t382 + t424 - t568) * t381;
t75 = -t265 + (-t307 + t445) * t379 + t446;
t33 = t75 * t148 + t448;
t372 = t379 * rSges(5,3);
t468 = t413 * t381;
t224 = t381 * (-t372 + t468);
t469 = t413 * t379;
t255 = rSges(5,3) * t381 + t469;
t103 = -t224 - t265 + (-t255 - t307) * t379;
t303 = t333 * t381;
t181 = -t302 * t379 - t381 * t303;
t52 = t103 * t181 - t237 * t469 - t415 * t468;
t581 = -m(5) * t52 - m(6) * t33;
t580 = t601 * t381;
t579 = -(t593 * t379 + t592 * t381) * t359 + (t595 * t379 + t594 * t381) * t360;
t361 = t381 * qJ(2);
t171 = -t379 * pkin(1) + t361 + t568;
t172 = -t356 + (rSges(6,3) + pkin(1)) * t381 + (qJ(2) + t569) * t379;
t578 = m(6) * (-t171 * t379 + t381 * t172);
t574 = -t379 / 0.2e1;
t517 = t379 / 0.2e1;
t515 = t381 / 0.2e1;
t475 = Icges(4,4) * t380;
t340 = -Icges(4,2) * t378 + t475;
t411 = Icges(4,1) * t378 + t475;
t572 = (t411 / 0.2e1 + t340 / 0.2e1) * t380;
t571 = t597 * t379;
t570 = t597 * t381;
t228 = (t332 + t491) * t381;
t182 = t361 - t372 + (-pkin(1) + t382) * t379 + t238;
t183 = -t358 + (rSges(5,3) + pkin(1)) * t381 + (qJ(2) + t566) * t379;
t416 = -t182 * t469 + t183 * t468;
t452 = t225 * t171 + t227 * t172;
t489 = (-t198 * t228 + t199 * t226 + t452) * t537 + (t164 * t333 + t416) * t538;
t215 = t301 + t352;
t216 = t351 + (-t485 + t491) * t381;
t490 = (t195 * t216 - t197 * t215 + t452) * t537 + (t237 * t303 - t302 * t415 + t416) * t538;
t565 = t489 - t490;
t561 = (t596 - t598) * t360 + t590 * t359;
t287 = -Icges(4,5) * t379 + t411 * t381;
t433 = t340 * t381 + t287;
t354 = Icges(4,4) * t455;
t286 = Icges(4,1) * t456 + Icges(4,5) * t381 + t354;
t434 = -Icges(4,2) * t456 + t286 + t354;
t476 = Icges(4,4) * t378;
t408 = Icges(4,2) * t380 + t476;
t285 = -Icges(4,6) * t379 + t408 * t381;
t342 = Icges(4,1) * t380 - t476;
t435 = -t342 * t381 + t285;
t284 = Icges(4,6) * t381 + t408 * t379;
t436 = -t342 * t379 + t284;
t560 = -(t379 * t435 - t381 * t436) * t378 + (t379 * t433 - t381 * t434) * t380;
t558 = qJD(1) * t578;
t556 = (t285 * t380 + t287 * t378) * t381;
t552 = t375 / 0.2e1 + t376 / 0.2e1;
t180 = -t302 * t381 + t303 * t379;
t449 = (-t215 * t381 + t216 * t379) * t537 + t180 * t538;
t386 = -t590 * t360 / 0.2e1 + (-t598 / 0.2e1 + t596 / 0.2e1) * t359;
t390 = (t587 * t379 + t588 * t381) * t574 + ((t586 + t588 + t591) * t381 + ((t601 + t585) * t381 - t580 + t587) * t379) * t517 + (t586 * t379 + t585 * t375 + ((-t601 + t585) * t381 + t580 + t587) * t381) * t515;
t543 = 0.4e1 * qJD(1);
t542 = 2 * qJD(3);
t540 = 2 * qJD(4);
t535 = -pkin(1) - pkin(6);
t170 = -t255 * t379 - t224;
t78 = t170 * t181 - t333 * t589;
t76 = m(5) * t78;
t529 = m(5) * (t182 * t415 + t183 * t233);
t528 = m(5) * (t182 * t303 + t183 * t302);
t104 = t376 * (-t337 + t492) + t262 + (t357 - t198) * t379;
t447 = t226 * t225 - t228 * t227;
t84 = t445 * t379 + t446;
t524 = m(6) * (t104 * t84 + t447);
t521 = m(6) * (t171 * t199 + t172 * t198);
t520 = m(6) * (t171 * t216 + t172 * t215);
t519 = m(6) * (t171 * t381 + t379 * t172);
t516 = -t381 / 0.2e1;
t514 = m(3) * ((rSges(3,3) * t381 + t361) * t381 + (rSges(3,3) + qJ(2)) * t375);
t414 = rSges(4,1) * t378 + rSges(4,2) * t380;
t387 = -t379 * rSges(4,3) + t414 * t381;
t200 = t535 * t379 + t361 + t387;
t201 = (rSges(4,3) - t535) * t381 + (qJ(2) + t414) * t379;
t513 = m(4) * (t200 * t317 + t201 * t316);
t512 = m(4) * (t200 * t381 + t201 * t379);
t510 = m(5) * (t182 * t381 + t183 * t379);
t502 = m(6) * (t198 * t379 + t199 * t381);
t501 = m(6) * (-t195 * t379 - t381 * t197);
t498 = m(6) * (t215 * t379 + t216 * t381);
t497 = m(6) * (t226 * t381 - t228 * t379);
t495 = m(6) * (-t226 * t379 - t381 * t228);
t487 = m(6) * qJD(3);
t240 = 0.2e1 * t552 * m(6);
t426 = t240 * qJD(1);
t423 = t84 * t148 + t447;
t404 = Icges(4,5) * t378 + Icges(4,6) * t380;
t282 = Icges(4,3) * t381 + t404 * t379;
t121 = t381 * t282 + t284 * t455 + t286 * t456;
t283 = -Icges(4,3) * t379 + t404 * t381;
t122 = -t381 * t283 - t285 * t455 - t287 * t456;
t420 = ((t571 * t379 + t579) * t381 - t570 * t375) * t517 + ((-t570 * t381 - t579) * t379 + t571 * t376) * t515;
t419 = t427 * t414;
t418 = t427 * t492;
t417 = t76 + t420;
t405 = Icges(4,5) * t380 - Icges(4,6) * t378;
t236 = t566 * t379;
t399 = -t236 * t379 - t238 * t381;
t392 = -t284 * t380 - t286 * t378;
t384 = -t390 + (t595 * t359 + t593 * t360 + t561 * t381 - t608) * t517 + (t594 * t359 + t592 * t360 - t561 * t379 - t609) * t515;
t383 = -t386 + (t359 * t604 + t360 * t602) * (t515 + t516);
t311 = t405 * t381;
t310 = t379 * t405;
t258 = t379 * t282;
t196 = t353 + (t412 + t493) * t381;
t194 = -t336 * t379 - t308;
t163 = t181 - t418;
t158 = t225 * t381 + t227 * t379;
t141 = t495 / 0.2e1;
t139 = t497 / 0.2e1;
t134 = m(6) * t158 * qJD(4);
t133 = t498 / 0.2e1;
t124 = -t283 * t379 + t556;
t123 = t392 * t381 + t258;
t119 = t501 / 0.2e1;
t118 = t502 / 0.2e1;
t98 = m(6) * (t225 * t379 - t227 * t381) - m(5) * t589;
t97 = t98 * qJD(4);
t92 = -t418 + t104;
t80 = t123 * t381 + t124 * t379;
t79 = t121 * t381 + t122 * t379;
t59 = t141 - t498 / 0.2e1;
t58 = t141 + t133;
t57 = t133 - t495 / 0.2e1;
t51 = t119 + t118;
t50 = t119 - t502 / 0.2e1;
t49 = t118 - t501 / 0.2e1;
t47 = t139 - t449;
t46 = -t497 / 0.2e1 + t449;
t45 = t139 + t449;
t34 = t510 + t512 + t514 + t519;
t27 = t422 + t451;
t26 = t386 + t520 + t528;
t25 = t283 * t375 + (t122 - t258 + (t283 - t392) * t381) * t381;
t24 = (-t123 + t258 + t122) * t379 + (t124 - t556 + (t283 + t392) * t379 + t121) * t381;
t10 = -t572 + (-t342 / 0.2e1 + t408 / 0.2e1) * t378 + t513 + t529 + t521 + t386;
t7 = t417 + t524;
t6 = t420 - t581;
t4 = (t25 / 0.2e1 + t80 / 0.2e1) * t381 + (-t79 / 0.2e1 + t24 / 0.2e1) * t379 + t390;
t3 = t390 + t565;
t2 = t390 - t565;
t1 = t384 + t489 + t490;
t5 = [t34 * qJD(2) + t10 * qJD(3) + t26 * qJD(4) + qJD(5) * t578, qJD(1) * t34 + qJD(3) * t27 + qJD(4) * t45, t10 * qJD(1) + t27 * qJD(2) + t1 * qJD(4) + t51 * qJD(5) + ((t188 * t345 - (t200 * t379 - t201 * t381) * t414) * t539 + (t171 * t194 + t172 * t196 + t195 * t199 - t197 * t198) * t537 + (-t182 * t236 + t183 * t238 + (-t233 + t237) * t415) * t538) * t542 + (t384 + (-t378 * t434 - t380 * t436) * t515 + t24 * t574 + (t378 * t433 + t380 * t435 + t79) * t517 + (t25 + t80) * t516 - t552 * t404) * qJD(3), t26 * qJD(1) + t45 * qJD(2) + t1 * qJD(3) + t384 * qJD(4) + t58 * qJD(5) + ((-t215 * t228 + t216 * t226 + t452) * t537 + (t180 * t333 + t416) * t538) * t540, t51 * qJD(3) + t58 * qJD(4) + t558; -t29 * qJD(3) + t46 * qJD(4) - t240 * qJD(5) + (-t519 / 0.4e1 - t510 / 0.4e1 - t512 / 0.4e1 - t514 / 0.4e1) * t543, 0, -t599 + t97 + (t399 * t538 + (t194 * t379 - t196 * t381) * t537 - t419 * t539) * t542, qJD(1) * t46 + qJD(3) * t98 + t97, -t426; (t383 + (-t408 + t342) * t378 / 0.2e1 + t572) * qJD(1) + t29 * qJD(2) + t4 * qJD(3) + t2 * qJD(4) + t50 * qJD(5) + (-t521 / 0.4e1 - t529 / 0.4e1 - t513 / 0.4e1) * t543, t599, t4 * qJD(1) + (m(4) * ((-t381 * t387 + (-t381 * rSges(4,3) - t414 * t379) * t379) * (-t316 * t379 - t317 * t381) - t345 * t419) + (t376 * t310 + (-t381 * t311 - t560) * t379) * t515 + (-t375 * t311 + (t379 * t310 + t560) * t381) * t517 + m(6) * (t194 * t195 - t196 * t197 + t75 * t92) + m(5) * (t103 * t163 - t236 * t237 - t238 * t415) + t420) * qJD(3) + t6 * qJD(4), t2 * qJD(1) + t6 * qJD(3) + ((t423 + t33) * t537 + (t52 + t78) * t538) * t540 + (t420 - t76 - t524) * qJD(4), t50 * qJD(1); t383 * qJD(1) + t47 * qJD(2) + t3 * qJD(3) + t390 * qJD(4) + t59 * qJD(5) + (-t520 / 0.4e1 - t528 / 0.4e1) * t543, t47 * qJD(1), t3 * qJD(1) + t7 * qJD(4) + ((t104 * t75 + t194 * t226 - t196 * t228 + t84 * t92 + t448) * t537 + (t163 * t170 + t333 * t399 + t52) * t538) * t542 + (t420 + t581) * qJD(3), t390 * qJD(1) + t7 * qJD(3) + (m(6) * t423 + t417) * qJD(4), t59 * qJD(1); t240 * qJD(2) + t49 * qJD(3) + t57 * qJD(4) - t558, t426, t49 * qJD(1) + (t194 * t381 + t196 * t379) * t487 + t134, t57 * qJD(1) + t158 * t487 + t134, 0;];
Cq = t5;
