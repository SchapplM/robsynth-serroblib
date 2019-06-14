% Calculate time derivative of joint inertia matrix for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14V3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:29
% EndTime: 2019-04-12 15:04:50
% DurationCPUTime: 36.76s
% Computational Cost: add. (75688->1306), mult. (209703->1828), div. (0->0), fcn. (233794->10), ass. (0->530)
t428 = sin(qJ(2));
t432 = cos(qJ(2));
t618 = Icges(4,5) * t432;
t624 = Icges(3,4) * t432;
t659 = -t618 + t624 + (Icges(3,1) + Icges(4,1)) * t428;
t429 = sin(qJ(1));
t423 = t429 ^ 2;
t433 = cos(qJ(1));
t424 = t433 ^ 2;
t597 = t423 + t424;
t658 = t659 * qJD(2);
t657 = 2 * m(3);
t656 = -t429 / 0.2e1;
t641 = t429 / 0.2e1;
t640 = -t433 / 0.2e1;
t639 = t433 / 0.2e1;
t655 = -qJD(1) / 0.2e1;
t654 = qJD(1) / 0.2e1;
t431 = cos(qJ(4));
t582 = qJD(4) * t431;
t544 = t428 * t582;
t427 = sin(qJ(4));
t589 = qJD(2) * t432;
t549 = t427 * t589;
t444 = t544 + t549;
t425 = sin(qJ(6));
t430 = cos(qJ(6));
t638 = cos(qJ(5));
t557 = t428 * t638;
t426 = sin(qJ(5));
t603 = t432 * t426;
t362 = t431 * t557 - t603;
t437 = -qJD(6) * t362 + t444;
t537 = t638 * qJD(5);
t538 = qJD(2) * t638;
t584 = qJD(4) * t427;
t285 = (t431 * t538 - t537) * t432 + (-t638 * t584 + (-qJD(5) * t431 + qJD(2)) * t426) * t428;
t609 = t427 * t428;
t518 = qJD(6) * t609 + t285;
t222 = -t425 * t518 + t430 * t437;
t223 = t425 * t437 + t430 * t518;
t522 = t428 * t537;
t583 = qJD(4) * t428;
t545 = t427 * t583;
t284 = -t426 * t545 - qJD(5) * t603 - t428 * t538 + (t426 * t589 + t522) * t431;
t144 = Icges(7,5) * t223 + Icges(7,6) * t222 + Icges(7,3) * t284;
t147 = Icges(7,4) * t223 + Icges(7,2) * t222 + Icges(7,6) * t284;
t150 = Icges(7,1) * t223 + Icges(7,4) * t222 + Icges(7,5) * t284;
t326 = -t362 * t425 + t430 * t609;
t327 = t362 * t430 + t425 * t609;
t607 = t428 * t431;
t361 = t426 * t607 + t432 * t638;
t240 = Icges(7,5) * t327 + Icges(7,6) * t326 + Icges(7,3) * t361;
t243 = Icges(7,4) * t327 + Icges(7,2) * t326 + Icges(7,6) * t361;
t246 = Icges(7,1) * t327 + Icges(7,4) * t326 + Icges(7,5) * t361;
t47 = t361 * t144 + t326 * t147 + t327 * t150 + t222 * t243 + t223 * t246 + t284 * t240;
t206 = Icges(6,5) * t285 - Icges(6,6) * t284 + Icges(6,3) * t444;
t209 = Icges(6,4) * t285 - Icges(6,2) * t284 + Icges(6,6) * t444;
t212 = Icges(6,1) * t285 - Icges(6,4) * t284 + Icges(6,5) * t444;
t290 = Icges(6,5) * t362 - Icges(6,6) * t361 + Icges(6,3) * t609;
t293 = Icges(6,4) * t362 - Icges(6,2) * t361 + Icges(6,6) * t609;
t296 = Icges(6,1) * t362 - Icges(6,4) * t361 + Icges(6,5) * t609;
t77 = t206 * t609 - t361 * t209 + t362 * t212 - t284 * t293 + t285 * t296 + t290 * t444;
t651 = -t47 - t77;
t505 = -Icges(3,2) * t428 + t624;
t347 = Icges(3,6) * t429 + t433 * t505;
t625 = Icges(3,4) * t428;
t510 = Icges(3,1) * t432 - t625;
t352 = Icges(3,5) * t429 + t433 * t510;
t470 = t347 * t428 - t352 * t432;
t454 = t470 * t429;
t346 = -Icges(3,6) * t433 + t429 * t505;
t351 = -Icges(3,5) * t433 + t429 * t510;
t471 = t346 * t428 - t351 * t432;
t455 = t471 * t433;
t500 = Icges(4,3) * t428 + t618;
t340 = Icges(4,6) * t429 + t433 * t500;
t619 = Icges(4,5) * t428;
t508 = Icges(4,1) * t432 + t619;
t350 = Icges(4,4) * t429 + t433 * t508;
t472 = t340 * t428 + t350 * t432;
t456 = t472 * t429;
t339 = -Icges(4,6) * t433 + t429 * t500;
t349 = -Icges(4,4) * t433 + t429 * t508;
t473 = t339 * t428 + t349 * t432;
t457 = t473 * t433;
t547 = t429 * t589;
t592 = qJD(1) * t433;
t446 = t428 * t592 + t547;
t501 = Icges(3,5) * t432 - Icges(3,6) * t428;
t341 = -Icges(3,3) * t433 + t429 * t501;
t503 = Icges(4,4) * t432 + Icges(4,6) * t428;
t344 = -Icges(4,2) * t433 + t429 * t503;
t650 = m(7) / 0.2e1;
t602 = t432 * t433;
t559 = t431 * t602;
t366 = t429 * t427 + t559;
t606 = t428 * t433;
t331 = t366 * t638 + t426 * t606;
t365 = t427 * t602 - t429 * t431;
t274 = -t331 * t425 + t365 * t430;
t275 = t331 * t430 + t365 * t425;
t330 = t366 * t426 - t433 * t557;
t111 = t240 * t330 + t243 * t274 + t246 * t275;
t605 = t429 * t432;
t364 = -t427 * t433 + t431 * t605;
t328 = t364 * t426 - t429 * t557;
t608 = t428 * t429;
t329 = t364 * t638 + t426 * t608;
t604 = t431 * t433;
t363 = t427 * t605 + t604;
t272 = -t329 * t425 + t363 * t430;
t273 = t329 * t430 + t363 * t425;
t192 = Icges(7,5) * t273 + Icges(7,6) * t272 + Icges(7,3) * t328;
t194 = Icges(7,4) * t273 + Icges(7,2) * t272 + Icges(7,6) * t328;
t196 = Icges(7,1) * t273 + Icges(7,4) * t272 + Icges(7,5) * t328;
t95 = t192 * t330 + t194 * t274 + t196 * t275;
t193 = Icges(7,5) * t275 + Icges(7,6) * t274 + Icges(7,3) * t330;
t195 = Icges(7,4) * t275 + Icges(7,2) * t274 + Icges(7,6) * t330;
t197 = Icges(7,1) * t275 + Icges(7,4) * t274 + Icges(7,5) * t330;
t96 = t193 * t330 + t195 * t274 + t197 * t275;
t37 = t111 * t361 + t328 * t95 + t330 * t96;
t649 = -t37 / 0.2e1;
t105 = t192 * t361 + t194 * t326 + t196 * t327;
t106 = t193 * t361 + t195 * t326 + t197 * t327;
t124 = t240 * t361 + t243 * t326 + t246 * t327;
t44 = t105 * t328 + t106 * t330 + t124 * t361;
t648 = t44 / 0.2e1;
t467 = (-qJD(4) * t432 + qJD(1)) * t427;
t593 = qJD(1) * t432;
t526 = -qJD(4) + t593;
t588 = qJD(2) * t433;
t546 = t428 * t588;
t287 = t433 * t467 + (-t429 * t526 - t546) * t431;
t523 = t432 * t538;
t524 = qJD(1) * t557;
t224 = qJD(5) * t331 + t287 * t426 + t429 * t524 - t433 * t523;
t647 = t224 / 0.2e1;
t591 = qJD(2) * t428;
t289 = t526 * t604 + (-t431 * t591 + t467) * t429;
t226 = qJD(5) * t329 + t289 * t426 - t429 * t523 - t433 * t524;
t646 = t226 / 0.2e1;
t645 = t284 / 0.2e1;
t644 = t328 / 0.2e1;
t643 = t330 / 0.2e1;
t642 = t361 / 0.2e1;
t388 = rSges(3,1) * t428 + rSges(3,2) * t432;
t637 = m(3) * t388;
t52 = t105 * t363 + t106 * t365 + t124 * t609;
t241 = Icges(6,5) * t329 - Icges(6,6) * t328 + Icges(6,3) * t363;
t244 = Icges(6,4) * t329 - Icges(6,2) * t328 + Icges(6,6) * t363;
t247 = Icges(6,1) * t329 - Icges(6,4) * t328 + Icges(6,5) * t363;
t137 = t241 * t609 - t244 * t361 + t247 * t362;
t242 = Icges(6,5) * t331 - Icges(6,6) * t330 + Icges(6,3) * t365;
t245 = Icges(6,4) * t331 - Icges(6,2) * t330 + Icges(6,6) * t365;
t248 = Icges(6,1) * t331 - Icges(6,4) * t330 + Icges(6,5) * t365;
t138 = t242 * t609 - t245 * t361 + t248 * t362;
t179 = t290 * t609 - t293 * t361 + t296 * t362;
t73 = t137 * t363 + t138 * t365 + t179 * t609;
t636 = t52 + t73;
t635 = rSges(4,1) * t428;
t634 = rSges(4,1) * t432;
t633 = rSges(4,2) * t433;
t421 = t429 * rSges(4,2);
t632 = t429 * rSges(3,3);
t590 = qJD(2) * t429;
t548 = t428 * t590;
t594 = qJD(1) * t429;
t288 = -t427 * t548 - t433 * t584 - t431 * t594 + (t427 * t592 + t429 * t582) * t432;
t208 = Icges(5,5) * t289 - Icges(5,6) * t288 + Icges(5,3) * t446;
t211 = Icges(5,4) * t289 - Icges(5,2) * t288 + Icges(5,6) * t446;
t214 = Icges(5,1) * t289 - Icges(5,4) * t288 + Icges(5,5) * t446;
t291 = Icges(5,5) * t364 - Icges(5,6) * t363 + Icges(5,3) * t608;
t294 = Icges(5,4) * t364 - Icges(5,2) * t363 + Icges(5,6) * t608;
t297 = Icges(5,1) * t364 - Icges(5,4) * t363 + Icges(5,5) * t608;
t481 = -t294 * t427 + t297 * t431;
t86 = (qJD(2) * t481 - t208) * t432 + (qJD(2) * t291 - t211 * t427 + t214 * t431 + (-t294 * t431 - t297 * t427) * qJD(4)) * t428;
t631 = t432 * t86;
t286 = qJD(1) * t363 - qJD(4) * t559 + t427 * t546 - t429 * t584;
t543 = t432 * t588;
t551 = t428 * t594;
t445 = t543 - t551;
t207 = Icges(5,5) * t287 + Icges(5,6) * t286 + Icges(5,3) * t445;
t210 = Icges(5,4) * t287 + Icges(5,2) * t286 + Icges(5,6) * t445;
t213 = Icges(5,1) * t287 + Icges(5,4) * t286 + Icges(5,5) * t445;
t292 = Icges(5,5) * t366 - Icges(5,6) * t365 + Icges(5,3) * t606;
t295 = Icges(5,4) * t366 - Icges(5,2) * t365 + Icges(5,6) * t606;
t298 = Icges(5,1) * t366 - Icges(5,4) * t365 + Icges(5,5) * t606;
t480 = -t295 * t427 + t298 * t431;
t87 = (qJD(2) * t480 - t207) * t432 + (qJD(2) * t292 - t210 * t427 + t213 * t431 + (-t295 * t431 - t298 * t427) * qJD(4)) * t428;
t630 = t432 * t87;
t629 = t433 * rSges(3,3);
t628 = -rSges(4,3) - qJ(3);
t627 = -rSges(5,3) - qJ(3);
t626 = t124 * t284 + t47 * t361;
t622 = Icges(5,4) * t427;
t621 = Icges(5,4) * t431;
t184 = -t291 * t432 + t428 * t481;
t612 = t184 * t432;
t185 = -t292 * t432 + t428 * t480;
t611 = t185 * t432;
t506 = Icges(5,1) * t431 - t622;
t348 = -Icges(5,5) * t432 + t428 * t506;
t610 = t348 * t431;
t585 = qJD(3) * t433;
t601 = qJ(3) * t543 + t428 * t585;
t414 = qJ(3) * t602;
t586 = qJD(3) * t429;
t600 = qJD(1) * t414 + t432 * t586;
t599 = rSges(4,2) * t592 + rSges(4,3) * t543;
t598 = t597 * qJ(3) * t428;
t342 = Icges(3,3) * t429 + t433 * t501;
t596 = qJD(1) * t342;
t345 = Icges(4,2) * t429 + t433 * t503;
t595 = qJD(1) * t345;
t587 = qJD(3) * t428;
t227 = t429 * t522 + t289 * t638 + (-qJD(5) * t364 + t446) * t426;
t164 = -qJD(6) * t273 - t227 * t425 + t288 * t430;
t165 = qJD(6) * t272 + t227 * t430 + t288 * t425;
t100 = Icges(7,4) * t165 + Icges(7,2) * t164 + Icges(7,6) * t226;
t102 = Icges(7,1) * t165 + Icges(7,4) * t164 + Icges(7,5) * t226;
t225 = t433 * t522 + t287 * t638 + (-qJD(5) * t366 + t445) * t426;
t162 = -qJD(6) * t275 - t225 * t425 - t286 * t430;
t163 = qJD(6) * t274 + t225 * t430 - t286 * t425;
t98 = Icges(7,5) * t165 + Icges(7,6) * t164 + Icges(7,3) * t226;
t22 = t100 * t274 + t102 * t275 + t162 * t194 + t163 * t196 + t192 * t224 + t330 * t98;
t101 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t224;
t97 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t224;
t99 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t224;
t23 = t101 * t275 + t162 * t195 + t163 * t197 + t193 * t224 + t274 * t99 + t330 * t97;
t32 = t144 * t330 + t147 * t274 + t150 * t275 + t162 * t243 + t163 * t246 + t224 * t240;
t511 = t429 * t95 + t433 * t96;
t63 = t96 * t429 - t433 * t95;
t1 = (qJD(2) * t511 - t32) * t432 + (-qJD(1) * t63 + qJD(2) * t111 + t22 * t429 + t23 * t433) * t428;
t167 = t290 * t365 - t293 * t330 + t296 * t331;
t146 = Icges(6,5) * t227 - Icges(6,6) * t226 + Icges(6,3) * t288;
t149 = Icges(6,4) * t227 - Icges(6,2) * t226 + Icges(6,6) * t288;
t152 = Icges(6,1) * t227 - Icges(6,4) * t226 + Icges(6,5) * t288;
t48 = t146 * t365 - t149 * t330 + t152 * t331 - t224 * t244 + t225 * t247 - t241 * t286;
t145 = Icges(6,5) * t225 - Icges(6,6) * t224 - Icges(6,3) * t286;
t148 = Icges(6,4) * t225 - Icges(6,2) * t224 - Icges(6,6) * t286;
t151 = Icges(6,1) * t225 - Icges(6,4) * t224 - Icges(6,5) * t286;
t49 = t145 * t365 - t148 * t330 + t151 * t331 - t224 * t245 + t225 * t248 - t242 * t286;
t127 = t241 * t365 - t244 * t330 + t247 * t331;
t128 = t242 * t365 - t245 * t330 + t248 * t331;
t492 = t127 * t429 + t128 * t433;
t493 = t127 * t433 - t128 * t429;
t70 = t206 * t365 - t209 * t330 + t212 * t331 - t224 * t293 + t225 * t296 - t286 * t290;
t13 = (qJD(2) * t492 - t70) * t432 + (qJD(1) * t493 + qJD(2) * t167 + t429 * t48 + t433 * t49) * t428;
t579 = t13 / 0.2e1 + t1 / 0.2e1;
t166 = t290 * t363 - t293 * t328 + t296 * t329;
t125 = t241 * t363 - t244 * t328 + t247 * t329;
t126 = t242 * t363 - t245 * t328 + t248 * t329;
t494 = t125 * t429 + t126 * t433;
t495 = t125 * t433 - t126 * t429;
t50 = t146 * t363 - t149 * t328 + t152 * t329 - t226 * t244 + t227 * t247 + t241 * t288;
t51 = t145 * t363 - t148 * t328 + t151 * t329 - t226 * t245 + t227 * t248 + t242 * t288;
t71 = t206 * t363 - t209 * t328 + t212 * t329 - t226 * t293 + t227 * t296 + t288 * t290;
t14 = (qJD(2) * t494 - t71) * t432 + (qJD(1) * t495 + qJD(2) * t166 + t429 * t50 + t433 * t51) * t428;
t110 = t240 * t328 + t243 * t272 + t246 * t273;
t24 = t100 * t272 + t102 * t273 + t164 * t194 + t165 * t196 + t192 * t226 + t328 * t98;
t25 = t101 * t273 + t164 * t195 + t165 * t197 + t193 * t226 + t272 * t99 + t328 * t97;
t33 = t144 * t328 + t147 * t272 + t150 * t273 + t164 * t243 + t165 * t246 + t226 * t240;
t93 = t192 * t328 + t194 * t272 + t196 * t273;
t94 = t193 * t328 + t195 * t272 + t197 * t273;
t512 = t429 * t93 + t433 * t94;
t62 = t94 * t429 - t433 * t93;
t2 = (qJD(2) * t512 - t33) * t432 + (-qJD(1) * t62 + qJD(2) * t110 + t24 * t429 + t25 * t433) * t428;
t578 = t14 / 0.2e1 + t2 / 0.2e1;
t15 = t127 * t288 - t128 * t286 + t167 * t444 + t363 * t48 + t365 * t49 + t609 * t70;
t3 = t111 * t444 + t22 * t363 + t23 * t365 - t286 * t96 + t288 * t95 + t32 * t609;
t577 = t3 / 0.2e1 + t15 / 0.2e1;
t16 = t125 * t288 - t126 * t286 + t166 * t444 + t363 * t50 + t365 * t51 + t609 * t71;
t4 = t110 * t444 + t24 * t363 + t25 * t365 - t286 * t94 + t288 * t93 + t33 * t609;
t576 = -t4 / 0.2e1 - t16 / 0.2e1;
t10 = qJD(1) * t511 - t22 * t433 + t23 * t429;
t19 = qJD(1) * t492 + t49 * t429 - t433 * t48;
t575 = t10 / 0.2e1 + t19 / 0.2e1;
t11 = qJD(1) * t512 - t24 * t433 + t25 * t429;
t20 = qJD(1) * t494 + t51 * t429 - t433 * t50;
t574 = t11 / 0.2e1 + t20 / 0.2e1;
t28 = t100 * t326 + t102 * t327 + t192 * t284 + t194 * t222 + t196 * t223 + t361 * t98;
t29 = t101 * t327 + t193 * t284 + t195 * t222 + t197 * t223 + t326 * t99 + t361 * t97;
t496 = t105 * t429 + t106 * t433;
t12 = qJD(1) * t496 - t28 * t433 + t29 * t429;
t490 = t137 * t429 + t138 * t433;
t54 = t146 * t609 - t149 * t361 + t152 * t362 + t241 * t444 - t244 * t284 + t247 * t285;
t55 = t145 * t609 - t148 * t361 + t151 * t362 + t242 * t444 - t245 * t284 + t248 * t285;
t573 = t490 * t655 + t55 * t656 + t54 * t639 - t12 / 0.2e1;
t572 = t32 / 0.2e1 + t29 / 0.2e1;
t571 = t33 / 0.2e1 + t28 / 0.2e1;
t38 = t110 * t609 + t363 * t93 + t365 * t94;
t65 = t125 * t363 + t126 * t365 + t166 * t609;
t570 = t38 / 0.2e1 + t65 / 0.2e1;
t39 = t111 * t609 + t363 * t95 + t365 * t96;
t66 = t127 * t363 + t128 * t365 + t167 * t609;
t569 = t39 / 0.2e1 + t66 / 0.2e1;
t53 = -t124 * t432 + t428 * t496;
t74 = -t179 * t432 + t428 * t490;
t568 = t53 / 0.2e1 + t74 / 0.2e1;
t567 = t62 / 0.2e1 - t495 / 0.2e1;
t566 = -t63 / 0.2e1 + t493 / 0.2e1;
t491 = t137 * t433 - t138 * t429;
t497 = t105 * t433 - t106 * t429;
t565 = -t497 / 0.2e1 - t491 / 0.2e1;
t40 = -t110 * t432 + t428 * t512;
t68 = -t166 * t432 + t428 * t494;
t564 = t68 / 0.2e1 + t40 / 0.2e1;
t41 = -t111 * t432 + t428 * t511;
t69 = -t167 * t432 + t428 * t492;
t563 = -t69 / 0.2e1 - t41 / 0.2e1;
t562 = t124 * t444 + t47 * t609;
t561 = t179 * t444 + t77 * t609;
t560 = qJ(3) * t608;
t103 = t163 * rSges(7,1) + t162 * rSges(7,2) + t224 * rSges(7,3);
t154 = t225 * rSges(6,1) - t224 * rSges(6,2) - t286 * rSges(6,3);
t199 = t275 * rSges(7,1) + t274 * rSges(7,2) + t330 * rSges(7,3);
t558 = t287 * rSges(5,1) + t286 * rSges(5,2) + rSges(5,3) * t543;
t251 = t331 * rSges(6,1) - t330 * rSges(6,2) + t365 * rSges(6,3);
t301 = t366 * rSges(5,1) - t365 * rSges(5,2) + rSges(5,3) * t606;
t356 = rSges(4,1) * t602 + rSges(4,3) * t606 + t421;
t556 = qJ(3) * t591;
t249 = rSges(7,1) * t327 + rSges(7,2) * t326 + rSges(7,3) * t361;
t554 = t249 * t592;
t299 = rSges(6,1) * t362 - rSges(6,2) * t361 + rSges(6,3) * t609;
t553 = t299 * t592;
t513 = rSges(5,1) * t431 - rSges(5,2) * t427;
t353 = -rSges(5,3) * t432 + t428 * t513;
t552 = t353 * t592;
t542 = qJD(2) * t648;
t541 = t110 / 0.2e1 + t105 / 0.2e1;
t540 = t111 / 0.2e1 + t106 / 0.2e1;
t539 = (t501 + t503) * qJD(2) / 0.2e1;
t536 = t429 * t627;
t533 = 2 * m(4);
t531 = 2 * m(5);
t529 = 2 * m(6);
t527 = 0.2e1 * m(7);
t525 = t597 * (qJ(3) * t589 + t587);
t521 = t428 * t536;
t182 = t291 * t606 - t365 * t294 + t366 * t297;
t183 = t292 * t606 - t365 * t295 + t366 * t298;
t487 = t182 * t433 - t183 * t429;
t520 = -t487 / 0.2e1 - t566;
t180 = t291 * t608 - t294 * t363 + t297 * t364;
t181 = t292 * t608 - t295 * t363 + t298 * t364;
t489 = t180 * t433 - t181 * t429;
t519 = -t489 / 0.2e1 + t567;
t517 = rSges(3,1) * t432 - rSges(3,2) * t428;
t516 = rSges(4,3) * t428 + t634;
t387 = -rSges(4,3) * t432 + t635;
t515 = t289 * rSges(5,1) - t288 * rSges(5,2);
t514 = -rSges(5,1) * t364 + rSges(5,2) * t363;
t504 = Icges(3,2) * t432 + t625;
t502 = -Icges(5,2) * t427 + t621;
t499 = -Icges(4,3) * t432 + t619;
t498 = Icges(5,5) * t431 - Icges(5,6) * t427;
t488 = t180 * t429 + t181 * t433;
t486 = t182 * t429 + t183 * t433;
t198 = rSges(7,1) * t273 + rSges(7,2) * t272 + rSges(7,3) * t328;
t485 = t198 * t433 - t199 * t429;
t484 = -t429 * t198 - t199 * t433;
t250 = rSges(6,1) * t329 - rSges(6,2) * t328 + rSges(6,3) * t363;
t483 = t250 * t433 - t251 * t429;
t482 = -t429 * t250 - t251 * t433;
t300 = rSges(5,3) * t608 - t514;
t479 = t300 * t433 - t301 * t429;
t478 = -t429 * t300 - t301 * t433;
t338 = -Icges(5,3) * t432 + t428 * t498;
t343 = -Icges(5,6) * t432 + t428 * t502;
t228 = t338 * t608 - t343 * t363 + t348 * t364;
t107 = -t228 * t432 + t428 * t488;
t469 = t107 + t40 + t68 - t612;
t229 = t338 * t606 - t365 * t343 + t366 * t348;
t108 = -t229 * t432 + t428 * t486;
t468 = -t108 - t41 - t69 + t611;
t466 = t55 / 0.2e1 + t70 / 0.2e1 + t572;
t465 = t71 / 0.2e1 + t54 / 0.2e1 + t571;
t36 = t110 * t361 + t328 * t93 + t330 * t94;
t464 = t36 * t641 + t37 * t639;
t463 = t166 / 0.2e1 + t137 / 0.2e1 + t541;
t462 = -t138 / 0.2e1 - t167 / 0.2e1 - t540;
t153 = rSges(7,1) * t223 + rSges(7,2) * t222 + rSges(7,3) * t284;
t461 = -t153 * t433 + t249 * t594;
t215 = rSges(6,1) * t285 - rSges(6,2) * t284 + rSges(6,3) * t444;
t460 = -t215 * t433 + t299 * t594;
t317 = (-rSges(5,1) * t427 - rSges(5,2) * t431) * t583 + (rSges(5,3) * t428 + t432 * t513) * qJD(2);
t459 = -t317 * t433 + t353 * t594;
t458 = t428 * t628 - t634;
t453 = -qJ(3) * t551 + t601;
t450 = qJD(2) * t504;
t449 = qJD(2) * (-Icges(4,4) * t428 + Icges(4,6) * t432);
t448 = qJD(2) * (-Icges(3,5) * t428 - Icges(3,6) * t432);
t447 = qJD(2) * t499;
t155 = t227 * rSges(6,1) - t226 * rSges(6,2) + t288 * rSges(6,3);
t104 = t165 * rSges(7,1) + t164 * rSges(7,2) + t226 * rSges(7,3);
t443 = qJD(1) * t458;
t302 = (-Icges(5,5) * t427 - Icges(5,6) * t431) * t583 + (Icges(5,3) * t428 + t432 * t498) * qJD(2);
t307 = (-Icges(5,2) * t431 - t622) * t583 + (Icges(5,6) * t428 + t432 * t502) * qJD(2);
t312 = (-Icges(5,1) * t427 - t621) * t583 + (Icges(5,5) * t428 + t432 * t506) * qJD(2);
t116 = t286 * t343 + t287 * t348 + t302 * t606 - t365 * t307 + t366 * t312 + t338 * t445;
t442 = t116 / 0.2e1 + t87 / 0.2e1 + t466;
t117 = -t288 * t343 + t289 * t348 + t302 * t608 - t363 * t307 + t364 * t312 + t338 * t446;
t441 = t86 / 0.2e1 + t117 / 0.2e1 + t465;
t440 = t228 / 0.2e1 + t184 / 0.2e1 + t463;
t439 = t229 / 0.2e1 + t185 / 0.2e1 - t462;
t438 = t429 * t570 + t433 * t569;
t436 = -t432 * t302 + t312 * t607 + t338 * t591 - t343 * t544 + t589 * t610;
t435 = t432 * t585 + (-t429 * t593 - t546) * qJ(3);
t434 = -qJ(3) * t446 - t428 * t586;
t413 = qJ(3) * t606;
t411 = qJ(3) * t605;
t374 = t517 * qJD(2);
t373 = t516 * qJD(2);
t357 = t433 * t517 + t632;
t355 = rSges(3,1) * t605 - rSges(3,2) * t608 - t629;
t354 = t429 * t516 - t633;
t337 = -t387 * t433 + t414;
t336 = -t387 * t429 + t411;
t334 = t413 + t356;
t333 = t429 * t458 + t633;
t321 = -t353 * t433 + t414;
t320 = -t353 * t429 + t411;
t319 = -rSges(3,1) * t548 + (rSges(3,1) * t602 + t632) * qJD(1) - t446 * rSges(3,2);
t318 = -t388 * t588 + (-t429 * t517 + t629) * qJD(1);
t309 = t429 * t449 + t595;
t308 = -qJD(1) * t344 + t433 * t449;
t306 = t429 * t448 + t596;
t305 = -qJD(1) * t341 + t433 * t448;
t279 = t413 + t301;
t278 = t521 + t514;
t269 = -t299 * t433 + t414;
t268 = -t299 * t429 + t411;
t267 = -t387 * t592 + (-t373 - t556) * t429 + t600;
t266 = -t373 * t433 + t387 * t594 + t435;
t265 = t429 * t354 + t356 * t433 + t598;
t264 = t433 * t443 + (-rSges(4,2) * qJD(1) - t587 + (t432 * t628 + t635) * qJD(2)) * t429;
t263 = -rSges(4,1) * t546 + t429 * t443 + t599 + t601;
t262 = -t432 * t301 - t353 * t606;
t261 = t300 * t432 + t353 * t608;
t260 = t429 * t342 - t433 * t470;
t259 = t429 * t341 - t455;
t258 = t429 * t345 + t433 * t472;
t257 = t429 * t344 + t457;
t256 = -t342 * t433 - t454;
t255 = -t341 * t433 - t429 * t471;
t254 = -t345 * t433 + t456;
t253 = -t344 * t433 + t429 * t473;
t252 = -t338 * t432 + (-t343 * t427 + t610) * t428;
t239 = t252 * t591;
t238 = t413 + t251;
t237 = -t250 - t560;
t236 = -t249 * t433 + t414;
t235 = -t249 * t429 + t411;
t232 = -t552 + (-t317 - t556) * t429 + t600;
t231 = t435 + t459;
t230 = t479 * t428;
t218 = -t478 + t598;
t217 = rSges(5,3) * t446 + t515;
t216 = -rSges(5,3) * t551 + t558;
t201 = -t432 * t251 - t299 * t606;
t200 = t250 * t432 + t299 * t608;
t191 = t536 * t589 + (t592 * t627 - t586) * t428 - t515;
t190 = qJD(1) * t521 + t558 + t601;
t189 = t251 * t609 - t299 * t365;
t188 = -t250 * t609 + t299 * t363;
t187 = t413 + t199;
t186 = -t198 - t560;
t178 = t483 * t428;
t177 = t433 * t599 + (-t387 * t423 - t424 * t635) * qJD(2) + (t433 * t354 + (-t356 + t421) * t429) * qJD(1) + t525;
t176 = t179 * t591;
t173 = -t482 + t598;
t170 = t250 * t365 - t251 * t363;
t169 = -t553 + (-t215 - t556) * t429 + t600;
t168 = t435 + t460;
t159 = -t432 * t199 - t249 * t606;
t158 = t198 * t432 + t249 * t608;
t157 = t199 * t609 - t249 * t365;
t156 = -t198 * t609 + t249 * t363;
t143 = t199 * t361 - t249 * t330;
t142 = -t198 * t361 + t249 * t328;
t141 = (t353 * t590 + t217) * t432 + (-qJD(2) * t300 + t429 * t317 + t552) * t428;
t140 = (-t353 * t588 - t216) * t432 + (qJD(2) * t301 + t459) * t428;
t136 = (-t343 * t589 + (-qJD(4) * t348 - t307) * t428) * t427 + t436;
t135 = -t155 + t434;
t134 = t453 + t154;
t131 = t485 * t428;
t130 = -t484 + t598;
t129 = t198 * t365 - t199 * t363;
t123 = t124 * t591;
t120 = t198 * t330 - t199 * t328;
t114 = qJD(1) * t479 + t216 * t433 + t429 * t217 + t525;
t113 = -t554 + (-t153 - t556) * t429 + t600;
t112 = t435 + t461;
t109 = t479 * t589 + (qJD(1) * t478 - t216 * t429 + t217 * t433) * t428;
t92 = -t104 + t434;
t91 = t453 + t103;
t89 = (t299 * t590 + t155) * t432 + (-qJD(2) * t250 + t429 * t215 + t553) * t428;
t88 = (-t299 * t588 - t154) * t432 + (qJD(2) * t251 + t460) * t428;
t83 = -t155 * t609 + t215 * t363 - t250 * t444 + t288 * t299;
t82 = t154 * t609 - t215 * t365 + t251 * t444 + t286 * t299;
t81 = t207 * t608 - t363 * t210 + t364 * t213 - t288 * t295 + t289 * t298 + t292 * t446;
t80 = t208 * t608 - t363 * t211 + t364 * t214 - t288 * t294 + t289 * t297 + t291 * t446;
t79 = t207 * t606 - t365 * t210 + t366 * t213 + t286 * t295 + t287 * t298 + t292 * t445;
t78 = t208 * t606 - t365 * t211 + t366 * t214 + t286 * t294 + t287 * t297 + t291 * t445;
t75 = qJD(1) * t483 + t154 * t433 + t429 * t155 + t525;
t72 = -t154 * t363 + t155 * t365 - t250 * t286 - t251 * t288;
t67 = t483 * t589 + (qJD(1) * t482 - t154 * t429 + t155 * t433) * t428;
t61 = (t249 * t590 + t104) * t432 + (-qJD(2) * t198 + t429 * t153 + t554) * t428;
t60 = (-t249 * t588 - t103) * t432 + (qJD(2) * t199 + t461) * t428;
t59 = -t104 * t361 + t153 * t328 - t198 * t284 + t226 * t249;
t58 = t103 * t361 - t153 * t330 + t199 * t284 - t224 * t249;
t57 = -t104 * t609 + t153 * t363 - t198 * t444 + t249 * t288;
t56 = t103 * t609 - t153 * t365 + t199 * t444 + t249 * t286;
t43 = qJD(1) * t485 + t103 * t433 + t429 * t104 + t525;
t42 = -t103 * t363 + t104 * t365 - t198 * t286 - t199 * t288;
t35 = -t103 * t328 + t104 * t330 + t198 * t224 - t199 * t226;
t34 = t485 * t589 + (qJD(1) * t484 - t103 * t429 + t104 * t433) * t428;
t31 = qJD(1) * t488 + t81 * t429 - t433 * t80;
t30 = qJD(1) * t486 + t79 * t429 - t433 * t78;
t27 = (qJD(2) * t488 - t117) * t432 + (qJD(1) * t489 + qJD(2) * t228 + t429 * t80 + t433 * t81) * t428;
t26 = (qJD(2) * t486 - t116) * t432 + (qJD(1) * t487 + qJD(2) * t229 + t429 * t78 + t433 * t79) * t428;
t18 = t137 * t288 - t138 * t286 + t54 * t363 + t55 * t365 + t561;
t17 = t176 + (qJD(2) * t490 - t77) * t432 + (qJD(1) * t491 + t54 * t429 + t55 * t433) * t428;
t9 = t105 * t226 + t106 * t224 + t28 * t328 + t29 * t330 + t626;
t8 = t105 * t288 - t106 * t286 + t28 * t363 + t29 * t365 + t562;
t7 = t123 + (qJD(2) * t496 - t47) * t432 + (qJD(1) * t497 + t28 * t429 + t29 * t433) * t428;
t6 = t110 * t284 + t224 * t94 + t226 * t93 + t24 * t328 + t25 * t330 + t33 * t361;
t5 = t111 * t284 + t22 * t328 + t224 * t96 + t226 * t95 + t23 * t330 + t32 * t361;
t21 = [(t134 * t238 + t135 * t237) * t529 + (t186 * t92 + t187 * t91) * t527 + (t190 * t279 + t191 * t278) * t531 + (t263 * t334 + t264 * t333) * t533 + (t318 * t357 + t319 * t355) * t657 + t436 - t343 * t549 - t348 * t545 - t307 * t609 + (-t504 + t499 + t510 + t508) * t591 + (t505 - t500 + t659) * t589 - t651; (t263 * t336 + t264 * t337 + t266 * t333 + t267 * t334) * m(4) + (t190 * t320 + t191 * t321 + t231 * t278 + t232 * t279) * m(5) + (t134 * t268 + t135 * t269 + t168 * t237 + t169 * t238) * m(6) + (t112 * t186 + t113 * t187 + t235 * t91 + t236 * t92) * m(7) + ((t319 * t388 + t355 * t374) * m(3) + t539 * t433 + (t340 * t654 + t347 * t655 + t447 * t656 + t450 * t641) * t432 + ((t350 + t352) * t655 + t658 * t641) * t428 - t441) * t433 + ((-t318 * t388 - t357 * t374) * m(3) + t539 * t429 + (t339 * t654 + t346 * t655 + t447 * t639 + t450 * t640) * t432 + ((t349 + t351) * t655 + t658 * t640) * t428 + t442) * t429 + (-t457 / 0.2e1 + t455 / 0.2e1 + t456 / 0.2e1 - t454 / 0.2e1) * qJD(2) + ((-t357 * t637 + (t347 / 0.2e1 - t340 / 0.2e1) * t432 + (t352 / 0.2e1 + t350 / 0.2e1) * t428 + t439) * t433 + (-t355 * t637 + (-t339 / 0.2e1 + t346 / 0.2e1) * t432 + (t349 / 0.2e1 + t351 / 0.2e1) * t428 + t440) * t429) * qJD(1); (t114 * t218 + t231 * t321 + t232 * t320) * t531 + (t177 * t265 + t266 * t337 + t267 * t336) * t533 + (t168 * t269 + t169 * t268 + t173 * t75) * t529 + (t112 * t236 + t113 * t235 + t130 * t43) * t527 - t433 * t31 - t433 * t20 - t433 * t11 + t429 * t30 + t429 * t19 + t429 * t10 - t433 * ((t433 * t309 + (t254 - t457) * qJD(1)) * t433 + (t253 * qJD(1) + (t340 * t589 - t350 * t591 + t595) * t429 + (-t308 + (-t339 * t432 + t349 * t428) * qJD(2) + t472 * qJD(1)) * t433) * t429) - t433 * ((t433 * t306 + (t256 + t455) * qJD(1)) * t433 + (t255 * qJD(1) + (-t347 * t589 - t352 * t591 + t596) * t429 + (-t305 + (t346 * t432 + t351 * t428) * qJD(2) - t470 * qJD(1)) * t433) * t429) + ((t429 * t355 + t357 * t433) * (t433 * t318 + t429 * t319 + (t355 * t433 - t357 * t429) * qJD(1)) + t597 * t388 * t374) * t657 + t429 * ((t429 * t305 + (t259 + t454) * qJD(1)) * t429 + (t260 * qJD(1) + (t346 * t589 + t351 * t591) * t433 + (-t306 + (-t347 * t432 - t352 * t428) * qJD(2) + (t342 - t471) * qJD(1)) * t429) * t433) + t429 * ((t429 * t308 + (t257 - t456) * qJD(1)) * t429 + (t258 * qJD(1) + (-t339 * t589 + t349 * t591) * t433 + (-t309 + (t340 * t432 - t350 * t428) * qJD(2) + (t345 + t473) * qJD(1)) * t429) * t433) + (-t489 + t62 - t495 + (-t253 - t255) * t433 + (t254 + t256) * t429) * t594 + (-t487 + t63 - t493 + (-t257 - t259) * t433 + (t258 + t260) * t429) * t592; (m(4) * t333 + m(5) * t278 + m(6) * t237 + m(7) * t186) * t543 + (m(4) * t334 + m(5) * t279 + m(6) * t238 + m(7) * t187) * t547 + ((-t186 * t594 + t187 * t592 + t429 * t91 + t433 * t92) * m(7) + (t134 * t429 + t135 * t433 - t237 * t594 + t238 * t592) * m(6) + (t190 * t429 + t191 * t433 - t278 * t594 + t279 * t592) * m(5) + (t263 * t429 + t264 * t433 - t333 * t594 + t334 * t592) * m(4)) * t428; (-m(4) * t177 - m(5) * t114 - m(6) * t75 - m(7) * t43 + ((m(4) * t337 + m(5) * t321 + m(6) * t269 + m(7) * t236) * t433 + (m(4) * t336 + m(5) * t320 + m(6) * t268 + m(7) * t235) * t429) * qJD(2)) * t432 + ((t112 * t433 + t113 * t429 + t235 * t592 - t236 * t594) * m(7) + (t168 * t433 + t169 * t429 + t268 * t592 - t269 * t594) * m(6) + (t231 * t433 + t232 * t429 + t320 * t592 - t321 * t594) * m(5) + (t266 * t433 + t267 * t429 + t336 * t592 - t337 * t594) * m(4) + (m(4) * t265 + m(5) * t218 + m(6) * t173 + m(7) * t130) * qJD(2)) * t428; 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + t650) * (-0.1e1 + t597) * t428 * t589; (t140 * t279 + t141 * t278 + t190 * t262 + t191 * t261) * m(5) + (t134 * t201 + t135 * t200 + t237 * t89 + t238 * t88) * m(6) + (t158 * t92 + t159 * t91 + t186 * t61 + t187 * t60) * m(7) + t239 + t123 + t176 + (-t136 + (t429 * t440 + t433 * t439) * qJD(2) + t651) * t432 + (t442 * t433 + t441 * t429 + (-t429 * t439 + t433 * t440) * qJD(1)) * t428; (t109 * t218 + t114 * t230 + t140 * t320 + t141 * t321 + t231 * t261 + t232 * t262) * m(5) + (t168 * t200 + t169 * t201 + t173 * t67 + t178 * t75 + t268 * t88 + t269 * t89) * m(6) + (t112 * t158 + t113 * t159 + t130 * t34 + t131 * t43 + t235 * t60 + t236 * t61) * m(7) + t573 * t432 + (t631 / 0.2e1 - t27 / 0.2e1 + (-t611 / 0.2e1 + t108 / 0.2e1 - t563) * qJD(1) + t520 * t589 - t578) * t433 + (-t630 / 0.2e1 + t26 / 0.2e1 + (-t612 / 0.2e1 + t107 / 0.2e1 + t564) * qJD(1) + t519 * t589 + t579) * t429 + ((t30 / 0.2e1 + t575) * t433 + (t31 / 0.2e1 + t574) * t429 + (t184 * t640 + t185 * t641 + t565) * qJD(2) + (-t429 * t520 + t433 * t519) * qJD(1)) * t428; (-m(5) * t109 - m(6) * t67 - m(7) * t34 + ((m(5) * t261 + m(6) * t200 + m(7) * t158) * t433 + (m(5) * t262 + m(6) * t201 + m(7) * t159) * t429) * qJD(2)) * t432 + ((t140 * t429 + t141 * t433 - t261 * t594 + t262 * t592) * m(5) + (-t200 * t594 + t201 * t592 + t429 * t88 + t433 * t89) * m(6) + (-t158 * t594 + t159 * t592 + t429 * t60 + t433 * t61) * m(7) + (m(5) * t230 + m(6) * t178 + m(7) * t131) * qJD(2)) * t428; (t109 * t230 + t140 * t262 + t141 * t261) * t531 + (t178 * t67 + t200 * t89 + t201 * t88) * t529 + (t131 * t34 + t158 * t61 + t159 * t60) * t527 + (t136 * t432 - t17 - t239 - t7 + (t429 * t469 - t433 * t468) * qJD(2)) * t432 + ((t1 + t13 + t26 - t630) * t433 + (t14 + t2 + t27 - t631) * t429 + (-t252 * t432 + t74 + t53 + (t184 * t429 + t185 * t433) * t428) * qJD(2) + (t429 * t468 + t433 * t469) * qJD(1)) * t428; (t156 * t92 + t157 * t91 + t186 * t57 + t187 * t56) * m(7) + (t134 * t189 + t135 * t188 + t237 * t83 + t238 * t82) * m(6) + t466 * t365 + t465 * t363 + t463 * t288 + t462 * t286 + t561 + t562; (t112 * t156 + t113 * t157 + t129 * t43 + t130 * t42 + t235 * t56 + t236 * t57) * m(7) + (t168 * t188 + t169 * t189 + t170 * t75 + t173 * t72 + t268 * t82 + t269 * t83) * m(6) + t576 * t433 + t577 * t429 + t575 * t365 + t574 * t363 + t567 * t288 + t566 * t286 + t565 * t544 + (-t428 * t573 + t565 * t589) * t427 + t438 * qJD(1); (-m(6) * t72 - m(7) * t42 + ((m(6) * t188 + m(7) * t156) * t433 + (m(6) * t189 + m(7) * t157) * t429) * qJD(2)) * t432 + ((-t188 * t594 + t189 * t592 + t429 * t82 + t433 * t83) * m(6) + (-t156 * t594 + t157 * t592 + t429 * t56 + t433 * t57) * m(7) + (m(6) * t170 + m(7) * t129) * qJD(2)) * t428; (t170 * t67 + t178 * t72 + t188 * t89 + t189 * t88 + t200 * t83 + t201 * t82) * m(6) + (t129 * t34 + t131 * t42 + t156 * t61 + t157 * t60 + t158 * t57 + t159 * t56) * m(7) + t579 * t365 + t578 * t363 + t564 * t288 + t563 * t286 + (-t18 / 0.2e1 - t8 / 0.2e1 + (t427 * t568 + t438) * qJD(2)) * t432 + (t577 * t433 - t576 * t429 + (t17 / 0.2e1 + t7 / 0.2e1) * t427 + t568 * t582 + (t73 / 0.2e1 + t52 / 0.2e1) * qJD(2) + (-t429 * t569 + t433 * t570) * qJD(1)) * t428; (t129 * t42 + t156 * t57 + t157 * t56) * t527 + (t170 * t72 + t188 * t83 + t189 * t82) * t529 + (t3 + t15) * t365 + (t4 + t16) * t363 + (t38 + t65) * t288 + (-t39 - t66) * t286 + t636 * t544 + ((t18 + t8) * t428 + t636 * t589) * t427; (t142 * t92 + t143 * t91 + t186 * t59 + t187 * t58) * m(7) + t572 * t330 + t571 * t328 + t541 * t226 + t540 * t224 + t626; t62 * t646 + t11 * t644 - t497 * t645 + t12 * t642 + t6 * t640 + (t112 * t142 + t113 * t143 + t120 * t43 + t130 * t35 + t235 * t58 + t236 * t59) * m(7) + t5 * t641 + t63 * t647 + t10 * t643 + t464 * qJD(1); 0.2e1 * ((-t35 + (t142 * t433 + t143 * t429) * qJD(2)) * t432 + (qJD(2) * t120 + t429 * t58 + t433 * t59 + (-t142 * t429 + t143 * t433) * qJD(1)) * t428) * t650; (t120 * t34 + t131 * t35 + t142 * t61 + t143 * t60 + t158 * t59 + t159 * t58) * m(7) + t41 * t647 + t1 * t643 + t53 * t645 + t7 * t642 + t40 * t646 + t2 * t644 + (-t9 / 0.2e1 + t464 * qJD(2)) * t432 + (t6 * t641 + t5 * t639 + t542 + (t36 * t639 + t429 * t649) * qJD(1)) * t428; t544 * t648 + (t120 * t42 + t129 * t35 + t142 * t57 + t143 * t56 + t156 * t59 + t157 * t58) * m(7) + t286 * t649 + t365 * t5 / 0.2e1 + t38 * t646 + t4 * t644 + t39 * t647 + t3 * t643 + t52 * t645 + t8 * t642 + t288 * t36 / 0.2e1 + t363 * t6 / 0.2e1 + (t432 * t542 + t428 * t9 / 0.2e1) * t427; t224 * t37 + t330 * t5 + t226 * t36 + t328 * t6 + t284 * t44 + t361 * t9 + (t120 * t35 + t142 * t59 + t143 * t58) * t527;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;