% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP3
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:56
% EndTime: 2019-12-31 19:51:05
% DurationCPUTime: 5.84s
% Computational Cost: add. (24827->348), mult. (18716->440), div. (0->0), fcn. (16467->8), ass. (0->242)
t360 = qJ(1) + qJ(2);
t357 = sin(t360);
t359 = pkin(8) + qJ(4);
t356 = cos(t359);
t355 = sin(t359);
t491 = rSges(6,1) + pkin(4);
t395 = t491 * t355;
t465 = rSges(6,3) + qJ(5);
t567 = -t465 * t356 + t395;
t573 = t567 * t357;
t575 = t357 * t573;
t358 = cos(t360);
t435 = t355 * t358;
t574 = t435 * t573;
t340 = Icges(5,4) * t356;
t314 = Icges(5,1) * t355 + t340;
t458 = Icges(6,5) * t356;
t572 = Icges(6,1) * t355 + t314 - t458;
t339 = Icges(6,5) * t355;
t377 = Icges(6,3) * t356 - t339;
t459 = Icges(5,4) * t355;
t568 = Icges(5,2) * t356 + t377 + t459;
t401 = qJD(1) + qJD(2);
t243 = t567 * t358;
t524 = m(6) / 0.2e1;
t318 = rSges(5,1) * t355 + rSges(5,2) * t356;
t353 = t357 ^ 2;
t354 = t358 ^ 2;
t552 = t353 + t354;
t543 = t552 * t318;
t423 = (-t243 * t358 - t575) * t524 - m(5) * t543 / 0.2e1;
t283 = t318 * t357;
t285 = t318 * t358;
t199 = t283 * t357 + t285 * t358;
t429 = t356 * t358;
t222 = -t358 * t395 + t465 * t429;
t451 = t222 * t358;
t525 = m(5) / 0.2e1;
t424 = (-t451 + t575) * t524 + t199 * t525;
t43 = t424 - t423;
t571 = t401 * t43;
t351 = t358 * rSges(6,2);
t472 = -pkin(7) - qJ(3);
t392 = t358 * t472;
t361 = cos(pkin(8));
t352 = pkin(3) * t361 + pkin(2);
t538 = t465 * t355 + t491 * t356;
t533 = t352 + t538;
t189 = -t533 * t357 + t351 - t392;
t338 = t357 * t472;
t466 = t357 * rSges(6,2);
t190 = t533 * t358 - t338 + t466;
t104 = t189 * t358 + t190 * t357;
t550 = rSges(4,2) * sin(pkin(8)) - rSges(4,1) * t361 - pkin(2);
t554 = rSges(4,3) + qJ(3);
t245 = t550 * t357 + t554 * t358;
t473 = sin(qJ(1)) * pkin(1);
t228 = t245 - t473;
t246 = t554 * t357 - t550 * t358;
t474 = cos(qJ(1)) * pkin(1);
t229 = t246 + t474;
t149 = t228 * t358 + t229 * t357;
t156 = t245 * t358 + t246 * t357;
t467 = rSges(5,1) * t356;
t391 = t352 + t467;
t436 = t355 * t357;
t404 = rSges(5,2) * t436 + t358 * rSges(5,3);
t217 = -t357 * t391 - t392 + t404;
t390 = -rSges(5,2) * t435 + t357 * rSges(5,3);
t218 = t358 * t391 - t338 + t390;
t536 = -t217 * t358 - t218 * t357;
t210 = t217 - t473;
t211 = t218 + t474;
t541 = -t210 * t358 - t211 * t357;
t559 = m(4) / 0.2e1;
t185 = t189 - t473;
t186 = t190 + t474;
t99 = t185 * t358 + t186 * t357;
t399 = (t104 + t99) * t524 + (-t536 - t541) * t525 + (t156 + t149) * t559;
t400 = (t99 - t104) * t524 + (-t541 + t536) * t525 + (t149 - t156) * t559;
t6 = t400 - t399;
t570 = t6 * qJD(1);
t311 = -Icges(5,2) * t355 + t340;
t569 = t311 + t572;
t315 = Icges(5,1) * t356 - t459;
t537 = Icges(6,1) * t356 + t339;
t566 = t537 + t315;
t241 = t538 * t357;
t564 = t185 - t189;
t562 = (Icges(5,6) - Icges(6,6)) * t356 + (Icges(6,4) + Icges(5,5)) * t355;
t262 = Icges(6,4) * t357 + t358 * t537;
t264 = Icges(5,5) * t357 + t315 * t358;
t561 = -t568 * t358 + t262 + t264;
t327 = Icges(6,5) * t429;
t254 = Icges(6,6) * t357 + Icges(6,3) * t435 + t327;
t260 = Icges(5,6) * t357 + t311 * t358;
t560 = -Icges(6,1) * t435 - t314 * t358 + t254 - t260 + t327;
t521 = m(4) * (-t246 * t228 + t229 * t245);
t489 = m(3) * (t474 * (-rSges(3,1) * t357 - rSges(3,2) * t358) + (t358 * rSges(3,1) - t357 * rSges(3,2)) * t473);
t397 = t190 * t435;
t97 = -t189 * t436 + t397;
t555 = m(6) * qJD(2) * t97;
t309 = Icges(6,4) * t356 + Icges(6,6) * t355;
t444 = t309 * t357;
t257 = -Icges(6,2) * t358 + t444;
t249 = t357 * t257;
t307 = Icges(6,3) * t355 + t458;
t253 = -Icges(6,6) * t358 + t307 * t357;
t261 = -Icges(6,4) * t358 + t357 * t537;
t131 = t253 * t435 + t261 * t429 + t249;
t553 = t131 * t358;
t551 = (t566 - t568) * t356 + (t307 - t569) * t355;
t328 = Icges(5,4) * t436;
t430 = t356 * t357;
t263 = Icges(5,1) * t430 - Icges(5,5) * t358 - t328;
t548 = -Icges(5,2) * t430 - t377 * t357 + t261 + t263 - t328;
t259 = Icges(5,4) * t430 - Icges(5,2) * t436 - Icges(5,6) * t358;
t547 = t572 * t357 - t253 + t259;
t545 = (t253 * t355 + t261 * t356) * t357;
t332 = t353 * t355;
t333 = t354 * t355;
t212 = 0.2e1 * (t332 / 0.2e1 + t333 / 0.2e1) * m(6);
t544 = t401 * t212;
t244 = t538 * t358;
t58 = -t190 * t185 + t186 * t189;
t88 = -t218 * t210 + t211 * t217;
t540 = t562 * t357;
t539 = t562 * t358;
t471 = ((-t186 + t190) * t243 + t564 * t573) * t524 + ((-t211 + t218) * t358 + (t210 - t217) * t357) * t318 * t525;
t102 = t210 * t283 - t211 * t285;
t107 = t217 * t283 - t218 * t285;
t75 = t185 * t573 + t186 * t222;
t82 = t189 * t573 + t190 * t222;
t534 = (t82 + t75) * t524 + (t107 + t102) * t525;
t370 = (-t307 / 0.2e1 + t569 / 0.2e1) * t356 + (-t568 / 0.2e1 + t566 / 0.2e1) * t355;
t532 = (t548 * t355 + t547 * t356) * t358 + (-t561 * t355 + t560 * t356) * t357;
t531 = 0.4e1 * qJD(1);
t529 = 0.4e1 * qJD(2);
t528 = 2 * qJD(4);
t516 = m(5) * t88;
t159 = t186 * t435;
t510 = m(6) * (-t564 * t436 + t159 - t397);
t509 = m(6) * (t159 + t397 + (-t185 - t189) * t436);
t506 = m(6) * t58;
t419 = t222 * t436 + t574;
t422 = t185 * t429 + t186 * t430;
t504 = m(6) * (t419 + t422);
t420 = t189 * t429 + t190 * t430;
t503 = m(6) * (t419 + t420);
t383 = t243 * t436 - t574;
t502 = m(6) * (t383 + t422);
t501 = m(6) * (t383 + t420);
t500 = m(6) * t75;
t499 = m(6) * t82;
t496 = m(6) * t99;
t495 = -t357 / 0.2e1;
t494 = t357 / 0.2e1;
t493 = -t358 / 0.2e1;
t443 = t309 * t358;
t258 = Icges(6,2) * t357 + t443;
t132 = t254 * t435 + t357 * t258 + t262 * t429;
t447 = t257 * t358;
t372 = t132 + t447;
t382 = -t254 * t436 + t258 * t358 - t262 * t430;
t23 = (t132 - t372) * t358 + (t131 + t382 - t249) * t357;
t255 = Icges(5,5) * t430 - Icges(5,6) * t436 - Icges(5,3) * t358;
t417 = -t357 * t255 - t263 * t429;
t133 = -t259 * t435 - t417;
t308 = Icges(5,5) * t356 - Icges(5,6) * t355;
t445 = t308 * t358;
t256 = Icges(5,3) * t357 + t445;
t416 = t357 * t256 + t264 * t429;
t134 = -t260 * t435 + t416;
t384 = t260 * t355 - t255;
t233 = t264 * t430;
t385 = t256 * t358 - t233;
t24 = (t358 * t384 + t134 - t416) * t358 + (t357 * t384 + t133 + t385) * t357;
t127 = -t447 + t545;
t25 = -t553 + (t127 + t372 - t545) * t357;
t130 = -t260 * t436 - t385;
t446 = t259 * t355;
t26 = (t130 - t233 + (t256 + t446) * t358 + t417) * t358 + t416 * t357;
t76 = -t127 * t358 - t357 * t382;
t77 = -(-t357 * (-t263 * t356 + t446) - t255 * t358) * t358 + t130 * t357;
t78 = t132 * t357 - t553;
t79 = -t133 * t358 + t134 * t357;
t2 = (t79 / 0.2e1 - t26 / 0.2e1 + t78 / 0.2e1 - t25 / 0.2e1) * t358 + (t24 / 0.2e1 + t77 / 0.2e1 + t23 / 0.2e1 + t76 / 0.2e1) * t357;
t490 = -qJD(3) * t43 + t2 * qJD(4);
t486 = m(4) * t149;
t485 = m(4) * t156;
t484 = m(5) * t102;
t483 = m(5) * t107;
t482 = m(5) * t541;
t481 = m(5) * t536;
t478 = m(6) * t104;
t469 = m(6) * qJD(4);
t468 = m(6) * qJD(5);
t44 = t423 + t424;
t463 = t44 * qJD(4);
t462 = t43 * qJD(4) - t212 * qJD(5);
t437 = t355 * t356;
t421 = (-t222 - t243) * t573;
t418 = -t243 * t429 - t430 * t573;
t407 = t552 * t437;
t403 = t332 + t333;
t94 = -t185 * t436 + t159;
t398 = m(6) * t94 * qJD(1);
t371 = (-t283 * t358 + t285 * t357) * t318;
t366 = t370 + t534;
t365 = -t370 + ((t253 + t259) * t356 + (-t261 + t263) * t355) * (t494 + t495);
t364 = t44 * qJD(3) + ((t25 + t26) * t358 / 0.2e1 + (t23 + t24 + t76 + t77) * t495 + (t308 * t357 + t560 * t355 + t561 * t356 + t551 * t358 + t444) * t494 + (-t547 * t355 + t548 * t356 + t551 * t357 - t443 - t445 + t78 + t79) * t493) * qJD(4);
t321 = -rSges(5,2) * t355 + t467;
t230 = t407 - t437;
t206 = t212 * qJD(3);
t137 = -t353 * t567 + t451;
t116 = (t466 + t244) * t358 + (-t351 + t241) * t357;
t74 = t116 * t403 + t418;
t68 = t501 / 0.2e1;
t65 = t502 / 0.2e1;
t64 = t503 / 0.2e1;
t59 = t504 / 0.2e1;
t49 = t509 / 0.2e1;
t48 = t510 / 0.2e1;
t47 = t478 - t481 + t485;
t46 = -t482 + t486 + t496;
t30 = t370 + t483 + t499;
t29 = t370 + t484 + t500;
t18 = t489 + t506 + t516 + t521;
t17 = t68 - t503 / 0.2e1;
t16 = t68 + t64;
t15 = t64 - t501 / 0.2e1;
t14 = t65 - t504 / 0.2e1;
t13 = t65 + t59;
t12 = t59 - t502 / 0.2e1;
t11 = t49 - t510 / 0.2e1;
t10 = t49 + t48;
t9 = t48 - t509 / 0.2e1;
t7 = t399 + t400;
t5 = t366 + t471;
t4 = t366 - t471;
t3 = t365 + t471 - t534;
t1 = [t18 * qJD(2) + t46 * qJD(3) + t29 * qJD(4) + t94 * t468, t18 * qJD(1) + t7 * qJD(3) + t5 * qJD(4) + t10 * qJD(5) + 0.2e1 * (t521 / 0.2e1 + t489 / 0.2e1 + t58 * t524 + t88 * t525) * qJD(2), qJD(1) * t46 + qJD(2) * t7 + t463, t29 * qJD(1) + t5 * qJD(2) + t13 * qJD(5) + ((t541 * t321 + t371) * t525 + (-t185 * t244 - t186 * t241 + t421) * t524) * t528 + t364, t10 * qJD(2) + t13 * qJD(4) + t398; -t6 * qJD(3) + t4 * qJD(4) + t11 * qJD(5) + (-t506 / 0.4e1 - t516 / 0.4e1 - t521 / 0.4e1 - t489 / 0.4e1) * t531, t47 * qJD(3) + t30 * qJD(4) + t468 * t97, qJD(2) * t47 + t463 - t570, t4 * qJD(1) + t30 * qJD(2) + t16 * qJD(5) + ((-t189 * t244 - t190 * t241 + t421) * t524 + (t536 * t321 + t371) * t525) * t528 + t364, t11 * qJD(1) + t16 * qJD(4) + t555; t6 * qJD(2) + (-t496 / 0.4e1 + t482 / 0.4e1 - t486 / 0.4e1) * t531 + t462, t570 + (-t478 / 0.4e1 + t481 / 0.4e1 - t485 / 0.4e1) * t529 + t462, 0, (t241 * t358 - t244 * t357) * t469 + t571, -t544; t365 * qJD(1) + t3 * qJD(2) + t14 * qJD(5) + (-t500 / 0.4e1 - t484 / 0.4e1) * t531 + t490, t3 * qJD(1) + t365 * qJD(2) + t17 * qJD(5) + (-t499 / 0.4e1 - t483 / 0.4e1) * t529 + t490, -t571, (m(5) * (t321 * t543 - (t357 * (rSges(5,1) * t430 - t404) + t358 * (rSges(5,1) * t429 + t390)) * t199) + m(6) * (t116 * t137 + t241 * t573 + t243 * t244) + ((t540 * t357 + t532) * t358 - t539 * t353) * t494 + ((t539 * t358 + t532) * t357 - t540 * t354) * t493) * qJD(4) + t74 * t468 + t401 * t2, t14 * qJD(1) + t17 * qJD(2) + t74 * t469 + (-t356 * t403 - t230 + t407) * t468; t9 * qJD(2) + t12 * qJD(4) + t206 - t398, t9 * qJD(1) + t15 * qJD(4) + t206 - t555, t544, t12 * qJD(1) + t15 * qJD(2) + (-t137 * t356 + (-t241 * t357 - t244 * t358 + t116) * t355 - t74 + t418) * t469 + t230 * t468, t230 * t469;];
Cq = t1;
