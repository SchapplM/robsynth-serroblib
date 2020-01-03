% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:49
% EndTime: 2019-12-31 17:08:03
% DurationCPUTime: 9.96s
% Computational Cost: add. (13374->465), mult. (33923->615), div. (0->0), fcn. (37553->6), ass. (0->292)
t320 = sin(qJ(1));
t498 = -t320 / 0.4e1;
t319 = sin(qJ(2));
t321 = cos(qJ(2));
t490 = sin(qJ(4));
t491 = cos(qJ(4));
t267 = t319 * t491 - t321 * t490;
t337 = t319 * t490 + t321 * t491;
t183 = -Icges(5,5) * t337 - Icges(5,6) * t267;
t322 = cos(qJ(1));
t247 = t267 * t322;
t248 = t337 * t322;
t265 = Icges(5,4) * t267;
t187 = -Icges(5,2) * t337 + t265;
t191 = Icges(5,1) * t337 + t265;
t427 = t187 + t191;
t264 = Icges(5,4) * t337;
t188 = Icges(5,2) * t267 + t264;
t190 = Icges(5,1) * t267 - t264;
t429 = t188 - t190;
t79 = t183 * t320 + t429 * t247 + t427 * t248;
t220 = Icges(5,4) * t248;
t159 = Icges(5,2) * t247 - Icges(5,6) * t320 + t220;
t435 = -Icges(5,1) * t247 + t159 + t220;
t219 = Icges(5,4) * t247;
t161 = Icges(5,1) * t248 - Icges(5,5) * t320 + t219;
t526 = -Icges(5,2) * t248 + t161 + t219;
t84 = t435 * t267 + t337 * t526;
t593 = (t79 + t84) * t498;
t581 = (t187 / 0.2e1 + t191 / 0.2e1) * t267;
t592 = (-t190 / 0.2e1 + t188 / 0.2e1) * t337 - t581;
t245 = t267 * t320;
t246 = t337 * t320;
t77 = -t183 * t322 + t429 * t245 + t427 * t246;
t591 = -t77 / 0.2e1;
t590 = -t77 / 0.4e1;
t588 = t84 / 0.2e1 + t79 / 0.2e1;
t218 = Icges(5,4) * t246;
t461 = Icges(5,2) * t245;
t361 = t461 + t218;
t458 = Icges(5,6) * t322;
t334 = t458 + t361;
t217 = Icges(5,4) * t245;
t468 = Icges(5,1) * t246;
t365 = t217 + t468;
t463 = Icges(5,5) * t322;
t335 = t463 + t365;
t350 = -t245 * t159 - t246 * t161;
t357 = Icges(5,5) * t246 + Icges(5,6) * t245;
t332 = Icges(5,3) * t322 + t357;
t533 = t332 - t357;
t318 = t322 ^ 2;
t407 = 0.2e1 * t218;
t360 = t407 + 0.2e1 * t458;
t406 = 0.2e1 * t463;
t323 = (t360 + t461) * t245 + (t406 + t468) * t246 + Icges(5,3) * t318;
t538 = t323 * t322;
t15 = (t533 * t320 + (-t335 + t365) * t248 - (t334 - t361) * t247 + t350) * t320 + t538;
t331 = Icges(5,5) * t248 + Icges(5,6) * t247 - Icges(5,3) * t320;
t328 = t322 * t331;
t91 = t328 - t350;
t92 = t247 * t159 + t248 * t161 - t320 * t331;
t16 = (-t247 * t361 - t248 * t365 - 0.2e1 * t328 + t350 + t91) * t322 + (-t245 * t334 - t246 * t335 - t322 * t533 + t323 + t92) * t320;
t47 = t320 * t91 - t538;
t48 = t320 * t92 - (t247 * t334 + t248 * t335 - t320 * t332) * t322;
t493 = t322 / 0.4e1;
t494 = -t322 / 0.4e1;
t587 = 0.2e1 * t16 * t493 + 0.2e1 * t48 * t494 + 0.2e1 * (t15 + t47) * t498 + (t84 / 0.4e1 + t79 / 0.4e1) * t320;
t193 = rSges(5,1) * t267 - rSges(5,2) * t337;
t454 = qJ(3) * t321;
t280 = pkin(2) * t319 - t454;
t482 = pkin(3) * t319;
t371 = t193 + t280 + t482;
t146 = t371 * t320;
t148 = t371 * t322;
t478 = rSges(4,1) * t319;
t281 = -rSges(4,3) * t321 + t478;
t414 = t280 + t281;
t199 = t414 * t320;
t201 = t414 * t322;
t314 = t322 * rSges(4,2);
t316 = t322 * pkin(5);
t472 = rSges(4,3) + qJ(3);
t492 = rSges(4,1) + pkin(2);
t523 = t472 * t319 + t492 * t321 + pkin(1);
t165 = -t320 * t523 + t314 + t316;
t166 = (rSges(4,2) + pkin(5)) * t320 + t523 * t322;
t444 = t321 * t322;
t445 = t320 * t321;
t437 = t165 * t444 + t166 * t445;
t455 = qJ(3) * t319;
t517 = pkin(2) + pkin(3);
t525 = t517 * t321 + pkin(1) + t455;
t531 = -t246 * rSges(5,1) - t245 * rSges(5,2);
t129 = t316 + (-rSges(5,3) - pkin(6)) * t322 - t525 * t320 + t531;
t164 = t248 * rSges(5,1) + t247 * rSges(5,2) - t320 * rSges(5,3);
t532 = -t320 * pkin(6) + t164;
t130 = t320 * pkin(5) + t322 * t525 + t532;
t440 = t129 * t444 + t130 * t445;
t446 = t319 * t322;
t448 = t319 * t320;
t518 = m(5) / 0.2e1;
t519 = m(4) / 0.2e1;
t480 = (-t146 * t446 + t148 * t448 + t440) * t518 + (-t199 * t446 + t201 * t448 + t437) * t519;
t300 = pkin(2) * t448;
t195 = t300 + (-t472 * t321 + t478) * t320;
t292 = qJ(3) * t444;
t299 = rSges(4,3) * t444;
t196 = -t492 * t446 + t292 + t299;
t173 = -rSges(5,1) * t245 + rSges(5,2) * t246;
t140 = t300 + (-t454 + t482) * t320 - t173;
t175 = -t247 * rSges(5,1) + t248 * rSges(5,2);
t141 = -t517 * t446 + t175 + t292;
t353 = t140 * t322 + t141 * t320;
t481 = (t353 * t319 + t440) * t518 + ((t195 * t322 + t196 * t320) * t319 + t437) * t519;
t552 = t480 - t481;
t586 = t552 * qJD(1);
t308 = Icges(4,5) * t319;
t469 = Icges(4,1) * t321;
t366 = t308 + t469;
t230 = Icges(4,4) * t320 + t366 * t322;
t465 = Icges(3,4) * t319;
t279 = Icges(3,1) * t321 - t465;
t232 = Icges(3,5) * t320 + t279 * t322;
t582 = t230 + t232;
t441 = Icges(5,1) - Icges(5,2);
t389 = t441 * t246;
t339 = t389 + t463;
t548 = -0.2e1 * t217;
t572 = t339 - t548;
t580 = t572 * t247 + t320 * (-Icges(5,5) * t245 + Icges(5,6) * t246);
t162 = rSges(5,3) * t322 - t531;
t117 = -t320 * t162 - t322 * t164;
t568 = -t320 * t173 - t175 * t322;
t578 = t117 * t568;
t167 = -Icges(5,5) * t247 + Icges(5,6) * t248;
t577 = (-t167 * t320 - t247 * t526 + t435 * t248) * t320;
t576 = (t167 * t322 - t245 * t526 + t435 * t246) * t320;
t575 = t568 * t321;
t192 = -rSges(5,1) * t337 - rSges(5,2) * t267;
t317 = t320 ^ 2;
t410 = t317 + t318;
t535 = t410 * t319;
t483 = m(5) * (t192 * t535 + t575);
t295 = Icges(4,5) * t444;
t222 = Icges(4,6) * t320 + Icges(4,3) * t446 + t295;
t272 = Icges(3,5) * t321 - Icges(3,6) * t319;
t224 = Icges(3,3) * t320 + t272 * t322;
t273 = Icges(4,4) * t321 + Icges(4,6) * t319;
t226 = Icges(4,2) * t320 + t273 * t322;
t573 = t222 * t446 + t582 * t444 + (t224 + t226) * t320;
t571 = (-Icges(3,6) + Icges(4,6)) * t321 + (-Icges(4,4) - Icges(3,5)) * t319;
t274 = Icges(3,2) * t321 + t465;
t457 = Icges(4,3) * t321;
t358 = t457 - t308;
t570 = (-t274 - t358) * t322 + t582;
t565 = -t429 * t337 / 0.2e1 + t581;
t497 = t320 / 0.2e1;
t542 = -t322 / 0.2e1;
t540 = t193 * t410;
t558 = (t129 * t322 + t130 * t320) * t192;
t227 = Icges(3,4) * t445 - Icges(3,2) * t448 - Icges(3,6) * t322;
t311 = Icges(3,4) * t321;
t462 = Icges(3,2) * t319;
t228 = Icges(3,6) * t320 + (t311 - t462) * t322;
t205 = t232 * t445;
t375 = t224 * t322 - t205;
t223 = Icges(3,5) * t445 - Icges(3,6) * t448 - Icges(3,3) * t322;
t296 = Icges(3,4) * t448;
t231 = Icges(3,1) * t445 - Icges(3,5) * t322 - t296;
t425 = -t320 * t223 - t231 * t444;
t557 = -t227 * t446 - t228 * t448 - t375 - t425;
t556 = -t228 * t446 + t573;
t450 = (-Icges(4,2) * t322 + t320 * t273) * t322;
t554 = t450 + t573;
t388 = t441 * t245;
t553 = t388 - t458;
t543 = -t320 / 0.2e1;
t464 = Icges(4,5) * t321;
t271 = Icges(4,3) * t319 + t464;
t221 = -Icges(4,6) * t322 + t271 * t320;
t229 = -Icges(4,4) * t322 + t366 * t320;
t537 = (t221 * t319 + t229 * t321) * t320;
t447 = t319 * t321;
t412 = t410 * t447;
t536 = (m(4) / 0.4e1 + m(5) / 0.4e1) * (t412 - t447);
t530 = t571 * t320;
t529 = t571 * t322;
t528 = t570 * t320;
t470 = Icges(3,1) * t319;
t367 = -t311 - t470;
t384 = (t367 * t322 - t228) * t320;
t386 = (-Icges(4,1) * t446 + t222 + t295) * t320;
t527 = t384 + t386;
t376 = -0.2e1 * t218;
t325 = t376 + t553;
t345 = (-t577 - (t325 * t248 + t580) * t322) * t497 + (-t576 - ((t376 - 0.2e1 * t458) * t246 - (-0.2e1 * t389 + t548 - 0.2e1 * t463) * t245) * t322) * t542;
t327 = t407 - t553;
t346 = (t577 - (t327 * t248 - t580) * t322) * t497 + (t576 - (-(-t548 + t406) * t245 + (t360 - 0.2e1 * t388) * t246) * t322) * t542;
t276 = Icges(4,1) * t319 - t464;
t524 = (t279 / 0.2e1 - t274 / 0.2e1 + t308 + t469 / 0.2e1 - t457 / 0.2e1) * t319 + (t311 + t470 / 0.2e1 - t462 / 0.2e1 + t276 / 0.2e1 - t271 / 0.2e1) * t321;
t522 = 0.4e1 * qJD(1);
t521 = 2 * qJD(2);
t512 = m(5) * (t146 * t175 - t148 * t173 - t558);
t511 = m(5) * (t353 * t193 - t558);
t378 = t410 * t192;
t431 = t321 * t540;
t510 = m(5) * (-t575 + (t117 - t378) * t319 + t431);
t438 = -t146 * t445 - t148 * t444;
t284 = pkin(2) * t321 + t455;
t415 = t410 * t284;
t98 = (pkin(3) * t444 + t532) * t322 + (pkin(3) * t445 + pkin(6) * t322 + t162) * t320 + t415;
t508 = m(5) * (t535 * t98 + t438);
t505 = m(5) * (t129 * t140 + t130 * t141);
t504 = m(5) * (t129 * t173 - t130 * t175);
t503 = m(5) * (t117 * t535 + t431);
t502 = m(5) * (-t129 * t448 + t130 * t446);
t479 = rSges(3,1) * t321;
t393 = pkin(1) + t479;
t411 = rSges(3,2) * t448 + t322 * rSges(3,3);
t197 = -t393 * t320 + t316 + t411;
t298 = rSges(3,2) * t446;
t198 = -t298 + t393 * t322 + (rSges(3,3) + pkin(5)) * t320;
t282 = rSges(3,1) * t319 + rSges(3,2) * t321;
t261 = t282 * t320;
t263 = t282 * t322;
t489 = m(3) * (t197 * t261 - t198 * t263);
t285 = rSges(4,1) * t321 + rSges(4,3) * t319;
t126 = t320 * (t285 * t320 - t314) + (t320 * rSges(4,2) + t285 * t322) * t322 + t415;
t430 = -t199 * t445 - t201 * t444;
t487 = m(4) * (t126 * t535 + t430);
t485 = m(4) * (t165 * t195 + t166 * t196);
t484 = m(4) * (-t165 * t448 + t166 * t446);
t449 = t227 * t319;
t340 = (t173 * t322 - t175 * t320) * t319;
t85 = -t340 * m(5) / 0.2e1;
t442 = t85 * qJD(3);
t423 = -t276 * t320 + t221;
t422 = -t367 * t320 + t227;
t421 = -t358 * t320 + t229;
t419 = -Icges(3,2) * t445 + t231 - t296;
t416 = t320 * (qJ(3) * t445 - t300) + t322 * (-pkin(2) * t446 + t292);
t413 = -t284 - t285;
t403 = t15 / 0.2e1 + t47 / 0.2e1;
t402 = t48 / 0.2e1 - t16 / 0.2e1;
t394 = t272 / 0.2e1 + t273 / 0.2e1;
t390 = -t346 - t345;
t387 = t423 * t322;
t385 = t422 * t322;
t383 = t421 * t322;
t381 = t419 * t322;
t374 = t228 * t319 - t223;
t81 = t325 * t267 - (t339 + 0.2e1 * t217) * t337;
t373 = t512 / 0.2e1 + t593 + (-t77 + t81) * t494;
t83 = t327 * t267 + t337 * t572;
t372 = t511 / 0.2e1 + t593 + (t77 + t83) * t493;
t370 = -pkin(3) * t321 + t192 - t284;
t369 = -t222 * t448 + t226 * t322 - t230 * t445;
t352 = t146 * t320 + t148 * t322;
t147 = t370 * t320;
t149 = t370 * t322;
t351 = t147 * t320 + t149 * t322;
t131 = -t450 + t537;
t330 = t402 + (t131 - t537 + t554) * t543 + t556 * t497 + (-t205 + (t224 + t449) * t322 + t425 + t557) * t542;
t329 = t369 * t543 + t403 + (-(-t231 * t321 + t449) * t320 + t131 - t223 * t322) * t542 + (t374 * t322 - t554 + t556) * t322 / 0.2e1 + (t221 * t446 + t229 * t444 + t374 * t320 + t369 + t375 + t557) * t497;
t326 = t403 * t320 + t402 * t322;
t286 = -rSges(3,2) * t319 + t479;
t202 = t413 * t322;
t200 = t413 * t320;
t145 = 0.4e1 * t536;
t139 = t322 * (-rSges(4,1) * t446 + t299) - t281 * t317 + t416;
t110 = -pkin(3) * t535 + t416 - t568;
t95 = -t483 / 0.2e1;
t87 = t503 / 0.2e1;
t86 = t340 * t518;
t66 = t484 + t502;
t53 = -t193 * t378 + t578;
t51 = t510 / 0.2e1;
t38 = t504 + t592;
t32 = t487 + t508;
t22 = t87 + t51 + t483 / 0.2e1;
t21 = t95 + t87 - t510 / 0.2e1;
t20 = t95 + t51 - t503 / 0.2e1;
t17 = t489 + t485 + t505 + t524 - t592;
t9 = t480 + t481;
t7 = m(5) * t53 + t346;
t6 = m(5) * (t352 * t192 + t568 * t98) + t345;
t4 = t329 * t320 + t330 * t322;
t3 = (t81 / 0.4e1 + t590) * t322 - t512 / 0.2e1 + t372 + t587;
t2 = t326 + t372 + t373;
t1 = -t511 / 0.2e1 + (-t83 / 0.4e1 + t590) * t322 + t373 + t587;
t5 = [t17 * qJD(2) + t66 * qJD(3) + t38 * qJD(4), t17 * qJD(1) + t9 * qJD(3) + t2 * qJD(4) + ((t165 * t202 + t166 * t200 - t195 * t201 - t196 * t199) * t519 + (t129 * t149 + t130 * t147 - t140 * t148 - t141 * t146) * t518) * t521 + ((m(3) * (-t197 * t286 - t261 * t282) - t83 / 0.2e1 + t591 + t394 * t322 - t330) * t322 + (m(3) * (-t198 * t286 + t263 * t282) + t394 * t320 - t329 + t588) * t320 + (-t387 / 0.2e1 + t385 / 0.2e1 + t386 / 0.2e1 + t384 / 0.2e1) * t319 + (-t383 / 0.2e1 - t381 / 0.2e1 + t570 * t497) * t321) * qJD(2), qJD(1) * t66 + t9 * qJD(2) - qJD(4) * t85, t38 * qJD(1) + t2 * qJD(2) - t442 + ((m(5) * (t129 * t192 + t173 * t193) + t81 / 0.2e1 + t591 - t402) * t322 + (m(5) * (t130 * t192 - t175 * t193) - t403 + t588) * t320) * qJD(4); t4 * qJD(2) + t552 * qJD(3) + t1 * qJD(4) + (-t489 / 0.4e1 - t485 / 0.4e1 - t505 / 0.4e1) * t522 + (-t524 - t565) * qJD(1), t4 * qJD(1) + (m(5) * (t110 * t98 - t146 * t147 - t148 * t149) + m(4) * (t126 * t139 - t199 * t200 - t201 * t202) + m(3) * ((t320 * (rSges(3,1) * t445 - t411) + t322 * (rSges(3,1) * t444 + t320 * rSges(3,3) - t298)) * (-t320 * t261 - t263 * t322) + t410 * t286 * t282) + t346 + ((-t530 * t320 + (t381 + t383 - t528) * t319 + ((t422 - t423) * t322 + t527) * t321) * t322 + t529 * t317) * t497 + ((-t529 * t322 + (-t387 + t385 + t527) * t321 + ((t419 + t421) * t322 - t528) * t319) * t320 + t530 * t318) * t542) * qJD(2) + t32 * qJD(3) + t6 * qJD(4), t586 + t32 * qJD(2) + t21 * qJD(4) + (-0.4e1 * t536 + 0.2e1 * (t518 + t519) * (-t321 * t535 + t412)) * qJD(3), t1 * qJD(1) + t6 * qJD(2) + t21 * qJD(3) + ((-t53 + (t117 - t98) * t568 + (-t352 - t540) * t192) * m(5) - t345 + t390) * qJD(4); -t552 * qJD(2) + t86 * qJD(4) + (-t502 / 0.4e1 - t484 / 0.4e1) * t522, -t586 + t145 * qJD(3) + t20 * qJD(4) + 0.4e1 * (-t508 / 0.4e1 - t487 / 0.4e1) * qJD(2) + ((-t321 * t110 + t438) * t518 + (-t321 * t139 + t430) * t519 + ((t351 + t98) * t518 + (t200 * t320 + t202 * t322 + t126) * t519) * t319) * t521, t145 * qJD(2), t86 * qJD(1) + t20 * qJD(2) + qJD(4) * t483; (-t504 + t565) * qJD(1) + t3 * qJD(2) + t442 + t326 * qJD(4), t3 * qJD(1) + ((t117 * t110 + t193 * t351) * m(5) - t346 + t390) * qJD(2) + t22 * qJD(3) + t7 * qJD(4), qJD(1) * t85 + t22 * qJD(2), t326 * qJD(1) + t7 * qJD(2) + (m(5) * (t192 * t540 - t578) + t345) * qJD(4);];
Cq = t5;
