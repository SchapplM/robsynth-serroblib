% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR12_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR12_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:58
% EndTime: 2019-12-31 18:07:13
% DurationCPUTime: 11.06s
% Computational Cost: add. (28902->548), mult. (38740->817), div. (0->0), fcn. (41893->8), ass. (0->324)
t331 = pkin(8) + qJ(4);
t323 = sin(t331);
t338 = cos(qJ(5));
t339 = cos(qJ(1));
t422 = t338 * t339;
t336 = sin(qJ(5));
t337 = sin(qJ(1));
t426 = t337 * t336;
t289 = -t323 * t426 + t422;
t425 = t337 * t338;
t428 = t336 * t339;
t290 = t323 * t425 + t428;
t324 = cos(t331);
t432 = t324 * t337;
t200 = Icges(6,5) * t290 + Icges(6,6) * t289 - Icges(6,3) * t432;
t447 = Icges(6,4) * t290;
t203 = Icges(6,2) * t289 - Icges(6,6) * t432 + t447;
t278 = Icges(6,4) * t289;
t206 = Icges(6,1) * t290 - Icges(6,5) * t432 + t278;
t291 = t323 * t428 + t425;
t292 = t323 * t422 - t426;
t431 = t324 * t339;
t111 = t200 * t431 + t291 * t203 - t292 * t206;
t542 = t111 * t339;
t541 = t337 * t111;
t202 = -Icges(6,5) * t292 + Icges(6,6) * t291 + Icges(6,3) * t431;
t438 = t202 * t323;
t280 = Icges(6,4) * t292;
t205 = Icges(6,2) * t291 + Icges(6,6) * t431 - t280;
t279 = Icges(6,4) * t291;
t207 = Icges(6,1) * t292 - Icges(6,5) * t431 - t279;
t536 = t205 * t336 + t207 * t338;
t122 = t536 * t324 - t438;
t299 = rSges(5,1) * t324 - rSges(5,2) * t323;
t332 = t337 ^ 2;
t333 = t339 ^ 2;
t318 = t332 + t333;
t504 = -m(6) / 0.2e1;
t377 = rSges(6,1) * t338 - rSges(6,2) * t336;
t260 = rSges(6,3) * t323 + t377 * t324;
t301 = pkin(4) * t324 + pkin(7) * t323;
t403 = t260 + t301;
t215 = t403 * t337;
t217 = t403 * t339;
t517 = t215 * t337 + t217 * t339;
t525 = -m(5) / 0.2e1;
t415 = t318 * t299 * t525 + t517 * t504;
t194 = ((rSges(6,3) + pkin(7)) * t323 + (pkin(4) + t377) * t324) * t339;
t281 = t299 * t337;
t283 = t299 * t339;
t219 = t337 * t281 + t283 * t339;
t434 = t323 * t337;
t242 = rSges(6,3) * t434 + (rSges(6,1) * t425 - rSges(6,2) * t426) * t324;
t390 = -pkin(4) * t432 - pkin(7) * t434 - t242;
t530 = t390 * t337;
t534 = m(6) / 0.2e1;
t535 = m(5) / 0.2e1;
t416 = (t194 * t339 - t530) * t534 + t219 * t535;
t49 = t416 - t415;
t540 = t49 * qJD(1);
t364 = -t291 * t205 - t292 * t207;
t412 = t289 * t203 + t290 * t206;
t539 = t364 + t412 + (-t200 * t337 - t202 * t339) * t324;
t515 = t292 * rSges(6,1) - t291 * rSges(6,2);
t213 = rSges(6,3) * t431 - t515;
t538 = -t213 * t323 + t260 * t431;
t110 = -t202 * t432 + t289 * t205 - t207 * t290;
t369 = Icges(6,5) * t338 - Icges(6,6) * t336;
t254 = Icges(6,3) * t323 + t324 * t369;
t445 = Icges(6,4) * t338;
t372 = -Icges(6,2) * t336 + t445;
t256 = Icges(6,6) * t323 + t324 * t372;
t446 = Icges(6,4) * t336;
t374 = Icges(6,1) * t338 - t446;
t258 = Icges(6,5) * t323 + t324 * t374;
t135 = t254 * t431 + t291 * t256 - t292 * t258;
t532 = t135 * t323;
t449 = Icges(5,4) * t323;
t373 = Icges(5,2) * t324 + t449;
t263 = Icges(5,6) * t339 + t337 * t373;
t311 = Icges(5,4) * t432;
t265 = Icges(5,1) * t434 + Icges(5,5) * t339 + t311;
t360 = -t324 * t263 - t323 * t265;
t529 = t339 * t360;
t317 = pkin(7) * t431;
t325 = t339 * qJ(2);
t461 = -pkin(6) - qJ(3);
t334 = sin(pkin(8));
t464 = pkin(3) * t334;
t357 = t339 * t464 + t325 + (-pkin(1) + t461) * t337;
t457 = rSges(6,3) * t324;
t463 = pkin(4) * t323;
t528 = (-t457 + t463) * t339 - t317 + t357 + t515;
t264 = -Icges(5,6) * t337 + t339 * t373;
t448 = Icges(5,4) * t324;
t375 = Icges(5,1) * t323 + t448;
t266 = -Icges(5,5) * t337 + t339 * t375;
t273 = -Icges(5,2) * t434 + t311;
t295 = -Icges(5,2) * t323 + t448;
t274 = t295 * t339;
t297 = Icges(5,1) * t324 - t449;
t276 = t297 * t337;
t277 = t297 * t339;
t527 = -t323 * (t337 * (t264 - t277) - t339 * (t263 - t276)) + t324 * (t337 * (t266 + t274) - t339 * (t265 + t273));
t491 = -t337 / 0.2e1;
t490 = t337 / 0.2e1;
t488 = t339 / 0.2e1;
t524 = t110 * t337;
t523 = t110 * t339;
t522 = t260 * t339;
t433 = t323 * t339;
t253 = Icges(6,3) * t324 - t323 * t369;
t423 = t338 * t258;
t429 = t336 * t256;
t521 = t323 * (t423 / 0.2e1 - t429 / 0.2e1 + t297 / 0.2e1 - t373 / 0.2e1 - t253 / 0.2e1);
t520 = (t324 * t264 + t323 * t266) * t339;
t321 = t339 * t461;
t383 = qJ(2) + t464;
t465 = pkin(1) * t339;
t211 = t290 * rSges(6,1) + t289 * rSges(6,2) - rSges(6,3) * t432;
t516 = -pkin(7) * t432 + t211;
t173 = t465 - t321 + (t383 + t463) * t337 + t516;
t518 = t339 * t173 - t337 * t528;
t168 = t211 * t323 + t260 * t432;
t514 = -t168 * t339 + t337 * t538;
t417 = (t337 * t194 + t339 * t390) * t534 + (-t281 * t339 + t337 * t283) * t535;
t361 = t423 - t429;
t348 = t253 - t361;
t436 = t254 * t323;
t512 = t324 * t348 - t436;
t349 = -t254 * t339 + t536;
t511 = t324 * t349 - t438;
t365 = -t203 * t336 + t206 * t338;
t350 = t254 * t337 - t365;
t440 = t200 * t323;
t510 = t324 * t350 - t440;
t272 = (-Icges(6,2) * t338 - t446) * t324;
t275 = (-Icges(6,1) * t336 - t445) * t324;
t509 = -(t258 / 0.2e1 + t272 / 0.2e1) * t336 + (t275 / 0.2e1 - t256 / 0.2e1) * t338;
t506 = 0.2e1 * t318;
t505 = 0.4e1 * qJD(1);
t133 = -t254 * t432 + t256 * t289 + t258 * t290;
t132 = t133 * t323;
t109 = -t200 * t432 + t412;
t368 = -t109 * t337 + t523;
t42 = t324 * t368 + t132;
t503 = -t42 / 0.2e1;
t225 = Icges(6,5) * t291 + Icges(6,6) * t292;
t408 = Icges(6,2) * t292 - t207 + t279;
t410 = -Icges(6,1) * t291 + t205 - t280;
t95 = t225 * t323 + (-t408 * t336 - t410 * t338) * t324;
t502 = t95 / 0.2e1;
t259 = -t377 * t323 + t457;
t124 = (t259 * t337 + t211) * t324 + (-t260 * t337 + t242) * t323;
t125 = (t259 * t339 - t213) * t324;
t500 = m(6) * (t124 * t173 + t125 * t528 - t168 * t390 + t194 * t538);
t230 = rSges(6,1) * t289 - rSges(6,2) * t290;
t231 = rSges(6,1) * t291 + rSges(6,2) * t292;
t282 = (-rSges(6,1) * t336 - rSges(6,2) * t338) * t324;
t498 = m(6) * (-t215 * t231 - t217 * t230 - t518 * t282);
t495 = m(6) * (-t173 * t390 + t194 * t528);
t494 = t318 / 0.2e1;
t493 = t323 / 0.2e1;
t492 = t324 / 0.2e1;
t489 = t337 / 0.4e1;
t487 = t339 / 0.4e1;
t486 = m(3) * ((rSges(3,3) * t339 - t337 * pkin(1) + t325) * t339 + (t465 + (rSges(3,3) + qJ(2)) * t337) * t337);
t380 = rSges(4,1) * t334 + rSges(4,2) * cos(pkin(8));
t394 = rSges(4,3) + pkin(1) + qJ(3);
t244 = -t394 * t337 + t380 * t339 + t325;
t245 = t394 * t339 + (qJ(2) + t380) * t337;
t485 = m(4) * (-t244 * t337 + t339 * t245);
t484 = m(4) * (t244 * t339 + t337 * t245);
t379 = rSges(5,1) * t323 + rSges(5,2) * t324;
t345 = -t337 * rSges(5,3) + t379 * t339;
t220 = t345 + t357;
t221 = -t321 + (rSges(5,3) + pkin(1)) * t339 + (t379 + t383) * t337;
t483 = m(5) * (t220 * t283 + t221 * t281);
t482 = m(5) * (-t220 * t337 + t339 * t221);
t481 = m(5) * (t220 * t339 + t337 * t221);
t477 = m(6) * t514;
t476 = m(6) * (t168 * t337 + t339 * t538);
t475 = m(6) * t518;
t474 = m(6) * (t337 * t173 + t339 * t528);
t470 = m(6) * (t215 * t339 - t217 * t337);
t174 = -t230 * t339 - t337 * t231;
t467 = m(6) * t174;
t176 = t337 * t230 - t231 * t339;
t466 = m(6) * t176;
t462 = qJD(4) / 0.2e1;
t459 = m(6) * qJD(5);
t455 = t337 * t42;
t112 = t202 * t431 - t364;
t367 = t112 * t339 - t541;
t43 = t324 * t367 + t532;
t454 = t339 * t43;
t453 = -t109 + t539;
t450 = t112 + t539;
t269 = (-Icges(6,5) * t336 - Icges(6,6) * t338) * t324;
t435 = t323 * t269;
t255 = Icges(6,6) * t324 - t323 * t372;
t430 = t336 * t255;
t257 = Icges(6,5) * t324 - t323 * t374;
t424 = t338 * t257;
t346 = m(6) * (-t124 * t339 + t125 * t337);
t356 = m(6) * t282 * t494;
t58 = -t346 / 0.2e1 + t356;
t421 = t58 * qJD(2);
t77 = t124 * t337 + t125 * t339;
t420 = t77 * qJD(3);
t419 = t77 * qJD(5);
t411 = -Icges(6,1) * t289 + t203 + t447;
t409 = -Icges(6,2) * t290 + t206 + t278;
t406 = t256 - t275;
t405 = t258 + t272;
t404 = pkin(7) * t324 + t259 - t463;
t397 = qJD(1) * t324;
t396 = qJD(5) * t324;
t182 = (t504 + t525 - m(4) / 0.2e1) * t506;
t395 = t182 * qJD(1);
t13 = -t532 + (t453 * t339 + t541) * t324;
t392 = -t13 / 0.2e1 - t43 / 0.2e1;
t12 = t132 + (-t450 * t337 + t523) * t324;
t391 = t503 + t12 / 0.2e1;
t370 = Icges(5,5) * t323 + Icges(5,6) * t324;
t156 = t339 * (Icges(5,3) * t339 + t337 * t370) + t263 * t432 + t265 * t434;
t262 = -Icges(5,3) * t337 + t339 * t370;
t157 = -t339 * t262 - t264 * t432 - t266 * t434;
t387 = -t432 / 0.4e1;
t386 = t431 / 0.4e1;
t382 = t318 * t379;
t105 = -t269 * t432 + t405 * t289 - t406 * t290;
t106 = t269 * t431 + t405 * t291 + t406 * t292;
t224 = Icges(6,5) * t289 - Icges(6,6) * t290;
t94 = t224 * t323 + (-t409 * t336 - t411 * t338) * t324;
t381 = t498 / 0.2e1 + (t106 + t95) * t489 + (t105 + t94) * t487;
t371 = Icges(5,5) * t324 - Icges(5,6) * t323;
t121 = t324 * t365 + t440;
t366 = -t121 * t337 - t122 * t339;
t362 = t211 * t339 + t213 * t337;
t358 = -m(6) * (t173 * t230 - t231 * t528) - t435 / 0.2e1;
t159 = -t337 * t262 + t520;
t17 = t450 * t339 + t524;
t55 = t109 * t339 + t524;
t355 = -t156 * t339 / 0.2e1 + t157 * t491 + (t157 - t529) * t490 + (t159 - t520 + (t262 + t360) * t337 + t156) * t488 - t55 / 0.2e1 + t17 / 0.2e1;
t18 = t453 * t337 - t542;
t56 = t112 * t337 + t542;
t354 = t332 * t262 / 0.2e1 + t159 * t490 + t18 / 0.2e1 + t56 / 0.2e1 + (t157 + (t262 - t360) * t339 + t529) * t488;
t84 = -t224 * t432 + t409 * t289 - t411 * t290;
t85 = -t225 * t432 + t408 * t289 - t410 * t290;
t35 = t85 * t337 + t339 * t84;
t86 = t224 * t431 + t409 * t291 + t411 * t292;
t87 = t225 * t431 + t408 * t291 + t410 * t292;
t36 = t87 * t337 + t339 * t86;
t353 = t35 * t488 + t36 * t490;
t344 = t12 * t489 + t13 * t487 + t17 * t386 - t455 / 0.4e1 + t454 / 0.4e1 - t55 * t431 / 0.4e1 + (t18 + t56) * t387;
t101 = (t254 + t424 - t430) * t324 + t348 * t323;
t148 = t324 * t361 + t436;
t237 = t256 * t337;
t239 = t258 * t337;
t88 = (-t237 * t336 + t239 * t338 + t200) * t324 + t350 * t323;
t238 = t256 * t339;
t240 = t258 * t339;
t89 = (t238 * t336 - t240 * t338 + t202) * t324 + t349 * t323;
t98 = t255 * t289 + t257 * t290 - t512 * t337;
t99 = t291 * t255 - t292 * t257 + t512 * t339;
t343 = t101 * t493 + t148 * t492 + t500 / 0.2e1 + (t121 + t133) * t434 / 0.4e1 - (-t122 + t135) * t433 / 0.4e1 + (t88 + t98) * t387 + (t89 + t99) * t386;
t342 = t424 / 0.2e1 - t430 / 0.2e1 + t254 / 0.2e1 - t375 / 0.2e1 - t295 / 0.2e1;
t271 = t371 * t339;
t270 = t337 * t371;
t216 = t404 * t339;
t214 = t404 * t337;
t181 = (-m(6) / 0.4e1 - m(5) / 0.4e1 - m(4) / 0.4e1) * t506 + (m(4) + m(5) + m(6)) * t494;
t180 = -t323 * t231 + t282 * t431;
t179 = t230 * t323 + t282 * t432;
t167 = t466 / 0.2e1;
t166 = t467 / 0.2e1;
t155 = t174 * t324;
t145 = t470 / 0.2e1;
t139 = t362 * t324;
t136 = t530 + (-t301 * t339 - t522) * t339;
t128 = (-pkin(4) * t433 + t213 + t317) * t339 + (-pkin(4) * t434 - t516) * t337;
t126 = (t435 + (-t405 * t336 - t406 * t338) * t324) * t323;
t114 = t476 / 0.2e1;
t113 = -t477 / 0.2e1;
t102 = (-t242 * t339 + t337 * t522) * t324 + t362 * t323;
t81 = -t291 * t238 + t292 * t240 + t511 * t339;
t80 = t291 * t237 - t292 * t239 + t510 * t339;
t79 = -t238 * t289 - t240 * t290 - t511 * t337;
t78 = t237 * t289 + t239 * t290 - t510 * t337;
t75 = m(6) * t77 * t462;
t69 = t475 + t482 + t485;
t66 = -t128 * t176 + t517 * t282;
t65 = t145 - t417;
t64 = -t470 / 0.2e1 + t417;
t63 = t145 + t417;
t57 = t346 / 0.2e1 + t356;
t48 = t415 + t416;
t47 = t324 * t509 - t358;
t46 = t474 + t481 + t484 + t486;
t45 = t148 * t323 + t324 * t366;
t34 = t114 - t467 / 0.2e1;
t33 = t113 - t466 / 0.2e1;
t32 = t167 + t113;
t31 = t167 + t477 / 0.2e1;
t30 = t166 + t114;
t29 = t166 - t476 / 0.2e1;
t28 = t81 * t337 + t339 * t80;
t27 = t79 * t337 + t339 * t78;
t26 = t324 * t342 + t483 + t495 - t521;
t23 = t106 * t323 + (-t337 * t86 + t339 * t87) * t324;
t22 = t105 * t323 + (-t337 * t84 + t339 * t85) * t324;
t21 = -t102 * t139 + t124 * t168 + t125 * t538;
t14 = (-t88 * t337 + t89 * t339 + t148) * t324 + (t101 - t366) * t323;
t9 = (-t337 * t80 + t339 * t81 + t135) * t324 + (-t367 + t99) * t323;
t8 = (-t337 * t78 + t339 * t79 + t133) * t324 + (-t368 + t98) * t323;
t7 = m(6) * t66 + t353;
t6 = (t392 * t337 + t391 * t339) * t324;
t5 = t337 * t355 + t339 * t354;
t4 = m(6) * t21 + (t8 * t491 + t9 * t488 + t45 / 0.2e1) * t324 + (t455 / 0.2e1 - t454 / 0.2e1 + t14 / 0.2e1) * t323;
t3 = t343 + ((-t17 / 0.4e1 + t55 / 0.4e1) * t339 + (t18 / 0.4e1 + t56 / 0.4e1) * t337) * t324 + (t42 / 0.4e1 - t12 / 0.4e1) * t337 + (-t43 / 0.4e1 - t13 / 0.4e1) * t339 + t381;
t2 = t344 + t343 - t498 / 0.2e1 + (-t106 / 0.4e1 - t95 / 0.4e1) * t337 + (-t105 / 0.4e1 - t94 / 0.4e1) * t339;
t1 = t344 + (-t148 / 0.2e1 + (-t99 / 0.4e1 - t89 / 0.4e1) * t339 + (t98 / 0.4e1 + t88 / 0.4e1) * t337) * t324 - t500 / 0.2e1 + (-t101 / 0.2e1 + (t135 / 0.4e1 - t122 / 0.4e1) * t339 + (-t133 / 0.4e1 - t121 / 0.4e1) * t337) * t323 + t381;
t10 = [t46 * qJD(2) + t69 * qJD(3) + t26 * qJD(4) + t47 * qJD(5), qJD(1) * t46 + qJD(3) * t181 + qJD(4) * t63 + qJD(5) * t30, qJD(1) * t69 + qJD(2) * t181 + qJD(4) * t48 + qJD(5) * t32, t26 * qJD(1) + t63 * qJD(2) + t48 * qJD(3) + t3 * qJD(5) + (m(6) * (-t173 * t216 + t194 * t215 + t214 * t528 + t217 * t390) + (t98 / 0.2e1 + m(5) * (t221 * t379 - t281 * t299) + t88 / 0.2e1 - t370 * t488 + (-t263 / 0.2e1 + t276 / 0.2e1) * t324 + (-t265 / 0.2e1 - t273 / 0.2e1) * t323 - t354) * t339 + (t89 / 0.2e1 + t99 / 0.2e1 + m(5) * (-t220 * t379 + t283 * t299) - t370 * t490 + (t264 / 0.2e1 - t277 / 0.2e1) * t324 + (t266 / 0.2e1 + t274 / 0.2e1) * t323 - t355) * t337) * qJD(4), t47 * qJD(1) + t30 * qJD(2) + t32 * qJD(3) + t3 * qJD(4) + (t126 + m(6) * (t168 * t230 + t173 * t179 + t180 * t528 - t231 * t538)) * qJD(5) + ((t502 + t106 / 0.2e1 - t391) * t339 + (-t94 / 0.2e1 - t105 / 0.2e1 - t392) * t337) * t396; t182 * qJD(3) + t64 * qJD(4) + t29 * qJD(5) + (-t474 / 0.4e1 - t481 / 0.4e1 - t484 / 0.4e1 - t486 / 0.4e1) * t505, 0, t395, t64 * qJD(1) + 0.2e1 * ((t214 * t337 + t216 * t339) * t534 + t382 * t525) * qJD(4) + t57 * qJD(5), t29 * qJD(1) + t57 * qJD(4) + (-t179 * t339 + t180 * t337) * t459; -t182 * qJD(2) + t49 * qJD(4) + t31 * qJD(5) + (-t475 / 0.4e1 - t482 / 0.4e1 - t485 / 0.4e1) * t505, -t395, 0, t540 + 0.2e1 * ((t214 * t339 - t216 * t337) * t462 + t419 / 0.4e1) * m(6), t31 * qJD(1) + t75 + (t179 * t337 + t180 * t339) * t459; t65 * qJD(2) - t49 * qJD(3) + t5 * qJD(4) + t1 * qJD(5) + (-t495 / 0.4e1 - t483 / 0.4e1) * t505 - t342 * t397 + t521 * qJD(1), qJD(1) * t65 + qJD(5) * t58, t419 * t504 - t540, t5 * qJD(1) + (m(5) * (-(-t339 * t345 + (-t339 * rSges(5,3) - t379 * t337) * t337) * t219 - t299 * t382) + m(6) * (t128 * t136 + t214 * t215 + t216 * t217) + (-t332 * t271 + (t337 * t270 + t527) * t339 + t28) * t490 + (t333 * t270 + (-t339 * t271 - t527) * t337 + t27) * t488) * qJD(4) + t7 * qJD(5), t1 * qJD(1) + t421 + t7 * qJD(4) + t420 * t504 + (-t45 / 0.2e1 + (t36 / 0.2e1 - t9 / 0.2e1) * t339 + (-t35 / 0.2e1 + t8 / 0.2e1) * t337) * t396 + (t22 * t488 + t23 * t490 + (-t14 / 0.2e1 + (t94 / 0.2e1 + t43 / 0.2e1) * t339 + (t502 + t503) * t337) * t323 + (t155 * t128 + t139 * t176 - t179 * t217 + t180 * t215 + t514 * t282 - t21) * m(6)) * qJD(5); t358 * qJD(1) + t34 * qJD(2) + t33 * qJD(3) + t2 * qJD(4) + t6 * qJD(5) - t509 * t397, qJD(1) * t34 - qJD(4) * t58, qJD(1) * t33 + t75, t2 * qJD(1) - t421 + (t8 * t488 + t55 * t434 / 0.2e1 - t27 * t432 / 0.2e1 + t9 * t490 - t56 * t433 / 0.2e1 + t28 * t431 / 0.2e1 + (t121 * t339 - t122 * t337) * t492 + (t89 * t337 + t88 * t339) * t493 - t353) * qJD(4) + t4 * qJD(5) + (t420 / 0.2e1 + (t102 * t128 - t124 * t217 + t125 * t215 - t136 * t139 - t168 * t216 + t214 * t538 - t66) * qJD(4)) * m(6), t6 * qJD(1) + t4 * qJD(4) + (m(6) * (-t139 * t155 + t168 * t179 + t180 * t538) + t126 * t493 + (t22 * t491 + t23 * t488 + (-t94 * t337 + t95 * t339) * t493) * t324) * qJD(5);];
Cq = t10;
