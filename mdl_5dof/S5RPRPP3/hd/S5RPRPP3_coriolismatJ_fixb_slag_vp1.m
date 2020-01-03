% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:21
% EndTime: 2019-12-31 18:12:33
% DurationCPUTime: 7.46s
% Computational Cost: add. (10233->393), mult. (12386->522), div. (0->0), fcn. (10933->5), ass. (0->237)
t496 = Icges(5,1) + Icges(4,3);
t495 = Icges(5,4) - Icges(4,5);
t494 = -Icges(5,5) + Icges(4,6);
t297 = pkin(7) + qJ(3);
t291 = sin(t297);
t292 = cos(t297);
t238 = Icges(4,5) * t292 - Icges(4,6) * t291;
t240 = -Icges(5,4) * t292 + Icges(5,5) * t291;
t301 = sin(qJ(1));
t302 = cos(qJ(1));
t493 = (t238 + t240) * t302 + t496 * t301;
t402 = Icges(4,4) * t291;
t244 = Icges(4,1) * t292 - t402;
t178 = Icges(4,5) * t301 + t244 * t302;
t285 = Icges(6,6) * t291;
t396 = Icges(6,3) * t292;
t179 = Icges(6,5) * t301 + (t285 + t396) * t302;
t492 = t179 + t178;
t384 = t292 * t301;
t386 = t291 * t301;
t488 = t302 * t496 + t384 * t495 + t386 * t494;
t245 = pkin(3) * t291 - qJ(4) * t292;
t320 = rSges(6,2) * t292 - rSges(6,3) * t291;
t395 = qJ(5) * t291;
t324 = t245 - t320 + t395;
t126 = t324 * t301;
t128 = t324 * t302;
t409 = rSges(5,2) * t291;
t321 = rSges(5,3) * t292 + t409;
t361 = t245 - t321;
t147 = t361 * t301;
t150 = t361 * t302;
t272 = rSges(5,3) * t386;
t287 = cos(pkin(7)) * pkin(2) + pkin(1);
t387 = t291 * qJ(4);
t418 = pkin(6) + qJ(2);
t110 = -t272 + (rSges(5,1) + t418) * t302 + (-t387 - t287 + (rSges(5,2) - pkin(3)) * t292) * t301;
t288 = t301 * t418;
t383 = t292 * t302;
t343 = t301 * rSges(5,1) - rSges(5,2) * t383;
t405 = rSges(5,3) + qJ(4);
t419 = pkin(3) * t292;
t111 = t288 + (t291 * t405 + t287 + t419) * t302 + t343;
t380 = t110 * t383 + t111 * t384;
t385 = t291 * t302;
t249 = t387 + t419;
t256 = qJ(5) * t384;
t348 = -rSges(6,2) * t386 - rSges(6,3) * t384 - t256;
t426 = rSges(6,1) + pkin(4);
t100 = (t418 + t426) * t302 + (-t249 - t287) * t301 + t348;
t404 = rSges(6,3) + qJ(5);
t352 = pkin(3) + t404;
t407 = rSges(6,2) + qJ(4);
t101 = t288 + t426 * t301 + (t291 * t407 + t292 * t352 + t287) * t302;
t413 = t100 * t383 + t101 * t384;
t457 = m(6) / 0.2e1;
t458 = m(5) / 0.2e1;
t416 = (-t126 * t385 + t128 * t386 + t413) * t457 + (-t147 * t385 + t150 * t386 + t380) * t458;
t281 = pkin(3) * t386;
t144 = t281 + (-t292 * t405 - t409) * t301;
t257 = qJ(4) * t383;
t325 = -pkin(3) * t385 + t257;
t357 = rSges(5,2) * t385 + rSges(5,3) * t383;
t145 = t325 + t357;
t122 = t281 + (t291 * t404 - t292 * t407) * t301;
t279 = rSges(6,2) * t383;
t123 = -t352 * t385 + t257 + t279;
t314 = t122 * t302 + t123 * t301;
t417 = (t291 * t314 + t413) * t457 + ((t144 * t302 + t145 * t301) * t291 + t380) * t458;
t2 = t417 - t416;
t491 = t2 * qJD(1);
t313 = t126 * t301 + t128 * t302;
t248 = rSges(4,1) * t291 + rSges(4,2) * t292;
t298 = t301 ^ 2;
t299 = t302 ^ 2;
t355 = t298 + t299;
t470 = t355 * t248;
t349 = -m(4) * t470 / 0.2e1 - m(6) * t313 / 0.2e1 + (-t147 * t301 - t150 * t302) * t458;
t224 = t248 * t301;
t226 = t248 * t302;
t124 = t224 * t301 + t226 * t302;
t350 = m(4) * t124 / 0.2e1 + (t122 * t301 - t123 * t302) * t457 + (t144 * t301 - t145 * t302) * t458;
t9 = t350 - t349;
t490 = t9 * qJD(1);
t398 = Icges(5,6) * t292;
t233 = Icges(5,3) * t291 - t398;
t181 = Icges(5,5) * t301 + t233 * t302;
t489 = -t178 * t384 - t181 * t386;
t487 = t302 * t493 + t489;
t241 = Icges(4,2) * t292 + t402;
t400 = Icges(6,2) * t292;
t486 = (t285 - t400 - t241) * t302 + t492;
t271 = Icges(4,4) * t386;
t177 = Icges(4,1) * t384 - Icges(4,5) * t302 - t271;
t182 = Icges(5,5) * t302 + Icges(5,6) * t384 - Icges(5,3) * t386;
t485 = -t177 * t383 + t182 * t385 + t301 * t488;
t175 = Icges(4,4) * t384 - Icges(4,2) * t386 - Icges(4,6) * t302;
t263 = Icges(5,6) * t386;
t186 = Icges(5,4) * t302 + Icges(5,2) * t384 - t263;
t484 = t175 * t291 - t186 * t292;
t262 = Icges(6,6) * t383;
t183 = Icges(6,4) * t301 + Icges(6,2) * t385 + t262;
t239 = Icges(6,4) * t291 + Icges(6,5) * t292;
t187 = Icges(6,1) * t301 + t239 * t302;
t483 = (t181 + t183) * t385 + t492 * t383 + (t187 + t493) * t301;
t482 = (Icges(6,4) - t494) * t292 + (-Icges(6,5) + t495) * t291;
t228 = t355 * t292;
t146 = (t228 - t292) * t291;
t411 = m(6) * qJD(5);
t481 = t146 * t411;
t286 = Icges(4,4) * t292;
t401 = Icges(4,2) * t291;
t176 = Icges(4,6) * t301 + (t286 - t401) * t302;
t264 = Icges(5,6) * t385;
t185 = Icges(5,4) * t301 - Icges(5,2) * t383 + t264;
t480 = t176 * t291 + t185 * t292 + t488;
t479 = -t176 * t385 - t185 * t383 + t483;
t389 = (Icges(6,1) * t302 - Icges(6,4) * t386 - Icges(6,5) * t384) * t302;
t478 = -t389 + t483;
t477 = -t175 * t385 - t176 * t386 - t185 * t384 + t186 * t383 - t485 - t487;
t476 = -t301 / 0.2e1;
t428 = t301 / 0.2e1;
t475 = -t302 / 0.2e1;
t473 = t100 * t302 + t101 * t301;
t260 = Icges(6,6) * t386;
t180 = Icges(6,5) * t302 - Icges(6,3) * t384 - t260;
t261 = Icges(6,6) * t384;
t184 = Icges(6,4) * t302 - Icges(6,2) * t386 - t261;
t472 = (t180 * t292 + t184 * t291) * t301;
t353 = m(5) / 0.4e1 + m(6) / 0.4e1;
t471 = t353 * t146;
t331 = (Icges(5,3) * t383 + t185 + t264) * t301;
t468 = qJD(3) * t302;
t467 = t482 * t301;
t466 = t482 * t302;
t465 = -t301 * t486 + t331;
t333 = (-Icges(6,3) * t385 + t183 + t262) * t301;
t315 = Icges(5,2) * t291 + t398;
t335 = (-t302 * t315 + t181) * t301;
t403 = Icges(4,1) * t291;
t319 = -t286 - t403;
t340 = (t302 * t319 - t176) * t301;
t464 = t333 + t335 + t340;
t463 = -t291 * (t244 / 0.2e1 - t241 / 0.2e1 - Icges(5,6) * t291 + t285 - t400 / 0.2e1 + t396 / 0.2e1 + (Icges(5,2) / 0.2e1 - Icges(5,3) / 0.2e1) * t292) - t292 * (t286 + t403 / 0.2e1 - t401 / 0.2e1 + t315 / 0.2e1 - t233 / 0.2e1 - Icges(6,6) * t292 + (-Icges(6,2) / 0.2e1 + Icges(6,3) / 0.2e1) * t291);
t227 = t355 * t291;
t462 = 0.2e1 * t227;
t461 = 0.4e1 * qJD(1);
t460 = 0.2e1 * qJD(3);
t410 = rSges(4,1) * t292;
t344 = t287 + t410;
t358 = rSges(4,2) * t386 + rSges(4,3) * t302;
t141 = -t301 * t344 + t302 * t418 + t358;
t342 = -rSges(4,2) * t385 + t301 * rSges(4,3);
t142 = t302 * t344 + t288 + t342;
t456 = m(4) * (t141 * t224 - t142 * t226);
t455 = m(4) * (t141 * t302 + t142 * t301);
t378 = -t147 * t384 - t150 * t383;
t363 = t355 * t249;
t68 = -t301 * (rSges(5,1) * t302 + rSges(5,2) * t384 - t272) + t302 * (rSges(5,3) * t385 + t343) + t363;
t451 = m(5) * (t227 * t68 + t378);
t449 = m(5) * (t110 * t144 + t111 * t145);
t448 = m(5) * (-t110 * t386 + t111 * t385);
t447 = m(5) * (t110 * t302 + t111 * t301);
t307 = t473 * t291;
t442 = m(6) * (t292 * t314 - t307);
t440 = m(6) * (-t126 * t383 + t128 * t384 - t307);
t379 = -t126 * t384 - t128 * t383;
t408 = rSges(6,2) * t291;
t58 = -t348 * t301 + t363 + (t292 * t404 + t408) * t299;
t438 = m(6) * (t227 * t58 + t379);
t435 = m(6) * (t100 * t122 + t101 * t123);
t434 = m(6) * (-t100 * t386 + t101 * t385);
t433 = m(6) * t473;
t425 = m(3) * t355 * (rSges(3,3) + qJ(2));
t289 = t291 ^ 2;
t290 = t292 ^ 2;
t356 = t355 * t290;
t422 = m(6) * (-t290 + (0.1e1 - t355) * t289 + t356);
t421 = m(6) * (-t228 * t292 - t289 * t355);
t420 = m(6) * (t227 * t291 + t356);
t412 = m(6) * qJD(3);
t86 = (t458 + t457) * t462;
t381 = t86 * qJD(1);
t373 = -t301 * t319 + t175;
t372 = -Icges(4,2) * t384 + t177 - t271;
t369 = Icges(6,2) * t384 + t180 - t260;
t368 = t301 * t315 + t182;
t367 = Icges(6,3) * t386 + t184 - t261;
t365 = -Icges(5,3) * t384 + t186 - t263;
t364 = t301 * (qJ(4) * t384 - t281) + t302 * t325;
t360 = -rSges(6,3) * t292 - t249 - t408;
t359 = rSges(5,2) * t292 - rSges(5,3) * t291 - t249;
t135 = m(6) * t228;
t354 = t135 * qJD(1);
t53 = -t100 * t384 + t101 * t383;
t351 = m(6) * t53 * qJD(1);
t341 = t373 * t302;
t339 = t372 * t302;
t336 = t369 * t302;
t334 = t368 * t302;
t332 = t367 * t302;
t330 = t365 * t302;
t323 = -t179 * t384 - t183 * t386 + t187 * t302;
t322 = t239 / 0.2e1 + t240 / 0.2e1 + t238 / 0.2e1;
t127 = t301 * t360 - t256;
t129 = (-qJ(5) * t292 + t360) * t302;
t308 = t127 * t301 + t129 * t302 + t58;
t81 = t389 - t472;
t306 = (t81 + t472 + t478) * t476 + t479 * t428 + ((t484 + t493) * t302 + t477 + t485 + t489) * t475;
t305 = t323 * t476 + (t81 + t488 * t302 + (t177 * t292 - t182 * t291 - t484) * t301) * t475 + (t480 * t302 - t478 + t479) * t302 / 0.2e1 + (-t180 * t383 - t184 * t385 + t301 * t480 + t323 + t477 + t487) * t428;
t252 = -rSges(4,2) * t291 + t410;
t151 = t359 * t302;
t148 = t359 * t301;
t131 = t420 / 0.2e1;
t130 = t421 / 0.2e1;
t120 = t422 / 0.2e1;
t89 = 0.4e1 * t471;
t87 = t353 * t462 - (m(5) + m(6)) * t227 / 0.2e1;
t84 = t298 * t321 + t302 * t357 + t364;
t63 = t302 * (-rSges(6,3) * t385 + t279) + t320 * t298 - t355 * t395 + t364;
t57 = t131 + t120 - t421 / 0.2e1;
t56 = t130 + t131 - t422 / 0.2e1;
t55 = t130 + t120 - t420 / 0.2e1;
t31 = t58 * t228 + t291 * t313;
t27 = t440 / 0.2e1;
t23 = t442 / 0.2e1;
t22 = t434 + t448;
t20 = t425 + t433 + t447 + t455;
t18 = t438 + t451;
t10 = t349 + t350;
t8 = t27 - t442 / 0.2e1;
t7 = t27 + t23;
t6 = t23 - t440 / 0.2e1;
t5 = t435 + t449 + t456 - t463;
t3 = t416 + t417;
t1 = t301 * t305 + t302 * t306;
t4 = [qJD(2) * t20 + qJD(3) * t5 + qJD(4) * t22 + t411 * t53, qJD(1) * t20 + qJD(3) * t10 + qJD(4) * t87, t5 * qJD(1) + t10 * qJD(2) + t3 * qJD(4) + t7 * qJD(5) + ((t110 * t151 + t111 * t148 - t144 * t150 - t145 * t147) * t458 + (t100 * t129 + t101 * t127 - t122 * t128 - t123 * t126) * t457) * t460 + (m(4) * (-t141 * t252 - t224 * t248) + t322 * t302 - t306) * t468 + ((m(4) * (-t142 * t252 + t226 * t248) + t322 * t301 - t305) * t301 + (t333 / 0.2e1 + t335 / 0.2e1 + t340 / 0.2e1 + t332 / 0.2e1 + t334 / 0.2e1 + t341 / 0.2e1) * t291 + (-t331 / 0.2e1 + t336 / 0.2e1 - t330 / 0.2e1 - t339 / 0.2e1 + t486 * t428) * t292) * qJD(3), qJD(1) * t22 + qJD(2) * t87 + qJD(3) * t3, qJD(3) * t7 + t351; t9 * qJD(3) - t86 * qJD(4) - t135 * qJD(5) + (-t425 / 0.4e1 - t455 / 0.4e1 - t447 / 0.4e1 - t433 / 0.4e1) * t461, 0, t490 + ((-t148 * t302 + t151 * t301) * t458 + (-t127 * t302 + t129 * t301) * t457) * t460, -t381, -t354; -t9 * qJD(2) + t1 * qJD(3) - t2 * qJD(4) + t8 * qJD(5) + (-t456 / 0.4e1 - t435 / 0.4e1 - t449 / 0.4e1) * t461 + t463 * qJD(1), -t490, t1 * qJD(1) + t18 * qJD(4) + t31 * t411 - ((-t466 * t302 + (t334 + t341 + t332 + t464) * t292 + ((t365 - t369 + t372) * t302 + t465) * t291) * t301 + t467 * t299) * t468 / 0.2e1 + (m(4) * (t252 * t470 - (t301 * (rSges(4,1) * t384 - t358) + t302 * (rSges(4,1) * t383 + t342)) * t124) + m(5) * (-t147 * t148 - t150 * t151 + t68 * t84) + m(6) * (-t126 * t127 - t128 * t129 + t58 * t63) + ((-t467 * t301 + (-t336 + t339 + t330 + t465) * t291 + ((t367 + t368 + t373) * t302 + t464) * t292) * t302 + t466 * t298) * t428) * qJD(3), qJD(3) * t18 - 0.4e1 * qJD(4) * t471 + qJD(5) * t56 - t491, t8 * qJD(1) + t56 * qJD(4) + t31 * t412 + t481; t86 * qJD(2) + t2 * qJD(3) + (-t448 / 0.4e1 - t434 / 0.4e1) * t461, t381, t491 + t89 * qJD(4) + t55 * qJD(5) + 0.4e1 * (-t451 / 0.4e1 - t438 / 0.4e1) * qJD(3) + ((-t292 * t84 + t378) * t458 + (-t292 * t63 + t379) * t457 + ((t148 * t301 + t151 * t302 + t68) * t458 + t308 * t457) * t291) * t460, t89 * qJD(3), t55 * qJD(3); t135 * qJD(2) + t6 * qJD(3) - t351, t354, t6 * qJD(1) + (t308 * t292 + (t313 + t63) * t291 - t31) * t412 + t57 * qJD(4) - t481, t57 * qJD(3), -t146 * t412;];
Cq = t4;
