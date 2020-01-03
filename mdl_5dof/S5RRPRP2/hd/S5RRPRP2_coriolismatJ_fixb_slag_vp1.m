% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:27
% EndTime: 2019-12-31 19:49:34
% DurationCPUTime: 4.36s
% Computational Cost: add. (22208->301), mult. (15156->386), div. (0->0), fcn. (13513->8), ass. (0->202)
t299 = qJ(1) + qJ(2);
t294 = pkin(8) + t299;
t292 = cos(t294);
t300 = sin(qJ(4));
t369 = t292 * t300;
t291 = sin(t294);
t302 = cos(qJ(4));
t407 = rSges(6,1) + pkin(4);
t331 = t407 * t300;
t390 = rSges(6,3) + qJ(5);
t470 = -t390 * t302 + t331;
t474 = t470 * t291;
t475 = t369 * t474;
t298 = Icges(5,4) * t302;
t268 = Icges(5,1) * t300 + t298;
t386 = Icges(6,5) * t302;
t473 = Icges(6,1) * t300 + t268 - t386;
t297 = Icges(6,5) * t300;
t320 = Icges(6,3) * t302 - t297;
t387 = Icges(5,4) * t300;
t471 = Icges(5,2) * t302 + t320 + t387;
t265 = -Icges(5,2) * t300 + t298;
t472 = t265 + t473;
t288 = t291 ^ 2;
t289 = t292 ^ 2;
t335 = t288 + t289;
t269 = Icges(5,1) * t302 - t387;
t446 = Icges(6,1) * t302 + t297;
t469 = t446 + t269;
t447 = t390 * t300 + t407 * t302;
t295 = sin(t299);
t399 = pkin(2) * t295;
t329 = t292 * pkin(7) - t399;
t442 = pkin(3) + t447;
t145 = t292 * rSges(6,2) - t442 * t291 + t329;
t401 = sin(qJ(1)) * pkin(1);
t135 = t145 - t401;
t467 = t135 - t145;
t465 = (-Icges(5,6) + Icges(6,6)) * t302 + (-Icges(6,4) - Icges(5,5)) * t300;
t200 = Icges(6,4) * t291 + t292 * t446;
t202 = Icges(5,5) * t291 + t269 * t292;
t464 = -t471 * t292 + t200 + t202;
t199 = -Icges(6,4) * t292 + t291 * t446;
t371 = t291 * t300;
t256 = Icges(5,4) * t371;
t370 = t291 * t302;
t201 = Icges(5,1) * t370 - Icges(5,5) * t292 - t256;
t463 = -Icges(5,2) * t370 - t320 * t291 + t199 + t201 - t256;
t368 = t292 * t302;
t255 = Icges(6,5) * t368;
t192 = Icges(6,6) * t291 + Icges(6,3) * t369 + t255;
t198 = Icges(5,6) * t291 + t265 * t292;
t462 = -Icges(6,1) * t369 - t268 * t292 + t192 - t198 + t255;
t261 = Icges(6,3) * t300 + t386;
t191 = -Icges(6,6) * t292 + t261 * t291;
t197 = Icges(5,4) * t370 - Icges(5,2) * t371 - Icges(5,6) * t292;
t461 = t473 * t291 - t191 + t197;
t296 = cos(t299);
t398 = pkin(2) * t296;
t146 = t398 + (rSges(6,2) + pkin(7)) * t291 + t442 * t292;
t167 = -t292 * t331 + t390 * t368;
t64 = t145 * t474 + t146 * t167;
t391 = rSges(5,1) * t302;
t330 = pkin(3) + t391;
t338 = rSges(5,2) * t371 + t292 * rSges(5,3);
t160 = -t330 * t291 + t329 + t338;
t258 = rSges(5,2) * t369;
t161 = t398 - t258 + t330 * t292 + (rSges(5,3) + pkin(7)) * t291;
t272 = rSges(5,1) * t300 + rSges(5,2) * t302;
t226 = t272 * t291;
t228 = t272 * t292;
t84 = t160 * t226 - t161 * t228;
t460 = -m(5) * t84 - m(6) * t64;
t436 = m(5) / 0.2e1;
t435 = m(6) / 0.2e1;
t459 = -t300 / 0.2e1;
t400 = cos(qJ(1)) * pkin(1);
t405 = m(3) * (t400 * (-rSges(3,1) * t295 - rSges(3,2) * t296) + (t296 * rSges(3,1) - t295 * rSges(3,2)) * t401);
t403 = m(4) * (t400 * (-rSges(4,1) * t291 - rSges(4,2) * t292 - t399) + (t292 * rSges(4,1) - t291 * rSges(4,2) + t398) * t401);
t332 = t146 * t369;
t80 = -t145 * t371 + t332;
t456 = t80 * m(6) * qJD(2);
t263 = Icges(6,4) * t302 + Icges(6,6) * t300;
t375 = t263 * t291;
t195 = -Icges(6,2) * t292 + t375;
t182 = t291 * t195;
t102 = t191 * t369 + t199 * t368 + t182;
t455 = t102 * t292;
t454 = (t469 - t471) * t302 + (t261 - t472) * t300;
t451 = (t191 * t300 + t199 * t302) * t291;
t155 = t160 - t401;
t156 = t161 + t400;
t74 = -t161 * t155 + t156 * t160;
t449 = t465 * t291;
t448 = t465 * t292;
t445 = -t464 * t300 + t462 * t302;
t444 = t463 * t300 + t461 * t302;
t136 = t146 + t400;
t186 = t470 * t292;
t397 = ((-t136 + t146) * t186 + t467 * t474) * t435 + ((-t156 + t161) * t292 + (t155 - t160) * t291) * t272 * t436;
t62 = t135 * t474 + t136 * t167;
t81 = t155 * t226 - t156 * t228;
t443 = (t64 + t62) * t435 + (t84 + t81) * t436;
t313 = t471 * t459 + t469 * t300 / 0.2e1 + (-t261 / 0.2e1 + t472 / 0.2e1) * t302;
t441 = 4 * qJD(1);
t439 = 2 * qJD(4);
t432 = m(5) * t74;
t430 = m(5) * t81;
t125 = t136 * t369;
t424 = m(6) * (-t467 * t371 + t125 - t332);
t423 = m(6) * (t125 + t332 + (-t135 - t145) * t371);
t42 = -t146 * t135 + t136 * t145;
t422 = m(6) * t42;
t351 = t167 * t371 + t475;
t354 = t135 * t368 + t136 * t370;
t420 = m(6) * (t351 + t354);
t353 = t145 * t368 + t146 * t370;
t419 = m(6) * (t351 + t353);
t326 = t186 * t371 - t475;
t418 = m(6) * (t326 + t354);
t417 = m(6) * (t326 + t353);
t416 = m(6) * t62;
t229 = t335 * t300;
t412 = t229 / 0.2e1;
t411 = -t291 / 0.2e1;
t410 = t291 / 0.2e1;
t409 = -t292 / 0.2e1;
t374 = t263 * t292;
t196 = Icges(6,2) * t291 + t374;
t103 = t192 * t369 + t291 * t196 + t200 * t368;
t380 = t195 * t292;
t315 = t103 + t380;
t325 = -t192 * t371 + t196 * t292 - t200 * t370;
t19 = (t103 - t315) * t292 + (t102 + t325 - t182) * t291;
t193 = Icges(5,5) * t370 - Icges(5,6) * t371 - Icges(5,3) * t292;
t349 = -t291 * t193 - t201 * t368;
t104 = -t197 * t369 - t349;
t262 = Icges(5,5) * t302 - Icges(5,6) * t300;
t376 = t262 * t292;
t194 = Icges(5,3) * t291 + t376;
t348 = t291 * t194 + t202 * t368;
t105 = -t198 * t369 + t348;
t327 = t198 * t300 - t193;
t173 = t202 * t370;
t328 = t194 * t292 - t173;
t20 = (t327 * t292 + t105 - t348) * t292 + (t327 * t291 + t104 + t328) * t291;
t98 = -t380 + t451;
t21 = -t455 + (t315 + t98 - t451) * t291;
t101 = -t198 * t371 - t328;
t379 = t197 * t300;
t22 = (t101 - t173 + (t194 + t379) * t292 + t349) * t292 + t348 * t291;
t58 = -t291 * t325 - t292 * t98;
t59 = -(-(-t201 * t302 + t379) * t291 - t193 * t292) * t292 + t101 * t291;
t60 = t103 * t291 - t455;
t61 = -t104 * t292 + t105 * t291;
t2 = (t61 / 0.2e1 - t22 / 0.2e1 + t60 / 0.2e1 - t21 / 0.2e1) * t292 + (t20 / 0.2e1 + t59 / 0.2e1 + t19 / 0.2e1 + t58 / 0.2e1) * t291;
t406 = t2 * qJD(4);
t393 = m(6) * qJD(4);
t392 = m(6) * qJD(5);
t363 = t300 * t302;
t352 = (-t167 - t186) * t474;
t350 = -t186 * t368 - t370 * t474;
t339 = t335 * t363;
t206 = (t412 + t459) * m(6);
t334 = t206 * qJD(3);
t77 = -t135 * t371 + t125;
t333 = m(6) * t77 * qJD(1);
t314 = (-t226 * t292 + t228 * t291) * t272;
t306 = t313 + t443;
t305 = -t313 + ((t191 + t197) * t302 + (-t199 + t201) * t300) * (t410 + t411);
t304 = ((t21 + t22) * t292 / 0.2e1 + (t19 + t20 + t58 + t59) * t411 + (t262 * t291 + t454 * t292 + t462 * t300 + t464 * t302 + t375) * t410 + (t454 * t291 - t461 * t300 + t463 * t302 - t374 - t376 + t60 + t61) * t409) * qJD(4);
t276 = -rSges(5,2) * t300 + t391;
t205 = m(6) * t412 + t300 * t435;
t190 = t339 - t363;
t187 = t447 * t292;
t185 = t447 * t291;
t152 = -t226 * t291 - t228 * t292;
t108 = t167 * t292 - t288 * t470;
t90 = t335 * t447;
t63 = t229 * t90 + t350;
t56 = t417 / 0.2e1;
t50 = t418 / 0.2e1;
t48 = t419 / 0.2e1;
t45 = t420 / 0.2e1;
t37 = t423 / 0.2e1;
t36 = t424 / 0.2e1;
t27 = t313 - t460;
t26 = t313 + t416 + t430;
t23 = t403 + t405 + t422 + t432;
t14 = t56 - t419 / 0.2e1;
t13 = t56 + t48;
t12 = t48 - t417 / 0.2e1;
t11 = t50 - t420 / 0.2e1;
t10 = t50 + t45;
t9 = t45 - t418 / 0.2e1;
t8 = t37 - t424 / 0.2e1;
t7 = t37 + t36;
t6 = t36 - t423 / 0.2e1;
t5 = t306 + t397;
t4 = t306 - t397;
t3 = t305 + t397 - t443;
t1 = [t23 * qJD(2) + t26 * qJD(4) + t77 * t392, t23 * qJD(1) + t5 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t403 / 0.2e1 + t405 / 0.2e1 + t42 * t435 + t74 * t436) * qJD(2), 0, t26 * qJD(1) + t5 * qJD(2) + t10 * qJD(5) + ((-t135 * t187 - t136 * t185 + t352) * t435 + ((-t155 * t292 - t156 * t291) * t276 + t314) * t436) * t439 + t304, t7 * qJD(2) + t10 * qJD(4) + t333; t4 * qJD(4) + t8 * qJD(5) + (-t422 / 0.4e1 - t432 / 0.4e1 - t403 / 0.4e1 - t405 / 0.4e1) * t441, t27 * qJD(4) + t80 * t392, 0, t4 * qJD(1) + t27 * qJD(2) + t13 * qJD(5) + ((-t145 * t187 - t146 * t185 + t352) * t435 + ((-t160 * t292 - t161 * t291) * t276 + t314) * t436) * t439 + t304, t8 * qJD(1) + t13 * qJD(4) + t456; 0, 0, 0, (t108 * t435 + t152 * t436) * t439 + t205 * qJD(5), t205 * qJD(4); t305 * qJD(1) + t3 * qJD(2) + t11 * qJD(5) + (-t416 / 0.4e1 - t430 / 0.4e1) * t441 + t406, t3 * qJD(1) + t14 * qJD(5) + t406 + (t305 + t460) * qJD(2), qJD(5) * t206, (m(5) * (t272 * t276 * t335 + (t291 * (rSges(5,1) * t370 - t338) + t292 * (rSges(5,1) * t368 + t291 * rSges(5,3) - t258)) * t152) + m(6) * (t108 * t90 + t185 * t474 + t186 * t187) + (t448 * t288 + (t444 * t292 + (t445 - t449) * t291) * t292) * t410 + (t449 * t289 + (t445 * t291 + (t444 - t448) * t292) * t291) * t409) * qJD(4) + t63 * t392 + (qJD(1) + qJD(2)) * t2, t11 * qJD(1) + t14 * qJD(2) + t334 + t63 * t393 + (-t229 * t302 - t190 + t339) * t392; t6 * qJD(2) + t9 * qJD(4) - t333, t6 * qJD(1) + t12 * qJD(4) - t456, -t206 * qJD(4), t9 * qJD(1) + t12 * qJD(2) - t334 + (-t108 * t302 + (-t185 * t291 - t187 * t292 + t90) * t300 - t63 + t350) * t393 + t190 * t392, t190 * t393;];
Cq = t1;
