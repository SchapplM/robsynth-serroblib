% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:37
% DurationCPUTime: 13.28s
% Computational Cost: add. (6019->445), mult. (16434->575), div. (0->0), fcn. (16278->6), ass. (0->214)
t480 = Icges(5,4) - Icges(6,5);
t469 = Icges(5,1) + Icges(6,1);
t468 = Icges(6,4) + Icges(5,5);
t467 = Icges(5,2) + Icges(6,3);
t466 = Icges(5,6) - Icges(6,6);
t264 = sin(qJ(1));
t263 = cos(pkin(7));
t383 = sin(pkin(7));
t397 = sin(qJ(4));
t398 = cos(qJ(4));
t220 = -t263 * t397 + t383 * t398;
t265 = cos(qJ(1));
t206 = t220 * t265;
t219 = t263 * t398 + t383 * t397;
t207 = t219 * t265;
t204 = t220 * t264;
t205 = t219 * t264;
t492 = t480 * t205;
t448 = t467 * t204 + t466 * t265 + t492;
t493 = t480 * t204;
t464 = -t469 * t205 - t468 * t265 - t493;
t478 = t206 * t448 - t207 * t464;
t494 = Icges(6,2) + Icges(5,3);
t482 = t466 * t204 + t468 * t205 + t494 * t265;
t473 = -t482 * t264 + t478;
t491 = t480 * t206;
t490 = t480 * t207;
t447 = -t467 * t206 + t466 * t264 - t490;
t465 = t469 * t207 - t468 * t264 + t491;
t487 = t480 * t220;
t486 = t480 * t219;
t485 = t204 * t448 - t205 * t464;
t446 = t466 * t206 + t468 * t207 - t494 * t264;
t443 = t467 * t219 - t487;
t481 = -t466 * t219 + t468 * t220;
t442 = t469 * t220 - t486;
t458 = rSges(6,1) + pkin(4);
t449 = rSges(6,3) + qJ(5);
t479 = -t447 * t206 + t465 * t207;
t477 = t447 * t204 - t465 * t205;
t476 = t448 * t219 + t220 * t464;
t475 = t482 * t265 + t485;
t474 = -t446 * t265 + t477;
t472 = -t443 * t204 + t442 * t205 + t481 * t265;
t471 = -t443 * t206 + t442 * t207 - t481 * t264;
t439 = t446 * t264 - t479;
t347 = qJD(4) * t265;
t411 = t220 * qJD(1);
t130 = -t219 * t347 - t264 * t411;
t351 = qJD(1) * t264;
t437 = qJD(4) * t220;
t131 = -t219 * t351 + t265 * t437;
t350 = qJD(1) * t265;
t455 = -t467 * t130 - t480 * t131 + t466 * t350;
t132 = qJD(4) * t205 - t265 * t411;
t133 = qJD(1) * t207 + t264 * t437;
t454 = t467 * t132 - t480 * t133 + t466 * t351;
t451 = -t480 * t130 - t469 * t131 + t468 * t350;
t450 = -t480 * t132 + t469 * t133 - t468 * t351;
t211 = t219 * qJD(4);
t463 = t480 * t211 + t467 * t437;
t462 = -t468 * t211 - t466 * t437;
t461 = -t469 * t211 - t480 * t437;
t362 = t449 * t219 + t458 * t220;
t459 = -(t130 * t466 + t131 * t468 - t350 * t494) * t264 + (-t132 * t466 + t133 * t468 - t351 * t494) * t265;
t456 = t471 * qJD(1);
t441 = (t264 * t439 + t265 * t473) * qJD(4);
t440 = (t474 * t264 + t265 * t475) * qJD(4);
t438 = t472 * qJD(1);
t348 = qJD(4) * t264;
t436 = -t475 + t439;
t435 = 0.2e1 * qJD(4);
t434 = t438 + t440;
t433 = t441 + t456;
t432 = t211 * t464 + t219 * t454 + t220 * t450 - t437 * t448;
t431 = t211 * t465 - t219 * t455 + t220 * t451 - t437 * t447;
t430 = -t130 * t443 + t131 * t442 - t206 * t463 + t207 * t461 - t264 * t462 - t350 * t481;
t429 = t132 * t443 + t133 * t442 - t204 * t463 + t205 * t461 + t265 * t462 - t351 * t481;
t427 = -t219 * t447 - t220 * t465;
t188 = qJD(5) * t204;
t424 = -t449 * t132 - t133 * t458 + t188;
t189 = qJD(5) * t206;
t423 = -t449 * t130 + t131 * t458 - t189;
t419 = -t265 * rSges(6,2) + t449 * t204 - t205 * t458;
t418 = -t449 * t206 + t207 * t458;
t417 = -t219 * t468 - t220 * t466;
t320 = qJD(3) * t383;
t238 = t264 * t320;
t255 = qJD(2) * t265;
t356 = t238 - t255;
t416 = -t188 + t356;
t231 = t265 * pkin(1) + t264 * qJ(2);
t209 = qJD(1) * t231 - t255;
t323 = t383 * qJ(3);
t396 = pkin(2) * t263;
t286 = t323 + t396;
t415 = -t286 * t350 - t209 - t238;
t325 = t265 * t383;
t375 = t263 * t265;
t280 = rSges(3,1) * t375 - rSges(3,2) * t325 + t264 * rSges(3,3);
t414 = t231 + t280;
t245 = pkin(3) * t375;
t315 = pkin(2) * t375 + t265 * t323 + t231;
t296 = t245 + t315;
t413 = (t205 * t467 + t464 - t493) * t265 + (-t207 * t467 + t465 + t491) * t264;
t412 = (t204 * t469 - t448 - t492) * t265 + (-t206 * t469 - t447 + t490) * t264;
t410 = (t204 * t468 - t205 * t466) * t265 + (-t206 * t468 + t207 * t466) * t264;
t259 = t264 * rSges(4,2);
t409 = rSges(4,1) * t375 + rSges(4,3) * t325 + t259 + t315;
t408 = -t219 * t469 + t443 - t487;
t407 = t220 * t467 - t442 + t486;
t318 = qJD(4) * t362;
t406 = t347 * t362 - t189;
t405 = qJD(1) ^ 2;
t403 = t264 / 0.2e1;
t402 = -t265 / 0.2e1;
t400 = -rSges(6,2) - pkin(6);
t399 = -rSges(5,3) - pkin(6);
t395 = t265 * pkin(6);
t394 = -qJD(1) / 0.2e1;
t392 = rSges(6,2) * t350 - t423;
t391 = rSges(6,2) * t351 + t424;
t386 = rSges(4,1) * t263;
t345 = qJD(1) * qJD(2);
t249 = t265 * t345;
t333 = t383 * rSges(4,3);
t287 = t333 + t386;
t304 = qJD(1) * t320;
t76 = t249 - t264 * t304 - t405 * (t265 * t287 + t259) + t415 * qJD(1);
t385 = t76 * t264;
t346 = qJD(5) * t219;
t384 = -t458 * t211 + t437 * t449 + t346;
t376 = t263 * t264;
t372 = -rSges(6,2) * t264 + t418;
t371 = t131 * rSges(5,1) + t130 * rSges(5,2);
t369 = -t458 * t204 - t449 * t205;
t368 = t458 * t206 + t449 * t207;
t363 = -t458 * t219 + t449 * t220;
t360 = t207 * rSges(5,1) + t206 * rSges(5,2);
t254 = qJD(2) * t264;
t353 = qJ(2) * t350 + t254;
t359 = qJD(1) * (-pkin(1) * t351 + t353) + t264 * t345;
t213 = t286 * t264;
t257 = t265 * qJ(2);
t229 = t264 * pkin(1) - t257;
t358 = -t213 - t229;
t321 = qJD(1) * t383;
t306 = t264 * t321;
t357 = rSges(3,2) * t306 + rSges(3,3) * t350;
t239 = t265 * t320;
t355 = t239 + t254;
t334 = t383 * rSges(3,2);
t354 = t265 * rSges(3,3) + t264 * t334;
t352 = -qJD(1) * t229 + t254;
t349 = qJD(3) * t263;
t344 = rSges(3,1) * t376;
t253 = pkin(6) * t351;
t343 = -qJD(1) * t245 + t253 + t415;
t222 = pkin(3) * t376 + t395;
t342 = -t222 + t358;
t341 = -pkin(6) * t264 + t296;
t340 = t239 + t353;
t338 = t253 - t356;
t185 = rSges(5,1) * t220 - rSges(5,2) * t219;
t167 = t185 * t347;
t335 = -rSges(3,1) * t263 - pkin(1);
t329 = t348 / 0.2e1;
t328 = -t347 / 0.2e1;
t326 = t257 - t395;
t30 = (t342 + t419) * qJD(1) + t355 + t406;
t324 = t30 * t362;
t319 = qJD(4) * t384;
t317 = qJD(1) * (-t286 * t351 + t239) + t265 * t304 + t359;
t316 = -qJD(1) * t213 + t239 + t352;
t302 = t133 * rSges(5,1) - t132 * rSges(5,2);
t118 = t205 * rSges(5,1) + t204 * rSges(5,2) + t265 * rSges(5,3);
t44 = t167 + (-t118 + t342) * qJD(1) + t355;
t122 = -t264 * rSges(5,3) + t360;
t45 = t185 * t348 + (t122 + t341) * qJD(1) + t356;
t300 = t45 * t264 + t44 * t265;
t299 = t341 + t372;
t298 = -t405 * t222 + t317;
t297 = -qJD(1) * t222 + t316;
t291 = -t118 * t264 - t122 * t265;
t275 = -pkin(1) - t333 - t323;
t274 = -pkin(3) * t263 - pkin(1) - t286;
t72 = -rSges(5,3) * t350 + t371;
t74 = -rSges(5,3) * t351 + t302;
t272 = -t264 * t74 - t265 * t72 + (-t118 * t265 + t122 * t264) * qJD(1);
t271 = (-rSges(4,1) - pkin(2)) * t263 + t275;
t261 = t265 * rSges(4,2);
t252 = rSges(4,2) * t350;
t202 = t344 - t354;
t201 = t264 * t287 - t261;
t164 = -rSges(5,1) * t211 - rSges(5,2) * t437;
t156 = qJD(1) * t414 - t255;
t155 = t254 + (-t202 - t229) * qJD(1);
t153 = rSges(5,1) * t206 - rSges(5,2) * t207;
t148 = rSges(5,1) * t204 - rSges(5,2) * t205;
t95 = t249 + (-qJD(1) * t280 - t209) * qJD(1);
t94 = qJD(1) * (-qJD(1) * t344 + t357) + t359;
t93 = qJD(1) * t409 + t356;
t92 = (-t201 + t358) * qJD(1) + t355;
t75 = qJD(1) * (-t287 * t351 + t252) + t317;
t58 = qJD(4) * t291 - t349;
t32 = -t349 + t346 + (t264 * t419 - t265 * t372) * qJD(4);
t31 = t299 * qJD(1) + t264 * t318 + t416;
t26 = t164 * t347 + t249 + (-t74 + (-qJD(4) * t185 - t320) * t264 + t343) * qJD(1);
t25 = t164 * t348 + (t72 + t167) * qJD(1) + t298;
t24 = t272 * qJD(4);
t15 = qJD(5) * t132 + t264 * t319 + (t265 * t318 - t392) * qJD(1) + t298;
t14 = -qJD(5) * t130 + t249 + t265 * t319 + ((-t320 - t318) * t264 + t343 + t391) * qJD(1);
t1 = qJD(5) * t437 + (t392 * t265 + t391 * t264 + (t264 * t372 + t265 * t419) * qJD(1)) * qJD(4);
t2 = [((t478 * t265 + (t436 + t485) * t264) * qJD(4) + t456) * t328 + (-t442 * t211 + t463 * t219 + t461 * t220 + t443 * t437) * qJD(1) + (t324 * t348 + t14 * (t326 + t419) + t15 * (t296 + t418) + (t14 * t274 + t15 * t400) * t264 + (-t297 + t340 - t406 + t423) * t31 + (t338 + t416 + t424) * t30 + (t299 * t30 - t31 * t419 + (t274 * t30 + t31 * t400) * t265 + (t30 * (rSges(6,2) - qJ(2)) + t31 * t274) * t264) * qJD(1)) * m(6) + (t26 * (-t118 + t326) + t44 * (-t302 + t338) + t25 * (t296 + t360) + t45 * (t340 + t371) + (t25 * t399 + t26 * t274) * t264 + ((t274 * t44 + t399 * t45) * t265 + (t44 * (rSges(5,3) - qJ(2)) + t45 * t274) * t264) * qJD(1) - (-qJD(1) * t118 + t167 + t297 - t44) * t45) * m(5) + (-(-qJD(1) * t201 + t316 - t92) * t93 + t76 * (t257 + t261) - t92 * t356 + t75 * t409 + t93 * (t252 + t340) + t271 * t385 + (t92 * t271 * t265 + (t92 * (-rSges(4,2) - qJ(2)) + t93 * (t275 - t386 - t396)) * t264) * qJD(1)) * m(4) + (t95 * (t264 * t335 + t257 + t354) + t155 * t255 + t94 * t414 + t156 * (t353 + t357) + (t155 * (t334 + t335) * t265 + (t155 * (-rSges(3,3) - qJ(2)) + t156 * t335) * t264) * qJD(1) - (-qJD(1) * t202 - t155 + t352) * t156) * m(3) - (t430 - t431) * t348 / 0.2e1 + (((t436 + t479) * t265 - t477 * t264) * qJD(4) + t434 - t438) * t329 + (t429 + t432 + t433) * t347 / 0.2e1 + ((-t427 + t471) * t265 + (t472 - t476) * t264) * qJD(4) * t394; 0.2e1 * (t14 * t403 + t15 * t402) * m(6) + 0.2e1 * (t25 * t402 + t26 * t403) * m(5) + 0.2e1 * (t385 / 0.2e1 + t75 * t402) * m(4) + 0.2e1 * (t402 * t94 + t403 * t95) * m(3); (t14 * t325 - t1 * t263 + (t15 * t383 - t30 * t321) * t264 + t30 * t306) * m(6) + (t26 * t325 - t24 * t263 + (t25 * t383 - t321 * t44) * t264 + t306 * t44) * m(5) + (t76 * t325 + (-t321 * t92 + t383 * t75) * t264 + t306 * t92) * m(4); ((t219 * t413 + t220 * t412) * qJD(4) + (t219 * t407 + t220 * t408) * qJD(1)) * t394 + (t432 * t265 + t431 * t264 + (t476 * t264 + t427 * t265) * qJD(1)) * qJD(1) / 0.2e1 + ((-t206 * t413 + t207 * t412 - t264 * t410) * qJD(4) + (-t206 * t407 + t207 * t408 - t264 * t417) * qJD(1)) * t329 + ((-t204 * t413 + t205 * t412 + t265 * t410) * qJD(4) + (-t204 * t407 + t205 * t408 + t265 * t417) * qJD(1)) * t328 + ((t14 * t362 + t30 * t384 - t1 * t372 + t32 * t392 + (t31 * t362 + t32 * t419) * qJD(1)) * t265 + (t15 * t362 + t31 * t384 + t1 * t419 + t32 * t391 + (t32 * t372 - t324) * qJD(1)) * t264 - (t205 * t31 + t207 * t30 + t220 * t32) * qJD(5) - (t30 * t369 + t31 * t368) * qJD(1) - ((t30 * t363 - t32 * t368) * t265 + (t31 * t363 + t32 * t369) * t264) * qJD(4)) * m(6) + (-(-t148 * t44 + t153 * t45) * qJD(1) - (t58 * (-t148 * t264 - t153 * t265) + t300 * (-rSges(5,1) * t219 - rSges(5,2) * t220)) * qJD(4) + t24 * t291 + t58 * t272 + t300 * t164 + (t25 * t264 + t26 * t265 + (-t44 * t264 + t45 * t265) * qJD(1)) * t185) * m(5) - (t430 * qJD(1) + ((t439 * qJD(1) + t448 * t130 - t131 * t464 - t454 * t206 + t450 * t207 - t350 * t482) * t265 + (-qJD(1) * t473 + t447 * t130 - t131 * t465 + t455 * t206 + t451 * t207 + t446 * t350 - t459) * t264) * t435) * t264 / 0.2e1 + (t429 * qJD(1) + ((-qJD(1) * t475 - t447 * t132 - t133 * t465 + t455 * t204 + t451 * t205 + t446 * t351) * t264 + (qJD(1) * t474 - t448 * t132 - t133 * t464 - t454 * t204 + t450 * t205 - t351 * t482 + t459) * t265) * t435) * t265 / 0.2e1 - (t434 + t440) * t351 / 0.2e1 - (t433 + t441) * t350 / 0.2e1; (t1 * t219 - t130 * t30 + t132 * t31 - t14 * t206 - t15 * t204 + t437 * t32 - (t204 * t30 - t206 * t31) * qJD(1) - (t32 * (t204 * t264 + t206 * t265) + (t31 * t264 + t30 * t265) * t219) * qJD(4)) * m(6);];
tauc = t2(:);
