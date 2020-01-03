% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:03
% EndTime: 2019-12-31 17:52:19
% DurationCPUTime: 12.21s
% Computational Cost: add. (6189->427), mult. (14277->519), div. (0->0), fcn. (14850->6), ass. (0->242)
t381 = sin(pkin(7));
t382 = cos(pkin(7));
t404 = sin(qJ(1));
t405 = cos(qJ(1));
t172 = -t381 * t404 - t382 * t405;
t173 = t405 * t381 - t404 * t382;
t469 = Icges(5,3) + Icges(6,3);
t225 = sin(qJ(4));
t226 = cos(qJ(4));
t282 = Icges(6,5) * t226 - Icges(6,6) * t225;
t284 = Icges(5,5) * t226 - Icges(5,6) * t225;
t492 = t282 + t284;
t453 = t469 * t172 + t492 * t173;
t377 = Icges(6,4) * t226;
t286 = -Icges(6,2) * t225 + t377;
t100 = Icges(6,6) * t172 + t173 * t286;
t379 = Icges(5,4) * t226;
t288 = -Icges(5,2) * t225 + t379;
t103 = Icges(5,6) * t172 + t173 * t288;
t449 = t100 + t103;
t378 = Icges(6,4) * t225;
t290 = Icges(6,1) * t226 - t378;
t106 = Icges(6,5) * t172 + t173 * t290;
t380 = Icges(5,4) * t225;
t292 = Icges(5,1) * t226 - t380;
t109 = Icges(5,5) * t172 + t173 * t292;
t477 = t106 + t109;
t496 = t449 * t225 - t226 * t477;
t459 = t453 * t172 - t173 * t496;
t452 = -t492 * t172 + t173 * t469;
t493 = t452 * t173;
t285 = Icges(6,2) * t226 + t378;
t287 = Icges(5,2) * t226 + t380;
t491 = -t285 - t287;
t289 = Icges(6,1) * t225 + t377;
t291 = Icges(5,1) * t225 + t379;
t490 = t289 + t291;
t108 = Icges(6,5) * t173 - t172 * t290;
t111 = Icges(5,5) * t173 - t172 * t292;
t476 = t108 + t111;
t488 = t476 * t226;
t423 = t491 * t225 + t490 * t226;
t283 = Icges(5,5) * t225 + Icges(5,6) * t226;
t450 = t283 * t172;
t281 = Icges(6,5) * t225 + Icges(6,6) * t226;
t451 = t281 * t172;
t466 = t173 * t423 + t450 + t451;
t487 = t466 * qJD(1);
t433 = t477 * t225 + t226 * t449;
t102 = Icges(6,6) * t173 - t172 * t286;
t105 = Icges(5,6) * t173 - t172 * t288;
t486 = t102 + t105;
t470 = Icges(5,6) + Icges(6,6);
t471 = Icges(5,5) + Icges(6,5);
t485 = t471 * t225 + t470 * t226;
t479 = Icges(5,4) + Icges(6,4);
t484 = t479 * t225 + (Icges(5,2) + Icges(6,2)) * t226;
t483 = (Icges(5,1) + Icges(6,1)) * t225 + t479 * t226;
t475 = (t286 + t288) * qJD(4);
t474 = (t290 + t292) * qJD(4);
t362 = t172 * t226;
t417 = -t476 * t362 + t493;
t473 = t452 * t172 + t173 * t488;
t361 = t173 * t225;
t458 = t486 * t361 - t473;
t363 = t172 * t225;
t434 = -t453 * t173 + t477 * t362 - t449 * t363;
t457 = t486 * t363 + t417;
t121 = t281 * t173;
t123 = t283 * t173;
t472 = t423 * t172 - t121 - t123;
t468 = t486 * t225;
t224 = -qJ(5) - pkin(6);
t465 = rSges(6,3) - t224;
t162 = t172 * qJD(1);
t333 = qJD(4) * t225;
t321 = t172 * t333;
t163 = t173 * qJD(1);
t367 = t163 * t226;
t261 = t321 + t367;
t332 = qJD(4) * t226;
t368 = t163 * t225;
t262 = t172 * t332 - t368;
t464 = (t469 * t162 + t471 * t261 + t470 * t262) * t173;
t463 = t492 * qJD(4);
t462 = t474 * t226 - t475 * t225 + (-t490 * t225 + t491 * t226) * qJD(4);
t444 = t468 - t488;
t461 = -t281 - t283;
t460 = rSges(6,1) + pkin(4);
t318 = t173 * t333;
t263 = -t162 * t226 + t318;
t264 = t162 * t225 + t173 * t332;
t456 = t469 * t163 + t471 * t263 + t470 * t264;
t454 = t472 * qJD(1);
t66 = t261 * rSges(5,1) + rSges(5,2) * t262 + t162 * rSges(5,3);
t446 = (-t434 * t172 + t457 * t173) * qJD(4);
t445 = (-t459 * t172 + t458 * t173) * qJD(4);
t218 = t405 * qJ(2);
t326 = t404 * pkin(1);
t199 = t326 - t218;
t215 = qJD(2) * t404;
t317 = qJD(1) * t405;
t339 = qJ(2) * t317 + t215;
t419 = qJD(1) * t199 - t215 + t339;
t222 = t405 * pkin(2);
t337 = t405 * pkin(1) + t404 * qJ(2);
t425 = t222 + t337;
t442 = -t172 * rSges(4,1) - t173 * rSges(4,2) + t425;
t441 = 0.2e1 * qJD(4);
t440 = t445 + t487;
t439 = t446 + t454;
t438 = qJD(4) * t496 + t485 * t163 + t483 * t263 + t484 * t264;
t437 = t444 * qJD(4) - t485 * t162 - t483 * t261 - t484 * t262;
t436 = t423 * t162 + t461 * t163 + t463 * t172 + t462 * t173;
t435 = t461 * t162 - t423 * t163 + t462 * t172 - t463 * t173;
t432 = t476 * t225 + t486 * t226;
t403 = pkin(4) * t226;
t213 = pkin(3) + t403;
t427 = -rSges(6,1) * t362 + rSges(6,2) * t363 - t172 * t213 + t465 * t173;
t325 = t404 * pkin(2);
t244 = t173 * rSges(4,1) - t172 * rSges(4,2) - t325;
t426 = -t199 + t244;
t254 = -t326 - t325;
t424 = t254 + t325;
t164 = t172 * qJD(5);
t422 = -t465 * t163 + t164;
t252 = t405 * rSges(3,1) + t404 * rSges(3,3);
t421 = t337 + t252;
t170 = t172 * pkin(6);
t416 = t173 * pkin(3) + t170;
t420 = -qJD(1) * t416 + t419;
t418 = rSges(6,1) * t367 + t465 * t162 + t163 * t213;
t335 = qJD(4) * t173;
t165 = qJD(5) * t173;
t308 = t163 * pkin(3) + t162 * pkin(6);
t298 = rSges(5,1) * t226 - rSges(5,2) * t225;
t388 = t172 * rSges(5,3);
t113 = t173 * t298 + t388;
t392 = rSges(6,2) * t226;
t415 = -t225 * t460 - t392;
t216 = qJD(2) * t405;
t414 = -qJD(1) * t425 + t216;
t353 = t287 * t173 - t109;
t357 = -t291 * t173 - t103;
t413 = t225 * t353 + t226 * t357;
t355 = t285 * t173 - t106;
t359 = -t289 * t173 - t100;
t412 = t225 * t355 + t226 * t359;
t296 = rSges(6,1) * t226 - rSges(6,2) * t225;
t387 = t172 * rSges(6,3);
t399 = pkin(3) - t213;
t411 = t387 + (t296 - t399) * t173;
t228 = qJD(1) ^ 2;
t410 = t162 / 0.2e1;
t409 = t163 / 0.2e1;
t159 = t163 * pkin(6);
t398 = rSges(6,1) * t263 + rSges(6,2) * t264 + pkin(4) * t318 + t162 * t399 - t159 - t422;
t397 = rSges(6,2) * t262 + t460 * t321 + t165 - t308 + t418;
t268 = -t199 - t325 + t416;
t295 = rSges(6,1) * t225 + t392;
t309 = pkin(4) * t225 + t295;
t293 = qJD(4) * t309;
t364 = t172 * t224;
t384 = t170 + t364 - t411;
t24 = t165 + t215 + t172 * t293 + (t268 - t384) * qJD(1);
t391 = t162 * t24;
t390 = t163 * rSges(5,3);
t385 = (-pkin(6) - t224) * t172 + t411;
t171 = t172 * pkin(3);
t142 = pkin(6) * t173 - t171;
t383 = -t142 + t427;
t358 = -t289 * t172 + t102;
t356 = -t291 * t172 + t105;
t354 = t285 * t172 + t108;
t352 = t287 * t172 + t111;
t160 = qJD(1) * t337 - t216;
t351 = t162 * pkin(3) - t159 - t160;
t316 = qJD(1) * t404;
t347 = (-pkin(1) * t316 + t215 + t339) * qJD(1);
t346 = t163 * rSges(4,1) - t162 * rSges(4,2);
t345 = t164 + t216;
t343 = -t285 + t290;
t342 = -t286 - t289;
t341 = -t287 + t292;
t340 = -t288 - t291;
t336 = qJD(4) * t172;
t297 = rSges(5,1) * t225 + rSges(5,2) * t226;
t334 = qJD(4) * t297;
t331 = t282 * qJD(1);
t330 = t284 * qJD(1);
t329 = qJD(4) ^ 2 * t403;
t328 = -t405 / 0.2e1;
t327 = t404 / 0.2e1;
t324 = t404 * rSges(3,1);
t312 = t336 / 0.2e1;
t311 = -t335 / 0.2e1;
t306 = pkin(4) * t361 + t295 * t173;
t305 = pkin(4) * t363 + t295 * t172;
t300 = t24 * t309;
t299 = t162 * rSges(4,1) + t163 * rSges(4,2);
t43 = t172 * t334 + t215 + (t113 + t268) * qJD(1);
t117 = -rSges(5,1) * t362 + rSges(5,2) * t363 + t173 * rSges(5,3);
t267 = t142 + t425;
t44 = t173 * t334 - t216 + (t117 + t267) * qJD(1);
t294 = -t172 * t43 - t173 * t44;
t272 = -t113 * t173 + t117 * t172;
t211 = qJD(1) * t216;
t269 = -t222 * t228 + t211;
t266 = pkin(3) + t298;
t265 = t213 + t296;
t260 = -t228 * t325 + t347;
t253 = -t324 - t326;
t249 = t267 + t383;
t248 = qJD(1) * t308 + t260;
t247 = t218 + t254;
t246 = t225 * t354 + t226 * t358;
t245 = t225 * t352 + t226 * t356;
t242 = (t225 * t342 + t226 * t343) * qJD(1);
t241 = (t225 * t340 + t226 * t341) * qJD(1);
t64 = rSges(5,1) * t263 + rSges(5,2) * t264 + t390;
t240 = -t113 * t162 - t117 * t163 + t172 * t66 + t173 * t64;
t220 = t405 * rSges(3,3);
t214 = rSges(3,3) * t317;
t182 = t298 * qJD(4);
t181 = t296 * qJD(4);
t136 = t297 * t172;
t134 = t297 * t173;
t119 = -qJD(1) * t160 - t228 * t252 + t211;
t118 = qJD(1) * (-rSges(3,1) * t316 + t214) + t347;
t92 = qJD(1) * t426 + t215;
t74 = (-t160 + t299) * qJD(1) + t269;
t73 = qJD(1) * t346 + t260;
t42 = qJD(4) * t272 - qJD(3);
t31 = (-t163 * t297 + t172 * t182) * qJD(4) + (-t64 + t351) * qJD(1) + t269;
t30 = qJD(1) * t66 + (t162 * t297 + t173 * t182) * qJD(4) + t248;
t25 = qJD(1) * t249 + t173 * t293 - t345;
t21 = -qJD(3) + (t172 * t383 + t173 * t384) * qJD(4);
t16 = t240 * qJD(4);
t15 = t172 * t329 + qJD(5) * t162 + (-t163 * t309 + t172 * t181) * qJD(4) + (t351 - t398) * qJD(1) + t269;
t14 = t173 * t329 + qJD(5) * t163 + t397 * qJD(1) + (t162 * t309 + t173 * t181) * qJD(4) + t248;
t1 = (t162 * t384 - t163 * t383 + t172 * t397 + t173 * t398) * qJD(4);
t2 = [t454 * t312 + (t474 * t225 + t475 * t226) * qJD(1) + (t14 * (t425 + t427) + t265 * t391 + (t173 * t265 + t247 - t364 + t387) * t15 + (-rSges(6,2) * t368 + t418 + t420) * t25 + (-t345 + t216 + t422) * t24 + ((-t385 + t424) * t25 + (t249 - t425) * t24) * qJD(1)) * m(6) + (t31 * (t170 + t247 + t388) + t30 * (t117 - t171 + t425) + (t30 * pkin(6) + t266 * t31) * t173 + (t162 * t266 - t297 * t335 - t159 - t390 + t414) * t43 + (-t297 * t336 + t43 + (-t113 + t424) * qJD(1) + t308 + t420 + t66) * t44) * m(5) + (t74 * t426 + t73 * t442 + (t299 + t414) * t92 + (t346 + t92 + (-t244 + t254) * qJD(1) + t419) * (t442 * qJD(1) - t216)) * m(4) + (t119 * (t218 + t220 + t253) + t118 * t421 + (t214 + (t324 - t220 + t253) * qJD(1) + t419) * (qJD(1) * t421 - t216)) * m(3) + (t435 + t437) * t335 / 0.2e1 - (t436 - t438 + t439) * t336 / 0.2e1 + (t440 - t487) * t311 + (((-t434 + t458 + t473) * t172 + t417 * t173) * t312 + t423 * qJD(1) + ((t24 * t415 + t300) * t173 + (-t21 * (t384 + t385) + (-t309 - t415) * t25) * t172) * m(6) + (-t432 + t472) * t410 + (t433 + t466) * t409 + (((t453 - t444) * t173 + t434) * t173 + ((t453 - t468) * t172 + t493 - t417 + t457) * t172) * t311) * qJD(4); 0.2e1 * (t14 * t328 + t15 * t327) * m(6) + 0.2e1 * (t30 * t328 + t31 * t327) * m(5) + 0.2e1 * (t327 * t74 + t328 * t73) * m(4) + 0.2e1 * (t118 * t328 + t119 * t327) * m(3); -m(5) * t16 - m(6) * t1; -((((-t352 - t354) * t173 + (t353 + t355) * t172) * t226 + ((t356 + t358) * t173 + (-t357 - t359) * t172) * t225) * qJD(4) + ((-t340 - t342) * t226 + (t341 + t343) * t225) * qJD(1)) * qJD(1) / 0.2e1 + (-t432 * t162 + t433 * t163 + t438 * t172 + t437 * t173) * qJD(1) / 0.2e1 + ((t121 * t336 + t331) * t172 + (t242 + (t246 * t173 + (-t451 - t412) * t172) * qJD(4)) * t173 + (t123 * t336 + t330) * t172 + (t241 + (t245 * t173 + (-t450 - t413) * t172) * qJD(4)) * t173) * t312 + ((t335 * t451 - t331) * t173 + (t242 + (-t412 * t172 + (-t121 + t246) * t173) * qJD(4)) * t172 + (t335 * t450 - t330) * t173 + (t241 + (-t413 * t172 + (-t123 + t245) * t173) * qJD(4)) * t172) * t311 + (t16 * t272 + t42 * t240 - t294 * t182 - (-t162 * t44 + t163 * t43 - t172 * t31 - t173 * t30) * t297 - (-t134 * t43 + t136 * t44) * qJD(1) - (t42 * (t134 * t173 + t136 * t172) - t294 * t298) * qJD(4)) * m(5) - (t436 * qJD(1) + (t459 * t163 + t458 * t162 + (t444 * t162 + t452 * t163) * t173 + (t162 * t496 + t453 * t163 + t456 * t172 - t464) * t172) * t441) * t172 / 0.2e1 + (t435 * qJD(1) + (t434 * t163 + t457 * t162 + (t452 * t162 - t444 * t163 + t464) * t173 + (t453 * t162 - t163 * t496 - t456 * t173) * t172) * t441) * t173 / 0.2e1 + (-(-t24 * t306 + t25 * t305) * qJD(1) - ((t21 * t306 + t25 * t296) * t173 + (t24 * (t296 + t403) + t305 * t21) * t172) * qJD(4) - (-t21 * t384 - t25 * t309) * t162 + (-t21 * t383 - t300) * t163 + (t1 * t384 + t14 * t309 + t25 * t181 + t21 * t398) * t173 + (t15 * t309 + t24 * (pkin(4) * t332 + t181) + t1 * t383 + t21 * t397) * t172) * m(6) + (t439 + t446) * t410 + (t440 + t445) * t409; (-t14 * t172 + t15 * t173 + t163 * t25 + t391 - (t172 * t24 + t173 * t25) * qJD(1)) * m(6);];
tauc = t2(:);
