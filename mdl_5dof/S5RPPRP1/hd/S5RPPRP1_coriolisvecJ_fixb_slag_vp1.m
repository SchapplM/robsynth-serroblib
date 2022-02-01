% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:13
% EndTime: 2022-01-23 09:12:38
% DurationCPUTime: 21.00s
% Computational Cost: add. (11133->445), mult. (15000->608), div. (0->0), fcn. (14318->8), ass. (0->225)
t505 = Icges(5,4) + Icges(6,4);
t463 = Icges(5,1) + Icges(6,1);
t496 = Icges(5,5) + Icges(6,5);
t483 = Icges(5,2) + Icges(6,2);
t504 = Icges(5,6) + Icges(6,6);
t503 = Icges(5,3) + Icges(6,3);
t246 = qJ(1) + pkin(7);
t242 = sin(t246);
t243 = cos(t246);
t252 = cos(qJ(4));
t248 = cos(pkin(8));
t250 = sin(qJ(4));
t371 = t248 * t250;
t183 = t242 * t371 + t243 * t252;
t516 = t505 * t183;
t370 = t248 * t252;
t373 = t243 * t250;
t184 = t242 * t370 - t373;
t515 = t505 * t184;
t247 = sin(pkin(8));
t378 = t242 * t247;
t460 = -t483 * t183 + t504 * t378 + t515;
t458 = t463 * t184 + t496 * t378 - t516;
t462 = -t504 * t183 + t496 * t184 + t503 * t378;
t185 = t242 * t252 - t243 * t371;
t514 = t505 * t185;
t513 = t505 * t252;
t512 = t505 * t250;
t376 = t242 * t250;
t186 = t243 * t370 + t376;
t511 = t505 * t186;
t510 = t460 * t183 - t458 * t184;
t476 = t462 * t378 - t510;
t375 = t243 * t247;
t461 = t504 * t185 + t496 * t186 + t503 * t375;
t459 = t483 * t185 + t504 * t375 + t511;
t457 = t463 * t186 + t496 * t375 + t514;
t494 = -t503 * t248 + (-t504 * t250 + t496 * t252) * t247;
t493 = -t504 * t248 + (-t483 * t250 + t513) * t247;
t492 = -t496 * t248 + (t463 * t252 - t512) * t247;
t509 = t462 * t375;
t508 = t460 * t185 + t458 * t186;
t507 = -t459 * t183 + t457 * t184 + t461 * t378;
t485 = t508 + t509;
t506 = t459 * t185 + t457 * t186 + t461 * t375;
t502 = (-t496 * t250 - t504 * t252) * t247;
t501 = (-t483 * t252 - t512) * t247;
t500 = (-t463 * t250 - t513) * t247;
t499 = t476 * t242;
t470 = t493 * t183 - t492 * t184 - t378 * t494;
t498 = t493 * t185 + t492 * t186 + t375 * t494;
t473 = -(t460 * t250 - t458 * t252) * t247 - t462 * t248;
t133 = qJD(1) * t185 - qJD(4) * t184;
t266 = t183 * qJD(4);
t134 = qJD(1) * t186 - t266;
t346 = qJD(1) * t247;
t327 = t243 * t346;
t468 = t504 * t133 + t496 * t134 + t503 * t327;
t131 = qJD(1) * t183 - qJD(4) * t186;
t132 = -qJD(1) * t184 + qJD(4) * t185;
t328 = t242 * t346;
t497 = t504 * t131 + t496 * t132 - t503 * t328;
t467 = t483 * t131 + t505 * t132 - t504 * t328;
t466 = t483 * t133 + t505 * t134 + t504 * t327;
t465 = t505 * t131 + t463 * t132 - t496 * t328;
t464 = t505 * t133 + t463 * t134 + t496 * t327;
t491 = t502 * qJD(4);
t490 = t501 * qJD(4);
t489 = t500 * qJD(4);
t488 = (t485 * t242 + t506 * t243) * t247;
t487 = (t507 * t243 + t499) * t247;
t232 = -qJD(4) * t248 + qJD(1);
t484 = t498 * t232;
t482 = t470 * t232;
t480 = t487 * qJD(4) - t482;
t479 = t488 * qJD(4) + t484;
t478 = -t468 * t248 + (t464 * t252 - t466 * t250 + (-t458 * t250 - t460 * t252) * qJD(4)) * t247;
t477 = -t497 * t248 + (t465 * t252 - t467 * t250 + (-t457 * t250 - t459 * t252) * qJD(4)) * t247;
t348 = qJD(1) * t242;
t444 = (t491 * t243 - t494 * t348) * t247 + t489 * t186 + t490 * t185 + t492 * t132 + t493 * t131;
t347 = qJD(1) * t243;
t443 = (t491 * t242 + t494 * t347) * t247 + t489 * t184 - t490 * t183 + t492 * t134 + t493 * t133;
t472 = t461 * t248 + (t459 * t250 - t457 * t252) * t247;
t471 = -t491 * t248 + (t489 * t252 - t490 * t250 + (-t492 * t250 - t493 * t252) * qJD(4)) * t247;
t456 = (t496 * t185 - t186 * t504) * t243 + (-t496 * t183 - t184 * t504) * t242;
t455 = t485 - t509;
t454 = t497 * t243;
t249 = -qJ(5) - pkin(6);
t400 = pkin(6) + t249;
t241 = pkin(4) * t252 + pkin(3);
t401 = pkin(3) - t241;
t359 = (t400 - rSges(6,3)) * t248 + (rSges(6,1) * t252 - rSges(6,2) * t250 - t401) * t247;
t338 = pkin(4) * t373;
t377 = t242 * t248;
t453 = -t184 * rSges(6,1) + t183 * rSges(6,2) - t241 * t377 + t338;
t430 = (t483 * t186 - t457 - t514) * t243 + (t483 * t184 - t458 + t516) * t242;
t341 = qJD(5) * t247;
t217 = t243 * t341;
t344 = qJD(4) * t247;
t326 = t242 * t344;
t404 = pkin(3) * t248;
t421 = t247 * t400;
t365 = (t404 + t421) * t242 - rSges(6,3) * t378 + t453;
t452 = t232 * t365 + t326 * t359 + t217;
t254 = qJD(1) ^ 2;
t445 = rSges(6,1) + pkin(4);
t442 = t471 * t232;
t436 = t247 * (qJD(4) * t243 * t359 - qJD(5) * t242);
t374 = t243 * t248;
t188 = pkin(3) * t374 + pkin(6) * t375;
t211 = t243 * pkin(2) + t242 * qJ(3);
t253 = cos(qJ(1));
t244 = t253 * pkin(1);
t419 = t244 + t211;
t435 = t188 + t419;
t433 = (t477 * t243 + t478 * t242 + (t242 * t472 + t243 * t473) * qJD(1)) * t247;
t432 = ((t247 * t347 * t461 + qJD(1) * t476 + t133 * t459 + t134 * t457 - t183 * t467 + t184 * t465) * t243 + ((t242 * t468 + t347 * t462 + t454) * t247 + t464 * t184 - t466 * t183 + t458 * t134 + t460 * t133 - t507 * qJD(1)) * t242) * t247;
t431 = (((-t348 * t461 + t454) * t247 + t465 * t186 + t467 * t185 + t457 * t132 + t459 * t131 + t485 * qJD(1)) * t243 + ((t243 * t468 - t348 * t462) * t247 + t464 * t186 + t466 * t185 + t458 * t132 + t460 * t131 - t506 * qJD(1)) * t242) * t247;
t429 = (t185 * t463 - t459 - t511) * t243 + (-t183 * t463 - t460 - t515) * t242;
t428 = t456 * t247;
t287 = rSges(4,1) * t374 - rSges(4,2) * t375 + t242 * rSges(4,3);
t427 = t287 + t419;
t426 = t492 + t501;
t425 = -t493 + t500;
t300 = t184 * rSges(5,1) - t183 * rSges(5,2);
t122 = rSges(5,3) * t378 + t300;
t182 = -rSges(5,3) * t248 + (rSges(5,1) * t252 - rSges(5,2) * t250) * t247;
t424 = -t122 * t232 + t182 * t326;
t420 = t248 * t401;
t402 = pkin(6) * t247;
t304 = t402 + t404;
t187 = t304 * t242;
t316 = t243 * rSges(3,1) - rSges(3,2) * t242;
t418 = t244 + t316;
t395 = pkin(4) * qJD(4);
t333 = t250 * t395;
t314 = t248 * t333;
t332 = t252 * t395;
t417 = t132 * rSges(6,1) + t131 * rSges(6,2) + qJD(1) * t338 + t242 * t332 - t243 * t314 + t249 * t328 + t217;
t226 = pkin(4) * t376;
t372 = t247 * t249;
t416 = t186 * rSges(6,1) + t185 * rSges(6,2) + rSges(6,3) * t375 + t241 * t374 - t243 * t372 + t226;
t408 = t242 / 0.2e1;
t407 = -t243 / 0.2e1;
t251 = sin(qJ(1));
t405 = pkin(1) * t251;
t403 = pkin(4) * t250;
t399 = -rSges(6,3) * t328 + (t402 + t420) * t348 + t417;
t299 = t134 * rSges(6,1) + t133 * rSges(6,2);
t398 = rSges(6,3) * t327 + t299 + t242 * t341 - pkin(4) * t266 + (t226 + (-t420 - t421) * t243) * qJD(1);
t396 = rSges(6,3) * t247;
t233 = qJD(3) * t242;
t236 = t243 * qJ(3);
t209 = pkin(2) * t242 - t236;
t320 = -t209 - t405;
t268 = t233 + (-t187 + t320) * qJD(1);
t46 = t268 + t424;
t388 = t243 * t46;
t280 = (-rSges(5,1) * t250 - rSges(5,2) * t252) * t247;
t196 = qJD(4) * t280;
t337 = t254 * t405;
t339 = qJD(1) * qJD(3);
t350 = qJ(3) * t347 + t233;
t290 = qJD(1) * (-pkin(2) * t348 + t350) + t242 * t339 - t337;
t281 = -t187 * t254 + t290;
t362 = t132 * rSges(5,1) + t131 * rSges(5,2);
t72 = -rSges(5,3) * t328 + t362;
t30 = t232 * t72 + (t182 * t348 - t196 * t243) * t344 + t281;
t387 = t30 * t243;
t336 = t254 * t244;
t306 = t243 * t339 - t336;
t313 = t182 * t243 * t344;
t234 = qJD(3) * t243;
t158 = qJD(1) * t211 - t234;
t358 = -t304 * t347 - t158;
t301 = -t134 * rSges(5,1) - t133 * rSges(5,2);
t74 = rSges(5,3) * t327 - t301;
t31 = t196 * t326 - t232 * t74 + (t313 + t358) * qJD(1) + t306;
t386 = t31 * t242;
t364 = -t188 + t416;
t361 = -rSges(6,2) * t184 - t183 * t445;
t360 = -rSges(6,2) * t186 + t185 * t445;
t279 = (-rSges(6,1) * t250 - rSges(6,2) * t252) * t247;
t340 = qJD(5) * t248;
t353 = qJD(4) * t279 - t247 * t333 - t340;
t352 = rSges(4,2) * t328 + rSges(4,3) * t347;
t351 = rSges(4,2) * t378 + t243 * rSges(4,3);
t349 = -qJD(1) * t209 + t233;
t335 = rSges(4,1) * t377;
t126 = t186 * rSges(5,1) + t185 * rSges(5,2) + rSges(5,3) * t375;
t323 = -rSges(4,1) * t248 - pkin(2);
t322 = -t344 / 0.2e1;
t321 = t344 / 0.2e1;
t318 = t236 - t405;
t317 = -t241 * t248 - pkin(2);
t310 = t242 * t322;
t309 = t242 * t321;
t308 = t243 * t322;
t307 = t243 * t321;
t210 = rSges(3,1) * t242 + rSges(3,2) * t243;
t267 = qJD(1) * t435 - t234;
t47 = t126 * t232 + t267 - t313;
t296 = t242 * t46 - t243 * t47;
t295 = -t242 * t72 + t243 * t74;
t294 = t122 * t243 - t126 * t242;
t289 = qJD(1) * t310;
t288 = qJD(1) * t307;
t274 = t349 + (-t187 - t405) * qJD(1);
t270 = -t404 - pkin(2) + (-rSges(5,3) - pkin(6)) * t247;
t154 = t335 - t351;
t150 = rSges(5,1) * t185 - rSges(5,2) * t186;
t148 = -rSges(5,1) * t183 - rSges(5,2) * t184;
t100 = qJD(1) * t427 - t234;
t99 = t233 + (-t154 + t320) * qJD(1);
t90 = (-qJD(1) * t287 - t158) * qJD(1) + t306;
t89 = qJD(1) * (-qJD(1) * t335 + t352) + t290;
t50 = t294 * t344 + qJD(2);
t29 = t364 * t232 + t267 - t436;
t28 = t268 + t452;
t25 = -t340 + qJD(2) + (-t242 * t364 - t243 * t365) * t344;
t20 = ((-t122 * t242 - t126 * t243) * qJD(1) + t295) * t344;
t19 = -t398 * t232 + t353 * t326 + (t358 + t436) * qJD(1) + t306;
t18 = t399 * t232 + (qJD(5) * t347 + (-t243 * t353 + t348 * t359) * qJD(4)) * t247 + t281;
t1 = (t398 * t243 - t399 * t242 + (t242 * t365 - t243 * t364) * qJD(1)) * t344;
t2 = [m(3) * ((-t210 * t254 - t337) * t418 + (-t336 + (-0.2e1 * t316 - t244 + t418) * t254) * (-t210 - t405)) + t442 + (((t506 + t476 + t510) * t243 + t455 * t242) * t344 + t484) * t310 + (t19 * (t318 + t453) + t28 * (t243 * t332 + t234 - t299) + t18 * (t419 + t416) + t29 * (t350 + t417) + (t19 * (-pkin(2) + t372 - t396) + t28 * (t314 - t341)) * t242 + ((-t251 * t29 - t253 * t28) * pkin(1) + t28 * ((-rSges(6,3) + t249) * t247 + t317) * t243 + (t28 * (-qJ(3) - t403) + t29 * (t317 - t396)) * t242) * qJD(1) - (-t28 + t274 + t452) * t29) * m(6) + (t31 * (-t300 + t318) + t46 * (t234 + t301) + t30 * (t126 + t435) + t47 * (t350 + t362) + t270 * t386 + ((-t251 * t47 - t253 * t46) * pkin(1) + t270 * t388 + (-t46 * qJ(3) + t47 * (-rSges(5,3) * t247 - pkin(2) - t304)) * t242) * qJD(1) - (t274 + t424 - t46) * t47) * m(5) + (t90 * (t242 * t323 + t318 + t351) + t99 * t234 + t89 * t427 + t100 * (t350 + t352) + ((-t100 * t251 - t253 * t99) * pkin(1) + t99 * (rSges(4,2) * t247 + t323) * t243 + (t99 * (-rSges(4,3) - qJ(3)) + t100 * t323) * t242) * qJD(1) - (-t99 + (-t154 - t405) * qJD(1) + t349) * t100) * m(4) + (t444 + t477) * t307 + (t498 - t472) * t289 + (-t470 + t473) * t288 + (t443 + t478 + t479) * t309 + (((t455 - t507 - t508) * t243 - t499) * t344 + t480 + t482) * t308; m(5) * t20 + m(6) * t1; 0.2e1 * (t18 * t407 + t19 * t408) * m(6) + 0.2e1 * (t386 / 0.2e1 - t387 / 0.2e1) * m(5) + 0.2e1 * (t407 * t89 + t408 * t90) * m(4); -(((-t250 * t426 + t252 * t425) * t232 + ((t250 * t430 + t252 * t429) * t247 - t456 * t248) * qJD(4)) * t247 - t502 * t232 * t248) * t232 / 0.2e1 + (-t248 * t471 + t433) * t232 / 0.2e1 - (t433 * qJD(4) + t442) * t248 / 0.2e1 + (t432 * qJD(4) + t443 * t232) * t378 / 0.2e1 - t479 * t328 / 0.2e1 + ((t183 * t430 + t184 * t429 + t242 * t428) * t344 + (-t183 * t426 + t184 * t425 + t378 * t502) * t232) * t310 + (-t248 * t443 + t432) * t309 + ((-t185 * t430 + t429 * t186 + t428 * t243) * t344 + (t185 * t426 + t186 * t425 + t375 * t502) * t232) * t308 + (-t248 * t444 + t431) * t307 + (-t248 * t498 + t488) * t289 + (t248 * t470 + t487) * t288 + (-(-t28 * t361 + t29 * t360) * t232 - (t25 * (-t242 * t360 + t243 * t361) + (-t247 * t403 + t279) * (t242 * t28 - t243 * t29)) * t344 + (-t18 * t364 - t19 * t365 + t28 * t398 - t29 * t399) * t248 + ((-t18 * t359 - t29 * t353 - t1 * t365 + t25 * t398 + (-t25 * t364 + t28 * t359) * qJD(1)) * t243 + (t19 * t359 + t28 * t353 - t1 * t364 - t25 * t399 + (t25 * t365 + t29 * t359) * qJD(1)) * t242) * t247) * m(6) + ((t122 * t31 - t126 * t30 + t46 * t74 - t47 * t72) * t248 + (t20 * t294 + t50 * (-t122 * t348 - t126 * t347 + t295) + t296 * t196 + (t386 - t387 + (t242 * t47 + t388) * qJD(1)) * t182) * t247 - (-t148 * t46 + t150 * t47) * t232 - (t50 * (t148 * t243 - t150 * t242) + t296 * t280) * t344) * m(5) + (qJD(1) * t480 + t431 * qJD(4) + t444 * t232) * t375 / 0.2e1; (-t1 * t248 + (t18 * t242 + t19 * t243) * t247) * m(6);];
tauc = t2(:);
