% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:59
% DurationCPUTime: 17.82s
% Computational Cost: add. (11027->422), mult. (14764->585), div. (0->0), fcn. (14210->6), ass. (0->212)
t480 = Icges(5,4) + Icges(6,4);
t447 = Icges(5,1) + Icges(6,1);
t466 = Icges(5,5) + Icges(6,5);
t462 = Icges(5,2) + Icges(6,2);
t479 = Icges(5,6) + Icges(6,6);
t478 = Icges(5,3) + Icges(6,3);
t241 = pkin(7) + qJ(2);
t238 = sin(t241);
t239 = cos(t241);
t246 = cos(qJ(4));
t243 = cos(pkin(8));
t245 = sin(qJ(4));
t355 = t243 * t245;
t181 = t238 * t355 + t239 * t246;
t492 = t480 * t181;
t354 = t243 * t246;
t357 = t239 * t245;
t182 = t238 * t354 - t357;
t491 = t480 * t182;
t242 = sin(pkin(8));
t362 = t238 * t242;
t443 = t447 * t182 + t466 * t362 - t492;
t445 = -t462 * t181 + t479 * t362 + t491;
t490 = t181 * t445 - t182 * t443;
t448 = -t479 * t181 + t466 * t182 + t478 * t362;
t487 = t480 * t246;
t486 = t480 * t245;
t183 = t238 * t246 - t239 * t355;
t485 = t480 * t183;
t360 = t238 * t245;
t184 = t239 * t354 + t360;
t484 = t480 * t184;
t483 = t183 * t445 + t184 * t443;
t459 = t448 * t362 - t490;
t359 = t239 * t242;
t446 = t479 * t183 + t466 * t184 + t478 * t359;
t444 = t462 * t183 + t479 * t359 + t484;
t442 = t447 * t184 + t466 * t359 + t485;
t477 = (-t466 * t245 - t479 * t246) * t242;
t476 = (-t462 * t246 - t486) * t242;
t475 = (-t447 * t245 - t487) * t242;
t405 = (-t245 * t445 + t246 * t443) * t242 - t448 * t243;
t441 = -t478 * t243 + (-t479 * t245 + t466 * t246) * t242;
t440 = -t479 * t243 + (-t462 * t245 + t487) * t242;
t439 = -t466 * t243 + (t447 * t246 - t486) * t242;
t473 = t448 * t359;
t471 = -t181 * t444 + t182 * t442 + t362 * t446;
t457 = t473 + t483;
t470 = t183 * t444 + t184 * t442 + t359 * t446;
t133 = qJD(2) * t183 - qJD(4) * t182;
t259 = t181 * qJD(4);
t134 = qJD(2) * t184 - t259;
t328 = qJD(2) * t242;
t310 = t239 * t328;
t453 = t479 * t133 + t466 * t134 + t478 * t310;
t131 = qJD(2) * t181 - qJD(4) * t184;
t132 = -qJD(2) * t182 + qJD(4) * t183;
t311 = t238 * t328;
t469 = t479 * t131 + t466 * t132 - t478 * t311;
t452 = t462 * t131 + t480 * t132 - t479 * t311;
t451 = t462 * t133 + t480 * t134 + t479 * t310;
t450 = t480 * t131 + t447 * t132 - t466 * t311;
t449 = t480 * t133 + t447 * t134 + t466 * t310;
t438 = t477 * qJD(4);
t437 = t476 * qJD(4);
t436 = t475 * qJD(4);
t464 = t459 * t238;
t421 = t440 * t181 - t439 * t182 - t441 * t362;
t463 = t440 * t183 + t439 * t184 + t441 * t359;
t461 = -t453 * t243 + (t449 * t246 - t451 * t245 + (-t443 * t245 - t445 * t246) * qJD(4)) * t242;
t460 = -t469 * t243 + (t450 * t246 - t452 * t245 + (-t442 * t245 - t444 * t246) * qJD(4)) * t242;
t455 = t446 * t243 + (t444 * t245 - t442 * t246) * t242;
t454 = -t438 * t243 + (t436 * t246 - t437 * t245 + (-t439 * t245 - t440 * t246) * qJD(4)) * t242;
t435 = (t457 * t238 + t470 * t239) * t242;
t434 = (t471 * t239 + t464) * t242;
t433 = (t466 * t183 - t184 * t479) * t239 + (-t466 * t181 - t182 * t479) * t238;
t432 = t469 * t239;
t229 = -qJD(4) * t243 + qJD(2);
t431 = t463 * t229;
t244 = -qJ(5) - pkin(6);
t383 = pkin(6) + t244;
t237 = pkin(4) * t246 + pkin(3);
t384 = pkin(3) - t237;
t342 = (t383 - rSges(6,3)) * t243 + (rSges(6,1) * t246 - rSges(6,2) * t245 - t384) * t242;
t320 = pkin(4) * t357;
t361 = t238 * t243;
t430 = -rSges(6,1) * t182 + rSges(6,2) * t181 - t237 * t361 + t320;
t429 = t421 * t229;
t411 = (t462 * t184 - t442 - t485) * t239 + (t462 * t182 - t443 + t492) * t238;
t323 = qJD(5) * t242;
t216 = t239 * t323;
t326 = qJD(4) * t242;
t309 = t238 * t326;
t387 = pkin(3) * t243;
t404 = t242 * t383;
t349 = (t387 + t404) * t238 - rSges(6,3) * t362 + t430;
t428 = t229 * t349 + t309 * t342 + t216;
t427 = rSges(6,1) + pkin(4);
t426 = t434 * qJD(4) - t429;
t425 = t435 * qJD(4) + t431;
t330 = qJD(2) * t238;
t424 = (t438 * t239 - t441 * t330) * t242 + t436 * t184 + t437 * t183 + t439 * t132 + t440 * t131;
t329 = qJD(2) * t239;
t423 = (t438 * t238 + t441 * t329) * t242 + t436 * t182 - t437 * t181 + t439 * t134 + t440 * t133;
t422 = t454 * t229;
t416 = t242 * (qJD(4) * t239 * t342 - qJD(5) * t238);
t414 = (t460 * t239 + t461 * t238 + (t455 * t238 + t405 * t239) * qJD(2)) * t242;
t413 = ((t446 * t242 * t329 + t459 * qJD(2) + t444 * t133 + t442 * t134 - t452 * t181 + t450 * t182) * t239 + ((t453 * t238 + t448 * t329 + t432) * t242 + t449 * t182 - t451 * t181 + t443 * t134 + t445 * t133 - t471 * qJD(2)) * t238) * t242;
t412 = (((-t446 * t330 + t432) * t242 + t450 * t184 + t452 * t183 + t442 * t132 + t444 * t131 + t457 * qJD(2)) * t239 + ((t453 * t239 - t448 * t330) * t242 + t449 * t184 + t451 * t183 + t443 * t132 + t445 * t131 - t470 * qJD(2)) * t238) * t242;
t410 = (t447 * t183 - t444 - t484) * t239 + (-t447 * t181 - t445 - t491) * t238;
t409 = t433 * t242;
t408 = t439 + t476;
t407 = -t440 + t475;
t289 = -rSges(5,1) * t182 + rSges(5,2) * t181;
t120 = rSges(5,3) * t362 - t289;
t180 = -rSges(5,3) * t243 + (rSges(5,1) * t246 - rSges(5,2) * t245) * t242;
t406 = -t120 * t229 + t180 * t309;
t403 = t243 * t384;
t385 = pkin(6) * t242;
t292 = t385 + t387;
t185 = t292 * t238;
t210 = t239 * pkin(2) + t238 * qJ(3);
t358 = t239 * t243;
t275 = rSges(4,1) * t358 - rSges(4,2) * t359 + t238 * rSges(4,3);
t400 = t210 + t275;
t379 = pkin(4) * qJD(4);
t317 = t245 * t379;
t301 = t243 * t317;
t316 = t246 * t379;
t399 = t132 * rSges(6,1) + t131 * rSges(6,2) + qJD(2) * t320 + t238 * t316 - t239 * t301 + t244 * t311 + t216;
t223 = pkin(4) * t360;
t356 = t242 * t244;
t398 = t184 * rSges(6,1) + t183 * rSges(6,2) + rSges(6,3) * t359 + t237 * t358 - t239 * t356 + t223;
t395 = t457 - t473;
t390 = t238 / 0.2e1;
t389 = -t239 / 0.2e1;
t386 = pkin(4) * t245;
t382 = -rSges(6,3) * t311 + (t385 + t403) * t330 + t399;
t288 = rSges(6,1) * t134 + rSges(6,2) * t133;
t381 = rSges(6,3) * t310 + t288 + t238 * t323 - pkin(4) * t259 + (t223 + (-t403 - t404) * t239) * qJD(2);
t380 = rSges(6,3) * t242;
t233 = t239 * qJ(3);
t208 = pkin(2) * t238 - t233;
t230 = qJD(3) * t238;
t279 = t230 + (-t185 - t208) * qJD(2);
t46 = t279 + t406;
t372 = t239 * t46;
t269 = (-rSges(5,1) * t245 - rSges(5,2) * t246) * t242;
t194 = qJD(4) * t269;
t321 = qJD(2) * qJD(3);
t332 = qJ(3) * t329 + t230;
t343 = qJD(2) * (-pkin(2) * t330 + t332) + t238 * t321;
t313 = -qJD(2) ^ 2 * t185 + t343;
t346 = t132 * rSges(5,1) + t131 * rSges(5,2);
t72 = -rSges(5,3) * t311 + t346;
t30 = t229 * t72 + (t180 * t330 - t194 * t239) * t326 + t313;
t371 = t30 * t239;
t225 = t239 * t321;
t300 = t180 * t239 * t326;
t231 = qJD(3) * t239;
t158 = qJD(2) * t210 - t231;
t341 = -t292 * t329 - t158;
t290 = -rSges(5,1) * t134 - rSges(5,2) * t133;
t74 = rSges(5,3) * t310 - t290;
t31 = t194 * t309 - t229 * t74 + t225 + (t300 + t341) * qJD(2);
t370 = t31 * t238;
t186 = pkin(3) * t358 + pkin(6) * t359;
t348 = -t186 + t398;
t345 = -rSges(6,2) * t182 - t427 * t181;
t344 = -rSges(6,2) * t184 + t427 * t183;
t336 = t186 + t210;
t268 = (-rSges(6,1) * t245 - rSges(6,2) * t246) * t242;
t322 = qJD(5) * t243;
t335 = qJD(4) * t268 - t242 * t317 - t322;
t334 = rSges(4,2) * t311 + rSges(4,3) * t329;
t333 = rSges(4,2) * t362 + t239 * rSges(4,3);
t331 = -qJD(2) * t208 + t230;
t319 = rSges(4,1) * t361;
t124 = t184 * rSges(5,1) + t183 * rSges(5,2) + rSges(5,3) * t359;
t312 = -qJD(2) * t185 + t331;
t306 = -rSges(4,1) * t243 - pkin(2);
t305 = -t326 / 0.2e1;
t304 = t326 / 0.2e1;
t303 = -t237 * t243 - pkin(2);
t297 = t238 * t305;
t296 = t238 * t304;
t295 = t239 * t305;
t294 = t239 * t304;
t278 = qJD(2) * t336 - t231;
t47 = t124 * t229 + t278 - t300;
t285 = t238 * t46 - t239 * t47;
t284 = -t238 * t72 + t239 * t74;
t283 = t120 * t239 - t124 * t238;
t277 = qJD(2) * t297;
t276 = qJD(2) * t294;
t260 = -t387 - pkin(2) + (-rSges(5,3) - pkin(6)) * t242;
t154 = t319 - t333;
t150 = rSges(5,1) * t183 - rSges(5,2) * t184;
t148 = -rSges(5,1) * t181 - rSges(5,2) * t182;
t126 = qJD(2) * t400 - t231;
t125 = t230 + (-t154 - t208) * qJD(2);
t90 = t225 + (-qJD(2) * t275 - t158) * qJD(2);
t89 = qJD(2) * (-qJD(2) * t319 + t334) + t343;
t50 = t283 * t326 + qJD(1);
t29 = t348 * t229 + t278 - t416;
t28 = t279 + t428;
t25 = -t322 + qJD(1) + (-t238 * t348 - t239 * t349) * t326;
t20 = ((-t120 * t238 - t124 * t239) * qJD(2) + t284) * t326;
t19 = t225 - t381 * t229 + t335 * t309 + (t341 + t416) * qJD(2);
t18 = t382 * t229 + (qJD(5) * t329 + (-t239 * t335 + t330 * t342) * qJD(4)) * t242 + t313;
t1 = (t381 * t239 - t382 * t238 + (t238 * t349 - t239 * t348) * qJD(2)) * t326;
t2 = [m(5) * t20 + m(6) * t1; (((t459 + t470 + t490) * t239 + t395 * t238) * t326 + t431) * t297 + (t19 * (t233 + t430) + t28 * (t239 * t316 + t231 - t288) + t18 * (t210 + t398) + t29 * (t332 + t399) + (t19 * (-pkin(2) + t356 - t380) + t28 * (t301 - t323)) * t238 + (t28 * ((-rSges(6,3) + t244) * t242 + t303) * t239 + (t28 * (-qJ(3) - t386) + t29 * (t303 - t380)) * t238) * qJD(2) - (-t28 + t312 + t428) * t29) * m(6) + (t31 * (t233 + t289) + t46 * (t231 + t290) + t30 * (t124 + t336) + t47 * (t332 + t346) + t260 * t370 + (t260 * t372 + (-t46 * qJ(3) + t47 * (-rSges(5,3) * t242 - pkin(2) - t292)) * t238) * qJD(2) - (t312 + t406 - t46) * t47) * m(5) + (t90 * (t306 * t238 + t233 + t333) + t125 * t231 + t89 * t400 + t126 * (t332 + t334) + (t125 * (rSges(4,2) * t242 + t306) * t239 + (t125 * (-rSges(4,3) - qJ(3)) + t126 * t306) * t238) * qJD(2) - (-qJD(2) * t154 - t125 + t331) * t126) * m(4) + (t424 + t460) * t294 + (t463 - t455) * t277 + (t405 - t421) * t276 + (t423 + t425 + t461) * t296 + (((t395 - t471 - t483) * t239 - t464) * t326 + t426 + t429) * t295 + t422; 0.2e1 * (t18 * t389 + t19 * t390) * m(6) + 0.2e1 * (t370 / 0.2e1 - t371 / 0.2e1) * m(5) + 0.2e1 * (t389 * t89 + t390 * t90) * m(4); -(((-t408 * t245 + t407 * t246) * t229 + ((t411 * t245 + t410 * t246) * t242 - t433 * t243) * qJD(4)) * t242 - t477 * t229 * t243) * t229 / 0.2e1 + (-t454 * t243 + t414) * t229 / 0.2e1 - (t414 * qJD(4) + t422) * t243 / 0.2e1 + (t413 * qJD(4) + t423 * t229) * t362 / 0.2e1 - t425 * t311 / 0.2e1 + ((t411 * t181 + t410 * t182 + t409 * t238) * t326 + (-t408 * t181 + t407 * t182 + t362 * t477) * t229) * t297 + (-t423 * t243 + t413) * t296 + ((-t183 * t411 + t410 * t184 + t409 * t239) * t326 + (t408 * t183 + t407 * t184 + t359 * t477) * t229) * t295 + (-t424 * t243 + t412) * t294 + (-t243 * t463 + t435) * t277 + (t421 * t243 + t434) * t276 + ((-t18 * t348 - t19 * t349 + t28 * t381 - t29 * t382) * t243 + ((-t18 * t342 - t29 * t335 - t1 * t349 + t25 * t381 + (-t25 * t348 + t28 * t342) * qJD(2)) * t239 + (t19 * t342 + t28 * t335 - t1 * t348 - t25 * t382 + (t25 * t349 + t29 * t342) * qJD(2)) * t238) * t242 - (-t28 * t345 + t29 * t344) * t229 - (t25 * (-t238 * t344 + t239 * t345) + (-t242 * t386 + t268) * (t238 * t28 - t239 * t29)) * t326) * m(6) + ((t120 * t31 - t124 * t30 + t46 * t74 - t47 * t72) * t243 + (t20 * t283 + t50 * (-t120 * t330 - t124 * t329 + t284) + t285 * t194 + (t370 - t371 + (t238 * t47 + t372) * qJD(2)) * t180) * t242 - (-t148 * t46 + t150 * t47) * t229 - (t50 * (t148 * t239 - t150 * t238) + t285 * t269) * t326) * m(5) + (t426 * qJD(2) + t412 * qJD(4) + t424 * t229) * t359 / 0.2e1; (-t1 * t243 + (t18 * t238 + t19 * t239) * t242) * m(6);];
tauc = t2(:);
