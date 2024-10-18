% Calculate time derivative of joint inertia matrix for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR12_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR12_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR12_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:21
% EndTime: 2024-09-28 18:07:34
% DurationCPUTime: 7.42s
% Computational Cost: add. (73605->769), mult. (94984->1130), div. (0->0), fcn. (96906->26), ass. (0->445)
t368 = cos(pkin(11));
t370 = cos(pkin(5));
t462 = t368 * t370;
t365 = sin(pkin(11));
t366 = sin(pkin(6));
t466 = t366 * t368;
t433 = pkin(10) * t466;
t316 = -t365 * pkin(4) + t370 * t433;
t504 = pkin(10) * t366;
t434 = t365 * t504;
t318 = pkin(4) * t462 + t434;
t372 = sin(qJ(4));
t376 = cos(qJ(4));
t398 = t316 * t372 + t318 * t376;
t220 = pkin(3) * t462 + t398;
t399 = t316 * t376 - t318 * t372;
t231 = pkin(3) * t365 - t399;
t373 = sin(qJ(3));
t377 = cos(qJ(3));
t515 = t220 * t377 - t231 * t373;
t139 = pkin(2) * t462 + t515;
t514 = -t220 * t373 - t231 * t377;
t143 = pkin(2) * t365 - t514;
t245 = t399 * qJD(4);
t374 = sin(qJ(2));
t378 = cos(qJ(2));
t351 = t377 * pkin(3) + pkin(2);
t396 = t373 * t374 - t377 * t378;
t439 = qJD(2) * t374;
t380 = (-qJD(3) * t396 - t373 * t439) * pkin(3);
t438 = qJD(2) * t378;
t212 = (t351 * t438 + t380) * t370;
t364 = qJ(2) + qJ(3);
t356 = sin(t364);
t362 = qJD(2) + qJD(3);
t314 = -pkin(3) * t356 * t362 - pkin(2) * t439;
t403 = t368 * t212 + t365 * t314;
t436 = qJD(4) * t377;
t437 = qJD(4) * t373;
t371 = sin(qJ(5));
t369 = cos(pkin(6));
t375 = cos(qJ(5));
t458 = t370 * t375;
t425 = t369 * t458;
t292 = -t365 * t371 + t368 * t425;
t459 = t370 * t371;
t460 = t369 * t375;
t294 = -t365 * t460 - t368 * t459;
t359 = qJ(4) + t364;
t349 = sin(t359);
t350 = cos(t359);
t353 = qJD(4) + t362;
t426 = t369 * t459;
t291 = t365 * t375 + t368 * t426;
t461 = t369 * t371;
t386 = t365 * t461 - t368 * t458;
t367 = sin(pkin(5));
t467 = t366 * t367;
t428 = t371 * t467;
t513 = -t291 * t350 + t349 * t386 + t368 * t428;
t131 = (-t292 * t349 + t294 * t350) * t353 + t513 * qJD(5);
t427 = t375 * t467;
t202 = t292 * t350 + t294 * t349 - t368 * t427;
t132 = (-t291 * t349 - t350 * t386) * t353 + t202 * qJD(5);
t475 = t353 * t366;
t480 = t350 * t365;
t253 = (t349 * t462 + t480) * t475;
t91 = rSges(6,1) * t132 + rSges(6,2) * t131 + rSges(6,3) * t253;
t523 = (qJD(3) * t515 + t245 * t373 + t398 * t436) * t378 - t143 * t439 + (qJD(3) * t514 + t245 * t377 - t398 * t437) * t374 + t139 * t438 - t403 + t91;
t468 = t365 * t370;
t315 = pkin(4) * t368 + t370 * t434;
t317 = pkin(4) * t468 - t433;
t400 = t315 * t372 + t317 * t376;
t219 = -pkin(3) * t468 - t400;
t401 = -t315 * t376 + t317 * t372;
t228 = pkin(3) * t368 - t401;
t517 = t219 * t377 - t228 * t373;
t138 = -pkin(2) * t468 + t517;
t516 = -t219 * t373 - t228 * t377;
t142 = pkin(2) * t368 - t516;
t244 = t401 * qJD(4);
t402 = -t365 * t212 + t368 * t314;
t290 = -t365 * t425 - t368 * t371;
t296 = t365 * t459 - t368 * t460;
t385 = t365 * t426 - t368 * t375;
t387 = t365 * t458 + t368 * t461;
t512 = t349 * t387 + t350 * t385 - t365 * t428;
t129 = (-t290 * t349 + t296 * t350) * t353 + t512 * qJD(5);
t200 = t290 * t350 + t296 * t349 + t365 * t427;
t130 = (t349 * t385 - t350 * t387) * t353 + t200 * qJD(5);
t479 = t350 * t368;
t252 = (-t349 * t468 + t479) * t475;
t90 = rSges(6,1) * t130 + rSges(6,2) * t129 + rSges(6,3) * t252;
t522 = (qJD(3) * t517 + t244 * t373 - t400 * t436) * t378 - t142 * t439 + (qJD(3) * t516 + t244 * t377 + t400 * t437) * t374 + t138 * t438 - t402 + t90;
t463 = t367 * t369;
t465 = t366 * t370;
t254 = t349 * t466 + (t350 * t465 + t463) * t365;
t121 = -rSges(6,1) * t512 + rSges(6,2) * t200 + rSges(6,3) * t254;
t379 = pkin(7) + pkin(8);
t363 = pkin(9) + t379;
t503 = t373 * pkin(3);
t435 = t378 * t503;
t267 = -t363 * t367 + (t351 * t374 + t435) * t370;
t335 = pkin(10) * t369 + t363;
t414 = t335 * t367 + t267;
t352 = t378 * pkin(2) + pkin(1);
t357 = cos(t364);
t330 = pkin(3) * t357 + t352;
t501 = pkin(1) - t330;
t521 = t138 * t374 + t142 * t378 + t414 * t365 + t501 * t368 + t121;
t481 = t349 * t365;
t255 = -t368 * t463 + (-t350 * t462 + t481) * t366;
t122 = -rSges(6,1) * t513 + rSges(6,2) * t202 + rSges(6,3) * t255;
t520 = t139 * t374 + t143 * t378 + t501 * t365 - t414 * t368 + t122;
t354 = pkin(5) + t364;
t347 = qJ(4) + t354;
t341 = cos(t347);
t507 = t353 / 0.2e1;
t325 = t341 * t507;
t355 = pkin(5) - t364;
t348 = -qJ(4) + t355;
t342 = cos(t348);
t477 = t353 * t342;
t269 = t325 + t477 / 0.2e1;
t476 = t353 * t365;
t213 = t269 * t368 - t349 * t476;
t340 = sin(t348);
t326 = t340 * t507;
t339 = sin(t347);
t478 = t353 * t339;
t271 = t326 - t478 / 0.2e1;
t214 = t271 * t368 - t350 * t476;
t152 = Icges(5,5) * t213 + Icges(5,6) * t214;
t464 = t367 * t368;
t519 = t152 * t464;
t268 = t325 - t477 / 0.2e1;
t270 = t326 + t478 / 0.2e1;
t203 = Icges(5,5) * t270 + Icges(5,6) * t268;
t518 = t203 * t370;
t500 = -pkin(2) + t351;
t510 = 2 * m(4);
t509 = 2 * m(5);
t508 = 2 * m(6);
t506 = t362 / 0.2e1;
t505 = pkin(2) * t374;
t502 = -pkin(1) + t352;
t115 = -Icges(6,5) * t512 + Icges(6,6) * t200 + Icges(6,3) * t254;
t117 = -Icges(6,4) * t512 + Icges(6,2) * t200 + Icges(6,6) * t254;
t119 = -Icges(6,1) * t512 + Icges(6,4) * t200 + Icges(6,5) * t254;
t84 = Icges(6,5) * t130 + Icges(6,6) * t129 + Icges(6,3) * t252;
t86 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t252;
t88 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t252;
t20 = t115 * t253 + t117 * t131 + t119 * t132 + t202 * t86 + t255 * t84 - t513 * t88;
t116 = -Icges(6,5) * t513 + Icges(6,6) * t202 + Icges(6,3) * t255;
t118 = -Icges(6,4) * t513 + Icges(6,2) * t202 + Icges(6,6) * t255;
t120 = -Icges(6,1) * t513 + Icges(6,4) * t202 + Icges(6,5) * t255;
t85 = Icges(6,5) * t132 + Icges(6,6) * t131 + Icges(6,3) * t253;
t87 = Icges(6,4) * t132 + Icges(6,2) * t131 + Icges(6,6) * t253;
t89 = Icges(6,1) * t132 + Icges(6,4) * t131 + Icges(6,5) * t253;
t21 = t116 * t253 + t118 * t131 + t120 * t132 + t202 * t87 + t255 * t85 - t513 * t89;
t392 = t349 * t375 + t350 * t461;
t423 = qJD(5) * t465;
t179 = -t371 * t423 + ((-t349 * t460 - t350 * t371) * t353 - t392 * qJD(5)) * t367;
t393 = -t349 * t371 + t350 * t460;
t180 = t375 * t423 + ((-t349 * t461 + t350 * t375) * t353 + t393 * qJD(5)) * t367;
t412 = t349 * t353 * t467;
t105 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t412;
t106 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t412;
t107 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t412;
t259 = t366 * t458 + t367 * t393;
t260 = t366 * t459 + t367 * t392;
t288 = -t350 * t467 + t369 * t370;
t145 = Icges(6,5) * t260 + Icges(6,6) * t259 + Icges(6,3) * t288;
t146 = Icges(6,4) * t260 + Icges(6,2) * t259 + Icges(6,6) * t288;
t147 = Icges(6,1) * t260 + Icges(6,4) * t259 + Icges(6,5) * t288;
t39 = t105 * t255 + t106 * t202 - t107 * t513 + t131 * t146 + t132 * t147 + t145 * t253;
t11 = t370 * t39 + (t20 * t365 - t21 * t368) * t367;
t474 = t353 * t368;
t215 = -t269 * t365 - t349 * t474;
t216 = -t271 * t365 - t350 * t474;
t153 = Icges(5,5) * t215 + Icges(5,6) * t216;
t154 = Icges(5,4) * t213 + Icges(5,2) * t214;
t155 = Icges(5,4) * t215 + Icges(5,2) * t216;
t156 = Icges(5,1) * t213 + Icges(5,4) * t214;
t157 = Icges(5,1) * t215 + Icges(5,4) * t216;
t334 = t342 / 0.2e1;
t308 = t334 + t341 / 0.2e1;
t263 = t308 * t368 - t481;
t333 = t339 / 0.2e1;
t307 = t333 - t340 / 0.2e1;
t264 = t307 * t368 + t480;
t184 = Icges(5,4) * t264 + Icges(5,2) * t263 - Icges(5,6) * t464;
t265 = -t308 * t365 - t349 * t368;
t266 = -t307 * t365 + t479;
t469 = t365 * t367;
t185 = Icges(5,4) * t266 + Icges(5,2) * t265 + Icges(5,6) * t469;
t186 = Icges(5,1) * t264 + Icges(5,4) * t263 - Icges(5,5) * t464;
t187 = Icges(5,1) * t266 + Icges(5,4) * t265 + Icges(5,5) * t469;
t204 = Icges(5,4) * t270 + Icges(5,2) * t268;
t205 = Icges(5,1) * t270 + Icges(5,4) * t268;
t306 = t333 + t340 / 0.2e1;
t309 = t334 - t341 / 0.2e1;
t232 = Icges(5,4) * t309 + Icges(5,2) * t306 + Icges(5,6) * t370;
t233 = Icges(5,1) * t309 + Icges(5,4) * t306 + Icges(5,5) * t370;
t30 = (-t153 * t464 + t155 * t263 + t157 * t264 + t185 * t214 + t187 * t213) * t469 - (t154 * t263 + t156 * t264 + t184 * t214 + t186 * t213 - t519) * t464 + (-t203 * t464 + t204 * t263 + t205 * t264 + t213 * t233 + t214 * t232) * t370;
t499 = -t11 - t30;
t498 = t522 * t370;
t496 = pkin(2) * qJD(2);
t110 = rSges(6,1) * t180 + rSges(6,2) * t179 + rSges(6,3) * t412;
t394 = -pkin(4) * t376 - t372 * t504;
t328 = pkin(3) - t394;
t329 = -pkin(4) * t372 + t376 * t504;
t262 = -t328 * t373 + t329 * t377;
t319 = t394 * qJD(4);
t320 = t329 * qJD(4);
t397 = -t328 * t377 - t329 * t373;
t421 = t262 + t503;
t444 = -t397 - t500;
t495 = -t110 - ((t319 * t373 + t320 * t377) * t374 - (t319 * t377 - t320 * t373) * t378 + (pkin(3) * t396 + t262 * t374 - t378 * t397) * qJD(3) + (t374 * t421 + t378 * t444) * qJD(2)) * t367;
t494 = t521 * t370;
t492 = Icges(3,4) * t374;
t491 = Icges(3,4) * t378;
t343 = sin(t354);
t473 = t362 * t343;
t346 = cos(t355);
t472 = t362 * t346;
t471 = t362 * t365;
t470 = t362 * t368;
t457 = t374 * t365;
t456 = t374 * t368;
t455 = t378 * t365;
t454 = t378 * t368;
t148 = rSges(6,1) * t260 + rSges(6,2) * t259 + rSges(6,3) * t288;
t453 = -(t335 - t363) * t370 - (t374 * t444 - t378 * t421) * t367 - t148;
t158 = rSges(5,1) * t213 + rSges(5,2) * t214;
t159 = rSges(5,1) * t215 + rSges(5,2) * t216;
t102 = t158 * t469 + t159 * t464;
t144 = t370 * t159;
t391 = t370 * t455 + t456;
t162 = t391 * t496 + t402;
t160 = t370 * t162;
t452 = t144 + t160;
t310 = t370 * t454 - t457;
t161 = -t310 * t496 + t403;
t451 = t161 * t469 + t162 * t464;
t450 = -t158 - t161;
t345 = cos(t354);
t331 = t345 * t506;
t280 = t331 + t472 / 0.2e1;
t248 = t280 * t368 - t356 * t471;
t344 = sin(t355);
t332 = t344 * t506;
t282 = t332 - t473 / 0.2e1;
t249 = t282 * t368 - t357 * t471;
t174 = rSges(4,1) * t248 + rSges(4,2) * t249;
t250 = -t280 * t365 - t356 * t470;
t251 = -t282 * t365 - t357 * t470;
t175 = rSges(4,1) * t250 + rSges(4,2) * t251;
t113 = t174 * t469 + t175 * t464;
t441 = t330 - t352;
t327 = -t367 * t379 + t370 * t505;
t443 = t267 - t327;
t182 = t441 * t365 + t368 * t443;
t183 = -t365 * t443 + t368 * t441;
t449 = t182 * t469 + t183 * t464;
t176 = t370 * t183;
t189 = rSges(5,1) * t266 + rSges(5,2) * t265 + rSges(5,3) * t469;
t181 = t370 * t189;
t448 = t176 + t181;
t188 = rSges(5,1) * t264 + rSges(5,2) * t263 - rSges(5,3) * t464;
t125 = t188 * t469 + t189 * t464;
t447 = -t182 - t188;
t337 = t346 / 0.2e1;
t323 = t337 + t345 / 0.2e1;
t274 = t323 * t368 - t356 * t365;
t336 = t343 / 0.2e1;
t322 = t336 - t344 / 0.2e1;
t275 = t322 * t368 + t357 * t365;
t197 = rSges(4,1) * t275 + rSges(4,2) * t274 - rSges(4,3) * t464;
t276 = -t323 * t365 - t356 * t368;
t277 = -t322 * t365 + t357 * t368;
t198 = rSges(4,1) * t277 + rSges(4,2) * t276 + rSges(4,3) * t469;
t126 = t197 * t469 + t198 * t464;
t420 = pkin(7) * t367 + t327;
t246 = t502 * t365 + t420 * t368;
t247 = -t420 * t365 + t502 * t368;
t446 = t246 * t469 + t247 * t464;
t234 = (t363 - t379) * t370 + (t500 * t374 + t435) * t367;
t241 = rSges(5,1) * t309 + rSges(5,2) * t306 + rSges(5,3) * t370;
t445 = -t234 - t241;
t297 = t310 * qJD(2);
t283 = pkin(2) * t297;
t299 = t391 * qJD(2);
t284 = pkin(2) * t299;
t442 = t283 * t469 - t284 * t464;
t440 = qJD(2) * t367;
t432 = t160 + t498;
t431 = -t161 - t523;
t430 = t176 + t494;
t429 = -t182 - t520;
t424 = -t234 + t453;
t18 = t115 * t252 + t117 * t129 + t119 * t130 + t200 * t86 + t254 * t84 - t512 * t88;
t19 = t116 * t252 + t118 * t129 + t120 * t130 + t200 * t87 + t254 * t85 - t512 * t89;
t38 = t105 * t254 + t106 * t200 - t107 * t512 + t129 * t146 + t130 * t147 + t145 * t252;
t10 = t370 * t38 + (t18 * t365 - t19 * t368) * t367;
t24 = t115 * t412 + t117 * t179 + t119 * t180 + t259 * t86 + t260 * t88 + t288 * t84;
t25 = t116 * t412 + t118 * t179 + t120 * t180 + t259 * t87 + t260 * t89 + t288 * t85;
t43 = t105 * t288 + t106 * t259 + t107 * t260 + t145 * t412 + t146 * t179 + t147 * t180;
t14 = t370 * t43 + (t24 * t365 - t25 * t368) * t367;
t422 = (t14 + ((t155 * t306 + t157 * t309 + t185 * t268 + t187 * t270) * t365 - (t154 * t306 + t156 * t309 + t184 * t268 + t186 * t270) * t368) * t367 + (t204 * t306 + t205 * t309 + t232 * t268 + t233 * t270 + (-t152 * t368 + t153 * t365) * t367 + t518) * t370) * t370 + (t10 - (t154 * t265 + t156 * t266 + t184 * t216 + t186 * t215) * t464 + (t204 * t265 + t205 * t266 + t215 * t233 + t216 * t232) * t370 + (t153 * t469 + t155 * t265 + t157 * t266 + t185 * t216 + t187 * t215 + t518 - t519) * t469) * t469;
t419 = t495 * t367;
t17 = t522 * t464 + t469 * t523;
t418 = t453 * t367;
t206 = rSges(5,1) * t270 + rSges(5,2) * t268;
t207 = (t500 * t438 + t380) * t367;
t417 = (-t206 - t207) * t367;
t416 = t445 * t367;
t321 = t336 + t344 / 0.2e1;
t324 = t337 - t345 / 0.2e1;
t258 = rSges(4,1) * t324 + rSges(4,2) * t321 + rSges(4,3) * t370;
t301 = t367 * t505 + (-pkin(7) + t379) * t370;
t415 = (-t258 - t301) * t367;
t45 = t464 * t521 + t469 * t520;
t413 = pkin(2) * t367 ^ 2 * t438;
t70 = t451 + t102;
t78 = t125 + t449;
t407 = (-t207 + t495) * t367;
t406 = t424 * t367;
t405 = (-t301 + t445) * t367;
t168 = Icges(4,5) * t248 + Icges(4,6) * t249;
t169 = Icges(4,5) * t250 + Icges(4,6) * t251;
t170 = Icges(4,4) * t248 + Icges(4,2) * t249;
t171 = Icges(4,4) * t250 + Icges(4,2) * t251;
t172 = Icges(4,1) * t248 + Icges(4,4) * t249;
t173 = Icges(4,1) * t250 + Icges(4,4) * t251;
t193 = Icges(4,4) * t275 + Icges(4,2) * t274 - Icges(4,6) * t464;
t194 = Icges(4,4) * t277 + Icges(4,2) * t276 + Icges(4,6) * t469;
t195 = Icges(4,1) * t275 + Icges(4,4) * t274 - Icges(4,5) * t464;
t196 = Icges(4,1) * t277 + Icges(4,4) * t276 + Icges(4,5) * t469;
t279 = t331 - t472 / 0.2e1;
t281 = t332 + t473 / 0.2e1;
t208 = Icges(4,5) * t281 + Icges(4,6) * t279;
t209 = Icges(4,4) * t281 + Icges(4,2) * t279;
t210 = Icges(4,1) * t281 + Icges(4,4) * t279;
t256 = Icges(4,4) * t324 + Icges(4,2) * t321 + Icges(4,6) * t370;
t257 = Icges(4,1) * t324 + Icges(4,4) * t321 + Icges(4,5) * t370;
t404 = ((t169 * t469 + t171 * t276 + t173 * t277 + t194 * t251 + t196 * t250) * t469 - (t168 * t469 + t170 * t276 + t172 * t277 + t193 * t251 + t195 * t250) * t464 + (t208 * t469 + t209 * t276 + t210 * t277 + t250 * t257 + t251 * t256) * t370) * t469 + t370 * ((t208 * t370 + t209 * t321 + t210 * t324 + t256 * t279 + t257 * t281) * t370 + ((t169 * t370 + t171 * t321 + t173 * t324 + t194 * t279 + t196 * t281) * t365 - (t168 * t370 + t170 * t321 + t172 * t324 + t193 * t279 + t195 * t281) * t368) * t367) + t422;
t16 = t17 + t451;
t395 = (-t301 + t424) * t367;
t44 = t45 + t449;
t311 = t370 * t456 + t455;
t390 = t370 * t457 - t454;
t54 = t115 * t254 + t117 * t200 - t119 * t512;
t55 = t116 * t254 + t118 * t200 - t120 * t512;
t65 = t145 * t254 + t146 * t200 - t147 * t512;
t3 = t18 * t254 + t19 * t255 + t252 * t54 + t253 * t55 + t288 * t38 + t412 * t65;
t56 = t115 * t255 + t117 * t202 - t119 * t513;
t57 = t116 * t255 + t118 * t202 - t120 * t513;
t66 = t145 * t255 + t146 * t202 - t147 * t513;
t4 = t20 * t254 + t21 * t255 + t252 * t56 + t253 * t57 + t288 * t39 + t412 * t66;
t63 = t115 * t288 + t117 * t259 + t119 * t260;
t64 = t116 * t288 + t118 * t259 + t120 * t260;
t71 = t145 * t288 + t146 * t259 + t147 * t260;
t6 = t24 * t254 + t25 * t255 + t252 * t63 + t253 * t64 + t288 * t43 + t412 * t71;
t389 = t3 * t469 / 0.2e1 + t288 * t14 / 0.2e1 - t4 * t464 / 0.2e1 + t252 * (t370 * t65 + (t365 * t54 - t368 * t55) * t367) / 0.2e1 + t253 * (t370 * t66 + (t365 * t56 - t368 * t57) * t367) / 0.2e1 + (t370 * t71 + (t365 * t63 - t368 * t64) * t367) * t412 / 0.2e1 + t370 * t6 / 0.2e1 + t254 * t10 / 0.2e1 + t255 * t11 / 0.2e1;
t211 = rSges(4,1) * t281 + rSges(4,2) * t279;
t388 = -t211 * t367 - t413;
t384 = t499 * t464 + t422;
t383 = -t413 + t417;
t382 = t407 - t413;
t37 = (-t169 * t464 + t171 * t274 + t173 * t275 + t194 * t249 + t196 * t248) * t469 - (-t168 * t464 + t170 * t274 + t172 * t275 + t193 * t249 + t195 * t248) * t464 + (-t208 * t464 + t209 * t274 + t210 * t275 + t248 * t257 + t249 * t256) * t370;
t381 = (-t37 + t499) * t464 + t404;
t305 = (rSges(3,1) * t378 - rSges(3,2) * t374) * t440;
t304 = (Icges(3,1) * t378 - t492) * t440;
t303 = (-Icges(3,2) * t374 + t491) * t440;
t302 = (Icges(3,5) * t378 - Icges(3,6) * t374) * t440;
t300 = t390 * qJD(2);
t298 = t311 * qJD(2);
t287 = rSges(3,3) * t370 + (rSges(3,1) * t374 + rSges(3,2) * t378) * t367;
t286 = Icges(3,5) * t370 + (Icges(3,1) * t374 + t491) * t367;
t285 = Icges(3,6) * t370 + (Icges(3,2) * t378 + t492) * t367;
t278 = t370 * t284;
t243 = -rSges(3,1) * t299 + rSges(3,2) * t300;
t242 = rSges(3,1) * t297 - rSges(3,2) * t298;
t240 = -Icges(3,1) * t299 + Icges(3,4) * t300;
t239 = Icges(3,1) * t297 - Icges(3,4) * t298;
t238 = -Icges(3,4) * t299 + Icges(3,2) * t300;
t237 = Icges(3,4) * t297 - Icges(3,2) * t298;
t236 = -Icges(3,5) * t299 + Icges(3,6) * t300;
t235 = Icges(3,5) * t297 - Icges(3,6) * t298;
t227 = t370 * t247;
t226 = -rSges(3,1) * t390 - rSges(3,2) * t391 + rSges(3,3) * t469;
t225 = rSges(3,1) * t311 + rSges(3,2) * t310 - rSges(3,3) * t464;
t224 = -Icges(3,1) * t390 - Icges(3,4) * t391 + Icges(3,5) * t469;
t223 = Icges(3,1) * t311 + Icges(3,4) * t310 - Icges(3,5) * t464;
t222 = -Icges(3,4) * t390 - Icges(3,2) * t391 + Icges(3,6) * t469;
t221 = Icges(3,4) * t311 + Icges(3,2) * t310 - Icges(3,6) * t464;
t192 = t370 * t198;
t165 = t370 * t175;
t149 = (t242 * t365 + t243 * t368) * t367;
t137 = -t197 * t370 - t258 * t464;
t136 = -t258 * t469 + t192;
t134 = -t188 * t370 - t241 * t464;
t133 = -t241 * t469 + t181;
t128 = -t174 * t370 - t211 * t464;
t127 = -t211 * t469 + t165;
t124 = -t158 * t370 - t206 * t464;
t123 = -t206 * t469 + t144;
t112 = (-t174 - t283) * t370 + t388 * t368;
t111 = t365 * t388 + t165 - t278;
t104 = (-t197 - t246) * t370 + t368 * t415;
t103 = t365 * t415 + t192 + t227;
t101 = t442 + t113;
t100 = t446 + t126;
t99 = t368 * t416 + t370 * t447;
t98 = t365 * t416 + t448;
t80 = (-t246 + t447) * t370 + t368 * t405;
t79 = t365 * t405 + t227 + t448;
t77 = t368 * t417 + t370 * t450;
t76 = t365 * t417 + t452;
t75 = -t122 * t288 + t148 * t255;
t74 = t121 * t288 - t148 * t254;
t73 = (-t283 + t450) * t370 + t383 * t368;
t72 = t365 * t383 - t278 + t452;
t69 = t78 + t446;
t68 = -t121 * t255 + t122 * t254;
t67 = t70 + t442;
t53 = t368 * t418 - t370 * t520;
t52 = t365 * t418 + t494;
t51 = t368 * t406 + t370 * t429;
t50 = t365 * t406 + t430;
t49 = t110 * t255 - t122 * t412 + t148 * t253 - t288 * t91;
t48 = -t110 * t254 + t121 * t412 - t148 * t252 + t288 * t90;
t47 = (-t246 + t429) * t370 + t368 * t395;
t46 = t365 * t395 + t227 + t430;
t42 = -t121 * t253 + t122 * t252 + t254 * t91 - t255 * t90;
t41 = t44 + t446;
t34 = t368 * t419 - t370 * t523;
t33 = t365 * t419 + t498;
t32 = t368 * t407 + t370 * t431;
t31 = t365 * t407 + t432;
t29 = (-t283 + t431) * t370 + t382 * t368;
t28 = t365 * t382 - t278 + t432;
t15 = t16 + t442;
t1 = [0; m(3) * t149 + m(4) * t101 + m(5) * t67 + m(6) * t15; -t11 * t464 - t30 * t464 - t37 * t464 - ((-t222 * t298 + t224 * t297 - t236 * t464 + t238 * t310 + t240 * t311) * t469 - (-t221 * t298 + t223 * t297 - t235 * t464 + t237 * t310 + t239 * t311) * t464 + (-t285 * t298 + t286 * t297 - t302 * t464 + t303 * t310 + t304 * t311) * t370) * t464 + ((t222 * t300 - t224 * t299 + t236 * t469 - t238 * t391 - t240 * t390) * t469 - (t221 * t300 - t223 * t299 + t235 * t469 - t237 * t391 - t239 * t390) * t464 + (t285 * t300 - t286 * t299 + t302 * t469 - t303 * t391 - t304 * t390) * t370) * t469 + (t15 * t41 + t28 * t46 + t29 * t47) * t508 + (t67 * t69 + t72 * t79 + t73 * t80) * t509 + (t100 * t101 + t103 * t111 + t104 * t112) * t510 + 0.2e1 * m(3) * ((-t225 * t370 - t287 * t464) * (-t242 * t370 - t305 * t464) + (t226 * t370 - t287 * t469) * (t243 * t370 - t305 * t469) + (t225 * t365 + t226 * t368) * t367 * t149) + t370 * (t370 ^ 2 * t302 + (((t238 * t378 + t240 * t374) * t365 - (t237 * t378 + t239 * t374) * t368 + ((-t222 * t374 + t224 * t378) * t365 - (-t221 * t374 + t223 * t378) * t368) * qJD(2)) * t367 + (-t235 * t368 + t236 * t365 + t303 * t378 + t304 * t374 + (-t285 * t374 + t286 * t378) * qJD(2)) * t370) * t367) + t404; m(4) * t113 + m(5) * t70 + m(6) * t16; m(6) * (t15 * t44 + t16 * t41 + t28 * t50 + t29 * t51 + t31 * t46 + t32 * t47) + m(5) * (t67 * t78 + t69 * t70 + t72 * t98 + t73 * t99 + t76 * t79 + t77 * t80) + m(4) * (t100 * t113 + t101 * t126 + t103 * t127 + t104 * t128 + t111 * t136 + t112 * t137) + t381; (t16 * t44 + t31 * t50 + t32 * t51) * t508 + (t70 * t78 + t76 * t98 + t77 * t99) * t509 + (t113 * t126 + t127 * t136 + t128 * t137) * t510 + t381; m(5) * t102 + m(6) * t17; m(6) * (t15 * t45 + t17 * t41 + t28 * t52 + t29 * t53 + t33 * t46 + t34 * t47) + m(5) * (t102 * t69 + t123 * t79 + t124 * t80 + t125 * t67 + t133 * t72 + t134 * t73) + t384; m(6) * (t16 * t45 + t17 * t44 + t31 * t52 + t32 * t53 + t33 * t50 + t34 * t51) + m(5) * (t102 * t78 + t123 * t98 + t124 * t99 + t125 * t70 + t133 * t76 + t134 * t77) + t384; (t102 * t125 + t123 * t133 + t124 * t134) * t509 + (t17 * t45 + t33 * t52 + t34 * t53) * t508 + t384; m(6) * t42; m(6) * (t15 * t68 + t28 * t74 + t29 * t75 + t41 * t42 + t46 * t48 + t47 * t49) + t389; m(6) * (t16 * t68 + t31 * t74 + t32 * t75 + t42 * t44 + t48 * t50 + t49 * t51) + t389; m(6) * (t17 * t68 + t33 * t74 + t34 * t75 + t42 * t45 + t48 * t52 + t49 * t53) + t389; (t42 * t68 + t48 * t74 + t49 * t75) * t508 + t252 * (t254 * t54 + t255 * t55 + t288 * t65) + t254 * t3 + t253 * (t254 * t56 + t255 * t57 + t288 * t66) + t255 * t4 + (t254 * t63 + t255 * t64 + t288 * t71) * t412 + t288 * t6;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
