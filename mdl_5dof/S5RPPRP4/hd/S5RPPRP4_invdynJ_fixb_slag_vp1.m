% Calculate vector of inverse dynamics joint torques for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:03
% EndTime: 2019-12-31 17:52:19
% DurationCPUTime: 13.41s
% Computational Cost: add. (6874->514), mult. (15426->595), div. (0->0), fcn. (16109->6), ass. (0->268)
t251 = sin(qJ(4));
t414 = sin(pkin(7));
t415 = cos(pkin(7));
t436 = sin(qJ(1));
t437 = cos(qJ(1));
t191 = -t414 * t436 - t415 * t437;
t192 = t437 * t414 - t436 * t415;
t252 = cos(qJ(4));
t410 = Icges(6,4) * t252;
t308 = -Icges(6,2) * t251 + t410;
t104 = Icges(6,6) * t191 + t192 * t308;
t412 = Icges(5,4) * t252;
t310 = -Icges(5,2) * t251 + t412;
t107 = Icges(5,6) * t191 + t192 * t310;
t488 = t104 + t107;
t536 = t251 * t488;
t513 = Icges(5,3) + Icges(6,3);
t304 = Icges(6,5) * t252 - Icges(6,6) * t251;
t306 = Icges(5,5) * t252 - Icges(5,6) * t251;
t535 = t304 + t306;
t311 = Icges(6,1) * t251 + t410;
t313 = Icges(5,1) * t251 + t412;
t529 = t311 + t313;
t411 = Icges(6,4) * t251;
t307 = Icges(6,2) * t252 + t411;
t413 = Icges(5,4) * t251;
t309 = Icges(5,2) * t252 + t413;
t530 = -t307 - t309;
t456 = t530 * t251 + t529 * t252;
t305 = Icges(5,5) * t251 + Icges(5,6) * t252;
t490 = t305 * t191;
t303 = Icges(6,5) * t251 + Icges(6,6) * t252;
t491 = t303 * t191;
t472 = t192 * t456 + t490 + t491;
t534 = t472 * qJD(1);
t503 = -t191 * t535 + t192 * t513;
t533 = t503 * t192;
t504 = t513 * t191 + t535 * t192;
t532 = t504 * t191;
t312 = Icges(6,1) * t252 - t411;
t112 = Icges(6,5) * t192 - t191 * t312;
t314 = Icges(5,1) * t252 - t413;
t115 = Icges(5,5) * t192 - t191 * t314;
t522 = t112 + t115;
t113 = Icges(5,5) * t191 + t192 * t314;
t531 = t113 * t252 - t536;
t110 = Icges(6,5) * t191 + t192 * t312;
t465 = -t110 - t113;
t474 = t251 * t465 - t252 * t488;
t526 = Icges(5,1) + Icges(6,1);
t525 = Icges(5,4) + Icges(6,4);
t515 = Icges(5,5) + Icges(6,5);
t524 = Icges(5,2) + Icges(6,2);
t514 = Icges(5,6) + Icges(6,6);
t106 = Icges(6,6) * t192 - t191 * t308;
t109 = Icges(5,6) * t192 - t191 * t310;
t523 = t106 + t109;
t496 = t110 * t252 + t531;
t396 = t191 * t252;
t448 = -t522 * t396 + t533;
t394 = t192 * t252;
t521 = t503 * t191 + t522 * t394;
t478 = t496 * t192 + t532;
t395 = t192 * t251;
t477 = t523 * t395 - t521;
t476 = t110 * t396 + t531 * t191 - t504 * t192;
t397 = t191 * t251;
t475 = t523 * t397 + t448;
t177 = t192 * qJD(1);
t176 = t191 * qJD(1);
t360 = t251 * qJD(4);
t348 = t192 * t360;
t282 = -t176 * t252 + t348;
t364 = qJD(4) * t252;
t283 = t176 * t251 + t192 * t364;
t520 = t514 * t177 + t525 * t282 + t524 * t283;
t351 = t191 * t360;
t400 = t177 * t252;
t280 = t351 + t400;
t363 = t191 * qJD(4);
t401 = t177 * t251;
t281 = t252 * t363 - t401;
t519 = t514 * t176 + t525 * t280 + t524 * t281;
t518 = t515 * t177 + t526 * t282 + t525 * t283;
t517 = -t515 * t176 - t526 * t280 - t525 * t281;
t512 = (t308 + t310) * qJD(4);
t511 = (t312 + t314) * qJD(4);
t510 = -t529 * t251 + t530 * t252;
t473 = t522 * t251 + t523 * t252;
t124 = t303 * t192;
t126 = t305 * t192;
t471 = t456 * t191 - t124 - t126;
t505 = rSges(6,3) + qJ(5) + pkin(6);
t509 = t523 * t251;
t507 = t513 * t177 + t515 * t282 + t514 * t283;
t506 = t513 * t176 + t515 * t280 + t514 * t281;
t502 = t535 * qJD(4);
t501 = t510 * qJD(4) - t512 * t251 + t511 * t252;
t500 = t473 * qJD(4) + t519 * t251 + t517 * t252;
t499 = -t474 * qJD(4) - t520 * t251 + t518 * t252;
t497 = -t522 * t252 + t509;
t495 = -t303 - t305;
t494 = -t476 * t191 + t475 * t192;
t493 = -t478 * t191 + t477 * t192;
t438 = rSges(6,1) + pkin(4);
t492 = t471 * qJD(1);
t324 = rSges(5,1) * t251 + rSges(5,2) * t252;
t365 = qJD(4) * t192;
t489 = t324 * t365;
t121 = -rSges(5,1) * t396 + rSges(5,2) * t397 + t192 * rSges(5,3);
t150 = -t191 * pkin(3) + pkin(6) * t192;
t248 = t437 * pkin(2);
t366 = t437 * pkin(1) + t436 * qJ(2);
t458 = t248 + t366;
t288 = t150 + t458;
t487 = t121 + t288;
t148 = -t191 * rSges(4,1) - t192 * rSges(4,2);
t486 = t148 + t458;
t462 = t505 * t191;
t69 = rSges(5,1) * t280 + rSges(5,2) * t281 + t176 * rSges(5,3);
t243 = t437 * qJ(2);
t356 = t436 * pkin(1);
t220 = t356 - t243;
t240 = qJD(2) * t436;
t347 = qJD(1) * t437;
t369 = qJ(2) * t347 + t240;
t450 = qJD(1) * t220 - t240 + t369;
t141 = qJD(4) * t176 + qJDD(4) * t192;
t244 = t251 * rSges(6,2);
t459 = -rSges(6,1) * t252 + t244;
t200 = t459 * qJD(4);
t254 = qJD(1) ^ 2;
t346 = qJD(1) * t436;
t274 = qJDD(1) * t366 - qJDD(2) * t437 + (-pkin(1) * t346 + t240 + t369) * qJD(1);
t355 = t436 * pkin(2);
t258 = qJDD(1) * t248 - t254 * t355 + t274;
t340 = t177 * pkin(3) + t176 * pkin(6);
t255 = qJD(1) * t340 + qJDD(1) * t150 + t258;
t421 = t252 * rSges(6,2);
t322 = rSges(6,1) * t251 + t421;
t390 = t252 * qJD(4) ^ 2;
t435 = pkin(4) * t252;
t236 = pkin(3) + t435;
t463 = -rSges(6,1) * t396 + rSges(6,2) * t397 - t191 * t236 + t505 * t192;
t416 = -t150 + t463;
t179 = qJD(5) * t192;
t449 = rSges(6,1) * t400 + t505 * t176 + t177 * t236;
t430 = rSges(6,2) * t281 + t351 * t438 + t179 - t340 + t449;
t11 = -t200 * t365 + qJD(5) * t177 - qJDD(5) * t191 + t141 * t322 + t416 * qJDD(1) + t430 * qJD(1) + (t141 * t251 + t192 * t390) * pkin(4) + t255;
t485 = t11 - g(2);
t484 = t493 * qJD(4) + t534;
t483 = t494 * qJD(4) + t492;
t482 = -t496 * qJD(4) + t518 * t251 + t520 * t252;
t481 = t497 * qJD(4) + t517 * t251 - t519 * t252;
t480 = t456 * t176 + t495 * t177 + t502 * t191 + t501 * t192;
t479 = t495 * t176 - t456 * t177 + t501 * t191 - t502 * t192;
t271 = t288 + t416;
t341 = pkin(4) * t251 + t322;
t315 = qJD(4) * t341;
t178 = t191 * qJD(5);
t241 = qJD(2) * t437;
t375 = t178 + t241;
t27 = t271 * qJD(1) + t192 * t315 - t375;
t469 = t27 * t341;
t467 = t324 * t363;
t464 = rSges(6,1) * t394 - rSges(6,2) * t395 + t192 * t236 + t462;
t461 = t192 * rSges(4,1) - t191 * rSges(4,2);
t188 = t191 * pkin(6);
t339 = t192 * pkin(3) + t188;
t276 = -t356 - t355;
t457 = t276 + t355;
t455 = (t503 * t176 - t497 * t177 + t506 * t192) * t192 + (t499 * t191 + t496 * t177 + t504 * t176 + (t500 - t507) * t192) * t191;
t454 = (t497 * t176 + t503 * t177 + t500 * t192) * t192 + (t507 * t191 + t504 * t177 - t496 * t176 + (t499 - t506) * t192) * t191;
t453 = -t505 * t177 + t178;
t367 = t437 * rSges(3,1) + t436 * rSges(3,3);
t452 = t366 + t367;
t451 = -qJD(1) * t339 + t450;
t284 = t251 * t438 + t421;
t447 = -qJD(1) * t458 + t241;
t383 = t309 * t192 - t113;
t387 = -t313 * t192 - t107;
t446 = t251 * t383 + t252 * t387;
t385 = t307 * t192 - t110;
t389 = -t311 * t192 - t104;
t445 = t251 * t385 + t252 * t389;
t444 = t141 / 0.2e1;
t142 = qJD(4) * t177 - qJDD(4) * t191;
t443 = t142 / 0.2e1;
t432 = pkin(3) - t236;
t172 = t177 * pkin(6);
t431 = rSges(6,1) * t282 + rSges(6,2) * t283 + pkin(4) * t348 + t176 * t432 - t172 - t453;
t334 = -t220 - t355;
t289 = t334 + t339;
t417 = t188 - t462 + (t459 + t432) * t192;
t270 = t289 - t417;
t26 = qJD(1) * t270 + t191 * t315 + t179 + t240;
t424 = t176 * t26;
t423 = t177 * rSges(5,3);
t181 = t191 * rSges(5,3);
t418 = -t339 + t464;
t388 = -t311 * t191 + t106;
t386 = -t313 * t191 + t109;
t384 = t307 * t191 + t112;
t382 = t309 * t191 + t115;
t173 = qJD(1) * t366 - t241;
t381 = t176 * pkin(3) - t172 - t173;
t376 = t177 * rSges(4,1) - t176 * rSges(4,2);
t374 = -t307 + t312;
t373 = -t308 - t311;
t372 = -t309 + t314;
t371 = -t310 - t313;
t370 = qJD(1) * t241 + qJDD(2) * t436;
t362 = t304 * qJD(1);
t361 = t306 * qJD(1);
t359 = m(4) + m(5) + m(6);
t358 = -t437 / 0.2e1;
t357 = t436 / 0.2e1;
t354 = t436 * rSges(3,1);
t345 = -t365 / 0.2e1;
t344 = t365 / 0.2e1;
t343 = -t363 / 0.2e1;
t342 = t363 / 0.2e1;
t337 = pkin(4) * t395 + t322 * t192;
t336 = pkin(4) * t397 + t322 * t191;
t326 = t26 * t341;
t325 = t176 * rSges(4,1) + t177 * rSges(4,2);
t224 = -rSges(5,1) * t252 + rSges(5,2) * t251;
t119 = t192 * t224 - t181;
t278 = -t119 + t289;
t43 = qJD(1) * t278 + t240 + t467;
t44 = qJD(1) * t487 - t241 + t489;
t317 = -t191 * t43 - t192 * t44;
t67 = rSges(5,1) * t282 + rSges(5,2) * t283 + t423;
t316 = -t191 * t69 - t192 * t67;
t294 = t119 * t192 + t121 * t191;
t287 = t461 + t334;
t117 = rSges(5,1) * t394 - rSges(5,2) * t395 + t181;
t286 = pkin(3) - t224;
t285 = t236 - t459;
t279 = -t248 * t254 + t370;
t275 = -t354 - t356;
t227 = rSges(2,1) * t437 - rSges(2,2) * t436;
t222 = rSges(2,1) * t436 + rSges(2,2) * t437;
t269 = t284 * t191;
t268 = t243 + t276;
t267 = t251 * t384 + t252 * t388;
t266 = t251 * t382 + t252 * t386;
t265 = (t251 * t373 + t252 * t374) * qJD(1);
t264 = (t251 * t371 + t252 * t372) * qJD(1);
t246 = t437 * rSges(3,3);
t238 = rSges(3,3) * t347;
t221 = t354 - t246;
t201 = t224 * qJD(4);
t139 = t324 * t191;
t137 = t324 * t192;
t96 = qJD(1) * t287 + t240;
t78 = qJDD(1) * t367 + qJD(1) * (-rSges(3,1) * t346 + t238) + t274;
t77 = -qJD(1) * t173 - t254 * t367 + t370 + (-t220 - t221) * qJDD(1);
t46 = qJD(1) * t376 + qJDD(1) * t148 + t258;
t45 = (-t173 + t325) * qJD(1) + t287 * qJDD(1) + t279;
t42 = qJD(4) * t294 - qJD(3);
t23 = -qJD(3) + (t191 * t416 + t192 * t417) * qJD(4);
t18 = qJD(1) * t69 + qJDD(1) * t121 + t141 * t324 - t201 * t365 + t255;
t17 = -t201 * t363 - t142 * t324 + (-t67 + t381) * qJD(1) + t278 * qJDD(1) + t279;
t16 = qJD(4) * t316 - t119 * t141 + t121 * t142 + qJDD(3);
t10 = -t200 * t363 + qJD(5) * t176 + qJDD(5) * t192 - t142 * t322 + (-t142 * t251 + t191 * t390) * pkin(4) + (t381 - t431) * qJD(1) + t270 * qJDD(1) + t279;
t1 = qJDD(3) + t416 * t142 - t417 * t141 + (-t430 * t191 - t431 * t192) * qJD(4);
t2 = [-m(2) * (-g(1) * t222 + g(2) * t227) + (((-t476 + t477 + t521) * t191 + ((t252 * t465 + t536) * t192 - t532 + t448 + t478) * t192) * qJD(4) + t492) * t342 + (t456 * qJD(4) + t511 * t251 + t512 * t252) * qJD(1) + (t285 * t424 - g(1) * (t268 + t464) + t485 * (t458 + t463) + (t285 * t192 + t268 + t462) * t10 + (-rSges(6,2) * t401 + t449 + t451) * t27 + (t241 - t375 + t453) * t26 + (t27 * t269 - (t469 + t23 * (t417 + t418)) * t191 + (-t26 * t284 + t326) * t192) * qJD(4) + ((-t418 + t457) * t27 + (-t458 + t271) * t26) * qJD(1)) * m(6) + (-t42 * (t117 + t119) * t363 - g(1) * (t117 + t268 + t339) + (t18 - g(2)) * t487 + (t286 * t192 + t181 + t188 + t268) * t17 + (t286 * t176 - t172 - t423 + t447 - t489) * t43 + (-t467 + t43 + (-t117 + t457) * qJD(1) + t340 + t451 + t69) * t44) * m(5) + ((t325 + t447) * t96 + (t376 + t96 + (t457 - t461) * qJD(1) + t450) * (qJD(1) * t486 - t241) + (t46 - g(2)) * t486 + (t45 - g(1)) * (t268 + t461)) * m(4) + ((t78 - g(2)) * t452 + (t77 - g(1)) * (t243 + t246 + t275) + (t238 + (t221 + t275) * qJD(1) + t450) * (qJD(1) * t452 - t241)) * m(3) + (t471 - t473) * t444 + (t472 - t474) * t443 + (t479 + t481) * t344 + (t480 - t482 + t483) * t343 + (m(2) * (t222 ^ 2 + t227 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3) - t510) * qJDD(1) + ((((t504 - t497) * t192 + t476) * t192 + ((t504 - t509) * t191 + t533 - t448 + t475) * t191) * qJD(4) + t484 - t534) * t345; (-m(3) - t359) * (g(1) * t436 - g(2) * t437) + 0.2e1 * (t10 * t357 + t11 * t358) * m(6) + 0.2e1 * (t17 * t357 + t18 * t358) * m(5) + 0.2e1 * (t357 * t45 + t358 * t46) * m(4) + 0.2e1 * (t357 * t77 + t358 * t78) * m(3); m(4) * qJDD(3) + m(5) * t16 + m(6) * t1 + g(3) * t359; t494 * t444 + t493 * t443 + t483 * t176 / 0.2e1 + t484 * t177 / 0.2e1 - (t480 * qJD(1) + t454 * qJD(4) + t472 * qJDD(1) + t477 * t141 + t478 * t142) * t191 / 0.2e1 + (t479 * qJD(1) + t455 * qJD(4) + t471 * qJDD(1) + t475 * t141 + t476 * t142) * t192 / 0.2e1 - ((((-t382 - t384) * t192 + (t383 + t385) * t191) * t252 + ((t386 + t388) * t192 + (-t387 - t389) * t191) * t251) * qJD(4) + ((-t371 - t373) * t252 + (t372 + t374) * t251) * qJD(1)) * qJD(1) / 0.2e1 + (-t473 * t176 - t474 * t177 + t482 * t191 + t481 * t192) * qJD(1) / 0.2e1 + (t474 * t191 - t473 * t192) * qJDD(1) / 0.2e1 + ((t365 * t490 - t361) * t192 + (t264 + (-t446 * t191 + (-t126 + t266) * t192) * qJD(4)) * t191 + (t365 * t491 - t362) * t192 + (t265 + (-t445 * t191 + (-t124 + t267) * t192) * qJD(4)) * t191) * t345 + (t176 * t475 + t177 * t476 + t455) * t344 + (t176 * t477 + t177 * t478 + t454) * t343 + ((t126 * t363 + t361) * t191 + (t264 + (t266 * t192 + (-t490 - t446) * t191) * qJD(4)) * t192 + (t124 * t363 + t362) * t191 + (t265 + (t267 * t192 + (-t491 - t445) * t191) * qJD(4)) * t192) * t342 + (-(-t137 * t43 + t139 * t44) * qJD(1) - (t42 * (t137 * t192 + t139 * t191) + t317 * t224) * qJD(4) - t16 * t294 + t42 * (t119 * t176 - t121 * t177 - t316) + t317 * t201 - (-t17 * t191 - t176 * t44 + t177 * t43 - t18 * t192) * t324 - g(1) * t139 - g(2) * t137 - g(3) * t224) * m(5) + (-g(3) * (-t252 * t438 + t244) - g(1) * t269 - (-t23 * t417 - t469) * t176 + (-t23 * t416 - t326) * t177 - (-t26 * t337 + t27 * t336) * qJD(1) + (t10 * t341 + t26 * (pkin(4) * t364 - t200) - t1 * t416 + t23 * t430 - (t26 * (-t459 + t435) + t336 * t23) * qJD(4)) * t191 + (-g(2) * t284 - t1 * t417 + t11 * t341 - t27 * t200 + t23 * t431 - (t23 * t337 - t27 * t459) * qJD(4)) * t192) * m(6); (t177 * t27 + t424 + (-qJD(1) * t27 - g(1) + t10) * t192 + (-t26 * qJD(1) - t485) * t191) * m(6);];
tau = t2;
