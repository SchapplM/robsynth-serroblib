% Calculate vector of inverse dynamics joint torques for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:13
% DurationCPUTime: 11.69s
% Computational Cost: add. (3103->392), mult. (7523->499), div. (0->0), fcn. (5785->4), ass. (0->211)
t484 = Icges(4,4) + Icges(5,4);
t481 = Icges(4,1) + Icges(5,1);
t479 = Icges(4,2) + Icges(5,2);
t480 = Icges(4,5) + Icges(5,5);
t478 = Icges(4,6) + Icges(5,6);
t225 = cos(qJ(3));
t483 = t484 * t225;
t223 = sin(qJ(3));
t482 = t484 * t223;
t460 = t225 * t479 + t482;
t459 = t223 * t481 + t483;
t477 = Icges(4,3) + Icges(5,3);
t224 = sin(qJ(1));
t339 = t224 * t225;
t476 = t484 * t339;
t474 = t223 * t480 + t225 * t478;
t226 = cos(qJ(1));
t475 = t480 * t226;
t470 = t224 * t460 + t226 * t478;
t469 = -t224 * t478 + t226 * t460;
t342 = t223 * t224;
t468 = t342 * t481 + t475 + t476;
t467 = -t224 * t480 + t226 * t459;
t462 = -t223 * t479 + t483;
t473 = t225 * t481 - t482;
t472 = t224 * t474 + t226 * t477;
t403 = -t224 * t477 + t226 * t474;
t438 = t223 * t467 + t225 * t469;
t471 = t438 * t226;
t466 = t462 * t226;
t465 = t473 * t224;
t464 = t473 * t226;
t463 = -t223 * t478 + t225 * t480;
t450 = t223 * t468 + t225 * t470;
t404 = t472 * t224;
t451 = t223 * t473 + t225 * t462;
t412 = t226 * t472 + t339 * t470 + t342 * t468;
t411 = -t226 * t403 - t339 * t469 - t342 * t467;
t410 = -t226 * t450 + t404;
t409 = -t224 * t403 + t471;
t407 = t223 * t469 - t225 * t467;
t458 = qJD(1) * t470 - qJD(3) * t466;
t309 = qJD(3) * t224;
t457 = -qJD(1) * t469 - t309 * t462;
t456 = -t464 * qJD(3) + (t224 * t459 + t475) * qJD(1);
t455 = qJD(1) * t467 + qJD(3) * t465;
t430 = t463 * t224;
t454 = t460 * qJD(3);
t453 = t459 * qJD(3);
t452 = t223 * t462 - t225 * t473;
t408 = t223 * t470 - t225 * t468;
t424 = t463 * t226;
t406 = t224 * t451 + t424;
t405 = t226 * t451 - t430;
t448 = (t342 * t479 - t468 - t476) * t226 + (t466 + t467) * t224;
t447 = (t465 - t470) * t226 + (-t464 + t469) * t224;
t446 = t459 + t462;
t445 = t460 - t473;
t444 = t403 * qJD(1);
t443 = t472 * qJD(1);
t442 = t226 ^ 2;
t441 = t451 * qJD(1) - qJD(3) * t474;
t440 = qJD(1) * t450 + qJD(3) * t430 + t444;
t439 = -qJD(1) * t438 - qJD(3) * t424 + t443;
t437 = t224 * t409 + t226 * t410;
t436 = t224 * t411 + t412 * t226;
t435 = qJD(3) * t408 - t223 * t455 + t225 * t457 + t443;
t434 = qJD(1) * t463 + qJD(3) * t452 + t223 * t453 + t225 * t454;
t433 = qJD(3) * t407 + t223 * t456 + t225 * t458 + t444;
t432 = (t223 * t446 + t225 * t445) * qJD(1);
t431 = t406 * qJD(1);
t362 = rSges(5,2) * t225;
t270 = rSges(5,1) * t223 + t362;
t429 = -pkin(3) * t223 - t270;
t222 = -qJ(4) - pkin(5);
t311 = qJD(1) * t226;
t299 = t225 * t311;
t300 = t223 * t311;
t312 = qJD(1) * t224;
t295 = t225 * t309;
t320 = pkin(3) * t295 + qJD(4) * t226;
t427 = t295 + t300;
t428 = rSges(5,1) * t427 + rSges(5,2) * t299 + pkin(3) * t300 + t222 * t312 + t320;
t426 = t474 * qJD(1);
t425 = t405 * qJD(1);
t422 = -t223 * t447 + t225 * t448;
t310 = qJD(3) * t223;
t364 = (-rSges(5,2) * t310 - rSges(5,3) * qJD(1)) * t224 + pkin(5) * t312 + t428;
t219 = t226 * rSges(5,3);
t308 = qJD(3) * t226;
t294 = t225 * t308;
t307 = qJD(4) * t224;
t282 = -pkin(3) * t294 + t307;
t169 = rSges(5,1) * t225 - rSges(5,2) * t223;
t298 = t169 * t308;
t196 = pkin(3) * t342;
t366 = pkin(5) + t222;
t396 = -t226 * t366 + t196;
t365 = -t298 + t282 + (t224 * t270 + t219 + t396) * qJD(1);
t421 = -t224 * t364 + t226 * t365;
t171 = pkin(1) * t226 + qJ(2) * t224;
t207 = qJD(2) * t226;
t117 = qJD(1) * t171 - t207;
t144 = t270 * qJD(3);
t305 = qJD(1) * qJD(3);
t148 = qJDD(3) * t224 + t226 * t305;
t209 = t226 * qJ(2);
t166 = pkin(1) * t224 - t209;
t216 = t224 * rSges(5,3);
t341 = t223 * t226;
t326 = pkin(3) * t341 + t224 * t366 + t226 * t270 - t216;
t372 = pkin(5) * t224;
t278 = t326 - t372;
t246 = -t166 + t278;
t370 = pkin(5) * qJD(1) ^ 2;
t283 = qJDD(4) - t370;
t306 = qJD(1) * qJD(2);
t318 = qJDD(2) * t224 + t226 * t306;
t340 = t223 * qJD(3) ^ 2;
t9 = -t144 * t309 + t148 * t169 + t283 * t226 + (t148 * t225 - t224 * t340) * pkin(3) + t246 * qJDD(1) + (-t117 - t307 - t365) * qJD(1) + t318;
t420 = -g(1) + t9;
t149 = qJDD(3) * t226 - t224 * t305;
t206 = qJD(2) * t224;
t317 = qJ(2) * t311 + t206;
t303 = qJD(1) * (-pkin(1) * t312 + t317) + qJDD(1) * t171 + t224 * t306;
t371 = pkin(5) * t226;
t281 = qJDD(1) * t371 + t303;
t286 = -pkin(3) * t225 - t169;
t108 = rSges(5,1) * t342 + rSges(5,2) * t339 + t219;
t327 = t108 + t396;
t10 = t283 * t224 + t286 * t149 + t327 * qJDD(1) + t364 * qJD(1) + (pkin(3) * t340 + qJD(1) * qJD(4) + qJD(3) * t144 - qJDD(2)) * t226 + t281;
t419 = -g(2) + t10;
t418 = qJD(3) * t436 + t431;
t417 = t437 * qJD(3) - t425;
t416 = -qJD(3) * t450 + t223 * t457 + t225 * t455;
t415 = qJD(3) * t438 - t223 * t458 + t225 * t456;
t414 = t224 * t441 + t226 * t434;
t413 = -t224 * t434 + t226 * t441;
t220 = t226 * rSges(4,3);
t109 = rSges(4,1) * t342 + rSges(4,2) * t339 + t220;
t284 = t171 + t371;
t398 = t109 + t284;
t130 = rSges(5,1) * t339 - rSges(5,2) * t342;
t197 = pkin(3) * t339;
t397 = t130 + t197;
t172 = -rSges(3,2) * t226 + rSges(3,3) * t224;
t394 = t440 * t442 + (t433 * t224 + (-t435 + t439) * t226) * t224;
t393 = t435 * t442 + (t439 * t224 + (-t433 + t440) * t226) * t224;
t271 = rSges(4,1) * t223 + rSges(4,2) * t225;
t145 = t271 * qJD(3);
t170 = rSges(4,1) * t225 - rSges(4,2) * t223;
t111 = -rSges(4,3) * t224 + t226 * t271;
t285 = t111 - t372;
t277 = -t166 + t285;
t133 = t170 * t226;
t71 = -qJD(3) * t133 + (t224 * t271 + t220) * qJD(1);
t15 = -t226 * t370 - t145 * t309 + t148 * t170 + (-t117 - t71) * qJD(1) + t277 * qJDD(1) + t318;
t301 = rSges(4,1) * t427 + rSges(4,2) * t299;
t73 = (-rSges(4,2) * t310 - rSges(4,3) * qJD(1)) * t224 + t301;
t16 = -t224 * t370 + qJD(1) * t73 + qJDD(1) * t109 - t149 * t170 + (qJD(3) * t145 - qJDD(2)) * t226 + t281;
t392 = t15 * t224 - t16 * t226;
t375 = rSges(5,1) + pkin(3);
t389 = t223 * t375 + t362;
t384 = g(1) * t224 - g(2) * t226;
t381 = -pkin(1) - pkin(5);
t379 = t148 / 0.2e1;
t378 = t149 / 0.2e1;
t377 = t224 / 0.2e1;
t374 = rSges(3,2) - pkin(1);
t367 = -pkin(1) + t222;
t361 = rSges(3,3) * t226;
t360 = t10 * t226;
t280 = t169 * t309 + t206 + t320;
t29 = qJD(1) * t246 + t280;
t357 = t226 * t29;
t135 = t170 * t309;
t43 = qJD(1) * t277 + t135 + t206;
t356 = t226 * t43;
t338 = t225 * t226;
t167 = rSges(3,2) * t224 + t361;
t321 = -t166 + t167;
t114 = t171 + t172;
t316 = rSges(3,2) * t312 + rSges(3,3) * t311;
t151 = qJD(1) * t166;
t315 = t206 - t151;
t304 = -rSges(4,3) + t381;
t297 = t223 * t309;
t296 = t223 * t308;
t290 = -t309 / 0.2e1;
t289 = t309 / 0.2e1;
t288 = -t308 / 0.2e1;
t287 = t308 / 0.2e1;
t272 = -t207 + t282;
t173 = rSges(2,1) * t226 - rSges(2,2) * t224;
t168 = rSges(2,1) * t224 + rSges(2,2) * t226;
t44 = qJD(1) * t398 - t170 * t308 - t207;
t265 = t224 * t43 - t226 * t44;
t264 = -t224 * t73 + t226 * t71;
t251 = -t109 * t224 - t111 * t226;
t231 = -t224 * t327 - t226 * t326;
t195 = rSges(5,2) * t341;
t132 = -rSges(5,1) * t338 + t195;
t131 = t170 * t224;
t81 = qJD(1) * t114 - t207;
t80 = qJD(1) * t321 + t206;
t45 = t251 * qJD(3);
t40 = qJD(1) * t316 + qJDD(1) * t172 - qJDD(2) * t226 + t303;
t39 = t321 * qJDD(1) + (-qJD(1) * t172 - t117) * qJD(1) + t318;
t30 = -t298 + (t284 + t327) * qJD(1) + t272;
t28 = t231 * qJD(3);
t1 = [-m(2) * (-g(1) * t168 + g(2) * t173) - t405 * t148 / 0.2e1 + t407 * t379 + (((t404 - t410 + t411) * t224 + (-t471 + (-t450 + t403) * t224 + t409 + t412) * t226) * qJD(3) + t431) * t290 + (-t451 * qJD(3) + t454 * t223 - t453 * t225) * qJD(1) + (-(qJD(1) * t278 - t151 + t280 - t29) * t30 + t29 * (rSges(5,1) * t294 - rSges(5,2) * t296 - t272) + t30 * (-rSges(5,2) * t297 + t317 + t428) + ((-rSges(5,3) + t367) * t357 + (t29 * (-qJ(2) + t429) + t30 * (-rSges(5,3) - pkin(1))) * t224) * qJD(1) + t419 * (-t222 * t226 + t108 + t171 + t196) + t420 * (t224 * t367 + t226 * t389 + t209 - t216)) * m(5) + (-(qJD(1) * t285 + t135 + t315 - t43) * t44 + t43 * (rSges(4,1) * t294 - rSges(4,2) * t296 + t207) + t44 * (-rSges(4,2) * t297 + t301 + t317) + (t304 * t356 + (t43 * (-qJ(2) - t271) + t44 * t304) * t224) * qJD(1) + (-g(2) + t16) * t398 + (-g(1) + t15) * (t224 * t381 + t111 + t209)) * m(4) + (-(qJD(1) * t167 + t315 - t80) * t81 + t80 * t207 + t81 * (t316 + t317) + (t80 * t374 * t226 + (t80 * (-rSges(3,3) - qJ(2)) - t81 * pkin(1)) * t224) * qJD(1) + (-g(2) + t40) * t114 + (-g(1) + t39) * (t224 * t374 + t209 + t361)) * m(3) + (t406 - t408) * t378 + ((((t450 + t403) * t226 - t404 + t411) * t226 + t403 * t224 ^ 2) * qJD(3) + t417 + t425) * t288 + (t413 + t416) * t287 + (m(2) * (t168 ^ 2 + t173 ^ 2) + Icges(2,3) + Icges(3,1) - t452) * qJDD(1) + (t414 + t415 + t418) * t289; (-m(3) - m(5)) * t384 - t360 * m(5) - t226 * t40 * m(3) + (-t384 + t392) * m(4) + 0.2e1 * (m(3) * t39 + m(5) * t9) * t377; t437 * t379 + t436 * t378 + (qJD(1) * t414 + qJD(3) * t393 - qJDD(1) * t405 + t148 * t409 + t149 * t410) * t377 + (qJD(1) * t413 + qJD(3) * t394 + qJDD(1) * t406 + t148 * t411 + t149 * t412) * t226 / 0.2e1 - ((t223 * t448 + t225 * t447) * qJD(3) + (t445 * t223 - t446 * t225) * qJD(1)) * qJD(1) / 0.2e1 + (t416 * t226 + t415 * t224 + (t224 * t408 + t226 * t407) * qJD(1)) * qJD(1) / 0.2e1 + (t224 * t407 - t226 * t408) * qJDD(1) / 0.2e1 - t418 * t312 / 0.2e1 + t417 * t311 / 0.2e1 + ((-t309 * t424 - t426) * t224 + ((t224 * t430 + t422) * qJD(3) + t432) * t226) * t290 + ((-t224 * t410 + t226 * t409) * qJD(1) + t393) * t289 + ((t308 * t430 - t426) * t226 + ((-t226 * t424 - t422) * qJD(3) - t432) * t224) * t288 + ((-t224 * t412 + t226 * t411) * qJD(1) + t394) * t287 + (-((t30 * t270 + (-pkin(3) * t338 + t132) * t28) * t226 + (-t397 * t28 + t29 * t429) * t224) * qJD(3) - g(1) * t397 - g(2) * (-t338 * t375 + t195) + g(3) * t389 + t9 * (t169 * t224 + t197) - t29 * pkin(3) * t297 + t286 * t360 + (qJD(3) * t421 - t327 * t148 - t326 * t149) * t231 + t28 * t421 - (t224 * t29 - t226 * t30) * t144 + (t29 * t132 - t30 * t130 + t28 * (t224 * t326 - t226 * t327) + (t224 * t30 + t357) * t169) * qJD(1)) * m(5) + (-g(1) * t131 + g(2) * t133 + g(3) * t271 + (qJD(3) * t264 - t109 * t148 - t111 * t149) * t251 + t45 * ((-t109 * t226 + t111 * t224) * qJD(1) + t264) - t265 * t145 + ((t224 * t44 + t356) * qJD(1) + t392) * t170 - (t131 * t44 + t133 * t43) * qJD(1) - (t45 * (-t131 * t224 - t133 * t226) - t265 * t271) * qJD(3)) * m(4); (t224 * t419 + t226 * t420) * m(5);];
tau = t1;
