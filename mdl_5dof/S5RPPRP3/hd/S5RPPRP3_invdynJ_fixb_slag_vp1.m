% Calculate vector of inverse dynamics joint torques for
% S5RPPRP3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:51:01
% DurationCPUTime: 15.33s
% Computational Cost: add. (6795->483), mult. (7991->576), div. (0->0), fcn. (6027->6), ass. (0->262)
t485 = Icges(5,3) + Icges(6,3);
t235 = qJ(1) + pkin(7);
t232 = sin(t235);
t233 = cos(t235);
t239 = cos(qJ(4));
t237 = sin(qJ(4));
t373 = Icges(6,4) * t237;
t275 = Icges(6,2) * t239 + t373;
t105 = -Icges(6,6) * t232 + t233 * t275;
t375 = Icges(5,4) * t237;
t276 = Icges(5,2) * t239 + t375;
t107 = -Icges(5,6) * t232 + t233 * t276;
t476 = t105 + t107;
t372 = Icges(6,4) * t239;
t277 = Icges(6,1) * t237 + t372;
t109 = -Icges(6,5) * t232 + t233 * t277;
t374 = Icges(5,4) * t239;
t278 = Icges(5,1) * t237 + t374;
t111 = -Icges(5,5) * t232 + t233 * t278;
t475 = t109 + t111;
t273 = Icges(6,5) * t237 + Icges(6,6) * t239;
t274 = Icges(5,5) * t237 + Icges(5,6) * t239;
t484 = t273 + t274;
t104 = Icges(6,6) * t233 + t232 * t275;
t106 = Icges(5,6) * t233 + t232 * t276;
t477 = t104 + t106;
t366 = t232 * t239;
t196 = Icges(6,4) * t366;
t367 = t232 * t237;
t370 = Icges(6,5) * t233;
t108 = Icges(6,1) * t367 + t196 + t370;
t197 = Icges(5,4) * t366;
t371 = Icges(5,5) * t233;
t110 = Icges(5,1) * t367 + t197 + t371;
t481 = t108 + t110;
t482 = t484 * t232 + t485 * t233;
t182 = -Icges(6,2) * t237 + t372;
t184 = -Icges(5,2) * t237 + t374;
t480 = t182 + t184;
t186 = Icges(6,1) * t239 - t373;
t188 = Icges(5,1) * t239 - t375;
t483 = t186 + t188;
t463 = t475 * t237 + t476 * t239;
t444 = -t485 * t232 + t484 * t233;
t458 = t481 * t237 + t477 * t239;
t479 = t277 + t278;
t445 = t482 * t232;
t470 = t483 * t237 + t480 * t239;
t478 = t463 * t233;
t430 = t482 * t233 + t477 * t366 + t481 * t367;
t429 = -t444 * t233 - t476 * t366 - t475 * t367;
t428 = -t458 * t233 + t445;
t427 = -t444 * t232 + t478;
t178 = Icges(6,5) * t239 - Icges(6,6) * t237;
t122 = t232 * t178;
t180 = Icges(5,5) * t239 - Icges(5,6) * t237;
t124 = t232 * t180;
t474 = t122 + t124;
t473 = (t275 + t276) * qJD(4);
t472 = t479 * qJD(4);
t471 = t480 * t237 - t483 * t239;
t368 = t180 * t233;
t369 = t178 * t233;
t469 = -t368 - t369;
t447 = t470 * t232 - t469;
t446 = t470 * t233 - t474;
t127 = t182 * t233;
t129 = t184 * t233;
t468 = (-t127 - t129) * qJD(4) + t477 * qJD(1);
t331 = qJD(4) * t232;
t467 = -t476 * qJD(1) - t480 * t331;
t131 = t186 * t233;
t133 = t188 * t233;
t466 = (-t131 - t133) * qJD(4) + (t479 * t232 + t370 + t371) * qJD(1);
t130 = t186 * t232;
t132 = t188 * t232;
t465 = (t130 + t132) * qJD(4) + t475 * qJD(1);
t464 = t470 * qJD(1) - t484 * qJD(4);
t462 = t427 * t232 + t428 * t233;
t461 = t429 * t232 + t430 * t233;
t460 = t473 * t239 + t472 * t237 + t471 * qJD(4) + (t178 + t180) * qJD(1);
t240 = cos(qJ(1));
t234 = t240 * pkin(1);
t425 = t476 * t237 - t475 * t239;
t459 = t447 * qJD(1);
t150 = rSges(3,1) * t232 + rSges(3,2) * t233;
t238 = sin(qJ(1));
t393 = pkin(1) * t238;
t138 = -t150 - t393;
t426 = t477 * t237 - t481 * t239;
t457 = t446 * qJD(1);
t215 = qJD(3) * t233;
t420 = t233 * pkin(2) + t232 * qJ(3);
t119 = qJD(1) * t420 - t215;
t324 = qJD(1) * qJD(4);
t142 = qJDD(4) * t232 + t233 * t324;
t381 = rSges(6,2) * t239;
t287 = rSges(6,1) * t237 + t381;
t161 = t287 * qJD(4);
t193 = rSges(6,1) * t239 - rSges(6,2) * t237;
t217 = t233 * qJ(3);
t392 = pkin(2) * t232;
t148 = -t217 + t392;
t296 = -pkin(6) * t232 - t393;
t224 = t232 * rSges(6,3);
t365 = t233 * t237;
t236 = -qJ(5) - pkin(6);
t386 = pkin(6) + t236;
t350 = pkin(4) * t365 + t232 * t386 + t233 * t287 - t224;
t258 = t296 + t350;
t251 = -t148 + t258;
t325 = qJD(1) * qJD(3);
t343 = qJDD(3) * t232 + t233 * t325;
t242 = qJD(1) ^ 2;
t454 = t242 * t234;
t261 = t343 - t454;
t302 = -pkin(6) * t242 + qJDD(5);
t327 = qJD(5) * t232;
t363 = t237 * qJD(4) ^ 2;
t227 = t233 * rSges(6,3);
t328 = qJD(4) * t239;
t313 = t233 * t328;
t301 = -pkin(4) * t313 + t327;
t330 = qJD(4) * t233;
t317 = t193 * t330;
t204 = pkin(4) * t367;
t421 = -t386 * t233 + t204;
t385 = -t317 + t301 + (t232 * t287 + t227 + t421) * qJD(1);
t10 = -t161 * t331 + t142 * t193 + t302 * t233 + (t142 * t239 - t232 * t363) * pkin(4) + (-t119 - t327 - t385) * qJD(1) + t251 * qJDD(1) + t261;
t456 = t10 - g(1);
t143 = qJDD(4) * t233 - t232 * t324;
t390 = pkin(6) * t233;
t209 = qJDD(1) * t390;
t231 = qJDD(1) * t234;
t298 = -t242 * t393 + t231;
t333 = qJD(1) * t232;
t214 = qJD(3) * t232;
t332 = qJD(1) * t233;
t342 = qJ(3) * t332 + t214;
t415 = qJD(1) * (-pkin(2) * t333 + t342) + qJDD(1) * t420 + t232 * t325;
t256 = t298 + t415;
t303 = -pkin(4) * t239 - t193;
t112 = rSges(6,1) * t367 + rSges(6,2) * t366 + t227;
t351 = t112 + t421;
t329 = qJD(4) * t237;
t318 = t239 * t332;
t319 = t237 * t332;
t315 = t232 * t328;
t349 = pkin(4) * t315 + qJD(5) * t233;
t434 = t315 + t319;
t417 = rSges(6,1) * t434 + rSges(6,2) * t318 + pkin(4) * t319 + t236 * t333 + t349;
t384 = (-rSges(6,2) * t329 - rSges(6,3) * qJD(1)) * t232 + pkin(6) * t333 + t417;
t11 = t209 + t302 * t232 + t303 * t143 + t351 * qJDD(1) + t384 * qJD(1) + (pkin(4) * t363 + qJD(1) * qJD(5) + qJD(4) * t161 - qJDD(3)) * t233 + t256;
t455 = t11 - g(2);
t453 = t461 * qJD(4) + t459;
t452 = t462 * qJD(4) - t457;
t451 = -t458 * qJD(4) + t467 * t237 + t465 * t239;
t450 = t463 * qJD(4) - t468 * t237 + t466 * t239;
t449 = t464 * t232 + t460 * t233;
t448 = -t460 * t232 + t464 * t233;
t442 = t444 * qJD(1);
t441 = t482 * qJD(1);
t440 = t233 ^ 2;
t152 = -rSges(4,2) * t233 + t232 * rSges(4,3);
t320 = t234 + t420;
t99 = t152 + t320;
t439 = t458 * qJD(1) + t474 * qJD(4) + t442;
t438 = -t463 * qJD(1) + t469 * qJD(4) + t441;
t437 = t426 * qJD(4) - t465 * t237 + t467 * t239 + t441;
t436 = t425 * qJD(4) + t466 * t237 + t468 * t239 + t442;
t435 = -pkin(4) * t237 - t287;
t346 = t182 + t277;
t347 = -t275 + t186;
t433 = (t237 * t346 - t239 * t347) * qJD(1);
t344 = t184 + t278;
t345 = -t276 + t188;
t432 = (t237 * t344 - t239 * t345) * qJD(1);
t134 = rSges(6,1) * t366 - rSges(6,2) * t367;
t205 = pkin(4) * t366;
t422 = t134 + t205;
t153 = t233 * rSges(3,1) - rSges(3,2) * t232;
t139 = t153 + t234;
t419 = t439 * t440 + (t436 * t232 + (-t437 + t438) * t233) * t232;
t418 = t437 * t440 + (t438 * t232 + (-t436 + t439) * t233) * t232;
t394 = rSges(6,1) + pkin(4);
t412 = t237 * t394 + t381;
t352 = t111 + t129;
t356 = t107 - t133;
t405 = t237 * t356 - t239 * t352;
t353 = -Icges(5,2) * t367 + t110 + t197;
t357 = t106 - t132;
t404 = t237 * t357 - t239 * t353;
t354 = t109 + t127;
t358 = t105 - t131;
t403 = t237 * t358 - t239 * t354;
t355 = -Icges(6,2) * t367 + t108 + t196;
t359 = t104 - t130;
t402 = t237 * t359 - t239 * t355;
t401 = -pkin(2) - pkin(6);
t399 = t142 / 0.2e1;
t398 = t143 / 0.2e1;
t397 = t232 / 0.2e1;
t396 = -t233 / 0.2e1;
t387 = -pkin(2) + t236;
t380 = rSges(4,3) * t233;
t288 = rSges(5,1) * t237 + rSges(5,2) * t239;
t162 = t288 * qJD(4);
t194 = rSges(5,1) * t239 - rSges(5,2) * t237;
t115 = -t232 * rSges(5,3) + t233 * t288;
t260 = t115 + t296;
t257 = -t148 + t260;
t295 = -t234 - t390;
t137 = t194 * t233;
t228 = t233 * rSges(5,3);
t75 = -qJD(4) * t137 + (t232 * t288 + t228) * qJD(1);
t16 = -t162 * t331 + t142 * t194 + t295 * t242 + (-t119 - t75) * qJD(1) + t257 * qJDD(1) + t343;
t379 = t16 * t232;
t113 = rSges(5,1) * t367 + rSges(5,2) * t366 + t228;
t321 = rSges(5,1) * t434 + rSges(5,2) * t318;
t77 = (-rSges(5,2) * t329 - rSges(5,3) * qJD(1)) * t232 + t321;
t17 = qJD(1) * t77 + qJDD(1) * t113 - t143 * t194 + t209 + t231 + t296 * t242 + (qJD(4) * t162 - qJDD(3)) * t233 + t415;
t378 = t17 * t233;
t141 = t194 * t331;
t45 = qJD(1) * t257 + t141 + t214;
t377 = t233 * t45;
t364 = t233 * t239;
t341 = rSges(4,2) * t333 + rSges(4,3) * t332;
t145 = qJD(1) * t148;
t340 = t214 - t145;
t335 = qJD(1) * t273;
t334 = qJD(1) * t274;
t326 = -m(4) - m(5) - m(6);
t323 = -rSges(5,3) + t401;
t316 = t232 * t329;
t314 = t233 * t329;
t309 = -t331 / 0.2e1;
t308 = t331 / 0.2e1;
t307 = -t330 / 0.2e1;
t306 = t330 / 0.2e1;
t305 = rSges(4,2) * t232 + t380 - t393;
t304 = t217 - t393;
t300 = t193 * t331 + t214 + t349;
t290 = -t215 + t301;
t195 = rSges(2,1) * t240 - rSges(2,2) * t238;
t192 = rSges(2,1) * t238 + rSges(2,2) * t240;
t259 = t420 - t295;
t46 = -t194 * t330 - t215 + (t113 + t259) * qJD(1);
t280 = t232 * t45 - t233 * t46;
t279 = -t232 * t77 + t233 * t75;
t266 = -t113 * t232 - t115 * t233;
t203 = rSges(6,2) * t365;
t136 = -rSges(6,1) * t364 + t203;
t135 = t194 * t232;
t47 = qJD(4) * t266 + qJD(2);
t34 = qJD(1) * t341 + qJDD(1) * t152 - qJDD(3) * t233 + t256;
t33 = (-t148 + t305) * qJDD(1) + (-qJD(1) * t152 - t119) * qJD(1) + t261;
t32 = -t317 + (t259 + t351) * qJD(1) + t290;
t31 = qJD(1) * t251 + t300;
t30 = qJD(2) + (-t232 * t351 - t233 * t350) * qJD(4);
t18 = qJD(4) * t279 - t113 * t142 - t115 * t143 + qJDD(2);
t1 = qJDD(2) - t350 * t143 - t351 * t142 + (-t384 * t232 + t385 * t233) * qJD(4);
t2 = [-m(2) * (-g(1) * t192 + g(2) * t195) - t446 * t142 / 0.2e1 + t425 * t399 + (((t427 + t430 - t478) * t233 + ((-t458 + t444) * t233 + t445 - t428 + t429) * t232) * qJD(4) + t459) * t309 + (-t470 * qJD(4) + t473 * t237 - t472 * t239) * qJD(1) + ((-t150 * t242 - g(2) + t298) * t139 + (-t454 + (-0.2e1 * t153 - t234 + t139) * t242 - g(1)) * t138) * m(3) + (t31 * (rSges(6,1) * t313 - rSges(6,2) * t314 - t290) + t32 * (-rSges(6,2) * t316 + t342 + t417) + ((-t238 * t32 - t240 * t31) * pkin(1) + t31 * (-rSges(6,3) + t387) * t233 + (t31 * (-qJ(3) + t435) + t32 * (-rSges(6,3) - pkin(2))) * t232) * qJD(1) - (t258 * qJD(1) - t145 + t300 - t31) * t32 + t455 * (-t233 * t236 + t112 + t204 + t320) + t456 * (t387 * t232 + t233 * t412 - t224 + t304)) * m(6) + (-(qJD(1) * t260 + t141 + t340 - t45) * t46 + t45 * (rSges(5,1) * t313 - rSges(5,2) * t314 + t215) + t46 * (-rSges(5,2) * t316 + t321 + t342) + ((-t238 * t46 - t240 * t45) * pkin(1) + t323 * t377 + (t45 * (-qJ(3) - t288) + t46 * t323) * t232) * qJD(1) + (t17 - g(2)) * (t113 + t320 + t390) + (t16 - g(1)) * (t232 * t401 + t115 + t304)) * m(5) + ((t33 - g(1)) * (t380 + (rSges(4,2) - pkin(2)) * t232 + t304) + (t34 - g(2)) * t99 + (-t340 + t341 + t342 + (-t305 - t392 - t393) * qJD(1)) * (qJD(1) * t99 - t215)) * m(4) + (-t426 + t447) * t398 + ((t444 * t232 ^ 2 + ((t458 + t444) * t233 - t445 + t429) * t233) * qJD(4) + t452 + t457) * t307 + (t448 + t451) * t306 + (m(3) * (t138 ^ 2 + t153 * t139) + m(2) * (t192 ^ 2 + t195 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1) - t471) * qJDD(1) + (t449 + t450 + t453) * t308; (m(3) + m(4)) * qJDD(2) + m(5) * t18 + m(6) * t1 + (-m(3) + t326) * g(3); t326 * (g(1) * t232 - g(2) * t233) + 0.2e1 * (t10 * t397 + t11 * t396) * m(6) + 0.2e1 * (t379 / 0.2e1 - t378 / 0.2e1) * m(5) + 0.2e1 * (t33 * t397 + t34 * t396) * m(4); t462 * t399 + t461 * t398 + (qJD(1) * t449 + t418 * qJD(4) - qJDD(1) * t446 + t427 * t142 + t428 * t143) * t397 + (qJD(1) * t448 + t419 * qJD(4) + qJDD(1) * t447 + t429 * t142 + t430 * t143) * t233 / 0.2e1 - ((((-t357 - t359) * t233 + (t356 + t358) * t232) * t239 + ((-t353 - t355) * t233 + (t352 + t354) * t232) * t237) * qJD(4) + ((-t344 - t346) * t239 + (-t345 - t347) * t237) * qJD(1)) * qJD(1) / 0.2e1 + (t451 * t233 + t450 * t232 + (t426 * t232 + t425 * t233) * qJD(1)) * qJD(1) / 0.2e1 + (t232 * t425 - t233 * t426) * qJDD(1) / 0.2e1 - t453 * t333 / 0.2e1 + t452 * t332 / 0.2e1 + ((-t331 * t369 - t335) * t232 + (t433 + (t402 * t233 + (t122 - t403) * t232) * qJD(4)) * t233 + (-t331 * t368 - t334) * t232 + (t432 + (t404 * t233 + (t124 - t405) * t232) * qJD(4)) * t233) * t309 + ((-t232 * t428 + t233 * t427) * qJD(1) + t418) * t308 + ((t122 * t330 - t335) * t233 + (-t433 + (t403 * t232 + (-t369 - t402) * t233) * qJD(4)) * t232 + (t124 * t330 - t334) * t233 + (-t432 + (t405 * t232 + (-t368 - t404) * t233) * qJD(4)) * t232) * t307 + ((-t232 * t430 + t429 * t233) * qJD(1) + t419) * t306 + (-g(1) * t422 - g(2) * (-t364 * t394 + t203) + g(3) * t412 - ((t32 * t287 + (-pkin(4) * t364 + t136) * t30) * t233 + (-t422 * t30 + t31 * t435) * t232) * qJD(4) + t10 * t205 + (t10 * t193 + t31 * (-pkin(4) * t329 - t161) - t1 * t351 - t30 * t384) * t232 + (-t1 * t350 + t11 * t303 + t32 * t161 + t30 * t385) * t233 + (t31 * t136 - t32 * t134 + (t32 * t193 + t30 * t350) * t232 + (t31 * t193 - t30 * t351) * t233) * qJD(1)) * m(6) + (t18 * t266 + t47 * ((-t113 * t233 + t115 * t232) * qJD(1) + t279) - t280 * t162 + (t379 - t378 + (t232 * t46 + t377) * qJD(1)) * t194 - g(1) * t135 + g(2) * t137 + g(3) * t288 - (t135 * t46 + t137 * t45) * qJD(1) - (t47 * (-t135 * t232 - t137 * t233) - t280 * t288) * qJD(4)) * m(5); (t455 * t232 + t233 * t456) * m(6);];
tau = t2;
