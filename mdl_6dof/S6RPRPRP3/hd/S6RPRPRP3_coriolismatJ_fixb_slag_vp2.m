% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:08:03
% EndTime: 2019-03-09 03:08:13
% DurationCPUTime: 5.89s
% Computational Cost: add. (11537->532), mult. (24523->700), div. (0->0), fcn. (24152->8), ass. (0->271)
t482 = Ifges(7,4) + Ifges(6,5);
t488 = Ifges(6,6) - Ifges(7,6);
t319 = sin(qJ(3));
t316 = sin(pkin(10));
t317 = cos(pkin(10));
t434 = sin(qJ(5));
t435 = cos(qJ(5));
t332 = t316 * t435 + t317 * t434;
t251 = t332 * t319;
t277 = t434 * t316 - t435 * t317;
t253 = t277 * t319;
t145 = -pkin(5) * t253 + qJ(6) * t251;
t456 = -t145 / 0.2e1;
t320 = cos(qJ(3));
t389 = t320 * t277;
t443 = -t389 / 0.2e1;
t487 = mrSges(6,3) + mrSges(7,2);
t375 = -cos(pkin(9)) * pkin(1) - pkin(2);
t275 = -pkin(3) * t320 - t319 * qJ(4) + t375;
t256 = t317 * t275;
t303 = sin(pkin(9)) * pkin(1) + pkin(7);
t393 = t317 * t319;
t140 = -pkin(8) * t393 + t256 + (-t303 * t316 - pkin(4)) * t320;
t294 = t320 * t303;
t184 = t316 * t275 + t317 * t294;
t396 = t316 * t319;
t160 = -pkin(8) * t396 + t184;
t54 = t140 * t435 - t160 * t434;
t51 = t320 * pkin(5) - t54;
t486 = t51 + t54;
t430 = t319 * pkin(3);
t292 = -qJ(4) * t320 + t430;
t293 = t319 * t303;
t202 = t317 * t292 + t316 * t293;
t203 = t316 * t292 - t293 * t317;
t340 = -t202 * t316 + t203 * t317;
t485 = m(7) + m(6);
t408 = t320 * mrSges(6,1);
t479 = t253 * mrSges(6,3);
t208 = -t408 + t479;
t484 = -t208 / 0.2e1;
t252 = t332 * t320;
t445 = t252 / 0.2e1;
t315 = t317 ^ 2;
t483 = t315 / 0.2e1;
t438 = t319 / 0.2e1;
t398 = t277 * qJ(6);
t431 = t332 * pkin(5);
t350 = -t398 - t431;
t432 = m(7) * t350;
t481 = Ifges(7,2) + Ifges(6,3);
t480 = t251 * mrSges(6,3);
t390 = t320 * qJ(6);
t372 = t434 * t140;
t374 = t435 * t160;
t55 = t374 + t372;
t50 = t55 - t390;
t459 = mrSges(6,3) / 0.2e1;
t460 = mrSges(7,2) / 0.2e1;
t378 = t459 + t460;
t478 = t378 * t251;
t186 = mrSges(7,1) * t332 + t277 * mrSges(7,3);
t187 = mrSges(6,1) * t332 - t277 * mrSges(6,2);
t477 = t186 + t187;
t312 = t320 * mrSges(7,3);
t205 = -t251 * mrSges(7,2) - t312;
t206 = mrSges(6,2) * t320 - t480;
t476 = t206 + t205;
t209 = mrSges(7,1) * t320 - t253 * mrSges(7,2);
t475 = -t208 + t209;
t474 = -t277 * t482 - t488 * t332;
t473 = -t251 * t482 + t488 * t253;
t311 = m(7) * qJ(6) + mrSges(7,3);
t471 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t470 = -mrSges(6,2) + t311;
t411 = t389 * mrSges(7,3);
t413 = t389 * mrSges(6,2);
t414 = t252 * mrSges(7,1);
t415 = t252 * mrSges(6,1);
t469 = -t411 / 0.2e1 + t413 / 0.2e1 - t414 / 0.2e1 - t415 / 0.2e1;
t468 = 0.2e1 * m(7);
t467 = 2 * qJD(3);
t466 = m(5) / 0.2e1;
t465 = -m(6) / 0.2e1;
t464 = m(6) / 0.2e1;
t463 = -m(7) / 0.2e1;
t462 = m(7) / 0.2e1;
t461 = mrSges(6,1) / 0.2e1;
t392 = t317 * t320;
t163 = t319 * pkin(4) - pkin(8) * t392 + t202;
t395 = t316 * t320;
t180 = -pkin(8) * t395 + t203;
t66 = t163 * t435 - t180 * t434;
t60 = -t319 * pkin(5) - t66;
t458 = t60 / 0.2e1;
t457 = Ifges(7,5) * t443 + Ifges(7,6) * t438 + Ifges(7,3) * t445;
t455 = t145 / 0.2e1;
t146 = mrSges(7,1) * t251 + mrSges(7,3) * t253;
t454 = t146 / 0.2e1;
t453 = t186 / 0.2e1;
t452 = t187 / 0.2e1;
t428 = pkin(8) + qJ(4);
t291 = t428 * t317;
t363 = t428 * t316;
t200 = t291 * t434 + t363 * t435;
t451 = -t200 / 0.2e1;
t201 = t291 * t435 - t363 * t434;
t450 = -t201 / 0.2e1;
t409 = t319 * mrSges(7,1);
t412 = t389 * mrSges(7,2);
t211 = -t409 - t412;
t449 = t211 / 0.2e1;
t233 = t251 * mrSges(6,2);
t448 = -t233 / 0.2e1;
t234 = t253 * mrSges(7,1);
t447 = -t234 / 0.2e1;
t446 = -t252 / 0.2e1;
t444 = -t253 / 0.2e1;
t442 = -t277 / 0.2e1;
t441 = t332 / 0.2e1;
t440 = -t316 / 0.2e1;
t439 = t317 / 0.2e1;
t437 = -t320 / 0.2e1;
t433 = m(5) * t303;
t108 = m(7) * t253;
t427 = -t50 + t55;
t426 = mrSges(5,1) * t316;
t425 = mrSges(5,2) * t317;
t424 = Ifges(5,4) * t316;
t423 = Ifges(5,4) * t317;
t422 = Ifges(6,4) * t253;
t421 = Ifges(6,4) * t332;
t420 = Ifges(5,5) * t317;
t419 = Ifges(7,5) * t251;
t418 = Ifges(7,5) * t277;
t417 = Ifges(5,2) * t316;
t416 = Ifges(5,6) * t316;
t410 = t277 * mrSges(7,2);
t407 = -mrSges(5,1) * t317 + mrSges(5,2) * t316 - mrSges(4,1);
t262 = pkin(4) * t396 + t293;
t352 = pkin(5) * t251 + qJ(6) * t253;
t80 = t352 + t262;
t25 = m(7) * (t80 * t253 - t320 * t50) + t253 * t146 - t320 * t205;
t403 = t25 * qJD(1);
t402 = t251 * t277;
t400 = t253 * t332;
t286 = t319 * mrSges(5,1) - mrSges(5,3) * t392;
t397 = t316 * t286;
t284 = -t319 * mrSges(5,2) - mrSges(5,3) * t395;
t394 = t317 * t284;
t391 = t319 * t320;
t67 = t434 * t163 + t435 * t180;
t369 = -t400 / 0.2e1;
t370 = t402 / 0.2e1;
t388 = mrSges(7,2) * t370 + mrSges(6,3) * t369;
t263 = pkin(4) * t395 + t294;
t314 = t316 ^ 2;
t383 = t314 + t315;
t382 = qJD(3) * t320;
t381 = m(5) * t438;
t380 = m(7) * t455;
t379 = -t432 / 0.2e1;
t377 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t376 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t304 = -pkin(4) * t317 - pkin(3);
t368 = t332 * t437;
t237 = Ifges(7,5) * t253;
t131 = -Ifges(7,6) * t320 + Ifges(7,3) * t251 - t237;
t133 = -Ifges(6,2) * t251 - Ifges(6,6) * t320 - t422;
t367 = t131 / 0.2e1 - t133 / 0.2e1;
t135 = -Ifges(7,1) * t253 - t320 * Ifges(7,4) + t419;
t240 = Ifges(6,4) * t251;
t137 = -Ifges(6,1) * t253 - t320 * Ifges(6,5) - t240;
t366 = t135 / 0.2e1 + t137 / 0.2e1;
t365 = t206 / 0.2e1 + t205 / 0.2e1;
t364 = t484 + t209 / 0.2e1;
t362 = t383 * mrSges(5,3);
t360 = t383 * qJ(4);
t359 = t200 * t252 - t201 * t389;
t355 = t378 * t253;
t351 = -pkin(5) * t252 - qJ(6) * t389;
t349 = pkin(5) * t277 - qJ(6) * t332;
t134 = -Ifges(6,4) * t389 - Ifges(6,2) * t252 + Ifges(6,6) * t319;
t136 = -Ifges(7,1) * t389 + Ifges(7,4) * t319 + Ifges(7,5) * t252;
t138 = -Ifges(6,1) * t389 - Ifges(6,4) * t252 + Ifges(6,5) * t319;
t147 = t411 + t414;
t148 = -t413 + t415;
t183 = -t294 * t316 + t256;
t204 = -mrSges(7,2) * t252 + mrSges(7,3) * t319;
t207 = -mrSges(6,2) * t319 - mrSges(6,3) * t252;
t210 = t319 * mrSges(6,1) + mrSges(6,3) * t389;
t247 = Ifges(5,6) * t319 + (-t417 + t423) * t320;
t248 = Ifges(5,5) * t319 + (t317 * Ifges(5,1) - t424) * t320;
t266 = (t425 + t426) * t320;
t283 = t320 * mrSges(5,2) - mrSges(5,3) * t396;
t285 = -t320 * mrSges(5,1) - mrSges(5,3) * t393;
t328 = t252 * t376 + t377 * t389;
t59 = qJ(6) * t319 + t67;
t81 = -t351 + t263;
t4 = t81 * t146 + t80 * t147 + t262 * t148 + t183 * t286 + t184 * t284 + t202 * t285 + t203 * t283 + t50 * t204 + t59 * t205 + t67 * t206 + t55 * t207 + t66 * t208 + t60 * t209 + t54 * t210 + t51 * t211 - t366 * t389 + t367 * t252 - (t263 * mrSges(6,2) + t136 / 0.2e1 + t138 / 0.2e1) * t253 + (t263 * mrSges(6,1) - t134 / 0.2e1 + t457) * t251 + m(6) * (t262 * t263 + t54 * t66 + t55 * t67) + m(7) * (t50 * t59 + t51 * t60 + t80 * t81) + m(5) * (t183 * t202 + t184 * t203) + (t303 * t266 + t248 * t439 + t247 * t440 + t375 * mrSges(4,1) + (-Ifges(4,4) + t420 / 0.2e1 - t416 / 0.2e1) * t319 - t377 * t253 - t376 * t251) * t319 + (t375 * mrSges(4,2) + (Ifges(5,1) * t483 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) + (t425 + t433) * t303 + (t303 * mrSges(5,1) - t423 + t417 / 0.2e1) * t316 - t481) * t319 + t328 + (Ifges(4,4) + t416 - t420) * t320) * t320;
t334 = t283 * t439 + t285 * t440;
t336 = t426 / 0.2e1 + t425 / 0.2e1;
t342 = -t183 * t316 + t184 * t317;
t6 = -t365 * t389 - (t207 / 0.2e1 + t204 / 0.2e1) * t253 + t364 * t252 + (-t210 / 0.2e1 + t449) * t251 + (t394 / 0.2e1 - t397 / 0.2e1 + t454 + t251 * t461 + mrSges(6,2) * t444 + t336 * t319) * t319 + (-t147 / 0.2e1 - t148 / 0.2e1 - t266 / 0.2e1 + t334) * t320 + (t60 * t251 + t51 * t252 - t59 * t253 + t80 * t319 - t320 * t81 - t389 * t50) * t462 + (-t66 * t251 - t54 * t252 - t67 * t253 + t262 * t319 - t263 * t320 - t389 * t55) * t464 + ((t342 - t294) * t320 + (t340 + t293) * t319) * t466;
t348 = t4 * qJD(1) + t6 * qJD(2);
t149 = -Ifges(7,3) * t253 - t419;
t150 = Ifges(6,2) * t253 - t240;
t151 = -Ifges(7,1) * t251 - t237;
t152 = -Ifges(6,1) * t251 + t422;
t5 = -t80 * t234 - t262 * t233 - (-t50 * mrSges(7,2) + t262 * mrSges(6,1) + t152 / 0.2e1 + t151 / 0.2e1 + t367) * t253 + (t80 * mrSges(7,3) - t51 * mrSges(7,2) + t149 / 0.2e1 - t150 / 0.2e1 - t366) * t251 + (m(7) * t80 + t146) * t145 + (m(7) * t51 + t475 + t479) * t55 + (m(7) * t50 + t476 + t480) * t54 + t473 * t437;
t7 = (t447 + t380 + t448) * t320 - (t486 * t463 + t408 / 0.2e1 - t355 - t364) * t253 + (t312 / 0.2e1 + t427 * t463 + t478 + t365) * t251;
t347 = t5 * qJD(1) - t7 * qJD(2);
t29 = m(5) * (-0.1e1 + t383) * t391 + t485 * (t251 * t252 + t253 * t389 - t391);
t346 = t6 * qJD(1) + t29 * qJD(2);
t345 = t7 * qJD(1);
t13 = -t475 * t253 - t476 * t251 + m(7) * (-t50 * t251 - t51 * t253) + m(6) * (-t55 * t251 + t54 * t253) + (m(5) * (-t183 * t317 - t184 * t316) - t317 * t285 - t316 * t283) * t319;
t344 = qJD(1) * t13;
t333 = t253 * mrSges(6,1) - t251 * mrSges(7,3) + t233 + t234;
t32 = t456 * t468 + t333;
t39 = (-t431 / 0.4e1 - t398 / 0.4e1 + t350 / 0.4e1) * t468 - t477;
t343 = qJD(1) * t32 + qJD(3) * t39;
t341 = -t200 * t253 - t201 * t251;
t27 = -t312 + (t374 / 0.4e1 + t372 / 0.4e1 - t390 / 0.2e1 - t55 / 0.4e1) * t468;
t339 = qJD(1) * t27 + qJD(5) * t311;
t167 = m(7) * t332;
t338 = -qJD(1) * t108 + qJD(3) * t167;
t337 = m(7) * t458 - t409 / 0.2e1;
t165 = t304 + t349;
t188 = mrSges(7,1) * t277 - mrSges(7,3) * t332;
t190 = Ifges(7,3) * t332 - t418;
t271 = Ifges(7,5) * t332;
t191 = Ifges(7,3) * t277 + t271;
t274 = Ifges(6,4) * t277;
t192 = -Ifges(6,2) * t332 - t274;
t193 = -Ifges(6,2) * t277 + t421;
t194 = -Ifges(7,1) * t277 + t271;
t195 = Ifges(7,1) * t332 + t418;
t196 = -Ifges(6,1) * t277 - t421;
t197 = Ifges(6,1) * t332 - t274;
t10 = -t350 * t188 + t304 * t187 - (t193 / 0.2e1 - t194 / 0.2e1 - t196 / 0.2e1 - t191 / 0.2e1) * t332 + (-t192 / 0.2e1 - t195 / 0.2e1 - t197 / 0.2e1 + t190 / 0.2e1) * t277 + (t186 - t432) * t165;
t327 = t351 * t462 + t469;
t11 = (t453 + t379 + t452) * t320 + t327;
t321 = -t365 * t200 - (t133 / 0.4e1 - t152 / 0.4e1 - t151 / 0.4e1 - t131 / 0.4e1) * t332 + (t149 / 0.4e1 - t150 / 0.4e1 - t137 / 0.4e1 - t135 / 0.4e1) * t277 + (-t253 * t450 + t251 * t451 - (-t55 / 0.2e1 + t50 / 0.2e1) * t332 + (-t54 / 0.2e1 - t51 / 0.2e1) * t277) * mrSges(7,2) + (mrSges(6,3) * t451 - t197 / 0.4e1 - t195 / 0.4e1 - t192 / 0.4e1 + t165 * mrSges(7,3) / 0.2e1 + t190 / 0.4e1) * t251 - (mrSges(6,3) * t450 + t304 * t461 - t193 / 0.4e1 + t196 / 0.4e1 + t194 / 0.4e1 + t191 / 0.4e1) * t253 + (t145 * t165 + t200 * t427 - t350 * t80) * t462 + t188 * t455 + t165 * t447 - t350 * t454 + t262 * t452 + t304 * t448 + t80 * t453 - t474 * t320 / 0.4e1 + (t462 * t486 + t364) * t201;
t323 = (-pkin(5) * t60 + qJ(6) * t59) * t463 + pkin(5) * t449 - qJ(6) * t204 / 0.2e1 - t59 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t458 - t66 * mrSges(6,1) / 0.2e1 + t67 * mrSges(6,2) / 0.2e1;
t2 = t321 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t319 + t323 + t328;
t331 = t2 * qJD(1) - t11 * qJD(2) + t10 * qJD(3);
t24 = m(5) * t360 + t362 + (t277 ^ 2 + t332 ^ 2) * t487 + t485 * (t200 * t332 - t201 * t277);
t326 = (t462 + t464) * (t251 * t332 + t277 * t253);
t34 = (t463 + t465 + (t314 / 0.2e1 + t483 - 0.1e1 / 0.2e1) * m(5)) * t319 + t326;
t322 = -m(5) * t342 / 0.2e1 + (-t277 * t55 - t332 * t54 + t341) * t465 + (-t277 * t50 + t332 * t51 + t341) * t463 - t334;
t324 = (t433 / 0.2e1 + t336) * t320 + t263 * t464 + t81 * t462 - t469;
t8 = -(-t355 + t364) * t332 + (t365 - t478) * t277 + t322 + t324;
t330 = -qJD(1) * t8 + qJD(2) * t34 + qJD(3) * t24;
t142 = (-t368 + t446) * m(7);
t325 = (t165 * t253 - t201 * t320 - t332 * t80) * t463 + t188 * t444 + t146 * t441;
t20 = 0.2e1 * mrSges(7,2) * t443 + t325 + t337;
t48 = (m(7) * t165 + t188) * t332;
t329 = qJD(1) * t20 - qJD(2) * t142 + qJD(3) * t48;
t189 = mrSges(6,1) * t277 + mrSges(6,2) * t332;
t141 = (-t368 + t445) * m(7);
t74 = m(7) * t201 - t410;
t71 = t350 * t462 + t379;
t52 = m(7) * t456 + t380;
t33 = t381 * t383 + t485 * t438 + t326 + t381;
t26 = 0.2e1 * t462 * t50 + t205;
t21 = t389 * t460 - t412 / 0.2e1 - t325 + t337;
t12 = t400 * t459 - mrSges(7,2) * t402 / 0.2e1 + t327 + t388 + (-t432 + t477) * t437;
t9 = mrSges(7,2) * t369 + mrSges(6,3) * t370 + t209 * t441 + t332 * t484 + t442 * t476 - t322 + t324 + t388;
t3 = t6 * qJD(3) - t7 * qJD(5);
t1 = Ifges(6,6) * t446 + Ifges(7,6) * t445 + t438 * t481 + t443 * t482 + t321 - t323;
t14 = [qJD(3) * t4 + qJD(4) * t13 + qJD(5) * t5 + qJD(6) * t25, t3, t9 * qJD(4) + t1 * qJD(5) + t21 * qJD(6) + ((Ifges(5,1) * t316 + t423) * t439 + (Ifges(5,2) * t317 + t424) * t440 + Ifges(4,5) + t407 * t303) * t382 + ((-pkin(3) * t294 + qJ(4) * t340) * t466 + (-t200 * t66 + t201 * t67 + t263 * t304) * t464 + (t165 * t81 + t200 * t60 + t201 * t59) * t462) * t467 + t348 + (t60 * t332 * mrSges(7,2) - Ifges(4,6) * t319 + t316 * t248 / 0.2e1 + t304 * t148 - pkin(3) * t266 + t263 * t189 + t81 * t188 + t165 * t147 + t247 * t439 + t134 * t442 + t191 * t445 + t193 * t446 + mrSges(4,2) * t293 + t277 * t457 + (Ifges(5,5) * t316 + Ifges(5,6) * t317 - t277 * t488 + t482 * t332) * t438 - t59 * t410 + (t197 + t195) * t443 + (t138 + t136) * t441 + (t207 + t204) * t201 + (-t210 + t211) * t200 + (t394 - t397) * qJ(4) + (-t277 * t67 - t332 * t66) * mrSges(6,3) + t340 * mrSges(5,3)) * qJD(3), qJD(3) * t9 + qJD(5) * t52 + t344, t1 * qJD(3) + t52 * qJD(4) + (t352 * mrSges(7,2) + t470 * t54 + t471 * t55 + t473) * qJD(5) + t26 * qJD(6) + t347, t21 * qJD(3) + t26 * qJD(5) + t403; t3, qJD(3) * t29, t33 * qJD(4) + t12 * qJD(5) + t141 * qJD(6) + (-mrSges(4,2) + t362) * t382 + ((t165 * t319 + t359) * t462 + (t304 * t319 + t359) * t464 + (t320 * t360 - t430) * t466) * t467 + t346 + ((t188 + t189 + t407) * t319 + t487 * (t252 * t332 + t277 * t389)) * qJD(3), qJD(3) * t33, t12 * qJD(3) + t333 * qJD(5) + (qJD(5) * t456 + qJD(6) * t444) * t468 - t345, t141 * qJD(3) - qJD(5) * t108; -qJD(4) * t8 + qJD(5) * t2 - qJD(6) * t20 - t348, qJD(4) * t34 - qJD(5) * t11 + qJD(6) * t142 - t346, qJD(4) * t24 + qJD(5) * t10 - qJD(6) * t48, qJD(5) * t71 + t330, t71 * qJD(4) + (t349 * mrSges(7,2) - t200 * t470 + t201 * t471 + t474) * qJD(5) + t74 * qJD(6) + t331, qJD(5) * t74 - t329; qJD(3) * t8 - qJD(5) * t32 + qJD(6) * t108 - t344, -qJD(3) * t34, -qJD(5) * t39 - qJD(6) * t167 - t330, 0, -t343, -t338; -qJD(3) * t2 + qJD(4) * t32 + qJD(6) * t27 - t347, qJD(3) * t11 + t345, qJD(4) * t39 - t331, t343, t311 * qJD(6), t339; t20 * qJD(3) - t108 * qJD(4) - t27 * qJD(5) - t403, -t142 * qJD(3), qJD(4) * t167 + t329, t338, -t339, 0;];
Cq  = t14;
