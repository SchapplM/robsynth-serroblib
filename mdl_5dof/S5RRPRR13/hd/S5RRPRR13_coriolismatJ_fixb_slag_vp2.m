% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:13
% EndTime: 2019-12-31 20:32:30
% DurationCPUTime: 7.95s
% Computational Cost: add. (16622->484), mult. (35629->683), div. (0->0), fcn. (38344->8), ass. (0->235)
t332 = sin(pkin(9));
t333 = cos(pkin(9));
t446 = sin(qJ(4));
t447 = cos(qJ(4));
t301 = -t446 * t332 + t447 * t333;
t302 = -t332 * t447 - t333 * t446;
t334 = sin(qJ(5));
t336 = cos(qJ(5));
t241 = t301 * t334 - t302 * t336;
t369 = t336 * t301 + t302 * t334;
t394 = Ifges(6,5) * t369 - Ifges(6,6) * t241;
t438 = pkin(7) + qJ(3);
t312 = t438 * t332;
t313 = t438 * t333;
t255 = -t312 * t446 + t313 * t447;
t205 = t301 * pkin(8) + t255;
t254 = -t447 * t312 - t313 * t446;
t346 = t302 * pkin(8) + t254;
t117 = t205 * t336 + t334 * t346;
t488 = -t205 * t334 + t336 * t346;
t523 = -t117 * mrSges(6,1) - t488 * mrSges(6,2);
t19 = t394 + t523;
t532 = t19 * qJD(5);
t432 = Ifges(6,4) * t241;
t136 = Ifges(6,2) * t369 + t432;
t518 = Ifges(6,1) * t369 - t432;
t531 = -t136 / 0.4e1 + t518 / 0.4e1;
t337 = cos(qJ(2));
t335 = sin(qJ(2));
t281 = t302 * t335;
t283 = t301 * t335;
t370 = t336 * t281 - t283 * t334;
t167 = mrSges(6,2) * t337 + mrSges(6,3) * t370;
t201 = t281 * t334 + t283 * t336;
t168 = -mrSges(6,1) * t337 - t201 * mrSges(6,3);
t400 = t332 * t335;
t308 = pkin(3) * t400 + t335 * pkin(6);
t244 = -pkin(4) * t281 + t308;
t325 = -pkin(3) * t333 - pkin(2);
t272 = -pkin(4) * t301 + t325;
t477 = -t117 / 0.2e1;
t193 = Ifges(6,4) * t370;
t503 = -Ifges(6,2) * t201 + t193;
t508 = t201 * mrSges(6,1);
t515 = t370 * mrSges(6,2) + t508;
t507 = t241 * mrSges(6,1);
t516 = t369 * mrSges(6,2) + t507;
t433 = Ifges(6,4) * t201;
t517 = Ifges(6,1) * t370 - t433;
t520 = t488 / 0.2e1;
t95 = Ifges(6,2) * t370 - t337 * Ifges(6,6) + t433;
t96 = Ifges(6,1) * t201 - t337 * Ifges(6,5) + t193;
t530 = t168 * t477 + t272 * t515 / 0.2e1 + t516 * t244 / 0.2e1 + t167 * t520 + (t96 / 0.4e1 + t503 / 0.4e1) * t369 + (t517 / 0.4e1 - t95 / 0.4e1) * t241;
t529 = t517 / 0.2e1;
t311 = -pkin(2) * t337 - t335 * qJ(3) - pkin(1);
t299 = t333 * t311;
t398 = t333 * t335;
t242 = -pkin(7) * t398 + t299 + (-pkin(6) * t332 - pkin(3)) * t337;
t397 = t333 * t337;
t266 = pkin(6) * t397 + t332 * t311;
t252 = -pkin(7) * t400 + t266;
t142 = t447 * t242 - t252 * t446;
t119 = -t283 * pkin(8) + t142;
t103 = -t337 * pkin(4) + t119;
t143 = t242 * t446 + t252 * t447;
t120 = t281 * pkin(8) + t143;
t409 = t120 * t336;
t43 = t103 * t334 + t409;
t50 = -t119 * t334 - t409;
t526 = t50 + t43;
t525 = t244 * t515;
t314 = t335 * pkin(2) - qJ(3) * t337;
t273 = pkin(6) * t400 + t333 * t314;
t243 = t335 * pkin(3) - pkin(7) * t397 + t273;
t274 = -pkin(6) * t398 + t332 * t314;
t399 = t332 * t337;
t253 = -pkin(7) * t399 + t274;
t144 = t447 * t243 - t253 * t446;
t284 = t301 * t337;
t118 = t335 * pkin(4) - t284 * pkin(8) + t144;
t145 = t446 * t243 + t447 * t253;
t456 = t284 / 0.2e1;
t282 = t302 * t337;
t459 = t282 / 0.2e1;
t350 = Ifges(5,5) * t456 + Ifges(5,6) * t459;
t199 = t282 * t336 - t284 * t334;
t363 = -t335 * mrSges(6,2) + t199 * mrSges(6,3);
t202 = t282 * t334 + t284 * t336;
t364 = t335 * mrSges(6,1) - t202 * mrSges(6,3);
t469 = t202 / 0.2e1;
t471 = t199 / 0.2e1;
t349 = Ifges(6,5) * t469 + Ifges(6,6) * t471;
t121 = pkin(8) * t282 + t145;
t352 = t118 * t334 + t121 * t336;
t353 = t118 * t336 - t121 * t334;
t521 = t335 / 0.2e1;
t368 = -t352 * mrSges(6,2) / 0.2e1 + t353 * mrSges(6,1) / 0.2e1 + Ifges(6,3) * t521 + t349;
t482 = m(6) * pkin(4);
t388 = t482 / 0.2e1;
t450 = t336 / 0.2e1;
t524 = Ifges(5,3) * t521 + t144 * mrSges(5,1) / 0.2e1 - t145 * mrSges(5,2) / 0.2e1 + (t334 ^ 2 + t336 ^ 2) * t118 * t388 + t350 + (t334 * t363 / 0.2e1 + t364 * t450) * pkin(4) + t368;
t384 = t507 / 0.2e1;
t385 = t508 / 0.2e1;
t494 = Ifges(6,5) * t370;
t510 = Ifges(6,6) * t201;
t396 = t494 - t510;
t230 = Ifges(6,4) * t369;
t139 = Ifges(6,1) * t241 + t230;
t502 = -Ifges(6,2) * t241 + t230;
t514 = t502 / 0.4e1 + t139 / 0.4e1;
t373 = t510 / 0.2e1 - t494 / 0.2e1;
t512 = 0.2e1 * mrSges(6,2);
t474 = -t201 / 0.2e1;
t511 = -t241 / 0.2e1;
t505 = m(4) * qJ(3) + mrSges(4,3);
t499 = t369 / 0.2e1;
t498 = t370 / 0.2e1;
t491 = Ifges(5,5) * t301 + Ifges(5,6) * t302 + t394;
t490 = Ifges(5,5) * t281 - Ifges(5,6) * t283 + t396;
t109 = -mrSges(6,1) * t370 + mrSges(6,2) * t201;
t133 = -mrSges(6,1) * t369 + mrSges(6,2) * t241;
t435 = Ifges(5,4) * t283;
t194 = Ifges(5,2) * t281 - t337 * Ifges(5,6) + t435;
t206 = t283 * mrSges(5,1) + t281 * mrSges(5,2);
t208 = Ifges(5,1) * t281 - t435;
t247 = -t302 * mrSges(5,1) + t301 * mrSges(5,2);
t434 = Ifges(5,4) * t302;
t249 = Ifges(5,2) * t301 - t434;
t250 = Ifges(5,1) * t301 + t434;
t260 = mrSges(5,2) * t337 + t281 * mrSges(5,3);
t416 = t283 * mrSges(5,3);
t261 = -mrSges(5,1) * t337 - t416;
t443 = pkin(4) * t302;
t386 = -t443 / 0.2e1;
t444 = pkin(4) * t283;
t387 = t444 / 0.2e1;
t410 = t120 * t334;
t42 = t103 * t336 - t410;
t448 = -t337 / 0.4e1;
t478 = -t488 / 0.2e1;
t481 = -t42 / 0.2e1;
t483 = m(6) / 0.2e1;
t51 = t119 * t336 - t410;
t486 = t109 * t386 + t133 * t387 + t491 * t448 + (t117 * t474 + t369 * t481 + t370 * t478 + t499 * t51 + t526 * t511) * mrSges(6,3) + ((-t244 * t302 + t272 * t283) * pkin(4) + t526 * t488 + (-t42 + t51) * t117) * t483 + t514 * t370 + t254 * t260 / 0.2e1 + t283 * t250 / 0.4e1 - t283 * t249 / 0.4e1 + t302 * t194 / 0.4e1 - t302 * t208 / 0.4e1 + t308 * t247 / 0.2e1 + t325 * t206 / 0.2e1 + (-t416 / 0.2e1 - t261 / 0.2e1) * t255 + t531 * t201 + t530;
t330 = t332 ^ 2;
t331 = t333 ^ 2;
t485 = m(4) / 0.2e1;
t484 = m(5) / 0.2e1;
t480 = -t95 / 0.2e1;
t473 = -t370 / 0.2e1;
t464 = t241 / 0.2e1;
t458 = -t283 / 0.2e1;
t457 = t283 / 0.2e1;
t455 = t301 / 0.2e1;
t453 = -t302 / 0.2e1;
t452 = t302 / 0.2e1;
t451 = -t334 / 0.2e1;
t449 = -t337 / 0.2e1;
t329 = t337 * pkin(6);
t442 = t42 * mrSges(6,2);
t441 = t43 * mrSges(6,1);
t440 = t50 * mrSges(6,1);
t439 = t51 * mrSges(6,2);
t437 = mrSges(4,1) * t332;
t436 = mrSges(4,2) * t333;
t430 = pkin(4) * qJD(4);
t424 = t199 * mrSges(6,1);
t421 = t202 * mrSges(6,2);
t417 = t282 * mrSges(5,1);
t415 = t284 * mrSges(5,2);
t278 = Ifges(5,4) * t281;
t195 = Ifges(5,1) * t283 - t337 * Ifges(5,5) + t278;
t309 = pkin(3) * t399 + t329;
t245 = -pkin(4) * t282 + t309;
t265 = -pkin(6) * t399 + t299;
t306 = t337 * mrSges(4,2) - mrSges(4,3) * t400;
t307 = -t337 * mrSges(4,1) - mrSges(4,3) * t398;
t351 = Ifges(4,5) * t333 - Ifges(4,6) * t332 - Ifges(3,4);
t359 = Ifges(6,4) * t202 + Ifges(6,2) * t199;
t360 = Ifges(5,4) * t284 + Ifges(5,2) * t282;
t361 = Ifges(6,1) * t202 + Ifges(6,4) * t199;
t362 = Ifges(5,1) * t284 + Ifges(5,4) * t282;
t365 = t421 - t424;
t366 = t415 - t417;
t367 = t436 + t437;
t3 = -t244 * t365 - t308 * t366 - t281 * t360 / 0.2e1 - t43 * t363 - t42 * t364 + (t265 * t397 + t266 * t399) * mrSges(4,3) + t361 * t474 + (t142 * t284 - t143 * t282) * mrSges(5,3) + t359 * t473 + t199 * t480 + (-0.2e1 * t367 * t329 - t265 * mrSges(4,1) + t266 * mrSges(4,2) + t143 * mrSges(5,2) - t142 * mrSges(5,1) - Ifges(5,5) * t283 - Ifges(5,6) * t281 - Ifges(6,5) * t201 - Ifges(6,6) * t370 + pkin(1) * mrSges(3,1) + (-Ifges(4,1) * t331 - m(4) * pkin(6) ^ 2 - Ifges(3,1) + Ifges(3,2) + Ifges(4,3) + Ifges(5,3) + Ifges(6,3) + (0.2e1 * Ifges(4,4) * t333 - Ifges(4,2) * t332) * t332) * t337 - t351 * t335) * t335 - m(4) * (t265 * t273 + t266 * t274) - m(5) * (t142 * t144 + t143 * t145 + t308 * t309) + t362 * t458 - t202 * t96 / 0.2e1 - t245 * t109 - t145 * t260 - t144 * t261 - t282 * t194 / 0.2e1 - t284 * t195 / 0.2e1 - t274 * t306 - t273 * t307 - t309 * (-mrSges(5,1) * t281 + mrSges(5,2) * t283) + (pkin(1) * mrSges(3,2) + t337 * t351 + t349 + t350) * t337 - t352 * t167 - t353 * t168 - m(6) * (t244 * t245 + t352 * t43 + t353 * t42);
t414 = t3 * qJD(1);
t207 = -Ifges(5,2) * t283 + t278;
t6 = m(6) * (t244 * t444 + t42 * t50 + t43 * t51) + t51 * t167 + t50 * t168 + t201 * t529 + t109 * t444 + t525 + t95 * t474 + t142 * t260 + t208 * t457 + t194 * t458 + t308 * t206 + (-t261 - t416) * t143 + (-t201 * t43 - t370 * t42) * mrSges(6,3) + (-t142 * mrSges(5,3) + t195 / 0.2e1 + t207 / 0.2e1) * t281 + t490 * t449 + (t96 + t503) * t498;
t413 = t6 * qJD(1);
t7 = t42 * t167 - t43 * t168 + t525 + t396 * t449 + (-t43 * mrSges(6,3) + t480 + t529) * t201 + (-t42 * mrSges(6,3) + t96 / 0.2e1 + t503 / 0.2e1) * t370;
t412 = t7 * qJD(1);
t16 = t370 * t167 - t201 * t168 + t281 * t260 - t283 * t261 + m(6) * (-t201 * t42 + t370 * t43) + m(5) * (-t142 * t283 + t143 * t281) + (m(4) * (-t265 * t333 - t266 * t332) - t332 * t306 - t333 * t307) * t335;
t411 = qJD(1) * t16;
t408 = t201 * t336;
t407 = t370 * t334;
t404 = t241 * t336;
t403 = t369 * t334;
t383 = -t254 * mrSges(5,3) / 0.2e1;
t293 = Ifges(5,4) * t301;
t248 = Ifges(5,2) * t302 + t293;
t251 = -Ifges(5,1) * t302 + t293;
t2 = t281 * t383 + (t248 + t251) * t281 / 0.4e1 + (t195 + t207) * t301 / 0.4e1 + t486 - t524;
t8 = t518 * t464 + t136 * t511 - t133 * t443 + t325 * t247 + t250 * t453 + t249 * t452 + (t251 / 0.2e1 + t248 / 0.2e1) * t301 + (t139 + t502) * t499 + (-m(6) * t443 + t516) * t272;
t358 = t2 * qJD(1) + t8 * qJD(2);
t14 = t272 * t516 + (t518 / 0.2e1 - t136 / 0.2e1) * t241 + (t139 / 0.2e1 + t502 / 0.2e1) * t369;
t340 = (mrSges(6,3) * t478 + t514) * t370 + (mrSges(6,3) * t477 + t531) * t201 + t394 * t448 + t530;
t5 = t340 - t368;
t357 = t5 * qJD(1) + t14 * qJD(2);
t339 = (-t201 * t511 + t370 * t499) * mrSges(6,3) + (t281 * t455 + t283 * t453) * mrSges(5,3) + (-t265 * t332 + t266 * t333) * t485 + (t142 * t302 + t143 * t301 - t254 * t283 + t255 * t281) * t484 + (t117 * t370 - t201 * t488 - t241 * t42 + t369 * t43) * t483 + t168 * t511 + t167 * t499 + t260 * t455 + t261 * t452 - t332 * t307 / 0.2e1 + t333 * t306 / 0.2e1;
t344 = -m(5) * t309 / 0.2e1 - m(6) * t245 / 0.2e1 + t424 / 0.2e1 - t421 / 0.2e1 + t417 / 0.2e1 - t415 / 0.2e1;
t10 = (-m(4) * pkin(6) / 0.2e1 - t436 / 0.2e1 - t437 / 0.2e1) * t337 + t339 + t344;
t21 = (t241 ^ 2 + t369 ^ 2) * mrSges(6,3) + (t301 ^ 2 + t302 ^ 2) * mrSges(5,3) + m(6) * (t117 * t369 - t241 * t488) + m(5) * (t254 * t302 + t255 * t301) + t505 * (t330 + t331);
t356 = -qJD(1) * t10 - qJD(2) * t21;
t27 = (t457 + t408 / 0.2e1 - t407 / 0.2e1) * t482 + t515 + t206;
t38 = (t453 + t404 / 0.2e1 - t403 / 0.2e1) * t482 + t516 + t247;
t355 = qJD(1) * t27 + qJD(2) * t38;
t48 = t498 * t512 + 0.2e1 * t385;
t52 = t499 * t512 + 0.2e1 * t384;
t354 = qJD(1) * t48 + qJD(2) * t52;
t343 = (t167 * t450 + t168 * t451 + (t201 * t451 + t336 * t473) * mrSges(6,3)) * pkin(4) - t373;
t12 = (t481 + t51 / 0.2e1) * mrSges(6,2) + (-t43 / 0.2e1 - t50 / 0.2e1) * mrSges(6,1) + t343 + t373;
t18 = (t478 + t520) * mrSges(6,2) + (t477 + t117 / 0.2e1) * mrSges(6,1);
t310 = (t334 * mrSges(6,1) + t336 * mrSges(6,2)) * pkin(4);
t345 = -qJD(1) * t12 - qJD(2) * t18 + qJD(4) * t310;
t303 = t310 * qJD(5);
t94 = m(6) * t386 + (t403 - t404) * t388;
t76 = m(6) * t387 + (t407 - t408) * t388;
t53 = t384 - t507 / 0.2e1;
t49 = t385 - t508 / 0.2e1;
t11 = -t441 / 0.2e1 - t442 / 0.2e1 + t440 / 0.2e1 - t439 / 0.2e1 + t343 - t373;
t9 = t339 + t329 * t485 + mrSges(4,2) * t397 / 0.2e1 + mrSges(4,1) * t399 / 0.2e1 - t344;
t4 = t340 + t368;
t1 = (t207 / 0.4e1 + t195 / 0.4e1) * t301 + (t383 + t251 / 0.4e1 + t248 / 0.4e1) * t281 + t486 + t524;
t13 = [-qJD(2) * t3 + qJD(3) * t16 + qJD(4) * t6 + qJD(5) * t7, t9 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) - t414 + ((-pkin(2) * t436 + Ifges(3,5) + (-pkin(2) * mrSges(4,1) + (Ifges(4,1) - Ifges(4,2)) * t333) * t332 + (-m(4) * pkin(2) - mrSges(4,1) * t333 + mrSges(4,2) * t332 - mrSges(3,1)) * pkin(6) + (t331 - t330) * Ifges(4,4)) * t337 + t359 * t499 + t360 * t455 + t361 * t464 + t362 * t453 + t136 * t471 + t139 * t469 + t245 * t133 + t272 * t365 + t249 * t459 + t251 * t456 + t309 * (-t301 * mrSges(5,1) - t302 * mrSges(5,2)) + t325 * t366 + 0.2e1 * (t117 * t352 + t272 * t245 + t353 * t488) * t483 + 0.2e1 * (t144 * t254 + t145 * t255 + t309 * t325) * t484 + (t117 * t199 - t202 * t488 - t241 * t353 + t352 * t369) * mrSges(6,3) + (t144 * t302 + t145 * t301 - t254 * t284 + t255 * t282) * mrSges(5,3) + (mrSges(5,1) * t254 + mrSges(6,1) * t488 + mrSges(3,2) * pkin(6) - mrSges(5,2) * t255 - mrSges(6,2) * t117 + Ifges(4,5) * t332 - Ifges(5,5) * t302 + Ifges(6,5) * t241 + Ifges(4,6) * t333 + Ifges(5,6) * t301 + Ifges(6,6) * t369 - qJ(3) * t367 - Ifges(3,6)) * t335 + t505 * (-t273 * t332 + t274 * t333)) * qJD(2), qJD(2) * t9 + qJD(4) * t76 + qJD(5) * t49 + t411, t413 + t1 * qJD(2) + t76 * qJD(3) + (-t143 * mrSges(5,1) - t142 * mrSges(5,2) - t439 + t440 + t490) * qJD(4) + t11 * qJD(5) + (m(6) * (t334 * t51 + t336 * t50) + (-t201 * t334 - t336 * t370) * mrSges(6,3)) * t430, t412 + t4 * qJD(2) + t49 * qJD(3) + t11 * qJD(4) + (t396 - t441 - t442) * qJD(5); qJD(3) * t10 + qJD(4) * t2 + qJD(5) * t5 + t414, qJD(3) * t21 + qJD(4) * t8 + qJD(5) * t14, qJD(4) * t94 + qJD(5) * t53 - t356, t94 * qJD(3) + (-t255 * mrSges(5,1) - t254 * mrSges(5,2) + t491 + t523) * qJD(4) + t532 + (m(6) * (-t117 * t336 + t334 * t488) + (-t241 * t334 - t336 * t369) * mrSges(6,3)) * t430 + t358, t53 * qJD(3) + t19 * qJD(4) + t357 + t532; -qJD(2) * t10 + qJD(4) * t27 + qJD(5) * t48 - t411, qJD(4) * t38 + qJD(5) * t52 + t356, 0, t355, t354; -qJD(2) * t2 - qJD(3) * t27 + qJD(5) * t12 - t413, -qJD(3) * t38 + qJD(5) * t18 - t358, -t355, -t303, -t303 - t345; -qJD(2) * t5 - qJD(3) * t48 - qJD(4) * t12 - t412, -qJD(3) * t52 - qJD(4) * t18 - t357, -t354, t345, 0;];
Cq = t13;
