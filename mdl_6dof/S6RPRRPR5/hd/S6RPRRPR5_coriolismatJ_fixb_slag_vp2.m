% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:27
% EndTime: 2019-03-09 05:12:41
% DurationCPUTime: 8.31s
% Computational Cost: add. (21103->519), mult. (40840->670), div. (0->0), fcn. (48555->8), ass. (0->299)
t290 = sin(qJ(6));
t281 = t290 * mrSges(7,1);
t293 = cos(qJ(6));
t282 = t293 * mrSges(7,2);
t506 = t282 + t281;
t288 = sin(pkin(10));
t457 = pkin(7) + qJ(2);
t261 = t457 * t288;
t289 = cos(pkin(10));
t262 = t457 * t289;
t292 = sin(qJ(3));
t294 = cos(qJ(3));
t241 = -t261 * t294 - t262 * t292;
t259 = t288 * t294 + t289 * t292;
t185 = -pkin(8) * t259 + t241;
t242 = -t261 * t292 + t262 * t294;
t258 = -t288 * t292 + t289 * t294;
t186 = pkin(8) * t258 + t242;
t291 = sin(qJ(4));
t463 = cos(qJ(4));
t107 = t463 * t185 - t291 * t186;
t314 = t291 * t258 + t259 * t463;
t70 = -pkin(5) * t314 + t107;
t553 = t70 * t506;
t557 = t553 / 0.2e1;
t507 = t185 * t291 + t463 * t186;
t528 = t507 * mrSges(6,2);
t529 = t507 * mrSges(5,1);
t543 = t107 * mrSges(6,3);
t544 = t107 * mrSges(5,2);
t556 = t528 - t529 + t543 - t544 + t553;
t555 = qJ(5) * t70;
t237 = -t463 * t258 + t259 * t291;
t521 = -t237 * pkin(5) + t507;
t554 = t521 * t70;
t517 = -t314 / 0.2e1;
t477 = t237 / 0.2e1;
t535 = Ifges(6,6) * t477;
t518 = Ifges(6,2) * t517 + t535;
t426 = t293 * mrSges(7,1);
t431 = t290 * mrSges(7,2);
t263 = -t426 + t431;
t145 = t263 * t237;
t511 = t263 * t314;
t552 = t70 * t145 + t521 * t511;
t361 = -pkin(2) * t289 - pkin(1);
t244 = -pkin(3) * t258 + t361;
t309 = -qJ(5) * t314 + t244;
t104 = pkin(4) * t237 + t309;
t157 = -mrSges(6,2) * t237 - mrSges(6,3) * t314;
t551 = m(6) * t104 + t157;
t424 = t293 * Ifges(7,4);
t336 = t290 * Ifges(7,1) + t424;
t520 = -Ifges(7,5) * t237 + t314 * t336;
t421 = t293 * t520;
t265 = Ifges(7,5) * t293 - Ifges(7,6) * t290;
t524 = t237 * t265;
t549 = -t524 / 0.4e1 + t528 / 0.2e1 - t529 / 0.2e1 + t421 / 0.4e1 + t543 / 0.2e1 - t544 / 0.2e1;
t547 = t520 / 0.2e1;
t542 = t290 * t521;
t541 = t293 * t521;
t286 = t290 ^ 2;
t287 = t293 ^ 2;
t383 = t286 + t287;
t484 = pkin(4) + pkin(9);
t525 = qJ(5) * t237;
t91 = t314 * t484 + t525;
t156 = pkin(4) * t314 + t525;
t539 = t107 * t291 - t463 * t507;
t495 = m(6) / 0.2e1;
t537 = (-pkin(4) * t507 + qJ(5) * t107) * t495;
t494 = m(6) / 0.4e1;
t516 = m(7) * t383;
t536 = 0.2e1 * (t516 / 0.4e1 + t494) * t314;
t534 = -t237 / 0.2e1;
t462 = m(6) * t507;
t532 = Ifges(6,3) * t477;
t531 = t237 * mrSges(7,1);
t530 = t237 * mrSges(7,2);
t526 = mrSges(7,3) * t383 + mrSges(5,1);
t375 = mrSges(6,3) + t506;
t459 = t259 * pkin(3);
t72 = t459 + t91;
t46 = -t290 * t72 + t541;
t423 = t293 * t46;
t47 = t293 * t72 + t542;
t429 = t290 * t47;
t522 = -t429 - t423;
t449 = Ifges(7,2) * t293;
t454 = Ifges(7,4) * t290;
t335 = t449 + t454;
t519 = -Ifges(7,6) * t237 + t314 * t335;
t458 = mrSges(6,1) + mrSges(5,3);
t512 = t314 * Ifges(6,6);
t349 = t286 / 0.2e1 + t287 / 0.2e1;
t339 = t349 * mrSges(7,3);
t450 = Ifges(7,2) * t290;
t267 = t424 - t450;
t386 = t293 * t267;
t456 = Ifges(7,1) * t293;
t269 = -t454 + t456;
t393 = t290 * t269;
t505 = t386 / 0.4e1 + t393 / 0.4e1;
t504 = t386 / 0.2e1 + t393 / 0.2e1;
t489 = -mrSges(7,3) / 0.2e1;
t503 = t383 * t489;
t447 = Ifges(7,6) * t293;
t451 = Ifges(7,5) * t290;
t502 = t447 / 0.2e1 + t451 / 0.2e1;
t405 = t237 * t290;
t151 = mrSges(7,1) * t314 - mrSges(7,3) * t405;
t390 = t293 * t151;
t404 = t237 * t293;
t154 = -mrSges(7,2) * t314 + mrSges(7,3) * t404;
t396 = t290 * t154;
t501 = t390 + t396;
t500 = m(7) / 0.4e1 + t494;
t445 = qJ(5) * mrSges(6,1);
t499 = -t445 / 0.2e1 + t505;
t250 = -(-m(7) - m(6)) * qJ(5) + t375;
t498 = 0.2e1 * m(7);
t497 = t259 ^ 2;
t496 = 2 * qJD(4);
t493 = m(7) / 0.2e1;
t491 = m(5) * pkin(3);
t490 = mrSges(7,2) / 0.2e1;
t488 = -t46 / 0.2e1;
t487 = t521 / 0.4e1;
t433 = t314 * Ifges(7,6);
t86 = t237 * t335 + t433;
t486 = -t86 / 0.4e1;
t483 = pkin(4) * mrSges(6,1);
t482 = qJ(5) / 0.2e1;
t481 = -t151 / 0.2e1;
t480 = -t512 / 0.2e1 + t532;
t474 = -t314 / 0.4e1;
t473 = t314 / 0.2e1;
t472 = -t265 / 0.2e1;
t379 = t463 * pkin(3);
t277 = -t379 - pkin(4);
t274 = -pkin(9) + t277;
t471 = t274 / 0.2e1;
t460 = pkin(3) * t291;
t276 = qJ(5) + t460;
t470 = t276 / 0.2e1;
t469 = -t290 / 0.2e1;
t468 = -t290 / 0.4e1;
t467 = t290 / 0.2e1;
t466 = -t293 / 0.2e1;
t465 = t293 / 0.2e1;
t464 = -t484 / 0.2e1;
t461 = m(6) * t156;
t455 = Ifges(6,4) * t237;
t453 = Ifges(5,5) * t237;
t452 = Ifges(6,5) * t314;
t448 = Ifges(5,6) * t314;
t446 = Ifges(7,3) * t237;
t110 = t156 + t459;
t394 = t290 * t314;
t150 = -mrSges(7,3) * t394 - t531;
t387 = t293 * t314;
t153 = mrSges(7,3) * t387 + t530;
t160 = Ifges(5,4) * t314 - t237 * Ifges(5,2);
t235 = Ifges(5,4) * t237;
t161 = Ifges(5,1) * t314 - t235;
t334 = t447 + t451;
t233 = t237 * mrSges(6,3);
t343 = -mrSges(6,2) * t314 + t233;
t232 = t237 * mrSges(5,2);
t344 = mrSges(5,1) * t314 - t232;
t345 = t259 * mrSges(4,1) + t258 * mrSges(4,2);
t355 = t387 / 0.2e1;
t358 = t394 / 0.2e1;
t66 = t237 * t484 + t309;
t42 = -t290 * t66 - t293 * t70;
t43 = -t290 * t70 + t293 * t66;
t84 = Ifges(7,3) * t314 + t237 * t334;
t210 = Ifges(7,4) * t404;
t434 = t314 * Ifges(7,5);
t89 = Ifges(7,1) * t405 + t210 + t434;
t1 = t86 * t355 + t89 * t358 - t446 * t473 + t552 + t405 * t547 + (Ifges(4,4) * t258 + (Ifges(4,1) - Ifges(4,2)) * t259) * t258 + (m(5) * t459 + t344) * t244 - t497 * Ifges(4,4) + t42 * t150 + t46 * t151 + t43 * t153 + t47 * t154 + (mrSges(5,1) * t459 - Ifges(5,1) * t473 - Ifges(5,4) * t534 + 0.2e1 * t518) * t237 + t551 * t110 + (t161 + t84) * t534 + ((-Ifges(5,4) + t334) * t473 + t532 - Ifges(5,2) * t534 + t480 + mrSges(5,2) * t459) * t314 + m(7) * (t42 * t46 + t43 * t47 + t554) + t519 * t404 / 0.2e1 + t361 * t345 + t104 * t343 + (t512 + t160) * t517;
t444 = t1 * qJD(1);
t435 = t237 * mrSges(6,1);
t432 = t276 * mrSges(6,1);
t430 = t290 * mrSges(7,3);
t428 = t290 * t519;
t427 = t290 * t89;
t425 = t293 * mrSges(7,3);
t422 = t293 * t86;
t152 = -t314 * t430 - t531;
t155 = t314 * t425 + t530;
t52 = -t290 * t91 + t541;
t53 = t293 * t91 + t542;
t4 = t52 * t151 + t42 * t152 + t53 * t154 + t43 * t155 + t156 * t157 + t104 * t461 + m(7) * (t42 * t52 + t43 * t53 + t554) + (-t244 * mrSges(5,2) + t535 + t104 * mrSges(6,3) + t235 / 0.2e1 + t518 - t161 / 0.2e1 - t84 / 0.2e1 + t519 * t465 + t520 * t467) * t237 + (t244 * mrSges(5,1) - t104 * mrSges(6,2) + t480 - t160 / 0.2e1 + t422 / 0.2e1 + t427 / 0.2e1 + (-Ifges(6,6) / 0.2e1 - Ifges(5,4) / 0.2e1 + t502) * t314 + (-Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t237) * t314 + t552;
t420 = t4 * qJD(1);
t147 = t506 * t237;
t254 = -t393 / 0.2e1;
t332 = t290 * t43 + t293 * t42;
t5 = -t42 * t154 + t43 * t151 - t521 * t147 + (t86 * t467 + t314 * t472 + (t449 * t467 + t254) * t237 + t332 * mrSges(7,3) + (t210 + t89) * t466) * t237;
t419 = t5 * qJD(1);
t418 = t52 * t293;
t417 = t53 * t290;
t413 = qJ(5) * t263;
t10 = (t258 ^ 2 + t497) * mrSges(4,3) - (-t237 * t458 + t145) * t237 + (t314 * t458 + t501) * t314 + m(7) * (-t237 * t521 + t314 * t332) + m(4) * (-t241 * t259 + t242 * t258) + (m(3) * qJ(2) + mrSges(3,3)) * (t288 ^ 2 + t289 ^ 2) + (m(6) + m(5)) * (-t107 * t314 - t237 * t507);
t412 = qJD(1) * t10;
t388 = t293 * t154;
t398 = t290 * t151;
t18 = (m(7) * (t290 * t42 - t293 * t43) - t388 + t398 - t551) * t314;
t411 = qJD(1) * t18;
t342 = t383 * t274;
t359 = t506 * t534;
t381 = t491 / 0.2e1;
t297 = t359 + (t277 * t495 + t342 * t493 - t381 * t463 + t503) * t314 + (-t291 * t381 - (t493 + t495) * t276) * t237;
t301 = t110 * t495 + (-t290 * t46 + t293 * t47) * t493 + t150 * t469 + t153 * t465 + t259 * t381;
t11 = t297 - t301 - t343 - t344 - t345;
t410 = t11 * qJD(1);
t341 = t383 * t484;
t377 = -mrSges(5,1) / 0.2e1 + mrSges(6,2) / 0.2e1;
t299 = (-t339 + t377) * t314 + t232 / 0.2e1 - t233 / 0.2e1 - t156 * t495 + (-t314 * t341 - t525) * t493 + t359;
t303 = -t461 / 0.2e1 - m(7) * (-t290 * t52 + t293 * t53) / 0.2e1 + t152 * t467 + t155 * t466;
t376 = mrSges(6,3) / 0.2e1 - mrSges(5,2) / 0.2e1;
t12 = -t237 * t376 + t314 * t377 + t299 + t303;
t409 = t12 * qJD(1);
t312 = (t281 / 0.2e1 + t282 / 0.2e1) * t314;
t315 = t388 / 0.2e1 - t398 / 0.2e1;
t22 = t237 * t339 + t312 - t315;
t407 = t22 * qJD(1);
t322 = t431 / 0.2e1 - t426 / 0.2e1;
t311 = t322 * t314;
t317 = -t396 / 0.2e1 - t390 / 0.2e1;
t24 = -t311 - t317;
t400 = t24 * qJD(1);
t399 = t276 * t263;
t397 = t290 * t153;
t395 = t290 * t155;
t391 = t293 * t150;
t389 = t293 * t152;
t54 = 0.2e1 * t536;
t385 = t54 * qJD(1);
t380 = t521 * t493;
t378 = t460 / 0.2e1;
t374 = t463 * mrSges(5,2);
t369 = t435 / 0.2e1;
t366 = -t430 / 0.2e1;
t365 = -t425 / 0.2e1;
t364 = -t521 * t263 / 0.2e1;
t363 = -t210 / 0.4e1 - t89 / 0.4e1;
t362 = mrSges(6,2) * t460 + t375 * t379;
t360 = t463 * t276;
t352 = t150 * t464;
t351 = t153 * t464;
t350 = t474 + t517;
t340 = -t335 * t469 + t254 - t386 / 0.2e1 - t336 * t465;
t338 = t154 * t464 + t486;
t330 = t417 + t418;
t296 = t52 * t365 + t53 * t366 + Ifges(6,5) * t473 + Ifges(6,4) * t477 + (t395 + t389) * t471 + t501 * t378 + (mrSges(6,1) * t378 + t505) * t314 + t145 * t379 / 0.2e1 + (-t539 * pkin(3) + t276 * t107 + t277 * t507) * t495 + (Ifges(5,5) + (t277 + t379) * mrSges(6,1)) * t534 + (t70 * t276 + t330 * t274 + (t291 * t332 + t463 * t521) * pkin(3)) * t493 + t519 * t468 + (t432 + Ifges(5,6)) * t517 + t511 * t470 + t557 + t549;
t300 = t537 + (t484 * t522 + t555) * t493 + t511 * t482 + t557;
t3 = -t296 - t428 / 0.4e1 - t453 / 0.2e1 + t455 / 0.2e1 + t452 / 0.2e1 - t448 / 0.2e1 + t300 + t46 * t365 + t293 * t352 + t290 * t351 + t47 * t366 + pkin(4) * t369 + t499 * t314 + t549;
t74 = -t362 + (t374 - m(7) * (t291 * t342 + t360) - m(6) * (t277 * t291 + t360)) * pkin(3) + t526 * t460;
t329 = -t3 * qJD(1) - t74 * qJD(3);
t162 = t340 - t399;
t318 = -t267 / 0.4e1 - t336 / 0.4e1 + t450 / 0.4e1;
t319 = -t335 / 0.4e1 + t269 / 0.4e1 + t456 / 0.4e1;
t302 = t319 * t293 + (-t424 / 0.4e1 + t318) * t290;
t298 = (-t274 * t339 + t302) * t237 + t147 * t470 + t364;
t307 = t446 / 0.2e1 + mrSges(7,1) * t488 + t47 * t490;
t7 = (Ifges(7,6) * t350 + t154 * t471 + t486) * t293 + (Ifges(7,5) * t350 + t274 * t481 + t363) * t290 + t298 + t307;
t328 = t7 * qJD(1) + t162 * qJD(3);
t305 = t369 - t397 / 0.2e1 - t391 / 0.2e1;
t306 = (-mrSges(6,1) / 0.2e1 + t322) * t237;
t16 = t306 + (t487 - t429 / 0.4e1 - t423 / 0.4e1) * t498 + t305;
t240 = 0.4e1 * t276 * t500 + t375;
t325 = -qJD(1) * t16 - qJD(3) * t240;
t324 = t484 * t339;
t323 = -t52 * mrSges(7,1) / 0.2e1 + t53 * t490;
t321 = t147 * t482 + t364;
t320 = -t481 * t484 + t363;
t316 = -t389 / 0.2e1 - t395 / 0.2e1;
t163 = t340 - t413;
t9 = (Ifges(7,3) / 0.2e1 + t324) * t237 + (-0.3e1 / 0.4e1 * t434 + t318 * t237 + t320) * t290 + (-0.3e1 / 0.4e1 * t433 + (-t454 / 0.4e1 + t319) * t237 + t338) * t293 + t321 + t323;
t92 = (t482 + t470) * t263 + (mrSges(7,1) * t378 + t336 / 0.2e1 + t267 / 0.2e1) * t293 + (-mrSges(7,2) * t460 / 0.2e1 + t269 / 0.2e1 - t335 / 0.2e1) * t290;
t310 = t9 * qJD(1) - t92 * qJD(3) + t163 * qJD(4);
t164 = (-0.1e1 / 0.2e1 + t349) * m(7) * t460 - t250;
t21 = t322 * t237 + (t487 - t417 / 0.4e1 - t418 / 0.4e1) * t498 + t316;
t308 = qJD(1) * t21 - qJD(3) * t164 + qJD(4) * t250;
t278 = qJ(5) * t379;
t165 = (t516 / 0.2e1 + t495) * t460 + t375 + t500 * (0.4e1 * qJ(5) + 0.2e1 * t460);
t93 = -t413 / 0.2e1 - t399 / 0.2e1 - t322 * t460 + t340;
t55 = m(6) * t517 - t473 * t516 + t536;
t25 = -t311 + t317;
t23 = t237 * t503 + t312 + t315;
t20 = t330 * t493 + t380 + t462 + (-mrSges(6,1) + t322) * t237 - t316;
t17 = -t522 * t493 + t507 * t495 + t380 + t462 / 0.2e1 + t306 - t305;
t15 = t297 + t301;
t13 = mrSges(5,1) * t473 + mrSges(5,2) * t534 + mrSges(6,2) * t517 + mrSges(6,3) * t477 + t299 - t303;
t8 = Ifges(7,3) * t534 + (-t433 / 0.4e1 + t338) * t293 + (-t434 / 0.4e1 + t320) * t290 + (t324 + t302) * t237 + t321 - t323 + t502 * t314;
t6 = -t427 / 0.4e1 - t422 / 0.4e1 + t210 * t468 + t334 * t474 + Ifges(7,6) * t355 + Ifges(7,5) * t358 + t315 * t274 + t298 - t307;
t2 = t296 + t377 * t507 - (t265 / 0.4e1 + Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1 - t483 / 0.2e1) * t237 + (Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1 + t499) * t314 + (-t519 / 0.4e1 + t351 + t47 * t489) * t290 + (t520 / 0.4e1 + mrSges(7,3) * t488 + t352) * t293 + t376 * t107 + t300;
t14 = [qJD(2) * t10 + qJD(3) * t1 + qJD(4) * t4 + qJD(5) * t18 - qJD(6) * t5, qJD(3) * t15 + qJD(4) * t13 + qJD(5) * t55 + qJD(6) * t25 + t412, t15 * qJD(2) + t2 * qJD(4) + t17 * qJD(5) + t6 * qJD(6) + t444 + (-t428 / 0.2e1 + t455 - t524 / 0.2e1 + t539 * t491 + t237 * mrSges(5,3) * t379 + (-mrSges(5,3) * t460 - t432 + t504) * t314 - t448 - Ifges(4,6) * t259 + Ifges(4,5) * t258 - t241 * mrSges(4,2) - t242 * mrSges(4,1) + t452 - t453 + (m(6) * t107 + m(7) * t70 + t511) * t276 + t522 * mrSges(7,3) + t421 / 0.2e1 + (-t435 + t462) * t277 + (-m(7) * t522 + t391 + t397) * t274 + t556) * qJD(3), t420 + t13 * qJD(2) + t2 * qJD(3) + t20 * qJD(5) + t8 * qJD(6) + ((-t330 * t484 + t555) * t493 + t537) * t496 + (qJ(5) * t511 + (-t52 * mrSges(7,3) - t484 * t152 + t547) * t293 + (-t53 * mrSges(7,3) - t484 * t155 - t519 / 0.2e1) * t290 + (Ifges(6,4) - Ifges(5,5) + t472 + t483) * t237 + (Ifges(6,5) - Ifges(5,6) - t445 + t504) * t314 + t556) * qJD(4), qJD(2) * t55 + qJD(3) * t17 + qJD(4) * t20 + qJD(6) * t23 + t411, -t419 + t25 * qJD(2) + t6 * qJD(3) + t8 * qJD(4) + t23 * qJD(5) + (-mrSges(7,1) * t43 - mrSges(7,2) * t42 + t524) * qJD(6); -qJD(3) * t11 - qJD(4) * t12 - qJD(5) * t54 - qJD(6) * t24 - t412, 0, -t410, -t409, -t385, qJD(6) * t263 - t400; qJD(2) * t11 - qJD(4) * t3 + qJD(5) * t16 + qJD(6) * t7 - t444, t410, -qJD(4) * t74 + qJD(5) * t240 + qJD(6) * t162, t165 * qJD(5) + t93 * qJD(6) + ((-t341 * t460 + t278) * t493 + (-pkin(4) * t460 + t278) * t495) * t496 + t329 + (t362 + (-t291 * t526 - t374) * pkin(3)) * qJD(4), qJD(4) * t165 - t325, t93 * qJD(4) + (-t274 * t506 - t334) * qJD(6) + t328; qJD(2) * t12 + qJD(3) * t3 + qJD(5) * t21 + qJD(6) * t9 - t420, t409, -qJD(5) * t164 - qJD(6) * t92 - t329, qJD(5) * t250 + qJD(6) * t163, t308 ((mrSges(7,2) * t484 - Ifges(7,6)) * t293 + (mrSges(7,1) * t484 - Ifges(7,5)) * t290) * qJD(6) + t310; qJD(2) * t54 - qJD(3) * t16 - qJD(4) * t21 - qJD(6) * t22 - t411, t385, qJD(4) * t164 + t325, -t308, 0, -qJD(6) * t506 - t407; qJD(2) * t24 - qJD(3) * t7 - qJD(4) * t9 + qJD(5) * t22 + t419, t400, qJD(4) * t92 - t328, -t310, t407, 0;];
Cq  = t14;
