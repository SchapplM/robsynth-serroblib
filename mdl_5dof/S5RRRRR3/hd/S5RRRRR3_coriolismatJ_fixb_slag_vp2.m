% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:11
% EndTime: 2019-12-05 18:55:29
% DurationCPUTime: 8.30s
% Computational Cost: add. (14946->502), mult. (33604->675), div. (0->0), fcn. (37325->8), ass. (0->309)
t312 = sin(qJ(5));
t314 = cos(qJ(4));
t313 = sin(qJ(4));
t517 = cos(qJ(5));
t406 = t517 * t313;
t352 = t312 * t314 + t406;
t451 = t312 * t313;
t576 = t517 * t314 - t451;
t222 = -mrSges(6,1) * t576 + mrSges(6,2) * t352;
t498 = Ifges(5,4) * t313;
t292 = Ifges(5,2) * t314 + t498;
t308 = Ifges(5,4) * t314;
t294 = Ifges(5,1) * t313 + t308;
t509 = pkin(3) * t313;
t604 = -t314 / 0.2e1;
t605 = t313 / 0.2e1;
t614 = -t222 * t509 + t292 * t605 + t294 * t604;
t295 = Ifges(5,1) * t314 - t498;
t365 = Ifges(5,2) * t313 - t308;
t562 = t295 * t605 + t365 * t604 - t614;
t483 = t352 * mrSges(6,3);
t515 = sin(qJ(3));
t427 = t515 * pkin(1);
t304 = t427 + pkin(5);
t241 = t576 * t304;
t267 = t576 * pkin(5);
t603 = t267 + t241;
t613 = t603 * t483 / 0.2e1;
t518 = cos(qJ(3));
t429 = t518 * pkin(1);
t305 = -t429 - pkin(2);
t508 = t314 * pkin(3);
t289 = t305 - t508;
t306 = -pkin(2) - t508;
t424 = t509 / 0.2e1;
t501 = mrSges(5,1) * t313;
t291 = mrSges(5,2) * t314 + t501;
t454 = t305 * t291;
t513 = pkin(2) * t291;
t612 = -t513 / 0.2e1 + t454 / 0.2e1 + m(6) * (t289 + t306) * t424 + t562;
t516 = sin(qJ(2));
t519 = cos(qJ(2));
t283 = t515 * t516 - t518 * t519;
t284 = -t515 * t519 - t518 * t516;
t230 = -t284 * pkin(2) + t283 * pkin(5);
t555 = m(6) * pkin(3);
t436 = t555 / 0.2e1;
t479 = t314 * mrSges(5,1);
t480 = t313 * mrSges(5,2);
t575 = -t480 / 0.2e1 + t479 / 0.2e1;
t279 = t284 * pkin(3);
t153 = t230 * t314 - t279;
t89 = t517 * t153 - t230 * t451;
t606 = t89 / 0.2e1;
t608 = -mrSges(6,2) / 0.2e1;
t90 = t312 * t153 + t230 * t406;
t585 = mrSges(6,1) * t606 + t90 * t608;
t611 = -(t312 * t90 + t517 * t89) * t436 - t575 * t230 - t585;
t428 = t516 * pkin(1);
t213 = t428 + t230;
t143 = t213 * t314 - t279;
t74 = t312 * t143 + t213 * t406;
t607 = -t74 / 0.2e1;
t609 = mrSges(6,1) / 0.2e1;
t72 = t517 * t143 - t213 * t451;
t586 = mrSges(6,2) * t607 + t72 * t609;
t610 = -(t312 * t74 + t517 * t72) * t436 - t575 * t213 - t586;
t559 = -m(6) / 0.2e1;
t527 = -t284 / 0.2e1;
t310 = t313 ^ 2;
t311 = t314 ^ 2;
t580 = t310 + t311;
t195 = t352 * t283;
t134 = mrSges(6,2) * t284 + mrSges(6,3) * t195;
t197 = t576 * t283;
t136 = -mrSges(6,1) * t284 + mrSges(6,3) * t197;
t431 = pkin(3) * t517;
t510 = pkin(3) * t312;
t542 = -t197 / 0.2e1;
t544 = t195 / 0.2e1;
t355 = Ifges(6,5) * t542 + Ifges(6,6) * t544;
t571 = Ifges(6,3) * t527 + t355;
t602 = -t136 * t431 / 0.2e1 - t134 * t510 / 0.2e1 - Ifges(5,3) * t527 - t571;
t430 = t519 * pkin(1);
t214 = t283 * pkin(2) + t284 * pkin(5) - t430;
t105 = t352 * t214;
t472 = t214 * t314;
t144 = pkin(3) * t283 + t472;
t75 = t312 * t144 + t214 * t406;
t599 = t105 - t75;
t106 = t576 * t214;
t73 = t517 * t144 - t214 * t451;
t598 = -t106 + t73;
t273 = Ifges(6,4) * t576;
t223 = -Ifges(6,2) * t352 + t273;
t496 = Ifges(6,4) * t352;
t224 = Ifges(6,2) * t576 + t496;
t225 = Ifges(6,1) * t576 - t496;
t226 = Ifges(6,1) * t352 + t273;
t372 = (t223 + t226) * t576 / 0.2e1 + (-t224 / 0.2e1 + t225 / 0.2e1) * t352;
t221 = mrSges(6,1) * t352 + mrSges(6,2) * t576;
t452 = t306 * t221;
t456 = t289 * t221;
t596 = t372 + t452 / 0.2e1 + t456 / 0.2e1 + t613;
t194 = t576 * t284;
t196 = t352 * t284;
t384 = t226 / 0.4e1 + t223 / 0.4e1;
t385 = -t225 / 0.4e1 + t224 / 0.4e1;
t594 = t385 * t194 + t384 * t196;
t438 = Ifges(6,5) * t576 - Ifges(6,6) * t352;
t268 = t352 * pkin(5);
t581 = -t267 * mrSges(6,1) + t268 * mrSges(6,2);
t357 = t438 + t581;
t593 = qJD(5) * t357;
t242 = t352 * t304;
t582 = -t241 * mrSges(6,1) + t242 * mrSges(6,2);
t358 = t438 + t582;
t592 = qJD(5) * t358;
t461 = t283 * t313;
t215 = mrSges(5,2) * t284 + mrSges(5,3) * t461;
t446 = t314 * t215;
t217 = mrSges(5,3) * t283 * t314 - mrSges(5,1) * t284;
t449 = t313 * t217;
t353 = t446 / 0.2e1 - t449 / 0.2e1;
t433 = t306 * t509;
t378 = t283 * t433;
t591 = (t267 * t74 - t268 * t72 - t378) * t559 - t353 * pkin(5);
t589 = t283 / 0.2e1;
t588 = Ifges(5,6) * t527;
t425 = -pkin(3) * t221 / 0.2e1;
t375 = t313 * t425;
t579 = (t375 + Ifges(6,3) / 0.2e1) * t284 - t355;
t530 = t268 / 0.2e1;
t532 = t267 / 0.2e1;
t578 = t194 * t532 + t196 * t530;
t534 = t242 / 0.2e1;
t535 = t241 / 0.2e1;
t577 = t194 * t535 + t196 * t534;
t383 = t292 / 0.4e1 - t295 / 0.4e1;
t481 = t283 * Ifges(5,6);
t146 = t365 * t284 + t481;
t210 = t284 * t294;
t572 = -t146 / 0.4e1 + t210 / 0.4e1;
t100 = -mrSges(6,1) * t194 + mrSges(6,2) * t196;
t137 = mrSges(6,1) * t283 + t194 * mrSges(6,3);
t531 = -t267 / 0.2e1;
t489 = t196 * mrSges(6,3);
t135 = -mrSges(6,2) * t283 + t489;
t549 = t135 / 0.2e1;
t570 = t306 * t100 / 0.2e1 + t137 * t531 - t268 * t549;
t526 = t289 / 0.2e1;
t533 = -t242 / 0.2e1;
t548 = t137 / 0.2e1;
t569 = t100 * t526 + t135 * t533 - t241 * t548;
t370 = t518 * t517;
t407 = t314 * t518;
t249 = (-t312 * t407 - t313 * t370) * pkin(1);
t408 = t313 * t518;
t250 = (-t312 * t408 + t314 * t370) * pkin(1);
t439 = t249 * t609 + t250 * t608;
t366 = t284 * t375;
t568 = t366 + t571;
t459 = t284 * t314;
t218 = mrSges(5,1) * t283 + mrSges(5,3) * t459;
t445 = t314 * t218;
t460 = t284 * t313;
t216 = -mrSges(5,2) * t283 + mrSges(5,3) * t460;
t450 = t313 * t216;
t500 = mrSges(5,3) * t284;
t566 = -t450 / 0.2e1 - t445 / 0.2e1 + t580 * t500 / 0.2e1;
t290 = -t479 + t480;
t207 = t290 * t284;
t565 = -pkin(2) * t207 / 0.2e1 + t570;
t522 = t305 / 0.2e1;
t564 = t207 * t522 + t569;
t497 = Ifges(6,4) * t194;
t104 = Ifges(6,1) * t196 + t497;
t528 = t283 / 0.4e1;
t80 = Ifges(6,2) * t196 + t283 * Ifges(6,6) - t497;
t561 = t438 * t528 - (-t104 / 0.4e1 + t80 / 0.4e1) * t352;
t558 = m(6) / 0.2e1;
t557 = -pkin(5) / 0.2e1;
t556 = m(5) * pkin(1);
t554 = -t73 / 0.2e1;
t553 = Ifges(6,4) * t542 + Ifges(6,2) * t544 + Ifges(6,6) * t527;
t551 = Ifges(6,1) * t197 / 0.2e1 - Ifges(6,4) * t195 / 0.2e1 + t284 * Ifges(6,5) / 0.2e1;
t547 = t365 * t589 + t588;
t545 = t194 / 0.2e1;
t543 = -t196 / 0.2e1;
t540 = -t222 / 0.2e1;
t523 = -t304 / 0.2e1;
t521 = -t306 / 0.2e1;
t208 = t291 * t283;
t514 = pkin(2) * t208;
t512 = pkin(3) * t100;
t102 = -mrSges(6,1) * t196 - mrSges(6,2) * t194;
t511 = pkin(3) * t102;
t506 = t73 * mrSges(6,2);
t504 = t75 * mrSges(6,1);
t307 = Ifges(5,5) * t314;
t491 = t105 * mrSges(6,1);
t490 = t106 * mrSges(6,2);
t484 = t576 * mrSges(6,3);
t482 = t283 * Ifges(5,5);
t458 = t284 * pkin(3) ^ 2;
t263 = t310 * t458;
t227 = t283 * t263;
t101 = -mrSges(6,1) * t195 - mrSges(6,2) * t197;
t147 = -Ifges(5,5) * t284 - t295 * t283;
t418 = t307 / 0.2e1;
t321 = Ifges(6,5) * t545 + Ifges(6,6) * t543 + t147 * t604 + (t418 - Ifges(4,4)) * t284 + (-pkin(3) * t101 + t547 + t588) * t313 + (-Ifges(5,3) - Ifges(4,2) + Ifges(4,1) - Ifges(6,3)) * t283;
t413 = t481 / 0.2e1;
t349 = -t511 + t146 / 0.2e1 + t413;
t148 = -t284 * t295 + t482;
t414 = -t482 / 0.2e1;
t368 = -t148 / 0.2e1 + t414;
t323 = Ifges(4,4) * t283 + t349 * t313 + t368 * t314 + t355;
t188 = Ifges(6,4) * t196;
t82 = -Ifges(6,1) * t194 + t283 * Ifges(6,5) + t188;
t326 = t75 * t134 + t73 * t136 + t194 * t551 + t196 * t553 + t82 * t542 + t80 * t544 - (-mrSges(4,1) * t284 - mrSges(4,2) * t283) * t430;
t359 = t445 + t450;
t360 = t313 * t215 + t314 * t217;
t382 = m(5) * t580;
t5 = (t214 * t382 + t359) * t213 + t360 * t214 + t74 * t135 + t72 * t137 + (mrSges(4,1) * t428 + t323) * t283 + (-mrSges(4,2) * t428 + t321) * t284 + t326 + m(6) * (t72 * t73 + t74 * t75 + t227) + (-t516 ^ 2 + t519 ^ 2) * Ifges(3,4) + (-m(4) * pkin(1) ^ 2 + Ifges(3,1) - Ifges(3,2)) * t519 * t516;
t478 = t5 * qJD(1);
t7 = t90 * t135 + t89 * t137 + m(6) * (t73 * t89 + t75 * t90 + t227) + t359 * t230 + (t230 * t382 + t360) * t214 + t323 * t283 + t321 * t284 + t326;
t477 = t7 * qJD(1);
t476 = t89 * t352;
t475 = t90 * t576;
t209 = t284 * t292;
t103 = Ifges(6,2) * t194 + t188;
t441 = Ifges(6,5) * t196 + Ifges(6,6) * t194;
t324 = (t103 / 0.2e1 + t82 / 0.2e1) * t196 + (t80 / 0.2e1 + t75 * mrSges(6,3) - t104 / 0.2e1) * t194 - t73 * t489 + t441 * t589;
t11 = t106 * t135 - t105 * t137 + m(6) * (-t105 * t73 + t106 * t75) + (t216 * t314 - t218 * t313) * t214 + ((-t210 / 0.2e1 + t349) * t314 + (-t512 + t209 / 0.2e1 + m(6) * t314 * t458 - t368) * t313) * t284 + t324;
t474 = t11 * qJD(1);
t14 = t73 * t135 - t75 * t137 - t460 * t512 + t324;
t473 = t14 * qJD(1);
t471 = t241 * t134;
t470 = t242 * t136;
t468 = t267 * t134;
t467 = t268 * t136;
t457 = t289 * t101;
t455 = t305 * t208;
t453 = t306 * t101;
t437 = mrSges(6,3) * t510;
t435 = pkin(3) * t459;
t434 = t289 * t509;
t417 = t484 / 0.2e1;
t415 = -t483 / 0.2e1;
t410 = t554 + t106 / 0.2e1;
t409 = t75 / 0.2e1 - t105 / 0.2e1;
t386 = t148 / 0.4e1 + t209 / 0.4e1;
t381 = -Ifges(5,6) * t313 + t307;
t379 = mrSges(6,3) * t431;
t377 = t515 * t509;
t369 = -t408 / 0.2e1;
t316 = t561 + t383 * t459 + t572 * t313 + t106 * t417 + t283 * t418 + t102 * t424 + (t294 - t365) * t460 / 0.4e1 + (t209 + t148) * t314 / 0.4e1 + (t103 + t82) * t576 / 0.4e1 - t599 * t415 - Ifges(5,6) * t461 / 0.2e1 + t484 * t554 + t381 * t528 + t435 * t540 + t366 + t594 + t602;
t345 = -t598 * t241 + t599 * t242 - t263;
t2 = (-t289 * t435 + t345) * t558 + t316 + t577 * mrSges(6,3) + t566 * t304 + t564 + t610;
t337 = t372 + t456;
t21 = m(6) * t434 + t337 + t454 + t562;
t364 = t2 * qJD(1) + t21 * qJD(2);
t320 = -t249 * t483 + t250 * t484 + (-mrSges(4,1) + t222 + t290) * t427 + (t580 * mrSges(5,3) - mrSges(4,2)) * t429;
t348 = t580 * t518;
t30 = m(6) * (t241 * t250 - t242 * t249 + t289 * t427) + (t348 * t304 + t515 * t305) * t556 + t320;
t331 = t515 * t102 / 0.2e1 + t218 * t369 + t216 * t407 / 0.2e1;
t332 = t249 * t548 + t250 * t549 + t353 * t304;
t239 = t283 * t434;
t338 = t241 * t90 - t242 * t89 + t249 * t73 + t250 * t75 - t239;
t343 = -t515 * t291 / 0.2e1;
t8 = -(pkin(2) / 0.2e1 + t522) * t208 + (t530 + t533) * t136 + (t531 + t535) * t134 + (t521 + t526) * t101 + (t284 * t343 + t331) * pkin(1) + (-pkin(1) * t284 * t377 + t338) * t558 + (-(t606 - t72 / 0.2e1) * t352 + (t90 / 0.2e1 + t607) * t576) * mrSges(6,3) + t332 + t591;
t363 = t8 * qJD(1) + t30 * qJD(2);
t330 = (t82 / 0.4e1 + t103 / 0.4e1) * t576 + t561;
t318 = (mrSges(6,3) * t534 + t384) * t196 + (mrSges(6,3) * t535 + t385) * t194 + t330 + t569;
t10 = t318 - t586 + t579;
t362 = t10 * qJD(1) + qJD(2) * t337;
t361 = t446 - t449;
t356 = (t310 / 0.2e1 + t311 / 0.2e1) * t500;
t322 = (t517 * t249 + t250 * t312) * t436 + t439 + (mrSges(5,1) * t369 - mrSges(5,2) * t407 / 0.2e1) * pkin(1);
t16 = -t415 * t603 + t322 - t596 - t612;
t336 = t372 + t452;
t27 = m(6) * t433 + t336 - t513 + t562;
t344 = -t598 * t267 + t599 * t268 - t263;
t4 = (-t306 * t435 + t344) * t558 + t316 + t578 * mrSges(6,3) + t566 * pkin(5) + t565 + t611;
t347 = t4 * qJD(1) - t16 * qJD(2) + t27 * qJD(3);
t317 = (mrSges(6,3) * t530 + t384) * t196 + (mrSges(6,3) * t532 + t385) * t194 + t330 + t570;
t13 = t317 - t585 + t579;
t342 = t596 - t613;
t19 = t342 - t439;
t346 = t13 * qJD(1) + t19 * qJD(2) + qJD(3) * t336;
t341 = -t352 * t409 + t410 * t576;
t325 = (-t312 * t137 / 0.2e1 + t517 * t549 + (t312 * t545 + t517 * t543) * mrSges(6,3)) * pkin(3);
t18 = -t409 * mrSges(6,1) + t410 * mrSges(6,2) + t325;
t288 = (mrSges(6,1) * t312 + mrSges(6,2) * t517) * pkin(3);
t339 = -qJD(1) * t18 + qJD(4) * t288;
t335 = -t352 * t437 - t379 * t576 + t381 + t438;
t334 = t147 * t605 + t224 * t544 + t226 * t542 + t576 * t553 - t352 * t551 + t314 * t547 + Ifges(4,6) * t284 + (Ifges(5,5) * t313 + Ifges(6,5) * t352 + Ifges(5,6) * t314 + Ifges(6,6) * t576) * t527 + (-Ifges(4,5) + t614) * t283;
t327 = (t425 + t294 / 0.4e1 - t365 / 0.4e1) * t284 + t511 / 0.2e1 - t481 / 0.4e1 + t572;
t319 = t307 * t528 + t313 * t413 + t314 * t414 + t330 + t594 - t602;
t280 = t288 * qJD(5);
t20 = t342 + t439;
t17 = t342 + t322 + t612;
t15 = -t506 / 0.2e1 - t504 / 0.2e1 - t490 / 0.2e1 - t491 / 0.2e1 + t325 + t441;
t12 = t317 + t568 + t585;
t9 = t318 + t568 + t586;
t6 = t334 + t338 * t558 + (t475 / 0.2e1 - t476 / 0.2e1) * mrSges(6,3) + t332 - t455 / 0.2e1 + t453 / 0.2e1 + t457 / 0.2e1 - t467 / 0.2e1 + t468 / 0.2e1 + t471 / 0.2e1 - t470 / 0.2e1 + t514 / 0.2e1 + ((t377 * t559 + t343) * t284 + t331) * pkin(1) + t74 * t417 + t72 * t415 - t591;
t3 = (t341 + t578) * mrSges(6,3) + (t216 * t557 + t327) * t313 + (t218 * t557 + ((m(6) * t521 + t540) * pkin(3) + t383) * t284 + t386) * t314 + t344 * t558 + pkin(5) * t356 + t319 + t565 - t611;
t1 = (t218 * t523 + ((t289 * t559 + t540) * pkin(3) + t383) * t284 + t386) * t314 + t345 * t558 + t304 * t356 + (t216 * t523 + t327) * t313 + (t341 + t577) * mrSges(6,3) + t319 + t564 - t610;
t22 = [qJD(2) * t5 + qJD(3) * t7 + qJD(4) * t11 + qJD(5) * t14, t478 + (t334 + m(6) * (t241 * t74 - t242 * t72 - t239) + t361 * t304 + Ifges(3,5) * t519 - Ifges(3,6) * t516 + (t518 * t283 + t515 * t284) * pkin(1) * mrSges(4,3) - t455 + t457 + t471 - t470 + (-t352 * t72 + t576 * t74) * mrSges(6,3)) * qJD(2) + t6 * qJD(3) + t1 * qJD(4) + t9 * qJD(5), t477 + t6 * qJD(2) + (t334 + m(6) * (t267 * t90 - t268 * t89 - t378) + t453 - t467 + t468 + t514 + t361 * pkin(5) + (t475 - t476) * mrSges(6,3)) * qJD(3) + t3 * qJD(4) + t12 * qJD(5), t474 + t1 * qJD(2) + t3 * qJD(3) + (-t214 * t501 + Ifges(5,5) * t460 + Ifges(5,6) * t459 - mrSges(5,2) * t472 - t490 - t491 - t196 * t379 + (-t517 * t105 + t106 * t312) * t555 + t194 * t437 + t441) * qJD(4) + t15 * qJD(5), t473 + t9 * qJD(2) + t12 * qJD(3) + t15 * qJD(4) + (t441 - t504 - t506) * qJD(5); qJD(3) * t8 + qJD(4) * t2 + qJD(5) * t10 - t478, qJD(3) * t30 + qJD(4) * t21 + qJD(5) * t337, ((-t515 * pkin(2) + t348 * pkin(5)) * t556 + m(6) * (-t268 * t249 + t267 * t250 + t306 * t427) + t320) * qJD(3) + t17 * qJD(4) + t20 * qJD(5) + t363, t17 * qJD(3) + (m(6) * (-t241 * t431 - t242 * t510) + t290 * t304 + t335 + t582) * qJD(4) + t592 + t364, t20 * qJD(3) + qJD(4) * t358 + t362 + t592; -qJD(2) * t8 + qJD(4) * t4 + qJD(5) * t13 - t477, -qJD(4) * t16 + qJD(5) * t19 - t363, qJD(4) * t27 + qJD(5) * t336, (m(6) * (-t267 * t431 - t268 * t510) + t290 * pkin(5) + t335 + t581) * qJD(4) + t593 + t347, qJD(4) * t357 + t346 + t593; -qJD(2) * t2 - qJD(3) * t4 + qJD(5) * t18 - t474, qJD(3) * t16 - t364, -t347, -t280, -t280 - t339; -qJD(2) * t10 - qJD(3) * t13 - qJD(4) * t18 - t473, -qJD(3) * t19 - t362, -t346, t339, 0;];
Cq = t22;
