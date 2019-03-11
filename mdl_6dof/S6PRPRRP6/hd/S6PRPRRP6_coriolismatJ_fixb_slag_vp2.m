% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:19
% EndTime: 2019-03-08 20:18:31
% DurationCPUTime: 5.73s
% Computational Cost: add. (6013->513), mult. (14562->709), div. (0->0), fcn. (12911->8), ass. (0->272)
t317 = sin(qJ(5));
t309 = t317 ^ 2;
t320 = cos(qJ(5));
t311 = t320 ^ 2;
t518 = t309 + t311;
t503 = -m(7) / 0.2e1;
t530 = Ifges(7,4) + Ifges(6,5);
t531 = mrSges(7,2) + mrSges(6,3);
t544 = t518 * t531;
t465 = Ifges(7,5) * t320;
t262 = Ifges(7,3) * t317 + t465;
t308 = Ifges(6,4) * t320;
t519 = -Ifges(6,2) * t317 + t308;
t543 = t519 / 0.2e1 - t262 / 0.2e1;
t370 = pkin(5) * t320 + qJ(6) * t317;
t253 = -pkin(4) - t370;
t313 = sin(pkin(6));
t318 = sin(qJ(4));
t322 = cos(qJ(2));
t319 = sin(qJ(2));
t424 = t319 * t320;
t156 = (t317 * t322 + t318 * t424) * t313;
t438 = t156 * t320;
t434 = t313 * t319;
t401 = t317 * t434;
t433 = t313 * t322;
t155 = t318 * t401 - t320 * t433;
t439 = t155 * t317;
t349 = (t438 + t439) * pkin(9);
t321 = cos(qJ(4));
t400 = t321 * t434;
t506 = -m(6) / 0.2e1;
t542 = -t531 * (t438 / 0.2e1 + t439 / 0.2e1) + (pkin(4) * t400 + t349) * t506 + (-t253 * t400 + t349) * t503;
t505 = m(6) / 0.2e1;
t481 = -t321 / 0.2e1;
t480 = t321 / 0.2e1;
t314 = cos(pkin(6));
t213 = t314 * t321 - t318 * t433;
t106 = t213 * t320 + t401;
t445 = t106 * t320;
t539 = t213 - t445;
t405 = mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t502 = m(7) / 0.2e1;
t538 = (-pkin(5) * t155 + qJ(6) * t156) * t502 + t405 * t156;
t454 = t320 * mrSges(7,1);
t459 = t317 * mrSges(7,3);
t373 = t454 + t459;
t455 = t320 * mrSges(6,1);
t460 = t317 * mrSges(6,2);
t374 = t455 - t460;
t524 = -t373 - t374;
t537 = t518 * pkin(9);
t536 = m(6) / 0.4e1 + m(7) / 0.4e1;
t105 = t213 * t317 - t313 * t424;
t212 = t314 * t318 + t321 * t433;
t220 = t370 * t321;
t419 = t320 * t321;
t241 = t318 * mrSges(6,1) - mrSges(6,3) * t419;
t242 = -t318 * mrSges(7,1) + mrSges(7,2) * t419;
t386 = t242 / 0.2e1 - t241 / 0.2e1;
t323 = -pkin(2) - pkin(8);
t379 = t317 * t323 - pkin(5);
t476 = pkin(4) * t318;
t252 = -pkin(9) * t321 + qJ(3) + t476;
t421 = t320 * t252;
t119 = t318 * t379 - t421;
t425 = t318 * t323;
t141 = -t317 * t425 + t421;
t415 = t119 + t141;
t418 = t320 * t323;
t430 = t317 * t252;
t116 = t430 + (qJ(6) + t418) * t318;
t142 = t318 * t418 + t430;
t416 = -t116 + t142;
t535 = (t416 * t105 + t415 * t106 + t212 * t220) * t503 - t386 * t106;
t533 = m(6) + m(7);
t532 = mrSges(6,1) + mrSges(7,1);
t529 = Ifges(7,2) + Ifges(6,3);
t528 = t318 * (-0.1e1 + t518);
t404 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t526 = t404 * t318;
t523 = m(7) * t253 - t373;
t453 = t320 * mrSges(6,2);
t462 = t317 * mrSges(6,1);
t260 = t453 + t462;
t447 = qJ(6) * t320;
t475 = pkin(5) * t317;
t258 = -t447 + t475;
t478 = m(7) * t258;
t522 = t260 + t478;
t305 = Ifges(7,5) * t317;
t521 = Ifges(7,1) * t320 + t305;
t520 = Ifges(7,6) * t317 + t530 * t320;
t267 = Ifges(6,1) * t317 + t308;
t417 = t321 * t323;
t473 = pkin(9) * t318;
t273 = pkin(4) * t321 + t473;
t420 = t320 * t273;
t151 = -t317 * t417 + t420;
t152 = t317 * t273 + t320 * t417;
t360 = -t151 * t317 + t152 * t320;
t126 = qJ(6) * t321 + t152;
t128 = t321 * t379 - t420;
t363 = t126 * t320 + t128 * t317;
t515 = Ifges(7,3) * t320 - t305;
t452 = t320 * mrSges(7,3);
t461 = t317 * mrSges(7,1);
t514 = -t478 / 0.2e1 - t462 / 0.2e1 - t461 / 0.2e1 - t453 / 0.2e1 + t452 / 0.2e1;
t513 = Ifges(6,6) * t481 + Ifges(7,6) * t480 + t318 * t543;
t301 = m(7) * qJ(6) + mrSges(7,3);
t512 = t321 * t536;
t259 = -t452 + t461;
t406 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t484 = -t260 / 0.2e1;
t511 = t406 * t317 - t405 * t320 - t259 / 0.2e1 + t484;
t223 = t259 * t318;
t224 = t260 * t318;
t428 = t317 * t321;
t238 = -mrSges(6,2) * t318 - mrSges(6,3) * t428;
t457 = t318 * mrSges(7,3);
t243 = -mrSges(7,2) * t428 + t457;
t385 = t243 / 0.2e1 + t238 / 0.2e1;
t337 = -t386 * t317 - t385 * t320 - t223 / 0.2e1 - t224 / 0.2e1;
t354 = t258 - t323;
t158 = t354 * t321;
t157 = t354 * t318;
t352 = t116 * t320 + t119 * t317 + t157;
t361 = -t141 * t317 + t142 * t320;
t429 = t317 * t318;
t237 = -t321 * mrSges(6,2) + mrSges(6,3) * t429;
t244 = mrSges(7,2) * t429 + t321 * mrSges(7,3);
t384 = t244 / 0.2e1 + t237 / 0.2e1;
t426 = t318 * t320;
t239 = t321 * mrSges(6,1) + mrSges(6,3) * t426;
t450 = t321 * mrSges(7,1);
t240 = -mrSges(7,2) * t426 - t450;
t490 = t240 / 0.2e1;
t387 = t490 - t239 / 0.2e1;
t225 = t321 * t259;
t226 = t321 * t260;
t388 = t225 / 0.2e1 + t226 / 0.2e1;
t510 = t388 * t213 + t387 * t105 + t384 * t106 + (-t151 * t105 + t152 * t106 - t213 * t417) * t505 + (t105 * t128 + t106 * t126 + t158 * t213) * t502 + ((-t361 + t425) * t505 + t352 * t503 + t337) * t212;
t509 = 0.2e1 * m(7);
t508 = 2 * qJD(4);
t507 = m(5) / 0.2e1;
t498 = Ifges(7,6) / 0.2e1;
t497 = -t105 / 0.2e1;
t496 = t128 / 0.2e1;
t495 = t142 / 0.2e1;
t493 = t212 / 0.2e1;
t221 = t373 * t321;
t492 = t221 / 0.2e1;
t222 = t374 * t321;
t491 = -t222 / 0.2e1;
t487 = -t373 / 0.2e1;
t486 = t258 / 0.2e1;
t485 = t259 / 0.2e1;
t477 = m(7) * t320;
t446 = t105 * t317;
t16 = -0.2e1 * t536 * t212 * t528 + 0.2e1 * (t446 - t539) * t512;
t470 = t16 * qJD(4);
t469 = m(7) * qJD(6);
t466 = Ifges(6,4) * t317;
t458 = t318 * mrSges(5,2);
t456 = t318 * Ifges(6,6);
t451 = t321 * mrSges(5,1);
t449 = -t374 - mrSges(5,1);
t340 = t531 * (-t309 / 0.2e1 - t311 / 0.2e1);
t362 = t141 * t320 + t142 * t317;
t364 = t116 * t317 - t119 * t320;
t390 = -t221 / 0.2e1 + t491;
t422 = t320 * t242;
t423 = t320 * t241;
t431 = t317 * t243;
t432 = t317 * t238;
t324 = (t220 * t503 + t390) * t321 + (t422 / 0.2e1 + (t362 - t364) * t502 - t431 / 0.2e1 - t432 / 0.2e1 - t423 / 0.2e1 + t340 * t321) * t318;
t333 = t370 * t503 + t460 / 0.2e1 - t459 / 0.2e1 - t455 / 0.2e1 - t454 / 0.2e1;
t14 = t324 + t333;
t442 = t14 * qJD(2);
t18 = m(5) * (-t212 * t321 + t213 * t318 + t433) * t434 + t533 * (t105 * t155 + t106 * t156 - t212 * t400);
t437 = t18 * qJD(1);
t436 = t212 * t317;
t427 = t318 * t106;
t414 = t537 * t212;
t413 = t537 * t321;
t412 = t318 * mrSges(5,1) + t321 * mrSges(5,2);
t310 = t318 ^ 2;
t312 = t321 ^ 2;
t409 = t310 + t312;
t143 = (t310 / 0.2e1 + t312 / 0.2e1 + 0.1e1 / 0.2e1) * t477;
t408 = t143 * qJD(2);
t403 = t498 - Ifges(6,6) / 0.2e1;
t402 = mrSges(4,3) + t412;
t399 = t323 * t434;
t398 = t212 * t419;
t395 = -t434 / 0.2e1;
t392 = t373 * t481;
t290 = Ifges(7,5) * t419;
t188 = Ifges(7,6) * t318 + Ifges(7,3) * t428 + t290;
t190 = t321 * t519 + t456;
t391 = -t188 / 0.2e1 + t190 / 0.2e1;
t263 = Ifges(6,2) * t320 + t466;
t382 = t515 / 0.2e1 + t263 / 0.2e1;
t265 = Ifges(7,1) * t317 - t465;
t381 = -t265 / 0.2e1 - t267 / 0.2e1;
t268 = Ifges(6,1) * t320 - t466;
t6 = t390 * t212 - t406 * t155 + t385 * t105 + t531 * t321 * (t446 / 0.2e1 + t445 / 0.2e1) + t535 + t538;
t227 = t515 * t321;
t228 = t321 * t263;
t229 = -Ifges(7,1) * t428 + t290;
t230 = t321 * t267;
t289 = Ifges(7,6) * t419;
t192 = t318 * Ifges(7,4) + t321 * t521;
t194 = t318 * Ifges(6,5) + t321 * t268;
t343 = -t192 / 0.2e1 - t194 / 0.2e1 + t526;
t8 = t142 * t242 + m(7) * (t116 * t141 + t119 * t142 + t158 * t220) + t141 * t243 + t220 * t225 + t158 * t221 + t141 * t238 - t142 * t241 + t318 * t289 / 0.2e1 + (-t323 * t222 + (-t456 / 0.2e1 - t116 * mrSges(7,2) - t142 * mrSges(6,3) + t229 / 0.2e1 - t230 / 0.2e1 - t391) * t320 + (-t119 * mrSges(7,2) + t141 * mrSges(6,3) + t228 / 0.2e1 + t227 / 0.2e1 + t343) * t317) * t321;
t369 = -t6 * qJD(1) + t8 * qJD(2);
t277 = t312 * t434;
t330 = (t310 * t434 + t277) * t507 + (t505 + t502) * (t155 * t429 + t156 * t426 + t277);
t334 = m(5) * t395 + (t506 + t503) * (-t105 * t320 + t106 * t317);
t20 = t330 + t334;
t26 = -m(6) * t362 - m(7) * t364 - t402 + t422 - t423 - t431 - t432 + (-m(4) - m(5)) * qJ(3);
t368 = -qJD(1) * t20 - qJD(2) * t26;
t38 = m(7) * (t116 * t318 - t158 * t419) + t318 * t243 - t225 * t419;
t42 = (t155 / 0.4e1 + t398 / 0.4e1 - t427 / 0.4e1) * t509;
t367 = qJD(1) * t42 - qJD(2) * t38;
t17 = t533 * (-t105 * t436 + t539 * t212);
t366 = t17 * qJD(1) + t16 * qJD(3);
t112 = t523 * t317;
t336 = (-t158 * t317 + (-t253 * t321 + t473) * t320) * t502 - t317 * t225 / 0.2e1;
t355 = m(7) * t496 - t450 / 0.2e1;
t35 = (t318 * mrSges(7,2) - t392) * t320 + t336 - t355;
t359 = -qJD(2) * t35 + qJD(4) * t112;
t55 = t457 + (t430 / 0.4e1 - t142 / 0.4e1 + (t418 / 0.4e1 + qJ(6) / 0.2e1) * t318) * t509;
t358 = qJD(2) * t55 + qJD(5) * t301;
t10 = (t352 * t502 + t361 * t505 - t337) * t321 + (t384 * t320 + t387 * t317 + (t158 + t363) * t502 + (t360 - 0.2e1 * t417) * t505 + t388) * t318;
t257 = t451 - t458;
t2 = (t257 / 0.2e1 + t458 / 0.2e1 + (-mrSges(5,1) / 0.2e1 - t374 / 0.2e1 + t487) * t321) * t434 + t510 + t542;
t191 = t321 * Ifges(7,4) - t318 * t521;
t193 = t321 * Ifges(6,5) - t268 * t318;
t5 = qJ(3) * t257 + t116 * t244 + t119 * t240 + t126 * t243 + t128 * t242 + t141 * t239 + t142 * t237 + t151 * t241 + t152 * t238 - t157 * t225 - t158 * t223 + m(6) * (t141 * t151 + t142 * t152) + m(7) * (t116 * t126 + t119 * t128 - t157 * t158) + (-Ifges(5,4) * t321 + t323 * t224 + (t191 / 0.2e1 + t193 / 0.2e1 - t404 * t321) * t320 + (t403 * t321 + t513) * t317) * t321 + (Ifges(5,4) * t318 + t323 * t226 + t343 * t320 + (-t403 * t318 + t391) * t317 + (-m(6) * t323 ^ 2 - Ifges(5,1) + Ifges(5,2) + t529) * t321) * t318;
t346 = t2 * qJD(1) + t5 * qJD(2) + t10 * qJD(3);
t62 = 0.4e1 * t512 * t528;
t345 = t16 * qJD(1) + t10 * qJD(2) + t62 * qJD(3);
t11 = ((t475 / 0.4e1 - t447 / 0.4e1 - t258 / 0.4e1) * t509 + t511) * t212;
t28 = -pkin(4) * t260 - t258 * t373 + (t259 + t478) * t253 + (-t381 + t543) * t320 + (t268 / 0.2e1 + t521 / 0.2e1 - t382) * t317;
t326 = t323 * t484 + (t268 / 0.4e1 + t521 / 0.4e1 - t263 / 0.4e1 - t515 / 0.4e1) * t320 + (-t519 / 0.4e1 + t262 / 0.4e1 - t267 / 0.4e1 - t265 / 0.4e1) * t317 + t340 * pkin(9);
t327 = (t158 * t258 + t220 * t253) * t502 + pkin(4) * t491 + t158 * t485 + t220 * t487 + t253 * t492 + t225 * t486 + t520 * t318 / 0.4e1;
t329 = (-pkin(5) * t128 + qJ(6) * t126) * t503 + pkin(5) * t490 - qJ(6) * t244 / 0.2e1 - t126 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t496 - t151 * mrSges(6,1) / 0.2e1 + t152 * mrSges(6,2) / 0.2e1;
t331 = (t495 - t116 / 0.2e1) * mrSges(7,2) + (t416 * t502 - t385) * pkin(9) + t188 / 0.4e1 - t190 / 0.4e1 + t229 / 0.4e1 - t230 / 0.4e1;
t332 = -t228 / 0.4e1 - t227 / 0.4e1 + t194 / 0.4e1 + t192 / 0.4e1 + (t141 / 0.2e1 + t119 / 0.2e1) * mrSges(7,2) + (t415 * t502 + t386) * pkin(9);
t3 = (t332 - t526) * t320 + ((-0.3e1 / 0.4e1 * Ifges(6,6) + t498) * t318 + t331) * t317 + t327 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1 + t326) * t321 + t329;
t40 = ((-t475 / 0.2e1 + t447 / 0.2e1 + t486) * m(7) - t511) * t321;
t341 = t11 * qJD(1) - t3 * qJD(2) + t40 * qJD(3) - t28 * qJD(4);
t335 = (-m(7) * t370 + t524) * qJD(5);
t288 = qJ(3) * t433;
t254 = (m(7) * pkin(9) + mrSges(7,2)) * t320;
t247 = t312 * t399;
t144 = (-0.1e1 / 0.2e1 + t409 / 0.2e1) * t477;
t78 = m(7) * t436;
t48 = (t430 + (0.2e1 * qJ(6) + t418) * t318) * t502 + m(7) * t495 + t243;
t43 = (-t398 + t427 + t155) * t502;
t41 = t514 * t321 + (t259 + t522) * t481;
t37 = -t320 * t392 + t336 + t355;
t19 = m(4) * t434 + t330 - t334;
t13 = t324 - t333;
t12 = t522 * t493 + (t485 - t514) * t212;
t9 = t10 * qJD(4);
t7 = t212 * t492 + t222 * t493 + (t238 + t243) * t497 - t532 * t155 / 0.2e1 + t531 * (t428 * t497 - t106 * t419 / 0.2e1) - t535 + t538;
t4 = (-t456 / 0.4e1 + t331) * t317 + t332 * t320 + t327 - t329 + t326 * t321 + t529 * t480 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t429 - t530 * t426 / 0.2e1;
t1 = (t451 + t257) * t434 / 0.2e1 + t510 + (t524 * t321 + t458) * t395 - t542;
t15 = [qJD(2) * t18 + qJD(4) * t17, t19 * qJD(3) + t1 * qJD(4) + t7 * qJD(5) + t43 * qJD(6) + t437 + (-t155 * t241 + t155 * t242 + t156 * t238 + t156 * t243 + ((-mrSges(3,2) + t402) * t322 + (-mrSges(3,1) + mrSges(4,2) + (-t225 - t226) * t321 - t409 * mrSges(5,3)) * t319) * t313 + 0.2e1 * (-t141 * t155 + t142 * t156 + t247) * t505 + 0.2e1 * (t116 * t156 + t119 * t155 - t158 * t400) * t502 + 0.2e1 * (t310 * t399 + t247 + t288) * t507 + m(4) * (-pkin(2) * t434 + t288)) * qJD(2), qJD(2) * t19 + t470, t1 * qJD(2) + t12 * qJD(5) - t78 * qJD(6) + ((-pkin(4) * t213 - t414) * t505 + (t213 * t253 - t414) * t502) * t508 + t366 + ((-t373 + t449) * t213 + (mrSges(5,2) - t544) * t212) * qJD(4), t7 * qJD(2) + t12 * qJD(4) + (-t532 * t106 + (mrSges(6,2) - mrSges(7,3)) * t105) * qJD(5) + ((-pkin(5) * t106 - qJ(6) * t105) * qJD(5) / 0.2e1 + t106 * qJD(6) / 0.2e1) * t509, m(7) * t106 * qJD(5) + t43 * qJD(2) - t78 * qJD(4); -qJD(3) * t20 + qJD(4) * t2 - qJD(5) * t6 - qJD(6) * t42 - t437, -qJD(3) * t26 + qJD(4) * t5 + qJD(5) * t8 + qJD(6) * t38, qJD(5) * t13 + qJD(6) * t144 + t368 + t9, t4 * qJD(5) + t37 * qJD(6) + t346 + (-mrSges(5,2) * t417 + pkin(4) * t224 - t253 * t223 - Ifges(5,6) * t321 + (m(6) * t360 + m(7) * t363) * pkin(9) + (-Ifges(5,5) + (-m(6) * pkin(4) + t449) * t323) * t318 + (t191 + t193) * t317 / 0.2e1 + ((-t239 + t240) * pkin(9) + t382 * t318 + t530 * t480) * t317 + ((t237 + t244) * pkin(9) + t381 * t318 + (-Ifges(7,6) + Ifges(6,6)) * t480 - t513) * t320 - t523 * t157 + t360 * mrSges(6,3) + t363 * mrSges(7,2)) * qJD(4), t13 * qJD(3) + t4 * qJD(4) + t48 * qJD(6) + t369 + (t289 + (-m(7) * pkin(5) - t532) * t142 + (-mrSges(6,2) + t301) * t141 + ((-mrSges(7,2) * qJ(6) - Ifges(6,6)) * t320 + (mrSges(7,2) * pkin(5) - t530) * t317) * t321) * qJD(5), qJD(3) * t144 + qJD(4) * t37 + qJD(5) * t48 - t367; qJD(2) * t20 + t470, qJD(5) * t14 + qJD(6) * t143 - t368 + t9, t62 * qJD(4) (t524 * t318 - t412) * qJD(4) + t41 * qJD(5) + ((t253 * t318 + t413) * t502 + (t413 - t476) * t505) * t508 + (qJD(4) * t544 + t317 * t469) * t321 + t345, t442 + t41 * qJD(4) + (t320 * t469 + t335) * t318, t408 + (qJD(4) * t428 + qJD(5) * t426) * m(7); -qJD(2) * t2 - qJD(5) * t11 - t366, qJD(5) * t3 + qJD(6) * t35 - t346, -qJD(5) * t40 - t345, qJD(5) * t28 - qJD(6) * t112 (-t370 * mrSges(7,2) - Ifges(6,6) * t317 + t520) * qJD(5) + t254 * qJD(6) + pkin(9) * t335 - t341, qJD(5) * t254 - t359; qJD(2) * t6 + qJD(4) * t11, -qJD(3) * t14 - qJD(4) * t3 + qJD(6) * t55 - t369, qJD(4) * t40 - t442, t341, t301 * qJD(6), t358; qJD(2) * t42, -qJD(3) * t143 - qJD(4) * t35 - qJD(5) * t55 + t367, -t408, t359, -t358, 0;];
Cq  = t15;
