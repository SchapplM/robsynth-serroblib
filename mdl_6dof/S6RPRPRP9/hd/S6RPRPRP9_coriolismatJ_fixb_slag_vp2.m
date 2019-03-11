% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:38
% EndTime: 2019-03-09 03:27:49
% DurationCPUTime: 6.62s
% Computational Cost: add. (11417->592), mult. (23929->777), div. (0->0), fcn. (23546->6), ass. (0->296)
t526 = Ifges(7,4) + Ifges(6,5);
t531 = Ifges(6,6) - Ifges(7,6);
t530 = mrSges(7,2) + mrSges(6,3);
t333 = sin(qJ(3));
t334 = cos(qJ(3));
t469 = pkin(3) * t333;
t305 = -qJ(4) * t334 + qJ(2) + t469;
t332 = cos(pkin(9));
t291 = t332 * t305;
t331 = sin(pkin(9));
t335 = -pkin(1) - pkin(7);
t381 = -t331 * t335 + pkin(4);
t414 = t332 * t334;
t183 = -pkin(8) * t414 + t333 * t381 + t291;
t324 = t333 * t335;
t225 = t331 * t305 + t332 * t324;
t417 = t331 * t334;
t197 = -pkin(8) * t417 + t225;
t471 = sin(qJ(5));
t472 = cos(qJ(5));
t71 = t183 * t472 - t197 * t471;
t64 = -t333 * pkin(5) - t71;
t463 = t64 + t71;
t348 = t331 * t472 + t332 * t471;
t264 = t348 * t333;
t217 = t264 * mrSges(7,2) + mrSges(7,3) * t334;
t218 = -mrSges(6,2) * t334 + t264 * mrSges(6,3);
t529 = t217 + t218;
t308 = pkin(3) * t334 + qJ(4) * t333;
t300 = t332 * t308;
t410 = t334 * t335;
t229 = -t331 * t410 + t300;
t230 = t331 * t308 + t332 * t410;
t359 = -t229 * t331 + t230 * t332;
t528 = m(6) + m(7);
t295 = t471 * t331 - t472 * t332;
t420 = t295 * qJ(6);
t467 = t348 * pkin(5);
t367 = -t420 - t467;
t470 = m(7) * t367;
t527 = -mrSges(6,1) - mrSges(7,1);
t525 = Ifges(7,2) + Ifges(6,3);
t265 = t348 * t334;
t524 = mrSges(6,3) * t265;
t267 = t295 * t334;
t523 = mrSges(6,3) * t267;
t413 = t333 * qJ(6);
t393 = t471 * t183;
t395 = t472 * t197;
t72 = t395 + t393;
t63 = t72 + t413;
t498 = -mrSges(6,3) / 0.2e1;
t398 = t498 - mrSges(7,2) / 0.2e1;
t522 = t267 * t398;
t438 = t332 * mrSges(5,2);
t441 = t331 * mrSges(5,1);
t521 = t333 * (t438 + t441);
t153 = -t267 * mrSges(7,1) + t265 * mrSges(7,3);
t154 = -t267 * mrSges(6,1) - t265 * mrSges(6,2);
t520 = t154 + t153;
t199 = mrSges(7,1) * t348 + t295 * mrSges(7,3);
t200 = mrSges(6,1) * t348 - t295 * mrSges(6,2);
t519 = t199 + t200;
t435 = t333 * mrSges(7,3);
t216 = -t265 * mrSges(7,2) + t435;
t219 = -t333 * mrSges(6,2) - t524;
t518 = t216 + t219;
t222 = t333 * mrSges(6,1) + t523;
t223 = -mrSges(7,1) * t333 - t267 * mrSges(7,2);
t407 = -t222 + t223;
t516 = -t295 * t526 - t531 * t348;
t515 = -t265 * t526 + t531 * t267;
t501 = m(7) / 0.2e1;
t503 = m(6) / 0.2e1;
t513 = t501 + t503;
t328 = m(7) * qJ(6) + mrSges(7,3);
t224 = -t324 * t331 + t291;
t302 = -t333 * mrSges(5,2) - mrSges(5,3) * t417;
t304 = t333 * mrSges(5,1) - mrSges(5,3) * t414;
t512 = -m(5) * (t224 * t332 + t225 * t331) - t331 * t302 - t332 * t304;
t511 = -m(7) * pkin(5) + t527;
t510 = -mrSges(6,2) + t328;
t266 = t295 * t333;
t396 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t397 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t509 = -t396 * t264 + t397 * t266;
t508 = 0.2e1 * m(7);
t507 = 2 * qJD(3);
t506 = -m(5) / 0.2e1;
t505 = m(5) / 0.2e1;
t504 = -m(6) / 0.2e1;
t502 = -m(7) / 0.2e1;
t500 = mrSges(7,1) / 0.2e1;
t499 = mrSges(6,2) / 0.2e1;
t497 = -mrSges(7,3) / 0.2e1;
t424 = t265 * qJ(6);
t468 = pkin(5) * t267;
t152 = t424 - t468;
t496 = t152 / 0.2e1;
t495 = t153 / 0.2e1;
t494 = t154 / 0.2e1;
t448 = t267 * mrSges(7,3);
t453 = t265 * mrSges(7,1);
t157 = t448 + t453;
t493 = t157 / 0.2e1;
t492 = t199 / 0.2e1;
t491 = t200 / 0.2e1;
t490 = -t217 / 0.2e1;
t489 = -t264 / 0.2e1;
t488 = t264 / 0.2e1;
t487 = -t265 / 0.2e1;
t486 = t265 / 0.2e1;
t485 = t266 / 0.2e1;
t484 = -t266 / 0.2e1;
t482 = -t295 / 0.2e1;
t481 = t295 / 0.2e1;
t479 = t331 / 0.2e1;
t478 = -t332 / 0.2e1;
t477 = t332 / 0.2e1;
t476 = t333 / 0.2e1;
t474 = -t334 / 0.2e1;
t473 = t334 / 0.2e1;
t465 = pkin(8) + qJ(4);
t464 = -t63 + t72;
t462 = Ifges(5,4) * t331;
t461 = Ifges(5,4) * t332;
t460 = Ifges(6,4) * t267;
t459 = Ifges(6,4) * t348;
t458 = Ifges(7,5) * t265;
t457 = Ifges(7,5) * t295;
t456 = t264 * mrSges(6,1);
t455 = t264 * mrSges(7,1);
t454 = t265 * mrSges(6,1);
t452 = t266 * mrSges(6,2);
t451 = t266 * mrSges(7,2);
t450 = t266 * mrSges(7,3);
t449 = t267 * mrSges(6,2);
t447 = t295 * mrSges(6,1);
t446 = t295 * mrSges(7,1);
t445 = t295 * mrSges(7,2);
t444 = t348 * mrSges(6,2);
t443 = t348 * mrSges(7,3);
t140 = Ifges(7,5) * t266 + Ifges(7,6) * t334 - Ifges(7,3) * t264;
t142 = Ifges(6,4) * t266 + Ifges(6,2) * t264 + Ifges(6,6) * t334;
t144 = Ifges(7,1) * t266 + Ifges(7,4) * t334 - Ifges(7,5) * t264;
t146 = Ifges(6,1) * t266 + Ifges(6,4) * t264 + Ifges(6,5) * t334;
t155 = -t450 - t455;
t156 = t452 - t456;
t220 = mrSges(6,1) * t334 - mrSges(6,3) * t266;
t434 = t334 * mrSges(7,1);
t221 = -t434 + t451;
t440 = t331 * Ifges(5,2);
t260 = Ifges(5,6) * t334 + (t440 - t461) * t333;
t437 = t332 * Ifges(5,1);
t261 = Ifges(5,5) * t334 + (-t437 + t462) * t333;
t418 = t331 * t333;
t292 = -pkin(4) * t418 + t324;
t293 = pkin(4) * t417 - t410;
t301 = -mrSges(5,2) * t334 + mrSges(5,3) * t418;
t415 = t332 * t333;
t303 = mrSges(5,1) * t334 + mrSges(5,3) * t415;
t145 = -Ifges(7,1) * t267 + t333 * Ifges(7,4) + t458;
t253 = Ifges(6,4) * t265;
t147 = -Ifges(6,1) * t267 + t333 * Ifges(6,5) - t253;
t385 = -t145 / 0.2e1 - t147 / 0.2e1;
t250 = Ifges(7,5) * t267;
t141 = t333 * Ifges(7,6) + Ifges(7,3) * t265 - t250;
t143 = -Ifges(6,2) * t265 + t333 * Ifges(6,6) - t460;
t386 = t143 / 0.2e1 - t141 / 0.2e1;
t436 = t332 * Ifges(5,5);
t439 = t331 * Ifges(5,6);
t192 = pkin(8) * t415 + t334 * t381 + t300;
t211 = pkin(8) * t418 + t230;
t74 = t471 * t192 + t472 * t211;
t65 = qJ(6) * t334 + t74;
t73 = t192 * t472 - t211 * t471;
t66 = -t334 * pkin(5) - t73;
t97 = -pkin(5) * t264 - qJ(6) * t266 + t292;
t368 = pkin(5) * t265 + t267 * qJ(6);
t98 = t293 + t368;
t3 = t98 * t155 + t293 * t156 + t97 * t157 + t65 * t216 + t63 * t217 + t72 * t218 + t74 * t219 + t71 * t220 + t64 * t221 + t73 * t222 + t66 * t223 + t224 * t303 + t225 * t301 + t229 * t304 + t230 * t302 - t385 * t266 - (t144 / 0.2e1 + t146 / 0.2e1 + t292 * mrSges(6,2)) * t267 + (t140 / 0.2e1 - t142 / 0.2e1 + t292 * mrSges(6,1)) * t265 + t386 * t264 + m(6) * (t292 * t293 + t71 * t73 + t72 * t74) + m(7) * (t63 * t65 + t64 * t66 + t97 * t98) + m(5) * (t224 * t229 + t225 * t230) + (t335 * t521 + t261 * t477 - t331 * t260 / 0.2e1 + qJ(2) * mrSges(4,1) + (-Ifges(4,4) + t436 / 0.2e1 - t439 / 0.2e1) * t334 - t397 * t267 + t396 * t265) * t334 + (-qJ(2) * mrSges(4,2) + (Ifges(4,4) - t436 + t439) * t333 + (-m(5) * t335 ^ 2 + Ifges(4,2) - Ifges(4,1) + Ifges(5,3) + (t335 * mrSges(5,2) - t437 / 0.2e1) * t332 + (t335 * mrSges(5,1) + t461 - t440 / 0.2e1) * t331 + t525) * t334 + t509) * t333;
t442 = t3 * qJD(1);
t158 = -Ifges(7,3) * t267 - t458;
t159 = Ifges(6,2) * t267 - t253;
t160 = -Ifges(7,1) * t265 - t250;
t161 = -Ifges(6,1) * t265 + t460;
t4 = t293 * t154 + t98 * t153 - (t160 / 0.2e1 + t161 / 0.2e1 - t63 * mrSges(7,2) - t386) * t267 + (t158 / 0.2e1 - t159 / 0.2e1 - t64 * mrSges(7,2) + t385) * t265 + (m(7) * t98 + t157) * t152 + (m(7) * t64 + t407 + t523) * t72 + (m(7) * t63 + t518 + t524) * t71 + t515 * t476;
t433 = t4 * qJD(1);
t432 = -mrSges(5,1) * t332 + mrSges(5,2) * t331 - mrSges(4,1);
t27 = t267 * t157 + m(7) * (t267 * t98 + t333 * t63) + t333 * t216;
t431 = qJD(1) * t27;
t37 = t513 * (-t264 * t267 + t265 * t266);
t430 = qJD(1) * t37;
t15 = t333 * mrSges(4,1) + t334 * mrSges(4,2) + mrSges(3,3) + t518 * t348 + t407 * t295 + (m(4) + m(3)) * qJ(2) + m(7) * (t64 * t295 + t348 * t63) + m(6) * (-t71 * t295 + t348 * t72) - t512;
t429 = t15 * qJD(1);
t307 = t465 * t332;
t382 = t465 * t331;
t214 = t307 * t471 + t382 * t472;
t428 = t214 * t265;
t215 = t307 * t472 - t382 * t471;
t115 = t215 * t267;
t425 = t264 * t265;
t171 = t266 * t267;
t422 = t266 * t333;
t419 = t331 * t303;
t416 = t332 * t301;
t412 = t333 * t334;
t411 = t334 * t267;
t78 = (t295 / 0.4e1 + t422 / 0.4e1 + t411 / 0.4e1) * t508;
t409 = t78 * qJD(1);
t329 = t331 ^ 2;
t330 = t332 ^ 2;
t402 = t329 + t330;
t401 = qJD(3) * t333;
t400 = m(5) * t476;
t399 = -t470 / 0.2e1;
t320 = -pkin(4) * t332 - pkin(3);
t389 = t348 * t474;
t388 = -t418 / 0.2e1;
t387 = t333 * t482;
t384 = -t216 / 0.2e1 - t219 / 0.2e1;
t383 = -t223 / 0.2e1 + t222 / 0.2e1;
t380 = t402 * mrSges(5,3);
t378 = t402 * qJ(4);
t377 = -t214 * t267 - t215 * t265;
t376 = -t115 + t428;
t370 = t398 * t265;
t366 = pkin(5) * t295 - qJ(6) * t348;
t31 = m(5) * (-0.1e1 + t402) * t412 + t528 * (t171 - t412 + t425);
t360 = -t224 * t331 + t225 * t332;
t338 = (t360 * t334 + (t359 - 0.2e1 * t410) * t333) * t506 + (-t264 * t73 - t265 * t71 - t266 * t74 - t267 * t72 - t292 * t334 + t293 * t333) * t504 + (t264 * t66 + t265 * t64 - t266 * t65 - t267 * t63 + t333 * t98 - t334 * t97) * t502;
t344 = t513 * (t214 * t295 + t215 * t348);
t350 = t302 * t478 + t304 * t479;
t352 = -t454 / 0.2e1 + t449 / 0.2e1;
t5 = -t384 * t267 - (t490 - t218 / 0.2e1) * t266 + t383 * t265 + (-t221 / 0.2e1 + t220 / 0.2e1) * t264 + (-t157 / 0.2e1 - t416 / 0.2e1 + t419 / 0.2e1 + t352) * t333 + (t156 / 0.2e1 + t155 / 0.2e1 - t521 + t350) * t334 + t338 + t344;
t365 = -t5 * qJD(1) + t31 * qJD(2);
t340 = t366 * t502 - t447 / 0.2e1 - t446 / 0.2e1 - t444 / 0.2e1 + t443 / 0.2e1;
t7 = (m(7) * t496 + t494 + t495) * t334 - (t463 * t502 + t383 + t522) * t266 + (t464 * t502 - t370 - t384) * t264 + t340;
t364 = t7 * qJD(1);
t14 = -t407 * t267 - t518 * t265 + m(7) * (-t265 * t63 - t267 * t64) + m(6) * (-t265 * t72 + t267 * t71) + t512 * t334;
t363 = qJD(1) * t14 + qJD(2) * t37;
t33 = (t468 / 0.4e1 - t424 / 0.4e1 - t152 / 0.4e1) * t508 - t520;
t42 = (-t467 / 0.4e1 - t420 / 0.4e1 + t367 / 0.4e1) * t508 - t519;
t362 = qJD(1) * t33 + qJD(3) * t42;
t30 = t435 + (t395 / 0.4e1 + t393 / 0.4e1 + t413 / 0.2e1 - t72 / 0.4e1) * t508;
t356 = qJD(1) * t30 + qJD(5) * t328;
t117 = m(7) * t267;
t175 = m(7) * t348;
t355 = -qJD(1) * t117 + qJD(3) * t175;
t354 = t292 * t504 + t502 * t97;
t353 = t66 * t502 + t434 / 0.2e1;
t173 = t320 + t366;
t201 = -t443 + t446;
t203 = Ifges(7,3) * t348 - t457;
t286 = Ifges(7,5) * t348;
t204 = Ifges(7,3) * t295 + t286;
t289 = Ifges(6,4) * t295;
t205 = -Ifges(6,2) * t348 - t289;
t206 = -Ifges(6,2) * t295 + t459;
t207 = -Ifges(7,1) * t295 + t286;
t208 = Ifges(7,1) * t348 + t457;
t209 = -Ifges(6,1) * t295 - t459;
t210 = Ifges(6,1) * t348 - t289;
t336 = t384 * t214 - (t143 / 0.4e1 - t161 / 0.4e1 - t160 / 0.4e1 - t141 / 0.4e1) * t348 + (t158 / 0.4e1 - t159 / 0.4e1 - t147 / 0.4e1 - t145 / 0.4e1) * t295 + (t115 / 0.2e1 - t428 / 0.2e1 - (-t72 / 0.2e1 + t63 / 0.2e1) * t348 + (-t64 / 0.2e1 - t71 / 0.2e1) * t295) * mrSges(7,2) + (t214 * t498 + t203 / 0.4e1 - t210 / 0.4e1 - t208 / 0.4e1 - t205 / 0.4e1) * t265 - (t209 / 0.4e1 + t207 / 0.4e1 + t204 / 0.4e1 - t206 / 0.4e1) * t267 + (t152 * t173 + t214 * t464 - t367 * t98) * t501 + t201 * t496 + t173 * t495 - t367 * t493 + t293 * t491 + t320 * t494 + t98 * t492 + t516 * t333 / 0.4e1 + (-t498 * t267 + t463 * t501 - t383) * t215;
t339 = (-pkin(5) * t66 + qJ(6) * t65) * t502 + pkin(5) * t221 / 0.2e1 + qJ(6) * t490 + t65 * t497 + t66 * t500 - t73 * mrSges(6,1) / 0.2e1 + t74 * t499;
t1 = (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t334 + t336 + t339 - t509;
t11 = -t367 * t201 + t320 * t200 - (-t204 / 0.2e1 - t209 / 0.2e1 + t206 / 0.2e1 - t207 / 0.2e1) * t348 + (t203 / 0.2e1 - t210 / 0.2e1 - t205 / 0.2e1 - t208 / 0.2e1) * t295 + (t199 - t470) * t173;
t341 = t368 * t502 - t453 / 0.2e1 - t448 / 0.2e1 + t352;
t12 = (t492 + t491 + t399) * t334 + t341;
t347 = t1 * qJD(1) - t12 * qJD(2) + t11 * qJD(3);
t337 = (-t370 + t384) * t295 - (t383 - t522) * t348 + t360 * t505 + (-t295 * t72 - t348 * t71 + t377) * t503 + (-t295 * t63 + t348 * t64 + t377) * t501 - t350;
t10 = -(t499 + t497) * t266 + (mrSges(6,1) / 0.2e1 + t500) * t264 + (t335 * t506 + t438 / 0.2e1 + t441 / 0.2e1) * t333 + t337 + t354;
t24 = m(5) * t378 + t380 + (t295 ^ 2 + t348 ^ 2) * t530 + t528 * (t214 * t348 - t215 * t295);
t343 = t513 * (t264 * t348 + t295 * t266);
t35 = (t502 + t504 + (t329 / 0.2e1 + t330 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t333 + t343;
t346 = qJD(1) * t10 + qJD(2) * t35 + qJD(3) * t24;
t149 = (-t389 + t487) * m(7);
t342 = (t173 * t267 + t215 * t333 - t348 * t98) * t501 + t267 * t201 / 0.2e1 - t348 * t493;
t22 = (t387 + t484) * mrSges(7,2) + t342 + t353;
t44 = (m(7) * t173 + t201) * t348;
t345 = -qJD(1) * t22 - qJD(2) * t149 + qJD(3) * t44;
t202 = t444 + t447;
t148 = (-t389 + t486) * m(7);
t80 = m(7) * t215 - t445;
t79 = (-t411 - t422) * t501 + m(7) * t481;
t75 = t367 * t501 + t399;
t36 = t37 * qJD(4);
t34 = t400 * t402 + t528 * t476 + t343 + t400;
t29 = 0.2e1 * t501 * t63 + t216;
t21 = mrSges(7,2) * t387 + t451 / 0.2e1 + t342 - t353;
t13 = t341 + (-t470 + t519) * t474;
t9 = t452 / 0.2e1 - t456 / 0.2e1 - t450 / 0.2e1 - t455 / 0.2e1 + t335 * t400 - mrSges(5,2) * t415 / 0.2e1 + mrSges(5,1) * t388 + t337 - t354;
t8 = t223 * t484 + (-t152 * t334 + t264 * t464 - t266 * t463) * t501 + t222 * t485 + t340 + t518 * t489 + t520 * t474 + t530 * (-t425 / 0.2e1 - t171 / 0.2e1);
t6 = -t338 + t521 * t473 + t302 * t414 / 0.2e1 - t304 * t417 / 0.2e1 + t301 * t415 / 0.2e1 + t303 * t388 + t221 * t488 + t220 * t489 + t223 * t486 + t222 * t487 + t344 + t529 * t484 - t518 * t267 / 0.2e1 + (-t449 + t454 + t157) * t476 + (-t521 + t156 + t155) * t474;
t2 = Ifges(6,6) * t488 + Ifges(7,6) * t489 + t473 * t525 + t485 * t526 + t336 - t339;
t16 = [qJD(2) * t15 + qJD(3) * t3 + qJD(4) * t14 + qJD(5) * t4 + qJD(6) * t27, t6 * qJD(3) + t8 * qJD(5) + t79 * qJD(6) + t36 + t429 + 0.2e1 * t513 * qJD(2) * (t264 * t295 - t266 * t348) t442 + t6 * qJD(2) + t9 * qJD(4) + t2 * qJD(5) + t21 * qJD(6) + ((Ifges(5,2) * t332 + t462) * t479 + (Ifges(5,1) * t331 + t461) * t478 - Ifges(4,5) + t432 * t335) * t401 + ((-pkin(3) * t324 + qJ(4) * t359) * t505 + (-t214 * t73 + t215 * t74 + t292 * t320) * t503 + (t173 * t97 + t214 * t66 + t215 * t65) * t501) * t507 + (t66 * t348 * mrSges(7,2) + (Ifges(5,5) * t331 + Ifges(5,6) * t332 - t295 * t531 + t526 * t348) * t473 + t206 * t488 + t204 * t489 + t260 * t477 + t261 * t479 + t140 * t481 + t142 * t482 - t65 * t445 - mrSges(4,2) * t410 + t173 * t155 + t97 * t201 + pkin(3) * t521 + t292 * t202 + t320 * t156 - Ifges(4,6) * t334 + (t210 + t208) * t485 + (t146 + t144) * t348 / 0.2e1 + t529 * t215 + (-t220 + t221) * t214 + (t416 - t419) * qJ(4) + (-t295 * t74 - t348 * t73) * mrSges(6,3) + t359 * mrSges(5,3)) * qJD(3), qJD(3) * t9 + t363, t433 + t8 * qJD(2) + t2 * qJD(3) + (t368 * mrSges(7,2) + t510 * t71 + t511 * t72 + t515) * qJD(5) + t29 * qJD(6), qJD(2) * t79 + qJD(3) * t21 + qJD(5) * t29 + t431; -qJD(3) * t5 - qJD(5) * t7 - qJD(6) * t78 + t36 - t429, qJD(3) * t31, t34 * qJD(4) + t13 * qJD(5) + t148 * qJD(6) + (t201 + t202 + t432) * t401 + ((t173 * t333 + t376) * t501 + (t320 * t333 + t376) * t503 + (t334 * t378 - t469) * t505) * t507 + t365 + ((-mrSges(4,2) + t380) * t334 + t530 * (t265 * t348 + t267 * t295)) * qJD(3), qJD(3) * t34 + t430, t13 * qJD(3) + (-t527 * t266 + (mrSges(6,2) - mrSges(7,3)) * t264) * qJD(5) + ((pkin(5) * t266 - qJ(6) * t264) * qJD(5) / 0.2e1 + qJD(6) * t484) * t508 - t364, -m(7) * t266 * qJD(5) + t148 * qJD(3) - t409; qJD(2) * t5 + qJD(4) * t10 + qJD(5) * t1 + qJD(6) * t22 - t442, qJD(4) * t35 - qJD(5) * t12 + qJD(6) * t149 - t365, qJD(4) * t24 + qJD(5) * t11 - qJD(6) * t44, qJD(5) * t75 + t346, t75 * qJD(4) + (t366 * mrSges(7,2) - t214 * t510 + t215 * t511 + t516) * qJD(5) + t80 * qJD(6) + t347, qJD(5) * t80 - t345; -qJD(3) * t10 - qJD(5) * t33 + qJD(6) * t117 - t363, -qJD(3) * t35 - t430, -qJD(5) * t42 - qJD(6) * t175 - t346, 0, -t362, -t355; qJD(2) * t7 - qJD(3) * t1 + qJD(4) * t33 + qJD(6) * t30 - t433, qJD(3) * t12 + t364, qJD(4) * t42 - t347, t362, t328 * qJD(6), t356; qJD(2) * t78 - qJD(3) * t22 - qJD(4) * t117 - qJD(5) * t30 - t431, -qJD(3) * t149 + t409, qJD(4) * t175 + t345, t355, -t356, 0;];
Cq  = t16;
