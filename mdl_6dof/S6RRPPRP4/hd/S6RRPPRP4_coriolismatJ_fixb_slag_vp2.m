% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:05
% EndTime: 2019-03-09 08:37:19
% DurationCPUTime: 7.00s
% Computational Cost: add. (10674->616), mult. (22838->819), div. (0->0), fcn. (21911->6), ass. (0->296)
t550 = mrSges(7,2) + mrSges(6,3);
t545 = Ifges(7,4) + Ifges(6,5);
t552 = Ifges(6,6) - Ifges(7,6);
t551 = -m(7) - m(6);
t362 = cos(qJ(2));
t361 = sin(qJ(2));
t317 = -pkin(2) * t362 - t361 * qJ(3) - pkin(1);
t358 = sin(pkin(9));
t440 = t358 * t362;
t336 = pkin(7) * t440;
t357 = t362 * pkin(3);
t359 = cos(pkin(9));
t162 = pkin(4) * t362 + t336 + t357 + (-pkin(8) * t361 - t317) * t359;
t437 = t359 * t362;
t223 = pkin(7) * t437 + t358 * t317;
t203 = -qJ(4) * t362 + t223;
t441 = t358 * t361;
t175 = pkin(8) * t441 + t203;
t360 = sin(qJ(5));
t490 = cos(qJ(5));
t68 = t162 * t490 - t360 * t175;
t63 = -t362 * pkin(5) - t68;
t549 = t63 + t68;
t304 = t360 * t358 + t359 * t490;
t271 = t362 * t304;
t464 = t271 * mrSges(7,2);
t484 = mrSges(6,1) + mrSges(7,1);
t483 = mrSges(6,2) - mrSges(7,3);
t546 = Ifges(5,4) + Ifges(4,5);
t544 = Ifges(7,2) + Ifges(6,3);
t543 = Ifges(4,6) - Ifges(5,6);
t542 = -t464 / 0.2e1;
t412 = t490 * t358;
t438 = t359 * t361;
t268 = t360 * t438 - t361 * t412;
t472 = t268 * mrSges(6,3);
t269 = t304 * t361;
t468 = t269 * mrSges(6,3);
t434 = t362 * qJ(6);
t415 = t490 * t175;
t436 = t360 * t162;
t69 = t415 + t436;
t62 = t69 + t434;
t316 = -t359 * pkin(3) - t358 * qJ(4) - pkin(2);
t290 = t359 * pkin(4) - t316;
t305 = t359 * t360 - t412;
t391 = pkin(5) * t304 + qJ(6) * t305;
t131 = t391 + t290;
t180 = mrSges(7,1) * t304 + mrSges(7,3) * t305;
t400 = m(7) * t131 + t180;
t474 = Ifges(7,5) * t304;
t184 = -Ifges(7,1) * t305 + t474;
t301 = Ifges(6,4) * t304;
t185 = -Ifges(6,1) * t305 - t301;
t541 = t185 + t184;
t435 = t360 * t362;
t270 = t359 * t435 - t362 * t412;
t206 = -mrSges(7,2) * t270 - mrSges(7,3) * t361;
t211 = mrSges(6,2) * t361 - mrSges(6,3) * t270;
t540 = t206 + t211;
t454 = t362 * mrSges(7,3);
t207 = -t268 * mrSges(7,2) + t454;
t208 = -mrSges(6,2) * t362 - t472;
t539 = t207 + t208;
t209 = mrSges(6,1) * t362 - t468;
t210 = -mrSges(7,1) * t362 + t269 * mrSges(7,2);
t538 = t209 - t210;
t535 = -t545 * t304 + t552 * t305;
t534 = -t545 * t268 - t552 * t269;
t452 = qJ(3) * t362;
t320 = t361 * pkin(2) - t452;
t439 = t359 * t320;
t247 = pkin(7) * t441 + t439;
t248 = -pkin(7) * t438 + t358 * t320;
t533 = -t247 * t358 + t248 * t359;
t353 = t361 * qJ(4);
t205 = t248 + t353;
t418 = -pkin(7) * t358 - pkin(3);
t214 = t361 * t418 - t439;
t532 = t205 * t359 + t214 * t358;
t514 = m(7) / 0.2e1;
t516 = m(6) / 0.2e1;
t531 = t516 + t514;
t530 = mrSges(4,2) / 0.2e1 - mrSges(5,3) / 0.2e1;
t352 = m(7) * qJ(6) + mrSges(7,3);
t145 = mrSges(7,1) * t268 - mrSges(7,3) * t269;
t333 = t359 * t353;
t528 = pkin(7) + (pkin(3) + pkin(4)) * t358;
t201 = -t361 * t528 + t333;
t392 = pkin(5) * t268 - qJ(6) * t269;
t78 = t201 + t392;
t529 = -m(7) * t78 - t145;
t527 = -m(7) * pkin(5) - t484;
t526 = -mrSges(6,2) + t352;
t421 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t422 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t525 = t421 * t270 + t422 * t271;
t524 = 0.2e1 * m(7);
t519 = m(4) / 0.2e1;
t518 = m(5) / 0.2e1;
t517 = -m(6) / 0.2e1;
t515 = -m(7) / 0.2e1;
t513 = mrSges(7,1) / 0.2e1;
t512 = mrSges(6,2) / 0.2e1;
t511 = -mrSges(7,3) / 0.2e1;
t482 = -pkin(8) + qJ(3);
t319 = t482 * t359;
t401 = t482 * t358;
t199 = t319 * t360 - t401 * t490;
t509 = -t199 / 0.2e1;
t508 = -t201 / 0.2e1;
t507 = -t269 / 0.2e1;
t506 = -t270 / 0.2e1;
t505 = t270 / 0.2e1;
t504 = t271 / 0.2e1;
t503 = -t304 / 0.2e1;
t502 = t304 / 0.2e1;
t501 = -t305 / 0.2e1;
t500 = t305 / 0.2e1;
t499 = -t358 / 0.2e1;
t498 = t358 / 0.2e1;
t497 = t359 / 0.2e1;
t496 = -t360 / 0.2e1;
t495 = t360 / 0.2e1;
t494 = -t361 / 0.2e1;
t488 = m(7) * t360;
t487 = t269 * pkin(5);
t486 = t305 * pkin(5);
t485 = t362 * pkin(7);
t481 = Ifges(4,4) * t358;
t480 = Ifges(4,4) * t359;
t479 = Ifges(6,4) * t269;
t478 = Ifges(6,4) * t305;
t477 = Ifges(5,5) * t358;
t476 = Ifges(5,5) * t359;
t475 = Ifges(7,5) * t268;
t473 = t268 * mrSges(6,2);
t471 = t268 * mrSges(7,3);
t470 = t269 * mrSges(6,1);
t469 = t269 * mrSges(7,1);
t467 = t270 * mrSges(6,1);
t466 = t270 * mrSges(7,1);
t465 = t271 * mrSges(6,2);
t463 = t271 * mrSges(7,3);
t133 = Ifges(7,5) * t271 - Ifges(7,6) * t361 + Ifges(7,3) * t270;
t135 = Ifges(6,4) * t271 - Ifges(6,2) * t270 - Ifges(6,6) * t361;
t137 = Ifges(7,1) * t271 - Ifges(7,4) * t361 + Ifges(7,5) * t270;
t139 = Ifges(6,1) * t271 - Ifges(6,4) * t270 - Ifges(6,5) * t361;
t146 = mrSges(6,1) * t268 + mrSges(6,2) * t269;
t147 = -t463 + t466;
t148 = t465 + t467;
t334 = qJ(4) * t437;
t202 = t362 * t528 - t334;
t222 = t359 * t317 - t336;
t204 = -t222 + t357;
t212 = -mrSges(6,1) * t361 - mrSges(6,3) * t271;
t455 = t361 * mrSges(7,1);
t213 = t455 + t464;
t256 = Ifges(5,6) * t361 + (Ifges(5,3) * t358 + t476) * t362;
t257 = Ifges(4,6) * t361 + (-Ifges(4,2) * t358 + t480) * t362;
t258 = Ifges(5,4) * t361 + (Ifges(5,1) * t359 + t477) * t362;
t259 = Ifges(4,5) * t361 + (Ifges(4,1) * t359 - t481) * t362;
t417 = pkin(3) * t358 + pkin(7);
t260 = t361 * t417 - t333;
t261 = t362 * t417 - t334;
t394 = t358 * mrSges(5,1) - t359 * mrSges(5,3);
t291 = t394 * t361;
t292 = t394 * t362;
t293 = (t358 * mrSges(4,1) + t359 * mrSges(4,2)) * t362;
t308 = t362 * mrSges(4,2) - mrSges(4,3) * t441;
t309 = -t361 * mrSges(4,2) - mrSges(4,3) * t440;
t310 = -mrSges(4,1) * t362 - mrSges(4,3) * t438;
t311 = mrSges(5,1) * t362 + mrSges(5,2) * t438;
t312 = t361 * mrSges(4,1) - mrSges(4,3) * t437;
t335 = mrSges(5,2) * t437;
t456 = t361 * mrSges(5,1);
t313 = t335 - t456;
t314 = -mrSges(5,2) * t440 + t361 * mrSges(5,3);
t315 = -mrSges(5,2) * t441 - t362 * mrSges(5,3);
t136 = Ifges(7,1) * t269 + t362 * Ifges(7,4) + t475;
t255 = Ifges(6,4) * t268;
t138 = Ifges(6,1) * t269 + t362 * Ifges(6,5) - t255;
t405 = t136 / 0.2e1 + t138 / 0.2e1;
t252 = Ifges(7,5) * t269;
t132 = t362 * Ifges(7,6) + Ifges(7,3) * t268 + t252;
t134 = -Ifges(6,2) * t268 + t362 * Ifges(6,6) + t479;
t406 = t132 / 0.2e1 - t134 / 0.2e1;
t163 = (-pkin(8) * t362 - t320) * t359 + (-pkin(4) + t418) * t361;
t176 = pkin(8) * t440 + t205;
t72 = t360 * t163 + t490 * t176;
t64 = -qJ(6) * t361 + t72;
t71 = t163 * t490 - t360 * t176;
t65 = t361 * pkin(5) - t71;
t79 = -t270 * pkin(5) + t271 * qJ(6) + t202;
t3 = t214 * t311 + t222 * t312 + t204 * t313 + t203 * t314 + t205 * t315 + t248 * t308 + t223 * t309 + t247 * t310 + t261 * t291 + t260 * t292 + t68 * t212 + t63 * t213 + t62 * t206 + t64 * t207 + t72 * t208 + t71 * t209 + t65 * t210 + t69 * t211 + t201 * t148 - t202 * t146 + (t137 / 0.2e1 + t139 / 0.2e1) * t269 + (-t135 / 0.2e1 + t133 / 0.2e1) * t268 - t79 * t145 + t78 * t147 + m(4) * (t222 * t247 + t223 * t248) + m(5) * (t203 * t205 + t204 * t214 + t260 * t261) + m(6) * (-t201 * t202 + t68 * t71 + t69 * t72) + m(7) * (t62 * t64 + t63 * t65 - t78 * t79) + (-pkin(1) * mrSges(3,2) + (-Ifges(5,2) - Ifges(4,3) + m(4) * pkin(7) ^ 2 + Ifges(3,1) - Ifges(3,2) + (pkin(7) * mrSges(4,2) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t359) * t359 + (pkin(7) * mrSges(4,1) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t358 + (-Ifges(4,4) + Ifges(5,5)) * t359) * t358 - t544) * t361 + (t358 * t543 - t359 * t546 + Ifges(3,4)) * t362 + t525) * t362 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t361 + pkin(7) * t293 - t422 * t269 - t421 * t268 + (t259 / 0.2e1 + t258 / 0.2e1 + (Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1) * t361) * t359 + (t256 / 0.2e1 - t257 / 0.2e1 + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t361) * t358) * t361 + t405 * t271 + t406 * t270;
t462 = t3 * qJD(1);
t461 = t304 * mrSges(6,2);
t460 = t304 * mrSges(7,2);
t459 = t304 * mrSges(7,3);
t458 = t305 * mrSges(6,1);
t457 = t305 * mrSges(7,1);
t143 = t469 + t471;
t144 = t470 - t473;
t149 = Ifges(7,3) * t269 - t475;
t150 = -Ifges(6,2) * t269 - t255;
t151 = -Ifges(7,1) * t268 + t252;
t152 = -Ifges(6,1) * t268 - t479;
t451 = qJ(6) * t268;
t393 = t451 + t487;
t4 = t201 * t144 + t78 * t143 + (t152 / 0.2e1 + t151 / 0.2e1 - t62 * mrSges(7,2) + t406) * t269 + (-t63 * mrSges(7,2) - t150 / 0.2e1 + t149 / 0.2e1 - t405) * t268 - t529 * t393 + (m(7) * t63 - t468 - t538) * t69 + (m(7) * t62 + t472 + t539) * t68 + t534 * t362 / 0.2e1;
t453 = t4 * qJD(1);
t450 = qJ(6) * t304;
t10 = t538 * t269 + t539 * t268 + m(7) * (t268 * t62 - t269 * t63) + m(6) * (t268 * t69 + t269 * t68) + ((-t310 + t311) * t359 + (-t308 - t315) * t358 + m(5) * (-t203 * t358 + t204 * t359) + m(4) * (-t222 * t359 - t223 * t358)) * t361;
t449 = qJD(1) * t10;
t413 = t490 * t208;
t414 = t490 * t207;
t419 = t490 * t69;
t420 = t490 * t62;
t11 = t538 * t435 + (-m(5) * t260 + m(6) * t201 + t146 - t291 - t529) * t438 + (-t413 + m(7) * (-t360 * t63 - t420) - t414 + m(6) * (t360 * t68 - t419) - t315 - m(5) * t203) * t362;
t448 = t11 * qJD(1);
t427 = t490 / 0.2e1;
t366 = (t549 * t360 - t419 + t420) * t514 + t209 * t496 + t210 * t495 + t414 / 0.2e1 + t413 / 0.2e1 + t550 * (t268 * t427 + t269 * t496);
t378 = -pkin(5) * t360 + qJ(6) * t490;
t416 = t362 * t490;
t396 = t416 / 0.2e1;
t408 = t435 / 0.2e1;
t368 = t362 * t378 * t515 + mrSges(6,2) * t396 + t408 * t484 + t416 * t511;
t12 = -t366 + t368;
t447 = t12 * qJD(1);
t22 = m(7) * (-t78 * t269 + t362 * t62) - t269 * t145 + t362 * t207;
t444 = t22 * qJD(1);
t429 = t515 + t517;
t428 = m(7) * t435;
t426 = qJD(6) * t488;
t425 = t513 + mrSges(6,1) / 0.2e1;
t424 = t512 + t511;
t423 = mrSges(7,2) / 0.2e1 + mrSges(6,3) / 0.2e1;
t404 = -t207 / 0.2e1 - t208 / 0.2e1;
t403 = -t209 / 0.2e1 + t210 / 0.2e1;
t200 = t319 * t490 + t360 * t401;
t399 = -t199 * t269 + t200 * t268;
t179 = -t450 + t486;
t364 = (t200 * t507 + t268 * t509 + (-t69 / 0.2e1 + t62 / 0.2e1) * t305 + (-t63 / 0.2e1 - t68 / 0.2e1) * t304) * mrSges(7,2) + (-t152 / 0.4e1 - t151 / 0.4e1 - t132 / 0.4e1 - t78 * mrSges(7,1) / 0.2e1 + mrSges(6,1) * t508 + t134 / 0.4e1) * t305 + (t78 * mrSges(7,3) / 0.2e1 + mrSges(6,2) * t508 + t149 / 0.4e1 - t150 / 0.4e1 - t138 / 0.4e1 - t136 / 0.4e1) * t304 + (t131 * t393 - t179 * t78 + (-t62 + t69) * t199 + t549 * t200) * t514 + t131 * t143 / 0.2e1 + t393 * t180 / 0.2e1 - t179 * t145 / 0.2e1 + t290 * t144 / 0.2e1 + t535 * t362 / 0.4e1;
t369 = (-pkin(5) * t65 + qJ(6) * t64) * t515 + pkin(5) * t213 / 0.2e1 - qJ(6) * t206 / 0.2e1 + t64 * t511 + t65 * t513 - t71 * mrSges(6,1) / 0.2e1 + t72 * t512;
t373 = -t474 / 0.4e1 + t301 / 0.4e1 - t185 / 0.4e1 - t184 / 0.4e1 + (-Ifges(7,3) / 0.4e1 - Ifges(6,2) / 0.4e1) * t305;
t298 = Ifges(7,5) * t305;
t182 = Ifges(7,3) * t304 - t298;
t183 = -Ifges(6,2) * t304 - t478;
t374 = t478 / 0.4e1 - t298 / 0.4e1 - t183 / 0.4e1 + t182 / 0.4e1 + (-Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1) * t304;
t2 = (-t472 / 0.2e1 + t404) * t199 + t364 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t361 + t374 * t269 + t373 * t268 + (-t468 / 0.2e1 + t403) * t200 + t369 - t525;
t9 = -t131 * (-t457 + t459) - t290 * (-t458 - t461) + (-Ifges(7,3) * t305 - t474) * t503 + t183 * t501 + t400 * t179 + (Ifges(6,2) * t305 - t301 + t541) * t502 + (t182 - t298 + t478 + (-Ifges(6,1) - Ifges(7,1)) * t304) * t500;
t390 = t2 * qJD(1) - t9 * qJD(2);
t19 = t551 * (t199 * t305 + t200 * t304) + t550 * (t304 ^ 2 + t305 ^ 2) + (-0.4e1 * (-m(5) / 0.4e1 - m(4) / 0.4e1) * qJ(3) + mrSges(4,3) + mrSges(5,2)) * (-t358 ^ 2 - t359 ^ 2);
t365 = (t311 / 0.2e1 - t310 / 0.2e1) * t358 + (t315 / 0.2e1 + t308 / 0.2e1) * t359 + (-t268 * t423 - t404) * t304 + (t269 * t423 + t403) * t305 + (-t222 * t358 + t223 * t359) * t519 + (t203 * t359 + t204 * t358) * t518 + (t304 * t69 - t305 * t68 + t399) * t516 + (t304 * t62 + t305 * t63 + t399) * t514;
t375 = -m(5) * t261 / 0.2e1 + t202 * t517 + t79 * t515;
t6 = t424 * t271 + t425 * t270 + (-m(4) * pkin(7) / 0.2e1 - t530 * t359 + (-mrSges(4,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t358) * t362 + t365 + t375;
t389 = qJD(1) * t6 - qJD(2) * t19;
t181 = mrSges(6,1) * t304 - mrSges(6,2) * t305;
t318 = -mrSges(5,1) * t359 - mrSges(5,3) * t358;
t36 = (-m(5) * t316 + m(6) * t290 + t181 - t318 + t400) * t358;
t376 = (-t199 * t360 - t200 * t490) * t362;
t363 = (-t358 * t260 + (-t316 * t361 - t452) * t359) * t518 + (t358 * t201 + t290 * t438 + t376) * t516 + (t131 * t438 + t358 * t78 + t376) * t514 + t291 * t499 - t318 * t438 / 0.2e1 + (t145 + t146) * t498 + (t180 + t181) * t438 / 0.2e1 + t550 * (t304 * t396 + t305 * t408);
t367 = t214 * t518 + (t360 * t72 + t490 * t71) * t516 + (t360 * t64 - t490 * t65) * t514 - t456 / 0.2e1 + t212 * t427 - t490 * t213 / 0.2e1 + t540 * t495;
t7 = t363 - t335 - t367;
t388 = -qJD(1) * t7 - qJD(2) * t36;
t387 = t429 * t438;
t371 = (-t131 * t269 + t200 * t362 + t305 * t78) * t514 + t180 * t507 + t145 * t500;
t379 = t65 * t515 - t455 / 0.2e1;
t17 = t371 + t379 + 0.2e1 * t542;
t44 = t400 * t305;
t386 = -qJD(1) * t17 - qJD(2) * t44;
t24 = -t484 * t269 + t483 * t268 + (-t393 / 0.4e1 - t487 / 0.4e1 - t451 / 0.4e1) * t524;
t30 = t484 * t305 + t483 * t304 + (t179 / 0.4e1 + t486 / 0.4e1 - t450 / 0.4e1) * t524;
t385 = qJD(1) * t24 + qJD(2) * t30;
t370 = t531 * (t360 * t268 + t269 * t490);
t40 = -m(5) * t438 - t370 + t387;
t372 = t531 * (t360 * t304 - t305 * t490);
t52 = (m(5) - t429) * t358 + t372;
t384 = qJD(1) * t40 - qJD(2) * t52;
t26 = t454 + (t415 / 0.4e1 + t434 / 0.2e1 + t436 / 0.4e1 - t69 / 0.4e1) * t524;
t381 = qJD(1) * t26 + qJD(5) * t352;
t120 = m(7) * t269;
t165 = m(7) * t305;
t380 = qJD(1) * t120 - qJD(2) * t165;
t81 = m(7) * t200 - t460;
t57 = -t551 * t499 + t372;
t42 = t387 + t370;
t29 = t459 / 0.2e1 - t461 / 0.2e1 - t458 / 0.2e1 - t457 / 0.2e1 + t425 * t305 + t424 * t304;
t25 = 0.2e1 * t514 * t62 + t207;
t23 = t470 / 0.2e1 - t473 / 0.2e1 + t471 / 0.2e1 + t469 / 0.2e1 - t425 * t269 + t424 * t268;
t16 = t542 + t464 / 0.2e1 + t371 - t379;
t13 = t366 + t368;
t8 = t363 + t367;
t5 = t365 - t466 / 0.2e1 - t467 / 0.2e1 - t465 / 0.2e1 + t463 / 0.2e1 + t485 * t519 - t375 + t530 * t437 + (mrSges(4,1) + mrSges(5,1)) * t440 / 0.2e1;
t1 = t364 + t403 * t200 + t404 * t199 + Ifges(7,6) * t505 + Ifges(6,6) * t506 + (-t200 * mrSges(6,3) / 0.2e1 + t374) * t269 + (mrSges(6,3) * t509 + t373) * t268 - t369 + t545 * t504 + t544 * t494;
t14 = [qJD(2) * t3 + qJD(3) * t10 + qJD(4) * t11 + qJD(5) * t4 + qJD(6) * t22, t462 + ((t139 + t137) * t501 + (t259 + t258) * t498 + (pkin(7) * mrSges(3,2) - Ifges(3,6)) * t361 + t182 * t505 + t183 * t506 + t257 * t497 + t133 * t502 + t135 * t503 - t359 * t256 / 0.2e1 + t316 * t292 + t261 * t318 + t290 * t148 - pkin(2) * t293 - t202 * t181 - t79 * t180 + t131 * t147 + (-t212 + t213) * t199 + (-t552 * t304 - t545 * t305) * t494 - t65 * t305 * mrSges(7,2) + t540 * t200 + t541 * t504 + (t358 * t546 + t359 * t543) * t361 / 0.2e1 + t532 * mrSges(5,2) + t533 * mrSges(4,3) - t64 * t460 + (-t304 * t72 + t305 * t71) * mrSges(6,3)) * qJD(2) + t5 * qJD(3) + t8 * qJD(4) + t1 * qJD(5) + t16 * qJD(6) + ((t309 + t314) * t359 + (-t312 + t313) * t358) * qJD(2) * qJ(3) + (Ifges(3,5) + (-Ifges(5,3) * t359 + t477) * t498 + (Ifges(4,2) * t359 + t481) * t499 + (-mrSges(4,1) * t359 + mrSges(4,2) * t358 - mrSges(3,1)) * pkin(7) + (-t476 + t480 + (Ifges(4,1) + Ifges(5,1)) * t358) * t497) * qJD(2) * t362 + 0.2e1 * ((-pkin(2) * t485 + qJ(3) * t533) * t519 + (qJ(3) * t532 + t261 * t316) * t518 + (-t199 * t71 + t200 * t72 - t202 * t290) * t516 + (-t131 * t79 + t199 * t65 + t200 * t64) * t514) * qJD(2), qJD(2) * t5 + qJD(4) * t42 + qJD(5) * t23 + t449, qJD(2) * t8 + qJD(3) * t42 + qJD(5) * t13 + t448, t453 + t1 * qJD(2) + t23 * qJD(3) + t13 * qJD(4) + (t392 * mrSges(7,2) + t526 * t68 + t527 * t69 + t534) * qJD(5) + t25 * qJD(6), qJD(2) * t16 + qJD(5) * t25 + t444; qJD(3) * t6 + qJD(4) * t7 + qJD(5) * t2 + qJD(6) * t17 - t462, -qJD(3) * t19 + qJD(4) * t36 - qJD(5) * t9 + qJD(6) * t44, qJD(4) * t57 + qJD(5) * t29 + t389, qJD(3) * t57 - t388, t29 * qJD(3) + (t391 * mrSges(7,2) - t199 * t526 + t200 * t527 + t535) * qJD(5) + t81 * qJD(6) + t390, qJD(5) * t81 - t386; -qJD(2) * t6 + qJD(4) * t40 + qJD(5) * t24 + qJD(6) * t120 - t449, -qJD(4) * t52 + qJD(5) * t30 - qJD(6) * t165 - t389, 0, t384, t385, t380; -t7 * qJD(2) - t40 * qJD(3) - t12 * qJD(5) + t362 * t426 - t448, qJD(3) * t52 + t388, -t384, 0, -t447 + (m(7) * t378 - t360 * t484 - t483 * t490) * qJD(5) + t426 (qJD(1) * t362 + qJD(5)) * t488; -qJD(2) * t2 - qJD(3) * t24 + qJD(4) * t12 + qJD(6) * t26 - t453, -qJD(3) * t30 - t390, -t385, t447, t352 * qJD(6), t381; -t17 * qJD(2) - t120 * qJD(3) - qJD(4) * t428 - t26 * qJD(5) - t444, qJD(3) * t165 + t386, -t380, -qJD(1) * t428, -t381, 0;];
Cq  = t14;
