% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:39
% EndTime: 2019-03-09 08:41:48
% DurationCPUTime: 5.66s
% Computational Cost: add. (11954->550), mult. (23348->718), div. (0->0), fcn. (22832->6), ass. (0->272)
t320 = sin(pkin(9));
t321 = cos(pkin(9));
t446 = sin(qJ(5));
t447 = cos(qJ(5));
t280 = t447 * t320 + t446 * t321;
t322 = sin(qJ(2));
t406 = t280 * t322;
t505 = t406 / 0.2e1;
t501 = Ifges(7,4) + Ifges(6,5);
t504 = Ifges(6,6) - Ifges(7,6);
t472 = m(7) / 0.2e1;
t503 = 0.2e1 * t472;
t466 = pkin(3) + pkin(7);
t502 = mrSges(7,2) + mrSges(6,3);
t323 = cos(qJ(2));
t414 = qJ(3) * t322;
t441 = pkin(2) + qJ(4);
t275 = -t323 * t441 - pkin(1) - t414;
t296 = t466 * t322;
t286 = t321 * t296;
t129 = pkin(4) * t322 + t286 + (pkin(8) * t323 - t275) * t320;
t180 = t321 * t275 + t320 * t296;
t402 = t321 * t323;
t156 = -pkin(8) * t402 + t180;
t64 = t129 * t447 - t156 * t446;
t58 = -t322 * pkin(5) - t64;
t439 = t58 + t64;
t279 = t446 * t320 - t447 * t321;
t188 = t280 * mrSges(7,1) + t279 * mrSges(7,3);
t307 = t320 * pkin(4) + qJ(3);
t353 = pkin(5) * t280 + qJ(6) * t279;
t158 = t353 + t307;
t445 = m(7) * t158;
t500 = t188 + t445;
t247 = t279 * t322;
t461 = -t247 / 0.2e1;
t460 = t247 / 0.2e1;
t499 = -t406 / 0.2e1;
t498 = -t323 / 0.2e1;
t497 = t323 / 0.2e1;
t496 = mrSges(6,1) + mrSges(7,1);
t495 = -Ifges(3,4) - Ifges(4,6);
t494 = Ifges(7,2) + Ifges(6,3);
t248 = t279 * t323;
t493 = t248 * mrSges(6,3);
t250 = t280 * t323;
t492 = t250 * mrSges(6,3);
t401 = t322 * qJ(6);
t383 = t446 * t129;
t386 = t447 * t156;
t65 = t386 + t383;
t57 = t65 + t401;
t468 = mrSges(6,3) / 0.2e1;
t389 = mrSges(7,2) / 0.2e1 + t468;
t491 = t250 * t389;
t419 = t322 * mrSges(7,3);
t204 = t248 * mrSges(7,2) + t419;
t207 = -mrSges(6,2) * t322 + t493;
t400 = t204 + t207;
t210 = mrSges(6,1) * t322 + t492;
t211 = -mrSges(7,1) * t322 - t250 * mrSges(7,2);
t399 = t210 - t211;
t319 = t321 ^ 2;
t394 = t320 ^ 2 + t319;
t488 = t279 * t504 - t280 * t501;
t487 = t248 * t501 + t250 * t504;
t366 = t322 * pkin(2) - qJ(3) * t323;
t276 = qJ(4) * t322 + t366;
t297 = t466 * t323;
t287 = t321 * t297;
t181 = -t276 * t320 + t287;
t182 = t321 * t276 + t320 * t297;
t486 = t181 * t321 + t182 * t320;
t474 = m(6) / 0.2e1;
t485 = t474 + t472;
t484 = t279 ^ 2 + t280 ^ 2;
t312 = m(7) * qJ(6) + mrSges(7,3);
t483 = -m(7) * pkin(5) - t496;
t482 = -mrSges(6,2) + t312;
t387 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t388 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t481 = t387 * t247 + t388 * t406;
t480 = 0.2e1 * m(7);
t477 = -m(5) / 0.2e1;
t476 = m(5) / 0.2e1;
t475 = -m(6) / 0.2e1;
t473 = -m(7) / 0.2e1;
t471 = -mrSges(6,1) / 0.2e1;
t470 = mrSges(7,1) / 0.2e1;
t469 = mrSges(6,2) / 0.2e1;
t467 = -mrSges(7,3) / 0.2e1;
t413 = qJ(6) * t248;
t444 = pkin(5) * t250;
t355 = t413 + t444;
t465 = -t355 / 0.2e1;
t392 = -pkin(8) - t441;
t292 = t392 * t320;
t360 = t321 * t392;
t183 = t292 * t446 - t360 * t447;
t464 = t183 / 0.2e1;
t412 = qJ(6) * t280;
t443 = pkin(5) * t279;
t352 = -t412 + t443;
t463 = -t352 / 0.2e1;
t417 = t323 * mrSges(7,1);
t209 = mrSges(7,2) * t406 - t417;
t462 = t209 / 0.2e1;
t458 = t250 / 0.2e1;
t456 = t279 / 0.2e1;
t455 = t280 / 0.2e1;
t454 = -t320 / 0.2e1;
t453 = t320 / 0.2e1;
t452 = -t321 / 0.2e1;
t451 = t321 / 0.2e1;
t440 = -t57 + t65;
t438 = Ifges(5,1) * t320;
t437 = Ifges(5,4) * t320;
t436 = Ifges(5,4) * t321;
t435 = Ifges(6,4) * t250;
t434 = Ifges(6,4) * t279;
t433 = Ifges(7,5) * t248;
t432 = Ifges(7,5) * t280;
t431 = t247 * mrSges(6,1);
t430 = t247 * mrSges(7,1);
t429 = t248 * mrSges(7,1);
t428 = t406 * mrSges(6,2);
t426 = t406 * mrSges(7,3);
t425 = t250 * mrSges(7,3);
t137 = -t426 + t430;
t138 = t428 + t431;
t139 = t425 - t429;
t179 = -t275 * t320 + t286;
t205 = -t247 * mrSges(7,2) + mrSges(7,3) * t323;
t206 = -mrSges(6,2) * t323 - t247 * mrSges(6,3);
t208 = mrSges(6,1) * t323 - mrSges(6,3) * t406;
t243 = Ifges(5,6) * t323 + (Ifges(5,2) * t321 + t437) * t322;
t244 = Ifges(5,5) * t323 + (t436 + t438) * t322;
t261 = (-pkin(4) * t321 - t466) * t322;
t262 = pkin(4) * t402 + t297;
t357 = t321 * mrSges(5,1) - t320 * mrSges(5,2);
t263 = t357 * t322;
t405 = t320 * t322;
t418 = t323 * mrSges(5,1);
t288 = -mrSges(5,3) * t405 + t418;
t404 = t320 * t323;
t421 = t322 * mrSges(5,1);
t289 = mrSges(5,3) * t404 + t421;
t403 = t321 * t322;
t416 = t323 * mrSges(5,2);
t290 = mrSges(5,3) * t403 - t416;
t420 = t322 * mrSges(5,2);
t291 = -mrSges(5,3) * t402 - t420;
t356 = -pkin(2) * t323 - t414;
t293 = -pkin(1) + t356;
t340 = -t261 * mrSges(6,2) + Ifges(6,4) * t460 + Ifges(7,5) * t461 + (Ifges(6,1) + Ifges(7,1)) * t499 + t501 * t498;
t294 = t323 * mrSges(4,2) - t322 * mrSges(4,3);
t368 = m(4) * t293 + t294;
t125 = -Ifges(7,1) * t250 + t322 * Ifges(7,4) - t433;
t242 = Ifges(6,4) * t248;
t127 = -Ifges(6,1) * t250 + t322 * Ifges(6,5) + t242;
t373 = -t125 / 0.2e1 - t127 / 0.2e1;
t239 = Ifges(7,5) * t250;
t121 = t322 * Ifges(7,6) - Ifges(7,3) * t248 - t239;
t123 = Ifges(6,2) * t248 + t322 * Ifges(6,6) - t435;
t374 = -t123 / 0.2e1 + t121 / 0.2e1;
t375 = Ifges(6,4) * t505 + Ifges(7,5) * t499 + Ifges(6,6) * t497 + Ifges(7,6) * t498 + (Ifges(6,2) + Ifges(7,3)) * t461;
t422 = t321 * Ifges(5,6);
t423 = t320 * Ifges(5,5);
t130 = pkin(4) * t323 + t287 + (-pkin(8) * t322 - t276) * t320;
t157 = pkin(8) * t403 + t182;
t68 = t446 * t130 + t447 * t157;
t59 = qJ(6) * t323 + t68;
t67 = t130 * t447 - t157 * t446;
t60 = -t323 * pkin(5) - t67;
t90 = pkin(5) * t247 - qJ(6) * t406 + t261;
t354 = -pkin(5) * t248 + qJ(6) * t250;
t91 = t354 + t262;
t3 = -t373 * t406 + t374 * t247 + (-t261 * mrSges(6,1) + t375) * t248 + t368 * t366 + m(5) * (t179 * t181 + t180 * t182 - t296 * t297) + m(6) * (t261 * t262 + t64 * t67 + t65 * t68) + m(7) * (t57 * t59 + t58 * t60 + t90 * t91) + t340 * t250 + (-pkin(1) * mrSges(3,2) - t293 * mrSges(4,3) + t243 * t452 + t244 * t454 + t387 * t248 + t388 * t250 - t296 * t357 + (-t423 / 0.2e1 - t422 / 0.2e1 - t495) * t323) * t323 + (-pkin(1) * mrSges(3,1) - t293 * mrSges(4,2) + (t422 + t423 + t495) * t322 + (-t319 * Ifges(5,2) / 0.2e1 + Ifges(4,2) - Ifges(4,3) + Ifges(3,1) + Ifges(5,3) - Ifges(3,2) + (-t436 - t438 / 0.2e1) * t320 + t494) * t323 - t481) * t322 - t297 * t263 + t179 * t288 + t181 * t289 + t180 * t290 + t182 * t291 + t262 * t138 + t58 * t209 + t67 * t210 + t60 * t211 + t59 * t204 + t57 * t205 + t65 * t206 + t68 * t207 + t64 * t208 + t91 * t137 + t90 * t139;
t424 = t3 * qJD(1);
t135 = -t250 * mrSges(7,1) - t248 * mrSges(7,3);
t136 = -t250 * mrSges(6,1) + t248 * mrSges(6,2);
t140 = -Ifges(7,3) * t250 + t433;
t141 = Ifges(6,2) * t250 + t242;
t142 = Ifges(7,1) * t248 - t239;
t143 = Ifges(6,1) * t248 + t435;
t4 = t262 * t136 + t91 * t135 + (-t142 / 0.2e1 - t143 / 0.2e1 + t57 * mrSges(7,2) - t374) * t250 + (t58 * mrSges(7,2) - t140 / 0.2e1 + t141 / 0.2e1 - t373) * t248 - (m(7) * t91 + t139) * t355 + (m(7) * t58 - t399 + t492) * t65 + (m(7) * t57 + t400 - t493) * t64 + t487 * t322 / 0.2e1;
t415 = t4 * qJD(1);
t25 = t322 * t204 + m(7) * (t250 * t91 + t322 * t57) + t250 * t139;
t411 = qJD(1) * t25;
t334 = m(5) * (t179 * t320 - t180 * t321) - t321 * t291 + t320 * t289;
t12 = t399 * t406 + t400 * t247 + m(7) * (t247 * t57 - t406 * t58) + m(6) * (t247 * t65 + t406 * t64) + (t334 - t368) * t322;
t410 = t12 * qJD(1);
t407 = t406 * t279;
t177 = t280 * t247;
t358 = 0.2e1 * t505;
t132 = t358 * m(7);
t393 = t132 * qJD(1);
t391 = mrSges(6,1) / 0.2e1 + t470;
t390 = t469 + t467;
t372 = -t204 / 0.2e1 - t207 / 0.2e1;
t371 = t205 / 0.2e1 + t206 / 0.2e1;
t370 = t462 - t208 / 0.2e1;
t369 = t210 / 0.2e1 - t211 / 0.2e1;
t187 = -t279 * mrSges(6,1) - t280 * mrSges(6,2);
t186 = -t279 * mrSges(7,1) + t280 * mrSges(7,3);
t184 = t292 * t447 + t360 * t446;
t365 = -t183 * t250 + t184 * t248;
t361 = t477 + t475 + t473;
t359 = t389 * t248;
t189 = -Ifges(7,3) * t279 - t432;
t271 = Ifges(7,5) * t279;
t190 = Ifges(7,3) * t280 - t271;
t274 = Ifges(6,4) * t280;
t191 = Ifges(6,2) * t279 - t274;
t192 = -Ifges(6,2) * t280 - t434;
t193 = -Ifges(7,1) * t280 - t271;
t194 = -Ifges(7,1) * t279 + t432;
t195 = -Ifges(6,1) * t280 + t434;
t196 = -Ifges(6,1) * t279 - t274;
t327 = t372 * t183 + (-t127 / 0.4e1 - t125 / 0.4e1 + t140 / 0.4e1 - t141 / 0.4e1) * t280 + (-t143 / 0.4e1 - t142 / 0.4e1 - t121 / 0.4e1 + t123 / 0.4e1) * t279 + (t248 * t464 + (-t58 / 0.2e1 - t64 / 0.2e1) * t280 + (-t65 / 0.2e1 + t57 / 0.2e1) * t279) * mrSges(7,2) + (mrSges(6,3) * t464 + t196 / 0.4e1 + t194 / 0.4e1 + t191 / 0.4e1 - t189 / 0.4e1) * t248 + (t192 / 0.4e1 - t195 / 0.4e1 - t193 / 0.4e1 - t190 / 0.4e1) * t250 + (-t158 * t355 + t183 * t440 - t352 * t91) * t472 + t188 * t465 + t158 * t135 / 0.2e1 + t139 * t463 + t262 * t187 / 0.2e1 + t307 * t136 / 0.2e1 + t91 * t186 / 0.2e1 + t488 * t322 / 0.4e1 + (t458 * mrSges(7,2) + t468 * t250 + t439 * t472 - t369) * t184;
t330 = (-pkin(5) * t60 + qJ(6) * t59) * t473 + pkin(5) * t462 - qJ(6) * t205 / 0.2e1 + t59 * t467 + t60 * t470 + t67 * t471 + t68 * t469;
t1 = (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t323 + t327 + t330 + t481;
t11 = t158 * t186 + t307 * t187 + (t189 / 0.2e1 - t196 / 0.2e1 - t191 / 0.2e1 - t194 / 0.2e1) * t280 + (-t190 / 0.2e1 - t195 / 0.2e1 + t192 / 0.2e1 - t193 / 0.2e1) * t279 - t500 * t352;
t351 = t1 * qJD(1) + t11 * qJD(2);
t24 = t484 * t502 + (m(7) + m(6)) * (-t183 * t279 - t184 * t280) + (m(5) * t441 + mrSges(5,3)) * t394;
t328 = (-t359 + t372) * t280 + (t369 + t491) * t279 + (-t179 * t321 - t180 * t320) * t476 + (t279 * t64 - t280 * t65 + t365) * t474 + (-t279 * t58 - t280 * t57 + t365) * t472;
t336 = t261 * t475 + t296 * t476 + t473 * t90;
t9 = (-t289 / 0.2e1 + t421 / 0.2e1) * t321 + (-t291 / 0.2e1 - t420 / 0.2e1) * t320 - t390 * t406 - t391 * t247 + t328 + t336;
t350 = qJD(1) * t9 + qJD(2) * t24;
t332 = (pkin(5) * t406 + qJ(6) * t247) * t472 + mrSges(6,2) * t461 + mrSges(7,3) * t460 + t496 * t505;
t5 = (t439 * t473 + t369 - t491) * t280 + (t440 * t473 - t359 - t372) * t279 + t332;
t349 = t5 * qJD(1);
t337 = -t183 * t406 + t184 * t247 + t262;
t329 = t297 * t477 + t337 * t475 + (t337 + t354) * t473;
t331 = (t60 * t279 + t59 * t280) * t472 + (-t67 * t279 + t68 * t280) * t474 + t486 * t476;
t10 = (t288 / 0.2e1 - t418 / 0.2e1) * t321 + (t290 / 0.2e1 + t416 / 0.2e1) * t320 + t390 * t250 + t391 * t248 + (t247 * t389 + t371) * t280 + (-t389 * t406 + t370) * t279 + t329 + t331;
t266 = t280 * mrSges(6,1);
t338 = -t279 * mrSges(6,2) + t188 + t266;
t41 = t320 * mrSges(5,1) + t321 * mrSges(5,2) + mrSges(4,3) + (m(5) + m(4)) * qJ(3) + m(6) * t307 + t445 + t338;
t348 = -qJD(1) * t10 + qJD(2) * t41;
t333 = (t158 * t250 + t184 * t322 + t279 * t91) * t472 + t188 * t458 + t139 * t456;
t339 = t60 * t473 + t417 / 0.2e1;
t21 = -mrSges(7,2) * t358 + t333 + t339;
t53 = t500 * t279;
t347 = -qJD(1) * t21 - qJD(2) * t53;
t33 = (t444 / 0.4e1 + t413 / 0.4e1 + t355 / 0.4e1) * t480 - t135 - t136;
t40 = (t443 / 0.4e1 - t412 / 0.4e1 + t352 / 0.4e1) * t480 - t186 - t187;
t346 = qJD(1) * t33 + qJD(2) * t40;
t37 = t485 * (t248 * t280 - t250 * t279);
t335 = t394 * t477 - 0.2e1 * (m(6) / 0.4e1 + m(7) / 0.4e1) * t484;
t39 = t335 + t361;
t345 = qJD(1) * t37 + qJD(2) * t39;
t13 = t399 * t250 + t400 * t248 + m(7) * (t248 * t57 - t250 * t58) + m(6) * (t248 * t65 + t250 * t64) + t334 * t323;
t344 = qJD(1) * t13 + qJD(3) * t37;
t27 = t419 + (t386 / 0.4e1 + t383 / 0.4e1 + t401 / 0.2e1 - t65 / 0.4e1) * t480;
t342 = qJD(1) * t27 + qJD(5) * t312;
t115 = m(7) * t250;
t160 = m(7) * t279;
t341 = qJD(1) * t115 + qJD(2) * t160;
t131 = m(7) * t505 - t406 * t472;
t77 = m(7) * t184 - t280 * mrSges(7,2);
t74 = m(7) * t463 + t352 * t472;
t61 = m(7) * t465 + t355 * t472;
t38 = t335 - t361;
t36 = t37 * qJD(4);
t26 = t503 * t57 + t204;
t20 = t333 - t339;
t8 = t291 * t454 + t289 * t452 + t428 / 0.2e1 + t431 / 0.2e1 - t426 / 0.2e1 + t430 / 0.2e1 + mrSges(5,2) * t405 / 0.2e1 - mrSges(5,1) * t403 / 0.2e1 + t328 - t336;
t7 = mrSges(5,1) * t402 / 0.2e1 - mrSges(5,2) * t404 / 0.2e1 - t329 + t290 * t453 + t288 * t451 - t250 * mrSges(6,2) / 0.2e1 + t425 / 0.2e1 + t248 * t471 - t429 / 0.2e1 + (m(4) * pkin(7) + mrSges(4,1)) * t323 + t370 * t279 + t371 * t280 + t331 + t502 * (-t177 / 0.2e1 + t407 / 0.2e1);
t6 = t211 * t455 + (t279 * t440 + t280 * t439) * t472 - t280 * t210 / 0.2e1 + t332 - t400 * t279 / 0.2e1 + t502 * (t248 * t456 + t250 * t455);
t2 = Ifges(6,6) * t461 + Ifges(7,6) * t460 + t494 * t497 + t501 * t505 + t327 - t330;
t14 = [qJD(2) * t3 + qJD(3) * t12 + qJD(4) * t13 + qJD(5) * t4 + qJD(6) * t25, t7 * qJD(3) + t8 * qJD(4) + t2 * qJD(5) + t20 * qJD(6) + t424 + (-qJ(3) * t263 + t158 * t137 + t307 * t138 - t183 * t208 + t183 * t209 + t184 * t205 + t184 * t206 + t90 * t188 + t190 * t460 + t192 * t461 + t261 * t266 + (-t441 * t288 - t181 * mrSges(5,3) - t296 * mrSges(5,2) + t244 / 0.2e1) * t321 + (-t441 * t290 - t182 * mrSges(5,3) - t296 * mrSges(5,1) - t243 / 0.2e1) * t320 + (-t59 * mrSges(7,2) - t68 * mrSges(6,3) - t375) * t280 + (-t60 * mrSges(7,2) + t67 * mrSges(6,3) + t340) * t279 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t451 + Ifges(5,6) * t454 + t279 * t388 - t280 * t387 - Ifges(4,4) + Ifges(3,5)) * t323 + 0.2e1 * (-qJ(3) * t296 - t486 * t441) * t476 + 0.2e1 * (-t183 * t67 + t184 * t68 + t261 * t307) * t474 + (t158 * t90 + t183 * t60 + t184 * t59) * t503 + (-qJ(3) * mrSges(4,1) + (Ifges(5,1) * t321 - t437) * t453 + (-Ifges(5,2) * t320 + t436) * t451 + Ifges(4,5) - Ifges(3,6)) * t322 + (m(4) * t356 - t323 * mrSges(3,1) + t322 * mrSges(3,2) + t294) * pkin(7) + (t194 + t196) * t505) * qJD(2), t7 * qJD(2) + t6 * qJD(5) + t131 * qJD(6) + t36 + t410 + 0.2e1 * t485 * qJD(3) * (t177 - t407) qJD(2) * t8 + qJD(5) * t61 + t344, t415 + t2 * qJD(2) + t6 * qJD(3) + t61 * qJD(4) + (t354 * mrSges(7,2) + t482 * t64 + t483 * t65 + t487) * qJD(5) + t26 * qJD(6), qJD(2) * t20 + qJD(3) * t131 + qJD(5) * t26 + t411; -qJD(3) * t10 + qJD(4) * t9 + qJD(5) * t1 + qJD(6) * t21 - t424, qJD(3) * t41 + qJD(4) * t24 + qJD(5) * t11 + qJD(6) * t53, qJD(4) * t38 + t348, qJD(3) * t38 + qJD(5) * t74 + t350, t74 * qJD(4) + (t353 * mrSges(7,2) - t482 * t183 + t483 * t184 + t488) * qJD(5) + t77 * qJD(6) + t351, qJD(5) * t77 - t347; qJD(2) * t10 - qJD(5) * t5 + qJD(6) * t132 + t36 - t410, qJD(4) * t39 - t348, 0, t345, -t338 * qJD(5) + (-t353 * qJD(5) / 0.2e1 + qJD(6) * t455) * t480 - t349, m(7) * t280 * qJD(5) + t393; -qJD(2) * t9 - qJD(5) * t33 + qJD(6) * t115 - t344, -qJD(3) * t39 - qJD(5) * t40 + qJD(6) * t160 - t350, -t345, 0, -t346, t341; -qJD(2) * t1 + qJD(3) * t5 + qJD(4) * t33 + qJD(6) * t27 - t415, qJD(4) * t40 - t351, t349, t346, t312 * qJD(6), t342; -qJD(2) * t21 - qJD(3) * t132 - qJD(4) * t115 - qJD(5) * t27 - t411, -qJD(4) * t160 + t347, -t393, -t341, -t342, 0;];
Cq  = t14;
