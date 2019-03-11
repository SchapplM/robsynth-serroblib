% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:32
% EndTime: 2019-03-09 03:44:43
% DurationCPUTime: 5.69s
% Computational Cost: add. (14065->529), mult. (26292->721), div. (0->0), fcn. (25104->8), ass. (0->278)
t327 = cos(qJ(3));
t328 = -pkin(3) - pkin(8);
t378 = -cos(pkin(10)) * pkin(1) - pkin(2);
t324 = sin(qJ(3));
t426 = qJ(4) * t324;
t359 = t378 - t426;
t229 = t327 * t328 + t359;
t305 = sin(pkin(10)) * pkin(1) + pkin(7);
t466 = pkin(4) + t305;
t257 = t466 * t324;
t323 = sin(qJ(5));
t326 = cos(qJ(5));
t142 = -t229 * t323 + t326 * t257;
t409 = t323 * t327;
t116 = pkin(9) * t409 + t142;
t106 = pkin(5) * t324 + t116;
t325 = cos(qJ(6));
t143 = t229 * t326 + t257 * t323;
t401 = t326 * t327;
t117 = -pkin(9) * t401 + t143;
t322 = sin(qJ(6));
t424 = t117 * t322;
t65 = t106 * t325 - t424;
t70 = t116 * t325 - t424;
t516 = t65 - t70;
t429 = t326 * mrSges(6,2);
t436 = t323 * mrSges(6,1);
t281 = t429 + t436;
t404 = t325 * t326;
t268 = t322 * t323 - t404;
t360 = t322 * t326 + t325 * t323;
t372 = -mrSges(7,1) * t360 + t268 * mrSges(7,2);
t343 = -t281 + t372;
t525 = mrSges(5,3) - t343;
t408 = t323 * t328;
t277 = -t323 * pkin(9) + t408;
t278 = (-pkin(9) + t328) * t326;
t184 = t277 * t325 + t278 * t322;
t371 = -t277 * t322 + t325 * t278;
t395 = -Ifges(7,5) * t360 + Ifges(7,6) * t268;
t34 = -t184 * mrSges(7,1) - t371 * mrSges(7,2) + t395;
t524 = t34 * qJD(6);
t238 = t360 * t324;
t487 = t238 / 0.2e1;
t410 = t323 * t324;
t235 = -t322 * t410 + t324 * t404;
t491 = t235 / 0.2e1;
t358 = Ifges(7,5) * t487 + Ifges(7,6) * t491;
t258 = t466 * t327;
t241 = t326 * t258;
t400 = t327 * qJ(4);
t282 = -t324 * pkin(3) + t400;
t266 = pkin(8) * t324 - t282;
t120 = pkin(5) * t327 + t241 + (-pkin(9) * t324 - t266) * t323;
t158 = t323 * t258 + t326 * t266;
t407 = t324 * t326;
t141 = pkin(9) * t407 + t158;
t76 = t120 * t325 - t141 * t322;
t77 = t120 * t322 + t141 * t325;
t510 = t77 * mrSges(7,2) / 0.2e1 - t76 * mrSges(7,1) / 0.2e1 - t358;
t507 = Ifges(7,3) * t327 / 0.2e1 - t510;
t523 = m(7) * pkin(5);
t484 = t268 / 0.2e1;
t306 = pkin(5) * t323 + qJ(4);
t476 = m(7) * t306;
t423 = t117 * t325;
t66 = t106 * t322 + t423;
t69 = -t116 * t322 - t423;
t515 = t66 + t69;
t520 = t268 * t322 + t360 * t325;
t236 = t268 * t327;
t237 = t360 * t327;
t144 = -t237 * mrSges(7,1) + t236 * mrSges(7,2);
t519 = qJD(6) * t144;
t518 = t372 * qJD(6);
t454 = Ifges(7,4) * t237;
t134 = Ifges(7,2) * t236 + t324 * Ifges(7,6) - t454;
t212 = Ifges(7,4) * t236;
t136 = -Ifges(7,1) * t237 + t324 * Ifges(7,5) + t212;
t147 = Ifges(7,2) * t237 + t212;
t148 = Ifges(7,1) * t236 + t454;
t176 = -mrSges(7,1) * t268 - mrSges(7,2) * t360;
t199 = pkin(5) * t401 + t258;
t439 = t237 * mrSges(7,3);
t192 = mrSges(7,1) * t324 + t439;
t492 = -t192 / 0.2e1;
t441 = t236 * mrSges(7,3);
t190 = -mrSges(7,2) * t324 + t441;
t493 = t190 / 0.2e1;
t517 = t184 * t492 + t371 * t493 + (t134 / 0.4e1 - t148 / 0.4e1) * t268 - (t147 / 0.4e1 + t136 / 0.4e1) * t360 + t199 * t176 / 0.2e1 + t306 * t144 / 0.2e1 + t324 * t395 / 0.4e1;
t450 = Ifges(6,6) * t326;
t452 = Ifges(6,5) * t323;
t513 = Ifges(4,4) - t452 / 0.2e1 - t450 / 0.2e1 + Ifges(5,6);
t319 = t323 ^ 2;
t320 = t326 ^ 2;
t394 = t319 + t320;
t511 = -t76 * t268 + t360 * t77;
t509 = t360 * t192 / 0.2e1 + t190 * t484;
t438 = t238 * mrSges(7,2);
t443 = t235 * mrSges(7,1);
t397 = t443 / 0.2e1 - t438 / 0.2e1;
t506 = 0.2e1 * m(7);
t505 = t236 ^ 2;
t504 = t237 ^ 2;
t503 = 2 * qJD(3);
t502 = m(6) / 0.2e1;
t501 = -m(7) / 0.2e1;
t500 = m(7) / 0.2e1;
t499 = pkin(5) / 0.2e1;
t85 = t235 * t237 + t236 * t238;
t497 = m(7) * t85;
t496 = -qJ(4) / 0.2e1;
t440 = t237 * mrSges(7,2);
t442 = t236 * mrSges(7,1);
t146 = -t440 - t442;
t495 = t146 / 0.2e1;
t494 = -t371 / 0.2e1;
t490 = t236 / 0.2e1;
t489 = -t237 / 0.2e1;
t488 = t237 / 0.2e1;
t486 = t258 / 0.2e1;
t485 = -t268 / 0.2e1;
t483 = -t360 / 0.2e1;
t482 = -t323 / 0.2e1;
t481 = t323 / 0.2e1;
t480 = t324 / 0.2e1;
t479 = -t326 / 0.2e1;
t478 = t326 / 0.2e1;
t477 = m(7) * t199;
t475 = pkin(3) * t327;
t474 = pkin(5) * t326;
t473 = t65 * mrSges(7,2);
t472 = t66 * mrSges(7,1);
t471 = t69 * mrSges(7,1);
t470 = t70 * mrSges(7,2);
t459 = mrSges(6,1) * t327;
t458 = mrSges(6,2) * t327;
t457 = mrSges(6,3) * t327;
t456 = Ifges(6,4) * t323;
t455 = Ifges(6,4) * t326;
t453 = Ifges(7,4) * t268;
t451 = Ifges(6,6) * t324;
t448 = pkin(5) * qJD(5);
t437 = t322 * t77;
t435 = t323 * mrSges(6,2);
t434 = t324 * mrSges(6,1);
t433 = t324 * mrSges(6,2);
t432 = t324 * Ifges(6,5);
t431 = t325 * t76;
t430 = t326 * mrSges(6,1);
t341 = t144 * t480 + t190 * t490 + t192 * t488;
t273 = -mrSges(6,3) * t401 - t433;
t403 = t326 * t273;
t271 = mrSges(6,3) * t409 + t434;
t411 = t323 * t271;
t353 = t403 / 0.2e1 - t411 / 0.2e1;
t374 = t319 / 0.2e1 + t320 / 0.2e1;
t370 = mrSges(6,3) * t374;
t377 = t410 / 0.2e1;
t10 = (t236 * t515 + t237 * t516) * t501 + (t504 / 0.2e1 + t505 / 0.2e1) * mrSges(7,3) + (t281 * t480 + t327 * t370 + t377 * t523 + t353) * t327 - t341;
t425 = t10 * qJD(1);
t15 = t341 - (t504 + t505) * mrSges(7,3) / 0.2e1;
t422 = t15 * qJD(1);
t356 = mrSges(7,1) * t487 + mrSges(7,2) * t491;
t18 = (t236 * t485 + t237 * t483) * mrSges(7,3) + t356 + t509;
t421 = t18 * qJD(1);
t420 = t235 * t360;
t419 = t238 * t268;
t253 = t359 - t475;
t280 = t327 * mrSges(5,2) - t324 * mrSges(5,3);
t373 = -m(5) * t253 - t280;
t24 = -t235 * t190 + t238 * t192 + m(7) * (-t235 * t66 + t238 * t65) + (-t403 + t411 + m(6) * (t142 * t323 - t143 * t326) + t373) * t324;
t418 = t24 * qJD(1);
t189 = -mrSges(7,2) * t327 + mrSges(7,3) * t235;
t413 = t322 * t189;
t412 = t322 * t237;
t308 = t324 * t327;
t191 = mrSges(7,1) * t327 - t238 * mrSges(7,3);
t406 = t325 * t191;
t405 = t325 * t236;
t284 = -Ifges(6,2) * t323 + t455;
t402 = t326 * t284;
t396 = Ifges(7,5) * t236 + Ifges(7,6) * t237;
t393 = qJD(3) * t324;
t392 = qJD(3) * t327;
t391 = m(7) * t499;
t390 = -t497 / 0.2e1;
t389 = t497 / 0.2e1;
t386 = -t65 / 0.2e1 + t70 / 0.2e1;
t385 = t66 / 0.2e1 + t69 / 0.2e1;
t381 = -t439 / 0.2e1;
t380 = mrSges(7,3) * t484;
t379 = t176 * t480 + t237 * t380 + t268 * t381;
t179 = -Ifges(7,2) * t360 - t453;
t180 = -Ifges(7,1) * t360 + t453;
t376 = t179 / 0.4e1 - t180 / 0.4e1;
t262 = Ifges(7,4) * t360;
t178 = Ifges(7,2) * t268 - t262;
t181 = -Ifges(7,1) * t268 - t262;
t375 = t181 / 0.4e1 + t178 / 0.4e1;
t279 = t430 - t435;
t286 = Ifges(6,1) * t326 - t456;
t369 = Ifges(6,1) * t323 + t455;
t368 = Ifges(6,2) * t326 + t456;
t133 = Ifges(7,4) * t238 + Ifges(7,2) * t235 + Ifges(7,6) * t327;
t135 = Ifges(7,1) * t238 + Ifges(7,4) * t235 + Ifges(7,5) * t327;
t145 = t438 - t443;
t157 = -t266 * t323 + t241;
t198 = (-t466 - t474) * t324;
t230 = Ifges(6,6) * t327 + t324 * t368;
t231 = -t327 * t368 + t451;
t232 = Ifges(6,5) * t327 + t324 * t369;
t233 = -t327 * t369 + t432;
t254 = t324 * t279;
t270 = -mrSges(6,3) * t410 + t459;
t272 = mrSges(6,3) * t407 - t458;
t4 = t142 * t270 + t157 * t271 + t143 * t272 + t158 * t273 - t258 * t254 + t135 * t489 + t136 * t487 + t134 * t491 + t133 * t490 + t198 * t146 + t199 * t145 + t66 * t189 + t77 * t190 + t65 * t191 + t76 * t192 + t373 * t282 + m(6) * (t142 * t157 + t143 * t158 - t257 * t258) + m(7) * (t198 * t199 + t65 * t76 + t66 * t77) + (t378 * mrSges(4,1) - t253 * mrSges(5,2) + t231 * t478 + t233 * t481 - t513 * t324 + t358) * t324 + (t230 * t479 + t232 * t482 - t257 * t279 + t378 * mrSges(4,2) - t253 * mrSges(5,3) + Ifges(7,5) * t489 + Ifges(7,6) * t490 + (Ifges(6,3) + Ifges(5,2) + Ifges(4,1) - Ifges(5,3) - Ifges(4,2) + Ifges(7,3)) * t324 + t513 * t327) * t327;
t354 = t270 * t479 + t272 * t482;
t355 = -t435 / 0.2e1 + t430 / 0.2e1;
t362 = t157 * t326 + t158 * t323;
t9 = t190 * t487 + t189 * t489 + t192 * t491 + t191 * t490 + (t355 * t327 + t354 + t495) * t327 + (t273 * t481 + t271 * t478 + t145 / 0.2e1 - t254 / 0.2e1) * t324 + (t198 * t324 + t199 * t327 + t235 * t65 + t236 * t76 - t237 * t77 + t238 * t66) * t500 + ((t258 - t362) * t327 + (t142 * t326 + t143 * t323 - t257) * t324) * t502;
t367 = t4 * qJD(1) + t9 * qJD(2);
t255 = t327 * t284;
t256 = t327 * t286;
t301 = Ifges(6,6) * t409;
t334 = t199 * t144 + (t136 / 0.2e1 + t147 / 0.2e1) * t236 + (t66 * mrSges(7,3) - t148 / 0.2e1 + t134 / 0.2e1) * t237 - t65 * t441 + t396 * t480;
t5 = t70 * t190 + t69 * t192 + m(7) * (t65 * t69 + t66 * t70) + t142 * t273 - t143 * t271 + t301 * t480 + ((t142 * mrSges(6,3) - t432 / 0.2e1 - t258 * mrSges(6,2) - t233 / 0.2e1 + t255 / 0.2e1) * t326 + (t143 * mrSges(6,3) - t258 * mrSges(6,1) + t256 / 0.2e1 + t231 / 0.2e1 + (-t146 - t477) * pkin(5)) * t323) * t327 + t334;
t366 = t5 * qJD(1) - t10 * qJD(2);
t8 = t65 * t190 - t66 * t192 + t334;
t365 = t8 * qJD(1) + t15 * qJD(2);
t62 = m(7) * (t235 * t236 - t237 * t238 + t308) + m(6) * (-t394 + 0.1e1) * t308;
t364 = t9 * qJD(1) + t62 * qJD(2);
t342 = m(7) * (-t268 * t515 - t360 * t516);
t344 = (-t235 * t322 + t238 * t325) * t391;
t11 = (t433 / 0.2e1 - t273 / 0.2e1) * t326 + (t434 / 0.2e1 + t271 / 0.2e1) * t323 - t374 * t457 - t342 / 0.2e1 + t344 + t18;
t363 = t11 * qJD(1);
t361 = t235 * t268 - t238 * t360;
t357 = -t157 * mrSges(6,1) / 0.2e1 + t158 * mrSges(6,2) / 0.2e1;
t352 = -t231 / 0.4e1 - t256 / 0.4e1 + t328 * t273 / 0.2e1;
t351 = -t233 / 0.4e1 + t255 / 0.4e1 - t328 * t271 / 0.2e1;
t348 = t361 * t501;
t31 = t306 * t176 - (t178 / 0.2e1 + t181 / 0.2e1) * t360 + (t179 / 0.2e1 - t180 / 0.2e1) * t268;
t32 = t379 - t397;
t330 = (t184 * mrSges(7,3) / 0.2e1 + t376) * t237 + (mrSges(7,3) * t494 + t375) * t236 + t517;
t6 = t330 - t507;
t347 = t6 * qJD(1) + t32 * qJD(2) + t31 * qJD(3);
t101 = t476 + (m(6) + m(5)) * qJ(4) + t525;
t331 = (t420 / 0.2e1 + t419 / 0.2e1) * mrSges(7,3) + m(6) * t486 + (-t184 * t235 + t238 * t371 + t199) * t500 - t442 / 0.2e1 - t440 / 0.2e1;
t333 = -m(6) * t362 / 0.2e1 + t511 * t501 + t191 * t484 + t189 * t483;
t17 = (t459 / 0.2e1 - t270 / 0.2e1) * t326 + (-t458 / 0.2e1 - t272 / 0.2e1) * t323 + t331 + t333;
t78 = t348 + (t501 + (-0.1e1 / 0.2e1 + t374) * m(6)) * t324;
t346 = -qJD(1) * t17 + qJD(2) * t78 - qJD(3) * t101;
t335 = (t322 * t492 + t325 * t493 + (-t405 / 0.2e1 + t412 / 0.2e1) * mrSges(7,3)) * pkin(5);
t14 = -mrSges(7,1) * t385 + mrSges(7,2) * t386 + t335;
t274 = (mrSges(7,1) * t322 + mrSges(7,2) * t325) * pkin(5);
t35 = (t494 + t371 / 0.2e1) * mrSges(7,2);
t345 = -t14 * qJD(1) - t35 * qJD(3) + t274 * qJD(5);
t19 = t236 * t380 - t360 * t381 + t356 - t509;
t329 = t376 * t237 + t375 * t236 + (t184 * t488 + t236 * t494 + t268 * t385 - t360 * t386) * mrSges(7,3) + (-t516 * t184 + t515 * t371) * t500 + t279 * t486 + t517;
t336 = (-t286 / 0.4e1 + t368 / 0.4e1 + mrSges(6,2) * t496) * t326 + t328 * t370;
t338 = t369 / 0.4e1 + t284 / 0.4e1 + mrSges(6,1) * t496 + (-t476 / 0.2e1 + t372 / 0.2e1) * pkin(5);
t1 = t329 + (t146 * t478 - t406 / 0.2e1 - t413 / 0.2e1 + (-t437 / 0.4e1 - t431 / 0.4e1 + t199 * t326 / 0.4e1) * t506) * pkin(5) + (-0.3e1 / 0.4e1 * t432 + t338 * t327 + t351) * t323 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t336) * t327 + (-0.3e1 / 0.4e1 * t451 + t352) * t326 + t357 + t510;
t20 = qJ(4) * t279 - t369 * t478 - t402 / 0.2e1 + (-t286 / 0.2e1 + t368 / 0.2e1) * t323 + t31 + (-t372 + t476) * t474;
t337 = t407 * t523;
t339 = (t235 * t325 + t238 * t322) * t391 + t397;
t29 = (-t176 / 0.2e1 - t279 / 0.2e1 + t355) * t324 - t337 / 0.2e1 + t339;
t340 = t1 * qJD(1) - t29 * qJD(2) + t20 * qJD(3);
t267 = t274 * qJD(6);
t74 = t348 + (t394 * t502 + m(5)) * t324 + (m(7) + m(6)) * t480;
t33 = t379 + t397;
t30 = t337 / 0.2e1 + t279 * t480 + t355 * t324 + t339 + t379;
t16 = (m(5) * t305 + mrSges(5,1) + t355) * t327 + t331 - t333 - t354;
t13 = -t473 / 0.2e1 - t472 / 0.2e1 - t470 / 0.2e1 + t471 / 0.2e1 + t335 + t396;
t12 = t342 / 0.2e1 + (t429 / 0.2e1 + t436 / 0.2e1) * t324 + t344 + t19 + t353 + t394 * t457 / 0.2e1;
t7 = t330 + t507;
t3 = t9 * qJD(3) + qJD(4) * t389 - t10 * qJD(5) + t15 * qJD(6);
t2 = (t431 + t437) * t391 + (-t432 / 0.4e1 + t351) * t323 + Ifges(6,6) * t407 / 0.2e1 + Ifges(6,5) * t377 + t329 + (-t451 / 0.4e1 + (t495 + t477 / 0.2e1) * pkin(5) + t352) * t326 - t357 + (t406 + t413) * t499 + (t323 * t338 + t336 + Ifges(6,3) / 0.2e1) * t327 + t507;
t21 = [qJD(3) * t4 + qJD(4) * t24 + qJD(5) * t5 + qJD(6) * t8, t3, t16 * qJD(4) + t2 * qJD(5) + t7 * qJD(6) + ((-qJ(4) * t257 + t328 * t362) * t502 + (t184 * t77 + t198 * t306 + t371 * t76) * t500) * t503 + (-pkin(3) * mrSges(5,1) + Ifges(6,5) * t478 + Ifges(7,5) * t485 + Ifges(6,6) * t482 + Ifges(7,6) * t483 - Ifges(5,4) + Ifges(4,5)) * t392 + (t402 / 0.2e1 - qJ(4) * mrSges(5,1) + t286 * t481 + Ifges(5,5) - Ifges(4,6)) * t393 + t367 + (-qJ(4) * t254 + t133 * t483 + t135 * t485 + t306 * t145 - t198 * t372 + t179 * t491 + t181 * t487 + t371 * t191 + t184 * t189 - t257 * t281 + (t328 * t270 - t157 * mrSges(6,3) + t232 / 0.2e1) * t326 + (t328 * t272 - t158 * mrSges(6,3) - t230 / 0.2e1) * t323 + (m(5) * (-t426 - t475) - t327 * mrSges(4,1) + t324 * mrSges(4,2) + t280) * t305 - t511 * mrSges(7,3)) * qJD(3), t418 + t16 * qJD(3) + t12 * qJD(5) + t19 * qJD(6) + (t85 * qJD(2) / 0.4e1 + (-t419 - t420) * qJD(4) / 0.2e1) * t506, t2 * qJD(3) + t12 * qJD(4) + (-t143 * mrSges(6,1) - t142 * mrSges(6,2) - Ifges(6,5) * t401 + t301 + t396 - t470 + t471) * qJD(5) + t13 * qJD(6) + (m(7) * (t322 * t70 + t325 * t69) + (-t405 + t412) * mrSges(7,3)) * t448 + t366, t7 * qJD(3) + t19 * qJD(4) + t13 * qJD(5) + (t396 - t472 - t473) * qJD(6) + t365; t3, t62 * qJD(3), t74 * qJD(4) + t30 * qJD(5) + t33 * qJD(6) + t361 * qJD(3) * mrSges(7,3) + (-mrSges(4,2) + t525) * t392 + (-mrSges(6,3) * t394 - mrSges(4,1) + mrSges(5,2)) * t393 + ((t184 * t238 + t235 * t371 + t306 * t327) * t500 + (t324 * t328 * t394 + t400) * t502 + m(5) * t282 / 0.2e1) * t503 + t364, qJD(1) * t389 + t74 * qJD(3), -t425 + t30 * qJD(3) + (t327 * t281 - t144 + (t236 * t322 + t237 * t325) * t523) * qJD(5) - t519, qJD(3) * t33 - qJD(5) * t144 + t422 - t519; qJD(4) * t17 + qJD(5) * t1 + qJD(6) * t6 - t367, -qJD(4) * t78 - qJD(5) * t29 + qJD(6) * t32 - t364, qJD(4) * t101 + qJD(5) * t20 + qJD(6) * t31, -t346 (-mrSges(6,1) * t408 - t328 * t429 + t34 - t450 - t452) * qJD(5) + t524 + (m(7) * (-t184 * t325 + t322 * t371) + t520 * mrSges(7,3)) * t448 + t340, t34 * qJD(5) + t347 + t524; qJD(2) * t390 - t17 * qJD(3) - t11 * qJD(5) - t18 * qJD(6) - t418, qJD(1) * t390 + t78 * qJD(3), t346, 0 (-t520 * t523 + t343) * qJD(5) + t518 - t363, qJD(5) * t372 - t421 + t518; -qJD(3) * t1 + qJD(4) * t11 + qJD(6) * t14 - t366, t29 * qJD(3) + t425, qJD(6) * t35 - t340, t363, -t267, -t267 - t345; -qJD(3) * t6 + qJD(4) * t18 - qJD(5) * t14 - t365, -t32 * qJD(3) - t422, -qJD(5) * t35 - t347, t421, t345, 0;];
Cq  = t21;
