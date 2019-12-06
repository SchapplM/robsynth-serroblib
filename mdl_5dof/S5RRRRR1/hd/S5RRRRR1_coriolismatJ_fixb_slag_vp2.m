% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:49:51
% EndTime: 2019-12-05 18:50:06
% DurationCPUTime: 6.11s
% Computational Cost: add. (9863->386), mult. (21102->532), div. (0->0), fcn. (24558->8), ass. (0->240)
t246 = sin(qJ(4));
t425 = pkin(3) * t246;
t239 = pkin(6) + t425;
t247 = cos(qJ(5));
t248 = cos(qJ(4));
t426 = sin(qJ(3));
t427 = sin(qJ(2));
t428 = cos(qJ(3));
t429 = cos(qJ(2));
t281 = t426 * t427 - t428 * t429;
t282 = t426 * t429 + t428 * t427;
t455 = t246 * t281 - t248 * t282;
t188 = t246 * t282 + t248 * t281;
t245 = sin(qJ(5));
t492 = t245 * t188;
t502 = -mrSges(6,2) * t455 - mrSges(6,3) * t492;
t523 = t247 * t502;
t491 = t247 * t188;
t501 = mrSges(6,1) * t455 - mrSges(6,3) * t491;
t524 = t245 * t501;
t471 = t523 / 0.2e1 - t524 / 0.2e1;
t527 = t471 * t239;
t317 = t523 - t524;
t331 = t429 * pkin(2) + pkin(1);
t269 = -t281 * pkin(3) + t331;
t263 = -pkin(4) * t188 + t269;
t261 = -pkin(6) * t455 + t263;
t255 = t247 * t261;
t386 = t455 * t247;
t357 = Ifges(6,4) * t386;
t475 = t455 * Ifges(6,6);
t421 = Ifges(6,4) * t247;
t322 = -Ifges(6,2) * t245 + t421;
t486 = t322 * t188 + t475;
t496 = t188 ^ 2;
t514 = -t492 / 0.2e1;
t526 = -(t269 * mrSges(5,2) - t245 * t357 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2)) * t455) * t188 + (t455 ^ 2 - t496) * Ifges(5,4) - t269 * mrSges(5,1) * t455 - t501 * t255 + (-Ifges(6,6) * t496 + (Ifges(6,2) * t514 + t475 / 0.2e1 + t486 / 0.2e1) * t455 - t502 * t261) * t245;
t525 = t501 / 0.2e1;
t232 = t245 * Ifges(6,1) + t421;
t360 = t247 * t232;
t422 = Ifges(6,4) * t245;
t231 = t247 * Ifges(6,2) + t422;
t366 = t245 * t231;
t305 = t366 / 0.2e1 - t360 / 0.2e1;
t323 = Ifges(6,1) * t247 - t422;
t431 = -t245 / 0.2e1;
t257 = t305 - t247 * t322 / 0.2e1 + t323 * t431;
t242 = Ifges(6,5) * t247;
t414 = Ifges(6,6) * t245;
t522 = Ifges(6,6) * t514 + t491 * Ifges(6,5) / 0.2e1 - (t242 / 0.2e1 - t414 / 0.2e1) * t188;
t230 = Ifges(6,5) * t245 + Ifges(6,6) * t247;
t473 = t455 * t230;
t477 = Ifges(5,6) * t455;
t507 = t247 * t486;
t476 = t455 * Ifges(6,5);
t485 = t323 * t188 + t476;
t508 = t245 * t485;
t520 = t473 / 0.2e1 - t477 + t507 / 0.2e1 + t508 / 0.2e1;
t441 = t188 / 0.4e1;
t442 = -t188 / 0.4e1;
t519 = t360 * t441 + t366 * t442 + t507 / 0.4e1 + t508 / 0.4e1;
t241 = t428 * pkin(2) + pkin(3);
t351 = t426 * t248;
t211 = pkin(2) * t351 + t246 * t241;
t206 = pkin(6) + t211;
t430 = t247 / 0.2e1;
t352 = t426 * t246;
t210 = -pkin(2) * t352 + t248 * t241;
t205 = -pkin(4) - t210;
t399 = t247 * mrSges(6,2);
t403 = t245 * mrSges(6,1);
t324 = t399 + t403;
t490 = t324 * t188;
t509 = t205 * t490;
t518 = t509 / 0.2e1 + (t502 * t430 + t501 * t431) * t206;
t516 = -t485 / 0.4e1;
t515 = -t486 / 0.4e1;
t493 = Ifges(5,5) * t188;
t513 = -t493 / 0.2e1;
t512 = Ifges(6,5) * t496;
t506 = -t473 / 0.4e1 + t477 / 0.2e1 + t513;
t505 = -t476 / 0.2e1 - t485 / 0.2e1;
t243 = t245 ^ 2;
t244 = t247 ^ 2;
t359 = t243 + t244;
t489 = t359 * t248;
t500 = t205 * t246 + t206 * t489;
t416 = Ifges(6,6) * t188;
t327 = pkin(4) * t455 - t188 * pkin(6);
t483 = pkin(4) / 0.2e1;
t456 = t359 * mrSges(6,3);
t474 = -mrSges(5,2) + t456;
t424 = pkin(3) * t248;
t240 = -pkin(4) - t424;
t336 = t359 * t210;
t470 = t240 * t211 + t239 * t336;
t469 = -mrSges(6,1) / 0.2e1;
t468 = mrSges(6,1) / 0.2e1;
t466 = mrSges(6,3) * (t244 / 0.2e1 + t243 / 0.2e1);
t325 = mrSges(6,1) * t247 - mrSges(6,2) * t245;
t465 = -mrSges(5,1) - t325;
t216 = (t428 * t248 - t352) * pkin(2);
t408 = t188 * mrSges(6,2);
t130 = -mrSges(6,3) * t245 * t455 + t408;
t409 = t188 * mrSges(6,1);
t133 = -mrSges(6,3) * t386 - t409;
t369 = t245 * t133;
t304 = t130 * t430 - t369 / 0.2e1;
t463 = t304 * t216;
t461 = m(5) * t269 - mrSges(5,1) * t188 + mrSges(5,2) * t455;
t380 = t211 * t455;
t382 = t210 * t188;
t459 = -t380 - t382;
t454 = m(6) * t261 * t359;
t453 = t474 * t216 + (-t426 * mrSges(4,1) - t428 * mrSges(4,2)) * pkin(2);
t452 = t490 * t483 + t506;
t450 = m(6) / 0.2e1;
t449 = -pkin(4) / 0.2e1;
t448 = -mrSges(6,2) / 0.2e1;
t446 = Ifges(6,1) / 0.2e1;
t444 = -t130 / 0.2e1;
t443 = -t133 / 0.2e1;
t440 = -t206 / 0.2e1;
t439 = t210 / 0.2e1;
t215 = (t428 * t246 + t351) * pkin(2);
t438 = -t215 / 0.2e1;
t437 = -t216 / 0.2e1;
t436 = -t231 / 0.4e1;
t435 = -t232 / 0.4e1;
t434 = -t239 / 0.2e1;
t433 = t239 / 0.2e1;
t432 = t240 / 0.2e1;
t419 = Ifges(6,5) * t188;
t413 = Ifges(6,3) * t455;
t249 = t282 ^ 2 * Ifges(4,4) + (t512 + (-Ifges(6,1) * t491 / 0.2e1 + t505) * t455) * t247 - t331 * (-t282 * mrSges(4,1) + t281 * mrSges(4,2)) + (-Ifges(4,4) * t281 + (Ifges(4,1) - Ifges(4,2)) * t282) * t281 + t526;
t278 = t282 * pkin(3);
t268 = t278 - t327;
t358 = t427 * pkin(2);
t262 = t358 + t268;
t258 = t247 * t262;
t259 = t245 * t262;
t1 = -pkin(1) * (-t427 * mrSges(3,1) - t429 * mrSges(3,2)) + (t427 ^ 2 - t429 ^ 2) * Ifges(3,4) + t133 * t258 + t130 * t259 + t249 + (m(4) * t331 - t281 * mrSges(4,1) - t282 * mrSges(4,2)) * t358 + (Ifges(3,2) - Ifges(3,1)) * t429 * t427 + t461 * (t358 + t278) + t262 * t454;
t412 = t1 * qJD(1);
t264 = t247 * t268;
t265 = t245 * t268;
t2 = t130 * t265 + t133 * t264 + t268 * t454 + t461 * t278 + t249;
t407 = t2 * qJD(1);
t406 = t210 * mrSges(5,2);
t405 = t211 * mrSges(5,1);
t404 = t215 * mrSges(5,1);
t400 = t246 * mrSges(5,1);
t396 = t248 * mrSges(5,2);
t308 = t247 * t327;
t309 = t245 * t327;
t3 = -t130 * t309 - t133 * t308 + (t512 + (-t491 * t446 + t505) * t455) * t247 - t327 * t454 + t526;
t395 = t3 * qJD(1);
t16 = -t130 * t255 + t263 * t369 + ((t357 - t416) * t247 + (-t419 - pkin(6) * t133 + (-t422 + (Ifges(6,1) - Ifges(6,2)) * t247) * t455) * t245) * t455;
t394 = t16 * qJD(1);
t391 = t455 * t246;
t388 = t188 * t248;
t379 = t211 * t325;
t378 = t215 * t455;
t377 = t215 * t325;
t376 = t216 * t188;
t375 = t240 * t490;
t365 = t246 * t325;
t125 = t325 * t455;
t356 = t125 * t449;
t350 = t242 * t442;
t349 = t205 * t125 / 0.2e1;
t348 = t125 * t432;
t337 = t242 - t414;
t335 = t359 * t216;
t310 = pkin(4) * t324;
t321 = -t310 / 0.2e1 - t257;
t274 = mrSges(6,3) * t336 - t379 - t405 - t406;
t22 = m(6) * (t205 * t211 + t206 * t336) + t274;
t280 = t471 * pkin(6) - t452 + t519;
t4 = (t133 * t439 - t188 * t436 + t206 * t525 + t516) * t245 + (t210 * t444 + t232 * t442 + t440 * t502 + t515) * t247 + (-t211 * t324 / 0.2e1 + Ifges(5,6) / 0.2e1 - t230 / 0.4e1) * t455 - t509 / 0.2e1 + t513 + t280;
t320 = -t4 * qJD(1) + t22 * qJD(2);
t303 = t205 * t324;
t110 = t303 - t257;
t283 = -0.3e1 / 0.2e1 * t421 + t435 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.4e1) * t245;
t267 = t283 * t455 + t416;
t297 = t436 + (t446 - Ifges(6,2) / 0.4e1) * t247;
t271 = (t442 - t188 / 0.2e1) * Ifges(6,5) + t297 * t455;
t307 = -t413 / 0.2e1 + t350;
t315 = t206 * t466;
t9 = t349 - t455 * t315 + (t133 * t440 + t262 * t468 + t271) * t247 + (t130 * t440 + t262 * t448 + t267) * t245 + t307;
t319 = t9 * qJD(1) + t110 * qJD(2);
t294 = (-t391 / 0.2e1 - t388 / 0.2e1) * pkin(3);
t298 = t324 * t455;
t15 = -t463 + t298 * t438 + (t432 - t205 / 0.2e1) * t490 + (t380 / 0.2e1 + t382 / 0.2e1 - t376 / 0.2e1 - t378 / 0.2e1 + t294) * mrSges(5,3) + t317 * (t433 + t440);
t23 = t465 * t215 + m(6) * (t205 * t215 + t206 * t335) + m(5) * (-t210 * t215 + t211 * t216) + t453;
t318 = -t15 * qJD(1) + t23 * qJD(2);
t316 = pkin(6) * t466;
t314 = t239 * t466;
t313 = -t399 / 0.2e1 - t403 / 0.2e1;
t306 = t130 * t431 + t247 * t443;
t302 = t240 * t324;
t296 = -t506 + t519;
t295 = t310 / 0.2e1;
t293 = t313 * t210;
t292 = t313 * t216;
t291 = -t303 / 0.2e1;
t290 = t298 / 0.2e1;
t273 = mrSges(6,3) * t489 - t365 - t396 - t400;
t111 = (m(6) * (t239 * t489 + t240 * t246) + t273) * pkin(3);
t286 = m(6) * (-pkin(4) * t215 + pkin(6) * t335);
t18 = (t500 * pkin(3) + t470) * t450 - t286 / 0.2e1 + t465 * (t425 / 0.2e1 + t211 / 0.2e1 + t438) + t474 * (t424 / 0.2e1 + t439 + t437);
t253 = (t246 * t290 + t304 * t248) * pkin(3) + t490 * t432 + t296;
t7 = t253 + (t188 * t435 + t515 + (-pkin(6) / 0.2e1 + t433) * t502) * t247 + (pkin(6) * t525 + t231 * t441 + t501 * t434 + t516) * t245 + t452;
t289 = t7 * qJD(1) + t18 * qJD(2) + t111 * qJD(3);
t11 = t348 - t455 * t314 + (t133 * t434 + t268 * t468 + t271) * t247 + (t130 * t434 + t268 * t448 + t267) * t245 + t307;
t134 = t302 - t257;
t252 = -t302 / 0.2e1 + t257;
t50 = t291 + t292 + t252;
t288 = t11 * qJD(1) - t50 * qJD(2) + t134 * qJD(3);
t287 = t313 * t424;
t285 = Ifges(4,5) * t281 + Ifges(4,6) * t282 - t305 * t188 + t493 + t520;
t13 = t356 + t350 + (-Ifges(6,3) / 0.2e1 - t316) * t455 + (-0.3e1 / 0.4e1 * t419 + (t443 + t409 / 0.2e1) * pkin(6) + (mrSges(6,1) * t449 + t297) * t455) * t247 + (t416 + (t444 - t408 / 0.2e1) * pkin(6) + (mrSges(6,2) * t483 + t283) * t455) * t245;
t135 = t310 + t257;
t52 = t295 + t291 + t293 + t257;
t77 = t295 + t287 + t252;
t279 = t13 * qJD(1) - t52 * qJD(2) - t77 * qJD(3) - t135 * qJD(4);
t272 = t413 / 0.2e1 + t522;
t260 = t283 * t245 + t297 * t247;
t209 = t302 / 0.2e1;
t190 = t303 / 0.2e1;
t78 = t209 + t287 + t321;
t53 = t190 + t293 + t321;
t51 = t209 + t190 + t292 - t257;
t17 = -t379 / 0.2e1 + t470 * t450 - t406 / 0.2e1 - t405 / 0.2e1 - t377 / 0.2e1 + t286 / 0.2e1 + mrSges(5,2) * t437 - t404 / 0.2e1 + t210 * t466 + (-t365 / 0.2e1 + t500 * t450 - t396 / 0.2e1 - t400 / 0.2e1 + t248 * t466) * pkin(3) + t216 * t456 / 0.2e1;
t14 = t356 + t309 * t448 + t308 * t468 + t306 * pkin(6) + (Ifges(6,3) / 0.2e1 - t316 + t260) * t455 + t522;
t12 = t348 - t265 * t448 + t264 * t469 + t306 * t239 + (-t314 + t260) * t455 + t272;
t10 = t349 - t259 * t448 + t258 * t469 + t306 * t206 + (-t315 + t260) * t455 + t272;
t8 = t375 / 0.2e1 + t285 + t527 + t215 * t290 + t463 + (t294 + t378 / 0.2e1) * mrSges(5,3) - (-t376 - t459) * mrSges(5,3) / 0.2e1 + t518;
t6 = t253 + t280 + t527;
t5 = t304 * t210 + t211 * t290 + t280 + t296 + t518;
t19 = [-qJD(2) * t1 - qJD(3) * t2 - qJD(4) * t3 - qJD(5) * t16, -t412 + (t459 * mrSges(5,3) + (-t428 * t281 + t426 * t282) * mrSges(4,3) * pkin(2) - Ifges(3,5) * t429 + Ifges(3,6) * t427 + t509 + t285 + t317 * t206) * qJD(2) + t8 * qJD(3) + t5 * qJD(4) + t10 * qJD(5), -t407 + t8 * qJD(2) + (t375 + t317 * t239 + (-t388 - t391) * pkin(3) * mrSges(5,3) + t285) * qJD(3) + t6 * qJD(4) + t12 * qJD(5), t5 * qJD(2) + t6 * qJD(3) + t14 * qJD(5) - t395 + (-pkin(4) * t490 - (-Ifges(5,5) + t305) * t188 + t317 * pkin(6) + t520) * qJD(4), -t394 + t10 * qJD(2) + t12 * qJD(3) + t14 * qJD(4) + ((-mrSges(6,2) * t261 - t475) * t247 + (-mrSges(6,1) * t261 - t476) * t245) * qJD(5); -qJD(3) * t15 - qJD(4) * t4 + qJD(5) * t9 + t412, qJD(3) * t23 + qJD(4) * t22 + qJD(5) * t110, t17 * qJD(4) + t51 * qJD(5) + t318 + (-t377 - t404 + 0.2e1 * (t215 * t240 + t239 * t335) * t450 + m(5) * (-t215 * t248 + t216 * t246) * pkin(3) + t453) * qJD(3), t17 * qJD(3) + (m(6) * (-pkin(4) * t211 + pkin(6) * t336) + t274) * qJD(4) + t53 * qJD(5) + t320, t51 * qJD(3) + t53 * qJD(4) + (-t206 * t325 + t337) * qJD(5) + t319; qJD(2) * t15 + qJD(4) * t7 + qJD(5) * t11 + t407, qJD(4) * t18 - qJD(5) * t50 - t318, qJD(4) * t111 + qJD(5) * t134, t78 * qJD(5) + (m(6) * (-pkin(4) * t246 + pkin(6) * t489) + t273) * qJD(4) * pkin(3) + t289, t78 * qJD(4) + (-t239 * t325 + t337) * qJD(5) + t288; qJD(2) * t4 - qJD(3) * t7 + qJD(5) * t13 + t395, -qJD(3) * t18 - qJD(5) * t52 - t320, -qJD(5) * t77 - t289, -t135 * qJD(5), (-pkin(6) * t325 + t337) * qJD(5) + t279; -qJD(2) * t9 - qJD(3) * t11 - qJD(4) * t13 + t394, qJD(3) * t50 + qJD(4) * t52 - t319, qJD(4) * t77 - t288, -t279, 0;];
Cq = t19;
