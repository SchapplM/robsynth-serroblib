% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:20
% EndTime: 2019-03-08 19:55:28
% DurationCPUTime: 4.84s
% Computational Cost: add. (7536->468), mult. (19149->660), div. (0->0), fcn. (18979->10), ass. (0->261)
t293 = cos(qJ(5));
t468 = Ifges(6,6) + Ifges(7,6);
t478 = t293 * t468;
t289 = sin(pkin(11));
t395 = sin(pkin(6));
t396 = cos(pkin(11));
t334 = t396 * t395;
t292 = sin(qJ(2));
t343 = t292 * t395;
t427 = cos(qJ(2));
t178 = t289 * t343 - t334 * t427;
t290 = sin(qJ(5));
t417 = -qJ(6) - pkin(9);
t236 = t417 * t290;
t239 = t417 * t293;
t335 = t427 * t395;
t179 = t289 * t335 + t292 * t334;
t294 = cos(qJ(4));
t378 = t293 * t294;
t89 = -t178 * t378 + t179 * t290;
t399 = t89 * t293;
t381 = t290 * t294;
t88 = t178 * t381 + t179 * t293;
t400 = t88 * t290;
t328 = t399 - t400;
t421 = pkin(5) * t293;
t273 = -pkin(4) - t421;
t291 = sin(qJ(4));
t383 = t273 * t291;
t387 = t178 * t291;
t453 = -m(7) / 0.2e1;
t470 = mrSges(7,3) + mrSges(6,3);
t477 = -t470 * (t399 / 0.2e1 - t400 / 0.2e1) - m(6) * (pkin(4) * t387 + pkin(9) * t328) / 0.2e1 + (-t178 * t383 + t236 * t88 - t239 * t89) * t453;
t455 = m(6) / 0.2e1;
t452 = m(7) / 0.2e1;
t465 = t290 * mrSges(7,1) + t293 * mrSges(7,2);
t211 = t465 * t291;
t271 = pkin(2) * t289 + pkin(8);
t422 = pkin(5) * t290;
t344 = t271 + t422;
t196 = t344 * t291;
t426 = m(7) * t196;
t476 = t211 + t426;
t411 = Ifges(7,4) * t290;
t248 = Ifges(7,1) * t293 - t411;
t184 = -t294 * Ifges(7,5) + t291 * t248;
t412 = Ifges(6,4) * t290;
t250 = Ifges(6,1) * t293 - t412;
t186 = -t294 * Ifges(6,5) + t291 * t250;
t366 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t338 = t366 * t294;
t475 = -t338 + t184 / 0.2e1 + t186 / 0.2e1;
t241 = t290 * mrSges(6,1) + t293 * mrSges(6,2);
t449 = mrSges(7,1) / 0.2e1;
t369 = mrSges(6,1) / 0.2e1 + t449;
t447 = mrSges(7,2) / 0.2e1;
t448 = mrSges(6,2) / 0.2e1;
t368 = t448 + t447;
t466 = t368 * t293;
t474 = t466 + t290 * t369 - t465 / 0.2e1 - t241 / 0.2e1;
t372 = m(6) / 0.4e1 + m(7) / 0.4e1;
t472 = 0.4e1 * t372;
t471 = mrSges(6,2) + mrSges(7,2);
t469 = Ifges(6,5) + Ifges(7,5);
t467 = Ifges(7,3) + Ifges(6,3);
t450 = m(7) * pkin(5);
t364 = mrSges(7,1) + t450;
t282 = Ifges(7,4) * t293;
t464 = -Ifges(7,2) * t290 + t282;
t247 = Ifges(7,1) * t290 + t282;
t283 = Ifges(6,4) * t293;
t463 = -Ifges(6,2) * t290 + t283;
t249 = Ifges(6,1) * t290 + t283;
t419 = t291 * pkin(4);
t420 = pkin(9) * t294;
t254 = t419 - t420;
t382 = t290 * t291;
t137 = t293 * t254 + t271 * t382;
t379 = t291 * t293;
t138 = t290 * t254 - t271 * t379;
t461 = -t137 * t290 + t138 * t293;
t113 = t291 * pkin(5) - qJ(6) * t378 + t137;
t127 = -qJ(6) * t381 + t138;
t460 = -t113 * t290 + t127 * t293;
t213 = t465 * t294;
t214 = t294 * t241;
t403 = t294 * mrSges(7,1);
t231 = -mrSges(7,3) * t379 - t403;
t232 = -mrSges(6,1) * t294 - mrSges(6,3) * t379;
t351 = t232 / 0.2e1 + t231 / 0.2e1;
t228 = t294 * mrSges(6,2) - mrSges(6,3) * t382;
t401 = t294 * mrSges(7,2);
t227 = -mrSges(7,3) * t382 + t401;
t440 = t227 / 0.2e1;
t353 = t440 + t228 / 0.2e1;
t303 = t290 * t351 - t293 * t353 + t213 / 0.2e1 + t214 / 0.2e1;
t397 = cos(pkin(6));
t131 = t179 * t291 - t294 * t397;
t132 = t179 * t294 + t291 * t397;
t272 = -pkin(2) * t396 - pkin(3);
t219 = -t294 * pkin(4) - t291 * pkin(9) + t272;
t190 = t293 * t219;
t363 = qJ(6) * t379;
t105 = -t363 + t190 + (-t271 * t290 - pkin(5)) * t294;
t130 = t290 * t219 + t271 * t378;
t112 = -qJ(6) * t382 + t130;
t197 = t344 * t294;
t317 = t105 * t290 - t112 * t293 + t197;
t129 = -t271 * t381 + t190;
t324 = -t129 * t290 + t130 * t293;
t234 = t291 * mrSges(6,1) - mrSges(6,3) * t378;
t233 = t291 * mrSges(7,1) - mrSges(7,3) * t378;
t437 = t233 / 0.2e1;
t350 = t437 + t234 / 0.2e1;
t229 = -t291 * mrSges(7,2) - mrSges(7,3) * t381;
t230 = -t291 * mrSges(6,2) - mrSges(6,3) * t381;
t352 = t229 / 0.2e1 + t230 / 0.2e1;
t212 = t241 * t291;
t441 = t211 / 0.2e1;
t355 = t212 / 0.2e1 + t441;
t384 = t271 * t291;
t77 = -t132 * t290 + t178 * t293;
t78 = t132 * t293 + t178 * t290;
t457 = t355 * t132 + t352 * t78 + t350 * t77 + (t132 * t384 + t137 * t77 + t138 * t78) * t455 + (t113 * t77 + t127 * t78 + t132 * t196) * t452 + ((t271 * t294 - t324) * t455 + t317 * t452 + t303) * t131;
t456 = 2 * qJD(4);
t446 = -mrSges(7,3) / 0.2e1;
t443 = -t113 / 0.2e1;
t408 = t290 * mrSges(6,2);
t415 = mrSges(6,1) * t293;
t333 = -t408 + t415;
t210 = t333 * t291;
t442 = -t210 / 0.2e1;
t439 = -t231 / 0.2e1;
t438 = -t232 / 0.2e1;
t436 = t236 / 0.2e1;
t237 = -mrSges(7,1) * t293 + mrSges(7,2) * t290;
t435 = t237 / 0.2e1;
t434 = t465 / 0.2e1;
t402 = t294 * mrSges(5,2);
t242 = t291 * mrSges(5,1) + t402;
t433 = t242 / 0.2e1;
t432 = t271 / 0.2e1;
t430 = t291 / 0.2e1;
t429 = t293 / 0.2e1;
t425 = m(7) * t197;
t424 = m(7) * t273;
t423 = m(7) * t291;
t286 = t291 ^ 2;
t288 = t294 ^ 2;
t316 = m(6) * t461;
t16 = m(6) * (t286 - t288) * t432 + (t317 * t453 + t324 * t455 - t303) * t294 + (t352 * t293 - t350 * t290 + (t196 + t460) * t452 + t316 / 0.2e1 + t355) * t291;
t285 = t290 ^ 2;
t287 = t293 ^ 2;
t346 = -t287 / 0.2e1 - t285 / 0.2e1;
t336 = t346 * mrSges(6,3);
t269 = mrSges(7,2) * t382;
t209 = mrSges(7,1) * t379 - t269;
t356 = t442 - t209 / 0.2e1;
t111 = t129 - t363;
t377 = -t105 + t111;
t21 = t356 * t294 + (mrSges(7,3) * t346 + t336) * t286 + (-t353 * t290 + ((-pkin(5) * t294 + t377) * t452 - t351) * t293) * t291;
t416 = t16 * qJD(4) + t21 * qJD(5);
t407 = t290 * t77;
t406 = t290 * t78;
t405 = t293 * t77;
t404 = t293 * t78;
t398 = -t333 - mrSges(5,1);
t386 = t178 * t294;
t25 = 0.2e1 * t372 * (t328 + t386) * t291;
t394 = qJD(1) * t25;
t389 = t178 * t286;
t388 = t178 * t288;
t385 = t21 * qJD(2);
t380 = t291 * t131;
t376 = t285 + t287;
t375 = qJD(4) * t294;
t374 = qJD(5) * t291;
t373 = t450 / 0.2e1;
t371 = -0.1e1 + t376;
t370 = t423 / 0.2e1;
t367 = t446 - mrSges(6,3) / 0.2e1;
t365 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t360 = -t387 / 0.2e1;
t359 = -t381 / 0.2e1;
t358 = t378 / 0.2e1;
t243 = Ifges(7,2) * t293 + t411;
t245 = Ifges(6,2) * t293 + t412;
t348 = -t245 / 0.2e1 - t243 / 0.2e1;
t347 = t247 / 0.2e1 + t249 / 0.2e1;
t345 = m(7) * t377;
t342 = t376 * mrSges(7,3);
t341 = -mrSges(6,1) - t364;
t337 = qJD(4) * (t237 + t398);
t329 = t404 - t407;
t41 = (-t178 / 0.2e1 + t406 / 0.2e1 + t405 / 0.2e1) * t423;
t47 = (m(7) * (-t105 * t293 - t112 * t290) - t290 * t227 - t293 * t231) * t291;
t327 = -qJD(1) * t41 + qJD(2) * t47;
t12 = (-t178 * t380 + t77 * t88 + t78 * t89) * t472 + m(5) * (-t132 * t294 + t179 - t380) * t178;
t326 = t12 * qJD(1) + t25 * qJD(3);
t320 = t132 - t329;
t13 = t320 * t131 * t472;
t15 = -0.2e1 * (t320 * t294 + t371 * t380) * t372;
t325 = t13 * qJD(1) + t15 * qJD(3);
t323 = t236 * t290 + t239 * t293;
t322 = t373 + t369;
t150 = t364 * t379 - t269;
t220 = -m(7) * t422 - t465;
t321 = qJD(2) * t150 - qJD(4) * t220;
t319 = t439 + t345 / 0.2e1;
t1 = (-t402 / 0.2e1 + t433 + (-mrSges(5,1) / 0.2e1 - t333 / 0.2e1 + t435) * t291) * t178 + t457 + t477;
t181 = Ifges(7,6) * t291 + t294 * t464;
t183 = Ifges(6,6) * t291 + t294 * t463;
t185 = Ifges(7,5) * t291 + t248 * t294;
t187 = Ifges(6,5) * t291 + t250 * t294;
t180 = -t294 * Ifges(7,6) + t291 * t464;
t182 = -t294 * Ifges(6,6) + t291 * t463;
t310 = -t180 / 0.2e1 - t182 / 0.2e1 + t365 * t294;
t7 = t105 * t233 + t112 * t229 + t113 * t231 + t127 * t227 + t129 * t234 + t130 * t230 + t137 * t232 + t138 * t228 + t196 * t213 + t197 * t211 + t272 * t242 + m(6) * (t129 * t137 + t130 * t138) + m(7) * (t105 * t113 + t112 * t127 + t196 * t197) + (Ifges(5,4) * t294 + t271 * t212 + t310 * t290 + t475 * t293) * t294 + (-Ifges(5,4) * t291 + t271 * t214 + (t185 / 0.2e1 + t187 / 0.2e1 + t366 * t291) * t293 + (-t181 / 0.2e1 - t183 / 0.2e1 - t365 * t291) * t290 + (m(6) * t271 ^ 2 + Ifges(5,1) - Ifges(5,2) - t467) * t294) * t291;
t314 = t1 * qJD(1) + t7 * qJD(2) + t16 * qJD(3);
t295 = (t370 * t421 - t356) * t131 + (-t367 * t382 + t353) * t77 + (t367 * t379 + t319 + t438) * t78;
t3 = -t322 * t88 + t368 * t89 + t295;
t215 = t291 * t243;
t216 = t291 * t245;
t217 = t291 * t247;
t218 = t291 * t249;
t9 = t111 * t227 + t129 * t228 - t130 * t232 + t196 * t209 + (t345 - t231) * t112 + (t271 * t210 + (t105 * mrSges(7,3) + t129 * mrSges(6,3) + t216 / 0.2e1 + t215 / 0.2e1 - t475) * t290 + (-t112 * mrSges(7,3) - t130 * mrSges(6,3) - t217 / 0.2e1 - t218 / 0.2e1 + t476 * pkin(5) + t310) * t293) * t291;
t313 = t3 * qJD(1) + t9 * qJD(2) + t21 * qJD(3);
t76 = t371 * t294 * t291 * t472;
t311 = t15 * qJD(1) + t16 * qJD(2) + t76 * qJD(3);
t308 = pkin(9) * t336 + t241 * t432;
t307 = -t249 / 0.4e1 - t247 / 0.4e1 - t463 / 0.4e1 - t464 / 0.4e1 + mrSges(7,3) * t436;
t10 = t474 * t131;
t32 = -pkin(4) * t241 + t273 * t465 + (t463 / 0.2e1 + t464 / 0.2e1 + t347) * t293 + (t248 / 0.2e1 + t250 / 0.2e1 + (t237 + t424) * pkin(5) + t348) * t290;
t280 = Ifges(7,5) * t293;
t281 = Ifges(6,5) * t293;
t297 = -t319 * t239 + (-t281 / 0.4e1 - t280 / 0.4e1) * t294 + pkin(4) * t442 + t196 * t434 + t227 * t436 + t273 * t209 / 0.2e1;
t299 = t250 / 0.4e1 + t248 / 0.4e1 - t245 / 0.4e1 - t243 / 0.4e1 - t239 * t446 + (t424 / 0.2e1 + t435) * pkin(5);
t300 = (t441 + t426 / 0.2e1) * pkin(5) - t180 / 0.4e1 - t182 / 0.4e1 - t217 / 0.4e1 - t218 / 0.4e1 - pkin(9) * t228 / 0.2e1;
t301 = pkin(9) * t438 - t216 / 0.4e1 - t215 / 0.4e1 + t186 / 0.4e1 + t184 / 0.4e1 + (t111 / 0.2e1 - t105 / 0.2e1) * mrSges(7,3);
t304 = mrSges(7,1) * t443 + t127 * t447 - t137 * mrSges(6,1) / 0.2e1 + t138 * t448;
t5 = (m(7) * t443 - t233 / 0.2e1) * pkin(5) + (-Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t308) * t291 + (t291 * t299 + t301 - t338) * t293 + ((0.3e1 / 0.4e1 * Ifges(6,6) + 0.3e1 / 0.4e1 * Ifges(7,6)) * t294 + t307 * t291 + t300) * t290 + t297 + t304;
t54 = t474 * t294;
t306 = -t10 * qJD(1) + t5 * qJD(2) + t54 * qJD(3) + t32 * qJD(4);
t109 = -m(7) * t323 + t342;
t152 = (-0.1e1 / 0.2e1 - t346) * t423;
t302 = ((-t236 * t291 + t112) * t293 + (t239 * t291 - t105) * t290) * t452;
t37 = (t440 - t401 / 0.2e1) * t293 + (t439 - t403 / 0.2e1) * t290 + t302 - t425 / 0.2e1;
t40 = 0.2e1 * (-t407 / 0.4e1 + t404 / 0.4e1 - t132 / 0.4e1) * m(7);
t305 = qJD(1) * t40 + qJD(2) * t37 + qJD(3) * t152 + qJD(4) * t109;
t151 = (t376 + 0.1e1) * t370;
t133 = t271 * t389;
t55 = t359 * t450 + (-t290 * t322 - t466) * t294 - (t241 + t465) * t294 / 0.2e1;
t42 = (-t405 - t406) * t370 + m(7) * t360;
t39 = (t132 + t329) * t452;
t36 = t302 + t290 * t439 + t227 * t429 + t425 / 0.2e1 + mrSges(7,2) * t358 + t381 * t449;
t11 = (t434 + t241 / 0.2e1 + t466 + (t373 + t322) * t290) * t131;
t8 = qJD(2) * t25 + qJD(4) * t15;
t6 = t297 + t113 * t373 + t301 * t293 + (t307 * t290 + t299 * t293 + t308) * t291 + ((Ifges(6,6) / 0.4e1 + Ifges(7,6) / 0.4e1) * t294 + t300) * t290 + pkin(5) * t437 - t304 + t467 * t430 + t468 * t359 + t469 * t358;
t4 = t373 * t88 + t295 + (mrSges(6,1) + mrSges(7,1)) * t88 / 0.2e1 - t471 * t89 / 0.2e1;
t2 = mrSges(5,1) * t387 / 0.2e1 + mrSges(5,2) * t386 / 0.2e1 + t178 * t433 + (-t333 + t237) * t360 + t457 - t477;
t14 = [qJD(2) * t12 + qJD(4) * t13 (-mrSges(3,2) * t335 - mrSges(3,1) * t343 - m(6) * t133 + m(5) * (t179 * t272 - t271 * t388 - t133) + t179 * (-t294 * mrSges(5,1) + t291 * mrSges(5,2)) + m(4) * (-t178 * t289 - t179 * t396) * pkin(2) + t178 * mrSges(4,2) - t179 * mrSges(4,1) + (m(6) * t130 + m(7) * t112 + t227 + t228) * t89 + (m(6) * t129 + m(7) * t105 + t231 + t232) * t88 + (-t212 - t476) * t387 + (-t388 - t389) * mrSges(5,3)) * qJD(2) + t2 * qJD(4) + t4 * qJD(5) + t42 * qJD(6) + t326, t8, t2 * qJD(2) + t11 * qJD(5) + t39 * qJD(6) + t132 * t337 + (-t376 * t470 + mrSges(5,2)) * qJD(4) * t131 + ((-pkin(9) * t131 * t376 - pkin(4) * t132) * t455 + (t131 * t323 + t132 * t273) * t452) * t456 + t325, t4 * qJD(2) + t11 * qJD(4) + (t341 * t78 - t471 * t77) * qJD(5), qJD(2) * t42 + qJD(4) * t39; qJD(4) * t1 + qJD(5) * t3 - qJD(6) * t41 - t326, qJD(4) * t7 + qJD(5) * t9 + qJD(6) * t47, -t394 + t416, t6 * qJD(5) + t36 * qJD(6) + (Ifges(5,5) + t347 * t293 + t348 * t290 + (-m(6) * pkin(4) + t398) * t271) * t375 + t314 + (mrSges(5,2) * t384 + m(7) * (t113 * t236 - t127 * t239 + t197 * t273) - Ifges(5,6) * t291 + t273 * t213 - t239 * t229 + t236 * t233 + t197 * t237 - pkin(4) * t214 + (t293 * t230 - t290 * t234 + t316) * pkin(9) + (t185 + t187) * t290 / 0.2e1 + (t290 * t469 + t478) * t430 + (t181 + t183) * t429 + t460 * mrSges(7,3) + t461 * mrSges(6,3)) * qJD(4), t6 * qJD(4) + (-mrSges(6,1) * t130 - mrSges(6,2) * t129 - mrSges(7,2) * t111 - t112 * t364) * qJD(5) + (-t478 + (mrSges(7,3) * pkin(5) - t469) * t290) * t374 + t313, qJD(4) * t36 + t327; t8, t394 + t416, t76 * qJD(4), t55 * qJD(5) + t151 * qJD(6) + t291 * t337 + (mrSges(6,3) * t376 - mrSges(5,2) + t342) * t375 + ((t376 * t420 - t419) * t455 + (-t294 * t323 + t383) * t452) * t456 + t311, t385 + t55 * qJD(4) + t269 * qJD(5) + (t293 * t341 + t408) * t374, t151 * qJD(4); -qJD(2) * t1 - qJD(5) * t10 + qJD(6) * t40 - t325, qJD(5) * t5 + qJD(6) * t37 - t314, qJD(5) * t54 + qJD(6) * t152 - t311, qJD(5) * t32 + qJD(6) * t109, t306 + (-mrSges(7,2) * t236 - mrSges(7,3) * t421 - pkin(9) * t415 + t280 + t281 + (mrSges(6,2) * pkin(9) - t468) * t290 + t364 * t239) * qJD(5), t305; -qJD(2) * t3 + qJD(4) * t10, -qJD(4) * t5 - qJD(6) * t150 - t313, -qJD(4) * t54 - t385, qJD(6) * t220 - t306, 0, -t321; qJD(2) * t41 - qJD(4) * t40, -qJD(4) * t37 + qJD(5) * t150 - t327, -t152 * qJD(4), -qJD(5) * t220 - t305, t321, 0;];
Cq  = t14;
