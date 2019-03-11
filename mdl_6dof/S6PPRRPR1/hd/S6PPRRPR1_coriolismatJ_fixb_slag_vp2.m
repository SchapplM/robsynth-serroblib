% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPRRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:38
% EndTime: 2019-03-08 18:45:50
% DurationCPUTime: 7.06s
% Computational Cost: add. (14447->539), mult. (39707->812), div. (0->0), fcn. (46580->14), ass. (0->286)
t304 = sin(qJ(4));
t307 = cos(qJ(4));
t274 = -pkin(4) * t307 - t304 * qJ(5) - pkin(3);
t302 = cos(pkin(13));
t258 = t302 * t274;
t299 = sin(pkin(13));
t383 = t302 * t304;
t195 = -pkin(10) * t383 + t258 + (-pkin(9) * t299 - pkin(5)) * t307;
t381 = t302 * t307;
t222 = pkin(9) * t381 + t299 * t274;
t386 = t299 * t304;
t208 = -pkin(10) * t386 + t222;
t303 = sin(qJ(6));
t306 = cos(qJ(6));
t123 = t195 * t306 - t208 * t303;
t124 = t195 * t303 + t208 * t306;
t280 = t304 * pkin(4) - qJ(5) * t307;
t229 = pkin(9) * t386 + t302 * t280;
t197 = t304 * pkin(5) - pkin(10) * t381 + t229;
t230 = -pkin(9) * t383 + t299 * t280;
t385 = t299 * t307;
t209 = -pkin(10) * t385 + t230;
t126 = t197 * t306 - t209 * t303;
t127 = t197 * t303 + t209 * t306;
t305 = sin(qJ(3));
t397 = sin(pkin(6));
t363 = sin(pkin(12)) * t397;
t424 = cos(qJ(3));
t301 = sin(pkin(7));
t352 = cos(pkin(12)) * t397;
t398 = cos(pkin(7));
t399 = cos(pkin(6));
t465 = t399 * t301 + t398 * t352;
t199 = t305 * t465 + t424 * t363;
t319 = -t301 * t352 + t398 * t399;
t145 = t199 * t304 - t307 * t319;
t146 = t199 * t307 + t304 * t319;
t343 = t306 * t299 + t303 * t302;
t243 = t343 * t307;
t215 = -mrSges(7,2) * t304 - mrSges(7,3) * t243;
t382 = t302 * t306;
t342 = t299 * t303 - t382;
t244 = t342 * t307;
t217 = mrSges(7,1) * t304 + mrSges(7,3) * t244;
t264 = -t304 * mrSges(6,2) - mrSges(6,3) * t385;
t266 = t304 * mrSges(6,1) - mrSges(6,3) * t381;
t372 = pkin(5) * t299 + pkin(9);
t268 = t372 * t304;
t269 = t372 * t307;
t221 = -pkin(9) * t385 + t258;
t344 = -t221 * t299 + t222 * t302;
t329 = pkin(9) * t307 - t344;
t241 = t303 * t386 - t304 * t382;
t242 = t343 * t304;
t176 = mrSges(7,1) * t242 - mrSges(7,2) * t241;
t407 = t302 * mrSges(6,2);
t410 = t299 * mrSges(6,1);
t351 = t407 + t410;
t251 = t351 * t304;
t364 = -t251 / 0.2e1 - t176 / 0.2e1;
t422 = pkin(9) * t304;
t216 = -mrSges(7,1) * t307 + t241 * mrSges(7,3);
t440 = -t216 / 0.2e1;
t214 = mrSges(7,2) * t307 - t242 * mrSges(7,3);
t442 = -t214 / 0.2e1;
t457 = -m(7) / 0.2e1;
t459 = -m(6) / 0.2e1;
t198 = t305 * t363 - t424 * t465;
t94 = -t146 * t299 + t198 * t302;
t95 = t146 * t302 + t198 * t299;
t69 = -t303 * t95 + t306 * t94;
t70 = t303 * t94 + t306 * t95;
t86 = t343 * t145;
t87 = t342 * t145;
t473 = t364 * t146 + (t145 * t329 + t146 * t422 + t229 * t94 + t230 * t95) * t459 + (t123 * t86 + t124 * t87 + t126 * t69 + t127 * t70 + t145 * t269 + t146 * t268) * t457 - t69 * t217 / 0.2e1 - t70 * t215 / 0.2e1 + t86 * t440 + t87 * t442 - t94 * t266 / 0.2e1 - t95 * t264 / 0.2e1;
t384 = t301 * t305;
t250 = t304 * t398 + t307 * t384;
t369 = t424 * t301;
t206 = -t299 * t250 - t302 * t369;
t207 = t302 * t250 - t299 * t369;
t131 = t206 * t306 - t207 * t303;
t132 = t206 * t303 + t207 * t306;
t249 = t304 * t384 - t307 * t398;
t169 = t343 * t249;
t170 = t342 * t249;
t439 = t216 / 0.2e1;
t441 = t214 / 0.2e1;
t448 = t132 / 0.2e1;
t456 = m(7) / 0.2e1;
t458 = m(6) / 0.2e1;
t403 = t307 * mrSges(5,2);
t281 = t304 * mrSges(5,1) + t403;
t471 = -t281 / 0.2e1;
t467 = t369 * t471;
t469 = t131 / 0.2e1;
t472 = -t250 * t364 + (t229 * t206 + t230 * t207 + t249 * t329 + t250 * t422) * t458 + (t123 * t169 + t124 * t170 + t126 * t131 + t127 * t132 + t249 * t269 + t250 * t268) * t456 + t217 * t469 + t215 * t448 + t169 * t439 + t170 * t441 + t206 * t266 / 0.2e1 + t207 * t264 / 0.2e1 + t467;
t252 = t351 * t307;
t263 = mrSges(6,2) * t307 - mrSges(6,3) * t386;
t265 = -t307 * mrSges(6,1) - mrSges(6,3) * t383;
t427 = t299 / 0.2e1;
t333 = t265 * t427 - t302 * t263 / 0.2e1;
t462 = t252 / 0.2e1 + t333;
t470 = t456 + t458;
t404 = t304 * mrSges(5,2);
t468 = -t307 * mrSges(5,1) - mrSges(4,1) + t404;
t359 = t304 * t369;
t466 = -t176 - t251;
t464 = -t229 * t299 + t230 * t302;
t296 = t302 ^ 2;
t461 = 2 * qJD(4);
t460 = m(5) / 0.2e1;
t455 = -mrSges(7,1) / 0.2e1;
t454 = mrSges(7,1) / 0.2e1;
t453 = -mrSges(7,2) / 0.2e1;
t452 = mrSges(7,2) / 0.2e1;
t451 = -mrSges(7,3) / 0.2e1;
t450 = mrSges(7,3) / 0.2e1;
t112 = t198 * t385 + t199 * t302;
t113 = -t198 * t381 + t199 * t299;
t79 = t112 * t306 - t113 * t303;
t449 = t79 / 0.2e1;
t447 = t145 / 0.2e1;
t175 = mrSges(7,1) * t241 + t242 * mrSges(7,2);
t446 = -t175 / 0.2e1;
t444 = -t198 / 0.2e1;
t200 = mrSges(7,1) * t343 - mrSges(7,2) * t342;
t443 = t200 / 0.2e1;
t438 = -t241 / 0.2e1;
t437 = -t242 / 0.2e1;
t436 = -t243 / 0.2e1;
t435 = -t244 / 0.2e1;
t434 = t249 / 0.2e1;
t431 = t342 / 0.2e1;
t430 = -t342 / 0.2e1;
t429 = t343 / 0.2e1;
t428 = -t299 / 0.2e1;
t426 = t302 / 0.2e1;
t425 = t304 / 0.2e1;
t298 = t307 ^ 2;
t423 = pkin(9) * t298;
t421 = pkin(10) + qJ(5);
t420 = Ifges(6,4) * t299;
t419 = Ifges(6,4) * t302;
t418 = Ifges(7,4) * t241;
t417 = Ifges(7,4) * t343;
t416 = Ifges(6,5) * t302;
t415 = Ifges(6,6) * t299;
t414 = t243 * mrSges(7,1);
t413 = t244 * mrSges(7,2);
t412 = t342 * mrSges(7,3);
t411 = t343 * mrSges(7,3);
t409 = t299 * Ifges(6,2);
t408 = t299 * t95;
t406 = t302 * t94;
t402 = t69 * t242;
t401 = t70 * t241;
t276 = -mrSges(6,1) * t302 + mrSges(6,2) * t299;
t400 = t276 - mrSges(5,1);
t396 = mrSges(7,3) * qJD(4);
t395 = t112 * t299;
t394 = t113 * t302;
t393 = t198 * t304;
t392 = t206 * t302;
t391 = t207 * t299;
t370 = t307 * t424;
t219 = (-t299 * t370 + t302 * t305) * t301;
t390 = t219 * t299;
t220 = (t299 * t305 + t302 * t370) * t301;
t389 = t220 * t302;
t380 = -Ifges(7,5) * t242 + Ifges(7,6) * t241;
t379 = -Ifges(7,5) * t342 - Ifges(7,6) * t343;
t378 = t299 ^ 2 + t296;
t375 = m(6) * t425;
t374 = -t412 / 0.2e1;
t373 = -t411 / 0.2e1;
t371 = t304 * t424;
t368 = t145 * t443;
t367 = t200 * t434;
t361 = qJ(5) * t378;
t297 = t304 ^ 2;
t360 = t297 * t369;
t358 = t249 * t371;
t356 = t369 / 0.2e1;
t201 = mrSges(7,1) * t342 + mrSges(7,2) * t343;
t355 = qJD(4) * (t201 + t400);
t354 = t146 * t470;
t353 = t250 * t470;
t350 = t299 * t94 - t302 * t95;
t100 = t145 * t393;
t80 = t112 * t303 + t113 * t306;
t13 = m(7) * (t69 * t79 + t70 * t80 - t100) + m(6) * (t112 * t94 + t113 * t95 - t100) + m(5) * (-t145 * t304 - t146 * t307 + t199) * t198;
t155 = t219 * t306 - t220 * t303;
t156 = t219 * t303 + t220 * t306;
t325 = (t145 * t369 - t198 * t249) * t304;
t8 = (t131 * t79 + t132 * t80 + t155 * t69 + t156 * t70 + t325) * t456 + (t206 * t112 + t207 * t113 + t219 * t94 + t220 * t95 + t325) * t458 + ((-t249 * t304 - t250 * t307) * t198 + (t145 * t371 + t146 * t370 + t198 * t305 - t199 * t424) * t301) * t460;
t349 = t13 * qJD(1) + t8 * qJD(2);
t88 = t145 * t146;
t16 = m(7) * (t69 * t86 + t70 * t87 + t88) + m(6) * (t145 * t350 + t88);
t345 = t206 * t299 - t207 * t302;
t9 = (t131 * t86 + t132 * t87 + t145 * t250 + t146 * t249 + t169 * t69 + t170 * t70) * t456 + ((t146 + t350) * t249 + (t250 + t345) * t145) * t458;
t348 = t16 * qJD(1) + t9 * qJD(2);
t196 = t249 * t250;
t40 = m(7) * (t131 * t169 + t132 * t170 + t196) + m(6) * (t249 * t345 + t196);
t347 = t9 * qJD(1) + t40 * qJD(2);
t218 = t301 * t358;
t41 = m(7) * (t131 * t155 + t132 * t156 + t218) + m(5) * (t250 * t370 - t305 * t369 + t358) * t301 + m(6) * (t206 * t219 + t207 * t220 + t218);
t346 = t8 * qJD(1) + t41 * qJD(2);
t340 = qJD(3) * t175 - qJD(4) * t200;
t339 = mrSges(7,1) * t449 + t453 * t80;
t338 = t452 * t87 + t455 * t86;
t337 = Ifges(7,5) * t244 / 0.2e1 + Ifges(7,6) * t243 / 0.2e1;
t336 = t155 * t454 + t156 * t453;
t335 = t169 * t454 + t170 * t453;
t334 = m(7) * (t241 * t69 - t242 * t70);
t332 = m(7) * (t131 * t241 - t132 * t242);
t331 = qJD(4) * (-t378 * mrSges(6,3) + mrSges(5,2));
t165 = -Ifges(7,2) * t242 - t307 * Ifges(7,6) - t418;
t166 = -Ifges(7,4) * t244 - Ifges(7,2) * t243 + Ifges(7,6) * t304;
t238 = Ifges(7,4) * t242;
t167 = -Ifges(7,1) * t241 - t307 * Ifges(7,5) - t238;
t168 = -Ifges(7,1) * t244 - Ifges(7,4) * t243 + Ifges(7,5) * t304;
t177 = -t413 + t414;
t239 = Ifges(6,6) * t304 + (-t409 + t419) * t307;
t240 = Ifges(6,5) * t304 + (t302 * Ifges(6,1) - t420) * t307;
t12 = -pkin(3) * t281 + t268 * t177 + t269 * t176 + t230 * t263 + t222 * t264 + t229 * t265 + t221 * t266 + t166 * t437 + t165 * t436 + t167 * t435 + t168 * t438 + t127 * t214 + t124 * t215 + t126 * t216 + t123 * t217 + m(7) * (t123 * t126 + t124 * t127 + t268 * t269) + m(6) * (t221 * t229 + t222 * t230) + (Ifges(7,5) * t438 + Ifges(7,6) * t437 + pkin(9) * t252 + t239 * t428 + t240 * t426 + (-Ifges(5,4) + t416 / 0.2e1 - t415 / 0.2e1) * t304) * t304 + (pkin(9) * t251 + (-Ifges(5,2) + Ifges(5,1) - Ifges(6,3) - Ifges(7,3) + m(6) * pkin(9) ^ 2 + t296 * Ifges(6,1) / 0.2e1 + (-t419 + t409 / 0.2e1) * t299) * t304 + t337 + (Ifges(5,4) + t415 - t416) * t307) * t307;
t275 = t421 * t299;
t277 = t421 * t302;
t212 = -t275 * t306 - t277 * t303;
t213 = -t275 * t303 + t277 * t306;
t292 = -pkin(5) * t302 - pkin(4);
t312 = (-t395 / 0.2e1 + t394 / 0.2e1) * mrSges(6,3) + (-t343 * t449 + t430 * t80) * mrSges(7,3) + (pkin(4) * t393 + (t394 - t395) * qJ(5)) * t458 + (t212 * t79 + t213 * t80 - t292 * t393) * t456;
t322 = t403 / 0.2e1 + (mrSges(5,1) / 0.2e1 - t276 / 0.2e1 - t201 / 0.2e1) * t304;
t323 = -t177 / 0.2e1 - t462;
t2 = t323 * t145 + (t471 + t322) * t198 + t312 + t473;
t308 = (-pkin(4) * t359 + (t389 - t390) * qJ(5)) * t458 + (t212 * t155 + t213 * t156 + t292 * t359) * t456 + t155 * t373 + t156 * t374 + t467 + (t201 + t276) * t304 * t356 + (-t390 / 0.2e1 + t389 / 0.2e1) * mrSges(6,3);
t6 = t177 * t434 + t249 * t462 - t308 + t472;
t328 = -t2 * qJD(1) + t6 * qJD(2) + t12 * qJD(3);
t315 = (t241 * t448 + t242 * t469) * mrSges(7,3) + t131 * t441 + t132 * t440 + t249 * t446;
t17 = t315 - t336;
t178 = Ifges(7,2) * t241 - t238;
t179 = -Ifges(7,1) * t242 + t418;
t24 = -t268 * t175 - t307 * t380 / 0.2e1 + t123 * t214 - t124 * t216 - (t167 / 0.2e1 + t178 / 0.2e1 - t123 * mrSges(7,3)) * t242 + (-t179 / 0.2e1 + t165 / 0.2e1 + t124 * mrSges(7,3)) * t241;
t324 = t175 * t447 + t439 * t70 + t442 * t69;
t4 = (-t402 / 0.2e1 - t401 / 0.2e1) * mrSges(7,3) + t324 + t339;
t327 = -t4 * qJD(1) + t17 * qJD(2) + t24 * qJD(3);
t29 = -t334 / 0.2e1 + (m(7) * t444 + (t444 + t408 / 0.2e1 + t406 / 0.2e1) * m(6)) * t304;
t48 = -t242 * t214 + t241 * t216 + m(7) * (t123 * t241 - t124 * t242) + (-t299 * t263 - t302 * t265 + m(6) * (-t221 * t302 - t222 * t299)) * t304;
t53 = -t332 / 0.2e1 + (m(7) * t356 + (t356 + t392 / 0.2e1 + t391 / 0.2e1) * m(6)) * t304;
t326 = -qJD(1) * t29 - qJD(2) * t53 + qJD(3) * t48;
t11 = t368 + t338;
t256 = Ifges(7,4) * t342;
t202 = -Ifges(7,2) * t343 - t256;
t203 = -Ifges(7,2) * t342 + t417;
t204 = -Ifges(7,1) * t342 - t417;
t205 = Ifges(7,1) * t343 - t256;
t311 = -(t167 / 0.4e1 + t178 / 0.4e1) * t342 - (-t179 / 0.4e1 + t165 / 0.4e1) * t343 - (t212 * t451 + t205 / 0.4e1 + t202 / 0.4e1) * t242 + (t213 * t450 - t204 / 0.4e1 + t203 / 0.4e1) * t241 + t212 * t441 + t213 * t440 + t268 * t443 + t292 * t446 - t307 * t379 / 0.4e1;
t317 = -Ifges(7,3) * t304 / 0.2e1 + t126 * t455 + t127 * t452 + t337;
t15 = t311 + t317;
t31 = t367 - t335;
t44 = t292 * t200 - (t203 / 0.2e1 - t204 / 0.2e1) * t343 - (t202 / 0.2e1 + t205 / 0.2e1) * t342;
t321 = t11 * qJD(1) + t31 * qJD(2) + t15 * qJD(3) + t44 * qJD(4);
t318 = t350 * t458 + (-t342 * t70 - t343 * t69) * t457;
t27 = t354 + t318;
t313 = t344 * t459 + (-t123 * t343 - t124 * t342 + t212 * t241 - t213 * t242) * t457 + t214 * t431 + t216 * t429 + t333;
t314 = (pkin(9) * t458 + t407 / 0.2e1 + t410 / 0.2e1) * t307 + t269 * t456 + t414 / 0.2e1 - t413 / 0.2e1;
t33 = (t241 * t429 - t242 * t431) * mrSges(7,3) + t313 + t314;
t316 = t345 * t458 + (-t131 * t343 - t132 * t342) * t457;
t51 = t353 + t316;
t78 = (t342 ^ 2 + t343 ^ 2) * mrSges(7,3) + m(7) * (-t212 * t343 - t213 * t342) + (m(6) * qJ(5) + mrSges(6,3)) * t378;
t320 = -qJD(1) * t27 - qJD(2) * t51 - qJD(3) * t33 + qJD(4) * t78;
t278 = pkin(9) * t360;
t187 = t297 * pkin(9) * t198;
t54 = t332 / 0.2e1 + (-t391 - t392) * t375 + (m(6) + m(7)) * t359 / 0.2e1;
t52 = t353 - t316;
t34 = t241 * t373 - t242 * t374 - t313 + t314;
t32 = t367 + t335;
t30 = t334 / 0.2e1 + (-t406 - t408) * t375 + (t457 + t459) * t393;
t28 = t354 - t318;
t18 = t315 + t336;
t14 = t311 - t317;
t10 = t368 - t338;
t7 = -t249 * t323 + t308 + t472;
t5 = t401 * t450 - t402 * t451 - t324 + t339;
t3 = t177 * t447 + t312 + (t281 / 0.2e1 + t322) * t198 + t462 * t145 - t473;
t1 = t8 * qJD(3) + t9 * qJD(4);
t19 = [t13 * qJD(3) + t16 * qJD(4), t1 (t112 * t265 + t113 * t263 + t468 * t199 + t80 * t214 + t79 * t216) * qJD(3) + t3 * qJD(4) + t30 * qJD(5) + t5 * qJD(6) + (mrSges(4,2) + t466 * t304 + (-t297 - t298) * mrSges(5,3)) * qJD(3) * t198 + 0.2e1 * ((t123 * t79 + t124 * t80 - t268 * t393) * t456 + (t112 * t221 + t113 * t222 - t187) * t458 + (-pkin(3) * t199 - t198 * t423 - t187) * t460) * qJD(3) + t349, t3 * qJD(3) + t28 * qJD(5) + t10 * qJD(6) + (-t342 * t87 - t343 * t86) * t396 + t146 * t355 + t145 * t331 + ((t146 * t292 + t212 * t86 + t213 * t87) * t456 + (-pkin(4) * t146 - t145 * t361) * t458) * t461 + t348, qJD(3) * t30 + qJD(4) * t28, t5 * qJD(3) + t10 * qJD(4) + (-mrSges(7,1) * t70 - mrSges(7,2) * t69) * qJD(6); t1, t41 * qJD(3) + t40 * qJD(4) (-mrSges(4,2) * t369 + m(7) * (t123 * t155 + t124 * t156) + t156 * t214 + t155 * t216 + m(5) * (t278 + (-pkin(3) * t305 + t423 * t424) * t301) + m(6) * (t219 * t221 + t220 * t222 + t278) + t219 * t265 + t220 * t263 + t468 * t384 + (m(7) * t268 - t466) * t359 + (t298 * t369 + t360) * mrSges(5,3)) * qJD(3) + t7 * qJD(4) + t54 * qJD(5) + t18 * qJD(6) + t346, t7 * qJD(3) + t52 * qJD(5) + t32 * qJD(6) + (-t169 * t343 - t170 * t342) * t396 + t250 * t355 + t249 * t331 + ((t169 * t212 + t170 * t213 + t250 * t292) * t456 + (-pkin(4) * t250 - t249 * t361) * t458) * t461 + t347, qJD(3) * t54 + qJD(4) * t52, t18 * qJD(3) + t32 * qJD(4) + (-mrSges(7,1) * t132 - mrSges(7,2) * t131) * qJD(6); -qJD(4) * t2 - qJD(5) * t29 - qJD(6) * t4 - t349, qJD(4) * t6 - qJD(5) * t53 + qJD(6) * t17 - t346, qJD(4) * t12 + qJD(5) * t48 + qJD(6) * t24, t34 * qJD(5) + t14 * qJD(6) + t328 + (-Ifges(5,6) * t304 + t239 * t426 + t240 * t427 + t292 * t177 + t269 * t201 + t168 * t429 + t166 * t430 - pkin(4) * t252 + t203 * t436 + t205 * t435 + t213 * t215 + t212 * t217 - t127 * t412 - t126 * t411 + m(7) * (t126 * t212 + t127 * t213 + t269 * t292) + pkin(9) * t404 + (m(6) * t464 + t302 * t264 - t299 * t266) * qJ(5) + (Ifges(5,5) + (Ifges(6,2) * t302 + t420) * t428 + (Ifges(6,1) * t299 + t419) * t426 + (-m(6) * pkin(4) + t400) * pkin(9)) * t307 + (Ifges(6,5) * t299 + Ifges(7,5) * t343 + Ifges(6,6) * t302 - Ifges(7,6) * t342) * t425 + t464 * mrSges(6,3)) * qJD(4), qJD(4) * t34 + t326, t14 * qJD(4) + (-mrSges(7,1) * t124 - mrSges(7,2) * t123 + t380) * qJD(6) + t327; qJD(3) * t2 - qJD(5) * t27 + qJD(6) * t11 - t348, -qJD(3) * t6 - qJD(5) * t51 + qJD(6) * t31 - t347, -qJD(5) * t33 + qJD(6) * t15 - t328, qJD(5) * t78 + qJD(6) * t44, t320 (-mrSges(7,1) * t213 - mrSges(7,2) * t212 + t379) * qJD(6) + t321; qJD(3) * t29 + qJD(4) * t27, qJD(3) * t53 + qJD(4) * t51, qJD(4) * t33 - qJD(6) * t175 - t326, qJD(6) * t200 - t320, 0, -t340; t4 * qJD(3) - t11 * qJD(4), -t17 * qJD(3) - t31 * qJD(4), -qJD(4) * t15 + qJD(5) * t175 - t327, -qJD(5) * t200 - t321, t340, 0;];
Cq  = t19;
