% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:04
% EndTime: 2019-03-09 02:18:13
% DurationCPUTime: 6.00s
% Computational Cost: add. (20797->405), mult. (39135->538), div. (0->0), fcn. (45701->10), ass. (0->248)
t251 = cos(qJ(6));
t246 = sin(pkin(11));
t247 = cos(pkin(11));
t406 = sin(qJ(4));
t408 = cos(qJ(4));
t225 = -t246 * t406 + t247 * t408;
t226 = t246 * t408 + t247 * t406;
t250 = sin(qJ(5));
t407 = cos(qJ(5));
t274 = t250 * t225 + t226 * t407;
t210 = -t407 * t225 + t226 * t250;
t249 = sin(qJ(6));
t463 = t249 * t210;
t479 = -mrSges(7,2) * t274 + mrSges(7,3) * t463;
t497 = t251 * t479;
t462 = t251 * t210;
t478 = mrSges(7,1) * t274 + mrSges(7,3) * t462;
t498 = t249 * t478;
t500 = t497 - t498;
t499 = t497 / 0.2e1 - t498 / 0.2e1;
t376 = t251 * mrSges(7,2);
t378 = t249 * mrSges(7,1);
t229 = t376 + t378;
t467 = t210 * t229;
t171 = t467 / 0.2e1;
t277 = t376 / 0.2e1 + t378 / 0.2e1;
t269 = t210 * t277;
t484 = t171 + t269;
t496 = t484 * qJD(6);
t230 = Ifges(7,5) * t249 + Ifges(7,6) * t251;
t441 = t274 * t230;
t446 = Ifges(6,6) * t274;
t241 = Ifges(7,4) * t251;
t234 = Ifges(7,1) * t249 + t241;
t333 = t251 * t234;
t391 = Ifges(7,4) * t249;
t232 = Ifges(7,2) * t251 + t391;
t341 = t249 * t232;
t457 = t341 / 0.4e1 - t333 / 0.4e1;
t495 = t457 * t210 + t441 / 0.4e1 - t446 / 0.2e1;
t228 = -mrSges(7,1) * t251 + t249 * mrSges(7,2);
t296 = sin(pkin(10)) * pkin(1) + qJ(3);
t282 = pkin(7) + t296;
t271 = t247 * t282;
t272 = t246 * t282;
t195 = -t271 * t406 - t272 * t408;
t157 = -t226 * pkin(8) + t195;
t196 = t271 * t408 - t272 * t406;
t257 = t225 * pkin(8) + t196;
t431 = t157 * t250 + t257 * t407;
t460 = t431 * t228;
t468 = t431 * mrSges(6,1);
t437 = -Ifges(7,2) * t249 + t241;
t454 = Ifges(7,6) * t274 - t210 * t437;
t485 = t251 * t454;
t235 = Ifges(7,1) * t251 - t391;
t455 = Ifges(7,5) * t274 - t210 * t235;
t486 = t249 * t455;
t80 = -t407 * t157 + t250 * t257;
t490 = t80 * mrSges(6,2);
t494 = t441 / 0.2e1 - t446 + t460 - t468 + t490 + t485 / 0.2e1 + t486 / 0.2e1;
t493 = t460 / 0.2e1 - t468 / 0.2e1 + t485 / 0.4e1 + t486 / 0.4e1 + t490 / 0.2e1;
t139 = t229 * t274;
t283 = -cos(pkin(10)) * pkin(1) - pkin(3) * t247 - pkin(2);
t213 = -pkin(4) * t225 + t283;
t108 = pkin(5) * t210 - pkin(9) * t274 + t213;
t42 = t108 * t251 - t249 * t431;
t43 = t249 * t108 + t251 * t431;
t492 = t431 * t139 + t42 * t478 + t43 * t479 - t80 * t467;
t491 = -t462 / 0.2e1;
t489 = t249 * t80;
t488 = t250 * t80;
t487 = t251 * t80;
t398 = t431 * t80;
t367 = t80 * t274;
t483 = -m(7) * t431 + t467;
t444 = t274 * mrSges(6,1);
t469 = t210 * mrSges(6,2);
t306 = -t469 + t444;
t278 = t469 / 0.2e1 - t444 / 0.2e1;
t451 = pkin(5) * t274;
t244 = t249 ^ 2;
t245 = t251 ^ 2;
t331 = t244 + t245;
t465 = t210 * t331;
t480 = -pkin(9) * t465 - t451;
t379 = t245 * mrSges(7,3);
t380 = t244 * mrSges(7,3);
t476 = (-t379 / 0.2e1 - t380 / 0.2e1) * t210;
t404 = pkin(4) * t250;
t238 = pkin(9) + t404;
t329 = t407 * pkin(4);
t239 = -t329 - pkin(5);
t475 = -t238 * t465 + t239 * t274;
t426 = m(7) / 0.2e1;
t474 = -t210 / 0.2e1;
t473 = -t229 / 0.2e1;
t420 = -t467 / 0.2e1;
t442 = t274 * t228;
t421 = t442 / 0.2e1;
t470 = mrSges(7,3) * t331;
t208 = Ifges(6,4) * t210;
t466 = t210 * t250;
t240 = Ifges(7,5) * t251;
t386 = Ifges(7,6) * t249;
t438 = t240 - t386;
t464 = t210 * t438;
t360 = t274 * t210;
t440 = (-t386 / 0.2e1 + t240 / 0.2e1) * t210;
t369 = t431 * t210;
t151 = pkin(9) * t210 + t451;
t400 = t226 * pkin(4);
t109 = t151 + t400;
t51 = t249 * t109 - t487;
t372 = t51 * t251;
t50 = t109 * t251 + t489;
t373 = t50 * t249;
t289 = t372 - t373;
t458 = t232 / 0.4e1 - t235 / 0.4e1;
t439 = -t341 / 0.2e1 + t333 / 0.2e1;
t453 = Ifges(7,6) / 0.2e1;
t415 = t210 / 0.2e1;
t417 = t274 / 0.2e1;
t392 = Ifges(6,4) * t274;
t452 = -t392 / 0.2e1;
t448 = mrSges(6,3) * t274;
t385 = Ifges(7,3) * t274;
t354 = t274 * t249;
t145 = -mrSges(7,2) * t210 - mrSges(7,3) * t354;
t353 = t274 * t251;
t148 = t210 * mrSges(7,1) - mrSges(7,3) * t353;
t410 = -t251 / 0.2e1;
t412 = -t249 / 0.2e1;
t276 = t145 * t412 + t148 * t410;
t432 = -t464 / 0.4e1 + t80 * t473;
t430 = t480 * t426 + t278;
t429 = (-t431 * t426 + t171) * pkin(5) + (t289 * t426 + t499) * pkin(9) + (t372 / 0.2e1 - t373 / 0.2e1) * mrSges(7,3) + t493;
t428 = t226 ^ 2;
t427 = m(6) / 0.2e1;
t425 = m(6) * pkin(4);
t424 = m(7) * pkin(4);
t423 = -mrSges(7,1) / 0.2e1;
t422 = -mrSges(7,2) / 0.2e1;
t419 = t139 / 0.2e1;
t418 = Ifges(6,2) * t415 + t452;
t413 = t238 / 0.2e1;
t411 = t249 / 0.2e1;
t409 = t251 / 0.2e1;
t401 = pkin(5) * t229;
t394 = mrSges(7,3) * t274;
t390 = Ifges(6,5) * t210;
t381 = t226 * mrSges(5,1);
t377 = t249 * mrSges(7,3);
t375 = t251 * t43;
t56 = t151 * t251 + t489;
t371 = t56 * t249;
t57 = t249 * t151 - t487;
t370 = t57 * t251;
t307 = t245 / 0.2e1 + t244 / 0.2e1;
t297 = mrSges(7,3) * t307;
t21 = t442 * t415 + (t274 * t297 - t276) * t274;
t364 = qJD(1) * t21;
t259 = (t249 * t51 + t251 * t50) * t426 + t479 * t411 + t478 * t409;
t281 = (-t274 * t407 - t466) * pkin(4) * t427 + t475 * t426;
t322 = t470 * t474 + t421;
t266 = t281 + t322;
t218 = t225 * mrSges(5,2);
t298 = -t218 - t306;
t17 = (-mrSges(5,1) - t425 / 0.2e1) * t226 - t259 + t266 + t298;
t362 = t17 * qJD(1);
t258 = -m(7) * (t249 * t57 + t251 * t56) / 0.2e1 + t479 * t412 + t478 * t410 + t278;
t280 = t322 + t430;
t18 = t258 + t280;
t361 = t18 * qJD(1);
t352 = t239 * t442;
t351 = t239 * t229;
t336 = t251 * t145;
t312 = t336 / 0.2e1;
t344 = t249 * t148;
t275 = t312 - t344 / 0.2e1;
t24 = t269 - t275;
t350 = t24 * qJD(1);
t101 = Ifges(7,6) * t210 + t274 * t437;
t349 = t249 * t101;
t104 = Ifges(7,5) * t210 + t235 * t274;
t339 = t251 * t104;
t263 = m(7) * (-t274 * t465 + t360);
t35 = t263 / 0.2e1;
t332 = t35 * qJD(1);
t330 = mrSges(7,3) * t370;
t325 = -t377 / 0.2e1;
t321 = t249 * t407;
t309 = t232 / 0.2e1 - t235 / 0.2e1;
t308 = -t437 / 0.2e1 - t234 / 0.2e1;
t301 = t235 * t411 + t409 * t437 + t439;
t300 = t407 * t422;
t299 = -t321 / 0.2e1;
t150 = Ifges(6,1) * t274 - t208;
t98 = Ifges(7,3) * t210 + t274 * t438;
t2 = (Ifges(5,4) * t225 + (Ifges(5,1) - Ifges(5,2)) * t226) * t225 - t428 * Ifges(5,4) + m(7) * (t42 * t50 + t43 * t51 + t398) + t283 * (t218 + t381) + (mrSges(6,1) * t400 - Ifges(6,1) * t417) * t210 - t454 * t354 / 0.2e1 + t455 * t353 / 0.2e1 + t104 * t491 + t101 * t463 / 0.2e1 + (mrSges(6,2) * t400 + t418) * t274 + (-t392 + t98) * t417 + (t385 - t464) * t415 + t51 * t145 + t50 * t148 + (-Ifges(6,2) * t274 + t150 - t208) * t474 + (m(6) * t400 + t306) * t213 + t492;
t290 = t249 * t42 - t375;
t262 = t210 * t290 + t367;
t7 = t139 * t417 - t467 * t415 + t499 * t274 - t275 * t210 + (t274 * t289 + t262 + t369) * t426;
t293 = t2 * qJD(1) + t7 * qJD(2);
t5 = m(7) * (t42 * t56 + t43 * t57 + t398) + t57 * t145 + t56 * t148 + (t98 / 0.2e1 + t418 + t213 * mrSges(6,1) + t454 * t412 + t455 * t409 + t452) * t274 + (t208 / 0.2e1 - t213 * mrSges(6,2) - t150 / 0.2e1 + t349 / 0.2e1 - t339 / 0.2e1 - t440 + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t274) * t210 + t492;
t288 = t370 - t371;
t8 = ((t288 + t80) * t426 + t419 + t499) * t274 + ((t290 + t431) * t426 + t420 - t275) * t210;
t292 = t5 * qJD(1) + t8 * qJD(2);
t140 = t232 * t274;
t141 = t234 * t274;
t6 = -t80 * t442 + t42 * t145 - t43 * t148 + (t290 * mrSges(7,3) + t101 * t410 - t141 * t409 + t230 * t474 + (t104 - t140) * t412) * t274;
t287 = t6 * qJD(1) - t21 * qJD(2);
t286 = t7 * qJD(1) + qJD(2) * t263;
t37 = m(7) * (0.1e1 - t331) * t360;
t285 = t8 * qJD(1) + t37 * qJD(2);
t13 = (t139 + t448) * t274 + (t225 ^ 2 + t428) * mrSges(5,3) - (-t210 * mrSges(6,3) + t336 - t344) * t210 + m(7) * t262 + m(6) * (t367 - t369) + m(5) * (-t195 * t226 + t196 * t225) + (m(4) * t296 + mrSges(4,3)) * (t246 ^ 2 + t247 ^ 2);
t284 = -qJD(1) * t13 - qJD(2) * t35;
t270 = t331 * t407;
t256 = (t228 - mrSges(6,1)) * t404 + (-mrSges(6,2) + t470) * t329;
t123 = (t238 * t270 + t239 * t250) * t424 + t256;
t253 = ((t270 * t274 + t466) * pkin(4) + t475) * t426 + t421 + t278 + t476;
t23 = -t442 / 0.2e1 + t253 - t476 - t430;
t252 = (t239 * t431 + (-t321 * t42 + t375 * t407 + t488) * pkin(4)) * t426 + Ifges(6,5) * t474 + t239 * t420 + t404 * t419 + t56 * t325 + t330 / 0.2e1 + pkin(4) * t148 * t299 + t312 * t329 + (t288 * t426 + t499) * t238 + t493 + t495;
t4 = -t252 - t390 / 0.2e1 + t429 + t495;
t268 = -t4 * qJD(1) + t23 * qJD(2) + t123 * qJD(4);
t152 = t301 + t351;
t65 = t420 + t269;
t265 = t385 / 0.2e1 + t50 * mrSges(7,1) / 0.2e1 + t51 * t422;
t9 = t352 / 0.2e1 + (Ifges(7,5) * t474 + t140 / 0.4e1 - t104 / 0.4e1 + t148 * t413) * t251 + (t210 * t453 + t101 / 0.4e1 + t141 / 0.4e1 + t145 * t413) * t249 + (t458 * t251 + (t437 / 0.4e1 + t234 / 0.4e1) * t249 + t238 * t297) * t274 + t265 + t432;
t267 = -t9 * qJD(1) - t65 * qJD(2) + t152 * qJD(4);
t264 = -t385 / 0.2e1 + t56 * t423 + t57 * mrSges(7,2) / 0.2e1;
t106 = (pkin(5) / 0.2e1 - t239 / 0.2e1) * t229 + (pkin(4) * t300 + t308) * t251 + (t329 * t423 + t309) * t249;
t260 = -t249 * t141 / 0.4e1 - t251 * t140 / 0.4e1 - t349 / 0.4e1 + t339 / 0.4e1 - t432 + (t325 + t377 / 0.2e1) * t43 - t458 * t353 - (t437 + t234) * t354 / 0.4e1;
t255 = t260 + (-t307 * t394 + t276) * pkin(9) + pkin(5) * t421;
t11 = t255 + t264 + t440;
t153 = t249 * t309 + t251 * t308 + t401;
t64 = (t473 + t277) * t210;
t261 = t11 * qJD(1) - t64 * qJD(2) - t106 * qJD(4) - t153 * qJD(5);
t107 = -t401 / 0.2e1 + t351 / 0.2e1 + (mrSges(7,1) * t299 + t251 * t300) * pkin(4) + t301;
t25 = t269 + t275;
t22 = t253 + t280;
t20 = t400 * t427 + t259 + t266;
t19 = -t258 + t280;
t12 = Ifges(7,5) * t491 + t463 * t453 + t255 - t264;
t10 = t260 - t440 - t352 / 0.2e1 + t265 + (-t331 * t394 / 0.2e1 + t276) * t238;
t3 = t252 - (Ifges(6,5) / 0.2e1 - t457) * t210 + (-Ifges(6,6) / 0.2e1 + t230 / 0.4e1) * t274 + t429;
t1 = qJD(3) * t35 + qJD(4) * t7 + qJD(5) * t8 - qJD(6) * t21;
t14 = [qJD(3) * t13 + qJD(4) * t2 + qJD(5) * t5 + qJD(6) * t6, t1, qJD(4) * t20 + qJD(5) * t19 + qJD(6) * t25 - t284, t20 * qJD(3) + t3 * qJD(5) + t10 * qJD(6) + t293 + (-t404 * t448 + (-t407 * t431 - t488) * t425 - t195 * mrSges(5,2) - t196 * mrSges(5,1) - t390 + Ifges(5,5) * t225 - Ifges(5,6) * t226 - (-mrSges(6,3) * t329 + t439) * t210 - t483 * t239 + (m(7) * t289 + t500) * t238 + t289 * mrSges(7,3) + t494) * qJD(4), t19 * qJD(3) + t3 * qJD(4) + t12 * qJD(6) + t292 + (-mrSges(7,3) * t371 + t330 + (-Ifges(6,5) - t439) * t210 + t483 * pkin(5) + (m(7) * t288 + t500) * pkin(9) + t494) * qJD(5), t25 * qJD(3) + t10 * qJD(4) + t12 * qJD(5) + (-t43 * mrSges(7,1) - t42 * mrSges(7,2) - t441) * qJD(6) + t287; t1, qJD(4) * t263 + qJD(5) * t37, t332, t22 * qJD(5) + t496 + t286 + (0.2e1 * t281 + t298 + t442 - t381 - (t379 + t380) * t210) * qJD(4), t22 * qJD(4) + (m(7) * t480 - mrSges(7,3) * t465 - t306 + t442) * qJD(5) + t496 + t285, qJD(6) * t442 - t364 + (qJD(4) + qJD(5)) * t484; -qJD(4) * t17 - qJD(5) * t18 - qJD(6) * t24 + t284, -t332, 0, -t362, -t361, -qJD(6) * t229 - t350; qJD(3) * t17 - qJD(5) * t4 - qJD(6) * t9 - t293, qJD(5) * t23 - qJD(6) * t65 - t286, t362, qJD(5) * t123 + qJD(6) * t152 ((-pkin(5) * t250 + pkin(9) * t270) * t424 + t256) * qJD(5) + t107 * qJD(6) + t268, t107 * qJD(5) + (t228 * t238 + t438) * qJD(6) + t267; qJD(3) * t18 + qJD(4) * t4 + qJD(6) * t11 - t292, -qJD(4) * t23 - qJD(6) * t64 - t285, t361, -qJD(6) * t106 - t268, -t153 * qJD(6) (pkin(9) * t228 + t438) * qJD(6) + t261; qJD(3) * t24 + qJD(4) * t9 - qJD(5) * t11 - t287, qJD(4) * t65 + qJD(5) * t64 + t364, t350, qJD(5) * t106 - t267, -t261, 0;];
Cq  = t14;
