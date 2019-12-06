% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:15:02
% EndTime: 2019-12-05 17:15:23
% DurationCPUTime: 5.50s
% Computational Cost: add. (8973->370), mult. (21189->524), div. (0->0), fcn. (22659->10), ass. (0->241)
t250 = sin(qJ(5));
t252 = cos(qJ(5));
t251 = sin(qJ(3));
t253 = cos(qJ(3));
t249 = sin(pkin(5));
t415 = sin(qJ(2));
t330 = t249 * t415;
t380 = cos(pkin(5));
t210 = t380 * t251 + t253 * t330;
t279 = t251 * t330 - t380 * t253;
t414 = sin(qJ(4));
t416 = cos(qJ(4));
t270 = t416 * t210 - t414 * t279;
t254 = cos(qJ(2));
t360 = t249 * t254;
t104 = -t250 * t360 + t252 * t270;
t148 = t414 * t210 + t416 * t279;
t388 = t252 * mrSges(6,2);
t392 = t250 * mrSges(6,1);
t291 = t392 / 0.2e1 + t388 / 0.2e1;
t282 = t291 * t148;
t391 = t250 * mrSges(6,3);
t334 = t391 / 0.2e1;
t335 = -t391 / 0.2e1;
t227 = t388 + t392;
t374 = t148 * t227;
t474 = t282 + t374 / 0.2e1 + (t334 + t335) * t104;
t378 = t104 * t252;
t103 = -t250 * t270 - t252 * t360;
t379 = t103 * t250;
t482 = m(6) * (t270 - t378 + t379) * t148;
t486 = t482 * qJD(1);
t491 = t474 * qJD(5) + t486;
t490 = qJD(3) + qJD(4);
t327 = t416 * t251;
t220 = -t414 * t253 - t327;
t196 = t220 * t360;
t325 = t414 * t251;
t219 = -t416 * t253 + t325;
t197 = t219 * t360;
t226 = -mrSges(6,1) * t252 + mrSges(6,2) * t250;
t489 = -t197 * mrSges(5,2) / 0.2e1 - (-t226 / 0.2e1 + mrSges(5,1) / 0.2e1) * t196;
t165 = t227 * t219;
t346 = t220 * t391;
t170 = -mrSges(6,2) * t219 + t346;
t387 = t252 * mrSges(6,3);
t172 = mrSges(6,1) * t219 + t220 * t387;
t417 = t252 / 0.2e1;
t420 = -t250 / 0.2e1;
t488 = -t170 * t417 - t172 * t420 - t165 / 0.2e1;
t162 = -t252 * t197 + t250 * t330;
t371 = t162 * t252;
t161 = t250 * t197 + t252 * t330;
t372 = t161 * t250;
t487 = (t372 / 0.2e1 - t371 / 0.2e1) * mrSges(6,3) + t489;
t322 = -t360 / 0.2e1;
t240 = Ifges(6,4) * t252;
t308 = Ifges(6,2) * t250 - t240;
t111 = -Ifges(6,6) * t220 + t308 * t219;
t409 = Ifges(6,4) * t250;
t309 = Ifges(6,1) * t252 - t409;
t113 = -Ifges(6,5) * t220 - t309 * t219;
t307 = Ifges(6,5) * t250 + Ifges(6,6) * t252;
t286 = t220 * t307;
t419 = t250 / 0.2e1;
t230 = Ifges(6,1) * t250 + t240;
t354 = t252 * t230;
t229 = Ifges(6,2) * t252 + t409;
t357 = t250 * t229;
t454 = -t354 / 0.2e1 + t357 / 0.2e1;
t284 = t111 * t417 + t113 * t419 + Ifges(5,6) * t220 - t286 / 0.2e1 + (-Ifges(5,5) + t454) * t219;
t464 = pkin(8) + pkin(7);
t353 = t464 * t253;
t439 = -t325 * t464 + t416 * t353;
t456 = t439 * t226;
t461 = t439 * mrSges(5,1);
t189 = t327 * t464 + t414 * t353;
t478 = t189 * mrSges(5,2);
t485 = t284 + t456 - t461 + t478;
t483 = t456 / 0.2e1 - t461 / 0.2e1 + t478 / 0.2e1;
t238 = -pkin(3) * t253 - pkin(2);
t479 = m(5) * t238;
t477 = t189 * t250;
t476 = t189 * t252;
t447 = t189 * t270;
t475 = t414 * t189;
t370 = t439 * t189;
t239 = Ifges(6,5) * t252;
t407 = Ifges(6,6) * t250;
t473 = -t239 / 0.2e1 + t407 / 0.2e1;
t444 = t270 * t226;
t448 = t270 * mrSges(5,1);
t462 = t148 * mrSges(5,2);
t472 = t444 - t448 + t462;
t471 = -t444 / 0.2e1 + t448 / 0.2e1 - t462 / 0.2e1;
t245 = t250 ^ 2;
t247 = t252 ^ 2;
t352 = t245 + t247;
t459 = t148 * t352;
t468 = -pkin(4) * t270 - pkin(9) * t459;
t167 = pkin(4) * t219 + pkin(9) * t220 + t238;
t73 = t167 * t250 + t252 * t439;
t386 = t252 * t73;
t72 = t167 * t252 - t250 * t439;
t304 = -t250 * t72 + t386;
t413 = pkin(3) * t251;
t433 = m(6) / 0.2e1;
t166 = t227 * t220;
t169 = mrSges(6,2) * t220 + t219 * t391;
t171 = -mrSges(6,1) * t220 + t219 * t387;
t440 = -t270 * t166 / 0.2e1 + t103 * t171 / 0.2e1 + t104 * t169 / 0.2e1;
t175 = -pkin(4) * t220 + pkin(9) * t219;
t168 = t175 + t413;
t80 = t168 * t252 + t477;
t81 = t168 * t250 - t476;
t467 = m(5) * t413 * t322 + (t103 * t80 + t104 * t81 + t148 * t439 + t447) * t433 + t440 + (-t304 * t433 + t488) * t148;
t173 = -mrSges(5,1) * t220 - mrSges(5,2) * t219;
t114 = t219 * Ifges(6,5) - t309 * t220;
t356 = t252 * t114;
t112 = t219 * Ifges(6,6) + t308 * t220;
t359 = t250 * t112;
t418 = -t252 / 0.2e1;
t460 = Ifges(5,4) + t473;
t466 = (-t356 / 0.2e1 + t359 / 0.2e1 + t460 * t219) * t219 - t439 * t166 + (t113 * t418 + t111 * t419 - t460 * t220 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2)) * t219) * t220 - t189 * t165 + t73 * t169 + t72 * t171 + t238 * t173;
t347 = t414 * pkin(3);
t236 = t347 + pkin(9);
t348 = t416 * pkin(3);
t237 = -t348 - pkin(4);
t465 = -t236 * t459 + t237 * t270;
t463 = pkin(4) * t439;
t457 = t237 * t439;
t455 = t251 ^ 2 + t253 ^ 2;
t228 = t251 * mrSges(4,1) + t253 * mrSges(4,2);
t446 = t228 * t322;
t382 = t81 * t252;
t383 = t80 * t250;
t303 = t382 - t383;
t443 = -t253 * mrSges(4,1) + t251 * mrSges(4,2);
t441 = t308 * t418 + t309 * t419;
t317 = t239 - t407;
t438 = -Ifges(6,3) * t220 / 0.2e1 + (t317 / 0.4e1 + t473) * t219;
t355 = t252 * t169;
t350 = pkin(9) * t355;
t358 = t250 * t171;
t351 = pkin(9) * t358;
t412 = pkin(4) * t165;
t437 = (t303 * pkin(9) - t463) * t433 + t412 / 0.2e1 - t351 / 0.2e1 + t350 / 0.2e1 + (t382 / 0.2e1 - t383 / 0.2e1) * mrSges(6,3) + t483;
t333 = t387 / 0.2e1;
t435 = t161 * t335 + t162 * t333 - t489;
t434 = -m(6) / 0.2e1;
t432 = m(5) * pkin(3);
t431 = m(6) * pkin(3);
t430 = -mrSges(6,1) / 0.2e1;
t429 = mrSges(6,2) / 0.2e1;
t428 = Ifges(6,3) / 0.2e1;
t164 = t226 * t220;
t424 = -t164 / 0.2e1;
t421 = t172 / 0.2e1;
t411 = pkin(4) * t227;
t396 = t219 * mrSges(5,3);
t395 = t220 * mrSges(5,3);
t394 = t245 * mrSges(6,3);
t393 = t247 * mrSges(6,3);
t89 = t175 * t250 - t476;
t381 = t89 * t252;
t375 = t148 * t196;
t368 = t189 * t196;
t326 = t415 * t249 ^ 2;
t22 = m(6) * (t103 * t161 + t104 * t162 - t375) + m(4) * (-t326 + (t210 * t253 + t251 * t279) * t249) * t254 + (-t197 * t270 - t326 * t254 - t375) * m(5);
t364 = t22 * qJD(1);
t362 = t237 * t165;
t361 = t237 * t227;
t344 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t343 = t236 * t358;
t342 = t236 * t355;
t337 = -t394 / 0.2e1;
t336 = -t393 / 0.2e1;
t329 = t250 * t416;
t328 = t252 * t416;
t215 = -t357 / 0.2e1;
t313 = t215 + t354 / 0.2e1 + t441;
t312 = t173 * t322;
t311 = -t328 / 0.2e1;
t310 = mrSges(6,3) * (t247 / 0.2e1 + t245 / 0.2e1);
t296 = t371 - t372;
t259 = (-t196 * t237 + t296 * t236) * t433 + (t416 * t196 - t414 * t197) * t432 / 0.2e1 + t446;
t2 = -t259 + t312 + t446 + t467 + t487;
t174 = mrSges(5,1) * t219 - mrSges(5,2) * t220;
t3 = -pkin(2) * t228 + t81 * t170 + t80 * t172 + (-Ifges(4,4) * t251 + pkin(3) * t174) * t251 + t413 * t479 + m(6) * (t72 * t80 + t73 * t81 + t370) + (Ifges(4,4) * t253 + (Ifges(4,1) - Ifges(4,2)) * t251) * t253 + t466;
t306 = t2 * qJD(1) + t3 * qJD(2);
t88 = t175 * t252 + t477;
t258 = ((t439 - t304) * t433 + t488) * t148 + (t103 * t88 + t104 * t89 + t447) * t433 + t312 + t440;
t273 = m(6) * (pkin(4) * t196 + t296 * pkin(9));
t4 = -t273 / 0.2e1 + t258 + t487;
t7 = m(6) * (t72 * t88 + t73 * t89 + t370) + t89 * t170 + t88 * t172 + t466;
t305 = t4 * qJD(1) + t7 * qJD(2);
t302 = -t88 * t250 + t381;
t11 = t73 * t172 - t189 * t164 + (t114 * t420 + t112 * t418 - mrSges(6,3) * t386 - t219 * t307 / 0.2e1 + (t230 * t417 + t215) * t220) * t220 + (-t170 + t346) * t72;
t277 = -t103 * t170 / 0.2e1 + t104 * t421 + t148 * t424;
t292 = t161 * mrSges(6,1) / 0.2e1 - t162 * mrSges(6,2) / 0.2e1;
t16 = (t379 / 0.2e1 - t378 / 0.2e1) * t220 * mrSges(6,3) + t277 + t292;
t299 = -t16 * qJD(1) - t11 * qJD(2);
t295 = t81 * t429 + t80 * t430;
t294 = t89 * t429 + t88 * t430;
t288 = t170 * t420 + t172 * t418;
t283 = t352 * t416;
t257 = ((-t103 * t329 + t104 * t328 + t414 * t148) * pkin(3) + t465) * t433 + t148 * t337 + t148 * t336 - t471;
t263 = t468 * t434 - (t337 + t336) * t148 + t471;
t13 = t257 + t263;
t264 = (t226 - mrSges(5,1)) * t347 + (t352 * mrSges(6,3) - mrSges(5,2)) * t348;
t52 = (t283 * t236 + t414 * t237) * t431 + t264;
t255 = (t302 * t236 + t457) * t434 + t362 / 0.2e1 + t343 / 0.2e1 - t342 / 0.2e1 + t88 * t334 - mrSges(6,3) * t381 / 0.2e1 + t166 * t347 / 0.2e1 + ((t73 * t328 - t72 * t329 + t475) * t434 + t329 * t421 + t170 * t311) * pkin(3) - t483;
t8 = t255 + t437;
t281 = t13 * qJD(1) - t8 * qJD(2) + t52 * qJD(3);
t265 = (t230 / 0.4e1 + t240 / 0.4e1 + t344 * t250) * t250 + (0.3e1 / 0.4e1 * t409 + t229 / 0.4e1 - t344 * t252) * t252;
t261 = t236 * t310 + t265;
t276 = t189 * t227 / 0.2e1 - t359 / 0.4e1 + t356 / 0.4e1;
t266 = t288 * t236 + t237 * t164 / 0.2e1 + t276;
t278 = (-0.3e1 / 0.4e1 * t407 + 0.3e1 / 0.4e1 * t239) * t219;
t10 = t278 + (t428 + t261) * t220 + t266 + t295;
t24 = -t374 / 0.2e1 + t282;
t65 = t313 + t361;
t280 = -t24 * qJD(1) + t10 * qJD(2) + t65 * qJD(3);
t262 = pkin(9) * t310 + t265;
t267 = pkin(4) * t424 + t288 * pkin(9) + t276;
t15 = t278 + (t428 + t262) * t220 + t267 + t294;
t26 = (-t227 / 0.2e1 + t291) * t148;
t268 = -t441 + t454;
t272 = (mrSges(6,2) * t311 + t329 * t430) * pkin(3);
t42 = (pkin(4) / 0.2e1 - t237 / 0.2e1) * t227 + t272 + t268;
t82 = t268 + t411;
t274 = t26 * qJD(1) - t15 * qJD(2) + t42 * qJD(3) + t82 * qJD(4);
t43 = -t411 / 0.2e1 + t361 / 0.2e1 + t272 + t313;
t17 = -t277 + t292 + (t103 * t335 + t104 * t333) * t220;
t14 = t262 * t220 + t267 - t294 + t438;
t12 = t257 - t263;
t9 = t261 * t220 + t266 - t295 + t438;
t6 = -t255 + t284 + t437;
t5 = t273 / 0.2e1 + t258 + t435;
t1 = (-t228 / 0.2e1 - t173 / 0.2e1) * t360 + t259 + t435 + t467;
t18 = [t22 * qJD(2) + t490 * t482, t1 * qJD(3) + t5 * qJD(4) + t17 * qJD(5) + t364 + (t161 * t172 + m(6) * (t161 * t72 + t162 * t73 - t368) + m(5) * (-t197 * t439 - t368) + t162 * t170 + t197 * t396 + m(4) * (t455 * t254 * pkin(7) - t415 * pkin(2)) * t249 + (-mrSges(3,1) + t174 + t443 + t479) * t330 + (t166 + t395) * t196 + (t455 * mrSges(4,3) - mrSges(3,2)) * t360) * qJD(2), t1 * qJD(2) + (m(6) * t465 - t416 * t270 * t432 + t279 * mrSges(4,2) - t210 * mrSges(4,1) - (t414 * t432 + t393 + t394) * t148 + t472) * qJD(3) + t12 * qJD(4) + t491, t5 * qJD(2) + t12 * qJD(3) + (m(6) * t468 - mrSges(6,3) * t459 + t472) * qJD(4) + t491, t17 * qJD(2) + (-mrSges(6,1) * t104 - mrSges(6,2) * t103) * qJD(5) + t490 * t474; qJD(3) * t2 + qJD(4) * t4 - qJD(5) * t16 - t364, qJD(3) * t3 + qJD(4) * t7 - qJD(5) * t11, (m(6) * (t236 * t303 + t457) + t342 - t343 + (-t416 * t439 - t475) * t432 + Ifges(4,5) * t253 - Ifges(4,6) * t251 - t362 + t348 * t396 + t347 * t395 + t443 * pkin(7) + t303 * mrSges(6,3) + t485) * qJD(3) + t6 * qJD(4) + t9 * qJD(5) + t306, t6 * qJD(3) + (m(6) * (pkin(9) * t302 - t463) + t350 - t351 + t412 + t302 * mrSges(6,3) + t485) * qJD(4) + t14 * qJD(5) + t305, t9 * qJD(3) + t14 * qJD(4) + (-mrSges(6,1) * t73 - mrSges(6,2) * t72 + t286) * qJD(5) + t299; -qJD(2) * t2 + qJD(4) * t13 - qJD(5) * t24 - t486, -qJD(4) * t8 + qJD(5) * t10 - t306, qJD(4) * t52 + qJD(5) * t65, ((-t414 * pkin(4) + t283 * pkin(9)) * t431 + t264) * qJD(4) + t43 * qJD(5) + t281, t43 * qJD(4) + (t226 * t236 + t317) * qJD(5) + t280; -qJD(2) * t4 - qJD(3) * t13 - qJD(5) * t26 - t486, qJD(3) * t8 + qJD(5) * t15 - t305, -qJD(5) * t42 - t281, -t82 * qJD(5), (pkin(9) * t226 + t317) * qJD(5) - t274; t16 * qJD(2) + t24 * qJD(3) + t26 * qJD(4), -qJD(3) * t10 - qJD(4) * t15 - t299, qJD(4) * t42 - t280, t274, 0;];
Cq = t18;
