% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR14_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:26
% EndTime: 2024-09-27 18:42:34
% DurationCPUTime: 5.34s
% Computational Cost: add. (28296->512), mult. (68962->654), div. (0->0), fcn. (70726->10), ass. (0->294)
t321 = cos(pkin(5));
t482 = t321 / 0.2e1;
t324 = sin(qJ(3));
t325 = sin(qJ(2));
t329 = cos(qJ(2));
t328 = cos(qJ(3));
t402 = t321 * t328;
t289 = (-t324 * t329 - t325 * t402) * pkin(1);
t403 = t321 * t324;
t290 = (-t325 * t403 + t328 * t329) * pkin(1);
t323 = sin(qJ(4));
t327 = cos(qJ(4));
t226 = t289 * t327 - t290 * t323;
t227 = t289 * t323 + t290 * t327;
t322 = sin(qJ(5));
t326 = cos(qJ(5));
t144 = t226 * t326 - t227 * t322;
t145 = t226 * t322 + t227 * t326;
t478 = pkin(3) * t327;
t315 = pkin(4) + t478;
t401 = t322 * t323;
t292 = -pkin(3) * t401 + t315 * t326;
t399 = t323 * t326;
t293 = pkin(3) * t399 + t315 * t322;
t499 = -pkin(3) / 0.2e1;
t398 = t144 * mrSges(6,1) / 0.2e1 - t145 * mrSges(6,2) / 0.2e1;
t505 = t398 - t227 * mrSges(5,2) / 0.2e1 + t226 * mrSges(5,1) / 0.2e1;
t530 = -m(6) / 0.2e1;
t531 = t505 - (t144 * t292 + t145 * t293) * t530 + t289 * mrSges(4,1) / 0.2e1 - t290 * mrSges(4,2) / 0.2e1 - m(5) * (t226 * t327 + t227 * t323) * t499;
t316 = pkin(1) * t329 + pkin(2);
t320 = sin(pkin(5));
t477 = pkin(3) * t328;
t291 = (-t316 - t477) * t320;
t285 = (-t323 * t324 + t327 * t328) * t320;
t476 = pkin(4) * t285;
t236 = t291 - t476;
t406 = t320 * t328;
t407 = t320 * t324;
t286 = -t323 * t406 - t327 * t407;
t225 = t285 * t322 - t286 * t326;
t369 = t326 * t285 + t286 * t322;
t517 = mrSges(6,1) * t225 + mrSges(6,2) * t369;
t107 = t236 * t517;
t304 = (-pkin(2) - t477) * t320;
t245 = t304 - t476;
t416 = t245 * t517;
t512 = Ifges(6,5) * t369;
t522 = Ifges(6,6) * t225;
t529 = t512 - t522;
t368 = Ifges(5,5) * t285 + Ifges(5,6) * t286 + t529;
t376 = -t522 / 0.2e1 + t512 / 0.2e1;
t527 = (-Ifges(6,6) * t321 / 0.2e1 - Ifges(6,4) * t225) * t225 + (Ifges(6,5) * t482 + Ifges(6,4) * t369 + (Ifges(6,1) - Ifges(6,2)) * t225) * t369;
t525 = t285 / 0.2e1;
t300 = (t326 * t327 - t401) * pkin(3);
t524 = t300 / 0.2e1;
t442 = t225 * mrSges(6,3);
t312 = pkin(2) * t402;
t491 = pkin(8) + pkin(9);
t271 = -t491 * t407 + t312;
t318 = t321 * pkin(3);
t254 = t318 + t271;
t392 = pkin(2) * t403;
t272 = t491 * t406 + t392;
t260 = t323 * t272;
t192 = t327 * t254 - t260;
t281 = t286 * pkin(10);
t372 = t321 * pkin(4) + t281;
t148 = t192 + t372;
t262 = t327 * t272;
t193 = t254 * t323 + t262;
t473 = pkin(10) * t285;
t159 = t193 + t473;
t424 = t159 * t326;
t84 = t148 * t322 + t424;
t432 = t84 * t225;
t305 = t316 * t402;
t479 = pkin(1) * t325;
t307 = pkin(8) * t320 + t479;
t373 = pkin(9) * t320 + t307;
t243 = -t324 * t373 + t305;
t237 = t318 + t243;
t382 = t316 * t403;
t244 = t328 * t373 + t382;
t240 = t327 * t244;
t156 = t237 * t323 + t240;
t437 = t286 * mrSges(5,3);
t131 = t156 * t437;
t228 = -mrSges(5,1) * t286 + mrSges(5,2) * t285;
t189 = t291 * t228;
t351 = (0.2e1 * Ifges(5,4) * t285 + Ifges(5,5) * t321) * t525 + (-Ifges(5,4) * t286 + Ifges(5,6) * t482 + (Ifges(5,1) - Ifges(5,2)) * (-t525 - t285 / 0.2e1)) * t286 + t368 * t482 + t527;
t238 = t323 * t244;
t155 = t327 * t237 - t238;
t426 = t155 * t285;
t120 = t155 + t372;
t133 = t156 + t473;
t428 = t133 * t322;
t66 = t120 * t326 - t428;
t435 = t66 * t369;
t427 = t133 * t326;
t67 = t120 * t322 + t427;
t51 = t67 * t442;
t521 = mrSges(5,3) * t426 + mrSges(6,3) * t435 - t107 - t131 - t189 - t351 + t51;
t299 = (-t322 * t327 - t399) * pkin(3);
t158 = t192 + t281;
t425 = t159 * t322;
t88 = t158 * t326 - t425;
t463 = t88 * mrSges(6,2);
t87 = -t158 * t322 - t424;
t464 = t87 * mrSges(6,1);
t458 = t464 / 0.2e1 - t463 / 0.2e1;
t83 = t148 * t326 - t425;
t520 = (t292 * t87 + t293 * t88 + t299 * t83 + t300 * t84) * t530 - t458;
t515 = pkin(4) / 0.2e1;
t441 = t369 * mrSges(6,3);
t138 = -mrSges(6,1) * t369 + mrSges(6,2) * t225;
t199 = -mrSges(6,2) * t321 + t441;
t200 = mrSges(6,1) * t321 - t442;
t229 = -mrSges(5,1) * t285 - mrSges(5,2) * t286;
t438 = t285 * mrSges(5,3);
t256 = -mrSges(5,2) * t321 + t438;
t257 = mrSges(5,1) * t321 + t437;
t302 = mrSges(4,1) * t321 - mrSges(4,3) * t407;
t303 = -mrSges(4,2) * t321 + mrSges(4,3) * t406;
t510 = t144 * t200 + t145 * t199 + t226 * t257 + t227 * t256 + t289 * t302 + t290 * t303 + (-t329 * mrSges(3,2) + (-mrSges(3,1) + (t138 + t229 + (-mrSges(4,1) * t328 + mrSges(4,2) * t324) * t320) * t320) * t325) * pkin(1);
t391 = pkin(3) * t407;
t475 = pkin(4) * t286;
t251 = t391 - t475;
t310 = Ifges(4,5) * t406;
t509 = t251 * t138 + t310 * t482;
t465 = t84 * mrSges(6,1);
t466 = t83 * mrSges(6,2);
t508 = -t465 / 0.2e1 - t466 / 0.2e1;
t471 = t67 * mrSges(6,1);
t472 = t66 * mrSges(6,2);
t507 = -t471 / 0.2e1 - t472 / 0.2e1;
t162 = -t243 * t323 - t240;
t142 = t162 - t473;
t163 = t327 * t243 - t238;
t143 = t281 + t163;
t78 = t142 * t322 + t143 * t326;
t467 = t78 * mrSges(6,2);
t77 = t142 * t326 - t143 * t322;
t468 = t77 * mrSges(6,1);
t459 = t468 / 0.2e1 - t467 / 0.2e1;
t294 = t299 * mrSges(6,1);
t436 = t300 * mrSges(6,2);
t503 = (mrSges(5,1) * t323 + mrSges(5,2) * t327) * pkin(3) - t294 + t436;
t180 = t293 * t442;
t268 = t323 * pkin(3) * t437;
t502 = -t180 / 0.2e1 + t268 / 0.2e1 + t199 * t524 + t299 * t200 / 0.2e1;
t501 = m(5) / 0.2e1;
t500 = m(6) / 0.2e1;
t498 = m(6) * pkin(4);
t497 = -t66 / 0.2e1;
t496 = -t67 / 0.2e1;
t495 = t78 / 0.2e1;
t494 = -t83 / 0.2e1;
t493 = -t84 / 0.2e1;
t197 = -t271 * t323 - t262;
t164 = t197 - t473;
t198 = t327 * t271 - t260;
t165 = t281 + t198;
t92 = t164 * t322 + t165 * t326;
t492 = t92 / 0.2e1;
t490 = t163 / 0.2e1;
t489 = t197 / 0.2e1;
t488 = t199 / 0.2e1;
t484 = -t292 / 0.2e1;
t483 = -t293 / 0.2e1;
t481 = -t322 / 0.2e1;
t480 = -t326 / 0.2e1;
t474 = pkin(4) * t322;
t132 = t155 + t281;
t72 = -t132 * t322 - t427;
t470 = t72 * mrSges(6,1);
t73 = t132 * t326 - t428;
t469 = t73 * mrSges(6,2);
t91 = t164 * t326 - t165 * t322;
t462 = t91 * mrSges(6,1);
t461 = t92 * mrSges(6,2);
t460 = t470 / 0.2e1 - t469 / 0.2e1;
t457 = t462 / 0.2e1 - t461 / 0.2e1;
t454 = pkin(3) * qJD(3);
t453 = pkin(4) * qJD(4);
t450 = t155 * mrSges(5,2);
t449 = t156 * mrSges(5,1);
t448 = t162 * mrSges(5,1);
t447 = t163 * mrSges(5,2);
t446 = t192 * mrSges(5,2);
t445 = t193 * mrSges(5,1);
t444 = t197 * mrSges(5,1);
t443 = t198 * mrSges(5,2);
t175 = t236 * t251;
t258 = -t307 * t324 + t305;
t259 = t307 * t328 + t382;
t265 = t291 * t391;
t296 = (mrSges(4,1) * t324 + mrSges(4,2) * t328) * t320;
t353 = Ifges(4,4) * t406 + Ifges(4,5) * t482 + (Ifges(4,1) - Ifges(4,2)) * t407;
t354 = -Ifges(4,4) * t407 - Ifges(4,6) * t321 + pkin(3) * t229;
t7 = m(6) * (t66 * t77 + t67 * t78 + t175) + m(5) * (t155 * t162 + t156 * t163 + t265) - t259 * t302 + t258 * t303 + t163 * t256 + t162 * t257 + (-t316 * t296 + (-t259 * mrSges(4,3) + t354) * t324 + (-t258 * mrSges(4,3) + t353) * t328) * t320 + t78 * t199 + t77 * t200 + t509 - t521;
t434 = t7 * qJD(1);
t113 = t138 * t475;
t118 = t155 * t256;
t119 = t156 * t257;
t45 = t72 * t200;
t46 = t73 * t199;
t8 = -t46 + t119 - m(6) * (-t236 * t475 + t66 * t72 + t67 * t73) - t118 + t113 - t45 + t521;
t433 = t8 * qJD(1);
t431 = t87 * t200;
t430 = t88 * t199;
t334 = t482 * t529 + t527;
t13 = t107 + t66 * t199 - t67 * t200 + (-t67 * t225 - t435) * mrSges(6,3) + t334;
t429 = t13 * qJD(1);
t423 = t192 * t256;
t422 = t193 * t286;
t393 = t320 * t479;
t394 = t320 ^ 2 * t479;
t21 = m(6) * (t144 * t66 + t145 * t67 + t236 * t393) + m(5) * (t155 * t226 + t156 * t227 + t291 * t393) + m(4) * (t258 * t289 + t259 * t290 - t316 * t394) + t510;
t421 = t21 * qJD(1);
t408 = t304 * t228;
t400 = t323 * t257;
t205 = t442 * t474;
t381 = -t441 / 0.2e1;
t397 = -t205 / 0.2e1 + t326 * pkin(4) * t381;
t115 = -t293 * mrSges(6,1) - t292 * mrSges(6,2);
t395 = qJD(5) * t115;
t390 = t326 * t441;
t389 = t327 * t438;
t384 = t66 / 0.2e1 + t83 / 0.2e1;
t383 = t496 + t493;
t380 = t441 / 0.2e1;
t297 = -pkin(8) * t407 + t312;
t375 = t297 / 0.2e1 + t258 / 0.2e1;
t298 = pkin(8) * t406 + t392;
t374 = -t298 / 0.2e1 - t259 / 0.2e1;
t370 = (-t236 - t245) * t286;
t367 = pkin(4) * t480 + t484;
t185 = t245 * t251;
t284 = t304 * t391;
t157 = t192 * t438;
t63 = t83 * t441;
t333 = t107 / 0.2e1 + t351 + t189 / 0.2e1 - t63 / 0.2e1 - t157 / 0.2e1 + (-t432 / 0.2e1 - t435 / 0.2e1) * mrSges(6,3) + t131 / 0.2e1 + (-t426 / 0.2e1 + t422 / 0.2e1) * mrSges(5,3) - t51 / 0.2e1 + t408 / 0.2e1 + t416 / 0.2e1;
t330 = (t489 + t162 / 0.2e1) * t257 + ((-pkin(2) / 0.2e1 - t316 / 0.2e1) * t296 + (mrSges(4,3) * t374 + t354) * t324 + (-mrSges(4,3) * t375 + t353) * t328) * t320 + (t91 / 0.2e1 + t77 / 0.2e1) * t200 + t333 + (t66 * t91 + t67 * t92 + t77 * t83 + t78 * t84 + t175 + t185) * t500 + (t155 * t197 + t156 * t198 + t162 * t192 + t163 * t193 + t265 + t284) * t501 + t374 * t302 + (t198 / 0.2e1 + t490) * t256 + (t492 + t495) * t199 + t375 * t303 + t509;
t2 = t330 - t531;
t340 = -mrSges(6,3) * t432 - t157 + t351 + t408 + t416 - t63;
t9 = mrSges(5,3) * t422 + m(6) * (t83 * t91 + t84 * t92 + t185) + m(5) * (t192 * t197 + t193 * t198 + t284) + t340 - t298 * t302 + t297 * t303 + t198 * t256 + t197 * t257 + t92 * t199 + t91 * t200 + (-pkin(2) * t296 + (-t298 * mrSges(4,3) + t354) * t324 + (-t297 * mrSges(4,3) + t353) * t328) * t320 + t509;
t365 = t2 * qJD(1) + t9 * qJD(2);
t364 = -t205 + t368;
t10 = -t113 + (-t257 + t437) * t193 + t340 + t423 + t430 + t431 + m(6) * (-t245 * t475 + t83 * t87 + t84 * t88);
t331 = t46 / 0.2e1 - t119 / 0.2e1 + t333 + t118 / 0.2e1 - t113 + t45 / 0.2e1 + t423 / 0.2e1 - t193 * t257 / 0.2e1 + t430 / 0.2e1 + t431 / 0.2e1;
t352 = t66 * t87 + t67 * t88 + t72 * t83 + t73 * t84;
t358 = m(6) * (t144 * t326 + t145 * t322);
t4 = t331 + 0.2e1 * (-t358 / 0.4e1 + m(6) * t370 / 0.4e1) * pkin(4) + t352 * t500 - t505;
t363 = t4 * qJD(1) + t10 * qJD(2);
t16 = t416 + t83 * t199 - t84 * t200 + (-t369 * t83 - t432) * mrSges(6,3) + t334;
t332 = (t236 / 0.2e1 + t245 / 0.2e1) * t517 + t384 * t199 + t383 * t200 + (t225 * t383 - t369 * t384) * mrSges(6,3) + t334;
t5 = t332 - t398;
t362 = t5 * qJD(1) + t16 * qJD(2);
t361 = (-t474 / 0.2e1 + t483) * mrSges(6,1);
t359 = m(6) * (t322 * t78 + t326 * t77);
t114 = -m(6) * (t292 * t299 + t293 * t300) + t503;
t355 = t292 * t381 + t256 * t478 / 0.2e1 + t502 + (t400 + t389) * t499;
t336 = (t292 * t72 + t293 * t73 + t299 * t66 + t300 * t67) * t500 + t355 + t460;
t12 = (t326 * t380 - t359 / 0.2e1) * pkin(4) + t205 / 0.2e1 + (t490 - t155 / 0.2e1) * mrSges(5,2) + (-t162 / 0.2e1 - t156 / 0.2e1) * mrSges(5,1) + t336 - t459;
t347 = (t322 * t92 + t326 * t91) * t498 / 0.2e1 + t397 + t457;
t15 = t292 * t380 + (t192 / 0.2e1 - t198 / 0.2e1) * mrSges(5,2) + (t193 / 0.2e1 + t489) * mrSges(5,1) + (t400 / 0.2e1 + (-t256 / 0.2e1 + t438 / 0.2e1) * t327) * pkin(3) + t347 - t502 + t520;
t357 = t12 * qJD(1) - t15 * qJD(2) - t114 * qJD(3);
t343 = (t225 * t483 + t369 * t484) * mrSges(6,3) + t292 * t488 + t200 * t483 + t376;
t339 = t343 - t376;
t18 = (t497 + t495) * mrSges(6,2) + (t496 - t77 / 0.2e1) * mrSges(6,1) + t339;
t20 = (t494 + t492) * mrSges(6,2) + (t493 - t91 / 0.2e1) * mrSges(6,1) + t339;
t356 = t18 * qJD(1) + t20 * qJD(2) + t115 * qJD(3);
t108 = -t294 / 0.2e1 + t361 + (t524 + t367) * mrSges(6,2);
t342 = (t326 * t488 + t200 * t481 + (t225 * t481 + t369 * t480) * mrSges(6,3)) * pkin(4) + t376;
t335 = t342 - t376;
t23 = (t497 + t73 / 0.2e1) * mrSges(6,2) + (t496 - t72 / 0.2e1) * mrSges(6,1) + t335;
t28 = (t494 + t88 / 0.2e1) * mrSges(6,2) + (t493 - t87 / 0.2e1) * mrSges(6,1) + t335;
t306 = (t322 * mrSges(6,1) + t326 * mrSges(6,2)) * pkin(4);
t350 = -qJD(1) * t23 - qJD(2) * t28 - qJD(3) * t108 + qJD(4) * t306;
t345 = -Ifges(4,6) * t407 - t292 * t441 - t180 + t268 + t310 + t368;
t341 = t343 + t376;
t338 = t342 + t376;
t301 = t306 * qJD(5);
t109 = -t436 / 0.2e1 + t294 / 0.2e1 + t367 * mrSges(6,2) + t361;
t27 = t338 + t458 + t508;
t22 = t338 + t460 + t507;
t19 = t341 + t457 + t508;
t17 = t341 + t459 + t507;
t14 = t347 - t443 / 0.2e1 + t444 / 0.2e1 - t446 / 0.2e1 - t445 / 0.2e1 + t355 + t368 - t520;
t11 = t336 + t359 * t515 - t447 / 0.2e1 + t448 / 0.2e1 - t450 / 0.2e1 - t449 / 0.2e1 + t368 + t397 + t459;
t6 = t332 + t398;
t3 = (pkin(4) * t370 + t352) * t500 + t331 + t358 * t515 + t505;
t1 = t330 + t531;
t24 = [qJD(2) * t21 + qJD(3) * t7 - qJD(4) * t8 + qJD(5) * t13, t1 * qJD(3) + t3 * qJD(4) + t6 * qJD(5) + t421 + (0.2e1 * (t144 * t83 + t145 * t84 + t245 * t393) * t500 + 0.2e1 * (t192 * t226 + t193 * t227 + t304 * t393) * t501 + m(4) * (-pkin(2) * t394 + t289 * t297 + t290 * t298) + t510) * qJD(2), t434 + t1 * qJD(2) + (t345 - t258 * mrSges(4,2) - t259 * mrSges(4,1) - t447 + t448 - t467 + t468 + m(6) * (t292 * t77 + t293 * t78)) * qJD(3) + t11 * qJD(4) + t17 * qJD(5) + (-t389 + m(5) * (t162 * t327 + t163 * t323)) * t454, -t433 + t3 * qJD(2) + t11 * qJD(3) + (t364 - t449 - t450 - t469 + t470) * qJD(4) + t22 * qJD(5) + (-t390 + m(6) * (t322 * t73 + t326 * t72)) * t453, t429 + t6 * qJD(2) + t17 * qJD(3) + t22 * qJD(4) + (t529 - t471 - t472) * qJD(5); qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t5 - t421, qJD(3) * t9 + qJD(4) * t10 + qJD(5) * t16, (m(6) * (t292 * t91 + t293 * t92) + t345 - t297 * mrSges(4,2) - t298 * mrSges(4,1) - t443 + t444 + t462 - t461) * qJD(3) + t14 * qJD(4) + t19 * qJD(5) + (-t389 + m(5) * (t197 * t327 + t198 * t323)) * t454 + t365, t14 * qJD(3) + (t364 - t445 - t446 - t463 + t464) * qJD(4) + t27 * qJD(5) + (-t390 + m(6) * (t322 * t88 + t326 * t87)) * t453 + t363, t19 * qJD(3) + t27 * qJD(4) + (t529 - t465 - t466) * qJD(5) + t362; -qJD(2) * t2 + qJD(4) * t12 + qJD(5) * t18 - t434, -qJD(4) * t15 + qJD(5) * t20 - t365, -qJD(4) * t114 + t395, ((t299 * t326 + t300 * t322) * t498 - t503) * qJD(4) + t109 * qJD(5) + t357, t109 * qJD(4) + t356 + t395; -qJD(2) * t4 - qJD(3) * t12 + qJD(5) * t23 + t433, qJD(3) * t15 + qJD(5) * t28 - t363, qJD(5) * t108 - t357, -t301, -t301 - t350; -qJD(2) * t5 - qJD(3) * t18 - qJD(4) * t23 - t429, -qJD(3) * t20 - qJD(4) * t28 - t362, -qJD(4) * t108 - t356, t350, 0;];
Cq = t24;
