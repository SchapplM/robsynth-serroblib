% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR9
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:14
% EndTime: 2019-12-05 17:19:32
% DurationCPUTime: 7.22s
% Computational Cost: add. (10108->508), mult. (24578->725), div. (0->0), fcn. (25365->10), ass. (0->276)
t316 = sin(qJ(3));
t472 = t316 / 0.2e1;
t313 = sin(pkin(5));
t315 = sin(qJ(4));
t317 = sin(qJ(2));
t319 = cos(qJ(4));
t320 = cos(qJ(3));
t321 = cos(qJ(2));
t386 = t320 * t321;
t219 = (-t315 * t386 + t317 * t319) * t313;
t220 = (t315 * t317 + t319 * t386) * t313;
t314 = sin(qJ(5));
t318 = cos(qJ(5));
t116 = t219 * t318 - t220 * t314;
t117 = t219 * t314 + t220 * t318;
t496 = mrSges(5,2) / 0.2e1;
t497 = -mrSges(5,1) / 0.2e1;
t498 = m(6) * pkin(4);
t512 = t116 * mrSges(6,1) / 0.2e1 - t117 * mrSges(6,2) / 0.2e1;
t530 = t219 * t497 + t220 * t496 - (t116 * t318 + t117 * t314) * t498 / 0.2e1 - t512;
t346 = t314 * t319 + t318 * t315;
t243 = t346 * t320;
t389 = t318 * t319;
t345 = t314 * t315 - t389;
t244 = t345 * t320;
t343 = Ifges(6,5) * t244 / 0.2e1 + Ifges(6,6) * t243 / 0.2e1;
t524 = -mrSges(6,1) / 0.2e1;
t529 = mrSges(6,2) / 0.2e1;
t293 = pkin(3) * t316 - pkin(8) * t320;
t394 = t315 * t316;
t225 = pkin(7) * t394 + t319 * t293;
t387 = t319 * t320;
t173 = pkin(4) * t316 - pkin(9) * t387 + t225;
t392 = t316 * t319;
t226 = -pkin(7) * t392 + t315 * t293;
t393 = t315 * t320;
t191 = -pkin(9) * t393 + t226;
t88 = t173 * t318 - t191 * t314;
t91 = t173 * t314 + t191 * t318;
t507 = t88 * t524 + t91 * t529 + t343;
t502 = Ifges(6,3) * t472 - t507;
t282 = -pkin(3) * t320 - pkin(8) * t316 - pkin(2);
t266 = t319 * t282;
t356 = -pkin(9) * t392 + t266;
t172 = (-pkin(7) * t315 - pkin(4)) * t320 + t356;
t224 = pkin(7) * t387 + t282 * t315;
t187 = -pkin(9) * t394 + t224;
t411 = t187 * t314;
t81 = t172 * t318 - t411;
t380 = pkin(7) * t393;
t186 = t356 - t380;
t98 = t186 * t318 - t411;
t528 = -t81 + t98;
t399 = t313 * t317;
t416 = cos(pkin(5));
t254 = t416 * t316 + t320 * t399;
t398 = t313 * t321;
t179 = -t254 * t315 - t319 * t398;
t180 = t254 * t319 - t315 * t398;
t359 = t318 * t179 - t180 * t314;
t93 = t179 * t314 + t180 * t318;
t27 = -t93 * mrSges(6,1) - t359 * mrSges(6,2);
t527 = t27 * qJD(5);
t493 = -pkin(9) - pkin(8);
t291 = t493 * t315;
t292 = t493 * t319;
t199 = t291 * t314 - t292 * t318;
t358 = t318 * t291 + t292 * t314;
t259 = Ifges(6,6) * t346;
t260 = Ifges(6,5) * t345;
t383 = -t260 - t259;
t36 = -t199 * mrSges(6,1) - t358 * mrSges(6,2) + t383;
t526 = t36 * qJD(5);
t431 = t199 * mrSges(6,3);
t447 = Ifges(6,4) * t346;
t525 = t431 - t447;
t410 = t187 * t318;
t82 = t172 * t314 + t410;
t97 = -t186 * t314 - t410;
t522 = t82 + t97;
t466 = pkin(4) * t319;
t303 = -pkin(3) - t466;
t370 = t316 * t398;
t408 = t220 * t319;
t409 = t219 * t315;
t518 = -(t408 / 0.2e1 - t409 / 0.2e1) * mrSges(5,3) - m(5) * (-pkin(3) * t370 + (t408 - t409) * pkin(8)) / 0.2e1 - m(6) * (t116 * t358 + t199 * t117 + t303 * t370) / 0.2e1;
t500 = m(6) / 0.2e1;
t517 = 0.2e1 * t500;
t501 = m(5) / 0.2e1;
t516 = 0.2e1 * t501;
t425 = t346 * mrSges(6,3);
t306 = Ifges(5,5) * t319;
t369 = -t306 / 0.2e1;
t444 = Ifges(5,6) * t315;
t513 = Ifges(4,4) + t369 + t444 / 0.2e1;
t242 = t346 * t316;
t230 = Ifges(6,4) * t242;
t241 = t314 * t394 - t316 * t389;
t136 = -Ifges(6,1) * t241 - t320 * Ifges(6,5) - t230;
t150 = Ifges(6,2) * t241 - t230;
t511 = t136 + t150;
t307 = Ifges(5,4) * t319;
t351 = -Ifges(5,2) * t315 + t307;
t288 = Ifges(5,1) * t315 + t307;
t509 = -t225 * t315 + t226 * t319;
t418 = t320 * mrSges(4,2);
t286 = t316 * mrSges(4,1) + t418;
t508 = -mrSges(5,1) * t319 + mrSges(5,2) * t315;
t253 = t316 * t399 - t416 * t320;
t132 = t346 * t253;
t133 = t345 * t253;
t342 = t132 * t524 + t133 * t529;
t484 = t242 / 0.2e1;
t485 = t241 / 0.2e1;
t429 = t242 * mrSges(6,3);
t201 = mrSges(6,2) * t320 - t429;
t489 = t201 / 0.2e1;
t505 = t359 * t489 + (t359 * t484 + t93 * t485) * mrSges(6,3);
t147 = -mrSges(6,1) * t241 - mrSges(6,2) * t242;
t475 = t303 / 0.2e1;
t504 = t147 * t475 + t358 * t489;
t449 = Ifges(6,4) * t241;
t134 = -Ifges(6,2) * t242 - t320 * Ifges(6,6) - t449;
t151 = -Ifges(6,1) * t242 + t449;
t376 = -Ifges(6,1) / 0.4e1 + Ifges(6,2) / 0.4e1;
t503 = -t376 * t242 - t151 / 0.4e1 + t134 / 0.4e1;
t494 = t93 / 0.2e1;
t490 = t179 / 0.2e1;
t430 = t241 * mrSges(6,3);
t203 = -mrSges(6,1) * t320 + t430;
t488 = -t203 / 0.2e1;
t487 = t203 / 0.2e1;
t486 = -t241 / 0.2e1;
t483 = -t242 / 0.2e1;
t482 = -t243 / 0.2e1;
t481 = -t244 / 0.2e1;
t480 = t253 / 0.2e1;
t256 = t508 * t316;
t479 = -t256 / 0.2e1;
t467 = pkin(4) * t315;
t366 = pkin(7) + t467;
t277 = t366 * t316;
t477 = t277 / 0.2e1;
t419 = t319 * mrSges(5,2);
t423 = t315 * mrSges(5,1);
t285 = t419 + t423;
t476 = t285 / 0.2e1;
t474 = -t315 / 0.2e1;
t473 = t315 / 0.2e1;
t471 = -t319 / 0.2e1;
t470 = t319 / 0.2e1;
t469 = -t320 / 0.2e1;
t468 = -t320 / 0.4e1;
t465 = pkin(7) * t316;
t464 = t81 * mrSges(6,2);
t463 = t82 * mrSges(6,1);
t456 = t97 * mrSges(6,1);
t455 = t98 * mrSges(6,2);
t454 = Ifges(6,1) - Ifges(6,2);
t450 = Ifges(5,4) * t315;
t448 = Ifges(6,4) * t345;
t441 = pkin(4) * qJD(4);
t434 = t358 * mrSges(6,3);
t428 = t345 * mrSges(6,2);
t427 = t345 * mrSges(6,3);
t426 = t346 * mrSges(6,1);
t424 = t314 * t91;
t421 = t316 * mrSges(4,2);
t420 = t318 * t88;
t417 = t508 - mrSges(4,1);
t415 = t116 * t346;
t414 = t117 * t345;
t413 = t179 * t315;
t412 = t180 * t319;
t407 = t224 * t319;
t170 = t253 * t254;
t24 = m(6) * (t132 * t359 + t133 * t93 + t170) + m(5) * (t170 + (-t412 + t413) * t253);
t404 = t24 * qJD(1);
t200 = t253 * t370;
t25 = m(6) * (t116 * t359 + t117 * t93 + t200) + m(5) * (t179 * t219 + t180 * t220 + t200) + m(4) * (t253 * t316 + t254 * t320 - t399) * t398;
t403 = t25 * qJD(1);
t402 = t253 * t315;
t401 = t277 * t315;
t202 = -mrSges(6,2) * t316 - mrSges(6,3) * t243;
t397 = t314 * t202;
t396 = t314 * t241;
t236 = -Ifges(5,6) * t320 + t351 * t316;
t395 = t315 * t236;
t204 = mrSges(6,1) * t316 + mrSges(6,3) * t244;
t391 = t318 * t204;
t390 = t318 * t242;
t352 = Ifges(5,1) * t319 - t450;
t238 = -Ifges(5,5) * t320 + t352 * t316;
t388 = t319 * t238;
t353 = t426 - t428;
t385 = t199 * t425 - t303 * t353;
t384 = -Ifges(6,5) * t242 + Ifges(6,6) * t241;
t309 = t315 ^ 2;
t311 = t319 ^ 2;
t382 = -t309 - t311;
t381 = pkin(4) * t392;
t379 = mrSges(5,3) * t394;
t375 = Ifges(5,2) / 0.4e1 - Ifges(5,1) / 0.4e1;
t374 = -t81 / 0.2e1 + t98 / 0.2e1;
t373 = t82 / 0.2e1 + t97 / 0.2e1;
t371 = t494 - t93 / 0.2e1;
t365 = t147 * t480;
t148 = mrSges(6,1) * t242 - mrSges(6,2) * t241;
t362 = t148 * t473;
t287 = Ifges(5,2) * t319 + t450;
t361 = t287 * t474;
t360 = t306 - t444;
t355 = t383 * t468;
t350 = Ifges(5,5) * t315 + Ifges(5,6) * t319;
t178 = mrSges(6,1) * t345 + mrSges(6,2) * t346;
t149 = mrSges(6,1) * t243 - mrSges(6,2) * t244;
t223 = t266 - t380;
t257 = t285 * t316;
t258 = t285 * t320;
t273 = mrSges(5,2) * t320 - t379;
t274 = -mrSges(5,2) * t316 - mrSges(5,3) * t393;
t275 = -mrSges(5,1) * t320 - mrSges(5,3) * t392;
t276 = mrSges(5,1) * t316 - mrSges(5,3) * t387;
t278 = t366 * t320;
t323 = (t148 / 0.2e1 + t257 / 0.2e1) * t254 + (t273 * t471 + t275 * t473 + t149 / 0.2e1 + t258 / 0.2e1) * t253 + (t254 * t465 + t179 * t225 + t180 * t226 + (pkin(7) * t320 + t223 * t315 - t407) * t253) * t501 + (t132 * t81 + t133 * t82 + t253 * t278 + t254 * t277 + t359 * t88 + t91 * t93) * t500 + t132 * t487 + t133 * t489 + t276 * t490 + t180 * t274 / 0.2e1 + t359 * t204 / 0.2e1 + t202 * t494;
t4 = (t414 / 0.2e1 + t415 / 0.2e1) * mrSges(6,3) + (-t286 / 0.2e1 + t418 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t508 / 0.2e1 - t178 / 0.2e1) * t316) * t398 + t323 + t518;
t135 = -Ifges(6,4) * t244 - Ifges(6,2) * t243 + Ifges(6,6) * t316;
t137 = -Ifges(6,1) * t244 - Ifges(6,4) * t243 + Ifges(6,5) * t316;
t237 = Ifges(5,6) * t316 + t351 * t320;
t239 = Ifges(5,5) * t316 + t352 * t320;
t7 = -pkin(2) * t286 + t226 * t273 + t224 * t274 + t225 * t275 + t223 * t276 + t277 * t149 + t278 * t148 + t135 * t483 + t134 * t482 + t136 * t481 + t137 * t486 + t81 * t204 + t91 * t201 + t82 * t202 + t88 * t203 + m(5) * (t223 * t225 + t224 * t226) + m(6) * (t277 * t278 + t81 * t88 + t82 * t91) + (t388 / 0.2e1 - t395 / 0.2e1 + pkin(7) * t257 + t343 + t513 * t320) * t320 + (t239 * t470 + t237 * t474 + pkin(7) * t258 + Ifges(6,5) * t486 + Ifges(6,6) * t483 - t513 * t316 + (m(5) * pkin(7) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) - Ifges(6,3)) * t320) * t316;
t349 = t4 * qJD(1) + t7 * qJD(2);
t339 = t277 * t147 + t384 * t469 + t82 * t430;
t10 = t256 * t465 - t98 * t201 - t97 * t203 - m(6) * (t81 * t97 + t82 * t98) - t148 * t381 + t224 * t275 + (t151 / 0.2e1 - t134 / 0.2e1) * t241 - (-t136 / 0.2e1 - t150 / 0.2e1 + t81 * mrSges(6,3)) * t242 + (t238 * t473 + t236 * t470 - m(6) * t277 * t466 + t350 * t469 + mrSges(5,3) * t407 + (-t288 * t471 + t361) * t316) * t316 - t339 + (-t379 - t273) * t223;
t325 = (-t412 / 0.2e1 + t413 / 0.2e1) * t316 * mrSges(5,3) + (t253 * t381 + t522 * t359) * t500 + t273 * t490 - t180 * t275 / 0.2e1 + t505 + (t528 * t500 - t487) * t93;
t6 = (t147 / 0.2e1 + t479) * t253 + t325 + t530;
t348 = t6 * qJD(1) - t10 * qJD(2);
t329 = t93 * t488 + t365 + t505;
t13 = t329 - t512;
t15 = -t203 * t82 + t151 * t486 + t134 * t485 + (t201 + t429) * t81 + t339 + t511 * t483;
t347 = t13 * qJD(1) + t15 * qJD(2);
t341 = t225 * t497 + t226 * t496;
t340 = t288 * t470 + t361;
t338 = (t132 * t318 + t133 * t314) * t498;
t336 = t353 * t480;
t322 = (t273 * t474 + t275 * t471) * pkin(8) + (t199 * t485 - t345 * t374 - t346 * t373 + t358 * t484) * mrSges(6,3) - (-t449 / 0.2e1 + t277 * t524 + t503) * t346 - (-t230 / 0.2e1 + t136 / 0.4e1 + t150 / 0.4e1 + mrSges(6,2) * t477 + t376 * t241) * t345 + pkin(3) * t256 / 0.2e1 - t199 * t487 - t395 / 0.4e1 + t388 / 0.4e1 + t504;
t326 = pkin(7) * t476 + (-t311 / 0.2e1 - t309 / 0.2e1) * pkin(8) * mrSges(5,3) + (-t288 / 0.4e1 - t307 / 0.4e1 + t375 * t315) * t315 + (-0.3e1 / 0.4e1 * t450 - t287 / 0.4e1 - t375 * t319 + (m(6) * t475 + t178 / 0.2e1) * pkin(4)) * t319;
t333 = t528 * t199 + t522 * t358;
t1 = t333 * t500 + (-t391 / 0.2e1 + t362 - t397 / 0.2e1 + 0.2e1 * (-t424 / 0.4e1 - t420 / 0.4e1 + t401 / 0.4e1) * m(6)) * pkin(4) + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t326) * t316 + t322 + (t369 + 0.3e1 / 0.4e1 * t444 - t306 / 0.4e1 + t260 / 0.4e1 + t259 / 0.4e1) * t320 + t341 + t507;
t327 = pkin(4) * t402 * t500 - t371 * t425;
t12 = (t426 / 0.2e1 - t428 / 0.2e1 + t476 - t419 / 0.2e1 - t423 / 0.2e1) * t253 - t338 / 0.2e1 + t327 + t342;
t18 = pkin(3) * t285 + t352 * t474 + t351 * t471 - t525 * t346 - (-t346 * t454 + t434 + t448) * t345 - t340 + t385 + (-m(6) * t303 - t178) * t467 + t358 * t427;
t335 = t12 * qJD(1) + t1 * qJD(2) - t18 * qJD(3);
t16 = t336 + t342;
t28 = -Ifges(6,4) * t345 ^ 2 - (-t345 * t454 + t525) * t346 + t385;
t324 = -(-t448 / 0.2e1 - t434 / 0.2e1) * t242 + (t447 / 0.2e1 + t431 / 0.2e1 - t376 * t345) * t241 + t199 * t488 + t353 * t477 + t355 - t511 * t345 / 0.4e1 - t503 * t346 + t504;
t9 = t324 - t502;
t334 = t16 * qJD(1) + t9 * qJD(2) - t28 * qJD(3);
t17 = t336 - t342;
t330 = (t318 * t489 + t314 * t488 + (t396 / 0.2e1 + t390 / 0.2e1) * mrSges(6,3)) * pkin(4);
t23 = -t373 * mrSges(6,1) + t374 * mrSges(6,2) + t330;
t26 = t371 * mrSges(6,1);
t281 = (mrSges(6,1) * t314 + mrSges(6,2) * t318) * pkin(4);
t332 = qJD(1) * t26 - qJD(2) * t23 + qJD(4) * t281;
t312 = t320 ^ 2;
t310 = t316 ^ 2;
t283 = t310 * pkin(7) * t398;
t267 = t281 * qJD(5);
t19 = -t464 / 0.2e1 - t463 / 0.2e1 - t455 / 0.2e1 + t456 / 0.2e1 + t330 + t384;
t14 = t329 + t512;
t11 = t253 * t476 + t338 / 0.2e1 + t419 * t480 + mrSges(5,1) * t402 / 0.2e1 + t327 + t17;
t8 = t324 + t502;
t5 = t253 * t479 + t325 + t365 - t530;
t3 = t323 + (-t414 - t415) * mrSges(6,3) / 0.2e1 + (t508 + t178) * t370 / 0.2e1 - t286 * t398 - t518;
t2 = Ifges(5,5) * t387 / 0.2e1 - Ifges(5,6) * t393 / 0.2e1 + pkin(4) * t362 + t326 * t316 + t360 * t468 + t355 + Ifges(5,3) * t472 + (t420 + t424) * t498 / 0.2e1 + t322 + (pkin(4) * t401 + t333) * t500 - t341 + (t391 + t397) * pkin(4) / 0.2e1 + t502;
t20 = [t25 * qJD(2) + t24 * qJD(3), t3 * qJD(3) + t5 * qJD(4) + t14 * qJD(5) + t403 + (t116 * t203 + t117 * t201 + t219 * t275 + t220 * t273 + ((-t320 * mrSges(4,1) - mrSges(3,1) + t421) * t317 + (-mrSges(3,2) + (t148 + t257) * t316 + (t310 + t312) * mrSges(4,3)) * t321) * t313 + (t81 * t116 + t82 * t117 + t277 * t370) * t517 + (t219 * t223 + t220 * t224 + t283) * t516 + m(4) * (t283 + (pkin(7) * t312 * t321 - pkin(2) * t317) * t313)) * qJD(2), t3 * qJD(2) + t11 * qJD(4) + t17 * qJD(5) + t404 + ((-t132 * t346 - t133 * t345) * mrSges(6,3) + (t178 + t417) * t254 + (t382 * mrSges(5,3) + mrSges(4,2)) * t253 + (t132 * t358 + t133 * t199 + t254 * t303) * t517 + (t382 * t253 * pkin(8) - pkin(3) * t254) * t516) * qJD(3), t5 * qJD(2) + t11 * qJD(3) + (-t180 * mrSges(5,1) - t179 * mrSges(5,2) + (t314 * t359 - t318 * t93) * t498 + t27) * qJD(4) + t527, t14 * qJD(2) + t17 * qJD(3) + t27 * qJD(4) + t527; qJD(3) * t4 + qJD(4) * t6 + qJD(5) * t13 - t403, qJD(3) * t7 - qJD(4) * t10 + qJD(5) * t15, t2 * qJD(4) + t8 * qJD(5) + t349 + (pkin(7) * t421 + m(6) * (t199 * t91 + t278 * t303 + t358 * t88) - t91 * t427 - t88 * t425 + t237 * t470 + t239 * t473 - Ifges(4,6) * t316 + t303 * t149 + t278 * t178 - t345 * t135 / 0.2e1 + (-Ifges(6,2) * t345 + t447) * t482 + (Ifges(6,1) * t346 - t448) * t481 + t346 * t137 / 0.2e1 - pkin(3) * t258 + t358 * t204 + t199 * t202 + (m(5) * t509 + t319 * t274 - t315 * t276) * pkin(8) + (Ifges(4,5) + (-m(5) * pkin(3) + t417) * pkin(7) + t340) * t320 + (Ifges(6,5) * t346 - Ifges(6,6) * t345 + t350) * t472 + t509 * mrSges(5,3)) * qJD(3), t2 * qJD(3) + (-t224 * mrSges(5,1) - t223 * mrSges(5,2) - Ifges(5,5) * t394 - Ifges(5,6) * t392 + t384 - t455 + t456) * qJD(4) + t19 * qJD(5) + (m(6) * (t314 * t98 + t318 * t97) + (t390 + t396) * mrSges(6,3)) * t441 + t348, t8 * qJD(3) + t19 * qJD(4) + (t384 - t463 - t464) * qJD(5) + t347; -qJD(2) * t4 + qJD(4) * t12 + qJD(5) * t16 - t404, qJD(4) * t1 + qJD(5) * t9 - t349, -qJD(4) * t18 - qJD(5) * t28, (t508 * pkin(8) + t36 + t360) * qJD(4) + t526 + (m(6) * (-t199 * t318 + t314 * t358) + (-t314 * t346 + t318 * t345) * mrSges(6,3)) * t441 + t335, t36 * qJD(4) + t334 + t526; -qJD(2) * t6 - qJD(3) * t12 - qJD(5) * t26, -qJD(3) * t1 + qJD(5) * t23 - t348, -t335, -t267, -t267 - t332; -t13 * qJD(2) - t16 * qJD(3) + t26 * qJD(4), -qJD(3) * t9 - qJD(4) * t23 - t347, -t334, t332, 0;];
Cq = t20;
