% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:41
% EndTime: 2019-03-08 21:12:52
% DurationCPUTime: 6.21s
% Computational Cost: add. (8716->566), mult. (21261->813), div. (0->0), fcn. (21604->10), ass. (0->282)
t317 = cos(pkin(11));
t322 = cos(qJ(3));
t315 = sin(pkin(11));
t316 = sin(pkin(6));
t323 = cos(qJ(2));
t397 = t316 * t323;
t382 = t315 * t397;
t320 = sin(qJ(2));
t398 = t316 * t320;
t195 = -t317 * t398 + t322 * t382;
t394 = t317 * t323;
t196 = (t315 * t320 + t322 * t394) * t316;
t318 = sin(qJ(6));
t321 = cos(qJ(6));
t110 = t195 * t321 - t196 * t318;
t111 = t195 * t318 + t196 * t321;
t445 = -pkin(9) + qJ(4);
t279 = t445 * t315;
t282 = t445 * t317;
t178 = t279 * t321 - t282 * t318;
t179 = t279 * t318 + t282 * t321;
t367 = qJ(5) * t315 + pkin(3);
t466 = pkin(4) + pkin(5);
t251 = t317 * t466 + t367;
t276 = -pkin(4) * t317 - t367;
t413 = t196 * t317;
t414 = t195 * t315;
t340 = (t413 + t414) * qJ(4);
t319 = sin(qJ(3));
t381 = t319 * t397;
t472 = -m(7) / 0.2e1;
t474 = -m(6) / 0.2e1;
t476 = -m(5) / 0.2e1;
t487 = mrSges(6,2) + mrSges(5,3);
t491 = -t487 * (t413 / 0.2e1 + t414 / 0.2e1) + (-pkin(3) * t381 + t340) * t476 + (t276 * t381 + t340) * t474 + (t178 * t110 + t179 * t111 - t251 * t381) * t472;
t312 = t317 ^ 2;
t385 = t315 ^ 2 + t312;
t467 = -mrSges(7,3) / 0.2e1;
t489 = m(6) + m(5);
t261 = t315 * t321 - t318 * t317;
t236 = t261 * t322;
t460 = t236 / 0.2e1;
t350 = t318 * t315 + t321 * t317;
t238 = t350 * t322;
t458 = t238 / 0.2e1;
t488 = -t319 / 0.2e1;
t486 = Ifges(6,4) + Ifges(5,5);
t485 = Ifges(5,6) - Ifges(6,6);
t484 = t315 * t466;
t401 = t315 * t319;
t424 = qJ(4) * t322;
t284 = pkin(3) * t319 - t424;
t403 = t284 * t317;
t214 = pkin(8) * t401 + t403;
t396 = t317 * t319;
t215 = -pkin(8) * t396 + t315 * t284;
t482 = -t214 * t315 + t215 * t317;
t185 = t319 * qJ(5) + t215;
t379 = -pkin(8) * t315 - pkin(4);
t190 = t319 * t379 - t403;
t481 = t185 * t317 + t190 * t315;
t480 = t474 + t476;
t479 = mrSges(6,3) / 0.2e1 - mrSges(5,2) / 0.2e1;
t263 = t322 * mrSges(5,2) - mrSges(5,3) * t401;
t265 = -t322 * mrSges(5,1) - mrSges(5,3) * t396;
t266 = t322 * mrSges(6,1) + mrSges(6,2) * t396;
t270 = -mrSges(6,2) * t401 - t322 * mrSges(6,3);
t478 = -(t266 / 0.2e1 - t265 / 0.2e1) * t315 - (t270 / 0.2e1 + t263 / 0.2e1) * t317;
t477 = 2 * qJD(3);
t475 = m(5) / 0.2e1;
t473 = m(6) / 0.2e1;
t471 = m(7) / 0.2e1;
t470 = m(7) / 0.4e1;
t469 = mrSges(7,1) / 0.2e1;
t468 = -mrSges(7,2) / 0.2e1;
t465 = Ifges(7,4) * t458 + Ifges(7,2) * t460 + Ifges(7,6) * t488;
t428 = t322 * mrSges(7,2);
t235 = t261 * t319;
t438 = t235 * mrSges(7,3);
t186 = -t428 + t438;
t464 = t186 / 0.2e1;
t430 = t322 * mrSges(7,1);
t237 = t350 * t319;
t436 = t237 * mrSges(7,3);
t188 = t430 - t436;
t463 = -t188 / 0.2e1;
t462 = t188 / 0.2e1;
t461 = -t235 / 0.2e1;
t459 = -t237 / 0.2e1;
t426 = cos(pkin(6));
t249 = t319 * t398 - t322 * t426;
t457 = -t249 / 0.2e1;
t456 = t350 / 0.2e1;
t455 = t261 / 0.2e1;
t454 = -t315 / 0.2e1;
t453 = t315 / 0.2e1;
t452 = t317 / 0.2e1;
t451 = -t318 / 0.2e1;
t450 = t319 / 0.2e1;
t449 = -t321 / 0.2e1;
t448 = m(6) * t276;
t447 = pkin(4) * t315;
t446 = pkin(8) * t322;
t444 = Ifges(5,4) * t315;
t443 = Ifges(5,4) * t317;
t442 = Ifges(7,4) * t237;
t441 = Ifges(7,4) * t261;
t440 = Ifges(6,5) * t315;
t439 = Ifges(6,5) * t317;
t437 = t236 * mrSges(7,1);
t435 = t238 * mrSges(7,2);
t434 = t318 * mrSges(7,1);
t433 = t319 * mrSges(6,1);
t432 = t319 * mrSges(4,2);
t431 = t321 * mrSges(7,2);
t429 = t322 * mrSges(4,2);
t281 = -mrSges(5,1) * t317 + mrSges(5,2) * t315;
t427 = t281 - mrSges(4,1);
t423 = t110 * t261;
t422 = t111 * t350;
t126 = t261 * t249;
t421 = t126 * t321;
t127 = t350 * t249;
t420 = t127 * t318;
t250 = t319 * t426 + t322 * t398;
t172 = t250 * t315 + t316 * t394;
t407 = t249 * t315;
t408 = t249 * t250;
t173 = t250 * t317 - t382;
t417 = t173 * t317;
t82 = t172 * t321 - t173 * t318;
t83 = t172 * t318 + t173 * t321;
t15 = m(7) * (-t126 * t82 - t127 * t83 + t408) + t489 * (-t172 * t407 - t249 * t417 + t408);
t419 = t15 * qJD(1);
t406 = t249 * t319;
t17 = m(7) * (t82 * t110 + t83 * t111) + 0.4e1 * (t406 * t470 + m(4) * (t250 * t322 + t406) / 0.4e1 - m(4) * t398 / 0.4e1) * t397 + t489 * (t172 * t195 + t173 * t196 + t249 * t381);
t418 = t17 * qJD(1);
t410 = t235 * t318;
t409 = t237 * t321;
t405 = t350 * t318;
t404 = t261 * t321;
t30 = (t428 / 0.2e1 - t186 / 0.2e1 + t438 / 0.2e1) * t321 + (t430 / 0.2e1 + t436 / 0.2e1 + t462) * t318;
t402 = t30 * qJD(2);
t399 = t315 * t322;
t395 = t317 * t322;
t393 = t318 * t188;
t391 = t321 * t186;
t134 = -mrSges(7,1) * t235 + mrSges(7,2) * t237;
t360 = t315 * mrSges(6,1) - t317 * mrSges(6,3);
t252 = t360 * t319;
t390 = t134 - t252;
t389 = t385 * qJ(4) * t249;
t388 = Ifges(7,5) * t235 - Ifges(7,6) * t237;
t387 = -Ifges(7,5) * t350 - Ifges(7,6) * t261;
t162 = mrSges(7,1) * t350 + mrSges(7,2) * t261;
t280 = -mrSges(6,1) * t317 - mrSges(6,3) * t315;
t386 = t280 - t162;
t278 = -pkin(3) * t322 - qJ(4) * t319 - pkin(2);
t198 = pkin(8) * t395 + t315 * t278;
t384 = t474 + t472;
t383 = m(6) / 0.4e1 + m(5) / 0.4e1;
t380 = t249 * t396;
t161 = t261 * mrSges(7,1) - mrSges(7,2) * t350;
t376 = t161 * t457;
t375 = -t407 / 0.2e1;
t372 = -t396 / 0.2e1;
t371 = t134 / 0.2e1 - t252 / 0.2e1;
t370 = t162 / 0.2e1 - t280 / 0.2e1;
t366 = qJ(5) * t317 - pkin(8);
t297 = pkin(8) * t399;
t197 = t278 * t317 - t297;
t183 = -qJ(5) * t322 + t198;
t362 = t370 * t317;
t361 = t315 * mrSges(5,1) + t317 * mrSges(5,2);
t133 = t237 * mrSges(7,1) + t235 * mrSges(7,2);
t285 = t319 * mrSges(4,1) + t429;
t135 = t435 - t437;
t180 = (t366 - t484) * t319;
t295 = qJ(5) * t395;
t181 = -t295 + (pkin(8) + t484) * t322;
t187 = mrSges(7,2) * t319 + mrSges(7,3) * t236;
t189 = -mrSges(7,1) * t319 - mrSges(7,3) * t238;
t231 = (-t366 + t447) * t319;
t232 = -t295 + (pkin(8) + t447) * t322;
t253 = t361 * t319;
t254 = t360 * t322;
t255 = t361 * t322;
t264 = -mrSges(5,2) * t319 - mrSges(5,3) * t399;
t267 = mrSges(5,1) * t319 - mrSges(5,3) * t395;
t296 = mrSges(6,2) * t395;
t268 = t296 - t433;
t269 = -mrSges(6,2) * t399 + mrSges(6,3) * t319;
t351 = -t197 * t315 + t198 * t317;
t310 = t322 * pkin(4);
t184 = -t197 + t310;
t352 = t183 * t317 + t184 * t315;
t143 = pkin(5) * t322 + t297 + t310 + (-pkin(9) * t319 - t278) * t317;
t148 = pkin(9) * t401 + t183;
t68 = t143 * t321 - t148 * t318;
t69 = t143 * t318 + t148 * t321;
t144 = (-pkin(9) * t322 - t284) * t317 + (-pkin(5) + t379) * t319;
t157 = pkin(9) * t399 + t185;
t72 = t144 * t321 - t157 * t318;
t73 = t144 * t318 + t157 * t321;
t324 = (t264 / 0.2e1 + t269 / 0.2e1) * t173 + (-t267 / 0.2e1 + t268 / 0.2e1) * t172 + (t253 / 0.2e1 - t371) * t250 + (t254 / 0.2e1 + t255 / 0.2e1 - t135 / 0.2e1 + t478) * t249 + (pkin(8) * t250 * t319 - t172 * t214 + t173 * t215 + (-t351 + t446) * t249) * t475 + (t172 * t190 + t173 * t185 + t231 * t250 + (t232 - t352) * t249) * t473 + (-t126 * t68 - t127 * t69 - t180 * t250 + t181 * t249 + t72 * t82 + t73 * t83) * t471 - t126 * t462 - t127 * t464 + t82 * t189 / 0.2e1 + t83 * t187 / 0.2e1;
t2 = t324 + (-t285 / 0.2e1 + t429 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t281 / 0.2e1 + t370) * t319) * t397 + (t422 / 0.2e1 + t423 / 0.2e1) * mrSges(7,3) + t491;
t122 = Ifges(7,2) * t235 + t322 * Ifges(7,6) + t442;
t226 = Ifges(7,4) * t235;
t124 = Ifges(7,1) * t237 + t322 * Ifges(7,5) + t226;
t125 = Ifges(7,1) * t238 + Ifges(7,4) * t236 - t319 * Ifges(7,5);
t227 = t319 * Ifges(6,6) + (t315 * Ifges(6,3) + t439) * t322;
t228 = t319 * Ifges(5,6) + (-t315 * Ifges(5,2) + t443) * t322;
t229 = t319 * Ifges(6,4) + (Ifges(6,1) * t317 + t440) * t322;
t230 = t319 * Ifges(5,5) + (Ifges(5,1) * t317 - t444) * t322;
t346 = Ifges(7,5) * t458 + Ifges(7,6) * t460;
t3 = m(5) * (t197 * t214 + t198 * t215) + (Ifges(7,5) * t459 + Ifges(7,6) * t461 + pkin(8) * t255 - Ifges(4,4) * t319 + (t229 / 0.2e1 + t230 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t319) * t317 + (-t228 / 0.2e1 + t227 / 0.2e1 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t319) * t315) * t319 + t124 * t458 + t122 * t460 + t235 * t465 - pkin(2) * t285 + t185 * t270 + t198 * t264 + t214 * t265 + t190 * t266 + t197 * t267 + t184 * t268 + t183 * t269 + t215 * t263 + t232 * t252 + t231 * t254 + t237 * t125 / 0.2e1 + t68 * t189 + t73 * t186 + t69 * t187 + t72 * t188 + t180 * t135 - t181 * t134 + m(6) * (t183 * t185 + t184 * t190 + t231 * t232) + m(7) * (-t180 * t181 + t68 * t72 + t69 * t73) + (pkin(8) * t253 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(7,3) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t312 + ((Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t315 + (-Ifges(5,4) + Ifges(6,5)) * t317) * t315) * t319 + t346 + (t315 * t485 - t317 * t486 + Ifges(4,4)) * t322) * t322;
t359 = t2 * qJD(1) + t3 * qJD(2);
t331 = (t459 * t83 + t461 * t82) * mrSges(7,3) + t133 * t457 + t82 * t464 + t83 * t463;
t345 = t110 * t469 + t111 * t468;
t6 = t331 - t345;
t136 = -Ifges(7,2) * t237 + t226;
t137 = Ifges(7,1) * t235 - t442;
t8 = t180 * t133 - t69 * t188 + t68 * t186 + t322 * t388 / 0.2e1 + (-t69 * mrSges(7,3) + t137 / 0.2e1 - t122 / 0.2e1) * t237 + (-t68 * mrSges(7,3) + t124 / 0.2e1 + t136 / 0.2e1) * t235;
t358 = t6 * qJD(1) + t8 * qJD(2);
t19 = m(7) * (-t235 * t69 + t237 * t68) + t237 * t188 - t235 * t186 + ((-t265 + t266) * t317 + (-t263 - t270) * t315 + m(6) * (-t183 * t315 + t184 * t317) + m(5) * (-t197 * t317 - t198 * t315)) * t319;
t341 = (t475 - t384) * t397;
t343 = m(7) * (-t235 * t83 + t237 * t82);
t354 = t172 * t317 - t173 * t315;
t28 = -t343 / 0.2e1 + (t354 * t480 + t341) * t319;
t357 = -qJD(1) * t28 + qJD(2) * t19;
t22 = t390 * t396 + (-t270 - t391 + t393) * t322 + m(7) * (t180 * t396 + (t318 * t68 - t321 * t69) * t322) + m(6) * (-t183 * t322 - t231 * t396);
t330 = (-t173 * t322 - t380) * t474 + (-t380 + (t318 * t82 - t321 * t83) * t322) * t472;
t335 = t195 * t473 + (t321 * t110 + t318 * t111) * t471;
t24 = t330 + t335;
t356 = qJD(1) * t24 - qJD(2) * t22;
t355 = -qJD(2) * t133 - qJD(3) * t161;
t117 = m(6) * t315 + (t405 / 0.2e1 + t404 / 0.2e1 + t453) * m(7);
t85 = -m(6) * t396 + (t372 + t410 / 0.2e1 - t409 / 0.2e1) * m(7);
t349 = qJD(2) * t85 - qJD(3) * t117;
t347 = 0.2e1 * (t470 + t383) * t250;
t344 = -t126 * t469 - t127 * t468;
t11 = t376 - t344;
t258 = Ifges(7,4) * t350;
t163 = -Ifges(7,2) * t261 - t258;
t164 = -Ifges(7,2) * t350 + t441;
t165 = -Ifges(7,1) * t350 - t441;
t166 = Ifges(7,1) * t261 - t258;
t21 = t251 * t161 + (t165 / 0.2e1 - t164 / 0.2e1) * t261 - (t166 / 0.2e1 + t163 / 0.2e1) * t350;
t326 = -(t124 / 0.4e1 + t136 / 0.4e1) * t350 + (t137 / 0.4e1 - t122 / 0.4e1) * t261 + (t163 / 0.4e1 + t166 / 0.4e1 + t178 * t467) * t235 + (-t164 / 0.4e1 + t165 / 0.4e1 + t179 * t467) * t237 + t178 * t464 + t179 * t463 + t180 * t161 / 0.2e1 + t251 * t133 / 0.2e1 + t322 * t387 / 0.4e1;
t333 = Ifges(7,3) * t450 - t72 * mrSges(7,1) / 0.2e1 + t73 * mrSges(7,2) / 0.2e1 - t346;
t5 = t326 + t333;
t339 = t11 * qJD(1) + t5 * qJD(2) + t21 * qJD(3);
t325 = (t235 * t456 + t261 * t459) * mrSges(7,3) + t351 * t475 + t352 * t473 + (t178 * t237 - t179 * t235 + t261 * t68 + t350 * t69) * t471 + t186 * t456 + t188 * t455 - t478;
t334 = t232 * t474 + t181 * t472 - t437 / 0.2e1 + t435 / 0.2e1;
t10 = (pkin(8) * t476 + t479 * t317 + (-mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t315) * t322 + t325 + t334;
t329 = (t261 * t82 + t350 * t83) * t472 + t480 * (t172 * t315 + t417);
t26 = t347 + t329;
t37 = (-t261 ^ 2 - t350 ^ 2) * mrSges(7,3) + m(7) * (t178 * t261 + t179 * t350) + (0.4e1 * qJ(4) * t383 + t487) * t385;
t338 = -qJD(1) * t26 + qJD(2) * t10 + qJD(3) * t37;
t327 = t371 * t315 + (t261 * t451 - t350 * t449) * t322 * mrSges(7,3) + (-t231 * t315 + (-t276 * t319 - t424) * t317) * t473 + (t251 * t396 + t180 * t315 + (t178 * t318 - t179 * t321) * t322) * t471;
t332 = t190 * t474 + (t318 * t73 + t321 * t72) * t472 + t187 * t451 + t189 * t449;
t14 = -t296 + (mrSges(6,1) / 0.2e1 + t362) * t319 + t327 + t332;
t43 = (t375 + t421 / 0.2e1 + t420 / 0.2e1) * m(7);
t60 = (m(7) * t251 - t386 - t448) * t315;
t337 = -qJD(1) * t43 - qJD(2) * t14 - qJD(3) * t60;
t314 = t322 ^ 2;
t313 = t319 ^ 2;
t283 = t313 * pkin(8) * t397;
t128 = m(7) * t454 + (t404 + t405) * t471;
t96 = m(7) * t372 + (t409 - t410) * t471;
t42 = m(6) * t375 + (-t420 - t421) * t471 + t384 * t407;
t31 = t391 / 0.2e1 + t436 * t451 - t393 / 0.2e1 + t438 * t449 + (t434 / 0.2e1 + t431 / 0.2e1) * t322;
t29 = t343 / 0.2e1 + t319 * t341 + t489 * t354 * t450;
t27 = t347 - t329;
t25 = -t330 + t335;
t16 = -t433 / 0.2e1 + t319 * t362 + t327 - t332;
t12 = t376 + t344;
t9 = t446 * t475 + t325 - t334 - t479 * t395 + (mrSges(5,1) + mrSges(6,1)) * t399 / 0.2e1;
t7 = t331 + t345;
t4 = t326 - t333;
t1 = t324 + (t281 + t280) * t381 / 0.2e1 + (t422 + t423) * t467 - (t285 + t429 + (mrSges(4,1) + t162) * t319) * t397 / 0.2e1 - t491;
t13 = [t17 * qJD(2) + t15 * qJD(3), t1 * qJD(3) + t29 * qJD(4) + t25 * qJD(5) + t7 * qJD(6) + t418 + (t110 * t188 + t111 * t186 - t195 * t265 + t195 * t266 + t196 * t263 + t196 * t270 + ((-t322 * mrSges(4,1) - mrSges(3,1) + t432) * t320 + (-mrSges(3,2) + (t313 + t314) * mrSges(4,3) + (t253 - t390) * t319) * t323) * t316 + 0.2e1 * (t183 * t196 + t184 * t195 + t231 * t381) * t473 + 0.2e1 * (t68 * t110 + t69 * t111 - t180 * t381) * t471 + 0.2e1 * (-t195 * t197 + t196 * t198 + t283) * t475 + m(4) * (t283 + (pkin(8) * t314 * t323 - pkin(2) * t320) * t316)) * qJD(2), t419 + t1 * qJD(2) + t27 * qJD(4) + t42 * qJD(5) + t12 * qJD(6) + ((t250 * t276 - t389) * t473 + (-t126 * t178 - t127 * t179 - t250 * t251) * t471 + (-pkin(3) * t250 - t389) * t475) * t477 + ((t126 * t261 + t127 * t350) * mrSges(7,3) + (t386 + t427) * t250 + (-t385 * t487 + mrSges(4,2)) * t249) * qJD(3), qJD(2) * t29 + qJD(3) * t27, qJD(2) * t25 + qJD(3) * t42, t7 * qJD(2) + t12 * qJD(3) + (-mrSges(7,1) * t83 - mrSges(7,2) * t82) * qJD(6); qJD(3) * t2 - qJD(4) * t28 - qJD(5) * t24 + qJD(6) * t6 - t418, qJD(3) * t3 + qJD(4) * t19 + qJD(5) * t22 + qJD(6) * t8, t9 * qJD(4) + t16 * qJD(5) + t4 * qJD(6) + ((t178 * t72 + t179 * t73 - t181 * t251) * t471 + t232 * t448 / 0.2e1) * t477 + t359 + (t164 * t460 - t350 * t465 + t228 * t452 + t125 * t455 + t166 * t458 + pkin(8) * t432 - Ifges(4,6) * t319 + (Ifges(7,5) * t261 - Ifges(7,6) * t350) * t488 - t317 * t227 / 0.2e1 + t232 * t280 + t276 * t254 + t251 * t135 - pkin(3) * t255 + t178 * t189 + t179 * t187 - t181 * t162 + (t230 + t229) * t453 + (t315 * t486 + t317 * t485) * t450 + (-t72 * t261 - t350 * t73) * mrSges(7,3) + t482 * mrSges(5,3) + t481 * mrSges(6,2) + ((t264 + t269) * t317 + (-t267 + t268) * t315 + m(5) * t482 + m(6) * t481) * qJ(4) + ((Ifges(5,2) * t317 + t444) * t454 + (-Ifges(6,3) * t317 + t440) * t453 + Ifges(4,5) + (-m(5) * pkin(3) + t427) * pkin(8) + (-t439 + t443 + (Ifges(5,1) + Ifges(6,1)) * t315) * t452) * t322) * qJD(3), qJD(3) * t9 + qJD(5) * t96 + t357, qJD(3) * t16 + qJD(4) * t96 + qJD(6) * t31 - t356, t4 * qJD(3) + t31 * qJD(5) + (-mrSges(7,1) * t69 - mrSges(7,2) * t68 + t388) * qJD(6) + t358; -qJD(2) * t2 - qJD(4) * t26 + qJD(5) * t43 + qJD(6) * t11 - t419, qJD(4) * t10 + qJD(5) * t14 + qJD(6) * t5 - t359, qJD(4) * t37 + qJD(5) * t60 + qJD(6) * t21, qJD(5) * t128 + t338, qJD(4) * t128 - t337 (-mrSges(7,1) * t179 - mrSges(7,2) * t178 + t387) * qJD(6) + t339; qJD(2) * t28 + qJD(3) * t26, -qJD(3) * t10 + qJD(5) * t85 - qJD(6) * t133 - t357, -qJD(5) * t117 - qJD(6) * t161 - t338, 0, t349, t355; qJD(2) * t24 - qJD(3) * t43, -qJD(3) * t14 - qJD(4) * t85 - qJD(6) * t30 + t356, qJD(4) * t117 + t337, -t349, 0, -t402 + (-t431 - t434) * qJD(6); -t6 * qJD(2) - t11 * qJD(3), -qJD(3) * t5 + qJD(4) * t133 + qJD(5) * t30 - t358, qJD(4) * t161 - t339, -t355, t402, 0;];
Cq  = t13;
