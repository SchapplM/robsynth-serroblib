% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:31:55
% EndTime: 2018-11-23 17:32:20
% DurationCPUTime: 25.35s
% Computational Cost: add. (11000->757), mult. (30149->986), div. (0->0), fcn. (21545->8), ass. (0->359)
t300 = sin(qJ(3));
t301 = sin(qJ(2));
t302 = cos(qJ(3));
t391 = qJD(3) * t302;
t375 = t301 * t391;
t384 = qJD(2) * qJD(3);
t303 = cos(qJ(2));
t393 = qJD(2) * t303;
t207 = -t300 * t384 + (-t300 * t393 - t375) * qJD(1);
t297 = sin(pkin(6));
t299 = cos(pkin(6));
t395 = qJD(2) * t301;
t367 = qJD(1) * t395;
t171 = -t207 * t297 + t299 * t367;
t466 = t171 / 0.2e1;
t392 = qJD(3) * t300;
t376 = t301 * t392;
t206 = t302 * t384 + (t302 * t393 - t376) * qJD(1);
t296 = sin(pkin(10));
t298 = cos(pkin(10));
t360 = t297 * t367;
t316 = t207 * t299 + t360;
t119 = t206 * t298 + t296 * t316;
t474 = t119 / 0.2e1;
t475 = -t119 / 0.2e1;
t413 = t298 * t299;
t118 = t206 * t296 - t207 * t413 - t298 * t360;
t476 = t118 / 0.2e1;
t477 = -t118 / 0.2e1;
t491 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t512 = Ifges(5,5) + Ifges(7,5);
t513 = Ifges(7,4) + Ifges(6,5);
t544 = Ifges(6,4) * t475 + Ifges(5,6) * t477 + t491 * t466 + t512 * t474 + t513 * t476;
t515 = Ifges(5,1) + Ifges(7,3);
t514 = Ifges(5,4) - Ifges(7,6);
t511 = Ifges(7,2) + Ifges(6,3);
t510 = -Ifges(6,6) + Ifges(7,6);
t543 = Ifges(6,6) + Ifges(5,4);
t355 = pkin(2) * t301 - pkin(9) * t303;
t268 = t355 * qJD(1);
t397 = qJD(1) * t301;
t378 = t300 * t397;
t208 = pkin(8) * t378 + t302 * t268;
t410 = t299 * t302;
t383 = qJ(4) * t410;
t320 = pkin(3) * t301 - t303 * t383;
t162 = qJD(1) * t320 + t208;
t428 = qJ(4) * t299;
t365 = pkin(9) + t428;
t344 = qJD(3) * t365;
t387 = qJD(4) * t300;
t374 = t299 * t387;
t214 = -t302 * t344 - t374;
t541 = t162 - t214;
t396 = qJD(1) * t303;
t295 = pkin(8) * t396;
t416 = t297 * t302;
t332 = pkin(3) * t300 - qJ(4) * t416;
t202 = t332 * t396 + t295;
t212 = qJD(3) * t332 - t297 * t387;
t540 = -t202 + t212;
t246 = t300 * t268;
t386 = qJD(4) * t302;
t373 = t299 * t386;
t407 = t301 * t302;
t408 = t300 * t303;
t417 = t297 * t301;
t527 = -t299 * t408 + t417;
t539 = t300 * t344 - t373 + t246 + (-pkin(8) * t407 + qJ(4) * t527) * qJD(1);
t538 = -t396 + qJD(3);
t536 = -t171 / 0.2e1;
t490 = -Ifges(6,4) + t512;
t489 = -Ifges(5,6) + t513;
t388 = qJD(4) * t298;
t377 = t302 * t397;
t264 = qJD(2) * t300 + t377;
t421 = t264 * t297;
t277 = -pkin(2) * t303 - t301 * pkin(9) - pkin(1);
t250 = t277 * qJD(1);
t282 = qJD(2) * pkin(9) + t295;
t184 = t302 * t250 - t300 * t282;
t149 = -t264 * t428 + t184;
t185 = t300 * t250 + t302 * t282;
t394 = qJD(2) * t302;
t263 = -t378 + t394;
t423 = t263 * t299;
t150 = -qJ(4) * t423 - t185;
t424 = t263 * t297;
t182 = pkin(3) * t264 - qJ(4) * t424;
t419 = t296 * t299;
t420 = t296 * t297;
t53 = t298 * t149 + t150 * t419 + t182 * t420;
t529 = qJ(5) * t421 - qJD(5) * t299 - t297 * t388 + t53;
t406 = t302 * t303;
t493 = t296 * t527 + t298 * t406;
t177 = t493 * qJD(1);
t228 = t298 * t391 - t392 * t419;
t401 = t177 - t228;
t414 = t298 * t264;
t333 = t538 * t297;
t198 = t333 + t423;
t426 = t198 * t296;
t148 = t414 + t426;
t197 = -t299 * t538 + t424;
t132 = qJ(4) * t198 + t185;
t138 = pkin(3) * t538 + t149;
t281 = -qJD(2) * pkin(2) + pkin(8) * t397;
t163 = -pkin(3) * t263 - qJ(4) * t421 + t281;
t40 = -t296 * t132 + t298 * (t138 * t299 + t163 * t297);
t307 = qJD(5) - t40;
t443 = pkin(4) + qJ(6);
t14 = pkin(5) * t148 + t197 * t443 + t307;
t422 = t264 * t296;
t147 = -t263 * t413 - t298 * t333 + t422;
t87 = -t138 * t297 + t299 * t163 + qJD(4);
t318 = -qJ(5) * t148 + t87;
t16 = t147 * t443 + t318;
t21 = pkin(4) * t147 + t318;
t22 = pkin(4) * t197 + t307;
t505 = -t514 * t147 + t148 * t515 - t512 * t197;
t63 = -t197 * Ifges(6,4) - t148 * Ifges(6,2) + t147 * Ifges(6,6);
t535 = t16 * mrSges(7,2) + t21 * mrSges(6,3) + t40 * mrSges(5,3) - t505 / 0.2e1 + t63 / 0.2e1 - t14 * mrSges(7,1) - t22 * mrSges(6,1) - t87 * mrSges(5,2);
t41 = t298 * t132 + t138 * t419 + t163 * t420;
t24 = qJ(5) * t197 - t41;
t17 = -pkin(5) * t147 + qJD(6) - t24;
t506 = t147 * t511 + t148 * t510 - t197 * t513;
t64 = t148 * Ifges(5,4) - t147 * Ifges(5,2) - t197 * Ifges(5,6);
t534 = t17 * mrSges(7,1) + t21 * mrSges(6,2) + t41 * mrSges(5,3) - t506 / 0.2e1 + t64 / 0.2e1 - t16 * mrSges(7,3) - t24 * mrSges(6,1) - t87 * mrSges(5,1);
t533 = -Ifges(4,6) / 0.2e1;
t447 = -t300 / 0.2e1;
t498 = -t298 * t539 - t419 * t541 + t420 * t540;
t532 = -t298 * (t162 * t299 + t202 * t297) + t539 * t296;
t172 = t263 * t296 + t264 * t413;
t173 = t263 * t298 - t264 * t419;
t97 = -t150 * t297 + t299 * t182;
t326 = -qJ(5) * t173 + t97;
t385 = qJD(5) * t296;
t531 = -t172 * t443 + (-qJD(6) * t298 - t385) * t297 - t326;
t389 = qJD(4) * t297;
t145 = t296 * t149;
t403 = -pkin(4) * t421 + t145;
t530 = -qJD(6) * t299 + t296 * t389 + t150 * t413 - pkin(5) * t173 - (-qJ(6) * t264 - t182 * t298) * t297 - t403;
t528 = pkin(5) * t172 - t529;
t495 = t297 * t541 + t299 * t540;
t412 = t298 * t300;
t329 = t296 * t302 + t299 * t412;
t494 = -t298 * t417 + t329 * t303;
t176 = t494 * qJD(1);
t227 = t329 * qJD(3);
t402 = t176 - t227;
t334 = -qJ(4) * t206 - qJD(4) * t264;
t112 = -t207 * pkin(3) + qJD(2) * t295 + t297 * t334;
t269 = t355 * qJD(2);
t251 = qJD(1) * t269;
t363 = pkin(8) * t367;
t131 = -qJD(3) * t185 + t302 * t251 + t300 * t363;
t77 = pkin(3) * t367 + t299 * t334 + t131;
t26 = t299 * t112 - t297 * t77;
t310 = -qJ(5) * t119 - qJD(5) * t148 + t26;
t1 = qJD(6) * t147 + t118 * t443 + t310;
t130 = t250 * t391 + t300 * t251 - t282 * t392 - t302 * t363;
t69 = qJ(4) * t316 + qJD(4) * t198 + t130;
t11 = t112 * t420 + t298 * t69 + t77 * t419;
t4 = -qJ(5) * t171 + qJD(5) * t197 - t11;
t3 = -pkin(5) * t118 - t4;
t5 = pkin(4) * t118 + t310;
t517 = t26 * mrSges(5,1) + t4 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) - t11 * mrSges(5,3) + t1 * mrSges(7,3) + 0.2e1 * t511 * t476 + (t510 - t514) * t474 + (t489 + t513) * t466 + Ifges(5,6) * t536 + t543 * t475 + (-t477 + t476) * Ifges(5,2);
t10 = -t296 * t69 + t298 * (t112 * t297 + t299 * t77);
t2 = pkin(5) * t119 + qJD(6) * t197 - t171 * t443 - t10;
t8 = -pkin(4) * t171 - t10;
t526 = t10 * mrSges(5,1) - t11 * mrSges(5,2) + t8 * mrSges(6,2) + t3 * mrSges(7,2) - t4 * mrSges(6,3) - t2 * mrSges(7,3) + t544;
t411 = t299 * t301;
t327 = t297 * t408 + t411;
t230 = t327 * qJD(1);
t501 = qJ(5) * t230 - t297 * (qJ(5) * t392 - qJD(5) * t302) - t498;
t238 = t296 * t410 + t412;
t525 = qJ(5) * t401 - qJD(5) * t238 + t495;
t524 = -t214 * t413 - t532;
t523 = Ifges(4,5) * t447 + t302 * t533;
t521 = -t40 * mrSges(5,1) + t41 * mrSges(5,2) - t22 * mrSges(6,2) - t17 * mrSges(7,2) + t24 * mrSges(6,3) + t14 * mrSges(7,3);
t519 = t526 + t544;
t518 = t510 * t476 + 0.2e1 * t515 * t474 + Ifges(6,4) * t536 - t1 * mrSges(7,2) + t26 * mrSges(5,2) - t5 * mrSges(6,3) + t8 * mrSges(6,1) - t10 * mrSges(5,3) + t2 * mrSges(7,1) + (t512 + t490) * t466 + (t474 - t475) * Ifges(6,2) + (t514 + t543) * t477;
t371 = -Ifges(3,6) * qJD(2) / 0.2e1;
t372 = Ifges(3,5) * qJD(2) / 0.2e1;
t237 = t296 * t300 - t298 * t410;
t507 = qJD(6) * t237 - t402 * t443 + t525;
t425 = t212 * t298;
t504 = t230 * t443 + (qJD(6) * t302 - t392 * t443 - t425) * t297 + t524 - t401 * pkin(5);
t503 = pkin(5) * t402 - t501;
t502 = -pkin(4) * t402 + t525;
t500 = pkin(4) * t230 + (-pkin(4) * t392 - t425) * t297 + t524;
t499 = (t212 * t297 + t214 * t299) * t298 + t532;
t331 = qJ(4) * t297 * t300 + pkin(3) * t302;
t257 = -pkin(2) - t331;
t265 = t365 * t300;
t497 = t298 * (t257 * t297 - t265 * t299);
t496 = t298 * (t150 * t299 + t182 * t297);
t292 = pkin(8) * t406;
t218 = t300 * t277 + t292;
t409 = t300 * t301;
t328 = t297 * t303 + t299 * t409;
t492 = t297 ^ 2 + t299 ^ 2;
t61 = t148 * Ifges(5,5) - t147 * Ifges(5,6) - t197 * Ifges(5,3);
t65 = -t197 * Ifges(7,1) + t147 * Ifges(7,4) + t148 * Ifges(7,5);
t66 = -t197 * Ifges(6,1) - t148 * Ifges(6,4) + t147 * Ifges(6,5);
t487 = t66 + t65 + t61;
t486 = -t131 * mrSges(4,1) + t130 * mrSges(4,2) - Ifges(4,5) * t206 - Ifges(4,6) * t207;
t434 = Ifges(4,2) * t263;
t440 = Ifges(4,4) * t264;
t169 = Ifges(4,6) * t538 + t434 + t440;
t253 = Ifges(4,4) * t263;
t442 = Ifges(4,1) * t264;
t170 = Ifges(4,5) * t538 + t253 + t442;
t293 = Ifges(3,4) * t396;
t336 = t184 * t302 + t185 * t300;
t438 = Ifges(4,4) * t302;
t348 = -Ifges(4,2) * t300 + t438;
t439 = Ifges(4,4) * t300;
t350 = Ifges(4,1) * t302 - t439;
t351 = mrSges(4,1) * t300 + mrSges(4,2) * t302;
t446 = t302 / 0.2e1;
t450 = t264 / 0.2e1;
t485 = -t336 * mrSges(4,3) + t281 * t351 + Ifges(3,1) * t397 / 0.2e1 + t293 / 0.2e1 + t372 + t263 * t348 / 0.2e1 + t350 * t450 + t169 * t447 + t170 * t446;
t445 = pkin(8) * t300;
t398 = t302 * t269 + t395 * t445;
t103 = -t301 * t373 + t320 * qJD(2) + (-t292 + (qJ(4) * t411 - t277) * t300) * qJD(3) + t398;
t321 = pkin(8) + t332;
t139 = (qJD(3) * t331 - t297 * t386) * t301 + t321 * t393;
t399 = t300 * t269 + t277 * t391;
t94 = (-t389 + (-pkin(8) * qJD(3) - qJD(2) * t428) * t300) * t303 + (-pkin(8) * t394 - t374 + (qJD(2) * t297 - t299 * t391) * qJ(4)) * t301 + t399;
t18 = -t296 * t94 + t298 * (t103 * t299 + t139 * t297);
t167 = -qJ(4) * t328 + t218;
t262 = t302 * t277;
t175 = -t301 * t383 + t262 + (-pkin(3) - t445) * t303;
t215 = t321 * t301;
t85 = -t296 * t167 + t298 * (t175 * t299 + t215 * t297);
t357 = Ifges(5,3) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1;
t358 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1 + Ifges(5,6) / 0.2e1;
t359 = Ifges(6,4) / 0.2e1 - Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t484 = t147 * t358 + t148 * t359 + t197 * t357 - t61 / 0.2e1 - t65 / 0.2e1 - t66 / 0.2e1 + t521;
t479 = pkin(1) * mrSges(3,1);
t478 = pkin(1) * mrSges(3,2);
t471 = -t147 / 0.2e1;
t470 = t147 / 0.2e1;
t469 = -t148 / 0.2e1;
t468 = t148 / 0.2e1;
t467 = -t169 / 0.2e1;
t460 = -t197 / 0.2e1;
t459 = t197 / 0.2e1;
t458 = t206 / 0.2e1;
t457 = t207 / 0.2e1;
t452 = -t263 / 0.2e1;
t451 = -t264 / 0.2e1;
t444 = qJD(3) / 0.2e1;
t441 = Ifges(3,4) * t301;
t435 = Ifges(4,5) * t302;
t432 = Ifges(4,6) * t300;
t427 = qJD(2) * mrSges(3,2);
t418 = t297 * t298;
t106 = -mrSges(5,1) * t197 - mrSges(5,3) * t148;
t110 = mrSges(6,1) * t148 - mrSges(6,2) * t197;
t405 = -t106 + t110;
t108 = mrSges(6,1) * t147 + mrSges(6,3) * t197;
t109 = -mrSges(7,1) * t147 - mrSges(7,2) * t197;
t404 = -t108 + t109;
t266 = t365 * t302;
t241 = t296 * t266;
t400 = pkin(4) * t416 + t241;
t240 = pkin(3) * t419 + qJ(4) * t418;
t390 = qJD(4) * t296;
t19 = t103 * t419 + t139 * t420 + t298 * t94;
t86 = t298 * t167 + t175 * t419 + t215 * t420;
t141 = t257 * t420 - t265 * t419 + t298 * t266;
t379 = -pkin(3) * t298 - pkin(4);
t366 = -qJ(5) * t296 - pkin(3);
t84 = t119 * mrSges(6,1) + t171 * mrSges(6,2);
t83 = -t118 * mrSges(7,1) + t171 * mrSges(7,2);
t54 = -t119 * mrSges(7,2) + t118 * mrSges(7,3);
t81 = t119 * mrSges(7,1) - t171 * mrSges(7,3);
t49 = -t103 * t297 + t299 * t139;
t123 = -t175 * t297 + t299 * t215;
t180 = t299 * t257 + t265 * t297;
t361 = t299 * t375;
t356 = m(4) * t281 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t263 + mrSges(4,2) * t264 + mrSges(3,3) * t397;
t216 = -qJ(5) * t299 - t240;
t352 = mrSges(4,1) * t302 - mrSges(4,2) * t300;
t349 = Ifges(4,1) * t300 + t438;
t347 = Ifges(4,2) * t302 + t439;
t346 = -t432 + t435;
t341 = t130 * t302 - t131 * t300;
t243 = t297 * t409 - t299 * t303;
t70 = -qJ(5) * t243 - t86;
t330 = t435 / 0.2e1 - t432 / 0.2e1;
t188 = -t296 * t328 + t298 * t407;
t324 = -qJ(5) * t188 + t123;
t323 = -qJ(5) * t238 + t180;
t127 = qJ(5) * t416 - t141;
t317 = 0.3e1 / 0.2e1 * Ifges(3,4) - t330;
t193 = qJD(2) * t327 + t297 * t375;
t12 = -qJ(5) * t193 - qJD(5) * t243 - t19;
t137 = qJD(2) * t493 - t296 * t361 - t298 * t376;
t309 = -qJ(5) * t137 - qJD(5) * t188 + t49;
t306 = t184 * mrSges(4,1) + Ifges(4,6) * t263 + Ifges(4,5) * t264 + t371 - (t303 * Ifges(3,2) + t441) * qJD(1) / 0.2e1 - t185 * mrSges(4,2) + (t538 / 0.2e1 + t444) * Ifges(4,3);
t290 = Ifges(4,3) * t367;
t287 = qJ(4) * t420;
t279 = mrSges(3,3) * t396 - t427;
t239 = pkin(3) * t413 - t287;
t220 = (-pkin(4) * t298 + t366) * t297;
t219 = t299 * t379 + t287;
t217 = -pkin(8) * t408 + t262;
t211 = mrSges(4,1) * t538 - t264 * mrSges(4,3);
t210 = -mrSges(4,2) * t538 + t263 * mrSges(4,3);
t209 = -pkin(8) * t377 + t246;
t195 = (-t298 * t443 + t366) * t297;
t194 = pkin(5) * t418 - t216;
t187 = t296 * t407 + t298 * t328;
t181 = pkin(5) * t420 + t287 + (-qJ(6) + t379) * t299;
t179 = -mrSges(4,2) * t367 + mrSges(4,3) * t207;
t178 = mrSges(4,1) * t367 - mrSges(4,3) * t206;
t156 = -qJD(3) * t218 + t398;
t155 = (-t301 * t394 - t303 * t392) * pkin(8) + t399;
t144 = -mrSges(4,1) * t207 + mrSges(4,2) * t206;
t140 = -t241 + t497;
t136 = qJD(2) * t494 - t296 * t376 + t298 * t361;
t134 = t206 * Ifges(4,1) + t207 * Ifges(4,4) + Ifges(4,5) * t367;
t133 = t206 * Ifges(4,4) + t207 * Ifges(4,2) + Ifges(4,6) * t367;
t128 = t400 - t497;
t122 = pkin(4) * t237 + t323;
t121 = t198 * t298 - t422 * t492;
t120 = t414 * t492 + t426;
t117 = t119 * mrSges(6,3);
t115 = t119 * mrSges(5,2);
t111 = -pkin(5) * t237 - t127;
t107 = mrSges(7,1) * t148 + mrSges(7,3) * t197;
t105 = mrSges(5,2) * t197 - mrSges(5,3) * t147;
t99 = t237 * t443 + t323;
t98 = t265 * t413 + pkin(5) * t238 + (qJ(6) * t302 - t257 * t298) * t297 + t400;
t90 = -mrSges(6,2) * t147 - mrSges(6,3) * t148;
t89 = mrSges(5,1) * t147 + mrSges(5,2) * t148;
t88 = -mrSges(7,2) * t148 + mrSges(7,3) * t147;
t82 = mrSges(6,1) * t118 - mrSges(6,3) * t171;
t80 = mrSges(5,1) * t171 - mrSges(5,3) * t119;
t79 = -mrSges(5,2) * t171 - mrSges(5,3) * t118;
t72 = pkin(4) * t187 + t324;
t71 = -pkin(4) * t243 - t85;
t56 = -t118 * mrSges(6,2) - t117;
t55 = t118 * mrSges(5,1) + t115;
t52 = -t145 + t496;
t46 = t187 * t443 + t324;
t45 = -pkin(5) * t187 - t70;
t44 = t403 - t496;
t42 = pkin(4) * t172 + t326;
t39 = pkin(5) * t188 - t243 * t443 - t85;
t15 = -pkin(4) * t193 - t18;
t13 = pkin(4) * t136 + t309;
t9 = -pkin(5) * t136 - t12;
t7 = pkin(5) * t137 - qJD(6) * t243 - t193 * t443 - t18;
t6 = qJD(6) * t187 + t136 * t443 + t309;
t20 = [(t346 * t444 + t372 + (t303 * t317 - 0.2e1 * t478) * qJD(1) + t356 * pkin(8) + t485) * t393 + m(4) * (t130 * t218 + t131 * t217 + t185 * t155 + t184 * t156) + t517 * t187 + t518 * t188 + t519 * t243 + (-Ifges(5,2) * t471 + Ifges(6,6) * t469 + t489 * t460 - t514 * t468 + t511 * t470 - t534) * t136 + (Ifges(5,4) * t471 - Ifges(6,2) * t469 + t490 * t460 + t515 * t468 + t510 * t470 - t535) * t137 + (-t290 / 0.2e1 + t486) * t303 + t155 * t210 + t156 * t211 + t217 * t178 + t218 * t179 + t123 * t55 + t19 * t105 + t18 * t106 + t7 * t107 + t12 * t108 + t9 * t109 + t15 * t110 + t6 * t88 + t49 * t89 + t13 * t90 + t39 * t81 + t70 * t82 + t45 * t83 + t71 * t84 + t85 * t80 + t86 * t79 + t72 * t56 + t46 * t54 + m(7) * (t1 * t46 + t14 * t7 + t16 * t6 + t17 * t9 + t2 * t39 + t3 * t45) + m(6) * (t12 * t24 + t13 * t21 + t15 * t22 + t4 * t70 + t5 * t72 + t71 * t8) + m(5) * (t10 * t85 + t11 * t86 + t123 * t26 + t18 * t40 + t19 * t41 + t49 * t87) + (t513 * t470 + t487 / 0.2e1 + t491 * t460 + t512 * t468 + Ifges(5,6) * t471 + Ifges(6,4) * t469 - t521) * t193 + (t350 * t458 + t348 * t457 + t134 * t446 + t133 * t447 + pkin(8) * t144 + (-t130 * t300 - t131 * t302) * mrSges(4,3) + (t302 * t467 + t170 * t447 + t349 * t451 + t281 * t352 + t347 * t452 + (t184 * t300 - t185 * t302) * mrSges(4,3) + t538 * t523) * qJD(3) + (-pkin(8) * t279 + t371 + (-0.2e1 * t479 - t317 * t301 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) + (m(4) * pkin(8) + t351) * pkin(8)) * t303) * qJD(1) + t306) * qJD(2)) * t301; (-t14 * t401 + t17 * t402) * mrSges(7,1) + (t14 * t230 - t16 * t402) * mrSges(7,3) + (t230 * t41 - t401 * t87) * mrSges(5,2) + (t21 * t402 - t22 * t230) * mrSges(6,2) + t495 * t89 + t498 * t105 + (t10 * t140 + t11 * t141 + t180 * t26 + t40 * t499 + t41 * t498 + t495 * t87) * m(5) + t499 * t106 + (t21 * t401 + t230 * t24) * mrSges(6,3) + t517 * t237 + t518 * t238 - t519 * t416 + (-t178 * t300 + t179 * t302) * pkin(9) - m(4) * (t184 * t208 + t185 * t209) + (Ifges(5,4) * t177 - Ifges(5,2) * t176 + Ifges(5,6) * t230 + t227 * t511 + t228 * t510) * t470 + (Ifges(5,4) * t228 - Ifges(5,2) * t227 + t176 * t511 + t177 * t510 + t230 * t513) * t471 + t500 * t110 + t501 * t108 + (t122 * t5 + t127 * t4 + t128 * t8 + t21 * t502 + t22 * t500 + t24 * t501) * m(6) + t502 * t90 + t503 * t109 + t504 * t107 + (-t63 + t505) * (t228 / 0.2e1 - t177 / 0.2e1) + (-t64 + t506) * (t227 / 0.2e1 - t176 / 0.2e1) + t507 * t88 + (t1 * t99 + t111 * t3 + t14 * t504 + t16 * t507 + t17 * t503 + t2 * t98) * m(7) - t487 * t230 / 0.2e1 + (t227 * t489 + t228 * t490) * t460 + (t176 * t489 + t177 * t490 + t230 * t491) * t459 + (Ifges(6,4) * t230 - Ifges(6,2) * t177 + Ifges(6,6) * t176 - t227 * t514 + t228 * t515) * t468 + (-Ifges(6,2) * t228 + Ifges(6,6) * t227 - t176 * t514 + t177 * t515 + t230 * t512) * t469 + (-t230 * t40 - t402 * t87) * mrSges(5,1) + (t40 * t401 + t402 * t41) * mrSges(5,3) + (t16 * t401 - t17 * t230) * mrSges(7,2) + ((-qJD(2) * t523 + t371 + (t479 + t441 / 0.2e1) * qJD(1) + (t279 + t427) * pkin(8) - t306) * t301 + (t372 - t293 / 0.2e1 + (t478 + t330 * t303 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t301) * qJD(1) + ((-m(4) * pkin(2) - mrSges(3,1) - t352) * qJD(2) - t356) * pkin(8) - qJD(3) * t346 - t485) * t303) * qJD(1) + t341 * mrSges(4,3) + t300 * t134 / 0.2e1 - t209 * t210 - t208 * t211 + m(4) * (-qJD(3) * t336 + t341) * pkin(9) + t180 * t55 + t140 * t80 + t141 * t79 - pkin(2) * t144 + t127 * t82 + t128 * t84 + t122 * t56 + t111 * t83 + t98 * t81 + t99 * t54 + t347 * t457 + t349 * t458 + t133 * t446 + (-t22 * t401 - t24 * t402) * mrSges(6,1) + ((-pkin(9) * t211 + Ifges(4,5) * t444 + t281 * mrSges(4,2) + t253 / 0.2e1 + t442 / 0.2e1 + t170 / 0.2e1 - t184 * mrSges(4,3)) * t302 + (t467 - pkin(9) * t210 + qJD(3) * t533 + t281 * mrSges(4,1) - t185 * mrSges(4,3) - t434 / 0.2e1 - t440 / 0.2e1 - t484 * t297) * t300) * qJD(3); t528 * t109 + (-t21 * t42 - t22 * t44 + t216 * t4 + t219 * t8 + t220 * t5 + (-t21 * t385 + t22 * t390) * t297 + t529 * t24) * m(6) + t529 * t108 + t530 * t107 + (t1 * t195 + t14 * t530 + t16 * t531 + t17 * t528 + t181 * t2 + t194 * t3) * m(7) + t531 * t88 - t486 + (-t358 * t118 - t359 * t119 + t357 * t171 + t526) * t299 + (-t90 * t385 + t484 * t264 - pkin(3) * t55 + (t105 * qJD(4) - t517) * t298 + (t405 * qJD(4) + t518) * t296) * t297 + (-Ifges(4,2) * t264 + t170 + t253) * t452 + (-t40 * t52 - t41 * t53 - t87 * t97 + t10 * t239 + t11 * t240 + (-pkin(3) * t26 + t388 * t41 - t390 * t40) * t297) * m(5) + (t184 * t263 + t185 * t264) * mrSges(4,3) + t290 - t538 * (Ifges(4,5) * t263 - Ifges(4,6) * t264) / 0.2e1 - t281 * (mrSges(4,1) * t264 + mrSges(4,2) * t263) + t239 * t80 + t240 * t79 + t219 * t84 + t220 * t56 - t184 * t210 + t185 * t211 + t216 * t82 + t194 * t83 + t195 * t54 + t181 * t81 - t53 * t105 - t52 * t106 - t44 * t110 - t97 * t89 - t42 * t90 + (Ifges(4,1) * t263 - t440) * t451 + t169 * t450 + (-Ifges(5,2) * t470 + Ifges(6,6) * t468 + t489 * t459 - t514 * t469 + t511 * t471 + t534) * t172 + (Ifges(5,4) * t470 - Ifges(6,2) * t468 + t490 * t459 + t515 * t469 + t510 * t471 + t535) * t173; t115 - t117 + (mrSges(5,1) - mrSges(6,2)) * t118 + (-t105 - t404) * t121 + (-t107 - t405) * t120 + t54 + (-t120 * t14 - t121 * t17 + t1) * m(7) + (-t120 * t22 + t121 * t24 + t5) * m(6) + (t120 * t40 - t121 * t41 + t26) * m(5); t404 * t197 + (t88 + t90) * t148 + t81 + t84 + (t148 * t16 + t17 * t197 + t2) * m(7) + (t148 * t21 - t197 * t24 + t8) * m(6); -t197 * t107 - t147 * t88 + 0.2e1 * (t3 / 0.2e1 + t14 * t460 + t16 * t471) * m(7) + t83;];
tauc  = t20(:);
