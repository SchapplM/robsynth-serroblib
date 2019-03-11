% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:13
% EndTime: 2019-03-09 04:50:45
% DurationCPUTime: 21.16s
% Computational Cost: add. (4513->666), mult. (8824->806), div. (0->0), fcn. (4937->6), ass. (0->297)
t486 = Ifges(7,4) + Ifges(6,5);
t487 = Ifges(6,4) + Ifges(5,5);
t484 = Ifges(7,2) + Ifges(6,3);
t483 = Ifges(5,6) - Ifges(6,6);
t480 = -Ifges(7,5) + t487;
t481 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t193 = cos(qJ(4));
t489 = t486 * t193;
t190 = sin(qJ(4));
t488 = t486 * t190;
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t321 = qJD(1) * qJD(3);
t146 = qJDD(1) * t194 - t191 * t321;
t323 = t193 * qJD(3);
t338 = qJD(1) * t194;
t138 = t190 * t338 - t323;
t333 = qJD(4) * t138;
t64 = qJDD(3) * t190 + t146 * t193 - t333;
t414 = t64 / 0.2e1;
t337 = qJD(3) * t190;
t139 = t193 * t338 + t337;
t65 = qJD(4) * t139 - t193 * qJDD(3) + t146 * t190;
t412 = t65 / 0.2e1;
t147 = -qJDD(1) * t191 - t194 * t321;
t136 = qJDD(4) - t147;
t408 = t136 / 0.2e1;
t406 = t138 / 0.2e1;
t405 = -t139 / 0.2e1;
t404 = t139 / 0.2e1;
t339 = qJD(1) * t191;
t170 = qJD(4) + t339;
t403 = -t170 / 0.2e1;
t485 = Ifges(6,2) + Ifges(5,3);
t482 = Ifges(6,6) - Ifges(7,6);
t444 = t190 * t484 + t489;
t443 = -t483 * t190 + t193 * t487;
t374 = Ifges(5,4) * t190;
t424 = t193 * t481 - t374 + t488;
t479 = t146 / 0.2e1;
t478 = t147 / 0.2e1;
t477 = (-Ifges(5,4) + t486) * t412 + t481 * t414 + t480 * t408;
t476 = t486 * t139;
t197 = -pkin(1) - pkin(7);
t166 = qJD(1) * t197 + qJD(2);
t148 = t191 * t166;
t475 = qJD(5) * t190 + t148;
t187 = t194 * pkin(8);
t149 = pkin(3) * t191 + qJ(2) - t187;
t114 = t149 * qJD(1);
t121 = qJD(3) * pkin(8) + t148;
t330 = qJD(4) * t193;
t331 = qJD(4) * t190;
t322 = qJD(1) * qJD(2);
t167 = qJDD(1) * qJ(2) + t322;
t68 = -pkin(3) * t147 - pkin(8) * t146 + t167;
t163 = qJDD(1) * t197 + qJDD(2);
t335 = qJD(3) * t194;
t84 = t191 * t163 + t166 * t335;
t79 = qJDD(3) * pkin(8) + t84;
t7 = t114 * t330 - t121 * t331 + t190 * t68 + t193 * t79;
t8 = -t114 * t331 - t121 * t330 - t190 * t79 + t193 * t68;
t259 = -t190 * t8 + t193 * t7;
t53 = t193 * t114 - t190 * t121;
t54 = t190 * t114 + t193 * t121;
t474 = -t53 * t330 - t54 * t331 + t259;
t4 = t136 * qJ(5) + t170 * qJD(5) + t7;
t223 = qJDD(5) - t8;
t5 = -pkin(4) * t136 + t223;
t260 = t190 * t5 + t193 * t4;
t38 = -pkin(4) * t170 + qJD(5) - t53;
t162 = t170 * qJ(5);
t39 = t162 + t54;
t473 = t38 * t330 - t39 * t331 + t260;
t410 = pkin(4) + pkin(5);
t312 = t410 * t190;
t359 = qJ(5) * t193;
t224 = t312 - t359;
t472 = t486 * t138;
t253 = t193 * mrSges(6,1) + t190 * mrSges(6,3);
t255 = mrSges(5,1) * t193 - mrSges(5,2) * t190;
t366 = t190 * mrSges(7,2);
t471 = -t253 - t255 - t366;
t377 = Ifges(4,4) * t191;
t250 = t194 * Ifges(4,1) - t377;
t132 = Ifges(5,4) * t138;
t429 = t139 * t481 + t170 * t480 - t132 + t472;
t461 = t138 * t484 + t170 * t482 + t476;
t470 = Ifges(4,5) * qJD(3) + qJD(1) * t250 + t190 * t461 + t193 * t429;
t376 = Ifges(4,4) * t194;
t243 = -t191 * Ifges(4,2) + t376;
t469 = Ifges(7,5) * t404 + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t243 / 0.2e1 + t487 * t405 + (Ifges(7,6) + t483) * t406 + (Ifges(7,3) + t485) * t403;
t1 = -qJ(6) * t64 - qJD(6) * t139 - t136 * t410 + t223;
t2 = qJ(6) * t65 + qJD(6) * t138 + t4;
t468 = t8 * mrSges(5,1) - t5 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3);
t467 = -m(6) - m(7);
t466 = t136 * t482 + t484 * t65 + t486 * t64;
t465 = mrSges(5,3) + mrSges(6,2);
t27 = mrSges(5,1) * t136 - mrSges(5,3) * t64;
t28 = -t136 * mrSges(6,1) + t64 * mrSges(6,2);
t463 = t28 - t27;
t30 = -mrSges(5,2) * t136 - mrSges(5,3) * t65;
t31 = -mrSges(6,2) * t65 + mrSges(6,3) * t136;
t462 = t31 + t30;
t379 = mrSges(7,3) * t138;
t85 = mrSges(7,2) * t170 + t379;
t383 = mrSges(6,2) * t138;
t90 = mrSges(6,3) * t170 - t383;
t459 = t85 + t90;
t380 = mrSges(5,3) * t139;
t88 = mrSges(5,1) * t170 - t380;
t382 = mrSges(6,2) * t139;
t89 = -mrSges(6,1) * t170 + t382;
t458 = -t88 + t89;
t381 = mrSges(5,3) * t138;
t86 = -mrSges(5,2) * t170 - t381;
t457 = -t90 - t86;
t456 = -t170 * t224 + t475;
t302 = t190 * t339;
t325 = qJD(6) * t193;
t386 = pkin(8) - qJ(6);
t261 = pkin(3) * t194 + pkin(8) * t191;
t143 = t261 * qJD(1);
t347 = t193 * t194;
t72 = t190 * t143 + t166 * t347;
t63 = qJ(5) * t338 + t72;
t455 = qJ(6) * t302 - t331 * t386 - t325 - t63;
t397 = pkin(4) * t190;
t228 = -t359 + t397;
t454 = t170 * t228 - t475;
t159 = t386 * t193;
t353 = t191 * t193;
t355 = t190 * t194;
t71 = t143 * t193 - t166 * t355;
t453 = qJD(4) * t159 - qJD(6) * t190 - (qJ(6) * t353 - t194 * t410) * qJD(1) + t71;
t19 = mrSges(5,1) * t65 + mrSges(5,2) * t64;
t451 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t146 - t19;
t450 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t138 + mrSges(5,2) * t139 + mrSges(4,3) * t338;
t192 = sin(qJ(1));
t354 = t191 * t192;
t449 = t193 * t410;
t349 = t192 * t193;
t195 = cos(qJ(1));
t352 = t191 * t195;
t117 = t190 * t352 + t349;
t344 = t195 * t193;
t350 = t192 * t190;
t118 = t191 * t344 - t350;
t448 = t118 * pkin(4) + qJ(5) * t117;
t447 = qJ(5) * t335 + t191 * qJD(5);
t446 = -t191 * t443 + t194 * t485;
t445 = -t191 * t444 + t194 * t482;
t442 = -t193 * t484 + t488;
t441 = t190 * t487 + t193 * t483;
t438 = t136 * t485 - t483 * t65 + t487 * t64;
t32 = qJ(6) * t139 + t53;
t437 = -t32 + qJD(5);
t336 = qJD(3) * t191;
t83 = t194 * t163 - t166 * t336;
t436 = -t191 * t84 - t194 * t83;
t393 = g(2) * t195;
t394 = g(1) * t192;
t435 = t393 - t394;
t434 = -m(5) + t467;
t257 = mrSges(4,1) * t194 - mrSges(4,2) * t191;
t433 = t194 * (-Ifges(4,1) * t191 - t376) / 0.2e1 + qJ(2) * t257;
t432 = mrSges(3,2) - mrSges(4,3) - mrSges(2,1);
t431 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t428 = g(1) * t354;
t427 = m(7) * pkin(5) + mrSges(7,1);
t360 = qJ(5) * t190;
t229 = pkin(4) * t193 + t360;
t150 = -pkin(3) - t229;
t269 = -m(7) * t410 - mrSges(7,1);
t277 = pkin(3) + t360;
t426 = m(5) * pkin(3) - m(6) * t150 + m(7) * t277 - t193 * t269 - t471;
t425 = -t191 * t424 + t194 * t480;
t373 = Ifges(5,4) * t193;
t423 = t190 * t481 + t373 - t489;
t256 = mrSges(4,1) * t191 + mrSges(4,2) * t194;
t262 = -m(7) * t386 + mrSges(7,3);
t421 = -t194 * t262 + mrSges(2,2) - mrSges(3,3) - t256;
t20 = -t170 * t410 + t437;
t33 = qJ(6) * t138 + t54;
t23 = t162 + t33;
t420 = -t1 * t190 - t193 * t2 - t20 * t330 + t23 * t331;
t418 = mrSges(5,1) + mrSges(6,1) + t427;
t357 = t166 * t194;
t122 = -qJD(3) * pkin(3) - t357;
t215 = qJ(5) * t139 - t122;
t24 = -t138 * t410 + qJD(6) + t215;
t364 = t193 * mrSges(6,3);
t252 = t190 * mrSges(6,1) - t364;
t254 = mrSges(5,1) * t190 + mrSges(5,2) * t193;
t41 = pkin(4) * t138 - t215;
t417 = t122 * t254 - t24 * (t190 * mrSges(7,1) - mrSges(7,2) * t193) + t252 * t41;
t416 = qJD(1) ^ 2;
t413 = -t65 / 0.2e1;
t411 = -m(3) - m(4);
t409 = -t136 / 0.2e1;
t407 = -t138 / 0.2e1;
t402 = t170 / 0.2e1;
t74 = mrSges(6,1) * t138 - mrSges(6,3) * t139;
t75 = -mrSges(7,1) * t138 + mrSges(7,2) * t139;
t385 = t74 - t75;
t378 = mrSges(7,3) * t139;
t375 = Ifges(5,4) * t139;
t361 = qJ(5) * t138;
t358 = qJ(6) * t194;
t356 = t190 * t191;
t351 = t191 * t197;
t348 = t192 * t194;
t346 = t194 * t195;
t92 = t190 * t149 + t193 * t351;
t343 = pkin(3) * t348 + pkin(8) * t354;
t342 = t195 * pkin(1) + t192 * qJ(2);
t334 = qJD(3) * t197;
t329 = qJD(4) * t194;
t328 = qJD(4) * t197;
t326 = qJD(5) * t193;
t324 = qJDD(1) * mrSges(3,2);
t319 = t86 + t459;
t87 = -mrSges(7,1) * t170 - t378;
t318 = t87 + t458;
t317 = pkin(8) * t348;
t137 = qJD(3) * t261 + qJD(2);
t299 = t194 * t334;
t305 = t190 * t137 + t149 * t330 + t193 * t299;
t298 = t191 * t328;
t304 = t149 * t331 + t190 * t299 + t193 * t298;
t81 = t191 * qJ(5) + t92;
t303 = t195 * pkin(7) + t342;
t301 = t190 * t336;
t297 = t193 * t329;
t18 = -t65 * mrSges(7,1) + t64 * mrSges(7,2);
t281 = t330 / 0.2e1;
t280 = -t329 / 0.2e1;
t26 = -t136 * mrSges(7,1) - t64 * mrSges(7,3);
t186 = t195 * qJ(2);
t279 = -pkin(1) * t192 + t186;
t276 = -t321 / 0.2e1;
t274 = (t167 + t322) * qJ(2);
t164 = t190 * t351;
t91 = t149 * t193 - t164;
t270 = pkin(3) * t354 + t303;
t258 = qJDD(3) * pkin(3) + t83;
t242 = -Ifges(5,2) * t190 + t373;
t241 = Ifges(5,2) * t193 + t374;
t236 = -Ifges(4,5) * t191 - Ifges(4,6) * t194;
t231 = Ifges(7,5) * t193 + Ifges(7,6) * t190;
t230 = Ifges(7,5) * t190 - Ifges(7,6) * t193;
t227 = t190 * t39 - t193 * t38;
t226 = t190 * t54 + t193 * t53;
t35 = t137 * t193 - t304;
t221 = t191 * (-Ifges(4,2) * t194 - t377);
t218 = Ifges(7,5) * t64 + Ifges(7,6) * t65 - Ifges(7,3) * t136;
t115 = t191 * t350 - t344;
t116 = t190 * t195 + t191 * t349;
t214 = t116 * pkin(4) + qJ(5) * t115 + t270;
t34 = -t190 * t298 + t305;
t213 = -g(1) * t115 + g(2) * t117 - g(3) * t355;
t174 = pkin(3) * t352;
t212 = -pkin(8) * t346 + t192 * t197 + t174 + t186;
t207 = Ifges(5,6) * t194 - t191 * t242;
t203 = -Ifges(7,3) * t194 - t191 * t231;
t202 = -t190 * t319 + t193 * t318;
t201 = qJ(5) * t64 + qJD(5) * t139 + t258;
t179 = -qJDD(1) * pkin(1) + qJDD(2);
t169 = mrSges(7,2) * t347;
t168 = qJ(5) * t347;
t157 = t386 * t190;
t155 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t339;
t142 = t256 * qJD(1);
t129 = t277 + t449;
t128 = t254 * t194;
t103 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t147;
t94 = -t168 + (-t197 + t397) * t194;
t82 = -pkin(4) * t191 - t91;
t80 = t168 + (t197 - t312) * t194;
t73 = pkin(4) * t139 + t361;
t67 = qJ(6) * t355 + t81;
t66 = -pkin(4) * t338 - t71;
t49 = -t138 * Ifges(5,2) + t170 * Ifges(5,6) + t375;
t43 = t164 + (-t149 - t358) * t193 - t410 * t191;
t42 = -t139 * t410 - t361;
t37 = (qJD(4) * t229 - t326) * t194 + (t197 - t228) * t336;
t29 = mrSges(7,2) * t136 + mrSges(7,3) * t65;
t25 = -pkin(4) * t335 - t35;
t22 = t34 + t447;
t21 = (t326 + (-t360 - t449) * qJD(4)) * t194 + (-t197 + t224) * t336;
t17 = mrSges(6,1) * t65 - mrSges(6,3) * t64;
t16 = qJ(6) * t297 + (qJD(6) * t194 + (-qJ(6) * qJD(3) - t328) * t191) * t190 + t305 + t447;
t15 = (qJ(6) * t336 - t137) * t193 + (qJ(6) * t331 - qJD(3) * t410 - t325) * t194 + t304;
t11 = t64 * Ifges(5,4) - t65 * Ifges(5,2) + t136 * Ifges(5,6);
t6 = pkin(4) * t65 - t201;
t3 = -t410 * t65 + qJDD(6) + t201;
t9 = [m(5) * (t122 * t336 * t197 + t34 * t54 + t35 * t53 + t7 * t92 + t8 * t91) + (qJD(3) * t445 - t329 * t442) * t406 + (qJD(3) * t446 - t329 * t441) * t402 + t436 * mrSges(4,3) + m(4) * (-t197 * t436 + t274) + t433 * t321 + (qJD(3) * t425 - t329 * t423) * t404 + t429 * t190 * t280 + t3 * t169 + t16 * t85 + t34 * t86 + t15 * t87 + t35 * t88 + t25 * t89 + t22 * t90 + t91 * t27 + t92 * t30 + t94 * t17 + t37 * t74 + t21 * t75 + t80 * t18 + t81 * t31 + t82 * t28 + t67 * t29 - pkin(1) * t324 + t43 * t26 + (t193 * t280 + t301 / 0.2e1) * t49 - t470 * t336 / 0.2e1 + (t443 * t408 + t444 * t412 + t461 * t281 + t231 * t409 + t242 * t413 + Ifges(4,1) * t479 + Ifges(4,4) * t478 + t424 * t414 + Ifges(4,5) * qJDD(3) + t6 * t252 + (m(5) * t258 + t451) * t197) * t194 + t243 * t478 + t250 * t479 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-mrSges(5,1) * t122 - mrSges(6,1) * t41 + mrSges(7,1) * t24 + mrSges(6,2) * t39 + mrSges(5,3) * t54 - mrSges(7,3) * t23) * (-t297 + t301) + (-mrSges(5,2) * t122 - mrSges(6,2) * t38 - mrSges(7,2) * t24 + mrSges(5,3) * t53 + mrSges(6,3) * t41 + mrSges(7,3) * t20) * (t190 * t329 + t191 * t323) - (Ifges(4,4) * t146 + Ifges(4,2) * t147 + 0.2e1 * Ifges(4,6) * qJDD(3) + t218) * t191 / 0.2e1 - t258 * t128 + (t256 + 0.2e1 * mrSges(3,3)) * t167 + (qJD(3) * t203 - t230 * t329) * t403 + (qJD(3) * t207 - t241 * t329) * t407 + t103 * t351 + m(7) * (t1 * t43 + t15 * t20 + t16 * t23 + t2 * t67 + t21 * t24 + t3 * t80) + m(6) * (t22 * t39 + t25 * t38 + t37 * t41 + t4 * t81 + t5 * t82 + t6 * t94) + (-t11 / 0.2e1 - t3 * mrSges(7,1) - t7 * mrSges(5,3) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + t466 / 0.2e1) * t355 + t155 * t299 + (t485 * t408 + t482 * t412 + t438 / 0.2e1 + t450 * t334 - Ifges(7,3) * t409 + Ifges(5,6) * t413 + t480 * t414 + t468) * t191 + (t53 * mrSges(5,1) - t38 * mrSges(6,1) - t20 * mrSges(7,1) - t54 * mrSges(5,2) + t23 * mrSges(7,2) + t39 * mrSges(6,3) - t469) * t335 + (t5 * mrSges(6,2) - t8 * mrSges(5,3) - t1 * mrSges(7,3) + t477) * t347 + qJD(3) ^ 2 * t236 / 0.2e1 + qJD(2) * t142 + qJ(2) * (-mrSges(4,1) * t147 + mrSges(4,2) * t146) + t179 * mrSges(3,2) + (-m(4) * t186 - m(6) * (t212 + t448) - m(7) * (t174 + t279 + t448) - m(5) * t212 - m(3) * t279 + t465 * t346 + (-m(4) * t197 + m(7) * pkin(7) - t432) * t192 - t418 * t118 + t431 * t117 + t421 * t195) * g(1) + (-m(3) * t342 - m(7) * t214 - m(5) * (t270 - t317) - m(6) * (t214 - t317) - m(4) * t303 + t465 * t348 + t432 * t195 - t418 * t116 + t431 * t115 + t421 * t192) * g(2) + t221 * t276 + m(3) * (-pkin(1) * t179 + t274); t324 + m(3) * t179 + (qJ(2) * t411 - mrSges(3,3)) * t416 + (-t142 - m(5) * t226 - m(6) * t227 - m(7) * (t190 * t23 - t193 * t20) + t202) * qJD(1) + (-t17 + t18 + (t190 * t318 + t193 * t319 + t155) * qJD(3) + m(4) * t83 + m(6) * (t323 * t39 + t337 * t38 - t6) + m(7) * (t20 * t337 + t23 * t323 + t3) + m(5) * (t323 * t54 - t337 * t53 + t258) + t451) * t194 + (t103 + (t29 + t462) * t193 + (t26 + t463) * t190 + (t385 + t450) * qJD(3) + t202 * qJD(4) + m(4) * t84 + m(6) * (qJD(3) * t41 + t473) + m(7) * (-qJD(3) * t24 - t420) + m(5) * (qJD(3) * t122 + t474)) * t191 + t435 * (-t411 - t434); t441 * t408 + t442 * t412 + t429 * t281 + (t221 / 0.2e1 - t433) * t416 + t423 * t414 + t420 * mrSges(7,3) + (t443 / 0.2e1 - t231 / 0.2e1) * qJD(4) * t170 + (t417 + t470 / 0.2e1) * t339 + (t404 * t424 + t417) * qJD(4) + t129 * t18 - t72 * t86 - t71 * t88 - t66 * t89 - t63 * t90 + t83 * mrSges(4,1) - t84 * mrSges(4,2) - pkin(3) * t19 + (-m(5) * t343 + t467 * (t192 * pkin(4) * t347 + t348 * t360 + t343) + (-(-m(7) * qJ(6) - mrSges(7,3)) * t191 + (-t193 * t427 + t471) * t194) * t192) * g(1) + t190 * t477 + (t444 / 0.2e1 - t242 / 0.2e1) * t333 + t258 * t255 + (pkin(3) * t258 - t122 * t148 - t53 * t71 - t54 * t72) * m(5) + (m(7) * t358 + t256 + t434 * t187 + (mrSges(7,3) - t465) * t194 + t426 * t191) * g(3) + (-t428 + t473) * mrSges(6,2) + (-t428 + t474) * mrSges(5,3) + t230 * t409 + t241 * t413 - t257 * t394 + t469 * t338 + t3 * (t193 * mrSges(7,1) + t366) - t155 * t357 - (t331 + t302) * t49 / 0.2e1 - t6 * t253 + ((-t445 / 0.2e1 + t207 / 0.2e1) * t138 + (-t446 / 0.2e1 + t203 / 0.2e1) * t170 + t425 * t405 - t38 * (-mrSges(6,1) * t194 - mrSges(6,2) * t353) - t20 * (-mrSges(7,1) * t194 + mrSges(7,3) * t353) - t53 * (mrSges(5,1) * t194 + mrSges(5,3) * t353) - t23 * (mrSges(7,2) * t194 - mrSges(7,3) * t356) - t54 * (-mrSges(5,2) * t194 + mrSges(5,3) * t356) - t39 * (mrSges(6,2) * t356 + mrSges(6,3) * t194)) * qJD(1) - t450 * t148 + t453 * t87 + t454 * t74 + t455 * t85 + Ifges(4,5) * t146 + Ifges(4,6) * t147 + t150 * t17 + t157 * t26 + t159 * t29 + t456 * t75 + (t1 * t157 + t129 * t3 + t159 * t2 + t20 * t453 + t23 * t455 + t24 * t456) * m(7) + t461 * t331 / 0.2e1 - t466 * t193 / 0.2e1 + t236 * t276 + t193 * t11 / 0.2e1 + (t150 * t6 - t38 * t66 - t39 * t63 + t454 * t41) * m(6) + (t257 + t426 * t194 + (-t262 + t465) * t191) * t393 + (t462 * t193 + t463 * t190 + t457 * t331 + t458 * t330 + m(5) * (-qJD(4) * t226 + t259) + m(6) * (-qJD(4) * t227 + t260) + (m(5) + m(6)) * t191 * t393) * pkin(8) + Ifges(4,3) * qJDD(3); (-Ifges(5,2) * t139 - t132 + t429) * t406 + (-t138 * t481 - t375 + t461 + t476) * t405 + t468 + (-t138 * t487 - t139 * t483) * t403 + (t139 * t484 - t472) * t407 - t32 * t85 - t33 * t87 - t73 * t74 - t42 * t75 - pkin(4) * t28 + (t29 + t31) * qJ(5) + (-pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t39 - t41 * t73) * m(6) - t218 + (qJ(5) * t2 - t1 * t410 - t20 * t33 + t23 * t437 - t24 * t42) * m(7) - t410 * t26 + (-Ifges(7,5) * t138 + Ifges(7,6) * t139) * t402 + t49 * t404 + t39 * t382 + t38 * t383 - t23 * t378 - t20 * t379 + t438 - t41 * (mrSges(6,1) * t139 + mrSges(6,3) * t138) - t24 * (-mrSges(7,1) * t139 - mrSges(7,2) * t138) - t122 * (mrSges(5,1) * t139 - mrSges(5,2) * t138) + (-m(6) * t39 - t381 + t457) * t53 + (-m(6) * t38 + t380 - t458) * t54 + t459 * qJD(5) + (t467 * (-t115 * pkin(4) + qJ(5) * t116) + t431 * t116 + t418 * t115) * g(1) + (t467 * (t117 * pkin(4) - qJ(5) * t118) - t431 * t118 - t418 * t117) * g(2) + (t128 - (t364 + (-m(6) * pkin(4) - mrSges(6,1)) * t190) * t194 - t269 * t355 - t169 + t467 * t168) * g(3); t385 * t139 - t459 * t170 + t26 + t28 + (-t139 * t24 - t170 * t23 + t1 + t213) * m(7) + (t139 * t41 - t170 * t39 + t213 + t5) * m(6); -t138 * t85 + t139 * t87 + (g(3) * t191 - t23 * t138 + t20 * t139 + t194 * t435 + t3) * m(7) + t18;];
tau  = t9;
