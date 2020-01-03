% Calculate vector of inverse dynamics joint torques for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:50
% EndTime: 2019-12-31 21:53:21
% DurationCPUTime: 15.33s
% Computational Cost: add. (6978->565), mult. (15426->724), div. (0->0), fcn. (10456->10), ass. (0->277)
t262 = cos(qJ(4));
t512 = mrSges(5,1) + mrSges(6,1);
t513 = t262 * t512;
t511 = Ifges(5,4) + Ifges(6,4);
t260 = sin(qJ(2));
t264 = cos(qJ(2));
t345 = qJD(1) * qJD(2);
t212 = qJDD(1) * t264 - t260 * t345;
t213 = qJDD(1) * t260 + t264 * t345;
t259 = sin(qJ(3));
t263 = cos(qJ(3));
t204 = t259 * t260 - t263 * t264;
t274 = t204 * qJD(3);
t115 = -qJD(1) * t274 + t212 * t259 + t213 * t263;
t205 = t259 * t264 + t260 * t263;
t195 = t205 * qJD(1);
t255 = qJD(2) + qJD(3);
t258 = sin(qJ(4));
t172 = -t195 * t258 + t255 * t262;
t254 = qJDD(2) + qJDD(3);
t63 = qJD(4) * t172 + t115 * t262 + t254 * t258;
t429 = t63 / 0.2e1;
t173 = t195 * t262 + t255 * t258;
t64 = -qJD(4) * t173 - t115 * t258 + t254 * t262;
t428 = t64 / 0.2e1;
t275 = t205 * qJD(3);
t116 = -qJD(1) * t275 + t212 * t263 - t213 * t259;
t114 = qJDD(4) - t116;
t427 = t114 / 0.2e1;
t508 = -mrSges(5,3) - mrSges(6,3);
t475 = Ifges(5,1) + Ifges(6,1);
t473 = Ifges(5,5) + Ifges(6,5);
t499 = Ifges(5,2) + Ifges(6,2);
t471 = Ifges(6,6) + Ifges(5,6);
t470 = Ifges(6,3) + Ifges(5,3);
t347 = qJD(4) * t258;
t194 = t204 * qJD(1);
t376 = t194 * t258;
t510 = t347 + t376;
t509 = t473 * t427 + t511 * t428 + t475 * t429;
t476 = mrSges(5,2) + mrSges(6,2);
t507 = t511 * t172;
t506 = t511 * t173;
t505 = t511 * t262;
t504 = t511 * t258;
t266 = -pkin(7) - pkin(6);
t227 = t266 * t264;
t208 = qJD(1) * t227;
t196 = t259 * t208;
t226 = t266 * t260;
t207 = qJD(1) * t226;
t199 = qJD(2) * pkin(2) + t207;
t154 = t199 * t263 + t196;
t409 = pkin(2) * t264;
t246 = pkin(1) + t409;
t225 = t246 * qJD(1);
t465 = t255 * Ifges(4,5);
t503 = -t225 * mrSges(4,2) - t154 * mrSges(4,3) + t465 / 0.2e1;
t190 = qJD(4) + t194;
t125 = pkin(3) * t194 - pkin(8) * t195 - t225;
t197 = t263 * t208;
t155 = t199 * t259 - t197;
t138 = pkin(8) * t255 + t155;
t56 = t262 * t125 - t138 * t258;
t30 = -qJ(5) * t173 + t56;
t27 = pkin(4) * t190 + t30;
t57 = t125 * t258 + t138 * t262;
t31 = qJ(5) * t172 + t57;
t464 = t255 * Ifges(4,6);
t502 = t225 * mrSges(4,1) - t56 * mrSges(5,1) - t27 * mrSges(6,1) + t57 * mrSges(5,2) + t31 * mrSges(6,2) + t464 / 0.2e1;
t256 = qJ(2) + qJ(3);
t251 = sin(t256);
t501 = t251 * t476;
t498 = t114 * t471 + t499 * t64 + t511 * t63;
t496 = t172 * t471 + t173 * t473 + t190 * t470;
t469 = t172 * t499 + t190 * t471 + t506;
t468 = t173 * t475 + t190 * t473 + t507;
t348 = qJD(3) * t263;
t335 = pkin(2) * t348;
t151 = pkin(3) * t195 + pkin(8) * t194;
t352 = qJD(1) * t260;
t338 = pkin(2) * t352;
t132 = t151 + t338;
t158 = t207 * t263 + t196;
t76 = t258 * t132 + t262 * t158;
t495 = t262 * t335 - t76;
t382 = t195 * mrSges(4,3);
t494 = -mrSges(4,1) * t255 - mrSges(5,1) * t172 + mrSges(5,2) * t173 + t382;
t493 = t510 * pkin(4);
t492 = -qJ(5) * t376 + t262 * qJD(5);
t491 = -t258 * t471 + t262 * t473;
t490 = -t258 * t499 + t505;
t489 = t262 * t475 - t504;
t488 = t251 * t513;
t157 = t207 * t259 - t197;
t349 = qJD(3) * t259;
t487 = pkin(2) * t349 - t157;
t346 = qJD(4) * t262;
t375 = t194 * t262;
t486 = -t346 - t375;
t379 = qJDD(1) * pkin(1);
t185 = -pkin(2) * t212 - t379;
t32 = -pkin(3) * t116 - pkin(8) * t115 + t185;
t203 = t213 * pkin(6);
t167 = qJDD(2) * pkin(2) - pkin(7) * t213 - t203;
t202 = t212 * pkin(6);
t171 = pkin(7) * t212 + t202;
t41 = t259 * t167 + t263 * t171 + t199 * t348 + t208 * t349;
t38 = pkin(8) * t254 + t41;
t5 = t125 * t346 - t138 * t347 + t258 * t32 + t262 * t38;
t6 = -qJD(4) * t57 - t258 * t38 + t262 * t32;
t484 = -t258 * t6 + t262 * t5;
t137 = -pkin(3) * t255 - t154;
t396 = mrSges(6,2) * t262;
t294 = mrSges(6,1) * t258 + t396;
t295 = mrSges(5,1) * t258 + mrSges(5,2) * t262;
t90 = -pkin(4) * t172 + qJD(5) + t137;
t483 = t137 * t295 + t294 * t90;
t252 = cos(t256);
t482 = t252 * mrSges(4,1) + (-mrSges(4,2) - t508) * t251;
t450 = t252 * pkin(3) + t251 * pkin(8);
t244 = pkin(4) * t262 + pkin(3);
t257 = -qJ(5) - pkin(8);
t452 = t252 * t244 - t251 * t257;
t480 = -m(5) * t450 - m(6) * t452;
t430 = m(6) * pkin(4);
t478 = t212 / 0.2e1;
t413 = t264 / 0.2e1;
t467 = Ifges(4,4) * t194;
t466 = t194 * Ifges(4,2);
t463 = t264 * Ifges(3,2);
t253 = t262 * qJ(5);
t299 = t195 * pkin(4) + t194 * t253;
t412 = pkin(2) * t259;
t243 = pkin(8) + t412;
t358 = -qJ(5) - t243;
t304 = qJD(4) * t358;
t75 = t262 * t132 - t158 * t258;
t461 = -t299 - t75 + (-qJD(5) - t335) * t258 + t262 * t304;
t309 = qJD(4) * t257;
t79 = t262 * t151 - t154 * t258;
t460 = -qJD(5) * t258 + t262 * t309 - t299 - t79;
t459 = t430 + mrSges(6,1);
t458 = t258 * t304 + t492 + t495;
t80 = t258 * t151 + t262 * t154;
t457 = t258 * t309 + t492 - t80;
t454 = t487 + t493;
t153 = pkin(3) * t204 - pkin(8) * t205 - t246;
t175 = t226 * t259 - t227 * t263;
t168 = t262 * t175;
t92 = t258 * t153 + t168;
t453 = t263 * t226 + t227 * t259;
t283 = -t244 * t251 - t252 * t257;
t408 = pkin(3) * t251;
t411 = pkin(2) * t260;
t449 = -m(6) * (t283 - t411) - m(5) * (-t408 - t411) + t488;
t448 = -m(6) * t283 + t488;
t447 = -t155 + t493;
t445 = t114 * t470 + t471 * t64 + t473 * t63;
t351 = qJD(1) * t264;
t405 = pkin(6) * t264;
t406 = pkin(6) * t260;
t444 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t352) * t405 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t351) * t406;
t265 = cos(qJ(1));
t361 = t258 * t265;
t363 = t252 * t265;
t443 = -t361 * t501 + t363 * t508;
t261 = sin(qJ(1));
t362 = t258 * t261;
t365 = t252 * t261;
t442 = -t362 * t501 + t365 * t508;
t441 = t202 * t264 + t203 * t260;
t440 = g(1) * t265 + g(2) * t261;
t438 = m(5) + m(6) + m(4);
t437 = mrSges(5,1) + t459;
t436 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t435 = -t482 + (t258 * t476 - t513) * t252;
t223 = -mrSges(3,1) * t264 + mrSges(3,2) * t260;
t434 = m(3) * pkin(1) + mrSges(2,1) - t223 + t482;
t395 = mrSges(5,3) * t172;
t119 = -mrSges(5,2) * t190 + t395;
t394 = mrSges(5,3) * t173;
t121 = mrSges(5,1) * t190 - t394;
t24 = mrSges(5,1) * t114 - mrSges(5,3) * t63;
t433 = m(5) * ((-t258 * t57 - t262 * t56) * qJD(4) + t484) - t121 * t346 - t119 * t347 - t258 * t24;
t3 = qJ(5) * t64 + qJD(5) * t172 + t5;
t432 = t6 * mrSges(5,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2);
t425 = -t172 / 0.2e1;
t424 = t172 / 0.2e1;
t423 = -t173 / 0.2e1;
t422 = t173 / 0.2e1;
t421 = -t190 / 0.2e1;
t420 = t190 / 0.2e1;
t419 = t194 / 0.2e1;
t417 = t195 / 0.2e1;
t410 = pkin(2) * t263;
t404 = pkin(8) * t262;
t393 = mrSges(6,3) * t172;
t392 = mrSges(6,3) * t173;
t391 = Ifges(3,4) * t260;
t390 = Ifges(3,4) * t264;
t381 = t195 * Ifges(4,4);
t163 = -qJD(2) * t204 - t274;
t378 = t163 * t258;
t377 = t163 * t262;
t372 = t205 * t258;
t371 = t205 * t262;
t369 = t243 * t262;
t360 = t261 * t262;
t359 = t262 * t265;
t350 = qJD(2) * t260;
t164 = qJD(2) * t205 + t275;
t337 = pkin(2) * t350;
t89 = pkin(3) * t164 - pkin(8) * t163 + t337;
t327 = qJD(2) * t266;
t210 = t260 * t327;
t211 = t264 * t327;
t96 = qJD(3) * t453 + t210 * t263 + t211 * t259;
t341 = t153 * t346 + t258 * t89 + t262 * t96;
t326 = t205 * t346;
t21 = -t64 * mrSges(6,1) + t63 * mrSges(6,2);
t312 = -t347 / 0.2e1;
t310 = -t258 * t96 + t262 * t89;
t308 = t345 / 0.2e1;
t91 = t262 * t153 - t175 * t258;
t298 = mrSges(3,1) * t260 + mrSges(3,2) * t264;
t296 = mrSges(4,1) * t251 + mrSges(4,2) * t252;
t291 = t391 + t463;
t288 = Ifges(3,5) * t264 - Ifges(3,6) * t260;
t282 = -qJ(5) * t163 - qJD(5) * t205;
t42 = t167 * t263 - t259 * t171 - t199 * t349 + t208 * t348;
t279 = pkin(1) * t298;
t183 = -t252 * t361 + t360;
t181 = t252 * t362 + t359;
t278 = t326 + t378;
t277 = t205 * t347 - t377;
t276 = t260 * (Ifges(3,1) * t264 - t391);
t39 = -pkin(3) * t254 - t42;
t97 = qJD(3) * t175 + t210 * t259 - t263 * t211;
t1 = pkin(4) * t114 - qJ(5) * t63 - qJD(5) * t173 + t6;
t133 = t381 + t464 - t466;
t134 = t195 * Ifges(4,1) + t465 - t467;
t20 = -pkin(4) * t64 + qJDD(5) + t39;
t267 = (t486 * t56 - t510 * t57 + t484) * mrSges(5,3) + (-t1 * t258 + t3 * t262 + t27 * t486 - t31 * t510) * mrSges(6,3) + t483 * qJD(4) + (-t467 + t134) * t419 - (-Ifges(4,1) * t194 - t381 + t496) * t195 / 0.2e1 + t258 * t509 + t498 * t262 / 0.2e1 + (t262 * t499 + t504) * t428 + (-Ifges(4,2) * t419 + t470 * t421 + t473 * t423 + t471 * t425 + t502) * t195 + (t258 * t473 + t262 * t471) * t427 + (t258 * t475 + t505) * t429 + t133 * t417 + t155 * t382 + t39 * (-t262 * mrSges(5,1) + mrSges(5,2) * t258) + t20 * (-t262 * mrSges(6,1) + mrSges(6,2) * t258) + (-t491 * t421 - t489 * t423 - t490 * t425 + t483 + t503) * t194 - t41 * mrSges(4,2) + t42 * mrSges(4,1) + Ifges(4,5) * t115 + Ifges(4,6) * t116 + Ifges(4,3) * t254 + (-t376 / 0.2e1 + t312) * t469 + (t346 / 0.2e1 + t375 / 0.2e1) * t468 + (t172 * t490 + t173 * t489 + t190 * t491) * qJD(4) / 0.2e1;
t248 = Ifges(3,4) * t351;
t245 = -pkin(3) - t410;
t236 = pkin(8) * t363;
t235 = pkin(8) * t365;
t224 = t253 + t404;
t222 = t257 * t258;
t219 = -t244 - t410;
t201 = t253 + t369;
t200 = t358 * t258;
t193 = Ifges(3,1) * t352 + Ifges(3,5) * qJD(2) + t248;
t192 = Ifges(3,6) * qJD(2) + qJD(1) * t291;
t184 = t252 * t359 + t362;
t182 = -t252 * t360 + t361;
t176 = -mrSges(4,2) * t255 - mrSges(4,3) * t194;
t150 = mrSges(4,1) * t194 + mrSges(4,2) * t195;
t126 = pkin(4) * t372 - t453;
t120 = mrSges(6,1) * t190 - t392;
t118 = -mrSges(6,2) * t190 + t393;
t103 = -mrSges(4,2) * t254 + mrSges(4,3) * t116;
t102 = mrSges(4,1) * t254 - mrSges(4,3) * t115;
t100 = -mrSges(6,1) * t172 + mrSges(6,2) * t173;
t69 = -qJ(5) * t372 + t92;
t48 = pkin(4) * t204 - t205 * t253 + t91;
t36 = pkin(4) * t278 + t97;
t26 = -mrSges(5,2) * t114 + mrSges(5,3) * t64;
t25 = -mrSges(6,2) * t114 + mrSges(6,3) * t64;
t23 = mrSges(6,1) * t114 - mrSges(6,3) * t63;
t22 = -mrSges(5,1) * t64 + mrSges(5,2) * t63;
t19 = -qJD(4) * t92 + t310;
t18 = -t175 * t347 + t341;
t8 = -qJ(5) * t326 + (-qJD(4) * t175 + t282) * t258 + t341;
t7 = pkin(4) * t164 + t282 * t262 + (-t168 + (qJ(5) * t205 - t153) * t258) * qJD(4) + t310;
t2 = [(t185 * mrSges(4,2) - t42 * mrSges(4,3) + Ifges(4,1) * t115 + Ifges(4,4) * t116 + Ifges(4,5) * t254 + t20 * t294 + t295 * t39 + t468 * t312 + t427 * t491 + t428 * t490 + t429 * t489) * t205 - (-m(4) * t42 + m(5) * t39 - t102 + t22) * t453 + (t445 / 0.2e1 + mrSges(4,1) * t185 + t1 * mrSges(6,1) - mrSges(4,3) * t41 - Ifges(4,4) * t115 - Ifges(4,2) * t116 - Ifges(4,6) * t254 + t427 * t470 + t428 * t471 + t429 * t473 + t432) * t204 + t371 * t509 + t291 * t478 - t192 * t350 / 0.2e1 + (t264 * t390 + t276) * t308 + (Ifges(4,1) * t417 + t134 / 0.2e1 - t467 / 0.2e1 + t503) * t163 + (-Ifges(4,4) * t417 - t155 * mrSges(4,3) - t133 / 0.2e1 + t466 / 0.2e1 + t471 * t424 + t473 * t422 + t470 * t420 + t496 / 0.2e1 - t502) * t164 - t279 * t345 + (-m(4) * t154 + m(5) * t137 + t494) * t97 + m(6) * (t1 * t48 + t126 * t20 + t27 * t7 + t3 * t69 + t31 * t8 + t36 * t90) + (-t326 / 0.2e1 - t378 / 0.2e1) * t469 - t498 * t372 / 0.2e1 + (t212 * t405 + t213 * t406 + t441) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t441) + (t193 * t413 + t288 * qJD(2) / 0.2e1 - t444) * qJD(2) - t223 * t379 + m(5) * (t18 * t57 + t19 * t56 + t5 * t92 + t6 * t91) + t213 * t390 / 0.2e1 + (-t277 * t511 - t278 * t499) * t424 + (-t277 * t475 - t278 * t511) * t422 + (-t362 * t430 - t438 * (t265 * t246 - t261 * t266) - t512 * t184 - t476 * t183 + t436 * t261 + (-t434 + t480) * t265) * g(2) + (-t512 * t182 - t476 * t181 + (-t258 * t430 + t266 * t438 + t436) * t265 + (-m(5) * (-t246 - t450) - m(6) * (-t246 - t452) + m(4) * t246 + t434) * t261) * g(1) + (Ifges(3,4) * t213 + Ifges(3,2) * t212) * t413 + m(4) * (t155 * t96 + t175 * t41 - t185 * t246 - t225 * t337) + (-mrSges(3,1) * t406 - mrSges(3,2) * t405 + 0.2e1 * Ifges(3,6) * t413) * qJDD(2) + (Ifges(3,1) * t213 + Ifges(3,4) * t478 + Ifges(3,5) * qJDD(2) - t308 * t463) * t260 + (-t1 * t371 + t27 * t277 - t278 * t31 - t3 * t372) * mrSges(6,3) + t150 * t337 + t468 * t377 / 0.2e1 + (-t277 * t473 - t278 * t471) * t420 + (t277 * t56 - t278 * t57 - t371 * t6 - t372 * t5) * mrSges(5,3) + t90 * (mrSges(6,1) * t278 - mrSges(6,2) * t277) + t137 * (mrSges(5,1) * t278 - mrSges(5,2) * t277) + t48 * t23 + t69 * t25 + t91 * t24 + t92 * t26 + t36 * t100 + t8 * t118 + t18 * t119 + t7 * t120 + t19 * t121 + t126 * t21 + Ifges(2,3) * qJDD(1) + t175 * t103 + t96 * t176 - pkin(1) * (-mrSges(3,1) * t212 + mrSges(3,2) * t213) - t246 * (-mrSges(4,1) * t116 + mrSges(4,2) * t115); t495 * t119 + t494 * t487 + (-m(4) * t409 - m(5) * (t450 + t409) - m(6) * (t452 + t409) + t223 + t435) * g(3) + t433 * t243 + t192 * t352 / 0.2e1 - t288 * t345 / 0.2e1 + (t444 + (t279 - t276 / 0.2e1) * qJD(1)) * qJD(1) + (t335 - t158) * t176 - t150 * t338 + t458 * t118 + (t1 * t200 + t20 * t219 + t201 * t3 + t461 * t27 + t458 * t31 + t454 * t90) * m(6) + t461 * t120 + t267 + (-g(2) * t235 - g(1) * t236 + t245 * t39 + (t137 * t259 + (-t258 * t56 + t262 * t57) * t263) * qJD(3) * pkin(2) - t137 * t157 - t56 * t75 - t57 * t76) * m(5) + (-t258 * t335 - t75) * t121 + (t261 * t449 + t442) * g(2) + (t265 * t449 + t443) * g(1) + t454 * t100 + t102 * t410 + t103 * t412 + t26 * t369 + (t154 * t157 - t155 * t158 + t225 * t338 + (t259 * t41 + t263 * t42 + (-t154 * t259 + t155 * t263) * qJD(3)) * pkin(2)) * m(4) - (-Ifges(3,2) * t352 + t193 + t248) * t351 / 0.2e1 + Ifges(3,3) * qJDD(2) + t200 * t23 + t201 * t25 - t202 * mrSges(3,2) - t203 * mrSges(3,1) + (m(4) * t411 + t296 + t298) * t440 + Ifges(3,6) * t212 + Ifges(3,5) * t213 + t219 * t21 + t245 * t22; -pkin(3) * t22 - t80 * t119 - t79 * t121 - t154 * t176 - t244 * t21 + t222 * t23 + t224 * t25 + t26 * t404 + t267 + t440 * t296 - t494 * t155 + t460 * t120 + t457 * t118 + t447 * t100 + (t261 * t448 + t442) * g(2) + (t448 * t265 + t443) * g(1) + (t1 * t222 - t20 * t244 + t224 * t3 + t27 * t460 + t31 * t457 + t447 * t90) * m(6) + (-g(1) * (-t265 * t408 + t236) - g(2) * (-t261 * t408 + t235) - t137 * t155 - t56 * t79 - t57 * t80 - pkin(3) * t39) * m(5) + (t435 + t480) * g(3) + t433 * pkin(8); t432 + (t258 * t459 + t295 + t396) * g(3) * t251 + (-t173 * t499 + t468 + t507) * t425 + (t172 * t473 - t173 * t471) * t421 + (-t183 * t437 + t184 * t476) * g(1) + (t395 - t119) * t56 + t459 * t1 + t27 * t393 + (-m(6) * (-t27 + t30) + t392 + t120) * t31 + (t394 + t121) * t57 + (t181 * t437 - t182 * t476) * g(2) + t469 * t422 + (t172 * t475 - t506) * t423 - t30 * t118 - t90 * (mrSges(6,1) * t173 + mrSges(6,2) * t172) - t137 * (mrSges(5,1) * t173 + mrSges(5,2) * t172) + t445 + ((-m(6) * t90 - t100) * t173 + t23) * pkin(4); -t172 * t118 + t173 * t120 + (g(3) * t252 - t31 * t172 + t27 * t173 - t440 * t251 + t20) * m(6) + t21;];
tau = t2;
