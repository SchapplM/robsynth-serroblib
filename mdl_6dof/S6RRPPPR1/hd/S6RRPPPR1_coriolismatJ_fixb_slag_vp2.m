% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:14
% EndTime: 2019-03-09 08:06:24
% DurationCPUTime: 6.64s
% Computational Cost: add. (13780->497), mult. (27508->700), div. (0->0), fcn. (30720->8), ass. (0->255)
t410 = sin(qJ(2));
t346 = t410 * pkin(2);
t461 = m(4) * t346;
t378 = sin(pkin(9));
t379 = cos(pkin(9));
t411 = cos(qJ(2));
t270 = t378 * t410 - t379 * t411;
t271 = -t378 * t411 - t379 * t410;
t460 = -t271 * mrSges(4,1) - t270 * mrSges(4,2);
t296 = cos(pkin(10));
t359 = t296 * qJ(5);
t345 = t410 * pkin(7);
t278 = -qJ(3) * t410 - t345;
t292 = t411 * pkin(7);
t349 = t411 * qJ(3) + t292;
t441 = t378 * t278 + t379 * t349;
t319 = t270 * t359 + t441;
t295 = sin(pkin(10));
t367 = t270 * t295;
t118 = -pkin(4) * t367 + t319;
t297 = sin(qJ(6));
t298 = cos(qJ(6));
t272 = -t297 * t295 - t296 * t298;
t172 = t272 * t270;
t394 = t172 * mrSges(7,2);
t273 = t295 * t298 - t297 * t296;
t169 = t273 * t270;
t397 = t169 * mrSges(7,1);
t320 = t397 / 0.2e1 + t394 / 0.2e1;
t435 = -m(7) / 0.2e1;
t437 = -m(6) / 0.2e1;
t439 = -m(5) / 0.2e1;
t430 = pkin(4) + pkin(5);
t343 = t430 * t270;
t97 = -t295 * t343 + t319;
t459 = t118 * t437 + t435 * t97 + t439 * t441 + t320;
t339 = t379 * pkin(2);
t287 = -t339 - pkin(3);
t313 = -t295 * qJ(5) + t287;
t259 = -t296 * pkin(4) + t313;
t276 = -mrSges(6,1) * t296 - mrSges(6,3) * t295;
t458 = m(6) * t259 + t276;
t450 = t273 * t271;
t427 = t450 / 0.2e1;
t456 = mrSges(5,3) + mrSges(6,2);
t207 = t273 * mrSges(7,1) + t272 * mrSges(7,2);
t455 = t207 * qJD(6);
t454 = 0.2e1 * t271;
t453 = Ifges(6,4) + Ifges(5,5);
t452 = Ifges(5,6) - Ifges(6,6);
t314 = t272 * t271;
t393 = t314 * mrSges(7,1);
t327 = -t295 * mrSges(6,1) + t296 * mrSges(6,3);
t191 = t327 * t271;
t91 = mrSges(7,1) * t450 + mrSges(7,2) * t314;
t451 = t191 - t91;
t419 = -t273 / 0.2e1;
t449 = t314 * t419;
t365 = t271 * t295;
t197 = -t270 * mrSges(5,2) + mrSges(5,3) * t365;
t202 = mrSges(6,2) * t365 + t270 * mrSges(6,3);
t448 = t197 + t202;
t358 = t296 * t271;
t200 = t270 * mrSges(5,1) + mrSges(5,3) * t358;
t201 = -t270 * mrSges(6,1) - mrSges(6,2) * t358;
t447 = t200 - t201;
t196 = mrSges(5,2) * t271 + mrSges(5,3) * t367;
t203 = mrSges(6,2) * t367 - mrSges(6,3) * t271;
t446 = t203 + t196;
t193 = -t271 * pkin(3) + t270 * qJ(4) + t346;
t220 = -t379 * t278 + t349 * t378;
t213 = t295 * t220;
t112 = t193 * t296 + t213;
t113 = t295 * t193 - t220 * t296;
t444 = -t112 * t295 + t113 * t296;
t87 = -t271 * qJ(5) + t113;
t88 = pkin(4) * t271 - t112;
t443 = t295 * t88 + t296 * t87;
t390 = t270 * mrSges(7,2);
t395 = t450 * mrSges(7,3);
t126 = t390 - t395;
t391 = t270 * mrSges(7,1);
t392 = t314 * mrSges(7,3);
t128 = -t391 - t392;
t420 = t272 / 0.2e1;
t318 = t126 * t420 + t128 * t419;
t304 = (t272 * t427 + t449) * mrSges(7,3) + t318;
t13 = t304 + t320;
t442 = t13 * qJD(1);
t252 = t298 * t273;
t364 = t272 * t297;
t434 = m(7) / 0.2e1;
t312 = (t252 - t364) * t434;
t440 = t271 ^ 2;
t294 = t296 ^ 2;
t438 = m(5) / 0.2e1;
t436 = m(6) / 0.2e1;
t433 = m(4) * pkin(2);
t432 = mrSges(6,1) / 0.2e1;
t165 = Ifges(7,4) * t450;
t71 = Ifges(7,1) * t314 - t270 * Ifges(7,5) - t165;
t431 = t71 / 0.2e1;
t429 = t169 / 0.2e1;
t428 = -t169 / 0.2e1;
t426 = -t450 / 0.2e1;
t425 = t314 / 0.2e1;
t338 = t378 * pkin(2);
t279 = t338 + qJ(4);
t408 = -pkin(8) + t279;
t260 = t408 * t295;
t261 = t408 * t296;
t188 = t297 * t260 + t261 * t298;
t424 = -t188 / 0.2e1;
t268 = Ifges(7,4) * t272;
t212 = Ifges(7,1) * t273 + t268;
t423 = t212 / 0.2e1;
t422 = -t271 / 0.2e1;
t421 = t271 / 0.2e1;
t418 = t273 / 0.2e1;
t417 = -t276 / 0.2e1;
t416 = -t295 / 0.2e1;
t415 = t295 / 0.2e1;
t414 = -t296 / 0.2e1;
t413 = t296 / 0.2e1;
t412 = -t298 / 0.2e1;
t409 = m(6) * t295;
t405 = Ifges(5,4) * t295;
t404 = Ifges(5,4) * t296;
t403 = Ifges(7,4) * t314;
t402 = Ifges(7,4) * t273;
t401 = Ifges(6,5) * t295;
t400 = Ifges(6,5) * t296;
t187 = t260 * t298 - t297 * t261;
t208 = -mrSges(7,1) * t272 + mrSges(7,2) * t273;
t228 = t296 * t430 - t313;
t277 = -mrSges(5,1) * t296 + mrSges(5,2) * t295;
t388 = t272 * mrSges(7,3);
t340 = t388 / 0.2e1;
t360 = t294 * t270;
t293 = t295 ^ 2;
t361 = t293 * t270;
t351 = (-t360 - t361) * t279;
t368 = t259 * t271;
t387 = t273 * mrSges(7,3);
t299 = (-t271 * t287 + t351) * t438 + (t351 - t368) * t436 + (-t169 * t187 + t172 * t188 + t228 * t271) * t434 + t208 * t421 + t271 * t417 + t277 * t422 + (-t270 * t378 + t271 * t379) * t433 / 0.2e1 + t387 * t429 + t172 * t340 + t456 * (-t361 / 0.2e1 - t360 / 0.2e1);
t125 = -mrSges(7,2) * t271 - t169 * mrSges(7,3);
t127 = mrSges(7,1) * t271 - t172 * mrSges(7,3);
t366 = t270 * t296;
t198 = -mrSges(5,1) * t271 + mrSges(5,3) * t366;
t389 = t271 * mrSges(6,1);
t199 = -mrSges(6,2) * t366 + t389;
t58 = -t213 + (pkin(8) * t270 - t193) * t296 + t430 * t271;
t65 = -pkin(8) * t367 + t87;
t37 = -t297 * t65 + t298 * t58;
t38 = t297 * t58 + t298 * t65;
t301 = (t112 * t296 + t113 * t295) * t438 + (t295 * t87 - t296 * t88) * t436 + (t272 * t37 + t273 * t38) * t434 + t127 * t420 + t125 * t418 + t198 * t413 + t199 * t414 + t461 / 0.2e1 + t446 * t415;
t6 = t299 - t301 - t460;
t399 = qJD(1) * t6;
t289 = -pkin(2) * t411 - pkin(1);
t192 = t270 * pkin(3) + t271 * qJ(4) + t289;
t214 = t295 * t441;
t110 = t192 * t296 - t214;
t111 = t295 * t192 + t296 * t441;
t119 = (-pkin(4) * t295 + t359) * t271 + t220;
t138 = -Ifges(6,6) * t271 + (-t295 * Ifges(6,3) - t400) * t270;
t139 = -Ifges(5,6) * t271 + (t295 * Ifges(5,2) - t404) * t270;
t140 = -Ifges(6,4) * t271 + (-t296 * Ifges(6,1) - t401) * t270;
t141 = -Ifges(5,5) * t271 + (-t296 * Ifges(5,1) + t405) * t270;
t189 = t327 * t270;
t328 = -t295 * mrSges(5,1) - t296 * mrSges(5,2);
t190 = t328 * t270;
t321 = -Ifges(7,5) * t172 / 0.2e1 + Ifges(7,6) * t429;
t57 = t214 + (pkin(8) * t271 - t192) * t296 - t343;
t84 = t270 * qJ(5) + t111;
t64 = -pkin(8) * t365 + t84;
t35 = -t297 * t64 + t298 * t57;
t36 = t297 * t57 + t298 * t64;
t68 = Ifges(7,4) * t172 - Ifges(7,2) * t169 + Ifges(7,6) * t271;
t69 = -Ifges(7,2) * t450 - t270 * Ifges(7,6) + t403;
t70 = Ifges(7,1) * t172 - Ifges(7,4) * t169 + Ifges(7,5) * t271;
t85 = -pkin(4) * t270 - t110;
t90 = t394 + t397;
t98 = (t295 * t430 - t359) * t271 - t220;
t1 = (-Ifges(3,2) + Ifges(3,1)) * t411 * t410 + (mrSges(4,1) * t346 + (-Ifges(6,2) - Ifges(5,3) + Ifges(4,1) - Ifges(7,3) - Ifges(4,2) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t294 + ((Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t295 + (-Ifges(5,4) + Ifges(6,5)) * t296) * t295) * t271 + t321 + (t295 * t452 - t296 * t453 + Ifges(4,4)) * t270) * t270 + m(5) * (t110 * t112 + t111 * t113 + t220 * t441) + (-t410 ^ 2 + t411 ^ 2) * Ifges(3,4) + t68 * t426 + t69 * t428 + t172 * t431 + t70 * t425 + m(6) * (t118 * t119 + t84 * t87 + t85 * t88) + m(7) * (t35 * t37 + t36 * t38 - t97 * t98) + (t460 + t461) * t289 - pkin(1) * (mrSges(3,1) * t410 + mrSges(3,2) * t411) + (-mrSges(4,2) * t346 + Ifges(7,5) * t425 + Ifges(7,6) * t426 - Ifges(4,4) * t271 + (-t140 / 0.2e1 - t141 / 0.2e1 - t441 * mrSges(5,2) + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t271) * t296 + (-t138 / 0.2e1 + t139 / 0.2e1 - t441 * mrSges(5,1) + (-Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t271) * t295) * t271 - t97 * t91 + t98 * t90 + t36 * t125 + t38 * t126 + t35 * t127 + t37 * t128 + t119 * t189 + t118 * t191 + t111 * t196 + t113 * t197 + t110 * t198 + t85 * t199 + t112 * t200 + t88 * t201 + t87 * t202 + t84 * t203 + t220 * t190;
t398 = t1 * qJD(1);
t384 = t297 * mrSges(7,1);
t383 = t298 * mrSges(7,2);
t352 = -Ifges(7,5) * t450 - Ifges(7,6) * t314;
t89 = -t450 * mrSges(7,2) + t393;
t92 = -Ifges(7,2) * t314 - t165;
t93 = -Ifges(7,1) * t450 - t403;
t4 = t98 * t89 - t270 * t352 / 0.2e1 + t35 * t126 - t36 * t128 + (-t36 * mrSges(7,3) + t93 / 0.2e1 - t69 / 0.2e1) * t314 - (-t35 * mrSges(7,3) + t431 + t92 / 0.2e1) * t450;
t382 = t4 * qJD(1);
t324 = t295 * t110 - t296 * t111;
t326 = t295 * t85 + t296 * t84;
t369 = t220 * t271;
t5 = t440 * t328 + t169 * t128 - t172 * t126 - m(7) * (-t169 * t35 + t172 * t36) + m(6) * t326 * t270 - m(5) * (t270 * t324 - t369) - m(4) * (-t270 * t441 - t369) - t447 * t367 + t448 * t366 + (m(6) * t119 - m(7) * t98 + t451) * t271 + (-t270 ^ 2 - t440) * mrSges(4,3);
t381 = t5 * qJD(1);
t12 = m(7) * (t314 * t35 + t36 * t450) + t450 * t126 + t314 * t128 + (t447 * t296 + t448 * t295 + m(6) * (t295 * t84 - t296 * t85) + m(5) * (t110 * t296 + t111 * t295)) * t271;
t377 = qJD(1) * t12;
t354 = t298 * t126;
t357 = t297 * t128;
t15 = t451 * t358 + (t202 + t354 - t357) * t270 + m(7) * (-t98 * t358 + (-t297 * t35 + t298 * t36) * t270) + m(6) * (t119 * t358 + t270 * t84);
t376 = qJD(1) * t15;
t317 = m(7) * (-t169 * t298 + t297 * t172);
t45 = -t317 / 0.2e1 + (t312 + t409) * t270;
t375 = qJD(1) * t45;
t371 = t314 * t298;
t356 = t297 * t450;
t355 = t297 * t273;
t350 = Ifges(7,5) * t272 - Ifges(7,6) * t273;
t348 = t293 + t294;
t347 = m(6) / 0.4e1 + m(5) / 0.4e1;
t337 = -t367 / 0.2e1;
t336 = -t366 / 0.2e1;
t333 = -t358 / 0.2e1;
t332 = t358 / 0.2e1;
t303 = (-t119 * t295 + (t270 * t279 + t368) * t296) * t437 + (-t228 * t358 + t295 * t98 + (-t187 * t297 + t188 * t298) * t270) * t435;
t305 = t88 * t436 + (t297 * t38 + t298 * t37) * t434 + t297 * t125 / 0.2e1 + t298 * t127 / 0.2e1;
t10 = (-t91 / 0.2e1 + t191 / 0.2e1) * t295 + (t432 + (t208 / 0.2e1 + t417) * t296) * t271 + (-t296 * mrSges(6,2) + (t272 * t412 - t355 / 0.2e1) * mrSges(7,3)) * t270 + t303 + t305;
t77 = (m(7) * t228 + t208 - t458) * t295;
t325 = qJD(1) * t10 - qJD(2) * t77;
t48 = 0.2e1 * mrSges(7,2) * t427 - t393;
t323 = qJD(1) * t48 - qJD(2) * t207;
t148 = t409 + (-t364 / 0.2e1 + t252 / 0.2e1 + t415) * m(7);
t53 = m(6) * t358 + (t332 - t371 / 0.2e1 - t356 / 0.2e1) * m(7);
t322 = qJD(1) * t53 - qJD(2) * t148;
t18 = (-t390 / 0.2e1 - t126 / 0.2e1 - t395 / 0.2e1) * t298 + (-t391 / 0.2e1 + t392 / 0.2e1 + t128 / 0.2e1) * t297;
t316 = t18 * qJD(1);
t307 = (t272 * t314 + t273 * t450) * t434 + t347 * t348 * t454;
t31 = (m(7) / 0.4e1 + t347) * t454 + t307;
t315 = t31 * qJD(1);
t209 = -Ifges(7,2) * t273 + t268;
t210 = Ifges(7,2) * t272 + t402;
t211 = Ifges(7,1) * t272 - t402;
t22 = t228 * t207 + (t211 / 0.2e1 - t210 / 0.2e1) * t273 + (t423 + t209 / 0.2e1) * t272;
t302 = (-t69 / 0.4e1 + t93 / 0.4e1) * t273 + (t92 / 0.4e1 + t71 / 0.4e1) * t272 - (-t187 * mrSges(7,3) / 0.2e1 + t212 / 0.4e1 + t209 / 0.4e1) * t450 + (mrSges(7,3) * t424 + t211 / 0.4e1 - t210 / 0.4e1) * t314 + t187 * t126 / 0.2e1 + t128 * t424 + t228 * t89 / 0.2e1 - t270 * t350 / 0.4e1 + t98 * t207 / 0.2e1;
t306 = Ifges(7,3) * t422 - t37 * mrSges(7,1) / 0.2e1 + t38 * mrSges(7,2) / 0.2e1 + t321;
t3 = t302 + t306;
t310 = -t3 * qJD(1) - t22 * qJD(2);
t39 = (-t272 ^ 2 - t273 ^ 2) * mrSges(7,3) + m(7) * (t187 * t273 - t188 * t272) + (0.4e1 * t279 * t347 + t456) * t348;
t300 = (-t200 / 0.2e1 + t201 / 0.2e1) * t295 + (t202 / 0.2e1 + t197 / 0.2e1) * t296 + (t420 * t450 + t449) * mrSges(7,3) + t324 * t439 + t326 * t436 + (t187 * t314 + t188 * t450 - t272 * t36 + t273 * t35) * t434 - t318;
t8 = ((-mrSges(6,3) / 0.2e1 + mrSges(5,2) / 0.2e1) * t296 + (t432 + mrSges(5,1) / 0.2e1) * t295) * t270 + t300 + t459;
t309 = t8 * qJD(1) + t39 * qJD(2);
t153 = m(7) * t416 + t312;
t52 = (t356 + t371) * t434 + m(6) * t333 + (t434 + t436) * t358;
t44 = m(6) * t337 + t317 / 0.2e1 + (t312 + t409 / 0.2e1) * t270;
t30 = t307 + (m(5) + m(6) + m(7)) * t422;
t19 = t354 / 0.2e1 - t297 * t392 / 0.2e1 - t357 / 0.2e1 - t395 * t412 + (-t383 / 0.2e1 - t384 / 0.2e1) * t270;
t14 = t304 - t320;
t11 = t91 * t415 + t208 * t333 + t191 * t416 + t276 * t332 + t389 / 0.2e1 - t303 + t305 + (t298 * t340 + mrSges(7,3) * t355 / 0.2e1) * t270;
t9 = mrSges(5,2) * t336 + mrSges(6,3) * t366 / 0.2e1 + t300 + (mrSges(5,1) + mrSges(6,1)) * t337 - t459;
t7 = t299 + t301;
t2 = t302 - t306;
t16 = [qJD(2) * t1 - qJD(3) * t5 + qJD(4) * t12 + qJD(5) * t15 + qJD(6) * t4, t7 * qJD(3) + t9 * qJD(4) + t11 * qJD(5) + t2 * qJD(6) + t398 + ((m(5) * t287 - t379 * t433 - mrSges(4,1) + t277) * t441 + (-t400 + t404 + (Ifges(5,1) + Ifges(6,1)) * t295) * t336 + ((-t198 + t199) * t295 + m(5) * t444 + t446 * t296 + m(6) * t443) * t279 + (t295 * t453 + t296 * t452) * t422 + t443 * mrSges(6,2) + t444 * mrSges(5,3) + (t140 + t141) * t415 + t210 * t428 + t172 * t423 + t70 * t418 + t68 * t420 + (Ifges(7,5) * t273 + Ifges(7,6) * t272) * t421 + t139 * t413 + t138 * t414 + t38 * t388 + mrSges(3,2) * t345 + m(7) * (t187 * t37 + t188 * t38 - t228 * t97) + (t270 * t339 + t271 * t338) * mrSges(4,3) - (t378 * t433 - mrSges(4,2)) * t220 - Ifges(3,6) * t410 + Ifges(3,5) * t411 + (-Ifges(6,3) * t296 + t401) * t337 + (Ifges(5,2) * t296 + t405) * t367 / 0.2e1 - t37 * t387 - mrSges(3,1) * t292 + t287 * t190 + Ifges(4,6) * t271 - Ifges(4,5) * t270 + t187 * t127 + t188 * t125 - t97 * t208 + t228 * t90 + t458 * t118 + t259 * t189) * qJD(2), -t381 + t7 * qJD(2) + (-t169 * t272 + t172 * t273) * m(7) * qJD(3) + t30 * qJD(4) + t44 * qJD(5) + t14 * qJD(6), qJD(2) * t9 + qJD(3) * t30 + qJD(5) * t52 + t377, qJD(2) * t11 + qJD(3) * t44 + qJD(4) * t52 + qJD(6) * t19 + t376, t382 + t2 * qJD(2) + t14 * qJD(3) + t19 * qJD(5) + (-mrSges(7,1) * t36 - mrSges(7,2) * t35 + t352) * qJD(6); qJD(3) * t6 + qJD(4) * t8 - qJD(5) * t10 + qJD(6) * t3 - t398, qJD(4) * t39 + qJD(5) * t77 + qJD(6) * t22, t399, t153 * qJD(5) + t309, qJD(4) * t153 - t325 (-mrSges(7,1) * t188 - mrSges(7,2) * t187 + t350) * qJD(6) - t310; -qJD(2) * t6 + qJD(4) * t31 + qJD(5) * t45 + qJD(6) * t13 + t381, -t399, 0, t315, t375, t442 - t455; -qJD(2) * t8 - qJD(3) * t31 + qJD(5) * t53 + qJD(6) * t48 - t377, -t148 * qJD(5) - t309 - t455, -t315, 0, t322, t323; qJD(2) * t10 - qJD(3) * t45 - qJD(4) * t53 - qJD(6) * t18 - t376, qJD(4) * t148 + t325, -t375, -t322, 0 (-t383 - t384) * qJD(6) - t316; -qJD(2) * t3 - qJD(3) * t13 - qJD(4) * t48 + qJD(5) * t18 - t382, qJD(4) * t207 + t310, -t442, -t323, t316, 0;];
Cq  = t16;
