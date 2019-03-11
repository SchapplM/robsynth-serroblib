% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:18
% EndTime: 2019-03-09 08:17:53
% DurationCPUTime: 25.41s
% Computational Cost: add. (5194->703), mult. (10951->901), div. (0->0), fcn. (6682->8), ass. (0->319)
t492 = Ifges(5,1) + Ifges(6,1);
t238 = sin(qJ(6));
t236 = cos(pkin(9));
t239 = sin(qJ(2));
t358 = qJD(1) * t239;
t335 = t236 * t358;
t235 = sin(pkin(9));
t241 = cos(qJ(6));
t375 = t235 * t241;
t119 = t238 * t335 - t358 * t375;
t162 = -t236 * t238 + t375;
t138 = t162 * qJD(6);
t469 = t138 + t119;
t279 = t238 * t235 + t241 * t236;
t268 = t279 * t239;
t120 = qJD(1) * t268;
t139 = t279 * qJD(6);
t468 = -t139 + t120;
t493 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t476 = -Ifges(4,4) + Ifges(3,5);
t491 = -Ifges(5,4) + Ifges(6,5);
t475 = Ifges(6,4) + Ifges(5,5);
t485 = Ifges(4,5) - Ifges(3,6);
t474 = -Ifges(5,6) + Ifges(6,6);
t397 = Ifges(6,5) * t236;
t399 = Ifges(5,4) * t236;
t490 = t235 * t492 - t397 + t399;
t242 = cos(qJ(2));
t357 = qJD(1) * t242;
t214 = pkin(7) * t357;
t167 = pkin(3) * t357 + t214;
t234 = qJD(2) * qJ(3);
t137 = qJD(4) + t234 + t167;
t227 = t242 * pkin(2);
t221 = t239 * qJ(3);
t321 = -pkin(1) - t221;
t277 = t321 - t227;
t152 = t277 * qJD(1);
t158 = qJD(2) * t235 + t236 * t357;
t159 = qJD(2) * t236 - t235 * t357;
t211 = pkin(7) * t358;
t351 = qJD(3) + t211;
t175 = -qJD(2) * pkin(2) + t351;
t181 = -t214 - t234;
t398 = Ifges(6,5) * t235;
t293 = -Ifges(6,3) * t236 + t398;
t400 = Ifges(5,4) * t235;
t297 = Ifges(5,2) * t236 + t400;
t304 = mrSges(6,1) * t236 + mrSges(6,3) * t235;
t306 = mrSges(5,1) * t236 - mrSges(5,2) * t235;
t371 = t236 * t239;
t377 = t235 * t239;
t388 = t242 * mrSges(4,3);
t406 = pkin(2) + qJ(4);
t123 = (-t242 * t406 + t321) * qJD(1);
t212 = pkin(3) * t358;
t125 = -qJD(2) * t406 + t212 + t351;
t52 = -t235 * t123 + t125 * t236;
t309 = qJD(5) - t52;
t46 = -pkin(4) * t358 + t309;
t53 = t236 * t123 + t235 * t125;
t47 = qJ(5) * t358 + t53;
t267 = qJ(5) * t159 - t137;
t57 = pkin(4) * t158 - t267;
t489 = -t158 * (Ifges(6,6) * t242 + t239 * t293) / 0.2e1 - (t239 * t490 + t242 * t475) * t159 / 0.2e1 - m(4) * (t175 * t242 + t181 * t239) * pkin(7) + t158 * (Ifges(5,6) * t242 + t239 * t297) / 0.2e1 - t53 * (-mrSges(5,2) * t242 + mrSges(5,3) * t371) - t52 * (mrSges(5,1) * t242 - mrSges(5,3) * t377) - t47 * (mrSges(6,2) * t371 + mrSges(6,3) * t242) - t46 * (-mrSges(6,1) * t242 + mrSges(6,2) * t377) - t152 * (-mrSges(4,2) * t239 - t388) + t137 * t239 * t306 + t239 * t57 * t304;
t456 = Ifges(3,1) + Ifges(6,2) + Ifges(5,3);
t405 = pkin(8) - t406;
t488 = -m(7) * t405 - t493 + (m(5) + m(6)) * t406;
t305 = mrSges(5,1) * t235 + mrSges(5,2) * t236;
t384 = qJ(5) * t236;
t389 = t235 * mrSges(6,1);
t480 = -mrSges(7,1) * t162 + t279 * mrSges(7,2);
t487 = -t305 - m(7) * (pkin(5) * t235 - t384) + t480 - t389 - (-m(6) * qJ(5) - mrSges(6,3)) * t236;
t396 = Ifges(4,6) * t239;
t291 = -t242 * Ifges(4,3) - t396;
t486 = t236 * (t159 * Ifges(5,4) - t158 * Ifges(5,2) + Ifges(5,6) * t358) + Ifges(4,5) * qJD(2) + qJD(1) * t291 + t235 * (t491 * t158 + t159 * t492 + t475 * t358);
t200 = qJD(6) - t358;
t281 = t158 * t238 + t159 * t241;
t415 = Ifges(7,4) * t281;
t92 = t158 * t241 - t159 * t238;
t29 = Ifges(7,2) * t92 + Ifges(7,6) * t200 + t415;
t435 = t29 / 0.2e1;
t349 = qJD(1) * qJD(2);
t170 = -t242 * qJDD(1) + t239 * t349;
t121 = qJDD(2) * t235 - t236 * t170;
t424 = t121 / 0.2e1;
t484 = Ifges(6,3) * t424;
t302 = t242 * mrSges(4,2) - t239 * mrSges(4,3);
t308 = mrSges(3,1) * t242 - mrSges(3,2) * t239;
t483 = t302 - t308;
t171 = qJDD(1) * t239 + t242 * t349;
t354 = qJD(3) * t239;
t278 = -qJD(4) * t242 - t354;
t383 = qJDD(1) * pkin(1);
t285 = -qJ(3) * t171 - t383;
t45 = qJD(1) * t278 + t170 * t406 + t285;
t157 = t171 * pkin(7);
t322 = qJDD(3) + t157;
t71 = pkin(3) * t171 - qJD(2) * qJD(4) - qJDD(2) * t406 + t322;
t20 = t235 * t71 + t236 * t45;
t219 = t239 * qJD(5);
t12 = t171 * qJ(5) + qJD(1) * t219 + t20;
t11 = pkin(8) * t121 + t12;
t122 = qJDD(2) * t236 + t170 * t235;
t19 = -t235 * t45 + t236 * t71;
t311 = qJDD(5) - t19;
t426 = pkin(4) + pkin(5);
t8 = -pkin(8) * t122 - t171 * t426 + t311;
t339 = t426 * t239;
t31 = -pkin(8) * t159 - qJD(1) * t339 + t309;
t33 = pkin(8) * t158 + t47;
t9 = -t238 * t33 + t241 * t31;
t1 = qJD(6) * t9 + t11 * t241 + t238 * t8;
t10 = t238 * t31 + t241 * t33;
t2 = -qJD(6) * t10 - t11 * t238 + t241 * t8;
t482 = t1 * t162 + t468 * t10 - t2 * t279 - t469 * t9;
t75 = -mrSges(6,2) * t121 + mrSges(6,3) * t171;
t76 = -mrSges(5,2) * t171 - mrSges(5,3) * t121;
t77 = mrSges(5,1) * t171 - mrSges(5,3) * t122;
t78 = -t171 * mrSges(6,1) + t122 * mrSges(6,2);
t481 = (t77 - t78) * t236 + (t75 + t76) * t235;
t26 = qJD(6) * t92 + t121 * t238 + t122 * t241;
t437 = t26 / 0.2e1;
t27 = -qJD(6) * t281 + t121 * t241 - t122 * t238;
t436 = t27 / 0.2e1;
t479 = -m(4) - m(5);
t478 = -m(7) - m(6);
t423 = t122 / 0.2e1;
t160 = qJDD(6) - t171;
t422 = t160 / 0.2e1;
t421 = t171 / 0.2e1;
t477 = mrSges(5,2) - mrSges(6,3);
t172 = t405 * t235;
t173 = t405 * t236;
t103 = t172 * t241 - t173 * t238;
t266 = -pkin(8) * t377 - t242 * t426;
t213 = pkin(2) * t358;
t386 = qJ(3) * t242;
t284 = qJ(4) * t239 - t386;
t131 = qJD(1) * t284 + t213;
t72 = -t235 * t131 + t167 * t236;
t44 = qJD(1) * t266 - t72;
t73 = t236 * t131 + t235 * t167;
t63 = qJ(5) * t357 + t73;
t49 = -pkin(8) * t335 + t63;
t473 = qJD(4) * t279 - qJD(6) * t103 + t238 * t49 - t241 * t44;
t102 = -t172 * t238 - t173 * t241;
t472 = -qJD(4) * t162 + qJD(6) * t102 - t238 * t44 - t241 * t49;
t32 = -mrSges(7,1) * t92 + mrSges(7,2) * t281;
t99 = mrSges(6,1) * t158 - mrSges(6,3) * t159;
t404 = -t32 + t99;
t471 = t491 * t121 + t122 * t492 + t475 * t171;
t387 = qJDD(2) / 0.2e1;
t395 = Ifges(4,6) * t242;
t292 = -t239 * Ifges(4,2) - t395;
t470 = Ifges(4,4) * qJD(2) + Ifges(7,5) * t281 + t92 * Ifges(7,6) + t200 * Ifges(7,3) + qJD(1) * t292;
t100 = mrSges(5,1) * t158 + mrSges(5,2) * t159;
t337 = mrSges(4,1) * t357;
t183 = -qJD(2) * mrSges(4,3) - t337;
t467 = -t183 + t100;
t466 = qJD(2) * mrSges(3,2) - mrSges(3,3) * t357 + t183;
t338 = mrSges(4,1) * t358;
t465 = mrSges(3,3) * t358 + t338 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2);
t464 = t239 * (-Ifges(4,2) * t242 + t396) + t242 * (Ifges(4,3) * t239 - t395);
t463 = t239 * t485 + t242 * t476;
t156 = t170 * pkin(7);
t462 = -t156 * t242 + t157 * t239;
t124 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t156;
t130 = -qJDD(2) * pkin(2) + t322;
t461 = -t124 * t242 + t130 * t239;
t288 = t19 * t236 + t20 * t235;
t14 = -pkin(4) * t171 + t311;
t460 = -t12 * t235 + t14 * t236;
t240 = sin(qJ(1));
t243 = cos(qJ(1));
t459 = g(1) * t243 + g(2) * t240;
t458 = -t122 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t424 - t171 * Ifges(5,6) / 0.2e1 + Ifges(6,5) * t423 + Ifges(6,6) * t421 + t484;
t350 = -m(5) + t478;
t457 = -mrSges(3,3) - mrSges(4,1) + mrSges(2,2);
t210 = Ifges(3,4) * t357;
t455 = Ifges(3,5) * qJD(2) + t158 * t474 + t159 * t475 + t358 * t456 + t210;
t453 = -mrSges(2,1) + t483;
t385 = qJ(5) * t235;
t273 = t236 * t426 + t385;
t448 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1);
t447 = -t388 + t487 * t242 + (m(4) * pkin(2) - mrSges(4,2) + t488) * t239;
t445 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t444 = -t9 * mrSges(7,1) + t10 * mrSges(7,2);
t443 = m(7) * pkin(8) + t493;
t439 = Ifges(7,4) * t437 + Ifges(7,2) * t436 + Ifges(7,6) * t422;
t438 = -t26 * Ifges(7,1) / 0.2e1 - t27 * Ifges(7,4) / 0.2e1 - t160 * Ifges(7,5) / 0.2e1;
t89 = Ifges(7,4) * t92;
t30 = Ifges(7,1) * t281 + Ifges(7,5) * t200 + t89;
t434 = t30 / 0.2e1;
t431 = -t92 / 0.2e1;
t430 = t92 / 0.2e1;
t429 = -t281 / 0.2e1;
t428 = t281 / 0.2e1;
t427 = pkin(3) + pkin(7);
t425 = -t121 / 0.2e1;
t420 = -t200 / 0.2e1;
t419 = t200 / 0.2e1;
t417 = -t239 / 0.2e1;
t414 = pkin(7) * t239;
t225 = t242 * pkin(7);
t402 = Ifges(3,4) * t239;
t401 = Ifges(3,4) * t242;
t356 = qJD(2) * t239;
t216 = pkin(2) * t356;
t105 = qJD(2) * t284 + t216 + t278;
t355 = qJD(2) * t242;
t169 = t427 * t355;
t56 = t236 * t105 + t235 * t169;
t376 = t235 * t240;
t374 = t235 * t242;
t370 = t236 * t242;
t369 = t239 * t243;
t368 = t240 * t236;
t367 = t242 * t243;
t115 = -mrSges(6,2) * t158 + mrSges(6,3) * t358;
t116 = -mrSges(5,2) * t358 - mrSges(5,3) * t158;
t366 = t115 + t116;
t117 = mrSges(5,1) * t358 - mrSges(5,3) * t159;
t118 = -mrSges(6,1) * t358 + mrSges(6,2) * t159;
t365 = t117 - t118;
t361 = t227 + t221;
t336 = t242 * qJ(4) + t361;
t154 = -pkin(1) - t336;
t186 = t427 * t239;
t96 = t236 * t154 + t235 * t186;
t362 = t211 + t212;
t187 = t242 * pkin(3) + t225;
t228 = t243 * pkin(7);
t360 = t243 * pkin(3) + t228;
t359 = t243 * pkin(1) + t240 * pkin(7);
t353 = qJD(4) * t235;
t352 = qJD(4) * t236;
t344 = pkin(4) * t374;
t343 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t160;
t340 = m(4) - t350;
t87 = t239 * qJ(5) + t96;
t334 = qJD(5) * t374;
t329 = t358 / 0.2e1;
t320 = -t349 / 0.2e1;
t133 = t171 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t59 = t121 * mrSges(5,1) + t122 * mrSges(5,2);
t58 = t121 * mrSges(6,1) - t122 * mrSges(6,3);
t55 = -t235 * t105 + t169 * t236;
t318 = -qJ(3) + t384;
t95 = -t235 * t154 + t186 * t236;
t43 = qJ(5) * t355 + t219 + t56;
t316 = pkin(2) * t367 + qJ(3) * t369 + t359;
t7 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t307 = mrSges(3,1) * t239 + mrSges(3,2) * t242;
t126 = t279 * t242;
t127 = -t238 * t370 + t241 * t374;
t303 = mrSges(7,1) * t126 + mrSges(7,2) * t127;
t299 = t242 * Ifges(3,2) + t402;
t296 = Ifges(6,4) * t235 - Ifges(6,6) * t236;
t294 = Ifges(5,5) * t235 + Ifges(5,6) * t236;
t290 = pkin(4) * t236 + t385;
t54 = pkin(8) * t374 - t339 - t95;
t64 = pkin(8) * t370 + t87;
t21 = -t238 * t64 + t241 * t54;
t22 = t238 * t54 + t241 * t64;
t61 = -mrSges(7,2) * t200 + mrSges(7,3) * t92;
t62 = mrSges(7,1) * t200 - mrSges(7,3) * t281;
t287 = -t238 * t62 + t241 * t61;
t143 = t235 * t243 + t239 * t368;
t144 = t236 * t243 - t239 * t376;
t283 = t143 * t241 - t144 * t238;
t282 = t143 * t238 + t144 * t241;
t275 = t240 * pkin(3) + qJ(4) * t367 + t316;
t274 = pkin(1) * t307;
t271 = t239 * (Ifges(3,1) * t242 - t402);
t141 = -t236 * t369 + t376;
t257 = -g(1) * t141 + g(2) * t143 - g(3) * t370;
t251 = t239 * (Ifges(6,2) * t242 + t239 * t296);
t250 = t239 * (Ifges(5,3) * t242 + t239 * t294);
t249 = pkin(3) * t170 - qJDD(4) + t124;
t246 = qJ(5) * t122 + qJD(5) * t159 + t249;
t196 = qJ(3) * t367;
t195 = t240 * t386;
t193 = -t236 * qJD(5) + qJD(3);
t176 = -pkin(1) - t361;
t174 = pkin(4) * t235 - t318;
t168 = t427 * t356;
t165 = -qJ(3) * t357 + t213;
t164 = t302 * qJD(1);
t147 = Ifges(3,6) * qJD(2) + qJD(1) * t299;
t142 = t235 * t369 + t368;
t140 = -t235 * t426 + t318;
t136 = -qJ(3) * t355 + t216 - t354;
t132 = mrSges(4,1) * t170 - qJDD(2) * mrSges(4,3);
t108 = t242 * t290 + t187;
t104 = -t290 * t358 - t362;
t101 = -t242 * t273 - t187;
t88 = -pkin(4) * t239 - t95;
t81 = t159 * Ifges(6,5) + Ifges(6,6) * t358 + t158 * Ifges(6,3);
t80 = pkin(2) * t170 - qJD(1) * t354 + t285;
t79 = t334 + (-t290 - t427) * t356;
t74 = t273 * t358 + t362;
t70 = t141 * t238 + t142 * t241;
t69 = t141 * t241 - t142 * t238;
t68 = qJD(6) * t126 + t162 * t356;
t67 = -qJD(2) * t268 + t138 * t242;
t65 = -pkin(4) * t357 - t72;
t60 = -t334 + (t273 + t427) * t356;
t48 = -pkin(4) * t355 - t55;
t36 = -t158 * t426 + t267;
t35 = -pkin(8) * t236 * t356 + t43;
t34 = qJD(2) * t266 - t55;
t23 = pkin(4) * t121 - t246;
t18 = -mrSges(7,2) * t160 + mrSges(7,3) * t27;
t17 = mrSges(7,1) * t160 - mrSges(7,3) * t26;
t13 = -t121 * t426 + t246;
t4 = -qJD(6) * t22 - t238 * t35 + t241 * t34;
t3 = qJD(6) * t21 + t238 * t34 + t241 * t35;
t5 = [(-Ifges(7,1) * t127 + Ifges(7,4) * t126) * t437 + (Ifges(7,1) * t68 + Ifges(7,4) * t67) * t428 + (-t176 * mrSges(4,3) - pkin(1) * mrSges(3,2) - t292 / 0.2e1 - Ifges(4,2) * t417) * t171 + (-qJDD(2) * mrSges(3,2) - t132) * t225 + (t465 * pkin(7) - Ifges(7,5) * t428 - Ifges(7,6) * t430 - Ifges(7,3) * t419 + t444 + t455 / 0.2e1 + t175 * mrSges(4,1) - t470 / 0.2e1) * t355 - t274 * t349 + t461 * mrSges(4,1) + (-t236 * t81 / 0.2e1 - t147 / 0.2e1 + t466 * pkin(7) + t181 * mrSges(4,1) + t486 / 0.2e1) * t356 + (t23 * t304 - t293 * t424 - t297 * t425 - t249 * t306 - t490 * t423 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t171 + (-Ifges(4,3) / 0.2e1 - Ifges(3,2) / 0.2e1) * t170 - 0.2e1 * t485 * t387) * t242 + (t463 * qJD(2) / 0.2e1 - t489) * qJD(2) + (-m(5) * t360 - t282 * mrSges(7,1) - t283 * mrSges(7,2) + t478 * (t144 * pkin(4) + qJ(5) * t143 + t360) + (-m(3) - m(4)) * t228 + t448 * t144 + t477 * t143 + t457 * t243 + (m(3) * pkin(1) - m(4) * t277 + t242 * t488 + t350 * t321 - t453) * t240) * g(1) + (-mrSges(3,1) * t414 + Ifges(4,4) * t417) * qJDD(2) + (t242 * (-Ifges(3,2) * t239 + t401) + t271 + t251 + t250) * t349 / 0.2e1 + m(7) * (t1 * t22 + t10 * t3 + t101 * t13 + t2 * t21 + t36 * t60 + t4 * t9) + m(6) * (t108 * t23 + t12 * t87 + t14 * t88 + t43 * t47 + t46 * t48 + t57 * t79) + m(5) * (-t137 * t168 - t187 * t249 + t19 * t95 + t20 * t96 + t52 * t55 + t53 * t56) + t308 * t383 + (t401 + (-t294 - t296) * t242 + t456 * t239) * t421 + (-mrSges(6,2) * t12 - mrSges(5,3) * t20 + t458) * t370 + m(4) * (t461 * pkin(7) + t136 * t152 + t176 * t80) + (-t170 * t225 + t171 * t414 + t462) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t462 * pkin(7)) + (Ifges(7,5) * t68 + Ifges(7,6) * t67) * t419 + (-Ifges(7,5) * t127 + Ifges(7,6) * t126) * t422 + Ifges(2,3) * qJDD(1) + t68 * t434 + t67 * t435 + t127 * t438 + t126 * t439 + (-t176 * mrSges(4,2) - pkin(1) * mrSges(3,1) - t299 / 0.2e1 + t291 / 0.2e1 + Ifges(4,6) * t417) * t170 + (-Ifges(7,4) * t127 + Ifges(7,2) * t126) * t436 + (Ifges(7,4) * t68 + Ifges(7,2) * t67) * t430 + t187 * t59 + t80 * t302 - t13 * t303 - t168 * t100 + t136 * t164 + t464 * t320 + (t1 * t126 + t10 * t67 + t127 * t2 - t68 * t9) * mrSges(7,3) + (-Ifges(3,4) * t170 + Ifges(3,5) * qJDD(2) + t121 * t474 + t122 * t475 + t171 * t456) * t239 / 0.2e1 + (t19 * mrSges(5,1) - t14 * mrSges(6,1) - t20 * mrSges(5,2) + t12 * mrSges(6,3) - Ifges(7,5) * t437 + Ifges(5,6) * t425 + Ifges(6,6) * t424 - Ifges(7,6) * t436 - Ifges(7,3) * t422 + t387 * t476 + t423 * t475 + t445) * t239 + t21 * t17 + t22 * t18 + (-m(3) * t359 - m(4) * t316 - m(5) * t275 - t70 * mrSges(7,1) - t69 * mrSges(7,2) + t478 * (t142 * pkin(4) + t141 * qJ(5) + t275) + t448 * t142 + t477 * t141 + t443 * t367 + t453 * t243 + t457 * t240) * g(2) + t60 * t32 + t3 * t61 + t4 * t62 + t133 * t414 + t343 * t417 + t36 * (-mrSges(7,1) * t67 + mrSges(7,2) * t68) + t87 * t75 + t88 * t78 + t95 * t77 + t96 * t76 + t79 * t99 + t101 * t7 + t108 * t58 + t43 * t115 + t56 * t116 + t55 * t117 + t48 * t118 + (-mrSges(6,2) * t14 + mrSges(5,3) * t19 - t471 / 0.2e1) * t374; t23 * t389 + (t398 - t400) * t423 + (-t353 - t63) * t115 + (-m(5) * t288 + m(6) * t460 - t481) * t406 + t468 * t435 + ((-t250 / 0.2e1 - t251 / 0.2e1 + t274 - t271 / 0.2e1 + t464 / 0.2e1) * qJD(1) + t489) * qJD(1) + (-m(4) * t361 - m(5) * t336 + t478 * (pkin(4) * t377 + t336) + t443 * t242 + t487 * t239 + t483) * g(3) + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) - t175 * t337 - t181 * t338 + (-t353 - t73) * t116 + (Ifges(7,5) * t138 - Ifges(7,6) * t139) * t419 + (Ifges(7,1) * t138 - Ifges(7,4) * t139) * t428 + (Ifges(7,4) * t138 - Ifges(7,2) * t139) * t430 + (-mrSges(6,3) * t23 + t329 * t81 + t475 * t421 + t423 * t492) * t236 + (t59 - t132) * qJ(3) - t249 * t305 + (Ifges(7,4) * t279 + Ifges(7,2) * t162) * t436 + (Ifges(7,1) * t279 + Ifges(7,4) * t162) * t437 + (Ifges(7,5) * t279 + Ifges(7,6) * t162) * t422 - t279 * t438 + (-qJ(3) * t249 + (-t235 * t53 - t236 * t52) * qJD(4) - t52 * t72 - t53 * t73 + (t362 + qJD(3)) * t137) * m(5) + t362 * t100 + (-Ifges(7,4) * t119 - Ifges(7,2) * t120) * t431 + t13 * t480 + (-Ifges(7,5) * t429 - Ifges(7,6) * t431 - Ifges(7,3) * t420 - t444) * t357 + t397 * t424 + t399 * t425 - (-Ifges(3,2) * t358 + t210 + t455) * t357 / 0.2e1 + (-Ifges(5,2) * t425 + t421 * t474 + t458 + t484) * t235 + t459 * t307 + t460 * mrSges(6,2) - t288 * mrSges(5,3) + t463 * t320 + (t352 - t65) * t118 + t485 * t170 + (-Ifges(7,1) * t119 - Ifges(7,4) * t120) * t429 + (-Ifges(7,5) * t119 - Ifges(7,6) * t120) * t420 - t486 * t358 / 0.2e1 + t162 * t439 + t174 * t58 - t165 * t164 + t156 * mrSges(3,2) - t157 * mrSges(3,1) + t140 * t7 - t465 * t214 - t466 * t211 + t467 * qJD(3) + t469 * t434 + (-mrSges(7,1) * t468 + mrSges(7,2) * t469) * t36 + t470 * t357 / 0.2e1 + t471 * t236 / 0.2e1 + t404 * t193 + t472 * t61 + t473 * t62 + (t1 * t103 + t102 * t2 + t13 * t140 + t473 * t9 + (-t193 - t74) * t36 + t472 * t10) * m(7) + t476 * t171 + (t478 * (t243 * t344 + t196) + t479 * t196 + t447 * t243) * g(1) + (t478 * (t240 * t344 + t195) + t479 * t195 + t447 * t240) * g(2) + t147 * t329 - t74 * t32 + t102 * t17 + t103 * t18 - t104 * t99 - t124 * mrSges(4,3) + t130 * mrSges(4,2) - pkin(2) * t133 + t482 * mrSges(7,3) + (t174 * t23 + (-t235 * t47 + t236 * t46) * qJD(4) - t46 * t65 - t47 * t63 + (t193 - t104) * t57) * m(6) + (-t352 - t72) * t117 + (-pkin(2) * t130 - qJ(3) * t124 - qJD(3) * t181 - t152 * t165) * m(4); -t279 * t17 + t162 * t18 - t469 * t62 + t468 * t61 + (-t404 - t467) * qJD(2) + t340 * t242 * g(3) + ((-t235 * t365 + t236 * t366 + t164) * qJD(1) - t459 * t340) * t239 + t133 + (qJD(2) * t36 + t482) * m(7) + (-t460 - qJD(2) * t57 - (-t235 * t46 - t236 * t47) * t358) * m(6) + (t288 - qJD(2) * t137 - (t235 * t52 - t236 * t53) * t358) * m(5) + (qJD(2) * t181 + t152 * t358 + t130) * m(4) + t481; t366 * t158 + t365 * t159 + t92 * t61 - t281 * t62 + t58 + t59 - t7 + (t10 * t92 - t281 * t9 - t13) * m(7) + (t158 * t47 - t159 * t46 + t23) * m(6) + (t158 * t53 + t159 * t52 - t249) * m(5) + (t239 * g(3) + t242 * t459) * t350; t241 * t17 + t238 * t18 + t404 * t159 + t287 * qJD(6) + (-t115 - t287) * t358 + t78 + (t1 * t238 - t159 * t36 + t2 * t241 + t257 + t200 * (t10 * t241 - t238 * t9)) * m(7) + (t159 * t57 - t358 * t47 + t14 + t257) * m(6); -t36 * (mrSges(7,1) * t281 + mrSges(7,2) * t92) + (Ifges(7,1) * t92 - t415) * t429 + t29 * t428 + (Ifges(7,5) * t92 - Ifges(7,6) * t281) * t420 - t9 * t61 + t10 * t62 - g(1) * (mrSges(7,1) * t69 - mrSges(7,2) * t70) - g(2) * (-mrSges(7,1) * t283 + mrSges(7,2) * t282) - g(3) * t303 + (t10 * t281 + t9 * t92) * mrSges(7,3) + t343 + (-Ifges(7,2) * t281 + t30 + t89) * t431 - t445;];
tau  = t5;
