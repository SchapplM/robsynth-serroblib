% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:36
% EndTime: 2019-03-09 08:42:15
% DurationCPUTime: 27.90s
% Computational Cost: add. (7034->720), mult. (14666->903), div. (0->0), fcn. (9359->10), ass. (0->312)
t471 = mrSges(6,2) - mrSges(7,3);
t438 = -m(7) * qJ(6) + t471;
t472 = mrSges(6,1) + mrSges(7,1);
t442 = -m(7) * pkin(5) - t472;
t255 = sin(qJ(2));
t250 = sin(pkin(9));
t251 = cos(pkin(9));
t254 = sin(qJ(5));
t406 = cos(qJ(5));
t279 = -t250 * t406 - t254 * t251;
t269 = t255 * t279;
t131 = qJD(1) * t269;
t321 = qJD(5) * t406;
t344 = qJD(5) * t254;
t155 = -t250 * t321 - t251 * t344;
t489 = t131 + t155;
t322 = t406 * t251;
t363 = t250 * t255;
t452 = -t254 * t363 + t255 * t322;
t130 = t452 * qJD(1);
t156 = -t250 * t344 + t251 * t321;
t488 = t156 + t130;
t257 = cos(qJ(2));
t350 = qJD(1) * t257;
t172 = -t251 * qJD(2) + t250 * t350;
t349 = qJD(2) * t250;
t274 = t251 * t350 + t349;
t112 = t254 * t172 - t406 * t274;
t230 = pkin(7) * t350;
t181 = pkin(3) * t350 + t230;
t249 = qJD(2) * qJ(3);
t153 = qJD(4) + t249 + t181;
t116 = pkin(4) * t274 + t153;
t262 = -t172 * t406 - t254 * t274;
t33 = -pkin(5) * t112 - qJ(6) * t262 + t116;
t107 = Ifges(7,5) * t262;
t351 = qJD(1) * t255;
t220 = qJD(5) + t351;
t49 = t220 * Ifges(7,6) - Ifges(7,3) * t112 + t107;
t384 = Ifges(6,4) * t262;
t52 = Ifges(6,2) * t112 + t220 * Ifges(6,6) + t384;
t496 = t33 * mrSges(7,1) - t52 / 0.2e1 + t49 / 0.2e1;
t476 = -m(6) - m(7);
t495 = qJD(2) / 0.2e1;
t494 = mrSges(7,2) + mrSges(6,3);
t470 = Ifges(3,1) + Ifges(5,3);
t469 = Ifges(6,1) + Ifges(7,1);
t468 = -Ifges(4,4) + Ifges(3,5);
t467 = Ifges(7,4) + Ifges(6,5);
t466 = Ifges(4,5) - Ifges(3,6);
t465 = Ifges(7,5) - Ifges(6,4);
t464 = Ifges(7,2) + Ifges(6,3);
t463 = Ifges(6,6) - Ifges(7,6);
t301 = t257 * mrSges(4,2) - t255 * mrSges(4,3);
t306 = mrSges(3,1) * t257 - mrSges(3,2) * t255;
t493 = t301 - t306;
t408 = t220 / 0.2e1;
t418 = t262 / 0.2e1;
t420 = -t112 / 0.2e1;
t421 = t112 / 0.2e1;
t474 = mrSges(7,3) * t33;
t492 = Ifges(6,4) * t421 + Ifges(7,5) * t420 + t467 * t408 + t469 * t418 - t474;
t491 = -Ifges(6,2) * t421 + Ifges(7,3) * t420 - t463 * t408 + t465 * t418 + t496;
t245 = pkin(9) + qJ(5);
t234 = sin(t245);
t235 = cos(t245);
t457 = mrSges(5,2) * t251;
t303 = mrSges(5,1) * t250 + t457;
t490 = t442 * t234 - t438 * t235 - t303;
t341 = qJD(1) * qJD(2);
t184 = -t257 * qJDD(1) + t255 * t341;
t132 = -qJDD(2) * t250 + t184 * t251;
t133 = qJDD(2) * t251 + t184 * t250;
t47 = qJD(5) * t112 + t254 * t132 + t133 * t406;
t429 = t47 / 0.2e1;
t48 = qJD(5) * t262 - t132 * t406 + t254 * t133;
t427 = t48 / 0.2e1;
t185 = qJDD(1) * t255 + t257 * t341;
t173 = qJDD(5) + t185;
t412 = t173 / 0.2e1;
t304 = mrSges(5,1) * t251 - mrSges(5,2) * t250;
t487 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3) - t304;
t174 = t250 * t254 - t322;
t238 = t255 * qJ(3);
t313 = -pkin(1) - t238;
t397 = pkin(2) + qJ(4);
t134 = (-t257 * t397 + t313) * qJD(1);
t228 = pkin(7) * t351;
t180 = -pkin(3) * t351 - t228;
t142 = -qJD(2) * t397 + qJD(3) - t180;
t78 = -t134 * t250 + t251 * t142;
t59 = pkin(4) * t351 + pkin(8) * t172 + t78;
t79 = t251 * t134 + t250 * t142;
t62 = -pkin(8) * t274 + t79;
t19 = -t254 * t62 + t406 * t59;
t20 = t254 * t59 + t406 * t62;
t346 = qJD(3) * t255;
t288 = -qJD(4) * t257 - t346;
t370 = qJDD(1) * pkin(1);
t291 = -qJ(3) * t185 - t370;
t70 = qJD(1) * t288 + t184 * t397 + t291;
t171 = t185 * pkin(7);
t314 = qJDD(3) + t171;
t97 = pkin(3) * t185 - qJD(2) * qJD(4) - qJDD(2) * t397 + t314;
t34 = -t250 * t70 + t251 * t97;
t21 = pkin(4) * t185 - pkin(8) * t133 + t34;
t35 = t250 * t97 + t251 * t70;
t24 = pkin(8) * t132 + t35;
t3 = t254 * t21 + t406 * t24 + t59 * t321 - t344 * t62;
t4 = -qJD(5) * t20 + t21 * t406 - t254 * t24;
t485 = t174 * t4 - t19 * t489 - t488 * t20 + t279 * t3;
t1 = qJ(6) * t173 + qJD(6) * t220 + t3;
t16 = -t220 * pkin(5) + qJD(6) - t19;
t17 = t220 * qJ(6) + t20;
t2 = -t173 * pkin(5) + qJDD(6) - t4;
t484 = t1 * t279 + t16 * t489 - t488 * t17 - t174 * t2;
t483 = -t257 * t494 + t493;
t428 = -t48 / 0.2e1;
t170 = t184 * pkin(7);
t135 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t170;
t109 = -pkin(3) * t184 + qJDD(4) - t135;
t68 = -pkin(4) * t132 + t109;
t5 = pkin(5) * t48 - qJ(6) * t47 - qJD(6) * t262 + t68;
t481 = 0.2e1 * Ifges(7,3) * t427 + t5 * mrSges(7,1) + t68 * mrSges(6,1) - t47 * Ifges(6,4) / 0.2e1 - t173 * Ifges(6,6) / 0.2e1 + (t427 - t428) * Ifges(6,2) + (t465 + Ifges(7,5)) * t429 + (-t463 + Ifges(7,6)) * t412;
t381 = Ifges(4,6) * t257;
t294 = -t255 * Ifges(4,2) - t381;
t479 = t16 * mrSges(7,1) + t20 * mrSges(6,2) + Ifges(4,4) * t495 + qJD(1) * t294 / 0.2e1 - t17 * mrSges(7,3) - t19 * mrSges(6,1);
t478 = -m(4) - m(5);
t477 = m(5) + m(3);
t475 = -t185 / 0.2e1;
t462 = t173 * t467 + t465 * t48 + t469 * t47;
t29 = mrSges(6,1) * t173 - mrSges(6,3) * t47;
t30 = -t173 * mrSges(7,1) + t47 * mrSges(7,2);
t461 = t30 - t29;
t31 = -mrSges(6,2) * t173 - mrSges(6,3) * t48;
t32 = -mrSges(7,2) * t48 + mrSges(7,3) * t173;
t460 = t31 + t32;
t108 = Ifges(6,4) * t112;
t383 = Ifges(7,5) * t112;
t459 = t220 * t467 + t262 * t469 + t108 - t383;
t224 = pkin(4) * t251 + pkin(3);
t136 = -t224 * t351 - t228;
t458 = pkin(5) * t488 - qJ(6) * t489 + qJD(6) * t174 + qJD(3) - t136;
t390 = mrSges(6,3) * t112;
t86 = -mrSges(6,2) * t220 + t390;
t392 = mrSges(7,2) * t112;
t87 = mrSges(7,3) * t220 + t392;
t395 = t86 + t87;
t389 = mrSges(6,3) * t262;
t88 = mrSges(6,1) * t220 - t389;
t391 = mrSges(7,2) * t262;
t89 = -mrSges(7,1) * t220 + t391;
t394 = t89 - t88;
t328 = mrSges(4,1) * t350;
t199 = -qJD(2) * mrSges(4,3) - t328;
t456 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t350 - t199;
t329 = mrSges(4,1) * t351;
t455 = -mrSges(3,3) * t351 - t329 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t382 = Ifges(4,6) * t255;
t454 = t255 * (-Ifges(4,2) * t257 + t382) + t257 * (Ifges(4,3) * t255 - t381);
t453 = t255 * t466 + t257 * t468;
t451 = t173 * t464 - t463 * t48 + t467 * t47;
t100 = -mrSges(5,2) * t185 + mrSges(5,3) * t132;
t101 = mrSges(5,1) * t185 - mrSges(5,3) * t133;
t450 = t250 * t100 + t251 * t101;
t449 = -t170 * t257 + t171 * t255;
t146 = -qJDD(2) * pkin(2) + t314;
t448 = -t135 * t257 + t146 * t255;
t447 = -t250 * t35 - t251 * t34;
t256 = sin(qJ(1));
t258 = cos(qJ(1));
t446 = g(1) * t258 + g(2) * t256;
t118 = mrSges(5,1) * t274 - t172 * mrSges(5,2);
t58 = -mrSges(6,1) * t112 + mrSges(6,2) * t262;
t445 = t199 - t118 - t58;
t227 = Ifges(3,4) * t350;
t443 = Ifges(3,5) * qJD(2) - Ifges(5,5) * t172 - Ifges(5,6) * t274 + t112 * t463 + t220 * t464 + t262 * t467 + t351 * t470 + t227;
t242 = t257 * pkin(2);
t353 = t242 + t238;
t307 = qJ(4) * t257 + t353;
t169 = -pkin(1) - t307;
t422 = pkin(3) + pkin(7);
t204 = t422 * t255;
t177 = t251 * t204;
t90 = pkin(4) * t255 + t177 + (pkin(8) * t257 - t169) * t250;
t115 = t251 * t169 + t250 * t204;
t359 = t251 * t257;
t96 = -pkin(8) * t359 + t115;
t393 = t254 * t90 + t406 * t96;
t287 = pkin(4) * t257 - pkin(8) * t363;
t348 = qJD(2) * t255;
t232 = pkin(2) * t348;
t290 = -qJ(3) * t257 + qJ(4) * t255;
t121 = qJD(2) * t290 + t232 + t288;
t347 = qJD(2) * t257;
t183 = t422 * t347;
t81 = -t121 * t250 + t251 * t183;
t63 = qJD(2) * t287 + t81;
t360 = t251 * t255;
t335 = pkin(8) * t360;
t82 = t251 * t121 + t250 * t183;
t71 = qJD(2) * t335 + t82;
t9 = -qJD(5) * t393 - t254 * t71 + t406 * t63;
t308 = -m(5) * t397 - mrSges(5,3);
t373 = t257 * mrSges(4,3);
t441 = -t373 + t490 * t257 + (-mrSges(4,2) - t308 + (m(4) - t476) * pkin(2) + t494) * t255;
t293 = -t257 * Ifges(4,3) - t382;
t434 = t251 * (-Ifges(5,4) * t172 - Ifges(5,2) * t274 + Ifges(5,6) * t351) + t250 * (-Ifges(5,1) * t172 - Ifges(5,4) * t274 + Ifges(5,5) * t351) + Ifges(4,5) * qJD(2) + qJD(1) * t293;
t286 = t313 - t242;
t165 = t286 * qJD(1);
t189 = -qJD(2) * pkin(2) + qJD(3) + t228;
t197 = -t230 - t249;
t386 = Ifges(5,4) * t250;
t297 = Ifges(5,2) * t251 + t386;
t385 = Ifges(5,4) * t251;
t300 = Ifges(5,1) * t250 + t385;
t432 = -m(4) * (t189 * t257 + t197 * t255) * pkin(7) - t165 * (-mrSges(4,2) * t255 - t373) - t79 * (-mrSges(5,2) * t257 + mrSges(5,3) * t360) - t78 * (mrSges(5,1) * t257 - mrSges(5,3) * t363) + t153 * t255 * t304 + (Ifges(5,6) * t257 + t255 * t297) * t274 / 0.2e1 + t172 * (Ifges(5,5) * t257 + t255 * t300) / 0.2e1;
t424 = -t133 * Ifges(5,4) / 0.2e1 - t132 * Ifges(5,2) / 0.2e1 + Ifges(5,6) * t475;
t419 = -t262 / 0.2e1;
t416 = t132 / 0.2e1;
t415 = t133 / 0.2e1;
t410 = t185 / 0.2e1;
t409 = -t220 / 0.2e1;
t405 = pkin(7) * t255;
t402 = g(3) * t257;
t240 = t257 * pkin(7);
t396 = -pkin(8) - t397;
t229 = pkin(2) * t351;
t147 = qJD(1) * t290 + t229;
t98 = -t147 * t250 + t251 * t181;
t74 = qJD(1) * t287 + t98;
t99 = t251 * t147 + t250 * t181;
t84 = qJD(1) * t335 + t99;
t28 = t254 * t74 + t406 * t84;
t388 = Ifges(3,4) * t255;
t387 = Ifges(3,4) * t257;
t128 = -mrSges(5,2) * t351 - mrSges(5,3) * t274;
t369 = t128 * t251;
t362 = t250 * t257;
t358 = t255 * t256;
t357 = t255 * t258;
t356 = t256 * t235;
t355 = t256 * t257;
t354 = t257 * t258;
t222 = t250 * pkin(4) + qJ(3);
t205 = t257 * pkin(3) + t240;
t352 = t258 * pkin(1) + t256 * pkin(7);
t345 = qJD(4) * t250;
t343 = qJD(5) * t257;
t342 = -m(5) + t476;
t336 = pkin(4) * t362;
t214 = pkin(4) * t363;
t330 = m(4) - t342;
t326 = t250 * t357;
t154 = pkin(4) * t359 + t205;
t15 = t48 * mrSges(6,1) + t47 * mrSges(6,2);
t14 = t48 * mrSges(7,1) - t47 * mrSges(7,3);
t312 = -t341 / 0.2e1;
t149 = t185 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t83 = -t132 * mrSges(5,1) + t133 * mrSges(5,2);
t310 = pkin(2) * t354 + qJ(3) * t357 + t352;
t305 = mrSges(3,1) * t255 + mrSges(3,2) * t257;
t299 = t257 * Ifges(3,2) + t388;
t295 = Ifges(5,5) * t250 + Ifges(5,6) * t251;
t27 = -t254 * t84 + t406 * t74;
t38 = -t254 * t96 + t406 * t90;
t281 = pkin(1) * t305;
t8 = t254 * t63 + t90 * t321 - t344 * t96 + t406 * t71;
t186 = t396 * t250;
t187 = t396 * t251;
t120 = t186 * t406 + t254 * t187;
t280 = -t186 * t254 + t187 * t406;
t277 = t255 * (Ifges(3,1) * t257 - t388);
t137 = (-pkin(7) - t224) * t348;
t265 = t255 * (Ifges(5,3) * t257 + t255 * t295);
t252 = -pkin(8) - qJ(4);
t243 = t258 * pkin(7);
t217 = qJ(3) * t354;
t216 = qJ(3) * t355;
t190 = -pkin(1) - t353;
t182 = t422 * t348;
t179 = -qJ(3) * t350 + t229;
t178 = t301 * qJD(1);
t161 = Ifges(3,6) * qJD(2) + qJD(1) * t299;
t150 = -qJ(3) * t347 + t232 - t346;
t148 = mrSges(4,1) * t184 - qJDD(2) * mrSges(4,3);
t144 = t279 * t257;
t143 = t254 * t362 - t257 * t322;
t141 = -t234 * t358 + t235 * t258;
t140 = t234 * t258 + t255 * t356;
t139 = t234 * t357 + t356;
t138 = t234 * t256 - t235 * t357;
t129 = mrSges(5,1) * t351 + mrSges(5,3) * t172;
t114 = -t169 * t250 + t177;
t106 = -pkin(5) * t279 + qJ(6) * t174 + t222;
t102 = pkin(2) * t184 - qJD(1) * t346 + t291;
t95 = qJD(2) * t452 - t279 * t343;
t94 = -qJD(2) * t269 + t174 * t343;
t77 = qJD(4) * t322 + qJD(5) * t120 - t254 * t345;
t76 = qJD(4) * t279 + qJD(5) * t280;
t69 = -pkin(5) * t143 - qJ(6) * t144 + t154;
t65 = t133 * Ifges(5,1) + t132 * Ifges(5,4) + t185 * Ifges(5,5);
t57 = -mrSges(7,1) * t112 - mrSges(7,3) * t262;
t56 = pkin(5) * t262 - qJ(6) * t112;
t37 = -t255 * pkin(5) - t38;
t36 = qJ(6) * t255 + t393;
t26 = -pkin(5) * t350 - t27;
t25 = qJ(6) * t350 + t28;
t23 = -pkin(5) * t95 - qJ(6) * t94 - qJD(6) * t144 + t137;
t7 = -pkin(5) * t347 - t9;
t6 = qJ(6) * t347 + qJD(6) * t255 + t8;
t10 = [(t1 * t255 - t144 * t5) * mrSges(7,3) + t294 * t475 + t448 * mrSges(4,1) + (-t116 * mrSges(6,1) + mrSges(7,2) * t17 + mrSges(6,3) * t20 - t491) * t95 + (t116 * mrSges(6,2) - t19 * mrSges(6,3) + t16 * mrSges(7,2) + t459 / 0.2e1 + t492) * t94 + m(5) * (t109 * t205 + t114 * t34 + t115 * t35 - t153 * t182 + t78 * t81 + t79 * t82) + (t144 * t68 - t255 * t3) * mrSges(6,2) + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t481) * t143 + (-t455 * pkin(7) + t189 * mrSges(4,1) + t443 / 0.2e1 + Ifges(6,6) * t421 + Ifges(7,6) * t420 + t464 * t408 + t467 * t418 - t479) * t347 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t449) + t454 * t312 + m(4) * (pkin(7) * t448 + t102 * t190 + t150 * t165) + (-t184 * t240 + t185 * t405 + t449) * mrSges(3,3) + (t453 * t495 - t432) * qJD(2) + t109 * t304 * t257 + (-t357 * t457 - m(3) * t352 - t326 * mrSges(5,1) + t478 * t310 + t476 * (pkin(4) * t326 + t256 * t224 - t252 * t354 + t310) + t442 * t139 + t438 * t138 + (-m(5) * qJ(4) - mrSges(5,3) - t494) * t354 + (-mrSges(2,1) + t493) * t258 + t487 * t256) * g(2) + t306 * t370 - t184 * t299 / 0.2e1 + t102 * t301 + (Ifges(6,4) * t144 + Ifges(6,6) * t255) * t428 + t184 * t293 / 0.2e1 + (-qJDD(2) * mrSges(3,2) - t148) * t240 + (t476 * (t258 * t224 + t252 * t355 + t243) + t442 * t141 + t438 * t140 + (-m(4) - t477) * t243 + t487 * t258 + (-t308 * t257 - (-m(5) * qJ(3) - t303) * t255 - m(4) * t286 + mrSges(2,1) + t476 * (t286 - t214) + t477 * pkin(1) - t483) * t256) * g(1) + (-qJDD(2) * mrSges(3,1) + t149) * t405 + (Ifges(7,5) * t144 + Ifges(7,6) * t255) * t427 - t281 * t341 + t35 * (-mrSges(5,2) * t255 - mrSges(5,3) * t359) - t65 * t362 / 0.2e1 + t34 * (mrSges(5,1) * t255 + mrSges(5,3) * t362) + t393 * t31 + m(6) * (t116 * t137 + t154 * t68 + t19 * t9 + t20 * t8 + t3 * t393 + t38 * t4) + (t257 * (-Ifges(3,2) * t255 + t387) + t277 + t265) * t341 / 0.2e1 + (-t456 * pkin(7) + t434 / 0.2e1 + t197 * mrSges(4,1) - t161 / 0.2e1) * t348 + t257 * (Ifges(3,4) * t185 - Ifges(3,2) * t184 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t257 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t185 + Ifges(4,3) * t184) / 0.2e1 + Ifges(2,3) * qJDD(1) + t2 * (-mrSges(7,1) * t255 + mrSges(7,2) * t144) + t4 * (mrSges(6,1) * t255 - mrSges(6,3) * t144) - t255 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t185 + Ifges(4,6) * t184) / 0.2e1 + t205 * t83 + t190 * (-mrSges(4,2) * t184 - mrSges(4,3) * t185) - pkin(1) * (mrSges(3,1) * t184 + mrSges(3,2) * t185) - t182 * t118 + t150 * t178 + t154 * t15 + t137 * t58 + t82 * t128 + t81 * t129 + t114 * t101 + t115 * t100 + t9 * t88 + t7 * t89 + t8 * t86 + t6 * t87 + t69 * t14 + t23 * t57 + t38 * t29 + t36 * t32 + t37 * t30 + t462 * t144 / 0.2e1 + m(7) * (t1 * t36 + t16 * t7 + t17 * t6 + t2 * t37 + t23 * t33 + t5 * t69) + (t144 * t467 + t255 * t464) * t412 + (t255 * t468 - t257 * t466) * qJDD(2) / 0.2e1 + (t144 * t469 + t255 * t467) * t429 + (-Ifges(3,4) * t184 + Ifges(3,5) * qJDD(2) + Ifges(5,5) * t133 + Ifges(5,6) * t132 + t185 * t470 + t451) * t255 / 0.2e1 + (t255 * t470 - t257 * t295 + t387) * t410 + (Ifges(5,5) * t255 - t257 * t300) * t415 + (Ifges(5,6) * t255 - t257 * t297) * t416 + t359 * t424; (qJ(3) * t109 + t447 * t397 + (-t250 * t79 - t251 * t78) * qJD(4) - t153 * t180 - t78 * t98 - t79 * t99) * m(5) - t434 * t351 / 0.2e1 + (-m(4) * t353 - m(5) * t307 - t257 * mrSges(5,3) + t476 * (-t252 * t257 + t214 + t353) + t490 * t255 + t483) * g(3) + t491 * t156 + t492 * t155 + (Ifges(6,2) * t420 - Ifges(7,3) * t421 + t409 * t463 - t419 * t465 + t496) * t130 - t481 * t279 + (Ifges(6,6) * t420 + Ifges(7,6) * t421 + t464 * t409 + t467 * t419 + t479) * t350 + t453 * t312 + t455 * t230 + t456 * t228 + t446 * t305 + t447 * mrSges(5,3) - (-Ifges(3,2) * t351 + t227 + t443) * t350 / 0.2e1 + t161 * t351 / 0.2e1 + (t155 / 0.2e1 + t131 / 0.2e1) * t459 + (-t116 * t136 + t280 * t4 + t120 * t3 + t222 * t68 + (t76 - t28) * t20 + (-t77 - t27) * t19) * m(6) + (t1 * t120 + t106 * t5 - t280 * t2 + t458 * t33 + (t76 - t25) * t17 + (t77 - t26) * t16) * m(7) - t461 * t280 - t450 * t397 + (-m(4) * t197 + m(5) * t153 + m(6) * t116 - t445) * qJD(3) + t109 * t303 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + t484 * mrSges(7,2) + t485 * mrSges(6,3) + (mrSges(6,1) * t488 + mrSges(6,2) * t489) * t116 + (-pkin(2) * t146 - qJ(3) * t135 - t165 * t179) * m(4) + (-t462 / 0.2e1 - t68 * mrSges(6,2) + t5 * mrSges(7,3) - Ifges(6,4) * t428 - Ifges(7,5) * t427 - t467 * t412 - t469 * t429) * t174 + (-t148 + t83) * qJ(3) - t189 * t328 - t197 * t329 + (-qJD(4) * t251 - t98) * t129 + t251 * t65 / 0.2e1 + t222 * t15 - t179 * t178 - t180 * t118 + t170 * mrSges(3,2) - t171 * mrSges(3,1) - pkin(2) * t149 + t146 * mrSges(4,2) - t135 * mrSges(4,3) - t136 * t58 + t106 * t14 - t27 * t88 - t26 * t89 - t28 * t86 - t25 * t87 + (-t345 - t99) * t128 + t394 * t77 + t395 * t76 + t458 * t57 + t460 * t120 + t466 * t184 + t468 * t185 - (Ifges(6,4) * t420 + Ifges(7,5) * t421 + t467 * t409 + t469 * t419 + t474) * t131 + (t476 * (t252 * t357 + t258 * t336 + t217) + t478 * t217 + t441 * t258) * g(1) + (t476 * (t252 * t358 + t256 * t336 + t216) + t478 * t216 + t441 * t256) * g(2) + (Ifges(5,5) * t251 - t250 * Ifges(5,6)) * t410 + (Ifges(5,1) * t251 - t386) * t415 + (-Ifges(5,2) * t250 + t385) * t416 + t250 * t424 + (t432 + (t454 / 0.2e1 + t281 - t277 / 0.2e1 - t265 / 0.2e1) * qJD(1)) * qJD(1); -t460 * t279 + t461 * t174 + (-t57 + t445) * qJD(2) + t330 * t402 + ((-t129 * t250 + t178 + t369) * qJD(1) - t446 * t330) * t255 + t149 + t450 + t488 * t395 - t489 * t394 + (-qJD(2) * t33 - t484) * m(7) + (-qJD(2) * t116 - t485) * m(6) + (-qJD(2) * t153 - (t250 * t78 - t251 * t79) * t351 - t447) * m(5) + (qJD(2) * t197 + t165 * t351 + t146) * m(4); t369 * t350 - t394 * t262 - t395 * t112 + t128 * t349 - t172 * t129 + t14 + t15 + t83 + (t255 * g(3) + t257 * t446) * t342 + (-t112 * t17 - t16 * t262 + t5) * m(7) + (-t112 * t20 + t19 * t262 + t68) * m(6) + (-t78 * t172 + t274 * t79 + t109) * m(5); -t116 * (mrSges(6,1) * t262 + mrSges(6,2) * t112) - t33 * (mrSges(7,1) * t262 - mrSges(7,3) * t112) + (t112 * t469 + t107 - t384 + t49) * t419 + (Ifges(7,3) * t262 + t383) * t421 + (-m(7) * t17 + t390 - t395) * t19 + (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t17 - (-pkin(5) * t235 - qJ(6) * t234) * t402 - t33 * t56 - g(2) * (pkin(5) * t140 - qJ(6) * t141) - g(1) * (-pkin(5) * t138 + qJ(6) * t139)) * m(7) + (-t234 * t471 + t235 * t472) * t402 + t451 + (-m(7) * t16 + t389 - t394) * t20 + (t138 * t472 + t139 * t471) * g(1) + (-Ifges(6,2) * t262 + t108 + t459) * t420 + (t112 * t467 - t262 * t463) * t409 + qJD(6) * t87 - t56 * t57 + qJ(6) * t32 - pkin(5) * t30 - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t4 * mrSges(6,1) + t1 * mrSges(7,3) + t17 * t391 - t16 * t392 + t52 * t418 + (-t472 * t140 - t471 * t141) * g(2); t262 * t57 - t220 * t87 + (-g(1) * t138 + g(2) * t140 - t17 * t220 - t235 * t402 + t262 * t33 + t2) * m(7) + t30;];
tau  = t10;
