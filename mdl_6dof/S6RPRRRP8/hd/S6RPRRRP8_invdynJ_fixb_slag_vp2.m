% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:57
% EndTime: 2019-03-09 06:23:26
% DurationCPUTime: 16.36s
% Computational Cost: add. (9591->679), mult. (18348->858), div. (0->0), fcn. (11828->10), ass. (0->311)
t506 = mrSges(6,1) + mrSges(7,1);
t505 = mrSges(6,2) - mrSges(7,3);
t504 = mrSges(6,3) + mrSges(7,2);
t503 = Ifges(6,1) + Ifges(7,1);
t535 = Ifges(6,4) - Ifges(7,5);
t502 = Ifges(7,4) + Ifges(6,5);
t501 = Ifges(7,2) + Ifges(6,3);
t500 = Ifges(6,6) - Ifges(7,6);
t278 = sin(qJ(3));
t456 = cos(qJ(4));
t360 = t456 * t278;
t277 = sin(qJ(4));
t457 = cos(qJ(3));
t363 = t277 * t457;
t298 = -t360 - t363;
t198 = t298 * qJD(1);
t194 = qJD(5) - t198;
t336 = t456 * t457;
t392 = qJD(1) * t278;
t199 = qJD(1) * t336 - t277 * t392;
t276 = sin(qJ(5));
t280 = cos(qJ(5));
t381 = qJD(3) + qJD(4);
t168 = t276 * t199 - t280 * t381;
t438 = mrSges(7,2) * t168;
t122 = mrSges(7,3) * t194 - t438;
t436 = mrSges(6,3) * t168;
t123 = -mrSges(6,2) * t194 - t436;
t396 = t122 + t123;
t169 = t280 * t199 + t276 * t381;
t435 = mrSges(6,3) * t169;
t124 = mrSges(6,1) * t194 - t435;
t437 = mrSges(7,2) * t169;
t125 = -mrSges(7,1) * t194 + t437;
t395 = -t124 + t125;
t283 = -pkin(1) - pkin(7);
t236 = qJD(1) * t283 + qJD(2);
t183 = -pkin(8) * t392 + t236 * t278;
t180 = t456 * t183;
t218 = t457 * t236;
t352 = qJD(1) * t457;
t184 = -pkin(8) * t352 + t218;
t181 = qJD(3) * pkin(3) + t184;
t128 = t277 * t181 + t180;
t114 = pkin(9) * t381 + t128;
t224 = pkin(3) * t392 + qJD(1) * qJ(2);
t130 = -pkin(4) * t198 - pkin(9) * t199 + t224;
t47 = -t114 * t276 + t130 * t280;
t48 = t114 * t280 + t130 * t276;
t317 = t276 * t48 + t280 * t47;
t273 = qJDD(3) + qJDD(4);
t235 = qJDD(1) * t283 + qJDD(2);
t391 = qJD(3) * t278;
t170 = t457 * t235 - t236 * t391;
t384 = qJD(1) * qJD(3);
t210 = qJDD(1) * t457 - t278 * t384;
t136 = qJDD(3) * pkin(3) - pkin(8) * t210 + t170;
t351 = qJD(3) * t457;
t171 = t278 * t235 + t236 * t351;
t211 = -qJD(1) * t351 - t278 * qJDD(1);
t148 = pkin(8) * t211 + t171;
t350 = qJD(4) * t456;
t390 = qJD(4) * t277;
t34 = t277 * t136 + t456 * t148 + t181 * t350 - t183 * t390;
t31 = pkin(9) * t273 + t34;
t387 = qJD(5) * t280;
t388 = qJD(5) * t276;
t119 = qJD(4) * t198 + t210 * t456 + t277 * t211;
t297 = -t277 * t278 + t336;
t120 = -qJD(1) * qJD(4) * t297 - t277 * t210 + t211 * t456;
t385 = qJD(1) * qJD(2);
t237 = qJDD(1) * qJ(2) + t385;
t178 = -pkin(3) * t211 + t237;
t39 = -pkin(4) * t120 - pkin(9) * t119 + t178;
t6 = -t114 * t388 + t130 * t387 + t276 * t39 + t280 * t31;
t7 = -qJD(5) * t48 - t276 * t31 + t280 * t39;
t332 = -t276 * t7 + t280 * t6;
t117 = qJDD(5) - t120;
t67 = qJD(5) * t169 + t276 * t119 - t280 * t273;
t27 = -mrSges(7,2) * t67 + mrSges(7,3) * t117;
t389 = qJD(5) * t168;
t66 = t280 * t119 + t276 * t273 - t389;
t28 = mrSges(6,1) * t117 - mrSges(6,3) * t66;
t29 = -t117 * mrSges(7,1) + t66 * mrSges(7,2);
t30 = -mrSges(6,2) * t117 - mrSges(6,3) * t67;
t514 = (-t28 + t29) * t276 + (t27 + t30) * t280;
t534 = m(6) * (-qJD(5) * t317 + t332) - t396 * t388 + t395 * t387 + t514;
t533 = t457 / 0.2e1;
t424 = t199 * mrSges(5,3);
t494 = mrSges(5,1) * t381 - mrSges(6,1) * t168 - mrSges(6,2) * t169 - t424;
t412 = t198 * t280;
t413 = t198 * t276;
t532 = -qJD(6) * t276 + (-t387 + t412) * qJ(6) + (t388 - t413) * pkin(5);
t2 = qJ(6) * t117 + qJD(6) * t194 + t6;
t4 = -pkin(5) * t117 + qJDD(6) - t7;
t333 = t2 * t280 + t276 * t4;
t38 = -pkin(5) * t194 + qJD(6) - t47;
t40 = qJ(6) * t194 + t48;
t531 = t38 * t387 - t40 * t388 + t333;
t132 = t184 * t277 + t180;
t530 = pkin(3) * t390 - t132;
t529 = t117 * t502 + t503 * t66 - t535 * t67;
t528 = -t168 * t500 + t169 * t502 + t501 * t194;
t167 = Ifges(6,4) * t168;
t427 = t168 * Ifges(7,5);
t499 = t169 * t503 + t502 * t194 - t167 + t427;
t473 = -t67 / 0.2e1;
t527 = Ifges(6,2) * t473;
t526 = t381 * Ifges(5,5);
t525 = t381 * Ifges(5,6);
t524 = t530 + t532;
t275 = qJ(3) + qJ(4);
t268 = cos(t275);
t281 = cos(qJ(1));
t402 = t276 * t281;
t267 = sin(t275);
t406 = t267 * t281;
t523 = -t268 * mrSges(6,2) * t402 - mrSges(5,2) * t406;
t522 = -t276 * t500 + t280 * t502;
t430 = Ifges(7,5) * t276;
t432 = Ifges(6,4) * t276;
t521 = t280 * t503 + t430 - t432;
t520 = -t47 * t387 - t48 * t388 + t332;
t421 = t280 * mrSges(6,1);
t517 = mrSges(6,2) * t276 - t421;
t279 = sin(qJ(1));
t448 = g(2) * t281;
t516 = -g(1) * t279 + t448;
t470 = t117 / 0.2e1;
t474 = t66 / 0.2e1;
t515 = Ifges(6,4) * t474 + Ifges(6,6) * t470 - t66 * Ifges(7,5) / 0.2e1 - t117 * Ifges(7,6) / 0.2e1 + Ifges(7,3) * t473 + t527;
t179 = t277 * t183;
t127 = t181 * t456 - t179;
t113 = -pkin(4) * t381 - t127;
t328 = t276 * mrSges(7,1) - t280 * mrSges(7,3);
t330 = mrSges(6,1) * t276 + mrSges(6,2) * t280;
t44 = t168 * pkin(5) - t169 * qJ(6) + t113;
t513 = t113 * t330 + t328 * t44;
t308 = mrSges(4,1) * t457 - mrSges(4,2) * t278;
t374 = Ifges(4,4) * t457;
t512 = (-Ifges(4,1) * t278 - t374) * t533 + qJ(2) * t308;
t398 = t279 * t280;
t368 = t268 * t398;
t403 = t276 * t279;
t369 = t268 * t403;
t405 = t268 * t279;
t407 = t267 * t279;
t511 = -mrSges(5,1) * t405 - t368 * t506 + t505 * t369 - t504 * t407;
t321 = -t280 * pkin(5) - t276 * qJ(6);
t221 = -pkin(4) + t321;
t329 = -t280 * mrSges(7,1) - t276 * mrSges(7,3);
t288 = m(7) * t221 + t329;
t331 = mrSges(5,1) * t267 + mrSges(5,2) * t268;
t510 = t331 - t504 * t268 + (-t288 - t517) * t267;
t433 = Ifges(5,4) * t199;
t143 = Ifges(5,2) * t198 + t433 + t525;
t192 = Ifges(5,4) * t198;
t144 = Ifges(5,1) * t199 + t192 + t526;
t35 = t136 * t456 - t277 * t148 - t181 * t390 - t183 * t350;
t32 = -t273 * pkin(4) - t35;
t429 = Ifges(7,5) * t280;
t322 = Ifges(7,3) * t276 + t429;
t431 = Ifges(6,4) * t280;
t325 = -Ifges(6,2) * t276 + t431;
t346 = t387 / 0.2e1;
t347 = t388 / 0.2e1;
t426 = t169 * Ifges(6,4);
t87 = -t168 * Ifges(6,2) + t194 * Ifges(6,6) + t426;
t367 = -t276 * t87 / 0.2e1;
t425 = t198 * mrSges(5,3);
t460 = t199 / 0.2e1;
t462 = -t198 / 0.2e1;
t464 = -t194 / 0.2e1;
t466 = -t169 / 0.2e1;
t467 = t168 / 0.2e1;
t468 = -t168 / 0.2e1;
t472 = t67 / 0.2e1;
t166 = Ifges(7,5) * t169;
t84 = t194 * Ifges(7,6) + t168 * Ifges(7,3) + t166;
t9 = t67 * pkin(5) - t66 * qJ(6) - t169 * qJD(6) + t32;
t509 = (-t429 + t431) * t474 + (-t38 * t412 + t40 * t413 + t531) * mrSges(7,2) + (t87 / 0.2e1 - t84 / 0.2e1) * t413 - (Ifges(5,1) * t198 - t433 + t528) * t199 / 0.2e1 + (-Ifges(7,3) * t472 + t470 * t500 + t515 + t527) * t280 + (t48 * mrSges(6,2) - t47 * mrSges(6,1) + t38 * mrSges(7,1) - t40 * mrSges(7,3) - t224 * mrSges(5,1) + Ifges(6,6) * t467 + Ifges(7,6) * t468 - Ifges(5,2) * t462 + t525 / 0.2e1 + t502 * t466 + t501 * t464) * t199 + (-t224 * mrSges(5,2) + t325 * t467 + t322 * t468 - t526 / 0.2e1 + t521 * t466 + t522 * t464 - t513) * t198 + (t412 * t47 + t413 * t48 + t520) * mrSges(6,3) + (t367 + t513) * qJD(5) + t32 * t517 + t143 * t460 + t128 * t424 + t127 * t425 + (-t325 / 0.2e1 + t322 / 0.2e1) * t389 + t430 * t472 + t432 * t473 + (t144 + t192) * t462 + (t169 * t521 + t194 * t522) * qJD(5) / 0.2e1 - (-t267 * mrSges(7,2) + t288 * t268) * t448 + (t529 / 0.2e1 + t502 * t470 + t503 * t474) * t276 + t9 * t329 + Ifges(5,5) * t119 + Ifges(5,6) * t120 - t34 * mrSges(5,2) + t35 * mrSges(5,1) + (-t412 / 0.2e1 + t346) * t499 + t84 * t347 + Ifges(5,3) * t273;
t471 = -m(3) - m(4);
t507 = -m(6) - m(7);
t497 = -t128 + t532;
t271 = t278 * pkin(3);
t253 = qJ(2) + t271;
t158 = -pkin(4) * t298 - pkin(9) * t297 + t253;
t440 = pkin(8) - t283;
t219 = t440 * t278;
t361 = t457 * t283;
t220 = -pkin(8) * t457 + t361;
t165 = -t219 * t456 + t277 * t220;
t495 = t276 * t158 + t280 * t165;
t493 = t277 * t219 + t456 * t220;
t491 = t117 * t501 - t500 * t67 + t502 * t66;
t490 = -t170 * t457 - t171 * t278;
t162 = -qJD(3) * t360 - qJD(4) * t363 - t277 * t351 - t278 * t350;
t163 = -t277 * t391 - t278 * t390 + t336 * t381;
t489 = -t127 * t162 - t128 * t163 - t297 * t35;
t24 = mrSges(6,1) * t67 + mrSges(6,2) * t66;
t488 = -m(6) * t32 - t24;
t485 = -m(4) * t490 + t278 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t211);
t484 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - mrSges(5,3);
t305 = t278 * mrSges(4,1) + mrSges(4,2) * t457;
t482 = -t305 - mrSges(3,3) - t331 + mrSges(2,2);
t238 = pkin(3) * t351 + qJD(2);
t91 = pkin(4) * t163 - pkin(9) * t162 + t238;
t202 = t440 * t391;
t203 = t220 * qJD(3);
t95 = qJD(4) * t493 + t277 * t202 + t456 * t203;
t21 = -qJD(5) * t495 - t276 * t95 + t280 * t91;
t481 = -m(7) * pkin(5) - t506;
t480 = -m(7) * qJ(6) + t505;
t479 = -m(6) * t113 + t494;
t478 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t477 = qJD(1) ^ 2;
t465 = t169 / 0.2e1;
t455 = pkin(3) * t277;
t454 = pkin(4) * t268;
t453 = pkin(5) * t199;
t452 = pkin(9) * t267;
t447 = g(3) * t268;
t445 = t297 * t9;
t434 = Ifges(4,4) * t278;
t428 = pkin(3) * qJD(4);
t423 = t297 * t32;
t418 = t162 * t276;
t417 = t162 * t280;
t416 = t163 * t276;
t415 = t163 * t280;
t410 = t297 * t280;
t404 = t268 * t281;
t223 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t352;
t399 = t278 * t223;
t397 = t281 * t280;
t157 = t199 * pkin(4) - t198 * pkin(9);
t69 = t280 * t127 + t276 * t157;
t133 = t184 * t456 - t179;
t338 = pkin(3) * t352;
t142 = t338 + t157;
t60 = t280 * t133 + t276 * t142;
t394 = pkin(4) * t405 + pkin(9) * t407;
t393 = t281 * pkin(1) + t279 * qJ(2);
t386 = qJDD(1) * mrSges(3,2);
t378 = t457 * pkin(3);
t377 = t456 * pkin(3);
t375 = t268 * t421;
t364 = t276 * t456;
t362 = t280 * t456;
t270 = t281 * qJ(2);
t345 = -pkin(1) * t279 + t270;
t254 = t268 * pkin(9);
t344 = -pkin(4) * t267 + t254;
t343 = -t384 / 0.2e1;
t342 = m(5) * t378;
t340 = (t237 + t385) * qJ(2);
t339 = pkin(3) * t350;
t337 = pkin(5) * t368 + qJ(6) * t369 + t394;
t320 = pkin(5) * t276 - qJ(6) * t280;
t318 = t276 * t40 - t280 * t38;
t315 = t276 * t339;
t314 = t280 * t339;
t68 = -t127 * t276 + t157 * t280;
t59 = -t133 * t276 + t142 * t280;
t92 = t158 * t280 - t165 * t276;
t282 = -pkin(8) - pkin(7);
t312 = t281 * t271 + t279 * t282 + t345;
t310 = t279 * t271 - t281 * t282 + t393;
t309 = -t378 - t452;
t307 = t457 * Ifges(4,1) - t434;
t306 = -Ifges(4,2) * t278 + t374;
t304 = -Ifges(4,5) * t278 - Ifges(4,6) * t457;
t301 = -t297 * t387 - t418;
t300 = -t297 * t388 + t417;
t20 = t158 * t387 - t165 * t388 + t276 * t91 + t280 * t95;
t295 = t278 * (-Ifges(4,2) * t457 - t434);
t294 = -t276 * t396 + t280 * t395;
t290 = -qJD(5) * t318 + t333;
t96 = qJD(4) * t165 - t456 * t202 + t277 * t203;
t265 = -pkin(1) * qJDD(1) + qJDD(2);
t258 = -t377 - pkin(4);
t248 = t279 * t378;
t222 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t392;
t207 = t305 * qJD(1);
t200 = -t377 + t221;
t197 = Ifges(4,5) * qJD(3) + qJD(1) * t307;
t196 = Ifges(4,6) * qJD(3) + qJD(1) * t306;
t191 = t267 * t397 - t403;
t190 = t267 * t402 + t398;
t189 = t267 * t398 + t402;
t188 = t267 * t403 - t397;
t187 = t199 * qJ(6);
t185 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t210;
t174 = -mrSges(5,2) * t381 + t425;
t156 = -mrSges(5,1) * t198 + mrSges(5,2) * t199;
t108 = -mrSges(5,2) * t273 + mrSges(5,3) * t120;
t107 = mrSges(5,1) * t273 - mrSges(5,3) * t119;
t105 = mrSges(7,1) * t168 - mrSges(7,3) * t169;
t104 = pkin(5) * t169 + qJ(6) * t168;
t99 = t297 * t320 - t493;
t76 = pkin(5) * t298 - t92;
t75 = -qJ(6) * t298 + t495;
t46 = -t68 - t453;
t45 = t187 + t69;
t43 = -t59 - t453;
t42 = t187 + t60;
t23 = mrSges(7,1) * t67 - mrSges(7,3) * t66;
t22 = t320 * t162 - (qJD(5) * t321 + qJD(6) * t280) * t297 + t96;
t11 = -pkin(5) * t163 - t21;
t10 = qJ(6) * t163 - qJD(6) * t298 + t20;
t1 = [(t502 * t163 + t503 * t300 + t301 * t535) * t465 + m(6) * (t20 * t48 + t21 * t47 + t495 * t6 + t7 * t92) + t495 * t30 + (-t491 / 0.2e1 - t178 * mrSges(5,1) + t34 * mrSges(5,3) + Ifges(5,4) * t119 + Ifges(5,2) * t120 + Ifges(5,6) * t273 - Ifges(6,6) * t473 - Ifges(7,6) * t472 - t470 * t501 - t474 * t502 - t478) * t298 + (t305 + 0.2e1 * mrSges(3,3)) * t237 + t490 * mrSges(4,3) + t489 * mrSges(5,3) + (t222 * t351 - t223 * t391 + t485) * t283 + (-m(5) * t127 - t479) * t96 + t528 * t163 / 0.2e1 + (qJD(5) * t84 + t529) * t410 / 0.2e1 + (-(mrSges(7,2) * t2 + mrSges(6,3) * t6 + t515) * t276 - t499 * t347 + t178 * mrSges(5,2) + Ifges(5,1) * t119 + Ifges(5,4) * t120 + Ifges(5,5) * t273 + t322 * t472 + t325 * t473 - t346 * t87 + t521 * t474 + t522 * t470) * t297 + t512 * t384 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + (Ifges(7,5) * t300 + Ifges(7,6) * t163 - Ifges(7,3) * t301) * t467 + (Ifges(6,4) * t300 + Ifges(6,2) * t301 + Ifges(6,6) * t163) * t468 + (Ifges(5,1) * t162 - Ifges(5,4) * t163) * t460 + m(5) * (t128 * t95 + t165 * t34 + t178 * t253 + t224 * t238) + t381 * (Ifges(5,5) * t162 - Ifges(5,6) * t163) / 0.2e1 + qJDD(3) * (Ifges(4,5) * t457 - Ifges(4,6) * t278) + t328 * t445 + t330 * t423 + (Ifges(4,1) * t210 + Ifges(4,4) * t211) * t533 - t196 * t351 / 0.2e1 + t162 * t367 + t84 * t418 / 0.2e1 - (-m(5) * t35 - t107 - t488) * t493 + m(4) * t340 - t197 * t391 / 0.2e1 - pkin(1) * t386 - t278 * (Ifges(4,4) * t210 + Ifges(4,2) * t211) / 0.2e1 + qJD(3) ^ 2 * t304 / 0.2e1 + t211 * t306 / 0.2e1 + t210 * t307 / 0.2e1 + t113 * (-mrSges(6,1) * t301 + mrSges(6,2) * t300) + t44 * (-mrSges(7,1) * t301 - mrSges(7,3) * t300) + t38 * (-mrSges(7,1) * t163 + mrSges(7,2) * t300) + t47 * (mrSges(6,1) * t163 - mrSges(6,3) * t300) + t48 * (-mrSges(6,2) * t163 + mrSges(6,3) * t301) + t40 * (mrSges(7,2) * t301 + mrSges(7,3) * t163) + t99 * t23 + t22 * t105 + t10 * t122 + t20 * t123 + t21 * t124 + t11 * t125 + m(7) * (t10 * t40 + t11 * t38 + t2 * t75 + t22 * t44 + t4 * t76 + t9 * t99) + t499 * t417 / 0.2e1 + t162 * t144 / 0.2e1 + (t163 * t501 + t300 * t502 + t301 * t500) * t194 / 0.2e1 - t163 * t143 / 0.2e1 + t165 * t108 + t95 * t174 + (-m(3) * t345 - m(4) * t270 - m(5) * t312 + t504 * t404 + t507 * (pkin(4) * t406 - pkin(9) * t404 + t312) + t481 * t191 + t480 * t190 + t482 * t281 + (-m(4) * t283 - t484) * t279) * g(1) + (-m(5) * t310 + t504 * t405 + t471 * t393 + t507 * (pkin(4) * t407 - pkin(9) * t405 + t310) + t481 * t189 + t480 * t188 + (-m(4) * pkin(7) + t484) * t281 + t482 * t279) * g(2) + t295 * t343 + t185 * t361 + t198 * (Ifges(5,4) * t162 - Ifges(5,2) * t163) / 0.2e1 + qJD(2) * t207 + qJ(2) * (-mrSges(4,1) * t211 + mrSges(4,2) * t210) + t224 * (mrSges(5,1) * t163 + mrSges(5,2) * t162) + t238 * t156 + t253 * (-mrSges(5,1) * t120 + mrSges(5,2) * t119) + t265 * mrSges(3,2) - (-mrSges(7,2) * t4 + mrSges(6,3) * t7) * t410 + t75 * t27 + t76 * t29 + t92 * t28 + m(3) * (-pkin(1) * t265 + t340); t386 + t457 * t185 + (t222 * t457 - t399) * qJD(3) + (qJ(2) * t471 - mrSges(3,3)) * t477 - (t23 + t24 - t107) * t297 + (-t105 + t494) * t162 + (t276 * t395 + t280 * t396 + t174) * t163 + m(3) * t265 - m(5) * t489 + m(7) * (-t162 * t44 + t38 * t416 + t40 * t415 - t445) + m(6) * (-t113 * t162 + t415 * t48 - t416 * t47 - t423) - (m(5) * t34 + m(6) * t520 + m(7) * t531 + t294 * qJD(5) + t108 + t514) * t298 + (-m(5) * t224 - m(6) * t317 - m(7) * t318 - t156 - t207 + t294) * qJD(1) + t516 * (m(5) - t471 - t507) + t485; -t494 * t530 + (t32 * t258 + (t113 * t277 + t362 * t48 - t364 * t47) * t428 - t113 * t132 - t47 * t59 - t48 * t60) * m(6) + (-t309 * t448 + t9 * t200 + (t362 * t40 + t364 * t38) * t428 - t38 * t43 - t40 * t42 + t524 * t44) * m(7) + t524 * t105 + (mrSges(6,3) * t406 + (-m(6) * (t309 - t454) + t375 + mrSges(5,1) * t268 + t342) * t281 + t523) * g(2) + t516 * t308 + (-m(7) * (t254 - t271) + m(5) * t271 - m(6) * (-t271 + t344) + t305 + t510) * g(3) + (-m(6) * (t248 + t394) - m(7) * (t248 + t337) - (-mrSges(5,2) * t267 + t342) * t279 + t511) * g(1) + (t295 / 0.2e1 - t512) * t477 - t156 * t338 + t108 * t455 + t196 * t352 / 0.2e1 + t197 * t392 / 0.2e1 + (m(7) * t290 + t534) * (pkin(9) + t455) + (t127 * t132 - t128 * t133 - t224 * t338 + (t456 * t35 + t277 * t34 + (-t127 * t277 + t128 * t456) * qJD(4)) * pkin(3)) * m(5) - t222 * t218 + t509 + t236 * t399 + (-t43 + t315) * t125 + (-t315 - t59) * t124 + (-t60 + t314) * t123 + Ifges(4,3) * qJDD(3) + (-t42 + t314) * t122 + t170 * mrSges(4,1) - t171 * mrSges(4,2) + t304 * t343 + t107 * t377 + t200 * t23 + Ifges(4,5) * t210 + Ifges(4,6) * t211 + t258 * t24 + (t339 - t133) * t174; -m(6) * (t47 * t68 + t48 * t69) - t45 * t122 - t69 * t123 - t68 * t124 - t46 * t125 - t127 * t174 + t221 * t23 + t479 * t128 + t497 * t105 + (-(m(6) * (-t452 - t454) - t375 - t267 * mrSges(6,3)) * t281 + mrSges(5,1) * t404 + t523) * g(2) + t488 * pkin(4) + (t221 * t9 - t38 * t46 - t40 * t45 + t44 * t497) * m(7) + (-m(6) * t344 - m(7) * t254 + t510) * g(3) + (-m(6) * t394 - m(7) * t337 + mrSges(5,2) * t407 + t511) * g(1) + ((t267 * t448 + t290) * m(7) + t534) * pkin(9) + t509; t478 + (Ifges(7,3) * t169 - t427) * t468 + t87 * t465 + t40 * t437 + t38 * t438 + (-t168 * t502 - t169 * t500) * t464 + (t188 * t506 + t189 * t505) * g(1) + (-Ifges(6,2) * t169 - t167 + t499) * t467 + (-t168 * t503 + t166 - t426 + t84) * t466 + (-m(7) * t38 - t395 + t435) * t48 + (-t190 * t506 - t191 * t505) * g(2) + (t330 + t328) * t447 + t491 + (-t104 * t44 - pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t40 + t320 * t447 - g(1) * (-pkin(5) * t188 + qJ(6) * t189) - g(2) * (pkin(5) * t190 - qJ(6) * t191)) * m(7) + (-m(7) * t40 - t396 - t436) * t47 - t104 * t105 + qJD(6) * t122 + qJ(6) * t27 - pkin(5) * t29 - t44 * (mrSges(7,1) * t169 + mrSges(7,3) * t168) - t113 * (mrSges(6,1) * t169 - mrSges(6,2) * t168); t169 * t105 - t194 * t122 + (-g(1) * t188 + g(2) * t190 + t44 * t169 - t40 * t194 - t276 * t447 + t4) * m(7) + t29;];
tau  = t1;
