% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:39
% EndTime: 2019-03-09 03:41:16
% DurationCPUTime: 23.56s
% Computational Cost: add. (11757->758), mult. (25232->1022), div. (0->0), fcn. (17829->18), ass. (0->326)
t287 = cos(qJ(3));
t276 = sin(pkin(11));
t278 = cos(pkin(11));
t282 = sin(qJ(5));
t286 = cos(qJ(5));
t309 = t276 * t282 - t278 * t286;
t303 = t309 * t287;
t186 = qJD(1) * t303;
t206 = t309 * qJD(5);
t479 = t206 - t186;
t225 = t276 * t286 + t278 * t282;
t304 = t225 * t287;
t185 = qJD(1) * t304;
t207 = t225 * qJD(5);
t474 = -t207 + t185;
t478 = mrSges(6,3) + mrSges(7,3);
t277 = sin(pkin(10));
t258 = pkin(1) * t277 + pkin(7);
t241 = t258 * qJD(1);
t283 = sin(qJ(3));
t228 = t283 * t241;
t198 = qJD(2) * t287 - t228;
t314 = pkin(3) * t283 - qJ(4) * t287;
t232 = t314 * qJD(1);
t139 = -t198 * t276 + t278 * t232;
t361 = t278 * t287;
t308 = pkin(4) * t283 - pkin(8) * t361;
t101 = qJD(1) * t308 + t139;
t140 = t278 * t198 + t276 * t232;
t350 = t287 * qJD(1);
t340 = t276 * t350;
t113 = -pkin(8) * t340 + t140;
t389 = pkin(8) + qJ(4);
t237 = t389 * t276;
t238 = t389 * t278;
t351 = qJD(5) * t286;
t354 = qJD(4) * t278;
t355 = qJD(4) * t276;
t462 = -t237 * t351 + (-t113 + t354) * t286 + (-qJD(5) * t238 - t101 - t355) * t282;
t166 = -t282 * t237 + t286 * t238;
t461 = -t225 * qJD(4) - qJD(5) * t166 - t286 * t101 + t113 * t282;
t358 = qJD(1) * t283;
t477 = -pkin(5) * t358 + t479 * pkin(9) + t461;
t476 = -t474 * pkin(9) - t462;
t259 = pkin(4) * t278 + pkin(3);
t274 = pkin(11) + qJ(5);
t266 = cos(t274);
t230 = pkin(5) * t266 + t259;
t270 = qJ(6) + t274;
t256 = sin(t270);
t257 = cos(t270);
t264 = sin(t274);
t473 = -m(6) * t259 - m(7) * t230 - t266 * mrSges(6,1) - t257 * mrSges(7,1) + t264 * mrSges(6,2) + t256 * mrSges(7,2);
t273 = -pkin(9) - t389;
t472 = -m(6) * t389 + m(7) * t273 - t478;
t252 = qJD(5) - t350;
t218 = qJD(3) * t278 - t276 * t358;
t219 = qJD(3) * t276 + t278 * t358;
t325 = t286 * t218 - t219 * t282;
t155 = t218 * t282 + t219 * t286;
t382 = Ifges(6,4) * t155;
t76 = Ifges(6,2) * t325 + t252 * Ifges(6,6) + t382;
t424 = t76 / 0.2e1;
t149 = Ifges(6,4) * t325;
t77 = t155 * Ifges(6,1) + t252 * Ifges(6,5) + t149;
t423 = t77 / 0.2e1;
t349 = m(5) + m(6) + m(7);
t449 = -m(4) - t349;
t471 = qJD(2) * qJD(3) + t258 * qJDD(1);
t275 = qJ(1) + pkin(10);
t265 = sin(t275);
t267 = cos(t275);
t470 = g(1) * t267 + g(2) * t265;
t281 = sin(qJ(6));
t285 = cos(qJ(6));
t469 = -t155 * t281 + t285 * t325;
t85 = t155 * t285 + t281 * t325;
t348 = qJD(1) * qJD(3);
t234 = qJDD(1) * t283 + t287 * t348;
t187 = qJDD(3) * t278 - t234 * t276;
t188 = qJDD(3) * t276 + t234 * t278;
t73 = qJD(5) * t325 + t187 * t282 + t188 * t286;
t74 = -qJD(5) * t155 + t187 * t286 - t188 * t282;
t24 = qJD(6) * t469 + t281 * t74 + t285 * t73;
t434 = t24 / 0.2e1;
t25 = -qJD(6) * t85 - t281 * t73 + t285 * t74;
t433 = t25 / 0.2e1;
t426 = t73 / 0.2e1;
t425 = t74 / 0.2e1;
t413 = t187 / 0.2e1;
t412 = t188 / 0.2e1;
t233 = -t287 * qJDD(1) + t283 * t348;
t223 = qJDD(5) + t233;
t213 = qJDD(6) + t223;
t411 = t213 / 0.2e1;
t410 = t223 / 0.2e1;
t468 = -t233 / 0.2e1;
t409 = t233 / 0.2e1;
t467 = t234 / 0.2e1;
t466 = qJD(3) / 0.2e1;
t165 = -t286 * t237 - t238 * t282;
t131 = -pkin(9) * t225 + t165;
t132 = -pkin(9) * t309 + t166;
t64 = t131 * t285 - t132 * t281;
t465 = qJD(6) * t64 + t281 * t477 - t476 * t285;
t65 = t131 * t281 + t132 * t285;
t464 = -qJD(6) * t65 + t476 * t281 + t285 * t477;
t435 = m(7) * pkin(5);
t463 = mrSges(6,1) + t435;
t322 = -mrSges(5,1) * t278 + mrSges(5,2) * t276;
t302 = m(5) * pkin(3) - t322;
t459 = t287 * t302;
t329 = -qJ(4) * t283 - pkin(2);
t279 = cos(pkin(10));
t401 = pkin(1) * t279;
t214 = -pkin(3) * t287 + t329 - t401;
t197 = t278 * t214;
t362 = t278 * t283;
t130 = -pkin(8) * t362 + t197 + (-t258 * t276 - pkin(4)) * t287;
t159 = t276 * t214 + t258 * t361;
t364 = t276 * t283;
t138 = -pkin(8) * t364 + t159;
t67 = t282 * t130 + t286 * t138;
t344 = mrSges(4,3) * t358;
t456 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t218 + mrSges(5,2) * t219 + t344;
t383 = Ifges(5,4) * t278;
t317 = -Ifges(5,2) * t276 + t383;
t384 = Ifges(5,4) * t276;
t319 = Ifges(5,1) * t278 - t384;
t455 = t218 * (Ifges(5,6) * t283 + t287 * t317) + t219 * (Ifges(5,5) * t283 + t287 * t319);
t144 = -mrSges(5,2) * t233 + mrSges(5,3) * t187;
t145 = mrSges(5,1) * t233 - mrSges(5,3) * t188;
t454 = t144 * t278 - t145 * t276;
t342 = t283 * qJDD(2) + t287 * t471;
t357 = qJD(3) * t283;
t136 = -t241 * t357 + t342;
t356 = qJD(3) * t287;
t137 = qJDD(2) * t287 - t241 * t356 - t283 * t471;
t453 = t136 * t287 - t137 * t283;
t118 = qJDD(3) * qJ(4) + (qJD(4) - t228) * qJD(3) + t342;
t260 = -pkin(2) - t401;
t240 = t260 * qJDD(1);
t353 = qJD(4) * t283;
t129 = pkin(3) * t233 - qJ(4) * t234 - qJD(1) * t353 + t240;
t62 = -t118 * t276 + t278 * t129;
t63 = t278 * t118 + t276 * t129;
t313 = -t276 * t62 + t278 * t63;
t199 = t283 * qJD(2) + t287 * t241;
t167 = pkin(4) * t340 + t199;
t452 = -pkin(5) * t474 - t167;
t115 = -t187 * mrSges(5,1) + t188 * mrSges(5,2);
t39 = -t74 * mrSges(6,1) + t73 * mrSges(6,2);
t8 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t451 = -t115 - t39 - t8;
t249 = qJD(6) + t252;
t450 = t219 * Ifges(5,5) + t155 * Ifges(6,5) + t85 * Ifges(7,5) + t218 * Ifges(5,6) + Ifges(6,6) * t325 + Ifges(7,6) * t469 - Ifges(5,3) * t350 + t252 * Ifges(6,3) + t249 * Ifges(7,3);
t182 = qJD(3) * qJ(4) + t199;
t189 = t214 * qJD(1);
t104 = -t182 * t276 + t278 * t189;
t91 = -pkin(4) * t350 - pkin(8) * t219 + t104;
t105 = t278 * t182 + t276 * t189;
t93 = pkin(8) * t218 + t105;
t45 = -t282 * t93 + t286 * t91;
t37 = -pkin(9) * t155 + t45;
t36 = pkin(5) * t252 + t37;
t46 = t282 * t91 + t286 * t93;
t38 = pkin(9) * t325 + t46;
t375 = t285 * t38;
t14 = t281 * t36 + t375;
t180 = -qJD(3) * pkin(3) + qJD(4) - t198;
t141 = -pkin(4) * t218 + t180;
t89 = -pkin(5) * t325 + t141;
t448 = -mrSges(7,1) * t89 + mrSges(7,3) * t14;
t378 = t281 * t38;
t13 = t285 * t36 - t378;
t447 = mrSges(7,2) * t89 - t13 * mrSges(7,3);
t263 = Ifges(4,4) * t350;
t446 = t278 * (t219 * Ifges(5,1) + t218 * Ifges(5,4) - Ifges(5,5) * t350) + Ifges(4,1) * t358 + Ifges(4,5) * qJD(3) + t263;
t324 = mrSges(4,1) * t287 - mrSges(4,2) * t283;
t445 = t283 * t478 + mrSges(3,1) + t324;
t321 = mrSges(5,1) * t276 + t278 * mrSges(5,2);
t323 = mrSges(4,1) * t283 + mrSges(4,2) * t287;
t363 = t276 * t287;
t444 = -t180 * t287 * t321 - t105 * (-mrSges(5,2) * t283 - mrSges(5,3) * t363) - t104 * (mrSges(5,1) * t283 - mrSges(5,3) * t361) - t260 * qJD(1) * t323;
t352 = qJD(5) * t282;
t52 = pkin(4) * t233 - pkin(8) * t188 + t62;
t55 = pkin(8) * t187 + t63;
t11 = t282 * t52 + t286 * t55 + t91 * t351 - t352 * t93;
t10 = pkin(9) * t74 + t11;
t12 = -qJD(5) * t46 - t282 * t55 + t286 * t52;
t9 = pkin(5) * t223 - pkin(9) * t73 + t12;
t2 = qJD(6) * t13 + t10 * t285 + t281 * t9;
t3 = -qJD(6) * t14 - t10 * t281 + t285 * t9;
t442 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t441 = -t12 * mrSges(6,1) + t11 * mrSges(6,2);
t396 = pkin(5) * t264;
t399 = pkin(4) * t276;
t440 = -m(6) * t399 - m(7) * (t396 + t399) + mrSges(3,2) - mrSges(4,3) - t321;
t386 = Ifges(4,4) * t283;
t318 = Ifges(4,2) * t287 + t386;
t439 = t13 * mrSges(7,1) + t45 * mrSges(6,1) - Ifges(4,6) * qJD(3) / 0.2e1 - qJD(1) * t318 / 0.2e1 - t14 * mrSges(7,2) - t46 * mrSges(6,2);
t437 = Ifges(7,4) * t434 + Ifges(7,2) * t433 + Ifges(7,6) * t411;
t436 = Ifges(7,1) * t434 + Ifges(7,4) * t433 + Ifges(7,5) * t411;
t432 = Ifges(6,4) * t426 + Ifges(6,2) * t425 + Ifges(6,6) * t410;
t431 = Ifges(6,1) * t426 + Ifges(6,4) * t425 + Ifges(6,5) * t410;
t402 = Ifges(7,4) * t85;
t41 = Ifges(7,2) * t469 + t249 * Ifges(7,6) + t402;
t430 = -t41 / 0.2e1;
t429 = t41 / 0.2e1;
t78 = Ifges(7,4) * t469;
t42 = t85 * Ifges(7,1) + t249 * Ifges(7,5) + t78;
t428 = -t42 / 0.2e1;
t427 = t42 / 0.2e1;
t422 = -t469 / 0.2e1;
t421 = t469 / 0.2e1;
t420 = -t85 / 0.2e1;
t419 = t85 / 0.2e1;
t418 = Ifges(5,1) * t412 + Ifges(5,4) * t413 + Ifges(5,5) * t409;
t417 = -t325 / 0.2e1;
t416 = t325 / 0.2e1;
t415 = -t155 / 0.2e1;
t414 = t155 / 0.2e1;
t408 = -t249 / 0.2e1;
t407 = t249 / 0.2e1;
t406 = -t252 / 0.2e1;
t405 = t252 / 0.2e1;
t284 = sin(qJ(1));
t400 = pkin(1) * t284;
t398 = pkin(5) * t155;
t393 = g(3) * t283;
t288 = cos(qJ(1));
t272 = t288 * pkin(1);
t388 = mrSges(6,3) * t325;
t387 = mrSges(6,3) * t155;
t385 = Ifges(4,4) * t287;
t127 = -qJDD(3) * pkin(3) + qJDD(4) - t137;
t373 = t127 * t283;
t366 = t265 * t287;
t365 = t267 * t287;
t246 = t283 * t258;
t168 = t256 * t366 + t257 * t267;
t169 = t256 * t267 - t257 * t366;
t360 = -t168 * mrSges(7,1) + t169 * mrSges(7,2);
t170 = -t256 * t365 + t257 * t265;
t171 = t256 * t265 + t257 * t365;
t359 = t170 * mrSges(7,1) - t171 * mrSges(7,2);
t203 = qJD(3) * t314 - t353;
t339 = t258 * t357;
t150 = t278 * t203 + t276 * t339;
t236 = t258 * t356;
t338 = t276 * t356;
t193 = pkin(4) * t338 + t236;
t202 = pkin(4) * t364 + t246;
t346 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t213;
t345 = Ifges(6,5) * t73 + Ifges(6,6) * t74 + Ifges(6,3) * t223;
t343 = mrSges(4,3) * t350;
t337 = m(5) * qJ(4) + mrSges(5,3);
t328 = -t348 / 0.2e1;
t327 = t348 / 0.2e1;
t66 = t286 * t130 - t138 * t282;
t320 = -mrSges(7,1) * t256 - mrSges(7,2) * t257;
t316 = Ifges(4,5) * t287 - Ifges(4,6) * t283;
t315 = Ifges(5,5) * t278 - Ifges(5,6) * t276;
t195 = t309 * t283;
t56 = -pkin(5) * t287 + pkin(9) * t195 + t66;
t194 = t225 * t283;
t57 = -pkin(9) * t194 + t67;
t28 = -t281 * t57 + t285 * t56;
t29 = t281 * t56 + t285 * t57;
t312 = -t104 * t276 + t105 * t278;
t125 = -t194 * t285 + t195 * t281;
t126 = -t194 * t281 - t195 * t285;
t156 = -t225 * t281 - t285 * t309;
t157 = t225 * t285 - t281 * t309;
t311 = t230 * t287 - t273 * t283;
t310 = t259 * t287 + t283 * t389;
t307 = t346 - t442;
t177 = -t264 * t365 + t265 * t266;
t175 = t264 * t366 + t266 * t267;
t305 = t283 * (Ifges(4,1) * t287 - t386);
t111 = qJD(3) * t308 + t150;
t191 = t276 * t203;
t128 = t191 + (-pkin(8) * t363 - t246 * t278) * qJD(3);
t34 = t282 * t111 + t286 * t128 + t130 * t351 - t138 * t352;
t294 = t287 * (Ifges(5,3) * t283 + t287 * t315);
t35 = -qJD(5) * t67 + t286 * t111 - t128 * t282;
t92 = -pkin(4) * t187 + t127;
t293 = t283 * t337 + t459;
t244 = -qJD(3) * mrSges(4,2) + t343;
t201 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t234;
t200 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t233;
t190 = pkin(5) * t309 - t259;
t184 = -mrSges(5,1) * t350 - mrSges(5,3) * t219;
t183 = mrSges(5,2) * t350 + mrSges(5,3) * t218;
t178 = t264 * t265 + t266 * t365;
t176 = t264 * t267 - t266 * t366;
t158 = -t258 * t363 + t197;
t151 = -t278 * t339 + t191;
t147 = t219 * Ifges(5,4) + t218 * Ifges(5,2) - Ifges(5,6) * t350;
t143 = pkin(5) * t194 + t202;
t135 = -qJD(3) * t304 + t206 * t283;
t134 = -qJD(3) * t303 - t207 * t283;
t124 = mrSges(6,1) * t252 - t387;
t123 = -mrSges(6,2) * t252 + t388;
t108 = -t185 * t281 - t186 * t285;
t107 = -t185 * t285 + t186 * t281;
t97 = t188 * Ifges(5,4) + t187 * Ifges(5,2) + t233 * Ifges(5,6);
t94 = -pkin(5) * t135 + t193;
t90 = -mrSges(6,1) * t325 + mrSges(6,2) * t155;
t69 = mrSges(7,1) * t249 - mrSges(7,3) * t85;
t68 = -mrSges(7,2) * t249 + mrSges(7,3) * t469;
t61 = -mrSges(6,2) * t223 + mrSges(6,3) * t74;
t60 = mrSges(6,1) * t223 - mrSges(6,3) * t73;
t48 = -qJD(6) * t126 - t134 * t281 + t135 * t285;
t47 = qJD(6) * t125 + t134 * t285 + t135 * t281;
t44 = -pkin(5) * t74 + t92;
t43 = -mrSges(7,1) * t469 + mrSges(7,2) * t85;
t27 = pkin(9) * t135 + t34;
t26 = pkin(5) * t357 - pkin(9) * t134 + t35;
t18 = -mrSges(7,2) * t213 + mrSges(7,3) * t25;
t17 = mrSges(7,1) * t213 - mrSges(7,3) * t24;
t16 = t285 * t37 - t378;
t15 = -t281 * t37 - t375;
t5 = -qJD(6) * t29 + t26 * t285 - t27 * t281;
t4 = qJD(6) * t28 + t26 * t281 + t27 * t285;
t1 = [(m(4) * (-t198 * t287 - t199 * t283) * t258 + t316 * t466 - t444) * qJD(3) + (-m(3) * t272 - mrSges(2,1) * t288 - t178 * mrSges(6,1) - t171 * mrSges(7,1) + mrSges(2,2) * t284 - t177 * mrSges(6,2) - t170 * mrSges(7,2) + t449 * (t267 * pkin(2) + t265 * pkin(7) + t272) + t440 * t265 + (-m(6) * t310 - m(7) * t311 - t293 - t445) * t267) * g(2) + (Ifges(6,4) * t134 + Ifges(6,2) * t135) * t416 + (Ifges(6,1) * t134 + Ifges(6,4) * t135) * t414 + (-Ifges(6,6) * t425 - Ifges(6,5) * t426 - Ifges(7,6) * t433 - Ifges(7,5) * t434 + Ifges(4,4) * t467 + Ifges(4,2) * t468 - Ifges(5,3) * t409 - Ifges(6,3) * t410 - Ifges(7,3) * t411 - Ifges(5,5) * t412 - Ifges(5,6) * t413 - t62 * mrSges(5,1) + t63 * mrSges(5,2) + (-Ifges(4,2) * t283 + t385) * t327 + t258 * t200 + t441 + t442) * t287 + (Ifges(4,1) * t234 + Ifges(4,4) * t468 + t315 * t409 + t317 * t413 + t319 * t412) * t283 + (-t11 * t194 + t12 * t195 - t134 * t45 + t135 * t46) * mrSges(6,3) + t92 * (mrSges(6,1) * t194 - mrSges(6,2) * t195) + (-Ifges(6,4) * t195 - Ifges(6,2) * t194) * t425 + (-Ifges(6,1) * t195 - Ifges(6,4) * t194) * t426 + (-Ifges(6,5) * t195 - Ifges(6,6) * t194) * t410 + t126 * t436 + t125 * t437 + t134 * t423 + t135 * t424 + t47 * t427 + t48 * t429 - t195 * t431 - t194 * t432 + (Ifges(7,1) * t126 + Ifges(7,4) * t125) * t434 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t419 + (m(3) * t400 + mrSges(2,1) * t284 - t176 * mrSges(6,1) - t169 * mrSges(7,1) + mrSges(2,2) * t288 - t175 * mrSges(6,2) - t168 * mrSges(7,2) + t449 * (t267 * pkin(7) - t400) + t440 * t267 + (-m(6) * (-pkin(2) - t310) - m(7) * (-pkin(2) - t311) - m(5) * t329 + t283 * mrSges(5,3) + t459 + m(4) * pkin(2) + t445) * t265) * g(1) + (-t362 * t62 - t364 * t63) * mrSges(5,3) + t305 * t327 + t294 * t328 + t455 * t466 + t385 * t467 + t318 * t468 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t421 + (Ifges(7,4) * t126 + Ifges(7,2) * t125) * t433 + m(7) * (t13 * t5 + t14 * t4 + t143 * t44 + t2 * t29 + t28 * t3 + t89 * t94) + m(6) * (t11 * t67 + t12 * t66 + t141 * t193 + t202 * t92 + t34 * t46 + t35 * t45) - (Ifges(5,5) * t188 + Ifges(5,6) * t187 + Ifges(5,3) * t233 + t345 + t346) * t287 / 0.2e1 + (t125 * t2 - t126 * t3 - t13 * t47 + t14 * t48) * mrSges(7,3) + (t450 / 0.2e1 - t199 * mrSges(4,3) + Ifges(6,5) * t414 + Ifges(7,5) * t419 + Ifges(6,6) * t416 + Ifges(7,6) * t421 + Ifges(6,3) * t405 + Ifges(7,3) * t407 + t439) * t357 + (-t198 * t356 + t453) * mrSges(4,3) + t362 * t418 + t446 * t356 / 0.2e1 + (Ifges(6,5) * t134 + Ifges(6,6) * t135) * t405 + m(4) * (t240 * t260 + t453 * t258) + t456 * t236 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t407 + (Ifges(7,5) * t126 + Ifges(7,6) * t125) * t411 + m(5) * (t104 * t150 + t105 * t151 + t158 * t62 + t159 * t63 + (t180 * t356 + t373) * t258) + (-t201 + t115) * t246 + qJDD(3) * (Ifges(4,5) * t283 + Ifges(4,6) * t287) - t97 * t364 / 0.2e1 + t260 * (mrSges(4,1) * t233 + mrSges(4,2) * t234) - t147 * t338 / 0.2e1 - t244 * t339 + t202 * t39 + t321 * t373 + t193 * t90 + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t279 - 0.2e1 * mrSges(3,2) * t277 + m(3) * (t277 ^ 2 + t279 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t151 * t183 + t150 * t184 - t240 * t324 + t158 * t145 + t159 * t144 + t141 * (-mrSges(6,1) * t135 + mrSges(6,2) * t134) + t143 * t8 + t44 * (-mrSges(7,1) * t125 + mrSges(7,2) * t126) + t34 * t123 + t35 * t124 + t94 * t43 + t89 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + t66 * t60 + t67 * t61 + t4 * t68 + t5 * t69 + t28 * t17 + t29 * t18; m(3) * qJDD(2) + t134 * t123 + t135 * t124 + t125 * t17 + t126 * t18 - t194 * t60 - t195 * t61 + t47 * t68 + t48 * t69 + (t200 + t454) * t283 + (t201 + t451) * t287 + (-m(3) + t449) * g(3) + ((t183 * t278 - t184 * t276 + t244) * t287 + (t43 + t90 + t456) * t283) * qJD(3) + m(7) * (t125 * t3 + t126 * t2 + t13 * t48 + t14 * t47 - t287 * t44 + t357 * t89) + m(6) * (-t11 * t195 - t12 * t194 + t134 * t46 + t135 * t45 + t141 * t357 - t287 * t92) + m(4) * (t136 * t283 + t137 * t287 + (-t198 * t283 + t199 * t287) * qJD(3)) + m(5) * (-t127 * t287 + t313 * t283 + (t180 * t283 + t287 * t312) * qJD(3)); (-mrSges(6,1) * t474 - mrSges(6,2) * t479) * t141 + (-t11 * t309 - t12 * t225 + t45 * t479 + t46 * t474) * mrSges(6,3) - t479 * t423 + (t444 - t455 / 0.2e1 + (t294 / 0.2e1 - t305 / 0.2e1) * qJD(1)) * qJD(1) + t464 * t69 + (t13 * t464 + t14 * t465 + t190 * t44 + t2 * t65 + t3 * t64 + t452 * t89) * m(7) + t465 * t68 - t450 * t358 / 0.2e1 + t461 * t124 + (t11 * t166 + t12 * t165 - t141 * t167 - t259 * t92 + t45 * t461 + t46 * t462) * m(6) + t462 * t123 + (Ifges(7,5) * t108 + Ifges(7,6) * t107) * t408 + (t354 - t140) * t183 + t156 * t437 + t108 * t428 + t107 * t430 + t225 * t431 + (Ifges(7,4) * t157 + Ifges(7,2) * t156) * t433 + (Ifges(7,1) * t157 + Ifges(7,4) * t156) * t434 + t157 * t436 + t470 * t323 + (-t324 - t293) * g(3) + (-t107 * t14 + t108 * t13 + t156 * t2 - t157 * t3) * mrSges(7,3) + t316 * t328 + (Ifges(7,1) * t419 + Ifges(7,4) * t421 + Ifges(7,5) * t407 + t427 + t447) * (qJD(6) * t156 - t206 * t285 - t207 * t281) + (-pkin(3) * t127 + qJ(4) * t313 + qJD(4) * t312 - t104 * t139 - t105 * t140) * m(5) + (Ifges(7,1) * t108 + Ifges(7,4) * t107) * t420 + t474 * t424 + (t470 * (t302 - t473) + t472 * g(3)) * t283 + (t470 * (-t337 + t472) + t473 * g(3)) * t287 + (Ifges(7,4) * t419 + Ifges(7,2) * t421 + Ifges(7,6) * t407 + t429 + t448) * (-qJD(6) * t157 + t206 * t281 - t207 * t285) + (Ifges(5,5) * t276 + Ifges(5,6) * t278) * t409 + (Ifges(7,5) * t157 + Ifges(7,6) * t156) * t411 + (Ifges(5,1) * t276 + t383) * t412 + (Ifges(5,2) * t278 + t384) * t413 + t276 * t418 - (-Ifges(4,2) * t358 + t263 + t446) * t350 / 0.2e1 + (-t355 - t139) * t184 + t313 * mrSges(5,3) + t454 * qJ(4) + (-m(5) * t180 + t344 - t456) * t199 + (Ifges(6,5) * t415 + Ifges(7,5) * t420 + Ifges(6,6) * t417 + Ifges(7,6) * t422 + Ifges(6,3) * t406 + Ifges(7,3) * t408 - t439) * t358 + (Ifges(7,4) * t108 + Ifges(7,2) * t107) * t422 + t452 * t43 + t278 * t97 / 0.2e1 - t259 * t39 + t147 * t340 / 0.2e1 - Ifges(4,6) * t233 + Ifges(4,5) * t234 + t190 * t8 + (-Ifges(6,1) * t186 - Ifges(6,4) * t185) * t415 + (-Ifges(6,5) * t186 - Ifges(6,6) * t185) * t406 + (-Ifges(6,4) * t186 - Ifges(6,2) * t185) * t417 - t309 * t432 + (Ifges(6,4) * t225 - Ifges(6,2) * t309) * t425 + (Ifges(6,1) * t225 - Ifges(6,4) * t309) * t426 + (Ifges(6,5) * t225 - Ifges(6,6) * t309) * t410 + t165 * t60 + t166 * t61 + t92 * (mrSges(6,1) * t309 + mrSges(6,2) * t225) - t167 * t90 + t44 * (-mrSges(7,1) * t156 + mrSges(7,2) * t157) - t136 * mrSges(4,2) + t137 * mrSges(4,1) - t89 * (-mrSges(7,1) * t107 + mrSges(7,2) * t108) - pkin(3) * t115 + t64 * t17 + t65 * t18 + Ifges(4,3) * qJDD(3) + t127 * t322 + (-t244 + t343) * t198 + (-Ifges(6,5) * t206 - Ifges(6,6) * t207) * t405 + (-Ifges(6,1) * t206 - Ifges(6,4) * t207) * t414 + (-Ifges(6,4) * t206 - Ifges(6,2) * t207) * t416; -t325 * t123 + t155 * t124 - t218 * t183 + t219 * t184 - t469 * t68 + t85 * t69 + (t13 * t85 - t14 * t469 + t44) * m(7) + (t155 * t45 - t325 * t46 + t92) * m(6) + (t104 * t219 - t105 * t218 + t127) * m(5) + (t287 * g(3) - t283 * t470) * t349 - t451; (t387 + t124) * t46 + t345 - (Ifges(7,4) * t420 + Ifges(7,2) * t422 + Ifges(7,6) * t408 + t430 - t448) * t85 + (Ifges(7,1) * t420 + Ifges(7,4) * t422 + Ifges(7,5) * t408 + t428 - t447) * t469 + (t2 * t281 + t285 * t3 + (-t13 * t281 + t14 * t285) * qJD(6)) * t435 + (m(7) * t396 + mrSges(6,1) * t264 + mrSges(6,2) * t266 - t320) * t393 + (t388 - t123) * t45 + (Ifges(6,5) * t325 - Ifges(6,6) * t155) * t406 + (Ifges(6,1) * t325 - t382) * t415 - t141 * (mrSges(6,1) * t155 + mrSges(6,2) * t325) + (-Ifges(6,2) * t155 + t149 + t77) * t417 + t307 + t76 * t414 - t43 * t398 - m(7) * (t13 * t15 + t14 * t16 + t398 * t89) + (mrSges(6,2) * t178 - t177 * t463 - t359) * g(1) + (-mrSges(6,2) * t176 + t175 * t463 - t360) * g(2) - t441 - t16 * t68 - t15 * t69 + ((-t281 * t69 + t285 * t68) * qJD(6) + t17 * t285 + t18 * t281) * pkin(5); -t89 * (mrSges(7,1) * t85 + mrSges(7,2) * t469) + (Ifges(7,1) * t469 - t402) * t420 + t41 * t419 + (Ifges(7,5) * t469 - Ifges(7,6) * t85) * t408 - t13 * t68 + t14 * t69 - g(1) * t359 - g(2) * t360 - t320 * t393 + (t13 * t469 + t14 * t85) * mrSges(7,3) + t307 + (-Ifges(7,2) * t85 + t42 + t78) * t422;];
tau  = t1;
