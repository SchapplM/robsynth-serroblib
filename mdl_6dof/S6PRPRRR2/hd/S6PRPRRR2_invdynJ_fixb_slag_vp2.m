% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:46
% EndTime: 2019-03-08 20:27:19
% DurationCPUTime: 18.39s
% Computational Cost: add. (8179->674), mult. (18623->943), div. (0->0), fcn. (14536->16), ass. (0->308)
t229 = sin(pkin(12));
t237 = sin(qJ(2));
t241 = cos(qJ(2));
t358 = cos(pkin(12));
t186 = -t241 * t229 - t237 * t358;
t231 = sin(pkin(6));
t159 = t186 * t231;
t148 = qJD(1) * t159;
t236 = sin(qJ(4));
t240 = cos(qJ(4));
t288 = pkin(4) * t236 - pkin(9) * t240;
t459 = t288 * qJD(4) + t148;
t285 = mrSges(5,1) * t240 - mrSges(5,2) * t236;
t239 = cos(qJ(5));
t224 = pkin(5) * t239 + pkin(4);
t228 = qJ(5) + qJ(6);
t226 = sin(t228);
t227 = cos(t228);
t235 = sin(qJ(5));
t284 = -mrSges(6,1) * t239 + mrSges(6,2) * t235;
t430 = m(6) * pkin(4) + m(7) * t224 + mrSges(7,1) * t227 - mrSges(7,2) * t226 - t284;
t242 = -pkin(10) - pkin(9);
t449 = -m(6) * pkin(9) + m(7) * t242 - mrSges(6,3) - mrSges(7,3);
t458 = -t449 * t236 + t240 * t430 + t285;
t334 = qJD(1) * t231;
t307 = t237 * t334;
t203 = t229 * t307;
t289 = t358 * t334;
t151 = t241 * t289 - t203;
t268 = -t240 * pkin(4) - t236 * pkin(9) - pkin(3);
t308 = t358 * pkin(2);
t181 = -t308 + t268;
t221 = pkin(2) * t229 + pkin(8);
t324 = qJD(5) * t239;
t326 = qJD(5) * t235;
t328 = qJD(4) * t239;
t339 = t239 * t240;
t436 = -t151 * t339 + t181 * t324 + (-t236 * t328 - t240 * t326) * t221 + t459 * t235;
t329 = qJD(4) * t236;
t305 = t221 * t329;
t343 = t235 * t240;
t457 = t151 * t343 + t235 * t305 + t239 * t459;
t455 = m(5) + m(6);
t192 = t221 * t339;
t267 = pkin(5) * t236 - pkin(10) * t339;
t454 = t267 * qJD(4) + (-t192 + (pkin(10) * t236 - t181) * t235) * qJD(5) + t457;
t327 = qJD(4) * t240;
t252 = t235 * t327 + t236 * t324;
t453 = -pkin(10) * t252 + t436;
t331 = qJD(2) * t240;
t306 = t235 * t331;
t309 = qJD(5) * t242;
t333 = qJD(1) * t241;
t201 = qJD(2) * pkin(2) + t231 * t333;
t140 = t229 * t201 + t237 * t289;
t136 = qJD(2) * pkin(8) + t140;
t233 = cos(pkin(6));
t216 = qJD(1) * t233 + qJD(3);
t100 = -t236 * t136 + t216 * t240;
t193 = t288 * qJD(2);
t65 = t239 * t100 + t235 * t193;
t452 = pkin(10) * t306 + t235 * t309 - t65;
t64 = -t100 * t235 + t239 * t193;
t451 = -qJD(2) * t267 + t239 * t309 - t64;
t323 = qJD(2) * qJD(4);
t199 = qJDD(2) * t240 - t236 * t323;
t184 = qJDD(5) - t199;
t180 = qJDD(6) + t184;
t392 = t180 / 0.2e1;
t332 = qJD(2) * t236;
t187 = -t235 * t332 + t328;
t200 = qJDD(2) * t236 + t240 * t323;
t119 = qJD(5) * t187 + qJDD(4) * t235 + t200 * t239;
t330 = qJD(4) * t235;
t188 = t239 * t332 + t330;
t120 = -qJD(5) * t188 + qJDD(4) * t239 - t200 * t235;
t234 = sin(qJ(6));
t238 = cos(qJ(6));
t126 = t187 * t234 + t188 * t238;
t37 = -qJD(6) * t126 - t119 * t234 + t120 * t238;
t404 = t37 / 0.2e1;
t292 = t238 * t187 - t188 * t234;
t36 = qJD(6) * t292 + t119 * t238 + t120 * t234;
t405 = t36 / 0.2e1;
t406 = Ifges(7,1) * t405 + Ifges(7,4) * t404 + Ifges(7,5) * t392;
t407 = Ifges(7,4) * t405 + Ifges(7,2) * t404 + Ifges(7,6) * t392;
t260 = -t237 * t229 + t241 * t358;
t158 = t260 * t231;
t283 = t235 * mrSges(6,1) + t239 * mrSges(6,2);
t384 = pkin(5) * t235;
t448 = -m(7) * (pkin(8) + t384) - t226 * mrSges(7,1) - t227 * mrSges(7,2) - mrSges(5,3) - t283 - t455 * pkin(8);
t273 = -t201 * t358 + t203;
t106 = qJD(2) * t268 + t273;
t215 = qJDD(1) * t233 + qJDD(3);
t291 = qJD(2) * t307;
t348 = t231 * t241;
t165 = qJDD(1) * t348 - t291;
t357 = qJDD(2) * pkin(2);
t154 = t165 + t357;
t301 = qJD(2) * t333;
t166 = (qJDD(1) * t237 + t301) * t231;
t96 = t229 * t154 + t358 * t166;
t91 = qJDD(2) * pkin(8) + t96;
t38 = -t136 * t329 + t236 * t215 + t216 * t327 + t240 * t91;
t34 = qJDD(4) * pkin(9) + t38;
t95 = t154 * t358 - t229 * t166;
t90 = -qJDD(2) * pkin(3) - t95;
t60 = -t199 * pkin(4) - t200 * pkin(9) + t90;
t354 = t216 * t236;
t101 = t136 * t240 + t354;
t94 = qJD(4) * pkin(9) + t101;
t8 = t106 * t324 + t235 * t60 + t239 * t34 - t326 * t94;
t47 = t106 * t235 + t239 * t94;
t9 = -qJD(5) * t47 - t235 * t34 + t239 * t60;
t287 = -t235 * t9 + t239 * t8;
t46 = t239 * t106 - t235 * t94;
t447 = -t46 * t324 - t47 * t326 + t287;
t230 = sin(pkin(11));
t232 = cos(pkin(11));
t345 = t233 * t241;
t446 = -t230 * t237 + t232 * t345;
t316 = mrSges(5,3) * t332;
t431 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t187 + mrSges(6,2) * t188 + t316;
t63 = -mrSges(7,1) * t292 + mrSges(7,2) * t126;
t319 = t63 + t431;
t93 = -qJD(4) * pkin(4) - t100;
t444 = m(6) * t93;
t70 = -pkin(5) * t187 + t93;
t445 = -m(7) * t70 - t319 - t444;
t398 = t119 / 0.2e1;
t397 = t120 / 0.2e1;
t391 = t184 / 0.2e1;
t443 = t199 / 0.2e1;
t442 = t200 / 0.2e1;
t164 = t239 * t181;
t342 = t236 * t239;
t105 = -pkin(10) * t342 + t164 + (-t221 * t235 - pkin(5)) * t240;
t129 = t235 * t181 + t192;
t344 = t235 * t236;
t116 = -pkin(10) * t344 + t129;
t50 = t105 * t238 - t116 * t234;
t441 = qJD(6) * t50 + t454 * t234 + t453 * t238;
t51 = t105 * t234 + t116 * t238;
t440 = -qJD(6) * t51 - t453 * t234 + t454 * t238;
t211 = t242 * t235;
t212 = t242 * t239;
t141 = t211 * t238 + t212 * t234;
t439 = qJD(6) * t141 + t451 * t234 + t452 * t238;
t142 = t211 * t234 - t212 * t238;
t438 = -qJD(6) * t142 - t452 * t234 + t451 * t238;
t437 = -qJD(5) * t129 + t457;
t408 = m(7) * pkin(5);
t435 = -mrSges(6,1) - t408;
t218 = qJD(5) - t331;
t213 = qJD(6) + t218;
t434 = t188 * Ifges(6,5) + t126 * Ifges(7,5) + t187 * Ifges(6,6) + Ifges(7,6) * t292 + t218 * Ifges(6,3) + t213 * Ifges(7,3);
t61 = -mrSges(6,1) * t120 + mrSges(6,2) * t119;
t433 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t200 + t61;
t432 = -t354 - (qJD(2) * t384 + t136) * t240 + pkin(5) * t326;
t270 = t234 * t235 - t238 * t239;
t162 = t270 * t236;
t429 = (mrSges(3,1) * t241 - mrSges(3,2) * t237) * t231 + t158 * mrSges(4,1) + t159 * mrSges(4,2);
t182 = Ifges(6,4) * t187;
t115 = t188 * Ifges(6,1) + t218 * Ifges(6,5) + t182;
t225 = Ifges(5,4) * t331;
t426 = Ifges(5,1) * t332 + Ifges(5,5) * qJD(4) + t239 * t115 + t225;
t85 = mrSges(6,1) * t184 - mrSges(6,3) * t119;
t86 = -mrSges(6,2) * t184 + mrSges(6,3) * t120;
t425 = -t235 * t85 + t239 * t86;
t39 = -t136 * t327 + t215 * t240 - t216 * t329 - t236 * t91;
t424 = -t236 * t39 + t240 * t38;
t423 = -m(7) - t455;
t422 = qJD(5) + qJD(6);
t421 = mrSges(5,1) + t430;
t420 = mrSges(5,2) + t449;
t40 = -pkin(10) * t188 + t46;
t28 = pkin(5) * t218 + t40;
t41 = pkin(10) * t187 + t47;
t361 = t238 * t41;
t15 = t234 * t28 + t361;
t419 = -t70 * mrSges(7,1) + t15 * mrSges(7,3);
t365 = t234 * t41;
t14 = t238 * t28 - t365;
t418 = t70 * mrSges(7,2) - t14 * mrSges(7,3);
t335 = t186 * t233;
t112 = t230 * t335 + t232 * t260;
t417 = -t230 * t260 + t232 * t335;
t416 = -mrSges(4,1) - t458;
t413 = -t9 * mrSges(6,1) + t8 * mrSges(6,2);
t6 = pkin(5) * t184 - pkin(10) * t119 + t9;
t7 = pkin(10) * t120 + t8;
t2 = qJD(6) * t14 + t234 * t6 + t238 * t7;
t3 = -qJD(6) * t15 - t234 * t7 + t238 * t6;
t412 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t411 = -mrSges(4,2) - t448;
t372 = Ifges(5,4) * t236;
t279 = t240 * Ifges(5,2) + t372;
t410 = t15 * mrSges(7,2) + Ifges(5,6) * qJD(4) / 0.2e1 + qJD(2) * t279 / 0.2e1 - t14 * mrSges(7,1);
t409 = qJD(2) ^ 2;
t403 = Ifges(6,1) * t398 + Ifges(6,4) * t397 + Ifges(6,5) * t391;
t367 = Ifges(7,4) * t126;
t54 = Ifges(7,2) * t292 + Ifges(7,6) * t213 + t367;
t402 = -t54 / 0.2e1;
t401 = t54 / 0.2e1;
t121 = Ifges(7,4) * t292;
t55 = Ifges(7,1) * t126 + Ifges(7,5) * t213 + t121;
t400 = -t55 / 0.2e1;
t399 = t55 / 0.2e1;
t396 = -t292 / 0.2e1;
t395 = t292 / 0.2e1;
t394 = -t126 / 0.2e1;
t393 = t126 / 0.2e1;
t389 = t188 / 0.2e1;
t388 = -t213 / 0.2e1;
t387 = t213 / 0.2e1;
t385 = pkin(5) * t188;
t251 = t233 * t260;
t108 = t230 * t186 + t232 * t251;
t350 = t231 * t236;
t82 = -t232 * t350 - t240 * t417;
t377 = (-t108 * t227 - t226 * t82) * mrSges(7,1) + (t108 * t226 - t227 * t82) * mrSges(7,2);
t111 = t186 * t232 - t230 * t251;
t84 = t112 * t240 + t230 * t350;
t376 = (-t111 * t227 - t226 * t84) * mrSges(7,1) + (t111 * t226 - t227 * t84) * mrSges(7,2);
t134 = -t159 * t240 + t233 * t236;
t375 = (-t134 * t226 - t158 * t227) * mrSges(7,1) + (-t134 * t227 + t158 * t226) * mrSges(7,2);
t374 = mrSges(6,3) * t187;
t373 = mrSges(6,3) * t188;
t371 = Ifges(5,4) * t240;
t370 = Ifges(6,4) * t188;
t369 = Ifges(6,4) * t235;
t368 = Ifges(6,4) * t239;
t35 = -qJDD(4) * pkin(4) - t39;
t363 = t236 * t35;
t349 = t231 * t240;
t346 = t233 * t237;
t325 = qJD(5) * t236;
t217 = pkin(2) * t348;
t322 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t180;
t16 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t321 = t16 + t433;
t318 = m(4) - t423;
t315 = mrSges(5,3) * t331;
t311 = Ifges(6,5) * t119 + Ifges(6,6) * t120 + Ifges(6,3) * t184;
t304 = t221 * t327;
t114 = t187 * Ifges(6,2) + t218 * Ifges(6,6) + t370;
t303 = -t235 * t114 / 0.2e1;
t293 = t323 / 0.2e1;
t290 = t446 * pkin(2);
t281 = Ifges(6,1) * t239 - t369;
t280 = Ifges(6,1) * t235 + t368;
t278 = -Ifges(6,2) * t235 + t368;
t277 = Ifges(6,2) * t239 + t369;
t276 = Ifges(5,5) * t240 - Ifges(5,6) * t236;
t275 = Ifges(6,5) * t239 - Ifges(6,6) * t235;
t274 = Ifges(6,5) * t235 + Ifges(6,6) * t239;
t77 = -t134 * t235 - t158 * t239;
t78 = t134 * t239 - t158 * t235;
t26 = -t234 * t78 + t238 * t77;
t27 = t234 * t77 + t238 * t78;
t133 = -t159 * t236 - t233 * t240;
t190 = t234 * t239 + t235 * t238;
t266 = t322 - t412;
t265 = -t230 * t345 - t232 * t237;
t264 = t93 * t283;
t135 = -qJD(2) * pkin(3) + t273;
t263 = t135 * (mrSges(5,1) * t236 + mrSges(5,2) * t240);
t262 = t236 * (Ifges(5,1) * t240 - t372);
t259 = t190 * t240;
t258 = t270 * t240;
t254 = t265 * pkin(2);
t253 = -t235 * t325 + t239 * t327;
t249 = Ifges(6,5) * t236 + t240 * t281;
t248 = Ifges(6,6) * t236 + t240 * t278;
t247 = Ifges(6,3) * t236 + t240 * t275;
t131 = t422 * t190;
t222 = -t308 - pkin(3);
t206 = -qJD(4) * mrSges(5,2) + t315;
t191 = t285 * qJD(2);
t169 = (t221 + t384) * t236;
t167 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t199;
t161 = t190 * t236;
t157 = qJD(2) * t258;
t156 = qJD(2) * t259;
t150 = qJD(2) * t158;
t149 = qJD(2) * t159;
t145 = mrSges(6,1) * t218 - t373;
t144 = -mrSges(6,2) * t218 + t374;
t138 = -mrSges(5,1) * t199 + mrSges(5,2) * t200;
t137 = pkin(5) * t252 + t304;
t128 = -t221 * t343 + t164;
t98 = mrSges(7,1) * t213 - mrSges(7,3) * t126;
t97 = -mrSges(7,2) * t213 + mrSges(7,3) * t292;
t76 = -qJD(4) * t133 + t150 * t240;
t75 = qJD(4) * t134 + t150 * t236;
t72 = -qJD(4) * t259 + t162 * t422;
t71 = -qJD(4) * t258 - t131 * t236;
t48 = t119 * Ifges(6,4) + t120 * Ifges(6,2) + t184 * Ifges(6,6);
t25 = -mrSges(7,2) * t180 + mrSges(7,3) * t37;
t24 = mrSges(7,1) * t180 - mrSges(7,3) * t36;
t23 = qJD(5) * t77 - t149 * t235 + t239 * t76;
t22 = -qJD(5) * t78 - t149 * t239 - t235 * t76;
t21 = -pkin(5) * t120 + t35;
t18 = t238 * t40 - t365;
t17 = -t234 * t40 - t361;
t5 = -qJD(6) * t27 + t22 * t238 - t23 * t234;
t4 = qJD(6) * t26 + t22 * t234 + t23 * t238;
t1 = [m(2) * qJDD(1) + t134 * t167 - t158 * t138 + t23 * t144 + t22 * t145 + t149 * t191 + t76 * t206 + t26 * t24 + t27 * t25 + t4 * t97 + t5 * t98 + t77 * t85 + t78 * t86 + (-mrSges(3,1) * t237 - mrSges(3,2) * t241) * t409 * t231 + (mrSges(4,1) * t149 - mrSges(4,2) * t150) * qJD(2) + t319 * t75 + t321 * t133 + t429 * qJDD(2) + (-m(2) - m(3) - t318) * g(3) + m(3) * (qJDD(1) * t233 ^ 2 + (t165 * t241 + t166 * t237) * t231) + m(4) * (t140 * t150 - t149 * t273 + t158 * t95 - t159 * t96 + t215 * t233) + m(5) * (-t100 * t75 + t101 * t76 - t133 * t39 + t134 * t38 - t135 * t149 - t158 * t90) + m(7) * (t133 * t21 + t14 * t5 + t15 * t4 + t2 * t27 + t26 * t3 + t70 * t75) + m(6) * (t133 * t35 + t22 * t46 + t23 * t47 + t75 * t93 + t77 * t9 + t78 * t8); (qJD(2) * t151 - t229 * t357 - t96) * mrSges(4,2) + (Ifges(7,5) * t71 + Ifges(7,6) * t72) * t387 + (t231 * t301 - t166) * mrSges(3,2) + (Ifges(7,4) * t71 + Ifges(7,2) * t72) * t395 + (Ifges(7,1) * t71 + Ifges(7,4) * t72) * t393 + (t273 * t148 - t140 * t151 + (t229 * t96 + t358 * t95) * pkin(2)) * m(4) + (-qJD(2) * t148 + qJDD(2) * t308 + t95) * mrSges(4,1) + (t135 * t148 - (-t100 * t236 + t101 * t240) * t151 + t222 * t90 + ((-t100 * t240 - t101 * t236) * qJD(4) + t424) * t221) * m(5) + (-t252 * t47 - t253 * t46 - t342 * t9 - t344 * t8) * mrSges(6,3) + (t46 * mrSges(6,1) - t47 * mrSges(6,2) + Ifges(7,5) * t393 + Ifges(7,6) * t395 + Ifges(7,3) * t387 - t410 + t434 / 0.2e1 - t101 * mrSges(5,3)) * t329 - (t239 * t114 + t235 * t115) * t325 / 0.2e1 - (t322 + t311) * t240 / 0.2e1 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (-t100 * t327 + t424) * mrSges(5,3) + t433 * t221 * t236 + t431 * t304 + t426 * t327 / 0.2e1 + (-t265 * mrSges(3,1) - (t230 * t346 - t232 * t241) * mrSges(3,2) - m(4) * t254 + t423 * (t111 * pkin(3) + t254) + t416 * t111 - t411 * t112) * g(1) + t21 * (mrSges(7,1) * t161 - mrSges(7,2) * t162) + (-Ifges(7,5) * t162 - Ifges(7,6) * t161) * t392 + (-Ifges(7,4) * t162 - Ifges(7,2) * t161) * t404 + (-Ifges(7,1) * t162 - Ifges(7,4) * t161) * t405 + (-t14 * t71 + t15 * t72 - t161 * t2 + t162 * t3) * mrSges(7,3) + t371 * t442 + (t291 + t165) * mrSges(3,1) + t279 * t443 + (-m(4) * t217 + t423 * (t158 * pkin(3) + t217) - t448 * t159 - t458 * t158 - t429) * g(3) + (Ifges(5,4) * t442 + Ifges(5,2) * t443 + t221 * t167 - t151 * t206 - Ifges(6,6) * t397 - Ifges(6,5) * t398 - Ifges(7,6) * t404 - Ifges(7,5) * t405 - Ifges(6,3) * t391 - Ifges(7,3) * t392 + (-Ifges(5,2) * t236 + t371) * t293 + t412 + t413) * t240 + (Ifges(5,1) * t200 + Ifges(5,4) * t443 + t275 * t391 + t278 * t397 + t281 * t398) * t236 + t440 * t98 + t441 * t97 + (t137 * t70 + t14 * t440 + t15 * t441 + t169 * t21 + t2 * t51 + t3 * t50) * m(7) + t436 * t144 + t437 * t145 + (t128 * t9 + t129 * t8 + (t327 * t93 + t363) * t221 + t436 * t47 + t437 * t46) * m(6) + (-t446 * mrSges(3,1) - (-t230 * t241 - t232 * t346) * mrSges(3,2) - m(4) * t290 + t423 * (t108 * pkin(3) + t290) + t416 * t108 + t411 * t417) * g(2) - t48 * t344 / 0.2e1 + t187 * (qJD(4) * t248 - t277 * t325) / 0.2e1 + t218 * (qJD(4) * t247 - t274 * t325) / 0.2e1 + t70 * (-mrSges(7,1) * t72 + mrSges(7,2) * t71) + t50 * t24 + t51 * t25 - t206 * t305 + t445 * t151 * t236 - t90 * t285 + qJD(4) ^ 2 * t276 / 0.2e1 + t71 * t399 + t72 * t401 + t342 * t403 - t162 * t406 - t161 * t407 + (qJD(4) * t249 - t280 * t325) * t389 + t283 * t363 + t303 * t327 + qJD(4) * t263 + t93 * (mrSges(6,1) * t252 + mrSges(6,2) * t253) + qJDD(4) * (Ifges(5,5) * t236 + Ifges(5,6) * t240) + t128 * t85 + t129 * t86 + t137 * t63 + t169 * t16 - t148 * t191 + t262 * t293 + t222 * t138; -t161 * t24 - t162 * t25 + t71 * t97 + t72 * t98 + m(7) * (t14 * t72 + t15 * t71 - t161 * t3 - t162 * t2) + m(4) * t215 + ((t144 * t239 - t145 * t235 + t206) * qJD(4) + m(5) * (qJD(4) * t101 + t39) + m(6) * (t328 * t47 - t330 * t46 - t35) - m(7) * t21 - t321) * t240 + (t167 + (-t235 * t144 - t239 * t145) * qJD(5) + m(5) * t38 + m(6) * t447 + (-m(5) * t100 - t445) * qJD(4) + t425) * t236 + (-t233 * g(3) + (-g(1) * t230 + g(2) * t232) * t231) * t318; (-Ifges(7,5) * t157 - Ifges(7,6) * t156) * t388 - t70 * (mrSges(7,1) * t156 - mrSges(7,2) * t157) + (-Ifges(7,4) * t157 - Ifges(7,2) * t156) * t396 + (-Ifges(7,1) * t157 - Ifges(7,4) * t156) * t394 + (-pkin(4) * t35 - t46 * t64 - t47 * t65) * m(6) + (-t14 * t157 + t15 * t156) * mrSges(7,3) + (t264 + t303) * qJD(5) + (t420 * t82 - t421 * (-t232 * t349 + t236 * t417)) * g(2) + (t187 * t278 + t188 * t281 + t218 * t275) * qJD(5) / 0.2e1 - (t187 * t248 + t188 * t249 + t218 * t247) * qJD(2) / 0.2e1 - (Ifges(7,4) * t393 + Ifges(7,2) * t395 + Ifges(7,6) * t387 + t401 + t419) * t131 + (Ifges(7,5) * t394 + Ifges(7,6) * t396 + Ifges(7,3) * t388 + t410) * t332 + (-t2 * mrSges(7,3) + t21 * mrSges(7,1) - 0.2e1 * t407 - (Ifges(7,1) * t393 + Ifges(7,4) * t395 + Ifges(7,5) * t387 + t399 + t418) * t422) * t270 + (mrSges(7,2) * t21 - mrSges(7,3) * t3 + 0.2e1 * t406) * t190 + t432 * t63 - t434 * t332 / 0.2e1 + (m(6) * ((-t235 * t47 - t239 * t46) * qJD(5) + t287) - t145 * t324 - t144 * t326 + t425) * pkin(9) - (-Ifges(5,2) * t332 + t225 + t426) * t331 / 0.2e1 + (t133 * t421 + t134 * t420) * g(3) + (t420 * t84 - t421 * (-t112 * t236 + t230 * t349)) * g(1) + (t316 - t431 - t444) * t101 + t439 * t97 + (t14 * t438 + t141 * t3 + t142 * t2 + t15 * t439 - t21 * t224 + t432 * t70) * m(7) + t438 * t98 + t447 * mrSges(6,3) - t264 * t331 - t276 * t323 / 0.2e1 + t115 * t324 / 0.2e1 - pkin(4) * t61 - t38 * mrSges(5,2) + t39 * mrSges(5,1) + t114 * t306 / 0.2e1 + (-t263 - t46 * (mrSges(6,1) * t236 - mrSges(6,3) * t339) - t47 * (-mrSges(6,2) * t236 - mrSges(6,3) * t343)) * qJD(2) + Ifges(5,3) * qJDD(4) + t35 * t284 + t277 * t397 + t280 * t398 - t157 * t400 - t156 * t402 + t235 * t403 + t274 * t391 + t141 * t24 + t142 * t25 - t65 * t144 - t64 * t145 + Ifges(5,6) * t199 + Ifges(5,5) * t200 + (t315 - t206) * t100 - t224 * t16 - t409 * t262 / 0.2e1 + t239 * t48 / 0.2e1; -t413 - (Ifges(7,4) * t394 + Ifges(7,2) * t396 + Ifges(7,6) * t388 + t402 - t419) * t126 + (Ifges(7,1) * t394 + Ifges(7,4) * t396 + Ifges(7,5) * t388 + t400 - t418) * t292 + (-(t111 * t235 - t239 * t84) * mrSges(6,2) - t376 + t435 * (-t111 * t239 - t235 * t84)) * g(1) + (-(t108 * t235 - t239 * t82) * mrSges(6,2) - t377 + t435 * (-t108 * t239 - t235 * t82)) * g(2) - (-Ifges(6,2) * t188 + t115 + t182) * t187 / 0.2e1 + (mrSges(6,2) * t78 + t435 * t77 - t375) * g(3) - t188 * (Ifges(6,1) * t187 - t370) / 0.2e1 - t17 * t98 - t18 * t97 + (t373 + t145) * t47 + t311 + t266 + (t2 * t234 + t238 * t3 + (-t14 * t234 + t15 * t238) * qJD(6)) * t408 + t114 * t389 - t93 * (mrSges(6,1) * t188 + mrSges(6,2) * t187) - t218 * (Ifges(6,5) * t187 - Ifges(6,6) * t188) / 0.2e1 - t63 * t385 - m(7) * (t14 * t17 + t15 * t18 + t385 * t70) + (t374 - t144) * t46 + ((-t234 * t98 + t238 * t97) * qJD(6) + t234 * t25 + t238 * t24) * pkin(5); -t70 * (mrSges(7,1) * t126 + mrSges(7,2) * t292) + (Ifges(7,1) * t292 - t367) * t394 + t54 * t393 + (Ifges(7,5) * t292 - Ifges(7,6) * t126) * t388 - t14 * t97 + t15 * t98 - g(1) * t376 - g(2) * t377 - g(3) * t375 + (t126 * t15 + t14 * t292) * mrSges(7,3) + t266 + (-Ifges(7,2) * t126 + t121 + t55) * t396;];
tau  = t1;
