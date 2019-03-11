% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:27
% EndTime: 2019-03-09 08:26:07
% DurationCPUTime: 25.14s
% Computational Cost: add. (11201->726), mult. (25632->897), div. (0->0), fcn. (19092->14), ass. (0->324)
t449 = mrSges(6,1) + mrSges(7,1);
t448 = mrSges(6,2) - mrSges(7,3);
t257 = sin(pkin(9));
t262 = sin(qJ(2));
t265 = cos(qJ(2));
t349 = cos(pkin(9));
t219 = t257 * t265 + t262 * t349;
t201 = t219 * qJD(1);
t256 = sin(pkin(10));
t258 = cos(pkin(10));
t179 = qJD(2) * t256 + t201 * t258;
t180 = qJD(2) * t258 - t201 * t256;
t261 = sin(qJ(5));
t264 = cos(qJ(5));
t113 = t179 * t261 - t180 * t264;
t320 = qJD(1) * qJD(2);
t307 = t262 * t320;
t319 = qJDD(1) * t265;
t224 = -t307 + t319;
t225 = qJDD(1) * t262 + t265 * t320;
t172 = t257 * t224 + t225 * t349;
t140 = qJDD(2) * t258 - t172 * t256;
t141 = qJDD(2) * t256 + t172 * t258;
t47 = -qJD(5) * t113 + t261 * t140 + t264 * t141;
t467 = -t47 / 0.2e1;
t456 = t264 * t179 + t180 * t261;
t48 = qJD(5) * t456 - t264 * t140 + t261 * t141;
t408 = -t48 / 0.2e1;
t171 = -t349 * t224 + t225 * t257;
t168 = qJDD(5) + t171;
t466 = -t168 / 0.2e1;
t453 = -m(6) - m(7);
t465 = -mrSges(6,3) - mrSges(7,2);
t447 = Ifges(6,1) + Ifges(7,1);
t446 = -Ifges(6,4) + Ifges(7,5);
t445 = Ifges(7,4) + Ifges(6,5);
t444 = -Ifges(6,6) + Ifges(7,6);
t443 = Ifges(6,3) + Ifges(7,2);
t303 = t349 * t265;
t328 = qJD(1) * t262;
t199 = -qJD(1) * t303 + t257 * t328;
t220 = t256 * t264 + t258 * t261;
t131 = t220 * t199;
t204 = t220 * qJD(5);
t433 = t131 + t204;
t218 = t256 * t261 - t264 * t258;
t132 = t218 * t199;
t203 = t218 * qJD(5);
t432 = t132 + t203;
t233 = -mrSges(3,1) * t265 + mrSges(3,2) * t262;
t255 = qJ(2) + pkin(9);
t252 = cos(t255);
t294 = -t258 * mrSges(5,1) + t256 * mrSges(5,2);
t274 = m(5) * pkin(3) - t294;
t464 = -t274 * t252 + t233;
t254 = pkin(10) + qJ(5);
t249 = sin(t254);
t251 = cos(t254);
t463 = -t448 * t249 + t251 * t449;
t347 = qJDD(1) * pkin(1);
t192 = -pkin(2) * t224 + qJDD(3) - t347;
t80 = pkin(3) * t171 - qJ(4) * t172 - qJD(4) * t201 + t192;
t216 = t225 * pkin(7);
t325 = qJD(3) * t262;
t158 = qJDD(2) * pkin(2) - qJ(3) * t225 - qJD(1) * t325 - t216;
t247 = pkin(7) * t319;
t326 = qJD(2) * t262;
t313 = pkin(7) * t326;
t324 = qJD(3) * t265;
t166 = qJ(3) * t224 + t247 + (-t313 + t324) * qJD(1);
t99 = t257 * t158 + t349 * t166;
t90 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t99;
t33 = -t256 * t90 + t258 * t80;
t20 = pkin(4) * t171 - pkin(8) * t141 + t33;
t253 = t265 * pkin(2);
t246 = t253 + pkin(1);
t229 = -qJD(1) * t246 + qJD(3);
t125 = pkin(3) * t199 - qJ(4) * t201 + t229;
t259 = -qJ(3) - pkin(7);
t232 = t259 * t262;
t222 = qJD(1) * t232;
t213 = qJD(2) * pkin(2) + t222;
t234 = t259 * t265;
t223 = qJD(1) * t234;
t304 = t349 * t223;
t163 = t257 * t213 - t304;
t151 = qJD(2) * qJ(4) + t163;
t82 = t258 * t125 - t151 * t256;
t59 = pkin(4) * t199 - pkin(8) * t179 + t82;
t83 = t256 * t125 + t258 * t151;
t69 = pkin(8) * t180 + t83;
t22 = t261 * t59 + t264 * t69;
t34 = t256 * t80 + t258 * t90;
t26 = pkin(8) * t140 + t34;
t4 = -qJD(5) * t22 + t20 * t264 - t26 * t261;
t2 = -pkin(5) * t168 + qJDD(6) - t4;
t394 = t168 / 0.2e1;
t407 = t48 / 0.2e1;
t409 = t47 / 0.2e1;
t98 = t158 * t349 - t257 * t166;
t93 = -qJDD(2) * pkin(3) + qJDD(4) - t98;
t63 = -t140 * pkin(4) + t93;
t7 = t48 * pkin(5) - t47 * qJ(6) - qJD(6) * t456 + t63;
t462 = -mrSges(6,2) * t63 - mrSges(7,2) * t2 + mrSges(6,3) * t4 + mrSges(7,3) * t7 - Ifges(7,5) * t407 + (-t409 + t467) * t447 + (-t394 + t466) * t445 + (-Ifges(6,4) + t446) * t408;
t461 = m(5) * qJ(4) + mrSges(5,3);
t112 = Ifges(6,4) * t113;
t196 = qJD(5) + t199;
t356 = Ifges(7,5) * t113;
t439 = t445 * t196 + t447 * t456 - t112 + t356;
t460 = -t439 / 0.2e1;
t374 = t258 * pkin(4);
t244 = pkin(3) + t374;
t250 = sin(t255);
t260 = -pkin(8) - qJ(4);
t459 = t252 * t244 - t250 * t260;
t263 = sin(qJ(1));
t266 = cos(qJ(1));
t458 = g(1) * t266 + g(2) * t263;
t457 = -m(3) * pkin(1) - mrSges(2,1) + t464;
t205 = t257 * t223;
t162 = t213 * t349 + t205;
t144 = -qJD(2) * pkin(3) + qJD(4) - t162;
t107 = -pkin(4) * t180 + t144;
t35 = t113 * pkin(5) - qJ(6) * t456 + t107;
t111 = Ifges(7,5) * t456;
t53 = Ifges(7,6) * t196 + Ifges(7,3) * t113 + t111;
t357 = Ifges(6,4) * t456;
t56 = -Ifges(6,2) * t113 + Ifges(6,6) * t196 + t357;
t455 = -mrSges(7,1) * t35 - t53 / 0.2e1 + t56 / 0.2e1;
t397 = t140 / 0.2e1;
t396 = t141 / 0.2e1;
t393 = t171 / 0.2e1;
t452 = t224 / 0.2e1;
t450 = mrSges(7,3) * t35;
t373 = qJD(2) / 0.2e1;
t36 = mrSges(6,1) * t168 - mrSges(6,3) * t47;
t37 = -t168 * mrSges(7,1) + t47 * mrSges(7,2);
t441 = t37 - t36;
t38 = -mrSges(6,2) * t168 - mrSges(6,3) * t48;
t39 = -mrSges(7,2) * t48 + mrSges(7,3) * t168;
t440 = t38 + t39;
t368 = mrSges(7,2) * t113;
t94 = mrSges(7,3) * t196 - t368;
t364 = mrSges(6,3) * t113;
t95 = -mrSges(6,2) * t196 - t364;
t370 = t94 + t95;
t363 = mrSges(6,3) * t456;
t96 = mrSges(6,1) * t196 - t363;
t367 = mrSges(7,2) * t456;
t97 = -mrSges(7,1) * t196 + t367;
t369 = t97 - t96;
t438 = Ifges(5,5) * t179;
t437 = t180 * Ifges(5,6);
t350 = qJDD(2) / 0.2e1;
t169 = t222 * t257 - t304;
t342 = t199 * t256;
t121 = -pkin(4) * t342 + t169;
t436 = t433 * pkin(5) + t432 * qJ(6) - qJD(6) * t220 - t121;
t435 = Ifges(4,5) * qJD(2);
t434 = Ifges(4,6) * qJD(2);
t365 = mrSges(4,3) * t201;
t431 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t180 - t179 * mrSges(5,2) - t365;
t430 = t443 * t168 + t444 * t48 + t445 * t47;
t215 = -pkin(7) * t307 + t247;
t429 = t215 * t265 + t216 * t262;
t126 = -t199 * mrSges(5,2) + mrSges(5,3) * t180;
t127 = mrSges(5,1) * t199 - mrSges(5,3) * t179;
t428 = t126 * t258 - t127 * t256;
t427 = -t256 * t33 + t258 * t34;
t426 = Ifges(5,3) * t199 + t444 * t113 + t443 * t196 + t445 * t456 + t437 + t438;
t425 = 0.2e1 * t350;
t293 = mrSges(5,1) * t256 + mrSges(5,2) * t258;
t424 = -m(3) * pkin(7) + m(5) * t259 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t293;
t276 = -t257 * t262 + t303;
t161 = -pkin(3) * t276 - qJ(4) * t219 - t246;
t177 = t257 * t232 - t234 * t349;
t105 = t258 * t161 - t177 * t256;
t335 = t219 * t258;
t79 = -pkin(4) * t276 - pkin(8) * t335 + t105;
t106 = t256 * t161 + t258 * t177;
t336 = t219 * t256;
t85 = -pkin(8) * t336 + t106;
t371 = t261 * t79 + t264 * t85;
t200 = t219 * qJD(2);
t202 = t276 * qJD(2);
t339 = t202 * t258;
t314 = pkin(2) * t326;
t117 = pkin(3) * t200 - qJ(4) * t202 - qJD(4) * t219 + t314;
t305 = qJD(2) * t259;
t197 = t262 * t305 + t324;
t198 = t265 * t305 - t325;
t139 = t197 * t349 + t257 * t198;
t73 = t258 * t117 - t139 * t256;
t51 = pkin(4) * t200 - pkin(8) * t339 + t73;
t340 = t202 * t256;
t74 = t256 * t117 + t258 * t139;
t62 = -pkin(8) * t340 + t74;
t9 = -qJD(5) * t371 - t261 * t62 + t264 * t51;
t423 = m(7) * pkin(5) + t449;
t295 = t252 * mrSges(4,1) - t250 * mrSges(4,2);
t422 = t295 + (mrSges(5,3) - t465) * t250;
t421 = m(7) * qJ(6) - t448;
t322 = qJD(5) * t264;
t323 = qJD(5) * t261;
t3 = t261 * t20 + t264 * t26 + t59 * t322 - t323 * t69;
t1 = qJ(6) * t168 + qJD(6) * t196 + t3;
t419 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t358 = Ifges(5,4) * t258;
t288 = -Ifges(5,2) * t256 + t358;
t359 = Ifges(5,4) * t256;
t290 = Ifges(5,1) * t258 - t359;
t382 = t258 / 0.2e1;
t418 = (Ifges(5,1) * t179 + Ifges(5,4) * t180 + Ifges(5,5) * t199) * t382 - t256 * (Ifges(5,4) * t179 + Ifges(5,2) * t180 + Ifges(5,6) * t199) / 0.2e1 + t144 * t293 + t229 * mrSges(4,2) + t179 * t290 / 0.2e1 + t180 * t288 / 0.2e1;
t391 = t196 / 0.2e1;
t399 = t456 / 0.2e1;
t401 = t113 / 0.2e1;
t402 = -t113 / 0.2e1;
t417 = Ifges(6,4) * t402 + Ifges(7,5) * t401 + t391 * t445 + t399 * t447 - t450;
t416 = -Ifges(6,2) * t402 + Ifges(7,3) * t401 + t391 * t444 + t399 * t446 - t455;
t414 = mrSges(6,1) * t63 + mrSges(7,1) * t7 + 0.2e1 * Ifges(7,3) * t407 + Ifges(6,4) * t467 + Ifges(6,6) * t466 + (t446 + Ifges(7,5)) * t409 + (t444 + Ifges(7,6)) * t394 + (-t408 + t407) * Ifges(6,2);
t21 = -t261 * t69 + t264 * t59;
t16 = -pkin(5) * t196 + qJD(6) - t21;
t17 = qJ(6) * t196 + t22;
t413 = t16 * mrSges(7,1) + t22 * mrSges(6,2) + t83 * mrSges(5,2) - t438 / 0.2e1 - t17 * mrSges(7,3) - t437 / 0.2e1 - t21 * mrSges(6,1) - t229 * mrSges(4,1) - t82 * mrSges(5,1);
t404 = Ifges(5,1) * t396 + Ifges(5,4) * t397 + Ifges(5,5) * t393;
t400 = -t456 / 0.2e1;
t392 = -t196 / 0.2e1;
t390 = -t199 / 0.2e1;
t389 = t199 / 0.2e1;
t386 = t201 / 0.2e1;
t381 = pkin(2) * t257;
t380 = pkin(2) * t262;
t379 = pkin(4) * t256;
t378 = pkin(7) * t265;
t375 = g(3) * t250;
t241 = qJ(4) + t381;
t372 = pkin(8) + t241;
t341 = t199 * t258;
t315 = pkin(2) * t328;
t137 = pkin(3) * t201 + qJ(4) * t199 + t315;
t170 = t222 * t349 + t205;
t91 = t258 * t137 - t170 * t256;
t68 = pkin(4) * t201 + pkin(8) * t341 + t91;
t92 = t256 * t137 + t258 * t170;
t81 = pkin(8) * t342 + t92;
t28 = t261 * t68 + t264 * t81;
t366 = mrSges(4,3) * t199;
t362 = Ifges(3,4) * t262;
t361 = Ifges(3,4) * t265;
t360 = Ifges(4,4) * t201;
t348 = qJ(4) * t250;
t333 = t250 * t266;
t332 = t252 * t266;
t331 = t263 * t249;
t330 = t263 * t251;
t329 = t266 * t249;
t327 = qJD(1) * t265;
t321 = m(5) - t453;
t318 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t328) * t378;
t310 = t349 * pkin(2);
t15 = t48 * mrSges(6,1) + t47 * mrSges(6,2);
t14 = t48 * mrSges(7,1) - t47 * mrSges(7,3);
t301 = t171 * mrSges(4,1) + t172 * mrSges(4,2);
t89 = -t140 * mrSges(5,1) + t141 * mrSges(5,2);
t138 = t197 * t257 - t349 * t198;
t176 = -t349 * t232 - t234 * t257;
t235 = t266 * t246;
t300 = -t263 * t259 + t235;
t245 = -t310 - pkin(3);
t108 = pkin(4) * t340 + t138;
t133 = pkin(4) * t336 + t176;
t296 = mrSges(3,1) * t262 + mrSges(3,2) * t265;
t289 = t265 * Ifges(3,2) + t362;
t287 = Ifges(3,5) * t265 - Ifges(3,6) * t262;
t286 = Ifges(5,5) * t258 - Ifges(5,6) * t256;
t284 = pkin(5) * t251 + qJ(6) * t249;
t282 = t256 * t82 - t258 * t83;
t27 = -t261 * t81 + t264 * t68;
t31 = -t261 * t85 + t264 * t79;
t210 = t372 * t256;
t211 = t372 * t258;
t279 = -t264 * t210 - t211 * t261;
t153 = -t210 * t261 + t211 * t264;
t278 = pkin(1) * t296;
t8 = t261 * t51 + t264 * t62 + t79 * t322 - t323 * t85;
t277 = t262 * (Ifges(3,1) * t265 - t362);
t228 = t245 - t374;
t248 = Ifges(3,4) * t327;
t231 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t327;
t209 = Ifges(3,1) * t328 + Ifges(3,5) * qJD(2) + t248;
t208 = Ifges(3,6) * qJD(2) + qJD(1) * t289;
t195 = Ifges(4,4) * t199;
t185 = -qJD(2) * mrSges(4,2) - t366;
t184 = t251 * t332 + t331;
t183 = t252 * t329 - t330;
t182 = t252 * t330 - t329;
t181 = t251 * t266 + t252 * t331;
t159 = mrSges(4,1) * t199 + mrSges(4,2) * t201;
t155 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t172;
t154 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t171;
t148 = t201 * Ifges(4,1) - t195 + t435;
t147 = -t199 * Ifges(4,2) + t360 + t434;
t146 = t218 * t219;
t145 = t220 * t219;
t142 = t218 * pkin(5) - t220 * qJ(6) + t228;
t119 = qJD(4) * t220 + qJD(5) * t153;
t118 = -qJD(4) * t218 + qJD(5) * t279;
t101 = mrSges(5,1) * t171 - mrSges(5,3) * t141;
t100 = -mrSges(5,2) * t171 + mrSges(5,3) * t140;
t88 = t202 * t220 + t322 * t335 - t323 * t336;
t87 = -t202 * t218 - t204 * t219;
t70 = t141 * Ifges(5,4) + t140 * Ifges(5,2) + t171 * Ifges(5,6);
t67 = mrSges(6,1) * t113 + mrSges(6,2) * t456;
t66 = mrSges(7,1) * t113 - mrSges(7,3) * t456;
t65 = pkin(5) * t456 + qJ(6) * t113;
t61 = pkin(5) * t145 + qJ(6) * t146 + t133;
t30 = pkin(5) * t276 - t31;
t29 = -qJ(6) * t276 + t371;
t24 = -pkin(5) * t201 - t27;
t23 = qJ(6) * t201 + t28;
t18 = pkin(5) * t88 - qJ(6) * t87 + qJD(6) * t146 + t108;
t6 = -pkin(5) * t200 - t9;
t5 = qJ(6) * t200 - qJD(6) * t276 + t8;
t10 = [t289 * t452 + (t421 * t181 + t423 * t182 + (-m(5) * (-t246 - t348) + m(4) * t246 + t422 + (-t246 - t459) * t453 - t457) * t263 + (t453 * (-t259 + t379) + m(4) * t259 + t424) * t266) * g(1) + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t414) * t145 + (mrSges(6,1) * t107 - mrSges(7,2) * t17 - mrSges(6,3) * t22 + t416) * t88 + (Ifges(4,5) * t373 + t148 / 0.2e1 + Ifges(4,1) * t386 + t286 * t389 - mrSges(4,3) * t162 + Ifges(4,4) * t390 + t418) * t202 + (-t413 - mrSges(4,3) * t163 - t147 / 0.2e1 + t443 * t391 + t445 * t399 - Ifges(4,6) * t373 - Ifges(4,4) * t386 + Ifges(5,3) * t389 - Ifges(4,2) * t390 + Ifges(7,6) * t401 + Ifges(6,6) * t402 + t426 / 0.2e1) * t200 + (mrSges(6,2) * t107 + t16 * mrSges(7,2) - t21 * mrSges(6,3) + t417 + t439 / 0.2e1) * t87 + (t287 * t373 - t318) * qJD(2) + (t224 * t378 + t429) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t429) + (-m(4) * t162 + m(5) * t144 - t431) * t138 + (-m(4) * t98 + m(5) * t93 - t155 + t89) * t176 + m(7) * (t1 * t29 + t16 * t6 + t17 * t5 + t18 * t35 + t2 * t30 + t61 * t7) - (Ifges(5,5) * t141 + Ifges(5,6) * t140 + Ifges(5,3) * t171 + t430) * t276 / 0.2e1 - (t192 * mrSges(4,1) + t33 * mrSges(5,1) - t34 * mrSges(5,2) - t99 * mrSges(4,3) - Ifges(4,4) * t172 + Ifges(5,5) * t396 + Ifges(4,2) * t171 - t425 * Ifges(4,6) + Ifges(5,6) * t397 + Ifges(6,6) * t408 + Ifges(7,6) * t407 + Ifges(5,3) * t393 + t443 * t394 + t445 * t409 + t419) * t276 + (t265 * (-Ifges(3,2) * t262 + t361) + t277) * t320 / 0.2e1 + m(6) * (t107 * t108 + t133 * t63 + t21 * t9 + t22 * t8 + t3 * t371 + t31 * t4) + t371 * t38 + t462 * t146 + (-m(4) * t300 - m(5) * t235 + t465 * t333 + t453 * (t244 * t332 - t260 * t333 + t263 * t379 + t300) - t423 * t184 - t421 * t183 + (-t461 * t250 - t295 + t457) * t266 + t424 * t263) * g(2) + t225 * t361 / 0.2e1 - qJDD(2) * mrSges(3,2) * t378 + t159 * t314 + (t192 * mrSges(4,2) - t98 * mrSges(4,3) + Ifges(4,1) * t172 - Ifges(4,4) * t171 + Ifges(4,5) * t425 + t286 * t393 + t288 * t397 + t290 * t396 + t93 * t293) * t219 + Ifges(2,3) * qJDD(1) + t265 * (Ifges(3,4) * t225 + Ifges(3,2) * t224 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t233 * t347 - t70 * t336 / 0.2e1 - t208 * t326 / 0.2e1 - t278 * t320 - pkin(1) * (-mrSges(3,1) * t224 + mrSges(3,2) * t225) - t231 * t313 + t139 * t185 + t177 * t154 - t246 * t301 + t133 * t15 + (-t33 * t335 - t336 * t34 - t339 * t82 - t340 * t83) * mrSges(5,3) + t74 * t126 + t73 * t127 + t105 * t101 + t106 * t100 + t108 * t67 + t5 * t94 + t8 * t95 + t9 * t96 + t6 * t97 + t18 * t66 + t61 * t14 + t31 * t36 + t30 * t37 + t29 * t39 + m(5) * (t105 * t33 + t106 * t34 + t73 * t82 + t74 * t83) + m(4) * (t139 * t163 + t177 * t99 - t192 * t246 + t229 * t314) + t265 * t209 * t373 + Ifges(3,6) * t265 * t350 + (Ifges(3,1) * t225 + Ifges(3,4) * t452 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t225) + t425 * Ifges(3,5)) * t262 + t335 * t404; (t208 / 0.2e1 + pkin(7) * t231) * t328 + (t318 + (t278 - t277 / 0.2e1) * qJD(1)) * qJD(1) + (-t195 + t148) * t389 + (t21 * t432 - t218 * t3 - t22 * t433) * mrSges(6,3) + t414 * t218 + t416 * t204 + (-t1 * t218 - t16 * t432 - t17 * t433) * mrSges(7,2) - (-Ifges(6,2) * t401 + Ifges(7,3) * t402 + t444 * t392 + t446 * t400 + t455) * t131 + (t435 / 0.2e1 - t286 * t390 + t418) * t199 + t436 * t66 + t369 * t119 + t370 * t118 + t440 * t153 + (-t341 * t82 - t342 * t83 + t427) * mrSges(5,3) + (-t282 * qJD(4) - t144 * t169 + t241 * t427 + t245 * t93 - t82 * t91 - t83 * t92) * m(5) + t428 * qJD(4) + t431 * t169 + (t162 * t169 - t163 * t170 - t229 * t315 + (t257 * t99 + t349 * t98) * pkin(2)) * m(4) + (mrSges(6,1) * t433 - mrSges(6,2) * t432) * t107 + (t1 * t153 + t142 * t7 - t279 * t2 + t436 * t35 + (t118 - t23) * t17 + (t119 - t24) * t16) * m(7) - t441 * t279 + (-t107 * t121 + t279 * t4 + t153 * t3 + t228 * t63 + (t118 - t28) * t22 + (-t119 - t27) * t21) * m(6) + t458 * (t296 + (m(4) + m(5)) * t380 + t453 * (-t252 * t260 - t380) + (mrSges(4,2) - t461 + t465) * t252 + (mrSges(4,1) + t274 + m(6) * t244 - m(7) * (-t244 - t284) + t463) * t250) - (-Ifges(3,2) * t328 + t209 + t248) * t327 / 0.2e1 + (-m(4) * t253 - m(5) * (t253 + t348) + t453 * (t253 + t459) + (-m(7) * t284 - t463) * t252 - t422 + t464) * g(3) + (t413 - Ifges(4,2) * t389 + Ifges(5,3) * t390 + Ifges(6,6) * t401 + Ifges(7,6) * t402 + t443 * t392 + t445 * t400 + t434 / 0.2e1) * t201 + t163 * t365 + t155 * t310 - t462 * t220 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) - (-Ifges(4,1) * t199 - t360 + t426) * t201 / 0.2e1 - t162 * t366 + t245 * t89 - t287 * t320 / 0.2e1 + Ifges(3,6) * t224 + Ifges(3,5) * t225 + t228 * t15 - t215 * mrSges(3,2) - t216 * mrSges(3,1) - t159 * t315 - t170 * t185 + Ifges(4,5) * t172 - Ifges(4,6) * t171 + t142 * t14 - t92 * t126 - t91 * t127 - t121 * t67 - t23 * t94 - t28 * t95 - t27 * t96 - t24 * t97 + t98 * mrSges(4,1) - t99 * mrSges(4,2) + t93 * t294 + t154 * t381 + t70 * t382 + t147 * t386 + (Ifges(6,4) * t401 + Ifges(7,5) * t402 + t392 * t445 + t400 * t447 + t450 + t460) * t132 + (-t417 + t460) * t203 + (t258 * t100 - t256 * t101) * t241 + (Ifges(5,5) * t256 + Ifges(5,6) * t258) * t393 + (Ifges(5,1) * t256 + t358) * t396 + (Ifges(5,2) * t258 + t359) * t397 + t256 * t404; (-t66 - t67 + t431) * t201 + (t185 + t428) * t199 + t301 + t440 * t220 + t441 * t218 + t258 * t101 + t256 * t100 - t432 * t370 + t433 * t369 + (-g(1) * t263 + g(2) * t266) * (m(4) + t321) + (t1 * t220 + t16 * t433 - t17 * t432 + t2 * t218 - t201 * t35) * m(7) + (-t107 * t201 - t21 * t433 - t218 * t4 - t22 * t432 + t220 * t3) * m(6) + (-t144 * t201 - t199 * t282 + t256 * t34 + t258 * t33) * m(5) + (t162 * t201 + t163 * t199 + t192) * m(4); -t369 * t456 + t370 * t113 - t180 * t126 + t179 * t127 + t14 + t15 + t89 + (t113 * t17 - t16 * t456 + t7) * m(7) + (t113 * t22 + t21 * t456 + t63) * m(6) + (t179 * t82 - t180 * t83 + t93) * m(5) + (t252 * g(3) - t250 * t458) * t321; t419 + (t249 * t449 + t251 * t448) * t375 + (-m(7) * t17 - t364 - t370) * t21 + (-Ifges(6,2) * t456 - t112 + t439) * t401 + t430 + t17 * t367 + t16 * t368 + (-m(7) * t16 + t363 - t369) * t22 - t35 * (mrSges(7,1) * t456 + mrSges(7,3) * t113) - t107 * (mrSges(6,1) * t456 - mrSges(6,2) * t113) + qJD(6) * t94 - t65 * t66 - pkin(5) * t37 + qJ(6) * t39 + (t183 * t449 + t184 * t448) * g(1) + (t181 * t449 + t182 * t448) * g(2) + (-t113 * t445 + t444 * t456) * t392 + t56 * t399 + (Ifges(7,3) * t456 - t356) * t402 + (-t35 * t65 - (-pkin(5) * t249 + qJ(6) * t251) * t375 - g(1) * (-pkin(5) * t183 + qJ(6) * t184) - g(2) * (-pkin(5) * t181 + qJ(6) * t182) - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t17) * m(7) + (-t113 * t447 + t111 - t357 + t53) * t400; t456 * t66 - t196 * t94 + (-g(1) * t183 - g(2) * t181 - t17 * t196 - t249 * t375 + t35 * t456 + t2) * m(7) + t37;];
tau  = t10;
