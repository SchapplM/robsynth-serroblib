% Calculate vector of inverse dynamics joint torques for
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:48
% EndTime: 2019-12-31 22:04:29
% DurationCPUTime: 21.74s
% Computational Cost: add. (6985->644), mult. (15437->837), div. (0->0), fcn. (10226->10), ass. (0->312)
t478 = mrSges(5,1) + mrSges(6,1);
t477 = mrSges(5,2) - mrSges(6,3);
t454 = -Ifges(5,4) + Ifges(6,5);
t484 = t454 + Ifges(6,5);
t242 = sin(qJ(2));
t246 = cos(qJ(2));
t289 = t246 * mrSges(3,1) - t242 * mrSges(3,2);
t476 = -mrSges(5,3) - mrSges(6,2);
t483 = t476 * t242 - t289;
t241 = sin(qJ(3));
t245 = cos(qJ(3));
t334 = qJD(1) * t242;
t309 = t245 * t334;
t181 = qJD(2) * t241 + t309;
t324 = qJD(1) * qJD(2);
t189 = qJDD(1) * t242 + t246 * t324;
t100 = -qJD(3) * t181 + qJDD(2) * t245 - t189 * t241;
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t306 = t241 * t334;
t331 = qJD(2) * t245;
t265 = t306 - t331;
t108 = t240 * t181 + t244 * t265;
t257 = t265 * qJD(3);
t99 = qJDD(2) * t241 + t189 * t245 - t257;
t34 = -qJD(4) * t108 + t240 * t100 + t244 * t99;
t414 = t34 / 0.2e1;
t253 = t244 * t181 - t240 * t265;
t35 = qJD(4) * t253 - t244 * t100 + t240 * t99;
t412 = t35 / 0.2e1;
t188 = qJDD(1) * t246 - t242 * t324;
t178 = qJDD(3) - t188;
t170 = qJDD(4) + t178;
t394 = t170 / 0.2e1;
t176 = t188 * pkin(6);
t482 = t176 * t246;
t481 = t242 * mrSges(4,3) + mrSges(2,1) - t483;
t230 = pkin(6) * t334;
t198 = -qJD(2) * pkin(2) + t230;
t127 = pkin(3) * t265 + t198;
t333 = qJD(1) * t246;
t216 = qJD(3) - t333;
t206 = qJD(4) + t216;
t291 = pkin(2) * t246 + pkin(7) * t242;
t192 = -pkin(1) - t291;
t171 = t192 * qJD(1);
t231 = pkin(6) * t333;
t199 = qJD(2) * pkin(7) + t231;
t114 = t241 * t171 + t245 * t199;
t84 = -pkin(8) * t265 + t114;
t365 = t240 * t84;
t113 = t245 * t171 - t199 * t241;
t83 = -pkin(8) * t181 + t113;
t75 = pkin(3) * t216 + t83;
t36 = t244 * t75 - t365;
t470 = qJD(5) - t36;
t27 = -pkin(4) * t206 + t470;
t40 = t108 * pkin(4) - qJ(5) * t253 + t127;
t103 = Ifges(5,4) * t108;
t367 = Ifges(6,5) * t108;
t453 = Ifges(6,4) + Ifges(5,5);
t455 = Ifges(5,1) + Ifges(6,1);
t448 = t206 * t453 + t253 * t455 - t103 + t367;
t480 = t127 * mrSges(5,2) - t40 * mrSges(6,3) + mrSges(6,2) * t27 - mrSges(5,3) * t36 + t448 / 0.2e1;
t479 = m(6) + m(5);
t452 = -Ifges(5,6) + Ifges(6,6);
t451 = Ifges(5,3) + Ifges(6,2);
t389 = t206 / 0.2e1;
t398 = t253 / 0.2e1;
t401 = t108 / 0.2e1;
t402 = -t108 / 0.2e1;
t475 = Ifges(5,4) * t402 + Ifges(6,5) * t401 + t389 * t453 + t398 * t455 + t480;
t360 = t244 * t84;
t37 = t240 * t75 + t360;
t28 = qJ(5) * t206 + t37;
t102 = Ifges(6,5) * t253;
t49 = t206 * Ifges(6,6) + t108 * Ifges(6,3) + t102;
t368 = Ifges(5,4) * t253;
t52 = -t108 * Ifges(5,2) + t206 * Ifges(5,6) + t368;
t474 = -mrSges(6,2) * t28 - mrSges(5,3) * t37 + mrSges(5,1) * t127 + mrSges(6,1) * t40 + t49 / 0.2e1 - t52 / 0.2e1;
t357 = qJDD(1) * pkin(1);
t116 = -pkin(2) * t188 - pkin(7) * t189 - t357;
t152 = qJDD(2) * pkin(7) + t176;
t44 = -qJD(3) * t114 + t245 * t116 - t152 * t241;
t20 = pkin(3) * t178 - pkin(8) * t99 + t44;
t327 = qJD(3) * t245;
t329 = qJD(3) * t241;
t43 = t241 * t116 + t245 * t152 + t171 * t327 - t199 * t329;
t26 = pkin(8) * t100 + t43;
t6 = -qJD(4) * t37 + t20 * t244 - t240 * t26;
t3 = -pkin(4) * t170 + qJDD(5) - t6;
t413 = -t35 / 0.2e1;
t177 = t189 * pkin(6);
t153 = -qJDD(2) * pkin(2) + t177;
t76 = -pkin(3) * t100 + t153;
t7 = pkin(4) * t35 - qJ(5) * t34 - qJD(5) * t253 + t76;
t473 = mrSges(5,2) * t76 + mrSges(6,2) * t3 - mrSges(5,3) * t6 - mrSges(6,3) * t7 + Ifges(5,4) * t413 + 0.2e1 * t394 * t453 + t412 * t484 + 0.2e1 * t414 * t455;
t239 = qJ(3) + qJ(4);
t233 = sin(t239);
t234 = cos(t239);
t287 = -mrSges(4,1) * t245 + mrSges(4,2) * t241;
t472 = m(4) * pkin(2) - t477 * t233 + t234 * t478 - t287;
t330 = qJD(2) * t246;
t305 = t241 * t330;
t260 = t242 * t327 + t305;
t243 = sin(qJ(1));
t247 = cos(qJ(1));
t342 = t246 * t247;
t162 = -t241 * t342 + t243 * t245;
t310 = t241 * t333;
t442 = -t231 + (-t310 + t329) * pkin(3);
t467 = -Ifges(5,2) * t402 + Ifges(6,3) * t401 + t389 * t452 + t398 * t454 + t474;
t390 = -t206 / 0.2e1;
t399 = -t253 / 0.2e1;
t466 = -Ifges(5,2) * t401 + Ifges(6,3) * t402 + t390 * t452 + t399 * t454 - t474;
t61 = pkin(4) * t253 + qJ(5) * t108;
t463 = Ifges(5,4) * t401 + Ifges(6,5) * t402 + t390 * t453 + t399 * t455 - t480;
t404 = t99 / 0.2e1;
t462 = -m(3) - m(4);
t403 = t100 / 0.2e1;
t393 = t178 / 0.2e1;
t460 = t188 / 0.2e1;
t459 = t189 / 0.2e1;
t458 = t265 / 0.2e1;
t232 = pkin(6) * t330;
t456 = -mrSges(3,3) + mrSges(2,2);
t182 = t240 * t241 - t244 * t245;
t435 = qJD(3) + qJD(4);
t117 = t435 * t182;
t183 = t240 * t245 + t241 * t244;
t118 = t435 * t183;
t139 = t183 * t333;
t264 = t182 * t246;
t140 = qJD(1) * t264;
t449 = -qJD(5) * t183 + t442 + (t117 - t140) * qJ(5) + (t118 - t139) * pkin(4);
t85 = -mrSges(6,2) * t108 + mrSges(6,3) * t206;
t375 = mrSges(5,3) * t108;
t86 = -mrSges(5,2) * t206 - t375;
t447 = t85 + t86;
t374 = mrSges(5,3) * t253;
t87 = mrSges(5,1) * t206 - t374;
t88 = -mrSges(6,1) * t206 + mrSges(6,2) * t253;
t446 = t87 - t88;
t150 = t183 * t242;
t180 = t245 * t192;
t348 = t242 * t245;
t112 = -pkin(8) * t348 + t180 + (-pkin(6) * t241 - pkin(3)) * t246;
t343 = t245 * t246;
t220 = pkin(6) * t343;
t130 = t241 * t192 + t220;
t352 = t241 * t242;
t120 = -pkin(8) * t352 + t130;
t445 = t240 * t112 + t244 * t120;
t444 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t265 + t181 * mrSges(4,2) + mrSges(3,3) * t334;
t229 = Ifges(3,4) * t333;
t175 = Ifges(4,4) * t265;
t93 = t181 * Ifges(4,1) + t216 * Ifges(4,5) - t175;
t443 = Ifges(3,1) * t334 + Ifges(3,5) * qJD(2) + t245 * t93 + t229;
t328 = qJD(3) * t242;
t261 = -t241 * t328 + t245 * t330;
t441 = t170 * t451 + t34 * t453 + t35 * t452;
t341 = t247 * t233;
t147 = -t243 * t234 + t246 * t341;
t148 = t233 * t243 + t234 * t342;
t440 = t147 * t478 + t148 * t477;
t344 = t243 * t246;
t145 = t233 * t344 + t234 * t247;
t146 = t234 * t344 - t341;
t439 = t145 * t478 + t146 * t477;
t438 = t177 * t242 + t482;
t437 = -t241 * t44 + t245 * t43;
t434 = Ifges(4,5) * t181 - Ifges(4,6) * t265 + Ifges(4,3) * t216 + t108 * t452 + t206 * t451 + t253 * t453;
t269 = pkin(3) * t242 - pkin(8) * t343;
t290 = pkin(2) * t242 - pkin(7) * t246;
t187 = t290 * qJD(2);
t332 = qJD(2) * t242;
t315 = pkin(6) * t332;
t336 = t245 * t187 + t241 * t315;
t64 = t269 * qJD(2) + (-t220 + (pkin(8) * t242 - t192) * t241) * qJD(3) + t336;
t81 = t241 * t187 + t192 * t327 + (-t242 * t331 - t246 * t329) * pkin(6);
t68 = -pkin(8) * t260 + t81;
t17 = -qJD(4) * t445 - t240 * t68 + t244 * t64;
t432 = m(6) * pkin(4) + t478;
t430 = m(6) * qJ(5) - t477;
t429 = -t44 * mrSges(4,1) + t43 * mrSges(4,2);
t325 = qJD(4) * t244;
t326 = qJD(4) * t240;
t5 = t240 * t20 + t244 * t26 + t75 * t325 - t326 * t84;
t2 = qJ(5) * t170 + qJD(5) * t206 + t5;
t426 = -t6 * mrSges(5,1) + t3 * mrSges(6,1) + t5 * mrSges(5,2) - t2 * mrSges(6,3);
t373 = Ifges(3,4) * t242;
t281 = t246 * Ifges(3,2) + t373;
t425 = t27 * mrSges(6,1) + t37 * mrSges(5,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t281 / 0.2e1 - t28 * mrSges(6,3) - t36 * mrSges(5,1);
t419 = mrSges(5,1) * t76 + mrSges(6,1) * t7 - mrSges(6,2) * t2 - mrSges(5,3) * t5 + 0.2e1 * Ifges(6,3) * t412 - t34 * Ifges(5,4) / 0.2e1 - t170 * Ifges(5,6) / 0.2e1 + t484 * t414 + (t452 + Ifges(6,6)) * t394 + (-t413 + t412) * Ifges(5,2);
t415 = m(5) * pkin(3);
t411 = Ifges(4,1) * t404 + Ifges(4,4) * t403 + Ifges(4,5) * t393;
t371 = Ifges(4,4) * t181;
t92 = -Ifges(4,2) * t265 + Ifges(4,6) * t216 + t371;
t405 = -t92 / 0.2e1;
t248 = -pkin(8) - pkin(7);
t392 = t181 / 0.2e1;
t383 = pkin(3) * t181;
t382 = pkin(3) * t240;
t381 = pkin(3) * t244;
t235 = t242 * pkin(6);
t378 = -qJD(1) / 0.2e1;
t377 = qJD(3) / 0.2e1;
t185 = t290 * qJD(1);
t164 = t241 * t185;
t350 = t241 * t246;
t111 = t164 + (-pkin(6) * t348 - pkin(8) * t350) * qJD(1);
t123 = pkin(6) * t306 + t245 * t185;
t94 = qJD(1) * t269 + t123;
t56 = t244 * t111 + t240 * t94;
t376 = mrSges(5,2) * t234;
t372 = Ifges(3,4) * t246;
t370 = Ifges(4,4) * t241;
t369 = Ifges(4,4) * t245;
t366 = t114 * mrSges(4,3);
t354 = t233 * t242;
t353 = t234 * t242;
t351 = t241 * t243;
t349 = t241 * t247;
t346 = t242 * t248;
t226 = pkin(3) * t245 + pkin(2);
t203 = t246 * t226;
t217 = pkin(3) * t352;
t190 = t235 + t217;
t335 = t247 * pkin(1) + t243 * pkin(6);
t319 = Ifges(4,5) * t99 + Ifges(4,6) * t100 + Ifges(4,3) * t178;
t317 = pkin(3) * t326;
t316 = pkin(3) * t325;
t314 = m(4) * pkin(7) + mrSges(4,3);
t128 = pkin(3) * t260 + t232;
t311 = qJD(3) * t248;
t22 = -t170 * mrSges(6,1) + t34 * mrSges(6,2);
t296 = t324 / 0.2e1;
t295 = -t145 * pkin(4) + qJ(5) * t146;
t294 = -t147 * pkin(4) + qJ(5) * t148;
t292 = t245 * t311;
t288 = mrSges(3,1) * t242 + mrSges(3,2) * t246;
t286 = mrSges(4,1) * t241 + mrSges(4,2) * t245;
t283 = Ifges(4,1) * t245 - t370;
t282 = Ifges(4,1) * t241 + t369;
t280 = -Ifges(4,2) * t241 + t369;
t279 = Ifges(4,2) * t245 + t370;
t278 = Ifges(3,5) * t246 - Ifges(3,6) * t242;
t277 = Ifges(4,5) * t245 - Ifges(4,6) * t241;
t276 = Ifges(4,5) * t241 + Ifges(4,6) * t245;
t275 = pkin(4) * t234 + qJ(5) * t233;
t55 = -t111 * t240 + t244 * t94;
t65 = t112 * t244 - t120 * t240;
t200 = t248 * t241;
t201 = t248 * t245;
t271 = t244 * t200 + t201 * t240;
t122 = t200 * t240 - t201 * t244;
t270 = t162 * pkin(3);
t268 = pkin(1) * t288;
t160 = t241 * t344 + t245 * t247;
t16 = t112 * t325 - t120 * t326 + t240 * t64 + t244 * t68;
t266 = t242 * (Ifges(3,1) * t246 - t373);
t262 = t160 * pkin(3);
t258 = t265 * mrSges(4,3);
t256 = Ifges(4,5) * t242 + t246 * t283;
t255 = Ifges(4,6) * t242 + t246 * t280;
t254 = Ifges(4,3) * t242 + t246 * t277;
t252 = -t426 + t441;
t237 = t247 * pkin(6);
t225 = -pkin(4) - t381;
t222 = qJ(5) + t382;
t214 = qJD(5) + t316;
t204 = mrSges(6,3) * t353;
t202 = qJ(5) * t353;
t196 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t333;
t186 = t241 * t311;
t172 = t286 * t242;
t163 = t245 * t342 + t351;
t161 = -t243 * t343 + t349;
t151 = t182 * t242;
t129 = -pkin(6) * t350 + t180;
t126 = mrSges(4,1) * t216 - mrSges(4,3) * t181;
t125 = -t216 * mrSges(4,2) - t258;
t124 = -pkin(6) * t309 + t164;
t104 = pkin(4) * t182 - qJ(5) * t183 - t226;
t82 = -qJD(3) * t130 + t336;
t79 = pkin(4) * t150 + qJ(5) * t151 + t190;
t78 = -mrSges(4,2) * t178 + mrSges(4,3) * t100;
t77 = mrSges(4,1) * t178 - mrSges(4,3) * t99;
t74 = qJD(4) * t122 + t186 * t240 - t244 * t292;
t73 = qJD(4) * t271 + t244 * t186 + t240 * t292;
t70 = -t326 * t352 + (t348 * t435 + t305) * t244 + t261 * t240;
t69 = -qJD(2) * t264 - t150 * t435;
t63 = mrSges(5,1) * t108 + mrSges(5,2) * t253;
t62 = mrSges(6,1) * t108 - mrSges(6,3) * t253;
t60 = pkin(4) * t246 - t65;
t59 = -qJ(5) * t246 + t445;
t57 = -mrSges(4,1) * t100 + mrSges(4,2) * t99;
t48 = -pkin(4) * t334 - t55;
t47 = qJ(5) * t334 + t56;
t46 = t383 + t61;
t41 = t99 * Ifges(4,4) + t100 * Ifges(4,2) + t178 * Ifges(4,6);
t39 = t244 * t83 - t365;
t38 = t240 * t83 + t360;
t24 = -mrSges(6,2) * t35 + mrSges(6,3) * t170;
t23 = -mrSges(5,2) * t170 - mrSges(5,3) * t35;
t21 = mrSges(5,1) * t170 - mrSges(5,3) * t34;
t18 = pkin(4) * t70 - qJ(5) * t69 + qJD(5) * t151 + t128;
t15 = -pkin(4) * t332 - t17;
t14 = qJ(5) * t332 - qJD(5) * t246 + t16;
t13 = mrSges(5,1) * t35 + mrSges(5,2) * t34;
t12 = mrSges(6,1) * t35 - mrSges(6,3) * t34;
t1 = [(-mrSges(3,1) * t235 + Ifges(3,5) * t242 + (-mrSges(3,2) * pkin(6) + Ifges(3,6)) * t246) * qJDD(2) + t475 * t69 + (t189 * t235 + t438 + t482) * mrSges(3,3) + m(5) * (t127 * t128 + t16 * t37 + t17 * t36 + t190 * t76 + t445 * t5 + t6 * t65) + t445 * t23 + t289 * t357 + (-t161 * mrSges(4,1) - t160 * mrSges(4,2) - t479 * (pkin(3) * t349 + t243 * t346 + t237) + t456 * t247 + t462 * t237 + t432 * t146 + t430 * t145 + (m(3) * pkin(1) - m(4) * t192 - t479 * (-pkin(1) - t203) + t481) * t243) * g(1) + (-t163 * mrSges(4,1) - t162 * mrSges(4,2) + t462 * t335 - t479 * (pkin(3) * t351 + t226 * t342 + t335) + t456 * t243 - t432 * t148 - t430 * t147 + (-m(4) * t291 + t346 * t479 - t481) * t247) * g(2) + t216 * (qJD(2) * t254 - t276 * t328) / 0.2e1 - t265 * (qJD(2) * t255 - t279 * t328) / 0.2e1 - t268 * t324 + t198 * (mrSges(4,1) * t260 + mrSges(4,2) * t261) + t419 * t150 + t372 * t459 + t281 * t460 - t41 * t352 / 0.2e1 + t266 * t296 + qJD(2) ^ 2 * t278 / 0.2e1 - t473 * t151 + m(6) * (t14 * t28 + t15 * t27 + t18 * t40 + t2 * t59 + t3 * t60 + t7 * t79) - t196 * t315 + t467 * t70 + (-t113 * t261 - t114 * t260 - t348 * t44 - t352 * t43) * mrSges(4,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t438) - (t319 + t441) * t246 / 0.2e1 + t443 * t330 / 0.2e1 + t444 * t232 + m(4) * (t113 * t82 + t114 * t81 + t129 * t44 + t130 * t43 + t198 * t232) + ((-Ifges(3,2) * t242 + t372) * t296 + Ifges(3,4) * t459 + Ifges(3,2) * t460 - Ifges(4,3) * t393 - Ifges(4,6) * t403 - Ifges(4,5) * t404 - Ifges(6,6) * t412 - Ifges(5,6) * t413 - t453 * t414 - t451 * t394 + t426 + t429) * t246 + (m(4) * t153 * pkin(6) + Ifges(3,1) * t189 + Ifges(3,4) * t460 + t277 * t393 + t280 * t403 + t283 * t404) * t242 + t57 * t235 + (qJD(2) * t256 - t282 * t328) * t392 + t305 * t405 + t348 * t411 + (t434 / 0.2e1 + t113 * mrSges(4,1) - t114 * mrSges(4,2) + Ifges(5,6) * t402 + Ifges(6,6) * t401 + t389 * t451 + t398 * t453 - t425) * t332 - (t241 * t93 + t245 * t92) * t328 / 0.2e1 + t59 * t24 + t60 * t22 + t18 * t62 + t65 * t21 + t79 * t12 + t14 * t85 + t16 * t86 + t17 * t87 + t15 * t88 + t81 * t125 + t82 * t126 + t128 * t63 + t129 * t77 + t130 * t78 + t153 * t172 - pkin(1) * (-mrSges(3,1) * t188 + mrSges(3,2) * t189) + t190 * t13 + Ifges(2,3) * qJDD(1); (t256 * t378 + t283 * t377) * t181 - t475 * t117 + (-pkin(2) * t153 - t113 * t123 - t114 * t124 - t198 * t231) * m(4) - (t22 - t21) * t271 + (t271 * t6 + t122 * t5 - t226 * t76 + (-t56 + t73) * t37 + (-t55 - t74) * t36 + t442 * t127) * m(5) + (t104 * t7 - t271 * t3 + t122 * t2 + t449 * t40 + (-t47 + t73) * t28 + (-t48 + t74) * t27) * m(6) + (-t366 + t405) * t329 + (t255 * t458 - t113 * (mrSges(4,1) * t242 - mrSges(4,3) * t343) - t114 * (-mrSges(4,2) * t242 - mrSges(4,3) * t350) + (-t266 / 0.2e1 + t268) * qJD(1)) * qJD(1) + t93 * t327 / 0.2e1 - t278 * t324 / 0.2e1 + (t23 + t24) * t122 + t419 * t182 + t196 * t230 + (-t242 * t314 - t479 * (t203 - t346) + (-m(6) * t275 - t472) * t246 + t483) * g(3) + t153 * t287 - t280 * t257 / 0.2e1 + (t198 * t286 + t254 * t378 + t277 * t377) * t216 + t473 * t183 + t92 * t310 / 0.2e1 - t463 * t140 + t466 * t139 + t467 * t118 - t434 * t334 / 0.2e1 + (-t126 * t327 - t125 * t329 + m(4) * ((-t113 * t245 - t114 * t241) * qJD(3) + t437) - t241 * t77 + t245 * t78) * pkin(7) + (-t113 * t327 + t437) * mrSges(4,3) + t442 * t63 - (-Ifges(3,2) * t334 + t229 + t443) * t333 / 0.2e1 - t444 * t231 - t446 * t74 + t447 * t73 + t449 * t62 + (g(1) * t247 + g(2) * t243) * (t288 + (t248 * t479 - t314 + t476) * t246 + (m(5) * t226 - m(6) * (-t226 - t275) + t472) * t242) + (Ifges(5,6) * t401 + Ifges(6,6) * t402 + t390 * t451 + t399 * t453 + t425) * t334 + t276 * t393 + t279 * t403 + t282 * t404 + t241 * t411 + Ifges(3,3) * qJDD(2) - pkin(2) * t57 - t47 * t85 - t56 * t86 - t55 * t87 - t48 * t88 + t104 * t12 - t124 * t125 - t123 * t126 - t176 * mrSges(3,2) - t177 * mrSges(3,1) + Ifges(3,6) * t188 + Ifges(3,5) * t189 - t226 * t13 + t245 * t41 / 0.2e1; (-m(5) * t127 - t63) * t383 + (-Ifges(4,2) * t181 - t175 + t93) * t458 + t86 * t316 + t319 + t252 - t216 * (-Ifges(4,5) * t265 - Ifges(4,6) * t181) / 0.2e1 - t198 * (t181 * mrSges(4,1) - mrSges(4,2) * t265) - t463 * t108 + (t2 * t222 + t225 * t3 - t40 * t46 + (-t39 + t214) * t28 + (t317 - t38) * t27) * m(6) - t181 * (-Ifges(4,1) * t265 - t371) / 0.2e1 + t466 * t253 + (-m(6) * (-t262 + t295) + m(5) * t262 + mrSges(4,1) * t160 - mrSges(4,2) * t161 + t439) * g(2) + (-m(6) * (t270 + t294) - m(5) * t270 - mrSges(4,1) * t162 + mrSges(4,2) * t163 + t440) * g(1) - t429 + (m(5) * t36 + t446) * t38 - t446 * t317 + (-m(5) * t37 - t447) * t39 + t181 * t366 + t21 * t381 + t23 * t382 + t92 * t392 + (t240 * t5 + t244 * t6 + (-t240 * t36 + t244 * t37) * qJD(4)) * t415 + (-(-mrSges(5,1) * t233 - t241 * t415 - t376) * t242 - m(6) * (-pkin(4) * t354 + t202 - t217) + mrSges(6,1) * t354 - t204 + t172) * g(3) + (-t258 - t125) * t113 - t46 * t62 + t114 * t126 + t214 * t85 + t222 * t24 + t225 * t22; (t108 * t27 + t253 * t28) * mrSges(6,2) + (-t204 + (t233 * t432 + t376) * t242) * g(3) + (t374 + t446) * t37 + (-t375 - t447) * t36 + t252 + t440 * g(1) + t439 * g(2) + t52 * t398 + (Ifges(6,3) * t253 - t367) * t402 - pkin(4) * t22 + qJ(5) * t24 - t61 * t62 + qJD(5) * t85 - t40 * (mrSges(6,1) * t253 + mrSges(6,3) * t108) - t127 * (mrSges(5,1) * t253 - mrSges(5,2) * t108) + (-t108 * t453 + t253 * t452) * t390 + (-Ifges(5,2) * t253 - t103 + t448) * t401 + (-t108 * t455 + t102 - t368 + t49) * t399 + (-pkin(4) * t3 - t294 * g(1) - t295 * g(2) - t202 * g(3) + qJ(5) * t2 - t27 * t37 + t28 * t470 - t40 * t61) * m(6); t253 * t62 - t206 * t85 + (-g(1) * t147 - g(2) * t145 - g(3) * t354 - t28 * t206 + t253 * t40 + t3) * m(6) + t22;];
tau = t1;
