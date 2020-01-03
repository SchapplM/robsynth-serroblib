% Calculate vector of inverse dynamics joint torques for
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:31
% EndTime: 2019-12-31 21:08:00
% DurationCPUTime: 18.60s
% Computational Cost: add. (3276->606), mult. (7305->749), div. (0->0), fcn. (4299->6), ass. (0->276)
t421 = -Ifges(6,5) - Ifges(4,5);
t441 = Ifges(5,4) + t421;
t422 = -Ifges(5,5) - Ifges(6,4);
t440 = -Ifges(4,6) - t422;
t443 = Ifges(4,1) + Ifges(6,3);
t442 = Ifges(6,2) + Ifges(5,3);
t183 = sin(qJ(3));
t186 = cos(qJ(3));
t184 = sin(qJ(2));
t187 = cos(qJ(2));
t304 = qJD(1) * qJD(2);
t209 = qJDD(1) * t184 + t187 * t304;
t307 = qJD(3) * t184;
t290 = t183 * t307;
t305 = t186 * qJD(2);
t60 = qJD(1) * t290 - qJD(3) * t305 - t183 * qJDD(2) - t186 * t209;
t382 = -t60 / 0.2e1;
t314 = qJD(1) * t184;
t292 = t186 * t314;
t130 = qJD(2) * t183 + t292;
t309 = qJD(3) * t130;
t61 = -t186 * qJDD(2) + t183 * t209 + t309;
t379 = t61 / 0.2e1;
t293 = t183 * t314;
t129 = t293 - t305;
t376 = -t129 / 0.2e1;
t373 = t130 / 0.2e1;
t313 = qJD(1) * t187;
t153 = -qJD(3) + t313;
t372 = -t153 / 0.2e1;
t444 = -mrSges(5,1) - mrSges(4,3);
t347 = Ifges(6,6) * t183;
t352 = Ifges(4,4) * t183;
t406 = t443 * t186 + t347 - t352;
t346 = Ifges(6,6) * t186;
t349 = Ifges(5,6) * t186;
t404 = t442 * t183 + t346 - t349;
t439 = -Ifges(4,3) - Ifges(6,1) - Ifges(5,1);
t438 = -m(6) * pkin(4) - mrSges(6,1) + t444;
t392 = t440 * t183 - t441 * t186;
t243 = t186 * mrSges(5,2) - t183 * mrSges(5,3);
t245 = mrSges(4,1) * t186 - mrSges(4,2) * t183;
t341 = t183 * mrSges(6,2);
t437 = t245 + t341 - t243;
t169 = pkin(6) * t313;
t308 = qJD(3) * t183;
t436 = -qJD(4) * t183 - t169 + (-t183 * t313 + t308) * pkin(3);
t135 = qJDD(1) * t187 - t184 * t304;
t306 = qJD(3) * t186;
t311 = qJD(2) * t187;
t435 = t183 * t311 + t184 * t306;
t434 = -pkin(4) * t129 + qJD(5);
t175 = t184 * pkin(7);
t179 = t187 * pkin(2);
t295 = -pkin(1) - t179;
t219 = t295 - t175;
t117 = t219 * qJD(1);
t142 = qJD(2) * pkin(7) + t169;
t65 = -t186 * t117 + t142 * t183;
t218 = pkin(4) * t130 + t65;
t433 = t218 + qJD(4);
t167 = Ifges(3,4) * t313;
t123 = Ifges(4,4) * t129;
t348 = Ifges(6,6) * t129;
t419 = t130 * t443 + t153 * t421 - t123 + t348;
t121 = Ifges(6,6) * t130;
t343 = t130 * Ifges(5,6);
t420 = t129 * t442 + t153 * t422 + t121 - t343;
t432 = Ifges(3,1) * t314 + Ifges(3,5) * qJD(2) + t183 * t420 + t186 * t419 + t167;
t354 = Ifges(3,4) * t184;
t413 = t187 * Ifges(3,2);
t238 = t354 + t413;
t431 = t439 * t372 + t441 * t373 + t440 * t376 + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t238 / 0.2e1;
t127 = qJDD(3) - t135;
t125 = t135 * pkin(6);
t100 = qJDD(2) * pkin(7) + t125;
t333 = qJDD(1) * pkin(1);
t67 = -t135 * pkin(2) - pkin(7) * t209 - t333;
t14 = -t183 * t100 - t117 * t308 - t142 * t306 + t186 * t67;
t211 = qJDD(4) - t14;
t361 = pkin(3) + qJ(5);
t1 = -pkin(4) * t60 + qJD(5) * t153 - t127 * t361 + t211;
t380 = -t61 / 0.2e1;
t430 = Ifges(5,6) * t380 + mrSges(6,1) * t1 + (-Ifges(4,4) + Ifges(6,6)) * t379 + (Ifges(5,2) + t443) * t382 + (-Ifges(5,4) / 0.2e1 - t421 / 0.2e1) * t127;
t429 = -m(5) - m(6);
t428 = t442 * t61 + (Ifges(5,6) - Ifges(6,6)) * t60 - t422 * t127;
t427 = t135 / 0.2e1;
t426 = t209 / 0.2e1;
t424 = pkin(6) * (t184 * t305 + t187 * t308);
t171 = pkin(6) * t311;
t423 = -mrSges(3,3) + mrSges(2,2);
t334 = qJ(4) * t186;
t220 = qJ(5) * t183 - t334;
t210 = t220 * t187;
t418 = -qJD(1) * t210 + qJD(3) * t220 - qJD(5) * t186 + t436;
t356 = mrSges(4,3) * t129;
t80 = mrSges(4,2) * t153 - t356;
t360 = mrSges(5,1) * t129;
t83 = mrSges(5,3) * t153 + t360;
t417 = -t80 + t83;
t355 = mrSges(4,3) * t130;
t81 = -mrSges(4,1) * t153 - t355;
t359 = mrSges(5,1) * t130;
t85 = -mrSges(5,2) * t153 + t359;
t416 = -t81 + t85;
t358 = mrSges(6,1) * t129;
t84 = -mrSges(6,2) * t153 - t358;
t415 = -t83 + t84;
t414 = -qJ(4) * t306 + t313 * t334 + t436;
t250 = pkin(2) * t184 - pkin(7) * t187;
t131 = t250 * qJD(1);
t108 = t183 * t131;
t266 = pkin(6) * t186 - qJ(4);
t328 = t183 * t187;
t203 = -pkin(4) * t328 - t184 * t266;
t378 = pkin(4) + pkin(7);
t412 = -qJD(1) * t203 - t378 * t308 - t108;
t144 = t378 * t186;
t294 = -pkin(6) * t183 - pkin(3);
t324 = t186 * t187;
t193 = pkin(4) * t324 + (-qJ(5) + t294) * t184;
t330 = t131 * t186;
t411 = -qJD(1) * t193 + qJD(3) * t144 + t330;
t247 = mrSges(3,1) * t187 - mrSges(3,2) * t184;
t410 = -t247 - mrSges(2,1);
t409 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t129 - mrSges(4,2) * t130 - mrSges(3,3) * t314;
t408 = -t184 * t421 + t187 * t406;
t407 = -t184 * t422 + t187 * t404;
t350 = Ifges(5,6) * t183;
t405 = t186 * t442 - t347 + t350;
t126 = t209 * pkin(6);
t403 = t125 * t187 + t126 * t184;
t122 = Ifges(5,6) * t129;
t47 = -t153 * Ifges(5,4) - t130 * Ifges(5,2) + t122;
t344 = t130 * Ifges(4,4);
t48 = -t129 * Ifges(4,2) - t153 * Ifges(4,6) + t344;
t402 = t183 * t48 + t186 * t47;
t401 = t444 * t184;
t13 = t186 * t100 + t117 * t306 - t142 * t308 + t183 * t67;
t400 = t13 * t186 - t14 * t183;
t66 = t117 * t183 + t142 * t186;
t399 = -t66 - t434;
t4 = -qJ(4) * t127 + qJD(4) * t153 - t13;
t6 = -pkin(3) * t127 + t211;
t398 = t183 * t6 - t186 * t4;
t397 = mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t394 = -t184 * t439 + t187 * t392;
t168 = pkin(6) * t314;
t141 = -qJD(2) * pkin(2) + t168;
t217 = -qJ(4) * t130 + t141;
t24 = t129 * t361 + t217;
t337 = t186 * mrSges(6,2);
t241 = mrSges(6,3) * t183 - t337;
t242 = mrSges(5,2) * t183 + mrSges(5,3) * t186;
t244 = mrSges(4,1) * t183 + mrSges(4,2) * t186;
t42 = pkin(3) * t129 + t217;
t393 = -t141 * t244 - t24 * t241 + t42 * t242;
t391 = -t127 * t439 + t440 * t61 + t441 * t60;
t286 = m(6) * qJ(5) + mrSges(6,3);
t388 = pkin(7) * (-m(4) + t429);
t267 = -qJ(4) * t183 - pkin(2);
t136 = -pkin(3) * t186 + t267;
t246 = mrSges(3,1) * t184 + mrSges(3,2) * t187;
t387 = t246 + t438 * t187 + (-m(6) * t267 - (-m(6) * t361 - mrSges(6,3)) * t186 - m(5) * t136 + m(4) * pkin(2) + t437) * t184;
t386 = mrSges(4,1) - mrSges(5,2) + t286;
t3 = -pkin(4) * t61 + qJDD(5) - t4;
t385 = -t14 * mrSges(4,1) + t13 * mrSges(4,2) - t6 * mrSges(5,2) - t3 * mrSges(6,2) + t4 * mrSges(5,3) + t1 * mrSges(6,3);
t381 = t60 / 0.2e1;
t377 = t127 / 0.2e1;
t375 = t129 / 0.2e1;
t374 = -t130 / 0.2e1;
t371 = t153 / 0.2e1;
t366 = pkin(4) * t184;
t176 = t184 * pkin(6);
t357 = mrSges(6,1) * t130;
t353 = Ifges(3,4) * t187;
t351 = Ifges(4,4) * t186;
t335 = qJ(4) * t129;
t329 = t183 * t184;
t327 = t184 * t186;
t188 = cos(qJ(1));
t326 = t184 * t188;
t185 = sin(qJ(1));
t325 = t185 * t187;
t323 = t187 * t188;
t322 = t188 * t183;
t134 = t250 * qJD(2);
t318 = t179 + t175;
t137 = -pkin(1) - t318;
t321 = t183 * t134 + t137 * t306;
t319 = -pkin(3) * t329 + qJ(4) * t327;
t156 = pkin(6) * t324;
t88 = t183 * t137 + t156;
t317 = t188 * pkin(1) + t185 * pkin(6);
t312 = qJD(2) * t184;
t155 = pkin(6) * t328;
t298 = pkin(6) * t312;
t93 = t176 - t319;
t271 = -t308 / 0.2e1;
t270 = t308 / 0.2e1;
t269 = -t306 / 0.2e1;
t268 = t306 / 0.2e1;
t30 = -t60 * mrSges(5,1) + t127 * mrSges(5,2);
t29 = -t61 * mrSges(6,1) + t127 * mrSges(6,2);
t27 = -t60 * mrSges(6,1) - t127 * mrSges(6,3);
t265 = t304 / 0.2e1;
t87 = t137 * t186 - t155;
t261 = pkin(3) * t435 + qJ(4) * t290 + t171;
t260 = pkin(3) * t324 + qJ(4) * t328 + t318;
t259 = pkin(2) * t323 + pkin(7) * t326 + t317;
t251 = t294 * t184;
t101 = -qJDD(2) * pkin(2) + t126;
t76 = qJ(4) * t187 - t88;
t249 = qJD(3) * t156 - t134 * t186 + t137 * t308;
t248 = -qJD(4) * t187 + t321;
t239 = Ifges(4,1) * t183 + t351;
t237 = -Ifges(4,2) * t183 + t351;
t236 = Ifges(4,2) * t186 + t352;
t234 = Ifges(5,4) * t183 + Ifges(5,5) * t186;
t233 = Ifges(6,4) * t186 - Ifges(6,5) * t183;
t231 = Ifges(3,5) * t187 - Ifges(3,6) * t184;
t229 = Ifges(4,5) * t183 + Ifges(4,6) * t186;
t228 = -Ifges(5,2) * t186 + t350;
t227 = Ifges(5,2) * t183 + t349;
t222 = -Ifges(6,3) * t183 + t346;
t216 = pkin(1) * t246;
t212 = t184 * (Ifges(3,1) * t187 - t354);
t40 = qJ(4) * t153 - t66;
t104 = t183 * t325 + t186 * t188;
t106 = -t185 * t186 + t187 * t322;
t204 = -g(1) * t106 - g(2) * t104 - g(3) * t329;
t200 = Ifges(5,4) * t184 + t187 * t228;
t195 = Ifges(4,6) * t184 + t187 * t237;
t192 = qJ(4) * t60 - qJD(4) * t130 + t101;
t180 = t188 * pkin(6);
t178 = t187 * pkin(3);
t143 = t378 * t183;
t139 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t313;
t118 = t244 * t184;
t116 = -t186 * t361 + t267;
t107 = t185 * t183 + t186 * t323;
t105 = t185 * t324 - t322;
t82 = mrSges(6,3) * t153 + t357;
t79 = -pkin(6) * t292 + t108;
t78 = pkin(6) * t293 + t330;
t77 = t178 - t87;
t75 = qJ(5) * t329 + t93;
t73 = -mrSges(5,2) * t129 - mrSges(5,3) * t130;
t71 = pkin(3) * t130 + t335;
t70 = -mrSges(6,2) * t130 + mrSges(6,3) * t129;
t69 = qJD(1) * t251 - t330;
t68 = t266 * t314 - t108;
t63 = -pkin(4) * t329 - t76;
t59 = qJ(5) * t187 + t155 + t178 + (-t137 + t366) * t186;
t39 = pkin(3) * t153 + qJD(4) + t65;
t37 = t130 * t361 + t335;
t34 = t183 * t298 - t249;
t33 = t321 - t424;
t32 = (-qJ(4) * t311 - qJD(4) * t184) * t186 + t261;
t31 = qJD(2) * t251 + t249;
t28 = mrSges(5,1) * t61 - mrSges(5,3) * t127;
t26 = -mrSges(4,2) * t127 - mrSges(4,3) * t61;
t25 = mrSges(4,1) * t127 + mrSges(4,3) * t60;
t23 = -t40 + t434;
t22 = -qJ(4) * t312 - t248 + t424;
t21 = t153 * t361 + t433;
t20 = qJD(2) * t210 + (qJD(5) * t183 + (qJ(5) * qJD(3) - qJD(4)) * t186) * t184 + t261;
t19 = (-pkin(4) * t327 - t155) * qJD(3) + t203 * qJD(2) + t248;
t18 = mrSges(6,2) * t60 + mrSges(6,3) * t61;
t17 = mrSges(4,1) * t61 - mrSges(4,2) * t60;
t16 = -mrSges(5,2) * t61 + mrSges(5,3) * t60;
t15 = -pkin(4) * t290 + qJD(2) * t193 + qJD(5) * t187 + t249;
t11 = -t60 * Ifges(4,4) - t61 * Ifges(4,2) + t127 * Ifges(4,6);
t5 = pkin(3) * t61 + t192;
t2 = qJD(5) * t129 + t361 * t61 + t192;
t7 = [(t39 * mrSges(5,1) + t21 * mrSges(6,1) + mrSges(4,2) * t141 - mrSges(6,2) * t24 + t65 * mrSges(4,3) - mrSges(5,3) * t42) * (t187 * t305 - t290) + t247 * t333 + t101 * t118 + t87 * t25 + t88 * t26 + t93 * t16 + t75 * t18 + t76 * t28 + t77 * t30 + t33 * t80 + t34 * t81 + t15 * t82 + t22 * t83 + t19 * t84 + t31 * t85 + t20 * t70 + t32 * t73 - t139 * t298 + (mrSges(4,1) * t141 + t40 * mrSges(5,1) - t23 * mrSges(6,1) - mrSges(5,2) * t42 - t66 * mrSges(4,3) + mrSges(6,3) * t24) * t435 + t63 * t29 + t59 * t27 + (qJD(2) * t200 + t227 * t307) * t374 + (-t11 / 0.2e1 - t13 * mrSges(4,3) - t3 * mrSges(6,1) + t4 * mrSges(5,1) + t428 / 0.2e1) * t329 - t216 * t304 + ((t222 - t239) * t307 + t408 * qJD(2)) * t373 - t409 * t171 + (t432 / 0.2e1 - t402 / 0.2e1) * t311 + m(5) * (t22 * t40 + t31 * t39 + t32 * t42 + t4 * t76 + t5 * t93 + t6 * t77) + m(6) * (t1 * t59 + t15 * t21 + t19 * t23 + t2 * t75 + t20 * t24 + t3 * t63) + (qJD(2) * t195 - t236 * t307) * t376 + (t176 * t209 + t403) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t403 * pkin(6)) + (qJD(2) * t407 + t307 * t405) * t375 + (t429 * (-t105 * pkin(3) - qJ(4) * t104 + t180) + t423 * t188 + (-m(4) - m(3)) * t180 + t386 * t105 + t397 * t104 + (-m(6) * t295 + m(3) * pkin(1) + (-m(4) - m(5)) * t219 - t401 - t410) * t185) * g(1) + (t419 * t271 + t420 * t268 + m(4) * t101 * pkin(6) + Ifges(3,1) * t209 + Ifges(3,4) * t427 + Ifges(3,5) * qJDD(2) + t2 * t241 + t228 * t381 + t237 * t380 - t5 * t242 - t265 * t413 + t48 * t269 + t47 * t270 + t392 * t377 + t404 * t379 + t406 * t382 + (m(6) * t378 + mrSges(6,1)) * t185 * g(1)) * t184 + (-qJDD(2) * mrSges(3,1) + t17) * t176 + (-t65 * mrSges(4,1) - t66 * mrSges(4,2) + t39 * mrSges(5,2) + t23 * mrSges(6,2) - t40 * mrSges(5,3) - t21 * mrSges(6,3) - t431) * t312 + ((-t229 + t233 + t234) * t307 + t394 * qJD(2)) * t372 + (mrSges(5,1) * t6 - mrSges(4,3) * t14 + t430) * t327 + m(4) * (t13 * t88 + t14 * t87 + t141 * t171 + t33 * t66 - t34 * t65) - pkin(1) * (-t135 * mrSges(3,1) + mrSges(3,2) * t209) + t212 * t265 + t353 * t426 + t238 * t427 + Ifges(2,3) * qJDD(1) + (-m(3) * t317 - m(4) * t259 + t429 * (t107 * pkin(3) + t106 * qJ(4) + t259) + t410 * t188 + t423 * t185 + t438 * t326 - t386 * t107 - t397 * t106) * g(2) + (-t391 / 0.2e1 + pkin(6) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t135) + t353 * t265 - Ifges(4,6) * t380 - Ifges(5,4) * t381 + Ifges(3,4) * t426 + Ifges(3,2) * t427 + Ifges(3,6) * qJDD(2) + t421 * t382 + t422 * t379 + t439 * t377 + t385) * t187 + qJD(2) ^ 2 * t231 / 0.2e1; t186 * t11 / 0.2e1 + t414 * t73 + t418 * t70 + (t1 * t143 + t116 * t2 + t144 * t3 + t21 * t411 + t23 * t412 + t24 * t418) * m(6) + t419 * t268 + t420 * t270 + Ifges(3,6) * t135 + t136 * t16 + t143 * t27 + t144 * t29 - t125 * mrSges(3,2) - t126 * mrSges(3,1) + t116 * t18 - t79 * t80 - t78 * t81 - t68 * t83 - t69 * t85 + (t65 * (mrSges(4,1) * t184 - mrSges(4,3) * t324) - t21 * (mrSges(6,1) * t324 - mrSges(6,3) * t184) - t39 * (mrSges(5,1) * t324 + mrSges(5,2) * t184) - t66 * (-mrSges(4,2) * t184 - mrSges(4,3) * t328) - t23 * (-mrSges(6,1) * t328 + mrSges(6,2) * t184) - t40 * (mrSges(5,1) * t328 - mrSges(5,3) * t184) + t394 * t371 + (-t408 / 0.2e1 + t200 / 0.2e1) * t130 + (-t407 / 0.2e1 + t195 / 0.2e1) * t129 + (t216 - t212 / 0.2e1) * qJD(1)) * qJD(1) - (-Ifges(3,2) * t314 + t167 + t432) * t313 / 0.2e1 + (t392 * t372 - t393) * qJD(3) - pkin(2) * t17 + (-t233 / 0.2e1 - t234 / 0.2e1) * t127 - t231 * t304 / 0.2e1 + t409 * t169 + t229 * t377 + t222 * t381 + (t306 * t39 + t308 * t40 + t398) * mrSges(5,1) + t139 * t168 + (t239 + t227) * t382 + (t306 * t65 - t308 * t66 + t400) * mrSges(4,3) + (t236 + t405) * t380 + (t136 * t5 - t39 * t69 - t40 * t68 + t414 * t42) * m(5) + (-pkin(2) * t101 - t141 * t169 + t65 * t78 - t66 * t79) * m(4) + (((t183 * t40 + t186 * t39) * qJD(3) + t398) * m(5) + (-t28 + t26) * t186 + (-t25 + t30) * t183 + t416 * t306 + t417 * t308 + ((-t183 * t66 + t186 * t65) * qJD(3) + t400) * m(4)) * pkin(7) + t431 * t314 + (t185 * t387 + t325 * t388) * g(2) + (t188 * t387 + t323 * t388) * g(1) + (-g(3) * t184 + t186 * t3 + t21 * t306 - t23 * t308) * mrSges(6,1) + t430 * t183 + t2 * (-t186 * mrSges(6,3) - t341) - t428 * t186 / 0.2e1 + (-t237 / 0.2e1 + t404 / 0.2e1) * qJD(3) * t129 + Ifges(3,5) * t209 + t47 * t269 + t48 * t271 + (-t247 - m(6) * (t260 + t366) - m(4) * t318 - m(5) * t260 + (-t186 * t286 - t437) * t187 + t401) * g(3) + (-t228 / 0.2e1 + t406 / 0.2e1) * t309 + t5 * t243 - t101 * t245 + (t402 / 0.2e1 + t393) * t313 + t411 * t82 + t412 * t84 + Ifges(3,3) * qJDD(2); t218 * t84 - t361 * t27 + t415 * qJD(4) + (-m(5) * t39 + t355 - t416) * t66 + (-m(5) * t40 + t356 - t417) * t65 + (-Ifges(4,2) * t130 - t123 + t419) * t375 + (t441 * t129 + t440 * t130) * t371 - t141 * (mrSges(4,1) * t130 - mrSges(4,2) * t129) - t24 * (mrSges(6,2) * t129 + mrSges(6,3) * t130) - t42 * (-mrSges(5,2) * t130 + mrSges(5,3) * t129) - t37 * t70 - t71 * t73 + t391 - t385 + (qJ(4) * t3 - t1 * t361 + t399 * t21 + t23 * t433 - t24 * t37) * m(6) + (Ifges(5,2) * t129 + t343 + t48) * t373 - pkin(3) * t30 + t23 * t357 + t21 * t358 + t39 * t360 - t40 * t359 + (-pkin(3) * t6 - qJ(4) * t4 - qJD(4) * t40 - t42 * t71) * m(5) + t399 * t82 + (t118 + t429 * t319 + (t183 * t286 - t242 - t337) * t184) * g(3) + (t429 * (-t106 * pkin(3) + qJ(4) * t107) - t397 * t107 + t386 * t106) * g(1) + (t429 * (-t104 * pkin(3) + qJ(4) * t105) - t397 * t105 + t386 * t104) * g(2) + (t29 - t28) * qJ(4) + (t130 * t442 + t122 - t348 + t47) * t376 + (-t129 * t443 + t121 - t344 + t420) * t374; t415 * t153 + (t70 + t73) * t130 + t27 + t30 + (t130 * t24 + t153 * t23 + t1 + t204) * m(6) + (t130 * t42 - t153 * t40 + t204 + t6) * m(5); -t129 * t70 - t153 * t82 + (-g(1) * t107 - g(2) * t105 - g(3) * t327 - t129 * t24 - t153 * t21 + t3) * m(6) + t29;];
tau = t7;
