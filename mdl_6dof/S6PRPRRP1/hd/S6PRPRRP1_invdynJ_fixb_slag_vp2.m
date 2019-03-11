% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:27
% EndTime: 2019-03-08 19:55:49
% DurationCPUTime: 15.72s
% Computational Cost: add. (4931->571), mult. (11489->766), div. (0->0), fcn. (8874->12), ass. (0->265)
t445 = Ifges(6,4) + Ifges(7,4);
t446 = Ifges(6,1) + Ifges(7,1);
t423 = -Ifges(7,5) - Ifges(6,5);
t444 = Ifges(6,2) + Ifges(7,2);
t443 = Ifges(6,6) + Ifges(7,6);
t199 = cos(qJ(5));
t187 = pkin(5) * t199 + pkin(4);
t453 = m(7) * t187;
t191 = sin(pkin(6));
t189 = sin(pkin(11));
t192 = cos(pkin(11));
t198 = sin(qJ(2));
t201 = cos(qJ(2));
t227 = t189 * t201 + t192 * t198;
t128 = t227 * t191;
t119 = qJD(1) * t128;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t250 = pkin(4) * t197 - pkin(9) * t200;
t452 = t250 * qJD(4) - t119;
t451 = t445 * t199;
t196 = sin(qJ(5));
t450 = t445 * t196;
t195 = -qJ(6) - pkin(9);
t449 = -m(7) * t195 + mrSges(6,3) + mrSges(7,3);
t308 = qJD(4) * t199;
t312 = qJD(2) * t197;
t155 = -t196 * t312 + t308;
t302 = qJD(2) * qJD(4);
t164 = qJDD(2) * t197 + t200 * t302;
t96 = qJD(5) * t155 + qJDD(4) * t196 + t164 * t199;
t380 = t96 / 0.2e1;
t310 = qJD(4) * t196;
t156 = t199 * t312 + t310;
t97 = -qJD(5) * t156 + qJDD(4) * t199 - t164 * t196;
t379 = t97 / 0.2e1;
t448 = -t155 / 0.2e1;
t374 = t156 / 0.2e1;
t311 = qJD(2) * t200;
t432 = t311 - qJD(5);
t372 = -t432 / 0.2e1;
t447 = qJD(4) / 0.2e1;
t442 = -Ifges(7,3) - Ifges(6,3);
t400 = -t196 * t443 - t199 * t423;
t398 = -t196 * t444 + t451;
t396 = t199 * t446 - t450;
t163 = qJDD(2) * t200 - t197 * t302;
t152 = qJDD(5) - t163;
t441 = t445 * t379 + t446 * t380 - t423 * t152 / 0.2e1;
t330 = t191 * t198;
t285 = qJD(1) * t330;
t167 = t189 * t285;
t313 = qJD(1) * t201;
t284 = t191 * t313;
t122 = t192 * t284 - t167;
t185 = pkin(2) * t189 + pkin(8);
t309 = qJD(4) * t197;
t282 = t185 * t309;
t322 = t196 * t200;
t440 = t122 * t322 + t196 * t282 + t199 * t452;
t251 = pkin(4) * t200 + pkin(9) * t197;
t226 = -pkin(3) - t251;
t368 = pkin(2) * t192;
t148 = t226 - t368;
t304 = qJD(5) * t199;
t320 = t199 * t200;
t439 = -t122 * t320 + t148 * t304 + t196 * t452;
t188 = Ifges(5,4) * t311;
t437 = t445 * t155;
t413 = t446 * t156 + t423 * t432 + t437;
t438 = Ifges(5,1) * t312 + Ifges(5,5) * qJD(4) + t199 * t413 + t188;
t190 = sin(pkin(10));
t193 = cos(pkin(10));
t194 = cos(pkin(6));
t324 = t194 * t201;
t224 = -t190 * t324 - t193 * t198;
t214 = t224 * pkin(2);
t319 = t201 * t192;
t153 = t189 * t198 - t319;
t218 = t153 * t194;
t82 = t190 * t218 - t193 * t227;
t325 = t194 * t198;
t314 = -t189 * t324 - t192 * t325;
t83 = -t193 * t153 + t190 * t314;
t436 = t82 * pkin(3) + pkin(8) * t83 + t214;
t165 = qJD(2) * pkin(2) + t284;
t111 = t189 * t165 + t192 * t285;
t107 = qJD(2) * pkin(8) + t111;
t180 = qJD(1) * t194 + qJD(3);
t335 = t180 * t197;
t72 = t107 * t200 + t335;
t67 = qJD(4) * pkin(9) + t72;
t110 = t165 * t192 - t167;
t77 = qJD(2) * t226 - t110;
t20 = -t196 * t67 + t199 * t77;
t21 = t196 * t77 + t199 * t67;
t179 = qJDD(1) * t194 + qJDD(3);
t307 = qJD(4) * t200;
t260 = qJD(2) * t285;
t328 = t191 * t201;
t132 = qJDD(1) * t328 - t260;
t338 = qJDD(2) * pkin(2);
t125 = t132 + t338;
t276 = qJD(2) * t313;
t133 = (qJDD(1) * t198 + t276) * t191;
t69 = t189 * t125 + t192 * t133;
t62 = qJDD(2) * pkin(8) + t69;
t12 = -t107 * t309 + t197 * t179 + t180 * t307 + t200 * t62;
t10 = qJDD(4) * pkin(9) + t12;
t306 = qJD(5) * t196;
t68 = t125 * t192 - t189 * t133;
t61 = -qJDD(2) * pkin(3) - t68;
t35 = -pkin(4) * t163 - pkin(9) * t164 + t61;
t3 = t199 * t10 + t196 * t35 + t77 * t304 - t306 * t67;
t4 = -qJD(5) * t21 - t10 * t196 + t199 * t35;
t249 = -t196 * t4 + t199 * t3;
t435 = -t20 * t304 - t21 * t306 + t249;
t434 = -t190 * t198 + t193 * t324;
t433 = t445 * t156;
t247 = mrSges(5,1) * t200 - mrSges(5,2) * t197;
t431 = m(6) * t251 + t197 * t449 + t200 * t453 + t247;
t353 = Ifges(5,4) * t197;
t238 = Ifges(5,2) * t200 + t353;
t430 = Ifges(5,6) * t447 + qJD(2) * t238 / 0.2e1 + t442 * t372 + t423 * t374 + t443 * t448;
t381 = m(7) * pkin(5);
t429 = -m(6) - m(5);
t428 = t163 / 0.2e1;
t427 = t164 / 0.2e1;
t426 = -mrSges(7,1) - mrSges(6,1);
t425 = mrSges(6,2) + mrSges(7,2);
t424 = mrSges(5,3) - mrSges(4,2);
t158 = t185 * t320;
t225 = pkin(5) * t197 - qJ(6) * t320;
t303 = qJD(6) * t199;
t421 = -t197 * t303 + t225 * qJD(4) + (-t158 + (qJ(6) * t197 - t148) * t196) * qJD(5) + t440;
t420 = t152 * t443 + t444 * t97 + t445 * t96;
t321 = t197 * t199;
t418 = (-qJ(6) * qJD(5) - qJD(4) * t185) * t321 + (-qJD(6) * t197 + (-qJ(6) * qJD(4) - qJD(5) * t185) * t200) * t196 + t439;
t417 = (-t197 * t308 - t200 * t306) * t185 + t439;
t101 = t196 * t148 + t158;
t416 = -qJD(5) * t101 + t440;
t414 = t155 * t444 - t432 * t443 + t433;
t37 = -mrSges(6,1) * t97 + mrSges(6,2) * t96;
t412 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t164 + t37;
t262 = qJD(5) * t195;
t283 = t196 * t311;
t159 = t250 * qJD(2);
t71 = -t197 * t107 + t180 * t200;
t40 = t196 * t159 + t199 * t71;
t411 = qJ(6) * t283 + t196 * t262 + t303 - t40;
t39 = t199 * t159 - t196 * t71;
t410 = -qJD(2) * t225 - qJD(6) * t196 + t199 * t262 - t39;
t366 = pkin(5) * t196;
t409 = pkin(5) * t306 - t335 - (qJD(2) * t366 + t107) * t200;
t408 = t381 + mrSges(7,1);
t294 = mrSges(5,3) * t312;
t407 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t155 + mrSges(6,2) * t156 + t294;
t406 = -t197 * t442 + t200 * t400;
t405 = t197 * t443 + t200 * t398;
t404 = -t197 * t423 + t200 * t396;
t127 = t189 * t330 - t191 * t319;
t403 = (mrSges(3,1) * t201 - mrSges(3,2) * t198) * t191 - mrSges(4,1) * t127 - mrSges(4,2) * t128;
t401 = -t196 * t423 + t199 * t443;
t399 = t199 * t444 + t450;
t397 = t196 * t446 + t451;
t394 = -t152 * t442 - t423 * t96 + t443 * t97;
t13 = -t107 * t307 + t179 * t200 - t180 * t309 - t197 * t62;
t393 = t12 * t200 - t13 * t197;
t392 = -m(7) + t429;
t391 = -mrSges(6,1) - t408;
t244 = -mrSges(7,1) * t199 + mrSges(7,2) * t196;
t246 = -mrSges(6,1) * t199 + mrSges(6,2) * t196;
t390 = m(6) * pkin(4) + mrSges(5,1) - t244 - t246 + t453;
t389 = -m(6) * pkin(9) + mrSges(5,2) - t449;
t102 = -mrSges(7,1) * t155 + mrSges(7,2) * t156;
t66 = -qJD(4) * pkin(4) - t71;
t41 = -pkin(5) * t155 + qJD(6) + t66;
t388 = -m(7) * t41 - t102;
t78 = t190 * t153 + t193 * t314;
t386 = -mrSges(4,1) - t431;
t1 = pkin(5) * t152 - qJ(6) * t96 - qJD(6) * t156 + t4;
t15 = qJ(6) * t155 + t21;
t2 = qJ(6) * t97 + qJD(6) * t155 + t3;
t14 = -qJ(6) * t156 + t20;
t8 = -pkin(5) * t432 + t14;
t385 = -t1 * t196 - t15 * t306 + t199 * t2 - t8 * t304;
t384 = -m(6) * t66 - t407;
t383 = -t4 * mrSges(6,1) + t3 * mrSges(6,2) + t2 * mrSges(7,2);
t382 = qJD(2) ^ 2;
t378 = t152 / 0.2e1;
t376 = t155 / 0.2e1;
t54 = mrSges(7,1) * t152 - mrSges(7,3) * t96;
t55 = mrSges(6,1) * t152 - mrSges(6,3) * t96;
t359 = t54 + t55;
t56 = -mrSges(7,2) * t152 + mrSges(7,3) * t97;
t57 = -mrSges(6,2) * t152 + mrSges(6,3) * t97;
t358 = t56 + t57;
t357 = mrSges(6,3) * t155;
t356 = mrSges(6,3) * t156;
t355 = mrSges(7,3) * t155;
t354 = mrSges(7,3) * t156;
t352 = Ifges(5,4) * t200;
t11 = -qJDD(4) * pkin(4) - t13;
t347 = t11 * t197;
t342 = t196 * t78;
t341 = t196 * t83;
t336 = t128 * t196;
t331 = t191 * t197;
t329 = t191 * t200;
t323 = t196 * t197;
t113 = mrSges(7,2) * t432 + t355;
t114 = mrSges(6,2) * t432 + t357;
t318 = t113 + t114;
t115 = -mrSges(7,1) * t432 - t354;
t116 = -mrSges(6,1) * t432 - t356;
t317 = t115 + t116;
t305 = qJD(5) * t197;
t181 = pkin(2) * t328;
t36 = -t97 * mrSges(7,1) + t96 * mrSges(7,2);
t301 = t36 + t412;
t296 = m(4) - t392;
t293 = mrSges(5,3) * t311;
t287 = t102 + t407;
t281 = t185 * t307;
t280 = t196 * t307;
t261 = t302 / 0.2e1;
t253 = t434 * pkin(2);
t245 = mrSges(6,1) * t196 + mrSges(6,2) * t199;
t243 = mrSges(7,1) * t196 + mrSges(7,2) * t199;
t233 = Ifges(5,5) * t200 - Ifges(5,6) * t197;
t105 = t128 * t200 + t194 * t197;
t47 = t105 * t199 + t127 * t196;
t46 = -t105 * t196 + t127 * t199;
t104 = t128 * t197 - t194 * t200;
t106 = -qJD(2) * pkin(3) - t110;
t221 = t106 * (mrSges(5,1) * t197 + mrSges(5,2) * t200);
t220 = t197 * (Ifges(5,1) * t200 - t353);
t213 = -t196 * t305 + t199 * t307;
t212 = t197 * t304 + t280;
t186 = -pkin(3) - t368;
t172 = t195 * t199;
t171 = t195 * t196;
t170 = -qJD(4) * mrSges(5,2) + t293;
t157 = t247 * qJD(2);
t138 = (t185 + t366) * t197;
t134 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t163;
t131 = t199 * t148;
t121 = t153 * t191 * qJD(2);
t120 = qJD(2) * t128;
t109 = -mrSges(5,1) * t163 + mrSges(5,2) * t164;
t108 = pkin(5) * t212 + t281;
t100 = -t185 * t322 + t131;
t90 = -qJ(6) * t323 + t101;
t79 = -t190 * t227 - t193 * t218;
t76 = -qJ(6) * t321 + t131 + (-t185 * t196 - pkin(5)) * t200;
t53 = t190 * t331 + t200 * t83;
t52 = -t190 * t329 + t197 * t83;
t51 = -t193 * t331 - t200 * t78;
t50 = t193 * t329 - t197 * t78;
t45 = -qJD(4) * t104 - t121 * t200;
t44 = qJD(4) * t105 - t121 * t197;
t7 = qJD(5) * t46 + t120 * t196 + t199 * t45;
t6 = -qJD(5) * t47 + t120 * t199 - t196 * t45;
t5 = -pkin(5) * t97 + qJDD(6) + t11;
t9 = [m(2) * qJDD(1) + t105 * t134 + t127 * t109 - t120 * t157 + t45 * t170 + t318 * t7 + t317 * t6 + t358 * t47 + t359 * t46 + (-mrSges(3,1) * t198 - mrSges(3,2) * t201) * t382 * t191 + (-mrSges(4,1) * t120 + mrSges(4,2) * t121) * qJD(2) + t287 * t44 + t301 * t104 + t403 * qJDD(2) + (-m(2) - m(3) - t296) * g(3) + m(6) * (t104 * t11 + t20 * t6 + t21 * t7 + t3 * t47 + t4 * t46 + t44 * t66) + m(7) * (t1 * t46 + t104 * t5 + t15 * t7 + t2 * t47 + t41 * t44 + t6 * t8) + m(5) * (-t104 * t13 + t105 * t12 + t106 * t120 + t127 * t61 - t44 * t71 + t45 * t72) + m(4) * (-t110 * t120 - t111 * t121 - t127 * t68 + t128 * t69 + t179 * t194) + m(3) * (qJDD(1) * t194 ^ 2 + (t132 * t201 + t133 * t198) * t191); (qJD(2) * t119 + t192 * t338 + t68) * mrSges(4,1) + (-t20 * t213 - t21 * t212 - t3 * t323 - t321 * t4) * mrSges(6,3) + (qJD(2) * t122 - t189 * t338 - t69) * mrSges(4,2) + (t20 * mrSges(6,1) + t8 * mrSges(7,1) - t21 * mrSges(6,2) - t15 * mrSges(7,2) - t72 * mrSges(5,3) - t430) * t309 + t238 * t428 + (t191 * t276 - t133) * mrSges(3,2) + (t384 + t388) * t122 * t197 + (-t224 * mrSges(3,1) - (t190 * t325 - t193 * t201) * mrSges(3,2) - m(4) * t214 - m(7) * (pkin(5) * t341 + t436) - t424 * t83 + t426 * (t320 * t82 + t341) - t425 * (t199 * t83 - t322 * t82) + t429 * t436 + t386 * t82) * g(1) + (t233 * t447 + t406 * t372 + t404 * t374 + t405 * t376 + t221) * qJD(4) + (-t372 * t401 - t374 * t397 - t376 * t399) * t305 + (-t1 * t321 - t15 * t212 - t2 * t323 - t213 * t8) * mrSges(7,3) + (t132 + t260) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (-t434 * mrSges(3,1) - (-t190 * t201 - t193 * t325) * mrSges(3,2) + t342 * t381 - m(4) * t253 + t424 * t78 + t426 * (t320 * t79 - t342) - t425 * (-t199 * t78 - t322 * t79) + t392 * (t79 * pkin(3) - pkin(8) * t78 + t253) + t386 * t79) * g(2) + t438 * t307 / 0.2e1 + (-t307 * t71 + t393) * mrSges(5,3) + t407 * t281 - t414 * t280 / 0.2e1 + t352 * t427 + t245 * t347 + t41 * (mrSges(7,1) * t212 + mrSges(7,2) * t213) + t66 * (mrSges(6,1) * t212 + mrSges(6,2) * t213) + t321 * t441 + t412 * t185 * t197 + (-t106 * t119 - (-t197 * t71 + t200 * t72) * t122 + t186 * t61 + ((-t197 * t72 - t200 * t71) * qJD(4) + t393) * t185) * m(5) - t394 * t200 / 0.2e1 + (Ifges(5,1) * t164 + Ifges(5,4) * t428 + t5 * t243 + t378 * t400 + t379 * t398 + t380 * t396) * t197 - (t196 * t413 + t199 * t414) * t305 / 0.2e1 + t416 * t116 + t417 * t114 + (t100 * t4 + t101 * t3 + (t307 * t66 + t347) * t185 + t417 * t21 + t416 * t20) * m(6) + t418 * t113 - t420 * t323 / 0.2e1 + (t1 * t76 + t108 * t41 + t138 * t5 + t15 * t418 + t2 * t90 + t421 * t8) * m(7) + t421 * t115 - t170 * t282 + ((t189 * t69 + t192 * t68) * pkin(2) + t110 * t119 - t111 * t122) * m(4) + qJDD(4) * (Ifges(5,5) * t197 + Ifges(5,6) * t200) + t186 * t109 + t119 * t157 + t138 * t36 + (Ifges(5,4) * t427 + Ifges(5,2) * t428 - t1 * mrSges(7,1) + (-Ifges(5,2) * t197 + t352) * t261 - t122 * t170 + t185 * t134 + t423 * t380 - t443 * t379 + t442 * t378 + t383) * t200 + t220 * t261 - t61 * t247 + t76 * t54 + (-t336 * t381 - m(4) * t181 - t128 * mrSges(5,3) + t426 * (-t127 * t320 + t336) - t425 * (t127 * t322 + t128 * t199) + t392 * (-t127 * pkin(3) + pkin(8) * t128 + t181) + t431 * t127 - t403) * g(3) + t90 * t56 + t100 * t55 + t101 * t57 + t108 * t102; m(4) * t179 + ((-t196 * t317 + t199 * t318 + t170) * qJD(4) + m(5) * (qJD(4) * t72 + t13) + m(6) * (-t20 * t310 + t21 * t308 - t11) + m(7) * (t15 * t308 - t310 * t8 - t5) - t301) * t200 + (t134 + t358 * t199 - t359 * t196 + t287 * qJD(4) + (-t196 * t318 - t199 * t317) * qJD(5) + m(5) * (-qJD(4) * t71 + t12) + m(6) * (qJD(4) * t66 + t435) + m(7) * (qJD(4) * t41 + t385)) * t197 + (-t194 * g(3) + (-g(1) * t190 + g(2) * t193) * t191) * t296; (-t116 * t304 - t114 * t306 - t196 * t55 + m(6) * ((-t21 * t196 - t20 * t199) * qJD(5) + t249) + t199 * t57) * pkin(9) + t409 * t102 + t410 * t115 + (t1 * t171 + t15 * t411 - t172 * t2 - t187 * t5 + t409 * t41 + t410 * t8) * m(7) + t411 * t113 + t430 * t312 + t385 * mrSges(7,3) + (t283 / 0.2e1 - t306 / 0.2e1) * t414 + (t389 * t51 + t390 * t50) * g(2) + (t104 * t390 + t105 * t389) * g(3) + (t389 * t53 + t390 * t52) * g(1) + t397 * t380 + t399 * t379 + t401 * t378 + (-pkin(4) * t11 - t20 * t39 - t21 * t40) * m(6) - t233 * t302 / 0.2e1 + t435 * mrSges(6,3) - (-Ifges(5,2) * t312 + t188 + t438) * t311 / 0.2e1 + t413 * t304 / 0.2e1 + (t293 - t170) * t71 + (-t8 * (mrSges(7,1) * t197 - mrSges(7,3) * t320) - t20 * (mrSges(6,1) * t197 - mrSges(6,3) * t320) - t15 * (-mrSges(7,2) * t197 - mrSges(7,3) * t322) - t21 * (-mrSges(6,2) * t197 - mrSges(6,3) * t322) - t221) * qJD(2) + (t294 + t384) * t72 + t196 * t441 - t382 * t220 / 0.2e1 + t420 * t199 / 0.2e1 + Ifges(5,3) * qJDD(4) - t187 * t36 + t171 * t54 - t172 * t56 + Ifges(5,6) * t163 + Ifges(5,5) * t164 - t12 * mrSges(5,2) + t13 * mrSges(5,1) + (t155 * t398 + t156 * t396 - t400 * t432) * qJD(5) / 0.2e1 - (t155 * t405 + t156 * t404 - t406 * t432) * qJD(2) / 0.2e1 - pkin(4) * t37 + t5 * t244 + t11 * t246 + t432 * (-t41 * t243 - t66 * t245) - t40 * t114 - t39 * t116; (t354 - m(7) * (t14 - t8) + t115) * t15 + t414 * t374 - t383 + (t357 - t114) * t20 + t8 * t355 + t408 * t1 + (t391 * t46 + t425 * t47) * g(3) + (-t155 * t423 - t156 * t443) * t432 / 0.2e1 - (t446 * t155 - t433) * t156 / 0.2e1 + (-t156 * t444 + t413 + t437) * t448 + (-t425 * (t196 * t79 - t199 * t51) + t391 * (-t196 * t51 - t199 * t79)) * g(2) + (-t425 * (t196 * t82 - t199 * t53) + t391 * (-t196 * t53 - t199 * t82)) * g(1) - t41 * (mrSges(7,1) * t156 + mrSges(7,2) * t155) - t66 * (mrSges(6,1) * t156 + mrSges(6,2) * t155) + (t356 + t116) * t21 - t14 * t113 + t394 + (t156 * t388 + t54) * pkin(5); -t155 * t113 + t156 * t115 + (-g(1) * t52 - g(2) * t50 - g(3) * t104 - t15 * t155 + t8 * t156 + t5) * m(7) + t36;];
tau  = t9;
