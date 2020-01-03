% Calculate vector of inverse dynamics joint torques for
% S5RRRPP7
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:57
% EndTime: 2019-12-31 21:04:29
% DurationCPUTime: 19.02s
% Computational Cost: add. (3261->596), mult. (7306->744), div. (0->0), fcn. (4311->6), ass. (0->275)
t445 = Ifges(6,4) + Ifges(5,5);
t446 = Ifges(5,4) + Ifges(4,5);
t444 = Ifges(6,2) + Ifges(5,3);
t443 = Ifges(4,6) - Ifges(5,6);
t450 = Ifges(6,5) - t446;
t439 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t182 = cos(qJ(3));
t449 = t445 * t182;
t179 = sin(qJ(3));
t448 = t445 * t179;
t180 = sin(qJ(2));
t306 = qJD(1) * t180;
t281 = t179 * t306;
t293 = t182 * qJD(2);
t129 = t281 - t293;
t369 = t129 / 0.2e1;
t280 = t182 * t306;
t130 = qJD(2) * t179 + t280;
t368 = -t130 / 0.2e1;
t367 = t130 / 0.2e1;
t183 = cos(qJ(2));
t305 = qJD(1) * t183;
t440 = -t305 + qJD(3);
t366 = t440 / 0.2e1;
t365 = -t440 / 0.2e1;
t447 = -mrSges(5,2) - mrSges(4,3);
t442 = Ifges(6,6) - Ifges(5,6);
t441 = -Ifges(4,3) - Ifges(5,2);
t405 = -t443 * t179 + t446 * t182;
t403 = t444 * t179 + t449;
t437 = m(6) * qJ(5) + mrSges(6,3) + t447;
t342 = Ifges(4,4) * t179;
t386 = t439 * t182 - t342 + t448;
t236 = t182 * mrSges(5,1) + t179 * mrSges(5,3);
t238 = mrSges(4,1) * t182 - mrSges(4,2) * t179;
t332 = t179 * mrSges(6,2);
t436 = -t236 - t238 - t332;
t292 = qJD(1) * qJD(2);
t135 = t183 * qJDD(1) - t180 * t292;
t426 = t135 / 0.2e1;
t136 = qJDD(1) * t180 + t183 * t292;
t435 = t136 / 0.2e1;
t434 = t445 * t130;
t169 = pkin(6) * t305;
t433 = qJD(4) * t179 + t169;
t432 = t445 * t129;
t167 = Ifges(3,4) * t305;
t123 = Ifges(4,4) * t129;
t391 = t130 * t439 - t440 * t450 - t123 + t432;
t422 = t129 * t444 - t440 * t442 + t434;
t431 = Ifges(3,1) * t306 + Ifges(3,5) * qJD(2) + t179 * t422 + t182 * t391 + t167;
t345 = Ifges(3,4) * t180;
t413 = t183 * Ifges(3,2);
t228 = t345 + t413;
t430 = Ifges(6,5) * t367 + Ifges(3,6) * qJD(2) / 0.2e1 + Ifges(6,3) * t365 + qJD(1) * t228 / 0.2e1 + t441 * t366 + t446 * t368 + (Ifges(6,6) + t443) * t369;
t147 = t440 * qJ(4);
t173 = t180 * pkin(7);
t176 = t183 * pkin(2);
t283 = -pkin(1) - t176;
t212 = t283 - t173;
t116 = t212 * qJD(1);
t145 = qJD(2) * pkin(7) + t169;
t64 = t179 * t116 + t182 * t145;
t36 = qJ(5) * t129 + t64;
t23 = t147 + t36;
t40 = t147 + t64;
t429 = -mrSges(5,2) * t40 - mrSges(4,3) * t64 + mrSges(6,3) * t23;
t428 = -m(6) - m(5);
t127 = qJDD(3) - t135;
t302 = qJD(3) * t129;
t59 = qJDD(2) * t179 + t136 * t182 - t302;
t60 = qJD(3) * t130 - t182 * qJDD(2) + t136 * t179;
t427 = -t127 * t442 + t444 * t60 + t445 * t59;
t303 = qJD(2) * t183;
t425 = pkin(6) * t303;
t424 = mrSges(2,2) - mrSges(3,3);
t347 = mrSges(6,3) * t129;
t78 = mrSges(6,2) * t440 + t347;
t351 = mrSges(5,2) * t129;
t83 = mrSges(5,3) * t440 - t351;
t420 = t78 + t83;
t348 = mrSges(4,3) * t130;
t81 = mrSges(4,1) * t440 - t348;
t350 = mrSges(5,2) * t130;
t82 = -mrSges(5,1) * t440 + t350;
t419 = -t81 + t82;
t349 = mrSges(4,3) * t129;
t79 = -mrSges(4,2) * t440 - t349;
t418 = -t83 - t79;
t373 = pkin(3) + pkin(4);
t284 = t373 * t179;
t326 = qJ(4) * t182;
t209 = -t284 + t326;
t417 = t209 * t440 + t433;
t243 = pkin(2) * t180 - pkin(7) * t183;
t133 = t243 * qJD(1);
t110 = t179 * t133;
t161 = qJ(4) * t306;
t294 = qJD(5) * t182;
t300 = qJD(3) * t179;
t318 = t180 * t182;
t319 = t179 * t183;
t352 = pkin(7) - qJ(5);
t416 = -t300 * t352 - t294 - t110 - t161 - (-pkin(6) * t318 + qJ(5) * t319) * qJD(1);
t360 = pkin(3) * t179;
t213 = -t326 + t360;
t415 = t213 * t440 - t433;
t143 = t352 * t182;
t282 = -pkin(6) * t179 - pkin(3);
t255 = -pkin(4) + t282;
t315 = t182 * t183;
t321 = t133 * t182;
t414 = qJD(3) * t143 - qJD(5) * t179 + t321 - (-qJ(5) * t315 + t180 * t255) * qJD(1);
t240 = mrSges(3,1) * t183 - mrSges(3,2) * t180;
t412 = -mrSges(2,1) - t240;
t410 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t129 - mrSges(4,2) * t130 - mrSges(3,3) * t306;
t409 = t182 * t373;
t304 = qJD(2) * t180;
t408 = qJ(4) * t304 - qJD(4) * t183;
t407 = -t180 * t441 + t183 * t405;
t406 = -t180 * t442 + t183 * t403;
t404 = -t182 * t444 + t448;
t402 = t179 * t446 + t182 * t443;
t401 = -t127 * t441 - t443 * t60 + t446 * t59;
t125 = t135 * pkin(6);
t126 = t136 * pkin(6);
t400 = t125 * t183 + t126 * t180;
t399 = t447 * t180;
t102 = qJDD(2) * pkin(7) + t125;
t298 = qJD(3) * t182;
t324 = qJDD(1) * pkin(1);
t66 = -pkin(2) * t135 - pkin(7) * t136 - t324;
t13 = t182 * t102 + t116 * t298 - t145 * t300 + t179 * t66;
t14 = -t179 * t102 - t116 * t300 - t145 * t298 + t182 * t66;
t398 = t13 * t182 - t14 * t179;
t63 = t182 * t116 - t179 * t145;
t35 = qJ(5) * t130 + t63;
t397 = -t35 + qJD(4);
t4 = t127 * qJ(4) + qJD(4) * t440 + t13;
t204 = qJDD(4) - t14;
t6 = -pkin(3) * t127 + t204;
t396 = t179 * t6 + t182 * t4;
t181 = sin(qJ(1));
t184 = cos(qJ(1));
t395 = g(1) * t184 + g(2) * t181;
t393 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t392 = (-Ifges(4,4) + t445) * t60 + t439 * t59 - t450 * t127;
t390 = m(6) * pkin(4) + mrSges(6,1);
t327 = qJ(4) * t179;
t214 = pkin(3) * t182 + t327;
t137 = -pkin(2) - t214;
t252 = -m(6) * t373 - mrSges(6,1);
t259 = pkin(2) + t327;
t389 = t437 * t183 + (m(4) * pkin(2) - m(5) * t137 + m(6) * t259 - t182 * t252 - t436) * t180;
t388 = -t180 * t450 + t183 * t386;
t168 = pkin(6) * t306;
t144 = -qJD(2) * pkin(2) + t168;
t329 = t182 * mrSges(5,3);
t235 = t179 * mrSges(5,1) - t329;
t237 = mrSges(4,1) * t179 + mrSges(4,2) * t182;
t211 = qJ(4) * t130 - t144;
t24 = -t129 * t373 + qJD(5) + t211;
t42 = pkin(3) * t129 - t211;
t387 = t144 * t237 + t42 * t235 + t24 * (-t179 * mrSges(6,1) + mrSges(6,2) * t182);
t341 = Ifges(4,4) * t182;
t385 = t179 * t439 + t341 - t449;
t382 = pkin(7) * (-m(4) + t428);
t380 = mrSges(4,1) + mrSges(5,1) + t390;
t1 = -qJ(5) * t59 - qJD(5) * t130 - t127 * t373 + t204;
t2 = qJ(5) * t60 + qJD(5) * t129 + t4;
t379 = -t14 * mrSges(4,1) + t6 * mrSges(5,1) + t1 * mrSges(6,1) + t13 * mrSges(4,2) - t2 * mrSges(6,2) - t4 * mrSges(5,3);
t343 = Ifges(4,4) * t130;
t48 = -t129 * Ifges(4,2) + Ifges(4,6) * t440 + t343;
t377 = -t48 / 0.2e1;
t376 = t59 / 0.2e1;
t375 = -t60 / 0.2e1;
t374 = t60 / 0.2e1;
t372 = -t127 / 0.2e1;
t371 = t127 / 0.2e1;
t370 = -t129 / 0.2e1;
t364 = t179 / 0.2e1;
t359 = pkin(6) * t180;
t346 = mrSges(6,3) * t130;
t344 = Ifges(3,4) * t183;
t328 = qJ(4) * t129;
t325 = qJ(5) * t180;
t320 = t179 * t180;
t317 = t180 * t184;
t316 = t181 * t183;
t314 = t183 * t184;
t313 = t184 * t179;
t134 = t243 * qJD(2);
t310 = t176 + t173;
t138 = -pkin(1) - t310;
t312 = t179 * t134 + t138 * t298;
t155 = pkin(6) * t315;
t311 = qJD(3) * t155 + t138 * t300;
t87 = t179 * t138 + t155;
t309 = t184 * pkin(1) + t181 * pkin(6);
t299 = qJD(3) * t180;
t296 = qJD(4) * t182;
t288 = pkin(6) * t304;
t279 = t179 * t303;
t18 = -t60 * mrSges(6,1) + t59 * mrSges(6,2);
t260 = t298 / 0.2e1;
t27 = -t127 * mrSges(5,1) + t59 * mrSges(5,2);
t25 = -t127 * mrSges(6,1) - t59 * mrSges(6,3);
t258 = t292 / 0.2e1;
t154 = pkin(6) * t319;
t86 = t138 * t182 - t154;
t254 = pkin(3) * t315 + qJ(4) * t319 + t310;
t253 = pkin(2) * t314 + pkin(7) * t317 + t309;
t244 = t282 * t180;
t242 = qJDD(2) * pkin(2) - t126;
t74 = -qJ(4) * t183 + t87;
t241 = -t134 * t182 + t311;
t239 = mrSges(3,1) * t180 + mrSges(3,2) * t183;
t227 = -Ifges(4,2) * t179 + t341;
t226 = Ifges(4,2) * t182 + t342;
t221 = Ifges(3,5) * t183 - Ifges(3,6) * t180;
t216 = Ifges(6,5) * t182 + Ifges(6,6) * t179;
t215 = Ifges(6,5) * t179 - Ifges(6,6) * t182;
t77 = -pkin(6) * t280 + t110;
t210 = pkin(1) * t239;
t205 = t180 * (Ifges(3,1) * t183 - t345);
t203 = Ifges(6,5) * t59 + Ifges(6,6) * t60 - Ifges(6,3) * t127;
t106 = t179 * t316 + t182 * t184;
t108 = -t181 * t182 + t183 * t313;
t199 = -g(1) * t108 - g(2) * t106 - g(3) * t320;
t194 = Ifges(4,6) * t180 + t183 * t227;
t190 = -Ifges(6,3) * t180 + t183 * t216;
t189 = qJ(4) * t59 + qJD(4) * t130 + t242;
t33 = (-t180 * t293 - t183 * t300) * pkin(6) + t312;
t177 = t184 * pkin(6);
t175 = t183 * pkin(3);
t152 = mrSges(6,2) * t318;
t149 = qJ(4) * t318;
t141 = t352 * t179;
t140 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t305;
t120 = t259 + t409;
t117 = t237 * t180;
t109 = t179 * t181 + t182 * t314;
t107 = t181 * t315 - t313;
t92 = -t149 + (pkin(6) + t360) * t180;
t80 = -mrSges(6,1) * t440 - t346;
t76 = pkin(6) * t281 + t321;
t75 = t175 - t86;
t73 = t149 + (-pkin(6) - t284) * t180;
t71 = -mrSges(6,1) * t129 + mrSges(6,2) * t130;
t70 = mrSges(5,1) * t129 - mrSges(5,3) * t130;
t69 = pkin(3) * t130 + t328;
t68 = qJD(1) * t244 - t321;
t67 = t161 + t77;
t62 = qJ(5) * t320 + t74;
t58 = pkin(4) * t183 + t154 + t175 + (-t138 - t325) * t182;
t39 = -pkin(3) * t440 + qJD(4) - t63;
t37 = -t130 * t373 - t328;
t34 = t179 * t288 - t241;
t32 = (qJD(3) * t214 - t296) * t180 + (pkin(6) + t213) * t303;
t31 = qJD(2) * t244 + t241;
t30 = -mrSges(5,2) * t60 + mrSges(5,3) * t127;
t29 = -mrSges(4,2) * t127 - mrSges(4,3) * t60;
t28 = mrSges(6,2) * t127 + mrSges(6,3) * t60;
t26 = mrSges(4,1) * t127 - mrSges(4,3) * t59;
t22 = t33 + t408;
t21 = (t296 + (-t327 - t409) * qJD(3)) * t180 + (-pkin(6) + t209) * t303;
t20 = -t373 * t440 + t397;
t19 = mrSges(4,1) * t60 + mrSges(4,2) * t59;
t17 = mrSges(5,1) * t60 - mrSges(5,3) * t59;
t16 = (-pkin(6) * qJD(2) + qJ(5) * qJD(3)) * t318 + (qJD(5) * t180 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t183) * t179 + t312 + t408;
t15 = (-qJ(5) * t303 - t134) * t182 + (qJ(5) * t300 + qJD(2) * t255 - t294) * t180 + t311;
t9 = t59 * Ifges(4,4) - t60 * Ifges(4,2) + t127 * Ifges(4,6);
t5 = pkin(3) * t60 - t189;
t3 = -t373 * t60 + qJDD(5) + t189;
t7 = [(-t9 / 0.2e1 - t3 * mrSges(6,1) - t13 * mrSges(4,3) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + t427 / 0.2e1) * t320 + (mrSges(5,2) * t6 - mrSges(4,3) * t14 - mrSges(6,3) * t1 + t392 / 0.2e1) * t318 + t431 * t303 / 0.2e1 + (t63 * mrSges(4,1) - t39 * mrSges(5,1) - t20 * mrSges(6,1) - t64 * mrSges(4,2) + t23 * mrSges(6,2) + t40 * mrSges(5,3) - t430) * t304 + t422 * t180 * t260 + (qJD(2) * t388 - t299 * t385) * t367 + (Ifges(3,4) * t435 + Ifges(3,2) * t426 + t203 / 0.2e1 + Ifges(3,6) * qJDD(2) - t401 / 0.2e1 + pkin(6) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t135) - Ifges(4,6) * t375 + Ifges(6,3) * t372 + t344 * t258 + t442 * t374 + t441 * t371 + t450 * t376 + t379) * t183 + (qJD(2) * t406 - t299 * t404) * t369 + (qJD(2) * t407 - t299 * t402) * t366 + (mrSges(4,1) * t144 + mrSges(5,1) * t42 - mrSges(6,1) * t24 + t429) * (t180 * t298 + t279) - t140 * t288 - pkin(1) * (-mrSges(3,1) * t135 + mrSges(3,2) * t136) + t344 * t435 + t3 * t152 + (-qJDD(2) * mrSges(3,1) + t19) * t359 + (t136 * t359 + t400) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t400) + t16 * t78 + t33 * t79 + t15 * t80 + t34 * t81 + t31 * t82 + t22 * t83 + t86 * t26 + t87 * t29 + t92 * t17 + t32 * t70 + t21 * t71 + t73 * t18 + t74 * t30 + t75 * t27 + t62 * t28 + t58 * t25 - t210 * t292 + (mrSges(4,2) * t144 + t39 * mrSges(5,2) + mrSges(6,2) * t24 - t63 * mrSges(4,3) - mrSges(5,3) * t42 - t20 * mrSges(6,3)) * (-t179 * t299 + t183 * t293) - (t179 * t391 + t182 * t48) * t299 / 0.2e1 - t242 * t117 + (-m(4) * pkin(6) * t242 + Ifges(3,1) * t136 + Ifges(3,4) * t426 + Ifges(3,5) * qJDD(2) + t216 * t372 + t227 * t375 + t5 * t235 - t258 * t413 + t405 * t371 + t403 * t374 + t386 * t376) * t180 + t240 * t324 + (-m(3) * t309 - m(4) * t253 + t428 * (t109 * pkin(3) + qJ(4) * t108 + t253) + t412 * t184 + t424 * t181 + t437 * t317 - t380 * t109 + t393 * t108) * g(2) + (t428 * (-t107 * pkin(3) - qJ(4) * t106 + t177) + t424 * t184 + (-m(4) - m(3)) * t177 + t380 * t107 - t393 * t106 + (-m(6) * t283 - (-m(6) * t352 + mrSges(6,3)) * t180 + m(3) * pkin(1) + (-m(4) - m(5)) * t212 - t399 - t412) * t181) * g(1) + t228 * t426 - t410 * t425 + m(4) * (t13 * t87 + t14 * t86 + t144 * t425 + t33 * t64 + t34 * t63) + t205 * t258 + (qJD(2) * t190 - t215 * t299) * t365 + (qJD(2) * t194 - t226 * t299) * t370 + t279 * t377 + Ifges(2,3) * qJDD(1) + m(5) * (t22 * t40 + t31 * t39 + t32 * t42 + t4 * t74 + t5 * t92 + t6 * t75) + m(6) * (t1 * t58 + t15 * t20 + t16 * t23 + t2 * t62 + t21 * t24 + t3 * t73) + qJD(2) ^ 2 * t221 / 0.2e1; (t298 * t39 + t396) * mrSges(5,2) + (-t298 * t63 + t398) * mrSges(4,3) - (-Ifges(3,2) * t306 + t167 + t431) * t305 / 0.2e1 + t430 * t306 + t385 * t376 + (t364 * t48 - t387) * t305 + (t181 * t389 + t316 * t382) * g(2) + (t184 * t389 + t314 * t382) * g(1) + (t403 / 0.2e1 - t227 / 0.2e1) * t302 + (g(3) * t180 - t1 * t179 - t182 * t2 - t20 * t298) * mrSges(6,3) + t404 * t374 + (t422 / 0.2e1 + t377 + t429) * t300 + (-(t407 / 0.2e1 - t190 / 0.2e1) * t440 + (-t406 / 0.2e1 + t194 / 0.2e1) * t129 + t388 * t368 - t20 * (-t180 * mrSges(6,1) - mrSges(6,3) * t315) - t63 * (mrSges(4,1) * t180 - mrSges(4,3) * t315) - t39 * (-mrSges(5,1) * t180 + mrSges(5,2) * t315) - t64 * (-mrSges(4,2) * t180 - mrSges(4,3) * t319) - t40 * (-mrSges(5,2) * t319 + mrSges(5,3) * t180) - t23 * (mrSges(6,2) * t180 + mrSges(6,3) * t319) + (t210 - t205 / 0.2e1) * qJD(1)) * qJD(1) + t410 * t169 + t182 * t9 / 0.2e1 + t137 * t17 + t141 * t25 + t143 * t28 + Ifges(3,6) * t135 + Ifges(3,5) * t136 + (t137 * t5 - t39 * t68 - t40 * t67 + t415 * t42) * m(5) + (t418 * t300 + t419 * t298 + (t30 + t29) * t182 + (-t26 + t27) * t179 + ((-t179 * t64 - t182 * t63) * qJD(3) + t398) * m(4) + ((-t179 * t40 + t182 * t39) * qJD(3) + t396) * m(5)) * pkin(7) + t140 * t168 + t395 * t239 + t402 * t371 + t391 * t260 + t392 * t364 - t125 * mrSges(3,2) - t126 * mrSges(3,1) + t120 * t18 - t77 * t79 - t76 * t81 - t68 * t82 - t67 * t83 - pkin(2) * t19 - t221 * t292 / 0.2e1 + (pkin(2) * t242 - t144 * t169 - t63 * t76 - t64 * t77) * m(4) + t242 * t238 + t414 * t80 + t415 * t70 + t416 * t78 + t417 * t71 + (t1 * t141 + t120 * t3 + t143 * t2 + t20 * t414 + t23 * t416 + t24 * t417) * m(6) + (-t240 - m(6) * (t254 - t325) - m(4) * t310 - m(5) * t254 + (-t182 * t390 + t436) * t183 + t399) * g(3) + t3 * (t182 * mrSges(6,1) + t332) - t427 * t182 / 0.2e1 + (t386 * t367 + t387) * qJD(3) + t215 * t372 + t226 * t375 - (-t405 / 0.2e1 + t216 / 0.2e1) * qJD(3) * t440 - t5 * t236 + Ifges(3,3) * qJDD(2); -t379 + (-t129 * t446 - t130 * t443) * t365 + (t130 * t444 - t432) * t370 - t144 * (mrSges(4,1) * t130 - mrSges(4,2) * t129) - t24 * (-mrSges(6,1) * t130 - mrSges(6,2) * t129) - t42 * (mrSges(5,1) * t130 + mrSges(5,3) * t129) - t23 * t346 - t20 * t347 + (-Ifges(4,2) * t130 - t123 + t391) * t369 - t35 * t78 - t36 * t80 - t69 * t70 - t37 * t71 - pkin(3) * t27 - t373 * t25 + (t2 * qJ(4) - t1 * t373 - t20 * t36 + t23 * t397 - t24 * t37) * m(6) + (-m(5) * t40 - t349 + t418) * t63 + (-m(5) * t39 + t348 - t419) * t64 + t420 * qJD(4) + t401 + (-pkin(3) * t6 + qJ(4) * t4 + qJD(4) * t40 - t42 * t69) * m(5) + (-t129 * t439 - t343 + t422 + t434) * t368 + (t428 * (-t108 * pkin(3) + qJ(4) * t109) + t393 * t109 + t380 * t108) * g(1) + (t428 * (-t106 * pkin(3) + qJ(4) * t107) + t393 * t107 + t380 * t106) * g(2) + (-t252 * t320 - t152 - (t329 + (-m(5) * pkin(3) - mrSges(5,1)) * t179) * t180 + t117 + t428 * t149) * g(3) + t40 * t350 + t39 * t351 + (-Ifges(6,5) * t129 + Ifges(6,6) * t130) * t366 + t48 * t367 - t203 + (t28 + t30) * qJ(4); -t420 * t440 + (t70 - t71) * t130 + t25 + t27 + (-t130 * t24 - t23 * t440 + t1 + t199) * m(6) + (t130 * t42 - t40 * t440 + t199 + t6) * m(5); -t129 * t78 + t130 * t80 + (-g(3) * t183 - t129 * t23 + t130 * t20 + t180 * t395 + t3) * m(6) + t18;];
tau = t7;
