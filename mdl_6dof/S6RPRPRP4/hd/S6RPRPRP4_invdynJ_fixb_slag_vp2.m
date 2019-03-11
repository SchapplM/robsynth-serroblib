% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:13
% EndTime: 2019-03-09 03:11:38
% DurationCPUTime: 17.67s
% Computational Cost: add. (4194->598), mult. (8176->746), div. (0->0), fcn. (4490->10), ass. (0->265)
t434 = m(6) + m(7);
t430 = Ifges(6,1) + Ifges(7,1);
t407 = Ifges(7,4) + Ifges(6,5);
t433 = Ifges(6,6) - Ifges(7,6);
t432 = -mrSges(6,3) - mrSges(7,2);
t174 = sin(qJ(5));
t178 = cos(qJ(3));
t299 = qJD(1) * t178;
t265 = t174 * t299;
t177 = cos(qJ(5));
t296 = qJD(3) * t177;
t118 = -t265 + t296;
t431 = -t118 / 0.2e1;
t408 = -Ifges(5,4) + Ifges(4,5);
t406 = Ifges(5,5) - Ifges(4,6);
t429 = Ifges(7,2) + Ifges(6,3);
t214 = pkin(5) * t174 - qJ(6) * t177;
t234 = t174 * mrSges(7,1) - t177 * mrSges(7,3);
t236 = mrSges(6,1) * t174 + mrSges(6,2) * t177;
t427 = -m(7) * t214 - t234 - t236;
t393 = t407 * t174 + t433 * t177;
t326 = Ifges(7,5) * t177;
t329 = Ifges(6,4) * t177;
t390 = t430 * t174 - t326 + t329;
t364 = pkin(3) + pkin(8);
t426 = t434 * t364 - t432;
t175 = sin(qJ(3));
t300 = qJD(1) * t175;
t142 = qJD(5) + t300;
t173 = cos(pkin(9));
t354 = pkin(1) * t173;
t149 = -pkin(2) - t354;
t238 = mrSges(4,1) * t175 + mrSges(4,2) * t178;
t319 = t178 * mrSges(5,3);
t165 = t178 * pkin(3);
t162 = t175 * qJ(4);
t252 = -pkin(2) - t162;
t206 = t252 - t165;
t78 = (t206 - t354) * qJD(1);
t425 = t149 * qJD(1) * t238 + t78 * (-mrSges(5,2) * t175 - t319) + (t175 * t393 + t178 * t429) * t142 / 0.2e1;
t287 = qJD(1) * qJD(3);
t122 = qJDD(1) * t175 + t178 * t287;
t113 = qJDD(5) + t122;
t121 = -t178 * qJDD(1) + t175 * t287;
t291 = qJD(5) * t177;
t58 = qJD(3) * t291 - qJD(5) * t265 + qJDD(3) * t174 - t177 * t121;
t31 = -mrSges(6,2) * t113 - mrSges(6,3) * t58;
t32 = -mrSges(7,2) * t58 + mrSges(7,3) * t113;
t341 = t31 + t32;
t298 = qJD(3) * t174;
t117 = t177 * t299 + t298;
t293 = qJD(5) * t117;
t57 = qJDD(3) * t177 + t121 * t174 - t293;
t29 = mrSges(6,1) * t113 - mrSges(6,3) * t57;
t30 = -t113 * mrSges(7,1) + t57 * mrSges(7,2);
t342 = t29 - t30;
t334 = mrSges(6,3) * t118;
t72 = mrSges(6,1) * t142 - t334;
t336 = mrSges(7,2) * t118;
t73 = -mrSges(7,1) * t142 + t336;
t339 = t72 - t73;
t335 = mrSges(6,3) * t117;
t70 = -mrSges(6,2) * t142 - t335;
t337 = mrSges(7,2) * t117;
t71 = mrSges(7,3) * t142 - t337;
t340 = t70 + t71;
t372 = t174 * t339 - t177 * t340;
t95 = t122 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t424 = t372 * qJD(5) - t341 * t174 - t342 * t177 - t95;
t423 = -t142 / 0.2e1;
t167 = qJ(1) + pkin(9);
t156 = sin(t167);
t422 = g(2) * t156;
t172 = sin(pkin(9));
t148 = pkin(1) * t172 + pkin(7);
t127 = t148 * qJD(1);
t295 = qJD(3) * t178;
t419 = qJD(2) * qJD(3) + t148 * qJDD(1);
t44 = qJDD(2) * t178 - t127 * t295 - t175 * t419;
t198 = qJDD(4) - t44;
t38 = -qJDD(3) * pkin(3) + t198;
t171 = qJD(3) * qJ(4);
t89 = t175 * qJD(2) + t178 * t127;
t77 = -t171 - t89;
t421 = -qJD(3) * t77 - t38;
t233 = t178 * mrSges(5,2) - t175 * mrSges(5,3);
t239 = mrSges(4,1) * t178 - mrSges(4,2) * t175;
t417 = t233 - t239;
t161 = t178 * qJD(2);
t68 = -t175 * (pkin(4) * qJD(1) + t127) + t161;
t381 = -t68 + qJD(4);
t59 = -qJD(3) * t364 + t381;
t66 = (-t178 * t364 + t252 - t354) * qJD(1);
t18 = -t174 * t66 + t177 * t59;
t19 = t174 * t59 + t177 * t66;
t21 = pkin(4) * t122 - qJDD(3) * t364 + t198;
t126 = t149 * qJDD(1);
t294 = qJD(4) * t175;
t184 = -qJ(4) * t122 - qJD(1) * t294 + t126;
t25 = t121 * t364 + t184;
t292 = qJD(5) * t174;
t3 = t174 * t21 + t177 * t25 + t59 * t291 - t292 * t66;
t4 = -qJD(5) * t19 - t174 * t25 + t177 * t21;
t240 = t174 * t3 + t177 * t4;
t416 = t18 * t292 - t19 * t291 - t240;
t14 = -pkin(5) * t142 + qJD(6) - t18;
t17 = qJ(6) * t142 + t19;
t1 = qJ(6) * t113 + qJD(6) * t142 + t3;
t2 = -pkin(5) * t113 + qJDD(6) - t4;
t241 = t1 * t174 - t177 * t2;
t415 = -t14 * t292 - t17 * t291 - t241;
t325 = Ifges(5,6) * t175;
t216 = -t178 * Ifges(5,3) - t325;
t109 = Ifges(6,4) * t117;
t328 = Ifges(7,5) * t117;
t404 = t118 * t430 + t142 * t407 - t109 + t328;
t331 = Ifges(6,4) * t118;
t48 = -Ifges(6,2) * t117 + Ifges(6,6) * t142 + t331;
t414 = Ifges(5,5) * qJD(3) + qJD(1) * t216 + t174 * t404 + t177 * t48;
t367 = t57 / 0.2e1;
t365 = t58 / 0.2e1;
t412 = t1 * mrSges(7,3);
t411 = t3 * mrSges(6,2);
t363 = t113 / 0.2e1;
t350 = g(3) * t178;
t410 = mrSges(6,1) + mrSges(7,1);
t409 = mrSges(6,2) - mrSges(7,3);
t405 = (-Ifges(6,4) + Ifges(7,5)) * t58 + t430 * t57 + t407 * t113;
t215 = pkin(5) * t177 + qJ(6) * t174;
t205 = -pkin(4) - t215;
t403 = qJD(5) * t215 - qJD(6) * t177 + qJD(4) - t161 - (qJD(1) * t205 - t127) * t175;
t94 = mrSges(5,1) * t121 - qJDD(3) * mrSges(5,3);
t402 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t121 - t94;
t280 = mrSges(5,1) * t299;
t131 = -qJD(3) * mrSges(5,3) - t280;
t62 = mrSges(6,1) * t117 + mrSges(6,2) * t118;
t401 = -t131 + t62;
t333 = Ifges(4,4) * t175;
t228 = t178 * Ifges(4,2) + t333;
t108 = Ifges(7,5) * t118;
t45 = Ifges(7,6) * t142 + Ifges(7,3) * t117 + t108;
t400 = Ifges(4,6) * qJD(3) + qJD(1) * t228 + t177 * t45;
t279 = mrSges(4,3) * t299;
t130 = -qJD(3) * mrSges(4,2) + t279;
t398 = t131 - t130;
t273 = mrSges(4,3) * t300;
t274 = mrSges(5,1) * t300;
t397 = t273 + t274 + (-mrSges(4,1) + mrSges(5,2)) * qJD(3);
t395 = t175 * t390 + t178 * t407;
t324 = Ifges(5,6) * t178;
t394 = t175 * (-Ifges(5,2) * t178 + t325) + t178 * (Ifges(5,3) * t175 - t324);
t392 = -t174 * t433 + t177 * t407;
t391 = t175 * t406 + t178 * t408;
t327 = Ifges(7,5) * t174;
t330 = Ifges(6,4) * t174;
t389 = t177 * t430 + t327 - t330;
t386 = t113 * t429 + t407 * t57 - t433 * t58;
t297 = qJD(3) * t175;
t43 = t175 * qJDD(2) - t127 * t297 + t178 * t419;
t385 = -t175 * t44 + t178 * t43;
t35 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t43;
t384 = t175 * t38 - t178 * t35;
t157 = cos(t167);
t383 = g(1) * t157 + t422;
t382 = -t57 * Ifges(6,4) / 0.2e1 - t113 * Ifges(6,6) / 0.2e1 + Ifges(7,5) * t367 + Ifges(7,6) * t363 + (Ifges(6,2) + Ifges(7,3)) * t365;
t288 = m(5) + t434;
t380 = mrSges(3,2) - mrSges(4,3) - mrSges(5,1);
t152 = Ifges(4,4) * t299;
t379 = Ifges(4,1) * t300 + Ifges(4,5) * qJD(3) - t117 * t433 + t118 * t407 + t142 * t429 + t152;
t378 = -mrSges(3,1) + t417;
t377 = -t319 + t427 * t178 + (m(5) * pkin(3) - mrSges(5,2) + t426) * t175;
t343 = pkin(4) + t148;
t110 = t343 * t175;
t302 = t165 + t162;
t106 = t149 - t302;
t163 = t178 * pkin(8);
t85 = t106 - t163;
t338 = t174 * t110 + t177 * t85;
t314 = qJ(4) * t178;
t213 = pkin(8) * t175 - t314;
t248 = pkin(3) * t297 - t294;
t75 = qJD(3) * t213 + t248;
t97 = t343 * t295;
t13 = -qJD(5) * t338 - t174 * t75 + t177 * t97;
t375 = -m(7) * pkin(5) - t410;
t374 = -m(7) * qJ(6) + t409;
t366 = -t58 / 0.2e1;
t362 = -t117 / 0.2e1;
t361 = t117 / 0.2e1;
t359 = t118 / 0.2e1;
t176 = sin(qJ(1));
t353 = pkin(1) * t176;
t179 = cos(qJ(1));
t166 = t179 * pkin(1);
t69 = pkin(4) * t299 + t89;
t153 = pkin(3) * t300;
t90 = qJD(1) * t213 + t153;
t34 = t174 * t69 + t177 * t90;
t332 = Ifges(4,4) * t178;
t64 = t171 + t69;
t26 = pkin(5) * t117 - qJ(6) * t118 + t64;
t313 = qJD(3) * t26;
t312 = qJD(3) * t64;
t309 = t157 * t178;
t308 = t174 * t175;
t307 = t174 * t178;
t305 = t175 * t177;
t304 = t177 * t178;
t137 = t178 * t148;
t111 = t178 * pkin(4) + t137;
t290 = qJD(5) * t178;
t289 = qJD(5) * t364;
t61 = mrSges(7,1) * t117 - mrSges(7,3) * t118;
t283 = -t61 - t401;
t272 = t157 * pkin(2) + t156 * pkin(7) + t166;
t268 = t174 * t290;
t253 = t157 * pkin(7) - t353;
t251 = -t287 / 0.2e1;
t88 = t127 * t175 - t161;
t237 = mrSges(6,1) * t177 - mrSges(6,2) * t174;
t235 = t177 * mrSges(7,1) + t174 * mrSges(7,3);
t226 = -Ifges(6,2) * t174 + t329;
t225 = Ifges(6,2) * t177 + t330;
t219 = Ifges(7,3) * t174 + t326;
t218 = -Ifges(7,3) * t177 + t327;
t217 = -t175 * Ifges(5,2) - t324;
t212 = t14 * t174 + t17 * t177;
t211 = t18 * t174 - t19 * t177;
t33 = -t174 * t90 + t177 * t69;
t208 = pkin(3) * t309 + t157 * t162 + t272;
t41 = t110 * t177 - t174 * t85;
t12 = t110 * t291 + t174 * t97 + t177 * t75 - t292 * t85;
t202 = t175 * (Ifges(4,1) * t178 - t333);
t190 = Ifges(6,6) * t178 + t175 * t225;
t189 = Ifges(7,6) * t178 + t175 * t218;
t22 = -pkin(4) * t121 - t35;
t186 = qJD(5) * t212 + t241;
t185 = -qJD(5) * t211 + t240;
t124 = qJ(4) + t214;
t120 = -qJ(4) * t299 + t153;
t119 = t233 * qJD(1);
t107 = t237 * t178;
t102 = Ifges(5,4) * qJD(3) + qJD(1) * t217;
t98 = -qJ(4) * t295 + t248;
t96 = t343 * t297;
t93 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t122;
t84 = -t156 * t308 + t157 * t177;
t83 = t156 * t305 + t157 * t174;
t82 = t156 * t177 + t157 * t308;
t81 = t156 * t174 - t157 * t305;
t74 = -qJD(3) * pkin(3) + qJD(4) + t88;
t63 = t178 * t215 + t111;
t60 = pkin(5) * t118 + qJ(6) * t117;
t39 = pkin(3) * t121 + t184;
t37 = -pkin(5) * t175 - t41;
t36 = qJ(6) * t175 + t338;
t28 = -pkin(5) * t299 - t33;
t27 = qJ(6) * t299 + t34;
t24 = (-qJD(5) * t214 + qJD(6) * t174) * t178 + (-t148 + t205) * t297;
t16 = mrSges(6,1) * t58 + mrSges(6,2) * t57;
t15 = mrSges(7,1) * t58 - mrSges(7,3) * t57;
t7 = -pkin(5) * t295 - t13;
t6 = qJ(6) * t295 + qJD(6) * t175 + t12;
t5 = pkin(5) * t58 - qJ(6) * t57 - qJD(6) * t118 + t22;
t8 = [t122 * t175 * Ifges(4,1) + (-m(3) * t166 - m(4) * t272 - m(5) * t208 - mrSges(2,1) * t179 + mrSges(2,2) * t176 - t434 * (t156 * pkin(4) + pkin(8) * t309 + t208) + t375 * t82 + t374 * t81 + t432 * t309 + t378 * t157 + t380 * t156) * g(2) + (m(3) * t353 + mrSges(2,1) * t176 + mrSges(2,2) * t179 - t434 * (t157 * pkin(4) + t253) + t375 * t84 + t374 * t83 + (-m(5) - m(4)) * t253 + t380 * t157 + (m(4) * pkin(2) - m(5) * t206 + t178 * t426 - t252 * t434 - t378) * t156) * g(1) + (t395 * t359 + t189 * t361 + t190 * t362 + t391 * qJD(3) / 0.2e1 + t425) * qJD(3) + t338 * t31 + m(6) * (t111 * t22 + t12 * t19 + t13 * t18 + t3 * t338 + t4 * t41 - t64 * t96) + (t178 * (-Ifges(4,2) * t175 + t332) + t202) * t287 / 0.2e1 + (t175 * t429 - t393 * t178) * t363 + (-t219 * t361 - t226 * t362 - t389 * t359 + t392 * t423) * t290 - t175 * t411 - t405 * t307 / 0.2e1 + (t175 * t407 - t178 * t390) * t367 + (t175 * t408 - t178 * t406) * qJDD(3) / 0.2e1 + t402 * t137 - t400 * t297 / 0.2e1 + t394 * t251 + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t382) * t304 + m(5) * (t106 * t39 + t78 * t98) + (t397 * t295 + t398 * t297 + m(4) * ((-t175 * t89 + t178 * t88) * qJD(3) + t385) + m(5) * ((t175 * t77 + t178 * t74) * qJD(3) + t384) + (-t93 + t95) * t175) * t148 + (-Ifges(4,4) * t121 + Ifges(4,5) * qJDD(3) + t386) * t175 / 0.2e1 + t122 * t332 / 0.2e1 + (t88 * mrSges(4,3) + t74 * mrSges(5,1) + t379 / 0.2e1 - t102 / 0.2e1 - t19 * mrSges(6,2) - t14 * mrSges(7,1) + t17 * mrSges(7,3) + t18 * mrSges(6,1)) * t295 + t106 * (-mrSges(5,2) * t121 - mrSges(5,3) * t122) + t98 * t119 + t111 * t16 - t96 * t62 + t22 * t107 + t12 * t70 + t6 * t71 + t13 * t72 + t7 * t73 + t24 * t61 + t63 * t15 + t41 * t29 + t36 * t32 + t37 * t30 + t2 * (-mrSges(7,1) * t175 - mrSges(7,2) * t307) + t4 * (mrSges(6,1) * t175 + mrSges(6,3) * t307) + (t297 * t77 + t384) * mrSges(5,1) + (-t297 * t89 + t385) * mrSges(4,3) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t173 - 0.2e1 * mrSges(3,2) * t172 + m(3) * (t172 ^ 2 + t173 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t48 * t268 / 0.2e1 + m(7) * (t1 * t36 + t14 * t7 + t17 * t6 + t2 * t37 + t24 * t26 + t5 * t63) + t414 * t297 / 0.2e1 + (mrSges(6,2) * t64 + mrSges(7,2) * t14 - mrSges(6,3) * t18 - mrSges(7,3) * t26) * (t174 * t297 - t177 * t290) + (-mrSges(6,1) * t64 - mrSges(7,1) * t26 + mrSges(7,2) * t17 + mrSges(6,3) * t19) * (t175 * t296 + t268) - (t174 * t45 + t177 * t404) * t290 / 0.2e1 + t5 * t235 * t178 + (m(4) * t149 - t239) * t126 + t178 * (Ifges(4,4) * t122 - Ifges(4,2) * t121 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t178 * (Ifges(5,5) * qJDD(3) - Ifges(5,6) * t122 + Ifges(5,3) * t121) / 0.2e1 - t175 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t122 + Ifges(5,6) * t121) / 0.2e1 + t149 * (mrSges(4,1) * t121 + mrSges(4,2) * t122) + t121 * t216 / 0.2e1 - t122 * t217 / 0.2e1 - t121 * t228 / 0.2e1 + t39 * t233 + (Ifges(7,6) * t175 - t178 * t218) * t365 + (Ifges(6,6) * t175 - t178 * t225) * t366 + t175 * t412; m(3) * qJDD(2) + (-m(3) - m(4) - t288) * g(3) + (t15 + t16 + (t174 * t340 + t177 * t339 + t397) * qJD(3) + m(4) * (qJD(3) * t88 + t43) + m(5) * (qJD(3) * t74 - t35) + m(6) * (t18 * t296 + t19 * t298 + t22) + m(7) * (-t14 * t296 + t17 * t298 + t5) + t402) * t175 + (t93 + (t130 - t283) * qJD(3) + m(4) * (qJD(3) * t89 + t44) + m(5) * t421 + m(6) * (t312 + t416) + m(7) * (t313 + t415) + t424) * t178; (-t94 + t16) * qJ(4) + (-m(5) * t302 - t434 * (t163 + t302) + t427 * t175 + t417) * g(3) + (t225 / 0.2e1 - t218 / 0.2e1) * t293 + (-t48 / 0.2e1 + t45 / 0.2e1) * t291 - t404 * t292 / 0.2e1 + t406 * t121 + t408 * t122 + t403 * t61 + (-t19 * (-mrSges(6,2) * t178 + mrSges(6,3) * t305) - t17 * (mrSges(7,2) * t305 + mrSges(7,3) * t178) - t18 * (mrSges(6,1) * t178 - mrSges(6,3) * t308) - t14 * (-mrSges(7,1) * t178 + mrSges(7,2) * t308) + (t190 / 0.2e1 - t189 / 0.2e1) * t117 + t395 * t431 + (t394 / 0.2e1 - t202 / 0.2e1) * qJD(1) - t425) * qJD(1) + t400 * t300 / 0.2e1 + t401 * qJD(4) + (-m(5) * t74 + t273 - t397) * t89 + (-m(5) * t77 - t279 - t398) * t88 + t389 * t367 + t391 * t251 + t392 * t363 - (t390 * t118 + t393 * t142) * qJD(5) / 0.2e1 + t383 * t238 - (-Ifges(4,2) * t300 + t152 + t379) * t299 / 0.2e1 + (-qJ(4) * t288 * t309 + t157 * t377) * g(1) + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + t124 * t15 - t120 * t119 - pkin(3) * t95 - t68 * t62 - t34 * t70 - t27 * t71 - t33 * t72 - t28 * t73 + t38 * mrSges(5,2) - t43 * mrSges(4,2) + t44 * mrSges(4,1) - t35 * mrSges(5,3) + t102 * t299 / 0.2e1 - t74 * t280 + (-t288 * t314 + t377) * t422 - t77 * t274 - t414 * t300 / 0.2e1 + (-t350 + t415) * mrSges(7,2) + (-t350 + t416) * mrSges(6,3) + t142 * (t26 * t235 + t64 * t237) + (t289 * t339 - t341 * t364 + t382) * t174 + (-t340 * t289 - t342 * t364 + t405 / 0.2e1) * t177 + (t124 * t5 - t14 * t28 - t17 * t27 - t186 * t364 + t26 * t403) * m(7) + (qJ(4) * t22 - t18 * t33 - t185 * t364 - t19 * t34 + t381 * t64) * m(6) + (-pkin(3) * t38 - qJ(4) * t35 - qJD(4) * t77 - t120 * t78) * m(5) + t5 * t234 + t22 * t236 + t219 * t365 + t226 * t366; t283 * qJD(3) + t288 * t350 + ((t119 - t372) * qJD(1) - t383 * t288) * t175 + (t212 * t300 + t186 - t313) * m(7) + (-t211 * t300 + t185 - t312) * m(6) + (t300 * t78 - t421) * m(5) - t424; (-t407 * t117 - t118 * t433) * t423 + (t409 * t82 + t410 * t81) * g(1) + (-t409 * t84 - t410 * t83) * g(2) + (-m(7) * t14 + t334 + t339) * t19 + (-m(7) * t17 - t335 - t340) * t18 + (-Ifges(6,2) * t118 - t109 + t404) * t361 - t411 + t412 + t386 - t26 * (mrSges(7,1) * t118 + mrSges(7,3) * t117) - t64 * (mrSges(6,1) * t118 - mrSges(6,2) * t117) + g(3) * t107 + qJD(6) * t71 - t60 * t61 - pkin(5) * t30 + qJ(6) * t32 + t4 * mrSges(6,1) - t2 * mrSges(7,1) + (-t117 * t430 + t108 - t331 + t45) * t431 + t235 * t350 + (-g(2) * (pkin(5) * t83 - qJ(6) * t84) - g(1) * (-pkin(5) * t81 + qJ(6) * t82) - t26 * t60 - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t17 + t215 * t350) * m(7) + t17 * t336 + t14 * t337 + t48 * t359 + (Ifges(7,3) * t118 - t328) * t362; t118 * t61 - t142 * t71 + (-g(1) * t81 + g(2) * t83 - g(3) * t304 + t118 * t26 - t142 * t17 + t2) * m(7) + t30;];
tau  = t8;
