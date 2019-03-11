% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:10:08
% EndTime: 2019-03-09 08:10:35
% DurationCPUTime: 20.66s
% Computational Cost: add. (8526->718), mult. (19084->906), div. (0->0), fcn. (13726->14), ass. (0->322)
t425 = -m(7) - m(5);
t435 = m(6) - t425;
t240 = sin(pkin(10));
t242 = cos(pkin(10));
t245 = sin(qJ(6));
t248 = cos(qJ(6));
t408 = -t240 * t245 + t242 * t248;
t179 = t408 * qJD(6);
t241 = sin(pkin(9));
t246 = sin(qJ(2));
t249 = cos(qJ(2));
t331 = cos(pkin(9));
t192 = t241 * t249 + t246 * t331;
t177 = t192 * qJD(1);
t415 = t408 * t177;
t413 = t179 + t415;
t434 = -mrSges(5,1) - mrSges(4,3);
t433 = -mrSges(5,2) + mrSges(4,1);
t423 = Ifges(4,1) + Ifges(6,3);
t268 = t248 * t240 + t245 * t242;
t180 = t268 * qJD(6);
t259 = t268 * t177;
t410 = -t259 - t180;
t234 = pkin(10) + qJ(6);
t229 = sin(t234);
t231 = cos(t234);
t277 = mrSges(6,1) * t240 + mrSges(6,2) * t242;
t362 = pkin(5) * t240;
t432 = -m(7) * t362 - mrSges(7,1) * t229 - mrSges(7,2) * t231 - t277;
t244 = -pkin(8) - qJ(5);
t350 = pkin(3) + qJ(5);
t431 = -m(7) * (-pkin(3) + t244) + mrSges(7,3) + m(6) * t350 + mrSges(6,3);
t247 = sin(qJ(1));
t430 = g(2) * t247;
t421 = -Ifges(4,6) + Ifges(5,5);
t289 = t331 * t249;
t310 = qJD(1) * t246;
t175 = -qJD(1) * t289 + t241 * t310;
t154 = -t242 * qJD(2) - t240 * t175;
t155 = -qJD(2) * t240 + t175 * t242;
t286 = t154 * t245 + t248 * t155;
t87 = t154 * t248 - t155 * t245;
t45 = -mrSges(7,1) * t286 - mrSges(7,2) * t87;
t95 = -mrSges(6,1) * t155 - mrSges(6,2) * t154;
t429 = -t45 - t95;
t204 = -mrSges(3,1) * t249 + mrSges(3,2) * t246;
t235 = qJ(2) + pkin(9);
t230 = sin(t235);
t232 = cos(t235);
t428 = t204 - t433 * t232 - (-mrSges(4,2) + mrSges(5,3)) * t230;
t221 = pkin(5) * t242 + pkin(4);
t243 = -qJ(3) - pkin(7);
t278 = -t242 * mrSges(6,1) + t240 * mrSges(6,2);
t427 = t278 - m(3) * pkin(7) - m(6) * (pkin(4) - t243) - m(7) * t221 + mrSges(2,2) - mrSges(3,3) + t434;
t304 = qJD(1) * qJD(2);
t294 = t246 * t304;
t303 = qJDD(1) * t249;
t197 = -t294 + t303;
t198 = qJDD(1) * t246 + t249 * t304;
t149 = -t197 * t331 + t198 * t241;
t117 = -qJDD(2) * t240 + t149 * t242;
t150 = t241 * t197 + t198 * t331;
t330 = qJDD(1) * pkin(1);
t167 = -pkin(2) * t197 + qJDD(3) - t330;
t253 = -qJ(4) * t150 - qJD(4) * t177 + t167;
t29 = qJD(5) * t175 + t149 * t350 + t253;
t189 = t198 * pkin(7);
t307 = qJD(3) * t246;
t136 = qJDD(2) * pkin(2) - qJ(3) * t198 - qJD(1) * t307 - t189;
t224 = pkin(7) * t303;
t308 = qJD(2) * t246;
t299 = pkin(7) * t308;
t306 = qJD(3) * t249;
t143 = qJ(3) * t197 + t224 + (-t299 + t306) * qJD(1);
t66 = t136 * t331 - t241 * t143;
t262 = qJDD(4) - t66;
t44 = t150 * pkin(4) - qJD(2) * qJD(5) - qJDD(2) * t350 + t262;
t15 = t240 * t44 + t242 * t29;
t10 = pkin(8) * t117 + t15;
t118 = qJDD(2) * t242 + t149 * t240;
t14 = -t240 * t29 + t242 * t44;
t5 = pkin(5) * t150 - pkin(8) * t118 + t14;
t233 = t249 * pkin(2);
t223 = t233 + pkin(1);
t200 = -qJD(1) * t223 + qJD(3);
t255 = -qJ(4) * t177 + t200;
t77 = t175 * t350 + t255;
t205 = t243 * t249;
t196 = qJD(1) * t205;
t181 = t241 * t196;
t203 = t243 * t246;
t195 = qJD(1) * t203;
t187 = qJD(2) * pkin(2) + t195;
t140 = t187 * t331 + t181;
t264 = qJD(4) - t140;
t356 = t177 * pkin(4);
t81 = -qJD(2) * t350 + t264 + t356;
t40 = -t240 * t77 + t242 * t81;
t22 = pkin(5) * t177 + pkin(8) * t154 + t40;
t41 = t240 * t81 + t242 * t77;
t24 = pkin(8) * t155 + t41;
t6 = t22 * t248 - t24 * t245;
t1 = qJD(6) * t6 + t10 * t248 + t245 * t5;
t7 = t22 * t245 + t24 * t248;
t2 = -qJD(6) * t7 - t10 * t245 + t248 * t5;
t426 = -t1 * t268 - t2 * t408 - t410 * t6 - t413 * t7;
t33 = qJD(6) * t286 + t117 * t245 + t118 * t248;
t388 = t33 / 0.2e1;
t34 = qJD(6) * t87 + t117 * t248 - t118 * t245;
t387 = t34 / 0.2e1;
t378 = t117 / 0.2e1;
t377 = t118 / 0.2e1;
t146 = qJDD(6) + t150;
t376 = t146 / 0.2e1;
t375 = t150 / 0.2e1;
t424 = t197 / 0.2e1;
t351 = qJD(2) / 0.2e1;
t422 = Ifges(4,5) - Ifges(5,4);
t295 = t331 * pkin(2);
t222 = -t295 - pkin(3);
t214 = -qJ(5) + t222;
t348 = -pkin(8) + t214;
t173 = t348 * t240;
t174 = t348 * t242;
t124 = t173 * t248 + t174 * t245;
t226 = pkin(2) * t310;
t287 = qJ(4) * t175 + t226;
t83 = t177 * t350 + t287;
t290 = t331 * t196;
t147 = t195 * t241 - t290;
t363 = pkin(4) * t175;
t100 = t147 - t363;
t98 = t242 * t100;
t28 = -pkin(5) * t175 + t98 + (-pkin(8) * t177 - t83) * t240;
t323 = t177 * t242;
t51 = t240 * t100 + t242 * t83;
t38 = pkin(8) * t323 + t51;
t420 = -qJD(5) * t408 - qJD(6) * t124 + t245 * t38 - t248 * t28;
t123 = -t173 * t245 + t174 * t248;
t419 = -qJD(5) * t268 + qJD(6) * t123 - t245 * t28 - t248 * t38;
t418 = t154 * Ifges(6,5);
t417 = t155 * Ifges(6,6);
t332 = qJDD(2) / 0.2e1;
t250 = cos(qJ(1));
t404 = g(1) * t250 + t430;
t416 = t232 * t404;
t414 = t433 * qJD(2) + t434 * t177;
t315 = t232 * t250;
t318 = t230 * t250;
t412 = pkin(3) * t315 + qJ(4) * t318;
t411 = -t232 * pkin(3) - t230 * qJ(4);
t188 = -pkin(7) * t294 + t224;
t407 = t188 * t249 + t189 * t246;
t103 = -mrSges(6,2) * t177 + mrSges(6,3) * t155;
t104 = mrSges(6,1) * t177 + mrSges(6,3) * t154;
t406 = -t103 * t240 - t104 * t242;
t68 = -mrSges(6,2) * t150 + mrSges(6,3) * t117;
t69 = mrSges(6,1) * t150 - mrSges(6,3) * t118;
t405 = t240 * t68 + t242 * t69;
t270 = t14 * t242 + t15 * t240;
t403 = m(4) - t425;
t169 = Ifges(4,4) * t175;
t170 = qJD(6) + t177;
t402 = Ifges(4,5) * qJD(2) - Ifges(7,5) * t87 + Ifges(7,6) * t286 + t170 * Ifges(7,3) + t423 * t177 - t169 + t417 - t418;
t401 = 0.2e1 * t332;
t364 = pkin(2) * t246;
t400 = t435 * t364 + (-mrSges(5,3) + t432) * t232 + (m(5) * pkin(3) - mrSges(5,2) + t431) * t230;
t141 = t241 * t187 - t290;
t129 = -qJD(2) * qJ(4) - t141;
t162 = mrSges(5,1) * t175 - qJD(2) * mrSges(5,3);
t399 = -m(5) * t129 - t162;
t398 = qJ(4) * t435;
t397 = -m(3) * pkin(1) - mrSges(2,1) + t428;
t396 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t102 = pkin(3) * t175 + t255;
t345 = Ifges(6,4) * t240;
t273 = Ifges(6,2) * t242 + t345;
t344 = Ifges(6,4) * t242;
t275 = Ifges(6,1) * t240 + t344;
t93 = qJD(5) - t129 - t363;
t393 = t141 * mrSges(4,3) + t102 * mrSges(5,2) - t93 * t278 - t129 * mrSges(5,1) + t154 * t275 / 0.2e1 - t155 * t273 / 0.2e1 - t200 * mrSges(4,1);
t120 = -qJD(2) * pkin(3) + t264;
t392 = t120 * mrSges(5,1) + t200 * mrSges(4,2) + t40 * mrSges(6,1) + t6 * mrSges(7,1) - t140 * mrSges(4,3) - t102 * mrSges(5,3) - t418 / 0.2e1 + t417 / 0.2e1 - t41 * mrSges(6,2) - t7 * mrSges(7,2);
t390 = Ifges(7,4) * t388 + Ifges(7,2) * t387 + Ifges(7,6) * t376;
t389 = Ifges(7,1) * t388 + Ifges(7,4) * t387 + Ifges(7,5) * t376;
t82 = Ifges(7,4) * t286;
t37 = -Ifges(7,1) * t87 + Ifges(7,5) * t170 + t82;
t385 = t37 / 0.2e1;
t384 = Ifges(6,1) * t377 + Ifges(6,4) * t378 + Ifges(6,5) * t375;
t383 = -t286 / 0.2e1;
t382 = t286 / 0.2e1;
t381 = t87 / 0.2e1;
t380 = -t87 / 0.2e1;
t379 = -m(6) - m(7);
t374 = -t170 / 0.2e1;
t373 = t170 / 0.2e1;
t372 = -t175 / 0.2e1;
t371 = t175 / 0.2e1;
t370 = -t177 / 0.2e1;
t369 = t177 / 0.2e1;
t366 = Ifges(7,4) * t87;
t365 = pkin(2) * t241;
t361 = pkin(7) * t249;
t358 = g(3) * t232;
t352 = -qJD(2) / 0.2e1;
t176 = t192 * qJD(2);
t190 = t241 * t246 - t289;
t178 = qJD(2) * t289 - t241 * t308;
t227 = pkin(2) * t308;
t261 = -qJ(4) * t178 - qJD(4) * t192 + t227;
t56 = qJD(5) * t190 + t176 * t350 + t261;
t291 = qJD(2) * t243;
t171 = t246 * t291 + t306;
t172 = t249 * t291 - t307;
t113 = t171 * t241 - t331 * t172;
t85 = pkin(4) * t178 + t113;
t27 = t240 * t85 + t242 * t56;
t347 = Ifges(3,4) * t246;
t346 = Ifges(3,4) * t249;
t337 = t177 * Ifges(4,4);
t336 = t177 * Ifges(5,6);
t335 = t232 * mrSges(7,3);
t152 = -t331 * t203 - t205 * t241;
t115 = pkin(4) * t192 + t152;
t284 = -qJ(4) * t192 - t223;
t99 = t190 * t350 + t284;
t54 = t240 * t115 + t242 * t99;
t326 = t176 * t240;
t325 = t176 * t242;
t324 = t177 * t240;
t320 = t190 * t240;
t319 = t190 * t242;
t317 = t231 * t250;
t316 = t232 * t244;
t312 = t247 * t229;
t311 = t247 * t231;
t67 = t241 * t136 + t331 * t143;
t309 = qJD(1) * t249;
t305 = m(5) - t379;
t302 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t310) * t361;
t301 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t146;
t300 = t162 + t429;
t296 = t233 - t411;
t11 = -t34 * mrSges(7,1) + t33 * mrSges(7,2);
t292 = -qJ(4) - t362;
t133 = t150 * mrSges(5,1) + qJDD(2) * mrSges(5,2);
t59 = -t117 * mrSges(6,1) + t118 * mrSges(6,2);
t210 = t250 * t223;
t285 = -t247 * t243 + t210;
t61 = -qJDD(2) * qJ(4) - qJD(2) * qJD(4) - t67;
t280 = mrSges(3,1) * t246 + mrSges(3,2) * t249;
t274 = t249 * Ifges(3,2) + t347;
t272 = Ifges(3,5) * t249 - Ifges(3,6) * t246;
t271 = Ifges(6,5) * t240 + Ifges(6,6) * t242;
t269 = t240 * t41 + t242 * t40;
t110 = t242 * t115;
t39 = pkin(5) * t192 + t110 + (-pkin(8) * t190 - t99) * t240;
t49 = pkin(8) * t319 + t54;
t16 = -t245 * t49 + t248 * t39;
t17 = t245 * t39 + t248 * t49;
t265 = pkin(1) * t280;
t260 = t246 * (Ifges(3,1) * t249 - t347);
t114 = t171 * t331 + t241 * t172;
t148 = t195 * t331 + t181;
t153 = t241 * t203 - t205 * t331;
t121 = t408 * t190;
t48 = -pkin(4) * t149 + qJDD(5) - t61;
t225 = Ifges(3,4) * t309;
t216 = qJ(4) + t365;
t202 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t309;
t199 = -t292 + t365;
t185 = Ifges(3,1) * t310 + Ifges(3,5) * qJD(2) + t225;
t184 = Ifges(3,6) * qJD(2) + qJD(1) * t274;
t168 = Ifges(5,6) * t175;
t160 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t175;
t159 = -t230 * t312 + t317;
t158 = t229 * t250 + t230 * t311;
t157 = t229 * t318 + t311;
t156 = t230 * t317 - t312;
t145 = t150 * mrSges(5,3);
t144 = t150 * mrSges(4,2);
t139 = pkin(3) * t190 + t284;
t138 = -mrSges(5,2) * t175 - mrSges(5,3) * t177;
t137 = mrSges(4,1) * t175 + mrSges(4,2) * t177;
t132 = mrSges(5,1) * t149 - qJDD(2) * mrSges(5,3);
t131 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t150;
t130 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t149;
t127 = -t175 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t337;
t126 = Ifges(5,4) * qJD(2) - t177 * Ifges(5,2) + t168;
t125 = Ifges(5,5) * qJD(2) + t175 * Ifges(5,3) - t336;
t122 = t268 * t190;
t116 = -t190 * pkin(4) + t153;
t112 = pkin(3) * t177 + t287;
t101 = t148 - t356;
t94 = pkin(3) * t176 + t261;
t86 = -t176 * pkin(4) + t114;
t84 = -t190 * t221 + t153;
t80 = t242 * t85;
t73 = -t177 * t221 + t148;
t72 = -t154 * Ifges(6,1) + t155 * Ifges(6,4) + Ifges(6,5) * t177;
t71 = -t154 * Ifges(6,4) + t155 * Ifges(6,2) + Ifges(6,6) * t177;
t65 = mrSges(7,1) * t170 + mrSges(7,3) * t87;
t64 = -mrSges(7,2) * t170 + mrSges(7,3) * t286;
t63 = -t176 * t221 + t114;
t62 = -qJDD(2) * pkin(3) + t262;
t60 = -pkin(5) * t155 + t93;
t58 = t176 * t408 - t180 * t190;
t57 = qJD(6) * t121 + t176 * t268;
t53 = -t240 * t99 + t110;
t52 = pkin(3) * t149 + t253;
t50 = -t240 * t83 + t98;
t46 = t118 * Ifges(6,4) + t117 * Ifges(6,2) + t150 * Ifges(6,6);
t36 = Ifges(7,2) * t286 + Ifges(7,6) * t170 - t366;
t26 = -t240 * t56 + t80;
t23 = -pkin(5) * t117 + t48;
t21 = pkin(8) * t325 + t27;
t20 = -mrSges(7,2) * t146 + mrSges(7,3) * t34;
t19 = mrSges(7,1) * t146 - mrSges(7,3) * t33;
t18 = pkin(5) * t178 + t80 + (-pkin(8) * t176 - t56) * t240;
t4 = -qJD(6) * t17 + t18 * t248 - t21 * t245;
t3 = qJD(6) * t16 + t18 * t245 + t21 * t248;
t8 = [(t402 / 0.2e1 + t423 * t369 + t422 * t351 - t126 / 0.2e1 - Ifges(5,2) * t370 - Ifges(5,6) * t371 + Ifges(4,4) * t372 + Ifges(7,3) * t373 + Ifges(7,5) * t380 + Ifges(7,6) * t382 + t392) * t178 + (t167 * mrSges(4,1) + t61 * mrSges(5,1) - t52 * mrSges(5,2) - t67 * mrSges(4,3) - Ifges(5,6) * t150 + t271 * t375 + t273 * t378 + t275 * t377 + t48 * t278 + (Ifges(5,3) + Ifges(4,2)) * t149 + (-t150 / 0.2e1 - t375) * Ifges(4,4) + t421 * t401) * t190 + (-m(4) * t285 - t157 * mrSges(7,1) - t156 * mrSges(7,2) - m(6) * (qJ(5) * t315 + t210 + t412) - mrSges(6,3) * t315 - t277 * t318 + t425 * (t285 + t412) + t427 * t247 + (-m(7) * (t230 * t362 - t316) - t335 + t397) * t250) * g(2) + (-t159 * mrSges(7,1) + t158 * mrSges(7,2) + (t243 * t403 + t427) * t250 + (-m(5) * t411 + t431 * t232 + (m(6) * qJ(4) - m(7) * t292 + t277) * t230 + (m(6) + t403) * t223 - t397) * t247) * g(1) + t139 * (-t149 * mrSges(5,2) - t145) + t71 * t325 / 0.2e1 + t72 * t326 / 0.2e1 + (-m(4) * t66 + m(5) * t62 - t131 + t133) * t152 + (t1 * t121 - t122 * t2 - t57 * t6 + t58 * t7) * mrSges(7,3) + (Ifges(7,5) * t57 + Ifges(7,6) * t58) * t373 + (Ifges(7,5) * t122 + Ifges(7,6) * t121) * t376 - t184 * t308 / 0.2e1 + (-t14 * t320 + t15 * t319 + t325 * t41 - t326 * t40) * mrSges(6,3) + (Ifges(7,1) * t57 + Ifges(7,4) * t58) * t380 + (Ifges(7,1) * t122 + Ifges(7,4) * t121) * t388 + (m(4) * t67 - m(5) * t61 + t130 - t132) * t153 + t46 * t319 / 0.2e1 + Ifges(3,6) * t249 * t332 - t202 * t299 - t265 * t304 + m(6) * (t116 * t48 + t14 * t53 + t15 * t54 + t26 * t40 + t27 * t41 + t86 * t93) + m(7) * (t1 * t17 + t16 * t2 + t23 * t84 + t3 * t7 + t4 * t6 + t60 * t63) + (Ifges(7,4) * t57 + Ifges(7,2) * t58) * t382 + (Ifges(7,4) * t122 + Ifges(7,2) * t121) * t387 + t17 * t20 + t16 * t19 + (t260 + t249 * (-Ifges(3,2) * t246 + t346)) * t304 / 0.2e1 - t204 * t330 + (m(4) * t141 + t160 + t399) * t114 + t274 * t424 + Ifges(2,3) * qJDD(1) + t198 * t346 / 0.2e1 + m(4) * (-t167 * t223 + t200 * t227) + m(5) * (t102 * t94 + t139 * t52) + t249 * (Ifges(3,4) * t198 + Ifges(3,2) * t197 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t223 * (t149 * mrSges(4,1) + t144) - pkin(1) * (-mrSges(3,1) * t197 + mrSges(3,2) * t198) + t58 * t36 / 0.2e1 + t60 * (-mrSges(7,1) * t58 + mrSges(7,2) * t57) + t63 * t45 + t3 * t64 + t4 * t65 + t54 * t68 + t53 * t69 + (Ifges(5,6) * t370 + Ifges(5,3) * t371 - Ifges(4,2) * t372 + t125 / 0.2e1 - t127 / 0.2e1 + (t271 - Ifges(4,4)) * t369 + t421 * t351 - t393) * t176 + t84 * t11 + (t62 * mrSges(5,1) + t14 * mrSges(6,1) + t167 * mrSges(4,2) - t15 * mrSges(6,2) - t66 * mrSges(4,3) - t52 * mrSges(5,3) + Ifges(4,5) * t332 + Ifges(6,5) * t377 + Ifges(7,5) * t388 + Ifges(5,2) * t150 + Ifges(6,6) * t378 + Ifges(7,6) * t387 + Ifges(7,3) * t376 + t423 * t375 + (-Ifges(5,6) - Ifges(4,4) / 0.2e1) * t149 - t401 * Ifges(5,4) + t396) * t192 + (-Ifges(4,4) * t149 + Ifges(4,5) * qJDD(2) + Ifges(6,5) * t118 + Ifges(6,6) * t117 + t423 * t150 + t301) * t192 / 0.2e1 + (t272 * t351 - t302) * qJD(2) + t137 * t227 - qJDD(2) * mrSges(3,2) * t361 + t86 * t95 + (Ifges(3,1) * t198 + Ifges(3,4) * t424 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t198) + t401 * Ifges(3,5)) * t246 + t27 * t103 + t26 * t104 + t116 * t59 + t23 * (-mrSges(7,1) * t121 + mrSges(7,2) * t122) + t94 * t138 + t249 * t185 * t351 + t320 * t384 + t57 * t385 + t122 * t389 + t121 * t390 + (t197 * t361 + t407) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t407) + (-m(4) * t140 + m(5) * t120 - t414) * t113; (Ifges(7,5) * t259 + Ifges(7,6) * t415) * t374 + (Ifges(7,1) * t259 + Ifges(7,4) * t415) * t381 + (Ifges(7,4) * t259 + Ifges(7,2) * t415) * t383 - t268 * t390 + (-t232 * t398 + t400) * t430 + (t250 * t400 - t315 * t398) * g(1) + (m(6) * t93 + m(7) * t60 + t399 - t429) * qJD(4) + t426 * mrSges(7,3) + (-m(6) * (qJ(5) * t232 + t296) - t232 * mrSges(6,3) - m(5) * t296 - m(7) * (t296 - t316) - t335 - m(4) * t233 + t432 * t230 + t428) * g(3) + (t302 + (t265 - t260 / 0.2e1) * qJD(1)) * qJD(1) - t71 * t323 / 0.2e1 - t72 * t324 / 0.2e1 + (-Ifges(7,5) * t180 - Ifges(7,6) * t179) * t373 + (-Ifges(7,1) * t180 - Ifges(7,4) * t179) * t380 - t137 * t226 + (t336 + t127) * t369 + (-Ifges(4,2) * t177 - t169 + t402) * t371 + (m(4) * t364 + mrSges(4,1) * t230 + mrSges(4,2) * t232 + t280) * t404 - t272 * t304 / 0.2e1 + (-t102 * t112 - t120 * t147 + t129 * t148 - t216 * t61 + t222 * t62) * m(5) + (t184 / 0.2e1 + pkin(7) * t202) * t310 + (-qJD(5) * t269 - t101 * t93 + t214 * t270 + t216 * t48 - t40 * t50 - t41 * t51) * m(6) + (-t160 + t162) * t148 + t131 * t295 - (-Ifges(3,2) * t310 + t185 + t225) * t309 / 0.2e1 - t259 * t37 / 0.2e1 + (Ifges(4,3) + Ifges(5,1) + Ifges(3,3)) * qJDD(2) + (t168 + t126) * t372 + (t59 - t132) * t216 + t48 * t277 + (t140 * t147 - t141 * t148 - t200 * t226 + (t241 * t67 + t331 * t66) * pkin(2)) * m(4) - t240 * t46 / 0.2e1 + t222 * t133 + Ifges(3,5) * t198 + t199 * t11 + Ifges(3,6) * t197 - t188 * mrSges(3,2) - t189 * mrSges(3,1) + (-Ifges(7,4) * t180 - Ifges(7,2) * t179) * t382 - t61 * mrSges(5,3) + t62 * mrSges(5,2) + t66 * mrSges(4,1) - t67 * mrSges(4,2) - t73 * t45 + t419 * t64 + t420 * t65 + (t1 * t124 + t123 * t2 + t199 * t23 + t419 * t7 + t420 * t6 - t60 * t73) * m(7) + (Ifges(5,3) * t372 + t352 * t421 + t393) * t177 + t421 * t149 + t422 * t150 + (-Ifges(7,5) * t381 + Ifges(5,2) * t369 - Ifges(7,6) * t383 - Ifges(7,3) * t374 - t352 * t422 + t392) * t175 + (-t175 * t423 + t177 * t271 + t125 - t337) * t370 - t101 * t95 - t51 * t103 - t50 * t104 + t123 * t19 + t124 * t20 - t112 * t138 + t130 * t365 + (Ifges(6,5) * t242 - Ifges(6,6) * t240) * t375 + (Ifges(6,1) * t242 - t345) * t377 + (-Ifges(6,2) * t240 + t344) * t378 + t242 * t384 - t180 * t385 - t413 * t36 / 0.2e1 + (mrSges(7,1) * t413 + mrSges(7,2) * t410) * t60 + t414 * t147 + (-t323 * t41 + t324 * t40 - t270) * mrSges(6,3) + t405 * t214 + t406 * qJD(5) + t408 * t389 + t23 * (mrSges(7,1) * t268 + mrSges(7,2) * t408) + (Ifges(7,5) * t408 - Ifges(7,6) * t268) * t376 + (Ifges(7,4) * t408 - Ifges(7,2) * t268) * t387 + (Ifges(7,1) * t408 - Ifges(7,4) * t268) * t388; -t268 * t19 + t408 * t20 - t240 * t69 + t242 * t68 + t144 - t145 - t413 * t65 + t410 * t64 + t433 * t149 + (t406 + t414) * t177 + (t160 - t300) * t175 + (-g(1) * t247 + g(2) * t250) * (m(4) + t305) + (t1 * t408 + t175 * t60 - t2 * t268 + t410 * t7 - t413 * t6) * m(7) + (-t14 * t240 + t15 * t242 + t175 * t93 - t177 * t269) * m(6) + (-t120 * t177 - t129 * t175 + t52) * m(5) + (t140 * t177 + t141 * t175 + t167) * m(4); t408 * t19 + t268 * t20 + t410 * t65 + t413 * t64 + (t103 * t242 - t104 * t240 + t138) * t177 + t300 * qJD(2) + t133 + (-qJD(2) * t60 - t426) * m(7) + (t270 - qJD(2) * t93 - (t240 * t40 - t242 * t41) * t177) * m(6) + (qJD(2) * t129 + t102 * t177 + t62) * m(5) + (-t230 * t404 + t358) * t305 + t405; t379 * t230 * g(3) - t155 * t103 - t154 * t104 - t286 * t64 - t87 * t65 + t11 + t59 + (-t286 * t7 - t6 * t87 + t23 - t416) * m(7) + (-t154 * t40 - t155 * t41 - t416 + t48) * m(6); -t60 * (-mrSges(7,1) * t87 + mrSges(7,2) * t286) + (Ifges(7,1) * t286 + t366) * t381 + t36 * t380 + (Ifges(7,5) * t286 + Ifges(7,6) * t87) * t374 - t6 * t64 + t7 * t65 - g(1) * (mrSges(7,1) * t156 - mrSges(7,2) * t157) - g(2) * (mrSges(7,1) * t158 + mrSges(7,2) * t159) - (-mrSges(7,1) * t231 + mrSges(7,2) * t229) * t358 + (t286 * t6 - t7 * t87) * mrSges(7,3) + t301 + (Ifges(7,2) * t87 + t37 + t82) * t383 + t396;];
tau  = t8;
