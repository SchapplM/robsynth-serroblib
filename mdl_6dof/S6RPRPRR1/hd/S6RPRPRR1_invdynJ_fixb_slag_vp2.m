% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR1
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:53
% EndTime: 2019-03-09 03:34:14
% DurationCPUTime: 16.34s
% Computational Cost: add. (11731->635), mult. (25586->840), div. (0->0), fcn. (18567->18), ass. (0->299)
t237 = qJ(3) + pkin(11);
t230 = qJ(5) + t237;
t213 = sin(t230);
t214 = cos(t230);
t248 = cos(qJ(6));
t374 = mrSges(7,1) * t248;
t402 = m(7) * pkin(5);
t255 = mrSges(6,2) * t214 + (t402 + mrSges(6,1) + t374) * t213;
t431 = -m(7) - m(6);
t246 = sin(qJ(3));
t367 = Ifges(4,4) * t246;
t244 = sin(qJ(6));
t372 = mrSges(7,2) * t244;
t447 = -t213 * t372 + t214 * (-m(7) * pkin(9) - mrSges(7,3));
t236 = qJD(3) + qJD(5);
t239 = sin(pkin(11));
t241 = cos(pkin(11));
t250 = cos(qJ(3));
t184 = -t239 * t246 + t241 * t250;
t174 = t184 * qJD(1);
t322 = t250 * qJD(1);
t323 = t246 * qJD(1);
t175 = -t239 * t322 - t241 * t323;
t245 = sin(qJ(5));
t249 = cos(qJ(5));
t275 = t174 * t245 - t249 * t175;
t105 = t236 * t248 - t244 * t275;
t106 = t236 * t244 + t248 * t275;
t345 = mrSges(6,1) * t236 + mrSges(7,1) * t105 - mrSges(7,2) * t106 - mrSges(6,3) * t275;
t382 = pkin(8) * t175;
t240 = sin(pkin(10));
t215 = pkin(1) * t240 + pkin(7);
t199 = t215 * qJD(1);
t297 = qJ(4) * qJD(1) + t199;
t328 = qJD(2) * t246;
t143 = t250 * t297 + t328;
t132 = t239 * t143;
t231 = t250 * qJD(2);
t142 = -t246 * t297 + t231;
t136 = qJD(3) * pkin(3) + t142;
t92 = t241 * t136 - t132;
t78 = qJD(3) * pkin(4) + t382 + t92;
t383 = pkin(8) * t174;
t332 = t241 * t143;
t93 = t239 * t136 + t332;
t80 = t93 + t383;
t43 = -t245 * t80 + t249 * t78;
t38 = -pkin(5) * t236 - t43;
t446 = -m(6) * t43 + m(7) * t38 - t345;
t320 = qJD(1) * qJD(3);
t192 = qJDD(1) * t250 - t246 * t320;
t193 = qJDD(1) * t246 + t250 * t320;
t129 = t192 * t239 + t193 * t241;
t229 = t250 * qJDD(2);
t319 = qJD(1) * qJD(4);
t326 = qJD(3) * t250;
t197 = t215 * qJDD(1);
t443 = qJD(2) * qJD(3) + t197;
t91 = -t199 * t326 + qJDD(3) * pkin(3) - qJ(4) * t193 + t229 + (-t319 - t443) * t246;
t327 = qJD(3) * t246;
t117 = t246 * qJDD(2) - t199 * t327 + t250 * t443;
t94 = qJ(4) * t192 + t250 * t319 + t117;
t53 = -t239 * t94 + t241 * t91;
t42 = qJDD(3) * pkin(4) - pkin(8) * t129 + t53;
t44 = t245 * t78 + t249 * t80;
t128 = t192 * t241 - t193 * t239;
t54 = t239 * t91 + t241 * t94;
t45 = pkin(8) * t128 + t54;
t11 = -qJD(5) * t44 - t245 * t45 + t249 * t42;
t445 = m(6) * t11;
t414 = t214 * pkin(5) + t213 * pkin(9);
t444 = m(7) * t414;
t203 = -t250 * mrSges(4,1) + t246 * mrSges(4,2);
t224 = sin(t237);
t226 = cos(t237);
t442 = -mrSges(5,1) * t226 + mrSges(5,2) * t224 + t203;
t441 = -t214 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t213;
t292 = mrSges(4,1) * t246 + mrSges(4,2) * t250;
t384 = pkin(3) * t246;
t440 = m(5) * t384 + mrSges(5,1) * t224 + mrSges(5,2) * t226 + t292 + t431 * (-pkin(4) * t224 - t384) + t255;
t234 = qJDD(3) + qJDD(5);
t299 = t249 * t174 + t175 * t245;
t64 = qJD(5) * t299 + t128 * t245 + t129 * t249;
t36 = qJD(6) * t105 + t234 * t244 + t248 * t64;
t65 = -qJD(5) * t275 + t128 * t249 - t129 * t245;
t63 = qJDD(6) - t65;
t19 = mrSges(7,1) * t63 - mrSges(7,3) * t36;
t37 = -qJD(6) * t106 + t234 * t248 - t244 * t64;
t20 = -mrSges(7,2) * t63 + mrSges(7,3) * t37;
t280 = -t244 * t19 + t248 * t20;
t324 = qJD(6) * t248;
t325 = qJD(6) * t244;
t110 = qJD(6) - t299;
t71 = -mrSges(7,2) * t110 + mrSges(7,3) * t105;
t72 = mrSges(7,1) * t110 - mrSges(7,3) * t106;
t439 = -t72 * t324 - t71 * t325 + t280;
t438 = t372 - t374;
t242 = cos(pkin(10));
t217 = -pkin(1) * t242 - pkin(2);
t232 = t250 * pkin(3);
t196 = t217 - t232;
t172 = qJD(1) * t196 + qJD(4);
t124 = -pkin(4) * t174 + t172;
t355 = t236 * Ifges(6,6);
t364 = Ifges(6,4) * t275;
t376 = t44 * mrSges(6,3);
t39 = pkin(9) * t236 + t44;
t58 = -pkin(5) * t299 - pkin(9) * t275 + t124;
t18 = t244 * t58 + t248 * t39;
t423 = t18 * mrSges(7,2);
t17 = -t244 * t39 + t248 * t58;
t424 = t17 * mrSges(7,1);
t358 = t110 * Ifges(7,3);
t359 = t106 * Ifges(7,5);
t361 = t105 * Ifges(7,6);
t50 = t358 + t359 + t361;
t418 = t299 * Ifges(6,2);
t73 = t355 + t364 + t418;
t437 = -t124 * mrSges(6,1) + t376 + t73 / 0.2e1 - t50 / 0.2e1 + t423 - t424 + t364 / 0.2e1 + t355 / 0.2e1;
t109 = Ifges(6,4) * t299;
t289 = mrSges(7,1) * t244 + mrSges(7,2) * t248;
t267 = t38 * t289;
t104 = Ifges(7,4) * t105;
t52 = t106 * Ifges(7,1) + t110 * Ifges(7,5) + t104;
t349 = t248 * t52;
t356 = t236 * Ifges(6,5);
t377 = t43 * mrSges(6,3);
t388 = t244 / 0.2e1;
t360 = t106 * Ifges(7,4);
t51 = t105 * Ifges(7,2) + t110 * Ifges(7,6) + t360;
t419 = t275 * Ifges(6,1);
t74 = t109 + t356 + t419;
t436 = -t124 * mrSges(6,2) - t267 - t349 / 0.2e1 + t51 * t388 + t377 - t74 / 0.2e1 - t356 / 0.2e1 - t109 / 0.2e1;
t432 = -m(3) - m(4);
t185 = t239 * t250 + t241 * t246;
t176 = t185 * qJD(3);
t429 = -t176 / 0.2e1;
t177 = t184 * qJD(3);
t428 = t177 / 0.2e1;
t16 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t56 = mrSges(6,1) * t234 - mrSges(6,3) * t64;
t422 = t16 - t56;
t107 = -mrSges(6,2) * t236 + mrSges(6,3) * t299;
t279 = -t244 * t72 + t248 * t71;
t270 = -t107 - t279;
t385 = pkin(3) * t241;
t216 = pkin(4) + t385;
t386 = pkin(3) * t239;
t171 = t245 * t216 + t249 * t386;
t413 = t214 * t438 + t441;
t126 = t184 * t245 + t185 * t249;
t274 = t249 * t184 - t185 * t245;
t81 = qJD(5) * t274 - t176 * t245 + t177 * t249;
t269 = t126 * t324 + t244 * t81;
t411 = -t17 * t324 - t18 * t325;
t160 = t199 * t250 + t328;
t118 = -t160 * qJD(3) - t197 * t246 + t229;
t410 = t117 * t250 - t118 * t246;
t409 = m(5) - t432;
t198 = t217 * qJDD(1);
t141 = -pkin(3) * t192 + qJDD(4) + t198;
t95 = -pkin(4) * t128 + t141;
t21 = -pkin(5) * t65 - pkin(9) * t64 + t95;
t10 = qJD(5) * t43 + t245 * t42 + t249 * t45;
t7 = pkin(9) * t234 + t10;
t2 = qJD(6) * t17 + t21 * t244 + t248 * t7;
t3 = -qJD(6) * t18 + t21 * t248 - t244 * t7;
t407 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t77 = pkin(5) * t275 - pkin(9) * t299;
t406 = mrSges(3,1) + m(5) * (t232 + pkin(2)) + m(4) * pkin(2) - t441 - t442;
t243 = -qJ(4) - pkin(7);
t404 = -m(4) * pkin(7) + m(5) * t243 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t401 = t36 / 0.2e1;
t400 = t37 / 0.2e1;
t397 = t63 / 0.2e1;
t394 = -t105 / 0.2e1;
t393 = -t106 / 0.2e1;
t392 = t106 / 0.2e1;
t391 = -t110 / 0.2e1;
t389 = -t175 / 0.2e1;
t247 = sin(qJ(1));
t387 = pkin(1) * t247;
t379 = t2 * t248;
t378 = t244 * t3;
t251 = cos(qJ(1));
t233 = t251 * pkin(1);
t371 = mrSges(5,3) * t174;
t370 = mrSges(5,3) * t175;
t369 = mrSges(7,3) * t244;
t368 = mrSges(7,3) * t248;
t366 = Ifges(4,4) * t250;
t365 = Ifges(5,4) * t175;
t363 = Ifges(7,4) * t244;
t362 = Ifges(7,4) * t248;
t341 = t126 * t244;
t340 = t126 * t248;
t238 = qJ(1) + pkin(10);
t225 = sin(t238);
t336 = t225 * t244;
t335 = t225 * t248;
t227 = cos(t238);
t334 = t227 * t244;
t333 = t227 * t248;
t331 = qJ(4) + t215;
t97 = t241 * t142 - t132;
t298 = qJD(3) * t331;
t150 = qJD(4) * t250 - t246 * t298;
t151 = -qJD(4) * t246 - t250 * t298;
t101 = t241 * t150 + t239 * t151;
t182 = t331 * t246;
t183 = t331 * t250;
t120 = -t239 * t182 + t241 * t183;
t329 = pkin(4) * t226 + t232;
t321 = m(5) - t431;
t317 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t63;
t222 = pkin(3) * t327;
t221 = pkin(3) * t323;
t315 = mrSges(4,3) * t323;
t314 = mrSges(4,3) * t322;
t307 = t349 / 0.2e1;
t305 = -t65 * mrSges(6,1) + t64 * mrSges(6,2);
t304 = -t325 / 0.2e1;
t144 = -pkin(4) * t175 + t221;
t145 = pkin(4) * t176 + t222;
t301 = -t128 * mrSges(5,1) + t129 * mrSges(5,2);
t96 = -t142 * t239 - t332;
t100 = -t150 * t239 + t241 * t151;
t119 = -t241 * t182 - t183 * t239;
t296 = t447 * t225;
t295 = t447 * t227;
t288 = Ifges(7,1) * t248 - t363;
t287 = t250 * Ifges(4,2) + t367;
t286 = -Ifges(7,2) * t244 + t362;
t285 = Ifges(4,5) * t250 - Ifges(4,6) * t246;
t284 = Ifges(7,5) * t248 - Ifges(7,6) * t244;
t282 = t17 * t248 + t18 * t244;
t281 = -t17 * t244 + t18 * t248;
t102 = -pkin(8) * t185 + t119;
t103 = pkin(8) * t184 + t120;
t67 = t102 * t245 + t103 * t249;
t140 = -pkin(4) * t184 + t196;
t70 = -pkin(5) * t274 - pkin(9) * t126 + t140;
t29 = t244 * t70 + t248 * t67;
t28 = -t244 * t67 + t248 * t70;
t277 = t249 * t102 - t103 * t245;
t170 = t216 * t249 - t245 * t386;
t272 = t96 - t383;
t271 = -pkin(8) * t177 + t100;
t268 = t126 * t325 - t248 * t81;
t266 = t105 * t286;
t265 = t106 * t288;
t264 = t110 * t284;
t263 = t217 * qJD(1) * t292;
t262 = t246 * (Ifges(4,1) * t250 - t367);
t256 = -qJD(6) * t282 - t378;
t254 = t256 + t379;
t14 = t36 * Ifges(7,4) + t37 * Ifges(7,2) + t63 * Ifges(7,6);
t15 = t36 * Ifges(7,1) + t37 * Ifges(7,4) + t63 * Ifges(7,5);
t8 = -pkin(5) * t234 - t11;
t253 = -t10 * mrSges(6,2) + t2 * t368 + t15 * t388 + t248 * t14 / 0.2e1 + Ifges(6,3) * t234 + (Ifges(7,1) * t244 + t362) * t401 + (Ifges(7,2) * t248 + t363) * t400 + t51 * t304 + (Ifges(7,5) * t244 + Ifges(7,6) * t248) * t397 + t8 * t438 + Ifges(6,6) * t65 + Ifges(6,5) * t64 + t11 * mrSges(6,1) + (t267 + t307) * qJD(6) + (t266 + t265 + t264) * qJD(6) / 0.2e1;
t235 = -pkin(8) + t243;
t220 = Ifges(4,4) * t322;
t202 = -qJD(3) * mrSges(4,2) + t314;
t200 = qJD(3) * mrSges(4,1) - t315;
t191 = pkin(2) + t329;
t179 = Ifges(4,1) * t323 + Ifges(4,5) * qJD(3) + t220;
t178 = Ifges(4,6) * qJD(3) + qJD(1) * t287;
t169 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t193;
t168 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t192;
t167 = Ifges(5,4) * t174;
t164 = -pkin(5) - t170;
t159 = -t199 * t246 + t231;
t153 = qJD(3) * mrSges(5,1) + t370;
t152 = -qJD(3) * mrSges(5,2) + t371;
t149 = t214 * t333 + t336;
t148 = -t214 * t334 + t335;
t147 = -t214 * t335 + t334;
t146 = t214 * t336 + t333;
t123 = -mrSges(5,1) * t174 - mrSges(5,2) * t175;
t122 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t129;
t121 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t128;
t112 = -t175 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t167;
t111 = t174 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t365;
t85 = -pkin(8) * t176 + t101;
t83 = t97 + t382;
t82 = qJD(5) * t126 + t249 * t176 + t177 * t245;
t76 = -mrSges(6,1) * t299 + mrSges(6,2) * t275;
t68 = t144 + t77;
t57 = -mrSges(6,2) * t234 + mrSges(6,3) * t65;
t47 = t245 * t272 + t249 * t83;
t33 = pkin(5) * t82 - pkin(9) * t81 + t145;
t27 = t244 * t77 + t248 * t43;
t26 = -t244 * t43 + t248 * t77;
t24 = qJD(5) * t277 + t245 * t271 + t249 * t85;
t23 = t244 * t68 + t248 * t47;
t22 = -t244 * t47 + t248 * t68;
t5 = -qJD(6) * t29 - t24 * t244 + t248 * t33;
t4 = qJD(6) * t28 + t24 * t248 + t244 * t33;
t1 = [(t17 * t268 - t18 * t269 - t2 * t341 - t3 * t340) * mrSges(7,3) + t193 * t246 * Ifges(4,1) + (m(4) * t217 + t203) * t198 + t250 * (Ifges(4,4) * t193 + Ifges(4,2) * t192) / 0.2e1 + (t263 + Ifges(5,5) * t428 + Ifges(5,6) * t429 + t285 * qJD(3) / 0.2e1) * qJD(3) + (m(4) * ((-t159 * t250 - t160 * t246) * qJD(3) + t410) - t200 * t326 - t202 * t327 - t246 * t169 + t250 * t168) * t215 + (-t159 * t326 - t160 * t327 + t410) * mrSges(4,3) - t269 * t51 / 0.2e1 + t446 * (qJD(5) * t67 + t245 * t85 - t249 * t271) + (-Ifges(7,1) * t268 - Ifges(7,4) * t269 + Ifges(7,5) * t82) * t392 + m(7) * (t17 * t5 + t18 * t4 + t2 * t29 + t28 * t3) + m(6) * (t10 * t67 + t124 * t145 + t140 * t95 + t24 * t44) + t193 * t366 / 0.2e1 - t82 * t423 + (t367 + t287) * t192 / 0.2e1 + (Ifges(5,1) * t177 - Ifges(5,4) * t176) * t389 + (-t176 * t93 - t177 * t92 + t184 * t54 - t185 * t53) * mrSges(5,3) + t174 * (Ifges(5,4) * t177 - Ifges(5,2) * t176) / 0.2e1 + t172 * (mrSges(5,1) * t176 + mrSges(5,2) * t177) + (Ifges(4,5) * t246 + Ifges(5,5) * t185 + Ifges(4,6) * t250 + Ifges(5,6) * t184) * qJDD(3) + (t262 + t250 * (-Ifges(4,2) * t246 + t366)) * t320 / 0.2e1 + t299 * (Ifges(6,4) * t81 - Ifges(6,2) * t82) / 0.2e1 + t275 * (Ifges(6,1) * t81 - Ifges(6,4) * t82) / 0.2e1 + (t95 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t64 + Ifges(6,4) * t65 + Ifges(6,5) * t234 + t284 * t397 + t286 * t400 + t288 * t401 + t289 * t8 + t304 * t52) * t126 - t82 * t376 - t81 * t377 + t123 * t222 + t82 * t424 + t112 * t428 + t111 * t429 - t14 * t341 / 0.2e1 + t15 * t340 / 0.2e1 + t179 * t326 / 0.2e1 - t178 * t327 / 0.2e1 + m(5) * (t100 * t92 + t101 * t93 + t119 * t53 + t120 * t54 + t141 * t196 + t172 * t222) + t140 * t305 + t236 * (Ifges(6,5) * t81 - Ifges(6,6) * t82) / 0.2e1 + t196 * t301 + t217 * (-mrSges(4,1) * t192 + mrSges(4,2) * t193) + t141 * (-mrSges(5,1) * t184 + mrSges(5,2) * t185) + t128 * (Ifges(5,4) * t185 + Ifges(5,2) * t184) + t129 * (Ifges(5,1) * t185 + Ifges(5,4) * t184) + t101 * t152 + t100 * t153 + t145 * t76 + t124 * (mrSges(6,1) * t82 + mrSges(6,2) * t81) + t120 * t121 + t119 * t122 + t24 * t107 + t82 * t50 / 0.2e1 - t82 * t73 / 0.2e1 + t81 * t74 / 0.2e1 + t4 * t71 + t5 * t72 + t67 * t57 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t242 - 0.2e1 * mrSges(3,2) * t240 + m(3) * (t240 ^ 2 + t242 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) - (Ifges(7,3) * t397 + Ifges(7,6) * t400 + Ifges(7,5) * t401 + t317 / 0.2e1 - Ifges(6,4) * t64 - Ifges(6,2) * t65 - Ifges(6,6) * t234 + t95 * mrSges(6,1) - t10 * mrSges(6,3) + t407) * t274 + (mrSges(2,1) * t247 - t147 * mrSges(7,1) + mrSges(2,2) * t251 - t146 * mrSges(7,2) + t431 * (-t227 * t235 - t387) + t409 * t387 + t404 * t227 + (-m(7) * (-t191 - t414) + m(6) * t191 + t406) * t225) * g(1) + t28 * t19 + t29 * t20 + t81 * t307 + (-mrSges(2,1) * t251 - t149 * mrSges(7,1) + mrSges(2,2) * t247 - t148 * mrSges(7,2) + t431 * (t227 * t191 - t225 * t235 + t233) - t409 * t233 + t404 * t225 + (-t406 - t444) * t227) * g(2) - (m(7) * t8 + t422 - t445) * t277 + t38 * (mrSges(7,1) * t269 - mrSges(7,2) * t268) + t110 * (-Ifges(7,5) * t268 - Ifges(7,6) * t269 + Ifges(7,3) * t82) / 0.2e1 + t105 * (-Ifges(7,4) * t268 - Ifges(7,2) * t269 + Ifges(7,6) * t82) / 0.2e1; m(3) * qJDD(2) + t185 * t121 + t184 * t122 + t177 * t152 - t176 * t153 + t246 * t168 + t250 * t169 - t345 * t82 - t422 * t274 + (-t200 * t246 + t202 * t250) * qJD(3) - t270 * t81 + (t57 + (-t244 * t71 - t248 * t72) * qJD(6) + t280) * t126 + (-t321 + t432) * g(3) + m(6) * (t10 * t126 + t11 * t274 - t43 * t82 + t44 * t81) + m(5) * (-t176 * t92 + t177 * t93 + t184 * t53 + t185 * t54) + m(7) * (t126 * t254 - t274 * t8 + t281 * t81 + t38 * t82) + m(4) * (t117 * t246 + t118 * t250 + (-t159 * t246 + t160 * t250) * qJD(3)); (t314 - t202) * t159 + (-t263 - t262 * qJD(1) / 0.2e1) * qJD(1) + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t411 * mrSges(7,3) + t446 * (t171 * qJD(5) - t245 * t83 + t249 * t272) + t122 * t385 + t121 * t386 + t111 * t389 + (Ifges(7,3) * t391 + Ifges(7,5) * t393 + Ifges(7,6) * t394 + t418 / 0.2e1 + t437) * t275 + (t284 * t391 + t288 * t393 + t286 * t394 + t17 * t368 + t18 * t369 - t419 / 0.2e1 + t436) * t299 + ((t239 * t54 + t241 * t53) * pkin(3) - t172 * t221 - t92 * t96 - t93 * t97) * m(5) + (t10 * t171 - t124 * t144 - t44 * t47) * m(6) + (t315 + t200) * t160 - (-Ifges(4,2) * t323 + t179 + t220) * t322 / 0.2e1 - (Ifges(5,2) * t175 + t112 + t167) * t174 / 0.2e1 + t175 * (Ifges(5,1) * t174 + t365) / 0.2e1 - t3 * t369 - t93 * t370 + t92 * t371 + t178 * t323 / 0.2e1 - t285 * t320 / 0.2e1 - t123 * t221 + Ifges(4,6) * t192 + Ifges(4,5) * t193 - qJD(3) * (Ifges(5,5) * t174 + Ifges(5,6) * t175) / 0.2e1 - t172 * (-mrSges(5,1) * t175 + mrSges(5,2) * t174) + t171 * t57 + t164 * t16 - t97 * t152 - t96 * t153 - t144 * t76 + Ifges(5,6) * t128 + Ifges(5,5) * t129 - t117 * mrSges(4,2) + t118 * mrSges(4,1) - t47 * t107 - t23 * t71 - t22 * t72 + t53 * mrSges(5,1) - t54 * mrSges(5,2) + (t164 * t8 - t17 * t22 - t18 * t23) * m(7) + (m(7) * t254 + t439) * (pkin(9) + t171) + (t227 * t440 + t295) * g(1) + (t225 * t440 + t296) * g(2) + (-m(5) * t232 - m(6) * t329 - m(7) * (t329 + t414) + t413 + t442) * g(3) + t253 + (t445 + (m(6) * t44 + m(7) * t281 - t270) * qJD(5) + t56) * t170; t279 * qJD(6) + t345 * t275 + t270 * t299 - t174 * t152 - t175 * t153 + t248 * t19 + t244 * t20 + t301 + t305 + (-g(1) * t225 + g(2) * t227) * t321 + (t110 * t281 + t2 * t244 + t3 * t248 - t275 * t38) * m(7) + (t275 * t43 - t299 * t44 + t95) * m(6) + (-t174 * t93 - t175 * t92 + t141) * m(5); (m(7) * (-t378 + t379 + t411) + t439) * pkin(9) + (-t264 / 0.2e1 - t266 / 0.2e1 - t265 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t275 + t282 * mrSges(7,3) + t436) * t299 + t256 * mrSges(7,3) + (t227 * t255 + t295) * g(1) + (t225 * t255 + t296) * g(2) + t345 * t44 + (-t358 / 0.2e1 - t361 / 0.2e1 - t359 / 0.2e1 + t437) * t275 - t8 * t402 - m(7) * (t17 * t26 + t18 * t27 + t38 * t44) + (t413 - t444) * g(3) - t43 * t107 - t27 * t71 - t26 * t72 - pkin(5) * t16 + t253; -t38 * (mrSges(7,1) * t106 + mrSges(7,2) * t105) + (Ifges(7,1) * t105 - t360) * t393 + t51 * t392 + (Ifges(7,5) * t105 - Ifges(7,6) * t106) * t391 - t17 * t71 + t18 * t72 - g(1) * (mrSges(7,1) * t148 - mrSges(7,2) * t149) - g(2) * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) + g(3) * t289 * t213 + (t105 * t17 + t106 * t18) * mrSges(7,3) + t317 + (-Ifges(7,2) * t106 + t104 + t52) * t394 + t407;];
tau  = t1;
