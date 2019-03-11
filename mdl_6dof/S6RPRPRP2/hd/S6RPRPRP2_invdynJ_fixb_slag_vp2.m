% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:31
% EndTime: 2019-03-09 03:04:59
% DurationCPUTime: 18.27s
% Computational Cost: add. (6737->629), mult. (14246->803), div. (0->0), fcn. (9565->14), ass. (0->272)
t399 = m(6) + m(7);
t273 = t399 + m(5);
t398 = -mrSges(6,3) - mrSges(7,2);
t382 = Ifges(6,1) + Ifges(7,1);
t397 = -Ifges(6,4) + Ifges(7,5);
t381 = Ifges(7,4) + Ifges(6,5);
t380 = Ifges(6,6) - Ifges(7,6);
t379 = -Ifges(6,3) - Ifges(7,2);
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t228 = t193 * mrSges(7,1) + t190 * mrSges(7,3);
t230 = mrSges(6,1) * t193 - mrSges(6,2) * t190;
t396 = -t228 - t230;
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t298 = sin(pkin(10));
t299 = cos(pkin(10));
t142 = t191 * t298 - t194 * t299;
t133 = t142 * qJD(1);
t390 = qJD(5) + t133;
t272 = qJD(1) * qJD(3);
t145 = qJDD(1) * t194 - t191 * t272;
t146 = qJDD(1) * t191 + t194 * t272;
t102 = t145 * t298 + t146 * t299;
t143 = t299 * t191 + t298 * t194;
t134 = t143 * qJD(1);
t210 = t193 * qJD(3) - t134 * t190;
t54 = qJD(5) * t210 + qJDD(3) * t190 + t102 * t193;
t350 = t54 / 0.2e1;
t113 = qJD(3) * t190 + t134 * t193;
t55 = qJD(5) * t113 - t193 * qJDD(3) + t102 * t190;
t348 = t55 / 0.2e1;
t101 = t145 * t299 - t146 * t298;
t100 = qJDD(5) - t101;
t347 = t100 / 0.2e1;
t186 = qJ(1) + pkin(9);
t175 = sin(t186);
t331 = g(2) * t175;
t159 = -mrSges(4,1) * t194 + mrSges(4,2) * t191;
t395 = -m(4) * pkin(2) - mrSges(3,1) + t159;
t187 = sin(pkin(9));
t168 = pkin(1) * t187 + pkin(7);
t153 = t168 * qJDD(1);
t394 = qJD(2) * qJD(3) + t153;
t155 = t168 * qJD(1);
t180 = t194 * qJDD(2);
t271 = qJD(1) * qJD(4);
t277 = qJD(3) * t194;
t61 = -t155 * t277 + qJDD(3) * pkin(3) - qJ(4) * t146 + t180 + (-t271 - t394) * t191;
t278 = qJD(3) * t191;
t90 = t191 * qJDD(2) - t155 * t278 + t194 * t394;
t65 = qJ(4) * t145 + t194 * t271 + t90;
t23 = t298 * t61 + t299 * t65;
t21 = qJDD(3) * pkin(8) + t23;
t274 = qJD(5) * t193;
t275 = qJD(5) * t190;
t188 = cos(pkin(9));
t170 = -pkin(1) * t188 - pkin(2);
t154 = t170 * qJDD(1);
t111 = -pkin(3) * t145 + qJDD(4) + t154;
t39 = -pkin(4) * t101 - pkin(8) * t102 + t111;
t182 = t194 * qJD(2);
t237 = qJ(4) * qJD(1) + t155;
t114 = -t191 * t237 + t182;
t107 = qJD(3) * pkin(3) + t114;
t279 = qJD(2) * t191;
t115 = t194 * t237 + t279;
t242 = t299 * t115;
t64 = t298 * t107 + t242;
t60 = qJD(3) * pkin(8) + t64;
t183 = t194 * pkin(3);
t152 = t170 - t183;
t132 = qJD(1) * t152 + qJD(4);
t75 = pkin(4) * t133 - pkin(8) * t134 + t132;
t3 = t190 * t39 + t193 * t21 + t75 * t274 - t275 * t60;
t27 = t190 * t75 + t193 * t60;
t4 = -qJD(5) * t27 - t190 * t21 + t193 * t39;
t234 = -t190 * t4 + t193 * t3;
t26 = -t190 * t60 + t193 * t75;
t393 = -t26 * t274 - t27 * t275 + t234;
t15 = -pkin(5) * t390 + qJD(6) - t26;
t16 = qJ(6) * t390 + t27;
t1 = qJ(6) * t100 + qJD(6) * t390 + t3;
t2 = -pkin(5) * t100 + qJDD(6) - t4;
t235 = t1 * t193 + t190 * t2;
t392 = t15 * t274 - t16 * t275 + t235;
t33 = -mrSges(6,2) * t100 - mrSges(6,3) * t55;
t34 = -mrSges(7,2) * t55 + mrSges(7,3) * t100;
t323 = t33 + t34;
t31 = mrSges(6,1) * t100 - mrSges(6,3) * t54;
t32 = -t100 * mrSges(7,1) + t54 * mrSges(7,2);
t324 = -t31 + t32;
t389 = t324 * t190 + t323 * t193;
t388 = m(3) + m(4);
t135 = t143 * qJD(3);
t386 = -t135 / 0.2e1;
t136 = t142 * qJD(3);
t385 = -t136 / 0.2e1;
t384 = mrSges(6,1) + mrSges(7,1);
t383 = mrSges(6,2) - mrSges(7,3);
t378 = t100 * t381 + t382 * t54 + t397 * t55;
t18 = mrSges(6,1) * t55 + mrSges(6,2) * t54;
t95 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t102;
t377 = t18 - t95;
t376 = t113 * t381 + t210 * t380 - t379 * t390;
t109 = Ifges(6,4) * t210;
t309 = t210 * Ifges(7,5);
t375 = t113 * t382 + t381 * t390 + t109 - t309;
t319 = mrSges(7,2) * t210;
t76 = mrSges(7,3) * t390 + t319;
t317 = mrSges(6,3) * t210;
t77 = -mrSges(6,2) * t390 + t317;
t322 = t76 + t77;
t316 = mrSges(6,3) * t113;
t78 = mrSges(6,1) * t390 - t316;
t318 = mrSges(7,2) * t113;
t79 = -mrSges(7,1) * t390 + t318;
t321 = t78 - t79;
t300 = qJDD(3) / 0.2e1;
t216 = pkin(5) * t190 - qJ(6) * t193;
t66 = t114 * t298 + t242;
t374 = -qJD(6) * t190 + t216 * t390 - t66;
t306 = t134 * mrSges(5,3);
t373 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t210 + mrSges(6,2) * t113 + t306;
t372 = Ifges(5,5) * qJD(3);
t371 = Ifges(5,6) * qJD(3);
t185 = qJ(3) + pkin(10);
t174 = sin(t185);
t176 = cos(t185);
t370 = -t176 * pkin(4) - t174 * pkin(8);
t227 = t190 * mrSges(7,1) - t193 * mrSges(7,3);
t229 = mrSges(6,1) * t190 + mrSges(6,2) * t193;
t104 = t298 * t115;
t63 = t107 * t299 - t104;
t59 = -qJD(3) * pkin(4) - t63;
t28 = -pkin(5) * t210 - t113 * qJ(6) + t59;
t368 = t28 * t227 + t59 * t229;
t367 = -t190 * t380 + t193 * t381;
t311 = Ifges(7,5) * t190;
t313 = Ifges(6,4) * t190;
t366 = t193 * t382 + t311 - t313;
t294 = t136 * t190;
t205 = t143 * t274 - t294;
t363 = -t100 * t379 - t380 * t55 + t381 * t54;
t125 = t155 * t194 + t279;
t91 = -qJD(3) * t125 - t153 * t191 + t180;
t362 = -t191 * t91 + t194 * t90;
t359 = 0.2e1 * t300;
t358 = Ifges(7,5) * t350 + Ifges(7,6) * t347 - t54 * Ifges(6,4) / 0.2e1 - t100 * Ifges(6,6) / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t348;
t231 = t176 * mrSges(5,1) - t174 * mrSges(5,2);
t357 = -t174 * t398 + t231;
t87 = pkin(4) * t142 - pkin(8) * t143 + t152;
t282 = qJ(4) + t168;
t140 = t282 * t191;
t141 = t282 * t194;
t93 = -t140 * t298 + t141 * t299;
t320 = t190 * t87 + t193 * t93;
t238 = qJD(3) * t282;
t116 = qJD(4) * t194 - t191 * t238;
t117 = -qJD(4) * t191 - t194 * t238;
t74 = t116 * t299 + t117 * t298;
t265 = pkin(3) * t278;
t85 = pkin(4) * t135 + pkin(8) * t136 + t265;
t9 = -qJD(5) * t320 - t190 * t74 + t193 * t85;
t356 = m(7) * pkin(5) + t384;
t355 = m(7) * qJ(6) - t383;
t177 = cos(t186);
t285 = t176 * t177;
t354 = (-g(1) * t285 - t176 * t331) * pkin(8);
t353 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t352 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t349 = -t55 / 0.2e1;
t346 = t210 / 0.2e1;
t345 = -t210 / 0.2e1;
t344 = -t113 / 0.2e1;
t343 = t113 / 0.2e1;
t342 = -t390 / 0.2e1;
t340 = t133 / 0.2e1;
t339 = t134 / 0.2e1;
t338 = -t134 / 0.2e1;
t335 = t190 / 0.2e1;
t192 = sin(qJ(1));
t334 = pkin(1) * t192;
t330 = g(3) * t174;
t195 = cos(qJ(1));
t184 = t195 * pkin(1);
t67 = t114 * t299 - t104;
t281 = qJD(1) * t191;
t266 = pkin(3) * t281;
t84 = pkin(4) * t134 + pkin(8) * t133 + t266;
t30 = t190 * t84 + t193 * t67;
t315 = Ifges(4,4) * t191;
t314 = Ifges(4,4) * t194;
t312 = Ifges(6,4) * t193;
t310 = Ifges(7,5) * t193;
t308 = t113 * Ifges(6,4);
t307 = t133 * mrSges(5,3);
t305 = t134 * Ifges(5,4);
t296 = t133 * t190;
t295 = t133 * t193;
t293 = t136 * t193;
t288 = t174 * t177;
t287 = t175 * t190;
t286 = t175 * t193;
t284 = t177 * t190;
t283 = t177 * t193;
t280 = qJD(1) * t194;
t71 = -mrSges(7,1) * t210 - mrSges(7,3) * t113;
t267 = t71 + t373;
t264 = mrSges(4,3) * t281;
t263 = mrSges(4,3) * t280;
t108 = Ifges(7,5) * t113;
t42 = Ifges(7,6) * t390 - Ifges(7,3) * t210 + t108;
t258 = t42 * t335;
t256 = t299 * pkin(3);
t255 = t298 * pkin(3);
t245 = -t275 / 0.2e1;
t244 = t274 / 0.2e1;
t243 = -t101 * mrSges(5,1) + t102 * mrSges(5,2);
t172 = t183 + pkin(2);
t189 = -qJ(4) - pkin(7);
t233 = t177 * t172 - t175 * t189 + t184;
t22 = -t298 * t65 + t299 * t61;
t232 = mrSges(4,1) * t191 + mrSges(4,2) * t194;
t224 = t194 * Ifges(4,2) + t315;
t223 = -Ifges(6,2) * t190 + t312;
t221 = Ifges(4,5) * t194 - Ifges(4,6) * t191;
t219 = Ifges(7,3) * t190 + t310;
t217 = t193 * pkin(5) + t190 * qJ(6);
t29 = -t190 * t67 + t193 * t84;
t40 = -t190 * t93 + t193 * t87;
t73 = t116 * t298 - t299 * t117;
t92 = t299 * t140 + t141 * t298;
t209 = -pkin(4) - t217;
t8 = t190 * t85 + t193 * t74 + t87 * t274 - t275 * t93;
t204 = t143 * t275 + t293;
t203 = t170 * qJD(1) * t232;
t202 = t191 * (Ifges(4,1) * t194 - t315);
t20 = -qJDD(3) * pkin(4) - t22;
t173 = Ifges(4,4) * t280;
t169 = -t256 - pkin(4);
t158 = -qJD(3) * mrSges(4,2) + t263;
t156 = qJD(3) * mrSges(4,1) - t264;
t139 = -t256 + t209;
t138 = Ifges(4,1) * t281 + Ifges(4,5) * qJD(3) + t173;
t137 = Ifges(4,6) * qJD(3) + qJD(1) * t224;
t130 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t146;
t129 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t145;
t127 = Ifges(5,4) * t133;
t124 = -t155 * t191 + t182;
t122 = -qJD(3) * mrSges(5,2) - t307;
t121 = t176 * t283 + t287;
t120 = t176 * t284 - t286;
t119 = t176 * t286 - t284;
t118 = t176 * t287 + t283;
t96 = mrSges(5,1) * t133 + mrSges(5,2) * t134;
t94 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t101;
t89 = t134 * Ifges(5,1) - t127 + t372;
t88 = -t133 * Ifges(5,2) + t305 + t371;
t70 = pkin(5) * t113 - qJ(6) * t210;
t48 = t143 * t216 + t92;
t45 = Ifges(6,2) * t210 + Ifges(6,6) * t390 + t308;
t36 = -pkin(5) * t142 - t40;
t35 = qJ(6) * t142 + t320;
t25 = -pkin(5) * t134 - t29;
t24 = qJ(6) * t134 + t30;
t17 = mrSges(7,1) * t55 - mrSges(7,3) * t54;
t14 = -t216 * t136 + (qJD(5) * t217 - qJD(6) * t193) * t143 + t73;
t7 = -pkin(5) * t135 - t9;
t6 = qJ(6) * t135 + qJD(6) * t142 + t8;
t5 = t55 * pkin(5) - t54 * qJ(6) - t113 * qJD(6) + t20;
t10 = [(m(4) * t170 + t159) * t154 + (-m(5) * t22 + m(6) * t20 + t377) * t92 - t375 * t293 / 0.2e1 + t376 * t135 / 0.2e1 + (t381 * t135 - t382 * t204 + t205 * t397) * t343 + (t203 + Ifges(5,5) * t385 + Ifges(5,6) * t386 + t221 * qJD(3) / 0.2e1) * qJD(3) + (t111 * mrSges(5,1) - t23 * mrSges(5,3) - Ifges(5,4) * t102 - Ifges(5,2) * t101 - Ifges(5,6) * t359 + Ifges(6,6) * t349 + Ifges(7,6) * t348 - t347 * t379 + t350 * t381 + t352 + t363 / 0.2e1) * t142 + (-t135 * t379 - t204 * t381 - t205 * t380) * t390 / 0.2e1 + (t375 * t245 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) + t358) * t190 + t111 * mrSges(5,2) - t22 * mrSges(5,3) + Ifges(5,1) * t102 + Ifges(5,4) * t101 + Ifges(5,5) * t359 + t20 * t229 + t219 * t348 + t223 * t349 + t5 * t227 + t42 * t244 + t347 * t367 + t350 * t366 + (t378 / 0.2e1 + mrSges(7,2) * t2 - mrSges(6,3) * t4) * t193) * t143 - t133 * (-Ifges(5,4) * t136 - Ifges(5,2) * t135) / 0.2e1 + t132 * (mrSges(5,1) * t135 - mrSges(5,2) * t136) + (-Ifges(5,1) * t136 - Ifges(5,4) * t135) * t339 + (-t135 * t64 + t136 * t63) * mrSges(5,3) + (t202 + t194 * (-Ifges(4,2) * t191 + t314)) * t272 / 0.2e1 + m(6) * (t26 * t9 + t27 * t8 + t3 * t320 + t4 * t40) + t320 * t33 + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t188 - 0.2e1 * mrSges(3,2) * t187 + m(3) * (t187 ^ 2 + t188 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t146 * t314 / 0.2e1 + t191 * (Ifges(4,4) * t145 + Ifges(4,5) * qJDD(3)) / 0.2e1 + (-m(5) * t233 - mrSges(2,1) * t195 + mrSges(2,2) * t192 + t398 * t288 - t399 * (pkin(4) * t285 + pkin(8) * t288 + t233) - t388 * t184 - t356 * t121 - t355 * t120 + (-t231 + t395) * t177 + t353 * t175) * g(2) + (mrSges(2,1) * t192 + mrSges(2,2) * t195 + t388 * t334 - t273 * (-t177 * t189 - t334) + t356 * t119 + t355 * t118 + t353 * t177 + (m(5) * t172 - t399 * (-t172 + t370) + t357 - t395) * t175) * g(1) + t146 * t191 * Ifges(4,1) + m(5) * (t111 * t152 + t132 * t265 + t23 * t93 + t64 * t74) + t138 * t277 / 0.2e1 - t137 * t278 / 0.2e1 + (m(4) * ((-t124 * t194 - t125 * t191) * qJD(3) + t362) + t194 * t129 - t156 * t277 - t158 * t278 - t191 * t130) * t168 + (-t124 * t277 - t125 * t278 + t362) * mrSges(4,3) - t205 * t45 / 0.2e1 + m(7) * (t1 * t35 + t14 * t28 + t15 * t7 + t16 * t6 + t2 * t36 + t48 * t5) + (-m(5) * t63 + m(6) * t59 + t373) * t73 + t96 * t265 - t136 * t258 + t88 * t386 + t89 * t385 + t145 * t224 / 0.2e1 + t194 * (Ifges(4,4) * t146 + Ifges(4,2) * t145 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t170 * (-mrSges(4,1) * t145 + mrSges(4,2) * t146) + t15 * (-mrSges(7,1) * t135 - mrSges(7,2) * t204) + t59 * (mrSges(6,1) * t205 - mrSges(6,2) * t204) + t28 * (mrSges(7,1) * t205 + mrSges(7,3) * t204) + t27 * (-mrSges(6,2) * t135 - mrSges(6,3) * t205) + t16 * (-mrSges(7,2) * t205 + mrSges(7,3) * t135) + t26 * (mrSges(6,1) * t135 + mrSges(6,3) * t204) + t74 * t122 + t93 * t94 + t6 * t76 + t8 * t77 + t9 * t78 + t7 * t79 + t14 * t71 + t40 * t31 + t48 * t17 + t35 * t34 + t36 * t32 + (Ifges(4,5) * t191 + Ifges(4,6) * t194) * t300 + (-Ifges(7,5) * t204 + Ifges(7,6) * t135 + Ifges(7,3) * t205) * t345 + (-Ifges(6,4) * t204 - Ifges(6,2) * t205 + Ifges(6,6) * t135) * t346 + t152 * t243; m(3) * qJDD(2) + t191 * t129 + t194 * t130 + (-t156 * t191 + t158 * t194) * qJD(3) + (t17 + t377) * t142 + t267 * t135 - (-t190 * t321 + t193 * t322 + t122) * t136 + (-t273 - t388) * g(3) + m(4) * (t191 * t90 + t194 * t91 + (-t124 * t191 + t125 * t194) * qJD(3)) + m(5) * (-t135 * t63 - t136 * t64 - t142 * t22) + m(7) * (t135 * t28 + t142 * t5 - t15 * t294 - t16 * t293) + m(6) * (t135 * t59 + t142 * t20 + t26 * t294 - t27 * t293) + (t94 + (-t190 * t322 - t193 * t321) * qJD(5) + m(5) * t23 + m(7) * t392 + m(6) * t393 + t389) * t143; t378 * t335 + (Ifges(6,2) * t349 - Ifges(7,3) * t348 + t347 * t380 - t358) * t193 - (t26 * mrSges(6,1) - t15 * mrSges(7,1) - t27 * mrSges(6,2) + t16 * mrSges(7,3) - t371 / 0.2e1 + t132 * mrSges(5,1) + Ifges(5,2) * t340 - Ifges(6,6) * t345 - Ifges(7,6) * t346 - t381 * t344 + t379 * t342) * t134 + (t347 * t381 + t350 * t382) * t190 + (-t132 * t266 + t63 * t66 - t64 * t67 + (t22 * t299 + t23 * t298) * pkin(3)) * m(5) - (t45 / 0.2e1 - t42 / 0.2e1) * t296 + (t113 * t366 + t367 * t390) * qJD(5) / 0.2e1 + (-(-t223 / 0.2e1 + t219 / 0.2e1) * t210 + t258 + t368) * qJD(5) + (-t127 + t89) * t340 + (-t305 + t376) * t338 + (-t321 * t274 - t322 * t275 + m(6) * ((-t190 * t27 - t193 * t26) * qJD(5) + t234) + m(7) * ((t15 * t193 - t16 * t190) * qJD(5) + t235) + t389) * (t255 + pkin(8)) + (t15 * t295 - t16 * t296 + t392) * mrSges(7,2) + (-t26 * t295 - t27 * t296 + t393) * mrSges(6,3) + (t169 * t20 - t26 * t29 - t27 * t30 - t59 * t66 + t354) * m(6) + (t139 * t5 - t15 * t25 - t16 * t24 + t374 * t28 + t354) * m(7) - (-Ifges(4,2) * t281 + t138 + t173) * t280 / 0.2e1 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (t263 - t158) * t124 + t64 * t306 + (-m(5) * t183 + t159 - t399 * (t183 - t370) + (-m(7) * t217 + t396) * t176 - t357) * g(3) + (-t310 + t312) * t350 + t311 * t348 + t313 * t349 + (t264 + t156) * t125 + t94 * t255 + t95 * t256 + t137 * t281 / 0.2e1 - t221 * t272 / 0.2e1 - t96 * t266 - (-t372 / 0.2e1 - t132 * mrSges(5,2) + Ifges(5,1) * t338 + t223 * t345 + t219 * t346 + t366 * t344 + t367 * t342 - t368) * t133 - t373 * t66 + t374 * t71 + t45 * t245 - t5 * t228 - t20 * t230 + t169 * t18 + Ifges(4,5) * t146 + Ifges(4,6) * t145 + t139 * t17 - t67 * t122 + Ifges(5,6) * t101 + Ifges(5,5) * t102 - t24 * t76 - t30 * t77 - t29 * t78 - t25 * t79 - t90 * mrSges(4,2) + t91 * mrSges(4,1) + t22 * mrSges(5,1) - t23 * mrSges(5,2) - t63 * t307 + (t232 + t273 * pkin(3) * t191 + (t398 + mrSges(5,2)) * t176 + (m(6) * pkin(4) - m(7) * t209 + mrSges(5,1) - t396) * t174) * (g(1) * t177 + t331) + t88 * t339 + (-t203 - t202 * qJD(1) / 0.2e1) * qJD(1) + (t295 / 0.2e1 + t244) * t375; t133 * t122 - t267 * t134 + (t322 * t390 - t324) * t193 + (-t321 * t390 + t323) * t190 + t243 + (-g(1) * t175 + g(2) * t177) * t273 + (t1 * t190 - t134 * t28 - t193 * t2 + t390 * (t15 * t190 + t16 * t193)) * m(7) + (-t134 * t59 + t190 * t3 + t193 * t4 + t390 * (-t190 * t26 + t193 * t27)) * m(6) + (t133 * t64 + t134 * t63 + t111) * m(5); t363 - t59 * (mrSges(6,1) * t113 + mrSges(6,2) * t210) + (-t113 * t380 + t210 * t381) * t342 + (t210 * t382 + t108 - t308 + t42) * t344 - t28 * (mrSges(7,1) * t113 - mrSges(7,3) * t210) + t352 + (-m(7) * t15 + t316 + t321) * t27 + (t227 + t229) * t330 + (-Ifges(6,2) * t113 + t109 + t375) * t345 + (t118 * t384 + t119 * t383) * g(2) + (t120 * t384 + t121 * t383) * g(1) + (t216 * t330 - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t16 - t28 * t70 - g(1) * (-pkin(5) * t120 + qJ(6) * t121) - g(2) * (-pkin(5) * t118 + qJ(6) * t119)) * m(7) + (Ifges(7,3) * t113 + t309) * t346 + (-m(7) * t16 + t317 - t322) * t26 + qJD(6) * t76 - t70 * t71 - pkin(5) * t32 + qJ(6) * t34 + t16 * t318 - t15 * t319 + t45 * t343; t113 * t71 - t390 * t76 + (-g(1) * t120 - g(2) * t118 + t28 * t113 - t16 * t390 - t190 * t330 + t2) * m(7) + t32;];
tau  = t10;
