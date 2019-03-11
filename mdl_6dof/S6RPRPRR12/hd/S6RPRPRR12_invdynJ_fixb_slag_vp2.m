% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:18
% EndTime: 2019-03-09 04:18:49
% DurationCPUTime: 21.97s
% Computational Cost: add. (6831->695), mult. (12677->913), div. (0->0), fcn. (7304->10), ass. (0->329)
t444 = mrSges(6,3) + mrSges(7,3);
t463 = -Ifges(4,5) + Ifges(5,4);
t462 = Ifges(5,5) - Ifges(4,6);
t233 = cos(qJ(3));
t216 = t233 * qJ(4);
t229 = sin(qJ(3));
t390 = pkin(8) * t229;
t271 = -t216 + t390;
t332 = qJD(1) * t229;
t335 = pkin(3) * t332 + qJD(1) * qJ(2);
t100 = qJD(1) * t271 + t335;
t228 = sin(qJ(5));
t232 = cos(qJ(5));
t325 = qJD(5) * t232;
t327 = qJD(5) * t228;
t319 = qJD(1) * qJD(3);
t170 = qJDD(1) * t229 + t233 * t319;
t169 = -t233 * qJDD(1) + t229 * t319;
t292 = -qJD(4) * t233 + qJD(2);
t318 = qJDD(1) * qJ(2);
t242 = qJ(4) * t169 + qJD(1) * t292 + t318;
t403 = pkin(3) + pkin(8);
t50 = t170 * t403 + t242;
t237 = -pkin(1) - pkin(7);
t185 = qJDD(1) * t237 + qJDD(2);
t187 = qJD(1) * t237 + qJD(2);
t330 = qJD(3) * t229;
t94 = t185 * t233 - t187 * t330;
t265 = qJDD(4) - t94;
t58 = -pkin(4) * t169 - qJDD(3) * t403 + t265;
t172 = t233 * t187;
t212 = t233 * qJD(1);
t123 = -pkin(4) * t212 + t172;
t424 = -t123 + qJD(4);
t96 = -qJD(3) * t403 + t424;
t11 = -t100 * t327 + t228 * t58 + t232 * t50 + t96 * t325;
t160 = qJD(3) * t232 + t228 * t332;
t73 = -qJD(5) * t160 - qJDD(3) * t228 + t170 * t232;
t10 = pkin(9) * t73 + t11;
t231 = cos(qJ(6));
t227 = sin(qJ(6));
t159 = -qJD(3) * t228 + t232 * t332;
t52 = t100 * t232 + t228 * t96;
t41 = pkin(9) * t159 + t52;
t367 = t227 * t41;
t194 = t212 + qJD(5);
t51 = -t100 * t228 + t232 * t96;
t40 = -pkin(9) * t160 + t51;
t37 = pkin(5) * t194 + t40;
t13 = t231 * t37 - t367;
t12 = -qJD(5) * t52 - t228 * t50 + t232 * t58;
t158 = qJDD(5) - t169;
t72 = qJD(5) * t159 + qJDD(3) * t232 + t170 * t228;
t9 = pkin(5) * t158 - pkin(9) * t72 + t12;
t2 = qJD(6) * t13 + t10 * t231 + t227 * t9;
t475 = t2 * mrSges(7,2);
t358 = t231 * t41;
t14 = t227 * t37 + t358;
t3 = -qJD(6) * t14 - t10 * t227 + t231 * t9;
t474 = t3 * mrSges(7,1);
t473 = t11 * mrSges(6,2);
t472 = t12 * mrSges(6,1);
t445 = m(5) + m(6);
t321 = m(7) + t445;
t417 = m(7) * pkin(5);
t471 = -mrSges(6,1) - t417;
t348 = t228 * t233;
t261 = -pkin(5) * t229 - pkin(9) * t348;
t380 = pkin(9) + t403;
t165 = pkin(3) * t212 + qJ(4) * t332;
t117 = pkin(8) * t212 + t165;
t171 = t229 * t187;
t122 = -pkin(4) * t332 + t171;
t63 = -t117 * t228 + t232 * t122;
t470 = -qJD(1) * t261 + t380 * t327 - t63;
t176 = t380 * t232;
t308 = t232 * t212;
t64 = t232 * t117 + t228 * t122;
t469 = pkin(9) * t308 + qJD(5) * t176 + t64;
t350 = t227 * t228;
t109 = -t212 * t350 + t231 * t308;
t315 = qJD(5) + qJD(6);
t323 = qJD(6) * t227;
t343 = t231 * t232;
t83 = -t227 * t327 - t228 * t323 + t315 * t343;
t432 = t83 + t109;
t263 = t227 * t232 + t231 * t228;
t110 = t263 * t212;
t84 = t315 * t263;
t431 = -t84 - t110;
t468 = t233 * t315 + qJD(1);
t288 = mrSges(4,1) * t229 + mrSges(4,2) * t233;
t467 = t229 * t444 + t288;
t378 = mrSges(6,3) * t159;
t98 = -mrSges(6,2) * t194 + t378;
t377 = mrSges(6,3) * t160;
t99 = mrSges(6,1) * t194 - t377;
t426 = -t228 * t99 + t232 * t98;
t53 = mrSges(6,1) * t158 - mrSges(6,3) * t72;
t54 = -mrSges(6,2) * t158 + mrSges(6,3) * t73;
t466 = -qJD(5) * t426 - t228 * t54 - t232 * t53;
t465 = -t13 * mrSges(7,1) + t14 * mrSges(7,2);
t293 = t231 * t159 - t160 * t227;
t148 = qJDD(6) + t158;
t23 = qJD(6) * t293 + t227 * t73 + t231 * t72;
t80 = t159 * t227 + t160 * t231;
t24 = -qJD(6) * t80 - t227 * t72 + t231 * t73;
t313 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t148;
t395 = Ifges(7,4) * t80;
t186 = t212 + t315;
t398 = -t186 / 0.2e1;
t406 = -t80 / 0.2e1;
t224 = qJD(3) * qJ(4);
t102 = t122 + t224;
t74 = -pkin(5) * t159 + t102;
t464 = t474 - t475 + t313 + (Ifges(7,5) * t293 - Ifges(7,6) * t80) * t398 + (t13 * t293 + t14 * t80) * mrSges(7,3) - t74 * (mrSges(7,1) * t80 + mrSges(7,2) * t293) + (Ifges(7,1) * t293 - t395) * t406;
t404 = -m(3) - m(4);
t179 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t332;
t311 = mrSges(5,1) * t332;
t181 = -qJD(3) * mrSges(5,3) + t311;
t461 = t179 - t181;
t310 = mrSges(5,1) * t212;
t460 = mrSges(4,3) * t212 + t310 + (-mrSges(4,1) + mrSges(5,2)) * qJD(3);
t269 = t228 * t51 - t232 * t52;
t270 = t11 * t228 + t12 * t232;
t241 = -qJD(5) * t269 + t270;
t459 = qJD(3) * t102 - t241;
t371 = Ifges(5,6) * t233;
t376 = Ifges(4,4) * t229;
t458 = t229 * (-Ifges(4,2) * t233 - t376) + t233 * (Ifges(5,2) * t229 + t371);
t457 = t229 * t463 + t462 * t233;
t226 = qJ(5) + qJ(6);
t213 = sin(t226);
t214 = cos(t226);
t286 = mrSges(6,1) * t228 + mrSges(6,2) * t232;
t456 = -t213 * mrSges(7,1) - t214 * mrSges(7,2) - t286;
t320 = qJD(1) * qJD(2);
t188 = t318 + t320;
t329 = qJD(3) * t233;
t95 = t229 * t185 + t187 * t329;
t82 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t95;
t85 = -qJDD(3) * pkin(3) + t265;
t455 = -t229 * t82 - t233 * t85;
t454 = -t229 * t95 - t233 * t94;
t453 = qJD(4) - t172;
t282 = t233 * Ifges(4,1) - t376;
t451 = Ifges(4,5) * qJD(3) + t160 * Ifges(6,5) + t80 * Ifges(7,5) + t159 * Ifges(6,6) + Ifges(7,6) * t293 + t194 * Ifges(6,3) + t186 * Ifges(7,3) + qJD(1) * t282;
t272 = t229 * Ifges(5,3) - t371;
t368 = t160 * Ifges(6,4);
t68 = t159 * Ifges(6,2) + t194 * Ifges(6,6) + t368;
t151 = Ifges(6,4) * t159;
t69 = t160 * Ifges(6,1) + t194 * Ifges(6,5) + t151;
t450 = Ifges(5,5) * qJD(3) + qJD(1) * t272 + t228 * t69 + t232 * t68;
t75 = Ifges(7,4) * t293;
t449 = -Ifges(7,2) * t80 + t75;
t235 = -pkin(9) - pkin(8);
t345 = t229 * t235;
t363 = t229 * mrSges(5,2);
t448 = -mrSges(3,3) + mrSges(2,2) - m(7) * (-pkin(5) * t348 - t216 - t345) - m(6) * t271 + t363 - (-m(5) * qJ(4) - mrSges(5,3)) * t233 - t467;
t391 = pkin(5) * t232;
t203 = pkin(4) + t391;
t447 = m(7) * t203 + mrSges(2,1) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t416 = t23 / 0.2e1;
t415 = t24 / 0.2e1;
t410 = t72 / 0.2e1;
t409 = t73 / 0.2e1;
t402 = t148 / 0.2e1;
t401 = t158 / 0.2e1;
t381 = pkin(4) - t237;
t174 = t380 * t228;
t90 = -t174 * t231 - t176 * t227;
t443 = -qJD(6) * t90 + t227 * t469 + t231 * t470;
t89 = t174 * t227 - t176 * t231;
t442 = qJD(6) * t89 + t227 * t470 - t231 * t469;
t127 = mrSges(5,1) * t170 - qJDD(3) * mrSges(5,3);
t38 = -mrSges(6,1) * t73 + mrSges(6,2) * t72;
t437 = -t127 + t38;
t262 = -t343 + t350;
t251 = qJD(3) * t262;
t436 = -t229 * t251 + t263 * t468;
t252 = qJD(3) * t263;
t435 = t229 * t252 + t262 * t468;
t87 = -mrSges(6,1) * t159 + mrSges(6,2) * t160;
t434 = -t181 + t87;
t433 = pkin(5) * t325 + t203 * t212 + t453;
t118 = t262 * t229;
t218 = t229 * pkin(3);
t334 = t218 - t216;
t290 = -t334 - t390;
t149 = qJ(2) - t290;
t177 = t381 * t233;
t152 = t228 * t177;
t81 = t232 * t149 + t152;
t429 = pkin(3) * t329 + qJ(4) * t330;
t393 = pkin(5) * t228;
t428 = -m(7) * t393 + t456;
t234 = cos(qJ(1));
t388 = g(2) * t234;
t230 = sin(qJ(1));
t389 = g(1) * t230;
t425 = t388 - t389;
t128 = -t169 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t422 = -t128 + t466;
t421 = -t13 * t431 - t14 * t432 - t2 * t263 + t262 * t3;
t420 = qJD(1) ^ 2;
t419 = Ifges(7,4) * t416 + Ifges(7,2) * t415 + Ifges(7,6) * t402;
t418 = Ifges(7,1) * t416 + Ifges(7,4) * t415 + Ifges(7,5) * t402;
t414 = Ifges(6,1) * t410 + Ifges(6,4) * t409 + Ifges(6,5) * t401;
t33 = Ifges(7,2) * t293 + Ifges(7,6) * t186 + t395;
t413 = -t33 / 0.2e1;
t34 = Ifges(7,1) * t80 + Ifges(7,5) * t186 + t75;
t412 = -t34 / 0.2e1;
t411 = t34 / 0.2e1;
t408 = -t293 / 0.2e1;
t407 = t293 / 0.2e1;
t405 = t80 / 0.2e1;
t399 = t160 / 0.2e1;
t397 = t186 / 0.2e1;
t394 = pkin(5) * t160;
t379 = mrSges(7,2) * t213;
t375 = Ifges(4,4) * t233;
t374 = Ifges(6,4) * t228;
t373 = Ifges(6,4) * t232;
t372 = Ifges(5,6) * t229;
t349 = t228 * t229;
t347 = t229 * t230;
t346 = t229 * t232;
t205 = t229 * t237;
t344 = t230 * t233;
t342 = t232 * t233;
t341 = t233 * t234;
t111 = t213 * t230 - t214 * t341;
t112 = -t213 * t341 - t214 * t230;
t339 = -t111 * mrSges(7,1) + t112 * mrSges(7,2);
t113 = -t213 * t234 - t214 * t344;
t114 = -t213 * t344 + t214 * t234;
t338 = t113 * mrSges(7,1) - t114 * mrSges(7,2);
t333 = t234 * pkin(1) + t230 * qJ(2);
t328 = qJD(3) * t237;
t326 = qJD(5) * t229;
t322 = qJDD(1) * mrSges(3,2);
t39 = -mrSges(7,1) * t293 + mrSges(7,2) * t80;
t314 = -t39 - t434;
t312 = Ifges(6,5) * t72 + Ifges(6,6) * t73 + Ifges(6,3) * t158;
t309 = t234 * pkin(7) + t333;
t193 = t233 * t328;
t301 = -t327 / 0.2e1;
t198 = qJ(4) + t393;
t299 = -pkin(9) * t229 - t149;
t298 = -t319 / 0.2e1;
t156 = t381 * t330;
t97 = qJD(2) + (qJD(3) * pkin(8) - qJD(4)) * t233 + t429;
t296 = -t232 * t156 - t228 * t97;
t294 = (t188 + t320) * qJ(2);
t157 = -pkin(4) * t329 + t193;
t289 = mrSges(4,1) * t233 - mrSges(4,2) * t229;
t287 = mrSges(6,1) * t232 - mrSges(6,2) * t228;
t284 = -t233 * mrSges(5,2) + t229 * mrSges(5,3);
t283 = -t233 * mrSges(5,3) - t363;
t281 = Ifges(6,1) * t232 - t374;
t280 = Ifges(6,1) * t228 + t373;
t279 = -t229 * Ifges(4,2) + t375;
t277 = -Ifges(6,2) * t228 + t373;
t276 = Ifges(6,2) * t232 + t374;
t274 = Ifges(6,5) * t232 - Ifges(6,6) * t228;
t273 = Ifges(6,5) * t228 + Ifges(6,6) * t232;
t153 = t232 * t177;
t62 = pkin(5) * t233 + t228 * t299 + t153;
t66 = pkin(9) * t346 + t81;
t27 = -t227 * t66 + t231 * t62;
t28 = t227 * t62 + t231 * t66;
t132 = -qJD(3) * pkin(3) + t453;
t141 = -t171 - t224;
t264 = t132 * t229 - t141 * t233;
t166 = t283 * qJD(1);
t259 = t166 + t426;
t258 = -t228 * t230 + t232 * t341;
t145 = -t228 * t234 - t230 * t342;
t59 = -pkin(4) * t170 - t82;
t257 = qJ(2) * t289;
t134 = -qJ(4) * t212 + t335;
t256 = t134 * t284;
t254 = t233 * (-Ifges(4,1) * t229 - t375);
t35 = -t149 * t327 - t228 * t156 + t177 * t325 + t232 * t97;
t250 = -t228 * t326 + t232 * t329;
t249 = t228 * t329 + t229 * t325;
t245 = -Ifges(6,5) * t229 + t233 * t280;
t244 = -Ifges(6,6) * t229 + t233 * t276;
t243 = -Ifges(6,3) * t229 + t233 * t273;
t240 = qJD(3) * t264 + t455;
t217 = t234 * qJ(2);
t209 = -qJDD(1) * pkin(1) + qJDD(2);
t204 = Ifges(5,6) * t332;
t184 = t229 * t214 * mrSges(7,1);
t175 = -pkin(4) * t229 + t205;
t173 = qJ(2) + t334;
t167 = t288 * qJD(1);
t150 = t287 * t229;
t146 = -t228 * t344 + t232 * t234;
t144 = -t228 * t341 - t230 * t232;
t138 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t212 + t204;
t135 = Ifges(4,6) * qJD(3) + qJD(1) * t279;
t126 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t170;
t125 = qJDD(3) * mrSges(4,1) + mrSges(4,3) * t169;
t124 = -t203 * t229 + t205;
t121 = t263 * t233;
t120 = t263 * t229;
t119 = t262 * t233;
t115 = t292 + t429;
t86 = -pkin(5) * t250 + t157;
t78 = -t149 * t228 + t153;
t65 = pkin(3) * t170 + t242;
t61 = mrSges(7,1) * t186 - mrSges(7,3) * t80;
t60 = -mrSges(7,2) * t186 + mrSges(7,3) * t293;
t46 = -t229 * t84 - t233 * t251;
t44 = -t118 * t315 + t233 * t252;
t36 = -qJD(5) * t81 + t296;
t31 = -pkin(5) * t73 + t59;
t29 = t72 * Ifges(6,4) + t73 * Ifges(6,2) + t158 * Ifges(6,6);
t26 = pkin(9) * t250 + t35;
t25 = t261 * qJD(3) + (t232 * t299 - t152) * qJD(5) + t296;
t18 = -mrSges(7,2) * t148 + mrSges(7,3) * t24;
t17 = mrSges(7,1) * t148 - mrSges(7,3) * t23;
t16 = t231 * t40 - t367;
t15 = -t227 * t40 - t358;
t8 = -mrSges(7,1) * t24 + mrSges(7,2) * t23;
t5 = -qJD(6) * t28 - t227 * t26 + t231 * t25;
t4 = qJD(6) * t27 + t227 * t25 + t231 * t26;
t1 = [-t455 * mrSges(5,1) + (-m(3) * t333 - m(4) * t309 - t146 * mrSges(6,1) - t114 * mrSges(7,1) - t145 * mrSges(6,2) - t113 * mrSges(7,2) - t321 * (pkin(3) * t347 + t309) + (-m(6) * pkin(4) - t447) * t234 + t448 * t230) * g(2) + (-t144 * mrSges(6,1) - t112 * mrSges(7,1) + t258 * mrSges(6,2) - t111 * mrSges(7,2) + t404 * t217 - t321 * (t234 * t218 + t217) + (m(3) * pkin(1) + t381 * m(6) + (-m(7) - m(4) - m(5)) * t237 + t447) * t230 + t448 * t234) * g(1) + m(4) * (-t237 * t454 + t294) + (Ifges(7,5) * t44 + Ifges(7,6) * t46) * t397 + (Ifges(7,4) * t120 - Ifges(7,2) * t118 + Ifges(7,6) * t233) * t415 + (Ifges(7,1) * t120 - Ifges(7,4) * t118 + Ifges(7,5) * t233) * t416 + t31 * (mrSges(7,1) * t118 + mrSges(7,2) * t120) + (Ifges(7,5) * t120 - Ifges(7,6) * t118 + Ifges(7,3) * t233) * t402 + (t229 * (Ifges(5,3) * t233 + t372) + t254) * t319 / 0.2e1 + (qJD(5) * t69 + t29) * t346 / 0.2e1 + (-Ifges(4,1) * t169 - Ifges(4,4) * t170 + Ifges(4,5) * qJDD(3) + t312 + t313) * t233 / 0.2e1 + t460 * t229 * t328 + (t288 + 0.2e1 * mrSges(3,3)) * t188 + (Ifges(7,4) * t44 + Ifges(7,2) * t46) * t407 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) - t169 * t282 / 0.2e1 + t65 * t283 + t170 * t272 / 0.2e1 - t170 * t279 / 0.2e1 + (-t118 * t2 - t120 * t3 - t13 * t44 + t14 * t46) * mrSges(7,3) + (t138 / 0.2e1 - t51 * mrSges(6,1) + t52 * mrSges(6,2) - t451 / 0.2e1 - t132 * mrSges(5,1) - Ifges(7,3) * t397 - Ifges(7,5) * t405 - Ifges(7,6) * t407 + t465) * t330 + (t11 * t346 - t12 * t349 - t249 * t51 + t250 * t52) * mrSges(6,3) + m(5) * (t115 * t134 + t173 * t65 + t237 * t240) + t229 * t68 * t301 + (-t135 / 0.2e1 + t450 / 0.2e1 + t141 * mrSges(5,1)) * t329 + t233 * t472 + t233 * t474 + (-t128 + t125) * t233 * t237 + (Ifges(7,1) * t44 + Ifges(7,4) * t46) * t405 + t169 * t372 / 0.2e1 - t233 * (Ifges(5,4) * qJDD(3) + Ifges(5,6) * t170) / 0.2e1 + (-t127 + t126) * t205 - t169 * Ifges(5,2) * t233 + (Ifges(6,6) * t233 + t229 * t276) * t409 + (Ifges(6,5) * t233 + t229 * t280) * t410 + t44 * t411 + t349 * t414 + t120 * t418 - t118 * t419 + t159 * (qJD(3) * t244 + t277 * t326) / 0.2e1 + t194 * (qJD(3) * t243 + t274 * t326) / 0.2e1 - pkin(1) * t322 + t229 * (Ifges(5,5) * qJDD(3) + Ifges(5,6) * t169 + Ifges(5,3) * t170) / 0.2e1 - t229 * (-Ifges(4,4) * t169 - Ifges(4,2) * t170 + Ifges(4,6) * qJDD(3)) / 0.2e1 + m(7) * (t124 * t31 + t13 * t5 + t14 * t4 + t2 * t28 + t27 * t3 + t74 * t86) + m(6) * (t102 * t157 + t11 * t81 + t12 * t78 + t175 * t59 + t35 * t52 + t36 * t51) + t102 * (-mrSges(6,1) * t250 + mrSges(6,2) * t249) + t209 * mrSges(3,2) + t175 * t38 + qJ(2) * (mrSges(4,1) * t170 - mrSges(4,2) * t169) + t173 * (-mrSges(5,2) * t170 + mrSges(5,3) * t169) + t115 * t166 + qJD(2) * t167 + t157 * t87 - t59 * t150 + t124 * t8 - t233 * t473 - t233 * t475 + t454 * mrSges(4,3) + t457 * qJD(3) ^ 2 / 0.2e1 + t458 * t298 + t461 * t193 + qJD(3) * t256 + t257 * t319 + (qJD(3) * t245 + t281 * t326) * t399 + (Ifges(6,3) * t233 + t229 * t273) * t401 + (t229 * t462 - t233 * t463) * qJDD(3) / 0.2e1 + t27 * t17 + t28 * t18 + t46 * t33 / 0.2e1 + t4 * t60 + t5 * t61 + t74 * (-mrSges(7,1) * t46 + mrSges(7,2) * t44) + m(3) * (-pkin(1) * t209 + t294) + t78 * t53 + t81 * t54 + t86 * t39 + t35 * t98 + t36 * t99; t322 + t119 * t17 - t121 * t18 + t436 * t61 + t435 * t60 + (qJ(2) * t404 - mrSges(3,3)) * t420 + (-m(5) * t134 - t167 - t259) * qJD(1) + (t126 + t8 + (t228 * t98 + t232 * t99 + t460) * qJD(3) + t437) * t229 + (t125 + (t179 - t314) * qJD(3) + t422) * t233 + m(3) * t209 - m(4) * t454 + m(5) * t240 + t425 * (t321 - t404) + (t119 * t3 - t121 * t2 + t13 * t436 + t14 * t435 + t229 * t31 + t329 * t74) * m(7) + ((t59 + (t228 * t52 + t232 * t51) * qJD(3)) * t229 + t459 * t233 + t269 * qJD(1)) * m(6); t194 * t287 * t102 - t460 * t171 + (qJ(4) * t59 + t102 * t424 - t51 * t63 - t52 * t64) * m(6) + t421 * mrSges(7,3) + (Ifges(7,4) * t110 + Ifges(7,2) * t109) * t408 - (Ifges(5,3) * t212 + t138 + t204) * t332 / 0.2e1 - (t159 * t276 + t160 * t280 + t194 * t273) * qJD(5) / 0.2e1 - (t159 * t244 + t160 * t245 + t194 * t243) * qJD(1) / 0.2e1 - t263 * t419 + t31 * (mrSges(7,1) * t263 - mrSges(7,2) * t262) + (-Ifges(7,5) * t262 - Ifges(7,6) * t263) * t402 + (-Ifges(7,4) * t262 - Ifges(7,2) * t263) * t415 + (-Ifges(7,1) * t262 - Ifges(7,4) * t263) * t416 - t262 * t418 + (-Ifges(7,5) * t84 - Ifges(7,6) * t83) * t397 + (-Ifges(7,1) * t84 - Ifges(7,4) * t83) * t405 + (-Ifges(7,4) * t84 - Ifges(7,2) * t83) * t407 + (Ifges(7,1) * t110 + Ifges(7,4) * t109) * t406 + (-m(6) * t290 - m(7) * (-t334 + t345) + m(5) * t334 + t283 + t428 * t233 + t467) * g(3) + (-pkin(3) * t85 - qJ(4) * t82 - qJD(4) * t141 - t134 * t165 - t187 * t264) * m(5) + (-t51 * (-mrSges(6,1) * t229 - mrSges(6,3) * t348) - t52 * (mrSges(6,2) * t229 + mrSges(6,3) * t342) - t256) * qJD(1) + (-t325 * t52 + t327 * t51 - t270) * mrSges(6,3) + t59 * t286 + (-Ifges(7,5) * t406 - Ifges(7,6) * t408 - Ifges(7,3) * t398 - t465) * t332 + (mrSges(7,1) * t432 + mrSges(7,2) * t431) * t74 + t432 * t413 + t433 * t39 + t434 * qJD(4) + t437 * qJ(4) + (-t257 - t254 / 0.2e1 + t458 / 0.2e1) * t420 + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + t132 * t311 + (-m(6) * t241 + t466) * t403 + (-t321 * (pkin(3) * t344 + qJ(4) * t347) + (-t284 + (-m(6) * pkin(8) + m(7) * t235 - t444) * t233 + t428 * t229) * t230) * g(1) + t442 * t60 + (t13 * t443 + t14 * t442 + t198 * t31 + t2 * t90 + t3 * t89 + t433 * t74) * m(7) + t443 * t61 - t289 * t389 + t277 * t409 + t281 * t410 - t84 * t411 + t110 * t412 + t232 * t414 + (Ifges(7,5) * t110 + Ifges(7,6) * t109) * t398 - t68 * t325 / 0.2e1 + t135 * t212 / 0.2e1 - t141 * t310 - t228 * t29 / 0.2e1 + t198 * t8 - t165 * t166 - t123 * t87 - pkin(3) * t128 - t450 * t212 / 0.2e1 + t451 * t332 / 0.2e1 + (t284 + t289 + (m(6) * t403 - m(7) * (-pkin(3) + t235) + m(5) * pkin(3) + t444) * t233 + (m(7) * t198 + qJ(4) * t445 - t456) * t229) * t388 + t457 * t298 - t461 * t172 + t69 * t301 + t274 * t401 + t462 * t170 + t463 * t169 - t82 * mrSges(5,3) + t85 * mrSges(5,2) + t89 * t17 + t90 * t18 + t94 * mrSges(4,1) - t95 * mrSges(4,2) - t64 * t98 - t63 * t99; t263 * t18 - t262 * t17 + t431 * t61 + t432 * t60 + t314 * qJD(3) + t259 * t212 + (-g(3) * t229 - t233 * t425) * t321 + (-qJD(3) * t74 - t421) * m(7) + (-t212 * t269 - t459) * m(6) + (qJD(3) * t141 + t134 * t212 + t85) * m(5) - t422; (-mrSges(6,2) * t144 + t258 * t471 - t339) * g(2) + (-t184 - (m(7) * t391 - t379) * t229 - t150) * g(3) - (-Ifges(6,2) * t160 + t151 + t69) * t159 / 0.2e1 + t293 * t412 - t80 * t413 + (t377 + t99) * t52 + t464 + (t378 - t98) * t51 + t312 - t39 * t394 - m(7) * (t13 * t15 + t14 * t16 + t394 * t74) + (mrSges(6,2) * t146 + t145 * t471 - t338) * g(1) - t160 * (Ifges(6,1) * t159 - t368) / 0.2e1 + (t2 * t227 + t231 * t3 + (-t13 * t227 + t14 * t231) * qJD(6)) * t417 + ((qJD(6) * t60 + t17) * t231 + t227 * t18 - t323 * t61) * pkin(5) - t194 * (Ifges(6,5) * t159 - Ifges(6,6) * t160) / 0.2e1 - t102 * (mrSges(6,1) * t160 + mrSges(6,2) * t159) + t449 * t408 + t68 * t399 - t473 + t472 - t16 * t60 - t15 * t61; t33 * t405 - t13 * t60 + t14 * t61 - g(1) * t338 - g(2) * t339 - g(3) * (-t229 * t379 + t184) + (t34 + t449) * t408 + t464;];
tau  = t1;
