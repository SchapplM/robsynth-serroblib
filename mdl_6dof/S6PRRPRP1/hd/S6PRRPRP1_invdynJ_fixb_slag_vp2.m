% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:41
% EndTime: 2019-03-08 21:23:19
% DurationCPUTime: 22.69s
% Computational Cost: add. (7162->682), mult. (16643->894), div. (0->0), fcn. (12623->14), ass. (0->324)
t477 = Ifges(6,4) + Ifges(7,4);
t240 = cos(qJ(5));
t395 = t240 * pkin(5);
t478 = -m(6) * pkin(4) - m(7) * (pkin(4) + t395) - mrSges(5,1);
t238 = sin(qJ(3));
t241 = cos(qJ(3));
t335 = qJD(2) * qJD(3);
t198 = qJDD(2) * t241 - t238 * t335;
t199 = qJDD(2) * t238 + t241 * t335;
t232 = sin(pkin(11));
t367 = cos(pkin(11));
t136 = t232 * t198 + t199 * t367;
t291 = t367 * t238;
t340 = qJD(2) * t241;
t184 = -qJD(2) * t291 - t232 * t340;
t237 = sin(qJ(5));
t149 = qJD(3) * t240 + t184 * t237;
t67 = qJD(5) * t149 + qJDD(3) * t237 + t136 * t240;
t415 = t67 / 0.2e1;
t150 = qJD(3) * t237 - t184 * t240;
t68 = -qJD(5) * t150 + qJDD(3) * t240 - t136 * t237;
t414 = t68 / 0.2e1;
t461 = Ifges(6,1) + Ifges(7,1);
t459 = -Ifges(7,5) - Ifges(6,5);
t458 = Ifges(6,2) + Ifges(7,2);
t457 = Ifges(6,6) + Ifges(7,6);
t456 = Ifges(7,3) + Ifges(6,3);
t252 = -t232 * t238 + t241 * t367;
t186 = t252 * qJD(3);
t194 = t232 * t241 + t291;
t336 = qJD(5) * t240;
t312 = t194 * t336;
t256 = t237 * t186 + t312;
t337 = qJD(5) * t237;
t183 = t252 * qJD(2);
t361 = t183 * t237;
t476 = t337 - t361;
t135 = t198 * t367 - t232 * t199;
t129 = qJDD(5) - t135;
t474 = t477 * t414 + t461 * t415 - t459 * t129 / 0.2e1;
t473 = t477 * t149;
t242 = cos(qJ(2));
t233 = sin(pkin(6));
t343 = qJD(2) * t233;
t305 = qJD(1) * t343;
t213 = t242 * t305;
t239 = sin(qJ(2));
t333 = qJDD(1) * t233;
t173 = t239 * t333 + t213;
t163 = qJDD(2) * pkin(8) + t173;
t234 = cos(pkin(6));
t344 = qJD(1) * t234;
t472 = qJD(3) * t344 + t163;
t366 = sin(pkin(10));
t290 = t366 * t239;
t368 = cos(pkin(10));
t293 = t368 * t242;
t182 = -t234 * t290 + t293;
t295 = t233 * t366;
t471 = -t182 * t238 + t241 * t295;
t354 = t233 * t239;
t187 = t234 * t241 - t238 * t354;
t470 = t477 * t150;
t469 = t477 * t240;
t468 = t477 * t237;
t467 = qJD(5) - t183;
t231 = qJ(3) + pkin(11);
t229 = sin(t231);
t230 = cos(t231);
t235 = -qJ(6) - pkin(9);
t466 = t478 * t230 + (-m(6) * pkin(9) + m(7) * t235 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3)) * t229;
t416 = m(7) * pkin(5);
t465 = t198 / 0.2e1;
t464 = t129 * t457 + t458 * t68 + t477 * t67;
t463 = -mrSges(7,1) - mrSges(6,1);
t462 = -mrSges(6,2) - mrSges(7,2);
t454 = t149 * t457 - t150 * t459 + t456 * t467;
t453 = t149 * t458 + t457 * t467 + t470;
t452 = t150 * t461 - t459 * t467 + t473;
t384 = mrSges(7,3) * t149;
t87 = -mrSges(7,2) * t467 + t384;
t386 = mrSges(6,3) * t149;
t88 = -mrSges(6,2) * t467 + t386;
t390 = t87 + t88;
t383 = mrSges(7,3) * t150;
t89 = mrSges(7,1) * t467 - t383;
t385 = mrSges(6,3) * t150;
t90 = mrSges(6,1) * t467 - t385;
t389 = t89 + t90;
t171 = Ifges(5,4) * t183;
t451 = t183 * Ifges(5,2);
t369 = qJDD(3) / 0.2e1;
t399 = pkin(3) * t232;
t224 = pkin(9) + t399;
t348 = qJ(6) + t224;
t286 = qJD(5) * t348;
t342 = qJD(2) * t238;
t327 = pkin(3) * t342;
t103 = -pkin(4) * t184 - pkin(9) * t183 + t327;
t345 = qJD(1) * t233;
t317 = t239 * t345;
t202 = qJD(2) * pkin(8) + t317;
t285 = qJ(4) * qJD(2) + t202;
t315 = t238 * t344;
t141 = t241 * t285 + t315;
t130 = t232 * t141;
t221 = t241 * t344;
t140 = -t238 * t285 + t221;
t73 = t140 * t367 - t130;
t31 = t237 * t103 + t240 * t73;
t450 = qJ(6) * t361 + qJD(6) * t240 - t237 * t286 - t31;
t30 = t240 * t103 - t237 * t73;
t360 = t183 * t240;
t449 = pkin(5) * t184 + qJ(6) * t360 - qJD(6) * t237 - t240 * t286 - t30;
t387 = mrSges(5,3) * t184;
t448 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t149 + mrSges(6,2) * t150 - t387;
t115 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t136;
t24 = -mrSges(6,1) * t68 + mrSges(6,2) * t67;
t447 = t24 - t115;
t292 = t367 * t141;
t71 = t140 * t232 + t292;
t446 = pkin(5) * t476 - t71;
t445 = t416 + mrSges(7,1);
t444 = Ifges(5,5) * qJD(3);
t443 = Ifges(5,6) * qJD(3);
t227 = pkin(3) * t241 + pkin(2);
t119 = -pkin(4) * t252 - pkin(9) * t194 - t227;
t236 = -qJ(4) - pkin(8);
t208 = t236 * t238;
t209 = t236 * t241;
t146 = t232 * t208 - t209 * t367;
t137 = t240 * t146;
t58 = t237 * t119 + t137;
t297 = qJD(3) * t236;
t177 = qJD(4) * t241 + t238 * t297;
t178 = -qJD(4) * t238 + t241 * t297;
t102 = t177 * t367 + t232 * t178;
t316 = t242 * t345;
t144 = t252 * t316;
t442 = -t144 + t102;
t272 = mrSges(7,1) * t237 + mrSges(7,2) * t240;
t274 = mrSges(6,1) * t237 + mrSges(6,2) * t240;
t134 = qJD(3) * pkin(3) + t140;
t60 = t134 * t367 - t130;
t55 = -qJD(3) * pkin(4) - t60;
t40 = -t149 * pkin(5) + qJD(6) + t55;
t441 = t40 * t272 + t55 * t274;
t440 = -t237 * t457 - t240 * t459;
t439 = -t237 * t458 + t469;
t438 = t240 * t461 - t468;
t289 = t366 * t242;
t294 = t368 * t239;
t180 = t234 * t294 + t289;
t296 = t233 * t368;
t437 = -t180 * t238 - t241 * t296;
t339 = qJD(3) * t238;
t326 = pkin(3) * t339;
t436 = t326 - t317;
t435 = t129 * t456 + t457 * t68 - t459 * t67;
t434 = -t336 + t360;
t332 = qJDD(1) * t234;
t75 = -t202 * t339 + t238 * t332 + t241 * t472;
t155 = t202 * t241 + t315;
t218 = t241 * t332;
t76 = -t155 * qJD(3) - t163 * t238 + t218;
t432 = -t238 * t76 + t241 * t75;
t334 = qJD(2) * qJD(4);
t338 = qJD(3) * t241;
t47 = -t202 * t338 + qJDD(3) * pkin(3) - qJ(4) * t199 + t218 + (-t334 - t472) * t238;
t48 = qJ(4) * t198 + t241 * t334 + t75;
t17 = t232 * t47 + t367 * t48;
t15 = qJDD(3) * pkin(9) + t17;
t212 = t239 * t305;
t172 = t242 * t333 - t212;
t162 = -qJDD(2) * pkin(2) - t172;
t117 = -pkin(3) * t198 + qJDD(4) + t162;
t39 = -pkin(4) * t135 - pkin(9) * t136 + t117;
t61 = t232 * t134 + t292;
t56 = qJD(3) * pkin(9) + t61;
t170 = -qJD(2) * t227 + qJD(4) - t316;
t84 = -pkin(4) * t183 + pkin(9) * t184 + t170;
t3 = t240 * t15 + t237 * t39 + t84 * t336 - t337 * t56;
t29 = t237 * t84 + t240 * t56;
t4 = -qJD(5) * t29 - t15 * t237 + t240 * t39;
t431 = -t237 * t4 + t240 * t3;
t429 = m(5) + m(7) + m(6);
t428 = -mrSges(6,1) - t445;
t324 = m(4) * pkin(8) + mrSges(4,3);
t427 = mrSges(3,2) - mrSges(5,3) - t324;
t426 = 0.2e1 * t369;
t273 = -mrSges(7,1) * t240 + mrSges(7,2) * t237;
t275 = -mrSges(6,1) * t240 + mrSges(6,2) * t237;
t425 = -t273 - t275 - t478;
t82 = -mrSges(7,1) * t149 + mrSges(7,2) * t150;
t424 = -m(7) * t40 - t82;
t423 = m(5) * t60 - t448;
t277 = -t241 * mrSges(4,1) + t238 * mrSges(4,2);
t251 = m(4) * pkin(2) - t277;
t422 = mrSges(3,1) + t251 - t466;
t421 = t170 * mrSges(5,2) + t444 / 0.2e1;
t2 = qJ(6) * t68 + qJD(6) * t149 + t3;
t420 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t123 = t180 * t230 - t229 * t296;
t125 = t182 * t230 + t229 * t295;
t166 = t229 * t234 + t230 * t354;
t419 = -g(1) * t125 - g(2) * t123 - g(3) * t166;
t21 = qJ(6) * t149 + t29;
t28 = -t237 * t56 + t240 * t84;
t20 = -qJ(6) * t150 + t28;
t8 = pkin(5) * t467 + t20;
t418 = t170 * mrSges(5,1) + t28 * mrSges(6,1) + t8 * mrSges(7,1) - t443 / 0.2e1 - t21 * mrSges(7,2) - t29 * mrSges(6,2);
t243 = qJD(2) ^ 2;
t417 = m(5) * pkin(3);
t413 = t129 / 0.2e1;
t412 = -t149 / 0.2e1;
t411 = t149 / 0.2e1;
t410 = -t150 / 0.2e1;
t409 = t150 / 0.2e1;
t408 = -t467 / 0.2e1;
t407 = t467 / 0.2e1;
t406 = -t183 / 0.2e1;
t405 = -t184 / 0.2e1;
t404 = t184 / 0.2e1;
t397 = g(3) * t233;
t34 = mrSges(7,1) * t129 - mrSges(7,3) * t67;
t35 = mrSges(6,1) * t129 - mrSges(6,3) * t67;
t392 = t34 + t35;
t36 = -mrSges(7,2) * t129 + mrSges(7,3) * t68;
t37 = -mrSges(6,2) * t129 + mrSges(6,3) * t68;
t391 = t36 + t37;
t388 = mrSges(5,3) * t183;
t382 = Ifges(4,4) * t238;
t381 = Ifges(4,4) * t241;
t380 = Ifges(5,4) * t184;
t365 = t180 * t237;
t363 = t182 * t237;
t359 = t186 * t240;
t358 = t194 * t237;
t357 = t194 * t240;
t356 = t230 * t237;
t355 = t230 * t240;
t353 = t233 * t242;
t350 = t237 * t242;
t349 = t240 * t242;
t341 = qJD(2) * t239;
t185 = t194 * qJD(3);
t104 = pkin(4) * t185 - pkin(9) * t186 + t326;
t331 = t240 * t102 + t237 * t104 + t119 * t336;
t330 = t82 + t448;
t323 = mrSges(4,3) * t342;
t322 = mrSges(4,3) * t340;
t321 = t233 * t350;
t319 = t233 * t349;
t318 = t367 * pkin(3);
t314 = t233 * t341;
t313 = t242 * t343;
t23 = -t68 * mrSges(7,1) + t67 * mrSges(7,2);
t299 = -t337 / 0.2e1;
t69 = -t135 * mrSges(5,1) + t136 * mrSges(5,2);
t287 = -t102 * t237 + t240 * t104;
t57 = t240 * t119 - t146 * t237;
t101 = t177 * t232 - t367 * t178;
t145 = -t367 * t208 - t209 * t232;
t282 = t471 * pkin(3);
t225 = -t318 - pkin(4);
t16 = -t232 * t48 + t367 * t47;
t269 = t241 * Ifges(4,2) + t382;
t266 = Ifges(4,5) * t241 - Ifges(4,6) * t238;
t261 = -qJ(6) * t186 - qJD(6) * t194;
t260 = t187 * pkin(3);
t188 = t234 * t238 + t241 * t354;
t109 = t232 * t187 + t188 * t367;
t85 = -t109 * t237 - t319;
t259 = -t109 * t240 + t321;
t255 = t194 * t337 - t359;
t203 = -qJD(2) * pkin(2) - t316;
t254 = t203 * (mrSges(4,1) * t238 + mrSges(4,2) * t241);
t253 = t238 * (Ifges(4,1) * t241 - t382);
t14 = -qJDD(3) * pkin(4) - t16;
t179 = -t234 * t293 + t290;
t181 = t234 * t289 + t294;
t247 = -g(1) * t181 - g(2) * t179 + g(3) * t353;
t245 = t437 * pkin(3);
t228 = Ifges(4,4) * t340;
t207 = -qJD(3) * mrSges(4,2) + t322;
t206 = qJD(3) * mrSges(4,1) - t323;
t204 = t225 - t395;
t196 = t277 * qJD(2);
t192 = t348 * t240;
t191 = t348 * t237;
t190 = Ifges(4,1) * t342 + Ifges(4,5) * qJD(3) + t228;
t189 = Ifges(4,6) * qJD(3) + qJD(2) * t269;
t176 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t199;
t175 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t198;
t165 = t229 * t354 - t234 * t230;
t160 = -qJD(3) * mrSges(5,2) + t388;
t154 = -t202 * t238 + t221;
t142 = -mrSges(4,1) * t198 + mrSges(4,2) * t199;
t139 = qJD(3) * t187 + t241 * t313;
t138 = -qJD(3) * t188 - t238 * t313;
t124 = t182 * t229 - t230 * t295;
t122 = t180 * t229 + t230 * t296;
t116 = -mrSges(5,1) * t183 - mrSges(5,2) * t184;
t114 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t135;
t111 = -t184 * Ifges(5,1) + t171 + t444;
t110 = -t380 + t443 + t451;
t108 = -t187 * t367 + t188 * t232;
t106 = t144 * t240 + t237 * t317;
t105 = -t144 * t237 + t240 * t317;
t100 = pkin(5) * t358 + t145;
t72 = t232 * t138 + t139 * t367;
t70 = -t138 * t367 + t139 * t232;
t46 = pkin(5) * t256 + t101;
t41 = -qJ(6) * t358 + t58;
t38 = -pkin(5) * t252 - qJ(6) * t357 + t57;
t27 = qJD(5) * t259 - t237 * t72 + t240 * t314;
t26 = qJD(5) * t85 + t237 * t314 + t240 * t72;
t19 = -qJD(5) * t58 + t287;
t18 = -t146 * t337 + t331;
t7 = -qJ(6) * t312 + (-qJD(5) * t146 + t261) * t237 + t331;
t6 = pkin(5) * t185 + t261 * t240 + (-t137 + (qJ(6) * t194 - t119) * t237) * qJD(5) + t287;
t5 = -t68 * pkin(5) + qJDD(6) + t14;
t1 = pkin(5) * t129 - qJ(6) * t67 - qJD(6) * t150 + t4;
t9 = [m(2) * qJDD(1) + t109 * t114 + t138 * t206 + t139 * t207 + t72 * t160 + t188 * t175 + t187 * t176 - t391 * t259 + t392 * t85 + t389 * t27 + t390 * t26 + t330 * t70 + (t23 + t447) * t108 + (-m(2) - m(3) - m(4) - t429) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t243 - t142 - t69) * t242 + (-mrSges(3,1) * t243 - mrSges(3,2) * qJDD(2) + (t116 + t196) * qJD(2)) * t239) * t233 + m(5) * (-t108 * t16 + t109 * t17 - t60 * t70 + t61 * t72 + (-t117 * t242 + t170 * t341) * t233) + m(4) * (t138 * t154 + t139 * t155 + t187 * t76 + t188 * t75 + (-t162 * t242 + t203 * t341) * t233) + m(3) * (qJDD(1) * t234 ^ 2 + (t172 * t242 + t173 * t239) * t233) + m(7) * (t1 * t85 + t108 * t5 - t2 * t259 + t21 * t26 + t27 * t8 + t40 * t70) + m(6) * (t108 * t14 - t259 * t3 + t26 * t29 + t27 * t28 + t4 * t85 + t55 * t70); -t256 * t453 / 0.2e1 + (Ifges(5,1) * t405 + t171 / 0.2e1 + t111 / 0.2e1 + t421 - t60 * mrSges(5,3)) * t186 + (t454 / 0.2e1 - Ifges(5,4) * t405 - t451 / 0.2e1 - t110 / 0.2e1 + t457 * t411 - t459 * t409 + t456 * t407 + t418 - t61 * mrSges(5,3)) * t185 + (-t101 * t60 - t117 * t227 - t145 * t16 + t146 * t17 + t170 * t436 + t442 * t61) * m(5) + t442 * t160 + (t241 * (-Ifges(4,2) * t238 + t381) + t253) * t335 / 0.2e1 + (-t206 * t338 - t207 * t339 + t241 * t175 - t238 * t176 + m(4) * ((-t154 * t241 - t155 * t238) * qJD(3) + t432)) * pkin(8) + (-t154 * t338 - t155 * t339 + t432) * mrSges(4,3) + t436 * t116 + (t101 * t55 + t14 * t145 + t3 * t58 + t4 * t57 + (-t106 + t18) * t29 + (-t105 + t19) * t28) * m(6) + (t1 * t38 + t100 * t5 + t2 * t41 + t40 * t46 + (-t105 + t6) * t8 + (-t106 + t7) * t21) * m(7) + (-t429 * t227 * t353 + (t463 * (t230 * t349 + t237 * t239) + t462 * (-t230 * t350 + t239 * t240) + t466 * t242 + (t236 * t429 - t237 * t416 - mrSges(5,3)) * t239) * t233) * g(3) + (-t435 / 0.2e1 - t117 * mrSges(5,1) - t1 * mrSges(7,1) + t17 * mrSges(5,3) + Ifges(5,4) * t136 + Ifges(5,2) * t135 + Ifges(5,6) * t426 - t413 * t456 - t414 * t457 + t415 * t459 - t420) * t252 + (t255 * t459 - t256 * t457) * t407 - t464 * t358 / 0.2e1 + (Ifges(4,1) * t199 + Ifges(4,4) * t465 + Ifges(4,5) * t426 + t206 * t316) * t238 + (-t365 * t416 + t463 * (-t179 * t355 + t365) + t462 * (t179 * t356 + t180 * t240) - t429 * (-t179 * t227 - t180 * t236) + t427 * t180 + t422 * t179) * g(2) + t447 * t145 + t448 * t101 - t241 * t207 * t316 + (t254 + t266 * qJD(3) / 0.2e1) * qJD(3) + (t255 * t28 - t256 * t29 - t3 * t358 - t357 * t4) * mrSges(6,3) + t199 * t381 / 0.2e1 + t269 * t465 + t162 * t277 + (-t363 * t416 + t463 * (-t181 * t355 + t363) + t462 * (t181 * t356 + t182 * t240) - t429 * (-t181 * t227 - t182 * t236) + t427 * t182 + t422 * t181) * g(1) + ((-m(6) * t55 + t423 + t424) * t316 + t117 * mrSges(5,2) - t16 * mrSges(5,3) + Ifges(5,1) * t136 + Ifges(5,4) * t135 + Ifges(5,5) * t426 + t14 * t274 + t5 * t272 + t413 * t440 + t414 * t439 + t415 * t438 + t452 * t299) * t194 + t40 * (mrSges(7,1) * t256 - mrSges(7,2) * t255) + t55 * (mrSges(6,1) * t256 - mrSges(6,2) * t255) + (-(t203 * t239 + (-t154 * t238 + t155 * t241) * t242) * t345 - pkin(2) * t162) * m(4) + (-t255 * t477 - t256 * t458) * t411 + (-t255 * t461 - t256 * t477) * t409 - t389 * t105 - t390 * t106 + t452 * t359 / 0.2e1 + (-t242 * t397 + t172 + t212) * mrSges(3,1) + (-t1 * t357 - t2 * t358 - t21 * t256 + t255 * t8) * mrSges(7,3) + (t239 * t397 - t173 + t213) * mrSges(3,2) + Ifges(3,3) * qJDD(2) + t241 * (Ifges(4,4) * t199 + Ifges(4,2) * t198 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t227 * t69 - pkin(2) * t142 + t146 * t114 + t357 * t474 - (t239 * t324 + t242 * t251) * t397 + t100 * t23 + t7 * t87 + t18 * t88 + t6 * t89 + t19 * t90 + t46 * t82 + t57 * t35 + t58 * t37 + t41 * t36 + t38 * t34 + t190 * t338 / 0.2e1 - t189 * t339 / 0.2e1 - t196 * t317 + Ifges(4,6) * t241 * t369; (t149 * t439 + t150 * t438 + t440 * t467) * qJD(5) / 0.2e1 + (-m(6) * (t123 * pkin(9) + t245) - m(7) * (-t123 * t235 + t245) + t123 * mrSges(5,2) - (-t180 * t241 + t238 * t296) * mrSges(4,2) + (-t417 - mrSges(4,1)) * t437 + t425 * t122) * g(2) + (Ifges(5,1) * t404 + t408 * t440 + t410 * t438 + t412 * t439 - t421 - t441) * t183 + t441 * qJD(5) - (-Ifges(4,2) * t342 + t190 + t228) * t340 / 0.2e1 + (-t90 * t336 - t88 * t337 + t240 * t37 - t237 * t35 + m(6) * ((-t29 * t237 - t28 * t240) * qJD(5) + t431)) * t224 + (-m(6) * (pkin(9) * t166 + t260) - m(5) * t260 + t166 * mrSges(5,2) - m(7) * (-t166 * t235 + t260) - mrSges(4,1) * t187 + mrSges(4,2) * t188 + t425 * t165) * g(3) + t423 * t71 + (-m(6) * (pkin(9) * t125 + t282) - m(7) * (-t125 * t235 + t282) - t471 * mrSges(4,1) - (-t182 * t241 - t238 * t295) * mrSges(4,2) - m(5) * t282 + t125 * mrSges(5,2) + t425 * t124) * g(1) + (t240 * t458 + t468) * t414 + (t237 * t461 + t469) * t415 + t60 * t388 + t114 * t399 + (t380 + t454) * t404 + (Ifges(5,2) * t406 - t408 * t456 + t410 * t459 - t412 * t457 + t418) * t184 + (-t237 * t459 + t240 * t457) * t413 + t464 * t240 / 0.2e1 + t446 * t82 + t449 * t89 + (-t1 * t191 + t192 * t2 + t204 * t5 + t450 * t21 + t446 * t40 + t449 * t8) * m(7) + t450 * t87 + t5 * t273 + t14 * t275 + (t323 + t206) * t155 + (t171 + t111) * t406 - qJD(2) * t254 + (t28 * t434 - t29 * t476 + t419 + t431) * mrSges(6,3) + (-t1 * t237 + t2 * t240 - t21 * t476 + t434 * t8 + t419) * mrSges(7,3) + t115 * t318 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t110 * t405 + (t16 * t367 + t17 * t232) * t417 + (t14 * t225 - t28 * t30 - t29 * t31 - t55 * t71) * m(6) + (t322 - t207) * t154 - m(5) * (t170 * t327 + t61 * t73) + (t361 / 0.2e1 + t299) * t453 + (-t360 / 0.2e1 + t336 / 0.2e1) * t452 + t225 * t24 + Ifges(4,5) * t199 + t204 * t23 + Ifges(4,6) * t198 - t191 * t34 + t192 * t36 - t73 * t160 - t243 * t253 / 0.2e1 + Ifges(5,6) * t135 + Ifges(5,5) * t136 + t237 * t474 - t31 * t88 - t30 * t90 - t75 * mrSges(4,2) + t76 * mrSges(4,1) - t61 * t387 + t16 * mrSges(5,1) - t17 * mrSges(5,2) + t189 * t342 / 0.2e1 - t266 * t335 / 0.2e1 - t116 * t327; -t183 * t160 + t330 * t184 + (t390 * t467 + t392) * t240 + (-t389 * t467 + t391) * t237 + t69 + (t1 * t240 + t184 * t40 + t2 * t237 + t247 + t467 * (t21 * t240 - t8 * t237)) * m(7) + (t184 * t55 + t3 * t237 + t4 * t240 + t247 + t467 * (-t28 * t237 + t29 * t240)) * m(6) + (-t183 * t61 - t184 * t60 + t117 + t247) * m(5); (t462 * (-t123 * t240 - t179 * t237) + t428 * (-t123 * t237 + t179 * t240)) * g(2) + (t462 * (-t125 * t240 - t181 * t237) + t428 * (-t125 * t237 + t181 * t240)) * g(1) + (t462 * (-t166 * t240 + t321) + t428 * (-t166 * t237 - t319)) * g(3) + t8 * t384 + t420 + (t385 + t90) * t29 + (t383 - m(7) * (t20 - t8) + t89) * t21 + (-t150 * t458 + t452 + t473) * t412 + (t149 * t461 - t470) * t410 + (t386 - t88) * t28 + t453 * t409 + t445 * t1 - t55 * (mrSges(6,1) * t150 + mrSges(6,2) * t149) - t40 * (mrSges(7,1) * t150 + mrSges(7,2) * t149) - t20 * t87 + (-t149 * t459 - t150 * t457) * t408 + t435 + (t150 * t424 + t34) * pkin(5); -t149 * t87 + t150 * t89 + (-g(1) * t124 - g(2) * t122 - g(3) * t165 - t149 * t21 + t150 * t8 + t5) * m(7) + t23;];
tau  = t9;
