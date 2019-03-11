% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR2
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:59
% EndTime: 2019-03-09 03:37:25
% DurationCPUTime: 18.00s
% Computational Cost: add. (11181->723), mult. (23952->951), div. (0->0), fcn. (16917->18), ass. (0->334)
t238 = sin(qJ(5));
t335 = sin(pkin(11));
t291 = t335 * pkin(3);
t212 = t291 + pkin(8);
t357 = pkin(9) + t212;
t287 = qJD(5) * t357;
t239 = sin(qJ(3));
t243 = cos(qJ(3));
t282 = qJD(1) * t335;
t336 = cos(pkin(11));
t283 = qJD(1) * t336;
t180 = -t239 * t282 + t243 * t283;
t330 = t180 * t238;
t181 = -t239 * t283 - t243 * t282;
t312 = qJD(1) * t239;
t298 = pkin(3) * t312;
t114 = -pkin(4) * t181 - pkin(8) * t180 + t298;
t242 = cos(qJ(5));
t234 = sin(pkin(10));
t213 = pkin(1) * t234 + pkin(7);
t203 = t213 * qJD(1);
t278 = qJ(4) * qJD(1) + t203;
t310 = qJD(2) * t239;
t160 = t243 * t278 + t310;
t142 = t335 * t160;
t226 = t243 * qJD(2);
t159 = -t239 * t278 + t226;
t94 = t159 * t336 - t142;
t55 = t238 * t114 + t242 * t94;
t441 = -pkin(9) * t330 + t238 * t287 + t55;
t329 = t180 * t242;
t54 = t242 * t114 - t238 * t94;
t440 = pkin(5) * t181 + pkin(9) * t329 - t242 * t287 - t54;
t359 = t242 * pkin(5);
t217 = pkin(4) + t359;
t233 = qJ(5) + qJ(6);
t227 = sin(t233);
t228 = cos(t233);
t273 = -mrSges(6,1) * t242 + mrSges(6,2) * t238;
t439 = -m(6) * pkin(4) - m(7) * t217 - t228 * mrSges(7,1) + t227 * mrSges(7,2) + t273;
t307 = qJD(5) * t238;
t438 = t307 - t330;
t241 = cos(qJ(6));
t174 = qJD(5) - t180;
t158 = qJD(3) * t238 - t181 * t242;
t235 = cos(pkin(10));
t215 = -pkin(1) * t235 - pkin(2);
t229 = t243 * pkin(3);
t200 = t215 - t229;
t179 = qJD(1) * t200 + qJD(4);
t102 = -pkin(4) * t180 + pkin(8) * t181 + t179;
t149 = qJD(3) * pkin(3) + t159;
t286 = t336 * t160;
t86 = t335 * t149 + t286;
t82 = qJD(3) * pkin(8) + t86;
t50 = t242 * t102 - t238 * t82;
t35 = -pkin(9) * t158 + t50;
t34 = pkin(5) * t174 + t35;
t237 = sin(qJ(6));
t157 = qJD(3) * t242 + t181 * t238;
t51 = t102 * t238 + t242 * t82;
t36 = pkin(9) * t157 + t51;
t341 = t237 * t36;
t13 = t241 * t34 - t341;
t339 = t241 * t36;
t14 = t237 * t34 + t339;
t417 = Ifges(5,6) * qJD(3);
t437 = t179 * mrSges(5,1) + t50 * mrSges(6,1) + t13 * mrSges(7,1) - t51 * mrSges(6,2) - t14 * mrSges(7,2) - t417 / 0.2e1;
t418 = Ifges(5,5) * qJD(3);
t436 = -t179 * mrSges(5,2) - t418 / 0.2e1;
t305 = m(5) + m(6) + m(7);
t435 = -mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t201 = t213 * qJDD(1);
t432 = qJD(2) * qJD(3) + t201;
t207 = -mrSges(4,1) * t243 + mrSges(4,2) * t239;
t231 = qJ(3) + pkin(11);
t220 = sin(t231);
t222 = cos(t231);
t431 = -t222 * mrSges(5,1) - t220 * t435 + t207;
t398 = m(7) * pkin(5);
t280 = t241 * t157 - t158 * t237;
t304 = qJD(1) * qJD(3);
t196 = qJDD(1) * t243 - t239 * t304;
t197 = qJDD(1) * t239 + t243 * t304;
t137 = t196 * t335 + t197 * t336;
t76 = qJD(5) * t157 + qJDD(3) * t238 + t137 * t242;
t77 = -qJD(5) * t158 + qJDD(3) * t242 - t137 * t238;
t28 = qJD(6) * t280 + t237 * t77 + t241 * t76;
t397 = t28 / 0.2e1;
t92 = t157 * t237 + t158 * t241;
t29 = -qJD(6) * t92 - t237 * t76 + t241 * t77;
t396 = t29 / 0.2e1;
t389 = t76 / 0.2e1;
t388 = t77 / 0.2e1;
t430 = m(3) + m(4);
t136 = t196 * t336 - t335 * t197;
t135 = qJDD(5) - t136;
t133 = qJDD(6) + t135;
t383 = t133 / 0.2e1;
t382 = t135 / 0.2e1;
t189 = t357 * t238;
t190 = t357 * t242;
t131 = -t189 * t237 + t190 * t241;
t429 = -qJD(6) * t131 + t237 * t441 + t440 * t241;
t130 = -t189 * t241 - t190 * t237;
t428 = qJD(6) * t130 + t440 * t237 - t241 * t441;
t169 = qJD(6) + t174;
t425 = t174 * Ifges(6,3);
t426 = t157 * Ifges(6,6);
t427 = t158 * Ifges(6,5) + t92 * Ifges(7,5) + Ifges(7,6) * t280 + t169 * Ifges(7,3) + t425 + t426;
t173 = Ifges(5,4) * t180;
t423 = t180 * Ifges(5,2);
t337 = qJDD(3) / 0.2e1;
t422 = mrSges(6,1) + t398;
t93 = t159 * t335 + t286;
t421 = pkin(5) * t438 - t93;
t127 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t137;
t41 = -mrSges(6,1) * t77 + mrSges(6,2) * t76;
t420 = t41 - t127;
t345 = t181 * mrSges(5,3);
t419 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t157 + mrSges(6,2) * t158 - t345;
t192 = t239 * t336 + t243 * t335;
t261 = t237 * t238 - t241 * t242;
t123 = t261 * t192;
t315 = qJ(4) + t213;
t187 = t315 * t239;
t188 = t315 * t243;
t125 = -t187 * t335 + t188 * t336;
t116 = t242 * t125;
t191 = t239 * t335 - t243 * t336;
t117 = pkin(4) * t191 - pkin(8) * t192 + t200;
t65 = t238 * t117 + t116;
t110 = t261 * t180;
t411 = qJD(5) + qJD(6);
t139 = t411 * t261;
t416 = -t139 + t110;
t195 = t237 * t242 + t238 * t241;
t109 = t195 * t180;
t140 = t411 * t195;
t415 = -t140 + t109;
t183 = t191 * qJD(3);
t306 = qJD(5) * t242;
t257 = -t183 * t238 + t192 * t306;
t309 = qJD(3) * t239;
t120 = t239 * qJDD(2) - t203 * t309 + t243 * t432;
t171 = t203 * t243 + t310;
t225 = t243 * qJDD(2);
t121 = -qJD(3) * t171 - t201 * t239 + t225;
t414 = t120 * t243 - t121 * t239;
t56 = mrSges(6,1) * t135 - mrSges(6,3) * t76;
t57 = -mrSges(6,2) * t135 + mrSges(6,3) * t77;
t413 = -t238 * t56 + t242 * t57;
t303 = qJD(1) * qJD(4);
t308 = qJD(3) * t243;
t83 = -t203 * t308 + qJDD(3) * pkin(3) - qJ(4) * t197 + t225 + (-t303 - t432) * t239;
t87 = qJ(4) * t196 + t243 * t303 + t120;
t46 = t335 * t83 + t336 * t87;
t44 = qJDD(3) * pkin(8) + t46;
t202 = t215 * qJDD(1);
t156 = -pkin(3) * t196 + qJDD(4) + t202;
t62 = -pkin(4) * t136 - pkin(8) * t137 + t156;
t11 = t102 * t306 + t238 * t62 + t242 * t44 - t307 * t82;
t12 = -qJD(5) * t51 - t238 * t44 + t242 * t62;
t412 = t11 * t242 - t12 * t238;
t410 = 0.2e1 * t337;
t85 = t149 * t336 - t142;
t81 = -qJD(3) * pkin(4) - t85;
t63 = -t157 * pkin(5) + t81;
t408 = -mrSges(7,1) * t63 + t14 * mrSges(7,3);
t407 = mrSges(7,2) * t63 - mrSges(7,3) * t13;
t6 = pkin(5) * t135 - pkin(9) * t76 + t12;
t9 = pkin(9) * t77 + t11;
t2 = qJD(6) * t13 + t237 * t6 + t241 * t9;
t3 = -qJD(6) * t14 - t237 * t9 + t241 * t6;
t405 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t404 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t403 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t402 = m(4) * pkin(2) + mrSges(3,1) - t431;
t400 = Ifges(7,4) * t397 + Ifges(7,2) * t396 + Ifges(7,6) * t383;
t399 = Ifges(7,1) * t397 + Ifges(7,4) * t396 + Ifges(7,5) * t383;
t395 = Ifges(6,1) * t389 + Ifges(6,4) * t388 + Ifges(6,5) * t382;
t368 = Ifges(7,4) * t92;
t39 = Ifges(7,2) * t280 + Ifges(7,6) * t169 + t368;
t394 = -t39 / 0.2e1;
t393 = t39 / 0.2e1;
t88 = Ifges(7,4) * t280;
t40 = Ifges(7,1) * t92 + Ifges(7,5) * t169 + t88;
t392 = -t40 / 0.2e1;
t391 = t40 / 0.2e1;
t387 = -t280 / 0.2e1;
t386 = t280 / 0.2e1;
t385 = -t92 / 0.2e1;
t384 = t92 / 0.2e1;
t381 = -t157 / 0.2e1;
t380 = -t158 / 0.2e1;
t379 = t158 / 0.2e1;
t378 = -t169 / 0.2e1;
t377 = t169 / 0.2e1;
t376 = -t174 / 0.2e1;
t375 = -t180 / 0.2e1;
t374 = -t181 / 0.2e1;
t370 = t242 / 0.2e1;
t240 = sin(qJ(1));
t367 = pkin(1) * t240;
t365 = pkin(5) * t158;
t364 = pkin(8) * t220;
t361 = g(3) * t220;
t244 = cos(qJ(1));
t230 = t244 * pkin(1);
t356 = mrSges(6,3) * t157;
t355 = mrSges(6,3) * t158;
t354 = Ifges(4,4) * t239;
t353 = Ifges(4,4) * t243;
t352 = Ifges(6,4) * t238;
t351 = Ifges(6,4) * t242;
t347 = t158 * Ifges(6,4);
t346 = t180 * mrSges(5,3);
t344 = t181 * Ifges(5,4);
t327 = t183 * t242;
t326 = t192 * t238;
t325 = t192 * t242;
t245 = -pkin(9) - pkin(8);
t324 = t220 * t245;
t232 = qJ(1) + pkin(10);
t221 = sin(t232);
t323 = t221 * t227;
t322 = t221 * t228;
t321 = t221 * t238;
t320 = t221 * t242;
t223 = cos(t232);
t319 = t223 * t227;
t318 = t223 * t228;
t317 = t223 * t238;
t316 = t223 * t242;
t151 = t222 * t323 + t318;
t152 = -t222 * t322 + t319;
t314 = -t151 * mrSges(7,1) + t152 * mrSges(7,2);
t153 = -t222 * t319 + t322;
t154 = t222 * t318 + t323;
t313 = t153 * mrSges(7,1) - t154 * mrSges(7,2);
t311 = qJD(1) * t243;
t52 = -mrSges(7,1) * t280 + mrSges(7,2) * t92;
t301 = t52 + t419;
t300 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t133;
t299 = Ifges(6,5) * t76 + Ifges(6,6) * t77 + Ifges(6,3) * t135;
t297 = pkin(3) * t309;
t295 = mrSges(4,3) * t312;
t294 = mrSges(4,3) * t311;
t150 = Ifges(6,4) * t157;
t73 = t158 * Ifges(6,1) + t174 * Ifges(6,5) + t150;
t293 = t73 * t370;
t292 = t336 * pkin(3);
t288 = -t307 / 0.2e1;
t284 = -t136 * mrSges(5,1) + t137 * mrSges(5,2);
t279 = qJD(3) * t315;
t161 = qJD(4) * t243 - t239 * t279;
t162 = -qJD(4) * t239 - t243 * t279;
t101 = t161 * t336 + t162 * t335;
t182 = t192 * qJD(3);
t115 = pkin(4) * t182 + pkin(8) * t183 + t297;
t281 = -t101 * t238 + t242 * t115;
t64 = t242 * t117 - t125 * t238;
t214 = -t292 - pkin(4);
t277 = pkin(4) * t222 + t364;
t45 = -t335 * t87 + t336 * t83;
t275 = mrSges(4,1) * t239 + mrSges(4,2) * t243;
t272 = mrSges(6,1) * t238 + mrSges(6,2) * t242;
t271 = -mrSges(7,1) * t227 - mrSges(7,2) * t228;
t270 = Ifges(6,1) * t242 - t352;
t269 = t243 * Ifges(4,2) + t354;
t268 = -Ifges(6,2) * t238 + t351;
t267 = Ifges(4,5) * t243 - Ifges(4,6) * t239;
t266 = Ifges(6,5) * t242 - Ifges(6,6) * t238;
t53 = pkin(5) * t191 - pkin(9) * t325 + t64;
t58 = -pkin(9) * t326 + t65;
t23 = -t237 * t58 + t241 * t53;
t24 = t237 * t53 + t241 * t58;
t264 = -t238 * t50 + t242 * t51;
t100 = t161 * t335 - t336 * t162;
t124 = t336 * t187 + t188 * t335;
t103 = -mrSges(6,2) * t174 + t356;
t104 = mrSges(6,1) * t174 - t355;
t263 = t103 * t242 - t104 * t238;
t262 = t217 * t222 - t324;
t260 = t300 + t405;
t165 = -t222 * t317 + t320;
t163 = t222 * t321 + t316;
t167 = -qJD(3) * mrSges(5,2) + t346;
t259 = t167 + t263;
t258 = t81 * t272;
t256 = t192 * t307 + t327;
t255 = t215 * qJD(1) * t275;
t254 = t239 * (Ifges(4,1) * t243 - t354);
t30 = t242 * t101 + t238 * t115 + t117 * t306 - t125 * t307;
t43 = -qJDD(3) * pkin(4) - t45;
t249 = (-t51 * t238 - t50 * t242) * qJD(5) + t412;
t236 = -qJ(4) - pkin(7);
t219 = Ifges(4,4) * t311;
t218 = t229 + pkin(2);
t206 = -qJD(3) * mrSges(4,2) + t294;
t204 = qJD(3) * mrSges(4,1) - t295;
t199 = t214 - t359;
t185 = Ifges(4,1) * t312 + Ifges(4,5) * qJD(3) + t219;
t184 = Ifges(4,6) * qJD(3) + qJD(1) * t269;
t176 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t197;
t175 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t196;
t170 = -t203 * t239 + t226;
t166 = t222 * t316 + t321;
t164 = -t222 * t320 + t317;
t128 = -mrSges(5,1) * t180 - mrSges(5,2) * t181;
t126 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t136;
t122 = t195 * t192;
t119 = -t181 * Ifges(5,1) + t173 + t418;
t118 = -t344 + t417 + t423;
t98 = pkin(5) * t326 + t124;
t72 = t157 * Ifges(6,2) + t174 * Ifges(6,6) + t347;
t67 = mrSges(7,1) * t169 - mrSges(7,3) * t92;
t66 = -mrSges(7,2) * t169 + mrSges(7,3) * t280;
t61 = pkin(5) * t257 + t100;
t48 = t123 * t411 + t195 * t183;
t47 = -t140 * t192 + t183 * t261;
t32 = t76 * Ifges(6,4) + t77 * Ifges(6,2) + t135 * Ifges(6,6);
t31 = -qJD(5) * t65 + t281;
t27 = -t77 * pkin(5) + t43;
t22 = -pkin(9) * t257 + t30;
t21 = -mrSges(7,2) * t133 + mrSges(7,3) * t29;
t20 = mrSges(7,1) * t133 - mrSges(7,3) * t28;
t19 = pkin(9) * t327 + pkin(5) * t182 + (-t116 + (pkin(9) * t192 - t117) * t238) * qJD(5) + t281;
t16 = t241 * t35 - t341;
t15 = -t237 * t35 - t339;
t10 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t5 = -qJD(6) * t24 + t19 * t241 - t22 * t237;
t4 = qJD(6) * t23 + t19 * t237 + t22 * t241;
t1 = [(-m(5) * t85 + m(6) * t81 + t419) * t100 + (-m(5) * t45 + m(6) * t43 + t420) * t124 + (m(4) * ((-t170 * t243 - t171 * t239) * qJD(3) + t414) - t204 * t308 - t206 * t309 + t243 * t175 - t239 * t176) * t213 + (-t170 * t308 - t171 * t309 + t414) * mrSges(4,3) - t257 * t72 / 0.2e1 + (-Ifges(6,1) * t256 - Ifges(6,4) * t257) * t379 + (t243 * (-Ifges(4,2) * t239 + t353) + t254) * t304 / 0.2e1 + (t300 + t299) * t191 / 0.2e1 + t197 * t239 * Ifges(4,1) + (t156 * mrSges(5,2) - t45 * mrSges(5,3) + Ifges(5,1) * t137 + Ifges(5,4) * t136 + Ifges(5,5) * t410 + t266 * t382 + t268 * t388 + t270 * t389 + t43 * t272 + t73 * t288) * t192 + (t156 * mrSges(5,1) - t46 * mrSges(5,3) - Ifges(5,4) * t137 + Ifges(6,5) * t389 + Ifges(7,5) * t397 - Ifges(5,2) * t136 - t410 * Ifges(5,6) + Ifges(6,6) * t388 + Ifges(7,6) * t396 + Ifges(6,3) * t382 + Ifges(7,3) * t383 + t404 + t405) * t191 + t128 * t297 + (-t86 * mrSges(5,3) + t425 / 0.2e1 + t426 / 0.2e1 - t118 / 0.2e1 - t423 / 0.2e1 - Ifges(5,4) * t374 + Ifges(7,3) * t377 + Ifges(6,5) * t379 + Ifges(7,5) * t384 + Ifges(7,6) * t386 + t427 / 0.2e1 + t437) * t182 + (t85 * mrSges(5,3) - t173 / 0.2e1 - t119 / 0.2e1 - t293 - Ifges(5,1) * t374 + t436) * t183 + t197 * t353 / 0.2e1 + t239 * (Ifges(4,4) * t196 + Ifges(4,5) * qJDD(3)) / 0.2e1 + t200 * t284 + t196 * t269 / 0.2e1 + m(7) * (t13 * t5 + t14 * t4 + t2 * t24 + t23 * t3 + t27 * t98 + t61 * t63) - t32 * t326 / 0.2e1 + (t255 + t267 * qJD(3) / 0.2e1) * qJD(3) + t81 * (mrSges(6,1) * t257 - mrSges(6,2) * t256) + t185 * t308 / 0.2e1 - t184 * t309 / 0.2e1 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t384 + (-t11 * t326 - t12 * t325 + t256 * t50 - t257 * t51) * mrSges(6,3) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t235 - 0.2e1 * mrSges(3,2) * t234 + m(3) * (t234 ^ 2 + t235 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + m(6) * (t11 * t65 + t12 * t64 + t30 * t51 + t31 * t50) + m(5) * (t101 * t86 + t125 * t46 + t156 * t200 + t179 * t297) + (m(4) * t215 + t207) * t202 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t377 + t174 * (-Ifges(6,5) * t256 - Ifges(6,6) * t257) / 0.2e1 + t243 * (Ifges(4,4) * t197 + Ifges(4,2) * t196 + Ifges(4,6) * qJDD(3)) / 0.2e1 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t386 + t215 * (-mrSges(4,1) * t196 + mrSges(4,2) * t197) + t101 * t167 + t157 * (-Ifges(6,4) * t256 - Ifges(6,2) * t257) / 0.2e1 + t125 * t126 + t30 * t103 + t31 * t104 + t98 * t10 + t63 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + t64 * t56 + t65 * t57 + t4 * t66 + t5 * t67 + t61 * t52 + t23 * t20 + t24 * t21 + (-t321 * t398 - mrSges(2,1) * t244 - t166 * mrSges(6,1) - t154 * mrSges(7,1) + mrSges(2,2) * t240 - t165 * mrSges(6,2) - t153 * mrSges(7,2) - t305 * (t223 * t218 - t221 * t236 + t230) - t430 * t230 + t403 * t221 + (-m(6) * t277 - m(7) * t262 - t402) * t223) * g(2) + (-t317 * t398 + mrSges(2,1) * t240 - t164 * mrSges(6,1) - t152 * mrSges(7,1) + mrSges(2,2) * t244 - t163 * mrSges(6,2) - t151 * mrSges(7,2) + t430 * t367 - t305 * (-t223 * t236 - t367) + t403 * t223 + (-m(7) * (-t218 - t262) - m(6) * (-t218 - t277) + m(5) * t218 + t402) * t221) * g(1) + (Ifges(4,5) * t239 + Ifges(4,6) * t243) * t337 + t47 * t391 + t48 * t393 + t325 * t395 - t123 * t399 - t122 * t400 + (-t122 * t2 + t123 * t3 - t13 * t47 + t14 * t48) * mrSges(7,3) + (-Ifges(7,1) * t123 - Ifges(7,4) * t122) * t397 + (-Ifges(7,5) * t123 - Ifges(7,6) * t122) * t383 + (-Ifges(7,4) * t123 - Ifges(7,2) * t122) * t396 + t27 * (mrSges(7,1) * t122 - mrSges(7,2) * t123); m(3) * qJDD(2) - t122 * t20 - t123 * t21 + t239 * t175 + t243 * t176 + t47 * t66 + t48 * t67 + (-t204 * t239 + t206 * t243) * qJD(3) + (t10 + t420) * t191 - t259 * t183 + t301 * t182 + (t126 + (-t238 * t103 - t242 * t104) * qJD(5) + t413) * t192 + (-t305 - t430) * g(3) + m(4) * (t120 * t239 + t121 * t243 + (-t170 * t239 + t171 * t243) * qJD(3)) + m(6) * (t182 * t81 - t183 * t264 + t191 * t43 + t192 * t249) + m(7) * (-t122 * t3 - t123 * t2 + t13 * t48 + t14 * t47 + t182 * t63 + t191 * t27) + m(5) * (-t182 * t85 - t183 * t86 - t191 * t45 + t192 * t46); -t419 * t93 + t421 * t52 + (-t179 * t298 + t85 * t93 - t86 * t94 + (t335 * t46 + t336 * t45) * pkin(3)) * m(5) - (Ifges(7,1) * t384 + Ifges(7,4) * t386 + Ifges(7,5) * t377 + t391 + t407) * t139 - (Ifges(7,4) * t384 + Ifges(7,2) * t386 + Ifges(7,6) * t377 + t393 + t408) * t140 + (m(6) * t249 - t103 * t307 - t104 * t306 + t413) * t212 + (-t255 - t254 * qJD(1) / 0.2e1) * qJD(1) + (t295 + t204) * t171 + (t330 / 0.2e1 + t288) * t72 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) - (-Ifges(4,2) * t312 + t185 + t219) * t311 / 0.2e1 + (t157 * t268 + t158 * t270 + t174 * t266) * qJD(5) / 0.2e1 + (t214 * t43 - t50 * t54 - t51 * t55 - t81 * t93) * m(6) - t267 * t304 / 0.2e1 + (-Ifges(6,5) * t380 - Ifges(7,5) * t385 + Ifges(5,2) * t375 - Ifges(6,6) * t381 - Ifges(7,6) * t387 - Ifges(6,3) * t376 - Ifges(7,3) * t378 + t437) * t181 + (-t438 * t51 + (-t306 + t329) * t50 + t412) * mrSges(6,3) + (g(1) * t223 + g(2) * t221) * (t275 + t305 * pkin(3) * t239 + (-m(6) * pkin(8) + m(7) * t245 - t435) * t222 + (mrSges(5,1) - t439) * t220) + (-m(7) * (t229 - t324) - m(6) * (t229 + t364) - m(5) * t229 + t439 * t222 + t431) * g(3) + (t266 * t376 + t268 * t381 + t270 * t380 - t258 + t436) * t180 - t86 * t345 + t184 * t312 / 0.2e1 + t43 * t273 + (t294 - t206) * t170 - t73 * t329 / 0.2e1 + t27 * (mrSges(7,1) * t261 + mrSges(7,2) * t195) + (Ifges(7,5) * t195 - Ifges(7,6) * t261) * t383 + (Ifges(7,4) * t195 - Ifges(7,2) * t261) * t396 + (Ifges(7,1) * t195 - Ifges(7,4) * t261) * t397 + (t109 * t14 - t110 * t13 - t195 * t3 - t2 * t261) * mrSges(7,3) - t261 * t400 + (t258 + t293) * qJD(5) - t128 * t298 + (t173 + t119) * t375 + t214 * t41 + t199 * t10 + Ifges(4,6) * t196 + Ifges(4,5) * t197 - t94 * t167 + Ifges(5,6) * t136 + Ifges(5,5) * t137 + t130 * t20 + t131 * t21 - t120 * mrSges(4,2) + t121 * mrSges(4,1) - t55 * t103 - t54 * t104 + t45 * mrSges(5,1) - t46 * mrSges(5,2) + (Ifges(5,1) * t180 + t344 + t427) * t181 / 0.2e1 + t428 * t66 + t429 * t67 + (t13 * t429 + t130 * t3 + t131 * t2 + t14 * t428 + t199 * t27 + t421 * t63) * m(7) + t126 * t291 + t127 * t292 + t85 * t346 + t32 * t370 + t118 * t374 + (Ifges(6,5) * t238 + Ifges(6,6) * t242) * t382 + (Ifges(6,2) * t242 + t352) * t388 + (Ifges(6,1) * t238 + t351) * t389 - t110 * t392 - t109 * t394 + t238 * t395 + t195 * t399 + (-Ifges(7,1) * t110 - Ifges(7,4) * t109) * t385 + (-Ifges(7,4) * t110 - Ifges(7,2) * t109) * t387 + (-Ifges(7,5) * t110 - Ifges(7,6) * t109) * t378 - t63 * (mrSges(7,1) * t109 - mrSges(7,2) * t110); -t261 * t20 + t195 * t21 + t238 * t57 + t242 * t56 + t415 * t67 + t416 * t66 + t263 * qJD(5) + t301 * t181 - t259 * t180 + t284 + (-g(1) * t221 + g(2) * t223) * t305 + (t13 * t415 + t14 * t416 + t181 * t63 + t195 * t2 - t261 * t3) * m(7) + (t11 * t238 + t12 * t242 + t174 * t264 + t181 * t81) * m(6) + (-t180 * t86 - t181 * t85 + t156) * m(5); (mrSges(6,2) * t166 - t165 * t422 - t313) * g(1) + t404 + t299 - (Ifges(7,4) * t385 + Ifges(7,2) * t387 + Ifges(7,6) * t378 + t394 - t408) * t92 + (Ifges(7,1) * t385 + Ifges(7,4) * t387 + Ifges(7,5) * t378 + t392 - t407) * t280 - t52 * t365 - m(7) * (t13 * t15 + t14 * t16 + t365 * t63) + (-mrSges(6,2) * t164 + t163 * t422 - t314) * g(2) + (t104 + t355) * t51 + (-t103 + t356) * t50 + (t238 * t398 - t271 + t272) * t361 + t260 - t81 * (mrSges(6,1) * t158 + mrSges(6,2) * t157) - t16 * t66 - t15 * t67 + (Ifges(6,5) * t157 - Ifges(6,6) * t158) * t376 + t72 * t379 + (Ifges(6,1) * t157 - t347) * t380 + (t2 * t237 + t241 * t3 + (-t13 * t237 + t14 * t241) * qJD(6)) * t398 + (-Ifges(6,2) * t158 + t150 + t73) * t381 + ((-t237 * t67 + t241 * t66) * qJD(6) + t241 * t20 + t237 * t21) * pkin(5); -t63 * (mrSges(7,1) * t92 + mrSges(7,2) * t280) + (Ifges(7,1) * t280 - t368) * t385 + t39 * t384 + (Ifges(7,5) * t280 - Ifges(7,6) * t92) * t378 - t13 * t66 + t14 * t67 - g(1) * t313 - g(2) * t314 - t271 * t361 + (t13 * t280 + t14 * t92) * mrSges(7,3) + t260 + (-Ifges(7,2) * t92 + t40 + t88) * t387;];
tau  = t1;
