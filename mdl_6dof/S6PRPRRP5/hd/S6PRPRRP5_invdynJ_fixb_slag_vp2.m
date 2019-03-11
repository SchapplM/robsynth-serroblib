% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:04
% EndTime: 2019-03-08 20:14:28
% DurationCPUTime: 15.43s
% Computational Cost: add. (4098->566), mult. (8700->739), div. (0->0), fcn. (5895->10), ass. (0->271)
t423 = Ifges(6,4) + Ifges(7,4);
t424 = Ifges(6,1) + Ifges(7,1);
t402 = Ifges(6,5) + Ifges(7,5);
t422 = Ifges(6,2) + Ifges(7,2);
t401 = Ifges(6,6) + Ifges(7,6);
t171 = cos(qJ(5));
t430 = t423 * t171;
t168 = sin(qJ(5));
t429 = t423 * t168;
t296 = qJD(4) * t171;
t172 = cos(qJ(4));
t300 = qJD(2) * t172;
t134 = -t168 * t300 + t296;
t169 = sin(qJ(4));
t288 = qJD(2) * qJD(4);
t139 = qJDD(2) * t172 - t169 * t288;
t60 = qJD(5) * t134 + qJDD(4) * t168 + t139 * t171;
t367 = t60 / 0.2e1;
t298 = qJD(4) * t168;
t135 = t171 * t300 + t298;
t61 = -qJD(5) * t135 + qJDD(4) * t171 - t139 * t168;
t366 = t61 / 0.2e1;
t365 = -m(5) - m(4);
t140 = -qJDD(2) * t169 - t172 * t288;
t130 = qJDD(5) - t140;
t364 = t130 / 0.2e1;
t428 = -t134 / 0.2e1;
t427 = -t135 / 0.2e1;
t302 = qJD(2) * t169;
t161 = qJD(5) + t302;
t426 = -t161 / 0.2e1;
t425 = qJD(4) / 0.2e1;
t421 = Ifges(6,3) + Ifges(7,3);
t387 = -t168 * t401 + t171 * t402;
t385 = -t168 * t422 + t430;
t383 = t171 * t424 - t429;
t162 = pkin(5) * t171 + pkin(4);
t167 = -qJ(6) - pkin(9);
t220 = mrSges(5,1) * t169 + mrSges(5,2) * t172;
t404 = mrSges(3,2) - mrSges(4,3);
t420 = -m(7) * (t162 * t169 + t167 * t172) + t172 * mrSges(7,3) - t220 + t404;
t419 = t139 / 0.2e1;
t418 = t140 / 0.2e1;
t417 = t402 * t364 + t423 * t366 + t367 * t424;
t174 = -pkin(2) - pkin(8);
t165 = sin(pkin(6));
t170 = sin(qJ(2));
t301 = qJD(2) * t170;
t265 = t165 * t301;
t151 = qJD(1) * t265;
t173 = cos(qJ(2));
t315 = t165 * t173;
t104 = qJDD(1) * t315 - t151;
t195 = qJDD(3) - t104;
t166 = cos(pkin(6));
t303 = qJD(1) * t166;
t416 = -qJD(4) * t303 + qJDD(2) * t174 + t195;
t304 = qJD(1) * t165;
t268 = t170 * t304;
t415 = t423 * t134;
t267 = t173 * t304;
t200 = qJD(3) - t267;
t118 = qJD(2) * t174 + t200;
t287 = qJDD(1) * t166;
t295 = qJD(4) * t172;
t22 = t118 * t295 + t169 * t416 + t172 * t287;
t20 = qJDD(4) * pkin(9) + t22;
t292 = qJD(5) * t171;
t293 = qJD(5) * t168;
t285 = qJDD(2) * qJ(3);
t286 = qJDD(1) * t170;
t84 = t165 * t286 + t285 + (qJD(3) + t267) * qJD(2);
t40 = -pkin(4) * t140 - pkin(9) * t139 + t84;
t266 = t172 * t303;
t81 = t118 * t169 + t266;
t65 = qJD(4) * pkin(9) + t81;
t222 = pkin(4) * t169 - pkin(9) * t172;
t141 = qJ(3) + t222;
t92 = qJD(2) * t141 + t268;
t3 = t168 * t40 + t171 * t20 + t92 * t292 - t293 * t65;
t26 = t168 * t92 + t171 * t65;
t4 = -qJD(5) * t26 - t168 * t20 + t171 * t40;
t221 = -t168 * t4 + t171 * t3;
t25 = -t168 * t65 + t171 * t92;
t414 = -t25 * t292 - t26 * t293 + t221;
t413 = t423 * t135;
t336 = Ifges(5,4) * t169;
t215 = t172 * Ifges(5,1) - t336;
t398 = t135 * t424 + t402 * t161 + t415;
t412 = Ifges(5,5) * qJD(4) + qJD(2) * t215 + t171 * t398;
t294 = qJD(4) * t174;
t335 = Ifges(5,4) * t172;
t210 = -Ifges(5,2) * t169 + t335;
t410 = Ifges(5,6) * t425 + qJD(2) * t210 / 0.2e1 + t421 * t426 + t402 * t427 + t401 * t428;
t2 = qJ(6) * t61 + qJD(6) * t134 + t3;
t409 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t368 = m(7) * pkin(5);
t408 = m(6) + m(7);
t407 = t130 * t401 + t422 * t61 + t423 * t60;
t346 = mrSges(3,1) - mrSges(4,2);
t405 = -mrSges(6,1) - mrSges(7,1);
t403 = mrSges(6,2) + mrSges(7,2);
t399 = t134 * t422 + t401 * t161 + t413;
t338 = mrSges(7,3) * t134;
t88 = -mrSges(7,2) * t161 + t338;
t340 = mrSges(6,3) * t134;
t89 = -mrSges(6,2) * t161 + t340;
t342 = t88 + t89;
t337 = mrSges(7,3) * t135;
t90 = mrSges(7,1) * t161 - t337;
t339 = mrSges(6,3) * t135;
t91 = mrSges(6,1) * t161 - t339;
t341 = t90 + t91;
t397 = mrSges(7,1) + t368;
t277 = mrSges(5,3) * t300;
t396 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t134 + mrSges(6,2) * t135 + t277;
t16 = -mrSges(6,1) * t61 + mrSges(6,2) * t60;
t395 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t139 - t16;
t236 = qJD(5) * t167;
t310 = t169 * t171;
t223 = pkin(4) * t172 + pkin(9) * t169;
t137 = t223 * qJD(2);
t80 = t118 * t172 - t169 * t303;
t41 = t171 * t137 - t168 * t80;
t394 = -(pkin(5) * t172 + qJ(6) * t310) * qJD(2) - t41 - qJD(6) * t168 + t171 * t236;
t263 = t168 * t302;
t290 = qJD(6) * t171;
t42 = t168 * t137 + t171 * t80;
t393 = -qJ(6) * t263 + t168 * t236 + t290 - t42;
t353 = pkin(5) * t168;
t392 = -t266 - (-qJD(2) * t353 + t118) * t169 + pkin(5) * t293;
t309 = t169 * t174;
t94 = t168 * t141 + t171 * t309;
t391 = -t169 * t387 + t172 * t421;
t390 = -t169 * t385 + t172 * t401;
t389 = -t169 * t383 + t172 * t402;
t388 = t168 * t402 + t171 * t401;
t386 = t171 * t422 + t429;
t384 = t168 * t424 + t430;
t381 = t130 * t421 + t401 * t61 + t402 * t60;
t297 = qJD(4) * t169;
t23 = -t118 * t297 - t169 * t287 + t172 * t416;
t380 = t169 * t22 + t172 * t23;
t379 = -m(5) - t408;
t378 = -mrSges(5,3) - t346;
t377 = -mrSges(6,1) - t397;
t217 = -mrSges(7,1) * t171 + mrSges(7,2) * t168;
t219 = -mrSges(6,1) * t171 + mrSges(6,2) * t168;
t376 = m(6) * pkin(4) + m(7) * t162 + mrSges(5,1) - t217 - t219;
t375 = -m(6) * pkin(9) + m(7) * t167 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t64 = -qJD(4) * pkin(4) - t80;
t43 = -pkin(5) * t134 + qJD(6) + t64;
t67 = -mrSges(7,1) * t134 + mrSges(7,2) * t135;
t374 = -m(7) * t43 - t67;
t1 = pkin(5) * t130 - qJ(6) * t60 - qJD(6) * t135 + t4;
t13 = -qJ(6) * t135 + t25;
t10 = pkin(5) * t161 + t13;
t14 = qJ(6) * t134 + t26;
t372 = -t1 * t168 - t10 * t292 - t14 * t293 + t171 * t2;
t371 = -m(6) * t64 - t396;
t327 = t172 * mrSges(6,3);
t369 = -m(6) * t141 + t327 + (-m(7) + t365) * qJ(3) + t420;
t175 = qJD(2) ^ 2;
t362 = t134 / 0.2e1;
t360 = t135 / 0.2e1;
t358 = t161 / 0.2e1;
t34 = mrSges(7,1) * t130 - mrSges(7,3) * t60;
t35 = mrSges(6,1) * t130 - mrSges(6,3) * t60;
t344 = t34 + t35;
t36 = -mrSges(7,2) * t130 + mrSges(7,3) * t61;
t37 = -mrSges(6,2) * t130 + mrSges(6,3) * t61;
t343 = t36 + t37;
t21 = -qJDD(4) * pkin(4) - t23;
t325 = t172 * t21;
t323 = cos(pkin(10));
t233 = t323 * t173;
t164 = sin(pkin(10));
t318 = t164 * t170;
t112 = -t166 * t233 + t318;
t322 = t112 * t168;
t234 = t323 * t170;
t317 = t164 * t173;
t114 = t166 * t317 + t234;
t321 = t114 * t168;
t138 = qJD(2) * qJ(3) + t268;
t320 = t138 * t173;
t319 = t164 * t165;
t316 = t165 * t170;
t314 = t168 * t169;
t313 = t168 * t170;
t312 = t168 * t172;
t311 = t168 * t173;
t308 = t170 * t171;
t307 = t171 * t172;
t305 = pkin(2) * t315 + qJ(3) * t316;
t299 = qJD(2) * t173;
t291 = qJD(5) * t172;
t289 = qJDD(2) * mrSges(4,2);
t15 = -t61 * mrSges(7,1) + t60 * mrSges(7,2);
t284 = t15 - t395;
t281 = t67 + t396;
t278 = mrSges(5,3) * t302;
t272 = t169 * t315;
t271 = t168 * t309;
t131 = qJD(4) * t223 + qJD(3);
t261 = t172 * t294;
t270 = t168 * t131 + t141 * t292 + t171 * t261;
t264 = t165 * t299;
t262 = t169 * t294;
t260 = t171 * t291;
t259 = t168 * t297;
t258 = t169 * t296;
t254 = qJD(1) * t299;
t239 = -t168 * t174 + pkin(5);
t235 = t165 * t323;
t232 = -t288 / 0.2e1;
t230 = t165 * t254;
t218 = mrSges(6,1) * t168 + mrSges(6,2) * t171;
t216 = mrSges(7,1) * t168 + mrSges(7,2) * t171;
t205 = -Ifges(5,5) * t169 - Ifges(5,6) * t172;
t199 = t26 * t168 + t25 * t171;
t198 = qJ(3) * t84 + qJD(3) * t138;
t196 = t138 * t299 + t170 * t84;
t117 = t166 * t172 - t272;
t75 = -t117 * t168 + t165 * t308;
t76 = t117 * t171 + t165 * t313;
t116 = t166 * t169 + t172 * t315;
t194 = t138 * (mrSges(5,1) * t172 - mrSges(5,2) * t169);
t193 = t169 * (-Ifges(5,2) * t172 - t336);
t192 = t172 * (-Ifges(5,1) * t169 - t335);
t188 = t168 * t291 + t258;
t187 = t259 - t260;
t97 = (t169 * t308 + t311) * t165;
t96 = (-t169 * t313 + t171 * t173) * t165;
t186 = -t168 * t342 - t171 * t341;
t184 = -qJD(5) * t94 + t171 * t131;
t148 = t167 * t171;
t147 = t167 * t168;
t145 = -qJD(4) * mrSges(5,2) - t278;
t136 = t220 * qJD(2);
t133 = -qJD(2) * pkin(2) + t200;
t132 = (-t174 + t353) * t172;
t129 = t171 * t141;
t115 = -t166 * t318 + t233;
t113 = t166 * t234 + t317;
t107 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t140;
t105 = (t254 + t286) * t165;
t103 = t114 * pkin(2);
t102 = t112 * pkin(2);
t95 = -qJDD(2) * pkin(2) + t195;
t93 = t129 - t271;
t87 = qJD(1) * t97;
t86 = qJD(1) * t96;
t85 = -pkin(5) * t187 + t262;
t78 = -mrSges(5,1) * t140 + mrSges(5,2) * t139;
t74 = -qJD(4) * t272 + t166 * t295 - t172 * t265;
t73 = -qJD(4) * t116 + t169 * t265;
t72 = -t112 * t169 + t172 * t235;
t71 = t112 * t172 + t169 * t235;
t70 = t114 * t169 + t172 * t319;
t69 = -t114 * t172 + t169 * t319;
t66 = -qJ(6) * t312 + t94;
t62 = -qJ(6) * t307 + t169 * t239 + t129;
t39 = -t168 * t261 + t184;
t38 = -qJD(5) * t271 + t270;
t18 = qJD(5) * t75 + t168 * t264 + t171 * t73;
t17 = -qJD(5) * t76 - t168 * t73 + t171 * t264;
t12 = -qJ(6) * t260 + (-qJD(6) * t172 + (qJ(6) * qJD(4) - qJD(5) * t174) * t169) * t168 + t270;
t11 = qJ(6) * t258 + (qJ(6) * t293 + qJD(4) * t239 - t290) * t172 + t184;
t5 = -pkin(5) * t61 + qJDD(6) + t21;
t6 = [t117 * t107 + t73 * t145 + t343 * t76 + t344 * t75 + t342 * t18 + t341 * t17 + t281 * t74 + t284 * t116 + (-m(2) - m(3) + t365 - t408) * g(3) + m(6) * (t116 * t21 + t17 * t25 + t18 * t26 + t3 * t76 + t4 * t75 + t64 * t74) + m(7) * (t1 * t75 + t10 * t17 + t116 * t5 + t14 * t18 + t2 * t76 + t43 * t74) + m(5) * (-t116 * t23 + t117 * t22 + t73 * t81 - t74 * t80) + (t136 * t299 + t170 * t78 + m(5) * t196 + m(3) * (t104 * t173 + t105 * t170) + m(4) * (t133 * t301 - t173 * t95 + t196) + (-t170 * t346 - t173 * t404) * t175 + (-t170 * t404 + t173 * t346) * qJDD(2)) * t165 + (m(2) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t166 ^ 2) * qJDD(1); t396 * t262 + t218 * t325 + qJ(3) * t78 + t66 * t36 + t62 * t34 - (t168 * t398 + t171 * t399) * t291 / 0.2e1 + (-t1 * t307 + t10 * t188 + t14 * t187 - t2 * t312) * mrSges(7,3) + (t3 * t94 + t4 * t93 + (t297 * t64 - t325) * t174 + (t38 - t87) * t26 + (t39 - t86) * t25) * m(6) + (t1 * t62 + t132 * t5 + t2 * t66 + t43 * t85 + (t12 - t87) * t14 + (t11 - t86) * t10) * m(7) - t412 * t297 / 0.2e1 + (t187 * t26 + t188 * t25 - t3 * t312 - t307 * t4) * mrSges(6,3) + t145 * t261 + t107 * t309 + (-t358 * t388 - t360 * t384 - t362 * t386) * t291 + t193 * t232 + (t5 * t216 + (-t268 * t80 + t294 * t81) * m(5) + t387 * t364 + Ifges(5,1) * t419 + Ifges(5,4) * t418 + (-t371 - t374) * t268 + t395 * t174 + t385 * t366 + t383 * t367 + Ifges(5,5) * qJDD(4)) * t172 + (-m(4) * t305 - (m(6) * t222 - t327) * t316 + t405 * t97 - t403 * t96 + t379 * (pkin(8) * t315 + t305) + (t170 * t420 + t378 * t173 - t311 * t368) * t165) * g(3) + (t25 * mrSges(6,1) + t10 * mrSges(7,1) - t26 * mrSges(6,2) - t14 * mrSges(7,2) - t81 * mrSges(5,3) - t410) * t295 + t200 * t136 + (t297 * t80 - t380) * mrSges(5,3) + t307 * t417 + t210 * t418 + t215 * t419 + (-(t133 * t170 + t320) * t304 - pkin(2) * t95 + t198) * m(4) + t192 * t288 / 0.2e1 - pkin(2) * t289 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (qJD(2) * qJD(3) - t230 + t285 + t84) * mrSges(4,3) + t84 * t220 + (-t151 + t95) * mrSges(4,2) + t43 * (-mrSges(7,1) * t187 - mrSges(7,2) * t188) + t64 * (-mrSges(6,1) * t187 - mrSges(6,2) * t188) + (-t105 + t230) * mrSges(3,2) + (t380 * t174 - t320 * t304 + t198) * m(5) + t85 * t67 + t12 * t88 + t38 * t89 + t11 * t90 + t39 * t91 + t93 * t35 + t94 * t37 + (t1 * mrSges(7,1) + (-t268 * t81 - t294 * t80) * m(5) + t381 / 0.2e1 + t421 * t364 - t268 * t145 - Ifges(5,4) * t139 / 0.2e1 - Ifges(5,2) * t140 / 0.2e1 + t401 * t366 + t402 * t367 - Ifges(5,6) * qJDD(4) + t409) * t169 - t341 * t86 - t342 * t87 + t399 * t259 / 0.2e1 + (t205 * t425 + t391 * t358 + t389 * t360 + t390 * t362 + t194) * qJD(4) - t407 * t312 / 0.2e1 + (t322 * t368 + m(4) * t102 + t405 * (t113 * t310 - t322) - t403 * (-t112 * t171 - t113 * t314) + t379 * (-pkin(8) * t112 - t102) - t378 * t112 + t369 * t113) * g(2) + (t321 * t368 + m(4) * t103 + t405 * (t115 * t310 - t321) - t403 * (-t114 * t171 - t115 * t314) + t379 * (-pkin(8) * t114 - t103) - t378 * t114 + t369 * t115) * g(1) + t132 * t15 + (t104 + t151) * mrSges(3,1); t289 - t175 * mrSges(4,3) + m(4) * t95 + (-t136 + t365 * t138 - m(6) * t199 - m(7) * (t10 * t171 + t14 * t168) + t186) * qJD(2) + ((-t168 * t341 + t171 * t342 + t145) * qJD(4) + m(5) * (qJD(4) * t81 + t23) + m(6) * (-t25 * t298 + t26 * t296 - t21) + m(7) * (-t10 * t298 + t14 * t296 - t5) - t284) * t172 + (t107 + t343 * t171 - t344 * t168 + t281 * qJD(4) + t186 * qJD(5) + m(5) * (-qJD(4) * t80 + t22) + m(6) * (qJD(4) * t64 + t414) + m(7) * (qJD(4) * t43 + t372)) * t169 + (m(4) - t379) * (-g(1) * t114 - g(2) * t112 + g(3) * t315); -t22 * mrSges(5,2) + t23 * mrSges(5,1) - pkin(4) * t16 + (t134 * t385 + t135 * t383 + t161 * t387) * qJD(5) / 0.2e1 - (t134 * t390 + t135 * t389 + t161 * t391) * qJD(2) / 0.2e1 + t412 * t302 / 0.2e1 + t414 * mrSges(6,3) + t392 * t67 + t393 * t88 + t394 * t90 + (t1 * t147 + t10 * t394 + t14 * t393 - t148 * t2 - t162 * t5 + t392 * t43) * m(7) + t384 * t367 + t386 * t366 + t388 * t364 + t161 * (t216 * t43 + t218 * t64) + (t116 * t376 + t117 * t375) * g(3) + (t375 * t70 + t376 * t69) * g(1) + (-t375 * t72 - t376 * t71) * g(2) + t205 * t232 + (-t91 * t292 - t89 * t293 + m(6) * (-qJD(5) * t199 + t221) - t168 * t35 + t171 * t37) * pkin(9) + (-t263 / 0.2e1 - t293 / 0.2e1) * t399 + t410 * t300 + (t277 + t371) * t81 + t372 * mrSges(7,3) + t168 * t417 + Ifges(5,3) * qJDD(4) + t5 * t217 + t21 * t219 + (t193 / 0.2e1 - t192 / 0.2e1) * t175 + (-pkin(4) * t21 - t25 * t41 - t26 * t42) * m(6) + (-t278 - t145) * t80 + (-t10 * (mrSges(7,1) * t172 + mrSges(7,3) * t310) - t25 * (mrSges(6,1) * t172 + mrSges(6,3) * t310) - t14 * (-mrSges(7,2) * t172 + mrSges(7,3) * t314) - t26 * (-mrSges(6,2) * t172 + mrSges(6,3) * t314) - t194) * qJD(2) - t42 * t89 - t41 * t91 + t398 * t292 / 0.2e1 + t407 * t171 / 0.2e1 + Ifges(5,5) * t139 + Ifges(5,6) * t140 + t147 * t34 - t148 * t36 - t162 * t15; t399 * t360 + (t337 - m(7) * (-t10 + t13) + t90) * t14 + t10 * t338 + t397 * t1 + (t339 + t91) * t26 + (t340 - t89) * t25 + (t377 * t75 + t403 * t76) * g(3) + (-t403 * (-t113 * t168 + t171 * t72) + t377 * (t113 * t171 + t168 * t72)) * g(2) + (-t403 * (-t115 * t168 - t171 * t70) + t377 * (t115 * t171 - t168 * t70)) * g(1) - t13 * t88 - t64 * (mrSges(6,1) * t135 + mrSges(6,2) * t134) - t43 * (mrSges(7,1) * t135 + mrSges(7,2) * t134) + (t134 * t402 - t135 * t401) * t426 + (t134 * t424 - t413) * t427 + (-t135 * t422 + t398 + t415) * t428 + t381 + (t135 * t374 + t34) * pkin(5) + t409; -t134 * t88 + t135 * t90 + (-g(1) * t69 + g(2) * t71 - g(3) * t116 + t10 * t135 - t14 * t134 + t5) * m(7) + t15;];
tau  = t6;
