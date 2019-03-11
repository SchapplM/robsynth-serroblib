% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:56
% EndTime: 2019-03-09 18:11:16
% DurationCPUTime: 9.83s
% Computational Cost: add. (13133->566), mult. (31761->734), div. (0->0), fcn. (22408->8), ass. (0->264)
t276 = cos(qJ(6));
t392 = -t276 / 0.2e1;
t440 = Ifges(4,1) + Ifges(5,1);
t439 = Ifges(5,4) + Ifges(4,5);
t274 = sin(qJ(3));
t278 = cos(qJ(2));
t391 = cos(qJ(3));
t314 = t391 * t278;
t275 = sin(qJ(2));
t332 = qJD(1) * t275;
t211 = -qJD(1) * t314 + t274 * t332;
t227 = t274 * t278 + t275 * t391;
t212 = t227 * qJD(1);
t331 = qJD(1) * t278;
t251 = -qJD(1) * pkin(1) - pkin(2) * t331;
t131 = t211 * pkin(3) - t212 * qJ(4) + t251;
t104 = -pkin(4) * t211 - t131;
t273 = sin(qJ(5));
t277 = cos(qJ(5));
t150 = t211 * t277 - t273 * t212;
t270 = qJD(2) + qJD(3);
t267 = qJD(5) - t270;
t272 = sin(qJ(6));
t302 = mrSges(7,1) * t272 + mrSges(7,2) * t276;
t206 = t212 * pkin(9);
t279 = -pkin(3) - pkin(4);
t405 = -pkin(8) - pkin(7);
t252 = t405 * t275;
t237 = qJD(1) * t252;
t221 = qJD(2) * pkin(2) + t237;
t253 = t405 * t278;
t238 = qJD(1) * t253;
t335 = t274 * t238;
t167 = t391 * t221 + t335;
t313 = qJD(4) - t167;
t103 = t270 * t279 - t206 + t313;
t266 = t270 * qJ(4);
t315 = t391 * t238;
t168 = t274 * t221 - t315;
t386 = pkin(9) * t211;
t305 = t168 + t386;
t113 = t266 + t305;
t61 = t103 * t277 - t113 * t273;
t55 = -pkin(5) * t267 - t61;
t424 = t55 * t302;
t152 = t211 * t273 + t277 * t212;
t125 = -t152 * t272 + t267 * t276;
t146 = qJD(6) - t150;
t126 = t152 * t276 + t267 * t272;
t373 = Ifges(7,4) * t126;
t49 = Ifges(7,2) * t125 + Ifges(7,6) * t146 + t373;
t120 = Ifges(7,4) * t125;
t50 = t126 * Ifges(7,1) + t146 * Ifges(7,5) + t120;
t434 = t50 * t392 + t272 * t49 / 0.2e1;
t453 = t152 / 0.2e1;
t442 = Ifges(6,1) * t453;
t454 = t104 * mrSges(6,2) + Ifges(6,4) * t150 + Ifges(6,5) * t267 + t424 - t434 + t442;
t444 = -Ifges(6,2) * t150 / 0.2e1;
t451 = -mrSges(5,1) - mrSges(4,1);
t425 = mrSges(6,3) * t152;
t348 = mrSges(6,1) * t267 + mrSges(7,1) * t125 - mrSges(7,2) * t126 - t425;
t450 = -m(7) * t55 + t348;
t52 = -pkin(5) * t150 - pkin(10) * t152 + t104;
t62 = t103 * t273 + t113 * t277;
t56 = pkin(10) * t267 + t62;
t21 = -t272 * t56 + t276 * t52;
t22 = t272 * t52 + t276 * t56;
t401 = -t146 / 0.2e1;
t403 = -t126 / 0.2e1;
t404 = -t125 / 0.2e1;
t282 = -t104 * mrSges(6,1) - t21 * mrSges(7,1) + t22 * mrSges(7,2) + 0.2e1 * Ifges(6,4) * t453 + 0.2e1 * Ifges(7,5) * t403 + Ifges(6,6) * t267 + 0.2e1 * Ifges(7,6) * t404 + 0.2e1 * Ifges(7,3) * t401 - t444;
t449 = t444 - t282;
t316 = qJD(2) * t405;
t307 = qJD(1) * t316;
t287 = t391 * t307;
t292 = t274 * t307;
t100 = qJD(3) * t168 + t275 * t292 - t278 * t287;
t288 = -t274 * t275 + t314;
t176 = t270 * t288;
t159 = t176 * qJD(1);
t283 = -t159 * pkin(9) + t100;
t177 = t270 * t227;
t160 = t177 * qJD(1);
t312 = qJD(3) * t391;
t329 = qJD(3) * t274;
t99 = t221 * t312 + t238 * t329 + t275 * t287 + t278 * t292;
t97 = t270 * qJD(4) + t99;
t70 = pkin(9) * t160 + t97;
t13 = qJD(5) * t61 + t273 * t283 + t277 * t70;
t247 = Ifges(7,5) * t272 + Ifges(7,6) * t276;
t372 = Ifges(7,4) * t272;
t248 = Ifges(7,2) * t276 + t372;
t371 = Ifges(7,4) * t276;
t249 = Ifges(7,1) * t272 + t371;
t326 = qJD(1) * qJD(2);
t311 = t275 * t326;
t79 = pkin(2) * t311 + t160 * pkin(3) - t159 * qJ(4) - t212 * qJD(4);
t53 = -pkin(4) * t160 - t79;
t66 = qJD(5) * t150 + t159 * t277 + t160 * t273;
t67 = qJD(5) * t152 + t159 * t273 - t277 * t160;
t10 = pkin(5) * t67 - pkin(10) * t66 + t53;
t3 = -qJD(6) * t22 + t10 * t276 - t13 * t272;
t377 = mrSges(7,3) * t272;
t416 = qJD(6) * t126;
t417 = qJD(6) * t125;
t42 = t276 * t66 + t417;
t299 = Ifges(7,5) * t276 - Ifges(7,6) * t272;
t421 = t146 * t299;
t300 = -Ifges(7,2) * t272 + t371;
t426 = t300 / 0.2e1;
t43 = -t272 * t66 - t416;
t327 = qJD(6) * t276;
t328 = qJD(6) * t272;
t433 = -t21 * t327 - t22 * t328;
t441 = qJD(6) / 0.2e1;
t301 = Ifges(7,1) * t276 - t372;
t445 = t301 / 0.2e1;
t448 = -(Ifges(6,6) - t247 / 0.2e1) * t67 - t13 * mrSges(6,2) + Ifges(6,5) * t66 - t3 * t377 + t433 * mrSges(7,3) + qJD(6) * t424 + t421 * t441 + t416 * t445 + t417 * t426 + t43 * t248 / 0.2e1 + t42 * t249 / 0.2e1;
t380 = -Ifges(4,4) + Ifges(5,5);
t438 = -Ifges(4,6) + Ifges(5,6);
t303 = mrSges(7,1) * t276 - mrSges(7,2) * t272;
t437 = mrSges(6,1) + t303;
t204 = Ifges(4,4) * t211;
t361 = t211 * Ifges(5,5);
t436 = t440 * t212 + t439 * t270 - t204 + t361;
t17 = mrSges(7,1) * t67 - mrSges(7,3) * t42;
t18 = -mrSges(7,2) * t67 + mrSges(7,3) * t43;
t298 = -t272 * t17 + t276 * t18;
t80 = -mrSges(7,2) * t146 + mrSges(7,3) * t125;
t81 = mrSges(7,1) * t146 - mrSges(7,3) * t126;
t435 = -t81 * t327 - t80 * t328 + t298;
t174 = t237 * t274 - t315;
t323 = pkin(2) * t329;
t432 = t323 - t174;
t96 = pkin(5) * t152 - pkin(10) * t150;
t124 = t206 + t167;
t242 = t277 * qJ(4) + t273 * t279;
t423 = -qJD(5) * t242 - t277 * t305 + (-qJD(4) + t124) * t273;
t378 = mrSges(4,3) * t211;
t182 = -mrSges(4,2) * t270 - t378;
t185 = -t211 * mrSges(5,2) + mrSges(5,3) * t270;
t334 = t182 + t185;
t360 = t212 * mrSges(4,3);
t379 = mrSges(5,2) * t212;
t420 = -t270 * t451 - t360 - t379;
t181 = t274 * t252 - t391 * t253;
t325 = t391 * pkin(2);
t260 = -t325 - pkin(3);
t257 = -pkin(4) + t260;
t390 = pkin(2) * t274;
t258 = qJ(4) + t390;
t189 = t273 * t257 + t277 * t258;
t418 = -t21 * t272 + t22 * t276;
t345 = Ifges(3,6) * qJD(2);
t375 = Ifges(3,4) * t275;
t415 = t345 / 0.2e1 + (t278 * Ifges(3,2) + t375) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t331);
t2 = qJD(6) * t21 + t10 * t272 + t13 * t276;
t414 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t42 + Ifges(7,6) * t43;
t129 = -mrSges(6,2) * t267 + mrSges(6,3) * t150;
t289 = -t272 * t81 + t276 * t80 + t129;
t297 = t21 * t276 + t22 * t272;
t402 = t126 / 0.2e1;
t412 = t125 * t426 + t301 * t402 + t421 / 0.2e1 - t297 * mrSges(7,3) + t454;
t410 = -0.2e1 * pkin(1);
t9 = t42 * Ifges(7,1) + t43 * Ifges(7,4) + t67 * Ifges(7,5);
t408 = t9 / 0.2e1;
t399 = -t211 / 0.2e1;
t398 = t211 / 0.2e1;
t396 = t212 / 0.2e1;
t394 = t270 / 0.2e1;
t393 = -t272 / 0.2e1;
t389 = pkin(4) * t212;
t388 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t332);
t14 = qJD(5) * t62 + t273 * t70 - t277 * t283;
t180 = -t391 * t252 - t253 * t274;
t144 = -pkin(9) * t227 + t180;
t145 = -pkin(9) * t288 + t181;
t294 = t277 * t144 - t145 * t273;
t385 = t14 * t294;
t384 = t61 * mrSges(6,3);
t383 = t66 * mrSges(6,3);
t382 = t67 * mrSges(6,3);
t381 = mrSges(5,2) + mrSges(4,3);
t261 = -t278 * pkin(2) - pkin(1);
t364 = t152 * t62;
t363 = t159 * mrSges(5,2);
t359 = t212 * Ifges(4,4);
t354 = t276 * mrSges(7,3);
t346 = Ifges(3,5) * qJD(2);
t344 = t100 * t180;
t343 = t100 * t227;
t135 = -pkin(3) * t270 + t313;
t340 = t135 * t211;
t336 = t270 * t277;
t161 = t212 * pkin(3) + t211 * qJ(4);
t175 = t391 * t237 + t335;
t330 = qJD(2) * t275;
t324 = pkin(2) * t332;
t310 = t278 * t326;
t166 = -pkin(3) * t288 - t227 * qJ(4) + t261;
t308 = pkin(2) * t312;
t121 = -t161 - t389;
t304 = t2 * t276 - t272 * t3;
t127 = pkin(4) * t288 - t166;
t171 = t227 * t277 - t273 * t288;
t293 = -t227 * t273 - t277 * t288;
t68 = -pkin(5) * t293 - pkin(10) * t171 + t127;
t88 = t144 * t273 + t145 * t277;
t31 = -t272 * t88 + t276 * t68;
t32 = t272 * t68 + t276 * t88;
t241 = -t273 * qJ(4) + t277 * t279;
t188 = t257 * t277 - t258 * t273;
t136 = t324 + t161;
t290 = t174 + t386;
t90 = pkin(2) * t330 + t177 * pkin(3) - t176 * qJ(4) - t227 * qJD(4);
t239 = t275 * t316;
t240 = t278 * t316;
t105 = t391 * t239 + t274 * t240 + t252 * t312 + t253 * t329;
t69 = -pkin(4) * t177 - t90;
t57 = t121 - t96;
t286 = -qJD(6) * t297 + t304;
t285 = (-t272 * t80 - t276 * t81) * qJD(6) + t298;
t106 = qJD(3) * t181 + t274 * t239 - t391 * t240;
t284 = -t176 * pkin(9) + t106;
t203 = Ifges(5,5) * t212;
t139 = t270 * Ifges(5,6) + t211 * Ifges(5,3) + t203;
t140 = -t211 * Ifges(4,2) + t270 * Ifges(4,6) + t359;
t147 = t266 + t168;
t8 = t42 * Ifges(7,4) + t43 * Ifges(7,2) + t67 * Ifges(7,6);
t280 = -t448 + t97 * mrSges(5,3) - t99 * mrSges(4,2) + t449 * t152 + t451 * t100 - t2 * t354 + t9 * t393 + t140 * t396 + (Ifges(5,3) * t212 - t361) * t399 + t8 * t392 - t167 * t378 + t147 * t379 + t434 * qJD(6) + (-Ifges(4,2) * t212 - t204 + t436) * t398 + t437 * t14 + t438 * t160 - t131 * (mrSges(5,1) * t212 + mrSges(5,3) * t211) - t251 * (mrSges(4,1) * t212 - mrSges(4,2) * t211) + t439 * t159 - (-t211 * t439 + t212 * t438) * t270 / 0.2e1 - (-t440 * t211 + t139 + t203 - t359) * t212 / 0.2e1 + (-t21 * t354 - t22 * t377 - t299 * t401 - t300 * t404 - t301 * t403 - t384 + t442 + t454) * t150;
t262 = Ifges(3,4) * t331;
t255 = t308 + qJD(4);
t236 = -pkin(10) + t242;
t235 = pkin(5) - t241;
t210 = Ifges(3,1) * t332 + t262 + t346;
t193 = t277 * qJD(4) + qJD(5) * t241;
t186 = pkin(5) - t188;
t179 = t212 * t272 + t276 * t336;
t178 = t212 * t276 - t272 * t336;
t165 = mrSges(4,1) * t211 + mrSges(4,2) * t212;
t164 = mrSges(5,1) * t211 - mrSges(5,3) * t212;
t128 = t206 + t175;
t111 = -t136 - t389;
t94 = -mrSges(6,1) * t150 + mrSges(6,2) * t152;
t83 = pkin(9) * t177 + t105;
t78 = qJD(5) * t171 + t176 * t273 - t277 * t177;
t77 = qJD(5) * t293 + t176 * t277 + t177 * t273;
t76 = t277 * t128 + t273 * t290;
t73 = t277 * t124 + t273 * t305;
t65 = Ifges(7,3) * t67;
t54 = t57 - t324;
t30 = t272 * t96 + t276 * t61;
t29 = -t272 * t61 + t276 * t96;
t28 = t272 * t54 + t276 * t76;
t27 = -t272 * t76 + t276 * t54;
t26 = t272 * t57 + t276 * t73;
t25 = -t272 * t73 + t276 * t57;
t24 = qJD(5) * t88 + t273 * t83 - t277 * t284;
t23 = qJD(5) * t294 + t273 * t284 + t277 * t83;
t16 = pkin(5) * t78 - pkin(10) * t77 + t69;
t15 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t5 = -qJD(6) * t32 + t16 * t276 - t23 * t272;
t4 = qJD(6) * t31 + t16 * t272 + t23 * t276;
t1 = [t449 * t78 + t69 * t94 + (t436 / 0.2e1 + mrSges(4,2) * t251 + mrSges(5,2) * t135 - mrSges(4,3) * t167 - mrSges(5,3) * t131 + Ifges(4,4) * t399 + Ifges(5,5) * t398 + t394 * t439 + t396 * t440) * t176 - t294 * t15 + (-t294 * t66 - t61 * t77 - t62 * t78 - t67 * t88) * mrSges(6,3) - (-t13 * mrSges(6,3) + t65 / 0.2e1 - Ifges(6,4) * t66 + t53 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t67 + t414) * t293 + (t288 * t99 + t343) * mrSges(4,3) + (-t345 / 0.2e1 + (mrSges(3,1) * t410 - 0.3e1 / 0.2e1 * t375 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t278) * qJD(1) + (m(4) * (qJD(1) * t261 + t251) + qJD(1) * (-mrSges(4,1) * t288 + mrSges(4,2) * t227) + t165) * pkin(2) - t415) * t330 + (t288 * t97 + t343) * mrSges(5,2) + (mrSges(4,1) * t261 + mrSges(5,1) * t166 + t380 * t227 - (Ifges(4,2) + Ifges(5,3)) * t288 - t381 * t181) * t160 + t79 * (-mrSges(5,1) * t288 - mrSges(5,3) * t227) + t4 * t80 + t5 * t81 + t31 * t17 + t32 * t18 - t348 * t24 + m(5) * (t105 * t147 + t106 * t135 + t131 * t90 + t166 * t79 + t181 * t97 + t344) + m(4) * (t105 * t168 - t106 * t167 + t181 * t99 + t344) + t334 * t105 - t420 * t106 + (t276 * t408 + Ifges(6,1) * t66 + t53 * mrSges(6,2) + t8 * t393 + t42 * t445 + t43 * t426 + (mrSges(6,3) + t302) * t14 + (-t2 * t272 - t276 * t3) * mrSges(7,3) + (-mrSges(7,3) * t418 + t247 * t401 + t248 * t404 + t249 * t403 + t303 * t55 + t392 * t49 + t393 * t50) * qJD(6) + (-Ifges(6,4) + t299 / 0.2e1) * t67) * t171 + t127 * (mrSges(6,1) * t67 + mrSges(6,2) * t66) + t23 * t129 + t90 * t164 + (-t147 * mrSges(5,2) - t168 * mrSges(4,3) + Ifges(5,3) * t398 - Ifges(4,2) * t399 + t139 / 0.2e1 + t131 * mrSges(5,1) - t140 / 0.2e1 + t251 * mrSges(4,1) + t380 * t396 + t438 * t394) * t177 + (mrSges(4,2) * t261 - mrSges(5,3) * t166 + t381 * t180 + t227 * t440 - t380 * t288) * t159 + (t442 + t412) * t77 + m(7) * (t2 * t32 + t21 * t5 + t22 * t4 + t24 * t55 + t3 * t31 - t385) + m(6) * (t104 * t69 + t127 * t53 + t13 * t88 + t23 * t62 - t24 * t61 - t385) + (t210 / 0.2e1 - t388 + t346 / 0.2e1 + (mrSges(3,2) * t410 + 0.3e1 / 0.2e1 * Ifges(3,4) * t278) * qJD(1)) * t278 * qJD(2); (-t159 * t325 - t160 * t390) * mrSges(4,3) + Ifges(3,5) * t310 + t182 * t308 + t280 - t334 * t175 + (-t160 * t258 + t340) * mrSges(5,2) - t111 * t94 + t415 * t332 + (-mrSges(3,1) * t310 + mrSges(3,2) * t311) * pkin(7) - t28 * t80 - t27 * t81 + (-t104 * t111 + t13 * t189 - t14 * t188 - t62 * t76) * m(6) + (t14 * t186 - t21 * t27 - t22 * t28) * m(7) + (t167 * t174 - t168 * t175 - t251 * t324 + (-t391 * t100 + t274 * t99 + (-t167 * t274 + t168 * t391) * qJD(3)) * pkin(2)) * m(4) - t188 * t383 - t189 * t382 - (Ifges(3,5) * t278 - Ifges(3,6) * t275) * t326 / 0.2e1 - t165 * t324 + (m(6) * t62 + m(7) * t418 + t289) * (qJD(5) * t188 + t255 * t277 + t273 * t323) - Ifges(3,6) * t311 + t331 * t388 + (t100 * t260 - t131 * t136 + t258 * t97 + (-t175 + t255) * t147 + t432 * t135) * m(5) - t432 * t420 + (m(7) * t286 + t435) * (-pkin(10) + t189) - t76 * t129 + (m(6) * t61 + t450) * (-qJD(5) * t189 + (-t290 + t323) * t277 + (t128 - t255) * t273) - t136 * t164 - (-Ifges(3,2) * t332 + t210 + t262) * t331 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t275 + mrSges(3,2) * t278) - t275 * (Ifges(3,1) * t278 - t375) / 0.2e1) * qJD(1) ^ 2 + t168 * t360 + t260 * t363 - mrSges(6,3) * t364 + t186 * t15 + t255 * t185; (-t241 * t66 - t242 * t67 - t364) * mrSges(6,3) + t280 + t289 * t193 + (-pkin(3) * t159 - qJ(4) * t160 + t340) * mrSges(5,2) + t285 * t236 - t26 * t80 - t25 * t81 + (t420 + t360) * t168 - t334 * t167 - t121 * t94 - t73 * t129 - t161 * t164 + qJD(4) * t185 + t235 * t15 + t348 * t423 + (-t104 * t121 + t13 * t242 - t14 * t241 + (t193 - t73) * t62 + t423 * t61) * m(6) + (-pkin(3) * t100 + qJ(4) * t97 - t131 * t161 - t135 * t168 + t147 * t313) * m(5) + (t14 * t235 + t193 * t418 - t21 * t25 - t22 * t26 + t286 * t236 - t423 * t55) * m(7); t363 - t178 * t81 - t179 * t80 - t270 * t185 + (t164 - t94) * t212 + (qJD(5) * t289 - t270 * t129 - t15 - t383) * t277 + (-t267 * t348 + t285 - t382) * t273 + ((qJD(5) * t418 - t14) * t277 - t178 * t21 - t179 * t22 + (t267 * t55 + t286) * t273) * m(7) + (-t104 * t212 + t13 * t273 - t14 * t277 + t267 * (-t273 * t61 + t277 * t62)) * m(6) + (t131 * t212 - t147 * t270 + t100) * m(5); -t30 * t80 - t29 * t81 + t282 * t152 - m(7) * (t21 * t29 + t22 * t30) + (t425 + t450) * t62 - pkin(5) * t15 + (-qJD(6) * t49 / 0.2e1 + t408) * t272 + (t384 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t152 - t412) * t150 + (m(7) * (t304 + t433) + t435) * pkin(10) + (t50 * t441 + t2 * mrSges(7,3) + t8 / 0.2e1) * t276 - t61 * t129 + (-m(7) * pkin(5) - t437) * t14 + t448; t65 - t55 * (mrSges(7,1) * t126 + mrSges(7,2) * t125) + (Ifges(7,1) * t125 - t373) * t403 + t49 * t402 + (Ifges(7,5) * t125 - Ifges(7,6) * t126) * t401 - t21 * t80 + t22 * t81 + (t125 * t21 + t126 * t22) * mrSges(7,3) + (-Ifges(7,2) * t126 + t120 + t50) * t404 + t414;];
tauc  = t1(:);
