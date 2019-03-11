% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:35
% EndTime: 2019-03-09 03:28:08
% DurationCPUTime: 23.06s
% Computational Cost: add. (6482->636), mult. (12965->818), div. (0->0), fcn. (8427->10), ass. (0->287)
t223 = sin(qJ(3));
t306 = qJD(1) * t223;
t207 = qJD(5) + t306;
t353 = -t207 / 0.2e1;
t219 = sin(pkin(9));
t220 = cos(pkin(9));
t226 = cos(qJ(3));
t305 = qJD(1) * t226;
t168 = qJD(3) * t219 + t220 * t305;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t298 = t220 * qJD(3);
t243 = t219 * t305 - t298;
t233 = t225 * t168 - t222 * t243;
t363 = -t233 / 0.2e1;
t398 = Ifges(7,4) + Ifges(6,5);
t420 = Ifges(7,2) + Ifges(6,3);
t434 = t420 * t353 + t398 * t363;
t104 = t222 * t168 + t225 * t243;
t364 = t104 / 0.2e1;
t433 = t223 / 0.2e1;
t432 = t226 / 0.2e1;
t422 = qJD(3) / 0.2e1;
t402 = mrSges(6,1) + mrSges(7,1);
t401 = mrSges(6,2) - mrSges(7,3);
t400 = Ifges(6,1) + Ifges(7,1);
t399 = Ifges(6,4) - Ifges(7,5);
t397 = Ifges(6,6) - Ifges(7,6);
t352 = t207 / 0.2e1;
t362 = t233 / 0.2e1;
t365 = -t104 / 0.2e1;
t228 = -pkin(1) - pkin(7);
t201 = qJD(1) * t228 + qJD(2);
t312 = t226 * t201;
t150 = -qJD(3) * pkin(3) + qJD(4) - t312;
t111 = pkin(4) * t243 + t150;
t184 = pkin(3) * t223 - qJ(4) * t226 + qJ(2);
t158 = t184 * qJD(1);
t183 = t223 * t201;
t161 = qJD(3) * qJ(4) + t183;
t90 = t220 * t158 - t161 * t219;
t62 = pkin(4) * t306 - pkin(8) * t168 + t90;
t91 = t219 * t158 + t220 * t161;
t64 = -pkin(8) * t243 + t91;
t19 = -t222 * t64 + t225 * t62;
t16 = -pkin(5) * t207 + qJD(6) - t19;
t29 = t104 * pkin(5) - qJ(6) * t233 + t111;
t101 = Ifges(6,4) * t104;
t332 = Ifges(7,5) * t104;
t394 = t207 * t398 + t233 * t400 - t101 + t332;
t413 = -mrSges(6,2) * t111 - mrSges(7,2) * t16 + mrSges(6,3) * t19 + mrSges(7,3) * t29 - t394 / 0.2e1;
t431 = -Ifges(6,4) * t365 - Ifges(7,5) * t364 - t398 * t352 - t400 * t362 + t413;
t424 = -m(6) - m(7);
t258 = Ifges(5,5) * t220 - Ifges(5,6) * t219;
t269 = mrSges(4,1) * t226 - mrSges(4,2) * t223;
t336 = Ifges(4,4) * t226;
t430 = qJ(2) * t269 + (-Ifges(4,1) * t223 - t336) * t432 + (Ifges(5,3) * t226 - t223 * t258) * t433;
t295 = qJD(1) * qJD(3);
t182 = qJDD(1) * t223 + t226 * t295;
t170 = qJDD(5) + t182;
t181 = qJDD(1) * t226 - t223 * t295;
t132 = qJDD(3) * t219 + t181 * t220;
t199 = qJDD(1) * t228 + qJDD(2);
t303 = qJD(3) * t226;
t120 = t223 * t199 + t201 * t303;
t110 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t120;
t270 = -qJD(4) * t226 + qJD(2);
t294 = qJDD(1) * qJ(2);
t87 = pkin(3) * t182 - qJ(4) * t181 + qJD(1) * t270 + t294;
t49 = -t110 * t219 + t220 * t87;
t22 = pkin(4) * t182 - pkin(8) * t132 + t49;
t131 = qJDD(3) * t220 - t181 * t219;
t50 = t220 * t110 + t219 * t87;
t24 = pkin(8) * t131 + t50;
t300 = qJD(5) * t225;
t301 = qJD(5) * t222;
t3 = t222 * t22 + t225 * t24 + t62 * t300 - t301 * t64;
t1 = qJ(6) * t170 + qJD(6) * t207 + t3;
t356 = t170 / 0.2e1;
t42 = qJD(5) * t233 - t225 * t131 + t222 * t132;
t371 = t42 / 0.2e1;
t372 = -t42 / 0.2e1;
t41 = -qJD(5) * t104 + t222 * t131 + t225 * t132;
t373 = t41 / 0.2e1;
t304 = qJD(3) * t223;
t119 = t199 * t226 - t201 * t304;
t116 = -qJDD(3) * pkin(3) + qJDD(4) - t119;
t65 = -pkin(4) * t131 + t116;
t5 = pkin(5) * t42 - qJ(6) * t41 - qJD(6) * t233 + t65;
t429 = mrSges(6,1) * t65 + mrSges(7,1) * t5 - mrSges(7,2) * t1 - mrSges(6,3) * t3 - t41 * Ifges(6,4) / 0.2e1 - t170 * Ifges(6,6) / 0.2e1 + 0.2e1 * Ifges(7,3) * t371 + (-t372 + t371) * Ifges(6,2) + (-t399 + Ifges(7,5)) * t373 + (-t397 + Ifges(7,6)) * t356;
t20 = t222 * t62 + t225 * t64;
t17 = qJ(6) * t207 + t20;
t100 = Ifges(7,5) * t233;
t43 = Ifges(7,6) * t207 + Ifges(7,3) * t104 + t100;
t333 = Ifges(6,4) * t233;
t46 = -Ifges(6,2) * t104 + Ifges(6,6) * t207 + t333;
t377 = mrSges(6,1) * t111 + mrSges(7,1) * t29 - mrSges(7,2) * t17 - mrSges(6,3) * t20 - t46 / 0.2e1 + t43 / 0.2e1;
t428 = -Ifges(6,2) * t365 + Ifges(7,3) * t364 - t352 * t397 - t362 * t399 + t377;
t261 = -Ifges(4,2) * t223 + t336;
t426 = t16 * mrSges(7,1) + t20 * mrSges(6,2) + Ifges(4,6) * t422 + qJD(1) * t261 / 0.2e1 - Ifges(5,5) * t168 / 0.2e1 + Ifges(5,6) * t243 / 0.2e1 - Ifges(5,3) * t306 / 0.2e1 + t397 * t364 - t17 * mrSges(7,3) - t19 * mrSges(6,1) + t434;
t425 = -m(5) - m(4);
t423 = t170 * t398 - t399 * t42 + t400 * t41;
t421 = mrSges(6,3) + mrSges(7,2);
t418 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t243 + t168 * mrSges(5,2) + mrSges(4,3) * t305;
t334 = Ifges(5,4) * t220;
t260 = -Ifges(5,2) * t219 + t334;
t335 = Ifges(5,4) * t219;
t262 = Ifges(5,1) * t220 - t335;
t417 = t168 * (Ifges(5,5) * t226 - t223 * t262) - (Ifges(5,6) * t226 - t223 * t260) * t243;
t218 = pkin(9) + qJ(5);
t211 = sin(t218);
t212 = cos(t218);
t415 = -t211 * t401 + t212 * t402;
t414 = t170 * t420 - t397 * t42 + t398 * t41;
t296 = qJD(1) * qJD(2);
t202 = t294 + t296;
t250 = t119 * t226 + t120 * t223;
t266 = t219 * mrSges(5,1) + t220 * mrSges(5,2);
t319 = t220 * t223;
t321 = t150 * t223;
t411 = t266 * t321 - t91 * (mrSges(5,3) * t219 * t223 - mrSges(5,2) * t226) - t90 * (mrSges(5,1) * t226 + mrSges(5,3) * t319);
t410 = mrSges(3,2) - t266 - mrSges(4,3) - mrSges(2,1);
t267 = -mrSges(5,1) * t220 + mrSges(5,2) * t219;
t240 = m(5) * pkin(3) - t267;
t268 = mrSges(4,1) * t223 + mrSges(4,2) * t226;
t282 = m(5) * qJ(4) + mrSges(5,3);
t409 = -t223 * t240 + t226 * t282 + mrSges(2,2) - mrSges(3,3) - t268;
t408 = -m(7) * pkin(5) - t402;
t407 = -m(7) * qJ(6) + t401;
t354 = t182 / 0.2e1;
t359 = t132 / 0.2e1;
t406 = Ifges(5,1) * t359 + Ifges(5,5) * t354;
t360 = t131 / 0.2e1;
t227 = cos(qJ(1));
t316 = t223 * t227;
t344 = pkin(8) + qJ(4);
t403 = g(2) * t344 * t316;
t348 = g(2) * t227;
t25 = mrSges(6,1) * t170 - mrSges(6,3) * t41;
t26 = -t170 * mrSges(7,1) + t41 * mrSges(7,2);
t396 = t26 - t25;
t27 = -mrSges(6,2) * t170 - mrSges(6,3) * t42;
t28 = -mrSges(7,2) * t42 + mrSges(7,3) * t170;
t395 = t27 + t28;
t172 = t219 * t225 + t220 * t222;
t153 = t172 * qJD(1);
t129 = t223 * t153;
t313 = t225 * t220;
t171 = t219 * t222 - t313;
t390 = t171 * t223;
t130 = qJD(1) * t390;
t287 = t219 * t306;
t133 = -pkin(4) * t287 + t183;
t154 = t171 * qJD(5);
t381 = t172 * qJD(5);
t393 = -qJD(6) * t172 - t133 + (t154 + t130) * qJ(6) + (t381 + t129) * pkin(5);
t339 = mrSges(6,3) * t104;
t73 = -mrSges(6,2) * t207 - t339;
t341 = mrSges(7,2) * t104;
t74 = mrSges(7,3) * t207 - t341;
t343 = t74 + t73;
t338 = mrSges(6,3) * t233;
t75 = mrSges(6,1) * t207 - t338;
t340 = mrSges(7,2) * t233;
t76 = -mrSges(7,1) * t207 + t340;
t342 = t76 - t75;
t144 = t171 * t226;
t392 = -qJD(3) * t144 - t223 * t381 - t153;
t142 = t172 * t226;
t307 = qJD(1) * t219;
t385 = -t219 * t301 + t220 * t300;
t391 = qJD(1) * t313 + qJD(3) * t142 - t222 * t307 + t223 * t385;
t209 = pkin(4) * t220 + pkin(3);
t256 = pkin(5) * t212 + qJ(6) * t211;
t388 = -m(7) * (-t209 - t256) + m(6) * t209 + t415;
t386 = t344 * t424 - t421;
t127 = -mrSges(5,2) * t306 - mrSges(5,3) * t243;
t128 = mrSges(5,1) * t306 - mrSges(5,3) * t168;
t384 = t127 * t220 - t128 * t219;
t255 = -t219 * t49 + t220 * t50;
t93 = -mrSges(5,2) * t182 + mrSges(5,3) * t131;
t94 = mrSges(5,1) * t182 - mrSges(5,3) * t132;
t383 = -t219 * t94 + t220 * t93;
t224 = sin(qJ(1));
t349 = g(1) * t224;
t382 = t348 - t349;
t12 = t42 * mrSges(7,1) - t41 * mrSges(7,3);
t13 = t42 * mrSges(6,1) + t41 * mrSges(6,2);
t70 = -t131 * mrSges(5,1) + t132 * mrSges(5,2);
t380 = -t12 - t13 - t70;
t4 = -qJD(5) * t20 + t22 * t225 - t222 * t24;
t166 = t220 * t184;
t275 = -t219 * t228 + pkin(4);
t318 = t220 * t226;
t102 = -pkin(8) * t318 + t223 * t275 + t166;
t315 = t223 * t228;
t123 = t219 * t184 + t220 * t315;
t320 = t219 * t226;
t113 = -pkin(8) * t320 + t123;
t327 = t222 * t102 + t225 * t113;
t257 = pkin(3) * t226 + qJ(4) * t223;
t145 = qJD(3) * t257 + t270;
t126 = t220 * t145;
t292 = pkin(8) * t319;
t69 = t126 + (t226 * t275 + t292) * qJD(3);
t302 = qJD(3) * t228;
t285 = t226 * t302;
t109 = t219 * t145 + t220 * t285;
t286 = t219 * t304;
t86 = pkin(8) * t286 + t109;
t15 = -qJD(5) * t327 - t222 * t86 + t225 * t69;
t379 = t223 * t282 + t226 * t240;
t376 = qJD(1) ^ 2;
t368 = Ifges(5,4) * t360 + t406;
t366 = -m(3) - m(4);
t350 = pkin(4) * t219;
t347 = g(3) * t226;
t175 = t257 * qJD(1);
t114 = t220 * t175 - t219 * t312;
t77 = (pkin(4) * t226 + t292) * qJD(1) + t114;
t115 = t219 * t175 + t220 * t312;
t89 = pkin(8) * t287 + t115;
t33 = t222 * t77 + t225 * t89;
t337 = Ifges(4,4) * t223;
t326 = t116 * t226;
t317 = t223 * t224;
t314 = t224 * t226;
t311 = t226 * t227;
t310 = t226 * t228;
t309 = t227 * t212;
t308 = t227 * pkin(1) + t224 * qJ(2);
t299 = qJDD(1) * mrSges(3,2);
t297 = m(5) - t424;
t288 = t227 * pkin(7) + t308;
t206 = t223 * t302;
t279 = t306 / 0.2e1;
t274 = -t295 / 0.2e1;
t271 = (t202 + t296) * qJ(2);
t167 = pkin(4) * t320 - t310;
t263 = t226 * Ifges(4,1) - t337;
t259 = -Ifges(4,5) * t223 - Ifges(4,6) * t226;
t254 = -t219 * t90 + t220 * t91;
t32 = -t222 * t89 + t225 * t77;
t53 = t102 * t225 - t113 * t222;
t189 = t344 * t219;
t190 = t344 * t220;
t249 = -t225 * t189 - t190 * t222;
t118 = -t189 * t222 + t190 * t225;
t151 = -pkin(4) * t286 + t206;
t14 = t102 * t300 - t113 * t301 + t222 * t69 + t225 * t86;
t245 = t223 * (-Ifges(4,2) * t226 - t337);
t215 = t227 * qJ(2);
t210 = -qJDD(1) * pkin(1) + qJDD(2);
t192 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t306;
t178 = t209 * t314;
t176 = t268 * qJD(1);
t160 = Ifges(4,5) * qJD(3) + qJD(1) * t263;
t147 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t182;
t146 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t181;
t141 = t172 * t223;
t137 = -t211 * t224 + t223 * t309;
t136 = t211 * t316 + t212 * t224;
t135 = t211 * t227 + t212 * t317;
t134 = t211 * t317 - t309;
t122 = -t219 * t315 + t166;
t108 = -t219 * t285 + t126;
t98 = pkin(5) * t171 - qJ(6) * t172 - t209;
t97 = Ifges(5,1) * t168 - Ifges(5,4) * t243 + Ifges(5,5) * t306;
t96 = Ifges(5,4) * t168 - Ifges(5,2) * t243 + Ifges(5,6) * t306;
t85 = -t222 * t223 * t298 - t225 * t286 + t226 * t385;
t83 = qJD(3) * t390 - t226 * t381;
t72 = qJD(4) * t172 + qJD(5) * t118;
t71 = -qJD(4) * t171 + qJD(5) * t249;
t63 = pkin(5) * t142 + qJ(6) * t144 + t167;
t60 = t132 * Ifges(5,4) + t131 * Ifges(5,2) + t182 * Ifges(5,6);
t57 = mrSges(6,1) * t104 + mrSges(6,2) * t233;
t56 = mrSges(7,1) * t104 - mrSges(7,3) * t233;
t55 = pkin(5) * t233 + qJ(6) * t104;
t52 = -pkin(5) * t223 - t53;
t51 = qJ(6) * t223 + t327;
t31 = -pkin(5) * t305 - t32;
t30 = qJ(6) * t305 + t33;
t18 = pkin(5) * t85 - qJ(6) * t83 + qJD(6) * t144 + t151;
t11 = -pkin(5) * t303 - t15;
t10 = qJ(6) * t303 + qJD(6) * t223 + t14;
t2 = -pkin(5) * t170 + qJDD(6) - t4;
t6 = [t428 * t85 + t430 * t295 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (t1 * t223 + t144 * t5) * mrSges(7,3) - t431 * t83 + t429 * t142 + (-t144 * t65 - t223 * t3) * mrSges(6,2) + t417 * t422 - t250 * mrSges(4,3) + (-t70 + t146) * t310 + (Ifges(4,1) * t181 - Ifges(4,4) * t182) * t432 + (Ifges(5,5) * t132 + Ifges(5,6) * t131 + Ifges(5,3) * t182 + t414) * t433 - t223 * (Ifges(4,4) * t181 - Ifges(4,2) * t182) / 0.2e1 + m(5) * (t108 * t90 + t109 * t91 + t122 * t49 + t123 * t50 + (t150 * t304 - t326) * t228) + t50 * (-mrSges(5,2) * t223 - mrSges(5,3) * t320) + t49 * (mrSges(5,1) * t223 - mrSges(5,3) * t318) - t60 * t320 / 0.2e1 + m(4) * (t228 * t250 + t271) + m(3) * (-pkin(1) * t210 + t271) - pkin(1) * t299 + t96 * t286 / 0.2e1 + (t268 + 0.2e1 * mrSges(3,3)) * t202 + (Ifges(6,6) * t365 + Ifges(7,6) * t364 + t420 * t352 + t398 * t362 - t426) * t303 + (-Ifges(7,5) * t144 + Ifges(7,6) * t223) * t371 + (-Ifges(6,4) * t144 + Ifges(6,6) * t223) * t372 - (t220 * t97 + t160) * t304 / 0.2e1 + t266 * t326 + t147 * t315 + t327 * t27 + m(6) * (t111 * t151 + t14 * t20 + t15 * t19 + t167 * t65 + t3 * t327 + t4 * t53) + t4 * (mrSges(6,1) * t223 + mrSges(6,3) * t144) + t2 * (-mrSges(7,1) * t223 - mrSges(7,2) * t144) + qJDD(3) * (Ifges(4,5) * t226 - Ifges(4,6) * t223) + t210 * mrSges(3,2) + qJ(2) * (mrSges(4,1) * t182 + mrSges(4,2) * t181) + qJD(2) * t176 + t167 * t13 + t151 * t57 + t109 * t127 + t108 * t128 + t122 * t94 + t123 * t93 + t10 * t74 + t15 * t75 + t11 * t76 - t182 * t261 / 0.2e1 + t181 * t263 / 0.2e1 + t14 * t73 + t63 * t12 + t51 * t28 + t52 * t26 + t53 * t25 + t18 * t56 + t418 * t206 + t192 * t285 + (-t144 * t398 + t223 * t420) * t356 + (-t144 * t400 + t223 * t398) * t373 + (t259 * t422 - t411) * qJD(3) - t423 * t144 / 0.2e1 + (-m(3) * t308 + t421 * t314 + t425 * t288 + t424 * (t209 * t317 + t227 * t350 - t314 * t344 + t288) + t408 * t135 + t407 * t134 + t410 * t227 + t409 * t224) * g(2) + (t424 * ((-t350 + t228) * t224 + t209 * t316 - t344 * t311 + t215) + t421 * t311 + t408 * t137 + t407 * t136 + (-m(3) + t425) * t215 + t409 * t227 + (m(3) * pkin(1) + t228 * t425 - t410) * t224) * g(1) + t245 * t274 + (Ifges(5,3) * t223 + t226 * t258) * t354 + (Ifges(5,5) * t223 + t226 * t262) * t359 + (Ifges(5,6) * t223 + t226 * t260) * t360 + t318 * t368 + m(7) * (t1 * t51 + t10 * t17 + t11 * t16 + t18 * t29 + t2 * t52 + t5 * t63); t299 - t395 * t390 + t396 * t141 + (qJ(2) * t366 - mrSges(3,3)) * t376 + (t147 + t383) * t223 + (-t127 * t219 - t128 * t220 - t176) * qJD(1) + (t146 + t380) * t226 + ((t192 + t384) * t226 + (t56 + t57 + t418) * t223) * qJD(3) + m(4) * t250 + m(3) * t210 + t392 * t343 + t391 * t342 + t382 * (t297 - t366) + (-t1 * t390 + t141 * t2 + t16 * t391 + t17 * t392 - t226 * t5 + t29 * t304) * m(7) + (t111 * t304 - t141 * t4 - t19 * t391 + t20 * t392 - t226 * t65 - t3 * t390) * m(6) + (-t326 + t255 * t223 + (t226 * t254 + t321) * qJD(3) - (t219 * t91 + t220 * t90) * qJD(1)) * m(5); t428 * t381 + (t245 / 0.2e1 - t430) * t376 + t431 * t154 + t429 * t171 + (t368 + t406) * t219 + (Ifges(6,4) * t364 + Ifges(7,5) * t365 + t398 * t353 + t400 * t363 + t413) * t130 + (-m(7) * t178 + ((-m(7) * t256 - t415) * t226 + t386 * t223) * t224) * g(1) + t382 * t269 + t383 * qJ(4) + t255 * mrSges(5,3) + t384 * qJD(4) + (Ifges(6,6) * t364 + Ifges(7,6) * t365 + t426 + t434) * t305 - t192 * t312 - t96 * t287 / 0.2e1 + t335 * t360 + (t60 / 0.2e1 + t97 * t279 + Ifges(5,6) * t354 + Ifges(5,2) * t360) * t220 + (t411 - t417 / 0.2e1) * qJD(1) - t396 * t249 + (t268 + (-t282 + t386) * t226 + (t240 + t388) * t223) * g(3) + (t403 + t1 * t118 - t249 * t2 + t5 * t98 + t393 * t29 + (t71 - t30) * t17 + (t72 - t31) * t16) * m(7) + (-g(1) * t178 + t403 - t111 * t133 + t249 * t4 + t118 * t3 - t209 * t65 + (t71 - t33) * t20 + (-t72 - t32) * t19) * m(6) - t418 * t183 + t342 * t72 + t343 * t71 + t393 * t56 + t395 * t118 + (-pkin(3) * t116 + qJ(4) * t255 + qJD(4) * t254 - t114 * t90 - t115 * t91 - t150 * t183) * m(5) - t379 * t349 + (mrSges(6,2) * t65 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t5 + Ifges(6,4) * t372 + Ifges(7,5) * t371 + t356 * t398 + t373 * t400 + t423 / 0.2e1) * t172 + (Ifges(6,2) * t364 - Ifges(7,3) * t365 + t353 * t397 + t363 * t399 + t377) * t129 - t209 * t13 + Ifges(4,5) * t181 - Ifges(4,6) * t182 - t133 * t57 - t115 * t127 - t114 * t128 - t120 * mrSges(4,2) + Ifges(4,3) * qJDD(3) + t119 * mrSges(4,1) + t98 * t12 - t30 * t74 - t32 * t75 - t31 * t76 + t116 * t267 - pkin(3) * t70 - t33 * t73 + (t223 * t421 + t388 * t226 + t379) * t348 + t259 * t274 + t160 * t279 + t334 * t359; -t342 * t233 + t343 * t104 + t168 * t128 + (t226 * t307 - t298) * t127 + (-t223 * g(3) - t226 * t382) * t297 + (t104 * t17 - t16 * t233 + t5) * m(7) + (t104 * t20 + t19 * t233 + t65) * m(6) + (t90 * t168 + t243 * t91 + t116) * m(5) - t380; (-(-pkin(5) * t211 + qJ(6) * t212) * t347 - t29 * t55 - g(2) * (pkin(5) * t136 - qJ(6) * t137) - g(1) * (-pkin(5) * t134 + qJ(6) * t135) - pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t17) * m(7) + (-m(7) * t16 + t338 - t342) * t20 + (-m(7) * t17 - t339 - t343) * t19 + (-Ifges(6,2) * t233 - t101 + t394) * t364 - t111 * (mrSges(6,1) * t233 - mrSges(6,2) * t104) - t29 * (mrSges(7,1) * t233 + mrSges(7,3) * t104) + (Ifges(7,3) * t233 - t332) * t365 + (t134 * t402 + t135 * t401) * g(1) + (t211 * t402 + t212 * t401) * t347 + (-t136 * t402 - t137 * t401) * g(2) + (-t104 * t398 - t233 * t397) * t353 + (-t104 * t400 + t100 - t333 + t43) * t363 + t414 + qJD(6) * t74 - t55 * t56 - pkin(5) * t26 + qJ(6) * t28 - t3 * mrSges(6,2) + t4 * mrSges(6,1) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t17 * t340 + t16 * t341 + t46 * t362; t233 * t56 - t207 * t74 + (-g(1) * t134 + g(2) * t136 - t17 * t207 - t211 * t347 + t233 * t29 + t2) * m(7) + t26;];
tau  = t6;
