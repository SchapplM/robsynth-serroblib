% Calculate vector of inverse dynamics joint torques for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:09
% EndTime: 2019-12-31 20:32:48
% DurationCPUTime: 21.03s
% Computational Cost: add. (8704->691), mult. (19780->946), div. (0->0), fcn. (14322->14), ass. (0->312)
t258 = sin(pkin(9));
t259 = cos(pkin(9));
t262 = sin(qJ(4));
t266 = cos(qJ(4));
t213 = t258 * t266 + t259 * t262;
t267 = cos(qJ(2));
t283 = t213 * t267;
t170 = qJD(1) * t283;
t191 = t213 * qJD(4);
t448 = t170 - t191;
t289 = t258 * t262 - t259 * t266;
t282 = t289 * t267;
t171 = qJD(1) * t282;
t190 = t289 * qJD(4);
t447 = t171 - t190;
t451 = mrSges(5,3) + mrSges(6,3);
t263 = sin(qJ(2));
t292 = pkin(2) * t263 - qJ(3) * t267;
t216 = t292 * qJD(1);
t335 = qJD(1) * t263;
t320 = t258 * t335;
t154 = pkin(6) * t320 + t259 * t216;
t341 = t259 * t267;
t288 = pkin(3) * t263 - pkin(7) * t341;
t120 = qJD(1) * t288 + t154;
t192 = t258 * t216;
t342 = t259 * t263;
t343 = t258 * t267;
t284 = -pkin(6) * t342 - pkin(7) * t343;
t141 = qJD(1) * t284 + t192;
t365 = pkin(7) + qJ(3);
t223 = t365 * t258;
t224 = t365 * t259;
t152 = -t262 * t223 + t266 * t224;
t436 = -t213 * qJD(3) - qJD(4) * t152 - t266 * t120 + t141 * t262;
t327 = qJD(4) * t266;
t330 = qJD(3) * t259;
t331 = qJD(3) * t258;
t435 = -t223 * t327 + (-t141 + t330) * t266 + (-qJD(4) * t224 - t120 - t331) * t262;
t450 = -pkin(4) * t335 - t447 * pkin(8) + t436;
t449 = -t448 * pkin(8) - t435;
t241 = pkin(3) * t259 + pkin(2);
t257 = pkin(9) + qJ(4);
t249 = cos(t257);
t215 = pkin(4) * t249 + t241;
t251 = qJ(5) + t257;
t239 = sin(t251);
t240 = cos(t251);
t248 = sin(t257);
t446 = -m(5) * t241 - m(6) * t215 - t249 * mrSges(5,1) - t240 * mrSges(6,1) + t248 * mrSges(5,2) + t239 * mrSges(6,2);
t256 = -pkin(8) - t365;
t445 = -m(5) * t365 + m(6) * t256 - t451;
t334 = qJD(1) * t267;
t238 = qJD(4) - t334;
t206 = qJD(2) * t259 - t320;
t318 = t259 * t335;
t207 = qJD(2) * t258 + t318;
t303 = t266 * t206 - t207 * t262;
t138 = t206 * t262 + t207 * t266;
t358 = Ifges(5,4) * t138;
t70 = Ifges(5,2) * t303 + Ifges(5,6) * t238 + t358;
t398 = t70 / 0.2e1;
t134 = Ifges(5,4) * t303;
t71 = Ifges(5,1) * t138 + Ifges(5,5) * t238 + t134;
t397 = t71 / 0.2e1;
t325 = qJD(1) * qJD(2);
t219 = qJDD(1) * t263 + t267 * t325;
t264 = sin(qJ(1));
t268 = cos(qJ(1));
t444 = g(1) * t268 + g(2) * t264;
t261 = sin(qJ(5));
t265 = cos(qJ(5));
t443 = -t138 * t261 + t265 * t303;
t80 = t138 * t265 + t261 * t303;
t172 = qJDD(2) * t259 - t219 * t258;
t173 = qJDD(2) * t258 + t219 * t259;
t67 = qJD(4) * t303 + t172 * t262 + t173 * t266;
t68 = -qJD(4) * t138 + t172 * t266 - t173 * t262;
t22 = qJD(5) * t443 + t261 * t68 + t265 * t67;
t408 = t22 / 0.2e1;
t23 = -qJD(5) * t80 - t261 * t67 + t265 * t68;
t407 = t23 / 0.2e1;
t400 = t67 / 0.2e1;
t399 = t68 / 0.2e1;
t387 = t172 / 0.2e1;
t386 = t173 / 0.2e1;
t250 = t267 * qJDD(1);
t313 = t263 * t325;
t218 = -t250 + t313;
t211 = qJDD(4) + t218;
t198 = qJDD(5) + t211;
t385 = t198 / 0.2e1;
t384 = t211 / 0.2e1;
t442 = -t218 / 0.2e1;
t383 = t218 / 0.2e1;
t441 = t219 / 0.2e1;
t440 = qJD(2) / 0.2e1;
t151 = -t266 * t223 - t224 * t262;
t112 = -pkin(8) * t213 + t151;
t113 = -pkin(8) * t289 + t152;
t57 = t112 * t265 - t113 * t261;
t439 = qJD(5) * t57 + t261 * t450 - t449 * t265;
t58 = t112 * t261 + t113 * t265;
t438 = -qJD(5) * t58 + t449 * t261 + t265 * t450;
t409 = m(6) * pkin(4);
t437 = mrSges(5,1) + t409;
t300 = -mrSges(4,1) * t259 + mrSges(4,2) * t258;
t281 = m(4) * pkin(2) - t300;
t433 = t267 * t281;
t307 = -qJ(3) * t263 - pkin(1);
t222 = -pkin(2) * t267 + t307;
t203 = t259 * t222;
t142 = -pkin(7) * t342 + t203 + (-pkin(6) * t258 - pkin(3)) * t267;
t167 = pkin(6) * t341 + t258 * t222;
t344 = t258 * t263;
t150 = -pkin(7) * t344 + t167;
t85 = t262 * t142 + t266 * t150;
t430 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t206 + mrSges(4,2) * t207 + mrSges(3,3) * t335;
t359 = Ifges(4,4) * t259;
t295 = -Ifges(4,2) * t258 + t359;
t360 = Ifges(4,4) * t258;
t297 = Ifges(4,1) * t259 - t360;
t429 = t206 * (Ifges(4,6) * t263 + t267 * t295) + t207 * (Ifges(4,5) * t263 + t267 * t297);
t243 = pkin(6) * t250;
t204 = -pkin(6) * t313 + t243;
t205 = t219 * pkin(6);
t428 = t204 * t267 + t205 * t263;
t329 = qJD(3) * t263;
t350 = qJDD(1) * pkin(1);
t128 = pkin(2) * t218 - qJ(3) * t219 - qJD(1) * t329 - t350;
t245 = pkin(6) * t335;
t177 = qJDD(2) * qJ(3) + t243 + (qJD(3) - t245) * qJD(2);
t86 = t259 * t128 - t177 * t258;
t87 = t258 * t128 + t259 * t177;
t427 = -t258 * t86 + t259 * t87;
t246 = pkin(6) * t334;
t319 = t258 * t334;
t199 = pkin(3) * t319 + t246;
t426 = -pkin(4) * t448 - t199;
t425 = m(4) + m(5) + m(6);
t233 = qJD(5) + t238;
t424 = t207 * Ifges(4,5) + t138 * Ifges(5,5) + t80 * Ifges(6,5) + t206 * Ifges(4,6) + Ifges(5,6) * t303 + Ifges(6,6) * t443 - Ifges(4,3) * t334 + t238 * Ifges(5,3) + t233 * Ifges(6,3);
t423 = -m(3) - t425;
t196 = t222 * qJD(1);
t226 = qJD(2) * qJ(3) + t246;
t143 = t259 * t196 - t226 * t258;
t144 = t258 * t196 + t259 * t226;
t221 = -qJD(2) * pkin(2) + qJD(3) + t245;
t299 = mrSges(4,1) * t258 + t259 * mrSges(4,2);
t422 = -t221 * t267 * t299 - t144 * (-mrSges(4,2) * t263 - mrSges(4,3) * t343) - t143 * (mrSges(4,1) * t263 - mrSges(4,3) * t341);
t96 = -pkin(3) * t334 - pkin(7) * t207 + t143;
t99 = pkin(7) * t206 + t144;
t49 = t262 * t96 + t266 * t99;
t43 = pkin(8) * t303 + t49;
t351 = t265 * t43;
t48 = -t262 * t99 + t266 * t96;
t42 = -pkin(8) * t138 + t48;
t41 = pkin(4) * t238 + t42;
t14 = t261 * t41 + t351;
t153 = -pkin(3) * t206 + t221;
t90 = -pkin(4) * t303 + t153;
t421 = -t90 * mrSges(6,1) + t14 * mrSges(6,3);
t354 = t261 * t43;
t13 = t265 * t41 - t354;
t420 = t90 * mrSges(6,2) - mrSges(6,3) * t13;
t244 = Ifges(3,4) * t334;
t419 = t259 * (t207 * Ifges(4,1) + t206 * Ifges(4,4) - Ifges(4,5) * t334) + Ifges(3,1) * t335 + Ifges(3,5) * qJD(2) + t244;
t302 = mrSges(3,1) * t267 - mrSges(3,2) * t263;
t418 = t451 * t263 + mrSges(2,1) + t302;
t328 = qJD(4) * t262;
t55 = pkin(3) * t218 - pkin(7) * t173 + t86;
t61 = pkin(7) * t172 + t87;
t11 = t262 * t55 + t266 * t61 + t96 * t327 - t328 * t99;
t10 = pkin(8) * t68 + t11;
t12 = -qJD(4) * t49 - t262 * t61 + t266 * t55;
t9 = pkin(4) * t211 - pkin(8) * t67 + t12;
t2 = qJD(5) * t13 + t10 * t265 + t261 * t9;
t3 = -qJD(5) * t14 - t10 * t261 + t265 * t9;
t416 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t415 = -t12 * mrSges(5,1) + t11 * mrSges(5,2);
t372 = pkin(4) * t248;
t375 = pkin(3) * t258;
t414 = -m(5) * t375 - m(6) * (t372 + t375) + mrSges(2,2) - mrSges(3,3) - t299;
t362 = Ifges(3,4) * t263;
t296 = Ifges(3,2) * t267 + t362;
t413 = t13 * mrSges(6,1) + t48 * mrSges(5,1) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t296 / 0.2e1 - t14 * mrSges(6,2) - t49 * mrSges(5,2);
t411 = Ifges(6,4) * t408 + Ifges(6,2) * t407 + Ifges(6,6) * t385;
t410 = Ifges(6,1) * t408 + Ifges(6,4) * t407 + Ifges(6,5) * t385;
t406 = Ifges(5,4) * t400 + Ifges(5,2) * t399 + Ifges(5,6) * t384;
t405 = Ifges(5,1) * t400 + Ifges(5,4) * t399 + Ifges(5,5) * t384;
t376 = Ifges(6,4) * t80;
t36 = Ifges(6,2) * t443 + Ifges(6,6) * t233 + t376;
t404 = -t36 / 0.2e1;
t403 = t36 / 0.2e1;
t74 = Ifges(6,4) * t443;
t37 = Ifges(6,1) * t80 + Ifges(6,5) * t233 + t74;
t402 = -t37 / 0.2e1;
t401 = t37 / 0.2e1;
t396 = -t443 / 0.2e1;
t395 = t443 / 0.2e1;
t394 = -t80 / 0.2e1;
t393 = t80 / 0.2e1;
t392 = Ifges(4,1) * t386 + Ifges(4,4) * t387 + Ifges(4,5) * t383;
t391 = -t303 / 0.2e1;
t390 = t303 / 0.2e1;
t389 = -t138 / 0.2e1;
t388 = t138 / 0.2e1;
t382 = -t233 / 0.2e1;
t381 = t233 / 0.2e1;
t380 = -t238 / 0.2e1;
t379 = t238 / 0.2e1;
t374 = pkin(4) * t138;
t369 = g(3) * t263;
t252 = t263 * pkin(6);
t364 = mrSges(5,3) * t303;
t363 = mrSges(5,3) * t138;
t361 = Ifges(3,4) * t267;
t186 = -qJDD(2) * pkin(2) + qJDD(3) + t205;
t348 = t186 * t263;
t340 = t264 * t267;
t339 = t267 * t268;
t162 = t239 * t340 + t240 * t268;
t163 = t239 * t268 - t240 * t340;
t338 = -t162 * mrSges(6,1) + t163 * mrSges(6,2);
t164 = -t239 * t339 + t264 * t240;
t165 = t264 * t239 + t240 * t339;
t337 = t164 * mrSges(6,1) - t165 * mrSges(6,2);
t187 = qJD(2) * t292 - t329;
t333 = qJD(2) * t263;
t321 = pkin(6) * t333;
t148 = t259 * t187 + t258 * t321;
t332 = qJD(2) * t267;
t247 = pkin(6) * t332;
t317 = t258 * t332;
t200 = pkin(3) * t317 + t247;
t217 = pkin(3) * t344 + t252;
t323 = Ifges(6,5) * t22 + Ifges(6,6) * t23 + Ifges(6,3) * t198;
t322 = Ifges(5,5) * t67 + Ifges(5,6) * t68 + Ifges(5,3) * t211;
t316 = m(4) * qJ(3) + mrSges(4,3);
t34 = -t68 * mrSges(5,1) + t67 * mrSges(5,2);
t8 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t306 = -t325 / 0.2e1;
t102 = -t172 * mrSges(4,1) + t173 * mrSges(4,2);
t84 = t266 * t142 - t150 * t262;
t301 = mrSges(3,1) * t263 + mrSges(3,2) * t267;
t298 = -mrSges(6,1) * t239 - mrSges(6,2) * t240;
t294 = Ifges(3,5) * t267 - Ifges(3,6) * t263;
t293 = Ifges(4,5) * t259 - Ifges(4,6) * t258;
t185 = t289 * t263;
t56 = -pkin(4) * t267 + pkin(8) * t185 + t84;
t184 = t213 * t263;
t60 = -pkin(8) * t184 + t85;
t32 = -t261 * t60 + t265 * t56;
t33 = t261 * t56 + t265 * t60;
t107 = -t184 * t265 + t185 * t261;
t108 = -t184 * t261 - t185 * t265;
t139 = -t213 * t261 - t265 * t289;
t140 = t213 * t265 - t261 * t289;
t291 = t215 * t267 - t256 * t263;
t290 = t241 * t267 + t263 * t365;
t287 = t323 - t416;
t286 = pkin(1) * t301;
t180 = -t248 * t339 + t264 * t249;
t178 = t248 * t340 + t249 * t268;
t285 = t263 * (Ifges(3,1) * t267 - t362);
t109 = qJD(2) * t288 + t148;
t175 = t258 * t187;
t121 = qJD(2) * t284 + t175;
t38 = t262 * t109 + t266 * t121 + t142 * t327 - t150 * t328;
t117 = -pkin(3) * t172 + t186;
t273 = t267 * (Ifges(4,3) * t263 + t267 * t293);
t39 = -qJD(4) * t85 + t266 * t109 - t121 * t262;
t272 = t263 * t316 + t433;
t227 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t334;
t181 = t264 * t248 + t249 * t339;
t179 = t248 * t268 - t249 * t340;
t174 = pkin(4) * t289 - t241;
t169 = -mrSges(4,1) * t334 - mrSges(4,3) * t207;
t168 = mrSges(4,2) * t334 + mrSges(4,3) * t206;
t166 = -pkin(6) * t343 + t203;
t155 = -pkin(6) * t318 + t192;
t149 = -t259 * t321 + t175;
t146 = pkin(4) * t184 + t217;
t130 = t207 * Ifges(4,4) + t206 * Ifges(4,2) - Ifges(4,6) * t334;
t126 = mrSges(4,1) * t218 - mrSges(4,3) * t173;
t125 = -mrSges(4,2) * t218 + mrSges(4,3) * t172;
t115 = -qJD(2) * t283 + t190 * t263;
t114 = -qJD(2) * t282 - t191 * t263;
t106 = mrSges(5,1) * t238 - t363;
t105 = -mrSges(5,2) * t238 + t364;
t98 = -t170 * t261 - t171 * t265;
t97 = -t170 * t265 + t171 * t261;
t93 = -pkin(4) * t115 + t200;
t91 = t173 * Ifges(4,4) + t172 * Ifges(4,2) + t218 * Ifges(4,6);
t83 = -mrSges(5,1) * t303 + mrSges(5,2) * t138;
t63 = mrSges(6,1) * t233 - mrSges(6,3) * t80;
t62 = -mrSges(6,2) * t233 + mrSges(6,3) * t443;
t53 = -mrSges(5,2) * t211 + mrSges(5,3) * t68;
t52 = mrSges(5,1) * t211 - mrSges(5,3) * t67;
t46 = -pkin(4) * t68 + t117;
t45 = -qJD(5) * t108 - t114 * t261 + t115 * t265;
t44 = qJD(5) * t107 + t114 * t265 + t115 * t261;
t40 = -mrSges(6,1) * t443 + mrSges(6,2) * t80;
t29 = pkin(8) * t115 + t38;
t28 = pkin(4) * t333 - pkin(8) * t114 + t39;
t18 = -mrSges(6,2) * t198 + mrSges(6,3) * t23;
t17 = mrSges(6,1) * t198 - mrSges(6,3) * t22;
t16 = t265 * t42 - t354;
t15 = -t261 * t42 - t351;
t5 = -qJD(5) * t33 - t261 * t29 + t265 * t28;
t4 = qJD(5) * t32 + t261 * t28 + t265 * t29;
t1 = [(Ifges(5,5) * t388 + Ifges(6,5) * t393 + Ifges(5,6) * t390 + Ifges(6,6) * t395 + Ifges(5,3) * t379 + Ifges(6,3) * t381 + t413 + t424 / 0.2e1) * t333 + (-t342 * t86 - t344 * t87) * mrSges(4,3) + t419 * t332 / 0.2e1 + t302 * t350 - t227 * t321 + (Ifges(5,5) * t114 + Ifges(5,6) * t115) * t379 - t286 * t325 + t102 * t252 + (Ifges(6,4) * t108 + Ifges(6,2) * t107) * t407 + (Ifges(6,4) * t44 + Ifges(6,2) * t45) * t395 + t299 * t348 + t45 * t403 - t185 * t405 - t184 * t406 + t108 * t410 + t107 * t411 + (Ifges(5,4) * t114 + Ifges(5,2) * t115) * t390 + (Ifges(6,5) * t44 + Ifges(6,6) * t45) * t381 + (Ifges(6,5) * t108 + Ifges(6,6) * t107) * t385 + (-mrSges(3,1) * t252 + Ifges(3,5) * t263 + (-mrSges(3,2) * pkin(6) + Ifges(3,6)) * t267) * qJDD(2) - (Ifges(4,5) * t173 + Ifges(4,6) * t172 + Ifges(4,3) * t218 + t322 + t323) * t267 / 0.2e1 + (-Ifges(5,5) * t185 - Ifges(5,6) * t184) * t384 + (-Ifges(5,4) * t185 - Ifges(5,2) * t184) * t399 + (-t11 * t184 - t114 * t48 + t115 * t49 + t12 * t185) * mrSges(5,3) + t117 * (mrSges(5,1) * t184 - mrSges(5,2) * t185) + (-Ifges(5,1) * t185 - Ifges(5,4) * t184) * t400 + m(6) * (t13 * t5 + t14 * t4 + t146 * t46 + t2 * t33 + t3 * t32 + t90 * t93) + m(5) * (t11 * t85 + t117 * t217 + t12 * t84 + t153 * t200 + t38 * t49 + t39 * t48) + m(4) * (t143 * t148 + t144 * t149 + t166 * t86 + t167 * t87 + (t221 * t332 + t348) * pkin(6)) + (t107 * t2 - t108 * t3 - t13 * t44 + t14 * t45) * mrSges(6,3) + (t285 + t267 * (-Ifges(3,2) * t263 + t361)) * t325 / 0.2e1 - t91 * t344 / 0.2e1 + (-pkin(6) * t218 * t267 + t219 * t252 + t428) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t428) + t430 * t247 - t130 * t317 / 0.2e1 + t429 * t440 + t361 * t441 + t296 * t442 + Ifges(2,3) * qJDD(1) + (t294 * t440 - t422) * qJD(2) + (-t86 * mrSges(4,1) + t87 * mrSges(4,2) + Ifges(3,4) * t441 - Ifges(4,5) * t386 - Ifges(5,5) * t400 - Ifges(6,5) * t408 + Ifges(3,2) * t442 - Ifges(4,6) * t387 - Ifges(5,6) * t399 - Ifges(6,6) * t407 - Ifges(4,3) * t383 - Ifges(5,3) * t384 - Ifges(6,3) * t385 + t415 + t416) * t267 + (Ifges(3,1) * t219 + Ifges(3,4) * t442 + t293 * t383 + t295 * t387 + t297 * t386) * t263 + t32 * t17 + t33 * t18 + (-t181 * mrSges(5,1) - t165 * mrSges(6,1) - t180 * mrSges(5,2) - t164 * mrSges(6,2) + t423 * (t268 * pkin(1) + t264 * pkin(6)) + t414 * t264 + (-m(5) * t290 - m(6) * t291 - t272 - t418) * t268) * g(2) + t4 * t62 + t5 * t63 + t273 * t306 + t84 * t52 + t85 * t53 + t90 * (-mrSges(6,1) * t45 + mrSges(6,2) * t44) + t93 * t40 + t38 * t105 + t39 * t106 + (-t179 * mrSges(5,1) - t163 * mrSges(6,1) - t178 * mrSges(5,2) - t162 * mrSges(6,2) + (-m(5) * (-pkin(1) - t290) - m(4) * t307 + t263 * mrSges(4,3) + t433 + m(3) * pkin(1) - m(6) * (-pkin(1) - t291) + t418) * t264 + (pkin(6) * t423 + t414) * t268) * g(1) + t46 * (-mrSges(6,1) * t107 + mrSges(6,2) * t108) + t146 * t8 + t153 * (-mrSges(5,1) * t115 + mrSges(5,2) * t114) + t166 * t126 + t167 * t125 + t149 * t168 + t148 * t169 + t200 * t83 + t217 * t34 - pkin(1) * (mrSges(3,1) * t218 + mrSges(3,2) * t219) + (Ifges(5,1) * t114 + Ifges(5,4) * t115) * t388 + t342 * t392 + t114 * t397 + t115 * t398 + t44 * t401 + (Ifges(6,1) * t108 + Ifges(6,4) * t107) * t408 + (Ifges(6,1) * t44 + Ifges(6,4) * t45) * t393; (t259 * t125 - t258 * t126) * qJ(3) + (t422 - t429 / 0.2e1 + (t273 / 0.2e1 - t285 / 0.2e1 + t286) * qJD(1)) * qJD(1) + (t13 * t98 + t139 * t2 - t14 * t97 - t140 * t3) * mrSges(6,3) + (Ifges(5,5) * t389 + Ifges(6,5) * t394 + Ifges(5,6) * t391 + Ifges(6,6) * t396 + Ifges(5,3) * t380 + Ifges(6,3) * t382 - t413) * t335 + (Ifges(6,1) * t393 + Ifges(6,4) * t395 + Ifges(6,5) * t381 + t401 + t420) * (qJD(5) * t139 - t190 * t265 - t191 * t261) + (Ifges(6,4) * t393 + Ifges(6,2) * t395 + Ifges(6,6) * t381 + t403 + t421) * (-qJD(5) * t140 + t190 * t261 - t191 * t265) - (-Ifges(3,2) * t335 + t244 + t419) * t334 / 0.2e1 + t444 * t301 + (-t302 - t272) * g(3) - t289 * t406 + (Ifges(5,5) * t213 - Ifges(5,6) * t289) * t384 + t117 * (mrSges(5,1) * t289 + mrSges(5,2) * t213) + (Ifges(5,4) * t213 - Ifges(5,2) * t289) * t399 + (Ifges(5,1) * t213 - Ifges(5,4) * t289) * t400 + (-t331 - t154) * t169 + (Ifges(6,5) * t98 + Ifges(6,6) * t97) * t382 + (Ifges(4,5) * t258 + Ifges(4,6) * t259) * t383 + (Ifges(6,5) * t140 + Ifges(6,6) * t139) * t385 + t435 * t105 + (t11 * t152 - t117 * t241 + t12 * t151 - t153 * t199 + t435 * t49 + t436 * t48) * m(5) + t436 * t106 + t186 * t300 + t98 * t402 + t97 * t404 + t213 * t405 + (Ifges(6,4) * t140 + Ifges(6,2) * t139) * t407 + (Ifges(6,1) * t140 + Ifges(6,4) * t139) * t408 + t140 * t410 + t139 * t411 + t130 * t319 / 0.2e1 + (t444 * (t281 - t446) + t445 * g(3)) * t263 + (t444 * (-t316 + t445) + t446 * g(3)) * t267 + t447 * t397 + t448 * t398 + t426 * t40 + (-t143 * t154 - t144 * t155 - t221 * t246 - pkin(2) * t186 + (-t143 * t258 + t144 * t259) * qJD(3) + t427 * qJ(3)) * m(4) + t427 * mrSges(4,3) - t430 * t246 + (t330 - t155) * t168 + (-Ifges(5,5) * t190 - Ifges(5,6) * t191) * t379 + (-Ifges(5,1) * t190 - Ifges(5,4) * t191) * t388 + (-Ifges(5,4) * t190 - Ifges(5,2) * t191) * t390 + (-t11 * t289 - t12 * t213 - t447 * t48 + t448 * t49) * mrSges(5,3) + (-mrSges(5,1) * t448 + mrSges(5,2) * t447) * t153 + t227 * t245 + (Ifges(6,4) * t98 + Ifges(6,2) * t97) * t396 + t294 * t306 + t438 * t63 + (t13 * t438 + t14 * t439 + t174 * t46 + t2 * t58 + t3 * t57 + t426 * t90) * m(6) + t439 * t62 - t424 * t335 / 0.2e1 + (-Ifges(5,5) * t171 - Ifges(5,6) * t170) * t380 + (-Ifges(5,4) * t171 - Ifges(5,2) * t170) * t391 + (-Ifges(5,1) * t171 - Ifges(5,4) * t170) * t389 + t57 * t17 + t58 * t18 - t90 * (-mrSges(6,1) * t97 + mrSges(6,2) * t98) - pkin(2) * t102 + t46 * (-mrSges(6,1) * t139 + mrSges(6,2) * t140) + t151 * t52 + t152 * t53 + t174 * t8 - t199 * t83 - t204 * mrSges(3,2) - t205 * mrSges(3,1) - Ifges(3,6) * t218 + Ifges(3,5) * t219 - t241 * t34 + (Ifges(6,1) * t98 + Ifges(6,4) * t97) * t394 + t259 * t91 / 0.2e1 + Ifges(3,3) * qJDD(2) + (Ifges(4,1) * t258 + t359) * t386 + (Ifges(4,2) * t259 + t360) * t387 + t258 * t392; -t303 * t105 + t138 * t106 - t206 * t168 + t207 * t169 - t443 * t62 + t80 * t63 + t102 + t34 + t8 + (t13 * t80 - t14 * t443 + t46) * m(6) + (t138 * t48 - t303 * t49 + t117) * m(5) + (t143 * t207 - t144 * t206 + t186) * m(4) + (t267 * g(3) - t263 * t444) * t425; -t415 + (-Ifges(5,2) * t138 + t134 + t71) * t391 + (mrSges(5,2) * t181 - t180 * t437 - t337) * g(1) - t40 * t374 - m(6) * (t13 * t15 + t14 * t16 + t374 * t90) + (t2 * t261 + t265 * t3 + (-t13 * t261 + t14 * t265) * qJD(5)) * t409 + (t363 + t106) * t49 + (Ifges(5,5) * t303 - Ifges(5,6) * t138) * t380 - t153 * (mrSges(5,1) * t138 + mrSges(5,2) * t303) + (Ifges(5,1) * t303 - t358) * t389 + (m(6) * t372 + mrSges(5,1) * t248 + mrSges(5,2) * t249 - t298) * t369 + (-mrSges(5,2) * t179 + t178 * t437 - t338) * g(2) + t287 - (Ifges(6,4) * t394 + Ifges(6,2) * t396 + Ifges(6,6) * t382 + t404 - t421) * t80 + (Ifges(6,1) * t394 + Ifges(6,4) * t396 + Ifges(6,5) * t382 + t402 - t420) * t443 + t322 + (t364 - t105) * t48 - t16 * t62 - t15 * t63 + t70 * t388 + ((-t261 * t63 + t265 * t62) * qJD(5) + t17 * t265 + t18 * t261) * pkin(4); -t90 * (mrSges(6,1) * t80 + mrSges(6,2) * t443) + (Ifges(6,1) * t443 - t376) * t394 + t36 * t393 + (Ifges(6,5) * t443 - Ifges(6,6) * t80) * t382 - t13 * t62 + t14 * t63 - g(1) * t337 - g(2) * t338 - t298 * t369 + (t13 * t443 + t14 * t80) * mrSges(6,3) + t287 + (-Ifges(6,2) * t80 + t37 + t74) * t396;];
tau = t1;
