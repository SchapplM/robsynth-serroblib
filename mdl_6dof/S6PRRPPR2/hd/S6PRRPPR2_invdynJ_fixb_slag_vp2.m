% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:09
% EndTime: 2019-03-08 21:04:36
% DurationCPUTime: 17.97s
% Computational Cost: add. (5700->668), mult. (13101->873), div. (0->0), fcn. (9764->14), ass. (0->318)
t402 = -m(6) - m(7);
t214 = sin(qJ(3));
t217 = cos(qJ(3));
t212 = -qJ(4) - pkin(8);
t269 = qJD(3) * t212;
t148 = qJD(4) * t217 + t214 * t269;
t149 = -qJD(4) * t214 + t217 * t269;
t208 = sin(pkin(11));
t330 = cos(pkin(11));
t164 = t208 * t217 + t214 * t330;
t218 = cos(qJ(2));
t210 = sin(pkin(6));
t304 = qJD(1) * t210;
t278 = t218 * t304;
t412 = -t148 * t208 + t330 * t149 + t164 * t278;
t411 = mrSges(5,2) - mrSges(6,3);
t399 = -mrSges(6,2) + mrSges(5,1);
t264 = t330 * t217;
t298 = qJD(3) * t214;
t157 = qJD(3) * t264 - t208 * t298;
t410 = pkin(5) * t157 - t412;
t155 = t164 * qJD(3);
t203 = pkin(3) * t298;
t232 = -qJ(5) * t157 - qJD(5) * t164 + t203;
t215 = sin(qJ(2));
t279 = t215 * t304;
t369 = pkin(4) + pkin(9);
t408 = -t155 * t369 - t232 + t279;
t213 = sin(qJ(6));
t216 = cos(qJ(6));
t251 = mrSges(7,1) * t213 + mrSges(7,2) * t216;
t211 = cos(pkin(6));
t331 = cos(pkin(10));
t267 = t331 * t215;
t209 = sin(pkin(10));
t314 = t209 * t218;
t151 = t211 * t267 + t314;
t268 = t210 * t331;
t225 = -t151 * t214 - t217 * t268;
t223 = t225 * pkin(3);
t397 = Ifges(6,5) - Ifges(5,6);
t302 = qJD(2) * t210;
t273 = qJD(1) * t302;
t185 = t218 * t273;
t292 = qJDD(1) * t210;
t144 = t215 * t292 + t185;
t133 = qJDD(2) * pkin(8) + t144;
t303 = qJD(1) * t211;
t405 = qJD(3) * t303 + t133;
t313 = t210 * t215;
t158 = t211 * t217 - t214 * t313;
t266 = t331 * t218;
t315 = t209 * t215;
t153 = -t211 * t315 + t266;
t312 = t210 * t217;
t404 = -t153 * t214 + t209 * t312;
t403 = m(7) * pkin(9);
t294 = qJD(2) * qJD(3);
t171 = qJDD(2) * t217 - t214 * t294;
t172 = qJDD(2) * t214 + t217 * t294;
t104 = -t171 * t330 + t172 * t208;
t301 = qJD(2) * t214;
t154 = -qJD(2) * t264 + t208 * t301;
t116 = -qJD(3) * t213 + t154 * t216;
t46 = qJD(6) * t116 + qJDD(3) * t216 + t104 * t213;
t372 = t46 / 0.2e1;
t117 = qJD(3) * t216 + t154 * t213;
t47 = -qJD(6) * t117 - qJDD(3) * t213 + t104 * t216;
t371 = t47 / 0.2e1;
t105 = t208 * t171 + t172 * t330;
t98 = qJDD(6) + t105;
t370 = t98 / 0.2e1;
t163 = t208 * t214 - t264;
t200 = pkin(3) * t217 + pkin(2);
t238 = -qJ(5) * t164 - t200;
t64 = t163 * t369 + t238;
t180 = t212 * t214;
t181 = t212 * t217;
t113 = -t330 * t180 - t181 * t208;
t70 = pkin(5) * t164 + t113;
t25 = t213 * t70 + t216 * t64;
t401 = -qJD(6) * t25 + t213 * t408 + t216 * t410;
t24 = -t213 * t64 + t216 * t70;
t400 = qJD(6) * t24 + t213 * t410 - t216 * t408;
t356 = g(3) * t210;
t352 = qJD(3) / 0.2e1;
t398 = Ifges(5,5) - Ifges(6,4);
t142 = Ifges(5,4) * t154;
t156 = t164 * qJD(2);
t145 = qJD(6) + t156;
t392 = t145 * Ifges(7,3);
t393 = t116 * Ifges(7,6);
t396 = t156 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t117 * Ifges(7,5) - t142 + t392 + t393;
t80 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t104;
t82 = mrSges(6,1) * t104 - qJDD(3) * mrSges(6,3);
t395 = t80 - t82;
t81 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t105;
t83 = t105 * mrSges(6,1) + qJDD(3) * mrSges(6,2);
t394 = -t81 + t83;
t333 = qJDD(3) / 0.2e1;
t350 = mrSges(6,1) * t154;
t130 = -qJD(3) * mrSges(6,3) + t350;
t61 = -mrSges(7,1) * t116 + mrSges(7,2) * t117;
t332 = t130 - t61;
t348 = mrSges(5,3) * t154;
t128 = -qJD(3) * mrSges(5,2) - t348;
t391 = t130 - t128;
t347 = mrSges(5,3) * t156;
t349 = mrSges(6,1) * t156;
t305 = -qJD(3) * t399 + t347 + t349;
t309 = t213 * t218;
t187 = t210 * t309;
t159 = t211 * t214 + t215 * t312;
t74 = -t158 * t330 + t159 * t208;
t62 = t216 * t74 + t187;
t206 = qJ(3) + pkin(11);
t204 = sin(t206);
t205 = cos(t206);
t390 = t204 * t411 - t205 * t399;
t175 = qJD(2) * pkin(8) + t279;
t291 = qJDD(1) * t211;
t57 = -t175 * t298 + t214 * t291 + t217 * t405;
t277 = t214 * t303;
t120 = t175 * t217 + t277;
t191 = t217 * t291;
t58 = -t120 * qJD(3) - t133 * t214 + t191;
t389 = -t214 * t58 + t217 * t57;
t22 = mrSges(7,1) * t98 - mrSges(7,3) * t46;
t23 = -mrSges(7,2) * t98 + mrSges(7,3) * t47;
t388 = t213 * t23 + t216 * t22;
t387 = m(5) - t402;
t254 = -mrSges(4,1) * t217 + mrSges(4,2) * t214;
t84 = mrSges(5,1) * t154 + mrSges(5,2) * t156;
t85 = -mrSges(6,2) * t154 - mrSges(6,3) * t156;
t386 = t254 * qJD(2) + t84 + t85;
t385 = 0.2e1 * t333;
t252 = t216 * mrSges(7,1) - t213 * mrSges(7,2);
t357 = pkin(5) * t154;
t194 = t217 * t303;
t261 = qJ(4) * qJD(2) + t175;
t108 = -t214 * t261 + t194;
t103 = qJD(3) * pkin(3) + t108;
t109 = t217 * t261 + t277;
t265 = t330 * t109;
t43 = t208 * t103 + t265;
t40 = -qJD(3) * qJ(5) - t43;
t27 = -t40 - t357;
t341 = t117 * Ifges(7,4);
t38 = t116 * Ifges(7,2) + t145 * Ifges(7,6) + t341;
t384 = t27 * t252 - t216 * t38 / 0.2e1;
t99 = t208 * t109;
t42 = t103 * t330 - t99;
t237 = qJD(5) - t42;
t354 = t156 * pkin(5);
t26 = -qJD(3) * t369 + t237 + t354;
t140 = -qJD(2) * t200 + qJD(4) - t278;
t221 = -qJ(5) * t156 + t140;
t41 = t154 * t369 + t221;
t10 = -t213 * t41 + t216 * t26;
t184 = t215 * t273;
t143 = t218 * t292 - t184;
t132 = -qJDD(2) * pkin(2) - t143;
t87 = -pkin(3) * t171 + qJDD(4) + t132;
t220 = -qJ(5) * t105 - qJD(5) * t156 + t87;
t15 = t104 * t369 + t220;
t293 = qJD(2) * qJD(4);
t297 = qJD(3) * t217;
t33 = -t175 * t297 + qJDD(3) * pkin(3) - qJ(4) * t172 + t191 + (-t293 - t405) * t214;
t34 = qJ(4) * t171 + t217 * t293 + t57;
t13 = -t208 * t34 + t33 * t330;
t236 = qJDD(5) - t13;
t5 = t105 * pkin(5) - qJDD(3) * t369 + t236;
t1 = qJD(6) * t10 + t15 * t216 + t213 * t5;
t11 = t213 * t26 + t216 * t41;
t2 = -qJD(6) * t11 - t15 * t213 + t216 * t5;
t383 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t382 = mrSges(7,3) + t399 + t403;
t60 = pkin(4) * t154 + t221;
t381 = -t140 * mrSges(5,1) + t60 * mrSges(6,2);
t380 = qJ(5) * t402 - t251 + t411;
t379 = m(5) * t43 - m(6) * t40 - t391;
t377 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(5,3) - t252;
t227 = m(4) * pkin(2) - t254;
t376 = t204 * t251 + mrSges(3,1) + t227 - t390;
t375 = t10 * mrSges(7,1) + t140 * mrSges(5,2) - t11 * mrSges(7,2) - t60 * mrSges(6,3);
t288 = m(4) * pkin(8) + mrSges(4,3);
t374 = mrSges(3,2) - t288 + t377;
t219 = qJD(2) ^ 2;
t373 = Ifges(7,1) * t372 + Ifges(7,4) * t371 + Ifges(7,5) * t370;
t368 = -t116 / 0.2e1;
t367 = -t117 / 0.2e1;
t366 = t117 / 0.2e1;
t365 = -t145 / 0.2e1;
t364 = -t154 / 0.2e1;
t363 = t154 / 0.2e1;
t362 = -t156 / 0.2e1;
t361 = t156 / 0.2e1;
t358 = pkin(3) * t208;
t355 = t1 * t213;
t353 = -qJD(3) / 0.2e1;
t14 = t208 * t33 + t330 * t34;
t346 = mrSges(7,3) * t216;
t345 = Ifges(4,4) * t214;
t344 = Ifges(4,4) * t217;
t343 = Ifges(7,4) * t213;
t342 = Ifges(7,4) * t216;
t340 = t156 * Ifges(5,4);
t339 = t156 * Ifges(6,6);
t329 = qJ(5) * t204;
t150 = -t211 * t266 + t315;
t328 = t150 * t205;
t152 = t211 * t314 + t267;
t326 = t152 * t205;
t324 = t155 * t213;
t323 = t155 * t216;
t322 = t156 * t213;
t321 = t163 * t213;
t320 = t163 * t216;
t317 = t205 * t218;
t316 = t209 * t210;
t311 = t210 * t218;
t308 = t216 * t218;
t307 = -t150 * t200 - t151 * t212;
t306 = -t152 * t200 - t153 * t212;
t300 = qJD(2) * t215;
t299 = qJD(2) * t217;
t296 = qJD(6) * t213;
t295 = qJD(6) * t216;
t290 = Ifges(7,5) * t46 + Ifges(7,6) * t47 + Ifges(7,3) * t98;
t202 = pkin(3) * t301;
t289 = t128 - t332;
t287 = mrSges(4,3) * t301;
t286 = mrSges(4,3) * t299;
t283 = t210 * t308;
t280 = t330 * pkin(3);
t276 = t210 * t300;
t275 = t218 * t302;
t270 = -t296 / 0.2e1;
t262 = qJ(5) * t154 + t202;
t51 = t108 * t208 + t265;
t260 = -pkin(4) * t328 - t150 * t329 + t307;
t259 = -pkin(4) * t326 - t152 * t329 + t306;
t258 = t214 * t278;
t257 = t217 * t278;
t256 = t404 * pkin(3);
t199 = -t280 - pkin(4);
t249 = Ifges(7,1) * t213 + t342;
t248 = t217 * Ifges(4,2) + t345;
t247 = Ifges(7,2) * t216 + t343;
t246 = Ifges(4,5) * t217 - Ifges(4,6) * t214;
t245 = Ifges(7,5) * t213 + Ifges(7,6) * t216;
t243 = t10 * t213 - t11 * t216;
t65 = -mrSges(7,2) * t145 + mrSges(7,3) * t116;
t66 = mrSges(7,1) * t145 - mrSges(7,3) * t117;
t242 = -t213 * t66 + t216 * t65;
t241 = -t213 * t65 - t216 * t66;
t239 = t158 * pkin(3);
t235 = -t213 * t74 + t283;
t234 = t163 * t295 + t324;
t233 = t163 * t296 - t323;
t176 = -qJD(2) * pkin(2) - t278;
t230 = t176 * (mrSges(4,1) * t214 + mrSges(4,2) * t217);
t229 = t214 * (Ifges(4,1) * t217 - t345);
t53 = t108 * t330 - t99;
t69 = t148 * t330 + t208 * t149;
t114 = t208 * t180 - t181 * t330;
t135 = t204 * t313 - t211 * t205;
t92 = t151 * t204 + t205 * t268;
t94 = t153 * t204 - t205 * t316;
t228 = -g(1) * t94 - g(2) * t92 - g(3) * t135;
t9 = -qJDD(3) * qJ(5) - qJD(3) * qJD(5) - t14;
t224 = -g(1) * t152 - g(2) * t150 + g(3) * t311;
t222 = -qJD(6) * t243 + t2 * t216 + t355;
t201 = Ifges(4,4) * t299;
t197 = qJ(5) + t358;
t179 = -qJD(3) * mrSges(4,2) + t286;
t178 = qJD(3) * mrSges(4,1) - t287;
t170 = t200 * t311;
t161 = Ifges(4,1) * t301 + Ifges(4,5) * qJD(3) + t201;
t160 = Ifges(4,6) * qJD(3) + qJD(2) * t248;
t147 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t172;
t146 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t171;
t141 = Ifges(6,6) * t154;
t119 = -t175 * t214 + t194;
t115 = Ifges(7,4) * t116;
t110 = -mrSges(4,1) * t171 + mrSges(4,2) * t172;
t107 = qJD(3) * t158 + t217 * t275;
t106 = -qJD(3) * t159 - t214 * t275;
t97 = t105 * mrSges(6,3);
t96 = t105 * mrSges(5,2);
t86 = pkin(4) * t163 + t238;
t78 = -t154 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t340;
t77 = Ifges(6,4) * qJD(3) - t156 * Ifges(6,2) + t141;
t76 = Ifges(6,5) * qJD(3) + t154 * Ifges(6,3) - t339;
t75 = t208 * t158 + t159 * t330;
t71 = -t163 * pkin(5) + t114;
t67 = pkin(4) * t156 + t262;
t59 = pkin(4) * t155 + t232;
t56 = -t155 * pkin(5) + t69;
t54 = t156 * t369 + t262;
t52 = t208 * t106 + t107 * t330;
t50 = -t106 * t330 + t107 * t208;
t49 = -t104 * mrSges(6,2) - t97;
t48 = t104 * mrSges(5,1) + t96;
t39 = t117 * Ifges(7,1) + t145 * Ifges(7,5) + t115;
t36 = -qJD(3) * pkin(4) + t237;
t30 = t53 - t354;
t29 = t51 - t357;
t21 = pkin(4) * t104 + t220;
t20 = qJD(6) * t62 + t213 * t50 + t216 * t276;
t19 = qJD(6) * t235 - t213 * t276 + t216 * t50;
t18 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t17 = t213 * t29 + t216 * t54;
t16 = -t213 * t54 + t216 * t29;
t12 = -qJDD(3) * pkin(4) + t236;
t7 = t46 * Ifges(7,4) + t47 * Ifges(7,2) + t98 * Ifges(7,6);
t6 = -pkin(5) * t104 - t9;
t3 = [m(2) * qJDD(1) + t106 * t178 + t107 * t179 + t159 * t146 + t158 * t147 + t19 * t66 + t20 * t65 + t62 * t22 - t235 * t23 + t394 * t74 + t305 * t50 + (t18 + t395) * t75 + t289 * t52 + (-m(2) - m(3) - m(4) - t387) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t219 - t110 - t48 - t49) * t218 + (-mrSges(3,1) * t219 - mrSges(3,2) * qJDD(2) + qJD(2) * t386) * t215) * t210 + m(3) * (qJDD(1) * t211 ^ 2 + (t143 * t218 + t144 * t215) * t210) + m(4) * (t106 * t119 + t107 * t120 + t158 * t58 + t159 * t57 + (-t132 * t218 + t176 * t300) * t210) + m(6) * (t12 * t74 + t36 * t50 - t40 * t52 - t75 * t9 + (-t21 * t218 + t300 * t60) * t210) + m(5) * (-t13 * t74 + t14 * t75 - t42 * t50 + t43 * t52 + (t140 * t300 - t218 * t87) * t210) + m(7) * (-t1 * t235 + t10 * t19 + t11 * t20 + t2 * t62 + t27 * t52 + t6 * t75); (-t179 * t298 - t178 * t297 + m(4) * ((-t119 * t217 - t120 * t214) * qJD(3) + t389) + t217 * t146 - t214 * t147) * pkin(8) + (-t119 * t297 - t120 * t298 + t389) * mrSges(4,3) + (-m(5) * t140 - m(6) * t60 - t386) * t279 + t172 * t344 / 0.2e1 + t214 * (Ifges(4,4) * t171 + Ifges(4,5) * qJDD(3)) / 0.2e1 + t400 * t65 + (t1 * t25 + t10 * t401 + t11 * t400 + t2 * t24 + t27 * t56 + t6 * t71) * m(7) + t401 * t66 + (t1 * t320 - t10 * t234 - t11 * t233 - t2 * t321 - t317 * t356) * mrSges(7,3) + t172 * t214 * Ifges(4,1) + (t215 * t356 - t144 + t185) * mrSges(3,2) + (-(t176 * t215 + (-t119 * t214 + t120 * t217) * t218) * t304 - pkin(2) * t132) * m(4) + (-m(7) * (-pkin(9) * t326 + t259) + mrSges(7,3) * t326 - m(5) * t306 - m(6) * t259 + t374 * t153 + t376 * t152) * g(1) + (-m(7) * (-pkin(9) * t328 + t260) + mrSges(7,3) * t328 - m(5) * t307 - m(6) * t260 + t374 * t151 + t376 * t150) * g(2) + (t396 / 0.2e1 + t392 / 0.2e1 + t393 / 0.2e1 - t77 / 0.2e1 + Ifges(5,1) * t361 - Ifges(6,2) * t362 - Ifges(6,6) * t363 + Ifges(5,4) * t364 + Ifges(7,5) * t366 - t42 * mrSges(5,3) + t36 * mrSges(6,1) + t398 * t352 + t375) * t157 + (-m(5) * t170 + t402 * (t170 + (pkin(4) * t205 + t329) * t311) + (-t317 * t403 + t390 * t218 + (-t309 * mrSges(7,1) - t308 * mrSges(7,2)) * t204 + (t212 * t387 + t377) * t215) * t210) * g(3) + t105 * Ifges(5,1) * t164 + (t12 * mrSges(6,1) + t87 * mrSges(5,2) - t13 * mrSges(5,3) - t21 * mrSges(6,3) + Ifges(5,5) * t333 + Ifges(7,5) * t372 + Ifges(7,6) * t371 + Ifges(7,3) * t370 + Ifges(6,2) * t105 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t104 - t385 * Ifges(6,4) + t383) * t164 + t145 * (Ifges(7,5) * t234 - Ifges(7,6) * t233) / 0.2e1 + t412 * (m(5) * t42 - m(6) * t36 - t305) + (-t218 * t356 + t143 + t184) * mrSges(3,1) + m(6) * (t21 * t86 + t59 * t60) + m(5) * (t140 * t203 - t200 * t87) + t84 * t203 + (Ifges(7,1) * t234 - Ifges(7,4) * t233) * t366 + (-Ifges(5,4) * t104 + Ifges(5,5) * qJDD(3) + t290) * t164 / 0.2e1 + (-m(5) * t13 + m(6) * t12 + t394) * t113 + (m(5) * t14 - m(6) * t9 + t395) * t114 + (t76 / 0.2e1 - t78 / 0.2e1 - Ifges(5,4) * t361 + Ifges(6,6) * t362 + Ifges(6,3) * t363 - Ifges(5,2) * t364 + t40 * mrSges(6,1) - t43 * mrSges(5,3) + t397 * t352 - t381) * t155 + (t246 * t352 + t230) * qJD(3) + (t217 * (-Ifges(4,2) * t214 + t344) + t229) * t294 / 0.2e1 + (qJD(6) * t39 + t7) * t320 / 0.2e1 + t178 * t258 + t27 * (mrSges(7,1) * t233 + mrSges(7,2) * t234) + t116 * (Ifges(7,4) * t234 - Ifges(7,2) * t233) / 0.2e1 + (-m(7) * t27 - t379 - t61) * (-t208 * t258 + t257 * t330) + t379 * t69 + Ifges(3,3) * qJDD(2) + t132 * t254 + t171 * t248 / 0.2e1 - (t215 * t288 + t218 * t227) * t356 + t38 * t323 / 0.2e1 + t39 * t324 / 0.2e1 + t217 * (Ifges(4,4) * t172 + Ifges(4,2) * t171 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t160 * t298 / 0.2e1 + t161 * t297 / 0.2e1 - t200 * t48 - pkin(2) * t110 + t59 * t85 + t86 * t49 + t71 * t18 + t56 * t61 + t24 * t22 + t25 * t23 + t321 * t373 + (Ifges(4,5) * t214 + Ifges(4,6) * t217) * t333 + (t87 * mrSges(5,1) + t9 * mrSges(6,1) - t21 * mrSges(6,2) - t14 * mrSges(5,3) + t245 * t370 + t247 * t371 + t249 * t372 - t6 * t252 + t38 * t270 + (-Ifges(5,4) - Ifges(6,6)) * t105 + (Ifges(5,2) + Ifges(6,3)) * t104 + t397 * t385) * t163 - t179 * t257; -t332 * qJD(5) + t391 * t53 + (m(7) * t222 + t295 * t65 - t296 * t66 + t388) * (-pkin(9) + t199) - t305 * t51 + (-m(5) * t239 - mrSges(4,1) * t158 + mrSges(4,2) * t159 + t402 * (-t135 * pkin(4) + t239) + t380 * (t204 * t211 + t205 * t313) + t382 * t135) * g(3) + t384 * qJD(6) + (-t11 * t295 - t355 + (t296 + t322) * t10) * mrSges(7,3) + (t12 * t199 - t197 * t9 - t36 * t51 - t60 * t67 + (-qJD(5) + t53) * t40) * m(6) + (-t10 * t16 - t11 * t17 + t197 * t6 + (-t30 + qJD(5)) * t27) * m(7) + (-t340 + t76) * t362 + (t339 + t78) * t361 + (t141 + t77) * t364 + (-m(5) * t223 - t225 * mrSges(4,1) - (-t151 * t217 + t214 * t268) * mrSges(4,2) + t380 * (t151 * t205 - t204 * t268) + t382 * t92 + t402 * (-t92 * pkin(4) + t223)) * g(2) + (-m(5) * t256 - t404 * mrSges(4,1) - (-t153 * t217 - t214 * t316) * mrSges(4,2) + t402 * (-t94 * pkin(4) + t256) + t380 * (t153 * t205 + t204 * t316) + t382 * t94) * g(1) + (t287 + t178) * t120 + (-t322 / 0.2e1 + t270) * t39 + (-Ifges(5,2) * t156 - t142 + t396) * t363 + (Ifges(6,3) * t364 - t11 * t346 + t245 * t365 + t247 * t368 + t249 * t367 + t353 * t397 + t381 + t384) * t156 + t397 * t104 + t398 * t105 + (-Ifges(5,1) * t362 - Ifges(7,5) * t367 + Ifges(6,2) * t361 - Ifges(7,6) * t368 - Ifges(7,3) * t365 - t353 * t398 + t375) * t154 - (-Ifges(4,2) * t301 + t161 + t201) * t299 / 0.2e1 - (t116 * t247 + t117 * t249 + t145 * t245) * qJD(6) / 0.2e1 - qJD(2) * t230 + (t18 - t82) * t197 + ((t13 * t330 + t14 * t208) * pkin(3) - t140 * t202 + t42 * t51 - t43 * t53) * m(5) - t219 * t229 / 0.2e1 + t6 * t251 + (t286 - t179) * t119 - t40 * t349 - t42 * t348 - t2 * t346 + t160 * t301 / 0.2e1 - t246 * t294 / 0.2e1 - t213 * t7 / 0.2e1 - t84 * t202 + t199 * t83 + Ifges(4,6) * t171 + Ifges(4,5) * t172 - t67 * t85 - t30 * t61 - t17 * t65 - t16 * t66 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t13 * mrSges(5,1) - t14 * mrSges(5,2) + (-Ifges(7,2) * t213 + t342) * t371 + (Ifges(7,1) * t216 - t343) * t372 + t216 * t373 + (Ifges(7,5) * t216 - Ifges(7,6) * t213) * t370 + t43 * t347 + t36 * t350 + t80 * t358 + t81 * t280 + t12 * mrSges(6,2) - t9 * mrSges(6,3) + (Ifges(6,1) + Ifges(4,3) + Ifges(5,3)) * qJDD(3); -t213 * t22 + t216 * t23 + t96 - t97 + t399 * t104 + t241 * qJD(6) + t289 * t154 + (t241 - t305) * t156 + (t1 * t216 + t154 * t27 - t2 * t213 + t224 - t145 * (t10 * t216 + t11 * t213)) * m(7) + (-t154 * t40 - t156 * t36 + t21 + t224) * m(6) + (t154 * t43 + t156 * t42 + t224 + t87) * m(5); t242 * qJD(6) + t332 * qJD(3) + (t242 + t85) * t156 + t83 + (-qJD(3) * t27 - t156 * t243 + t222 + t228) * m(7) + (qJD(3) * t40 + t156 * t60 + t12 + t228) * m(6) + t388; -t27 * (mrSges(7,1) * t117 + mrSges(7,2) * t116) + (Ifges(7,1) * t116 - t341) * t367 + t38 * t366 + (Ifges(7,5) * t116 - Ifges(7,6) * t117) * t365 - t10 * t65 + t11 * t66 - g(1) * ((-t152 * t213 + t216 * t94) * mrSges(7,1) + (-t152 * t216 - t213 * t94) * mrSges(7,2)) - g(2) * ((-t150 * t213 + t216 * t92) * mrSges(7,1) + (-t150 * t216 - t213 * t92) * mrSges(7,2)) - g(3) * ((t135 * t216 + t187) * mrSges(7,1) + (-t135 * t213 + t283) * mrSges(7,2)) + (t10 * t116 + t11 * t117) * mrSges(7,3) + t290 + (-Ifges(7,2) * t117 + t115 + t39) * t368 + t383;];
tau  = t3;
