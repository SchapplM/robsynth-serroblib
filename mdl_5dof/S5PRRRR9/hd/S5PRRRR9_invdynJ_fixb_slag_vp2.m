% Calculate vector of inverse dynamics joint torques for
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:03
% EndTime: 2019-12-05 17:19:44
% DurationCPUTime: 15.47s
% Computational Cost: add. (5287->602), mult. (12275->860), div. (0->0), fcn. (9031->14), ass. (0->282)
t206 = sin(qJ(3));
t210 = cos(qJ(3));
t247 = mrSges(4,1) * t210 - mrSges(4,2) * t206;
t212 = -pkin(9) - pkin(8);
t374 = -m(5) * pkin(8) + m(6) * t212 - mrSges(5,3) - mrSges(6,3);
t209 = cos(qJ(4));
t196 = pkin(4) * t209 + pkin(3);
t200 = qJ(4) + qJ(5);
t198 = sin(t200);
t199 = cos(t200);
t205 = sin(qJ(4));
t246 = -mrSges(5,1) * t209 + mrSges(5,2) * t205;
t394 = m(5) * pkin(3) + m(6) * t196 + mrSges(6,1) * t199 - mrSges(6,2) * t198 - t246;
t403 = t206 * t374 - t394 * t210 - mrSges(3,1) - t247;
t248 = pkin(3) * t206 - pkin(8) * t210;
t168 = t248 * qJD(3);
t174 = -pkin(3) * t210 - pkin(8) * t206 - pkin(2);
t207 = sin(qJ(2));
t283 = qJD(4) * t209;
t285 = qJD(4) * t205;
t287 = qJD(3) * t209;
t202 = sin(pkin(5));
t294 = qJD(1) * t202;
t211 = cos(qJ(2));
t297 = t210 * t211;
t378 = t205 * t168 + t174 * t283 + (-t206 * t287 - t210 * t285) * pkin(7) - (t205 * t207 + t209 * t297) * t294;
t288 = qJD(3) * t206;
t276 = pkin(7) * t288;
t402 = t209 * t168 + t205 * t276 - (-t205 * t297 + t207 * t209) * t294;
t400 = -m(4) - m(5);
t298 = t209 * t210;
t193 = pkin(7) * t298;
t232 = pkin(4) * t206 - pkin(9) * t298;
t399 = t232 * qJD(3) + (-t193 + (pkin(9) * t206 - t174) * t205) * qJD(4) + t402;
t286 = qJD(3) * t210;
t263 = t205 * t286;
t220 = t206 * t283 + t263;
t398 = pkin(9) * t220 - t378;
t289 = qJD(2) * t210;
t264 = t205 * t289;
t270 = qJD(4) * t212;
t269 = t207 * t294;
t172 = qJD(2) * pkin(7) + t269;
t203 = cos(pkin(5));
t293 = qJD(1) * t203;
t119 = -t206 * t172 + t210 * t293;
t165 = t248 * qJD(2);
t75 = t209 * t119 + t205 * t165;
t397 = pkin(9) * t264 + t205 * t270 - t75;
t74 = -t119 * t205 + t209 * t165;
t396 = -qJD(2) * t232 + t209 * t270 - t74;
t292 = qJD(2) * t202;
t261 = qJD(1) * t292;
t183 = t211 * t261;
t281 = qJDD(1) * t202;
t134 = t207 * t281 + t183;
t393 = qJDD(2) * pkin(7) + qJD(3) * t293 + t134;
t370 = m(6) - t400;
t392 = pkin(2) * t370 - t403;
t245 = t205 * mrSges(5,1) + t209 * mrSges(5,2);
t391 = -t198 * mrSges(6,1) - t199 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - t245;
t204 = sin(qJ(5));
t208 = cos(qJ(5));
t291 = qJD(2) * t206;
t160 = -t205 * t291 + t287;
t161 = qJD(3) * t205 + t209 * t291;
t250 = t208 * t160 - t161 * t204;
t282 = qJD(2) * qJD(3);
t170 = qJDD(2) * t206 + t210 * t282;
t82 = qJD(4) * t160 + qJDD(3) * t205 + t170 * t209;
t83 = -qJD(4) * t161 + qJDD(3) * t209 - t170 * t205;
t25 = qJD(5) * t250 + t204 * t83 + t208 * t82;
t356 = t25 / 0.2e1;
t90 = t160 * t204 + t161 * t208;
t26 = -qJD(5) * t90 - t204 * t82 + t208 * t83;
t355 = t26 / 0.2e1;
t348 = t82 / 0.2e1;
t347 = t83 / 0.2e1;
t169 = qJDD(2) * t210 - t206 * t282;
t157 = qJDD(4) - t169;
t153 = qJDD(5) + t157;
t342 = t153 / 0.2e1;
t341 = t157 / 0.2e1;
t390 = t169 / 0.2e1;
t389 = t170 / 0.2e1;
t126 = t205 * t174 + t193;
t301 = t205 * t206;
t104 = -pkin(9) * t301 + t126;
t159 = t209 * t174;
t299 = t206 * t209;
t91 = -pkin(9) * t299 + t159 + (-pkin(7) * t205 - pkin(4)) * t210;
t41 = -t104 * t204 + t208 * t91;
t387 = qJD(5) * t41 + t399 * t204 - t398 * t208;
t42 = t104 * t208 + t204 * t91;
t386 = -qJD(5) * t42 + t398 * t204 + t399 * t208;
t180 = t212 * t205;
t181 = t212 * t209;
t112 = t180 * t204 - t181 * t208;
t385 = -qJD(5) * t112 - t397 * t204 + t396 * t208;
t111 = t180 * t208 + t181 * t204;
t384 = qJD(5) * t111 + t396 * t204 + t397 * t208;
t192 = qJD(4) - t289;
t186 = qJD(5) + t192;
t383 = t161 * Ifges(5,5) + t90 * Ifges(6,5) + t160 * Ifges(5,6) + Ifges(6,6) * t250 + t192 * Ifges(5,3) + t186 * Ifges(6,3);
t357 = m(6) * pkin(4);
t381 = -mrSges(5,1) - t357;
t267 = t206 * t293;
t332 = pkin(4) * t205;
t380 = pkin(4) * t285 - t267 - (qJD(2) * t332 + t172) * t210;
t34 = -mrSges(5,1) * t83 + mrSges(5,2) * t82;
t379 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t170 + t34;
t377 = -qJD(4) * t126 + t402;
t273 = mrSges(4,3) * t291;
t376 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t160 + mrSges(5,2) * t161 + t273;
t233 = t204 * t205 - t208 * t209;
t130 = t233 * t206;
t197 = Ifges(4,4) * t289;
t156 = Ifges(5,4) * t160;
t79 = Ifges(5,1) * t161 + Ifges(5,5) * t192 + t156;
t373 = Ifges(4,1) * t291 + Ifges(4,5) * qJD(3) + t209 * t79 + t197;
t280 = qJDD(1) * t203;
t57 = -t172 * t288 + t206 * t280 + t393 * t210;
t58 = -t172 * t286 - t393 * t206 + t210 * t280;
t372 = -t206 * t58 + t210 * t57;
t120 = t172 * t210 + t267;
t109 = qJD(3) * pkin(8) + t120;
t268 = t211 * t294;
t122 = qJD(2) * t174 - t268;
t47 = qJDD(3) * pkin(8) + t57;
t182 = t207 * t261;
t133 = t211 * t281 - t182;
t123 = -qJDD(2) * pkin(2) - t133;
t69 = -pkin(3) * t169 - pkin(8) * t170 + t123;
t13 = -t109 * t285 + t122 * t283 + t205 * t69 + t209 * t47;
t61 = t109 * t209 + t122 * t205;
t14 = -qJD(4) * t61 - t205 * t47 + t209 * t69;
t371 = t13 * t209 - t14 * t205;
t369 = qJD(4) + qJD(5);
t368 = mrSges(4,1) + t394;
t367 = mrSges(4,2) + t374;
t108 = -qJD(3) * pkin(3) - t119;
t364 = -m(5) * t108 - t376;
t10 = pkin(9) * t83 + t13;
t38 = pkin(9) * t160 + t61;
t314 = t204 * t38;
t60 = -t109 * t205 + t209 * t122;
t37 = -pkin(9) * t161 + t60;
t33 = pkin(4) * t192 + t37;
t15 = t208 * t33 - t314;
t7 = pkin(4) * t157 - pkin(9) * t82 + t14;
t2 = qJD(5) * t15 + t10 * t208 + t204 * t7;
t311 = t208 * t38;
t16 = t204 * t33 + t311;
t3 = -qJD(5) * t16 - t10 * t204 + t208 * t7;
t363 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t362 = -t14 * mrSges(5,1) + t13 * mrSges(5,2);
t271 = pkin(7) + t332;
t361 = -m(6) * t271 + t400 * pkin(7) + t391;
t322 = Ifges(4,4) * t206;
t241 = t210 * Ifges(4,2) + t322;
t360 = t16 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t241 / 0.2e1 - t15 * mrSges(6,1);
t213 = qJD(2) ^ 2;
t359 = Ifges(6,4) * t356 + Ifges(6,2) * t355 + Ifges(6,6) * t342;
t358 = Ifges(6,1) * t356 + Ifges(6,4) * t355 + Ifges(6,5) * t342;
t354 = Ifges(5,1) * t348 + Ifges(5,4) * t347 + Ifges(5,5) * t341;
t334 = Ifges(6,4) * t90;
t31 = Ifges(6,2) * t250 + Ifges(6,6) * t186 + t334;
t353 = -t31 / 0.2e1;
t352 = t31 / 0.2e1;
t84 = Ifges(6,4) * t250;
t32 = Ifges(6,1) * t90 + Ifges(6,5) * t186 + t84;
t351 = -t32 / 0.2e1;
t350 = t32 / 0.2e1;
t320 = Ifges(5,4) * t161;
t78 = Ifges(5,2) * t160 + Ifges(5,6) * t192 + t320;
t349 = -t78 / 0.2e1;
t346 = -t250 / 0.2e1;
t345 = t250 / 0.2e1;
t344 = -t90 / 0.2e1;
t343 = t90 / 0.2e1;
t339 = t161 / 0.2e1;
t338 = -t186 / 0.2e1;
t337 = t186 / 0.2e1;
t335 = mrSges(6,3) * t15;
t333 = pkin(4) * t161;
t330 = t16 * mrSges(6,3);
t308 = cos(pkin(10));
t252 = t308 * t211;
t201 = sin(pkin(10));
t306 = t201 * t207;
t138 = -t203 * t252 + t306;
t253 = t308 * t207;
t305 = t201 * t211;
t139 = t203 * t253 + t305;
t254 = t202 * t308;
t97 = t139 * t210 - t206 * t254;
t327 = (t138 * t199 - t198 * t97) * mrSges(6,1) + (-t138 * t198 - t199 * t97) * mrSges(6,2);
t140 = t203 * t305 + t253;
t141 = -t203 * t306 + t252;
t99 = t201 * t202 * t206 + t141 * t210;
t326 = (t140 * t199 - t198 * t99) * mrSges(6,1) + (-t140 * t198 - t199 * t99) * mrSges(6,2);
t303 = t202 * t210;
t145 = t203 * t206 + t207 * t303;
t302 = t202 * t211;
t325 = (-t145 * t198 - t199 * t302) * mrSges(6,1) + (-t145 * t199 + t198 * t302) * mrSges(6,2);
t324 = mrSges(5,3) * t160;
t323 = mrSges(5,3) * t161;
t321 = Ifges(4,4) * t210;
t319 = Ifges(5,4) * t205;
t318 = Ifges(5,4) * t209;
t48 = -qJDD(3) * pkin(3) - t58;
t313 = t206 * t48;
t304 = t202 * t207;
t300 = t205 * t210;
t290 = qJD(2) * t207;
t284 = qJD(4) * t206;
t279 = Ifges(6,5) * t25 + Ifges(6,6) * t26 + Ifges(6,3) * t153;
t278 = Ifges(5,5) * t82 + Ifges(5,6) * t83 + Ifges(5,3) * t157;
t275 = pkin(7) * t286;
t272 = mrSges(4,3) * t289;
t266 = t202 * t290;
t265 = t211 * t292;
t251 = t282 / 0.2e1;
t243 = Ifges(5,1) * t209 - t319;
t242 = Ifges(5,1) * t205 + t318;
t240 = -Ifges(5,2) * t205 + t318;
t239 = Ifges(5,2) * t209 + t319;
t238 = Ifges(4,5) * t210 - Ifges(4,6) * t206;
t237 = Ifges(5,5) * t209 - Ifges(5,6) * t205;
t236 = Ifges(5,5) * t205 + Ifges(5,6) * t209;
t102 = -t145 * t205 - t209 * t302;
t230 = -t145 * t209 + t205 * t302;
t43 = t102 * t208 + t204 * t230;
t44 = t102 * t204 - t208 * t230;
t163 = t204 * t209 + t205 * t208;
t231 = t279 - t363;
t144 = -t203 * t210 + t206 * t304;
t173 = -qJD(2) * pkin(2) - t268;
t228 = t173 * (mrSges(4,1) * t206 + mrSges(4,2) * t210);
t227 = t206 * (Ifges(4,1) * t210 - t322);
t226 = t163 * t210;
t225 = t233 * t210;
t221 = -t205 * t284 + t209 * t286;
t219 = Ifges(5,5) * t206 + t210 * t243;
t218 = Ifges(5,6) * t206 + t210 * t240;
t217 = Ifges(5,3) * t206 + t210 * t237;
t94 = t369 * t163;
t176 = -qJD(3) * mrSges(4,2) + t272;
t171 = t271 * t206;
t164 = t247 * qJD(2);
t135 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t169;
t129 = t163 * t206;
t128 = qJD(2) * t225;
t127 = qJD(2) * t226;
t125 = -pkin(7) * t300 + t159;
t121 = pkin(4) * t220 + t275;
t118 = mrSges(5,1) * t192 - t323;
t117 = -mrSges(5,2) * t192 + t324;
t105 = -mrSges(4,1) * t169 + mrSges(4,2) * t170;
t101 = -qJD(3) * t144 + t210 * t265;
t100 = qJD(3) * t145 + t206 * t265;
t93 = t369 * t233;
t76 = -pkin(4) * t160 + t108;
t73 = mrSges(6,1) * t186 - mrSges(6,3) * t90;
t72 = -mrSges(6,2) * t186 + mrSges(6,3) * t250;
t65 = -mrSges(5,2) * t157 + mrSges(5,3) * t83;
t64 = mrSges(5,1) * t157 - mrSges(5,3) * t82;
t50 = -qJD(3) * t226 + t130 * t369;
t49 = -qJD(3) * t225 - t206 * t94;
t39 = -mrSges(6,1) * t250 + mrSges(6,2) * t90;
t36 = qJD(4) * t102 + t101 * t209 + t205 * t266;
t35 = qJD(4) * t230 - t101 * t205 + t209 * t266;
t28 = t82 * Ifges(5,4) + t83 * Ifges(5,2) + t157 * Ifges(5,6);
t27 = -pkin(4) * t83 + t48;
t20 = -mrSges(6,2) * t153 + mrSges(6,3) * t26;
t19 = mrSges(6,1) * t153 - mrSges(6,3) * t25;
t18 = t208 * t37 - t314;
t17 = -t204 * t37 - t311;
t9 = -qJD(5) * t44 - t204 * t36 + t208 * t35;
t8 = qJD(5) * t43 + t204 * t35 + t208 * t36;
t6 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t1 = [m(2) * qJDD(1) + t101 * t176 + t102 * t64 - t230 * t65 + t36 * t117 + t35 * t118 + t145 * t135 + t43 * t19 + t44 * t20 + t8 * t72 + t9 * t73 + (t6 + t379) * t144 + (t39 + t376) * t100 + (-m(2) - m(3) - t370) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t213 - t105) * t211 + (-mrSges(3,1) * t213 - mrSges(3,2) * qJDD(2) - qJD(2) * t164) * t207) * t202 + m(6) * (t100 * t76 + t144 * t27 + t15 * t9 + t16 * t8 + t2 * t44 + t3 * t43) + m(5) * (t100 * t108 + t102 * t14 - t13 * t230 + t144 * t48 + t35 * t60 + t36 * t61) + m(4) * (-t100 * t119 + t101 * t120 - t144 * t58 + t145 * t57 + (-t123 * t211 + t173 * t290) * t202) + m(3) * (qJDD(1) * t203 ^ 2 + (t133 * t211 + t134 * t207) * t202); (Ifges(6,5) * t49 + Ifges(6,6) * t50) * t337 + (-m(6) * t76 + t364 - t39) * t206 * t268 + (t392 * t140 + t141 * t361) * g(1) + (t392 * t138 + t139 * t361) * g(2) + (t60 * mrSges(5,1) - t61 * mrSges(5,2) + Ifges(6,5) * t343 + Ifges(6,6) * t345 + Ifges(6,3) * t337 - t360 - t120 * mrSges(4,3) + t383 / 0.2e1) * t288 + (-t119 * t286 + t372) * mrSges(4,3) - (t205 * t79 + t209 * t78) * t284 / 0.2e1 - (t278 + t279) * t210 / 0.2e1 + (-Ifges(6,5) * t130 - Ifges(6,6) * t129) * t342 + (-t2 * t129 + t3 * t130 - t15 * t49 + t16 * t50) * mrSges(6,3) + (-Ifges(6,1) * t130 - Ifges(6,4) * t129) * t356 + (-Ifges(6,4) * t130 - Ifges(6,2) * t129) * t355 + t27 * (mrSges(6,1) * t129 - mrSges(6,2) * t130) + t321 * t389 + t241 * t390 + (Ifges(6,1) * t49 + Ifges(6,4) * t50) * t343 + (-t13 * t301 - t14 * t299 - t220 * t61 - t221 * t60) * mrSges(5,3) - t28 * t301 / 0.2e1 + t299 * t354 - t130 * t358 - t129 * t359 + (qJD(3) * t219 - t242 * t284) * t339 + t263 * t349 + t49 * t350 + t50 * t352 + t245 * t313 + (Ifges(6,4) * t49 + Ifges(6,2) * t50) * t345 + t108 * (mrSges(5,1) * t220 + mrSges(5,2) * t221) + (t182 + t133) * mrSges(3,1) + (-t370 * (pkin(2) * t302 + pkin(7) * t304) + (t403 * t211 + (-m(6) * t332 + t391) * t207) * t202) * g(3) + t160 * (qJD(3) * t218 - t239 * t284) / 0.2e1 + t192 * (qJD(3) * t217 - t236 * t284) / 0.2e1 - t176 * t276 + Ifges(3,3) * qJDD(2) + t164 * t269 + t379 * pkin(7) * t206 + (-pkin(2) * t123 + ((-t119 * t210 - t120 * t206) * qJD(3) + t372) * pkin(7) - (t173 * t207 + (-t119 * t206 + t120 * t210) * t211) * t294) * m(4) + t373 * t286 / 0.2e1 + t376 * t275 + t377 * t118 + (t125 * t14 + t126 * t13 + (t108 * t286 + t313) * pkin(7) + t378 * t61 + t377 * t60) * m(5) + t378 * t117 + t386 * t73 + t387 * t72 + (t121 * t76 + t15 * t386 + t16 * t387 + t171 * t27 + t2 * t42 + t3 * t41) * m(6) + (Ifges(4,4) * t389 + Ifges(4,2) * t390 - Ifges(6,6) * t355 - Ifges(6,5) * t356 - Ifges(5,3) * t341 - Ifges(6,3) * t342 - Ifges(5,6) * t347 - Ifges(5,5) * t348 - t176 * t268 + pkin(7) * t135 + (-Ifges(4,2) * t206 + t321) * t251 + t362 + t363) * t210 + (Ifges(4,1) * t170 + Ifges(4,4) * t390 + t237 * t341 + t240 * t347 + t243 * t348) * t206 + t227 * t251 + qJD(3) * t228 + qJD(3) ^ 2 * t238 / 0.2e1 - t123 * t247 + t41 * t19 + t42 * t20 + (t183 - t134) * mrSges(3,2) + t76 * (-mrSges(6,1) * t50 + mrSges(6,2) * t49) - pkin(2) * t105 + t121 * t39 + t125 * t64 + t126 * t65 + t171 * t6 + qJDD(3) * (Ifges(4,5) * t206 + Ifges(4,6) * t210); (Ifges(6,4) * t163 - Ifges(6,2) * t233) * t355 + (Ifges(6,1) * t163 - Ifges(6,4) * t233) * t356 + (Ifges(6,5) * t163 - Ifges(6,6) * t233) * t342 + t27 * (mrSges(6,1) * t233 + mrSges(6,2) * t163) - t233 * t359 + (t127 * t16 - t128 * t15 - t163 * t3 - t2 * t233) * mrSges(6,3) + ((t128 - t93) * mrSges(6,2) + (-t127 + t94) * mrSges(6,1)) * t76 + (-Ifges(6,4) * t128 - Ifges(6,2) * t127) * t346 + (-Ifges(6,1) * t128 - Ifges(6,4) * t127) * t344 + (-Ifges(6,5) * t128 - Ifges(6,6) * t127) * t338 + (-t228 - t60 * (mrSges(5,1) * t206 - mrSges(5,3) * t298) - t61 * (-mrSges(5,2) * t206 - mrSges(5,3) * t300)) * qJD(2) + (t144 * t368 + t145 * t367) * g(3) + (t367 * t97 - t368 * (-t139 * t206 - t210 * t254)) * g(2) + (t367 * t99 - t368 * (-t141 * t206 + t201 * t303)) * g(1) + (t273 + t364) * t120 + (Ifges(6,5) * t344 + Ifges(6,6) * t346 + Ifges(6,3) * t338 + t360) * t291 + (t160 * t240 + t161 * t243 + t192 * t237) * qJD(4) / 0.2e1 - (t160 * t218 + t161 * t219 + t192 * t217) * qJD(2) / 0.2e1 + t205 * t354 + t163 * t358 + t236 * t341 + t239 * t347 + t242 * t348 + t285 * t349 - t128 * t351 - t127 * t353 + (-pkin(3) * t48 - t60 * t74 - t61 * t75) * m(5) - (Ifges(6,1) * t343 + Ifges(6,4) * t345 + Ifges(6,5) * t337 - t335 + t350) * t93 - t213 * t227 / 0.2e1 + (t272 - t176) * t119 - (Ifges(6,4) * t343 + Ifges(6,2) * t345 + Ifges(6,6) * t337 + t330 + t352) * t94 + t79 * t283 / 0.2e1 - t238 * t282 / 0.2e1 + Ifges(4,3) * qJDD(3) + t78 * t264 / 0.2e1 + t192 * t108 * t245 + (m(5) * ((-t205 * t61 - t209 * t60) * qJD(4) + t371) - t205 * t64 + t209 * t65 - t118 * t283 - t117 * t285) * pkin(8) + (-t283 * t60 - t285 * t61 + t371) * mrSges(5,3) - (-Ifges(4,2) * t291 + t197 + t373) * t289 / 0.2e1 + t380 * t39 - t383 * t291 / 0.2e1 + t384 * t72 + (t111 * t3 + t112 * t2 + t15 * t385 + t16 * t384 - t196 * t27 + t380 * t76) * m(6) + t385 * t73 + t48 * t246 - pkin(3) * t34 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t111 * t19 + t112 * t20 - t75 * t117 - t74 * t118 + Ifges(4,6) * t169 + Ifges(4,5) * t170 - t196 * t6 + t209 * t28 / 0.2e1; (-t326 - (-t140 * t205 - t209 * t99) * mrSges(5,2) + t381 * (t140 * t209 - t205 * t99)) * g(1) + (-t327 - (-t138 * t205 - t209 * t97) * mrSges(5,2) + t381 * (t138 * t209 - t205 * t97)) * g(2) - (-Ifges(5,2) * t161 + t156 + t79) * t160 / 0.2e1 + (-mrSges(5,2) * t230 + t102 * t381 - t325) * g(3) - (t76 * mrSges(6,1) + Ifges(6,4) * t344 + Ifges(6,2) * t346 + Ifges(6,6) * t338 - t330 + t353) * t90 + (-t76 * mrSges(6,2) + Ifges(6,1) * t344 + Ifges(6,4) * t346 + Ifges(6,5) * t338 + t335 + t351) * t250 + (t2 * t204 + t208 * t3 + (-t15 * t204 + t16 * t208) * qJD(5)) * t357 + t78 * t339 + (t324 - t117) * t60 - t39 * t333 - m(6) * (t15 * t17 + t16 * t18 + t333 * t76) - t161 * (Ifges(5,1) * t160 - t320) / 0.2e1 + (t323 + t118) * t61 + t278 + t231 - t362 - t18 * t72 - t17 * t73 - t108 * (mrSges(5,1) * t161 + mrSges(5,2) * t160) - t192 * (Ifges(5,5) * t160 - Ifges(5,6) * t161) / 0.2e1 + ((-t204 * t73 + t208 * t72) * qJD(5) + t208 * t19 + t204 * t20) * pkin(4); -t76 * (mrSges(6,1) * t90 + mrSges(6,2) * t250) + (Ifges(6,1) * t250 - t334) * t344 + t31 * t343 + (Ifges(6,5) * t250 - Ifges(6,6) * t90) * t338 - t15 * t72 + t16 * t73 - g(1) * t326 - g(2) * t327 - g(3) * t325 + (t15 * t250 + t16 * t90) * mrSges(6,3) + t231 + (-Ifges(6,2) * t90 + t32 + t84) * t346;];
tau = t1;
