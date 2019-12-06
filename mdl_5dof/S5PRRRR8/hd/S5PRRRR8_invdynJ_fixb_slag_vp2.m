% Calculate vector of inverse dynamics joint torques for
% S5PRRRR8
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:50
% EndTime: 2019-12-05 17:15:31
% DurationCPUTime: 11.72s
% Computational Cost: add. (5310->536), mult. (12179->773), div. (0->0), fcn. (9087->14), ass. (0->263)
t388 = mrSges(5,2) - mrSges(6,3);
t214 = sin(qJ(5));
t218 = cos(qJ(5));
t383 = t218 * mrSges(6,1) - t214 * mrSges(6,2);
t387 = -mrSges(5,1) - t383;
t217 = sin(qJ(2));
t212 = sin(pkin(5));
t302 = qJD(1) * t212;
t272 = t217 * t302;
t216 = sin(qJ(3));
t298 = qJD(3) * t216;
t366 = pkin(3) * t298 - t272;
t215 = sin(qJ(4));
t219 = cos(qJ(4));
t220 = cos(qJ(3));
t174 = t215 * t216 - t219 * t220;
t232 = t174 * qJD(4);
t123 = -qJD(3) * t174 - t232;
t175 = t215 * t220 + t216 * t219;
t233 = t175 * qJD(4);
t124 = qJD(3) * t175 + t233;
t386 = pkin(4) * t124 - pkin(9) * t123 + t366;
t222 = -pkin(8) - pkin(7);
t273 = qJD(3) * t222;
t179 = t216 * t273;
t190 = t222 * t216;
t191 = t222 * t220;
t243 = t219 * t190 + t191 * t215;
t257 = t220 * t273;
t221 = cos(qJ(2));
t271 = t221 * t302;
t371 = qJD(4) * t243 + t174 * t271 + t219 * t179 + t215 * t257;
t136 = t190 * t215 - t191 * t219;
t370 = qJD(4) * t136 - t175 * t271 + t179 * t215 - t219 * t257;
t172 = t175 * qJD(2);
t209 = qJD(3) + qJD(4);
t133 = -t172 * t214 + t209 * t218;
t134 = t172 * t218 + t209 * t214;
t332 = t172 * mrSges(5,3);
t324 = -mrSges(5,1) * t209 - mrSges(6,1) * t133 + mrSges(6,2) * t134 + t332;
t210 = qJ(3) + qJ(4);
t206 = sin(t210);
t207 = cos(t210);
t357 = m(6) * pkin(9);
t358 = m(6) * pkin(4);
t385 = (-t358 + t387) * t207 + (-t357 + t388) * t206;
t290 = qJD(2) * qJD(3);
t180 = qJDD(2) * t220 - t216 * t290;
t181 = qJDD(2) * t216 + t220 * t290;
t83 = -qJD(2) * t233 + t180 * t219 - t181 * t215;
t213 = cos(pkin(5));
t312 = t212 * t217;
t167 = t213 * t220 - t216 * t312;
t323 = cos(pkin(10));
t260 = t323 * t221;
t211 = sin(pkin(10));
t314 = t211 * t217;
t166 = -t213 * t314 + t260;
t311 = t212 * t220;
t384 = -t166 * t216 + t211 * t311;
t208 = qJDD(3) + qJDD(4);
t82 = -qJD(2) * t232 + t180 * t215 + t181 * t219;
t38 = qJD(5) * t133 + t208 * t214 + t218 * t82;
t356 = t38 / 0.2e1;
t39 = -qJD(5) * t134 + t208 * t218 - t214 * t82;
t355 = t39 / 0.2e1;
t171 = t174 * qJD(2);
t162 = qJD(5) + t171;
t333 = t134 * Ifges(6,4);
t49 = t133 * Ifges(6,2) + t162 * Ifges(6,6) + t333;
t382 = -t49 / 0.2e1;
t81 = qJDD(5) - t83;
t353 = t81 / 0.2e1;
t381 = m(5) + m(6);
t380 = t180 / 0.2e1;
t379 = t181 / 0.2e1;
t203 = pkin(3) * t220 + pkin(2);
t114 = pkin(4) * t174 - pkin(9) * t175 - t203;
t59 = t114 * t218 - t136 * t214;
t378 = qJD(5) * t59 + t386 * t214 + t371 * t218;
t60 = t114 * t214 + t136 * t218;
t377 = -qJD(5) * t60 - t371 * t214 + t386 * t218;
t16 = -mrSges(6,1) * t39 + mrSges(6,2) * t38;
t71 = mrSges(5,1) * t208 - mrSges(5,3) * t82;
t376 = t16 - t71;
t375 = t220 * Ifges(4,2);
t252 = mrSges(6,1) * t214 + mrSges(6,2) * t218;
t301 = qJD(1) * t213;
t198 = t220 * t301;
t182 = qJD(2) * pkin(7) + t272;
t258 = pkin(8) * qJD(2) + t182;
t127 = -t216 * t258 + t198;
t121 = qJD(3) * pkin(3) + t127;
t270 = t216 * t301;
t128 = t220 * t258 + t270;
t321 = t128 * t215;
t57 = t121 * t219 - t321;
t53 = -pkin(4) * t209 - t57;
t374 = t252 * t53;
t373 = -mrSges(5,3) - t252;
t151 = -t206 * t312 + t207 * t213;
t152 = t206 * t213 + t207 * t312;
t369 = t387 * t151 + t388 * t152;
t315 = t211 * t212;
t117 = -t166 * t206 + t207 * t315;
t118 = t166 * t207 + t206 * t315;
t368 = t387 * t117 + t388 * t118;
t261 = t323 * t217;
t313 = t211 * t221;
t164 = t213 * t261 + t313;
t262 = t212 * t323;
t115 = -t164 * t206 - t207 * t262;
t116 = t164 * t207 - t206 * t262;
t367 = t387 * t115 + t388 * t116;
t293 = qJD(5) * t218;
t237 = t123 * t214 + t175 * t293;
t19 = mrSges(6,1) * t81 - mrSges(6,3) * t38;
t20 = -mrSges(6,2) * t81 + mrSges(6,3) * t39;
t365 = -t214 * t19 + t218 * t20;
t300 = qJD(2) * t212;
t266 = qJD(1) * t300;
t193 = t221 * t266;
t289 = qJDD(1) * t212;
t157 = t217 * t289 + t193;
t150 = qJDD(2) * pkin(7) + t157;
t288 = qJDD(1) * t213;
t66 = qJD(3) * t198 + t220 * t150 - t182 * t298 + t216 * t288;
t140 = t182 * t220 + t270;
t67 = -t140 * qJD(3) - t150 * t216 + t220 * t288;
t364 = -t216 * t67 + t220 * t66;
t52 = qJDD(3) * pkin(3) - pkin(8) * t181 + t67;
t55 = pkin(8) * t180 + t66;
t308 = t219 * t128;
t58 = t121 * t215 + t308;
t13 = -qJD(4) * t58 - t215 * t55 + t219 * t52;
t282 = m(4) * pkin(7) + mrSges(4,3);
t362 = mrSges(3,2) - t282 + t373;
t255 = -mrSges(4,1) * t220 + mrSges(4,2) * t216;
t231 = m(4) * pkin(2) - t255;
t361 = mrSges(3,1) + t231 - t385;
t192 = t217 * t266;
t156 = t221 * t289 - t192;
t149 = -qJDD(2) * pkin(2) - t156;
t113 = -pkin(3) * t180 + t149;
t21 = -pkin(4) * t83 - pkin(9) * t82 + t113;
t54 = pkin(9) * t209 + t58;
t161 = -qJD(2) * t203 - t271;
t78 = pkin(4) * t171 - pkin(9) * t172 + t161;
t22 = -t214 * t54 + t218 * t78;
t295 = qJD(4) * t219;
t296 = qJD(4) * t215;
t12 = t121 * t295 - t128 * t296 + t215 * t52 + t219 * t55;
t9 = pkin(9) * t208 + t12;
t2 = qJD(5) * t22 + t21 * t214 + t218 * t9;
t23 = t214 * t78 + t218 * t54;
t3 = -qJD(5) * t23 + t21 * t218 - t214 * t9;
t360 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t223 = qJD(2) ^ 2;
t359 = Ifges(6,1) * t356 + Ifges(6,4) * t355 + Ifges(6,5) * t353;
t352 = -t133 / 0.2e1;
t351 = -t134 / 0.2e1;
t350 = t134 / 0.2e1;
t349 = -t162 / 0.2e1;
t346 = t172 / 0.2e1;
t345 = t218 / 0.2e1;
t344 = pkin(3) * t215;
t343 = pkin(3) * t219;
t342 = g(3) * t212;
t341 = t2 * t218;
t340 = t214 * t3;
t338 = mrSges(5,3) * t171;
t337 = Ifges(4,4) * t216;
t336 = Ifges(4,4) * t220;
t335 = Ifges(6,4) * t214;
t334 = Ifges(6,4) * t218;
t331 = t172 * Ifges(5,4);
t319 = t171 * t214;
t318 = t171 * t218;
t317 = t175 * t214;
t316 = t175 * t218;
t310 = t212 * t221;
t299 = qJD(2) * t217;
t297 = qJD(3) * t220;
t294 = qJD(5) * t214;
t292 = t216 * qJD(2);
t291 = t220 * qJD(2);
t287 = Ifges(6,5) * t38 + Ifges(6,6) * t39 + Ifges(6,3) * t81;
t286 = pkin(3) * t292;
t280 = mrSges(4,3) * t292;
t279 = mrSges(4,3) * t291;
t277 = t214 * t310;
t275 = t218 * t310;
t129 = Ifges(6,4) * t133;
t50 = t134 * Ifges(6,1) + t162 * Ifges(6,5) + t129;
t274 = t50 * t345;
t269 = t212 * t299;
t268 = t221 * t300;
t265 = -t294 / 0.2e1;
t263 = t117 * pkin(4) + pkin(9) * t118;
t259 = t290 / 0.2e1;
t256 = t384 * pkin(3);
t112 = pkin(4) * t172 + pkin(9) * t171;
t251 = Ifges(6,1) * t218 - t335;
t250 = t337 + t375;
t249 = -Ifges(6,2) * t214 + t334;
t248 = Ifges(4,5) * t220 - Ifges(4,6) * t216;
t247 = Ifges(6,5) * t218 - Ifges(6,6) * t214;
t168 = t213 * t216 + t217 * t311;
t245 = t219 * t167 - t168 * t215;
t96 = t167 * t215 + t168 * t219;
t242 = t167 * pkin(3);
t75 = -t214 * t96 - t275;
t238 = -t218 * t96 + t277;
t236 = -t218 * t123 + t175 * t294;
t183 = -qJD(2) * pkin(2) - t271;
t235 = t183 * (mrSges(4,1) * t216 + mrSges(4,2) * t220);
t234 = t216 * (Ifges(4,1) * t220 - t337);
t229 = -t164 * t216 - t220 * t262;
t228 = t229 * pkin(3);
t227 = -t340 + (-t23 * t214 - t22 * t218) * qJD(5);
t226 = t227 + t341;
t10 = -pkin(4) * t208 - t13;
t160 = Ifges(5,4) * t171;
t48 = t134 * Ifges(6,5) + t133 * Ifges(6,6) + t162 * Ifges(6,3);
t7 = t38 * Ifges(6,4) + t39 * Ifges(6,2) + t81 * Ifges(6,6);
t92 = -t171 * Ifges(5,2) + t209 * Ifges(5,6) + t331;
t93 = t172 * Ifges(5,1) + t209 * Ifges(5,5) - t160;
t224 = -t161 * (mrSges(5,1) * t172 - mrSges(5,2) * t171) + (Ifges(6,3) * t172 - t171 * t247) * t349 + (Ifges(6,5) * t172 - t171 * t251) * t351 + (Ifges(6,6) * t172 - t171 * t249) * t352 - t209 * (-Ifges(5,5) * t171 - Ifges(5,6) * t172) / 0.2e1 - t22 * (mrSges(6,1) * t172 + mrSges(6,3) * t318) - t23 * (-mrSges(6,2) * t172 + mrSges(6,3) * t319) + (t133 * t249 + t134 * t251 + t162 * t247) * qJD(5) / 0.2e1 - (-Ifges(5,1) * t171 - t331 + t48) * t172 / 0.2e1 + (-Ifges(5,2) * t172 - t160 + t93) * t171 / 0.2e1 + (t374 + t274) * qJD(5) + t171 * t374 + (Ifges(6,2) * t218 + t335) * t355 + (Ifges(6,1) * t214 + t334) * t356 + t214 * t359 + t49 * t265 + t50 * t318 / 0.2e1 + t319 * t382 + Ifges(5,5) * t82 + Ifges(5,6) * t83 - t57 * t338 + mrSges(6,3) * t341 + t7 * t345 + t92 * t346 + (Ifges(6,5) * t214 + Ifges(6,6) * t218) * t353 - t12 * mrSges(5,2) + t13 * mrSges(5,1) + Ifges(5,3) * t208 - t10 * t383;
t204 = Ifges(4,4) * t291;
t202 = -pkin(4) - t343;
t187 = -qJD(3) * mrSges(4,2) + t279;
t186 = qJD(3) * mrSges(4,1) - t280;
t177 = t255 * qJD(2);
t170 = Ifges(4,1) * t292 + Ifges(4,5) * qJD(3) + t204;
t169 = Ifges(4,6) * qJD(3) + qJD(2) * t250;
t165 = t213 * t313 + t261;
t163 = -t213 * t260 + t314;
t159 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t181;
t158 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t180;
t148 = t151 * pkin(4);
t141 = -mrSges(5,2) * t209 - t338;
t139 = -t182 * t216 + t198;
t130 = -mrSges(4,1) * t180 + mrSges(4,2) * t181;
t126 = qJD(3) * t167 + t220 * t268;
t125 = -qJD(3) * t168 - t216 * t268;
t110 = t115 * pkin(4);
t109 = mrSges(5,1) * t171 + mrSges(5,2) * t172;
t91 = t112 + t286;
t85 = mrSges(6,1) * t162 - mrSges(6,3) * t134;
t84 = -mrSges(6,2) * t162 + mrSges(6,3) * t133;
t72 = -mrSges(5,2) * t208 + mrSges(5,3) * t83;
t62 = t127 * t219 - t321;
t61 = t127 * t215 + t308;
t30 = -mrSges(5,1) * t83 + mrSges(5,2) * t82;
t29 = qJD(4) * t96 - t219 * t125 + t126 * t215;
t28 = qJD(4) * t245 + t125 * t215 + t126 * t219;
t27 = t112 * t214 + t218 * t57;
t26 = t112 * t218 - t214 * t57;
t25 = t214 * t91 + t218 * t62;
t24 = -t214 * t62 + t218 * t91;
t18 = qJD(5) * t238 - t214 * t28 + t218 * t269;
t17 = qJD(5) * t75 + t214 * t269 + t218 * t28;
t1 = [m(2) * qJDD(1) + t125 * t186 + t126 * t187 + t28 * t141 + t168 * t158 + t167 * t159 + t17 * t84 + t18 * t85 + t75 * t19 - t238 * t20 + t96 * t72 - t376 * t245 + t324 * t29 + (-m(2) - m(3) - m(4) - t381) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t223 - t130 - t30) * t221 + (-mrSges(3,1) * t223 - mrSges(3,2) * qJDD(2) + (t109 + t177) * qJD(2)) * t217) * t212 + m(6) * (-t10 * t245 + t17 * t23 + t18 * t22 - t2 * t238 + t29 * t53 + t3 * t75) + m(5) * (t12 * t96 + t13 * t245 + t28 * t58 - t29 * t57 + (-t113 * t221 + t161 * t299) * t212) + m(4) * (t125 * t139 + t126 * t140 + t167 * t67 + t168 * t66 + (-t149 * t221 + t183 * t299) * t212) + m(3) * (qJDD(1) * t213 ^ 2 + (t156 * t221 + t157 * t217) * t212); (-pkin(2) * t149 - (t183 * t217 + (-t139 * t216 + t140 * t220) * t221) * t302) * m(4) + (-t2 * t317 + t22 * t236 - t23 * t237 - t3 * t316) * mrSges(6,3) + (t235 + t248 * qJD(3) / 0.2e1) * qJD(3) - t376 * t243 + (-t10 * t243 + t2 * t60 + t22 * t377 + t23 * t378 + t3 * t59 + t370 * t53) * m(6) + (-t113 * t203 + t12 * t136 + t13 * t243 + t161 * t366 - t370 * t57 + t371 * t58) * m(5) + (t217 * t342 - t157 + t193) * mrSges(3,2) + (-t221 * t342 + t156 + t192) * mrSges(3,1) + (-t57 * t123 - t58 * t124) * mrSges(5,3) + (m(4) * ((-t139 * t220 - t140 * t216) * qJD(3) + t364) - t186 * t297 - t187 * t298 - t216 * t159 + t220 * t158) * pkin(7) + (-t139 * t297 - t140 * t298 + t364) * mrSges(4,3) + t366 * t109 + (Ifges(6,5) * t356 + t287 / 0.2e1 - t12 * mrSges(5,3) + Ifges(6,3) * t353 + Ifges(6,6) * t355 - Ifges(5,4) * t82 - Ifges(5,2) * t83 + t113 * mrSges(5,1) - Ifges(5,6) * t208 + t360) * t174 + (t113 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,1) * t82 + Ifges(5,4) * t83 + Ifges(5,5) * t208 + t10 * t252 + t247 * t353 + t249 * t355 + t251 * t356 + t265 * t50) * t175 + t316 * t359 + ((t385 * t221 + (t222 * t381 + t373) * t217) * t212 - t381 * t203 * t310) * g(3) + t377 * t85 + t378 * t84 + t371 * t141 + t123 * t274 + t234 * t259 + t370 * t324 - t7 * t317 / 0.2e1 + t170 * t297 / 0.2e1 - t169 * t298 / 0.2e1 + Ifges(3,3) * qJDD(2) + t336 * t379 + t250 * t380 + t237 * t382 + t22 * mrSges(6,1) * t124 - t23 * mrSges(6,2) * t124 + t59 * t19 + t60 * t20 + t53 * (mrSges(6,1) * t237 - mrSges(6,2) * t236) + t162 * (-Ifges(6,5) * t236 - Ifges(6,6) * t237 + Ifges(6,3) * t124) / 0.2e1 + t133 * (-Ifges(6,4) * t236 - Ifges(6,2) * t237 + Ifges(6,6) * t124) / 0.2e1 + t123 * t93 / 0.2e1 + t124 * t48 / 0.2e1 - t124 * t92 / 0.2e1 - pkin(2) * t130 + t136 * t72 + t161 * (mrSges(5,1) * t124 + mrSges(5,2) * t123) + t149 * t255 + (Ifges(5,1) * t123 - Ifges(5,4) * t124) * t346 + (-Ifges(6,1) * t236 - Ifges(6,4) * t237 + Ifges(6,5) * t124) * t350 - t171 * (Ifges(5,4) * t123 - Ifges(5,2) * t124) / 0.2e1 - (t217 * t282 + t221 * t231) * t342 - t203 * t30 + t209 * (Ifges(5,5) * t123 - Ifges(5,6) * t124) / 0.2e1 + qJDD(3) * (Ifges(4,5) * t216 + Ifges(4,6) * t220) + (Ifges(4,4) * t379 + Ifges(4,2) * t380 - t187 * t271 + t336 * t259) * t220 + (Ifges(4,1) * t181 + Ifges(4,4) * t380 + t186 * t271 - t259 * t375) * t216 + (-t381 * (-t165 * t203 - t166 * t222) + t362 * t166 + t361 * t165) * g(1) + (-t381 * (-t163 * t203 - t164 * t222) + t362 * t164 + t361 * t163) * g(2) - t177 * t272; (-t214 * t85 + t218 * t84 + t141) * pkin(3) * t295 - (-Ifges(4,2) * t292 + t170 + t204) * t291 / 0.2e1 + ((t12 * t215 + t13 * t219 + (-t215 * t57 + t219 * t58) * qJD(4)) * pkin(3) - t161 * t286 + t57 * t61 - t58 * t62) * m(5) + t224 + (t280 + t186) * t140 + (t279 - t187) * t139 + (-t22 * t293 - t23 * t294 - t340) * mrSges(6,3) + (m(6) * t226 - t293 * t85 - t294 * t84 + t365) * (pkin(9) + t344) + (-m(5) * t228 - m(6) * (t116 * pkin(9) + t110 + t228) - t229 * mrSges(4,1) - (-t164 * t220 + t216 * t262) * mrSges(4,2) + t367) * g(2) + (-m(5) * t242 - m(6) * (pkin(9) * t152 + t148 + t242) - mrSges(4,1) * t167 + mrSges(4,2) * t168 + t369) * g(3) + (t10 * t202 + (t215 * t53 + (-t22 * t214 + t23 * t218) * t219) * qJD(4) * pkin(3) - t22 * t24 - t23 * t25 - t53 * t61) * m(6) + (-t384 * mrSges(4,1) - (-t166 * t220 - t216 * t315) * mrSges(4,2) - m(5) * t256 - m(6) * (t256 + t263) + t368) * g(1) + t324 * (pkin(3) * t296 - t61) - t248 * t290 / 0.2e1 + t169 * t292 / 0.2e1 + Ifges(4,3) * qJDD(3) - t109 * t286 + t58 * t332 - t66 * mrSges(4,2) + t67 * mrSges(4,1) - t25 * t84 - t24 * t85 - qJD(2) * t235 - t62 * t141 + t71 * t343 + t72 * t344 + Ifges(4,6) * t180 + Ifges(4,5) * t181 + t202 * t16 - t223 * t234 / 0.2e1; t224 + (-m(6) * t110 + t367) * g(2) + (-m(6) * t263 + t368) * g(1) + (-m(6) * t148 + t369) * g(3) + (-t324 + t332) * t58 - m(6) * (t22 * t26 + t23 * t27 + t53 * t58) + t227 * mrSges(6,3) - t10 * t358 + t226 * t357 - t27 * t84 - t26 * t85 - t57 * t141 - pkin(4) * t16 + ((-t214 * t84 - t218 * t85) * qJD(5) + (-g(2) * t116 - g(3) * t152) * m(6) + t365) * pkin(9); -t53 * (mrSges(6,1) * t134 + mrSges(6,2) * t133) + (Ifges(6,1) * t133 - t333) * t351 + t49 * t350 + (Ifges(6,5) * t133 - Ifges(6,6) * t134) * t349 - t22 * t84 + t23 * t85 - g(1) * ((-t118 * t214 + t165 * t218) * mrSges(6,1) + (-t118 * t218 - t165 * t214) * mrSges(6,2)) - g(2) * ((-t116 * t214 + t163 * t218) * mrSges(6,1) + (-t116 * t218 - t163 * t214) * mrSges(6,2)) - g(3) * ((-t152 * t214 - t275) * mrSges(6,1) + (-t152 * t218 + t277) * mrSges(6,2)) + (t133 * t22 + t134 * t23) * mrSges(6,3) + t287 + (-Ifges(6,2) * t134 + t129 + t50) * t352 + t360;];
tau = t1;
