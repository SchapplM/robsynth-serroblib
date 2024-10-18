% Calculate vector of inverse dynamics joint torques for
% S5RRRRR13
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
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:24
% EndTime: 2024-09-27 17:30:32
% DurationCPUTime: 4.92s
% Computational Cost: add. (8869->471), mult. (13746->647), div. (0->0), fcn. (8742->22), ass. (0->272)
t259 = qJD(1) + qJD(2);
t272 = cos(qJ(2));
t359 = pkin(1) * qJD(1);
t325 = t272 * t359;
t208 = pkin(2) * t259 + t325;
t271 = cos(qJ(3));
t266 = sin(qJ(3));
t267 = sin(qJ(2));
t326 = t267 * t359;
t304 = t266 * t326;
t156 = t271 * t208 - t304;
t157 = t208 * t266 + t271 * t326;
t270 = cos(qJ(4));
t263 = cos(pkin(5));
t330 = qJD(4) * t270;
t315 = t263 * t330;
t265 = sin(qJ(4));
t343 = t263 * t265;
t401 = pkin(3) * t315 - t156 * t270 + t157 * t343;
t340 = t267 * t271;
t290 = -t266 * t272 - t340;
t179 = t290 * t359;
t341 = t266 * t267;
t289 = t271 * t272 - t341;
t180 = t289 * t359;
t242 = pkin(2) * t271 + pkin(3);
t332 = qJD(3) * t271;
t400 = t179 * t343 - t242 * t315 + (-pkin(2) * t332 + t180) * t270;
t342 = t263 * t270;
t399 = -t179 * t342 + t180 * t265 + (-t265 * t271 - t266 * t342) * qJD(3) * pkin(2);
t372 = m(6) * pkin(4);
t398 = -mrSges(5,1) - t372;
t206 = t242 * t343;
t262 = sin(pkin(5));
t253 = t262 * pkin(9);
t217 = pkin(2) * t266 + t253;
t365 = pkin(10) * t262;
t310 = -t217 - t365;
t397 = (t270 * t310 - t206) * qJD(4) + t399;
t333 = qJD(3) * t266;
t324 = pkin(2) * t333;
t303 = t263 * t324;
t396 = -(qJD(4) * t310 - t303) * t265 + t400;
t323 = t262 * (pkin(9) + pkin(10));
t302 = t265 * t323;
t395 = -qJD(4) * t302 + t401;
t232 = pkin(3) * t343;
t83 = -t156 * t265 - t157 * t342;
t394 = (-t270 * t323 - t232) * qJD(4) - t83;
t264 = sin(qJ(5));
t269 = cos(qJ(5));
t171 = (-t264 * t265 + t269 * t270) * t262;
t247 = qJD(3) + t259;
t138 = t247 * t171;
t257 = qJDD(1) + qJDD(2);
t246 = qJDD(3) + t257;
t331 = qJD(4) * t265;
t151 = (t246 * t270 - t247 * t331) * t262;
t152 = (t246 * t265 + t247 * t330) * t262;
t57 = qJD(5) * t138 + t151 * t264 + t152 * t269;
t393 = Ifges(6,5) * t57;
t292 = t264 * t270 + t265 * t269;
t285 = t292 * qJD(5);
t347 = t247 * t262;
t58 = t151 * t269 - t152 * t264 - t285 * t347;
t392 = Ifges(6,6) * t58;
t391 = -mrSges(5,3) - mrSges(6,3);
t390 = Ifges(5,5) * t152;
t389 = Ifges(5,6) * t151;
t209 = t263 * t246 + qJDD(4);
t388 = Ifges(5,3) * t209;
t200 = qJDD(5) + t209;
t387 = Ifges(6,3) * t200;
t258 = t262 ^ 2;
t346 = t247 * t258;
t386 = t179 + t324;
t128 = pkin(9) * t347 + t157;
t300 = pkin(10) * t347 + t128;
t135 = pkin(3) * t247 + t156;
t322 = t135 * t343;
t60 = t270 * t300 + t322;
t353 = t264 * t60;
t210 = t263 * t247 + qJD(4);
t116 = t135 * t342;
t59 = -t265 * t300 + t116;
t50 = pkin(4) * t210 + t59;
t21 = t269 * t50 - t353;
t367 = pkin(1) * t272;
t196 = -qJD(2) * t326 + qJDD(1) * t367;
t170 = pkin(2) * t257 + t196;
t334 = qJD(2) * t272;
t197 = (qJD(1) * t334 + qJDD(1) * t267) * pkin(1);
t77 = -qJD(3) * t304 + t266 * t170 + t271 * t197 + t208 * t332;
t66 = t246 * t253 + t77;
t73 = t128 * t270 + t322;
t78 = -qJD(3) * t157 + t271 * t170 - t197 * t266;
t74 = pkin(3) * t246 + t78;
t14 = -qJD(4) * t73 - t265 * t66 + t74 * t342;
t8 = pkin(4) * t209 - pkin(10) * t152 + t14;
t13 = -t128 * t331 + t135 * t315 + t270 * t66 + t74 * t343;
t9 = pkin(10) * t151 + t13;
t4 = t21 * qJD(5) + t264 * t8 + t269 * t9;
t352 = t269 * t60;
t22 = t264 * t50 + t352;
t5 = -t22 * qJD(5) - t264 * t9 + t269 * t8;
t385 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t384 = t14 * mrSges(5,1) - t13 * mrSges(5,2);
t207 = t242 * t342;
t254 = t263 * pkin(4);
t115 = t265 * t310 + t207 + t254;
t146 = t270 * t217 + t206;
t344 = t262 * t270;
t230 = pkin(10) * t344;
t129 = t230 + t146;
t64 = t115 * t269 - t129 * t264;
t383 = qJD(5) * t64 + t397 * t264 - t396 * t269;
t65 = t115 * t264 + t129 * t269;
t382 = -qJD(5) * t65 + t396 * t264 + t397 * t269;
t233 = pkin(3) * t342;
t140 = t233 + t254 - t302;
t190 = pkin(9) * t344 + t232;
t162 = t230 + t190;
t88 = t140 * t264 + t162 * t269;
t381 = -qJD(5) * t88 - t395 * t264 + t394 * t269;
t87 = t140 * t269 - t162 * t264;
t380 = qJD(5) * t87 + t394 * t264 + t395 * t269;
t379 = -t190 * qJD(4) - t83;
t316 = t262 * t331;
t378 = -pkin(9) * t316 + t401;
t377 = (-qJD(4) * t217 - t303) * t265 - t400;
t376 = -qJD(4) * t146 + t399;
t375 = t135 * (mrSges(5,1) * t265 + mrSges(5,2) * t270);
t362 = Ifges(5,4) * t265;
t374 = t265 * (Ifges(5,1) * t270 - t362) * t346;
t243 = pkin(2) + t367;
t182 = -pkin(1) * t341 + t271 * t243;
t177 = pkin(3) + t182;
t154 = t177 * t343;
t183 = pkin(1) * t340 + t266 * t243;
t169 = t253 + t183;
t92 = t270 * t169 + t154;
t373 = m(3) * pkin(1);
t172 = t292 * t262;
t139 = t247 * t172;
t369 = t139 / 0.2e1;
t261 = qJ(1) + qJ(2);
t252 = cos(t261);
t238 = pkin(2) * t252;
t366 = pkin(4) * t270;
t363 = mrSges(6,3) * t138;
t361 = Ifges(5,4) * t270;
t360 = Ifges(6,4) * t139;
t358 = pkin(4) * qJD(5);
t357 = t139 * mrSges(6,3);
t255 = qJ(3) + t261;
t239 = sin(t255);
t356 = t239 * mrSges(4,2);
t240 = cos(t255);
t355 = t240 * mrSges(4,2);
t354 = t258 * t74;
t100 = (-t247 * t366 - t135) * t262;
t351 = t100 * t262;
t350 = t135 * t258;
t349 = t239 * t262;
t345 = t262 * t265;
t260 = qJ(4) + qJ(5);
t248 = pkin(5) + t260;
t227 = cos(t248) / 0.2e1;
t312 = pkin(5) - t260;
t237 = cos(t312);
t194 = t237 / 0.2e1 + t227;
t249 = sin(t260);
t130 = -t194 * t240 + t239 * t249;
t299 = sin(t312);
t328 = sin(t248) / 0.2e1;
t193 = t328 - t299 / 0.2e1;
t251 = cos(t260);
t294 = -t193 * t240 - t239 * t251;
t339 = -t130 * mrSges(6,1) + t294 * mrSges(6,2);
t131 = -t194 * t239 - t240 * t249;
t293 = t193 * t239 - t240 * t251;
t338 = t131 * mrSges(6,1) + t293 * mrSges(6,2);
t337 = (t328 + t299 / 0.2e1) * mrSges(6,1) + (t227 - t237 / 0.2e1) * mrSges(6,2);
t335 = t240 * pkin(3) + pkin(9) * t349;
t329 = m(4) + m(5) + m(6);
t327 = t387 + t392 + t393;
t321 = t247 * t345;
t320 = t247 * t344;
t117 = t243 * t332 + (qJD(2) * t289 - t267 * t333) * pkin(1);
t118 = -t243 * t333 + (qJD(2) * t290 - t267 * t332) * pkin(1);
t319 = t270 * t117 + t118 * t343 + t177 * t315;
t318 = t388 + t389 + t390;
t317 = t238 + t335;
t241 = pkin(3) + t366;
t314 = t345 / 0.2e1;
t313 = t330 / 0.2e1;
t311 = -t169 - t365;
t309 = -t117 * t265 + t118 * t342;
t195 = pkin(4) * t343 - t323;
t308 = -t195 * t239 + t240 * t241;
t307 = pkin(4) * t321;
t306 = mrSges(5,3) * t321;
t305 = mrSges(5,3) * t320;
t221 = pkin(4) * t316;
t298 = t238 + t308;
t297 = t311 * t265;
t155 = t177 * t342;
t82 = t155 + t254 + t297;
t89 = t230 + t92;
t34 = -t264 * t89 + t269 * t82;
t35 = t264 * t82 + t269 * t89;
t187 = (-mrSges(5,1) * t270 + mrSges(5,2) * t265) * t262;
t295 = -qJD(4) * t262 * t375 - t74 * t187;
t288 = -t239 * t265 + t240 * t342;
t167 = -t239 * t342 - t240 * t265;
t72 = -t128 * t265 + t116;
t287 = (-t265 * t73 - t270 * t72) * mrSges(5,3);
t286 = (t270 * Ifges(5,2) + t362) * t262;
t168 = -t239 * t343 + t240 * t270;
t284 = -t240 * mrSges(4,1) - t168 * mrSges(5,1) + t293 * mrSges(6,1) - t167 * mrSges(5,2) - t131 * mrSges(6,2) + t391 * t349;
t283 = t210 * t262 * (Ifges(5,5) * t270 - Ifges(5,6) * t265);
t250 = sin(t261);
t280 = -t252 * mrSges(3,1) + t250 * mrSges(3,2) + t284;
t279 = qJD(4) * t287 + t295;
t136 = Ifges(6,4) * t138;
t205 = qJD(5) + t210;
t70 = Ifges(6,2) * t138 + Ifges(6,6) * t205 + t360;
t71 = Ifges(6,1) * t139 + Ifges(6,5) * t205 + t136;
t278 = -t100 * (mrSges(6,1) * t139 + mrSges(6,2) * t138) + t21 * t363 + t70 * t369 + t327 - t139 * (Ifges(6,1) * t138 - t360) / 0.2e1 - t205 * (Ifges(6,5) * t138 - Ifges(6,6) * t139) / 0.2e1 - (-Ifges(6,2) * t139 + t136 + t71) * t138 / 0.2e1 + t385;
t166 = -t239 * t270 - t240 * t343;
t277 = -t166 * mrSges(5,1) - t294 * mrSges(6,1) + t288 * mrSges(5,2) - t130 * mrSges(6,2) + (m(5) * pkin(3) + m(6) * t241 + mrSges(4,1)) * t239 + (m(6) * t195 + (-m(5) * pkin(9) + t391) * t262) * t240;
t276 = mrSges(3,2) * t252 + (pkin(2) * t329 + mrSges(3,1)) * t250 + t277;
t107 = (qJD(4) + qJD(5)) * t171;
t108 = (-qJD(4) * t292 - t285) * t262;
t113 = Ifges(5,6) * t210 + t247 * t286;
t203 = Ifges(5,4) * t320;
t114 = Ifges(5,1) * t321 + Ifges(5,5) * t210 + t203;
t48 = -pkin(4) * t151 - t262 * t74;
t275 = (Ifges(5,1) * t152 + Ifges(5,4) * t151 + Ifges(5,5) * t209) * t314 + t78 * mrSges(4,1) + t107 * t71 / 0.2e1 + t108 * t70 / 0.2e1 + t100 * (-mrSges(6,1) * t108 + mrSges(6,2) * t107) + t138 * (Ifges(6,4) * t107 + Ifges(6,2) * t108) / 0.2e1 + t205 * (Ifges(6,5) * t107 + Ifges(6,6) * t108) / 0.2e1 + Ifges(4,3) * t246 + (-Ifges(5,2) * t265 + t361) * t313 * t346 + t151 * t286 / 0.2e1 - t113 * t316 / 0.2e1 + (Ifges(6,1) * t107 + Ifges(6,4) * t108) * t369 + (Ifges(5,4) * t152 + Ifges(5,2) * t151 + Ifges(5,6) * t209) * t344 / 0.2e1 + (t327 + t318) * t263 / 0.2e1 + (t283 + t374) * qJD(4) / 0.2e1 + (t13 * t344 - t14 * t345) * mrSges(5,3) + (t152 * (t265 * Ifges(5,1) + t361) / 0.2e1 + t209 * (Ifges(5,5) * t265 + Ifges(5,6) * t270) / 0.2e1 + t114 * t313) * t262 + (t387 / 0.2e1 + t392 / 0.2e1 + t393 / 0.2e1 + t389 / 0.2e1 + t390 / 0.2e1 + t388 / 0.2e1 + t384 + t385) * t263 + (-t21 * t107 + t22 * t108) * mrSges(6,3) + (t48 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t57 + Ifges(6,4) * t58 + Ifges(6,5) * t200) * t172 + (-t48 * mrSges(6,1) + t4 * mrSges(6,3) + Ifges(6,4) * t57 + Ifges(6,2) * t58 + Ifges(6,6) * t200) * t171;
t274 = t196 * mrSges(3,1) + Ifges(3,3) * t257 + t275;
t273 = cos(qJ(1));
t268 = sin(qJ(1));
t256 = t273 * pkin(1);
t202 = t241 * t262;
t189 = -pkin(9) * t345 + t233;
t178 = (-t242 - t366) * t262;
t174 = t262 * t324 + t221;
t149 = t247 * t187;
t148 = -mrSges(5,2) * t210 + t305;
t147 = mrSges(5,1) * t210 - t306;
t145 = -t217 * t265 + t207;
t137 = (-t177 - t366) * t262;
t110 = mrSges(5,1) * t209 - mrSges(5,3) * t152;
t109 = -mrSges(5,2) * t209 + mrSges(5,3) * t151;
t104 = mrSges(6,1) * t205 - t357;
t103 = -mrSges(6,2) * t205 + t363;
t93 = -t118 * t262 + t221;
t91 = -t169 * t265 + t155;
t90 = -mrSges(5,1) * t151 + mrSges(5,2) * t152;
t81 = -mrSges(6,1) * t138 + mrSges(6,2) * t139;
t47 = -mrSges(6,2) * t200 + mrSges(6,3) * t58;
t46 = mrSges(6,1) * t200 - mrSges(6,3) * t57;
t33 = -qJD(4) * t92 + t309;
t32 = -t169 * t331 + t319;
t28 = (t270 * t311 - t154) * qJD(4) + t309;
t27 = qJD(4) * t297 + t319;
t25 = t269 * t59 - t353;
t24 = -t264 * t59 - t352;
t23 = -mrSges(6,1) * t58 + mrSges(6,2) * t57;
t7 = -qJD(5) * t35 - t264 * t27 + t269 * t28;
t6 = qJD(5) * t34 + t264 * t28 + t269 * t27;
t1 = [t34 * t46 + t35 * t47 + t93 * t81 + t6 * t103 + t7 * t104 + t92 * t109 + t91 * t110 + t137 * t23 + (-t117 * t247 - t183 * t246 - t77) * mrSges(4,2) + t33 * t147 + t32 * t148 - t197 * mrSges(3,2) + Ifges(2,3) * qJDD(1) + (t118 * t247 + t182 * t246) * mrSges(4,1) + ((-t257 * t267 - t259 * t334) * mrSges(3,2) + (-qJD(2) * t259 * t267 + t257 * t272) * mrSges(3,1)) * pkin(1) + m(5) * (t13 * t92 + t14 * t91 + t32 * t73 + t33 * t72 + (t118 * t135 + t177 * t74) * t258) + (t280 + t356 - m(5) * (t256 + t317) - m(6) * (t256 + t298) - m(4) * (t238 + t256) + t268 * mrSges(2,2) + (-mrSges(2,1) - t373) * t273) * g(2) + m(6) * (t100 * t93 + t137 * t48 + t21 * t7 + t22 * t6 + t34 * t5 + t35 * t4) + m(4) * (t117 * t157 + t118 * t156 + t182 * t78 + t183 * t77) + (t196 * t272 + t197 * t267) * t373 + t274 + (-t118 * t149 - t177 * t90 + t279) * t262 + (t355 + mrSges(2,2) * t273 + (mrSges(2,1) + (m(3) + t329) * pkin(1)) * t268 + t276) * g(1); t64 * t46 + t65 * t47 + t145 * t110 + t146 * t109 + t174 * t81 + t178 * t23 + t259 * mrSges(3,1) * t326 + t377 * t148 + t376 * t147 + t280 * g(2) + t276 * g(1) + t382 * t104 + t383 * t103 + (t386 * t149 + t179 * t81 - t242 * t90 + t279) * t262 + (g(1) * t240 + g(2) * t239 + t180 * t247 - t77 + (-t246 * t266 - t247 * t332) * pkin(2)) * mrSges(4,2) + (t259 * t325 - t197) * mrSges(3,2) + (-t179 * t247 + (t246 * t271 - t247 * t333) * pkin(2)) * mrSges(4,1) + t274 + (-t298 * g(2) + t100 * t174 + t178 * t48 + t179 * t351 + t21 * t382 + t22 * t383 + t4 * t65 + t5 * t64) * m(6) + (-t238 * g(2) - t156 * t179 - t157 * t180 + (-t156 * t333 + t157 * t332 + t266 * t77 + t271 * t78) * pkin(2)) * m(4) + (-t317 * g(2) + t13 * t146 + t14 * t145 + t242 * t354 - t386 * t350 + t376 * t72 + t377 * t73) * m(5); t157 * t247 * mrSges(4,1) + t87 * t46 + t88 * t47 + t189 * t110 + t190 * t109 - t202 * t23 + t378 * t148 + t379 * t147 + t381 * t104 + t380 * t103 + (t284 + t356) * g(2) + (t277 + t355) * g(1) + (-pkin(3) * t90 + (-t149 - t81) * t157 + (pkin(4) * t265 * t81 + t287) * qJD(4) + t295) * t262 + (t156 * t247 - t77) * mrSges(4,2) + t275 + (-t308 * g(2) + t100 * t221 - t157 * t351 - t202 * t48 + t21 * t381 + t22 * t380 + t4 * t88 + t5 * t87) * m(6) + (pkin(3) * t354 - t335 * g(2) + t13 * t190 + t14 * t189 + t157 * t350 + t378 * t73 + t379 * t72) * m(5); t278 + t346 * t375 - t81 * t307 - m(6) * (t100 * t307 + t21 * t24 + t22 * t25) + (t264 * t4 + t269 * t5 + (-t21 * t264 + t22 * t269) * qJD(5)) * t372 + t22 * t357 + t318 + (t147 + t306) * t73 + (-t148 + t305) * t72 - (-Ifges(5,2) * t321 + t114 + t203) * t320 / 0.2e1 + (-t264 * t358 - t24) * t104 + (t269 * t358 - t25) * t103 + (-t344 * t372 + t187 - t337) * g(3) + (-mrSges(5,2) * t166 + t398 * t288 - t339) * g(2) + (t168 * mrSges(5,2) + t398 * t167 - t338) * g(1) + (t264 * t47 + t269 * t46) * pkin(4) + (-t283 / 0.2e1 + t113 * t314 - t374 / 0.2e1) * t247 + t384; t278 - t21 * t103 - g(1) * t338 - g(3) * t337 + (t104 + t357) * t22 - g(2) * t339;];
tau = t1;
