% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:23
% EndTime: 2019-03-08 23:31:40
% DurationCPUTime: 10.00s
% Computational Cost: add. (6024->643), mult. (14780->853), div. (0->0), fcn. (9990->10), ass. (0->312)
t217 = cos(qJ(3));
t290 = qJD(2) * t217;
t206 = Ifges(4,4) * t290;
t412 = -t206 / 0.2e1;
t408 = Ifges(5,1) + Ifges(6,1);
t400 = Ifges(6,4) + Ifges(5,5);
t212 = sin(qJ(4));
t287 = qJD(4) * t212;
t360 = pkin(9) - pkin(10);
t213 = sin(qJ(3));
t292 = qJD(2) * t213;
t214 = sin(qJ(2));
t209 = sin(pkin(6));
t294 = qJD(1) * t209;
t272 = t214 * t294;
t179 = qJD(2) * pkin(8) + t272;
t210 = cos(pkin(6));
t304 = t210 * t217;
t142 = qJD(1) * t304 - t213 * t179;
t258 = pkin(3) * t213 - pkin(9) * t217;
t175 = t258 * qJD(2);
t216 = cos(qJ(4));
t79 = t216 * t142 + t212 * t175;
t69 = qJ(5) * t292 + t79;
t411 = pkin(10) * t212 * t290 + t360 * t287 + t69;
t191 = t360 * t216;
t302 = t216 * t217;
t279 = pkin(10) * t302;
t361 = pkin(4) + pkin(5);
t78 = -t212 * t142 + t175 * t216;
t410 = qJD(4) * t191 - (-t213 * t361 - t279) * qJD(2) + t78;
t293 = qJD(1) * t213;
t270 = t210 * t293;
t143 = t179 * t217 + t270;
t120 = qJD(3) * pkin(9) + t143;
t187 = -pkin(3) * t217 - pkin(9) * t213 - pkin(2);
t218 = cos(qJ(2));
t271 = t218 * t294;
t144 = qJD(2) * t187 - t271;
t61 = -t212 * t120 + t216 * t144;
t378 = qJD(5) - t61;
t409 = qJD(3) / 0.2e1;
t407 = Ifges(5,6) - Ifges(6,6);
t282 = t216 * qJD(3);
t170 = t212 * t292 - t282;
t171 = qJD(3) * t212 + t216 * t292;
t211 = sin(qJ(6));
t215 = cos(qJ(6));
t236 = t170 * t211 + t171 * t215;
t98 = t170 * t215 - t171 * t211;
t92 = Ifges(7,4) * t98;
t406 = Ifges(7,2) * t236 - t92;
t405 = -pkin(10) * t171 + t378;
t323 = qJD(2) * pkin(2);
t180 = -t271 - t323;
t265 = Ifges(4,5) * t409;
t62 = t216 * t120 + t212 * t144;
t237 = t212 * t62 + t216 * t61;
t200 = qJD(4) - t290;
t45 = -pkin(4) * t200 + t378;
t193 = t200 * qJ(5);
t48 = t193 + t62;
t238 = t212 * t48 - t216 * t45;
t326 = Ifges(6,5) * t216;
t244 = Ifges(6,3) * t212 + t326;
t332 = Ifges(5,4) * t216;
t248 = -Ifges(5,2) * t212 + t332;
t253 = mrSges(6,1) * t212 - mrSges(6,3) * t216;
t255 = mrSges(5,1) * t212 + mrSges(5,2) * t216;
t259 = qJD(3) * pkin(3) + t142;
t324 = Ifges(6,6) * t212;
t325 = Ifges(5,6) * t212;
t329 = Ifges(5,5) * t216;
t331 = Ifges(6,4) * t216;
t343 = t216 / 0.2e1;
t345 = t212 / 0.2e1;
t346 = -t212 / 0.2e1;
t351 = t171 / 0.2e1;
t353 = t170 / 0.2e1;
t354 = -t170 / 0.2e1;
t168 = Ifges(5,4) * t170;
t328 = Ifges(6,5) * t170;
t385 = t171 * t408 + t200 * t400 - t168 + t328;
t388 = t200 / 0.2e1;
t327 = Ifges(6,5) * t212;
t333 = Ifges(5,4) * t212;
t394 = t408 * t216 + t327 - t333;
t226 = qJ(5) * t171 + t259;
t65 = pkin(4) * t170 - t226;
t167 = Ifges(6,5) * t171;
t83 = t200 * Ifges(6,6) + t170 * Ifges(6,3) + t167;
t334 = Ifges(5,4) * t171;
t86 = -t170 * Ifges(5,2) + t200 * Ifges(5,6) + t334;
t369 = t238 * mrSges(6,2) + t237 * mrSges(5,3) + t259 * t255 - t244 * t353 - t248 * t354 - t65 * t253 - t83 * t345 - t86 * t346 - t394 * t351 - (-t325 + t329 + t324 + t331) * t388 - t385 * t343;
t402 = t292 / 0.2e1;
t404 = -t180 * mrSges(4,2) + t142 * mrSges(4,3) - Ifges(4,1) * t402 - t265 + t369 + t412;
t30 = -t200 * t361 + t405;
t42 = pkin(10) * t170 + t62;
t34 = t193 + t42;
t10 = -t211 * t34 + t215 * t30;
t11 = t211 * t30 + t215 * t34;
t281 = qJD(2) * qJD(3);
t263 = t213 * t281;
t194 = qJD(6) - t200;
t350 = -t194 / 0.2e1;
t266 = t213 * t287;
t280 = qJD(3) * qJD(4);
t133 = t216 * t280 + (t217 * t282 - t266) * qJD(2);
t286 = qJD(4) * t216;
t288 = qJD(3) * t217;
t134 = t212 * t280 + (t212 * t288 + t213 * t286) * qJD(2);
t28 = qJD(6) * t98 + t133 * t215 + t134 * t211;
t29 = -qJD(6) * t236 - t133 * t211 + t134 * t215;
t386 = Ifges(7,5) * t28 + Ifges(7,6) * t29;
t44 = -t170 * t361 + t226;
t403 = (t10 * t98 + t11 * t236) * mrSges(7,3) - Ifges(7,3) * t263 + (Ifges(7,5) * t98 - Ifges(7,6) * t236) * t350 - t44 * (mrSges(7,1) * t236 + mrSges(7,2) * t98) + t386;
t401 = -qJD(3) / 0.2e1;
t399 = t400 * t263 + (-Ifges(5,4) + Ifges(6,5)) * t134 + t408 * t133;
t190 = t360 * t212;
t129 = t190 * t215 - t191 * t211;
t398 = qJD(6) * t129 + t410 * t211 - t411 * t215;
t130 = t190 * t211 + t191 * t215;
t397 = -qJD(6) * t130 + t411 * t211 + t410 * t215;
t307 = qJ(5) * t216;
t232 = -t212 * t361 + t307;
t284 = qJD(5) * t212;
t396 = qJD(4) * t232 + t284 + t270 - (qJD(2) * t232 - t179) * t217;
t395 = t408 * t212 - t326 + t332;
t330 = Ifges(7,4) * t236;
t392 = Ifges(7,1) * t98 - t330;
t365 = t28 / 0.2e1;
t364 = t29 / 0.2e1;
t33 = Ifges(7,1) * t236 + Ifges(7,5) * t194 + t92;
t389 = t33 / 0.2e1;
t359 = -t236 / 0.2e1;
t387 = -t263 / 0.2e1;
t264 = Ifges(4,6) * t401;
t181 = -qJ(5) * t211 - t215 * t361;
t380 = qJD(6) * t181 - t211 * t42 + t215 * t405;
t182 = qJ(5) * t215 - t211 * t361;
t379 = -qJD(6) * t182 - t211 * t405 - t215 * t42;
t235 = t211 * t216 - t212 * t215;
t153 = t235 * t213;
t377 = t212 * t400 + t216 * t407;
t376 = -qJ(5) * t133 - qJD(5) * t171;
t178 = t258 * qJD(3);
t137 = (t178 + t272) * qJD(2);
t305 = t209 * t218;
t268 = qJD(2) * t305;
t231 = qJD(3) * t210 + t268;
t225 = qJD(1) * t231;
t289 = qJD(3) * t213;
t89 = -t179 * t289 + t217 * t225;
t19 = -t120 * t287 + t212 * t137 + t144 * t286 + t216 * t89;
t20 = -t120 * t286 + t137 * t216 - t144 * t287 - t212 * t89;
t375 = t19 * t216 - t20 * t212;
t16 = qJ(5) * t263 + t200 * qJD(5) + t19;
t17 = -pkin(4) * t263 - t20;
t374 = t16 * t216 + t17 * t212;
t373 = qJD(4) - qJD(6);
t308 = qJ(5) * t212;
t370 = -t216 * t361 - t308;
t276 = Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t277 = -Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t278 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t368 = -t277 * t170 + t278 * t171 + t276 * t200 + t11 * mrSges(7,2) + t180 * mrSges(4,1) + t48 * mrSges(6,3) + t61 * mrSges(5,1) + t264 - (Ifges(4,4) * t213 + t217 * Ifges(4,2)) * qJD(2) / 0.2e1 - t194 * Ifges(7,3) - t236 * Ifges(7,5) - t98 * Ifges(7,6) + Ifges(5,6) * t354 + Ifges(6,6) * t353 - t10 * mrSges(7,1) - t143 * mrSges(4,3) - t45 * mrSges(6,1) - t62 * mrSges(5,2) + (Ifges(5,3) + Ifges(6,2)) * t388 + t400 * t351;
t367 = Ifges(7,4) * t365 + Ifges(7,2) * t364 + Ifges(7,6) * t387;
t366 = Ifges(7,1) * t365 + Ifges(7,4) * t364 + Ifges(7,5) * t387;
t363 = -t98 / 0.2e1;
t362 = t98 / 0.2e1;
t358 = t236 / 0.2e1;
t357 = t133 / 0.2e1;
t356 = -t134 / 0.2e1;
t355 = t134 / 0.2e1;
t352 = -t171 / 0.2e1;
t349 = t194 / 0.2e1;
t348 = -t200 / 0.2e1;
t344 = -t216 / 0.2e1;
t342 = pkin(8) * t212;
t341 = pkin(10) * t213;
t93 = -mrSges(6,2) * t134 + mrSges(6,3) * t263;
t96 = -mrSges(5,2) * t263 - mrSges(5,3) * t134;
t338 = t93 + t96;
t94 = mrSges(5,1) * t263 - mrSges(5,3) * t133;
t95 = -mrSges(6,1) * t263 + t133 * mrSges(6,2);
t337 = -t94 + t95;
t336 = mrSges(5,3) * t170;
t335 = mrSges(5,3) * t171;
t306 = t209 * t214;
t156 = t213 * t306 - t304;
t267 = t179 * t288;
t90 = t213 * t225 + t267;
t320 = t156 * t90;
t106 = mrSges(6,1) * t170 - mrSges(6,3) * t171;
t43 = -mrSges(7,1) * t98 + mrSges(7,2) * t236;
t313 = t106 - t43;
t309 = qJ(5) * t170;
t303 = t212 * t217;
t301 = t216 * t218;
t234 = t211 * t212 + t215 * t216;
t103 = t373 * t234;
t229 = t217 * t234;
t149 = qJD(2) * t229;
t300 = t103 - t149;
t104 = t373 * t235;
t228 = t235 * t217;
t148 = qJD(2) * t228;
t299 = -t104 + t148;
t138 = -mrSges(5,2) * t200 - t336;
t141 = -mrSges(6,2) * t170 + mrSges(6,3) * t200;
t298 = t138 + t141;
t139 = mrSges(5,1) * t200 - t335;
t140 = -mrSges(6,1) * t200 + mrSges(6,2) * t171;
t297 = t139 - t140;
t296 = t212 * t178 + t187 * t286;
t295 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t170 - mrSges(5,2) * t171 - mrSges(4,3) * t292;
t202 = pkin(8) * t302;
t147 = t212 * t187 + t202;
t291 = qJD(2) * t214;
t283 = qJD(5) * t216;
t275 = mrSges(4,3) * t290;
t274 = t212 * t305;
t273 = -pkin(4) - t342;
t269 = t209 * t291;
t189 = -qJD(3) * mrSges(4,2) + t275;
t262 = -m(4) * t143 - t189;
t201 = pkin(8) * t303;
t146 = t187 * t216 - t201;
t12 = pkin(10) * t134 + t16;
t13 = -pkin(10) * t133 - t263 * t361 - t20;
t1 = qJD(6) * t10 + t12 * t215 + t13 * t211;
t2 = -qJD(6) * t11 - t12 * t211 + t13 * t215;
t260 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t131 = -qJ(5) * t217 + t147;
t257 = qJD(4) * t202 - t178 * t216 + t187 * t287;
t256 = mrSges(5,1) * t216 - mrSges(5,2) * t212;
t254 = mrSges(6,1) * t216 + mrSges(6,3) * t212;
t247 = Ifges(5,2) * t216 + t333;
t243 = -Ifges(6,3) * t216 + t327;
t242 = pkin(4) * t216 + t308;
t241 = pkin(4) * t212 - t307;
t75 = -mrSges(7,2) * t194 + mrSges(7,3) * t98;
t76 = mrSges(7,1) * t194 - mrSges(7,3) * t236;
t239 = -t211 * t76 + t215 * t75;
t101 = t212 * t341 + t131;
t208 = t217 * pkin(4);
t91 = pkin(5) * t217 + t201 + t208 + (-t187 - t341) * t216;
t36 = t101 * t215 + t211 * t91;
t35 = -t101 * t211 + t215 * t91;
t157 = t210 * t213 + t217 * t306;
t112 = t157 * t212 + t209 * t301;
t113 = t157 * t216 - t274;
t46 = t112 * t215 - t113 * t211;
t47 = t112 * t211 + t113 * t215;
t233 = pkin(8) + t241;
t230 = -m(4) * t142 - m(5) * t259 - t295;
t227 = -pkin(8) + t232;
t73 = (-t213 * t282 - t217 * t287) * pkin(8) + t296;
t224 = -t20 * mrSges(5,1) + t17 * mrSges(6,1) + t19 * mrSges(5,2) - t16 * mrSges(6,3) + t260;
t220 = qJD(2) ^ 2;
t205 = qJ(5) * t289;
t199 = Ifges(6,2) * t263;
t198 = Ifges(5,3) * t263;
t183 = -pkin(3) - t242;
t174 = (-mrSges(4,1) * t217 + mrSges(4,2) * t213) * qJD(2);
t166 = pkin(3) - t370;
t162 = (mrSges(4,1) * t213 + mrSges(4,2) * t217) * t281;
t155 = qJD(4) * t241 - t284;
t154 = t234 * t213;
t150 = t233 * t213;
t136 = (t212 * t214 + t217 * t301) * t294;
t135 = -t216 * t272 + t271 * t303;
t132 = -t146 + t208;
t128 = t227 * t213;
t125 = Ifges(6,4) * t133;
t124 = Ifges(5,5) * t133;
t123 = Ifges(5,6) * t134;
t122 = Ifges(6,6) * t134;
t111 = -qJD(3) * t156 + t217 * t268;
t110 = qJD(3) * t157 + t213 * t268;
t105 = pkin(4) * t171 + t309;
t82 = t270 + (qJD(2) * t241 + t179) * t217;
t77 = -t171 * t361 - t309;
t74 = t289 * t342 - t257;
t72 = (qJD(4) * t242 - t283) * t213 + t233 * t288;
t71 = -pkin(4) * t292 - t78;
t68 = t273 * t289 + t257;
t67 = mrSges(5,1) * t134 + mrSges(5,2) * t133;
t66 = mrSges(6,1) * t134 - mrSges(6,3) * t133;
t64 = t135 * t211 + t136 * t215;
t63 = t135 * t215 - t136 * t211;
t59 = -qJD(5) * t217 + t205 + t73;
t54 = t133 * Ifges(5,4) - t134 * Ifges(5,2) + Ifges(5,6) * t263;
t53 = t133 * Ifges(6,5) + Ifges(6,6) * t263 + t134 * Ifges(6,3);
t52 = -qJD(3) * t228 + t103 * t213;
t51 = qJD(3) * t229 + t153 * t373;
t49 = (qJD(4) * t370 + t283) * t213 + t227 * t288;
t40 = t205 + (-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t216 * t213 + (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t212) * t217 + t296;
t39 = -qJD(4) * t112 + t111 * t216 + t212 * t269;
t38 = -qJD(4) * t274 + t111 * t212 + t157 * t286 - t216 * t269;
t37 = pkin(10) * t266 + (-t279 + (-pkin(5) + t273) * t213) * qJD(3) + t257;
t32 = Ifges(7,2) * t98 + Ifges(7,6) * t194 + t330;
t25 = mrSges(7,2) * t263 + mrSges(7,3) * t29;
t24 = -mrSges(7,1) * t263 - mrSges(7,3) * t28;
t23 = pkin(4) * t134 + t376 + t90;
t18 = -t134 * t361 - t231 * t293 - t267 - t376;
t9 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t6 = qJD(6) * t46 + t211 * t38 + t215 * t39;
t5 = -qJD(6) * t47 - t211 * t39 + t215 * t38;
t4 = -qJD(6) * t36 - t211 * t40 + t215 * t37;
t3 = qJD(6) * t35 + t211 * t37 + t215 * t40;
t7 = [-t157 * mrSges(4,3) * t263 + t111 * t189 + t46 * t24 + t47 * t25 + t5 * t76 + t6 * t75 + t298 * t39 - t297 * t38 + t338 * t113 + t337 * t112 + ((-mrSges(3,2) * t220 - t162) * t218 + (-mrSges(3,1) * t220 + qJD(2) * t174) * t214) * t209 + (qJD(3) * t275 + t66 + t67 - t9) * t156 + (-t295 + t313) * t110 + m(6) * (t110 * t65 + t112 * t17 + t113 * t16 + t156 * t23 + t38 * t45 + t39 * t48) + m(5) * (-t110 * t259 - t112 * t20 + t113 * t19 - t38 * t61 + t39 * t62 + t320) + m(7) * (t1 * t47 + t10 * t5 + t11 * t6 - t110 * t44 - t156 * t18 + t2 * t46) + m(4) * (-t110 * t142 + t111 * t143 + t320 + t157 * t89 + (t180 - t271) * t269); m(7) * (t1 * t36 + t10 * t4 + t11 * t3 + t128 * t18 + t2 * t35 + t44 * t49) + m(6) * (t131 * t16 + t132 * t17 + t150 * t23 + t45 * t68 + t48 * t59 + t65 * t72) + (-t174 + 0.2e1 * (-t180 / 0.2e1 - t323 / 0.2e1) * m(4)) * t272 - m(7) * (t10 * t63 + t11 * t64) - m(6) * (t135 * t45 + t136 * t48) - m(5) * (-t135 * t61 + t136 * t62) + (-t63 + t4) * t76 + (-t64 + t3) * t75 + (t244 * t355 + t248 * t356 + t23 * t253 + t53 * t345 + t54 * t346 + (-t19 * t212 - t20 * t216) * mrSges(5,3) + (-t16 * t212 + t17 * t216) * mrSges(6,2) + (mrSges(4,3) + t255) * t90 + (mrSges(4,2) * t291 + (-m(6) * t65 + m(7) * t44 - t230 - t313) * t218) * t294 + (t86 * t344 + t247 * t353 + t243 * t354 - t259 * t256 + t65 * t254 + (t212 * t61 - t216 * t62) * mrSges(5,3) + (-t212 * t45 - t216 * t48) * mrSges(6,2) + t395 * t352 + t377 * t348 + t385 * t346) * qJD(4) + ((-Ifges(7,5) * t154 / 0.2e1 + Ifges(7,6) * t153 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t329 / 0.2e1 - t325 / 0.2e1 + t331 / 0.2e1 + t324 / 0.2e1) * t213 + (0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(7,3) - 0.3e1 / 0.2e1 * Ifges(4,2) - t276) * t217) * qJD(2) + t264 + t368) * qJD(3) + t394 * t357 + (qJD(4) * t83 + t399) * t343 + (t67 + (m(4) + m(5)) * t90 + t262 * qJD(3)) * pkin(8)) * t213 + (-t1 * t153 - t10 * t51 + t11 * t52 - t154 * t2) * mrSges(7,3) + t18 * (mrSges(7,1) * t153 + mrSges(7,2) * t154) + (Ifges(7,4) * t154 - Ifges(7,2) * t153) * t364 + (Ifges(7,1) * t154 - Ifges(7,4) * t153) * t365 + m(5) * (t146 * t20 + t147 * t19 + t61 * t74 + t62 * t73) + (t224 - t198 / 0.2e1 - t199 / 0.2e1 + t123 / 0.2e1 - t124 / 0.2e1 - t125 / 0.2e1 - t122 / 0.2e1 + (-mrSges(4,1) * t291 + t218 * t262) * t294 + t277 * t134 + (m(4) * pkin(8) + mrSges(4,3)) * t89 - t278 * t133 + (t230 * pkin(8) + t265 + 0.3e1 / 0.2e1 * t206 - t404) * qJD(3) + t386) * t217 + t297 * t135 - t298 * t136 - pkin(2) * t162 + t150 * t66 + t146 * t94 + t147 * t96 + t73 * t138 + t74 * t139 + t68 * t140 + t59 * t141 + t128 * t9 + t131 * t93 + t132 * t95 + t72 * t106 + t44 * (-mrSges(7,1) * t52 + mrSges(7,2) * t51) + t52 * t32 / 0.2e1 + t49 * t43 + t35 * t24 + t36 * t25 + (Ifges(7,5) * t51 + Ifges(7,6) * t52) * t349 + (Ifges(7,1) * t51 + Ifges(7,4) * t52) * t358 + (Ifges(7,4) * t51 + Ifges(7,2) * t52) * t362 + t154 * t366 - t153 * t367 + t51 * t389; ((t264 + (-Ifges(7,5) * t235 - Ifges(7,6) * t234) * t401 + Ifges(4,4) * t402 + t377 * t409 - t368) * t213 + (t412 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t292 + t265 + t404) * t217) * qJD(2) - t234 * t367 + t18 * (mrSges(7,1) * t234 - mrSges(7,2) * t235) + (-Ifges(7,4) * t235 - Ifges(7,2) * t234) * t364 + (-Ifges(7,1) * t235 - Ifges(7,4) * t234) * t365 + (-t1 * t234 - t10 * t300 + t11 * t299 + t2 * t235) * mrSges(7,3) - t235 * t366 + (Ifges(7,5) * t103 - Ifges(7,6) * t104) * t349 + (Ifges(7,1) * t103 - Ifges(7,4) * t104) * t358 + (Ifges(7,4) * t103 - Ifges(7,2) * t104) * t362 + (t148 / 0.2e1 - t104 / 0.2e1) * t32 + (Ifges(7,5) * t149 - Ifges(7,6) * t148) * t350 + (Ifges(7,1) * t149 - Ifges(7,4) * t148) * t359 + (Ifges(7,4) * t149 - Ifges(7,2) * t148) * t363 - m(5) * (-t143 * t259 + t61 * t78 + t62 * t79) + (-t369 + (-m(5) * t237 - m(6) * t238 - t212 * t298 - t216 * t297) * pkin(9)) * qJD(4) - m(6) * (t45 * t71 + t48 * t69 + t65 * t82) + (-t149 / 0.2e1 + t103 / 0.2e1) * t33 + t395 * t357 + t397 * t76 + t398 * t75 + (t1 * t130 + t10 * t397 + t11 * t398 + t129 * t2 + t166 * t18 + t396 * t44) * m(7) + t399 * t345 + (t155 - t82) * t106 + t374 * mrSges(6,2) + m(6) * (pkin(9) * t374 + t155 * t65 + t183 * t23) + m(5) * (-pkin(3) * t90 + pkin(9) * t375) + t375 * mrSges(5,3) + (t212 * t337 + t216 * t338) * pkin(9) + (-mrSges(7,1) * t299 + mrSges(7,2) * t300) * t44 + t295 * t143 - t142 * t189 + t183 * t66 + t166 * t9 - t79 * t138 - t78 * t139 - t71 * t140 - t69 * t141 + t129 * t24 + t130 * t25 - t89 * mrSges(4,2) - pkin(3) * t67 + t396 * t43 + t54 * t343 + t53 * t344 + t243 * t355 + t247 * t356 - t23 * t254 + (-mrSges(4,1) - t256) * t90; (-t400 * t170 - t171 * t407) * t348 + (-t170 * t408 + t167 - t334 + t83) * t352 + t98 * t389 + t259 * (mrSges(5,1) * t171 - mrSges(5,2) * t170) - t403 + (-pkin(4) * t17 + qJ(5) * t16 - t105 * t65 + t378 * t48 - t45 * t62) * m(6) + t379 * t76 + (t1 * t182 + t10 * t379 + t11 * t380 + t181 * t2 - t44 * t77) * m(7) + t380 * t75 - t224 + t198 + t199 - t123 + t124 + t125 + t122 + t406 * t363 + (t32 - t392) * t359 + (-t298 - t336) * t61 + (t297 + t335) * t62 + t181 * t24 + t182 * t25 - t65 * (mrSges(6,1) * t171 + mrSges(6,3) * t170) + qJD(5) * t141 - t105 * t106 + qJ(5) * t93 - pkin(4) * t95 - t77 * t43 + (t170 * t45 + t171 * t48) * mrSges(6,2) + t86 * t351 + (Ifges(6,3) * t171 - t328) * t354 + (-Ifges(5,2) * t171 - t168 + t385) * t353; t211 * t25 + t215 * t24 + t313 * t171 + t239 * qJD(6) + (-t141 - t239) * t200 + t95 + (t1 * t211 - t171 * t44 + t2 * t215 + t194 * (-t10 * t211 + t11 * t215)) * m(7) + (t171 * t65 - t200 * t48 + t17) * m(6); t392 * t359 + t32 * t358 - t10 * t75 + t11 * t76 + t260 + (t33 - t406) * t363 + t403;];
tauc  = t7(:);
