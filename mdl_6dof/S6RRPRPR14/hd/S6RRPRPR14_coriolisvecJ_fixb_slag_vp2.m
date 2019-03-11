% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:35
% EndTime: 2019-03-09 11:33:08
% DurationCPUTime: 16.26s
% Computational Cost: add. (8154->697), mult. (20952->894), div. (0->0), fcn. (14690->8), ass. (0->324)
t235 = sin(qJ(2));
t231 = sin(pkin(6));
t310 = qJD(1) * t231;
t290 = t235 * t310;
t209 = qJD(4) + t290;
t344 = t209 / 0.2e1;
t234 = sin(qJ(4));
t238 = cos(qJ(2));
t289 = t238 * t310;
t198 = t234 * t289;
t232 = cos(pkin(6));
t309 = qJD(1) * t232;
t221 = qJD(2) + t309;
t237 = cos(qJ(4));
t168 = t221 * t237 - t198;
t346 = t168 / 0.2e1;
t347 = -t168 / 0.2e1;
t391 = Ifges(6,4) - Ifges(5,5);
t166 = qJD(6) + t168;
t350 = t166 / 0.2e1;
t167 = t221 * t234 + t237 * t289;
t233 = sin(qJ(6));
t236 = cos(qJ(6));
t104 = t167 * t233 + t209 * t236;
t354 = t104 / 0.2e1;
t103 = t167 * t236 - t209 * t233;
t356 = t103 / 0.2e1;
t425 = Ifges(7,5) * t354 + Ifges(7,6) * t356 + Ifges(7,3) * t350;
t432 = -Ifges(5,1) * t346 + Ifges(6,2) * t347 + t344 * t391 - t425;
t431 = Ifges(6,6) + Ifges(5,4);
t349 = -t167 / 0.2e1;
t301 = pkin(1) * t309;
t219 = t238 * t301;
t362 = pkin(3) + pkin(8);
t281 = -qJ(5) * t234 - t362;
t304 = qJD(4) * t237;
t305 = qJD(4) * t234;
t286 = pkin(4) * t304 + qJ(5) * t305 + qJD(3);
t361 = pkin(4) + pkin(10);
t429 = (pkin(10) * qJD(4) - qJD(5)) * t237 + t286 - t219 - (-t237 * t361 + t281) * t290;
t240 = -pkin(2) - pkin(9);
t335 = pkin(5) - t240;
t283 = qJD(4) * t335;
t214 = pkin(2) * t290;
t265 = pkin(9) * t235 - qJ(3) * t238;
t152 = t265 * t310 + t214;
t182 = pkin(8) * t289 + t235 * t301;
t154 = pkin(3) * t289 + t182;
t78 = -t234 * t152 + t154 * t237;
t428 = -t234 * t283 - (pkin(5) * t234 * t235 - t238 * t361) * t310 + t78;
t426 = t431 * t349;
t308 = qJD(2) * t231;
t285 = qJD(1) * t308;
t279 = t235 * t285;
t113 = qJD(4) * t167 - t234 * t279;
t353 = -t113 / 0.2e1;
t114 = -qJD(4) * t198 + t221 * t304 - t237 * t279;
t400 = -t114 / 0.2e1;
t195 = t234 * pkin(4) - qJ(5) * t237 + qJ(3);
t189 = pkin(10) * t234 + t195;
t197 = t335 * t237;
t125 = -t189 * t233 + t197 * t236;
t424 = qJD(6) * t125 + t233 * t428 + t236 * t429;
t126 = t189 * t236 + t197 * t233;
t423 = -qJD(6) * t126 - t233 * t429 + t236 * t428;
t315 = t235 * t237;
t79 = t237 * t152 + t234 * t154;
t422 = -(pkin(5) * t315 + qJ(5) * t238) * t310 - t79 - t237 * t283;
t208 = qJD(2) * t219;
t210 = t221 * qJD(3);
t318 = t231 * t235;
t282 = t362 * t318;
t260 = qJD(2) * t282;
t102 = -qJD(1) * t260 + t208 + t210;
t284 = -qJ(3) * t235 - pkin(1);
t146 = (t238 * t240 + t284) * t231;
t127 = qJD(1) * t146;
t225 = t232 * t235 * pkin(1);
t317 = t231 * t238;
t155 = (t317 * t362 + t225) * qJD(2);
t131 = qJD(1) * t155;
t207 = pkin(2) * t279;
t248 = t221 * t240 + t290 * t362 + qJD(3) - t219;
t306 = qJD(3) * t235;
t249 = (qJD(2) * t265 - t306) * t231;
t415 = qJD(1) * t249 + qJD(4) * t248 + t207;
t23 = -t127 * t304 + t131 * t237 - t234 * t415;
t278 = t238 * t285;
t21 = -pkin(4) * t278 - t23;
t307 = qJD(2) * t238;
t287 = t231 * t307;
t259 = t361 * t287;
t11 = -pkin(5) * t113 - qJD(1) * t259 - t23;
t243 = qJ(5) * t113 - qJD(5) * t168 + t102;
t20 = t114 * t361 + t243;
t56 = t127 * t234 - t237 * t248;
t258 = pkin(5) * t168 + t56;
t411 = qJD(5) + t258;
t34 = -t209 * t361 + t411;
t211 = t221 * qJ(3);
t117 = t211 + t154;
t256 = -qJ(5) * t168 + t117;
t42 = t167 * t361 + t256;
t5 = -t233 * t42 + t236 * t34;
t1 = qJD(6) * t5 + t11 * t233 + t20 * t236;
t6 = t233 * t34 + t236 * t42;
t2 = -qJD(6) * t6 + t11 * t236 - t20 * t233;
t277 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t29 = pkin(4) * t114 + t243;
t397 = -t278 / 0.2e1;
t49 = -qJD(6) * t104 + t114 * t236 - t233 * t278;
t46 = Ifges(7,6) * t49;
t48 = qJD(6) * t103 + t114 * t233 + t236 * t278;
t47 = Ifges(7,5) * t48;
t7 = -Ifges(7,3) * t113 + t46 + t47;
t421 = t277 + Ifges(6,4) * t397 + t102 * mrSges(5,2) + t21 * mrSges(6,1) + Ifges(5,5) * t278 / 0.2e1 + t7 / 0.2e1 - t23 * mrSges(5,3) - t29 * mrSges(6,3) + t431 * t400 + (Ifges(6,2) + Ifges(5,1)) * t353;
t181 = pkin(8) * t290 - t219;
t420 = -qJD(3) - t181;
t58 = pkin(4) * t167 + t256;
t419 = t5 * mrSges(7,1) + t117 * mrSges(5,2) - t6 * mrSges(7,2) - t58 * mrSges(6,3) + t426 - t432;
t398 = -t221 / 0.2e1;
t418 = mrSges(4,1) + mrSges(3,3);
t417 = mrSges(4,2) - mrSges(3,1);
t276 = t1 * t233 + t2 * t236;
t338 = t167 * pkin(5);
t57 = t237 * t127 + t234 * t248;
t45 = -t209 * qJ(5) - t57;
t35 = -t45 - t338;
t414 = t209 * t35 - t276;
t387 = -qJD(5) - t56;
t43 = -pkin(4) * t209 - t387;
t413 = -t43 * mrSges(6,1) - t56 * mrSges(5,3);
t412 = t45 * mrSges(6,1) - t57 * mrSges(5,3);
t348 = t167 / 0.2e1;
t410 = -Ifges(5,4) * t349 + Ifges(6,6) * t348 - t419 + t432;
t409 = Ifges(6,6) * t347 + Ifges(6,3) * t348;
t399 = t114 / 0.2e1;
t401 = t113 / 0.2e1;
t408 = Ifges(6,6) * t401 + Ifges(6,3) * t399;
t22 = -t127 * t305 + t234 * t131 + t237 * t415;
t13 = -qJ(5) * t278 - qJD(5) * t209 - t22;
t405 = Ifges(6,5) / 0.2e1;
t407 = t102 * mrSges(5,1) + t13 * mrSges(6,1) - t29 * mrSges(6,2) - t22 * mrSges(5,3) + Ifges(5,4) * t401 + Ifges(5,2) * t399 + Ifges(5,6) * t397 + t278 * t405 + t408;
t404 = -Ifges(5,6) / 0.2e1;
t396 = -t310 / 0.2e1;
t390 = Ifges(6,5) - Ifges(5,6);
t386 = -Ifges(3,4) * t289 / 0.2e1 + Ifges(3,5) * t398;
t298 = -Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t385 = -t298 * t209 - t413 + t419 + t425;
t288 = t235 * t308;
t216 = pkin(2) * t288;
t123 = t216 + t249;
t222 = pkin(8) * t318;
t340 = pkin(1) * t238;
t291 = -pkin(2) - t340;
t130 = pkin(3) * t318 + t222 + (-pkin(9) + t291) * t232;
t254 = -t237 * t123 - t130 * t304 + t146 * t305 - t234 * t155;
t26 = -t231 * (qJ(5) * t307 + qJD(5) * t235) + t254;
t75 = -mrSges(7,2) * t166 + mrSges(7,3) * t103;
t76 = mrSges(7,1) * t166 - mrSges(7,3) * t104;
t263 = t233 * t76 - t236 * t75;
t32 = -mrSges(7,1) * t113 - mrSges(7,3) * t48;
t33 = mrSges(7,2) * t113 + mrSges(7,3) * t49;
t264 = t233 * t33 + t236 * t32;
t384 = t263 * qJD(6) - t264;
t294 = Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1;
t295 = t404 + t405;
t382 = t117 * mrSges(5,1) - t58 * mrSges(6,2) + Ifges(5,4) * t347 + Ifges(6,5) * t344 + Ifges(5,2) * t348 + t209 * t404 + t409;
t383 = t294 * t168 - t295 * t209 - t382 - t412;
t172 = (-pkin(2) * t238 + t284) * t231;
t157 = qJD(1) * t172;
t299 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t327 = Ifges(4,6) * t238;
t381 = t295 * t167 - t298 * t168 + t299 * t209 + t181 * mrSges(3,3) + t43 * mrSges(6,2) + Ifges(3,1) * t290 / 0.2e1 + Ifges(4,4) * t398 + (-Ifges(4,2) * t235 - t327) * t396 + Ifges(5,5) * t346 + Ifges(5,6) * t349 + Ifges(6,4) * t347 + Ifges(6,5) * t348 - t157 * mrSges(4,3) - t45 * mrSges(6,3) - t56 * mrSges(5,1) - t57 * mrSges(5,2) - t386 + (Ifges(5,3) + Ifges(6,1)) * t344;
t331 = Ifges(7,4) * t233;
t268 = Ifges(7,2) * t236 + t331;
t330 = Ifges(7,4) * t236;
t270 = Ifges(7,1) * t233 + t330;
t273 = mrSges(7,1) * t236 - mrSges(7,2) * t233;
t274 = t233 * t5 - t236 * t6;
t325 = Ifges(7,6) * t236;
t328 = Ifges(7,5) * t233;
t343 = -t233 / 0.2e1;
t351 = -t166 / 0.2e1;
t355 = -t104 / 0.2e1;
t357 = -t103 / 0.2e1;
t332 = Ifges(7,4) * t104;
t40 = Ifges(7,2) * t103 + Ifges(7,6) * t166 + t332;
t373 = -t40 / 0.2e1;
t100 = Ifges(7,4) * t103;
t41 = Ifges(7,1) * t104 + Ifges(7,5) * t166 + t100;
t379 = mrSges(7,3) * t274 + t236 * t373 + (t325 + t328) * t351 + t268 * t357 + t270 * t355 + t35 * t273 + t343 * t41;
t378 = -Ifges(5,4) * t346 - Ifges(5,2) * t349 + t344 * t390 + t382 + t409;
t9 = Ifges(7,1) * t48 + Ifges(7,4) * t49 - Ifges(7,5) * t113;
t376 = t9 / 0.2e1;
t375 = Ifges(3,5) / 0.2e1;
t374 = Ifges(4,5) / 0.2e1;
t372 = t40 / 0.2e1;
t371 = t41 / 0.2e1;
t370 = t48 / 0.2e1;
t369 = t49 / 0.2e1;
t360 = m(6) * t21;
t359 = pkin(1) * mrSges(3,1);
t358 = pkin(1) * mrSges(3,2);
t341 = t236 / 0.2e1;
t91 = -t113 * mrSges(6,1) + mrSges(6,2) * t278;
t92 = mrSges(5,1) * t278 + mrSges(5,3) * t113;
t334 = t92 - t91;
t90 = mrSges(6,1) * t114 - mrSges(6,3) * t278;
t93 = -mrSges(5,2) * t278 - mrSges(5,3) * t114;
t333 = t93 - t90;
t322 = t221 * Ifges(4,5);
t118 = mrSges(6,1) * t167 - mrSges(6,3) * t209;
t59 = -mrSges(7,1) * t103 + mrSges(7,2) * t104;
t321 = t59 - t118;
t177 = -mrSges(4,1) * t289 - mrSges(4,3) * t221;
t96 = mrSges(5,1) * t167 + mrSges(5,2) * t168;
t320 = t96 - t177;
t319 = qJ(5) * t167;
t316 = t233 * t237;
t314 = t236 * t237;
t120 = -mrSges(5,2) * t209 - mrSges(5,3) * t167;
t313 = -t118 + t120;
t119 = mrSges(6,1) * t168 + mrSges(6,2) * t209;
t121 = mrSges(5,1) * t209 - mrSges(5,3) * t168;
t312 = t119 - t121;
t74 = t234 * t130 + t237 * t146;
t311 = t221 * t417 + t290 * t418;
t188 = pkin(8) * t317 + t225;
t303 = t232 * t340;
t302 = -qJD(2) + t221 / 0.2e1;
t300 = -Ifges(5,1) / 0.2e1 - Ifges(6,2) / 0.2e1;
t297 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t296 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t292 = t234 * t317;
t171 = -t232 * qJ(3) - t188;
t73 = t130 * t237 - t234 * t146;
t145 = pkin(3) * t317 - t171;
t280 = t234 * t290;
t275 = t233 * t6 + t236 * t5;
t272 = mrSges(7,1) * t233 + mrSges(7,2) * t236;
t271 = Ifges(7,1) * t236 - t331;
t269 = -Ifges(7,2) * t233 + t330;
t267 = Ifges(7,5) * t236 - Ifges(7,6) * t233;
t186 = t232 * t237 - t292;
t44 = pkin(5) * t186 - t318 * t361 - t73;
t185 = t232 * t234 + t237 * t317;
t255 = -qJ(5) * t186 + t145;
t61 = t185 * t361 + t255;
t16 = -t233 * t61 + t236 * t44;
t17 = t233 * t44 + t236 * t61;
t262 = t234 * t43 - t237 * t45;
t261 = t234 * t56 + t237 * t57;
t220 = qJD(2) * t303;
t183 = -pkin(8) * t288 + t220;
t68 = -qJ(5) * t318 - t74;
t31 = -t234 * t123 - t130 * t305 - t146 * t304 + t155 * t237;
t143 = t185 * t233 + t236 * t318;
t142 = t185 * t236 - t233 * t318;
t169 = -pkin(8) * t279 + t208;
t228 = t232 * qJD(3);
t129 = t220 + t228 - t260;
t253 = (-qJ(3) * t307 - t306) * t231;
t252 = Ifges(3,6) * t398 + (Ifges(3,4) * t235 + Ifges(3,2) * t238) * t396 + t322 / 0.2e1 + (-Ifges(4,6) * t235 - Ifges(4,3) * t238) * t310 / 0.2e1 - t157 * mrSges(4,2) - t182 * mrSges(3,3);
t184 = t188 * qJD(2);
t137 = -t169 - t210;
t170 = qJD(1) * t184;
t251 = -t169 * mrSges(3,2) - t137 * mrSges(4,3) + t170 * t417;
t250 = t23 * mrSges(5,1) - t22 * mrSges(5,2) + t21 * mrSges(6,2) - t13 * mrSges(6,3);
t140 = qJD(4) * t185 - t234 * t288;
t245 = qJ(5) * t140 - qJD(5) * t186 + t129;
t206 = Ifges(6,1) * t278;
t205 = Ifges(3,5) * t278;
t204 = Ifges(4,5) * t279;
t203 = Ifges(5,3) * t278;
t196 = t335 * t234;
t187 = -t222 + t303;
t180 = -qJ(3) * t289 + t214;
t179 = (mrSges(4,2) * t238 - mrSges(4,3) * t235) * t310;
t176 = -mrSges(3,2) * t221 + mrSges(3,3) * t289;
t174 = -qJD(5) * t237 + t286;
t173 = t232 * t291 + t222;
t165 = -t183 - t228;
t162 = (-t233 * t315 + t236 * t238) * t310;
t161 = (-t233 * t238 - t235 * t314) * t310;
t159 = t221 * t236 - t233 * t280;
t158 = -t221 * t233 - t236 * t280;
t156 = t216 + t253;
t153 = -qJD(1) * t282 + t219;
t151 = -t211 - t182;
t144 = -pkin(2) * t221 - t420;
t141 = -qJD(4) * t292 + t232 * t304 - t237 * t288;
t134 = qJD(1) * t253 + t207;
t111 = Ifges(6,4) * t113;
t110 = Ifges(5,5) * t113;
t109 = Ifges(6,5) * t114;
t108 = Ifges(5,6) * t114;
t97 = -mrSges(6,2) * t167 - mrSges(6,3) * t168;
t95 = pkin(4) * t168 + t319;
t89 = t219 + (-pkin(4) * t237 + t281) * t290;
t77 = pkin(4) * t185 + t255;
t72 = -pkin(4) * t289 - t78;
t71 = -qJ(5) * t289 - t79;
t70 = t168 * t361 + t319;
t69 = -pkin(4) * t318 - t73;
t67 = qJD(6) * t142 + t141 * t233 + t236 * t287;
t66 = -qJD(6) * t143 + t141 * t236 - t233 * t287;
t65 = mrSges(5,1) * t114 - mrSges(5,2) * t113;
t64 = -mrSges(6,2) * t114 + mrSges(6,3) * t113;
t50 = -pkin(5) * t185 - t68;
t38 = t57 - t338;
t36 = pkin(4) * t141 + t245;
t28 = t141 * t361 + t245;
t27 = -pkin(4) * t287 - t31;
t19 = -pkin(5) * t141 - t26;
t18 = -pkin(5) * t140 - t259 - t31;
t15 = t233 * t38 + t236 * t70;
t14 = -t233 * t70 + t236 * t38;
t12 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t10 = -pkin(5) * t114 - t13;
t8 = Ifges(7,4) * t48 + Ifges(7,2) * t49 - Ifges(7,6) * t113;
t4 = -qJD(6) * t17 + t18 * t236 - t233 * t28;
t3 = qJD(6) * t16 + t18 * t233 + t236 * t28;
t24 = [(t378 + t412) * t141 + (t410 + t413) * t140 + ((Ifges(5,1) + Ifges(7,3)) * t353 - t298 * t278 + Ifges(7,6) * t369 + Ifges(7,5) * t370 - Ifges(6,6) * t399 + Ifges(5,4) * t400 - Ifges(6,2) * t401 + t421) * t186 + m(5) * (t102 * t145 + t117 * t129 + t22 * t74 + t23 * t73 - t254 * t57 - t31 * t56) - t254 * t120 + ((t151 * mrSges(4,1) + (-Ifges(3,6) / 0.2e1 + t374) * t221 + t252) * t235 + (t144 * mrSges(4,1) + (t375 - Ifges(4,4) / 0.2e1) * t221 + t381) * t238) * t308 + (Ifges(7,1) * t67 + Ifges(7,4) * t66) * t354 + (Ifges(7,1) * t143 + Ifges(7,4) * t142) * t370 + ((t173 * mrSges(4,1) - t187 * mrSges(3,3) - t172 * mrSges(4,3) + (-Ifges(4,4) + t375) * t232 + (t238 * t296 - 0.2e1 * t358) * t231) * t238 + (t171 * mrSges(4,1) - t172 * mrSges(4,2) - t188 * mrSges(3,3) + (-Ifges(3,6) + t374) * t232 + (-t235 * t296 - 0.2e1 * t359) * t231 + (-0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(4,2) + t299) * t317) * t235) * t285 + (Ifges(7,4) * t143 + Ifges(7,2) * t142) * t369 + t311 * t184 + (-Ifges(5,4) * t353 - Ifges(5,2) * t400 + t278 * t295 + t407 + t408) * t185 + (Ifges(7,5) * t67 + Ifges(7,6) * t66) * t350 + (Ifges(7,5) * t143 + Ifges(7,6) * t142) * t353 + (t204 / 0.2e1 + t205 / 0.2e1 + t251) * t232 + t77 * t64 + t68 * t90 + t3 * t75 + t4 * t76 + t35 * (-mrSges(7,1) * t66 + mrSges(7,2) * t67) + t19 * t59 + t50 * t12 + t16 * t32 + t17 * t33 + (Ifges(7,4) * t67 + Ifges(7,2) * t66) * t356 + m(4) * (t134 * t172 + t137 * t171 + t144 * t184 + t151 * t165 + t156 * t157 + t170 * t173) + m(3) * (t169 * t188 - t170 * t187 + t181 * t184 + t182 * t183) + m(6) * (t13 * t68 + t21 * t69 + t26 * t45 + t27 * t43 + t29 * t77 + t36 * t58) + m(7) * (t1 * t17 + t10 * t50 + t16 * t2 + t19 * t35 + t3 * t6 + t4 * t5) + (t1 * t142 - t143 * t2 - t5 * t67 + t6 * t66) * mrSges(7,3) + t69 * t91 + t73 * t92 + t74 * t93 + t36 * t97 + ((-mrSges(4,1) * t137 + mrSges(4,2) * t134 + mrSges(3,3) * t169) * t238 + (-t110 / 0.2e1 - t108 / 0.2e1 + t203 / 0.2e1 + t206 / 0.2e1 + t111 / 0.2e1 + t109 / 0.2e1 - t134 * mrSges(4,3) + t418 * t170 + t295 * t114 + t298 * t113 + t250) * t235) * t231 + t26 * t118 + t27 * t119 + t31 * t121 + t129 * t96 + t142 * t8 / 0.2e1 + t10 * (-mrSges(7,1) * t142 + mrSges(7,2) * t143) + t145 * t65 + t67 * t371 + t66 * t372 + t143 * t376 + t165 * t177 + t156 * t179 + t183 * t176; t422 * t59 + (t35 * (-mrSges(7,1) * t314 + mrSges(7,2) * t316) + (Ifges(7,4) * t316 + Ifges(7,2) * t314) * t356 + (Ifges(7,1) * t316 + Ifges(7,4) * t314) * t354 + (Ifges(7,5) * t316 + Ifges(7,6) * t314) * t350 + t314 * t372 + t316 * t371 - t261 * mrSges(5,3) - t262 * mrSges(6,1) + (m(5) * t261 + m(6) * t262) * t240 + (t314 * t6 - t316 * t5) * mrSges(7,3) + (t313 * t240 + t378) * t237 + (t312 * t240 + t410) * t234) * qJD(4) + (t233 * t376 + t8 * t341 - t10 * t273 + t270 * t370 + t268 * t369 + t297 * t114 + (-t328 / 0.2e1 - t325 / 0.2e1 + t294) * t113 + (t1 * t236 - t2 * t233) * mrSges(7,3) + (m(5) * t22 - m(6) * t13 + t333) * t240 + (-mrSges(7,3) * t275 + t267 * t350 + t269 * t356 + t271 * t354 + t272 * t35 + t341 * t41 + t343 * t40) * qJD(6) + t407) * t234 + (-pkin(2) * t170 - qJ(3) * t137 - t144 * t182 + t151 * t420 - t157 * t180) * m(4) + (t47 / 0.2e1 + t46 / 0.2e1 - t294 * t114 + (-Ifges(7,3) / 0.2e1 + t300) * t113 + (m(5) * t23 + t334 - t360) * t240 + t421) * t237 + t251 + (-t6 * t161 + t5 * t162) * mrSges(7,3) + t423 * t76 + (t1 * t126 - t10 * t196 + t125 * t2 + t35 * t422 + t423 * t5 + t424 * t6) * m(7) + t424 * t75 - t311 * t182 + t320 * qJD(3) + qJ(3) * t65 + t204 + t205 - m(5) * (t117 * t153 - t56 * t78 + t57 * t79) - m(6) * (t43 * t72 + t45 * t71 + t58 * t89) + (-t89 + t174) * t97 + (-t177 + t176) * t181 + m(6) * (t58 * t174 + t29 * t195) + m(5) * (t102 * qJ(3) + t117 * qJD(3)) + (Ifges(7,5) * t162 + Ifges(7,6) * t161) * t351 + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t355 - t71 * t118 - t72 * t119 - t79 * t120 - t78 * t121 + t125 * t32 + t126 * t33 + ((t302 * Ifges(4,4) + (-pkin(2) * qJD(2) - t144) * mrSges(4,1) + (-t327 / 0.2e1 + t358) * t310 + (t390 * t234 - t391 * t237) * qJD(2) / 0.2e1 - t381 + t386) * t238 + (-t322 / 0.2e1 + (((Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t235 + t359) * t231 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t317) * qJD(1) + t302 * Ifges(3,6) + (-qJ(3) * qJD(2) - t151) * mrSges(4,1) + (t167 * t297 - t383) * t237 + (t167 * t294 + t168 * t300 - t385) * t234 - t252) * t235) * t310 - t153 * t96 - t162 * t41 / 0.2e1 - t35 * (-mrSges(7,1) * t161 + mrSges(7,2) * t162) + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t357 + t161 * t373 - t180 * t179 + t195 * t64 - t196 * t12; -t158 * t76 - t159 * t75 + (mrSges(4,1) * t307 + t179 * t235) * t310 + (-t97 - t320) * t221 + (t12 + t312 * t290 + (t233 * t75 + t236 * t76 + t312) * qJD(4) + t333) * t234 + (t334 + t209 * (t59 + t313) + t384) * t237 + (-t158 * t5 - t159 * t6 + (qJD(4) * t275 + t10) * t234 + (qJD(6) * t274 + t414) * t237) * m(7) + (-t13 * t234 + t209 * t262 - t21 * t237 - t221 * t58) * m(6) + (-t117 * t221 + t209 * t261 + t22 * t234 + t23 * t237) * m(5) + (t151 * t221 + t157 * t290 + t170) * m(4); (qJ(5) * t10 - t14 * t5 - t15 * t6 + t411 * t35) * m(7) + (t385 + t426) * t167 + t258 * t59 - ((-m(7) * t274 - t263) * qJD(6) + m(7) * t276 + t264) * t361 + (-pkin(4) * t21 - qJ(5) * t13 + t387 * t45 - t43 * t57 - t58 * t95) * m(6) + ((-t297 - t300) * t167 + t379 + t383) * t168 + t379 * qJD(6) - t276 * mrSges(7,3) + t250 - t312 * t57 + t313 * t56 + t9 * t341 + t8 * t343 + t10 * t272 + t321 * qJD(5) + (-t90 + t12) * qJ(5) + t111 + t109 - t110 - t15 * t75 - t14 * t76 + t203 + t206 + t267 * t353 - pkin(4) * t91 - t95 * t97 - t108 + t269 * t369 + t271 * t370; -t321 * t209 + (-t263 + t97) * t168 + t360 - m(6) * (-t168 * t58 - t209 * t45) + t91 + (-t166 * t274 - t414) * m(7) - t384; -t35 * (mrSges(7,1) * t104 + mrSges(7,2) * t103) + (Ifges(7,1) * t103 - t332) * t355 + t40 * t354 + (Ifges(7,5) * t103 - Ifges(7,6) * t104) * t351 - t5 * t75 + t6 * t76 + (t103 * t5 + t104 * t6) * mrSges(7,3) + t277 + t7 + (-Ifges(7,2) * t104 + t100 + t41) * t357;];
tauc  = t24(:);
