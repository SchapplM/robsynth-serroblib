% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:19
% EndTime: 2019-03-09 09:01:50
% DurationCPUTime: 16.98s
% Computational Cost: add. (11387->724), mult. (34305->987), div. (0->0), fcn. (26711->10), ass. (0->342)
t323 = cos(pkin(11));
t302 = t323 * pkin(2);
t228 = -t302 - pkin(3);
t226 = -pkin(9) + t228;
t238 = cos(qJ(5));
t307 = qJD(5) * t238;
t235 = sin(qJ(5));
t239 = cos(qJ(2));
t233 = cos(pkin(6));
t363 = pkin(1) * t233;
t225 = t239 * t363;
t220 = qJD(1) * t225;
t232 = sin(pkin(6));
t236 = sin(qJ(2));
t349 = pkin(8) + qJ(3);
t295 = t349 * t236;
t286 = t232 * t295;
t169 = -qJD(1) * t286 + t220;
t231 = sin(pkin(11));
t224 = t236 * t363;
t318 = t232 * t239;
t397 = t318 * t349 + t224;
t170 = t397 * qJD(1);
t290 = t323 * t170;
t118 = t169 * t231 + t290;
t282 = t323 * t318;
t257 = qJD(1) * t282;
t311 = qJD(1) * t232;
t301 = t236 * t311;
t188 = t231 * t301 - t257;
t361 = pkin(4) * t188;
t86 = t118 - t361;
t195 = (t231 * t239 + t236 * t323) * t232;
t191 = qJD(1) * t195;
t218 = pkin(2) * t301;
t289 = qJ(4) * t188 + t218;
t384 = pkin(3) + pkin(9);
t91 = t191 * t384 + t289;
t42 = t235 * t86 + t238 * t91;
t424 = t226 * t307 - t42;
t160 = t231 * t170;
t119 = t169 * t323 - t160;
t423 = t119 - qJD(4);
t351 = -mrSges(5,2) + mrSges(4,1);
t281 = pkin(5) * t238 + pkin(10) * t235;
t422 = qJD(5) * t281 - (-pkin(4) - t281) * t191 - t423;
t421 = pkin(10) * t188 + t424;
t213 = qJD(2) * t220;
t309 = qJD(2) * t232;
t294 = qJD(1) * t309;
t285 = t236 * t294;
t192 = -pkin(8) * t285 + t213;
t312 = pkin(8) * t318 + t224;
t201 = t312 * qJD(2);
t193 = qJD(1) * t201;
t222 = qJD(1) * t233 + qJD(2);
t247 = (-qJD(2) * t295 + qJD(3) * t239) * t232;
t139 = qJD(1) * t247 + t213;
t319 = t232 * t236;
t152 = -qJD(2) * t397 - qJD(3) * t319;
t242 = t152 * qJD(1);
t82 = t323 * t139 + t231 * t242;
t75 = -t222 * qJD(4) - t82;
t81 = t139 * t231 - t323 * t242;
t420 = -t193 * mrSges(3,1) - t192 * mrSges(3,2) - t82 * mrSges(4,2) - t75 * mrSges(5,3) - t351 * t81;
t148 = t188 * t238 - t222 * t235;
t189 = qJD(2) * t195;
t181 = qJD(1) * t189;
t107 = qJD(5) * t148 + t181 * t235;
t149 = t188 * t235 + t222 * t238;
t185 = qJD(5) + t191;
t234 = sin(qJ(6));
t237 = cos(qJ(6));
t109 = -t149 * t234 + t185 * t237;
t182 = qJD(2) * t257 - t231 * t285;
t47 = qJD(6) * t109 + t107 * t237 + t182 * t234;
t388 = t47 / 0.2e1;
t110 = t149 * t237 + t185 * t234;
t48 = -qJD(6) * t110 - t107 * t234 + t182 * t237;
t387 = t48 / 0.2e1;
t308 = qJD(5) * t235;
t61 = pkin(4) * t182 + t81;
t212 = pkin(2) * t285;
t254 = -qJ(4) * t182 - qJD(4) * t191 + t212;
t62 = t181 * t384 + t254;
t150 = pkin(2) * t222 + t169;
t101 = t150 * t323 - t160;
t255 = qJD(4) - t101;
t359 = t191 * pkin(4);
t63 = -t222 * t384 + t255 + t359;
t207 = (-pkin(2) * t239 - pkin(1)) * t232;
t202 = qJD(1) * t207 + qJD(3);
t245 = -t191 * qJ(4) + t202;
t84 = t188 * t384 + t245;
t7 = t235 * t61 + t238 * t62 + t63 * t307 - t308 * t84;
t419 = t7 * mrSges(6,2);
t32 = t235 * t63 + t238 * t84;
t8 = -qJD(5) * t32 - t235 * t62 + t238 * t61;
t418 = t8 * mrSges(6,1);
t108 = qJD(5) * t149 - t238 * t181;
t381 = t108 / 0.2e1;
t31 = -t235 * t84 + t238 * t63;
t417 = t31 * mrSges(6,1);
t416 = t32 * mrSges(6,2);
t415 = Ifges(4,5) - Ifges(5,4);
t414 = -Ifges(4,6) + Ifges(5,5);
t20 = -mrSges(7,1) * t48 + mrSges(7,2) * t47;
t76 = mrSges(6,1) * t182 - mrSges(6,3) * t107;
t413 = t76 - t20;
t412 = t148 * Ifges(6,6);
t411 = t149 * Ifges(6,5);
t410 = t182 * Ifges(6,5);
t409 = t182 * Ifges(6,6);
t408 = t185 * Ifges(6,3);
t102 = t231 * t150 + t290;
t93 = -t222 * qJ(4) - t102;
t69 = -t93 - t361;
t407 = t69 * (mrSges(6,1) * t238 - mrSges(6,2) * t235);
t362 = pkin(2) * t231;
t227 = qJ(4) + t362;
t205 = pkin(5) * t235 - pkin(10) * t238 + t227;
t317 = t234 * t235;
t157 = t205 * t237 - t226 * t317;
t406 = qJD(6) * t157 + t422 * t234 + t421 * t237;
t315 = t235 * t237;
t158 = t205 * t234 + t226 * t315;
t405 = -qJD(6) * t158 - t421 * t234 + t422 * t237;
t184 = Ifges(4,4) * t188;
t404 = t191 * Ifges(4,1) + t222 * Ifges(4,5) - t184 + t408 + t411 + t412;
t331 = t188 * mrSges(5,1);
t155 = -mrSges(5,3) * t222 + t331;
t96 = -mrSges(6,1) * t148 + mrSges(6,2) * t149;
t324 = t155 - t96;
t132 = -t188 * t237 - t191 * t317;
t296 = t234 * t308;
t306 = qJD(6) * t238;
t403 = t237 * t306 + t132 - t296;
t133 = -t188 * t234 + t191 * t315;
t402 = t234 * t306 + t237 * t308 + t133;
t330 = t191 * mrSges(5,1);
t346 = mrSges(4,3) * t191;
t313 = t351 * t222 - t330 - t346;
t321 = t191 * t238;
t401 = t307 + t321;
t400 = t235 * t7 + t238 * t8;
t57 = -pkin(4) * t181 - t75;
t25 = pkin(5) * t108 - pkin(10) * t107 + t57;
t5 = pkin(10) * t182 + t7;
t29 = pkin(10) * t185 + t32;
t37 = -pkin(5) * t148 - pkin(10) * t149 + t69;
t9 = -t234 * t29 + t237 * t37;
t1 = qJD(6) * t9 + t234 * t25 + t237 * t5;
t10 = t234 * t37 + t237 * t29;
t2 = -qJD(6) * t10 - t234 * t5 + t237 * t25;
t280 = t1 * t237 - t2 * t234;
t268 = Ifges(7,5) * t237 - Ifges(7,6) * t234;
t339 = Ifges(7,4) * t237;
t271 = -Ifges(7,2) * t234 + t339;
t340 = Ifges(7,4) * t234;
t274 = Ifges(7,1) * t237 - t340;
t276 = mrSges(7,1) * t234 + mrSges(7,2) * t237;
t278 = t10 * t234 + t9 * t237;
t28 = -pkin(5) * t185 - t31;
t106 = Ifges(7,4) * t109;
t147 = qJD(6) - t148;
t40 = Ifges(7,1) * t110 + Ifges(7,5) * t147 + t106;
t326 = t237 * t40;
t375 = t147 / 0.2e1;
t377 = t110 / 0.2e1;
t379 = t109 / 0.2e1;
t341 = Ifges(7,4) * t110;
t39 = Ifges(7,2) * t109 + Ifges(7,6) * t147 + t341;
t390 = -t39 / 0.2e1;
t399 = -t278 * mrSges(7,3) + t326 / 0.2e1 + t234 * t390 + t268 * t375 + t271 * t379 + t274 * t377 + t28 * t276;
t6 = -pkin(5) * t182 - t8;
t398 = qJD(5) * (t10 * t237 - t234 * t9) - t6;
t114 = -mrSges(6,2) * t185 + mrSges(6,3) * t148;
t64 = -mrSges(7,2) * t147 + mrSges(7,3) * t109;
t65 = mrSges(7,1) * t147 - mrSges(7,3) * t110;
t396 = qJD(5) * (t234 * t65 - t237 * t64 - t114) - t191 * t114 - t413;
t194 = t231 * t319 - t282;
t248 = -t195 * qJ(4) + t207;
t100 = t194 * t384 + t248;
t166 = pkin(2) * t233 + t225 - t286;
t186 = qJ(3) * t318 + t312;
t260 = -t166 * t323 + t231 * t186;
t85 = t195 * pkin(4) - t233 * t384 + t260;
t348 = t238 * t100 + t235 * t85;
t299 = t236 * t309;
t190 = qJD(2) * t282 - t231 * t299;
t219 = pkin(2) * t299;
t253 = -qJ(4) * t190 - qJD(4) * t195 + t219;
t73 = t189 * t384 + t253;
t221 = qJD(2) * t225;
t151 = t221 + t247;
t97 = t151 * t231 - t323 * t152;
t74 = pkin(4) * t190 + t97;
t17 = -qJD(5) * t348 - t235 * t73 + t238 * t74;
t395 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t394 = -0.2e1 * pkin(1);
t11 = Ifges(7,5) * t47 + Ifges(7,6) * t48 + Ifges(7,3) * t108;
t393 = t11 / 0.2e1;
t12 = t47 * Ifges(7,4) + t48 * Ifges(7,2) + t108 * Ifges(7,6);
t392 = t12 / 0.2e1;
t391 = Ifges(7,1) * t388 + Ifges(7,4) * t387 + Ifges(7,5) * t381;
t389 = t39 / 0.2e1;
t332 = t185 * Ifges(6,6);
t344 = Ifges(6,4) * t149;
t71 = t148 * Ifges(6,2) + t332 + t344;
t386 = -t71 / 0.2e1;
t143 = Ifges(6,4) * t148;
t333 = t185 * Ifges(6,5);
t72 = t149 * Ifges(6,1) + t143 + t333;
t385 = -t72 / 0.2e1;
t383 = t107 / 0.2e1;
t382 = -t108 / 0.2e1;
t380 = -t109 / 0.2e1;
t378 = -t110 / 0.2e1;
t376 = -t147 / 0.2e1;
t258 = t194 * t238 - t233 * t235;
t374 = t258 / 0.2e1;
t164 = t194 * t235 + t233 * t238;
t373 = t164 / 0.2e1;
t372 = -t188 / 0.2e1;
t371 = t188 / 0.2e1;
t369 = -t191 / 0.2e1;
t368 = t191 / 0.2e1;
t366 = t222 / 0.2e1;
t365 = -t236 / 0.2e1;
t364 = t238 / 0.2e1;
t356 = t238 * t6;
t354 = t31 * mrSges(6,3);
t353 = t32 * mrSges(6,3);
t350 = -Ifges(4,4) - Ifges(5,6);
t347 = mrSges(4,3) * t188;
t345 = Ifges(3,4) * t236;
t343 = Ifges(6,4) * t235;
t342 = Ifges(6,4) * t238;
t338 = Ifges(3,5) * t239;
t337 = t109 * Ifges(7,6);
t336 = t110 * Ifges(7,5);
t335 = t147 * Ifges(7,3);
t334 = t182 * mrSges(5,1);
t329 = t191 * Ifges(4,4);
t328 = t191 * Ifges(5,6);
t327 = t195 * t81;
t115 = mrSges(6,1) * t185 - mrSges(6,3) * t149;
t55 = -mrSges(7,1) * t109 + mrSges(7,2) * t110;
t325 = t115 - t55;
t322 = t191 * t235;
t316 = t234 * t238;
t314 = t237 * t238;
t98 = t323 * t151 + t231 * t152;
t123 = t231 * t166 + t323 * t186;
t38 = t335 + t336 + t337;
t305 = t38 * t364;
t304 = t238 * t386;
t303 = Ifges(6,5) * t107 - Ifges(6,6) * t108 + Ifges(6,3) * t182;
t90 = -t233 * qJD(4) - t98;
t116 = -t233 * qJ(4) - t123;
t300 = t239 * t311;
t298 = t226 * t308;
t288 = mrSges(3,3) * t301;
t287 = mrSges(3,3) * t300;
t275 = Ifges(6,1) * t235 + t342;
t273 = Ifges(7,1) * t234 + t339;
t272 = Ifges(6,2) * t238 + t343;
t270 = Ifges(7,2) * t237 + t340;
t269 = Ifges(6,5) * t235 + Ifges(6,6) * t238;
t267 = Ifges(7,5) * t234 + Ifges(7,6) * t237;
t26 = mrSges(7,1) * t108 - mrSges(7,3) * t47;
t27 = -mrSges(7,2) * t108 + mrSges(7,3) * t48;
t266 = -t234 * t26 + t237 * t27;
t36 = pkin(10) * t195 + t348;
t89 = -pkin(4) * t194 - t116;
t51 = -pkin(5) * t258 - pkin(10) * t164 + t89;
t19 = t234 * t51 + t237 * t36;
t18 = -t234 * t36 + t237 * t51;
t265 = -t234 * t64 - t237 * t65;
t263 = t31 * t235 - t32 * t238;
t41 = -t235 * t91 + t238 * t86;
t45 = -t100 * t235 + t238 * t85;
t130 = t164 * t237 + t195 * t234;
t129 = -t164 * t234 + t195 * t237;
t66 = -pkin(4) * t189 - t90;
t16 = -t100 * t308 + t235 * t74 + t238 * t73 + t85 * t307;
t252 = t222 * (-Ifges(3,6) * t236 + t338);
t244 = -qJD(5) * t263 + t400;
t243 = qJD(5) * t28 - qJD(6) * t278 + t280;
t77 = -mrSges(6,2) * t182 - mrSges(6,3) * t108;
t241 = qJD(6) * t265 - t185 * t325 + t266 + t77;
t217 = Ifges(3,4) * t300;
t211 = t294 * t338;
t203 = -pkin(8) * t319 + t225;
t200 = -pkin(8) * t299 + t221;
t199 = t312 * qJD(1);
t198 = -pkin(8) * t301 + t220;
t197 = -t222 * mrSges(3,2) + t287;
t196 = mrSges(3,1) * t222 - t288;
t183 = Ifges(5,6) * t188;
t177 = Ifges(5,4) * t182;
t176 = Ifges(4,5) * t182;
t175 = Ifges(5,5) * t181;
t174 = Ifges(4,6) * t181;
t172 = t182 * mrSges(5,3);
t171 = t182 * mrSges(4,2);
t168 = Ifges(3,1) * t301 + t222 * Ifges(3,5) + t217;
t167 = Ifges(3,6) * t222 + (Ifges(3,2) * t239 + t345) * t311;
t153 = -mrSges(4,2) * t222 - t347;
t141 = -t191 * t314 + t222 * t234;
t140 = t191 * t316 + t222 * t237;
t135 = -mrSges(5,2) * t188 - mrSges(5,3) * t191;
t134 = mrSges(4,1) * t188 + mrSges(4,2) * t191;
t131 = t194 * pkin(3) + t248;
t128 = pkin(3) * t191 + t289;
t126 = -t188 * Ifges(4,2) + t222 * Ifges(4,6) + t329;
t125 = t222 * Ifges(5,4) - t191 * Ifges(5,2) + t183;
t124 = t222 * Ifges(5,5) + t188 * Ifges(5,3) - t328;
t121 = qJD(5) * t164 - t189 * t238;
t120 = qJD(5) * t258 + t189 * t235;
t117 = -t233 * pkin(3) + t260;
t113 = t188 * pkin(3) + t245;
t99 = pkin(5) * t149 - pkin(10) * t148;
t95 = pkin(3) * t189 + t253;
t92 = -t222 * pkin(3) + t255;
t87 = t119 - t359;
t83 = pkin(3) * t181 + t254;
t54 = mrSges(6,1) * t108 + mrSges(6,2) * t107;
t53 = qJD(6) * t129 + t120 * t237 + t190 * t234;
t52 = -qJD(6) * t130 - t120 * t234 + t190 * t237;
t50 = t107 * Ifges(6,1) - t108 * Ifges(6,4) + t410;
t49 = t107 * Ifges(6,4) - t108 * Ifges(6,2) + t409;
t35 = -pkin(5) * t195 - t45;
t33 = pkin(5) * t188 - t41;
t30 = pkin(5) * t121 - pkin(10) * t120 + t66;
t24 = t234 * t99 + t237 * t31;
t23 = -t234 * t31 + t237 * t99;
t15 = -pkin(5) * t190 - t17;
t14 = pkin(10) * t190 + t16;
t4 = -qJD(6) * t19 - t14 * t234 + t237 * t30;
t3 = qJD(6) * t18 + t14 * t237 + t234 * t30;
t13 = [t130 * t391 + t129 * t392 + (t194 * t75 + t327) * mrSges(5,1) + (-t194 * t82 + t327) * mrSges(4,3) + (Ifges(7,5) * t53 + Ifges(7,6) * t52 + Ifges(7,3) * t121) * t375 + (Ifges(7,1) * t53 + Ifges(7,4) * t52 + Ifges(7,5) * t121) * t377 + (Ifges(7,4) * t53 + Ifges(7,2) * t52 + Ifges(7,6) * t121) * t379 + t121 * t386 + t52 * t389 + t185 * (Ifges(6,5) * t120 - Ifges(6,6) * t121) / 0.2e1 + (t404 / 0.2e1 + Ifges(4,1) * t368 - Ifges(5,2) * t369 - Ifges(5,6) * t371 + Ifges(4,4) * t372 + t202 * mrSges(4,2) + t411 / 0.2e1 + t412 / 0.2e1 + t417 - t416 - t113 * mrSges(5,3) - t125 / 0.2e1 + t408 / 0.2e1 - t101 * mrSges(4,3) + t92 * mrSges(5,1) + t415 * t366) * t190 + m(6) * (t16 * t32 + t17 * t31 + t348 * t7 + t45 * t8 + t57 * t89 + t66 * t69) + t348 * t77 + ((Ifges(3,5) * t233 / 0.2e1 - t203 * mrSges(3,3) + (mrSges(3,2) * t394 + 0.3e1 / 0.2e1 * Ifges(3,4) * t239) * t232) * t239 + (-t312 * mrSges(3,3) - Ifges(3,6) * t233 + (mrSges(3,1) * t394 - 0.3e1 / 0.2e1 * t345 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t239) * t232 + (m(4) * t207 + mrSges(4,1) * t194 + mrSges(4,2) * t195) * pkin(2)) * t236) * t294 + t148 * (Ifges(6,4) * t120 - Ifges(6,2) * t121) / 0.2e1 + (Ifges(6,1) * t164 + Ifges(6,4) * t258 + Ifges(6,5) * t195) * t383 + (Ifges(7,4) * t130 + Ifges(7,2) * t129 - Ifges(7,6) * t258) * t387 + (Ifges(7,1) * t130 + Ifges(7,4) * t129 - Ifges(7,5) * t258) * t388 + t1 * (mrSges(7,2) * t258 + mrSges(7,3) * t129) + t2 * (-mrSges(7,1) * t258 - mrSges(7,3) * t130) + t57 * (-mrSges(6,1) * t258 + mrSges(6,2) * t164) + m(3) * (t192 * t312 - t193 * t203 - t198 * t201 + t199 * t200) + (t192 * t239 + t193 * t236) * t232 * mrSges(3,3) + (Ifges(6,5) * t373 + Ifges(6,6) * t374 + t260 * mrSges(4,3) + t117 * mrSges(5,1) + (Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t233 + t350 * t194 + (Ifges(4,1) + Ifges(5,2) + Ifges(6,3) / 0.2e1) * t195) * t182 + m(4) * (-t101 * t97 + t102 * t98 + t123 * t82 + t202 * t219 + t260 * t81) - t258 * t393 + (-t120 * t31 - t121 * t32 - t164 * t8 + t258 * t7) * mrSges(6,3) + (Ifges(7,5) * t130 + Ifges(7,6) * t129 - Ifges(7,3) * t258) * t381 + (Ifges(6,4) * t164 + Ifges(6,2) * t258 + Ifges(6,6) * t195) * t382 + t50 * t373 + t49 * t374 + t149 * (Ifges(6,1) * t120 - Ifges(6,4) * t121) / 0.2e1 + (t211 / 0.2e1 + t176 / 0.2e1 - t174 / 0.2e1 - t177 / 0.2e1 + t175 / 0.2e1 + t420) * t233 + t207 * t171 + t83 * (-mrSges(5,2) * t194 - mrSges(5,3) * t195) + t200 * t197 - t201 * t196 + t98 * t153 + t90 * t155 + t95 * t135 + t6 * (-mrSges(7,1) * t129 + mrSges(7,2) * t130) + t120 * t72 / 0.2e1 + t10 * (-mrSges(7,2) * t121 + mrSges(7,3) * t52) + t9 * (mrSges(7,1) * t121 - mrSges(7,3) * t53) + t121 * t38 / 0.2e1 + (t239 * t168 / 0.2e1 + t167 * t365 + t236 * pkin(2) * t134 + t252 / 0.2e1 + (-t198 * t239 - t199 * t236) * mrSges(3,3)) * t309 + t69 * (mrSges(6,1) * t121 + mrSges(6,2) * t120) + t16 * t114 + t17 * t115 + (t207 * mrSges(4,1) + t116 * mrSges(5,1) - t131 * mrSges(5,2) - t123 * mrSges(4,3) + (-Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1) * t233 + t350 * t195 + (Ifges(4,2) + Ifges(5,3)) * t194) * t181 + m(7) * (t1 * t19 + t10 * t3 + t15 * t28 + t18 * t2 + t35 * t6 + t4 * t9) + t66 * t96 + m(5) * (t113 * t95 + t116 * t75 + t117 * t81 + t131 * t83 + t90 * t93 + t92 * t97) + t89 * t54 + t45 * t76 + t4 * t65 + t3 * t64 + t15 * t55 + t28 * (-mrSges(7,1) * t52 + mrSges(7,2) * t53) + t53 * t40 / 0.2e1 + t35 * t20 - t313 * t97 + t18 * t26 + t19 * t27 + t195 * t303 / 0.2e1 + t195 * t418 + (-Ifges(4,4) * t368 + Ifges(5,6) * t369 + Ifges(5,3) * t371 - Ifges(4,2) * t372 + t202 * mrSges(4,1) - t113 * mrSges(5,2) - t126 / 0.2e1 + t124 / 0.2e1 - t102 * mrSges(4,3) + t93 * mrSges(5,1) + t414 * t366) * t189 - t195 * t419 - t131 * t172; -(t234 * t40 + t237 * t39) * t306 / 0.2e1 - (t72 + t326) * t308 / 0.2e1 + (-t181 * t362 - t182 * t302) * mrSges(4,3) + t132 * t390 + t314 * t391 + t313 * t118 + (-t1 * t316 - t2 * t314 + t402 * t9) * mrSges(7,3) + (t28 * t403 + t401 * t9) * mrSges(7,1) + (-mrSges(7,2) * t401 - mrSges(7,3) * t403) * t10 + (-t33 + t298) * t55 + (-t41 - t298) * t115 + (-t153 + t155) * t119 - t402 * t28 * mrSges(7,2) + t424 * t114 + (t288 + t196) * t199 + t102 * t346 + t308 * t354 + t276 * t356 + (Ifges(7,5) * t133 + Ifges(7,6) * t132 - Ifges(7,3) * t321) * t376 + (Ifges(7,1) * t133 + Ifges(7,4) * t132 - Ifges(7,5) * t321) * t378 + (Ifges(7,4) * t133 + Ifges(7,2) * t132 - Ifges(7,6) * t321) * t380 + t322 * t385 + t296 * t389 + ((Ifges(7,3) * t238 - t235 * t268) * t375 + (Ifges(7,5) * t238 - t235 * t274) * t377 + (Ifges(7,6) * t238 - t235 * t271) * t379 + t304 + t305 + t407) * qJD(5) - t324 * qJD(4) + (-Ifges(4,2) * t191 - t184 + t404) * t371 + t405 * t65 + (-t28 * t33 + t1 * t158 + t157 * t2 + (t28 * t308 - t356) * t226 + t405 * t9 + t406 * t10) * m(7) + t406 * t64 + t413 * t226 * t238 + (-Ifges(4,1) * t188 + t124 - t329) * t369 + (t101 * t118 - t102 * t119 - t202 * t218 + (t231 * t82 - t323 * t81) * pkin(2)) * m(4) + t420 - (t148 * t272 + t149 * t275 + t185 * t269) * qJD(5) / 0.2e1 + (Ifges(5,3) * t191 + t125 + t183) * t372 + t191 * t304 + t191 * t305 - ((-Ifges(3,2) * t301 + t168 + t217) * t239 + t252) * t311 / 0.2e1 + t176 - t177 + t50 * t364 - t174 + t175 + (-t267 * t375 - t270 * t379 - t273 * t377) * t306 + t92 * t331 + t228 * t334 + t211 + (t226 * t244 + t227 * t57 - t31 * t41 - t32 * t42 + (qJD(4) - t87) * t69) * m(6) + (pkin(1) * (mrSges(3,1) * t236 + mrSges(3,2) * t239) + (Ifges(3,1) * t239 - t345) * t365) * qJD(1) ^ 2 * t232 ^ 2 + (Ifges(5,2) * t188 + t126 + t328) * t368 - t202 * (mrSges(4,1) * t191 - mrSges(4,2) * t188) - t113 * (-mrSges(5,2) * t191 + mrSges(5,3) * t188) + t157 * t26 + t158 * t27 - t128 * t135 - t133 * t40 / 0.2e1 - t307 * t353 - t101 * t347 - t87 * t96 - t93 * t330 - t12 * t316 / 0.2e1 + t167 * t301 / 0.2e1 - t134 * t218 + t188 * t417 + t191 * t407 + (-mrSges(5,1) * t181 + t54) * t227 + (t287 - t197) * t198 + (t31 * t322 - t32 * t321 - t400) * mrSges(6,3) + (-t113 * t128 - t118 * t92 - t227 * t75 + t228 * t81 + t423 * t93) * m(5) + t342 * t382 + (t393 + Ifges(7,3) * t381 - Ifges(6,2) * t382 + Ifges(7,6) * t387 + Ifges(7,5) * t388 + t226 * t77 - t49 / 0.2e1 - t409 / 0.2e1 + t57 * mrSges(6,1) + t395) * t235 + (t268 * t381 + Ifges(6,1) * t383 + t271 * t387 + t274 * t388 + t410 / 0.2e1 + t57 * mrSges(6,2)) * t238 - (-t188 * t415 + t191 * t414) * t222 / 0.2e1 - t188 * t416 - t343 * t383 - t185 * (-Ifges(6,3) * t188 + t191 * t269) / 0.2e1 - t148 * (-Ifges(6,6) * t188 + t191 * t272) / 0.2e1 - t149 * (-Ifges(6,5) * t188 + t191 * t275) / 0.2e1 - Ifges(3,6) * t285; -t132 * t65 - t133 * t64 + t171 - t172 + t313 * t191 + t351 * t181 + (t153 - t324) * t188 + t396 * t235 + t241 * t238 + (-t10 * t133 - t132 * t9 - t235 * t398 + t243 * t238 + t28 * t321) * m(7) + (t188 * t69 - t235 * t8 + t238 * t7 - t185 * (t235 * t32 + t238 * t31)) * m(6) + (-t188 * t93 - t191 * t92 + t83) * m(5) + (t101 * t191 + t102 * t188 + t212) * m(4); t334 + t191 * t135 - t140 * t65 - t141 * t64 + t324 * t222 - t396 * t238 + t241 * t235 + (-t10 * t141 - t140 * t9 + t243 * t235 + t238 * t398 + t28 * t322) * m(7) + (-t191 * t263 - t222 * t69 + t244) * m(6) + (t113 * t191 + t222 * t93 + t81) * m(5); t303 + t234 * t391 + t325 * t32 - t31 * t114 - t24 * t64 - t23 * t65 - pkin(5) * t20 - t419 + t418 + (-t38 / 0.2e1 + t71 / 0.2e1 + t353 + t332 / 0.2e1 - t335 / 0.2e1 - t337 / 0.2e1 - t336 / 0.2e1 - t69 * mrSges(6,1) + t10 * mrSges(7,2) - t9 * mrSges(7,1) + t344 / 0.2e1) * t149 + t6 * (-mrSges(7,1) * t237 + mrSges(7,2) * t234) + t237 * t392 + t267 * t381 + t270 * t387 + t273 * t388 + (-t333 / 0.2e1 - t143 / 0.2e1 - t69 * mrSges(6,2) + t385 + t354 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t149 - t399) * t148 + t399 * qJD(6) + t280 * mrSges(7,3) + (-pkin(5) * t6 - t10 * t24 - t23 * t9 - t28 * t32) * m(7) + (m(7) * t280 + t266 + (-m(7) * t278 + t265) * qJD(6)) * pkin(10); -t28 * (mrSges(7,1) * t110 + mrSges(7,2) * t109) + (Ifges(7,1) * t109 - t341) * t378 + t39 * t377 + (Ifges(7,5) * t109 - Ifges(7,6) * t110) * t376 - t9 * t64 + t10 * t65 + (t10 * t110 + t109 * t9) * mrSges(7,3) + t11 + (-Ifges(7,2) * t110 + t106 + t40) * t380 + t395;];
tauc  = t13(:);
