% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:51
% EndTime: 2019-03-09 08:50:20
% DurationCPUTime: 16.66s
% Computational Cost: add. (16068->678), mult. (41339->940), div. (0->0), fcn. (31889->10), ass. (0->300)
t269 = sin(pkin(10));
t276 = cos(qJ(2));
t324 = cos(pkin(10));
t291 = t324 * t276;
t273 = sin(qJ(2));
t310 = qJD(1) * t273;
t233 = -qJD(1) * t291 + t269 * t310;
t268 = sin(pkin(11));
t270 = cos(pkin(11));
t272 = sin(qJ(5));
t275 = cos(qJ(5));
t280 = t268 * t272 - t270 * t275;
t170 = t280 * t233;
t237 = t280 * qJD(5);
t398 = -t237 - t170;
t250 = t268 * t275 + t270 * t272;
t169 = t250 * t233;
t238 = t250 * qJD(5);
t397 = -t238 - t169;
t249 = t269 * t276 + t273 * t324;
t235 = t249 * qJD(1);
t300 = pkin(2) * t310;
t175 = pkin(3) * t235 + qJ(4) * t233 + t300;
t344 = -qJ(3) - pkin(7);
t259 = t344 * t276;
t253 = qJD(1) * t259;
t239 = t269 * t253;
t258 = t344 * t273;
t252 = qJD(1) * t258;
t196 = t252 * t324 + t239;
t121 = t268 * t175 + t270 * t196;
t316 = t233 * t268;
t107 = pkin(8) * t316 + t121;
t354 = pkin(2) * t269;
t262 = qJ(4) + t354;
t345 = pkin(8) + t262;
t244 = t345 * t268;
t245 = t345 * t270;
t304 = qJD(5) * t275;
t306 = qJD(4) * t270;
t307 = qJD(4) * t268;
t120 = t270 * t175 - t196 * t268;
t315 = t233 * t270;
t91 = pkin(4) * t235 + pkin(8) * t315 + t120;
t408 = -t244 * t304 + (-t107 + t306) * t275 + (-qJD(5) * t245 - t307 - t91) * t272;
t187 = -t272 * t244 + t275 * t245;
t407 = -t250 * qJD(4) - qJD(5) * t187 + t107 * t272 - t275 * t91;
t203 = qJD(2) * t268 + t235 * t270;
t288 = t270 * qJD(2) - t235 * t268;
t147 = t203 * t275 + t272 * t288;
t271 = sin(qJ(6));
t274 = cos(qJ(6));
t420 = -t203 * t272 + t275 * t288;
t428 = -t147 * t271 + t274 * t420;
t86 = t147 * t274 + t271 * t420;
t427 = pkin(9) * t397 + t408;
t426 = -pkin(5) * t235 - pkin(9) * t398 + t407;
t230 = qJD(5) + t233;
t265 = -pkin(2) * t276 - pkin(1);
t311 = qJD(1) * t265;
t255 = qJD(3) + t311;
t163 = t233 * pkin(3) - t235 * qJ(4) + t255;
t246 = qJD(2) * pkin(2) + t252;
t292 = t324 * t253;
t192 = t269 * t246 - t292;
t185 = qJD(2) * qJ(4) + t192;
t111 = t270 * t163 - t185 * t268;
t73 = pkin(4) * t233 - pkin(8) * t203 + t111;
t112 = t268 * t163 + t270 * t185;
t92 = pkin(8) * t288 + t112;
t44 = -t272 * t92 + t275 * t73;
t33 = -pkin(9) * t147 + t44;
t32 = pkin(5) * t230 + t33;
t45 = t272 * t73 + t275 * t92;
t34 = pkin(9) * t420 + t45;
t326 = t271 * t34;
t11 = t274 * t32 - t326;
t325 = t274 * t34;
t12 = t271 * t32 + t325;
t234 = t249 * qJD(2);
t222 = qJD(1) * t234;
t278 = -t269 * t273 + t291;
t236 = t278 * qJD(2);
t223 = qJD(1) * t236;
t98 = qJD(5) * t420 - t223 * t280;
t99 = -qJD(5) * t147 - t223 * t250;
t30 = qJD(6) * t428 + t271 * t99 + t274 * t98;
t31 = -qJD(6) * t86 - t271 * t98 + t274 * t99;
t302 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t222;
t355 = Ifges(7,4) * t86;
t221 = qJD(6) + t230;
t366 = -t221 / 0.2e1;
t375 = -t86 / 0.2e1;
t305 = qJD(5) * t272;
t317 = t223 * t270;
t303 = qJD(1) * qJD(2);
t296 = t273 * t303;
t287 = pkin(2) * t296;
t138 = pkin(3) * t222 - qJ(4) * t223 - qJD(4) * t235 + t287;
t293 = qJD(2) * t344;
t231 = qJD(3) * t276 + t273 * t293;
t213 = t231 * qJD(1);
t232 = -t273 * qJD(3) + t276 * t293;
t277 = qJD(1) * t232;
t159 = t324 * t213 + t269 * t277;
t154 = qJD(2) * qJD(4) + t159;
t77 = t270 * t138 - t154 * t268;
t59 = pkin(4) * t222 - pkin(8) * t317 + t77;
t318 = t223 * t268;
t78 = t268 * t138 + t270 * t154;
t61 = -pkin(8) * t318 + t78;
t15 = t272 * t59 + t275 * t61 + t73 * t304 - t305 * t92;
t10 = pkin(9) * t99 + t15;
t16 = -qJD(5) * t45 - t272 * t61 + t275 * t59;
t8 = pkin(5) * t222 - pkin(9) * t98 + t16;
t2 = qJD(6) * t11 + t10 * t274 + t271 * t8;
t3 = -qJD(6) * t12 - t10 * t271 + t274 * t8;
t416 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t191 = t246 * t324 + t239;
t178 = -qJD(2) * pkin(3) + qJD(4) - t191;
t139 = -pkin(4) * t288 + t178;
t80 = -pkin(5) * t420 + t139;
t423 = t302 + t416 + (Ifges(7,5) * t428 - Ifges(7,6) * t86) * t366 + (t11 * t428 + t12 * t86) * mrSges(7,3) - t80 * (mrSges(7,1) * t86 + mrSges(7,2) * t428) + (Ifges(7,1) * t428 - t355) * t375;
t361 = -t233 / 0.2e1;
t422 = Ifges(4,6) * qJD(2) / 0.2e1;
t81 = Ifges(7,4) * t428;
t421 = -Ifges(7,2) * t86 + t81;
t332 = t203 * Ifges(5,5);
t333 = t288 * Ifges(5,6);
t340 = Ifges(4,4) * t235;
t419 = t112 * mrSges(5,2) + t12 * mrSges(7,2) + t45 * mrSges(6,2) + Ifges(4,2) * t361 + t340 / 0.2e1 + t422 - t11 * mrSges(7,1) - t111 * mrSges(5,1) - t255 * mrSges(4,1) - t44 * mrSges(6,1) - t332 / 0.2e1 - t333 / 0.2e1;
t387 = t30 / 0.2e1;
t386 = t31 / 0.2e1;
t373 = t98 / 0.2e1;
t372 = t99 / 0.2e1;
t364 = t222 / 0.2e1;
t308 = qJD(2) * t273;
t417 = pkin(2) * t308;
t186 = -t275 * t244 - t245 * t272;
t156 = -pkin(9) * t250 + t186;
t157 = -pkin(9) * t280 + t187;
t93 = t156 * t274 - t157 * t271;
t415 = qJD(6) * t93 + t426 * t271 + t274 * t427;
t94 = t156 * t271 + t157 * t274;
t414 = -qJD(6) * t94 - t271 * t427 + t426 * t274;
t413 = mrSges(4,2) * t235;
t406 = Ifges(4,5) * qJD(2);
t195 = t252 * t269 - t292;
t155 = -pkin(4) * t316 + t195;
t402 = -pkin(5) * t397 - t155;
t190 = -pkin(3) * t278 - t249 * qJ(4) + t265;
t199 = t269 * t258 - t259 * t324;
t132 = t270 * t190 - t199 * t268;
t350 = pkin(8) * t270;
t106 = -pkin(4) * t278 - t249 * t350 + t132;
t133 = t268 * t190 + t270 * t199;
t313 = t249 * t268;
t114 = -pkin(8) * t313 + t133;
t55 = t272 * t106 + t275 * t114;
t110 = t169 * t271 + t170 * t274;
t193 = -t250 * t271 - t274 * t280;
t124 = qJD(6) * t193 - t237 * t274 - t238 * t271;
t401 = t124 - t110;
t109 = t169 * t274 - t170 * t271;
t194 = t250 * t274 - t271 * t280;
t125 = -qJD(6) * t194 + t237 * t271 - t238 * t274;
t400 = t125 - t109;
t342 = mrSges(4,3) * t235;
t399 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t288 - mrSges(5,2) * t203 - t342;
t396 = -t268 * t77 + t270 * t78;
t395 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t309 = qJD(1) * t276;
t322 = Ifges(3,6) * qJD(2);
t341 = Ifges(3,4) * t273;
t394 = pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t309) + t322 / 0.2e1 + (t276 * Ifges(3,2) + t341) * qJD(1) / 0.2e1;
t336 = Ifges(5,2) * t268;
t338 = Ifges(5,4) * t270;
t284 = -t336 + t338;
t339 = Ifges(5,4) * t268;
t285 = Ifges(5,1) * t270 - t339;
t286 = mrSges(5,1) * t268 + mrSges(5,2) * t270;
t356 = t270 / 0.2e1;
t357 = -t268 / 0.2e1;
t393 = (t203 * Ifges(5,1) + Ifges(5,4) * t288 + Ifges(5,5) * t233) * t356 + (t203 * Ifges(5,4) + Ifges(5,2) * t288 + Ifges(5,6) * t233) * t357 + t178 * t286 + t255 * mrSges(4,2) + t203 * t285 / 0.2e1 + t288 * t284 / 0.2e1;
t391 = -0.2e1 * pkin(1);
t389 = Ifges(7,4) * t387 + Ifges(7,2) * t386 + Ifges(7,6) * t364;
t388 = Ifges(7,1) * t387 + Ifges(7,4) * t386 + Ifges(7,5) * t364;
t38 = Ifges(7,2) * t428 + Ifges(7,6) * t221 + t355;
t385 = -t38 / 0.2e1;
t384 = t38 / 0.2e1;
t39 = Ifges(7,1) * t86 + Ifges(7,5) * t221 + t81;
t383 = -t39 / 0.2e1;
t382 = t39 / 0.2e1;
t381 = Ifges(6,4) * t373 + Ifges(6,2) * t372 + Ifges(6,6) * t364;
t380 = Ifges(6,1) * t373 + Ifges(6,4) * t372 + Ifges(6,5) * t364;
t337 = Ifges(6,4) * t147;
t71 = Ifges(6,2) * t420 + Ifges(6,6) * t230 + t337;
t379 = t71 / 0.2e1;
t144 = Ifges(6,4) * t420;
t72 = Ifges(6,1) * t147 + Ifges(6,5) * t230 + t144;
t378 = t72 / 0.2e1;
t377 = -t428 / 0.2e1;
t376 = t428 / 0.2e1;
t374 = t86 / 0.2e1;
t370 = -t420 / 0.2e1;
t369 = t420 / 0.2e1;
t368 = -t147 / 0.2e1;
t367 = t147 / 0.2e1;
t365 = t221 / 0.2e1;
t363 = -t230 / 0.2e1;
t362 = t230 / 0.2e1;
t360 = t233 / 0.2e1;
t359 = -t235 / 0.2e1;
t352 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t310);
t347 = t428 * Ifges(7,6);
t346 = t86 * Ifges(7,5);
t343 = mrSges(4,3) * t233;
t335 = t420 * Ifges(6,6);
t334 = t147 * Ifges(6,5);
t331 = t221 * Ifges(7,3);
t330 = t230 * Ifges(6,3);
t329 = t235 * Ifges(4,1);
t323 = Ifges(3,5) * qJD(2);
t158 = t213 * t269 - t324 * t277;
t198 = -t324 * t258 - t259 * t269;
t320 = t158 * t198;
t314 = t236 * t268;
t150 = pkin(3) * t234 - qJ(4) * t236 - qJD(4) * t249 + t417;
t177 = t231 * t324 + t269 * t232;
t102 = t268 * t150 + t270 * t177;
t162 = mrSges(5,1) * t318 + mrSges(5,2) * t317;
t301 = Ifges(6,5) * t98 + Ifges(6,6) * t99 + Ifges(6,3) * t222;
t299 = t324 * pkin(2);
t53 = -t99 * mrSges(6,1) + t98 * mrSges(6,2);
t9 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t295 = t276 * t303;
t290 = t222 * mrSges(4,1) + t223 * mrSges(4,2);
t54 = t275 * t106 - t114 * t272;
t101 = t270 * t150 - t177 * t268;
t176 = t231 * t269 - t324 * t232;
t264 = -t299 - pkin(3);
t128 = pkin(4) * t318 + t158;
t140 = pkin(4) * t314 + t176;
t173 = pkin(4) * t313 + t198;
t283 = Ifges(5,5) * t270 - Ifges(5,6) * t268;
t282 = t268 * t78 + t270 * t77;
t180 = t280 * t249;
t46 = -pkin(5) * t278 + pkin(9) * t180 + t54;
t179 = t250 * t249;
t50 = -pkin(9) * t179 + t55;
t21 = -t271 * t50 + t274 * t46;
t22 = t271 * t46 + t274 * t50;
t281 = t111 * t268 - t112 * t270;
t117 = -t179 * t274 + t180 * t271;
t118 = -t179 * t271 - t180 * t274;
t67 = pkin(4) * t234 - t236 * t350 + t101;
t79 = -pkin(8) * t314 + t102;
t23 = t106 * t304 - t114 * t305 + t272 * t67 + t275 * t79;
t254 = -t270 * pkin(4) + t264;
t24 = -qJD(5) * t55 - t272 * t79 + t275 * t67;
t266 = Ifges(3,4) * t309;
t243 = Ifges(3,1) * t310 + t266 + t323;
t229 = Ifges(4,4) * t233;
t211 = -qJD(2) * mrSges(4,2) - t343;
t200 = pkin(5) * t280 + t254;
t188 = mrSges(4,1) * t233 + t413;
t182 = -t229 + t329 + t406;
t172 = mrSges(5,1) * t222 - mrSges(5,3) * t317;
t171 = -mrSges(5,2) * t222 - mrSges(5,3) * t318;
t165 = mrSges(5,1) * t233 - mrSges(5,3) * t203;
t164 = -mrSges(5,2) * t233 + mrSges(5,3) * t288;
t137 = t222 * Ifges(5,5) + t223 * t285;
t136 = t222 * Ifges(5,6) + t223 * t284;
t129 = t233 * Ifges(5,3) + t332 + t333;
t123 = mrSges(6,1) * t230 - mrSges(6,3) * t147;
t122 = -mrSges(6,2) * t230 + mrSges(6,3) * t420;
t119 = pkin(5) * t179 + t173;
t116 = -t236 * t250 + t237 * t249;
t115 = -t236 * t280 - t238 * t249;
t90 = -mrSges(6,1) * t420 + mrSges(6,2) * t147;
t75 = -mrSges(6,2) * t222 + mrSges(6,3) * t99;
t74 = mrSges(6,1) * t222 - mrSges(6,3) * t98;
t70 = t330 + t334 + t335;
t66 = mrSges(7,1) * t221 - mrSges(7,3) * t86;
t65 = -mrSges(7,2) * t221 + mrSges(7,3) * t428;
t62 = -pkin(5) * t116 + t140;
t56 = -pkin(5) * t99 + t128;
t47 = -mrSges(7,1) * t428 + mrSges(7,2) * t86;
t41 = -qJD(6) * t118 - t115 * t271 + t116 * t274;
t40 = qJD(6) * t117 + t115 * t274 + t116 * t271;
t37 = t331 + t346 + t347;
t26 = -mrSges(7,2) * t222 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t222 - mrSges(7,3) * t30;
t20 = pkin(9) * t116 + t23;
t17 = pkin(5) * t234 - pkin(9) * t115 + t24;
t14 = t274 * t33 - t326;
t13 = -t271 * t33 - t325;
t5 = -qJD(6) * t22 + t17 * t274 - t20 * t271;
t4 = qJD(6) * t21 + t17 * t271 + t20 * t274;
t1 = [(-t191 * t236 - t192 * t234 + t198 * t223 - t199 * t222) * mrSges(4,3) + (t329 / 0.2e1 + t283 * t360 + t182 / 0.2e1 + (-t111 * t270 - t112 * t268) * mrSges(5,3) + t393) * t236 + (-mrSges(5,3) * t282 + t136 * t357 + t137 * t356 + (Ifges(4,1) + Ifges(5,1) * t270 ^ 2 / 0.2e1 + (-t338 + t336 / 0.2e1) * t268) * t223 - t222 * Ifges(4,4) + t283 * t364 + (t286 + mrSges(4,3)) * t158) * t249 + t265 * t290 + (-t322 / 0.2e1 + (mrSges(3,1) * t391 - 0.3e1 / 0.2e1 * t341 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t276) * qJD(1) + (t413 + t188 + m(4) * (t255 + t311)) * pkin(2) - t394) * t308 + (t234 * t359 + t236 * t361) * Ifges(4,4) + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t386 + m(6) * (t128 * t173 + t139 * t140 + t15 * t55 + t16 * t54 + t23 * t45 + t24 * t44) + m(7) * (t11 * t5 + t119 * t56 + t12 * t4 + t2 * t22 + t21 * t3 + t62 * t80) + (Ifges(7,1) * t118 + Ifges(7,4) * t117) * t387 + (-Ifges(6,5) * t180 + Ifges(7,5) * t118 - Ifges(6,6) * t179 + Ifges(7,6) * t117) * t364 + (-(Ifges(5,3) + Ifges(4,2)) * t222 - mrSges(4,1) * qJD(1) * t417 - mrSges(5,1) * t77 + mrSges(5,2) * t78 + mrSges(4,3) * t159 - Ifges(6,5) * t373 - Ifges(7,5) * t387 - Ifges(6,6) * t372 - Ifges(7,6) * t386 - (Ifges(6,3) + Ifges(7,3)) * t364 - (-Ifges(4,4) + t283) * t223 - t395 - t416 - t302 / 0.2e1 - t301 / 0.2e1) * t278 + (-t11 * t40 + t117 * t2 - t118 * t3 + t12 * t41) * mrSges(7,3) + t177 * t211 + t198 * t162 + t132 * t172 + t173 * t53 + t133 * t171 + t102 * t164 + t101 * t165 + t139 * (-mrSges(6,1) * t116 + mrSges(6,2) * t115) + t140 * t90 + t23 * t122 + t24 * t123 + t56 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) + t119 * t9 + t80 * (-mrSges(7,1) * t41 + mrSges(7,2) * t40) + t54 * t74 + t55 * t75 + t4 * t65 + t5 * t66 + t62 * t47 + t21 * t25 + t22 * t26 + t118 * t388 + t117 * t389 + (Ifges(7,4) * t40 + Ifges(7,2) * t41) * t376 + t115 * t378 + t116 * t379 - t180 * t380 - t179 * t381 + t40 * t382 + t41 * t384 + (Ifges(7,5) * t40 + Ifges(7,6) * t41) * t365 + (Ifges(6,1) * t115 + Ifges(6,4) * t116) * t367 + (Ifges(6,4) * t115 + Ifges(6,2) * t116) * t369 + (Ifges(7,1) * t40 + Ifges(7,4) * t41) * t374 + (Ifges(6,5) * t115 + Ifges(6,6) * t116) * t362 + (t335 / 0.2e1 + t346 / 0.2e1 + t347 / 0.2e1 + t129 / 0.2e1 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t233 + t37 / 0.2e1 + t70 / 0.2e1 + t334 / 0.2e1 + t331 / 0.2e1 + t330 / 0.2e1 - t419) * t234 + m(4) * (t159 * t199 - t176 * t191 + t177 * t192 + t320) + m(5) * (t101 * t111 + t102 * t112 + t132 * t77 + t133 * t78 + t176 * t178 + t320) - t399 * t176 + (Ifges(4,5) * t236 / 0.2e1 - Ifges(4,6) * t234 / 0.2e1 + (-t352 + t243 / 0.2e1 + t323 / 0.2e1 + (mrSges(3,2) * t391 + 0.3e1 / 0.2e1 * Ifges(3,4) * t276) * qJD(1)) * t276) * qJD(2) + (-Ifges(6,1) * t180 - Ifges(6,4) * t179) * t373 + (-Ifges(6,4) * t180 - Ifges(6,2) * t179) * t372 + (-t115 * t44 + t116 * t45 - t15 * t179 + t16 * t180) * mrSges(6,3) + t128 * (mrSges(6,1) * t179 - mrSges(6,2) * t180); t414 * t66 + (t11 * t414 + t12 * t415 + t2 * t94 + t200 * t56 + t3 * t93 + t402 * t80) * m(7) + t415 * t65 + (Ifges(6,4) * t170 + Ifges(6,2) * t169) * t370 - Ifges(3,6) * t296 + (Ifges(6,1) * t170 + Ifges(6,4) * t169) * t368 + (-mrSges(7,1) * t400 + mrSges(7,2) * t401) * t80 + (-t11 * t401 + t12 * t400 + t193 * t2 - t194 * t3) * mrSges(7,3) + t402 * t47 + (-Ifges(6,1) * t237 - Ifges(6,4) * t238) * t367 + (-Ifges(6,4) * t237 - Ifges(6,2) * t238) * t369 + (-Ifges(6,5) * t237 - Ifges(6,6) * t238) * t362 + (Ifges(7,5) * t110 + Ifges(7,6) * t109) * t366 + (-t121 + t306) * t164 + (Ifges(6,5) * t170 + Ifges(6,6) * t169) * t363 + (Ifges(7,1) * t110 + Ifges(7,4) * t109) * t375 + (Ifges(7,4) * t110 + Ifges(7,2) * t109) * t377 + (t171 * t270 - t172 * t268) * t262 + (-mrSges(3,1) * t295 + mrSges(3,2) * t296) * pkin(7) + t394 * t310 + (t406 / 0.2e1 - t283 * t361 + t393) * t233 + t407 * t123 + (t128 * t254 - t139 * t155 + t15 * t187 + t16 * t186 + t407 * t44 + t408 * t45) * m(6) + t408 * t122 + (-t111 * t315 - t112 * t316 + t396) * mrSges(5,3) + (-t281 * qJD(4) - t111 * t120 - t112 * t121 + t158 * t264 - t178 * t195 + t262 * t396) * m(5) + (-mrSges(6,1) * t397 + mrSges(6,2) * t398) * t139 + t399 * t195 + ((-t158 * t324 + t159 * t269) * pkin(2) + t191 * t195 - t192 * t196 - t255 * t300) * m(4) + (-Ifges(4,1) * t233 + t129 - t340 + t37 + t70) * t359 + (Ifges(5,5) * t268 + Ifges(6,5) * t250 + Ifges(7,5) * t194 + Ifges(5,6) * t270 - Ifges(6,6) * t280 + Ifges(7,6) * t193) * t364 + t128 * (mrSges(6,1) * t280 + mrSges(6,2) * t250) + (Ifges(6,4) * t250 - Ifges(6,2) * t280) * t372 + (Ifges(6,1) * t250 - Ifges(6,4) * t280) * t373 + (-t15 * t280 - t16 * t250 + t397 * t45 - t398 * t44) * mrSges(6,3) - t280 * t381 + t268 * t137 / 0.2e1 + t264 * t162 + t254 * t53 - Ifges(4,6) * t222 + Ifges(4,5) * t223 - t196 * t211 + t200 * t9 + t56 * (-mrSges(7,1) * t193 + mrSges(7,2) * t194) + t186 * t74 + t187 * t75 - t170 * t72 / 0.2e1 - t169 * t71 / 0.2e1 - t159 * mrSges(4,2) - t155 * t90 + t93 * t25 + t94 * t26 + (-t222 * t354 - t223 * t299) * mrSges(4,3) + (-mrSges(5,1) * t270 + mrSges(5,2) * t268 - mrSges(4,1)) * t158 + (Ifges(7,1) * t194 + Ifges(7,4) * t193) * t387 + t194 * t388 + t193 * t389 + (Ifges(7,4) * t124 + Ifges(7,2) * t125) * t376 - t237 * t378 - t238 * t379 + t250 * t380 + t124 * t382 + t110 * t383 + t125 * t384 + t109 * t385 + (Ifges(7,4) * t194 + Ifges(7,2) * t193) * t386 + (Ifges(7,5) * t124 + Ifges(7,6) * t125) * t365 + (Ifges(7,1) * t124 + Ifges(7,4) * t125) * t374 + t309 * t352 + t136 * t356 + (-t120 - t307) * t165 - t188 * t300 - (Ifges(3,5) * t276 - Ifges(3,6) * t273) * t303 / 0.2e1 + Ifges(3,5) * t295 + t192 * t342 + (Ifges(5,1) * t268 + t338) * t317 / 0.2e1 + (Ifges(6,5) * t368 + Ifges(7,5) * t375 - Ifges(4,2) * t360 + Ifges(6,6) * t370 + Ifges(7,6) * t377 + Ifges(5,3) * t361 + Ifges(6,3) * t363 + Ifges(7,3) * t366 + t419 + t422) * t235 - (Ifges(5,2) * t270 + t339) * t318 / 0.2e1 - t191 * t343 + (-t229 + t182) * t360 - (-Ifges(3,2) * t310 + t243 + t266) * t309 / 0.2e1 + (-t273 * (Ifges(3,1) * t276 - t341) / 0.2e1 + pkin(1) * (mrSges(3,1) * t273 + mrSges(3,2) * t276)) * qJD(1) ^ 2; (-t90 - t47 + t399) * t235 + (t164 * t270 - t165 * t268 + t211) * t233 + t400 * t66 + t401 * t65 + t397 * t123 + t398 * t122 + t290 + t270 * t172 + t268 * t171 + t250 * t75 - t280 * t74 + t194 * t26 + t193 * t25 + (t11 * t400 + t12 * t401 + t193 * t3 + t194 * t2 - t235 * t80) * m(7) + (-t139 * t235 + t15 * t250 - t16 * t280 + t397 * t44 + t398 * t45) * m(6) + (-t178 * t235 - t233 * t281 + t282) * m(5) + (t191 * t235 + t192 * t233 + t287) * m(4); -t420 * t122 + t147 * t123 - t288 * t164 + t203 * t165 - t428 * t65 + t86 * t66 + t162 + t53 + t9 + (t11 * t86 - t12 * t428 + t56) * m(7) + (t147 * t44 - t420 * t45 + t128) * m(6) + (t111 * t203 - t112 * t288 + t158) * m(5); t301 + t395 - m(7) * (t11 * t13 + t12 * t14) + (t147 * t45 + t420 * t44) * mrSges(6,3) - t86 * t385 - t139 * (mrSges(6,1) * t147 + mrSges(6,2) * t420) + t45 * t123 - t44 * t122 - t14 * t65 - t13 * t66 + t421 * t377 + t71 * t367 + (Ifges(6,1) * t420 - t337) * t368 + (Ifges(6,5) * t420 - Ifges(6,6) * t147) * t363 + (-t147 * t47 + t274 * t25 + t271 * t26 + (-t271 * t66 + t274 * t65) * qJD(6) + (-t147 * t80 + t2 * t271 + t274 * t3 + (-t11 * t271 + t12 * t274) * qJD(6)) * m(7)) * pkin(5) + (-Ifges(6,2) * t147 + t144 + t72) * t370 + t428 * t383 + t423; t38 * t374 - t11 * t65 + t12 * t66 + (t39 + t421) * t377 + t423;];
tauc  = t1(:);
