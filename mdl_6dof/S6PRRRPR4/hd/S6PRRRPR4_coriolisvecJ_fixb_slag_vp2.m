% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:24:12
% EndTime: 2018-11-23 15:24:26
% DurationCPUTime: 14.71s
% Computational Cost: add. (10552->689), mult. (26702->975), div. (0->0), fcn. (19745->12), ass. (0->318)
t277 = sin(qJ(2));
t271 = sin(pkin(6));
t332 = qJD(1) * t271;
t317 = t277 * t332;
t245 = qJD(2) * pkin(8) + t317;
t276 = sin(qJ(3));
t280 = cos(qJ(3));
t273 = cos(pkin(6));
t331 = qJD(1) * t273;
t196 = -t245 * t276 + t280 * t331;
t303 = pkin(3) * t276 - pkin(9) * t280;
t242 = t303 * qJD(2);
t275 = sin(qJ(4));
t279 = cos(qJ(4));
t142 = -t196 * t275 + t242 * t279;
t339 = t279 * t280;
t292 = pkin(4) * t276 - qJ(5) * t339;
t364 = -qJ(5) - pkin(9);
t305 = qJD(4) * t364;
t442 = -qJD(2) * t292 - qJD(5) * t275 + t279 * t305 - t142;
t143 = t196 * t279 + t242 * t275;
t328 = qJD(2) * t280;
t313 = t275 * t328;
t321 = qJD(5) * t279;
t441 = -qJ(5) * t313 - t275 * t305 + t143 - t321;
t270 = sin(pkin(12));
t272 = cos(pkin(12));
t233 = t270 * t279 + t272 * t275;
t289 = t233 * t280;
t200 = qJD(2) * t289;
t217 = t233 * qJD(4);
t440 = t200 - t217;
t293 = t270 * t275 - t272 * t279;
t288 = t293 * t280;
t201 = qJD(2) * t288;
t218 = t293 * qJD(4);
t336 = -t201 + t218;
t421 = t270 * t441 + t272 * t442;
t420 = t270 * t442 - t272 * t441;
t261 = qJD(4) - t328;
t326 = qJD(3) * t279;
t330 = qJD(2) * t276;
t239 = -t275 * t330 + t326;
t240 = qJD(3) * t275 + t279 * t330;
t166 = t239 * t270 + t240 * t272;
t425 = pkin(10) * t166;
t197 = t245 * t280 + t276 * t331;
t184 = qJD(3) * pkin(9) + t197;
t247 = -pkin(3) * t280 - pkin(9) * t276 - pkin(2);
t281 = cos(qJ(2));
t316 = t281 * t332;
t199 = qJD(2) * t247 - t316;
t126 = -t184 * t275 + t199 * t279;
t103 = -qJ(5) * t240 + t126;
t84 = pkin(4) * t261 + t103;
t128 = t184 * t279 + t199 * t275;
t104 = qJ(5) * t239 + t128;
t97 = t270 * t104;
t46 = t272 * t84 - t97;
t26 = pkin(5) * t261 - t425 + t46;
t274 = sin(qJ(6));
t278 = cos(qJ(6));
t304 = t239 * t272 - t240 * t270;
t412 = pkin(10) * t304;
t342 = t272 * t104;
t47 = t270 * t84 + t342;
t33 = t47 + t412;
t12 = t26 * t278 - t274 * t33;
t13 = t26 * t274 + t278 * t33;
t320 = qJD(2) * qJD(3);
t308 = t276 * t320;
t258 = Ifges(7,3) * t308;
t94 = t166 * t278 + t274 * t304;
t373 = Ifges(7,4) * t94;
t252 = qJD(6) + t261;
t379 = -t252 / 0.2e1;
t393 = -t94 / 0.2e1;
t415 = -t166 * t274 + t278 * t304;
t395 = -t415 / 0.2e1;
t87 = Ifges(7,4) * t415;
t45 = Ifges(7,1) * t94 + Ifges(7,5) * t252 + t87;
t183 = -qJD(3) * pkin(3) - t196;
t148 = -pkin(4) * t239 + qJD(5) + t183;
t86 = -pkin(5) * t304 + t148;
t439 = t258 + (Ifges(7,5) * t415 - Ifges(7,6) * t94) * t379 + (t12 * t415 + t13 * t94) * mrSges(7,3) + (-Ifges(7,2) * t94 + t45 + t87) * t395 - t86 * (mrSges(7,1) * t94 + mrSges(7,2) * t415) + (Ifges(7,1) * t415 - t373) * t393;
t438 = -pkin(5) * t330 + pkin(10) * t336 + t421;
t437 = pkin(10) * t440 + t420;
t44 = Ifges(7,2) * t415 + Ifges(7,6) * t252 + t373;
t435 = t44 / 0.2e1;
t434 = -t328 / 0.2e1;
t338 = t280 * t281;
t190 = (-t275 * t338 + t277 * t279) * t332;
t191 = (t275 * t277 + t279 * t338) * t332;
t125 = t190 * t272 - t191 * t270;
t243 = t303 * qJD(3);
t322 = qJD(4) * t279;
t335 = t243 * t275 + t247 * t322;
t340 = t276 * t279;
t101 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t340 + (-qJD(5) * t276 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t280) * t275 + t335;
t263 = pkin(8) * t339;
t327 = qJD(3) * t276;
t370 = pkin(8) * t275;
t334 = t243 * t279 + t327 * t370;
t85 = -t276 * t321 + t292 * qJD(3) + (-t263 + (qJ(5) * t276 - t247) * t275) * qJD(4) + t334;
t48 = -t101 * t270 + t272 * t85;
t429 = -t125 + t48;
t127 = t190 * t270 + t191 * t272;
t50 = t101 * t272 + t270 * t85;
t428 = -t127 + t50;
t427 = -Ifges(4,1) / 0.2e1;
t426 = Ifges(4,4) * t434;
t250 = t364 * t275;
t251 = t364 * t279;
t180 = t250 * t272 + t251 * t270;
t149 = -pkin(10) * t233 + t180;
t181 = t250 * t270 - t251 * t272;
t150 = -pkin(10) * t293 + t181;
t72 = t149 * t278 - t150 * t274;
t424 = qJD(6) * t72 + t274 * t438 + t278 * t437;
t73 = t149 * t274 + t150 * t278;
t423 = -qJD(6) * t73 - t274 * t437 + t278 * t438;
t422 = t166 * Ifges(6,4);
t170 = pkin(4) * t313 + t197;
t324 = qJD(4) * t275;
t419 = pkin(4) * t324 - pkin(5) * t440 - t170;
t325 = qJD(3) * t280;
t418 = t275 * t325 + t276 * t322;
t319 = qJD(3) * qJD(4);
t323 = qJD(4) * t276;
t188 = t279 * t319 + (-t275 * t323 + t279 * t325) * qJD(2);
t189 = -qJD(2) * t418 - t275 * t319;
t123 = t188 * t272 + t189 * t270;
t343 = t271 * t281;
t314 = qJD(2) * t343;
t287 = qJD(1) * (qJD(3) * t273 + t314);
t155 = -t245 * t327 + t280 * t287;
t193 = (t243 + t317) * qJD(2);
t62 = -qJD(4) * t128 - t155 * t275 + t193 * t279;
t40 = pkin(4) * t308 - qJ(5) * t188 - qJD(5) * t240 + t62;
t61 = t155 * t279 - t184 * t324 + t193 * t275 + t199 * t322;
t42 = qJ(5) * t189 + qJD(5) * t239 + t61;
t17 = -t270 * t42 + t272 * t40;
t8 = pkin(5) * t308 - pkin(10) * t123 + t17;
t122 = -t188 * t270 + t189 * t272;
t18 = t270 * t40 + t272 * t42;
t9 = pkin(10) * t122 + t18;
t2 = qJD(6) * t12 + t274 * t8 + t278 * t9;
t3 = -qJD(6) * t13 - t274 * t9 + t278 * t8;
t31 = qJD(6) * t415 + t122 * t274 + t123 * t278;
t32 = -qJD(6) * t94 + t122 * t278 - t123 * t274;
t417 = mrSges(7,1) * t3 - t2 * mrSges(7,2) + Ifges(7,5) * t31 + Ifges(7,6) * t32;
t246 = -qJD(2) * pkin(2) - t316;
t356 = t240 * Ifges(5,4);
t153 = Ifges(5,2) * t239 + Ifges(5,6) * t261 + t356;
t231 = Ifges(5,4) * t239;
t154 = Ifges(5,1) * t240 + Ifges(5,5) * t261 + t231;
t294 = t126 * t279 + t128 * t275;
t296 = Ifges(5,5) * t279 - Ifges(5,6) * t275;
t361 = Ifges(5,4) * t279;
t298 = -Ifges(5,2) * t275 + t361;
t362 = Ifges(5,4) * t275;
t300 = Ifges(5,1) * t279 - t362;
t301 = mrSges(5,1) * t275 + mrSges(5,2) * t279;
t374 = t279 / 0.2e1;
t375 = -t275 / 0.2e1;
t376 = t261 / 0.2e1;
t380 = t240 / 0.2e1;
t284 = -t294 * mrSges(5,3) + t296 * t376 + t239 * t298 / 0.2e1 + t300 * t380 + t183 * t301 + t153 * t375 + t154 * t374;
t347 = Ifges(4,5) * qJD(3);
t416 = t196 * mrSges(4,3) + t330 * t427 + t426 - t347 / 0.2e1 - t246 * mrSges(4,2) - t284;
t399 = t31 / 0.2e1;
t398 = t32 / 0.2e1;
t82 = Ifges(6,2) * t304 + Ifges(6,6) * t261 + t422;
t414 = t82 / 0.2e1;
t390 = t122 / 0.2e1;
t389 = t123 / 0.2e1;
t377 = -t261 / 0.2e1;
t388 = -t304 / 0.2e1;
t413 = t308 / 0.2e1;
t366 = qJD(2) / 0.2e1;
t14 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t68 = -mrSges(6,1) * t122 + mrSges(6,2) * t123;
t411 = t14 + t68;
t410 = Ifges(6,4) * t304;
t265 = pkin(4) * t272 + pkin(5);
t371 = pkin(4) * t270;
t212 = t265 * t278 - t274 * t371;
t53 = -t103 * t270 - t342;
t36 = t53 - t412;
t54 = t103 * t272 - t97;
t37 = t54 - t425;
t409 = qJD(6) * t212 - t274 * t36 - t278 * t37;
t213 = t265 * t274 + t278 * t371;
t408 = -qJD(6) * t213 + t274 * t37 - t278 * t36;
t203 = t247 * t275 + t263;
t405 = -t275 * t62 + t279 * t61;
t402 = -t62 * mrSges(5,1) - t17 * mrSges(6,1) + t61 * mrSges(5,2) + t18 * mrSges(6,2) - Ifges(5,5) * t188 - Ifges(6,5) * t123 - Ifges(5,6) * t189 - Ifges(6,6) * t122 - t417;
t401 = Ifges(7,4) * t399 + Ifges(7,2) * t398 + Ifges(7,6) * t413;
t400 = Ifges(7,1) * t399 + Ifges(7,4) * t398 + Ifges(7,5) * t413;
t397 = Ifges(6,4) * t389 + Ifges(6,2) * t390 + Ifges(6,6) * t413;
t396 = Ifges(6,1) * t389 + Ifges(6,4) * t390 + Ifges(6,5) * t413;
t394 = t415 / 0.2e1;
t392 = t94 / 0.2e1;
t387 = t304 / 0.2e1;
t386 = -t166 / 0.2e1;
t385 = t166 / 0.2e1;
t384 = t188 / 0.2e1;
t383 = t189 / 0.2e1;
t382 = -t239 / 0.2e1;
t381 = -t240 / 0.2e1;
t378 = t252 / 0.2e1;
t372 = pkin(4) * t240;
t369 = pkin(8) * t280;
t363 = Ifges(4,4) * t276;
t352 = t261 * Ifges(6,3);
t134 = -t200 * t278 + t201 * t274;
t160 = t233 * t278 - t274 * t293;
t89 = -qJD(6) * t160 - t217 * t278 + t218 * t274;
t349 = t134 - t89;
t135 = -t200 * t274 - t201 * t278;
t159 = -t233 * t274 - t278 * t293;
t88 = qJD(6) * t159 - t217 * t274 - t218 * t278;
t348 = t135 - t88;
t346 = Ifges(4,6) * qJD(3);
t156 = t245 * t325 + t276 * t287;
t344 = t271 * t277;
t219 = -t273 * t280 + t276 * t344;
t345 = t156 * t219;
t341 = t275 * t276;
t235 = t279 * t247;
t167 = -qJ(5) * t340 + t235 + (-pkin(4) - t370) * t280;
t176 = -qJ(5) * t341 + t203;
t106 = t167 * t270 + t176 * t272;
t333 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t239 - mrSges(5,2) * t240 - mrSges(4,3) * t330;
t244 = pkin(4) * t341 + pkin(8) * t276;
t329 = qJD(2) * t277;
t318 = mrSges(4,3) * t328;
t198 = pkin(4) * t418 + pkin(8) * t325;
t266 = -pkin(4) * t279 - pkin(3);
t315 = t271 * t329;
t310 = t347 / 0.2e1;
t309 = -t346 / 0.2e1;
t307 = 0.3e1 / 0.2e1 * t328;
t102 = -mrSges(6,1) * t304 + mrSges(6,2) * t166;
t306 = m(6) * t148 + t102;
t105 = t167 * t272 - t176 * t270;
t302 = mrSges(5,1) * t279 - mrSges(5,2) * t275;
t299 = Ifges(5,1) * t275 + t361;
t297 = Ifges(5,2) * t279 + t362;
t295 = Ifges(5,5) * t275 + Ifges(5,6) * t279;
t211 = t293 * t276;
t71 = -pkin(5) * t280 + pkin(10) * t211 + t105;
t210 = t233 * t276;
t74 = -pkin(10) * t210 + t106;
t34 = -t274 * t74 + t278 * t71;
t35 = t274 * t71 + t278 * t74;
t220 = t273 * t276 + t280 * t344;
t174 = -t220 * t275 - t279 * t343;
t291 = -t220 * t279 + t275 * t343;
t110 = t174 * t272 + t270 * t291;
t111 = t174 * t270 - t272 * t291;
t56 = t110 * t278 - t111 * t274;
t57 = t110 * t274 + t111 * t278;
t140 = -t210 * t278 + t211 * t274;
t141 = -t210 * t274 - t211 * t278;
t290 = -m(4) * t196 + m(5) * t183 - t333;
t112 = -pkin(4) * t189 + t156;
t283 = t128 * mrSges(5,2) + t13 * mrSges(7,2) + t197 * mrSges(4,3) + t47 * mrSges(6,2) + Ifges(5,3) * t377 - t240 * Ifges(5,5) - t239 * Ifges(5,6) + t346 / 0.2e1 + (t280 * Ifges(4,2) + t363) * t366 - t252 * Ifges(7,3) - t94 * Ifges(7,5) - t415 * Ifges(7,6) - t352 / 0.2e1 - t166 * Ifges(6,5) - t304 * Ifges(6,6) - t12 * mrSges(7,1) - t126 * mrSges(5,1) - t246 * mrSges(4,1) - t46 * mrSges(6,1);
t282 = qJD(2) ^ 2;
t260 = Ifges(5,3) * t308;
t259 = Ifges(6,3) * t308;
t249 = -qJD(3) * mrSges(4,2) + t318;
t241 = (-mrSges(4,1) * t280 + mrSges(4,2) * t276) * qJD(2);
t227 = (mrSges(4,1) * t276 + mrSges(4,2) * t280) * t320;
t204 = pkin(5) * t293 + t266;
t202 = -t275 * t369 + t235;
t195 = mrSges(5,1) * t261 - mrSges(5,3) * t240;
t194 = -mrSges(5,2) * t261 + mrSges(5,3) * t239;
t173 = qJD(3) * t220 + t276 * t314;
t172 = -qJD(3) * t219 + t280 * t314;
t168 = pkin(5) * t210 + t244;
t162 = -mrSges(5,2) * t308 + mrSges(5,3) * t189;
t161 = mrSges(5,1) * t308 - mrSges(5,3) * t188;
t145 = -qJD(3) * t288 - t217 * t276;
t144 = -qJD(3) * t289 + t293 * t323;
t139 = mrSges(6,1) * t261 - mrSges(6,3) * t166;
t138 = -mrSges(6,2) * t261 + mrSges(6,3) * t304;
t137 = -qJD(4) * t203 + t334;
t136 = (-t276 * t326 - t280 * t324) * pkin(8) + t335;
t133 = pkin(5) * t166 + t372;
t131 = -mrSges(5,1) * t189 + mrSges(5,2) * t188;
t116 = Ifges(5,1) * t188 + Ifges(5,4) * t189 + Ifges(5,5) * t308;
t115 = Ifges(5,4) * t188 + Ifges(5,2) * t189 + Ifges(5,6) * t308;
t109 = mrSges(6,1) * t308 - mrSges(6,3) * t123;
t108 = -mrSges(6,2) * t308 + mrSges(6,3) * t122;
t107 = -pkin(5) * t144 + t198;
t96 = qJD(4) * t291 - t172 * t275 + t279 * t315;
t95 = qJD(4) * t174 + t172 * t279 + t275 * t315;
t83 = Ifges(6,1) * t166 + Ifges(6,5) * t261 + t410;
t76 = mrSges(7,1) * t252 - mrSges(7,3) * t94;
t75 = -mrSges(7,2) * t252 + mrSges(7,3) * t415;
t67 = t125 * t274 + t127 * t278;
t66 = t125 * t278 - t127 * t274;
t65 = -pkin(5) * t122 + t112;
t60 = -qJD(6) * t141 + t144 * t278 - t145 * t274;
t59 = qJD(6) * t140 + t144 * t274 + t145 * t278;
t52 = -mrSges(7,1) * t415 + mrSges(7,2) * t94;
t51 = t270 * t96 + t272 * t95;
t49 = -t270 * t95 + t272 * t96;
t27 = pkin(10) * t144 + t50;
t25 = pkin(5) * t327 - pkin(10) * t145 + t48;
t24 = -mrSges(7,2) * t308 + mrSges(7,3) * t32;
t23 = mrSges(7,1) * t308 - mrSges(7,3) * t31;
t7 = -qJD(6) * t57 - t274 * t51 + t278 * t49;
t6 = qJD(6) * t56 + t274 * t49 + t278 * t51;
t5 = -qJD(6) * t35 + t25 * t278 - t27 * t274;
t4 = qJD(6) * t34 + t25 * t274 + t27 * t278;
t1 = [t172 * t249 + t95 * t194 + t96 * t195 - t291 * t162 + t174 * t161 + t51 * t138 + t49 * t139 + t110 * t109 + t111 * t108 + t6 * t75 + t7 * t76 + t56 * t23 + t57 * t24 - t220 * mrSges(4,3) * t308 + ((-mrSges(3,2) * t282 - t227) * t281 + (-mrSges(3,1) * t282 + qJD(2) * t241) * t277) * t271 + (qJD(3) * t318 + t131 + t411) * t219 + (t102 + t52 - t333) * t173 + m(4) * (t155 * t220 + t345 + t172 * t197 - t173 * t196 + (t246 - t316) * t315) + m(7) * (t12 * t7 + t13 * t6 + t173 * t86 + t2 * t57 + t219 * t65 + t3 * t56) + m(6) * (t110 * t17 + t111 * t18 + t112 * t219 + t148 * t173 + t46 * t49 + t47 * t51) + m(5) * (t126 * t96 + t128 * t95 + t173 * t183 + t174 * t62 - t291 * t61 + t345); (t4 - t67) * t75 + (-t190 + t137) * t195 + (-t249 * t338 + (-mrSges(4,1) * t328 - t241) * t277) * t332 + (Ifges(4,4) * t307 + t290 * pkin(8) + t310 - t416) * t325 - m(4) * (t197 * t280 * t316 + t246 * t317) + m(5) * (t126 * t137 + t128 * t136 + t202 * t62 + t203 * t61) - m(7) * (t12 * t66 + t13 * t67) - m(5) * (t126 * t190 + t128 * t191) + (-t258 / 0.2e1 - t259 / 0.2e1 - t260 / 0.2e1 + t155 * mrSges(4,3) + t402) * t280 + m(7) * (t107 * t86 + t12 * t5 + t13 * t4 + t168 * t65 + t2 * t35 + t3 * t34) + m(4) * (-pkin(2) * qJD(1) * t315 + t155 * t369) + t428 * t138 + t429 * t139 + (t105 * t17 + t106 * t18 + t112 * t244 + t148 * t198 + t428 * t47 + t429 * t46) * m(6) + (-t66 + t5) * t76 + (t115 * t375 + t116 * t374 + t300 * t384 + t298 * t383 + (-t275 * t61 - t279 * t62) * mrSges(5,3) + (mrSges(4,3) + t301) * t156 + (mrSges(4,2) * t329 + (-m(7) * t86 - t290 - t306 - t52) * t281) * t332 + (t295 * t377 + t297 * t382 + t299 * t381 + t183 * t302 + t154 * t375 - t279 * t153 / 0.2e1 + (t126 * t275 - t128 * t279) * mrSges(5,3)) * qJD(4) + (Ifges(4,1) * t307 - 0.3e1 / 0.2e1 * Ifges(4,2) * t328 + (t434 + t376) * Ifges(5,3) + t352 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4) * t330 + t309 - t283 + (-Ifges(6,5) * t211 + Ifges(7,5) * t141 - Ifges(6,6) * t210 + Ifges(7,6) * t140 + t296 * t276 + (-Ifges(6,3) - Ifges(7,3)) * t280) * t366) * qJD(3) + (t131 + (m(4) + m(5)) * t156 + (-m(4) * t197 - t249) * qJD(3)) * pkin(8)) * t276 + (-t191 + t136) * t194 + (-t12 * t59 + t13 * t60 + t140 * t2 - t141 * t3) * mrSges(7,3) + t244 * t68 - pkin(2) * t227 + t202 * t161 + t203 * t162 + t198 * t102 + t168 * t14 + t145 * t83 / 0.2e1 + t148 * (-mrSges(6,1) * t144 + mrSges(6,2) * t145) + t65 * (-mrSges(7,1) * t140 + mrSges(7,2) * t141) + t107 * t52 + t106 * t108 + t105 * t109 + t86 * (-mrSges(7,1) * t60 + mrSges(7,2) * t59) + t59 * t45 / 0.2e1 + t35 * t24 + t34 * t23 + t60 * t435 + t140 * t401 + t144 * t414 + (t144 * t47 - t145 * t46 + t17 * t211 - t18 * t210) * mrSges(6,3) + t112 * (mrSges(6,1) * t210 - mrSges(6,2) * t211) + (-Ifges(6,1) * t211 - Ifges(6,4) * t210) * t389 + (-Ifges(6,4) * t211 - Ifges(6,2) * t210) * t390 + (Ifges(6,5) * t145 + Ifges(6,6) * t144) * t376 + (Ifges(7,5) * t59 + Ifges(7,6) * t60) * t378 + (Ifges(6,1) * t145 + Ifges(6,4) * t144) * t385 + (Ifges(6,4) * t145 + Ifges(6,2) * t144) * t387 + (Ifges(7,1) * t59 + Ifges(7,4) * t60) * t392 + (Ifges(7,4) * t59 + Ifges(7,2) * t60) * t394 - t211 * t396 - t210 * t397 + (Ifges(7,4) * t141 + Ifges(7,2) * t140) * t398 + (Ifges(7,1) * t141 + Ifges(7,4) * t140) * t399 + t141 * t400; t419 * t52 - m(5) * (t126 * t142 + t128 * t143 + t183 * t197) + (-Ifges(6,5) * t201 - Ifges(6,6) * t200) * t377 + (-Ifges(6,1) * t201 - Ifges(6,4) * t200) * t386 + (-Ifges(6,4) * t201 - Ifges(6,2) * t200) * t388 + (-t218 / 0.2e1 + t201 / 0.2e1) * t83 + t420 * t138 + t421 * t139 + (t112 * t266 - t148 * t170 + t17 * t180 + t18 * t181 + t420 * t47 + t421 * t46) * m(6) + (-t134 / 0.2e1 + t89 / 0.2e1) * t44 + (-t135 / 0.2e1 + t88 / 0.2e1) * t45 + (mrSges(7,1) * t349 - mrSges(7,2) * t348) * t86 + (t12 * t348 - t13 * t349 + t159 * t2 - t160 * t3) * mrSges(7,3) + t112 * (mrSges(6,1) * t293 + mrSges(6,2) * t233) + (Ifges(6,1) * t233 - Ifges(6,4) * t293) * t389 + (Ifges(6,4) * t233 - Ifges(6,2) * t293) * t390 - t293 * t397 + (-t217 / 0.2e1 + t200 / 0.2e1) * t82 + (-t161 * t275 + t162 * t279) * pkin(9) + t405 * mrSges(5,3) + m(5) * (-pkin(3) * t156 + pkin(9) * t405) + (t306 * t275 * pkin(4) + (-m(5) * t294 - t194 * t275 - t195 * t279) * pkin(9) + t284) * qJD(4) + (-mrSges(4,1) - t302) * t156 + t333 * t197 + t275 * t116 / 0.2e1 + t266 * t68 - t196 * t249 + t204 * t14 - t143 * t194 - t142 * t195 + t180 * t109 + t181 * t108 - t170 * t102 + t65 * (-mrSges(7,1) * t159 + mrSges(7,2) * t160) - t155 * mrSges(4,2) - pkin(3) * t131 + t72 * t23 + t73 * t24 + t159 * t401 + (-mrSges(6,1) * t440 - mrSges(6,2) * t336) * t148 + (-t17 * t233 - t18 * t293 + t336 * t46 + t440 * t47) * mrSges(6,3) + (-Ifges(6,5) * t218 - Ifges(6,6) * t217) * t376 + (-Ifges(6,1) * t218 - Ifges(6,4) * t217) * t385 + (-Ifges(6,4) * t218 - Ifges(6,2) * t217) * t387 + t423 * t76 + t424 * t75 + (t12 * t423 + t13 * t424 + t2 * t73 + t204 * t65 + t3 * t72 + t419 * t86) * m(7) + ((t310 + t426 + t416) * t280 + ((-Ifges(6,3) / 0.2e1 - Ifges(5,3) / 0.2e1) * t261 + t309 + (t363 / 0.2e1 + (t427 + Ifges(4,2) / 0.2e1) * t280) * qJD(2) + t283) * t276 + (Ifges(6,5) * t233 + Ifges(7,5) * t160 - Ifges(6,6) * t293 + Ifges(7,6) * t159 + t295) * t327 / 0.2e1) * qJD(2) + t115 * t374 + (Ifges(7,5) * t88 + Ifges(7,6) * t89) * t378 + (Ifges(7,5) * t135 + Ifges(7,6) * t134) * t379 + t297 * t383 + t299 * t384 + (Ifges(7,1) * t88 + Ifges(7,4) * t89) * t392 + (Ifges(7,1) * t135 + Ifges(7,4) * t134) * t393 + (Ifges(7,4) * t88 + Ifges(7,2) * t89) * t394 + (Ifges(7,4) * t135 + Ifges(7,2) * t134) * t395 + t233 * t396 + (Ifges(7,4) * t160 + Ifges(7,2) * t159) * t398 + (Ifges(7,1) * t160 + Ifges(7,4) * t159) * t399 + t160 * t400; (-t102 * t240 + t108 * t270 + t109 * t272) * pkin(4) + t259 + t260 + t408 * t76 + (t12 * t408 + t13 * t409 - t133 * t86 + t2 * t213 + t212 * t3) * m(7) + t409 * t75 + (-Ifges(6,2) * t166 + t410 + t83) * t388 + t166 * t414 + (t166 * t47 + t304 * t46) * mrSges(6,3) + (Ifges(5,5) * t239 + Ifges(6,5) * t304 - Ifges(5,6) * t240 - Ifges(6,6) * t166) * t377 - t148 * (mrSges(6,1) * t166 + mrSges(6,2) * t304) - t402 + t94 * t435 + (-Ifges(5,2) * t240 + t154 + t231) * t382 + t439 + (t126 * t239 + t128 * t240) * mrSges(5,3) - t183 * (mrSges(5,1) * t240 + mrSges(5,2) * t239) + t212 * t23 + t213 * t24 - t126 * t194 + t128 * t195 - t53 * t139 - t54 * t138 - t133 * t52 + (-t148 * t372 - t46 * t53 - t47 * t54 + (t17 * t272 + t18 * t270) * pkin(4)) * m(6) + (Ifges(6,1) * t304 - t422) * t386 + t153 * t380 + (Ifges(5,1) * t239 - t356) * t381; -t304 * t138 + t166 * t139 - t415 * t75 + t94 * t76 + (t12 * t94 - t13 * t415 + t65) * m(7) + (t166 * t46 - t304 * t47 + t112) * m(6) + t411; -t12 * t75 + t13 * t76 + t44 * t392 + t417 + t439;];
tauc  = t1(:);
