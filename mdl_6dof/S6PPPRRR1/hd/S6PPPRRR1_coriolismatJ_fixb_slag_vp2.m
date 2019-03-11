% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPPRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:32
% EndTime: 2019-03-08 18:39:38
% DurationCPUTime: 4.37s
% Computational Cost: add. (16671->385), mult. (49958->620), div. (0->0), fcn. (61878->16), ass. (0->233)
t240 = sin(qJ(6));
t235 = t240 ^ 2;
t243 = cos(qJ(6));
t237 = t243 ^ 2;
t405 = t237 + t235;
t411 = mrSges(6,2) - t405 * (m(7) * pkin(11) + mrSges(7,3));
t242 = sin(qJ(4));
t344 = sin(pkin(14));
t345 = sin(pkin(7));
t298 = t345 * t344;
t375 = cos(qJ(4));
t239 = sin(pkin(8));
t347 = cos(pkin(14));
t300 = t347 * t345;
t348 = cos(pkin(8));
t349 = cos(pkin(7));
t403 = t349 * t239 + t348 * t300;
t161 = t242 * t403 + t375 * t298;
t241 = sin(qJ(5));
t244 = cos(qJ(5));
t256 = -t239 * t300 + t348 * t349;
t137 = t244 * t161 + t241 * t256;
t160 = t242 * t298 - t375 * t403;
t113 = -t137 * t240 + t160 * t243;
t114 = t137 * t243 + t160 * t240;
t136 = t161 * t241 - t244 * t256;
t356 = t243 * mrSges(7,2);
t359 = t240 * mrSges(7,1);
t296 = t356 + t359;
t206 = t296 * t241;
t332 = t240 * t244;
t215 = -t241 * mrSges(7,2) - mrSges(7,3) * t332;
t329 = t243 * t244;
t217 = t241 * mrSges(7,1) - mrSges(7,3) * t329;
t207 = t296 * t244;
t333 = t240 * t241;
t327 = mrSges(7,3) * t333;
t214 = mrSges(7,2) * t244 - t327;
t331 = t241 * t243;
t216 = -mrSges(7,1) * t244 - mrSges(7,3) * t331;
t220 = -pkin(5) * t244 - t241 * pkin(11) - pkin(4);
t182 = -pkin(10) * t332 + t243 * t220;
t183 = pkin(10) * t329 + t240 * t220;
t289 = t182 * t240 - t183 * t243;
t377 = -t243 / 0.2e1;
t379 = t240 / 0.2e1;
t393 = m(7) / 0.2e1;
t253 = (pkin(10) * t244 + t289) * t393 + t214 * t377 + t216 * t379 + t207 / 0.2e1;
t410 = -t113 * t217 / 0.2e1 - t114 * t215 / 0.2e1 - t137 * t206 / 0.2e1 - t253 * t136;
t378 = t241 / 0.2e1;
t357 = t241 * mrSges(6,2);
t409 = -t244 * mrSges(6,1) - mrSges(5,1) + t357;
t361 = Ifges(7,6) * t240;
t363 = Ifges(7,5) * t243;
t283 = t363 / 0.2e1 - t361 / 0.2e1;
t408 = Ifges(6,4) - t283;
t355 = t244 * mrSges(6,2);
t358 = t241 * mrSges(6,1);
t224 = t355 + t358;
t319 = t239 * t375;
t407 = -t224 * t319 / 0.2e1;
t346 = sin(pkin(6));
t301 = cos(pkin(13)) * t346;
t272 = t349 * t301;
t299 = t346 * sin(pkin(13));
t350 = cos(pkin(6));
t250 = -t272 * t347 + t299 * t344 - t300 * t350;
t255 = -t301 * t345 + t349 * t350;
t404 = -t255 * t239 + t250 * t348;
t368 = t241 * pkin(5);
t228 = -pkin(11) * t244 + t368;
t185 = pkin(10) * t333 + t228 * t243;
t186 = -pkin(10) * t331 + t228 * t240;
t402 = -t185 * t240 + t186 * t243;
t279 = t356 / 0.2e1 + t359 / 0.2e1;
t401 = t296 / 0.2e1 + t279;
t381 = t224 / 0.2e1;
t297 = mrSges(7,1) * t243 - mrSges(7,2) * t240;
t382 = -t297 / 0.2e1;
t400 = t241 * (t382 - mrSges(6,1) / 0.2e1) + t381 - t355 / 0.2e1;
t399 = -m(7) * pkin(5) - mrSges(6,1) - t297;
t398 = t297 * t378 + t355 / 0.2e1 + t358 / 0.2e1 + t381;
t396 = 2 * qJD(4);
t395 = m(6) / 0.2e1;
t394 = m(6) / 0.4e1;
t392 = mrSges(7,1) / 0.2e1;
t391 = -mrSges(7,2) / 0.2e1;
t159 = t272 * t344 + t298 * t350 + t299 * t347;
t112 = t159 * t375 - t242 * t404;
t246 = t239 * t250 + t255 * t348;
t94 = t112 * t241 - t244 * t246;
t390 = t94 / 0.2e1;
t205 = t297 * t241;
t389 = -t205 / 0.2e1;
t388 = t206 / 0.2e1;
t386 = -t214 / 0.2e1;
t385 = t215 / 0.2e1;
t384 = t216 / 0.2e1;
t383 = t217 / 0.2e1;
t380 = -t240 / 0.2e1;
t376 = t243 / 0.2e1;
t374 = pkin(5) * t205;
t373 = pkin(10) * t206;
t236 = t241 ^ 2;
t372 = pkin(10) * t236;
t238 = t244 ^ 2;
t371 = pkin(10) * t238;
t370 = pkin(11) * t214;
t369 = pkin(11) * t216;
t367 = t241 * pkin(10);
t366 = mrSges(7,3) * t241;
t365 = Ifges(7,4) * t240;
t364 = Ifges(7,4) * t243;
t362 = Ifges(7,5) * t244;
t360 = Ifges(7,6) * t244;
t111 = t159 * t242 + t375 * t404;
t95 = t112 * t244 + t241 * t246;
t67 = t111 * t243 - t240 * t95;
t354 = t67 * t240;
t68 = t111 * t240 + t243 * t95;
t353 = t68 * t243;
t77 = t111 * t332 + t112 * t243;
t352 = t77 * t240;
t78 = -t111 * t329 + t112 * t240;
t351 = t78 * t243;
t343 = t113 * t240;
t342 = t114 * t243;
t120 = t160 * t332 + t243 * t161;
t341 = t120 * t240;
t121 = -t160 * t329 + t240 * t161;
t340 = t121 * t243;
t335 = t239 * t242;
t203 = t241 * t348 + t244 * t335;
t162 = -t240 * t203 - t243 * t319;
t339 = t162 * t240;
t163 = t203 * t243 - t240 * t319;
t338 = t163 * t243;
t293 = -Ifges(7,2) * t240 + t364;
t274 = t293 * t241;
t190 = t274 - t360;
t334 = t240 * t190;
t295 = Ifges(7,1) * t243 - t365;
t275 = t295 * t241;
t192 = t275 - t362;
t330 = t243 * t192;
t326 = -m(7) * t241 / 0.4e1;
t325 = -t366 / 0.2e1;
t324 = mrSges(7,3) * t380;
t322 = mrSges(7,3) * t376;
t318 = t241 * t375;
t317 = t244 * t375;
t308 = t327 / 0.2e1;
t307 = t243 * t325;
t306 = t236 * t319;
t305 = t239 * t318;
t202 = t241 * t335 - t244 * t348;
t304 = t202 * t318;
t294 = Ifges(7,1) * t240 + t364;
t292 = Ifges(7,2) * t243 + t365;
t291 = -t361 + t363;
t290 = Ifges(7,5) * t240 + Ifges(7,6) * t243;
t288 = -t202 * t241 - t203 * t244;
t130 = (mrSges(7,2) * pkin(5) - t364) * t243 + (pkin(5) * mrSges(7,1) + t365 + (-Ifges(7,1) + Ifges(7,2)) * t243) * t240;
t280 = t185 * t392 + t186 * t391;
t52 = t374 / 0.2e1 + (Ifges(7,3) / 0.2e1 + (t237 / 0.2e1 + t235 / 0.2e1) * pkin(11) * mrSges(7,3)) * t241 + (0.3e1 / 0.4e1 * t362 + t369 / 0.2e1 - t192 / 0.4e1 + (pkin(10) * t391 + (Ifges(7,2) / 0.2e1 - Ifges(7,1) / 0.4e1) * t243) * t241) * t243 + (-0.3e1 / 0.4e1 * t360 + t370 / 0.2e1 + t190 / 0.4e1 + (0.3e1 / 0.2e1 * t364 - pkin(10) * mrSges(7,1) / 0.2e1 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.4e1) * t240) * t241) * t240 + t280;
t287 = t52 * qJD(4) + t130 * qJD(5);
t286 = t391 * t78 + t392 * t77;
t285 = -t353 + t95 + t354;
t284 = -t241 * t94 - t244 * t95 + t112;
t282 = t120 * t392 + t121 * t391;
t180 = (-t240 * t317 + t242 * t243) * t239;
t181 = (t240 * t242 + t243 * t317) * t239;
t281 = t180 * t392 + t181 * t391;
t278 = t137 - t342 + t343;
t277 = -t136 * t241 - t137 * t244 + t161;
t276 = t203 - t338 + t339;
t10 = m(7) * (t67 * t77 + t68 * t78) + 0.4e1 * (t284 * t394 + t326 * t94) * t111;
t5 = (t113 * t77 + t114 * t78 + t120 * t67 + t121 * t68 + (-t111 * t136 - t160 * t94) * t241) * t393 + (t111 * t277 + t160 * t284) * t395;
t9 = (t162 * t77 + t163 * t78 + t180 * t67 + t181 * t68 + (-t111 * t202 + t319 * t94) * t241) * t393 + (t288 * t111 + (t111 * t242 - t112 * t375 + t317 * t95 + t318 * t94) * t239) * t395;
t270 = t10 * qJD(1) + t5 * qJD(2) + t9 * qJD(3);
t12 = (t202 * t285 + t276 * t94) * t393;
t13 = m(7) * t285 * t94;
t8 = (t136 * t285 + t278 * t94) * t393;
t269 = t13 * qJD(1) + t8 * qJD(2) + t12 * qJD(3);
t30 = (t180 * t113 + t181 * t114 + t162 * t120 + t163 * t121 + (t136 * t319 - t160 * t202) * t241) * t393 + (t288 * t160 + (t136 * t318 + t137 * t317 + t160 * t242 - t161 * t375) * t239) * t395;
t33 = m(7) * (t113 * t120 + t114 * t121) + 0.4e1 * (t136 * t326 + t277 * t394) * t160;
t268 = t5 * qJD(1) + t33 * qJD(2) + t30 * qJD(3);
t83 = m(7) * (t162 * t180 + t163 * t181 + t239 * t304) + m(6) * (t203 * t317 - t242 * t319 + t304) * t239;
t267 = t9 * qJD(1) + t30 * qJD(2) + t83 * qJD(3);
t32 = (t136 * t276 + t202 * t278) * t393;
t38 = m(7) * t278 * t136;
t266 = t8 * qJD(1) + t38 * qJD(2) + t32 * qJD(3);
t80 = m(7) * t276 * t202;
t265 = t12 * qJD(1) + t32 * qJD(2) + t80 * qJD(3);
t264 = t185 * t113 + t186 * t114 + t137 * t367;
t263 = qJD(4) * (-t241 * t206 + mrSges(5,2) + (-t236 - t238) * mrSges(6,3));
t262 = t113 * t386 + t114 * t384 + t136 * t389;
t260 = t162 * t386 + t163 * t384 + t202 * t389;
t251 = (t340 / 0.2e1 - t341 / 0.2e1) * mrSges(7,3) + (t160 * t368 + (t340 - t341) * pkin(11)) * t393;
t14 = -t400 * t160 - m(7) * t264 / 0.2e1 + t251 + t410;
t245 = (t185 * t162 + t186 * t163 + t203 * t367) * t393 + t162 * t383 + t163 * t385 + t203 * t388 + t407 + t253 * t202;
t247 = (-pkin(5) * t305 + (-t180 * t240 + t181 * t243) * pkin(11)) * t393 + t180 * t324 + t181 * t322 + t305 * t382 + t407;
t36 = -t245 + t247;
t248 = t253 * t94 + (t185 * t67 + t186 * t68 + t367 * t95) * t393 + t67 * t383 + t68 * t385 + t95 * t388;
t257 = m(7) * (t111 * t368 + (t351 - t352) * pkin(11));
t4 = (-t351 / 0.2e1 + t352 / 0.2e1) * mrSges(7,3) + t400 * t111 - t257 / 0.2e1 + t248;
t191 = Ifges(7,6) * t241 + t244 * t293;
t193 = Ifges(7,5) * t241 + t244 * t295;
t48 = t186 * t214 + t183 * t215 + t185 * t216 + t182 * t217 + m(7) * (t182 * t185 + t183 * t186) - pkin(4) * t224 + (t373 + t330 / 0.2e1 - t334 / 0.2e1 + t408 * t244) * t244 + (pkin(10) * t207 + t193 * t376 + t191 * t380 - t408 * t241 + (m(7) * pkin(10) ^ 2 + Ifges(6,1) - Ifges(6,2) - Ifges(7,3)) * t244) * t241;
t259 = t4 * qJD(1) - t14 * qJD(2) - t36 * qJD(3) + t48 * qJD(4);
t26 = (t342 / 0.2e1 - t343 / 0.2e1) * t366 + t262 + t282;
t50 = (t338 / 0.2e1 - t339 / 0.2e1) * t366 + t260 + t281;
t55 = t182 * t214 - t183 * t216 + (pkin(10) * t205 + t192 * t380 + t190 * t377 + t244 * t290 / 0.2e1 + (t292 * t379 + t294 * t377) * t241 + t289 * mrSges(7,3)) * t241;
t252 = (-t353 / 0.2e1 + t354 / 0.2e1) * t366 + t67 * t214 / 0.2e1 - t68 * t216 / 0.2e1 + t205 * t390;
t6 = t252 - t286;
t258 = t6 * qJD(1) - t26 * qJD(2) - t50 * qJD(3) + t55 * qJD(4);
t221 = pkin(10) * t306;
t151 = t160 * t372;
t108 = t111 * t372;
t76 = t401 * t202;
t53 = -t374 / 0.2e1 + t370 * t380 + t369 * t377 + t330 / 0.4e1 - t294 * t333 / 0.2e1 - t334 / 0.4e1 - t292 * t331 / 0.2e1 + t243 * t275 / 0.4e1 - t240 * t274 / 0.4e1 + t373 / 0.2e1 + Ifges(7,3) * t378 + t280 + t405 * pkin(11) * t325 + (-t291 / 0.4e1 + t283) * t244;
t51 = t162 * t308 + t163 * t307 - t260 + t281;
t40 = t401 * t136;
t37 = t245 + t247;
t27 = t113 * t308 + t114 * t307 - t262 + t282;
t19 = t279 * t94 + t296 * t390;
t15 = t160 * t398 + t264 * t393 + t251 - t410;
t11 = t30 * qJD(4) + t32 * qJD(5);
t7 = t252 + t286;
t3 = t257 / 0.2e1 + t78 * t322 + t77 * t324 + t248 + t398 * t111;
t2 = t9 * qJD(4) + t12 * qJD(5);
t1 = t5 * qJD(4) + t8 * qJD(5);
t16 = [t10 * qJD(4) + t13 * qJD(5), t1, t2 (t112 * t409 + t78 * t214 + t77 * t216) * qJD(4) + t3 * qJD(5) + t7 * qJD(6) + t111 * t263 + ((t182 * t77 + t183 * t78 - t108) * t393 + (-pkin(4) * t112 - t111 * t371 - t108) * t395) * t396 + t270, t3 * qJD(4) + (t399 * t95 + t411 * t94) * qJD(5) + t19 * qJD(6) + t269, t7 * qJD(4) + t19 * qJD(5) + (-mrSges(7,1) * t68 - mrSges(7,2) * t67) * qJD(6); t1, t33 * qJD(4) + t38 * qJD(5), t11 (t120 * t216 + t121 * t214 + t161 * t409) * qJD(4) + t15 * qJD(5) + t27 * qJD(6) + t160 * t263 + ((t120 * t182 + t121 * t183 - t151) * t393 + (-pkin(4) * t161 - t160 * t371 - t151) * t395) * t396 + t268, t15 * qJD(4) + (t411 * t136 + t399 * t137) * qJD(5) + t40 * qJD(6) + t266, t27 * qJD(4) + t40 * qJD(5) + (-mrSges(7,1) * t114 - mrSges(7,2) * t113) * qJD(6); t2, t11, t83 * qJD(4) + t80 * qJD(5) (t206 * t305 + m(7) * (t180 * t182 + t181 * t183 + t221) + t181 * t214 + t180 * t216 + m(6) * (t221 + (-pkin(4) * t242 + t371 * t375) * t239) - mrSges(5,2) * t319 + t409 * t335 + (t238 * t319 + t306) * mrSges(6,3)) * qJD(4) + t37 * qJD(5) + t51 * qJD(6) + t267, t37 * qJD(4) + (t411 * t202 + t399 * t203) * qJD(5) + t76 * qJD(6) + t265, t51 * qJD(4) + t76 * qJD(5) + (-mrSges(7,1) * t163 - mrSges(7,2) * t162) * qJD(6); qJD(5) * t4 + qJD(6) * t6 - t270, -qJD(5) * t14 - qJD(6) * t26 - t268, -qJD(5) * t36 - qJD(6) * t50 - t267, qJD(5) * t48 + qJD(6) * t55, t53 * qJD(6) + t259 + (-Ifges(6,6) * t241 - pkin(5) * t207 + pkin(10) * t357 + t191 * t376 + t193 * t379 + t290 * t378 + (m(7) * t402 + t243 * t215 - t240 * t217) * pkin(11) + (pkin(10) * t399 + t292 * t380 + t294 * t376 + Ifges(6,5)) * t244 + t402 * mrSges(7,3)) * qJD(5), t53 * qJD(5) + (-mrSges(7,1) * t183 - mrSges(7,2) * t182 - t241 * t290) * qJD(6) + t258; -qJD(4) * t4 - t269, qJD(4) * t14 - t266, qJD(4) * t36 - t265, -qJD(6) * t52 - t259, -t130 * qJD(6) (-pkin(11) * t297 + t291) * qJD(6) - t287; -t6 * qJD(4), t26 * qJD(4), t50 * qJD(4), qJD(5) * t52 - t258, t287, 0;];
Cq  = t16;
