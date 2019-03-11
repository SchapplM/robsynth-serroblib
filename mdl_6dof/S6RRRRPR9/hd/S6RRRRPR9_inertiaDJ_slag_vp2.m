% Calculate time derivative of joint inertia matrix for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:53
% EndTime: 2019-03-09 22:48:09
% DurationCPUTime: 7.64s
% Computational Cost: add. (17955->735), mult. (46270->1075), div. (0->0), fcn. (46197->12), ass. (0->297)
t314 = sin(pkin(12));
t316 = cos(pkin(12));
t288 = -mrSges(6,1) * t316 + mrSges(6,2) * t314;
t415 = t288 - mrSges(5,1);
t387 = t316 / 0.2e1;
t388 = t314 / 0.2e1;
t318 = sin(qJ(6));
t322 = cos(qJ(6));
t278 = t314 * t322 + t316 * t318;
t392 = t278 / 0.2e1;
t328 = t314 * t318 - t316 * t322;
t393 = -t328 / 0.2e1;
t414 = Ifges(6,5) * t388 + Ifges(7,5) * t392 + Ifges(6,6) * t387 + Ifges(7,6) * t393;
t317 = cos(pkin(6));
t321 = sin(qJ(2));
t315 = sin(pkin(6));
t325 = cos(qJ(2));
t357 = t315 * t325;
t274 = t317 * t321 * pkin(1) + pkin(8) * t357;
t259 = pkin(9) * t317 + t274;
t260 = (-pkin(2) * t325 - pkin(9) * t321 - pkin(1)) * t315;
t320 = sin(qJ(3));
t324 = cos(qJ(3));
t197 = t324 * t259 + t320 * t260;
t403 = -pkin(10) - pkin(9);
t295 = t403 * t320;
t296 = t403 * t324;
t319 = sin(qJ(4));
t323 = cos(qJ(4));
t413 = t323 * t295 + t296 * t319;
t358 = t315 * t321;
t301 = pkin(8) * t358;
t384 = pkin(1) * t325;
t273 = t317 * t384 - t301;
t268 = t328 * qJD(6);
t412 = qJD(3) + qJD(4);
t271 = t317 * t320 + t324 * t358;
t354 = qJD(2) * t315;
t339 = t325 * t354;
t243 = -qJD(3) * t271 - t320 * t339;
t270 = t317 * t324 - t320 * t358;
t244 = qJD(3) * t270 + t324 * t339;
t353 = qJD(2) * t321;
t340 = t315 * t353;
t411 = Ifges(4,5) * t244 + Ifges(4,6) * t243 + Ifges(4,3) * t340;
t231 = mrSges(7,1) * t328 + mrSges(7,2) * t278;
t410 = ((t231 + t415) * t319 - mrSges(5,2) * t323) * pkin(3) * qJD(4);
t409 = 2 * m(5);
t408 = 2 * m(6);
t407 = 2 * m(7);
t406 = -2 * mrSges(3,3);
t269 = t278 * qJD(6);
t217 = mrSges(7,1) * t269 - t268 * mrSges(7,2);
t405 = 0.2e1 * t217;
t264 = t274 * qJD(2);
t404 = 0.2e1 * t264;
t329 = t323 * t270 - t271 * t319;
t128 = qJD(4) * t329 + t243 * t319 + t244 * t323;
t107 = -t128 * t314 + t316 * t340;
t402 = t107 / 0.2e1;
t218 = -Ifges(7,5) * t268 - Ifges(7,6) * t269;
t401 = t218 / 0.2e1;
t219 = -Ifges(7,4) * t268 - Ifges(7,2) * t269;
t400 = t219 / 0.2e1;
t220 = -Ifges(7,1) * t268 - Ifges(7,4) * t269;
t399 = t220 / 0.2e1;
t233 = Ifges(7,4) * t278 - Ifges(7,2) * t328;
t397 = t233 / 0.2e1;
t234 = Ifges(7,1) * t278 - Ifges(7,4) * t328;
t396 = t234 / 0.2e1;
t395 = -t268 / 0.2e1;
t394 = -t269 / 0.2e1;
t375 = Ifges(6,4) * t316;
t390 = Ifges(6,1) * t388 + t375 / 0.2e1;
t389 = -t314 / 0.2e1;
t386 = m(7) * t319;
t262 = (pkin(2) * t321 - pkin(9) * t325) * t354;
t263 = t273 * qJD(2);
t351 = qJD(3) * t324;
t352 = qJD(3) * t320;
t140 = -t259 * t352 + t260 * t351 + t320 * t262 + t324 * t263;
t106 = pkin(10) * t243 + t140;
t196 = -t320 * t259 + t324 * t260;
t157 = -pkin(3) * t357 - t271 * pkin(10) + t196;
t174 = pkin(10) * t270 + t197;
t349 = qJD(4) * t323;
t350 = qJD(4) * t319;
t141 = -t197 * qJD(3) + t324 * t262 - t263 * t320;
t99 = pkin(3) * t340 - pkin(10) * t244 + t141;
t30 = t323 * t106 + t157 * t349 - t174 * t350 + t319 * t99;
t27 = (qJ(5) * t353 - qJD(5) * t325) * t315 + t30;
t209 = t270 * t319 + t271 * t323;
t129 = qJD(4) * t209 - t323 * t243 + t244 * t319;
t193 = -t243 * pkin(3) + t264;
t48 = t129 * pkin(4) - t128 * qJ(5) - t209 * qJD(5) + t193;
t12 = -t27 * t314 + t316 * t48;
t385 = mrSges(6,3) * t12;
t383 = pkin(3) * t323;
t279 = t319 * t320 - t323 * t324;
t237 = t412 * t279;
t280 = t319 * t324 + t320 * t323;
t238 = t412 * t280;
t345 = pkin(3) * t352;
t147 = pkin(4) * t238 + qJ(5) * t237 - qJD(5) * t280 + t345;
t341 = qJD(3) * t403;
t285 = t320 * t341;
t335 = t324 * t341;
t179 = t413 * qJD(4) + t323 * t285 + t319 * t335;
t80 = t316 * t147 - t179 * t314;
t382 = t80 * mrSges(6,3);
t13 = t316 * t27 + t314 * t48;
t380 = mrSges(7,3) * t269;
t379 = mrSges(7,3) * t328;
t378 = Ifges(4,4) * t320;
t377 = Ifges(4,4) * t324;
t376 = Ifges(6,4) * t314;
t373 = t13 * t316;
t372 = t263 * mrSges(3,2);
t371 = t264 * mrSges(3,1);
t81 = t314 * t147 + t316 * t179;
t370 = t316 * t81;
t258 = t301 + (-pkin(2) - t384) * t317;
t216 = -t270 * pkin(3) + t258;
t109 = -pkin(4) * t329 - t209 * qJ(5) + t216;
t94 = t319 * t157 + t323 * t174;
t87 = -qJ(5) * t357 + t94;
t54 = t314 * t109 + t316 * t87;
t250 = t295 * t319 - t296 * t323;
t180 = qJD(4) * t250 + t285 * t319 - t323 * t335;
t369 = t180 * t413;
t368 = t237 * t314;
t367 = t237 * t316;
t366 = t413 * t319;
t365 = t280 * t314;
t364 = t280 * t316;
t300 = pkin(3) * t349 + qJD(5);
t362 = t300 * t314;
t361 = t300 * t316;
t305 = pkin(3) * t319 + qJ(5);
t360 = t305 * t314;
t359 = t305 * t316;
t183 = -t314 * t209 - t316 * t357;
t184 = t316 * t209 - t314 * t357;
t122 = -mrSges(6,1) * t183 + mrSges(6,2) * t184;
t192 = -mrSges(5,1) * t357 - t209 * mrSges(5,3);
t356 = -t192 + t122;
t160 = -mrSges(6,1) * t368 - mrSges(6,2) * t367;
t309 = -pkin(3) * t324 - pkin(2);
t224 = pkin(4) * t279 - qJ(5) * t280 + t309;
t166 = t314 * t224 + t316 * t250;
t355 = t314 ^ 2 + t316 ^ 2;
t348 = 0.2e1 * mrSges(7,3);
t347 = 0.2e1 * t315;
t108 = t128 * t316 + t314 * t340;
t120 = t183 * t322 - t184 * t318;
t38 = qJD(6) * t120 + t107 * t318 + t108 * t322;
t121 = t183 * t318 + t184 * t322;
t39 = -qJD(6) * t121 + t107 * t322 - t108 * t318;
t8 = Ifges(7,5) * t38 + Ifges(7,6) * t39 + Ifges(7,3) * t129;
t344 = pkin(3) * t350;
t343 = Ifges(4,6) * t357;
t116 = t237 * t328 - t269 * t280;
t117 = t237 * t278 + t280 * t268;
t58 = Ifges(7,5) * t116 + Ifges(7,6) * t117 + Ifges(7,3) * t238;
t342 = Ifges(5,5) * t128 - Ifges(5,6) * t129 + Ifges(5,3) * t340;
t306 = -pkin(5) * t316 - pkin(4);
t14 = -t39 * mrSges(7,1) + t38 * mrSges(7,2);
t61 = -t107 * mrSges(6,1) + t108 * mrSges(6,2);
t65 = -t117 * mrSges(7,1) + t116 * mrSges(7,2);
t53 = t316 * t109 - t314 * t87;
t338 = t355 * t300;
t93 = t157 * t323 - t319 * t174;
t165 = t316 * t224 - t250 * t314;
t337 = t355 * qJD(5);
t336 = -t219 * t328 + t278 * t220 - t269 * t233 - t268 * t234;
t88 = pkin(4) * t357 - t93;
t333 = Ifges(6,1) * t316 - t376;
t332 = -Ifges(6,2) * t314 + t375;
t331 = Ifges(6,5) * t316 - Ifges(6,6) * t314;
t34 = -pkin(5) * t329 - pkin(11) * t184 + t53;
t45 = pkin(11) * t183 + t54;
t16 = -t318 * t45 + t322 * t34;
t17 = t318 * t34 + t322 * t45;
t330 = 0.2e1 * t355 * mrSges(6,3);
t133 = pkin(5) * t279 - pkin(11) * t364 + t165;
t146 = -pkin(11) * t365 + t166;
t73 = t133 * t322 - t146 * t318;
t74 = t133 * t318 + t146 * t322;
t275 = (-pkin(11) - t305) * t314;
t311 = t316 * pkin(11);
t276 = t311 + t359;
t214 = t275 * t322 - t276 * t318;
t215 = t275 * t318 + t276 * t322;
t287 = (-pkin(11) - qJ(5)) * t314;
t289 = qJ(5) * t316 + t311;
t245 = t287 * t322 - t289 * t318;
t246 = t287 * t318 + t289 * t322;
t31 = -t319 * t106 - t157 * t350 - t174 * t349 + t323 * t99;
t28 = -pkin(4) * t340 - t31;
t10 = Ifges(7,1) * t38 + Ifges(7,4) * t39 + Ifges(7,5) * t129;
t4 = pkin(5) * t129 - pkin(11) * t108 + t12;
t7 = pkin(11) * t107 + t13;
t2 = qJD(6) * t16 + t318 * t4 + t322 * t7;
t22 = -pkin(5) * t107 + t28;
t291 = Ifges(6,2) * t316 + t376;
t3 = -qJD(6) * t17 - t318 * t7 + t322 * t4;
t43 = Ifges(6,4) * t108 + Ifges(6,2) * t107 + Ifges(6,6) * t129;
t44 = Ifges(6,1) * t108 + Ifges(6,4) * t107 + Ifges(6,5) * t129;
t56 = Ifges(7,4) * t121 + Ifges(7,2) * t120 - Ifges(7,6) * t329;
t57 = Ifges(7,1) * t121 + Ifges(7,4) * t120 - Ifges(7,5) * t329;
t72 = -pkin(5) * t183 + t88;
t9 = Ifges(7,4) * t38 + Ifges(7,2) * t39 + Ifges(7,6) * t129;
t327 = t43 * t387 + t44 * t388 + t108 * t390 + t10 * t392 + t9 * t393 + t56 * t394 + t57 * t395 + t38 * t396 + t39 * t397 + t121 * t399 + t120 * t400 - t329 * t401 + t291 * t402 - t2 * t379 - t17 * t380 + mrSges(6,3) * t373 + t342 + t28 * t288 + t22 * t231 + t72 * t217 + (t16 * t268 - t278 * t3) * mrSges(7,3) + t31 * mrSges(5,1) - t30 * mrSges(5,2) + t414 * t129;
t131 = Ifges(6,6) * t238 - t237 * t332;
t132 = Ifges(6,5) * t238 - t237 * t333;
t210 = t278 * t280;
t211 = t328 * t280;
t137 = -Ifges(7,4) * t211 - Ifges(7,2) * t210 + Ifges(7,6) * t279;
t138 = -Ifges(7,1) * t211 - Ifges(7,4) * t210 + Ifges(7,5) * t279;
t139 = -pkin(5) * t368 + t180;
t69 = pkin(5) * t238 + pkin(11) * t367 + t80;
t75 = pkin(11) * t368 + t81;
t19 = qJD(6) * t73 + t318 * t69 + t322 * t75;
t20 = -qJD(6) * t74 - t318 * t75 + t322 * t69;
t204 = pkin(5) * t365 - t413;
t229 = Ifges(5,6) * t238;
t230 = Ifges(5,5) * t237;
t59 = Ifges(7,4) * t116 + Ifges(7,2) * t117 + Ifges(7,6) * t238;
t60 = Ifges(7,1) * t116 + Ifges(7,4) * t117 + Ifges(7,5) * t238;
t326 = t131 * t387 + t132 * t388 - t367 * t390 + t60 * t392 + t59 * t393 + t137 * t394 + t138 * t395 + t116 * t396 + t117 * t397 - t211 * t399 - t210 * t400 + t279 * t401 - t19 * t379 - t74 * t380 + mrSges(6,3) * t370 - t230 - t229 + t139 * t231 + t204 * t217 - t179 * mrSges(5,2) + (-t20 * t278 + t268 * t73) * mrSges(7,3) + t291 * t368 / 0.2e1 + t414 * t238 + t415 * t180;
t310 = Ifges(4,5) * t351;
t308 = -pkin(4) - t383;
t299 = Ifges(3,5) * t339;
t294 = Ifges(4,1) * t320 + t377;
t293 = Ifges(4,2) * t324 + t378;
t286 = t306 - t383;
t284 = (Ifges(4,1) * t324 - t378) * qJD(3);
t283 = (-Ifges(4,2) * t320 + t377) * qJD(3);
t282 = (mrSges(4,1) * t320 + mrSges(4,2) * t324) * qJD(3);
t248 = -mrSges(4,1) * t357 - t271 * mrSges(4,3);
t247 = mrSges(4,2) * t357 + t270 * mrSges(4,3);
t241 = Ifges(5,1) * t280 - Ifges(5,4) * t279;
t240 = Ifges(5,4) * t280 - Ifges(5,2) * t279;
t239 = mrSges(5,1) * t279 + mrSges(5,2) * t280;
t226 = mrSges(6,1) * t279 - mrSges(6,3) * t364;
t225 = -mrSges(6,2) * t279 - mrSges(6,3) * t365;
t221 = (mrSges(6,1) * t314 + mrSges(6,2) * t316) * t280;
t207 = mrSges(4,1) * t340 - mrSges(4,3) * t244;
t206 = -mrSges(4,2) * t340 + mrSges(4,3) * t243;
t203 = Ifges(4,1) * t271 + Ifges(4,4) * t270 - Ifges(4,5) * t357;
t202 = Ifges(4,4) * t271 + Ifges(4,2) * t270 - t343;
t199 = -qJD(5) * t278 - qJD(6) * t246;
t198 = -qJD(5) * t328 + qJD(6) * t245;
t191 = mrSges(5,2) * t357 + mrSges(5,3) * t329;
t190 = Ifges(6,5) * t279 + t280 * t333;
t189 = Ifges(6,6) * t279 + t280 * t332;
t188 = Ifges(6,3) * t279 + t280 * t331;
t182 = mrSges(7,1) * t279 + mrSges(7,3) * t211;
t181 = -mrSges(7,2) * t279 - mrSges(7,3) * t210;
t177 = -mrSges(4,1) * t243 + mrSges(4,2) * t244;
t173 = -Ifges(5,1) * t237 - Ifges(5,4) * t238;
t172 = -Ifges(5,4) * t237 - Ifges(5,2) * t238;
t171 = mrSges(5,1) * t238 - mrSges(5,2) * t237;
t168 = mrSges(6,1) * t238 + mrSges(6,3) * t367;
t167 = -mrSges(6,2) * t238 + mrSges(6,3) * t368;
t164 = -qJD(6) * t215 - t278 * t300;
t163 = qJD(6) * t214 - t300 * t328;
t159 = Ifges(4,1) * t244 + Ifges(4,4) * t243 + Ifges(4,5) * t340;
t158 = Ifges(4,4) * t244 + Ifges(4,2) * t243 + Ifges(4,6) * t340;
t149 = mrSges(7,1) * t210 - mrSges(7,2) * t211;
t148 = -mrSges(5,1) * t329 + mrSges(5,2) * t209;
t143 = Ifges(5,1) * t209 + Ifges(5,4) * t329 - Ifges(5,5) * t357;
t142 = Ifges(5,4) * t209 + Ifges(5,2) * t329 - Ifges(5,6) * t357;
t136 = -Ifges(7,5) * t211 - Ifges(7,6) * t210 + Ifges(7,3) * t279;
t135 = -mrSges(6,1) * t329 - mrSges(6,3) * t184;
t134 = mrSges(6,2) * t329 + mrSges(6,3) * t183;
t130 = Ifges(6,3) * t238 - t237 * t331;
t112 = -mrSges(5,2) * t340 - mrSges(5,3) * t129;
t111 = mrSges(5,1) * t340 - mrSges(5,3) * t128;
t90 = -mrSges(7,2) * t238 + mrSges(7,3) * t117;
t89 = mrSges(7,1) * t238 - mrSges(7,3) * t116;
t86 = Ifges(6,1) * t184 + Ifges(6,4) * t183 - Ifges(6,5) * t329;
t85 = Ifges(6,4) * t184 + Ifges(6,2) * t183 - Ifges(6,6) * t329;
t84 = Ifges(6,5) * t184 + Ifges(6,6) * t183 - Ifges(6,3) * t329;
t77 = -mrSges(7,1) * t329 - mrSges(7,3) * t121;
t76 = mrSges(7,2) * t329 + mrSges(7,3) * t120;
t70 = mrSges(5,1) * t129 + mrSges(5,2) * t128;
t68 = Ifges(5,1) * t128 - Ifges(5,4) * t129 + Ifges(5,5) * t340;
t67 = Ifges(5,4) * t128 - Ifges(5,2) * t129 + Ifges(5,6) * t340;
t66 = -mrSges(7,1) * t120 + mrSges(7,2) * t121;
t64 = mrSges(6,1) * t129 - mrSges(6,3) * t108;
t63 = -mrSges(6,2) * t129 + mrSges(6,3) * t107;
t55 = Ifges(7,5) * t121 + Ifges(7,6) * t120 - Ifges(7,3) * t329;
t42 = Ifges(6,5) * t108 + Ifges(6,6) * t107 + Ifges(6,3) * t129;
t24 = -mrSges(7,2) * t129 + mrSges(7,3) * t39;
t23 = mrSges(7,1) * t129 - mrSges(7,3) * t38;
t1 = [(mrSges(3,3) * t321 * t404 + (0.2e1 * t263 * mrSges(3,3) - t342 - t411) * t325 + ((t273 * t406 + Ifges(3,5) * t317 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t325) * t347) * t325 + (t274 * t406 + Ifges(4,5) * t271 + Ifges(5,5) * t209 - 0.2e1 * Ifges(3,6) * t317 + Ifges(4,6) * t270 + Ifges(5,6) * t329 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t321) * t347 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3)) * t357) * t321) * qJD(2)) * t315 - (t42 + t8 - t67) * t329 + (-t142 + t84 + t55) * t129 + (-mrSges(4,1) * t270 + mrSges(4,2) * t271) * t404 + (t16 * t3 + t17 * t2 + t22 * t72) * t407 + (t12 * t53 + t13 * t54 + t28 * t88) * t408 + (t193 * t216 + t30 * t94 + t31 * t93) * t409 + t271 * t159 + t270 * t158 + 0.2e1 * t258 * t177 + t243 * t202 + t244 * t203 + 0.2e1 * t140 * t247 + 0.2e1 * t141 * t248 + t209 * t68 + 0.2e1 * t216 * t70 + 0.2e1 * t197 * t206 + 0.2e1 * t196 * t207 + t183 * t43 + t184 * t44 + 0.2e1 * t30 * t191 + 0.2e1 * t31 * t192 + 0.2e1 * t193 * t148 + t128 * t143 + 0.2e1 * t13 * t134 + 0.2e1 * t12 * t135 + t120 * t9 + t121 * t10 + 0.2e1 * t28 * t122 + t107 * t85 + t108 * t86 + 0.2e1 * t93 * t111 + 0.2e1 * t94 * t112 + 0.2e1 * t88 * t61 + 0.2e1 * t2 * t76 + 0.2e1 * t3 * t77 + 0.2e1 * t72 * t14 + 0.2e1 * t54 * t63 + 0.2e1 * t53 * t64 + 0.2e1 * t22 * t66 + t39 * t56 + t38 * t57 + 0.2e1 * t16 * t23 + 0.2e1 * t17 * t24 + 0.2e1 * m(3) * (t263 * t274 - t264 * t273) + 0.2e1 * m(4) * (t140 * t197 + t141 * t196 + t258 * t264) + (t299 - 0.2e1 * t371 - 0.2e1 * t372) * t317; -(t61 - t111) * t413 + m(6) * (t12 * t165 + t13 * t166 + t180 * t88 - t28 * t413 + t53 * t80 + t54 * t81) + m(5) * (t179 * t94 - t180 * t93 + t193 * t309 + t216 * t345 + t250 * t30 + t31 * t413) + ((-t310 / 0.2e1 + t230 / 0.2e1 + t229 / 0.2e1) * t325 + (Ifges(5,5) * t280 / 0.2e1 - Ifges(5,6) * t279 / 0.2e1 - Ifges(3,6) + Ifges(4,5) * t320 / 0.2e1 + Ifges(4,6) * t324 / 0.2e1) * t353) * t315 - (t130 / 0.2e1 + t58 / 0.2e1 - t172 / 0.2e1) * t329 + ((t140 * t324 - t141 * t320 + (-t196 * t324 - t197 * t320) * qJD(3)) * pkin(9) - pkin(2) * t264) * m(4) + m(7) * (t139 * t72 + t16 * t20 + t17 * t19 + t2 * t74 + t204 * t22 + t3 * t73) + t189 * t402 - t371 - t372 + t299 + t309 * t70 + t258 * t282 + t270 * t283 / 0.2e1 + t271 * t284 / 0.2e1 + t243 * t293 / 0.2e1 + t244 * t294 / 0.2e1 + t128 * t241 / 0.2e1 + t250 * t112 + t193 * t239 + t28 * t221 + t13 * t225 + t12 * t226 + t209 * t173 / 0.2e1 - t210 * t9 / 0.2e1 - t211 * t10 / 0.2e1 + t216 * t171 + t204 * t14 + t183 * t131 / 0.2e1 + t184 * t132 / 0.2e1 + t108 * t190 / 0.2e1 + t179 * t191 + t2 * t181 + t3 * t182 - pkin(2) * t177 + t88 * t160 + t165 * t64 + t166 * t63 + t54 * t167 + t53 * t168 + t22 * t149 + t81 * t134 + t80 * t135 + t39 * t137 / 0.2e1 + t38 * t138 / 0.2e1 + t139 * t66 + t120 * t59 / 0.2e1 + t121 * t60 / 0.2e1 + t116 * t57 / 0.2e1 + t117 * t56 / 0.2e1 + t16 * t89 + t17 * t90 + t73 * t23 + t74 * t24 + t19 * t76 + t20 * t77 + t72 * t65 + (-t240 / 0.2e1 + t188 / 0.2e1 + t136 / 0.2e1) * t129 + (t237 * t93 - t238 * t94 - t279 * t30 - t280 * t31) * mrSges(5,3) + ((t203 / 0.2e1 - pkin(9) * t248 - t196 * mrSges(4,3)) * t324 + (t343 / 0.2e1 - t202 / 0.2e1 - pkin(9) * t247 - t197 * mrSges(4,3) + pkin(3) * t148) * t320) * qJD(3) + (t264 * mrSges(4,2) + t159 / 0.2e1 - pkin(9) * t207 - t141 * mrSges(4,3)) * t320 + t356 * t180 + (-t264 * mrSges(4,1) + t158 / 0.2e1 + pkin(9) * t206 + t140 * mrSges(4,3)) * t324 - (t143 / 0.2e1 + t86 * t387 + t85 * t389) * t237 + (t68 / 0.2e1 + t44 * t387 + t43 * t389) * t280 + (t84 / 0.2e1 + t55 / 0.2e1 - t142 / 0.2e1) * t238 + (t42 / 0.2e1 + t8 / 0.2e1 - t67 / 0.2e1) * t279; -0.2e1 * t413 * t160 + 0.2e1 * (-t179 * t279 + t180 * t280 + t237 * t413 - t238 * t250) * mrSges(5,3) - (-t189 * t314 + t190 * t316 + t241) * t237 + (t188 + t136 - t240) * t238 + (t139 * t204 + t19 * t74 + t20 * t73) * t407 + t324 * t283 + t320 * t284 + 0.2e1 * t309 * t171 - 0.2e1 * pkin(2) * t282 + 0.2e1 * t180 * t221 + 0.2e1 * t81 * t225 + 0.2e1 * t80 * t226 - t210 * t59 - t211 * t60 + 0.2e1 * t204 * t65 + 0.2e1 * t19 * t181 + 0.2e1 * t20 * t182 + 0.2e1 * t166 * t167 + 0.2e1 * t165 * t168 + 0.2e1 * t139 * t149 + t117 * t137 + t116 * t138 + 0.2e1 * t73 * t89 + 0.2e1 * t74 * t90 + (t165 * t80 + t166 * t81 - t369) * t408 + (t179 * t250 + t309 * t345 - t369) * t409 + (t130 + t58 - t172) * t279 + (-t131 * t314 + t132 * t316 + t173) * t280 + (t324 * t294 + (0.2e1 * pkin(3) * t239 - t293) * t320) * qJD(3); t327 + (t72 * t386 * qJD(4) + (m(5) * t31 + t111 + (m(5) * t94 + t191) * qJD(4)) * t323 + (m(5) * t30 + t112 + (-m(5) * t93 + m(6) * t88 + t356 + t66) * qJD(4)) * t319) * pkin(3) + (t300 * t134 + t305 * t63) * t316 + (-t135 * t300 - t305 * t64 - t385) * t314 + m(6) * (-t12 * t360 + t13 * t359 + t28 * t308 + t361 * t54 - t362 * t53) + m(7) * (t16 * t164 + t163 * t17 + t2 * t215 + t214 * t3 + t22 * t286) + t308 * t61 + t286 * t14 + t214 * t23 + t215 * t24 + t163 * t76 + t164 * t77 + t141 * mrSges(4,1) - t140 * mrSges(4,2) + t411; t326 + (m(5) * (t179 * t319 - t180 * t323) + (t237 * t323 - t238 * t319) * mrSges(5,3) + (-t323 * t279 * mrSges(5,3) + (t280 * mrSges(5,3) + t149 + t221) * t319 + m(5) * (t250 * t323 - t366) - m(6) * t366 + t204 * t386) * qJD(4)) * pkin(3) + t310 + t308 * t160 + t286 * t65 + (-Ifges(4,6) * t320 + (-mrSges(4,1) * t324 + mrSges(4,2) * t320) * pkin(9)) * qJD(3) + t214 * t89 + t215 * t90 + t163 * t181 + t164 * t182 + (t167 * t305 + t225 * t300) * t316 + (-t305 * t168 - t300 * t226 - t382) * t314 + m(6) * (-t165 * t362 + t166 * t361 + t180 * t308 + t359 * t81 - t360 * t80) + m(7) * (t139 * t286 + t163 * t74 + t164 * t73 + t19 * t215 + t20 * t214); t286 * t405 + t300 * t330 + 0.2e1 * t410 + (t163 * t215 + t164 * t214 + t286 * t344) * t407 + (t305 * t338 + t308 * t344) * t408 + (-t163 * t328 - t164 * t278 + t214 * t268 - t215 * t269) * t348 + t336; t327 + m(6) * (-pkin(4) * t28 + (-t314 * t53 + t316 * t54) * qJD(5) + (-t12 * t314 + t373) * qJ(5)) + (qJ(5) * t63 + qJD(5) * t134) * t316 + (-qJ(5) * t64 - qJD(5) * t135 - t385) * t314 + m(7) * (t16 * t199 + t17 * t198 + t2 * t246 + t22 * t306 + t245 * t3) + t306 * t14 + t245 * t23 + t246 * t24 + t198 * t76 + t199 * t77 - pkin(4) * t61; t326 + m(6) * (-pkin(4) * t180 + (-t165 * t314 + t166 * t316) * qJD(5) + (-t314 * t80 + t370) * qJ(5)) + (qJ(5) * t167 + qJD(5) * t225) * t316 + (-qJ(5) * t168 - qJD(5) * t226 - t382) * t314 + m(7) * (t139 * t306 + t19 * t246 + t198 * t74 + t199 * t73 + t20 * t245) + t306 * t65 + t245 * t89 + t246 * t90 + t198 * t181 + t199 * t182 - pkin(4) * t160; (t286 + t306) * t217 + t410 + m(7) * (t163 * t246 + t164 * t245 + t198 * t215 + t199 * t214 + t306 * t344) + m(6) * (-pkin(4) * t344 + qJ(5) * t338 + t305 * t337) + (t338 + t337) * mrSges(6,3) + ((-t164 - t199) * t278 - (t163 + t198) * t328 - (t215 + t246) * t269 - (-t214 - t245) * t268) * mrSges(7,3) + t336; t306 * t405 + (t198 * t246 + t199 * t245) * t407 + (qJ(5) * t355 * t408 + t330) * qJD(5) + (-t198 * t328 - t199 * t278 + t245 * t268 - t246 * t269) * t348 + t336; m(6) * t28 + m(7) * t22 + t14 + t61; m(6) * t180 + m(7) * t139 + t160 + t65; (m(6) + m(7)) * t344 + t217; t217; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t20 - mrSges(7,2) * t19 + t58; mrSges(7,1) * t164 - mrSges(7,2) * t163 + t218; mrSges(7,1) * t199 - mrSges(7,2) * t198 + t218; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
