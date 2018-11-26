% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2018-11-23 16:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:09:22
% EndTime: 2018-11-23 16:09:41
% DurationCPUTime: 19.55s
% Computational Cost: add. (29724->739), mult. (98867->1047), div. (0->0), fcn. (85744->14), ass. (0->351)
t293 = sin(qJ(3));
t286 = sin(pkin(12));
t288 = sin(pkin(6));
t363 = qJD(1) * t288;
t355 = t286 * t363;
t290 = cos(pkin(7));
t296 = cos(qJ(3));
t384 = cos(pkin(12));
t344 = t384 * t296;
t335 = t290 * t344;
t287 = sin(pkin(7));
t385 = cos(pkin(6));
t346 = t385 * t287;
t336 = t296 * t346;
t485 = -t288 * t335 - t336;
t235 = qJD(1) * t485 + t293 * t355;
t285 = sin(pkin(13));
t289 = cos(pkin(13));
t292 = sin(qJ(5));
t295 = cos(qJ(5));
t269 = t285 * t295 + t289 * t292;
t178 = t269 * t235;
t263 = t269 * qJD(5);
t503 = t178 + t263;
t268 = t285 * t292 - t295 * t289;
t179 = t268 * t235;
t262 = t268 * qJD(5);
t366 = t179 + t262;
t345 = t384 * t293;
t373 = t286 * t296;
t243 = t293 * t346 + (t290 * t345 + t373) * t288;
t236 = t243 * qJD(1);
t347 = t288 * t384;
t259 = -t287 * t347 + t290 * t385;
t255 = qJD(1) * t259 + qJD(3);
t194 = t236 * t289 + t255 * t285;
t476 = t236 / 0.2e1;
t501 = -Ifges(4,6) / 0.2e1;
t502 = Ifges(4,4) * t476 - t255 * t501 - t194 * Ifges(5,5) / 0.2e1;
t230 = qJD(5) + t235;
t500 = -t230 / 0.2e1;
t341 = -t236 * t285 + t289 * t255;
t480 = t194 * t295 + t292 * t341;
t499 = -t480 / 0.2e1;
t481 = -t194 * t292 + t295 * t341;
t498 = -t481 / 0.2e1;
t415 = pkin(10) + qJ(4);
t274 = t415 * t285;
t275 = t415 * t289;
t307 = -t295 * t274 - t275 * t292;
t356 = pkin(1) * t385;
t281 = t286 * t356;
t334 = qJD(1) * t347;
t258 = qJ(2) * t334 + qJD(1) * t281;
t299 = (t290 * t347 + t346) * pkin(9);
t224 = qJD(1) * t299 + t258;
t282 = t384 * t356;
t279 = qJD(1) * t282;
t375 = t286 * t288;
t300 = t385 * pkin(2) + (-pkin(9) * t290 - qJ(2)) * t375;
t231 = qJD(1) * t300 + t279;
t254 = (-pkin(9) * t286 * t287 - pkin(2) * t384 - pkin(1)) * t288;
t247 = qJD(1) * t254 + qJD(2);
t310 = t231 * t290 + t247 * t287;
t162 = -t293 * t224 + t296 * t310;
t183 = pkin(3) * t236 + qJ(4) * t235;
t103 = -t162 * t285 + t289 * t183;
t83 = pkin(10) * t235 * t289 + pkin(4) * t236 + t103;
t104 = t289 * t162 + t285 * t183;
t378 = t235 * t285;
t93 = pkin(10) * t378 + t104;
t464 = -qJD(4) * t268 + qJD(5) * t307 - t292 * t83 - t295 * t93;
t423 = t289 / 0.2e1;
t496 = -pkin(11) * t236 + t464;
t214 = t296 * t224;
t370 = t290 * t293;
t372 = t287 * t293;
t163 = t231 * t370 + t247 * t372 + t214;
t120 = -pkin(4) * t378 + t163;
t495 = pkin(5) * t503 + t366 * pkin(11) - t120;
t252 = (-t286 * t370 + t344) * t363;
t340 = t287 * t355;
t210 = -t252 * t285 + t289 * t340;
t211 = t252 * t289 + t285 * t340;
t260 = -t285 * t372 + t289 * t290;
t261 = t285 * t290 + t289 * t372;
t308 = t295 * t260 - t261 * t292;
t360 = qJD(3) * t296;
t352 = t287 * t360;
t494 = -qJD(5) * t308 + t210 * t292 + t211 * t295 + t268 * t352;
t301 = t288 * (t290 * t373 + t345);
t251 = qJD(1) * t301;
t361 = qJD(3) * t293;
t353 = t287 * t361;
t493 = -t251 + t353;
t136 = qJD(6) - t481;
t131 = -t255 * pkin(3) + qJD(4) - t162;
t185 = -t231 * t287 + t290 * t247;
t410 = Ifges(5,4) * t289;
t325 = -t285 * Ifges(5,2) + t410;
t424 = -t285 / 0.2e1;
t482 = t255 * Ifges(4,5) / 0.2e1;
t489 = -Ifges(4,4) / 0.2e1;
t452 = t235 * t489 + t482;
t473 = t341 / 0.2e1;
t411 = Ifges(5,4) * t285;
t474 = Ifges(5,1) * t423 - t411 / 0.2e1;
t486 = Ifges(4,1) * t476;
t491 = t194 * t474 + t325 * t473 + (t194 * Ifges(5,1) + Ifges(5,4) * t341 + Ifges(5,5) * t235) * t423 + (t194 * Ifges(5,4) + Ifges(5,2) * t341 + Ifges(5,6) * t235) * t424 - t162 * mrSges(4,3) + t452 + t486 + t131 * (mrSges(5,1) * t285 + mrSges(5,2) * t289) + t185 * mrSges(4,2);
t490 = Ifges(6,3) * t500 - t341 * Ifges(5,6) / 0.2e1 + Ifges(6,5) * t499 + Ifges(6,6) * t498 + t502;
t374 = t286 * t293;
t237 = (t336 + (t335 - t374) * t288) * qJD(3);
t222 = qJD(1) * t237;
t98 = qJD(5) * t481 - t222 * t268;
t488 = t98 / 0.2e1;
t99 = qJD(5) * t480 + t222 * t269;
t487 = -t99 / 0.2e1;
t238 = t243 * qJD(3);
t223 = qJD(1) * t238;
t428 = t223 / 0.2e1;
t112 = mrSges(6,1) * t230 - mrSges(6,3) * t480;
t291 = sin(qJ(6));
t294 = cos(qJ(6));
t108 = t230 * t294 - t291 * t480;
t109 = t230 * t291 + t294 * t480;
t60 = -mrSges(7,1) * t108 + mrSges(7,2) * t109;
t386 = t112 - t60;
t249 = -t274 * t292 + t275 * t295;
t465 = -qJD(4) * t269 - qJD(5) * t249 + t292 * t93 - t295 * t83;
t298 = qJD(2) * t301;
t125 = qJD(1) * t298 + (t293 * t310 + t214) * qJD(3);
t380 = t222 * t285;
t105 = pkin(4) * t380 + t125;
t51 = -qJD(6) * t109 + t223 * t294 - t291 * t98;
t48 = Ifges(7,6) * t51;
t50 = qJD(6) * t108 + t223 * t291 + t294 * t98;
t49 = Ifges(7,5) * t50;
t15 = Ifges(7,3) * t99 + t48 + t49;
t130 = pkin(3) * t235 - qJ(4) * t236 + t185;
t132 = qJ(4) * t255 + t163;
t80 = t289 * t130 - t132 * t285;
t57 = pkin(4) * t235 - pkin(10) * t194 + t80;
t81 = t285 * t130 + t289 * t132;
t59 = pkin(10) * t341 + t81;
t28 = t292 * t57 + t295 * t59;
t26 = pkin(11) * t230 + t28;
t100 = -pkin(4) * t341 + t131;
t46 = -pkin(5) * t481 - pkin(11) * t480 + t100;
t13 = -t26 * t291 + t294 * t46;
t35 = t99 * pkin(5) - t98 * pkin(11) + t105;
t358 = qJD(5) * t295;
t359 = qJD(5) * t292;
t379 = t222 * t289;
t333 = qJD(2) * t347;
t273 = t296 * t333;
t362 = qJD(2) * t288;
t354 = t286 * t362;
t338 = qJD(1) * t354;
t305 = t290 * t338;
t351 = t290 * t360;
t124 = t231 * t351 + t247 * t352 + qJD(1) * t273 + (-qJD(3) * t224 - t305) * t293;
t113 = qJD(4) * t255 + t124;
t306 = t287 * t338;
t143 = pkin(3) * t223 - qJ(4) * t222 - qJD(4) * t236 + t306;
t75 = -t113 * t285 + t289 * t143;
t58 = pkin(4) * t223 - pkin(10) * t379 + t75;
t76 = t289 * t113 + t285 * t143;
t62 = -pkin(10) * t380 + t76;
t7 = t292 * t58 + t295 * t62 + t57 * t358 - t359 * t59;
t5 = pkin(11) * t223 + t7;
t1 = qJD(6) * t13 + t291 * t35 + t294 * t5;
t14 = t26 * t294 + t291 * t46;
t2 = -qJD(6) * t14 - t291 * t5 + t294 * t35;
t337 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t472 = -t98 * Ifges(6,4) / 0.2e1;
t451 = -t223 * Ifges(6,6) / 0.2e1 + t472;
t439 = t99 / 0.2e1;
t453 = Ifges(6,2) * t439;
t484 = t337 + t105 * mrSges(6,1) - t7 * mrSges(6,3) + t15 / 0.2e1 + t451 + t453;
t483 = Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t27 = -t292 * t59 + t295 * t57;
t135 = Ifges(6,4) * t481;
t392 = t230 * Ifges(6,5);
t79 = Ifges(6,1) * t480 + t135 + t392;
t479 = t100 * mrSges(6,2) + t79 / 0.2e1 - t27 * mrSges(6,3);
t477 = -t185 * mrSges(4,1) - t80 * mrSges(5,1) - t27 * mrSges(6,1) + t81 * mrSges(5,2) + t28 * mrSges(6,2) + t163 * mrSges(4,3) + t490 - (Ifges(4,2) + Ifges(5,3)) * t235 / 0.2e1;
t475 = t259 / 0.2e1;
t22 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t88 = mrSges(6,1) * t223 - mrSges(6,3) * t98;
t471 = t22 - t88;
t175 = mrSges(5,1) * t380 + mrSges(5,2) * t379;
t47 = t99 * mrSges(6,1) + t98 * mrSges(6,2);
t469 = -t175 - t47;
t284 = -pkin(4) * t289 - pkin(3);
t234 = pkin(5) * t268 - pkin(11) * t269 + t284;
t187 = t234 * t291 + t249 * t294;
t468 = -qJD(6) * t187 - t496 * t291 + t495 * t294;
t186 = t234 * t294 - t249 * t291;
t467 = qJD(6) * t186 + t495 * t291 + t496 * t294;
t466 = pkin(5) * t236 - t465;
t219 = t260 * t292 + t261 * t295;
t371 = t287 * t296;
t201 = -t291 * t219 - t294 * t371;
t460 = qJD(6) * t201 + t493 * t291 - t494 * t294;
t303 = -t294 * t219 + t291 * t371;
t459 = qJD(6) * t303 + t494 * t291 + t493 * t294;
t458 = -mrSges(4,1) * t255 - mrSges(5,1) * t341 + mrSges(5,2) * t194 + mrSges(4,3) * t236;
t169 = -mrSges(5,2) * t235 + mrSges(5,3) * t341;
t170 = mrSges(5,1) * t235 - mrSges(5,3) * t194;
t455 = t169 * t289 - t170 * t285;
t454 = t1 * t294 - t2 * t291;
t25 = -pkin(5) * t230 - t27;
t319 = t13 * t294 + t14 * t291;
t321 = Ifges(7,5) * t294 - Ifges(7,6) * t291;
t406 = Ifges(7,4) * t294;
t324 = -Ifges(7,2) * t291 + t406;
t407 = Ifges(7,4) * t291;
t327 = Ifges(7,1) * t294 - t407;
t329 = mrSges(7,1) * t291 + mrSges(7,2) * t294;
t421 = t294 / 0.2e1;
t422 = -t291 / 0.2e1;
t433 = t136 / 0.2e1;
t435 = t109 / 0.2e1;
t437 = t108 / 0.2e1;
t408 = Ifges(7,4) * t109;
t44 = Ifges(7,2) * t108 + Ifges(7,6) * t136 + t408;
t106 = Ifges(7,4) * t108;
t45 = Ifges(7,1) * t109 + Ifges(7,5) * t136 + t106;
t450 = -t319 * mrSges(7,3) + t25 * t329 + t321 * t433 + t324 * t437 + t327 * t435 + t421 * t45 + t422 * t44;
t8 = -qJD(5) * t28 - t292 * t62 + t295 * t58;
t364 = qJ(2) * t347 + t281;
t240 = t299 + t364;
t244 = t282 + t300;
t309 = t244 * t290 + t254 * t287;
t167 = -t293 * t240 + t296 * t309;
t200 = t243 * t289 + t259 * t285;
t242 = t288 * t374 + t485;
t195 = -t244 * t287 + t290 * t254;
t151 = pkin(3) * t242 - qJ(4) * t243 + t195;
t228 = t296 * t240;
t168 = t244 * t370 + t254 * t372 + t228;
t158 = qJ(4) * t259 + t168;
t90 = t289 * t151 - t158 * t285;
t67 = pkin(4) * t242 - pkin(10) * t200 + t90;
t199 = -t243 * t285 + t259 * t289;
t91 = t285 * t151 + t289 * t158;
t74 = pkin(10) * t199 + t91;
t414 = t292 * t67 + t295 * t74;
t376 = t237 * t289;
t145 = t244 * t351 + t254 * t352 + t273 + (-qJD(3) * t240 - t290 * t354) * t293;
t126 = qJD(4) * t259 + t145;
t339 = t287 * t354;
t166 = pkin(3) * t238 - qJ(4) * t237 - qJD(4) * t243 + t339;
t84 = -t126 * t285 + t289 * t166;
t64 = pkin(4) * t238 - pkin(10) * t376 + t84;
t377 = t237 * t285;
t85 = t289 * t126 + t285 * t166;
t73 = -pkin(10) * t377 + t85;
t12 = -qJD(5) * t414 - t292 * t73 + t295 * t64;
t442 = Ifges(6,1) * t488 + Ifges(6,4) * t487 + Ifges(6,5) * t428;
t449 = -m(4) * t162 + m(5) * t131 + t458;
t401 = t136 * Ifges(7,3);
t402 = t109 * Ifges(7,5);
t403 = t108 * Ifges(7,6);
t43 = t401 + t402 + t403;
t391 = t230 * Ifges(6,6);
t409 = Ifges(6,4) * t480;
t78 = Ifges(6,2) * t481 + t391 + t409;
t447 = t14 * mrSges(7,2) + t28 * mrSges(6,3) - t43 / 0.2e1 + t78 / 0.2e1 - t100 * mrSges(6,1) - t13 * mrSges(7,1);
t446 = Ifges(6,2) / 0.2e1;
t17 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + t99 * Ifges(7,5);
t444 = t17 / 0.2e1;
t441 = t50 / 0.2e1;
t440 = t51 / 0.2e1;
t438 = -t108 / 0.2e1;
t436 = -t109 / 0.2e1;
t434 = -t136 / 0.2e1;
t432 = t481 / 0.2e1;
t431 = t480 / 0.2e1;
t430 = Ifges(5,5) * t428 + t222 * t474;
t427 = t230 / 0.2e1;
t425 = t242 / 0.2e1;
t413 = mrSges(4,3) * t222;
t412 = mrSges(4,3) * t223;
t405 = Ifges(5,5) * t289;
t404 = Ifges(5,6) * t285;
t383 = t125 * t296;
t369 = t291 * t262;
t368 = t294 * t262;
t365 = Ifges(4,5) * t222 - Ifges(4,6) * t223;
t357 = Ifges(6,5) * t98 - Ifges(6,6) * t99 + Ifges(6,3) * t223;
t152 = -t179 * t291 + t236 * t294;
t343 = t152 - t369;
t153 = t179 * t294 + t236 * t291;
t342 = -t153 - t368;
t332 = -t1 * t291 - t2 * t294;
t330 = mrSges(7,1) * t294 - mrSges(7,2) * t291;
t326 = Ifges(7,1) * t291 + t406;
t323 = Ifges(7,2) * t294 + t407;
t322 = -t404 + t405;
t320 = Ifges(7,5) * t291 + Ifges(7,6) * t294;
t318 = t13 * t291 - t14 * t294;
t317 = -t285 * t75 + t289 * t76;
t316 = -t285 * t80 + t289 * t81;
t30 = pkin(11) * t242 + t414;
t161 = -t259 * pkin(3) - t167;
t107 = -t199 * pkin(4) + t161;
t165 = t199 * t292 + t200 * t295;
t311 = t295 * t199 - t200 * t292;
t52 = -pkin(5) * t311 - t165 * pkin(11) + t107;
t19 = t291 * t52 + t294 * t30;
t18 = -t291 * t30 + t294 * t52;
t70 = -mrSges(7,2) * t136 + mrSges(7,3) * t108;
t71 = mrSges(7,1) * t136 - mrSges(7,3) * t109;
t315 = -t291 * t71 + t294 * t70;
t33 = -t292 * t74 + t295 * t67;
t122 = t165 * t294 + t242 * t291;
t121 = -t165 * t291 + t242 * t294;
t11 = t292 * t64 + t295 * t73 + t67 * t358 - t359 * t74;
t302 = -(-qJ(2) * t355 + t279) * t286 + t258 * t384;
t264 = (mrSges(3,1) * t385 - mrSges(3,3) * t375) * qJD(1);
t265 = (-mrSges(3,2) * t385 + mrSges(3,3) * t347) * qJD(1);
t146 = t298 + (t293 * t309 + t228) * qJD(3);
t114 = pkin(4) * t377 + t146;
t196 = -mrSges(4,2) * t255 - mrSges(4,3) * t235;
t184 = mrSges(4,1) * t235 + mrSges(4,2) * t236;
t180 = mrSges(4,1) * t223 + mrSges(4,2) * t222;
t177 = mrSges(5,1) * t223 - mrSges(5,3) * t379;
t176 = -mrSges(5,2) * t223 - mrSges(5,3) * t380;
t149 = t223 * Ifges(5,6) + t222 * t325;
t111 = -mrSges(6,2) * t230 + mrSges(6,3) * t481;
t102 = qJD(5) * t165 + t237 * t269;
t101 = qJD(5) * t311 - t237 * t268;
t89 = -mrSges(6,2) * t223 - mrSges(6,3) * t99;
t87 = pkin(5) * t480 - pkin(11) * t481;
t86 = -mrSges(6,1) * t481 + mrSges(6,2) * t480;
t54 = -qJD(6) * t122 - t101 * t291 + t238 * t294;
t53 = qJD(6) * t121 + t101 * t294 + t238 * t291;
t36 = t102 * pkin(5) - t101 * pkin(11) + t114;
t32 = -mrSges(7,2) * t99 + mrSges(7,3) * t51;
t31 = mrSges(7,1) * t99 - mrSges(7,3) * t50;
t29 = -pkin(5) * t242 - t33;
t21 = t27 * t294 + t291 * t87;
t20 = -t27 * t291 + t294 * t87;
t16 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + t99 * Ifges(7,6);
t10 = -pkin(5) * t238 - t12;
t9 = pkin(11) * t238 + t11;
t6 = -pkin(5) * t223 - t8;
t4 = -qJD(6) * t19 - t291 * t9 + t294 * t36;
t3 = qJD(6) * t18 + t291 * t36 + t294 * t9;
t23 = [(t482 + t486 + (t489 + t322 / 0.2e1) * t235 + t491) * t237 + (Ifges(6,5) * t431 + Ifges(5,6) * t473 + Ifges(6,6) * t432 + Ifges(6,3) * t427 + t483 * t235 - t477 - t502) * t238 + m(7) * (t1 * t19 + t10 * t25 + t13 * t4 + t14 * t3 + t18 * t2 + t29 * t6) + (Ifges(5,5) * t200 + Ifges(6,5) * t165 + Ifges(5,6) * t199 + (Ifges(5,3) + Ifges(6,3)) * t242) * t428 + (Ifges(6,1) * t431 + Ifges(6,4) * t432 + Ifges(6,5) * t427 + t479) * t101 + m(5) * (t125 * t161 + t75 * t90 + t76 * t91 + t80 * t84 + t81 * t85) + (Ifges(4,2) * t242 + Ifges(5,3) * t425 + t259 * t501) * t223 + (t105 * t165 - t242 * t7) * mrSges(6,2) - (Ifges(7,5) * t441 - Ifges(6,6) * t428 + Ifges(7,6) * t440 + Ifges(7,3) * t439 + t453 + t472 + t484) * t311 + t449 * t146 - t76 * mrSges(5,2) * t242 + t75 * mrSges(5,1) * t242 + (t1 * t121 - t122 * t2 - t13 * t53 + t14 * t54) * mrSges(7,3) + (Ifges(7,1) * t53 + Ifges(7,4) * t54) * t435 + (Ifges(7,1) * t122 + Ifges(7,4) * t121) * t441 + (Ifges(7,5) * t53 + Ifges(7,6) * t54) * t433 + (Ifges(7,5) * t122 + Ifges(7,6) * t121) * t439 + m(4) * (t124 * t168 - t125 * t167 + t145 * t163 + (qJD(1) * t195 + t185) * t339) + t184 * t339 + (-Ifges(6,4) * t431 + Ifges(7,5) * t435 - Ifges(6,2) * t432 - Ifges(6,6) * t427 + Ifges(7,6) * t437 + Ifges(7,3) * t433 - t447) * t102 + (Ifges(7,4) * t53 + Ifges(7,2) * t54) * t437 + (Ifges(7,4) * t122 + Ifges(7,2) * t121) * t440 + (t199 * t76 - t200 * t75 - t376 * t80 - t377 * t81) * mrSges(5,3) - t168 * t412 - t167 * t413 - (Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * t242) * t380 / 0.2e1 + (Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * t242) * t379 / 0.2e1 + (mrSges(4,2) * t306 + mrSges(4,3) * t125 + Ifges(4,1) * t222 - Ifges(4,4) * t223) * t243 - 0.2e1 * t264 * t354 + t124 * (-mrSges(4,2) * t259 - mrSges(4,3) * t242) + t8 * (mrSges(6,1) * t242 - mrSges(6,3) * t165) + t365 * t475 + (Ifges(6,4) * t165 + Ifges(6,6) * t242) * t487 + (Ifges(6,1) * t165 + Ifges(6,5) * t242) * t488 + t29 * t22 + t18 * t31 + t19 * t32 + t357 * t425 + t200 * t430 + m(6) * (t100 * t114 + t105 * t107 + t11 * t28 + t12 * t27 + t33 * t8 + t414 * t7) + t414 * t89 + t165 * t442 + t122 * t444 + t53 * t45 / 0.2e1 + t54 * t44 / 0.2e1 + t25 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + (-Ifges(4,4) * t242 + Ifges(4,5) * t475 + t322 * t425) * t222 + t10 * t60 + (-t125 * t259 + t242 * t306) * mrSges(4,1) + t3 * t70 + t4 * t71 + t33 * t88 + m(3) * ((t384 * t364 + (qJ(2) * t375 - t282) * t286) * qJD(1) + t302) * t362 + t107 * t47 + t11 * t111 + t12 * t112 + t114 * t86 + t121 * t16 / 0.2e1 + t6 * (-mrSges(7,1) * t121 + mrSges(7,2) * t122) + 0.2e1 * t265 * t333 + t85 * t169 + t84 * t170 + t161 * t175 + t91 * t176 + t90 * t177 + t195 * t180 + t145 * t196 + t199 * t149 / 0.2e1 + t125 * (-mrSges(5,1) * t199 + mrSges(5,2) * t200); -t303 * t32 + (-t1 * t303 + t13 * t459 + t14 * t460 + t2 * t201 - t308 * t6) * m(7) - t471 * t308 + (m(6) * t27 - m(7) * t25 + t386) * (-qJD(5) * t219 - t295 * t210 + t211 * t292 - t269 * t352) + (-m(6) * t100 - t449 - t86) * t251 + (t8 * t308 + t7 * t219 + (t100 * t361 - t105 * t296) * t287 - t494 * t28) * m(6) - t494 * t111 + t264 * t355 - t372 * t412 + t290 * t180 + (-t210 * t80 - t211 * t81 + t75 * t260 + t76 * t261 + (-t383 + (t131 * t293 + t296 * t316) * qJD(3)) * t287) * m(5) - t184 * t340 + t261 * t176 + t260 * t177 - t252 * t196 + (-t163 * t252 - t185 * t340 + (t305 + t124 * t293 - t383 + (-t162 * t293 + t163 * t296) * qJD(3)) * t287) * m(4) + (t196 + t455) * t352 + (t86 + t458) * t353 + t459 * t71 + t460 * t70 + (-t413 + t469) * t371 - m(3) * t302 * t363 - t265 * t334 + t201 * t31 - t210 * t170 - t211 * t169 + t219 * t89; ((t405 / 0.2e1 - t404 / 0.2e1) * t235 + (-t285 * t81 - t289 * t80) * mrSges(5,3) + (Ifges(4,1) / 0.2e1 - t483) * t236 + t452 + t491) * t235 + (t477 + t490) * t236 + (t176 * t289 - t177 * t285) * qJ(4) + (-t8 * mrSges(6,3) + t6 * t329 + t327 * t441 + t324 * t440 + t321 * t439 + 0.2e1 * t442 + t105 * mrSges(6,2) + t17 * t421 + t16 * t422 + t332 * mrSges(7,3) + (t323 * t438 + t326 * t436 + t320 * t434 + t25 * t330 - t294 * t44 / 0.2e1 + t45 * t422 + t318 * mrSges(7,3)) * qJD(6)) * t269 + (t49 / 0.2e1 + t48 / 0.2e1 + (t446 + Ifges(7,3) / 0.2e1) * t99 + t451 + t484) * t268 + (t27 * t366 - t28 * t503) * mrSges(6,3) + (-mrSges(7,2) * t503 - mrSges(7,3) * t343) * t14 + (mrSges(6,1) * t503 - mrSges(6,2) * t366) * t100 + (mrSges(7,1) * t503 - mrSges(7,3) * t342) * t13 + (-Ifges(7,5) * t368 + Ifges(7,6) * t369 + Ifges(7,3) * t263) * t433 + (-Ifges(7,1) * t368 + Ifges(7,4) * t369 + Ifges(7,5) * t263) * t435 + (-Ifges(7,4) * t368 + Ifges(7,2) * t369 + Ifges(7,6) * t263) * t437 + (t369 / 0.2e1 - t152 / 0.2e1) * t44 + (Ifges(7,1) * t153 + Ifges(7,4) * t152 - Ifges(7,5) * t178) * t436 + (Ifges(7,4) * t153 + Ifges(7,2) * t152 - Ifges(7,6) * t178) * t438 + (-t262 / 0.2e1 - t179 / 0.2e1) * t79 + (-Ifges(6,5) * t262 - Ifges(6,6) * t263) * t427 + (-Ifges(6,1) * t262 - Ifges(6,4) * t263) * t431 + (-Ifges(6,4) * t262 - Ifges(6,2) * t263) * t432 + (-t368 / 0.2e1 - t153 / 0.2e1) * t45 + (-mrSges(5,1) * t289 + mrSges(5,2) * t285 - mrSges(4,1)) * t125 + (-t100 * t120 + t105 * t284 + t249 * t7 + t27 * t465 + t28 * t464 + t307 * t8) * m(6) + (t1 * t187 + t13 * t468 + t14 * t467 + t186 * t2 + t25 * t466 - t307 * t6) * m(7) - t471 * t307 + ((Ifges(5,1) * t285 + t410) * t423 + (Ifges(5,2) * t289 + t411) * t424) * t222 + (mrSges(7,1) * t343 + mrSges(7,2) * t342) * t25 + t284 * t47 + t249 * t89 + (-pkin(3) * t125 + qJ(4) * t317 + qJD(4) * t316 - t103 * t80 - t104 * t81 - t131 * t163) * m(5) + t365 + t455 * qJD(4) - t458 * t163 + t464 * t111 + t465 * t112 + t466 * t60 + t467 * t70 + t468 * t71 + (Ifges(6,4) * t179 + Ifges(6,2) * t178) * t498 + (Ifges(6,1) * t179 + Ifges(6,4) * t178) * t499 + (Ifges(6,5) * t179 + Ifges(6,6) * t178) * t500 + t149 * t423 + (Ifges(5,5) * t285 + Ifges(5,6) * t289) * t428 + t285 * t430 + (t43 - t78) * (t263 / 0.2e1 + t178 / 0.2e1) + (Ifges(7,5) * t153 + Ifges(7,6) * t152 - Ifges(7,3) * t178) * t434 - t120 * t86 + t317 * mrSges(5,3) - t124 * mrSges(4,2) - t104 * t169 - t103 * t170 - pkin(3) * t175 + t186 * t31 + t187 * t32 - t162 * t196; -t341 * t169 + t194 * t170 + t291 * t32 + t294 * t31 + t386 * t480 + t315 * qJD(6) + (-t111 - t315) * t481 + (-t136 * t318 - t480 * t25 - t332) * m(7) + (t27 * t480 - t28 * t481 + t105) * m(6) + (t194 * t80 - t341 * t81 + t125) * m(5) - t469; -pkin(5) * t22 + t454 * mrSges(7,3) + t386 * t28 + t357 + (-t403 / 0.2e1 - t402 / 0.2e1 - t401 / 0.2e1 + t391 / 0.2e1 + t409 / 0.2e1 + t447) * t480 + t8 * mrSges(6,1) - t7 * mrSges(6,2) + t16 * t421 + t291 * t444 - t21 * t70 - t20 * t71 - t27 * t111 + (-t392 / 0.2e1 - t135 / 0.2e1 + (t446 - Ifges(6,1) / 0.2e1) * t480 - t450 - t479) * t481 + t450 * qJD(6) + t320 * t439 + t323 * t440 + t326 * t441 - t6 * t330 + (-pkin(5) * t6 - t13 * t20 - t14 * t21 - t25 * t28) * m(7) + (-t291 * t31 + t294 * t32 + m(7) * t454 + (-m(7) * t319 - t291 * t70 - t294 * t71) * qJD(6)) * pkin(11); -t25 * (mrSges(7,1) * t109 + mrSges(7,2) * t108) - t13 * t70 + t14 * t71 + (Ifges(7,1) * t108 - t408) * t436 + t44 * t435 + (Ifges(7,5) * t108 - Ifges(7,6) * t109) * t434 + (t108 * t13 + t109 * t14) * mrSges(7,3) + t337 + t15 + (-Ifges(7,2) * t109 + t106 + t45) * t438;];
tauc  = t23(:);
