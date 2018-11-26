% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:27:08
% EndTime: 2018-11-23 18:27:18
% DurationCPUTime: 10.68s
% Computational Cost: add. (19096->653), mult. (48903->843), div. (0->0), fcn. (35977->8), ass. (0->309)
t481 = Ifges(6,1) + Ifges(7,1);
t475 = Ifges(6,5) + Ifges(7,4);
t297 = sin(qJ(3));
t298 = sin(qJ(2));
t301 = cos(qJ(3));
t302 = cos(qJ(2));
t378 = t301 * t302;
t266 = -t297 * t298 + t378;
t255 = t266 * qJD(1);
t267 = t297 * t302 + t301 * t298;
t256 = t267 * qJD(1);
t296 = sin(qJ(4));
t300 = cos(qJ(4));
t322 = t255 * t296 + t300 * t256;
t482 = Ifges(5,1) * t322 / 0.2e1;
t480 = Ifges(6,6) - Ifges(7,6);
t295 = sin(qJ(5));
t299 = cos(qJ(5));
t294 = qJD(2) + qJD(3);
t351 = qJD(4) + t294;
t190 = t295 * t322 - t299 * t351;
t185 = Ifges(6,4) * t190;
t191 = t295 * t351 + t299 * t322;
t350 = t300 * t255 - t256 * t296;
t206 = qJD(5) - t350;
t403 = Ifges(7,5) * t190;
t459 = t191 * t481 + t475 * t206 - t185 + t403;
t440 = -t190 / 0.2e1;
t447 = -pkin(8) - pkin(7);
t280 = t447 * t302;
t272 = qJD(1) * t280;
t257 = t297 * t272;
t279 = t447 * t298;
t271 = qJD(1) * t279;
t263 = qJD(2) * pkin(2) + t271;
t221 = t301 * t263 + t257;
t249 = t256 * pkin(9);
t188 = t221 - t249;
t176 = pkin(3) * t294 + t188;
t260 = t301 * t272;
t222 = t263 * t297 - t260;
t422 = pkin(9) * t255;
t189 = t222 + t422;
t177 = t296 * t189;
t118 = t300 * t176 - t177;
t289 = -pkin(2) * t302 - pkin(1);
t278 = qJD(1) * t289;
t231 = -t255 * pkin(3) + t278;
t428 = -t299 / 0.2e1;
t429 = t295 / 0.2e1;
t430 = -t295 / 0.2e1;
t107 = -pkin(4) * t351 - t118;
t342 = mrSges(7,1) * t295 - mrSges(7,3) * t299;
t344 = mrSges(6,1) * t295 + mrSges(6,2) * t299;
t53 = t190 * pkin(5) - t191 * qJ(6) + t107;
t472 = t107 * t344 + t53 * t342;
t184 = Ifges(7,5) * t191;
t96 = Ifges(7,6) * t206 + Ifges(7,3) * t190 + t184;
t406 = Ifges(6,4) * t191;
t99 = -Ifges(6,2) * t190 + Ifges(6,6) * t206 + t406;
t479 = t231 * mrSges(5,2) - t118 * mrSges(5,3) + Ifges(5,4) * t350 + Ifges(5,5) * t351 - t459 * t428 - t99 * t429 - t96 * t430 + t472 + t482;
t478 = Ifges(5,2) / 0.2e1;
t404 = Ifges(6,4) * t299;
t337 = -Ifges(6,2) * t295 + t404;
t476 = t337 * t440;
t228 = t294 * t266;
t217 = t228 * qJD(1);
t229 = t294 * t267;
t218 = t229 * qJD(1);
t116 = qJD(4) * t322 + t217 * t296 + t300 * t218;
t115 = qJD(4) * t350 + t217 * t300 - t218 * t296;
t86 = -qJD(5) * t190 + t299 * t115;
t87 = qJD(5) * t191 + t295 * t115;
t474 = (-Ifges(6,4) + Ifges(7,5)) * t87 + t481 * t86 + t475 * t116;
t24 = mrSges(6,1) * t87 + mrSges(6,2) * t86;
t178 = t300 * t189;
t119 = t296 * t176 + t178;
t357 = qJD(2) * t447;
t369 = qJD(3) * t301;
t370 = qJD(3) * t297;
t166 = t357 * t256 + t263 * t369 + t272 * t370;
t123 = -pkin(9) * t218 + t166;
t268 = t297 * t279;
t310 = (t378 * t447 - t268) * qJD(2) * qJD(1);
t167 = -qJD(3) * t222 + t310;
t420 = t217 * pkin(9);
t29 = qJD(4) * t119 + t123 * t296 - t300 * (t167 - t420);
t473 = m(6) * t29 + t24;
t226 = -t271 * t297 + t260;
t192 = t226 - t422;
t227 = t301 * t271 + t257;
t193 = -t249 + t227;
t288 = pkin(2) * t301 + pkin(3);
t367 = qJD(4) * t300;
t368 = qJD(4) * t296;
t383 = t297 * t300;
t458 = -t300 * t192 + t193 * t296 - t288 * t368 - (t297 * t367 + (t296 * t301 + t383) * qJD(3)) * pkin(2);
t471 = t295 * t475 + t299 * t480;
t333 = Ifges(6,5) * t299 - Ifges(6,6) * t295;
t335 = Ifges(7,4) * t299 + Ifges(7,6) * t295;
t470 = t335 + t333;
t401 = Ifges(7,5) * t299;
t469 = t295 * t481 - t401 + t404;
t402 = Ifges(7,5) * t295;
t339 = Ifges(7,1) * t299 + t402;
t405 = Ifges(6,4) * t295;
t341 = Ifges(6,1) * t299 - t405;
t468 = t341 + t339;
t467 = Ifges(5,6) * t351 / 0.2e1;
t108 = pkin(10) * t351 + t119;
t128 = -pkin(4) * t350 - pkin(10) * t322 + t231;
t28 = t300 * t123 + t118 * qJD(4) + (-t263 * t370 + t272 * t369 + t310 - t420) * t296;
t373 = qJD(1) * t298;
t291 = pkin(2) * t373;
t196 = pkin(3) * t218 + qJD(2) * t291;
t35 = pkin(4) * t116 - pkin(10) * t115 + t196;
t365 = qJD(5) * t299;
t366 = qJD(5) * t295;
t6 = -t108 * t366 + t128 * t365 + t299 * t28 + t295 * t35;
t2 = qJ(6) * t116 + qJD(6) * t206 + t6;
t52 = t108 * t299 + t128 * t295;
t7 = -qJD(5) * t52 - t28 * t295 + t299 * t35;
t4 = -pkin(5) * t116 - t7;
t51 = -t108 * t295 + t128 * t299;
t460 = qJD(6) - t51;
t40 = -pkin(5) * t206 + t460;
t466 = t2 * t299 + t295 * t4 + t40 * t365;
t400 = Ifges(5,2) * t350;
t407 = Ifges(5,4) * t322;
t453 = t467 + t407 / 0.2e1;
t465 = -t231 * mrSges(5,1) - t51 * mrSges(6,1) + t40 * mrSges(7,1) + t52 * mrSges(6,2) + t400 / 0.2e1 + t453;
t435 = -t322 / 0.2e1;
t423 = pkin(5) * t322;
t248 = pkin(5) * t366 - qJ(6) * t365 - qJD(6) * t295;
t328 = pkin(5) * t295 - qJ(6) * t299;
t462 = -t328 * t350 - t119 + t248;
t389 = t350 * t299;
t390 = t350 * t295;
t348 = pkin(5) * t390 - qJ(6) * t389;
t461 = -t348 + t248 - t458;
t202 = t322 * qJ(6);
t377 = -mrSges(5,1) * t351 + mrSges(6,1) * t190 + mrSges(6,2) * t191 + mrSges(5,3) * t322;
t133 = t192 * t296 + t193 * t300;
t384 = t296 * t297;
t215 = t288 * t367 + (-t297 * t368 + (t300 * t301 - t384) * qJD(3)) * pkin(2);
t457 = -t133 + t215;
t225 = t266 * t296 + t267 * t300;
t238 = -t266 * pkin(3) + t289;
t321 = t300 * t266 - t267 * t296;
t149 = -pkin(4) * t321 - t225 * pkin(10) + t238;
t232 = t301 * t279 + t280 * t297;
t204 = -pkin(9) * t267 + t232;
t233 = -t301 * t280 + t268;
t205 = pkin(9) * t266 + t233;
t154 = t204 * t296 + t205 * t300;
t456 = t295 * t149 + t299 * t154;
t455 = t300 * t204 - t205 * t296;
t41 = qJ(6) * t206 + t52;
t361 = t41 * t366;
t417 = t7 * t295;
t418 = t299 * t6;
t454 = m(7) * (-t361 + t466) + m(6) * (-t365 * t51 - t366 * t52 - t417 + t418);
t273 = t298 * t357;
t274 = t302 * t357;
t172 = t301 * t273 + t297 * t274 + t279 * t369 + t280 * t370;
t146 = -pkin(9) * t229 + t172;
t173 = -qJD(3) * t233 - t273 * t297 + t301 * t274;
t147 = -pkin(9) * t228 + t173;
t38 = qJD(4) * t455 + t146 * t300 + t147 * t296;
t135 = qJD(4) * t321 + t228 * t300 - t229 * t296;
t136 = qJD(4) * t225 + t228 * t296 + t300 * t229;
t371 = qJD(2) * t298;
t212 = pkin(2) * t371 + pkin(3) * t229;
t47 = pkin(4) * t136 - pkin(10) * t135 + t212;
t13 = -qJD(5) * t456 - t295 * t38 + t299 * t47;
t164 = pkin(4) * t322 - pkin(10) * t350;
t313 = t206 * t333;
t314 = t206 * t335;
t315 = t191 * t339;
t316 = t191 * t341;
t331 = Ifges(7,3) * t295 + t401;
t317 = t190 * t331;
t324 = t295 * t52 + t299 * t51;
t452 = -t324 * mrSges(6,3) + t476 + t316 / 0.2e1 + t317 / 0.2e1 + t314 / 0.2e1 + t315 / 0.2e1 + t313 / 0.2e1 - (t295 * t41 - t299 * t40) * mrSges(7,2) + t479;
t362 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t363 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t364 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t97 = Ifges(6,5) * t191 - Ifges(6,6) * t190 + Ifges(6,3) * t206;
t98 = Ifges(7,4) * t191 + Ifges(7,2) * t206 + Ifges(7,6) * t190;
t451 = t362 * t190 + t364 * t191 + t363 * t206 - t119 * mrSges(5,3) + t97 / 0.2e1 + t98 / 0.2e1 + t41 * mrSges(7,3) - t453 - t465;
t450 = t86 / 0.2e1;
t449 = -t87 / 0.2e1;
t448 = t87 / 0.2e1;
t445 = pkin(1) * mrSges(3,1);
t444 = pkin(1) * mrSges(3,2);
t443 = t116 / 0.2e1;
t439 = t190 / 0.2e1;
t438 = -t191 / 0.2e1;
t437 = t191 / 0.2e1;
t436 = -t206 / 0.2e1;
t433 = t255 / 0.2e1;
t432 = -t256 / 0.2e1;
t431 = t256 / 0.2e1;
t427 = t299 / 0.2e1;
t426 = m(4) * t278;
t425 = pkin(3) * t256;
t424 = pkin(3) * t300;
t30 = -mrSges(7,2) * t87 + mrSges(7,3) * t116;
t33 = -mrSges(6,2) * t116 - mrSges(6,3) * t87;
t415 = t30 + t33;
t31 = mrSges(6,1) * t116 - mrSges(6,3) * t86;
t32 = -t116 * mrSges(7,1) + t86 * mrSges(7,2);
t414 = -t31 + t32;
t413 = mrSges(4,3) * t255;
t412 = mrSges(6,3) * t190;
t411 = mrSges(6,3) * t191;
t409 = Ifges(3,4) * t298;
t408 = Ifges(4,4) * t256;
t398 = t455 * t29;
t397 = t256 * mrSges(4,3);
t396 = Ifges(3,5) * qJD(2);
t395 = Ifges(3,6) * qJD(2);
t394 = qJD(2) * mrSges(3,1);
t393 = qJD(2) * mrSges(3,2);
t392 = t119 * t322;
t388 = t215 * t295;
t387 = t215 * t299;
t385 = t295 * t300;
t381 = t299 * t300;
t61 = t299 * t118 + t295 * t164;
t125 = t188 * t300 - t177;
t138 = t164 + t425;
t57 = t299 * t125 + t295 * t138;
t134 = t138 + t291;
t59 = t299 * t133 + t295 * t134;
t139 = -mrSges(7,2) * t190 + mrSges(7,3) * t206;
t140 = -mrSges(6,2) * t206 - t412;
t375 = t139 + t140;
t141 = mrSges(6,1) * t206 - t411;
t142 = -mrSges(7,1) * t206 + mrSges(7,2) * t191;
t374 = -t141 + t142;
t251 = pkin(2) * t383 + t296 * t288;
t372 = qJD(1) * t302;
t356 = t396 / 0.2e1;
t355 = -t395 / 0.2e1;
t124 = t188 * t296 + t178;
t250 = -pkin(2) * t384 + t288 * t300;
t345 = mrSges(6,1) * t299 - mrSges(6,2) * t295;
t343 = mrSges(7,1) * t299 + mrSges(7,3) * t295;
t336 = Ifges(6,2) * t299 + t405;
t330 = -Ifges(7,3) * t299 + t402;
t329 = pkin(5) * t299 + qJ(6) * t295;
t60 = -t118 * t295 + t164 * t299;
t56 = -t125 * t295 + t138 * t299;
t58 = -t133 * t295 + t134 * t299;
t75 = t149 * t299 - t154 * t295;
t275 = -pkin(4) - t329;
t12 = t149 * t365 - t154 * t366 + t295 * t47 + t299 * t38;
t312 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t311 = (-qJD(5) * t324 - t417) * mrSges(6,3);
t39 = qJD(4) * t154 + t146 * t296 - t300 * t147;
t11 = pkin(5) * t87 - qJ(6) * t86 - qJD(6) * t191 + t29;
t19 = Ifges(7,5) * t86 + Ifges(7,6) * t116 + Ifges(7,3) * t87;
t20 = Ifges(6,4) * t86 - Ifges(6,2) * t87 + Ifges(6,6) * t116;
t306 = -t28 * mrSges(5,2) + mrSges(6,3) * t418 + Ifges(5,5) * t115 - Ifges(5,6) * t116 - t11 * t343 + t19 * t428 + t20 * t427 + t330 * t448 + t336 * t449 + t469 * t450 + t471 * t443 + t474 * t429 + (-t99 / 0.2e1 + t96 / 0.2e1) * t366 + t459 * t365 / 0.2e1 + (-t345 - mrSges(5,1)) * t29 + (t476 + t472) * qJD(5) + t466 * mrSges(7,2) + (t317 + t316 + t315 + t314 + t313) * qJD(5) / 0.2e1;
t305 = t415 * t299 + t414 * t295 + (-t295 * t375 + t299 * t374) * qJD(5) + t454;
t200 = t255 * Ifges(4,2) + t294 * Ifges(4,6) + t408;
t247 = Ifges(4,4) * t255;
t201 = Ifges(4,1) * t256 + Ifges(4,5) * t294 + t247;
t304 = t306 - t278 * (mrSges(4,1) * t256 + mrSges(4,2) * t255) - t40 * mrSges(7,2) * t389 - t294 * (Ifges(4,5) * t255 - Ifges(4,6) * t256) / 0.2e1 + t200 * t431 + (Ifges(4,1) * t255 - t408) * t432 - t166 * mrSges(4,2) + t167 * mrSges(4,1) + Ifges(4,5) * t217 - Ifges(4,6) * t218 + t221 * t413 - (-Ifges(4,2) * t256 + t201 + t247) * t255 / 0.2e1 + (t389 * t51 + t390 * t52) * mrSges(6,3) + (-t407 + t98 + t97) * t435 + (t467 + Ifges(6,6) * t439 + Ifges(7,6) * t440 + t475 * t438 + (Ifges(6,3) + Ifges(7,2)) * t436 + t465) * t322 + (Ifges(5,1) * t435 + t322 * t478 + t331 * t440 + t337 * t439 + t470 * t436 + t468 * t438 - t479) * t350;
t290 = Ifges(3,4) * t372;
t277 = mrSges(3,3) * t372 - t393;
t276 = -mrSges(3,3) * t373 + t394;
t264 = t275 - t424;
t254 = Ifges(3,1) * t373 + t290 + t396;
t253 = t395 + (Ifges(3,2) * t302 + t409) * qJD(1);
t245 = -pkin(4) - t250;
t237 = pkin(3) * t368 + t248;
t236 = mrSges(4,1) * t294 - t397;
t235 = -mrSges(4,2) * t294 + t413;
t234 = t291 + t425;
t230 = -t250 + t275;
t220 = -mrSges(4,1) * t255 + mrSges(4,2) * t256;
t194 = -mrSges(5,2) * t351 + mrSges(5,3) * t350;
t163 = -mrSges(5,1) * t350 + mrSges(5,2) * t322;
t161 = -mrSges(7,2) * t390 + mrSges(7,3) * t322;
t130 = mrSges(7,1) * t190 - mrSges(7,3) * t191;
t129 = pkin(5) * t191 + qJ(6) * t190;
t113 = Ifges(7,2) * t116;
t111 = Ifges(6,3) * t116;
t90 = t225 * t328 - t455;
t85 = Ifges(7,4) * t86;
t84 = Ifges(6,5) * t86;
t83 = Ifges(6,6) * t87;
t82 = Ifges(7,6) * t87;
t63 = t124 + t348;
t55 = pkin(5) * t321 - t75;
t54 = -qJ(6) * t321 + t456;
t49 = -t60 - t423;
t48 = t61 + t202;
t45 = -t58 - t423;
t44 = t202 + t59;
t43 = -t56 - t423;
t42 = t202 + t57;
t23 = mrSges(7,1) * t87 - mrSges(7,3) * t86;
t14 = t328 * t135 + (qJD(5) * t329 - qJD(6) * t299) * t225 + t39;
t9 = -pkin(5) * t136 - t13;
t8 = qJ(6) * t136 - qJD(6) * t321 + t12;
t1 = [t377 * t39 + (t19 * t429 + t20 * t430 + t196 * mrSges(5,2) + t11 * t342 + t331 * t448 + t337 * t449 + Ifges(5,1) * t115 - Ifges(5,4) * t116 + (mrSges(5,3) + t344) * t29 + (-t295 * t6 - t299 * t7) * mrSges(6,3) + (-t2 * t295 + t299 * t4) * mrSges(7,2) + (t53 * t343 + t107 * t345 + t330 * t440 + t336 * t439 + t99 * t428 + (t295 * t51 - t299 * t52) * mrSges(6,3) + (-t295 * t40 - t299 * t41) * mrSges(7,2) + t469 * t438 + t471 * t436 + t459 * t430) * qJD(5) + t468 * t450 + t470 * t443 + (qJD(5) * t96 + t474) * t427) * t225 + (-t400 / 0.2e1 + t451) * t136 + m(7) * (t11 * t90 + t14 * t53 + t2 * t54 + t4 * t55 + t40 * t9 + t41 * t8) + m(4) * (t166 * t233 + t167 * t232 + t172 * t222 + t173 * t221) + (-t115 * t455 - t116 * t154) * mrSges(5,3) - t455 * t24 + (t166 * t266 - t167 * t267 - t217 * t232 - t218 * t233 - t221 * t228 - t222 * t229) * mrSges(4,3) + (-t266 * t218 - t229 * t433) * Ifges(4,2) + (t266 * t217 - t218 * t267 + t228 * t433 - t229 * t431) * Ifges(4,4) + t289 * (mrSges(4,1) * t218 + mrSges(4,2) * t217) + (t482 + t452) * t135 + m(6) * (t107 * t39 + t12 * t52 + t13 * t51 + t456 * t6 + t7 * t75 - t398) + m(5) * (-t118 * t39 + t119 * t38 + t154 * t28 + t196 * t238 + t212 * t231 - t398) + t456 * t33 - (-t28 * mrSges(5,3) + t84 / 0.2e1 - t83 / 0.2e1 + t111 / 0.2e1 + t85 / 0.2e1 + t113 / 0.2e1 + t82 / 0.2e1 - Ifges(5,4) * t115 + t196 * mrSges(5,1) + t362 * t87 + t364 * t86 + (Ifges(5,2) + t363) * t116 + t312) * t321 + t294 * (Ifges(4,5) * t228 - Ifges(4,6) * t229) / 0.2e1 + t278 * (mrSges(4,1) * t229 + mrSges(4,2) * t228) + (-pkin(7) * t276 + t254 / 0.2e1 + t356 + (-0.2e1 * t444 + 0.3e1 / 0.2e1 * Ifges(3,4) * t302) * qJD(1)) * t302 * qJD(2) + (t217 * t267 + t228 * t431) * Ifges(4,1) + (-pkin(7) * t277 - t253 / 0.2e1 + t355 + (-0.2e1 * t445 - 0.3e1 / 0.2e1 * t409 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t302) * qJD(1) + (t220 + qJD(1) * (-mrSges(4,1) * t266 + mrSges(4,2) * t267) + 0.2e1 * t426) * pkin(2)) * t371 + t54 * t30 + t55 * t32 + t75 * t31 + t90 * t23 + t14 * t130 + t8 * t139 + t12 * t140 + t13 * t141 + t9 * t142 + t38 * t194 + t212 * t163 + t228 * t201 / 0.2e1 - t229 * t200 / 0.2e1 + t172 * t235 + t173 * t236 + t238 * (mrSges(5,1) * t116 + mrSges(5,2) * t115); (t11 * t230 + t461 * t53 + (t387 - t44) * t41 + (t388 - t45) * t40) * m(7) + t461 * t130 + t457 * t194 + (t118 * t458 + t119 * t457 - t231 * t234 - t250 * t29 + t251 * t28) * m(5) + (t245 * t29 + (t387 - t59) * t52 + (-t388 - t58) * t51 - t458 * t107) * m(6) - t377 * t458 + ((-qJD(5) * t375 + t414) * t295 + (qJD(5) * t374 + t415) * t299 + t454) * (pkin(10) + t251) + (-t7 * mrSges(6,3) + t374 * t215 + (-t41 * mrSges(7,2) - t52 * mrSges(6,3)) * qJD(5)) * t295 + (-t51 * qJD(5) * mrSges(6,3) + t375 * t215) * t299 + t304 + ((t235 * t301 - t236 * t297) * qJD(3) + (-t217 * t301 - t218 * t297) * mrSges(4,3)) * pkin(2) + ((t166 * t297 + t167 * t301 + (-t221 * t297 + t222 * t301) * qJD(3)) * pkin(2) - t221 * t226 - t222 * t227) * m(4) + (-t115 * t250 - t116 * t251 + t392) * mrSges(5,3) + ((t356 - t254 / 0.2e1 - t290 / 0.2e1 + qJD(1) * t444 + (t276 - t394) * pkin(7)) * t302 + (t355 + t253 / 0.2e1 + (t445 + t409 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t302) * qJD(1) + (t277 + t393) * pkin(7) + (-t220 - t426) * pkin(2)) * t298) * qJD(1) - t44 * t139 - t59 * t140 - t58 * t141 - t45 * t142 - t41 * t161 + t230 * t23 - t234 * t163 - t227 * t235 - t226 * t236 + t245 * t24 + t222 * t397; t304 + (-t63 + t237) * t130 + (t236 + t397) * t222 + m(7) * (t11 * t264 + t237 * t53) + (-t256 * t163 + (-t115 * t300 - t116 * t296) * mrSges(5,3) + (t377 * t296 + (t295 * t374 + t299 * t375 + t194) * t300 + m(7) * (t381 * t41 + t385 * t40) + m(6) * (t107 * t296 + t381 * t52 - t385 * t51)) * qJD(4) + (0.2e1 * t231 * t432 + t28 * t296 - t29 * t300 + (-t118 * t296 + t119 * t300) * qJD(4)) * m(5)) * pkin(3) + t305 * (pkin(3) * t296 + pkin(10)) - t377 * t124 + (-mrSges(7,2) * t366 - t161) * t41 - m(5) * (-t118 * t124 + t119 * t125) + t311 + t264 * t23 + mrSges(5,3) * t392 - m(7) * (t40 * t43 + t41 * t42 + t53 * t63) - m(6) * (t107 * t124 + t51 * t56 + t52 * t57) - t42 * t139 - t57 * t140 - t56 * t141 - t43 * t142 - t125 * t194 - t221 * t235 + t473 * (-pkin(4) - t424); t306 + ((t478 - Ifges(5,1) / 0.2e1) * t322 - t452) * t350 + t305 * pkin(10) - mrSges(7,2) * t361 - m(6) * (t107 * t119 + t51 * t60 + t52 * t61) + t462 * t130 - t451 * t322 + t311 - t377 * t119 + t275 * t23 - t48 * t139 - t61 * t140 - t60 * t141 - t49 * t142 - t118 * t194 - t473 * pkin(4) + (t11 * t275 - t40 * t49 - t41 * t48 + t462 * t53) * m(7); t312 + t82 - t83 + (-t374 + t411) * t52 + (-t375 - t412) * t51 + qJ(6) * t30 + t99 * t437 + (Ifges(7,3) * t191 - t403) * t440 - pkin(5) * t32 + t84 + t85 - t129 * t130 + qJD(6) * t139 + t111 + t113 - t53 * (mrSges(7,1) * t191 + mrSges(7,3) * t190) - t107 * (mrSges(6,1) * t191 - mrSges(6,2) * t190) + (t190 * t40 + t191 * t41) * mrSges(7,2) + (-t475 * t190 - t191 * t480) * t436 + (-pkin(5) * t4 + qJ(6) * t2 - t129 * t53 - t40 * t52 + t41 * t460) * m(7) + (-Ifges(6,2) * t191 - t185 + t459) * t439 + (-t190 * t481 + t184 - t406 + t96) * t438; t191 * t130 - t206 * t139 + 0.2e1 * (t4 / 0.2e1 + t53 * t437 + t41 * t436) * m(7) + t32;];
tauc  = t1(:);
