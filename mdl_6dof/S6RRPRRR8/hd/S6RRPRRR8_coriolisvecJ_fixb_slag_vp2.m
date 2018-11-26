% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:25:45
% EndTime: 2018-11-23 17:26:05
% DurationCPUTime: 20.00s
% Computational Cost: add. (24786->771), mult. (61542->1095), div. (0->0), fcn. (46885->10), ass. (0->338)
t335 = sin(pkin(11));
t336 = cos(pkin(11));
t339 = sin(qJ(4));
t343 = cos(qJ(4));
t299 = t335 * t343 + t336 * t339;
t344 = cos(qJ(2));
t352 = t299 * t344;
t262 = qJD(1) * t352;
t281 = t299 * qJD(4);
t391 = -t262 + t281;
t357 = t335 * t339 - t336 * t343;
t351 = t357 * t344;
t263 = qJD(1) * t351;
t280 = t357 * qJD(4);
t390 = -t263 + t280;
t340 = sin(qJ(2));
t359 = pkin(2) * t340 - qJ(3) * t344;
t301 = t359 * qJD(1);
t389 = qJD(1) * t340;
t371 = t335 * t389;
t254 = pkin(7) * t371 + t336 * t301;
t396 = t336 * t344;
t356 = pkin(3) * t340 - pkin(8) * t396;
t219 = qJD(1) * t356 + t254;
t282 = t335 * t301;
t397 = t336 * t340;
t398 = t335 * t344;
t353 = -pkin(7) * t397 - pkin(8) * t398;
t239 = qJD(1) * t353 + t282;
t431 = pkin(8) + qJ(3);
t306 = t431 * t335;
t307 = t431 * t336;
t252 = -t339 * t306 + t343 * t307;
t499 = -t299 * qJD(3) - qJD(4) * t252 - t343 * t219 + t239 * t339;
t383 = qJD(4) * t343;
t385 = qJD(3) * t336;
t386 = qJD(3) * t335;
t498 = -t306 * t383 + (-t239 + t385) * t343 + (-qJD(4) * t307 - t219 - t386) * t339;
t527 = -pkin(4) * t389 + pkin(9) * t390 + t499;
t526 = pkin(9) * t391 - t498;
t293 = t336 * qJD(2) - t371;
t370 = t336 * t389;
t378 = t335 * qJD(2);
t294 = t370 + t378;
t236 = t293 * t339 + t294 * t343;
t338 = sin(qJ(5));
t342 = cos(qJ(5));
t364 = t343 * t293 - t294 * t339;
t163 = t236 * t342 + t338 * t364;
t337 = sin(qJ(6));
t341 = cos(qJ(6));
t494 = -t236 * t338 + t342 * t364;
t107 = t163 * t341 + t337 * t494;
t331 = pkin(7) * t389;
t304 = -qJD(2) * pkin(2) + qJD(3) + t331;
t253 = -pkin(3) * t293 + t304;
t181 = -pkin(4) * t364 + t253;
t117 = -pkin(5) * t494 + t181;
t388 = qJD(1) * t344;
t327 = qJD(4) - t388;
t319 = qJD(5) + t327;
t507 = pkin(10) * t163;
t305 = -pkin(2) * t344 - t340 * qJ(3) - pkin(1);
t286 = t305 * qJD(1);
t332 = pkin(7) * t388;
t313 = qJD(2) * qJ(3) + t332;
t241 = t336 * t286 - t335 * t313;
t374 = pkin(3) * t388;
t192 = -t294 * pkin(8) + t241 - t374;
t242 = t335 * t286 + t336 * t313;
t196 = pkin(8) * t293 + t242;
t129 = t343 * t192 - t196 * t339;
t115 = -pkin(9) * t236 + t129;
t111 = pkin(4) * t327 + t115;
t130 = t192 * t339 + t196 * t343;
t116 = pkin(9) * t364 + t130;
t112 = t338 * t116;
t53 = t342 * t111 - t112;
t41 = t53 - t507;
t40 = pkin(5) * t319 + t41;
t489 = pkin(10) * t494;
t114 = t342 * t116;
t54 = t111 * t338 + t114;
t42 = t54 + t489;
t409 = t337 * t42;
t16 = t341 * t40 - t409;
t408 = t341 * t42;
t17 = t337 * t40 + t408;
t517 = -t163 * t337 + t341 * t494;
t348 = qJD(2) * t351;
t184 = -qJD(1) * t348 + qJD(4) * t364;
t349 = qJD(2) * t352;
t185 = -qJD(1) * t349 - qJD(4) * t236;
t86 = qJD(5) * t494 + t184 * t342 + t185 * t338;
t87 = -qJD(5) * t163 - t184 * t338 + t185 * t342;
t31 = qJD(6) * t517 + t337 * t87 + t341 * t86;
t32 = -qJD(6) * t107 - t337 * t86 + t341 * t87;
t377 = qJD(1) * qJD(2);
t368 = t340 * t377;
t376 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t368;
t424 = Ifges(7,4) * t107;
t311 = qJD(6) + t319;
t445 = -t311 / 0.2e1;
t462 = t107 / 0.2e1;
t463 = -t107 / 0.2e1;
t465 = -t517 / 0.2e1;
t48 = Ifges(7,2) * t517 + Ifges(7,6) * t311 + t424;
t381 = qJD(5) * t342;
t382 = qJD(5) * t338;
t276 = qJD(2) * t359 - t340 * qJD(3);
t264 = t276 * qJD(1);
t303 = (qJD(3) - t331) * qJD(2);
t217 = t336 * t264 - t335 * t303;
t350 = t356 * qJD(2);
t188 = qJD(1) * t350 + t217;
t218 = t335 * t264 + t336 * t303;
t367 = t344 * t377;
t363 = t335 * t367;
t193 = -pkin(8) * t363 + t218;
t74 = -qJD(4) * t130 + t343 * t188 - t193 * t339;
t59 = pkin(4) * t368 - pkin(9) * t184 + t74;
t384 = qJD(4) * t339;
t73 = t339 * t188 + t192 * t383 + t343 * t193 - t196 * t384;
t61 = pkin(9) * t185 + t73;
t12 = t111 * t381 - t116 * t382 + t338 * t59 + t342 * t61;
t10 = pkin(10) * t87 + t12;
t13 = -qJD(5) * t54 - t338 * t61 + t342 * t59;
t9 = pkin(5) * t368 - pkin(10) * t86 + t13;
t2 = qJD(6) * t16 + t10 * t341 + t337 * t9;
t3 = -qJD(6) * t17 - t10 * t337 + t341 * t9;
t488 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t99 = Ifges(7,4) * t517;
t49 = Ifges(7,1) * t107 + Ifges(7,5) * t311 + t99;
t525 = (Ifges(7,1) * t517 - t424) * t463 + (Ifges(7,5) * t517 - Ifges(7,6) * t107) * t445 + (t107 * t17 + t16 * t517) * mrSges(7,3) - t117 * (mrSges(7,1) * t107 + mrSges(7,2) * t517) + t48 * t462 + t376 + t488 + (-Ifges(7,2) * t107 + t49 + t99) * t465;
t237 = -t299 * t338 - t342 * t357;
t164 = qJD(5) * t237 - t280 * t342 - t281 * t338;
t195 = -t262 * t338 - t263 * t342;
t393 = t164 - t195;
t238 = t299 * t342 - t338 * t357;
t165 = -qJD(5) * t238 + t280 * t338 - t281 * t342;
t194 = -t262 * t342 + t263 * t338;
t392 = t165 - t194;
t157 = Ifges(6,4) * t494;
t375 = Ifges(6,5) * t86 + Ifges(6,6) * t87 + Ifges(6,3) * t368;
t425 = Ifges(6,4) * t163;
t443 = -t319 / 0.2e1;
t457 = -t163 / 0.2e1;
t459 = -t494 / 0.2e1;
t479 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t96 = Ifges(6,1) * t163 + Ifges(6,5) * t319 + t157;
t524 = t375 + t479 + (Ifges(6,5) * t494 - Ifges(6,6) * t163) * t443 + (t163 * t54 + t494 * t53) * mrSges(6,3) + (-Ifges(6,2) * t163 + t157 + t96) * t459 - t181 * (mrSges(6,1) * t163 + mrSges(6,2) * t494) + (Ifges(6,1) * t494 - t425) * t457 + t525;
t251 = -t343 * t306 - t307 * t339;
t212 = -pkin(9) * t299 + t251;
t213 = -pkin(9) * t357 + t252;
t504 = t212 * t381 - t213 * t382 + t338 * t527 - t526 * t342;
t144 = t338 * t212 + t342 * t213;
t503 = -qJD(5) * t144 + t526 * t338 + t342 * t527;
t521 = t392 * pkin(10) + t504;
t520 = -pkin(5) * t389 - t393 * pkin(10) + t503;
t287 = t335 * t374 + t332;
t497 = t391 * pkin(4) - t287;
t508 = Ifges(3,5) / 0.2e1;
t143 = t342 * t212 - t213 * t338;
t122 = -pkin(10) * t238 + t143;
t123 = pkin(10) * t237 + t144;
t66 = t122 * t337 + t123 * t341;
t506 = -qJD(6) * t66 - t521 * t337 + t520 * t341;
t65 = t122 * t341 - t123 * t337;
t505 = qJD(6) * t65 + t520 * t337 + t521 * t341;
t329 = pkin(4) * t342 + pkin(5);
t379 = qJD(6) * t341;
t380 = qJD(6) * t337;
t394 = t338 * t341;
t57 = -t115 * t338 - t114;
t45 = t57 - t489;
t58 = t342 * t115 - t112;
t46 = t58 - t507;
t502 = t337 * t46 - t341 * t45 - t329 * t380 + (-t338 * t379 + (-t337 * t342 - t394) * qJD(5)) * pkin(4);
t395 = t337 * t338;
t501 = -t337 * t45 - t341 * t46 + t329 * t379 + (-t338 * t380 + (t341 * t342 - t395) * qJD(5)) * pkin(4);
t500 = -t392 * pkin(5) + t497;
t95 = Ifges(6,2) * t494 + Ifges(6,6) * t319 + t425;
t491 = t95 / 0.2e1;
t369 = -Ifges(3,6) * qJD(2) / 0.2e1;
t490 = qJD(2) * t508;
t292 = t336 * t305;
t240 = -pkin(8) * t397 + t292 + (-pkin(7) * t335 - pkin(3)) * t344;
t259 = pkin(7) * t396 + t335 * t305;
t399 = t335 * t340;
t248 = -pkin(8) * t399 + t259;
t171 = t343 * t240 - t339 * t248;
t272 = t357 * t340;
t141 = -pkin(4) * t344 + t272 * pkin(9) + t171;
t172 = t339 * t240 + t343 * t248;
t271 = t299 * t340;
t145 = -pkin(9) * t271 + t172;
t91 = t338 * t141 + t342 * t145;
t476 = -t74 * mrSges(5,1) + t73 * mrSges(5,2);
t387 = qJD(2) * t340;
t330 = Ifges(3,4) * t388;
t421 = Ifges(4,2) * t335;
t427 = Ifges(4,4) * t336;
t360 = -t421 + t427;
t428 = Ifges(4,4) * t335;
t361 = Ifges(4,1) * t336 - t428;
t430 = mrSges(4,2) * t336;
t438 = t336 / 0.2e1;
t439 = -t335 / 0.2e1;
t475 = -(t241 * t336 + t242 * t335) * mrSges(4,3) + t304 * (mrSges(4,1) * t335 + t430) + Ifges(3,1) * t389 / 0.2e1 + t330 / 0.2e1 + t490 + t293 * t360 / 0.2e1 + t294 * t361 / 0.2e1 + (Ifges(4,4) * t294 + Ifges(4,2) * t293 - Ifges(4,6) * t388) * t439 + (Ifges(4,1) * t294 + Ifges(4,4) * t293 - Ifges(4,5) * t388) * t438;
t474 = t31 / 0.2e1;
t473 = t32 / 0.2e1;
t470 = t86 / 0.2e1;
t469 = t87 / 0.2e1;
t468 = pkin(1) * mrSges(3,1);
t467 = pkin(1) * mrSges(3,2);
t464 = t517 / 0.2e1;
t207 = -t271 * t342 + t272 * t338;
t208 = -t271 * t338 - t272 * t342;
t138 = t207 * t341 - t208 * t337;
t461 = t138 / 0.2e1;
t139 = t207 * t337 + t208 * t341;
t460 = t139 / 0.2e1;
t458 = t494 / 0.2e1;
t456 = t163 / 0.2e1;
t455 = t184 / 0.2e1;
t454 = t185 / 0.2e1;
t453 = t207 / 0.2e1;
t452 = t208 / 0.2e1;
t451 = -t364 / 0.2e1;
t450 = t364 / 0.2e1;
t449 = -t236 / 0.2e1;
t448 = t236 / 0.2e1;
t447 = -t271 / 0.2e1;
t446 = -t272 / 0.2e1;
t444 = t311 / 0.2e1;
t442 = t319 / 0.2e1;
t441 = -t327 / 0.2e1;
t440 = t327 / 0.2e1;
t429 = Ifges(3,4) * t340;
t426 = Ifges(5,4) * t236;
t422 = Ifges(4,5) * t336;
t419 = Ifges(4,6) * t335;
t131 = t194 * t341 - t195 * t337;
t167 = t237 * t337 + t238 * t341;
t70 = -qJD(6) * t167 - t164 * t337 + t165 * t341;
t407 = t131 - t70;
t132 = t194 * t337 + t195 * t341;
t166 = t237 * t341 - t238 * t337;
t69 = qJD(6) * t166 + t164 * t341 + t165 * t337;
t406 = t132 - t69;
t403 = qJD(2) * mrSges(3,2);
t266 = mrSges(4,1) * t363 + t367 * t430;
t373 = pkin(7) * t387;
t246 = t336 * t276 + t335 * t373;
t326 = pkin(7) * t367;
t275 = pkin(3) * t363 + t326;
t288 = (pkin(3) * t378 + pkin(7) * qJD(2)) * t344;
t302 = pkin(3) * t399 + t340 * pkin(7);
t372 = Ifges(5,5) * t184 + Ifges(5,6) * t185 + Ifges(5,3) * t368;
t328 = -pkin(3) * t336 - pkin(2);
t39 = -t87 * mrSges(6,1) + t86 * mrSges(6,2);
t8 = -t32 * mrSges(7,1) + t31 * mrSges(7,2);
t124 = -t185 * mrSges(5,1) + t184 * mrSges(5,2);
t90 = t342 * t141 - t338 * t145;
t362 = m(4) * t304 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t293 + mrSges(4,2) * t294 + mrSges(3,3) * t389;
t154 = -pkin(4) * t185 + t275;
t215 = t280 * t340 - t349;
t186 = -pkin(4) * t215 + t288;
t244 = pkin(4) * t271 + t302;
t67 = -pkin(5) * t344 - t208 * pkin(10) + t90;
t68 = pkin(10) * t207 + t91;
t35 = -t337 * t68 + t341 * t67;
t36 = t337 * t67 + t341 * t68;
t265 = pkin(4) * t357 + t328;
t354 = t422 / 0.2e1 - t419 / 0.2e1;
t214 = -t281 * t340 - t348;
t209 = t350 + t246;
t267 = t335 * t276;
t220 = qJD(2) * t353 + t267;
t98 = -qJD(4) * t172 + t343 * t209 - t220 * t339;
t81 = pkin(4) * t387 - pkin(9) * t214 + t98;
t97 = t339 * t209 + t343 * t220 + t240 * t383 - t248 * t384;
t88 = pkin(9) * t215 + t97;
t24 = t141 * t381 - t145 * t382 + t338 * t81 + t342 * t88;
t25 = -qJD(5) * t91 - t338 * t88 + t342 * t81;
t345 = t494 * Ifges(6,6) + t236 * Ifges(5,5) + t163 * Ifges(6,5) + Ifges(4,6) * t293 - Ifges(4,3) * t388 / 0.2e1 + Ifges(4,5) * t294 + t369 - (Ifges(3,2) * t344 + t429) * qJD(1) / 0.2e1 + t517 * Ifges(7,6) + t327 * Ifges(5,3) + t364 * Ifges(5,6) - t242 * mrSges(4,2) + t241 * mrSges(4,1) + t53 * mrSges(6,1) - t54 * mrSges(6,2) + t129 * mrSges(5,1) - t130 * mrSges(5,2) + t16 * mrSges(7,1) - t17 * mrSges(7,2) + t107 * Ifges(7,5) + t319 * Ifges(6,3) + t311 * Ifges(7,3);
t314 = mrSges(3,3) * t388 - t403;
t279 = pkin(4) * t394 + t329 * t337;
t278 = -pkin(4) * t395 + t329 * t341;
t274 = (mrSges(4,1) * t340 - mrSges(4,3) * t396) * t377;
t273 = (-mrSges(4,2) * t340 - mrSges(4,3) * t398) * t377;
t261 = -mrSges(4,1) * t388 - t294 * mrSges(4,3);
t260 = mrSges(4,2) * t388 + t293 * mrSges(4,3);
t258 = -pkin(7) * t398 + t292;
t255 = -pkin(7) * t370 + t282;
t250 = (Ifges(4,5) * t340 + t344 * t361) * t377;
t249 = (Ifges(4,6) * t340 + t344 * t360) * t377;
t247 = -t336 * t373 + t267;
t230 = Ifges(5,4) * t364;
t206 = mrSges(5,1) * t327 - mrSges(5,3) * t236;
t205 = -mrSges(5,2) * t327 + mrSges(5,3) * t364;
t190 = -pkin(5) * t237 + t265;
t170 = -mrSges(5,2) * t368 + mrSges(5,3) * t185;
t169 = mrSges(5,1) * t368 - mrSges(5,3) * t184;
t168 = -mrSges(5,1) * t364 + mrSges(5,2) * t236;
t156 = -pkin(5) * t207 + t244;
t151 = t236 * Ifges(5,1) + t327 * Ifges(5,5) + t230;
t150 = Ifges(5,2) * t364 + t327 * Ifges(5,6) + t426;
t148 = mrSges(6,1) * t319 - mrSges(6,3) * t163;
t147 = -mrSges(6,2) * t319 + mrSges(6,3) * t494;
t125 = pkin(4) * t236 + pkin(5) * t163;
t121 = t184 * Ifges(5,1) + t185 * Ifges(5,4) + Ifges(5,5) * t368;
t120 = t184 * Ifges(5,4) + t185 * Ifges(5,2) + Ifges(5,6) * t368;
t119 = -qJD(5) * t208 - t214 * t338 + t215 * t342;
t118 = qJD(5) * t207 + t214 * t342 + t215 * t338;
t110 = -mrSges(6,1) * t494 + mrSges(6,2) * t163;
t93 = mrSges(7,1) * t311 - mrSges(7,3) * t107;
t92 = -mrSges(7,2) * t311 + mrSges(7,3) * t517;
t89 = -pkin(5) * t119 + t186;
t78 = -mrSges(6,2) * t368 + mrSges(6,3) * t87;
t77 = mrSges(6,1) * t368 - mrSges(6,3) * t86;
t62 = -pkin(5) * t87 + t154;
t50 = -mrSges(7,1) * t517 + mrSges(7,2) * t107;
t44 = -qJD(6) * t139 - t118 * t337 + t119 * t341;
t43 = qJD(6) * t138 + t118 * t341 + t119 * t337;
t38 = t86 * Ifges(6,1) + t87 * Ifges(6,4) + Ifges(6,5) * t368;
t37 = t86 * Ifges(6,4) + t87 * Ifges(6,2) + Ifges(6,6) * t368;
t27 = -mrSges(7,2) * t368 + mrSges(7,3) * t32;
t26 = mrSges(7,1) * t368 - mrSges(7,3) * t31;
t21 = pkin(10) * t119 + t24;
t20 = pkin(5) * t387 - pkin(10) * t118 + t25;
t19 = t341 * t41 - t409;
t18 = -t337 * t41 - t408;
t7 = t31 * Ifges(7,1) + t32 * Ifges(7,4) + Ifges(7,5) * t368;
t6 = t31 * Ifges(7,4) + t32 * Ifges(7,2) + Ifges(7,6) * t368;
t5 = -qJD(6) * t36 + t20 * t341 - t21 * t337;
t4 = qJD(6) * t35 + t20 * t337 + t21 * t341;
t1 = [(t218 * mrSges(4,2) - t217 * mrSges(4,1) - Ifges(5,6) * t454 - Ifges(5,5) * t455 - Ifges(6,6) * t469 - Ifges(6,5) * t470 - Ifges(7,6) * t473 - Ifges(7,5) * t474 + (t362 * pkin(7) + t475 + t490) * qJD(2) + (-0.2e1 * t467 + (0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * t422 + 0.3e1 / 0.2e1 * t419) * t344 + (-0.3e1 / 0.2e1 * Ifges(4,3) - Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(4,1) * t336 ^ 2 / 0.2e1 + (m(4) * pkin(7) + t430) * pkin(7) + (pkin(7) * mrSges(4,1) - t427 + t421 / 0.2e1) * t335) * t340) * t377 + t476 - t479 - t488) * t344 + ((-pkin(7) * t314 + t345 + t369) * qJD(2) + t250 * t438 + t249 * t439 + pkin(7) * t266 + (-t217 * t336 - t218 * t335) * mrSges(4,3) + (Ifges(5,5) * t446 + Ifges(5,6) * t447 + Ifges(6,5) * t452 + Ifges(6,6) * t453 + Ifges(7,5) * t460 + Ifges(7,6) * t461 - 0.2e1 * t468 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t354) * t340) * t377) * t340 + t119 * t491 + (Ifges(7,4) * t139 + Ifges(7,2) * t138) * t473 + (Ifges(6,4) * t208 + Ifges(6,2) * t207) * t469 + (Ifges(7,1) * t139 + Ifges(7,4) * t138) * t474 + (Ifges(6,1) * t208 + Ifges(6,4) * t207) * t470 - (t376 + t375 + t372) * t344 / 0.2e1 + m(5) * (t129 * t98 + t130 * t97 + t171 * t74 + t172 * t73 + t253 * t288 + t275 * t302) + m(6) * (t12 * t91 + t13 * t90 + t154 * t244 + t181 * t186 + t24 * t54 + t25 * t53) + m(7) * (t117 * t89 + t156 * t62 + t16 * t5 + t17 * t4 + t2 * t36 + t3 * t35) + (-Ifges(5,4) * t272 - Ifges(5,2) * t271) * t454 + (-Ifges(5,1) * t272 - Ifges(5,4) * t271) * t455 + (-t129 * t214 + t130 * t215 - t271 * t73 + t272 * t74) * mrSges(5,3) + t275 * (mrSges(5,1) * t271 - mrSges(5,2) * t272) + m(4) * (t217 * t258 + t218 * t259 + t241 * t246 + t242 * t247) + t44 * t48 / 0.2e1 + t43 * t49 / 0.2e1 + t35 * t26 + t36 * t27 + (t138 * t2 - t139 * t3 - t16 * t43 + t17 * t44) * mrSges(7,3) + (-t118 * t53 + t119 * t54 + t12 * t207 - t13 * t208) * mrSges(6,3) + (Ifges(5,5) * t214 + Ifges(5,6) * t215) * t440 + t89 * t50 + t90 * t77 + t91 * t78 + t4 * t92 + t5 * t93 + (Ifges(6,5) * t118 + Ifges(6,6) * t119) * t442 + (Ifges(7,5) * t43 + Ifges(7,6) * t44) * t444 + t121 * t446 + t120 * t447 + (Ifges(5,1) * t214 + Ifges(5,4) * t215) * t448 + (Ifges(5,4) * t214 + Ifges(5,2) * t215) * t450 + t38 * t452 + t37 * t453 + (Ifges(6,1) * t118 + Ifges(6,4) * t119) * t456 + (Ifges(6,4) * t118 + Ifges(6,2) * t119) * t458 + t7 * t460 + t6 * t461 + (Ifges(7,1) * t43 + Ifges(7,4) * t44) * t462 + (Ifges(7,4) * t43 + Ifges(7,2) * t44) * t464 + t117 * (-mrSges(7,1) * t44 + mrSges(7,2) * t43) + t118 * t96 / 0.2e1 + t62 * (-mrSges(7,1) * t138 + mrSges(7,2) * t139) + t24 * t147 + t25 * t148 + t156 * t8 + t171 * t169 + t172 * t170 + t181 * (-mrSges(6,1) * t119 + mrSges(6,2) * t118) + t186 * t110 + t97 * t205 + t98 * t206 + t154 * (-mrSges(6,1) * t207 + mrSges(6,2) * t208) + t214 * t151 / 0.2e1 + t215 * t150 / 0.2e1 + t244 * t39 + t253 * (-mrSges(5,1) * t215 + mrSges(5,2) * t214) + t247 * t260 + t246 * t261 + t259 * t273 + t258 * t274 + t288 * t168 + t302 * t124; t497 * t110 + t498 * t205 + (t129 * t499 + t130 * t498 + t251 * t74 + t252 * t73 - t253 * t287 + t275 * t328) * m(5) + t499 * t206 + t500 * t50 + ((-t330 / 0.2e1 + (t344 * t354 + t467) * qJD(1) + (t508 + (Ifges(4,1) * t335 + t427) * t438 + (Ifges(4,2) * t336 + t428) * t439) * qJD(2) + ((-m(4) * pkin(2) - mrSges(4,1) * t336 + mrSges(4,2) * t335 - mrSges(3,1)) * qJD(2) - t362) * pkin(7) - t475) * t344 + (-t345 + (t314 + t403) * pkin(7) + t369 + (t429 / 0.2e1 + t468 + (-Ifges(3,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t344) * qJD(1)) * t340 + (Ifges(4,5) * t335 + Ifges(5,5) * t299 + Ifges(6,5) * t238 + Ifges(7,5) * t167 + Ifges(4,6) * t336 - Ifges(5,6) * t357 + Ifges(6,6) * t237 + Ifges(7,6) * t166) * t387 / 0.2e1) * qJD(1) + (t70 / 0.2e1 - t131 / 0.2e1) * t48 + (t69 / 0.2e1 - t132 / 0.2e1) * t49 + m(4) * (-t241 * t386 + t242 * t385 + (-t217 * t335 + t218 * t336) * qJ(3)) - m(4) * (t241 * t254 + t242 * t255) + t275 * (mrSges(5,1) * t357 + mrSges(5,2) * t299) - t357 * t120 / 0.2e1 + (t262 / 0.2e1 - t281 / 0.2e1) * t150 + (-Ifges(5,5) * t263 - Ifges(5,6) * t262) * t441 + (-Ifges(5,1) * t263 - Ifges(5,4) * t262) * t449 + (-Ifges(5,4) * t263 - Ifges(5,2) * t262) * t451 + (t263 / 0.2e1 - t280 / 0.2e1) * t151 + t503 * t148 + (t12 * t144 + t13 * t143 + t154 * t265 + t181 * t497 + t503 * t53 + t504 * t54) * m(6) + t504 * t147 + t505 * t92 + (t117 * t500 + t16 * t506 + t17 * t505 + t190 * t62 + t2 * t66 + t3 * t65) * m(7) + t506 * t93 + (t129 * t390 - t130 * t391 - t299 * t74 - t357 * t73) * mrSges(5,3) + (Ifges(5,4) * t299 - Ifges(5,2) * t357) * t454 + (Ifges(5,1) * t299 - Ifges(5,4) * t357) * t455 + (-Ifges(5,5) * t280 - Ifges(5,6) * t281) * t440 + (-Ifges(5,1) * t280 - Ifges(5,4) * t281) * t448 + (-Ifges(5,4) * t280 - Ifges(5,2) * t281) * t450 + (qJ(3) * t273 + qJD(3) * t260 + t218 * mrSges(4,3) + t249 / 0.2e1) * t336 + (-qJ(3) * t274 - qJD(3) * t261 - t217 * mrSges(4,3) + t250 / 0.2e1) * t335 + (t164 / 0.2e1 - t195 / 0.2e1) * t96 + (t165 / 0.2e1 - t194 / 0.2e1) * t95 + t65 * t26 + t66 * t27 + (Ifges(6,5) * t164 + Ifges(6,6) * t165) * t442 + (Ifges(6,5) * t195 + Ifges(6,6) * t194) * t443 + (Ifges(7,5) * t69 + Ifges(7,6) * t70) * t444 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t445 + (Ifges(6,1) * t164 + Ifges(6,4) * t165) * t456 + (Ifges(6,1) * t195 + Ifges(6,4) * t194) * t457 + (Ifges(6,4) * t164 + Ifges(6,2) * t165) * t458 + (Ifges(6,4) * t195 + Ifges(6,2) * t194) * t459 + (Ifges(7,1) * t69 + Ifges(7,4) * t70) * t462 + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t463 + (Ifges(7,4) * t69 + Ifges(7,2) * t70) * t464 + (Ifges(7,4) * t132 + Ifges(7,2) * t131) * t465 + (Ifges(6,4) * t238 + Ifges(6,2) * t237) * t469 + (Ifges(6,1) * t238 + Ifges(6,4) * t237) * t470 + (Ifges(7,4) * t167 + Ifges(7,2) * t166) * t473 + (Ifges(7,1) * t167 + Ifges(7,4) * t166) * t474 + t143 * t77 + t144 * t78 + t166 * t6 / 0.2e1 + t167 * t7 / 0.2e1 + t62 * (-mrSges(7,1) * t166 + mrSges(7,2) * t167) + t190 * t8 + t237 * t37 / 0.2e1 + t238 * t38 / 0.2e1 + t154 * (-mrSges(6,1) * t237 + mrSges(6,2) * t238) + t251 * t169 + t252 * t170 - t255 * t260 - t254 * t261 + t265 * t39 - pkin(2) * t266 - t287 * t168 + t299 * t121 / 0.2e1 + (mrSges(5,1) * t391 - mrSges(5,2) * t390) * t253 + (-mrSges(6,1) * t392 + mrSges(6,2) * t393) * t181 + (t12 * t237 - t13 * t238 + t392 * t54 - t393 * t53) * mrSges(6,3) + t328 * t124 + (mrSges(7,1) * t407 - mrSges(7,2) * t406) * t117 + (t16 * t406 + t166 * t2 - t167 * t3 - t17 * t407) * mrSges(7,3); t107 * t93 - t517 * t92 - t494 * t147 + t163 * t148 - t364 * t205 + t236 * t206 - t293 * t260 + t294 * t261 + t124 + t266 + t39 + t8 + (t107 * t16 - t17 * t517 + t62) * m(7) + (t163 * t53 - t494 * t54 + t154) * m(6) + (t129 * t236 - t130 * t364 + t275) * m(5) + (t241 * t294 - t242 * t293 + t326) * m(4); -t476 + t524 + t501 * t92 + (-t117 * t125 + t16 * t502 + t17 * t501 + t2 * t279 + t278 * t3) * m(7) + t502 * t93 + (-Ifges(5,2) * t236 + t151 + t230) * t451 + t372 + t163 * t491 + (-t236 * t110 + t338 * t78 + t342 * t77 + (t147 * t342 - t148 * t338) * qJD(5) + (t12 * t338 + t13 * t342 - t181 * t236 + t381 * t54 - t382 * t53) * m(6)) * pkin(4) + (t129 * t364 + t130 * t236) * mrSges(5,3) + (Ifges(5,5) * t364 - Ifges(5,6) * t236) * t441 + (Ifges(5,1) * t364 - t426) * t449 - t253 * (mrSges(5,1) * t236 + mrSges(5,2) * t364) + t150 * t448 - t125 * t50 - t58 * t147 - t57 * t148 - m(6) * (t53 * t57 + t54 * t58) - t129 * t205 + t130 * t206 + t278 * t26 + t279 * t27; (-t163 * t50 + t341 * t26 + t337 * t27 + (-t337 * t93 + t341 * t92) * qJD(6) + (-t117 * t163 - t16 * t380 + t17 * t379 + t2 * t337 + t3 * t341) * m(7)) * pkin(5) - m(7) * (t16 * t18 + t17 * t19) - t19 * t92 - t18 * t93 + t95 * t456 - t53 * t147 + t54 * t148 + t524; -t16 * t92 + t17 * t93 + t525;];
tauc  = t1(:);
