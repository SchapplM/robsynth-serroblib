% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:52
% EndTime: 2019-03-09 22:39:50
% DurationCPUTime: 28.52s
% Computational Cost: add. (14506->830), mult. (34744->1134), div. (0->0), fcn. (23790->8), ass. (0->351)
t310 = sin(qJ(3));
t314 = cos(qJ(3));
t311 = sin(qJ(2));
t388 = qJD(1) * t311;
t365 = t314 * t388;
t253 = qJD(2) * t310 + t365;
t309 = sin(qJ(4));
t313 = cos(qJ(4));
t366 = t310 * t388;
t385 = qJD(2) * t314;
t327 = t366 - t385;
t191 = t309 * t253 + t313 * t327;
t308 = sin(qJ(6));
t312 = cos(qJ(6));
t320 = t313 * t253 - t309 * t327;
t118 = t191 * t312 - t308 * t320;
t315 = cos(qJ(2));
t387 = qJD(1) * t315;
t292 = qJD(3) - t387;
t279 = qJD(4) + t292;
t316 = -pkin(4) - pkin(5);
t268 = -pkin(2) * t315 - pkin(8) * t311 - pkin(1);
t245 = t268 * qJD(1);
t304 = pkin(7) * t387;
t276 = qJD(2) * pkin(8) + t304;
t199 = t314 * t245 - t276 * t310;
t162 = -pkin(9) * t253 + t199;
t148 = pkin(3) * t292 + t162;
t200 = t310 * t245 + t314 * t276;
t163 = -pkin(9) * t327 + t200;
t400 = t309 * t163;
t69 = t313 * t148 - t400;
t463 = qJD(5) - t69;
t485 = pkin(10) * t320;
t503 = -t485 + t463;
t53 = t279 * t316 + t503;
t267 = t279 * qJ(5);
t497 = pkin(10) * t191;
t396 = t313 * t163;
t70 = t309 * t148 + t396;
t60 = t70 + t497;
t57 = t267 + t60;
t12 = -t308 * t57 + t312 * t53;
t13 = t308 * t53 + t312 * t57;
t379 = qJD(4) * t313;
t380 = qJD(4) * t309;
t347 = pkin(2) * t311 - pkin(8) * t315;
t262 = t347 * qJD(2);
t246 = qJD(1) * t262;
t378 = qJD(1) * qJD(2);
t358 = t311 * t378;
t351 = pkin(7) * t358;
t136 = -qJD(3) * t200 + t314 * t246 + t310 * t351;
t382 = qJD(3) * t311;
t384 = qJD(2) * t315;
t325 = -t310 * t382 + t314 * t384;
t377 = qJD(2) * qJD(3);
t212 = qJD(1) * t325 + t314 * t377;
t82 = pkin(3) * t358 - pkin(9) * t212 + t136;
t381 = qJD(3) * t314;
t383 = qJD(3) * t310;
t135 = t245 * t381 + t310 * t246 - t276 * t383 - t314 * t351;
t362 = t310 * t384;
t324 = t311 * t381 + t362;
t213 = -qJD(1) * t324 - t310 * t377;
t97 = pkin(9) * t213 + t135;
t19 = -t148 * t380 - t163 * t379 - t309 * t97 + t313 * t82;
t369 = t316 * t311;
t350 = qJD(2) * t369;
t92 = -qJD(4) * t191 + t313 * t212 + t309 * t213;
t10 = -pkin(10) * t92 + qJD(1) * t350 - t19;
t18 = t148 * t379 - t163 * t380 + t309 * t82 + t313 * t97;
t15 = qJ(5) * t358 + t279 * qJD(5) + t18;
t93 = qJD(4) * t320 + t309 * t212 - t313 * t213;
t9 = pkin(10) * t93 + t15;
t2 = qJD(6) * t12 + t10 * t308 + t312 * t9;
t3 = -qJD(6) * t13 + t10 * t312 - t308 * t9;
t32 = qJD(6) * t118 + t308 * t93 + t312 * t92;
t490 = t191 * t308 + t312 * t320;
t33 = -qJD(6) * t490 - t308 * t92 + t312 * t93;
t376 = Ifges(7,5) * t32 + Ifges(7,6) * t33 - Ifges(7,3) * t358;
t412 = Ifges(7,4) * t490;
t271 = qJD(6) - t279;
t427 = -t271 / 0.2e1;
t442 = -t490 / 0.2e1;
t303 = pkin(7) * t388;
t275 = -qJD(2) * pkin(2) + t303;
t218 = pkin(3) * t327 + t275;
t319 = qJ(5) * t320 - t218;
t63 = t191 * t316 + t319;
t521 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t376 + (t118 * t12 + t13 * t490) * mrSges(7,3) + (Ifges(7,5) * t118 - Ifges(7,6) * t490) * t427 + (Ifges(7,1) * t118 - t412) * t442 - t63 * (mrSges(7,1) * t490 + mrSges(7,2) * t118);
t16 = -pkin(4) * t358 - t19;
t443 = t118 / 0.2e1;
t444 = -t118 / 0.2e1;
t480 = -Ifges(5,6) + Ifges(6,6);
t481 = Ifges(6,2) + Ifges(5,3);
t483 = Ifges(6,4) + Ifges(5,5);
t460 = t358 * t481 + t480 * t93 + t483 * t92;
t115 = Ifges(7,4) * t118;
t501 = Ifges(7,2) * t490 - t115;
t55 = Ifges(7,2) * t118 + Ifges(7,6) * t271 + t412;
t56 = Ifges(7,1) * t490 + Ifges(7,5) * t271 + t115;
t520 = t19 * mrSges(5,1) - t16 * mrSges(6,1) - t18 * mrSges(5,2) + t15 * mrSges(6,3) + t442 * t55 + t443 * t56 + t501 * t444 + t460 - t521;
t254 = t309 * t310 - t313 * t314;
t459 = qJD(3) + qJD(4);
t203 = t459 * t254;
t326 = t254 * t315;
t223 = qJD(1) * t326;
t516 = t203 - t223;
t255 = t309 * t314 + t310 * t313;
t204 = t459 * t255;
t222 = t255 * t387;
t391 = t204 - t222;
t80 = t313 * t162 - t400;
t467 = pkin(3) * t379 + qJD(5) - t80;
t259 = t347 * qJD(1);
t214 = pkin(7) * t366 + t314 * t259;
t395 = t314 * t315;
t332 = pkin(3) * t311 - pkin(9) * t395;
t178 = qJD(1) * t332 + t214;
t239 = t310 * t259;
t397 = t311 * t314;
t398 = t310 * t315;
t196 = t239 + (-pkin(7) * t397 - pkin(9) * t398) * qJD(1);
t114 = t309 * t178 + t313 * t196;
t104 = qJ(5) * t388 + t114;
t445 = -pkin(9) - pkin(8);
t367 = qJD(3) * t445;
t260 = t310 * t367;
t261 = t314 * t367;
t277 = t445 * t310;
t278 = t445 * t314;
t144 = t313 * t260 + t309 * t261 + t277 * t379 + t278 * t380;
t518 = -t104 + t144;
t113 = t178 * t313 - t309 * t196;
t211 = t309 * t277 - t313 * t278;
t145 = qJD(4) * t211 + t260 * t309 - t313 * t261;
t517 = -t113 - t145;
t484 = Ifges(5,1) + Ifges(6,1);
t482 = Ifges(6,5) - Ifges(5,4);
t510 = pkin(10) * t516 - qJD(1) * t369 - t517;
t509 = -pkin(10) * t391 - t518;
t505 = -t485 + t467;
t79 = t162 * t309 + t396;
t504 = pkin(3) * t380 - t497 - t79;
t420 = pkin(3) * t310;
t248 = t387 * t420 + t304;
t502 = pkin(3) * t383 - t248;
t185 = Ifges(5,4) * t191;
t409 = Ifges(6,5) * t191;
t461 = t483 * t279 + t320 * t484 - t185 + t409;
t498 = t461 / 0.2e1;
t446 = t93 / 0.2e1;
t496 = Ifges(6,3) * t446;
t401 = qJ(5) * t191;
t495 = qJ(5) * t516 - qJD(5) * t255 + t502;
t454 = t32 / 0.2e1;
t453 = t33 / 0.2e1;
t448 = t92 / 0.2e1;
t433 = t212 / 0.2e1;
t432 = t213 / 0.2e1;
t488 = -t358 / 0.2e1;
t487 = t358 / 0.2e1;
t486 = pkin(4) * t320;
t210 = -t313 * t277 - t278 * t309;
t173 = -pkin(10) * t255 + t210;
t174 = pkin(10) * t254 + t211;
t101 = t173 * t308 + t174 * t312;
t479 = -qJD(6) * t101 + t308 * t509 + t312 * t510;
t100 = t173 * t312 - t174 * t308;
t478 = qJD(6) * t100 + t308 * t510 - t312 * t509;
t477 = t483 * t358 + t482 * t93 + t484 * t92;
t265 = t312 * qJ(5) + t308 * t316;
t472 = -qJD(6) * t265 - t308 * t503 - t312 * t60;
t264 = -t308 * qJ(5) + t312 * t316;
t471 = qJD(6) * t264 - t308 * t60 + t312 * t503;
t300 = -pkin(3) * t313 - pkin(4);
t295 = -pkin(5) + t300;
t296 = pkin(3) * t309 + qJ(5);
t225 = t295 * t308 + t296 * t312;
t470 = -qJD(6) * t225 - t308 * t505 + t312 * t504;
t224 = t295 * t312 - t296 * t308;
t469 = qJD(6) * t224 + t308 * t504 + t312 * t505;
t468 = t316 * t391 - t495;
t232 = t255 * t311;
t464 = t316 * t320;
t462 = pkin(4) * t391 + t495;
t294 = pkin(7) * t395;
t221 = t310 * t268 + t294;
t416 = Ifges(4,4) * t253;
t176 = -Ifges(4,2) * t327 + Ifges(4,6) * t292 + t416;
t249 = Ifges(4,4) * t327;
t177 = t253 * Ifges(4,1) + t292 * Ifges(4,5) - t249;
t333 = t199 * t314 + t200 * t310;
t408 = Ifges(4,6) * t310;
t339 = Ifges(4,5) * t314 - t408;
t415 = Ifges(4,4) * t310;
t343 = Ifges(4,1) * t314 - t415;
t344 = mrSges(4,1) * t310 + mrSges(4,2) * t314;
t422 = t314 / 0.2e1;
t423 = t292 / 0.2e1;
t429 = t253 / 0.2e1;
t458 = -t333 * mrSges(4,3) + t275 * t344 + t343 * t429 + t339 * t423 - t310 * t176 / 0.2e1 + t177 * t422;
t184 = Ifges(6,5) * t320;
t107 = t279 * Ifges(6,6) + t191 * Ifges(6,3) + t184;
t413 = Ifges(5,4) * t320;
t110 = -t191 * Ifges(5,2) + t279 * Ifges(5,6) + t413;
t66 = t267 + t70;
t84 = t191 * pkin(4) - t319;
t457 = t218 * mrSges(5,1) + t84 * mrSges(6,1) + t107 / 0.2e1 - t110 / 0.2e1 - t66 * mrSges(6,2) - t70 * mrSges(5,3);
t456 = Ifges(7,4) * t454 + Ifges(7,2) * t453 + Ifges(7,6) * t488;
t455 = Ifges(7,1) * t454 + Ifges(7,4) * t453 + Ifges(7,5) * t488;
t452 = Ifges(6,5) * t448 + Ifges(6,6) * t487 + t496;
t451 = -t92 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t446 + Ifges(5,6) * t488;
t447 = -t93 / 0.2e1;
t441 = t490 / 0.2e1;
t439 = Ifges(4,1) * t433 + Ifges(4,4) * t432 + Ifges(4,5) * t487;
t438 = -t191 / 0.2e1;
t437 = t191 / 0.2e1;
t435 = -t320 / 0.2e1;
t434 = t320 / 0.2e1;
t430 = -t253 / 0.2e1;
t426 = t271 / 0.2e1;
t425 = -t279 / 0.2e1;
t424 = t279 / 0.2e1;
t419 = pkin(7) * t310;
t418 = mrSges(5,3) * t191;
t417 = mrSges(5,3) * t320;
t414 = Ifges(4,4) * t314;
t411 = Ifges(4,5) * t253;
t410 = Ifges(4,5) * t310;
t407 = Ifges(4,3) * t292;
t406 = t118 * Ifges(7,6);
t405 = t490 * Ifges(7,5);
t404 = t271 * Ifges(7,3);
t158 = t222 * t312 + t223 * t308;
t195 = t254 * t308 + t255 * t312;
t74 = -qJD(6) * t195 + t203 * t308 + t204 * t312;
t403 = t158 - t74;
t159 = t222 * t308 - t223 * t312;
t194 = t254 * t312 - t255 * t308;
t73 = qJD(6) * t194 - t203 * t312 + t204 * t308;
t402 = t159 - t73;
t399 = t310 * t311;
t164 = -mrSges(6,2) * t191 + mrSges(6,3) * t279;
t165 = -mrSges(5,2) * t279 - t418;
t394 = t164 + t165;
t166 = mrSges(5,1) * t279 - t417;
t167 = -mrSges(6,1) * t279 + mrSges(6,2) * t320;
t393 = t166 - t167;
t252 = t314 * t268;
t198 = -pkin(9) * t397 + t252 + (-pkin(3) - t419) * t315;
t206 = -pkin(9) * t399 + t221;
t131 = t309 * t198 + t313 * t206;
t390 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t327 + t253 * mrSges(4,2) + mrSges(3,3) * t388;
t386 = qJD(2) * t311;
t389 = t314 * t262 + t386 * t419;
t263 = pkin(3) * t399 + t311 * pkin(7);
t305 = pkin(7) * t384;
t371 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t370 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t368 = Ifges(4,5) * t212 + Ifges(4,6) * t213 + Ifges(4,3) * t358;
t219 = pkin(3) * t324 + t305;
t301 = -pkin(3) * t314 - pkin(2);
t357 = t315 * t378;
t354 = t384 / 0.2e1;
t353 = -t382 / 0.2e1;
t189 = -pkin(3) * t213 + pkin(7) * t357;
t130 = t198 * t313 - t309 * t206;
t123 = -qJ(5) * t315 + t131;
t233 = t254 * t311;
t345 = -qJ(5) * t233 - t263;
t124 = t315 * pkin(4) - t130;
t342 = Ifges(4,1) * t310 + t414;
t341 = -Ifges(4,2) * t310 + t414;
t340 = Ifges(4,2) * t314 + t415;
t338 = -pkin(3) * t253 - t401;
t83 = pkin(5) * t315 + pkin(10) * t233 + t124;
t94 = pkin(10) * t232 + t123;
t46 = -t308 * t94 + t312 * t83;
t47 = t308 * t83 + t312 * t94;
t98 = -mrSges(7,2) * t271 + mrSges(7,3) * t118;
t99 = mrSges(7,1) * t271 - mrSges(7,3) * t490;
t336 = -t308 * t99 + t312 * t98;
t335 = t135 * t314 - t136 * t310;
t168 = t232 * t312 + t233 * t308;
t169 = t232 * t308 - t233 * t312;
t77 = -mrSges(6,1) * t358 + t92 * mrSges(6,2);
t331 = qJ(5) * t255 - t301;
t129 = t332 * qJD(2) + (-t294 + (pkin(9) * t311 - t268) * t310) * qJD(3) + t389;
t160 = t310 * t262 + t268 * t381 + (-t311 * t385 - t315 * t383) * pkin(7);
t134 = -pkin(9) * t324 + t160;
t41 = t129 * t313 - t309 * t134 - t198 * t380 - t206 * t379;
t329 = t310 * t341;
t328 = t314 * t341;
t40 = t309 * t129 + t313 * t134 + t198 * t379 - t206 * t380;
t137 = -qJD(2) * t326 - t232 * t459;
t323 = qJ(5) * t137 - qJD(5) * t233 - t219;
t322 = qJ(5) * t92 + qJD(5) * t320 - t189;
t38 = qJ(5) * t386 - qJD(5) * t315 + t40;
t302 = Ifges(3,4) * t387;
t273 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t387;
t237 = Ifges(3,1) * t388 + Ifges(3,5) * qJD(2) + t302;
t236 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t311 + t315 * Ifges(3,2)) * qJD(1);
t220 = -pkin(7) * t398 + t252;
t217 = mrSges(4,1) * t292 - mrSges(4,3) * t253;
t216 = -t292 * mrSges(4,2) - mrSges(4,3) * t327;
t215 = -pkin(7) * t365 + t239;
t188 = -mrSges(4,2) * t358 + mrSges(4,3) * t213;
t187 = mrSges(4,1) * t358 - mrSges(4,3) * t212;
t186 = pkin(4) * t254 - t331;
t175 = -Ifges(4,6) * t327 + t407 + t411;
t161 = -qJD(3) * t221 + t389;
t157 = t254 * t316 + t331;
t150 = pkin(4) * t232 - t345;
t149 = -mrSges(4,1) * t213 + mrSges(4,2) * t212;
t139 = t212 * Ifges(4,4) + t213 * Ifges(4,2) + Ifges(4,6) * t358;
t138 = -t380 * t399 + (t397 * t459 + t362) * t313 + t325 * t309;
t128 = mrSges(5,1) * t191 + mrSges(5,2) * t320;
t127 = mrSges(6,1) * t191 - mrSges(6,3) * t320;
t126 = t401 + t486;
t125 = t232 * t316 + t345;
t109 = Ifges(6,4) * t320 + t279 * Ifges(6,2) + t191 * Ifges(6,6);
t108 = Ifges(5,5) * t320 - t191 * Ifges(5,6) + t279 * Ifges(5,3);
t105 = -pkin(4) * t388 - t113;
t103 = -t338 + t486;
t78 = -mrSges(5,2) * t358 - mrSges(5,3) * t93;
t76 = mrSges(5,1) * t358 - mrSges(5,3) * t92;
t75 = -mrSges(6,2) * t93 + mrSges(6,3) * t358;
t72 = -t401 + t464;
t65 = -pkin(4) * t279 + t463;
t64 = t338 + t464;
t58 = -mrSges(7,1) * t118 + mrSges(7,2) * t490;
t54 = t404 + t405 + t406;
t52 = pkin(4) * t138 - t323;
t51 = -qJD(6) * t169 - t137 * t308 + t138 * t312;
t50 = qJD(6) * t168 + t137 * t312 + t138 * t308;
t49 = mrSges(5,1) * t93 + mrSges(5,2) * t92;
t48 = mrSges(6,1) * t93 - mrSges(6,3) * t92;
t39 = -pkin(4) * t386 - t41;
t35 = t138 * t316 + t323;
t34 = pkin(4) * t93 - t322;
t29 = mrSges(7,2) * t358 + mrSges(7,3) * t33;
t28 = -mrSges(7,1) * t358 - mrSges(7,3) * t32;
t23 = pkin(10) * t138 + t38;
t22 = -pkin(10) * t137 + t350 - t41;
t11 = t316 * t93 + t322;
t8 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t5 = -qJD(6) * t47 + t22 * t312 - t23 * t308;
t4 = qJD(6) * t46 + t22 * t308 + t23 * t312;
t1 = [-t139 * t399 / 0.2e1 + (t314 * t353 - t362 / 0.2e1) * t176 + (-t333 * t384 + (-t135 * t310 - t136 * t314 + (t199 * t310 - t200 * t314) * qJD(3)) * t311) * mrSges(4,3) + (t135 * t315 - t200 * t386 + t275 * t325) * mrSges(4,2) + (t311 * t149 + (t390 * t315 + (t344 * t387 - t273) * t311) * qJD(2)) * pkin(7) + m(7) * (t11 * t125 + t12 * t5 + t13 * t4 + t2 * t47 + t3 * t46 + t35 * t63) + m(6) * (t123 * t15 + t124 * t16 + t150 * t34 + t38 * t66 + t39 * t65 + t52 * t84) + m(5) * (t130 * t19 + t131 * t18 + t189 * t263 + t218 * t219 + t40 * t70 + t41 * t69) + t46 * t28 + t47 * t29 + (-t137 * t84 - t15 * t315 + t233 * t34 + t386 * t66) * mrSges(6,3) + (t137 * t218 + t18 * t315 - t189 * t233 - t386 * t70) * mrSges(5,2) + (-Ifges(5,4) * t233 - Ifges(5,6) * t315) * t447 + (-Ifges(6,5) * t233 - Ifges(6,6) * t315) * t446 + t16 * (mrSges(6,1) * t315 - mrSges(6,2) * t233) + t19 * (-mrSges(5,1) * t315 + mrSges(5,3) * t233) + (t310 * t353 + t314 * t354) * t177 + m(4) * (t135 * t221 + t136 * t220 + t160 * t200 + t161 * t199 + (t275 + t303) * t305) + t315 * t376 / 0.2e1 - t327 * (-t340 * t382 + (Ifges(4,6) * t311 + t315 * t341) * qJD(2)) / 0.2e1 + (mrSges(5,1) * t189 + mrSges(6,1) * t34 - mrSges(6,2) * t15 - mrSges(5,3) * t18 - Ifges(5,2) * t447 + t448 * t482 + t451 + t452 + t496) * t232 + t237 * t354 - (t236 + t54 + qJD(1) * (Ifges(7,5) * t169 + Ifges(7,6) * t168 + Ifges(7,3) * t315)) * t386 / 0.2e1 + t12 * (-mrSges(7,1) * t386 - mrSges(7,3) * t50) + t65 * (-mrSges(6,1) * t386 + mrSges(6,2) * t137) + (-t136 * t315 + t199 * t386 + t275 * t324) * mrSges(4,1) + t13 * (mrSges(7,2) * t386 + mrSges(7,3) * t51) + t69 * (mrSges(5,1) * t386 - mrSges(5,3) * t137) + (Ifges(5,4) * t137 + Ifges(5,6) * t386) * t438 + (-0.2e1 * pkin(1) * (mrSges(3,1) * t311 + mrSges(3,2) * t315) + (-0.3e1 / 0.2e1 * t311 ^ 2 + 0.3e1 / 0.2e1 * t315 ^ 2) * Ifges(3,4)) * t378 + ((t311 * t339 - t483 * t233 + t480 * t232 + (-Ifges(4,3) - t481) * t315) * qJD(1) + t175 + t109 + t108) * t386 / 0.2e1 + (Ifges(6,5) * t137 + Ifges(6,6) * t386) * t437 + t137 * t498 + qJD(2) ^ 2 * (Ifges(3,5) * t315 - Ifges(3,6) * t311) / 0.2e1 + t51 * t55 / 0.2e1 + t50 * t56 / 0.2e1 + t35 * t58 + t63 * (-mrSges(7,1) * t51 + mrSges(7,2) * t50) - (t368 + t460) * t315 / 0.2e1 + t4 * t98 + t5 * t99 + ((-Ifges(4,6) * t314 - t410) * t382 + (Ifges(4,3) * t311 + t315 * t339) * qJD(2)) * t423 + (Ifges(7,5) * t50 + Ifges(7,6) * t51 - Ifges(7,3) * t386) * t426 + (-t342 * t382 + (Ifges(4,5) * t311 + t315 * t343) * qJD(2)) * t429 + (-Ifges(4,6) * t315 + t311 * t341) * t432 + (-Ifges(4,5) * t315 + t311 * t343) * t433 + t397 * t439 + (Ifges(7,1) * t50 + Ifges(7,4) * t51 - Ifges(7,5) * t386) * t441 + (Ifges(7,4) * t50 + Ifges(7,2) * t51 - Ifges(7,6) * t386) * t443 + (Ifges(7,4) * t169 + Ifges(7,2) * t168 + Ifges(7,6) * t315) * t453 + (Ifges(7,1) * t169 + Ifges(7,4) * t168 + Ifges(7,5) * t315) * t454 + t169 * t455 + t168 * t456 + t123 * t75 + t124 * t77 + t125 * t8 + t52 * t127 + t130 * t76 + t131 * t78 - t477 * t233 / 0.2e1 + (-Ifges(5,2) * t438 + Ifges(6,3) * t437 + t480 * t424 + t482 * t434 + t457) * t138 + (t137 * t483 + t386 * t481) * t424 + t150 * t48 + (t137 * t484 + t386 * t483) * t434 + (-t233 * t484 - t315 * t483) * t448 + t38 * t164 + t40 * t165 + t41 * t166 + t39 * t167 + t11 * (-mrSges(7,1) * t168 + mrSges(7,2) * t169) + t160 * t216 + t161 * t217 + t219 * t128 + t220 * t187 + t221 * t188 + t263 * t49 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t311 * t357 + t2 * (-mrSges(7,2) * t315 + mrSges(7,3) * t168) + t3 * (mrSges(7,1) * t315 - mrSges(7,3) * t169); (mrSges(5,1) * t391 - mrSges(5,2) * t516) * t218 + (-t15 * t254 + t16 * t255 - t391 * t66 - t516 * t65) * mrSges(6,2) + (mrSges(6,1) * t391 + mrSges(6,3) * t516) * t84 + (-t18 * t254 - t19 * t255 - t391 * t70 + t516 * t69) * mrSges(5,3) + (t73 / 0.2e1 - t159 / 0.2e1) * t56 + (-t187 * t310 + t188 * t314 + (-t216 * t310 - t217 * t314) * qJD(3) + m(4) * (-qJD(3) * t333 + t335)) * pkin(8) + (t77 - t76) * t210 + (t75 + t78) * t211 + (mrSges(7,1) * t403 - mrSges(7,2) * t402) * t63 + (t12 * t402 - t13 * t403 + t194 * t2 - t195 * t3) * mrSges(7,3) - t393 * t145 + t394 * t144 + (t107 - t110) * (t204 / 0.2e1 - t222 / 0.2e1) - m(4) * (t199 * t214 + t200 * t215) + (t74 / 0.2e1 - t158 / 0.2e1) * t55 + (-Ifges(5,4) * t223 - Ifges(6,5) * t203 - Ifges(5,2) * t222 + Ifges(6,3) * t204) * t437 + t461 * (-t203 / 0.2e1 + t223 / 0.2e1) + (t222 * t480 - t223 * t483) * t425 + (t222 * t482 - t223 * t484) * t435 + (-t203 * t483 + t204 * t480) * t424 + (-t203 * t484 + t204 * t482) * t434 + (-Ifges(5,4) * t203 - Ifges(6,5) * t223 - Ifges(5,2) * t204 + Ifges(6,3) * t222) * t438 + t335 * mrSges(4,3) + ((-t237 / 0.2e1 + pkin(1) * mrSges(3,2) * qJD(1) - t302 / 0.2e1 + (Ifges(3,5) / 0.2e1 - t328 / 0.2e1) * qJD(2) + (-m(4) * t275 + (-m(4) * pkin(2) - mrSges(4,1) * t314 + mrSges(4,2) * t310 - mrSges(3,1)) * qJD(2) - t390) * pkin(7) - t458) * t315 + ((-Ifges(7,5) * t195 / 0.2e1 - Ifges(7,6) * t194 / 0.2e1 + t410 / 0.2e1 - Ifges(3,6) / 0.2e1 + pkin(7) * mrSges(3,2) + t371 * t255 + t370 * t254) * qJD(2) + (pkin(1) * mrSges(3,1) + (Ifges(3,4) / 0.2e1 + t408 / 0.2e1) * t311) * qJD(1) - t175 / 0.2e1 - t411 / 0.2e1 - t407 / 0.2e1 + t406 / 0.2e1 + t405 / 0.2e1 - t109 / 0.2e1 + t54 / 0.2e1 - t370 * t191 - t108 / 0.2e1 + t404 / 0.2e1 - t199 * mrSges(4,1) + t200 * mrSges(4,2) + pkin(7) * t273 + t12 * mrSges(7,1) - t13 * mrSges(7,2) + t65 * mrSges(6,1) - t69 * mrSges(5,1) - t66 * mrSges(6,3) + t70 * mrSges(5,2) - t371 * t320 + t236 / 0.2e1 - qJD(3) * t329 / 0.2e1 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t279 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + t329 / 0.2e1) * t387) * t311) * qJD(1) + (qJD(2) * t328 / 0.2e1 + t128 * t420 + t458) * qJD(3) + (t18 * t211 + t189 * t301 - t19 * t210 + (-t114 + t144) * t70 + t517 * t69 + t502 * t218) * m(5) + (t15 * t211 + t16 * t210 + t186 * t34 + t462 * t84 + t518 * t66 + (-t105 + t145) * t65) * m(6) + t462 * t127 + t100 * t28 + t101 * t29 + t468 * t58 + t139 * t422 + (Ifges(7,5) * t73 + Ifges(7,6) * t74) * t426 + (Ifges(7,5) * t159 + Ifges(7,6) * t158) * t427 + t340 * t432 + t342 * t433 + t310 * t439 + (Ifges(7,1) * t73 + Ifges(7,4) * t74) * t441 + (Ifges(7,1) * t159 + Ifges(7,4) * t158) * t442 + (Ifges(7,4) * t73 + Ifges(7,2) * t74) * t443 + (Ifges(7,4) * t159 + Ifges(7,2) * t158) * t444 + (Ifges(6,5) * t255 + Ifges(6,3) * t254) * t446 + (Ifges(5,4) * t255 - Ifges(5,2) * t254) * t447 + t254 * t451 + t254 * t452 + (Ifges(7,4) * t195 + Ifges(7,2) * t194) * t453 + (Ifges(7,1) * t195 + Ifges(7,4) * t194) * t454 + t195 * t455 + t194 * t456 + t477 * t255 / 0.2e1 + t478 * t98 + t479 * t99 + (t100 * t3 + t101 * t2 + t11 * t157 + t12 * t479 + t13 * t478 + t468 * t63) * m(7) - pkin(2) * t149 + (t254 * t482 + t255 * t484) * t448 + t157 * t8 - t104 * t164 - t114 * t165 - t113 * t166 - t105 * t167 + t186 * t48 + t11 * (-mrSges(7,1) * t194 + mrSges(7,2) * t195) - t215 * t216 - t214 * t217 - t248 * t128 + t34 * (mrSges(6,1) * t254 - mrSges(6,3) * t255) + t189 * (mrSges(5,1) * t254 + mrSges(5,2) * t255) + t301 * t49; t520 + t393 * t79 + (mrSges(5,2) * t218 + mrSges(6,2) * t65 - mrSges(5,3) * t69 - mrSges(6,3) * t84 - Ifges(5,4) * t437 - Ifges(6,5) * t438 - t425 * t483 - t435 * t484 + t498) * t191 + t368 - m(5) * (-t69 * t79 + t70 * t80) + (-Ifges(4,2) * t253 + t177 - t249) * t327 / 0.2e1 + (-t199 * t327 + t200 * t253) * mrSges(4,3) - t275 * (t253 * mrSges(4,1) - mrSges(4,2) * t327) - t292 * (-Ifges(4,5) * t327 - Ifges(4,6) * t253) / 0.2e1 - t64 * t58 + t467 * t164 + (-t103 * t84 + t15 * t296 + t16 * t300 + t467 * t66 - t65 * t79) * m(6) + t469 * t98 + t470 * t99 + (t12 * t470 + t13 * t469 + t2 * t225 + t224 * t3 - t63 * t64) * m(7) + t176 * t429 + (-Ifges(4,1) * t327 - t416) * t430 - t103 * t127 + (m(6) * t65 * t380 - t253 * t128 + t309 * t78 + t313 * t76 + (t165 * t313 - t309 * t393) * qJD(4) + (t18 * t309 + t19 * t313 + 0.2e1 * t218 * t430 + t379 * t70 - t380 * t69) * m(5)) * pkin(3) - t135 * mrSges(4,2) + t136 * mrSges(4,1) + (-Ifges(5,2) * t437 + Ifges(6,3) * t438 + t480 * t425 + t482 * t435 - t457) * t320 - t80 * t165 - t199 * t216 + t200 * t217 + t224 * t28 + t225 * t29 + t296 * t75 + t300 * t77; t520 + (t191 * t65 + t320 * t66) * mrSges(6,2) + (Ifges(6,3) * t320 - t409) * t438 - t84 * (mrSges(6,1) * t320 + mrSges(6,3) * t191) - t218 * (mrSges(5,1) * t320 - mrSges(5,2) * t191) + (t393 + t417) * t70 + (-t394 - t418) * t69 + (-Ifges(5,2) * t320 - t185 + t461) * t437 + (-pkin(4) * t16 + qJ(5) * t15 - t126 * t84 + t463 * t66 - t65 * t70) * m(6) - t72 * t58 + qJ(5) * t75 - pkin(4) * t77 + t110 * t434 + t471 * t98 + t472 * t99 + (t12 * t472 + t13 * t471 + t2 * t265 + t264 * t3 - t63 * t72) * m(7) - t126 * t127 + (-t191 * t483 + t320 * t480) * t425 + (-t191 * t484 + t107 + t184 - t413) * t435 + qJD(5) * t164 + t264 * t28 + t265 * t29; t312 * t28 + t308 * t29 + (t127 - t58) * t320 + t336 * qJD(6) + (-t164 - t336) * t279 + t77 + (-t320 * t63 + t2 * t308 + t3 * t312 + t271 * (-t12 * t308 + t13 * t312)) * m(7) + (-t279 * t66 + t320 * t84 + t16) * m(6); t55 * t441 - t12 * t98 + t13 * t99 + (t56 - t501) * t444 + t521;];
tauc  = t1(:);
