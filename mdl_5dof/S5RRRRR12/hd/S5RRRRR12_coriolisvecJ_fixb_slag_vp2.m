% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:37
% EndTime: 2019-12-31 22:47:34
% DurationCPUTime: 24.36s
% Computational Cost: add. (20055->800), mult. (60895->1137), div. (0->0), fcn. (49647->12), ass. (0->391)
t303 = cos(qJ(2));
t295 = cos(pkin(5));
t436 = pkin(1) * t295;
t290 = t303 * t436;
t282 = qJD(1) * t290;
t299 = sin(qJ(2));
t293 = sin(pkin(5));
t294 = cos(pkin(6));
t360 = pkin(9) * t294 + pkin(8);
t344 = t293 * t360;
t322 = t299 * t344;
t217 = -qJD(1) * t322 + t282;
t289 = t299 * t436;
t310 = -t303 * t344 - t289;
t218 = t310 * qJD(1);
t292 = sin(pkin(6));
t435 = pkin(9) * t292;
t312 = (pkin(2) * t299 - t303 * t435) * t293;
t245 = qJD(1) * t312;
t302 = cos(qJ(3));
t298 = sin(qJ(3));
t386 = t294 * t298;
t390 = t292 * t298;
t134 = t302 * t217 + t218 * t386 + t245 * t390;
t285 = pkin(9) * t390;
t385 = t294 * t302;
t267 = pkin(2) * t385 - t285;
t260 = t267 * qJD(3);
t528 = -t134 + t260;
t167 = -t218 * t292 + t294 * t245;
t379 = t299 * t302;
t380 = t298 * t303;
t319 = t294 * t379 + t380;
t371 = qJD(1) * t293;
t235 = t319 * t371;
t377 = t302 * t303;
t381 = t298 * t299;
t317 = -t294 * t381 + t377;
t236 = t317 * t371;
t527 = -pkin(3) * t235 + pkin(10) * t236 - t167 + (pkin(3) * t298 - pkin(10) * t302) * t292 * qJD(3);
t359 = t299 * t371;
t347 = t292 * t359;
t526 = pkin(10) * t347 - t528;
t284 = qJD(1) * t295 + qJD(2);
t320 = t294 * t377 - t381;
t389 = t292 * t302;
t189 = t284 * t389 + t320 * t371;
t187 = qJD(4) - t189;
t318 = t294 * t380 + t379;
t311 = t318 * t293;
t190 = qJD(1) * t311 + t284 * t390;
t358 = t303 * t371;
t238 = t294 * t284 - t292 * t358 + qJD(3);
t297 = sin(qJ(4));
t301 = cos(qJ(4));
t157 = -t190 * t297 + t238 * t301;
t308 = (qJD(2) * t317 + qJD(3) * t320) * t293;
t368 = qJD(3) * t302;
t356 = t292 * t368;
t159 = qJD(1) * t308 + t284 * t356;
t370 = qJD(2) * t293;
t351 = qJD(1) * t370;
t345 = t299 * t351;
t323 = t292 * t345;
t91 = qJD(4) * t157 + t159 * t301 + t297 * t323;
t466 = t91 / 0.2e1;
t525 = Ifges(5,4) * t466;
t182 = t236 * t297 - t301 * t347;
t266 = t294 * t297 + t301 * t390;
t222 = qJD(4) * t266 + t297 * t356;
t524 = t182 - t222;
t183 = t236 * t301 + t297 * t347;
t265 = -t301 * t294 + t297 * t390;
t221 = -qJD(4) * t265 + t301 * t356;
t523 = t183 - t221;
t369 = qJD(3) * t298;
t357 = t292 * t369;
t522 = t235 - t357;
t269 = pkin(2) * t386 + pkin(9) * t389;
t253 = pkin(10) * t294 + t269;
t254 = (-pkin(3) * t302 - pkin(10) * t298 - pkin(2)) * t292;
t366 = qJD(4) * t301;
t367 = qJD(4) * t297;
t492 = -t253 * t367 + t254 * t366 + t527 * t297 - t526 * t301;
t206 = t298 * t217;
t261 = t269 * qJD(3);
t490 = t218 * t385 - t206 - (-pkin(3) * t359 - t245 * t302) * t292 + t261;
t188 = pkin(2) * t284 + t217;
t242 = (-pkin(2) * t303 - t299 * t435 - pkin(1)) * t293;
t231 = qJD(1) * t242;
t151 = -t188 * t292 + t294 * t231;
t97 = -pkin(3) * t189 - pkin(10) * t190 + t151;
t387 = t293 * t303;
t184 = t284 * t435 + (t360 * t387 + t289) * qJD(1);
t327 = t188 * t294 + t231 * t292;
t111 = t184 * t302 + t298 * t327;
t99 = pkin(10) * t238 + t111;
t42 = t297 * t97 + t301 * t99;
t296 = sin(qJ(5));
t300 = cos(qJ(5));
t39 = pkin(11) * t187 + t42;
t158 = t190 * t301 + t238 * t297;
t110 = -t298 * t184 + t302 * t327;
t98 = -pkin(3) * t238 - t110;
t50 = -pkin(4) * t157 - pkin(11) * t158 + t98;
t17 = -t296 * t39 + t300 * t50;
t18 = t296 * t50 + t300 * t39;
t504 = t17 * mrSges(6,1) - t18 * mrSges(6,2);
t480 = t98 * mrSges(5,1) - t42 * mrSges(5,3) + t504;
t153 = qJD(5) - t157;
t409 = Ifges(6,3) * t153;
t117 = -t158 * t296 + t187 * t300;
t410 = Ifges(6,6) * t117;
t118 = t158 * t300 + t187 * t296;
t414 = Ifges(6,5) * t118;
t43 = t409 + t410 + t414;
t411 = Ifges(5,6) * t187;
t413 = Ifges(5,2) * t157;
t421 = Ifges(5,4) * t158;
t84 = t411 + t413 + t421;
t506 = t84 / 0.2e1 - t43 / 0.2e1;
t521 = t506 - t480;
t520 = -pkin(10) * qJD(5) * t301 - t111 + t187 * (pkin(4) * t297 - pkin(11) * t301);
t278 = -pkin(4) * t301 - pkin(11) * t297 - pkin(3);
t145 = pkin(3) * t190 - pkin(10) * t189;
t74 = t301 * t110 + t297 * t145;
t519 = t367 * pkin(10) + pkin(11) * t190 - qJD(5) * t278 + t74;
t92 = qJD(4) * t158 + t159 * t297 - t301 * t323;
t464 = t92 / 0.2e1;
t307 = (qJD(2) * t319 + qJD(3) * t318) * t293;
t160 = qJD(1) * t307 + t284 * t357;
t513 = -t160 / 0.2e1;
t518 = -pkin(11) * t522 + t492;
t517 = -pkin(4) * t524 + pkin(11) * t523 + t490;
t332 = t17 * t300 + t18 * t296;
t334 = Ifges(6,5) * t300 - Ifges(6,6) * t296;
t416 = Ifges(6,4) * t300;
t336 = -Ifges(6,2) * t296 + t416;
t417 = Ifges(6,4) * t296;
t338 = Ifges(6,1) * t300 - t417;
t339 = mrSges(6,1) * t296 + mrSges(6,2) * t300;
t41 = -t297 * t99 + t301 * t97;
t38 = -pkin(4) * t187 - t41;
t438 = t300 / 0.2e1;
t439 = -t296 / 0.2e1;
t418 = Ifges(6,4) * t118;
t44 = Ifges(6,2) * t117 + Ifges(6,6) * t153 + t418;
t114 = Ifges(6,4) * t117;
t45 = Ifges(6,1) * t118 + Ifges(6,5) * t153 + t114;
t452 = t153 / 0.2e1;
t456 = t118 / 0.2e1;
t458 = t117 / 0.2e1;
t516 = -t332 * mrSges(6,3) + t334 * t452 + t336 * t458 + t338 * t456 + t339 * t38 + t438 * t45 + t439 * t44;
t447 = t160 / 0.2e1;
t465 = -t92 / 0.2e1;
t37 = -qJD(5) * t118 + t160 * t300 - t296 * t91;
t472 = t37 / 0.2e1;
t36 = qJD(5) * t117 + t160 * t296 + t300 * t91;
t473 = t36 / 0.2e1;
t220 = t310 * qJD(2);
t204 = qJD(1) * t220;
t246 = qJD(2) * t312;
t239 = qJD(1) * t246;
t279 = qJD(2) * t282;
t314 = qJD(2) * t322;
t203 = -qJD(1) * t314 + t279;
t355 = t294 * t369;
t349 = -t184 * t368 - t188 * t355 - t298 * t203 - t231 * t357;
t61 = -t204 * t385 + (-pkin(3) * t345 - t239 * t302) * t292 - t349;
t19 = pkin(4) * t92 - pkin(11) * t91 + t61;
t354 = t294 * t368;
t63 = -t184 * t369 + t188 * t354 + t302 * t203 + t204 * t386 + t231 * t356 + t239 * t390;
t60 = pkin(10) * t323 + t63;
t163 = -t204 * t292 + t294 * t239;
t77 = pkin(3) * t160 - pkin(10) * t159 + t163;
t10 = t297 * t77 + t301 * t60 + t97 * t366 - t367 * t99;
t8 = pkin(11) * t160 + t10;
t1 = qJD(5) * t17 + t19 * t296 + t300 * t8;
t2 = -qJD(5) * t18 + t19 * t300 - t296 * t8;
t484 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t515 = -t484 - mrSges(5,1) * t61 - Ifges(6,5) * t473 - Ifges(6,6) * t472 - Ifges(6,3) * t464 + t525 - (t513 - t447) * Ifges(5,6) - (t464 - t465) * Ifges(5,2);
t514 = t159 / 0.2e1;
t512 = t238 / 0.2e1;
t511 = Ifges(5,4) * t465;
t510 = t151 * mrSges(4,2);
t509 = t296 * t519 + t300 * t520;
t508 = t296 * t520 - t300 * t519;
t489 = t301 * t253 + t297 * t254;
t493 = -qJD(4) * t489 + t526 * t297 + t527 * t301;
t152 = Ifges(5,4) * t157;
t415 = Ifges(5,5) * t187;
t424 = Ifges(5,1) * t158;
t85 = t152 + t415 + t424;
t467 = t85 / 0.2e1;
t483 = -t98 * mrSges(5,2) + t41 * mrSges(5,3);
t507 = -t152 / 0.2e1 - t467 - t415 / 0.2e1 + t483 - t516;
t505 = -t151 * mrSges(4,1) - t41 * mrSges(5,1) + t42 * mrSges(5,2);
t441 = -t238 / 0.2e1;
t443 = -t190 / 0.2e1;
t503 = Ifges(4,1) * t443 + Ifges(4,5) * t441;
t502 = Ifges(5,1) * t466 + Ifges(5,5) * t447;
t474 = t502 + t511;
t501 = t61 * mrSges(5,2) + t474 + t502;
t444 = -t189 / 0.2e1;
t446 = -t187 / 0.2e1;
t449 = -t158 / 0.2e1;
t451 = -t157 / 0.2e1;
t500 = Ifges(5,5) * t449 - Ifges(4,2) * t444 - Ifges(4,6) * t441 + Ifges(5,6) * t451 + Ifges(5,3) * t446;
t498 = Ifges(4,3) * t512;
t497 = -Ifges(3,6) * t284 / 0.2e1;
t252 = t285 + (-pkin(2) * t302 - pkin(3)) * t294;
t171 = pkin(4) * t265 - pkin(11) * t266 + t252;
t173 = -pkin(11) * t389 + t489;
t115 = t171 * t300 - t173 * t296;
t496 = qJD(5) * t115 + t517 * t296 + t518 * t300;
t116 = t171 * t296 + t173 * t300;
t495 = -qJD(5) * t116 - t518 * t296 + t517 * t300;
t494 = pkin(4) * t522 - t493;
t216 = pkin(2) * t295 + t290 - t322;
t164 = -t216 * t292 + t294 * t242;
t210 = -t293 * t320 - t295 * t389;
t211 = t295 * t390 + t311;
t108 = pkin(3) * t210 - pkin(10) * t211 + t164;
t372 = pkin(8) * t387 + t289;
t205 = (t292 * t295 + t294 * t387) * pkin(9) + t372;
t130 = t302 * t205 + t216 * t386 + t242 * t390;
t264 = -t292 * t387 + t295 * t294;
t113 = pkin(10) * t264 + t130;
t491 = t297 * t108 + t301 * t113;
t11 = -qJD(4) * t42 - t297 * t60 + t301 * t77;
t488 = t10 * t301 - t11 * t297;
t487 = t1 * t300 - t2 * t296;
t485 = Ifges(4,1) * t514 + Ifges(4,4) * t513;
t129 = -t298 * t205 + t302 * (t216 * t294 + t242 * t292);
t353 = t299 * t370;
t346 = t292 * t353;
t283 = qJD(2) * t290;
t219 = t283 - t314;
t78 = -t205 * t369 + t216 * t354 + t302 * t219 + t220 * t386 + t242 * t356 + t246 * t390;
t71 = pkin(10) * t346 + t78;
t165 = t295 * t357 + t307;
t166 = t295 * t356 + t308;
t168 = -t220 * t292 + t294 * t246;
t87 = pkin(3) * t165 - pkin(10) * t166 + t168;
t16 = -qJD(4) * t491 - t297 * t71 + t301 * t87;
t482 = -t110 * mrSges(4,3) + t510;
t481 = t111 * mrSges(4,3) + t505;
t5 = Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t92;
t479 = -mrSges(5,3) * t10 - t525 + t5 / 0.2e1 - t515;
t6 = Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t92;
t477 = t6 / 0.2e1;
t476 = Ifges(6,1) * t473 + Ifges(6,4) * t472 + Ifges(6,5) * t464;
t470 = -t44 / 0.2e1;
t401 = t187 * Ifges(5,3);
t405 = t158 * Ifges(5,5);
t406 = t157 * Ifges(5,6);
t83 = t401 + t405 + t406;
t469 = t83 / 0.2e1;
t155 = Ifges(4,6) * t160;
t156 = Ifges(4,5) * t159;
t94 = Ifges(4,3) * t323 - t155 + t156;
t463 = t94 / 0.2e1;
t462 = Ifges(4,5) * t323 / 0.2e1 + t485;
t461 = pkin(1) * mrSges(3,1);
t460 = pkin(1) * mrSges(3,2);
t459 = -t117 / 0.2e1;
t457 = -t118 / 0.2e1;
t395 = t238 * Ifges(4,6);
t400 = t189 * Ifges(4,2);
t422 = Ifges(4,4) * t190;
t127 = t395 + t400 + t422;
t455 = -t127 / 0.2e1;
t186 = Ifges(4,4) * t189;
t396 = t238 * Ifges(4,5);
t398 = t190 * Ifges(4,1);
t128 = t186 + t396 + t398;
t454 = t128 / 0.2e1;
t453 = -t153 / 0.2e1;
t450 = t157 / 0.2e1;
t448 = t158 / 0.2e1;
t445 = t187 / 0.2e1;
t442 = t190 / 0.2e1;
t440 = t295 / 0.2e1;
t437 = -t301 / 0.2e1;
t90 = Ifges(5,5) * t91;
t89 = Ifges(5,6) * t92;
t9 = -pkin(4) * t160 - t11;
t432 = t297 * t9;
t429 = t63 * mrSges(4,2);
t64 = (t204 * t294 + t239 * t292) * t302 + t349;
t428 = t64 * mrSges(4,1);
t427 = qJD(2) / 0.2e1;
t426 = mrSges(4,3) * t189;
t425 = mrSges(4,3) * t190;
t423 = Ifges(3,4) * t299;
t420 = Ifges(5,4) * t297;
t419 = Ifges(5,4) * t301;
t403 = t159 * Ifges(4,4);
t121 = mrSges(5,1) * t187 - mrSges(5,3) * t158;
t62 = -mrSges(6,1) * t117 + mrSges(6,2) * t118;
t392 = -t62 + t121;
t388 = t293 * t299;
t384 = t296 * t297;
t383 = t296 * t301;
t382 = t297 * t300;
t378 = t300 * t301;
t223 = -t266 * t296 - t300 * t389;
t147 = qJD(5) * t223 + t221 * t300 + t296 * t357;
t150 = t183 * t300 + t235 * t296;
t376 = t147 - t150;
t321 = -t266 * t300 + t296 * t389;
t148 = qJD(5) * t321 - t221 * t296 + t300 * t357;
t149 = -t183 * t296 + t235 * t300;
t375 = t148 - t149;
t31 = Ifges(5,3) * t160 - t89 + t90;
t259 = t372 * qJD(1);
t350 = t497 - (Ifges(3,2) * t303 + t423) * t371 / 0.2e1 - t259 * mrSges(3,3);
t348 = -t205 * t368 - t216 * t355 - t298 * t219 - t242 * t357;
t343 = t428 - t429;
t342 = -t11 * mrSges(5,1) + t10 * mrSges(5,2);
t340 = mrSges(6,1) * t300 - mrSges(6,2) * t296;
t337 = Ifges(6,1) * t296 + t416;
t335 = Ifges(6,2) * t300 + t417;
t333 = Ifges(6,5) * t296 + Ifges(6,6) * t300;
t52 = pkin(11) * t210 + t491;
t112 = -pkin(3) * t264 - t129;
t169 = t211 * t297 - t264 * t301;
t170 = t211 * t301 + t264 * t297;
t65 = pkin(4) * t169 - pkin(11) * t170 + t112;
t21 = t296 * t65 + t300 * t52;
t20 = -t296 * t52 + t300 * t65;
t331 = t297 * t42 + t301 * t41;
t53 = t108 * t301 - t113 * t297;
t73 = -t110 * t297 + t145 * t301;
t136 = t170 * t300 + t210 * t296;
t135 = -t170 * t296 + t210 * t300;
t176 = -t253 * t297 + t254 * t301;
t15 = t108 * t366 - t113 * t367 + t297 * t87 + t301 * t71;
t257 = -pkin(8) * t359 + t282;
t281 = Ifges(3,4) * t358;
t313 = t284 * Ifges(3,5) - t257 * mrSges(3,3) + Ifges(3,1) * t359 / 0.2e1 + t281 / 0.2e1;
t263 = t372 * qJD(2);
t309 = t110 * mrSges(4,1) - t111 * mrSges(4,2) + t190 * Ifges(4,5) + t189 * Ifges(4,6) + t498;
t72 = -t220 * t385 + (-pkin(3) * t353 - t246 * t302) * t292 - t348;
t306 = t411 / 0.2e1 - t410 / 0.2e1 - t414 / 0.2e1 - t409 / 0.2e1 + t421 / 0.2e1 + t521;
t277 = Ifges(3,5) * t303 * t351;
t268 = -pkin(8) * t388 + t290;
t262 = -pkin(8) * t353 + t283;
t256 = -t284 * mrSges(3,2) + mrSges(3,3) * t358;
t255 = mrSges(3,1) * t284 - mrSges(3,3) * t359;
t250 = qJD(1) * t263;
t249 = -pkin(8) * t345 + t279;
t248 = pkin(10) * t378 + t278 * t296;
t247 = -pkin(10) * t383 + t278 * t300;
t172 = pkin(4) * t389 - t176;
t162 = mrSges(4,1) * t238 - t425;
t161 = -mrSges(4,2) * t238 + t426;
t144 = -mrSges(4,1) * t189 + mrSges(4,2) * t190;
t142 = -mrSges(4,2) * t323 - mrSges(4,3) * t160;
t141 = mrSges(4,1) * t323 - mrSges(4,3) * t159;
t140 = t189 * t378 + t190 * t296;
t139 = -t189 * t383 + t190 * t300;
t133 = -t206 + (t218 * t294 + t245 * t292) * t302;
t120 = -mrSges(5,2) * t187 + mrSges(5,3) * t157;
t104 = -qJD(4) * t169 + t166 * t301 + t297 * t346;
t103 = qJD(4) * t170 + t166 * t297 - t301 * t346;
t102 = mrSges(4,1) * t160 + mrSges(4,2) * t159;
t101 = pkin(4) * t158 - pkin(11) * t157;
t100 = -mrSges(5,1) * t157 + mrSges(5,2) * t158;
t95 = -t160 * Ifges(4,2) + Ifges(4,6) * t323 + t403;
t82 = mrSges(6,1) * t153 - mrSges(6,3) * t118;
t81 = -mrSges(6,2) * t153 + mrSges(6,3) * t117;
t79 = (t220 * t294 + t246 * t292) * t302 + t348;
t67 = -mrSges(5,2) * t160 - mrSges(5,3) * t92;
t66 = mrSges(5,1) * t160 - mrSges(5,3) * t91;
t58 = -pkin(4) * t190 - t73;
t51 = -pkin(4) * t210 - t53;
t49 = qJD(5) * t135 + t104 * t300 + t165 * t296;
t48 = -qJD(5) * t136 - t104 * t296 + t165 * t300;
t40 = mrSges(5,1) * t92 + mrSges(5,2) * t91;
t28 = t101 * t296 + t300 * t41;
t27 = t101 * t300 - t296 * t41;
t24 = pkin(4) * t103 - pkin(11) * t104 + t72;
t23 = -mrSges(6,2) * t92 + mrSges(6,3) * t37;
t22 = mrSges(6,1) * t92 - mrSges(6,3) * t36;
t14 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t13 = -pkin(4) * t165 - t16;
t12 = pkin(11) * t165 + t15;
t4 = -qJD(5) * t21 - t12 * t296 + t24 * t300;
t3 = qJD(5) * t20 + t12 * t300 + t24 * t296;
t7 = [m(5) * (t10 * t491 + t11 * t53 + t112 * t61 + t15 * t42 + t16 * t41 + t72 * t98) + t491 * t67 + m(3) * (t249 * t372 - t250 * t268 - t257 * t263 + t259 * t262) + (Ifges(6,5) * t136 + Ifges(6,6) * t135) * t464 + (Ifges(6,5) * t49 + Ifges(6,6) * t48) * t452 + (Ifges(6,4) * t136 + Ifges(6,2) * t135) * t472 + (Ifges(6,4) * t49 + Ifges(6,2) * t48) * t458 + (Ifges(5,5) * t104 + Ifges(5,3) * t165) * t445 + (Ifges(5,5) * t170 + Ifges(5,3) * t210) * t447 + (Ifges(5,4) * t170 + Ifges(5,6) * t210) * t465 + (Ifges(5,4) * t104 + Ifges(5,6) * t165) * t450 + (Ifges(6,1) * t136 + Ifges(6,4) * t135) * t473 + (Ifges(6,1) * t49 + Ifges(6,4) * t48) * t456 + (Ifges(5,1) * t170 + Ifges(5,5) * t210) * t466 + (Ifges(5,1) * t104 + Ifges(5,5) * t165) * t448 + (t1 * t135 - t136 * t2 - t17 * t49 + t18 * t48) * mrSges(6,3) + (-t10 * t210 + t104 * t98 - t165 * t42 + t170 * t61) * mrSges(5,2) + ((Ifges(3,5) * t440 - t268 * mrSges(3,3) + (-0.2e1 * t460 + 0.3e1 / 0.2e1 * Ifges(3,4) * t303) * t293) * t303 + (-Ifges(3,6) * t295 + t292 * (Ifges(4,5) * t211 - Ifges(4,6) * t210 + Ifges(4,3) * t264) / 0.2e1 - t372 * mrSges(3,3) + (-0.2e1 * t461 - 0.3e1 / 0.2e1 * t423 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t303) * t293) * t299) * t351 + t211 * t462 + t264 * t463 + t104 * t467 + t165 * t469 + t170 * t474 + t136 * t476 + t135 * t477 + (-t110 * t166 - t111 * t165 - t210 * t63 - t211 * t64) * mrSges(4,3) + t166 * t454 + t165 * t455 + (-Ifges(5,4) * t448 + Ifges(6,5) * t456 - Ifges(5,2) * t450 - Ifges(5,6) * t445 + Ifges(6,6) * t458 + Ifges(6,3) * t452 - t521) * t103 + m(4) * (t110 * t79 + t111 * t78 + t129 * t64 + t130 * t63 + t151 * t168 + t163 * t164) + m(6) * (t1 * t21 + t13 * t38 + t17 * t4 + t18 * t3 + t2 * t20 + t51 * t9) + t479 * t169 + (t313 * t303 + (t497 + (t498 + t309) * t292 + t350) * t299) * t370 + t189 * (Ifges(4,4) * t166 - Ifges(4,2) * t165) / 0.2e1 - t264 * t429 + (Ifges(4,5) * t166 - Ifges(4,6) * t165) * t512 + (Ifges(4,4) * t211 - Ifges(4,2) * t210 + Ifges(4,6) * t264) * t513 + (Ifges(4,1) * t211 - Ifges(4,4) * t210 + Ifges(4,5) * t264) * t514 + t249 * (-t295 * mrSges(3,2) + mrSges(3,3) * t387) - t250 * (mrSges(3,1) * t295 - mrSges(3,3) * t388) + t277 * t440 + (Ifges(4,1) * t166 - Ifges(4,4) * t165) * t442 + t264 * t428 + t3 * t81 + t4 * t82 + t53 * t66 + t13 * t62 + t51 * t14 + t48 * t44 / 0.2e1 + t38 * (-mrSges(6,1) * t48 + mrSges(6,2) * t49) + t49 * t45 / 0.2e1 + t72 * t100 + t112 * t40 + t15 * t120 + t16 * t121 + t9 * (-mrSges(6,1) * t135 + mrSges(6,2) * t136) + t129 * t141 + t130 * t142 + t78 * t161 + t79 * t162 + t164 * t102 + t41 * (mrSges(5,1) * t165 - mrSges(5,3) * t104) + t151 * (mrSges(4,1) * t165 + mrSges(4,2) * t166) + t168 * t144 + t210 * t31 / 0.2e1 + t11 * (mrSges(5,1) * t210 - mrSges(5,3) * t170) - t210 * t95 / 0.2e1 + t163 * (mrSges(4,1) * t210 + mrSges(4,2) * t211) + t262 * t256 - t263 * t255 + t20 * t22 + t21 * t23; t483 * t523 - t480 * t524 + (Ifges(5,5) * t183 - Ifges(5,6) * t182) * t446 + (Ifges(5,4) * t183 - Ifges(5,2) * t182) * t451 + (Ifges(5,1) * t183 - Ifges(5,4) * t182) * t449 + (Ifges(4,4) * t444 - t128 / 0.2e1 - t482 + t503) * t236 + (-t11 * mrSges(5,3) + t501 + t511) * t266 + (t10 * t489 + t11 * t176 + t252 * t61 + t41 * t493 + t42 * t492 + t490 * t98) * m(5) + t489 * t67 + (-Ifges(6,5) * t321 + Ifges(6,6) * t223) * t464 + (-Ifges(6,4) * t321 + Ifges(6,2) * t223) * t472 + (-Ifges(6,1) * t321 + Ifges(6,4) * t223) * t473 + (t1 * t223 - t17 * t376 + t18 * t375 + t2 * t321) * mrSges(6,3) + t9 * (-mrSges(6,1) * t223 - mrSges(6,2) * t321) - t321 * t476 + (t84 - t43) * (t182 / 0.2e1 - t222 / 0.2e1) + (-Ifges(4,4) * t443 - t83 / 0.2e1 + t127 / 0.2e1 + t481 + t500) * t235 + t223 * t477 + (t147 / 0.2e1 - t150 / 0.2e1) * t45 + ((-t281 / 0.2e1 + t371 * t460 - t313) * t303 + ((t461 + t423 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t303) * t371 + (-qJD(2) + t284 / 0.2e1) * Ifges(3,6) + ((Ifges(4,5) * t298 + Ifges(4,6) * t302) * t292 * t427 + (t294 * t427 + t441) * Ifges(4,3) - t309) * t292 - t350) * t299) * t371 + (t463 + t156 / 0.2e1 - t155 / 0.2e1 + t343) * t294 + t492 * t120 + t493 * t121 + t494 * t62 + (-t133 - t261) * t162 + t490 * t100 - m(4) * (t110 * t133 + t111 * t134 + t151 * t167) + (Ifges(5,4) * t221 - Ifges(5,2) * t222) * t450 + (Ifges(6,5) * t147 + Ifges(6,6) * t148 + Ifges(6,3) * t222) * t452 + (Ifges(6,5) * t150 + Ifges(6,6) * t149 + Ifges(6,3) * t182) * t453 + (Ifges(6,1) * t147 + Ifges(6,4) * t148 + Ifges(6,5) * t222) * t456 + (Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t182) * t457 + (Ifges(6,4) * t147 + Ifges(6,2) * t148 + Ifges(6,6) * t222) * t458 + (Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t182) * t459 + ((-m(4) * t163 - t102) * pkin(2) + (t163 * mrSges(4,2) - t64 * mrSges(4,3) + t462 + t485) * t298 + (-t163 * mrSges(4,1) + t403 / 0.2e1 + t95 / 0.2e1 - t31 / 0.2e1 - t90 / 0.2e1 + t89 / 0.2e1 + t63 * mrSges(4,3) + (-Ifges(4,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t160 + t342) * t302 + ((t186 / 0.2e1 + t398 / 0.2e1 + t396 / 0.2e1 + t454 + t482) * t302 + (-t400 / 0.2e1 - t422 / 0.2e1 - t395 / 0.2e1 + t455 + t469 + t406 / 0.2e1 + t405 / 0.2e1 + t401 / 0.2e1 - t481) * t298) * qJD(3)) * t292 + m(4) * (-t110 * t261 + t111 * t260 + t267 * t64 + t269 * t63) + t479 * t265 + t277 + t528 * t161 + (-mrSges(6,1) * t375 + mrSges(6,2) * t376) * t38 + (Ifges(5,5) * t221 - Ifges(5,6) * t222) * t445 + (Ifges(5,1) * t221 - Ifges(5,4) * t222) * t448 + t115 * t22 + t116 * t23 - t167 * t144 + t172 * t14 + t176 * t66 + t495 * t82 + t496 * t81 + (t1 * t116 + t115 * t2 + t17 * t495 + t172 * t9 + t18 * t496 + t38 * t494) * m(6) - t249 * mrSges(3,2) - t250 * mrSges(3,1) + t252 * t40 - t257 * t256 + t259 * t255 + t267 * t141 + t269 * t142 + (t148 / 0.2e1 - t149 / 0.2e1) * t44 + (-t183 / 0.2e1 + t221 / 0.2e1) * t85; ((t424 / 0.2e1 + (m(6) * t38 - t392) * pkin(10) - t507) * qJD(4) + t67 * pkin(10) + t515) * t301 + (t128 + t186) * t444 + t508 * t81 + (pkin(10) * t432 + t1 * t248 + t17 * t509 + t18 * t508 + t2 * t247 - t38 * t58) * m(6) + t509 * t82 + (t500 + t505) * t190 + (t83 - t422) * t443 + t488 * mrSges(5,3) + (t334 * t464 + t336 * t472 + t338 * t473 + (-t413 / 0.2e1 - t306 - t120 * pkin(10)) * qJD(4) + (t14 - t66) * pkin(10) + (t45 * t439 + t300 * t470 + t38 * t340 + t335 * t459 + t337 * t457 + t333 * t453 + (t17 * t296 - t18 * t300) * mrSges(6,3)) * qJD(5) + t501) * t297 + t139 * t470 + t382 * t476 + (Ifges(6,1) * t140 + Ifges(6,4) * t139) * t457 + (Ifges(6,5) * t140 + Ifges(6,6) * t139) * t453 + t94 + ((-Ifges(5,2) * t297 + t419) * t451 + t331 * mrSges(5,3) + (Ifges(5,5) * t301 - Ifges(5,6) * t297) * t446 + (Ifges(5,1) * t301 - t420) * t449 - t98 * (mrSges(5,1) * t297 + mrSges(5,2) * t301) + t85 * t437 - t510 + t503 + (Ifges(6,5) * t457 + Ifges(6,6) * t459 + Ifges(6,3) * t453 - t504 + t506) * t297) * t189 + t343 + (Ifges(6,4) * t140 + Ifges(6,2) * t139) * t459 + t419 * t466 + t420 * t465 + (-t1 * t384 - t18 * t139 + t17 * t140 - t2 * t382) * mrSges(6,3) + (-t100 + t162 + t425) * t111 + (-t161 + t426) * t110 - t6 * t384 / 0.2e1 + t127 * t442 + t339 * t432 + t5 * t437 - t58 * t62 - pkin(3) * t40 - t74 * t120 - t73 * t121 - t140 * t45 / 0.2e1 - t38 * (-mrSges(6,1) * t139 + mrSges(6,2) * t140) + t247 * t22 + t248 * t23 + (-pkin(3) * t61 - t111 * t98 - t41 * t73 - t42 * t74 + (-qJD(4) * t331 + t488) * pkin(10)) * m(5); -t342 + t31 + t306 * t158 + t392 * t42 + t487 * mrSges(6,3) - t28 * t81 - t27 * t82 + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t158 + t507) * t157 + t516 * qJD(5) - t41 * t120 + t333 * t464 + t335 * t472 + t337 * t473 - t9 * t340 + t296 * t476 + t6 * t438 - pkin(4) * t14 + (-pkin(4) * t9 - t17 * t27 - t18 * t28 - t38 * t42) * m(6) + (-t22 * t296 + t23 * t300 + m(6) * t487 + (-m(6) * t332 - t296 * t81 - t300 * t82) * qJD(5)) * pkin(11); -t38 * (mrSges(6,1) * t118 + mrSges(6,2) * t117) + (Ifges(6,1) * t117 - t418) * t457 + t44 * t456 + (Ifges(6,5) * t117 - Ifges(6,6) * t118) * t453 - t17 * t81 + t18 * t82 + (t117 * t17 + t118 * t18) * mrSges(6,3) + t5 + (-Ifges(6,2) * t118 + t114 + t45) * t459 + t484;];
tauc = t7(:);
