% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:29:50
% EndTime: 2018-11-23 17:30:08
% DurationCPUTime: 18.75s
% Computational Cost: add. (17364->868), mult. (43517->1195), div. (0->0), fcn. (32472->10), ass. (0->389)
t316 = sin(qJ(4));
t319 = cos(qJ(5));
t321 = cos(qJ(2));
t312 = sin(pkin(6));
t391 = qJD(1) * t312;
t315 = sin(qJ(5));
t317 = sin(qJ(2));
t401 = t315 * t317;
t233 = (-t316 * t401 + t319 * t321) * t391;
t320 = cos(qJ(4));
t382 = qJD(5) * t319;
t386 = qJD(4) * t316;
t512 = -t315 * t386 + t320 * t382 + t233;
t355 = pkin(4) * t320 + pkin(10) * t316;
t270 = qJD(4) * t355 + qJD(3);
t279 = pkin(4) * t316 - pkin(10) * t320 + qJ(3);
t322 = -pkin(2) - pkin(9);
t384 = qJD(4) * t322;
t367 = t320 * t384;
t400 = t316 * t322;
t375 = t315 * t400;
t154 = -qJD(5) * t375 + t315 * t270 + t279 * t382 + t319 * t367;
t372 = t317 * t391;
t296 = pkin(2) * t372;
t346 = pkin(9) * t317 - qJ(3) * t321;
t221 = t346 * t391 + t296;
t371 = t321 * t391;
t313 = cos(pkin(6));
t390 = qJD(1) * t313;
t379 = pkin(1) * t390;
t254 = pkin(8) * t371 + t317 * t379;
t223 = pkin(3) * t371 + t254;
t140 = t320 * t221 + t316 * t223;
t123 = pkin(10) * t371 + t140;
t301 = t321 * t379;
t464 = pkin(3) + pkin(8);
t160 = t301 + (-t355 - t464) * t372;
t77 = t319 * t123 + t315 * t160;
t511 = t154 - t77;
t399 = t317 * t319;
t234 = (t315 * t321 + t316 * t399) * t391;
t262 = t319 * t270;
t303 = t319 * t400;
t360 = t320 * t372;
t364 = -t315 * t322 + pkin(5);
t435 = pkin(11) * t320;
t436 = pkin(11) * t319;
t76 = -t123 * t315 + t319 * t160;
t510 = pkin(5) * t360 + pkin(11) * t234 + t262 + (-t303 + (-t279 + t435) * t315) * qJD(5) + (t316 * t436 + t320 * t364) * qJD(4) - t76;
t509 = -t512 * pkin(11) + t511;
t463 = -pkin(11) - pkin(10);
t373 = qJD(5) * t463;
t304 = qJD(2) + t390;
t240 = -t316 * t304 - t320 * t371;
t404 = t240 * t315;
t165 = t304 * t322 + t372 * t464 + qJD(3) - t301;
t365 = -qJ(3) * t317 - pkin(1);
t215 = (t321 * t322 + t365) * t312;
t193 = qJD(1) * t215;
t103 = t165 * t320 - t316 * t193;
t359 = t316 * t371;
t241 = t304 * t320 - t359;
t164 = pkin(4) * t241 - pkin(10) * t240;
t70 = t319 * t103 + t315 * t164;
t508 = pkin(11) * t404 + t315 * t373 - t70;
t69 = -t103 * t315 + t319 * t164;
t507 = -pkin(5) * t241 + t240 * t436 + t319 * t373 - t69;
t294 = t304 * qJ(3);
t186 = t294 + t223;
t110 = -pkin(4) * t240 - pkin(10) * t241 + t186;
t104 = t165 * t316 + t193 * t320;
t289 = qJD(4) + t372;
t89 = pkin(10) * t289 + t104;
t52 = t319 * t110 - t315 * t89;
t53 = t110 * t315 + t319 * t89;
t344 = t315 * t53 + t319 * t52;
t353 = mrSges(6,1) * t315 + mrSges(6,2) * t319;
t439 = t319 / 0.2e1;
t440 = -t315 / 0.2e1;
t172 = -t241 * t315 + t289 * t319;
t237 = qJD(5) - t240;
t173 = t241 * t319 + t289 * t315;
t433 = Ifges(6,4) * t173;
t85 = Ifges(6,2) * t172 + Ifges(6,6) * t237 + t433;
t169 = Ifges(6,4) * t172;
t86 = Ifges(6,1) * t173 + Ifges(6,5) * t237 + t169;
t88 = -pkin(4) * t289 - t103;
t506 = t344 * mrSges(6,3) - t88 * t353 - t439 * t86 - t440 * t85;
t314 = sin(qJ(6));
t318 = cos(qJ(6));
t363 = t318 * t172 - t173 * t314;
t102 = Ifges(7,4) * t363;
t109 = t172 * t314 + t173 * t318;
t39 = -pkin(11) * t173 + t52;
t35 = pkin(5) * t237 + t39;
t40 = pkin(11) * t172 + t53;
t408 = t314 * t40;
t15 = t318 * t35 - t408;
t407 = t318 * t40;
t16 = t314 * t35 + t407;
t430 = Ifges(7,4) * t109;
t228 = qJD(6) + t237;
t448 = -t228 / 0.2e1;
t456 = -t109 / 0.2e1;
t458 = -t363 / 0.2e1;
t48 = Ifges(7,1) * t109 + Ifges(7,5) * t228 + t102;
t71 = -pkin(5) * t172 + t88;
t505 = (Ifges(7,5) * t363 - Ifges(7,6) * t109) * t448 + (t109 * t16 + t15 * t363) * mrSges(7,3) + (-Ifges(7,2) * t109 + t102 + t48) * t458 - t71 * (mrSges(7,1) * t109 + mrSges(7,2) * t363) + (Ifges(7,1) * t363 - t430) * t456;
t493 = -t304 / 0.2e1;
t504 = mrSges(4,1) + mrSges(3,3);
t503 = mrSges(4,2) - mrSges(3,1);
t383 = qJD(5) * t315;
t366 = qJD(2) * t391;
t357 = t321 * t366;
t358 = t317 * t366;
t287 = pkin(2) * t358;
t387 = qJD(3) * t317;
t329 = (qJD(2) * t346 - t387) * t312;
t170 = qJD(1) * t329 + t287;
t308 = t313 * t317 * pkin(1);
t402 = t312 * t321;
t224 = (t402 * t464 + t308) * qJD(2);
t197 = qJD(1) * t224;
t385 = qJD(4) * t320;
t59 = t165 * t385 + t320 * t170 - t193 * t386 + t316 * t197;
t55 = pkin(10) * t357 + t59;
t288 = qJD(2) * t301;
t291 = t304 * qJD(3);
t403 = t312 * t317;
t362 = t464 * t403;
t341 = qJD(2) * t362;
t171 = -qJD(1) * t341 + t288 + t291;
t389 = qJD(2) * t317;
t368 = t316 * t389;
t181 = -t304 * t386 + (-t321 * t385 + t368) * t391;
t182 = -qJD(4) * t359 + t304 * t385 - t320 * t358;
t80 = pkin(4) * t182 - pkin(10) * t181 + t171;
t13 = t110 * t382 + t315 * t80 + t319 * t55 - t383 * t89;
t97 = -qJD(5) * t173 - t181 * t315 + t319 * t357;
t10 = pkin(11) * t97 + t13;
t14 = -qJD(5) * t53 - t315 * t55 + t319 * t80;
t96 = qJD(5) * t172 + t181 * t319 + t315 * t357;
t6 = pkin(5) * t182 - pkin(11) * t96 + t14;
t2 = qJD(6) * t15 + t10 * t318 + t314 * t6;
t3 = -qJD(6) * t16 - t10 * t314 + t318 * t6;
t502 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t269 = t319 * t279;
t192 = t316 * t364 - t319 * t435 + t269;
t239 = t315 * t279 + t303;
t208 = -t315 * t435 + t239;
t120 = t192 * t318 - t208 * t314;
t501 = qJD(6) * t120 + t510 * t314 + t509 * t318;
t121 = t192 * t314 + t208 * t318;
t500 = -qJD(6) * t121 - t509 * t314 + t510 * t318;
t139 = -t316 * t221 + t223 * t320;
t122 = -pkin(4) * t371 - t139;
t499 = t512 * pkin(5) + t316 * t384 - t122;
t338 = t314 * t315 - t318 * t319;
t479 = qJD(5) + qJD(6);
t498 = t479 * t338;
t348 = Ifges(6,5) * t319 - Ifges(6,6) * t315;
t431 = Ifges(6,4) * t319;
t350 = -Ifges(6,2) * t315 + t431;
t432 = Ifges(6,4) * t315;
t352 = Ifges(6,1) * t319 - t432;
t444 = t237 / 0.2e1;
t451 = t173 / 0.2e1;
t453 = t172 / 0.2e1;
t497 = -t348 * t444 - t350 * t453 - t352 * t451 + t506;
t253 = pkin(8) * t372 - t301;
t496 = -qJD(3) - t253;
t32 = qJD(6) * t363 + t314 * t97 + t318 * t96;
t470 = t32 / 0.2e1;
t33 = -qJD(6) * t109 - t314 * t96 + t318 * t97;
t469 = t33 / 0.2e1;
t47 = Ifges(7,2) * t363 + Ifges(7,6) * t228 + t430;
t494 = t47 / 0.2e1;
t450 = t182 / 0.2e1;
t492 = -t391 / 0.2e1;
t93 = Ifges(6,6) * t97;
t94 = Ifges(6,5) * t96;
t36 = Ifges(6,3) * t182 + t93 + t94;
t30 = Ifges(7,6) * t33;
t31 = Ifges(7,5) * t32;
t7 = Ifges(7,3) * t182 + t30 + t31;
t491 = t36 + t7;
t490 = t14 * mrSges(6,1) - t13 * mrSges(6,2) + t502;
t292 = t463 * t315;
t293 = t463 * t319;
t226 = t292 * t318 + t293 * t314;
t485 = qJD(6) * t226 + t507 * t314 + t318 * t508;
t227 = t292 * t314 - t293 * t318;
t484 = -qJD(6) * t227 - t314 * t508 + t507 * t318;
t57 = -mrSges(7,1) * t363 + mrSges(7,2) * t109;
t483 = m(7) * t71 + t57;
t305 = pkin(8) * t403;
t438 = pkin(1) * t321;
t374 = -pkin(2) - t438;
t196 = pkin(3) * t403 + t305 + (-pkin(9) + t374) * t313;
t127 = t316 * t196 + t320 * t215;
t118 = pkin(10) * t403 + t127;
t266 = pkin(8) * t402 + t308;
t244 = -t313 * qJ(3) - t266;
t214 = pkin(3) * t402 - t244;
t263 = t313 * t316 + t320 * t402;
t376 = t316 * t402;
t264 = t313 * t320 - t376;
t134 = pkin(4) * t263 - pkin(10) * t264 + t214;
t68 = t319 * t118 + t315 * t134;
t231 = t304 * t319 + t315 * t360;
t232 = t304 * t315 - t319 * t360;
t273 = t314 * t319 + t315 * t318;
t258 = t273 * t320;
t482 = -qJD(4) * t258 - t231 * t318 + t232 * t314 + t316 * t498;
t201 = t479 * t273;
t260 = t338 * t320;
t481 = -qJD(4) * t260 - t201 * t316 - t231 * t314 - t232 * t318;
t235 = Ifges(5,4) * t240;
t413 = t289 * Ifges(5,5);
t415 = t241 * Ifges(5,1);
t145 = t235 + t413 + t415;
t480 = t103 * mrSges(5,3) - t145 / 0.2e1 - t186 * mrSges(5,2) - t413 / 0.2e1 - t235 / 0.2e1;
t345 = t13 * t319 - t14 * t315;
t60 = -t165 * t386 - t316 * t170 - t193 * t385 + t197 * t320;
t478 = t60 * mrSges(5,1) - t59 * mrSges(5,2) + Ifges(5,5) * t181 - Ifges(5,6) * t182;
t477 = -Ifges(3,4) * t371 / 0.2e1 + Ifges(3,5) * t493;
t475 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t450;
t474 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t450;
t473 = Ifges(3,5) / 0.2e1;
t472 = Ifges(4,5) / 0.2e1;
t471 = Ifges(5,2) / 0.2e1;
t38 = t96 * Ifges(6,1) + t97 * Ifges(6,4) + t182 * Ifges(6,5);
t468 = t38 / 0.2e1;
t467 = -t85 / 0.2e1;
t466 = t96 / 0.2e1;
t465 = t97 / 0.2e1;
t461 = pkin(1) * mrSges(3,1);
t460 = pkin(1) * mrSges(3,2);
t457 = t363 / 0.2e1;
t455 = t109 / 0.2e1;
t454 = -t172 / 0.2e1;
t452 = -t173 / 0.2e1;
t447 = t228 / 0.2e1;
t445 = -t237 / 0.2e1;
t443 = -t263 / 0.2e1;
t441 = t264 / 0.2e1;
t437 = pkin(5) * t315;
t434 = Ifges(5,4) * t241;
t429 = Ifges(4,6) * t321;
t428 = t363 * Ifges(7,6);
t427 = t109 * Ifges(7,5);
t424 = t172 * Ifges(6,6);
t423 = t173 * Ifges(6,5);
t422 = t181 * Ifges(5,1);
t421 = t181 * Ifges(5,4);
t420 = t182 * Ifges(5,4);
t419 = t228 * Ifges(7,3);
t418 = t237 * Ifges(6,3);
t417 = t240 * Ifges(5,2);
t412 = t289 * Ifges(5,6);
t409 = t304 * Ifges(4,5);
t161 = mrSges(5,1) * t357 - mrSges(5,3) * t181;
t50 = -mrSges(6,1) * t97 + mrSges(6,2) * t96;
t406 = t161 - t50;
t136 = -t201 * t320 + t338 * t386;
t151 = t233 * t314 + t234 * t318;
t398 = t136 - t151;
t138 = t273 * t386 + t320 * t498;
t150 = t233 * t318 - t234 * t314;
t397 = t138 - t150;
t146 = t273 * t240;
t396 = -t146 + t201;
t147 = t338 * t240;
t395 = -t147 + t498;
t111 = -mrSges(6,1) * t172 + mrSges(6,2) * t173;
t188 = mrSges(5,1) * t289 - mrSges(5,3) * t241;
t394 = t188 - t111;
t163 = -mrSges(5,1) * t240 + mrSges(5,2) * t241;
t249 = -mrSges(4,1) * t371 - mrSges(4,3) * t304;
t393 = -t249 + t163;
t392 = -t304 * t503 - t372 * t504;
t388 = qJD(2) * t321;
t381 = t313 * t438;
t380 = t304 / 0.2e1 - qJD(2);
t378 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t370 = t312 * t389;
t369 = t312 * t388;
t67 = -t118 * t315 + t319 * t134;
t126 = t196 * t320 - t316 * t215;
t361 = t316 * t372;
t354 = mrSges(6,1) * t319 - mrSges(6,2) * t315;
t351 = Ifges(6,1) * t315 + t431;
t349 = Ifges(6,2) * t319 + t432;
t347 = Ifges(6,5) * t315 + Ifges(6,6) * t319;
t212 = t264 * t319 + t312 * t401;
t49 = pkin(5) * t263 - pkin(11) * t212 + t67;
t211 = -t264 * t315 + t312 * t399;
t58 = pkin(11) * t211 + t68;
t22 = -t314 * t58 + t318 * t49;
t23 = t314 * t49 + t318 * t58;
t343 = t315 * t52 - t319 * t53;
t74 = mrSges(6,1) * t182 - mrSges(6,3) * t96;
t75 = -mrSges(6,2) * t182 + mrSges(6,3) * t97;
t342 = -t315 * t74 + t319 * t75;
t340 = t103 * t316 - t104 * t320;
t128 = -mrSges(6,2) * t237 + mrSges(6,3) * t172;
t129 = mrSges(6,1) * t237 - mrSges(6,3) * t173;
t339 = -t315 * t128 - t319 * t129;
t132 = t211 * t318 - t212 * t314;
t133 = t211 * t314 + t212 * t318;
t302 = qJD(2) * t381;
t255 = -pkin(8) * t370 + t302;
t298 = pkin(2) * t370;
t190 = t298 + t329;
t73 = -t316 * t190 - t196 * t386 - t215 * t385 + t224 * t320;
t72 = t320 * t190 + t196 * t385 - t215 * t386 + t316 * t224;
t65 = pkin(10) * t369 + t72;
t310 = t313 * qJD(3);
t195 = t302 + t310 - t341;
t209 = -qJD(4) * t263 + t312 * t368;
t210 = -qJD(4) * t376 + t313 * t385 - t320 * t370;
t95 = pkin(4) * t210 - pkin(10) * t209 + t195;
t20 = -t118 * t383 + t134 * t382 + t315 * t95 + t319 * t65;
t242 = -pkin(8) * t358 + t288;
t117 = -pkin(4) * t403 - t126;
t245 = (-pkin(2) * t321 + t365) * t312;
t333 = (-qJ(3) * t388 - t387) * t312;
t230 = qJD(1) * t245;
t332 = Ifges(3,6) * t493 + (Ifges(3,4) * t317 + Ifges(3,2) * t321) * t492 + t409 / 0.2e1 + (-Ifges(4,6) * t317 - Ifges(4,3) * t321) * t391 / 0.2e1 - t230 * mrSges(4,2) - t254 * mrSges(3,3);
t256 = t266 * qJD(2);
t205 = -t242 - t291;
t243 = qJD(1) * t256;
t331 = -t242 * mrSges(3,2) - t205 * mrSges(4,3) + t243 * t503;
t66 = -pkin(4) * t369 - t73;
t21 = -qJD(5) * t68 - t315 * t65 + t319 * t95;
t56 = -pkin(4) * t357 - t60;
t327 = -t415 / 0.2e1 + t480;
t326 = t103 * mrSges(5,1) + t253 * mrSges(3,3) + t289 * Ifges(5,3) + t241 * Ifges(5,5) + t240 * Ifges(5,6) + Ifges(3,1) * t372 / 0.2e1 + Ifges(4,4) * t493 + (-t317 * Ifges(4,2) - t429) * t492 - t104 * mrSges(5,2) - t230 * mrSges(4,3) - t477;
t144 = t412 + t417 + t434;
t46 = t419 + t427 + t428;
t84 = t418 + t423 + t424;
t324 = t412 / 0.2e1 - t418 / 0.2e1 - t423 / 0.2e1 - t424 / 0.2e1 - t419 / 0.2e1 - t427 / 0.2e1 - t428 / 0.2e1 - t84 / 0.2e1 - t46 / 0.2e1 + t144 / 0.2e1 - t15 * mrSges(7,1) + t16 * mrSges(7,2) - t52 * mrSges(6,1) + t53 * mrSges(6,2) + t104 * mrSges(5,3) - t186 * mrSges(5,1) + t434 / 0.2e1;
t323 = (-t417 / 0.2e1 - t324) * t320;
t309 = -pkin(5) * t319 - pkin(4);
t286 = Ifges(3,5) * t357;
t285 = Ifges(4,5) * t358;
t284 = Ifges(5,3) * t357;
t271 = (-t322 + t437) * t320;
t265 = -t305 + t381;
t259 = t338 * t316;
t257 = t273 * t316;
t252 = -qJ(3) * t371 + t296;
t251 = (mrSges(4,2) * t321 - mrSges(4,3) * t317) * t391;
t248 = -mrSges(3,2) * t304 + mrSges(3,3) * t371;
t246 = t313 * t374 + t305;
t238 = t269 - t375;
t236 = -t255 - t310;
t225 = t298 + t333;
t222 = -qJD(1) * t362 + t301;
t220 = -t294 - t254;
t213 = -pkin(2) * t304 - t496;
t202 = qJD(1) * t333 + t287;
t187 = -mrSges(5,2) * t289 + mrSges(5,3) * t240;
t162 = -mrSges(5,2) * t357 - mrSges(5,3) * t182;
t155 = -qJD(5) * t239 - t315 * t367 + t262;
t115 = -qJD(5) * t212 - t209 * t315 + t319 * t369;
t114 = qJD(5) * t211 + t209 * t319 + t315 * t369;
t112 = mrSges(5,1) * t182 + mrSges(5,2) * t181;
t100 = Ifges(5,5) * t357 - t420 + t422;
t99 = -t182 * Ifges(5,2) + Ifges(5,6) * t357 + t421;
t87 = -pkin(5) * t211 + t117;
t83 = pkin(5) * t404 + t104;
t82 = mrSges(7,1) * t228 - mrSges(7,3) * t109;
t81 = -mrSges(7,2) * t228 + mrSges(7,3) * t363;
t45 = -pkin(5) * t115 + t66;
t42 = -qJD(6) * t133 - t114 * t314 + t115 * t318;
t41 = qJD(6) * t132 + t114 * t318 + t115 * t314;
t37 = t96 * Ifges(6,4) + t97 * Ifges(6,2) + t182 * Ifges(6,6);
t34 = -pkin(5) * t97 + t56;
t29 = -mrSges(7,2) * t182 + mrSges(7,3) * t33;
t28 = mrSges(7,1) * t182 - mrSges(7,3) * t32;
t19 = t318 * t39 - t408;
t18 = -t314 * t39 - t407;
t17 = pkin(11) * t115 + t20;
t12 = pkin(5) * t210 - pkin(11) * t114 + t21;
t11 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t5 = -qJD(6) * t23 + t12 * t318 - t17 * t314;
t4 = qJD(6) * t22 + t12 * t314 + t17 * t318;
t1 = [(Ifges(6,5) * t212 + Ifges(7,5) * t133 + Ifges(6,6) * t211 + Ifges(7,6) * t132 + (Ifges(6,3) + Ifges(7,3)) * t263) * t450 + (t84 + t46) * t210 / 0.2e1 + (-t103 * t209 - t104 * t210 - t263 * t59 - t264 * t60) * mrSges(5,3) + t241 * (Ifges(5,1) * t209 - Ifges(5,4) * t210) / 0.2e1 + t181 * (Ifges(5,1) * t264 - Ifges(5,4) * t263) / 0.2e1 + t289 * (Ifges(5,5) * t209 - Ifges(5,6) * t210) / 0.2e1 + m(5) * (t103 * t73 + t104 * t72 + t126 * t60 + t127 * t59 + t171 * t214 + t186 * t195) + m(4) * (t202 * t245 + t205 * t244 + t213 * t256 + t220 * t236 + t225 * t230 + t243 * t246) + m(3) * (t242 * t266 - t243 * t265 + t253 * t256 + t254 * t255) + m(7) * (t15 * t5 + t16 * t4 + t2 * t23 + t22 * t3 + t34 * t87 + t45 * t71) + m(6) * (t117 * t56 + t13 * t68 + t14 * t67 + t20 * t53 + t21 * t52 + t66 * t88) + ((-t245 * mrSges(4,3) + Ifges(5,5) * t441 + Ifges(5,6) * t443 + t246 * mrSges(4,1) - t265 * mrSges(3,3) + (-Ifges(4,4) + t473) * t313 + (t321 * t378 - 0.2e1 * t460) * t312) * t321 + (t244 * mrSges(4,1) - t245 * mrSges(4,2) - t266 * mrSges(3,3) + (-Ifges(3,6) + t472) * t313 + (-t317 * t378 - 0.2e1 * t461) * t312 + (0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + Ifges(5,3) / 0.2e1) * t402) * t317) * t366 + t171 * (mrSges(5,1) * t263 + mrSges(5,2) * t264) + t14 * (mrSges(6,1) * t263 - mrSges(6,3) * t212) + t13 * (-mrSges(6,2) * t263 + mrSges(6,3) * t211) + t2 * (-mrSges(7,2) * t263 + mrSges(7,3) * t132) + t3 * (mrSges(7,1) * t263 - mrSges(7,3) * t133) + t255 * t248 + t236 * t249 + t225 * t251 + t56 * (-mrSges(6,1) * t211 + mrSges(6,2) * t212) + t214 * t112 + t52 * (mrSges(6,1) * t210 - mrSges(6,3) * t114) + t53 * (-mrSges(6,2) * t210 + mrSges(6,3) * t115) + t15 * (mrSges(7,1) * t210 - mrSges(7,3) * t41) + t16 * (-mrSges(7,2) * t210 + mrSges(7,3) * t42) + t186 * (mrSges(5,1) * t210 + mrSges(5,2) * t209) - t210 * t144 / 0.2e1 + t211 * t37 / 0.2e1 + t209 * t145 / 0.2e1 + t72 * t187 + t73 * t188 + t195 * t163 - t392 * t256 + t126 * t161 + t127 * t162 + t20 * t128 + t21 * t129 + t34 * (-mrSges(7,1) * t132 + mrSges(7,2) * t133) + t114 * t86 / 0.2e1 + t88 * (-mrSges(6,1) * t115 + mrSges(6,2) * t114) + t115 * t85 / 0.2e1 + t117 * t50 + t66 * t111 + t4 * t81 + t5 * t82 + t87 * t11 + t71 * (-mrSges(7,1) * t42 + mrSges(7,2) * t41) + t67 * t74 + t68 * t75 + t45 * t57 + t41 * t48 / 0.2e1 + t22 * t28 + t23 * t29 + t240 * (Ifges(5,4) * t209 - Ifges(5,2) * t210) / 0.2e1 - t182 * (Ifges(5,4) * t264 - Ifges(5,2) * t263) / 0.2e1 + (t285 / 0.2e1 + t286 / 0.2e1 + t331) * t313 + t491 * t263 / 0.2e1 + t100 * t441 + t99 * t443 + (Ifges(6,5) * t114 + Ifges(6,6) * t115 + Ifges(6,3) * t210) * t444 + (Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t210) * t447 + (Ifges(6,1) * t114 + Ifges(6,4) * t115 + Ifges(6,5) * t210) * t451 + ((-mrSges(4,1) * t205 + mrSges(4,2) * t202 + mrSges(3,3) * t242) * t321 + (t284 / 0.2e1 - t202 * mrSges(4,3) + t504 * t243 + t478) * t317 + ((t220 * mrSges(4,1) + (-Ifges(3,6) / 0.2e1 + t472) * t304 + t332) * t317 + (t213 * mrSges(4,1) + (t473 - Ifges(4,4) / 0.2e1) * t304 + t326) * t321) * qJD(2)) * t312 + (Ifges(6,4) * t114 + Ifges(6,2) * t115 + Ifges(6,6) * t210) * t453 + (Ifges(7,1) * t41 + Ifges(7,4) * t42 + Ifges(7,5) * t210) * t455 + (Ifges(7,4) * t41 + Ifges(7,2) * t42 + Ifges(7,6) * t210) * t457 + (Ifges(6,4) * t212 + Ifges(6,2) * t211 + Ifges(6,6) * t263) * t465 + (Ifges(6,1) * t212 + Ifges(6,4) * t211 + Ifges(6,5) * t263) * t466 + t212 * t468 + (Ifges(7,4) * t133 + Ifges(7,2) * t132 + Ifges(7,6) * t263) * t469 + (Ifges(7,1) * t133 + Ifges(7,4) * t132 + Ifges(7,5) * t263) * t470 + t133 * t474 + t132 * t475 + t42 * t494; (t323 + (-m(5) * t340 + t320 * t187) * t322 + (t348 * t445 + t350 * t454 + t352 * t452 + t327 + (m(6) * t88 - t394) * t322 + t506) * t316) * qJD(4) + t34 * (mrSges(7,1) * t258 - mrSges(7,2) * t260) + (-t15 * t398 + t16 * t397 - t2 * t258 + t260 * t3) * mrSges(7,3) + (-Ifges(7,5) * t260 - Ifges(7,6) * t258) * t450 + (-Ifges(7,4) * t260 - Ifges(7,2) * t258) * t469 + (-Ifges(7,1) * t260 - Ifges(7,4) * t258) * t470 + t511 * t128 + t331 + m(6) * (t13 * t239 + t14 * t238 + t53 * t154 + t52 * t155) + (t155 - t76) * t129 + (-t249 + t248) * t253 + ((-t326 + (-pkin(2) * qJD(2) - t213) * mrSges(4,1) + qJD(2) * (Ifges(5,5) * t320 - Ifges(5,6) * t316) / 0.2e1 + t380 * Ifges(4,4) + (t460 - t429 / 0.2e1) * t391 + t477) * t321 + (-t409 / 0.2e1 + (((Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t317 + t461) * t312 + (Ifges(3,2) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1) * t402) * qJD(1) + t380 * Ifges(3,6) + (-qJ(3) * qJD(2) - t220) * mrSges(4,1) + t327 * t316 + t323 - t332) * t317) * t391 + t285 + t286 + (t38 * t439 + t37 * t440 + t348 * t450 + t350 * t465 + t352 * t466 + t56 * t353 + t100 / 0.2e1 + t422 / 0.2e1 - t420 / 0.2e1 - t60 * mrSges(5,3) + t171 * mrSges(5,2) + (-t13 * t315 - t14 * t319) * mrSges(6,3) + (m(5) * t60 - m(6) * t56 + t406) * t322 + (mrSges(6,3) * t343 + t319 * t467 + t347 * t445 + t349 * t454 + t351 * t452 + t354 * t88 + t440 * t86) * qJD(5)) * t320 + (-t233 * t53 + t234 * t52) * mrSges(6,3) + (-t151 / 0.2e1 + t136 / 0.2e1) * t48 + t271 * t11 - t252 * t251 + t238 * t74 + t239 * t75 - t88 * (-mrSges(6,1) * t233 + mrSges(6,2) * t234) - t234 * t86 / 0.2e1 - t222 * t163 - t140 * t187 - t139 * t188 + m(5) * (t171 * qJ(3) + t186 * qJD(3) + t400 * t59) + (-mrSges(7,1) * t397 + mrSges(7,2) * t398) * t71 + t392 * t254 + t393 * qJD(3) + (-t150 / 0.2e1 + t138 / 0.2e1) * t47 + t120 * t28 + t121 * t29 - t122 * t111 + qJ(3) * t112 - m(6) * (t122 * t88 + t52 * t76 + t53 * t77) - m(5) * (t103 * t139 + t104 * t140 + t186 * t222) + (t93 / 0.2e1 + t94 / 0.2e1 - t421 / 0.2e1 + t171 * mrSges(5,1) + t7 / 0.2e1 + t30 / 0.2e1 + t31 / 0.2e1 - t99 / 0.2e1 + t36 / 0.2e1 + t322 * t162 - t59 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t471) * t182 + t490) * t316 + (Ifges(6,5) * t234 + Ifges(6,6) * t233) * t445 + (Ifges(7,5) * t136 + Ifges(7,6) * t138) * t447 + (Ifges(7,5) * t151 + Ifges(7,6) * t150) * t448 + (Ifges(6,1) * t234 + Ifges(6,4) * t233) * t452 + (-pkin(2) * t243 - qJ(3) * t205 - t213 * t254 + t220 * t496 - t230 * t252) * m(4) + t499 * t57 + t500 * t82 + t501 * t81 + (t120 * t3 + t121 * t2 + t15 * t500 + t16 * t501 + t271 * t34 + t499 * t71) * m(7) + (Ifges(6,4) * t234 + Ifges(6,2) * t233) * t454 + (Ifges(7,1) * t136 + Ifges(7,4) * t138) * t455 + (Ifges(7,1) * t151 + Ifges(7,4) * t150) * t456 + (Ifges(7,4) * t136 + Ifges(7,2) * t138) * t457 + (Ifges(7,4) * t151 + Ifges(7,2) * t150) * t458 + t233 * t467 - t260 * t474 - t258 * t475; -t232 * t128 - t231 * t129 - t257 * t28 - t259 * t29 + t482 * t82 + t481 * t81 - t393 * t304 + (mrSges(4,1) * t388 + t251 * t317) * t391 + (t187 * t372 - t11 + (t128 * t319 - t129 * t315 + t187) * qJD(4) + t406) * t320 + (t339 * qJD(5) + t162 + t342 + t289 * (t57 - t394)) * t316 + (-t2 * t259 - t257 * t3 - t320 * t34 + (t361 + t386) * t71 + t481 * t16 + t482 * t15) * m(7) + ((-qJD(4) * t343 - t56) * t320 + (qJD(4) * t88 - qJD(5) * t344 + t345) * t316 - t231 * t52 - t232 * t53 + t88 * t361) * m(6) + (-t186 * t304 - t289 * t340 + t316 * t59 + t320 * t60) * m(5) + (t220 * t304 + t230 * t372 + t243) * m(4); (Ifges(7,5) * t273 - Ifges(7,6) * t338 + t347) * t450 + t34 * (mrSges(7,1) * t338 + mrSges(7,2) * t273) + (t15 * t395 - t16 * t396 - t2 * t338 - t273 * t3) * mrSges(7,3) + (Ifges(7,4) * t273 - Ifges(7,2) * t338) * t469 + (Ifges(7,1) * t273 - Ifges(7,4) * t338) * t470 - t338 * t475 + (-Ifges(7,5) * t147 - Ifges(7,6) * t146) * t448 + (-Ifges(7,1) * t147 - Ifges(7,4) * t146) * t456 + (-Ifges(7,4) * t147 - Ifges(7,2) * t146) * t458 + (-t201 / 0.2e1 + t146 / 0.2e1) * t47 + ((-m(6) * t344 + t339) * qJD(5) + m(6) * t345 + t342) * pkin(10) + (-t498 / 0.2e1 + t147 / 0.2e1) * t48 + (-Ifges(7,5) * t498 - Ifges(7,6) * t201) * t447 + (-Ifges(7,1) * t498 - Ifges(7,4) * t201) * t455 + (-Ifges(7,4) * t498 - Ifges(7,2) * t201) * t457 + t324 * t241 + t484 * t82 + t485 * t81 + (t15 * t484 + t16 * t485 + t2 * t227 + t226 * t3 + t309 * t34 - t71 * t83) * m(7) + t284 + t309 * t11 + t226 * t28 + t227 * t29 - t103 * t187 + (mrSges(7,1) * t396 - mrSges(7,2) * t395) * t71 + t394 * t104 - t70 * t128 - t69 * t129 - t83 * t57 - pkin(4) * t50 + t478 + (-pkin(4) * t56 - t104 * t88 - t52 * t69 - t53 * t70) * m(6) + t37 * t439 + ((-Ifges(5,1) / 0.2e1 + t471) * t241 + t480 + t497) * t240 + (t437 * t483 - t497) * qJD(5) + t345 * mrSges(6,3) - t56 * t354 + t349 * t465 + t351 * t466 + t315 * t468 + t273 * t474; t109 * t494 + (-Ifges(6,2) * t173 + t169 + t86) * t454 + t490 + t491 + (t318 * t28 + t314 * t29 + m(7) * (t2 * t314 + t3 * t318) - t483 * t173 + (-t314 * t82 + t318 * t81 + m(7) * (-t15 * t314 + t16 * t318)) * qJD(6)) * pkin(5) - t88 * (mrSges(6,1) * t173 + mrSges(6,2) * t172) - t52 * t128 + t53 * t129 - t19 * t81 - t18 * t82 - m(7) * (t15 * t18 + t16 * t19) + (Ifges(6,5) * t172 - Ifges(6,6) * t173) * t445 + t85 * t451 + (t172 * t52 + t173 * t53) * mrSges(6,3) + (Ifges(6,1) * t172 - t433) * t452 + t505; -t15 * t81 + t16 * t82 + t47 * t455 + t502 + t505 + t7;];
tauc  = t1(:);
