% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:38:52
% EndTime: 2018-11-23 16:39:13
% DurationCPUTime: 22.06s
% Computational Cost: add. (41529->805), mult. (135677->1137), div. (0->0), fcn. (116932->14), ass. (0->383)
t308 = sin(pkin(13));
t310 = sin(pkin(6));
t311 = cos(pkin(13));
t321 = cos(qJ(3));
t312 = cos(pkin(7));
t317 = sin(qJ(3));
t402 = t312 * t317;
t332 = t310 * (-t308 * t402 + t311 * t321);
t277 = qJD(1) * t332;
t309 = sin(pkin(7));
t389 = qJD(3) * t321;
t548 = -t309 * t389 + t277;
t313 = cos(pkin(6));
t404 = t310 * t311;
t339 = t309 * t313 + t312 * t404;
t330 = t339 * qJD(1);
t392 = qJD(1) * t310;
t374 = t308 * t392;
t255 = -t317 * t374 + t321 * t330;
t315 = sin(qJ(5));
t316 = sin(qJ(4));
t319 = cos(qJ(5));
t320 = cos(qJ(4));
t290 = t315 * t320 + t316 * t319;
t206 = t290 * t255;
t503 = qJD(4) + qJD(5);
t268 = t503 * t290;
t397 = t206 - t268;
t289 = t315 * t316 - t319 * t320;
t207 = t289 * t255;
t267 = t503 * t289;
t396 = -t207 + t267;
t406 = t309 * t317;
t286 = t312 * t320 - t316 * t406;
t365 = t309 * t374;
t509 = qJD(4) * t286 - t316 * t365 - t548 * t320;
t287 = t312 * t316 + t320 * t406;
t508 = -qJD(4) * t287 + t548 * t316 - t320 * t365;
t483 = -pkin(11) - pkin(10);
t375 = qJD(4) * t483;
t293 = t316 * t375;
t297 = t483 * t316;
t298 = t483 * t320;
t342 = t319 * t297 + t298 * t315;
t366 = t320 * t375;
t226 = qJD(5) * t342 + t319 * t293 + t315 * t366;
t301 = qJ(2) * t404;
t456 = pkin(1) * t313;
t380 = qJD(1) * t456;
t282 = qJD(1) * t301 + t308 * t380;
t245 = pkin(9) * t330 + t282;
t300 = t311 * t380;
t408 = t308 * t310;
t329 = pkin(2) * t313 + (-pkin(9) * t312 - qJ(2)) * t408;
t254 = qJD(1) * t329 + t300;
t278 = (-pkin(9) * t308 * t309 - pkin(2) * t311 - pkin(1)) * t310;
t272 = qJD(1) * t278 + qJD(2);
t346 = t254 * t312 + t272 * t309;
t178 = -t317 * t245 + t346 * t321;
t265 = t313 * t406 + (t308 * t321 + t311 * t402) * t310;
t256 = t265 * qJD(1);
t212 = pkin(3) * t256 - pkin(10) * t255;
t127 = -t178 * t316 + t320 * t212;
t112 = -pkin(11) * t255 * t320 + pkin(4) * t256 + t127;
t128 = t320 * t178 + t316 * t212;
t409 = t255 * t316;
t121 = -pkin(11) * t409 + t128;
t69 = t315 * t112 + t319 * t121;
t547 = t226 - t69;
t235 = t321 * t245;
t401 = t312 * t321;
t331 = (t308 * t401 + t311 * t317) * t310;
t328 = qJD(2) * t331;
t149 = qJD(1) * t328 + (t317 * t346 + t235) * qJD(3);
t283 = -t309 * t404 + t312 * t313;
t279 = qJD(1) * t283 + qJD(3);
t221 = t256 * t320 + t279 * t316;
t405 = t309 * t321;
t532 = t310 * (-t308 * t317 + t311 * t401) + t313 * t405;
t257 = t532 * qJD(3);
t242 = qJD(1) * t257;
t186 = -qJD(4) * t221 - t242 * t316;
t117 = -t186 * pkin(4) + t149;
t258 = t265 * qJD(3);
t243 = qJD(1) * t258;
t466 = t243 / 0.2e1;
t538 = -Ifges(6,4) / 0.2e1;
t220 = -t256 * t316 + t279 * t320;
t185 = qJD(4) * t220 + t242 * t320;
t367 = t319 * t220 - t221 * t315;
t93 = qJD(5) * t367 + t185 * t319 + t186 * t315;
t348 = t220 * t315 + t319 * t221;
t94 = qJD(5) * t348 + t185 * t315 - t319 * t186;
t546 = t117 * mrSges(6,2) + t93 * Ifges(6,1) + 0.2e1 * Ifges(6,5) * t466 + 0.2e1 * t94 * t538;
t544 = -t255 / 0.2e1;
t543 = -t279 / 0.2e1;
t179 = t254 * t402 + t272 * t406 + t235;
t142 = pkin(4) * t409 + t179;
t388 = qJD(4) * t316;
t542 = pkin(4) * t388 - t397 * pkin(5) + t396 * pkin(12) - t142;
t541 = -pkin(12) * t256 + t547;
t253 = -t255 + qJD(4);
t248 = qJD(5) + t253;
t314 = sin(qJ(6));
t318 = cos(qJ(6));
t132 = t248 * t318 - t314 * t348;
t60 = qJD(6) * t132 + t243 * t314 + t318 * t93;
t133 = t248 * t314 + t318 * t348;
t61 = -qJD(6) * t133 + t243 * t318 - t314 * t93;
t23 = -mrSges(7,1) * t61 + mrSges(7,2) * t60;
t213 = -t254 * t309 + t312 * t272;
t153 = -pkin(3) * t255 - pkin(10) * t256 + t213;
t155 = pkin(10) * t279 + t179;
t110 = t153 * t316 + t155 * t320;
t92 = pkin(11) * t220 + t110;
t414 = t319 * t92;
t109 = t320 * t153 - t155 * t316;
t91 = -pkin(11) * t221 + t109;
t84 = pkin(4) * t253 + t91;
t40 = t315 * t84 + t414;
t327 = qJD(2) * t332;
t148 = qJD(1) * t327 + qJD(3) * t178;
t391 = qJD(2) * t310;
t373 = t308 * t391;
t363 = qJD(1) * t373;
t341 = t309 * t363;
t200 = pkin(3) * t243 - pkin(10) * t242 + t341;
t67 = -qJD(4) * t110 - t148 * t316 + t320 * t200;
t47 = pkin(4) * t243 - pkin(11) * t185 + t67;
t387 = qJD(4) * t320;
t66 = t320 * t148 + t153 * t387 - t155 * t388 + t316 * t200;
t53 = pkin(11) * t186 + t66;
t11 = -qJD(5) * t40 - t315 * t53 + t319 * t47;
t8 = -pkin(5) * t243 - t11;
t540 = m(7) * t8 + t23;
t343 = t319 * t286 - t287 * t315;
t510 = qJD(5) * t343 + t508 * t315 + t509 * t319;
t276 = qJD(1) * t331;
t390 = qJD(3) * t317;
t539 = t309 * t390 - t276;
t275 = t297 * t315 - t298 * t319;
t227 = qJD(5) * t275 + t293 * t315 - t319 * t366;
t68 = t112 * t319 - t121 * t315;
t537 = -t227 - t68;
t393 = t308 * t456 + t301;
t260 = pkin(9) * t339 + t393;
t303 = t311 * t456;
t266 = t303 + t329;
t345 = t266 * t312 + t278 * t309;
t191 = -t317 * t260 + t345 * t321;
t58 = Ifges(7,6) * t61;
t59 = Ifges(7,5) * t60;
t18 = Ifges(7,3) * t94 + t58 + t59;
t38 = pkin(12) * t248 + t40;
t154 = -t279 * pkin(3) - t178;
t124 = -t220 * pkin(4) + t154;
t79 = -pkin(5) * t367 - pkin(12) * t348 + t124;
t21 = -t314 * t38 + t318 * t79;
t30 = t94 * pkin(5) - t93 * pkin(12) + t117;
t385 = qJD(5) * t319;
t386 = qJD(5) * t315;
t10 = t315 * t47 + t319 * t53 + t84 * t385 - t386 * t92;
t7 = pkin(12) * t243 + t10;
t2 = qJD(6) * t21 + t30 * t314 + t318 * t7;
t22 = t314 * t79 + t318 * t38;
t3 = -qJD(6) * t22 + t30 * t318 - t314 * t7;
t362 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t528 = t93 * t538;
t501 = -t243 * Ifges(6,6) / 0.2e1 + t528;
t484 = t94 / 0.2e1;
t502 = Ifges(6,2) * t484;
t536 = t362 + t117 * mrSges(6,1) + t18 / 0.2e1 + t501 + t502;
t534 = Ifges(4,4) * t544 + Ifges(4,5) * t543;
t533 = Ifges(4,6) * t543 - t256 * Ifges(4,4) / 0.2e1;
t531 = t256 * Ifges(4,1) / 0.2e1 - t534;
t530 = Ifges(4,2) * t544 + t533;
t472 = t185 / 0.2e1;
t471 = t186 / 0.2e1;
t526 = Ifges(6,5) * t93;
t525 = Ifges(6,6) * t94;
t524 = Ifges(4,5) * t242;
t522 = Ifges(4,6) * t243;
t520 = t148 * mrSges(4,2);
t222 = -t266 * t309 + t312 * t278;
t169 = -pkin(3) * t532 - pkin(10) * t265 + t222;
t251 = t321 * t260;
t192 = t266 * t402 + t278 * t406 + t251;
t177 = pkin(10) * t283 + t192;
t119 = t316 * t169 + t320 * t177;
t228 = -t265 * t316 + t283 * t320;
t105 = pkin(11) * t228 + t119;
t118 = t320 * t169 - t177 * t316;
t229 = t265 * t320 + t283 * t316;
t97 = -pkin(4) * t532 - pkin(11) * t229 + t118;
t518 = t319 * t105 + t315 * t97;
t517 = pkin(5) * t256 - t537;
t307 = -pkin(4) * t320 - pkin(3);
t261 = pkin(5) * t289 - pkin(12) * t290 + t307;
t218 = t261 * t314 + t275 * t318;
t516 = -qJD(6) * t218 - t314 * t541 + t318 * t542;
t217 = t261 * t318 - t275 * t314;
t515 = qJD(6) * t217 + t314 * t542 + t318 * t541;
t136 = mrSges(6,1) * t248 - mrSges(6,3) * t348;
t87 = -mrSges(7,1) * t132 + mrSges(7,2) * t133;
t412 = t87 - t136;
t247 = t286 * t315 + t287 * t319;
t340 = -t318 * t247 + t314 * t405;
t513 = qJD(6) * t340 - t314 * t510 + t318 * t539;
t230 = -t314 * t247 - t318 * t405;
t512 = qJD(6) * t230 + t314 * t539 + t318 * t510;
t511 = -qJD(5) * t247 - t315 * t509 + t319 * t508;
t395 = mrSges(4,1) * t279 + mrSges(5,1) * t220 - mrSges(5,2) * t221 - mrSges(4,3) * t256;
t156 = Ifges(6,4) * t367;
t157 = qJD(6) - t367;
t355 = Ifges(7,5) * t318 - Ifges(7,6) * t314;
t335 = t157 * t355;
t441 = Ifges(7,4) * t314;
t359 = Ifges(7,1) * t318 - t441;
t336 = t133 * t359;
t440 = Ifges(7,4) * t318;
t357 = -Ifges(7,2) * t314 + t440;
t337 = t132 * t357;
t360 = mrSges(7,1) * t314 + mrSges(7,2) * t318;
t416 = t315 * t92;
t39 = t319 * t84 - t416;
t37 = -pkin(5) * t248 - t39;
t338 = t37 * t360;
t439 = Ifges(6,5) * t248;
t459 = -t318 / 0.2e1;
t460 = t314 / 0.2e1;
t446 = Ifges(6,1) * t348;
t108 = t156 + t439 + t446;
t499 = t124 * mrSges(6,2) + t108 / 0.2e1;
t442 = Ifges(7,4) * t133;
t75 = Ifges(7,2) * t132 + Ifges(7,6) * t157 + t442;
t131 = Ifges(7,4) * t132;
t76 = Ifges(7,1) * t133 + Ifges(7,5) * t157 + t131;
t507 = t39 * mrSges(6,3) - t439 / 0.2e1 - t337 / 0.2e1 - t336 / 0.2e1 - t335 / 0.2e1 + t75 * t460 + t76 * t459 - t338 - t499 - t156 / 0.2e1;
t505 = -t316 * t67 + t320 * t66;
t504 = -t21 * t314 + t22 * t318;
t425 = t221 * Ifges(5,4);
t140 = t220 * Ifges(5,2) + t253 * Ifges(5,6) + t425;
t219 = Ifges(5,4) * t220;
t141 = t221 * Ifges(5,1) + t253 * Ifges(5,5) + t219;
t349 = t109 * t320 + t110 * t316;
t444 = Ifges(5,4) * t320;
t445 = Ifges(5,4) * t316;
t457 = t320 / 0.2e1;
t464 = t253 / 0.2e1;
t467 = t221 / 0.2e1;
t469 = t220 / 0.2e1;
t500 = -t349 * mrSges(5,3) + t154 * (mrSges(5,1) * t316 + mrSges(5,2) * t320) + (Ifges(5,5) * t320 - Ifges(5,6) * t316) * t464 + (-Ifges(5,2) * t316 + t444) * t469 + (Ifges(5,1) * t320 - t445) * t467 - t316 * t140 / 0.2e1 + t141 * t457;
t198 = qJD(4) * t228 + t257 * t320;
t163 = qJD(3) * t191 + t327;
t364 = t309 * t373;
t205 = pkin(3) * t258 - pkin(10) * t257 + t364;
t78 = -qJD(4) * t119 - t163 * t316 + t320 * t205;
t55 = pkin(4) * t258 - pkin(11) * t198 + t78;
t199 = -qJD(4) * t229 - t257 * t316;
t77 = t320 * t163 + t169 * t387 - t177 * t388 + t316 * t205;
t63 = pkin(11) * t199 + t77;
t15 = -qJD(5) * t518 - t315 * t63 + t319 * t55;
t497 = t67 * mrSges(5,1) - t66 * mrSges(5,2);
t496 = t11 * mrSges(6,1) - t10 * mrSges(6,2);
t114 = pkin(5) * t348 - pkin(12) * t367;
t495 = t213 * mrSges(4,2) - t178 * mrSges(4,3) + t531;
t436 = Ifges(6,6) * t248;
t437 = Ifges(6,2) * t367;
t443 = Ifges(6,4) * t348;
t107 = t436 + t437 + t443;
t434 = Ifges(7,3) * t157;
t435 = Ifges(7,6) * t132;
t438 = Ifges(7,5) * t133;
t74 = t434 + t435 + t438;
t494 = t22 * mrSges(7,2) + t107 / 0.2e1 - t74 / 0.2e1 - t124 * mrSges(6,1) - t21 * mrSges(7,1);
t493 = t213 * mrSges(4,1) + t109 * mrSges(5,1) + t39 * mrSges(6,1) - t110 * mrSges(5,2) - t40 * mrSges(6,2) - t179 * mrSges(4,3) + t530;
t491 = Ifges(6,2) / 0.2e1;
t487 = t60 / 0.2e1;
t486 = t61 / 0.2e1;
t485 = t76 / 0.2e1;
t482 = Ifges(5,1) * t472 + Ifges(5,4) * t471 + Ifges(5,5) * t466;
t481 = -t132 / 0.2e1;
t480 = t132 / 0.2e1;
t479 = -t133 / 0.2e1;
t478 = t133 / 0.2e1;
t476 = -t157 / 0.2e1;
t475 = t157 / 0.2e1;
t474 = t367 / 0.2e1;
t473 = t348 / 0.2e1;
t468 = -t221 / 0.2e1;
t465 = t248 / 0.2e1;
t461 = -t314 / 0.2e1;
t458 = t318 / 0.2e1;
t455 = t2 * t318;
t454 = t3 * t314;
t81 = mrSges(6,1) * t243 - mrSges(6,3) * t93;
t449 = t23 - t81;
t448 = mrSges(4,3) * t242;
t447 = mrSges(4,3) * t243;
t433 = t367 * Ifges(6,6);
t432 = t348 * Ifges(6,5);
t428 = t21 * t318;
t426 = t220 * Ifges(5,6);
t424 = t221 * Ifges(5,5);
t421 = t248 * Ifges(6,3);
t420 = t253 * Ifges(5,3);
t411 = t149 * t321;
t403 = t311 * (-mrSges(3,2) * t313 + mrSges(3,3) * t404) * qJD(1);
t400 = t314 * t267;
t399 = t318 * t267;
t394 = -t522 + t524;
t384 = qJD(6) * t314;
t383 = qJD(6) * t318;
t381 = Ifges(6,3) * t243 - t525 + t526;
t113 = -mrSges(6,1) * t367 + mrSges(6,2) * t348;
t377 = -t113 + t395;
t376 = Ifges(5,5) * t185 + Ifges(5,6) * t186 + Ifges(5,3) * t243;
t170 = t207 * t314 + t256 * t318;
t369 = t170 - t400;
t171 = -t207 * t318 + t256 * t314;
t368 = -t171 - t399;
t361 = mrSges(7,1) * t318 - mrSges(7,2) * t314;
t358 = Ifges(7,1) * t314 + t440;
t356 = Ifges(7,2) * t318 + t441;
t354 = Ifges(7,5) * t314 + Ifges(7,6) * t318;
t353 = t22 * t314 + t428;
t49 = -pkin(12) * t532 + t518;
t176 = -t283 * pkin(3) - t191;
t134 = -t228 * pkin(4) + t176;
t181 = t228 * t315 + t229 * t319;
t347 = t319 * t228 - t229 * t315;
t80 = -pkin(5) * t347 - t181 * pkin(12) + t134;
t25 = t314 * t80 + t318 * t49;
t24 = -t314 * t49 + t318 * t80;
t50 = -t105 * t315 + t319 * t97;
t144 = t181 * t318 - t314 * t532;
t143 = -t181 * t314 - t318 * t532;
t344 = -(-qJ(2) * t374 + t300) * t308 + t282 * t311;
t14 = -t105 * t386 + t315 * t55 + t319 * t63 + t97 * t385;
t284 = (mrSges(3,1) * t313 - mrSges(3,3) * t408) * qJD(1);
t19 = t60 * Ifges(7,4) + t61 * Ifges(7,2) + t94 * Ifges(7,6);
t20 = t60 * Ifges(7,1) + t61 * Ifges(7,4) + t94 * Ifges(7,5);
t326 = mrSges(7,3) * t455 + qJD(6) * t338 + t19 * t458 + t20 * t460 + t356 * t486 + t358 * t487 - t8 * t361 + t381 - t75 * t384 / 0.2e1 + t383 * t485 + t354 * t484 + t496 + (t337 + t336 + t335) * qJD(6) / 0.2e1;
t164 = t328 + (t317 * t345 + t251) * qJD(3);
t122 = -t199 * pkin(4) + t164;
t102 = -mrSges(7,2) * t157 + mrSges(7,3) * t132;
t103 = mrSges(7,1) * t157 - mrSges(7,3) * t133;
t31 = mrSges(7,1) * t94 - mrSges(7,3) * t60;
t32 = -mrSges(7,2) * t94 + mrSges(7,3) * t61;
t325 = -t102 * t384 - t103 * t383 - t314 * t31 + t318 * t32 + m(7) * (-t21 * t383 - t22 * t384 - t454 + t455);
t324 = -t434 / 0.2e1 - t435 / 0.2e1 - t438 / 0.2e1 + t436 / 0.2e1 + t40 * mrSges(6,3) + t443 / 0.2e1 + t494;
t223 = -mrSges(4,2) * t279 + mrSges(4,3) * t255;
t211 = -mrSges(4,1) * t255 + mrSges(4,2) * t256;
t208 = mrSges(4,1) * t243 + mrSges(4,2) * t242;
t190 = mrSges(5,1) * t253 - mrSges(5,3) * t221;
t189 = -mrSges(5,2) * t253 + mrSges(5,3) * t220;
t147 = -mrSges(5,2) * t243 + mrSges(5,3) * t186;
t146 = mrSges(5,1) * t243 - mrSges(5,3) * t185;
t139 = t420 + t424 + t426;
t135 = -mrSges(6,2) * t248 + mrSges(6,3) * t367;
t123 = -mrSges(5,1) * t186 + mrSges(5,2) * t185;
t115 = t185 * Ifges(5,4) + t186 * Ifges(5,2) + t243 * Ifges(5,6);
t106 = t421 + t432 + t433;
t101 = qJD(5) * t181 + t198 * t315 - t319 * t199;
t100 = qJD(5) * t347 + t198 * t319 + t199 * t315;
t98 = pkin(4) * t221 + t114;
t82 = -mrSges(6,2) * t243 - mrSges(6,3) * t94;
t71 = -qJD(6) * t144 - t100 * t314 + t258 * t318;
t70 = qJD(6) * t143 + t100 * t318 + t258 * t314;
t48 = pkin(5) * t532 - t50;
t45 = mrSges(6,1) * t94 + mrSges(6,2) * t93;
t44 = t319 * t91 - t416;
t43 = t315 * t91 + t414;
t35 = t101 * pkin(5) - t100 * pkin(12) + t122;
t29 = t114 * t314 + t318 * t39;
t28 = t114 * t318 - t314 * t39;
t27 = t314 * t98 + t318 * t44;
t26 = -t314 * t44 + t318 * t98;
t13 = -pkin(5) * t258 - t15;
t12 = pkin(12) * t258 + t14;
t5 = -qJD(6) * t25 - t12 * t314 + t318 * t35;
t4 = qJD(6) * t24 + t12 * t318 + t314 * t35;
t1 = [(Ifges(5,1) * t198 + Ifges(5,4) * t199) * t467 + (Ifges(5,1) * t229 + Ifges(5,4) * t228) * t472 + (Ifges(7,5) * t70 + Ifges(7,6) * t71) * t475 + (Ifges(7,5) * t144 + Ifges(7,6) * t143) * t484 + (Ifges(7,4) * t70 + Ifges(7,2) * t71) * t480 + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t486 + (-m(4) * t178 + m(5) * t154 - t395) * t164 + (Ifges(6,1) * t473 + Ifges(6,4) * t474 + Ifges(6,5) * t465 + t499) * t100 - (t381 + t376) * t532 / 0.2e1 - (mrSges(4,1) * t341 - t148 * mrSges(4,3) + Ifges(4,2) * t243 - Ifges(4,4) * t242 + t526 / 0.2e1 - t525 / 0.2e1 + Ifges(5,6) * t471 + Ifges(5,5) * t472 + (Ifges(6,3) + Ifges(5,3)) * t466 + t496 + t497) * t532 + (t143 * t2 - t144 * t3 - t21 * t70 + t22 * t71) * mrSges(7,3) + (-t109 * t198 + t110 * t199 + t228 * t66 - t229 * t67) * mrSges(5,3) + m(6) * (t10 * t518 + t11 * t50 + t117 * t134 + t122 * t124 + t14 * t40 + t15 * t39) + t518 * t82 + (t10 * t347 - t100 * t39 - t101 * t40 - t11 * t181) * mrSges(6,3) + (Ifges(5,4) * t198 + Ifges(5,2) * t199) * t469 + (Ifges(5,4) * t229 + Ifges(5,2) * t228) * t471 + (Ifges(5,5) * t467 + Ifges(6,5) * t473 + Ifges(5,6) * t469 + Ifges(6,6) * t474 + Ifges(5,3) * t464 + Ifges(6,3) * t465 + t493 + t530) * t258 + (t495 + t531) * t257 + t211 * t364 + (t139 + t106) * t258 / 0.2e1 + (m(3) * ((t311 * t393 + (qJ(2) * t408 - t303) * t308) * qJD(1) + t344) + 0.2e1 * t403) * t391 + t546 * t181 + (Ifges(7,1) * t70 + Ifges(7,4) * t71) * t478 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t487 + (Ifges(5,5) * t198 + Ifges(5,6) * t199) * t464 + (Ifges(5,5) * t229 + Ifges(5,6) * t228) * t466 + (-Ifges(6,4) * t473 + Ifges(7,5) * t478 - Ifges(6,2) * t474 - Ifges(6,6) * t465 + Ifges(7,6) * t480 + Ifges(7,3) * t475 - t494) * t101 + (mrSges(4,2) * t341 + mrSges(4,3) * t149 + Ifges(4,1) * t242 - Ifges(4,4) * t243) * t265 + t149 * (-mrSges(5,1) * t228 + mrSges(5,2) * t229) + t163 * t223 + t228 * t115 / 0.2e1 + t222 * t208 + m(7) * (t13 * t37 + t2 * t25 + t21 * t5 + t22 * t4 + t24 * t3 + t48 * t8) - t192 * t447 - t191 * t448 + m(4) * (t148 * t192 - t149 * t191 + t163 * t179 + (qJD(1) * t222 + t213) * t364) + t154 * (-mrSges(5,1) * t199 + mrSges(5,2) * t198) + t199 * t140 / 0.2e1 + t198 * t141 / 0.2e1 + t77 * t189 + t78 * t190 + t176 * t123 - (Ifges(7,5) * t487 - Ifges(6,6) * t466 + Ifges(7,6) * t486 + Ifges(7,3) * t484 + t502 + t528 + t536) * t347 + (-t520 - t522 / 0.2e1 + t524 / 0.2e1 - t149 * mrSges(4,1) + t394 / 0.2e1) * t283 + t24 * t31 + t25 * t32 + t48 * t23 + t37 * (-mrSges(7,1) * t71 + mrSges(7,2) * t70) + t71 * t75 / 0.2e1 + t50 * t81 + m(5) * (t109 * t78 + t110 * t77 + t118 * t67 + t119 * t66 + t149 * t176) + t13 * t87 + t4 * t102 + t5 * t103 - 0.2e1 * t284 * t373 + t122 * t113 + t134 * t45 + t14 * t135 + t15 * t136 + t143 * t19 / 0.2e1 + t144 * t20 / 0.2e1 + t8 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + t118 * t146 + t119 * t147 + t229 * t482 + t70 * t485; -m(4) * (-t178 * t276 + t179 * t277) + t508 * t190 + t509 * t189 + t312 * t208 + t286 * t146 + t287 * t147 - t277 * t223 + t247 * t82 + t230 * t31 - t340 * t32 + t510 * t135 + t513 * t103 - t449 * t343 + t512 * t102 + t377 * t276 - t412 * t511 + (-m(3) * t344 + t308 * t284 - t403) * t392 + (-t2 * t340 + t21 * t513 + t22 * t512 + t230 * t3 - t343 * t8 - t37 * t511) * m(7) + (t10 * t247 + t11 * t343 - t124 * t276 + t39 * t511 + t40 * t510) * m(6) + (t508 * t109 + t509 * t110 - t154 * t276 + t67 * t286 + t66 * t287) * m(5) + (-t211 * t374 - t317 * t447 + (-t123 - t45 - t448) * t321 + (t223 * t321 - t317 * t377) * qJD(3) + m(6) * (-t117 * t321 + t124 * t390) + m(5) * (t154 * t390 - t411) + (t148 * t317 - t178 * t390 + t179 * t389 - t213 * t374 + t312 * t363 - t411) * m(4)) * t309; (m(5) * t505 - t146 * t316 + t147 * t320 + (-m(5) * t349 - t316 * t189 - t320 * t190) * qJD(4)) * pkin(10) + t505 * mrSges(5,3) + ((m(6) * t124 + t113) * t316 * pkin(4) + t500) * qJD(4) + (-pkin(3) * t149 - t109 * t127 - t110 * t128 - t154 * t179) * m(5) + (-Ifges(7,5) * t399 + Ifges(7,6) * t400 + Ifges(7,3) * t268) * t475 + (-Ifges(7,1) * t399 + Ifges(7,4) * t400 + Ifges(7,5) * t268) * t478 + (-Ifges(7,4) * t399 + Ifges(7,2) * t400 + Ifges(7,6) * t268) * t480 + (-t170 / 0.2e1 + t400 / 0.2e1) * t75 - t449 * t342 + m(6) * (t10 * t275 + t11 * t342 + t117 * t307 + t226 * t40 - t227 * t39) + (t2 * t218 + t21 * t516 + t217 * t3 + t22 * t515 - t342 * t8 + t37 * t517) * m(7) + (-t267 / 0.2e1 + t207 / 0.2e1) * t108 + (-Ifges(6,5) * t267 - Ifges(6,6) * t268) * t465 + (-Ifges(6,1) * t267 - Ifges(6,4) * t268) * t473 + (-Ifges(6,4) * t267 - Ifges(6,2) * t268) * t474 + (-t171 / 0.2e1 - t399 / 0.2e1) * t76 - t520 + t547 * t135 + (t74 - t107) * (t268 / 0.2e1 - t206 / 0.2e1) + (t359 * t487 + t357 * t486 + t355 * t484 + t8 * t360 - t11 * mrSges(6,3) + t19 * t461 + t20 * t458 + (-t2 * t314 - t3 * t318) * mrSges(7,3) + (-mrSges(7,3) * t504 + t354 * t476 + t356 * t481 + t358 * t479 + t37 * t361 + t75 * t459 + t76 * t461) * qJD(6) + t546) * t290 - m(6) * (t124 * t142 + t39 * t68 + t40 * t69) + (-mrSges(5,1) * t320 + mrSges(5,2) * t316 - mrSges(4,1)) * t149 + t394 + (-t493 - t420 / 0.2e1 - t139 / 0.2e1 - t424 / 0.2e1 - t426 / 0.2e1 - t421 / 0.2e1 - t106 / 0.2e1 - t432 / 0.2e1 - t433 / 0.2e1 - t533) * t256 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t256 - t495 - t500 + t534) * t255 + t307 * t45 + t275 * t82 - t178 * t223 + t217 * t31 + t218 * t32 - t128 * t189 - t127 * t190 - t367 * (-Ifges(6,4) * t207 - Ifges(6,2) * t206) / 0.2e1 - t348 * (-Ifges(6,1) * t207 - Ifges(6,4) * t206) / 0.2e1 - t248 * (-Ifges(6,5) * t207 - Ifges(6,6) * t206) / 0.2e1 + (-t10 * mrSges(6,3) + t59 / 0.2e1 + t58 / 0.2e1 + (t491 + Ifges(7,3) / 0.2e1) * t94 + t501 + t536) * t289 + t515 * t102 + t516 * t103 + t517 * t87 + t537 * t136 + t115 * t457 + (mrSges(7,1) * t369 + mrSges(7,2) * t368) * t37 - pkin(3) * t123 - t142 * t113 + t395 * t179 + (mrSges(7,2) * t397 - mrSges(7,3) * t369) * t22 + (t39 * t396 + t397 * t40) * mrSges(6,3) + (-mrSges(7,1) * t397 - mrSges(7,3) * t368) * t21 + (-mrSges(6,1) * t397 - mrSges(6,2) * t396) * t124 + (Ifges(5,5) * t316 + Ifges(5,6) * t320) * t466 + (Ifges(5,2) * t320 + t445) * t471 + (Ifges(5,1) * t316 + t444) * t472 + (Ifges(7,5) * t171 + Ifges(7,6) * t170 + Ifges(7,3) * t206) * t476 + (Ifges(7,1) * t171 + Ifges(7,4) * t170 + Ifges(7,5) * t206) * t479 + (Ifges(7,4) * t171 + Ifges(7,2) * t170 + Ifges(7,6) * t206) * t481 + t316 * t482; (-t446 / 0.2e1 + t507) * t367 + (-t157 * t428 + (-t157 * t22 - t3) * t314) * mrSges(7,3) + (t437 / 0.2e1 + t324) * t348 + t325 * (pkin(4) * t315 + pkin(12)) + t376 + (t109 * t220 + t110 * t221) * mrSges(5,3) - m(6) * (-t39 * t43 + t40 * t44) + t497 + t326 - m(7) * (t21 * t26 + t22 * t27 + t37 * t43) - t253 * (Ifges(5,5) * t220 - Ifges(5,6) * t221) / 0.2e1 - t154 * (mrSges(5,1) * t221 + mrSges(5,2) * t220) - t109 * t189 + t110 * t190 - t412 * t43 + (-t221 * t113 + t315 * t82 + t319 * t81 + ((m(7) * t37 + t412) * t315 + (m(7) * t504 + t102 * t318 - t103 * t314 + t135) * t319) * qJD(5) + (0.2e1 * t124 * t468 + t10 * t315 + t11 * t319 + (-t315 * t39 + t319 * t40) * qJD(5)) * m(6)) * pkin(4) - t27 * t102 - t26 * t103 - t44 * t135 + t140 * t467 + (Ifges(5,1) * t220 - t425) * t468 - (-Ifges(5,2) * t221 + t141 + t219) * t220 / 0.2e1 + t540 * (-pkin(4) * t319 - pkin(5)); ((-Ifges(6,1) / 0.2e1 + t491) * t348 + t353 * mrSges(7,3) + t507) * t367 + t324 * t348 + t325 * pkin(12) - t412 * t40 + t326 - m(7) * (t21 * t28 + t22 * t29 + t37 * t40) + (-qJD(6) * t353 - t454) * mrSges(7,3) - t29 * t102 - t28 * t103 - t39 * t135 - t540 * pkin(5); -t21 * t102 + (Ifges(7,1) * t132 - t442) * t479 + t75 * t478 + (Ifges(7,5) * t132 - Ifges(7,6) * t133) * t476 + t22 * t103 - t37 * (mrSges(7,1) * t133 + mrSges(7,2) * t132) + (t132 * t21 + t133 * t22) * mrSges(7,3) + t362 + t18 + (-Ifges(7,2) * t133 + t131 + t76) * t481;];
tauc  = t1(:);
