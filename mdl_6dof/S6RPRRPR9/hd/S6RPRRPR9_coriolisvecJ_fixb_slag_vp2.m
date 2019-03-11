% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:26:57
% EndTime: 2019-03-09 05:27:41
% DurationCPUTime: 21.82s
% Computational Cost: add. (31319->772), mult. (103872->1088), div. (0->0), fcn. (89399->14), ass. (0->360)
t289 = sin(pkin(12));
t291 = sin(pkin(6));
t293 = cos(pkin(12));
t301 = cos(qJ(3));
t294 = cos(pkin(7));
t298 = sin(qJ(3));
t369 = t294 * t298;
t309 = (-t289 * t369 + t293 * t301) * t291;
t257 = qJD(1) * t309;
t290 = sin(pkin(7));
t356 = qJD(3) * t301;
t500 = -t290 * t356 + t257;
t295 = cos(pkin(6));
t368 = t294 * t301;
t372 = t290 * t301;
t304 = t295 * t372 + t291 * (-t289 * t298 + t293 * t368);
t237 = t304 * qJD(1);
t288 = sin(pkin(13));
t292 = cos(pkin(13));
t297 = sin(qJ(4));
t300 = cos(qJ(4));
t273 = t288 * t300 + t292 * t297;
t187 = t273 * t237;
t266 = t273 * qJD(4);
t365 = t187 - t266;
t272 = t288 * t297 - t292 * t300;
t188 = t272 * t237;
t267 = t272 * qJD(4);
t364 = -t188 + t267;
t373 = t290 * t298;
t270 = t294 * t300 - t297 * t373;
t359 = qJD(1) * t291;
t347 = t289 * t359;
t337 = t290 * t347;
t462 = qJD(4) * t270 - t297 * t337 - t500 * t300;
t271 = t294 * t297 + t300 * t373;
t461 = -qJD(4) * t271 + t500 * t297 - t300 * t337;
t414 = -qJ(5) - pkin(10);
t341 = qJD(4) * t414;
t264 = qJD(5) * t300 + t297 * t341;
t308 = -qJD(5) * t297 + t300 * t341;
t219 = t292 * t264 + t288 * t308;
t371 = t291 * t293;
t281 = qJ(2) * t371;
t418 = pkin(1) * t295;
t353 = qJD(1) * t418;
t263 = qJD(1) * t281 + t289 * t353;
t311 = (t290 * t295 + t294 * t371) * pkin(9);
t230 = qJD(1) * t311 + t263;
t280 = t293 * t353;
t375 = t289 * t291;
t307 = pkin(2) * t295 + (-pkin(9) * t294 - qJ(2)) * t375;
t236 = qJD(1) * t307 + t280;
t258 = (-pkin(9) * t289 * t290 - pkin(2) * t293 - pkin(1)) * t291;
t252 = qJD(1) * t258 + qJD(2);
t319 = t236 * t294 + t252 * t290;
t161 = -t298 * t230 + t319 * t301;
t247 = t295 * t373 + (t289 * t301 + t293 * t369) * t291;
t240 = t247 * qJD(1);
t193 = pkin(3) * t240 - pkin(10) * t237;
t111 = t300 * t161 + t297 * t193;
t376 = t237 * t297;
t100 = -qJ(5) * t376 + t111;
t110 = -t161 * t297 + t300 * t193;
t90 = -qJ(5) * t237 * t300 + pkin(4) * t240 + t110;
t52 = t292 * t100 + t288 * t90;
t499 = t219 - t52;
t265 = -t290 * t371 + t294 * t295;
t259 = qJD(1) * t265 + qJD(3);
t202 = -t240 * t297 + t259 * t300;
t238 = t304 * qJD(3);
t225 = qJD(1) * t238;
t168 = qJD(4) * t202 + t225 * t300;
t203 = t240 * t300 + t259 * t297;
t169 = -qJD(4) * t203 - t225 * t297;
t106 = t168 * t288 - t292 * t169;
t107 = t168 * t292 + t169 * t288;
t239 = t247 * qJD(3);
t226 = qJD(1) * t239;
t426 = t226 / 0.2e1;
t491 = -Ifges(6,4) / 0.2e1;
t215 = t301 * t230;
t310 = (t289 * t368 + t293 * t298) * t291;
t306 = qJD(2) * t310;
t132 = qJD(1) * t306 + (t298 * t319 + t215) * qJD(3);
t96 = -t169 * pkin(4) + t132;
t498 = t96 * mrSges(6,2) + t107 * Ifges(6,1) + 0.2e1 * Ifges(6,5) * t426 + 0.2e1 * t106 * t491;
t496 = -t237 / 0.2e1;
t495 = -t259 / 0.2e1;
t479 = Ifges(6,3) + Ifges(5,3);
t494 = -pkin(11) * t240 + t499;
t162 = t236 * t369 + t252 * t373 + t215;
t125 = pkin(4) * t376 + t162;
t355 = qJD(4) * t297;
t493 = pkin(4) * t355 - t365 * pkin(5) + t364 * pkin(11) - t125;
t463 = t461 * t288 + t462 * t292;
t256 = qJD(1) * t310;
t357 = qJD(3) * t298;
t492 = t290 * t357 - t256;
t218 = t264 * t288 - t292 * t308;
t51 = -t100 * t288 + t292 * t90;
t490 = -t218 - t51;
t360 = t289 * t418 + t281;
t243 = t311 + t360;
t283 = t293 * t418;
t248 = t283 + t307;
t318 = t248 * t294 + t258 * t290;
t172 = -t298 * t243 + t318 * t301;
t235 = qJD(4) - t237;
t296 = sin(qJ(6));
t299 = cos(qJ(6));
t320 = t202 * t288 + t292 * t203;
t117 = t235 * t296 + t299 * t320;
t65 = -qJD(6) * t117 - t107 * t296 + t226 * t299;
t62 = Ifges(7,6) * t65;
t116 = t235 * t299 - t296 * t320;
t64 = qJD(6) * t116 + t107 * t299 + t226 * t296;
t63 = Ifges(7,5) * t64;
t15 = Ifges(7,3) * t106 + t62 + t63;
t194 = -t236 * t290 + t294 * t252;
t135 = -pkin(3) * t237 - pkin(10) * t240 + t194;
t137 = pkin(10) * t259 + t162;
t87 = t300 * t135 - t137 * t297;
t73 = -qJ(5) * t203 + t87;
t70 = pkin(4) * t235 + t73;
t88 = t135 * t297 + t137 * t300;
t74 = qJ(5) * t202 + t88;
t71 = t292 * t74;
t29 = t288 * t70 + t71;
t26 = pkin(11) * t235 + t29;
t136 = -t259 * pkin(3) - t161;
t109 = -t202 * pkin(4) + qJD(5) + t136;
t338 = t292 * t202 - t203 * t288;
t60 = -pkin(5) * t338 - pkin(11) * t320 + t109;
t13 = -t26 * t296 + t299 * t60;
t35 = t106 * pkin(5) - t107 * pkin(11) + t96;
t305 = qJD(2) * t309;
t131 = qJD(1) * t305 + qJD(3) * t161;
t358 = qJD(2) * t291;
t346 = t289 * t358;
t335 = qJD(1) * t346;
t316 = t290 * t335;
t181 = pkin(3) * t226 - pkin(10) * t225 + t316;
t50 = -qJD(4) * t88 - t131 * t297 + t300 * t181;
t30 = pkin(4) * t226 - qJ(5) * t168 - qJD(5) * t203 + t50;
t354 = qJD(4) * t300;
t49 = t300 * t131 + t135 * t354 - t137 * t355 + t297 * t181;
t34 = qJ(5) * t169 + qJD(5) * t202 + t49;
t8 = t288 * t30 + t292 * t34;
t6 = pkin(11) * t226 + t8;
t1 = qJD(6) * t13 + t296 * t35 + t299 * t6;
t14 = t26 * t299 + t296 * t60;
t2 = -qJD(6) * t14 - t296 * t6 + t299 * t35;
t334 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t483 = t107 * t491;
t456 = -t226 * Ifges(6,6) / 0.2e1 + t483;
t440 = t106 / 0.2e1;
t457 = Ifges(6,2) * t440;
t489 = t334 + t96 * mrSges(6,1) + t15 / 0.2e1 + t456 + t457;
t487 = Ifges(4,4) * t496 + Ifges(4,5) * t495;
t486 = Ifges(4,6) * t495 - t240 * Ifges(4,4) / 0.2e1;
t485 = t240 * Ifges(4,1) / 0.2e1 - t487;
t484 = Ifges(4,2) * t496 + t486;
t431 = t168 / 0.2e1;
t430 = t169 / 0.2e1;
t478 = Ifges(4,5) * t225;
t477 = Ifges(6,5) * t107;
t475 = Ifges(4,6) * t226;
t474 = Ifges(6,6) * t106;
t473 = t131 * mrSges(4,2);
t287 = -pkin(4) * t300 - pkin(3);
t241 = pkin(5) * t272 - pkin(11) * t273 + t287;
t276 = t414 * t300;
t342 = t414 * t297;
t254 = -t292 * t276 + t288 * t342;
t195 = t241 * t299 - t254 * t296;
t471 = qJD(6) * t195 + t493 * t296 + t494 * t299;
t196 = t241 * t296 + t254 * t299;
t470 = -qJD(6) * t196 - t494 * t296 + t493 * t299;
t469 = pkin(5) * t240 - t490;
t119 = mrSges(6,1) * t235 - mrSges(6,3) * t320;
t75 = -mrSges(7,1) * t116 + mrSges(7,2) * t117;
t468 = t75 - t119;
t204 = -t248 * t290 + t294 * t258;
t152 = -pkin(3) * t304 - pkin(10) * t247 + t204;
t233 = t301 * t243;
t173 = t248 * t369 + t258 * t373 + t233;
t160 = pkin(10) * t265 + t173;
t98 = t297 * t152 + t300 * t160;
t228 = t270 * t288 + t271 * t292;
t315 = -t299 * t228 + t296 * t372;
t466 = qJD(6) * t315 - t463 * t296 + t492 * t299;
t210 = -t296 * t228 - t299 * t372;
t465 = qJD(6) * t210 + t492 * t296 + t463 * t299;
t464 = t462 * t288 - t461 * t292;
t362 = mrSges(4,1) * t259 + mrSges(5,1) * t202 - mrSges(5,2) * t203 - mrSges(4,3) * t240;
t460 = Ifges(5,5) * t168 + Ifges(5,6) * t169 + t479 * t226 - t474 + t477;
t459 = -t297 * t50 + t300 * t49;
t458 = t1 * t299 - t2 * t296;
t138 = qJD(6) - t338;
t390 = t203 * Ifges(5,4);
t121 = t202 * Ifges(5,2) + t235 * Ifges(5,6) + t390;
t201 = Ifges(5,4) * t202;
t122 = t203 * Ifges(5,1) + t235 * Ifges(5,5) + t201;
t321 = t297 * t88 + t300 * t87;
t409 = Ifges(5,4) * t300;
t410 = Ifges(5,4) * t297;
t419 = t300 / 0.2e1;
t424 = t235 / 0.2e1;
t427 = t203 / 0.2e1;
t428 = t202 / 0.2e1;
t455 = -t321 * mrSges(5,3) + (Ifges(5,5) * t300 - Ifges(5,6) * t297) * t424 + (-Ifges(5,2) * t297 + t409) * t428 + (Ifges(5,1) * t300 - t410) * t427 + t136 * (mrSges(5,1) * t297 + mrSges(5,2) * t300) - t297 * t121 / 0.2e1 + t122 * t419;
t381 = t288 * t74;
t28 = t292 * t70 - t381;
t25 = -pkin(5) * t235 - t28;
t324 = t13 * t299 + t14 * t296;
t326 = Ifges(7,5) * t299 - Ifges(7,6) * t296;
t406 = Ifges(7,4) * t299;
t328 = -Ifges(7,2) * t296 + t406;
t407 = Ifges(7,4) * t296;
t330 = Ifges(7,1) * t299 - t407;
t331 = mrSges(7,1) * t296 + mrSges(7,2) * t299;
t420 = t299 / 0.2e1;
t421 = -t296 / 0.2e1;
t434 = t138 / 0.2e1;
t436 = t117 / 0.2e1;
t438 = t116 / 0.2e1;
t408 = Ifges(7,4) * t117;
t54 = Ifges(7,2) * t116 + Ifges(7,6) * t138 + t408;
t114 = Ifges(7,4) * t116;
t55 = Ifges(7,1) * t117 + Ifges(7,5) * t138 + t114;
t454 = -t324 * mrSges(7,3) + t25 * t331 + t326 * t434 + t328 * t438 + t330 * t436 + t420 * t55 + t421 * t54;
t386 = t235 * Ifges(6,5);
t396 = t320 * Ifges(6,1);
t399 = t338 * Ifges(6,4);
t86 = t386 + t396 + t399;
t453 = t109 * mrSges(6,2) + t86 / 0.2e1;
t451 = t194 * mrSges(4,2) - t161 * mrSges(4,3) + t485;
t7 = -t288 * t34 + t292 * t30;
t450 = t50 * mrSges(5,1) + t7 * mrSges(6,1) - t49 * mrSges(5,2) - t8 * mrSges(6,2);
t400 = t138 * Ifges(7,3);
t401 = t117 * Ifges(7,5);
t402 = t116 * Ifges(7,6);
t53 = t400 + t401 + t402;
t385 = t235 * Ifges(6,6);
t395 = t320 * Ifges(6,4);
t398 = t338 * Ifges(6,2);
t85 = t385 + t395 + t398;
t449 = t109 * mrSges(6,1) + t13 * mrSges(7,1) + t53 / 0.2e1 - t85 / 0.2e1 - t14 * mrSges(7,2);
t448 = t194 * mrSges(4,1) + t87 * mrSges(5,1) + t28 * mrSges(6,1) - t88 * mrSges(5,2) - t29 * mrSges(6,2) - t162 * mrSges(4,3) + t484;
t17 = Ifges(7,1) * t64 + Ifges(7,4) * t65 + Ifges(7,5) * t106;
t446 = t17 / 0.2e1;
t443 = t64 / 0.2e1;
t442 = t65 / 0.2e1;
t441 = Ifges(5,1) * t431 + Ifges(5,4) * t430 + Ifges(5,5) * t426;
t439 = -t116 / 0.2e1;
t437 = -t117 / 0.2e1;
t435 = -t138 / 0.2e1;
t433 = t338 / 0.2e1;
t432 = t320 / 0.2e1;
t425 = -t235 / 0.2e1;
t417 = pkin(4) * t203;
t22 = -mrSges(7,1) * t65 + mrSges(7,2) * t64;
t93 = mrSges(6,1) * t226 - mrSges(6,3) * t107;
t413 = t22 - t93;
t208 = -t247 * t297 + t265 * t300;
t180 = qJD(4) * t208 + t238 * t300;
t209 = t247 * t300 + t265 * t297;
t146 = qJD(3) * t172 + t305;
t336 = t290 * t346;
t186 = pkin(3) * t239 - pkin(10) * t238 + t336;
t59 = -qJD(4) * t98 - t146 * t297 + t300 * t186;
t41 = pkin(4) * t239 - qJ(5) * t180 - qJD(5) * t209 + t59;
t179 = -qJD(4) * t209 - t238 * t297;
t58 = t300 * t146 + t152 * t354 - t160 * t355 + t297 * t186;
t45 = qJ(5) * t179 + qJD(5) * t208 + t58;
t12 = t288 * t41 + t292 * t45;
t97 = t300 * t152 - t160 * t297;
t77 = -pkin(4) * t304 - qJ(5) * t209 + t97;
t83 = qJ(5) * t208 + t98;
t43 = t288 * t77 + t292 * t83;
t412 = mrSges(4,3) * t225;
t411 = mrSges(4,3) * t226;
t397 = t338 * Ifges(6,6);
t394 = t320 * Ifges(6,5);
t391 = t202 * Ifges(5,6);
t389 = t203 * Ifges(5,5);
t377 = t132 * t301;
t370 = t293 * (-mrSges(3,2) * t295 + mrSges(3,3) * t371) * qJD(1);
t367 = t296 * t267;
t366 = t299 * t267;
t361 = -t475 + t478;
t91 = -mrSges(6,1) * t338 + mrSges(6,2) * t320;
t352 = t91 - t362;
t61 = t106 * mrSges(6,1) + t107 * mrSges(6,2);
t153 = t188 * t296 + t240 * t299;
t340 = t153 - t367;
t154 = -t188 * t299 + t240 * t296;
t339 = -t154 - t366;
t333 = -t1 * t296 - t2 * t299;
t332 = mrSges(7,1) * t299 - mrSges(7,2) * t296;
t329 = Ifges(7,1) * t296 + t406;
t327 = Ifges(7,2) * t299 + t407;
t325 = Ifges(7,5) * t296 + Ifges(7,6) * t299;
t323 = t13 * t296 - t14 * t299;
t11 = -t288 * t45 + t292 * t41;
t42 = -t288 * t83 + t292 * t77;
t37 = -pkin(11) * t304 + t43;
t159 = -t265 * pkin(3) - t172;
t115 = -t208 * pkin(4) + t159;
t163 = -t292 * t208 + t209 * t288;
t164 = t208 * t288 + t209 * t292;
t66 = t163 * pkin(5) - t164 * pkin(11) + t115;
t19 = t296 * t66 + t299 * t37;
t18 = -t296 * t37 + t299 * t66;
t80 = -mrSges(7,2) * t138 + mrSges(7,3) * t116;
t81 = mrSges(7,1) * t138 - mrSges(7,3) * t117;
t322 = -t296 * t81 + t299 * t80;
t127 = t164 * t299 - t296 * t304;
t126 = -t164 * t296 - t299 * t304;
t317 = -(-qJ(2) * t347 + t280) * t289 + t263 * t293;
t268 = (mrSges(3,1) * t295 - mrSges(3,3) * t375) * qJD(1);
t147 = t306 + (t298 * t318 + t233) * qJD(3);
t101 = -t179 * pkin(4) + t147;
t286 = -pkin(4) * t292 - pkin(5);
t253 = -t276 * t288 - t292 * t342;
t227 = -t292 * t270 + t271 * t288;
t205 = -mrSges(4,2) * t259 + mrSges(4,3) * t237;
t192 = -mrSges(4,1) * t237 + mrSges(4,2) * t240;
t189 = mrSges(4,1) * t226 + mrSges(4,2) * t225;
t171 = mrSges(5,1) * t235 - mrSges(5,3) * t203;
t170 = -mrSges(5,2) * t235 + mrSges(5,3) * t202;
t130 = -mrSges(5,2) * t226 + mrSges(5,3) * t169;
t129 = mrSges(5,1) * t226 - mrSges(5,3) * t168;
t120 = t235 * Ifges(5,3) + t389 + t391;
t118 = -mrSges(6,2) * t235 + mrSges(6,3) * t338;
t113 = t179 * t288 + t180 * t292;
t112 = -t292 * t179 + t180 * t288;
t108 = -mrSges(5,1) * t169 + mrSges(5,2) * t168;
t94 = t168 * Ifges(5,4) + t169 * Ifges(5,2) + t226 * Ifges(5,6);
t92 = -mrSges(6,2) * t226 - mrSges(6,3) * t106;
t84 = t235 * Ifges(6,3) + t394 + t397;
t78 = pkin(5) * t320 - pkin(11) * t338 + t417;
t69 = -qJD(6) * t127 - t113 * t296 + t239 * t299;
t68 = qJD(6) * t126 + t113 * t299 + t239 * t296;
t46 = t112 * pkin(5) - t113 * pkin(11) + t101;
t39 = -mrSges(7,2) * t106 + mrSges(7,3) * t65;
t38 = mrSges(7,1) * t106 - mrSges(7,3) * t64;
t36 = pkin(5) * t304 - t42;
t33 = t292 * t73 - t381;
t32 = t288 * t73 + t71;
t21 = t296 * t78 + t299 * t33;
t20 = -t296 * t33 + t299 * t78;
t16 = t64 * Ifges(7,4) + t65 * Ifges(7,2) + t106 * Ifges(7,6);
t10 = pkin(11) * t239 + t12;
t9 = -pkin(5) * t239 - t11;
t5 = -pkin(5) * t226 - t7;
t4 = -qJD(6) * t19 - t10 * t296 + t299 * t46;
t3 = qJD(6) * t18 + t10 * t299 + t296 * t46;
t23 = [(m(3) * ((t293 * t360 + (qJ(2) * t375 - t283) * t289) * qJD(1) + t317) + 0.2e1 * t370) * t358 + (mrSges(4,2) * t316 + mrSges(4,3) * t132 + Ifges(4,1) * t225 - Ifges(4,4) * t226) * t247 + (t120 + t84) * t239 / 0.2e1 + m(4) * (t131 * t173 - t132 * t172 + t146 * t162 + (qJD(1) * t204 + t194) * t336) + (Ifges(5,5) * t427 + Ifges(6,5) * t432 + Ifges(5,6) * t428 + Ifges(6,6) * t433 + t479 * t424 + t448 + t484) * t239 + (t451 + t485) * t238 + m(5) * (t132 * t159 + t49 * t98 + t50 * t97 + t58 * t88 + t59 * t87) + (Ifges(7,5) * t68 + Ifges(7,6) * t69) * t434 + (Ifges(7,5) * t127 + Ifges(7,6) * t126) * t440 + (Ifges(7,5) * t443 - Ifges(6,6) * t426 + Ifges(7,6) * t442 + Ifges(7,3) * t440 + t457 + t483 + t489) * t163 + (Ifges(5,5) * t180 + Ifges(5,6) * t179) * t424 + (Ifges(5,5) * t209 + Ifges(5,6) * t208) * t426 + (Ifges(7,4) * t68 + Ifges(7,2) * t69) * t438 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t442 + (Ifges(5,4) * t180 + Ifges(5,2) * t179) * t428 + (Ifges(5,4) * t209 + Ifges(5,2) * t208) * t430 + (Ifges(7,1) * t68 + Ifges(7,4) * t69) * t436 + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t443 + (t1 * t126 - t127 * t2 - t13 * t68 + t14 * t69) * mrSges(7,3) + (Ifges(6,1) * t432 + Ifges(6,4) * t433 + Ifges(6,5) * t424 + t453) * t113 + t498 * t164 + m(6) * (t101 * t109 + t11 * t28 + t115 * t96 + t12 * t29 + t42 * t7 + t43 * t8) + m(7) * (t1 * t19 + t13 * t4 + t14 * t3 + t18 * t2 + t25 * t9 + t36 * t5) + (Ifges(5,1) * t180 + Ifges(5,4) * t179) * t427 + (Ifges(5,1) * t209 + Ifges(5,4) * t208) * t431 + (-Ifges(6,4) * t432 + Ifges(7,5) * t436 - Ifges(6,2) * t433 - Ifges(6,6) * t424 + Ifges(7,6) * t438 + Ifges(7,3) * t434 + t449) * t112 + (t179 * t88 - t180 * t87 + t208 * t49 - t209 * t50) * mrSges(5,3) + (-t112 * t29 - t113 * t28 - t163 * t8 - t164 * t7) * mrSges(6,3) + t132 * (-mrSges(5,1) * t208 + mrSges(5,2) * t209) + t208 * t94 / 0.2e1 + t204 * t189 + t146 * t205 + t179 * t121 / 0.2e1 + t136 * (-mrSges(5,1) * t179 + mrSges(5,2) * t180) + t180 * t122 / 0.2e1 + t58 * t170 + t59 * t171 + t159 * t108 + t126 * t16 / 0.2e1 + t5 * (-mrSges(7,1) * t126 + mrSges(7,2) * t127) + t97 * t129 + t98 * t130 + t12 * t118 + t11 * t119 + t115 * t61 + t101 * t91 + t43 * t92 + t42 * t93 + t3 * t80 + t4 * t81 + t9 * t75 + t25 * (-mrSges(7,1) * t69 + mrSges(7,2) * t68) + t69 * t54 / 0.2e1 + t68 * t55 / 0.2e1 + t192 * t336 + (-m(4) * t161 + m(5) * t136 - t362) * t147 - (t460 / 0.2e1 - t131 * mrSges(4,3) - Ifges(4,4) * t225 + Ifges(4,2) * t226 - t474 / 0.2e1 + t477 / 0.2e1 + mrSges(4,1) * t316 + Ifges(5,6) * t430 + Ifges(5,5) * t431 + t479 * t426 + t450) * t304 - 0.2e1 * t268 * t346 + t209 * t441 + t127 * t446 + (-t132 * mrSges(4,1) - t473 + t478 / 0.2e1 - t475 / 0.2e1 + t361 / 0.2e1) * t265 - t173 * t411 - t172 * t412 + t36 * t22 + t18 * t38 + t19 * t39; t466 * t81 + t465 * t80 + t462 * t170 + t461 * t171 + t463 * t118 - m(4) * (-t161 * t256 + t162 * t257) + t294 * t189 + t270 * t129 + t271 * t130 - t257 * t205 + t228 * t92 + t210 * t38 - t315 * t39 - t352 * t256 + t413 * t227 + t468 * t464 + (-m(3) * t317 + t289 * t268 - t370) * t359 + (-t1 * t315 + t13 * t466 + t14 * t465 + t2 * t210 + t227 * t5 + t25 * t464) * m(7) + (-t109 * t256 - t7 * t227 + t8 * t228 - t28 * t464 + t29 * t463) * m(6) + (-t136 * t256 + t50 * t270 + t49 * t271 + t461 * t87 + t462 * t88) * m(5) + (m(6) * (t109 * t357 - t301 * t96) + m(5) * (t136 * t357 - t377) - t192 * t347 - t298 * t411 + (-t108 - t61 - t412) * t301 + (t205 * t301 + t298 * t352) * qJD(3) + (t131 * t298 - t161 * t357 + t162 * t356 - t194 * t347 + t294 * t335 - t377) * m(4)) * t290; ((-m(5) * t321 - t297 * t170 - t300 * t171) * qJD(4) - t129 * t297 + t130 * t300 + m(5) * t459) * pkin(10) + t459 * mrSges(5,3) + ((m(6) * t109 + t91) * t297 * pkin(4) + t455) * qJD(4) + (-Ifges(6,4) * t267 - Ifges(6,2) * t266) * t433 + (-t154 / 0.2e1 - t366 / 0.2e1) * t55 + (-Ifges(7,5) * t366 + Ifges(7,6) * t367 + Ifges(7,3) * t266) * t434 + (-Ifges(7,1) * t366 + Ifges(7,4) * t367 + Ifges(7,5) * t266) * t436 + (-Ifges(7,4) * t366 + Ifges(7,2) * t367 + Ifges(7,6) * t266) * t438 + (-t153 / 0.2e1 + t367 / 0.2e1) * t54 + (t53 - t85) * (t266 / 0.2e1 - t187 / 0.2e1) - t338 * (-Ifges(6,4) * t188 - Ifges(6,2) * t187) / 0.2e1 - t320 * (-Ifges(6,1) * t188 - Ifges(6,4) * t187) / 0.2e1 + (-Ifges(6,5) * t188 - Ifges(6,6) * t187) * t425 + (-t267 / 0.2e1 + t188 / 0.2e1) * t86 + (-Ifges(6,5) * t267 - Ifges(6,6) * t266) * t424 + (-Ifges(6,1) * t267 - Ifges(6,4) * t266) * t432 + t490 * t119 + (-t8 * mrSges(6,3) + t62 / 0.2e1 + t63 / 0.2e1 + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t106 + t456 + t489) * t272 + (-t391 / 0.2e1 - t389 / 0.2e1 - t120 / 0.2e1 - t84 / 0.2e1 - t394 / 0.2e1 - t397 / 0.2e1 + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1) * t235 - t448 - t486) * t240 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t240 - t451 - t455 + t487) * t237 + t361 + t499 * t118 + (t326 * t440 + t330 * t443 + t328 * t442 + t5 * t331 - t7 * mrSges(6,3) + t16 * t421 + t17 * t420 + t333 * mrSges(7,3) + (t325 * t435 + t327 * t439 + t329 * t437 + t25 * t332 + t55 * t421 - t299 * t54 / 0.2e1 + t323 * mrSges(7,3)) * qJD(6) + t498) * t273 - t473 + (-pkin(3) * t132 - t110 * t87 - t111 * t88 - t136 * t162) * m(5) - m(6) * (t109 * t125 + t28 * t51 + t29 * t52) + t287 * t61 + t254 * t92 - t161 * t205 + t195 * t38 + t196 * t39 - t111 * t170 - t110 * t171 - t125 * t91 - pkin(3) * t108 + m(6) * (-t218 * t28 + t219 * t29 - t253 * t7 + t254 * t8 + t287 * t96) + (-mrSges(5,1) * t300 + mrSges(5,2) * t297 - mrSges(4,1)) * t132 + (mrSges(7,1) * t340 + mrSges(7,2) * t339) * t25 + t469 * t75 + t470 * t81 + t471 * t80 + (t1 * t196 + t13 * t470 + t14 * t471 + t195 * t2 + t25 * t469 + t253 * t5) * m(7) + t94 * t419 + (Ifges(5,5) * t297 + Ifges(5,6) * t300) * t426 + (Ifges(5,2) * t300 + t410) * t430 + (Ifges(5,1) * t297 + t409) * t431 + (Ifges(7,5) * t154 + Ifges(7,6) * t153 + Ifges(7,3) * t187) * t435 + (Ifges(7,1) * t154 + Ifges(7,4) * t153 + Ifges(7,5) * t187) * t437 + (Ifges(7,4) * t154 + Ifges(7,2) * t153 + Ifges(7,6) * t187) * t439 + t297 * t441 + t362 * t162 + (mrSges(7,2) * t365 - mrSges(7,3) * t340) * t14 + (t28 * t364 + t29 * t365) * mrSges(6,3) + (-mrSges(7,1) * t365 - mrSges(7,3) * t339) * t13 + (-mrSges(6,1) * t365 - mrSges(6,2) * t364) * t109 + t413 * t253; -t468 * t32 + ((-m(7) * t324 - t296 * t80 - t299 * t81) * qJD(6) - t296 * t38 + t299 * t39 + m(7) * t458) * (pkin(4) * t288 + pkin(11)) + t458 * mrSges(7,3) + (t28 * mrSges(6,3) - t386 / 0.2e1 - t399 / 0.2e1 - t396 / 0.2e1 - t453 - t454) * t338 + t454 * qJD(6) - (-Ifges(5,2) * t203 + t122 + t201) * t202 / 0.2e1 + t460 + ((t288 * t8 + t292 * t7) * pkin(4) - t109 * t417 + t28 * t32 - t29 * t33) * m(6) + (t202 * t87 + t203 * t88) * mrSges(5,3) + t450 - t5 * t332 + (t29 * mrSges(6,3) + t385 / 0.2e1 + t398 / 0.2e1 + t395 / 0.2e1 - t402 / 0.2e1 - t401 / 0.2e1 - t400 / 0.2e1 - t449) * t320 + t286 * t22 + (-t13 * t20 - t14 * t21 - t25 * t32 + t286 * t5) * m(7) - t136 * (mrSges(5,1) * t203 + mrSges(5,2) * t202) - t87 * t170 + t88 * t171 - t33 * t118 - t21 * t80 - t20 * t81 + t16 * t420 + (Ifges(5,5) * t202 - Ifges(5,6) * t203) * t425 + t121 * t427 + t325 * t440 + t327 * t442 + t329 * t443 + t296 * t446 - t203 * (Ifges(5,1) * t202 - t390) / 0.2e1 + (-t203 * t91 + t288 * t92 + t292 * t93) * pkin(4); t296 * t39 + t299 * t38 - t468 * t320 + t322 * qJD(6) + (-t118 - t322) * t338 + t61 + (-t138 * t323 - t320 * t25 - t333) * m(7) + (t28 * t320 - t29 * t338 + t96) * m(6); -t13 * t80 + t14 * t81 - t25 * (mrSges(7,1) * t117 + mrSges(7,2) * t116) + (Ifges(7,1) * t116 - t408) * t437 + t54 * t436 + (Ifges(7,5) * t116 - Ifges(7,6) * t117) * t435 + (t116 * t13 + t117 * t14) * mrSges(7,3) + t334 + t15 + (-Ifges(7,2) * t117 + t114 + t55) * t439;];
tauc  = t23(:);
