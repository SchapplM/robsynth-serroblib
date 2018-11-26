% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2018-11-23 17:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:38:30
% EndTime: 2018-11-23 17:38:52
% DurationCPUTime: 22.52s
% Computational Cost: add. (8521->760), mult. (22263->971), div. (0->0), fcn. (16044->8), ass. (0->324)
t264 = cos(qJ(2));
t257 = sin(pkin(6));
t349 = qJD(1) * t257;
t329 = t264 * t349;
t231 = -qJD(3) + t329;
t391 = -t231 / 0.2e1;
t260 = sin(qJ(3));
t263 = cos(qJ(3));
t364 = cos(pkin(6));
t285 = qJD(1) * t364 + qJD(2);
t261 = sin(qJ(2));
t330 = t261 * t349;
t186 = t260 * t285 + t263 * t330;
t394 = t186 / 0.2e1;
t432 = Ifges(5,4) + Ifges(4,5);
t435 = Ifges(4,1) + Ifges(5,1);
t180 = qJD(6) + t186;
t398 = t180 / 0.2e1;
t279 = t263 * t285;
t185 = t260 * t330 - t279;
t259 = sin(qJ(6));
t262 = cos(qJ(6));
t124 = t185 * t262 + t231 * t259;
t404 = t124 / 0.2e1;
t123 = -t185 * t259 + t231 * t262;
t406 = t123 / 0.2e1;
t464 = Ifges(7,5) * t404 + Ifges(7,6) * t406 + Ifges(7,3) * t398;
t475 = t432 * t391 + t435 * t394 + t464;
t397 = -t185 / 0.2e1;
t474 = Ifges(4,4) * t397;
t332 = pkin(1) * t364;
t246 = t264 * t332;
t202 = -pkin(8) * t330 + qJD(1) * t246;
t282 = (pkin(2) * t261 - pkin(9) * t264) * t257;
t203 = qJD(1) * t282;
t117 = -t260 * t202 + t203 * t263;
t345 = qJD(3) * t263;
t381 = pkin(9) - qJ(5);
t473 = -qJD(5) * t260 + t345 * t381 + t117;
t222 = t263 * qJ(4) * t329;
t350 = qJ(4) * t345 + t260 * qJD(4);
t472 = t350 - t222;
t396 = t185 / 0.2e1;
t395 = -t186 / 0.2e1;
t434 = Ifges(6,1) + Ifges(5,3);
t245 = t261 * t332;
t409 = -pkin(4) - pkin(10);
t343 = pkin(3) - t409;
t280 = pkin(5) * t263 - t260 * t343;
t360 = t257 * t264;
t470 = qJD(3) * t280 - (-t245 + (-pkin(8) + t280) * t360) * qJD(1) + t472;
t317 = t343 * t261;
t355 = t263 * qJ(5);
t338 = t264 * t355;
t469 = -(-t317 - t338) * t349 + t473;
t468 = Ifges(5,5) * t395;
t390 = t231 / 0.2e1;
t467 = -qJD(1) / 0.2e1;
t383 = -qJD(2) / 0.2e1;
t118 = t263 * t202 + t260 * t203;
t102 = qJ(4) * t330 + t118;
t346 = qJD(3) * t260;
t466 = -qJD(5) * t263 - t346 * t381 - t102;
t465 = Ifges(5,5) * t396 + t474;
t322 = Ifges(3,6) * t364;
t380 = Ifges(3,4) * t261;
t463 = (t322 + (Ifges(3,2) * t264 + t380) * t257) * t467;
t462 = Ifges(3,6) * t383;
t429 = Ifges(5,2) + Ifges(4,3);
t428 = Ifges(4,6) - Ifges(5,6);
t224 = -t263 * pkin(3) - t260 * qJ(4) - pkin(2);
t217 = t263 * pkin(4) - t224;
t179 = pkin(5) * t260 + pkin(10) * t263 + t217;
t232 = t381 * t260;
t122 = t179 * t259 + t232 * t262;
t461 = -qJD(6) * t122 - t259 * t469 + t262 * t470;
t121 = t179 * t262 - t232 * t259;
t460 = qJD(6) * t121 + t259 * t470 + t262 * t469;
t357 = t260 * t264;
t459 = -(pkin(5) * t261 + qJ(5) * t357) * t349 + t466;
t410 = pkin(3) + pkin(4);
t458 = -(-t261 * t410 - t338) * t349 + t473;
t457 = qJ(5) * t260 * t329 - t466;
t339 = t410 * t260;
t456 = -(-t245 + (-pkin(8) - t339) * t360) * qJD(1) - qJD(3) * t339 + t472;
t221 = t231 * qJ(4);
t351 = pkin(8) * t360 + t245;
t197 = t364 * pkin(9) + t351;
t160 = qJD(2) * pkin(9) + qJD(1) * t197;
t198 = (-pkin(2) * t264 - pkin(9) * t261 - pkin(1)) * t257;
t169 = qJD(1) * t198;
t98 = t263 * t160 + t260 * t169;
t84 = -t221 + t98;
t455 = -t84 * mrSges(5,2) - t98 * mrSges(4,3);
t352 = t260 * t160 - t263 * t169;
t324 = qJD(4) + t352;
t82 = pkin(3) * t231 + t324;
t175 = Ifges(6,4) * t185;
t91 = -t186 * Ifges(6,2) + t231 * Ifges(6,6) + t175;
t454 = t91 / 0.2e1 - t352 * mrSges(4,3) - t82 * mrSges(5,2);
t159 = -t285 * pkin(2) - t202;
t172 = t186 * qJ(5);
t315 = -t172 + t324;
t44 = t231 * t410 + t315;
t77 = t185 * pkin(3) - t186 * qJ(4) + t159;
t281 = qJD(5) - t77;
t49 = -pkin(4) * t185 + t281;
t28 = pkin(5) * t186 + t185 * t409 + t281;
t35 = t231 * t343 + t315;
t5 = -t259 * t35 + t262 * t28;
t6 = t259 * t28 + t262 * t35;
t453 = t49 * mrSges(6,1) + t5 * mrSges(7,1) + t159 * mrSges(4,2) - t6 * mrSges(7,2) - t77 * mrSges(5,3) - t44 * mrSges(6,3) + t465 + t475;
t67 = qJ(5) * t185 + t98;
t53 = t221 - t67;
t452 = -t159 * mrSges(4,1) - t77 * mrSges(5,1) - t49 * mrSges(6,2) + t53 * mrSges(6,3) + Ifges(5,6) * t390 + t468 + (Ifges(4,4) + Ifges(6,4)) * t394 + (Ifges(6,5) + Ifges(4,6)) * t391 + (Ifges(4,2) + t434) * t397;
t342 = qJD(1) * qJD(2);
t321 = t264 * t342;
t313 = t257 * t321;
t361 = t257 * t261;
t337 = t260 * t361;
t316 = qJD(3) * t337;
t135 = qJD(1) * t316 - qJD(3) * t279 - t263 * t313;
t204 = qJD(2) * t282;
t193 = qJD(1) * t204;
t215 = -pkin(8) * t361 + t246;
t206 = t215 * qJD(2);
t194 = qJD(1) * t206;
t38 = -t160 * t345 - t169 * t346 + t193 * t263 - t260 * t194;
t274 = qJ(5) * t135 - qJD(5) * t186 - t38;
t348 = qJD(2) * t257;
t283 = t317 * t348;
t16 = -qJD(1) * t283 + t274;
t136 = qJD(3) * t186 + t260 * t313;
t422 = t351 * qJD(2);
t195 = qJD(1) * t422;
t33 = t136 * pkin(3) + t135 * qJ(4) - t186 * qJD(4) + t195;
t7 = -t135 * pkin(5) + t136 * t409 - t33;
t1 = qJD(6) * t5 + t16 * t262 + t259 * t7;
t451 = t1 * mrSges(7,2);
t2 = -qJD(6) * t6 - t16 * t259 + t262 * t7;
t450 = t2 * mrSges(7,1);
t318 = t364 * t263;
t327 = t264 * t348;
t149 = -qJD(3) * t318 - t263 * t327 + t316;
t214 = t260 * t364 + t263 * t361;
t150 = qJD(3) * t214 + t260 * t327;
t50 = t150 * pkin(3) + t149 * qJ(4) - t214 * qJD(4) + t422;
t205 = t351 * qJD(1);
t443 = -t352 * mrSges(4,1) - t82 * mrSges(5,1) - t53 * mrSges(6,1) - t98 * mrSges(4,2) + t44 * mrSges(6,2) - t205 * mrSges(3,3) + t84 * mrSges(5,3) + Ifges(6,5) * t397 + Ifges(6,6) * t394 + Ifges(6,3) * t391 + t462 + t463;
t433 = -Ifges(4,4) + Ifges(5,5);
t442 = Ifges(6,4) * t395 + Ifges(6,5) * t390 - Ifges(4,2) * t397 - t428 * t391 + t433 * t394 + t434 * t396 - t452;
t431 = Ifges(6,4) - Ifges(5,5);
t441 = -Ifges(6,2) * t395 - Ifges(6,6) * t390 - t431 * t396 + t453 + t474 + t475;
t320 = t261 * t342;
t314 = t257 * t320;
t54 = qJD(6) * t123 + t136 * t262 - t259 * t314;
t440 = -t54 / 0.2e1;
t55 = -qJD(6) * t124 - t136 * t259 - t262 * t314;
t439 = -t55 / 0.2e1;
t438 = -t135 / 0.2e1;
t437 = t135 / 0.2e1;
t436 = -t136 / 0.2e1;
t430 = Ifges(6,5) - Ifges(5,6);
t426 = t194 * mrSges(3,2);
t66 = -t172 + t352;
t425 = qJD(4) + t66;
t108 = -mrSges(5,1) * t314 - t135 * mrSges(5,2);
t34 = -pkin(3) * t314 - t38;
t424 = m(5) * t34 + t108;
t423 = mrSges(3,1) * t285 - mrSges(4,1) * t185 - mrSges(4,2) * t186 - mrSges(3,3) * t330;
t303 = -t1 * t262 + t2 * t259;
t420 = (Ifges(6,3) + t429) * t314 + (-Ifges(6,5) - t428) * t136 + (-Ifges(6,6) - t432) * t135;
t307 = -Ifges(6,6) / 0.2e1 - Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t419 = t307 * t231 + t453 - t454 + t464;
t376 = Ifges(7,4) * t262;
t295 = -Ifges(7,2) * t259 + t376;
t377 = Ifges(7,4) * t259;
t297 = Ifges(7,1) * t262 - t377;
t298 = mrSges(7,1) * t259 + mrSges(7,2) * t262;
t301 = t259 * t6 + t262 * t5;
t372 = Ifges(7,6) * t259;
t374 = Ifges(7,5) * t262;
t387 = -t262 / 0.2e1;
t389 = t259 / 0.2e1;
t399 = -t180 / 0.2e1;
t405 = -t124 / 0.2e1;
t407 = -t123 / 0.2e1;
t378 = Ifges(7,4) * t124;
t42 = Ifges(7,2) * t123 + Ifges(7,6) * t180 + t378;
t120 = Ifges(7,4) * t123;
t43 = Ifges(7,1) * t124 + Ifges(7,5) * t180 + t120;
t45 = -pkin(5) * t231 - t53;
t418 = mrSges(7,3) * t301 + (-t372 + t374) * t399 + t295 * t407 + t297 * t405 - t45 * t298 + t387 * t43 + t389 * t42;
t416 = -t42 / 0.2e1;
t415 = t42 / 0.2e1;
t414 = t43 / 0.2e1;
t413 = Ifges(6,4) * t436 + Ifges(6,2) * t438 + Ifges(6,6) * t314 / 0.2e1;
t52 = Ifges(7,5) * t54;
t51 = Ifges(7,6) * t55;
t258 = qJ(4) + pkin(5);
t379 = Ifges(3,4) * t264;
t373 = Ifges(3,2) * t261;
t368 = t202 * mrSges(3,3);
t366 = t261 * Ifges(3,1);
t137 = mrSges(6,1) * t231 - mrSges(6,3) * t185;
t65 = -mrSges(7,1) * t123 + mrSges(7,2) * t124;
t365 = t65 - t137;
t213 = -t318 + t337;
t362 = qJ(5) * t213;
t359 = t259 * t260;
t358 = t260 * t262;
t356 = t262 * t264;
t106 = mrSges(6,2) * t314 + t135 * mrSges(6,3);
t138 = mrSges(4,2) * t231 - mrSges(4,3) * t185;
t142 = -mrSges(5,2) * t185 - mrSges(5,3) * t231;
t354 = -t138 - t142;
t139 = -mrSges(4,1) * t231 - mrSges(4,3) * t186;
t140 = mrSges(5,1) * t231 + mrSges(5,2) * t186;
t353 = t139 - t140;
t112 = t186 * pkin(3) + t185 * qJ(4);
t116 = t263 * t197 + t260 * t198;
t347 = qJD(2) * t261;
t344 = qJD(4) * t264;
t8 = -Ifges(7,3) * t135 + t51 + t52;
t341 = Ifges(6,4) / 0.2e1 + Ifges(4,4) / 0.2e1;
t340 = t142 + t365;
t196 = -t364 * pkin(2) - t215;
t328 = t257 * t347;
t326 = t257 * t344;
t323 = Ifges(3,5) * t364;
t73 = -t135 * mrSges(6,1) + t136 * mrSges(6,2);
t115 = -t260 * t197 + t198 * t263;
t310 = -Ifges(4,1) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(6,2) / 0.2e1;
t309 = Ifges(5,5) / 0.2e1 - t341;
t308 = -Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1 + Ifges(4,6) / 0.2e1;
t306 = -Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1 - Ifges(6,1) / 0.2e1;
t305 = t450 - t451;
t302 = t1 * t259 + t2 * t262;
t300 = t259 * t5 - t262 * t6;
t101 = pkin(3) * t360 - t115;
t299 = mrSges(7,1) * t262 - mrSges(7,2) * t259;
t296 = Ifges(7,1) * t259 + t376;
t294 = Ifges(7,2) * t262 + t377;
t292 = Ifges(7,5) * t259 + Ifges(7,6) * t262;
t29 = -mrSges(7,1) * t135 - mrSges(7,3) * t54;
t30 = mrSges(7,2) * t135 + mrSges(7,3) * t55;
t291 = -t259 * t29 + t262 * t30;
t99 = t213 * pkin(3) - t214 * qJ(4) + t196;
t40 = pkin(5) * t214 + t213 * t409 - t99;
t70 = pkin(4) * t360 - qJ(5) * t214 + t101;
t64 = pkin(10) * t360 + t70;
t14 = -t259 * t64 + t262 * t40;
t15 = t259 * t40 + t262 * t64;
t85 = -mrSges(7,2) * t180 + mrSges(7,3) * t123;
t86 = mrSges(7,1) * t180 - mrSges(7,3) * t124;
t290 = -t259 * t86 + t262 * t85;
t289 = -t259 * t85 - t262 * t86;
t288 = -t260 * t84 + t263 * t82;
t287 = -t260 * t98 + t263 * t352;
t286 = t410 * t328;
t100 = -qJ(4) * t360 + t116;
t57 = -t197 * t345 - t198 * t346 + t204 * t263 - t260 * t206;
t151 = -t213 * t259 + t257 * t356;
t152 = t213 * t262 + t259 * t360;
t37 = -t160 * t346 + t169 * t345 + t260 * t193 + t263 * t194;
t56 = -t197 * t346 + t198 * t345 + t260 * t204 + t263 * t206;
t278 = qJ(4) * t328 + t56;
t31 = qJ(4) * t314 - t231 * qJD(4) + t37;
t273 = qJ(5) * t149 - qJD(5) * t214 - t57;
t270 = qJ(5) * t150 + qJD(5) * t213 + t278;
t17 = -qJ(5) * t136 - qJD(5) * t185 - t31;
t267 = -t308 * t231 + t452 - t455;
t236 = Ifges(3,4) * t329;
t233 = pkin(9) * t263 - t355;
t230 = Ifges(3,5) * t313;
t211 = pkin(3) * t346 - t350;
t200 = -mrSges(3,2) * t285 + mrSges(3,3) * t329;
t165 = (-t259 * t261 + t260 * t356) * t349;
t164 = (-t259 * t357 - t261 * t262) * t349;
t156 = Ifges(3,1) * t330 + Ifges(3,5) * t285 + t236;
t141 = -mrSges(6,2) * t231 - mrSges(6,3) * t186;
t119 = -t222 + (t245 + (pkin(3) * t260 + pkin(8)) * t360) * qJD(1);
t113 = mrSges(5,1) * t185 - mrSges(5,3) * t186;
t111 = mrSges(6,1) * t186 + mrSges(6,2) * t185;
t110 = -mrSges(4,2) * t314 - mrSges(4,3) * t136;
t109 = -mrSges(6,1) * t314 - mrSges(6,3) * t136;
t107 = mrSges(4,1) * t314 + mrSges(4,3) * t135;
t105 = -mrSges(5,2) * t136 + mrSges(5,3) * t314;
t104 = -pkin(3) * t330 - t117;
t92 = t186 * Ifges(5,4) - t231 * Ifges(5,2) + t185 * Ifges(5,6);
t90 = t186 * Ifges(4,5) - t185 * Ifges(4,6) - t231 * Ifges(4,3);
t83 = -pkin(4) * t186 - t112;
t80 = -t100 - t362;
t78 = -pkin(4) * t213 - t99;
t75 = -qJD(6) * t152 - t150 * t259 - t262 * t328;
t74 = qJD(6) * t151 + t150 * t262 - t259 * t328;
t72 = mrSges(4,1) * t136 - mrSges(4,2) * t135;
t71 = mrSges(5,1) * t136 + mrSges(5,3) * t135;
t69 = -t258 * t360 + t116 + t362;
t63 = -t135 * Ifges(4,1) - t136 * Ifges(4,4) + Ifges(4,5) * t314;
t62 = -t135 * Ifges(5,1) + Ifges(5,4) * t314 + t136 * Ifges(5,5);
t61 = t136 * Ifges(6,1) + t135 * Ifges(6,4) - Ifges(6,5) * t314;
t60 = -t135 * Ifges(4,4) - t136 * Ifges(4,2) + Ifges(4,6) * t314;
t58 = -t135 * Ifges(5,5) + Ifges(5,6) * t314 + t136 * Ifges(5,3);
t46 = -pkin(3) * t328 - t57;
t39 = t278 - t326;
t36 = -pkin(5) * t185 + t186 * t409 - t112;
t32 = -t150 * pkin(4) - t50;
t27 = -t270 + t326;
t26 = -t136 * pkin(4) - t33;
t25 = -t286 + t273;
t22 = (pkin(5) * t347 - t344) * t257 + t270;
t21 = -t283 + t273;
t20 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t19 = -t149 * pkin(5) + t150 * t409 - t50;
t18 = -qJD(1) * t286 + t274;
t13 = pkin(5) * t314 - t17;
t12 = t259 * t36 + t262 * t67;
t11 = -t259 * t67 + t262 * t36;
t10 = t54 * Ifges(7,1) + t55 * Ifges(7,4) - t135 * Ifges(7,5);
t9 = t54 * Ifges(7,4) + t55 * Ifges(7,2) - t135 * Ifges(7,6);
t4 = -qJD(6) * t15 + t19 * t262 - t21 * t259;
t3 = qJD(6) * t14 + t19 * t259 + t21 * t262;
t23 = [(-m(3) * t215 + m(4) * t196 - mrSges(3,1) * t364 + mrSges(4,1) * t213 + mrSges(4,2) * t214) * t195 - t420 * t360 / 0.2e1 - t364 * t426 + (t1 * t151 - t152 * t2 - t5 * t74 + t6 * t75) * mrSges(7,3) + ((t264 * (t323 + (t366 + t379) * t257) + (-t213 * t428 + t214 * t432 - t360 * t429) * t261) * qJD(1) + t285 * (Ifges(3,5) * t264 - Ifges(3,6) * t261) + t264 * t156 + (t92 + t90) * t261) * t348 / 0.2e1 + m(4) * (t115 * t38 + t116 * t37 - t352 * t57 + t56 * t98) + (Ifges(7,1) * t74 + Ifges(7,4) * t75) * t404 + (t63 + t62 + t8) * t214 / 0.2e1 + (t61 + t58) * t213 / 0.2e1 + ((-t373 + t379) * t321 + (Ifges(3,1) * t264 - t380) * t320 - 0.2e1 * pkin(1) * (mrSges(3,1) * t261 + mrSges(3,2) * t264) * t342) * t257 ^ 2 + m(3) * (t194 * t351 + t205 * t206) + (t194 * t360 + t195 * t361 - t215 * t313 - t314 * t351) * mrSges(3,3) + (t213 * t434 - t214 * t431 + t360 * t430) * t136 / 0.2e1 + (Ifges(7,4) * t74 + Ifges(7,2) * t75) * t406 + (t442 + t455) * t150 + t214 * t450 + (t443 - Ifges(6,3) * t390 - Ifges(6,6) * t395 + Ifges(4,6) * t397 + t429 * t391 - t430 * t396 + t432 * t394 + (Ifges(6,5) * t213 - Ifges(6,6) * t214 + Ifges(6,3) * t360) * t467 + t463) * t328 + (Ifges(7,5) * t74 + Ifges(7,6) * t75) * t398 + t214 * t413 + t74 * t414 + t75 * t415 + (-m(3) * t202 + m(4) * t159 - t423) * t422 + t14 * t29 + t15 * t30 + t22 * t65 - t214 * t451 + t69 * t20 + t45 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + t78 * t73 + t3 * t85 + t4 * t86 + t99 * t71 + t100 * t105 + t70 * t106 + t101 * t108 + t80 * t109 + (Ifges(4,4) * t214 - Ifges(4,2) * t213 - Ifges(4,6) * t360) * t436 + (Ifges(6,4) * t213 - Ifges(6,2) * t214 + Ifges(6,6) * t360) * t437 + (Ifges(7,5) * t152 + Ifges(7,6) * t151 - t432 * t360 + t433 * t213 + (Ifges(7,3) + t435) * t214) * t438 + t32 * t111 + t50 * t113 + t115 * t107 + t116 * t110 + t27 * t137 + t56 * t138 + t57 * t139 + t46 * t140 + t25 * t141 + t39 * t142 + t151 * t9 / 0.2e1 + t152 * t10 / 0.2e1 + t13 * (-mrSges(7,1) * t151 + mrSges(7,2) * t152) + (-t441 + t454) * t149 + t196 * t72 + t206 * t200 - t213 * t60 / 0.2e1 + t54 * (Ifges(7,1) * t152 + Ifges(7,4) * t151 + Ifges(7,5) * t214) / 0.2e1 + t55 * (Ifges(7,4) * t152 + Ifges(7,2) * t151 + Ifges(7,6) * t214) / 0.2e1 + t33 * (mrSges(5,1) * t213 - mrSges(5,3) * t214) + t26 * (mrSges(6,1) * t214 + mrSges(6,2) * t213) + m(7) * (t1 * t15 + t13 * t69 + t14 * t2 + t22 * t45 + t3 * t6 + t4 * t5) + m(6) * (t17 * t80 + t18 * t70 + t25 * t44 + t26 * t78 + t27 * t53 + t32 * t49) + m(5) * (t100 * t31 + t101 * t34 + t33 * t99 + t39 * t84 + t46 * t82 + t50 * t77) + t31 * (-mrSges(5,2) * t213 - mrSges(5,3) * t360) + t18 * (-mrSges(6,2) * t360 - mrSges(6,3) * t214) + t38 * (-mrSges(4,1) * t360 - mrSges(4,3) * t214) + t37 * (mrSges(4,2) * t360 - mrSges(4,3) * t213) + t17 * (mrSges(6,1) * t360 - mrSges(6,3) * t213) + t34 * (mrSges(5,1) * t360 + mrSges(5,2) * t214) + t364 * (-Ifges(3,6) * t314 + t230) / 0.2e1 - t327 * t368; t230 + (-t164 * t6 + t165 * t5) * mrSges(7,3) + (t31 * mrSges(5,2) + t37 * mrSges(4,3) - t33 * mrSges(5,1) - t195 * mrSges(4,1) - t26 * mrSges(6,2) + t17 * mrSges(6,3) + t9 * t389 + t10 * t387 - t61 / 0.2e1 + t60 / 0.2e1 - t58 / 0.2e1 - t13 * t298 + t297 * t440 + t295 * t439 + t302 * mrSges(7,3) + t306 * t136 + (t374 / 0.2e1 - t372 / 0.2e1 + t309) * t135 + (m(4) * t37 + m(5) * t31 + t105 + t110) * pkin(9) + (-mrSges(7,3) * t300 + t262 * t415 + t292 * t398 + t294 * t406 + t296 * t404 - t299 * t45 + t389 * t43) * qJD(6)) * t263 - m(4) * (-t117 * t352 + t118 * t98 + t159 * t205) + t423 * t205 + (t8 / 0.2e1 + t62 / 0.2e1 + t63 / 0.2e1 + t52 / 0.2e1 + t51 / 0.2e1 - t33 * mrSges(5,3) + t195 * mrSges(4,2) + t413 + t26 * mrSges(6,1) + t34 * mrSges(5,2) - t38 * mrSges(4,3) - t18 * mrSges(6,3) + t309 * t136 + (-Ifges(7,3) / 0.2e1 + t310) * t135 + (-m(4) * t38 - t107 + t424) * pkin(9) + t305) * t260 - m(5) * (t102 * t84 + t104 * t82 + t119 * t77) + m(5) * (t211 * t77 + t224 * t33) + (Ifges(7,5) * t165 + Ifges(7,6) * t164) * t399 + (Ifges(7,1) * t165 + Ifges(7,4) * t164) * t405 + (-t119 + t211) * t113 + (t287 * mrSges(4,3) + t288 * mrSges(5,2) + (Ifges(7,5) * t358 - Ifges(7,6) * t359) * t398 + (Ifges(7,1) * t358 - Ifges(7,4) * t359) * t404 + (Ifges(7,4) * t358 - Ifges(7,2) * t359) * t406 + t358 * t414 + t359 * t416 + (m(4) * t287 + m(5) * t288) * pkin(9) + t45 * (mrSges(7,1) * t359 + mrSges(7,2) * t358) + (-t358 * t5 - t359 * t6) * mrSges(7,3) + (-t353 * pkin(9) + t441 - t91 / 0.2e1) * t263 + (t354 * pkin(9) + t442) * t260) * qJD(3) + t456 * t111 - t426 + (Ifges(7,4) * t165 + Ifges(7,2) * t164) * t407 + t164 * t416 + t457 * t137 + t458 * t141 + (-t17 * t233 + t18 * t232 + t217 * t26 + t44 * t458 + t456 * t49 + t457 * t53) * m(6) + t459 * t65 + t460 * t85 + (t1 * t122 + t121 * t2 + t13 * t233 + t45 * t459 + t460 * t6 + t461 * t5) * m(7) + t461 * t86 + ((-t443 + (-Ifges(6,5) * t263 - Ifges(6,6) * t260) * t383 + t462 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t231 + (t260 * t432 + t263 * t428) * qJD(2) / 0.2e1 - t92 / 0.2e1 - t90 / 0.2e1 + t307 * t186 + t308 * t185 + (t322 / 0.2e1 + (pkin(1) * mrSges(3,1) + t380 / 0.2e1) * t257) * qJD(1)) * t261 + (-t236 / 0.2e1 + Ifges(3,5) * t383 - t156 / 0.2e1 + ((pkin(1) * mrSges(3,2) - t366 / 0.2e1 + t373 / 0.2e1) * t257 - t323 / 0.2e1) * qJD(1) + t368 + (t185 * t306 - t186 * t309 + t267) * t260 + (-t185 * t309 + t186 * t310 - t419) * t263) * t264) * t349 + t121 * t29 + t122 * t30 - t118 * t138 - t117 * t139 - t104 * t140 - t102 * t142 - t165 * t43 / 0.2e1 - t45 * (-mrSges(7,1) * t164 + mrSges(7,2) * t165) + (-m(4) * t195 - t72) * pkin(2) - t195 * mrSges(3,1) - t202 * t200 + t217 * t73 + t224 * t71 + t232 * t106 + (t20 - t109) * t233; t13 * t299 + (-t11 * t5 - t12 * t6 + t13 * t258 + t425 * t45) * m(7) + t303 * mrSges(7,3) + t9 * t387 - t354 * t352 - ((-m(7) * t301 + t289) * qJD(6) - m(7) * t303 + t291) * t343 + (t105 - t109) * qJ(4) + t418 * qJD(6) + t420 - t17 * mrSges(6,1) + t18 * mrSges(6,2) + (-pkin(3) * t34 + qJ(4) * t31 - t112 * t77 + t324 * t84 - t82 * t98) * m(5) + ((t306 - t310) * t185 + t468 + t341 * t186 + t267 + t418) * t186 + (-t175 / 0.2e1 + t419 + t465) * t185 + (-t17 * qJ(4) - t18 * t410 - t425 * t53 - t44 * t67 - t49 * t83) * m(6) - t410 * t106 + t31 * mrSges(5,3) - t34 * mrSges(5,1) - t37 * mrSges(4,2) + t38 * mrSges(4,1) - t12 * t85 - t11 * t86 - pkin(3) * t108 + t292 * t437 + t294 * t439 + t296 * t440 - t83 * t111 - t112 * t113 - t67 * t141 + t340 * qJD(4) + t258 * t20 - t259 * t10 / 0.2e1 + t353 * t98 + t365 * t66; t289 * qJD(6) + t340 * t231 + (-t111 + t113 + t289) * t186 - m(5) * (-t186 * t77 - t231 * t84) + t291 + t106 + (-t180 * t301 + t231 * t45 - t303) * m(7) + (-t186 * t49 - t231 * t53 + t18) * m(6) + t424; t259 * t30 + t262 * t29 - t365 * t185 + t290 * qJD(6) + (t141 + t290) * t186 + t73 + (-t180 * t300 - t185 * t45 + t302) * m(7) + (t185 * t53 + t186 * t44 + t26) * m(6); -t45 * (mrSges(7,1) * t124 + mrSges(7,2) * t123) + (Ifges(7,1) * t123 - t378) * t405 + t42 * t404 + (Ifges(7,5) * t123 - Ifges(7,6) * t124) * t399 - t5 * t85 + t6 * t86 + (t123 * t5 + t124 * t6) * mrSges(7,3) + t305 + t8 + (-Ifges(7,2) * t124 + t120 + t43) * t407;];
tauc  = t23(:);
