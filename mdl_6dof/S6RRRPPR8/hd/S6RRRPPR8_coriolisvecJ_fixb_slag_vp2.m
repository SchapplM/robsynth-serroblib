% Calculate vector of centrifugal and Coriolis load on the joints for
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:37
% EndTime: 2019-03-09 16:04:24
% DurationCPUTime: 23.95s
% Computational Cost: add. (8521->756), mult. (22263->971), div. (0->0), fcn. (16044->8), ass. (0->320)
t263 = cos(qJ(2));
t362 = cos(pkin(6));
t330 = pkin(1) * t362;
t245 = t263 * t330;
t260 = sin(qJ(2));
t256 = sin(pkin(6));
t347 = qJD(1) * t256;
t328 = t260 * t347;
t202 = -pkin(8) * t328 + qJD(1) * t245;
t284 = qJD(1) * t362 + qJD(2);
t159 = -t284 * pkin(2) - t202;
t327 = t263 * t347;
t231 = -qJD(3) + t327;
t389 = -t231 / 0.2e1;
t259 = sin(qJ(3));
t262 = cos(qJ(3));
t186 = t259 * t284 + t262 * t328;
t392 = t186 / 0.2e1;
t180 = qJD(6) + t186;
t396 = t180 / 0.2e1;
t278 = t262 * t284;
t185 = t259 * t328 - t278;
t258 = sin(qJ(6));
t261 = cos(qJ(6));
t124 = t185 * t261 + t231 * t258;
t402 = t124 / 0.2e1;
t123 = -t185 * t258 + t231 * t261;
t404 = t123 / 0.2e1;
t433 = Ifges(4,1) + Ifges(5,1);
t172 = t186 * qJ(5);
t244 = t260 * t330;
t358 = t256 * t263;
t349 = pkin(8) * t358 + t244;
t197 = t362 * pkin(9) + t349;
t160 = qJD(2) * pkin(9) + qJD(1) * t197;
t198 = (-pkin(2) * t263 - pkin(9) * t260 - pkin(1)) * t256;
t169 = qJD(1) * t198;
t350 = t259 * t160 - t262 * t169;
t322 = qJD(4) + t350;
t314 = -t172 + t322;
t408 = pkin(3) + pkin(4);
t44 = t231 * t408 + t314;
t394 = t185 / 0.2e1;
t395 = -t185 / 0.2e1;
t473 = Ifges(4,4) + Ifges(6,4);
t464 = Ifges(5,5) * t394 + t473 * t395;
t430 = Ifges(5,4) + Ifges(4,5);
t472 = Ifges(6,6) + t430;
t77 = t185 * pkin(3) - t186 * qJ(4) + t159;
t280 = qJD(5) - t77;
t49 = -pkin(4) * t185 + t280;
t407 = -pkin(4) - pkin(10);
t28 = pkin(5) * t186 + t185 * t407 + t280;
t341 = pkin(3) - t407;
t35 = t231 * t341 + t314;
t5 = -t258 * t35 + t261 * t28;
t6 = t258 * t28 + t261 * t35;
t474 = t49 * mrSges(6,1) + t5 * mrSges(7,1) + t159 * mrSges(4,2) - t6 * mrSges(7,2) - t77 * mrSges(5,3) - t44 * mrSges(6,3) + 0.2e1 * Ifges(7,5) * t402 + 0.2e1 * Ifges(7,6) * t404 + 0.2e1 * Ifges(7,3) * t396 + t464 + (Ifges(6,2) + t433) * t392 + t472 * t389;
t281 = (pkin(2) * t260 - pkin(9) * t263) * t256;
t203 = qJD(1) * t281;
t117 = -t259 * t202 + t203 * t262;
t343 = qJD(3) * t262;
t379 = pkin(9) - qJ(5);
t471 = -qJD(5) * t259 + t343 * t379 + t117;
t222 = t262 * qJ(4) * t327;
t348 = qJ(4) * t343 + t259 * qJD(4);
t470 = t348 - t222;
t393 = -t186 / 0.2e1;
t432 = Ifges(6,1) + Ifges(5,3);
t279 = pkin(5) * t262 - t259 * t341;
t469 = qJD(3) * t279 - (-t244 + (-pkin(8) + t279) * t358) * qJD(1) + t470;
t316 = t341 * t260;
t353 = t262 * qJ(5);
t336 = t263 * t353;
t468 = -(-t316 - t336) * t347 + t471;
t467 = Ifges(5,5) * t393;
t388 = t231 / 0.2e1;
t466 = -qJD(1) / 0.2e1;
t381 = -qJD(2) / 0.2e1;
t118 = t262 * t202 + t259 * t203;
t102 = qJ(4) * t328 + t118;
t344 = qJD(3) * t259;
t465 = -qJD(5) * t262 - t344 * t379 - t102;
t320 = Ifges(3,6) * t362;
t378 = Ifges(3,4) * t260;
t461 = (t320 + (Ifges(3,2) * t263 + t378) * t256) * t466;
t460 = Ifges(3,6) * t381;
t427 = Ifges(5,2) + Ifges(4,3);
t426 = Ifges(4,6) - Ifges(5,6);
t224 = -t262 * pkin(3) - t259 * qJ(4) - pkin(2);
t217 = t262 * pkin(4) - t224;
t179 = pkin(5) * t259 + pkin(10) * t262 + t217;
t232 = t379 * t259;
t122 = t179 * t258 + t232 * t261;
t459 = -qJD(6) * t122 - t258 * t468 + t261 * t469;
t121 = t179 * t261 - t232 * t258;
t458 = qJD(6) * t121 + t258 * t469 + t261 * t468;
t355 = t259 * t263;
t457 = -(pkin(5) * t260 + qJ(5) * t355) * t347 + t465;
t456 = -(-t260 * t408 - t336) * t347 + t471;
t455 = qJ(5) * t259 * t327 - t465;
t337 = t408 * t259;
t454 = -(-t244 + (-pkin(8) - t337) * t358) * qJD(1) - qJD(3) * t337 + t470;
t82 = pkin(3) * t231 + t322;
t453 = t82 * mrSges(5,2) + t350 * mrSges(4,3);
t221 = t231 * qJ(4);
t98 = t262 * t160 + t259 * t169;
t84 = -t221 + t98;
t452 = -t84 * mrSges(5,2) - t98 * mrSges(4,3);
t67 = qJ(5) * t185 + t98;
t53 = t221 - t67;
t451 = -t159 * mrSges(4,1) - t77 * mrSges(5,1) - t49 * mrSges(6,2) + t53 * mrSges(6,3) + Ifges(5,6) * t388 + t467 + t473 * t392 + (Ifges(6,5) + Ifges(4,6)) * t389 + (Ifges(4,2) + t432) * t395;
t340 = qJD(1) * qJD(2);
t319 = t263 * t340;
t312 = t256 * t319;
t359 = t256 * t260;
t335 = t259 * t359;
t315 = qJD(3) * t335;
t135 = qJD(1) * t315 - qJD(3) * t278 - t262 * t312;
t204 = qJD(2) * t281;
t193 = qJD(1) * t204;
t215 = -pkin(8) * t359 + t245;
t206 = t215 * qJD(2);
t194 = qJD(1) * t206;
t38 = -t160 * t343 - t169 * t344 + t193 * t262 - t259 * t194;
t273 = qJ(5) * t135 - qJD(5) * t186 - t38;
t346 = qJD(2) * t256;
t282 = t316 * t346;
t16 = -qJD(1) * t282 + t273;
t136 = qJD(3) * t186 + t259 * t312;
t420 = t349 * qJD(2);
t195 = qJD(1) * t420;
t33 = t136 * pkin(3) + t135 * qJ(4) - t186 * qJD(4) + t195;
t7 = -t135 * pkin(5) + t136 * t407 - t33;
t1 = qJD(6) * t5 + t16 * t261 + t258 * t7;
t449 = t1 * mrSges(7,2);
t2 = -qJD(6) * t6 - t16 * t258 + t261 * t7;
t448 = t2 * mrSges(7,1);
t214 = t259 * t362 + t262 * t359;
t325 = t263 * t346;
t149 = qJD(3) * t214 + t259 * t325;
t150 = -t315 + (qJD(3) * t362 + t325) * t262;
t50 = t149 * pkin(3) - t150 * qJ(4) - t214 * qJD(4) + t420;
t205 = t349 * qJD(1);
t441 = -t350 * mrSges(4,1) - t82 * mrSges(5,1) - t53 * mrSges(6,1) - t98 * mrSges(4,2) + t44 * mrSges(6,2) - t205 * mrSges(3,3) + t84 * mrSges(5,3) + Ifges(6,5) * t395 + Ifges(6,6) * t392 + Ifges(6,3) * t389 + t460 + t461;
t431 = -Ifges(4,4) + Ifges(5,5);
t440 = Ifges(6,4) * t393 + Ifges(6,5) * t388 - Ifges(4,2) * t395 - t426 * t389 + t431 * t392 + t432 * t394 - t451;
t429 = -Ifges(6,4) + Ifges(5,5);
t439 = Ifges(4,4) * t395 - Ifges(6,2) * t393 - Ifges(6,6) * t388 + t430 * t389 + t433 * t392 + t429 * t394 + t474;
t318 = t260 * t340;
t313 = t256 * t318;
t54 = qJD(6) * t123 + t136 * t261 - t258 * t313;
t438 = -t54 / 0.2e1;
t55 = -qJD(6) * t124 - t136 * t258 - t261 * t313;
t437 = -t55 / 0.2e1;
t436 = -t135 / 0.2e1;
t435 = t135 / 0.2e1;
t434 = -t136 / 0.2e1;
t428 = Ifges(6,5) - Ifges(5,6);
t424 = t194 * mrSges(3,2);
t66 = -t172 + t350;
t423 = qJD(4) + t66;
t108 = -mrSges(5,1) * t313 - t135 * mrSges(5,2);
t34 = -pkin(3) * t313 - t38;
t422 = m(5) * t34 + t108;
t421 = mrSges(3,1) * t284 - mrSges(4,1) * t185 - mrSges(4,2) * t186 - mrSges(3,3) * t328;
t302 = -t1 * t261 + t2 * t258;
t418 = (Ifges(6,3) + t427) * t313 + (-Ifges(6,5) - t426) * t136 - t472 * t135;
t306 = -Ifges(6,6) / 0.2e1 - Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t417 = t306 * t231 + t453 + t474;
t374 = Ifges(7,4) * t261;
t294 = -Ifges(7,2) * t258 + t374;
t375 = Ifges(7,4) * t258;
t296 = Ifges(7,1) * t261 - t375;
t297 = mrSges(7,1) * t258 + mrSges(7,2) * t261;
t300 = t258 * t6 + t261 * t5;
t370 = Ifges(7,6) * t258;
t372 = Ifges(7,5) * t261;
t385 = -t261 / 0.2e1;
t387 = t258 / 0.2e1;
t397 = -t180 / 0.2e1;
t403 = -t124 / 0.2e1;
t405 = -t123 / 0.2e1;
t376 = Ifges(7,4) * t124;
t42 = Ifges(7,2) * t123 + Ifges(7,6) * t180 + t376;
t120 = Ifges(7,4) * t123;
t43 = Ifges(7,1) * t124 + Ifges(7,5) * t180 + t120;
t45 = -pkin(5) * t231 - t53;
t416 = mrSges(7,3) * t300 + (-t370 + t372) * t397 + t294 * t405 + t296 * t403 - t45 * t297 + t385 * t43 + t387 * t42;
t414 = -t42 / 0.2e1;
t413 = t42 / 0.2e1;
t412 = t43 / 0.2e1;
t411 = Ifges(6,4) * t434 + Ifges(6,2) * t436 + Ifges(6,6) * t313 / 0.2e1;
t52 = Ifges(7,5) * t54;
t51 = Ifges(7,6) * t55;
t257 = qJ(4) + pkin(5);
t377 = Ifges(3,4) * t263;
t371 = Ifges(3,2) * t260;
t366 = t202 * mrSges(3,3);
t364 = t260 * Ifges(3,1);
t137 = mrSges(6,1) * t231 - mrSges(6,3) * t185;
t65 = -mrSges(7,1) * t123 + mrSges(7,2) * t124;
t363 = t65 - t137;
t213 = -t262 * t362 + t335;
t360 = qJ(5) * t213;
t357 = t258 * t259;
t356 = t259 * t261;
t354 = t261 * t263;
t106 = mrSges(6,2) * t313 + t135 * mrSges(6,3);
t138 = mrSges(4,2) * t231 - mrSges(4,3) * t185;
t142 = -mrSges(5,2) * t185 - mrSges(5,3) * t231;
t352 = -t138 - t142;
t139 = -mrSges(4,1) * t231 - mrSges(4,3) * t186;
t140 = mrSges(5,1) * t231 + mrSges(5,2) * t186;
t351 = t139 - t140;
t112 = t186 * pkin(3) + t185 * qJ(4);
t116 = t262 * t197 + t259 * t198;
t345 = qJD(2) * t260;
t342 = qJD(4) * t263;
t8 = -Ifges(7,3) * t135 + t51 + t52;
t339 = Ifges(6,4) / 0.2e1 + Ifges(4,4) / 0.2e1;
t338 = t142 + t363;
t196 = -t362 * pkin(2) - t215;
t326 = t256 * t345;
t324 = t256 * t342;
t321 = Ifges(3,5) * t362;
t73 = -t135 * mrSges(6,1) + t136 * mrSges(6,2);
t115 = -t259 * t197 + t198 * t262;
t309 = -Ifges(4,1) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(6,2) / 0.2e1;
t308 = Ifges(5,5) / 0.2e1 - t339;
t307 = -Ifges(4,6) / 0.2e1 - Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1;
t305 = -Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1 - Ifges(6,1) / 0.2e1;
t304 = t448 - t449;
t301 = t1 * t258 + t2 * t261;
t299 = t258 * t5 - t261 * t6;
t101 = pkin(3) * t358 - t115;
t298 = mrSges(7,1) * t261 - mrSges(7,2) * t258;
t295 = Ifges(7,1) * t258 + t374;
t293 = Ifges(7,2) * t261 + t375;
t291 = Ifges(7,5) * t258 + Ifges(7,6) * t261;
t29 = -mrSges(7,1) * t135 - mrSges(7,3) * t54;
t30 = mrSges(7,2) * t135 + mrSges(7,3) * t55;
t290 = -t258 * t29 + t261 * t30;
t99 = t213 * pkin(3) - t214 * qJ(4) + t196;
t40 = pkin(5) * t214 + t213 * t407 - t99;
t70 = pkin(4) * t358 - qJ(5) * t214 + t101;
t64 = pkin(10) * t358 + t70;
t14 = -t258 * t64 + t261 * t40;
t15 = t258 * t40 + t261 * t64;
t85 = -mrSges(7,2) * t180 + mrSges(7,3) * t123;
t86 = mrSges(7,1) * t180 - mrSges(7,3) * t124;
t289 = -t258 * t86 + t261 * t85;
t288 = -t258 * t85 - t261 * t86;
t287 = -t259 * t84 + t262 * t82;
t286 = -t259 * t98 + t262 * t350;
t285 = t408 * t326;
t100 = -qJ(4) * t358 + t116;
t57 = -t197 * t343 - t198 * t344 + t204 * t262 - t259 * t206;
t151 = -t213 * t258 + t256 * t354;
t152 = t213 * t261 + t258 * t358;
t37 = -t160 * t344 + t169 * t343 + t259 * t193 + t262 * t194;
t56 = -t197 * t344 + t198 * t343 + t259 * t204 + t262 * t206;
t277 = qJ(4) * t326 + t56;
t31 = qJ(4) * t313 - t231 * qJD(4) + t37;
t272 = -qJ(5) * t150 - qJD(5) * t214 - t57;
t269 = qJ(5) * t149 + qJD(5) * t213 + t277;
t17 = -qJ(5) * t136 - qJD(5) * t185 - t31;
t266 = t307 * t231 + t451 - t452;
t236 = Ifges(3,4) * t327;
t233 = pkin(9) * t262 - t353;
t230 = Ifges(3,5) * t312;
t211 = pkin(3) * t344 - t348;
t200 = -mrSges(3,2) * t284 + mrSges(3,3) * t327;
t165 = (-t258 * t260 + t259 * t354) * t347;
t164 = (-t258 * t355 - t260 * t261) * t347;
t156 = Ifges(3,1) * t328 + t284 * Ifges(3,5) + t236;
t141 = -mrSges(6,2) * t231 - mrSges(6,3) * t186;
t119 = -t222 + (t244 + (pkin(3) * t259 + pkin(8)) * t358) * qJD(1);
t113 = mrSges(5,1) * t185 - mrSges(5,3) * t186;
t111 = mrSges(6,1) * t186 + mrSges(6,2) * t185;
t110 = -mrSges(4,2) * t313 - mrSges(4,3) * t136;
t109 = -mrSges(6,1) * t313 - mrSges(6,3) * t136;
t107 = mrSges(4,1) * t313 + mrSges(4,3) * t135;
t105 = -mrSges(5,2) * t136 + mrSges(5,3) * t313;
t104 = -pkin(3) * t328 - t117;
t92 = t186 * Ifges(5,4) - t231 * Ifges(5,2) + t185 * Ifges(5,6);
t90 = t186 * Ifges(4,5) - t185 * Ifges(4,6) - t231 * Ifges(4,3);
t83 = -pkin(4) * t186 - t112;
t80 = -t100 - t360;
t78 = -pkin(4) * t213 - t99;
t75 = qJD(6) * t151 + t149 * t261 - t258 * t326;
t74 = -qJD(6) * t152 - t149 * t258 - t261 * t326;
t72 = mrSges(4,1) * t136 - mrSges(4,2) * t135;
t71 = mrSges(5,1) * t136 + mrSges(5,3) * t135;
t69 = -t257 * t358 + t116 + t360;
t63 = -t135 * Ifges(4,1) - t136 * Ifges(4,4) + Ifges(4,5) * t313;
t62 = -t135 * Ifges(5,1) + Ifges(5,4) * t313 + t136 * Ifges(5,5);
t61 = t136 * Ifges(6,1) + t135 * Ifges(6,4) - Ifges(6,5) * t313;
t60 = -t135 * Ifges(4,4) - t136 * Ifges(4,2) + Ifges(4,6) * t313;
t58 = -t135 * Ifges(5,5) + Ifges(5,6) * t313 + t136 * Ifges(5,3);
t46 = -pkin(3) * t326 - t57;
t39 = t277 - t324;
t36 = -pkin(5) * t185 + t186 * t407 - t112;
t32 = -t149 * pkin(4) - t50;
t27 = -t269 + t324;
t26 = -t136 * pkin(4) - t33;
t25 = -t285 + t272;
t22 = (pkin(5) * t345 - t342) * t256 + t269;
t21 = -t282 + t272;
t20 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t19 = t150 * pkin(5) + t149 * t407 - t50;
t18 = -qJD(1) * t285 + t273;
t13 = pkin(5) * t313 - t17;
t12 = t258 * t36 + t261 * t67;
t11 = -t258 * t67 + t261 * t36;
t10 = Ifges(7,1) * t54 + Ifges(7,4) * t55 - Ifges(7,5) * t135;
t9 = t54 * Ifges(7,4) + t55 * Ifges(7,2) - t135 * Ifges(7,6);
t4 = -qJD(6) * t15 + t19 * t261 - t21 * t258;
t3 = qJD(6) * t14 + t19 * t258 + t21 * t261;
t23 = [m(3) * (t194 * t349 + t205 * t206) + (t63 + t62 + t8) * t214 / 0.2e1 + (t61 + t58) * t213 / 0.2e1 + (-0.2e1 * pkin(1) * (mrSges(3,1) * t260 + mrSges(3,2) * t263) * t340 + (-t371 + t377) * t319 + (Ifges(3,1) * t263 - t378) * t318) * t256 ^ 2 + (t194 * t358 + t195 * t359 - t215 * t312 - t313 * t349) * mrSges(3,3) + m(4) * (t115 * t38 + t116 * t37 - t350 * t57 + t56 * t98) + (-m(3) * t202 + m(4) * t159 - t421) * t420 + (-m(3) * t215 + m(4) * t196 - mrSges(3,1) * t362 + mrSges(4,1) * t213 + mrSges(4,2) * t214) * t195 - t418 * t358 / 0.2e1 + t18 * (-mrSges(6,2) * t358 - mrSges(6,3) * t214) + t38 * (-mrSges(4,1) * t358 - mrSges(4,3) * t214) + t37 * (mrSges(4,2) * t358 - mrSges(4,3) * t213) + t17 * (mrSges(6,1) * t358 - mrSges(6,3) * t213) + t34 * (mrSges(5,1) * t358 + mrSges(5,2) * t214) + t31 * (-mrSges(5,2) * t213 - mrSges(5,3) * t358) + (Ifges(7,4) * t75 + Ifges(7,2) * t74) * t404 - t325 * t366 + (Ifges(7,1) * t75 + Ifges(7,4) * t74) * t402 + t214 * t411 + t75 * t412 + t74 * t413 + t362 * (-Ifges(3,6) * t313 + t230) / 0.2e1 + (t1 * t151 - t152 * t2 - t5 * t75 + t6 * t74) * mrSges(7,3) + (Ifges(4,4) * t214 - Ifges(4,2) * t213 - Ifges(4,6) * t358) * t434 + (Ifges(6,4) * t213 - Ifges(6,2) * t214 + Ifges(6,6) * t358) * t435 + (Ifges(7,5) * t152 + Ifges(7,6) * t151 - t430 * t358 + t431 * t213 + (Ifges(7,3) + t433) * t214) * t436 - t362 * t424 + (t441 + t461 + t427 * t389 - t428 * t394 + t430 * t392 - Ifges(6,6) * t393 + Ifges(4,6) * t395 - Ifges(6,3) * t388 + (Ifges(6,5) * t213 - Ifges(6,6) * t214 + Ifges(6,3) * t358) * t466) * t326 + t14 * t29 + t15 * t30 + (((-t213 * t426 + t214 * t430 - t358 * t427) * t260 + t263 * (t321 + (t364 + t377) * t256)) * qJD(1) + (t92 + t90) * t260 + t284 * (Ifges(3,5) * t263 - Ifges(3,6) * t260) + t263 * t156) * t346 / 0.2e1 + (t213 * t432 + t214 * t429 + t358 * t428) * t136 / 0.2e1 + t22 * t65 + t69 * t20 + t45 * (-mrSges(7,1) * t74 + mrSges(7,2) * t75) + t78 * t73 + t3 * t85 + t4 * t86 + t99 * t71 + t100 * t105 + t70 * t106 + t101 * t108 + t80 * t109 + t32 * t111 + t50 * t113 + t115 * t107 + t116 * t110 + t27 * t137 + t56 * t138 + t57 * t139 + t46 * t140 + t25 * t141 + t39 * t142 + t151 * t9 / 0.2e1 + t152 * t10 / 0.2e1 + t13 * (-mrSges(7,1) * t151 + mrSges(7,2) * t152) - t214 * t449 + m(7) * (t1 * t15 + t13 * t69 + t14 * t2 + t22 * t45 + t3 * t6 + t4 * t5) + m(6) * (t17 * t80 + t18 * t70 + t25 * t44 + t26 * t78 + t27 * t53 + t32 * t49) + m(5) * (t100 * t31 + t101 * t34 + t33 * t99 + t39 * t84 + t46 * t82 + t50 * t77) + t196 * t72 + t206 * t200 - t213 * t60 / 0.2e1 + t54 * (Ifges(7,1) * t152 + Ifges(7,4) * t151 + Ifges(7,5) * t214) / 0.2e1 + t55 * (Ifges(7,4) * t152 + Ifges(7,2) * t151 + Ifges(7,6) * t214) / 0.2e1 + t33 * (mrSges(5,1) * t213 - mrSges(5,3) * t214) + t26 * (mrSges(6,1) * t214 + mrSges(6,2) * t213) + (t440 + t452) * t149 + (t439 + t453) * t150 + (Ifges(7,5) * t75 + Ifges(7,6) * t74) * t396 + t214 * t448; -m(4) * (-t117 * t350 + t118 * t98 + t159 * t205) + (t286 * mrSges(4,3) + t287 * mrSges(5,2) + t45 * (mrSges(7,1) * t357 + mrSges(7,2) * t356) + (m(4) * t286 + m(5) * t287) * pkin(9) + (Ifges(7,4) * t356 - Ifges(7,2) * t357) * t404 + t356 * t412 + t357 * t414 + (Ifges(7,5) * t356 - Ifges(7,6) * t357) * t396 + (Ifges(7,1) * t356 - Ifges(7,4) * t357) * t402 + (-t356 * t5 - t357 * t6) * mrSges(7,3) + (-t351 * pkin(9) + t439) * t262 + (t352 * pkin(9) + t440) * t259) * qJD(3) + m(5) * (t211 * t77 + t224 * t33) + (-t164 * t6 + t165 * t5) * mrSges(7,3) + (-m(4) * t195 - t72) * pkin(2) + (Ifges(7,1) * t165 + Ifges(7,4) * t164) * t403 + (Ifges(7,4) * t165 + Ifges(7,2) * t164) * t405 + t164 * t414 + (Ifges(7,5) * t165 + Ifges(7,6) * t164) * t397 + t230 + (t20 - t109) * t233 + t421 * t205 + (-t18 * mrSges(6,3) - t33 * mrSges(5,3) + t195 * mrSges(4,2) + t26 * mrSges(6,1) + t34 * mrSges(5,2) - t38 * mrSges(4,3) + t63 / 0.2e1 + t62 / 0.2e1 + t411 + t8 / 0.2e1 + t52 / 0.2e1 + t51 / 0.2e1 + t308 * t136 + (-Ifges(7,3) / 0.2e1 + t309) * t135 + (-m(4) * t38 - t107 + t422) * pkin(9) + t304) * t259 - m(5) * (t102 * t84 + t104 * t82 + t119 * t77) + (-t119 + t211) * t113 + ((-t441 + (Ifges(6,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t231 + t306 * t186 - t307 * t185 + (t320 / 0.2e1 + (pkin(1) * mrSges(3,1) + t378 / 0.2e1) * t256) * qJD(1) - t92 / 0.2e1 + (-Ifges(6,5) * t262 - Ifges(6,6) * t259) * t381 - t90 / 0.2e1 + t460 + (t259 * t430 + t262 * t426) * qJD(2) / 0.2e1) * t260 + (-t156 / 0.2e1 - t236 / 0.2e1 + (-t321 / 0.2e1 + (-t364 / 0.2e1 + t371 / 0.2e1 + pkin(1) * mrSges(3,2)) * t256) * qJD(1) + t366 + Ifges(3,5) * t381 + (t185 * t305 - t186 * t308 + t266) * t259 + (-t185 * t308 + t186 * t309 - t417) * t262) * t263) * t347 + (-t33 * mrSges(5,1) - t195 * mrSges(4,1) - t26 * mrSges(6,2) + t17 * mrSges(6,3) + t31 * mrSges(5,2) + t37 * mrSges(4,3) + t9 * t387 + t10 * t385 - t13 * t297 + t296 * t438 + t294 * t437 - t61 / 0.2e1 + t60 / 0.2e1 - t58 / 0.2e1 + t301 * mrSges(7,3) + t305 * t136 + (t372 / 0.2e1 - t370 / 0.2e1 + t308) * t135 + (m(4) * t37 + m(5) * t31 + t105 + t110) * pkin(9) + (-mrSges(7,3) * t299 + t261 * t413 + t291 * t396 + t293 * t404 + t295 * t402 - t298 * t45 + t387 * t43) * qJD(6)) * t262 + t121 * t29 + t122 * t30 - t424 - t118 * t138 - t117 * t139 - t104 * t140 - t102 * t142 - t165 * t43 / 0.2e1 - t45 * (-mrSges(7,1) * t164 + mrSges(7,2) * t165) - t195 * mrSges(3,1) - t202 * t200 + t217 * t73 + t224 * t71 + t232 * t106 + t454 * t111 + t455 * t137 + t456 * t141 + (-t17 * t233 + t18 * t232 + t217 * t26 + t44 * t456 + t454 * t49 + t455 * t53) * m(6) + t457 * t65 + t458 * t85 + t459 * t86 + (t1 * t122 + t121 * t2 + t13 * t233 + t45 * t457 + t458 * t6 + t459 * t5) * m(7); t416 * qJD(6) - t352 * t350 + (-t17 * qJ(4) - t18 * t408 - t423 * t53 - t44 * t67 - t49 * t83) * m(6) - t408 * t106 + t418 + t338 * qJD(4) + (t417 + t464) * t185 + t351 * t98 + t13 * t298 + (t105 - t109) * qJ(4) + t302 * mrSges(7,3) + t9 * t385 + t363 * t66 + (-t11 * t5 - t12 * t6 + t13 * t257 + t423 * t45) * m(7) - (-m(7) * t302 + t290 + (-m(7) * t300 + t288) * qJD(6)) * t341 + t291 * t435 + t293 * t437 + t295 * t438 - t17 * mrSges(6,1) + t18 * mrSges(6,2) + (t467 + t266 + t339 * t186 + (t305 - t309) * t185 + t416) * t186 + t31 * mrSges(5,3) - t34 * mrSges(5,1) - t37 * mrSges(4,2) + t38 * mrSges(4,1) + (-pkin(3) * t34 + qJ(4) * t31 - t112 * t77 + t322 * t84 - t82 * t98) * m(5) - t12 * t85 - t11 * t86 - pkin(3) * t108 - t83 * t111 - t112 * t113 - t67 * t141 + t257 * t20 - t258 * t10 / 0.2e1; t288 * qJD(6) + t338 * t231 + (-t111 + t113 + t288) * t186 - m(5) * (-t186 * t77 - t231 * t84) + t290 + t106 + (-t180 * t300 + t231 * t45 - t302) * m(7) + (-t186 * t49 - t231 * t53 + t18) * m(6) + t422; t258 * t30 + t261 * t29 - t363 * t185 + t289 * qJD(6) + (t141 + t289) * t186 + t73 + (-t180 * t299 - t185 * t45 + t301) * m(7) + (t185 * t53 + t186 * t44 + t26) * m(6); -t45 * (mrSges(7,1) * t124 + mrSges(7,2) * t123) + (Ifges(7,1) * t123 - t376) * t403 + t42 * t402 + (Ifges(7,5) * t123 - Ifges(7,6) * t124) * t397 - t5 * t85 + t6 * t86 + (t123 * t5 + t124 * t6) * mrSges(7,3) + t304 + t8 + (-Ifges(7,2) * t124 + t120 + t43) * t405;];
tauc  = t23(:);
