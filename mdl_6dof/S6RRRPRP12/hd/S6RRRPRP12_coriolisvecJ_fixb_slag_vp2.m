% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:43
% EndTime: 2019-03-09 17:52:42
% DurationCPUTime: 30.07s
% Computational Cost: add. (10843->806), mult. (28195->1050), div. (0->0), fcn. (20695->8), ass. (0->323)
t273 = cos(qJ(2));
t267 = sin(pkin(6));
t369 = qJD(1) * t267;
t347 = t273 * t369;
t250 = -qJD(3) + t347;
t381 = cos(pkin(6));
t337 = t381 * qJD(1);
t331 = pkin(1) * t337;
t270 = sin(qJ(2));
t348 = t270 * t369;
t216 = -pkin(8) * t348 + t273 * t331;
t289 = (pkin(2) * t270 - pkin(9) * t273) * t267;
t217 = qJD(1) * t289;
t269 = sin(qJ(3));
t272 = cos(qJ(3));
t141 = -t269 * t216 + t217 * t272;
t266 = t272 * pkin(4);
t366 = qJD(3) * t272;
t419 = pkin(4) + pkin(9);
t420 = pkin(3) + pkin(10);
t473 = -(t266 * t273 - t270 * t420) * t369 + t141 + t419 * t366;
t367 = qJD(3) * t269;
t336 = pkin(3) * t367 - qJD(4) * t269;
t219 = pkin(8) * t347 + t270 * t331;
t334 = t269 * t347;
t351 = pkin(3) * t334 + t219;
t472 = -t336 + t351 + t250 * (pkin(10) * t269 - qJ(4) * t272);
t445 = Ifges(6,1) + Ifges(7,1);
t443 = Ifges(7,4) + Ifges(6,5);
t333 = t272 * t347;
t471 = t333 - t366;
t293 = t337 + qJD(2);
t286 = t272 * t293;
t199 = t269 * t348 - t286;
t470 = -t199 / 0.2e1;
t441 = Ifges(7,5) - Ifges(6,4);
t439 = Ifges(6,6) - Ifges(7,6);
t338 = -qJ(4) * t269 - pkin(2);
t229 = -t272 * t420 + t338;
t251 = t419 * t269;
t268 = sin(qJ(5));
t271 = cos(qJ(5));
t364 = qJD(5) * t271;
t365 = qJD(5) * t268;
t462 = -t229 * t365 + t251 * t364 + t268 * t473 - t472 * t271;
t431 = t271 * t229 + t268 * t251;
t461 = -qJD(5) * t431 + t472 * t268 + t271 * t473;
t178 = -pkin(2) * t293 - t216;
t200 = t269 * t293 + t272 * t348;
t279 = -t200 * qJ(4) + t178;
t104 = t199 * pkin(3) + t279;
t179 = pkin(9) * t293 + t219;
t212 = (-pkin(2) * t273 - pkin(9) * t270 - pkin(1)) * t267;
t190 = qJD(1) * t212;
t123 = t272 * t179 + t269 * t190;
t238 = t250 * qJ(4);
t110 = t238 - t123;
t116 = -t250 * Ifges(5,5) - t200 * Ifges(5,6) + t199 * Ifges(5,3);
t119 = t200 * Ifges(4,4) - t199 * Ifges(4,2) - t250 * Ifges(4,6);
t357 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1;
t360 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t469 = t104 * mrSges(5,2) + t123 * mrSges(4,3) - t116 / 0.2e1 + t119 / 0.2e1 - t110 * mrSges(5,1) - t178 * mrSges(4,1) + t357 * t200 + t360 * t250;
t148 = t199 * t268 - t250 * t271;
t295 = t271 * t199 + t250 * t268;
t92 = -pkin(4) * t199 + t123;
t77 = -t238 + t92;
t27 = -pkin(5) * t295 - qJ(6) * t148 + t77;
t122 = t179 * t269 - t272 * t190;
t290 = pkin(4) * t200 + t122;
t66 = t250 * t420 + qJD(4) + t290;
t69 = t199 * t420 + t279;
t19 = -t268 * t69 + t271 * t66;
t20 = t268 * t66 + t271 * t69;
t301 = t19 * t268 - t20 * t271;
t194 = qJD(5) + t200;
t434 = qJD(6) - t19;
t11 = -pkin(5) * t194 + t434;
t12 = qJ(6) * t194 + t20;
t302 = t11 * t268 + t12 * t271;
t388 = Ifges(7,5) * t268;
t306 = -Ifges(7,3) * t271 + t388;
t391 = Ifges(6,4) * t268;
t312 = Ifges(6,2) * t271 + t391;
t319 = mrSges(7,1) * t271 + mrSges(7,3) * t268;
t321 = mrSges(6,1) * t271 - mrSges(6,2) * t268;
t401 = t271 / 0.2e1;
t402 = -t271 / 0.2e1;
t404 = -t268 / 0.2e1;
t408 = -t194 / 0.2e1;
t414 = -t148 / 0.2e1;
t415 = -t295 / 0.2e1;
t416 = t295 / 0.2e1;
t145 = Ifges(6,4) * t295;
t389 = Ifges(7,5) * t295;
t437 = t148 * t445 + t443 * t194 + t145 - t389;
t387 = Ifges(7,5) * t271;
t390 = Ifges(6,4) * t271;
t452 = t268 * t445 - t387 + t390;
t144 = Ifges(7,5) * t148;
t60 = Ifges(7,6) * t194 - Ifges(7,3) * t295 + t144;
t392 = Ifges(6,4) * t148;
t63 = Ifges(6,2) * t295 + Ifges(6,6) * t194 + t392;
t468 = t401 * t60 + t402 * t63 + t404 * t437 + t27 * t319 + t306 * t416 + t312 * t415 + t77 * t321 + (t268 * t443 + t271 * t439) * t408 + t452 * t414 - t302 * mrSges(7,2) + t301 * mrSges(6,3);
t467 = Ifges(5,6) * t470;
t466 = t200 / 0.2e1;
t446 = Ifges(5,1) + Ifges(4,3);
t465 = Ifges(5,4) - Ifges(4,5);
t442 = Ifges(5,5) - Ifges(4,6);
t464 = -qJ(6) * t471 + qJD(6) * t269 + t462;
t463 = pkin(5) * t471 - t461;
t142 = t272 * t216 + t269 * t217;
t129 = -qJ(4) * t348 - t142;
t115 = -pkin(4) * t334 - t129;
t185 = t268 * t348 - t271 * t334;
t376 = t268 * t273;
t186 = (t269 * t376 + t270 * t271) * t369;
t305 = pkin(5) * t271 + qJ(6) * t268;
t291 = -pkin(4) - t305;
t304 = -pkin(5) * t268 + qJ(6) * t271;
t460 = -pkin(5) * t185 + qJ(6) * t186 - t115 + (qJD(5) * t304 + qJD(6) * t268) * t272 + (-pkin(9) + t291) * t367;
t459 = -t419 * t367 - t115;
t458 = -t468 - t469;
t457 = t27 * mrSges(7,1) + t77 * mrSges(6,1) + t60 / 0.2e1 - t63 / 0.2e1 - t12 * mrSges(7,2) - t20 * mrSges(6,3);
t456 = t250 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t466 + t467;
t363 = qJD(1) * qJD(2);
t340 = t273 * t363;
t329 = t267 * t340;
t378 = t267 * t270;
t354 = t269 * t378;
t332 = qJD(3) * t354;
t155 = qJD(1) * t332 - qJD(3) * t286 - t272 * t329;
t156 = qJD(3) * t200 + t269 * t329;
t341 = t270 * t363;
t330 = t267 * t341;
t75 = qJD(5) * t295 + t156 * t268 + t271 * t330;
t76 = qJD(5) * t148 - t271 * t156 + t268 * t330;
t455 = -t155 * t443 + t441 * t76 + t445 * t75;
t450 = -qJD(4) - t122;
t108 = pkin(3) * t250 - t450;
t192 = Ifges(4,4) * t199;
t121 = t200 * Ifges(4,1) - t250 * Ifges(4,5) - t192;
t356 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t358 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t359 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t449 = -Ifges(5,4) / 0.2e1;
t361 = Ifges(4,5) / 0.2e1 + t449;
t61 = t148 * Ifges(6,5) + Ifges(6,6) * t295 + t194 * Ifges(6,3);
t62 = t148 * Ifges(7,4) + t194 * Ifges(7,2) - Ifges(7,6) * t295;
t276 = t295 * t356 - t148 * t359 + t194 * t358 - t250 * t361 - t104 * mrSges(5,3) - t11 * mrSges(7,1) - t20 * mrSges(6,2) + t456 + t121 / 0.2e1 + t61 / 0.2e1 + t62 / 0.2e1 + t108 * mrSges(5,1) + t12 * mrSges(7,3) + t122 * mrSges(4,3) + t178 * mrSges(4,2) + t19 * mrSges(6,1);
t362 = Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1;
t454 = (-t199 * t357 + t200 * t362 + t276) * t272;
t453 = -t268 * t439 + t271 * t443;
t451 = t271 * t445 + t388 - t391;
t412 = -t155 / 0.2e1;
t447 = -t156 / 0.2e1;
t440 = Ifges(7,2) + Ifges(6,3);
t14 = Ifges(6,5) * t75 - Ifges(6,6) * t76 - Ifges(6,3) * t155;
t15 = Ifges(7,4) * t75 - Ifges(7,2) * t155 + Ifges(7,6) * t76;
t438 = t14 + t15;
t349 = pkin(1) * t381;
t227 = -pkin(8) * t378 + t273 * t349;
t220 = t227 * qJD(2);
t208 = qJD(1) * t220;
t436 = t208 * mrSges(3,2);
t355 = Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1;
t435 = t199 * t355;
t133 = -t155 * mrSges(5,1) + mrSges(5,2) * t330;
t218 = qJD(2) * t289;
t207 = qJD(1) * t218;
t58 = -t179 * t366 - t190 * t367 + t207 * t272 - t269 * t208;
t54 = -pkin(3) * t330 - t58;
t433 = m(5) * t54 + t133;
t432 = mrSges(3,1) * t293 - mrSges(4,1) * t199 - mrSges(4,2) * t200 - mrSges(3,3) * t348;
t428 = t155 * t465 + t156 * t442 + t330 * t446;
t368 = qJD(2) * t267;
t346 = t270 * t368;
t297 = t420 * t346;
t30 = -pkin(4) * t155 - qJD(1) * t297 - t58;
t377 = t267 * t273;
t228 = pkin(8) * t377 + t270 * t349;
t221 = t228 * qJD(2);
t209 = qJD(1) * t221;
t278 = t155 * qJ(4) - t200 * qJD(4) + t209;
t36 = t156 * t420 + t278;
t4 = -qJD(5) * t20 - t268 * t36 + t271 * t30;
t224 = -t272 * t381 + t354;
t210 = -pkin(2) * t381 - t227;
t225 = t269 * t381 + t272 * t378;
t284 = -t225 * qJ(4) + t210;
t105 = t224 * t420 + t284;
t211 = pkin(9) * t381 + t228;
t139 = -t269 * t211 + t212 * t272;
t127 = pkin(3) * t377 - t139;
t94 = pkin(4) * t225 + pkin(10) * t377 + t127;
t397 = t271 * t105 + t268 * t94;
t345 = t273 * t368;
t169 = -t332 + (qJD(3) * t381 + t345) * t272;
t79 = -t211 * t366 - t212 * t367 + t218 * t272 - t269 * t220;
t40 = pkin(4) * t169 - t297 - t79;
t168 = qJD(3) * t225 + t269 * t345;
t280 = -t169 * qJ(4) - t225 * qJD(4) + t221;
t52 = t168 * t420 + t280;
t8 = -qJD(5) * t397 - t268 * t52 + t271 * t40;
t424 = t75 / 0.2e1;
t423 = -t76 / 0.2e1;
t422 = t76 / 0.2e1;
t421 = Ifges(5,2) * t412 + Ifges(5,6) * t447 + t330 * t449;
t413 = t148 / 0.2e1;
t407 = t194 / 0.2e1;
t403 = t268 / 0.2e1;
t43 = -mrSges(6,1) * t155 - mrSges(6,3) * t75;
t44 = t155 * mrSges(7,1) + t75 * mrSges(7,2);
t399 = t43 - t44;
t45 = mrSges(6,2) * t155 - mrSges(6,3) * t76;
t46 = -mrSges(7,2) * t76 - mrSges(7,3) * t155;
t398 = t45 + t46;
t396 = mrSges(6,3) * t295;
t395 = mrSges(6,3) * t148;
t394 = Ifges(3,4) * t270;
t393 = Ifges(3,4) * t273;
t386 = Ifges(3,2) * t270;
t385 = t216 * mrSges(3,3);
t384 = t219 * mrSges(3,3);
t383 = t270 * Ifges(3,1);
t379 = qJ(4) * t199;
t109 = t200 * t420 + t379;
t33 = t271 * t109 + t268 * t92;
t160 = mrSges(5,1) * t199 + mrSges(5,3) * t250;
t90 = -mrSges(6,1) * t295 + mrSges(6,2) * t148;
t382 = -t160 + t90;
t380 = Ifges(3,6) * qJD(2);
t375 = t268 * t420;
t374 = t271 * t420;
t111 = mrSges(7,2) * t295 + mrSges(7,3) * t194;
t112 = -mrSges(6,2) * t194 + t396;
t373 = t111 + t112;
t113 = mrSges(6,1) * t194 - t395;
t114 = -mrSges(7,1) * t194 + mrSges(7,2) * t148;
t372 = -t113 + t114;
t159 = -mrSges(4,1) * t250 - mrSges(4,3) * t200;
t161 = mrSges(5,1) * t200 - mrSges(5,2) * t250;
t371 = t159 - t161;
t158 = mrSges(4,2) * t250 - mrSges(4,3) * t199;
t370 = -t160 + t158;
t140 = t272 * t211 + t269 * t212;
t252 = t272 * pkin(9) + t266;
t353 = t271 * t377;
t343 = Ifges(3,5) * t381;
t342 = Ifges(3,6) * t381;
t3 = t268 * t30 + t271 * t36 + t66 * t364 - t365 * t69;
t1 = -qJ(6) * t155 + qJD(6) * t194 + t3;
t2 = pkin(5) * t155 - t4;
t323 = -t1 * t268 + t2 * t271;
t322 = -t268 * t3 - t271 * t4;
t320 = mrSges(6,1) * t268 + mrSges(6,2) * t271;
t318 = mrSges(7,1) * t268 - mrSges(7,3) * t271;
t313 = -Ifges(6,2) * t268 + t390;
t307 = Ifges(7,3) * t268 + t387;
t34 = -t105 * t268 + t271 * t94;
t32 = -t109 * t268 + t271 * t92;
t165 = -t229 * t268 + t251 * t271;
t126 = qJ(4) * t377 - t140;
t170 = t224 * t271 + t267 * t376;
t7 = -t105 * t365 + t268 * t40 + t271 * t52 + t94 * t364;
t57 = -t179 * t367 + t190 * t366 + t269 * t207 + t272 * t208;
t78 = -t211 * t367 + t212 * t366 + t269 * t218 + t272 * t220;
t106 = -pkin(4) * t224 - t126;
t285 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t49 = -qJ(4) * t330 + t250 * qJD(4) - t57;
t29 = -pkin(4) * t156 - t49;
t59 = -qJ(4) * t346 + qJD(4) * t377 - t78;
t282 = (t342 + (Ifges(3,2) * t273 + t394) * t267) * qJD(1);
t41 = -pkin(4) * t168 - t59;
t255 = Ifges(3,4) * t347;
t248 = Ifges(3,5) * t329;
t242 = -pkin(3) * t272 + t338;
t240 = qJ(4) - t304;
t223 = -qJ(4) * t366 + t336;
t215 = -mrSges(3,2) * t293 + mrSges(3,3) * t347;
t213 = qJD(5) * t305 - qJD(6) * t271 + qJD(4);
t193 = t272 * t305 + t252;
t175 = Ifges(3,1) * t348 + Ifges(3,5) * t293 + t255;
t174 = t282 + t380;
t171 = t224 * t268 - t353;
t162 = -pkin(5) * t269 - t165;
t157 = qJ(6) * t269 + t431;
t143 = -qJ(4) * t333 + t351;
t138 = -mrSges(5,2) * t199 - mrSges(5,3) * t200;
t136 = pkin(3) * t200 + t379;
t135 = -mrSges(4,2) * t330 - mrSges(4,3) * t156;
t134 = mrSges(4,1) * t330 + mrSges(4,3) * t155;
t132 = mrSges(5,1) * t156 - mrSges(5,3) * t330;
t131 = -pkin(3) * t348 - t141;
t125 = t224 * pkin(3) + t284;
t120 = -t250 * Ifges(5,1) - t200 * Ifges(5,4) + t199 * Ifges(5,5);
t117 = t200 * Ifges(4,5) - t199 * Ifges(4,6) - t250 * Ifges(4,3);
t102 = qJD(5) * t170 + t168 * t268 + t271 * t346;
t101 = -t168 * t271 - qJD(5) * t353 + (qJD(5) * t224 + t346) * t268;
t96 = mrSges(4,1) * t156 - mrSges(4,2) * t155;
t95 = -mrSges(5,2) * t156 + mrSges(5,3) * t155;
t89 = -mrSges(7,1) * t295 - mrSges(7,3) * t148;
t88 = pkin(5) * t148 - qJ(6) * t295;
t83 = -t155 * Ifges(4,1) - t156 * Ifges(4,4) + Ifges(4,5) * t330;
t82 = -t155 * Ifges(4,4) - t156 * Ifges(4,2) + Ifges(4,6) * t330;
t80 = Ifges(5,5) * t330 + t155 * Ifges(5,6) + t156 * Ifges(5,3);
t68 = t168 * pkin(3) + t280;
t67 = -pkin(3) * t346 - t79;
t53 = t156 * pkin(3) + t278;
t47 = t200 * t291 - t122;
t42 = -pkin(5) * t170 - qJ(6) * t171 + t106;
t26 = -pkin(5) * t225 - t34;
t25 = qJ(6) * t225 + t397;
t24 = pkin(5) * t199 - t32;
t23 = -qJ(6) * t199 + t33;
t22 = mrSges(6,1) * t76 + mrSges(6,2) * t75;
t21 = mrSges(7,1) * t76 - mrSges(7,3) * t75;
t16 = Ifges(6,4) * t75 - Ifges(6,2) * t76 - Ifges(6,6) * t155;
t13 = Ifges(7,5) * t75 - Ifges(7,6) * t155 + Ifges(7,3) * t76;
t10 = pkin(5) * t101 - qJ(6) * t102 - qJD(6) * t171 + t41;
t9 = pkin(5) * t76 - qJ(6) * t75 - qJD(6) * t148 + t29;
t6 = -pkin(5) * t169 - t8;
t5 = qJ(6) * t169 + qJD(6) * t225 + t7;
t17 = [(Ifges(4,4) * t225 - Ifges(4,2) * t224 - Ifges(4,6) * t377) * t447 + m(3) * (t208 * t228 + t219 * t220) + (Ifges(7,5) * t171 + Ifges(7,6) * t225) * t422 - (t282 + t174) * t346 / 0.2e1 + (t121 + t62 + t61) * t169 / 0.2e1 + (-0.2e1 * pkin(1) * (mrSges(3,1) * t270 + mrSges(3,2) * t273) * t363 + (-t386 + t393) * t340 + (Ifges(3,1) * t273 - t394) * t341) * t267 ^ 2 + t397 * t45 + m(6) * (t106 * t29 + t19 * t8 + t20 * t7 + t3 * t397 + t34 * t4 + t41 * t77) + (t102 * t77 - t169 * t20 + t171 * t29 - t225 * t3) * mrSges(6,2) + (Ifges(6,4) * t102 + Ifges(6,6) * t169) * t416 + (Ifges(6,4) * t171 + Ifges(6,6) * t225) * t423 + t455 * t171 / 0.2e1 + m(4) * (-t122 * t79 + t123 * t78 + t139 * t58 + t140 * t57) + (t208 * t377 + t209 * t378 - t227 * t329 - t228 * t330) * mrSges(3,3) + (t1 * t225 - t102 * t27 + t12 * t169 - t171 * t9) * mrSges(7,3) + (t293 * (Ifges(3,5) * t273 - Ifges(3,6) * t270) + t273 * t175 + (t120 + t117) * t270 + ((t224 * t442 - t225 * t465 - t377 * t446) * t270 + t273 * (t343 + (t383 + t393) * t267)) * qJD(1)) * t368 / 0.2e1 - (t168 * t442 - t169 * t465 + t346 * t446) * t250 / 0.2e1 + (Ifges(7,5) * t102 + Ifges(7,6) * t169) * t415 + (-m(3) * t227 + m(4) * t210 - mrSges(3,1) * t381 + mrSges(4,1) * t224 + mrSges(4,2) * t225) * t209 + (-Ifges(6,2) * t416 + Ifges(7,3) * t415 - t439 * t407 + t441 * t413 + t457) * t101 + t4 * (mrSges(6,1) * t225 - mrSges(6,3) * t171) + t2 * (-mrSges(7,1) * t225 + mrSges(7,2) * t171) + t53 * (-mrSges(5,2) * t224 - mrSges(5,3) * t225) + t224 * t80 / 0.2e1 - t224 * t82 / 0.2e1 + t220 * t215 + t210 * t96 + t178 * (mrSges(4,1) * t168 + mrSges(4,2) * t169) + t11 * (-mrSges(7,1) * t169 + mrSges(7,2) * t102) + t104 * (-mrSges(5,2) * t168 - mrSges(5,3) * t169) + t19 * (mrSges(6,1) * t169 - mrSges(6,3) * t102) + t168 * t116 / 0.2e1 - t168 * t119 / 0.2e1 + t79 * t159 + t59 * t160 + t67 * t161 + t78 * t158 + t140 * t135 + t126 * t132 + t127 * t133 + t68 * t138 + t139 * t134 + t125 * t95 + t6 * t114 + t5 * t111 + t7 * t112 + t8 * t113 + t106 * t22 + t10 * t89 + t41 * t90 + t42 * t21 + t34 * t43 + t26 * t44 + t25 * t46 + (Ifges(4,4) * t169 - Ifges(4,2) * t168 + Ifges(4,6) * t346) * t470 + (Ifges(4,1) * t169 - Ifges(4,4) * t168 + Ifges(4,5) * t346) * t466 + t169 * t456 - t428 * t377 / 0.2e1 + (-m(3) * t216 + m(4) * t178 - t432) * t221 - t381 * t436 + t437 * t102 / 0.2e1 + (t83 + t438) * t225 / 0.2e1 + (t3 * mrSges(6,3) + t1 * mrSges(7,2) - t13 / 0.2e1 - t9 * mrSges(7,1) - t29 * mrSges(6,1) + t16 / 0.2e1 - Ifges(7,3) * t422 + Ifges(6,2) * t423 - t441 * t424) * t170 + (-Ifges(4,4) * t224 - Ifges(4,5) * t377 + t443 * t171 + t439 * t170 + (Ifges(4,1) + t440) * t225) * t412 + (t102 * t443 + t169 * t440) * t407 + (t171 * t445 + t225 * t443) * t424 + (t102 * t445 + t169 * t443) * t413 + t110 * (mrSges(5,1) * t168 - mrSges(5,3) * t346) + t123 * (-mrSges(4,2) * t346 - mrSges(4,3) * t168) + t199 * (Ifges(5,5) * t346 - Ifges(5,6) * t169 + Ifges(5,3) * t168) / 0.2e1 - t200 * (Ifges(5,4) * t346 - Ifges(5,2) * t169 + Ifges(5,6) * t168) / 0.2e1 + t108 * (mrSges(5,1) * t169 + mrSges(5,2) * t346) - t122 * (mrSges(4,1) * t346 - mrSges(4,3) * t169) + t225 * t421 + t156 * (-Ifges(5,5) * t377 - Ifges(5,6) * t225 + Ifges(5,3) * t224) / 0.2e1 + t155 * (-Ifges(5,4) * t377 - Ifges(5,2) * t225 + Ifges(5,6) * t224) / 0.2e1 + t54 * (mrSges(5,1) * t225 - mrSges(5,2) * t377) + t58 * (-mrSges(4,1) * t377 - mrSges(4,3) * t225) + t49 * (mrSges(5,1) * t224 + mrSges(5,3) * t377) + t57 * (mrSges(4,2) * t377 - mrSges(4,3) * t224) + t381 * (-Ifges(3,6) * t330 + t248) / 0.2e1 - t346 * t384 - t345 * t385 + m(7) * (t1 * t25 + t10 * t27 + t11 * t6 + t12 * t5 + t2 * t26 + t42 * t9) + m(5) * (t104 * t68 + t108 * t67 + t110 * t59 + t125 * t53 + t126 * t49 + t127 * t54); -t436 + t461 * t113 + (t165 * t4 + t19 * t461 + t20 * t462 + t252 * t29 + t3 * t431 + t459 * t77) * m(6) + t462 * t112 + t463 * t114 + (t1 * t157 + t11 * t463 + t12 * t464 + t162 * t2 + t193 * t9 + t27 * t460) * m(7) + t464 * t111 + (t454 + (t435 + t458) * t269 + (-t371 * t272 - t370 * t269 + m(5) * (t108 * t272 + t110 * t269) + m(4) * (t122 * t272 - t123 * t269)) * pkin(9)) * qJD(3) + (-Ifges(6,2) * t415 + Ifges(7,3) * t416 - t408 * t439 + t414 * t441 - t457) * t185 + (t29 * t321 + t9 * t319 + t306 * t423 + t312 * t422 + t13 * t401 + t16 * t402 + t82 / 0.2e1 - t80 / 0.2e1 - t49 * mrSges(5,1) + t57 * mrSges(4,3) + t53 * mrSges(5,2) - t209 * mrSges(4,1) - t355 * t156 + (-t268 * t359 + t271 * t356 - t357) * t155 + (t268 * t4 - t271 * t3) * mrSges(6,3) + (-t1 * t271 - t2 * t268) * mrSges(7,2) + (m(4) * t57 - m(5) * t49 - t132 + t135) * pkin(9) + (-t77 * t320 - t27 * t318 + t313 * t415 + t307 * t416 + t63 * t403 + (t19 * t271 + t20 * t268) * mrSges(6,3) + (-t11 * t271 + t12 * t268) * mrSges(7,2) + t451 * t414 + t453 * t408 + t437 * t402) * qJD(5) - t452 * t75 / 0.2e1 + (qJD(5) * t60 + t455) * t404) * t272 + t460 * t89 + ((-t380 / 0.2e1 + t384 + t122 * mrSges(4,1) - t108 * mrSges(5,2) + t110 * mrSges(5,3) + t123 * mrSges(4,2) - t117 / 0.2e1 - t120 / 0.2e1 + t174 / 0.2e1 + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t250 - t361 * t200 - t360 * t199 + (t342 / 0.2e1 + (pkin(1) * mrSges(3,1) + t394 / 0.2e1) * t267) * qJD(1) + (-t269 * t465 - t442 * t272) * qJD(2) / 0.2e1) * t270 + (-t175 / 0.2e1 + (-t343 / 0.2e1 + (t386 / 0.2e1 + pkin(1) * mrSges(3,2) - t383 / 0.2e1) * t267) * qJD(1) + t385 - qJD(2) * Ifges(3,5) / 0.2e1 - t255 / 0.2e1 + (-t435 + t469) * t269 - t454) * t273) * t369 + t459 * t90 + t431 * t45 + t252 * t22 + (-m(4) * t209 - t96) * pkin(2) + t242 * t95 - t216 * t215 - t209 * mrSges(3,1) + t193 * t21 + t165 * t43 - t141 * t159 - t129 * t160 - t131 * t161 + t162 * t44 + t157 * t46 - t142 * t158 + (-t437 / 0.2e1 - mrSges(6,2) * t77 - t11 * mrSges(7,2) + t19 * mrSges(6,3) + mrSges(7,3) * t27 + Ifges(6,4) * t415 + Ifges(7,5) * t416 + t408 * t443 + t414 * t445) * t186 + t248 + t432 * t219 + (t15 / 0.2e1 + t14 / 0.2e1 + t83 / 0.2e1 + t421 + t54 * mrSges(5,1) - t58 * mrSges(4,3) - t53 * mrSges(5,3) + t209 * mrSges(4,2) - t356 * t76 - t359 * t75 - t357 * t156 + (-t358 - t362) * t155 + (-m(4) * t58 - t134 + t433) * pkin(9) + t285) * t269 + m(5) * (t104 * t223 + t242 * t53) - m(4) * (-t122 * t141 + t123 * t142 + t178 * t219) - m(5) * (t104 * t143 + t108 * t131 + t110 * t129) + (t223 - t143) * t138; (-t132 + t22) * qJ(4) - m(7) * (t11 * t24 + t12 * t23 + t27 * t47) + (-pkin(3) * t54 - qJ(4) * t49 - t104 * t136 - t108 * t123 + t110 * t450) * m(5) + ((-t355 + t362) * t199 - t458) * t200 + t9 * t318 + t29 * t320 + t455 * t401 + m(7) * (-t1 * t375 + t2 * t374 + t27 * t213 + t9 * t240) + m(6) * (t29 * qJ(4) + t77 * qJD(4) - t3 * t375 - t374 * t4) + t428 + (t213 - t47) * t89 + t322 * mrSges(6,3) + t323 * mrSges(7,2) + (t276 + t467 - t192 / 0.2e1) * t199 + t451 * t424 + t453 * t412 + (-(-m(6) * t301 + m(7) * t302 + t268 * t372 + t271 * t373) * t420 + t468) * qJD(5) - m(6) * (t19 * t32 + t20 * t33 - t290 * t77) + t290 * t90 - (t268 * t398 + t271 * t399) * t420 + t240 * t21 - pkin(3) * t133 - t136 * t138 - t24 * t114 - t23 * t111 - t33 * t112 - t32 * t113 - t57 * mrSges(4,2) + t58 * mrSges(4,1) - t49 * mrSges(5,3) + t54 * mrSges(5,2) + t13 * t403 + t16 * t404 + t307 * t422 + t313 * t423 + t370 * t122 + t371 * t123 + t382 * qJD(4); t200 * t138 + (t89 + t382) * t250 + (t194 * t373 + t399) * t271 + (t194 * t372 + t398) * t268 - m(5) * (-t104 * t200 + t110 * t250) + (t194 * t302 + t250 * t27 - t323) * m(7) + (-t194 * t301 + t250 * t77 - t322) * m(6) + t433; t285 + (-t11 * t295 + t12 * t148) * mrSges(7,2) + (-t372 + t395) * t20 + (-t373 + t396) * t19 + (Ifges(7,3) * t148 + t389) * t416 - t77 * (mrSges(6,1) * t148 + mrSges(6,2) * t295) + t63 * t413 - t27 * (mrSges(7,1) * t148 - mrSges(7,3) * t295) + qJD(6) * t111 - t88 * t89 - pkin(5) * t44 + qJ(6) * t46 + (-t148 * t439 + t295 * t443) * t408 + (-pkin(5) * t2 + qJ(6) * t1 - t11 * t20 + t12 * t434 - t27 * t88) * m(7) + (-Ifges(6,2) * t148 + t145 + t437) * t415 + (t295 * t445 + t144 - t392 + t60) * t414 + t438; -t194 * t111 + t148 * t89 + 0.2e1 * (t2 / 0.2e1 + t12 * t408 + t27 * t413) * m(7) + t44;];
tauc  = t17(:);
