% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:49:49
% EndTime: 2018-11-23 17:50:14
% DurationCPUTime: 25.72s
% Computational Cost: add. (10843->810), mult. (28195->1055), div. (0->0), fcn. (20695->8), ass. (0->324)
t274 = cos(qJ(2));
t268 = sin(pkin(6));
t370 = qJD(1) * t268;
t349 = t274 * t370;
t250 = -qJD(3) + t349;
t382 = cos(pkin(6));
t338 = t382 * qJD(1);
t332 = pkin(1) * t338;
t271 = sin(qJ(2));
t350 = t271 * t370;
t216 = -pkin(8) * t350 + t274 * t332;
t290 = (pkin(2) * t271 - pkin(9) * t274) * t268;
t217 = qJD(1) * t290;
t270 = sin(qJ(3));
t273 = cos(qJ(3));
t141 = -t270 * t216 + t217 * t273;
t267 = t273 * pkin(4);
t367 = qJD(3) * t273;
t421 = pkin(4) + pkin(9);
t422 = pkin(3) + pkin(10);
t475 = -(t267 * t274 - t271 * t422) * t370 + t141 + t421 * t367;
t368 = qJD(3) * t270;
t337 = pkin(3) * t368 - qJD(4) * t270;
t219 = pkin(8) * t349 + t271 * t332;
t335 = t270 * t349;
t353 = pkin(3) * t335 + t219;
t474 = -t337 + t353 + t250 * (pkin(10) * t270 - qJ(4) * t273);
t451 = Ifges(6,1) + Ifges(7,1);
t449 = Ifges(7,4) + Ifges(6,5);
t334 = t273 * t349;
t473 = t334 - t367;
t447 = -Ifges(7,5) + Ifges(6,4);
t445 = Ifges(6,6) - Ifges(7,6);
t340 = -qJ(4) * t270 - pkin(2);
t229 = -t273 * t422 + t340;
t251 = t421 * t270;
t269 = sin(qJ(5));
t272 = cos(qJ(5));
t365 = qJD(5) * t272;
t366 = qJD(5) * t269;
t468 = -t229 * t366 + t251 * t365 + t269 * t475 - t474 * t272;
t437 = t272 * t229 + t269 * t251;
t467 = -qJD(5) * t437 + t474 * t269 + t272 * t475;
t294 = t338 + qJD(2);
t287 = t273 * t294;
t198 = t270 * t350 - t287;
t177 = -pkin(2) * t294 - t216;
t199 = t270 * t294 + t273 * t350;
t280 = -t199 * qJ(4) + t177;
t104 = t198 * pkin(3) + t280;
t178 = pkin(9) * t294 + t219;
t212 = (-pkin(2) * t274 - pkin(9) * t271 - pkin(1)) * t268;
t189 = qJD(1) * t212;
t123 = t273 * t178 + t270 * t189;
t238 = t250 * qJ(4);
t110 = t238 - t123;
t116 = -t250 * Ifges(5,5) - t199 * Ifges(5,6) + t198 * Ifges(5,3);
t119 = t199 * Ifges(4,4) - t198 * Ifges(4,2) - t250 * Ifges(4,6);
t359 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1;
t360 = Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1;
t472 = t104 * mrSges(5,2) + t123 * mrSges(4,3) - t116 / 0.2e1 + t119 / 0.2e1 - t110 * mrSges(5,1) - t177 * mrSges(4,1) + t359 * t199 - t360 * t250;
t148 = t198 * t269 - t250 * t272;
t296 = t272 * t198 + t250 * t269;
t92 = -pkin(4) * t198 + t123;
t77 = -t238 + t92;
t27 = -pkin(5) * t296 - qJ(6) * t148 + t77;
t122 = t178 * t270 - t273 * t189;
t291 = pkin(4) * t199 + t122;
t66 = t250 * t422 + qJD(4) + t291;
t69 = t198 * t422 + t280;
t19 = -t269 * t69 + t272 * t66;
t20 = t269 * t66 + t272 * t69;
t302 = t19 * t269 - t20 * t272;
t193 = qJD(5) + t199;
t440 = qJD(6) - t19;
t11 = -pkin(5) * t193 + t440;
t12 = qJ(6) * t193 + t20;
t303 = t11 * t269 + t12 * t272;
t320 = mrSges(7,1) * t272 + mrSges(7,3) * t269;
t322 = mrSges(6,1) * t272 - mrSges(6,2) * t269;
t402 = t272 / 0.2e1;
t403 = -t272 / 0.2e1;
t405 = -t269 / 0.2e1;
t145 = Ifges(6,4) * t296;
t390 = Ifges(7,5) * t296;
t443 = t148 * t451 + t449 * t193 + t145 - t390;
t144 = Ifges(7,5) * t148;
t60 = Ifges(7,6) * t193 - Ifges(7,3) * t296 + t144;
t393 = Ifges(6,4) * t148;
t63 = Ifges(6,2) * t296 + Ifges(6,6) * t193 + t393;
t471 = -t303 * mrSges(7,2) + t302 * mrSges(6,3) + t27 * t320 + t77 * t322 + t402 * t60 + t403 * t63 + t405 * t443;
t452 = Ifges(5,1) + Ifges(4,3);
t450 = Ifges(5,4) - Ifges(4,5);
t448 = Ifges(5,5) - Ifges(4,6);
t470 = -qJ(6) * t473 + qJD(6) * t270 + t468;
t469 = pkin(5) * t473 - t467;
t142 = t273 * t216 + t270 * t217;
t129 = -qJ(4) * t350 - t142;
t115 = -pkin(4) * t335 - t129;
t184 = t269 * t350 - t272 * t335;
t377 = t269 * t274;
t185 = (t270 * t377 + t271 * t272) * t370;
t306 = pkin(5) * t272 + qJ(6) * t269;
t292 = -pkin(4) - t306;
t305 = -pkin(5) * t269 + qJ(6) * t272;
t466 = -pkin(5) * t184 + qJ(6) * t185 - t115 + (qJD(5) * t305 + qJD(6) * t269) * t273 + (-pkin(9) + t292) * t368;
t465 = -t421 * t368 - t115;
t389 = Ifges(7,5) * t269;
t307 = -Ifges(7,3) * t272 + t389;
t392 = Ifges(6,4) * t269;
t313 = Ifges(6,2) * t272 + t392;
t408 = t193 / 0.2e1;
t414 = t148 / 0.2e1;
t416 = -t296 / 0.2e1;
t417 = t296 / 0.2e1;
t388 = Ifges(7,5) * t272;
t391 = Ifges(6,4) * t272;
t434 = t269 * t451 - t388 + t391;
t435 = t269 * t449 + t272 * t445;
t464 = -t307 * t416 - t313 * t417 - t408 * t435 - t414 * t434 + t471 + t472;
t463 = -mrSges(6,1) * t77 - mrSges(7,1) * t27 + t12 * mrSges(7,2) + t20 * mrSges(6,3) + t63 / 0.2e1 - t60 / 0.2e1;
t364 = qJD(1) * qJD(2);
t342 = t274 * t364;
t330 = t268 * t342;
t379 = t268 * t271;
t355 = t270 * t379;
t333 = qJD(3) * t355;
t155 = qJD(1) * t333 - qJD(3) * t287 - t273 * t330;
t156 = qJD(3) * t199 + t270 * t330;
t343 = t271 * t364;
t331 = t268 * t343;
t75 = qJD(5) * t296 + t156 * t269 + t272 * t331;
t76 = qJD(5) * t148 - t272 * t156 + t269 * t331;
t462 = -t449 * t155 - t447 * t76 + t451 * t75;
t458 = -qJD(4) - t122;
t108 = pkin(3) * t250 - t458;
t191 = Ifges(4,4) * t198;
t121 = t199 * Ifges(4,1) - t250 * Ifges(4,5) - t191;
t356 = -Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1;
t358 = -Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t457 = -Ifges(5,4) / 0.2e1;
t361 = Ifges(4,5) / 0.2e1 + t457;
t362 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t190 = Ifges(5,6) * t198;
t453 = -t250 / 0.2e1;
t454 = -t199 / 0.2e1;
t418 = Ifges(5,4) * t453 + Ifges(5,2) * t454 + t190 / 0.2e1;
t61 = t148 * Ifges(6,5) + Ifges(6,6) * t296 + t193 * Ifges(6,3);
t62 = t148 * Ifges(7,4) + t193 * Ifges(7,2) - Ifges(7,6) * t296;
t277 = t296 * t358 + t148 * t362 - t193 * t356 - t250 * t361 - t104 * mrSges(5,3) - t11 * mrSges(7,1) - t20 * mrSges(6,2) - t418 + t121 / 0.2e1 + t61 / 0.2e1 + t62 / 0.2e1 + t108 * mrSges(5,1) + t12 * mrSges(7,3) + t122 * mrSges(4,3) + t177 * mrSges(4,2) + t19 * mrSges(6,1);
t363 = Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1;
t461 = (-t198 * t359 + t199 * t363 + t277) * t273;
t460 = -t269 * t445 + t272 * t449;
t459 = t272 * t451 + t389 - t392;
t413 = -t155 / 0.2e1;
t455 = -t156 / 0.2e1;
t446 = -Ifges(7,2) - Ifges(6,3);
t14 = Ifges(6,5) * t75 - Ifges(6,6) * t76 - Ifges(6,3) * t155;
t15 = Ifges(7,4) * t75 - Ifges(7,2) * t155 + Ifges(7,6) * t76;
t444 = t14 + t15;
t351 = pkin(1) * t382;
t227 = -pkin(8) * t379 + t274 * t351;
t220 = t227 * qJD(2);
t207 = qJD(1) * t220;
t442 = t207 * mrSges(3,2);
t357 = Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1;
t441 = t198 * t357;
t133 = -t155 * mrSges(5,1) + mrSges(5,2) * t331;
t218 = qJD(2) * t290;
t206 = qJD(1) * t218;
t58 = -t178 * t367 - t189 * t368 + t206 * t273 - t270 * t207;
t54 = -pkin(3) * t331 - t58;
t439 = m(5) * t54 + t133;
t438 = mrSges(3,1) * t294 - mrSges(4,1) * t198 - mrSges(4,2) * t199 - mrSges(3,3) * t350;
t432 = t450 * t155 + t448 * t156 + t331 * t452;
t369 = qJD(2) * t268;
t348 = t271 * t369;
t298 = t422 * t348;
t30 = -pkin(4) * t155 - qJD(1) * t298 - t58;
t378 = t268 * t274;
t228 = pkin(8) * t378 + t271 * t351;
t221 = t228 * qJD(2);
t208 = qJD(1) * t221;
t279 = t155 * qJ(4) - t199 * qJD(4) + t208;
t36 = t156 * t422 + t279;
t4 = -qJD(5) * t20 - t269 * t36 + t272 * t30;
t339 = t382 * t273;
t224 = -t339 + t355;
t210 = -pkin(2) * t382 - t227;
t225 = t270 * t382 + t273 * t379;
t285 = -t225 * qJ(4) + t210;
t105 = t224 * t422 + t285;
t211 = pkin(9) * t382 + t228;
t139 = -t270 * t211 + t212 * t273;
t127 = pkin(3) * t378 - t139;
t94 = pkin(4) * t225 + pkin(10) * t378 + t127;
t398 = t272 * t105 + t269 * t94;
t347 = t274 * t369;
t167 = -qJD(3) * t339 - t273 * t347 + t333;
t79 = -t211 * t367 - t212 * t368 + t218 * t273 - t270 * t220;
t40 = -pkin(4) * t167 - t298 - t79;
t168 = qJD(3) * t225 + t270 * t347;
t281 = t167 * qJ(4) - t225 * qJD(4) + t221;
t52 = t168 * t422 + t281;
t8 = -qJD(5) * t398 - t269 * t52 + t272 * t40;
t426 = t75 / 0.2e1;
t425 = -t76 / 0.2e1;
t424 = t76 / 0.2e1;
t423 = Ifges(5,2) * t413 + Ifges(5,6) * t455 + t331 * t457;
t415 = -t148 / 0.2e1;
t409 = -t193 / 0.2e1;
t404 = t269 / 0.2e1;
t43 = -mrSges(6,1) * t155 - mrSges(6,3) * t75;
t44 = t155 * mrSges(7,1) + t75 * mrSges(7,2);
t400 = t43 - t44;
t45 = mrSges(6,2) * t155 - mrSges(6,3) * t76;
t46 = -mrSges(7,2) * t76 - mrSges(7,3) * t155;
t399 = t45 + t46;
t397 = mrSges(6,3) * t296;
t396 = mrSges(6,3) * t148;
t395 = Ifges(3,4) * t271;
t394 = Ifges(3,4) * t274;
t387 = Ifges(3,2) * t271;
t386 = t216 * mrSges(3,3);
t385 = t219 * mrSges(3,3);
t384 = t271 * Ifges(3,1);
t380 = qJ(4) * t198;
t109 = t199 * t422 + t380;
t33 = t272 * t109 + t269 * t92;
t160 = mrSges(5,1) * t198 + mrSges(5,3) * t250;
t90 = -mrSges(6,1) * t296 + mrSges(6,2) * t148;
t383 = -t160 + t90;
t381 = Ifges(3,6) * qJD(2);
t376 = t269 * t422;
t375 = t272 * t422;
t111 = mrSges(7,2) * t296 + mrSges(7,3) * t193;
t112 = -mrSges(6,2) * t193 + t397;
t374 = t111 + t112;
t113 = mrSges(6,1) * t193 - t396;
t114 = -mrSges(7,1) * t193 + mrSges(7,2) * t148;
t373 = -t113 + t114;
t158 = mrSges(4,2) * t250 - mrSges(4,3) * t198;
t372 = -t160 + t158;
t159 = -mrSges(4,1) * t250 - mrSges(4,3) * t199;
t161 = mrSges(5,1) * t199 - mrSges(5,2) * t250;
t371 = -t161 + t159;
t140 = t273 * t211 + t270 * t212;
t252 = t273 * pkin(9) + t267;
t345 = Ifges(3,5) * t382;
t344 = Ifges(3,6) * t382;
t3 = t269 * t30 + t272 * t36 + t66 * t365 - t366 * t69;
t1 = -qJ(6) * t155 + qJD(6) * t193 + t3;
t2 = pkin(5) * t155 - t4;
t324 = -t1 * t269 + t2 * t272;
t323 = -t269 * t3 - t272 * t4;
t321 = mrSges(6,1) * t269 + mrSges(6,2) * t272;
t319 = mrSges(7,1) * t269 - mrSges(7,3) * t272;
t314 = -Ifges(6,2) * t269 + t391;
t308 = Ifges(7,3) * t269 + t388;
t34 = -t105 * t269 + t272 * t94;
t32 = -t109 * t269 + t272 * t92;
t164 = -t229 * t269 + t251 * t272;
t126 = qJ(4) * t378 - t140;
t169 = t224 * t272 + t268 * t377;
t7 = -t105 * t366 + t269 * t40 + t272 * t52 + t94 * t365;
t57 = -t178 * t368 + t189 * t367 + t270 * t206 + t273 * t207;
t78 = -t211 * t368 + t212 * t367 + t270 * t218 + t273 * t220;
t106 = -pkin(4) * t224 - t126;
t286 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t49 = -qJ(4) * t331 + t250 * qJD(4) - t57;
t29 = -pkin(4) * t156 - t49;
t59 = -qJ(4) * t348 + qJD(4) * t378 - t78;
t283 = (t344 + (Ifges(3,2) * t274 + t395) * t268) * qJD(1);
t41 = -pkin(4) * t168 - t59;
t255 = Ifges(3,4) * t349;
t248 = Ifges(3,5) * t330;
t242 = -pkin(3) * t273 + t340;
t240 = qJ(4) - t305;
t223 = -qJ(4) * t367 + t337;
t215 = -mrSges(3,2) * t294 + mrSges(3,3) * t349;
t213 = qJD(5) * t306 - qJD(6) * t272 + qJD(4);
t192 = t273 * t306 + t252;
t174 = Ifges(3,1) * t350 + Ifges(3,5) * t294 + t255;
t173 = t283 + t381;
t170 = t224 * t269 - t272 * t378;
t162 = -pkin(5) * t270 - t164;
t157 = qJ(6) * t270 + t437;
t143 = -qJ(4) * t334 + t353;
t138 = -mrSges(5,2) * t198 - mrSges(5,3) * t199;
t136 = pkin(3) * t199 + t380;
t135 = -mrSges(4,2) * t331 - mrSges(4,3) * t156;
t134 = mrSges(4,1) * t331 + mrSges(4,3) * t155;
t132 = mrSges(5,1) * t156 - mrSges(5,3) * t331;
t131 = -pkin(3) * t350 - t141;
t125 = t224 * pkin(3) + t285;
t120 = -t250 * Ifges(5,1) - t199 * Ifges(5,4) + t198 * Ifges(5,5);
t117 = t199 * Ifges(4,5) - t198 * Ifges(4,6) - t250 * Ifges(4,3);
t102 = -t224 * t366 - t269 * t348 + (qJD(5) * t378 + t168) * t272;
t101 = qJD(5) * t169 + t168 * t269 + t272 * t348;
t96 = mrSges(4,1) * t156 - mrSges(4,2) * t155;
t95 = -mrSges(5,2) * t156 + mrSges(5,3) * t155;
t89 = -mrSges(7,1) * t296 - mrSges(7,3) * t148;
t88 = pkin(5) * t148 - qJ(6) * t296;
t83 = -t155 * Ifges(4,1) - t156 * Ifges(4,4) + Ifges(4,5) * t331;
t82 = -t155 * Ifges(4,4) - t156 * Ifges(4,2) + Ifges(4,6) * t331;
t80 = Ifges(5,5) * t331 + t155 * Ifges(5,6) + t156 * Ifges(5,3);
t68 = t168 * pkin(3) + t281;
t67 = -pkin(3) * t348 - t79;
t53 = t156 * pkin(3) + t279;
t47 = t199 * t292 - t122;
t42 = -pkin(5) * t169 - qJ(6) * t170 + t106;
t26 = -pkin(5) * t225 - t34;
t25 = qJ(6) * t225 + t398;
t24 = pkin(5) * t198 - t32;
t23 = -qJ(6) * t198 + t33;
t22 = mrSges(6,1) * t76 + mrSges(6,2) * t75;
t21 = mrSges(7,1) * t76 - mrSges(7,3) * t75;
t16 = Ifges(6,4) * t75 - Ifges(6,2) * t76 - Ifges(6,6) * t155;
t13 = Ifges(7,5) * t75 - Ifges(7,6) * t155 + Ifges(7,3) * t76;
t10 = -pkin(5) * t102 - qJ(6) * t101 - qJD(6) * t170 + t41;
t9 = pkin(5) * t76 - qJ(6) * t75 - qJD(6) * t148 + t29;
t6 = pkin(5) * t167 - t8;
t5 = -qJ(6) * t167 + qJD(6) * t225 + t7;
t17 = [((t120 + t117) * t271 + t294 * (Ifges(3,5) * t274 - Ifges(3,6) * t271) + t274 * t174 + ((t224 * t448 - t225 * t450 - t378 * t452) * t271 + t274 * (t345 + (t384 + t394) * t268)) * qJD(1)) * t369 / 0.2e1 + (t1 * t225 - t101 * t27 - t12 * t167 - t170 * t9) * mrSges(7,3) + t462 * t170 / 0.2e1 + t443 * t101 / 0.2e1 + (t83 + t444) * t225 / 0.2e1 - (t283 + t173) * t348 / 0.2e1 - (t121 + t62 + t61) * t167 / 0.2e1 + ((-t387 + t394) * t342 + (Ifges(3,1) * t274 - t395) * t343 - 0.2e1 * pkin(1) * (mrSges(3,1) * t271 + mrSges(3,2) * t274) * t364) * t268 ^ 2 + m(6) * (t106 * t29 + t19 * t8 + t20 * t7 + t3 * t398 + t34 * t4 + t41 * t77) + t398 * t45 + (Ifges(7,5) * t101 - Ifges(7,6) * t167) * t416 + (Ifges(7,5) * t170 + Ifges(7,6) * t225) * t424 + (-m(3) * t216 + m(4) * t177 - t438) * t221 + m(4) * (-t122 * t79 + t123 * t78 + t139 * t58 + t140 * t57) - t382 * t442 + (t1 * mrSges(7,2) + t3 * mrSges(6,3) - t13 / 0.2e1 - t9 * mrSges(7,1) - t29 * mrSges(6,1) + t16 / 0.2e1 - Ifges(7,3) * t424 + Ifges(6,2) * t425 + t447 * t426) * t169 + (-Ifges(4,4) * t224 - Ifges(4,5) * t378 + t449 * t170 + t445 * t169 + (Ifges(4,1) - t446) * t225) * t413 + (t101 * t449 + t167 * t446) * t408 + (t170 * t451 + t225 * t449) * t426 + (t101 * t451 - t167 * t449) * t414 + (t101 * t77 + t167 * t20 + t170 * t29 - t225 * t3) * mrSges(6,2) + (-m(3) * t227 + m(4) * t210 - mrSges(3,1) * t382 + mrSges(4,1) * t224 + mrSges(4,2) * t225) * t208 + (t207 * t378 + t208 * t379 - t227 * t330 - t228 * t331) * mrSges(3,3) + m(3) * (t207 * t228 + t219 * t220) - t347 * t386 - t348 * t385 + t382 * (-Ifges(3,6) * t331 + t248) / 0.2e1 + (Ifges(6,2) * t417 - Ifges(7,3) * t416 + t408 * t445 + t414 * t447 + t463) * t102 + t156 * (-Ifges(5,5) * t378 - Ifges(5,6) * t225 + Ifges(5,3) * t224) / 0.2e1 + t155 * (-Ifges(5,4) * t378 - Ifges(5,2) * t225 + Ifges(5,6) * t224) / 0.2e1 + t54 * (mrSges(5,1) * t225 - mrSges(5,2) * t378) + t58 * (-mrSges(4,1) * t378 - mrSges(4,3) * t225) + t49 * (mrSges(5,1) * t224 + mrSges(5,3) * t378) + t57 * (mrSges(4,2) * t378 - mrSges(4,3) * t224) + t4 * (mrSges(6,1) * t225 - mrSges(6,3) * t170) + t2 * (-mrSges(7,1) * t225 + mrSges(7,2) * t170) + t53 * (-mrSges(5,2) * t224 - mrSges(5,3) * t225) + t224 * t80 / 0.2e1 - t224 * t82 / 0.2e1 + t220 * t215 + t210 * t96 + t110 * (mrSges(5,1) * t168 - mrSges(5,3) * t348) + t123 * (-mrSges(4,2) * t348 - mrSges(4,3) * t168) - t198 * (-Ifges(4,4) * t167 - Ifges(4,2) * t168 + Ifges(4,6) * t348) / 0.2e1 + t198 * (Ifges(5,5) * t348 + Ifges(5,6) * t167 + Ifges(5,3) * t168) / 0.2e1 + t199 * (-Ifges(4,1) * t167 - Ifges(4,4) * t168 + Ifges(4,5) * t348) / 0.2e1 + t108 * (-mrSges(5,1) * t167 + mrSges(5,2) * t348) - t122 * (mrSges(4,1) * t348 + mrSges(4,3) * t167) + t177 * (mrSges(4,1) * t168 - mrSges(4,2) * t167) - t168 * t119 / 0.2e1 + t19 * (-mrSges(6,1) * t167 - mrSges(6,3) * t101) + t11 * (mrSges(7,1) * t167 + mrSges(7,2) * t101) + t104 * (-mrSges(5,2) * t168 + mrSges(5,3) * t167) + t168 * t116 / 0.2e1 + t59 * t160 + t67 * t161 + t78 * t158 + t79 * t159 + (Ifges(6,4) * t101 - Ifges(6,6) * t167) * t417 + (Ifges(6,4) * t170 + Ifges(6,6) * t225) * t425 + t126 * t132 + t127 * t133 + t68 * t138 + t139 * t134 + t140 * t135 + t125 * t95 + t5 * t111 + t7 * t112 + t8 * t113 + t6 * t114 + m(7) * (t1 * t25 + t10 * t27 + t11 * t6 + t12 * t5 + t2 * t26 + t42 * t9) + m(5) * (t104 * t68 + t108 * t67 + t110 * t59 + t125 * t53 + t126 * t49 + t127 * t54) + t106 * t22 + t10 * t89 + t41 * t90 + t34 * t43 + t26 * t44 + t25 * t46 + t42 * t21 - t432 * t378 / 0.2e1 + t167 * t418 + t225 * t423 + (t167 * t450 + t168 * t448 + t348 * t452) * t453 + (Ifges(5,4) * t348 + Ifges(5,2) * t167 + Ifges(5,6) * t168) * t454 + (Ifges(4,4) * t225 - Ifges(4,2) * t224 - Ifges(4,6) * t378) * t455; (t57 * mrSges(4,3) - t49 * mrSges(5,1) - t208 * mrSges(4,1) + t53 * mrSges(5,2) + t82 / 0.2e1 - t80 / 0.2e1 + t29 * t322 + t9 * t320 + t13 * t402 + t16 * t403 + t307 * t425 + t313 * t424 - t357 * t156 + (t269 * t362 + t272 * t358 - t359) * t155 + (t269 * t4 - t272 * t3) * mrSges(6,3) + (-t1 * t272 - t2 * t269) * mrSges(7,2) + (m(4) * t57 - m(5) * t49 - t132 + t135) * pkin(9) + (t314 * t416 + t308 * t417 + t63 * t404 - t77 * t321 - t27 * t319 + (t19 * t272 + t20 * t269) * mrSges(6,3) + (-t11 * t272 + t12 * t269) * mrSges(7,2) + t459 * t415 + t460 * t409 + t443 * t403) * qJD(5) - t434 * t75 / 0.2e1 + (qJD(5) * t60 + t462) * t405) * t273 + t248 + t438 * t219 + (-t58 * mrSges(4,3) + t54 * mrSges(5,1) + t208 * mrSges(4,2) + t83 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1 - t53 * mrSges(5,3) + t423 - t358 * t76 + t362 * t75 - t359 * t156 + (t356 - t363) * t155 + (-m(4) * t58 - t134 + t439) * pkin(9) + t286) * t270 - t442 + t437 * t45 + ((t385 - t381 / 0.2e1 + t122 * mrSges(4,1) - t108 * mrSges(5,2) + t123 * mrSges(4,2) + t110 * mrSges(5,3) - t117 / 0.2e1 - t120 / 0.2e1 + t173 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t250 - t361 * t199 + t360 * t198 + (t344 / 0.2e1 + (pkin(1) * mrSges(3,1) + t395 / 0.2e1) * t268) * qJD(1) + (-t450 * t270 - t448 * t273) * qJD(2) / 0.2e1) * t271 + (-t174 / 0.2e1 - t255 / 0.2e1 + t386 - qJD(2) * Ifges(3,5) / 0.2e1 + (-t345 / 0.2e1 + (pkin(1) * mrSges(3,2) - t384 / 0.2e1 + t387 / 0.2e1) * t268) * qJD(1) + (-t441 + t472) * t270 - t461) * t274) * t370 + (-Ifges(6,2) * t416 + Ifges(7,3) * t417 - t409 * t445 - t415 * t447 + t463) * t184 + t468 * t112 + (t164 * t4 + t19 * t467 + t20 * t468 + t252 * t29 + t3 * t437 + t465 * t77) * m(6) + t469 * t114 + t470 * t111 + (t1 * t157 + t11 * t469 + t12 * t470 + t162 * t2 + t192 * t9 + t27 * t466) * m(7) + t467 * t113 + t466 * t89 + t465 * t90 + (t461 + (t441 - t464) * t270 + (-t371 * t273 - t372 * t270 + m(5) * (t108 * t273 + t110 * t270) + m(4) * (t122 * t273 - t123 * t270)) * pkin(9)) * qJD(3) - m(4) * (-t122 * t141 + t123 * t142 + t177 * t219) - m(5) * (t104 * t143 + t108 * t131 + t110 * t129) + t252 * t22 + t242 * t95 + (t223 - t143) * t138 - t208 * mrSges(3,1) - t216 * t215 + t192 * t21 - t129 * t160 - t131 * t161 + t162 * t44 + t164 * t43 + t157 * t46 - t142 * t158 - t141 * t159 + m(5) * (t104 * t223 + t242 * t53) + (-m(4) * t208 - t96) * pkin(2) + (-mrSges(6,2) * t77 - t11 * mrSges(7,2) + t19 * mrSges(6,3) + mrSges(7,3) * t27 + Ifges(6,4) * t416 + Ifges(7,5) * t417 + t409 * t449 + t415 * t451 - t443 / 0.2e1) * t185; t462 * t402 + (t277 - t190 / 0.2e1 - t191 / 0.2e1) * t198 + (-t132 + t22) * qJ(4) - m(7) * (t11 * t24 + t12 * t23 + t27 * t47) + m(6) * (t29 * qJ(4) + t77 * qJD(4) - t3 * t376 - t375 * t4) + m(7) * (-t1 * t376 + t2 * t375 + t27 * t213 + t9 * t240) - m(6) * (t19 * t32 + t20 * t33 - t291 * t77) + t291 * t90 - (t269 * t399 + t272 * t400) * t422 + (-pkin(3) * t54 - qJ(4) * t49 - t104 * t136 - t108 * t123 + t110 * t458) * m(5) + (t213 - t47) * t89 + t432 + (t313 * t416 + t307 * t417 - (-m(6) * t302 + m(7) * t303 + t269 * t373 + t272 * t374) * t422 + t434 * t415 + t435 * t409 + t471) * qJD(5) + ((-t357 + t363) * t198 + t464) * t199 + t383 * qJD(4) + t371 * t123 + t372 * t122 + t240 * t21 - pkin(3) * t133 - t136 * t138 - t23 * t111 - t33 * t112 - t32 * t113 - t24 * t114 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t54 * mrSges(5,2) - t49 * mrSges(5,3) + t9 * t319 + t29 * t321 + t323 * mrSges(6,3) + t324 * mrSges(7,2) + t13 * t404 + t16 * t405 + t308 * t424 + t314 * t425 + t459 * t426 + t460 * t413; t199 * t138 + (t89 + t383) * t250 + (t193 * t374 + t400) * t272 + (t193 * t373 + t399) * t269 - m(5) * (-t104 * t199 + t110 * t250) + (t193 * t303 + t250 * t27 - t324) * m(7) + (-t193 * t302 + t250 * t77 - t323) * m(6) + t439; t286 + (-t373 + t396) * t20 + (-t374 + t397) * t19 - t77 * (mrSges(6,1) * t148 + mrSges(6,2) * t296) + t63 * t414 + (Ifges(7,3) * t148 + t390) * t417 - t27 * (mrSges(7,1) * t148 - mrSges(7,3) * t296) + qJD(6) * t111 - t88 * t89 - pkin(5) * t44 + qJ(6) * t46 + (-t11 * t296 + t12 * t148) * mrSges(7,2) + (-t148 * t445 + t296 * t449) * t409 + (-pkin(5) * t2 + qJ(6) * t1 - t11 * t20 + t12 * t440 - t27 * t88) * m(7) + (-Ifges(6,2) * t148 + t145 + t443) * t416 + (t296 * t451 + t144 - t393 + t60) * t415 + t444; -t193 * t111 + t148 * t89 + 0.2e1 * (t2 / 0.2e1 + t12 * t409 + t27 * t414) * m(7) + t44;];
tauc  = t17(:);
