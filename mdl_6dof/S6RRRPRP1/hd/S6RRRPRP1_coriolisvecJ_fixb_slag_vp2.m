% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:41:10
% EndTime: 2018-11-23 17:41:20
% DurationCPUTime: 10.16s
% Computational Cost: add. (13869->604), mult. (36400->785), div. (0->0), fcn. (26783->8), ass. (0->288)
t471 = Ifges(6,4) + Ifges(7,4);
t472 = Ifges(6,1) + Ifges(7,1);
t465 = Ifges(6,5) + Ifges(7,5);
t470 = Ifges(6,2) + Ifges(7,2);
t464 = Ifges(7,6) + Ifges(6,6);
t303 = cos(qJ(5));
t475 = t471 * t303;
t300 = sin(qJ(5));
t474 = t471 * t300;
t297 = qJD(2) + qJD(3);
t473 = t297 / 0.2e1;
t455 = -t464 * t300 + t465 * t303;
t453 = -t470 * t300 + t475;
t451 = t472 * t303 - t474;
t301 = sin(qJ(3));
t302 = sin(qJ(2));
t304 = cos(qJ(3));
t305 = cos(qJ(2));
t378 = t304 * t305;
t269 = -t301 * t302 + t378;
t257 = t269 * qJD(1);
t270 = t301 * t305 + t304 * t302;
t258 = t270 * qJD(1);
t298 = sin(pkin(10));
t299 = cos(pkin(10));
t320 = t257 * t298 + t299 * t258;
t190 = t297 * t303 - t300 * t320;
t469 = t471 * t190;
t191 = t297 * t300 + t303 * t320;
t468 = t471 * t191;
t437 = -pkin(8) - pkin(7);
t284 = t437 * t305;
t275 = qJD(1) * t284;
t259 = t301 * t275;
t283 = t437 * t302;
t274 = qJD(1) * t283;
t265 = qJD(2) * pkin(2) + t274;
t218 = t304 * t265 + t259;
t252 = t258 * qJ(4);
t185 = t218 - t252;
t174 = pkin(3) * t297 + t185;
t262 = t304 * t275;
t219 = t265 * t301 - t262;
t392 = qJ(4) * t257;
t186 = t219 + t392;
t177 = t298 * t186;
t107 = t174 * t299 - t177;
t291 = -pkin(2) * t305 - pkin(1);
t282 = qJD(1) * t291;
t227 = -t257 * pkin(3) + qJD(4) + t282;
t467 = t227 * mrSges(5,2) - t107 * mrSges(5,3);
t349 = t299 * t257 - t258 * t298;
t398 = t320 * Ifges(5,4);
t466 = Ifges(5,6) * t473 + t398 / 0.2e1 + t349 * Ifges(5,2) / 0.2e1;
t225 = t297 * t269;
t214 = t225 * qJD(1);
t226 = t297 * t270;
t215 = t226 * qJD(1);
t158 = t214 * t298 + t299 * t215;
t159 = t214 * t299 - t215 * t298;
t95 = qJD(5) * t190 + t159 * t303;
t96 = -qJD(5) * t191 - t159 * t300;
t463 = t158 * t464 + t470 * t96 + t471 * t95;
t462 = t158 * t465 + t471 * t96 + t472 * t95;
t201 = qJD(5) - t349;
t461 = t190 * t470 + t201 * t464 + t468;
t447 = t191 * t472 + t201 * t465 + t469;
t383 = t298 * t301;
t406 = pkin(2) * qJD(3);
t249 = (t299 * t304 - t383) * t406;
t223 = -t274 * t301 + t262;
t192 = t223 - t392;
t224 = t304 * t274 + t259;
t193 = -t252 + t224;
t122 = t192 * t298 + t193 * t299;
t417 = pkin(3) * t258;
t128 = pkin(4) * t320 - pkin(9) * t349 + t417;
t372 = qJD(1) * t302;
t293 = pkin(2) * t372;
t123 = t128 + t293;
t50 = t303 * t122 + t300 * t123;
t459 = t249 * t303 - t50;
t381 = t299 * t301;
t441 = -t299 * t192 + t193 * t298 - (t298 * t304 + t381) * t406;
t388 = t349 * t300;
t458 = qJ(6) * t388 + t303 * qJD(6);
t105 = -pkin(4) * t297 - t107;
t336 = mrSges(7,1) * t300 + mrSges(7,2) * t303;
t338 = mrSges(6,1) * t300 + mrSges(6,2) * t303;
t70 = -pkin(5) * t190 + qJD(6) + t105;
t457 = t105 * t338 + t70 * t336;
t456 = t300 * t465 + t303 * t464;
t454 = t303 * t470 + t474;
t452 = t300 * t472 + t475;
t430 = t191 / 0.2e1;
t450 = t457 + t453 * t190 / 0.2e1 + t451 * t430 + t455 * t201 / 0.2e1;
t382 = t299 * t186;
t108 = t298 * t174 + t382;
t106 = pkin(9) * t297 + t108;
t115 = -pkin(4) * t349 - pkin(9) * t320 + t227;
t37 = -t106 * t300 + t303 * t115;
t25 = -qJ(6) * t191 + t37;
t22 = pkin(5) * t201 + t25;
t38 = t106 * t303 + t115 * t300;
t26 = qJ(6) * t190 + t38;
t449 = -t227 * mrSges(5,1) - t37 * mrSges(6,1) - t22 * mrSges(7,1) + t38 * mrSges(6,2) + t26 * mrSges(7,2) + t466;
t427 = -t320 / 0.2e1;
t448 = -t349 / 0.2e1;
t401 = t349 * Ifges(5,4);
t296 = t303 * qJ(6);
t342 = pkin(5) * t320 - t296 * t349;
t288 = pkin(3) * t298 + pkin(9);
t376 = -qJ(6) - t288;
t347 = qJD(5) * t376;
t113 = t185 * t299 - t177;
t47 = -t113 * t300 + t303 * t128;
t446 = -qJD(6) * t300 + t303 * t347 - t342 - t47;
t290 = pkin(2) * t304 + pkin(3);
t251 = pkin(2) * t381 + t298 * t290;
t244 = pkin(9) + t251;
t377 = -qJ(6) - t244;
t348 = qJD(5) * t377;
t49 = -t122 * t300 + t303 * t123;
t445 = -t342 - t49 + (-qJD(6) - t249) * t300 + t303 * t348;
t48 = t303 * t113 + t300 * t128;
t444 = t300 * t347 + t458 - t48;
t443 = t300 * t348 + t458 + t459;
t198 = pkin(5) * t388;
t367 = qJD(5) * t300;
t364 = pkin(5) * t367;
t442 = -t198 + t364 - t441;
t440 = -t122 + t249;
t230 = t304 * t283 + t284 * t301;
t204 = -qJ(4) * t270 + t230;
t271 = t301 * t283;
t231 = -t304 * t284 + t271;
t205 = qJ(4) * t269 + t231;
t141 = t204 * t298 + t205 * t299;
t136 = t303 * t141;
t221 = -t299 * t269 + t270 * t298;
t222 = t269 * t298 + t270 * t299;
t239 = -t269 * pkin(3) + t291;
t139 = t221 * pkin(4) - t222 * pkin(9) + t239;
t62 = t300 * t139 + t136;
t439 = t95 / 0.2e1;
t438 = t96 / 0.2e1;
t436 = pkin(1) * mrSges(3,1);
t435 = pkin(1) * mrSges(3,2);
t434 = t158 / 0.2e1;
t433 = -t190 / 0.2e1;
t431 = -t191 / 0.2e1;
t429 = -t201 / 0.2e1;
t425 = t257 / 0.2e1;
t424 = t258 / 0.2e1;
t423 = -t297 / 0.2e1;
t422 = -t300 / 0.2e1;
t419 = t303 / 0.2e1;
t418 = m(4) * t282;
t416 = pkin(5) * t303;
t308 = (t378 * t437 - t271) * qJD(2) * qJD(1);
t319 = -t214 * qJ(4) - t258 * qJD(4);
t368 = qJD(3) * t304;
t369 = qJD(3) * t301;
t359 = qJD(2) * t437;
t161 = t359 * t258 + t265 * t368 + t275 * t369;
t99 = -qJ(4) * t215 + qJD(4) * t257 + t161;
t33 = t299 * t99 + (-t265 * t369 + t275 * t368 + t308 + t319) * t298;
t366 = qJD(5) * t303;
t196 = pkin(3) * t215 + qJD(2) * t293;
t53 = pkin(4) * t158 - pkin(9) * t159 + t196;
t7 = -t106 * t367 + t115 * t366 + t300 * t53 + t303 * t33;
t415 = t303 * t7;
t414 = mrSges(4,3) * t257;
t413 = mrSges(4,3) * t258;
t412 = Ifges(3,4) * t302;
t411 = Ifges(4,4) * t258;
t140 = -t299 * t204 + t205 * t298;
t162 = -qJD(3) * t219 + t308;
t32 = t298 * t99 - t299 * (t162 + t319);
t404 = t140 * t32;
t399 = t320 * Ifges(5,1);
t397 = t297 * Ifges(5,5);
t395 = t300 * t37;
t394 = Ifges(3,5) * qJD(2);
t393 = Ifges(3,6) * qJD(2);
t391 = qJD(2) * mrSges(3,1);
t390 = qJD(2) * mrSges(3,2);
t389 = t108 * t320;
t387 = t349 * t303;
t386 = t222 * t300;
t384 = t288 * t303;
t375 = -mrSges(5,1) * t297 - mrSges(6,1) * t190 + mrSges(6,2) * t191 + mrSges(5,3) * t320;
t371 = qJD(1) * t305;
t370 = qJD(2) * t302;
t276 = t302 * t359;
t277 = t305 * t359;
t170 = t304 * t276 + t301 * t277 + t283 * t368 + t284 * t369;
t119 = -qJ(4) * t226 + qJD(4) * t269 + t170;
t171 = -qJD(3) * t231 - t276 * t301 + t304 * t277;
t120 = -qJ(4) * t225 - qJD(4) * t270 + t171;
t42 = t119 * t299 + t120 * t298;
t163 = t225 * t298 + t299 * t226;
t164 = t225 * t299 - t226 * t298;
t211 = pkin(2) * t370 + pkin(3) * t226;
t69 = pkin(4) * t163 - pkin(9) * t164 + t211;
t365 = t139 * t366 + t300 * t69 + t303 * t42;
t363 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t362 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t361 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t289 = -pkin(3) * t299 - pkin(4);
t358 = t222 * t366;
t357 = t394 / 0.2e1;
t356 = -t393 / 0.2e1;
t27 = -t96 * mrSges(7,1) + t95 * mrSges(7,2);
t351 = -t300 * t42 + t303 * t69;
t350 = t158 * mrSges(5,1) + t159 * mrSges(5,2);
t41 = t119 * t298 - t299 * t120;
t61 = t303 * t139 - t141 * t300;
t112 = t185 * t298 + t382;
t250 = -pkin(2) * t383 + t290 * t299;
t8 = -qJD(5) * t38 - t300 * t33 + t303 * t53;
t1 = pkin(5) * t158 - qJ(6) * t95 - qJD(6) * t191 + t8;
t345 = -t8 * mrSges(6,3) - t1 * mrSges(7,3);
t243 = -pkin(4) - t250;
t344 = -mrSges(6,3) * t37 - mrSges(7,3) * t22;
t343 = -mrSges(6,3) * t38 - t26 * mrSges(7,3);
t3 = qJ(6) * t96 + qJD(6) * t190 + t7;
t341 = -t1 * t303 - t3 * t300;
t340 = -t300 * t7 - t303 * t8;
t339 = mrSges(6,1) * t303 - mrSges(6,2) * t300;
t337 = mrSges(7,1) * t303 - mrSges(7,2) * t300;
t323 = t22 * t300 - t26 * t303;
t322 = -t300 * t38 - t303 * t37;
t321 = -t303 * t38 + t395;
t318 = -qJ(6) * t164 - qJD(6) * t222;
t309 = t8 * mrSges(6,1) + t1 * mrSges(7,1) - t7 * mrSges(6,2) - t3 * mrSges(7,2);
t307 = m(6) * (qJD(5) * t322 - t8 * t300 + t415);
t12 = -pkin(5) * t96 + t32;
t144 = t397 + t399 + t401;
t202 = t257 * Ifges(4,2) + t297 * Ifges(4,6) + t411;
t253 = Ifges(4,4) * t257;
t203 = Ifges(4,1) * t258 + t297 * Ifges(4,5) + t253;
t84 = t191 * Ifges(7,5) + t190 * Ifges(7,6) + t201 * Ifges(7,3);
t85 = t191 * Ifges(6,5) + t190 * Ifges(6,6) + t201 * Ifges(6,3);
t306 = (-t398 + t85 + t84) * t427 + (Ifges(5,1) * t427 + Ifges(5,5) * t423 + t429 * t455 + t431 * t451 + t433 * t453 - t457 - t467) * t349 + (-mrSges(5,1) - t339) * t32 - t282 * (mrSges(4,1) * t258 + mrSges(4,2) * t257) + (-t367 / 0.2e1 + t388 / 0.2e1) * t461 + (t366 / 0.2e1 - t387 / 0.2e1) * t447 + (t37 * t387 + t38 * t388 + t415) * mrSges(6,3) - t33 * mrSges(5,2) + t452 * t439 + t454 * t438 + t456 * t434 + (t22 * t387 + t26 * t388 + t3 * t303) * mrSges(7,3) + t450 * qJD(5) - (-Ifges(4,2) * t258 + t203 + t253) * t257 / 0.2e1 + (t401 + t144) * t448 + (Ifges(4,5) * t257 - Ifges(4,6) * t258) * t423 + t202 * t424 + t218 * t414 + t462 * t300 / 0.2e1 + t463 * t419 + (-Ifges(5,6) * t423 - Ifges(5,2) * t448 + t464 * t433 + t465 * t431 + (Ifges(7,3) + Ifges(6,3)) * t429 + t449) * t320 - Ifges(5,6) * t158 - t12 * t337 + Ifges(5,5) * t159 - t161 * mrSges(4,2) + t162 * mrSges(4,1) + Ifges(4,5) * t214 - Ifges(4,6) * t215 - t258 * (Ifges(4,1) * t257 - t411) / 0.2e1;
t292 = Ifges(3,4) * t371;
t280 = mrSges(3,3) * t371 - t390;
t279 = -mrSges(3,3) * t372 + t391;
t278 = t289 - t416;
t267 = t296 + t384;
t266 = t376 * t300;
t256 = Ifges(3,1) * t372 + t292 + t394;
t255 = t393 + (t305 * Ifges(3,2) + t412) * qJD(1);
t236 = t243 - t416;
t235 = mrSges(4,1) * t297 - t413;
t234 = -mrSges(4,2) * t297 + t414;
t233 = t293 + t417;
t229 = t244 * t303 + t296;
t228 = t377 * t300;
t217 = -mrSges(4,1) * t257 + mrSges(4,2) * t258;
t194 = -mrSges(5,2) * t297 + mrSges(5,3) * t349;
t155 = Ifges(6,3) * t158;
t154 = Ifges(7,3) * t158;
t152 = -mrSges(5,1) * t349 + mrSges(5,2) * t320;
t133 = mrSges(6,1) * t201 - mrSges(6,3) * t191;
t132 = mrSges(7,1) * t201 - mrSges(7,3) * t191;
t131 = -mrSges(6,2) * t201 + mrSges(6,3) * t190;
t130 = -mrSges(7,2) * t201 + mrSges(7,3) * t190;
t124 = -mrSges(7,1) * t190 + mrSges(7,2) * t191;
t102 = pkin(5) * t386 + t140;
t94 = Ifges(6,5) * t95;
t93 = Ifges(7,5) * t95;
t92 = Ifges(6,6) * t96;
t91 = Ifges(7,6) * t96;
t74 = t112 + t198;
t46 = -mrSges(6,2) * t158 + mrSges(6,3) * t96;
t45 = -mrSges(7,2) * t158 + mrSges(7,3) * t96;
t44 = mrSges(6,1) * t158 - mrSges(6,3) * t95;
t43 = mrSges(7,1) * t158 - mrSges(7,3) * t95;
t39 = -qJ(6) * t386 + t62;
t36 = pkin(5) * t221 - t222 * t296 + t61;
t28 = -mrSges(6,1) * t96 + mrSges(6,2) * t95;
t21 = (t164 * t300 + t358) * pkin(5) + t41;
t10 = -qJD(5) * t62 + t351;
t9 = -t141 * t367 + t365;
t5 = -qJ(6) * t358 + (-qJD(5) * t141 + t318) * t300 + t365;
t4 = pkin(5) * t163 + t318 * t303 + (-t136 + (qJ(6) * t222 - t139) * t300) * qJD(5) + t351;
t2 = [(t84 / 0.2e1 + t85 / 0.2e1 - t108 * mrSges(5,3) + t361 * t201 + t363 * t191 + t362 * t190 - t449 - t466) * t163 + (t461 * t422 + t397 / 0.2e1 + t401 / 0.2e1 + t399 / 0.2e1 + (-t22 * t303 - t26 * t300) * mrSges(7,3) + t144 / 0.2e1 + t447 * t419 + t322 * mrSges(6,3) + t450 + t467) * t164 + t61 * t44 + t62 * t46 + t36 * t43 + t39 * t45 + (t161 * t269 - t162 * t270 - t214 * t230 - t215 * t231 - t218 * t225 - t219 * t226) * mrSges(4,3) + (-t269 * t215 - t226 * t425) * Ifges(4,2) + (t269 * t214 - t270 * t215 + t225 * t425 - t226 * t424) * Ifges(4,4) + t291 * (mrSges(4,1) * t215 + mrSges(4,2) * t214) + m(4) * (t161 * t231 + t162 * t230 + t170 * t219 + t171 * t218) + m(7) * (t1 * t36 + t102 * t12 + t21 * t70 + t22 * t4 + t26 * t5 + t3 * t39) + t282 * (mrSges(4,1) * t226 + mrSges(4,2) * t225) + (t196 * mrSges(5,2) + t12 * t336 + Ifges(5,1) * t159 - Ifges(5,4) * t158 + (mrSges(5,3) + t338) * t32 + t341 * mrSges(7,3) + t340 * mrSges(6,3) + (mrSges(6,3) * t321 + mrSges(7,3) * t323 + t105 * t339 + t337 * t70 + t454 * t433 + t452 * t431 + t456 * t429 - t461 * t303 / 0.2e1) * qJD(5) + t451 * t439 + t453 * t438 + t455 * t434 + t462 * t419 + (qJD(5) * t447 + t463) * t422) * t222 + (t256 / 0.2e1 - pkin(7) * t279 + t357 + (-0.2e1 * t435 + 0.3e1 / 0.2e1 * Ifges(3,4) * t305) * qJD(1)) * t305 * qJD(2) + t102 * t27 + t21 * t124 + t5 * t130 + t9 * t131 + t4 * t132 + t10 * t133 + t140 * t28 + (Ifges(4,5) * t225 - Ifges(4,6) * t226) * t473 + t42 * t194 + t211 * t152 + t225 * t203 / 0.2e1 - t226 * t202 / 0.2e1 + t170 * t234 + t171 * t235 + t239 * t350 + (-t33 * mrSges(5,3) + t93 / 0.2e1 + t91 / 0.2e1 + t154 / 0.2e1 + t94 / 0.2e1 + t92 / 0.2e1 + t155 / 0.2e1 - Ifges(5,4) * t159 + t196 * mrSges(5,1) + t362 * t96 + t363 * t95 + (Ifges(5,2) + t361) * t158 + t309) * t221 + t375 * t41 + (t140 * t159 - t141 * t158) * mrSges(5,3) + m(5) * (-t107 * t41 + t108 * t42 + t141 * t33 + t196 * t239 + t211 * t227 + t404) + m(6) * (t10 * t37 + t105 * t41 + t38 * t9 + t61 * t8 + t62 * t7 + t404) + (t270 * t214 + t225 * t424) * Ifges(4,1) + (-t255 / 0.2e1 - pkin(7) * t280 + t356 + (-0.2e1 * t436 - 0.3e1 / 0.2e1 * t412 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t305) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t269 + mrSges(4,2) * t270) + 0.2e1 * t418 + t217) * pkin(2)) * t370; t244 * t307 + t443 * t130 + (t1 * t228 + t12 * t236 + t22 * t445 + t229 * t3 + t26 * t443 + t442 * t70) * m(7) + t445 * t132 + t442 * t124 + ((t234 * t304 - t235 * t301) * qJD(3) + (-t214 * t304 - t215 * t301) * mrSges(4,3)) * pkin(2) + t440 * t194 + t306 + ((t161 * t301 + t162 * t304 + (-t218 * t301 + t219 * t304) * qJD(3)) * pkin(2) - t218 * t223 - t219 * t224) * m(4) + t219 * t413 + (t107 * t441 + t108 * t440 - t227 * t233 - t250 * t32 + t251 * t33) * m(5) - t441 * t375 + (-t441 * t105 + t243 * t32 - t249 * t395 - t37 * t49 + t38 * t459) * m(6) - t50 * t131 - t49 * t133 + (t249 * t131 + t244 * t46 + (-t133 * t244 + t344) * qJD(5)) * t303 + (-t249 * t133 - t244 * t44 + (-t131 * t244 + t343) * qJD(5) + t345) * t300 + t228 * t43 + t229 * t45 - t233 * t152 - t224 * t234 - t223 * t235 + t236 * t27 + t243 * t28 + (-t158 * t251 - t159 * t250 + t389) * mrSges(5,3) + ((-t292 / 0.2e1 - t256 / 0.2e1 + t357 + qJD(1) * t435 + (t279 - t391) * pkin(7)) * t305 + (t255 / 0.2e1 + t356 + (t436 + t412 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t305) * qJD(1) + (t280 + t390) * pkin(7) + (-t217 - t418) * pkin(2)) * t302) * qJD(1); t288 * t307 + t46 * t384 + t444 * t130 + t306 - t74 * t124 - t48 * t131 - t47 * t133 + ((-t133 * t288 + t344) * t303 + (pkin(5) * t124 - t131 * t288 + t343) * t300) * qJD(5) + (-t288 * t44 + t345) * t300 + t446 * t132 - t113 * t194 - t218 * t234 + t266 * t43 + t267 * t45 + t278 * t27 + t289 * t28 - t375 * t112 + (t389 + (-t158 * t298 - t159 * t299) * pkin(3)) * mrSges(5,3) + (t235 + t413) * t219 - t152 * t417 + (t1 * t266 + t12 * t278 + t267 * t3 + (t364 - t74) * t70 + t444 * t26 + t446 * t22) * m(7) + (-t105 * t112 + t289 * t32 - t37 * t47 - t38 * t48) * m(6) + ((t298 * t33 - t299 * t32) * pkin(3) + t107 * t112 - t108 * t113 - t227 * t417) * m(5); -t349 * t194 + (-t124 - t375) * t320 + (t43 + t44 + t201 * (t130 + t131)) * t303 + (t45 + t46 - t201 * (t132 + t133)) * t300 + t350 + (-t201 * t323 - t320 * t70 - t341) * m(7) + (-t105 * t320 - t201 * t321 - t340) * m(6) + (t107 * t320 - t108 * t349 + t196) * m(5); t155 + t154 + t309 + (-t124 * t191 + t43) * pkin(5) + (-(-t22 + t25) * t26 + (-t191 * t70 + t1) * pkin(5)) * m(7) - t25 * t130 - t37 * t131 + t26 * t132 + t38 * t133 - t70 * (mrSges(7,1) * t191 + mrSges(7,2) * t190) - t105 * (mrSges(6,1) * t191 + mrSges(6,2) * t190) + t92 + t93 + t94 + (t190 * t37 + t191 * t38) * mrSges(6,3) + (t190 * t22 + t191 * t26) * mrSges(7,3) + t91 + (t190 * t472 - t468) * t431 + t461 * t430 + (t190 * t465 - t464 * t191) * t429 + (-t191 * t470 + t447 + t469) * t433; -t190 * t130 + t191 * t132 + 0.2e1 * (t12 / 0.2e1 + t26 * t433 + t22 * t430) * m(7) + t27;];
tauc  = t2(:);
