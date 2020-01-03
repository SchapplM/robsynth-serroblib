% Calculate vector of inverse dynamics joint torques for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:07
% EndTime: 2019-12-31 21:16:33
% DurationCPUTime: 13.18s
% Computational Cost: add. (8791->639), mult. (19494->859), div. (0->0), fcn. (13914->14), ass. (0->304)
t482 = -mrSges(5,3) - mrSges(6,3);
t286 = qJ(2) + qJ(3);
t280 = sin(t286);
t287 = sin(pkin(9));
t386 = t287 * mrSges(5,2);
t284 = pkin(9) + qJ(5);
t278 = sin(t284);
t398 = mrSges(6,2) * t278;
t481 = (-t386 - t398) * t280;
t288 = cos(pkin(9));
t400 = t288 * pkin(4);
t270 = pkin(3) + t400;
t281 = cos(t286);
t289 = -pkin(8) - qJ(4);
t452 = t281 * t270 - t280 * t289;
t480 = m(6) * t452;
t290 = sin(qJ(5));
t294 = cos(qJ(5));
t291 = sin(qJ(3));
t292 = sin(qJ(2));
t295 = cos(qJ(2));
t410 = cos(qJ(3));
t233 = t291 * t295 + t292 * t410;
t218 = t233 * qJD(1);
t285 = qJD(2) + qJD(3);
t194 = t218 * t288 + t285 * t287;
t341 = t410 * t295;
t361 = qJD(1) * t292;
t217 = -qJD(1) * t341 + t291 * t361;
t407 = pkin(2) * t295;
t274 = pkin(1) + t407;
t251 = t274 * qJD(1);
t153 = pkin(3) * t217 - qJ(4) * t218 - t251;
t297 = -pkin(7) - pkin(6);
t254 = t297 * t295;
t236 = qJD(1) * t254;
t220 = t410 * t236;
t253 = t297 * t292;
t235 = qJD(1) * t253;
t225 = qJD(2) * pkin(2) + t235;
t180 = t291 * t225 - t220;
t162 = t285 * qJ(4) + t180;
t86 = t288 * t153 - t162 * t287;
t56 = pkin(4) * t217 - pkin(8) * t194 + t86;
t333 = -t218 * t287 + t288 * t285;
t87 = t287 * t153 + t288 * t162;
t61 = pkin(8) * t333 + t87;
t14 = -t290 * t61 + t294 * t56;
t479 = t14 * mrSges(6,1);
t15 = t290 * t56 + t294 * t61;
t478 = t15 * mrSges(6,2);
t477 = t86 * mrSges(5,1);
t476 = t87 * mrSges(5,2);
t119 = t194 * t294 + t290 * t333;
t211 = qJD(5) + t217;
t462 = t333 * Ifges(5,6);
t466 = t194 * Ifges(5,5);
t469 = -t194 * t290 + t294 * t333;
t475 = t119 * Ifges(6,5) + Ifges(6,6) * t469 + t217 * Ifges(5,3) + t211 * Ifges(6,3) + t462 + t466;
t231 = t287 * t294 + t288 * t290;
t151 = t231 * t217;
t214 = t231 * qJD(5);
t474 = t151 + t214;
t313 = t287 * t290 - t288 * t294;
t152 = t313 * t217;
t213 = t313 * qJD(5);
t473 = t152 + t213;
t472 = -t281 * mrSges(4,1) + (mrSges(4,2) - mrSges(6,3)) * t280;
t357 = qJD(1) * qJD(2);
t241 = qJDD(1) * t295 - t292 * t357;
t384 = t288 * mrSges(5,1);
t279 = cos(t284);
t399 = mrSges(6,1) * t279;
t471 = (t384 + t399) * t280;
t470 = t384 - t386;
t242 = qJDD(1) * t292 + t295 * t357;
t309 = -t291 * t292 + t341;
t302 = t309 * qJD(3);
t141 = qJD(1) * t302 + t291 * t241 + t242 * t410;
t303 = t233 * qJD(3);
t142 = qJD(1) * t303 - t410 * t241 + t291 * t242;
t382 = qJDD(1) * pkin(1);
t208 = -pkin(2) * t241 - t382;
t57 = pkin(3) * t142 - qJ(4) * t141 - qJD(4) * t218 + t208;
t283 = qJDD(2) + qJDD(3);
t229 = t242 * pkin(6);
t188 = qJDD(2) * pkin(2) - pkin(7) * t242 - t229;
t228 = t241 * pkin(6);
t196 = pkin(7) * t241 + t228;
t338 = qJD(3) * t410;
t358 = qJD(3) * t291;
t74 = t291 * t188 + t410 * t196 + t225 * t338 + t236 * t358;
t70 = qJ(4) * t283 + qJD(4) * t285 + t74;
t17 = -t287 * t70 + t288 * t57;
t18 = t287 * t57 + t288 * t70;
t317 = -t17 * t287 + t18 * t288;
t122 = -t141 * t287 + t283 * t288;
t123 = t141 * t288 + t283 * t287;
t36 = qJD(5) * t469 + t122 * t290 + t123 * t294;
t433 = t36 / 0.2e1;
t37 = -qJD(5) * t119 + t122 * t294 - t123 * t290;
t432 = t37 / 0.2e1;
t424 = t122 / 0.2e1;
t423 = t123 / 0.2e1;
t140 = qJDD(5) + t142;
t422 = t140 / 0.2e1;
t421 = t142 / 0.2e1;
t467 = t241 / 0.2e1;
t411 = t295 / 0.2e1;
t465 = t285 * Ifges(4,5);
t464 = t285 * Ifges(4,6);
t463 = t295 * Ifges(3,2);
t409 = pkin(2) * t291;
t269 = qJ(4) + t409;
t222 = (-pkin(8) - t269) * t287;
t282 = t288 * pkin(8);
t373 = t269 * t288;
t223 = t282 + t373;
t173 = t222 * t290 + t223 * t294;
t330 = pkin(2) * t338;
t262 = t330 + qJD(4);
t378 = t217 * t288;
t329 = t218 * pkin(4) + pkin(8) * t378;
t174 = pkin(3) * t218 + qJ(4) * t217;
t348 = pkin(2) * t361;
t158 = t174 + t348;
t219 = t291 * t236;
t183 = t235 * t410 + t219;
t97 = t288 * t158 - t183 * t287;
t63 = t329 + t97;
t379 = t217 * t287;
t355 = pkin(8) * t379;
t98 = t287 * t158 + t288 * t183;
t76 = t355 + t98;
t461 = -qJD(5) * t173 - t231 * t262 + t290 * t76 - t294 * t63;
t172 = t222 * t294 - t223 * t290;
t460 = qJD(5) * t172 - t262 * t313 - t290 * t63 - t294 * t76;
t246 = t289 * t287;
t383 = qJ(4) * t288;
t247 = t282 + t383;
t191 = t246 * t290 + t247 * t294;
t179 = t225 * t410 + t219;
t99 = t288 * t174 - t179 * t287;
t67 = t329 + t99;
t100 = t287 * t174 + t288 * t179;
t80 = t355 + t100;
t459 = -qJD(4) * t231 - qJD(5) * t191 + t290 * t80 - t294 * t67;
t190 = t246 * t294 - t247 * t290;
t458 = -qJD(4) * t313 + qJD(5) * t190 - t290 * t67 - t294 * t80;
t157 = -t285 * pkin(3) + qJD(4) - t179;
t323 = mrSges(5,1) * t287 + mrSges(5,2) * t288;
t457 = t157 * t323;
t389 = t218 * mrSges(4,3);
t454 = mrSges(4,1) * t285 + mrSges(5,1) * t333 - mrSges(5,2) * t194 - t389;
t453 = t410 * t253 + t291 * t254;
t314 = -t270 * t280 - t281 * t289;
t406 = pkin(3) * t280;
t408 = pkin(2) * t292;
t450 = -m(6) * (t314 - t408) - m(5) * (-t406 - t408) + t471;
t449 = -m(6) * t314 + t471;
t360 = qJD(1) * t295;
t403 = pkin(6) * t295;
t404 = pkin(6) * t292;
t448 = (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t360) * t404 + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t361) * t403;
t296 = cos(qJ(1));
t370 = t281 * t296;
t447 = t481 * t296 + t370 * t482;
t293 = sin(qJ(1));
t371 = t281 * t293;
t446 = t481 * t293 + t371 * t482;
t145 = -mrSges(5,2) * t217 + mrSges(5,3) * t333;
t146 = mrSges(5,1) * t217 - mrSges(5,3) * t194;
t445 = t288 * t145 - t287 * t146;
t444 = t228 * t295 + t229 * t292;
t443 = g(1) * t296 + g(2) * t293;
t267 = t280 * mrSges(5,3);
t442 = -t267 + t472 + (t398 - t399 - t470) * t281;
t60 = -t122 * mrSges(5,1) + t123 * mrSges(5,2);
t75 = t188 * t410 - t291 * t196 - t225 * t358 + t236 * t338;
t71 = -t283 * pkin(3) + qJDD(4) - t75;
t441 = m(5) * t71 + t60;
t440 = -m(3) * pkin(6) + m(5) * t297 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t323;
t11 = pkin(4) * t142 - pkin(8) * t123 + t17;
t12 = pkin(8) * t122 + t18;
t2 = qJD(5) * t14 + t11 * t290 + t12 * t294;
t3 = -qJD(5) * t15 + t11 * t294 - t12 * t290;
t439 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t438 = m(5) * t157 - t454;
t250 = -mrSges(3,1) * t295 + mrSges(3,2) * t292;
t437 = -(m(5) * pkin(3) + t470) * t281 - mrSges(2,1) - m(3) * pkin(1) + t250 + t472;
t435 = Ifges(6,4) * t433 + Ifges(6,2) * t432 + Ifges(6,6) * t422;
t434 = Ifges(6,1) * t433 + Ifges(6,4) * t432 + Ifges(6,5) * t422;
t431 = Ifges(5,1) * t423 + Ifges(5,4) * t424 + Ifges(5,5) * t421;
t393 = Ifges(6,4) * t119;
t52 = Ifges(6,2) * t469 + Ifges(6,6) * t211 + t393;
t430 = t52 / 0.2e1;
t114 = Ifges(6,4) * t469;
t53 = Ifges(6,1) * t119 + Ifges(6,5) * t211 + t114;
t429 = t53 / 0.2e1;
t428 = -t469 / 0.2e1;
t427 = t469 / 0.2e1;
t426 = -t119 / 0.2e1;
t425 = t119 / 0.2e1;
t419 = -t211 / 0.2e1;
t418 = t211 / 0.2e1;
t417 = -t217 / 0.2e1;
t416 = t217 / 0.2e1;
t414 = t218 / 0.2e1;
t412 = t288 / 0.2e1;
t405 = pkin(4) * t287;
t397 = Ifges(3,4) * t292;
t396 = Ifges(3,4) * t295;
t395 = Ifges(5,4) * t287;
t394 = Ifges(5,4) * t288;
t391 = t179 * mrSges(4,3);
t388 = t218 * Ifges(4,4);
t69 = mrSges(5,1) * t142 - mrSges(5,3) * t123;
t385 = t287 * t69;
t342 = qJD(2) * t297;
t237 = t292 * t342;
t238 = t295 * t342;
t124 = qJD(3) * t453 + t410 * t237 + t291 * t238;
t184 = qJD(2) * t309 + t302;
t185 = qJD(2) * t233 + t303;
t359 = qJD(2) * t292;
t347 = pkin(2) * t359;
t93 = pkin(3) * t185 - qJ(4) * t184 - qJD(4) * t233 + t347;
t50 = t288 * t124 + t287 * t93;
t381 = t184 * t287;
t380 = t184 * t288;
t375 = t233 * t287;
t374 = t233 * t288;
t265 = t280 * qJ(4);
t178 = -pkin(3) * t309 - qJ(4) * t233 - t274;
t198 = t291 * t253 - t254 * t410;
t111 = t287 * t178 + t288 * t198;
t362 = t281 * pkin(3) + t265;
t352 = Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t140;
t351 = t410 * pkin(2);
t340 = -t287 * (t194 * Ifges(5,4) + Ifges(5,2) * t333 + t217 * Ifges(5,6)) / 0.2e1;
t339 = (t194 * Ifges(5,1) + Ifges(5,4) * t333 + t217 * Ifges(5,5)) * t412;
t10 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t335 = t357 / 0.2e1;
t49 = -t124 * t287 + t288 * t93;
t110 = t288 * t178 - t198 * t287;
t182 = t235 * t291 - t220;
t273 = -t351 - pkin(3);
t327 = mrSges(3,1) * t292 + mrSges(3,2) * t295;
t325 = mrSges(4,1) * t280 + mrSges(4,2) * t281;
t322 = Ifges(5,1) * t288 - t395;
t321 = t397 + t463;
t320 = -Ifges(5,2) * t287 + t394;
t319 = Ifges(3,5) * t295 - Ifges(3,6) * t292;
t318 = Ifges(5,5) * t288 - Ifges(5,6) * t287;
t316 = -t287 * t86 + t288 * t87;
t79 = -pkin(4) * t309 - pkin(8) * t374 + t110;
t92 = -pkin(8) * t375 + t111;
t32 = -t290 * t92 + t294 * t79;
t33 = t290 * t79 + t294 * t92;
t310 = pkin(1) * t327;
t308 = t292 * (Ifges(3,1) * t295 - t397);
t125 = qJD(3) * t198 + t291 * t237 - t410 * t238;
t107 = -pkin(4) * t333 + t157;
t160 = -t217 * Ifges(4,2) + t388 + t464;
t209 = Ifges(4,4) * t217;
t161 = t218 * Ifges(4,1) - t209 + t465;
t41 = t123 * Ifges(5,4) + t122 * Ifges(5,2) + t142 * Ifges(5,6);
t45 = -t122 * pkin(4) + t71;
t298 = -t218 * t479 - t333 * (Ifges(5,6) * t218 - t217 * t320) / 0.2e1 + (-Ifges(6,5) * t213 - Ifges(6,6) * t214) * t418 + (-Ifges(6,1) * t213 - Ifges(6,4) * t214) * t425 + (-Ifges(6,4) * t213 - Ifges(6,2) * t214) * t427 + t251 * (mrSges(4,1) * t218 - mrSges(4,2) * t217) + t45 * (mrSges(6,1) * t313 + mrSges(6,2) * t231) + (Ifges(6,5) * t231 - Ifges(6,6) * t313) * t422 + (Ifges(6,4) * t231 - Ifges(6,2) * t313) * t432 + (Ifges(6,1) * t231 - Ifges(6,4) * t313) * t433 - t313 * t435 + (-Ifges(4,2) * t218 + t161 - t209) * t416 - t71 * t470 + t217 * t340 + t217 * t339 - t217 * t391 + t180 * t389 + t218 * t476 + t218 * t478 - (-Ifges(4,1) * t217 - t388 + t475) * t218 / 0.2e1 - t218 * t477 + (-t378 * t86 - t379 * t87 + t317) * mrSges(5,3) + (mrSges(6,1) * t474 - mrSges(6,2) * t473) * t107 + (t14 * t473 - t15 * t474 - t2 * t313 - t231 * t3) * mrSges(6,3) - t74 * mrSges(4,2) + t75 * mrSges(4,1) + Ifges(4,5) * t141 - Ifges(4,6) * t142 - t151 * t52 / 0.2e1 - t152 * t53 / 0.2e1 + t217 * t457 - t194 * (Ifges(5,5) * t218 - t217 * t322) / 0.2e1 + Ifges(4,3) * t283 - t285 * (-Ifges(4,5) * t217 - Ifges(4,6) * t218) / 0.2e1 + t41 * t412 + t160 * t414 + (Ifges(5,3) * t218 - t217 * t318) * t417 + (Ifges(6,5) * t152 + Ifges(6,6) * t151 + Ifges(6,3) * t218) * t419 + (Ifges(5,5) * t287 + Ifges(5,6) * t288) * t421 + (Ifges(5,1) * t287 + t394) * t423 + (Ifges(5,2) * t288 + t395) * t424 + (Ifges(6,1) * t152 + Ifges(6,4) * t151 + Ifges(6,5) * t218) * t426 + (Ifges(6,4) * t152 + Ifges(6,2) * t151 + Ifges(6,6) * t218) * t428 - t213 * t429 - t214 * t430 + t287 * t431 + t231 * t434;
t276 = Ifges(3,4) * t360;
t257 = t296 * t274;
t256 = qJ(4) * t370;
t255 = qJ(4) * t371;
t245 = t273 - t400;
t216 = Ifges(3,1) * t361 + Ifges(3,5) * qJD(2) + t276;
t215 = Ifges(3,6) * qJD(2) + qJD(1) * t321;
t205 = pkin(4) * t379;
t204 = t278 * t293 + t279 * t370;
t203 = -t278 * t370 + t279 * t293;
t202 = t278 * t296 - t279 * t371;
t201 = t278 * t371 + t279 * t296;
t199 = -mrSges(4,2) * t285 - mrSges(4,3) * t217;
t176 = mrSges(4,1) * t217 + mrSges(4,2) * t218;
t166 = t313 * t233;
t165 = t231 * t233;
t155 = pkin(4) * t375 - t453;
t137 = t182 - t205;
t132 = -t205 + t180;
t130 = -mrSges(4,2) * t283 - mrSges(4,3) * t142;
t129 = mrSges(4,1) * t283 - mrSges(4,3) * t141;
t96 = mrSges(6,1) * t211 - mrSges(6,3) * t119;
t95 = -mrSges(6,2) * t211 + mrSges(6,3) * t469;
t83 = pkin(4) * t381 + t125;
t68 = -mrSges(5,2) * t142 + mrSges(5,3) * t122;
t66 = -t184 * t231 + t213 * t233;
t65 = -t184 * t313 - t214 * t233;
t59 = -mrSges(6,1) * t469 + mrSges(6,2) * t119;
t38 = -pkin(8) * t381 + t50;
t27 = pkin(4) * t185 - pkin(8) * t380 + t49;
t20 = -mrSges(6,2) * t140 + mrSges(6,3) * t37;
t19 = mrSges(6,1) * t140 - mrSges(6,3) * t36;
t5 = -qJD(5) * t33 + t27 * t294 - t290 * t38;
t4 = qJD(5) * t32 + t27 * t290 + t294 * t38;
t1 = [t242 * t396 / 0.2e1 + (t295 * t396 + t308) * t335 + (-m(5) * t257 - t204 * mrSges(6,1) - t203 * mrSges(6,2) + (-m(6) - m(4)) * (-t293 * t297 + t257) + (-m(6) * t405 + t440) * t293 + (-t480 - (m(5) * qJ(4) + mrSges(5,3)) * t280 + t437) * t296) * g(2) + m(4) * (t124 * t180 + t198 * t74 - t208 * t274 - t251 * t347) + (-t14 * t65 + t15 * t66 - t165 * t2 + t166 * t3) * mrSges(6,3) + (-Ifges(6,5) * t166 - Ifges(6,6) * t165) * t422 + (-Ifges(6,4) * t166 - Ifges(6,2) * t165) * t432 + (-Ifges(6,1) * t166 - Ifges(6,4) * t165) * t433 + t45 * (mrSges(6,1) * t165 - mrSges(6,2) * t166) + (-t17 * t374 - t18 * t375 - t380 * t86 - t381 * t87) * mrSges(5,3) + (t479 - t478 - t160 / 0.2e1 + t462 / 0.2e1 + t466 / 0.2e1 - t251 * mrSges(4,1) - t464 / 0.2e1 - Ifges(4,4) * t414 + Ifges(5,3) * t416 - Ifges(4,2) * t417 + Ifges(6,3) * t418 + Ifges(6,5) * t425 + Ifges(6,6) * t427 - t180 * mrSges(4,3) + t477 - t476 + t475 / 0.2e1) * t185 - (Ifges(5,5) * t123 + Ifges(5,6) * t122 + Ifges(5,3) * t142 + t352) * t309 / 0.2e1 - (t208 * mrSges(4,1) + t17 * mrSges(5,1) - t18 * mrSges(5,2) - t74 * mrSges(4,3) - Ifges(4,4) * t141 + Ifges(5,5) * t423 + Ifges(6,5) * t433 + Ifges(4,2) * t142 - Ifges(4,6) * t283 + Ifges(5,6) * t424 + Ifges(6,6) * t432 + Ifges(5,3) * t421 + Ifges(6,3) * t422 + t439) * t309 + (-t202 * mrSges(6,1) - t201 * mrSges(6,2) + (-m(6) * (-t297 + t405) + m(4) * t297 + t440) * t296 + (-m(6) * (-t274 - t452) + m(4) * t274 - m(5) * (-t274 - t265) + t267 - t437) * t293) * g(1) + (t208 * mrSges(4,2) - t75 * mrSges(4,3) + Ifges(4,1) * t141 - Ifges(4,4) * t142 + Ifges(4,5) * t283 + t318 * t421 + t320 * t424 + t322 * t423 + t323 * t71) * t233 + (t241 * t403 + t242 * t404 + t444) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t444) + (t216 * t411 + t319 * qJD(2) / 0.2e1 - t448) * qJD(2) + (t340 + t339 + t457 + t161 / 0.2e1 + t333 * t320 / 0.2e1 + t194 * t322 / 0.2e1 - t251 * mrSges(4,2) + t465 / 0.2e1 + Ifges(4,1) * t414 + t318 * t416 + Ifges(4,4) * t417 - t391) * t184 + (-mrSges(3,1) * t404 - mrSges(3,2) * t403 + 0.2e1 * Ifges(3,6) * t411) * qJDD(2) + (Ifges(3,1) * t242 + Ifges(3,4) * t467 + Ifges(3,5) * qJDD(2) - t335 * t463) * t292 + (-m(4) * t179 + t438) * t125 - (-m(4) * t75 - t129 + t441) * t453 + m(6) * (t107 * t83 + t14 * t5 + t15 * t4 + t155 * t45 + t2 * t33 + t3 * t32) + (Ifges(6,5) * t65 + Ifges(6,6) * t66) * t418 + t32 * t19 + t33 * t20 + m(5) * (t110 * t17 + t111 * t18 + t49 * t86 + t50 * t87) + (Ifges(6,4) * t65 + Ifges(6,2) * t66) * t427 + t176 * t347 + (Ifges(3,4) * t242 + Ifges(3,2) * t241) * t411 + (Ifges(6,1) * t65 + Ifges(6,4) * t66) * t425 + t83 * t59 + t4 * t95 + t5 * t96 + t107 * (-mrSges(6,1) * t66 + mrSges(6,2) * t65) + t110 * t69 + t111 * t68 + t50 * t145 + t49 * t146 + t155 * t10 + t321 * t467 + t198 * t130 + t124 * t199 - pkin(1) * (-mrSges(3,1) * t241 + mrSges(3,2) * t242) - t274 * (mrSges(4,1) * t142 + mrSges(4,2) * t141) + Ifges(2,3) * qJDD(1) + t65 * t429 + t66 * t430 + t374 * t431 - t166 * t434 - t165 * t435 - t310 * t357 - t215 * t359 / 0.2e1 - t41 * t375 / 0.2e1 - t250 * t382; ((t410 * t75 + t291 * t74 + (-t179 * t291 + t180 * t410) * qJD(3)) * pkin(2) + t179 * t182 - t180 * t183 + t251 * t348) * m(4) - (-Ifges(3,2) * t361 + t216 + t276) * t360 / 0.2e1 + (-m(4) * t407 - m(6) * (t452 + t407) - m(5) * (t362 + t407) + t250 + t442) * g(3) + t445 * t262 + (t293 * t450 + t446) * g(2) + (t296 * t450 + t447) * g(1) + t454 * t182 + t460 * t95 + (-t107 * t137 + t461 * t14 + t460 * t15 + t172 * t3 + t173 * t2 + t245 * t45) * m(6) + t461 * t96 + (m(6) * t107 + t438 + t59) * pkin(2) * t358 + t130 * t409 + t68 * t373 + t129 * t351 + (t330 - t183) * t199 + (t448 + (t310 - t308 / 0.2e1) * qJD(1)) * qJD(1) + t298 + (-g(1) * t256 - g(2) * t255 - t157 * t182 + t262 * t316 + t269 * t317 + t273 * t71 - t86 * t97 - t87 * t98) * m(5) + Ifges(3,3) * qJDD(2) - t137 * t59 + (m(4) * t408 + t325 + t327) * t443 - t98 * t145 - t97 * t146 + t172 * t19 + t173 * t20 - t228 * mrSges(3,2) - t229 * mrSges(3,1) + Ifges(3,6) * t241 + Ifges(3,5) * t242 + t245 * t10 + t273 * t60 - t176 * t348 - t319 * t357 / 0.2e1 + t215 * t361 / 0.2e1 - t269 * t385; -pkin(3) * t60 - qJ(4) * t385 - t270 * t10 - t100 * t145 - t132 * t59 - t99 * t146 - t179 * t199 + t190 * t19 + t191 * t20 + t68 * t383 + t298 + t459 * t96 + t458 * t95 + t443 * t325 + t454 * t180 + t445 * qJD(4) + (t449 * t293 + t446) * g(2) + (t449 * t296 + t447) * g(1) + (-t107 * t132 + t14 * t459 + t15 * t458 + t190 * t3 + t191 * t2 - t270 * t45) * m(6) + (-g(1) * (-t296 * t406 + t256) - g(2) * (-t293 * t406 + t255) - pkin(3) * t71 + t317 * qJ(4) + t316 * qJD(4) - t100 * t87 - t157 * t180 - t86 * t99) * m(5) + (-m(5) * t362 + t442 - t480) * g(3); t119 * t96 - t469 * t95 - t333 * t145 + t194 * t146 - m(5) * (-t194 * t86 + t333 * t87) + t10 + t441 + (g(3) * t281 - t280 * t443) * (m(5) + m(6)) + (t119 * t14 - t15 * t469 + t45) * m(6); -t107 * (mrSges(6,1) * t119 + mrSges(6,2) * t469) + (Ifges(6,1) * t469 - t393) * t426 + t52 * t425 + (Ifges(6,5) * t469 - Ifges(6,6) * t119) * t419 - t14 * t95 + t15 * t96 - g(1) * (mrSges(6,1) * t203 - mrSges(6,2) * t204) - g(2) * (-mrSges(6,1) * t201 + mrSges(6,2) * t202) - g(3) * (-mrSges(6,1) * t278 - mrSges(6,2) * t279) * t280 + (t119 * t15 + t14 * t469) * mrSges(6,3) + t352 + (-Ifges(6,2) * t119 + t114 + t53) * t428 + t439;];
tau = t1;
