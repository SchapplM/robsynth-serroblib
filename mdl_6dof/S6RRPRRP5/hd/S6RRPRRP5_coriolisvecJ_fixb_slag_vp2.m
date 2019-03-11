% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:50
% EndTime: 2019-03-09 11:58:51
% DurationCPUTime: 30.95s
% Computational Cost: add. (15732->793), mult. (47356->1124), div. (0->0), fcn. (37597->10), ass. (0->330)
t270 = sin(pkin(6));
t269 = sin(pkin(11));
t274 = sin(qJ(2));
t277 = cos(qJ(2));
t372 = cos(pkin(11));
t323 = t372 * t277;
t286 = -t269 * t274 + t323;
t282 = t270 * t286;
t221 = qJD(1) * t282;
t445 = t221 - qJD(4);
t474 = Ifges(6,4) + Ifges(7,4);
t475 = Ifges(6,1) + Ifges(7,1);
t473 = Ifges(6,5) + Ifges(7,5);
t472 = Ifges(6,2) + Ifges(7,2);
t471 = Ifges(6,6) + Ifges(7,6);
t271 = cos(pkin(6));
t402 = pkin(1) * t271;
t264 = t277 * t402;
t259 = qJD(1) * t264;
t394 = pkin(8) + qJ(3);
t332 = t274 * t394;
t318 = t270 * t332;
t208 = -qJD(1) * t318 + t259;
t263 = t274 * t402;
t366 = t270 * t277;
t442 = t366 * t394 + t263;
t209 = t442 * qJD(1);
t325 = t372 * t209;
t160 = t208 * t269 + t325;
t273 = sin(qJ(4));
t276 = cos(qJ(4));
t486 = -t160 - t445 * (pkin(4) * t273 - pkin(10) * t276);
t324 = t372 * t274;
t355 = qJD(1) * t270;
t340 = t277 * t355;
t222 = -t269 * t340 - t324 * t355;
t261 = qJD(1) * t271 + qJD(2);
t189 = t222 * t273 + t261 * t276;
t224 = qJD(2) * t282;
t217 = qJD(1) * t224;
t142 = qJD(4) * t189 + t217 * t276;
t190 = -t222 * t276 + t261 * t273;
t219 = -t286 * t355 + qJD(4);
t272 = sin(qJ(5));
t275 = cos(qJ(5));
t144 = -t190 * t272 + t219 * t275;
t228 = (t269 * t277 + t324) * t270;
t223 = qJD(2) * t228;
t216 = qJD(1) * t223;
t73 = qJD(5) * t144 + t142 * t275 + t216 * t272;
t436 = t73 / 0.2e1;
t145 = t190 * t275 + t219 * t272;
t74 = -qJD(5) * t145 - t142 * t272 + t216 * t275;
t435 = t74 / 0.2e1;
t143 = qJD(4) * t190 + t217 * t273;
t426 = t143 / 0.2e1;
t470 = Ifges(6,3) + Ifges(7,3);
t485 = t474 * t144;
t364 = t272 * t276;
t178 = -t221 * t364 - t222 * t275;
t349 = qJD(4) * t276;
t337 = t272 * t349;
t346 = qJD(5) * t275;
t455 = t273 * t346 + t178 + t337;
t484 = t474 * t145;
t378 = t219 * Ifges(5,6);
t389 = Ifges(5,4) * t190;
t110 = t189 * Ifges(5,2) + t378 + t389;
t191 = pkin(2) * t261 + t208;
t133 = t269 * t191 + t325;
t127 = pkin(9) * t261 + t133;
t248 = (-pkin(2) * t277 - pkin(1)) * t270;
t241 = qJD(1) * t248 + qJD(3);
t154 = -t221 * pkin(3) + t222 * pkin(9) + t241;
t83 = t127 * t276 + t154 * t273;
t396 = t83 * mrSges(5,3);
t483 = t396 + t110 / 0.2e1;
t187 = qJD(5) - t189;
t60 = t145 * Ifges(7,5) + t144 * Ifges(7,6) + t187 * Ifges(7,3);
t61 = t145 * Ifges(6,5) + t144 * Ifges(6,6) + t187 * Ifges(6,3);
t466 = t61 + t60;
t482 = t466 / 0.2e1;
t481 = t473 * t426 + t474 * t435 + t436 * t475;
t465 = t472 * t144 + t471 * t187 + t484;
t464 = t475 * t145 + t473 * t187 + t485;
t401 = pkin(2) * t269;
t266 = pkin(9) + t401;
t338 = t266 * t349;
t199 = t269 * t209;
t161 = t208 * t372 - t199;
t341 = t274 * t355;
t322 = pkin(2) * t341;
t173 = -pkin(3) * t222 - pkin(9) * t221 + t322;
t97 = -t273 * t161 + t173 * t276;
t84 = pkin(4) * t222 - t97;
t480 = t338 - t84;
t351 = qJD(4) * t273;
t339 = t266 * t351;
t98 = t276 * t161 + t273 * t173;
t85 = -pkin(10) * t222 + t98;
t479 = t486 * t275 + (t339 + t85) * t272;
t342 = t372 * pkin(2);
t267 = -t342 - pkin(3);
t246 = -t276 * pkin(4) - t273 * pkin(10) + t267;
t478 = t246 * t346 + t486 * t272 - t275 * t85;
t477 = t474 * t275;
t476 = t474 * t272;
t469 = t143 * t470 + t471 * t74 + t473 * t73;
t468 = t143 * t471 + t472 * t74 + t474 * t73;
t112 = mrSges(5,1) * t216 - mrSges(5,3) * t142;
t23 = -mrSges(6,1) * t74 + mrSges(6,2) * t73;
t463 = -t112 + t23;
t362 = t275 * t276;
t179 = t221 * t362 - t222 * t272;
t249 = t266 * t362;
t345 = qJD(6) * t275;
t370 = t221 * t273;
t462 = -pkin(5) * t370 + qJ(6) * t179 - t273 * t345 + (pkin(5) * t273 - qJ(6) * t362) * qJD(4) + (-t249 + (qJ(6) * t273 - t246) * t272) * qJD(5) + t479;
t363 = t273 * t275;
t461 = -qJ(6) * t178 + (-qJ(6) * qJD(5) - qJD(4) * t266) * t363 + (-qJD(6) * t273 + (-qJ(6) * qJD(4) - qJD(5) * t266) * t276) * t272 + t478;
t393 = -qJ(6) - pkin(10);
t326 = qJD(5) * t393;
t131 = pkin(4) * t190 - pkin(10) * t189;
t82 = -t273 * t127 + t154 * t276;
t49 = t275 * t131 - t272 * t82;
t460 = -pkin(5) * t190 - qJD(6) * t272 - t49 + (qJ(6) * t189 + t326) * t275;
t371 = t189 * t272;
t50 = t272 * t131 + t275 * t82;
t459 = qJ(6) * t371 + t272 * t326 + t345 - t50;
t198 = t272 * t246 + t249;
t458 = -qJD(5) * t198 + t479;
t348 = qJD(5) * t272;
t350 = qJD(4) * t275;
t457 = (-t273 * t350 - t276 * t348) * t266 + t478;
t456 = pkin(5) * t455 + t480;
t347 = qJD(5) * t273;
t454 = t272 * t347 - t275 * t349 + t179;
t391 = mrSges(4,3) * t222;
t359 = -mrSges(4,1) * t261 - mrSges(5,1) * t189 + mrSges(5,2) * t190 - t391;
t453 = -t272 * t471 + t275 * t473;
t452 = t272 * t473 + t275 * t471;
t451 = -t272 * t472 + t477;
t450 = t275 * t472 + t476;
t449 = t275 * t475 - t476;
t448 = t272 * t475 + t477;
t447 = t351 - t370;
t254 = qJD(2) * t259;
t281 = (-qJD(2) * t332 + qJD(3) * t277) * t270;
t183 = qJD(1) * t281 + t254;
t367 = t270 * t274;
t193 = -qJD(2) * t442 - qJD(3) * t367;
t279 = qJD(1) * t193;
t118 = t183 * t372 + t269 * t279;
t353 = qJD(2) * t270;
t331 = qJD(1) * t353;
t317 = t274 * t331;
t289 = pkin(2) * t317;
t155 = pkin(3) * t216 - pkin(9) * t217 + t289;
t32 = t276 * t118 - t127 * t351 + t154 * t349 + t273 * t155;
t33 = -t273 * t118 - t127 * t349 - t154 * t351 + t155 * t276;
t446 = -t273 * t33 + t276 * t32;
t27 = pkin(10) * t216 + t32;
t117 = t183 * t269 - t372 * t279;
t53 = pkin(4) * t143 - pkin(10) * t142 + t117;
t67 = pkin(10) * t219 + t83;
t132 = t191 * t372 - t199;
t126 = -t261 * pkin(3) - t132;
t79 = -t189 * pkin(4) - t190 * pkin(10) + t126;
t5 = t275 * t27 + t272 * t53 + t79 * t346 - t348 * t67;
t26 = t272 * t79 + t275 * t67;
t6 = -qJD(5) * t26 - t27 * t272 + t275 * t53;
t309 = -t272 * t6 + t275 * t5;
t25 = -t272 * t67 + t275 * t79;
t18 = -qJ(6) * t145 + t25;
t10 = pkin(5) * t187 + t18;
t19 = qJ(6) * t144 + t26;
t290 = t25 * t275 + t26 * t272;
t307 = mrSges(7,1) * t272 + mrSges(7,2) * t275;
t308 = mrSges(6,1) * t272 + mrSges(6,2) * t275;
t404 = t275 / 0.2e1;
t415 = t187 / 0.2e1;
t422 = t145 / 0.2e1;
t424 = t144 / 0.2e1;
t66 = -pkin(4) * t219 - t82;
t44 = -pkin(5) * t144 + qJD(6) + t66;
t441 = -t290 * mrSges(6,3) - (t10 * t275 + t19 * t272) * mrSges(7,3) + t307 * t44 + t308 * t66 + t451 * t424 + t449 * t422 + t453 * t415 - t465 * t272 / 0.2e1 + t464 * t404;
t440 = -0.2e1 * pkin(1);
t439 = m(5) / 0.2e1;
t438 = m(6) / 0.2e1;
t437 = m(7) / 0.2e1;
t432 = m(7) * t44;
t431 = -t110 / 0.2e1;
t185 = Ifges(5,4) * t189;
t379 = t219 * Ifges(5,5);
t111 = t190 * Ifges(5,1) + t185 + t379;
t429 = -t111 / 0.2e1;
t428 = t142 / 0.2e1;
t427 = -t143 / 0.2e1;
t425 = -t144 / 0.2e1;
t423 = -t145 / 0.2e1;
t416 = -t187 / 0.2e1;
t202 = t228 * t273 - t271 * t276;
t414 = -t202 / 0.2e1;
t203 = t228 * t276 + t271 * t273;
t412 = t203 / 0.2e1;
t411 = -t221 / 0.2e1;
t410 = -t222 / 0.2e1;
t408 = t271 / 0.2e1;
t405 = -t274 / 0.2e1;
t400 = pkin(5) * t272;
t397 = t82 * mrSges(5,3);
t205 = pkin(2) * t271 + t264 - t318;
t356 = pkin(8) * t366 + t263;
t220 = qJ(3) * t366 + t356;
t169 = t269 * t205 + t372 * t220;
t157 = pkin(9) * t271 + t169;
t227 = t269 * t367 - t270 * t323;
t177 = t227 * pkin(3) - t228 * pkin(9) + t248;
t101 = t276 * t157 + t273 * t177;
t90 = pkin(10) * t227 + t101;
t168 = t205 * t372 - t269 * t220;
t156 = -t271 * pkin(3) - t168;
t99 = t202 * pkin(4) - t203 * pkin(10) + t156;
t38 = t272 * t99 + t275 * t90;
t392 = mrSges(4,3) * t221;
t390 = Ifges(3,4) * t274;
t388 = Ifges(5,4) * t273;
t387 = Ifges(5,4) * t276;
t382 = Ifges(3,5) * t277;
t377 = t222 * Ifges(4,4);
t28 = -pkin(4) * t216 - t33;
t376 = t273 * t28;
t153 = mrSges(5,1) * t219 - mrSges(5,3) * t190;
t93 = -mrSges(6,1) * t144 + mrSges(6,2) * t145;
t373 = -t93 + t153;
t369 = t221 * t276;
t365 = t272 * t273;
t105 = -mrSges(7,2) * t187 + mrSges(7,3) * t144;
t106 = -mrSges(6,2) * t187 + mrSges(6,3) * t144;
t361 = t105 + t106;
t107 = mrSges(7,1) * t187 - mrSges(7,3) * t145;
t108 = mrSges(6,1) * t187 - mrSges(6,3) * t145;
t360 = -t107 - t108;
t352 = qJD(4) * t272;
t343 = Ifges(5,5) * t142 - Ifges(5,6) * t143 + Ifges(5,3) * t216;
t336 = t274 * t353;
t22 = -t74 * mrSges(7,1) + t73 * mrSges(7,2);
t328 = t349 / 0.2e1;
t327 = -t347 / 0.2e1;
t37 = -t272 * t90 + t275 * t99;
t100 = -t273 * t157 + t177 * t276;
t260 = qJD(2) * t264;
t192 = t260 + t281;
t129 = t192 * t269 - t372 * t193;
t321 = pkin(2) * t336;
t320 = mrSges(3,3) * t341;
t319 = mrSges(3,3) * t340;
t1 = pkin(5) * t143 - qJ(6) * t73 - qJD(6) * t145 + t6;
t2 = qJ(6) * t74 + qJD(6) * t144 + t5;
t310 = -t1 * t272 + t2 * t275;
t306 = Ifges(5,1) * t276 - t388;
t301 = -Ifges(5,2) * t273 + t387;
t296 = Ifges(5,5) * t276 - Ifges(5,6) * t273;
t176 = t203 * t275 + t227 * t272;
t175 = -t203 * t272 + t227 * t275;
t130 = t192 * t372 + t269 * t193;
t174 = pkin(3) * t223 - pkin(9) * t224 + t321;
t43 = -t273 * t130 - t157 * t349 + t174 * t276 - t177 * t351;
t89 = -pkin(4) * t227 - t100;
t42 = t276 * t130 - t157 * t351 + t273 * t174 + t177 * t349;
t35 = pkin(10) * t223 + t42;
t166 = qJD(4) * t203 + t224 * t273;
t167 = -qJD(4) * t202 + t224 * t276;
t58 = pkin(4) * t166 - pkin(10) * t167 + t129;
t7 = t272 * t58 + t275 * t35 + t99 * t346 - t348 * t90;
t287 = t261 * (-Ifges(3,6) * t274 + t382);
t36 = -pkin(4) * t223 - t43;
t235 = t356 * qJD(2);
t225 = -pkin(8) * t317 + t254;
t226 = qJD(1) * t235;
t280 = -t226 * mrSges(3,1) - t117 * mrSges(4,1) - t225 * mrSges(3,2) - t118 * mrSges(4,2);
t8 = -qJD(5) * t38 - t272 * t35 + t275 * t58;
t268 = -pkin(5) * t275 - pkin(4);
t258 = Ifges(3,4) * t340;
t256 = t393 * t275;
t255 = t393 * t272;
t253 = t331 * t382;
t244 = -pkin(8) * t367 + t264;
t240 = (t266 + t400) * t273;
t237 = t275 * t246;
t234 = -pkin(8) * t336 + t260;
t233 = t356 * qJD(1);
t232 = -pkin(8) * t341 + t259;
t231 = -t261 * mrSges(3,2) + t319;
t230 = mrSges(3,1) * t261 - t320;
t218 = Ifges(4,4) * t221;
t213 = Ifges(4,5) * t217;
t212 = Ifges(4,6) * t216;
t210 = t217 * mrSges(4,2);
t207 = Ifges(3,1) * t341 + t261 * Ifges(3,5) + t258;
t206 = Ifges(3,6) * t261 + (Ifges(3,2) * t277 + t390) * t355;
t197 = -t266 * t364 + t237;
t194 = -mrSges(4,2) * t261 + t392;
t188 = -qJ(6) * t365 + t198;
t182 = -qJ(6) * t363 + t237 + (-t266 * t272 - pkin(5)) * t276;
t180 = -mrSges(4,1) * t221 - mrSges(4,2) * t222;
t172 = -t222 * Ifges(4,1) + t261 * Ifges(4,5) + t218;
t171 = t221 * Ifges(4,2) + t261 * Ifges(4,6) - t377;
t152 = -mrSges(5,2) * t219 + mrSges(5,3) * t189;
t113 = -mrSges(5,2) * t216 - mrSges(5,3) * t143;
t109 = t190 * Ifges(5,5) + t189 * Ifges(5,6) + t219 * Ifges(5,3);
t92 = -mrSges(7,1) * t144 + mrSges(7,2) * t145;
t91 = mrSges(5,1) * t143 + mrSges(5,2) * t142;
t88 = qJD(5) * t175 + t167 * t275 + t223 * t272;
t87 = -qJD(5) * t176 - t167 * t272 + t223 * t275;
t76 = t142 * Ifges(5,1) - t143 * Ifges(5,4) + t216 * Ifges(5,5);
t75 = t142 * Ifges(5,4) - t143 * Ifges(5,2) + t216 * Ifges(5,6);
t59 = pkin(5) * t371 + t83;
t55 = -pkin(5) * t175 + t89;
t48 = -mrSges(6,2) * t143 + mrSges(6,3) * t74;
t47 = -mrSges(7,2) * t143 + mrSges(7,3) * t74;
t46 = mrSges(6,1) * t143 - mrSges(6,3) * t73;
t45 = mrSges(7,1) * t143 - mrSges(7,3) * t73;
t29 = qJ(6) * t175 + t38;
t20 = pkin(5) * t202 - qJ(6) * t176 + t37;
t11 = -pkin(5) * t87 + t36;
t9 = -pkin(5) * t74 + t28;
t4 = qJ(6) * t87 + qJD(6) * t175 + t7;
t3 = pkin(5) * t166 - qJ(6) * t88 - qJD(6) * t176 + t8;
t12 = [(t225 * t277 + t226 * t274) * t270 * mrSges(3,3) + t464 * t88 / 0.2e1 + t465 * t87 / 0.2e1 + t468 * t175 / 0.2e1 + t469 * t202 / 0.2e1 + (t175 * t472 + t176 * t474 + t202 * t471) * t435 + (t474 * t175 + t176 * t475 + t473 * t202) * t436 + (t473 * t166 + t474 * t87 + t475 * t88) * t422 + (t175 * t471 + t176 * t473 + t202 * t470) * t426 + (t166 * t470 + t471 * t87 + t473 * t88) * t415 + (t166 * t471 + t472 * t87 + t474 * t88) * t424 + m(3) * (t225 * t356 - t226 * t244 - t232 * t235 + t233 * t234) + ((Ifges(3,5) * t408 - t244 * mrSges(3,3) + (mrSges(3,2) * t440 + 0.3e1 / 0.2e1 * Ifges(3,4) * t277) * t270) * t277 + (-t356 * mrSges(3,3) - Ifges(3,6) * t271 + (mrSges(3,1) * t440 - 0.3e1 / 0.2e1 * t390 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t277) * t270 + (m(4) * t248 + mrSges(4,1) * t227 + mrSges(4,2) * t228) * pkin(2)) * t274) * t331 + t248 * t210 + (t117 * t228 - t118 * t227 - t132 * t224 - t133 * t223) * mrSges(4,3) + t359 * t129 + m(7) * (t1 * t20 + t10 * t3 + t11 * t44 + t19 * t4 + t2 * t29 + t55 * t9) + m(6) * (t25 * t8 + t26 * t7 + t28 * t89 + t36 * t66 + t37 * t6 + t38 * t5) + m(5) * (t100 * t33 + t101 * t32 + t117 * t156 + t126 * t129 + t42 * t83 + t43 * t82) + (t253 / 0.2e1 + t213 / 0.2e1 - t212 / 0.2e1 + t280) * t271 + t89 * t23 + t11 * t92 + t36 * t93 + t66 * (-mrSges(6,1) * t87 + mrSges(6,2) * t88) + t44 * (-mrSges(7,1) * t87 + mrSges(7,2) * t88) + t227 * t343 / 0.2e1 + t55 * t22 + m(4) * (-t117 * t168 + t118 * t169 - t129 * t132 + t130 * t133 + t241 * t321) + t38 * t48 + t20 * t45 + t37 * t46 + t29 * t47 + (t274 * pkin(2) * t180 + t277 * t207 / 0.2e1 + t206 * t405 + t287 / 0.2e1 + (-t232 * t277 - t233 * t274) * mrSges(3,3)) * t353 + (-t168 * mrSges(4,3) + Ifges(4,1) * t228 - Ifges(4,4) * t227 + Ifges(4,5) * t408) * t217 + (-t169 * mrSges(4,3) + Ifges(5,5) * t412 + Ifges(5,6) * t414 - Ifges(4,4) * t228 + t248 * mrSges(4,1) - Ifges(4,6) * t271 / 0.2e1 + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t227) * t216 + t4 * t105 + t7 * t106 + t3 * t107 + t8 * t108 + t100 * t112 + t101 * t113 + t42 * t152 + t43 * t153 + t156 * t91 + t25 * (mrSges(6,1) * t166 - mrSges(6,3) * t88) + t10 * (mrSges(7,1) * t166 - mrSges(7,3) * t88) + t26 * (-mrSges(6,2) * t166 + mrSges(6,3) * t87) + t19 * (-mrSges(7,2) * t166 + mrSges(7,3) * t87) + t167 * t111 / 0.2e1 + t126 * (mrSges(5,1) * t166 + mrSges(5,2) * t167) + t28 * (-mrSges(6,1) * t175 + mrSges(6,2) * t176) + t9 * (-mrSges(7,1) * t175 + mrSges(7,2) * t176) + t130 * t194 + t5 * (-mrSges(6,2) * t202 + mrSges(6,3) * t175) + t2 * (-mrSges(7,2) * t202 + mrSges(7,3) * t175) + t6 * (mrSges(6,1) * t202 - mrSges(6,3) * t176) + t1 * (mrSges(7,1) * t202 - mrSges(7,3) * t176) + t117 * (mrSges(5,1) * t202 + mrSges(5,2) * t203) + t223 * t109 / 0.2e1 + t83 * (-mrSges(5,2) * t223 - mrSges(5,3) * t166) + t82 * (mrSges(5,1) * t223 - mrSges(5,3) * t167) + t189 * (Ifges(5,4) * t167 - Ifges(5,2) * t166 + Ifges(5,6) * t223) / 0.2e1 + t190 * (Ifges(5,1) * t167 - Ifges(5,4) * t166 + Ifges(5,5) * t223) / 0.2e1 + t219 * (Ifges(5,5) * t167 - Ifges(5,6) * t166 + Ifges(5,3) * t223) / 0.2e1 - t223 * t171 / 0.2e1 + t224 * t172 / 0.2e1 + t176 * t481 + t166 * t482 + (Ifges(4,1) * t224 - Ifges(4,4) * t223) * t410 + t76 * t412 + t75 * t414 + (Ifges(5,4) * t203 - Ifges(5,2) * t202 + Ifges(5,6) * t227) * t427 + (Ifges(5,1) * t203 - Ifges(5,4) * t202 + Ifges(5,5) * t227) * t428 + t166 * t431 + t221 * (Ifges(4,4) * t224 - Ifges(4,2) * t223) / 0.2e1 + t32 * (-mrSges(5,2) * t227 - mrSges(5,3) * t202) + t33 * (mrSges(5,1) * t227 - mrSges(5,3) * t203) + t234 * t231 - t235 * t230 + t241 * (mrSges(4,1) * t223 + mrSges(4,2) * t224) + t261 * (Ifges(4,5) * t224 - Ifges(4,6) * t223) / 0.2e1; (pkin(1) * (mrSges(3,1) * t274 + mrSges(3,2) * t277) + (Ifges(3,1) * t277 - t390) * t405) * qJD(1) ^ 2 * t270 ^ 2 - t359 * t160 + (mrSges(7,1) * t447 + mrSges(7,3) * t454) * t10 + (mrSges(6,1) * t447 + mrSges(6,3) * t454) * t25 + (mrSges(6,1) * t455 - mrSges(6,2) * t454) * t66 + (mrSges(7,1) * t455 - mrSges(7,2) * t454) * t44 + (-mrSges(7,2) * t447 - mrSges(7,3) * t455) * t19 + (-mrSges(6,2) * t447 - mrSges(6,3) * t455) * t26 + t480 * t93 - t468 * t365 / 0.2e1 - t469 * t276 / 0.2e1 + (t369 * t82 + t446) * mrSges(5,3) + (t474 * t178 + t179 * t475) * t423 + (t178 * t471 + t179 * t473) * t416 + (t178 * t472 + t179 * t474) * t425 + (t117 * t267 - t126 * t160 - t82 * t97 - t83 * t98) * m(5) + (t197 * t6 + t198 * t5 + t458 * t25 + t457 * t26 - t66 * t84) * m(6) + (((-t273 * t83 - t276 * t82) * qJD(4) + t446) * m(5) + t463 * t273 + (t349 * t66 + t376) * m(6) + t276 * t113) * t266 + t456 * t92 + t457 * t106 + t458 * t108 + t461 * t105 + (t1 * t182 + t10 * t462 + t188 * t2 + t19 * t461 + t240 * t9 + t44 * t456) * m(7) + t462 * t107 + (-t452 * t347 + (t273 * t470 + t276 * t453) * qJD(4)) * t415 + (t273 * t453 - t276 * t470) * t426 + (-t450 * t347 + (t273 * t471 + t276 * t451) * qJD(4)) * t424 + (t273 * t451 - t276 * t471) * t435 + (-t448 * t347 + (t273 * t473 + t276 * t449) * qJD(4)) * t422 + (t273 * t449 - t276 * t473) * t436 + (-t216 * t401 - t217 * t342) * mrSges(4,3) + (t132 * t160 - t133 * t161 - t241 * t322 + (-t117 * t372 + t118 * t269) * pkin(2)) * m(4) - ((-Ifges(3,2) * t341 + t207 + t258) * t277 + t287) * t355 / 0.2e1 + (t320 + t230) * t233 + t280 + (t189 * t301 + t190 * t306 + t219 * t296) * qJD(4) / 0.2e1 + (Ifges(4,1) * t221 + t109 + t377) * t222 / 0.2e1 - t83 * mrSges(5,2) * t222 + t82 * mrSges(5,1) * t222 + (Ifges(4,2) * t222 + t172 + t218) * t411 - t133 * t391 + t1 * (-mrSges(7,1) * t276 - mrSges(7,3) * t363) + t6 * (-mrSges(6,1) * t276 - mrSges(6,3) * t363) + t2 * (mrSges(7,2) * t276 - mrSges(7,3) * t365) - t445 * t126 * (mrSges(5,1) * t273 + mrSges(5,2) * t276) + (t319 - t231) * t232 + t253 + t111 * t328 + (-t339 - t98) * t152 - t212 + t213 + t9 * t307 * t273 - t349 * t397 + (-t338 - t97) * t153 + t5 * (mrSges(6,2) * t276 - mrSges(6,3) * t365) + t206 * t341 / 0.2e1 - t180 * t322 - Ifges(3,6) * t317 + (-t466 / 0.2e1 + t473 * t423 + t470 * t416 + t471 * t425 + t483) * t370 + t464 * (t272 * t327 + t275 * t328 - t179 / 0.2e1) + t465 * (t275 * t327 - t337 / 0.2e1 - t178 / 0.2e1) + (t482 - t396 + t431) * t351 + t182 * t45 + t188 * t47 - t161 * t194 + t197 * t46 + t198 * t48 + t308 * t376 + t132 * t392 + t363 * t481 + t171 * t410 + (Ifges(5,2) * t276 + t388) * t427 + (Ifges(5,1) * t273 + t387) * t428 + t369 * t429 + t240 * t22 - t241 * (-mrSges(4,1) * t222 + mrSges(4,2) * t221) - t261 * (Ifges(4,5) * t221 + Ifges(4,6) * t222) / 0.2e1 + t267 * t91 - t219 * (-Ifges(5,3) * t222 + t221 * t296) / 0.2e1 - t189 * (-Ifges(5,6) * t222 + t221 * t301) / 0.2e1 - t190 * (-Ifges(5,5) * t222 + t221 * t306) / 0.2e1 + t273 * t76 / 0.2e1 + t276 * t75 / 0.2e1 + t117 * (-mrSges(5,1) * t276 + mrSges(5,2) * t273) + t216 * (Ifges(5,5) * t273 + Ifges(5,6) * t276) / 0.2e1; t216 * mrSges(4,1) - t221 * t194 + t210 - t361 * t179 + t360 * t178 + (-t221 * t152 - t22 + (t272 * t360 + t275 * t361 + t152) * qJD(4) - t463) * t276 + (t113 + (t47 + t48) * t275 + (-t45 - t46) * t272 + (-t272 * t361 + t275 * t360) * qJD(5) + t445 * (-t92 + t373)) * t273 - m(7) * (t10 * t178 + t179 * t19) - m(6) * (t178 * t25 + t179 * t26) + 0.2e1 * ((-t10 * t352 + t19 * t350 - t9) * t437 + (-t25 * t352 + t26 * t350 - t28) * t438 + m(5) * t83 * t411 + (qJD(4) * t83 + t33) * t439) * t276 + 0.2e1 * ((qJD(4) * t44 - t10 * t346 - t19 * t348 + t310) * t437 + (qJD(4) * t66 - t25 * t346 - t26 * t348 + t309) * t438 + (-qJD(4) * t82 + t32) * t439 + (t82 * t439 - t432 / 0.2e1 - m(6) * t66 / 0.2e1) * t221) * t273 + (m(5) * t126 + t359) * t222 + (-t132 * t222 - t133 * t221 + t289) * m(4); (-t185 / 0.2e1 + t429 - t126 * mrSges(5,2) + t397 - t379 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t190 - t441) * t189 + ((t92 + t432) * t400 + t441) * qJD(5) + t468 * t404 + t448 * t436 + t450 * t435 + t452 * t426 + ((-m(6) * t290 - t106 * t272 - t275 * t108) * qJD(5) - t272 * t46 + t275 * t48 + m(6) * t309) * pkin(10) + t459 * t105 + (t1 * t255 + t10 * t460 + t19 * t459 - t2 * t256 + t268 * t9 - t44 * t59) * m(7) + t460 * t107 + t343 + t373 * t83 - t59 * t92 + t33 * mrSges(5,1) - t32 * mrSges(5,2) - pkin(4) * t23 + (-pkin(4) * t28 - t25 * t49 - t26 * t50 - t66 * t83) * m(6) - t50 * t106 - t49 * t108 + (t378 / 0.2e1 + t19 * mrSges(7,2) - t10 * mrSges(7,1) + t26 * mrSges(6,2) - t25 * mrSges(6,1) - t126 * mrSges(5,1) - t60 / 0.2e1 - t61 / 0.2e1 + t389 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t187 + (-Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1) * t145 + (-Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t144 + t483) * t190 - t82 * t152 + t272 * t481 + t255 * t45 - t256 * t47 + t268 * t22 + t28 * (-mrSges(6,1) * t275 + mrSges(6,2) * t272) + t9 * (-mrSges(7,1) * t275 + mrSges(7,2) * t272) + t309 * mrSges(6,3) + t310 * mrSges(7,3); (-t145 * t92 + t45) * pkin(5) + t469 + (-(-t10 + t18) * t19 + (-t145 * t44 + t1) * pkin(5)) * m(7) + t6 * mrSges(6,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t1 * mrSges(7,1) + (t10 * t144 + t145 * t19) * mrSges(7,3) + (t144 * t25 + t145 * t26) * mrSges(6,3) - t18 * t105 - t25 * t106 + t19 * t107 + t26 * t108 - t66 * (mrSges(6,1) * t145 + mrSges(6,2) * t144) - t44 * (mrSges(7,1) * t145 + mrSges(7,2) * t144) + (t144 * t475 - t484) * t423 + t465 * t422 + (t144 * t473 - t145 * t471) * t416 + (-t145 * t472 + t464 + t485) * t425; -t144 * t105 + t145 * t107 + 0.2e1 * (t9 / 0.2e1 + t10 * t422 + t19 * t425) * m(7) + t22;];
tauc  = t12(:);
