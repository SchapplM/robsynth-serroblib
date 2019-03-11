% Calculate vector of inverse dynamics joint torques for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:21
% EndTime: 2019-03-08 18:52:48
% DurationCPUTime: 14.66s
% Computational Cost: add. (6095->558), mult. (15526->773), div. (0->0), fcn. (13603->14), ass. (0->269)
t449 = Ifges(6,4) + Ifges(7,4);
t450 = Ifges(6,1) + Ifges(7,1);
t432 = Ifges(6,5) + Ifges(7,5);
t448 = Ifges(6,2) + Ifges(7,2);
t447 = Ifges(6,6) + Ifges(7,6);
t199 = sin(qJ(4));
t202 = cos(qJ(4));
t257 = pkin(4) * t199 - pkin(10) * t202;
t342 = cos(pkin(6));
t186 = qJD(1) * t342 + qJD(2);
t193 = sin(pkin(12));
t195 = sin(pkin(6));
t200 = sin(qJ(3));
t196 = cos(pkin(12));
t341 = cos(pkin(7));
t273 = t196 * t341;
t260 = t200 * t273;
t368 = cos(qJ(3));
t211 = (t193 * t368 + t260) * t195;
t194 = sin(pkin(7));
t332 = t194 * t200;
t90 = qJD(1) * t211 + t186 * t332;
t457 = t257 * qJD(4) - t90;
t201 = cos(qJ(5));
t456 = t449 * t201;
t198 = sin(qJ(5));
t455 = t449 * t198;
t319 = qJD(4) * t201;
t323 = qJD(3) * t199;
t163 = -t198 * t323 + t319;
t313 = qJD(3) * qJD(4);
t173 = qJDD(3) * t199 + t202 * t313;
t104 = qJD(5) * t163 + qJDD(4) * t198 + t173 * t201;
t380 = t104 / 0.2e1;
t164 = qJD(4) * t198 + t201 * t323;
t105 = -qJD(5) * t164 + qJDD(4) * t201 - t173 * t198;
t379 = t105 / 0.2e1;
t172 = qJDD(3) * t202 - t199 * t313;
t160 = qJDD(5) - t172;
t378 = t160 / 0.2e1;
t454 = -t163 / 0.2e1;
t453 = -t164 / 0.2e1;
t321 = qJD(3) * t202;
t438 = t321 - qJD(5);
t372 = -t438 / 0.2e1;
t190 = pkin(5) * t201 + pkin(4);
t452 = m(7) * t190;
t451 = qJD(4) / 0.2e1;
t446 = -Ifges(7,3) - Ifges(6,3);
t408 = -t198 * t447 + t201 * t432;
t406 = -t198 * t448 + t456;
t404 = t201 * t450 - t455;
t197 = -qJ(6) - pkin(10);
t445 = -m(7) * t197 + mrSges(6,3) + mrSges(7,3);
t444 = t432 * t378 + t449 * t379 + t380 * t450;
t233 = t368 * t273;
t227 = t195 * t233;
t221 = qJD(1) * t227;
t334 = t193 * t195;
t299 = qJD(1) * t334;
t300 = t194 * t368;
t206 = -t186 * t300 + t200 * t299 - t221;
t320 = qJD(4) * t199;
t311 = pkin(9) * t320;
t330 = t198 * t202;
t443 = t198 * t311 + t457 * t201 - t206 * t330;
t258 = -pkin(4) * t202 - pkin(10) * t199;
t178 = -pkin(3) + t258;
t315 = qJD(5) * t201;
t328 = t201 * t202;
t442 = t178 * t315 + t457 * t198 + t206 * t328;
t441 = t449 * t163;
t440 = t449 * t164;
t191 = Ifges(5,4) * t321;
t422 = t164 * t450 - t432 * t438 + t441;
t439 = Ifges(5,1) * t323 + Ifges(5,5) * qJD(4) + t201 * t422 + t191;
t357 = Ifges(5,4) * t199;
t245 = Ifges(5,2) * t202 + t357;
t437 = t446 * t372 + t432 * t453 + t447 * t454 + Ifges(5,6) * t451 + qJD(3) * t245 / 0.2e1;
t381 = m(7) * pkin(5);
t436 = t172 / 0.2e1;
t435 = t173 / 0.2e1;
t434 = -mrSges(6,1) - mrSges(7,1);
t433 = mrSges(6,2) + mrSges(7,2);
t430 = t104 * t449 + t105 * t448 + t160 * t447;
t188 = pkin(9) * t328;
t232 = pkin(5) * t199 - qJ(6) * t328;
t314 = qJD(6) * t201;
t428 = -t199 * t314 + t232 * qJD(4) + (-t188 + (qJ(6) * t199 - t178) * t198) * qJD(5) + t443;
t329 = t199 * t201;
t427 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t329 + (-qJD(6) * t199 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t202) * t198 + t442;
t317 = qJD(5) * t198;
t426 = (-t199 * t319 - t202 * t317) * pkin(9) + t442;
t133 = t198 * t178 + t188;
t425 = -qJD(5) * t133 + t443;
t423 = t163 * t448 - t438 * t447 + t440;
t277 = qJD(5) * t197;
t297 = t198 * t321;
t170 = t257 * qJD(3);
t304 = t195 * t196 * t194;
t131 = -qJD(1) * t304 + t186 * t341;
t88 = qJD(3) * pkin(9) + t90;
t66 = t131 * t202 - t199 * t88;
t40 = t198 * t170 + t201 * t66;
t421 = qJ(6) * t297 + t198 * t277 + t314 - t40;
t39 = t201 * t170 - t198 * t66;
t420 = -qJD(3) * t232 - qJD(6) * t198 + t201 * t277 - t39;
t335 = t131 * t199;
t366 = pkin(5) * t198;
t419 = pkin(5) * t317 - t335 - (qJD(3) * t366 + t88) * t202;
t418 = t381 + mrSges(7,1);
t57 = -mrSges(6,1) * t105 + mrSges(6,2) * t104;
t417 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t173 + t57;
t306 = mrSges(5,3) * t323;
t416 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t163 + mrSges(6,2) * t164 + t306;
t340 = cos(pkin(11));
t256 = t342 * t340;
t339 = sin(pkin(11));
t209 = t193 * t339 - t196 * t256;
t275 = t195 * t340;
t415 = t194 * t275 + t209 * t341;
t255 = t342 * t339;
t210 = t193 * t340 + t196 * t255;
t274 = t195 * t339;
t414 = -t194 * t274 + t210 * t341;
t413 = -t199 * t446 + t202 * t408;
t412 = t199 * t447 + t202 * t406;
t411 = t199 * t432 + t202 * t404;
t409 = t198 * t432 + t201 * t447;
t407 = t201 * t448 + t455;
t405 = t198 * t450 + t456;
t403 = t104 * t432 + t105 * t447 - t160 * t446;
t119 = -mrSges(5,1) * t172 + mrSges(5,2) * t173;
t402 = mrSges(4,1) * qJDD(3) - t119;
t254 = -mrSges(5,1) * t202 + mrSges(5,2) * t199;
t165 = t254 * qJD(3);
t401 = mrSges(4,1) * qJD(3) - t165;
t185 = qJDD(1) * t342 + qJDD(2);
t129 = -qJDD(1) * t304 + t185 * t341;
t318 = qJD(4) * t202;
t235 = t195 * t260;
t268 = qJD(3) * t299;
t269 = qJD(3) * t300;
t292 = qJDD(1) * t334;
t52 = qJD(3) * t221 + qJDD(1) * t235 + t185 * t332 + t186 * t269 - t200 * t268 + t368 * t292;
t50 = qJDD(3) * pkin(9) + t52;
t11 = t199 * t129 + t131 * t318 + t202 * t50 - t320 * t88;
t12 = t129 * t202 - t131 * t320 - t199 * t50 - t88 * t318;
t400 = t11 * t202 - t12 * t199;
t322 = qJD(3) * t200;
t298 = t194 * t322;
t53 = -qJD(3) * qJD(1) * t235 + qJDD(1) * t227 + t185 * t300 - t186 * t298 - t200 * t292 - t368 * t268;
t51 = -qJDD(3) * pkin(3) - t53;
t28 = -t172 * pkin(4) - t173 * pkin(10) + t51;
t67 = t202 * t88 + t335;
t62 = qJD(4) * pkin(10) + t67;
t82 = qJD(3) * t178 + t206;
t9 = qJDD(4) * pkin(10) + t11;
t3 = t198 * t28 + t201 * t9 + t82 * t315 - t317 * t62;
t21 = t198 * t82 + t201 * t62;
t4 = -qJD(5) * t21 - t198 * t9 + t201 * t28;
t399 = -t198 * t4 + t201 * t3;
t398 = -m(5) - m(6) - m(7);
t397 = -mrSges(6,1) - t418;
t251 = -mrSges(7,1) * t201 + mrSges(7,2) * t198;
t253 = -mrSges(6,1) * t201 + mrSges(6,2) * t198;
t396 = m(6) * pkin(4) + mrSges(5,1) - t251 - t253 + t452;
t395 = -m(6) * pkin(10) + mrSges(5,2) - t445;
t112 = -mrSges(7,1) * t163 + mrSges(7,2) * t164;
t61 = -qJD(4) * pkin(4) - t66;
t45 = -pkin(5) * t163 + qJD(6) + t61;
t394 = -m(7) * t45 - t112;
t392 = m(4) + m(3) - t398;
t391 = -m(6) * t258 + t199 * t445 + t202 * t452 + mrSges(4,1) - t254;
t390 = -m(6) * t61 - t416;
t1 = pkin(5) * t160 - qJ(6) * t104 - qJD(6) * t164 + t4;
t78 = mrSges(7,1) * t160 - mrSges(7,3) * t104;
t79 = mrSges(6,1) * t160 - mrSges(6,3) * t104;
t389 = m(6) * t4 + m(7) * t1 + t78 + t79;
t2 = qJ(6) * t105 + qJD(6) * t163 + t3;
t80 = -mrSges(7,2) * t160 + mrSges(7,3) * t105;
t81 = -mrSges(6,2) * t160 + mrSges(6,3) * t105;
t388 = m(6) * t3 + m(7) * t2 + t80 + t81;
t358 = mrSges(7,3) * t163;
t123 = mrSges(7,2) * t438 + t358;
t360 = mrSges(6,3) * t163;
t124 = mrSges(6,2) * t438 + t360;
t19 = qJ(6) * t163 + t21;
t387 = m(6) * t21 + m(7) * t19 + t123 + t124;
t349 = t164 * mrSges(7,3);
t125 = -mrSges(7,1) * t438 - t349;
t359 = mrSges(6,3) * t164;
t126 = -mrSges(6,1) * t438 - t359;
t20 = -t198 * t62 + t201 * t82;
t18 = -qJ(6) * t164 + t20;
t17 = -pkin(5) * t438 + t18;
t386 = m(6) * t20 + m(7) * t17 + t125 + t126;
t385 = -t4 * mrSges(6,1) + t3 * mrSges(6,2) + t2 * mrSges(7,2);
t384 = t390 + t394;
t10 = -qJDD(4) * pkin(4) - t12;
t56 = -t105 * mrSges(7,1) + t104 * mrSges(7,2);
t7 = -pkin(5) * t105 + qJDD(6) + t10;
t383 = -m(5) * t12 + m(6) * t10 + m(7) * t7 + t417 + t56;
t382 = -m(5) * t66 - t384;
t203 = qJD(3) ^ 2;
t376 = t163 / 0.2e1;
t374 = t164 / 0.2e1;
t352 = t10 * t199;
t139 = t193 * t256 + t196 * t339;
t72 = t139 * t368 - t200 * t415;
t346 = t198 * t72;
t140 = -t193 * t255 + t196 * t340;
t74 = t140 * t368 - t200 * t414;
t345 = t198 * t74;
t343 = t202 * Ifges(5,4);
t337 = qJD(3) * mrSges(4,2);
t276 = t194 * t342;
t111 = t200 * t276 + t211;
t336 = t111 * t198;
t333 = t193 * t200;
t331 = t198 * t199;
t324 = mrSges(4,2) * qJDD(3);
t316 = qJD(5) * t199;
t310 = pkin(9) * t318;
t305 = mrSges(5,3) * t321;
t296 = t198 * t318;
t272 = t313 / 0.2e1;
t252 = mrSges(6,1) * t198 + mrSges(6,2) * t201;
t250 = mrSges(7,1) * t198 + mrSges(7,2) * t201;
t240 = Ifges(5,5) * t202 - Ifges(5,6) * t199;
t234 = t368 * t276;
t110 = t195 * t333 - t227 - t234;
t219 = t341 * t342 - t304;
t77 = t111 * t202 + t199 * t219;
t37 = t110 * t201 - t198 * t77;
t38 = t110 * t198 + t201 * t77;
t87 = -qJD(3) * pkin(3) + t206;
t229 = t87 * (mrSges(5,1) * t199 + mrSges(5,2) * t202);
t228 = t199 * (Ifges(5,1) * t202 - t357);
t142 = t199 * t341 + t202 * t332;
t116 = -t198 * t142 - t201 * t300;
t224 = -t201 * t142 + t198 * t300;
t223 = -t198 * t316 + t201 * t318;
t222 = t199 * t315 + t296;
t141 = t199 * t332 - t202 * t341;
t76 = t111 * t199 - t202 * t219;
t205 = t194 * t210 + t274 * t341;
t204 = t194 * t209 - t275 * t341;
t182 = t197 * t201;
t181 = t197 * t198;
t180 = -qJD(4) * mrSges(5,2) + t305;
t174 = (pkin(9) + t366) * t199;
t162 = t201 * t178;
t135 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t172;
t132 = -pkin(9) * t330 + t162;
t130 = pkin(5) * t222 + t310;
t118 = -qJ(6) * t331 + t133;
t114 = -qJD(4) * t141 + t202 * t269;
t109 = -qJ(6) * t329 + t162 + (-pkin(9) * t198 - pkin(5)) * t202;
t103 = t111 * qJD(3);
t102 = (t234 + (t233 - t333) * t195) * qJD(3);
t73 = t140 * t200 + t368 * t414;
t71 = t139 * t200 + t368 * t415;
t36 = t199 * t205 + t74 * t202;
t35 = t199 * t74 - t202 * t205;
t34 = t199 * t204 + t72 * t202;
t33 = t199 * t72 - t202 * t204;
t32 = -qJD(4) * t76 + t102 * t202;
t5 = [m(4) * (t90 * t102 + t52 * t111 + t129 * t219) + m(5) * (t11 * t77 + t32 * t67) + t32 * t180 + t77 * t135 - t111 * t324 - t102 * t337 + m(3) * (t185 * t342 + (t193 ^ 2 + t196 ^ 2) * t195 ^ 2 * qJDD(1)) + m(2) * qJDD(1) + t387 * (qJD(5) * t37 + t103 * t198 + t201 * t32) + t386 * (-qJD(5) * t38 + t103 * t201 - t198 * t32) + t388 * t38 + t389 * t37 + (-m(4) * t53 + m(5) * t51 - t402) * t110 + (m(4) * t206 + m(5) * t87 - t401) * t103 + t383 * t76 + t382 * (qJD(4) * t77 + t102 * t199) + (-m(2) - t392) * g(3); m(3) * t185 + t165 * t298 + t114 * t180 + t142 * t135 + m(5) * (t11 * t142 + t67 * t114 + (t322 * t87 - t368 * t51) * t194) + m(4) * (t129 * t341 + (t368 * t53 + t200 * t52 + (t200 * t206 + t368 * t90) * qJD(3)) * t194) + t386 * (qJD(5) * t224 - t198 * t114 + t201 * t298) + t387 * (qJD(5) * t116 + t201 * t114 + t198 * t298) + (-mrSges(4,1) * t203 - t324) * t332 - t388 * t224 + t389 * t116 + (-mrSges(4,2) * t203 + t402) * t300 + t383 * t141 + t382 * (qJD(4) * t142 + t199 * t269) + (-g(1) * t274 + g(2) * t275 - g(3) * t342) * t392; (t20 * mrSges(6,1) + t17 * mrSges(7,1) - t21 * mrSges(6,2) - t19 * mrSges(7,2) - t67 * mrSges(5,3) - t437) * t320 + t45 * (mrSges(7,1) * t222 + mrSges(7,2) * t223) + t61 * (mrSges(6,1) * t222 + mrSges(6,2) * t223) + t51 * t254 - (t198 * t422 + t201 * t423) * t316 / 0.2e1 + (t180 * t202 + (-t199 * t66 + t202 * t67) * m(5) - t384 * t199 - t337) * t206 + t343 * t435 + t245 * t436 + (-t20 * t223 - t21 * t222 - t3 * t331 - t329 * t4) * mrSges(6,3) + t417 * pkin(9) * t199 + (t240 * t451 + t413 * t372 + t411 * t374 + t412 * t376 + t229) * qJD(4) + (Ifges(5,4) * t435 + Ifges(5,2) * t436 + (-Ifges(5,2) * t199 + t343) * t272 + pkin(9) * t135 - t1 * mrSges(7,1) - t432 * t380 - t447 * t379 + t446 * t378 + t385) * t202 + (Ifges(5,1) * t173 + Ifges(5,4) * t436 + t7 * t250 + t378 * t408 + t379 * t406 + t380 * t404) * t199 + (-t346 * t381 + mrSges(4,2) * t72 + t398 * (-t71 * pkin(3) + pkin(9) * t72) + t434 * (-t328 * t71 + t346) - t433 * (t201 * t72 + t330 * t71) + t391 * t71) * g(2) + (-t345 * t381 + mrSges(4,2) * t74 + t398 * (-t73 * pkin(3) + pkin(9) * t74) + t434 * (-t328 * t73 + t345) - t433 * (t201 * t74 + t330 * t73) + t391 * t73) * g(1) + (-t336 * t381 + mrSges(4,2) * t111 + t434 * (-t110 * t328 + t336) - t433 * (t110 * t330 + t111 * t201) + t398 * (-t110 * pkin(3) + pkin(9) * t111) + t391 * t110) * g(3) + t401 * t90 - t403 * t202 / 0.2e1 + (-t372 * t409 - t374 * t405 - t376 * t407) * t316 - t423 * t296 / 0.2e1 + t425 * t126 + t426 * t124 + (t132 * t4 + t133 * t3 + (t318 * t61 + t352) * pkin(9) + t426 * t21 + t425 * t20) * m(6) + t427 * t123 + t428 * t125 + (t1 * t109 + t118 * t2 + t130 * t45 + t17 * t428 + t174 * t7 + t19 * t427) * m(7) - t430 * t331 / 0.2e1 + (-g(1) * t74 - g(2) * t72 - g(3) * t111 - t318 * t66 + t400) * mrSges(5,3) + t439 * t318 / 0.2e1 + (-t87 * t90 - pkin(3) * t51 + ((-t199 * t67 - t202 * t66) * qJD(4) + t400) * pkin(9)) * m(5) + (-t1 * t329 - t17 * t223 - t19 * t222 - t2 * t331) * mrSges(7,3) + t416 * t310 + t329 * t444 + Ifges(4,3) * qJDD(3) + t252 * t352 + qJDD(4) * (Ifges(5,5) * t199 + Ifges(5,6) * t202) + t174 * t56 + t130 * t112 + t132 * t79 + t133 * t81 + t118 * t80 - pkin(3) * t119 + t109 * t78 - t180 * t311 + t228 * t272 - t52 * mrSges(4,2) + t53 * mrSges(4,1); (-t1 * t198 - t17 * t315 - t19 * t317 + t2 * t201) * mrSges(7,3) + t437 * t323 + (t395 * t77 + t396 * t76) * g(3) + (t35 * t396 + t36 * t395) * g(1) + (t33 * t396 + t34 * t395) * g(2) + t7 * t251 + t10 * t253 + (t297 / 0.2e1 - t317 / 0.2e1) * t423 + (-pkin(4) * t10 - t20 * t39 - t21 * t40) * m(6) + (t306 + t390) * t67 + (-t124 * t317 - t126 * t315 + t201 * t81 - t198 * t79 + m(6) * ((-t21 * t198 - t20 * t201) * qJD(5) + t399)) * pkin(10) + (-t20 * t315 - t21 * t317 + t399) * mrSges(6,3) + t405 * t380 + t407 * t379 + t409 * t378 + (t305 - t180) * t66 + t422 * t315 / 0.2e1 + t430 * t201 / 0.2e1 - (-Ifges(5,2) * t323 + t191 + t439) * t321 / 0.2e1 + t438 * (-t45 * t250 - t61 * t252) + t198 * t444 + t419 * t112 + t420 * t125 + (t1 * t181 + t17 * t420 - t182 * t2 + t19 * t421 - t190 * t7 + t419 * t45) * m(7) + t421 * t123 + (t163 * t406 + t164 * t404 - t408 * t438) * qJD(5) / 0.2e1 + Ifges(5,3) * qJDD(4) - t190 * t56 + Ifges(5,5) * t173 + t181 * t78 - t182 * t80 + Ifges(5,6) * t172 - t39 * t126 - t40 * t124 - (t163 * t412 + t164 * t411 - t413 * t438) * qJD(3) / 0.2e1 - t240 * t313 / 0.2e1 - t11 * mrSges(5,2) + t12 * mrSges(5,1) + (-t229 - t19 * (-mrSges(7,2) * t199 - mrSges(7,3) * t330) - t21 * (-mrSges(6,2) * t199 - mrSges(6,3) * t330) - t17 * (mrSges(7,1) * t199 - mrSges(7,3) * t328) - t20 * (mrSges(6,1) * t199 - mrSges(6,3) * t328)) * qJD(3) - t203 * t228 / 0.2e1 - pkin(4) * t57; (-t433 * (-t198 * t71 - t201 * t34) + t397 * (-t198 * t34 + t201 * t71)) * g(2) + (t349 + t125 - m(7) * (-t17 + t18)) * t19 + t418 * t1 - t385 + (t360 - t124) * t20 + (t37 * t397 + t38 * t433) * g(3) + t423 * t374 + (-t433 * (-t198 * t73 - t201 * t36) + t397 * (-t198 * t36 + t201 * t73)) * g(1) + (t163 * t432 - t164 * t447) * t438 / 0.2e1 + (t163 * t450 - t440) * t453 + (-t164 * t448 + t422 + t441) * t454 + (t359 + t126) * t21 + t17 * t358 - t45 * (mrSges(7,1) * t164 + mrSges(7,2) * t163) - t61 * (mrSges(6,1) * t164 + mrSges(6,2) * t163) - t18 * t123 + t403 + (t164 * t394 + t78) * pkin(5); -t163 * t123 + t164 * t125 + (-g(1) * t35 - g(2) * t33 - g(3) * t76 - t19 * t163 + t17 * t164 + t7) * m(7) + t56;];
tau  = t5;
