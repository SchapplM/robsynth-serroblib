% Calculate vector of inverse dynamics joint torques for
% S5PRRRR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:17
% EndTime: 2019-12-05 17:24:06
% DurationCPUTime: 22.46s
% Computational Cost: add. (7481->683), mult. (19360->1009), div. (0->0), fcn. (16188->14), ass. (0->327)
t236 = sin(pkin(5));
t237 = cos(pkin(6));
t245 = cos(qJ(3));
t246 = cos(qJ(2));
t351 = t245 * t246;
t241 = sin(qJ(3));
t242 = sin(qJ(2));
t355 = t241 * t242;
t269 = -t237 * t355 + t351;
t166 = t269 * t236;
t148 = qJD(1) * t166;
t235 = sin(pkin(6));
t365 = t235 * t241;
t230 = pkin(8) * t365;
t358 = t237 * t245;
t202 = pkin(2) * t358 - t230;
t189 = t202 * qJD(3);
t454 = -t148 + t189;
t359 = t237 * t241;
t362 = t235 * t245;
t203 = pkin(2) * t359 + pkin(8) * t362;
t176 = pkin(9) * t237 + t203;
t292 = -pkin(3) * t245 - pkin(9) * t241;
t177 = (-pkin(2) + t292) * t235;
t274 = (pkin(3) * t241 - pkin(9) * t245) * t235;
t188 = qJD(3) * t274;
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t346 = qJD(1) * t236;
t322 = t242 * t346;
t300 = t235 * t322;
t339 = qJD(4) * t244;
t340 = qJD(4) * t240;
t456 = -t176 * t340 + t177 * t339 + t454 * t244 + (t188 - t300) * t240;
t353 = t242 * t245;
t354 = t241 * t246;
t271 = t237 * t353 + t354;
t257 = t271 * t236;
t455 = qJD(1) * t257 - t203 * qJD(3);
t239 = sin(qJ(5));
t243 = cos(qJ(5));
t343 = qJD(2) * t245;
t320 = t235 * t343;
t218 = qJD(4) - t320;
t214 = qJD(2) * pkin(2) + t246 * t346;
t238 = cos(pkin(5));
t345 = qJD(1) * t238;
t225 = t237 * t345;
t121 = t225 + (qJD(2) * t292 - t214) * t235;
t266 = t214 * t237 + t235 * t345;
t344 = qJD(2) * t235;
t204 = pkin(8) * t344 + t322;
t368 = t204 * t245;
t109 = t241 * t266 + t368;
t229 = qJD(2) * t237 + qJD(3);
t93 = pkin(9) * t229 + t109;
t46 = t121 * t240 + t244 * t93;
t41 = pkin(10) * t218 + t46;
t321 = t241 * t344;
t163 = t229 * t244 - t240 * t321;
t164 = t229 * t240 + t244 * t321;
t108 = -t241 * t204 + t245 * t266;
t92 = -pkin(3) * t229 - t108;
t47 = -pkin(4) * t163 - pkin(10) * t164 + t92;
t17 = t239 * t47 + t243 * t41;
t442 = t17 * mrSges(6,2);
t16 = -t239 * t41 + t243 * t47;
t443 = t16 * mrSges(6,1);
t384 = Ifges(5,4) * t164;
t457 = Ifges(5,6) * t218;
t458 = Ifges(5,2) * t163;
t81 = t384 + t457 + t458;
t469 = -t442 + t443 - t81 / 0.2e1;
t398 = pkin(9) * t244;
t342 = qJD(3) * t241;
t318 = t235 * t342;
t468 = -pkin(10) * t318 - t456;
t200 = -t244 * t237 + t240 * t365;
t341 = qJD(3) * t245;
t317 = t235 * t341;
t138 = -qJD(4) * t200 + t244 * t317;
t363 = t235 * t244;
t201 = t237 * t240 + t241 * t363;
t139 = qJD(4) * t201 + t240 * t317;
t467 = pkin(4) * t139 - pkin(10) * t138 - t455;
t124 = -t164 * t239 + t218 * t243;
t125 = t164 * t243 + t218 * t239;
t156 = qJD(5) - t163;
t42 = Ifges(6,5) * t125 + Ifges(6,6) * t124 + Ifges(6,3) * t156;
t466 = -mrSges(5,3) * t46 + t42 / 0.2e1 + t469;
t288 = mrSges(6,1) * t239 + mrSges(6,2) * t243;
t446 = m(5) + m(6);
t465 = -pkin(9) * t446 + mrSges(4,2) - mrSges(5,3) - t288;
t228 = qJDD(2) * t237 + qJDD(3);
t315 = qJD(2) * t346;
t332 = qJDD(1) * t236;
t195 = -t242 * t315 + t246 * t332;
t171 = qJDD(2) * pkin(2) + t195;
t331 = qJDD(1) * t238;
t314 = t235 * t331;
t196 = t242 * t332 + t246 * t315;
t448 = pkin(8) * qJDD(2) * t235 + qJD(3) * t266 + t196;
t37 = t245 * (t171 * t237 + t314) - t204 * t341 - t448 * t241;
t32 = -pkin(3) * t228 - t37;
t333 = qJD(2) * qJD(3);
t194 = (qJDD(2) * t241 + t245 * t333) * t235;
t90 = qJD(4) * t163 + t194 * t244 + t228 * t240;
t91 = -qJD(4) * t164 - t194 * t240 + t228 * t244;
t15 = -pkin(4) * t91 - pkin(10) * t90 + t32;
t193 = (-qJDD(2) * t245 + t241 * t333) * t235;
t180 = qJDD(4) + t193;
t36 = t171 * t359 - t204 * t342 + t241 * t314 + t245 * t448;
t31 = pkin(9) * t228 + t36;
t130 = -t171 * t235 + t237 * t331;
t69 = pkin(3) * t193 - pkin(9) * t194 + t130;
t8 = t121 * t339 + t240 * t69 + t244 * t31 - t340 * t93;
t3 = pkin(10) * t180 + t8;
t2 = -qJD(5) * t17 + t15 * t243 - t239 * t3;
t444 = t2 * mrSges(6,1);
t1 = qJD(5) * t16 + t15 * t239 + t243 * t3;
t445 = t1 * mrSges(6,2);
t464 = t444 - t445;
t187 = qJD(2) * t274;
t71 = t244 * t108 + t240 * t187;
t463 = t340 * pkin(9) + pkin(10) * t321 + t71;
t291 = pkin(4) * t240 - pkin(10) * t244;
t462 = -qJD(5) * t398 + t291 * qJD(4) - t214 * t359 - t368 - (t241 * t345 + t291 * t343) * t235;
t175 = t230 + (-pkin(2) * t245 - pkin(3)) * t237;
t100 = pkin(4) * t200 - pkin(10) * t201 + t175;
t433 = t244 * t176 + t240 * t177;
t102 = -pkin(10) * t362 + t433;
t38 = t100 * t243 - t102 * t239;
t461 = qJD(5) * t38 + t239 * t467 - t243 * t468;
t39 = t100 * t239 + t102 * t243;
t460 = -qJD(5) * t39 + t239 * t468 + t243 * t467;
t154 = Ifges(5,4) * t163;
t459 = Ifges(5,5) * t218;
t303 = mrSges(4,3) * t321;
t453 = mrSges(4,1) * t229 + mrSges(5,1) * t163 - mrSges(5,2) * t164 - t303;
t352 = t244 * t245;
t150 = (-t239 * t352 + t241 * t243) * t344;
t452 = t239 * t339 + t150;
t451 = t37 * mrSges(4,1) - t36 * mrSges(4,2);
t289 = -mrSges(6,1) * t243 + mrSges(6,2) * t239;
t261 = m(6) * pkin(4) - t289;
t290 = -mrSges(5,1) * t244 + mrSges(5,2) * t240;
t328 = m(6) * pkin(10) + mrSges(6,3);
t450 = pkin(3) * t446 + t240 * t328 + t244 * t261 + mrSges(4,1) - t290;
t401 = t180 / 0.2e1;
t412 = t91 / 0.2e1;
t413 = t90 / 0.2e1;
t423 = Ifges(5,1) * t413 + Ifges(5,4) * t412 + Ifges(5,5) * t401;
t33 = qJD(5) * t124 + t180 * t239 + t243 * t90;
t422 = t33 / 0.2e1;
t34 = -qJD(5) * t125 + t180 * t243 - t239 * t90;
t421 = t34 / 0.2e1;
t87 = qJDD(5) - t91;
t414 = t87 / 0.2e1;
t14 = -mrSges(6,1) * t34 + mrSges(6,2) * t33;
t61 = mrSges(5,1) * t180 - mrSges(5,3) * t90;
t440 = t14 - t61;
t439 = mrSges(5,1) + t261;
t304 = mrSges(5,2) - t328;
t217 = -pkin(4) * t244 - pkin(10) * t240 - pkin(3);
t336 = qJD(5) * t243;
t438 = t217 * t336 + t239 * t462 - t243 * t463;
t338 = qJD(5) * t239;
t437 = -t217 * t338 + t239 * t463 + t243 * t462;
t388 = mrSges(5,3) * t164;
t129 = mrSges(5,1) * t218 - t388;
t58 = -mrSges(6,1) * t124 + mrSges(6,2) * t125;
t373 = t129 - t58;
t435 = t240 * t336 + t452;
t151 = (t239 * t241 + t243 * t352) * t344;
t337 = qJD(5) * t240;
t434 = t239 * t337 - t243 * t339 + t151;
t432 = t237 * t351 - t355;
t431 = t1 * t243 - t2 * t239;
t9 = -qJD(4) * t46 - t240 * t31 + t244 * t69;
t428 = t9 * mrSges(5,1) - t8 * mrSges(5,2) + Ifges(5,5) * t90 + Ifges(5,6) * t91 + Ifges(5,3) * t180;
t57 = -qJD(4) * t433 + t188 * t244 - t189 * t240;
t247 = qJD(2) ^ 2;
t6 = Ifges(6,4) * t33 + Ifges(6,2) * t34 + Ifges(6,6) * t87;
t425 = t6 / 0.2e1;
t424 = Ifges(6,1) * t422 + Ifges(6,4) * t421 + Ifges(6,5) * t414;
t381 = Ifges(6,4) * t125;
t43 = Ifges(6,2) * t124 + Ifges(6,6) * t156 + t381;
t419 = -t43 / 0.2e1;
t418 = t43 / 0.2e1;
t122 = Ifges(6,4) * t124;
t44 = Ifges(6,1) * t125 + Ifges(6,5) * t156 + t122;
t417 = -t44 / 0.2e1;
t416 = t44 / 0.2e1;
t410 = -t124 / 0.2e1;
t409 = t124 / 0.2e1;
t408 = -t125 / 0.2e1;
t407 = t125 / 0.2e1;
t406 = -t156 / 0.2e1;
t405 = t156 / 0.2e1;
t402 = t164 / 0.2e1;
t400 = t237 / 0.2e1;
t4 = -pkin(4) * t180 - t9;
t395 = t240 * t4;
t394 = t244 * t8;
t389 = mrSges(5,3) * t163;
t387 = mrSges(5,3) * t240;
t386 = Ifges(4,4) * t241;
t385 = Ifges(4,4) * t245;
t383 = Ifges(5,4) * t240;
t382 = Ifges(5,4) * t244;
t380 = Ifges(6,4) * t239;
t379 = Ifges(6,4) * t243;
t378 = t163 * Ifges(5,6);
t377 = t164 * Ifges(5,5);
t376 = t218 * Ifges(5,3);
t375 = t229 * Ifges(4,5);
t374 = t229 * Ifges(4,6);
t372 = cos(pkin(11));
t371 = sin(pkin(11));
t370 = t163 * t239;
t369 = t163 * t243;
t366 = t235 * t240;
t364 = t235 * t242;
t361 = t236 * t246;
t360 = t237 * t238;
t357 = t239 * t240;
t356 = t240 * t243;
t327 = t236 * t364;
t347 = pkin(2) * t361 + pkin(8) * t327;
t334 = -m(4) - t446;
t5 = Ifges(6,5) * t33 + Ifges(6,6) * t34 + Ifges(6,3) * t87;
t329 = pkin(9) * t339;
t324 = Ifges(4,5) * t194 - Ifges(4,6) * t193 + Ifges(4,3) * t228;
t312 = t339 / 0.2e1;
t311 = -t337 / 0.2e1;
t310 = t236 * t372;
t309 = t236 * t371;
t308 = t372 * t242;
t307 = t372 * t246;
t306 = t371 * t242;
t305 = t371 * t246;
t302 = mrSges(4,3) * t320;
t299 = qJD(2) * t327;
t294 = t235 * t310;
t293 = t235 * t309;
t287 = Ifges(5,1) * t244 - t383;
t286 = Ifges(6,1) * t243 - t380;
t285 = Ifges(6,1) * t239 + t379;
t284 = -Ifges(5,2) * t240 + t382;
t283 = -Ifges(6,2) * t239 + t379;
t282 = Ifges(6,2) * t243 + t380;
t281 = Ifges(5,5) * t244 - Ifges(5,6) * t240;
t280 = Ifges(6,5) * t243 - Ifges(6,6) * t239;
t279 = Ifges(6,5) * t239 + Ifges(6,6) * t243;
t45 = t121 * t244 - t240 * t93;
t132 = -t236 * t432 - t238 * t362;
t270 = t237 * t354 + t353;
t133 = t236 * t270 + t238 * t365;
t272 = -t235 * t361 + t360;
t99 = t133 * t244 + t240 * t272;
t50 = t132 * t243 - t239 * t99;
t51 = t132 * t239 + t243 * t99;
t70 = -t108 * t240 + t187 * t244;
t110 = -t176 * t240 + t177 * t244;
t140 = -t201 * t239 - t243 * t362;
t273 = -t201 * t243 + t239 * t362;
t265 = t241 * (Ifges(4,1) * t245 - t386);
t264 = (Ifges(4,2) * t245 + t386) * t235;
t255 = mrSges(3,2) + (pkin(8) * t334 - mrSges(4,3)) * t235;
t254 = t238 * t305 + t308;
t253 = -t238 * t307 + t306;
t251 = t254 * t245;
t250 = t253 * t241;
t249 = t253 * t245;
t98 = t133 * t240 - t244 * t272;
t223 = Ifges(4,4) * t320;
t199 = -t238 * t306 + t307;
t198 = t238 * t308 + t305;
t192 = t254 * pkin(2);
t191 = t253 * pkin(2);
t186 = (-mrSges(4,1) * t245 + mrSges(4,2) * t241) * t344;
t185 = -mrSges(4,2) * t229 + t302;
t173 = t217 * t239 + t243 * t398;
t172 = t217 * t243 - t239 * t398;
t153 = -t214 * t235 + t225;
t143 = Ifges(4,1) * t321 + t223 + t375;
t142 = qJD(2) * t264 + t374;
t137 = mrSges(4,1) * t228 - mrSges(4,3) * t194;
t136 = -mrSges(4,2) * t228 - mrSges(4,3) * t193;
t135 = t235 * t254 + t237 * t309;
t134 = t235 * t253 - t237 * t310;
t128 = -mrSges(5,2) * t218 + t389;
t123 = mrSges(4,1) * t193 + mrSges(4,2) * t194;
t119 = t148 * t240 - t244 * t300;
t118 = -t199 * t359 - t251;
t116 = -t198 * t359 - t249;
t107 = pkin(4) * t164 - pkin(10) * t163;
t101 = pkin(4) * t362 - t110;
t97 = t199 * t245 + (-t237 * t254 + t293) * t241;
t96 = t199 * t241 + t237 * t251 - t245 * t293;
t95 = t198 * t245 - t237 * t250 - t241 * t294;
t94 = t198 * t241 + t237 * t249 + t245 * t294;
t84 = t238 * t317 + (t269 * qJD(2) + qJD(3) * t432) * t236;
t83 = t238 * t318 + (qJD(2) * t271 + qJD(3) * t270) * t236;
t82 = Ifges(5,1) * t164 + t154 + t459;
t80 = t376 + t377 + t378;
t78 = mrSges(6,1) * t156 - mrSges(6,3) * t125;
t77 = -mrSges(6,2) * t156 + mrSges(6,3) * t124;
t66 = qJD(5) * t273 - t138 * t239 + t243 * t318;
t65 = qJD(5) * t140 + t138 * t243 + t239 * t318;
t62 = -mrSges(5,2) * t180 + mrSges(5,3) * t91;
t59 = -pkin(4) * t321 - t70;
t55 = t135 * t240 + t244 * t97;
t53 = t134 * t240 + t244 * t95;
t49 = -pkin(4) * t318 - t57;
t40 = -pkin(4) * t218 - t45;
t35 = -mrSges(5,1) * t91 + mrSges(5,2) * t90;
t27 = -qJD(4) * t98 + t240 * t299 + t84 * t244;
t26 = qJD(4) * t99 + t84 * t240 - t244 * t299;
t24 = t90 * Ifges(5,4) + t91 * Ifges(5,2) + t180 * Ifges(5,6);
t21 = t107 * t239 + t243 * t45;
t20 = t107 * t243 - t239 * t45;
t19 = -mrSges(6,2) * t87 + mrSges(6,3) * t34;
t18 = mrSges(6,1) * t87 - mrSges(6,3) * t33;
t13 = qJD(5) * t50 + t239 * t83 + t243 * t27;
t12 = -qJD(5) * t51 - t239 * t27 + t243 * t83;
t7 = [t123 * t360 + m(2) * qJDD(1) + t12 * t78 + t27 * t128 + t13 * t77 + t133 * t136 + t50 * t18 + t84 * t185 + t51 * t19 + t99 * t62 + t440 * t98 - t453 * t83 - t373 * t26 + (-t137 + t35) * t132 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t247 - t123 * t235) * t246 + (-mrSges(3,1) * t247 - mrSges(3,2) * qJDD(2) + t186 * t344) * t242) * t236 + (-m(2) - m(3) + t334) * g(3) + m(6) * (t1 * t51 + t12 * t16 + t13 * t17 + t2 * t50 + t26 * t40 + t4 * t98) + m(5) * (t132 * t32 - t26 * t45 + t27 * t46 + t8 * t99 + t83 * t92 - t9 * t98) + m(3) * (qJDD(1) * t238 ^ 2 + (t195 * t246 + t196 * t242) * t236) + m(4) * (-t108 * t83 + t109 * t84 + t130 * t272 - t37 * t132 + t36 * t133 + t153 * t299); (-pkin(2) * t123 + (mrSges(4,2) * t130 - mrSges(4,3) * t37) * t241 + (-t130 * mrSges(4,1) + t36 * mrSges(4,3) - t428) * t245) * t235 + (t110 * t9 + t433 * t8 + t175 * t32 - t455 * t92 + t456 * t46 + (t119 + t57) * t45) * m(5) + t456 * t128 + t451 * t237 + t454 * t185 + t453 * t455 + (-pkin(2) * t130 * t235 + g(1) * t192 + g(2) * t191 - g(3) * t347 + t108 * t455 + t109 * t454 - t153 * t300 + t202 * t37 + t203 * t36) * m(4) + (mrSges(5,2) * t32 - mrSges(5,3) * t9 + 0.2e1 * t423) * t201 + (Ifges(6,3) * t414 + Ifges(6,6) * t421 + Ifges(6,5) * t422 - t8 * mrSges(5,3) + t5 / 0.2e1 - t24 / 0.2e1 + t32 * mrSges(5,1) - Ifges(5,6) * t401 - Ifges(5,2) * t412 - Ifges(5,4) * t413 + t464) * t200 + (t82 / 0.2e1 + t92 * mrSges(5,2) - t45 * mrSges(5,3) + t154 / 0.2e1 + Ifges(5,1) * t402 + t459 / 0.2e1) * t138 + t460 * t78 + (t1 * t39 + t101 * t4 + t2 * t38 + (-t119 + t49) * t40 + t461 * t17 + t460 * t16) * m(6) + t461 * t77 + (Ifges(4,3) * t400 + (Ifges(4,5) * t241 + Ifges(4,6) * t245) * t235) * t228 - (Ifges(4,6) * t400 + t264) * t193 + (Ifges(4,5) * t400 + (t241 * Ifges(4,1) + t385) * t235) * t194 + t57 * t129 + t110 * t61 + t101 * t14 + ((t143 / 0.2e1 + t153 * mrSges(4,2) + t375 / 0.2e1) * t245 + (t45 * mrSges(5,1) - t46 * mrSges(5,2) + t378 / 0.2e1 + t377 / 0.2e1 + t376 / 0.2e1 + t153 * mrSges(4,1) - t374 / 0.2e1 + t80 / 0.2e1 - t142 / 0.2e1) * t241 + (t265 / 0.2e1 + t245 * (-Ifges(4,2) * t241 + t385) / 0.2e1) * t344 + (-t108 * t245 - t109 * t241) * mrSges(4,3)) * t235 * qJD(3) + t49 * t58 + t40 * (-mrSges(6,1) * t66 + mrSges(6,2) * t65) + t38 * t18 + t39 * t19 + (Ifges(6,5) * t65 + Ifges(6,6) * t66) * t405 + (Ifges(6,1) * t65 + Ifges(6,4) * t66) * t407 + (-t166 * mrSges(4,1) - t439 * (t166 * t244 + t240 * t327) + t465 * t257 + t304 * (t166 * t240 - t244 * t327) - t446 * (t166 * pkin(3) + t347)) * g(3) + (t254 * mrSges(3,1) - t118 * mrSges(4,1) + t304 * (t118 * t240 - t199 * t363) + t255 * t199 - t439 * (t118 * t244 + t199 * t366) + t465 * (t199 * t358 - t241 * t254) - t446 * (t118 * pkin(3) - t192)) * g(1) + (t253 * mrSges(3,1) - t116 * mrSges(4,1) + t304 * (t116 * t240 - t198 * t363) + t255 * t198 - t439 * (t116 * t244 + t198 * t366) + t465 * (t198 * t358 - t250) - t446 * (t116 * pkin(3) - t191)) * g(2) + (Ifges(6,4) * t65 + Ifges(6,2) * t66) * t409 + t373 * t119 + ((-mrSges(3,1) * t246 + (-t235 * mrSges(4,3) + mrSges(3,2)) * t242) * g(3) + (-t186 * t364 + (mrSges(3,1) * t242 + mrSges(3,2) * t246) * qJD(2)) * qJD(1)) * t236 + Ifges(3,3) * qJDD(2) + (t92 * mrSges(5,1) - t458 / 0.2e1 - Ifges(5,4) * t402 + Ifges(6,3) * t405 + Ifges(6,5) * t407 + Ifges(6,6) * t409 - t457 / 0.2e1 + t466) * t139 + t175 * t35 + t195 * mrSges(3,1) - t196 * mrSges(3,2) + t202 * t137 + t203 * t136 + t324 * t400 + t65 * t416 + t66 * t418 + t140 * t425 + t433 * t62 + (-Ifges(6,1) * t273 + Ifges(6,4) * t140) * t422 + (-Ifges(6,5) * t273 + Ifges(6,6) * t140) * t414 + (-Ifges(6,4) * t273 + Ifges(6,2) * t140) * t421 + t4 * (-mrSges(6,1) * t140 - mrSges(6,2) * t273) + (t1 * t140 - t16 * t65 + t17 * t66 + t2 * t273) * mrSges(6,3) - t273 * t424; t437 * t78 + (-t40 * t59 + t1 * t173 + t172 * t2 + (t339 * t40 + t395) * pkin(9) + t438 * t17 + t437 * t16) * m(6) + t438 * t77 + (mrSges(6,1) * t435 - mrSges(6,2) * t434) * t40 + (-t1 * t357 + t16 * t434 - t17 * t435 - t2 * t356) * mrSges(6,3) + t440 * pkin(9) * t240 + t218 * t92 * (mrSges(5,1) * t240 + mrSges(5,2) * t244) + t452 * t419 + t451 + (Ifges(6,5) * t408 + Ifges(6,6) * t410 + Ifges(6,3) * t406 - t469) * t240 * t320 + (-m(5) * t92 + t303 + t453) * t109 - t6 * t357 / 0.2e1 + (Ifges(6,1) * t151 + Ifges(6,4) * t150) * t408 + (Ifges(6,4) * t151 + Ifges(6,2) * t150) * t410 - t71 * t128 - t235 ^ 2 * t247 * t265 / 0.2e1 - pkin(3) * t35 + (-t70 - t329) * t129 + (t163 * t284 + t164 * t287 + t218 * t281) * qJD(4) / 0.2e1 + t324 + (-pkin(3) * t32 + (-t240 * t9 + t394 + (-t240 * t46 - t244 * t45) * qJD(4)) * pkin(9) - t45 * t70 - t46 * t71) * m(5) + (Ifges(6,5) * t151 + Ifges(6,6) * t150) * t406 + (t132 * t450 + t133 * t465) * g(3) + (t450 * t94 + t465 * t95) * g(2) + (t450 * t96 + t465 * t97) * g(1) + (t302 - t185) * t108 + t82 * t312 + t32 * t290 + (-t339 * t45 + t394) * mrSges(5,3) + t239 * t44 * t311 + (-t46 * (-mrSges(5,2) * t241 - t245 * t387) - t153 * (mrSges(4,1) * t241 + mrSges(4,2) * t245) - t45 * (mrSges(5,1) * t241 - mrSges(5,3) * t352)) * t344 + (t311 * t43 + t312 * t44) * t243 - t9 * t387 + (-t59 + t329) * t58 + t142 * t321 / 0.2e1 + t244 * t445 + (-pkin(9) * t128 + t466) * t340 + t172 * t18 + t173 * t19 - t244 * t444 + t288 * t395 + t62 * t398 + (Ifges(5,5) * t240 + Ifges(5,6) * t244) * t401 + (-t279 * t337 + (Ifges(6,3) * t240 + t244 * t280) * qJD(4)) * t405 + (-t285 * t337 + (Ifges(6,5) * t240 + t244 * t286) * qJD(4)) * t407 + (-t282 * t337 + (Ifges(6,6) * t240 + t244 * t283) * qJD(4)) * t409 + (Ifges(5,2) * t244 + t383) * t412 + (Ifges(5,1) * t240 + t382) * t413 + (-Ifges(6,3) * t244 + t240 * t280) * t414 + t151 * t417 + (-Ifges(6,6) * t244 + t240 * t283) * t421 + (-Ifges(6,5) * t244 + t240 * t286) * t422 + t240 * t423 + t356 * t424 + t244 * t24 / 0.2e1 - t244 * t5 / 0.2e1 - (t164 * (Ifges(5,5) * t241 + t245 * t287) + t163 * (Ifges(5,6) * t241 + t245 * t284) + t229 * (Ifges(4,5) * t245 - Ifges(4,6) * t241) + t241 * t80 + t218 * (Ifges(5,3) * t241 + t245 * t281) + (-Ifges(4,2) * t321 + t240 * t42 + t244 * t82 + t143 + t223) * t245) * t344 / 0.2e1; (-m(6) * t40 + t373 + t388) * t46 + (t304 * t99 + t439 * t98) * g(3) + (t304 * t53 - t439 * (t134 * t244 - t240 * t95)) * g(2) + ((-t338 + t370) * t17 + (-t336 + t369) * t16 + t431) * mrSges(6,3) + (-t77 * t338 - t78 * t336 + m(6) * ((-t16 * t243 - t17 * t239) * qJD(5) + t431) - t239 * t18 + t243 * t19) * pkin(10) + t156 * t40 * t288 - pkin(4) * t14 - t21 * t77 - t20 * t78 + (t124 * t283 + t125 * t286 + t156 * t280) * qJD(5) / 0.2e1 - (Ifges(5,1) * t163 - t384 + t42) * t164 / 0.2e1 - (-Ifges(5,2) * t164 + t154 + t82) * t163 / 0.2e1 + t4 * t289 + (-pkin(4) * t4 - t16 * t20 - t17 * t21) * m(6) + t428 + t164 * t442 + (t304 * t55 - t439 * (t135 * t244 - t240 * t97)) * g(1) - t92 * (mrSges(5,1) * t164 + mrSges(5,2) * t163) - t164 * t443 - t218 * (Ifges(5,5) * t163 - Ifges(5,6) * t164) / 0.2e1 + t81 * t402 + (Ifges(6,3) * t164 + t163 * t280) * t406 + (Ifges(6,5) * t164 + t163 * t286) * t408 + (Ifges(6,6) * t164 + t163 * t283) * t410 + t279 * t414 + t336 * t416 + t369 * t417 + t370 * t418 + t338 * t419 + t282 * t421 + t285 * t422 + t239 * t424 + t243 * t425 + (-t128 + t389) * t45; -t40 * (mrSges(6,1) * t125 + mrSges(6,2) * t124) + (Ifges(6,1) * t124 - t381) * t408 + t43 * t407 + (Ifges(6,5) * t124 - Ifges(6,6) * t125) * t406 - t16 * t77 + t17 * t78 - g(1) * ((-t239 * t55 + t243 * t96) * mrSges(6,1) + (-t239 * t96 - t243 * t55) * mrSges(6,2)) - g(2) * ((-t239 * t53 + t243 * t94) * mrSges(6,1) + (-t239 * t94 - t243 * t53) * mrSges(6,2)) - g(3) * (mrSges(6,1) * t50 - mrSges(6,2) * t51) + (t124 * t16 + t125 * t17) * mrSges(6,3) + t5 + (-Ifges(6,2) * t125 + t122 + t44) * t410 + t464;];
tau = t7;
