% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:13
% EndTime: 2019-03-08 20:45:43
% DurationCPUTime: 17.36s
% Computational Cost: add. (6899->677), mult. (14379->939), div. (0->0), fcn. (10176->14), ass. (0->319)
t217 = sin(qJ(4));
t458 = -t217 / 0.2e1;
t432 = m(6) + m(7);
t457 = m(5) + t432;
t221 = cos(qJ(4));
t269 = pkin(4) * t221 + pkin(9) * t217;
t168 = qJD(4) * t269 + qJD(3);
t277 = -pkin(9) * t221 + qJ(3);
t182 = pkin(4) * t217 + t277;
t216 = sin(qJ(5));
t220 = cos(qJ(5));
t222 = cos(qJ(2));
t224 = -pkin(2) - pkin(8);
t320 = qJD(4) * t224;
t287 = t221 * t320;
t334 = t217 * t224;
t298 = t216 * t334;
t318 = qJD(5) * t220;
t213 = sin(pkin(6));
t329 = qJD(1) * t213;
t218 = sin(qJ(2));
t333 = t218 * t220;
t423 = -qJD(5) * t298 + t216 * t168 + t182 * t318 + t220 * t287 - (t216 * t222 + t217 * t333) * t329;
t337 = t216 * t218;
t456 = t220 * t168 - (-t217 * t337 + t220 * t222) * t329;
t205 = pkin(5) * t220 + pkin(4);
t211 = qJ(5) + qJ(6);
t209 = sin(t211);
t210 = cos(t211);
t267 = -mrSges(6,1) * t220 + mrSges(6,2) * t216;
t455 = -m(6) * pkin(4) - m(7) * t205 - mrSges(7,1) * t210 + mrSges(7,2) * t209 + t267;
t208 = t217 * qJD(2);
t203 = t208 + qJD(5);
t196 = qJD(6) + t203;
t378 = -t196 / 0.2e1;
t322 = qJD(4) * t220;
t326 = qJD(2) * t221;
t171 = -t216 * t326 + t322;
t324 = qJD(4) * t216;
t172 = t220 * t326 + t324;
t215 = sin(qJ(6));
t219 = cos(qJ(6));
t94 = t171 * t215 + t172 * t219;
t385 = -t94 / 0.2e1;
t272 = t219 * t171 - t172 * t215;
t387 = -t272 / 0.2e1;
t454 = Ifges(7,5) * t385 + Ifges(7,6) * t387 + Ifges(7,3) * t378;
t383 = -m(4) - m(5);
t453 = -t171 / 0.2e1;
t452 = -t172 / 0.2e1;
t451 = -t203 / 0.2e1;
t360 = Ifges(5,4) * t221;
t450 = Ifges(5,2) * t458 + t360 / 0.2e1;
t195 = t220 * t334;
t278 = -t216 * t224 + pkin(5);
t335 = t217 * t220;
t310 = pkin(10) * t335;
t449 = (-t195 + (pkin(10) * t221 - t182) * t216) * qJD(5) + (t221 * t278 + t310) * qJD(4) + t456;
t323 = qJD(4) * t217;
t289 = t216 * t323;
t317 = qJD(5) * t221;
t235 = -t220 * t317 + t289;
t448 = -pkin(10) * t235 - t423;
t290 = t216 * t208;
t223 = -pkin(10) - pkin(9);
t296 = qJD(5) * t223;
t294 = t222 * t329;
t255 = qJD(3) - t294;
t154 = qJD(2) * t224 + t255;
t214 = cos(pkin(6));
t328 = qJD(1) * t214;
t113 = t154 * t221 - t217 * t328;
t176 = t269 * qJD(2);
t75 = t220 * t113 + t216 * t176;
t447 = -pkin(10) * t290 + t216 * t296 - t75;
t74 = -t113 * t216 + t220 * t176;
t446 = t220 * t296 - (pkin(5) * t221 + t310) * qJD(2) - t74;
t266 = t216 * mrSges(6,1) + t220 * mrSges(6,2);
t368 = mrSges(3,1) - mrSges(4,2);
t373 = pkin(5) * t216;
t445 = -m(7) * t373 - t209 * mrSges(7,1) - t210 * mrSges(7,2) - mrSges(5,3) - t266 - t368;
t444 = -m(7) * t223 + mrSges(6,3) + mrSges(7,3);
t295 = t218 * t329;
t125 = qJD(2) * t182 + t295;
t293 = t221 * t328;
t114 = t154 * t217 + t293;
t97 = qJD(4) * pkin(9) + t114;
t51 = t220 * t125 - t216 * t97;
t35 = -pkin(10) * t172 + t51;
t30 = pkin(5) * t203 + t35;
t52 = t125 * t216 + t220 * t97;
t36 = pkin(10) * t171 + t52;
t353 = t215 * t36;
t15 = t219 * t30 - t353;
t314 = qJDD(1) * t214;
t321 = qJD(4) * t221;
t327 = qJD(2) * t218;
t292 = t213 * t327;
t192 = qJD(1) * t292;
t338 = t213 * t222;
t140 = qJDD(1) * t338 - t192;
t246 = qJDD(3) - t140;
t438 = -qJD(4) * t328 + qJDD(2) * t224 + t246;
t47 = t154 * t321 + t217 * t438 + t221 * t314;
t43 = qJDD(4) * pkin(9) + t47;
t312 = qJDD(2) * qJ(3);
t313 = qJDD(1) * t218;
t119 = t213 * t313 + t312 + (qJD(3) + t294) * qJD(2);
t315 = qJD(2) * qJD(4);
t180 = qJDD(2) * t221 - t217 * t315;
t181 = -t217 * qJDD(2) - t221 * t315;
t73 = -pkin(4) * t181 - pkin(9) * t180 + t119;
t14 = -qJD(5) * t52 - t216 * t43 + t220 * t73;
t167 = qJDD(5) - t181;
t85 = qJD(5) * t171 + qJDD(4) * t216 + t180 * t220;
t6 = pkin(5) * t167 - pkin(10) * t85 + t14;
t319 = qJD(5) * t216;
t13 = t125 * t318 + t216 * t73 + t220 * t43 - t319 * t97;
t86 = -qJD(5) * t172 + qJDD(4) * t220 - t180 * t216;
t8 = pkin(10) * t86 + t13;
t2 = qJD(6) * t15 + t215 * t6 + t219 * t8;
t443 = t2 * mrSges(7,2);
t350 = t219 * t36;
t16 = t215 * t30 + t350;
t3 = -qJD(6) * t16 - t215 * t8 + t219 * t6;
t442 = t3 * mrSges(7,1);
t441 = t13 * mrSges(6,2);
t440 = t14 * mrSges(6,1);
t249 = t215 * t216 - t219 * t220;
t410 = qJD(5) + qJD(6);
t439 = t410 * t249;
t254 = t13 * t220 - t14 * t216;
t437 = -t51 * t318 - t52 * t319 + t254;
t436 = t457 * pkin(8) + pkin(2) * (m(4) + t457) - t445;
t268 = mrSges(5,1) * t217 + mrSges(5,2) * t221;
t431 = mrSges(3,2) - mrSges(4,3);
t434 = t217 * t455 - t268 + t431;
t433 = t16 * mrSges(7,2) + Ifges(5,6) * qJD(4) / 0.2e1 + qJD(2) * t450 + Ifges(6,5) * t452 + Ifges(6,6) * t453 + Ifges(6,3) * t451 - t15 * mrSges(7,1) + t454;
t25 = qJD(6) * t272 + t215 * t86 + t219 * t85;
t397 = t25 / 0.2e1;
t26 = -qJD(6) * t94 - t215 * t85 + t219 * t86;
t396 = t26 / 0.2e1;
t389 = t85 / 0.2e1;
t388 = t86 / 0.2e1;
t163 = qJDD(6) + t167;
t382 = t163 / 0.2e1;
t381 = t167 / 0.2e1;
t190 = t223 * t216;
t191 = t223 * t220;
t117 = t190 * t219 + t191 * t215;
t430 = qJD(6) * t117 + t215 * t446 + t219 * t447;
t118 = t190 * t215 - t191 * t219;
t429 = -qJD(6) * t118 - t215 * t447 + t219 * t446;
t127 = t216 * t182 + t195;
t336 = t216 * t221;
t101 = -pkin(10) * t336 + t127;
t166 = t220 * t182;
t332 = t220 * t221;
t90 = -pkin(10) * t332 + t217 * t278 + t166;
t42 = t101 * t219 + t215 * t90;
t428 = -qJD(6) * t42 + t215 * t448 + t219 * t449;
t41 = -t101 * t215 + t219 * t90;
t427 = qJD(6) * t41 + t215 * t449 - t219 * t448;
t398 = m(7) * pkin(5);
t425 = -mrSges(6,1) - t398;
t424 = -qJD(5) * t127 - t216 * t287 + t456;
t34 = -mrSges(6,1) * t86 + mrSges(6,2) * t85;
t422 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t180 - t34;
t174 = t215 * t220 + t216 * t219;
t135 = t174 * t221;
t421 = t249 * qJD(2) - qJD(4) * t135 + t217 * t439;
t137 = t249 * t221;
t158 = t174 * qJD(2);
t99 = t410 * t174;
t420 = -qJD(4) * t137 - t217 * t99 - t158;
t419 = -t293 - (-qJD(2) * t373 + t154) * t217 + pkin(5) * t319;
t418 = t217 * t249;
t302 = mrSges(5,3) * t326;
t417 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t171 + mrSges(6,2) * t172 + t302;
t416 = m(6) * pkin(9) + t444;
t48 = -t154 * t323 - t217 * t314 + t221 * t438;
t414 = t217 * t47 + t221 * t48;
t361 = Ifges(5,4) * t217;
t264 = t221 * Ifges(5,1) - t361;
t164 = Ifges(6,4) * t171;
t81 = t172 * Ifges(6,1) + t203 * Ifges(6,5) + t164;
t413 = Ifges(5,5) * qJD(4) + qJD(2) * t264 + t220 * t81;
t69 = mrSges(6,1) * t167 - mrSges(6,3) * t85;
t70 = -mrSges(6,2) * t167 + mrSges(6,3) * t86;
t412 = -t216 * t69 + t220 * t70;
t408 = mrSges(5,1) - t455;
t407 = mrSges(5,2) - t416;
t96 = -qJD(4) * pkin(4) - t113;
t405 = -m(6) * t96 - t417;
t403 = -m(6) * t277 + t434 + t444 * t221 + (-m(7) + t383) * qJ(3);
t225 = qJD(2) ^ 2;
t402 = m(5) / 0.2e1;
t401 = m(6) / 0.2e1;
t400 = Ifges(7,4) * t397 + Ifges(7,2) * t396 + Ifges(7,6) * t382;
t399 = Ifges(7,1) * t397 + Ifges(7,4) * t396 + Ifges(7,5) * t382;
t395 = Ifges(6,1) * t389 + Ifges(6,4) * t388 + Ifges(6,5) * t381;
t375 = Ifges(7,4) * t94;
t32 = Ifges(7,2) * t272 + Ifges(7,6) * t196 + t375;
t394 = -t32 / 0.2e1;
t393 = t32 / 0.2e1;
t87 = Ifges(7,4) * t272;
t33 = Ifges(7,1) * t94 + Ifges(7,5) * t196 + t87;
t392 = -t33 / 0.2e1;
t391 = t33 / 0.2e1;
t386 = t272 / 0.2e1;
t384 = t94 / 0.2e1;
t379 = t172 / 0.2e1;
t377 = t196 / 0.2e1;
t374 = pkin(5) * t172;
t372 = t15 * mrSges(7,3);
t371 = t16 * mrSges(7,3);
t345 = cos(pkin(11));
t275 = t345 * t218;
t212 = sin(pkin(11));
t340 = t212 * t222;
t149 = t214 * t340 + t275;
t342 = t212 * t213;
t103 = t149 * t217 + t221 * t342;
t274 = t345 * t222;
t341 = t212 * t218;
t150 = -t214 * t341 + t274;
t366 = (-t103 * t209 + t150 * t210) * mrSges(7,1) + (-t103 * t210 - t150 * t209) * mrSges(7,2);
t147 = -t214 * t274 + t341;
t276 = t213 * t345;
t105 = -t147 * t217 + t221 * t276;
t148 = t214 * t275 + t340;
t365 = (t105 * t209 + t148 * t210) * mrSges(7,1) + (t105 * t210 - t148 * t209) * mrSges(7,2);
t299 = t217 * t338;
t153 = t214 * t221 - t299;
t339 = t213 * t218;
t364 = (-t153 * t209 + t210 * t339) * mrSges(7,1) + (-t153 * t210 - t209 * t339) * mrSges(7,2);
t363 = mrSges(6,3) * t171;
t362 = mrSges(6,3) * t172;
t359 = Ifges(6,4) * t216;
t358 = Ifges(6,4) * t220;
t354 = t172 * Ifges(6,4);
t44 = -qJDD(4) * pkin(4) - t48;
t347 = t221 * t44;
t179 = qJD(2) * qJ(3) + t295;
t343 = t179 * t222;
t330 = pkin(2) * t338 + qJ(3) * t339;
t325 = qJD(2) * t222;
t316 = qJDD(2) * mrSges(4,2);
t7 = -mrSges(7,1) * t26 + mrSges(7,2) * t25;
t311 = t7 - t422;
t309 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t163;
t308 = Ifges(6,5) * t85 + Ifges(6,6) * t86 + Ifges(6,3) * t167;
t45 = -mrSges(7,1) * t272 + mrSges(7,2) * t94;
t306 = t45 + t417;
t305 = -t383 + t432;
t303 = mrSges(5,3) * t208;
t291 = t213 * t325;
t288 = t217 * t320;
t285 = qJD(1) * t325;
t273 = -t315 / 0.2e1;
t270 = t213 * t285;
t263 = Ifges(6,1) * t220 - t359;
t262 = Ifges(6,1) * t216 + t358;
t260 = -Ifges(6,2) * t216 + t358;
t259 = Ifges(6,2) * t220 + t359;
t258 = -Ifges(5,5) * t217 - Ifges(5,6) * t221;
t257 = Ifges(6,5) * t220 - Ifges(6,6) * t216;
t256 = Ifges(6,5) * t216 + Ifges(6,6) * t220;
t253 = t216 * t52 + t220 * t51;
t108 = -t153 * t216 + t213 * t333;
t109 = t153 * t220 + t213 * t337;
t49 = t108 * t219 - t109 * t215;
t50 = t108 * t215 + t109 * t219;
t123 = -mrSges(6,2) * t203 + t363;
t124 = mrSges(6,1) * t203 - t362;
t250 = -t216 * t123 - t220 * t124;
t248 = qJ(3) * t119 + qJD(3) * t179;
t247 = t309 + t442 - t443;
t152 = t214 * t217 + t221 * t338;
t245 = t119 * t218 + t179 * t325;
t244 = t179 * (mrSges(5,1) * t221 - mrSges(5,2) * t217);
t243 = t217 * (-Ifges(5,2) * t221 - t361);
t242 = t221 * (-Ifges(5,1) * t217 - t360);
t236 = t216 * t317 + t217 * t322;
t232 = Ifges(6,5) * t221 - t217 * t263;
t231 = Ifges(6,6) * t221 - t217 * t260;
t230 = Ifges(6,3) * t221 - t217 * t257;
t186 = -qJD(4) * mrSges(5,2) - t303;
t175 = t268 * qJD(2);
t170 = -qJD(2) * pkin(2) + t255;
t169 = (-t224 + t373) * t221;
t143 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t181;
t141 = (t285 + t313) * t213;
t134 = t174 * t217;
t130 = qJD(2) * t418;
t129 = t217 * t158;
t128 = -qJDD(2) * pkin(2) + t246;
t126 = t166 - t298;
t120 = -pkin(5) * t235 + t288;
t111 = -mrSges(5,1) * t181 + mrSges(5,2) * t180;
t107 = -qJD(4) * t299 + t214 * t321 - t221 * t292;
t106 = -qJD(4) * t152 + t217 * t292;
t80 = t171 * Ifges(6,2) + t203 * Ifges(6,6) + t354;
t78 = -pkin(5) * t171 + t96;
t77 = mrSges(7,1) * t196 - mrSges(7,3) * t94;
t76 = -mrSges(7,2) * t196 + mrSges(7,3) * t272;
t56 = t174 * t323 + t221 * t439;
t54 = qJD(4) * t418 - t221 * t99;
t38 = qJD(5) * t108 + t106 * t220 + t216 * t291;
t37 = -qJD(5) * t109 - t106 * t216 + t220 * t291;
t28 = t85 * Ifges(6,4) + t86 * Ifges(6,2) + t167 * Ifges(6,6);
t27 = -pkin(5) * t86 + t44;
t20 = -mrSges(7,2) * t163 + mrSges(7,3) * t26;
t19 = mrSges(7,1) * t163 - mrSges(7,3) * t25;
t18 = t219 * t35 - t353;
t17 = -t215 * t35 - t350;
t10 = -qJD(6) * t50 - t215 * t38 + t219 * t37;
t9 = qJD(6) * t49 + t215 * t37 + t219 * t38;
t1 = [t10 * t77 + t106 * t186 + t108 * t69 + t109 * t70 + t38 * t123 + t37 * t124 + t153 * t143 + t49 * t19 + t50 * t20 + t9 * t76 + t311 * t152 + t306 * t107 + (-m(2) - m(3) - t305) * g(3) + m(5) * (t106 * t114 - t107 * t113 - t152 * t48 + t153 * t47) + m(7) * (t10 * t15 + t107 * t78 + t152 * t27 + t16 * t9 + t2 * t50 + t3 * t49) + m(6) * (t107 * t96 + t108 * t14 + t109 * t13 + t152 * t44 + t37 * t51 + t38 * t52) + (t175 * t325 + t218 * t111 + m(5) * t245 + m(3) * (t140 * t222 + t141 * t218) + m(4) * (-t128 * t222 + t170 * t327 + t245) + (-t218 * t368 - t222 * t431) * t225 + (-t218 * t431 + t222 * t368) * qJDD(2)) * t213 + (m(2) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t214 ^ 2) * qJDD(1); (t113 * t323 - t414) * mrSges(5,3) + (-t135 * t2 + t137 * t3 - t15 * t54 + t16 * t56) * mrSges(7,3) - (t216 * t81 + t220 * t80) * t317 / 0.2e1 + (t309 + t308) * t217 / 0.2e1 + t27 * (mrSges(7,1) * t135 - mrSges(7,2) * t137) + (-Ifges(7,5) * t137 - Ifges(7,6) * t135 + Ifges(7,3) * t217) * t382 + (-Ifges(7,4) * t137 - Ifges(7,2) * t135 + Ifges(7,6) * t217) * t396 + (-Ifges(7,1) * t137 - Ifges(7,4) * t135 + Ifges(7,5) * t217) * t397 + t427 * t76 + (t120 * t78 + t15 * t428 + t16 * t427 + t169 * t27 + t2 * t42 + t3 * t41) * m(7) + t428 * t77 + t423 * t123 + (t126 * t14 + t127 * t13 + (t323 * t96 - t347) * t224 + t423 * t52 + t424 * t51) * m(6) + t424 * t124 - t413 * t323 / 0.2e1 + (((-t113 * t217 + t114 * t221) * qJD(4) + t414) * t224 + t248 - (t343 + (t113 * t221 + t114 * t217) * t218) * t329) * m(5) + t417 * t288 + t255 * t175 + (-(t170 * t218 + t343) * t329 - pkin(2) * t128 + t248) * m(4) + (Ifges(7,4) * t54 + Ifges(7,2) * t56) * t386 + (-t217 * t295 + t287) * t186 + (Ifges(7,1) * t54 + Ifges(7,4) * t56) * t384 + (Ifges(5,4) * t180 + Ifges(5,2) * t181) * t458 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + (Ifges(7,5) * t54 + Ifges(7,6) * t56) * t377 + t217 * t440 + t217 * t442 + (-t13 * t336 - t14 * t332 + t235 * t52 + t236 * t51) * mrSges(6,3) + t181 * t450 + t143 * t334 + qJD(4) * t244 + t96 * (-mrSges(6,1) * t235 - mrSges(6,2) * t236) + t422 * t221 * t224 + t266 * t347 + (m(7) * t78 - t405 + t45) * t221 * t295 + (qJD(2) * qJD(3) + t119 - t270 + t312) * mrSges(4,3) + t221 * (Ifges(5,1) * t180 + Ifges(5,4) * t181) / 0.2e1 + t80 * t289 / 0.2e1 + qJDD(4) * (Ifges(5,5) * t221 - Ifges(5,6) * t217) + (-m(4) * t330 - t457 * (pkin(8) * t338 + t330) + (t445 * t222 + (t221 * t416 + t434) * t218) * t213) * g(3) - t217 * t443 - t217 * t441 + t169 * t7 + (t128 - t192) * mrSges(4,2) + (t149 * t436 + t150 * t403) * g(1) + (t147 * t436 + t148 * t403) * g(2) + t126 * t69 + t127 * t70 + t120 * t45 + qJ(3) * t111 + (-t141 + t270) * mrSges(3,2) + t78 * (-mrSges(7,1) * t56 + mrSges(7,2) * t54) + (t51 * mrSges(6,1) - t52 * mrSges(6,2) - t114 * mrSges(5,3) + Ifges(7,5) * t384 + Ifges(7,6) * t386 + Ifges(7,3) * t377 - t433) * t321 + t41 * t19 + t42 * t20 + (t140 + t192) * mrSges(3,1) + t242 * t315 / 0.2e1 - pkin(2) * t316 + t171 * (qJD(4) * t231 - t259 * t317) / 0.2e1 + t203 * (qJD(4) * t230 - t256 * t317) / 0.2e1 + (qJD(4) * t232 - t262 * t317) * t379 + (Ifges(6,3) * t217 + t221 * t257) * t381 + (Ifges(6,6) * t217 + t221 * t260) * t388 + (Ifges(6,5) * t217 + t221 * t263) * t389 + t54 * t391 + t56 * t393 + t332 * t395 - t137 * t399 - t135 * t400 + t243 * t273 - t28 * t336 / 0.2e1 + qJD(4) ^ 2 * t258 / 0.2e1 + t180 * t264 / 0.2e1 + t119 * t268; t316 - t225 * mrSges(4,3) - t134 * t19 - t418 * t20 + t421 * t77 + t420 * t76 + ((t123 * t220 - t124 * t216 + t186) * qJD(4) - t311) * t221 + (qJD(4) * t306 + qJD(5) * t250 + t143 + t412) * t217 + m(4) * t128 + 0.2e1 * ((qJD(4) * t114 + t48) * t402 + (t322 * t52 - t324 * t51 - t44) * t401) * t221 + 0.2e1 * ((-qJD(4) * t113 + t47) * t402 + (qJD(4) * t96 + t437) * t401) * t217 + (-m(6) * t253 + t179 * t383 - t175 + t250) * qJD(2) + (-g(1) * t149 - g(2) * t147 + g(3) * t338) * t305 + (-t134 * t3 + t15 * t421 + t16 * t420 - t2 * t418 - t27 * t221 + t323 * t78) * m(7); (t171 * t260 + t172 * t263 + t203 * t257) * qJD(5) / 0.2e1 - (t171 * t231 + t172 * t232 + t203 * t230) * qJD(2) / 0.2e1 - (t319 + t290) * t80 / 0.2e1 + t429 * t77 + (t117 * t3 + t118 * t2 + t15 * t429 + t16 * t430 - t205 * t27 + t419 * t78) * m(7) + t430 * t76 + t419 * t45 + (-t124 * t318 - t123 * t319 + m(6) * (-qJD(5) * t253 + t254) + t412) * pkin(9) + t413 * t208 / 0.2e1 + (t152 * t408 + t153 * t407) * g(3) + (t407 * t103 - t408 * (t149 * t221 - t217 * t342)) * g(1) + (-t407 * t105 - t408 * (t147 * t221 + t217 * t276)) * g(2) + (t302 + t405) * t114 - (Ifges(7,1) * t384 + Ifges(7,4) * t386 + Ifges(7,5) * t377 - t372 + t391) * t439 + ((-t130 - t439) * mrSges(7,2) + (t129 + t99) * mrSges(7,1)) * t78 + (-pkin(4) * t44 - t51 * t74 - t52 * t75) * m(6) + (-t242 / 0.2e1 + t243 / 0.2e1) * t225 + (t433 + t454) * t326 + t203 * t96 * t266 + (Ifges(7,1) * t130 + Ifges(7,4) * t129) * t385 - (Ifges(7,4) * t384 + Ifges(7,2) * t386 + Ifges(7,6) * t377 + t371 + t393) * t99 + (Ifges(7,4) * t130 + Ifges(7,2) * t129) * t387 + (-t51 * (mrSges(6,1) * t221 + mrSges(6,3) * t335) - t244 - t52 * (mrSges(6,3) * t216 * t217 - mrSges(6,2) * t221)) * qJD(2) + (-t303 - t186) * t113 + Ifges(5,3) * qJDD(4) + t220 * t28 / 0.2e1 - t205 * t7 + Ifges(5,5) * t180 + Ifges(5,6) * t181 + t437 * mrSges(6,3) - t74 * t124 + t117 * t19 + t118 * t20 - t75 * t123 - t47 * mrSges(5,2) + t48 * mrSges(5,1) - pkin(4) * t34 + (-t129 * t16 + t130 * t15 - t174 * t3 - t2 * t249) * mrSges(7,3) + t27 * (mrSges(7,1) * t249 + mrSges(7,2) * t174) + (Ifges(7,5) * t174 - Ifges(7,6) * t249) * t382 + (Ifges(7,4) * t174 - Ifges(7,2) * t249) * t396 + (Ifges(7,1) * t174 - Ifges(7,4) * t249) * t397 - t249 * t400 + t81 * t318 / 0.2e1 + t256 * t381 + t259 * t388 + t262 * t389 + t130 * t392 + t129 * t394 + t216 * t395 + t174 * t399 + t258 * t273 + t44 * t267 + (Ifges(7,5) * t130 + Ifges(7,6) * t129) * t378; t308 + t247 - t45 * t374 - m(7) * (t15 * t17 + t16 * t18 + t374 * t78) + (Ifges(6,1) * t171 - t354) * t452 + (Ifges(6,5) * t171 - Ifges(6,6) * t172) * t451 - t96 * (mrSges(6,1) * t172 + mrSges(6,2) * t171) - t18 * t76 - t17 * t77 - t441 + t440 + t80 * t379 + (t2 * t215 + t219 * t3 + (-t15 * t215 + t16 * t219) * qJD(6)) * t398 + (-mrSges(7,2) * t78 + Ifges(7,1) * t385 + Ifges(7,4) * t387 + Ifges(7,5) * t378 + t372 + t392) * t272 - (mrSges(7,1) * t78 + Ifges(7,4) * t385 + Ifges(7,2) * t387 + Ifges(7,6) * t378 - t371 + t394) * t94 + (t362 + t124) * t52 + (t363 - t123) * t51 + (-Ifges(6,2) * t172 + t164 + t81) * t453 + (mrSges(6,2) * t109 + t108 * t425 - t364) * g(3) + (-t365 - (t105 * t220 - t148 * t216) * mrSges(6,2) + t425 * (t105 * t216 + t148 * t220)) * g(2) + (-t366 - (-t103 * t220 - t150 * t216) * mrSges(6,2) + t425 * (-t103 * t216 + t150 * t220)) * g(1) + ((-t215 * t77 + t219 * t76) * qJD(6) + t19 * t219 + t20 * t215) * pkin(5); -t78 * (mrSges(7,1) * t94 + mrSges(7,2) * t272) + (Ifges(7,1) * t272 - t375) * t385 + t32 * t384 + (Ifges(7,5) * t272 - Ifges(7,6) * t94) * t378 - t15 * t76 + t16 * t77 - g(1) * t366 - g(2) * t365 - g(3) * t364 + (t15 * t272 + t16 * t94) * mrSges(7,3) + t247 + (-Ifges(7,2) * t94 + t33 + t87) * t387;];
tau  = t1;
