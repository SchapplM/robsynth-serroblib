% Calculate time derivative of joint inertia matrix for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:26
% EndTime: 2019-03-10 02:00:47
% DurationCPUTime: 8.73s
% Computational Cost: add. (13957->802), mult. (36912->1123), div. (0->0), fcn. (35434->10), ass. (0->308)
t300 = sin(qJ(4));
t304 = cos(qJ(4));
t301 = sin(qJ(3));
t354 = qJD(4) * t301;
t305 = cos(qJ(3));
t356 = qJD(3) * t305;
t312 = -t300 * t354 + t304 * t356;
t353 = qJD(4) * t304;
t311 = t300 * t356 + t301 * t353;
t266 = -pkin(3) * t305 - pkin(10) * t301 - pkin(2);
t361 = t304 * t305;
t286 = pkin(9) * t361;
t221 = t300 * t266 + t286;
t404 = qJD(4) * t221;
t298 = cos(pkin(6));
t297 = sin(pkin(6));
t302 = sin(qJ(2));
t365 = t297 * t302;
t235 = t298 * t301 + t305 * t365;
t306 = cos(qJ(2));
t359 = qJD(2) * t297;
t336 = t306 * t359;
t195 = qJD(3) * t235 + t301 * t336;
t234 = -t298 * t305 + t301 * t365;
t196 = -qJD(3) * t234 + t305 * t336;
t364 = t297 * t306;
t313 = -t235 * t304 + t300 * t364;
t358 = qJD(2) * t302;
t337 = t297 * t358;
t100 = qJD(4) * t313 - t196 * t300 + t304 * t337;
t197 = -t235 * t300 - t304 * t364;
t101 = qJD(4) * t197 + t196 * t304 + t300 * t337;
t299 = sin(qJ(5));
t303 = cos(qJ(5));
t114 = t197 * t303 + t299 * t313;
t39 = qJD(5) * t114 + t100 * t299 + t101 * t303;
t115 = t197 * t299 - t303 * t313;
t40 = -qJD(5) * t115 + t100 * t303 - t101 * t299;
t7 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t195;
t8 = Ifges(6,5) * t39 + Ifges(6,6) * t40 + Ifges(6,3) * t195;
t403 = t7 + t8;
t400 = t300 / 0.2e1;
t379 = t304 / 0.2e1;
t399 = pkin(4) * t303;
t398 = -mrSges(6,1) - mrSges(7,1);
t249 = t299 * t304 + t300 * t303;
t394 = qJD(4) + qJD(5);
t186 = t394 * t249;
t314 = t299 * t300 - t303 * t304;
t132 = -t186 * t301 - t314 * t356;
t233 = t314 * t301;
t133 = t233 * t394 - t249 * t356;
t357 = qJD(3) * t301;
t64 = Ifges(7,5) * t132 + Ifges(7,6) * t133 + Ifges(7,3) * t357;
t65 = Ifges(6,5) * t132 + Ifges(6,6) * t133 + Ifges(6,3) * t357;
t397 = -t64 - t65;
t396 = pkin(4) * qJD(5);
t281 = pkin(8) * t365;
t377 = pkin(1) * t306;
t223 = t281 + (-pkin(2) - t377) * t298;
t141 = t234 * pkin(3) - t235 * pkin(10) + t223;
t242 = t298 * t302 * pkin(1) + pkin(8) * t364;
t224 = pkin(9) * t298 + t242;
t225 = (-pkin(2) * t306 - pkin(9) * t302 - pkin(1)) * t297;
t149 = t305 * t224 + t301 * t225;
t143 = -pkin(10) * t364 + t149;
t73 = t300 * t141 + t304 * t143;
t247 = t304 * t266;
t362 = t301 * t304;
t376 = pkin(9) * t300;
t172 = -pkin(11) * t362 + t247 + (-pkin(4) - t376) * t305;
t363 = t300 * t301;
t199 = -pkin(11) * t363 + t221;
t112 = t299 * t172 + t303 * t199;
t387 = -pkin(11) - pkin(10);
t273 = t387 * t300;
t274 = t387 * t304;
t203 = t299 * t273 - t303 * t274;
t241 = t298 * t377 - t281;
t267 = -mrSges(5,1) * t304 + mrSges(5,2) * t300;
t395 = -m(5) * pkin(3) + t267;
t393 = 0.2e1 * m(5);
t392 = 2 * m(6);
t391 = 2 * m(7);
t390 = 0.2e1 * pkin(9);
t389 = -2 * mrSges(3,3);
t386 = t100 / 0.2e1;
t385 = t197 / 0.2e1;
t384 = -t313 / 0.2e1;
t372 = Ifges(5,4) * t300;
t317 = Ifges(5,1) * t304 - t372;
t228 = -Ifges(5,5) * t305 + t301 * t317;
t383 = t228 / 0.2e1;
t371 = Ifges(5,4) * t304;
t271 = Ifges(5,1) * t300 + t371;
t382 = t271 / 0.2e1;
t381 = -t300 / 0.2e1;
t380 = -t304 / 0.2e1;
t378 = m(5) * t305;
t375 = pkin(9) * t305;
t296 = t301 * pkin(9);
t72 = t304 * t141 - t143 * t300;
t51 = pkin(4) * t234 + pkin(11) * t313 + t72;
t61 = pkin(11) * t197 + t73;
t24 = t299 * t51 + t303 * t61;
t374 = Ifges(4,4) * t301;
t373 = Ifges(4,4) * t305;
t370 = Ifges(5,6) * t300;
t230 = t241 * qJD(2);
t369 = t230 * mrSges(3,2);
t231 = t242 * qJD(2);
t368 = t231 * mrSges(3,1);
t367 = t231 * mrSges(4,1);
t366 = t231 * mrSges(4,2);
t185 = t394 * t314;
t118 = -Ifges(7,5) * t185 - Ifges(7,6) * t186;
t119 = -Ifges(6,5) * t185 - Ifges(6,6) * t186;
t264 = (pkin(3) * t301 - pkin(10) * t305) * qJD(3);
t360 = t304 * t264 + t357 * t376;
t265 = pkin(4) * t363 + t296;
t355 = qJD(4) * t300;
t352 = qJD(5) * t299;
t351 = qJD(5) * t303;
t10 = Ifges(6,4) * t39 + Ifges(6,2) * t40 + Ifges(6,6) * t195;
t9 = Ifges(7,4) * t39 + Ifges(7,2) * t40 + Ifges(7,6) * t195;
t349 = t9 / 0.2e1 + t10 / 0.2e1;
t44 = Ifges(5,5) * t101 + Ifges(5,6) * t100 + Ifges(5,3) * t195;
t348 = pkin(4) * t355;
t347 = Ifges(4,6) * t364;
t11 = Ifges(7,1) * t39 + Ifges(7,4) * t40 + Ifges(7,5) * t195;
t12 = Ifges(6,1) * t39 + Ifges(6,4) * t40 + Ifges(6,5) * t195;
t346 = t11 / 0.2e1 + t12 / 0.2e1;
t55 = Ifges(7,4) * t115 + Ifges(7,2) * t114 + Ifges(7,6) * t234;
t56 = Ifges(6,4) * t115 + Ifges(6,2) * t114 + Ifges(6,6) * t234;
t345 = t55 / 0.2e1 + t56 / 0.2e1;
t57 = Ifges(7,1) * t115 + Ifges(7,4) * t114 + Ifges(7,5) * t234;
t58 = Ifges(6,1) * t115 + Ifges(6,4) * t114 + Ifges(6,5) * t234;
t344 = t57 / 0.2e1 + t58 / 0.2e1;
t66 = Ifges(7,4) * t132 + Ifges(7,2) * t133 + Ifges(7,6) * t357;
t67 = Ifges(6,4) * t132 + Ifges(6,2) * t133 + Ifges(6,6) * t357;
t343 = t67 / 0.2e1 + t66 / 0.2e1;
t68 = Ifges(7,1) * t132 + Ifges(7,4) * t133 + Ifges(7,5) * t357;
t69 = Ifges(6,1) * t132 + Ifges(6,4) * t133 + Ifges(6,5) * t357;
t342 = t68 / 0.2e1 + t69 / 0.2e1;
t229 = (pkin(2) * t302 - pkin(9) * t306) * t359;
t82 = -t224 * t357 + t225 * t356 + t301 * t229 + t305 * t230;
t80 = pkin(10) * t337 + t82;
t97 = t195 * pkin(3) - t196 * pkin(10) + t231;
t26 = -qJD(4) * t73 - t300 * t80 + t304 * t97;
t17 = pkin(4) * t195 - pkin(11) * t101 + t26;
t25 = t141 * t353 - t143 * t355 + t300 * t97 + t304 * t80;
t22 = pkin(11) * t100 + t25;
t6 = -qJD(5) * t24 + t303 * t17 - t22 * t299;
t2 = pkin(5) * t195 - qJ(6) * t39 - qJD(6) * t115 + t6;
t30 = mrSges(7,1) * t195 - mrSges(7,3) * t39;
t341 = m(7) * t2 + t30;
t96 = -Ifges(5,1) * t313 + Ifges(5,4) * t197 + Ifges(5,5) * t234;
t340 = t96 * t379;
t339 = Ifges(4,5) * t196 - Ifges(4,6) * t195 + Ifges(4,3) * t337;
t216 = pkin(4) * t311 + pkin(9) * t356;
t289 = -pkin(4) * t304 - pkin(3);
t338 = qJD(4) * t387;
t120 = -Ifges(7,4) * t185 - Ifges(7,2) * t186;
t121 = -Ifges(6,4) * t185 - Ifges(6,2) * t186;
t331 = t120 / 0.2e1 + t121 / 0.2e1;
t122 = -Ifges(7,1) * t185 - Ifges(7,4) * t186;
t123 = -Ifges(6,1) * t185 - Ifges(6,4) * t186;
t330 = t122 / 0.2e1 + t123 / 0.2e1;
t232 = t249 * t301;
t152 = -Ifges(7,4) * t233 - Ifges(7,2) * t232 - Ifges(7,6) * t305;
t153 = -Ifges(6,4) * t233 - Ifges(6,2) * t232 - Ifges(6,6) * t305;
t329 = t152 / 0.2e1 + t153 / 0.2e1;
t154 = -Ifges(7,1) * t233 - Ifges(7,4) * t232 - Ifges(7,5) * t305;
t155 = -Ifges(6,1) * t233 - Ifges(6,4) * t232 - Ifges(6,5) * t305;
t328 = t154 / 0.2e1 + t155 / 0.2e1;
t191 = Ifges(7,4) * t249 - Ifges(7,2) * t314;
t192 = Ifges(6,4) * t249 - Ifges(6,2) * t314;
t327 = t191 / 0.2e1 + t192 / 0.2e1;
t193 = Ifges(7,1) * t249 - Ifges(7,4) * t314;
t194 = Ifges(6,1) * t249 - Ifges(6,4) * t314;
t326 = t193 / 0.2e1 + t194 / 0.2e1;
t107 = mrSges(7,1) * t357 - mrSges(7,3) * t132;
t106 = (pkin(4) * t301 - pkin(11) * t361) * qJD(3) + (-t286 + (pkin(11) * t301 - t266) * t300) * qJD(4) + t360;
t146 = t300 * t264 + t266 * t353 + (-t304 * t357 - t305 * t355) * pkin(9);
t124 = -pkin(11) * t311 + t146;
t43 = -qJD(5) * t112 + t303 * t106 - t124 * t299;
t28 = pkin(5) * t357 - qJ(6) * t132 + qJD(6) * t233 + t43;
t325 = m(7) * t28 + t107;
t13 = -t40 * mrSges(7,1) + t39 * mrSges(7,2);
t324 = (-mrSges(6,2) - mrSges(7,2)) * t303;
t23 = -t299 * t61 + t303 * t51;
t70 = -t133 * mrSges(7,1) + t132 * mrSges(7,2);
t116 = t186 * mrSges(7,1) - t185 * mrSges(7,2);
t111 = t303 * t172 - t199 * t299;
t148 = -t301 * t224 + t225 * t305;
t202 = t303 * t273 + t274 * t299;
t220 = -t300 * t375 + t247;
t323 = -qJD(4) * t220 + t146;
t322 = t337 / 0.2e1;
t262 = t300 * t338;
t263 = t304 * t338;
t137 = -qJD(5) * t203 - t262 * t299 + t303 * t263;
t78 = qJ(6) * t185 - qJD(6) * t249 + t137;
t321 = m(7) * t78 + t185 * mrSges(7,3);
t142 = pkin(3) * t364 - t148;
t293 = Ifges(5,5) * t353;
t320 = -Ifges(5,6) * t355 / 0.2e1 + t293 / 0.2e1 + t118 / 0.2e1 + t119 / 0.2e1;
t319 = Ifges(5,5) * t400 + Ifges(5,6) * t379 - (Ifges(6,6) + Ifges(7,6)) * t314 / 0.2e1 + (Ifges(6,5) + Ifges(7,5)) * t249 / 0.2e1;
t318 = mrSges(5,1) * t300 + mrSges(5,2) * t304;
t316 = -Ifges(5,2) * t300 + t371;
t269 = Ifges(5,2) * t304 + t372;
t315 = t25 * t304 - t26 * t300;
t83 = -t224 * t356 - t225 * t357 + t229 * t305 - t301 * t230;
t5 = t299 * t17 + t303 * t22 + t51 * t351 - t352 * t61;
t42 = t299 * t106 + t303 * t124 + t172 * t351 - t199 * t352;
t136 = t303 * t262 + t299 * t263 + t273 * t351 + t274 * t352;
t91 = -pkin(4) * t197 + t142;
t310 = -t299 * t186 + (t249 * t299 - t303 * t314) * qJD(5);
t81 = -pkin(3) * t337 - t83;
t77 = -qJ(6) * t186 - qJD(6) * t314 + t136;
t309 = t137 * mrSges(6,1) + t78 * mrSges(7,1) - t136 * mrSges(6,2) - t77 * mrSges(7,2) + t118 + t119;
t3 = qJ(6) * t40 + qJD(6) * t114 + t5;
t308 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t403;
t47 = -pkin(4) * t100 + t81;
t29 = qJ(6) * t133 - qJD(6) * t232 + t42;
t307 = t43 * mrSges(6,1) + t28 * mrSges(7,1) - t42 * mrSges(6,2) - t29 * mrSges(7,2) - t397;
t160 = Ifges(5,5) * t312 - Ifges(5,6) * t311 + Ifges(5,3) * t357;
t294 = Ifges(4,5) * t356;
t288 = pkin(5) + t399;
t276 = Ifges(3,5) * t336;
t272 = Ifges(4,1) * t301 + t373;
t270 = Ifges(4,2) * t305 + t374;
t261 = -mrSges(5,1) * t305 - mrSges(5,3) * t362;
t260 = mrSges(5,2) * t305 - mrSges(5,3) * t363;
t259 = (Ifges(4,1) * t305 - t374) * qJD(3);
t258 = t317 * qJD(4);
t257 = (-Ifges(4,2) * t301 + t373) * qJD(3);
t256 = t316 * qJD(4);
t254 = (mrSges(4,1) * t301 + mrSges(4,2) * t305) * qJD(3);
t253 = t318 * qJD(4);
t243 = t318 * t301;
t227 = -Ifges(5,6) * t305 + t301 * t316;
t226 = -Ifges(5,3) * t305 + (Ifges(5,5) * t304 - t370) * t301;
t222 = pkin(5) * t314 + t289;
t212 = -mrSges(5,2) * t357 - mrSges(5,3) * t311;
t211 = mrSges(5,1) * t357 - mrSges(5,3) * t312;
t207 = -mrSges(6,1) * t305 + mrSges(6,3) * t233;
t206 = -mrSges(7,1) * t305 + mrSges(7,3) * t233;
t205 = mrSges(6,2) * t305 - mrSges(6,3) * t232;
t204 = mrSges(7,2) * t305 - mrSges(7,3) * t232;
t201 = -mrSges(4,1) * t364 - t235 * mrSges(4,3);
t200 = mrSges(4,2) * t364 - t234 * mrSges(4,3);
t188 = mrSges(6,1) * t314 + mrSges(6,2) * t249;
t187 = mrSges(7,1) * t314 + mrSges(7,2) * t249;
t178 = pkin(5) * t232 + t265;
t171 = mrSges(5,1) * t311 + mrSges(5,2) * t312;
t167 = pkin(5) * t186 + t348;
t166 = mrSges(6,1) * t232 - mrSges(6,2) * t233;
t165 = mrSges(7,1) * t232 - mrSges(7,2) * t233;
t164 = -qJ(6) * t314 + t203;
t163 = -qJ(6) * t249 + t202;
t162 = -t271 * t354 + (Ifges(5,5) * t301 + t305 * t317) * qJD(3);
t161 = -t269 * t354 + (Ifges(5,6) * t301 + t305 * t316) * qJD(3);
t159 = mrSges(4,1) * t337 - mrSges(4,3) * t196;
t158 = -mrSges(4,2) * t337 - mrSges(4,3) * t195;
t157 = Ifges(4,1) * t235 - Ifges(4,4) * t234 - Ifges(4,5) * t364;
t156 = Ifges(4,4) * t235 - Ifges(4,2) * t234 - t347;
t151 = -Ifges(6,5) * t233 - Ifges(6,6) * t232 - Ifges(6,3) * t305;
t150 = -Ifges(7,5) * t233 - Ifges(7,6) * t232 - Ifges(7,3) * t305;
t147 = t360 - t404;
t145 = mrSges(5,1) * t234 + mrSges(5,3) * t313;
t144 = -mrSges(5,2) * t234 + mrSges(5,3) * t197;
t131 = -mrSges(5,1) * t197 - mrSges(5,2) * t313;
t125 = mrSges(4,1) * t195 + mrSges(4,2) * t196;
t117 = mrSges(6,1) * t186 - mrSges(6,2) * t185;
t110 = -mrSges(6,2) * t357 + mrSges(6,3) * t133;
t109 = -mrSges(7,2) * t357 + mrSges(7,3) * t133;
t108 = mrSges(6,1) * t357 - mrSges(6,3) * t132;
t105 = Ifges(4,1) * t196 - Ifges(4,4) * t195 + Ifges(4,5) * t337;
t104 = Ifges(4,4) * t196 - Ifges(4,2) * t195 + Ifges(4,6) * t337;
t95 = -Ifges(5,4) * t313 + Ifges(5,2) * t197 + Ifges(5,6) * t234;
t94 = -Ifges(5,5) * t313 + Ifges(5,6) * t197 + Ifges(5,3) * t234;
t90 = -pkin(5) * t133 + t216;
t89 = mrSges(6,1) * t234 - mrSges(6,3) * t115;
t88 = mrSges(7,1) * t234 - mrSges(7,3) * t115;
t87 = -mrSges(6,2) * t234 + mrSges(6,3) * t114;
t86 = -mrSges(7,2) * t234 + mrSges(7,3) * t114;
t85 = -qJ(6) * t232 + t112;
t84 = -pkin(5) * t305 + qJ(6) * t233 + t111;
t75 = mrSges(5,1) * t195 - mrSges(5,3) * t101;
t74 = -mrSges(5,2) * t195 + mrSges(5,3) * t100;
t71 = -mrSges(6,1) * t133 + mrSges(6,2) * t132;
t63 = -mrSges(6,1) * t114 + mrSges(6,2) * t115;
t62 = -mrSges(7,1) * t114 + mrSges(7,2) * t115;
t59 = -pkin(5) * t114 + t91;
t54 = Ifges(6,5) * t115 + Ifges(6,6) * t114 + Ifges(6,3) * t234;
t53 = Ifges(7,5) * t115 + Ifges(7,6) * t114 + Ifges(7,3) * t234;
t52 = -mrSges(5,1) * t100 + mrSges(5,2) * t101;
t46 = Ifges(5,1) * t101 + Ifges(5,4) * t100 + Ifges(5,5) * t195;
t45 = Ifges(5,4) * t101 + Ifges(5,2) * t100 + Ifges(5,6) * t195;
t33 = -mrSges(6,2) * t195 + mrSges(6,3) * t40;
t32 = -mrSges(7,2) * t195 + mrSges(7,3) * t40;
t31 = mrSges(6,1) * t195 - mrSges(6,3) * t39;
t20 = qJ(6) * t114 + t24;
t19 = -pkin(5) * t40 + t47;
t18 = pkin(5) * t234 - qJ(6) * t115 + t23;
t14 = -mrSges(6,1) * t40 + mrSges(6,2) * t39;
t1 = [(t18 * t2 + t19 * t59 + t20 * t3) * t391 + (t23 * t6 + t24 * t5 + t47 * t91) * t392 + (t142 * t81 + t25 * t73 + t26 * t72) * t393 + (-t104 + t44 + 0.2e1 * t367 + t403) * t234 + (t55 + t56) * t40 + (t105 + 0.2e1 * t366) * t235 + (t276 - 0.2e1 * t368 - 0.2e1 * t369) * t298 - t313 * t46 + (t57 + t58) * t39 + (-t306 * t339 + 0.2e1 * (t230 * t306 + t231 * t302) * mrSges(3,3) + ((t241 * t389 + Ifges(3,5) * t298 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t306) * t297) * t306 + (t242 * t389 + Ifges(4,5) * t235 - 0.2e1 * Ifges(3,6) * t298 - Ifges(4,6) * t234 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t302 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t306) * t297) * t302) * qJD(2)) * t297 + 0.2e1 * t223 * t125 + 0.2e1 * t83 * t201 + t196 * t157 + t197 * t45 + 0.2e1 * t82 * t200 + 0.2e1 * t149 * t158 + 0.2e1 * t148 * t159 + 0.2e1 * t25 * t144 + 0.2e1 * t26 * t145 + 0.2e1 * t142 * t52 + 0.2e1 * t81 * t131 + t100 * t95 + t101 * t96 + 0.2e1 * t3 * t86 + 0.2e1 * t5 * t87 + 0.2e1 * t2 * t88 + 0.2e1 * t6 * t89 + 0.2e1 * t91 * t14 + 0.2e1 * t72 * t75 + 0.2e1 * t73 * t74 + 0.2e1 * t19 * t62 + 0.2e1 * t47 * t63 + (-t156 + t94 + t53 + t54) * t195 + 0.2e1 * t59 * t13 + 0.2e1 * t18 * t30 + 0.2e1 * t23 * t31 + 0.2e1 * t20 * t32 + 0.2e1 * t24 * t33 + (t12 + t11) * t115 + (t9 + t10) * t114 + 0.2e1 * m(3) * (t230 * t242 - t231 * t241) + 0.2e1 * m(4) * (t148 * t83 + t149 * t82 + t223 * t231); -t368 + m(6) * (t111 * t6 + t112 * t5 + t216 * t91 + t23 * t43 + t24 * t42 + t265 * t47) + m(7) * (t178 * t19 + t18 * t28 + t2 * t84 + t20 * t29 + t3 * t85 + t59 * t90) + t101 * t383 + t162 * t384 + t161 * t385 + t227 * t386 + (-t257 / 0.2e1 + t160 / 0.2e1 + t64 / 0.2e1 + t65 / 0.2e1) * t234 + ((-t148 * mrSges(4,3) + t340 + t95 * t381 + t157 / 0.2e1) * t305 + (t347 / 0.2e1 - t156 / 0.2e1 + t94 / 0.2e1 + t53 / 0.2e1 + t54 / 0.2e1 - t149 * mrSges(4,3)) * t301 + (-t301 * t200 + (t131 - t201) * t305 + m(4) * (-t148 * t305 - t149 * t301) + t142 * t378) * pkin(9)) * qJD(3) + (Ifges(4,5) * t322 + t46 * t379 + t366 + t105 / 0.2e1 + t45 * t381 - t83 * mrSges(4,3) + (t380 * t95 + t381 * t96) * qJD(4) + (-t159 + t52) * pkin(9)) * t301 + m(4) * (-pkin(2) * t231 - t296 * t83 + t375 * t82) + m(5) * (t146 * t73 + t147 * t72 + t220 * t26 + t221 * t25 + t296 * t81) + (Ifges(4,6) * t322 + pkin(9) * t158 + t82 * mrSges(4,3) - t367 + t104 / 0.2e1 - t44 / 0.2e1 - t7 / 0.2e1 - t8 / 0.2e1) * t305 + (-t306 * t294 / 0.2e1 - Ifges(3,6) * t358) * t297 + (-t270 / 0.2e1 + t226 / 0.2e1 + t150 / 0.2e1 + t151 / 0.2e1) * t195 + t342 * t115 + t343 * t114 + t344 * t132 + t345 * t133 - t346 * t233 - t349 * t232 - t369 + t328 * t39 + t329 * t40 + t265 * t14 + t196 * t272 / 0.2e1 + t223 * t254 + t235 * t259 / 0.2e1 + t25 * t260 + t26 * t261 + t81 * t243 + t220 * t75 + t221 * t74 + t3 * t204 + t5 * t205 + t2 * t206 + t6 * t207 + t72 * t211 + t73 * t212 + t216 * t63 + t142 * t171 + t178 * t13 + t19 * t165 + t47 * t166 + t146 * t144 + t147 * t145 - pkin(2) * t125 + t112 * t33 + t18 * t107 + t23 * t108 + t20 * t109 + t24 * t110 + t111 * t31 + t29 * t86 + t42 * t87 + t28 * t88 + t43 * t89 + t90 * t62 + t91 * t71 + t84 * t30 + t85 * t32 + t59 * t70 + t276; (t178 * t90 + t28 * t84 + t29 * t85) * t391 + (t111 * t43 + t112 * t42 + t216 * t265) * t392 + (t146 * t221 + t147 * t220) * t393 - (t68 + t69) * t233 - (t67 + t66) * t232 + (t152 + t153) * t133 + (t154 + t155) * t132 + (-t160 + t257 + (-t300 * t227 + t304 * t228 + t243 * t390 + t272) * qJD(3) + t397) * t305 + (t171 * t390 - t300 * t161 + t304 * t162 + t259 + (-t227 * t304 - t228 * t300) * qJD(4) + (0.2e1 * pkin(9) ^ 2 * t378 + t150 + t151 + t226 - t270) * qJD(3)) * t301 + 0.2e1 * t265 * t71 - 0.2e1 * pkin(2) * t254 + 0.2e1 * t146 * t260 + 0.2e1 * t147 * t261 + 0.2e1 * t220 * t211 + 0.2e1 * t221 * t212 + 0.2e1 * t29 * t204 + 0.2e1 * t42 * t205 + 0.2e1 * t28 * t206 + 0.2e1 * t43 * t207 + 0.2e1 * t216 * t166 + 0.2e1 * t178 * t70 + 0.2e1 * t90 * t165 + 0.2e1 * t84 * t107 + 0.2e1 * t85 * t109 + 0.2e1 * t111 * t108 + 0.2e1 * t112 * t110; t258 * t384 + t256 * t385 + t269 * t386 + t45 * t379 + t101 * t382 + m(7) * (t163 * t2 + t164 * t3 + t167 * t59 + t18 * t78 + t19 * t222 + t20 * t77) + t395 * t81 + t46 * t400 + t339 + (t18 * t185 - t186 * t20 - t2 * t249 - t3 * t314) * mrSges(7,3) + (t185 * t23 - t186 * t24 - t249 * t6 - t314 * t5) * mrSges(6,3) - t349 * t314 + (m(5) * (-t353 * t72 - t355 * t73 + t315) + t304 * t74 - t300 * t75 - t145 * t353 - t144 * t355) * pkin(10) - t344 * t185 - t345 * t186 + t346 * t249 + m(6) * (t136 * t24 + t137 * t23 + t202 * t6 + t203 * t5 + t289 * t47 + t348 * t91) + (t340 + (-t95 / 0.2e1 + pkin(4) * t63) * t300) * qJD(4) + t326 * t39 + t327 * t40 + t330 * t115 + t331 * t114 + t319 * t195 + t320 * t234 + ((-t300 * t73 - t304 * t72) * qJD(4) + t315) * mrSges(5,3) + t289 * t14 + t142 * t253 + t222 * t13 + t202 * t31 + t203 * t33 + t19 * t187 + t47 * t188 + t163 * t30 + t164 * t32 + t167 * t62 + t136 * t87 + t137 * t89 + t91 * t117 + t59 * t116 + t78 * t88 - t82 * mrSges(4,2) + t83 * mrSges(4,1) + t77 * t86 - pkin(3) * t52; ((-mrSges(4,1) + t395) * qJD(3) * pkin(9) - t320) * t305 + (t258 * t379 + t256 * t381 + pkin(9) * t253 + (t269 * t380 + t271 * t381) * qJD(4) + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t319) * qJD(3)) * t301 + (qJD(4) * t383 + t161 / 0.2e1 + t356 * t382 + t323 * mrSges(5,3) + (m(5) * t323 - qJD(4) * t261 + t212) * pkin(10)) * t304 + (t185 * t84 - t186 * t85 - t249 * t28 - t29 * t314) * mrSges(7,3) + (t111 * t185 - t112 * t186 - t249 * t43 - t314 * t42) * mrSges(6,3) - t343 * t314 + t342 * t249 + t326 * t132 + t327 * t133 - t328 * t185 - t329 * t186 - t330 * t233 - t331 * t232 + (t162 / 0.2e1 - t147 * mrSges(5,3) - t269 * t356 / 0.2e1 + (-t227 / 0.2e1 - t221 * mrSges(5,3) + (m(6) * t265 + t166) * pkin(4)) * qJD(4) + (m(5) * (-t147 - t404) - t211 - qJD(4) * t260) * pkin(10)) * t300 + t289 * t71 + t265 * t117 + t222 * t70 + t202 * t108 + t203 * t110 + t77 * t204 + t136 * t205 + t78 * t206 + t137 * t207 + t216 * t188 + t90 * t187 - pkin(3) * t171 + t178 * t116 + t163 * t107 + t164 * t109 + t167 * t165 + t294 + m(6) * (t111 * t137 + t112 * t136 + t202 * t43 + t203 * t42 + t216 * t289) + m(7) * (t163 * t28 + t164 * t29 + t167 * t178 + t222 * t90 + t77 * t85 + t78 * t84); -0.2e1 * pkin(3) * t253 + 0.2e1 * t222 * t116 + 0.2e1 * t289 * t117 + 0.2e1 * t167 * t187 + t304 * t256 + t300 * t258 + (t122 + t123) * t249 - (t120 + t121) * t314 - (t191 + t192) * t186 - (t193 + t194) * t185 + (t304 * t271 + (0.2e1 * pkin(4) * t188 - t269) * t300) * qJD(4) + (t136 * t203 + t137 * t202 + t289 * t348) * t392 + (t163 * t78 + t164 * t77 + t167 * t222) * t391 + 0.2e1 * (t163 * t185 - t164 * t186 - t249 * t78 - t314 * t77) * mrSges(7,3) + 0.2e1 * (-t136 * t314 - t137 * t249 + t185 * t202 - t186 * t203) * mrSges(6,3); t308 - t25 * mrSges(5,2) + t26 * mrSges(5,1) + t341 * t288 + (t303 * t31 + (t32 + t33) * t299 + ((t86 + t87) * t303 + (-t88 - t89) * t299) * qJD(5) + m(7) * (-t18 * t352 + t20 * t351 + t299 * t3) + m(6) * (-t23 * t352 + t24 * t351 + t299 * t5 + t303 * t6)) * pkin(4) + t44; t307 + t160 - t146 * mrSges(5,2) + t147 * mrSges(5,1) + t325 * t288 + (t303 * t108 + (t109 + t110) * t299 + ((t204 + t205) * t303 + (-t206 - t207) * t299) * qJD(5) + m(7) * (t29 * t299 + t351 * t85 - t352 * t84) + m(6) * (-t111 * t352 + t112 * t351 + t299 * t42 + t303 * t43)) * pkin(4); t293 + t321 * t288 + (pkin(10) * t267 - t370) * qJD(4) + (m(7) * (-t163 * t352 + t164 * t351 + t299 * t77) + m(6) * (t136 * t299 + t137 * t303 - t202 * t352 + t203 * t351) + t310 * mrSges(7,3) + (t303 * t185 + t310) * mrSges(6,3)) * pkin(4) + t309; 0.2e1 * (t324 + ((-t288 + t399) * m(7) + t398) * t299) * t396; pkin(5) * t341 + t308; pkin(5) * t325 + t307; pkin(5) * t321 + t309; (t324 + (-m(7) * pkin(5) + t398) * t299) * t396; 0; m(7) * t19 + t13; m(7) * t90 + t70; m(7) * t167 + t116; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
