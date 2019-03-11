% Calculate time derivative of joint inertia matrix for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:20
% EndTime: 2019-03-09 23:51:40
% DurationCPUTime: 9.36s
% Computational Cost: add. (10368->828), mult. (27174->1155), div. (0->0), fcn. (25266->10), ass. (0->327)
t278 = sin(qJ(4));
t401 = (Ifges(6,4) + Ifges(5,5)) * t278;
t279 = sin(qJ(3));
t283 = cos(qJ(3));
t341 = qJD(3) * t283;
t321 = t278 * t341;
t282 = cos(qJ(4));
t338 = qJD(4) * t282;
t288 = t279 * t338 + t321;
t276 = cos(pkin(6));
t275 = sin(pkin(6));
t280 = sin(qJ(2));
t353 = t275 * t280;
t204 = t276 * t279 + t283 * t353;
t284 = cos(qJ(2));
t344 = qJD(2) * t275;
t324 = t284 * t344;
t152 = qJD(3) * t204 + t279 * t324;
t203 = -t276 * t283 + t279 * t353;
t153 = -qJD(3) * t203 + t283 * t324;
t343 = qJD(2) * t280;
t325 = t275 * t343;
t352 = t275 * t284;
t328 = t278 * t352;
t77 = -qJD(4) * t328 + t153 * t278 + t204 * t338 - t282 * t325;
t154 = t204 * t278 + t282 * t352;
t78 = -qJD(4) * t154 + t153 * t282 + t278 * t325;
t26 = Ifges(5,5) * t78 - Ifges(5,6) * t77 + Ifges(5,3) * t152;
t27 = Ifges(6,4) * t78 + Ifges(6,2) * t152 + Ifges(6,6) * t77;
t398 = t26 + t27;
t207 = t276 * t280 * pkin(1) + pkin(8) * t352;
t184 = pkin(9) * t276 + t207;
t185 = (-pkin(2) * t284 - pkin(9) * t280 - pkin(1)) * t275;
t196 = (pkin(2) * t280 - pkin(9) * t284) * t344;
t260 = pkin(8) * t353;
t374 = pkin(1) * t284;
t206 = t276 * t374 - t260;
t197 = t206 * qJD(2);
t342 = qJD(3) * t279;
t58 = -t184 * t341 - t185 * t342 + t283 * t196 - t279 * t197;
t287 = pkin(3) * t325 + t58;
t34 = mrSges(5,1) * t77 + mrSges(5,2) * t78;
t397 = m(5) * t287 - t34;
t277 = sin(qJ(6));
t281 = cos(qJ(6));
t298 = t277 * t282 - t278 * t281;
t199 = t298 * t279;
t396 = -mrSges(7,1) * t277 - mrSges(7,2) * t281;
t395 = qJD(4) - qJD(6);
t355 = qJ(5) * t278;
t385 = pkin(4) + pkin(5);
t394 = -t282 * t385 - t355;
t300 = pkin(4) * t282 + t355;
t336 = qJD(5) * t282;
t393 = qJD(4) * t300 - t336;
t392 = 0.2e1 * m(5);
t391 = 2 * m(6);
t390 = 2 * m(7);
t389 = 0.2e1 * pkin(9);
t388 = -2 * mrSges(3,3);
t155 = t204 * t282 - t328;
t84 = t154 * t281 - t155 * t277;
t387 = t84 / 0.2e1;
t85 = t154 * t277 + t155 * t281;
t386 = t85 / 0.2e1;
t384 = pkin(10) - pkin(11);
t297 = t277 * t278 + t281 * t282;
t146 = t395 * t297;
t382 = t146 / 0.2e1;
t147 = t395 * t298;
t381 = -t147 / 0.2e1;
t150 = -Ifges(7,4) * t298 - Ifges(7,2) * t297;
t380 = t150 / 0.2e1;
t151 = -Ifges(7,1) * t298 - Ifges(7,4) * t297;
t379 = t151 / 0.2e1;
t378 = -t199 / 0.2e1;
t200 = t297 * t279;
t377 = t200 / 0.2e1;
t376 = -t297 / 0.2e1;
t375 = -t298 / 0.2e1;
t373 = pkin(9) * t278;
t372 = pkin(9) * t283;
t96 = t199 * t395 + t297 * t341;
t97 = t146 * t279 - t298 * t341;
t371 = -Ifges(7,5) * t96 - Ifges(7,6) * t97;
t368 = Ifges(4,4) * t279;
t367 = Ifges(4,4) * t283;
t366 = Ifges(5,4) * t278;
t365 = Ifges(5,4) * t282;
t364 = Ifges(6,5) * t278;
t363 = Ifges(6,5) * t282;
t362 = Ifges(5,6) * t282;
t234 = -t277 * qJ(5) - t281 * t385;
t193 = t281 * qJD(5) + qJD(6) * t234;
t361 = t193 * mrSges(7,2);
t235 = t281 * qJ(5) - t277 * t385;
t194 = -t277 * qJD(5) - qJD(6) * t235;
t360 = t194 * mrSges(7,1);
t359 = t197 * mrSges(3,2);
t198 = t207 * qJD(2);
t358 = t198 * mrSges(3,1);
t357 = t198 * mrSges(4,1);
t356 = t198 * mrSges(4,2);
t354 = qJ(5) * t282;
t351 = t278 * t279;
t350 = t279 * t282;
t349 = t282 * t283;
t183 = t260 + (-pkin(2) - t374) * t276;
t107 = pkin(3) * t203 - pkin(10) * t204 + t183;
t119 = t283 * t184 + t279 * t185;
t109 = -pkin(10) * t352 + t119;
t47 = t278 * t107 + t282 * t109;
t87 = Ifges(7,5) * t146 - Ifges(7,6) * t147;
t118 = -t279 * t184 + t283 * t185;
t301 = Ifges(6,3) * t278 + t363;
t186 = -Ifges(6,6) * t283 + t279 * t301;
t302 = -Ifges(5,2) * t278 + t365;
t189 = -Ifges(5,6) * t283 + t279 * t302;
t348 = t186 - t189;
t303 = Ifges(6,1) * t282 + t364;
t190 = -Ifges(6,4) * t283 + t279 * t303;
t304 = Ifges(5,1) * t282 - t366;
t191 = -Ifges(5,5) * t283 + t279 * t304;
t347 = t190 + t191;
t233 = (pkin(3) * t279 - pkin(10) * t283) * qJD(3);
t239 = -pkin(3) * t283 - pkin(10) * t279 - pkin(2);
t346 = t278 * t233 + t239 * t338;
t320 = t282 * t341;
t345 = Ifges(5,5) * t320 + Ifges(5,3) * t342;
t265 = pkin(9) * t349;
t182 = t278 * t239 + t265;
t340 = qJD(4) * t278;
t221 = Ifges(6,4) * t338 + Ifges(6,6) * t340;
t339 = qJD(4) * t279;
t337 = qJD(5) * t278;
t23 = -qJD(6) * t85 - t277 * t78 + t281 * t77;
t24 = qJD(6) * t84 + t277 * t77 + t281 * t78;
t3 = Ifges(7,5) * t24 + Ifges(7,6) * t23 - Ifges(7,3) * t152;
t334 = Ifges(4,6) * t352;
t25 = Ifges(6,5) * t78 + Ifges(6,6) * t152 + Ifges(6,3) * t77;
t28 = Ifges(5,4) * t78 - Ifges(5,2) * t77 + Ifges(5,6) * t152;
t333 = t25 / 0.2e1 - t28 / 0.2e1;
t29 = Ifges(6,1) * t78 + Ifges(6,4) * t152 + Ifges(6,5) * t77;
t30 = Ifges(5,1) * t78 - Ifges(5,4) * t77 + Ifges(5,5) * t152;
t332 = t29 / 0.2e1 + t30 / 0.2e1;
t62 = Ifges(6,5) * t155 + Ifges(6,6) * t203 + Ifges(6,3) * t154;
t65 = Ifges(5,4) * t155 - Ifges(5,2) * t154 + Ifges(5,6) * t203;
t331 = t62 / 0.2e1 - t65 / 0.2e1;
t66 = Ifges(6,1) * t155 + Ifges(6,4) * t203 + Ifges(6,5) * t154;
t67 = Ifges(5,1) * t155 - Ifges(5,4) * t154 + Ifges(5,5) * t203;
t330 = -t67 / 0.2e1 - t66 / 0.2e1;
t251 = t384 * t282;
t329 = Ifges(7,3) * t342;
t40 = t203 * qJ(5) + t47;
t327 = Ifges(4,5) * t153 - Ifges(4,6) * t152 + Ifges(4,3) * t325;
t108 = pkin(3) * t352 - t118;
t326 = -pkin(4) - t373;
t323 = t278 * t339;
t242 = -Ifges(6,3) * t282 + t364;
t127 = -t242 * t339 + (Ifges(6,6) * t279 + t283 * t301) * qJD(3);
t245 = Ifges(5,2) * t282 + t366;
t130 = -t245 * t339 + (Ifges(5,6) * t279 + t283 * t302) * qJD(3);
t319 = t127 / 0.2e1 - t130 / 0.2e1;
t247 = Ifges(6,1) * t278 - t363;
t131 = -t247 * t339 + (Ifges(6,4) * t279 + t283 * t303) * qJD(3);
t248 = Ifges(5,1) * t278 + t365;
t132 = -t248 * t339 + (Ifges(5,5) * t279 + t283 * t304) * qJD(3);
t318 = t132 / 0.2e1 + t131 / 0.2e1;
t317 = t186 / 0.2e1 - t189 / 0.2e1;
t316 = t190 / 0.2e1 + t191 / 0.2e1;
t219 = t301 * qJD(4);
t222 = t302 * qJD(4);
t315 = t219 / 0.2e1 - t222 / 0.2e1;
t224 = t303 * qJD(4);
t225 = t304 * qJD(4);
t314 = t224 / 0.2e1 + t225 / 0.2e1;
t313 = t242 / 0.2e1 - t245 / 0.2e1;
t312 = t247 / 0.2e1 + t248 / 0.2e1;
t51 = -t152 * mrSges(6,1) + t78 * mrSges(6,2);
t46 = t107 * t282 - t278 * t109;
t264 = t278 * t372;
t181 = t239 * t282 - t264;
t311 = Ifges(6,4) * t320 + Ifges(6,2) * t342 + Ifges(6,6) * t288;
t310 = t325 / 0.2e1;
t220 = Ifges(5,5) * t338 - Ifges(5,6) * t340;
t309 = t220 / 0.2e1 + t221 / 0.2e1 - t87 / 0.2e1;
t161 = -qJ(5) * t283 + t182;
t308 = qJD(4) * t265 - t233 * t282 + t239 * t340;
t307 = t362 / 0.2e1 - Ifges(6,6) * t282 / 0.2e1 + Ifges(7,5) * t298 / 0.2e1 + Ifges(7,6) * t297 / 0.2e1 + t401 / 0.2e1;
t241 = -t282 * mrSges(5,1) + t278 * mrSges(5,2);
t306 = mrSges(5,1) * t278 + mrSges(5,2) * t282;
t240 = -t282 * mrSges(6,1) - t278 * mrSges(6,3);
t305 = mrSges(6,1) * t278 - mrSges(6,3) * t282;
t299 = pkin(4) * t278 - t354;
t31 = -pkin(11) * t155 - t203 * t385 - t46;
t32 = pkin(11) * t154 + t40;
t11 = -t277 * t32 + t281 * t31;
t12 = t277 * t31 + t281 * t32;
t274 = t283 * pkin(4);
t134 = pkin(5) * t283 + t264 + t274 + (-pkin(11) * t279 - t239) * t282;
t137 = pkin(11) * t351 + t161;
t74 = t134 * t281 - t137 * t277;
t75 = t134 * t277 + t137 * t281;
t250 = t384 * t278;
t159 = t250 * t281 - t251 * t277;
t160 = t250 * t277 + t251 * t281;
t57 = -t184 * t342 + t185 * t341 + t279 * t196 + t283 * t197;
t55 = pkin(10) * t325 + t57;
t68 = pkin(3) * t152 - pkin(10) * t153 + t198;
t15 = -t107 * t340 - t109 * t338 - t278 * t55 + t282 * t68;
t296 = pkin(9) + t299;
t295 = qJ(5) * t155 - t108;
t294 = -t278 * t385 + t354;
t76 = pkin(11) * t323 + (-pkin(11) * t349 + (-pkin(5) + t326) * t279) * qJD(3) + t308;
t267 = qJ(5) * t342;
t79 = t267 + (-pkin(9) * qJD(3) + pkin(11) * qJD(4)) * t350 + (-qJD(5) + (-pkin(9) * qJD(4) + pkin(11) * qJD(3)) * t278) * t283 + t346;
t19 = qJD(6) * t74 + t277 * t76 + t281 * t79;
t20 = -qJD(6) * t75 - t277 * t79 + t281 * t76;
t293 = -t20 * mrSges(7,1) + t19 * mrSges(7,2) + t371;
t14 = t107 * t338 - t109 * t340 + t278 * t68 + t282 * t55;
t292 = -pkin(9) + t294;
t231 = t384 * t340;
t232 = qJD(4) * t251;
t100 = qJD(6) * t159 - t231 * t281 + t232 * t277;
t101 = -qJD(6) * t160 + t231 * t277 + t232 * t281;
t291 = t101 * mrSges(7,1) - t100 * mrSges(7,2) + t87;
t7 = -pkin(11) * t78 - t152 * t385 - t15;
t10 = t152 * qJ(5) + t203 * qJD(5) + t14;
t8 = pkin(11) * t77 + t10;
t1 = qJD(6) * t11 + t277 * t7 + t281 * t8;
t2 = -qJD(6) * t12 - t277 * t8 + t281 * t7;
t290 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + t3;
t289 = t320 - t323;
t169 = mrSges(6,2) * t320 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t340) * t279;
t116 = (-t282 * t342 - t283 * t340) * pkin(9) + t346;
t286 = qJ(5) * t78 + qJD(5) * t155 + t287;
t272 = Ifges(4,5) * t341;
t253 = Ifges(3,5) * t324;
t249 = Ifges(4,1) * t279 + t367;
t246 = Ifges(4,2) * t283 + t368;
t236 = -pkin(3) - t300;
t230 = -mrSges(6,2) * t351 - mrSges(6,3) * t283;
t229 = mrSges(6,1) * t283 + mrSges(6,2) * t350;
t228 = -mrSges(5,1) * t283 - mrSges(5,3) * t350;
t227 = mrSges(5,2) * t283 - mrSges(5,3) * t351;
t226 = (Ifges(4,1) * t283 - t368) * qJD(3);
t223 = (-Ifges(4,2) * t279 + t367) * qJD(3);
t218 = (mrSges(4,1) * t279 + mrSges(4,2) * t283) * qJD(3);
t217 = t306 * qJD(4);
t216 = t305 * qJD(4);
t212 = pkin(3) - t394;
t209 = t306 * t279;
t208 = t305 * t279;
t202 = qJD(4) * t299 - t337;
t192 = t296 * t279;
t188 = -Ifges(6,2) * t283 + (Ifges(6,4) * t282 + Ifges(6,6) * t278) * t279;
t187 = -Ifges(5,3) * t283 + (Ifges(5,5) * t282 - Ifges(5,6) * t278) * t279;
t176 = qJD(4) * t294 + t337;
t171 = -mrSges(6,2) * t288 + mrSges(6,3) * t342;
t170 = -mrSges(5,2) * t342 - mrSges(5,3) * t288;
t168 = mrSges(5,1) * t342 - mrSges(5,3) * t289;
t164 = mrSges(7,1) * t283 - mrSges(7,3) * t200;
t163 = -mrSges(7,2) * t283 - mrSges(7,3) * t199;
t162 = -t181 + t274;
t158 = -mrSges(4,1) * t352 - mrSges(4,3) * t204;
t157 = mrSges(4,2) * t352 - mrSges(4,3) * t203;
t156 = t292 * t279;
t148 = mrSges(7,1) * t297 - mrSges(7,2) * t298;
t136 = mrSges(5,1) * t288 + mrSges(5,2) * t289;
t135 = mrSges(6,1) * t288 - mrSges(6,3) * t289;
t133 = mrSges(7,1) * t199 + mrSges(7,2) * t200;
t129 = -Ifges(6,4) * t323 + t311;
t128 = -Ifges(5,5) * t323 - Ifges(5,6) * t288 + t345;
t126 = mrSges(4,1) * t325 - mrSges(4,3) * t153;
t125 = -mrSges(4,2) * t325 - mrSges(4,3) * t152;
t124 = Ifges(4,1) * t204 - Ifges(4,4) * t203 - Ifges(4,5) * t352;
t123 = Ifges(4,4) * t204 - Ifges(4,2) * t203 - t334;
t122 = Ifges(7,1) * t200 - Ifges(7,4) * t199 + Ifges(7,5) * t283;
t121 = Ifges(7,4) * t200 - Ifges(7,2) * t199 + Ifges(7,6) * t283;
t120 = Ifges(7,5) * t200 - Ifges(7,6) * t199 + Ifges(7,3) * t283;
t117 = t342 * t373 - t308;
t115 = t279 * t393 + t296 * t341;
t114 = -mrSges(6,1) * t203 + mrSges(6,2) * t155;
t113 = mrSges(5,1) * t203 - mrSges(5,3) * t155;
t112 = -mrSges(5,2) * t203 - mrSges(5,3) * t154;
t111 = -mrSges(6,2) * t154 + mrSges(6,3) * t203;
t110 = t326 * t342 + t308;
t104 = -qJD(5) * t283 + t116 + t267;
t95 = mrSges(5,1) * t154 + mrSges(5,2) * t155;
t94 = mrSges(6,1) * t154 - mrSges(6,3) * t155;
t91 = mrSges(4,1) * t152 + mrSges(4,2) * t153;
t90 = (qJD(4) * t394 + t336) * t279 + t292 * t341;
t89 = Ifges(7,1) * t146 - Ifges(7,4) * t147;
t88 = Ifges(7,4) * t146 - Ifges(7,2) * t147;
t86 = mrSges(7,1) * t147 + mrSges(7,2) * t146;
t83 = mrSges(7,2) * t342 + mrSges(7,3) * t97;
t82 = -mrSges(7,1) * t342 - mrSges(7,3) * t96;
t81 = Ifges(4,1) * t153 - Ifges(4,4) * t152 + Ifges(4,5) * t325;
t80 = Ifges(4,4) * t153 - Ifges(4,2) * t152 + Ifges(4,6) * t325;
t64 = Ifges(6,4) * t155 + Ifges(6,2) * t203 + Ifges(6,6) * t154;
t63 = Ifges(5,5) * t155 - Ifges(5,6) * t154 + Ifges(5,3) * t203;
t60 = -mrSges(7,1) * t203 - mrSges(7,3) * t85;
t59 = mrSges(7,2) * t203 + mrSges(7,3) * t84;
t52 = -mrSges(6,2) * t77 + mrSges(6,3) * t152;
t50 = mrSges(5,1) * t152 - mrSges(5,3) * t78;
t49 = -mrSges(5,2) * t152 - mrSges(5,3) * t77;
t48 = pkin(4) * t154 - t295;
t45 = -mrSges(7,1) * t97 + mrSges(7,2) * t96;
t44 = Ifges(7,1) * t96 + Ifges(7,4) * t97 - Ifges(7,5) * t342;
t43 = Ifges(7,4) * t96 + Ifges(7,2) * t97 - Ifges(7,6) * t342;
t42 = -t329 - t371;
t41 = -pkin(4) * t203 - t46;
t39 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t38 = -t154 * t385 + t295;
t37 = Ifges(7,1) * t85 + Ifges(7,4) * t84 - Ifges(7,5) * t203;
t36 = Ifges(7,4) * t85 + Ifges(7,2) * t84 - Ifges(7,6) * t203;
t35 = Ifges(7,5) * t85 + Ifges(7,6) * t84 - Ifges(7,3) * t203;
t33 = mrSges(6,1) * t77 - mrSges(6,3) * t78;
t18 = -mrSges(7,1) * t152 - mrSges(7,3) * t24;
t17 = mrSges(7,2) * t152 + mrSges(7,3) * t23;
t16 = pkin(4) * t77 - t286;
t13 = -pkin(4) * t152 - t15;
t9 = -t385 * t77 + t286;
t6 = -mrSges(7,1) * t23 + mrSges(7,2) * t24;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t23 - Ifges(7,5) * t152;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t23 - Ifges(7,6) * t152;
t21 = [(t66 + t67) * t78 + 0.2e1 * m(3) * (t197 * t207 - t198 * t206) + 0.2e1 * m(4) * (t118 * t58 + t119 * t57 + t183 * t198) + (t62 - t65) * t77 + (t1 * t12 + t11 * t2 + t38 * t9) * t390 + (t10 * t40 + t13 * t41 + t16 * t48) * t391 + 0.2e1 * t183 * t91 + 0.2e1 * t58 * t158 + 0.2e1 * t57 * t157 + t153 * t124 + 0.2e1 * t119 * t125 + 0.2e1 * t118 * t126 + 0.2e1 * t108 * t34 + 0.2e1 * t10 * t111 + 0.2e1 * t14 * t112 + 0.2e1 * t15 * t113 + 0.2e1 * t13 * t114 + 0.2e1 * t16 * t94 + t84 * t4 + t85 * t5 + (t29 + t30) * t155 + 0.2e1 * t1 * t59 + 0.2e1 * t2 * t60 + 0.2e1 * t47 * t49 + 0.2e1 * t46 * t50 + 0.2e1 * t41 * t51 + 0.2e1 * t40 * t52 + 0.2e1 * t48 * t33 + t23 * t36 + t24 * t37 + 0.2e1 * t38 * t6 + 0.2e1 * t9 * t39 + (t25 - t28) * t154 + 0.2e1 * t12 * t17 + 0.2e1 * t11 * t18 + (t63 + t64 - t123 - t35) * t152 + (-t108 * t287 + t14 * t47 + t15 * t46) * t392 - 0.2e1 * t287 * t95 + (-t284 * t327 + 0.2e1 * (t197 * t284 + t198 * t280) * mrSges(3,3) + ((t206 * t388 + Ifges(3,5) * t276 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t284) * t275) * t284 + (t207 * t388 + Ifges(4,5) * t204 - 0.2e1 * Ifges(3,6) * t276 - Ifges(4,6) * t203 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t280 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t284) * t275) * t280) * qJD(2)) * t275 + (t81 + 0.2e1 * t356) * t204 + (t253 - 0.2e1 * t358 - 0.2e1 * t359) * t276 + (-t3 - t80 + 0.2e1 * t357 + t398) * t203; (-m(4) * t198 - t91) * pkin(2) + t43 * t387 + m(6) * (t10 * t161 + t104 * t40 + t110 * t41 + t115 * t48 + t13 * t162 + t16 * t192) + m(7) * (t1 * t75 + t11 * t20 + t12 * t19 + t156 * t9 + t2 * t74 + t38 * t90) + t5 * t377 + t4 * t378 + t44 * t386 + t153 * t249 / 0.2e1 + t183 * t218 + t204 * t226 / 0.2e1 + t14 * t227 + t15 * t228 + t13 * t229 + t10 * t230 + t16 * t208 + t181 * t50 + t182 * t49 + t192 * t33 + t161 * t52 + t162 * t51 + t1 * t163 + t2 * t164 + t46 * t168 + t41 * t169 + t47 * t170 + t40 * t171 + t156 * t6 + t9 * t133 + t48 * t135 + t108 * t136 + t117 * t113 + t23 * t121 / 0.2e1 + t24 * t122 / 0.2e1 + t104 * t111 + t110 * t114 + t115 * t94 + t116 * t112 + t96 * t37 / 0.2e1 + t97 * t36 / 0.2e1 + t90 * t39 + t11 * t82 + t12 * t83 + t74 * t18 + t75 * t17 + t19 * t59 + t20 * t60 + t38 * t45 - t287 * t209 + (-t223 / 0.2e1 + t128 / 0.2e1 + t129 / 0.2e1 - t42 / 0.2e1) * t203 + (-t246 / 0.2e1 + t187 / 0.2e1 + t188 / 0.2e1 - t120 / 0.2e1) * t152 + m(5) * (t116 * t47 + t117 * t46 + t14 * t182 + t15 * t181) + t253 + (Ifges(4,6) * t310 + t57 * mrSges(4,3) - t357 + t80 / 0.2e1 - t26 / 0.2e1 - t27 / 0.2e1 + t3 / 0.2e1 + (m(4) * t57 + t125) * pkin(9) + (-t118 * mrSges(4,3) + t124 / 0.2e1 - t330 * t282 + t331 * t278 + (-m(4) * t118 + m(5) * t108 - t158 + t95) * pkin(9)) * qJD(3)) * t283 + (-Ifges(3,6) * t343 - t284 * t272 / 0.2e1) * t275 - t358 - t359 + (Ifges(4,5) * t310 - t58 * mrSges(4,3) + t356 + t81 / 0.2e1 + t332 * t282 + t333 * t278 + (t278 * t330 + t282 * t331) * qJD(4) + (-t119 * mrSges(4,3) + t334 / 0.2e1 - t123 / 0.2e1 + t63 / 0.2e1 + t64 / 0.2e1 - t35 / 0.2e1) * qJD(3) + (-qJD(3) * t157 - t126 + m(4) * (-qJD(3) * t119 - t58) - t397) * pkin(9)) * t279 + t316 * t78 + t317 * t77 + t318 * t155 + t319 * t154; (t156 * t90 + t19 * t75 + t20 * t74) * t390 + (t104 * t161 + t110 * t162 + t115 * t192) * t391 + (t116 * t182 + t117 * t181) * t392 - 0.2e1 * pkin(2) * t218 + 0.2e1 * t116 * t227 + 0.2e1 * t117 * t228 + 0.2e1 * t110 * t229 + 0.2e1 * t104 * t230 + 0.2e1 * t115 * t208 - t199 * t43 + t200 * t44 + 0.2e1 * t181 * t168 + 0.2e1 * t182 * t170 + 0.2e1 * t192 * t135 + 0.2e1 * t19 * t163 + 0.2e1 * t20 * t164 + 0.2e1 * t162 * t169 + 0.2e1 * t161 * t171 + 0.2e1 * t156 * t45 + 0.2e1 * t90 * t133 + t97 * t121 + t96 * t122 + 0.2e1 * t74 * t82 + 0.2e1 * t75 * t83 + (t136 * t389 + t226 + (t131 + t132) * t282 + (t127 - t130) * t278 + (-t278 * t347 + t282 * t348) * qJD(4) + (pkin(9) ^ 2 * t283 * t392 - t120 + t187 + t188 - t246) * qJD(3)) * t279 + (-t128 - t129 + t223 + t42 + (t209 * t389 + t278 * t348 + t282 * t347 + t249) * qJD(3)) * t283; m(6) * (t16 * t236 + t202 * t48) + t5 * t375 + t4 * t376 + t24 * t379 + t23 * t380 + t36 * t381 + t37 * t382 + t89 * t386 + t88 * t387 + t236 * t33 + t16 * t240 + t48 * t216 + t108 * t217 + t212 * t6 + t202 * t94 + t159 * t18 + t160 * t17 + t176 * t39 + t9 * t148 + t100 * t59 + t101 * t60 + t38 * t86 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t327 + m(7) * (t1 * t160 + t100 * t12 + t101 * t11 + t159 * t2 + t176 * t38 + t212 * t9) + (-t1 * t297 - t11 * t146 - t12 * t147 + t2 * t298) * mrSges(7,3) - t287 * t241 + ((t49 + t52) * t282 + (-t50 + t51) * t278 + ((-t113 + t114) * t282 + (-t111 - t112) * t278) * qJD(4) + m(6) * (t10 * t282 + t13 * t278 + t338 * t41 - t340 * t40) + m(5) * (t14 * t282 - t15 * t278 - t338 * t46 - t340 * t47)) * pkin(10) + (t13 * mrSges(6,2) - t15 * mrSges(5,3) + t332) * t278 + (t10 * mrSges(6,2) + t14 * mrSges(5,3) - t333) * t282 + t307 * t152 + t309 * t203 + t397 * pkin(3) + t312 * t78 + t313 * t77 + t314 * t155 + t315 * t154 + ((t41 * mrSges(6,2) - t46 * mrSges(5,3) - t330) * t282 + (-t40 * mrSges(6,2) - t47 * mrSges(5,3) + t331) * t278) * qJD(4); t279 * pkin(9) * t217 + m(6) * (t115 * t236 + t192 * t202) + t44 * t375 + t43 * t376 + t89 * t377 + t88 * t378 + t96 * t379 + t97 * t380 + t121 * t381 + t122 * t382 + t236 * t135 + t115 * t240 + t192 * t216 + t202 * t208 + t212 * t45 + t159 * t82 + t160 * t83 + t100 * t163 + t101 * t164 + t176 * t133 + t156 * t86 + t90 * t148 - pkin(3) * t136 + m(7) * (t100 * t75 + t101 * t74 + t156 * t176 + t159 * t20 + t160 * t19 + t212 * t90) + (-t146 * t74 - t147 * t75 - t19 * t297 + t20 * t298) * mrSges(7,3) + t272 + ((-m(5) * pkin(3) - mrSges(4,1) + t241) * t372 + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t307) * t279) * qJD(3) + (-t117 * mrSges(5,3) + t110 * mrSges(6,2) + t315 * t279 + t313 * t341 + (-t161 * mrSges(6,2) - t182 * mrSges(5,3) - t279 * t312 + t317) * qJD(4) + (-t168 + t169 + (-t227 - t230) * qJD(4) + m(6) * (-qJD(4) * t161 + t110) + m(5) * (-qJD(4) * t182 - t117)) * pkin(10) + t318) * t278 + (t116 * mrSges(5,3) + t104 * mrSges(6,2) + t314 * t279 + t312 * t341 + (t162 * mrSges(6,2) - t181 * mrSges(5,3) + t279 * t313 + t316) * qJD(4) + (t170 + t171 + (-t228 + t229) * qJD(4) + m(6) * (qJD(4) * t162 + t104) + m(5) * (-qJD(4) * t181 + t116)) * pkin(10) - t319) * t282 - t309 * t283; (t100 * t160 + t101 * t159 + t176 * t212) * t390 + 0.2e1 * t176 * t148 + 0.2e1 * t212 * t86 - t147 * t150 - t297 * t88 + t146 * t151 - t298 * t89 + 0.2e1 * t236 * t216 - 0.2e1 * pkin(3) * t217 + (-t219 + t222) * t282 + (t224 + t225) * t278 + ((t247 + t248) * t282 + (t242 - t245) * t278) * qJD(4) + 0.2e1 * (m(6) * t236 + t240) * t202 + 0.2e1 * (-t100 * t297 + t101 * t298 - t146 * t159 - t147 * t160) * mrSges(7,3); -t290 + t234 * t18 + t235 * t17 + t194 * t60 + t193 * t59 + qJD(5) * t111 - pkin(4) * t51 + qJ(5) * t52 - t13 * mrSges(6,1) - t14 * mrSges(5,2) + t15 * mrSges(5,1) + t10 * mrSges(6,3) + m(7) * (t1 * t235 + t11 * t194 + t12 * t193 + t2 * t234) + m(6) * (-pkin(4) * t13 + qJ(5) * t10 + qJD(5) * t40) + t398; -Ifges(5,6) * t321 + t311 + t234 * t82 + t235 * t83 + qJD(5) * t230 + t194 * t164 + t193 * t163 - pkin(4) * t169 + qJ(5) * t171 - t110 * mrSges(6,1) - t116 * mrSges(5,2) + t117 * mrSges(5,1) + t104 * mrSges(6,3) + m(7) * (t19 * t235 + t193 * t75 + t194 * t74 + t20 * t234) + m(6) * (-pkin(4) * t110 + qJ(5) * t104 + qJD(5) * t161) + (Ifges(7,3) * qJD(3) + (-t362 - t401) * qJD(4)) * t279 + t293 + t345; m(7) * (t100 * t235 + t101 * t234 + t159 * t194 + t160 * t193) - t393 * mrSges(6,2) + (-t146 * t234 - t147 * t235 - t193 * t297 + t194 * t298) * mrSges(7,3) + (m(6) * t336 + (-m(6) * t300 + t240 + t241) * qJD(4)) * pkin(10) - t291 + t220 + t221; (t193 * t235 + t194 * t234) * t390 + 0.2e1 * t361 - 0.2e1 * t360 + 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t277 * t17 + t281 * t18 + (-t277 * t60 + t281 * t59) * qJD(6) + m(7) * (t1 * t277 + t2 * t281 + (-t11 * t277 + t12 * t281) * qJD(6)) + m(6) * t13 + t51; t277 * t83 + t281 * t82 + (t163 * t281 - t164 * t277) * qJD(6) + m(7) * (t19 * t277 + t20 * t281 + (-t277 * t74 + t281 * t75) * qJD(6)) + m(6) * t110 + t169; m(7) * (t277 * t100 + t281 * t101 + (-t159 * t277 + t160 * t281) * qJD(6)) + (m(6) * pkin(10) + mrSges(6,2)) * t338 + (-t281 * t146 - t277 * t147 + (-t277 * t298 - t281 * t297) * qJD(6)) * mrSges(7,3); m(7) * (t193 * t277 + t194 * t281) + (m(7) * (-t234 * t277 + t235 * t281) - t396) * qJD(6); 0; t290; -t293 - t329; t291; t360 - t361; t396 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
