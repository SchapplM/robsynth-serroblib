% Calculate time derivative of joint inertia matrix for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:55
% EndTime: 2019-03-09 21:42:12
% DurationCPUTime: 7.11s
% Computational Cost: add. (6344->770), mult. (17273->1030), div. (0->0), fcn. (15288->8), ass. (0->293)
t377 = Ifges(5,5) + Ifges(7,5);
t260 = sin(qJ(4));
t261 = sin(qJ(3));
t263 = cos(qJ(4));
t323 = qJD(4) * t263;
t264 = cos(qJ(3));
t326 = qJD(3) * t264;
t267 = t260 * t326 + t261 * t323;
t376 = -Ifges(6,4) / 0.2e1;
t325 = qJD(4) * t260;
t327 = qJD(3) * t261;
t373 = pkin(9) * (t263 * t327 + t264 * t325);
t258 = cos(pkin(6));
t257 = sin(pkin(6));
t262 = sin(qJ(2));
t340 = t257 * t262;
t162 = -t258 * t264 + t261 * t340;
t265 = cos(qJ(2));
t329 = qJD(2) * t257;
t313 = t265 * t329;
t118 = -qJD(3) * t162 + t264 * t313;
t339 = t257 * t265;
t372 = -qJD(4) * t339 + t118;
t235 = pkin(8) * t340;
t363 = pkin(1) * t265;
t165 = t258 * t363 - t235;
t163 = t258 * t261 + t264 * t340;
t117 = qJD(3) * t163 + t261 * t313;
t328 = qJD(2) * t262;
t314 = t257 * t328;
t60 = t163 * t323 + t260 * t372 - t263 * t314;
t61 = -t163 * t325 + t260 * t314 + t263 * t372;
t11 = Ifges(5,5) * t61 - Ifges(5,6) * t60 + Ifges(5,3) * t117;
t15 = Ifges(7,1) * t117 + Ifges(7,4) * t60 + Ifges(7,5) * t61;
t16 = Ifges(6,1) * t117 - Ifges(6,4) * t61 + Ifges(6,5) * t60;
t371 = t11 + t15 + t16;
t370 = 2 * m(5);
t369 = 2 * m(6);
t368 = 2 * m(7);
t367 = 0.2e1 * pkin(9);
t366 = 2 * mrSges(7,1);
t365 = -2 * mrSges(3,3);
t364 = pkin(5) + pkin(10);
t362 = pkin(9) * t260;
t361 = pkin(10) * t260;
t360 = pkin(10) * t263;
t255 = t261 * pkin(9);
t359 = -mrSges(6,1) - mrSges(7,1);
t358 = pkin(4) + qJ(6);
t145 = t235 + (-pkin(2) - t363) * t258;
t77 = t162 * pkin(3) - t163 * pkin(10) + t145;
t166 = t258 * t262 * pkin(1) + pkin(8) * t339;
t146 = pkin(9) * t258 + t166;
t147 = (-pkin(2) * t265 - pkin(9) * t262 - pkin(1)) * t257;
t91 = t264 * t146 + t261 * t147;
t79 = -pkin(10) * t339 + t91;
t27 = t260 * t77 + t263 * t79;
t357 = Ifges(4,4) * t261;
t356 = Ifges(4,4) * t264;
t355 = Ifges(5,4) * t260;
t354 = Ifges(5,4) * t263;
t353 = Ifges(6,4) * t263;
t352 = Ifges(5,6) * t260;
t351 = Ifges(5,6) * t263;
t350 = Ifges(6,6) * t260;
t349 = Ifges(6,6) * t263;
t348 = Ifges(7,6) * t260;
t347 = Ifges(7,6) * t263;
t159 = t165 * qJD(2);
t346 = t159 * mrSges(3,2);
t160 = t166 * qJD(2);
t345 = t160 * mrSges(3,1);
t344 = t160 * mrSges(4,1);
t343 = t160 * mrSges(4,2);
t342 = qJ(5) * t260;
t341 = qJD(3) * mrSges(7,3);
t338 = t260 * t261;
t337 = t260 * t264;
t336 = t261 * t263;
t335 = t263 * t264;
t197 = (pkin(3) * t261 - pkin(10) * t264) * qJD(3);
t202 = -pkin(3) * t264 - pkin(10) * t261 - pkin(2);
t334 = t260 * t197 + t202 * t323;
t308 = t263 * t326;
t333 = mrSges(6,1) * t308 + mrSges(6,2) * t327;
t332 = Ifges(5,5) * t308 + Ifges(5,3) * t327;
t331 = pkin(4) * t338 + t255;
t241 = pkin(9) * t335;
t143 = t260 * t202 + t241;
t183 = Ifges(7,4) * t325 + Ifges(7,5) * t323;
t330 = qJ(5) * qJD(5);
t324 = qJD(4) * t261;
t322 = qJD(6) * t260;
t240 = pkin(9) * t337;
t320 = Ifges(4,6) * t339;
t218 = t364 * t263;
t319 = mrSges(7,1) * t325;
t318 = Ifges(4,5) * t118 - Ifges(4,6) * t117 + Ifges(4,3) * t314;
t282 = -Ifges(5,2) * t260 + t354;
t149 = -Ifges(5,6) * t264 + t261 * t282;
t277 = Ifges(6,3) * t260 - t349;
t152 = -Ifges(6,5) * t264 + t261 * t277;
t278 = Ifges(7,2) * t260 + t347;
t153 = -Ifges(7,4) * t264 + t261 * t278;
t317 = -t149 + t152 + t153;
t283 = Ifges(5,1) * t263 - t355;
t150 = -Ifges(5,5) * t264 + t261 * t283;
t274 = Ifges(7,3) * t263 + t348;
t151 = -Ifges(7,5) * t264 + t261 * t274;
t281 = -Ifges(6,2) * t263 + t350;
t154 = -Ifges(6,4) * t264 + t261 * t281;
t316 = t150 + t151 - t154;
t315 = -pkin(4) - t362;
t311 = t260 * t324;
t34 = t61 * mrSges(6,1) + t117 * mrSges(6,2);
t33 = -t60 * mrSges(7,1) + t117 * mrSges(7,2);
t31 = t61 * mrSges(7,1) - t117 * mrSges(7,3);
t26 = -t260 * t79 + t263 * t77;
t90 = -t261 * t146 + t147 * t264;
t142 = t202 * t263 - t240;
t307 = pkin(4) * t325 - qJD(5) * t260;
t306 = pkin(4) * t267 + pkin(9) * t326 + qJ(5) * t311;
t305 = Ifges(7,1) * t327 + Ifges(7,4) * t267 + Ifges(7,5) * t308;
t304 = Ifges(6,1) * t327 + Ifges(6,4) * t311 + Ifges(6,5) * t267;
t13 = Ifges(6,4) * t117 - Ifges(6,2) * t61 + Ifges(6,6) * t60;
t17 = Ifges(5,1) * t61 - Ifges(5,4) * t60 + Ifges(5,5) * t117;
t9 = Ifges(7,5) * t117 + Ifges(7,6) * t60 + Ifges(7,3) * t61;
t303 = -t13 / 0.2e1 + t17 / 0.2e1 + t9 / 0.2e1;
t302 = t314 / 0.2e1;
t10 = Ifges(6,5) * t117 - Ifges(6,6) * t61 + Ifges(6,3) * t60;
t12 = Ifges(7,4) * t117 + Ifges(7,2) * t60 + Ifges(7,6) * t61;
t14 = Ifges(5,4) * t61 - Ifges(5,2) * t60 + Ifges(5,6) * t117;
t301 = t12 / 0.2e1 - t14 / 0.2e1 + t10 / 0.2e1;
t119 = t163 * t260 + t263 * t339;
t120 = t163 * t263 - t260 * t339;
t43 = Ifges(7,5) * t162 + Ifges(7,6) * t119 + Ifges(7,3) * t120;
t47 = Ifges(6,4) * t162 - Ifges(6,2) * t120 + Ifges(6,6) * t119;
t51 = Ifges(5,1) * t120 - Ifges(5,4) * t119 + Ifges(5,5) * t162;
t300 = t43 / 0.2e1 - t47 / 0.2e1 + t51 / 0.2e1;
t44 = Ifges(6,5) * t162 - Ifges(6,6) * t120 + Ifges(6,3) * t119;
t46 = Ifges(7,4) * t162 + Ifges(7,2) * t119 + Ifges(7,6) * t120;
t48 = Ifges(5,4) * t120 - Ifges(5,2) * t119 + Ifges(5,6) * t162;
t299 = t44 / 0.2e1 + t46 / 0.2e1 - t48 / 0.2e1;
t24 = -qJ(5) * t162 - t27;
t280 = Ifges(6,2) * t260 + t349;
t102 = t280 * t324 + (Ifges(6,4) * t261 + t264 * t281) * qJD(3);
t215 = Ifges(5,1) * t260 + t354;
t98 = -t215 * t324 + (Ifges(5,5) * t261 + t264 * t283) * qJD(3);
t275 = -Ifges(7,3) * t260 + t347;
t99 = t275 * t324 + (Ifges(7,5) * t261 + t264 * t274) * qJD(3);
t298 = t98 / 0.2e1 + t99 / 0.2e1 - t102 / 0.2e1;
t276 = Ifges(6,3) * t263 + t350;
t100 = t276 * t324 + (Ifges(6,5) * t261 + t264 * t277) * qJD(3);
t279 = Ifges(7,2) * t263 - t348;
t101 = t279 * t324 + (Ifges(7,4) * t261 + t264 * t278) * qJD(3);
t213 = Ifges(5,2) * t263 + t355;
t97 = -t213 * t324 + (Ifges(5,6) * t261 + t264 * t282) * qJD(3);
t297 = t100 / 0.2e1 + t101 / 0.2e1 - t97 / 0.2e1;
t125 = qJ(5) * t264 - t143;
t78 = pkin(3) * t339 - t90;
t296 = qJD(4) * t241 - t197 * t263 + t202 * t325;
t295 = -t149 / 0.2e1 + t152 / 0.2e1 + t153 / 0.2e1;
t294 = t150 / 0.2e1 + t151 / 0.2e1 - t154 / 0.2e1;
t178 = t274 * qJD(4);
t181 = t281 * qJD(4);
t187 = t283 * qJD(4);
t293 = t178 / 0.2e1 - t181 / 0.2e1 + t187 / 0.2e1;
t179 = t277 * qJD(4);
t180 = t278 * qJD(4);
t185 = t282 * qJD(4);
t292 = t179 / 0.2e1 + t180 / 0.2e1 - t185 / 0.2e1;
t247 = Ifges(6,5) * t325;
t248 = Ifges(5,5) * t323;
t291 = -Ifges(5,6) * t325 / 0.2e1 + t248 / 0.2e1 + t183 / 0.2e1 + t323 * t376 + t247 / 0.2e1;
t290 = -t275 / 0.2e1 + t280 / 0.2e1 + t215 / 0.2e1;
t289 = -t276 / 0.2e1 - t279 / 0.2e1 - t213 / 0.2e1;
t288 = t351 / 0.2e1 - (Ifges(7,4) + Ifges(6,5)) * t263 / 0.2e1 + (t376 + t377 / 0.2e1) * t260;
t287 = -qJD(5) * t264 + t334;
t204 = -t263 * mrSges(5,1) + t260 * mrSges(5,2);
t286 = mrSges(5,1) * t260 + mrSges(5,2) * t263;
t203 = t263 * mrSges(6,2) - t260 * mrSges(6,3);
t285 = -mrSges(6,2) * t260 - mrSges(6,3) * t263;
t284 = -mrSges(7,2) * t263 + mrSges(7,3) * t260;
t273 = -pkin(4) * t263 - t342;
t158 = (pkin(2) * t262 - pkin(9) * t265) * t329;
t39 = -t146 * t327 + t147 * t326 + t261 * t158 + t264 * t159;
t37 = pkin(10) * t314 + t39;
t52 = t117 * pkin(3) - t118 * pkin(10) + t160;
t7 = -t260 * t37 + t263 * t52 - t79 * t323 - t77 * t325;
t272 = -qJ(5) * t263 + qJ(6) * t260;
t271 = t263 * (m(6) * pkin(10) - t359);
t40 = -t146 * t326 - t147 * t327 + t158 * t264 - t261 * t159;
t6 = t260 * t52 + t263 * t37 + t77 * t323 - t325 * t79;
t270 = -qJ(5) * t120 + t78;
t269 = t308 - t311;
t38 = -pkin(3) * t314 - t40;
t134 = -mrSges(7,1) * t267 + mrSges(7,2) * t327;
t4 = -qJ(5) * t117 - qJD(5) * t162 - t6;
t266 = -qJ(5) * t61 - qJD(5) * t120 + t38;
t256 = t264 * pkin(4);
t249 = Ifges(4,5) * t326;
t223 = mrSges(7,1) * t308;
t221 = Ifges(3,5) * t313;
t217 = t364 * t260;
t216 = Ifges(4,1) * t261 + t356;
t214 = Ifges(4,2) * t264 + t357;
t205 = -mrSges(7,2) * t260 - mrSges(7,3) * t263;
t198 = -pkin(3) + t273;
t196 = qJD(4) * t218;
t195 = t364 * t325;
t194 = mrSges(6,1) * t336 - mrSges(6,2) * t264;
t193 = -mrSges(7,1) * t338 - mrSges(7,2) * t264;
t192 = mrSges(6,1) * t338 + mrSges(6,3) * t264;
t191 = mrSges(7,1) * t336 + mrSges(7,3) * t264;
t190 = -mrSges(5,1) * t264 - mrSges(5,3) * t336;
t189 = mrSges(5,2) * t264 - mrSges(5,3) * t338;
t188 = (Ifges(4,1) * t264 - t357) * qJD(3);
t186 = (-Ifges(4,2) * t261 + t356) * qJD(3);
t177 = (mrSges(4,1) * t261 + mrSges(4,2) * t264) * qJD(3);
t176 = t286 * qJD(4);
t175 = t285 * qJD(4);
t174 = t284 * qJD(4);
t170 = t286 * t261;
t169 = t285 * t261;
t168 = t284 * t261;
t167 = -t263 * t358 - pkin(3) - t342;
t161 = -qJ(5) * t323 + t307;
t157 = -qJ(5) * t336 + t331;
t156 = -Ifges(6,1) * t264 + (Ifges(6,5) * t260 - t353) * t261;
t155 = -Ifges(7,1) * t264 + (Ifges(7,4) * t260 + Ifges(7,5) * t263) * t261;
t148 = -Ifges(5,3) * t264 + (Ifges(5,5) * t263 - t352) * t261;
t135 = -mrSges(6,1) * t311 + t333;
t133 = mrSges(6,1) * t267 - mrSges(6,3) * t327;
t132 = t223 + (-t319 - t341) * t261;
t131 = -mrSges(5,2) * t327 - mrSges(5,3) * t267;
t130 = mrSges(5,1) * t327 - mrSges(5,3) * t269;
t126 = -t142 + t256;
t124 = -mrSges(4,1) * t339 - t163 * mrSges(4,3);
t123 = mrSges(4,2) * t339 - t162 * mrSges(4,3);
t122 = t261 * t272 + t331;
t121 = qJD(4) * t272 - qJD(6) * t263 + t307;
t109 = -pkin(5) * t338 - t125;
t108 = -mrSges(7,2) * t269 + mrSges(7,3) * t267;
t107 = mrSges(5,1) * t267 + mrSges(5,2) * t269;
t106 = -mrSges(6,2) * t267 - mrSges(6,3) * t269;
t105 = qJ(6) * t264 + t240 + t256 + (pkin(5) * t261 - t202) * t263;
t104 = -Ifges(6,4) * t308 + t304;
t103 = -Ifges(7,5) * t311 + t305;
t96 = -Ifges(5,5) * t311 - Ifges(5,6) * t267 + t332;
t95 = mrSges(4,1) * t314 - mrSges(4,3) * t118;
t94 = -mrSges(4,2) * t314 - mrSges(4,3) * t117;
t93 = Ifges(4,1) * t163 - Ifges(4,4) * t162 - Ifges(4,5) * t339;
t92 = Ifges(4,4) * t163 - Ifges(4,2) * t162 - t320;
t89 = t327 * t362 - t296;
t88 = t334 - t373;
t87 = (-qJ(5) * t326 - qJD(5) * t261) * t263 + t306;
t86 = mrSges(5,1) * t162 - mrSges(5,3) * t120;
t85 = -mrSges(5,2) * t162 - mrSges(5,3) * t119;
t84 = mrSges(6,1) * t120 + mrSges(6,2) * t162;
t83 = -mrSges(7,1) * t119 + mrSges(7,2) * t162;
t82 = mrSges(6,1) * t119 - mrSges(6,3) * t162;
t81 = mrSges(7,1) * t120 - mrSges(7,3) * t162;
t80 = t315 * t327 + t296;
t74 = -qJ(5) * t327 - t287 + t373;
t69 = -mrSges(6,2) * t119 - mrSges(6,3) * t120;
t68 = mrSges(5,1) * t119 + mrSges(5,2) * t120;
t67 = -mrSges(7,2) * t120 + mrSges(7,3) * t119;
t66 = mrSges(4,1) * t117 + mrSges(4,2) * t118;
t65 = Ifges(4,1) * t118 - Ifges(4,4) * t117 + Ifges(4,5) * t314;
t64 = Ifges(4,4) * t118 - Ifges(4,2) * t117 + Ifges(4,6) * t314;
t63 = t272 * t326 + (t322 + (qJ(6) * qJD(4) - qJD(5)) * t263) * t261 + t306;
t62 = (-pkin(5) * t336 - t240) * qJD(4) + (-pkin(5) * t337 + (-pkin(9) * t263 + qJ(5)) * t261) * qJD(3) + t287;
t50 = Ifges(6,1) * t162 - Ifges(6,4) * t120 + Ifges(6,5) * t119;
t49 = Ifges(7,1) * t162 + Ifges(7,4) * t119 + Ifges(7,5) * t120;
t45 = Ifges(5,5) * t120 - Ifges(5,6) * t119 + Ifges(5,3) * t162;
t42 = -pkin(5) * t311 + qJD(6) * t264 + (pkin(5) * t335 + (-qJ(6) + t315) * t261) * qJD(3) + t296;
t32 = mrSges(6,1) * t60 - mrSges(6,3) * t117;
t30 = mrSges(5,1) * t117 - mrSges(5,3) * t61;
t29 = -mrSges(5,2) * t117 - mrSges(5,3) * t60;
t28 = pkin(4) * t119 + t270;
t25 = -pkin(4) * t162 - t26;
t23 = t119 * t358 + t270;
t22 = -mrSges(6,2) * t60 - mrSges(6,3) * t61;
t21 = mrSges(5,1) * t60 + mrSges(5,2) * t61;
t20 = -mrSges(7,2) * t61 + mrSges(7,3) * t60;
t19 = -pkin(5) * t119 - t24;
t18 = pkin(5) * t120 - t162 * t358 - t26;
t8 = pkin(4) * t60 + t266;
t5 = -pkin(4) * t117 - t7;
t3 = qJD(6) * t119 + t358 * t60 + t266;
t2 = -pkin(5) * t60 - t4;
t1 = pkin(5) * t61 - qJD(6) * t162 - t117 * t358 - t7;
t35 = [(-t64 + 0.2e1 * t344 + t371) * t162 + (t1 * t18 + t19 * t2 + t23 * t3) * t368 + (t24 * t4 + t25 * t5 + t28 * t8) * t369 + (t26 * t7 + t27 * t6 + t38 * t78) * t370 + 0.2e1 * m(3) * (t159 * t166 - t160 * t165) + 0.2e1 * m(4) * (t145 * t160 + t39 * t91 + t40 * t90) + (t43 + t51 - t47) * t61 + (t44 + t46 - t48) * t60 + (t65 + 0.2e1 * t343) * t163 + (t221 - 0.2e1 * t345 - 0.2e1 * t346) * t258 + (t9 + t17 - t13) * t120 + 0.2e1 * t23 * t20 + 0.2e1 * t28 * t22 + 0.2e1 * t27 * t29 + 0.2e1 * t26 * t30 + 0.2e1 * t18 * t31 + 0.2e1 * t24 * t32 + 0.2e1 * t19 * t33 + 0.2e1 * t25 * t34 + 0.2e1 * t3 * t67 + 0.2e1 * t38 * t68 + 0.2e1 * t8 * t69 + 0.2e1 * t78 * t21 + 0.2e1 * t1 * t81 + 0.2e1 * t4 * t82 + 0.2e1 * t2 * t83 + 0.2e1 * t5 * t84 + 0.2e1 * t6 * t85 + 0.2e1 * t7 * t86 + 0.2e1 * t91 * t94 + 0.2e1 * t90 * t95 + t118 * t93 + 0.2e1 * t39 * t123 + 0.2e1 * t40 * t124 + 0.2e1 * t145 * t66 + (t45 + t49 + t50 - t92) * t117 + (-t265 * t318 + 0.2e1 * (t159 * t265 + t160 * t262) * mrSges(3,3) + ((t165 * t365 + Ifges(3,5) * t258 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t265) * t257) * t265 + (t166 * t365 + Ifges(4,5) * t163 - 0.2e1 * Ifges(3,6) * t258 - Ifges(4,6) * t162 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t262 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t265) * t257) * t262) * qJD(2)) * t257 + (t10 + t12 - t14) * t119; -t345 - t346 + (t96 / 0.2e1 + t103 / 0.2e1 + t104 / 0.2e1 - t186 / 0.2e1) * t162 + m(5) * (t142 * t7 + t143 * t6 + t26 * t89 + t27 * t88) + m(6) * (t125 * t4 + t126 * t5 + t157 * t8 + t24 * t74 + t25 * t80 + t28 * t87) + m(7) * (t1 * t105 + t109 * t2 + t122 * t3 + t18 * t42 + t19 * t62 + t23 * t63) + (t148 / 0.2e1 + t155 / 0.2e1 + t156 / 0.2e1 - t214 / 0.2e1) * t117 + (-m(4) * t160 - t66) * pkin(2) + (-t40 * mrSges(4,3) + t343 + t65 / 0.2e1 + Ifges(4,5) * t302 + t303 * t263 + t301 * t260 + (-t260 * t300 + t263 * t299) * qJD(4) + (-t92 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1 + t45 / 0.2e1 + t320 / 0.2e1 - t91 * mrSges(4,3)) * qJD(3) + (-qJD(3) * t123 + t21 - t95 + m(5) * t38 + m(4) * (-qJD(3) * t91 - t40)) * pkin(9)) * t261 + (t39 * mrSges(4,3) + Ifges(4,6) * t302 - t11 / 0.2e1 - t15 / 0.2e1 - t16 / 0.2e1 + t64 / 0.2e1 - t344 + (m(4) * t39 + t94) * pkin(9) + (t93 / 0.2e1 - t90 * mrSges(4,3) + t300 * t263 + t299 * t260 + (-m(4) * t90 + m(5) * t78 - t124 + t68) * pkin(9)) * qJD(3)) * t264 + (-t265 * t249 / 0.2e1 - Ifges(3,6) * t328) * t257 + t221 + t63 * t67 + t42 * t81 + t74 * t82 + t62 * t83 + t80 * t84 + t87 * t69 + t88 * t85 + t89 * t86 + t105 * t31 + t28 * t106 + t78 * t107 + t23 * t108 + t109 * t33 + t122 * t20 + t125 * t32 + t294 * t61 + t295 * t60 + t126 * t34 + t26 * t130 + t27 * t131 + t18 * t132 + t24 * t133 + t19 * t134 + t25 * t135 + t297 * t119 + t142 * t30 + t143 * t29 + t157 * t22 + t298 * t120 + t3 * t168 + t8 * t169 + t38 * t170 + t145 * t177 + t163 * t188 / 0.2e1 + t6 * t189 + t7 * t190 + t1 * t191 + t4 * t192 + t2 * t193 + t5 * t194 + t118 * t216 / 0.2e1; 0.2e1 * t122 * t108 + 0.2e1 * t105 * t132 + 0.2e1 * t125 * t133 + 0.2e1 * t109 * t134 + 0.2e1 * t126 * t135 + 0.2e1 * t142 * t130 + 0.2e1 * t143 * t131 + 0.2e1 * t157 * t106 + 0.2e1 * t63 * t168 + 0.2e1 * t87 * t169 - 0.2e1 * pkin(2) * t177 + 0.2e1 * t88 * t189 + 0.2e1 * t89 * t190 + 0.2e1 * t42 * t191 + 0.2e1 * t74 * t192 + 0.2e1 * t62 * t193 + 0.2e1 * t80 * t194 + (t142 * t89 + t143 * t88) * t370 + (t125 * t74 + t126 * t80 + t157 * t87) * t369 + (t105 * t42 + t109 * t62 + t122 * t63) * t368 + (-t103 - t104 + t186 - t96 + (t170 * t367 + t260 * t317 + t263 * t316 + t216) * qJD(3)) * t264 + (t107 * t367 + t188 + (t98 + t99 - t102) * t263 + (-t97 + t100 + t101) * t260 + (pkin(9) ^ 2 * t264 * t370 + t148 + t155 + t156 - t214) * qJD(3) + (-t260 * t316 + t263 * t317) * qJD(4)) * t261; t318 + m(5) * (-pkin(3) * t38 + t360 * t6 - t361 * t7) + m(6) * (t161 * t28 + t198 * t8 - t360 * t4 + t361 * t5) - pkin(3) * t21 - t39 * mrSges(4,2) + t40 * mrSges(4,1) + t121 * t67 + t288 * t117 + t289 * t60 + t290 * t61 + t291 * t162 + t292 * t119 + t293 * t120 + t161 * t69 + ((t25 * mrSges(6,1) + t18 * mrSges(7,1) - t26 * mrSges(5,3) + t300) * t263 + (t24 * mrSges(6,1) - t19 * mrSges(7,1) - t27 * mrSges(5,3) + t299) * t260 + ((t84 - t86) * t263 + (t82 - t85) * t260 + m(6) * (t24 * t260 + t25 * t263) + m(5) * (-t26 * t263 - t260 * t27)) * pkin(10)) * qJD(4) + (t2 * mrSges(7,1) - t4 * mrSges(6,1) + t6 * mrSges(5,3) + (t29 - t32) * pkin(10) - t301) * t263 + t167 * t20 + t23 * t174 + t28 * t175 + t78 * t176 - t195 * t83 + t196 * t81 + t198 * t22 + t8 * t203 + t38 * t204 + t3 * t205 + t217 * t31 + t218 * t33 + (-t7 * mrSges(5,3) + t5 * mrSges(6,1) + t1 * mrSges(7,1) + (-t30 + t34) * pkin(10) + t303) * t260 + m(7) * (t1 * t217 + t121 * t23 + t167 * t3 + t18 * t196 - t19 * t195 + t2 * t218); t176 * t255 - pkin(3) * t107 + t198 * t106 + t167 * t108 + t121 * t168 + t122 * t174 + t217 * t132 + t218 * t134 + t157 * t175 + t161 * t169 + t196 * t191 - t195 * t193 + t87 * t203 + t63 * t205 + t249 + m(7) * (t105 * t196 - t109 * t195 + t121 * t122 + t167 * t63 + t217 * t42 + t218 * t62) + m(6) * (t157 * t161 + t198 * t87) - t291 * t264 + ((-m(5) * pkin(3) - mrSges(4,1) + t204) * t264 * pkin(9) + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t288) * t261) * qJD(3) + (t88 * mrSges(5,3) - t74 * mrSges(6,1) + t62 * mrSges(7,1) + t293 * t261 + t290 * t326 + (t126 * mrSges(6,1) + t105 * mrSges(7,1) - t142 * mrSges(5,3) + t289 * t261 + t294) * qJD(4) + (t131 - t133 + (-t190 + t194) * qJD(4) + m(5) * (-qJD(4) * t142 + t88) + m(6) * (qJD(4) * t126 - t74)) * pkin(10) - t297) * t263 + (-t89 * mrSges(5,3) + t80 * mrSges(6,1) + t42 * mrSges(7,1) + t292 * t261 + t289 * t326 + (t125 * mrSges(6,1) - t109 * mrSges(7,1) - t143 * mrSges(5,3) - t290 * t261 + t295) * qJD(4) + (-t130 + t135 + (-t189 + t192) * qJD(4) + m(5) * (-qJD(4) * t143 - t89) + m(6) * (qJD(4) * t125 + t80)) * pkin(10) + t298) * t260; 0.2e1 * t198 * t175 + 0.2e1 * t121 * t205 + 0.2e1 * t167 * t174 + (t121 * t167 - t195 * t218 + t196 * t217) * t368 - 0.2e1 * pkin(3) * t176 + 0.2e1 * (m(6) * t198 + t203) * t161 + (-t195 * t366 - t179 - t180 + t185) * t263 + (t196 * t366 + t178 - t181 + t187) * t260 + ((t217 * t366 + t215 - t275 + t280) * t263 + (-0.2e1 * mrSges(7,1) * t218 - t213 - t276 - t279) * t260) * qJD(4); (-t82 + t83) * qJD(5) + (-t32 + t33) * qJ(5) + m(6) * (-pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t24) + m(7) * (qJ(5) * t2 + qJD(5) * t19 - qJD(6) * t18 - t1 * t358) - t1 * mrSges(7,3) + t2 * mrSges(7,2) - t4 * mrSges(6,3) + t5 * mrSges(6,2) - t6 * mrSges(5,2) + t7 * mrSges(5,1) - pkin(4) * t34 - qJD(6) * t81 - t358 * t31 + t371; t304 + t305 + (-t352 - t353) * t326 + (-t377 * t260 - t351) * t324 + (-t192 + t193) * qJD(5) + (-t133 + t134) * qJ(5) + m(7) * (qJ(5) * t62 + qJD(5) * t109 - qJD(6) * t105 - t358 * t42) + m(6) * (-pkin(4) * t80 - qJ(5) * t74 - qJD(5) * t125) - t42 * mrSges(7,3) + t62 * mrSges(7,2) - t74 * mrSges(6,3) + t80 * mrSges(6,2) - t88 * mrSges(5,2) + t89 * mrSges(5,1) - pkin(4) * t135 - qJD(6) * t191 - t358 * t132 + t332; m(7) * (-qJ(5) * t195 + qJD(5) * t218 - qJD(6) * t217 - t196 * t358) + t247 + t248 - t195 * mrSges(7,2) - t196 * mrSges(7,3) - mrSges(7,1) * t322 + qJD(5) * t271 + ((-mrSges(6,1) * pkin(4) - mrSges(7,1) * t358 - Ifges(6,4)) * t263 + (qJ(5) * t359 - Ifges(5,6)) * t260 + (m(6) * t273 + t203 + t204) * pkin(10)) * qJD(4) + t183; 0.2e1 * m(6) * t330 + 0.2e1 * qJD(6) * mrSges(7,3) + 0.2e1 * m(7) * (qJD(6) * t358 + t330) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJD(5); m(6) * t5 + m(7) * t1 + t31 + t34; t223 + m(6) * t80 + m(7) * t42 + (t325 * t359 - t341) * t261 + t333; m(7) * t196 + qJD(4) * t271; -m(7) * qJD(6); 0; m(7) * t2 + t33; m(7) * t62 + t134; -m(7) * t195 - t319; m(7) * qJD(5); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t35(1) t35(2) t35(4) t35(7) t35(11) t35(16); t35(2) t35(3) t35(5) t35(8) t35(12) t35(17); t35(4) t35(5) t35(6) t35(9) t35(13) t35(18); t35(7) t35(8) t35(9) t35(10) t35(14) t35(19); t35(11) t35(12) t35(13) t35(14) t35(15) t35(20); t35(16) t35(17) t35(18) t35(19) t35(20) t35(21);];
Mq  = res;
