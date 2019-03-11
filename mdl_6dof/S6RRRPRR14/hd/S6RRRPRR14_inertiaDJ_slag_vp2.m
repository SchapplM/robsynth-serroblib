% Calculate time derivative of joint inertia matrix for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:27
% EndTime: 2019-03-09 20:12:42
% DurationCPUTime: 6.97s
% Computational Cost: add. (10244->717), mult. (25822->1032), div. (0->0), fcn. (24271->10), ass. (0->293)
t351 = pkin(4) + pkin(9);
t365 = Ifges(5,1) + Ifges(4,3);
t364 = Ifges(4,5) - Ifges(5,4);
t363 = -Ifges(4,6) + Ifges(5,5);
t255 = sin(qJ(5));
t259 = cos(qJ(5));
t260 = cos(qJ(3));
t310 = qJD(5) * t260;
t256 = sin(qJ(3));
t314 = qJD(3) * t256;
t268 = t255 * t314 - t259 * t310;
t297 = t255 * t310;
t269 = t259 * t314 + t297;
t254 = sin(qJ(6));
t258 = cos(qJ(6));
t277 = t254 * t259 + t258 * t255;
t342 = -t277 / 0.2e1;
t320 = t258 * t259;
t276 = t254 * t255 - t320;
t341 = -t276 / 0.2e1;
t339 = -t255 / 0.2e1;
t325 = qJ(4) * t256;
t352 = pkin(3) + pkin(10);
t188 = -t260 * t352 - pkin(2) - t325;
t224 = t351 * t256;
t195 = t255 * t224;
t128 = t259 * t188 + t195;
t252 = sin(pkin(6));
t257 = sin(qJ(2));
t324 = t252 * t257;
t237 = pkin(8) * t324;
t253 = cos(pkin(6));
t261 = cos(qJ(2));
t337 = pkin(1) * t261;
t185 = t253 * t337 - t237;
t362 = qJD(5) + qJD(6);
t309 = qJD(6) * t254;
t312 = qJD(5) * t255;
t135 = -t254 * t312 - t255 * t309 + t320 * t362;
t136 = t362 * t277;
t361 = qJD(6) * (-t254 * t276 - t258 * t277) - t135 * t254 + t136 * t258;
t360 = t260 * t362;
t359 = 2 * m(6);
t358 = 2 * m(7);
t357 = -2 * mrSges(3,3);
t356 = 2 * mrSges(3,3);
t182 = t253 * t256 + t260 * t324;
t316 = qJD(2) * t252;
t300 = t261 * t316;
t141 = qJD(3) * t182 + t256 * t300;
t303 = t256 * t324;
t181 = -t253 * t260 + t303;
t323 = t252 * t261;
t143 = t181 * t259 + t255 * t323;
t315 = qJD(2) * t257;
t301 = t252 * t315;
t73 = qJD(5) * t143 + t141 * t255 + t259 * t301;
t355 = t73 / 0.2e1;
t270 = -t181 * t255 + t259 * t323;
t85 = t143 * t258 + t254 * t270;
t354 = t85 / 0.2e1;
t86 = t143 * t254 - t258 * t270;
t353 = t86 / 0.2e1;
t350 = -t135 / 0.2e1;
t349 = -t136 / 0.2e1;
t139 = -Ifges(7,4) * t276 - Ifges(7,2) * t277;
t348 = t139 / 0.2e1;
t140 = -Ifges(7,1) * t276 - Ifges(7,4) * t277;
t347 = t140 / 0.2e1;
t346 = t143 / 0.2e1;
t345 = -t270 / 0.2e1;
t176 = t276 * t260;
t344 = t176 / 0.2e1;
t177 = t277 * t260;
t343 = -t177 / 0.2e1;
t332 = Ifges(6,4) * t259;
t220 = -Ifges(6,2) * t255 + t332;
t340 = t220 / 0.2e1;
t338 = -t259 / 0.2e1;
t336 = pkin(11) + t352;
t186 = t253 * t257 * pkin(1) + pkin(8) * t323;
t168 = pkin(9) * t253 + t186;
t169 = (-pkin(2) * t261 - pkin(9) * t257 - pkin(1)) * t252;
t103 = -t256 * t168 + t169 * t260;
t100 = pkin(3) * t323 - t103;
t66 = pkin(4) * t182 + pkin(10) * t323 + t100;
t167 = t237 + (-pkin(2) - t337) * t253;
t267 = -qJ(4) * t182 + t167;
t74 = t181 * t352 + t267;
t30 = t255 * t66 + t259 * t74;
t335 = Ifges(4,4) * t256;
t334 = Ifges(4,4) * t260;
t333 = Ifges(6,4) * t255;
t331 = Ifges(5,6) * t256;
t330 = Ifges(5,6) * t260;
t174 = t185 * qJD(2);
t329 = t174 * mrSges(3,2);
t175 = t186 * qJD(2);
t328 = t175 * mrSges(3,1);
t327 = t175 * mrSges(4,1);
t326 = t175 * mrSges(4,2);
t222 = Ifges(6,1) * t259 - t333;
t322 = t255 * t222;
t321 = t255 * t352;
t319 = t259 * t220;
t318 = t259 * t260;
t317 = t259 * t352;
t88 = -Ifges(7,5) * t136 - Ifges(7,6) * t135;
t104 = t260 * t168 + t256 * t169;
t225 = t351 * t260;
t313 = qJD(3) * t260;
t311 = qJD(5) * t259;
t308 = qJD(6) * t258;
t307 = 2 * mrSges(7,3);
t306 = 0.2e1 * t252;
t142 = -qJD(3) * t303 + (qJD(3) * t253 + t300) * t260;
t72 = qJD(5) * t270 + t141 * t259 - t255 * t301;
t21 = qJD(6) * t85 + t254 * t72 + t258 * t73;
t22 = -qJD(6) * t86 - t254 * t73 + t258 * t72;
t6 = Ifges(7,5) * t21 + Ifges(7,6) * t22 + Ifges(7,3) * t142;
t26 = Ifges(6,5) * t73 + Ifges(6,6) * t72 + Ifges(6,3) * t142;
t96 = t276 * t360 + t277 * t314;
t97 = -t276 * t314 + t277 * t360;
t41 = Ifges(7,5) * t96 + Ifges(7,6) * t97 + Ifges(7,3) * t313;
t304 = m(5) * pkin(9) + mrSges(5,1);
t286 = -Ifges(6,5) * t255 - Ifges(6,6) * t259;
t302 = t286 * qJD(5) / 0.2e1 + t88 / 0.2e1;
t295 = t323 / 0.2e1;
t294 = Ifges(6,5) * t259 / 0.2e1 + Ifges(6,6) * t339 + Ifges(7,5) * t341 + Ifges(7,6) * t342;
t212 = t336 * t259;
t293 = pkin(11) * t260 - t188;
t29 = -t255 * t74 + t259 * t66;
t114 = t142 * mrSges(5,1) + mrSges(5,2) * t301;
t292 = -t136 * mrSges(7,1) - t135 * mrSges(7,2);
t290 = pkin(3) * t314 - qJD(4) * t256;
t163 = (pkin(10) * t256 - qJ(4) * t260) * qJD(3) + t290;
t210 = t351 * t313;
t291 = -t163 * t255 + t259 * t210;
t289 = mrSges(6,1) * t259 - mrSges(6,2) * t255;
t216 = mrSges(6,1) * t255 + mrSges(6,2) * t259;
t215 = t260 * mrSges(5,2) - t256 * mrSges(5,3);
t288 = Ifges(6,1) * t255 + t332;
t287 = Ifges(6,2) * t259 + t333;
t285 = -pkin(3) * t260 - t325;
t24 = pkin(5) * t182 + pkin(11) * t270 + t29;
t25 = pkin(11) * t143 + t30;
t12 = t24 * t258 - t25 * t254;
t13 = t24 * t254 + t25 * t258;
t284 = t255 * t29 - t259 * t30;
t101 = -mrSges(6,2) * t182 + mrSges(6,3) * t143;
t102 = mrSges(6,1) * t182 + mrSges(6,3) * t270;
t283 = t259 * t101 - t255 * t102;
t196 = t259 * t224;
t112 = pkin(5) * t256 + t255 * t293 + t196;
t121 = -pkin(11) * t318 + t128;
t56 = t112 * t258 - t121 * t254;
t57 = t112 * t254 + t121 * t258;
t127 = -t188 * t255 + t196;
t282 = -t127 * t255 + t128 * t259;
t281 = -t135 * t277 - t136 * t276;
t207 = mrSges(6,3) * t255 * t260 + mrSges(6,1) * t256;
t208 = -mrSges(6,2) * t256 - mrSges(6,3) * t318;
t278 = -t255 * t207 + t259 * t208;
t211 = t336 * t255;
t146 = -t211 * t258 - t212 * t254;
t145 = t211 * t254 - t212 * t258;
t191 = t336 * t312;
t192 = qJD(5) * t212;
t83 = qJD(6) * t145 + t191 * t254 - t192 * t258;
t84 = -qJD(6) * t146 + t191 * t258 + t192 * t254;
t275 = t84 * mrSges(7,1) - t83 * mrSges(7,2) + t88;
t99 = qJ(4) * t323 - t104;
t173 = (pkin(2) * t257 - pkin(9) * t261) * t316;
t54 = -t168 * t313 - t169 * t314 + t173 * t260 - t256 * t174;
t37 = pkin(4) * t142 - t301 * t352 - t54;
t263 = -qJ(4) * t142 - qJD(4) * t182 + t175;
t44 = t141 * t352 + t263;
t11 = -qJD(5) * t30 - t255 * t44 + t259 * t37;
t4 = pkin(5) * t142 - pkin(11) * t73 + t11;
t10 = t255 * t37 + t259 * t44 + t66 * t311 - t312 * t74;
t5 = pkin(11) * t72 + t10;
t2 = qJD(6) * t12 + t254 * t4 + t258 * t5;
t3 = -qJD(6) * t13 - t254 * t5 + t258 * t4;
t274 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t6;
t52 = (-pkin(11) * t255 * t256 + pkin(5) * t260) * qJD(3) + (t259 * t293 - t195) * qJD(5) + t291;
t70 = t259 * t163 - t188 * t312 + t255 * t210 + t224 * t311;
t55 = pkin(11) * t269 + t70;
t15 = qJD(6) * t56 + t254 * t52 + t258 * t55;
t16 = -qJD(6) * t57 - t254 * t55 + t258 * t52;
t273 = t16 * mrSges(7,1) - t15 * mrSges(7,2) + t41;
t64 = -Ifges(6,4) * t270 + Ifges(6,2) * t143 + Ifges(6,6) * t182;
t65 = -Ifges(6,1) * t270 + Ifges(6,4) * t143 + Ifges(6,5) * t182;
t272 = t338 * t64 + t339 * t65;
t271 = t363 * t141 + t364 * t142 + t365 * t301;
t53 = -t168 * t314 + t169 * t313 + t256 * t173 + t260 * t174;
t75 = -pkin(4) * t181 - t99;
t117 = t268 * Ifges(6,5) + t269 * Ifges(6,6) + Ifges(6,3) * t313;
t266 = -t12 * t136 + t13 * t135 + t2 * t277 - t276 * t3;
t265 = t135 * t57 - t136 * t56 + t15 * t277 - t16 * t276;
t264 = t135 * t146 - t136 * t145 - t276 * t84 + t277 * t83;
t46 = -qJ(4) * t301 + qJD(4) * t323 - t53;
t38 = -pkin(4) * t141 - t46;
t248 = Ifges(4,5) * t313;
t247 = Ifges(5,5) * t314;
t242 = pkin(5) * t255 + qJ(4);
t234 = pkin(5) * t311 + qJD(4);
t229 = Ifges(3,5) * t300;
t223 = Ifges(4,1) * t256 + t334;
t221 = Ifges(4,2) * t260 + t335;
t218 = -Ifges(5,2) * t256 - t330;
t217 = -Ifges(5,3) * t260 - t331;
t213 = -pkin(2) + t285;
t209 = t351 * t314;
t206 = (Ifges(4,1) * t260 - t335) * qJD(3);
t205 = t288 * qJD(5);
t204 = (-Ifges(4,2) * t256 + t334) * qJD(3);
t203 = t287 * qJD(5);
t201 = (-Ifges(5,2) * t260 + t331) * qJD(3);
t200 = (Ifges(5,3) * t256 - t330) * qJD(3);
t199 = (mrSges(4,1) * t256 + mrSges(4,2) * t260) * qJD(3);
t198 = (-mrSges(5,2) * t256 - mrSges(5,3) * t260) * qJD(3);
t197 = t289 * qJD(5);
t189 = (-mrSges(7,1) * t254 - mrSges(7,2) * t258) * qJD(6) * pkin(5);
t187 = t289 * t260;
t180 = pkin(5) * t318 + t225;
t179 = -qJ(4) * t313 + t290;
t172 = Ifges(6,5) * t256 - t260 * t288;
t171 = Ifges(6,6) * t256 - t260 * t287;
t170 = Ifges(6,3) * t256 + t260 * t286;
t159 = mrSges(6,1) * t313 - mrSges(6,3) * t268;
t158 = -mrSges(6,2) * t313 + mrSges(6,3) * t269;
t154 = mrSges(7,1) * t256 + mrSges(7,3) * t177;
t153 = -mrSges(7,2) * t256 + mrSges(7,3) * t176;
t152 = -mrSges(4,1) * t323 - mrSges(4,3) * t182;
t151 = mrSges(4,2) * t323 - mrSges(4,3) * t181;
t150 = mrSges(5,1) * t182 - mrSges(5,2) * t323;
t149 = mrSges(5,1) * t181 + mrSges(5,3) * t323;
t147 = -pkin(5) * t297 + (-pkin(5) * t259 - t351) * t314;
t137 = mrSges(7,1) * t277 - mrSges(7,2) * t276;
t123 = -mrSges(6,1) * t269 + mrSges(6,2) * t268;
t122 = -mrSges(5,2) * t181 - mrSges(5,3) * t182;
t120 = -mrSges(7,1) * t176 - mrSges(7,2) * t177;
t119 = -t222 * t310 + (Ifges(6,5) * t260 + t256 * t288) * qJD(3);
t118 = -t220 * t310 + (Ifges(6,6) * t260 + t256 * t287) * qJD(3);
t116 = mrSges(4,1) * t301 - mrSges(4,3) * t142;
t115 = -mrSges(4,2) * t301 - mrSges(4,3) * t141;
t113 = mrSges(5,1) * t141 - mrSges(5,3) * t301;
t111 = Ifges(4,1) * t182 - Ifges(4,4) * t181 - Ifges(4,5) * t323;
t110 = Ifges(4,4) * t182 - Ifges(4,2) * t181 - Ifges(4,6) * t323;
t109 = -Ifges(5,4) * t323 - Ifges(5,2) * t182 + Ifges(5,6) * t181;
t108 = -Ifges(5,5) * t323 - Ifges(5,6) * t182 + Ifges(5,3) * t181;
t107 = -Ifges(7,1) * t177 + Ifges(7,4) * t176 + Ifges(7,5) * t256;
t106 = -Ifges(7,4) * t177 + Ifges(7,2) * t176 + Ifges(7,6) * t256;
t105 = -Ifges(7,5) * t177 + Ifges(7,6) * t176 + Ifges(7,3) * t256;
t98 = pkin(3) * t181 + t267;
t95 = -mrSges(6,1) * t143 - mrSges(6,2) * t270;
t92 = -mrSges(5,2) * t141 - mrSges(5,3) * t142;
t91 = mrSges(4,1) * t141 + mrSges(4,2) * t142;
t90 = -Ifges(7,1) * t136 - Ifges(7,4) * t135;
t89 = -Ifges(7,4) * t136 - Ifges(7,2) * t135;
t87 = mrSges(7,1) * t135 - mrSges(7,2) * t136;
t81 = -mrSges(7,2) * t313 + mrSges(7,3) * t97;
t80 = mrSges(7,1) * t313 - mrSges(7,3) * t96;
t79 = Ifges(4,1) * t142 - Ifges(4,4) * t141 + Ifges(4,5) * t301;
t78 = Ifges(4,4) * t142 - Ifges(4,2) * t141 + Ifges(4,6) * t301;
t77 = Ifges(5,4) * t301 - Ifges(5,2) * t142 + Ifges(5,6) * t141;
t76 = Ifges(5,5) * t301 - Ifges(5,6) * t142 + Ifges(5,3) * t141;
t71 = -qJD(5) * t128 + t291;
t63 = -Ifges(6,5) * t270 + Ifges(6,6) * t143 + Ifges(6,3) * t182;
t59 = mrSges(7,1) * t182 - mrSges(7,3) * t86;
t58 = -mrSges(7,2) * t182 + mrSges(7,3) * t85;
t51 = pkin(3) * t141 + t263;
t50 = -pkin(3) * t301 - t54;
t49 = -pkin(5) * t143 + t75;
t48 = mrSges(6,1) * t142 - mrSges(6,3) * t73;
t47 = -mrSges(6,2) * t142 + mrSges(6,3) * t72;
t45 = -mrSges(7,1) * t97 + mrSges(7,2) * t96;
t43 = Ifges(7,1) * t96 + Ifges(7,4) * t97 + Ifges(7,5) * t313;
t42 = Ifges(7,4) * t96 + Ifges(7,2) * t97 + Ifges(7,6) * t313;
t39 = -mrSges(7,1) * t85 + mrSges(7,2) * t86;
t34 = Ifges(7,1) * t86 + Ifges(7,4) * t85 + Ifges(7,5) * t182;
t33 = Ifges(7,4) * t86 + Ifges(7,2) * t85 + Ifges(7,6) * t182;
t32 = Ifges(7,5) * t86 + Ifges(7,6) * t85 + Ifges(7,3) * t182;
t31 = -mrSges(6,1) * t72 + mrSges(6,2) * t73;
t28 = Ifges(6,1) * t73 + Ifges(6,4) * t72 + Ifges(6,5) * t142;
t27 = Ifges(6,4) * t73 + Ifges(6,2) * t72 + Ifges(6,6) * t142;
t23 = -pkin(5) * t72 + t38;
t18 = -mrSges(7,2) * t142 + mrSges(7,3) * t22;
t17 = mrSges(7,1) * t142 - mrSges(7,3) * t21;
t9 = -mrSges(7,1) * t22 + mrSges(7,2) * t21;
t8 = Ifges(7,1) * t21 + Ifges(7,4) * t22 + Ifges(7,5) * t142;
t7 = Ifges(7,4) * t21 + Ifges(7,2) * t22 + Ifges(7,6) * t142;
t1 = [(t12 * t3 + t13 * t2 + t23 * t49) * t358 + (t10 * t30 + t11 * t29 + t38 * t75) * t359 + 0.2e1 * m(3) * (t174 * t186 - t175 * t185) + 0.2e1 * m(4) * (t103 * t54 + t104 * t53 + t167 * t175) + 0.2e1 * m(5) * (t100 * t50 + t46 * t99 + t51 * t98) - t270 * t28 + (t175 * t257 * t356 + (t174 * t356 - t271) * t261 + ((t185 * t357 + Ifges(3,5) * t253 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t261) * t306) * t261 + (t186 * t357 - 0.2e1 * Ifges(3,6) * t253 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t257) * t306 + t364 * t182 + t363 * t181 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t365) * t323) * t257) * qJD(2)) * t252 + (t26 + t6 - t77 + t79 + 0.2e1 * t326) * t182 + (t76 - t78 + 0.2e1 * t327) * t181 + (t229 - 0.2e1 * t328 - 0.2e1 * t329) * t253 + (t111 + t63 + t32 - t109) * t142 + (-t110 + t108) * t141 + 0.2e1 * t167 * t91 + 0.2e1 * t46 * t149 + 0.2e1 * t50 * t150 + 0.2e1 * t53 * t151 + 0.2e1 * t54 * t152 + t143 * t27 + 0.2e1 * t103 * t116 + 0.2e1 * t51 * t122 + 0.2e1 * t11 * t102 + 0.2e1 * t99 * t113 + 0.2e1 * t100 * t114 + 0.2e1 * t104 * t115 + 0.2e1 * t38 * t95 + 0.2e1 * t98 * t92 + 0.2e1 * t10 * t101 + t86 * t8 + t85 * t7 + t72 * t64 + t73 * t65 + 0.2e1 * t75 * t31 + 0.2e1 * t2 * t58 + 0.2e1 * t3 * t59 + 0.2e1 * t29 * t48 + 0.2e1 * t49 * t9 + 0.2e1 * t30 * t47 + 0.2e1 * t23 * t39 + t22 * t33 + t21 * t34 + 0.2e1 * t13 * t18 + 0.2e1 * t12 * t17; t7 * t344 + t119 * t345 + t118 * t346 + t43 * t353 + t42 * t354 + t172 * t355 + t8 * t343 - t328 - t329 + (t28 * t339 - t46 * mrSges(5,1) + t53 * mrSges(4,3) + t27 * t338 - t327 - t76 / 0.2e1 + t78 / 0.2e1 + (t65 * t338 + t255 * t64 / 0.2e1) * qJD(5)) * t260 + (t50 * mrSges(5,1) - t54 * mrSges(4,3) + t326 + t79 / 0.2e1 - t77 / 0.2e1 + t26 / 0.2e1 + t6 / 0.2e1) * t256 + ((-t113 + t115) * t260 + (t114 - t116) * t256 + ((t150 - t152) * t260 + (t149 - t151) * t256) * qJD(3) + m(5) * (t100 * t313 + t50 * t256 - t46 * t260 + t314 * t99) + m(4) * (-t103 * t313 - t104 * t314 - t54 * t256 + t53 * t260)) * pkin(9) + ((-t247 / 0.2e1 - t248 / 0.2e1) * t261 + (-Ifges(3,6) + (Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1) * t260 + (Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t256) * t315) * t252 + ((-t109 / 0.2e1 + t111 / 0.2e1 + t63 / 0.2e1 + t32 / 0.2e1 + Ifges(5,4) * t295 + t100 * mrSges(5,1) - t103 * mrSges(4,3)) * t260 + (Ifges(4,6) * t295 + t108 / 0.2e1 - t110 / 0.2e1 + t99 * mrSges(5,1) - t104 * mrSges(4,3) - t272) * t256) * qJD(3) + m(7) * (t12 * t16 + t13 * t15 + t147 * t49 + t180 * t23 + t2 * t57 + t3 * t56) + (-m(4) * t175 - t91) * pkin(2) + t213 * t92 + t51 * t215 + t225 * t31 + t98 * t198 + t167 * t199 + t11 * t207 + t10 * t208 - t209 * t95 + t38 * t187 + t180 * t9 + t179 * t122 + t72 * t171 / 0.2e1 + t147 * t39 + t2 * t153 + t3 * t154 + t30 * t158 + t29 * t159 + t127 * t48 + t128 * t47 + t23 * t120 + t75 * t123 + t71 * t102 + t22 * t106 / 0.2e1 + t21 * t107 / 0.2e1 + t96 * t34 / 0.2e1 + t97 * t33 / 0.2e1 + t70 * t101 + t12 * t80 + t13 * t81 + t56 * t17 + t57 * t18 + t15 * t58 + t16 * t59 + t49 * t45 + t229 + (-t201 / 0.2e1 + t206 / 0.2e1 + t117 / 0.2e1 + t41 / 0.2e1) * t182 + (t200 / 0.2e1 - t204 / 0.2e1) * t181 + (-t218 / 0.2e1 + t223 / 0.2e1 + t170 / 0.2e1 + t105 / 0.2e1) * t142 + (t217 / 0.2e1 - t221 / 0.2e1) * t141 + m(6) * (t10 * t128 + t11 * t127 - t209 * t75 + t225 * t38 + t29 * t71 + t30 * t70) + m(5) * (t179 * t98 + t213 * t51); 0.2e1 * t213 * t198 + 0.2e1 * t225 * t123 - 0.2e1 * pkin(2) * t199 + 0.2e1 * t71 * t207 + 0.2e1 * t70 * t208 - 0.2e1 * t209 * t187 + t176 * t42 - t177 * t43 + 0.2e1 * t180 * t45 + 0.2e1 * t147 * t120 + 0.2e1 * t15 * t153 + 0.2e1 * t16 * t154 + 0.2e1 * t128 * t158 + 0.2e1 * t127 * t159 + t97 * t106 + t96 * t107 + 0.2e1 * t56 * t80 + 0.2e1 * t57 * t81 + 0.2e1 * (m(5) * t213 + t215) * t179 + (t127 * t71 + t128 * t70 - t209 * t225) * t359 + (t147 * t180 + t15 * t57 + t16 * t56) * t358 + (t117 - t201 + t206 + t41 + (t171 * t259 + t172 * t255 + t217 - t221) * qJD(3)) * t256 + (-t259 * t118 - t255 * t119 - t200 + t204 + (t171 * t255 - t172 * t259) * qJD(5) + (t105 + t170 - t218 + t223) * qJD(3)) * t260; -t205 * t345 - t203 * t346 + t21 * t347 + t22 * t348 + t34 * t349 + t33 * t350 + t90 * t353 + t89 * t354 + t222 * t355 + t72 * t340 + t8 * t341 + t7 * t342 + (-t352 * t48 - t11 * mrSges(6,3) + t28 / 0.2e1) * t259 + (-t352 * t47 - t10 * mrSges(6,3) - t27 / 0.2e1) * t255 + (t284 * mrSges(6,3) - (-m(6) * t284 + t283) * t352 + t272) * qJD(5) + t302 * t182 + t294 * t142 + (-t149 + t95) * qJD(4) + (-t113 + t31) * qJ(4) - t266 * mrSges(7,3) + t271 + t234 * t39 + t242 * t9 + t38 * t216 + t75 * t197 + t146 * t18 + t145 * t17 + t23 * t137 - pkin(3) * t114 + t49 * t87 + t83 * t58 + t84 * t59 - t53 * mrSges(4,2) + t54 * mrSges(4,1) + t50 * mrSges(5,2) - t46 * mrSges(5,3) + m(7) * (t12 * t84 + t13 * t83 + t145 * t3 + t146 * t2 + t23 * t242 + t234 * t49) + m(5) * (-pkin(3) * t50 - qJ(4) * t46 - qJD(4) * t99) + m(6) * (qJ(4) * t38 + qJD(4) * t75 - t10 * t321 - t11 * t317); t96 * t347 + t97 * t348 + t107 * t349 + t106 * t350 + t43 * t341 + t42 * t342 + t90 * t343 + t89 * t344 + ((-t171 / 0.2e1 - t128 * mrSges(6,3) - t260 * t222 / 0.2e1) * t259 + (-t172 / 0.2e1 + t260 * t340 + t127 * mrSges(6,3)) * t255 - (m(6) * t282 + t278) * t352) * qJD(5) + (-t352 * t159 - t71 * mrSges(6,3) + t119 / 0.2e1) * t259 + (-t352 * t158 - t70 * mrSges(6,3) - t118 / 0.2e1) * t255 + ((-pkin(3) * mrSges(5,1) - Ifges(5,4) + t294) * t260 + (-Ifges(4,6) + t319 / 0.2e1 - qJ(4) * mrSges(5,1) + t322 / 0.2e1) * t256 + (m(5) * t285 - t260 * mrSges(4,1) + t256 * mrSges(4,2) + t215) * pkin(9)) * qJD(3) + t302 * t256 - t265 * mrSges(7,3) + t234 * t120 + t242 * t45 - t209 * t216 + t225 * t197 + qJD(4) * t187 + t180 * t87 + t146 * t81 + t147 * t137 + t83 * t153 + t84 * t154 + t145 * t80 + qJ(4) * t123 + m(7) * (t145 * t16 + t146 * t15 + t147 * t242 + t180 * t234 + t56 * t84 + t57 * t83) + (t304 * qJD(4) - t203 * t338 - t205 * t339) * t260 + m(6) * (-qJ(4) * t209 + qJD(4) * t225 - t317 * t71 - t321 * t70) + t247 + t248; 0.2e1 * t234 * t137 + 0.2e1 * t242 * t87 - t135 * t139 - t277 * t89 - t136 * t140 - t276 * t90 + (t145 * t84 + t146 * t83 + t234 * t242) * t358 - t259 * t205 + t255 * t203 + 0.2e1 * qJ(4) * t197 + (-t319 - t322) * qJD(5) + 0.2e1 * (mrSges(5,3) + t216 + (m(5) + m(6)) * qJ(4)) * qJD(4) - t264 * t307; t135 * t58 - t136 * t59 - t276 * t17 + t277 * t18 + t255 * t47 + t259 * t48 + t283 * qJD(5) + m(7) * t266 + m(6) * (-qJD(5) * t284 + t10 * t255 + t11 * t259) + m(5) * t50 + t114; t135 * t153 - t136 * t154 + t255 * t158 + t259 * t159 + t277 * t81 - t276 * t80 + t278 * qJD(5) + t304 * t313 + m(7) * t265 + m(6) * (qJD(5) * t282 + t255 * t70 + t259 * t71); m(7) * t264 + t281 * t307; -0.2e1 * m(7) * t281; t11 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (-t12 * t309 + t13 * t308 + t2 * t254 + t258 * t3) - t59 * t309 + t258 * t17 + t58 * t308 + t254 * t18) * pkin(5) + t274 + t26; t71 * mrSges(6,1) - t70 * mrSges(6,2) + (t153 * t308 + t254 * t81 - t154 * t309 + t258 * t80 + m(7) * (t15 * t254 + t16 * t258 + t308 * t57 - t309 * t56)) * pkin(5) + t117 + t273; ((mrSges(6,2) * t352 - Ifges(6,6)) * t259 + (mrSges(6,1) * t352 - Ifges(6,5)) * t255) * qJD(5) + (m(7) * (t254 * t83 + t258 * t84 + (-t145 * t254 + t146 * t258) * qJD(6)) + t361 * mrSges(7,3)) * pkin(5) + t275; -m(7) * pkin(5) * t361 - t216 * qJD(5) + t292; 0.2e1 * t189; t274; t273; t275; t292; t189; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
