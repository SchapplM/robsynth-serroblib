% Calculate time derivative of joint inertia matrix for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:59
% EndTime: 2019-03-09 17:04:15
% DurationCPUTime: 8.11s
% Computational Cost: add. (9808->676), mult. (25838->968), div. (0->0), fcn. (24981->10), ass. (0->289)
t382 = Ifges(7,4) + Ifges(6,5);
t242 = cos(qJ(5));
t239 = sin(qJ(5));
t320 = Ifges(7,5) * t239;
t258 = Ifges(7,1) * t242 + t320;
t322 = Ifges(6,4) * t239;
t259 = Ifges(6,1) * t242 - t322;
t366 = (t258 + t259) * qJD(5);
t381 = t366 / 0.2e1;
t380 = Ifges(4,3) + Ifges(5,3);
t319 = Ifges(7,5) * t242;
t321 = Ifges(6,4) * t242;
t365 = -t319 + t321 + (Ifges(6,1) + Ifges(7,1)) * t239;
t236 = sin(pkin(11));
t243 = cos(qJ(3));
t240 = sin(qJ(3));
t311 = cos(pkin(11));
t266 = t311 * t240;
t194 = t236 * t243 + t266;
t292 = qJD(5) * t242;
t273 = t194 * t292;
t305 = t236 * t240;
t246 = t243 * t311 - t305;
t187 = t246 * qJD(3);
t302 = t239 * t187;
t248 = t273 + t302;
t293 = qJD(5) * t239;
t274 = t194 * t293;
t301 = t242 * t187;
t247 = t274 - t301;
t238 = cos(pkin(6));
t237 = sin(pkin(6));
t241 = sin(qJ(2));
t304 = t237 * t241;
t189 = t238 * t240 + t243 * t304;
t244 = cos(qJ(2));
t297 = qJD(2) * t237;
t277 = t244 * t297;
t153 = -qJD(3) * t189 - t240 * t277;
t188 = t238 * t243 - t240 * t304;
t154 = qJD(3) * t188 + t243 * t277;
t100 = -t153 * t311 + t154 * t236;
t101 = t236 * t153 + t154 * t311;
t133 = t236 * t188 + t189 * t311;
t303 = t237 * t244;
t109 = t239 * t133 + t242 * t303;
t296 = qJD(2) * t241;
t278 = t237 * t296;
t59 = -qJD(5) * t109 + t242 * t101 + t239 * t278;
t280 = t239 * t303;
t60 = -qJD(5) * t280 + t101 * t239 + t133 * t292 - t242 * t278;
t7 = Ifges(6,5) * t59 - Ifges(6,6) * t60 + Ifges(6,3) * t100;
t8 = Ifges(7,4) * t59 + Ifges(7,2) * t100 + Ifges(7,6) * t60;
t379 = t7 + t8;
t335 = t239 / 0.2e1;
t334 = -t242 / 0.2e1;
t10 = Ifges(7,1) * t59 + Ifges(7,4) * t100 + Ifges(7,5) * t60;
t11 = Ifges(6,1) * t59 - Ifges(6,4) * t60 + Ifges(6,5) * t100;
t378 = t10 + t11;
t24 = mrSges(6,1) * t100 - mrSges(6,3) * t59;
t25 = -t100 * mrSges(7,1) + t59 * mrSges(7,2);
t377 = t25 - t24;
t26 = -mrSges(6,2) * t100 - mrSges(6,3) * t60;
t27 = -mrSges(7,2) * t60 + mrSges(7,3) * t100;
t376 = t26 + t27;
t110 = t242 * t133 - t280;
t132 = -t188 * t311 + t189 * t236;
t40 = Ifges(7,1) * t110 + Ifges(7,4) * t132 + Ifges(7,5) * t109;
t41 = Ifges(6,1) * t110 - Ifges(6,4) * t109 + Ifges(6,5) * t132;
t375 = t40 + t41;
t186 = t194 * qJD(3);
t70 = -Ifges(7,1) * t247 + Ifges(7,4) * t186 + Ifges(7,5) * t248;
t71 = -Ifges(6,1) * t247 - Ifges(6,4) * t248 + Ifges(6,5) * t186;
t374 = t70 + t71;
t72 = -mrSges(7,2) * t109 + mrSges(7,3) * t132;
t73 = -mrSges(6,2) * t132 - mrSges(6,3) * t109;
t373 = -t73 - t72;
t74 = mrSges(6,1) * t132 - mrSges(6,3) * t110;
t75 = -mrSges(7,1) * t132 + mrSges(7,2) * t110;
t372 = -t74 + t75;
t105 = mrSges(6,1) * t186 + mrSges(6,3) * t247;
t106 = -t186 * mrSges(7,1) - t247 * mrSges(7,2);
t371 = t106 - t105;
t107 = -mrSges(6,2) * t186 - mrSges(6,3) * t248;
t108 = -mrSges(7,2) * t248 + mrSges(7,3) * t186;
t370 = t107 + t108;
t115 = -Ifges(7,4) * t246 + t194 * t258;
t116 = -Ifges(6,5) * t246 + t194 * t259;
t299 = t115 + t116;
t309 = t194 * t239;
t142 = mrSges(6,2) * t246 - mrSges(6,3) * t309;
t145 = -mrSges(7,2) * t309 - mrSges(7,3) * t246;
t369 = t142 + t145;
t308 = t194 * t242;
t143 = -mrSges(6,1) * t246 - mrSges(6,3) * t308;
t144 = mrSges(7,1) * t246 + mrSges(7,2) * t308;
t368 = -t143 + t144;
t230 = -pkin(3) * t243 - pkin(2);
t141 = -pkin(4) * t246 - pkin(10) * t194 + t230;
t326 = -qJ(4) - pkin(9);
t211 = t326 * t243;
t156 = -t211 * t311 + t305 * t326;
t367 = t239 * t141 + t242 * t156;
t191 = t238 * t241 * pkin(1) + pkin(8) * t303;
t172 = pkin(9) * t238 + t191;
t173 = (-pkin(2) * t244 - pkin(9) * t241 - pkin(1)) * t237;
t121 = t243 * t172 + t240 * t173;
t223 = pkin(8) * t304;
t332 = pkin(1) * t244;
t190 = t238 * t332 - t223;
t267 = qJD(3) * t326;
t185 = qJD(4) * t243 + t240 * t267;
t245 = -qJD(4) * t240 + t243 * t267;
t127 = t185 * t311 + t236 * t245;
t295 = qJD(3) * t240;
t288 = pkin(3) * t295;
t128 = pkin(4) * t186 - pkin(10) * t187 + t288;
t32 = t242 * t127 + t239 * t128 + t141 * t292 - t156 * t293;
t33 = -qJD(5) * t367 - t127 * t239 + t128 * t242;
t363 = -t239 * t33 + t242 * t32;
t28 = qJ(6) * t186 - qJD(6) * t246 + t32;
t31 = -pkin(5) * t186 - t33;
t362 = t239 * t31 + t242 * t28;
t176 = (pkin(2) * t241 - pkin(9) * t244) * t297;
t177 = t190 * qJD(2);
t77 = -qJD(3) * t121 + t243 * t176 - t177 * t240;
t44 = pkin(3) * t278 - qJ(4) * t154 - qJD(4) * t189 + t77;
t294 = qJD(3) * t243;
t76 = -t172 * t295 + t173 * t294 + t240 * t176 + t243 * t177;
t51 = qJ(4) * t153 + qJD(4) * t188 + t76;
t18 = t236 * t44 + t311 * t51;
t16 = pkin(10) * t278 + t18;
t178 = t191 * qJD(2);
t119 = -t153 * pkin(3) + t178;
t30 = t100 * pkin(4) - t101 * pkin(10) + t119;
t102 = qJ(4) * t188 + t121;
t120 = -t240 * t172 + t243 * t173;
t88 = -pkin(3) * t303 - t189 * qJ(4) + t120;
t50 = t311 * t102 + t236 * t88;
t43 = -pkin(10) * t303 + t50;
t171 = t223 + (-pkin(2) - t332) * t238;
t135 = -t188 * pkin(3) + t171;
t63 = t132 * pkin(4) - t133 * pkin(10) + t135;
t3 = t242 * t16 + t239 * t30 + t63 * t292 - t293 * t43;
t325 = t239 * t63 + t242 * t43;
t4 = -qJD(5) * t325 - t16 * t239 + t242 * t30;
t361 = -t239 * t4 + t242 * t3;
t1 = qJ(6) * t100 + qJD(6) * t132 + t3;
t2 = -pkin(5) * t100 - t4;
t360 = t1 * t242 + t2 * t239;
t212 = -Ifges(7,3) * t242 + t320;
t215 = Ifges(6,2) * t242 + t322;
t336 = -t215 / 0.2e1;
t359 = t336 + t212 / 0.2e1;
t318 = Ifges(6,6) * t242;
t358 = Ifges(7,6) * t334 + t318 / 0.2e1 + t382 * t335;
t256 = Ifges(7,3) * t239 + t319;
t199 = t256 * qJD(5);
t257 = -Ifges(6,2) * t239 + t321;
t202 = t257 * qJD(5);
t357 = -t202 / 0.2e1 + t199 / 0.2e1;
t200 = Ifges(6,5) * t292 - Ifges(6,6) * t293;
t201 = Ifges(7,4) * t292 + Ifges(7,6) * t293;
t356 = t201 / 0.2e1 + t200 / 0.2e1;
t255 = t242 * pkin(5) + t239 * qJ(6);
t291 = qJD(6) * t242;
t355 = qJD(5) * t255 - t291;
t354 = 2 * m(5);
t353 = 2 * m(6);
t352 = 2 * m(7);
t351 = -2 * mrSges(3,3);
t350 = -2 * mrSges(5,3);
t126 = t185 * t236 - t311 * t245;
t349 = 0.2e1 * t126;
t155 = -t211 * t236 - t326 * t266;
t348 = 0.2e1 * t155;
t347 = 0.2e1 * t178;
t346 = m(5) * pkin(3);
t333 = t242 / 0.2e1;
t331 = pkin(3) * t236;
t324 = Ifges(4,4) * t240;
t323 = Ifges(4,4) * t243;
t317 = t177 * mrSges(3,2);
t316 = t178 * mrSges(3,1);
t310 = t126 * t155;
t228 = pkin(10) + t331;
t307 = t228 * t239;
t306 = t228 * t242;
t111 = -Ifges(7,6) * t246 + t194 * t256;
t114 = -Ifges(6,6) * t246 + t194 * t257;
t300 = t111 - t114;
t298 = Ifges(6,5) * t301 + Ifges(6,3) * t186;
t290 = 0.2e1 * t237;
t287 = Ifges(4,6) * t303;
t36 = Ifges(7,5) * t110 + Ifges(7,6) * t132 + Ifges(7,3) * t109;
t39 = Ifges(6,4) * t110 - Ifges(6,2) * t109 + Ifges(6,6) * t132;
t286 = t36 / 0.2e1 - t39 / 0.2e1;
t285 = t40 / 0.2e1 + t41 / 0.2e1;
t284 = mrSges(7,2) * t293;
t283 = mrSges(6,3) * t293;
t282 = mrSges(6,3) * t292;
t281 = mrSges(7,2) * t292;
t279 = t311 * pkin(3);
t276 = t228 * t293;
t275 = t228 * t292;
t270 = -t293 / 0.2e1;
t269 = t293 / 0.2e1;
t268 = t292 / 0.2e1;
t58 = t100 * mrSges(5,1) + t101 * mrSges(5,2);
t136 = t186 * mrSges(5,1) + t187 * mrSges(5,2);
t265 = Ifges(5,5) * t278;
t264 = Ifges(5,6) * t278;
t263 = Ifges(7,4) * t301 + Ifges(7,2) * t186 + t248 * Ifges(7,6);
t229 = -t279 - pkin(4);
t17 = -t236 * t51 + t311 * t44;
t49 = -t236 * t102 + t311 * t88;
t210 = -t242 * mrSges(6,1) + t239 * mrSges(6,2);
t261 = mrSges(6,1) * t239 + mrSges(6,2) * t242;
t209 = -t242 * mrSges(7,1) - t239 * mrSges(7,3);
t260 = mrSges(7,1) * t239 - mrSges(7,3) * t242;
t254 = pkin(5) * t239 - qJ(6) * t242;
t19 = -t239 * t43 + t242 * t63;
t91 = t141 * t242 - t156 * t239;
t42 = pkin(4) * t303 - t49;
t249 = Ifges(4,5) * t154 + Ifges(5,5) * t101 + Ifges(4,6) * t153 - Ifges(5,6) * t100 + t380 * t278;
t15 = -pkin(4) * t278 - t17;
t233 = Ifges(4,5) * t294;
t222 = Ifges(3,5) * t277;
t219 = Ifges(4,1) * t240 + t323;
t216 = Ifges(4,2) * t243 + t324;
t206 = (Ifges(4,1) * t243 - t324) * qJD(3);
t203 = (-Ifges(4,2) * t240 + t323) * qJD(3);
t198 = (mrSges(4,1) * t240 + mrSges(4,2) * t243) * qJD(3);
t197 = t261 * qJD(5);
t196 = t260 * qJD(5);
t192 = -t255 + t229;
t184 = -pkin(5) * t293 + qJ(6) * t292 + qJD(6) * t239;
t183 = Ifges(5,5) * t187;
t181 = Ifges(5,6) * t186;
t158 = -mrSges(4,1) * t303 - t189 * mrSges(4,3);
t157 = mrSges(4,2) * t303 + t188 * mrSges(4,3);
t149 = Ifges(5,1) * t194 + Ifges(5,4) * t246;
t148 = Ifges(5,4) * t194 + Ifges(5,2) * t246;
t147 = -mrSges(5,1) * t246 + mrSges(5,2) * t194;
t140 = t261 * t194;
t139 = t260 * t194;
t138 = Ifges(5,1) * t187 - Ifges(5,4) * t186;
t137 = Ifges(5,4) * t187 - Ifges(5,2) * t186;
t130 = mrSges(4,1) * t278 - mrSges(4,3) * t154;
t129 = -mrSges(4,2) * t278 + mrSges(4,3) * t153;
t125 = Ifges(4,1) * t189 + Ifges(4,4) * t188 - Ifges(4,5) * t303;
t124 = Ifges(4,4) * t189 + Ifges(4,2) * t188 - t287;
t118 = -mrSges(5,1) * t303 - t133 * mrSges(5,3);
t117 = mrSges(5,2) * t303 - t132 * mrSges(5,3);
t113 = -Ifges(7,2) * t246 + (Ifges(7,4) * t242 + Ifges(7,6) * t239) * t194;
t112 = -Ifges(6,3) * t246 + (Ifges(6,5) * t242 - Ifges(6,6) * t239) * t194;
t104 = -mrSges(4,1) * t153 + mrSges(4,2) * t154;
t103 = t194 * t254 + t155;
t90 = Ifges(4,1) * t154 + Ifges(4,4) * t153 + Ifges(4,5) * t278;
t89 = Ifges(4,4) * t154 + Ifges(4,2) * t153 + Ifges(4,6) * t278;
t86 = mrSges(5,1) * t278 - mrSges(5,3) * t101;
t85 = -mrSges(5,2) * t278 - mrSges(5,3) * t100;
t84 = mrSges(6,1) * t248 - mrSges(6,2) * t247;
t83 = mrSges(7,1) * t248 + mrSges(7,3) * t247;
t82 = mrSges(5,1) * t132 + mrSges(5,2) * t133;
t81 = pkin(5) * t246 - t91;
t80 = -qJ(6) * t246 + t367;
t79 = Ifges(5,1) * t133 - Ifges(5,4) * t132 - Ifges(5,5) * t303;
t78 = Ifges(5,4) * t133 - Ifges(5,2) * t132 - Ifges(5,6) * t303;
t69 = -Ifges(6,4) * t247 - Ifges(6,2) * t248 + Ifges(6,6) * t186;
t68 = -Ifges(7,4) * t274 + t263;
t67 = -Ifges(6,5) * t274 - Ifges(6,6) * t248 + t298;
t66 = -Ifges(7,5) * t247 + Ifges(7,6) * t186 + Ifges(7,3) * t248;
t65 = mrSges(6,1) * t109 + mrSges(6,2) * t110;
t64 = mrSges(7,1) * t109 - mrSges(7,3) * t110;
t52 = t254 * t187 + t194 * t355 + t126;
t46 = Ifges(5,1) * t101 - Ifges(5,4) * t100 + t265;
t45 = Ifges(5,4) * t101 - Ifges(5,2) * t100 + t264;
t38 = Ifges(7,4) * t110 + Ifges(7,2) * t132 + Ifges(7,6) * t109;
t37 = Ifges(6,5) * t110 - Ifges(6,6) * t109 + Ifges(6,3) * t132;
t23 = t109 * pkin(5) - t110 * qJ(6) + t42;
t22 = mrSges(6,1) * t60 + mrSges(6,2) * t59;
t21 = mrSges(7,1) * t60 - mrSges(7,3) * t59;
t13 = -pkin(5) * t132 - t19;
t12 = qJ(6) * t132 + t325;
t9 = Ifges(6,4) * t59 - Ifges(6,2) * t60 + Ifges(6,6) * t100;
t6 = Ifges(7,5) * t59 + Ifges(7,6) * t100 + Ifges(7,3) * t60;
t5 = t60 * pkin(5) - t59 * qJ(6) - t110 * qJD(6) + t15;
t14 = [(t15 * t42 + t19 * t4 + t3 * t325) * t353 + 0.2e1 * t325 * t26 + (mrSges(3,3) * t241 * t347 + (0.2e1 * mrSges(3,3) * t177 - t249) * t244 + ((t190 * t351 + Ifges(3,5) * t238 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t244) * t290) * t244 + (t191 * t351 + Ifges(4,5) * t189 + Ifges(5,5) * t133 - 0.2e1 * Ifges(3,6) * t238 + Ifges(4,6) * t188 - Ifges(5,6) * t132 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t241) * t290 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t380) * t303) * t241) * qJD(2)) * t237 + 0.2e1 * m(3) * (t177 * t191 - t178 * t190) + 0.2e1 * m(4) * (t120 * t77 + t121 * t76 + t171 * t178) + (t36 - t39) * t60 + (t1 * t12 + t13 * t2 + t23 * t5) * t352 + (t119 * t135 + t17 * t49 + t18 * t50) * t354 + (-mrSges(4,1) * t188 + mrSges(4,2) * t189) * t347 + t375 * t59 + t378 * t110 + (-t9 + t6) * t109 + t188 * t89 + t189 * t90 + 0.2e1 * t171 * t104 + t153 * t124 + t154 * t125 + 0.2e1 * t76 * t157 + 0.2e1 * t77 * t158 + (t37 + t38 - t78) * t100 + 0.2e1 * t121 * t129 + 0.2e1 * t120 * t130 + t133 * t46 + 0.2e1 * t135 * t58 + 0.2e1 * t18 * t117 + 0.2e1 * t17 * t118 + 0.2e1 * t119 * t82 + (-t45 + t379) * t132 + t101 * t79 + 0.2e1 * t50 * t85 + 0.2e1 * t49 * t86 + 0.2e1 * t1 * t72 + 0.2e1 * t3 * t73 + 0.2e1 * t4 * t74 + 0.2e1 * t2 * t75 + 0.2e1 * t5 * t64 + 0.2e1 * t15 * t65 + 0.2e1 * t42 * t22 + 0.2e1 * t13 * t25 + 0.2e1 * t12 * t27 + 0.2e1 * t23 * t21 + 0.2e1 * t19 * t24 + (t222 - 0.2e1 * t316 - 0.2e1 * t317) * t238; t325 * t107 + ((-t77 * t240 + t76 * t243 + (-t120 * t243 - t121 * t240) * qJD(3)) * pkin(9) - pkin(2) * t178) * m(4) - (-t264 / 0.2e1 - t18 * mrSges(5,3) + t7 / 0.2e1 + t8 / 0.2e1 - t45 / 0.2e1) * t246 + m(6) * (t126 * t42 + t15 * t155 + t19 * t33 + t3 * t367 + t32 * t325 + t4 * t91) + t367 * t26 + (t22 - t86) * t155 + (t65 - t118) * t126 + m(7) * (t1 * t80 + t103 * t5 + t12 * t28 + t13 * t31 + t2 * t81 + t23 * t52) + (t111 / 0.2e1 - t114 / 0.2e1) * t60 + (t115 / 0.2e1 + t116 / 0.2e1) * t59 - t316 - t317 + t230 * t58 + t153 * t216 / 0.2e1 + t154 * t219 / 0.2e1 + t171 * t198 + t188 * t203 / 0.2e1 + t189 * t206 / 0.2e1 + t156 * t85 + t133 * t138 / 0.2e1 + t5 * t139 + t15 * t140 + t3 * t142 + t4 * t143 + t2 * t144 + t1 * t145 + t119 * t147 + t101 * t149 / 0.2e1 + t135 * t136 + t127 * t117 + t103 * t21 - pkin(2) * t104 + t19 * t105 + t13 * t106 + t12 * t108 + t91 * t24 + t80 * t27 + t81 * t25 + t23 * t83 + t42 * t84 + t28 * t72 + t32 * t73 + t33 * t74 + t31 * t75 + t52 * t64 + (-t50 * mrSges(5,3) + t37 / 0.2e1 + t38 / 0.2e1 - t78 / 0.2e1) * t186 + (-t137 / 0.2e1 + t67 / 0.2e1 + t68 / 0.2e1) * t132 + ((-t183 / 0.2e1 + t181 / 0.2e1 - t233 / 0.2e1) * t244 + (Ifges(4,5) * t240 / 0.2e1 + Ifges(4,6) * t243 / 0.2e1 - Ifges(3,6)) * t296) * t237 + ((-pkin(9) * t158 - t120 * mrSges(4,3) + t125 / 0.2e1) * t243 + (-pkin(9) * t157 - t121 * mrSges(4,3) + pkin(3) * t82 + t287 / 0.2e1 - t124 / 0.2e1) * t240) * qJD(3) + m(5) * (t119 * t230 - t126 * t49 + t127 * t50 + t135 * t288 - t155 * t17 + t156 * t18) + (-t49 * mrSges(5,3) + t79 / 0.2e1 + t285 * t242 + t286 * t239) * t187 + (t265 / 0.2e1 - t17 * mrSges(5,3) + t46 / 0.2e1 + (t10 / 0.2e1 + t11 / 0.2e1) * t242 + (t6 / 0.2e1 - t9 / 0.2e1) * t239 + (-t239 * t285 + t242 * t286) * qJD(5)) * t194 + (t70 / 0.2e1 + t71 / 0.2e1) * t110 + (-t69 / 0.2e1 + t66 / 0.2e1) * t109 + (-t148 / 0.2e1 + t112 / 0.2e1 + t113 / 0.2e1) * t100 + (pkin(9) * t129 + t76 * mrSges(4,3) - t178 * mrSges(4,1) + t89 / 0.2e1) * t243 + (-pkin(9) * t130 - t77 * mrSges(4,3) + t178 * mrSges(4,2) + t90 / 0.2e1) * t240 + t222; t243 * t203 + t240 * t206 + 0.2e1 * t230 * t136 - 0.2e1 * pkin(2) * t198 + t84 * t348 + 0.2e1 * t52 * t139 + t140 * t349 + 0.2e1 * t32 * t142 + 0.2e1 * t33 * t143 + 0.2e1 * t31 * t144 + 0.2e1 * t28 * t145 + 0.2e1 * t91 * t105 + 0.2e1 * t81 * t106 + 0.2e1 * t367 * t107 + 0.2e1 * t80 * t108 + 0.2e1 * t103 * t83 + (t243 * t219 + (0.2e1 * pkin(3) * t147 - t216) * t240) * qJD(3) + (t127 * t156 + t230 * t288 + t310) * t354 + (t32 * t367 + t33 * t91 + t310) * t353 + (t103 * t52 + t28 * t80 + t31 * t81) * t352 - (t127 * t350 - t137 + t67 + t68) * t246 + (t156 * t350 + t112 + t113 - t148) * t186 + (mrSges(5,3) * t348 + t239 * t300 + t242 * t299 + t149) * t187 + (mrSges(5,3) * t349 + t138 + t374 * t242 + (t66 - t69) * t239 + (-t239 * t299 + t242 * t300) * qJD(5)) * t194; -t325 * t283 + t356 * t132 + t357 * t109 + t358 * t100 + t359 * t60 + m(7) * (-t184 * t23 + t192 * t5 + ((-t12 * t239 + t13 * t242) * qJD(5) + t360) * t228) + t360 * mrSges(7,2) + m(6) * (t15 * t229 + ((-t19 * t242 - t239 * t325) * qJD(5) + t361) * t228) + t361 * mrSges(6,3) + t365 * t59 / 0.2e1 + (t17 * t311 + t18 * t236) * t346 + t85 * t331 + t9 * t333 + t6 * t334 + t36 * t269 + t39 * t270 + t86 * t279 + t372 * t275 + t373 * t276 + t375 * t268 + t376 * t306 + t377 * t307 + t378 * t335 + t110 * t381 + t15 * t210 + t229 * t22 + t23 * t196 + t42 * t197 + t5 * t209 + t192 * t21 - t184 * t64 + t249 - t76 * mrSges(4,2) + t77 * mrSges(4,1) - t18 * mrSges(5,2) + t17 * mrSges(5,1) - t19 * t282 - t12 * t284 + t13 * t281; m(6) * ((-t239 * t367 - t242 * t91) * qJD(5) + t363) * t228 + t365 * (t194 * t270 + t301 / 0.2e1) - t356 * t246 - t367 * t283 + t368 * t275 - t369 * t276 + (t194 * t212 + t299) * t268 + t370 * t306 + t371 * t307 + t357 * t309 + (-mrSges(5,3) * t331 + t358) * t186 + t359 * t302 + m(7) * (-t103 * t184 + t192 * t52 + ((-t239 * t80 + t242 * t81) * qJD(5) + t362) * t228) + t362 * mrSges(7,2) + t363 * mrSges(6,3) + t183 - t181 + (m(6) * t229 - t311 * t346 - mrSges(5,1) + t210) * t126 + t273 * t336 + t69 * t333 + t66 * t334 + t114 * t270 + t111 * t269 + (-mrSges(4,1) * t294 + mrSges(4,2) * t295) * pkin(9) + t374 * t335 + t308 * t381 + t229 * t84 + t103 * t196 + t155 * t197 + t52 * t209 + t192 * t83 - t184 * t139 - Ifges(4,6) * t295 - t91 * t282 - t80 * t284 + t233 + (t236 * t346 - mrSges(5,2)) * t127 - t187 * mrSges(5,3) * t279 + t81 * t281; 0.2e1 * t192 * t196 + 0.2e1 * t197 * t229 + (-t199 + t202) * t242 + t366 * t239 + 0.2e1 * (-m(7) * t192 - t209) * t184 + (t365 * t242 + (t212 - t215) * t239) * qJD(5); -t377 * t242 + t376 * t239 + (t239 * t372 - t242 * t373) * qJD(5) + m(6) * (t239 * t3 + t242 * t4 + (-t19 * t239 + t242 * t325) * qJD(5)) + m(7) * (t1 * t239 - t2 * t242 + (t12 * t242 + t13 * t239) * qJD(5)) + m(5) * t119 + t58; m(5) * t288 - t371 * t242 + t370 * t239 + (t239 * t368 + t242 * t369) * qJD(5) + m(6) * (t239 * t32 + t242 * t33 + (-t239 * t91 + t242 * t367) * qJD(5)) + m(7) * (t239 * t28 - t242 * t31 + (t239 * t81 + t242 * t80) * qJD(5)) + t136; 0; 0; -pkin(5) * t25 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) + qJD(6) * t72 + qJ(6) * t27 + t1 * mrSges(7,3) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t379; -Ifges(6,6) * t302 - pkin(5) * t106 + m(7) * (-pkin(5) * t31 + qJ(6) * t28 + qJD(6) * t80) + qJD(6) * t145 + qJ(6) * t108 + t28 * mrSges(7,3) - t32 * mrSges(6,2) + t33 * mrSges(6,1) - t31 * mrSges(7,1) + (-t382 * t239 - t318) * t194 * qJD(5) + t263 + t298; -t355 * mrSges(7,2) + (m(7) * t291 + (-m(7) * t255 + t209 + t210) * qJD(5)) * t228 + t200 + t201; m(7) * t184 + ((-mrSges(6,2) + mrSges(7,3)) * t242 + (-mrSges(6,1) - mrSges(7,1)) * t239) * qJD(5); 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t25; m(7) * t31 + t106; (m(7) * t228 + mrSges(7,2)) * t292; m(7) * t293; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
