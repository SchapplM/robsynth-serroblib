% Calculate time derivative of joint inertia matrix for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:35
% EndTime: 2019-03-09 14:08:51
% DurationCPUTime: 7.75s
% Computational Cost: add. (18282->604), mult. (46084->908), div. (0->0), fcn. (48852->12), ass. (0->254)
t253 = sin(qJ(6));
t257 = cos(qJ(6));
t258 = cos(qJ(5));
t360 = (t253 ^ 2 + t257 ^ 2) * t258;
t250 = sin(pkin(6));
t364 = 0.2e1 * t250;
t249 = sin(pkin(12));
t251 = cos(pkin(12));
t255 = sin(qJ(4));
t259 = cos(qJ(4));
t219 = -t249 * t255 + t251 * t259;
t209 = t219 * qJD(4);
t312 = t251 * t255;
t220 = t249 * t259 + t312;
t210 = t220 * qJD(4);
t254 = sin(qJ(5));
t272 = t258 * t219 - t220 * t254;
t132 = qJD(5) * t272 + t209 * t258 - t210 * t254;
t175 = t219 * t254 + t220 * t258;
t297 = qJD(6) * t257;
t270 = t253 * t132 + t175 * t297;
t298 = qJD(6) * t253;
t307 = t257 * t132;
t269 = t175 * t298 - t307;
t333 = pkin(9) + qJ(3);
t228 = t333 * t251;
t286 = t333 * t249;
t184 = t259 * t228 - t255 * t286;
t167 = pkin(10) * t219 + t184;
t296 = t259 * qJD(3);
t301 = qJD(4) * t259;
t302 = qJD(4) * t255;
t164 = -qJD(3) * t312 - t228 * t301 + (t302 * t333 - t296) * t249;
t265 = -t209 * pkin(10) + t164;
t300 = qJD(5) * t254;
t281 = t259 * t286;
t163 = -qJD(4) * t281 + t251 * t296 + (-qJD(3) * t249 - qJD(4) * t228) * t255;
t183 = -t255 * t228 - t281;
t267 = -t220 * pkin(10) + t183;
t361 = -pkin(10) * t210 + qJD(5) * t267 + t163;
t67 = -t167 * t300 + t254 * t265 + t258 * t361;
t110 = t258 * t167 + t254 * t267;
t241 = -pkin(3) * t251 - pkin(2);
t196 = -pkin(4) * t219 + t241;
t111 = -pkin(5) * t272 - pkin(11) * t175 + t196;
t78 = -t110 * t253 + t111 * t257;
t133 = qJD(5) * t175 + t209 * t254 + t258 * t210;
t336 = pkin(4) * t210;
t84 = pkin(5) * t133 - pkin(11) * t132 + t336;
t19 = qJD(6) * t78 + t253 * t84 + t257 * t67;
t79 = t110 * t257 + t111 * t253;
t20 = -qJD(6) * t79 - t253 * t67 + t257 * t84;
t363 = t19 * t257 - t20 * t253;
t260 = cos(qJ(2));
t313 = t250 * t260;
t252 = cos(pkin(6));
t256 = sin(qJ(2));
t314 = t250 * t256;
t206 = -t249 * t314 + t251 * t252;
t207 = t249 * t252 + t251 * t314;
t169 = t206 * t255 + t207 * t259;
t214 = t252 * t256 * pkin(1) + pkin(8) * t313;
t197 = qJ(3) * t252 + t214;
t198 = (-pkin(2) * t260 - qJ(3) * t256 - pkin(1)) * t250;
t161 = -t249 * t197 + t251 * t198;
t136 = -pkin(3) * t313 - t207 * pkin(9) + t161;
t162 = t251 * t197 + t249 * t198;
t150 = pkin(9) * t206 + t162;
t90 = t259 * t136 - t255 * t150;
t74 = -pkin(4) * t313 - t169 * pkin(10) + t90;
t168 = t206 * t259 - t207 * t255;
t91 = t255 * t136 + t259 * t150;
t80 = pkin(10) * t168 + t91;
t332 = t254 * t74 + t258 * t80;
t29 = -pkin(11) * t313 + t332;
t115 = t168 * t254 + t169 * t258;
t238 = pkin(8) * t314;
t337 = pkin(1) * t260;
t199 = t238 + (-pkin(2) - t337) * t252;
t170 = -t206 * pkin(3) + t199;
t121 = -t168 * pkin(4) + t170;
t273 = t258 * t168 - t169 * t254;
t54 = -pkin(5) * t273 - t115 * pkin(11) + t121;
t16 = -t253 * t29 + t257 * t54;
t303 = qJD(2) * t250;
t289 = t260 * t303;
t146 = -qJD(4) * t169 - t220 * t289;
t202 = t214 * qJD(2);
t282 = t249 * t289;
t182 = pkin(3) * t282 + t202;
t104 = -t146 * pkin(4) + t182;
t145 = qJD(4) * t168 + t219 * t289;
t64 = qJD(5) * t273 + t145 * t258 + t146 * t254;
t65 = qJD(5) * t115 + t145 * t254 - t258 * t146;
t23 = t65 * pkin(5) - t64 * pkin(11) + t104;
t290 = t256 * t303;
t299 = qJD(5) * t258;
t188 = (-qJD(3) * t256 + (pkin(2) * t256 - qJ(3) * t260) * qJD(2)) * t250;
t295 = t252 * t337;
t201 = -pkin(8) * t290 + qJD(2) * t295;
t193 = qJD(3) * t252 + t201;
t152 = t251 * t188 - t249 * t193;
t311 = t251 * t260;
t129 = (pkin(3) * t256 - pkin(9) * t311) * t303 + t152;
t153 = t249 * t188 + t251 * t193;
t148 = -pkin(9) * t282 + t153;
t50 = -qJD(4) * t91 + t259 * t129 - t148 * t255;
t33 = pkin(4) * t290 - pkin(10) * t145 + t50;
t49 = t255 * t129 + t136 * t301 + t259 * t148 - t150 * t302;
t35 = pkin(10) * t146 + t49;
t8 = t254 * t33 + t258 * t35 + t74 * t299 - t300 * t80;
t5 = pkin(11) * t290 + t8;
t2 = qJD(6) * t16 + t23 * t253 + t257 * t5;
t17 = t253 * t54 + t257 * t29;
t3 = -qJD(6) * t17 + t23 * t257 - t253 * t5;
t362 = t2 * t257 - t3 * t253;
t340 = t253 / 0.2e1;
t339 = t257 / 0.2e1;
t285 = -t298 / 0.2e1;
t9 = -qJD(5) * t332 - t254 * t35 + t258 * t33;
t359 = 2 * m(4);
t358 = 2 * m(5);
t357 = 2 * m(6);
t356 = 2 * m(7);
t355 = -2 * mrSges(3,3);
t354 = -2 * mrSges(6,3);
t68 = t167 * t299 + t254 * t361 - t258 * t265;
t353 = 0.2e1 * t68;
t109 = t254 * t167 - t258 * t267;
t352 = 0.2e1 * t109;
t271 = -t257 * t115 + t253 * t313;
t41 = qJD(6) * t271 - t253 * t64 + t257 * t290;
t351 = t41 / 0.2e1;
t102 = -t253 * t115 - t257 * t313;
t350 = t102 / 0.2e1;
t349 = t219 / 0.2e1;
t348 = t220 / 0.2e1;
t244 = Ifges(7,5) * t297;
t347 = Ifges(7,6) * t285 + t244 / 0.2e1;
t329 = Ifges(7,4) * t253;
t278 = Ifges(7,1) * t257 - t329;
t224 = t278 * qJD(6);
t346 = t224 / 0.2e1;
t345 = Ifges(7,5) * t340 + Ifges(7,6) * t339;
t328 = Ifges(7,4) * t257;
t232 = Ifges(7,1) * t253 + t328;
t343 = t232 / 0.2e1;
t342 = t251 / 0.2e1;
t341 = -t253 / 0.2e1;
t331 = Ifges(4,4) * t249;
t330 = Ifges(4,4) * t251;
t327 = Ifges(7,6) * t253;
t326 = pkin(4) * qJD(5);
t325 = t109 * t68;
t322 = t254 * mrSges(6,1);
t321 = t258 * mrSges(6,2);
t106 = -mrSges(6,1) * t313 - t115 * mrSges(6,3);
t69 = -mrSges(7,1) * t102 - mrSges(7,2) * t271;
t320 = -t106 + t69;
t319 = t109 * t254;
t318 = t175 * t253;
t317 = t175 * t257;
t309 = t253 * t258;
t229 = -mrSges(7,1) * t257 + mrSges(7,2) * t253;
t308 = t254 * t229;
t306 = t257 * t258;
t305 = Ifges(6,5) * t132 - Ifges(6,6) * t133;
t304 = Ifges(5,5) * t209 - Ifges(5,6) * t210;
t189 = t251 * mrSges(4,2) * t289 + mrSges(4,1) * t282;
t40 = qJD(6) * t102 + t253 * t290 + t257 * t64;
t12 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t65;
t294 = Ifges(6,5) * t64 - Ifges(6,6) * t65 + Ifges(6,3) * t290;
t15 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t6 = -pkin(5) * t290 - t9;
t293 = m(7) * t6 + t15;
t292 = Ifges(5,5) * t145 + Ifges(5,6) * t146 + Ifges(5,3) * t290;
t58 = mrSges(7,1) * t270 - mrSges(7,2) * t269;
t291 = m(7) * t68 + t58;
t26 = t65 * mrSges(6,1) + t64 * mrSges(6,2);
t284 = t297 / 0.2e1;
t92 = -t146 * mrSges(5,1) + t145 * mrSges(5,2);
t85 = t133 * mrSges(6,1) + t132 * mrSges(6,2);
t280 = mrSges(7,3) * t360;
t279 = mrSges(7,1) * t253 + mrSges(7,2) * t257;
t277 = -Ifges(7,2) * t253 + t328;
t75 = mrSges(7,2) * t273 + mrSges(7,3) * t102;
t76 = -mrSges(7,1) * t273 + mrSges(7,3) * t271;
t276 = -t253 * t76 + t257 * t75;
t30 = -t254 * t80 + t258 * t74;
t223 = t277 * qJD(6);
t231 = Ifges(7,2) * t257 + t329;
t268 = t257 * t223 + t253 * t224 - t231 * t298 + t232 * t297;
t44 = -Ifges(7,5) * t269 - Ifges(7,6) * t270 + Ifges(7,3) * t133;
t21 = mrSges(7,1) * t65 - mrSges(7,3) * t40;
t22 = -mrSges(7,2) * t65 + mrSges(7,3) * t41;
t264 = -t75 * t298 + m(7) * (-t16 * t297 - t17 * t298 + t362) + t257 * t22 - t253 * t21 - t76 * t297;
t134 = mrSges(7,2) * t272 - mrSges(7,3) * t318;
t135 = -mrSges(7,1) * t272 - mrSges(7,3) * t317;
t72 = mrSges(7,1) * t133 + mrSges(7,3) * t269;
t73 = -mrSges(7,2) * t133 - mrSges(7,3) * t270;
t263 = -t134 * t298 + m(7) * (-t297 * t78 - t298 * t79 + t363) + t257 * t73 - t253 * t72 - t135 * t297;
t13 = Ifges(7,4) * t40 + Ifges(7,2) * t41 + Ifges(7,6) * t65;
t14 = Ifges(7,1) * t40 + Ifges(7,4) * t41 + Ifges(7,5) * t65;
t221 = t279 * qJD(6);
t28 = pkin(5) * t313 - t30;
t52 = -Ifges(7,4) * t271 + Ifges(7,2) * t102 - Ifges(7,6) * t273;
t53 = -Ifges(7,1) * t271 + Ifges(7,4) * t102 - Ifges(7,5) * t273;
t262 = t9 * mrSges(6,1) - t8 * mrSges(6,2) - t271 * t346 - t273 * t347 + t13 * t339 + t14 * t340 + t28 * t221 + t223 * t350 + t6 * t229 + t231 * t351 + t53 * t284 + t52 * t285 + t40 * t343 + t65 * t345 + t294 + ((-t16 * t257 - t17 * t253) * qJD(6) + t362) * mrSges(7,3);
t45 = -Ifges(7,4) * t269 - Ifges(7,2) * t270 + Ifges(7,6) * t133;
t46 = -Ifges(7,1) * t269 - Ifges(7,4) * t270 + Ifges(7,5) * t133;
t98 = -Ifges(7,6) * t272 + t175 * t277;
t99 = -Ifges(7,5) * t272 + t175 * t278;
t261 = t109 * t221 + t307 * t343 + t133 * t345 - t223 * t318 / 0.2e1 + t317 * t346 - t272 * t347 + t46 * t340 + t45 * t339 + t99 * t284 - t67 * mrSges(6,2) + t305 + (t229 - mrSges(6,1)) * t68 - t270 * t231 / 0.2e1 + (t175 * t232 + t98) * t285 + ((-t253 * t79 - t257 * t78) * qJD(6) + t363) * mrSges(7,3);
t243 = -pkin(4) * t258 - pkin(5);
t242 = pkin(4) * t254 + pkin(11);
t235 = Ifges(3,5) * t289;
t213 = -t238 + t295;
t203 = t209 * mrSges(5,2);
t195 = (mrSges(4,1) * t256 - mrSges(4,3) * t311) * t303;
t194 = (-mrSges(4,3) * t249 * t260 - mrSges(4,2) * t256) * t303;
t187 = -mrSges(4,1) * t313 - t207 * mrSges(4,3);
t186 = mrSges(4,2) * t313 + t206 * mrSges(4,3);
t179 = Ifges(5,1) * t220 + Ifges(5,4) * t219;
t178 = Ifges(5,4) * t220 + Ifges(5,2) * t219;
t177 = (t256 * Ifges(4,5) + (t251 * Ifges(4,1) - t331) * t260) * t303;
t176 = (t256 * Ifges(4,6) + (-t249 * Ifges(4,2) + t330) * t260) * t303;
t173 = Ifges(5,1) * t209 - Ifges(5,4) * t210;
t172 = Ifges(5,4) * t209 - Ifges(5,2) * t210;
t171 = t210 * mrSges(5,1) + t203;
t157 = -mrSges(5,1) * t313 - t169 * mrSges(5,3);
t156 = mrSges(5,2) * t313 + t168 * mrSges(5,3);
t140 = Ifges(6,1) * t175 + Ifges(6,4) * t272;
t139 = Ifges(6,4) * t175 + Ifges(6,2) * t272;
t138 = -mrSges(6,1) * t272 + mrSges(6,2) * t175;
t124 = t279 * t175;
t118 = -mrSges(5,2) * t290 + mrSges(5,3) * t146;
t117 = mrSges(5,1) * t290 - mrSges(5,3) * t145;
t108 = Ifges(5,1) * t169 + Ifges(5,4) * t168 - Ifges(5,5) * t313;
t107 = Ifges(5,4) * t169 + Ifges(5,2) * t168 - Ifges(5,6) * t313;
t105 = mrSges(6,2) * t313 + mrSges(6,3) * t273;
t97 = -Ifges(7,3) * t272 + (Ifges(7,5) * t257 - t327) * t175;
t89 = Ifges(5,1) * t145 + Ifges(5,4) * t146 + Ifges(5,5) * t290;
t88 = Ifges(5,4) * t145 + Ifges(5,2) * t146 + Ifges(5,6) * t290;
t87 = Ifges(6,1) * t132 - Ifges(6,4) * t133;
t86 = Ifges(6,4) * t132 - Ifges(6,2) * t133;
t83 = -mrSges(6,1) * t273 + mrSges(6,2) * t115;
t82 = Ifges(6,1) * t115 + Ifges(6,4) * t273 - Ifges(6,5) * t313;
t81 = Ifges(6,4) * t115 + Ifges(6,2) * t273 - Ifges(6,6) * t313;
t57 = -mrSges(6,2) * t290 - mrSges(6,3) * t65;
t56 = mrSges(6,1) * t290 - mrSges(6,3) * t64;
t51 = -Ifges(7,5) * t271 + Ifges(7,6) * t102 - Ifges(7,3) * t273;
t25 = Ifges(6,1) * t64 - Ifges(6,4) * t65 + Ifges(6,5) * t290;
t24 = Ifges(6,4) * t64 - Ifges(6,2) * t65 + Ifges(6,6) * t290;
t1 = [0.2e1 * m(3) * (t201 * t214 - t202 * t213) + (-t81 + t51) * t65 + (t152 * t161 + t153 * t162 + t199 * t202) * t359 + (t16 * t3 + t17 * t2 + t28 * t6) * t356 + (t170 * t182 + t49 * t91 + t50 * t90) * t358 + (t104 * t121 + t30 * t9 + t332 * t8) * t357 + 0.2e1 * t332 * t57 + t252 * t235 + ((t214 * t355 + Ifges(4,5) * t207 + Ifges(5,5) * t169 + Ifges(6,5) * t115 - 0.2e1 * Ifges(3,6) * t252 + Ifges(4,6) * t206 + Ifges(5,6) * t168 + Ifges(6,6) * t273 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t256) * t364) * t256 + (t213 * t355 + t251 * (Ifges(4,1) * t207 + Ifges(4,4) * t206) - t249 * (Ifges(4,4) * t207 + Ifges(4,2) * t206) + Ifges(3,5) * t252 + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t251 + Ifges(4,6) * t249 + Ifges(3,4)) * t260) * t364 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3) - Ifges(6,3)) * t314) * t260) * t303 + t206 * t176 + 0.2e1 * t202 * (-mrSges(4,1) * t206 + mrSges(4,2) * t207) + t207 * t177 + 0.2e1 * t199 * t189 + 0.2e1 * t182 * (-mrSges(5,1) * t168 + mrSges(5,2) * t169) + 0.2e1 * t153 * t186 + 0.2e1 * t152 * t187 + 0.2e1 * t162 * t194 + 0.2e1 * t161 * t195 + t168 * t88 + t169 * t89 + 0.2e1 * t170 * t92 + 0.2e1 * t49 * t156 + 0.2e1 * t50 * t157 + t145 * t108 + t146 * t107 + 0.2e1 * t121 * t26 + t115 * t25 + 0.2e1 * t90 * t117 + 0.2e1 * t91 * t118 + 0.2e1 * t9 * t106 + t102 * t13 + 0.2e1 * t104 * t83 + 0.2e1 * t8 * t105 + t64 * t82 + 0.2e1 * t2 * t75 + 0.2e1 * t3 * t76 + 0.2e1 * t6 * t69 + t41 * t52 + t40 * t53 + 0.2e1 * t30 * t56 + 0.2e1 * t28 * t15 + 0.2e1 * t16 * t21 + 0.2e1 * t17 * t22 - t271 * t14 - (t12 - t24) * t273 - t294 * t313 - t292 * t313 + 0.2e1 * t201 * (-t252 * mrSges(3,2) + mrSges(3,3) * t313) - 0.2e1 * t202 * (mrSges(3,1) * t252 - mrSges(3,3) * t314); m(5) * (t163 * t91 + t164 * t90 + t182 * t241 + t183 * t50 + t184 * t49) + m(7) * (t109 * t6 + t16 * t20 + t17 * t19 + t2 * t79 + t28 * t68 + t3 * t78) - (-pkin(4) * t83 + t107 / 0.2e1) * t210 + (-t139 / 0.2e1 + t97 / 0.2e1) * t65 + (qJD(3) * t186 + qJ(3) * t194 + t153 * mrSges(4,3) - t202 * mrSges(4,1) + t176 / 0.2e1) * t251 + (-qJ(3) * t195 - t152 * mrSges(4,3) - qJD(3) * t187 + t202 * mrSges(4,2) + t177 / 0.2e1) * t249 + t89 * t348 + t88 * t349 + t45 * t350 + t98 * t351 + t235 + (t15 - t56) * t109 + m(6) * (t104 * t196 - t109 * t9 + t110 * t8 + t121 * t336 - t30 * t68 + t332 * t67) + (-t209 * t90 - t210 * t91 + t219 * t49 - t220 * t50) * mrSges(5,3) + (-t332 * mrSges(6,3) + t51 / 0.2e1 - t81 / 0.2e1) * t133 + m(4) * (-pkin(2) * t202 + (-t161 * t249 + t162 * t251) * qJD(3) + (-t152 * t249 + t153 * t251) * qJ(3)) + t241 * t92 + t182 * (-mrSges(5,1) * t219 + mrSges(5,2) * t220) + t209 * t108 / 0.2e1 + t196 * t26 - t201 * mrSges(3,2) - t202 * mrSges(3,1) + t145 * t179 / 0.2e1 + t183 * t117 + t184 * t118 - pkin(2) * t189 + t168 * t172 / 0.2e1 + t169 * t173 / 0.2e1 + t146 * t178 / 0.2e1 + t170 * t171 + t163 * t156 + t164 * t157 + t2 * t134 + t3 * t135 + t104 * t138 + t64 * t140 / 0.2e1 + t121 * t85 + t6 * t124 + t115 * t87 / 0.2e1 + t110 * t57 + t40 * t99 / 0.2e1 + t67 * t105 + t17 * t73 + t19 * t75 + t20 * t76 + t78 * t21 + t79 * t22 + t16 * t72 + t28 * t58 - t271 * t46 / 0.2e1 - (t44 / 0.2e1 - t86 / 0.2e1) * t273 + (((Ifges(6,5) * t175 / 0.2e1 + Ifges(6,6) * t272 / 0.2e1 + Ifges(5,5) * t348 + Ifges(5,6) * t349 - Ifges(3,6) + Ifges(4,5) * t249 / 0.2e1 + Ifges(4,6) * t342) * t256 + (-t249 * (Ifges(4,2) * t251 + t331) / 0.2e1 + (Ifges(4,1) * t249 + t330) * t342) * t260) * qJD(2) - (t304 + t305) * t260 / 0.2e1) * t250 - (-t8 * mrSges(6,3) + t12 / 0.2e1 - t24 / 0.2e1) * t272 + t320 * t68 + (t53 * t339 + t52 * t341 - t30 * mrSges(6,3) + t82 / 0.2e1) * t132 + (t14 * t339 + t13 * t341 - t9 * mrSges(6,3) + t25 / 0.2e1 + (t53 * t341 - t257 * t52 / 0.2e1) * qJD(6)) * t175; t58 * t352 + t124 * t353 + 0.2e1 * t19 * t134 + 0.2e1 * t20 * t135 + 0.2e1 * t241 * t171 + t219 * t172 + t220 * t173 + t209 * t179 + 0.2e1 * t196 * t85 + 0.2e1 * t78 * t72 + 0.2e1 * t79 * t73 - (-0.2e1 * pkin(4) * t138 + t178) * t210 - (t354 * t67 + t44 - t86) * t272 + (t110 * t354 - t139 + t97) * t133 + (mrSges(6,3) * t352 - t253 * t98 + t257 * t99 + t140) * t132 + (t110 * t67 + t196 * t336 + t325) * t357 + (t163 * t184 + t164 * t183) * t358 + (t19 * t79 + t20 * t78 + t325) * t356 + (mrSges(6,3) * t353 - t253 * t45 + t257 * t46 + t87 + (-t253 * t99 - t257 * t98) * qJD(6)) * t175 + 0.2e1 * (t163 * t219 - t164 * t220 - t183 * t209 - t184 * t210) * mrSges(5,3) + (qJ(3) * t359 + 0.2e1 * mrSges(4,3)) * (t249 ^ 2 + t251 ^ 2) * qJD(3); t257 * t21 + t253 * t22 + t276 * qJD(6) + m(7) * (t2 * t253 + t257 * t3 + (-t16 * t253 + t17 * t257) * qJD(6)) + m(6) * t104 + m(5) * t182 + m(4) * t202 + t92 + t26 + t189; m(7) * (t19 * t253 + t20 * t257 + (-t253 * t78 + t257 * t79) * qJD(6)) + t134 * t297 + t253 * t73 - t135 * t298 + t257 * t72 + t203 - (-m(6) * pkin(4) - mrSges(5,1)) * t210 + t85; 0; (m(6) * (t254 * t8 + t258 * t9) + t258 * t56 + t254 * t57 + (t320 * t254 + (t105 + t276) * t258 + m(7) * (-t16 * t309 + t17 * t306 + t254 * t28) + m(6) * (-t254 * t30 + t258 * t332)) * qJD(5)) * pkin(4) + t264 * t242 + t293 * t243 + t262 - t49 * mrSges(5,2) + t50 * mrSges(5,1) + t292; t261 + (m(6) * (t254 * t67 - t258 * t68) + (-t132 * t258 - t133 * t254) * mrSges(6,3) + ((t175 * mrSges(6,3) + t124) * t254 + (mrSges(6,3) * t272 + t134 * t257 - t135 * t253) * t258 + m(7) * (t306 * t79 - t309 * t78 + t319) + m(6) * (t110 * t258 + t319)) * qJD(5)) * pkin(4) + t291 * t243 + t263 * t242 - t163 * mrSges(5,2) + t164 * mrSges(5,1) + t304; 0; 0.2e1 * t243 * t221 + (0.2e1 * t308 + (t242 * t360 + t243 * t254) * t356 - 0.2e1 * t321 - 0.2e1 * t322 + 0.2e1 * t280) * t326 + t268; -pkin(5) * t293 + pkin(11) * t264 + t262; -pkin(5) * t291 + pkin(11) * t263 + t261; 0; (t243 - pkin(5)) * t221 + (t308 + m(7) * (-pkin(5) * t254 + pkin(11) * t360) - t321 - t322 + t280) * t326 + t268; -0.2e1 * pkin(5) * t221 + t268; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t20 - mrSges(7,2) * t19 + t44; -t221; t244 - t279 * pkin(4) * t299 + (t229 * t242 - t327) * qJD(6); t244 + (pkin(11) * t229 - t327) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
