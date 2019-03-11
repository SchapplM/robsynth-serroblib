% Calculate time derivative of joint inertia matrix for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:24
% EndTime: 2019-03-09 10:21:37
% DurationCPUTime: 5.48s
% Computational Cost: add. (10830->557), mult. (29611->845), div. (0->0), fcn. (29839->12), ass. (0->241)
t193 = sin(pkin(6));
t303 = 0.2e1 * t193;
t302 = Ifges(5,3) + Ifges(6,3);
t191 = sin(pkin(12));
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t253 = cos(pkin(12));
t160 = t191 * t200 + t197 * t253;
t196 = sin(qJ(6));
t235 = qJD(6) * t196;
t205 = -t191 * t197 + t200 * t253;
t155 = t205 * qJD(4);
t199 = cos(qJ(6));
t242 = t199 * t155;
t206 = t160 * t235 - t242;
t276 = t196 / 0.2e1;
t275 = t199 / 0.2e1;
t221 = -t235 / 0.2e1;
t198 = sin(qJ(2));
t195 = cos(pkin(6));
t274 = pkin(1) * t195;
t180 = t198 * t274;
t201 = cos(qJ(2));
t245 = t193 * t201;
t246 = t193 * t198;
t269 = -pkin(8) - qJ(3);
t301 = (t245 * t269 - t180) * qJD(2) - qJD(3) * t246;
t154 = t160 * qJD(4);
t237 = qJD(4) * t200;
t300 = Ifges(5,5) * t237 + Ifges(6,5) * t155 - Ifges(6,6) * t154;
t234 = qJD(6) * t199;
t299 = t196 * t155 + t160 * t234;
t192 = sin(pkin(11));
t184 = pkin(2) * t192 + pkin(9);
t241 = qJ(5) + t184;
t217 = qJD(4) * t241;
t138 = qJD(5) * t200 - t197 * t217;
t204 = -qJD(5) * t197 - t200 * t217;
t105 = t138 * t253 + t191 * t204;
t238 = qJD(4) * t197;
t233 = pkin(4) * t238;
t109 = pkin(5) * t154 - pkin(10) * t155 + t233;
t194 = cos(pkin(11));
t186 = -pkin(2) * t194 - pkin(3);
t169 = -pkin(4) * t200 + t186;
t110 = -pkin(5) * t205 - pkin(10) * t160 + t169;
t158 = t241 * t200;
t218 = t241 * t197;
t113 = t158 * t253 - t191 * t218;
t71 = t110 * t199 - t113 * t196;
t38 = qJD(6) * t71 + t105 * t199 + t109 * t196;
t72 = t110 * t196 + t113 * t199;
t39 = -qJD(6) * t72 - t105 * t196 + t109 * t199;
t298 = -t196 * t39 + t199 * t38;
t244 = t194 * t201;
t144 = t192 * t246 - t193 * t244;
t145 = (t192 * t201 + t194 * t198) * t193;
t127 = t145 * t200 + t195 * t197;
t161 = (-pkin(2) * t201 - pkin(1)) * t193;
t101 = t144 * pkin(3) - t145 * pkin(9) + t161;
t181 = t201 * t274;
t222 = t269 * t198;
t128 = pkin(2) * t195 + t193 * t222 + t181;
t157 = pkin(8) * t245 + t180;
t137 = qJ(3) * t245 + t157;
t97 = t192 * t128 + t194 * t137;
t85 = pkin(9) * t195 + t97;
t58 = t200 * t101 - t197 * t85;
t41 = pkin(4) * t144 - qJ(5) * t127 + t58;
t126 = -t145 * t197 + t195 * t200;
t59 = t197 * t101 + t200 * t85;
t47 = qJ(5) * t126 + t59;
t23 = t191 * t41 + t253 * t47;
t20 = pkin(10) * t144 + t23;
t96 = t128 * t194 - t192 * t137;
t84 = -pkin(3) * t195 - t96;
t68 = -pkin(4) * t126 + t84;
t81 = -t126 * t253 + t127 * t191;
t82 = t191 * t126 + t127 * t253;
t36 = pkin(5) * t81 - pkin(10) * t82 + t68;
t10 = -t196 * t20 + t199 * t36;
t178 = qJD(2) * t181;
t118 = t178 + (qJD(2) * t222 + qJD(3) * t201) * t193;
t79 = t118 * t192 - t194 * t301;
t239 = qJD(2) * t193;
t140 = (-t192 * t198 + t244) * t239;
t94 = -qJD(4) * t127 - t140 * t197;
t53 = -pkin(4) * t94 + t79;
t95 = qJD(4) * t126 + t140 * t200;
t60 = t191 * t95 - t253 * t94;
t61 = t191 * t94 + t253 * t95;
t13 = pkin(5) * t60 - pkin(10) * t61 + t53;
t139 = qJD(2) * t145;
t228 = t198 * t239;
t216 = pkin(2) * t228;
t100 = pkin(3) * t139 - pkin(9) * t140 + t216;
t80 = t194 * t118 + t192 * t301;
t27 = -qJD(4) * t59 + t200 * t100 - t197 * t80;
t15 = pkin(4) * t139 - qJ(5) * t95 - qJD(5) * t127 + t27;
t26 = t197 * t100 + t101 * t237 + t200 * t80 - t238 * t85;
t21 = qJ(5) * t94 + qJD(5) * t126 + t26;
t6 = t191 * t15 + t253 * t21;
t4 = pkin(10) * t139 + t6;
t1 = qJD(6) * t10 + t13 * t196 + t199 * t4;
t11 = t196 * t36 + t199 * t20;
t2 = -qJD(6) * t11 + t13 * t199 - t196 * t4;
t297 = t1 * t199 - t196 * t2;
t219 = t155 * (t196 ^ 2 + t199 ^ 2);
t296 = 2 * m(6);
t295 = 2 * m(7);
t294 = -2 * mrSges(3,3);
t293 = -2 * mrSges(4,3);
t292 = -2 * mrSges(6,3);
t291 = -2 * Ifges(4,4);
t104 = t138 * t191 - t204 * t253;
t290 = 0.2e1 * t104;
t289 = 0.2e1 * t145;
t288 = 0.2e1 * t161;
t287 = m(4) * pkin(2);
t286 = m(6) * pkin(4);
t67 = t144 * t196 + t199 * t82;
t35 = -qJD(6) * t67 + t139 * t199 - t196 * t61;
t285 = t35 / 0.2e1;
t66 = t144 * t199 - t196 * t82;
t284 = t66 / 0.2e1;
t187 = Ifges(7,5) * t234;
t282 = Ifges(7,6) * t221 + t187 / 0.2e1;
t263 = Ifges(7,4) * t196;
t214 = Ifges(7,1) * t199 - t263;
t167 = t214 * qJD(6);
t281 = t167 / 0.2e1;
t280 = Ifges(7,5) * t276 + Ifges(7,6) * t275;
t262 = Ifges(7,4) * t199;
t174 = Ifges(7,1) * t196 + t262;
t278 = t174 / 0.2e1;
t277 = -t196 / 0.2e1;
t273 = pkin(4) * t191;
t270 = t80 * mrSges(4,2);
t34 = qJD(6) * t66 + t139 * t196 + t199 * t61;
t12 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t49 = mrSges(6,1) * t139 - mrSges(6,3) * t61;
t268 = t12 - t49;
t37 = -mrSges(7,1) * t66 + mrSges(7,2) * t67;
t70 = mrSges(6,1) * t144 - mrSges(6,3) * t82;
t267 = t37 - t70;
t265 = Ifges(5,4) * t197;
t264 = Ifges(5,4) * t200;
t261 = Ifges(7,6) * t196;
t260 = t139 * Ifges(6,5);
t259 = t139 * Ifges(6,6);
t258 = t144 * Ifges(5,6);
t148 = -pkin(8) * t228 + t178;
t257 = t148 * mrSges(3,2);
t149 = t157 * qJD(2);
t256 = t149 * mrSges(3,1);
t112 = t158 * t191 + t218 * t253;
t252 = t104 * t112;
t251 = t154 * t205;
t250 = t160 * t196;
t249 = t160 * t199;
t183 = pkin(10) + t273;
t248 = t183 * t196;
t247 = t183 * t199;
t236 = qJD(6) * t160;
t7 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t60;
t232 = mrSges(7,3) * t234;
t231 = mrSges(7,3) * t235;
t230 = Ifges(3,5) * t201 * t239 + Ifges(4,5) * t140 - Ifges(4,6) * t139;
t229 = t253 * pkin(4);
t227 = t183 * t235;
t226 = t183 * t234;
t31 = t60 * mrSges(6,1) + t61 * mrSges(6,2);
t220 = t234 / 0.2e1;
t114 = t154 * mrSges(6,1) + t155 * mrSges(6,2);
t170 = -mrSges(7,1) * t199 + mrSges(7,2) * t196;
t215 = mrSges(7,1) * t196 + mrSges(7,2) * t199;
t213 = -Ifges(7,2) * t196 + t262;
t212 = -t10 * t196 + t11 * t199;
t42 = -mrSges(7,2) * t81 + mrSges(7,3) * t66;
t43 = mrSges(7,1) * t81 - mrSges(7,3) * t67;
t211 = -t196 * t43 + t199 * t42;
t210 = -t197 * t27 + t200 * t26;
t209 = -t104 * t205 + t112 * t154;
t208 = Ifges(5,5) * t95 + Ifges(6,5) * t61 + Ifges(5,6) * t94 - Ifges(6,6) * t60 + t139 * t302;
t5 = t15 * t253 - t191 * t21;
t22 = -t191 * t47 + t253 * t41;
t202 = (-t10 * t199 - t11 * t196) * qJD(6) + t297;
t63 = -Ifges(7,5) * t206 - Ifges(7,6) * t299 + Ifges(7,3) * t154;
t185 = -t229 - pkin(5);
t175 = Ifges(5,1) * t197 + t264;
t173 = Ifges(5,2) * t200 + t265;
t172 = Ifges(7,2) * t199 + t263;
t168 = (Ifges(5,1) * t200 - t265) * qJD(4);
t166 = (-Ifges(5,2) * t197 + t264) * qJD(4);
t165 = t213 * qJD(6);
t163 = (mrSges(5,1) * t197 + mrSges(5,2) * t200) * qJD(4);
t162 = t215 * qJD(6);
t156 = -pkin(8) * t246 + t181;
t132 = t140 * mrSges(4,2);
t124 = Ifges(6,1) * t160 + Ifges(6,4) * t205;
t123 = Ifges(6,4) * t160 + Ifges(6,2) * t205;
t122 = -mrSges(6,1) * t205 + mrSges(6,2) * t160;
t120 = -mrSges(7,1) * t205 - mrSges(7,3) * t249;
t119 = mrSges(7,2) * t205 - mrSges(7,3) * t250;
t117 = t215 * t160;
t116 = Ifges(6,1) * t155 - Ifges(6,4) * t154;
t115 = Ifges(6,4) * t155 - Ifges(6,2) * t154;
t108 = -Ifges(7,5) * t205 + t160 * t214;
t107 = -Ifges(7,6) * t205 + t160 * t213;
t106 = -Ifges(7,3) * t205 + (Ifges(7,5) * t199 - t261) * t160;
t103 = mrSges(5,1) * t144 - mrSges(5,3) * t127;
t102 = -mrSges(5,2) * t144 + mrSges(5,3) * t126;
t93 = -mrSges(7,2) * t154 - mrSges(7,3) * t299;
t92 = mrSges(7,1) * t154 + mrSges(7,3) * t206;
t78 = -mrSges(7,1) * t299 + mrSges(7,2) * t206;
t76 = mrSges(5,1) * t139 - mrSges(5,3) * t95;
t75 = -mrSges(5,2) * t139 + mrSges(5,3) * t94;
t74 = Ifges(5,1) * t127 + Ifges(5,4) * t126 + Ifges(5,5) * t144;
t73 = Ifges(5,4) * t127 + Ifges(5,2) * t126 + t258;
t69 = -mrSges(6,2) * t144 - mrSges(6,3) * t81;
t65 = -Ifges(7,1) * t206 - Ifges(7,4) * t299 + Ifges(7,5) * t154;
t64 = -Ifges(7,4) * t206 - Ifges(7,2) * t299 + Ifges(7,6) * t154;
t62 = -mrSges(5,1) * t94 + mrSges(5,2) * t95;
t52 = mrSges(6,1) * t81 + mrSges(6,2) * t82;
t51 = Ifges(5,1) * t95 + Ifges(5,4) * t94 + t139 * Ifges(5,5);
t50 = Ifges(5,4) * t95 + Ifges(5,2) * t94 + t139 * Ifges(5,6);
t48 = -mrSges(6,2) * t139 - mrSges(6,3) * t60;
t45 = Ifges(6,1) * t82 - Ifges(6,4) * t81 + Ifges(6,5) * t144;
t44 = Ifges(6,4) * t82 - Ifges(6,2) * t81 + Ifges(6,6) * t144;
t30 = Ifges(7,1) * t67 + Ifges(7,4) * t66 + Ifges(7,5) * t81;
t29 = Ifges(7,4) * t67 + Ifges(7,2) * t66 + Ifges(7,6) * t81;
t28 = Ifges(7,5) * t67 + Ifges(7,6) * t66 + Ifges(7,3) * t81;
t25 = Ifges(6,1) * t61 - Ifges(6,4) * t60 + t260;
t24 = Ifges(6,4) * t61 - Ifges(6,2) * t60 + t259;
t19 = -t144 * pkin(5) - t22;
t17 = -mrSges(7,2) * t60 + mrSges(7,3) * t35;
t16 = mrSges(7,1) * t60 - mrSges(7,3) * t34;
t9 = Ifges(7,1) * t34 + Ifges(7,4) * t35 + Ifges(7,5) * t60;
t8 = Ifges(7,4) * t34 + Ifges(7,2) * t35 + Ifges(7,6) * t60;
t3 = -t139 * pkin(5) - t5;
t14 = [0.2e1 * t6 * t69 + 0.2e1 * t5 * t70 + 0.2e1 * t59 * t75 + 0.2e1 * t58 * t76 + t66 * t8 + t67 * t9 + 0.2e1 * t68 * t31 + t61 * t45 + 0.2e1 * t22 * t49 + 0.2e1 * t53 * t52 + 0.2e1 * t1 * t42 + 0.2e1 * t2 * t43 + 0.2e1 * t23 * t48 + 0.2e1 * m(3) * (t148 * t157 - t149 * t156) + t35 * t29 + 0.2e1 * t3 * t37 + t34 * t30 + t79 * mrSges(4,3) * t289 + 0.2e1 * t19 * t12 + 0.2e1 * t10 * t16 + 0.2e1 * t11 * t17 + 0.2e1 * m(4) * (-t79 * t96 + t80 * t97) + t132 * t288 + (t1 * t11 + t10 * t2 + t19 * t3) * t295 + (t22 * t5 + t23 * t6 + t53 * t68) * t296 + (t7 - t24) * t81 + (-0.2e1 * t79 * mrSges(4,1) + t230 - 0.2e1 * t256 - 0.2e1 * t257 - 0.2e1 * t270) * t195 + 0.2e1 * m(5) * (t26 * t59 + t27 * t58 + t79 * t84) + (0.2e1 * (t148 * t201 + t149 * t198) * mrSges(3,3) + ((t156 * t294 + Ifges(3,5) * t195 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t201) * t303) * t201 + (0.2e1 * pkin(2) * (mrSges(4,1) * t144 + mrSges(4,2) * t145) + t157 * t294 - 0.2e1 * Ifges(3,6) * t195 + t287 * t288 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t198 + (Ifges(3,1) - Ifges(3,2)) * t201) * t303) * t198) * qJD(2)) * t193 + (mrSges(4,1) * t288 + t97 * t293 + t145 * t291 + Ifges(5,5) * t127 + Ifges(6,5) * t82 - Ifges(4,6) * t195 + Ifges(5,6) * t126 - Ifges(6,6) * t81 + ((2 * Ifges(4,2)) + t302) * t144) * t139 + (t28 - t44) * t60 + (Ifges(4,1) * t289 + Ifges(4,5) * t195 + t144 * t291 + t293 * t96) * t140 + (t80 * t293 + t208) * t144 + t82 * t25 + 0.2e1 * t84 * t62 + t94 * t73 + t95 * t74 + 0.2e1 * t26 * t102 + 0.2e1 * t27 * t103 + t126 * t50 + t127 * t51 + 0.2e1 * t79 * (-mrSges(5,1) * t126 + mrSges(5,2) * t127); t71 * t16 + t72 * t17 - t19 * t78 + t67 * t65 / 0.2e1 + t38 * t42 + t39 * t43 + t230 + (-t23 * mrSges(6,3) + t28 / 0.2e1 - t44 / 0.2e1) * t154 + (-t197 * t76 + t200 * t75) * t184 + t64 * t284 + t107 * t285 + (t192 * t80 - t194 * t79) * t287 + (-mrSges(5,1) * t200 + mrSges(5,2) * t197 - mrSges(4,1)) * t79 + m(7) * (t1 * t72 + t10 * t39 + t104 * t19 + t11 * t38 + t112 * t3 + t2 * t71) + (t63 / 0.2e1 - t115 / 0.2e1) * t81 + (t106 / 0.2e1 - t123 / 0.2e1) * t60 + t300 * t144 / 0.2e1 + (-t22 * mrSges(6,3) + t29 * t277 + t30 * t275 + t45 / 0.2e1) * t155 + (-t5 * mrSges(6,3) + t8 * t277 + t9 * t275 + t25 / 0.2e1 + t260 / 0.2e1 + (t30 * t277 - t199 * t29 / 0.2e1) * qJD(6)) * t160 + t267 * t104 + t268 * t112 + t10 * t92 + t11 * t93 + t105 * t69 + t34 * t108 / 0.2e1 + t113 * t48 + t68 * t114 + t82 * t116 / 0.2e1 + t3 * t117 + t1 * t119 + t2 * t120 - t256 - t257 + t53 * t122 + t61 * t124 / 0.2e1 + t84 * t163 + t126 * t166 / 0.2e1 + t127 * t168 / 0.2e1 + t169 * t31 + t94 * t173 / 0.2e1 + t95 * t175 / 0.2e1 + t186 * t62 + (-t139 * t192 - t140 * t194) * pkin(2) * mrSges(4,3) + t210 * mrSges(5,3) + t197 * t51 / 0.2e1 + (((-t59 * t197 - t58 * t200) * qJD(4) + t210) * t184 + t186 * t79) * m(5) + t200 * t50 / 0.2e1 + t139 * (Ifges(5,5) * t197 + Ifges(5,6) * t200) / 0.2e1 - t270 - (-t6 * mrSges(6,3) + t7 / 0.2e1 - t24 / 0.2e1 - t259 / 0.2e1) * t205 - Ifges(3,6) * t228 + m(6) * (-t104 * t22 + t105 * t23 - t112 * t5 + t113 * t6 + t169 * t53 + t233 * t68) + ((-t58 * mrSges(5,3) - t184 * t103 + t74 / 0.2e1) * t200 + (pkin(4) * t52 - t59 * mrSges(5,3) - t184 * t102 - t73 / 0.2e1 - t258 / 0.2e1) * t197) * qJD(4); t117 * t290 - 0.2e1 * t112 * t78 + 0.2e1 * t169 * t114 + 0.2e1 * t38 * t119 + 0.2e1 * t39 * t120 + 0.2e1 * t186 * t163 + t200 * t166 + t197 * t168 + 0.2e1 * t71 * t92 + 0.2e1 * t72 * t93 + (t200 * t175 + (0.2e1 * pkin(4) * t122 - t173) * t197) * qJD(4) + (t38 * t72 + t39 * t71 + t252) * t295 + (t105 * t113 + t169 * t233 + t252) * t296 - (t105 * t292 - t115 + t63) * t205 + (t113 * t292 + t106 - t123) * t154 + (0.2e1 * mrSges(6,3) * t112 - t107 * t196 + t199 * t108 + t124) * t155 + (mrSges(6,3) * t290 - t196 * t64 + t199 * t65 + t116 + (-t107 * t199 - t108 * t196) * qJD(6)) * t160; m(4) * t216 + t139 * mrSges(4,1) + t197 * t75 + t200 * t76 + t132 - t268 * t205 + t267 * t154 + (t102 * t200 - t103 * t197) * qJD(4) + (t211 + t69) * t155 + (-t196 * t16 + t199 * t17 + t48 + (-t196 * t42 - t199 * t43) * qJD(6)) * t160 + m(7) * (t154 * t19 + t155 * t212 + t160 * t202 - t205 * t3) + m(6) * (-t154 * t22 + t155 * t23 + t160 * t6 + t205 * t5) + m(5) * (t197 * t26 + t200 * t27 + (-t197 * t58 + t200 * t59) * qJD(4)); t154 * t117 + t205 * t78 + m(7) * t209 + m(6) * (t105 * t160 + t113 * t155 + t209) + (m(7) * (t155 * t72 + t160 * t38 - t236 * t71) + t155 * t119 + t160 * t93 - t120 * t236) * t199 + (m(7) * (-t155 * t71 - t160 * t39 - t236 * t72) - t119 * t236 - t155 * t120 - t160 * t92) * t196; 0.2e1 * m(6) * (t155 * t160 - t251) + 0.2e1 * m(7) * (t160 * t219 - t251); t30 * t220 + t29 * t221 + m(7) * (t183 * t202 + t185 * t3) + t27 * mrSges(5,1) + t208 - t26 * mrSges(5,2) - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t165 * t284 + t172 * t285 + (t191 * t6 + t253 * t5) * t286 + t8 * t275 + t9 * t276 + t34 * t278 + t60 * t280 + t67 * t281 + t81 * t282 + t48 * t273 + t17 * t247 + t49 * t229 + t297 * mrSges(7,3) + t19 * t162 + t3 * t170 + t185 * t12 - t43 * t226 - t42 * t227 - t11 * t231 - t10 * t232 - t16 * t248; t108 * t220 + t300 + t64 * t275 + t65 * t276 + t242 * t278 + t154 * t280 + t249 * t281 + t93 * t247 + m(7) * ((-t196 * t72 - t199 * t71) * qJD(6) + t298) * t183 + (t160 * t174 + t107) * t221 + (-mrSges(5,1) * t237 + mrSges(5,2) * t238) * t184 + (t191 * t286 - mrSges(6,2)) * t105 + (m(7) * t185 - t253 * t286 - mrSges(6,1) + t170) * t104 + (-t154 * t273 - t155 * t229) * mrSges(6,3) + t298 * mrSges(7,3) - t299 * t172 / 0.2e1 + t112 * t162 - t185 * t78 - t205 * t282 - t120 * t226 - t119 * t227 - t72 * t231 - t71 * t232 - Ifges(5,6) * t238 - t92 * t248 - t165 * t250 / 0.2e1; -mrSges(5,2) * t237 - mrSges(5,1) * t238 + (-t154 * t253 + t155 * t191) * t286 + t154 * t170 - t205 * t162 + m(7) * (t185 * t154 + t183 * t219) - t114 + mrSges(7,3) * t219; 0.2e1 * t162 * t185 + t165 * t199 + t167 * t196 + (-t172 * t196 + t174 * t199) * qJD(6); t199 * t16 + t196 * t17 + t211 * qJD(6) + m(7) * (qJD(6) * t212 + t1 * t196 + t199 * t2) + m(6) * t53 + t31; m(7) * (t196 * t38 + t199 * t39 + (-t196 * t71 + t199 * t72) * qJD(6)) + t119 * t234 + t196 * t93 - t120 * t235 + t199 * t92 + m(6) * t233 + t114; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t39 - mrSges(7,2) * t38 + t63; t78; t187 + (t170 * t183 - t261) * qJD(6); -t162; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
