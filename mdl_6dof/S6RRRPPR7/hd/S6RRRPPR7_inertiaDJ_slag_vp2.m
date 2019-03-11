% Calculate time derivative of joint inertia matrix for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:35
% EndTime: 2019-03-09 15:56:50
% DurationCPUTime: 6.94s
% Computational Cost: add. (5081->541), mult. (11884->778), div. (0->0), fcn. (10087->8), ass. (0->219)
t253 = Ifges(5,4) + Ifges(4,5);
t198 = cos(qJ(3));
t231 = qJD(3) * t198;
t195 = sin(qJ(3));
t233 = qJD(3) * t195;
t274 = Ifges(5,6) * t233 + t231 * t253;
t273 = -Ifges(5,2) - Ifges(4,3);
t196 = sin(qJ(2));
t199 = cos(qJ(2));
t234 = qJD(2) * t199;
t223 = t195 * t234;
t202 = t196 * t231 + t223;
t247 = Ifges(5,5) * t195;
t211 = Ifges(5,1) * t198 + t247;
t249 = Ifges(4,4) * t195;
t212 = Ifges(4,1) * t198 - t249;
t272 = (t211 + t212) * qJD(3);
t235 = qJD(2) * t196;
t271 = qJ(4) * t235 - qJD(4) * t199;
t188 = t195 * qJ(4);
t207 = -t198 * pkin(3) - t188;
t194 = sin(qJ(6));
t197 = cos(qJ(6));
t192 = sin(pkin(10));
t193 = cos(pkin(10));
t205 = t192 * t195 + t193 * t198;
t123 = t205 * t196;
t160 = -pkin(2) * t199 - pkin(8) * t196 - pkin(1);
t258 = pkin(7) * t199;
t177 = t195 * t258;
t191 = t199 * pkin(3);
t84 = pkin(4) * t199 + t177 + t191 + (-qJ(5) * t196 - t160) * t198;
t178 = t198 * t258;
t113 = t195 * t160 + t178;
t103 = -qJ(4) * t199 + t113;
t243 = t195 * t196;
t93 = qJ(5) * t243 + t103;
t43 = -t192 * t93 + t193 * t84;
t25 = pkin(5) * t199 - pkin(9) * t123 + t43;
t206 = t192 * t198 - t193 * t195;
t122 = t206 * t196;
t44 = t192 * t84 + t193 * t93;
t28 = -pkin(9) * t122 + t44;
t10 = -t194 * t28 + t197 * t25;
t156 = (pkin(2) * t196 - pkin(8) * t199) * qJD(2);
t259 = pkin(7) * t195;
t226 = -pkin(3) - t259;
t228 = qJD(5) * t198;
t239 = qJD(3) * t178 + t160 * t233;
t38 = (-qJ(5) * t234 - t156) * t198 + (qJ(5) * t233 - t228 + (-pkin(4) + t226) * qJD(2)) * t196 + t239;
t240 = t195 * t156 + t160 * t231;
t242 = t196 * t198;
t39 = (-pkin(7) * qJD(2) + qJ(5) * qJD(3)) * t242 + (qJD(5) * t196 + (-pkin(7) * qJD(3) + qJ(5) * qJD(2)) * t199) * t195 + t240 + t271;
t12 = -t192 * t39 + t193 * t38;
t232 = qJD(3) * t196;
t71 = t205 * t234 + t206 * t232;
t6 = -pkin(5) * t235 - pkin(9) * t71 + t12;
t13 = t192 * t38 + t193 * t39;
t129 = t205 * qJD(3);
t70 = t129 * t196 - t206 * t234;
t9 = pkin(9) * t70 + t13;
t1 = qJD(6) * t10 + t194 * t6 + t197 * t9;
t11 = t194 * t25 + t197 * t28;
t2 = -qJD(6) * t11 - t194 * t9 + t197 * t6;
t270 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t269 = 2 * m(4);
t268 = 2 * m(6);
t267 = 2 * m(7);
t266 = -0.2e1 * pkin(1);
t265 = 0.2e1 * pkin(7);
t200 = -pkin(3) - pkin(4);
t138 = -t192 * t194 + t193 * t197;
t157 = -qJ(4) * t192 + t193 * t200;
t145 = -pkin(5) + t157;
t158 = t193 * qJ(4) + t192 * t200;
t94 = t145 * t197 - t158 * t194;
t59 = qJD(4) * t138 + qJD(6) * t94;
t257 = t59 * mrSges(7,2);
t140 = t192 * t197 + t193 * t194;
t95 = t145 * t194 + t158 * t197;
t60 = -qJD(4) * t140 - qJD(6) * t95;
t256 = t60 * mrSges(7,1);
t254 = qJD(2) / 0.2e1;
t252 = pkin(8) - qJ(5);
t67 = -t122 * t197 - t123 * t194;
t22 = qJD(6) * t67 + t194 * t70 + t197 * t71;
t68 = -t122 * t194 + t123 * t197;
t23 = -qJD(6) * t68 - t194 * t71 + t197 * t70;
t251 = Ifges(7,5) * t22 + Ifges(7,6) * t23;
t128 = t206 * qJD(3);
t91 = t194 * t206 - t197 * t205;
t45 = qJD(6) * t91 - t128 * t194 + t129 * t197;
t92 = -t194 * t205 - t197 * t206;
t46 = -qJD(6) * t92 - t128 * t197 - t129 * t194;
t250 = Ifges(7,5) * t45 + Ifges(7,6) * t46;
t248 = Ifges(4,4) * t198;
t246 = Ifges(5,5) * t198;
t245 = Ifges(4,6) * t198;
t244 = t199 * Ifges(4,6);
t125 = -t233 * t252 - t228;
t163 = t252 * t198;
t127 = qJD(3) * t163 - qJD(5) * t195;
t73 = t193 * t125 + t192 * t127;
t119 = -Ifges(5,4) * t199 + t196 * t211;
t120 = -Ifges(4,5) * t199 + t196 * t212;
t241 = t119 + t120;
t161 = t252 * t195;
t101 = t192 * t161 + t193 * t163;
t222 = t198 * t234;
t230 = qJD(4) * t198;
t238 = qJ(4) * t222 + t196 * t230;
t237 = qJ(4) * t231 + t195 * qJD(4);
t159 = -pkin(2) + t207;
t227 = t200 * t195;
t225 = t195 * t232;
t164 = -Ifges(5,3) * t198 + t247;
t165 = Ifges(4,2) * t198 + t249;
t221 = t164 / 0.2e1 - t165 / 0.2e1;
t166 = Ifges(5,1) * t195 - t246;
t167 = Ifges(4,1) * t195 + t248;
t220 = t166 / 0.2e1 + t167 / 0.2e1;
t34 = -t70 * mrSges(6,1) + t71 * mrSges(6,2);
t5 = -t23 * mrSges(7,1) + t22 * mrSges(7,2);
t14 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t85 = t128 * mrSges(6,1) + t129 * mrSges(6,2);
t72 = -t125 * t192 + t193 * t127;
t100 = t193 * t161 - t163 * t192;
t112 = t160 * t198 - t177;
t136 = t198 * pkin(4) - t159;
t219 = -pkin(7) + t227;
t209 = Ifges(5,3) * t195 + t246;
t117 = -Ifges(5,6) * t199 + t196 * t209;
t210 = -Ifges(4,2) * t195 + t248;
t118 = t196 * t210 - t244;
t218 = t117 - t118 + t244;
t217 = -t156 * t198 + t239;
t216 = -t198 * mrSges(4,1) + t195 * mrSges(4,2);
t215 = mrSges(4,1) * t195 + mrSges(4,2) * t198;
t162 = -t198 * mrSges(5,1) - t195 * mrSges(5,3);
t214 = mrSges(5,1) * t195 - mrSges(5,3) * t198;
t130 = t138 * qJD(6);
t131 = t140 * qJD(6);
t213 = -t131 * mrSges(7,1) - t130 * mrSges(7,2);
t208 = Ifges(6,5) * t129 - Ifges(6,6) * t128;
t74 = pkin(9) * t206 + t100;
t75 = -pkin(9) * t205 + t101;
t32 = -t194 * t75 + t197 * t74;
t33 = t194 * t74 + t197 * t75;
t51 = -pkin(9) * t129 + t72;
t52 = -pkin(9) * t128 + t73;
t7 = qJD(6) * t32 + t194 * t51 + t197 * t52;
t8 = -qJD(6) * t33 - t194 * t52 + t197 * t51;
t204 = t8 * mrSges(7,1) - t7 * mrSges(7,2) + t250;
t111 = qJD(3) * t227 + t237;
t176 = qJ(4) * t242;
t102 = t196 * t219 + t176;
t203 = t222 - t225;
t108 = mrSges(5,2) * t222 + (-mrSges(5,1) * qJD(2) - mrSges(5,2) * t233) * t196;
t62 = (-t198 * t235 - t199 * t233) * pkin(7) + t240;
t201 = Ifges(6,5) * t71 - t202 * Ifges(5,6) + Ifges(6,6) * t70 - t222 * t253 + t273 * t235 + t251;
t50 = (t198 * t200 - t188) * t232 + t219 * t234 + t238;
t155 = -mrSges(5,2) * t243 - mrSges(5,3) * t199;
t154 = mrSges(5,1) * t199 + mrSges(5,2) * t242;
t153 = -mrSges(4,1) * t199 - mrSges(4,3) * t242;
t152 = mrSges(4,2) * t199 - mrSges(4,3) * t243;
t149 = t210 * qJD(3);
t148 = t209 * qJD(3);
t147 = t215 * qJD(3);
t146 = t214 * qJD(3);
t133 = t214 * t196;
t126 = pkin(3) * t233 - t237;
t121 = -t176 + (pkin(3) * t195 + pkin(7)) * t196;
t110 = -mrSges(5,2) * t202 + mrSges(5,3) * t235;
t109 = -mrSges(4,2) * t235 - mrSges(4,3) * t202;
t107 = mrSges(4,1) * t235 - mrSges(4,3) * t203;
t106 = -t112 + t191;
t105 = mrSges(6,1) * t199 - mrSges(6,3) * t123;
t104 = -mrSges(6,2) * t199 - mrSges(6,3) * t122;
t99 = pkin(5) * t205 + t136;
t98 = -Ifges(6,1) * t206 - Ifges(6,4) * t205;
t97 = -Ifges(6,4) * t206 - Ifges(6,2) * t205;
t96 = mrSges(6,1) * t205 - mrSges(6,2) * t206;
t90 = mrSges(4,1) * t202 + mrSges(4,2) * t203;
t89 = mrSges(5,1) * t202 - mrSges(5,3) * t203;
t87 = Ifges(6,1) * t129 - Ifges(6,4) * t128;
t86 = Ifges(6,4) * t129 - Ifges(6,2) * t128;
t81 = -t167 * t232 + (Ifges(4,5) * t196 + t199 * t212) * qJD(2);
t80 = -t166 * t232 + (Ifges(5,4) * t196 + t199 * t211) * qJD(2);
t79 = -t165 * t232 + (Ifges(4,6) * t196 + t199 * t210) * qJD(2);
t78 = -t164 * t232 + (Ifges(5,6) * t196 + t199 * t209) * qJD(2);
t77 = pkin(5) * t128 + t111;
t76 = mrSges(6,1) * t122 + mrSges(6,2) * t123;
t66 = Ifges(6,1) * t123 - Ifges(6,4) * t122 + t199 * Ifges(6,5);
t65 = Ifges(6,4) * t123 - Ifges(6,2) * t122 + t199 * Ifges(6,6);
t64 = pkin(5) * t122 + t102;
t63 = t235 * t259 - t217;
t61 = pkin(3) * t202 + pkin(7) * t234 + qJ(4) * t225 - t238;
t58 = -mrSges(6,1) * t235 - mrSges(6,3) * t71;
t57 = mrSges(6,2) * t235 + mrSges(6,3) * t70;
t56 = mrSges(7,1) * t199 - mrSges(7,3) * t68;
t55 = -mrSges(7,2) * t199 + mrSges(7,3) * t67;
t54 = t226 * t235 + t217;
t53 = t62 + t271;
t49 = Ifges(7,1) * t92 + Ifges(7,4) * t91;
t48 = Ifges(7,4) * t92 + Ifges(7,2) * t91;
t47 = -mrSges(7,1) * t91 + mrSges(7,2) * t92;
t31 = -mrSges(7,1) * t67 + mrSges(7,2) * t68;
t30 = Ifges(6,1) * t71 + Ifges(6,4) * t70 - Ifges(6,5) * t235;
t29 = Ifges(6,4) * t71 + Ifges(6,2) * t70 - Ifges(6,6) * t235;
t27 = Ifges(7,1) * t68 + Ifges(7,4) * t67 + Ifges(7,5) * t199;
t26 = Ifges(7,4) * t68 + Ifges(7,2) * t67 + Ifges(7,6) * t199;
t24 = -pkin(5) * t70 + t50;
t18 = mrSges(7,2) * t235 + mrSges(7,3) * t23;
t17 = -mrSges(7,1) * t235 - mrSges(7,3) * t22;
t16 = Ifges(7,1) * t45 + Ifges(7,4) * t46;
t15 = Ifges(7,4) * t45 + Ifges(7,2) * t46;
t4 = Ifges(7,1) * t22 + Ifges(7,4) * t23 - Ifges(7,5) * t235;
t3 = Ifges(7,4) * t22 + Ifges(7,2) * t23 - Ifges(7,6) * t235;
t19 = [0.2e1 * t62 * t152 + 0.2e1 * t63 * t153 + 0.2e1 * t54 * t154 + 0.2e1 * t53 * t155 + 0.2e1 * t61 * t133 + 0.2e1 * t121 * t89 - t122 * t29 + t123 * t30 + 0.2e1 * t102 * t34 + 0.2e1 * t13 * t104 + 0.2e1 * t12 * t105 + 0.2e1 * t106 * t108 + 0.2e1 * t103 * t110 + 0.2e1 * t112 * t107 + 0.2e1 * t113 * t109 + 0.2e1 * t50 * t76 + t67 * t3 + t68 * t4 + t70 * t65 + t71 * t66 + 0.2e1 * t64 * t5 + 0.2e1 * t1 * t55 + 0.2e1 * t2 * t56 + 0.2e1 * t44 * t57 + 0.2e1 * t43 * t58 + t23 * t26 + t22 * t27 + 0.2e1 * t24 * t31 + 0.2e1 * t10 * t17 + 0.2e1 * t11 * t18 + (t90 * t265 + (t80 + t81) * t198 + (t78 - t79) * t195 + (t218 * t198 + (t199 * t253 - t241) * t195) * qJD(3) + (mrSges(3,1) * t266 - Ifges(6,5) * t123 + Ifges(6,6) * t122 - Ifges(7,5) * t68 - Ifges(7,6) * t67 + (-(2 * Ifges(3,4)) + t253 * t198 + (-Ifges(4,6) + Ifges(5,6)) * t195) * t196 + (pkin(7) ^ 2 * t269 + t215 * t265 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(6,3)) - (2 * Ifges(7,3)) + t273) * t199) * qJD(2)) * t196 + (t1 * t11 + t10 * t2 + t24 * t64) * t267 + (t102 * t50 + t12 * t43 + t13 * t44) * t268 + (t112 * t63 + t113 * t62) * t269 + 0.2e1 * m(5) * (t103 * t53 + t106 * t54 + t121 * t61) + ((mrSges(3,2) * t266 + 0.2e1 * Ifges(3,4) * t199 + t195 * t218 + t198 * t241) * qJD(2) + t201) * t199; -t205 * t29 / 0.2e1 + (t12 * t206 - t128 * t44 - t129 * t43 - t13 * t205) * mrSges(6,3) + (t245 * t254 - Ifges(3,6) * qJD(2) + (mrSges(3,2) * qJD(2) + t147) * pkin(7) - (-Ifges(6,5) * t206 + Ifges(7,5) * t92 - Ifges(6,6) * t205 + Ifges(7,6) * t91) * qJD(2) / 0.2e1 + (t148 / 0.2e1 - t149 / 0.2e1 - t220 * qJD(3) + t253 * t254) * t195 + (t272 / 0.2e1 - Ifges(5,6) * t254 + qJD(3) * t221) * t198) * t196 - t206 * t30 / 0.2e1 - t274 * t199 / 0.2e1 + m(5) * (t121 * t126 + t159 * t61) + t61 * t162 + t121 * t146 + t159 * t89 + t126 * t133 + t136 * t34 - t128 * t65 / 0.2e1 + t129 * t66 / 0.2e1 - t122 * t86 / 0.2e1 + t123 * t87 / 0.2e1 + t100 * t58 + t101 * t57 + t102 * t85 + t73 * t104 + t72 * t105 + t111 * t76 - pkin(2) * t90 + t91 * t3 / 0.2e1 + t92 * t4 / 0.2e1 + t50 * t96 + t70 * t97 / 0.2e1 + t71 * t98 / 0.2e1 + t99 * t5 + t77 * t31 + t67 * t15 / 0.2e1 + t68 * t16 / 0.2e1 + t64 * t14 + t7 * t55 + t8 * t56 + t45 * t27 / 0.2e1 + t46 * t26 / 0.2e1 + t24 * t47 + t23 * t48 / 0.2e1 + t22 * t49 / 0.2e1 + t33 * t18 + t32 * t17 + (t1 * t91 - t10 * t45 + t11 * t46 - t2 * t92) * mrSges(7,3) + (t250 + t208) * t199 / 0.2e1 + (-t78 / 0.2e1 + t79 / 0.2e1 + t62 * mrSges(4,3) + t53 * mrSges(5,2) + (t106 * mrSges(5,2) - t112 * mrSges(4,3) + t120 / 0.2e1 + t119 / 0.2e1) * qJD(3)) * t198 + m(6) * (t100 * t12 + t101 * t13 + t102 * t111 + t136 * t50 + t43 * t72 + t44 * t73) + m(7) * (t1 * t33 + t10 * t8 + t11 * t7 + t2 * t32 + t24 * t99 + t64 * t77) + (Ifges(3,5) + t220 * t198 + t221 * t195 + (-m(4) * pkin(2) - mrSges(3,1) + t216) * pkin(7)) * t234 + ((t109 + t110) * t198 + (-t107 + t108) * t195 + ((-t153 + t154) * t198 + (-t152 - t155) * t195) * qJD(3) + m(4) * (-t112 * t231 - t113 * t233 - t195 * t63 + t198 * t62) + m(5) * (-t103 * t233 + t106 * t231 + t195 * t54 + t198 * t53)) * pkin(8) + (t80 / 0.2e1 + t81 / 0.2e1 - t63 * mrSges(4,3) + t54 * mrSges(5,2) + (-t113 * mrSges(4,3) - t103 * mrSges(5,2) + t244 / 0.2e1 + t117 / 0.2e1 - t118 / 0.2e1) * qJD(3)) * t195; -0.2e1 * pkin(2) * t147 + 0.2e1 * t111 * t96 - t128 * t97 + t129 * t98 + 0.2e1 * t136 * t85 - t205 * t86 - t206 * t87 + 0.2e1 * t99 * t14 + 0.2e1 * t159 * t146 + t91 * t15 + t92 * t16 + t45 * t49 + t46 * t48 + 0.2e1 * t77 * t47 + (-t148 + t149) * t198 + t272 * t195 + (t100 * t72 + t101 * t73 + t111 * t136) * t268 + (t32 * t8 + t33 * t7 + t77 * t99) * t267 + ((t166 + t167) * t198 + (t164 - t165) * t195) * qJD(3) + 0.2e1 * (m(5) * t159 + t162) * t126 + 0.2e1 * (-t32 * t45 + t33 * t46 + t7 * t91 - t8 * t92) * mrSges(7,3) + 0.2e1 * (-t100 * t129 - t101 * t128 - t205 * t73 + t206 * t72) * mrSges(6,3); -t201 + (t193 * t104 - t192 * t105 + t155) * qJD(4) + m(6) * (t12 * t157 + t13 * t158 + (-t192 * t43 + t193 * t44) * qJD(4)) + t157 * t58 + t158 * t57 - pkin(3) * t108 + qJ(4) * t110 + t94 * t17 + t95 * t18 + t59 * t55 + t60 * t56 - t62 * mrSges(4,2) + t63 * mrSges(4,1) + t53 * mrSges(5,3) - t54 * mrSges(5,1) - t12 * mrSges(6,1) + t13 * mrSges(6,2) + m(5) * (-pkin(3) * t54 + qJ(4) * t53 + qJD(4) * t103) + m(7) * (t1 * t95 + t10 * t60 + t11 * t59 + t2 * t94) + t270 - Ifges(4,6) * t223 + ((Ifges(6,3) + Ifges(7,3)) * qJD(2) + (-t195 * t253 - t245) * qJD(3)) * t196; -Ifges(4,6) * t233 - t72 * mrSges(6,1) + t73 * mrSges(6,2) + m(6) * (t157 * t72 + t158 * t73 + (-t100 * t192 + t101 * t193) * qJD(4)) + m(7) * (t32 * t60 + t33 * t59 + t7 * t95 + t8 * t94) + (qJD(3) * t207 + t230) * mrSges(5,2) + (-t45 * t94 + t46 * t95 + t59 * t91 - t60 * t92) * mrSges(7,3) + (-t158 * t128 - t157 * t129 + (-t192 * t206 - t193 * t205) * qJD(4)) * mrSges(6,3) + (m(5) * t230 + (m(5) * t207 + t162 + t216) * qJD(3)) * pkin(8) - t204 - t208 + t274; (t59 * t95 + t60 * t94) * t267 + 0.2e1 * t257 - 0.2e1 * t256 + 0.2e1 * (t193 * mrSges(6,2) + t192 * mrSges(6,1) + m(6) * (-t157 * t192 + t158 * t193) + m(5) * qJ(4) + mrSges(5,3)) * qJD(4); t130 * t55 - t131 * t56 + t138 * t17 + t140 * t18 + t192 * t57 + t193 * t58 + m(7) * (t1 * t140 - t10 * t131 + t11 * t130 + t138 * t2) + m(6) * (t12 * t193 + t13 * t192) + m(5) * t54 + t108; (m(5) * pkin(8) + mrSges(5,2)) * t231 + (-t128 * t192 - t129 * t193) * mrSges(6,3) + m(7) * (t130 * t33 - t131 * t32 + t138 * t8 + t140 * t7) + m(6) * (t192 * t73 + t193 * t72) + (t130 * t91 + t131 * t92 - t138 * t45 + t140 * t46) * mrSges(7,3); m(7) * (t130 * t95 - t131 * t94 + t138 * t60 + t140 * t59) - t213; (t130 * t140 - t131 * t138) * t267; m(6) * t50 + m(7) * t24 + t34 + t5; m(6) * t111 + m(7) * t77 + t14 + t85; 0; 0; 0; -Ifges(7,3) * t235 + t251 - t270; t204; t256 - t257; t213; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
