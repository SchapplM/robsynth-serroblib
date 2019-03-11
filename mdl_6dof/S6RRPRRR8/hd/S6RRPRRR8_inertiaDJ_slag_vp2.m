% Calculate time derivative of joint inertia matrix for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:23
% EndTime: 2019-03-09 14:02:34
% DurationCPUTime: 4.37s
% Computational Cost: add. (13573->514), mult. (31613->775), div. (0->0), fcn. (31051->10), ass. (0->213)
t212 = sin(pkin(11));
t259 = pkin(8) + qJ(3);
t194 = t259 * t212;
t213 = cos(pkin(11));
t195 = t259 * t213;
t216 = sin(qJ(4));
t220 = cos(qJ(4));
t156 = -t194 * t220 - t195 * t216;
t189 = t212 * t220 + t213 * t216;
t131 = -pkin(9) * t189 + t156;
t157 = -t194 * t216 + t195 * t220;
t226 = t212 * t216 - t213 * t220;
t132 = -pkin(9) * t226 + t157;
t215 = sin(qJ(5));
t219 = cos(qJ(5));
t88 = t131 * t215 + t132 * t219;
t217 = sin(qJ(2));
t221 = cos(qJ(2));
t193 = -pkin(2) * t221 - qJ(3) * t217 - pkin(1);
t184 = t213 * t193;
t250 = t213 * t217;
t148 = -pkin(8) * t250 + t184 + (-pkin(7) * t212 - pkin(3)) * t221;
t249 = t213 * t221;
t163 = pkin(7) * t249 + t193 * t212;
t252 = t212 * t217;
t155 = -pkin(8) * t252 + t163;
t103 = t148 * t216 + t155 * t220;
t176 = t226 * qJD(4);
t244 = qJD(2) * t217;
t167 = t189 * t217;
t168 = t226 * t217;
t125 = -t167 * t219 + t168 * t215;
t177 = t189 * qJD(4);
t243 = qJD(2) * t221;
t133 = -t177 * t217 - t226 * t243;
t134 = t176 * t217 - t189 * t243;
t68 = qJD(5) * t125 + t133 * t219 + t134 * t215;
t126 = -t167 * t215 - t168 * t219;
t69 = -qJD(5) * t126 - t133 * t215 + t134 * t219;
t276 = -Ifges(6,5) * t68 - Ifges(6,6) * t69 - Ifges(6,3) * t244;
t275 = -Ifges(5,5) * t133 - Ifges(5,6) * t134 - Ifges(5,3) * t244;
t274 = 2 * m(4);
t273 = 2 * m(5);
t272 = 2 * m(6);
t271 = 2 * m(7);
t270 = -0.2e1 * pkin(1);
t269 = 0.2e1 * pkin(7);
t146 = -t189 * t215 - t219 * t226;
t147 = t189 * t219 - t215 * t226;
t214 = sin(qJ(6));
t218 = cos(qJ(6));
t97 = t146 * t218 - t147 * t214;
t268 = t97 / 0.2e1;
t98 = t146 * t214 + t147 * t218;
t267 = t98 / 0.2e1;
t266 = t146 / 0.2e1;
t265 = t147 / 0.2e1;
t264 = -t226 / 0.2e1;
t263 = t189 / 0.2e1;
t262 = t213 / 0.2e1;
t260 = pkin(4) * t177;
t95 = qJD(5) * t146 - t176 * t219 - t177 * t215;
t96 = -qJD(5) * t147 + t176 * t215 - t177 * t219;
t39 = qJD(6) * t97 + t214 * t96 + t218 * t95;
t40 = -qJD(6) * t98 - t214 * t95 + t218 * t96;
t258 = Ifges(7,5) * t39 + Ifges(7,6) * t40;
t246 = t216 * t155;
t102 = t148 * t220 - t246;
t82 = -pkin(4) * t221 + pkin(9) * t168 + t102;
t89 = -pkin(9) * t167 + t103;
t54 = t215 * t82 + t219 * t89;
t257 = Ifges(6,5) * t95 + Ifges(6,6) * t96;
t256 = mrSges(4,2) * t213;
t255 = Ifges(4,4) * t212;
t254 = Ifges(4,4) * t213;
t204 = pkin(4) * t219 + pkin(5);
t236 = qJD(6) * t218;
t237 = qJD(6) * t214;
t248 = t214 * t215;
t144 = t204 * t236 + (-t215 * t237 + (t218 * t219 - t248) * qJD(5)) * pkin(4);
t253 = t144 * mrSges(7,2);
t251 = t212 * t221;
t247 = t215 * t218;
t245 = -Ifges(5,5) * t176 - Ifges(5,6) * t177;
t233 = t212 * t243;
t169 = mrSges(4,1) * t233 + t243 * t256;
t173 = -t217 * qJD(3) + (pkin(2) * t217 - qJ(3) * t221) * qJD(2);
t234 = pkin(7) * t244;
t153 = t173 * t213 + t212 * t234;
t208 = pkin(7) * t243;
t180 = pkin(3) * t233 + t208;
t192 = pkin(3) * t252 + pkin(7) * t217;
t242 = qJD(3) * t212;
t241 = qJD(3) * t213;
t240 = qJD(4) * t220;
t239 = qJD(5) * t215;
t238 = qJD(5) * t219;
t79 = t125 * t218 - t126 * t214;
t25 = qJD(6) * t79 + t214 * t69 + t218 * t68;
t80 = t125 * t214 + t126 * t218;
t26 = -qJD(6) * t80 - t214 * t68 + t218 * t69;
t235 = -Ifges(7,5) * t25 - Ifges(7,6) * t26 - Ifges(7,3) * t244;
t203 = -pkin(3) * t213 - pkin(2);
t31 = -t69 * mrSges(6,1) + mrSges(6,2) * t68;
t57 = -t96 * mrSges(6,1) + mrSges(6,2) * t95;
t9 = -t26 * mrSges(7,1) + mrSges(7,2) * t25;
t17 = -t40 * mrSges(7,1) + mrSges(7,2) * t39;
t53 = -t215 * t89 + t219 * t82;
t90 = -t134 * mrSges(5,1) + mrSges(5,2) * t133;
t145 = -t204 * t237 + (-t215 * t236 + (-t214 * t219 - t247) * qJD(5)) * pkin(4);
t138 = t145 * mrSges(7,1);
t232 = t138 - t253;
t87 = t131 * t219 - t132 * t215;
t109 = -pkin(4) * t134 + t180;
t150 = pkin(4) * t167 + t192;
t118 = -t194 * t240 + t220 * t241 + (-qJD(4) * t195 - t242) * t216;
t107 = -pkin(9) * t177 + t118;
t119 = -qJD(3) * t189 - qJD(4) * t157;
t108 = pkin(9) * t176 + t119;
t42 = t107 * t219 + t108 * t215 + t131 * t238 - t132 * t239;
t27 = pkin(10) * t96 + t42;
t43 = -qJD(5) * t88 - t107 * t215 + t108 * t219;
t28 = -pkin(10) * t95 + t43;
t70 = -pkin(10) * t147 + t87;
t71 = pkin(10) * t146 + t88;
t32 = -t214 * t71 + t218 * t70;
t5 = qJD(6) * t32 + t214 * t28 + t218 * t27;
t33 = t214 * t70 + t218 * t71;
t6 = -qJD(6) * t33 - t214 * t27 + t218 * t28;
t230 = mrSges(7,1) * t6 - t5 * mrSges(7,2) + t258;
t229 = Ifges(4,1) * t213 - t255;
t228 = -t212 * Ifges(4,2) + t254;
t227 = -Ifges(4,5) * t213 + Ifges(4,6) * t212;
t34 = -pkin(5) * t221 - pkin(10) * t126 + t53;
t38 = pkin(10) * t125 + t54;
t15 = -t214 * t38 + t218 * t34;
t16 = t214 * t34 + t218 * t38;
t164 = pkin(4) * t226 + t203;
t127 = (pkin(3) * t217 - pkin(8) * t249) * qJD(2) + t153;
t165 = t212 * t173;
t136 = t165 + (-pkin(7) * t250 - pkin(8) * t251) * qJD(2);
t56 = -qJD(4) * t103 + t127 * t220 - t136 * t216;
t46 = pkin(4) * t244 - pkin(9) * t133 + t56;
t55 = -qJD(4) * t246 + t127 * t216 + t136 * t220 + t148 * t240;
t48 = pkin(9) * t134 + t55;
t14 = -qJD(5) * t54 - t215 * t48 + t219 * t46;
t10 = pkin(5) * t244 - pkin(10) * t68 + t14;
t13 = t215 * t46 + t219 * t48 + t238 * t82 - t239 * t89;
t11 = pkin(10) * t69 + t13;
t2 = qJD(6) * t15 + t10 * t214 + t11 * t218;
t3 = -qJD(6) * t16 + t10 * t218 - t11 * t214;
t225 = mrSges(7,1) * t3 - t2 * mrSges(7,2) - t235;
t224 = (-mrSges(6,1) * t215 - mrSges(6,2) * t219) * qJD(5) * pkin(4);
t223 = mrSges(6,1) * t43 - t42 * mrSges(6,2) + t230 + t257;
t222 = mrSges(6,1) * t14 - t13 * mrSges(6,2) + t225 - t276;
t191 = -mrSges(4,1) * t221 - mrSges(4,3) * t250;
t190 = mrSges(4,2) * t221 - mrSges(4,3) * t252;
t182 = (-mrSges(7,1) * t214 - mrSges(7,2) * t218) * qJD(6) * pkin(5);
t179 = (mrSges(4,1) * t217 - mrSges(4,3) * t249) * qJD(2);
t178 = (-mrSges(4,2) * t217 - mrSges(4,3) * t251) * qJD(2);
t175 = pkin(4) * t247 + t204 * t214;
t174 = -pkin(4) * t248 + t204 * t218;
t170 = t176 * mrSges(5,2);
t162 = -pkin(7) * t251 + t184;
t161 = (t217 * Ifges(4,5) + t221 * t229) * qJD(2);
t160 = (t217 * Ifges(4,6) + t221 * t228) * qJD(2);
t159 = -mrSges(5,1) * t221 + mrSges(5,3) * t168;
t158 = mrSges(5,2) * t221 - mrSges(5,3) * t167;
t154 = -t213 * t234 + t165;
t152 = Ifges(5,1) * t189 - Ifges(5,4) * t226;
t151 = Ifges(5,4) * t189 - Ifges(5,2) * t226;
t143 = -Ifges(5,1) * t176 - Ifges(5,4) * t177;
t142 = -Ifges(5,4) * t176 - Ifges(5,2) * t177;
t141 = t177 * mrSges(5,1) - t170;
t124 = -Ifges(5,1) * t168 - Ifges(5,4) * t167 - Ifges(5,5) * t221;
t123 = -Ifges(5,4) * t168 - Ifges(5,2) * t167 - Ifges(5,6) * t221;
t114 = -mrSges(5,2) * t244 + mrSges(5,3) * t134;
t113 = mrSges(5,1) * t244 - mrSges(5,3) * t133;
t112 = -mrSges(6,1) * t221 - mrSges(6,3) * t126;
t111 = mrSges(6,2) * t221 + mrSges(6,3) * t125;
t110 = -pkin(5) * t146 + t164;
t101 = Ifges(6,1) * t147 + Ifges(6,4) * t146;
t100 = Ifges(6,4) * t147 + Ifges(6,2) * t146;
t99 = -mrSges(6,1) * t146 + mrSges(6,2) * t147;
t91 = -pkin(5) * t125 + t150;
t85 = -mrSges(6,1) * t125 + mrSges(6,2) * t126;
t84 = Ifges(5,1) * t133 + Ifges(5,4) * t134 + Ifges(5,5) * t244;
t83 = Ifges(5,4) * t133 + Ifges(5,2) * t134 + Ifges(5,6) * t244;
t81 = -pkin(5) * t96 + t260;
t76 = Ifges(6,1) * t126 + Ifges(6,4) * t125 - Ifges(6,5) * t221;
t75 = Ifges(6,4) * t126 + Ifges(6,2) * t125 - Ifges(6,6) * t221;
t73 = -mrSges(7,1) * t221 - mrSges(7,3) * t80;
t72 = mrSges(7,2) * t221 + mrSges(7,3) * t79;
t64 = -mrSges(6,2) * t244 + mrSges(6,3) * t69;
t63 = mrSges(6,1) * t244 - mrSges(6,3) * t68;
t62 = Ifges(7,1) * t98 + Ifges(7,4) * t97;
t61 = Ifges(7,4) * t98 + Ifges(7,2) * t97;
t60 = -mrSges(7,1) * t97 + mrSges(7,2) * t98;
t59 = Ifges(6,1) * t95 + Ifges(6,4) * t96;
t58 = Ifges(6,4) * t95 + Ifges(6,2) * t96;
t52 = -mrSges(7,1) * t79 + mrSges(7,2) * t80;
t51 = -pkin(5) * t69 + t109;
t50 = Ifges(7,1) * t80 + Ifges(7,4) * t79 - Ifges(7,5) * t221;
t49 = Ifges(7,4) * t80 + Ifges(7,2) * t79 - Ifges(7,6) * t221;
t30 = Ifges(6,1) * t68 + Ifges(6,4) * t69 + Ifges(6,5) * t244;
t29 = Ifges(6,4) * t68 + Ifges(6,2) * t69 + Ifges(6,6) * t244;
t21 = -mrSges(7,2) * t244 + mrSges(7,3) * t26;
t20 = mrSges(7,1) * t244 - mrSges(7,3) * t25;
t19 = Ifges(7,1) * t39 + Ifges(7,4) * t40;
t18 = Ifges(7,4) * t39 + Ifges(7,2) * t40;
t8 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t244;
t7 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t244;
t1 = [0.2e1 * t180 * (mrSges(5,1) * t167 - mrSges(5,2) * t168) + (-t212 * t160 + t213 * t161 + t169 * t269) * t217 + (t15 * t3 + t16 * t2 + t51 * t91) * t271 + (t109 * t150 + t13 * t54 + t14 * t53) * t272 + (t102 * t56 + t103 * t55 + t180 * t192) * t273 + (t153 * t162 + t154 * t163) * t274 + 0.2e1 * t154 * t190 + 0.2e1 * t153 * t191 + 0.2e1 * t192 * t90 + 0.2e1 * t163 * t178 + 0.2e1 * t162 * t179 - t167 * t83 - t168 * t84 + 0.2e1 * t55 * t158 + 0.2e1 * t56 * t159 + 0.2e1 * t150 * t31 + t133 * t124 + t134 * t123 + t125 * t29 + t126 * t30 + 0.2e1 * t13 * t111 + 0.2e1 * t14 * t112 + 0.2e1 * t102 * t113 + 0.2e1 * t103 * t114 + 0.2e1 * t109 * t85 + 0.2e1 * t91 * t9 + t79 * t7 + t80 * t8 + 0.2e1 * t2 * t72 + 0.2e1 * t3 * t73 + t69 * t75 + t68 * t76 + 0.2e1 * t53 * t63 + 0.2e1 * t54 * t64 + t26 * t49 + t25 * t50 + 0.2e1 * t51 * t52 + 0.2e1 * t15 * t20 + 0.2e1 * t16 * t21 + (t235 + t275 + t276) * t221 + ((mrSges(3,2) * t270 + 0.2e1 * (Ifges(3,4) + t227) * t221) * t221 + (mrSges(3,1) * t270 + Ifges(7,5) * t80 + Ifges(7,6) * t79 + Ifges(6,5) * t126 + Ifges(6,6) * t125 - Ifges(5,5) * t168 - Ifges(5,6) * t167 + (-0.2e1 * Ifges(3,4) - t227) * t217 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) + (mrSges(4,1) * t212 + t256) * t269 + t213 * t229 - t212 * t228 - Ifges(7,3) - Ifges(6,3) + pkin(7) ^ 2 * t274 - Ifges(5,3)) * t221) * t217) * qJD(2); -(t258 + t257 + t245) * t221 / 0.2e1 - (t123 / 0.2e1 - pkin(4) * t85) * t177 + ((pkin(7) * mrSges(3,2) - Ifges(3,6) + Ifges(7,5) * t267 + Ifges(7,6) * t268 + Ifges(6,5) * t265 + Ifges(6,6) * t266 + Ifges(5,5) * t263 + Ifges(5,6) * t264 + Ifges(4,5) * t212 / 0.2e1 + Ifges(4,6) * t262) * t217 + (-t212 * (Ifges(4,2) * t213 + t255) / 0.2e1 + (Ifges(4,1) * t212 + t254) * t262 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(4,1) * t213 + mrSges(4,2) * t212 - mrSges(3,1)) * pkin(7)) * t221) * qJD(2) + m(4) * (-t162 * t242 + t163 * t241 + (-t153 * t212 + t154 * t213) * qJ(3)) + (t160 / 0.2e1 + t154 * mrSges(4,3) + qJ(3) * t178 + qJD(3) * t190) * t213 + (t161 / 0.2e1 - t153 * mrSges(4,3) - qJ(3) * t179 - qJD(3) * t191) * t212 + (-t15 * t39 + t16 * t40 + t2 * t97 - t3 * t98) * mrSges(7,3) + (t13 * t146 - t14 * t147 - t53 * t95 + t54 * t96) * mrSges(6,3) + t84 * t263 + t83 * t264 + t30 * t265 + t29 * t266 + t8 * t267 + t7 * t268 + t203 * t90 + t192 * t141 - t176 * t124 / 0.2e1 - t167 * t142 / 0.2e1 - t168 * t143 / 0.2e1 - pkin(2) * t169 + t156 * t113 + t157 * t114 + t118 * t158 + t119 * t159 + t164 * t31 + t150 * t57 + t134 * t151 / 0.2e1 + t133 * t152 / 0.2e1 + t125 * t58 / 0.2e1 + t126 * t59 / 0.2e1 + t42 * t111 + t43 * t112 + t69 * t100 / 0.2e1 + t68 * t101 / 0.2e1 + t109 * t99 + t110 * t9 + t95 * t76 / 0.2e1 + t96 * t75 / 0.2e1 + t87 * t63 + t88 * t64 + t91 * t17 + t79 * t18 / 0.2e1 + t80 * t19 / 0.2e1 + t81 * t52 + t5 * t72 + t6 * t73 + t51 * t60 + t26 * t61 / 0.2e1 + t25 * t62 / 0.2e1 + t40 * t49 / 0.2e1 + t39 * t50 / 0.2e1 + t32 * t20 + t33 * t21 + m(6) * (t109 * t164 + t13 * t88 + t14 * t87 + t150 * t260 + t42 * t54 + t43 * t53) + t180 * (mrSges(5,1) * t226 + mrSges(5,2) * t189) + (t102 * t176 - t103 * t177 - t189 * t56 - t226 * t55) * mrSges(5,3) + m(7) * (t110 * t51 + t15 * t6 + t16 * t5 + t2 * t33 + t3 * t32 + t81 * t91) + m(5) * (t102 * t119 + t103 * t118 + t156 * t56 + t157 * t55 + t180 * t203); t96 * t100 + t95 * t101 + 0.2e1 * t110 * t17 + 0.2e1 * t203 * t141 - t226 * t142 + t189 * t143 + t146 * t58 + t147 * t59 - t176 * t152 + 0.2e1 * t164 * t57 + t97 * t18 + t98 * t19 + t39 * t62 + t40 * t61 + 0.2e1 * t81 * t60 - (-0.2e1 * pkin(4) * t99 + t151) * t177 + (t110 * t81 + t32 * t6 + t33 * t5) * t271 + (t118 * t157 + t119 * t156) * t273 + (t164 * t260 + t42 * t88 + t43 * t87) * t272 + 0.2e1 * (-t32 * t39 + t33 * t40 + t5 * t97 - t6 * t98) * mrSges(7,3) + 0.2e1 * (t146 * t42 - t147 * t43 - t87 * t95 + t88 * t96) * mrSges(6,3) + 0.2e1 * (-t118 * t226 - t119 * t189 + t156 * t176 - t157 * t177) * mrSges(5,3) + (qJ(3) * t274 + 0.2e1 * mrSges(4,3)) * (t212 ^ 2 + t213 ^ 2) * qJD(3); m(4) * t208 + m(5) * t180 + m(6) * t109 + m(7) * t51 + t169 + t31 + t9 + t90; m(7) * t81 - t170 - (-m(6) * pkin(4) - mrSges(5,1)) * t177 + t17 + t57; 0; t222 + m(7) * (t144 * t16 + t145 * t15 + t174 * t3 + t175 * t2) + t174 * t20 + t175 * t21 + t144 * t72 + t145 * t73 - t55 * mrSges(5,2) + t56 * mrSges(5,1) + (-t112 * t239 + t111 * t238 + t215 * t64 + t219 * t63 + m(6) * (t13 * t215 + t14 * t219 + t238 * t54 - t239 * t53)) * pkin(4) - t275; m(7) * (t144 * t33 + t145 * t32 + t174 * t6 + t175 * t5) + t119 * mrSges(5,1) - t118 * mrSges(5,2) + (t144 * t97 - t145 * t98 - t174 * t39 + t175 * t40) * mrSges(7,3) + (m(6) * (t215 * t42 + t219 * t43 + (-t215 * t87 + t219 * t88) * qJD(5)) + (t215 * t96 - t219 * t95 + (t146 * t219 + t147 * t215) * qJD(5)) * mrSges(6,3)) * pkin(4) + t223 + t245; 0; 0.2e1 * t138 + (t144 * t175 + t145 * t174) * t271 - 0.2e1 * t253 + 0.2e1 * t224; (-t73 * t237 + t218 * t20 + m(7) * (-t15 * t237 + t16 * t236 + t2 * t214 + t218 * t3) + t72 * t236 + t214 * t21) * pkin(5) + t222; (m(7) * (t214 * t5 + t218 * t6 + (-t214 * t32 + t218 * t33) * qJD(6)) + (t214 * t40 - t218 * t39 + (t214 * t98 + t218 * t97) * qJD(6)) * mrSges(7,3)) * pkin(5) + t223; 0; t224 + (m(7) * (t144 * t214 + t145 * t218 - t174 * t237 + t175 * t236) - mrSges(7,2) * t236 - mrSges(7,1) * t237) * pkin(5) + t232; 0.2e1 * t182; t225; t230; 0; t232; t182; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
