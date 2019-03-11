% Calculate time derivative of joint inertia matrix for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:17
% EndTime: 2019-03-09 09:17:25
% DurationCPUTime: 3.69s
% Computational Cost: add. (3792->456), mult. (9279->664), div. (0->0), fcn. (7963->8), ass. (0->201)
t145 = sin(qJ(6));
t148 = cos(qJ(6));
t200 = t145 ^ 2 + t148 ^ 2;
t253 = m(7) * (-0.1e1 + t200);
t252 = Ifges(4,4) + Ifges(3,5);
t142 = sin(pkin(6));
t147 = sin(qJ(2));
t206 = t142 * t147;
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t150 = cos(qJ(2));
t205 = t142 * t150;
t77 = -t142 * pkin(1) - pkin(2) * t205 - qJ(3) * t206;
t67 = pkin(3) * t205 - t77;
t42 = (pkin(4) * t147 + pkin(9) * t150) * t142 + t67;
t143 = cos(pkin(6));
t125 = pkin(8) * t206;
t172 = -qJ(4) * t206 + t125;
t229 = pkin(1) * t150;
t183 = -pkin(2) - t229;
t173 = -pkin(3) + t183;
t51 = (-pkin(9) + t173) * t143 + t172;
t220 = t146 * t42 + t149 * t51;
t17 = pkin(10) * t206 + t220;
t129 = t143 * t147 * pkin(1);
t87 = pkin(8) * t205 + t129;
t76 = t143 * qJ(3) + t87;
t66 = qJ(4) * t205 - t76;
t56 = t143 * pkin(4) - t66;
t84 = -t143 * t149 + t146 * t205;
t186 = t149 * t205;
t85 = -t143 * t146 - t186;
t29 = -pkin(5) * t84 - pkin(10) * t85 + t56;
t10 = -t145 * t17 + t148 * t29;
t11 = t145 * t29 + t148 * t17;
t198 = qJD(2) * t142;
t180 = t150 * t198;
t197 = qJD(2) * t147;
t181 = t142 * t197;
t62 = qJD(5) * t84 + t149 * t181;
t64 = t145 * t206 + t148 * t85;
t25 = -qJD(6) * t64 - t145 * t62 + t148 * t180;
t63 = -t145 * t85 + t148 * t206;
t26 = qJD(6) * t63 + t145 * t180 + t148 * t62;
t12 = -mrSges(7,1) * t25 + mrSges(7,2) * t26;
t248 = qJD(5) * t220;
t44 = mrSges(6,1) * t180 - mrSges(6,3) * t62;
t151 = -pkin(2) - pkin(3);
t137 = -pkin(9) + t151;
t202 = qJ(3) * t180 + qJD(3) * t206;
t35 = (pkin(4) * t150 + t137 * t147) * t198 + t202;
t54 = -qJD(4) * t206 + (t129 + (pkin(8) - qJ(4)) * t205) * qJD(2);
t9 = -t146 * t54 + t149 * t35 - t248;
t160 = -m(6) * (t9 + t248) - t44 + t12;
t195 = qJD(5) * t148;
t36 = mrSges(7,2) * t84 + mrSges(7,3) * t63;
t37 = -mrSges(7,1) * t84 - mrSges(7,3) * t64;
t4 = -pkin(5) * t180 - t9;
t69 = -mrSges(6,2) * t206 + mrSges(6,3) * t84;
t251 = t160 + m(7) * (qJD(5) * t10 * t145 - t11 * t195 + t4) + qJD(5) * (t145 * t37 - t148 * t36 - t69);
t191 = qJD(6) * t148;
t194 = qJD(5) * t149;
t158 = t145 * t194 + t146 * t191;
t249 = t146 ^ 2 - t149 ^ 2;
t188 = t143 * t229;
t118 = qJD(2) * t188;
t133 = t143 * qJD(3);
t48 = t142 * (pkin(8) * t197 + qJD(4) * t150) - qJ(4) * t181 - t118 - t133;
t246 = 0.2e1 * m(7);
t245 = -0.2e1 * pkin(1);
t244 = 2 * mrSges(4,2);
t243 = -2 * mrSges(3,3);
t242 = -2 * mrSges(5,3);
t241 = 0.2e1 * t137;
t144 = qJ(3) + pkin(4);
t240 = 0.2e1 * t144;
t239 = m(7) * t4;
t238 = t63 / 0.2e1;
t237 = t64 / 0.2e1;
t216 = Ifges(7,4) * t145;
t168 = Ifges(7,1) * t148 - t216;
t81 = Ifges(7,5) * t149 - t146 * t168;
t236 = t81 / 0.2e1;
t235 = t146 * t194 * t253;
t106 = Ifges(7,2) * t148 + t216;
t234 = t106 / 0.2e1;
t233 = t145 / 0.2e1;
t232 = -t146 / 0.2e1;
t231 = t146 / 0.2e1;
t230 = t148 / 0.2e1;
t228 = pkin(5) * t146;
t227 = pkin(5) * t149;
t226 = pkin(10) * t146;
t225 = pkin(10) * t149;
t224 = t54 * mrSges(5,2);
t82 = -pkin(8) * t181 + t118;
t223 = t82 * mrSges(3,2);
t222 = mrSges(4,2) - mrSges(5,3);
t219 = -mrSges(5,1) * t143 + mrSges(6,1) * t84 - mrSges(6,2) * t85 + mrSges(5,3) * t205;
t218 = Ifges(6,4) * t146;
t217 = Ifges(6,4) * t149;
t215 = Ifges(7,4) * t148;
t214 = Ifges(7,6) * t145;
t213 = t146 * mrSges(6,2);
t212 = t146 * mrSges(7,3);
t169 = mrSges(7,1) * t145 + mrSges(7,2) * t148;
t92 = t169 * qJD(6);
t211 = t146 * t92;
t83 = t87 * qJD(2);
t210 = t147 * t83;
t18 = -t146 * t51 + t149 * t42;
t16 = -pkin(5) * t206 - t18;
t209 = qJD(5) * t16;
t101 = t144 + t226 + t227;
t204 = t145 * t149;
t71 = t101 * t148 - t137 * t204;
t208 = qJD(6) * t71;
t203 = t148 * t149;
t72 = t101 * t145 + t137 * t203;
t207 = qJD(6) * t72;
t201 = mrSges(5,1) * t180 + mrSges(5,2) * t181;
t196 = qJD(5) * t146;
t193 = qJD(6) * t145;
t192 = qJD(6) * t146;
t190 = qJD(6) * t149;
t189 = Ifges(5,4) - Ifges(4,5) + Ifges(3,4);
t61 = -qJD(5) * t186 - t143 * t196 + t146 * t181;
t5 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t61;
t187 = Ifges(6,5) * t62 - Ifges(6,6) * t61 + Ifges(6,3) * t180;
t184 = Ifges(4,6) * t181 + t180 * t252;
t182 = (-mrSges(4,1) - mrSges(3,1)) * t83;
t179 = t145 * t196;
t178 = t146 * t195;
t177 = t145 * t192;
t15 = pkin(5) * t61 - pkin(10) * t62 - t48;
t8 = t146 * t35 + t149 * t54 + t42 * t194 - t196 * t51;
t3 = pkin(10) * t180 + t8;
t1 = qJD(6) * t10 + t145 * t15 + t148 * t3;
t2 = -qJD(6) * t11 - t145 * t3 + t148 * t15;
t171 = t1 * t148 - t2 * t145;
t104 = t149 * mrSges(6,1) - t213;
t170 = -t146 * mrSges(6,1) - t149 * mrSges(6,2);
t103 = -mrSges(7,1) * t148 + mrSges(7,2) * t145;
t108 = Ifges(7,1) * t145 + t215;
t167 = -Ifges(7,2) * t145 + t215;
t13 = -mrSges(7,2) * t61 + mrSges(7,3) * t25;
t14 = mrSges(7,1) * t61 - mrSges(7,3) * t26;
t166 = t148 * t13 - t145 * t14;
t21 = Ifges(7,4) * t64 + Ifges(7,2) * t63 - Ifges(7,6) * t84;
t22 = Ifges(7,1) * t64 + Ifges(7,4) * t63 - Ifges(7,5) * t84;
t162 = -t145 * t21 / 0.2e1 + t22 * t230;
t159 = t148 * t194 - t177;
t31 = -mrSges(7,1) * t63 + mrSges(7,2) * t64;
t43 = -mrSges(6,2) * t180 - mrSges(6,3) * t61;
t70 = mrSges(6,1) * t206 - mrSges(6,3) * t85;
t156 = t43 + m(6) * (-qJD(5) * t18 + t8) + (t31 - t70) * qJD(5);
t45 = Ifges(7,5) * t177 + (-Ifges(7,5) * t203 - Ifges(7,3) * t146) * qJD(5) + t158 * Ifges(7,6);
t155 = -t10 * t191 - t11 * t193 + t171;
t89 = qJD(3) + (t225 - t228) * qJD(5);
t33 = -t137 * t178 + t145 * t89 + t208;
t34 = t137 * t179 + t148 * t89 - t207;
t154 = -t145 * t34 + t148 * t33 - t191 * t71 - t193 * t72;
t100 = mrSges(7,1) * t149 + t148 * t212;
t73 = -mrSges(7,1) * t196 + mrSges(7,3) * t159;
t74 = mrSges(7,2) * t196 + mrSges(7,3) * t158;
t88 = t169 * t146;
t99 = -mrSges(7,2) * t149 + t145 * t212;
t153 = -qJD(5) * t88 - t100 * t191 - t145 * t73 + t148 * t74 - t193 * t99;
t152 = (-t145 * t36 - t148 * t37) * qJD(6) + m(7) * (t155 + t209) + t156 + t166;
t132 = Ifges(7,5) * t191;
t131 = Ifges(6,6) * t196;
t109 = -Ifges(6,1) * t146 - t217;
t107 = -Ifges(6,2) * t149 - t218;
t105 = Ifges(7,5) * t145 + Ifges(7,6) * t148;
t98 = (-Ifges(6,1) * t149 + t218) * qJD(5);
t97 = t168 * qJD(6);
t96 = (Ifges(6,2) * t146 - t217) * qJD(5);
t95 = t167 * qJD(6);
t94 = -Ifges(7,6) * t193 + t132;
t93 = t170 * qJD(5);
t91 = mrSges(4,2) * t205 + mrSges(4,3) * t143;
t86 = -t125 + t188;
t80 = Ifges(7,6) * t149 - t146 * t167;
t79 = Ifges(7,3) * t149 + (-Ifges(7,5) * t148 + t214) * t146;
t78 = t143 * t183 + t125;
t75 = t133 + t82;
t68 = pkin(2) * t181 - t202;
t57 = t143 * t173 + t172;
t55 = -mrSges(7,1) * t158 - mrSges(7,2) * t159;
t53 = t151 * t181 + t202;
t47 = t108 * t192 + (-Ifges(7,5) * t146 - t149 * t168) * qJD(5);
t46 = t106 * t192 + (-Ifges(7,6) * t146 - t149 * t167) * qJD(5);
t41 = Ifges(6,1) * t85 + Ifges(6,4) * t84 + Ifges(6,5) * t206;
t40 = Ifges(6,4) * t85 + Ifges(6,2) * t84 + Ifges(6,6) * t206;
t30 = mrSges(6,1) * t61 + mrSges(6,2) * t62;
t28 = Ifges(6,1) * t62 - Ifges(6,4) * t61 + Ifges(6,5) * t180;
t27 = Ifges(6,4) * t62 - Ifges(6,2) * t61 + Ifges(6,6) * t180;
t20 = Ifges(7,5) * t64 + Ifges(7,6) * t63 - Ifges(7,3) * t84;
t7 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t61;
t6 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t61;
t19 = [(t1 * t11 + t10 * t2 + t16 * t4) * t246 + (t20 - t40) * t61 + (t27 - t5) * t84 + 0.2e1 * t219 * t48 + 0.2e1 * m(5) * (t48 * t66 + t53 * t67 + t54 * t57) + 0.2e1 * m(4) * (t68 * t77 + t75 * t76 + t78 * t83) + 0.2e1 * m(3) * (t82 * t87 - t83 * t86) + 0.2e1 * t67 * t201 + (0.2e1 * t182 + t184 - 0.2e1 * t223 + 0.2e1 * t224) * t143 + 0.2e1 * m(6) * (t18 * t9 + t220 * t8 - t48 * t56) + 0.2e1 * t220 * t43 + (t147 * t187 + t54 * t147 * t242 + t210 * t244 + 0.2e1 * t53 * (mrSges(5,1) * t147 - mrSges(5,2) * t150) + 0.2e1 * t68 * (-mrSges(4,1) * t150 - mrSges(4,3) * t147) + 0.2e1 * (t150 * t82 + t210) * mrSges(3,3) + ((0.2e1 * t77 * mrSges(4,1) - 0.2e1 * t76 * mrSges(4,2) + t87 * t243 + t66 * t242 + (mrSges(3,1) * t245 - 0.2e1 * t147 * t189) * t142 + (-(2 * Ifges(5,5)) - (2 * Ifges(3,6)) + Ifges(4,6)) * t143) * t147 + (t78 * t244 + t86 * t243 - 0.2e1 * t77 * mrSges(4,3) + t57 * t242 + Ifges(6,5) * t85 + Ifges(6,6) * t84 + (mrSges(3,2) * t245 + 0.2e1 * t150 * t189) * t142 + ((2 * Ifges(5,6)) + t252) * t143 + ((2 * Ifges(3,1)) + (2 * Ifges(4,1)) - (2 * Ifges(5,1)) - (2 * Ifges(3,2)) + (2 * Ifges(5,2)) - (2 * Ifges(4,3)) + Ifges(6,3)) * t206) * t150) * qJD(2)) * t142 + 0.2e1 * t11 * t13 + 0.2e1 * t10 * t14 + 0.2e1 * t16 * t12 + t25 * t21 + t26 * t22 + 0.2e1 * t4 * t31 + 0.2e1 * t1 * t36 + 0.2e1 * t2 * t37 + 0.2e1 * t18 * t44 + 0.2e1 * t56 * t30 + t62 * t41 + t63 * t6 + t64 * t7 + 0.2e1 * t8 * t69 + 0.2e1 * t9 * t70 + t85 * t28 + 0.2e1 * t75 * t91; t184 + t26 * t236 + t47 * t237 + t46 * t238 + t182 + m(4) * (-pkin(2) * t83 + qJ(3) * t75 + qJD(3) * t76) + m(5) * (-qJ(3) * t48 - qJD(3) * t66 + t151 * t54) + (-t45 / 0.2e1 + t96 / 0.2e1) * t84 + (t147 * (-Ifges(6,5) * t194 + t131) / 0.2e1 + ((Ifges(6,5) * t232 - Ifges(6,6) * t149 / 0.2e1 + Ifges(5,6) - pkin(2) * mrSges(4,2) - t151 * mrSges(5,3)) * t150 + (-qJ(3) * t222 - Ifges(5,5) - Ifges(3,6)) * t147) * qJD(2)) * t142 + (t91 - t219) * qJD(3) + (-t8 * mrSges(6,3) + t5 / 0.2e1 - t27 / 0.2e1 + (t18 * mrSges(6,3) - t41 / 0.2e1 - t162) * qJD(5) + (m(7) * t209 + t156) * t137) * t149 + (t79 / 0.2e1 - t107 / 0.2e1) * t61 + (-mrSges(5,1) - t104) * t48 + m(7) * (t1 * t72 + t10 * t34 + t11 * t33 + t2 * t71) + m(6) * (qJD(3) * t56 - t144 * t48) + (t9 * mrSges(6,3) + t6 * t233 - t148 * t7 / 0.2e1 - t28 / 0.2e1 + (t21 * t230 + t22 * t233) * qJD(6) + (t220 * mrSges(6,3) - t20 / 0.2e1 + t40 / 0.2e1) * qJD(5) + (-qJD(5) * t69 + t160 + t239) * t137) * t146 - t223 + t224 + t33 * t36 + t34 * t37 + t16 * t55 + t71 * t14 + t72 * t13 + t10 * t73 + t11 * t74 + t75 * mrSges(4,3) + t25 * t80 / 0.2e1 - t4 * t88 + t56 * t93 + t85 * t98 / 0.2e1 + t1 * t99 + t2 * t100 + t62 * t109 / 0.2e1 + t144 * t30; (t33 * t72 + t34 * t71) * t246 + 0.2e1 * t71 * t73 + 0.2e1 * t72 * t74 + 0.2e1 * t33 * t99 + 0.2e1 * t34 * t100 + t93 * t240 + (t45 - t96 + (t145 * t80 - t148 * t81 - t241 * t88 - t109) * qJD(5)) * t149 + (t55 * t241 + t145 * t46 - t148 * t47 - t98 + (t145 * t81 + t148 * t80) * qJD(6) + (t137 ^ 2 * t149 * t246 + t107 - t79) * qJD(5)) * t146 + (m(6) * t240 + 0.2e1 * mrSges(5,1) + 0.2e1 * mrSges(4,3) + 0.2e1 * t104 + 0.2e1 * (m(4) + m(5)) * qJ(3)) * qJD(3); m(4) * t83 + m(5) * t54 + t146 * t251 + t152 * t149 + t222 * t180; (t55 + (m(7) * (t145 * t71 - t148 * t72) - t148 * t99 + t145 * t100) * qJD(5)) * t146 + (m(7) * (t196 * t241 + t154) + t153) * t149; -0.2e1 * t235; m(5) * t53 + t152 * t146 - t149 * t251 + t201; -t149 * t55 + (m(7) * (t137 * t249 + t203 * t72 - t204 * t71) + t99 * t203 - t100 * t204) * qJD(5) + (m(7) * t154 + t153) * t146; -t249 * qJD(5) * t253; 0.2e1 * t235; -t8 * mrSges(6,2) + t9 * mrSges(6,1) + t16 * t92 - t84 * t94 / 0.2e1 + t95 * t238 + t97 * t237 + t4 * t103 + t61 * t105 / 0.2e1 + t25 * t234 + t26 * t108 / 0.2e1 + t7 * t233 + t6 * t230 + t162 * qJD(6) + (-t12 - t239) * pkin(5) + ((-t10 * t148 - t11 * t145) * qJD(6) + t171) * mrSges(7,3) + (m(7) * t155 - t191 * t37 - t193 * t36 + t166) * pkin(10) + t187; t131 + t137 * t211 - pkin(5) * t55 + t149 * t94 / 0.2e1 + (t105 * t232 - Ifges(6,5) * t149 + (t213 + (-m(7) * pkin(5) - mrSges(6,1) + t103) * t149) * t137) * qJD(5) + (-t108 * t194 / 0.2e1 + t33 * mrSges(7,3) + t97 * t232 + t46 / 0.2e1 + (-t71 * mrSges(7,3) + t106 * t231 + t236) * qJD(6) + (m(7) * (t33 - t208) + t74 - qJD(6) * t100) * pkin(10)) * t148 + (t194 * t234 - t34 * mrSges(7,3) + t95 * t231 + t47 / 0.2e1 + (-t80 / 0.2e1 - t72 * mrSges(7,3) + t108 * t231) * qJD(6) + (m(7) * (-t34 - t207) - qJD(6) * t99 - t73) * pkin(10)) * t145; t211 + (t149 * t103 + m(7) * (-t200 * t226 - t227) - t200 * t212 - t104) * qJD(5); -t149 * t92 + (t146 * t103 + m(7) * (t200 * t225 - t228) + t200 * t149 * mrSges(7,3) + t170) * qJD(5); -0.2e1 * pkin(5) * t92 + t145 * t97 + t148 * t95 + (-t106 * t145 + t108 * t148) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t34 - mrSges(7,2) * t33 + t45; (t145 * t190 + t178) * mrSges(7,2) + (-t148 * t190 + t179) * mrSges(7,1); t55; t132 + (pkin(10) * t103 - t214) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
