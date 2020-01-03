% Calculate time derivative of joint inertia matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:41
% EndTime: 2019-12-31 21:26:51
% DurationCPUTime: 3.89s
% Computational Cost: add. (5007->454), mult. (13366->695), div. (0->0), fcn. (12665->10), ass. (0->206)
t252 = Ifges(4,3) + Ifges(5,3);
t170 = cos(qJ(3));
t163 = sin(pkin(10));
t167 = sin(qJ(3));
t211 = t163 * t167;
t216 = cos(pkin(10));
t173 = t170 * t216 - t211;
t129 = t173 * qJD(3);
t183 = t216 * t167;
t135 = t163 * t170 + t183;
t169 = cos(qJ(5));
t166 = sin(qJ(5));
t200 = qJD(5) * t166;
t174 = -t169 * t129 + t135 * t200;
t234 = t166 / 0.2e1;
t232 = t169 / 0.2e1;
t186 = -t200 / 0.2e1;
t165 = cos(pkin(5));
t168 = sin(qJ(2));
t164 = sin(pkin(5));
t171 = cos(qJ(2));
t209 = t164 * t171;
t133 = t165 * t168 * pkin(1) + pkin(7) * t209;
t117 = pkin(8) * t165 + t133;
t118 = (-pkin(2) * t171 - pkin(8) * t168 - pkin(1)) * t164;
t77 = t170 * t117 + t167 * t118;
t210 = t164 * t168;
t154 = pkin(7) * t210;
t229 = pkin(1) * t171;
t132 = t165 * t229 - t154;
t199 = qJD(5) * t169;
t175 = t166 * t129 + t135 * t199;
t128 = t135 * qJD(3);
t201 = qJD(3) * t170;
t202 = qJD(3) * t167;
t251 = Ifges(4,5) * t201 + Ifges(5,5) * t129 - Ifges(4,6) * t202 - Ifges(5,6) * t128;
t250 = 2 * m(5);
t249 = 2 * m(6);
t248 = -2 * mrSges(3,3);
t247 = -2 * mrSges(5,3);
t226 = -qJ(4) - pkin(8);
t184 = qJD(3) * t226;
t127 = qJD(4) * t170 + t167 * t184;
t172 = -qJD(4) * t167 + t170 * t184;
t80 = t127 * t163 - t172 * t216;
t246 = 0.2e1 * t80;
t145 = t226 * t170;
t103 = -t145 * t163 - t183 * t226;
t245 = 0.2e1 * t103;
t122 = t133 * qJD(2);
t244 = 0.2e1 * t122;
t243 = m(5) * pkin(3);
t203 = qJD(2) * t168;
t192 = t164 * t203;
t131 = t165 * t167 + t170 * t210;
t204 = qJD(2) * t164;
t191 = t171 * t204;
t101 = -qJD(3) * t131 - t167 * t191;
t130 = t165 * t170 - t167 * t210;
t102 = qJD(3) * t130 + t170 * t191;
t63 = t163 * t101 + t102 * t216;
t86 = t163 * t130 + t131 * t216;
t68 = -t166 * t86 - t169 * t209;
t34 = qJD(5) * t68 + t166 * t192 + t169 * t63;
t242 = t34 / 0.2e1;
t176 = t166 * t209 - t169 * t86;
t35 = qJD(5) * t176 - t166 * t63 + t169 * t192;
t241 = t35 / 0.2e1;
t240 = t68 / 0.2e1;
t161 = Ifges(6,5) * t199;
t239 = Ifges(6,6) * t186 + t161 / 0.2e1;
t223 = Ifges(6,4) * t166;
t180 = Ifges(6,1) * t169 - t223;
t142 = t180 * qJD(5);
t238 = t142 / 0.2e1;
t237 = Ifges(6,5) * t234 + Ifges(6,6) * t232;
t235 = -t166 / 0.2e1;
t233 = t167 / 0.2e1;
t231 = t170 / 0.2e1;
t228 = pkin(3) * t163;
t76 = -t167 * t117 + t170 * t118;
t52 = -pkin(3) * t209 - t131 * qJ(4) + t76;
t64 = qJ(4) * t130 + t77;
t29 = t163 * t52 + t216 * t64;
t23 = -pkin(9) * t209 + t29;
t85 = -t130 * t216 + t131 * t163;
t116 = t154 + (-pkin(2) - t229) * t165;
t87 = -t130 * pkin(3) + t116;
t36 = t85 * pkin(4) - t86 * pkin(9) + t87;
t11 = t166 * t36 + t169 * t23;
t62 = -t101 * t216 + t102 * t163;
t75 = -t101 * pkin(3) + t122;
t15 = t62 * pkin(4) - t63 * pkin(9) + t75;
t120 = (pkin(2) * t168 - pkin(8) * t171) * t204;
t121 = t132 * qJD(2);
t44 = -t77 * qJD(3) + t170 * t120 - t121 * t167;
t24 = pkin(3) * t192 - qJ(4) * t102 - qJD(4) * t131 + t44;
t43 = -t117 * t202 + t118 * t201 + t167 * t120 + t170 * t121;
t30 = qJ(4) * t101 + qJD(4) * t130 + t43;
t9 = t163 * t24 + t216 * t30;
t7 = pkin(9) * t192 + t9;
t2 = -qJD(5) * t11 + t15 * t169 - t166 * t7;
t227 = t166 * t2;
t225 = Ifges(4,4) * t167;
t224 = Ifges(4,4) * t170;
t222 = Ifges(6,4) * t169;
t221 = Ifges(6,6) * t166;
t220 = t103 * t80;
t219 = t121 * mrSges(3,2);
t104 = -t145 * t216 + t211 * t226;
t160 = -pkin(3) * t170 - pkin(2);
t92 = -pkin(4) * t173 - pkin(9) * t135 + t160;
t56 = t104 * t169 + t166 * t92;
t81 = t127 * t216 + t163 * t172;
t196 = pkin(3) * t202;
t82 = pkin(4) * t128 - pkin(9) * t129 + t196;
t17 = -qJD(5) * t56 - t166 * t81 + t169 * t82;
t218 = t166 * t17;
t217 = t169 * mrSges(6,3);
t215 = t135 * t166;
t214 = t135 * t169;
t158 = pkin(9) + t228;
t213 = t158 * t166;
t212 = t158 * t169;
t149 = Ifges(6,1) * t166 + t222;
t206 = t169 * t149;
t198 = 0.2e1 * t164;
t3 = Ifges(6,5) * t34 + Ifges(6,6) * t35 + Ifges(6,3) * t62;
t195 = mrSges(6,3) * t200;
t194 = mrSges(6,3) * t199;
t193 = t216 * pkin(3);
t190 = t158 * t200;
t189 = t158 * t199;
t33 = t62 * mrSges(5,1) + t63 * mrSges(5,2);
t185 = t199 / 0.2e1;
t88 = t128 * mrSges(5,1) + t129 * mrSges(5,2);
t144 = -mrSges(6,1) * t169 + mrSges(6,2) * t166;
t181 = mrSges(6,1) * t166 + mrSges(6,2) * t169;
t179 = -Ifges(6,2) * t166 + t222;
t10 = -t166 * t23 + t169 * t36;
t178 = -t167 * t44 + t170 * t43;
t55 = -t104 * t166 + t169 * t92;
t177 = Ifges(4,5) * t102 + Ifges(5,5) * t63 + Ifges(4,6) * t101 - Ifges(5,6) * t62 + t192 * t252;
t8 = -t163 * t30 + t216 * t24;
t28 = -t163 * t64 + t216 * t52;
t38 = -Ifges(6,5) * t174 - Ifges(6,6) * t175 + Ifges(6,3) * t128;
t159 = -t193 - pkin(4);
t153 = Ifges(3,5) * t191;
t150 = Ifges(4,1) * t167 + t224;
t148 = Ifges(4,2) * t170 + t225;
t147 = Ifges(6,2) * t169 + t223;
t143 = (Ifges(4,1) * t170 - t225) * qJD(3);
t141 = (-Ifges(4,2) * t167 + t224) * qJD(3);
t140 = t179 * qJD(5);
t138 = (mrSges(4,1) * t167 + mrSges(4,2) * t170) * qJD(3);
t137 = t181 * qJD(5);
t106 = -mrSges(4,1) * t209 - t131 * mrSges(4,3);
t105 = mrSges(4,2) * t209 + t130 * mrSges(4,3);
t98 = Ifges(5,1) * t135 + Ifges(5,4) * t173;
t97 = Ifges(5,4) * t135 + Ifges(5,2) * t173;
t96 = -mrSges(5,1) * t173 + mrSges(5,2) * t135;
t94 = -mrSges(6,1) * t173 - mrSges(6,3) * t214;
t93 = mrSges(6,2) * t173 - mrSges(6,3) * t215;
t91 = t181 * t135;
t90 = Ifges(5,1) * t129 - Ifges(5,4) * t128;
t89 = Ifges(5,4) * t129 - Ifges(5,2) * t128;
t84 = mrSges(4,1) * t192 - mrSges(4,3) * t102;
t83 = -mrSges(4,2) * t192 + mrSges(4,3) * t101;
t79 = Ifges(4,1) * t131 + Ifges(4,4) * t130 - Ifges(4,5) * t209;
t78 = Ifges(4,4) * t131 + Ifges(4,2) * t130 - Ifges(4,6) * t209;
t74 = -mrSges(5,1) * t209 - t86 * mrSges(5,3);
t73 = mrSges(5,2) * t209 - t85 * mrSges(5,3);
t72 = -Ifges(6,5) * t173 + t135 * t180;
t71 = -Ifges(6,6) * t173 + t135 * t179;
t70 = -Ifges(6,3) * t173 + (Ifges(6,5) * t169 - t221) * t135;
t67 = -mrSges(6,2) * t128 - mrSges(6,3) * t175;
t66 = mrSges(6,1) * t128 + mrSges(6,3) * t174;
t65 = -mrSges(4,1) * t101 + mrSges(4,2) * t102;
t54 = Ifges(4,1) * t102 + Ifges(4,4) * t101 + Ifges(4,5) * t192;
t53 = Ifges(4,4) * t102 + Ifges(4,2) * t101 + Ifges(4,6) * t192;
t50 = mrSges(5,1) * t192 - mrSges(5,3) * t63;
t49 = -mrSges(5,2) * t192 - mrSges(5,3) * t62;
t48 = mrSges(6,1) * t175 - mrSges(6,2) * t174;
t47 = mrSges(5,1) * t85 + mrSges(5,2) * t86;
t46 = Ifges(5,1) * t86 - Ifges(5,4) * t85 - Ifges(5,5) * t209;
t45 = Ifges(5,4) * t86 - Ifges(5,2) * t85 - Ifges(5,6) * t209;
t42 = mrSges(6,1) * t85 + mrSges(6,3) * t176;
t41 = -mrSges(6,2) * t85 + mrSges(6,3) * t68;
t40 = -Ifges(6,1) * t174 - Ifges(6,4) * t175 + Ifges(6,5) * t128;
t39 = -Ifges(6,4) * t174 - Ifges(6,2) * t175 + Ifges(6,6) * t128;
t37 = -mrSges(6,1) * t68 - mrSges(6,2) * t176;
t26 = Ifges(5,1) * t63 - Ifges(5,4) * t62 + Ifges(5,5) * t192;
t25 = Ifges(5,4) * t63 - Ifges(5,2) * t62 + Ifges(5,6) * t192;
t22 = pkin(4) * t209 - t28;
t21 = -Ifges(6,1) * t176 + Ifges(6,4) * t68 + Ifges(6,5) * t85;
t20 = -Ifges(6,4) * t176 + Ifges(6,2) * t68 + Ifges(6,6) * t85;
t19 = -Ifges(6,5) * t176 + Ifges(6,6) * t68 + Ifges(6,3) * t85;
t16 = qJD(5) * t55 + t166 * t82 + t169 * t81;
t14 = -mrSges(6,2) * t62 + mrSges(6,3) * t35;
t13 = mrSges(6,1) * t62 - mrSges(6,3) * t34;
t12 = -mrSges(6,1) * t35 + mrSges(6,2) * t34;
t6 = -pkin(4) * t192 - t8;
t5 = Ifges(6,1) * t34 + Ifges(6,4) * t35 + Ifges(6,5) * t62;
t4 = Ifges(6,4) * t34 + Ifges(6,2) * t35 + Ifges(6,6) * t62;
t1 = qJD(5) * t10 + t15 * t166 + t169 * t7;
t18 = [(-t45 + t19) * t62 + (-t25 + t3) * t85 + (mrSges(3,3) * t168 * t244 + (0.2e1 * t121 * mrSges(3,3) - t177) * t171 + ((t132 * t248 + Ifges(3,5) * t165 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t171) * t198) * t171 + (t133 * t248 + Ifges(4,5) * t131 + Ifges(5,5) * t86 - 0.2e1 * Ifges(3,6) * t165 + Ifges(4,6) * t130 - Ifges(5,6) * t85 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t168) * t198 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - t252) * t209) * t168) * qJD(2)) * t164 + 0.2e1 * m(4) * (t116 * t122 + t43 * t77 + t44 * t76) + 0.2e1 * m(3) * (t121 * t133 - t122 * t132) + (-0.2e1 * t122 * mrSges(3,1) + t153 - 0.2e1 * t219) * t165 + (-mrSges(4,1) * t130 + mrSges(4,2) * t131) * t244 + (t1 * t11 + t10 * t2 + t22 * t6) * t249 + (t28 * t8 + t29 * t9 + t75 * t87) * t250 + t131 * t54 + t130 * t53 + 0.2e1 * t116 * t65 + 0.2e1 * t43 * t105 + 0.2e1 * t44 * t106 + t101 * t78 + t102 * t79 + t86 * t26 + 0.2e1 * t87 * t33 + 0.2e1 * t77 * t83 + 0.2e1 * t76 * t84 + 0.2e1 * t9 * t73 + 0.2e1 * t8 * t74 + 0.2e1 * t75 * t47 + t68 * t4 + t63 * t46 + 0.2e1 * t29 * t49 + 0.2e1 * t28 * t50 + 0.2e1 * t6 * t37 + 0.2e1 * t1 * t41 + 0.2e1 * t2 * t42 + t34 * t21 + t35 * t20 + 0.2e1 * t22 * t12 + 0.2e1 * t10 * t13 + 0.2e1 * t11 * t14 - t176 * t5; t153 - t219 + (t79 * t231 + (pkin(3) * t47 - t78 / 0.2e1) * t167) * qJD(3) + (t46 / 0.2e1 - t28 * mrSges(5,3) + t20 * t235 + t21 * t232) * t129 + (t26 / 0.2e1 - t8 * mrSges(5,3) + t4 * t235 + t5 * t232 + (-t169 * t20 / 0.2e1 + t21 * t235) * qJD(5)) * t135 + (-t106 * t201 - t105 * t202 - t167 * t84 + t170 * t83 + m(4) * (-t201 * t76 - t202 * t77 + t178)) * pkin(8) + m(5) * (-t103 * t8 + t104 * t9 + t160 * t75 + t196 * t87 - t28 * t80 + t29 * t81) + t39 * t240 + t71 * t241 + t72 * t242 + t53 * t231 + t54 * t233 + ((-t167 * t77 - t170 * t76) * qJD(3) + t178) * mrSges(4,3) + (t37 - t74) * t80 + (-mrSges(4,1) * t170 + mrSges(4,2) * t167 - mrSges(3,1)) * t122 + (t12 - t50) * t103 + (-m(4) * t122 - t65) * pkin(2) + m(6) * (t1 * t56 + t10 * t17 + t103 * t6 + t11 * t16 + t2 * t55 + t22 * t80) + (-t89 / 0.2e1 + t38 / 0.2e1) * t85 + t160 * t33 + t101 * t148 / 0.2e1 + t102 * t150 / 0.2e1 + t130 * t141 / 0.2e1 + t131 * t143 / 0.2e1 + t116 * t138 + t104 * t49 + t86 * t90 / 0.2e1 + t6 * t91 + t1 * t93 + t2 * t94 + t75 * t96 + t63 * t98 / 0.2e1 + t87 * t88 + t81 * t73 + t11 * t67 + t10 * t66 + t55 * t13 + t56 * t14 + t22 * t48 + t16 * t41 + t17 * t42 + (-t97 / 0.2e1 + t70 / 0.2e1) * t62 - t176 * t40 / 0.2e1 + ((Ifges(4,5) * t233 + Ifges(4,6) * t231 + Ifges(5,5) * t135 / 0.2e1 + Ifges(5,6) * t173 / 0.2e1 - Ifges(3,6)) * t203 - t251 * t171 / 0.2e1) * t164 - (-t25 / 0.2e1 + t3 / 0.2e1 - t9 * mrSges(5,3)) * t173 + (-t45 / 0.2e1 + t19 / 0.2e1 - t29 * mrSges(5,3)) * t128; -0.2e1 * pkin(2) * t138 + t48 * t245 + t170 * t141 + t167 * t143 + 0.2e1 * t16 * t93 + 0.2e1 * t160 * t88 + 0.2e1 * t17 * t94 + 0.2e1 * t55 * t66 + 0.2e1 * t56 * t67 + t91 * t246 + (t170 * t150 + (0.2e1 * pkin(3) * t96 - t148) * t167) * qJD(3) + (t16 * t56 + t17 * t55 + t220) * t249 + (t104 * t81 + t160 * t196 + t220) * t250 - (t247 * t81 + t38 - t89) * t173 + (t104 * t247 + t70 - t97) * t128 + (mrSges(5,3) * t245 - t166 * t71 + t169 * t72 + t98) * t129 + (mrSges(5,3) * t246 - t166 * t39 + t169 * t40 + t90 + (-t166 * t72 - t169 * t71) * qJD(5)) * t135; t177 - mrSges(6,3) * t227 + m(6) * (t159 * t6 + (t1 * t169 - t227 + (-t10 * t169 - t11 * t166) * qJD(5)) * t158) - t13 * t213 + t21 * t185 + t20 * t186 - t10 * t194 - t11 * t195 - t42 * t189 - t41 * t190 + t5 * t234 + t62 * t237 + t85 * t239 + t140 * t240 + t147 * t241 + t149 * t242 + (t163 * t9 + t216 * t8) * t243 + t49 * t228 + t4 * t232 + t14 * t212 + t1 * t217 + t50 * t193 + t159 * t12 + t6 * t144 + t22 * t137 - t43 * mrSges(4,2) + t44 * mrSges(4,1) + t8 * mrSges(5,1) - t9 * mrSges(5,2) - t176 * t238; -t175 * t147 / 0.2e1 + (-mrSges(4,1) * t201 + mrSges(4,2) * t202) * pkin(8) + (-t128 * t228 - t129 * t193) * mrSges(5,3) + (t135 * t149 + t71) * t186 + (t163 * t243 - mrSges(5,2)) * t81 + (m(6) * t159 - t216 * t243 - mrSges(5,1) + t144) * t80 - mrSges(6,3) * t218 - t66 * t213 - t140 * t215 / 0.2e1 + t129 * t206 / 0.2e1 + t72 * t185 - t55 * t194 - t56 * t195 - t94 * t189 - t93 * t190 + t40 * t234 + t128 * t237 + t214 * t238 + t39 * t232 + t67 * t212 + t16 * t217 + m(6) * (t16 * t169 - t218 + (-t166 * t56 - t169 * t55) * qJD(5)) * t158 + t159 * t48 + t103 * t137 - t173 * t239 + t251; 0.2e1 * t137 * t159 + t140 * t169 + t142 * t166 + (-t147 * t166 + t206) * qJD(5); t169 * t13 + t166 * t14 + (-t166 * t42 + t169 * t41) * qJD(5) + m(6) * (t1 * t166 + t169 * t2 + (-t10 * t166 + t11 * t169) * qJD(5)) + m(5) * t75 + t33; t93 * t199 + t166 * t67 + m(6) * (t16 * t166 + t169 * t17 + (-t166 * t55 + t169 * t56) * qJD(5)) - t94 * t200 + t169 * t66 + m(5) * t196 + t88; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t17 - mrSges(6,2) * t16 + t38; t161 + (t144 * t158 - t221) * qJD(5); -t137; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t18(1), t18(2), t18(4), t18(7), t18(11); t18(2), t18(3), t18(5), t18(8), t18(12); t18(4), t18(5), t18(6), t18(9), t18(13); t18(7), t18(8), t18(9), t18(10), t18(14); t18(11), t18(12), t18(13), t18(14), t18(15);];
Mq = res;
