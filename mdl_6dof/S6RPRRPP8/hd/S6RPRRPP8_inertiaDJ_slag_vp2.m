% Calculate time derivative of joint inertia matrix for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:55
% EndTime: 2019-03-09 04:54:01
% DurationCPUTime: 3.17s
% Computational Cost: add. (1524->428), mult. (3473->572), div. (0->0), fcn. (2223->4), ass. (0->166)
t200 = 2 * qJD(2);
t208 = Ifges(7,4) + Ifges(6,5);
t113 = sin(qJ(4));
t115 = cos(qJ(4));
t204 = t113 ^ 2 + t115 ^ 2;
t116 = cos(qJ(3));
t165 = qJD(4) * t116;
t149 = t113 * t165;
t114 = sin(qJ(3));
t171 = qJD(3) * t114;
t120 = t115 * t171 + t149;
t147 = t113 * t171;
t148 = t115 * t165;
t118 = t147 - t148;
t207 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t206 = mrSges(7,2) + mrSges(6,3);
t178 = qJ(5) * t113;
t125 = -pkin(4) * t115 - t178;
t68 = -pkin(3) + t125;
t72 = mrSges(6,2) * t115 - mrSges(6,3) * t113;
t205 = m(6) * t68 + t72;
t73 = -mrSges(5,1) * t115 + mrSges(5,2) * t113;
t203 = -m(5) * pkin(3) - mrSges(4,1) + t73;
t202 = 2 * m(7);
t201 = 2 * mrSges(7,1);
t198 = pkin(5) + pkin(8);
t196 = mrSges(6,1) + mrSges(7,1);
t194 = -Ifges(5,5) - Ifges(7,5);
t112 = -pkin(4) - qJ(6);
t170 = qJD(3) * t116;
t28 = -mrSges(6,1) * t118 - mrSges(6,3) * t170;
t29 = mrSges(7,1) * t118 + mrSges(7,2) * t170;
t193 = -t28 + t29;
t25 = mrSges(5,1) * t170 + mrSges(5,3) * t120;
t101 = mrSges(6,2) * t170;
t30 = -mrSges(6,1) * t120 + t101;
t192 = t30 - t25;
t173 = t115 * t116;
t60 = mrSges(5,1) * t114 - mrSges(5,3) * t173;
t64 = mrSges(6,1) * t173 + mrSges(6,2) * t114;
t191 = -t60 + t64;
t176 = t113 * t116;
t62 = mrSges(6,1) * t176 - mrSges(6,3) * t114;
t63 = -mrSges(7,1) * t176 + mrSges(7,2) * t114;
t190 = -t62 + t63;
t151 = t115 * t170;
t163 = qJD(5) * t115;
t188 = qJ(5) * t151 + t114 * t163;
t117 = -pkin(1) - pkin(7);
t174 = t114 * t117;
t67 = pkin(3) * t114 - pkin(8) * t116 + qJ(2);
t32 = t113 * t67 + t115 * t174;
t186 = Ifges(5,4) * t113;
t185 = Ifges(5,4) * t115;
t184 = Ifges(6,6) * t113;
t183 = Ifges(6,6) * t115;
t182 = Ifges(7,6) * t113;
t181 = Ifges(7,6) * t115;
t180 = t114 * Ifges(6,4);
t179 = t114 * Ifges(5,6);
t177 = qJ(5) * t115;
t175 = t114 * t115;
t172 = qJ(5) * qJD(5);
t169 = qJD(3) * t117;
t168 = qJD(4) * t113;
t167 = qJD(4) * t114;
t166 = qJD(4) * t115;
t164 = qJD(4) * t117;
t162 = qJD(6) * t114;
t161 = -mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t159 = mrSges(5,2) - t206;
t134 = -Ifges(5,2) * t113 + t185;
t34 = t116 * t134 + t179;
t129 = -Ifges(6,3) * t113 + t183;
t37 = Ifges(6,5) * t114 - t116 * t129;
t130 = Ifges(7,2) * t113 + t181;
t38 = Ifges(7,4) * t114 + t116 * t130;
t158 = t34 - t37 - t38;
t150 = t116 * t169;
t49 = qJD(2) + (pkin(3) * t116 + pkin(8) * t114) * qJD(3);
t157 = t113 * t49 + t115 * t150 + t166 * t67;
t59 = -mrSges(5,2) * t114 - mrSges(5,3) * t176;
t156 = t59 + t190;
t61 = mrSges(7,1) * t173 - mrSges(7,3) * t114;
t155 = t61 + t191;
t154 = pkin(4) * t170;
t82 = t198 * t115;
t152 = qJ(5) * t167;
t146 = t113 * t167;
t145 = t114 * t166;
t143 = t112 * t114;
t92 = t113 * t174;
t31 = t115 * t67 - t92;
t142 = pkin(4) * t148 + qJ(5) * t120 + t114 * t169;
t141 = pkin(4) * t168 - qJD(5) * t113;
t127 = -Ifges(7,3) * t113 + t181;
t132 = Ifges(6,2) * t113 + t183;
t80 = Ifges(5,1) * t113 + t185;
t140 = t127 / 0.2e1 - t132 / 0.2e1 - t80 / 0.2e1;
t128 = Ifges(6,3) * t115 + t184;
t131 = Ifges(7,2) * t115 - t182;
t79 = Ifges(5,2) * t115 + t186;
t139 = -t128 / 0.2e1 - t131 / 0.2e1 - t79 / 0.2e1;
t23 = -qJ(5) * t114 - t32;
t138 = mrSges(5,1) * t113 + mrSges(5,2) * t115;
t137 = -mrSges(6,2) * t113 - mrSges(6,3) * t115;
t136 = -mrSges(7,2) * t115 + mrSges(7,3) * t113;
t135 = Ifges(5,1) * t115 - t186;
t133 = Ifges(6,2) * t115 - t184;
t126 = Ifges(7,3) * t115 + t182;
t7 = -t113 * t150 + t115 * t49 - t117 * t145 - t168 * t67;
t124 = qJ(6) * t113 - t177;
t123 = t115 * (m(6) * pkin(8) + t196);
t35 = Ifges(5,5) * t114 + t116 * t135;
t36 = Ifges(7,5) * t114 + t116 * t126;
t39 = -t116 * t133 + t180;
t122 = t114 * t194 - t35 - t36 + t39;
t121 = Ifges(6,4) * t120 + Ifges(5,6) * t147 + t148 * t208 + t170 * t207;
t106 = Ifges(7,4) * t168;
t105 = Ifges(5,5) * t166;
t104 = Ifges(6,5) * t168;
t103 = Ifges(7,5) * t166;
t99 = pkin(4) * t176;
t81 = t198 * t113;
t74 = -mrSges(7,2) * t113 - mrSges(7,3) * t115;
t66 = qJD(4) * t82;
t65 = t198 * t168;
t58 = t135 * qJD(4);
t57 = t134 * qJD(4);
t56 = t133 * qJD(4);
t55 = t130 * qJD(4);
t54 = t129 * qJD(4);
t53 = t126 * qJD(4);
t52 = t138 * qJD(4);
t51 = t137 * qJD(4);
t50 = t136 * qJD(4);
t47 = t138 * t116;
t46 = t137 * t116;
t45 = t136 * t116;
t44 = t112 * t115 - pkin(3) - t178;
t40 = -qJ(5) * t166 + t141;
t33 = t99 + (-t117 - t177) * t116;
t27 = -mrSges(7,1) * t120 - mrSges(7,3) * t170;
t26 = -mrSges(5,2) * t170 + mrSges(5,3) * t118;
t24 = -pkin(4) * t114 - t31;
t22 = t99 + (-t117 + t124) * t116;
t21 = qJD(4) * t124 - qJD(6) * t115 + t141;
t19 = -pkin(5) * t176 - t23;
t18 = mrSges(7,2) * t120 - mrSges(7,3) * t118;
t17 = -mrSges(5,1) * t118 - mrSges(5,2) * t120;
t16 = mrSges(6,2) * t118 + mrSges(6,3) * t120;
t15 = t92 + (pkin(5) * t116 - t67) * t115 + t143;
t14 = t132 * t165 + (Ifges(6,4) * t116 + t114 * t133) * qJD(3);
t13 = t131 * t165 + (Ifges(7,4) * t116 - t114 * t130) * qJD(3);
t12 = t128 * t165 + (Ifges(6,5) * t116 + t114 * t129) * qJD(3);
t11 = t127 * t165 + (Ifges(7,5) * t116 - t114 * t126) * qJD(3);
t10 = -t80 * t165 + (Ifges(5,5) * t116 - t114 * t135) * qJD(3);
t9 = -t79 * t165 + (Ifges(5,6) * t116 - t114 * t134) * qJD(3);
t8 = -pkin(4) * t147 - t116 * t163 + t142;
t6 = -t117 * t146 + t157;
t5 = -t7 - t154;
t4 = -qJ(5) * t170 + (t113 * t164 - qJD(5)) * t114 - t157;
t3 = (qJ(6) * qJD(4) - qJD(5)) * t173 + (qJD(3) * t143 + qJD(6) * t116) * t113 + t142;
t2 = (-pkin(5) * t166 + qJ(5) * qJD(3)) * t116 + (qJD(5) + (pkin(5) * qJD(3) - t164) * t113) * t114 + t157;
t1 = -pkin(5) * t149 - t162 + (-pkin(5) * t175 + t112 * t116) * qJD(3) - t7;
t20 = [0.2e1 * t6 * t59 + 0.2e1 * t7 * t60 + 0.2e1 * t1 * t61 + 0.2e1 * t4 * t62 + 0.2e1 * t2 * t63 + 0.2e1 * t5 * t64 + 0.2e1 * t3 * t45 + 0.2e1 * t8 * t46 + 0.2e1 * t15 * t27 + 0.2e1 * t23 * t28 + 0.2e1 * t19 * t29 + 0.2e1 * t24 * t30 + 0.2e1 * t31 * t25 + 0.2e1 * t32 * t26 + 0.2e1 * t33 * t16 + 0.2e1 * t22 * t18 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t200 + 0.2e1 * m(6) * (t23 * t4 + t24 * t5 + t33 * t8) + (t1 * t15 + t19 * t2 + t22 * t3) * t202 + 0.2e1 * m(5) * (t31 * t7 + t32 * t6) + (mrSges(4,1) * t200 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * Ifges(4,4) * t114 + 0.2e1 * t117 * t47 + t122 * t115 + (-t114 * t208 + t158) * t113) * qJD(3) + t121) * t114 + (mrSges(4,2) * t200 - 0.2e1 * t117 * t17 + (t10 + t11 - t14) * t115 + (t12 + t13 - t9) * t113 + ((-t158 - t179) * t115 + t122 * t113) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) - 0.2e1 * Ifges(4,4) * t116 + (-Ifges(6,4) - t194) * t173 + (-Ifges(5,6) + t208) * t176 + (-0.2e1 * m(5) * t117 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + t207) * t114) * qJD(3)) * t116; (-t16 - t17 - t18 - m(6) * t8 - m(7) * t3 + (t156 * t115 + t155 * t113 + m(5) * (-t113 * t31 + t115 * t32) + m(6) * (t113 * t24 - t115 * t23) + m(7) * (t113 * t15 + t115 * t19)) * qJD(3)) * t116 + ((t26 + t193) * t115 + (t27 + t192) * t113 + (t45 + t46 + t47) * qJD(3) + (-t113 * t156 + t115 * t155) * qJD(4) + m(5) * (-t113 * t7 + t115 * t6 - t166 * t31 - t168 * t32 - 0.2e1 * t150) + m(6) * (qJD(3) * t33 + t113 * t5 - t115 * t4 + t166 * t24 + t168 * t23) + m(7) * (qJD(3) * t22 + t1 * t113 + t115 * t2 + t15 * t166 - t168 * t19)) * t114; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t204) * t114 * t170; -pkin(3) * t17 + t68 * t16 + t44 * t18 + t21 * t45 + t22 * t50 + t81 * t27 + t82 * t29 + t3 * t74 + t33 * t51 + t40 * t46 + t66 * t61 - t65 * t63 + t8 * t72 + m(7) * (t1 * t81 + t15 * t66 - t19 * t65 + t2 * t82 + t21 * t22 + t3 * t44) + m(6) * (t33 * t40 + t68 * t8) + (t106 / 0.2e1 + t103 / 0.2e1 + t104 / 0.2e1 + t105 / 0.2e1 + (t117 * t203 - Ifges(4,5)) * qJD(3)) * t114 + (-t4 * mrSges(6,1) + t6 * mrSges(5,3) + t2 * mrSges(7,1) + t9 / 0.2e1 - t12 / 0.2e1 - t13 / 0.2e1 + t140 * t171 + (t24 * mrSges(6,1) - t31 * mrSges(5,3) + t15 * mrSges(7,1) + t35 / 0.2e1 + t36 / 0.2e1 - t39 / 0.2e1 - t180 / 0.2e1) * qJD(4) + (t26 - t28 + t191 * qJD(4) + m(5) * (-qJD(4) * t31 + t6) + m(6) * (qJD(4) * t24 - t4)) * pkin(8)) * t115 + (t5 * mrSges(6,1) - t7 * mrSges(5,3) + t1 * mrSges(7,1) + t10 / 0.2e1 + t11 / 0.2e1 - t14 / 0.2e1 - t139 * t171 + (t23 * mrSges(6,1) - t32 * mrSges(5,3) - t19 * mrSges(7,1) - t34 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1 - t179 / 0.2e1) * qJD(4) + ((-t59 + t62) * qJD(4) + m(5) * (-qJD(4) * t32 - t7) + m(6) * (qJD(4) * t23 + t5) + t192) * pkin(8)) * t113 + (-t117 * t52 + (t53 / 0.2e1 + t56 / 0.2e1 + t58 / 0.2e1) * t115 + (-t54 / 0.2e1 + t55 / 0.2e1 - t57 / 0.2e1) * t113 + (-t117 * mrSges(4,2) - Ifges(4,6) + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t115 + (Ifges(7,5) / 0.2e1 - Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t113) * qJD(3) + (t113 * t140 + t115 * t139) * qJD(4)) * t116; m(7) * (t113 * t66 - t115 * t65 + t166 * t81 - t168 * t82) * t114 + (-m(6) * t40 - m(7) * t21 - t50 - t51 - t52) * t116 + ((m(7) * t44 + t203 + t205 + t74) * t114 + (m(7) * (t113 * t81 + t115 * t82) - mrSges(4,2) + t204 * (mrSges(5,3) + t196)) * t116) * qJD(3) + (m(6) + m(5)) * t204 * pkin(8) * t170; 0.2e1 * t68 * t51 + 0.2e1 * t21 * t74 + 0.2e1 * t44 * t50 + (t21 * t44 - t65 * t82 + t66 * t81) * t202 - 0.2e1 * pkin(3) * t52 + 0.2e1 * t205 * t40 + (-t65 * t201 + t54 - t55 + t57) * t115 + (t66 * t201 + t53 + t56 + t58) * t113 + ((t201 * t81 - t127 + t132 + t80) * t115 + (-0.2e1 * mrSges(7,1) * t82 - t128 - t131 - t79) * t113) * qJD(4); t190 * qJD(5) + t193 * qJ(5) + m(6) * (-pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t23) + m(7) * (qJ(5) * t2 + qJD(5) * t19 - qJD(6) * t15 + t1 * t112) + t112 * t27 - qJD(6) * t61 - pkin(4) * t30 + t7 * mrSges(5,1) - t1 * mrSges(7,3) + t2 * mrSges(7,2) - t4 * mrSges(6,3) + t5 * mrSges(6,2) - t6 * mrSges(5,2) + (-Ifges(5,6) * t115 + t113 * t194) * t165 + (-t113 * t208 + t115 * t194) * t171 + t121; t161 * t145 + m(6) * (-pkin(4) * t145 + t188) + m(7) * (t112 * t145 + t188) - t159 * t151 + (t159 * t167 + t161 * t170 + m(6) * (-t152 - t154) + m(7) * (t112 * t170 - t152 - t162)) * t113; m(7) * (-qJ(5) * t65 + qJD(5) * t82 - qJD(6) * t81 + t112 * t66) + t106 + t103 + t104 + t105 - t66 * mrSges(7,3) - t65 * mrSges(7,2) - qJD(6) * t113 * mrSges(7,1) + qJD(5) * t123 + ((-mrSges(6,1) * pkin(4) + mrSges(7,1) * t112 - Ifges(6,4)) * t115 + (-qJ(5) * t196 - Ifges(5,6)) * t113 + (m(6) * t125 + t72 + t73) * pkin(8)) * qJD(4); 0.2e1 * qJD(6) * mrSges(7,3) + 0.2e1 * m(7) * (-qJD(6) * t112 + t172) + 0.2e1 * m(6) * t172 + 0.2e1 * t206 * qJD(5); t101 - t196 * t149 + m(6) * t5 + m(7) * t1 + (-mrSges(7,3) * t116 - t175 * t196) * qJD(3); (m(6) + m(7)) * (t113 * t170 + t145); m(7) * t66 + qJD(4) * t123; -m(7) * qJD(6); 0; m(7) * t2 + t29; (-t146 + t151) * m(7); -m(7) * t65 - mrSges(7,1) * t168; m(7) * qJD(5); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
