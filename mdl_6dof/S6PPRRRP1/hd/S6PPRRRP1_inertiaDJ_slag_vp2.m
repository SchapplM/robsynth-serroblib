% Calculate time derivative of joint inertia matrix for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:24
% EndTime: 2019-03-08 18:52:30
% DurationCPUTime: 2.33s
% Computational Cost: add. (2143->350), mult. (6522->517), div. (0->0), fcn. (6303->12), ass. (0->168)
t177 = Ifges(6,5) + Ifges(7,5);
t210 = -Ifges(6,3) - Ifges(7,3);
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t73 = -mrSges(7,1) * t106 + mrSges(7,2) * t103;
t90 = -pkin(5) * t106 - pkin(4);
t209 = m(7) * t90 + t73;
t107 = cos(qJ(4));
t148 = qJD(4) * t107;
t130 = t106 * t148;
t104 = sin(qJ(4));
t146 = qJD(5) * t104;
t112 = -t103 * t146 + t130;
t145 = qJD(5) * t106;
t111 = t103 * t148 + t104 * t145;
t208 = -m(6) * pkin(4) - mrSges(6,1) * t106 + mrSges(6,2) * t103 - mrSges(5,1);
t176 = Ifges(6,6) + Ifges(7,6);
t207 = t177 * t103 + t176 * t106;
t194 = m(6) / 0.2e1;
t143 = t194 + m(7) / 0.2e1;
t206 = 0.2e1 * t143;
t160 = Ifges(7,4) * t106;
t120 = -Ifges(7,2) * t103 + t160;
t162 = Ifges(6,4) * t106;
t121 = -Ifges(6,2) * t103 + t162;
t204 = (t120 + t121) * qJD(5);
t161 = Ifges(7,4) * t103;
t122 = Ifges(7,1) * t106 - t161;
t163 = Ifges(6,4) * t103;
t123 = Ifges(6,1) * t106 - t163;
t203 = (t122 + t123) * qJD(5);
t202 = -m(6) * pkin(10) - mrSges(6,3);
t201 = mrSges(7,3) - t202;
t200 = t208 + t209;
t199 = 0.2e1 * m(6);
t198 = 0.2e1 * m(7);
t197 = 2 * pkin(9);
t196 = -0.2e1 * mrSges(7,3);
t195 = m(5) / 0.2e1;
t193 = m(5) * pkin(3);
t190 = m(7) * pkin(5);
t101 = cos(pkin(7));
t102 = cos(pkin(6));
t99 = sin(pkin(6));
t159 = cos(pkin(12)) * t99;
t98 = sin(pkin(7));
t117 = t101 * t102 - t159 * t98;
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t139 = t101 * t159;
t116 = t102 * t98 + t139;
t181 = sin(pkin(12)) * t99;
t28 = t105 * t116 + t108 * t181;
t14 = t104 * t28 - t107 * t117;
t15 = t104 * t117 + t107 * t28;
t141 = t105 * t181;
t22 = (t108 * t116 - t141) * qJD(3);
t6 = qJD(4) * t15 + t22 * t104;
t189 = t14 * t6;
t186 = pkin(9) * t103;
t23 = t28 * qJD(3);
t157 = t108 * t98;
t27 = -t102 * t157 - t108 * t139 + t141;
t185 = t23 * t27;
t150 = qJD(3) * t108;
t134 = t98 * t150;
t158 = t105 * t98;
t49 = t101 * t104 + t107 * t158;
t30 = qJD(4) * t49 + t104 * t134;
t48 = -t101 * t107 + t104 * t158;
t184 = t48 * t30;
t183 = t6 * t104;
t7 = -qJD(4) * t14 + t22 * t107;
t182 = t7 * t107;
t179 = -mrSges(5,1) * t107 + mrSges(5,2) * t104 - mrSges(4,1);
t178 = -mrSges(6,2) - mrSges(7,2);
t175 = -qJ(6) - pkin(10);
t24 = mrSges(7,1) * t111 + mrSges(7,2) * t112;
t25 = mrSges(6,1) * t111 + mrSges(6,2) * t112;
t174 = t24 + t25;
t149 = qJD(4) * t104;
t34 = mrSges(7,1) * t149 - mrSges(7,3) * t112;
t35 = mrSges(6,1) * t149 - mrSges(6,3) * t112;
t173 = t34 + t35;
t36 = -mrSges(7,2) * t149 - mrSges(7,3) * t111;
t37 = -mrSges(6,2) * t149 - mrSges(6,3) * t111;
t172 = t36 + t37;
t44 = -Ifges(7,5) * t107 + t104 * t122;
t45 = -Ifges(6,5) * t107 + t104 * t123;
t171 = t44 + t45;
t69 = (pkin(4) * t104 - pkin(10) * t107) * qJD(4);
t71 = -pkin(4) * t107 - pkin(10) * t104 - pkin(3);
t170 = t103 * t69 + t145 * t71;
t169 = t106 * t69 + t149 * t186;
t52 = (mrSges(7,1) * t103 + mrSges(7,2) * t106) * t104;
t124 = mrSges(6,1) * t103 + mrSges(6,2) * t106;
t53 = t124 * t104;
t168 = t52 + t53;
t147 = qJD(5) * t103;
t57 = mrSges(7,1) * t147 + mrSges(7,2) * t145;
t58 = t124 * qJD(5);
t167 = t57 + t58;
t154 = t103 * t104;
t65 = mrSges(7,2) * t107 - mrSges(7,3) * t154;
t66 = mrSges(6,2) * t107 - mrSges(6,3) * t154;
t166 = t65 + t66;
t153 = t104 * t106;
t67 = -mrSges(7,1) * t107 - mrSges(7,3) * t153;
t68 = -mrSges(6,1) * t107 - mrSges(6,3) * t153;
t165 = t67 + t68;
t152 = t106 * t107;
t88 = pkin(9) * t152;
t40 = t103 * t71 + t88;
t29 = -qJD(4) * t48 + t107 * t134;
t156 = t29 * t107;
t155 = t30 * t104;
t151 = qJD(3) * t105;
t144 = qJD(6) * t106;
t140 = pkin(5) * t147;
t77 = Ifges(7,2) * t106 + t161;
t78 = Ifges(6,2) * t106 + t163;
t138 = -t77 / 0.2e1 - t78 / 0.2e1;
t79 = Ifges(7,1) * t103 + t160;
t80 = Ifges(6,1) * t103 + t162;
t137 = t79 / 0.2e1 + t80 / 0.2e1;
t136 = mrSges(7,1) + t190;
t135 = t98 * t151;
t129 = t176 * t103;
t128 = qJD(5) * t175;
t127 = -t130 * t177 + t149 * t210;
t126 = mrSges(6,1) + t136;
t125 = t30 * t14 + t48 * t6;
t9 = t103 * t27 + t106 * t15;
t8 = -t103 * t15 + t106 * t27;
t42 = -t107 * Ifges(7,6) + t104 * t120;
t43 = -t107 * Ifges(6,6) + t104 * t121;
t119 = t107 * t176 - t42 - t43;
t118 = t14 * t148 + t183;
t31 = -t103 * t49 - t106 * t157;
t115 = t103 * t157 - t106 * t49;
t114 = -t108 * t23 + t151 * t27;
t113 = t148 * t48 + t155;
t96 = Ifges(6,5) * t145;
t95 = Ifges(7,5) * t145;
t76 = t175 * t106;
t72 = t175 * t103;
t70 = (pkin(5) * t103 + pkin(9)) * t104;
t59 = (mrSges(5,1) * t104 + mrSges(5,2) * t107) * qJD(4);
t56 = t106 * t71;
t47 = -qJD(6) * t103 + t106 * t128;
t46 = t103 * t128 + t144;
t39 = -t107 * t186 + t56;
t38 = pkin(5) * t111 + pkin(9) * t148;
t33 = -qJ(6) * t154 + t40;
t26 = -qJ(6) * t153 + t56 + (-pkin(5) - t186) * t107;
t21 = -t80 * t146 + (Ifges(6,5) * t104 + t107 * t123) * qJD(4);
t20 = -t79 * t146 + (Ifges(7,5) * t104 + t107 * t122) * qJD(4);
t19 = -t78 * t146 + (Ifges(6,6) * t104 + t107 * t121) * qJD(4);
t18 = -t77 * t146 + (Ifges(7,6) * t104 + t107 * t120) * qJD(4);
t17 = -qJD(5) * t40 + t169;
t16 = (-t106 * t149 - t107 * t147) * pkin(9) + t170;
t13 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t153 + (-qJD(6) * t104 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t107) * t103 + t170;
t12 = qJD(5) * t115 - t103 * t29 + t106 * t135;
t11 = qJD(5) * t31 + t103 * t135 + t106 * t29;
t10 = -t104 * t144 + (pkin(5) * t104 - qJ(6) * t152) * qJD(4) + (-t88 + (qJ(6) * t104 - t71) * t103) * qJD(5) + t169;
t5 = qJD(5) * t8 + t103 * t23 + t106 * t7;
t4 = -qJD(5) * t9 - t103 * t7 + t106 * t23;
t1 = [0.2e1 * m(5) * (t15 * t7 + t185 + t189) + 0.2e1 * m(4) * (t22 * t28 + t185) + 0.2e1 * (t4 * t8 + t5 * t9 + t189) * t206; m(5) * (t29 * t15 + t49 * t7 + t125) + (t11 * t9 - t115 * t5 + t12 * t8 + t31 * t4 + t125) * t206 + 0.2e1 * (t114 * t195 + m(4) * (t105 * t22 + t150 * t28 + t114) / 0.2e1) * t98; 0.2e1 * m(5) * (-t105 * t150 * t98 ^ 2 + t29 * t49 + t184) + 0.4e1 * t143 * (-t11 * t115 + t12 * t31 + t184); -t22 * mrSges(4,2) + t27 * t59 + t172 * t9 + t173 * t8 + t168 * t6 + t166 * t5 + t165 * t4 + t174 * t14 + m(7) * (t10 * t8 + t13 * t9 + t14 * t38 + t26 * t4 + t33 * t5 + t6 * t70) + m(6) * (t16 * t9 + t17 * t8 + t39 * t4 + t40 * t5) + (t118 * t194 + (-t149 * t15 + t118 + t182) * t195) * t197 + (t183 + t182 + (-t104 * t15 + t107 * t14) * qJD(4)) * mrSges(5,3) + (t179 - t193) * t23; t174 * t48 - t172 * t115 + t173 * t31 + t168 * t30 + t165 * t12 + t166 * t11 + (-t108 * t59 + (-mrSges(4,2) * t108 + t105 * t179) * qJD(3)) * t98 + m(6) * (t11 * t40 - t115 * t16 + t12 * t39 + t17 * t31) + m(7) * (t10 * t31 + t11 * t33 - t115 * t13 + t12 * t26 + t30 * t70 + t38 * t48) - t135 * t193 + (t113 * t194 + (-t149 * t49 + t113 + t156) * t195) * t197 + (t155 + t156 + (-t104 * t49 + t107 * t48) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t59 + 0.2e1 * t10 * t67 + 0.2e1 * t13 * t65 + 0.2e1 * t16 * t66 + 0.2e1 * t17 * t68 + 0.2e1 * t70 * t24 + 0.2e1 * t26 * t34 + 0.2e1 * t33 * t36 + 0.2e1 * t39 * t35 + 0.2e1 * t40 * t37 + 0.2e1 * t38 * t52 + (t10 * t26 + t13 * t33 + t38 * t70) * t198 + (t16 * t40 + t17 * t39) * t199 + ((0.2e1 * Ifges(5,4) * t107 + t103 * t119 + t106 * t171 + t197 * t53) * qJD(4) + t127) * t107 + (t25 * t197 + (t20 + t21) * t106 + (-t18 - t19) * t103 + (t119 * t106 + (t107 * t177 - t171) * t103) * qJD(5) + ((t106 * t177 - 0.2e1 * Ifges(5,4) - t129) * t104 + ((pkin(9) ^ 2) * t199 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) + t210) * t107) * qJD(4)) * t104; -t7 * mrSges(5,2) + t167 * t14 + m(7) * (t14 * t140 + t4 * t72 + t46 * t9 + t47 * t8 - t5 * t76) + t200 * t6 + t201 * (-t4 * t103 + t5 * t106 + (-t103 * t9 - t106 * t8) * qJD(5)); -t29 * mrSges(5,2) + t167 * t48 + m(7) * (-t11 * t76 - t115 * t46 + t12 * t72 + t140 * t48 + t31 * t47) + t200 * t30 + t201 * (-t12 * t103 + t11 * t106 + (t103 * t115 - t106 * t31) * qJD(5)); t90 * t24 + t46 * t65 + t47 * t67 + t70 * t57 + t72 * t34 + t38 * t73 - t76 * t36 - pkin(4) * t25 + m(7) * (t10 * t72 - t13 * t76 + t26 * t47 + t33 * t46 + t38 * t90) + (-t95 / 0.2e1 - t96 / 0.2e1 + (pkin(9) * t208 + Ifges(5,5)) * qJD(4)) * t107 + (t20 / 0.2e1 + t21 / 0.2e1 - t10 * mrSges(7,3) - t17 * mrSges(6,3) + t138 * t148 + (-m(6) * t17 - t35) * pkin(10) + (pkin(5) * t52 - pkin(10) * t66 - t33 * mrSges(7,3) - t42 / 0.2e1 - t43 / 0.2e1 + (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t107 + t70 * t190 + t202 * t40) * qJD(5)) * t103 + (t18 / 0.2e1 + t19 / 0.2e1 + t13 * mrSges(7,3) + t16 * mrSges(6,3) + t137 * t148 + (-t39 * mrSges(6,3) - t26 * mrSges(7,3) + t44 / 0.2e1 + t45 / 0.2e1) * qJD(5) + (-qJD(5) * t68 + t37 + m(6) * (-qJD(5) * t39 + t16)) * pkin(10)) * t106 + (pkin(9) * t58 + (-t103 * t137 + t106 * t138) * qJD(5) - t204 * t103 / 0.2e1 + t203 * t106 / 0.2e1 + (-Ifges(5,6) + mrSges(5,2) * pkin(9) + t207 / 0.2e1) * qJD(4)) * t104; 0.2e1 * t90 * t57 - 0.2e1 * pkin(4) * t58 + (-t46 * t76 + t47 * t72) * t198 + (t47 * t196 + (0.2e1 * pkin(5) * t209 - t76 * t196 - t77 - t78) * qJD(5) + t203) * t103 + (0.2e1 * t46 * mrSges(7,3) + (t196 * t72 + t79 + t80) * qJD(5) + t204) * t106; t126 * t4 + t178 * t5; t11 * t178 + t12 * t126; mrSges(6,1) * t17 + mrSges(7,1) * t10 - mrSges(6,2) * t16 - mrSges(7,2) * t13 - t129 * t148 + (m(7) * t10 + t34) * pkin(5) - t207 * t146 - t127; -mrSges(7,2) * t46 + t95 + t96 + t136 * t47 + ((-mrSges(6,1) * pkin(10) - mrSges(7,3) * pkin(5)) * t106 + (mrSges(6,2) * pkin(10) - t176) * t103) * qJD(5); 0; m(7) * t6; m(7) * t30; m(7) * t38 + t24; m(7) * t140 + t57; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
