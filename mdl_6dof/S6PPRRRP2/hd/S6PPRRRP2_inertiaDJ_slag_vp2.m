% Calculate time derivative of joint inertia matrix for
% S6PPRRRP2
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:06
% EndTime: 2019-03-08 18:56:11
% DurationCPUTime: 2.57s
% Computational Cost: add. (2179->351), mult. (6636->516), div. (0->0), fcn. (6435->12), ass. (0->169)
t187 = Ifges(7,4) + Ifges(6,5);
t214 = -Ifges(7,2) - Ifges(6,3);
t108 = sin(qJ(4));
t107 = sin(qJ(5));
t111 = cos(qJ(4));
t155 = qJD(4) * t111;
t140 = t107 * t155;
t110 = cos(qJ(5));
t152 = qJD(5) * t110;
t119 = t108 * t152 + t140;
t101 = sin(pkin(12));
t103 = sin(pkin(6));
t106 = cos(pkin(6));
t109 = sin(qJ(3));
t112 = cos(qJ(3));
t104 = cos(pkin(12));
t105 = cos(pkin(7));
t161 = t104 * t105;
t102 = sin(pkin(7));
t162 = t102 * t112;
t213 = (-t101 * t109 + t112 * t161) * t103 + t106 * t162;
t212 = m(6) + m(7);
t121 = -t102 * t103 * t104 + t105 * t106;
t163 = t102 * t109;
t37 = t106 * t163 + (t101 * t112 + t109 * t161) * t103;
t20 = t108 * t121 + t37 * t111;
t31 = t213 * qJD(3);
t10 = qJD(4) * t20 + t31 * t108;
t19 = t37 * t108 - t111 * t121;
t157 = qJD(3) * t112;
t142 = t102 * t157;
t58 = t105 * t108 + t111 * t163;
t39 = qJD(4) * t58 + t108 * t142;
t57 = -t111 * t105 + t108 * t163;
t211 = t57 * t10 + t39 * t19;
t210 = mrSges(6,3) + mrSges(7,2);
t175 = Ifges(7,5) * t107;
t132 = Ifges(7,1) * t110 + t175;
t177 = Ifges(6,4) * t107;
t133 = Ifges(6,1) * t110 - t177;
t209 = (t132 + t133) * qJD(5);
t208 = m(7) * qJD(6);
t206 = m(7) * qJ(6) + mrSges(7,3);
t151 = qJD(5) * t111;
t154 = qJD(5) * t107;
t156 = qJD(4) * t108;
t76 = (pkin(4) * t108 - pkin(10) * t111) * qJD(4);
t80 = -pkin(4) * t111 - pkin(10) * t108 - pkin(3);
t24 = pkin(9) * (t107 * t156 - t110 * t151) + t110 * t76 - t154 * t80;
t197 = m(6) / 0.2e1;
t205 = 0.2e1 * pkin(10) * (m(7) / 0.2e1 + t197);
t129 = pkin(5) * t110 + qJ(6) * t107;
t150 = qJD(6) * t110;
t204 = qJD(5) * t129 - t150;
t203 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t202 = -mrSges(6,2) + t206;
t201 = 0.2e1 * m(6);
t200 = 0.2e1 * pkin(9);
t198 = m(5) / 0.2e1;
t195 = m(5) * pkin(3);
t193 = pkin(9) * t111;
t7 = t19 * t10;
t32 = t37 * qJD(3);
t192 = t32 * t213;
t11 = -qJD(4) * t19 + t31 * t111;
t13 = -t107 * t213 + t110 * t20;
t4 = qJD(5) * t13 + t107 * t11 - t32 * t110;
t191 = t4 * t107;
t12 = t107 * t20 + t110 * t213;
t5 = -qJD(5) * t12 + t107 * t32 + t11 * t110;
t190 = t5 * t110;
t25 = t57 * t39;
t189 = qJD(4) / 0.2e1;
t188 = -mrSges(5,1) * t111 + mrSges(5,2) * t108 - mrSges(4,1);
t139 = t110 * t155;
t153 = qJD(5) * t108;
t120 = -t107 * t153 + t139;
t34 = mrSges(7,1) * t119 - mrSges(7,3) * t120;
t35 = mrSges(6,1) * t119 + mrSges(6,2) * t120;
t186 = t34 + t35;
t44 = mrSges(6,1) * t156 - mrSges(6,3) * t120;
t45 = mrSges(7,2) * t139 + (-mrSges(7,1) * qJD(4) - mrSges(7,2) * t154) * t108;
t185 = -t44 + t45;
t46 = -mrSges(6,2) * t156 - mrSges(6,3) * t119;
t47 = -mrSges(7,2) * t119 + mrSges(7,3) * t156;
t184 = t46 + t47;
t53 = -t111 * Ifges(7,4) + t108 * t132;
t54 = -t111 * Ifges(6,5) + t108 * t133;
t183 = t53 + t54;
t134 = mrSges(7,1) * t107 - mrSges(7,3) * t110;
t60 = t134 * t108;
t135 = mrSges(6,1) * t107 + mrSges(6,2) * t110;
t61 = t135 * t108;
t182 = t60 + t61;
t64 = t134 * qJD(5);
t65 = t135 * qJD(5);
t181 = t64 + t65;
t160 = t107 * t108;
t72 = mrSges(6,2) * t111 - mrSges(6,3) * t160;
t75 = -mrSges(7,2) * t160 - mrSges(7,3) * t111;
t180 = t72 + t75;
t159 = t108 * t110;
t73 = -mrSges(6,1) * t111 - mrSges(6,3) * t159;
t74 = mrSges(7,1) * t111 + mrSges(7,2) * t159;
t179 = -t73 + t74;
t82 = -t110 * mrSges(6,1) + t107 * mrSges(6,2);
t178 = t82 - mrSges(5,1);
t49 = t107 * t80 + t110 * t193;
t176 = Ifges(6,4) * t110;
t174 = Ifges(7,5) * t110;
t173 = Ifges(6,6) * t110;
t172 = Ifges(6,6) * t111;
t171 = t10 * t108;
t170 = t11 * t111;
t169 = t110 * t80;
t158 = qJD(3) * t109;
t143 = t102 * t158;
t38 = -qJD(4) * t57 + t111 * t142;
t40 = t107 * t58 + t110 * t162;
t16 = -qJD(5) * t40 + t107 * t143 + t110 * t38;
t168 = t16 * t110;
t146 = t107 * t162;
t17 = -qJD(5) * t146 + t107 * t38 - t110 * t143 + t152 * t58;
t167 = t17 * t107;
t166 = t38 * t111;
t165 = t39 * t108;
t81 = -t110 * mrSges(7,1) - t107 * mrSges(7,3);
t149 = t81 + t178;
t84 = -Ifges(7,3) * t110 + t175;
t85 = Ifges(6,2) * t110 + t177;
t148 = t84 / 0.2e1 - t85 / 0.2e1;
t86 = Ifges(7,1) * t107 - t174;
t87 = Ifges(6,1) * t107 + t176;
t147 = t86 / 0.2e1 + t87 / 0.2e1;
t130 = Ifges(7,3) * t107 + t174;
t51 = -Ifges(7,6) * t111 + t108 * t130;
t131 = -Ifges(6,2) * t107 + t176;
t52 = t108 * t131 - t172;
t136 = t51 - t52 + t172;
t128 = pkin(5) * t107 - qJ(6) * t110;
t127 = -t119 * Ifges(7,6) - t187 * t139 + t214 * t156;
t126 = pkin(9) + t128;
t125 = -t112 * t32 - t158 * t213;
t124 = t155 * t19 + t171;
t123 = t155 * t57 + t165;
t23 = t107 * t76 + t80 * t152 + (-t107 * t151 - t110 * t156) * pkin(9);
t100 = Ifges(7,4) * t152;
t99 = Ifges(6,5) * t152;
t97 = Ifges(7,6) * t154;
t77 = -pkin(4) - t129;
t68 = t131 * qJD(5);
t67 = t130 * qJD(5);
t66 = (mrSges(5,1) * t108 + mrSges(5,2) * t111) * qJD(4);
t56 = qJD(5) * t128 - qJD(6) * t107;
t55 = t126 * t108;
t48 = -t107 * t193 + t169;
t43 = -t169 + (pkin(9) * t107 + pkin(5)) * t111;
t42 = -qJ(6) * t111 + t49;
t41 = t110 * t58 - t146;
t29 = -t87 * t153 + (Ifges(6,5) * t108 + t111 * t133) * qJD(4);
t28 = -t86 * t153 + (Ifges(7,4) * t108 + t111 * t132) * qJD(4);
t27 = -t85 * t153 + (Ifges(6,6) * t108 + t111 * t131) * qJD(4);
t26 = -t84 * t153 + (Ifges(7,6) * t108 + t111 * t130) * qJD(4);
t22 = t204 * t108 + t126 * t155;
t21 = -pkin(5) * t156 - t24;
t18 = qJ(6) * t156 - qJD(6) * t111 + t23;
t15 = pkin(10) * t168;
t3 = pkin(10) * t190;
t1 = [0.2e1 * m(5) * (t11 * t20 - t192 + t7) + 0.2e1 * m(4) * (t31 * t37 - t192) + 0.2e1 * t212 * (t12 * t4 + t13 * t5 + t7); m(5) * (t58 * t11 + t38 * t20 + t211) + 0.2e1 * (t125 * t198 + m(4) * (t109 * t31 + t157 * t37 + t125) / 0.2e1) * t102 + t212 * (t12 * t17 + t16 * t13 + t4 * t40 + t41 * t5 + t211); 0.2e1 * m(5) * (-t102 ^ 2 * t109 * t157 + t58 * t38 + t25) + 0.2e1 * t212 * (t41 * t16 + t17 * t40 + t25); -t31 * mrSges(4,2) - t213 * t66 + t180 * t5 + t179 * t4 + t186 * t19 + t184 * t13 + t185 * t12 + t182 * t10 + m(7) * (t10 * t55 + t12 * t21 + t13 * t18 + t19 * t22 + t4 * t43 + t42 * t5) + m(6) * (-t12 * t24 + t13 * t23 - t4 * t48 + t49 * t5) + (t124 * t197 + (-t156 * t20 + t124 + t170) * t198) * t200 + (t171 + t170 + (-t108 * t20 + t111 * t19) * qJD(4)) * mrSges(5,3) + (t188 - t195) * t32; t186 * t57 + t184 * t41 + t185 * t40 + t182 * t39 + t179 * t17 + t180 * t16 + (-t112 * t66 + (-mrSges(4,2) * t112 + t109 * t188) * qJD(3)) * t102 + m(6) * (t16 * t49 - t17 * t48 + t23 * t41 - t24 * t40) + m(7) * (t16 * t42 + t17 * t43 + t18 * t41 + t21 * t40 + t22 * t57 + t39 * t55) - t143 * t195 + (t123 * t197 + (-t156 * t58 + t123 + t166) * t198) * t200 + (t165 + t166 + (-t108 * t58 + t111 * t57) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t66 + 0.2e1 * t18 * t75 + 0.2e1 * t21 * t74 + 0.2e1 * t22 * t60 + 0.2e1 * t23 * t72 + 0.2e1 * t24 * t73 + 0.2e1 * t55 * t34 + 0.2e1 * t42 * t47 + 0.2e1 * t43 * t45 + 0.2e1 * t48 * t44 + 0.2e1 * t49 * t46 + 0.2e1 * m(7) * (t18 * t42 + t21 * t43 + t22 * t55) + (t23 * t49 + t24 * t48) * t201 + ((0.2e1 * Ifges(5,4) * t111 + t107 * t136 + t110 * t183 + t200 * t61) * qJD(4) + t127) * t111 + (t35 * t200 + (t28 + t29) * t110 + (t26 - t27) * t107 + (t136 * t110 + (t111 * t187 - t183) * t107) * qJD(5) + ((-0.2e1 * Ifges(5,4) + t187 * t110 + (-Ifges(6,6) + Ifges(7,6)) * t107) * t108 + (pkin(9) ^ 2 * t201 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) + t214) * t111) * qJD(4)) * t108; -t11 * mrSges(5,2) + t181 * t19 + t149 * t10 + m(7) * (t10 * t77 + t19 * t56 + t3) + m(6) * (-pkin(4) * t10 + t3) + (t12 * t152 - t13 * t154 + t191) * t205 + t210 * (t191 + t190 + (-t107 * t13 + t110 * t12) * qJD(5)); -t38 * mrSges(5,2) + t181 * t57 + t149 * t39 + m(6) * (-pkin(4) * t39 + t15) + m(7) * (t39 * t77 + t56 * t57 + t15) + (t152 * t40 - t154 * t41 + t167) * t205 + t210 * (t167 + t168 + (-t107 * t41 + t110 * t40) * qJD(5)); m(7) * (t22 * t77 + t55 * t56) - pkin(4) * t35 + t56 * t60 + t55 * t64 + t77 * t34 + t22 * t81 + (-t100 / 0.2e1 - t97 / 0.2e1 - t99 / 0.2e1 + (Ifges(5,5) + (-m(6) * pkin(4) + t178) * pkin(9)) * qJD(4)) * t111 + (t21 * mrSges(7,2) - t24 * mrSges(6,3) + t28 / 0.2e1 + t29 / 0.2e1 + t148 * t155 + (-t42 * mrSges(7,2) - t49 * mrSges(6,3) + t51 / 0.2e1 - t52 / 0.2e1 + t172 / 0.2e1) * qJD(5) + (-t180 * qJD(5) + m(7) * (-t42 * qJD(5) + t21) + m(6) * (-t49 * qJD(5) - t24) + t185) * pkin(10)) * t107 + (t23 * mrSges(6,3) + t18 * mrSges(7,2) - t26 / 0.2e1 + t27 / 0.2e1 + t147 * t155 + (t43 * mrSges(7,2) - t48 * mrSges(6,3) + t53 / 0.2e1 + t54 / 0.2e1) * qJD(5) + (t179 * qJD(5) + m(7) * (t43 * qJD(5) + t18) + m(6) * (-t48 * qJD(5) + t23) + t184) * pkin(10)) * t110 + (-Ifges(5,6) * qJD(4) + t173 * t189 + (qJD(4) * mrSges(5,2) + t65) * pkin(9) + (t67 / 0.2e1 - t68 / 0.2e1 - t147 * qJD(5) + t187 * t189) * t107 + (t209 / 0.2e1 - Ifges(7,6) * t189 + qJD(5) * t148) * t110) * t108; -0.2e1 * pkin(4) * t65 + 0.2e1 * t64 * t77 + 0.2e1 * (m(7) * t77 + t81) * t56 + (-t67 + t68) * t110 + t209 * t107 + ((t86 + t87) * t110 + (t84 - t85) * t107) * qJD(5); t13 * t208 + t202 * t5 + t203 * t4; t202 * t16 + t203 * t17 + t41 * t208; -Ifges(6,6) * t140 + t18 * mrSges(7,3) + qJD(6) * t75 + qJ(6) * t47 + m(7) * (-pkin(5) * t21 + qJ(6) * t18 + qJD(6) * t42) - t23 * mrSges(6,2) + t24 * mrSges(6,1) - t21 * mrSges(7,1) - pkin(5) * t45 + (-t107 * t187 - t173) * t153 - t127; -Ifges(6,6) * t154 + t100 + t97 + t99 - t204 * mrSges(7,2) + (m(7) * t150 + (-m(7) * t129 + t81 + t82) * qJD(5)) * pkin(10); 0.2e1 * t206 * qJD(6); m(7) * t4; m(7) * t17; m(7) * t21 + t45; (m(7) * pkin(10) + mrSges(7,2)) * t152; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
