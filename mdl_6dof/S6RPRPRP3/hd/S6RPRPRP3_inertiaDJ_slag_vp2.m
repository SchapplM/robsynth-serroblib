% Calculate time derivative of joint inertia matrix for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:55
% EndTime: 2019-03-09 03:07:59
% DurationCPUTime: 2.50s
% Computational Cost: add. (2412->350), mult. (5751->506), div. (0->0), fcn. (4844->8), ass. (0->147)
t197 = Ifges(7,4) + Ifges(6,5);
t195 = Ifges(7,6) - Ifges(6,6);
t196 = -Ifges(7,2) - Ifges(6,3);
t123 = cos(pkin(10));
t127 = cos(qJ(5));
t122 = sin(pkin(10));
t125 = sin(qJ(5));
t161 = t122 * t125;
t101 = -t127 * t123 + t161;
t102 = t122 * t127 + t123 * t125;
t126 = sin(qJ(3));
t152 = qJD(5) * t126;
t128 = cos(qJ(3));
t154 = qJD(3) * t128;
t47 = -t101 * t154 - t102 * t152;
t151 = qJD(5) * t127;
t158 = t123 * t126;
t48 = t102 * t154 + t151 * t158 - t152 * t161;
t14 = t48 * mrSges(7,1) - t47 * mrSges(7,3);
t15 = t48 * mrSges(6,1) + t47 * mrSges(6,2);
t134 = -t14 - t15;
t145 = t122 * t154;
t170 = mrSges(5,2) * t123;
t84 = mrSges(5,1) * t145 + t154 * t170;
t194 = -t84 + t134;
t147 = -cos(pkin(9)) * pkin(1) - pkin(2);
t193 = 0.2e1 * t147;
t192 = m(7) + m(6);
t191 = mrSges(6,3) + mrSges(7,2);
t116 = sin(pkin(9)) * pkin(1) + pkin(7);
t159 = t122 * t128;
t98 = -pkin(3) * t128 - t126 * qJ(4) + t147;
t83 = t123 * t98;
t61 = -t116 * t159 + t83;
t157 = t123 * t128;
t62 = t116 * t157 + t122 * t98;
t190 = -t122 * t61 + t123 * t62;
t155 = qJD(3) * t126;
t146 = t116 * t155;
t149 = t126 * qJD(4);
t180 = pkin(3) * t126;
t92 = -t149 + (-qJ(4) * t128 + t180) * qJD(3);
t59 = t122 * t146 + t123 * t92;
t77 = t122 * t92;
t60 = -t123 * t146 + t77;
t189 = -t122 * t59 + t123 * t60;
t93 = t101 * qJD(5);
t94 = t102 * qJD(5);
t188 = t195 * t94 - t197 * t93;
t187 = m(7) * qJ(6) + mrSges(7,3);
t45 = -pkin(8) * t158 + t83 + (-t116 * t122 - pkin(4)) * t128;
t160 = t122 * t126;
t51 = -pkin(8) * t160 + t62;
t176 = t125 * t45 + t127 * t51;
t27 = (pkin(4) * t126 - pkin(8) * t157) * qJD(3) + t59;
t110 = t126 * t116;
t37 = t77 + (-pkin(8) * t159 - t123 * t110) * qJD(3);
t4 = -qJD(5) * t176 - t125 * t37 + t127 * t27;
t186 = 2 * m(5);
t185 = 2 * m(6);
t184 = 0.2e1 * m(7);
t121 = t123 ^ 2;
t183 = 0.2e1 * t116;
t182 = t123 / 0.2e1;
t179 = pkin(8) + qJ(4);
t23 = -mrSges(7,2) * t48 + mrSges(7,3) * t155;
t26 = -mrSges(6,2) * t155 - mrSges(6,3) * t48;
t178 = t23 + t26;
t24 = mrSges(6,1) * t155 - mrSges(6,3) * t47;
t25 = -mrSges(7,1) * t155 + t47 * mrSges(7,2);
t177 = -t24 + t25;
t52 = t94 * mrSges(7,1) + t93 * mrSges(7,3);
t53 = t94 * mrSges(6,1) - t93 * mrSges(6,2);
t175 = -t52 - t53;
t80 = t102 * t126;
t71 = -t80 * mrSges(7,2) - mrSges(7,3) * t128;
t72 = mrSges(6,2) * t128 - t80 * mrSges(6,3);
t174 = t71 + t72;
t81 = t101 * t126;
t73 = -mrSges(6,1) * t128 + t81 * mrSges(6,3);
t74 = mrSges(7,1) * t128 - t81 * mrSges(7,2);
t173 = -t73 + t74;
t169 = Ifges(5,4) * t122;
t168 = Ifges(5,4) * t123;
t167 = t122 * Ifges(5,2);
t162 = -mrSges(5,1) * t123 + mrSges(5,2) * t122 - mrSges(4,1);
t106 = t116 * t154;
t79 = pkin(4) * t145 + t106;
t91 = pkin(4) * t160 + t110;
t156 = t122 ^ 2 + t121;
t153 = qJD(5) * t125;
t150 = t125 * qJD(4);
t148 = t127 * qJD(4);
t117 = -pkin(4) * t123 - pkin(3);
t144 = t126 * t154;
t109 = t179 * t123;
t141 = qJD(5) * t179;
t31 = t123 * t148 - t109 * t153 + (-t127 * t141 - t150) * t122;
t32 = t123 * t150 + t109 * t151 + (-t125 * t141 + t148) * t122;
t142 = t179 * t122;
t69 = t109 * t125 + t127 * t142;
t70 = t127 * t109 - t125 * t142;
t143 = t70 * t31 + t32 * t69;
t140 = t156 * mrSges(5,3);
t139 = t156 * qJ(4);
t138 = -Ifges(5,5) * t123 + Ifges(5,6) * t122;
t12 = -t125 * t51 + t127 * t45;
t135 = t196 * t155 - t195 * t48 - t197 * t47;
t3 = t125 * t27 + t127 * t37 + t45 * t151 - t51 * t153;
t133 = -t31 * t81 + t32 * t80 + t70 * t47 + t69 * t48;
t131 = -pkin(5) * t48 + qJ(6) * t47 - qJD(6) * t81;
t104 = -mrSges(5,1) * t128 - mrSges(5,3) * t158;
t103 = mrSges(5,2) * t128 - mrSges(5,3) * t160;
t97 = (mrSges(5,1) * t126 - mrSges(5,3) * t157) * qJD(3);
t96 = (-mrSges(5,2) * t126 - mrSges(5,3) * t159) * qJD(3);
t95 = (mrSges(5,1) * t122 + t170) * t126;
t76 = (t126 * Ifges(5,5) + (t123 * Ifges(5,1) - t169) * t128) * qJD(3);
t75 = (t126 * Ifges(5,6) + (-t167 + t168) * t128) * qJD(3);
t68 = Ifges(6,1) * t102 - Ifges(6,4) * t101;
t67 = Ifges(7,1) * t102 + Ifges(7,5) * t101;
t66 = Ifges(6,4) * t102 - Ifges(6,2) * t101;
t65 = Ifges(7,5) * t102 + Ifges(7,3) * t101;
t64 = mrSges(6,1) * t101 + mrSges(6,2) * t102;
t63 = mrSges(7,1) * t101 - mrSges(7,3) * t102;
t58 = pkin(5) * t101 - qJ(6) * t102 + t117;
t57 = -Ifges(6,1) * t93 - Ifges(6,4) * t94;
t56 = -Ifges(7,1) * t93 + Ifges(7,5) * t94;
t55 = -Ifges(6,4) * t93 - Ifges(6,2) * t94;
t54 = -Ifges(7,5) * t93 + Ifges(7,3) * t94;
t50 = mrSges(6,1) * t80 - mrSges(6,2) * t81;
t49 = mrSges(7,1) * t80 + mrSges(7,3) * t81;
t36 = -Ifges(6,1) * t81 - Ifges(6,4) * t80 - Ifges(6,5) * t128;
t35 = -Ifges(7,1) * t81 - Ifges(7,4) * t128 + Ifges(7,5) * t80;
t34 = -Ifges(6,4) * t81 - Ifges(6,2) * t80 - Ifges(6,6) * t128;
t33 = -Ifges(7,5) * t81 - Ifges(7,6) * t128 + Ifges(7,3) * t80;
t22 = pkin(5) * t94 + qJ(6) * t93 - qJD(6) * t102;
t20 = pkin(5) * t80 + qJ(6) * t81 + t91;
t11 = pkin(5) * t128 - t12;
t10 = -qJ(6) * t128 + t176;
t9 = Ifges(6,1) * t47 - Ifges(6,4) * t48 + Ifges(6,5) * t155;
t8 = Ifges(7,1) * t47 + Ifges(7,4) * t155 + Ifges(7,5) * t48;
t7 = Ifges(6,4) * t47 - Ifges(6,2) * t48 + Ifges(6,6) * t155;
t6 = Ifges(7,5) * t47 + Ifges(7,6) * t155 + Ifges(7,3) * t48;
t5 = -t131 + t79;
t2 = -pkin(5) * t155 - t4;
t1 = qJ(6) * t155 - qJD(6) * t128 + t3;
t13 = [t135 * t128 + (-t122 * t75 + t123 * t76 + t84 * t183) * t126 + 0.2e1 * t176 * t26 + (t12 * t4 + t176 * t3 + t79 * t91) * t185 - (t8 + t9) * t81 + (t6 - t7) * t80 + ((mrSges(4,1) * t193 + (-(2 * Ifges(4,4)) - t138) * t126 - t197 * t81 + t195 * t80) * t126 + (t95 * t183 + mrSges(4,2) * t193 + 0.2e1 * (Ifges(4,4) + t138) * t128 + (-(2 * Ifges(5,3)) + t116 ^ 2 * t186 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + Ifges(5,1) * t121 + (t167 - 0.2e1 * t168) * t122 + t196) * t126) * t128) * qJD(3) + 0.2e1 * t10 * t23 + 0.2e1 * t12 * t24 + 0.2e1 * t11 * t25 + 0.2e1 * t20 * t14 + (t1 * t10 + t11 * t2 + t20 * t5) * t184 + (t61 * t59 + t62 * t60) * t186 + (t33 - t34) * t48 + (t35 + t36) * t47 + 0.2e1 * t5 * t49 + 0.2e1 * t1 * t71 + 0.2e1 * t3 * t72 + 0.2e1 * t4 * t73 + 0.2e1 * t2 * t74 + 0.2e1 * t79 * t50 + 0.2e1 * t91 * t15 + 0.2e1 * t62 * t96 + 0.2e1 * t61 * t97 + 0.2e1 * t60 * t103 + 0.2e1 * t59 * t104; -t178 * t81 + t177 * t80 + t173 * t48 + t174 * t47 + (-t122 * t97 + t123 * t96) * t126 + t194 * t128 + ((t103 * t123 - t104 * t122) * t128 + (t49 + t50 + t95) * t126) * qJD(3) + m(7) * (-t1 * t81 + t10 * t47 + t11 * t48 - t128 * t5 + t20 * t155 + t2 * t80) + m(6) * (-t12 * t48 - t128 * t79 + t91 * t155 + t176 * t47 - t3 * t81 - t4 * t80) + m(5) * (t189 * t126 + (t116 * t126 ^ 2 + (-t116 * t128 + t190) * t128) * qJD(3)); 0.2e1 * m(5) * (-0.1e1 + t156) * t144 + 0.2e1 * t192 * (-t81 * t47 + t80 * t48 - t144); (t8 / 0.2e1 + t9 / 0.2e1) * t102 + (t6 / 0.2e1 - t7 / 0.2e1) * t101 + ((Ifges(5,5) * t122 / 0.2e1 + Ifges(5,6) * t182 - Ifges(4,6) + t116 * mrSges(4,2) + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t102 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t101) * t126 + (Ifges(4,5) + (Ifges(5,1) * t122 + t168) * t182 - t122 * (Ifges(5,2) * t123 + t169) / 0.2e1 + (-m(5) * pkin(3) + t162) * t116) * t128) * qJD(3) + t177 * t69 + t178 * t70 + t173 * t32 + t174 * t31 + (-t101 * t3 - t102 * t4 + t12 * t93 - t176 * t94) * mrSges(6,3) + (-t1 * t101 - t10 * t94 + t102 * t2 - t11 * t93) * mrSges(7,2) + m(7) * (t1 * t70 + t10 * t31 + t11 * t32 + t2 * t69 + t20 * t22 + t5 * t58) + m(6) * (t117 * t79 - t12 * t32 + t176 * t31 + t3 * t70 - t4 * t69) - (t56 / 0.2e1 + t57 / 0.2e1) * t81 + (t54 / 0.2e1 - t55 / 0.2e1) * t80 + (t65 / 0.2e1 - t66 / 0.2e1) * t48 + (t67 / 0.2e1 + t68 / 0.2e1) * t47 + (qJ(4) * t96 + qJD(4) * t103 + t60 * mrSges(5,3) + t75 / 0.2e1) * t123 + (-qJ(4) * t97 - qJD(4) * t104 - t59 * mrSges(5,3) + t76 / 0.2e1) * t122 + (t33 / 0.2e1 - t34 / 0.2e1) * t94 - (t35 / 0.2e1 + t36 / 0.2e1) * t93 + t22 * t49 - t188 * t128 / 0.2e1 + t20 * t52 + m(5) * (t189 * qJ(4) + t190 * qJD(4)) + t58 * t14 + t5 * t63 + t79 * t64 - pkin(3) * t84 + t91 * t53 + t117 * t15; t175 * t128 + ((-mrSges(4,2) + t140) * t128 + (t63 + t64 + t162) * t126) * qJD(3) + m(7) * (-t128 * t22 + t58 * t155 + t133) + m(6) * (t117 * t155 + t133) + m(5) * (t156 * t149 + (t128 * t139 - t180) * qJD(3)) + t191 * (-t101 * t47 + t102 * t48 - t80 * t93 + t81 * t94); 0.2e1 * t117 * t53 + 0.2e1 * t22 * t63 + 0.2e1 * t58 * t52 + (-t66 + t65) * t94 - (t67 + t68) * t93 + (t56 + t57) * t102 + (-t55 + t54) * t101 + t143 * t185 + (t22 * t58 + t143) * t184 + (t139 * t186 + 0.2e1 * t140) * qJD(4) + 0.2e1 * t191 * (-t101 * t31 + t102 * t32 - t69 * t93 - t70 * t94); m(5) * t106 + m(6) * t79 + m(7) * t5 - t194; (m(5) + t192) * t155; m(7) * t22 - t175; 0; m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t10) + t1 * mrSges(7,3) + qJD(6) * t71 + qJ(6) * t23 + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) - pkin(5) * t25 - t135; m(7) * t131 + t134; m(7) * qJD(6) * t70 + (pkin(5) * t93 - qJ(6) * t94 - qJD(6) * t101) * mrSges(7,2) + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t32 + (-mrSges(6,2) + t187) * t31 + t188; 0; 0.2e1 * t187 * qJD(6); m(7) * t2 + t25; m(7) * t48; m(7) * t32 - t93 * mrSges(7,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
