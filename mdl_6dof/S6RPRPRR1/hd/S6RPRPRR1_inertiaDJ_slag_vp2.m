% Calculate time derivative of joint inertia matrix for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:54
% EndTime: 2019-03-09 03:33:57
% DurationCPUTime: 2.00s
% Computational Cost: add. (4871->273), mult. (10028->408), div. (0->0), fcn. (9981->10), ass. (0->127)
t176 = -2 * mrSges(6,3);
t105 = sin(qJ(6));
t108 = cos(qJ(6));
t155 = mrSges(7,1) * t108;
t91 = mrSges(7,2) * t105 - t155;
t181 = -mrSges(6,1) + t91;
t141 = t105 ^ 2 + t108 ^ 2;
t106 = sin(qJ(5));
t169 = cos(qJ(5));
t102 = sin(pkin(11));
t107 = sin(qJ(3));
t103 = cos(pkin(11));
t109 = cos(qJ(3));
t143 = t103 * t109;
t85 = -t102 * t107 + t143;
t86 = t102 * t109 + t103 * t107;
t117 = -t106 * t86 + t169 * t85;
t78 = t86 * qJD(3);
t79 = t85 * qJD(3);
t41 = qJD(5) * t117 - t106 * t78 + t169 * t79;
t131 = t141 * t41;
t139 = qJD(6) * t108;
t61 = t106 * t85 + t169 * t86;
t183 = t105 * t41 + t139 * t61;
t167 = pkin(3) * t102;
t97 = pkin(3) * t103 + pkin(4);
t75 = -t106 * t167 + t169 * t97;
t66 = t75 * qJD(5);
t182 = t66 * mrSges(6,2);
t96 = sin(pkin(10)) * pkin(1) + pkin(7);
t146 = qJ(4) + t96;
t82 = t146 * t107;
t83 = t146 * t109;
t58 = -t102 * t83 - t103 * t82;
t118 = -pkin(8) * t86 + t58;
t69 = t102 * t82;
t59 = t103 * t83 - t69;
t50 = pkin(8) * t85 + t59;
t24 = t106 * t118 + t169 * t50;
t135 = -cos(pkin(10)) * pkin(1) - pkin(2);
t90 = -pkin(3) * t109 + t135;
t64 = -pkin(4) * t85 + t90;
t25 = -pkin(5) * t117 - pkin(9) * t61 + t64;
t14 = t105 * t25 + t108 * t24;
t142 = t14 * qJD(6);
t48 = -t86 * qJD(4) + (-t143 * t146 + t69) * qJD(3);
t110 = -t79 * pkin(8) + t48;
t179 = -t106 * t50 + t118 * t169;
t129 = qJD(3) * t146;
t49 = t103 * (qJD(4) * t109 - t107 * t129) + t102 * (-t107 * qJD(4) - t109 * t129);
t46 = -pkin(8) * t78 + t49;
t10 = qJD(5) * t179 + t106 * t110 + t169 * t46;
t42 = qJD(5) * t61 + t106 * t79 + t169 * t78;
t170 = pkin(5) * t42;
t99 = qJD(3) * t107 * pkin(3);
t65 = pkin(4) * t78 + t99;
t17 = -pkin(9) * t41 + t170 + t65;
t3 = -t10 * t105 + t108 * t17 - t142;
t180 = -t3 - t142;
t130 = t141 * t66;
t140 = qJD(6) * t105;
t137 = t61 * t140;
t148 = t108 * t41;
t115 = t137 - t148;
t178 = 2 * m(6);
t177 = 2 * m(7);
t11 = qJD(5) * t24 + t106 * t46 - t110 * t169;
t175 = 0.2e1 * t11;
t174 = 0.2e1 * t65;
t173 = 0.2e1 * t90;
t172 = m(5) * pkin(3);
t166 = t105 * t3;
t13 = -t105 * t24 + t108 * t25;
t145 = qJD(6) * t13;
t2 = t10 * t108 + t105 * t17 + t145;
t165 = t108 * t2;
t164 = t11 * t179;
t76 = t106 * t97 + t167 * t169;
t67 = t76 * qJD(5);
t163 = t179 * t67;
t162 = t42 * t117;
t161 = t117 * t67;
t160 = t61 * t66;
t159 = t78 * t85;
t158 = t79 * t86;
t157 = Ifges(7,5) * t148 + Ifges(7,3) * t42;
t156 = mrSges(6,1) * t42 + mrSges(6,2) * t41;
t154 = Ifges(7,4) * t105;
t153 = Ifges(7,4) * t108;
t152 = Ifges(7,6) * t105;
t150 = t105 * t61;
t147 = t108 * t61;
t73 = pkin(9) + t76;
t144 = qJD(6) * t73;
t138 = 0.2e1 * t109;
t134 = t78 * mrSges(5,1) + mrSges(5,2) * t79;
t132 = -(2 * Ifges(6,4)) - t152;
t128 = -t11 * t117 - t179 * t42;
t126 = mrSges(7,1) * t105 + mrSges(7,2) * t108;
t125 = Ifges(7,1) * t108 - t154;
t124 = -Ifges(7,2) * t105 + t153;
t123 = Ifges(7,5) * t105 + Ifges(7,6) * t108;
t122 = -t105 * t13 + t108 * t14;
t15 = mrSges(7,1) * t42 + mrSges(7,3) * t115;
t16 = -mrSges(7,2) * t42 - mrSges(7,3) * t183;
t121 = -t105 * t15 + t108 * t16;
t43 = mrSges(7,2) * t117 - mrSges(7,3) * t150;
t44 = -mrSges(7,1) * t117 - mrSges(7,3) * t147;
t120 = -t105 * t44 + t108 * t43;
t87 = t126 * qJD(6);
t119 = mrSges(7,3) * t131 - t117 * t87 + t42 * t91 - t156;
t88 = t124 * qJD(6);
t89 = t125 * qJD(6);
t92 = Ifges(7,2) * t108 + t154;
t93 = Ifges(7,1) * t105 + t153;
t114 = t105 * t89 + t108 * t88 + t139 * t93 - t140 * t92;
t112 = -t166 + (-t105 * t14 - t108 * t13) * qJD(6);
t21 = -Ifges(7,6) * t117 + t124 * t61;
t22 = -Ifges(7,5) * t117 + t125 * t61;
t6 = -Ifges(7,4) * t115 - Ifges(7,2) * t183 + Ifges(7,6) * t42;
t7 = -Ifges(7,1) * t115 - Ifges(7,4) * t183 + Ifges(7,5) * t42;
t98 = Ifges(7,5) * t139;
t111 = -t10 * mrSges(6,2) + mrSges(7,3) * t165 - t179 * t87 + t22 * t139 / 0.2e1 + t93 * t148 / 0.2e1 + t105 * t7 / 0.2e1 + Ifges(6,5) * t41 + t108 * t6 / 0.2e1 - t88 * t150 / 0.2e1 + t89 * t147 / 0.2e1 - t117 * (-Ifges(7,6) * t140 + t98) / 0.2e1 + (t123 / 0.2e1 - Ifges(6,6)) * t42 - t183 * t92 / 0.2e1 - (t61 * t93 + t21) * t140 / 0.2e1 + t181 * t11;
t72 = -pkin(5) - t75;
t35 = t126 * t61;
t12 = -mrSges(7,1) * t183 + mrSges(7,2) * t115;
t1 = [-0.2e1 * Ifges(5,2) * t159 + 0.2e1 * Ifges(5,1) * t158 + t134 * t173 + 0.2e1 * t64 * t156 + 0.2e1 * t2 * t43 + 0.2e1 * t3 * t44 + t35 * t175 + 0.2e1 * t14 * t16 + 0.2e1 * t179 * t12 + 0.2e1 * t13 * t15 + t24 * t42 * t176 + (-t105 * t21 + t108 * t22 + t176 * t179) * t41 + (t13 * t3 + t14 * t2 - t164) * t177 + (t10 * t24 + t64 * t65 - t164) * t178 + 0.2e1 * m(5) * (t48 * t58 + t49 * t59) + ((t135 * mrSges(4,2) + Ifges(4,4) * t109) * t138 + (0.2e1 * pkin(3) * (-mrSges(5,1) * t85 + mrSges(5,2) * t86) + 0.2e1 * t135 * mrSges(4,1) + t172 * t173 - 0.2e1 * Ifges(4,4) * t107 + (Ifges(4,1) - Ifges(4,2)) * t138) * t107) * qJD(3) - (mrSges(6,1) * t174 + t10 * t176 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t42 + t132 * t41 + t157) * t117 + (mrSges(6,2) * t174 + mrSges(6,3) * t175 + 0.2e1 * Ifges(6,1) * t41 - t105 * t6 + t108 * t7 + (Ifges(7,5) * t108 + t132) * t42 + (-t105 * t22 - t108 * t21 + t117 * t123) * qJD(6)) * t61 + 0.2e1 * (-t78 * t86 + t79 * t85) * Ifges(5,4) + 0.2e1 * (-t48 * t86 + t49 * t85 - t58 * t79 - t59 * t78) * mrSges(5,3); t117 * t12 + t42 * t35 + t120 * t41 + ((-t105 * t43 - t108 * t44) * qJD(6) + t121) * t61 + m(6) * (t10 * t61 + t24 * t41 + t128) + m(7) * (t122 * t41 + (t112 + t165) * t61 + t128) + m(5) * (t48 * t85 + t49 * t86 - t58 * t78 + t59 * t79); 0.2e1 * m(6) * (t41 * t61 - t162) + 0.2e1 * m(7) * (t131 * t61 - t162) + 0.2e1 * m(5) * (t158 - t159); t111 + m(6) * (t10 * t76 - t11 * t75 + t24 * t66 - t163) + (-mrSges(7,3) * t145 - t44 * t144 + t66 * t43 + t73 * t16 + m(7) * (-t13 * t144 + t14 * t66 + t2 * t73)) * t108 + (-t43 * t144 + t180 * mrSges(7,3) + (-m(7) * t13 - t44) * t66 + (m(7) * t180 - t15) * t73) * t105 + (m(5) * (t102 * t49 + t103 * t48) + (-t102 * t78 - t103 * t79) * mrSges(5,3)) * pkin(3) + (t117 * t66 - t41 * t75 - t42 * t76 + t61 * t67) * mrSges(6,3) - Ifges(5,6) * t78 + Ifges(5,5) * t79 - t72 * t12 + t67 * t35 + t48 * mrSges(5,1) - t49 * mrSges(5,2) + ((-mrSges(4,1) * t96 + Ifges(4,5)) * t109 + (mrSges(4,2) * t96 - Ifges(4,6)) * t107) * qJD(3) + m(7) * (t11 * t72 - t163); (-mrSges(4,1) * t107 - mrSges(4,2) * t109) * qJD(3) + m(6) * (t41 * t76 - t42 * t75 + t160 - t161) + m(7) * (t42 * t72 - t161 + t141 * (t41 * t73 + t160)) + (t102 * t79 - t103 * t78) * t172 + t119 - t134; 0.2e1 * t72 * t87 - 0.2e1 * t182 + (t130 * t73 + t67 * t72) * t177 + (t66 * t76 - t67 * t75) * t178 + t114 + 0.2e1 * t181 * t67 + 0.2e1 * mrSges(7,3) * t130; m(5) * t99 + t105 * t16 + t108 * t15 + t120 * qJD(6) + m(7) * (qJD(6) * t122 + t105 * t2 + t108 * t3) + m(6) * t65 + t134 + t156; 0; 0; 0; t111 + (-m(7) * t11 + t12) * pkin(5) + (-t44 * t139 - t43 * t140 + m(7) * (-t13 * t139 - t14 * t140 + t165 - t166) + t121) * pkin(9) + t112 * mrSges(7,3); m(7) * (pkin(9) * t131 - t170) + t119; -t182 + (t72 - pkin(5)) * t87 + t114 + (m(7) * pkin(9) + mrSges(7,3)) * t130 + (-m(7) * pkin(5) + t181) * t67; 0; -0.2e1 * pkin(5) * t87 + t114; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t137 - Ifges(7,6) * t183 + t157; t12; t98 - t126 * t66 + (-t73 * t155 + (mrSges(7,2) * t73 - Ifges(7,6)) * t105) * qJD(6); -t87; t98 + (pkin(9) * t91 - t152) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
