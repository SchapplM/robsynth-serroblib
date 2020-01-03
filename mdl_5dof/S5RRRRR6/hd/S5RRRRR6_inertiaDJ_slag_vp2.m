% Calculate time derivative of joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:14:48
% EndTime: 2020-01-03 12:14:52
% DurationCPUTime: 1.52s
% Computational Cost: add. (3591->248), mult. (7806->350), div. (0->0), fcn. (6753->8), ass. (0->130)
t115 = sin(qJ(3));
t119 = cos(qJ(3));
t184 = t115 ^ 2 + t119 ^ 2;
t113 = sin(qJ(5));
t117 = cos(qJ(5));
t114 = sin(qJ(4));
t118 = cos(qJ(4));
t140 = qJD(4) * t118;
t141 = qJD(4) * t114;
t116 = sin(qJ(2));
t103 = pkin(1) * t116 + pkin(7);
t158 = -pkin(8) - t103;
t128 = qJD(3) * t158;
t120 = cos(qJ(2));
t153 = pkin(1) * qJD(2);
t131 = t120 * t153;
t66 = t115 * t128 + t119 * t131;
t67 = -t115 * t131 + t119 * t128;
t86 = t158 * t115;
t110 = t119 * pkin(8);
t87 = t103 * t119 + t110;
t23 = t114 * t67 + t118 * t66 + t86 * t140 - t141 * t87;
t180 = qJD(3) + qJD(4);
t90 = t114 * t119 + t115 * t118;
t59 = t180 * t90;
t57 = t59 * pkin(9);
t12 = t23 - t57;
t89 = -t114 * t115 + t118 * t119;
t58 = t180 * t89;
t169 = pkin(9) * t58;
t51 = t114 * t86 + t118 * t87;
t24 = -t51 * qJD(4) - t114 * t66 + t118 * t67;
t13 = t24 - t169;
t168 = pkin(9) * t90;
t50 = -t114 * t87 + t118 * t86;
t39 = t50 - t168;
t84 = t89 * pkin(9);
t40 = t84 + t51;
t16 = -t113 * t40 + t117 * t39;
t2 = qJD(5) * t16 + t113 * t13 + t117 * t12;
t17 = t113 * t39 + t117 * t40;
t3 = -qJD(5) * t17 - t113 * t12 + t117 * t13;
t183 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t171 = -pkin(8) - pkin(7);
t101 = t171 * t115;
t102 = pkin(7) * t119 + t110;
t130 = qJD(3) * t171;
t95 = t115 * t130;
t96 = t119 * t130;
t37 = t101 * t140 - t102 * t141 + t114 * t96 + t118 * t95;
t25 = t37 - t57;
t65 = t114 * t101 + t118 * t102;
t38 = -t65 * qJD(4) - t114 * t95 + t118 * t96;
t26 = t38 - t169;
t64 = t118 * t101 - t102 * t114;
t42 = t64 - t168;
t43 = t84 + t65;
t27 = -t113 * t43 + t117 * t42;
t5 = qJD(5) * t27 + t113 * t26 + t117 * t25;
t28 = t113 * t42 + t117 * t43;
t6 = -qJD(5) * t28 - t113 * t25 + t117 * t26;
t182 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t181 = t184 * t120;
t179 = t24 * mrSges(5,1) - t23 * mrSges(5,2) + t183;
t178 = t38 * mrSges(5,1) - t37 * mrSges(5,2) + t182;
t99 = -mrSges(4,1) * t119 + mrSges(4,2) * t115;
t177 = (t184 * mrSges(4,3) - mrSges(3,2)) * t120 + (t99 - mrSges(3,1)) * t116;
t176 = 2 * m(5);
t175 = 2 * m(6);
t53 = -t113 * t90 + t117 * t89;
t20 = qJD(5) * t53 - t113 * t59 + t117 * t58;
t54 = t113 * t89 + t117 * t90;
t21 = -qJD(5) * t54 - t113 * t58 - t117 * t59;
t9 = -mrSges(6,1) * t21 + mrSges(6,2) * t20;
t174 = 0.2e1 * t9;
t32 = -mrSges(6,1) * t53 + mrSges(6,2) * t54;
t173 = 0.2e1 * t32;
t33 = mrSges(5,1) * t59 + mrSges(5,2) * t58;
t172 = 0.2e1 * t33;
t60 = -mrSges(5,1) * t89 + mrSges(5,2) * t90;
t170 = pkin(3) * t60;
t165 = t89 * pkin(4);
t164 = pkin(1) * t120;
t163 = t21 * mrSges(6,3);
t104 = pkin(3) * t118 + pkin(4);
t138 = qJD(5) * t117;
t139 = qJD(5) * t113;
t145 = t113 * t114;
t48 = t104 * t138 + (-t114 * t139 + (t117 * t118 - t145) * qJD(4)) * pkin(3);
t160 = t48 * mrSges(6,2);
t143 = qJD(3) * t115;
t108 = pkin(3) * t143;
t109 = t116 * t153;
t97 = t109 + t108;
t159 = t97 * t60;
t157 = Ifges(6,5) * t20 + Ifges(6,6) * t21;
t156 = Ifges(4,4) * t115;
t154 = Ifges(4,6) * t115;
t152 = pkin(3) * qJD(4);
t151 = pkin(4) * qJD(5);
t144 = t114 * t117;
t142 = qJD(3) * t119;
t137 = 2 * mrSges(5,3);
t136 = 0.2e1 * mrSges(6,3);
t135 = mrSges(5,3) * t152;
t134 = mrSges(6,3) * t151;
t133 = t117 * t20 * mrSges(6,3);
t132 = t118 * t58 * mrSges(5,3);
t106 = -t119 * pkin(3) - pkin(2);
t49 = -t104 * t139 + (-t114 * t138 + (-t113 * t118 - t144) * qJD(4)) * pkin(3);
t47 = t49 * mrSges(6,1);
t129 = t47 - t160;
t46 = pkin(4) * t59 + t108;
t127 = Ifges(5,5) * t58 - Ifges(5,6) * t59 + t157;
t125 = mrSges(4,1) * t115 + mrSges(4,2) * t119;
t98 = t106 - t164;
t124 = t117 * t53 * t134 + t127 + (pkin(4) * t163 + t134 * t54) * t113;
t123 = (-mrSges(5,1) * t114 - mrSges(5,2) * t118) * t152;
t122 = 0.2e1 * t20 * t54 * Ifges(6,1) + 0.2e1 * t21 * Ifges(6,2) * t53 - 0.2e1 * t89 * Ifges(5,2) * t59 + 0.2e1 * t90 * t58 * Ifges(5,1) + (Ifges(4,1) * t119 - t156) * t143 + (0.2e1 * Ifges(4,4) * t119 + (Ifges(4,1) - Ifges(4,2)) * t115) * t142 + 0.2e1 * (t20 * t53 + t21 * t54) * Ifges(6,4) + 0.2e1 * (t89 * t58 - t90 * t59) * Ifges(5,4);
t76 = -pkin(3) * t145 + t104 * t117;
t77 = pkin(3) * t144 + t104 * t113;
t121 = t118 * t89 * t135 + Ifges(4,5) * t142 + t77 * t163 + t127 + (-pkin(3) * t59 * mrSges(5,3) + t90 * t135) * t114 + (-t20 * t76 + t48 * t53 - t49 * t54) * mrSges(6,3);
t105 = -pkin(2) - t164;
t100 = Ifges(4,2) * t119 + t156;
t94 = t125 * qJD(3);
t85 = (-t113 * mrSges(6,1) - t117 * mrSges(6,2)) * t151;
t69 = t106 - t165;
t68 = t98 - t165;
t41 = t109 + t46;
t1 = [t122 + (t16 * t3 + t17 * t2 + t41 * t68) * t175 + (t23 * t51 + t24 * t50 + t97 * t98) * t176 + 0.2e1 * (m(4) * (t181 * t103 + t105 * t116) + t177) * t153 + (-t16 * t20 + t17 * t21 + t2 * t53 - t3 * t54) * t136 + (t23 * t89 - t24 * t90 - t50 * t58 - t51 * t59) * t137 - t100 * t143 + t41 * t173 + t68 * t174 + 0.2e1 * t159 + t98 * t172 + 0.2e1 * t105 * t94; t122 + m(6) * (t16 * t6 + t17 * t5 + t2 * t28 + t27 * t3 + t41 * t69 + t46 * t68) + m(5) * (t106 * t97 + t108 * t98 + t23 * t65 + t24 * t64 + t37 * t51 + t38 * t50) + (-pkin(2) + t105) * t94 + (t68 + t69) * t9 + (t98 + t106) * t33 + (t41 + t46) * t32 + (-t100 + t170) * t143 + ((-t24 - t38) * t90 + (t23 + t37) * t89 - (t51 + t65) * t59 + (-t50 - t64) * t58) * mrSges(5,3) + ((-t3 - t6) * t54 + (t2 + t5) * t53 + (t17 + t28) * t21 + (-t16 - t27) * t20) * mrSges(6,3) + (m(4) * (-pkin(2) * t116 + t181 * pkin(7)) + t177) * t153 + t159; t122 + (t27 * t6 + t28 * t5 + t46 * t69) * t175 + (t106 * t108 + t37 * t65 + t38 * t64) * t176 + (-t100 + 0.2e1 * t170) * t143 + (-t20 * t27 + t21 * t28 + t5 * t53 - t54 * t6) * t136 + (t37 * t89 - t38 * t90 - t58 * t64 - t59 * t65) * t137 + t46 * t173 + t69 * t174 - 0.2e1 * pkin(2) * t94 + t106 * t172; t121 + (-t132 + m(5) * (t114 * t23 + t118 * t24 + t140 * t51 - t141 * t50)) * pkin(3) - t125 * t131 + (t103 * t99 - t154) * qJD(3) + m(6) * (t16 * t49 + t17 * t48 + t2 * t77 + t3 * t76) + t179; t121 + (-t132 + m(5) * (t114 * t37 + t118 * t38 + t140 * t65 - t141 * t64)) * pkin(3) + (pkin(7) * t99 - t154) * qJD(3) + m(6) * (t27 * t49 + t28 * t48 + t5 * t77 + t6 * t76) + t178; -0.2e1 * t160 + 0.2e1 * t47 + (t48 * t77 + t49 * t76) * t175 + 0.2e1 * t123; (-t133 + m(6) * (t113 * t2 + t117 * t3 + t138 * t17 - t139 * t16)) * pkin(4) + t124 + t179; (-t133 + m(6) * (t113 * t5 + t117 * t6 + t138 * t28 - t139 * t27)) * pkin(4) + t124 + t178; t123 + (m(6) * (t113 * t48 + t117 * t49 + t138 * t77 - t76 * t139) - mrSges(6,2) * t138 - mrSges(6,1) * t139) * pkin(4) + t129; 0.2e1 * t85; t157 + t183; t157 + t182; t129; t85; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
