% Calculate time derivative of joint inertia matrix for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:08
% EndTime: 2019-12-31 21:16:13
% DurationCPUTime: 1.85s
% Computational Cost: add. (3507->290), mult. (8000->439), div. (0->0), fcn. (7270->8), ass. (0->127)
t128 = sin(pkin(9));
t129 = cos(pkin(9));
t182 = -mrSges(5,1) * t129 + mrSges(5,2) * t128 - mrSges(4,1);
t132 = sin(qJ(2));
t171 = -pkin(7) - pkin(6);
t117 = t171 * t132;
t135 = cos(qJ(2));
t118 = t171 * t135;
t131 = sin(qJ(3));
t134 = cos(qJ(3));
t181 = t134 * t117 + t118 * t131;
t109 = t131 * t132 - t134 * t135;
t110 = t131 * t135 + t134 * t132;
t124 = -pkin(2) * t135 - pkin(1);
t83 = t109 * pkin(3) - t110 * qJ(4) + t124;
t98 = t117 * t131 - t118 * t134;
t49 = -t128 * t98 + t129 * t83;
t50 = t128 * t83 + t129 * t98;
t180 = -t128 * t49 + t129 * t50;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t137 = t128 * t130 - t129 * t133;
t102 = t137 * qJD(5);
t179 = qJD(2) + qJD(3);
t163 = pkin(2) * qJD(3);
t108 = t128 * t133 + t129 * t130;
t89 = mrSges(6,1) * t137 + mrSges(6,2) * t108;
t178 = (-mrSges(4,2) * t134 + (t89 + t182) * t131) * t163;
t177 = 2 * m(5);
t176 = 2 * m(6);
t145 = qJD(2) * t171;
t112 = t132 * t145;
t141 = t135 * t145;
t57 = t98 * qJD(3) + t112 * t131 - t134 * t141;
t175 = 0.2e1 * t57;
t103 = t108 * qJD(5);
t77 = mrSges(6,1) * t103 - t102 * mrSges(6,2);
t174 = 0.2e1 * t77;
t127 = t129 ^ 2;
t173 = 0.2e1 * t124;
t170 = pkin(2) * t134;
t92 = t179 * t109;
t93 = t179 * t110;
t40 = pkin(2) * qJD(2) * t132 + pkin(3) * t93 + qJ(4) * t92 - qJD(4) * t110;
t56 = t181 * qJD(3) + t134 * t112 + t131 * t141;
t15 = -t128 * t56 + t129 * t40;
t169 = t15 * mrSges(5,3);
t168 = t57 * t181;
t16 = t128 * t40 + t129 * t56;
t157 = t129 * t92;
t160 = t128 * t92;
t45 = -mrSges(5,1) * t160 - mrSges(5,2) * t157;
t166 = Ifges(5,4) * t128;
t165 = Ifges(5,4) * t129;
t164 = Ifges(5,2) * t128;
t162 = t128 * t15;
t159 = t129 * t16;
t156 = t131 * t181;
t155 = t110 * t128;
t154 = t110 * t129;
t120 = pkin(2) * t131 + qJ(4);
t152 = t120 * t129;
t151 = -Ifges(6,5) * t102 - Ifges(6,6) * t103;
t150 = t128 ^ 2 + t127;
t149 = 2 * mrSges(6,3);
t148 = 0.2e1 * t135;
t25 = -t110 * t103 + t137 * t92;
t26 = t110 * t102 + t108 * t92;
t147 = Ifges(6,5) * t25 + Ifges(6,6) * t26 + Ifges(6,3) * t93;
t146 = t131 * t163;
t121 = -pkin(4) * t129 - pkin(3);
t8 = -t26 * mrSges(6,1) + t25 * mrSges(6,2);
t78 = -Ifges(6,4) * t102 - Ifges(6,2) * t103;
t79 = -Ifges(6,1) * t102 - Ifges(6,4) * t103;
t90 = Ifges(6,4) * t108 - Ifges(6,2) * t137;
t91 = Ifges(6,1) * t108 - Ifges(6,4) * t137;
t144 = -t102 * t91 - t103 * t90 + t108 * t79 - t137 * t78;
t119 = t134 * t163 + qJD(4);
t143 = t150 * t119;
t142 = t150 * qJD(4);
t139 = Ifges(5,5) * t129 - Ifges(5,6) * t128;
t33 = pkin(4) * t109 - pkin(8) * t154 + t49;
t39 = -pkin(8) * t155 + t50;
t11 = -t130 * t39 + t133 * t33;
t12 = t130 * t33 + t133 * t39;
t138 = 0.2e1 * t150 * mrSges(5,3);
t105 = (-pkin(8) - t120) * t128;
t125 = t129 * pkin(8);
t106 = t125 + t152;
t75 = t105 * t133 - t106 * t130;
t76 = t105 * t130 + t106 * t133;
t114 = (-pkin(8) - qJ(4)) * t128;
t116 = qJ(4) * t129 + t125;
t95 = t114 * t133 - t116 * t130;
t96 = t114 * t130 + t116 * t133;
t13 = pkin(8) * t160 + t16;
t9 = pkin(4) * t93 + pkin(8) * t157 + t15;
t2 = t11 * qJD(5) + t13 * t133 + t130 * t9;
t3 = -t12 * qJD(5) - t13 * t130 + t133 * t9;
t31 = t93 * Ifges(5,6) - (-t164 + t165) * t92;
t32 = t93 * Ifges(5,5) - (Ifges(5,1) * t129 - t166) * t92;
t71 = t108 * t110;
t72 = t137 * t110;
t34 = -Ifges(6,4) * t72 - Ifges(6,2) * t71 + Ifges(6,6) * t109;
t35 = -Ifges(6,1) * t72 - Ifges(6,4) * t71 + Ifges(6,5) * t109;
t36 = -pkin(4) * t160 + t57;
t6 = Ifges(6,4) * t25 + Ifges(6,2) * t26 + t93 * Ifges(6,6);
t69 = pkin(4) * t155 - t181;
t7 = Ifges(6,1) * t25 + Ifges(6,4) * t26 + t93 * Ifges(6,5);
t136 = mrSges(5,3) * t159 - t56 * mrSges(4,2) + t69 * t77 - t71 * t78 / 0.2e1 - t72 * t79 / 0.2e1 + t36 * t89 + t26 * t90 / 0.2e1 + t25 * t91 / 0.2e1 - Ifges(4,5) * t92 - Ifges(4,6) * t93 - t102 * t35 / 0.2e1 - t103 * t34 / 0.2e1 - t137 * t6 / 0.2e1 + t108 * t7 / 0.2e1 + t128 * t32 / 0.2e1 + t129 * t31 / 0.2e1 + t109 * t151 / 0.2e1 - (Ifges(5,1) * t128 + t165) * t157 / 0.2e1 + (Ifges(5,2) * t129 + t166) * t160 / 0.2e1 + t182 * t57 + (Ifges(5,5) * t128 + Ifges(6,5) * t108 + Ifges(5,6) * t129 - Ifges(6,6) * t137) * t93 / 0.2e1 + (t11 * t102 - t12 * t103 - t3 * t108 - t137 * t2) * mrSges(6,3);
t123 = -pkin(3) - t170;
t113 = t121 - t170;
t85 = mrSges(5,1) * t109 - mrSges(5,3) * t154;
t84 = -mrSges(5,2) * t109 - mrSges(5,3) * t155;
t80 = (mrSges(5,1) * t128 + mrSges(5,2) * t129) * t110;
t66 = -t108 * qJD(4) - t96 * qJD(5);
t65 = -t137 * qJD(4) + t95 * qJD(5);
t59 = mrSges(6,1) * t109 + mrSges(6,3) * t72;
t58 = -mrSges(6,2) * t109 - mrSges(6,3) * t71;
t52 = mrSges(5,1) * t93 + mrSges(5,3) * t157;
t51 = -mrSges(5,2) * t93 + mrSges(5,3) * t160;
t48 = -t76 * qJD(5) - t108 * t119;
t47 = t75 * qJD(5) - t137 * t119;
t41 = mrSges(6,1) * t71 - mrSges(6,2) * t72;
t18 = -mrSges(6,2) * t93 + mrSges(6,3) * t26;
t17 = mrSges(6,1) * t93 - mrSges(6,3) * t25;
t1 = [t93 * (-Ifges(6,5) * t72 - Ifges(6,6) * t71) + t80 * t175 + (t11 * t3 + t12 * t2 + t36 * t69) * t176 + (t93 * mrSges(4,1) - t92 * mrSges(4,2)) * t173 + (mrSges(4,3) * t175 - t128 * t31 + t129 * t32 + (-(2 * Ifges(4,4)) + t139) * t93 - (Ifges(5,1) * t127 + (2 * Ifges(4,1)) + (t164 - 0.2e1 * t165) * t128) * t92) * t110 + (t15 * t49 + t16 * t50 - t168) * t177 + 0.2e1 * m(4) * (t56 * t98 - t168) + (-0.2e1 * mrSges(4,3) * t56 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3)) * t93 + t147) * t109 + 0.2e1 * (t181 * t92 - t98 * t93) * mrSges(4,3) - 0.2e1 * t181 * t45 - 0.2e1 * (-Ifges(4,4) + t139) * t92 * t109 + 0.2e1 * t11 * t17 + 0.2e1 * t12 * t18 + t26 * t34 + t25 * t35 + 0.2e1 * t36 * t41 + 0.2e1 * t50 * t51 + 0.2e1 * t49 * t52 + 0.2e1 * t2 * t58 + 0.2e1 * t3 * t59 + 0.2e1 * t69 * t8 - t71 * t6 - t72 * t7 + 0.2e1 * t16 * t84 + 0.2e1 * t15 * t85 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t135) * t148 + (0.2e1 * pkin(2) * (mrSges(4,1) * t109 + mrSges(4,2) * t110) + m(4) * pkin(2) * t173 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t132 + (Ifges(3,1) - Ifges(3,2)) * t148) * t132) * qJD(2); (m(4) * (t131 * t56 - t134 * t57) + (-t131 * t93 + t134 * t92) * mrSges(4,3) + (-t134 * t109 * mrSges(4,3) - m(5) * t156 + m(4) * (t134 * t98 - t156) + (m(6) * t69 + t110 * mrSges(4,3) + t41 + t80) * t131) * qJD(3)) * pkin(2) + (Ifges(3,5) * t135 - Ifges(3,6) * t132 + (-mrSges(3,1) * t135 + mrSges(3,2) * t132) * pkin(6)) * qJD(2) + (t119 * t84 + t120 * t51) * t129 + (-t119 * t85 - t120 * t52 - t169) * t128 + m(5) * (t180 * t119 - t120 * t162 + t123 * t57 + t16 * t152) + m(6) * (t11 * t48 + t113 * t36 + t12 * t47 + t2 * t76 + t3 * t75) + t136 + t47 * t58 + t48 * t59 + t75 * t17 + t76 * t18 + t113 * t8 + t123 * t45; t113 * t174 + t119 * t138 + 0.2e1 * t178 + (t113 * t146 + t47 * t76 + t48 * t75) * t176 + (t120 * t143 + t123 * t146) * t177 + (t75 * t102 - t76 * t103 - t48 * t108 - t137 * t47) * t149 + t144; m(5) * (-pkin(3) * t57 + t180 * qJD(4) + (t159 - t162) * qJ(4)) + (qJ(4) * t51 + qJD(4) * t84) * t129 + (-qJ(4) * t52 - qJD(4) * t85 - t169) * t128 + m(6) * (t11 * t66 + t12 * t65 + t121 * t36 + t2 * t96 + t3 * t95) + t136 - pkin(3) * t45 + t65 * t58 + t66 * t59 + t95 * t17 + t96 * t18 + t121 * t8; (t113 + t121) * t77 + t178 + m(6) * (t121 * t146 + t47 * t96 + t48 * t95 + t65 * t76 + t66 * t75) + m(5) * (-pkin(3) * t146 + qJ(4) * t143 + t120 * t142) + (t143 + t142) * mrSges(5,3) + ((-t48 - t66) * t108 - (t47 + t65) * t137 - (t76 + t96) * t103 - (-t75 - t95) * t102) * mrSges(6,3) + t144; t121 * t174 + (t65 * t96 + t66 * t95) * t176 + (t150 * qJ(4) * t177 + t138) * qJD(4) + (t95 * t102 - t96 * t103 - t66 * t108 - t137 * t65) * t149 + t144; m(5) * t57 + m(6) * t36 + t45 + t8; (m(5) + m(6)) * t146 + t77; t77; 0; mrSges(6,1) * t3 - mrSges(6,2) * t2 + t147; mrSges(6,1) * t48 - mrSges(6,2) * t47 + t151; mrSges(6,1) * t66 - mrSges(6,2) * t65 + t151; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
