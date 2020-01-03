% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:42
% EndTime: 2019-12-31 19:40:48
% DurationCPUTime: 2.52s
% Computational Cost: add. (1592->329), mult. (3710->425), div. (0->0), fcn. (1705->4), ass. (0->150)
t101 = sin(qJ(5));
t103 = cos(qJ(5));
t121 = mrSges(6,1) * t101 + mrSges(6,2) * t103;
t102 = sin(qJ(2));
t104 = cos(qJ(2));
t127 = pkin(4) * t102 + pkin(7) * t104;
t147 = qJD(1) * t104;
t148 = qJD(1) * t102;
t51 = -qJD(1) * pkin(1) - pkin(2) * t147 - qJ(3) * t148;
t35 = pkin(3) * t147 + qJD(4) - t51;
t21 = qJD(1) * t127 + t35;
t87 = pkin(6) * t148;
t149 = qJD(3) + t87;
t80 = qJ(4) * t148;
t139 = -t80 + t149;
t105 = -pkin(2) - pkin(3);
t97 = -pkin(7) + t105;
t33 = qJD(2) * t97 + t139;
t5 = -t101 * t33 + t103 * t21;
t6 = t101 * t21 + t103 * t33;
t124 = t6 * t101 + t5 * t103;
t55 = qJD(2) * t101 + t103 * t147;
t169 = Ifges(6,4) * t55;
t145 = qJD(2) * t103;
t54 = t101 * t147 - t145;
t78 = qJD(5) + t148;
t16 = Ifges(6,2) * t54 + Ifges(6,6) * t78 - t169;
t53 = Ifges(6,4) * t54;
t17 = -Ifges(6,1) * t55 + Ifges(6,5) * t78 + t53;
t170 = -t103 / 0.2e1;
t171 = t101 / 0.2e1;
t88 = pkin(6) * t147;
t59 = -qJ(4) * t147 + t88;
t98 = qJD(2) * qJ(3);
t44 = -t59 - t98;
t39 = qJD(2) * pkin(4) - t44;
t196 = -t124 * mrSges(6,3) + t39 * t121 - t16 * t171 - t17 * t170;
t195 = -t147 / 0.2e1;
t194 = -mrSges(3,1) - mrSges(4,1);
t151 = Ifges(6,6) * t101;
t152 = Ifges(6,5) * t103;
t116 = -t151 + t152;
t154 = Ifges(6,4) * t103;
t118 = -Ifges(6,2) * t101 + t154;
t155 = Ifges(6,4) * t101;
t120 = Ifges(6,1) * t103 - t155;
t174 = t78 / 0.2e1;
t177 = -t55 / 0.2e1;
t178 = t54 / 0.2e1;
t185 = -qJD(2) / 0.2e1;
t186 = qJD(1) / 0.2e1;
t187 = -qJD(1) / 0.2e1;
t191 = Ifges(4,6) / 0.2e1;
t67 = t88 + t98;
t70 = mrSges(4,2) * t147 + qJD(2) * mrSges(4,3);
t85 = Ifges(4,5) * t148;
t193 = -(m(4) * t67 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t147 + t70) * pkin(6) + t35 * mrSges(5,2) + t51 * mrSges(4,1) + qJD(2) * t191 + Ifges(4,3) * t195 + t85 / 0.2e1 + (Ifges(3,4) * t102 + t104 * Ifges(3,2)) * t187 + (-t104 * Ifges(5,1) - Ifges(5,4) * t102) * t186 - t44 * mrSges(5,3) + t118 * t178 + t120 * t177 - t67 * mrSges(4,2) + t116 * t174 + (Ifges(3,6) + Ifges(5,5)) * t185 + t196;
t192 = -Ifges(3,1) / 0.2e1;
t141 = qJD(2) * qJD(5);
t143 = qJD(5) * t104;
t28 = -t103 * t141 + (t101 * t143 + t102 * t145) * qJD(1);
t190 = -t28 / 0.2e1;
t146 = qJD(2) * t102;
t29 = t101 * t141 + (-t101 * t146 + t103 * t143) * qJD(1);
t189 = -t29 / 0.2e1;
t188 = Ifges(3,4) * t195;
t160 = pkin(6) - qJ(4);
t109 = pkin(4) * t104 + t102 * t97;
t108 = t109 * qJD(2);
t142 = qJD(1) * qJD(2);
t135 = t104 * t142;
t91 = t102 * qJD(3);
t158 = qJ(3) * t135 + qJD(1) * t91;
t12 = qJD(1) * t108 + t158;
t144 = qJD(2) * t104;
t43 = -qJD(4) * t102 + t144 * t160;
t37 = t43 * qJD(1);
t1 = qJD(5) * t5 + t101 * t12 + t103 * t37;
t2 = -qJD(5) * t6 - t101 * t37 + t103 * t12;
t126 = -t1 * t103 + t2 * t101;
t123 = -t101 * t5 + t103 * t6;
t56 = t87 - t80;
t184 = -t56 - qJD(3);
t183 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t28 + Ifges(6,6) * t29;
t137 = t105 * qJD(2);
t38 = t137 + t139;
t62 = -qJD(2) * pkin(2) + t149;
t181 = (m(4) * t62 + (mrSges(4,2) + mrSges(3,3)) * t148 + t194 * qJD(2)) * pkin(6) - t38 * mrSges(5,3) - t51 * mrSges(4,3) - t6 * mrSges(6,2) + t78 * Ifges(6,3) - t55 * Ifges(6,5) + t54 * Ifges(6,6) - (-Ifges(5,4) * t104 - t102 * Ifges(5,2)) * t186 - (t102 * Ifges(4,1) - Ifges(4,5) * t104) * t187 - t148 * t192 - t188 + t35 * mrSges(5,1) + t5 * mrSges(6,1) + t62 * mrSges(4,2) - (Ifges(5,6) + Ifges(4,4) + Ifges(3,5)) * t185;
t179 = -t54 / 0.2e1;
t176 = t55 / 0.2e1;
t175 = -t78 / 0.2e1;
t173 = pkin(1) * mrSges(3,1);
t172 = pkin(1) * mrSges(3,2);
t110 = pkin(6) * t146 + qJD(4) * t104;
t136 = t102 * t142;
t96 = qJD(2) * qJD(3);
t27 = -qJ(4) * t136 + qJD(1) * t110 - t96;
t72 = t160 * t104;
t164 = t27 * t72;
t159 = -qJD(2) * mrSges(5,1) + mrSges(6,1) * t54 + mrSges(6,2) * t55 + mrSges(5,3) * t147;
t157 = mrSges(5,1) * t135 + mrSges(5,2) * t136;
t156 = qJ(3) * t144 + t91;
t153 = Ifges(6,5) * t101;
t150 = Ifges(6,6) * t103;
t63 = -t104 * pkin(2) - t102 * qJ(3) - pkin(1);
t140 = t70 - t159;
t138 = m(4) * pkin(6) + mrSges(4,2);
t52 = t104 * pkin(3) - t63;
t132 = t102 * t137;
t131 = -0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t130 = -Ifges(5,5) / 0.2e1 + t191 - Ifges(3,6) / 0.2e1;
t129 = Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t125 = t1 * t101 + t103 * t2;
t122 = -mrSges(6,1) * t103 + mrSges(6,2) * t101;
t119 = Ifges(6,1) * t101 + t154;
t117 = Ifges(6,2) * t103 + t155;
t18 = mrSges(6,1) * t135 - mrSges(6,3) * t28;
t19 = -mrSges(6,2) * t135 + mrSges(6,3) * t29;
t115 = -t101 * t18 + t103 * t19;
t30 = -mrSges(6,2) * t78 + mrSges(6,3) * t54;
t31 = mrSges(6,1) * t78 + mrSges(6,3) * t55;
t114 = -t101 * t31 + t103 * t30;
t113 = -t101 * t30 - t103 * t31;
t32 = t127 + t52;
t71 = t160 * t102;
t13 = -t101 * t71 + t103 * t32;
t14 = t101 * t32 + t103 * t71;
t100 = qJ(3) + pkin(4);
t82 = qJ(3) * t147;
t77 = Ifges(6,3) * t135;
t64 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t148;
t61 = -pkin(6) * t136 + t96;
t60 = (mrSges(5,1) * t102 - mrSges(5,2) * t104) * qJD(1);
t57 = (-mrSges(4,1) * t104 - mrSges(4,3) * t102) * qJD(1);
t42 = pkin(2) * t146 - t156;
t41 = -qJ(4) * t146 + t110;
t40 = t105 * t148 + t82;
t36 = pkin(2) * t136 - t158;
t34 = t132 + t156;
t24 = qJD(1) * t132 + t158;
t23 = qJD(1) * t109 + t82;
t20 = t108 + t156;
t11 = t101 * t23 + t103 * t59;
t10 = -t101 * t59 + t103 * t23;
t9 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t8 = t28 * Ifges(6,1) + t29 * Ifges(6,4) + Ifges(6,5) * t135;
t7 = t28 * Ifges(6,4) + t29 * Ifges(6,2) + Ifges(6,6) * t135;
t4 = -qJD(5) * t14 - t101 * t43 + t103 * t20;
t3 = qJD(5) * t13 + t101 * t20 + t103 * t43;
t15 = [t52 * t157 + t43 * t64 + t72 * t9 + t42 * t57 + t34 * t60 + t3 * t30 + t4 * t31 + t14 * t19 + t13 * t18 + t159 * t41 + m(6) * (t1 * t14 + t13 * t2 + t3 * t6 - t39 * t41 + t4 * t5 - t164) + m(5) * (t24 * t52 + t34 * t35 + t37 * t71 + t38 * t43 + t41 * t44 - t164) + m(4) * (t36 * t63 + t42 * t51) + (t77 / 0.2e1 - t36 * mrSges(4,3) + t24 * mrSges(5,1) - t37 * mrSges(5,3) + (t130 * qJD(2) + (t63 * mrSges(4,1) + t72 * mrSges(5,3) + t102 * t131 - 0.2e1 * t173) * qJD(1) + t193) * qJD(2) + t183) * t102 + (t120 * t190 + t118 * t189 - t36 * mrSges(4,1) - t24 * mrSges(5,2) + t8 * t170 + t7 * t171 + t138 * t61 + (mrSges(5,3) + t121) * t27 + t125 * mrSges(6,3) + (t103 * t16 / 0.2e1 + t17 * t171 + (t150 + t153) * t174 + t117 * t178 + t119 * t177 + t39 * t122 + t123 * mrSges(6,3)) * qJD(5) + (t129 * qJD(2) + ((-t152 / 0.2e1 + t151 / 0.2e1 - t131) * t104 - t63 * mrSges(4,3) - 0.2e1 * t172 - t71 * mrSges(5,3) + (0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + t138 * pkin(6)) * t102) * qJD(1) + t181) * qJD(2)) * t104; t100 * t9 - t101 * t8 / 0.2e1 + t61 * mrSges(4,3) - t59 * t64 - t40 * t60 + t37 * mrSges(5,2) - t11 * t30 - t10 * t31 + t126 * mrSges(6,3) + t117 * t189 + t119 * t190 + (-mrSges(5,1) + t122) * t27 + t140 * qJD(3) - t159 * t56 + t7 * t170 + ((t188 + (t172 + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t104) * qJD(1) - t181) * t104 + ((t192 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t147 + (t173 + (Ifges(3,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t102) * qJD(1) - t85 / 0.2e1 - t193) * t102 + ((pkin(6) * mrSges(3,2) + (-mrSges(4,2) + mrSges(5,3)) * qJ(3) + t130) * t102 + (-pkin(2) * mrSges(4,2) - t105 * mrSges(5,3) - t153 / 0.2e1 - t150 / 0.2e1 + (-m(4) * pkin(2) + t194) * pkin(6) + t129) * t104) * qJD(2)) * qJD(1) + (t116 * t175 + t118 * t179 + t120 * t176 - t196) * qJD(5) + m(4) * (qJ(3) * t61 + qJD(3) * t67) + (-t10 * t5 - t100 * t27 - t11 * t6 - t184 * t39) * m(6) + (t115 - m(6) * t126 + (-m(6) * t124 + t113) * qJD(5)) * t97 + (-m(4) * t51 - t57) * (pkin(2) * t148 - t82) + (-t27 * qJ(3) + t105 * t37 + t184 * t44 - t35 * t40 - t38 * t59) * m(5); t113 * qJD(5) - t140 * qJD(2) + ((-mrSges(5,3) + t138) * t144 + (t113 + t57 - t60) * t102) * qJD(1) - m(4) * (qJD(2) * t67 - t148 * t51) + t115 + (-t39 * qJD(2) - t124 * t78 - t126) * m(6) + (qJD(2) * t44 - t148 * t35 + t37) * m(5); t101 * t19 + t103 * t18 + t114 * qJD(5) + m(5) * t24 + m(6) * (qJD(5) * t123 + t125) + ((-m(5) * t44 + m(6) * t39 - t159) * t104 + (m(5) * t38 + m(6) * t123 + t114 + t64) * t102) * qJD(1) + t157; t77 - t39 * (-mrSges(6,1) * t55 + mrSges(6,2) * t54) + (Ifges(6,1) * t54 + t169) * t176 + t16 * t177 + (Ifges(6,5) * t54 + Ifges(6,6) * t55) * t175 - t5 * t30 + t6 * t31 + (t5 * t54 - t55 * t6) * mrSges(6,3) + (Ifges(6,2) * t55 + t17 + t53) * t179 + t183;];
tauc = t15(:);
