% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:21
% EndTime: 2019-12-31 17:09:27
% DurationCPUTime: 2.69s
% Computational Cost: add. (1734->317), mult. (4671->480), div. (0->0), fcn. (2993->6), ass. (0->157)
t204 = qJD(2) / 0.2e1;
t124 = sin(qJ(4));
t126 = cos(qJ(4));
t122 = sin(pkin(7));
t123 = cos(pkin(7));
t138 = t122 * t124 - t123 * t126;
t125 = sin(qJ(2));
t127 = cos(qJ(2));
t158 = t123 * t127;
t137 = pkin(3) * t125 - pkin(6) * t158;
t140 = pkin(2) * t125 - qJ(3) * t127;
t101 = t140 * qJD(1);
t157 = qJD(1) * t125;
t150 = t122 * t157;
t61 = pkin(5) * t150 + t123 * t101;
t38 = qJD(1) * t137 + t61;
t159 = t123 * t125;
t160 = t122 * t127;
t135 = -pkin(5) * t159 - pkin(6) * t160;
t86 = t122 * t101;
t48 = qJD(1) * t135 + t86;
t175 = pkin(6) + qJ(3);
t106 = t175 * t122;
t107 = t175 * t123;
t58 = -t106 * t126 - t107 * t124;
t203 = -qJD(3) * t138 + qJD(4) * t58 - t124 * t38 - t126 * t48;
t59 = -t106 * t124 + t107 * t126;
t99 = t122 * t126 + t123 * t124;
t202 = -qJD(3) * t99 - qJD(4) * t59 + t124 * t48 - t126 * t38;
t148 = -Ifges(3,6) * qJD(2) / 0.2e1;
t201 = Ifges(3,5) * t204;
t96 = t123 * qJD(2) - t150;
t149 = t123 * t157;
t97 = t122 * qJD(2) + t149;
t145 = -t124 * t97 + t126 * t96;
t132 = t137 * qJD(2);
t120 = pkin(5) * t157;
t103 = (qJD(3) - t120) * qJD(2);
t82 = qJD(2) * t140 - t125 * qJD(3);
t71 = t82 * qJD(1);
t36 = -t122 * t103 + t123 * t71;
t25 = qJD(1) * t132 + t36;
t155 = qJD(1) * qJD(2);
t146 = t127 * t155;
t143 = t122 * t146;
t37 = t123 * t103 + t122 * t71;
t27 = -pkin(6) * t143 + t37;
t156 = qJD(1) * t127;
t153 = pkin(3) * t156;
t121 = pkin(5) * t156;
t111 = qJD(2) * qJ(3) + t121;
t105 = -pkin(2) * t127 - t125 * qJ(3) - pkin(1);
t90 = t105 * qJD(1);
t50 = -t122 * t111 + t123 * t90;
t26 = -t97 * pkin(6) - t153 + t50;
t51 = t123 * t111 + t122 * t90;
t28 = pkin(6) * t96 + t51;
t8 = -t124 * t28 + t126 * t26;
t1 = qJD(4) * t8 + t124 * t25 + t126 * t27;
t9 = t124 * t26 + t126 * t28;
t2 = -qJD(4) * t9 - t124 * t27 + t126 * t25;
t200 = -t2 * mrSges(5,1) + t1 * mrSges(5,2);
t84 = t138 * qJD(4);
t104 = -qJD(2) * pkin(2) + qJD(3) + t120;
t119 = Ifges(3,4) * t156;
t167 = Ifges(4,2) * t122;
t169 = Ifges(4,4) * t123;
t141 = -t167 + t169;
t170 = Ifges(4,4) * t122;
t142 = Ifges(4,1) * t123 - t170;
t172 = mrSges(4,2) * t123;
t182 = t123 / 0.2e1;
t183 = -t122 / 0.2e1;
t199 = -(t122 * t51 + t123 * t50) * mrSges(4,3) + t104 * (mrSges(4,1) * t122 + t172) + Ifges(3,1) * t157 / 0.2e1 + t119 / 0.2e1 + t201 + (Ifges(4,4) * t97 + Ifges(4,2) * t96 - Ifges(4,6) * t156) * t183 + (Ifges(4,1) * t97 + Ifges(4,4) * t96 - Ifges(4,5) * t156) * t182 + t96 * t141 / 0.2e1 + t97 * t142 / 0.2e1;
t198 = m(4) * pkin(5);
t133 = t138 * t127;
t130 = qJD(2) * t133;
t23 = -qJD(1) * t130 + qJD(4) * t145;
t197 = t23 / 0.2e1;
t134 = t99 * t127;
t131 = qJD(2) * t134;
t47 = t124 * t96 + t126 * t97;
t24 = -qJD(1) * t131 - qJD(4) * t47;
t196 = t24 / 0.2e1;
t195 = -t145 / 0.2e1;
t194 = t145 / 0.2e1;
t193 = -t47 / 0.2e1;
t192 = t47 / 0.2e1;
t77 = t99 * t125;
t191 = -t77 / 0.2e1;
t78 = t138 * t125;
t190 = -t78 / 0.2e1;
t189 = pkin(1) * mrSges(3,1);
t188 = pkin(1) * mrSges(3,2);
t117 = qJD(4) - t156;
t185 = -t117 / 0.2e1;
t184 = t117 / 0.2e1;
t181 = Ifges(5,4) * t47;
t69 = qJD(1) * t134;
t85 = t99 * qJD(4);
t174 = -t69 + t85;
t70 = qJD(1) * t133;
t173 = -t70 + t84;
t171 = Ifges(3,4) * t125;
t168 = Ifges(4,5) * t123;
t166 = Ifges(4,6) * t122;
t152 = pkin(5) * qJD(2) * t125;
t53 = t122 * t152 + t123 * t82;
t66 = pkin(5) * t158 + t122 * t105;
t161 = qJD(2) * mrSges(3,2);
t72 = mrSges(4,1) * t143 + t146 * t172;
t147 = t125 * t155;
t154 = Ifges(5,5) * t23 + Ifges(5,6) * t24 + Ifges(5,3) * t147;
t151 = pkin(3) * t122 + pkin(5);
t7 = -t24 * mrSges(5,1) + t23 * mrSges(5,2);
t144 = m(4) * t104 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t96 + mrSges(4,2) * t97 + mrSges(3,3) * t157;
t95 = t123 * t105;
t49 = -pkin(6) * t159 + t95 + (-pkin(5) * t122 - pkin(3)) * t127;
t55 = -pkin(6) * t122 * t125 + t66;
t18 = -t124 * t55 + t126 * t49;
t19 = t124 * t49 + t126 * t55;
t92 = t151 * t127 * qJD(2);
t136 = t168 / 0.2e1 - t166 / 0.2e1;
t129 = t50 * mrSges(4,1) + t8 * mrSges(5,1) + t117 * Ifges(5,3) + t47 * Ifges(5,5) + t145 * Ifges(5,6) - Ifges(4,3) * t156 / 0.2e1 + Ifges(4,6) * t96 + Ifges(4,5) * t97 + t148 - (Ifges(3,2) * t127 + t171) * qJD(1) / 0.2e1 - t51 * mrSges(4,2) - t9 * mrSges(5,2);
t118 = -pkin(3) * t123 - pkin(2);
t112 = mrSges(3,3) * t156 - t161;
t102 = t151 * t125;
t91 = t122 * t153 + t121;
t81 = qJD(1) * t92;
t80 = (mrSges(4,1) * t125 - mrSges(4,3) * t158) * t155;
t79 = (-mrSges(4,2) * t125 - mrSges(4,3) * t160) * t155;
t73 = t122 * t82;
t68 = -mrSges(4,1) * t156 - t97 * mrSges(4,3);
t67 = mrSges(4,2) * t156 + t96 * mrSges(4,3);
t65 = -pkin(5) * t160 + t95;
t62 = -pkin(5) * t149 + t86;
t60 = -pkin(3) * t96 + t104;
t57 = (Ifges(4,5) * t125 + t127 * t142) * t155;
t56 = (Ifges(4,6) * t125 + t127 * t141) * t155;
t54 = -t123 * t152 + t73;
t43 = Ifges(5,4) * t145;
t39 = qJD(2) * t135 + t73;
t35 = t125 * t84 - t131;
t34 = -t125 * t85 - t130;
t33 = t132 + t53;
t32 = mrSges(5,1) * t117 - mrSges(5,3) * t47;
t31 = -mrSges(5,2) * t117 + mrSges(5,3) * t145;
t17 = -mrSges(5,2) * t147 + mrSges(5,3) * t24;
t16 = mrSges(5,1) * t147 - mrSges(5,3) * t23;
t15 = -mrSges(5,1) * t145 + mrSges(5,2) * t47;
t12 = Ifges(5,1) * t47 + Ifges(5,5) * t117 + t43;
t11 = Ifges(5,2) * t145 + Ifges(5,6) * t117 + t181;
t6 = t23 * Ifges(5,1) + t24 * Ifges(5,4) + Ifges(5,5) * t147;
t5 = t23 * Ifges(5,4) + t24 * Ifges(5,2) + Ifges(5,6) * t147;
t4 = -qJD(4) * t19 - t124 * t39 + t126 * t33;
t3 = qJD(4) * t18 + t124 * t33 + t126 * t39;
t10 = [((-pkin(5) * t112 + t129 + t148) * qJD(2) + t57 * t182 + t56 * t183 + pkin(5) * t72 + (-t122 * t37 - t123 * t36) * mrSges(4,3) + (-0.2e1 * t189 + Ifges(5,5) * t190 + Ifges(5,6) * t191 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t136) * t125) * t155) * t125 + m(4) * (t36 * t65 + t37 * t66 + t50 * t53 + t51 * t54) + (t37 * mrSges(4,2) - t36 * mrSges(4,1) - Ifges(5,6) * t196 - Ifges(5,5) * t197 - t154 / 0.2e1 + (t144 * pkin(5) + t199 + t201) * qJD(2) + (-0.2e1 * t188 + (-0.3e1 / 0.2e1 * t168 + 0.3e1 / 0.2e1 * t166 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t127 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(4,1) * t123 ^ 2 / 0.2e1 - Ifges(5,3) / 0.2e1 + (t172 + t198) * pkin(5) + (pkin(5) * mrSges(4,1) - t169 + t167 / 0.2e1) * t122) * t125) * t155 + t200) * t127 + t102 * t7 + t92 * t15 + t66 * t79 + t65 * t80 + t60 * (-mrSges(5,1) * t35 + mrSges(5,2) * t34) + t54 * t67 + t53 * t68 + m(5) * (t1 * t19 + t102 * t81 + t18 * t2 + t3 * t9 + t4 * t8 + t60 * t92) + t34 * t12 / 0.2e1 + t35 * t11 / 0.2e1 + t3 * t31 + t4 * t32 + t18 * t16 + t19 * t17 + (-t1 * t77 + t2 * t78 - t34 * t8 + t35 * t9) * mrSges(5,3) + (-Ifges(5,4) * t78 - Ifges(5,2) * t77) * t196 + t81 * (mrSges(5,1) * t77 - mrSges(5,2) * t78) + (-Ifges(5,1) * t78 - Ifges(5,4) * t77) * t197 + (Ifges(5,5) * t34 + Ifges(5,6) * t35) * t184 + t6 * t190 + t5 * t191 + (Ifges(5,1) * t34 + Ifges(5,4) * t35) * t192 + (Ifges(5,4) * t34 + Ifges(5,2) * t35) * t194; -t138 * t5 / 0.2e1 + t81 * (mrSges(5,1) * t138 + mrSges(5,2) * t99) + (Ifges(5,4) * t99 - Ifges(5,2) * t138) * t196 + (Ifges(5,1) * t99 - Ifges(5,4) * t138) * t197 + (-t1 * t138 + t173 * t8 - t174 * t9 - t2 * t99) * mrSges(5,3) + (-t85 / 0.2e1 + t69 / 0.2e1) * t11 + (-Ifges(5,5) * t70 - Ifges(5,6) * t69) * t185 + (-Ifges(5,1) * t70 - Ifges(5,4) * t69) * t193 + (-Ifges(5,4) * t70 - Ifges(5,2) * t69) * t195 + (-t84 / 0.2e1 + t70 / 0.2e1) * t12 + (t56 / 0.2e1 + qJD(3) * t67 + qJ(3) * t79 + t37 * mrSges(4,3)) * t123 - m(4) * (t50 * t61 + t51 * t62) + (t57 / 0.2e1 - qJD(3) * t68 - qJ(3) * t80 - t36 * mrSges(4,3)) * t122 + t118 * t7 + t99 * t6 / 0.2e1 - t91 * t15 - pkin(2) * t72 - t62 * t67 - t61 * t68 + t58 * t16 + t59 * t17 + t202 * t32 + t203 * t31 + (t1 * t59 + t118 * t81 + t2 * t58 + t202 * t8 + t203 * t9 - t60 * t91) * m(5) + (((t189 + t171 / 0.2e1) * qJD(1) + (t112 + t161) * pkin(5) - t129 + t148 + (Ifges(4,5) * t122 + Ifges(5,5) * t99 + Ifges(4,6) * t123 - Ifges(5,6) * t138) * t204) * t125 + (-t119 / 0.2e1 + (t188 + t136 * t127 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t125) * qJD(1) + (Ifges(3,5) / 0.2e1 + (Ifges(4,1) * t122 + t169) * t182 + (Ifges(4,2) * t123 + t170) * t183) * qJD(2) + ((-m(4) * pkin(2) - mrSges(4,1) * t123 + mrSges(4,2) * t122 - mrSges(3,1)) * qJD(2) - t144) * pkin(5) - t199) * t127) * qJD(1) + (-Ifges(5,5) * t84 - Ifges(5,6) * t85) * t184 + (-Ifges(5,1) * t84 - Ifges(5,4) * t85) * t192 + (-Ifges(5,4) * t84 - Ifges(5,2) * t85) * t194 + m(4) * ((-t122 * t50 + t123 * t51) * qJD(3) + (-t36 * t122 + t37 * t123) * qJ(3)) + (mrSges(5,1) * t174 - mrSges(5,2) * t173) * t60; t146 * t198 - t145 * t31 + t47 * t32 - t96 * t67 + t97 * t68 - m(4) * (-t50 * t97 + t51 * t96) + t7 + t72 + (-t145 * t9 + t47 * t8 + t81) * m(5); -t60 * (mrSges(5,1) * t47 + mrSges(5,2) * t145) + (Ifges(5,1) * t145 - t181) * t193 + t11 * t192 + (Ifges(5,5) * t145 - Ifges(5,6) * t47) * t185 - t8 * t31 + t9 * t32 + (t145 * t8 + t47 * t9) * mrSges(5,3) + t154 + (-Ifges(5,2) * t47 + t12 + t43) * t195 - t200;];
tauc = t10(:);
