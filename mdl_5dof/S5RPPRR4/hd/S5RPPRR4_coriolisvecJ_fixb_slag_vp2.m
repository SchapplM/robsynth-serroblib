% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:17
% EndTime: 2022-01-23 09:16:24
% DurationCPUTime: 2.82s
% Computational Cost: add. (3526->291), mult. (9758->433), div. (0->0), fcn. (7228->8), ass. (0->141)
t128 = sin(pkin(9));
t130 = cos(pkin(9));
t133 = sin(qJ(4));
t135 = cos(qJ(4));
t112 = t128 * t135 + t130 * t133;
t106 = t112 * qJD(4);
t131 = cos(pkin(8));
t138 = qJD(1) * t112;
t177 = t131 * t138 - t106;
t157 = t130 * t135;
t141 = t128 * t133 - t157;
t105 = t141 * qJD(4);
t153 = qJD(1) * t131;
t176 = t141 * t153 - t105;
t148 = qJD(4) * t135;
t149 = qJD(4) * t133;
t129 = sin(pkin(8));
t158 = t129 * t130;
t137 = -pkin(6) * t158 + (-qJ(2) * t128 - pkin(3)) * t131;
t114 = -pkin(2) * t131 - qJ(3) * t129 - pkin(1);
t101 = qJD(1) * t114 + qJD(2);
t92 = t130 * t101;
t55 = qJD(1) * t137 + t92;
t154 = qJD(1) * t129;
t145 = t128 * t154;
t146 = qJ(2) * t153;
t72 = t128 * t101 + t130 * t146;
t60 = -pkin(6) * t145 + t72;
t150 = qJD(3) * t129;
t151 = qJD(2) * t131;
t103 = -t128 * t151 - t130 * t150;
t96 = t103 * qJD(1);
t104 = -t128 * t150 + t130 * t151;
t97 = t104 * qJD(1);
t14 = t133 * t96 + t135 * t97 + t55 * t148 - t149 * t60;
t90 = t129 * t105;
t80 = qJD(1) * t90;
t12 = pkin(7) * t80 + t14;
t36 = t133 * t55 + t135 * t60;
t15 = -qJD(4) * t36 - t133 * t97 + t135 * t96;
t89 = t129 * t106;
t79 = qJD(1) * t89;
t13 = pkin(7) * t79 + t15;
t132 = sin(qJ(5));
t134 = cos(qJ(5));
t85 = t129 * t138;
t22 = -pkin(7) * t85 + t36;
t162 = t134 * t22;
t123 = qJD(4) - t153;
t35 = -t133 * t60 + t135 * t55;
t87 = -t133 * t145 + t154 * t157;
t21 = -pkin(7) * t87 + t35;
t18 = pkin(4) * t123 + t21;
t7 = t132 * t18 + t162;
t3 = -qJD(5) * t7 - t12 * t132 + t13 * t134;
t1 = t3 * mrSges(6,1);
t118 = qJD(5) + t123;
t143 = -t132 * t87 - t134 * t85;
t50 = -t132 * t85 + t134 * t87;
t165 = Ifges(6,4) * t50;
t45 = Ifges(6,4) * t143;
t20 = Ifges(6,1) * t50 + Ifges(6,5) * t118 + t45;
t29 = -qJD(5) * t50 + t132 * t79 + t134 * t80;
t26 = Ifges(6,6) * t29;
t28 = qJD(5) * t143 + t132 * t80 - t134 * t79;
t27 = Ifges(6,5) * t28;
t119 = qJ(2) * t154 + qJD(3);
t102 = pkin(3) * t145 + t119;
t58 = pkin(4) * t85 + t102;
t164 = t132 * t22;
t6 = t134 * t18 - t164;
t186 = t1 + t26 + t27 - (Ifges(6,5) * t143 - Ifges(6,6) * t50) * t118 / 0.2e1 + (t143 * t6 + t50 * t7) * mrSges(6,3) - (-Ifges(6,2) * t50 + t20 + t45) * t143 / 0.2e1 - t58 * (mrSges(6,1) * t50 + mrSges(6,2) * t143) - (Ifges(6,1) * t143 - t165) * t50 / 0.2e1;
t19 = Ifges(6,2) * t143 + Ifges(6,6) * t118 + t165;
t184 = t19 / 0.2e1;
t66 = -t112 * t132 - t134 * t141;
t182 = qJD(5) * t66 + t177 * t132 + t176 * t134;
t67 = t112 * t134 - t132 * t141;
t181 = -qJD(5) * t67 - t176 * t132 + t177 * t134;
t172 = t50 / 0.2e1;
t170 = t87 / 0.2e1;
t2 = qJD(5) * t6 + t12 * t134 + t13 * t132;
t169 = t2 * mrSges(6,2);
t167 = -t131 / 0.2e1;
t166 = Ifges(5,4) * t87;
t110 = t130 * t114;
t63 = t110 + t137;
t160 = qJ(2) * t131;
t156 = t128 * t114 + t130 * t160;
t159 = t128 * t129;
t73 = -pkin(6) * t159 + t156;
t40 = t133 * t63 + t135 * t73;
t140 = (mrSges(4,1) * t128 + mrSges(4,2) * t130) * t129;
t161 = -mrSges(5,1) * t85 - mrSges(5,2) * t87 - qJD(1) * t140;
t113 = pkin(3) * t159 + t129 * qJ(2);
t126 = t129 ^ 2;
t127 = t131 ^ 2;
t155 = t126 + t127;
t152 = qJD(2) * t129;
t147 = qJD(1) * qJD(2);
t144 = qJ(2) * t147;
t39 = -t133 * t73 + t135 * t63;
t142 = t97 * t128 + t96 * t130;
t99 = t141 * t129;
t33 = -pkin(4) * t131 + pkin(7) * t99 + t39;
t98 = t112 * t129;
t34 = -pkin(7) * t98 + t40;
t10 = -t132 * t34 + t134 * t33;
t11 = t132 * t33 + t134 * t34;
t56 = t132 * t99 - t134 * t98;
t57 = -t132 * t98 - t134 * t99;
t23 = t133 * t103 + t135 * t104 + t63 * t148 - t149 * t73;
t139 = -t15 * mrSges(5,1) + t14 * mrSges(5,2) + t169;
t24 = -qJD(4) * t40 + t135 * t103 - t104 * t133;
t121 = t126 * t144;
t108 = (-mrSges(4,1) * t131 - mrSges(4,3) * t158) * qJD(1);
t107 = (mrSges(4,2) * t131 - mrSges(4,3) * t159) * qJD(1);
t81 = Ifges(5,4) * t85;
t77 = Ifges(5,5) * t79;
t76 = Ifges(5,6) * t80;
t75 = t79 * mrSges(5,2);
t74 = -pkin(4) * t90 + t152;
t71 = -t128 * t146 + t92;
t70 = mrSges(5,1) * t123 - mrSges(5,3) * t87;
t69 = -mrSges(5,2) * t123 - mrSges(5,3) * t85;
t68 = -pkin(4) * t80 + t129 * t147;
t65 = pkin(4) * t98 + t113;
t44 = Ifges(5,1) * t87 + Ifges(5,5) * t123 - t81;
t43 = -Ifges(5,2) * t85 + Ifges(5,6) * t123 + t166;
t42 = mrSges(6,1) * t118 - mrSges(6,3) * t50;
t41 = -mrSges(6,2) * t118 + mrSges(6,3) * t143;
t32 = -qJD(5) * t57 + t132 * t89 + t134 * t90;
t31 = qJD(5) * t56 + t132 * t90 - t134 * t89;
t30 = -mrSges(6,1) * t143 + mrSges(6,2) * t50;
t25 = t28 * mrSges(6,2);
t17 = pkin(7) * t89 + t24;
t16 = pkin(7) * t90 + t23;
t9 = t134 * t21 - t164;
t8 = -t132 * t21 - t162;
t5 = -qJD(5) * t11 - t132 * t16 + t134 * t17;
t4 = qJD(5) * t10 + t132 * t17 + t134 * t16;
t37 = [t65 * t25 - t113 * t75 + t31 * t20 / 0.2e1 + (-Ifges(5,1) * t89 + Ifges(5,4) * t90) * t170 + (-t14 * t98 + t15 * t99 + t35 * t89 + t36 * t90) * mrSges(5,3) - t85 * (-Ifges(5,4) * t89 + Ifges(5,2) * t90) / 0.2e1 + t102 * (-t90 * mrSges(5,1) - t89 * mrSges(5,2)) + t123 * (-Ifges(5,5) * t89 + Ifges(5,6) * t90) / 0.2e1 + t32 * t184 + (Ifges(6,1) * t31 + Ifges(6,4) * t32) * t172 + m(4) * (t97 * t156 + t72 * t104 + t96 * (-t128 * t160 + t110) + t71 * t103 + t121 + t119 * t152) + m(5) * (t14 * t40 + t15 * t39 + t23 * t36 + t24 * t35 + (qJD(1) * t113 + t102) * t152) + 0.2e1 * m(3) * (t127 * t144 + t121) + (t2 * t56 - t3 * t57 - t6 * t31 + t7 * t32) * mrSges(6,3) + (-t27 / 0.2e1 - t26 / 0.2e1 - t1 + t77 / 0.2e1 - t76 / 0.2e1 + t97 * mrSges(4,2) - t96 * mrSges(4,1) + t139) * t131 + 0.2e1 * t155 * mrSges(3,3) * t147 + (-t142 * mrSges(4,3) + ((mrSges(5,1) * t98 - mrSges(5,2) * t99 + t140) * qJD(1) - t161) * qJD(2)) * t129 + (-t113 * mrSges(5,1) + t40 * mrSges(5,3) - Ifges(5,4) * t99 - Ifges(5,2) * t98 + Ifges(5,6) * t167) * t80 - (-t39 * mrSges(5,3) - Ifges(5,1) * t99 - Ifges(5,4) * t98 + Ifges(5,5) * t167) * t79 + t143 * (Ifges(6,4) * t31 + Ifges(6,2) * t32) / 0.2e1 + (-t65 * mrSges(6,1) + t11 * mrSges(6,3) + Ifges(6,4) * t57 + Ifges(6,2) * t56 + Ifges(6,6) * t167) * t29 + (-t10 * mrSges(6,3) + Ifges(6,1) * t57 + Ifges(6,4) * t56 + Ifges(6,5) * t167) * t28 + t4 * t41 + t5 * t42 + t58 * (-mrSges(6,1) * t32 + mrSges(6,2) * t31) + t68 * (-mrSges(6,1) * t56 + mrSges(6,2) * t57) + t23 * t69 + t24 * t70 + t74 * t30 - t89 * t44 / 0.2e1 + t90 * t43 / 0.2e1 + t104 * t107 + t103 * t108 + t118 * (Ifges(6,5) * t31 + Ifges(6,6) * t32) / 0.2e1 + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + t58 * t74 + t65 * t68); t177 * t70 + t176 * t69 + t181 * t42 + t182 * t41 + (-t28 * t66 + t29 * t67) * mrSges(6,3) + (t112 * t80 - t141 * t79) * mrSges(5,3) + (-t58 * t154 + t181 * t6 + t182 * t7 + t2 * t67 + t3 * t66) * m(6) + (-t102 * t154 + t112 * t14 - t141 * t15 + t176 * t36 + t177 * t35) * m(5) + m(4) * t142 + ((-t107 * t130 + t108 * t128) * t131 + (-t30 + t161) * t129 + (-t119 * t129 - (-t128 * t71 + t130 * t72) * t131) * m(4) + (-m(3) * qJ(2) - mrSges(3,3)) * t155 * qJD(1)) * qJD(1); -t80 * mrSges(5,1) - t29 * mrSges(6,1) - t143 * t41 + t50 * t42 + t85 * t69 + t87 * t70 + t25 - t75 - m(5) * (-t35 * t87 - t36 * t85) + (m(5) * qJD(2) + t128 * t107 + t130 * t108 + (t128 * t72 + t130 * t71 + qJD(2)) * m(4)) * t154 + (-t143 * t7 + t50 * t6 + t68) * m(6); -t139 - t123 * (-Ifges(5,5) * t85 - Ifges(5,6) * t87) / 0.2e1 + (-Ifges(5,2) * t87 + t44 - t81) * t85 / 0.2e1 + t50 * t184 + (-t35 * t85 + t36 * t87) * mrSges(5,3) - t87 * (-Ifges(5,1) * t85 - t166) / 0.2e1 - t102 * (mrSges(5,1) * t87 - mrSges(5,2) * t85) + t43 * t170 - m(6) * (t6 * t8 + t7 * t9) + t76 - t77 - t9 * t41 - t8 * t42 - t35 * t69 + t36 * t70 + (-t87 * t30 + (-t132 * t42 + t134 * t41) * qJD(5) + (t132 * t29 - t134 * t28) * mrSges(6,3) + (t132 * t2 + t134 * t3 - t58 * t87 + (-t132 * t6 + t134 * t7) * qJD(5)) * m(6)) * pkin(4) + t186; t19 * t172 - t6 * t41 + t7 * t42 - t169 + t186;];
tauc = t37(:);
