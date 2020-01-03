% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:25
% DurationCPUTime: 2.55s
% Computational Cost: add. (1692->276), mult. (3516->424), div. (0->0), fcn. (1815->6), ass. (0->140)
t121 = qJ(2) * qJD(1);
t74 = -pkin(1) - pkin(2);
t59 = qJD(1) * t74 + qJD(2);
t67 = sin(pkin(8));
t68 = cos(pkin(8));
t42 = t68 * t121 + t67 * t59;
t36 = -qJD(1) * pkin(6) + t42;
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t24 = qJD(3) * t73 - t36 * t71;
t21 = -qJD(4) * pkin(4) - t24;
t176 = m(6) * t21;
t133 = qJD(1) * t71;
t108 = mrSges(5,3) * t133;
t72 = cos(qJ(5));
t128 = qJD(4) * t72;
t70 = sin(qJ(5));
t50 = t133 * t70 + t128;
t122 = t70 * qJD(4);
t51 = t133 * t72 - t122;
t137 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t50 - mrSges(6,2) * t51 - t108;
t129 = qJD(4) * t71;
t110 = t67 * t129;
t141 = t70 * t73;
t138 = t72 * t73;
t45 = t138 * t67 - t68 * t70;
t175 = -qJD(5) * t45 + t110 * t70 - (-t141 * t68 + t67 * t72) * qJD(1);
t44 = -t141 * t67 - t68 * t72;
t174 = qJD(5) * t44 - t110 * t72 - (t138 * t68 + t67 * t70) * qJD(1);
t25 = qJD(3) * t71 + t36 * t73;
t123 = t25 * qJD(4);
t124 = t24 * qJD(4);
t120 = qJD(1) * qJD(2);
t106 = t68 * t120;
t14 = t106 * t73 + t124;
t15 = t106 * t71 + t123;
t149 = t15 * t71;
t87 = t14 * t73 + t149;
t173 = -t71 * t123 - t73 * t124 + t87;
t118 = qJD(4) * qJD(5);
t126 = qJD(5) * t71;
t127 = qJD(4) * t73;
t79 = -t126 * t70 + t127 * t72;
t31 = -qJD(1) * t79 + t118 * t72;
t165 = t31 / 0.2e1;
t109 = t73 * t122;
t125 = qJD(5) * t72;
t78 = t125 * t71 + t109;
t32 = qJD(1) * t78 - t118 * t70;
t164 = t32 / 0.2e1;
t172 = -m(3) * qJ(2) - mrSges(3,3);
t130 = qJD(2) * t68;
t136 = t68 * qJ(2) + t67 * t74;
t53 = -pkin(6) + t136;
t171 = -t53 * t129 + t73 * t130;
t131 = qJD(2) * t67;
t98 = -pkin(4) * t71 + pkin(7) * t73;
t43 = qJD(4) * t98 + t131;
t37 = t43 * qJD(1);
t22 = qJD(4) * pkin(7) + t25;
t41 = -t67 * t121 + t59 * t68;
t35 = qJD(1) * pkin(3) - t41;
t99 = pkin(4) * t73 + pkin(7) * t71;
t23 = qJD(1) * t99 + t35;
t5 = -t22 * t70 + t23 * t72;
t1 = qJD(5) * t5 + t14 * t72 + t37 * t70;
t6 = t22 * t72 + t23 * t70;
t2 = -qJD(5) * t6 - t14 * t70 + t37 * t72;
t96 = t1 * t72 - t2 * t70;
t49 = Ifges(6,4) * t50;
t132 = qJD(1) * t73;
t60 = qJD(5) + t132;
t18 = -Ifges(6,1) * t51 + Ifges(6,5) * t60 + t49;
t139 = t72 * t18;
t160 = t60 / 0.2e1;
t161 = -t51 / 0.2e1;
t162 = t50 / 0.2e1;
t153 = Ifges(6,4) * t51;
t17 = Ifges(6,2) * t50 + Ifges(6,6) * t60 - t153;
t89 = Ifges(6,5) * t72 - Ifges(6,6) * t70;
t151 = Ifges(6,4) * t72;
t91 = -Ifges(6,2) * t70 + t151;
t152 = Ifges(6,4) * t70;
t93 = Ifges(6,1) * t72 - t152;
t94 = mrSges(6,1) * t70 + mrSges(6,2) * t72;
t157 = t6 * t70;
t95 = -t5 * t72 - t157;
t168 = t89 * t160 + t93 * t161 + t91 * t162 + t21 * t94 + t95 * mrSges(6,3) - t70 * t17 / 0.2e1 + t139 / 0.2e1;
t119 = qJD(4) * qJD(1);
t105 = t71 * t119;
t166 = Ifges(6,4) * t165 + Ifges(6,2) * t164 - Ifges(6,6) * t105 / 0.2e1;
t155 = Ifges(5,4) * t71;
t154 = Ifges(5,4) * t73;
t148 = t25 * t73;
t147 = t50 * Ifges(6,6);
t146 = t51 * Ifges(6,5);
t145 = t53 * t71;
t144 = t60 * Ifges(6,3);
t143 = t68 * t71;
t142 = t70 * t71;
t140 = t71 * t72;
t135 = Ifges(5,5) * qJD(4);
t134 = Ifges(5,6) * qJD(4);
t117 = mrSges(5,3) * t132;
t100 = -t67 * qJ(2) + t68 * t74;
t52 = pkin(3) - t100;
t97 = -qJD(5) * t53 * t73 + t43;
t92 = Ifges(6,1) * t70 + t151;
t90 = Ifges(6,2) * t72 + t152;
t88 = Ifges(6,5) * t70 + Ifges(6,6) * t72;
t19 = -mrSges(6,1) * t105 - mrSges(6,3) * t31;
t20 = mrSges(6,2) * t105 + mrSges(6,3) * t32;
t86 = -t70 * t19 + t72 * t20;
t33 = -mrSges(6,2) * t60 + mrSges(6,3) * t50;
t34 = mrSges(6,1) * t60 + mrSges(6,3) * t51;
t85 = -t70 * t33 - t72 * t34;
t84 = -t41 * t67 + t42 * t68;
t54 = (mrSges(5,1) * t73 - mrSges(5,2) * t71) * qJD(1);
t83 = (-mrSges(5,1) * t71 - mrSges(5,2) * t73) * qJD(4);
t82 = (-t71 * Ifges(5,1) - t154) * qJD(1);
t81 = (-t73 * Ifges(5,2) - t155) * qJD(1);
t80 = Ifges(6,5) * t31 + Ifges(6,6) * t32 - Ifges(6,3) * t105;
t38 = t52 + t99;
t77 = qJD(5) * t38 + t171;
t75 = qJD(1) ^ 2;
t58 = -qJD(4) * mrSges(5,2) - t117;
t55 = t98 * qJD(1);
t48 = qJD(1) * t83;
t47 = t82 + t135;
t46 = t81 + t134;
t16 = t144 - t146 + t147;
t13 = t138 * t53 + t38 * t70;
t12 = -t141 * t53 + t38 * t72;
t11 = t24 * t72 + t55 * t70;
t10 = -t24 * t70 + t55 * t72;
t9 = -mrSges(6,1) * t32 + mrSges(6,2) * t31;
t8 = Ifges(6,1) * t31 + Ifges(6,4) * t32 - Ifges(6,5) * t105;
t4 = -t70 * t77 + t72 * t97;
t3 = t70 * t97 + t72 * t77;
t7 = [(-t137 - t176) * (-t53 * t127 - t71 * t130) - t173 * mrSges(5,3) - t73 * (Ifges(5,2) * t71 - t154) * t119 + m(4) * ((-t100 * t67 + t136 * t68) * qJD(1) + t84) * qJD(2) + t5 * (-mrSges(6,1) * t129 + mrSges(6,3) * t79) + t6 * (mrSges(6,2) * t129 + mrSges(6,3) * t78) + 0.2e1 * t54 * t131 + t17 * t109 / 0.2e1 + 0.2e1 * mrSges(4,2) * t106 + t73 * t80 / 0.2e1 + t21 * (-mrSges(6,1) * t78 - mrSges(6,2) * t79) + qJD(4) ^ 2 * (-Ifges(5,5) * t73 + Ifges(5,6) * t71) / 0.2e1 + t171 * t58 + 0.2e1 * (mrSges(4,1) * t67 - t172) * t120 - t94 * t149 - (-Ifges(5,1) * t73 + t155) * t105 + t52 * t48 + t3 * t33 + t4 * t34 + t12 * t19 + t13 * t20 + m(5) * (((-t24 * t73 - t25 * t71) * qJD(4) + t87) * t53 + ((-t24 * t71 + t148) * t68 + (qJD(1) * t52 + t35) * t67) * qJD(2)) - t8 * t140 / 0.2e1 + t2 * (mrSges(6,1) * t73 + mrSges(6,3) * t140) + t1 * (-mrSges(6,2) * t73 + mrSges(6,3) * t142) + m(6) * (t1 * t13 + t12 * t2 + t15 * t145 + t3 * t6 + t4 * t5) + (t72 * t17 + t70 * t18) * t126 / 0.2e1 - (t82 + t47 + t139) * t127 / 0.2e1 + (t81 + t46) * t129 / 0.2e1 - (qJD(1) * (t73 * Ifges(6,3) - t71 * t89) + t16) * t129 / 0.2e1 + t35 * t83 + t9 * t145 + (t88 * t126 + (-Ifges(6,3) * t71 - t73 * t89) * qJD(4)) * t160 + (t92 * t126 + (-Ifges(6,5) * t71 - t73 * t93) * qJD(4)) * t161 + (t90 * t126 + (-Ifges(6,6) * t71 - t73 * t91) * qJD(4)) * t162 + (Ifges(6,6) * t73 - t71 * t91) * t164 + (Ifges(6,5) * t73 - t71 * t93) * t165 + t142 * t166; t44 * t19 + t45 * t20 + t172 * t75 + t175 * t34 + t174 * t33 + (-t75 * mrSges(4,2) - t48) * t68 + (t1 * t45 + t174 * t6 + t175 * t5 + t2 * t44) * m(6) + (-t75 * mrSges(4,1) + t71 * t9 + (t137 * t73 - t58 * t71) * qJD(4) + m(5) * t173 + m(6) * (t127 * t21 + t149)) * t67 + ((-t137 * t71 - t73 * t58) * t68 - t67 * t54 - m(4) * t84 - t143 * t176 - (-t24 * t143 + t68 * t148 + (t130 + t35) * t67) * m(5)) * qJD(1); (-t9 + (t33 * t72 - t34 * t70 + t117 + t58) * qJD(4) + m(5) * (-t15 + t123) + m(6) * (-t122 * t5 + t128 * t6 - t15)) * t73 + (t85 * qJD(5) + (t108 + t137) * qJD(4) + m(5) * (t14 - t124) + m(6) * (qJD(4) * t21 - qJD(5) * t157 - t125 * t5 + t96) + t86) * t71; t92 * t165 + t90 * t164 + t72 * t166 + t70 * t8 / 0.2e1 - t24 * t58 - t11 * t33 - t10 * t34 - pkin(4) * t9 - t14 * mrSges(5,2) - t137 * t25 + (-t72 * mrSges(6,1) + t70 * mrSges(6,2) - mrSges(5,1)) * t15 + t96 * mrSges(6,3) + t168 * qJD(5) + ((-t46 / 0.2e1 + t16 / 0.2e1 + t144 / 0.2e1 + t147 / 0.2e1 - t146 / 0.2e1 + t5 * mrSges(6,1) - t6 * mrSges(6,2) + t35 * mrSges(5,1) + t134 / 0.2e1 - qJD(4) * t88 / 0.2e1 - t25 * mrSges(5,3) + Ifges(5,4) * t133 / 0.2e1) * t71 + (t47 / 0.2e1 + t35 * mrSges(5,2) - t135 / 0.2e1 - t24 * mrSges(5,3) + (-t154 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t71) * qJD(1) + t168) * t73) * qJD(1) + (-pkin(4) * t15 - t10 * t5 - t11 * t6 - t21 * t25) * m(6) + (t86 + m(6) * t96 + (m(6) * t95 + t85) * qJD(5)) * pkin(7); -t1 * mrSges(6,2) + t2 * mrSges(6,1) - t21 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + t51 * (Ifges(6,1) * t50 + t153) / 0.2e1 + t17 * t161 - t60 * (Ifges(6,5) * t50 + Ifges(6,6) * t51) / 0.2e1 - t5 * t33 + t6 * t34 + (t5 * t50 - t51 * t6) * mrSges(6,3) + t80 - (Ifges(6,2) * t51 + t18 + t49) * t50 / 0.2e1;];
tauc = t7(:);
