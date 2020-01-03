% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR7
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:57
% EndTime: 2019-12-31 17:06:04
% DurationCPUTime: 2.49s
% Computational Cost: add. (1721->286), mult. (4655->423), div. (0->0), fcn. (3033->6), ass. (0->148)
t134 = cos(pkin(7));
t87 = sin(pkin(7));
t89 = sin(qJ(2));
t91 = cos(qJ(2));
t98 = t134 * t91 - t87 * t89;
t62 = t98 * qJD(1);
t177 = -t62 / 0.2e1;
t60 = qJD(4) - t62;
t176 = -qJD(2) / 0.2e1;
t115 = t134 * t89;
t129 = qJD(1) * t91;
t63 = -qJD(1) * t115 - t129 * t87;
t142 = t63 * mrSges(4,3);
t88 = sin(qJ(4));
t90 = cos(qJ(4));
t47 = qJD(2) * t90 + t63 * t88;
t48 = qJD(2) * t88 - t63 * t90;
t135 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t47 + mrSges(5,2) * t48 - t142;
t132 = Ifges(3,6) * qJD(2);
t152 = Ifges(3,4) * t89;
t175 = t132 / 0.2e1 + (t91 * Ifges(3,2) + t152) * qJD(1) / 0.2e1 + pkin(5) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t129);
t136 = -qJ(3) - pkin(5);
t114 = qJD(2) * t136;
t61 = qJD(3) * t91 + t114 * t89;
t52 = t61 * qJD(1);
t94 = -t89 * qJD(3) + t114 * t91;
t93 = t87 * t94;
t23 = qJD(1) * t93 + t134 * t52;
t125 = qJD(1) * qJD(2);
t119 = t89 * t125;
t113 = pkin(2) * t119;
t72 = t87 * t91 + t115;
t64 = t72 * qJD(2);
t56 = qJD(1) * t64;
t65 = t98 * qJD(2);
t57 = qJD(1) * t65;
t26 = pkin(3) * t56 - pkin(6) * t57 + t113;
t85 = -pkin(2) * t91 - pkin(1);
t131 = qJD(1) * t85;
t76 = qJD(3) + t131;
t29 = -t62 * pkin(3) + t63 * pkin(6) + t76;
t79 = t136 * t91;
t75 = qJD(1) * t79;
t66 = t134 * t75;
t121 = t136 * t89;
t74 = qJD(1) * t121;
t70 = qJD(2) * pkin(2) + t74;
t41 = t87 * t70 - t66;
t37 = qJD(2) * pkin(6) + t41;
t8 = t29 * t90 - t37 * t88;
t1 = qJD(4) * t8 + t23 * t90 + t26 * t88;
t9 = t29 * t88 + t37 * t90;
t2 = -qJD(4) * t9 - t23 * t88 + t26 * t90;
t24 = qJD(4) * t47 + t57 * t90;
t25 = -qJD(4) * t48 - t57 * t88;
t174 = t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t24 + Ifges(5,6) * t25;
t140 = t63 * Ifges(4,4);
t173 = t140 / 0.2e1 + Ifges(4,2) * t177;
t102 = Ifges(5,5) * t90 - Ifges(5,6) * t88;
t149 = Ifges(5,4) * t90;
t104 = -Ifges(5,2) * t88 + t149;
t150 = Ifges(5,4) * t88;
t106 = Ifges(5,1) * t90 - t150;
t46 = Ifges(5,4) * t47;
t14 = Ifges(5,1) * t48 + Ifges(5,5) * t60 + t46;
t137 = t90 * t14;
t151 = Ifges(5,4) * t48;
t13 = Ifges(5,2) * t47 + Ifges(5,6) * t60 + t151;
t138 = t88 * t13;
t164 = t48 / 0.2e1;
t107 = mrSges(5,1) * t88 + mrSges(5,2) * t90;
t139 = t87 * t75;
t40 = t134 * t70 + t139;
t36 = -qJD(2) * pkin(3) - t40;
t99 = t36 * t107;
t172 = t137 / 0.2e1 - t138 / 0.2e1 + t60 * t102 / 0.2e1 + t106 * t164 + t47 * t104 / 0.2e1 + t99;
t171 = -0.2e1 * pkin(1);
t169 = t24 / 0.2e1;
t168 = t25 / 0.2e1;
t167 = Ifges(4,6) * t176 + t173;
t166 = -t47 / 0.2e1;
t165 = -t48 / 0.2e1;
t163 = t56 / 0.2e1;
t162 = -t60 / 0.2e1;
t159 = -t88 / 0.2e1;
t158 = t90 / 0.2e1;
t157 = pkin(2) * t87;
t130 = qJD(1) * t89;
t156 = pkin(5) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t130);
t154 = mrSges(5,3) * t88;
t153 = mrSges(5,3) * t90;
t59 = Ifges(4,4) * t62;
t92 = t134 * t94;
t22 = -qJD(1) * t92 + t52 * t87;
t44 = -t115 * t136 - t79 * t87;
t148 = t22 * t44;
t147 = t47 * Ifges(5,6);
t146 = t48 * Ifges(5,5);
t145 = t60 * Ifges(5,3);
t144 = t62 * mrSges(4,3);
t141 = t63 * Ifges(4,1);
t133 = Ifges(3,5) * qJD(2);
t128 = qJD(2) * t89;
t127 = qJD(4) * t88;
t126 = qJD(4) * t90;
t124 = pkin(2) * t130;
t120 = t134 * pkin(2);
t118 = t91 * t125;
t116 = t56 * mrSges(4,1) + t57 * mrSges(4,2);
t111 = -t1 * t88 - t2 * t90;
t110 = -t8 * t90 - t9 * t88;
t109 = t8 * t88 - t9 * t90;
t108 = mrSges(5,1) * t90 - mrSges(5,2) * t88;
t105 = Ifges(5,1) * t88 + t149;
t103 = Ifges(5,2) * t90 + t150;
t101 = Ifges(5,5) * t88 + Ifges(5,6) * t90;
t27 = -mrSges(5,2) * t60 + mrSges(5,3) * t47;
t28 = mrSges(5,1) * t60 - mrSges(5,3) * t48;
t100 = t27 * t90 - t28 * t88;
t39 = -pkin(3) * t98 - t72 * pkin(6) + t85;
t45 = t121 * t87 - t134 * t79;
t17 = t39 * t90 - t45 * t88;
t18 = t39 * t88 + t45 * t90;
t86 = Ifges(3,4) * t129;
t84 = -t120 - pkin(3);
t83 = pkin(6) + t157;
t69 = Ifges(3,1) * t130 + t133 + t86;
t55 = Ifges(5,3) * t56;
t50 = -qJD(2) * mrSges(4,2) + t144;
t43 = t134 * t74 + t139;
t42 = t74 * t87 - t66;
t38 = -mrSges(4,1) * t62 - mrSges(4,2) * t63;
t35 = Ifges(4,5) * qJD(2) - t141 + t59;
t33 = pkin(2) * t128 + pkin(3) * t64 - pkin(6) * t65;
t32 = -pkin(3) * t63 - pkin(6) * t62 + t124;
t31 = t134 * t61 + t93;
t30 = t61 * t87 - t92;
t16 = -mrSges(5,2) * t56 + mrSges(5,3) * t25;
t15 = mrSges(5,1) * t56 - mrSges(5,3) * t24;
t12 = t145 + t146 + t147;
t11 = t32 * t88 + t43 * t90;
t10 = t32 * t90 - t43 * t88;
t7 = -mrSges(5,1) * t25 + mrSges(5,2) * t24;
t6 = t24 * Ifges(5,1) + t25 * Ifges(5,4) + t56 * Ifges(5,5);
t5 = t24 * Ifges(5,4) + t25 * Ifges(5,2) + t56 * Ifges(5,6);
t4 = -qJD(4) * t18 - t31 * t88 + t33 * t90;
t3 = qJD(4) * t17 + t31 * t90 + t33 * t88;
t19 = [t85 * t116 + t31 * t50 + t44 * t7 + t3 * t27 + t4 * t28 + t17 * t15 + t18 * t16 + t135 * t30 + m(4) * (t23 * t45 - t30 * t40 + t31 * t41 + t148) + m(5) * (t1 * t18 + t17 * t2 + t3 * t9 + t30 * t36 + t4 * t8 + t148) + (t145 / 0.2e1 + t147 / 0.2e1 + t146 / 0.2e1 + t8 * mrSges(5,1) - t9 * mrSges(5,2) + t76 * mrSges(4,1) + t12 / 0.2e1 + t167 + t173) * t64 + (t76 * mrSges(4,2) - t141 / 0.2e1 + t59 / 0.2e1 + t35 / 0.2e1 + t110 * mrSges(5,3) + t172) * t65 + (-t40 * t65 - t41 * t64 + t44 * t57 - t45 * t56) * mrSges(4,3) + (Ifges(4,5) * t65 / 0.2e1 - Ifges(4,6) * t64 / 0.2e1 + (t69 / 0.2e1 - t156 + t133 / 0.2e1 + (mrSges(3,2) * t171 + 0.3e1 / 0.2e1 * Ifges(3,4) * t91) * qJD(1)) * t91) * qJD(2) - (t55 / 0.2e1 - Ifges(4,4) * t57 - t23 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t56 + t174) * t98 + (t102 * t163 + t106 * t169 + t104 * t168 - Ifges(4,4) * t56 + Ifges(4,1) * t57 + t6 * t158 + t5 * t159 + (mrSges(4,3) + t107) * t22 + t111 * mrSges(5,3) + (t101 * t162 + t103 * t166 + t105 * t165 + t36 * t108 - t90 * t13 / 0.2e1 + t14 * t159 + t109 * mrSges(5,3)) * qJD(4)) * t72 + (-t132 / 0.2e1 + (mrSges(3,1) * t171 - 0.3e1 / 0.2e1 * t152 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t91) * qJD(1) + (m(4) * (t76 + t131) + t38 + qJD(1) * (-mrSges(4,1) * t98 + mrSges(4,2) * t72)) * pkin(2) - t175) * t128; (-t76 * t124 + t40 * t42 - t41 * t43 + (-t134 * t22 + t23 * t87) * pkin(2)) * m(4) - t9 * (mrSges(5,2) * t63 - t154 * t62) - t2 * t154 - t8 * (-mrSges(5,1) * t63 - t153 * t62) - t41 * t142 + t62 * t138 / 0.2e1 - t38 * t124 - (Ifges(3,5) * t91 - Ifges(3,6) * t89) * t125 / 0.2e1 - m(5) * (t10 * t8 + t11 * t9) + (m(5) * t84 - mrSges(4,1) - t108) * t22 - Ifges(3,6) * t119 + (-t126 * t8 - t127 * t9) * mrSges(5,3) + (m(5) * (t1 * t90 - t2 * t88) + t90 * t16 - t88 * t15 - t28 * t126 - t27 * t127) * t83 + t88 * t6 / 0.2e1 - t76 * (-mrSges(4,1) * t63 + mrSges(4,2) * t62) + t84 * t7 - (-Ifges(3,2) * t130 + t69 + t86) * t129 / 0.2e1 + (Ifges(4,1) * t62 + t12 + t140) * t63 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t89 + mrSges(3,2) * t91) - t89 * (Ifges(3,1) * t91 - t152) / 0.2e1) * qJD(1) ^ 2 - Ifges(4,6) * t56 + Ifges(4,5) * t57 - t43 * t50 - t23 * mrSges(4,2) - t11 * t27 - t10 * t28 + (-m(5) * t36 - t135) * t42 + t175 * t130 + (m(5) * t110 * t83 + t172) * qJD(4) + (-t120 * t57 - t157 * t56) * mrSges(4,3) + t5 * t158 + (-Ifges(5,3) * t63 + t102 * t62) * t162 + t101 * t163 + (-Ifges(5,5) * t63 + t106 * t62) * t165 + (-Ifges(5,6) * t63 + t104 * t62) * t166 + t63 * t167 + t103 * t168 + t105 * t169 + (Ifges(4,2) * t63 + t137 + t35 + t59) * t177 + t1 * t153 + t129 * t156 + t40 * t144 - t62 * t99 + (-mrSges(3,1) * t118 + mrSges(3,2) * t119) * pkin(5) + (Ifges(4,5) * t62 + Ifges(4,6) * t63) * t176 + Ifges(3,5) * t118; t90 * t15 + t88 * t16 + t135 * t63 + t100 * qJD(4) + (-t100 - t50) * t62 + t116 + (-t60 * t109 + t36 * t63 - t111) * m(5) + (-t40 * t63 - t41 * t62 + t113) * m(4); t55 - t36 * (mrSges(5,1) * t48 + mrSges(5,2) * t47) + (Ifges(5,1) * t47 - t151) * t165 + t13 * t164 + (Ifges(5,5) * t47 - Ifges(5,6) * t48) * t162 - t8 * t27 + t9 * t28 + (t47 * t8 + t48 * t9) * mrSges(5,3) + (-Ifges(5,2) * t48 + t14 + t46) * t166 + t174;];
tauc = t19(:);
