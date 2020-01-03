% Calculate vector of inverse dynamics joint torques for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:31
% DurationCPUTime: 5.19s
% Computational Cost: add. (1131->294), mult. (2671->381), div. (0->0), fcn. (1608->8), ass. (0->129)
t171 = m(4) + m(5);
t176 = mrSges(4,1) + mrSges(5,1);
t175 = mrSges(4,2) - mrSges(5,3);
t169 = Ifges(4,1) + Ifges(5,1);
t167 = Ifges(5,4) + Ifges(4,5);
t168 = -Ifges(4,4) + Ifges(5,5);
t86 = sin(qJ(2));
t88 = cos(qJ(2));
t71 = -mrSges(3,1) * t88 + mrSges(3,2) * t86;
t83 = qJ(2) + pkin(6);
t81 = sin(t83);
t82 = cos(t83);
t174 = t175 * t81 - t176 * t82 + t71;
t87 = sin(qJ(1));
t142 = g(2) * t87;
t89 = cos(qJ(1));
t173 = g(1) * t89 + t142;
t127 = cos(pkin(6));
t107 = t127 * t86;
t84 = sin(pkin(6));
t62 = t84 * t88 + t107;
t148 = t62 / 0.2e1;
t120 = qJD(1) * qJD(2);
t112 = t86 * t120;
t121 = qJDD(1) * t88;
t66 = -t112 + t121;
t172 = t66 / 0.2e1;
t170 = qJD(2) / 0.2e1;
t166 = Ifges(5,6) - Ifges(4,6);
t106 = t127 * t88;
t123 = t86 * qJD(1);
t49 = -qJD(1) * t106 + t123 * t84;
t138 = Ifges(5,5) * t49;
t47 = Ifges(4,4) * t49;
t51 = t62 * qJD(1);
t165 = t167 * qJD(2) + t169 * t51 + t138 - t47;
t136 = t49 * mrSges(4,3);
t137 = t49 * mrSges(5,2);
t42 = qJD(2) * mrSges(5,3) - t137;
t164 = -qJD(2) * mrSges(4,2) - t136 + t42;
t134 = t51 * mrSges(4,3);
t135 = t51 * mrSges(5,2);
t163 = qJD(2) * t176 - t134 - t135;
t161 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(3,3);
t79 = pkin(5) * t121;
t59 = -pkin(5) * t112 + t79;
t67 = qJDD(1) * t86 + t120 * t88;
t60 = t67 * pkin(5);
t159 = t59 * t88 + t60 * t86;
t156 = m(3) * pkin(1) + mrSges(2,1) - t174;
t153 = -t49 / 0.2e1;
t152 = t49 / 0.2e1;
t150 = t51 / 0.2e1;
t147 = pkin(2) * t84;
t145 = pkin(2) * t88;
t144 = pkin(5) * t88;
t141 = Ifges(3,4) * t86;
t140 = Ifges(3,4) * t88;
t139 = Ifges(4,4) * t51;
t85 = -qJ(3) - pkin(5);
t72 = t85 * t88;
t65 = qJD(1) * t72;
t131 = t84 * t65;
t23 = qJDD(2) * pkin(2) - qJ(3) * t67 - qJD(3) * t123 - t60;
t126 = qJD(2) * t86;
t116 = pkin(5) * t126;
t125 = qJD(3) * t88;
t29 = qJ(3) * t66 + t79 + (-t116 + t125) * qJD(1);
t5 = t127 * t29 + t84 * t23;
t53 = t127 * t65;
t114 = t85 * t86;
t64 = qJD(1) * t114;
t58 = qJD(2) * pkin(2) + t64;
t28 = t84 * t58 - t53;
t124 = qJDD(1) * pkin(1);
t122 = t88 * qJD(1);
t119 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t123) * t144;
t118 = pkin(2) * t123;
t117 = pkin(2) * t126;
t78 = pkin(1) + t145;
t113 = t127 * pkin(2);
t35 = -t127 * t66 + t67 * t84;
t36 = t127 * t67 + t84 * t66;
t110 = t35 * mrSges(4,1) + t36 * mrSges(4,2);
t109 = t35 * mrSges(5,1) - t36 * mrSges(5,3);
t105 = qJD(2) * t85;
t20 = -qJDD(2) * mrSges(5,1) + t36 * mrSges(5,2);
t103 = -g(1) * t87 + g(2) * t89;
t102 = mrSges(3,1) * t86 + mrSges(3,2) * t88;
t99 = t88 * Ifges(3,2) + t141;
t98 = Ifges(3,5) * t88 - Ifges(3,6) * t86;
t97 = pkin(3) * t82 + qJ(4) * t81;
t96 = pkin(1) * t102;
t95 = t86 * (Ifges(3,1) * t88 - t141);
t45 = -pkin(2) * t66 + qJDD(3) - t124;
t4 = t127 * t23 - t84 * t29;
t27 = t127 * t58 + t131;
t94 = -t84 * t86 + t106;
t68 = -qJD(1) * t78 + qJD(3);
t92 = -qJD(3) * t86 + t105 * t88;
t80 = Ifges(3,4) * t122;
t77 = -t113 - pkin(3);
t75 = qJ(4) + t147;
t70 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t122;
t56 = Ifges(3,1) * t123 + Ifges(3,5) * qJD(2) + t80;
t55 = Ifges(3,6) * qJD(2) + qJD(1) * t99;
t52 = t94 * qJD(2);
t50 = t62 * qJD(2);
t48 = t105 * t86 + t125;
t46 = Ifges(5,5) * t51;
t34 = t127 * t64 + t131;
t33 = t64 * t84 - t53;
t26 = -pkin(3) * t94 - qJ(4) * t62 - t78;
t25 = mrSges(4,1) * t49 + mrSges(4,2) * t51;
t24 = mrSges(5,1) * t49 - mrSges(5,3) * t51;
t21 = -mrSges(5,2) * t35 + qJDD(2) * mrSges(5,3);
t19 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t36;
t18 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t35;
t17 = qJD(2) * qJ(4) + t28;
t14 = -Ifges(4,2) * t49 + Ifges(4,6) * qJD(2) + t139;
t13 = Ifges(5,6) * qJD(2) + Ifges(5,3) * t49 + t46;
t12 = -qJD(2) * pkin(3) + qJD(4) - t27;
t8 = pkin(3) * t51 + qJ(4) * t49 + t118;
t7 = pkin(3) * t49 - qJ(4) * t51 + t68;
t6 = pkin(3) * t50 - qJ(4) * t52 - qJD(4) * t62 + t117;
t3 = -qJDD(2) * pkin(3) + qJDD(4) - t4;
t2 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t5;
t1 = pkin(3) * t35 - qJ(4) * t36 - qJD(4) * t51 + t45;
t9 = [-qJDD(2) * mrSges(3,2) * t144 + t94 * (Ifges(4,4) * t36 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (-t168 * t94 + t169 * t62) * t36 / 0.2e1 - t94 * (Ifges(5,5) * t36 + Ifges(5,6) * qJDD(2)) / 0.2e1 + t1 * (-mrSges(5,1) * t94 - mrSges(5,3) * t62) + t45 * (-mrSges(4,1) * t94 + mrSges(4,2) * t62) + (Ifges(3,1) * t67 + Ifges(3,4) * t172 - pkin(5) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t67)) * t86 + t88 * (Ifges(3,4) * t67 + Ifges(3,2) * t66) / 0.2e1 + (t165 / 0.2e1 + t68 * mrSges(4,2) - t7 * mrSges(5,3) + Ifges(4,4) * t153 + Ifges(5,5) * t152 + t167 * t170 + t169 * t150 - t27 * mrSges(4,3) + t12 * mrSges(5,2)) * t52 + t26 * t109 - t78 * t110 + qJDD(2) * Ifges(3,6) * t88 + (-m(4) * t27 + m(5) * t12 - t163) * (-t127 * t92 + t48 * t84) + (m(4) * t28 + m(5) * t17 + t164) * (t127 * t48 + t84 * t92) + (t144 * t66 + t159) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t159) + (t2 * t94 + t3 * t62 - t142) * mrSges(5,2) + (0.2e1 * Ifges(3,5) * t86 - t166 * t94 + t167 * t62) * qJDD(2) / 0.2e1 + (-t4 * t62 + t5 * t94 - t142) * mrSges(4,3) - pkin(1) * (-mrSges(3,1) * t66 + mrSges(3,2) * t67) + t6 * t24 - t71 * t124 - t55 * t126 / 0.2e1 + (-(Ifges(5,3) + Ifges(4,2)) * t94 + 0.2e1 * t168 * t148) * t35 + t99 * t172 + (qJDD(2) * t167 + t169 * t36) * t148 + (t170 * t98 - t119) * qJD(2) + t25 * t117 + ((t171 * t85 - mrSges(5,2) - mrSges(4,3) + t161) * t89 + (-m(5) * (-t78 - t97) + m(4) * t78 + t156) * t87) * g(1) + (t161 * t87 - t171 * (t89 * t78 - t87 * t85) + (-m(5) * t97 - t156) * t89) * g(2) + (t95 + t88 * (-Ifges(3,2) * t86 + t140)) * t120 / 0.2e1 + (t166 * t170 + t168 * t150 + t68 * mrSges(4,1) + t7 * mrSges(5,1) + t13 / 0.2e1 - t14 / 0.2e1 + Ifges(5,3) * t152 - Ifges(4,2) * t153 - t28 * mrSges(4,3) - t17 * mrSges(5,2)) * t50 + t67 * t140 / 0.2e1 + Ifges(2,3) * qJDD(1) + m(5) * (t1 * t26 + t6 * t7) + m(4) * (t117 * t68 - t45 * t78) + t88 * t56 * t170 + (-m(4) * t4 + m(5) * t3 - t19 + t20) * (-t107 * t85 - t72 * t84) + (m(4) * t5 + m(5) * t2 + t18 + t21) * (t114 * t84 - t127 * t72) - t70 * t116 - t96 * t120; (pkin(5) * t70 + t55 / 0.2e1) * t123 + ((t127 * t4 + t5 * t84) * pkin(2) - t68 * t118 + t27 * t33 - t28 * t34) * m(4) + t163 * t33 - t164 * t34 + (-Ifges(4,2) * t51 + t165 - t47) * t152 + t166 * t35 + t167 * t36 - (t166 * t51 - t167 * t49) * qJD(2) / 0.2e1 + (t119 + (-t95 / 0.2e1 + t96) * qJD(1)) * qJD(1) + t75 * t21 + t77 * t20 + Ifges(3,6) * t66 + Ifges(3,5) * t67 - t68 * (mrSges(4,1) * t51 - mrSges(4,2) * t49) - t59 * mrSges(3,2) - t60 * mrSges(3,1) - t7 * (t51 * mrSges(5,1) + t49 * mrSges(5,3)) + qJD(4) * t42 - t8 * t24 + t2 * mrSges(5,3) - t3 * mrSges(5,1) + t4 * mrSges(4,1) - t5 * mrSges(4,2) - t27 * t136 + t28 * t134 + t17 * t135 + t12 * t137 + t19 * t113 - (-t169 * t49 + t13 - t139 + t46) * t51 / 0.2e1 + (-t12 * t33 + t2 * t75 + t3 * t77 - t7 * t8 + (-t34 + qJD(4)) * t17) * m(5) + (-m(4) * t145 - m(5) * (t97 + t145) + t174) * g(3) + t18 * t147 + t173 * (t102 + t171 * pkin(2) * t86 + (-m(5) * qJ(4) + t175) * t82 + (m(5) * pkin(3) + t176) * t81) + (Ifges(5,2) + Ifges(4,3) + Ifges(3,3)) * qJDD(2) - (-Ifges(3,2) * t123 + t56 + t80) * t122 / 0.2e1 + t14 * t150 + (Ifges(5,3) * t51 - t138) * t153 - t25 * t118 - t98 * t120 / 0.2e1; t163 * t51 + t164 * t49 + t109 + t110 + (-t12 * t51 + t17 * t49 + t1 + t103) * m(5) + (t27 * t51 + t28 * t49 + t103 + t45) * m(4); -qJD(2) * t42 + t51 * t24 + (g(3) * t82 - t17 * qJD(2) - t173 * t81 + t7 * t51 + t3) * m(5) + t20;];
tau = t9;
