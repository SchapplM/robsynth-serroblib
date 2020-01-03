% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:50:59
% DurationCPUTime: 1.71s
% Computational Cost: add. (2133->229), mult. (3850->290), div. (0->0), fcn. (2344->6), ass. (0->110)
t177 = mrSges(6,1) + mrSges(5,1);
t175 = Ifges(5,1) + Ifges(6,1);
t168 = Ifges(5,5) + Ifges(6,4);
t176 = -Ifges(6,5) + Ifges(5,4);
t108 = cos(qJ(2));
t135 = pkin(1) * qJD(1);
t116 = -t108 * t135 + qJD(3);
t103 = qJD(1) + qJD(2);
t105 = cos(pkin(8));
t152 = cos(qJ(4));
t123 = t152 * t105;
t104 = sin(pkin(8));
t106 = sin(qJ(4));
t132 = t106 * t104;
t110 = t123 - t132;
t81 = t110 * qJD(4);
t68 = t103 * t81;
t163 = t68 / 0.2e1;
t167 = Ifges(6,6) - Ifges(5,6);
t73 = t110 * t103;
t149 = Ifges(6,5) * t73;
t71 = Ifges(5,4) * t73;
t88 = t104 * t152 + t106 * t105;
t74 = t88 * t103;
t173 = qJD(4) * t168 + t175 * t74 - t149 + t71;
t113 = -mrSges(4,1) * t105 + mrSges(4,2) * t104;
t172 = mrSges(5,1) * t73 - mrSges(5,2) * t74 - t113 * t103;
t143 = t73 * mrSges(5,3);
t61 = t73 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t138 = -qJD(4) * mrSges(5,2) + t143 + t61;
t142 = t74 * mrSges(5,3);
t137 = -t74 * mrSges(6,2) + qJD(4) * t177 - t142;
t131 = t104 ^ 2 + t105 ^ 2;
t121 = t131 * mrSges(4,3);
t107 = sin(qJ(2));
t129 = t107 * t135;
t90 = qJ(3) * t103 + t129;
t122 = pkin(7) * t103 + t90;
t64 = t122 * t104;
t65 = t122 * t105;
t35 = -t106 * t64 + t152 * t65;
t134 = pkin(1) * qJD(2);
t124 = qJD(1) * t134;
t118 = t108 * t124;
t86 = t103 * qJD(3) + t118;
t7 = qJD(4) * t35 + t86 * t88;
t153 = t7 * t88;
t125 = t152 * t64;
t133 = t106 * t65;
t34 = -t125 - t133;
t171 = -t34 * t81 + t153;
t23 = qJD(4) * qJ(5) + t35;
t139 = -qJD(4) * t125 + t86 * t123;
t5 = -t86 * t132 + (qJD(5) - t133) * qJD(4) + t139;
t82 = t88 * qJD(4);
t170 = t5 * t110 - t23 * t82;
t92 = (-pkin(7) - qJ(3)) * t104;
t99 = t105 * pkin(7);
t93 = qJ(3) * t105 + t99;
t57 = t106 * t92 + t152 * t93;
t166 = -qJD(4) * t57 - t116 * t88;
t111 = -t106 * t93 + t152 * t92;
t165 = qJD(4) * t111 + t116 * t110;
t164 = -t34 + qJD(5);
t162 = t73 / 0.2e1;
t161 = -t73 / 0.2e1;
t159 = t74 / 0.2e1;
t97 = pkin(1) * t107 + qJ(3);
t83 = (-pkin(7) - t97) * t104;
t84 = t105 * t97 + t99;
t112 = -t106 * t84 + t152 * t83;
t155 = t112 * t7;
t154 = t111 * t7;
t150 = Ifges(5,4) * t74;
t148 = pkin(1) * t108;
t53 = t106 * t83 + t152 * t84;
t69 = t103 * t82;
t145 = t53 * t69;
t144 = t68 * mrSges(6,2);
t141 = -qJD(4) / 0.2e1;
t140 = qJD(4) / 0.2e1;
t136 = mrSges(4,3) * t103;
t127 = t107 * t134;
t126 = t108 * t134;
t98 = -t105 * pkin(3) - pkin(2);
t38 = t69 * mrSges(5,1) + t68 * mrSges(5,2);
t37 = t69 * mrSges(6,1) - t68 * mrSges(6,3);
t120 = t131 * t90;
t119 = qJD(3) * t131;
t117 = t107 * t124;
t115 = -t111 * t68 - t57 * t69;
t36 = pkin(4) * t82 - qJ(5) * t81 - qJD(5) * t88;
t54 = -pkin(4) * t110 - t88 * qJ(5) + t98;
t72 = t103 * t98 + t116;
t21 = -t73 * pkin(4) - t74 * qJ(5) + t72;
t22 = -qJD(4) * pkin(4) + t164;
t70 = Ifges(6,5) * t74;
t39 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t73 + t70;
t40 = Ifges(5,2) * t73 + Ifges(5,6) * qJD(4) + t150;
t6 = (-qJD(4) * t65 - t104 * t86) * t106 + t139;
t8 = pkin(4) * t69 - qJ(5) * t68 - qJD(5) * t74 + t117;
t109 = mrSges(6,2) * t153 + t113 * t117 + (t175 * t68 - t176 * t69) * t88 / 0.2e1 + t86 * t121 + (t39 / 0.2e1 + t21 * mrSges(6,1) + t72 * mrSges(5,1) - t40 / 0.2e1 - t35 * mrSges(5,3) + Ifges(6,3) * t161 - Ifges(5,2) * t162 - t176 * t159 + t167 * t140) * t82 + (mrSges(5,2) * t117 - t8 * mrSges(6,3) + (-Ifges(5,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t69 + t175 * t163) * t88 + (-mrSges(5,1) * t117 - t8 * mrSges(6,1) + t6 * mrSges(5,3) + (-Ifges(5,2) - Ifges(6,3)) * t69 + 0.2e1 * t176 * t163) * t110 + (t173 / 0.2e1 + t72 * mrSges(5,2) + t22 * mrSges(6,2) - t21 * mrSges(6,3) + Ifges(5,4) * t162 + Ifges(6,5) * t161 + t168 * t140 + t175 * t159) * t81;
t95 = qJD(3) + t126;
t91 = t98 - t148;
t89 = -t103 * pkin(2) + t116;
t48 = t54 - t148;
t46 = -mrSges(6,1) * t73 - mrSges(6,3) * t74;
t45 = pkin(4) * t74 - qJ(5) * t73;
t29 = t36 + t127;
t1 = [m(5) * (t53 * t6 - t155) + m(6) * (t21 * t29 + t48 * t8 + t5 * t53 - t155) + t91 * t38 + t48 * t37 + t29 * t46 - mrSges(3,1) * t117 + t109 - t112 * t144 + (-m(5) * t34 + m(6) * t22 - t137) * (qJD(4) * t53 + t88 * t95) + (m(5) * t35 + m(6) * t23 + t138) * (qJD(4) * t112 + t110 * t95) + (-t103 * t126 - t118) * mrSges(3,2) + (m(5) * (qJD(1) * t91 + t72) + m(4) * (t89 + (-pkin(2) - t148) * qJD(1)) - t103 * mrSges(3,1) - t172) * t127 + (-t145 + t170) * mrSges(6,2) + (-t112 * t68 - t145 + t171) * mrSges(5,3) + (m(4) * (t86 * t97 + t90 * t95) + t95 * t136) * t131; t119 * t136 + ((-mrSges(3,2) * qJD(2) + (mrSges(3,2) - t121) * t103) * t108 + (-t46 + (-qJD(2) + t103) * mrSges(3,1) + t172) * t107) * t135 + (t115 + t171) * mrSges(5,3) + (t115 + t170) * mrSges(6,2) + t98 * t38 + t54 * t37 + t36 * t46 + t109 + t165 * t138 + t166 * t137 + (t5 * t57 + t54 * t8 - t154 + t165 * t23 - t166 * t22 + (-t129 + t36) * t21) * m(6) + (t117 * t98 - t129 * t72 + t165 * t35 + t166 * t34 + t57 * t6 - t154) * m(5) + (-(t107 * t89 + t108 * t120) * t135 + t131 * t86 * qJ(3) - pkin(2) * t117 + t90 * t119) * m(4); t137 * t74 - t138 * t73 + (m(4) + m(5)) * t117 - m(5) * (-t34 * t74 + t35 * t73) + t37 + t38 + (-m(4) * t120 - t121 * t103) * t103 + (-t22 * t74 - t23 * t73 + t8) * m(6); (Ifges(6,3) * t74 + t149) * t162 - t21 * (t74 * mrSges(6,1) - t73 * mrSges(6,3)) - t72 * (mrSges(5,1) * t74 + mrSges(5,2) * t73) + t40 * t159 + qJD(5) * t61 - t45 * t46 + t5 * mrSges(6,3) - t6 * mrSges(5,2) - t177 * t7 + t167 * t69 + t168 * t68 + (t137 + t142) * t35 + (-t138 + t143) * t34 + (-pkin(4) * t68 - qJ(5) * t69 - t22 * t73 + t23 * t74) * mrSges(6,2) + (t167 * t74 + t168 * t73) * t141 + (-pkin(4) * t7 + qJ(5) * t5 + t164 * t23 - t21 * t45 - t22 * t35) * m(6) + (-Ifges(5,2) * t74 + t173 + t71) * t161 - (t175 * t73 - t150 + t39 + t70) * t74 / 0.2e1; t144 - qJD(4) * t61 + t74 * t46 + 0.2e1 * (t7 / 0.2e1 + t23 * t141 + t21 * t159) * m(6);];
tauc = t1(:);
