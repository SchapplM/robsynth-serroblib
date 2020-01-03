% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP1
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:58:57
% DurationCPUTime: 1.30s
% Computational Cost: add. (1286->220), mult. (2605->292), div. (0->0), fcn. (1292->6), ass. (0->113)
t170 = Ifges(5,1) + Ifges(6,1);
t95 = sin(qJ(4));
t97 = cos(qJ(4));
t169 = -mrSges(5,1) * t97 + mrSges(5,2) * t95 - mrSges(4,1);
t96 = sin(qJ(2));
t98 = cos(qJ(2));
t165 = mrSges(3,1) * t96 + mrSges(3,2) * t98;
t137 = qJD(1) * pkin(1);
t128 = t96 * t137;
t127 = t98 * t137;
t90 = qJD(1) + qJD(2);
t73 = t90 * pkin(2) + t127;
t93 = sin(pkin(8));
t94 = cos(pkin(8));
t29 = t94 * t128 + t93 * t73;
t21 = pkin(7) * t90 + t29;
t119 = qJ(5) * t90 + t21;
t134 = qJD(3) * t95;
t12 = t119 * t97 + t134;
t148 = t90 * t95;
t69 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t148;
t70 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t148;
t147 = t90 * t97;
t71 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t147;
t144 = t97 * mrSges(5,3);
t72 = -qJD(4) * mrSges(5,2) + t90 * t144;
t164 = (t71 + t72) * t97 - (t69 + t70) * t95;
t163 = pkin(2) * t94;
t162 = pkin(4) * t97;
t17 = t21 * t97 + t134;
t138 = pkin(1) * qJD(2);
t146 = t93 * t96;
t61 = (t94 * t98 - t146) * t138;
t41 = qJD(1) * t61;
t6 = -qJD(4) * t17 - t41 * t95;
t161 = t6 * t95;
t158 = Ifges(5,4) * t95;
t156 = Ifges(6,4) * t95;
t88 = t97 * qJD(3);
t16 = -t21 * t95 + t88;
t154 = t16 * mrSges(5,3);
t153 = t16 * t95;
t152 = t17 * mrSges(5,3);
t145 = t94 * t96;
t86 = pkin(1) * t98 + pkin(2);
t140 = pkin(1) * t145 + t93 * t86;
t57 = pkin(7) + t140;
t151 = t57 * t97;
t150 = t61 * t97;
t149 = t90 * mrSges(4,2);
t143 = qJD(4) * t88 + t97 * t41;
t133 = qJD(4) * t95;
t125 = t90 * t133;
t132 = qJD(4) * t97;
t50 = t90 * mrSges(6,2) * t132 + mrSges(6,1) * t125;
t136 = -qJ(5) - t57;
t83 = pkin(2) * t93 + pkin(7);
t135 = -qJ(5) - t83;
t87 = t97 * qJD(5);
t126 = pkin(4) * t133;
t124 = -pkin(3) - t162;
t121 = t169 * t90;
t78 = t93 * t128;
t28 = t73 * t94 - t78;
t120 = -pkin(1) * t146 + t86 * t94;
t118 = qJD(4) * t136;
t117 = qJD(4) * t135;
t56 = -pkin(3) - t120;
t113 = -mrSges(6,1) * t97 + mrSges(6,2) * t95;
t109 = t119 * t95;
t11 = t88 - t109;
t10 = qJD(4) * pkin(4) + t11;
t112 = -t10 * t95 + t12 * t97;
t111 = -t16 * t97 - t17 * t95;
t110 = t17 * t97 - t153;
t108 = m(6) * t112;
t107 = pkin(1) * (t93 * t98 + t145);
t106 = (Ifges(5,2) * t97 + t158) * t90;
t105 = (Ifges(6,2) * t97 + t156) * t90;
t104 = (mrSges(5,1) * t95 + mrSges(5,2) * t97) * qJD(4);
t59 = qJD(2) * t107;
t40 = qJD(1) * t59;
t15 = t124 * t90 + qJD(5) - t28;
t19 = pkin(4) * t125 + t40;
t2 = -qJD(4) * t109 + t90 * t87 + t143;
t20 = -pkin(3) * t90 - t28;
t3 = (-qJD(5) * t90 - t41) * t95 - t12 * qJD(4);
t46 = Ifges(6,6) * qJD(4) + t105;
t47 = Ifges(5,6) * qJD(4) + t106;
t79 = Ifges(6,4) * t147;
t48 = Ifges(6,1) * t148 + Ifges(6,5) * qJD(4) + t79;
t80 = Ifges(5,4) * t147;
t49 = Ifges(5,1) * t148 + Ifges(5,5) * qJD(4) + t80;
t5 = -t21 * t133 + t143;
t99 = t19 * t113 + (-t6 * mrSges(5,3) - t3 * mrSges(6,3)) * t95 + t5 * t144 + t20 * t104 + t2 * t97 * mrSges(6,3) - t41 * mrSges(4,2) + t15 * (mrSges(6,1) * t95 + mrSges(6,2) * t97) * qJD(4) + t169 * t40 + (t170 * t97 - t156 - t158) * t125 + ((Ifges(5,5) + Ifges(6,5)) * t97 + (-Ifges(5,6) - Ifges(6,6)) * t95) * qJD(4) ^ 2 / 0.2e1 - (t105 + t106 + t47 + t46) * t133 / 0.2e1 + (t49 + t48 + ((-0.2e1 * Ifges(5,2) - 0.2e1 * Ifges(6,2) + t170) * t95 + 0.3e1 * (Ifges(5,4) + Ifges(6,4)) * t97) * t90) * t132 / 0.2e1;
t89 = t97 * qJ(5);
t84 = -pkin(3) - t163;
t74 = t124 - t163;
t68 = t83 * t97 + t89;
t67 = t135 * t95;
t62 = t113 * t90;
t60 = t94 * t127 - t78;
t58 = qJD(1) * t107;
t51 = t90 * t104;
t39 = -qJD(5) * t95 + t97 * t117;
t38 = t95 * t117 + t87;
t36 = t56 - t162;
t35 = t59 + t126;
t23 = t89 + t151;
t22 = t136 * t95;
t8 = (-qJD(5) - t61) * t95 + t97 * t118;
t7 = t95 * t118 + t150 + t87;
t1 = [t99 + (-t95 * t70 + t97 * t72 - t149) * t61 + t121 * t59 + m(5) * (t17 * t150 + t5 * t151 - t61 * t153 - t57 * t161 + t20 * t59 + t40 * t56) + m(4) * (-t40 * t120 + t41 * t140 - t28 * t59 + t29 * t61) + m(6) * (t10 * t8 + t12 * t7 + t15 * t35 + t19 * t36 + t2 * t23 + t22 * t3) + (t111 * mrSges(5,3) + (m(5) * t111 - t97 * t70 - t95 * t72) * t57 + (-t10 * t97 - t12 * t95 + (-t22 * t97 - t23 * t95) * t90) * mrSges(6,3)) * qJD(4) + t7 * t71 + t56 * t51 + t35 * t62 + t8 * t69 + t36 * t50 + t165 * t138 * (-qJD(1) - t90); t99 + ((-t154 - t70 * t83 + (-t67 * t90 - t10) * mrSges(6,3)) * t97 + (-t152 + pkin(4) * t62 - t72 * t83 + (-t68 * t90 - t12) * mrSges(6,3)) * t95) * qJD(4) + (t149 - t164) * t60 + t38 * t71 + t74 * t50 + t84 * t51 + t39 * t69 + (-t121 - t62) * t58 + t165 * t137 * (-qJD(2) + t90) + (t10 * t39 - t112 * t60 + t12 * t38 + t19 * t74 + t2 * t68 + t3 * t67 + (t126 - t58) * t15) * m(6) + ((-t40 * t94 + t41 * t93) * pkin(2) + t28 * t58 - t29 * t60) * m(4) + (-t110 * t60 - t20 * t58 + t40 * t84 + (t111 * qJD(4) + t5 * t97 - t161) * t83) * m(5); m(5) * (t5 * t95 + t6 * t97) + m(6) * (t2 * t95 + t3 * t97) + (m(5) * t110 + t108 + (mrSges(5,3) + mrSges(6,3)) * t90 * (-t95 ^ 2 - t97 ^ 2) + t164) * qJD(4); t6 * mrSges(5,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2) - t11 * t71 - t16 * t72 + t17 * t70 + (m(6) * pkin(4) + mrSges(6,1)) * t3 + (-m(6) * (-t10 + t11) + t69) * t12 + ((t10 * mrSges(6,3) + t154 - t79 / 0.2e1 - t80 / 0.2e1 - t49 / 0.2e1 - t48 / 0.2e1 - t20 * mrSges(5,2) - t15 * mrSges(6,2) + (Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,3)) * qJD(4)) * t97 + (t152 + t12 * mrSges(6,3) + t46 / 0.2e1 + t47 / 0.2e1 - t20 * mrSges(5,1) - t15 * mrSges(6,1) + (Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t148 + (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t15 - t62) * pkin(4) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t147) * t95) * t90; m(6) * t19 + (t95 * t69 - t97 * t71 - t108) * t90 + t50;];
tauc = t1(:);
