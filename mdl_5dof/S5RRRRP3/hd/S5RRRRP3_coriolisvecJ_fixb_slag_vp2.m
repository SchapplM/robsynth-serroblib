% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:10
% DurationCPUTime: 1.32s
% Computational Cost: add. (2030->216), mult. (3632->295), div. (0->0), fcn. (1715->6), ass. (0->101)
t162 = mrSges(5,1) + mrSges(6,1);
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t161 = -mrSges(5,1) * t85 + mrSges(5,2) * t82 - mrSges(4,1);
t83 = sin(qJ(3));
t129 = qJD(3) * t83;
t130 = pkin(1) * qJD(1);
t87 = cos(qJ(2));
t117 = t87 * t130;
t79 = qJD(1) + qJD(2);
t66 = t79 * pkin(2) + t117;
t84 = sin(qJ(2));
t86 = cos(qJ(3));
t134 = t84 * t86;
t100 = t83 * t87 + t134;
t128 = qJD(3) * t86;
t90 = (qJD(2) * t100 + t128 * t84) * pkin(1);
t16 = qJD(1) * t90 + t129 * t66;
t78 = qJD(3) + t79;
t97 = (mrSges(5,1) * t82 + mrSges(5,2) * t85) * qJD(4);
t42 = t78 * t97;
t160 = m(5) * t16 + t42;
t126 = qJD(4) * t85;
t135 = t83 * t84;
t91 = (-t84 * t129 + (t86 * t87 - t135) * qJD(2)) * pkin(1);
t15 = qJD(1) * t91 + t128 * t66;
t118 = t84 * t130;
t35 = t118 * t86 + t66 * t83;
t24 = pkin(8) * t78 + t35;
t7 = t126 * t24 + t15 * t82;
t150 = t7 * t82;
t139 = t24 * t82;
t19 = -qJD(4) * pkin(4) + qJD(5) + t139;
t12 = t85 * t15;
t4 = t12 + (qJD(5) - t139) * qJD(4);
t159 = t19 * t126 + t4 * t85 + t150;
t121 = qJD(4) * qJ(5);
t20 = t24 * t85 + t121;
t125 = t20 * qJD(4);
t127 = qJD(4) * t82;
t6 = -t127 * t24 + t12;
t152 = t6 * t85;
t158 = m(6) * (-t125 * t82 + t159) + m(5) * (t150 + t152);
t137 = t78 * t82;
t133 = (-mrSges(6,2) - mrSges(5,3)) * t137 + t162 * qJD(4);
t136 = t78 * t85;
t119 = mrSges(6,2) * t136;
t65 = qJD(4) * mrSges(6,3) + t119;
t132 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t136 + t65;
t101 = t19 * t82 + t20 * t85;
t156 = mrSges(3,1) * t84 + mrSges(3,2) * t87;
t114 = (t82 ^ 2 + t85 ^ 2) * t24;
t154 = pkin(2) * t86;
t151 = t7 * mrSges(5,3);
t147 = Ifges(5,4) * t82;
t145 = Ifges(6,5) * t82;
t144 = Ifges(6,5) * t85;
t142 = t20 * mrSges(6,2);
t76 = pkin(1) * t87 + pkin(2);
t131 = pkin(1) * t134 + t83 * t76;
t124 = t82 * qJD(5);
t115 = t137 / 0.2e1;
t112 = t127 / 0.2e1;
t110 = t161 * t78;
t71 = t83 * t118;
t34 = t66 * t86 - t71;
t109 = -pkin(1) * t135 + t76 * t86;
t108 = t133 * qJD(4);
t107 = t78 * t112;
t103 = -mrSges(6,1) * t85 - mrSges(6,3) * t82;
t51 = t103 * t78;
t105 = -t110 - t51;
t102 = pkin(4) * t82 - qJ(5) * t85;
t67 = -pkin(4) * t85 - qJ(5) * t82 - pkin(3);
t99 = (Ifges(6,1) * t82 - t144) * t78;
t98 = (Ifges(5,2) * t85 + t147) * t78;
t96 = (mrSges(6,1) * t82 - mrSges(6,3) * t85) * qJD(4);
t58 = pkin(4) * t127 - t121 * t85 - t124;
t22 = t129 * t76 + t90;
t13 = t67 * t78 - t34;
t23 = -pkin(3) * t78 - t34;
t68 = Ifges(6,5) * t137;
t37 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t136 + t68;
t38 = Ifges(5,6) * qJD(4) + t98;
t39 = Ifges(6,4) * qJD(4) + t99;
t69 = Ifges(5,4) * t136;
t40 = Ifges(5,1) * t137 + Ifges(5,5) * qJD(4) + t69;
t8 = (qJD(4) * t102 - t124) * t78 + t16;
t88 = t8 * t103 - t78 * (Ifges(6,3) * t82 + t144) * t126 + (-Ifges(6,3) * t85 + t145) * t107 + t13 * t96 + t23 * t97 + t37 * t112 + mrSges(5,3) * t152 + t161 * t16 + 0.2e1 * (t145 - t147 + (Ifges(5,1) + Ifges(6,1)) * t85) * t107 + ((Ifges(6,4) + Ifges(5,5)) * t85 + (-Ifges(5,6) + Ifges(6,6)) * t82) * qJD(4) ^ 2 / 0.2e1 - (t98 + t38) * t127 / 0.2e1 + t159 * mrSges(6,2) + ((0.3e1 * Ifges(5,4) * t85 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t82) * t78 + t99 + t40 + t39) * t126 / 0.2e1;
t74 = pkin(2) * t83 + pkin(8);
t61 = t67 - t154;
t57 = t117 * t86 - t71;
t56 = t100 * t130;
t55 = pkin(8) + t131;
t54 = -pkin(3) - t109;
t41 = t78 * t96;
t36 = pkin(2) * t129 + t58;
t33 = -t109 + t67;
t21 = t128 * t76 + t91;
t11 = t22 + t58;
t1 = [t110 * t22 + t88 + m(4) * (-t109 * t16 + t131 * t15 + t35 * t21 - t34 * t22) + (-t21 * t78 - t15) * mrSges(4,2) + t33 * t41 + t11 * t51 + t54 * t42 + (t151 - t133 * t21 + (-t132 * t55 - t142) * qJD(4)) * t82 + (-t108 * t55 + t132 * t21) * t85 + m(5) * (t21 * t114 + t16 * t54 + t22 * t23) + m(6) * (t101 * t21 + t11 * t13 + t33 * t8) + t55 * t158 + t156 * qJD(2) * pkin(1) * (-qJD(1) - t79); t88 - m(6) * (t101 * t57 + t13 * t56) - m(5) * (t57 * t114 + t23 * t56) + (t57 * t78 - t15) * mrSges(4,2) - m(4) * (-t34 * t56 + t35 * t57) + m(6) * (t13 * t36 + t61 * t8) + t74 * t158 + (t151 + t133 * t57 + (-t132 * t74 - t142) * qJD(4)) * t82 + (-t108 * t74 - t132 * t57) * t85 + t105 * t56 + t156 * t130 * (-qJD(2) + t79) + t61 * t41 + t36 * t51 + (m(5) * t114 * t128 + (m(4) * t15 + (-m(4) * t34 + m(5) * t23 + t110) * qJD(3)) * t83 + (-m(4) * t16 + (m(4) * t35 + m(6) * t101 - t78 * mrSges(4,2) + t132 * t85 - t133 * t82) * qJD(3)) * t86) * pkin(2) + t160 * (-pkin(3) - t154); t88 - m(5) * (t114 * t34 + t23 * t35) + pkin(8) * t158 + (t151 + t133 * t34 + (-pkin(8) * t132 - t142) * qJD(4)) * t82 + (-pkin(8) * t108 - t132 * t34) * t85 + t105 * t35 + (t34 * t78 - t15) * mrSges(4,2) + t67 * t41 + t58 * t51 - t160 * pkin(3) + (-t101 * t34 + t8 * t67 + (-t35 + t58) * t13) * m(6); -t6 * mrSges(5,2) + t4 * mrSges(6,3) + qJD(5) * t65 - t162 * t7 + m(6) * (-t7 * pkin(4) + t4 * qJ(5) + t20 * qJD(5)) + ((-m(6) * t19 + t133) * t85 + (m(6) * t20 + t132) * t82) * t24 + ((-t23 * mrSges(5,2) + t13 * mrSges(6,3) - t69 / 0.2e1 - t39 / 0.2e1 - t40 / 0.2e1 - t19 * mrSges(6,2) + Ifges(6,5) * t136 / 0.2e1) * t85 + (-t23 * mrSges(5,1) - t13 * mrSges(6,1) - t68 / 0.2e1 - t37 / 0.2e1 + t38 / 0.2e1 + t142 + Ifges(5,4) * t115 + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t136) * t82 + ((Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,2)) * t85 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1 - qJ(5) * mrSges(6,2)) * t82) * qJD(4) + (-m(6) * t13 - t51) * t102) * t78; t51 * t137 + (-t65 + t119) * qJD(4) + 0.2e1 * (t7 / 0.2e1 + t13 * t115 - t125 / 0.2e1) * m(6);];
tauc = t1(:);
