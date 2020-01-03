% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:29
% DurationCPUTime: 1.09s
% Computational Cost: add. (1258->207), mult. (2591->285), div. (0->0), fcn. (1268->6), ass. (0->102)
t160 = mrSges(5,1) + mrSges(6,1);
t85 = sin(qJ(4));
t87 = cos(qJ(4));
t159 = -mrSges(5,1) * t87 + mrSges(5,2) * t85 - mrSges(4,1);
t128 = pkin(1) * qJD(1);
t86 = sin(qJ(2));
t116 = t86 * t128;
t88 = cos(qJ(2));
t115 = t88 * t128;
t80 = qJD(1) + qJD(2);
t68 = t80 * pkin(2) + t115;
t83 = sin(pkin(8));
t84 = cos(pkin(8));
t27 = t84 * t116 + t83 * t68;
t20 = pkin(7) * t80 + t27;
t139 = t20 * t85;
t15 = qJD(3) * t87 - t139;
t155 = -t15 + qJD(5);
t12 = -qJD(4) * pkin(4) + t155;
t124 = qJD(4) * t87;
t16 = qJD(3) * t85 + t20 * t87;
t127 = pkin(1) * qJD(2);
t136 = t83 * t86;
t55 = (t84 * t88 - t136) * t127;
t35 = qJD(1) * t55;
t6 = qJD(4) * t16 + t35 * t85;
t149 = t6 * t85;
t133 = qJD(3) * t124 + t87 * t35;
t4 = (qJD(5) - t139) * qJD(4) + t133;
t158 = t12 * t124 + t4 * t87 + t149;
t156 = mrSges(3,1) * t86 + mrSges(3,2) * t88;
t137 = t80 * t87;
t117 = mrSges(6,2) * t137;
t67 = qJD(4) * mrSges(6,3) + t117;
t131 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t137 + t67;
t134 = t85 * mrSges(5,3);
t138 = t80 * t85;
t132 = -mrSges(6,2) * t138 + t160 * qJD(4) - t80 * t134;
t154 = t131 * t87 - t132 * t85;
t153 = pkin(2) * t84;
t125 = qJD(4) * t85;
t5 = -t20 * t125 + t133;
t151 = t5 * t87;
t135 = t84 * t86;
t78 = pkin(1) * t88 + pkin(2);
t130 = pkin(1) * t135 + t83 * t78;
t51 = pkin(7) + t130;
t150 = t51 * t6;
t148 = t6 * t87;
t145 = Ifges(5,4) * t85;
t143 = Ifges(6,5) * t85;
t142 = Ifges(6,5) * t87;
t119 = qJD(4) * qJ(5);
t14 = t16 + t119;
t141 = t14 * mrSges(6,2);
t140 = t15 * mrSges(5,3);
t126 = qJD(4) * t51;
t123 = t14 * qJD(4);
t122 = t85 * qJD(5);
t113 = t138 / 0.2e1;
t111 = t125 / 0.2e1;
t109 = t159 * t80;
t70 = t83 * t116;
t26 = t68 * t84 - t70;
t108 = -pkin(1) * t136 + t78 * t84;
t107 = t80 * t111;
t104 = -mrSges(6,1) * t87 - mrSges(6,3) * t85;
t103 = -t16 * mrSges(5,3) - t141;
t102 = pkin(4) * t85 - qJ(5) * t87;
t101 = t12 * t85 + t14 * t87;
t100 = -t15 * t85 + t16 * t87;
t99 = -pkin(4) * t87 - qJ(5) * t85 - pkin(3);
t98 = pkin(1) * (t83 * t88 + t135);
t97 = (Ifges(6,1) * t85 - t142) * t80;
t96 = (Ifges(5,2) * t87 + t145) * t80;
t95 = (mrSges(5,1) * t85 + mrSges(5,2) * t87) * qJD(4);
t94 = (mrSges(6,1) * t85 - mrSges(6,3) * t87) * qJD(4);
t53 = qJD(2) * t98;
t56 = pkin(4) * t125 - t87 * t119 - t122;
t34 = qJD(1) * t53;
t10 = t34 + (t102 * qJD(4) - t122) * t80;
t11 = t99 * t80 - t26;
t19 = -pkin(3) * t80 - t26;
t71 = Ifges(6,5) * t138;
t40 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t137 + t71;
t41 = Ifges(5,6) * qJD(4) + t96;
t42 = Ifges(6,4) * qJD(4) + t97;
t72 = Ifges(5,4) * t137;
t43 = Ifges(5,1) * t138 + Ifges(5,5) * qJD(4) + t72;
t89 = t10 * t104 - t80 * (Ifges(6,3) * t85 + t142) * t124 + mrSges(5,3) * t151 + t11 * t94 + t19 * t95 + (-Ifges(6,3) * t87 + t143) * t107 + t40 * t111 + t159 * t34 + 0.2e1 * (t143 - t145 + (Ifges(5,1) + Ifges(6,1)) * t87) * t107 + ((Ifges(6,4) + Ifges(5,5)) * t87 + (-Ifges(5,6) + Ifges(6,6)) * t85) * qJD(4) ^ 2 / 0.2e1 - (t96 + t41) * t125 / 0.2e1 + t158 * mrSges(6,2) + ((0.3e1 * Ifges(5,4) * t87 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t85) * t80 + t97 + t43 + t42) * t124 / 0.2e1;
t76 = -pkin(3) - t153;
t63 = t99 - t153;
t59 = t102 * t80;
t57 = t104 * t80;
t54 = t84 * t115 - t70;
t52 = qJD(1) * t98;
t50 = -pkin(3) - t108;
t45 = t80 * t95;
t44 = t80 * t94;
t21 = -t108 + t99;
t17 = t53 + t56;
t1 = [(t6 * mrSges(5,3) - t132 * t55 + (-t131 * t51 + t103) * qJD(4) + m(6) * (t12 * t55 - t51 * t123 + t150) + m(5) * (-t16 * t126 - t15 * t55 + t150)) * t85 + (t131 * t55 + (-t132 * t51 - t140) * qJD(4) + m(6) * (t12 * t126 + t14 * t55 + t4 * t51) + m(5) * (-t15 * t126 + t16 * t55 + t5 * t51)) * t87 + t109 * t53 + m(6) * (t10 * t21 + t11 * t17) + m(5) * (t19 * t53 + t34 * t50) + t89 + m(4) * (-t34 * t108 + t35 * t130 - t26 * t53 + t27 * t55) + (-t55 * t80 - t35) * mrSges(4,2) + t50 * t45 + t17 * t57 + t21 * t44 + t156 * t127 * (-qJD(1) - t80); t89 + (-t85 * t141 + (-t15 * t87 - t16 * t85) * mrSges(5,3)) * qJD(4) + (t80 * mrSges(4,2) - t154) * t54 + t76 * t45 + t56 * t57 + t63 * t44 - t35 * mrSges(4,2) + t6 * t134 + (-t109 - t57) * t52 + t156 * t128 * (-qJD(2) + t80) + (t10 * t63 - t101 * t54 + (-t52 + t56) * t11) * m(6) + (-t100 * t54 - t19 * t52 + t34 * t76) * m(5) + ((-t34 * t84 + t35 * t83) * pkin(2) + t26 * t52 - t27 * t54) * m(4) + (m(6) * (-t85 * t123 + t158) + m(5) * (-t15 * t124 - t16 * t125 + t149 + t151) + (-t131 * t85 - t132 * t87) * qJD(4)) * (pkin(2) * t83 + pkin(7)); m(5) * (t5 * t85 - t148) + m(6) * (t4 * t85 - t148) + (m(5) * t100 + m(6) * t101 + (mrSges(6,2) + mrSges(5,3)) * t80 * (-t85 ^ 2 - t87 ^ 2) + t154) * qJD(4); -t5 * mrSges(5,2) + t4 * mrSges(6,3) + qJD(5) * t67 - t59 * t57 - t160 * t6 + t132 * t16 - t131 * t15 + ((t11 * mrSges(6,3) - t19 * mrSges(5,2) - t72 / 0.2e1 - t42 / 0.2e1 - t43 / 0.2e1 - t12 * mrSges(6,2) + t140 + Ifges(6,5) * t137 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,2)) * qJD(4)) * t87 + (-t11 * mrSges(6,1) - t19 * mrSges(5,1) - t71 / 0.2e1 - t40 / 0.2e1 + t41 / 0.2e1 + Ifges(5,4) * t113 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1 - qJ(5) * mrSges(6,2)) * qJD(4) + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t137 - t103) * t85) * t80 + (-t6 * pkin(4) + t4 * qJ(5) - t11 * t59 - t12 * t16 + t14 * t155) * m(6); t57 * t138 + (-t67 + t117) * qJD(4) + 0.2e1 * (t6 / 0.2e1 + t11 * t113 - t123 / 0.2e1) * m(6);];
tauc = t1(:);
