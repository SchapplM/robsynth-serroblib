% Calculate vector of inverse dynamics joint torques for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:39
% EndTime: 2019-12-05 15:42:47
% DurationCPUTime: 3.66s
% Computational Cost: add. (2653->298), mult. (6112->405), div. (0->0), fcn. (4504->12), ass. (0->126)
t94 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t112 = qJDD(4) + qJDD(5);
t118 = qJD(4) + qJD(5);
t122 = sin(qJ(5));
t124 = cos(qJ(5));
t119 = sin(pkin(9));
t123 = sin(qJ(4));
t120 = cos(pkin(9));
t125 = cos(qJ(4));
t153 = t120 * t125;
t83 = -t119 * t123 + t153;
t75 = t83 * qJD(2);
t84 = t119 * t125 + t120 * t123;
t76 = t84 * qJD(2);
t137 = -t122 * t76 + t124 * t75;
t50 = t122 * t75 + t124 * t76;
t164 = Ifges(6,4) * t50;
t77 = t83 * qJD(4);
t52 = qJD(2) * t77 + qJDD(2) * t84;
t78 = t84 * qJD(4);
t53 = -qJD(2) * t78 + qJDD(2) * t83;
t17 = qJD(5) * t137 + t122 * t53 + t124 * t52;
t18 = -qJD(5) * t50 - t122 * t52 + t124 * t53;
t106 = t120 * qJD(1);
t160 = pkin(6) + qJ(3);
t91 = t160 * t119;
t70 = -qJD(2) * t91 + t106;
t152 = qJD(2) * t120;
t87 = qJ(3) * t152 + t119 * qJD(1);
t71 = pkin(6) * t152 + t87;
t39 = t123 * t70 + t125 * t71;
t104 = t120 * qJDD(1);
t62 = t104 + (-pkin(6) * qJDD(2) - t94) * t119;
t146 = qJDD(2) * t120;
t74 = t119 * qJDD(1) + t120 * t94;
t63 = pkin(6) * t146 + t74;
t20 = -qJD(4) * t39 - t123 * t63 + t125 * t62;
t6 = qJDD(4) * pkin(4) - pkin(7) * t52 + t20;
t151 = qJD(4) * t125;
t155 = t123 * t71;
t19 = -qJD(4) * t155 + t123 * t62 + t125 * t63 + t70 * t151;
t7 = pkin(7) * t53 + t19;
t32 = pkin(7) * t75 + t39;
t157 = t122 * t32;
t38 = t125 * t70 - t155;
t31 = -pkin(7) * t76 + t38;
t28 = qJD(4) * pkin(4) + t31;
t8 = t124 * t28 - t157;
t2 = qJD(5) * t8 + t122 * t6 + t124 * t7;
t42 = Ifges(6,4) * t137;
t24 = Ifges(6,1) * t50 + Ifges(6,5) * t118 + t42;
t154 = t124 * t32;
t9 = t122 * t28 + t154;
t3 = -qJD(5) * t9 - t122 * t7 + t124 * t6;
t100 = pkin(3) * t120 + pkin(2);
t90 = -qJD(2) * t100 + qJD(3);
t56 = -pkin(4) * t75 + t90;
t187 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t112 - (Ifges(6,5) * t137 - Ifges(6,6) * t50) * t118 / 0.2e1 + (t137 * t8 + t50 * t9) * mrSges(6,3) - (-Ifges(6,2) * t50 + t24 + t42) * t137 / 0.2e1 - t56 * (mrSges(6,1) * t50 + mrSges(6,2) * t137) - (Ifges(6,1) * t137 - t164) * t50 / 0.2e1;
t186 = t119 ^ 2 + t120 ^ 2;
t116 = pkin(9) + qJ(4);
t107 = sin(t116);
t109 = cos(t116);
t111 = qJ(5) + t116;
t98 = sin(t111);
t99 = cos(t111);
t138 = t99 * mrSges(6,1) - t98 * mrSges(6,2);
t185 = -t109 * mrSges(5,1) + t107 * mrSges(5,2) - t138;
t117 = pkin(8) + qJ(2);
t108 = sin(t117);
t110 = cos(t117);
t176 = g(1) * t110 + g(2) * t108;
t23 = Ifges(6,2) * t137 + Ifges(6,6) * t118 + t164;
t183 = t23 / 0.2e1;
t73 = -t119 * t94 + t104;
t177 = -t119 * t73 + t120 * t74;
t134 = -mrSges(4,1) * t120 + mrSges(4,2) * t119;
t175 = m(4) * pkin(2) + m(5) * t100 + mrSges(3,1) + m(6) * (pkin(4) * t109 + t100) - t134 - t185;
t174 = mrSges(3,2) + m(6) * (-pkin(7) - t160) - mrSges(6,3) - m(5) * t160 - mrSges(5,3) - m(4) * qJ(3) - mrSges(4,3);
t171 = t50 / 0.2e1;
t169 = t76 / 0.2e1;
t168 = m(2) + m(3);
t167 = pkin(4) * t78;
t165 = Ifges(5,4) * t76;
t92 = t160 * t120;
t58 = -t123 * t91 + t125 * t92;
t150 = m(4) + m(5) + m(6);
t147 = qJDD(2) * t119;
t140 = -t53 * mrSges(5,1) + t52 * mrSges(5,2);
t139 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t57 = -t123 * t92 - t125 * t91;
t136 = -mrSges(4,1) * t146 + mrSges(4,2) * t147;
t135 = mrSges(6,1) * t98 + mrSges(6,2) * t99;
t132 = -t119 * (-qJ(3) * qJD(2) * t119 + t106) + t120 * t87;
t40 = -pkin(7) * t84 + t57;
t41 = pkin(7) * t83 + t58;
t21 = -t122 * t41 + t124 * t40;
t22 = t122 * t40 + t124 * t41;
t54 = -t122 * t84 + t124 * t83;
t55 = t122 * t83 + t124 * t84;
t88 = -qJDD(2) * t100 + qJDD(3);
t36 = -t91 * t151 + qJD(3) * t153 + (-qJD(3) * t119 - qJD(4) * t92) * t123;
t37 = -qJD(3) * t84 - qJD(4) * t58;
t102 = -qJDD(2) * pkin(2) + qJDD(3);
t72 = Ifges(5,4) * t75;
t68 = -pkin(4) * t83 - t100;
t66 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t76;
t65 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t75;
t46 = t76 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t72;
t45 = t75 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t165;
t44 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t53;
t43 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t52;
t35 = mrSges(6,1) * t118 - mrSges(6,3) * t50;
t34 = -mrSges(6,2) * t118 + mrSges(6,3) * t137;
t33 = -pkin(4) * t53 + t88;
t30 = -pkin(7) * t77 + t37;
t29 = -pkin(7) * t78 + t36;
t27 = -qJD(5) * t55 - t122 * t77 - t124 * t78;
t26 = qJD(5) * t54 - t122 * t78 + t124 * t77;
t25 = -mrSges(6,1) * t137 + mrSges(6,2) * t50;
t13 = -mrSges(6,2) * t112 + mrSges(6,3) * t18;
t12 = mrSges(6,1) * t112 - mrSges(6,3) * t17;
t11 = t124 * t31 - t157;
t10 = -t122 * t31 - t154;
t5 = -qJD(5) * t22 - t122 * t29 + t124 * t30;
t4 = qJD(5) * t21 + t122 * t30 + t124 * t29;
t1 = [t54 * t12 + t55 * t13 + t26 * t34 + t27 * t35 + t83 * t43 + t84 * t44 + t77 * t65 - t78 * t66 + t168 * qJDD(1) + m(4) * (t119 * t74 + t120 * t73) + m(5) * (t19 * t84 + t20 * t83 - t38 * t78 + t39 * t77) + m(6) * (t2 * t55 + t26 * t9 + t27 * t8 + t3 * t54) + (-t150 - t168) * g(3); t137 * (Ifges(6,4) * t26 + Ifges(6,2) * t27) / 0.2e1 + t90 * (mrSges(5,1) * t78 + mrSges(5,2) * t77) + qJD(4) * (Ifges(5,5) * t77 - Ifges(5,6) * t78) / 0.2e1 + t75 * (Ifges(5,4) * t77 - Ifges(5,2) * t78) / 0.2e1 + (-t38 * t77 - t39 * t78) * mrSges(5,3) + (Ifges(5,1) * t77 - Ifges(5,4) * t78) * t169 + m(6) * (t167 * t56 + t2 * t22 + t21 * t3 + t33 * t68 + t4 * t9 + t5 * t8) + (t186 * t94 + t177) * mrSges(4,3) + t118 * (Ifges(6,5) * t26 + Ifges(6,6) * t27) / 0.2e1 - t78 * t45 / 0.2e1 + t77 * t46 / 0.2e1 + t36 * t65 + t37 * t66 + t56 * (-mrSges(6,1) * t27 + mrSges(6,2) * t26) + t57 * t43 + t58 * t44 + t4 * t34 + t5 * t35 + t21 * t12 + t22 * t13 + t26 * t24 / 0.2e1 + t68 * t139 - t100 * t140 - pkin(2) * t136 + t102 * t134 + (-t26 * t8 + t27 * t9) * mrSges(6,3) + m(5) * (-t100 * t88 + t19 * t58 + t20 * t57 + t36 * t39 + t37 * t38) + (-mrSges(6,1) * t33 + mrSges(6,3) * t2 + Ifges(6,4) * t17 + Ifges(6,2) * t18 + Ifges(6,6) * t112) * t54 + m(4) * (-pkin(2) * t102 + qJ(3) * t177 + t132 * qJD(3)) + (t108 * t175 + t110 * t174) * g(1) + (t108 * t174 - t110 * t175) * g(2) + t27 * t183 + (Ifges(6,1) * t26 + Ifges(6,4) * t27) * t171 + (Ifges(4,4) * t119 + Ifges(4,2) * t120) * t146 + (Ifges(4,1) * t119 + Ifges(4,4) * t120) * t147 + (mrSges(5,2) * t88 - mrSges(5,3) * t20 + Ifges(5,1) * t52 + Ifges(5,4) * t53 + Ifges(5,5) * qJDD(4)) * t84 + (-mrSges(5,1) * t88 + mrSges(5,3) * t19 + Ifges(5,4) * t52 + Ifges(5,2) * t53 + Ifges(5,6) * qJDD(4)) * t83 + (mrSges(6,2) * t33 - mrSges(6,3) * t3 + Ifges(6,1) * t17 + Ifges(6,4) * t18 + Ifges(6,5) * t112) * t55 + t25 * t167 + Ifges(3,3) * qJDD(2); -t137 * t34 + t50 * t35 - t75 * t65 + t76 * t66 - t186 * qJD(2) ^ 2 * mrSges(4,3) + t136 + t139 + t140 + (-g(1) * t108 + g(2) * t110) * t150 + (-t137 * t9 + t50 * t8 + t33) * m(6) + (t38 * t76 - t39 * t75 + t88) * m(5) + (-qJD(2) * t132 + t102) * m(4); -m(6) * (t10 * t8 + t11 * t9) - (-Ifges(5,2) * t76 + t46 + t72) * t75 / 0.2e1 + t50 * t183 - t90 * (mrSges(5,1) * t76 + mrSges(5,2) * t75) - qJD(4) * (Ifges(5,5) * t75 - Ifges(5,6) * t76) / 0.2e1 - t38 * t65 + t39 * t66 + Ifges(5,6) * t53 + Ifges(5,5) * t52 + t176 * (mrSges(5,1) * t107 + mrSges(5,2) * t109 + t135) + t185 * g(3) - t11 * t34 - t10 * t35 - t19 * mrSges(5,2) + t20 * mrSges(5,1) + t187 + (t124 * t12 + t122 * t13 - t76 * t25 + (-g(3) * t109 + t176 * t107 + t122 * t2 + t124 * t3 - t56 * t76) * m(6) + (-t122 * t35 + t124 * t34 + (-t122 * t8 + t124 * t9) * m(6)) * qJD(5)) * pkin(4) + t45 * t169 + Ifges(5,3) * qJDD(4) - t76 * (Ifges(5,1) * t75 - t165) / 0.2e1 + (t38 * t75 + t39 * t76) * mrSges(5,3); -g(3) * t138 + t176 * t135 + t23 * t171 - t8 * t34 + t9 * t35 + t187;];
tau = t1;
