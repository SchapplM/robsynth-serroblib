% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:20
% EndTime: 2020-01-03 11:22:26
% DurationCPUTime: 1.63s
% Computational Cost: add. (1825->232), mult. (5005->368), div. (0->0), fcn. (3679->8), ass. (0->123)
t100 = cos(pkin(9));
t102 = cos(pkin(7));
t101 = cos(pkin(8));
t99 = sin(pkin(7));
t137 = t101 * t99;
t97 = sin(pkin(9));
t72 = t100 * t102 + t97 * t137;
t61 = t72 * qJD(1);
t152 = qJD(5) + t61;
t103 = sin(qJ(5));
t104 = cos(qJ(5));
t98 = sin(pkin(8));
t131 = qJD(1) * t98;
t120 = t104 * t131;
t128 = t101 * t102;
t64 = (t100 * t128 + t97 * t99) * qJD(1);
t134 = t104 * t98;
t77 = t100 * t134 - t101 * t103;
t151 = -t77 * qJD(5) - t102 * t120 + t103 * t64;
t135 = t103 * t98;
t108 = t100 * t135 + t101 * t104;
t121 = t103 * t131;
t150 = -t108 * qJD(5) - t102 * t121 - t104 * t64;
t130 = qJD(1) * t99;
t119 = t100 * t130;
t126 = qJD(1) * t102;
t63 = t101 * t119 - t126 * t97;
t38 = t104 * t63 + t121 * t99;
t148 = t38 / 0.2e1;
t147 = t72 / 0.2e1;
t146 = Ifges(6,4) * t38;
t85 = qJ(2) * t130 + qJD(3);
t145 = t85 * t99;
t144 = t97 * t98;
t143 = t98 * t99;
t123 = t98 * t130;
t36 = -t103 * t63 + t120 * t99;
t142 = -mrSges(5,1) * t123 - mrSges(6,1) * t36 + mrSges(6,2) * t38 + mrSges(5,3) * t63;
t141 = mrSges(5,1) * t61 + mrSges(5,2) * t63 - (-mrSges(4,1) * t102 - mrSges(4,3) * t137) * qJD(1);
t118 = t101 * t126;
t81 = -pkin(2) * t102 - qJ(3) * t99 - pkin(1);
t71 = t81 * qJD(1) + qJD(2);
t46 = qJ(2) * t118 + t98 * t71;
t32 = -qJ(4) * t126 + t46;
t112 = pkin(3) * t98 - qJ(4) * t101;
t51 = t112 * t130 + t85;
t17 = t100 * t32 + t97 * t51;
t139 = qJ(2) * t128 + t98 * t81;
t50 = -qJ(4) * t102 + t139;
t58 = (qJ(2) + t112) * t99;
t140 = t100 * t50 + t97 * t58;
t95 = t99 ^ 2;
t96 = t102 ^ 2;
t138 = t95 + t96;
t136 = t102 * t98;
t109 = (-qJD(4) * t101 + qJD(2)) * t99;
t106 = qJD(1) * t109;
t129 = qJD(3) * t99;
t122 = t98 * t129;
t107 = -qJD(4) * t102 - t122;
t125 = qJD(2) * t102;
t88 = t101 * t125;
t82 = qJD(1) * t88;
t52 = t107 * qJD(1) + t82;
t23 = -t100 * t106 + t52 * t97;
t133 = t23 * t100;
t74 = t101 * t129 + t125 * t98;
t65 = t74 * qJD(1);
t132 = t65 * t101;
t124 = qJD(1) * qJD(2);
t89 = qJ(2) * t136;
t117 = qJ(2) * t124;
t45 = -qJD(1) * t89 + t101 * t71;
t116 = t101 * t81 - t89;
t115 = t102 * pkin(3) - t116;
t24 = t100 * t52 + t106 * t97;
t31 = pkin(3) * t126 + qJD(4) - t45;
t11 = pkin(4) * t61 - pkin(6) * t63 + t31;
t13 = pkin(6) * t123 + t17;
t5 = -t103 * t13 + t104 * t11;
t1 = qJD(5) * t5 + t103 * t65 + t104 * t24;
t6 = t103 * t11 + t104 * t13;
t2 = -qJD(5) * t6 - t103 * t24 + t104 * t65;
t114 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t16 = t100 * t51 - t32 * t97;
t111 = t100 * t58 - t50 * t97;
t66 = -qJD(1) * t122 + t82;
t110 = -t66 * t98 + t132;
t73 = t100 * t137 - t102 * t97;
t18 = pkin(4) * t72 - pkin(6) * t73 + t115;
t20 = pkin(6) * t143 + t140;
t7 = -t103 * t20 + t104 * t18;
t8 = t103 * t18 + t104 * t20;
t47 = -t103 * t73 + t99 * t134;
t48 = t104 * t73 + t99 * t135;
t69 = (mrSges(4,1) * t98 + mrSges(4,2) * t101) * t130;
t91 = t95 * t117;
t78 = (mrSges(4,2) * t102 - mrSges(4,3) * t143) * qJD(1);
t75 = t88 - t122;
t62 = t118 * t97 - t119;
t56 = t107 + t88;
t54 = t77 * t130;
t53 = t108 * t130;
t43 = -mrSges(5,2) * t123 - mrSges(5,3) * t61;
t41 = t48 * qJD(5);
t40 = t47 * qJD(5);
t35 = Ifges(6,4) * t36;
t34 = t38 * qJD(5);
t33 = t36 * qJD(5);
t30 = Ifges(6,5) * t33;
t29 = Ifges(6,6) * t34;
t26 = t100 * t56 + t109 * t97;
t25 = -t100 * t109 + t56 * t97;
t22 = mrSges(6,1) * t152 - mrSges(6,3) * t38;
t21 = -mrSges(6,2) * t152 + mrSges(6,3) * t36;
t19 = -pkin(4) * t143 - t111;
t14 = mrSges(6,1) * t34 + mrSges(6,2) * t33;
t12 = -pkin(4) * t123 - t16;
t10 = Ifges(6,1) * t38 + Ifges(6,5) * t152 + t35;
t9 = Ifges(6,2) * t36 + Ifges(6,6) * t152 + t146;
t4 = -qJD(5) * t8 - t103 * t26 + t104 * t74;
t3 = qJD(5) * t7 + t103 * t74 + t104 * t26;
t15 = [0.2e1 * t138 * mrSges(3,3) * t124 + t152 * (Ifges(6,5) * t40 - Ifges(6,6) * t41) / 0.2e1 + (t65 * mrSges(4,1) + t66 * mrSges(4,2)) * t102 + t75 * t78 + t40 * t10 / 0.2e1 - t41 * t9 / 0.2e1 + t26 * t43 + t23 * (-mrSges(6,1) * t47 + mrSges(6,2) * t48) + t19 * t14 + t3 * t21 + t4 * t22 + (t30 / 0.2e1 - t29 / 0.2e1 + t65 * mrSges(5,1) + t114) * t72 - (t8 * mrSges(6,3) + Ifges(6,4) * t48 + Ifges(6,2) * t47 + Ifges(6,6) * t147) * t34 + (-t7 * mrSges(6,3) + Ifges(6,1) * t48 + Ifges(6,4) * t47 + Ifges(6,5) * t147) * t33 + m(5) * (-t23 * t111 + t65 * t115 + t24 * t140 - t16 * t25 + t17 * t26 + t31 * t74) + t141 * t74 + t142 * t25 + m(4) * (qJD(2) * t145 - t65 * t116 + t66 * t139 - t45 * t74 + t46 * t75 + t91) + 0.2e1 * m(3) * (t117 * t96 + t91) + m(6) * (t1 * t8 + t12 * t25 + t19 * t23 + t2 * t7 + t3 * t6 + t4 * t5) + (t1 * t47 - t2 * t48 - t5 * t40 - t6 * t41) * mrSges(6,3) + t12 * (mrSges(6,1) * t41 + mrSges(6,2) * t40) + t36 * (Ifges(6,4) * t40 - Ifges(6,2) * t41) / 0.2e1 + (Ifges(6,1) * t40 - Ifges(6,4) * t41) * t148 + ((-t23 * mrSges(5,1) - t24 * mrSges(5,2)) * t98 + 0.2e1 * qJD(2) * t69 + t110 * mrSges(4,3)) * t99 + (t23 * t73 - t24 * t72) * mrSges(5,3) + t65 * mrSges(5,2) * t73; t14 * t144 - t64 * t43 - t142 * t62 + t151 * t22 + t150 * t21 + (t108 * t33 - t34 * t77) * mrSges(6,3) + (t1 * t77 - t108 * t2 - t12 * t62 + t23 * t144 + t150 * t6 + t151 * t5) * m(6) - m(4) * t110 + (-t132 + (t100 * t24 + t23 * t97) * t98 + t16 * t62 - t17 * t64) * m(5) + (-t99 * t69 + (-t101 * t78 - t141 * t98) * t102 + (-t128 * t46 + t45 * t136 - t145) * m(4) - m(5) * t31 * t136 + (-m(3) * qJ(2) - mrSges(3,3)) * t138 * qJD(1)) * qJD(1); -t100 * t14 + t54 * t21 - t53 * t22 + ((-t103 * t21 - t104 * t22) * qJD(5) + (t103 * t33 - t104 * t34) * mrSges(6,3)) * t97 + m(5) * (t24 * t97 - t133) + (-t141 * t101 + (t100 * t43 + t142 * t97 + t78) * t98 - m(5) * (-t100 * t17 * t98 + t101 * t31 + t16 * t144) + (t101 * t45 + t46 * t98 + qJD(2)) * m(4)) * t130 + (-t133 + (t1 * t104 - t103 * t2 + (-t103 * t6 - t104 * t5) * qJD(5)) * t97 - t5 * t53 + t6 * t54 + t12 * t144 * t130) * m(6); t61 * t43 - t142 * t63 + (-t33 * mrSges(6,3) + t152 * t21) * t104 + (-t34 * mrSges(6,3) - t152 * t22) * t103 + (t1 * t103 + t104 * t2 - t12 * t63 + t152 * (-t103 * t5 + t104 * t6)) * m(6) + (t16 * t63 + t17 * t61 + t65) * m(5); t30 - t29 - t12 * (mrSges(6,1) * t38 + mrSges(6,2) * t36) - t38 * (Ifges(6,1) * t36 - t146) / 0.2e1 + t9 * t148 - t152 * (Ifges(6,5) * t36 - Ifges(6,6) * t38) / 0.2e1 - t5 * t21 + t6 * t22 + (t36 * t5 + t38 * t6) * mrSges(6,3) + t114 - (-Ifges(6,2) * t38 + t10 + t35) * t36 / 0.2e1;];
tauc = t15(:);
