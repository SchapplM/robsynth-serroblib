% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:05
% DurationCPUTime: 2.38s
% Computational Cost: add. (1937->250), mult. (5122->369), div. (0->0), fcn. (3643->6), ass. (0->119)
t107 = qJD(4) + qJD(5);
t111 = sin(qJ(5));
t113 = cos(qJ(5));
t108 = sin(pkin(8));
t109 = cos(pkin(8));
t114 = cos(qJ(4));
t136 = qJD(1) * t114;
t112 = sin(qJ(4));
t137 = qJD(1) * t112;
t70 = -t108 * t137 - t109 * t136;
t71 = t108 * t136 - t109 * t137;
t116 = t111 * t70 + t113 * t71;
t123 = -t111 * t71 + t113 * t70;
t36 = Ifges(6,4) * t123;
t13 = Ifges(6,1) * t116 + Ifges(6,5) * t107 + t36;
t146 = Ifges(6,4) * t116;
t78 = -t108 * t112 - t109 * t114;
t68 = t78 * qJD(4);
t63 = qJD(1) * t68;
t79 = t108 * t114 - t109 * t112;
t69 = t79 * qJD(4);
t64 = qJD(1) * t69;
t16 = qJD(5) * t123 - t111 * t64 + t113 * t63;
t17 = -qJD(5) * t116 - t111 * t63 - t113 * t64;
t130 = qJD(1) * qJD(2);
t125 = t114 * t130;
t126 = t112 * t130;
t131 = qJD(4) * t114;
t132 = qJD(4) * t112;
t139 = qJD(1) * t108;
t89 = qJ(2) * t139 + qJD(3);
t74 = -pkin(6) * t139 + t89;
t145 = -pkin(6) + qJ(2);
t85 = t145 * t109;
t81 = qJD(1) * t85;
t25 = t108 * t126 + t109 * t125 + t74 * t131 - t132 * t81;
t19 = -pkin(7) * t64 + t25;
t48 = t112 * t74 + t114 * t81;
t26 = -qJD(4) * t48 + t108 * t125 - t109 * t126;
t20 = -pkin(7) * t63 + t26;
t29 = pkin(7) * t70 + t48;
t144 = t111 * t29;
t47 = -t112 * t81 + t114 * t74;
t28 = -pkin(7) * t71 + t47;
t27 = qJD(4) * pkin(4) + t28;
t6 = t113 * t27 - t144;
t2 = qJD(5) * t6 + t111 * t20 + t113 * t19;
t142 = t113 * t29;
t7 = t111 * t27 + t142;
t3 = -qJD(5) * t7 - t111 * t19 + t113 * t20;
t138 = qJD(1) * t109;
t67 = -qJD(1) * pkin(1) - pkin(2) * t138 - qJ(3) * t139 + qJD(2);
t58 = pkin(3) * t138 - t67;
t37 = -pkin(4) * t70 + t58;
t166 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t16 + Ifges(6,6) * t17 - (Ifges(6,5) * t123 - Ifges(6,6) * t116) * t107 / 0.2e1 + (t116 * t7 + t123 * t6) * mrSges(6,3) - (-Ifges(6,2) * t116 + t13 + t36) * t123 / 0.2e1 - t37 * (mrSges(6,1) * t116 + mrSges(6,2) * t123) - (Ifges(6,1) * t123 - t146) * t116 / 0.2e1;
t165 = mrSges(4,2) + mrSges(3,3);
t12 = Ifges(6,2) * t123 + Ifges(6,6) * t107 + t146;
t163 = t12 / 0.2e1;
t162 = t68 / 0.2e1;
t151 = -t70 / 0.2e1;
t161 = t123 / 0.2e1;
t18 = -mrSges(6,1) * t123 + mrSges(6,2) * t116;
t156 = -m(6) * t37 - t18;
t152 = t116 / 0.2e1;
t150 = t71 / 0.2e1;
t147 = Ifges(5,4) * t71;
t84 = t145 * t108;
t56 = t112 * t84 + t114 * t85;
t115 = qJD(1) ^ 2;
t141 = qJ(2) * t115;
t135 = qJD(2) * t112;
t134 = qJD(2) * t114;
t133 = qJD(3) * t108;
t129 = -t109 * pkin(2) - t108 * qJ(3) - pkin(1);
t105 = t108 ^ 2;
t128 = t105 * t130;
t106 = t109 ^ 2;
t127 = t106 * t130;
t124 = qJD(1) * t133;
t55 = -t112 * t85 + t114 * t84;
t72 = t109 * pkin(3) - t129;
t120 = qJ(2) * t127;
t119 = t64 * mrSges(5,1) + t63 * mrSges(5,2);
t118 = -t17 * mrSges(6,1) + t16 * mrSges(6,2);
t34 = -pkin(7) * t79 + t55;
t35 = pkin(7) * t78 + t56;
t10 = -t111 * t35 + t113 * t34;
t11 = t111 * t34 + t113 * t35;
t49 = -t111 * t79 + t113 * t78;
t50 = t111 * t78 + t113 * t79;
t83 = t111 * t114 + t112 * t113;
t82 = -t111 * t112 + t113 * t114;
t32 = t108 * t135 + t109 * t134 + t84 * t131 - t132 * t85;
t80 = (-mrSges(4,1) * t109 - mrSges(4,3) * t108) * qJD(1);
t33 = -qJD(4) * t56 + t108 * t134 - t109 * t135;
t99 = t106 * t141;
t90 = qJ(2) * t128;
t65 = Ifges(5,4) * t70;
t60 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t71;
t59 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t70;
t57 = pkin(4) * t69 + t133;
t54 = pkin(4) * t64 + t124;
t53 = t107 * t83;
t52 = t107 * t82;
t51 = -pkin(4) * t78 + t72;
t46 = -mrSges(5,1) * t70 + mrSges(5,2) * t71;
t39 = Ifges(5,1) * t71 + Ifges(5,5) * qJD(4) + t65;
t38 = Ifges(5,2) * t70 + Ifges(5,6) * qJD(4) + t147;
t31 = mrSges(6,1) * t107 - mrSges(6,3) * t116;
t30 = -mrSges(6,2) * t107 + mrSges(6,3) * t123;
t24 = -pkin(7) * t68 + t33;
t23 = -pkin(7) * t69 + t32;
t22 = -qJD(5) * t50 - t111 * t68 - t113 * t69;
t21 = qJD(5) * t49 - t111 * t69 + t113 * t68;
t9 = t113 * t28 - t144;
t8 = -t111 * t28 - t142;
t5 = -qJD(5) * t11 - t111 * t23 + t113 * t24;
t4 = qJD(5) * t10 + t111 * t24 + t113 * t23;
t1 = [(m(5) * (qJD(1) * t72 + t58) - 0.2e1 * t80 + t46) * t133 + 0.2e1 * t165 * (t127 + t128) + 0.2e1 * m(3) * (t90 + t120) + t51 * t118 + t72 * t119 + (t150 * t68 + t63 * t79) * Ifges(5,1) + (t25 * t78 - t26 * t79 - t47 * t68 - t48 * t69 - t55 * t63 - t56 * t64) * mrSges(5,3) + (t151 * t69 - t78 * t64) * Ifges(5,2) + (-t69 * t150 + t162 * t70 + t78 * t63 - t64 * t79) * Ifges(5,4) + t107 * (Ifges(6,5) * t21 + Ifges(6,6) * t22) / 0.2e1 - t69 * t38 / 0.2e1 + t54 * (-mrSges(6,1) * t49 + mrSges(6,2) * t50) + t57 * t18 + t32 * t59 + t33 * t60 + t37 * (-mrSges(6,1) * t22 + mrSges(6,2) * t21) + t4 * t30 + t5 * t31 + t58 * (mrSges(5,1) * t69 + mrSges(5,2) * t68) + qJD(4) * (Ifges(5,5) * t68 - Ifges(5,6) * t69) / 0.2e1 + t21 * t13 / 0.2e1 + m(4) * (0.2e1 * t120 + t90 + (t89 * qJD(2) + (-qJD(1) * t129 - t67) * qJD(3)) * t108) + (t152 * t21 + t16 * t50) * Ifges(6,1) + m(6) * (t10 * t3 + t11 * t2 + t37 * t57 + t4 * t7 + t5 * t6 + t51 * t54) + (t22 * t152 + t16 * t49 + t161 * t21 + t17 * t50) * Ifges(6,4) + (t161 * t22 + t17 * t49) * Ifges(6,2) + t39 * t162 + t22 * t163 + (-t10 * t16 + t11 * t17 + t2 * t49 - t21 * t6 + t22 * t7 - t3 * t50) * mrSges(6,3) + m(5) * (t25 * t56 + t26 * t55 + t32 * t48 + t33 * t47) + (-mrSges(5,1) * t78 + mrSges(5,2) * t79) * t124; t123 * t30 - t116 * t31 + t70 * t59 - t71 * t60 + (-m(4) - m(5)) * t124 - m(3) * (t105 * t141 + t99) - m(4) * (t139 * t89 + t99) - m(5) * (t47 * t71 - t48 * t70) - t118 - t119 + t165 * t115 * (-t105 - t106) + (-t116 * t6 + t123 * t7 - t54) * m(6); t52 * t30 - t53 * t31 + (-t112 * t60 + t114 * t59) * qJD(4) + (-t16 * t82 + t17 * t83) * mrSges(6,3) + (-t112 * t64 - t114 * t63) * mrSges(5,3) + m(5) * (t112 * t25 + t114 * t26 + (-t112 * t47 + t114 * t48) * qJD(4)) + m(6) * (t2 * t83 + t3 * t82 + t52 * t7 - t53 * t6) + (-m(5) * t58 - t46 + t80 + (qJD(2) + t67) * m(4) + t156) * t139; (t47 * t70 + t48 * t71) * mrSges(5,3) + t116 * t163 + (-Ifges(5,2) * t71 + t39 + t65) * t151 - t58 * (mrSges(5,1) * t71 + mrSges(5,2) * t70) - qJD(4) * (Ifges(5,5) * t70 - Ifges(5,6) * t71) / 0.2e1 + Ifges(5,5) * t63 - Ifges(5,6) * t64 - t47 * t59 + t48 * t60 - t9 * t30 - t8 * t31 - t25 * mrSges(5,2) + t26 * mrSges(5,1) + ((t111 * t17 - t113 * t16) * mrSges(6,3) + m(6) * (t111 * t2 + t113 * t3) + t156 * t71 + (-t111 * t31 + t113 * t30 + m(6) * (-t111 * t6 + t113 * t7)) * qJD(5)) * pkin(4) - t71 * (Ifges(5,1) * t70 - t147) / 0.2e1 - m(6) * (t6 * t8 + t7 * t9) + t38 * t150 + t166; t12 * t152 - t6 * t30 + t7 * t31 + t166;];
tauc = t1(:);
