% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:40
% EndTime: 2019-12-05 15:49:50
% DurationCPUTime: 2.06s
% Computational Cost: add. (1710->274), mult. (4565->400), div. (0->0), fcn. (3195->10), ass. (0->145)
t86 = cos(qJ(4));
t125 = qJD(2) * t86;
t77 = Ifges(5,4) * t125;
t165 = -t77 / 0.2e1;
t164 = qJD(4) / 0.2e1;
t83 = sin(qJ(4));
t126 = qJD(2) * t83;
t163 = t126 / 0.2e1;
t108 = Ifges(5,5) * t164;
t79 = sin(pkin(5));
t127 = qJD(1) * t79;
t84 = sin(qJ(2));
t114 = t84 * t127;
t87 = cos(qJ(2));
t113 = t87 * t127;
t67 = qJD(2) * pkin(2) + t113;
t78 = sin(pkin(10));
t80 = cos(pkin(10));
t41 = t114 * t80 + t67 * t78;
t39 = qJD(2) * pkin(7) + t41;
t81 = cos(pkin(5));
t71 = qJD(1) * t81 + qJD(3);
t24 = -t39 * t83 + t71 * t86;
t68 = t78 * t114;
t40 = t67 * t80 - t68;
t38 = -qJD(2) * pkin(3) - t40;
t82 = sin(qJ(5));
t142 = Ifges(6,4) * t82;
t85 = cos(qJ(5));
t100 = Ifges(6,1) * t85 - t142;
t101 = mrSges(6,1) * t82 + mrSges(6,2) * t85;
t25 = t39 * t86 + t71 * t83;
t23 = qJD(4) * pkin(8) + t25;
t92 = -pkin(4) * t86 - pkin(8) * t83 - pkin(3);
t26 = qJD(2) * t92 - t40;
t5 = -t23 * t82 + t26 * t85;
t6 = t23 * t85 + t26 * t82;
t103 = t5 * t85 + t6 * t82;
t139 = Ifges(6,6) * t82;
t140 = Ifges(6,5) * t85;
t148 = t85 / 0.2e1;
t149 = -t82 / 0.2e1;
t124 = qJD(4) * t82;
t63 = t126 * t85 + t124;
t152 = t63 / 0.2e1;
t22 = -qJD(4) * pkin(4) - t24;
t143 = Ifges(6,4) * t63;
t123 = qJD(4) * t85;
t62 = -t126 * t82 + t123;
t73 = qJD(5) - t125;
t29 = Ifges(6,2) * t62 + Ifges(6,6) * t73 + t143;
t61 = Ifges(6,4) * t62;
t30 = Ifges(6,1) * t63 + Ifges(6,5) * t73 + t61;
t141 = Ifges(6,4) * t85;
t98 = -Ifges(6,2) * t82 + t141;
t89 = -t103 * mrSges(6,3) + t30 * t148 + t29 * t149 + t73 * (-t139 + t140) / 0.2e1 + t100 * t152 + t62 * t98 / 0.2e1 + t22 * t101;
t162 = -t38 * mrSges(5,2) + t24 * mrSges(5,3) - Ifges(5,1) * t163 - t108 + t165 - t89;
t161 = qJD(2) * t79;
t55 = (t78 * t87 + t80 * t84) * t79;
t50 = qJD(1) * t55;
t105 = pkin(4) * t83 - pkin(8) * t86;
t66 = t105 * qJD(4);
t27 = (t66 + t50) * qJD(2);
t119 = t24 * qJD(4);
t93 = t78 * t84 - t80 * t87;
t52 = t93 * t161;
t45 = qJD(1) * t52;
t7 = -t45 * t86 + t119;
t1 = qJD(5) * t5 + t27 * t82 + t7 * t85;
t2 = -qJD(5) * t6 + t27 * t85 - t7 * t82;
t104 = t1 * t85 - t2 * t82;
t116 = qJD(4) * qJD(5);
t121 = qJD(5) * t82;
t122 = qJD(4) * t86;
t46 = t85 * t116 + (-t121 * t83 + t122 * t85) * qJD(2);
t120 = qJD(5) * t85;
t47 = -t82 * t116 + (-t120 * t83 - t122 * t82) * qJD(2);
t160 = -t2 * mrSges(6,1) + t1 * mrSges(6,2) - Ifges(6,5) * t46 - Ifges(6,6) * t47;
t158 = 2 * m(5);
t157 = t46 / 0.2e1;
t156 = t47 / 0.2e1;
t53 = t113 * t80 - t68;
t155 = -t53 / 0.2e1;
t154 = -t62 / 0.2e1;
t153 = -t63 / 0.2e1;
t151 = -t73 / 0.2e1;
t75 = pkin(2) * t78 + pkin(7);
t150 = t75 / 0.2e1;
t147 = pkin(2) * t80;
t36 = t55 * t83 - t81 * t86;
t118 = t25 * qJD(4);
t8 = -t45 * t83 + t118;
t144 = t36 * t8;
t51 = qJD(2) * t55;
t44 = qJD(1) * t51;
t54 = t93 * t79;
t136 = t44 * t54;
t132 = t82 * t86;
t131 = t85 * t86;
t111 = mrSges(5,3) * t126;
t130 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t62 - mrSges(6,2) * t63 - t111;
t128 = Ifges(5,6) * qJD(4);
t117 = qJD(4) * qJD(2);
t115 = t8 * t150;
t112 = qJD(4) * t75 * t83;
t110 = mrSges(5,3) * t125;
t109 = t83 * t117;
t107 = -t128 / 0.2e1;
t102 = mrSges(6,1) * t85 - mrSges(6,2) * t82;
t99 = Ifges(6,1) * t82 + t141;
t97 = Ifges(6,2) * t85 + t142;
t96 = Ifges(6,5) * t82 + Ifges(6,6) * t85;
t31 = mrSges(6,1) * t109 - mrSges(6,3) * t46;
t32 = -mrSges(6,2) * t109 + mrSges(6,3) * t47;
t95 = -t82 * t31 + t85 * t32;
t37 = t55 * t86 + t81 * t83;
t18 = t37 * t85 + t54 * t82;
t17 = -t37 * t82 + t54 * t85;
t48 = -mrSges(6,2) * t73 + mrSges(6,3) * t62;
t49 = mrSges(6,1) * t73 - mrSges(6,3) * t63;
t94 = -t82 * t48 - t85 * t49;
t60 = t92 - t147;
t34 = t131 * t75 + t60 * t82;
t33 = -t132 * t75 + t60 * t85;
t91 = t6 * mrSges(6,2) - t73 * Ifges(6,3) - t63 * Ifges(6,5) - t62 * Ifges(6,6) + t128 / 0.2e1 + (Ifges(5,4) * t83 + t86 * Ifges(5,2)) * qJD(2) / 0.2e1 - t38 * mrSges(5,1) - t5 * mrSges(6,1);
t76 = -pkin(3) - t147;
t72 = Ifges(6,3) * t109;
t70 = -qJD(4) * mrSges(5,2) + t110;
t65 = t105 * qJD(2);
t64 = (-mrSges(5,1) * t86 + mrSges(5,2) * t83) * qJD(2);
t58 = (mrSges(5,1) * t83 + mrSges(5,2) * t86) * t117;
t21 = -mrSges(6,1) * t47 + mrSges(6,2) * t46;
t20 = t131 * t53 + t50 * t82;
t19 = -t132 * t53 + t50 * t85;
t16 = -qJD(4) * t36 - t52 * t86;
t15 = qJD(4) * t37 - t52 * t83;
t14 = Ifges(6,1) * t46 + Ifges(6,4) * t47 + Ifges(6,5) * t109;
t13 = Ifges(6,4) * t46 + Ifges(6,2) * t47 + Ifges(6,6) * t109;
t12 = -qJD(5) * t34 + t112 * t82 + t66 * t85;
t11 = qJD(5) * t33 - t112 * t85 + t66 * t82;
t10 = t24 * t85 + t65 * t82;
t9 = -t24 * t82 + t65 * t85;
t4 = qJD(5) * t17 + t16 * t85 + t51 * t82;
t3 = -qJD(5) * t18 - t16 * t82 + t51 * t85;
t28 = [t16 * t70 + t17 * t31 + t18 * t32 + t36 * t21 + t3 * t49 + t4 * t48 + t51 * t64 + t54 * t58 - t130 * t15 + m(5) * (-t15 * t24 + t16 * t25 + t37 * t7 + t38 * t51 + t136 + t144) + m(4) * (-t40 * t51 - t41 * t52 - t45 * t55 + t136) + m(6) * (t1 * t18 + t15 * t22 + t17 * t2 + t3 * t5 + t4 * t6 + t144) + (-t51 * mrSges(4,1) + t52 * mrSges(4,2) + (t36 * t86 - t37 * t83) * qJD(4) * mrSges(5,3) + (-mrSges(3,1) * t84 - mrSges(3,2) * t87) * t161) * qJD(2); t33 * t31 + t34 * t32 - t50 * t64 + t76 * t58 + (-t19 + t12) * t49 + (-t20 + t11) * t48 + (qJD(2) * t53 + t45) * mrSges(4,2) + (qJD(2) * t50 - t44) * mrSges(4,1) - m(6) * (t19 * t5 + t20 * t6) + m(6) * (t1 * t34 + t11 * t6 + t12 * t5 + t2 * t33) + (-t38 * t50 / 0.2e1 + t44 * t76 / 0.2e1) * t158 + (-t72 / 0.2e1 - t53 * t70 + t7 * mrSges(5,3) - t44 * mrSges(5,1) + (t150 * t7 + t155 * t25) * t158 + (0.3e1 / 0.2e1 * t77 + t108 + (-m(5) * t24 + m(6) * t22 - t130) * t75 - t162) * qJD(4) + t160) * t86 + (t14 * t148 + t13 * t149 + t75 * t21 + t44 * mrSges(5,2) + t100 * t157 + t98 * t156 + (mrSges(5,3) + t101) * t8 + t130 * t53 + (-t1 * t82 - t2 * t85) * mrSges(6,3) + 0.2e1 * (t155 * t22 + t115) * m(6) + (t24 * t53 / 0.2e1 + t115) * t158 + (t96 * t151 + t99 * t153 + t97 * t154 + t22 * t102 - t85 * t29 / 0.2e1 + t30 * t149 + (t5 * t82 - t6 * t85) * mrSges(6,3)) * qJD(5) + (-t75 * t70 + t107 + (-m(5) * t75 - mrSges(5,3)) * t25 + ((-0.3e1 / 0.2e1 * Ifges(5,4) + t140 / 0.2e1 - t139 / 0.2e1) * t83 + (-Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,1)) * t86) * qJD(2) - t91) * qJD(4)) * t83 + (t40 * t50 - t41 * t53 + (-t44 * t80 - t45 * t78) * pkin(2)) * m(4); (-t21 + (t48 * t85 - t49 * t82 - t110 + t70) * qJD(4) + m(5) * (-t8 + t118) + m(6) * (t123 * t6 - t124 * t5 - t8)) * t86 + (t94 * qJD(5) + (-t111 - t130) * qJD(4) + m(5) * (t7 - t119) + m(6) * (qJD(4) * t22 - t120 * t5 - t121 * t6 + t104) + t95) * t83; t99 * t157 + t97 * t156 + t13 * t148 + t82 * t14 / 0.2e1 - t24 * t70 - t10 * t48 - t9 * t49 - pkin(4) * t21 - t7 * mrSges(5,2) + (-mrSges(5,1) - t102) * t8 + t130 * t25 + t104 * mrSges(6,3) + t89 * qJD(5) + ((t25 * mrSges(5,3) + Ifges(5,4) * t163 + t164 * t96 + t107 + t91) * t83 + (t165 + t108 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t126 + t162) * t86) * qJD(2) + (-t8 * pkin(4) - t10 * t6 - t22 * t25 - t5 * t9) * m(6) + (t95 + m(6) * t104 + (-m(6) * t103 + t94) * qJD(5)) * pkin(8); t72 - t22 * (mrSges(6,1) * t63 + mrSges(6,2) * t62) + (Ifges(6,1) * t62 - t143) * t153 + t29 * t152 + (Ifges(6,5) * t62 - Ifges(6,6) * t63) * t151 - t5 * t48 + t6 * t49 + (t5 * t62 + t6 * t63) * mrSges(6,3) + (-Ifges(6,2) * t63 + t30 + t61) * t154 - t160;];
tauc = t28(:);
