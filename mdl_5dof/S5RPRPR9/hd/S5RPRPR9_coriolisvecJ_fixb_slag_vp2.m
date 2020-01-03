% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:23:48
% DurationCPUTime: 1.97s
% Computational Cost: add. (1488->278), mult. (3497->378), div. (0->0), fcn. (1740->6), ass. (0->145)
t84 = cos(qJ(5));
t139 = Ifges(6,4) * t84;
t82 = sin(qJ(5));
t101 = Ifges(6,1) * t82 + t139;
t104 = mrSges(6,1) * t84 - mrSges(6,2) * t82;
t85 = cos(qJ(3));
t77 = t85 * qJD(2);
t71 = sin(pkin(8)) * pkin(1) + pkin(6);
t61 = t71 * qJD(1);
t108 = pkin(4) * qJD(1) + t61;
t83 = sin(qJ(3));
t94 = t108 * t83;
t31 = t77 - t94;
t170 = -t31 + qJD(4);
t86 = -pkin(3) - pkin(7);
t20 = t86 * qJD(3) + t170;
t110 = -cos(pkin(8)) * pkin(1) - pkin(2);
t93 = -qJ(4) * t83 + t110;
t42 = t86 * t85 + t93;
t25 = t42 * qJD(1);
t5 = t20 * t84 - t25 * t82;
t6 = t20 * t82 + t25 * t84;
t105 = t5 * t82 - t6 * t84;
t135 = Ifges(6,6) * t84;
t138 = Ifges(6,5) * t82;
t147 = -t84 / 0.2e1;
t148 = -t82 / 0.2e1;
t118 = t83 * qJD(1);
t70 = qJD(5) + t118;
t149 = -t70 / 0.2e1;
t122 = qJD(3) * t84;
t125 = qJD(1) * t85;
t58 = -t82 * t125 + t122;
t141 = Ifges(6,4) * t58;
t57 = -qJD(3) * t82 - t84 * t125;
t15 = Ifges(6,2) * t57 + Ifges(6,6) * t70 + t141;
t150 = -t58 / 0.2e1;
t151 = -t57 / 0.2e1;
t54 = Ifges(6,4) * t57;
t16 = Ifges(6,1) * t58 + Ifges(6,5) * t70 + t54;
t117 = qJD(3) * qJ(4);
t124 = qJD(2) * t83;
t44 = t61 * t85 + t124;
t37 = -t44 - t117;
t75 = pkin(4) * t125;
t24 = -t37 + t75;
t140 = Ifges(6,4) * t82;
t99 = Ifges(6,2) * t84 + t140;
t175 = t105 * mrSges(6,3) + t101 * t150 + t24 * t104 + t15 * t147 + t16 * t148 + (t135 + t138) * t149 + t99 * t151;
t174 = -t125 / 0.2e1;
t74 = pkin(3) * t118;
t69 = qJD(3) * t74;
t120 = qJD(4) * t83;
t97 = pkin(7) * t83 - qJ(4) * t85;
t89 = t97 * qJD(3) - t120;
t22 = t89 * qJD(1) + t69;
t23 = (t108 * t85 + t124) * qJD(3);
t1 = qJD(5) * t5 + t22 * t84 + t23 * t82;
t2 = -qJD(5) * t6 - t22 * t82 + t23 * t84;
t161 = t1 * t82 + t2 * t84;
t173 = m(6) * t161;
t163 = mrSges(5,1) + mrSges(4,3);
t172 = mrSges(5,2) - mrSges(4,1);
t164 = qJD(3) / 0.2e1;
t165 = -qJD(3) / 0.2e1;
t166 = qJD(1) / 0.2e1;
t53 = -pkin(3) * t85 + t93;
t38 = t53 * qJD(1);
t63 = t110 * qJD(1);
t171 = t38 * mrSges(5,2) + t44 * mrSges(4,3) + Ifges(4,6) * t164 + (Ifges(4,4) * t83 + t85 * Ifges(4,2)) * t166 + Ifges(5,5) * t165 - (-Ifges(5,6) * t83 - t85 * Ifges(5,3)) * qJD(1) / 0.2e1 - t37 * mrSges(5,1) - t63 * mrSges(4,1) + t175;
t169 = -Ifges(4,1) / 0.2e1;
t168 = Ifges(4,4) * t174;
t131 = t61 * t83;
t72 = qJD(3) * t77;
t28 = -t72 + (-qJD(4) + t131) * qJD(3);
t167 = m(5) * t28;
t162 = t5 * t84 + t6 * t82;
t43 = -t77 + t131;
t160 = -t43 - qJD(4);
t119 = qJD(5) * t85;
t114 = t84 * t119;
t116 = qJD(3) * qJD(5);
t123 = qJD(3) * t83;
t29 = -t82 * t116 + (t82 * t123 - t114) * qJD(1);
t115 = t82 * t119;
t30 = -t84 * t116 + (t83 * t122 + t115) * qJD(1);
t159 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t29 + Ifges(6,6) * t30;
t33 = -mrSges(6,2) * t70 + mrSges(6,3) * t57;
t34 = mrSges(6,1) * t70 - mrSges(6,3) * t58;
t95 = t84 * t33 - t82 * t34;
t121 = qJD(3) * t85;
t109 = qJD(1) * t121;
t17 = mrSges(6,1) * t109 - mrSges(6,3) * t29;
t18 = -mrSges(6,2) * t109 + mrSges(6,3) * t30;
t96 = t84 * t17 + t82 * t18;
t158 = -t95 * qJD(5) - t96;
t157 = -m(4) * t44 + m(5) * t37;
t127 = t172 * qJD(3) + t163 * t118;
t35 = -qJD(3) * pkin(3) - t160;
t156 = m(4) * t43 + m(5) * t35 + t127;
t154 = 0.2e1 * t71;
t153 = m(4) / 0.2e1;
t152 = t15 / 0.2e1;
t142 = pkin(4) + t71;
t137 = Ifges(6,5) * t84;
t136 = Ifges(6,6) * t82;
t64 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t125;
t65 = -mrSges(5,1) * t125 - qJD(3) * mrSges(5,3);
t129 = t64 - t65;
t21 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t128 = -t65 + t21;
t113 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t112 = -Ifges(4,6) / 0.2e1 + Ifges(5,5) / 0.2e1;
t111 = -0.3e1 / 0.2e1 * Ifges(5,6) - 0.3e1 / 0.2e1 * Ifges(4,4);
t56 = t142 * t85;
t103 = mrSges(6,1) * t82 + mrSges(6,2) * t84;
t102 = Ifges(6,1) * t84 - t140;
t100 = -Ifges(6,2) * t82 + t139;
t55 = t142 * t83;
t13 = t42 * t84 + t55 * t82;
t12 = -t42 * t82 + t55 * t84;
t90 = -t85 * t117 - t120;
t88 = t38 * mrSges(5,3) + t6 * mrSges(6,2) - t70 * Ifges(6,3) - t58 * Ifges(6,5) - t57 * Ifges(6,6) + t118 * t169 + Ifges(4,5) * t165 + t168 + Ifges(5,4) * t164 + (-Ifges(5,2) * t83 - Ifges(5,6) * t85) * t166 - t35 * mrSges(5,1) - t43 * mrSges(4,3) - t5 * mrSges(6,1) - t63 * mrSges(4,2);
t76 = pkin(3) * t123;
t68 = Ifges(6,3) * t109;
t60 = -qJ(4) * t125 + t74;
t59 = (mrSges(5,2) * t85 - mrSges(5,3) * t83) * qJD(1);
t48 = t76 + t90;
t47 = qJD(3) * t56;
t46 = t142 * t123;
t45 = t97 * qJD(1) + t74;
t41 = t44 * qJD(3);
t40 = -t61 * t123 + t72;
t39 = t90 * qJD(1) + t69;
t36 = t76 + t89;
t32 = t44 + t75;
t19 = t72 + (qJD(4) - t94) * qJD(3);
t11 = t32 * t82 + t45 * t84;
t10 = t32 * t84 - t45 * t82;
t9 = -mrSges(6,1) * t30 + mrSges(6,2) * t29;
t8 = t29 * Ifges(6,1) + t30 * Ifges(6,4) + Ifges(6,5) * t109;
t7 = t29 * Ifges(6,4) + t30 * Ifges(6,2) + Ifges(6,6) * t109;
t4 = -t13 * qJD(5) - t36 * t82 + t47 * t84;
t3 = t12 * qJD(5) + t36 * t84 + t47 * t82;
t14 = [t12 * t17 + t13 * t18 - t46 * t21 + t3 * t33 + t4 * t34 + t48 * t59 + t56 * t9 + m(6) * (t1 * t13 + t12 * t2 + t19 * t56 - t24 * t46 + t3 * t6 + t4 * t5) + m(5) * (t38 * t48 + t39 * t53) + (-t39 * mrSges(5,3) + t68 / 0.2e1 + ((m(5) / 0.2e1 + t153) * t154 + t163) * t41 + ((t110 * mrSges(4,1) - t53 * mrSges(5,2) + t111 * t83) * qJD(1) + t112 * qJD(3) + (-t129 + t157) * t71 - t171) * qJD(3) + t159) * t83 + (t19 * t104 + t39 * mrSges(5,2) + t7 * t147 + t8 * t148 + t40 * mrSges(4,3) - t28 * mrSges(5,1) - t29 * t101 / 0.2e1 - t30 * t99 / 0.2e1 + (-t1 * t84 + t2 * t82) * mrSges(6,3) + (-t167 / 0.2e1 + t40 * t153) * t154 + (t70 * (t136 - t137) / 0.2e1 + t100 * t151 + t102 * t150 - t24 * t103 + t16 * t147 + t82 * t152 + t162 * mrSges(6,3)) * qJD(5) + (t113 * qJD(3) - t88 + t156 * t71 + ((-t138 / 0.2e1 - t135 / 0.2e1 - t111) * t85 + t110 * mrSges(4,2) - t53 * mrSges(5,3) + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + Ifges(6,3) / 0.2e1) * t83) * qJD(1)) * qJD(3)) * t85; m(6) * (-t6 * t114 + t5 * t115) + (-m(4) - m(5)) * t41 * t85 + t163 * qJD(1) * qJD(3) * (-t83 ^ 2 - t85 ^ 2) + (-t173 + (m(6) * t24 + t128 - t157 + t64) * qJD(3) + t158) * t85 + (t9 + m(4) * t40 - t167 + m(6) * t19 + (m(6) * t162 + t82 * t33 + t84 * t34 + t156) * qJD(3)) * t83; t29 * t102 / 0.2e1 + t30 * t100 / 0.2e1 + t19 * t103 + t84 * t8 / 0.2e1 + t7 * t148 - t60 * t59 - t11 * t33 - t10 * t34 - t40 * mrSges(4,2) - t28 * mrSges(5,3) - t31 * t21 + qJ(4) * t9 - t127 * t44 + t129 * t43 + t172 * t41 + t128 * qJD(4) - t161 * mrSges(6,3) + t175 * qJD(5) + ((t168 + Ifges(5,6) * t174 + (-pkin(3) * mrSges(5,1) + t137 / 0.2e1 - t136 / 0.2e1 + t113) * qJD(3) + t88) * t85 + ((Ifges(5,6) / 0.2e1 + Ifges(4,4) / 0.2e1) * t118 + (-qJ(4) * mrSges(5,1) + t112) * qJD(3) + (Ifges(4,2) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t169) * t125 + t171) * t83) * qJD(1) + (t19 * qJ(4) - t10 * t5 - t11 * t6 + t170 * t24) * m(6) + (t96 + t173 + (-m(6) * t105 + t95) * qJD(5)) * t86 + (-pkin(3) * t41 - qJ(4) * t28 + t160 * t37 - t35 * t44 - t38 * t60) * m(5); -t128 * qJD(3) + (mrSges(5,1) * t121 + (t59 + t95) * t83) * qJD(1) + (-t24 * qJD(3) - t70 * t105 + t161) * m(6) + (qJD(3) * t37 + t38 * t118 + t41) * m(5) - t158; t68 - t24 * (mrSges(6,1) * t58 + mrSges(6,2) * t57) + (Ifges(6,1) * t57 - t141) * t150 + t58 * t152 + (Ifges(6,5) * t57 - Ifges(6,6) * t58) * t149 - t5 * t33 + t6 * t34 + (t5 * t57 + t58 * t6) * mrSges(6,3) + (-Ifges(6,2) * t58 + t16 + t54) * t151 + t159;];
tauc = t14(:);
