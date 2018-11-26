% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:42:21
% EndTime: 2018-11-23 15:42:23
% DurationCPUTime: 2.23s
% Computational Cost: add. (1826->332), mult. (3657->439), div. (0->0), fcn. (1607->4), ass. (0->149)
t181 = -qJD(1) / 0.2e1;
t84 = sin(qJ(6));
t86 = cos(qJ(6));
t107 = mrSges(7,1) * t86 - mrSges(7,2) * t84;
t85 = sin(qJ(4));
t135 = qJD(1) * t85;
t127 = pkin(4) * t135 - qJD(2);
t83 = pkin(1) + qJ(3);
t87 = cos(qJ(4));
t113 = -qJ(5) * t87 + t83;
t94 = pkin(8) * t85 + t113;
t25 = t94 * qJD(1) + t127;
t162 = pkin(4) + pkin(8);
t77 = qJD(1) * qJ(2) + qJD(3);
t66 = -qJD(1) * pkin(7) + t77;
t116 = pkin(5) * qJD(1) - t66;
t42 = t116 * t87;
t182 = t42 + qJD(5);
t26 = -t162 * qJD(4) + t182;
t7 = -t25 * t84 + t26 * t86;
t8 = t25 * t86 + t26 * t84;
t108 = t7 * t84 - t8 * t86;
t54 = qJD(4) * t86 + t84 * t135;
t154 = Ifges(7,4) * t54;
t53 = -qJD(4) * t84 + t86 * t135;
t128 = t87 * qJD(1);
t68 = qJD(6) + t128;
t13 = Ifges(7,2) * t53 + Ifges(7,6) * t68 + t154;
t51 = Ifges(7,4) * t53;
t14 = Ifges(7,1) * t54 + Ifges(7,5) * t68 + t51;
t163 = t86 / 0.2e1;
t164 = -t84 / 0.2e1;
t126 = qJD(4) * qJ(5);
t58 = t85 * t66;
t41 = -pkin(5) * t135 + t58;
t34 = t41 + t126;
t185 = t108 * mrSges(7,3) + t34 * t107 - t13 * t163 + t14 * t164;
t140 = -mrSges(5,1) + mrSges(6,2);
t184 = mrSges(6,1) + mrSges(5,3);
t148 = Ifges(7,6) * t86;
t151 = Ifges(7,5) * t84;
t101 = t148 + t151;
t153 = Ifges(7,4) * t84;
t102 = Ifges(7,2) * t86 + t153;
t152 = Ifges(7,4) * t86;
t104 = Ifges(7,1) * t84 + t152;
t165 = t68 / 0.2e1;
t167 = t54 / 0.2e1;
t169 = t53 / 0.2e1;
t179 = qJD(4) / 0.2e1;
t180 = -qJD(4) / 0.2e1;
t35 = t113 * qJD(1) + t127;
t46 = -t58 - t126;
t67 = t83 * qJD(1) - qJD(2);
t183 = -t67 * mrSges(5,1) - t46 * mrSges(6,1) + t35 * mrSges(6,2) - Ifges(6,5) * t179 - Ifges(5,6) * t180 - t101 * t165 - t102 * t169 - t104 * t167 + t185 + ((-Ifges(5,4) - Ifges(6,6)) * t87 + (Ifges(5,2) + Ifges(6,3)) * t85) * t181;
t131 = qJD(4) * t87;
t178 = pkin(4) * t131 + t85 * t126;
t109 = t7 * t86 + t8 * t84;
t115 = -qJD(5) * t87 + qJD(3);
t136 = qJ(5) * t85;
t19 = ((t162 * t87 + t136) * qJD(4) + t115) * qJD(1);
t132 = qJD(4) * t85;
t125 = qJD(1) * qJD(2);
t70 = t87 * t125;
t22 = -t116 * t132 - t70;
t1 = t7 * qJD(6) + t19 * t86 + t22 * t84;
t2 = -t8 * qJD(6) - t19 * t84 + t22 * t86;
t177 = t1 * t84 + t2 * t86;
t123 = qJD(4) * qJD(6);
t130 = qJD(6) * t85;
t28 = -t84 * t123 + (t86 * t130 + t84 * t131) * qJD(1);
t29 = -t86 * t123 + (-t84 * t130 + t86 * t131) * qJD(1);
t176 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t28 + Ifges(7,6) * t29;
t175 = (t85 ^ 2 + t87 ^ 2) * t66;
t118 = qJD(1) * t132;
t17 = -mrSges(7,1) * t118 - mrSges(7,3) * t28;
t18 = mrSges(7,2) * t118 + mrSges(7,3) * t29;
t100 = t86 * t17 + t84 * t18;
t31 = -mrSges(7,2) * t68 + mrSges(7,3) * t53;
t32 = mrSges(7,1) * t68 - mrSges(7,3) * t54;
t99 = t86 * t31 - t84 * t32;
t173 = -t99 * qJD(6) - t100;
t172 = t28 / 0.2e1;
t171 = t29 / 0.2e1;
t170 = -t53 / 0.2e1;
t168 = -t54 / 0.2e1;
t166 = -t68 / 0.2e1;
t82 = qJ(2) - pkin(7);
t155 = pkin(5) - t82;
t150 = Ifges(7,5) * t86;
t149 = Ifges(7,6) * t84;
t36 = t66 * t132 - t70;
t147 = t36 * t82;
t146 = t36 * t87;
t142 = t84 * t87;
t141 = t86 * t87;
t37 = t85 * t125 + t66 * t131;
t62 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t135;
t64 = mrSges(6,1) * t135 - qJD(4) * mrSges(6,3);
t139 = t62 - t64;
t21 = -mrSges(7,1) * t53 + mrSges(7,2) * t54;
t138 = -t64 + t21;
t137 = qJD(4) * t140 + t128 * t184;
t55 = pkin(4) * t128 + qJ(5) * t135;
t133 = qJD(3) * t67;
t124 = qJD(1) * qJD(3);
t122 = t62 + t138;
t121 = -Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t120 = -Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1;
t119 = -0.3e1 / 0.2e1 * Ifges(6,6) - 0.3e1 / 0.2e1 * Ifges(5,4);
t117 = qJD(4) * t155;
t114 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t111 = t1 * t86 - t2 * t84;
t106 = mrSges(7,1) * t84 + mrSges(7,2) * t86;
t105 = Ifges(7,1) * t86 - t153;
t103 = -Ifges(7,2) * t84 + t152;
t98 = t84 * t31 + t86 * t32;
t44 = -qJD(4) * pkin(4) - t66 * t87 + qJD(5);
t97 = t44 * t85 - t46 * t87;
t78 = t85 * pkin(4);
t45 = t78 + t94;
t61 = t155 * t87;
t16 = t45 * t86 + t61 * t84;
t15 = -t45 * t84 + t61 * t86;
t93 = t98 + t137;
t92 = -m(7) * t108 + t99;
t73 = Ifges(6,6) * t135;
t91 = t35 * mrSges(6,3) + t8 * mrSges(7,2) - t68 * Ifges(7,3) - t54 * Ifges(7,5) - t53 * Ifges(7,6) + Ifges(5,5) * t180 + (t87 * Ifges(5,1) - Ifges(5,4) * t85) * t181 + Ifges(6,4) * t179 - Ifges(6,2) * t128 / 0.2e1 + t73 / 0.2e1 - t44 * mrSges(6,1) - t67 * mrSges(5,2) - t7 * mrSges(7,1);
t89 = qJD(1) ^ 2;
t60 = t155 * t85;
t59 = t113 + t78;
t57 = qJD(1) * (mrSges(5,1) * t85 + mrSges(5,2) * t87);
t56 = (-mrSges(6,2) * t85 - mrSges(6,3) * t87) * qJD(1);
t43 = pkin(8) * t128 + t55;
t40 = t115 + t178;
t39 = qJD(2) * t85 - t87 * t117;
t38 = -qJD(2) * t87 - t85 * t117;
t33 = -qJD(4) * qJD(5) - t37;
t30 = qJD(3) + (qJD(4) * pkin(8) - qJD(5)) * t87 + t178;
t27 = ((pkin(4) * t87 + t136) * qJD(4) + t115) * qJD(1);
t20 = (-pkin(5) * t128 + qJD(5)) * qJD(4) + t37;
t11 = t41 * t84 + t43 * t86;
t10 = t41 * t86 - t43 * t84;
t9 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t6 = t28 * Ifges(7,1) + t29 * Ifges(7,4) - Ifges(7,5) * t118;
t5 = t28 * Ifges(7,4) + t29 * Ifges(7,2) - Ifges(7,6) * t118;
t4 = -t16 * qJD(6) - t30 * t84 + t38 * t86;
t3 = t15 * qJD(6) + t30 * t86 + t38 * t84;
t12 = [qJD(3) * t57 + t15 * t17 + t16 * t18 + t39 * t21 + t3 * t31 + t4 * t32 + t40 * t56 - t60 * t9 + 0.2e1 * (qJD(3) * mrSges(4,3) + qJD(2) * t114) * qJD(1) + m(6) * (t27 * t59 + t35 * t40) + m(5) * (qJD(2) * t175 + t83 * t124 + t133) + m(4) * (qJD(2) * t77 + t133 + (qJ(2) * qJD(2) + qJD(3) * t83) * qJD(1)) + m(7) * (t1 * t16 + t15 * t2 - t20 * t60 + t3 * t8 + t34 * t39 + t4 * t7) + (-t27 * mrSges(6,3) + mrSges(5,2) * t124 + t184 * t36 - t137 * qJD(2) - m(5) * t147 + m(6) * (-qJD(2) * t44 - t147) + ((-m(6) * t46 + t139) * t82 + (t83 * mrSges(5,1) - t59 * mrSges(6,2) + t119 * t87) * qJD(1) + t120 * qJD(4) - t183) * qJD(4) + t176) * t87 + (t5 * t163 + t33 * mrSges(6,1) + t84 * t6 / 0.2e1 - t20 * t107 + t104 * t172 + t102 * t171 - t27 * mrSges(6,2) + mrSges(5,1) * t124 + t139 * qJD(2) + t111 * mrSges(7,3) + m(6) * (-qJD(2) * t46 - t33 * t82) + ((-t149 + t150) * t165 + t103 * t169 + t105 * t167 + t34 * t106 + t14 * t163 + t13 * t164 - t109 * mrSges(7,3)) * qJD(6) + (t121 * qJD(4) + t91 + (m(6) * t44 + t137) * t82 + (-t83 * mrSges(5,2) + t59 * mrSges(6,3) + (-t151 / 0.2e1 - t148 / 0.2e1 - t119) * t85 + (-Ifges(7,3) - 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(6,3) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t87) * qJD(1)) * qJD(4) + (m(5) * t82 - mrSges(5,3)) * t37) * t85; t84 * t17 - t86 * t18 + t98 * qJD(6) + m(7) * (t109 * qJD(6) - t111) - m(6) * t27 - t114 * t89 + ((-t77 - qJD(3)) * m(4) - t122 * t85 + t93 * t87 - m(7) * (-t7 * t141 - t8 * t142 + t34 * t85) - m(6) * (-t44 * t87 - t46 * t85) + (t140 * t87 + (mrSges(5,2) - mrSges(6,3)) * t85) * qJD(4) + (-qJD(3) - t175) * m(5)) * qJD(1); -t89 * mrSges(4,3) + t85 * t9 + t173 * t87 + m(7) * (t108 * qJD(6) * t87 - t1 * t142 - t2 * t141 + t20 * t85) + m(6) * (-t33 * t85 - t146) + m(5) * (t37 * t85 - t146) + (m(6) * t97 + (m(7) * t34 + t122) * t87 + (m(7) * t109 + t93) * t85) * qJD(4) + (-m(5) * t67 - m(6) * t35 - t56 - t57 + (-t67 + qJD(2)) * m(4) - t92) * qJD(1); t105 * t172 + t103 * t171 + t20 * t106 + t6 * t163 + t5 * t164 - t55 * t56 - t37 * mrSges(5,2) + t42 * t21 - t11 * t31 - t10 * t32 - t33 * mrSges(6,3) + qJ(5) * t9 - t100 * t162 + t140 * t36 + t138 * qJD(5) - t177 * mrSges(7,3) + (-t137 * t85 - t139 * t87) * t66 + (t101 * t166 + t102 * t170 + t104 * t168 - t162 * t92 + t185) * qJD(6) + ((-t91 - t73 / 0.2e1 + (pkin(4) * mrSges(6,1) - t150 / 0.2e1 + t149 / 0.2e1 + t121) * qJD(4) - Ifges(5,4) * t135 / 0.2e1) * t85 + ((-qJ(5) * mrSges(6,1) + t120) * qJD(4) + (-Ifges(6,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t135 + (Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t128 + t183) * t87) * qJD(1) + (qJ(5) * t20 - t10 * t7 - t11 * t8 - t177 * t162 + t182 * t34) * m(7) + (-pkin(4) * t36 - qJ(5) * t33 - qJD(5) * t46 - t35 * t55 - t97 * t66) * m(6); -t138 * qJD(4) + (-mrSges(6,1) * t132 + (t56 + t99) * t87) * qJD(1) + (-qJD(4) * t34 - t68 * t108 + t177) * m(7) + (t46 * qJD(4) + t35 * t128 + t36) * m(6) - t173; -Ifges(7,3) * t118 - t34 * (mrSges(7,1) * t54 + mrSges(7,2) * t53) + (Ifges(7,1) * t53 - t154) * t168 + t13 * t167 + (Ifges(7,5) * t53 - Ifges(7,6) * t54) * t166 - t7 * t31 + t8 * t32 + (t53 * t7 + t54 * t8) * mrSges(7,3) + (-Ifges(7,2) * t54 + t14 + t51) * t170 + t176;];
tauc  = t12(:);
