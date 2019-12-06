% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:10
% DurationCPUTime: 2.32s
% Computational Cost: add. (2226->281), mult. (5784->406), div. (0->0), fcn. (3981->6), ass. (0->135)
t116 = qJD(3) + qJD(5);
t119 = sin(qJ(5));
t121 = cos(qJ(5));
t117 = sin(pkin(9));
t118 = cos(pkin(9));
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t97 = -t117 * t120 + t118 * t122;
t89 = t97 * qJD(2);
t138 = qJD(2) * t122;
t139 = qJD(2) * t120;
t90 = -t117 * t138 - t118 * t139;
t128 = t119 * t90 + t121 * t89;
t42 = Ifges(6,4) * t128;
t48 = t119 * t89 - t121 * t90;
t13 = Ifges(6,1) * t48 + Ifges(6,5) * t116 + t42;
t152 = Ifges(6,4) * t48;
t98 = t117 * t122 + t118 * t120;
t91 = t98 * qJD(3);
t80 = qJD(2) * t91;
t92 = t97 * qJD(3);
t81 = qJD(2) * t92;
t17 = qJD(5) * t128 - t119 * t80 + t121 * t81;
t18 = -qJD(5) * t48 - t119 * t81 - t121 * t80;
t135 = qJD(1) * qJD(3);
t112 = t122 * t135;
t146 = -qJ(4) - pkin(6);
t127 = qJD(3) * t146;
t87 = qJD(4) * t122 + t120 * t127;
t57 = qJD(2) * t87 + t112;
t88 = -t120 * qJD(4) + t122 * t127;
t58 = qJD(2) * t88 - t120 * t135;
t26 = -t117 * t57 + t118 * t58;
t20 = -pkin(7) * t81 + t26;
t27 = t117 * t58 + t118 * t57;
t21 = -pkin(7) * t80 + t27;
t155 = pkin(7) * t90;
t108 = t146 * t122;
t136 = t120 * qJD(1);
t86 = -qJD(2) * t108 + t136;
t68 = t117 * t86;
t107 = t146 * t120;
t115 = t122 * qJD(1);
t85 = qJD(2) * t107 + t115;
t78 = qJD(3) * pkin(3) + t85;
t34 = t118 * t78 - t68;
t24 = qJD(3) * pkin(4) + t155 + t34;
t156 = pkin(7) * t89;
t144 = t118 * t86;
t35 = t117 * t78 + t144;
t25 = t35 + t156;
t6 = -t119 * t25 + t121 * t24;
t2 = qJD(5) * t6 + t119 * t20 + t121 * t21;
t7 = t119 * t24 + t121 * t25;
t3 = -qJD(5) * t7 - t119 * t21 + t121 * t20;
t113 = -pkin(3) * t122 - pkin(2);
t140 = qJD(2) * t113;
t104 = qJD(4) + t140;
t55 = -t89 * pkin(4) + t104;
t177 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 - (Ifges(6,5) * t128 - Ifges(6,6) * t48) * t116 / 0.2e1 - (-Ifges(6,2) * t48 + t13 + t42) * t128 / 0.2e1 - t55 * (mrSges(6,1) * t48 + mrSges(6,2) * t128) - (Ifges(6,1) * t128 - t152) * t48 / 0.2e1;
t176 = qJD(3) / 0.2e1;
t175 = t128 * t6 + t48 * t7;
t12 = Ifges(6,2) * t128 + Ifges(6,6) * t116 + t152;
t173 = t12 / 0.2e1;
t133 = Ifges(4,5) * t176;
t36 = -t117 * t85 - t144;
t28 = t36 - t156;
t38 = t118 * t85 - t68;
t30 = t38 + t155;
t111 = pkin(3) * t118 + pkin(4);
t151 = pkin(3) * t117;
t83 = t111 * t121 - t119 * t151;
t172 = t83 * qJD(5) - t119 * t28 - t121 * t30;
t84 = t111 * t119 + t121 * t151;
t171 = -t84 * qJD(5) + t119 * t30 - t121 * t28;
t165 = t128 / 0.2e1;
t163 = t48 / 0.2e1;
t161 = -t90 / 0.2e1;
t160 = -t91 / 0.2e1;
t159 = t92 / 0.2e1;
t158 = pkin(2) * mrSges(4,1);
t157 = pkin(2) * mrSges(4,2);
t153 = Ifges(5,4) * t90;
t51 = t119 * t97 + t121 * t98;
t150 = t18 * t51;
t50 = -t119 * t98 + t121 * t97;
t149 = t50 * t17;
t148 = t80 * t98;
t147 = t97 * t81;
t39 = t117 * t88 + t118 * t87;
t145 = Ifges(4,4) * t120;
t60 = t117 * t107 - t118 * t108;
t141 = Ifges(4,6) * qJD(3);
t137 = qJD(3) * t120;
t134 = m(4) * pkin(6) + mrSges(4,3);
t132 = -t141 / 0.2e1;
t131 = t80 * mrSges(5,1) + t81 * mrSges(5,2);
t130 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t129 = qJD(2) * t137;
t37 = -t117 * t87 + t118 * t88;
t59 = t118 * t107 + t108 * t117;
t102 = -pkin(6) * t139 + t115;
t114 = Ifges(4,4) * t138;
t126 = Ifges(4,1) * t139 / 0.2e1 + t114 / 0.2e1 + t133 - t102 * mrSges(4,3);
t103 = pkin(6) * t138 + t136;
t125 = t103 * mrSges(4,3) + t141 / 0.2e1 + (t122 * Ifges(4,2) + t145) * qJD(2) / 0.2e1;
t124 = pkin(3) * t129;
t40 = -pkin(7) * t98 + t59;
t41 = pkin(7) * t97 + t60;
t10 = -t119 * t41 + t121 * t40;
t11 = t119 * t40 + t121 * t41;
t106 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t138;
t105 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t139;
t96 = t103 * qJD(3);
t95 = -pkin(6) * t129 + t112;
t82 = Ifges(5,4) * t89;
t67 = -t97 * pkin(4) + t113;
t66 = qJD(3) * mrSges(5,1) + t90 * mrSges(5,3);
t65 = -qJD(3) * mrSges(5,2) + t89 * mrSges(5,3);
t62 = pkin(3) * t137 + pkin(4) * t91;
t61 = pkin(3) * t139 - pkin(4) * t90;
t56 = pkin(4) * t80 + t124;
t49 = -mrSges(5,1) * t89 - mrSges(5,2) * t90;
t44 = -Ifges(5,1) * t90 + Ifges(5,5) * qJD(3) + t82;
t43 = Ifges(5,2) * t89 + Ifges(5,6) * qJD(3) - t153;
t33 = mrSges(6,1) * t116 - mrSges(6,3) * t48;
t32 = -mrSges(6,2) * t116 + mrSges(6,3) * t128;
t31 = -pkin(7) * t91 + t39;
t29 = -pkin(7) * t92 + t37;
t23 = -qJD(5) * t51 - t119 * t92 - t121 * t91;
t22 = qJD(5) * t50 - t119 * t91 + t121 * t92;
t19 = -mrSges(6,1) * t128 + mrSges(6,2) * t48;
t5 = -qJD(5) * t11 - t119 * t31 + t121 * t29;
t4 = qJD(5) * t10 + t119 * t29 + t121 * t31;
t1 = [t22 * t32 + t23 * t33 + t92 * t65 - t91 * t66 + (-t149 + t150) * mrSges(6,3) + (-t147 - t148) * mrSges(5,3) + (-t120 * t105 + t122 * t106 + (-t120 ^ 2 - t122 ^ 2) * qJD(2) * mrSges(4,3)) * qJD(3) + m(4) * (t95 * t120 - t122 * t96 + (-t102 * t120 + t103 * t122) * qJD(3)) + m(5) * (t26 * t97 + t27 * t98 - t34 * t91 + t35 * t92) + m(6) * (t2 * t51 + t22 * t7 + t23 * t6 + t3 * t50); m(5) * (t26 * t59 + t27 * t60 + t34 * t37 + t35 * t39) + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + t55 * t62 + t56 * t67) + (-t26 * t98 + t27 * t97 - t34 * t92 - t35 * t91 - t59 * t81 - t60 * t80) * mrSges(5,3) + (t160 * t89 - t97 * t80) * Ifges(5,2) + (t159 * t89 - t161 * t91 + t147 - t148) * Ifges(5,4) + (-t10 * t17 + t11 * t18 + t2 * t50 - t22 * t6 + t23 * t7 - t3 * t51) * mrSges(6,3) + t116 * (Ifges(6,5) * t22 + Ifges(6,6) * t23) / 0.2e1 + t62 * t19 + t39 * t65 + t37 * t66 + t55 * (-mrSges(6,1) * t23 + mrSges(6,2) * t22) + t56 * (-mrSges(6,1) * t50 + mrSges(6,2) * t51) + t4 * t32 + t5 * t33 + t22 * t13 / 0.2e1 + t104 * (mrSges(5,1) * t91 + mrSges(5,2) * t92) + (t134 * t96 + (t132 + (-m(4) * t103 - t106) * pkin(6) + (-0.2e1 * t158 - 0.3e1 / 0.2e1 * t145 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t122) * qJD(2) + (m(5) * (t104 + t140) + t49 + qJD(2) * (-mrSges(5,1) * t97 + mrSges(5,2) * t98)) * pkin(3) - t125) * qJD(3)) * t120 + (t163 * t22 + t17 * t51) * Ifges(6,1) + (t165 * t23 + t50 * t18) * Ifges(6,2) + (t163 * t23 + t165 * t22 + t149 + t150) * Ifges(6,4) + (t134 * t95 + (t133 + (-0.2e1 * t157 + 0.3e1 / 0.2e1 * Ifges(4,4) * t122) * qJD(2) + (-m(4) * t102 - t105) * pkin(6) + t126) * qJD(3)) * t122 + (t161 * t92 + t81 * t98) * Ifges(5,1) + t67 * t130 + t113 * t131 + (Ifges(5,5) * t92 - Ifges(5,6) * t91) * t176 + t23 * t173 + t44 * t159 + t43 * t160; (-t17 * t83 + t18 * t84 + t175) * mrSges(6,3) - m(5) * (t34 * t36 + t35 * t38) + t177 + t48 * t173 + (t34 * t89 - t35 * t90 + (-t117 * t80 - t118 * t81) * pkin(3)) * mrSges(5,3) - t104 * (-mrSges(5,1) * t90 + mrSges(5,2) * t89) + t103 * t105 - t102 * t106 - t95 * mrSges(4,2) - t96 * mrSges(4,1) - qJD(3) * (Ifges(5,5) * t89 + Ifges(5,6) * t90) / 0.2e1 - Ifges(5,6) * t80 + Ifges(5,5) * t81 - t38 * t65 - t36 * t66 - t61 * t19 - t27 * mrSges(5,2) + t26 * mrSges(5,1) - (Ifges(5,2) * t90 + t44 + t82) * t89 / 0.2e1 + t171 * t33 + t172 * t32 + (t171 * t6 + t172 * t7 + t2 * t84 + t3 * t83 - t55 * t61) * m(6) + ((t133 - t114 / 0.2e1 + qJD(2) * t157 - t126) * t122 + (t132 + (t158 + t145 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t122) * qJD(2) + (-m(5) * t104 - t49) * pkin(3) + t125) * t120) * qJD(2) + t90 * (Ifges(5,1) * t89 + t153) / 0.2e1 + t43 * t161 + m(5) * (t117 * t27 + t118 * t26) * pkin(3); -t128 * t32 + t48 * t33 - t89 * t65 - t90 * t66 + t130 + t131 + (-t128 * t7 + t48 * t6 + t56) * m(6) + (-t34 * t90 - t35 * t89 + t124) * m(5); t175 * mrSges(6,3) + t12 * t163 - t6 * t32 + t7 * t33 + t177;];
tauc = t1(:);
