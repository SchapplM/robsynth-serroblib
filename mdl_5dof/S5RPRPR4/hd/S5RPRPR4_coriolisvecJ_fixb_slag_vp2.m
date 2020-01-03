% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:15
% EndTime: 2020-01-03 11:38:26
% DurationCPUTime: 2.54s
% Computational Cost: add. (2586->285), mult. (6429->407), div. (0->0), fcn. (4341->8), ass. (0->141)
t119 = qJD(3) + qJD(5);
t124 = sin(qJ(5));
t126 = cos(qJ(5));
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t125 = sin(qJ(3));
t127 = cos(qJ(3));
t105 = -t120 * t125 + t122 * t127;
t97 = t105 * qJD(1);
t146 = qJD(1) * t127;
t147 = qJD(1) * t125;
t98 = -t120 * t146 - t122 * t147;
t134 = t124 * t98 + t126 * t97;
t42 = Ifges(6,4) * t134;
t48 = t124 * t97 - t126 * t98;
t15 = Ifges(6,1) * t48 + Ifges(6,5) * t119 + t42;
t161 = Ifges(6,4) * t48;
t106 = t120 * t127 + t122 * t125;
t99 = t106 * qJD(3);
t86 = qJD(1) * t99;
t100 = t105 * qJD(3);
t87 = qJD(1) * t100;
t19 = t134 * qJD(5) - t124 * t86 + t126 * t87;
t144 = qJD(4) * t127;
t145 = qJD(3) * t125;
t114 = sin(pkin(8)) * pkin(1) + pkin(6);
t108 = t114 * qJD(1);
t118 = t127 * qJD(2);
t80 = qJD(3) * t118 - t108 * t145;
t55 = (-qJ(4) * t145 + t144) * qJD(1) + t80;
t142 = t125 * qJD(4);
t132 = qJ(4) * qJD(1) + t108;
t143 = t125 * qJD(2);
t73 = t127 * t132 + t143;
t56 = -qJD(1) * t142 - t73 * qJD(3);
t24 = -t120 * t55 + t122 * t56;
t12 = -pkin(7) * t87 + t24;
t25 = t120 * t56 + t122 * t55;
t13 = -pkin(7) * t86 + t25;
t165 = pkin(7) * t98;
t63 = t120 * t73;
t72 = -t125 * t132 + t118;
t67 = qJD(3) * pkin(3) + t72;
t32 = t122 * t67 - t63;
t22 = qJD(3) * pkin(4) + t165 + t32;
t166 = pkin(7) * t97;
t154 = t122 * t73;
t33 = t120 * t67 + t154;
t23 = t33 + t166;
t6 = -t124 * t23 + t126 * t22;
t2 = qJD(5) * t6 + t12 * t124 + t126 * t13;
t20 = -qJD(5) * t48 - t124 * t87 - t126 * t86;
t7 = t124 * t22 + t126 * t23;
t3 = -qJD(5) * t7 + t12 * t126 - t124 * t13;
t140 = -cos(pkin(8)) * pkin(1) - pkin(2);
t107 = -pkin(3) * t127 + t140;
t148 = qJD(1) * t107;
t96 = qJD(4) + t148;
t57 = -t97 * pkin(4) + t96;
t185 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 - (Ifges(6,5) * t134 - Ifges(6,6) * t48) * t119 / 0.2e1 - (-Ifges(6,2) * t48 + t15 + t42) * t134 / 0.2e1 - t57 * (mrSges(6,1) * t48 + mrSges(6,2) * t134) - (Ifges(6,1) * t134 - t161) * t48 / 0.2e1;
t184 = t134 * t6 + t48 * t7;
t182 = -Ifges(4,1) / 0.2e1;
t14 = Ifges(6,2) * t134 + Ifges(6,6) * t119 + t161;
t181 = t14 / 0.2e1;
t117 = Ifges(4,4) * t146;
t180 = -t117 / 0.2e1;
t34 = -t120 * t72 - t154;
t28 = t34 - t166;
t35 = t122 * t72 - t63;
t29 = t35 + t165;
t115 = pkin(3) * t122 + pkin(4);
t160 = pkin(3) * t120;
t94 = t115 * t126 - t124 * t160;
t179 = t94 * qJD(5) - t124 * t28 - t126 * t29;
t95 = t115 * t124 + t126 * t160;
t178 = -t95 * qJD(5) + t124 * t29 - t126 * t28;
t172 = t134 / 0.2e1;
t170 = t48 / 0.2e1;
t168 = -t98 / 0.2e1;
t167 = -t99 / 0.2e1;
t164 = t100 / 0.2e1;
t162 = Ifges(5,4) * t98;
t59 = t105 * t124 + t106 * t126;
t159 = t20 * t59;
t58 = t105 * t126 - t106 * t124;
t158 = t58 * t19;
t149 = qJ(4) + t114;
t133 = qJD(3) * t149;
t76 = -t125 * t133 + t144;
t77 = -t127 * t133 - t142;
t37 = t120 * t77 + t122 * t76;
t103 = t149 * t125;
t104 = t149 * t127;
t53 = -t120 * t103 + t122 * t104;
t157 = Ifges(4,4) * t125;
t156 = t105 * t87;
t110 = t140 * qJD(1);
t155 = t110 * mrSges(4,2);
t152 = t86 * t106;
t151 = Ifges(4,5) * qJD(3);
t150 = Ifges(4,6) * qJD(3);
t141 = pkin(3) * t145;
t139 = t151 / 0.2e1;
t138 = -t150 / 0.2e1;
t137 = m(4) * t114 + mrSges(4,3);
t136 = t86 * mrSges(5,1) + t87 * mrSges(5,2);
t135 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t36 = -t120 * t76 + t122 * t77;
t52 = -t122 * t103 - t104 * t120;
t88 = -t108 * t125 + t118;
t131 = t88 * mrSges(4,3) + t147 * t182 + t180 - t151 / 0.2e1;
t130 = qJD(1) * t141;
t38 = -pkin(7) * t106 + t52;
t39 = pkin(7) * t105 + t53;
t10 = -t124 * t39 + t126 * t38;
t11 = t124 * t38 + t126 * t39;
t89 = t108 * t127 + t143;
t128 = t89 * mrSges(4,3) + t150 / 0.2e1 + (t127 * Ifges(4,2) + t157) * qJD(1) / 0.2e1 - t110 * mrSges(4,1);
t111 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t146;
t109 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t147;
t93 = Ifges(5,4) * t97;
t81 = t89 * qJD(3);
t79 = qJD(3) * mrSges(5,1) + t98 * mrSges(5,3);
t78 = -qJD(3) * mrSges(5,2) + t97 * mrSges(5,3);
t75 = pkin(4) * t99 + t141;
t74 = pkin(3) * t147 - pkin(4) * t98;
t71 = -t105 * pkin(4) + t107;
t62 = pkin(4) * t86 + t130;
t54 = -mrSges(5,1) * t97 - mrSges(5,2) * t98;
t44 = -t98 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t93;
t43 = t97 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t162;
t41 = mrSges(6,1) * t119 - mrSges(6,3) * t48;
t40 = -mrSges(6,2) * t119 + mrSges(6,3) * t134;
t31 = -pkin(7) * t99 + t37;
t30 = -pkin(7) * t100 + t36;
t27 = -qJD(5) * t59 - t100 * t124 - t126 * t99;
t26 = qJD(5) * t58 + t100 * t126 - t124 * t99;
t21 = -mrSges(6,1) * t134 + mrSges(6,2) * t48;
t5 = -qJD(5) * t11 - t124 * t31 + t126 * t30;
t4 = qJD(5) * t10 + t124 * t30 + t126 * t31;
t1 = [(-t10 * t19 + t11 * t20 + t2 * t58 - t26 * t6 + t27 * t7 - t3 * t59) * mrSges(6,3) + (t137 * t81 + (t138 + (-m(4) * t89 - t111) * t114 + (t140 * mrSges(4,1) - 0.3e1 / 0.2e1 * t157 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t127) * qJD(1) + (m(5) * (t96 + t148) + t54 + qJD(1) * (-mrSges(5,1) * t105 + mrSges(5,2) * t106)) * pkin(3) - t128) * qJD(3)) * t125 + t96 * (mrSges(5,1) * t99 + mrSges(5,2) * t100) + qJD(3) * (Ifges(5,5) * t100 - Ifges(5,6) * t99) / 0.2e1 + t119 * (Ifges(6,5) * t26 + Ifges(6,6) * t27) / 0.2e1 + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + t57 * t75 + t62 * t71) + t37 * t78 + t36 * t79 + t62 * (-mrSges(6,1) * t58 + mrSges(6,2) * t59) + t75 * t21 + t57 * (-mrSges(6,1) * t27 + mrSges(6,2) * t26) + t4 * t40 + t5 * t41 + t26 * t15 / 0.2e1 + m(5) * (t24 * t52 + t25 * t53 + t32 * t36 + t33 * t37) + t44 * t164 + t43 * t167 + (-t100 * t32 + t105 * t25 - t106 * t24 - t33 * t99 - t52 * t87 - t53 * t86) * mrSges(5,3) + (-t105 * t86 + t167 * t97) * Ifges(5,2) + (t164 * t97 - t168 * t99 - t152 + t156) * Ifges(5,4) + t71 * t135 + t107 * t136 + (t137 * t80 + (0.3e1 / 0.2e1 * t117 + t139 + (-m(4) * t88 - t109) * t114 + 0.2e1 * t155 - t131) * qJD(3)) * t127 + (t100 * t168 + t87 * t106) * Ifges(5,1) + (t170 * t26 + t19 * t59) * Ifges(6,1) + (t172 * t27 + t58 * t20) * Ifges(6,2) + (t170 * t27 + t172 * t26 + t158 + t159) * Ifges(6,4) + t27 * t181; t100 * t78 + t26 * t40 + t27 * t41 - t99 * t79 + (-t158 + t159) * mrSges(6,3) + (-t152 - t156) * mrSges(5,3) + (-t125 * t109 + t127 * t111 + (-t125 ^ 2 - t127 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + m(4) * (t80 * t125 - t127 * t81 + (-t125 * t88 + t127 * t89) * qJD(3)) + m(5) * (t100 * t33 + t105 * t24 + t106 * t25 - t32 * t99) + m(6) * (t2 * t59 + t26 * t7 + t27 * t6 + t3 * t58); (-t19 * t94 + t20 * t95 + t184) * mrSges(6,3) - (Ifges(5,2) * t98 + t44 + t93) * t97 / 0.2e1 + t89 * t109 - t88 * t111 - t96 * (-mrSges(5,1) * t98 + mrSges(5,2) * t97) - qJD(3) * (Ifges(5,5) * t97 + Ifges(5,6) * t98) / 0.2e1 - t34 * t79 - t80 * mrSges(4,2) - t81 * mrSges(4,1) - Ifges(5,6) * t86 + Ifges(5,5) * t87 - t74 * t21 - t35 * t78 + t24 * mrSges(5,1) - t25 * mrSges(5,2) + m(5) * (t120 * t25 + t122 * t24) * pkin(3) + t185 + t43 * t168 - m(5) * (t32 * t34 + t33 * t35) + ((t180 - t155 + t139 + t131) * t127 + (t138 + (t157 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t182) * t127) * qJD(1) + (-m(5) * t96 - t54) * pkin(3) + t128) * t125) * qJD(1) + t178 * t41 + t179 * t40 + (t178 * t6 + t179 * t7 + t2 * t95 + t3 * t94 - t57 * t74) * m(6) + t48 * t181 + (t32 * t97 - t33 * t98 + (-t120 * t86 - t122 * t87) * pkin(3)) * mrSges(5,3) + t98 * (Ifges(5,1) * t97 + t162) / 0.2e1; -t134 * t40 + t48 * t41 - t97 * t78 - t98 * t79 + t135 + t136 + (-t134 * t7 + t48 * t6 + t62) * m(6) + (-t32 * t98 - t33 * t97 + t130) * m(5); t184 * mrSges(6,3) + t14 * t170 - t6 * t40 + t7 * t41 + t185;];
tauc = t1(:);
