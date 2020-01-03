% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:21
% DurationCPUTime: 2.02s
% Computational Cost: add. (2630->266), mult. (7084->394), div. (0->0), fcn. (4758->6), ass. (0->136)
t128 = qJD(2) + qJD(3);
t127 = qJD(4) + t128;
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t130 = sin(qJ(3));
t131 = sin(qJ(2));
t133 = cos(qJ(3));
t134 = cos(qJ(2));
t110 = -t130 * t131 + t133 * t134;
t98 = t110 * qJD(1);
t111 = t130 * t134 + t133 * t131;
t99 = t111 * qJD(1);
t138 = -t129 * t99 + t132 * t98;
t77 = t128 * t110;
t67 = t77 * qJD(1);
t78 = t128 * t111;
t68 = t78 * qJD(1);
t13 = qJD(4) * t138 - t129 * t68 + t132 * t67;
t61 = t129 * t98 + t132 * t99;
t14 = -qJD(4) * t61 - t129 * t67 - t132 * t68;
t163 = Ifges(5,4) * t61;
t167 = pkin(7) * t98;
t170 = -pkin(6) - pkin(5);
t123 = t170 * t134;
t116 = qJD(1) * t123;
t103 = t133 * t116;
t122 = t170 * t131;
t115 = qJD(1) * t122;
t106 = qJD(2) * pkin(2) + t115;
t72 = t106 * t130 - t103;
t46 = t72 + t167;
t158 = t129 * t46;
t100 = t130 * t116;
t71 = t133 * t106 + t100;
t93 = t99 * pkin(7);
t45 = t71 - t93;
t41 = pkin(3) * t128 + t45;
t15 = t132 * t41 - t158;
t141 = qJD(2) * t170;
t137 = qJD(1) * t141;
t107 = t131 * t137;
t108 = t134 * t137;
t145 = qJD(3) * t133;
t146 = qJD(3) * t130;
t35 = t106 * t145 + t133 * t107 + t130 * t108 + t116 * t146;
t17 = -pkin(7) * t68 + t35;
t36 = -qJD(3) * t72 - t107 * t130 + t133 * t108;
t18 = -pkin(7) * t67 + t36;
t2 = qJD(4) * t15 + t129 * t18 + t132 * t17;
t55 = Ifges(5,4) * t138;
t29 = Ifges(5,1) * t61 + Ifges(5,5) * t127 + t55;
t156 = t132 * t46;
t16 = t129 * t41 + t156;
t3 = -qJD(4) * t16 - t129 * t17 + t132 * t18;
t125 = -pkin(2) * t134 - pkin(1);
t121 = qJD(1) * t125;
t79 = -t98 * pkin(3) + t121;
t186 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t13 + Ifges(5,6) * t14 - (Ifges(5,5) * t138 - Ifges(5,6) * t61) * t127 / 0.2e1 - (-Ifges(5,2) * t61 + t29 + t55) * t138 / 0.2e1 - t79 * (mrSges(5,1) * t61 + mrSges(5,2) * t138) - (Ifges(5,1) * t138 - t163) * t61 / 0.2e1;
t162 = t16 * t61;
t28 = Ifges(5,2) * t138 + Ifges(5,6) * t127 + t163;
t184 = t28 / 0.2e1;
t124 = pkin(2) * t133 + pkin(3);
t143 = qJD(4) * t132;
t144 = qJD(4) * t129;
t151 = t129 * t130;
t75 = -t115 * t130 + t103;
t47 = t75 - t167;
t76 = t133 * t115 + t100;
t48 = -t93 + t76;
t183 = -t129 * t47 - t132 * t48 + t124 * t143 + (-t130 * t144 + (t132 * t133 - t151) * qJD(3)) * pkin(2);
t150 = t130 * t132;
t182 = t129 * t48 - t132 * t47 - t124 * t144 + (-t130 * t143 + (-t129 * t133 - t150) * qJD(3)) * pkin(2);
t178 = t138 * t15;
t81 = t130 * t122 - t133 * t123;
t176 = t138 / 0.2e1;
t174 = t61 / 0.2e1;
t172 = t98 / 0.2e1;
t171 = t99 / 0.2e1;
t169 = pkin(1) * mrSges(3,1);
t168 = pkin(1) * mrSges(3,2);
t165 = m(4) * t121;
t164 = mrSges(4,3) * t98;
t161 = t99 * mrSges(4,3);
t160 = t99 * Ifges(4,4);
t159 = Ifges(3,4) * t131;
t155 = Ifges(3,5) * qJD(2);
t154 = Ifges(3,6) * qJD(2);
t153 = qJD(2) * mrSges(3,1);
t152 = qJD(2) * mrSges(3,2);
t149 = qJD(1) * t131;
t148 = qJD(1) * t134;
t147 = qJD(2) * t131;
t142 = pkin(2) * t147;
t140 = t155 / 0.2e1;
t139 = -t154 / 0.2e1;
t80 = t133 * t122 + t123 * t130;
t56 = -pkin(7) * t111 + t80;
t57 = pkin(7) * t110 + t81;
t30 = -t129 * t57 + t132 * t56;
t31 = t129 * t56 + t132 * t57;
t73 = t110 * t132 - t111 * t129;
t74 = t110 * t129 + t111 * t132;
t117 = t131 * t141;
t118 = t134 * t141;
t37 = t133 * t117 + t130 * t118 + t122 * t145 + t123 * t146;
t38 = -qJD(3) * t81 - t117 * t130 + t133 * t118;
t53 = t98 * Ifges(4,2) + t128 * Ifges(4,6) + t160;
t92 = Ifges(4,4) * t98;
t54 = t99 * Ifges(4,1) + t128 * Ifges(4,5) + t92;
t135 = mrSges(5,3) * t178 - t35 * mrSges(4,2) - t121 * (mrSges(4,1) * t99 + mrSges(4,2) * t98) + t61 * t184 + t36 * mrSges(4,1) + t53 * t171 - t99 * (Ifges(4,1) * t98 - t160) / 0.2e1 + t71 * t164 - t128 * (Ifges(4,5) * t98 - Ifges(4,6) * t99) / 0.2e1 - Ifges(4,6) * t68 + Ifges(4,5) * t67 - (-Ifges(4,2) * t99 + t54 + t92) * t98 / 0.2e1 + t186;
t126 = Ifges(3,4) * t148;
t120 = mrSges(3,3) * t148 - t152;
t119 = -mrSges(3,3) * t149 + t153;
t97 = Ifges(3,1) * t149 + t126 + t155;
t96 = t154 + (t134 * Ifges(3,2) + t159) * qJD(1);
t95 = pkin(2) * t150 + t124 * t129;
t94 = -pkin(2) * t151 + t124 * t132;
t85 = -t110 * pkin(3) + t125;
t84 = mrSges(4,1) * t128 - t161;
t83 = -mrSges(4,2) * t128 + t164;
t82 = pkin(2) * t149 + pkin(3) * t99;
t70 = -mrSges(4,1) * t98 + mrSges(4,2) * t99;
t62 = pkin(3) * t78 + t142;
t51 = pkin(3) * t68 + qJD(1) * t142;
t50 = mrSges(5,1) * t127 - mrSges(5,3) * t61;
t49 = -mrSges(5,2) * t127 + mrSges(5,3) * t138;
t33 = -mrSges(5,1) * t138 + mrSges(5,2) * t61;
t26 = -pkin(7) * t77 + t38;
t25 = -pkin(7) * t78 + t37;
t24 = -qJD(4) * t74 - t129 * t77 - t132 * t78;
t23 = qJD(4) * t73 - t129 * t78 + t132 * t77;
t20 = t132 * t45 - t158;
t19 = -t129 * t45 - t156;
t5 = -qJD(4) * t31 - t129 * t25 + t132 * t26;
t4 = qJD(4) * t30 + t129 * t26 + t132 * t25;
t1 = [m(5) * (t15 * t5 + t16 * t4 + t2 * t31 + t3 * t30 + t51 * t85 + t62 * t79) + t125 * (mrSges(4,1) * t68 + mrSges(4,2) * t67) + (t110 * t35 - t111 * t36 - t67 * t80 - t68 * t81 - t71 * t77 - t72 * t78) * mrSges(4,3) + (-t110 * t68 - t172 * t78) * Ifges(4,2) + (t110 * t67 - t68 * t111 - t171 * t78 + t172 * t77) * Ifges(4,4) + t127 * (Ifges(5,5) * t23 + Ifges(5,6) * t24) / 0.2e1 - t78 * t53 / 0.2e1 + t79 * (-mrSges(5,1) * t24 + mrSges(5,2) * t23) + t37 * t83 + t38 * t84 + t85 * (-mrSges(5,1) * t14 + mrSges(5,2) * t13) + t51 * (-mrSges(5,1) * t73 + mrSges(5,2) * t74) + t77 * t54 / 0.2e1 + t62 * t33 + t4 * t49 + t5 * t50 + t23 * t29 / 0.2e1 + m(4) * (t35 * t81 + t36 * t80 + t37 * t72 + t38 * t71) + (t97 / 0.2e1 - pkin(5) * t119 + t140 + (-0.2e1 * t168 + 0.3e1 / 0.2e1 * Ifges(3,4) * t134) * qJD(1)) * t134 * qJD(2) + (-t13 * t30 + t14 * t31 - t15 * t23 + t16 * t24 + t2 * t73 - t3 * t74) * mrSges(5,3) + (t73 * t14 + t176 * t24) * Ifges(5,2) + (t73 * t13 + t14 * t74 + t174 * t24 + t176 * t23) * Ifges(5,4) + t128 * (Ifges(4,5) * t77 - Ifges(4,6) * t78) / 0.2e1 + t121 * (mrSges(4,1) * t78 + mrSges(4,2) * t77) + (-t96 / 0.2e1 - pkin(5) * t120 + t139 + (-0.2e1 * t169 - 0.3e1 / 0.2e1 * t159 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t134) * qJD(1) + (0.2e1 * t165 + t70 + qJD(1) * (-mrSges(4,1) * t110 + mrSges(4,2) * t111)) * pkin(2)) * t147 + (t67 * t111 + t171 * t77) * Ifges(4,1) + (t13 * t74 + t174 * t23) * Ifges(5,1) + t24 * t184; t182 * t50 + t183 * t49 + (-t13 * t94 + t14 * t95 + t162) * mrSges(5,3) - m(4) * (t71 * t75 + t72 * t76) + t72 * t161 - t82 * t33 - t76 * t83 - t75 * t84 + t135 + ((t140 - t97 / 0.2e1 - t126 / 0.2e1 + qJD(1) * t168 + (t119 - t153) * pkin(5)) * t134 + (t139 + t96 / 0.2e1 + (t169 + t159 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t134) * qJD(1) + (t120 + t152) * pkin(5) + (-t70 - t165) * pkin(2)) * t131) * qJD(1) + (m(4) * (t130 * t35 + t133 * t36 + t145 * t72 - t146 * t71) + t83 * t145 - t84 * t146 + (-t130 * t68 - t133 * t67) * mrSges(4,3)) * pkin(2) + (t15 * t182 + t16 * t183 + t2 * t95 + t3 * t94 - t79 * t82) * m(5); (t84 + t161) * t72 + mrSges(5,3) * t162 - t71 * t83 - t20 * t49 - t19 * t50 + t135 + (-t99 * t33 + (-t129 * t50 + t132 * t49) * qJD(4) + (t129 * t14 - t13 * t132) * mrSges(5,3) + (t129 * t2 + t132 * t3 + t143 * t16 - t144 * t15 - t79 * t99) * m(5)) * pkin(3) - m(5) * (t15 * t19 + t16 * t20); t28 * t174 - t15 * t49 + t16 * t50 + (t178 + t162) * mrSges(5,3) + t186;];
tauc = t1(:);
