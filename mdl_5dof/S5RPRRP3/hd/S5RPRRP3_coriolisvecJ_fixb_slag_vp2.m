% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:18
% DurationCPUTime: 2.36s
% Computational Cost: add. (1953->246), mult. (4844->338), div. (0->0), fcn. (2972->6), ass. (0->118)
t144 = Ifges(5,4) + Ifges(6,4);
t161 = Ifges(5,1) + Ifges(6,1);
t167 = Ifges(5,5) + Ifges(6,5);
t160 = Ifges(5,2) + Ifges(6,2);
t166 = Ifges(5,6) + Ifges(6,6);
t113 = sin(qJ(4));
t114 = sin(qJ(3));
t115 = cos(qJ(4));
t116 = cos(qJ(3));
t98 = -t113 * t114 + t115 * t116;
t93 = t98 * qJD(1);
t169 = t144 * t93;
t99 = t113 * t116 + t115 * t114;
t94 = t99 * qJD(1);
t168 = t144 * t94;
t110 = qJD(3) + qJD(4);
t165 = t166 * t110 + t160 * t93 + t168;
t164 = t167 * t110 + t161 * t94 + t169;
t163 = -Ifges(4,1) / 0.2e1;
t133 = qJD(1) * t116;
t108 = Ifges(4,4) * t133;
t162 = -t108 / 0.2e1;
t105 = sin(pkin(8)) * pkin(1) + pkin(6);
t101 = t105 * qJD(1);
t122 = pkin(7) * qJD(1) + t101;
t129 = t114 * qJD(2);
t69 = t116 * t122 + t129;
t63 = t113 * t69;
t109 = t116 * qJD(2);
t119 = t122 * t114;
t68 = t109 - t119;
t66 = qJD(3) * pkin(3) + t68;
t18 = t115 * t66 - t63;
t86 = t94 * qJ(5);
t13 = t18 - t86;
t157 = -t93 / 0.2e1;
t154 = t94 / 0.2e1;
t127 = -cos(pkin(8)) * pkin(1) - pkin(2);
t100 = -pkin(3) * t116 + t127;
t95 = qJD(1) * t100;
t153 = m(5) * t95;
t152 = pkin(4) * t94;
t149 = mrSges(5,3) * t93;
t148 = mrSges(6,3) * t93;
t147 = t94 * mrSges(5,3);
t143 = pkin(7) + t105;
t25 = t115 * t68 - t63;
t71 = -mrSges(6,2) * t110 + t148;
t72 = -mrSges(5,2) * t110 + t149;
t142 = t71 + t72;
t73 = mrSges(6,1) * t110 - t94 * mrSges(6,3);
t74 = mrSges(5,1) * t110 - t147;
t141 = t73 + t74;
t96 = t143 * t114;
t97 = t143 * t116;
t46 = -t113 * t96 + t115 * t97;
t140 = Ifges(4,4) * t114;
t139 = qJ(5) * t93;
t103 = t127 * qJD(1);
t138 = t103 * mrSges(4,2);
t60 = t110 * t99;
t48 = t60 * qJD(1);
t137 = t113 * t48;
t65 = t115 * t69;
t136 = Ifges(4,5) * qJD(3);
t135 = Ifges(4,6) * qJD(3);
t134 = qJD(1) * t114;
t132 = qJD(3) * t114;
t131 = qJD(4) * t113;
t130 = qJD(4) * t115;
t128 = pkin(3) * t132;
t126 = t136 / 0.2e1;
t125 = -t135 / 0.2e1;
t124 = m(4) * t105 + mrSges(4,3);
t24 = -t113 * t68 - t65;
t45 = -t113 * t97 - t115 * t96;
t123 = qJD(3) * t143;
t81 = -t101 * t114 + t109;
t121 = t81 * mrSges(4,3) + t134 * t163 + t162 - t136 / 0.2e1;
t19 = t113 * t66 + t65;
t82 = t101 * t116 + t129;
t118 = t82 * mrSges(4,3) + t135 / 0.2e1 + (t116 * Ifges(4,2) + t140) * qJD(1) / 0.2e1 - t103 * mrSges(4,1);
t106 = qJD(3) * t109;
t61 = -qJD(3) * t119 + t106;
t62 = t69 * qJD(3);
t7 = -t113 * t62 + t115 * t61 + t66 * t130 - t131 * t69;
t87 = t114 * t123;
t88 = t116 * t123;
t11 = -t113 * t88 - t115 * t87 - t96 * t130 - t131 * t97;
t8 = -qJD(4) * t19 - t113 * t61 - t115 * t62;
t12 = -qJD(4) * t46 + t113 * t87 - t115 * t88;
t59 = t110 * t98;
t10 = pkin(4) * t110 + t13;
t2 = -qJ(5) * t48 + qJD(5) * t93 + t7;
t47 = t59 * qJD(1);
t3 = -qJ(5) * t47 - qJD(5) * t94 + t8;
t53 = -t93 * pkin(4) + qJD(5) + t95;
t117 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t7 * mrSges(5,2) - t2 * mrSges(6,2) + t10 * t148 + t18 * t149 - t53 * (mrSges(6,1) * t94 + mrSges(6,2) * t93) - t95 * (mrSges(5,1) * t94 + mrSges(5,2) * t93) - t166 * t48 + t167 * t47 - (t161 * t93 - t168) * t94 / 0.2e1 + t165 * t154 - (-t166 * t94 + t167 * t93) * t110 / 0.2e1 + (-t160 * t94 + t164 + t169) * t157;
t107 = pkin(3) * t115 + pkin(4);
t104 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t133;
t102 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t134;
t76 = t82 * qJD(3);
t75 = -t101 * t132 + t106;
t70 = pkin(3) * t134 + t152;
t67 = -t98 * pkin(4) + t100;
t52 = -mrSges(5,1) * t93 + mrSges(5,2) * t94;
t51 = -mrSges(6,1) * t93 + mrSges(6,2) * t94;
t40 = t47 * mrSges(6,2);
t39 = pkin(4) * t60 + t128;
t32 = pkin(4) * t48 + qJD(1) * t128;
t31 = qJ(5) * t98 + t46;
t30 = -qJ(5) * t99 + t45;
t17 = -t86 + t25;
t16 = t24 - t139;
t14 = t19 + t139;
t5 = -qJ(5) * t59 - qJD(5) * t99 + t12;
t4 = -qJ(5) * t60 + qJD(5) * t98 + t11;
t1 = [m(6) * (t10 * t5 + t14 * t4 + t2 * t31 + t3 * t30 + t32 * t67 + t39 * t53) + (t124 * t76 + (t125 + (-m(4) * t82 - t104) * t105 + (t127 * mrSges(4,1) - 0.3e1 / 0.2e1 * t140 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t116) * qJD(1) + (0.2e1 * t153 + t52 + qJD(1) * (-mrSges(5,1) * t98 + mrSges(5,2) * t99)) * pkin(3) - t118) * qJD(3)) * t114 + (t124 * t75 + (0.3e1 / 0.2e1 * t108 + t126 + (-m(4) * t81 - t102) * t105 + 0.2e1 * t138 - t121) * qJD(3)) * t116 + (mrSges(5,2) * t100 - mrSges(5,3) * t45 - mrSges(6,3) * t30 + t144 * t98 + t161 * t99) * t47 - (-t100 * mrSges(5,1) - t67 * mrSges(6,1) + t46 * mrSges(5,3) + t31 * mrSges(6,3) + t144 * t99 + t160 * t98) * t48 + (-t10 * t59 - t14 * t60 + t2 * t98 - t3 * t99) * mrSges(6,3) + (-t18 * t59 - t19 * t60 + t7 * t98 - t8 * t99) * mrSges(5,3) + t95 * (mrSges(5,1) * t60 + mrSges(5,2) * t59) + t32 * (-mrSges(6,1) * t98 + mrSges(6,2) * t99) + t4 * t71 + t11 * t72 + t5 * t73 + t12 * t74 + m(5) * (t11 * t19 + t12 * t18 + t45 * t8 + t46 * t7) + t53 * (mrSges(6,1) * t60 + mrSges(6,2) * t59) + t39 * t51 + t67 * t40 + t164 * t59 / 0.2e1 - t165 * t60 / 0.2e1 + (t144 * t59 - t160 * t60) * t93 / 0.2e1 + (-t144 * t60 + t161 * t59) * t154 + (-t166 * t60 + t167 * t59) * t110 / 0.2e1; -t141 * t60 + t142 * t59 + (-t114 * t102 + t116 * t104 + (-t114 ^ 2 - t116 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + m(4) * (t75 * t114 - t116 * t76 + (-t114 * t81 + t116 * t82) * qJD(3)) + m(5) * (-t18 * t60 + t19 * t59 + t7 * t99 + t8 * t98) + m(6) * (-t10 * t60 + t14 * t59 + t2 * t99 + t3 * t98) + (mrSges(6,3) + mrSges(5,3)) * (-t98 * t47 - t48 * t99); (-t107 * t47 + t14 * t94) * mrSges(6,3) - m(5) * (t18 * t24 + t19 * t25) + t117 + t19 * t147 + t82 * t102 - t81 * t104 - t70 * t51 - t17 * t71 - t25 * t72 - t16 * t73 - t24 * t74 - t75 * mrSges(4,2) - t76 * mrSges(4,1) + ((t162 - t138 + t126 + t121) * t116 + (t125 + (t140 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t163) * t116) * qJD(1) + (-t52 - t153) * pkin(3) + t118) * t114) * qJD(1) + (-mrSges(6,3) * t137 + (-t115 * t47 - t137) * mrSges(5,3) + (-t113 * t141 + t115 * t142) * qJD(4) + m(5) * (t113 * t7 + t115 * t8 + t130 * t19 - t131 * t18)) * pkin(3) + (t3 * t107 - t10 * t16 - t14 * t17 - t53 * t70 + (-t10 * t131 + t113 * t2 + t130 * t14) * pkin(3)) * m(6); t117 + (t19 * mrSges(5,3) + t14 * mrSges(6,3) - pkin(4) * t51) * t94 - t13 * t71 - t18 * t72 + t14 * t73 + t19 * t74 - pkin(4) * t47 * mrSges(6,3) + (-t53 * t152 - (-t10 + t13) * t14 + t3 * pkin(4)) * m(6); t48 * mrSges(6,1) - t93 * t71 + t94 * t73 + t40 + 0.2e1 * (t32 / 0.2e1 + t10 * t154 + t14 * t157) * m(6);];
tauc = t1(:);
