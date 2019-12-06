% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:31
% EndTime: 2019-12-05 16:43:39
% DurationCPUTime: 2.33s
% Computational Cost: add. (1659->244), mult. (4315->339), div. (0->0), fcn. (2678->4), ass. (0->111)
t135 = Ifges(5,4) + Ifges(6,4);
t155 = Ifges(5,1) + Ifges(6,1);
t159 = Ifges(5,5) + Ifges(6,5);
t154 = Ifges(5,2) + Ifges(6,2);
t158 = Ifges(5,6) + Ifges(6,6);
t108 = sin(qJ(4));
t109 = sin(qJ(3));
t110 = cos(qJ(4));
t111 = cos(qJ(3));
t88 = -t108 * t109 + t110 * t111;
t80 = t88 * qJD(2);
t161 = t135 * t80;
t89 = t108 * t111 + t110 * t109;
t81 = t89 * qJD(2);
t160 = t135 * t81;
t107 = qJD(3) + qJD(4);
t157 = t158 * t107 + t154 * t80 + t160;
t156 = t159 * t107 + t155 * t81 + t161;
t118 = Ifges(4,5) * qJD(3) / 0.2e1;
t147 = -pkin(7) - pkin(6);
t101 = t147 * t111;
t122 = t109 * qJD(1);
t77 = -qJD(2) * t101 + t122;
t69 = t108 * t77;
t106 = t111 * qJD(1);
t120 = qJD(2) * t147;
t76 = t109 * t120 + t106;
t72 = qJD(3) * pkin(3) + t76;
t27 = t110 * t72 - t69;
t73 = t81 * qJ(5);
t11 = t27 - t73;
t151 = -t80 / 0.2e1;
t148 = t81 / 0.2e1;
t104 = -pkin(3) * t111 - pkin(2);
t99 = qJD(2) * t104;
t146 = m(5) * t99;
t145 = pkin(2) * mrSges(4,1);
t144 = pkin(2) * mrSges(4,2);
t143 = pkin(4) * t81;
t140 = mrSges(5,3) * t80;
t139 = mrSges(6,3) * t80;
t138 = t81 * mrSges(5,3);
t61 = -mrSges(6,2) * t107 + t139;
t62 = -mrSges(5,2) * t107 + t140;
t134 = t61 + t62;
t63 = mrSges(6,1) * t107 - t81 * mrSges(6,3);
t64 = mrSges(5,1) * t107 - t138;
t133 = t63 + t64;
t30 = t110 * t76 - t69;
t100 = t147 * t109;
t55 = t108 * t100 - t110 * t101;
t132 = Ifges(4,4) * t109;
t131 = qJ(5) * t80;
t52 = t107 * t89;
t46 = t52 * qJD(2);
t130 = t108 * t46;
t71 = t110 * t77;
t128 = Ifges(4,6) * qJD(3);
t127 = qJD(2) * t109;
t126 = qJD(2) * t111;
t125 = qJD(3) * t109;
t124 = qJD(4) * t108;
t123 = qJD(4) * t110;
t121 = m(4) * pkin(6) + mrSges(4,3);
t119 = qJD(3) * t147;
t117 = -t128 / 0.2e1;
t116 = qJD(2) * t125;
t29 = -t108 * t76 - t71;
t54 = t110 * t100 + t101 * t108;
t105 = Ifges(4,4) * t126;
t95 = -pkin(6) * t127 + t106;
t115 = Ifges(4,1) * t127 / 0.2e1 + t105 / 0.2e1 + t118 - t95 * mrSges(4,3);
t96 = pkin(6) * t126 + t122;
t114 = t96 * mrSges(4,3) + t128 / 0.2e1 + (Ifges(4,2) * t111 + t132) * qJD(2) / 0.2e1;
t93 = t109 * t119;
t28 = t108 * t72 + t71;
t102 = qJD(3) * t106;
t66 = qJD(2) * t93 + t102;
t67 = (t111 * t120 - t122) * qJD(3);
t7 = t108 * t67 + t110 * t66 + t72 * t123 - t124 * t77;
t94 = t111 * t119;
t14 = t100 * t123 + t101 * t124 + t108 * t94 + t110 * t93;
t8 = -qJD(4) * t28 - t108 * t66 + t110 * t67;
t15 = -qJD(4) * t55 - t108 * t93 + t110 * t94;
t51 = t107 * t88;
t10 = pkin(4) * t107 + t11;
t2 = -qJ(5) * t46 + qJD(5) * t80 + t7;
t45 = t51 * qJD(2);
t3 = -qJ(5) * t45 - qJD(5) * t81 + t8;
t53 = -t80 * pkin(4) + qJD(5) + t99;
t112 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t7 * mrSges(5,2) - t2 * mrSges(6,2) + t10 * t139 + t27 * t140 - t53 * (mrSges(6,1) * t81 + mrSges(6,2) * t80) - t99 * (mrSges(5,1) * t81 + mrSges(5,2) * t80) - t158 * t46 + t159 * t45 - (t155 * t80 - t160) * t81 / 0.2e1 + t157 * t148 - (-t158 * t81 + t159 * t80) * t107 / 0.2e1 + (-t154 * t81 + t156 + t161) * t151;
t103 = pkin(3) * t110 + pkin(4);
t98 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t126;
t97 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t127;
t86 = t96 * qJD(3);
t85 = -pkin(6) * t116 + t102;
t68 = -t88 * pkin(4) + t104;
t60 = pkin(3) * t127 + t143;
t50 = -mrSges(5,1) * t80 + mrSges(5,2) * t81;
t49 = -mrSges(6,1) * t80 + mrSges(6,2) * t81;
t40 = t45 * mrSges(6,2);
t39 = pkin(3) * t125 + pkin(4) * t52;
t38 = qJ(5) * t88 + t55;
t37 = -qJ(5) * t89 + t54;
t26 = pkin(3) * t116 + pkin(4) * t46;
t17 = -t73 + t30;
t16 = t29 - t131;
t12 = t28 + t131;
t5 = -qJ(5) * t51 - qJD(5) * t89 + t15;
t4 = -qJ(5) * t52 + qJD(5) * t88 + t14;
t1 = [-t133 * t52 + t134 * t51 + (-t109 * t97 + t111 * t98 + (-t109 ^ 2 - t111 ^ 2) * qJD(2) * mrSges(4,3)) * qJD(3) + m(4) * (t85 * t109 - t111 * t86 + (-t109 * t95 + t111 * t96) * qJD(3)) + m(5) * (-t27 * t52 + t28 * t51 + t7 * t89 + t8 * t88) + m(6) * (-t10 * t52 + t12 * t51 + t2 * t89 + t3 * t88) + (mrSges(6,3) + mrSges(5,3)) * (-t88 * t45 - t46 * t89); m(6) * (t10 * t5 + t12 * t4 + t2 * t38 + t26 * t68 + t3 * t37 + t39 * t53) + t99 * (mrSges(5,1) * t52 + mrSges(5,2) * t51) + t26 * (-mrSges(6,1) * t88 + mrSges(6,2) * t89) + t4 * t61 + t14 * t62 + t5 * t63 + t15 * t64 + t53 * (mrSges(6,1) * t52 + mrSges(6,2) * t51) + t39 * t49 + m(5) * (t14 * t28 + t15 * t27 + t54 * t8 + t55 * t7) + (-t27 * t51 - t28 * t52 + t7 * t88 - t8 * t89) * mrSges(5,3) + (mrSges(5,2) * t104 - mrSges(5,3) * t54 - mrSges(6,3) * t37 + t135 * t88 + t155 * t89) * t45 - (-mrSges(5,1) * t104 - mrSges(6,1) * t68 + mrSges(5,3) * t55 + mrSges(6,3) * t38 + t135 * t89 + t154 * t88) * t46 + (-t10 * t51 - t12 * t52 + t2 * t88 - t3 * t89) * mrSges(6,3) + (t121 * t85 + (t118 + (-0.2e1 * t144 + 0.3e1 / 0.2e1 * Ifges(4,4) * t111) * qJD(2) + (-m(4) * t95 - t97) * pkin(6) + t115) * qJD(3)) * t111 + (t121 * t86 + (t117 + (-m(4) * t96 - t98) * pkin(6) + (-0.2e1 * t145 - 0.3e1 / 0.2e1 * t132 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t111) * qJD(2) + (t50 + 0.2e1 * t146 + qJD(2) * (-mrSges(5,1) * t88 + mrSges(5,2) * t89)) * pkin(3) - t114) * qJD(3)) * t109 + t68 * t40 + t156 * t51 / 0.2e1 - t157 * t52 / 0.2e1 + (t135 * t51 - t154 * t52) * t80 / 0.2e1 + (-t135 * t52 + t155 * t51) * t148 + (-t158 * t52 + t159 * t51) * t107 / 0.2e1; t112 + (-t103 * t45 + t12 * t81) * mrSges(6,3) - m(5) * (t27 * t29 + t28 * t30) + (-mrSges(6,3) * t130 + (-t110 * t45 - t130) * mrSges(5,3) + (-t108 * t133 + t110 * t134) * qJD(4) + m(5) * (t108 * t7 + t110 * t8 + t123 * t28 - t124 * t27)) * pkin(3) + ((t118 - t105 / 0.2e1 + qJD(2) * t144 - t115) * t111 + (t117 + (t145 + t132 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t111) * qJD(2) + (-t50 - t146) * pkin(3) + t114) * t109) * qJD(2) + t96 * t97 - t95 * t98 - t85 * mrSges(4,2) - t86 * mrSges(4,1) - t17 * t61 - t30 * t62 - t16 * t63 - t29 * t64 - t60 * t49 + t28 * t138 + (t3 * t103 - t10 * t16 - t12 * t17 - t53 * t60 + (-t10 * t124 + t108 * t2 + t12 * t123) * pkin(3)) * m(6); t112 + (-t53 * t143 - (-t10 + t11) * t12 + t3 * pkin(4)) * m(6) + (t28 * mrSges(5,3) + t12 * mrSges(6,3) - pkin(4) * t49) * t81 - t11 * t61 - t27 * t62 + t12 * t63 + t28 * t64 - pkin(4) * t45 * mrSges(6,3); t46 * mrSges(6,1) - t80 * t61 + t81 * t63 + t40 + 0.2e1 * (t26 / 0.2e1 + t10 * t148 + t12 * t151) * m(6);];
tauc = t1(:);
