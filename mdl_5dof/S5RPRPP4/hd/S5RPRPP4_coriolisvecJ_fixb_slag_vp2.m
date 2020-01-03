% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:24
% DurationCPUTime: 2.51s
% Computational Cost: add. (1151->238), mult. (2642->319), div. (0->0), fcn. (1477->4), ass. (0->119)
t147 = Ifges(5,1) + Ifges(6,1);
t149 = Ifges(6,4) + Ifges(5,5);
t148 = mrSges(6,1) + mrSges(5,1);
t146 = Ifges(5,4) - Ifges(6,5);
t111 = cos(pkin(7));
t78 = sin(pkin(7));
t79 = sin(qJ(3));
t80 = cos(qJ(3));
t87 = -t111 * t79 - t78 * t80;
t47 = t87 * qJD(1);
t124 = Ifges(6,5) * t47;
t46 = Ifges(5,4) * t47;
t108 = qJD(1) * t79;
t92 = t111 * t80;
t91 = qJD(1) * t92;
t50 = -t108 * t78 + t91;
t145 = t149 * qJD(3) + t147 * t50 - t124 + t46;
t116 = mrSges(6,2) + mrSges(5,3);
t142 = -Ifges(5,6) + Ifges(6,6);
t128 = mrSges(5,3) * t47;
t130 = mrSges(6,2) * t47;
t34 = qJD(3) * mrSges(6,3) + t130;
t114 = -qJD(3) * mrSges(5,2) + t128 + t34;
t122 = t50 * mrSges(5,3);
t129 = mrSges(6,2) * t50;
t113 = t148 * qJD(3) - t122 - t129;
t141 = (qJ(2) * (m(3) + m(4)));
t140 = qJD(1) ^ 2;
t139 = t47 / 0.2e1;
t138 = -t47 / 0.2e1;
t135 = t50 / 0.2e1;
t81 = -pkin(1) - pkin(6);
t133 = pkin(3) * t78;
t112 = qJ(4) - t81;
t60 = t112 * t79;
t26 = t112 * t92 - t60 * t78;
t101 = qJ(4) * qJD(3);
t104 = qJD(4) * t79;
t105 = qJD(3) * t80;
t64 = qJD(1) * t81 + qJD(2);
t97 = t64 * t105;
t28 = t97 + (-t101 * t80 - t104) * qJD(1);
t103 = qJD(4) * t80;
t106 = qJD(3) * t79;
t98 = t64 * t106;
t83 = -t98 + (t101 * t79 - t103) * qJD(1);
t3 = -t111 * t83 + t28 * t78;
t132 = t26 * t3;
t54 = t78 * t79 - t92;
t131 = t3 * t54;
t127 = Ifges(4,4) * t79;
t126 = Ifges(4,4) * t80;
t125 = Ifges(5,4) * t50;
t100 = qJD(1) * qJD(3);
t41 = t87 * t100;
t123 = t41 * mrSges(6,2);
t102 = qJ(4) * qJD(1);
t43 = (t64 - t102) * t79;
t121 = t78 * t43;
t120 = t79 * Ifges(4,2);
t107 = qJD(1) * t80;
t62 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t107;
t119 = t79 * t62;
t118 = t80 * Ifges(4,1);
t57 = t80 * t64;
t117 = -qJD(3) / 0.2e1;
t4 = t111 * t28 + t78 * t83;
t35 = t111 * t43;
t44 = -t102 * t80 + t57;
t37 = qJD(3) * pkin(3) + t44;
t9 = t78 * t37 + t35;
t94 = t80 * t100;
t58 = pkin(3) * t94 + qJD(1) * qJD(2);
t72 = t79 * pkin(3) + qJ(2);
t110 = Ifges(4,5) * qJD(3);
t109 = Ifges(4,6) * qJD(3);
t65 = pkin(3) * t105 + qJD(2);
t99 = pkin(3) * t107;
t59 = pkin(3) * t108 + qJD(1) * qJ(2) + qJD(4);
t96 = t111 * pkin(3);
t95 = t79 * t100;
t93 = t112 * t80;
t90 = mrSges(4,1) * t79 + mrSges(4,2) * t80;
t88 = qJ(2) * (mrSges(4,1) * t80 - mrSges(4,2) * t79);
t8 = t111 * t37 - t121;
t86 = t106 * t112 - t103;
t2 = qJD(3) * qJD(5) + t4;
t48 = -qJD(3) * t92 + t106 * t78;
t49 = t87 * qJD(3);
t6 = -qJD(3) * pkin(4) + qJD(5) - t8;
t7 = qJD(3) * qJ(5) + t9;
t85 = -t2 * t87 - t48 * t7 - t49 * t6 + t131;
t84 = -t4 * t87 - t48 * t9 + t49 * t8 + t131;
t71 = -t96 - pkin(4);
t69 = qJ(5) + t133;
t61 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t108;
t56 = t90 * qJD(1);
t52 = t110 + (t118 - t127) * qJD(1);
t51 = t109 + (-t120 + t126) * qJD(1);
t45 = Ifges(6,5) * t50;
t42 = -qJD(3) * t93 - t104;
t40 = -qJD(3) * t91 + t78 * t95;
t39 = t41 * mrSges(5,2);
t38 = t40 * mrSges(6,1);
t27 = -t111 * t60 - t78 * t93;
t22 = -pkin(4) * t87 + qJ(5) * t54 + t72;
t21 = -mrSges(5,1) * t47 + mrSges(5,2) * t50;
t20 = -mrSges(6,1) * t47 - mrSges(6,3) * t50;
t17 = Ifges(5,2) * t47 + Ifges(5,6) * qJD(3) + t125;
t16 = Ifges(6,6) * qJD(3) - Ifges(6,3) * t47 + t45;
t15 = pkin(4) * t50 - qJ(5) * t47 + t99;
t14 = t111 * t44 - t121;
t13 = t44 * t78 + t35;
t12 = t111 * t42 + t78 * t86;
t11 = -t111 * t86 + t42 * t78;
t10 = -pkin(4) * t47 - qJ(5) * t50 + t59;
t5 = -pkin(4) * t48 - qJ(5) * t49 + qJD(5) * t54 + t65;
t1 = -pkin(4) * t40 - qJ(5) * t41 - qJD(5) * t50 + t58;
t18 = [m(5) * (-t11 * t8 + t12 * t9 + t27 * t4 + t58 * t72 + t59 * t65 + t132) + m(6) * (t1 * t22 + t10 * t5 + t11 * t6 + t12 * t7 + t2 * t27 + t132) + (-t22 * mrSges(6,3) + t116 * t26 + t146 * t87 - t147 * t54) * t41 + (((2 * mrSges(3,3)) + t90 + (2 * t141)) * qJD(2) + (0.2e1 * t88 + 0.3e1 / 0.2e1 * t80 * t120 - 0.3e1 / 0.2e1 * t79 * t118 + (0.3e1 / 0.2e1 * t79 ^ 2 - 0.3e1 / 0.2e1 * t80 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1) - t22 * t38 + t72 * t39 - t113 * t11 + t114 * t12 + qJD(2) * t56 + t65 * t21 + t5 * t20 - t84 * mrSges(5,3) - t85 * mrSges(6,2) + ((-t51 / 0.2e1 + t81 * t61 - t109 / 0.2e1) * t80 + (-t52 / 0.2e1 - t81 * t62 - t110 / 0.2e1) * t79) * qJD(3) + (-mrSges(5,1) * t72 - (-Ifges(5,2) - Ifges(6,3)) * t87 - t146 * t54 + t116 * t27) * t40 + t58 * (-mrSges(5,1) * t87 - mrSges(5,2) * t54) + t1 * (-mrSges(6,1) * t87 + mrSges(6,3) * t54) + (t17 / 0.2e1 - t16 / 0.2e1 + (Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * qJD(3) - t59 * mrSges(5,1) - t10 * mrSges(6,1) - Ifges(6,3) * t138 + Ifges(5,2) * t139 + t146 * t135) * t48 + ((Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * qJD(3) + t59 * mrSges(5,2) - t10 * mrSges(6,3) + Ifges(6,5) * t138 + Ifges(5,4) * t139 + t145 / 0.2e1 + t147 * t135) * t49; t113 * t49 - t114 * t48 + (t80 * t61 - t119) * qJD(3) + m(5) * t84 + m(6) * t85 + (-mrSges(3,3) - t141) * t140 + (-m(5) * t59 - m(6) * t10 - t20 - t21 - t56) * qJD(1) + t116 * (-t40 * t87 + t41 * t54); (-t10 * t15 - t13 * t6 + t2 * t69 + t3 * t71 + (-t14 + qJD(5)) * t7) * m(6) + t113 * t13 - t114 * t14 + (mrSges(6,2) * t69 + mrSges(5,3) * t133 - t142) * t40 + (t13 * t8 - t14 * t9 - t59 * t99 + (-t111 * t3 + t4 * t78) * pkin(3)) * m(5) + (-Ifges(5,2) * t50 + t145 + t46) * t138 + t8 * t128 - t61 * t57 + t51 * t107 / 0.2e1 + t52 * t108 / 0.2e1 - (-Ifges(4,5) * t79 - Ifges(4,6) * t80) * t100 / 0.2e1 - mrSges(4,2) * t97 - mrSges(4,1) * t98 - t21 * t99 - Ifges(4,5) * t95 - Ifges(4,6) * t94 + (-mrSges(5,3) * t96 + t149) * t41 + (t142 * t50 + t149 * t47) * t117 + (-t88 - t80 * (-Ifges(4,1) * t79 - t126) / 0.2e1 + t79 * (-Ifges(4,2) * t80 - t127) / 0.2e1) * t140 + qJD(5) * t34 - t15 * t20 + t2 * mrSges(6,3) - t4 * mrSges(5,2) + t17 * t135 + t64 * t119 + t9 * t122 + t71 * t123 + t7 * t129 - t6 * t130 - (t147 * t47 - t125 + t16 + t45) * t50 / 0.2e1 - t148 * t3 - t59 * (mrSges(5,1) * t50 + mrSges(5,2) * t47) - t10 * (mrSges(6,1) * t50 - mrSges(6,3) * t47) + (Ifges(6,3) * t50 + t124) * t139; -t40 * mrSges(5,1) - t41 * mrSges(6,3) + t113 * t50 - t114 * t47 - t38 + t39 + (-t7 * t47 - t6 * t50 + t1) * m(6) + (-t47 * t9 + t50 * t8 + t58) * m(5); t123 - qJD(3) * t34 + t50 * t20 + 0.2e1 * (t3 / 0.2e1 + t7 * t117 + t10 * t135) * m(6);];
tauc = t18(:);
