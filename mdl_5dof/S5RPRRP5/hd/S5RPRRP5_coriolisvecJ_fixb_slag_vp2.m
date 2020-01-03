% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:44
% DurationCPUTime: 0.97s
% Computational Cost: add. (1209->191), mult. (2576->262), div. (0->0), fcn. (1234->6), ass. (0->88)
t144 = mrSges(5,1) + mrSges(6,1);
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t143 = -mrSges(5,1) * t82 + mrSges(5,2) * t80 - mrSges(4,1);
t114 = qJD(4) * t82;
t138 = pkin(1) * sin(pkin(8));
t108 = qJD(1) * t138;
t72 = cos(pkin(8)) * pkin(1) + pkin(2);
t65 = t72 * qJD(1);
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t35 = t83 * t108 + t81 * t65;
t75 = qJD(1) + qJD(3);
t21 = t75 * pkin(7) + t35;
t126 = t21 * t80;
t16 = qJD(2) * t82 - t126;
t140 = -t16 + qJD(5);
t12 = -qJD(4) * pkin(4) + t140;
t17 = qJD(2) * t80 + t21 * t82;
t34 = -t81 * t108 + t65 * t83;
t28 = t34 * qJD(3);
t6 = qJD(4) * t17 + t28 * t80;
t134 = t6 * t80;
t121 = qJD(2) * t114 + t82 * t28;
t4 = (qJD(5) - t126) * qJD(4) + t121;
t142 = t12 * t114 + t4 * t82 + t134;
t118 = t83 * t138 + t81 * t72;
t101 = -t81 * t138 + t72 * t83;
t123 = t75 * t82;
t109 = mrSges(6,2) * t123;
t63 = qJD(4) * mrSges(6,3) + t109;
t119 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t123 + t63;
t122 = t80 * mrSges(5,3);
t124 = t75 * t80;
t120 = -mrSges(6,2) * t124 + t144 * qJD(4) - t75 * t122;
t139 = t119 * t82 - t120 * t80;
t115 = qJD(4) * t80;
t5 = -t21 * t115 + t121;
t136 = t5 * t82;
t53 = pkin(7) + t118;
t135 = t53 * t6;
t133 = t6 * t82;
t132 = Ifges(5,4) * t80;
t130 = Ifges(6,5) * t80;
t129 = Ifges(6,5) * t82;
t111 = qJD(4) * qJ(5);
t14 = t17 + t111;
t128 = t14 * mrSges(6,2);
t127 = t16 * mrSges(5,3);
t116 = qJD(4) * t53;
t113 = t14 * qJD(4);
t112 = t80 * qJD(5);
t106 = t124 / 0.2e1;
t104 = t115 / 0.2e1;
t102 = t143 * t75;
t100 = t75 * t104;
t97 = -mrSges(6,1) * t82 - mrSges(6,3) * t80;
t96 = -t17 * mrSges(5,3) - t128;
t95 = pkin(4) * t80 - qJ(5) * t82;
t94 = t12 * t80 + t14 * t82;
t93 = -t16 * t80 + t17 * t82;
t64 = -t82 * pkin(4) - t80 * qJ(5) - pkin(3);
t92 = (Ifges(6,1) * t80 - t129) * t75;
t91 = (Ifges(5,2) * t82 + t132) * t75;
t90 = (mrSges(5,1) * t80 + mrSges(5,2) * t82) * qJD(4);
t89 = (mrSges(6,1) * t80 - mrSges(6,3) * t82) * qJD(4);
t41 = t118 * qJD(3);
t54 = pkin(4) * t115 - t82 * t111 - t112;
t29 = t35 * qJD(3);
t11 = t64 * t75 - t34;
t20 = -t75 * pkin(3) - t34;
t68 = Ifges(6,5) * t124;
t42 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t123 + t68;
t43 = Ifges(5,6) * qJD(4) + t91;
t44 = Ifges(6,4) * qJD(4) + t92;
t69 = Ifges(5,4) * t123;
t45 = Ifges(5,1) * t124 + Ifges(5,5) * qJD(4) + t69;
t9 = t29 + (t95 * qJD(4) - t112) * t75;
t84 = -t75 * (Ifges(6,3) * t80 + t129) * t114 + t9 * t97 + mrSges(5,3) * t136 + t11 * t89 + t20 * t90 + (-Ifges(6,3) * t82 + t130) * t100 + t42 * t104 + t143 * t29 + 0.2e1 * (t130 - t132 + (Ifges(5,1) + Ifges(6,1)) * t82) * t100 + ((Ifges(6,4) + Ifges(5,5)) * t82 + (-Ifges(5,6) + Ifges(6,6)) * t80) * qJD(4) ^ 2 / 0.2e1 - (t91 + t43) * t115 / 0.2e1 + t142 * mrSges(6,2) + ((0.3e1 * Ifges(5,4) * t82 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t80) * t75 + t92 + t45 + t44) * t114 / 0.2e1;
t57 = t95 * t75;
t55 = t97 * t75;
t52 = -pkin(3) - t101;
t47 = t75 * t90;
t46 = t75 * t89;
t40 = t101 * qJD(3);
t22 = -t101 + t64;
t15 = t41 + t54;
t1 = [t84 + t102 * t41 + m(6) * (t11 * t15 + t22 * t9) + m(5) * (t20 * t41 + t29 * t52) + m(4) * (-t29 * t101 + t28 * t118 - t34 * t41 + t35 * t40) + (-t40 * t75 - t28) * mrSges(4,2) + t22 * t46 + t52 * t47 + t15 * t55 + (t6 * mrSges(5,3) - t120 * t40 + (-t119 * t53 + t96) * qJD(4) + m(6) * (-t53 * t113 + t12 * t40 + t135) + m(5) * (-t17 * t116 - t16 * t40 + t135)) * t80 + (t119 * t40 + (-t120 * t53 - t127) * qJD(4) + m(6) * (t12 * t116 + t14 * t40 + t4 * t53) + m(5) * (-t16 * t116 + t17 * t40 + t5 * t53)) * t82; m(5) * (t5 * t80 - t133) + m(6) * (t4 * t80 - t133) + (m(5) * t93 + m(6) * t94 + (mrSges(6,2) + mrSges(5,3)) * t75 * (-t80 ^ 2 - t82 ^ 2) + t139) * qJD(4); t84 + (-t80 * t128 + (-t16 * t82 - t17 * t80) * mrSges(5,3)) * qJD(4) + (t75 * mrSges(4,2) - t139) * t34 + (-t102 - t55) * t35 + t6 * t122 + t64 * t46 - pkin(3) * t47 + t54 * t55 - t28 * mrSges(4,2) + (-t94 * t34 + t9 * t64 + (t54 - t35) * t11) * m(6) + (-t29 * pkin(3) - t20 * t35 - t93 * t34) * m(5) + (m(6) * (-t80 * t113 + t142) + m(5) * (-t16 * t114 - t17 * t115 + t134 + t136) + (-t119 * t80 - t120 * t82) * qJD(4)) * pkin(7); -t5 * mrSges(5,2) + t4 * mrSges(6,3) + qJD(5) * t63 - t57 * t55 - t144 * t6 + t120 * t17 - t119 * t16 + ((-t69 / 0.2e1 - t44 / 0.2e1 - t45 / 0.2e1 - t20 * mrSges(5,2) + t11 * mrSges(6,3) - t12 * mrSges(6,2) + t127 + Ifges(6,5) * t123 / 0.2e1 + (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1 - pkin(4) * mrSges(6,2)) * qJD(4)) * t82 + (-t68 / 0.2e1 - t42 / 0.2e1 + t43 / 0.2e1 - t20 * mrSges(5,1) - t11 * mrSges(6,1) + Ifges(5,4) * t106 + (-Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1 - qJ(5) * mrSges(6,2)) * qJD(4) + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t123 - t96) * t80) * t75 + (-t6 * pkin(4) + t4 * qJ(5) - t11 * t57 - t12 * t17 + t14 * t140) * m(6); t55 * t124 + (-t63 + t109) * qJD(4) + 0.2e1 * (t6 / 0.2e1 + t11 * t106 - t113 / 0.2e1) * m(6);];
tauc = t1(:);
