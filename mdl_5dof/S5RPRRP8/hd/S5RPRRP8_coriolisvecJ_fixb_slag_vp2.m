% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:11
% DurationCPUTime: 1.28s
% Computational Cost: add. (1418->201), mult. (2381->272), div. (0->0), fcn. (923->4), ass. (0->96)
t138 = mrSges(5,1) + mrSges(6,1);
t137 = mrSges(6,2) + mrSges(5,3);
t136 = Ifges(6,4) + Ifges(5,5);
t64 = sin(qJ(4));
t112 = Ifges(6,5) * t64;
t66 = cos(qJ(4));
t113 = Ifges(5,4) * t66;
t131 = -qJD(4) / 0.2e1;
t134 = -Ifges(6,3) / 0.2e1;
t114 = Ifges(5,4) * t64;
t51 = Ifges(5,2) * t66 + t114;
t111 = Ifges(6,5) * t66;
t52 = Ifges(6,1) * t64 - t111;
t61 = -qJD(1) + qJD(3);
t130 = qJD(4) / 0.2e1;
t132 = t61 / 0.2e1;
t106 = t61 * t66;
t56 = Ifges(5,4) * t106;
t68 = -pkin(1) - pkin(2);
t54 = t68 * qJD(1) + qJD(2);
t65 = sin(qJ(3));
t67 = cos(qJ(3));
t91 = qJ(2) * qJD(1);
t28 = t54 * t65 + t67 * t91;
t13 = pkin(7) * t61 + t28;
t110 = t13 * t64;
t7 = -qJD(4) * pkin(4) + qJD(5) + t110;
t107 = t61 * t64;
t88 = t107 / 0.2e1;
t73 = t52 * t132 + Ifges(5,1) * t88 + t56 / 0.2e1 + t7 * mrSges(6,2) + t136 * t130;
t133 = -t61 / 0.2e1;
t55 = Ifges(6,5) * t107;
t8 = qJD(4) * qJ(5) + t13 * t66;
t74 = Ifges(6,6) * t130 + t106 * t134 + t55 / 0.2e1 + Ifges(5,6) * t131 + t51 * t133 - t8 * mrSges(6,2);
t135 = ((t66 * t134 + t112 / 0.2e1 - t51 / 0.2e1) * t61 + t74) * t64 + ((t52 / 0.2e1 + Ifges(5,1) * t64 / 0.2e1 + t113 / 0.2e1) * t61 + t73) * t66 - (t136 * t66 + (-Ifges(5,6) + Ifges(6,6)) * t64) * t131;
t102 = t138 * qJD(4) - t137 * t107;
t99 = t67 * qJ(2) + t65 * t68;
t79 = -t65 * qJ(2) + t67 * t68;
t129 = (mrSges(6,1) * t64 - mrSges(6,3) * t66) * qJD(4);
t128 = (m(3) * qJ(2) + mrSges(3,3)) * qJD(1);
t127 = t64 * t7 + t66 * t8;
t90 = qJ(2) * qJD(3);
t95 = qJD(3) * t54;
t97 = qJD(2) * t67;
t10 = t67 * t95 + (-t65 * t90 + t97) * qJD(1);
t9 = t66 * t10;
t1 = t9 + (qJD(5) - t110) * qJD(4);
t93 = qJD(4) * t66;
t3 = t10 * t64 + t13 * t93;
t117 = t3 * t64;
t126 = t1 * t66 + t117;
t89 = mrSges(6,2) * t106;
t46 = qJD(4) * mrSges(6,3) + t89;
t125 = -t8 * m(6) - t46;
t81 = t13 * (t64 ^ 2 + t66 ^ 2);
t39 = (Ifges(6,1) * t66 + t112) * qJD(4);
t40 = (Ifges(5,1) * t66 - t114) * qJD(4);
t124 = -t137 * t3 + (t39 + t40) * t133;
t98 = qJD(2) * t65;
t11 = t65 * t95 + (t67 * t90 + t98) * qJD(1);
t27 = t54 * t67 - t65 * t91;
t12 = -pkin(3) * t61 - t27;
t34 = (mrSges(5,1) * t64 + mrSges(5,2) * t66) * qJD(4);
t75 = pkin(4) * t64 - qJ(5) * t66;
t29 = t75 * qJD(4) - qJD(5) * t64;
t4 = t29 * t61 + t11;
t48 = mrSges(6,1) * t66 + mrSges(6,3) * t64;
t49 = -mrSges(5,1) * t66 + mrSges(5,2) * t64;
t47 = -pkin(4) * t66 - qJ(5) * t64 - pkin(3);
t5 = t47 * t61 - t27;
t121 = -t10 * mrSges(4,2) + (t49 - mrSges(4,1)) * t11 + t12 * t34 + t5 * t129 - t4 * t48;
t35 = (Ifges(6,3) * t64 + t111) * qJD(4);
t38 = (-Ifges(5,2) * t64 + t113) * qJD(4);
t120 = -(t39 / 0.2e1 + t40 / 0.2e1) * t64 + (t35 / 0.2e1 - t38 / 0.2e1) * t66;
t45 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t106;
t101 = t45 + t46;
t119 = m(4) * t28 + m(6) * t127 - t61 * mrSges(4,2) + t101 * t66 - t102 * t64;
t30 = t48 * t61;
t31 = t49 * t61;
t103 = t30 - t31;
t94 = qJD(4) * t64;
t78 = t61 * mrSges(4,1) + t103;
t2 = -t13 * t94 + t9;
t77 = t2 * t66 + t117;
t72 = -t1 * mrSges(6,2) - t2 * mrSges(5,3) + t35 * t132 + t38 * t133;
t26 = t99 * qJD(3) + t98;
t70 = (m(6) * t7 - t102) * t66 + (t125 - t45) * t64;
t42 = -pkin(7) + t99;
t41 = pkin(3) - t79;
t32 = t75 * t61;
t25 = t79 * qJD(3) + t97;
t20 = t61 * t34;
t19 = t61 * t129;
t14 = -t47 - t79;
t6 = t26 - t29;
t15 = [t14 * t19 + t41 * t20 + t26 * t31 - t6 * t30 + 0.2e1 * qJD(2) * t128 + (t101 * t25 + t72) * t66 + (-t102 * t25 + t124) * t64 + m(6) * (t126 * t42 + t127 * t25 + t14 * t4 + t5 * t6) + m(4) * (t10 * t99 - t11 * t79 + t28 * t25 - t27 * t26) + m(5) * (t11 * t41 + t12 * t26 + t25 * t81 + t77 * t42) + (-t26 * mrSges(4,1) - t25 * mrSges(4,2) + t120) * t61 + (t70 * t42 - t135) * qJD(4) - t121; (-t78 * qJD(3) + (-t101 * t64 - t102 * t66) * qJD(4) + m(5) * (qJD(3) * t12 + t77) + m(6) * (qJD(3) * t5 + t7 * t93 - t8 * t94 + t126) + m(4) * (-qJD(3) * t27 + t10)) * t65 + (-t128 + (m(4) * t27 - m(5) * t12 - m(6) * t5 + t78) * t65) * qJD(1) + (-t19 - t20 + m(5) * (qJD(3) * t81 - t11) - m(6) * t4 - m(4) * t11 + t119 * qJD(3) + (-m(5) * t81 - t119) * qJD(1)) * t67; -pkin(3) * t20 + t47 * t19 - t29 * t30 + t103 * t28 + (-t101 * t27 - t72) * t66 + (t102 * t27 - t124) * t64 + (t28 * mrSges(4,1) + t27 * mrSges(4,2) - t120) * t61 + (t70 * pkin(7) + t135) * qJD(4) + (-t127 * t27 + t4 * t47 + (-t28 + t29) * t5 + t126 * pkin(7)) * m(6) + (-pkin(3) * t11 + t77 * pkin(7) - t12 * t28 - t27 * t81) * m(5) + t121; -t2 * mrSges(5,2) + t1 * mrSges(6,3) + qJD(5) * t46 + t32 * t30 - t138 * t3 - t70 * t13 + ((-t56 / 0.2e1 - t12 * mrSges(5,2) + t5 * mrSges(6,3) + Ifges(6,5) * t106 / 0.2e1 - t73) * t66 + (-t55 / 0.2e1 - t12 * mrSges(5,1) - t5 * mrSges(6,1) + Ifges(5,4) * t88 + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t106 - t74) * t64 + ((Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,2)) * t66 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1 - qJ(5) * mrSges(6,2)) * t64) * qJD(4)) * t61 + (-t3 * pkin(4) + t1 * qJ(5) + t8 * qJD(5) - t5 * t32) * m(6); -t30 * t107 + 0.2e1 * (t3 / 0.2e1 + t5 * t88) * m(6) + (t125 + t89) * qJD(4);];
tauc = t15(:);
