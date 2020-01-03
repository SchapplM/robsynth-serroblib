% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:28
% EndTime: 2020-01-03 11:27:37
% DurationCPUTime: 1.72s
% Computational Cost: add. (2330->218), mult. (5808->318), div. (0->0), fcn. (4200->8), ass. (0->100)
t90 = cos(pkin(9));
t95 = cos(qJ(4));
t108 = t90 * t95;
t88 = sin(pkin(9));
t93 = sin(qJ(4));
t77 = -t88 * t93 + t108;
t70 = t77 * qJD(1);
t78 = t88 * t95 + t90 * t93;
t71 = t78 * qJD(1);
t92 = sin(qJ(5));
t94 = cos(qJ(5));
t102 = t70 * t94 - t71 * t92;
t45 = t70 * t92 + t71 * t94;
t116 = Ifges(6,4) * t45;
t39 = Ifges(6,4) * t102;
t87 = qJD(4) + qJD(5);
t15 = Ifges(6,1) * t45 + Ifges(6,5) * t87 + t39;
t72 = t77 * qJD(4);
t63 = qJD(1) * t72;
t73 = t78 * qJD(4);
t64 = qJD(1) * t73;
t19 = qJD(5) * t102 + t63 * t94 - t64 * t92;
t104 = qJD(4) * t95;
t105 = qJD(3) * t88;
t106 = pkin(6) * qJD(1);
t82 = sin(pkin(8)) * pkin(1) + qJ(3);
t79 = qJD(1) * t82;
t84 = t90 * qJD(2);
t54 = t84 + (-t79 - t106) * t88;
t61 = qJD(2) * t88 + t79 * t90;
t55 = t106 * t90 + t61;
t81 = qJD(3) * t108;
t22 = t54 * t104 + qJD(1) * t81 + (-qJD(1) * t105 - qJD(4) * t55) * t93;
t12 = -pkin(7) * t64 + t22;
t32 = t54 * t93 + t55 * t95;
t97 = t78 * qJD(3);
t23 = -qJD(1) * t97 - qJD(4) * t32;
t13 = -pkin(7) * t63 + t23;
t26 = pkin(7) * t70 + t32;
t114 = t26 * t92;
t31 = t54 * t95 - t55 * t93;
t25 = -pkin(7) * t71 + t31;
t24 = qJD(4) * pkin(4) + t25;
t6 = t24 * t94 - t114;
t2 = qJD(5) * t6 + t12 * t94 + t13 * t92;
t20 = -qJD(5) * t45 - t63 * t92 - t64 * t94;
t113 = t26 * t94;
t7 = t24 * t92 + t113;
t3 = -qJD(5) * t7 - t12 * t92 + t13 * t94;
t98 = -cos(pkin(8)) * pkin(1) - pkin(3) * t90 - pkin(2);
t69 = qJD(1) * t98 + qJD(3);
t48 = -pkin(4) * t70 + t69;
t134 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 - (Ifges(6,5) * t102 - Ifges(6,6) * t45) * t87 / 0.2e1 + (t102 * t6 + t45 * t7) * mrSges(6,3) - (-Ifges(6,2) * t45 + t15 + t39) * t102 / 0.2e1 - t48 * (mrSges(6,1) * t45 + mrSges(6,2) * t102) - (Ifges(6,1) * t102 - t116) * t45 / 0.2e1;
t101 = qJD(1) * (t88 ^ 2 + t90 ^ 2);
t133 = mrSges(4,3) * t101;
t14 = Ifges(6,2) * t102 + Ifges(6,6) * t87 + t116;
t131 = t14 / 0.2e1;
t118 = pkin(6) + t82;
t74 = t118 * t88;
t75 = t118 * t90;
t47 = -t74 * t93 + t75 * t95;
t125 = t102 / 0.2e1;
t123 = t45 / 0.2e1;
t121 = t72 / 0.2e1;
t120 = -t73 / 0.2e1;
t117 = Ifges(5,4) * t71;
t50 = t77 * t92 + t78 * t94;
t115 = t20 * t50;
t49 = t77 * t94 - t78 * t92;
t112 = t49 * t19;
t111 = t64 * t78;
t109 = t77 * t63;
t103 = -t20 * mrSges(6,1) + mrSges(6,2) * t19;
t46 = -t74 * t95 - t75 * t93;
t35 = -pkin(7) * t78 + t46;
t36 = pkin(7) * t77 + t47;
t10 = t35 * t94 - t36 * t92;
t11 = t35 * t92 + t36 * t94;
t100 = -(-t79 * t88 + t84) * t88 + t61 * t90;
t33 = -t74 * t104 + t81 + (-qJD(4) * t75 - t105) * t93;
t34 = -qJD(4) * t47 - t97;
t68 = Ifges(5,4) * t70;
t58 = t63 * mrSges(5,2);
t57 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t71;
t56 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t70;
t53 = -pkin(4) * t77 + t98;
t41 = Ifges(5,1) * t71 + Ifges(5,5) * qJD(4) + t68;
t40 = Ifges(5,2) * t70 + Ifges(5,6) * qJD(4) + t117;
t38 = mrSges(6,1) * t87 - mrSges(6,3) * t45;
t37 = -mrSges(6,2) * t87 + mrSges(6,3) * t102;
t30 = -pkin(7) * t72 + t34;
t29 = -pkin(7) * t73 + t33;
t28 = -qJD(5) * t50 - t72 * t92 - t73 * t94;
t27 = qJD(5) * t49 + t72 * t94 - t73 * t92;
t21 = -mrSges(6,1) * t102 + mrSges(6,2) * t45;
t9 = t25 * t94 - t114;
t8 = -t25 * t92 - t113;
t5 = -qJD(5) * t11 - t29 * t92 + t30 * t94;
t4 = qJD(5) * t10 + t29 * t94 + t30 * t92;
t1 = [(m(4) * (t101 * t82 + t100) + 0.2e1 * t133) * qJD(3) + t87 * (Ifges(6,5) * t27 + Ifges(6,6) * t28) / 0.2e1 + (t123 * t27 + t19 * t50) * Ifges(6,1) + (t125 * t28 + t49 * t20) * Ifges(6,2) + (t123 * t28 + t125 * t27 + t112 + t115) * Ifges(6,4) + (t121 * t71 + t63 * t78) * Ifges(5,1) + t53 * t103 + t33 * t56 + t34 * t57 + t48 * (-mrSges(6,1) * t28 + mrSges(6,2) * t27) + t4 * t37 + t5 * t38 + t27 * t15 / 0.2e1 + (t120 * t70 - t77 * t64) * Ifges(5,2) + t98 * (t64 * mrSges(5,1) + t58) + (t64 * (-mrSges(6,1) * t49 + mrSges(6,2) * t50) + t73 * t21) * pkin(4) + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + (t48 * t73 + t53 * t64) * pkin(4)) + (t22 * t77 - t23 * t78 - t31 * t72 - t32 * t73 - t46 * t63 - t47 * t64) * mrSges(5,3) + (t120 * t71 + t121 * t70 + t109 - t111) * Ifges(5,4) + m(5) * (t22 * t47 + t23 * t46 + t31 * t34 + t32 * t33) + t40 * t120 + t41 * t121 + t69 * (mrSges(5,1) * t73 + mrSges(5,2) * t72) + qJD(4) * (Ifges(5,5) * t72 - Ifges(5,6) * t73) / 0.2e1 + t28 * t131 + (-t10 * t19 + t11 * t20 + t2 * t49 - t27 * t6 + t28 * t7 - t3 * t50) * mrSges(6,3); t27 * t37 + t28 * t38 + t72 * t56 - t73 * t57 + (-t112 + t115) * mrSges(6,3) + (-t109 - t111) * mrSges(5,3) + m(5) * (t22 * t78 + t23 * t77 - t31 * t73 + t32 * t72) + m(6) * (t2 * t50 + t27 * t7 + t28 * t6 + t3 * t49); -t102 * t37 + t45 * t38 - t70 * t56 + t71 * t57 + t58 - (-m(6) * pkin(4) - mrSges(5,1)) * t64 - m(5) * (-t31 * t71 + t32 * t70) - m(6) * (t102 * t7 - t45 * t6) + t103 + (-m(4) * t100 - t133) * qJD(1); (t31 * t70 + t32 * t71) * mrSges(5,3) - t69 * (mrSges(5,1) * t71 + mrSges(5,2) * t70) - qJD(4) * (Ifges(5,5) * t70 - Ifges(5,6) * t71) / 0.2e1 + t71 * t40 / 0.2e1 - Ifges(5,6) * t64 - t71 * (Ifges(5,1) * t70 - t117) / 0.2e1 - t31 * t56 + t32 * t57 + Ifges(5,5) * t63 - t9 * t37 - t8 * t38 - t22 * mrSges(5,2) + t23 * mrSges(5,1) - (-Ifges(5,2) * t71 + t41 + t68) * t70 / 0.2e1 + t45 * t131 - m(6) * (t6 * t8 + t7 * t9) + (-t71 * t21 + (t94 * t37 - t92 * t38) * qJD(5) + (-t19 * t94 + t20 * t92) * mrSges(6,3) + (t2 * t92 + t3 * t94 - t48 * t71 + (-t6 * t92 + t7 * t94) * qJD(5)) * m(6)) * pkin(4) + t134; t14 * t123 - t6 * t37 + t7 * t38 + t134;];
tauc = t1(:);
