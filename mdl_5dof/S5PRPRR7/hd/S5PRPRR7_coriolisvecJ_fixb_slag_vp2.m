% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:45
% EndTime: 2019-12-05 15:59:48
% DurationCPUTime: 1.40s
% Computational Cost: add. (1258->209), mult. (2711->306), div. (0->0), fcn. (1554->6), ass. (0->112)
t86 = -pkin(2) - pkin(6);
t130 = pkin(7) - t86;
t81 = sin(qJ(4));
t60 = t130 * t81;
t84 = cos(qJ(4));
t61 = t130 * t84;
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t31 = t60 * t80 - t61 * t83;
t107 = qJD(4) * t81;
t54 = t130 * t107;
t55 = qJD(4) * t61;
t82 = sin(qJ(2));
t95 = t80 * t84 + t83 * t81;
t92 = t95 * t82;
t141 = -qJD(1) * t92 + qJD(5) * t31 + t54 * t80 - t55 * t83;
t32 = -t60 * t83 - t61 * t80;
t115 = t83 * t84;
t94 = t80 * t81 - t115;
t91 = t94 * t82;
t140 = qJD(1) * t91 - qJD(5) * t32 + t54 * t83 + t55 * t80;
t85 = cos(qJ(2));
t103 = t85 * qJD(1);
t99 = qJD(3) - t103;
t56 = qJD(2) * t86 + t99;
t109 = qJD(2) * t82;
t100 = qJD(1) * t109;
t68 = t84 * t100;
t33 = -t107 * t56 + t68;
t106 = qJD(4) * t84;
t34 = t81 * t100 + t56 * t106;
t97 = t33 * t84 + t34 * t81;
t43 = t95 * t85;
t110 = qJD(2) * t81;
t64 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t110;
t108 = qJD(2) * t84;
t65 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t108;
t96 = t64 * t84 - t65 * t81;
t139 = t96 * qJD(4);
t77 = qJD(4) + qJD(5);
t111 = qJD(1) * t82;
t66 = qJD(2) * qJ(3) + t111;
t116 = t66 * t85;
t138 = t116 + (t81 ^ 2 + t84 ^ 2) * t56 * t82;
t30 = t77 * t95;
t136 = -t30 / 0.2e1;
t50 = t95 * qJD(2);
t135 = t50 / 0.2e1;
t51 = t108 * t83 - t110 * t80;
t134 = -t51 / 0.2e1;
t133 = t51 / 0.2e1;
t132 = t81 / 0.2e1;
t131 = -t84 / 0.2e1;
t129 = mrSges(6,3) * t50;
t128 = Ifges(5,4) * t81;
t127 = Ifges(5,4) * t84;
t126 = Ifges(6,4) * t51;
t38 = -pkin(7) * t110 + t56 * t81;
t124 = t38 * t80;
t123 = t38 * t83;
t122 = t51 * mrSges(6,3);
t73 = pkin(4) * t81 + qJ(3);
t53 = qJD(2) * t73 + t111;
t121 = t53 * t85;
t88 = t77 * t94;
t24 = t88 * qJD(2);
t119 = t95 * t24;
t23 = t30 * qJD(2);
t118 = t94 * t23;
t62 = (qJD(3) + t103) * qJD(2);
t117 = t62 * t82;
t26 = mrSges(6,1) * t50 + mrSges(6,2) * t51;
t59 = (mrSges(5,1) * t81 + mrSges(5,2) * t84) * qJD(2);
t114 = t26 + t59;
t113 = Ifges(5,5) * qJD(4);
t112 = Ifges(5,6) * qJD(4);
t105 = qJD(5) * t80;
t104 = qJD(5) * t83;
t102 = pkin(7) * t108;
t101 = -m(6) * t53 - t26;
t39 = t84 * t56 - t102;
t69 = pkin(4) * t106 + qJD(3);
t98 = mrSges(5,1) * t84 - mrSges(5,2) * t81;
t35 = qJD(4) * pkin(4) + t39;
t12 = t35 * t83 - t124;
t13 = t35 * t80 + t123;
t93 = qJ(3) * t62 + qJD(3) * t66;
t27 = t68 + (pkin(7) * qJD(2) - t56) * t107;
t28 = -qJD(4) * t102 + t34;
t2 = qJD(5) * t12 + t27 * t80 + t28 * t83;
t29 = -t105 * t81 - t107 * t80 + t115 * t77;
t3 = -qJD(5) * t13 + t27 * t83 - t28 * t80;
t90 = -t12 * t30 + t13 * t29 + t2 * t95 - t3 * t94;
t19 = -Ifges(6,2) * t50 + Ifges(6,6) * t77 + t126;
t46 = Ifges(6,4) * t50;
t20 = Ifges(6,1) * t51 + Ifges(6,5) * t77 - t46;
t89 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t12 * t129 + t19 * t133 - t53 * (mrSges(6,1) * t51 - mrSges(6,2) * t50) + (-Ifges(6,1) * t50 - t126) * t134 - t77 * (-Ifges(6,5) * t50 - Ifges(6,6) * t51) / 0.2e1 + Ifges(6,6) * t24 - Ifges(6,5) * t23 + (-Ifges(6,2) * t51 + t20 - t46) * t135;
t87 = qJD(2) ^ 2;
t63 = -qJD(2) * pkin(2) + t99;
t52 = t98 * qJD(4) * qJD(2);
t49 = t113 + (Ifges(5,1) * t84 - t128) * qJD(2);
t48 = t112 + (-Ifges(5,2) * t81 + t127) * qJD(2);
t45 = (t69 + t103) * qJD(2);
t42 = t94 * t85;
t37 = mrSges(6,1) * t77 - t122;
t36 = -mrSges(6,2) * t77 - t129;
t17 = t39 * t83 - t124;
t16 = -t39 * t80 - t123;
t9 = -qJD(2) * t91 + t43 * t77;
t8 = qJD(2) * t92 + t85 * t88;
t4 = -mrSges(6,1) * t24 - mrSges(6,2) * t23;
t1 = [t8 * t36 + t9 * t37 + (t4 + t52) * t82 - t85 * t139 + (t23 * t42 - t24 * t43) * mrSges(6,3) + (t114 * t85 + (t64 * t81 + t65 * t84) * t82) * qJD(2) + m(4) * (t117 + (t116 + (t63 - t103) * t82) * qJD(2)) + m(5) * (qJD(2) * t138 - t97 * t85 + t117) + m(6) * (qJD(2) * t121 + t12 * t9 + t13 * t8 - t2 * t43 + t3 * t42 + t45 * t82) + ((-mrSges(3,2) + mrSges(4,3)) * t85 + (-mrSges(3,1) + mrSges(4,2)) * t82) * t87; t77 * (-Ifges(6,5) * t30 - Ifges(6,6) * t29) / 0.2e1 + t69 * t26 + t73 * t4 + qJ(3) * t52 + t53 * (t29 * mrSges(6,1) - t30 * mrSges(6,2)) + t45 * (mrSges(6,1) * t95 - mrSges(6,2) * t94) + qJD(3) * t59 - t29 * t19 / 0.2e1 + t20 * t136 + t140 * t37 + t141 * t36 - t114 * t103 + (t135 * t29 - t119) * Ifges(6,2) + (-t133 * t30 + t118) * Ifges(6,1) + (qJD(2) * t99 + t62) * mrSges(4,3) + (t23 * t31 + t24 * t32 - t90) * mrSges(6,3) + (-t65 * t111 + t62 * mrSges(5,2) - t33 * mrSges(5,3) + (t66 * mrSges(5,1) + t86 * t64 - t48 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4) * t108 - t112 / 0.2e1) * qJD(4)) * t84 + (t134 * t29 - t136 * t50 + t23 * t95 - t24 * t94) * Ifges(6,4) + (-t64 * t111 + t62 * mrSges(5,1) - t34 * mrSges(5,3) + (-t86 * t65 - t49 / 0.2e1 - t66 * mrSges(5,2) - t113 / 0.2e1 + (0.3e1 / 0.2e1 * t128 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t84) * qJD(2)) * qJD(4)) * t81 + (t93 - (pkin(2) * t109 + t63 * t82 + t116) * qJD(1)) * m(4) + (-qJD(1) * t138 + t97 * t86 + t93) * m(5) + (-qJD(1) * t121 + t140 * t12 + t141 * t13 + t2 * t32 + t3 * t31 + t45 * t73 + t53 * t69) * m(6); -t87 * mrSges(4,3) + t29 * t36 - t30 * t37 + t139 + (-t118 + t119) * mrSges(6,3) + m(5) * t97 + m(6) * t90 + (-m(5) * t66 - t59 + (-t66 + t111) * m(4) + t101) * qJD(2); t89 - m(6) * (t12 * t16 + t13 * t17) - t96 * t56 - t34 * mrSges(5,2) - t17 * t36 - t16 * t37 + t33 * mrSges(5,1) + t13 * t122 + (t36 * t104 + m(6) * (t104 * t13 - t105 * t12 + t2 * t80 + t3 * t83) - t37 * t105 + (t23 * t83 + t24 * t80) * mrSges(6,3)) * pkin(4) + (t84 * t48 / 0.2e1 + t49 * t132 - t66 * t98 + ((-Ifges(5,1) * t81 - t127) * t131 + (-Ifges(5,2) * t84 - t128) * t132) * qJD(2) + t101 * t84 * pkin(4) + (Ifges(5,6) * t131 - Ifges(5,5) * t81 / 0.2e1) * qJD(4)) * qJD(2); t89 - t12 * t36 + (t37 + t122) * t13;];
tauc = t1(:);
