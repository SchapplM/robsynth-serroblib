% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:12:56
% DurationCPUTime: 0.95s
% Computational Cost: add. (706->171), mult. (1358->237), div. (0->0), fcn. (550->4), ass. (0->87)
t132 = Ifges(4,1) + Ifges(5,1);
t111 = pkin(1) * qJD(1);
t68 = qJD(1) + qJD(2);
t72 = sin(qJ(2));
t53 = t68 * pkin(6) + t72 * t111;
t71 = sin(qJ(3));
t73 = cos(qJ(3));
t131 = (t71 ^ 2 + t73 ^ 2) * t53;
t74 = cos(qJ(2));
t100 = t74 * t111;
t65 = -t73 * pkin(3) - pkin(2);
t23 = t65 * t68 + qJD(4) - t100;
t86 = -t73 * mrSges(5,1) + t71 * mrSges(5,2);
t41 = t86 * t68;
t126 = -m(5) * t23 - t41;
t95 = qJ(4) * t68 + t53;
t22 = t95 * t73;
t21 = t95 * t71;
t9 = qJD(3) * pkin(3) - t21;
t125 = m(5) * (t22 * t73 - t71 * t9);
t123 = t74 * pkin(1);
t122 = Ifges(4,4) * t71;
t120 = Ifges(5,4) * t71;
t107 = qJD(3) * t71;
t110 = pkin(1) * qJD(2);
t98 = qJD(1) * t110;
t91 = t74 * t98;
t56 = t73 * t91;
t10 = -t53 * t107 + t56;
t118 = t10 * t73;
t115 = t68 * t71;
t114 = t68 * t73;
t113 = t71 * mrSges(4,3);
t112 = -qJ(4) - pkin(6);
t106 = qJD(3) * t73;
t99 = t68 * t107;
t32 = t68 * mrSges(5,2) * t106 + mrSges(5,1) * t99;
t63 = t72 * pkin(1) + pkin(6);
t109 = -qJ(4) - t63;
t108 = qJD(2) * t72;
t66 = t73 * qJD(4);
t105 = -qJD(1) - t68;
t104 = -qJD(2) + t68;
t102 = t74 * t110;
t101 = pkin(3) * t107;
t94 = qJD(3) * t112;
t93 = qJD(3) * t109;
t92 = t72 * t98;
t88 = -t22 * t71 - t9 * t73;
t87 = -t73 * mrSges(4,1) + t71 * mrSges(4,2);
t50 = qJD(3) * mrSges(4,1) - t68 * t113;
t52 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t114;
t85 = t73 * t50 + t71 * t52;
t84 = qJD(3) * t95;
t11 = -t53 * t106 - t71 * t91;
t83 = m(4) * (-t11 * t71 + t118);
t82 = (Ifges(4,2) * t73 + t122) * t68;
t81 = (Ifges(5,2) * t73 + t120) * t68;
t80 = (mrSges(4,1) * t71 + mrSges(4,2) * t73) * qJD(3);
t2 = t68 * t66 - t71 * t84 + t56;
t28 = Ifges(5,6) * qJD(3) + t81;
t29 = Ifges(4,6) * qJD(3) + t82;
t61 = Ifges(5,4) * t114;
t30 = Ifges(5,1) * t115 + Ifges(5,5) * qJD(3) + t61;
t62 = Ifges(4,4) * t114;
t31 = Ifges(4,1) * t115 + Ifges(4,5) * qJD(3) + t62;
t40 = pkin(3) * t99 + t92;
t54 = -t68 * pkin(2) - t100;
t75 = t2 * t73 * mrSges(5,3) + t54 * t80 + t23 * (mrSges(5,1) * t71 + mrSges(5,2) * t73) * qJD(3) + t87 * t92 + mrSges(4,3) * t118 + t40 * t86 + (t132 * t73 - t120 - t122) * t99 + ((Ifges(4,5) + Ifges(5,5)) * t73 + (-Ifges(4,6) - Ifges(5,6)) * t71) * qJD(3) ^ 2 / 0.2e1 - (t28 + t29 + t81 + t82) * t107 / 0.2e1 + (t30 + t31 + ((t132 - 0.2e1 * Ifges(4,2) - 0.2e1 * Ifges(5,2)) * t71 + 0.3e1 * (Ifges(4,4) + Ifges(5,4)) * t73) * t68) * t106 / 0.2e1;
t67 = t73 * qJ(4);
t64 = -pkin(2) - t123;
t58 = t73 * pkin(6) + t67;
t57 = t112 * t71;
t55 = t65 - t123;
t51 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t114;
t49 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t115;
t48 = pkin(1) * t108 + t101;
t47 = t73 * t63 + t67;
t46 = t109 * t71;
t42 = t87 * t68;
t39 = -t71 * qJD(4) + t73 * t94;
t38 = t71 * t94 + t66;
t33 = t68 * t80;
t7 = (-qJD(4) - t102) * t71 + t73 * t93;
t6 = t73 * t102 + t71 * t93 + t66;
t3 = (-qJD(4) * t68 - t91) * t71 - t73 * t84;
t1 = [((m(4) * (qJD(1) * t64 + t54) + t42 + t105 * mrSges(3,1)) * t72 + (m(4) * t131 + t105 * mrSges(3,2) - t71 * t50 + t73 * t52) * t74) * t110 + t75 + (-t85 * t63 + ((-t46 * t73 - t47 * t71) * t68 + t88) * mrSges(5,3)) * qJD(3) + (-t11 * mrSges(4,3) - t3 * mrSges(5,3)) * t71 + m(5) * (t2 * t47 + t22 * t6 + t23 * t48 + t3 * t46 + t40 * t55 + t9 * t7) + t63 * t83 + t55 * t32 + t64 * t33 + t48 * t41 + t7 * t49 + t6 * t51; (-t85 * qJD(3) + t83) * pkin(6) + ((t104 * mrSges(3,1) + t126 - t42) * t72 + ((-t51 - t52) * t73 + (t49 + t50) * t71 + t104 * mrSges(3,2) - t125) * t74 + (-pkin(2) * t108 - t74 * t131 - t54 * t72) * m(4)) * t111 + (-t3 * t71 + ((-t57 * t73 - t58 * t71) * t68 + t88) * qJD(3)) * mrSges(5,3) + m(5) * (t23 * t101 + t2 * t58 + t22 * t38 + t3 * t57 + t9 * t39 + t40 * t65) + t75 + t65 * t32 - t11 * t113 + t41 * t101 - pkin(2) * t33 + t39 * t49 + t38 * t51; t11 * mrSges(4,1) - t10 * mrSges(4,2) - t2 * mrSges(5,2) + t21 * t51 + t85 * t53 + (m(5) * pkin(3) + mrSges(5,1)) * t3 + (-m(5) * (-t21 - t9) + t49) * t22 + ((-t54 * mrSges(4,2) - t23 * mrSges(5,2) - t30 / 0.2e1 - t31 / 0.2e1 - t61 / 0.2e1 - t62 / 0.2e1 + t9 * mrSges(5,3) + (Ifges(5,5) / 0.2e1 + Ifges(4,5) / 0.2e1 - pkin(3) * mrSges(5,3)) * qJD(3)) * t73 + (-t54 * mrSges(4,1) - t23 * mrSges(5,1) + t28 / 0.2e1 + t29 / 0.2e1 + t22 * mrSges(5,3) + (Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t115 + (-Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(3) + t126 * pkin(3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(4,1) / 0.2e1) * t114) * t71) * t68; m(5) * t40 + (t71 * t49 - t73 * t51 - t125) * t68 + t32;];
tauc = t1(:);
