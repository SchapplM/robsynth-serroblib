% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:34
% DurationCPUTime: 1.32s
% Computational Cost: add. (1309->202), mult. (2830->291), div. (0->0), fcn. (1421->6), ass. (0->100)
t53 = sin(qJ(5));
t102 = Ifges(6,6) * t53;
t55 = cos(qJ(5));
t103 = Ifges(6,5) * t55;
t112 = t55 / 0.2e1;
t113 = -t53 / 0.2e1;
t89 = t53 * qJD(4);
t56 = cos(qJ(4));
t92 = qJD(1) * t56;
t41 = t55 * t92 + t89;
t116 = t41 / 0.2e1;
t106 = Ifges(6,4) * t41;
t88 = t55 * qJD(4);
t40 = -t53 * t92 + t88;
t54 = sin(qJ(4));
t93 = qJD(1) * t54;
t49 = qJD(5) + t93;
t13 = Ifges(6,2) * t40 + Ifges(6,6) * t49 + t106;
t37 = Ifges(6,4) * t40;
t14 = Ifges(6,1) * t41 + Ifges(6,5) * t49 + t37;
t48 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t38 = t48 * qJD(1) + qJD(3);
t30 = -t54 * qJD(2) + t38 * t56;
t22 = -qJD(4) * pkin(4) - t30;
t104 = Ifges(6,4) * t55;
t70 = -Ifges(6,2) * t53 + t104;
t105 = Ifges(6,4) * t53;
t72 = Ifges(6,1) * t55 - t105;
t73 = mrSges(6,1) * t53 + mrSges(6,2) * t55;
t31 = qJD(2) * t56 + t54 * t38;
t23 = qJD(4) * pkin(7) + t31;
t50 = sin(pkin(8)) * pkin(1) + qJ(3);
t36 = t54 * pkin(4) - pkin(7) * t56 + t50;
t32 = t36 * qJD(1);
t7 = -t23 * t53 + t32 * t55;
t8 = t23 * t55 + t32 * t53;
t75 = t53 * t8 + t55 * t7;
t131 = -t75 * mrSges(6,3) + t14 * t112 + t13 * t113 + (-t102 + t103) * t49 / 0.2e1 + t72 * t116 + t40 * t70 / 0.2e1 + t22 * t73;
t28 = -mrSges(6,2) * t49 + mrSges(6,3) * t40;
t29 = mrSges(6,1) * t49 - mrSges(6,3) * t41;
t123 = -m(6) * t75 - t53 * t28 - t55 * t29;
t87 = qJD(4) * qJD(5);
t90 = qJD(5) * t56;
t26 = t55 * t87 + (-t53 * t90 - t54 * t88) * qJD(1);
t91 = qJD(4) * t56;
t83 = qJD(1) * t91;
t15 = mrSges(6,1) * t83 - t26 * mrSges(6,3);
t27 = -t53 * t87 + (t54 * t89 - t55 * t90) * qJD(1);
t16 = -mrSges(6,2) * t83 + t27 * mrSges(6,3);
t130 = t123 * qJD(5) - t53 * t15 + t55 * t16;
t107 = Ifges(5,4) * t54;
t95 = Ifges(5,5) * qJD(4);
t129 = t95 / 0.2e1 + (t56 * Ifges(5,1) - t107) * qJD(1) / 0.2e1 - t30 * mrSges(5,3) + t131;
t25 = qJD(4) * t31;
t9 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t128 = t9 + m(6) * (t7 * t89 - t8 * t88 + t25);
t81 = -Ifges(5,6) * qJD(4) / 0.2e1;
t127 = m(4) + m(5);
t24 = qJD(4) * t30;
t77 = pkin(4) * t56 + pkin(7) * t54;
t39 = t77 * qJD(4) + qJD(3);
t33 = t39 * qJD(1);
t1 = qJD(5) * t7 + t24 * t55 + t33 * t53;
t2 = -qJD(5) * t8 - t24 * t53 + t33 * t55;
t76 = t1 * t55 - t2 * t53;
t124 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t26 + Ifges(6,6) * t27;
t121 = t26 / 0.2e1;
t120 = t27 / 0.2e1;
t119 = -t40 / 0.2e1;
t117 = -t41 / 0.2e1;
t115 = -t49 / 0.2e1;
t99 = t48 * t54;
t84 = mrSges(5,3) * t92;
t97 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t40 - mrSges(6,2) * t41 - t84;
t96 = mrSges(4,3) * qJD(1);
t44 = qJD(1) * t50;
t86 = mrSges(5,3) * t93;
t85 = t48 * t91;
t82 = -t95 / 0.2e1;
t80 = 0.2e1 * t44;
t74 = t55 * mrSges(6,1) - t53 * mrSges(6,2);
t71 = Ifges(6,1) * t53 + t104;
t69 = Ifges(6,2) * t55 + t105;
t67 = Ifges(6,5) * t53 + Ifges(6,6) * t55;
t18 = t36 * t53 + t55 * t99;
t17 = t36 * t55 - t53 * t99;
t45 = -qJD(4) * mrSges(5,2) - t86;
t64 = t28 * t55 - t29 * t53 + t45;
t59 = t44 * mrSges(5,1) - t31 * mrSges(5,3) + t49 * Ifges(6,3) + t40 * Ifges(6,6) + t41 * Ifges(6,5) + t7 * mrSges(6,1) - t8 * mrSges(6,2) + t81 - (Ifges(5,4) * t56 - Ifges(5,2) * t54) * qJD(1) / 0.2e1;
t58 = m(6) * (qJD(4) * t22 + t76) + t130;
t47 = Ifges(6,3) * t83;
t43 = t77 * qJD(1);
t42 = qJD(1) * (t54 * mrSges(5,1) + mrSges(5,2) * t56);
t11 = t30 * t55 + t43 * t53;
t10 = -t30 * t53 + t43 * t55;
t6 = Ifges(6,1) * t26 + Ifges(6,4) * t27 + Ifges(6,5) * t83;
t5 = Ifges(6,4) * t26 + Ifges(6,2) * t27 + Ifges(6,6) * t83;
t4 = -t18 * qJD(5) + t55 * t39 - t53 * t85;
t3 = t17 * qJD(5) + t53 * t39 + t55 * t85;
t12 = [m(6) * (t1 * t18 + t2 * t17 + t8 * t3 + t7 * t4) + t3 * t28 + t4 * t29 + t17 * t15 + t18 * t16 + (t127 * t80 + t42 + 0.2e1 * t96) * qJD(3) + (qJD(1) * qJD(3) * mrSges(5,1) + t47 / 0.2e1 + (m(5) * t48 - mrSges(5,3)) * t24 + (0.3e1 / 0.2e1 * Ifges(5,4) * t93 + t82 - t80 * mrSges(5,2) + (-m(5) * t30 + m(6) * t22 - t97) * t48 - t129) * qJD(4) + t124) * t54 + (t72 * t121 + t70 * t120 + t6 * t112 + t5 * t113 + (mrSges(5,3) + t73) * t25 + (-t1 * t53 - t2 * t55) * mrSges(6,3) + (t81 + t59) * qJD(4) + (-t55 * t13 / 0.2e1 + t14 * t113 + t67 * t115 + t69 * t119 + t71 * t117 + t22 * t74 + (t7 * t53 - t8 * t55) * mrSges(6,3)) * qJD(5) + (-m(6) * t25 + qJD(4) * t45 - t9) * t48 + (qJD(3) * mrSges(5,2) + ((-0.3e1 / 0.2e1 * Ifges(5,4) + t103 / 0.2e1 - t102 / 0.2e1) * t56 + t50 * mrSges(5,1) + (0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,1)) * t54) * qJD(4)) * qJD(1)) * t56; ((-t64 - t86) * qJD(4) + t128) * t54 + ((-t84 - t97) * qJD(4) + t58) * t56; (t64 * qJD(4) - t128) * t56 + (-t97 * qJD(4) + t58) * t54 + (-t127 * t44 + t123 - t42 - t96) * qJD(1); t71 * t121 + t69 * t120 + t5 * t112 + t53 * t6 / 0.2e1 - t30 * t45 - t24 * mrSges(5,2) - t11 * t28 - t10 * t29 - pkin(4) * t9 + t97 * t31 + (-mrSges(5,1) - t74) * t25 + t76 * mrSges(6,3) + t131 * qJD(5) + ((t81 + qJD(4) * t67 / 0.2e1 + Ifges(5,4) * t92 / 0.2e1 - t59) * t56 + (t44 * mrSges(5,2) + t82 + (-t107 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t56) * qJD(1) + t129) * t54) * qJD(1) + (-pkin(4) * t25 - t10 * t7 - t11 * t8 - t22 * t31) * m(6) + (m(6) * t76 + t130) * pkin(7); t47 - t22 * (mrSges(6,1) * t41 + mrSges(6,2) * t40) + (Ifges(6,1) * t40 - t106) * t117 + t13 * t116 + (Ifges(6,5) * t40 - Ifges(6,6) * t41) * t115 - t7 * t28 + t8 * t29 + (t40 * t7 + t41 * t8) * mrSges(6,3) + (-Ifges(6,2) * t41 + t14 + t37) * t119 + t124;];
tauc = t12(:);
