% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:07:52
% DurationCPUTime: 1.48s
% Computational Cost: add. (893->198), mult. (2254->280), div. (0->0), fcn. (1162->4), ass. (0->93)
t73 = sin(qJ(2));
t95 = t73 * qJD(1);
t67 = pkin(5) * t95;
t129 = qJD(3) + t67;
t128 = -pkin(6) * t95 + t129;
t113 = qJD(2) - qJD(4);
t127 = -t113 / 0.2e1;
t116 = -qJD(2) / 0.2e1;
t126 = -mrSges(3,1) - mrSges(4,1);
t75 = cos(qJ(2));
t94 = t75 * qJD(1);
t68 = pkin(5) * t94;
t47 = -pkin(6) * t94 + t68;
t105 = pkin(2) + pkin(3);
t72 = sin(qJ(4));
t74 = cos(qJ(4));
t51 = -t72 * qJ(3) - t105 * t74;
t125 = t51 * qJD(4) + t128 * t74 - t47 * t72;
t52 = t74 * qJ(3) - t105 * t72;
t124 = -t52 * qJD(4) - t128 * t72 - t47 * t74;
t100 = Ifges(3,4) * t73;
t118 = -Ifges(4,5) * t95 / 0.2e1;
t119 = Ifges(4,3) / 0.2e1;
t89 = qJ(3) * t73 + pkin(1);
t54 = -pkin(2) * t75 - t89;
t39 = t54 * qJD(1);
t71 = qJD(2) * qJ(3);
t57 = t68 + t71;
t59 = mrSges(4,2) * t94 + qJD(2) * mrSges(4,3);
t123 = -(m(4) * t57 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t94 + t59) * pkin(5) - t57 * mrSges(4,2) - t94 * t119 - t118 - (Ifges(3,2) * t75 + t100) * qJD(1) / 0.2e1 + t39 * mrSges(4,1) + (-Ifges(4,6) + Ifges(3,6)) * t116;
t117 = -Ifges(3,4) * t94 / 0.2e1;
t120 = -Ifges(3,1) / 0.2e1;
t53 = -qJD(2) * pkin(2) + t129;
t99 = Ifges(4,5) * t75;
t122 = (m(4) * t53 + (mrSges(4,2) + mrSges(3,3)) * t95 + t126 * qJD(2)) * pkin(5) - t39 * mrSges(4,3) + (Ifges(4,1) * t73 - t99) * qJD(1) / 0.2e1 - t95 * t120 - t117 + t53 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t116;
t104 = pkin(5) - pkin(6);
t60 = t104 * t73;
t48 = qJD(2) * t60;
t70 = qJD(2) * qJD(3);
t27 = -qJD(1) * t48 + t70;
t61 = t104 * t75;
t49 = qJD(2) * t61;
t40 = qJD(1) * t49;
t24 = -t105 * qJD(2) + t128;
t38 = t47 + t71;
t5 = t24 * t74 - t38 * t72;
t1 = t5 * qJD(4) + t27 * t74 + t40 * t72;
t82 = t72 * t75 - t73 * t74;
t14 = t113 * t82;
t10 = t14 * qJD(1);
t37 = -t72 * t94 + t74 * t95;
t106 = t37 / 0.2e1;
t6 = t24 * t72 + t38 * t74;
t2 = -t6 * qJD(4) - t27 * t72 + t40 * t74;
t41 = t105 * t75 + t89;
t23 = t41 * qJD(1);
t36 = -t72 * t95 - t74 * t94;
t29 = Ifges(5,4) * t36;
t30 = Ifges(5,4) * t37;
t7 = Ifges(5,2) * t36 - Ifges(5,6) * t113 + t30;
t8 = Ifges(5,1) * t37 - Ifges(5,5) * t113 + t29;
t81 = t72 * t73 + t74 * t75;
t15 = t113 * t81;
t9 = t15 * qJD(1);
t121 = -t37 * (-Ifges(5,1) * t36 + t30) / 0.2e1 - t2 * mrSges(5,1) + t1 * mrSges(5,2) - Ifges(5,5) * t9 + Ifges(5,6) * t10 - t7 * t106 + t23 * (t37 * mrSges(5,1) + t36 * mrSges(5,2)) + (Ifges(5,5) * t36 - Ifges(5,6) * t37) * t127 - (Ifges(5,2) * t37 - t29 - t8) * t36 / 0.2e1;
t108 = t36 / 0.2e1;
t103 = pkin(1) * mrSges(3,1);
t102 = pkin(1) * mrSges(3,2);
t97 = qJ(3) * t75;
t96 = qJD(3) * t73;
t93 = Ifges(3,5) / 0.2e1 + Ifges(4,4) / 0.2e1;
t92 = 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t91 = -Ifges(3,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t90 = m(4) * pkin(5) + mrSges(4,2);
t88 = t36 * t5 + t37 * t6;
t84 = pkin(2) * t73 - t97;
t19 = mrSges(5,2) * t113 + t36 * mrSges(5,3);
t20 = -mrSges(5,1) * t113 - t37 * mrSges(5,3);
t83 = t74 * t19 - t72 * t20;
t17 = t60 * t74 - t61 * t72;
t18 = t60 * t72 + t61 * t74;
t80 = -t105 * t73 + t97;
t31 = t84 * qJD(2) - t96;
t21 = t80 * qJD(2) + t96;
t50 = -qJD(2) * t67 + t70;
t44 = (-t75 * mrSges(4,1) - mrSges(4,3) * t73) * qJD(1);
t28 = t80 * qJD(1);
t22 = t31 * qJD(1);
t16 = t21 * qJD(1);
t11 = -mrSges(5,1) * t36 + mrSges(5,2) * t37;
t4 = -t18 * qJD(4) + t48 * t72 + t49 * t74;
t3 = t17 * qJD(4) - t48 * t74 + t49 * t72;
t12 = [(Ifges(5,5) * t15 - Ifges(5,6) * t14) * t127 + t41 * (mrSges(5,1) * t10 + mrSges(5,2) * t9) + t16 * (mrSges(5,1) * t81 - mrSges(5,2) * t82) + t31 * t44 + t23 * (mrSges(5,1) * t14 + mrSges(5,2) * t15) - t14 * t7 / 0.2e1 + t15 * t8 / 0.2e1 + t3 * t19 + t4 * t20 + t21 * t11 + (t10 * t81 - t14 * t108) * Ifges(5,2) + (t15 * t106 - t82 * t9) * Ifges(5,1) + m(5) * (t1 * t18 + t16 * t41 + t17 * t2 + t21 * t23 + t3 * t6 + t4 * t5) + m(4) * (t22 * t54 + t31 * t39) + (-t1 * t81 - t18 * t10 - t6 * t14 - t5 * t15 - t17 * t9 + t2 * t82) * mrSges(5,3) + (t10 * t82 - t14 * t106 + t15 * t108 - t81 * t9) * Ifges(5,4) + (-t22 * mrSges(4,3) + (t91 * qJD(2) + (t54 * mrSges(4,1) + t92 * t73 - 0.2e1 * t103) * qJD(1) + t123) * qJD(2)) * t73 + (-t22 * mrSges(4,1) + t90 * t50 + (t93 * qJD(2) + (-t54 * mrSges(4,3) - 0.2e1 * t102 - t92 * t75 + (-0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) + t90 * pkin(5)) * t73) * qJD(1) + t122) * qJD(2)) * t75; m(4) * (t50 * qJ(3) + t57 * qJD(3)) + t124 * t20 + t125 * t19 + (-t10 * t52 - t51 * t9 - t88) * mrSges(5,3) + ((t117 + (t102 + t99 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t126) * pkin(5) + t93) * qJD(2) - t122) * t75 + (t118 + (t103 + t100 / 0.2e1) * qJD(1) + (-Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + t119 + t120) * t94 + (pkin(5) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t91) * qJD(2) - t123) * t73) * qJD(1) + t50 * mrSges(4,3) + qJD(3) * t59 - t28 * t11 + (-m(4) * t39 - t44) * t84 * qJD(1) + (t1 * t52 + t124 * t5 + t125 * t6 + t2 * t51 - t23 * t28) * m(5) + t121; t83 * qJD(4) + (-t11 + t44) * t95 + (-t10 * t72 - t74 * t9) * mrSges(5,3) + (t90 * t94 - t59 - t83) * qJD(2) - m(4) * (qJD(2) * t57 - t39 * t95) + (t1 * t72 + t2 * t74 - t23 * t95 - t113 * (-t5 * t72 + t6 * t74)) * m(5); t88 * mrSges(5,3) - t5 * t19 + t6 * t20 - t121;];
tauc = t12(:);
