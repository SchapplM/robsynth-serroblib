% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:36
% DurationCPUTime: 1.30s
% Computational Cost: add. (740->173), mult. (1994->252), div. (0->0), fcn. (1128->4), ass. (0->98)
t122 = -Ifges(5,4) - Ifges(6,4);
t121 = Ifges(5,1) + Ifges(6,1);
t98 = Ifges(6,5) + Ifges(5,5);
t120 = -Ifges(5,2) - Ifges(6,2);
t97 = Ifges(5,6) + Ifges(6,6);
t58 = cos(qJ(4));
t119 = t122 * t58;
t57 = sin(qJ(4));
t118 = t122 * t57;
t55 = sin(pkin(8));
t51 = t55 ^ 2;
t56 = cos(pkin(8));
t52 = t56 ^ 2;
t117 = (t51 + t52) * mrSges(4,3);
t116 = t97 * t58;
t42 = -pkin(3) * t56 - pkin(6) * t55 - pkin(2);
t100 = t56 * t58;
t48 = qJ(3) * t100;
t21 = t57 * t42 + t48;
t49 = -qJD(2) * t56 + qJD(4);
t83 = t55 * qJD(2);
t74 = t58 * t83;
t29 = mrSges(6,1) * t49 - mrSges(6,3) * t74;
t30 = mrSges(5,1) * t49 - mrSges(5,3) * t74;
t95 = t29 + t30;
t76 = t57 * t83;
t27 = -mrSges(6,2) * t49 - mrSges(6,3) * t76;
t28 = -mrSges(5,2) * t49 - mrSges(5,3) * t76;
t96 = t27 + t28;
t115 = t95 * t57 - t96 * t58;
t107 = mrSges(6,2) * t57;
t108 = mrSges(6,1) * t58;
t110 = -t58 / 0.2e1;
t39 = (pkin(4) * t57 + qJ(3)) * t55;
t50 = t56 * qJD(1);
t22 = qJD(2) * t39 + qJD(5) - t50;
t82 = qJ(3) * qJD(2);
t40 = t55 * t82 - t50;
t34 = t42 * qJD(2) + qJD(3);
t41 = qJD(1) * t55 + t56 * t82;
t10 = t58 * t34 - t41 * t57;
t11 = t34 * t57 + t41 * t58;
t63 = t10 * t57 - t11 * t58;
t65 = mrSges(5,1) * t58 - mrSges(5,2) * t57;
t72 = qJ(5) * t83;
t8 = -t58 * t72 + t10;
t3 = pkin(4) * t49 + t8;
t9 = -t57 * t72 + t11;
t67 = t3 * t57 - t9 * t58;
t114 = t63 * mrSges(5,3) + t67 * mrSges(6,3) + t40 * t65 + t22 * (-t107 + t108) + (-(t121 * t58 + t118) * t57 / 0.2e1 + (t120 * t57 - t119) * t110) * t83 + (-t98 * t57 - t116 / 0.2e1 + t97 * t110) * t49;
t86 = qJD(4) * t58;
t36 = (pkin(4) * t86 + qJD(3)) * t55;
t33 = qJD(2) * t36;
t109 = m(6) * t33;
t102 = t40 * t55;
t101 = t56 * t57;
t99 = t58 * t55;
t88 = qJD(3) * t56;
t94 = t42 * t86 + t58 * t88;
t91 = qJ(3) * t57;
t90 = qJ(5) * t55;
t89 = qJD(3) * t55;
t87 = qJD(4) * t57;
t85 = qJD(5) * t57;
t84 = qJD(5) * t58;
t81 = qJ(5) * qJD(4);
t80 = qJD(2) * qJD(3);
t79 = qJD(2) * qJD(4);
t78 = t56 * t91;
t77 = t58 * t90;
t75 = t57 * t88;
t73 = 0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,4);
t71 = t56 * t80;
t70 = t55 * t79;
t69 = t51 * qJ(3) * t80 + t40 * t89;
t6 = t34 * t86 - t41 * t87 + t58 * t71;
t1 = (-t58 * t81 - t85) * t83 + t6;
t62 = t11 * qJD(4);
t2 = -t62 + (-t75 + (t57 * t81 - t84) * t55) * qJD(2);
t68 = -t1 * t57 - t2 * t58;
t7 = -t57 * t71 - t62;
t66 = -t6 * t57 - t7 * t58;
t64 = mrSges(6,1) * t57 + mrSges(6,2) * t58;
t61 = t7 * mrSges(5,1) - t6 * mrSges(5,2) - t1 * mrSges(6,2);
t43 = t70 * t108;
t38 = t58 * t42;
t32 = (mrSges(5,1) * t57 + mrSges(5,2) * t58) * t83;
t31 = t64 * t83;
t25 = t65 * t70;
t24 = -t70 * t107 + t43;
t20 = t38 - t78;
t15 = -t57 * t90 + t21;
t14 = -t21 * qJD(4) - t75;
t13 = -qJD(4) * t78 + t94;
t12 = -t77 + t38 + (-pkin(4) - t91) * t56;
t5 = -t75 - t55 * t84 + (-t48 + (-t42 + t90) * t57) * qJD(4);
t4 = -t55 * t85 + (-t77 - t78) * qJD(4) + t94;
t16 = [(-t24 - t25 - t109) * t56 + ((-t96 * t57 - t95 * t58) * qJD(4) + m(5) * (-t10 * t86 - t11 * t87 - t57 * t7 + t58 * t6 - t71) + m(6) * (t1 * t58 - t2 * t57 - t3 * t86 - t9 * t87)) * t55 + (mrSges(5,3) + mrSges(6,3)) * t51 * t79 * (-t57 ^ 2 - t58 ^ 2); t13 * t28 + t14 * t30 + t39 * t24 + t4 * t27 + t5 * t29 + t36 * t31 + (-t2 * mrSges(6,1) - t61) * t56 + (t66 * mrSges(5,3) + t68 * mrSges(6,3) + qJ(3) * t25 + qJD(3) * t32 + t33 * t64) * t55 + m(4) * (t41 * t88 + t69) + m(6) * (t1 * t15 + t12 * t2 + t22 * t36 + t3 * t5 + t33 * t39 + t4 * t9) + m(5) * (t10 * t14 + t11 * t13 + t20 * t7 + t21 * t6 + t69) + t114 * t55 * qJD(4) + ((m(4) * qJ(3) * t52 + 0.2e1 * t117) * qJD(3) + ((mrSges(5,2) * t89 + (-t21 * mrSges(5,3) - t15 * mrSges(6,3) + t97 * t56 - t73 * t99) * qJD(4)) * t58 + (mrSges(5,1) * t89 + (t12 * mrSges(6,3) + t20 * mrSges(5,3) + t73 * t55 * t57 + t98 * t56 + (-0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t99) * qJD(4)) * t57) * t55) * qJD(2); -t115 * qJD(4) + m(5) * (-t63 * qJD(4) - t66) + m(6) * (-t67 * qJD(4) - t68) + ((-t31 - t32) * t55 + t115 * t56 - m(4) * (t41 * t56 + t102) - m(5) * (-t10 * t101 + t11 * t100 + t102) - m(6) * (t9 * t100 - t3 * t101 + t22 * t55) - qJD(2) * t117) * qJD(2); -t10 * t28 + t11 * t30 - t8 * t27 + (m(6) * pkin(4) + mrSges(6,1)) * t2 + (t29 - m(6) * (-t3 + t8)) * t9 + (((t120 * t58 + t118) * t57 / 0.2e1 + (-t121 * t57 + t119) * t110) * t83 + (-m(6) * t22 - t31) * t58 * pkin(4) + (-t116 + (mrSges(6,3) * pkin(4) - t98) * t57) * qJD(4) - t114) * t83 + t61; t43 + t109 + (-mrSges(6,2) * t87 + t57 * t27 + t58 * t29 - m(6) * (-t3 * t58 - t57 * t9)) * t83;];
tauc = t16(:);
